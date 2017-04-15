extern crate bit_vec;
use self::bit_vec::BitVec;

use std::fs;
use std::fs::File;
use std::io::{Write, Read, Seek, SeekFrom};

#[derive(Debug)]
#[derive(Default)]
pub struct SeriesBlock {
    t0: u64,         // beginning timestamp of this block (
    bv: BitVec,      // the stream of bits
    length: usize,   // total number of values written to this block

    state: SeriesState,
}

#[derive(Debug)]
#[derive(Default)]
struct SeriesState {
    prev_bits: u64,                // the bits for the value we emitted previously
    prev_nibbles: u8,              // the number of zero nibbles for the previous value
    prev_timestamp: u64,
    prev_delta: u64,
}

#[derive(Debug)]
#[derive(PartialEq)]
pub struct Point {
    pub t: u64,
    pub v: f64
}

struct SeriesBlockIterator<'a> {
    series_block: &'a SeriesBlock, // a ref to a SeriesBlock
    state: SeriesState,
    index: usize                   // our index in the bit stream
}

impl SeriesBlock {
    pub fn new(t0: u64) -> SeriesBlock {
        SeriesBlock {
            t0: t0,
            state: SeriesState { prev_nibbles: 0xff, ..Default::default() },
            ..Default::default()
        }
    }

    pub fn start_time(&self) -> u64 {
        self.t0
    }

    pub fn length(&self) -> usize {
        self.length
    }

    pub fn append(&mut self, p: Point) {
        self.length += 1;

        // we need a bitwise representation of our float
        let bits = unsafe { std::mem::transmute::<f64,u64>(p.v) };

        // if this is the first point, write 12 bits of the timestamp
        if self.bv.len() == 0 {
            let tdelta = p.t - self.t0;

            self.state.prev_timestamp = p.t;
            self.state.prev_bits = bits;
            self.state.prev_delta = tdelta;

            self.append_bits(tdelta, 12);
            self.append_bits(bits, 64);
            return;
        }

        // if this is any other point, calculate the double delta of the timestamp
        let tdelta = p.t - self.state.prev_timestamp;
        let dd: i64 = tdelta as i64 - self.state.prev_delta as i64;
        self.state.prev_delta = tdelta;
        self.state.prev_timestamp = p.t;

        // check to see what range the delta fell into
        match dd {
            0 => {
                self.bv.push(false);
            },
            -63 ... 64 => {
                self.append_bits(0b10, 2);
                self.append_bits(dd as u64, 7);
            },
            -255 ... 256 => {
                self.append_bits(0b110, 3);
                self.append_bits(dd as u64, 9);
            },
            -2047 ... 2048 => {
                self.append_bits(0b1110, 4);
                self.append_bits(dd as u64, 12);
            },
            _ => {
                self.append_bits(0b1111, 4);
                self.append_bits(dd as u64, 32);
            }
        }

        // get the XOR of the current value and the previous value
        let xor = self.state.prev_bits ^ bits;
        self.state.prev_bits = bits;

        // control code 00 = repeated point
        if xor == 0 {
            self.append_bits(0b0, 1);
            //self.append_bits(0b00, 2);
            return;
        }

        let lz = xor.leading_zeros() as u8;
        let zero_nibbles = lz.wrapping_div(4);
        let remaining_bits = 64 - (4 * zero_nibbles);

        // if we haven't recorded the leading zeros yet, do that now
        if (self.state.prev_nibbles != zero_nibbles) || (self.state.prev_nibbles == 0xff) {
            // control code 11 = reset lz
            // TODO: decide if it's efficient to reset
            // maybe just if lz > s.state.prev_nibbles ?
            self.state.prev_nibbles = zero_nibbles;
            self.append_bits(0b11, 2);
            self.append_bits(zero_nibbles as u64, 4);
            self.append_bits(xor, remaining_bits);
        } else {
            // control code 10 = continue using current lz
            self.append_bits(0b10, 2);
            self.append_bits(xor, remaining_bits);
        }

        // control code 01 = null value
        // TODO: this needs to get handled outside of this function

        return;
    }

    fn append_bits(&mut self, value: u64, num_bits: u8) {
        // TODO: bail out if n > 64, otherwise we'll overflow
        for x in (0..num_bits).rev() {
            self.bv.push(value & (1 << x) > 0);
        }
    }

    pub fn append_null(&mut self) {
        self.length += 1;
        self.append_bits(0b01, 2);
    }

    fn get_bits(&self, index: &mut usize, num_bits: u8) -> u64 {
        let mut bits: u64 = 0;
        for x in (0..num_bits).rev() {
            let value = self.bv.get(*index).unwrap_or(false);
            if value { bits ^= 1 << x };
            *index += 1;
        }

        return bits;
    }

    fn get_bits_until_zero(&self, index: &mut usize, max: usize) -> u64 {
        let mut bits: u64 = 0;
        let mut n: usize = 0;

        while n < max {
            let value = self.bv.get(*index).unwrap_or(false);
            *index += 1;

            bits = bits << 1;
            if value { bits += 1; } else { break; }
            n += 1;
        }

        return bits;
    }

    fn extract_timestamp_delta(&self, index: &mut usize) -> Option<u64> {
        let mut timestamp: u64 = 0;

        if *index == 0 {
            timestamp = self.get_bits(index, 12);
        } else if *index >= self.bv.len() {
            return None;
        } else {
            let control_code = self.get_bits_until_zero(index, 4);
            let mut xor: u64 = 0;

            match control_code {
                0b0 => {
                    timestamp = 0;
                }
                0b10 => {
                    timestamp = self.get_bits(index, 7);
                    if timestamp & (1<<6) != 0 { timestamp = timestamp.wrapping_sub(1<<7); }
                }
                0b110 => {
                    timestamp = self.get_bits(index, 9);
                    if timestamp & (2^9) != 0 { timestamp.wrapping_sub(2^9); }
                }
                0b1110 => {
                    timestamp = self.get_bits(index, 12);
                    if timestamp & (2^12) != 0 { timestamp.wrapping_sub(2^12); }
                }
                0b1111 => {
                    timestamp = self.get_bits(index, 32);
                }
                _ => {
                    // this should theoretically be impossible
                    return None;
                }
            }
        }

        return Some(timestamp);
    }

    fn extract_value(&self, index: &mut usize, prev_nibbles: &mut u8) -> Option<u64> {
        let mut xor: u64 = 0;

        // TODO: this is hacky. we need a better start condition.
        if *index == 12 {
            xor = self.get_bits(index, 64);
        } else if *index >= self.bv.len() {
            return None;
        } else {
            // read the first bit; if it's `0` we have a duplicate value
            if self.get_bits(index, 1) == 0b0 {
                xor = 0;
            } else {
                // if it's a `1`, check the next bit to see if we can keep the same nibble count
                if self.get_bits(index, 1) == 0b0 {
                    xor = self.get_bits(index, 64 - 4 * *prev_nibbles);
                } else {
                    let nibbles = self.get_bits(index, 4) as u8;
                    xor = self.get_bits(index, 64 - 4 * nibbles);
                    *prev_nibbles = nibbles;
                }
            }
        }

        return Some(xor);
    }

    fn iter(&self) -> SeriesBlockIterator {
        SeriesBlockIterator {
            series_block: self,
            index: 0,
            state: SeriesState {
                prev_bits: 0,
                prev_nibbles: 0,
                prev_timestamp: 0,
                prev_delta: 0
            }
        }
    }

    pub fn write_to_disk(&self, file: &mut File) {
        file.write_all(&self.bv.to_bytes());
    }
}

impl<'a> Iterator for SeriesBlockIterator<'a> {
    type Item = Point;

    fn next(&mut self) -> Option<Point> {
        let mut t: u64 = 0;
        let mut v: f64 = 0f64;

        let delta = self.series_block.extract_timestamp_delta(&mut self.index);
        match delta {
            None => { return None; },
            Some(x) => {
                let new_delta = self.state.prev_delta.wrapping_add(x);
                t = self.state.prev_timestamp.wrapping_add(new_delta);
                self.state.prev_delta = new_delta;
                self.state.prev_timestamp = t;
            }
        }

        let value = self.series_block.extract_value(&mut self.index, &mut self.state.prev_nibbles);
        match value {
            None => { return None; },
            Some(x) => {
                let bits = x ^ self.state.prev_bits;
                self.state.prev_bits = bits;
                v = unsafe { std::mem::transmute::<u64,f64>(bits) };
            }
        }

        return Some(Point {t: t, v: v});
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    extern crate rand;
    extern crate num;

    use self::rand::distributions::{Normal, IndependentSample};

    fn bits(v: f64) -> u64 {
        unsafe { std::mem::transmute::<f64,u64>(v) }
    }

    #[test]
    fn test_new_series_block() {
        let s = SeriesBlock::new(10000);

        // ensure that everything was initialized correctly
        assert_eq!(s.t0, 10000u64);
        assert_eq!(s.state.prev_nibbles, 0xff);
        assert_eq!(s.bv.len(), 0);
        assert_eq!(s.length, 0);
        assert_eq!(s.state.prev_bits, 0);
        assert_eq!(s.state.prev_timestamp, 0);
        assert_eq!(s.state.prev_delta, 0);
    }

    #[test]
    fn test_append_series_block() {
        let mut s = SeriesBlock::new(0);

        /// Point #0
        /// t: 60
        /// v: 3.14 -> 0x40091eb851eb851f
        ///
        /// since this is the first point, we should encode the last 12 bits of the timestamp
        /// and the entire value, for a total of 12 + 64 = 72 bits
        s.append(Point{t: 60, v: 3.14});
        assert_eq!(s.length, 1);
        assert_eq!(s.state.prev_bits, bits(3.14));
        assert_eq!(s.state.prev_timestamp, 60);
        assert_eq!(s.state.prev_delta, 60);
        assert_eq!(s.bv.len(), 76);
        assert_eq!(s.state.prev_nibbles, 0xff);

        /// Point #1
        /// t: 120 (delta of 60 from the previous point)
        /// v: 3.15 -> 0x4009333333333333 xor with previous = 18 lz = ~4 nibbles
        ///
        /// this point should have the same delta as the previous point, so we should encode
        /// a zero for the timestamp indicating a double delta of 0. the value will have 4
        /// nibbles for the leading zeros (16 bits), leaving 48 significant bits plus the four
        /// bits storing the nibble count.
        s.append(Point{t: 120, v: 3.15});
        assert_eq!(s.length, 2);
        assert_eq!(s.state.prev_bits, bits(3.15));
        assert_eq!(s.state.prev_timestamp, 120);
        assert_eq!(s.state.prev_delta, 60);
        assert_eq!(s.bv.len(), 131); // previous + 1 + 2 + 4 + 48
        assert_eq!(s.state.prev_nibbles, 4);

        /// Point #2
        /// t: 180 (delta of 60 from the previous point)
        /// v: 3.15 -> 0x4009333333333333 xor with previous = 18 lz = ~4 nibbles
        ///
        /// this point should again have the same delta as the previous point, so we should encode
        /// a zero for the timestamp indicating a double delta of 0. the value is exactly the same
        /// as the previous value, so we should just encode 2 bits for a repeat value.
        s.append(Point{t: 180, v: 3.15});
        assert_eq!(s.length, 3);
        assert_eq!(s.state.prev_bits, bits(3.15));
        assert_eq!(s.state.prev_timestamp, 180);
        assert_eq!(s.state.prev_delta, 60);
        assert_eq!(s.bv.len(), 133); // previous + 1 + 1
        assert_eq!(s.state.prev_nibbles, 4);

        /// Point #3
        /// t: 181 (delta of 1 from the previous point)
        /// v: 3.1501 -> 0x40093367a0f9096c xor with previous = 25 lz = ~6 nibbles
        ///
        /// this point has a delta of 1 from the previous point, so we encode `0b10` to indicate
        /// the range -63..64 for the double delta. we'll write only the 7 LSB of the timestamp.
        /// the value is slightly different than the previous, so we have 6 nibbles of leading
        /// zero overlap, leaving 40 significant bits.
        s.append(Point{t: 181, v: 3.1501});
        assert_eq!(s.length, 4);
        assert_eq!(s.state.prev_bits, bits(3.1501));
        assert_eq!(s.state.prev_timestamp, 181);
        assert_eq!(s.state.prev_delta, 1);
        assert_eq!(s.bv.len(), 188); // previous + 9 + 1 + 4 + 40
        assert_eq!(s.state.prev_nibbles, 6);

        /// Point #4
        /// t: 182 (delta of 1 from the previous point)
        /// v: 3.1502 -> 0x4009339c0ebedfa4 xor with previous = 24 lz = ~6 nibbles
        ///
        /// this point should again have the same delta as the previous point, so we should encode
        /// a zero for the timestamp indicating a double delta of 0. the value retains the same
        /// number of leading zero nibbles as the previous point, so we write a `0b10` to indicate
        /// that the leading zero nibbles are unchanged, and then write our 40 significant bits.
        s.append(Point{t: 182, v: 3.1502});
        assert_eq!(s.length, 5);
        assert_eq!(s.state.prev_bits, bits(3.1502));
        assert_eq!(s.state.prev_timestamp, 182);
        assert_eq!(s.state.prev_delta, 1);
        assert_eq!(s.bv.len(), 231); // previous + 1 + 1 + 40
        assert_eq!(s.state.prev_nibbles, 6);

        /// Point #5
        /// t: 240 (delta of 58 from the previous point)
        /// v: 4.20 -> 0x4010cccccccccccd xor with previous = 11 lz = ~2 nibbles
        ///
        /// this point has a delta of 58 from the previous point so we encode `0b10` to indicate
        /// the range -63..64 for the double delta. we'll write only the 7 LSB of the timestamp.
        /// the value is much higher than the previous, we'll write a `0b11` plus the value `2` in
        /// four bits to indicate the new leading zero nibble count, then the 56 significant bits.
        s.append(Point{t: 240, v: 4.20});
        assert_eq!(s.length, 6);
        assert_eq!(s.state.prev_bits, bits(4.20));
        assert_eq!(s.state.prev_timestamp, 240);
        assert_eq!(s.state.prev_delta, 58);
        assert_eq!(s.bv.len(), 302); // previous + 2 + 7 + 2 + 4 + 56
        assert_eq!(s.state.prev_nibbles, 2);
    }

    #[test]
/*
    fn test_extract_series_block() {
        let mut s = SeriesBlock::new(0);
        let mut index = 0;

        s.append(Point{t: 60, v: 3.14}); // 0x40091eb851eb851f
        s.append(Point{t: 120, v: 3.15}); // 0x4009333333333333
        s.append(Point{t: 180, v: 3.15}); // 0x4009333333333333
        s.append(Point{t: 181, v: 3.1501}); // 0x40093367a0f9096c
        s.append(Point{t: 182, v: 3.1502}); // 0x4009339c0ebedfa4
        s.append(Point{t: 240, v: 4.20}); // 0x4010cccccccccccd

        let value = s.extract(&mut index, 0, 0);
        assert_eq!(value, 3.14);

        let previous = unsafe { std::mem::transmute::<f64, u64>(value) };
        let value = s.extract(&mut index, previous, 0);
        assert_eq!(value, 3.15);

        let previous = unsafe { std::mem::transmute::<f64, u64>(value) };
        let value = s.extract(&mut index, previous, 4);
        assert_eq!(value, 3.15);

        let previous = unsafe { std::mem::transmute::<f64, u64>(value) };
        let value = s.extract(&mut index, previous, 4);
        assert_eq!(value, 3.1501);

        let previous = unsafe { std::mem::transmute::<f64, u64>(value) };
        let value = s.extract(&mut index, previous, 6);
        assert_eq!(value, 3.1502);

        let previous = unsafe { std::mem::transmute::<f64, u64>(value) };
        let value = s.extract(&mut index, previous, 6);
        assert_eq!(value, 4.20);
    }
*/
    #[test]
    fn test_series_block_iterator() {
        let mut s = SeriesBlock::new(0);

        s.append(Point{t: 60, v: 3.14}); // 0x40091eb851eb851f
        s.append(Point{t: 120, v: 3.15}); // 0x4009333333333333
        s.append(Point{t: 180, v: 3.15}); // 0x4009333333333333
        s.append(Point{t: 181, v: 3.1501}); // 0x40093367a0f9096c
        s.append(Point{t: 182, v: 3.1502}); // 0x4009339c0ebedfa4
        s.append(Point{t: 240, v: 4.20}); // 0x4010cccccccccccd

        s.append(Point{t: 300, v: 4.19});
        s.append(Point{t: 330, v: 4.18});
        s.append(Point{t: 360, v: 4.20});
        s.append(Point{t: 420, v: 4.40});

        let mut iter = s.iter();

        assert_eq!(iter.next(), Some(Point{t: 60, v: 3.14}));
        assert_eq!(iter.next(), Some(Point{t: 120, v: 3.15}));
        assert_eq!(iter.next(), Some(Point{t: 180, v: 3.15}));
        assert_eq!(iter.next(), Some(Point{t: 181, v: 3.1501}));
        assert_eq!(iter.next(), Some(Point{t: 182, v: 3.1502}));
        assert_eq!(iter.next(), Some(Point{t: 240, v: 4.20}));

        assert_eq!(iter.next(), Some(Point{t: 300, v: 4.19}));
        assert_eq!(iter.next(), Some(Point{t: 330, v: 4.18}));
        assert_eq!(iter.next(), Some(Point{t: 360, v: 4.20}));
        assert_eq!(iter.next(), Some(Point{t: 420, v: 4.40}));

        assert_eq!(iter.next(), None);
    }

    #[test]
    fn test_series_block_massive() {
        use std::usize;

        const LENGTH: usize = 24*60;

        struct SeriesBlockTest<'a> {
            series_block: SeriesBlock,
            inputs: [f64; LENGTH],
            name: &'a str
        }

        impl<'a> SeriesBlockTest<'a> {
            pub fn new<F: FnOnce(&mut [f64])>(name: &'a str, assign: F) -> SeriesBlockTest<'a> {
                let mut s = SeriesBlockTest {
                    name: name,
                    series_block: SeriesBlock::new(0),
                    inputs: [0f64; LENGTH]
                };

                assign(&mut s.inputs);

                return s;
            }
        }

        struct TestContainer<'a> {
            tests: Vec<SeriesBlockTest<'a>>
        }

        impl<'a> TestContainer<'a> {
            pub fn new() -> TestContainer<'a> {
                TestContainer{ tests: Vec::new() }
            }

            pub fn add<F: FnOnce(&mut [f64])>(&mut self, name: &'a str, assign: F) {
                let s = SeriesBlockTest::new(name, assign);
                self.tests.push(s);
            }
        }

        // helper function to round floats to a fixed number of decimal places
        fn round_f64(value: f64, places: i32) -> f64 {
            if places <= 0 {
                return value.round();
            } else {
                return (value * 10.0_f64.powi(places)).round() / 10.0_f64.powi(places);
            }
        }

        // arrays to hold random number sequences
        let mut rand_f64: [f64; LENGTH] = [0f64; LENGTH];
        let mut normal_f64: [f64; LENGTH] = [0f64; LENGTH];
        let mut normal10_f64: [f64; LENGTH] = [0f64; LENGTH];

        // fill the array with pure random numbers
        for f in rand_f64.iter_mut() {
            *f = rand::random::<f64>();
        };

        // fill the array with numbers that follow a normal distribution (μ=1.0, σ=0.1)
        let normal = Normal::new(1.0, 0.1);
        for f in normal_f64.iter_mut() {
            *f = normal.ind_sample(&mut rand::thread_rng());
        };

        // fill the array with numbers that follow a normal distribution (μ=10.0, σ=0.1)
        let normal10 = Normal::new(10.0, 0.1);
        for f in normal10_f64.iter_mut() {
            *f = normal10.ind_sample(&mut rand::thread_rng());
        };

        let mut tc = TestContainer::new();

        tc.add("constant zero",              |ii| { for i in ii { *i = 0f64; } });
        tc.add("constant one",               |ii| { for i in ii { *i = 1f64; } });
        tc.add("constant +.3f",              |ii| { for i in ii { *i = 1234.567f64; } });
        tc.add("constant -.3f",              |ii| { for i in ii { *i = -1234.567f64; } });
        tc.add("constant +.0f",              |ii| { for i in ii { *i = 1234f64; } });
        tc.add("constant -.0f",              |ii| { for i in ii { *i = -1234f64; } });
        tc.add("constant f64::MAX",          |ii| { for i in ii { *i = std::f64::MAX; } });
        tc.add("constant f64::MIN",          |ii| { for i in ii { *i = std::f64::MIN; } });
        tc.add("constant f64::MIN_POSITIVE", |ii| { for i in ii { *i = std::f64::MIN_POSITIVE; } });
        tc.add("alternating 1 or 0 (n=1)",   |ii| { for (n, i) in ii.into_iter().enumerate() { if n.wrapping_rem(2) < 1 { *i = 0f64; } else { *i = 1f64; } } });
        tc.add("alternating 1 or 0 (n=100)", |ii| { for (n, i) in ii.into_iter().enumerate() { if n.wrapping_rem(200) < 100 { *i = 0f64; } else { *i = 1f64; } } });
        tc.add("random number",              |ii| { for (n, i) in ii.into_iter().enumerate() { *i = rand_f64[n]; } });
        tc.add("random (μ=1.0, σ=0.1)",      |ii| { for (n, i) in ii.into_iter().enumerate() { *i = normal_f64[n]; } });
        tc.add("random (μ=1.0, σ=0.1) .3f",  |ii| { for (n, i) in ii.into_iter().enumerate() { *i = round_f64(normal_f64[n], 3); } });
        tc.add("random (μ=1.0, σ=0.1) .2f",  |ii| { for (n, i) in ii.into_iter().enumerate() { *i = round_f64(normal_f64[n], 2); } });
        tc.add("random (μ=1.0, σ=0.1) .1f",  |ii| { for (n, i) in ii.into_iter().enumerate() { *i = round_f64(normal_f64[n], 1); } });
        tc.add("random (μ=1.0, σ=0.1) .0f",  |ii| { for (n, i) in ii.into_iter().enumerate() { *i = round_f64(normal_f64[n], 0); } });
        tc.add("random (μ=10.0, σ=0.1)",     |ii| { for (n, i) in ii.into_iter().enumerate() { *i = normal10_f64[n]; } });
        tc.add("random (μ=10.0, σ=0.1) .3f", |ii| { for (n, i) in ii.into_iter().enumerate() { *i = round_f64(normal10_f64[n], 3); } });
        tc.add("random (μ=10.0, σ=0.1) .2f", |ii| { for (n, i) in ii.into_iter().enumerate() { *i = round_f64(normal10_f64[n], 2); } });
        tc.add("random (μ=10.0, σ=0.1) .1f", |ii| { for (n, i) in ii.into_iter().enumerate() { *i = round_f64(normal10_f64[n], 1); } });
        tc.add("random (μ=10.0, σ=0.1) .0f", |ii| { for (n, i) in ii.into_iter().enumerate() { *i = round_f64(normal10_f64[n], 0); } });

        // display the test header
        println!("{:<32} {:>10} {:>8}", "test case", "length", "bpp");

        for mut test in tc.tests {
            // compress the test data
            let mut t = 0u64;
            for input in test.inputs.into_iter() {
                t += 60;
                test.series_block.append(Point{t: t, v: *input});
            }

            assert_eq!(test.series_block.length, LENGTH);

            let len = test.series_block.bv.len();
            let bpp = len as f64 / LENGTH as f64;
            println!("{:<32} {:>10} {:>8.3}", test.name, len, bpp/8.0);
        }
    }
}

