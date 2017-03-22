extern crate bit_vec;

//use std::collections::BTreeMap;
use bit_vec::BitVec;

struct Series {
    t0: u64,        // beginning timestamp of this block (
    vprev: u64,     // previous value, stored bitwise
    lz: u8,         // current number of leading zero nibbles (4-bits)
    bv: BitVec,     // the stream of bits
    finished: bool, // is this block complete?
    length: usize   // total number of values written to this block
}

struct SeriesIterator<'a> {
    series: &'a Series,
    prev_bits: u64,
    prev_nibbles: u8,
    index: usize
}

impl Series {
    pub fn new(t0: u64) -> Series {
        Series {
            t0: t0,
            vprev: 0,
            lz: 0xff,
            bv: BitVec::new(),
            finished: false,
            length: 0,
        }
    }

    pub fn append(&mut self, v: f64) {
        self.length += 1;

        // we need a bitwise representation of our float
        let bits = unsafe { std::mem::transmute::<f64,u64>(v) };

        // if this is the first point, write the entire value
        if self.bv.len() == 0 {
            self.vprev = bits;
            self.append_bits(bits, 64);
            return;
        }

        let xor = self.vprev ^ bits;
        self.vprev = bits;

        // control code 00 = repeated point
        if xor == 0 {
            self.append_bits(0b00, 2);
            return;
        }

        let lz = xor.leading_zeros() as u8;
        let zero_nibbles = lz.wrapping_div(4);
        let remaining_bits = 64 - (4 * zero_nibbles);

        // if we haven't recorded the leading zeros yet, do that now
        if (self.lz != zero_nibbles) || (self.lz == 0xff) {
            // control code 11 = reset lz
            // TODO: decide if it's efficient to reset
            // maybe just if lz > s.lz ?
            self.lz = zero_nibbles;
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

    fn extract(&self, index: &mut usize, previous: u64, previous_nibbles: u8) -> f64 {
        let mut value: u64 = 0;
        if *index == 0 {
            value = self.get_bits(index, 64);
            println!("value: {:015x}", value);
        } else {
            let control_code = self.get_bits(index, 2);
            match control_code {
                0b00 => {
                    value = previous;
                }
                0b01 => {

                }
                0b10 => {
                    value = self.get_bits(index, 64 - 4 * previous_nibbles);
                    println!("value: {:015x}", value);
                    value ^= previous;
                    println!("value: {:015x}", value);
                }
                0b11 => {
                    let nibbles = self.get_bits(index, 4) as u8;
                    println!("nibbles: {:?}", nibbles);
                    let remaining_bits = 64 - 4 * nibbles;
                    println!("rem bits: {:?}", remaining_bits);
                    value = self.get_bits(index, 64 - 4 * nibbles);
                    println!("value: {:015x}", value);
                    value ^= previous;
                    println!("value: {:015x}", value);
                }
                _ => {
                    // impossible!
                }
            }
        }

        unsafe { std::mem::transmute::<u64, f64>(value) }
    }

    fn extract_bits(&self, index: &mut usize, prev_nibbles: &mut u8) -> Option<u64> {
        let mut xor: u64 = 0;
        println!("vec length: {:?}, {:?}", self.bv.len(), index);
        if *index == 0 {
            xor = self.get_bits(index, 64);
        } else if *index >= self.bv.len() {
            return None;
        } else {
            let control_code = self.get_bits(index, 2);
            match control_code {
                0b00 => {
                    xor = 0;
                }
                0b01 => {
                    // TODO: figure out how to handle null values
                }
                0b10 => {
                    xor = self.get_bits(index, 64 - 4 * *prev_nibbles);
                    println!("value: {:015x}", xor);
                }
                0b11 => {
                    let nibbles = self.get_bits(index, 4) as u8;
                    println!("nibbles: {:?}", nibbles);
                    let remaining_bits = 64 - 4 * nibbles;
                    println!("rem bits: {:?}", remaining_bits);
                    xor = self.get_bits(index, 64 - 4 * nibbles);
                    println!("value: {:015x}", xor);
                    *prev_nibbles = nibbles;
                }
                _ => {
                    // this should theoretically be impossible
                    return None;
                }
            }
        }

        return Some(xor);
    }

    fn iter(&self) -> SeriesIterator {
        SeriesIterator {
            series: self,
            index: 0,
            prev_bits: 0,
            prev_nibbles: 0
        }
    }
}

impl<'a> Iterator for SeriesIterator<'a> {
    type Item = f64;

    fn next(&mut self) -> Option<f64> {
        let value = self.series.extract_bits(&mut self.index, &mut self.prev_nibbles);
        match value {
            None => None,
            Some(x) => {
                let bits = x ^ self.prev_bits;
                self.prev_bits = bits;
                return Some(unsafe { std::mem::transmute::<u64,f64>(bits) });
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bits(v: f64) -> u64 {
        unsafe { std::mem::transmute::<f64,u64>(v) }
    }

    #[test]
    fn test_new_series() {
        let s = Series::new(10000);

        // ensure that everything was initialized correctly
        assert_eq!(s.t0, 10000u64);
        assert_eq!(s.lz, 0xff);
        assert_eq!(s.bv.len(), 0);
        assert_eq!(s.length, 0);
        assert_eq!(s.vprev, 0);
        assert_eq!(s.finished, false);
    }

    #[test]
    fn test_append_series() {
        let mut s = Series::new(0);

        s.append(3.14); // 0x40091eb851eb851f
        assert_eq!(s.vprev, bits(3.14));
        assert_eq!(s.length, 1);
        assert_eq!(s.bv.len(), 64);
        assert_eq!(s.lz, 0xff);

        s.append(3.15); // 0x4009333333333333 xor with previous = 18 lz = ~4 nibbles
        assert_eq!(s.vprev, bits(3.15));
        assert_eq!(s.length, 2);
        assert_eq!(s.bv.len(), 118); // previous + 2 + 4 + 48
        assert_eq!(s.lz, 4);

        s.append(3.15); // repeat value
        assert_eq!(s.vprev, bits(3.15));
        assert_eq!(s.length, 3);
        assert_eq!(s.bv.len(), 120); // previous + 2
        assert_eq!(s.lz, 4);

        s.append(3.1501); // 0x40093367a0f9096c xor with previous = 25 lz = ~6 nibbles
        assert_eq!(s.vprev, bits(3.1501));
        assert_eq!(s.length, 4);
        assert_eq!(s.bv.len(), 166); // previous + 2 + 4 + 40
        assert_eq!(s.lz, 6);

        s.append(3.1502); // 0x4009339c0ebedfa4 xor with previous = 24 lz = ~6 nibbles
        assert_eq!(s.vprev, bits(3.1502));
        assert_eq!(s.length, 5);
        assert_eq!(s.bv.len(), 208); // previous + 2 + 40
        assert_eq!(s.lz, 6);

        s.append(4.20); // 0x4010cccccccccccd xor with previous = 11 lz = ~2 nibbles
        assert_eq!(s.vprev, bits(4.20));
        assert_eq!(s.length, 6);
        assert_eq!(s.bv.len(), 270); // previous + 2 + 4 + 56
        assert_eq!(s.lz, 2);
    }

    #[test]
    fn test_extract_series() {
        let mut s = Series::new(0);
        let mut index = 0;

        s.append(3.14); // 0x40091eb851eb851f
        s.append(3.15); // 0x4009333333333333
        s.append(3.15); // 0x4009333333333333
        s.append(3.1501); // 0x40093367a0f9096c
        s.append(3.1502); // 0x4009339c0ebedfa4
        s.append(4.20); // 0x4010cccccccccccd

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

    #[test]
    fn test_series_iterator() {
        let mut s = Series::new(0);
        let mut index = 0;

        s.append(3.14); // 0x40091eb851eb851f
        s.append(3.15); // 0x4009333333333333
        s.append(3.15); // 0x4009333333333333
        s.append(3.1501); // 0x40093367a0f9096c
        s.append(3.1502); // 0x4009339c0ebedfa4
        s.append(4.20); // 0x4010cccccccccccd

        let mut iter = s.iter();
        assert_eq!(iter.next(), Some(3.14));
        assert_eq!(iter.next(), Some(3.15));
        assert_eq!(iter.next(), Some(3.15));
        assert_eq!(iter.next(), Some(3.1501));
        assert_eq!(iter.next(), Some(3.1502));
        assert_eq!(iter.next(), Some(4.20));
        assert_eq!(iter.next(), None);
    }
}

