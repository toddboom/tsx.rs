extern crate bit_vec;

//use std::collections::BTreeMap;
use bit_vec::BitVec;

struct Series {
    t0: u64,        // beginning timestamp of this block (
    vprev: u64,     // previous value, stored bitwise
    lz: u8,         // current number of leading zero nibbles (4-bits)
    bv: BitVec,     // the stream of bits
    finished: bool, // is this block complete?
}

impl Series {
    pub fn new(t0: u64) -> Series {
        Series {
            t0: t0,
            vprev: 0,
            lz: 0xff,
            bv: BitVec::new(),
            finished: false,
        }
    }

    pub fn append(&mut self, v: f64) {
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
        if bits == 0 {
            self.append_bits(0b00, 2);
            return;
        }

        let lz = xor.leading_zeros() as u8;
        let zero_nibbles = lz % 4;
        let remaining_bits = 64 - (4 * zero_nibbles);

        // if we haven't recorded the leading zeros yet, do that now
        if self.lz == 0xff {
            self.lz = zero_nibbles;
        }

        // control code 10 = continue using current lz
        self.append_bits(0b10, 2);
        self.append_bits(bits, remaining_bits);

        // control code 01 = null value
        // TODO: this needs to get handled outside of this function

        // control code 11 = reset lz
        // TODO: decide if it's efficient to reset

        return;
    }

    fn append_bits(&mut self, value: u64, num_bits: u8) {
        // TODO: bail out if n > 64, otherwise we'll overflow
        for x in num_bits-1..0 {
            self.bv.push(value & (2 << x) > 0);
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
        let s = Series::new(0
        );

        assert_eq!(s.t0, 0u64);
    }

    #[test]
    fn test_push_series() {
        let mut s = Series::new(0);

        s.append(3.14);
        assert_eq!(s.vprev, bits(3.14));
    }
}

