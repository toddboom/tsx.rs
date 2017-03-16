extern crate bit_vec;

use std::collections::BTreeMap;
use bit_vec::BitVec;

struct Series {
    t0: u64,        // beginning timestamp of this block (
    vprev: u64,     // previous value, stored bitwise
    bv: BitVec,     // the stream of bits
    finished: bool, // is this block complete?
}

impl Series {
    pub fn new(t0: u64) -> Series {
        Series {
            t0: t0,
            vprev: 0,
            bv: BitVec::new(),
            finished: false,
        }
    }

    pub fn append(&mut self, v: f64) {
        // we need a bitwise representation of our float
        let bits = unsafe { std::mem::transmute::<f64,u64>(v) };

        // if this is the first point, write the entire value
        if bv.len() == 0 {
            self.append_bits(v, 64);

            return;
        }

        let xor = self.vprev ^ bits;

        self.vprev = bits;
    }

    fn append_bits(&mut self, v: u64, n: u8) {
        let bits = unsafe { std::mem::transmute::<f64,u64>(v) };

        // TODO: bail out if n > 64, otherwise we'll overflow
        for x in n-1..0 {
            self.bv.push(bits & (2 << n) > 0);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_series() {
        let s = Series::new(0
        );

        assert_eq!(s.t0, 0u64);
    }

    #[test]
    fn test_push_series() {
        let mut s = Series::new(0);

        s.push(999, 3.14);
        assert_eq!(s.tn, 999);
        assert_eq!(s.vn, 3.14);
    }
}

