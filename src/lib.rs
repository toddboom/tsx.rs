extern crate bit_vec;

use std::collections::BTreeMap;
use bit_vec::BitVec;

struct Series {
    t0: u64,
    tn: u64,
    vn: f64,
    bv: BitVec,

    leading: u8,
    trailing: u8,
    finished: bool,

    td: u64
}

impl Series {
    pub fn new() -> Series {
        Series {
            t0: 0,
            tn: 0,
            vn: 0f64,
            bv: BitVec::new(),
            leading: 0,
            trailing: 0,
            finished: false,
            td: 0
        }
    }

    pub fn push(&mut self, t: u64, v: f64) {

        // is this the first point?
        if self.tn == 0 {
            self.tn = t;
            self.vn = v;
            self.td = t;

            self.push_bits(t, 14);
        }
    }

    fn push_bits(&mut self, v: u64, n: u8) {
        for x in 0..n {
            self.bv.push(v && 2 << n == 1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_series() {
        let s = Series::new();

        assert_eq!(s.t0, 0u64);
    }

    #[test]
    fn test_push_series() {
        let s = Series::new();

        s.push(999, 3.14);
        assert_eq!(s.t0, 999);
        assert_eq!(s.tn, 999);
        assert_eq!(s.td, 999);
        assert_eq!(s.vn, 3.14);
    }
}
