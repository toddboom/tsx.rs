extern crate bit_vec;

use std::collections::BTreeMap;
use bit_vec::BitVec;

struct Series {
    t0: u64,        // first timestamp in this block
    tn: u64,        // latest timestamp in this block
    td: u64,        // the timestamp delta?
    vn: f64,        // latest value in this block
    bv: BitVec,     // the stream of bits

    leading: u8,
    trailing: u8,
    finished: bool,
}

impl Series {
    pub fn new() -> Series {
        Series {
            t0: 0,
            tn: 0,
            td: 0,
            vn: 0f64,
            bv: BitVec::new(),
            leading: 0,
            trailing: 0,
            finished: false,
        }
    }

    pub fn push(&mut self, t: u64, v: f64) {
        // is this the first point?
        if self.tn == 0 {
            self.tn = t.clone();
            self.vn = v.clone();
            self.td = t.clone();

            self.push_bits(v, 64);
        }
    }

    fn push_bits(&mut self, v: f64, n: u8) {
        let b64 = unsafe { std::mem::transmute::<f64, u64>(v) };
        for x in n-1..0 {
            self.bv.push(b64 & (1 << x) > 0);
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
        let mut s = Series::new();

        s.push(999, 3.14);
        assert_eq!(s.tn, 999);
        assert_eq!(s.td, 999);
        assert_eq!(s.vn, 3.14);
    }
}
