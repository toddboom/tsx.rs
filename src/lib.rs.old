extern crate bit_vec;

use std::collections::BTreeMap;
use bit_vec::BitVec;

struct Series {
    t0: u64,        // beginning timestamp of this block
    tn: u64,        // latest timestamp in this block
    vn: f64,        // latest value in this block
    td: u64,        // the timestamp delta?
    bv: BitVec,     // the stream of bits

    leading: u8,
    trailing: u8,
    finished: bool,
}

impl Series {
    pub fn new(t0: u64) -> Series {
        let mut s = Series {
            t0: 0,
            tn: 0,
            td: 0,
            vn: 0f64,
            bv: BitVec::new(),
            leading: 0xff,
            trailing: 0,
            finished: false,
        };

        s.push_timestamp(t0, 32);
        s
    }

    pub fn push(&mut self, t: u64, v: f64) {
        // is this the first point?
        if self.tn == 0 {
            self.tn = t;
            self.vn = v;
            self.td = t - self.t0;

            self.push_timestamp(t, 14);
            self.push_value(v, 64);

            return;
        }

        let td = t - self.tn;
        let dd = td - self.td;
        ;

        match dd {
            0 => {
                self.bv.push(false);
            },
            1 ... 127 => {
                self.push_timestamp(0b10, 2);
                self.push_timestamp(dd, 7);
            },
            128 ... 511 => {
                self.push_timestamp(0b110, 3);
                self.push_timestamp(dd, 9);
            },
            512 ... 4095 => {
                self.push_timestamp(0b1110, 4);
                self.push_timestamp(dd, 12
                );
            },
            _ => {
                self.push_timestamp(0b1111, 4);
                self.push_timestamp(dd, 32);
            }
        }

    }

    fn push_value(&mut self, v: f64, n: u8) {
        let bits = unsafe { std::mem::transmute::<f64,u64>(v) };

        // TODO: bail out if n > 64, otherwise we'll overflow
        for x in n-1..0 {
            self.bv.push(bits & (2 << n) > 0);
        }
    }

    fn push_timestamp(&mut self, v: u64, n: u8) {
        for x in n-1..0 {
            self.bv.push(v & (2 << n) > 0);
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

