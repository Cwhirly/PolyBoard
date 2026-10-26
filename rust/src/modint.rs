#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Mod {
    p: u32,
}

impl Mod {
    pub fn new(p: u32) -> Self {
        assert!(p >= 2, "modulus must be at least 2");
        Self { p }
    }

    pub const fn p(self) -> u32 {
        self.p
    }

    pub fn add(self, a: u32, b: u32) -> u32 {
        let s = a + b;
        if s >= self.p {
            s - self.p
        } else {
            s
        }
    }

    pub fn sub(self, a: u32, b: u32) -> u32 {
        if a >= b {
            a - b
        } else {
            a + self.p - b
        }
    }

    pub fn neg(self, a: u32) -> u32 {
        if a == 0 {
            0
        } else {
            self.p - a
        }
    }

    pub fn mul(self, a: u32, b: u32) -> u32 {
        (a as u64 * b as u64 % self.p as u64) as u32
    }

    pub fn pow(self, mut a: u32, mut e: u64) -> u32 {
        let mut r = 1u32;
        a %= self.p;
        while e != 0 {
            if e & 1 != 0 {
                r = self.mul(r, a);
            }
            a = self.mul(a, a);
            e >>= 1;
        }
        r
    }

    pub fn inv(self, a: u32) -> Option<u32> {
        let a = a % self.p;
        let (mut r0, mut r1) = (a as i64, self.p as i64);
        let (mut s0, mut s1) = (1i64, 0i64);
        while r1 != 0 {
            let q = r0 / r1;
            let r2 = r0 - q * r1;
            r0 = r1;
            r1 = r2;
            let s2 = s0 - q * s1;
            s0 = s1;
            s1 = s2;
        }
        if r0 != 1 {
            None
        } else {
            Some(((s0 % self.p as i64 + self.p as i64) % self.p as i64) as u32)
        }
    }

    pub fn sqrt(self, a: u32) -> Option<u32> {
        let p = self.p;
        let a = a % p;
        if a == 0 {
            return Some(0);
        }
        if p == 2 {
            return Some(a);
        }
        if self.pow(a, (p - 1) as u64 / 2) != 1 {
            return None;
        }
        if p % 4 == 3 {
            let r = self.pow(a, (p as u64 + 1) / 4);
            return Some(r.min(p - r));
        }
        let mut q = p - 1;
        let mut s = 0u32;
        while q & 1 == 0 {
            q >>= 1;
            s += 1;
        }
        let mut z = 2u32;
        while self.pow(z, (p - 1) as u64 / 2) != p - 1 {
            z += 1;
        }
        let mut c = self.pow(z, q as u64);
        let mut x = self.pow(a, (q as u64 + 1) / 2);
        let mut t = self.pow(a, q as u64);
        let mut m = s;
        while t != 1 {
            let mut i = 0u32;
            let mut tt = t;
            while tt != 1 {
                tt = self.mul(tt, tt);
                i += 1;
            }
            let b = self.pow(c, 1u64 << (m - i - 1));
            let b2 = self.mul(b, b);
            x = self.mul(x, b);
            t = self.mul(t, b2);
            c = b2;
            m = i;
        }
        Some(x.min(p - x))
    }

    pub fn from_i64(self, x: i64) -> u32 {
        ((x % self.p as i64 + self.p as i64) % self.p as i64) as u32
    }
}
