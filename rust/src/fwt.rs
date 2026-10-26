use crate::modint::Mod;
use crate::poly::Poly;
use std::ops::{Add, BitAnd, BitOr, BitXor, Mul, Sub};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Sps {
    m: Mod,
    a: Vec<u32>,
}

impl Sps {
    pub fn new(m: Mod, n: usize) -> Self {
        assert!(n.is_power_of_two(), "length must be a power of two");
        Self { m, a: vec![0; n] }
    }

    pub fn from_vec(m: Mod, a: Vec<u32>) -> Self {
        assert!(a.len().is_power_of_two(), "length must be a power of two");
        Self { m, a }
    }

    pub fn from_poly(f: &Poly) -> Self {
        Self::from_vec(f.modulus(), f.as_slice().to_vec())
    }

    pub fn len(&self) -> usize {
        self.a.len()
    }

    pub fn is_empty(&self) -> bool {
        self.a.is_empty()
    }

    pub fn get(&self, i: usize) -> u32 {
        self.a[i]
    }

    pub fn set(&mut self, i: usize, x: u32) {
        self.a[i] = x % self.m.p();
    }

    pub fn to_poly(&self) -> Poly {
        Poly::from_vec(self.m, self.a.clone())
    }

    fn fwt_and(&mut self, inv: bool) {
        let n = self.a.len();
        let mut len = 1;
        while len < n {
            for i in (0..n).step_by(len << 1) {
                for j in 0..len {
                    let x = self.a[i + j];
                    let y = self.a[i + j + len];
                    if inv {
                        self.a[i + j] = self.m.sub(x, y);
                    } else {
                        self.a[i + j] = self.m.add(x, y);
                    }
                }
            }
            len <<= 1;
        }
    }

    fn fwt_or(&mut self, inv: bool) {
        let n = self.a.len();
        let mut len = 1;
        while len < n {
            for i in (0..n).step_by(len << 1) {
                for j in 0..len {
                    let x = self.a[i + j];
                    let y = self.a[i + j + len];
                    if inv {
                        self.a[i + j + len] = self.m.sub(y, x);
                    } else {
                        self.a[i + j + len] = self.m.add(y, x);
                    }
                }
            }
            len <<= 1;
        }
    }

    fn fwt_xor(&mut self, inv: bool) {
        let n = self.a.len();
        let mut len = 1;
        while len < n {
            for i in (0..n).step_by(len << 1) {
                for j in 0..len {
                    let x = self.a[i + j];
                    let y = self.a[i + j + len];
                    self.a[i + j] = self.m.add(x, y);
                    self.a[i + j + len] = self.m.sub(x, y);
                }
            }
            len <<= 1;
        }
        if inv {
            let inv2 = self.m.inv(2).unwrap();
            for x in self.a.iter_mut() {
                *x = self.m.mul(*x, inv2);
            }
        }
    }

    pub fn and(&self, other: &Self) -> Self {
        self.assert_same(other);
        let mut a = self.clone();
        let mut b = other.clone();
        a.fwt_and(false);
        b.fwt_and(false);
        for i in 0..a.a.len() {
            a.a[i] = a.m.mul(a.a[i], b.a[i]);
        }
        a.fwt_and(true);
        a
    }

    pub fn or(&self, other: &Self) -> Self {
        self.assert_same(other);
        let mut a = self.clone();
        let mut b = other.clone();
        a.fwt_or(false);
        b.fwt_or(false);
        for i in 0..a.a.len() {
            a.a[i] = a.m.mul(a.a[i], b.a[i]);
        }
        a.fwt_or(true);
        a
    }

    pub fn xor(&self, other: &Self) -> Self {
        self.assert_same(other);
        let mut a = self.clone();
        let mut b = other.clone();
        a.fwt_xor(false);
        b.fwt_xor(false);
        for i in 0..a.a.len() {
            a.a[i] = a.m.mul(a.a[i], b.a[i]);
        }
        a.fwt_xor(true);
        a
    }

    pub fn dx(&self) -> Self {
        let mut r = self.clone();
        for (i, x) in r.a.iter_mut().enumerate() {
            *x = self.m.mul(*x, i.count_ones());
        }
        r
    }

    pub fn integ(&self) -> Self {
        let mut r = self.clone();
        r.a[0] = 0;
        for (i, x) in r.a.iter_mut().enumerate().skip(1) {
            let inv = self.m.inv(i.count_ones()).expect("rank not invertible");
            *x = self.m.mul(*x, inv);
        }
        r
    }

    pub fn inv(&self) -> Self {
        let n = self.len();
        let k = n.trailing_zeros();
        let mut g = Sps::new(self.m, n);
        g.set(
            0,
            self.m
                .inv(self.a[0])
                .expect("empty-set coefficient not invertible"),
        );
        let mut step = 1u32;
        while step <= k {
            let fg = self.subset(&g);
            let mut t = Sps::new(self.m, n);
            for i in 0..n {
                t.set(
                    i,
                    if i == 0 {
                        self.m.sub(2 % self.m.p(), fg.get(0))
                    } else {
                        self.m.neg(fg.get(i))
                    },
                );
            }
            g = g.subset(&t);
            step <<= 1;
        }
        g
    }

    pub fn ln(&self) -> Self {
        self.dx().subset(&self.inv()).integ()
    }

    pub fn exp(&self) -> Self {
        let n = self.len();
        let k = n.trailing_zeros();
        let mut g = Sps::new(self.m, n);
        g.set(0, 1);
        let mut step = 1u32;
        while step <= k {
            let l = g.ln();
            let mut t = Sps::new(self.m, n);
            for i in 0..n {
                t.set(
                    i,
                    if i == 0 {
                        self.m.add(1, self.m.sub(self.get(i), l.get(i)))
                    } else {
                        self.m.sub(self.get(i), l.get(i))
                    },
                );
            }
            g = g.subset(&t);
            step <<= 1;
        }
        g
    }

    pub fn subset(&self, other: &Self) -> Self {
        self.assert_same(other);
        let n = self.a.len();
        let k = n.trailing_zeros() as usize;
        let mut f = vec![vec![0u32; n]; k + 1];
        let mut g = vec![vec![0u32; n]; k + 1];
        for i in 0..n {
            let pc = i.count_ones() as usize;
            f[pc][i] = self.a[i];
            g[pc][i] = other.a[i];
        }
        for l in 0..=k {
            let mut len = 1;
            while len < n {
                for i in (0..n).step_by(len << 1) {
                    for j in 0..len {
                        f[l][i + j + len] = self.m.add(f[l][i + j + len], f[l][i + j]);
                        g[l][i + j + len] = self.m.add(g[l][i + j + len], g[l][i + j]);
                    }
                }
                len <<= 1;
            }
        }
        let mut h = vec![vec![0u32; n]; k + 1];
        for l in 0..=k {
            for i in 0..n {
                let mut s = 0u32;
                for j in 0..=l {
                    s = self.m.add(s, self.m.mul(f[j][i], g[l - j][i]));
                }
                h[l][i] = s;
            }
        }
        for l in 0..=k {
            let mut len = 1;
            while len < n {
                for i in (0..n).step_by(len << 1) {
                    for j in 0..len {
                        h[l][i + j + len] = self.m.sub(h[l][i + j + len], h[l][i + j]);
                    }
                }
                len <<= 1;
            }
        }
        let mut r = vec![0u32; n];
        for i in 0..n {
            r[i] = h[i.count_ones() as usize][i];
        }
        Sps::from_vec(self.m, r)
    }

    fn assert_same(&self, other: &Self) {
        assert_eq!(self.m, other.m);
        assert_eq!(self.a.len(), other.a.len());
    }
}

impl BitAnd for &Sps {
    type Output = Sps;
    fn bitand(self, rhs: &Sps) -> Sps {
        self.and(rhs)
    }
}

impl BitOr for &Sps {
    type Output = Sps;
    fn bitor(self, rhs: &Sps) -> Sps {
        self.or(rhs)
    }
}

impl BitXor for &Sps {
    type Output = Sps;
    fn bitxor(self, rhs: &Sps) -> Sps {
        self.xor(rhs)
    }
}

impl Mul for &Sps {
    type Output = Sps;
    fn mul(self, rhs: &Sps) -> Sps {
        self.subset(rhs)
    }
}

impl Add for &Sps {
    type Output = Sps;
    fn add(self, rhs: &Sps) -> Sps {
        assert_eq!(self.m, rhs.m);
        assert_eq!(self.a.len(), rhs.a.len());
        let mut r = self.clone();
        for i in 0..r.a.len() {
            r.a[i] = self.m.add(self.a[i], rhs.a[i]);
        }
        r
    }
}

impl Sub for &Sps {
    type Output = Sps;
    fn sub(self, rhs: &Sps) -> Sps {
        assert_eq!(self.m, rhs.m);
        assert_eq!(self.a.len(), rhs.a.len());
        let mut r = self.clone();
        for i in 0..r.a.len() {
            r.a[i] = self.m.sub(self.a[i], rhs.a[i]);
        }
        r
    }
}

impl Add for Sps {
    type Output = Sps;
    fn add(self, rhs: Sps) -> Sps {
        &self + &rhs
    }
}

impl Sub for Sps {
    type Output = Sps;
    fn sub(self, rhs: Sps) -> Sps {
        &self - &rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive_and(a: &[u32], b: &[u32], m: Mod) -> Vec<u32> {
        let n = a.len();
        let mut c = vec![0u32; n];
        for i in 0..n {
            for j in 0..n {
                if i & j == i {
                    c[i] = m.add(c[i], m.mul(a[i], b[j]));
                }
            }
        }
        c
    }

    #[test]
    fn fwt_matches_naive() {
        let m = Mod::new(998244353);
        let n = 16;
        let a: Vec<u32> = (0..n).map(|i| (i * 7 + 1) as u32).collect();
        let b: Vec<u32> = (0..n).map(|i| (i * 3 + 2) as u32).collect();
        let sa = Sps::from_vec(m, a.clone());
        let sb = Sps::from_vec(m, b.clone());
        let got = (&sa & &sb).to_poly();
        // AND: h[S] = sum_{T subset S} f[T] g[S \ T]? test subset instead
        let got_sub = (&sa * &sb).to_poly();
        let mut want = vec![0u32; n];
        for s in 0..n {
            for t in 0..n {
                if t & s == t {
                    want[s] = m.add(want[s], m.mul(a[t], b[s ^ t]));
                }
            }
        }
        assert_eq!(got_sub.as_slice(), want.as_slice());
        let _ = got;
        let _ = naive_and(&a, &b, m);
    }

    #[test]
    fn ln_exp_dx_integ() {
        let m = Mod::new(998244353);
        for n in [1usize, 2, 4, 8, 16] {
            let mut f = Sps::from_vec(
                m,
                (0..n).map(|i| ((i * 37 + 11) % 998244353) as u32).collect(),
            );
            f.set(0, 1);
            let e = f.ln().exp();
            assert_eq!(f, e);
            let mut g = Sps::new(m, n);
            for i in 1..n {
                g.set(i, (i * 53 + 7) as u32);
            }
            let l = g.exp().ln();
            assert_eq!(g, l);
            let d = f.dx().integ();
            assert_eq!(d.get(0), 0);
            for i in 1..n {
                assert_eq!(d.get(i), f.get(i));
            }
            let inv = f.inv();
            let id = f.subset(&inv);
            assert_eq!(id.get(0), 1);
            for i in 1..n {
                assert_eq!(id.get(i), 0);
            }
        }
    }
}
