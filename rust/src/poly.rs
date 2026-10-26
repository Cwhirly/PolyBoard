use crate::modint::Mod;
use crate::ntt::convolution;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Shl, Shr, Sub, SubAssign};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Poly {
    m: Mod,
    a: Vec<u32>,
}

impl Poly {
    pub fn new(m: Mod, n: usize) -> Self {
        Self { m, a: vec![0; n] }
    }

    pub fn from_vec(m: Mod, a: Vec<u32>) -> Self {
        Self { m, a }
    }

    pub fn zero(m: Mod) -> Self {
        Self { m, a: Vec::new() }
    }

    pub fn one(m: Mod) -> Self {
        Self { m, a: vec![1] }
    }

    pub fn modulus(&self) -> Mod {
        self.m
    }

    pub fn len(&self) -> usize {
        self.a.len()
    }

    pub fn is_empty(&self) -> bool {
        self.a.is_empty()
    }

    pub fn coeff(&self, i: usize) -> u32 {
        self.a.get(i).copied().unwrap_or(0)
    }

    pub fn set(&mut self, i: usize, x: u32) {
        if i >= self.a.len() {
            self.a.resize(i + 1, 0);
        }
        self.a[i] = x % self.m.p();
    }

    pub fn as_slice(&self) -> &[u32] {
        &self.a
    }

    pub fn resize(&mut self, n: usize) {
        self.a.resize(n, 0);
    }

    pub fn push(&mut self, x: u32) {
        self.a.push(x % self.m.p());
    }

    pub fn truncate(&mut self, n: usize) {
        self.a.truncate(n);
    }

    pub fn rev(&mut self) {
        self.a.reverse();
    }

    pub fn tp(&mut self) {
        let p = self.m.p();
        for x in self.a.iter_mut() {
            *x %= p;
        }
    }

    pub fn trunc(&self, n: usize) -> Self {
        let mut r = self.clone();
        r.resize(n);
        r
    }

    pub fn mul_trunc(&self, other: &Self, n: usize) -> Self {
        self.assert_same(other);
        let mut c = convolution(&self.a, &other.a, self.m);
        c.resize(n, 0);
        c.truncate(n);
        Poly::from_vec(self.m, c)
    }

    pub fn derivative(&self) -> Self {
        if self.a.len() <= 1 {
            return Poly::zero(self.m);
        }
        let mut r = vec![0u32; self.a.len() - 1];
        for i in 1..self.a.len() {
            r[i - 1] = self.m.mul(self.a[i], i as u32);
        }
        Poly::from_vec(self.m, r)
    }

    pub fn integral(&self) -> Self {
        let mut r = vec![0u32; self.a.len() + 1];
        for i in 0..self.a.len() {
            r[i + 1] = self
                .m
                .mul(self.a[i], self.m.inv((i + 1) as u32).unwrap_or(0));
        }
        Poly::from_vec(self.m, r)
    }

    pub fn inv(&self, n: usize) -> Self {
        assert!(
            !self.a.is_empty() && self.a[0] != 0,
            "constant term must be invertible"
        );
        let mut g = vec![self.m.inv(self.a[0]).expect("constant term not invertible")];
        let mut m = 1usize;
        while m < n {
            let m2 = (m << 1).min(n);
            let f = self.trunc(m2);
            let mut t = f.mul_trunc(&Poly::from_vec(self.m, g.clone()), m2);
            t.resize(m2);
            t.a[0] = self.m.sub(2 % self.m.p(), t.a[0]);
            for x in t.a.iter_mut().skip(1) {
                *x = self.m.neg(*x);
            }
            let g2 = Poly::from_vec(self.m, g).mul_trunc(&t, m2);
            g = g2.a;
            m = m2;
        }
        Poly::from_vec(self.m, g)
    }

    pub fn div_rem(&self, other: &Self) -> (Self, Self) {
        self.assert_same(other);
        assert!(
            !other.a.is_empty() && other.a[0] != 0,
            "divisor must have invertible constant term"
        );
        if self.a.len() < other.a.len() {
            return (Poly::zero(self.m), self.clone());
        }
        let n = self.a.len();
        let m = other.a.len();
        let mut fr = self.clone();
        fr.rev();
        let mut gr = other.clone();
        gr.rev();
        let qr = fr.mul_trunc(&gr.inv(n - m + 1), n - m + 1);
        let mut q = qr;
        q.rev();
        q.resize(n - m + 1);
        let r = self - &(&q * other);
        (q, r.trunc(m - 1))
    }

    pub fn pow_mod_x(&self, mut k: u64, n: usize) -> Self {
        assert!(
            !self.a.is_empty() && self.a[0] != 0,
            "constant term must be invertible"
        );
        let mut r = Poly::one(self.m);
        let mut b = self.trunc(n);
        while k != 0 {
            if k & 1 != 0 {
                r = r.mul_trunc(&b, n);
            }
            b = b.mul_trunc(&b, n);
            k >>= 1;
        }
        r.resize(n);
        r
    }

    fn assert_same(&self, other: &Self) {
        assert_eq!(self.m, other.m, "moduli differ");
    }
}

impl Add for &Poly {
    type Output = Poly;
    fn add(self, rhs: &Poly) -> Poly {
        self.assert_same(rhs);
        let n = self.len().max(rhs.len());
        let mut a = vec![0u32; n];
        for i in 0..n {
            a[i] = self.m.add(self.coeff(i), rhs.coeff(i));
        }
        Poly::from_vec(self.m, a)
    }
}

impl Sub for &Poly {
    type Output = Poly;
    fn sub(self, rhs: &Poly) -> Poly {
        self.assert_same(rhs);
        let n = self.len().max(rhs.len());
        let mut a = vec![0u32; n];
        for i in 0..n {
            a[i] = self.m.sub(self.coeff(i), rhs.coeff(i));
        }
        Poly::from_vec(self.m, a)
    }
}

impl Mul for &Poly {
    type Output = Poly;
    fn mul(self, rhs: &Poly) -> Poly {
        self.assert_same(rhs);
        let a = convolution(&self.a, &rhs.a, self.m);
        Poly::from_vec(self.m, a)
    }
}

impl Neg for &Poly {
    type Output = Poly;
    fn neg(self) -> Poly {
        Poly::from_vec(self.m, self.a.iter().map(|&x| self.m.neg(x)).collect())
    }
}

impl Add for Poly {
    type Output = Poly;
    fn add(self, rhs: Poly) -> Poly {
        &self + &rhs
    }
}

impl Sub for Poly {
    type Output = Poly;
    fn sub(self, rhs: Poly) -> Poly {
        &self - &rhs
    }
}

impl Mul for Poly {
    type Output = Poly;
    fn mul(self, rhs: Poly) -> Poly {
        &self * &rhs
    }
}

impl Neg for Poly {
    type Output = Poly;
    fn neg(self) -> Poly {
        -&self
    }
}

impl AddAssign<&Poly> for Poly {
    fn add_assign(&mut self, rhs: &Poly) {
        *self = &*self + rhs;
    }
}

impl SubAssign<&Poly> for Poly {
    fn sub_assign(&mut self, rhs: &Poly) {
        *self = &*self - rhs;
    }
}

impl MulAssign<&Poly> for Poly {
    fn mul_assign(&mut self, rhs: &Poly) {
        *self = &*self * rhs;
    }
}

impl Mul<u32> for &Poly {
    type Output = Poly;
    fn mul(self, rhs: u32) -> Poly {
        Poly::from_vec(
            self.m,
            self.a
                .iter()
                .map(|&x| self.m.mul(x, rhs % self.m.p()))
                .collect(),
        )
    }
}

impl Mul<u32> for Poly {
    type Output = Poly;
    fn mul(self, rhs: u32) -> Poly {
        &self * rhs
    }
}

impl MulAssign<u32> for Poly {
    fn mul_assign(&mut self, rhs: u32) {
        *self = &*self * rhs;
    }
}

impl Shl<usize> for &Poly {
    type Output = Poly;
    fn shl(self, k: usize) -> Poly {
        let mut a = vec![0u32; k];
        a.extend_from_slice(&self.a);
        Poly::from_vec(self.m, a)
    }
}

impl Shr<usize> for &Poly {
    type Output = Poly;
    fn shr(self, k: usize) -> Poly {
        if k >= self.a.len() {
            Poly::zero(self.m)
        } else {
            Poly::from_vec(self.m, self.a[k..].to_vec())
        }
    }
}

impl Shl<usize> for Poly {
    type Output = Poly;
    fn shl(self, k: usize) -> Poly {
        &self << k
    }
}

impl Shr<usize> for Poly {
    type Output = Poly;
    fn shr(self, k: usize) -> Poly {
        &self >> k
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arithmetic() {
        let m = Mod::new(998244353);
        let f = Poly::from_vec(m, vec![1, 2, 3]);
        let g = Poly::from_vec(m, vec![4, 5]);
        assert_eq!((&f + &g).as_slice(), &[5, 7, 3]);
        assert_eq!((&f - &g).as_slice(), &[m.sub(1, 4), m.sub(2, 5), 3]);
        assert_eq!((&f * &g).as_slice(), &[4, 13, 22, 15]);
    }

    #[test]
    fn inverse_and_div() {
        let m = Mod::new(998244353);
        let f = Poly::from_vec(m, vec![1, 2, 3, 4, 5]);
        let inv = f.inv(8);
        let prod = f.mul_trunc(&inv, 8);
        assert_eq!(prod.as_slice()[0], 1);
        assert!(prod.as_slice()[1..].iter().all(|&x| x == 0));
        let g = Poly::from_vec(m, vec![1, 1, 1]);
        let (q, r) = f.div_rem(&g);
        assert_eq!(&(&q * &g) + &r, f);
    }
}
