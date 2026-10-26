use crate::modint::Mod;
use std::ops::{Add, Mul, Sub};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Matrix {
    md: Mod,
    rows: usize,
    cols: usize,
    a: Vec<Vec<u32>>,
}

impl Matrix {
    pub fn new(md: Mod, rows: usize, cols: usize) -> Self {
        Self {
            md,
            rows,
            cols,
            a: vec![vec![0; cols]; rows],
        }
    }

    pub fn identity(md: Mod, n: usize) -> Self {
        let mut r = Self::new(md, n, n);
        for i in 0..n {
            r.a[i][i] = 1;
        }
        r
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn get(&self, i: usize, j: usize) -> u32 {
        self.a[i][j]
    }

    pub fn set(&mut self, i: usize, j: usize, x: u32) {
        self.a[i][j] = x % self.md.p();
    }

    pub fn transpose(&self) -> Self {
        let mut r = Self::new(self.md, self.cols, self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                r.a[j][i] = self.a[i][j];
            }
        }
        r
    }

    pub fn gauss(&mut self) -> bool {
        let mut rank = 0usize;
        for col in 0..self.cols {
            let mut piv = None;
            for i in rank..self.rows {
                if self.a[i][col] != 0 && self.md.inv(self.a[i][col]).is_some() {
                    piv = Some(i);
                    break;
                }
            }
            let Some(p) = piv else { continue };
            self.a.swap(rank, p);
            let inv = self.md.inv(self.a[rank][col]).unwrap();
            for j in col..self.cols {
                self.a[rank][j] = self.md.mul(self.a[rank][j], inv);
            }
            for i in 0..self.rows {
                if i != rank && self.a[i][col] != 0 {
                    let f = self.a[i][col];
                    for j in col..self.cols {
                        self.a[i][j] = self.md.sub(self.a[i][j], self.md.mul(f, self.a[rank][j]));
                    }
                }
            }
            rank += 1;
            if rank == self.rows {
                break;
            }
        }
        rank == self.rows
    }

    pub fn inv(&self) -> Option<Self> {
        if self.rows != self.cols {
            return None;
        }
        let n = self.rows;
        let mut aug = Self::new(self.md, n, n * 2);
        for i in 0..n {
            for j in 0..n {
                aug.a[i][j] = self.a[i][j];
            }
            aug.a[i][n + i] = 1;
        }
        let mut rank = 0usize;
        for col in 0..n {
            let mut piv = None;
            for i in rank..n {
                if aug.a[i][col] != 0 && self.md.inv(aug.a[i][col]).is_some() {
                    piv = Some(i);
                    break;
                }
            }
            let p = piv?;
            aug.a.swap(rank, p);
            let inv = self.md.inv(aug.a[rank][col]).unwrap();
            for j in 0..(n * 2) {
                aug.a[rank][j] = self.md.mul(aug.a[rank][j], inv);
            }
            for i in 0..n {
                if i != rank && aug.a[i][col] != 0 {
                    let f = aug.a[i][col];
                    for j in 0..(n * 2) {
                        aug.a[i][j] = self.md.sub(aug.a[i][j], self.md.mul(f, aug.a[rank][j]));
                    }
                }
            }
            rank += 1;
        }
        let mut res = Self::new(self.md, n, n);
        for i in 0..n {
            for j in 0..n {
                res.a[i][j] = aug.a[i][n + j];
            }
        }
        Some(res)
    }

    pub fn det(&self) -> Option<u32> {
        if self.rows != self.cols {
            return None;
        }
        let n = self.rows;
        let mut a = self.a.clone();
        let mut det = 1u32;
        for col in 0..n {
            let mut piv = None;
            for i in col..n {
                if a[i][col] != 0 && self.md.inv(a[i][col]).is_some() {
                    piv = Some(i);
                    break;
                }
            }
            let p = match piv {
                Some(p) => p,
                None => return Some(0),
            };
            if p != col {
                a.swap(p, col);
                det = self.md.neg(det);
            }
            det = self.md.mul(det, a[col][col]);
            let inv = self.md.inv(a[col][col]).unwrap();
            for i in (col + 1)..n {
                if a[i][col] != 0 {
                    let f = self.md.mul(a[i][col], inv);
                    for j in col..n {
                        a[i][j] = self.md.sub(a[i][j], self.md.mul(f, a[col][j]));
                    }
                }
            }
        }
        Some(det)
    }

    pub fn pow(&self, mut k: u64) -> Self {
        assert_eq!(self.rows, self.cols);
        let mut r = Self::identity(self.md, self.rows);
        let mut b = self.clone();
        while k != 0 {
            if k & 1 != 0 {
                r = &r * &b;
            }
            b = &b * &b;
            k >>= 1;
        }
        r
    }
}

impl Add for &Matrix {
    type Output = Matrix;
    fn add(self, rhs: &Matrix) -> Matrix {
        assert_eq!(
            (self.rows, self.cols, self.md),
            (rhs.rows, rhs.cols, rhs.md)
        );
        let mut r = self.clone();
        for i in 0..self.rows {
            for j in 0..self.cols {
                r.a[i][j] = self.md.add(self.a[i][j], rhs.a[i][j]);
            }
        }
        r
    }
}

impl Sub for &Matrix {
    type Output = Matrix;
    fn sub(self, rhs: &Matrix) -> Matrix {
        assert_eq!(
            (self.rows, self.cols, self.md),
            (rhs.rows, rhs.cols, rhs.md)
        );
        let mut r = self.clone();
        for i in 0..self.rows {
            for j in 0..self.cols {
                r.a[i][j] = self.md.sub(self.a[i][j], rhs.a[i][j]);
            }
        }
        r
    }
}

impl Mul for &Matrix {
    type Output = Matrix;
    fn mul(self, rhs: &Matrix) -> Matrix {
        assert_eq!(self.cols, rhs.rows);
        assert_eq!(self.md, rhs.md);
        let mut r = Matrix::new(self.md, self.rows, rhs.cols);
        for i in 0..self.rows {
            for k in 0..self.cols {
                if self.a[i][k] == 0 {
                    continue;
                }
                for j in 0..rhs.cols {
                    if rhs.a[k][j] == 0 {
                        continue;
                    }
                    r.a[i][j] = self
                        .md
                        .add(r.a[i][j], self.md.mul(self.a[i][k], rhs.a[k][j]));
                }
            }
        }
        r
    }
}

impl Add for Matrix {
    type Output = Matrix;
    fn add(self, rhs: Matrix) -> Matrix {
        &self + &rhs
    }
}

impl Sub for Matrix {
    type Output = Matrix;
    fn sub(self, rhs: Matrix) -> Matrix {
        &self - &rhs
    }
}

impl Mul for Matrix {
    type Output = Matrix;
    fn mul(self, rhs: Matrix) -> Matrix {
        &self * &rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inverse() {
        let md = Mod::new(998244353);
        let mut a = Matrix::new(md, 2, 2);
        a.set(0, 0, 1);
        a.set(0, 1, 2);
        a.set(1, 0, 3);
        a.set(1, 1, 4);
        let inv = a.inv().unwrap();
        let id = &a * &inv;
        assert_eq!(id.get(0, 0), 1);
        assert_eq!(id.get(0, 1), 0);
        assert_eq!(id.get(1, 0), 0);
        assert_eq!(id.get(1, 1), 1);
        assert_eq!(a.det(), Some(md.sub(0, 2)));
    }
}
