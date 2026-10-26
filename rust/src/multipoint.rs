use crate::modint::Mod;
use crate::poly::Poly;

pub struct Multipoint {
    xs: Vec<u32>,
    tree: Vec<Poly>,
}

impl Multipoint {
    pub fn new(m: Mod, xs: Vec<u32>) -> Self {
        let n = xs.len();
        let mut tree = vec![Poly::zero(m); (n << 2) + 5];
        if n != 0 {
            Self::build(m, &xs, 1, 0, n - 1, &mut tree);
        }
        Self { xs, tree }
    }

    fn build(m: Mod, xs: &[u32], p: usize, l: usize, r: usize, tree: &mut [Poly]) {
        if l == r {
            tree[p] = Poly::from_vec(m, vec![m.neg(xs[l]), 1]);
            return;
        }
        let mid = (l + r) >> 1;
        Self::build(m, xs, p << 1, l, mid, tree);
        Self::build(m, xs, p << 1 | 1, mid + 1, r, tree);
        tree[p] = &tree[p << 1] * &tree[p << 1 | 1];
    }

    pub fn eval(&self, f: &Poly) -> Vec<u32> {
        if self.xs.is_empty() {
            return Vec::new();
        }
        let mut r = vec![0u32; self.xs.len()];
        let mut f = f.clone();
        if f.len() >= self.tree[1].len() {
            f = f.div_rem(&self.tree[1]).1;
        }
        self.eval_at(1, 0, self.xs.len() - 1, &f, &mut r);
        r
    }

    fn eval_at(&self, p: usize, l: usize, r: usize, f: &Poly, out: &mut [u32]) {
        if l == r {
            out[l] = f.coeff(0);
            return;
        }
        let mid = (l + r) >> 1;
        let fl = f.div_rem(&self.tree[p << 1]).1;
        let fr = f.div_rem(&self.tree[p << 1 | 1]).1;
        self.eval_at(p << 1, l, mid, &fl, out);
        self.eval_at(p << 1 | 1, mid + 1, r, &fr, out);
    }
}

pub struct Interpolation {
    m: Mod,
    xs: Vec<u32>,
    tree: Vec<Poly>,
    w: Vec<u32>,
}

impl Interpolation {
    pub fn new(m: Mod, xs: Vec<u32>, ys: Vec<u32>) -> Self {
        assert_eq!(xs.len(), ys.len());
        let mp = Multipoint::new(m, xs.clone());
        let n = xs.len();
        let tree = mp.tree.clone();
        if n != 0 {
            let md = tree[1].derivative();
            let vals = mp.eval(&md);
            let mut w = vec![0u32; n];
            for i in 0..n {
                w[i] = m.mul(ys[i], m.inv(vals[i]).expect("derivative zero at point"));
            }
            Self { m, xs, tree, w }
        } else {
            Self {
                m,
                xs,
                tree: Vec::new(),
                w: Vec::new(),
            }
        }
    }

    pub fn interpolate(&self) -> Poly {
        if self.xs.is_empty() {
            return Poly::zero(self.m);
        }
        self.dfs(1, 0, self.xs.len() - 1)
    }

    fn dfs(&self, p: usize, l: usize, r: usize) -> Poly {
        if l == r {
            return Poly::from_vec(self.m, vec![self.w[l]]);
        }
        let mid = (l + r) >> 1;
        let left = self.dfs(p << 1, l, mid);
        let right = self.dfs(p << 1 | 1, mid + 1, r);
        &(&left * &self.tree[p << 1 | 1]) + &(&right * &self.tree[p << 1])
    }
}

pub fn pmpe(f: &Poly, xs: &[u32]) -> Vec<u32> {
    Multipoint::new(f.modulus(), xs.to_vec()).eval(f)
}

pub fn pfi(m: Mod, xs: &[u32], ys: &[u32]) -> Poly {
    Interpolation::new(m, xs.to_vec(), ys.to_vec()).interpolate()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eval_interp() {
        let m = Mod::new(998244353);
        let xs: Vec<u32> = (0..12).map(|i| (i * 7 + 3) as u32).collect();
        let f = Poly::from_vec(m, (0..8).map(|i| (i * 11 + 5) as u32).collect());
        let ys = pmpe(&f, &xs);
        for (i, &x) in xs.iter().enumerate() {
            let mut v = 0u32;
            for j in (0..f.len()).rev() {
                v = m.add(m.mul(v, x), f.coeff(j));
            }
            assert_eq!(ys[i], v);
        }
        let g = pfi(m, &xs, &ys);
        for &x in &xs {
            let mut a = 0u32;
            let mut b = 0u32;
            for j in (0..f.len()).rev() {
                a = m.add(m.mul(a, x), f.coeff(j));
            }
            for j in (0..g.len()).rev() {
                b = m.add(m.mul(b, x), g.coeff(j));
            }
            assert_eq!(a, b);
        }
    }
}
