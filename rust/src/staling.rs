use crate::fps::{exp, ln};
use crate::modint::Mod;
use crate::poly::Poly;
use crate::pre::factorials;

fn build(md: Mod, l: usize, r: usize) -> Poly {
    if l == r {
        return Poly::from_vec(md, vec![l as u32, 1]);
    }
    let mid = (l + r) >> 1;
    &build(md, l, mid) * &build(md, mid + 1, r)
}

pub fn row(md: Mod, n: usize) -> Poly {
    if n == 0 {
        return Poly::one(md);
    }
    build(md, 0, n - 1)
}

pub fn column(md: Mod, k: usize, n: usize) -> Poly {
    if k > n {
        return Poly::new(md, n + 1);
    }
    if k == 0 {
        let mut r = Poly::new(md, n + 1);
        r.set(0, 1);
        return r;
    }
    let (fac, ifac) = factorials(md, n);
    let mut f = Poly::new(md, n + 1);
    for i in 1..=n {
        f.set(i, md.inv(i as u32).expect("i not invertible"));
    }
    let shifted = &f >> 1;
    let l = ln(&shifted, n);
    let lk = &l * (k % md.p() as usize) as u32;
    let mut e = exp(&lk, n);
    e = e << k;
    e.resize(n + 1);
    for i in 0..=n {
        e.set(i, md.mul(md.mul(e.coeff(i), ifac[k]), fac[i]));
    }
    e
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_kind_small() {
        let md = Mod::new(998244353);
        let n = 12;
        let mut c = vec![vec![0u32; n + 1]; n + 1];
        c[0][0] = 1;
        for i in 1..=n {
            for j in 1..=i {
                c[i][j] = md.add(c[i - 1][j - 1], md.mul((i - 1) as u32, c[i - 1][j]));
            }
        }
        let r = row(md, n);
        for k in 0..=n {
            assert_eq!(r.coeff(k), c[n][k], "row k={k}");
        }
        for k in 1..=6 {
            let col = column(md, k, n);
            for i in 0..=n {
                assert_eq!(col.coeff(i), c[i][k], "col k={k} n={i}");
            }
        }
    }
}
