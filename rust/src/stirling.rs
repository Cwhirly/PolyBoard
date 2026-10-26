use crate::modint::Mod;
use crate::poly::Poly;
use crate::pre::factorials;

pub fn row(md: Mod, n: usize) -> Poly {
    let (_, ifac) = factorials(md, n);
    let mut f = Poly::new(md, n + 1);
    let mut g = Poly::new(md, n + 1);
    for i in 0..=n {
        f.set(i, md.mul(ifac[i], md.pow(i as u32, n as u64)));
        g.set(i, if i & 1 == 1 { md.neg(ifac[i]) } else { ifac[i] });
    }
    let mut h = &f * &g;
    h.resize(n + 1);
    h
}

fn dev_mul(l: usize, r: usize, md: Mod) -> Poly {
    if l == r {
        return Poly::from_vec(md, vec![1, md.neg(l as u32)]);
    }
    let mid = (l + r) >> 1;
    let a = dev_mul(l, mid, md);
    let b = dev_mul(mid + 1, r, md);
    &a * &b
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
    let p = dev_mul(1, k, md);
    let inv = p.inv(n + 1);
    let mut r = &inv << k;
    r.resize(n + 1);
    r
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modint::Mod;

    #[test]
    fn small_values() {
        let md = Mod::new(998244353);
        let r = row(md, 10);
        // S(10, k) recurrence check
        let mut s = vec![vec![0u32; 11]; 11];
        s[0][0] = 1;
        for n in 1..=10 {
            for k in 1..=n {
                s[n][k] = md.add(s[n - 1][k - 1], md.mul(k as u32, s[n - 1][k]));
            }
        }
        for k in 0..=10 {
            assert_eq!(r.coeff(k), s[10][k], "k={k}");
        }
        for k in 1..=5 {
            let c = column(md, k, 10);
            for n in 0..=10 {
                assert_eq!(c.coeff(n), s[n][k], "k={k} n={n}");
            }
        }
    }
}
