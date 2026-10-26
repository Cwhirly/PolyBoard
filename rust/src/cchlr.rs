use crate::ffp::ffptp;
use crate::modint::Mod;
use crate::poly::Poly;
use crate::pre::factorials;

fn calc(g: &Poly, mut n: u64) -> Poly {
    let md = g.modulus();
    let mut res = Poly::one(md);
    let mut a = Poly::from_vec(md, vec![0, 1]);
    while n != 0 {
        if n & 1 != 0 {
            res = &res * &a;
            res = res.div_rem(g).1;
            res.resize(g.len());
        }
        a = &a * &a;
        a = a.div_rem(g).1;
        a.resize(g.len());
        n >>= 1;
    }
    res.resize(g.len());
    res
}

pub fn fiduccia_with(md: Mod, c: &[u32], a: &[u32], k: usize, n: u64) -> u32 {
    let mut g = poly_from_c(md, c, k);
    g.set(k, 1);
    let t = calc(&g, n);
    let mut ans = 0u32;
    for i in 0..k {
        ans = md.add(ans, md.mul(t.coeff(i), a[i]));
    }
    ans
}

fn poly_from_c(md: Mod, c: &[u32], k: usize) -> Poly {
    let mut g = Poly::new(md, k + 1);
    for i in 1..=k {
        g.set(i - 1, md.neg(c[k - i + 1]));
    }
    g
}

pub fn bostan_mori(md: Mod, mut f: Poly, mut g: Poly, mut n: u64) -> u32 {
    while n != 0 {
        let mut t = g.clone();
        for i in (1..t.len()).step_by(2) {
            t.set(i, md.neg(t.coeff(i)));
        }
        f = &f * &t;
        g = &g * &t;
        let start = (n & 1) as usize;
        let mut nf = Poly::zero(md);
        let mut idx = 0;
        let mut i = start;
        while i < f.len() {
            nf.set(idx, f.coeff(i));
            idx += 1;
            i += 2;
        }
        f = nf;
        let mut ng = Poly::zero(md);
        idx = 0;
        i = 0;
        while i < g.len() {
            ng.set(idx, g.coeff(i));
            idx += 1;
            i += 2;
        }
        g = ng;
        n >>= 1;
    }
    if f.is_empty() {
        0
    } else {
        md.mul(
            f.coeff(0),
            md.inv(g.coeff(0)).expect("denominator not invertible"),
        )
    }
}

pub fn bostan_mori_cchlr(md: Mod, c: &[u32], a: &[u32], k: usize, n: u64) -> u32 {
    let mut g = Poly::new(md, k + 1);
    g.set(0, 1);
    for i in 1..=k {
        g.set(i, md.neg(c[i]));
    }
    let t = Poly::from_vec(md, a[..k].to_vec());
    let mut res = &t * &g;
    res.resize(k);
    bostan_mori(md, res, g, n)
}

pub fn bostam_mori_ccnhlr(md: Mod, c: &[u32], a: &[u32], k: usize, n: u64, p0: &Poly) -> u32 {
    let d = p0.len().saturating_sub(1);
    if k == 0 {
        let mut r = 0u32;
        for i in (0..=d).rev() {
            r = md.add(md.mul(r, n as u32), p0.coeff(i));
        }
        return r;
    }
    if (n as usize) < k {
        return a[n as usize];
    }
    let mut dd = Poly::new(md, k + 1);
    dd.set(0, 1);
    for i in 1..=k {
        dd.set(i, md.neg(c[i]));
    }
    let mut aa = Poly::new(md, k);
    for i in 0..k {
        aa.set(i, a[i]);
    }
    let mut cc = &aa * &dd;
    cc.resize(k);
    let fall = ffptp(p0);
    let mdmax = d + k + 2;
    let (fac, ifac) = factorials(md, mdmax);
    let mut u = Poly::new(md, d + 1);
    let mut v = Poly::new(md, d + 1);
    for j in 0..=d {
        u.set(j, md.mul(md.mul(fall.coeff(j), fac[j]), fac[d - j]));
        v.set(j, if j & 1 == 1 { md.neg(ifac[j]) } else { ifac[j] });
    }
    let mut r = &u * &v;
    for t in 0..=d {
        r.set(t, md.mul(r.coeff(t), ifac[d - t]));
    }
    r.resize(d + 1);
    let mut gg = Poly::new(md, k);
    for t in 0..k {
        gg.set(t, md.mul(md.mul(fac[d + t], ifac[d]), ifac[t]));
    }
    let mut kk = &r * &gg;
    kk.resize(k);
    let mut w = Poly::new(md, d + 2);
    for t in 0..=d + 1 {
        let b = md.mul(md.mul(fac[d + 1], ifac[t]), ifac[d + 1 - t]);
        w.set(t, if t & 1 == 1 { md.neg(b) } else { b });
    }
    let mut e = Poly::new(md, k);
    for i in 0..k {
        e.set(i, md.sub(cc.coeff(i), kk.coeff(i)));
    }
    let mut num = &e * &w;
    num.resize(k + d + 1);
    for i in 0..=d {
        num.set(i, md.add(num.coeff(i), r.coeff(i)));
    }
    let mut den = &w * &dd;
    den.resize(k + d + 2);
    bostan_mori(md, num, den, n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nonhomogeneous() {
        let md = Mod::new(998244353);
        let c = [0, 1, 1];
        let a = [0, 1];
        let p0 = Poly::from_vec(md, vec![0, 1]); // p(n) = n
        let mut w = vec![0u32; 40];
        w[0] = 0;
        w[1] = 1;
        for i in 2..40 {
            w[i] = md.add(md.add(w[i - 1], w[i - 2]), i as u32 % md.p());
        }
        for n in 2..40 {
            assert_eq!(
                bostam_mori_ccnhlr(md, &c, &a, 2, n as u64, &p0),
                w[n],
                "n={n}"
            );
        }
    }

    #[test]
    fn fib() {
        let md = Mod::new(998244353);
        // a_n = a_{n-1} + a_{n-2}, a0=0,a1=1
        let c = [0, 1, 1];
        let a = [0, 1];
        for n in 0..30 {
            let mut want = vec![0u32; 30];
            want[0] = 0;
            want[1] = 1;
            for i in 2..30 {
                want[i] = md.add(want[i - 1], want[i - 2]);
            }
            assert_eq!(fiduccia_with(md, &c, &a, 2, n as u64), want[n as usize]);
            assert_eq!(bostan_mori_cchlr(md, &c, &a, 2, n as u64), want[n as usize]);
        }
    }
}
