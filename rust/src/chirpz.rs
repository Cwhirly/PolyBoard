use crate::modint::Mod;
use crate::poly::Poly;

pub fn ci2(m: Mod, i: u32) -> u32 {
    let i = i as u64;
    (i * (i.wrapping_sub(1)) / 2 % (m.p() as u64 - 1)) as u32
}

pub fn chirp_z(f: &Poly, c: u32, m: usize) -> Poly {
    let md = f.modulus();
    let n = f.len();
    if n == 0 || m == 0 {
        return Poly::zero(md);
    }
    if c == 0 {
        let mut r = Poly::new(md, m);
        let mut s = 0u32;
        for i in 0..n {
            s = md.add(s, f.coeff(i));
        }
        r.set(0, s);
        for i in 1..m {
            r.set(i, f.coeff(0));
        }
        return r;
    }
    let l = n + m - 1;
    let mut f = f.clone();
    f.resize(l);
    let mut r = Poly::new(md, l);
    let mut g = Poly::new(md, l);
    let ic = md.inv(c).expect("c not invertible");
    let (mut w, mut iw, mut p, mut ip) = (1u32, 1u32, 1u32, 1u32);
    for i in 0..l {
        r.set(l - 1 - i, w);
        g.set(i, md.mul(iw, f.coeff(i)));
        w = md.mul(w, p);
        iw = md.mul(iw, ip);
        p = md.mul(p, c);
        ip = md.mul(ip, ic);
    }
    let t = &r * &g;
    let mut res = Poly::new(md, m);
    let (mut v, mut vp) = (1u32, 1u32);
    for i in 0..m {
        res.set(i, md.mul(t.coeff(l - 1 - i), v));
        v = md.mul(v, vp);
        vp = md.mul(vp, ic);
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modint::Mod;

    #[test]
    fn chirp_small() {
        for &pp in &[998244353u32, 1000000007, 12] {
            let md = Mod::new(pp);
            for c in [0u32, 1, 2, 3, 12345] {
                if c != 0 && md.inv(c).is_none() {
                    continue;
                }
                let n = 12;
                let m = 9;
                let f = Poly::from_vec(md, (0..n).map(|i| (i * 37 + 11) as u32 % pp).collect());
                let got = chirp_z(&f, c, m);
                let mut cp = 1u32;
                for i in 0..m {
                    let mut want = 0u32;
                    let mut pw = 1u32;
                    for j in 0..n {
                        want = md.add(want, md.mul(f.coeff(j), pw));
                        pw = md.mul(pw, cp);
                    }
                    assert_eq!(got.coeff(i), want, "p={pp} c={c} i={i}");
                    cp = md.mul(cp, c);
                }
            }
        }
    }
}
