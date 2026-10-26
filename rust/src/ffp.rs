use crate::modint::Mod;
use crate::multipoint::pmpe;
use crate::poly::Poly;
use crate::pre::factorials;

fn build(md: Mod, l: usize, r: usize, f: &Poly) -> (Poly, Poly) {
    if l == r {
        return (
            Poly::from_vec(md, vec![f.coeff(l)]),
            Poly::from_vec(md, vec![md.neg(l as u32), 1]),
        );
    }
    let mid = (l + r) >> 1;
    let (la, lt) = build(md, l, mid, f);
    let (ra, rt) = build(md, mid + 1, r, f);
    (&la + &(&lt * &ra), &lt * &rt)
}

pub fn ptffp(f: &Poly) -> Poly {
    let md = f.modulus();
    let n = f.len();
    if n == 0 {
        return Poly::zero(md);
    }
    let (ans, _) = build(md, 0, n - 1, f);
    ans
}

pub fn ffptp(f: &Poly) -> Poly {
    let md = f.modulus();
    let n = f.len();
    if n == 0 {
        return Poly::zero(md);
    }
    let xs: Vec<u32> = (0..=n).map(|i| i as u32).collect();
    let vals = pmpe(f, &xs);
    let (_, ifac) = factorials(md, n);
    let mut g = Poly::new(md, n + 1);
    let mut h = Poly::new(md, n + 1);
    for i in 0..=n {
        g.set(i, md.mul(vals[i], ifac[i]));
        h.set(i, if i & 1 == 1 { md.neg(ifac[i]) } else { ifac[i] });
    }
    let mut r = &g * &h;
    r.resize(n + 1);
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip() {
        let md = Mod::new(998244353);
        let f = Poly::from_vec(md, (0..8).map(|i| (i * 17 + 3) as u32).collect());
        let c = ffptp(&f);
        let back = ptffp(&c);
        for i in 0..f.len() {
            assert_eq!(back.coeff(i), f.coeff(i), "i={i}");
        }
    }
}
