use crate::poly::Poly;

pub fn ln(f: &Poly, n: usize) -> Poly {
    if n == 0 {
        return Poly::zero(f.modulus());
    }
    let inv = f.inv(n);
    let mut r = f
        .derivative()
        .mul_trunc(&inv, n.saturating_sub(1))
        .integral();
    r.resize(n);
    r
}

pub fn exp(f: &Poly, n: usize) -> Poly {
    if n == 0 {
        return Poly::zero(f.modulus());
    }
    let md = f.modulus();
    let mut g = Poly::one(md);
    let mut m = 1usize;
    while m < n {
        let m2 = (m << 1).min(n);
        let lg = ln(&g, m2);
        let mut t = f.trunc(m2) - lg;
        t.set(0, md.add(t.coeff(0), 1));
        g = g.mul_trunc(&t, m2);
        m = m2;
    }
    g.resize(n);
    g
}

pub fn sqrt(f: &Poly, n: usize) -> Option<Poly> {
    if n == 0 {
        return Some(Poly::zero(f.modulus()));
    }
    let md = f.modulus();
    let s = md.sqrt(f.coeff(0))?;
    let mut g = Poly::from_vec(md, vec![s]);
    let inv2 = md.inv(2)?;
    let mut m = 1usize;
    while m < n {
        let m2 = (m << 1).min(n);
        let ig = g.inv(m2);
        let t = f.trunc(m2).mul_trunc(&ig, m2);
        g = ((&g + &t) * inv2).trunc(m2);
        g.resize(m2);
        m = m2;
    }
    g.resize(n);
    Some(g)
}

pub fn pow(f: &Poly, k: u64, n: usize) -> Poly {
    let md = f.modulus();
    if n == 0 {
        return Poly::zero(md);
    }
    if k == 0 {
        let mut r = Poly::one(md);
        r.resize(n);
        return r;
    }
    let mut f = f.clone();
    if f.coeff(0) == 0 {
        let mut z = 0usize;
        while z < f.len() && f.coeff(z) == 0 {
            z += 1;
        }
        if z == f.len() {
            return Poly::zero(md);
        }
        let shift = z as u64 * k;
        if shift >= n as u64 {
            return Poly::zero(md);
        }
        f = f >> z;
        let r = pow(&f, k, n - shift as usize);
        return r << shift as usize;
    }
    let c = f.coeff(0);
    let invc = md.inv(c).expect("constant term not invertible");
    let h = &f * invc;
    let l = ln(&h, n);
    let mut e = &l * (k % md.p() as u64) as u32;
    e = exp(&e, n);
    e = &e * md.pow(c, k);
    e.resize(n);
    e
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modint::Mod;

    #[test]
    fn ln_exp_roundtrip() {
        let md = Mod::new(998244353);
        let f = Poly::from_vec(md, vec![1, 3, 5, 7, 9, 11]);
        let n = 6;
        let l = ln(&f, n);
        let e = exp(&l, n);
        assert_eq!(e.as_slice(), f.as_slice());
    }

    #[test]
    fn sqrt_square() {
        let md = Mod::new(998244353);
        let f = Poly::from_vec(md, vec![1, 2, 3, 4, 5]);
        let s = sqrt(&f, 5).unwrap();
        let sq = s.mul_trunc(&s, 5);
        assert_eq!(sq.as_slice(), f.as_slice());
    }

    #[test]
    fn pow_small() {
        let md = Mod::new(998244353);
        let f = Poly::from_vec(md, vec![1, 1]);
        let p = pow(&f, 3, 4);
        assert_eq!(p.as_slice(), &[1, 3, 3, 1]);
    }
}
