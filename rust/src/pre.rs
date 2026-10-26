use crate::modint::Mod;

pub fn factorials(m: Mod, n: usize) -> (Vec<u32>, Vec<u32>) {
    let mut fac = vec![1u32; n + 1];
    for i in 1..=n {
        fac[i] = m.mul(fac[i - 1], i as u32);
    }
    let mut ifac = vec![1u32; n + 1];
    ifac[n] = m.inv(fac[n]).expect("factorial not invertible");
    for i in (1..=n).rev() {
        ifac[i - 1] = m.mul(ifac[i], i as u32);
    }
    (fac, ifac)
}

pub fn inverse_table(m: Mod, n: usize) -> Vec<u32> {
    let mut inv = vec![0u32; n + 1];
    if n >= 1 {
        inv[1] = 1;
    }
    for i in 2..=n {
        let q = m.p() / i as u32;
        inv[i] = m.mul(m.sub(0, q), inv[(m.p() as usize) % i]);
    }
    inv
}
