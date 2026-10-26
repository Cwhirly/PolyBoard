use crate::modint::Mod;
use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};

const P1: u32 = 998244353;
const P2: u32 = 469762049;
const P3: u32 = 167772161;
const G: u32 = 3;
const NAIVE_LIMIT: usize = 4096;

fn is_prime(n: u32) -> bool {
    if n < 2 {
        return false;
    }
    for p in [2u32, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % p == 0 {
            return n == p;
        }
    }
    let mut d = n - 1;
    let mut s = 0;
    while d & 1 == 0 {
        d >>= 1;
        s += 1;
    }
    for a in [2u64, 3, 5, 7, 11, 13, 17, 19] {
        if a % n as u64 == 0 {
            continue;
        }
        let mut x = mod_pow(a as u32, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        let mut ok = false;
        for _ in 1..s {
            x = (x as u64 * x as u64 % n as u64) as u32;
            if x == n - 1 {
                ok = true;
                break;
            }
        }
        if !ok {
            return false;
        }
    }
    true
}

fn mod_pow(mut a: u32, mut e: u32, p: u32) -> u32 {
    let mut r = 1u32;
    while e != 0 {
        if e & 1 != 0 {
            r = (r as u64 * a as u64 % p as u64) as u32;
        }
        a = (a as u64 * a as u64 % p as u64) as u32;
        e >>= 1;
    }
    r
}

fn primitive_root(p: u32) -> u32 {
    static CACHE: OnceLock<Mutex<HashMap<u32, u32>>> = OnceLock::new();
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    if let Some(&g) = cache.lock().unwrap().get(&p) {
        return g;
    }
    let mut m = p - 1;
    let mut factors = Vec::new();
    let mut d = 2u32;
    while d as u64 * d as u64 <= m as u64 {
        if m % d == 0 {
            factors.push(d);
            while m % d == 0 {
                m /= d;
            }
        }
        d += 1;
    }
    if m > 1 {
        factors.push(m);
    }
    let mut g = 2u32;
    'outer: loop {
        for &q in &factors {
            if mod_pow(g, (p - 1) / q, p) == 1 {
                g += 1;
                continue 'outer;
            }
        }
        break;
    }
    cache.lock().unwrap().insert(p, g);
    g
}

fn ntt(a: &mut [u32], invert: bool, p: u32, g: u32) {
    let n = a.len();
    if n <= 1 {
        return;
    }
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            a.swap(i, j);
        }
    }
    let mut len = 2usize;
    while len <= n {
        let mut wlen = mod_pow(g, (p - 1) / len as u32, p);
        if invert {
            wlen = mod_pow(wlen, p - 2, p);
        }
        let half = len >> 1;
        for i in (0..n).step_by(len) {
            let mut w = 1u32;
            for j in 0..half {
                let u = a[i + j];
                let v = (a[i + j + half] as u64 * w as u64 % p as u64) as u32;
                let x = u as u64 + v as u64;
                a[i + j] = if x >= p as u64 {
                    (x - p as u64) as u32
                } else {
                    x as u32
                };
                let y = u as i64 - v as i64;
                a[i + j + half] = if y < 0 {
                    (y + p as i64) as u32
                } else {
                    y as u32
                };
                w = (w as u64 * wlen as u64 % p as u64) as u32;
            }
        }
        len <<= 1;
    }
    if invert {
        let inv_n = mod_pow(n as u32, p - 2, p);
        for x in a.iter_mut() {
            *x = (*x as u64 * inv_n as u64 % p as u64) as u32;
        }
    }
}

fn convolution_naive(a: &[u32], b: &[u32], m: Mod) -> Vec<u32> {
    let mut c = vec![0u32; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        if x == 0 {
            continue;
        }
        for (j, &y) in b.iter().enumerate() {
            c[i + j] = m.add(c[i + j], m.mul(x, y));
        }
    }
    c
}

fn convolution_prime(a: &[u32], b: &[u32], p: u32, g: u32) -> Vec<u32> {
    let need = a.len() + b.len() - 1;
    let n = need.next_power_of_two();
    let mut fa = vec![0u32; n];
    let mut fb = vec![0u32; n];
    for i in 0..a.len() {
        fa[i] = a[i] % p;
    }
    for i in 0..b.len() {
        fb[i] = b[i] % p;
    }
    ntt(&mut fa, false, p, g);
    ntt(&mut fb, false, p, g);
    for i in 0..n {
        fa[i] = (fa[i] as u64 * fb[i] as u64 % p as u64) as u32;
    }
    ntt(&mut fa, true, p, g);
    fa.truncate(need);
    fa
}

fn inv_mod(a: u32, p: u32) -> u32 {
    mod_pow(a, p - 2, p)
}

fn convolution_3(a: &[u32], b: &[u32], m: Mod) -> Vec<u32> {
    let c1 = convolution_prime(a, b, P1, G);
    let c2 = convolution_prime(a, b, P2, G);
    let c3 = convolution_prime(a, b, P3, G);
    let m1 = P1 as u128;
    let m2 = P2 as u128;
    let m3 = P3 as u128;
    let m12 = m1 * m2;
    let i12 = inv_mod((m1 % m2) as u32, P2) as u128;
    let i123 = inv_mod((m12 % m3) as u32, P3) as u128;
    let mut out = Vec::with_capacity(c1.len());
    for i in 0..c1.len() {
        let x1 = c1[i] as u128;
        let t2 = ((c2[i] as u128 + m2 - x1 % m2) % m2 * i12) % m2;
        let x12 = x1 + m1 * t2;
        let t3 = ((c3[i] as u128 + m3 - x12 % m3) % m3 * i123) % m3;
        let x = x12 + m12 * t3;
        out.push((x % m.p() as u128) as u32);
    }
    out
}

pub fn convolution(a: &[u32], b: &[u32], m: Mod) -> Vec<u32> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    if a.len() * b.len() <= NAIVE_LIMIT {
        return convolution_naive(a, b, m);
    }
    let need = a.len() + b.len() - 1;
    let n = need.next_power_of_two();
    let p = m.p();
    if is_prime(p) && (p - 1) as usize % n == 0 {
        let g = primitive_root(p);
        return convolution_prime(a, b, p, g);
    }
    convolution_3(a, b, m)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn naive_and_ntt_agree() {
        for &p in &[998244353u32, 1000000007, 12] {
            let m = Mod::new(p);
            let a: Vec<u32> = (0..50).map(|i| (i * 1234567 % p as usize) as u32).collect();
            let b: Vec<u32> = (0..60).map(|i| (i * 7654321 % p as usize) as u32).collect();
            let c = convolution(&a, &b, m);
            let d = convolution_naive(&a, &b, m);
            assert_eq!(c, d, "p={p}");
        }
    }
}
