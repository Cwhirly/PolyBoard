use polyboard::modint::Mod;
use polyboard::ntt::convolution;
use std::io::{self, Read};

fn main() {
    let mut s = String::new();
    io::stdin().read_to_string(&mut s).unwrap();
    let mut it = s.split_whitespace();
    let n: usize = it.next().unwrap().parse().unwrap();
    let p: u32 = it.next().unwrap().parse().unwrap();
    let a: Vec<u32> = (0..n)
        .map(|_| it.next().unwrap().parse().unwrap())
        .collect();
    let b: Vec<u32> = (0..n)
        .map(|_| it.next().unwrap().parse().unwrap())
        .collect();
    let c = convolution(&a, &b, Mod::new(p));
    let mut h = 0u64;
    for &x in &c {
        h = h.wrapping_mul(131).wrapping_add(x as u64);
    }
    println!("{}", h);
}
