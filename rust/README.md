# PolyBoard Rust 1.0.0

Rust port of the C++ PolyBoard template. The Rust version starts at 1.0.0 and mirrors C++ 1.19.0; later C++ updates are synced here.

## Status

| Area | C++ | Rust |
| --- | --- | --- |
| Modular arithmetic / Barrett | yes | `modint::Mod` |
| NTT / arbitrary-modulus convolution | yes | `ntt::convolution`（单 NTT 或三素数 CRT） |
| `poly` basics, inverse, div/rem, pow, high part | yes | `poly::Poly` |
| FPS `ln` / `exp` / `sqrt` / `pow` | yes | `fps` |
| Matrix add/sub/mul/pow/gauss/inv/det | yes | `matrix::Matrix` |
| FWT / subset convolution | yes | `fwt::Sps` |
| Stirling row/column | yes | `stirling` |
| Chirp-Z | yes | `chirpz::chirp_z` |
| FPS trig / `Trans` | yes | `fps` |
| Multipoint evaluation / interpolation | yes | `multipoint` |
| FFP / Staling / CCHLR / CCNHLR | yes | `ffp` / `staling` / `cchlr` |
| High part (`operator&`) | yes | `Poly::high_mul` |

## Build

```bash
cargo test
```

## Notes

- The modulus is runtime (`Mod::new(p)`), so one build supports all moduli.
- Convolution uses one NTT when `p` is NTT-friendly, otherwise three NTT primes plus CRT.
- C++ `P_FIXED` / `ATN` switches are not needed: the Rust path selects NTT or CRT at runtime.
- The Rust API is idiomatic rather than a line-by-line transliteration; algorithms and limits follow the C++ version.
- Convolution was cross-checked against C++ at `n = 10^5` for `p = 998244353` (NTT) and `p = 1000000007` (three-prime CRT); hashes match.
- `Chirp-Z` additionally handles `c = 0` directly, which the short C++ implementation does not.
