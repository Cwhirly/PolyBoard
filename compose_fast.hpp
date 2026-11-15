#include <algorithm>
#include <cassert>
#include <type_traits>
#include <vector>
namespace FastCompose {
template <unsigned Mod>
class ModInt {
    static_assert((Mod >> 31) == 0);
    static unsigned safe_mod(int v) {
        if ((v %= (int)Mod) < 0) v += (int)Mod;
        return v;
    }
    struct private_constructor_t {};
    static inline private_constructor_t private_constructor{};
    ModInt(private_constructor_t, unsigned v) : v_(v) {}
    unsigned v_;
public:
    static unsigned mod() { return Mod; }
    static ModInt from_raw(unsigned v) { return ModInt(private_constructor, v); }
    ModInt() : v_() {}
    ModInt(int v) : v_(safe_mod(v)) {}
    unsigned val() const { return v_; }
    ModInt operator-() const { return from_raw(v_ == 0 ? v_ : Mod - v_); }
    ModInt pow(int e) const {
        if (e < 0) return inv().pow(-e);
        for (ModInt x(*this), res(from_raw(1));; x *= x) {
            if (e & 1) res *= x;
            if ((e >>= 1) == 0) return res;
        }
    }
    ModInt inv() const {
        int x1 = 1, x3 = 0, a = v_, b = Mod;
        while (b) {
            int q = a / b, x1_old = x1, a_old = a;
            x1 = x3, x3 = x1_old - x3 * q, a = b, b = a_old - b * q;
        }
        return from_raw(x1 < 0 ? x1 + (int)Mod : x1);
    }
    ModInt &operator+=(const ModInt &a) {
        if ((v_ += a.v_) >= Mod) v_ -= Mod;
        return *this;
    }
    ModInt &operator-=(const ModInt &a) {
        if ((v_ += Mod - a.v_) >= Mod) v_ -= Mod;
        return *this;
    }
    ModInt &operator*=(const ModInt &a) {
        v_ = (unsigned long long)v_ * a.v_ % Mod;
        return *this;
    }
    friend ModInt operator+(const ModInt &a, const ModInt &b) { return ModInt(a) += b; }
    friend ModInt operator-(const ModInt &a, const ModInt &b) { return ModInt(a) -= b; }
    friend ModInt operator*(const ModInt &a, const ModInt &b) { return ModInt(a) *= b; }
    friend bool operator==(const ModInt &a, const ModInt &b) { return a.v_ == b.v_; }
    friend bool operator!=(const ModInt &a, const ModInt &b) { return a.v_ != b.v_; }
};
template <typename Tp>
class NttInfo {
    static Tp least_quadratic_nonresidue() {
        for (int i = 2;; ++i)
            if (Tp::from_raw(i).pow((Tp::mod() - 1) / 2) == -1) return Tp::from_raw(i);
    }
    const int ordlog2_;
    const Tp zeta_;
    const Tp invzeta_;
    mutable std::vector<Tp> root_;
    mutable std::vector<Tp> invroot_;
    NttInfo()
        : ordlog2_(__builtin_ctz(Tp::mod() - 1)),
          zeta_(least_quadratic_nonresidue().pow((Tp::mod() - 1) >> ordlog2_)),
          invzeta_(zeta_.inv()), root_{Tp::from_raw(1)}, invroot_{Tp::from_raw(1)} {}
public:
    static const NttInfo &get() { static NttInfo info; return info; }
    const std::vector<Tp> &root(int n) const {
        if (const int s = root_.size(); s < n) {
            root_.resize(n);
            for (int i = __builtin_ctz(s); (1 << i) < n; ++i) {
                const int j = 1 << i;
                root_[j] = zeta_.pow(1 << (ordlog2_ - i - 2));
                for (int k = j + 1; k < j * 2; ++k) root_[k] = root_[k - j] * root_[j];
            }
        }
        return root_;
    }
    const std::vector<Tp> &inv_root(int n) const {
        if (const int s = invroot_.size(); s < n) {
            invroot_.resize(n);
            for (int i = __builtin_ctz(s); (1 << i) < n; ++i) {
                const int j = 1 << i;
                invroot_[j] = invzeta_.pow(1 << (ordlog2_ - i - 2));
                for (int k = j + 1; k < j * 2; ++k) invroot_[k] = invroot_[k - j] * invroot_[j];
            }
        }
        return invroot_;
    }
};
template <typename Tp>
class Binomial {
    std::vector<Tp> factorial_, invfactorial_;
    Binomial() : factorial_{Tp::from_raw(1)}, invfactorial_{Tp::from_raw(1)} {}
    void preprocess(int n) {
        if (const int nn = factorial_.size(); nn < n) {
            int k = nn;
            while (k < n) k *= 2;
            factorial_.resize(k);
            invfactorial_.resize(k);
            for (int i = nn; i != k; ++i) factorial_[i] = factorial_[i - 1] * Tp::from_raw(i);
            invfactorial_.back() = factorial_.back().inv();
            for (int i = k - 2; i >= nn; --i)
                invfactorial_[i] = invfactorial_[i + 1] * Tp::from_raw(i + 1);
        }
    }
public:
    static const Binomial &get(int n) { static Binomial bin; bin.preprocess(n); return bin; }
    Tp binom(int n, int m) const { return n < m ? Tp() : factorial_[n] * invfactorial_[m] * invfactorial_[n - m]; }
};
int ntt_len(int n) {
    --n;
    n |= n >> 1, n |= n >> 2, n |= n >> 4, n |= n >> 8;
    return (n | n >> 16) + 1;
}
template <typename Tp>
void ntt(std::vector<Tp> &a) {
    const int n = a.size();
    for (int j = 0, l = n >> 1; j != l; ++j) {
        auto u = a[j], v = a[j + l];
        a[j] = u + v, a[j + l] = u - v;
    }
    auto &&root = NttInfo<Tp>::get().root(n / 2);
    for (int i = n >> 1; i >= 2; i >>= 1) {
        for (int j = 0, l = i >> 1; j != l; ++j) {
            auto u = a[j], v = a[j + l];
            a[j] = u + v, a[j + l] = u - v;
        }
        for (int j = i, l = i >> 1, m = 1; j != n; j += i, ++m)
            for (int k = j; k != j + l; ++k) {
                auto u = a[k], v = a[k + l] * root[m];
                a[k] = u + v, a[k + l] = u - v;
            }
    }
}
template <typename Tp>
void intt(std::vector<Tp> &a) {
    const int n = a.size();
    auto &&root = NttInfo<Tp>::get().inv_root(n / 2);
    for (int i = 2; i < n; i <<= 1) {
        for (int j = 0, l = i >> 1; j != l; ++j) {
            auto u = a[j], v = a[j + l];
            a[j] = u + v, a[j + l] = u - v;
        }
        for (int j = i, l = i >> 1, m = 1; j != n; j += i, ++m)
            for (int k = j; k != j + l; ++k) {
                auto u = a[k], v = a[k + l];
                a[k] = u + v, a[k + l] = (u - v) * root[m];
            }
    }
    const auto iv = Tp::from_raw(Tp::mod() - Tp::mod() / n);
    for (int j = 0, l = n >> 1; j != l; ++j) {
        auto u = a[j] * iv, v = a[j + l] * iv;
        a[j] = u + v, a[j + l] = u - v;
    }
}
template <typename Tp>
std::vector<Tp> convolution_ntt(std::vector<Tp> a, std::vector<Tp> b) {
    const int n = a.size(), m = b.size(), len = ntt_len(n + m - 1);
    a.resize(len), b.resize(len);
    ntt(a), ntt(b);
    for (int i = 0; i != len; ++i) a[i] *= b[i];
    intt(a);
    a.resize(n + m - 1);
    return a;
}
template <typename Tp>
std::vector<Tp> composition(const std::vector<Tp> &f, const std::vector<Tp> &g, int n) {
    if (n <= 0) return {};
    if (g.empty()) return std::vector<Tp>(n);
    struct rec {
        const std::vector<Tp> &f_;
        const Tp g0_;
        rec(const std::vector<Tp> &f, Tp g0) : f_(f), g0_(g0) {}
        std::vector<Tp> run(const std::vector<Tp> &Q, int d, int n) const {
            if (n == 0) {
                std::vector<Tp> P(d), invQ(d);
                for (int i = d - 1, j = 0; j < (int)f_.size() && i >= 0;) P[i--] = f_[j++];
                auto &&bin = Binomial<Tp>::get(d * 2);
                for (int i = 0; i < d; ++i) invQ[i] = bin.binom(d + i - 1, i) * g0_.pow(i);
                auto PinvQ = convolution_ntt(P, invQ);
                PinvQ.resize(d);
                return PinvQ;
            }
            const int len = ntt_len((d * 2 + 1) * (n * 2 + 2) - 1);
            std::vector<Tp> nttQ(len);
            for (int i = 0; i <= d; ++i)
                for (int j = 0; j <= n; ++j) nttQ[i * (n * 2 + 2) + j] = Q[i * (n + 1) + j];
            ntt(nttQ);
            std::vector<Tp> VV(len / 2);
            for (int i = 0; i != len; i += 2) VV[i / 2] = nttQ[i] * nttQ[i + 1];
            intt(VV);
            std::vector<Tp> V((d * 2 + 1) * (n / 2 + 1));
            for (int i = 0; i <= d * 2; ++i)
                for (int j = 0; j <= n / 2; ++j) V[i * (n / 2 + 1) + j] = VV[i * (n + 1) + j];
            const auto T = run(V, d * 2, n / 2);
            std::vector<Tp> nttT(len / 2);
            for (int i = 0; i < d * 2; ++i)
                for (int j = 0; j <= n / 2; ++j) nttT[i * (n + 1) + j] = T[i * (n / 2 + 1) + j];
            ntt(nttT);
            std::vector<Tp> UU(len);
            for (int i = 0; i != len; i += 2) {
                UU[i] = nttT[i / 2] * nttQ[i + 1];
                UU[i + 1] = nttT[i / 2] * nttQ[i];
            }
            intt(UU);
            std::vector<Tp> U(d * (n + 1));
            for (int i = 0; i < d; ++i)
                for (int j = 0; j <= n; ++j) U[i * (n + 1) + j] = UU[(i + d) * (n * 2 + 2) + j];
            return U;
        }
    } a(f, g[0]);
    std::vector<Tp> Q(n * 2);
    Q[0] = Tp::from_raw(1);
    for (int i = n, j = 0; j < (int)g.size() && i < n * 2;) Q[i++] = -g[j++];
    auto res = a.run(Q, 1, std::max(n - 1, (int)f.size() - 1));
    res.resize(n);
    return res;
}
}
