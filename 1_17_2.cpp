#include <bits/stdc++.h>

using namespace std;
namespace IO {
    const int __SIZE = (1 << 21) + 1;
    char ibuf[__SIZE], *iS, *iT, obuf[__SIZE],
        *oS = obuf, *oT = oS + __SIZE - 1, _c, qu[55];
    int __f, qr, _eof;
#define Gc()                                                                   \
    (iS == iT ? (iT = (iS = ibuf) + fread(ibuf, 1, __SIZE, stdin),             \
                 (iS == iT ? EOF : *iS++))                                     \
              : *iS++)

    void flush() {
        fwrite(obuf, 1, oS - obuf, stdout);
        oS = obuf;
    }

    void gc(char &x) { x = Gc(); }

    void pc(char x) {
        *oS++ = x;
        if (oS == oT)
            flush();
    }

    void pstr(const char *s) {
        int __len = strlen(s);
        for (__f = 0; __f < __len; ++__f)
            pc(s[__f]);
    }

    void gstr(char *s) {
        for (_c = Gc(); _c < 32 || _c > 126 || _c == ' ';)
            _c = Gc();
        for (; _c > 31 && _c < 127 && _c != ' ' && _c != '\n' && _c != '\r';
             ++s, _c = Gc())
            *s = _c;
        *s = 0;
    }

    template <class I> bool read(I &x) {
        _eof = 0;
        for (__f = 1, _c = Gc(); (_c < '0' || _c > '9') && !_eof; _c = Gc()) {
            if (_c == '-')
                __f = -1;
            _eof |= _c == EOF;
        }
        for (x = 0; _c <= '9' && _c >= '0' && !_eof; _c = Gc()) {
            x = x * 10 + (_c & 15), _eof |= _c == EOF;
        }
        x *= __f;
        return !_eof;
    }

    template <class I> void print(I x) {
        if (!x)
            pc('0');
        if (x < 0) {
            pc('-');
            x = -x;
        }
        while (x) {
            qu[++qr] = x % 10 + '0', x /= 10;
        }
        while (qr)
            pc(qu[qr--]);
    }

    struct Flusher_ {
        ~Flusher_() { flush(); }
    } io_flusher_;
}
using IO::gc;
using IO::gstr;
using IO::pc;
using IO::print;
using IO::pstr;
using IO::read;

namespace CONSTS {
    constexpr int __qpow(int __a, int __b, int __p) {
        long long __r = 1, __x = __a;
        while (__b) {
            if (__b & 1)
                __r = __r * __x % __p;
            __x = __x * __x % __p, __b >>= 1;
        }
        return __r;
    }

    constexpr bool __isroot(int __g, int __p) {
        int __m = __p - 1;
        for (int __i = 2; 1ll * __i * __i <= __m; ++__i)
            if (__m % __i == 0) {
                if (__qpow(__g, (__p - 1) / __i, __p) == 1)
                    return false;
                while (__m % __i == 0)
                    __m /= __i;
            }
        return __m == 1 || __qpow(__g, (__p - 1) / __m, __p) != 1;
    }

    constexpr int __root(int __p) {
        for (int __g = 2;; ++__g)
            if (__isroot(__g, __p))
                return __g;
    }

    constexpr int __v2(int __x) {
        int __r = 0;
        while (!(__x & 1))
            __x >>= 1, ++__r;
        return __r;
    }

    constexpr int __lgup(int __x) {
        int __r = 0;
        while ((1 << __r) < __x)
            ++__r;
        return __r;
    }

    const constexpr int P = 998244353, N = 600005, Y = __root(P),
                        I = __qpow(Y, P - 2, P), B = (P + 1) >> 1,
                        _I_ = __qpow(Y, (P - 1) / 4, P), LG = __lgup(2 * N - 1),
                        NL = 1 << LG;
    static_assert(P % 4 == 1, "P must be 1 mod 4");
    static_assert(LG <= __v2(P - 1), "N is too large for P");
}
using namespace CONSTS;

namespace Pre {
    int gp[NL], igp[NL], ny[N], fac[N], ifac[N];

    int Q(int __a, int __b) {
        int __res = 1;
        while (__b) {
            if (__b & 1)
                __res = 1ll * __res * __a % P;
            __a = 1ll * __a * __a % P, __b >>= 1;
        }
        return __res % P;
    }

    int Inv(int __x) { return Q(__x, P - 2); }

    void initYG() {
        int __off = 0;
        for (int __len = NL; __len >= 2; __len >>= 1) {
            long long __g = Q(Y, (P - 1) / __len), __ig = Q(I, (P - 1) / __len);
            int __w = 1, __iw = 1;
            for (int j = 0; j < (__len >> 1); ++j) {
                gp[__off + j] = __w, igp[__off + j] = __iw;
                __w = __w * __g % P, __iw = __iw * __ig % P;
            }
            __off += __len >> 1;
        }
        ny[1] = 1;
        for (int i = 2; i < N; ++i)
            ny[i] = 1ll * (P - P / i) * ny[P % i] % P;
        fac[0] = ifac[0] = 1;
        for (int i = 1; i < N; ++i)
            fac[i] = 1ll * fac[i - 1] * i % P;
        ifac[N - 1] = Inv(fac[N - 1]);
        for (int i = N - 1; i >= 1; --i)
            ifac[i - 1] = 1ll * ifac[i] * i % P;
    }
}

namespace Pint {
    template <class T> T addt(T &__a, T __b) {
        if ((__a += __b) >= P)
            __a -= P;
        return __a;
    }

    template <class T> T delt(T &__a, T __b) {
        if ((__a -= __b) < 0)
            __a += P;
        return __a;
    }

    template <class T> T add(T __a, T __b) { return addt(__a, __b); }

    template <class T> T del(T __a, T __b) { return delt(__a, __b); }

    template <class T> T tp(T x) {
        while (x < 0)
            x += P;
        return x;
    }

    template <class T> T Vmax(T a, T b) { return ((a > b) ? a : b); }
}

using namespace Pint;
namespace Quad {
    mt19937 rnd(time(0));

    struct NTC {
        int Re, Im;
    };

    class Solver {
        int ii;

        static int Q(int __a, int __b) {
            int __res = 1;
            while (__b) {
                if (__b & 1)
                    __res = 1ll * __res * __a % P;
                __a = 1ll * __a * __a % P, __b >>= 1;
            }
            return __res % P;
        }

        NTC mul(const NTC &__a, const NTC &__b) const {
            NTC __res;
            __res.Re =
                (1ll * __a.Re * __b.Re + 1ll * ii * __a.Im % P * __b.Im) % P;
            __res.Im = (1ll * __a.Re * __b.Im + 1ll * __a.Im * __b.Re) % P;
            return __res;
        }

        NTC q(NTC __a, int __b) const {
            NTC __res = {1, 0};
            while (__b) {
                if (__b & 1)
                    __res = mul(__res, __a);
                __a = mul(__a, __a), __b >>= 1;
            }
            return __res;
        }

        int Cipolla(int __n) {
            if (Q(__n, (P - 1) >> 1) == P - 1)
                return -1;
            int __a;
            do {
                __a = rnd() % P;
                ii = (1ll * __a * __a % P - __n + P) % P;
            } while (Q(ii, (P - 1) >> 1) != P - 1);
            return q({__a, 1}, B).Re;
        }

    public:
        int work(int __n) {
            if (!__n)
                return 0;
            int __u = Cipolla(__n % P);
            if (__u == -1)
                return P + 1;
            int __v = P - __u;
            if (__u > __v)
                swap(__u, __v);
            return __u;
        }
    };

    int work(int __n) { return Solver().work(__n); }
}

namespace POLY {
    class poly;

    namespace PIO {
        poly pin(int __n);
        void ppri(const poly &__f, int __n);
    }

    namespace UCPF {
        poly Ln(poly __f);
        poly Exp(poly __f);
        poly Sqrt(poly __f);
        poly Sin(poly __f);
        poly Cos(poly __f);
        poly Tan(poly __f);
        poly ArcSin(poly __f);
        poly ArcCos(poly __f);
        poly ArcTan(poly __f);
        poly AND(poly __f, poly __g);
        poly OR(poly __f, poly __g);
        poly XOR(poly __f, poly __g);
        poly Subset(poly __f, poly __g);
        poly Pow_For_Luogu(poly __f, int __k1, int __k2, int __kk);
        poly Trans(poly __f, int __c);
    }

    using namespace UCPF;
    using namespace PIO;

    const constexpr int MUL_LIM = 4096;

    class poly {
        vector<int> a;

    public:
        poly() {}
        poly(int __lim_siz) { resize(__lim_siz); }

        int size() const { return a.size(); }

        bool empty() const { return a.empty(); }

        int &operator[](int __id) { return a[__id]; }

        const int &operator[](int __id) const { return a[__id]; }

        void resize(int __lim_siz) { a.resize(__lim_siz); }

        void reserve(int __lim_siz) { a.reserve(__lim_siz); }

        void pb(int __x) { a.push_back(__x); }

        void rev() { reverse(a.begin(), a.end()); }

        void rev(int __l, int __r) {
            reverse(a.begin() + __l, a.begin() + __r);
        }

        void tp() {
            for (int i = 0; i < size(); ++i)
                a[i] = Pint::tp(a[i]);
        }

        poly operator+=(const poly &__t) {
            int __n = Vmax(size(), __t.size());
            resize(__n);
            for (int i = 0; i < __t.size(); ++i)
                a[i] = add(a[i], __t[i]);
            return *this;
        }

        poly operator-=(const poly &__t) {
            int __n = Vmax(size(), __t.size());
            resize(__n);
            for (int i = 0; i < __t.size(); ++i)
                a[i] = del(a[i], __t[i]);
            return *this;
        }

        poly operator+(poly __t) const {
            poly __res = *this;
            return __t += __res;
        }

        poly operator-(poly __t) const {
            poly __res = *this;
            return __res -= __t;
        }

        void NTT(bool __) {
            int __n = size();
            if (__) {
                for (int __len = __n; __len >= 2; __len >>= 1) {
                    int __half = __len >> 1;
                    const int *__w = Pre::gp + (NL - __len);
                    for (int i = 0; i < __n; i += __len) {
                        for (int j = 0; j < __half; ++j) {
                            int __u = a[i + j], __v = a[i + j + __half];
                            a[i + j] = add(__u, __v);
                            a[i + j + __half] =
                                1ll * __w[j] * del(__u, __v) % P;
                        }
                    }
                }
            } else {
                for (int __len = 2; __len <= __n; __len <<= 1) {
                    int __half = __len >> 1;
                    const int *__w = Pre::igp + (NL - __len);
                    for (int i = 0; i < __n; i += __len) {
                        for (int j = 0; j < __half; ++j) {
                            int __u = a[i + j],
                                __v = 1ll * __w[j] * a[i + j + __half] % P;
                            a[i + j] = add(__u, __v);
                            a[i + j + __half] = del(__u, __v);
                        }
                    }
                }
                long long __inv = Pre::Inv(__n);
                for (int i = 0; i < __n; ++i)
                    a[i] = a[i] * __inv % P;
            }
        }

        poly operator*=(const poly &__t) {
            int __n = size(), __m = __t.size();
            if (!__n || !__m) {
                a.clear();
                return *this;
            }
            if (1ll * __n * __m <= MUL_LIM) {
                poly __res(__n + __m - 1);
                for (int i = 0; i < __n; ++i) {
                    if (!a[i])
                        continue;
                    long long __ai = a[i];
                    for (int j = 0; j < __m; ++j)
                        __res[i + j] = (__res[i + j] + __ai * __t[j]) % P;
                }
                return *this = __res;
            }
            int __len = 1;
            while (__len < __n + __m - 1)
                __len <<= 1;
            poly __u(__len);
            for (int i = 0; i < __m; ++i)
                __u[i] = __t[i];
            resize(__len);
            NTT(1);
            __u.NTT(1);
            for (int i = 0; i < __len; ++i)
                a[i] = 1ll * a[i] * __u[i] % P;
            NTT(0);
            resize(__n + __m - 1);
            return *this;
        }

        poly operator*=(int __t) {
            if (__t < 0)
                __t += P;
            for (int i = 0; i < size(); ++i)
                a[i] = 1ll * a[i] * __t % P;
            return *this;
        }

        poly operator*(const poly &__t) const {
            poly __res = *this;
            return __res *= __t;
        }

        poly operator*(int __t) const {
            poly __res = *this;
            return __res *= __t;
        }

        poly operator&(const poly &__t) const {
            int __n = size(), __m = __t.size(), __len = 1;
            while (__len < __n)
                __len <<= 1;
            poly __u = *this, __v = __t;
            __u.resize(__len);
            __v.resize(__len);
            __u.NTT(1);
            __v.NTT(1);
            for (int i = 0; i < __len; ++i)
                __u[i] = 1ll * __u[i] * __v[i] % P;
            __u.NTT(0);
            poly __res(__n - __m + 1);
            for (int i = __m - 1; i < __n; ++i)
                __res[i - __m + 1] = __u[i];
            return __res;
        }

        poly operator~() const {
            int __n = size(), __dep = 1, __lim = 1;
            poly __res, __f;
            if (!__n)
                return __res;
            while (__lim < (__n << 1))
                __lim <<= 1;
            __res.reserve(__lim);
            __res.pb(Pre::Inv(a[0]));
            while (__dep < (__n << 1)) {
                __f = *this;
                __f.resize(__dep);
                __f.resize(__dep << 1);
                if (__dep != 1)
                    __res.resize(__dep >> 1);
                __res.resize(__dep << 1);
                __f.NTT(1);
                __res.NTT(1);
                for (int i = 0; i < (__dep << 1); ++i)
                    __res[i] = 1ll * __res[i] *
                               (2 + P - 1ll * __f[i] * __res[i] % P) % P;
                __res.NTT(0);
                __dep <<= 1;
            }
            __res.resize(__n);
            return __res;
        }

        poly operator/(const poly &__t) const {
            int __n = size(), __m = __t.size();
            poly __f = *this, __g = __t;
            __f.rev();
            __g.rev();
            __g.resize(__n - __m + 1);
            poly __res = ~__g * __f;
            __res.resize(__n - __m + 1);
            __res.rev();
            return __res;
        }

        poly operator/=(const poly &__t) { return *this = *this / __t; }

        poly operator/=(int __t) {
            if (__t < 0)
                __t += P;
            __t = Pre::Inv(__t);
            for (int i = 0; i < size(); ++i)
                a[i] = 1ll * a[i] * __t % P;
            return *this;
        }

        poly operator/(int __t) const {
            poly __res = *this;
            return __res /= __t;
        }

        poly operator%=(const poly &__t) {
            if (__t.size() > size())
                return *this;
            poly __q = *this / __t;
            return *this = *this - __q * __t;
        }

        poly operator%(const poly &__t) const {
            poly __res = *this;
            return __res %= __t;
        }

        poly operator<<=(int __t) {
            a.insert(a.begin(), __t, 0);
            return *this;
        }

        poly operator<<(int __t) const {
            poly __res = *this;
            return __res <<= __t;
        }

        poly operator>>=(int __t) {
            if (__t >= size()) {
                a.clear();
                return *this;
            }
            a.erase(a.begin(), a.begin() + __t);
            return *this;
        }

        poly operator>>(int __t) const {
            poly __res = *this;
            return __res >>= __t;
        }

        poly operator^=(int __t) {
            int __n = size(), __fir = 0;
            while (__fir < __n && !a[__fir])
                ++__fir;
            if (__fir == __n) {
                *this = poly(__n);
                if (!__t && __n)
                    a[0] = 1;
                return *this;
            }
            int __con = a[__fir];
            *this >>= __fir;
            *this /= __con;
            *this = Ln(*this);
            *this *= __t;
            *this = Exp(*this);
            *this *= Pre::Q(__con, __t);
            *this <<= (__fir * __t);
            return *this;
        }

        poly operator^(int __t) const {
            poly __res = *this;
            return __res ^= __t;
        }
    };

    using ffp = poly;

    class sps {
        vector<int> a;

    public:
        sps() {}
        sps(int __lim_siz) { resize(__lim_siz); }
        explicit sps(const poly &__f) {
            resize(__f.size());
            for (int i = 0; i < __f.size(); ++i)
                a[i] = __f[i];
        }

        int size() const { return a.size(); }

        bool empty() const { return a.empty(); }

        int &operator[](int __id) { return a[__id]; }

        const int &operator[](int __id) const { return a[__id]; }

        void resize(int __lim_siz) { a.resize(__lim_siz); }

        void reserve(int __lim_siz) { a.reserve(__lim_siz); }

        poly to_poly() const {
            poly __f(size());
            for (int i = 0; i < size(); ++i)
                __f[i] = a[i];
            return __f;
        }

        void Fwtand(int __) {
            int __n = size();
            for (int x = 2; x <= __n; x <<= 1) {
                int k = x >> 1;
                for (int i = 0; i < __n; i += x)
                    for (int j = 0; j < k; ++j)
                        a[i + j] = __ ? add(a[i + j], a[i + j + k])
                                      : del(a[i + j], a[i + j + k]);
            }
        }

        void Fwtor(int __) {
            int __n = size();
            for (int x = 2; x <= __n; x <<= 1) {
                int k = x >> 1;
                for (int i = 0; i < __n; i += x)
                    for (int j = 0; j < k; ++j)
                        a[i + j + k] = __ ? add(a[i + j + k], a[i + j])
                                          : del(a[i + j + k], a[i + j]);
            }
        }

        void Fmtor(int __) { Fwtor(__); }

        void Fwtxor(int __) {
            int __n = size();
            for (int x = 2; x <= __n; x <<= 1) {
                int k = x >> 1;
                for (int i = 0; i < __n; i += x)
                    for (int j = 0; j < k; ++j) {
                        int __u = a[i + j], __v = a[i + j + k];
                        a[i + j] = add(__u, __v);
                        a[i + j + k] = del(__u, __v);
                    }
            }
            if (!__) {
                long long __inv = Pre::Inv(__n);
                for (int i = 0; i < __n; ++i)
                    a[i] = a[i] * __inv % P;
            }
        }

        sps operator&(const sps &__g) const {
            sps __f = *this, __h = __g;
            __f.Fwtand(1);
            __h.Fwtand(1);
            for (int i = 0; i < __f.size(); ++i)
                __f[i] = 1ll * __f[i] * __h[i] % P;
            __f.Fwtand(0);
            return __f;
        }

        sps operator|(const sps &__g) const {
            sps __f = *this, __h = __g;
            __f.Fwtor(1);
            __h.Fwtor(1);
            for (int i = 0; i < __f.size(); ++i)
                __f[i] = 1ll * __f[i] * __h[i] % P;
            __f.Fwtor(0);
            return __f;
        }

        sps operator^(const sps &__g) const {
            sps __f = *this, __h = __g;
            __f.Fwtxor(1);
            __h.Fwtxor(1);
            for (int i = 0; i < __f.size(); ++i)
                __f[i] = 1ll * __f[i] * __h[i] % P;
            __f.Fwtxor(0);
            return __f;
        }

        sps operator*(const sps &__g) const {
            int __n = size(), __lgn = __lg(__n);
            vector<sps> __F(__lgn + 1, sps(__n)), __G(__lgn + 1, sps(__n)),
                __H(__lgn + 1, sps(__n));
            for (int i = 0; i < __n; ++i) {
                __F[__builtin_popcount(i)][i] = a[i];
                __G[__builtin_popcount(i)][i] = __g[i];
            }
            for (int i = 0; i <= __lgn; ++i) {
                __F[i].Fwtor(1);
                __G[i].Fwtor(1);
            }
            for (int i = 0; i <= __lgn; ++i)
                for (int k = 0; k <= i; ++k) {
                    int *__h = &__H[i][0], *__x = &__F[k][0],
                        *__y = &__G[i - k][0];
                    for (int S = 0; S < __n; ++S)
                        __h[S] = (__h[S] + 1ll * __x[S] * __y[S]) % P;
                }
            for (int i = 0; i <= __lgn; ++i)
                __H[i].Fwtor(0);
            sps __ans(__n);
            for (int i = 0; i < __n; ++i)
                __ans[i] = __H[__builtin_popcount(i)][i];
            return __ans;
        }

        sps operator*=(int __t) {
            if (__t < 0)
                __t += P;
            for (int i = 0; i < size(); ++i)
                a[i] = 1ll * a[i] * __t % P;
            return *this;
        }

        sps operator*(int __t) const {
            sps __res = *this;
            return __res *= __t;
        }

        sps And(const sps &__g) const { return *this & __g; }

        sps Or(const sps &__g) const { return *this | __g; }

        sps Xor(const sps &__g) const { return *this ^ __g; }

        sps Subset(const sps &__g) const { return *this * __g; }
    };
}

using namespace POLY;
poly Dx(const poly &__f) {
    poly __ans;
    __ans.reserve(__f.empty() ? 0 : __f.size() - 1);
    for (int i = 1; i < __f.size(); ++i)
        __ans.pb(1ll * i * __f[i] % P);
    return __ans;
}

poly Integ(const poly &__f) {
    poly __ans;
    __ans.reserve(__f.size() + 1);
    __ans.pb(0);
    for (int i = 0; i < __f.size(); ++i)
        __ans.pb(1ll * Pre::ny[i + 1] * __f[i] % P);
    return __ans;
}

poly UCPF::Ln(poly __f) {
    int __n = __f.size();
    poly __d = Dx(__f);
    __d *= ~__f;
    __d.resize(__n);
    __d = Integ(__d);
    __d.resize(__n);
    return __d;
}

poly UCPF::Exp(poly __f) {
    int __n = __f.size();
    poly __ans;
    int __dep = 1;
    __ans.pb(1);
    while (__dep < (__n << 1)) {
        poly __lnf = Ln(__ans);
        __lnf = __f - __lnf;
        __lnf.resize(__dep << 1);
        ++__lnf[0];
        __ans *= __lnf;
        __ans.resize(__dep << 1);
        __dep <<= 1;
    }
    __ans.resize(__n);
    return __ans;
}

poly UCPF::Sqrt(poly __f) {
    int __con = __f[0];
    __f *= Pre::Inv(__con);
    __f = Exp(Ln(__f) * B);
    __f *= (__con != 1) ? Quad::work(__con) : 1;
    return __f;
}

poly UCPF::Sin(poly __f) {
    return (Exp(__f * _I_) - Exp(__f * (P - _I_))) *
           (1ll * B * Pre::Inv(_I_) % P);
}

poly UCPF::Cos(poly __f) { return (Exp(__f * _I_) + Exp(__f * (P - _I_))) * B; }

poly UCPF::Tan(poly __f) { return Sin(__f) * ~Cos(__f); }

poly UCPF::ArcSin(poly __f) {
    poly __t = __f;
    int __n = __f.size();
    __f = __f * __f;
    __f.resize(__n);
    __f *= (P - 1);
    __f[0]++;
    __f = Sqrt(__f);
    __f = ~__f;
    __f *= Dx(__t);
    __f.resize(__n);
    __f = Integ(__f);
    __f.resize(__n);
    return __f;
}

poly UCPF::ArcCos(poly __f) { return ArcSin(__f) * (P - 1); }

poly UCPF::ArcTan(poly __f) {
    poly __t = __f;
    int __n = __f.size();
    __f = __f * __f;
    __f.resize(__n);
    __f[0]++;
    __f = ~__f;
    __f *= Dx(__t);
    __f.resize(__n);
    __f = Integ(__f);
    __f.resize(__n);
    return __f;
}

poly UCPF::AND(poly __f, poly __g) { return (sps(__f) & sps(__g)).to_poly(); }

poly UCPF::OR(poly __f, poly __g) { return (sps(__f) | sps(__g)).to_poly(); }

poly UCPF::XOR(poly __f, poly __g) { return (sps(__f) ^ sps(__g)).to_poly(); }

poly UCPF::Subset(poly __f, poly __g) {
    return (sps(__f) * sps(__g)).to_poly();
}

poly UCPF::Pow_For_Luogu(poly __f, int __k1, int __k2, int __kk) {
    int __n = __f.size(), __fir = 0;
    while (__fir < __n && !__f[__fir])
        ++__fir;
    if (__fir == __n) {
        __f = poly(__n);
        if (!__kk && __n)
            __f[0] = 1;
        return __f;
    }
    if (1ll * __kk * __fir >= __n)
        return __f = poly(__n);
    int __con = __f[__fir];
    __f >>= __fir;
    __f /= __con;
    __f = Exp(Ln(__f) * __k2);
    __f *= Pre::Q(__con, __k1);
    __f <<= (__fir * __k2);
    return __f;
}

poly UCPF::Trans(poly __f, int __c) {
    int __n = __f.size();
    poly __g(__n), __h(__n);
    for (int i = 0; i < __n; ++i)
        __g[i] = 1ll * __f[i] * Pre::fac[i] % P;
    long long __pw = 1;
    for (int i = 0; i < __n; ++i) {
        __h[i] = __pw * Pre::ifac[i] % P;
        __pw = __pw * __c % P;
    }
    __g.rev();
    __h *= __g;
    __h.resize(__n);
    __h.rev();
    for (int i = 0; i < __n; ++i)
        __h[i] = 1ll * __h[i] * Pre::ifac[i] % P;
    return __h;
}

poly PIO::pin(int __n) {
    poly __f;
    __f.reserve(__n);
    for (int i = 0; i < __n; ++i) {
        int x;
        read(x);
        __f.pb(x);
    }
    return __f;
}

void PIO::ppri(const poly &__f, int __n) {
    int __s = __f.size();
    for (int i = 0; i < __n; ++i) {
        print(i < __s ? Pint::tp(__f[i]) : 0);
        pc(' ');
    }
    pc('\n');
}

namespace Multipoint_Evel {
    struct Multipoint {
        vector<poly> t;
        vector<int> an;

        void build(int l, int r, int p, const vector<int> &__x) {
            if (p == 1) {
                t.assign((r - l + 1) << 2, poly());
                an.assign(r - l + 1, 0);
            }
            if (l == r) {
                t[p].reserve(2);
                t[p].pb(1);
                t[p].pb(P - __x[l]);
                return;
            }
            int mid = (l + r) >> 1;
            build(l, mid, p << 1, __x);
            build(mid + 1, r, p << 1 | 1, __x);
            t[p] = t[p << 1] * t[p << 1 | 1];
        }

        void query(int l, int r, int p, poly __f) {
            if (l == r) {
                an[l] = __f[0];
                return;
            }
            poly __u = __f & t[p << 1 | 1], __v = __f & t[p << 1];
            int mid = (l + r) >> 1;
            query(l, mid, p << 1, __u);
            query(mid + 1, r, p << 1 | 1, __v);
        }

        poly eval(const poly &__f, const vector<int> &__x) {
            int __n = __f.size(), __m = __x.size();
            if (!__m)
                return poly();
            if (!__n)
                return poly(__m);
            poly __g = __f;
            int __con = __g[0];
            __g.resize(__m + 1);
            build(0, __m - 1, 1, __x);
            __g.rev();
            __g.resize((__m + 1) << 1);
            poly __tmp = __g * ~t[1];
            __tmp.resize(__m);
            query(0, __m - 1, 1, __tmp);
            poly __ans(__m);
            for (int i = 0; i < __m; ++i)
                __ans[i] = (1ll * __x[i] * an[i] + __con) % P;
            return __ans;
        }
    };

    poly PMPE(const poly &__f, const vector<int> &__x) {
        return Multipoint().eval(__f, __x);
    }
}

using Multipoint_Evel::PMPE;

namespace Fast_Interpolation {
    struct Interpolation {
        vector<poly> q, ret;

        void build(int l, int r, int p, const vector<int> &__x) {
            if (p == 1) {
                q.assign((r - l + 1) << 2, poly());
                ret.assign((r - l + 1) << 2, poly());
            }
            if (l == r) {
                q[p].reserve(2);
                q[p].pb(P - __x[l]);
                q[p].pb(1);
                return;
            }
            int mid = (l + r) >> 1;
            build(l, mid, p << 1, __x);
            build(mid + 1, r, p << 1 | 1, __x);
            q[p] = q[p << 1] * q[p << 1 | 1];
        }

        void solve(int l, int r, int p, const vector<int> &__v) {
            if (l == r) {
                ret[p].pb(__v[l]);
                return;
            }
            int mid = (l + r) >> 1;
            solve(l, mid, p << 1, __v);
            solve(mid + 1, r, p << 1 | 1, __v);
            ret[p] = ret[p << 1] * q[p << 1 | 1] + ret[p << 1 | 1] * q[p << 1];
        }

        poly interpolate(const vector<int> &__x, const vector<int> &__y) {
            int __n = __x.size();
            if (!__n)
                return poly();
            build(0, __n - 1, 1, __x);
            poly __f = Dx(q[1]);
            int __con = __f[0];
            __f.resize(__n + 1);
            Multipoint_Evel::Multipoint __mp;
            __mp.build(0, __n - 1, 1, __x);
            __f.rev();
            __f.resize(__n << 1);
            poly __tmp = __f * ~__mp.t[1];
            __tmp.resize(__n);
            __mp.query(0, __n - 1, 1, __tmp);
            vector<int> __v(__n);
            for (int i = 0; i < __n; ++i)
                __v[i] = 1ll * __y[i] *
                         Pre::Inv((1ll * __x[i] * __mp.an[i] + __con) % P) % P;
            solve(0, __n - 1, 1, __v);
            return std::move(ret[1]);
        }
    };

    poly PFI(const vector<int> &__x, const vector<int> &__y) {
        return Interpolation().interpolate(__x, __y);
    }
}

using Fast_Interpolation::PFI;

namespace ChirpZ {
    int Ci2(int __i) { return 1ll * __i * (__i - 1) / 2 % (P - 1); }

    poly Chirp_Z(poly __f, int __c, int __m) {
        int __n = __f.size(), __l = __n + __m - 1;
        __f.resize(__l);
        poly __r(__l), __g(__l);
        long long __w = 1, __iw = 1, __p = 1, __ip = 1, __ic = Pre::Inv(__c);
        for (int i = 0; i < __l; ++i) {
            __r[__l - 1 - i] = __w;
            __g[i] = __iw * __f[i] % P;
            __w = __w * __p % P, __iw = __iw * __ip % P;
            __p = __p * __c % P, __ip = __ip * __ic % P;
        }
        poly __t = __r * __g;
        poly __res(__m);
        long long __v = 1, __vp = 1;
        for (int i = 0; i < __m; ++i) {
            __res[i] = __t[__l - 1 - i] * __v % P;
            __v = __v * __vp % P, __vp = __vp * __ic % P;
        }
        return __res;
    }
}

using ChirpZ::Chirp_Z;

namespace Stirling {
    poly row(int __n) {
        poly __f(__n + 1), __g(__n + 1);
        __f[0] = __g[0] = 1;
        for (int i = 1; i <= __n; ++i) {
            __f[i] = Pre::ifac[i];
            __g[i] = (i & 1) ? P - Pre::ifac[i] : Pre::ifac[i];
        }
        for (int i = 0; i <= __n; ++i)
            __f[i] = 1ll * __f[i] * Pre::Q(i, __n) % P;
        __f *= __g;
        __f.resize(__n + 1);
        return __f;
    }

    void Dev_Mul(int l, int r, int p, vector<poly> &__t) {
        if (l == r) {
            __t[p].reserve(2);
            __t[p].pb(1);
            __t[p].pb(P - l);
            return;
        }
        int mid = (l + r) >> 1;
        Dev_Mul(l, mid, p << 1, __t);
        Dev_Mul(mid + 1, r, p << 1 | 1, __t);
        __t[p] = __t[p << 1] * __t[p << 1 | 1];
    }

    poly column(int __k, int __n) {
        poly __ans;
        if (__k > __n) {
            __ans.resize(__n + 1);
            return __ans;
        }
        vector<poly> __t((__k << 2) + 5);
        Dev_Mul(1, __k, 1, __t);
        __t[1].resize(__n + 1);
        __ans = (~__t[1]) << __k;
        __ans.resize(__n + 1);
        return __ans;
    }
}

namespace Staling {
    poly row(int __n) {
        poly __ans;
        if (__n == 1) {
            __ans.pb(0);
            __ans.pb(1);
            return __ans;
        }
        poly __l = row(__n >> 1);
        __l *= Trans(__l, __n >> 1);
        if (((__n >> 1) << 1) == __n)
            return __l;
        return __l = __l * (__n - 1) + (__l << 1);
    }

    poly column(int __k, int __n) {
        poly __f;
        __f.pb(0);
        for (int i = 1; i <= __n; ++i)
            __f.pb(Pre::ny[i]);
        __f = Exp(Ln(__f >> 1) * __k) << __k;
        __f.resize(__n + 1);
        __f /= Pre::fac[__k];
        for (int i = 0; i <= __n; ++i)
            __f[i] = 1ll * __f[i] * Pre::fac[i] % P;
        return __f;
    }
}

namespace FFP {
    void Dev_Mul(int l, int r, int p, const poly &__a, vector<poly> &__t,
                 vector<poly> &__ans) {
        if (l == r) {
            __t[p].reserve(2);
            __t[p].pb(P - l);
            __t[p].pb(1);
            __ans[p].pb(__a[l]);
            return;
        }
        int mid = (l + r) >> 1;
        Dev_Mul(l, mid, p << 1, __a, __t, __ans);
        Dev_Mul(mid + 1, r, p << 1 | 1, __a, __t, __ans);
        __t[p] = __t[p << 1] * __t[p << 1 | 1];
        __ans[p] = __ans[p << 1] + __t[p << 1] * __ans[p << 1 | 1];
    }

    ffp PTFFP(const poly &__f) {
        int __n = __f.size();
        vector<poly> __t((__n << 2) + 5), __ans((__n << 2) + 5);
        Dev_Mul(0, __n - 1, 1, __f, __t, __ans);
        return std::move(__ans[1]);
    }

    poly FFPTP(ffp __f) {
        int __n = __f.size();
        vector<int> __x(__n + 1);
        for (int i = 0; i <= __n; ++i)
            __x[i] = i;
        poly __r = Multipoint_Evel::Multipoint().eval(__f, __x);
        poly __g(__n + 1), __h(__n + 1);
        for (int i = 0; i <= __n; ++i) {
            __g[i] = 1ll * __r[i] * Pre::ifac[i] % P;
            __h[i] = (i & 1) ? P - Pre::ifac[i] : Pre::ifac[i];
        }
        __g *= __h;
        __g.resize(__n + 1);
        return __g;
    }
}

namespace CCHLR {
    poly __calc(const poly &__g, int __b) {
        poly __res, __a;
        __res.pb(1);
        __a.pb(0), __a.pb(1);
        while (__b) {
            if (__b & 1) {
                __res = __res * __a;
                __res %= __g;
                __res.resize(__g.size());
            }
            __a = __a * __a;
            __a %= __g;
            __a.resize(__g.size());
            __b >>= 1;
        }
        __res.resize(__g.size());
        return __res;
    }

    int Fiduccia(int *c, int *a, int k, int n) {
        poly __g;
        for (int i = 1; i <= k; ++i)
            __g.pb((P - c[k - i + 1]) % P);
        __g.pb(1);
        poly __t = __calc(__g, n);
        long long __ans = 0;
        for (int i = 0; i < k; ++i)
            __ans = (__ans + 1ll * __t[i] * a[i]) % P;
        return __ans;
    }

    int Bostan_Mori(poly __f, poly __g, int __n) {
        while (__n) {
            poly __t = __g;
            for (int i = 1; i < __g.size(); i += 2)
                __t[i] = (P - __g[i]) % P;
            __f *= __t;
            __g *= __t;
            int i = (__n & 1);
            for (; i < __f.size(); i += 2)
                __f[i >> 1] = __f[i];
            __f.resize(i >> 1);
            for (i = 0; i < __g.size(); i += 2)
                __g[i >> 1] = __g[i];
            __g.resize(i >> 1);
            __n >>= 1;
        }
        if (!__f.size())
            return 0;
        return 1ll * __f[0] * Pre::Inv(__g[0]) % P;
    }

    int Bostan_Mori_CCHLR(int *c, int *a, int k, int n) {
        poly __g;
        __g.pb(1);
        for (int i = 1; i <= k; ++i)
            __g.pb((P - c[i]) % P);
        poly __t;
        for (int i = 0; i < k; ++i)
            __t.pb((a[i] % P + P) % P);
        poly __res = __t * __g;
        __res.resize(k);
        return Bostan_Mori(__res, __g, n);
    }
}

using CCHLR::Bostan_Mori;
using CCHLR::Bostan_Mori_CCHLR;
using CCHLR::Fiduccia;

namespace CCNHLR {
    poly __calc(const poly &__g, int __b) {
        poly __res, __a;
        __res.pb(1);
        __a.pb(0), __a.pb(1);
        while (__b) {
            if (__b & 1) {
                __res = __res * __a;
                __res %= __g;
                __res.resize(__g.size());
            }
            __a = __a * __a;
            __a %= __g;
            __a.resize(__g.size());
            __b >>= 1;
        }
        __res.resize(__g.size());
        return __res;
    }

    poly __ksm(int __b) {
        poly __res, __a;
        __res.pb(1);
        __a.pb(P - 1), __a.pb(1);
        while (__b) {
            if (__b & 1)
                __res *= __a;
            __a *= __a;
            __b >>= 1;
        }
        return __res;
    }

    int Bostam_Mori_CCNHLR(int *c, int *a, int k, int n, const poly &__p0) {
        int __d = __p0.size() - 1;
        if (!k) {
            long long __r = 0;
            for (int i = __d; i >= 0; --i)
                __r = (__r * n + __p0[i]) % P;
            return __r;
        }
        if (n < k)
            return a[n];
        poly __D(k + 1);
        __D[0] = 1;
        for (int i = 1; i <= k; ++i)
            __D[i] = (P - c[i]) % P;
        poly __A(k);
        for (int i = 0; i < k; ++i)
            __A[i] = a[i];
        poly __C = __A * __D;
        __C.resize(k);
        poly __fall = FFP::FFPTP(__p0);
        poly __u(__d + 1), __v(__d + 1);
        for (int j = 0; j <= __d; ++j) {
            __u[j] = 1ll * __fall[j] * Pre::fac[j] % P * Pre::fac[__d - j] % P;
            __v[j] = (j & 1) ? P - Pre::ifac[j] : Pre::ifac[j];
        }
        poly __R = __u * __v;
        for (int t = 0; t <= __d; ++t)
            __R[t] = 1ll * __R[t] * Pre::ifac[__d - t] % P;
        __R.resize(__d + 1);
        poly __G(k);
        for (int t = 0; t < k; ++t)
            __G[t] =
                1ll * Pre::fac[__d + t] * Pre::ifac[__d] % P * Pre::ifac[t] % P;
        poly __K = __R * __G;
        __K.resize(k);
        poly __W(__d + 2);
        for (int t = 0; t <= __d + 1; ++t) {
            int __b = 1ll * Pre::fac[__d + 1] * Pre::ifac[t] % P *
                      Pre::ifac[__d + 1 - t] % P;
            __W[t] = (t & 1) ? P - __b : __b;
        }
        poly __E(k);
        for (int i = 0; i < k; ++i)
            __E[i] = (__C[i] + P - __K[i]) % P;
        poly __Num = __E * __W;
        __Num.resize(k + __d + 1);
        for (int i = 0; i <= __d; ++i)
            __Num[i] = (__Num[i] + __R[i]) % P;
        poly __Den = __W * __D;
        __Den.resize(k + __d + 2);
        return Bostan_Mori(__Num, __Den, n);
    }
}

using CCNHLR::Bostam_Mori_CCNHLR;

int n, m, c;

inline void work()
{
    poly f;
    read(n);
    read(c);
    read(m);
    
    f=pin(n);
    ppri(Chirp_Z(f, c, m), m);
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(0), cout.tie(0);
    Pre::initYG();
    work();
    return 0;
}
