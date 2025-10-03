/*

POLYBOARD (1.15.2)

Menu:
IO ------------------------ Fastio (From:luuia)
Pre ----------------------- Preworks
    | Q ----------------------- Fastpow
    | Inv --------------------- Get number theory inv
    | initYG ------------------ Serve for NTT
Pint ---------------------- Integers Calculate in Modulo
    | addt -------------------- Plus to
    | delt -------------------- Minus to
    | add --------------------- Plus
    | del --------------------- Minus
    | tp ---------------------- To Positive
    | Vmax -------------------- Manumax
Quad ----------------------- Quadratic Residue
    | NTC --------------------- Manucomplex
    | Q ----------------------- Fastpow with Any Modulo
    | q ----------------------- Conplex Fastpow
    | Cipolla ----------------- Cipolla (
    | work -------------------- Get Quadratic Residue
POLY ----------------------- POLY
    | PIO --------------------- Poly Input/Output
        | pin --------------------- Statement of Input
        | ppri -------------------- Statement of Output
    | Poly -------------------- Class of Poly
        | resize ------------------ Resize ans Shrink
        | shrink ------------------ Update n
        | rev --------------------- Reverse
        | pb ---------------------- Push Back
        | tp ---------------------- To Positive
        | NTT --------------------- NTT (
        | Fwtand ------------------ Convolute And
        | Fwtor ------------------- Convolute Or
        | Fwtxor ------------------ Convolute Xor
UCPF ----------------------- Unclassable Polyfunctions
    | Dx ----------------------- Get Derivative
    | Integ -------------------- Get Integral
    | Ln ----------------------- Get Logarithm of e
    | Exp ---------------------- Get Power of e
    | Sqrt --------------------- Get Sqrt
    | Sin ---------------------- Get Sin
    | Cos ---------------------- Get Cos
    | Tan ---------------------- Get Tan
    | ArcSin ------------------- Get ArcSin
    | ArcCos ------------------- Get ArcCos
    | ArcTan ------------------- Get ArcTan
    | And ---------------------- And Together
    | Or ----------------------- Or Together
    | Xor ---------------------- Xor Together
    | Pow_For_Luogu ------------ For Extremely Large Pow Number
    | Trans -------------------- Translation
Multipoint_Evel ------------ Multipoint_Evel
    | build -------------------- Conquer Multiplication
    | multipoint_evel ---------- Calculate
    | PMPE --------------------- For Direct Calling
Fast_Interpolation --------- Fast_Interpolation
    | build1 ------------------- Conquer Multiplication
    | Inv_evel ----------------- Calculate
    | PFI ---------------------- For Direct Calling
Chirp_Z -------------------- Chirp-Z Transform
    | Ci2 ---------------------- Calculate \binom{i}{2}
    | Chirp_Z ------------------ Calculate
Stirling ------------------- Second Stirling Numbers
    | row ---------------------- For the Same Row
    | Dev_Mul ------------------ Conquer Multiplication
    | Column ------------------- For the Same Colume
Staling -------------------- First Stirling Numbers
    | row ---------------------- For the Same Row
    | Column ------------------- For the Same Column
FFP ------------------------ Descending Power Polynomial iff Normal Polynomial
    | Dev_Mul ------------------ Conquer Multiplication
    | PTFFP -------------------- Normal Polynomial to Descending Power Polynomial
    | FFPTP -------------------- Descending Power Polynomial to Normal Polynomial
CCHLR ---------------------- Constant Coefficient Homogeneous Linear Recursion
    | __calc ------------------- Polynomial Fast Pow
    | Fiduccia ----------------- Fiduccia Algorithm
CCNLR --------------------- Constant Coefficient Nonhomogeneous Linear Recursion
Update:





*/

#include <bits/stdc++.h>

using namespace std;
namespace IO
{
    const int __SIZE = (1 << 21) + 1;
    char ibuf[__SIZE], *iS, *iT, obuf[__SIZE], *oS = obuf, *oT = oS + __SIZE - 1, _c, qu[55];
    int __f, qr, _eof;
#define Gc() (iS == iT ? (iT = (iS = ibuf) + fread(ibuf, 1, __SIZE, stdin), (iS == iT ? EOF : *iS++)) : *iS++)

    inline void flush()
    {
        fwrite(obuf, 1, oS - obuf, stdout);
        oS = obuf;
    }

    inline void gc(char &x)
    {
        x = Gc();
    }

    inline void pc(char x)
    {
        *oS++ = x;
        if (oS == oT)
        {
            flush();
        }
    }

    inline void pstr(const char *s)
    {
        int __len = strlen(s);
        for (__f = 0; __f < __len; ++__f)
        {
            pc(s[__f]);
        }
    }

    inline void gstr(char *s)
    {
        for (_c = Gc(); _c < 32 || _c > 126 || _c == ' ';)
        {
            _c = Gc();
        }
        for (; _c > 31 && _c < 127 && _c != ' ' && _c != '\n' && _c != '\r'; ++s, _c = Gc())
        {
            *s = _c;
        }
        *s = 0;
    }

    template <class I>
    inline bool read(I &x)
    {
        _eof = 0;
        for (__f = 1, _c = Gc(); (_c < '0' || _c > '9') && !_eof; _c = Gc())
        {
            if (_c == '-')
            {
                __f = -1;
            }
            _eof |= _c == EOF;
        }
        for (x = 0; _c <= '9' && _c >= '0' && !_eof; _c = Gc())
        {
            x = x * 10 + (_c & 15), _eof |= _c == EOF;
        }
        x *= __f;
        return !_eof;
    }

    template <class I>
    inline void print(I x)
    {
        if (!x)
        {
            pc('0');
        }
        if (x < 0)
        {
            pc('-');
            x = -x;
        }
        while (x)
        {
            qu[++qr] = x % 10 + '0', x /= 10;
        }
        while (qr)
        {
            pc(qu[qr--]);
        }
    }

    struct Flusher_
    {
        ~Flusher_() { flush(); }
    } io_flusher_;
}
using IO::gc;
using IO::gstr;
using IO::pc;
using IO::print;
using IO::pstr;
using IO::read;

namespace CONSTS
{
    const constexpr int P = 998244353, Y = 3, I = 332748118, B = (P + 1) >> 1, N = 600005, _I_ = 86583718;
}
using namespace CONSTS;

namespace Pre
{
    int gp[N], igp[N], ny[N], inv[N];

    inline int Q(int a, int b)
    {
        int res = 1;
        while (b)
        {
            if (b & 1)
            {
                res = 1ll * res * a % P;
            }
            a = 1ll * a * a % P, b >>= 1;
        }
        return res % P;
    }

    inline int Inv(int __x)
    {
        return Q(__x, P - 2);
    }

    inline void initYG()
    {
        int tmpg = Q(Y, (P - 1) / (1 << 19)), tmpig = Q(I, (P - 1) / (1 << 19));
        gp[0] = igp[0] = 1;
        for (register int i = 1; i <= (1 << 19); ++i)
        {
            gp[i] = 1ll * tmpg * gp[i - 1] % P;
            igp[i] = 1ll * tmpig * igp[i - 1] % P;
            
        }
        ny[1] = 1;
        for (int i = 2; i <= 500010; ++i)
        {
            ny[i] = 1ll * (P - P / i) * ny[P % i] % P;
        }
    }
}

namespace Pint
{
    template <class T>
    inline T addt(T &__a, T __b)
    {
        if ((__a += __b) >= P)
        {
            __a -= P;
        }
        return __a;
    }

    template <class T>
    inline T delt(T &__a, T __b)
    {
        if ((__a -= __b) < 0)
        {
            __a += P;
        }
        return __a;
    }

    template <class T>
    inline T add(T __a, T __b)
    {
        return addt(__a, __b);
    }

    template <class T>
    inline T del(T __a, T __b)
    {
        return delt(__a, __b);
    }

    template <class T>
    inline T tp(T x)
    {
        while (x < 0)
        {
            x += P;
        }
        return x;
    }

    template <class T>
    inline T Vmax(T a, T b)
    {
        return ((a > b) ? a : b);
    }
}

using namespace Pint;
mt19937 rnd(time(0));
namespace Quad
{
    int t, n, p = 998244353, ii;

    struct NTC
    {
        int Re, Im;
        NTC operator*(NTC __A) const
        {
            NTC __res;
            __res.Re = (1ll * Re * __A.Re % p + 1ll * ii * Im % p * __A.Im % p + p) % p;
            __res.Im = (1ll * Re * __A.Im % p + 1ll * Im * __A.Re % p + p) % p;
            return __res;
        }
    };
    inline int Q(int __a, int __b, int __p)
    {
        int __res = 1;
        while (__b)
        {
            if (__b & 1)
            {
                __res = 1ll * __res * __a % __p;
            }
            __a = 1ll * __a * __a % __p;
            __b >>= 1;
        }
        return __res % __p;
    }

    inline NTC q(NTC __a, int __b, int &__p)
    {
        NTC __res = {1, 0};
        while (__b)
        {
            if (__b & 1)
            {
                __res = __res * __a;
            }
            __a = __a * __a;
            __b >>= 1;
        }
        return __res;
    }
    inline int Cipolla(int n, int p)
    {
        p = 998244353;
        n %= p;
        if (Q(n, (p - 1) >> 1, p) == p - 1)
            return -1;
        int a;
        while (true)
        {
            a = rnd() % p;
            ii = ((1ll * a * a % p - n + p) % p + p) % p;
            if (Q(ii, (p - 1) >> 1, 998244353) == p - 1)
                break;
        }
        NTC x = {a, 1};
        return (q(x, B, p).Re % p + p) % p;
    }
    inline int work(int n, int p)
    {
        if (!n)
        {
            return 0;
        }
        int u = Cipolla(n, p), v = p - u;
        if (u == -1)
        {
            return 998244354;
        }
        u = (u + p) % p;
        v = (v + p) % p;
        if (u > v)
            swap(u, v);
        return u;
    }
}

namespace POLY
{
    int __Binary_reverse[2000005];
    class poly;
    class ffp;
    namespace PIO
    {
        void pin(poly &f, int __n);
        void ppri(poly __f, int __n);
        void pin(ffp &f, int __n);
        void ppri(ffp __f, int __n);
    }

    namespace UCPF
    {
        inline poly Ln(poly __f);
        inline poly Exp(poly __f);
        inline poly Sqrt(poly __f);
        inline poly Sin(poly __f);
        inline poly Cos(poly __f);
        inline poly Tan(poly __f);
        inline poly ArcSin(poly __f);
        inline poly ArcCos(poly __f);
        inline poly ArcTan(poly __f);
        inline poly AND(poly __f, poly __g);
        inline poly OR(poly __f, poly __g);
        inline poly XOR(poly __f, poly __g);
        inline poly Pow_For_Luogu(poly __f, int __k1, int __k2, int kk);
        inline poly Trans(poly __f, int c);
    }
    using namespace UCPF;
    using namespace PIO;

    class poly
    {
    public:
        int n;
        vector<int> a;
        poly()
        {
            n = 0;
        }
        poly(int __n)
        {
            resize(__n);
        }

        int &operator[](int id)
        {
            return a[id];
        }

        inline void resize(int __lim_siz)
        {
            a.resize(__lim_siz);
            a.shrink_to_fit();
            n = __lim_siz;
        }

        inline void shrink()
        {
            resize(n);
        }

        vector<int>::iterator begin()
        {
            return a.begin();
        }

        vector<int>::iterator end()
        {
            return a.end();
        }

        inline void rev(vector<int>::iterator __st, vector<int>::iterator __en)
        {
            reverse(__st, __en);
        }

        inline void pb(int x)
        {
            a.push_back(x);
            ++n;
        }

        inline void tp(int __n)
        {
            for (register int i = 0; i < __n; ++i)
            {
                a[i] = Pint::tp(a[i]);
            }
        }

        poly operator+=(const poly &__tmpa)
        {
            int nm = Vmax(n, __tmpa.n);
            resize(nm);
            for (register int i = 0; i < nm; ++i)
            {
                a[i] = add((i < __tmpa.n) ? __tmpa.a[i] : 0, (i < n) ? a[i] : 0);
            }
            return *this;
        }

        poly operator-=(const poly &__tmpa)
        {
            int nm = Vmax(n, __tmpa.n);
            resize(nm);
            for (register int i = 0; i < nm; ++i)
            {
                a[i] = del((i < n) ? a[i] : 0, (i < __tmpa.n) ? __tmpa.a[i] : 0);
            }
            return *this;
        }

        poly operator+(poly __tmpa) const
        {
            poly __tmp = *this;
            return (__tmpa += __tmp);
        }

        poly operator-(poly __tmpa) const
        {
            poly __tmp = *this;
            return (__tmp -= __tmpa);
        }

        inline void NTT(bool __)
        {
            for (register int i = 0; i < n; ++i)
            {
                if (i < __Binary_reverse[i])
                {
                    swap(a[i], a[__Binary_reverse[i]]);
                }
            }
            for (register int r = 1, base = (1 << 18); r < n; r <<= 1, base >>= 1)
            {
                for (register int i = 0; i < n; i += (r << 1))
                {
                    for (register int j = 0, bj = 0; j < r; ++j, bj += base)
                    {
                        int __I = a[j | i], __II = 1ll * (__ ? Pre::gp[bj] : Pre::igp[bj]) * a[j | r | i] % P;
                        a[j | i] = add(__I, __II);
                        a[j | r | i] = del(__I, __II);
                    }
                }
            }
            if (!__)
            {
                long long invtmp = Pre::Inv(n);
                for (register int i = 0; i < n; ++i)
                {
                    a[i] = a[i] * invtmp % P;
                }
            }
            tp(n);
        }
        inline void Fwtand(int __)
        {
            if (!__)
                __ = -1;
            for (register int x = 2; x <= n; x <<= 1)
            {
                int k = x >> 1;
                for (register int i = 0; i < n; i += x)
                {
                    for (register int j = 0; j < k; ++j)
                    {
                        a[i + j] = 1ll * add(1ll * a[i + j], 1ll * a[i + j + k] * __ % P) % P;
                    }
                }
            }
        }
        inline void Fwtor(int __)
        {
            if (!__)
                __ = -1;
            for (register int x = 2; x <= n; x <<= 1)
            {
                int k = x >> 1;
                for (register int i = 0; i < n; i += x)
                {
                    for (register int j = 0; j < k; ++j)
                    {
                        a[i + j + k] = 1ll * add(1ll * a[i + j + k], 1ll * a[i + j] * __ % P) % P;
                    }
                }
            }
        }
        inline void Fwtxor(int __)
        {
            if (!__)
                __ = B;
            for (register int x = 2; x <= n; x <<= 1)
            {
                int k = x >> 1;
                for (register int i = 0; i < n; i += x)
                {
                    for (register int j = 0; j < k; ++j)
                    {
                        a[i + j] = add(a[i + j], a[i + j + k]);
                        a[i + j + k] = (a[i + j] - (a[i + j + k] << 1) % P + P) % P;
                        a[i + j] = 1ll * __ * a[i + j] % P;
                        a[i + j + k] = 1ll * __ * a[i + j + k] % P;
                    }
                }
            }
        }
        poly operator*=(const poly &__tmpa)
        {
            poly __tmpmuly = __tmpa;
            int __siz = n + __tmpa.n - 1;
            int __len = 1, __L = 0;
            while (__len < n + __tmpa.n)
            {
                __len <<= 1;
                __L++;
            }
            for (register int i = 1; i < __len; ++i)
            {
                __Binary_reverse[i] = (__Binary_reverse[i >> 1] >> 1) | ((i & 1) << (__L - 1));
            }
            resize(__len);
            __tmpmuly.resize(__len);
            NTT(1);
            __tmpmuly.NTT(1);
            for (register int i = 0; i < __len; ++i)
            {
                a[i] = 1ll * a[i] * __tmpmuly[i] % P;
            }
            NTT(0);
            resize(__siz);
            return *this;
        }

        poly operator*=(const int __tmpa)
        {
            for (register int i = 0; i < n; ++i)
            {
                a[i] = 1ll * a[i] * __tmpa % P;
            }
            return *this;
        }

        template <class T>
        poly operator*(const T &__tmpa) const
        {
            poly ans = *this;
            return ans *= __tmpa;
        }

        poly operator&(const poly &__tmpa) const
        {
            poly __tmpmulx = *this, __tmpmuly = __tmpa;
            int nnn = a.size(), mmm = __tmpa.a.size();
            int len = 1, L = 0;
            while (len < nnn)
            {
                len <<= 1;
                L++;
            }
            for (register int i = 1; i < len; ++i)
            {
                __Binary_reverse[i] = (__Binary_reverse[i >> 1] >> 1) | ((i & 1) << (L - 1));
            }
            __tmpmulx.resize(len);
            __tmpmuly.resize(len);
            __tmpmulx.NTT(1);
            __tmpmuly.NTT(1);
            poly res(len);
            for (register long long i = 0; i < len; ++i)
            {
                res[i] = 1ll * ((1ll * __tmpmulx[i] * __tmpmuly[i] % P) + P) % P;
            }
            res.NTT(0);
            poly ans(nnn - mmm + 1);
            for (register int __i = mmm - 1; __i < nnn; __i++)
            {
                ans[__i - mmm + 1] = res[__i];
            }
            return ans;
        }

        poly operator~() const
        {
            poly __f = *this, __tmpf;
            int cst = __f.a[0];
            int L = 0;
            poly ans;
            int dep = 1;
            ans.pb(Pre::Inv(cst));
            while (dep < (n << 1))
            {
                ++L;
                for (register int i = 1; i < (dep << 1); ++i)
                {
                    __Binary_reverse[i] = (__Binary_reverse[i >> 1] >> 1) | ((i & 1) << (L - 1));
                }
                __tmpf.a.clear();
                __tmpf = *this;
                __tmpf.resize(dep);
                __tmpf.resize(dep << 1);
                if (dep != 1)
                {
                    ans.resize(dep >> 1);
                }
                ans.resize(dep << 1);
                __tmpf.NTT(1);
                ans.NTT(1);
                for (register int i = 0; i < (dep << 1); ++i)
                {
                    ans[i] = 1ll * ans[i] * Pint::tp(2ll - 1ll * __tmpf[i] * ans[i] % P + P) % P;
                }
                ans.NTT(0);
                ans.resize(dep << 1);
                ans.tp(ans.n);
                dep <<= 1;
            }
            ans.resize(n);
            return ans;
        }

        poly operator/=(poly g)
        {
            rev(a.begin(), a.end());
            g.rev(g.a.begin(), g.a.end());
            int m = g.n;
            g.resize(n - m + 1);
            poly q = ~g * *this;
            q.rev(q.a.begin(), q.a.begin() + n - m + 1);
            return q;
        }

        poly operator/=(int g)
        {
            for (register int i = 0; i < n; ++i)
            {
                a[i] = 1ll * a[i] * Pre::Inv(g) % P;
            }
            tp(n);
            return *this;
        }

        template <class T>
        poly operator/(T &g) const
        {
            poly res = *this;
            return res /= g;
        }

        poly operator%=(poly g)
        {
            if (g.n > n)
                return *this;
            poly q = *this / g;
            q.resize(n - g.n + 1);
            return *this = *this - q * g;
        }

        poly operator%(poly g) const
        {
            poly res = *this;
            return res %= g;
        }

        poly operator<<=(int g)
        {
            poly res;
            for (register int i = 1; i <= g; ++i)
            {
                res.pb(0);
            }
            for (register int i = 0; i < n; ++i)
            {
                res.pb(a[i]);
            }
            return *this = res;
        }

        poly operator<<(int g) const
        {
            poly res = *this;
            return res <<= g;
        }

        poly operator>>=(int g)
        {
            poly res;
            for (register int i = g; i < n; ++i)
            {
                res.pb(a[i]);
            }
            return *this = res;
        }

        poly operator>>(int g) const
        {
            poly res = *this;
            return res >>= g;
        }

        poly operator^=(int g)
        {
            int __fir;
            for (register int i = 0; i < n; ++i)
            {
                if (a[i])
                {
                    __fir = i;
                    break;
                }
            }
            int con = a[__fir];
            *this >>= __fir;
            *this /= con;
            *this = Ln(*this);
            *this *= g;
            *this = Exp(*this);
            *this *= Pre::Q(con, g);
            *this <<= (__fir * g);
            return *this;
        }

        poly operator^(int g) const
        {
            poly res = *this;
            return res ^= g;
        }
    };
    class ffp
    {
    public:
        int n;
        vector<int> a;
        ffp()
        {
            n = 0;
        }
        ffp(int __n)
        {
            resize(__n);
        }

        int &operator[](int id)
        {
            return a[id];
        }

        inline void resize(int __lim_siz)
        {
            a.resize(__lim_siz);
            a.shrink_to_fit();
            n = __lim_siz;
        }

        inline void shrink()
        {
            resize(n);
        }

        inline void rev(vector<int>::iterator __st, vector<int>::iterator __en)
        {
            reverse(__st, __en);
        }

        inline void pb(int x)
        {
            a.push_back(x);
            ++n;
        }

        inline void tp(int __n)
        {
            for (register int i = 0; i < __n; ++i)
            {
                a[i] = Pint::tp(a[i]);
            }
        }

        ffp operator+=(const ffp &__tmpa)
        {
            int nm = Vmax(n, __tmpa.n);
            resize(nm);
            for (register int i = 0; i < nm; ++i)
            {
                a[i] = add((i < __tmpa.n) ? __tmpa.a[i] : 0, (i < n) ? a[i] : 0);
            }
            return *this;
        }

        ffp operator-=(const ffp &__tmpa)
        {
            int nm = Vmax(n, __tmpa.n);
            resize(nm);
            for (register int i = 0; i < nm; ++i)
            {
                a[i] = del((i < n) ? a[i] : 0, (i < __tmpa.n) ? __tmpa.a[i] : 0);
            }
            return *this;
        }

        ffp operator+(ffp __tmpa) const
        {
            ffp __tmp = *this;
            return (__tmpa += __tmp);
        }

        ffp operator-(ffp __tmpa) const
        {
            ffp __tmp = *this;
            return (__tmp -= __tmpa);
        }
    };
}

using namespace POLY;
inline poly Dx(poly &__f)
{
    poly ans;
    for (register int i = 1; i < __f.n; ++i)
    {
        ans.pb(1ll * i * __f[i] % P);
    }
    return ans;
}

inline poly Integ(poly &__f)
{
    poly ans;
    ans.pb(0);
    for (register int i = 0; i < __f.n; ++i)
    {
        ans.a.push_back(1ll * Pre::ny[i + 1] * __f.a[i] % P);
    }
    return ans;
}

inline poly UCPF::Ln(poly __F)
{
    poly dF = Dx(__F);
    dF *= ~__F;
    dF.resize(__F.n);
    dF = Integ(dF);
    dF.resize(__F.n);
    return dF;
}

inline poly UCPF::Exp(poly __f)
{
    poly ans;
    int dep = 1;
    ans.pb(1);
    poly lnf, tmpff;
    while (dep < __f.n << 1)
    {
        lnf = Ln(ans);
        lnf = __f - lnf;
        lnf.resize(dep << 1);
        ++lnf.a[0];
        ans *= lnf;
        ans.resize(dep << 1);
        dep <<= 1;
    }
    ans.resize(__f.n);
    return ans;
}

inline poly UCPF::Sqrt(poly __f)
{
    int tmp = __f[0];
    __f *= Pre::Inv(__f[0]);
    __f = Ln(__f);
    __f *= B;
    __f = Exp(__f);
    tmp = (tmp != 1) ? Quad::work(tmp, 998244353) : 1;
    __f *= tmp;
    return __f;
}

inline poly UCPF::Sin(poly __f)
{
    return (Exp(__f * _I_) - Exp(__f * (P - _I_))) * B * Pre::Inv(_I_);
}

inline poly UCPF::Cos(poly __f)
{
    return (Exp(__f * _I_) + Exp(__f * (P - _I_))) * B;
}

inline poly UCPF::Tan(poly __f)
{
    return Sin(__f) * ~Cos(__f);
}

inline poly UCPF::ArcSin(poly __f)
{
    poly __tmp = __f;
    int __n = __f.n;
    __f = __f * __f;
    __f.resize(__n);
    __f *= (P - 1);
    __f[0]++;
    __f = Sqrt(__f);
    __f = ~__f;
    __f *= Dx(__tmp);
    __f.resize(__n);
    __f = Integ(__f);
    return __f;
}

inline poly UCPF::ArcCos(poly __f)
{
    return ArcSin(__f) * (P - 1);
}

inline poly UCPF::ArcTan(poly __f)
{
    poly __tmp = __f;
    int __n = __f.n;
    __f = __f * __f;
    __f.resize(__n);
    __f[0]++;
    __f = ~__f;
    __f *= Dx(__tmp);
    __f.resize(__n);
    __f = Integ(__f);
    return __f;
}

inline poly UCPF::AND(poly __f, poly __g)
{
    __f.Fwtand(1);
    __g.Fwtand(1);
    for (register int i = 0; i < __f.n; ++i)
    {
        __f[i] = (1ll * __f[i] * __g[i] % P);
    }
    __f.Fwtand(0);
    __f.tp(__f.n);
    return __f;
}

inline poly UCPF::OR(poly __f, poly __g)
{
    __f.Fwtor(1);
    __g.Fwtor(1);
    for (register int i = 0; i < __f.n; ++i)
    {
        __f[i] = (1ll * __f[i] * __g[i] % P);
    }
    __f.Fwtor(0);
    __f.tp(__f.n);
    return __f;
}

inline poly UCPF::XOR(poly __f, poly __g)
{
    __f.Fwtxor(1);
    __g.Fwtxor(1);
    for (register int i = 0; i < __f.n; ++i)
    {
        __f[i] = 1ll * __f[i] * __g[i] % P;
    }
    __f.Fwtxor(0);
    __f.tp(__f.n);
    return __f;
}

inline poly UCPF::Pow_For_Luogu(poly __f, int __k1, int __k2, int kk)
{
    int __fir;
    for (register int i = 0; i < __f.n; ++i)
    {
        if (__f.a[i])
        {
            __fir = i;
            break;
        }
    }
    if (1ll * kk * __fir >= 1ll * __f.n)
    {
        __f.a.clear();
        __f.resize(__f.n);
        return __f;
    }
    int con = __f.a[__fir];
    __f >>= __fir;
    __f /= con;
    __f = Ln(__f);
    __f *= __k2;
    __f = Exp(__f);
    __f *= Pre::Q(con, __k1);
    __f <<= (__fir * __k2);
    return __f;
}

#ifdef TRANS
int f[N];
inline poly UCPF::Trans(poly __f, int c)
{
    poly __g(__f.n), __h(__f.n);
    f[0] = 1;
    for (int i = 1; i < __f.n; ++i)
    {
        f[i] = 1ll * f[i - 1] * i % P;
    }
    for (int i = 0; i < __f.n; ++i)
    {
        __g[i] = 1ll * __f[i] * f[i] % P;
        __h[i] = 1ll * Pre::Q(c, i) * Pre::Inv(f[i]) % P;
    }
    __g.rev(__g.a.begin(), __g.a.end());
    __h *= __g;
    __h.resize(__f.n);
    __h.rev(__h.a.begin(), __h.a.end());
    for (int i = 0; i < __h.n; ++i)
    {
        __h[i] = 1ll * __h[i] * Pre::Inv(f[i]) % P;
    }
    return __h;
}
#endif

void PIO::pin(poly &f, int __n)
{
    for (register int i = 0; i < __n; ++i)
    {
        int x;
        read(x);
        f.pb(x);
    }
}

void PIO::ppri(poly __f, int __n)
{
    __f.tp(__f.n);
    __f.resize(__n);
    for (register int i = 0; i < __n; ++i)
    {
        print(__f[i]);
        pc(' ');
    }
    pc('\n');
}

void PIO::pin(ffp &f, int __n)
{
    for (register int i = 0; i < __n; ++i)
    {
        int x;
        read(x);
        f.pb(x);
    }
}

void PIO::ppri(ffp __f, int __n)
{
    __f.tp(__f.n);
    __f.resize(__n);
    for (register int i = 0; i < __n; ++i)
    {
        print(__f[i]);
        pc(' ');
    }
    pc('\n');
}

#define PMPE PMPE
#ifdef PMPE

namespace Multipoint_Evel
{
    poly t[N << 1];
    int an[N];
    void build(int l, int r, int p, int *a)
    {
        if (l == r)
        {
            t[p].pb(1ll);
            t[p].pb(P - a[l]);
            return;
        }
        int mid = (l + r) >> 1;
        build(l, mid, p << 1, a);
        build(mid + 1, r, p << 1 | 1, a);
        t[p] = t[p << 1] * t[p << 1 | 1];
    }

    void multipoint_evel(int l, int r, int p, poly __f)
    {

        if (l == r)
        {
            an[l] = __f[0];
            return;
        }
        poly r1 = __f & t[p << 1 | 1], r2 = __f & t[p << 1];
        int mid = (l + r) >> 1;
        multipoint_evel(l, mid, p << 1, r1);
        multipoint_evel(mid + 1, r, p << 1 | 1, r2);
    }

    inline void PMPE(int *val, poly __f, int mm, int *ans)
    {
        int con = __f[0];
        int lm = mm, nn = __f.n;
        if (nn <= mm - 1)
        {
            __f.resize(mm + 1);
        }
        else if (nn >= mm)
        {
            mm = Vmax(mm, nn - 1);
        }
        build(1, mm, 1, val);
        __f.rev(__f.a.begin(), __f.a.end());
        __f.resize(nn << 1);
        poly tmp = __f * ~t[1];
        tmp.resize(mm);
        multipoint_evel(1, mm, 1, tmp);
        for (int i = 1; i <= lm; i++)
        {
            ans[i] = ((1ll * val[i] * an[i] % P + con) % P + P) % P;
        }
    }
}
using Multipoint_Evel::PMPE;
#endif
#ifdef FPI
namespace Fast_Interpolation
{
    poly q[N << 2], ret[N << 2];
    int a[N];
    void build1(int l, int r, int p, int *x)
    {
        if (l == r)
        {
            q[p].pb(P - x[l]);
            q[p].pb(1ll);
            return;
        }
        int mid = (l + r) >> 1;
        build1(l, mid, p << 1, x);
        build1(mid + 1, r, p << 1 | 1, x);
        q[p] = q[p << 1] * q[p << 1 | 1];
    }

    void Inv_evel(int l, int r, int p)
    {
        if (l == r)
        {
            ret[p].pb(a[l]);
            return;
        }
        int mid = (l + r) >> 1;
        Inv_evel(l, mid, p << 1);
        Inv_evel(mid + 1, r, p << 1 | 1);
        ret[p] = ret[p << 1] * q[p << 1 | 1] + ret[p << 1 | 1] * q[p << 1];
    }

    inline poly PFI(int *x, int *y, int __n)
    {
        build1(1, __n, 1, x);
        poly F = Dx(q[1]);
        int con = F[0];
        F.resize(__n + 1);
        Multipoint_Evel::build(1, __n, 1, x);
        F.rev(F.a.begin(), F.a.end());
        F.resize(__n << 1);
        poly tmp = F * ~Multipoint_Evel::t[1];

        tmp.resize(__n);
        Multipoint_Evel::multipoint_evel(1, __n, 1, tmp);
        for (int i = 1; i <= __n; ++i)
        {
            a[i] = 1ll * y[i] * Pre::Inv(((1ll * x[i] * Multipoint_Evel::an[i] % P + con) % P + P) % P) % P;
        }
        Inv_evel(1, __n, 1);
        return ret[1];
    }
}
using Fast_Interpolation::PFI;

#endif

namespace Chirp_Z
{
    inline int Ci2(int &i)
    {
        register long long k = 1ll * i * (i - 1) / 2;
        k %= (P - 1);
        return k;
    }
    inline poly Chirp_Z(poly __f, int c, int m)
    { // You can't pass P6800 by this code,because this PolyBoard is limit at 100000
        poly res;
        int __n = __f.n;
        poly R(__n + m), G(__n + m);
        __f.resize(__n + m - 1);
        for (register int i = 0; i < __n + m - 1; ++i)
        {
            R[__n + m - 2 - i] = 1ll * Pre::Q(c, Ci2(i)) % P;
            G[i] = 1ll * Pre::Q(c, (P - 1ll - Ci2(i)) % P) * __f[i] % P;
        }
        poly ans = R * G;
        for (int i = 0; i < m; i++)
        {
            res.pb(1ll * ans[__n + m - 2 - i] * Pre::Q(c, (P - 1ll - Ci2(i)) % P) % P);
        }
        return res;
    }
}

#ifdef STIRLING

namespace Stirling
{
    // In fact, if you want to pass the problems in luogu,you have to change the modulo.
    inline poly row(int __n)
    {
        poly f(__n + 1), _f(__n + 1);
        f[0] = _f[0] = 1;
        for (int i = 1; i <= __n; i++)
        {
            f[i] = 1ll * f[i - 1] * Pre::Inv(i) % P;
            _f[i] = f[i];
        }
        for (int i = 0; i <= __n; i++)
        {
            f[i] = 1ll * f[i] * Pre::Q(i, __n) % P;
            _f[i] = 1ll * _f[i] * ((i & 1) ? (P - 1) : 1) % P;
        }
        f *= _f;
        f.tp(f.n);
        f.resize(__n + 1);
        return f;
    }

    poly t[N << 1];
    void Dev_Mul(int l, int r, int p)
    {
        if (l == r)
        {
            t[p].pb(1ll);
            t[p].pb(P - l);
            return;
        }
        int mid = (l + r) >> 1;
        Dev_Mul(l, mid, p << 1);
        Dev_Mul(mid + 1, r, p << 1 | 1);
        t[p] = t[p << 1] * t[p << 1 | 1];
    }

    inline poly column(int __k, int __n)
    {
        poly ans;
        if (__k > __n)
        {
            ans.resize(__n + 1);
            return ans;
        }
        Dev_Mul(1, __k, 1);
        t[1].resize(__n + 1);
        ans = (~t[1]) << __k;
        ans.resize(__n + 1);
        return ans;
    }
}

namespace Staling
{
    int f[N];

    inline poly row(int __n)
    {
        poly ans;
        if (__n == 1)
        {
            ans.pb(0);
            ans.pb(1);
            return ans;
        }
        poly lst = row(__n >> 1);
        poly lst_c = Trans(lst, __n >> 1);
        lst *= lst_c;
        if (((__n >> 1) << 1) == __n)
        {
            return lst;
        }
        else
        {
            lst = lst * (__n - 1) + (lst << 1);
            return lst;
        }
    }

    inline poly column(int __k, int __n)
    {
        poly __f;
        __f.pb(0);
        f[0] = 1;
        for (int i = 1; i <= __n; i++)
        {
            f[i] = 1ll * f[i - 1] * i % P;
            __f.pb(Pre::Inv(i));
        }
        __f = Exp(Ln(__f >> 1) * __k) << __k;
        __f.resize(__n + 1);
        __f /= f[__k];
        for (int i = 0; i <= __n; i++)
        {
            __f[i] = 1ll * __f[i] * f[i] % P;
        }
        return __f;
    }
}

#endif

// #define FFP FFP
// #define PTFFP PTFFP
#ifdef FFP

namespace FFP
{
    int f[N], a[N], r[N];

#ifdef PTFFP

    poly t[N << 1], ans[N << 1];
    void Dev_Mul(int l, int r, int p, poly &a)
    {
        if (l == r)
        {
            t[p].pb(P - l);
            t[p].pb(1);
            ans[p].pb(a[l]);
            return;
        }
        int mid = (l + r) >> 1;
        Dev_Mul(l, mid, p << 1, a);
        Dev_Mul(mid + 1, r, p << 1 | 1, a);
        t[p] = t[p << 1] * t[p << 1 | 1];
        ans[p] = ans[p << 1] + t[p << 1] * ans[p << 1 | 1];
    }

    inline ffp PTFFP(poly __f)
    {
        Dev_Mul(0, __f.n - 1, 1, __f);
        ffp __res;
        for (int i = 0; i < ans[1].n; ++i)
            __res.pb(ans[1][i]);
        return __res;
    }

#endif
    inline poly FFPTP(ffp ___f)
    {
        poly __f;
        for (int i = 0; i < ___f.n; ++i)
        {
            __f.pb(___f[i]);
        }
        poly __g, __h;
        f[0] = 1;
        for (int i = 1; i <= __f.n + 1; i++)
        {
            f[i] = 1ll * f[i - 1] * Pre::Inv(i) % P;
            a[i] = i - 1;
        }
        PMPE(a, __f, __f.n + 1, r);
        for (int i = 0; i <= __f.n; i++)
        {
            __g.pb(1ll * r[i + 1] * f[i] % P);
            __h.pb(1ll * ((i & 1) ? (P - 1) : 1) * f[i] % P);
        }
        __g *= __h;
        __g.resize(__f.n + 1);

        return __g;
    }
}
#endif

namespace CCHLR
{
    poly gamma;
    poly __calc(int b)
    {
        poly res, a;
        res.pb(1);
        a.pb(0), a.pb(1);
        while (b)
        {
            if (b & 1)
            {
                res = res * a;
                res %= gamma;
                res.resize(gamma.n);
            }
            a = a * a;
            a %= gamma;
            a.resize(gamma.n);
            b >>= 1;
        }
        return res;
    }
    int Fiduccia(int *c, int *a, int k, int n)
    {

        for (register int i = 1; i <= k; ++i)
        {
            gamma.pb((P - c[k - i + 1]) % P);
        }
        gamma.pb(1);
        poly tmp = __calc(n);
        long long ans = 0;
        for (int i = 0; i < k; ++i)
        {
            ans += 1ll * tmp[i] * a[i] % P;
            ans %= P;
        }
        return tp(ans);
    }

    int Bostan_Mori(poly __f, poly __g, int __n)
    {
        while (__n)
        {
            poly t = __g;
            for (int i = 1; i < __g.n; i += 2)
            {
                t[i] = ((P - __g[i]) % P + P) % P;
            }
            __f *= t;
            __g *= t;
            int i = (__n & 1);
            for (; i < __f.n; i += 2)
            {
                __f[i >> 1] = __f[i];
            }
            __f.resize(i >> 1);
            for (i = 0; i < __g.n; i += 2)
            {
                __g[i >> 1] = __g[i];
            }
            __g.resize(i >> 1);
            __n >>= 1;
        }
        if (__f.n == 0)
            return 0;
        return 1ll * __f[0] * Pre::Inv(__g[0]) % P;
    }

    int Bostan_Mori_CCHLR(int *c, int *a, int k, int n)
    {
        gamma.pb(1);
        for (register int i = 1; i <= k; ++i)
        {
            gamma.pb(((P - c[i]) % P + P) % P);
        }
        poly __tmp;
        for (int i = 0; i < k; i++)
        {
            __tmp.pb((a[i] % P + P) % P);
        }
        poly res = __tmp * gamma;
        res.resize(k);
        return Bostan_Mori(res, gamma, n);
    }
}
using CCHLR::Bostan_Mori;
using CCHLR::Bostan_Mori_CCHLR;
using CCHLR::Fiduccia;

#define CCNHLR CCNHLR
#ifdef CCNHLR
namespace CCNHLR
{
    poly gamma;
    int __y[N], __x[N];
    poly __calc(int b)
    {
        poly res, a;
        res.pb(1);
        a.pb(0), a.pb(1);
        while (b)
        {
            if (b & 1)
            {
                res = res * a;
                res %= gamma;
                res.resize(gamma.n);
            }
            a = a * a;
            a %= gamma;
            a.resize(gamma.n);
            b >>= 1;
        }
        return res;
    }
    poly __ksm(int b)
    {
        poly res, a;
        res.pb(1);
        a.pb(P - 1), a.pb(1);
        while (b)
        {
            if (b & 1)
            {
                res *= a;
            }
            a *= a;
            b >>= 1;
        }
        return res;
    }
    inline int Bostam_Mori_CCNHLR(int *c, int *a, int k, int n, poly p)
    {// Very Slow
        gamma.pb(1);
        for (register int i = 1; i <= k; ++i)
        {
            gamma.pb(((P - c[i]) % P + P) % P);
        }
        poly __tmp, __tmp2;
        for (int i = 0; i < k; i++)
        {
            __tmp.pb(a[i]);
        }
        __tmp2.pb(0);
        for (int i = 1; i <= k; i++)
            __tmp2.pb(c[i]);
        poly __p = __tmp * __tmp2;
        __p.resize(k);
        for (int i = 0; i < k; i++)
        {
            __p[i] = P - __p[i];
            __p[i] =(__p[i] + a[i]) % P;
        }
        __p.resize(k + p.n + 1);
        for (int i = 1; i <= p.n + 1; i++)
            __x[i] = i + k - 1;
        PMPE(__x, p, p.n + 1, __y); // P __p
        for (int i = 1; i <= p.n + 1; i++)
        {
            __p[i + k - 1] = __y[i];
        }
        gamma.resize(p.n + k + 1);
        poly __a = __p * (~gamma);

        __a.resize(p.n + k + 1);
        gamma.resize(0);
        for (register int i = 1; i <= k; ++i)
        {
            gamma.pb((P - c[k - i + 1]) % P);
        }
        gamma.pb(1);
        gamma *= __ksm(p.n + 1);
        poly tmp = __calc(n);
        long long ans = 0;
        for (int i = 0; i < p.n + k + 1; ++i)
        {
            ans += 1ll * tmp[i] * __a[i] % P;
            ans %= P;
        }
        return tp(ans);
    }
}
using CCNHLR::Bostam_Mori_CCNHLR;
#endif

int n, k, m;
int a[N], c[N];
poly f;

inline void work()
{
    read(n), read(m), read(k);
    for (int i = 0; i < k; i++)
        read(a[i]);
    for (int i = 1; i <= k; i++)
        read(c[i]);
    pin(f, m + 1);
    print(Bostam_Mori_CCNHLR(c, a, k, n, f));
}

signed main()
{
    ios::sync_with_stdio(false);
    cin.tie(0), cout.tie(0);
    Pre::initYG();
    work();
    return 0;
}