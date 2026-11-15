#define POLYBOARD_P 998244353
#define POLYBOARD_ATN 1
#include "../matrix.cpp"
using namespace std;

template <class T> void out(const T &a, int n) {
    for (int i = 0; i < n; ++i)
        cout << (i ? " " : "") << a[i];
    cout << "\n";
}


static bool same(const poly &a, const poly &b, int n) {
    for (int i = 0; i < n; ++i) {
        int x = i < a.size() ? a[i] : 0;
        int y = i < b.size() ? b[i] : 0;
        if (x != y)
            return false;
    }
    return true;
}

static void rnd(poly &f, mt19937 &rng) {
    for (int i = 0; i < f.size(); ++i)
        f[i] = rng() % P;
}

static bool prop(string op, int n, unsigned seed) {
    mt19937 rng(seed);
    if (op == "sqrt") {
        poly f(n);
        rnd(f, rng);
        f[0] = 1;
        poly h = UCPF::Sqrt(f);
        poly p = h * h;
        return same(p, f, n);
    }
    if (op == "trig") {
        poly f(n);
        rnd(f, rng);
        f[0] = 0;
        poly s = UCPF::Sin(f), c = UCPF::Cos(f);
        poly p = s * s + c * c;
        poly one(1);
        one[0] = 1;
        return same(p, one, n);
    }
    if (op == "trans") {
        poly f(n);
        rnd(f, rng);
        int c = rng() % P;
        poly h = UCPF::Trans(f, c);
        poly g(n);
        vector<int> pw(n);
        pw[0] = 1;
        for (int i = 1; i < n; ++i)
            pw[i] = 1ll * pw[i - 1] * c % P;
        for (int k = 0; k < n; ++k) {
            long long z = 0;
            for (int i = k; i < n; ++i) {
                int b = 1ll * Pre::fac[i] * Pre::ifac[k] % P * Pre::ifac[i - k] % P;
                z += 1ll * f[i] * b % P * pw[i - k] % P;
                if (z >= (1ll << 61))
                    z %= P;
            }
            g[k] = z % P;
        }
        return same(h, g, n);
    }
    if (op == "div") {
        poly f(2 * n), g(n + 1);
        rnd(f, rng);
        rnd(g, rng);
        g[0] = 1;
        poly q = f / g;
        poly p = q * g + (f % g);
        return same(p, f, f.size());
    }
    if (op == "bit") {
        int len = 1 << 3;
        poly f(len), g(len);
        rnd(f, rng);
        rnd(g, rng);
        for (int t = 0; t < 3; ++t) {
            poly h = t == 0 ? UCPF::AND(f, g)
                     : t == 1 ? UCPF::OR(f, g)
                              : UCPF::XOR(f, g);
            poly e(len);
            for (int i = 0; i < len; ++i)
                for (int j = 0; j < len; ++j) {
                    int k = t == 0 ? (i & j) : t == 1 ? (i | j) : (i ^ j);
                    e[k] = mod(e[k] + 1ll * f[i] * g[j]);
                }
            if (!same(h, e, len))
                return false;
        }
        return true;
    }
    if (op == "subset") {
        int k = 3, len = 1 << k;
        sps f(len), g(len);
        for (int i = 0; i < len; ++i)
            f[i] = rng() % P, g[i] = rng() % P;
        sps h = f * g;
        sps e(len);
        for (int s = 0; s < len; ++s)
            for (int t = s;; t = (t - 1) & s) {
                e[s] = mod(e[s] + 1ll * f[t] * g[s ^ t]);
                if (!t)
                    break;
            }
        for (int i = 0; i < len; ++i)
            if (h[i] != e[i])
                return false;
        return true;
    }
    if (op == "shift") {
        poly f(n + 1);
        rnd(f, rng);
        poly y(n + 1);
        for (int i = 0; i <= n; ++i) {
            long long z = 0;
            for (int j = n; j >= 0; --j)
                z = (z * i + f[j]) % P;
            y[i] = z;
        }
        int m = n + 1 + rng() % 100;
        poly h = UCPF::Shift_Values(y, m);
        for (int i = 0; i <= n; ++i) {
            long long z = 0;
            for (int j = n; j >= 0; --j)
                z = (z * (m + i) + f[j]) % P;
            if (h[i] != z)
                return false;
        }
        return true;
    }
    if (op == "interp") {
        poly f(n);
        rnd(f, rng);
        vector<int> x(n);
        for (int i = 0; i < n; ++i)
            x[i] = i;
        poly y = PMPE(f, x);
        vector<int> yv(n);
        for (int i = 0; i < n; ++i)
            yv[i] = y[i];
        poly h = PFI(x, yv);
        return same(h, f, n);
    }
    return false;
}

static unsigned long long bench(string op, int n, unsigned seed, int rep) {
    mt19937 rng(seed);
    poly f(n), g(n);
    rnd(f, rng);
    rnd(g, rng);
    f[0] = 0;
    g[0] = 1;
    poly h;
    auto st = chrono::steady_clock::now();
    for (int t = 0; t < rep; ++t) {
        if (op == "conv")
            h = f * g;
        else if (op == "inv")
            h = ~g;
        else if (op == "ln")
            h = UCPF::Ln(g);
        else if (op == "exp")
            h = UCPF::Exp(f);
        else if (op == "compose")
            h = UCPF::Compose(f, g);
        else
            return 0;
    }
    auto en = chrono::steady_clock::now();
    unsigned long long z = 0;
    for (int i = 0; i < h.size(); ++i)
        z = z * 131 + h[i];
    fprintf(stderr, "%s n=%d rep=%d %.3fms hash=%llu\n", op.c_str(), n, rep,
            chrono::duration<double, milli>(en - st).count(), z);
    return z;
}

int main(int argc, char **argv) {
    if (argc < 2)
        return 2;
    string mode = argv[1];
    Pre::initYG();
    if (mode == "conv") {
        int n, m;
        cin >> n >> m;
        poly f = PIO::pin(n), g = PIO::pin(m);
        out(f * g, n + m - 1);
    } else if (mode == "inv") {
        int n;
        cin >> n;
        poly f = PIO::pin(n);
        out(~f, n);
    } else if (mode == "ln") {
        int n;
        cin >> n;
        poly f = PIO::pin(n);
        out(UCPF::Ln(f), n);
    } else if (mode == "exp") {
        int n;
        cin >> n;
        poly f = PIO::pin(n);
        out(UCPF::Exp(f), n);
    } else if (mode == "compose") {
        int n, m;
        cin >> n >> m;
        poly f = PIO::pin(n), g = PIO::pin(m);
        out(UCPF::Compose(f, g), n);
    } else if (mode == "high") {
        int n, m;
        cin >> n >> m;
        poly f = PIO::pin(n), g = PIO::pin(m);
        out(f & g, n - m + 1);
    } else if (mode == "spsinv") {
        int k;
        cin >> k;
        int n = 1 << k;
        sps f(n);
        for (int i = 0; i < n; ++i)
            cin >> f[i];
        out(~f, n);
    } else if (mode == "spsexp") {
        int k;
        cin >> k;
        int n = 1 << k;
        sps f(n);
        for (int i = 0; i < n; ++i)
            cin >> f[i];
        out(UCSF::Exp(f), n);
    } else if (mode == "spsln") {
        int k;
        cin >> k;
        int n = 1 << k;
        sps f(n);
        for (int i = 0; i < n; ++i)
            cin >> f[i];
        out(UCSF::Ln(f), n);
    } else if (mode == "matmul") {
        int n;
        cin >> n;
        MATRIX::matrix a(n, n), b(n, n);
        for (int i = 1; i <= n; ++i)
            for (int j = 1; j <= n; ++j)
                cin >> a(i, j);
        for (int i = 1; i <= n; ++i)
            for (int j = 1; j <= n; ++j)
                cin >> b(i, j);
        auto c = a * b;
        for (int i = 1; i <= n; ++i)
            for (int j = 1; j <= n; ++j)
                cout << c(i, j) << (j == n ? "\n" : " ");
    } else if (mode == "matinv") {
        int n;
        cin >> n;
        MATRIX::matrix a(n, n), b;
        for (int i = 1; i <= n; ++i)
            for (int j = 1; j <= n; ++j)
                cin >> a(i, j);
        if (!a.inv(b)) {
            cout << "-1\n";
        } else {
            for (int i = 1; i <= n; ++i)
                for (int j = 1; j <= n; ++j)
                    cout << b(i, j) << (j == n ? "\n" : " ");
        }
    } else if (mode == "matdet") {
        int n;
        cin >> n;
        MATRIX::matrix a(n, n);
        for (int i = 1; i <= n; ++i)
            for (int j = 1; j <= n; ++j)
                cin >> a(i, j);
        cout << a.det() << "\n";
    } else if (mode == "shift") {
        int n, m;
        cin >> n >> m;
        poly y = PIO::pin(n + 1);
        out(UCPF::Shift_Values(y, m), n + 1);
    } else if (mode == "badinv") {
        int n;
        cin >> n;
        poly f = PIO::pin(n);
        if (f[0] == 0)
            return 3;
        out(~f, n);
    } else if (mode == "badln") {
        int n;
        cin >> n;
        poly f = PIO::pin(n);
        if (f[0] == 0)
            return 3;
        out(UCPF::Ln(f), n);
    } else if (mode == "badexp") {
        int n;
        cin >> n;
        poly f = PIO::pin(n);
        if (f[0] != 0)
            return 3;
        out(UCPF::Exp(f), n);
    } else if (mode == "badmat") {
        int n;
        cin >> n;
        MATRIX::matrix a(n, n);
        for (int i = 1; i <= n; ++i)
            for (int j = 1; j <= n; ++j)
                cin >> a(i, j);
        if (a.det() == 0)
            return 3;
        cout << a.det() << "\n";
    } else if (mode == "prop") {
        string op;
        int n;
        unsigned seed;
        cin >> op >> n >> seed;
        return prop(op, n, seed) ? 0 : 1;
    } else if (mode == "bench") {
        string op;
        int n, rep;
        unsigned seed;
        cin >> op >> n >> seed >> rep;
        return bench(op, n, seed, rep) ? 0 : 1;
    } else {
        return 2;
    }
    return 0;
}
