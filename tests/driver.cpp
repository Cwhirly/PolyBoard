#define POLYBOARD_P 998244353
#define POLYBOARD_ATN 1
#include "../matrix.cpp"
using namespace std;

template <class T> void out(const T &a, int n) {
    for (int i = 0; i < n; ++i)
        cout << (i ? " " : "") << a[i];
    cout << "\n";
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
    } else {
        return 2;
    }
    return 0;
}
