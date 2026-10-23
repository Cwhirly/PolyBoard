#include "poly.cpp"

namespace MATRIX {
    class matrix {
        vector<vector<int>> a;
        int n, m;

        int __rref(int __lim) {
            int __r = 1;
            for (int k = 1; k <= __lim && __r <= n; ++k) {
                int __p = 0, __iv = 0;
                for (int i = __r; i <= n; ++i)
                    if (a[i][k] && (__iv = Pre::Inv(a[i][k]))) {
                        __p = i;
                        break;
                    }
                if (!__p)
                    continue;
                swap(a[__r], a[__p]);
                for (int j = k; j <= m; ++j)
                    a[__r][j] = 1ll * a[__r][j] * __iv % P;
                for (int i = 1; i <= n; ++i) {
                    if (i == __r || !a[i][k])
                        continue;
                    int __x = a[i][k];
                    for (int j = k; j <= m; ++j)
                        a[i][j] = (a[i][j] + 1ll * __x * (P - a[__r][j])) % P;
                }
                ++__r;
            }
            return __r - 1;
        }

    public:
        matrix() { n = m = 0; }
        matrix(int __n, int __m) { resize(__n, __m); }

        void resize(int __n, int __m) {
            a.resize(__n + 1);
            for (int i = 1; i <= __n; ++i)
                a[i].resize(__m + 1);
            n = __n, m = __m;
        }

        int row() const { return n; }

        int col() const { return m; }

        vector<int> &operator[](int __i) { return a[__i]; }

        const vector<int> &operator[](int __i) const { return a[__i]; }

        int &operator()(int __i, int __j) { return a[__i][__j]; }

        const int &operator()(int __i, int __j) const { return a[__i][__j]; }

        matrix &operator+=(const matrix &__t) {
            for (int i = 1; i <= n; ++i)
                for (int j = 1; j <= m; ++j)
                    addt(a[i][j], __t.a[i][j]);
            return *this;
        }

        matrix &operator-=(const matrix &__t) {
            for (int i = 1; i <= n; ++i)
                for (int j = 1; j <= m; ++j)
                    delt(a[i][j], __t.a[i][j]);
            return *this;
        }

        matrix &operator*=(const matrix &__t) {
            matrix __ans(n, __t.m);
            for (int i = 1; i <= n; ++i) {
                vector<int> &__ri = __ans.a[i];
                for (int k = 1; k <= m; ++k) {
                    int __x = a[i][k];
                    if (!__x)
                        continue;
                    const vector<int> &__tk = __t.a[k];
                    for (int j = 1; j <= __t.m; ++j)
                        __ri[j] = (__ri[j] + 1ll * __x * __tk[j]) % P;
                }
            }
            a.swap(__ans.a);
            n = __ans.n, m = __ans.m;
            return *this;
        }

        matrix operator+(const matrix &__t) const {
            matrix __res = *this;
            return __res += __t;
        }

        matrix operator-(const matrix &__t) const {
            matrix __res = *this;
            return __res -= __t;
        }

        matrix operator-() const {
            matrix __res(n, m);
            for (int i = 1; i <= n; ++i)
                for (int j = 1; j <= m; ++j)
                    __res.a[i][j] = a[i][j] ? P - a[i][j] : 0;
            return __res;
        }

        matrix operator*(const matrix &__t) const {
            matrix __res = *this;
            return __res *= __t;
        }

        bool gauss() { return __rref(m) == n; }

        bool inv(matrix &__res) const {
            if (n != m) {
                __res = matrix(n, m);
                return false;
            }
            matrix __t(n, m << 1);
            for (int i = 1; i <= n; ++i) {
                for (int j = 1; j <= m; ++j)
                    __t.a[i][j] = a[i][j];
                __t.a[i][m + i] = 1;
            }
            if (__t.__rref(n) != n) {
                __res = matrix(n, m);
                return false;
            }
            __res.resize(n, m);
            for (int i = 1; i <= n; ++i)
                for (int j = 1; j <= m; ++j)
                    __res.a[i][j] = __t.a[i][m + j];
            return true;
        }

        int det() const {
            matrix __t = *this;
            int __res = 1;
            for (int k = 1; k < n; ++k) {
                int __p = 0, __iv = 0;
                for (int i = k; i <= n; ++i)
                    if (__t.a[i][k] && (__iv = Pre::Inv(__t.a[i][k]))) {
                        __p = i;
                        break;
                    }
                if (!__p)
                    return 0;
                if (__p != k) {
                    swap(__t.a[k], __t.a[__p]);
                    __res = P - __res;
                }
                for (int i = k + 1; i <= n; ++i) {
                    if (!__t.a[i][k])
                        continue;
                    int __x = 1ll * __t.a[i][k] * __iv % P;
                    for (int j = k; j <= n; ++j)
                        __t.a[i][j] = (__t.a[i][j] + 1ll * __x * (P - __t.a[k][j])) % P;
                }
            }
            for (int i = 1; i <= n; ++i)
                __res = 1ll * __res * __t.a[i][i] % P;
            return __res;
        }

        matrix operator~() const {
            matrix __res;
            if (!inv(__res))
                __res = matrix(n, m);
            return __res;
        }
    };

    namespace MIO {
        matrix pin(int __n, int __m) {
            matrix __f(__n, __m);
            for (int i = 1; i <= __n; ++i)
                for (int j = 1; j <= __m; ++j)
                    read(__f(i, j));
            return __f;
        }

        void ppri(const matrix &__f, int __n, int __m) {
            for (int i = 1; i <= __n; ++i) {
                for (int j = 1; j <= __m; ++j) {
                    print((i <= __f.row() && j <= __f.col()) ? Pint::tp(__f(i, j)) : 0);
                    pc(' ');
                }
                pc('\n');
            }
        }
    }
}


#if __INCLUDE_LEVEL__ == 0

int n;

inline void work()
{
    read(n);
    MATRIX::matrix f = MATRIX::MIO::pin(n, n), g;
    if (!f.inv(g)) {
        print(-1);
        pc('\n');
        return;
    }
    MATRIX::MIO::ppri(g, n, n);
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(0), cout.tie(0);
    Pre::initYG();
    work();
    return 0;
}

#endif
