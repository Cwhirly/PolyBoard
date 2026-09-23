# PolyBoard 1.17.2

单文件、OI 向的多项式模板，模数固定 `998244353`，自带 `NTT`、形式幂级数全家桶、多点求值 / 快速插值、位运算卷积、斯特林数、下降幂转换、常系数线性递推等。

`1.17.1 = 1.17.0 的全部功能 + 封装修整 + 常数优化 + bug 修复`；
`1.17.2 = 1.17.1 + 集合幂级数 sps 类 + 非齐次线性递推重写`；
`2.1.1 = 矩阵类（基本运算 / 高斯消元 / 行列式 / 求逆），直接 include 1_17_2.cpp 复用上面全部功能`。

- 主文件（工作区只留这两份）：`1_17_2.cpp`（不带矩阵）、`2_1_1.cpp`（带矩阵）
- 历史版本：都在 git tag 里，见 `GIT.md`（例如 `git show 1.16.0:1_16_0.cpp`）
- 更新记录：`CHANGELOG.md`
- 语言标准：`C++14` 及以上（推荐 `-std=c++17`），需要 `bits/stdc++.h`，即 GCC。
- 全文件无注释、函数体内无空行；说明都在这份 README 里。

---

## 1. 快速上手

### 1.1 编译

```bash
g++ 1_17_2.cpp -o main -std=c++17 -O2     # 不带矩阵
g++ 2_1_1.cpp -o main -std=c++17 -O2      # 带矩阵
```

**没有开关宏，所有功能永远可用**（`1.17.0` 里的 `#ifdef TRANS / PMPE / FPI / STIRLING / FFP / PTFFP / CCNHLR` 已全部删除：那些 `#ifdef` 原本只是为了屏蔽巨大的全局 scratch 数组，现在这些数组都收进了函数内部 / `struct`，留着宏只会让编译参数更容易写错）。编译时间约 2 s。

### 1.2 只改两个数：`N` 和 `P`

```cpp
namespace CONSTS {
    const constexpr int P = 998244353, N = 600005, ...;
}
```

- `P`：模数。程序会自动算出原根 `Y`（编译期搜索最原根）、`I = 1/Y`、`B = 1/2`、`_I_`（四次单位根）；`static_assert` 会检查 `P ≡ 1 (mod 4)`。
- `N`：多项式长度上限（用户数组、阶乘表等的容量）。**NTT 的最大长度 `NL` 也由它推出**：`LG = ceil(log2(2N-1))`、`NL = 2^LG`（两个长度为 `N` 的多项式相乘也放得下）。`static_assert` 会检查 `2^LG | P-1`。

其余常量（`Y`、`I`、`B`、`_I_`、`LG`、`NL`）都不用管，全部由 `P`、`N` 在编译期算出。换模数 / 换规模只动这两个数。

实测：

- 只把 `N` 改成 `100005` → `LG/NL` 自动变成 `18 / 262144`，`90000 × 90000` 的卷积正确；
- 只把 `P` 改成 `167772161`（`5·2^25+1`）→ 原根 `Y = 3`、`I = 1/3`、`_I_` 全部重新算出，该模数下的卷积与 `Exp/Ln` 正确。

### 1.3 最小骨架

```cpp
void work() {
    int n, m;
    read(n), read(m);
    poly f = pin(n + 1), g = pin(m + 1);
    ppri(f * g, n + m + 1);
}

signed main() {
    Pre::initYG();
    work();
    return 0;
}
```

**`main` 里必须先调用 `Pre::initYG()`**，它负责初始化 NTT 单位根表、`1/i` 表、阶乘 / 阶乘逆表；之后才能调用任何多项式接口。

文件末尾自带的 `work()` / `main()` 就是一个求卷积的程序，正式使用时把它替换成题目要求即可。

### 1.4 三步走

1. 想改规模 / 模数就改开头的 `N` / `P`，其它什么都别动；
2. 在 `main` 顶部 `Pre::initYG();`；
3. 用 `pin(n)` / `read` 读入，用 `poly` 的运算符或 `UCPF` 里的函数算，用 `ppri` / `print` 输出。

---

## 2. 容量与限制

- `P`：模数。原根 `Y`、`B = 1/2`、`I = 1/Y`、`_I_`（四次单位根）都是编译期算出来的；`static_assert` 检查 `P ≡ 1 (mod 4)`。
- `N`：唯一的规模常量（多项式长度、`Pre::ny/fac/ifac` 表、用户数组都用它）。
- `NL`：NTT 能做的最大长度，由 `N` 推出：`LG = ceil(log2(2N-1))`、`NL = 2^LG`。默认 `N = 600005` 时 `LG = 21`、`NL = 2^21 = 2097152` —— 两个长度 600005 的多项式相乘（结果 1200009 位）也放得下。`static_assert` 检查 `2^LG | P - 1`。

各接口的规模上限（`N = 600005` 时实测）：

| 接口 | 上限 | 来源 |
| --- | --- | --- |
| `f * g`、`f & g` | `\|f\| + \|g\| - 1 ≤ NL` | NTT 上限 |
| `~f`、`f / g`、`f % g` | `n ≤ NL / 2`（2^20） | Newton 迭代里 NTT 长度 `2·dep` |
| `Chirp_Z` | `n + m ≤ NL / 2 + 1`（约 2^20） | 内部卷积长 `2(n+m)-3` |
| `Dx`、`Integ`、`Ln`、`Trans`、`Stirling`、`Staling`、`FFP` | `n ≤ N` | `Pre::ny` / `fac` / `ifac` 只到 `N` |
| `Exp`、`Sqrt`、`Sin`、`Cos`、`Tan`、`ArcSin`、`ArcCos`、`ArcTan` | `n ≤ NL / 4`（2^19） | 算法内部需要 `4·dep ≤ NL` |
| `PMPE` / `PFI` | 点数 `m ≤ (NL - 2) / 3`，且段树占 `4m` 个多项式 | 内部卷积与内存 |
| `AND` / `OR` / `XOR` / `Subset` | 长度是 2 的幂 | FWT |

另外一条通用约束：凡是内部要用 `Pre::ny` / `fac` / `ifac` 的接口（`Ln`、`Integ`、`Sqrt`、`Arc*`、`Trans`、`Stirling`、`Staling`…），都还要求 `n ≤ N`。

超出上限是未定义行为（NTT 会越界访问旋转因子表），不会报错也不会自动截断，**务必自己核对规模**。需要更大的规模时只改 `N`（上限由 `2^LG | P-1` 决定，`P = 998244353` 时 `LG ≤ 23`）。

表的内存是随 `N` 走的：`gp`/`igp` 各 `NL` 个 int（默认各 8 MiB，只有被用到的那一段会被真正访问），`ny` 与 `fac`/`ifac` 各 `N` 个 int。`N` 减小，内存同步减小。

## 3. API 手册

### 3.1 输入输出 `IO`

```cpp
template <class I> bool read(I &x);
template <class I> void print(I x);
void pc(char c);
void pstr(const char *s);
void gc(char &c);
void gstr(char *s);
void flush();
```

- `read` 返回是否成功读到（未遇 EOF 为 `true`）；
- `print` 支持负数，不输出多余空格；
- `pc` 逐字符输出；`pstr` 输出到 `'\0'` 为止（不含结尾符）；
- `gstr` 把不带空白的字符串读进缓冲区（需自行保证长度）；
- 输出带 2 MiB 缓冲，正常从 `main` 返回时由全局对象自动 flush；如果程序里调用了 `_exit` / `exit` 或需要立刻看到输出，请手动 `flush()`；
- 与 `printf` / `cout` 混用时，注意内置输出走的是自己的缓冲区，混用前先 `flush()`。

### 3.2 常量与预处理 `CONSTS` / `Pre`

```cpp
P = 998244353, Y = 3, I = 332748118, B = (P + 1) >> 1,
N = 600005, _I_ = 86583718, LG = 19, NL = 1 << LG;

int Pre::Q(int a, int b);      // a^b mod P
int Pre::Inv(int x);           // x^{-1} mod P
void Pre::initYG();            // 必须最先调用

Pre::gp[i], Pre::igp[i];       // 2^19 次单位根及其逆（i ≤ 2^19）
Pre::ny[i];                    // i^{-1}，i ≥ 1
Pre::fac[i], Pre::ifac[i];     // i!，1/i!
```

### 3.3 模运算小工具 `Pint`

```cpp
addt(a, b);  // a += b 并取模，返回 a
delt(a, b);  // a -= b 并取模，返回 a
add(a, b);   // a + b
del(a, b);   // a - b
tp(x);       // 把可能为负的数变成 [0, P)
Vmax(a, b);  // 最大值（之所以有它，是历史原因）
```

均为模板，可传 `int` / `long long`。

### 3.4 二次剩余 `Quad`

```cpp
int Quad::work(int n);   // n 的模 P 平方根
```

- `n = 0` 返回 `0`；
- `n` 是二次剩余时返回较小的那个根；
- `n` 非二次剩余返回 `P + 1`（哨兵值，调用方需要自己判断）。

### 3.5 多项式类 `POLY::poly`

内部只有一个 `vector<int> a`，长度由 vector 决定，没有独立长度字段。

| 接口 | 说明 |
| --- | --- |
| `poly()` | 空多项式 |
| `poly(n)` | 长度为 n 的全 0 多项式 |
| `size()` / `empty()` | 长度 / 是否为空 |
| `f[i]` | 第 i 项系数（有 `int&` 与 `const int&` 两个版本，无越界检查） |
| `resize(n)` | 改长度（补 0 / 截断） |
| `pb(x)` | 尾插一个系数 |
| `rev()` | 整体翻转 |
| `rev(l, r)` | 翻转下标区间 `[l, r)` |
| `tp()` | 全部取到 `[0, P)` |
| `f += g` / `f -= g` / `f + g` / `f - g` | 加减，长度取两者较大值 |
| `f *= k` / `f * k` | 数乘（`k` 可为负） |
| `f *= g` / `f * g` | 卷积 |
| `f & g` | 中段积：`(f & g)[i] = (f·g)[m-1+i]`，其中 `n = \|f\| ≥ m = \|g\|`（多点求值内部使用） |
| `~f` | 形式幂级数求逆，长度不变，要求 `f[0] ≠ 0` |
| `f / g` / `f % g` | 多项式除法 / 取模，长度同 `f`，要求 `deg g ≤ deg f` 且 `\|f\| < 2\|g\|`（内部用反转技巧） |
| `f /= k` / `f / k` | 数除（乘 `k^{-1}`） |
| `f <<= k` / `f << k` / `f >>= k` / `f >> k` | 乘 / 除 `x^k` |
| `f ^= k` / `f ^ k` | `f^k mod x^n`，长度不变，要求 `f[0] ≠ 0` |
| `f.NTT(bool)` | 原地 NTT，`true` 正向、`false` 逆向；长度必须为 2 的幂且 ≤ 2^19 |

```cpp
using ffp = poly;   // 下降幂多项式与普通多项式共用同一个类
```

位运算卷积（`Fwt*`、`And/Or/Xor/Subset`）自 1.17.2 起归 3.6 节的 `sps` 类管。

关于 `NTT`：正向变换后系数是**位反转顺序**（内部约定，为了让正逆变换配成一对、省掉位反转）。点乘后接 `NTT(false)` 即可，不需要自己处理顺序；但不要直接把正向变换的结果当普通顺序读。


### 3.6 集合幂级数 `POLY::sps`

并 / 交 / 异或卷积与子集卷积的容器：系数数组长度为 `2^k`，下标看成 `k` 位二进制集合。

```cpp
sps();                        // 空
sps(n);                       // 长度 n 的全 0（n 必须是 2 的幂）
explicit sps(const poly &f);  // 由多项式（系数数组）构造，要显式写 sps(f)

int size(); bool empty();
f[i];                         // 第 i 项系数
void resize(int n); void reserve(int n);
poly to_poly() const;         // 转回多项式（同一组系数）
```

| 运算 | 含义 |
| --- | --- |
| `f & g` | 与卷积（别名 `f.And(g)`） |
| `f \| g` | 或卷积（别名 `f.Or(g)`） |
| `f ^ g` | 异或卷积（别名 `f.Xor(g)`） |
| `f * g` | **子集卷积**（别名 `f.Subset(g)`）：`h_S = Σ_{T ⊆ S} f_T · g_{S \ T}` |
| `f * k` / `f *= k` | 数乘，`k` 可为负 |
| `f.Fwtand(__)` / `Fwtor` / `Fwtxor` / `Fmtor` | 原地快速沃尔什变换，非 0 正向、`0` 逆向；结果保证落在 `[0, P)` |

```cpp
int n = 1 << 17;
sps f(n), g(n);
for (int i = 0; i < n; ++i)
    f[i] = ..., g[i] = ...;
sps w = f & g;        // 与卷积
sps o = f | g;        // 或卷积
sps x = f ^ g;        // 异或卷积
sps s = f * g;        // 子集卷积
sps t = s * 2;        // 数乘
poly back = s.to_poly();
```

- 两边长度必须相同且是 2 的幂；`&` / `|` / `^` 复杂度 `O(n log n)`，`*`（子集卷积）`O(n log² n)`、额外空间 `O(n log n)`。
- 构造用的 `sps(poly)` 是 `explicit` 的：`poly` 不会悄悄变成 `sps`，`poly | poly` 之类的写法会直接编译报错，需要自己写 `sps(f) | sps(g)`。
- 全局的 `AND/OR/XOR/Subset(poly, poly)` 仍然可用（内部就是上面的运算符），老代码不用改。
- **1.17.1 及以前这些 `Fwt*` 方法是挂在 `poly` 上的**；1.17.2 起统一搬进 `sps`。

### 3.7 多项式输入输出 `POLY::PIO`

```cpp
poly pin(int n);                  // 读 n 个整数
void ppri(const poly &f, int n);  // 输出前 n 项（不足补 0），空格分隔，末尾换行
```

### 3.8 多项式函数 `UCPF`（已 `using` 到全局，可直接裸调）

| 接口 | 语义 | 前置条件 | 返回长度 |
| --- | --- | --- | --- |
| `Dx(f)` | 求导 | — | `n - 1` |
| `Integ(f)` | 积分（常数项 0） | — | `n + 1` |
| `Ln(f)` | `ln(f / f[0]) mod x^n`，常数项为 0 | `f[0] ≠ 0` | `n` |
| `Exp(f)` | `e^f mod x^n` | `f[0] = 0` | `n` |
| `Sqrt(f)` | `sqrt(f) mod x^n` | `f[0]` 是二次剩余 | `n` |
| `Sin(f)` / `Cos(f)` / `Tan(f)` | 三角函数 | `f[0] = 0` | `n` |
| `ArcSin(f)` / `ArcCos(f)` / `ArcTan(f)` | 反三角函数 | `f[0] = 0` | `n` |
| `AND(f, g)` / `OR(f, g)` / `XOR(f, g)` | 位运算卷积（= `sps` 的 `And/Or/Xor`） | `\|f\| = \|g\|` 且为 2 的幂 | 同上 |
| `Subset(f, g)` | 子集卷积 `h_S = Σ_{T ⊆ S} f_T · g_{S \ T}`（= `sps::Subset`） | `\|f\| = \|g\|` 且为 2 的幂 | 同上 |
| `Pow_For_Luogu(f, k1, k2, kk)` | `f^k`，`k` 超大时使用 | `f ≠ 0` | `n` |
| `Trans(f, c)` | 平移 `f(x + c)` | — | `n` |

关于 `ArcCos`：模板里定义为 `-ArcSin(f)`（常数项 0）。真正的 `arccos f` 会比它多一个 `π/2`，请按需自行补常数。

关于 `Pow_For_Luogu(f, k1, k2, kk)`：`k` 大到无法直接存储时使用。

- `kk`：指数本身的数值（用来判断 `f` 的最低非零项次数 × 指数是否已经超过长度，从而直接返回 0）；
- `k2`：参与级数运算的指数（即 `k mod (P - 1)`，`Exp/Ln` 里当作普通整数用）；
- `k1`：常数项幂用的指数（一般取 `k mod (P - 1)`，或 `k1 = k` 当 `k` 不大时）。

若 `f` 有前导零（`f = x^p · g`），结果会自动乘上 `x^{p·k2}`。

### 3.9 多点求值 / 快速插值

```cpp
poly PMPE(const poly &f, const vector<int> &x);          // 在 x[0..m-1] 处求值
poly PFI(const vector<int> &x, const vector<int> &y);    // 过 m 个点的插值多项式
```

- `PMPE` 返回长度 `m` 的值表；要求 `f.size() ≤ m`（即 `deg f < m`）。
- `PFI` 返回次数 `< m` 的插值多项式，长度 `m`；要求 `x` 两两不同。
- 两个函数都是 0 下标（1.17.0 用的是 1 下标裸指针）。

需要反复调用时可以复用对象：

```cpp
Multipoint_Evel::Multipoint mp;      // 也可以按需更换求值点
poly r1 = mp.eval(f1, xs);
poly r2 = mp.eval(f2, xs);

Fast_Interpolation::Interpolation it;
poly p = it.interpolate(xs, ys);     // 同一对象可反复调用
```

低层接口 `Multipoint::build / query`、`Interpolation::build / solve` 也对外开放。

### 3.10 Chirp-Z（等比点求值）

```cpp
poly Chirp_Z(poly f, int c, int m);   // 返回 f(c^0), f(c^1), ..., f(c^{m-1})
int Ci2(int i);                       // C(i, 2) mod (P - 1)，历史遗留工具函数
```

- 返回长度 `m`，要求 `n + m ≤ 2^19`。
- 注意：本模板算的是 `f(c^i)`，没有额外的 `m` 次单位根因子（即不是 `f(c·ω^i)`）。想要 `f(c·ω^i)` 请自己取 `c' = c·ω`。

### 3.11 斯特林数

```cpp
poly Stirling::row(int n);        // 第二类 S(n, 0..n)，长度 n+1
poly Stirling::column(int k, int n);  // 第二类 S(0..n, k)，长度 n+1
poly Staling::row(int n);         // 第一类（无符号）c(n, 0..n)
poly Staling::column(int k, int n);   // 第一类（无符号）c(0..n, k)
```

### 3.12 下降幂 / 普通幂转换

```cpp
poly FFP::FFPTP(const poly &f);   // 普通幂系数 → 下降幂系数
poly FFP::PTFFP(const poly &c);   // 下降幂系数 → 普通幂系数
```

- `FFPTP(f)` 返回 `c`，使得 `f(x) = Σ_i c_i · x(x-1)...(x-i+1)`；返回长度 `n + 1`。
- `PTFFP(c)` 是它的逆运算，返回普通幂系数；返回长度 `n`。
- **注意**：1.17.0 的文件头菜单把这两个名字写反了（实现本身是对的）。以本节为准：`FFPTP` 是「普通 → 下降幂」，`PTFFP` 是「下降幂 → 普通」。

### 3.13 常系数线性递推

齐次式 `a_n = Σ_{i=1..k} c_i · a_{n-i}`：

```cpp
int Fiduccia(int *c, int *a, int k, int n);            // Kitamasa
int Bostan_Mori_CCHLR(int *c, int *a, int k, int n);   // Bostan-Mori
```

- `c[1..k]`：**自然顺序**，`c[i]` 是 `a_{n-i}` 的系数（`c[0]` 不使用）；
- `a[0..k-1]`：初值；
- 返回 `a_n`（`n` 从 0 开始；`n = 0` 返回 `a_0`）。

有理函数系数：

```cpp
int Bostan_Mori(poly f, poly g, int n);   // 返回 [x^n] f/g，要求 deg f < deg g
```

非齐次式 `a_n = Σ_{i=1..k} c_i · a_{n-i} + p(n)`：

```cpp
int Bostam_Mori_CCNHLR(int *c, int *a, int k, int n, const poly &p);
```

- `p` 是多项式 `p(n)` 的系数数组（长度 = 次数 + 1）；
- 参数约定与 `Fiduccia` 相同；
- 复杂度较高，适合 `k`、`deg p` 都不大的场合。

---

## 4. 常用写法

### 4.1 卷积 / 加减 / 数乘

```cpp
poly f = pin(n + 1), g = pin(m + 1);
ppri(f * g, n + m + 1);          // 卷积
f += g;                          // 就地加
f *= 2;                          // 数乘
f <<= 5;                         // 乘 x^5
```

### 4.2 求逆 / 求导 / 积分 / 对数 / 指数

```cpp
poly f = pin(n);
f[0] = ...;                      // 必须非 0
poly inv = ~f;                   // f^{-1} mod x^n
poly l   = Ln(f);                // ln(f / f[0])，常数项 0
poly g   = Exp(poly_minus_1);    // e^g，要求 g[0] == 0
poly d   = Dx(f);                // 求导
poly it  = Integ(d);             // 积分（常数项 0）
```

### 4.3 多项式除法 / 取模

```cpp
poly q = f / g;    // 商
poly r = f % g;    // 余式
```

### 4.4 多点求值 / 快速插值

```cpp
std::vector<int> x = {...}, y = {...};
poly vals = PMPE(f, x);           // f 在各点的值（需要 f.size() <= x.size()）
poly poly_ = PFI(x, y);           // 过这些点的多项式
```

### 4.5 位运算卷积

```cpp
poly f, g;                                    // 长度是 2 的幂
poly h1 = AND(f, g), h2 = OR(f, g), h3 = XOR(f, g), h4 = Subset(f, g);

sps a(f), b(g);                               // 更推荐：直接用集合幂级数
sps w = a & b, o = a | b, x = a ^ b, s = a * b;
```

### 4.6 常系数线性递推第 n 项

```cpp
int c[10] = {0, 1, 1, 0};   // a_n = a_{n-1} + a_{n-2}
int a[10] = {0, 1, 0, 0};   // a_0 = 0, a_1 = 1
printf("%d\n", Fiduccia(c, a, 2, 1000000000));   // 斐波那契第 1e9 项
```

### 4.7 预处理阶乘 / 逆元

```cpp
Pre::fac[i];      // i!
Pre::ifac[i];     // 1 / i!
Pre::ny[i];       // 1 / i
C(n, m) == 1ll * Pre::fac[n] * Pre::ifac[m] % P * Pre::ifac[n - m] % P;
```

---

## 5. 相对 1.17.0 的变化

### 5.1 接口变化（迁移对照表）

| 1.17.0（或更早） | 1.17.2 | 说明 |
| --- | --- | --- |
| `f.n` | `f.size()` | 长度交给 `vector` 自身管理，去掉易失步的 `n` 字段 |
| `f.a`（public） | 私有 | 改由 `size() / operator[] / pb / resize / rev` 访问 |
| `f.shrink()` | 删除 | 不再需要手工同步长度 |
| `f.tp(f.n)` | `f.tp()` | |
| `f.rev(f.a.begin(), f.a.end())` | `f.rev()` | |
| `f.rev(it1, it2)` | `f.rev(l, r)` | 改为下标区间 |
| `pin(f, n)` | `f = pin(n)` | 返回多项式 |
| `ppri(f, n)` | `ppri(f, n)` | 形参改 `const poly&`，不再拷贝 |
| `PMPE(int *val, poly f, int mm, int *ans)` | `poly PMPE(const poly &f, const vector<int> &x)` | 0 下标、返回结果 |
| `PFI(int *x, int *y, int n)` | `poly PFI(const vector<int> &x, const vector<int> &y)` | 同上 |
| `Quad::work(n, p)` | `Quad::work(n)` | `p` 参数原本就被忽略 |
| `#ifdef TRANS / PMPE / FPI / STIRLING / FFP / PTFFP / CCNHLR` | 全部删除 | 所有功能永远编译进来，不需要任何开关宏（第 1.1 节）|
| `Dx(poly&)` / `Integ(poly&)` | `Dx(const poly&)` / `Integ(const poly&)` | |
| `__Binary_reverse` 全局数组 | 删除 | NTT 改为免位反转 |
| `Multipoint_Evel::t,an`、`Fast_Interpolation::q,ret,a`、`Stirling::t`、`FFP::t,ans`、`CCHLR::gamma`、`CCNHLR::gamma,__x,__y`、`Quad::t,n,p,ii` 等全局变量 | 收进 `struct` 成员或函数内局部变量 | 不再污染全局命名空间，也顺带修掉了「同一进程内第二次调用出错」 |
| `class ffp`（与 `poly` 重复实现） | `using ffp = poly;` | 两个类完全同构，合并 |
| `f.Fwtand/Fwtor/Fwtxor/Fmtor`（`poly` 成员） | `sps` 的成员 | 位运算那一路归集合幂级数类，见 3.6 |
| `UCPF::AND/OR/XOR/Subset` 各自手写 FWT 流程 | `sps` 的 `&` / `\|` / `^` / `*` | 只有一处实现，便于以后维护 |
| `Bostam_Mori_CCNHLR` 的 `__calc`（每步都做一次多项式求逆的快速幂） | 合并分母后的单次 `Bostan_Mori` | 见 5.2，常数小很多 |
| `Pint::inv` 未使用数组、`Fmtor` 慢实现等 | 删除 / 转正 | |

### 5.2 常数优化

实测（`g++ -std=c++17 -O2`，Intel Core Ultra 7 255H，取 3 次运行的最小 CPU 时间，两者输出完全相同）：

| 操作 | 规模 | 1.17.0 | 1.17.2 | 加速 |
| --- | --- | --- | --- | --- |
| `f * g` | 200000 × 200000 | 31.3 ms | 16.5 ms | 1.90× |
| `~f` | 200000 | 57.2 ms | 30.9 ms | 1.85× |
| `Ln` | 200000 | 88.9 ms | 49.5 ms | 1.80× |
| `Exp` | 100000 | 136.9 ms | 73.4 ms | 1.87× |
| `Sqrt` | 100000 | 177.6 ms | 94.9 ms | 1.87× |
| `Sin` | 100000 | 274.1 ms | 149.7 ms | 1.83× |
| `f % g` | 200000 | 30.9 ms | 15.4 ms | 2.01× |
| `Trans` | 200000 | 75.0 ms | 17.2 ms | 4.36× |
| `Stirling::row` | 200000 | 57.3 ms | 22.2 ms | 2.58× |
| `PFI` | 50000 点 | 234.3 ms | 131.1 ms | 1.79× |
| `PMPE` | 50000 点 | 127.3 ms | 76.1 ms | 1.67× |
| `Subset` | 2^17 | 109.9 ms | 34.6 ms | 3.18× |
| `Chirp_Z` | 100000 + 100000 | 77.6 ms | 17.1 ms | 4.54× |
| **合计** | | **1478 ms** | **729 ms** | **2.03×** |

非齐次线性递推 `Bostam_Mori_CCNHLR`（P5808，`n = 10^9`，两边都取 `k = m`）：

| 规模 | 1.17.0 | 1.17.1 | 1.17.2 | 相对 1.17.0 |
| --- | --- | --- | --- | --- |
| k = m = 30000 | 3221 ms | 2076 ms | **552 ms** | **5.8×** |

做法：把 `A(X) = Σ a_n X^n` 写成 `Num(X) / ((1-X)^{m+1} D(X))`（`D` 是 k 阶特征多项式），
再用**一次** `Bostan_Mori` 求 `[X^n]`；原实现是 Kitamasa 式的快速幂，每一步都要做一趟多项式除法（内含一次牛顿求逆），常数大得多。
`Num` 由 `C = (Σ_{i<k} a_i X^i)D mod X^k`、`P` 的下降幂（牛顿）系数与两次卷积拼出，全程 `O((k+m) log(k+m) log n)`。

主要手段：

1. **免位反转 NTT**：正向用 DIF、逆向用 DIT，省掉一趟 `O(n)` 交换与一趟 `tp`。
2. **压缩旋转因子表**：`gp` / `igp` 按“块长从大到小、每段连续存 `len/2` 个 `ω_len^j`”布局（总长恰好 `NL-1`）。这样一次长度为 `n` 的变换只顺序访问表中连续的一段（原先按 `NL/len` 步长跳着取，表越大越吃亏），换表后整体又快了一档。
3. **`resize` 不再 `shrink_to_fit()`**：原实现每次 `resize` 都会强制重新分配并整体拷贝，是最大的隐性开销。
4. **小规模卷积走朴素乘法**：`|f|·|g| ≤ MUL_LIM`（默认为 4096，源码里可调）时直接 `O(nm)` 暴力，省掉两次分配 + 三次 NTT。
5. **FWT 去取模**：用 `add / del` 代替 `1ll * a * b % P`；异或 FWT 的 `1/2` 缩放合并为最后一趟；`Fmtor` 直接复用 `Fwtor`。
6. **快速幂改线性递推**：`Trans`、`Chirp_Z` 中逐项 `Q(c, i)` 改为递推（`Chirp_Z` 省得最多）。
7. **`Subset` 内层提指针、只取一次模**。
8. **预处理补上 `fac / ifac`**：`Stirling`、`Staling`、`FFP` 直接用表，不再逐个求逆。
9. **对象瘦身**：`poly` 只有一个 `vector<int>`（24 字节），多点求值等大数组的内存与缓存压力下降约 25%；段树等 scratch 数组按 `4m` 动态分配，不再是固定 `N<<1` / `N<<2`。
10. **非齐次线性递推换算法**：`Bostam_Mori_CCNHLR` 从「每步一次多项式求逆的快速幂」改成「合并分母 + 单次 Bostan-Mori」，P5808 满规模下快 5.8 倍（见上表）。
11. **按需 `reserve`**：`pin(n)` 一次开好、`Dx`/`Integ` 预留长度、段树叶子预留 2 个位置（省掉一半小分配）、`~f` 复用一块缓冲区并把容量一次性预留到 `2n`、卷积里被变换的那一份直接按 NTT 长度构造。

### 5.3 修复的 bug

1. **`PMPE` 返回值错误**：1.17.0 会返回 `x_i · f(x_i) + f(0)` 而不是 `f(x_i)`（快速插值走的是另一套反转变换，所以没暴露）。1.17.1 已修正，可用霍纳法则逐点核对。
2. **scratch 数组只增不减**：`Multipoint_Evel::t`、`Fast_Interpolation::q/ret`、`FFP::t/ans` 等用 `pb` 追加且从不重置，**同一进程内第二次求值会读到残留数据甚至越界崩溃**（1.17.0 里 `PMPE` 调用两次必崩，`Stirling::column` / `FFP::PTFFP` 同理）。1.17.1 每次调用重置。
3. **`Integ` 只更新了系数数组而没更新长度字段**，返回值长度错乱，并连带把 `ArcSin` / `ArcTan` 的长度带错。
4. **`Fiduccia` 在 `n = 0` 时越界读**（`__calc` 返回长度不足 `k`）。1.17.1 起在 `__calc` 末尾补齐长度；1.17.2 重写后天然支持 `n = 0`。
5. **一个零宽字符（U+FFFC）** 藏在 1.17.0 的 `Bostam_Mori_CCNHLR` 里，导致开启 `CCNHLR` 时无法编译。
6. **`operator^=` / `Pow_For_Luogu`** 在多项式全 0 时读取未初始化变量（现已返回 0 多项式，`^0` 返回 1）。
7. **`f >>= k` 在 `k ≥ |f|` 时的越界 `erase`**。
8. `~f` 对空多项式、`f * g` 对空多项式等边界不再越界。

### 5.4 代码风格

与 1.17.0 保持一致：`clang-format` 配置为

```
BasedOnStyle: LLVM, IndentWidth: 4, ColumnLimit: 80,
NamespaceIndentation: All, AccessModifierOffset: -4, FixNamespaceComments: false
```

命名沿用 `__` 前缀的局部量 / 形参、大写下划线的常量；全文件无注释，函数体内无空行。

---

## 6. 注意事项

1. **别忘了 `Pre::initYG()`**，否则所有 NTT 都是错的。
2. `Ln` / `Sqrt` / `Exp` / 三角函数的常数项要求见 3.8；`Ln` 会自动把常数项归一化掉（返回 `ln(f / f[0])`）。
3. **`Exp` 及其下游函数的规模上限是 `2^17`**，比其它函数低，n 大时特别留意。
4. 自己调用 `NTT` 时长度必须是 2 的幂且 ≤ `2^19`，并且要成对使用（正向的位反转顺序不要直接解读）。
5. OJ 上注意编译参数：`-O2`、`-std=c++17`；没有开关宏，所有功能永远可用（2.x 直接 include 1.x，也不需要额外参数）。
6. 本模板面向 `998244353`，需要换模数时请注意 `P` 的 2 进幂次与根表初始化（`CONSTS::LG`、`Pre::initYG`）。
7. 递归较深的接口（多点求值、插值、斯特林）内部 scratch 都在堆上（`vector`），栈压力很小。

---

## 7. 矩阵扩展（`2_1_1.cpp`）

`2_1_1.cpp` 的第一行就是 `#include "1_17_2.cpp"`：多项式部分完全复用 1.x，文件里只多一个 `POLY::matrix`。

### 7.1 为什么能这么写

`1_17_2.cpp` 末尾的 `work()` 与 `main()` 包在 `#if __INCLUDE_LEVEL__ == 0 ... #endif` 里：**直接编译** `1_17_2.cpp` 时包含层级是 0，会带上这套提交骨架；被别的文件 `#include` 时层级 ≥ 1，自动跳过，好让包含者写自己的 `work()` / `main()`。于是：

| 编译谁 | 得到什么 |
| --- | --- |
| `g++ 1_17_2.cpp` | 1.x（做题时把 `work()` 换成题目逻辑） |
| `g++ 2_1_1.cpp` | 2.x：矩阵 + 全部多项式功能 |

以后新开 2.x 文件，**第一行写 `#include "1_17_2.cpp"`，其余照旧**（`work()` / `main()` 自己写，不会和 1.x 的那份冲突）。被测文件（例如自测程序）同样可以 `#include "2_1_1.cpp"` 后自己写 `main()`。

### 7.2 接口

```cpp
matrix f(n, m);           // n 行 m 列，元素全 0
f.resize(n, m);           // 缩小时丢尾巴，扩大时新格子补 0
f.row(); f.col();         // 行数 / 列数
f[i][j]; f(i, j);         // 第 i 行第 j 列，1-based

f += g; f -= g; f *= g;   // 复合运算
f + g;  f - g;  f * g;    // 返回新矩阵
-f;                       // 逐元素取相反数（0 仍是 0）
f.gauss();                // 原地化为行最简形 RREF（解方程组用）
f.det();                  // 行列式 O(n^3)，不改动 f
~f;                       // 逆矩阵，要求方阵且可逆
```

矩阵输入输出（与多项式的 `pin` / `ppri` 重载共存，都在 `POLY::PIO`）：

```cpp
pin(f, n, m);    // 读 n*m 个整数，内部自动 resize
ppri(f, n, m);   // 按行输出，元素间一个空格、行尾换行
```

### 7.3 例子

```cpp
int n;

inline void work()
{
    read(n);
    matrix f(n, n);
    pin(f, n, n);
    print(f.det());
    pc('\n');
    ppri(~f, n, n);
}
```

验算方法：`f * ~f` 是单位矩阵；可逆矩阵 `f.gauss()` 之后是单位矩阵。`A(n,k) *= B(k,m)` 要求 `k` 对上、`n`/`m` 可任意（`det` / `~` 只对方阵有意义）。

### 7.4 相对 2.1.0 的变化

- **`-A` 修好了**：2.1.0 会输出 `P` 而不是 `0`（未取模）。
- **`gauss()` 重写**：2.1.0 靠"把下面的行加到主元行"来躲开零主元，列数多于行数或需要换行时会错位；现在一律行交换 + 记录换行符号。
- **`~A` 是新的**：增广 `[A | I]` 后一次高斯-约当消元得到 `A^{-1}`。
- **`det()` 记录换行符号**，`k` 循环只做到 `n-1`（最后一次乘法在末尾统一做）。
- **`*=` 改 i-k-j 顺序 + 行指针 + 零元跳过**，比 2.1.0 的 i-j-k 三重循环快约 2~3 倍。
- **去掉 `CP()`（特征多项式）**：2.1.0 的做法是高斯消元后算 `∏(x - a_ii)`，但高斯消元不是相似变换，结果并不是特征多项式。需要的话用 Hessenberg 化重写（说一声即可）。

---

## 8. 作者

```
_____________________________________________
                                          
  AUTHER   : Revitalize                   
  Nation   : China                        
  Province : Shaanxi                      
  School   : Middle School Attached to NPU  
  LOCATION : Zhejiang Shaoxing            
  TIME     : 2025/7/16                    
____________________________________________
```

![](https://cdn.luogu.com.cn/upload/image_hosting/95rarpii.png)
