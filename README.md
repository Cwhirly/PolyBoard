# PolyBoard 1.19.0

PolyBoard 是面向 OI 的单文件 C++ 多项式与矩阵模板，默认模数 `998244353`。主要功能包括 NTT / MTT 卷积、形式幂级数运算、多点求值与快速插值、集合幂级数与位运算卷积、斯特林数、下降幂转换、常系数线性递推，以及矩阵基本运算、高斯消元、行列式与求逆。

- 主文件：`poly.cpp`（多项式）、`matrix.cpp`（矩阵，第一行 `#include "poly.cpp"`）。
- 版本管理：文件名固定，版本由 git tag / release 区分，历史文件见 `GIT.md`。
- 更新记录：`CHANGELOG.md`。
- 编译环境：GCC，`bits/stdc++.h`，C++14 及以上（推荐 `-std=c++17`）。
- 代码风格：文件内无注释，函数体内无空行；使用说明集中在本文件。

---

## 1. 快速开始

### 1.1 编译

```bash
g++ poly.cpp   -o main -std=c++17 -O2
g++ matrix.cpp -o main -std=c++17 -O2
```

模板不提供模块级开关宏，所有功能都会参与编译。可配置项为 `POLYBOARD_P`、`POLYBOARD_N`、`POLYBOARD_ATN`、`POLYBOARD_P_FIXED` 四个宏，见第 2 节。完整编译约 2 秒。

### 1.2 初始化

固定模数：

```cpp
signed main() {
    Pre::initYG();
    work();
    return 0;
}
```

运行时模数：

```cpp
void work() {
    int p;
    read(p);
    Pre::init(p);
    ...
}
```

`Pre::initYG()` / `Pre::init(p)` 必须在任何多项式或矩阵接口之前调用。前者初始化 NTT 单位根表（仅 `ATN=1`）、`ny`、`fac`、`ifac`；后者额外设置运行时的 `P`、Barrett 状态、`LG/NL` 和阶乘表。

### 1.3 最小示例

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

文件末尾的 `work()` / `main()` 是一个求卷积的示例，正式使用时替换为题目逻辑。

---

## 2. 配置

### 2.1 配置宏

四个宏都可以在 `#include "poly.cpp"` 之前定义，也可以直接修改文件开头的默认值。

```cpp
#ifndef POLYBOARD_P
#define POLYBOARD_P 998244353
#endif
#ifndef POLYBOARD_N
#define POLYBOARD_N 600005
#endif
#ifndef POLYBOARD_ATN
#define POLYBOARD_ATN 1
#endif
#ifndef POLYBOARD_P_FIXED
#define POLYBOARD_P_FIXED 1
#endif
```

| 宏 | 默认值 | 含义 |
| --- | --- | --- |
| `POLYBOARD_P` | `998244353` | 模数；`P_FIXED=1` 时为编译期常量，`P_FIXED=0` 时为运行时初值 |
| `POLYBOARD_N` | `600005` | 多项式长度、`ny/fac/ifac` 等表的容量上限 |
| `POLYBOARD_ATN` | `1` | `1` 使用 NTT，`0` 使用 MTT |
| `POLYBOARD_P_FIXED` | `1` | `1` 表示模数编译期确定，`0` 表示模数运行时输入 |

### 2.2 固定模数

```cpp
#define POLYBOARD_P 1000000007
#define POLYBOARD_N 600005
#define POLYBOARD_ATN 0
#define POLYBOARD_P_FIXED 1
#include "poly.cpp"
```

- `ATN=1`：`P` 必须是 NTT 素数。`LG = min(ceil(log2(2N-1)), v2(P-1))`，`NL = 2^LG`。编译期检查 `P ≡ 1 (mod 4)`、`2^LG | P-1` 和原根存在性。
- `ATN=0`：`LG = ceil(log2(2N-1))`，`NL = 2^LG`。`Pre::initYG` 跳过 NTT 旋转因子表，乘法走 5 次 FFT 的 `poly::MTT`。
- `P_FIXED=1`：`P` 为 `constexpr`，`mod(x)` 内联为 `x % P`，没有运行时模数开销。

固定模数下，`Y`、`I = 1/Y`、`B = 1/2`、`_I_`（四次单位根）均由编译期计算。

### 2.3 运行时模数

```cpp
#define POLYBOARD_P_FIXED 0
#define POLYBOARD_ATN 0
#include "poly.cpp"

void work() {
    int p;
    read(p);
    Pre::init(p);
    ...
}
```

- `P_FIXED=0` 时，`P` 是运行时 `int`，必须令 `ATN=0`；`static_assert` 会拒绝 `ATN=1`。
- 必须在 `main` 或 `work` 开头调用 `Pre::init(p)`。`Pre::init` 完成 `P`、`B`、`LG`、`NL`、Barrett 状态和 `ny/fac/ifac` 的初始化。
- 所有 `% P` 已统一替换为 `mod`。固定模数下 `mod` 编译为 `% P`；运行时模数下 `mod` 使用 Barrett。
- 运行时不计算 `Y`、`I`、`_I_`。依赖四次单位根的函数（`Sin`、`Cos` 等）建议只在固定模数 `ATN=1` 下使用。
- `matrix.cpp` 同样支持运行时模数，只需在使用矩阵前调用 `Pre::init(p)`。

### 2.4 可用性

| 模数 | `ATN` | 多项式乘法 | 形式幂级数 | 矩阵 |
| --- | --- | --- | --- | --- |
| `998244353`（默认） | `1` | NTT | 全功能 | 全功能 |
| 其它 NTT 素数 | `1` | NTT | 全功能 | 全功能 |
| 非 NTT 素数（如 `1000000007`） | `0` | MTT | 需要元素可逆；`Sin/Cos` 还要求 `P ≡ 1 (mod 4)` | 全功能 |
| 合数（如 `12`） | `0` | MTT | 保证乘法与加减；除法、求逆、`Ln/Exp` 可能失败 | 基本运算正确；`gauss/inv` 要求主元可逆，否则返回 `false` |

`Pre::Inv` 使用扩展欧几里得算法，只要元素本身可逆，任意模数下都能求出逆元；不可逆时返回 `0`。

### 2.5 `mod` 与 Barrett

`mod(long long x)` 返回 `[0, P)` 内的值：

- `P_FIXED=1`：`x % P`，由编译器按常量模数优化。
- `P_FIXED=0`：Barrett 归约。`bar` 在 `Pre::init(p)` 中初始化。

新增代码中的模运算应写成 `mod(...)`，不要直接写 `% P`。这样固定模数路径保持常量优化，运行时路径自动使用 Barrett。

---

## 3. 容量与限制

### 3.1 常量

默认配置下：

| 名称 | 值 | 说明 |
| --- | --- | --- |
| `P` | `998244353` | 模数 |
| `Y` | `3` | 原根 |
| `I` | `332748118` | `Y^{-1}` |
| `B` | `499122177` | `2^{-1}` |
| `_I_` | `911660635` | 四次单位根 |
| `N` | `600005` | 表与用户数组容量 |
| `LG` | `21` | `log2(NL)` |
| `NL` | `2097152` | NTT / MTT 最大变换长度 |
| `MUL_LIM` | `4096` | 小规模卷积转暴力的阈值 |

`Pre::gp`、`Pre::igp` 的长度为 `NL`；`Pre::ny`、`Pre::fac`、`Pre::ifac` 的长度为 `N`。

### 3.2 接口上限

下表为 `N = 600005`、`ATN = 1` 时的实测上限。超出上限属于未定义行为，模板不会自动截断或报错。

| 接口 | 上限 | 来源 |
| --- | --- | --- |
| `f * g`、`f & g` | `f.size() + g.size() - 1 ≤ NL` | NTT 长度 |
| `~f`、`f / g`、`f % g` | `n ≤ NL / 2`（约 `2^20`） | Newton 迭代内部 NTT 长度 |
| `Chirp_Z` | `n + m ≤ NL / 2 + 1` | 内部卷积长度 |
| `Dx`、`Integ`、`Ln`、`Trans`、`Stirling`、`Staling`、`FFP` | `n ≤ N` | `ny/fac/ifac` 表长 |
| `Exp`、`Sqrt`、`Sin`、`Cos`、`Tan`、`ArcSin`、`ArcCos`、`ArcTan` | `n ≤ NL / 4`（约 `2^19`） | 算法内部需要 `4 · dep ≤ NL` |
| `PMPE`、`PFI` | 点数 `m ≤ (NL - 2) / 3` | 内部卷积与段树内存 |
| `AND`、`OR`、`XOR`、`Subset` | 长度为 2 的幂 | FWT |

内部使用 `Pre::ny`、`Pre::fac`、`Pre::ifac` 的接口还要求 `n ≤ N`。需要更大规模时修改 `N`，并确认 `2^LG` 不超过 `P-1` 的 2 的幂次。

### 3.3 内存

- `Pre::gp`、`Pre::igp`：各 `NL` 个 `int`（默认各 8 MiB）。`ATN=0` 时不构建，但固定模数下仍按 `NL` 预留；运行时模数下按 `NLMAX` 预留。
- `Pre::ny`、`Pre::fac`、`Pre::ifac`：各 `N` 个 `int`。
- `poly` 只含一个 `vector<int>`，对象本身约 24 字节。
- 多点求值与快速插值的段树按 `4m` 动态分配。

减小 `N` 会同步减小上述静态数组。

---

## 4. 基础类型

### 4.1 `POLY::poly`

内部只有一个 `vector<int> a`，长度由 `vector` 维护。

| 接口 | 说明 |
| --- | --- |
| `poly()` | 空多项式 |
| `poly(n)` | 长度为 `n` 的全零多项式 |
| `size()`、`empty()` | 长度、是否为空 |
| `f[i]` | 第 `i` 项系数，提供 `int&` 与 `const int&` 两个版本，无越界检查 |
| `resize(n)` | 调整长度，补零或截断 |
| `pb(x)` | 尾部插入系数 |
| `rev()` | 整体翻转 |
| `rev(l, r)` | 翻转下标区间 `[l, r)` |
| `tp()` | 所有系数归一到 `[0, P)` |
| `f += g`、`f -= g`、`f + g`、`f - g` | 加减，长度取较大者 |
| `f *= k`、`f * k` | 数乘，`k` 可为负 |
| `f *= g`、`f * g` | 卷积 |
| `f & g` | 乘积高段：`(f & g)[i] = (f·g)[m-1+i]`，其中 `n = f.size() ≥ m = g.size()`，多点求值内部使用 |
| `~f` | 形式幂级数求逆，长度不变，要求 `f[0] ≠ 0` |
| `f / g`、`f % g` | 多项式除法、取模，要求 `deg g ≤ deg f` 且 `f.size() < 2 * g.size()` |
| `f /= k`、`f / k` | 数除，乘以 `k^{-1}` |
| `f <<= k`、`f << k`、`f >>= k`、`f >> k` | 乘或除 `x^k` |
| `f ^= k`、`f ^ k` | `f^k mod x^n`，长度不变，要求 `f[0] ≠ 0` |
| `f.NTT(bool)` | 原地 NTT，`true` 正向，`false` 逆向；长度必须是 2 的幂且不超过 `NL` |

`using ffp = poly;`。下降幂多项式与普通多项式共用 `poly`。

NTT 正向变换的输出为位反转顺序，这是内部约定；点乘后调用 `NTT(false)` 即可，不需要手动位反转。正向结果不能按自然顺序直接读取。

### 4.2 `POLY::sps`

集合幂级数，用于并、交、异或卷积与子集卷积。系数数组长度为 `2^k`，下标视为 `k` 位二进制集合。

```cpp
sps();                        // 空
sps(n);                       // 长度 n 的全零，n 必须是 2 的幂
explicit sps(const poly &f);  // 由系数数组构造

int size(); bool empty();
f[i];
void resize(int n); void reserve(int n);
poly to_poly() const;
```

| 运算 | 含义 |
| --- | --- |
| `f & g` | 与卷积，别名 `f.And(g)` |
| `f.Or(g)` | 或卷积 |
| `f ^ g` | 异或卷积，别名 `f.Xor(g)` |
| `f * g` | 子集卷积，别名 `f.Subset(g)`；`h_S = Σ_{T ⊆ S} f_T · g_{S \ T}` |
| `f * k`、`f *= k` | 数乘，`k` 可为负 |
| `f.Fwtand(x)`、`Fwtor`、`Fwtxor`、`Fmtor` | 原地沃尔什变换；非零为正向，`0` 为逆向；结果归一到 `[0, P)` |

- 两侧长度必须相同且为 2 的幂。
- `&`、`|`、`^` 复杂度 `O(n log n)`；子集卷积 `O(n log² n)`，额外空间 `O(n log n)`。
- `sps(poly)` 为 `explicit`，`poly` 不会隐式转换为 `sps`。
- 全局 `AND/OR/XOR/Subset(poly, poly)` 仍然可用，内部调用 `sps` 的对应运算。

### 4.3 `MATRIX::matrix`

矩阵类位于 `matrix.cpp` 的 `MATRIX` 命名空间，下标从 1 开始。

```cpp
MATRIX::matrix f(n, m);       // n 行 m 列
MATRIX::matrix g;             // 空矩阵
f.row(); f.col();
f(i, j);                      // 第 i 行第 j 列，1 下标
```

| 接口 | 说明 |
| --- | --- |
| `f += g`、`f -= g`、`f + g`、`f - g` | 矩阵加减，要求尺寸相同 |
| `f *= k`、`f * k` | 数乘 |
| `f *= g`、`f * g` | 矩阵乘法 |
| `f ^= k`、`f ^ k` | 矩阵快速幂，要求方阵 |
| `f.gauss()` | 高斯消元，返回是否满行秩 |
| `f.inv(g)` | 求逆并写入 `g`，返回是否可逆；失败时 `g` 为全零矩阵 |
| `~f` | 求逆；不可逆时返回全零矩阵 |
| `f.det()` | 行列式 |

### 4.4 `ffp`

`ffp` 是 `poly` 的别名，用于下降幂多项式。相关转换函数见 5.8 节。

---

## 5. API 参考

### 5.1 输入输出 `IO`

```cpp
template <class I> bool read(I &x);
template <class I> void print(I x);
void pc(char c);
void pstr(const char *s);
void gc(char &c);
void gstr(char *s);
void flush();
```

- `read` 返回是否成功读取；遇到 EOF 返回 `false`。
- `print` 支持负数，不输出多余空格。
- `pc` 输出单个字符；`pstr` 输出到 `'\0'` 为止。
- `gstr` 读入不含空白的字符串，调用方负责缓冲区长度。
- 输出使用 2 MiB 缓冲，正常从 `main` 返回时自动刷新。调用 `exit` / `_exit` 或需要立即输出时，手动 `flush()`。
- 与 `printf` / `cout` 混用时，先 `flush()`。

### 5.2 常量与预处理 `CONSTS` / `Pre`

```cpp
int Pre::Q(int a, int b);   // a^b mod P
int Pre::Inv(int x);        // x^{-1} mod P，扩展欧几里得；不可逆返回 0
void Pre::initYG();         // 固定模数初始化
void Pre::init(int p);      // 运行时模数初始化，仅 P_FIXED=0 可用

Pre::gp[i], Pre::igp[i];    // NTT 单位根及其逆
Pre::ny[i];                 // i^{-1}
Pre::fac[i], Pre::ifac[i];  // i!，1/i!
```

`Pre::gp/igp` 的布局为按块长从大到小、每段连续存放 `len/2` 个单位根，总长度不超过 `NL-1`。

### 5.3 模运算工具 `Pint`

```cpp
Pint::addt(a, b);  // a += b 并取模，返回 a
Pint::delt(a, b);  // a -= b 并取模，返回 a
Pint::add(a, b);   // a + b mod P
Pint::del(a, b);   // a - b mod P
Pint::tp(x);       // 将可能为负的数归一到 [0, P)
Pint::Vmax(a, b);  // 最大值
```

以上均为模板，可传 `int` 或 `long long`。

### 5.4 二次剩余 `Quad`

```cpp
int Quad::work(int n);   // n 在模 P 下的平方根
```

- `n = 0` 返回 `0`。
- `n` 为二次剩余时返回较小的根。
- `n` 非二次剩余时返回 `P + 1`，调用方需自行判断。

### 5.5 多项式输入输出 `POLY::PIO`

```cpp
poly pin(int n);                  // 读入 n 个整数
void ppri(const poly &f, int n);  // 输出前 n 项，不足补零；空格分隔，末尾换行
```

### 5.6 多项式函数 `UCPF`

`UCPF` 中的函数已通过 `using` 暴露到全局，可直接调用。

| 接口 | 语义 | 前置条件 | 返回长度 |
| --- | --- | --- | --- |
| `Dx(f)` | 求导 | 无 | `n - 1` |
| `Integ(f)` | 积分，常数项为 0 | 无 | `n + 1` |
| `Ln(f)` | `ln(f / f[0]) mod x^n`，常数项为 0 | `f[0] ≠ 0` | `n` |
| `Exp(f)` | `e^f mod x^n` | `f[0] = 0` | `n` |
| `Sqrt(f)` | `sqrt(f) mod x^n` | `f[0]` 为二次剩余 | `n` |
| `Sin(f)`、`Cos(f)`、`Tan(f)` | 三角函数 | `f[0] = 0` | `n` |
| `ArcSin(f)`、`ArcCos(f)`、`ArcTan(f)` | 反三角函数 | `f[0] = 0` | `n` |
| `AND(f, g)`、`OR(f, g)`、`XOR(f, g)` | 位运算卷积 | 长度相同且为 2 的幂 | 同输入 |
| `Subset(f, g)` | 子集卷积 | 长度相同且为 2 的幂 | 同输入 |
| `Pow_For_Luogu(f, k1, k2, kk)` | `f^k`，指数超大时使用 | `f ≠ 0` | `n` |
| `Trans(f, c)` | 平移 `f(x + c)` | 无 | `n` |

- `ArcCos(f)` 定义为 `-ArcSin(f)`，常数项为 0；真正的 `arccos` 会多一个 `π/2`，需要时自行补常数。
- `Pow_For_Luogu` 的三个指数参数：
  - `kk`：指数的真实值，用于判断 `f` 的最低非零项次数乘以指数是否已超过长度，从而直接返回 0；
  - `k2`：参与级数运算的指数，通常取 `k mod (P-1)`；
  - `k1`：常数项幂使用的指数，通常取 `k mod (P-1)`，`k` 较小时可直接取 `k`。
- 若 `f = x^p · g`，结果会自动乘上 `x^{p·k2}`。

### 5.7 多点求值与快速插值

```cpp
poly PMPE(const poly &f, const vector<int> &x);
poly PFI(const vector<int> &x, const vector<int> &y);
```

- `PMPE` 返回长度 `m` 的值表，要求 `f.size() ≤ m`。
- `PFI` 返回次数小于 `m` 的插值多项式，长度 `m`，要求 `x` 两两不同。
- 两者均为 0 下标。

需要重复调用时可复用对象：

```cpp
Multipoint_Evel::Multipoint mp;
poly r1 = mp.eval(f1, xs);
poly r2 = mp.eval(f2, xs);

Fast_Interpolation::Interpolation it;
poly p = it.interpolate(xs, ys);
```

低层接口 `Multipoint::build/query`、`Interpolation::build/solve` 同样对外开放。

### 5.8 Chirp-Z

```cpp
poly Chirp_Z(poly f, int c, int m);   // f(c^0), f(c^1), ..., f(c^{m-1})
int Ci2(int i);                       // C(i, 2) mod (P - 1)
```

- 返回长度 `m`，要求 `n + m ≤ NL / 2 + 1`。
- 计算的是 `f(c^i)`，不含额外的单位根因子。需要 `f(c·ω^i)` 时，令 `c' = c·ω` 后调用。

### 5.9 斯特林数与下降幂

```cpp
poly Stirling::row(int n);              // 第二类 S(n, 0..n)
poly Stirling::column(int k, int n);    // 第二类 S(0..n, k)
poly Staling::row(int n);               // 第一类无符号 c(n, 0..n)
poly Staling::column(int k, int n);     // 第一类无符号 c(0..n, k)

poly FFP::FFPTP(const poly &f);         // 普通幂系数 -> 下降幂系数
poly FFP::PTFFP(const poly &c);         // 下降幂系数 -> 普通幂系数
```

- `FFPTP(f)` 返回 `c`，使 `f(x) = Σ c_i · x(x-1)…(x-i+1)`，返回长度 `n + 1`。
- `PTFFP` 为其逆运算，返回普通幂系数，长度 `n`。

### 5.10 常系数线性递推

齐次递推 `a_n = Σ_{i=1..k} c_i · a_{n-i}`：

```cpp
int Fiduccia(int *c, int *a, int k, int n);
int Bostan_Mori_CCHLR(int *c, int *a, int k, int n);
```

- `c[1..k]` 使用自然顺序，`c[i]` 是 `a_{n-i}` 的系数，`c[0]` 不使用。
- `a[0..k-1]` 为初值。
- 返回 `a_n`，`n` 从 0 开始。

有理函数系数：

```cpp
int Bostan_Mori(poly f, poly g, int n);   // [x^n] f/g，要求 deg f < deg g
```

非齐次递推 `a_n = Σ_{i=1..k} c_i · a_{n-i} + p(n)`：

```cpp
int Bostam_Mori_CCNHLR(int *c, int *a, int k, int n, const poly &p);
```

- `p` 是 `p(n)` 的系数数组，长度为次数加一。
- 参数约定与 `Fiduccia` 相同。
- 复杂度适合 `k` 与 `deg p` 均不太大的场景。

### 5.11 矩阵输入输出与接口

```cpp
MATRIX::matrix f = MATRIX::MIO::pin(n, m);
MATRIX::MIO::ppri(f, n, m);
```

`MATRIX::MIO::pin` 按行读入 `n × m` 个整数；`ppri` 输出前 `n` 行、前 `m` 列，不足补零。矩阵接口见 4.3 节。

---

## 6. 内部算法

### 6.1 NTT

正向使用 DIF（自然输入、位反转输出），逆向使用 DIT（位反转输入、自然输出）。正逆变换共用同一张旋转因子表，不进行显式的位反转或 `tp` 扫描。

`ATN=1` 时，`LG = min(ceil(log2(2N-1)), v2(P-1))`，`NL = 2^LG`。变换长度必须是 2 的幂且不超过 `NL`。

### 6.2 MTT

`ATN=0` 时，卷积由 `poly::MTT` 完成：

- 将系数拆成 15 位高低两段；
- 使用 3 次正向 FFT、2 次逆向 FFT（共 5 次）；
- 正向 DIF、逆向 DIT，与 NTT 的结构对应；
- 逆变换后按 `ac + (ad + bc)·2^15 + bd·2^30` 合并，并归约到 `[0, P)`。

FFT 实现位于全局 `namespace FFT` 的 `run` 函数，模板参数为复数类型，旋转因子表为函数内静态表。

### 6.3 Newton 迭代

求逆、`Ln`、`Exp`、`Sqrt` 使用牛顿迭代，迭代过程中通过 NTT 或 MTT 完成卷积。`ATN=0` 时走通用牛顿迭代，`ATN=1` 时走原有的专用 NTT 路径。

### 6.4 FWT 与子集卷积

`sps` 的与、或、异或卷积使用原地沃尔什变换；子集卷积按集合大小分层，复杂度 `O(n log² n)`。变换结果在每层结束后归一到 `[0, P)`。

### 6.5 非齐次线性递推

将生成函数写成

\[
A(X) = \frac{\operatorname{Num}(X)}{(1-X)^{m+1} D(X)}
\]

其中 `D` 是特征多项式。随后用一次 Bostan-Mori 求 `[X^n]`，避免逐步多项式除法。`Num` 由初值、`D` 的截断和 `p(n)` 的下降幂系数通过两次卷积得到。

### 6.6 模运算

固定模数下 `mod(x)` 为 `x % P`；运行时模数下使用 Barrett 归约。NTT 蝶形、点乘、`Get` 及矩阵内层中的乘法均通过 `mod` 完成。

---

## 7. 示例

### 7.1 卷积、加减、数乘

```cpp
poly f = pin(n + 1), g = pin(m + 1);
ppri(f * g, n + m + 1);
f += g;
f *= 2;
f <<= 5;
```

### 7.2 求逆、求导、积分、对数、指数

```cpp
poly f = pin(n);
f[0] = ...;              // 非零
poly inv = ~f;           // f^{-1} mod x^n
poly l = Ln(f);          // ln(f / f[0])
poly g = Exp(poly_1);    // e^g，g[0] = 0
poly d = Dx(f);
poly it = Integ(d);
```

### 7.3 多项式除法与取模

```cpp
poly q = f / g;
poly r = f % g;
```

### 7.4 多点求值与快速插值

```cpp
std::vector<int> x = {...}, y = {...};
poly vals = PMPE(f, x);
poly p = PFI(x, y);
```

### 7.5 位运算卷积

```cpp
poly f, g;  // 长度为 2 的幂
poly h1 = AND(f, g), h2 = OR(f, g), h3 = XOR(f, g), h4 = Subset(f, g);

sps a(f), b(g);
sps w = a & b, o = a | b, x = a ^ b, s = a * b;
```

### 7.6 常系数线性递推

```cpp
int c[10] = {0, 1, 1, 0};   // a_n = a_{n-1} + a_{n-2}
int a[10] = {0, 1, 0, 0};   // a_0 = 0, a_1 = 1
printf("%d\n", Fiduccia(c, a, 2, 1000000000));
```

### 7.7 阶乘与逆元

```cpp
Pre::fac[i];      // i!
Pre::ifac[i];     // 1 / i!
Pre::ny[i];       // 1 / i
C(n, m) == mod(1ll * Pre::fac[n] * Pre::ifac[m] * Pre::ifac[n - m]);
```

### 7.8 矩阵

```cpp
MATRIX::matrix f = MATRIX::MIO::pin(n, n), g;
if (!f.inv(g)) {
    print(-1);
    pc('\n');
} else {
    MATRIX::MIO::ppri(g, n, n);
}
```

---

## 8. 性能

测试环境：`g++ -std=c++17 -O2`，Intel Core Ultra 7 255H；取三次运行的最小 CPU 时间。不同机器与编译器的结果会有差异。

### 8.1 1.17.0 与 1.17.2 对比

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
| `Subset` | `2^17` | 109.9 ms | 34.6 ms | 3.18× |
| `Chirp_Z` | 100000 + 100000 | 77.6 ms | 17.1 ms | 4.54× |

### 8.2 非齐次线性递推

`Bostam_Mori_CCNHLR`，`n = 10^9`，`k = m`：

| 规模 | 1.17.0 | 1.17.1 | 1.17.2 |
| --- | --- | --- | --- |
| `k = m = 30000` | 3221 ms | 2076 ms | 552 ms |

### 8.3 当前版本

`N = 600005`、`P = 998244353`：

| 操作 | 规模 | 时间 |
| --- | --- | --- |
| NTT 卷积 | `10^5` | 约 15 ms |
| MTT 卷积 | `10^5` | 约 28 ms |
| NTT 卷积 | `6 × 10^5` | 约 102 ms |
| MTT 卷积 | `6 × 10^5` | 约 260 ms |

运行时模数 `P = 1000000007`、MTT、`10^5` 规模约 26 ms。

---

## 9. 版本与迁移

### 9.1 1.17.0 到 1.19.0 的接口变化

| 旧接口 | 新接口 | 说明 |
| --- | --- | --- |
| `f.n` | `f.size()` | 长度由 `vector` 管理 |
| `f.a`（public） | 私有 | 通过 `size/operator[]/pb/resize/rev` 访问 |
| `f.shrink()` | 删除 | 不再需要手工同步长度 |
| `f.tp(f.n)` | `f.tp()` | |
| `f.rev(begin, end)` | `f.rev()` 或 `f.rev(l, r)` | 迭代器版本改为下标区间 |
| `pin(f, n)` | `f = pin(n)` | 返回多项式 |
| `PMPE(int *val, poly f, int mm, int *ans)` | `PMPE(const poly&, const vector<int>&)` | 0 下标，返回结果 |
| `PFI(int *x, int *y, int n)` | `PFI(const vector<int>&, const vector<int>&)` | 同上 |
| `Quad::work(n, p)` | `Quad::work(n)` | `p` 参数原本未使用 |
| `Dx(poly&)`、`Integ(poly&)` | `Dx(const poly&)`、`Integ(const poly&)` | |
| `#ifdef TRANS / PMPE / FPI / STIRLING / FFP / PTFFP / CCNHLR` | 删除 | 所有功能始终参与编译 |
| `class ffp` | `using ffp = poly` | 两个类同构，合并 |
| `f.Fwtand/Fwtor/Fwtxor/Fmtor`（`poly` 成员） | `sps` 成员 | 位运算卷积归入集合幂级数 |
| `POLY::MTT` 占位实现 | `poly::MTT` 五次 FFT | 1.19.0 实现 |

1.19.0 还新增了 `POLYBOARD_P_FIXED` 与 `mod`：`P_FIXED=0` 时 `P` 由输入决定，必须 `ATN=0`，并在使用前调用 `Pre::init(p)`；`poly` 与 `matrix` 的模运算统一走 `mod`。

### 9.2 矩阵 2.1.x

- `matrix` 位于 `MATRIX` 命名空间，下标从 1 开始。
- `gauss()` 返回是否满行秩；`inv(g)` 返回是否可逆，失败时 `g` 置全零；`~` 不可逆时返回全零矩阵。
- `det()` 记录换行符号，`k` 循环只到 `n-1`，最后一次乘法在末尾统一完成。
- `*=` 采用 i-k-j 顺序、行指针与零元跳过。
- 消元主元要求可逆，配合扩展欧几里得的 `Pre::Inv`，合数模数下结论仍然可靠。
- `pin` 改为返回值写法：`MATRIX::matrix f = MATRIX::MIO::pin(n, m);`。
- `CP()`（特征多项式）已移除；原实现基于高斯消元，结果并非特征多项式。

### 9.3 已修复问题

- `multipoint_eval_query` 在 `l == r` 时提前返回，消除了死循环与越界。
- `Pint::del` 的负数归位修正为 `x += P`。
- `multipoint_eval_build` 的 `4 * m + 5` 越界改为 `4 * m + 4`。
- `matrix::operator~` 在不可逆时返回全零矩阵。
- `FFP` 的阶乘下标越界修正为 `i <= n`。
- `Dx` / `Integ` 恢复为返回新对象。
- 1.17.1 中 `~` 的 Newton 起始长度从 `1 << 20` 改为 `2` 后逐轮倍增，修复小规模求逆错误。
- 全局 scratch 数组收进 `struct` 或函数局部，修复同一进程内第二次调用出错的问题。

---

## 10. 注意事项

1. 固定模数下先调用 `Pre::initYG()`；运行时模数下先调用 `Pre::init(p)`。
2. 超出第 3 节的上限属于未定义行为，模板不会自动截断。
3. `Ln` 要求 `f[0] ≠ 0`；`Exp`、`Sin`、`Cos`、`Tan`、`Arc*` 要求 `f[0] = 0`；`~f`、`f^k` 要求 `f[0] ≠ 0`。
4. NTT 正向输出为位反转顺序；点乘后接 `NTT(false)`。
5. `P_FIXED=0` 时必须 `ATN=0`，且运行时不计算 `Y/I/_I_`；`Sin/Cos` 建议只在固定模数下使用。
6. 新增代码中的模运算统一写成 `mod(...)`。
7. 矩阵接口位于 `MATRIX` 命名空间，输入输出使用 `MATRIX::MIO`。
8. 输出使用自定义缓冲；与 `printf` / `cout` 混用时先 `flush()`。

---

## 附录：作者

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
