# 更新记录

> 每条先保留原来的主线描述，再补文件层面的具体变化。
> 1.1.0–1.9.1 只有一份 `Code.cpp`；1.10.0–1.13.0、1.14.0–1.15.2 分别指向同一个提交，是同一提交里的多份快照，不能用 `git diff tag1 tag2` 直接比较，下面的“文件”行按快照文件本身统计。

## 1.1.0

- **修复了原始版本的 $O(n\log^2n)$ 乘法逆。**
- 提供了更快的取模意义下加减法。
- 封装成更加完整的 `vector` 实现。
- 更新了 Cipolla 命名空间。

文件：`Code.cpp` 759 行。`poly` 内部改为 `vector<int>`，提供 `resize/shrink/rev/pb/tp/NTT`；`Pint` 提供 `addt/delt/add/del/tp/Vmax`；`Quad` 提供 `NTC/Q/q/Cipolla/work`；乘法逆改用牛顿迭代。

## 1.2.0

- **修复了直接取模导致 RE 的问题。**
- 更新了多项式除法与多项式取模。

文件：752 行（1.1.0 是 759 行）。去掉全局 `#define int long long`，乘法处显式 `1ll`；新增 `operator%`，移除半成品 `prefix/modfix`；`rev` 改为迭代器区间，`tp` 改为引用版。

## 1.3.0

- 更新了形式导数，形式积分，对数函数，指数函数。

文件：749 行。和 1.2.0 的代码差异只有 13 行，主要是打开 `using namespace UCPF` 并把 demo 从除法换成 `Sqrt`；`Dx/Integ/Ln/Exp` 的实现其实在 1.1.0 就已存在，这一版是把它们接到统一命名空间下。

## 1.4.0

- **修复了基于指数函数的多项式开根。**
- 更新了基于牛顿迭代的多项式开根，支持非首一多项式的开根。

文件：778 行。新增 `Sqrt_slow`（基于 `Exp` 的慢速开根），`Sqrt` 改为牛顿迭代版。

## 1.4.1

- 更新了 Sqrt、Pow 和三角函数。

文件：793 行，+175/−160。新增 `Sin/Cos/Tan` 与多项式快速幂 `operator^/^=`；移除临时的 `Sqrt_slow`。

## 1.5.0

- 更新了反三角函数（$\arcsin,\arccos,\arctan$）。

文件：829 行。新增 `ArcSin/ArcCos/ArcTan`。

## 1.6.0

- 修复了多项式多点求值。

文件：896 行（+67）。新增 `PMPE`、乘积树 `build`、`multipoint_evel`。

## 1.7.0

- 更新了多项式快速插值。

文件：945 行（+49）。新增 `PFI`、`Inv_evel`、`build1`。

## 1.8.0

- 更新了 Chirp-Z 变换。

文件：971 行（+26）。新增 `Chirp_Z/Ci2`。

## 1.8.1

- 将二进制 FWT 模板化，并封装按位与、按位或、按位异或卷积。

文件：1067 行（+115）。新增 `Fwtand/Fwtor/Fwtxor/AND/OR/XOR`。

## 1.9.0 / 1.9.1

- 更新了第二类斯特林数（行）。

文件：1062 行。新增 `Stirling::row`。仓库里没有单独的 `1.9.0` tag，实际 tag 是 `1.9.1`。

## 1.10.0

- 更新了第二类斯特林数（列）。
- 更新了适用于极大指数的多项式快速幂。

文件：`1_10_0.cpp` 1171 行。新增 `Stirling::column`、`Dev_Mul`、`Pow`；补了位移运算符 `<< >> <<= >>=`。

## 1.11.0

- 更新了第一类斯特林数（行）。

文件：1224 行（+53）。新增 `Trans`（后来用于普通/下降幂转换）。

## 1.12.0

- 更新了第一类斯特林数（列）。

文件：1256 行（+32）。新增 `Staling::column`；`Pow` 改名为 `Pow_For_Luogu`，修极大指数快速幂的入口。

## 1.13.0 / 1.14.0

- 更新了普通多项式与下降幂多项式的相互转化。

文件：1.13.0 是 1324 行，1.14.0 是 1312 行。新增 `FFPTP/PTFFP/Dev_Mul2`；1.14.0 去掉 `Dev_Mul2` 并重整了文件。

## 1.15.0

- 更新了基于 Fiduccia 算法的常系数齐次线性递推。

文件：1400 行（+88）。新增 `CCHLR::Fiduccia` 与 `__calc`。

## 1.15.1

- 将下降幂多项式封装成单独类。

文件：1504 行（+104）。新增 `ffp` 类，`pin/ppri` 增加 `ffp` 重载。

## 1.15.2

- 更新了基于 Bostan-Mori 算法的常系数齐次线性递推。

文件：1521 行（+17）。菜单从 `FFP to Class!` 改成 `Bostan-Mori`，但文件里实际只多了一份 `Fiduccia`（重复实现）；真正的 `Bostan_Mori` / `Bostan_Mori_CCHLR` 到 1.16.0 才出现。

## 1.16.0

- 更新了基于分块矩阵与 Fiduccia 算法的常系数非齐次线性递推。

文件：1669 行（+148）。新增 `Bostan_Mori`、`Bostan_Mori_CCHLR`、`Bostam_Mori_CCNHLR`、`__ksm`。

## 1.17.0

- 更新了子集卷积，优化与或卷积为 FMT。

文件：1387 行（1669 → 1387，−282）。新增 `Subset`、`Fmtor`；把命名空间/类结构定为 `CONSTS/Pre/Pint/Quad/POLY/UCPF/CCHLR/CCNHLR/...`，为 1.17.1 的重构打底。

## 1.17.1

- 重构了代码，改进了常数。

文件：1249 行（1387 → 1249）。13 项基准合计约 1.5 s → 0.73 s（2.0×）：
- NTT 免位反转（正向 DIF、逆向 DIT，`1/n` 用 `ny` 预表）；`resize` 去掉 `shrink_to_fit`；小卷积走朴素（`MUL_LIM=4096`）。
- FWT 去取模；`Trans`/`Chirp_Z` 的逐项快速幂改线性递推；`Subset` 内层提指针；补 `fac/ifac` 表。
- `poly` 缩到 24 字节；段树 scratch 按 `4m` 动态分配。
- 封装：`poly` 私有 `vector`，删 `n`/`shrink`；全局 scratch 收进 `struct`/函数；`ffp` 合并为 `using ffp = poly`；删掉所有 `#ifdef`。
- 接口：`poly f = pin(n)`；`PMPE(f, x)`/`PFI(x, y)` 返回值；`Quad::work(n)`；`Dx/Integ` 收 `const&`；`rev()/tp()` 不带参。
- 修 bug：PMPE 返回 `x_i f(x_i)+f(0)`、scratch 复用崩溃、`Integ` 长度、`Fiduccia` n=0、`CCNHLR` 的 U+FFFC、全 0 多项式取幂、`>>=` 越界。

## 1.17.2

- 更新了集合幂级数类，仅支持各项卷积，优化了非齐次线性递推。

文件：1339 行（1249 → 1339）。
- 新增 `sps` 类：FWT 从 `poly` 搬入，重载 `& | ^ *`，`And/Or/Xor/Subset` 保留为一行别名；`sps(const poly&)` 改成 `explicit`。
- CCNHLR 重写为一次 Bostan-Mori：生成函数并分式，五个部件全是卷积；P5808（k=m=30000, n=1e9）3221 ms → 552 ms（相对 1.17.0 5.8×，相对 1.17.1 3.8×）。
- `FFPTP` 内部仍然是多点求值（`Multipoint().eval(f, {0..d})`），这部分在 CCNHLR 里约占 14.8%。

## 1.18.0

- 模数与规模改成可注入：`#define POLYBOARD_P` / `POLYBOARD_N` 放在 include 之前即可，不必改模板本身。
- `LG = min(ceil(log2(2N-1)), v2(P-1))`，换成非 NTT 模数时旋转因子表自动缩到合法大小，不会越界。
- `Pre::Inv` 改用扩展欧几里得，任何模数下只要元素可逆就正确（不可逆返回 `0`）。
- 原根搜索限制在 1000 以内，避免合数模数下编译期死循环；`static_assert` 只在默认模数时启用。
- 可用性：默认 998244353 / 其它 NTT 素数全功能；非 NTT 素数（如 1e9+7）多项式不可用、矩阵可用；合数模数矩阵的 `gauss/inv` 要求主元可逆。
- 文件：1396 行（1339 → 1396，含 clang-format 的展开）。

## 2.1.0

- 更新了矩阵类，内含基本运算与高斯消元，行列式求值。

文件：`2_1_0.cpp` 1920 行，自带完整多项式代码 + `matrix`；`CP()` 的特征多项式实现有误（2.1.1 删除）。

## 2.1.1

- 修复了矩阵类的 bug：取相反数会输出 P、高斯消元躲零主元时换行错位、行列式没记换行符号。
- 更新了矩阵求逆，用 `~` 运算符；`*=` 换成 i-k-j 顺序 + 零元跳过。
- 2.x 改为第一行 `#include "1_17_2.cpp"` 复用 1.x 的全部多项式功能，工作区只保留这两份代码。

文件：`2_1_1.cpp` 202 行（不含被 include 的 1.x）。`~` 用增广矩阵 `[A|I]` 一次高斯-约当；`*=` 比 2.1.0 快约 2–3×；删掉 `CP()`；`1_17_2.cpp` 末尾加 `#if __INCLUDE_LEVEL__ == 0` 守卫。

## 2.1.2

- `gauss` 返回是否满行秩，新增 `inv(matrix&)` 返回是否可逆（失败时结果置全 0），`~` 不可逆时返回全 0 矩阵。
- `pin` 统一成返回值写法：`matrix f = pin(n, m);`（对照 `poly f = pin(n)`）；`ppri` 与多项式版同口径（截断/补 0、输出前 `tp`）。
- 消元主元要求可逆，配合扩展欧几里得的 `Inv`，合数模数下 `gauss` / `inv` 结论仍然可靠。

文件：`2_1_2.cpp` 219 行（202 → 219）。命名空间改成 `MATRIX` / `MIO`；`gauss`/`inv` 共用 `__rref`。

## 1.19.0

- 文件固定为 `poly.cpp` / `matrix.cpp`，版本只由 tag / release 区分；`matrix.cpp` 第一行 `#include "poly.cpp"`，以后不用再改名。
- 新增 `POLYBOARD_ATN` 编译期开关：`1` 走 NTT，`0` 走 MTT；默认 `1`。
- 新增 `POLYBOARD_P_FIXED`：默认 `1`；`0` 时 `P` 是运行时 `int`，必须 `ATN=0`，在 `main` 里调 `Pre::init(P)`。所有 `% P` 改成 `mod`：固定 P 下内联为 `% P`，运行时走 Barrett；`matrix.cpp` 同步改用 `mod`。
- `Pre::initYG` 在 `ATN=0` 时跳过 NTT 旋转因子表；`LG/NL` 在 `ATN=0` 时按 FFT 长度取；`P_FIXED=0` 时由 `Pre::init(P)` 计算 `B/LG/NL`、Barrett 和阶乘/逆元表。
- `poly` 的 `operator*=`、`operator&`、`operator~` 在 `ATN=0` 时走 `poly::MTT` / 通用牛顿迭代；`ATN=1` 保持原 NTT 路径。
- `poly::MTT` 是 5 次 FFT（`namespace FFT` 的 `run`，DIF 正 + DIT 逆），不再占位。
- 验证：固定 P 的 ATN=0/1、运行时 P=998244353/1000000007/12、n=1e5/6e5、矩阵固定/运行时 1e9+7 均通过；运行时与固定路径结果一致。
- 文件：`poly.cpp` 1591 行、`matrix.cpp` 219 行。

## 附：版本与 tag 对照

- 1.1.0–1.9.1：只有 `Code.cpp`。
- 1.9.0 没有独立 tag，仓库用 `1.9.1`。
- 1.10.0–1.13.0：同一个提交 `83cb613`，四份快照 `1_10_0.cpp`…`1_13_0.cpp`。
- 1.14.0–1.15.2：同一个提交 `6d55043`，四份快照 `1_14_0.cpp`…`1_15_2.cpp`。
- 1.16.0 起每个版本一个提交。
- 2.1.0 在 1.17.0 的提交里；2.1.1 起改成 include 复用。
- 1.19.0 起文件名固定为 `poly.cpp` / `matrix.cpp`，版本只由 tag 区分；旧 tag 里仍是 `1_x_y.cpp`。
