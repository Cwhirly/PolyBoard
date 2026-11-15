# PolyBoard tests

测试数据不全部提交，`gen.py` 在 `build/cases` 下生成。`build/` 已被忽略。

## 功能测试

```bash
tests/run.sh              # 全部
tests/run.sh conv inv     # 只测指定功能
tests/run.sh -l           # 列出功能
tests/run.sh -r           # 强制重新编译
tests/run.sh -g           # 重新生成数据
tests/run.sh -q           # 只输出失败和汇总
```

脚本会按需运行 `python3 gen.py`，生成 `cases/<feature>/<name>.in/.out`。带 `.err` 的 case 期望程序返回非零。

| 名称 | 内容 |
| --- | --- |
| `conv` | NTT 卷积 |
| `inv` | 形式幂级数求逆 |
| `ln` | `Ln` |
| `exp` | `Exp` |
| `compose` | `UCPF::Compose` |
| `high` | `operator&` |
| `spsinv` / `spsexp` | 集合幂级数求逆、`Exp` |
| `matmul` / `matinv` / `matdet` | 矩阵 |
| `badinv` / `badln` / `badexp` / `badmat` | 期望错误 |
| `shift` | 连续点值平移 |
| `revert` | 复合逆 |
| `prop` | 恒等式/性质测试 |

## 性质测试

`prop` 不存固定答案，由 driver 内部检查：

- `Sqrt(f)^2 = f`
- `Sin^2 + Cos^2 = 1`
- `Trans` 对二项式展开
- `(f/g)*g + f%g = f`
- 位运算卷积对暴力
- 子集卷积对暴力
- 多点求值再插值回原多项式
- 连续点值平移对直接求值
- 复合逆的 `f(g)=x`

## 性能测试

```bash
tests/bench.sh
tests/bench.sh -n 20000 -r 5 conv exp compose
```

每个操作跑若干次，输出 `op n rep time hash`。时间依赖机器，只适合本机前后对比。

## 加 case

在 `gen.py` 里加生成逻辑，或在 `build/cases/<feature>/` 手动放 `.in/.out`；`.err` 表示期望失败。
