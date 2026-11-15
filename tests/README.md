# PolyBoard tests

数据在 `cases/<feature>/`，每个 case 是 `<name>.in` + `<name>.out`。带 `.err` 的 case 期望程序返回非零。

## 使用

```bash
tests/run.sh              # 全部
tests/run.sh conv inv     # 只测指定功能
tests/run.sh -l           # 列出功能
tests/run.sh -r           # 强制重新编译
tests/run.sh -q           # 只输出失败和汇总
```

脚本默认编译 `driver.cpp`，它包含 `../matrix.cpp`，所以 poly、sps、matrix 都可以测。默认模数是 `998244353`。

## 功能

| 名称 | 内容 |
| --- | --- |
| `conv` | NTT 卷积 |
| `inv` | 形式幂级数求逆 |
| `ln` | `Ln` |
| `exp` | `Exp` |
| `compose` | `UCPF::Compose`，含 `g[0] != 0` |
| `high` | `operator&` 的高半段 |
| `spsinv` | 集合幂级数求逆 |
| `spsexp` | 集合幂级数 `Exp` |
| `matmul` / `matinv` / `matdet` | 矩阵乘、逆、行列式 |
| `badinv` / `badln` / `badexp` / `badmat` | 期望错误 |

## 生成数据

```bash
python3 tests/gen.py
```

会把数据重新写入 `cases/`。随机种子固定，生成结果可复现。

## 加 case

在对应的 `cases/<feature>/` 放 `<name>.in` 和 `<name>.out` 即可；期望报错的 case 放一个空的 `<name>.err`。
