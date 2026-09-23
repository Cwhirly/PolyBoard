# 版本管理速查（PolyBoard）

仓库：<https://github.com/Cwhirly/PolyBoard>（public，主分支 `main`，SSH 推送）
本机已经配置好：`git` + `gh`（已登录 `Cwhirly`）。下面所有命令都在项目目录里执行。

---

## 0. 一句话流程

> 新建版本文件 → `git add` → `git commit` → `git tag` → `git push` + `git push --tags`

每个"发布"就是一个 **提交 + 一个 tag**。tag 名 = 文件名去掉下划线和 `.cpp`：`1_17_2.cpp` → `1.17.2`。

---

## 1. 每次写完一个新版本怎么做

假设你刚写完 `1_18_0.cpp`：

```bash
cd ~/Desktop/Codes/cpp/Templetes/PolyBoard

git status                  # 先看一眼改了哪些文件（红色=没进版本库，绿色=已暂存）
git diff                    # 看具体改了啥（可选）

git add 1_18_0.cpp          # 把这个文件加入本次提交（只写文件名，别用 git add . ）
git commit -m "release: 1.18.0 <一句话说明这次加了/改了什么>"

git tag 1.18.0              # 给这个提交打个标签

git push                    # 推提交
git push --tags             # 推标签（新标签必须单独推，push 不会带上）
```

完事。刷新 GitHub 页面就能看到 tag。

> 如果 `1_18_0.cpp` 是在 `1_17_2.cpp` 基础上复制改的，记得 `git add` 的只有新文件；
> 老文件没动就不用管。README 如果也改了，一起 `git add README.md` 即可。

### 提交信息怎么写

看以前的记录就行，用的是「类型: 说明」：

| 前缀 | 用途 | 例子 |
| --- | --- | --- |
| `release:` | 发布一个新版本 | `release: 1.17.2 集合幂级数 sps 类 + 非齐次递推重写` |
| `feat:` | 加功能（不发版） | `feat: 加矩阵类` |
| `fix:` | 修 bug | `fix: 多点求值返回 x*f(x)+f(0)` |
| `style:` | 只改格式 | `style: 1_17_0 用 clang-format 整理` |
| `chore:` | 杂事（忽略文件、文档） | `chore: 忽略编译产物` |
| `docs:` | 只改文档 | `docs: README 补 sps 一节` |

---

## 2. 平时常用的查看命令

```bash
git log --oneline --decorate        # 提交 + tag 一览（现在最常用的一个）
git log --oneline -5                # 最近 5 条
git tag -l                          # 所有 tag
git tag -l --sort=-v:refname        # 按版本号从新到旧列 tag

git show 1.16.0 --stat              # 某个版本改了哪些文件
git show 1.16.0:1_16_0.cpp > /tmp/old.cpp   # 把某个版本的文件导出到别处看

git diff 1.17.0 1.17.1              # 比较两个版本（tag 之间可以直接 diff）
git diff HEAD                       # 工作区和最新提交的差别（提交前检查用）

git status                          # 当前状态（最常用）
git remote -v                       # 看远程仓库地址
gh repo view --web                  # 用浏览器打开仓库
```

---

## 3. 提交了但还没 push，想改/想撤

```bash
# 只改提交信息
git commit --amend -m "新的说明"

# 忘了 add 某个文件
git add 漏掉的文件 && git commit --amend --no-edit

# 完全撤掉最近一次提交，但保留文件内容（最常用）
git reset --soft HEAD~1

# 打错了 tag（还没 push）
git tag -d 1.18.0
git tag 1.18.0                      # 重新打

# 已经 push 过的 tag 想改（必须删远程的，别人拉过的仓库会报警）
git push origin :refs/tags/1.18.0 && git tag -d 1.18.0
```

**规矩：push 之后不要改历史。** 已经推上去的提交/标签就用新的提交去修正（`release: 1.18.1 修正 ...`），否则以后自己 `git pull` 都会冲突。

---

## 4. 给"当初忘了打标签"的老版本补标签

先找到这个文件是在哪个提交被加进来的：

```bash
git log --diff-filter=A --format='%h %ad %s' --date=short -- 1_10_0.cpp
# 83cb613 2025-07-14 release      <- 提交号在最左边
```

然后拿提交号打标签：

```bash
git tag 1.10.0 83cb613               # 给老提交补 tag
git push --tags                      # 推上去
```

反过来，想知道某个 tag 指向哪个提交：`git rev-parse 1.10.0`。

### 本仓库的历史版本 → 提交对照（已按此补好）

| tag | 提交 | 说明 |
| --- | --- | --- |
| 1.1.0 / 1.2.0 / 1.3.0 / 1.4.0 / 1.4.1 / 1.5.0 / 1.6.0 / 1.7.0 / 1.8.0 / 1.8.1 / 1.9.1 | 早期提交 | 建仓时就打好的；那时文件还叫 `Code.cpp`，每个 tag 里的 `Code.cpp` 就是那一版 |
| 1.10.0 / 1.11.0 / 1.12.0 / 1.13.0 | `83cb613` | 从这一版起改名叫 `1_x_y.cpp`；四个文件是同一个提交加的，所以 tag 指向同一个提交 |
| 1.14.0 / 1.15.0 / 1.15.1 / 1.15.2 | `6d55043` | 同上 |
| 1.16.0 | `2daef72` | `1_16_0.cpp` 首次出现 |
| 1.17.0 | `536f298` | 这个提交才把「子集反演」补进 `1_17_0.cpp` |
| 1.17.1 | `release: 1.17.1` 那个提交 | |
| 1.17.2 | `release: 1.17.2` 那个提交 | |

> 一个提交里加了多个版本文件时，多个 tag 会指向同一个提交，这是正常的（文件内容都在，回退到那个 tag 也能拿到对应文件）。

---

## 5. 编译产物不要进版本库

`.gitignore` 已经忽略了 `*.exe`、`*.out`、`a.out`、`.vscode/`、`tempCodeRunnerFile.cpp` 以及本项目里两个没后缀的编译产物 `1_17_0`、`1_17_2`。

今后编译建议显式指定输出名，省得再往 `.gitignore` 里加名字：

```bash
g++ 1_17_2.cpp -o main -std=c++17 -O2     # 或者 -o main.exe / -o a.out
```

如果又有新的产物没被忽略（`git status` 里以 `??` 出现），要么删掉，要么把它加进 `.gitignore` 再提交一次。

---

## 6. 想让 GitHub 上也有"发布页"（可选）

```bash
gh release create 1.17.2 --title "1.17.2" --notes "集合幂级数 sps 类；非齐次递推重写（P5808 满规模 3221ms → 552ms）"
```

`gh release create` 会自动用同名 tag（没有就新建），生成一个带说明和源码包的页面。

---

## 7. 万一手滑了

| 现象 | 处理 |
| --- | --- |
| 提交里混进了编译产物 | `git rm --cached 文件名` 然后 `git commit --amend`（文件还在硬盘上） |
| 想放弃某个文件的改动 | `git checkout -- 文件名`（**注意：改动会丢**） |
| 想看看某个旧版本长什么样 | `git show 1.16.0:1_16_0.cpp \| less` |
| push 报 `rejected`（远程更新了） | `git pull --rebase` 再 `git push` |
| 不知道自己在哪 | `git log --oneline -3` + `git status` |
