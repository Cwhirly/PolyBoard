import os, random, itertools
P = 998244353
ROOT = os.path.dirname(os.path.abspath(__file__))
C = os.path.join(ROOT, "cases")

def wr(feature, name, inp, out=None, err=False):
    d = os.path.join(C, feature)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, name + ".in"), "w") as f:
        f.write(inp)
    if out is not None:
        with open(os.path.join(d, name + ".out"), "w") as f:
            f.write(out)
    if err:
        open(os.path.join(d, name + ".err"), "w").close()

def conv(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] = (c[i + j] + x * y) % P
    return c

def inv(f, n):
    g = [0] * n
    g[0] = pow(f[0], P - 2, P)
    for i in range(1, n):
        s = 0
        for j in range(1, i + 1):
            s = (s + f[j] * g[i - j]) % P
        g[i] = (-s * g[0]) % P
    return g

def ln(f, n):
    g = inv(f, n)
    d = [(i * f[i]) % P for i in range(1, n)]
    q = conv(d, g)[:n - 1]
    r = [0] * n
    for i in range(1, n):
        r[i] = q[i - 1] * pow(i, P - 2, P) % P
    return r

def exp(f, n):
    g = [0] * n
    g[0] = 1
    for i in range(1, n):
        s = 0
        for j in range(1, i + 1):
            s = (s + j * f[j] * g[i - j]) % P
        g[i] = s * pow(i, P - 2, P) % P
    return g

def compose(f, g, n):
    r = [0] * n
    for c in reversed(f):
        r = conv(r, g)[:n]
        r[0] = (r[0] + c) % P
    return r

def fmt(a):
    return " ".join(map(str, a)) + "\n"

def subset_inv(f, k):
    n = 1 << k
    g = [0] * n
    g[0] = pow(f[0], P - 2, P)
    for s in range(1, n):
        z = 0
        t = s
        while t:
            z = (z + f[t] * g[s ^ t]) % P
            t = (t - 1) & s
        g[s] = (-z * g[0]) % P
    return g

def subset_exp(f, k):
    n = 1 << k
    g = [0] * n
    g[0] = 1
    for s in range(1, n):
        m = s & -s
        z = 0
        t = s
        while t:
            if t & m:
                z = (z + f[t] * g[s ^ t]) % P
            t = (t - 1) & s
        g[s] = z
    return g

def mat_mul(a, b):
    n = len(a)
    c = [[0] * n for _ in range(n)]
    for i in range(n):
        for k in range(n):
            for j in range(n):
                c[i][j] = (c[i][j] + a[i][k] * b[k][j]) % P
    return c

def mat_det(a):
    a = [row[:] for row in a]
    n = len(a)
    r = 1
    for i in range(n):
        p = next((j for j in range(i, n) if a[j][i]), None)
        if p is None:
            return 0
        if p != i:
            a[i], a[p] = a[p], a[i]
            r = -r
        r = r * a[i][i] % P
        iv = pow(a[i][i], P - 2, P)
        for j in range(i + 1, n):
            if not a[j][i]:
                continue
            x = a[j][i] * iv % P
            for k in range(i, n):
                a[j][k] = (a[j][k] - x * a[i][k]) % P
    return r % P

def mat_inv(a):
    n = len(a)
    m = [[0] * (2 * n) for _ in range(n)]
    for i in range(n):
        for j in range(n):
            m[i][j] = a[i][j]
        m[i][n + i] = 1
    for i in range(n):
        p = next((j for j in range(i, n) if m[j][i]), None)
        if p is None:
            return None
        m[i], m[p] = m[p], m[i]
        iv = pow(m[i][i], P - 2, P)
        for j in range(2 * n):
            m[i][j] = m[i][j] * iv % P
        for j in range(n):
            if j == i or not m[j][i]:
                continue
            x = m[j][i]
            for k in range(2 * n):
                m[j][k] = (m[j][k] - x * m[i][k]) % P
    return [row[n:] for row in m]

def mat_fmt(a):
    return "".join(" ".join(map(str, row)) + "\n" for row in a)

random.seed(20261115)

for n, m, tag in [(1,1,"tiny"),(2,3,"small"),(17,19,"mid"),(257,300,"big"),(4096,4096,"large")]:
    f = [random.randrange(P) for _ in range(n)]
    g = [random.randrange(P) for _ in range(m)]
    wr("conv", f"{tag}_{n}_{m}", f"{n} {m}\n{fmt(f)}{fmt(g)}", fmt(conv(f,g)))

for n, tag in [(1,"tiny"),(2,"small"),(7,"small"),(64,"mid"),(300,"big")]:
    f = [random.randrange(P) for _ in range(n)]
    f[0] = random.randrange(1, P)
    wr("inv", f"{tag}_{n}", f"{n}\n{fmt(f)}", fmt(inv(f,n)))

for n, tag in [(1,"tiny"),(2,"small"),(8,"small"),(64,"mid")]:
    f = [1] + [random.randrange(P) for _ in range(n - 1)]
    wr("ln", f"{tag}_{n}", f"{n}\n{fmt(f)}", fmt(ln(f,n)))

for n, tag in [(1,"tiny"),(2,"small"),(8,"small"),(64,"mid")]:
    f = [0] + [random.randrange(P) for _ in range(n - 1)]
    wr("exp", f"{tag}_{n}", f"{n}\n{fmt(f)}", fmt(exp(f,n)))

for n, m, tag in [(1,1,"tiny"),(2,3,"small"),(5,7,"small"),(12,20,"mid"),(30,50,"big")]:
    f = [random.randrange(P) for _ in range(n)]
    g = [random.randrange(P) for _ in range(m)]
    wr("compose", f"{tag}_{n}_{m}", f"{n} {m}\n{fmt(f)}{fmt(g)}", fmt(compose(f,g,n)))
    f = [random.randrange(P) for _ in range(n)]
    g = [random.randrange(P) for _ in range(m)]
    g[0] = 1
    wr("compose", f"{tag}_{n}_{m}_g0", f"{n} {m}\n{fmt(f)}{fmt(g)}", fmt(compose(f,g,n)))

for n, m, tag in [(2,1,"tiny"),(8,3,"small"),(64,15,"mid")]:
    f = [random.randrange(P) for _ in range(n)]
    g = [random.randrange(P) for _ in range(m)]
    full = conv(f,g)
    wr("high", f"{tag}_{n}_{m}", f"{n} {m}\n{fmt(f)}{fmt(g)}", fmt(full[m-1:n]))

for k, tag in [(1,"small"),(2,"small"),(3,"mid")]:
    n = 1 << k
    f = [random.randrange(P) for _ in range(n)]
    f[0] = random.randrange(1, P)
    wr("spsinv", f"k{k}", f"{k}\n{fmt(f)}", fmt(subset_inv(f,k)))

for k, tag in [(1,"small"),(2,"small"),(3,"mid")]:
    n = 1 << k
    f = [0] + [random.randrange(P) for _ in range(n - 1)]
    wr("spsexp", f"k{k}", f"{k}\n{fmt(f)}", fmt(subset_exp(f,k)))

for n, tag in [(1,"tiny"),(2,"small"),(4,"small"),(6,"mid")]:
    a = [[random.randrange(P) for _ in range(n)] for _ in range(n)]
    b = [[random.randrange(P) for _ in range(n)] for _ in range(n)]
    wr("matmul", f"{tag}_{n}", f"{n}\n{mat_fmt(a)}{mat_fmt(b)}", mat_fmt(mat_mul(a,b)))
    wr("matdet", f"{tag}_{n}", f"{n}\n{mat_fmt(a)}", f"{mat_det(a)}\n")
    c = mat_inv(a)
    if c is not None:
        wr("matinv", f"{tag}_{n}", f"{n}\n{mat_fmt(a)}", mat_fmt(c))

# expected errors
wr("badinv", "zero", "3\n0 1 2\n", err=True)
wr("badln", "zero", "3\n0 1 2\n", err=True)
wr("badexp", "nonzero", "3\n1 0 0\n", err=True)
wr("badmat", "singular", "2\n1 2\n2 4\n", err=True)
print("generated")

# corner cases
wr("conv", "corner_zero", "2 3\n0 0\n0 1 2\n", fmt(conv([0,0],[0,1,2])))
wr("inv", "corner_const", "3\n5 0 0\n", fmt(inv([5,0,0],3)))
wr("ln", "corner_one", "1\n1\n", "0\n")
wr("exp", "corner_zero", "3\n0 0 0\n", fmt(exp([0,0,0],3)))
wr("compose", "corner_const_f", "3 4\n7 0 0\n1 2 3 4\n", fmt(compose([7,0,0],[1,2,3,4],3)))
wr("compose", "corner_zero_g", "3 1\n1 2 3\n0\n", fmt(compose([1,2,3],[0],3)))
wr("spsinv", "corner_one", "2\n1 0 0 0\n", fmt(subset_inv([1,0,0,0],2)))
wr("spsexp", "corner_zero", "2\n0 0 0 0\n", fmt(subset_exp([0,0,0,0],2)))
wr("high", "corner_equal", "8 5\n1 2 3 4 5 6 7 8\n6 7 8 9 10\n", fmt(conv([1,2,3,4,5,6,7,8],[6,7,8,9,10])[4:8]))
