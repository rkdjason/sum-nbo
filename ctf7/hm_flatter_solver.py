#!/usr/bin/env python3
import argparse
import itertools
import math
import os
import re
import subprocess
import tempfile
import time
from pathlib import Path

N = int(
    "e505004fb5d34eb712d48ff4bbe8d27fc388133c6c0e734001061c0ee0a4edc6"
    "37c04fe8dd376185de8ba04d0ccdbabb93ab7c371b88d92e865eec42b028c61d"
    "d7004ebf2ebb5d69d0a09142be5c9de4da16e514eea318172ecda6cd192073eb"
    "afb1e02d522ec05334590ea6d75960c4937bf64f9700db177a4aa3da6aae6807"
    "e5e32c0d0e428a0db68d299f20c235d84ef459b0cf11828659c31663c9ea8204"
    "4b28152c89a9c36c3ec4303bd36664fd77fb02c58340bdae21120326d83fc017"
    "34bc90048dec9fe35f08c8fdc523abf84a91ec430f49567237c3153a2035ff62"
    "5613b6dc3e6cb14d50e18b8a79b25d678465b3ad02f5b7d818a1e2d635a0baf1", 16)
E = 65537
CT = int(
    "8919342826ef38215af31e00c9290c4c50ef9ff9e1afc59147fab5b096361035"
    "e85f5fc95b73b0697813b57b831a807d41bcbecde5b9e6639e2845b14e395ed0"
    "e5d995e63709ac0c5ee2337228ee76bcbad857b14904aa2e8e9997671908a634"
    "d0d1dda1d062ce7f2e3293ddec8f5cce26029292d594a062dcf317d2a8380f43"
    "d72551889efceb876c8945a50382272e76ed6b6fcdff160344e9e948e2b6e740"
    "e78bedf25f30e2c7eeb5f74686c8eadc29cea04ff08cfd86dfd3d2a1632bf04a"
    "d5cfa369892a2da40f0dc0098ce6b731d841aab3d0c8b78eb69c4625c47c4ad7"
    "158d49bb5d879581e02bc525abe47f39f699864bc5ce1de719430dae7aa5480b", 16)
MASK = int(
    "fffffffffffffffffffffffff0ffffffffffffffffffffffc00000000000fffe"
    "0000000000000000000003ffe00000000000000000ffffffffffffffffffff"
    "fffffffffffffffffffffff000000000000003ffe00000000000000000001ff"
    "fffffffffffffffffffffffffc3fffffffffffffffffffffffffffffffffffff", 16)
PMASK = int(
    "ffa360d46885c534d186538170633fafc2c0548a2e24a2c1c0000000000039e2"
    "0000000000000000000000a52000000000000000003e2de4c436d2ca740a6246"
    "99e1a1af94045c63261323c000000000000003bba00000000000000000000e5"
    "0b0bc2461fcbac0726360c2c0809450a9a892cbf1d98ceee48827591ccc593c9", 16)

X = 1 << 154
Y = 1 << 229
SX = 1 << 265
SY = 1 << 600
ALPHA = 1 << 335
INV_SX = pow(SX, -1, N)


def candidate_constant(a: int, d: int) -> int:
    base = PMASK | (a << 150) | (d << 920)
    return base + SX * X + SY * Y


def normalized_beta(c: int) -> int:
    b = (c * INV_SX) % N
    return b - N if b > N // 2 else b


def monomials(m):
    return [(i, j) for i in range(m + 1) for j in range(m - i + 1)]


def polynomial_powers(m, alpha, beta):
    powers = [{(0, 0): 1}]
    for _ in range(m):
        prev = powers[-1]
        cur = {}
        for (i, j), c in prev.items():
            cur[(i + 1, j)] = cur.get((i + 1, j), 0) + c
            cur[(i, j + 1)] = cur.get((i, j + 1), 0) + c * alpha
            cur[(i, j)] = cur.get((i, j), 0) + c * beta
        powers.append({k: v for k, v in cur.items() if v})
    return powers


def write_lattice(path: Path, m: int, t: int, beta: int):
    mons = monomials(m)
    idx = {z: i for i, z in enumerate(mons)}
    xp = [pow(X, i) for i in range(m + 1)]
    yp = [pow(Y, i) for i in range(m + 1)]
    np = [pow(N, max(t - k, 0)) for k in range(m + 1)]
    fp = polynomial_powers(m, ALPHA, beta)
    with path.open("w") as out:
        out.write("[")
        first = True
        for k in range(m + 1):
            mult = np[k]
            for shift_y in range(m - k + 1):
                row = [0] * len(mons)
                for (i, j), c in fp[k].items():
                    jj = j + shift_y
                    row[idx[(i, jj)]] = c * mult * xp[i] * yp[jj]
                if not first:
                    out.write("\n")
                first = False
                out.write("[" + " ".join(map(str, row)) + "]")
        out.write("\n]\n")
    return mons, xp, yp


def parse_lattice(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s in ("[", "]"):
            continue
        s = s.strip("[] ")
        if s:
            rows.append([int(x) for x in s.split()])
    return rows


def primitive(poly):
    vals = [abs(v) for v in poly.values() if v]
    if not vals:
        return {}
    g = 0
    for v in vals:
        g = math.gcd(g, v)
    if g > 1:
        poly = {k: v // g for k, v in poly.items()}
    # Canonical sign.
    lead = max(poly)
    if poly[lead] < 0:
        poly = {k: -v for k, v in poly.items()}
    return {k: v for k, v in poly.items() if v}


def divide_by_linear_x(poly, alpha, beta):
    work = dict(poly)
    quo = {}
    maxx = max((i for i, _ in work), default=0)
    for i in range(maxx, 0, -1):
        ys = sorted([j for (ii, j), c in work.items() if ii == i and c])
        for j in ys:
            c = work.pop((i, j), 0)
            if not c:
                continue
            quo[(i - 1, j)] = quo.get((i - 1, j), 0) + c
            work[(i - 1, j + 1)] = work.get((i - 1, j + 1), 0) - c * alpha
            work[(i - 1, j)] = work.get((i - 1, j), 0) - c * beta
            if work.get((i - 1, j + 1)) == 0:
                work.pop((i - 1, j + 1), None)
            if work.get((i - 1, j)) == 0:
                work.pop((i - 1, j), None)
    rem = {k: v for k, v in work.items() if v}
    if rem:
        return None
    return primitive(quo)


def reconstruct(rows, mons, xp, yp, p_lower, t, beta, max_polys=14):
    ranked = sorted(rows, key=lambda r: sum(z * z for z in r))
    threshold = pow(p_lower, 2 * t)
    result = []
    signatures = set()
    for row in ranked:
        nz = [(j, z) for j, z in enumerate(row) if z]
        if not nz:
            continue
        norm2 = sum(z * z for _, z in nz)
        if norm2 * len(nz) >= threshold:
            continue
        poly = {}
        okay = True
        for j, z in nz:
            i, k = mons[j]
            scale = xp[i] * yp[k]
            if z % scale:
                okay = False
                break
            c = z // scale
            if c:
                poly[(i, k)] = c
        if not okay or not poly:
            continue
        poly = primitive(poly)
        # Remove exact powers of the original linear polynomial. The quotient
        # also vanishes at the desired integer root whenever the product does.
        while True:
            q = divide_by_linear_x(poly, ALPHA, beta)
            if q is None or not q:
                break
            poly = q
        if len(poly) <= 1:
            continue
        # Avoid duplicate scalar-equivalent rows.
        lead = max(poly)
        sig = tuple(sorted((k, v % 1000003) for k, v in poly.items()))
        if sig in signatures:
            continue
        signatures.add(sig)
        result.append(poly)
        if len(result) >= max_polys:
            break
    return result


def poly_singular(poly, prime):
    terms = []
    for (i, j), c in sorted(poly.items(), reverse=True):
        c %= prime
        if not c:
            continue
        factors = []
        if c != 1 or (i == 0 and j == 0):
            factors.append(str(c))
        if i:
            factors.append("x" if i == 1 else f"x^{i}")
        if j:
            factors.append("y" if j == 1 else f"y^{j}")
        terms.append("*".join(factors) if factors else "1")
    return "+".join(terms) if terms else "0"


def is_probable_prime(n):
    if n < 2:
        return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for p in small:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in [2, 3, 5, 7, 11]:
        if a >= n:
            continue
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def primes_30(count=14):
    ans = []
    n = 1000000007
    while len(ans) < count:
        if is_probable_prime(n):
            ans.append(n)
        n += 2
    return ans


def parse_constant(s, prime):
    s = s.strip().replace(" ", "")
    if not s or s == "0":
        return 0
    if re.fullmatch(r"-?\d+", s):
        return int(s) % prime
    return None


def groebner_residue(polys, prime, singular="Singular", timeout=180):
    take = min(10, len(polys))
    declarations = [f"poly h{i}={poly_singular(polys[i], prime)};" for i in range(take)]
    ideal = ",".join(f"h{i}" for i in range(take))
    script = f"""
option(redSB);
ring R={prime},(x,y),lp;
{os.linesep.join(declarations)}
ideal I={ideal};
ideal G=std(I);
if (size(G)>0 && G[1]==1) {{ print(\"ONE\"); exit; }}
poly rx=reduce(x,G);
poly ry=reduce(y,G);
print(\"RX_BEGIN\"); print(rx); print(\"RX_END\");
print(\"RY_BEGIN\"); print(ry); print(\"RY_END\");
"""
    try:
        cp = subprocess.run([singular, "-q"], input=script, text=True,
                            capture_output=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None
    out = cp.stdout + "\n" + cp.stderr
    if "ONE" in out:
        return False
    mx = re.search(r"RX_BEGIN\s*(.*?)\s*RX_END", out, re.S)
    my = re.search(r"RY_BEGIN\s*(.*?)\s*RY_END", out, re.S)
    if not mx or not my:
        return None
    xr = parse_constant(mx.group(1).strip().splitlines()[-1], prime)
    yr = parse_constant(my.group(1).strip().splitlines()[-1], prime)
    if xr is None or yr is None:
        return None
    return xr, yr


def crt_pair(value, modulus, residue, prime):
    if modulus == 1:
        return residue, prime
    k = ((residue - value) % prime) * pow(modulus % prime, -1, prime) % prime
    return value + modulus * k, modulus * prime


def centered(v, mod):
    return v - mod if v > mod // 2 else v


def solve_candidate(a, d, work: Path, flatter, singular, m, t):
    c = candidate_constant(a, d)
    beta = normalized_beta(c)
    p_lower = PMASK | (a << 150) | (d << 920)
    lat = work / f"a{a}_m{m}t{t}.lat"
    red = work / f"a{a}_m{m}t{t}.red"
    meta = work / f"a{a}_m{m}t{t}.log"
    start = time.time()
    mons, xp, yp = write_lattice(lat, m, t, beta)
    with meta.open("w") as log:
        log.write(f"a={a} d={d} m={m} t={t} dim={len(mons)} beta_bits={abs(beta).bit_length()}\n")
        log.flush()
        try:
            cp = subprocess.run([flatter, "-delta", "0.99", str(lat), str(red)],
                                stdout=log, stderr=subprocess.STDOUT, timeout=2700)
        except subprocess.TimeoutExpired:
            log.write("FLATTER_TIMEOUT\n")
            return None
        if cp.returncode != 0 or not red.exists():
            log.write(f"FLATTER_FAIL rc={cp.returncode}\n")
            return None
    rows = parse_lattice(red)
    polys = reconstruct(rows, mons, xp, yp, p_lower, t, beta)
    with meta.open("a") as log:
        log.write(f"rows={len(rows)} polys={len(polys)} reduce_elapsed={time.time()-start:.2f}\n")
    # Keep reduced matrix only when debugging; it is large.
    try:
        lat.unlink()
        red.unlink()
    except OSError:
        pass
    if len(polys) < 2:
        return None

    xv = yv = 0
    mod = 1
    useful = 0
    for prime in primes_30(18):
        rr = groebner_residue(polys, prime, singular=singular)
        with meta.open("a") as log:
            log.write(f"prime={prime} residue={rr}\n")
        if rr is False:
            return None
        if rr is None:
            continue
        xr, yr = rr
        xv, newmod = crt_pair(xv, mod, xr, prime)
        yv, _ = crt_pair(yv, mod, yr, prime)
        mod = newmod
        useful += 1
        x0, y0 = centered(xv, mod), centered(yv, mod)
        if abs(x0) < X and abs(y0) < Y:
            p = c + SX * x0 + SY * y0
            if 1 < p < N and N % p == 0 and (p & MASK) == PMASK:
                q = N // p
                return p, q
        if mod > 4 * Y:
            break
    return None


def write_solution(path, p, q):
    phi = (p - 1) * (q - 1)
    priv = pow(E, -1, phi)
    mm = pow(CT, priv, N)
    pt = mm.to_bytes((mm.bit_length() + 7) // 8, "big")
    text = (
        "status=VERIFIED\n"
        f"p={p}\nq={q}\n"
        f"p_hex=0x{p:x}\nq_hex=0x{q:x}\n"
        f"plaintext_hex={pt.hex()}\nplaintext_repr={pt!r}\n"
        "verify_product=True\nverify_mask=True\n"
    )
    path.write_text(text)
    print(text, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("d", type=int)
    ap.add_argument("work", type=Path)
    ap.add_argument("solution", type=Path)
    ap.add_argument("--flatter", default="flatter")
    ap.add_argument("--singular", default="Singular")
    ap.add_argument("--params", default="22,6;23,7")
    ap.add_argument("--a-start", type=int, default=0)
    ap.add_argument("--a-end", type=int, default=16)
    args = ap.parse_args()
    assert 0 <= args.d < 16
    args.work.mkdir(parents=True, exist_ok=True)
    params = [tuple(map(int, x.split(","))) for x in args.params.split(";") if x]
    for m, t in params:
        for a in range(args.a_start, args.a_end):
            print(f"TRY d={args.d} a={a} m={m} t={t}", flush=True)
            ans = solve_candidate(a, args.d, args.work, args.flatter, args.singular, m, t)
            if ans:
                write_solution(args.solution, *ans)
                return 0
    args.solution.write_text(f"status=NOT_FOUND\nd={args.d}\n")
    return 1

if __name__ == "__main__":
    raise SystemExit(main())
