#!/usr/bin/env sage -python
from sage.all import *
import argparse
import itertools
import os
import subprocess
import sys
import time

N = Integer('e505004fb5d34eb712d48ff4bbe8d27fc388133c6c0e734001061c0ee0a4edc637c04fe8dd376185de8ba04d0ccdbabb93ab7c371b88d92e865eec42b028c61dd7004ebf2ebb5d69d0a09142be5c9de4da16e514eea318172ecda6cd192073ebafb1e02d522ec05334590ea6d75960c4937bf64f9700db177a4aa3da6aae6807e5e32c0d0e428a0db68d299f20c235d84ef459b0cf11828659c31663c9ea82044b28152c89a9c36c3ec4303bd36664fd77fb02c58340bdae21120326d83fc01734bc90048dec9fe35f08c8fdc523abf84a91ec430f49567237c3153a2035ff625613b6dc3e6cb14d50e18b8a79b25d678465b3ad02f5b7d818a1e2d635a0baf1', 16)
e = Integer(65537)
ct = Integer('8919342826ef38215af31e00c9290c4c50ef9ff9e1afc59147fab5b096361035e85f5fc95b73b0697813b57b831a807d41bcbecde5b9e6639e2845b14e395ed0e5d995e63709ac0c5ee2337228ee76bcbad857b14904aa2e8e9997671908a634d0d1dda1d062ce7f2e3293ddec8f5cce26029292d594a062dcf317d2a8380f43d72551889efceb876c8945a50382272e76ed6b6fcdff160344e9e948e2b6e740e78bedf25f30e2c7eeb5f74686c8eadc29cea04ff08cfd86dfd3d2a1632bf04ad5cfa369892a2da40f0dc0098ce6b731d841aab3d0c8b78eb69c4625c47c4ad7158d49bb5d879581e02bc525abe47f39f699864bc5ce1de719430dae7aa5480b', 16)
MASK = Integer('fffffffffffffffffffffffff0ffffffffffffffffffffffc00000000000fffe0000000000000000000003ffe00000000000000000fffffffffffffffffffffffffffffffffffffffffffff000000000000003ffe000000000000000000001fffffffffffffffffffffffffffc3fffffffffffffffffffffffffffffffffffff', 16)
PM = Integer('ffa360d46885c534d186538170633fafc2c0548a2e24a2c1c0000000000039e20000000000000000000000a52000000000000000003e2de4c436d2ca740a624699e1a1af94045c63261323c000000000000003bba000000000000000000000e50b0bc2461fcbac0726360c2c0809450a9a892cbf1d98ceee48827591ccc593c9', 16)

OX, LX = 265, 155
OY, LY = 600, 230
CX, CY = Integer(1) << (LX - 1), Integer(1) << (LY - 1)
BX, BY = CX, CY
INV = inverse_mod(Integer(1) << OX, N)
CYCOEF = ((Integer(1) << OY) * INV) % N
if CYCOEF > N // 2:
    CYCOEF -= N

ZZxy = PolynomialRing(ZZ, names=('x', 'y'), order='invlex')
x, y = ZZxy.gens()


def centered_constant(a, d):
    base = PM | (Integer(a) << 150) | (Integer(d) << 920)
    shifted = base + (CX << OX) + (CY << OY)
    b = (shifted * INV) % N
    if b > N // 2:
        b -= N
    return base, b


def monomial_exponents(m):
    return [(i, j) for i in range(m + 1) for j in range(m - i + 1)]


def powers_of_linear(m, ay, b):
    powers = [{(0, 0): Integer(1)}]
    for _ in range(m):
        prev = powers[-1]
        out = {}
        for (i, j), c in prev.items():
            out[(i + 1, j)] = out.get((i + 1, j), 0) + c
            out[(i, j + 1)] = out.get((i, j + 1), 0) + c * ay
            out[(i, j)] = out.get((i, j), 0) + c * b
        powers.append(out)
    return powers


def build_lattice(a, d, m, t):
    base, b = centered_constant(a, d)
    exps = monomial_exponents(m)
    pos = {z: i for i, z in enumerate(exps)}
    scales = [(BX ** i) * (BY ** j) for i, j in exps]
    powers = powers_of_linear(m, CYCOEF, b)
    rows = []
    for k in range(m + 1):
        nk = N ** max(t - k, 0)
        for s in range(m - k + 1):
            row = [Integer(0)] * len(exps)
            for (i, j), c in powers[k].items():
                jj = j + s
                col = pos[(i, jj)]
                row[col] = c * nk * scales[col]
            rows.append(row)
    assert len(rows) == len(exps)
    return base, b, exps, scales, rows


def write_fplll_matrix(path, rows):
    with open(path, 'w', encoding='ascii', buffering=1024 * 1024) as f:
        f.write('[\n')
        for row in rows:
            f.write('[' + ' '.join(str(v) for v in row) + ']\n')
        f.write(']\n')


def read_fplll_matrix(path):
    rows = []
    with open(path, 'r', encoding='ascii', buffering=1024 * 1024) as f:
        for line in f:
            s = line.strip()
            if not s or s in ('[', ']'):
                continue
            if s.startswith('[['):
                s = s[1:]
            if s.endswith(']]'):
                s = s[:-1]
            if s.startswith('['):
                s = s[1:]
            if s.endswith(']'):
                s = s[:-1]
            if s.strip():
                rows.append([Integer(z) for z in s.split()])
    return rows


def reduce_with_flatter(rows, delta, workdir):
    inp = os.path.join(workdir, 'lattice.in')
    out = os.path.join(workdir, 'lattice.out')
    write_fplll_matrix(inp, rows)
    cmd = ['flatter', '-delta', str(delta), inp, out]
    start = time.time()
    cp = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
    elapsed = time.time() - start
    if cp.returncode != 0:
        raise RuntimeError('flatter failed (%d):\n%s' % (cp.returncode, cp.stdout[-8000:]))
    reduced = read_fplll_matrix(out)
    if len(reduced) != len(rows):
        raise RuntimeError('reduced matrix row mismatch: %d != %d' % (len(reduced), len(rows)))
    return reduced, elapsed


def primitive_poly(poly):
    if poly == 0:
        return poly
    cont = gcd([Integer(c) for c in poly.coefficients()])
    if cont not in (0, 1, -1):
        poly = poly // cont
    if poly.lc() < 0:
        poly = -poly
    return poly


def reconstruct_polynomials(reduced, exps, scales, fpoly, t, max_polys=18):
    threshold_sq = Integer(1) << (2 * 1023 * t)
    candidates = []
    for row in reduced:
        nz = [(v, exp, sc) for v, exp, sc in zip(row, exps, scales) if v]
        if not nz:
            continue
        norm2 = sum(v * v for v, _, _ in nz)
        weight = len(nz)
        if norm2 * weight >= threshold_sq:
            continue
        poly = ZZxy(0)
        valid = True
        for v, (i, j), sc in nz:
            if v % sc:
                valid = False
                break
            poly += (v // sc) * (x ** i) * (y ** j)
        if not valid or poly == 0:
            continue
        poly = primitive_poly(poly)
        while poly != 0:
            qpoly, rem = poly.quo_rem(fpoly)
            if rem != 0:
                break
            poly = primitive_poly(qpoly)
        if poly == 0 or poly.is_constant():
            continue
        if any(poly == old for _, old in candidates):
            continue
        candidates.append((norm2, poly))
    candidates.sort(key=lambda z: z[0])
    return [p for _, p in candidates[:max_polys]]


def linear_solution_from_groebner(G, prime):
    eqs = []
    for g in G:
        if g == 1:
            return None
        if g.total_degree() <= 1 and not g.is_constant():
            eqs.append(g)
    if len(eqs) < 2:
        return None
    F = GF(prime)
    for g1, g2 in itertools.combinations(eqs, 2):
        a1 = F(g1.monomial_coefficient(g1.parent().gen(0)))
        b1 = F(g1.monomial_coefficient(g1.parent().gen(1)))
        c1 = F(g1.constant_coefficient())
        a2 = F(g2.monomial_coefficient(g2.parent().gen(0)))
        b2 = F(g2.monomial_coefficient(g2.parent().gen(1)))
        c2 = F(g2.constant_coefficient())
        det = a1 * b2 - a2 * b1
        if det == 0:
            continue
        rx = (-c1 * b2 + c2 * b1) / det
        ry = (-a1 * c2 + a2 * c1) / det
        if all(g(rx, ry) == 0 for g in G):
            return int(rx), int(ry)
    return None


def modular_unique_root(polys, prime, take=8):
    F = GF(prime)
    R = PolynomialRing(F, names=('x', 'y'), order='lex')
    mapped = [R(p) for p in polys[:take]]
    selected = []
    for p in mapped:
        if p == 0 or p.is_constant():
            continue
        selected.append(p)
        if len(selected) < 2:
            continue
        try:
            G = ideal(selected).groebner_basis()
        except Exception:
            continue
        if any(g == 1 for g in G):
            return None
        sol = linear_solution_from_groebner(G, prime)
        if sol is not None:
            return sol
    return None


def crt_signed(residues, moduli, bound):
    value = Integer(crt(residues, moduli))
    modulus = prod(moduli)
    if modulus <= 2 * bound:
        return None
    if value > modulus // 2:
        value -= modulus
    if not (-bound <= value < bound):
        return None
    return value


def recover_root(polys, prime_count=10, take=10):
    rx, ry, mods = [], [], []
    prime = next_prime(Integer(1000000007))
    attempts = 0
    while len(mods) < prime_count and attempts < prime_count * 5:
        attempts += 1
        sol = modular_unique_root(polys, prime, take=take)
        if sol is not None:
            sx, sy = sol
            rx.append(Integer(sx))
            ry.append(Integer(sy))
            mods.append(Integer(prime))
            xv = crt_signed(rx, mods, BX)
            yv = crt_signed(ry, mods, BY)
            if xv is not None and yv is not None:
                if all(poly(xv, yv) == 0 for poly in polys):
                    return xv, yv, mods
        prime = next_prime(prime + 1000)
    return None


def verify_and_decrypt(base, xc, yc):
    ux = xc + CX
    uy = yc + CY
    if not (0 <= ux < (Integer(1) << LX) and 0 <= uy < (Integer(1) << LY)):
        return None
    p = base + (ux << OX) + (uy << OY)
    if p <= 1 or N % p:
        return None
    q = N // p
    if p * q != N or (p & MASK) != PM:
        return None
    phi = (p - 1) * (q - 1)
    priv = inverse_mod(e, phi)
    msg = power_mod(ct, priv, N)
    raw = int(msg).to_bytes((int(msg).bit_length() + 7) // 8, 'big')
    return p, q, raw


def solve_candidate(a, d, m, t, delta, work_root, max_polys, prime_count):
    tag = 'a%02d_d%02d' % (a, d)
    workdir = os.path.join(work_root, tag)
    os.makedirs(workdir, exist_ok=True)
    t0 = time.time()
    base, b, exps, scales, rows = build_lattice(a, d, m, t)
    print('CANDIDATE', tag, 'dim', len(rows), 'build_sec', round(time.time() - t0, 3), flush=True)
    reduced, redsec = reduce_with_flatter(rows, delta, workdir)
    print('REDUCED', tag, 'sec', round(redsec, 3), flush=True)
    fpoly = x + CYCOEF * y + b
    polys = reconstruct_polynomials(reduced, exps, scales, fpoly, t, max_polys=max_polys)
    print('POLYS', tag, len(polys), [(p.total_degree(), len(p.monomials())) for p in polys[:8]], flush=True)
    if len(polys) < 2:
        return None
    root = recover_root(polys, prime_count=prime_count, take=min(max_polys, 12))
    if root is None:
        print('NO_ROOT', tag, flush=True)
        return None
    xc, yc, mods = root
    print('ROOT', tag, 'xc', xc, 'yc', yc, 'crt_bits', prod(mods).nbits(), flush=True)
    ans = verify_and_decrypt(base, xc, yc)
    if ans is None:
        print('ROOT_FAILED_VERIFY', tag, flush=True)
        return None
    p, q, raw = ans
    print('FOUND_P=0x%x' % p, flush=True)
    print('FOUND_Q=0x%x' % q, flush=True)
    print('FOUND_P_DEC=%d' % p, flush=True)
    print('FOUND_Q_DEC=%d' % q, flush=True)
    print('PLAINTEXT_HEX=%s' % raw.hex(), flush=True)
    print('PLAINTEXT_REPR=%r' % raw, flush=True)
    with open(os.path.join(work_root, 'FOUND.txt'), 'w') as f:
        f.write('p=0x%x\nq=0x%x\nplaintext_hex=%s\nplaintext_repr=%r\n' % (p, q, raw.hex(), raw))
    return ans


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--d', type=int, required=True)
    ap.add_argument('--a-start', type=int, default=0)
    ap.add_argument('--a-end', type=int, default=16)
    ap.add_argument('--m', type=int, default=24)
    ap.add_argument('--t', type=int, default=7)
    ap.add_argument('--delta', type=float, default=0.99)
    ap.add_argument('--max-polys', type=int, default=18)
    ap.add_argument('--prime-count', type=int, default=10)
    ap.add_argument('--work', default='work')
    args = ap.parse_args()
    assert 0 <= args.d < 16 and 0 <= args.a_start <= args.a_end <= 16
    print('CONFIG', vars(args), flush=True)
    print('MASK_KNOWN', MASK.popcount(), 'MASK_UNKNOWN', 1024 - MASK.popcount(), flush=True)
    print('BOUNDS_BITS', BX.nbits() - 1, BY.nbits() - 1, 'sum', (BX.nbits() - 1) + (BY.nbits() - 1), flush=True)
    for a in range(args.a_start, args.a_end):
        try:
            ans = solve_candidate(a, args.d, args.m, args.t, args.delta, args.work, args.max_polys, args.prime_count)
            if ans is not None:
                return 0
        except Exception as exc:
            print('ERROR a=%d d=%d: %s: %s' % (a, args.d, type(exc).__name__, exc), flush=True)
            import traceback
            traceback.print_exc()
    print('NOT_FOUND d=%d range=%d:%d' % (args.d, args.a_start, args.a_end), flush=True)
    return 1


if __name__ == '__main__':
    sys.exit(main())
