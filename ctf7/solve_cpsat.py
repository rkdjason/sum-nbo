#!/usr/bin/env python3
import argparse
import math
import os
from pathlib import Path
from ortools.sat.python import cp_model

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


def ceildiv(a, b):
    return (a + b - 1) // b


def fixed_q_data(d):
    mod = 1 << 265
    qlow = []
    for a in range(16):
        plow = (PMASK | (a << 150)) & (mod - 1)
        qlow.append((N * pow(plow, -1, mod)) & (mod - 1))
    all1024 = (1 << 1024) - 1
    dmask = 15 << 920
    remaining = (all1024 ^ MASK) & (all1024 ^ dmask)
    pmin = PMASK | (d << 920)
    pmax = pmin | remaining
    qmin, qmax = ceildiv(N, pmax), N // pmin
    fixed = {}
    for i in range(265):
        vals = {(x >> i) & 1 for x in qlow}
        if len(vals) == 1:
            fixed[i] = vals.pop()
    common = 0
    for i in range(1023, -1, -1):
        if ((qmin >> i) & 1) != ((qmax >> i) & 1):
            break
        fixed[i] = (qmin >> i) & 1
        common += 1
    return qlow, fixed, common, qmin, qmax


def make_number(model, prefix, fixed_bits, width=1024, limb_bits=8):
    bits = []
    for i in range(width):
        if i in fixed_bits:
            bits.append(int(fixed_bits[i]))
        else:
            bits.append(model.NewBoolVar(f"{prefix}b{i}"))
    limbs = []
    base = 1 << limb_bits
    for j in range((width + limb_bits - 1) // limb_bits):
        chunk = bits[j*limb_bits:(j+1)*limb_bits]
        if all(isinstance(x, int) for x in chunk):
            limbs.append(sum(int(x) << k for k, x in enumerate(chunk)))
        else:
            v = model.NewIntVar(0, base - 1, f"{prefix}L{j}")
            model.Add(v == sum((int(x) if isinstance(x, int) else x) << k for k, x in enumerate(chunk)))
            limbs.append(v)
    return bits, limbs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("d", type=int)
    ap.add_argument("result", type=Path)
    ap.add_argument("--seconds", type=float, default=21000)
    ap.add_argument("--workers", type=int, default=max(1, os.cpu_count() or 1))
    ap.add_argument("--limb-bits", type=int, default=8, choices=(4, 8, 16))
    args = ap.parse_args()
    d = args.d
    assert 0 <= d < 16
    model = cp_model.CpModel()
    qlow, qfixed, common, qmin, qmax = fixed_q_data(d)

    pfixed = {}
    for i in range(1024):
        if 920 <= i < 924:
            pfixed[i] = (d >> (i-920)) & 1
        elif (MASK >> i) & 1:
            pfixed[i] = (PMASK >> i) & 1
    pbits, plimbs = make_number(model, "p", pfixed, limb_bits=args.limb_bits)
    qbits, qlimbs = make_number(model, "q", qfixed, limb_bits=args.limb_bits)

    selectors = [model.NewBoolVar(f"a{a}") for a in range(16)]
    model.AddExactlyOne(selectors)
    for a, sel in enumerate(selectors):
        for j in range(4):
            b = pbits[150+j]
            val = (a >> j) & 1
            if isinstance(b, int):
                if b != val:
                    model.Add(sel == 0)
            else:
                model.Add(b == val).OnlyEnforceIf(sel)
        # 265 exact low q bits. Whole-limb equalities are stronger than individual bits.
        full_limbs = 265 // args.limb_bits
        rem = 265 % args.limb_bits
        mask_limb = (1 << args.limb_bits) - 1
        for j in range(full_limbs):
            val = (qlow[a] >> (j*args.limb_bits)) & mask_limb
            limb = qlimbs[j]
            if isinstance(limb, int):
                if limb != val:
                    model.Add(sel == 0)
            else:
                model.Add(limb == val).OnlyEnforceIf(sel)
        if rem:
            j = full_limbs
            for k in range(rem):
                b = qbits[j*args.limb_bits+k]
                val = (qlow[a] >> (j*args.limb_bits+k)) & 1
                if isinstance(b, int):
                    if b != val:
                        model.Add(sel == 0)
                else:
                    model.Add(b == val).OnlyEnforceIf(sel)

    base = 1 << args.limb_bits
    nlimbs = (1024 + args.limb_bits - 1) // args.limb_bits
    outlimbs = 2 * nlimbs
    Nchunks = [(N >> (args.limb_bits*k)) & (base-1) for k in range(outlimbs)]

    products = [[] for _ in range(outlimbs)]
    nonlinear = 0
    for i, px in enumerate(plimbs):
        if isinstance(px, int) and px == 0:
            continue
        for j, qx in enumerate(qlimbs):
            if isinstance(qx, int) and qx == 0:
                continue
            k = i+j
            if isinstance(px, int) and isinstance(qx, int):
                products[k].append(px*qx)
            elif isinstance(px, int):
                products[k].append(px*qx)
            elif isinstance(qx, int):
                products[k].append(qx*px)
            else:
                z = model.NewIntVar(0, (base-1)*(base-1), f"mul_{i}_{j}")
                model.AddMultiplicationEquality(z, [px, qx])
                products[k].append(z)
                nonlinear += 1

    carry = 0
    max_carry = 0
    carry_vars = []
    for k in range(outlimbs):
        max_terms = len(products[k]) * (base-1)*(base-1)
        next_max = (max_carry + max_terms) // base + 1
        nxt = model.NewIntVar(0, max(1, next_max), f"carry{k+1}") if k+1 < outlimbs else 0
        model.Add(carry + sum(products[k]) == Nchunks[k] + base*nxt)
        carry_vars.append(nxt)
        carry = nxt
        max_carry = next_max
    if not isinstance(carry, int):
        model.Add(carry == 0)

    primary = [x for x in selectors]
    # Branch on unknown p bits from both ends inward; q and carries should propagate.
    pvars = [b for b in pbits if not isinstance(b, int)]
    pvars.sort(key=lambda v: min(int(v.Name()[2:]), 1023-int(v.Name()[2:])))
    primary.extend(pvars)
    model.AddDecisionStrategy(primary, cp_model.CHOOSE_FIRST, cp_model.SELECT_MIN_VALUE)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = args.seconds
    solver.parameters.num_search_workers = args.workers
    solver.parameters.log_search_progress = True
    solver.parameters.cp_model_presolve = True
    solver.parameters.linearization_level = 2
    solver.parameters.random_seed = d + 1
    status = solver.Solve(model)
    status_name = solver.StatusName(status)
    print(f"d={d} limb_bits={args.limb_bits} common_high={common} nonlinear={nonlinear} status={status_name}")
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        args.result.write_text(f"status={status_name}\nd={d}\ncommon_high={common}\nnonlinear={nonlinear}\n")
        return 1

    def val(bits):
        x = 0
        for i,b in enumerate(bits):
            bitv = b if isinstance(b,int) else solver.Value(b)
            x |= int(bitv) << i
        return x
    p, q = val(pbits), val(qbits)
    if p*q != N or (p&MASK) != PMASK:
        args.result.write_text(f"status=INVALID\np={p}\nq={q}\n")
        return 2
    priv = pow(E, -1, (p-1)*(q-1))
    m = pow(CT, priv, N)
    pt = m.to_bytes((m.bit_length()+7)//8, "big")
    out = (
        "status=VERIFIED\n"
        f"p={p}\nq={q}\n"
        f"p_hex=0x{p:x}\nq_hex=0x{q:x}\n"
        f"plaintext_hex={pt.hex()}\nplaintext_repr={pt!r}\n"
        "verify_product=True\nverify_mask=True\n"
    )
    args.result.write_text(out)
    print(out)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
