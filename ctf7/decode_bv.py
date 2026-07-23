#!/usr/bin/env python3
import argparse
import re
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


def candidates(text):
    vals = []
    for h in re.findall(r"#x([0-9a-fA-F]+)", text):
        if len(h) >= 250:
            vals.append(int(h, 16))
    for b in re.findall(r"#b([01]+)", text):
        if len(b) >= 1000:
            vals.append(int(b, 2))
    for line in text.splitlines():
        fields = line.strip().split()
        for f in fields:
            if len(f) >= 1000 and set(f) <= {"0", "1"}:
                vals.append(int(f, 2))
    seen = set()
    for x in vals:
        if x not in seen:
            seen.add(x)
            yield x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("solver_output", type=Path)
    ap.add_argument("result", type=Path)
    args = ap.parse_args()
    text = args.solver_output.read_text(errors="replace")
    vals = list(candidates(text))
    p = q = None
    for x in vals:
        if x > 1 and N % x == 0:
            y = N // x
            if (x & MASK) == PMASK:
                p, q = x, y
                break
            if (y & MASK) == PMASK:
                p, q = y, x
                break
    if p is None:
        status = "UNSAT" if re.search(r"\bunsat\b", text, re.I) else "UNKNOWN_OR_UNPARSED"
        args.result.write_text(f"status={status}\nvalues_seen={len(vals)}\n")
        return 1
    assert p * q == N and (p & MASK) == PMASK
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
