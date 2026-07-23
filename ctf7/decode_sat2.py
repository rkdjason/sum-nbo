#!/usr/bin/env python3
import argparse
import math
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


def parse_meta(path: Path):
    data = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        key, *vals = line.split()
        data[key] = vals
    return data


def parse_model(path: Path):
    text = path.read_text(errors="replace")
    if "UNSATISFIABLE" in text and "SATISFIABLE" not in text.replace("UNSATISFIABLE", ""):
        return None, "UNSAT"
    if "SATISFIABLE" not in text:
        return None, "UNKNOWN"
    model = set()
    for line in text.splitlines():
        if line.startswith("v ") or line == "v":
            for tok in line[1:].split():
                try:
                    x = int(tok)
                except ValueError:
                    continue
                if x > 0:
                    model.add(x)
    return model, "SAT"


def recover(bits, model):
    x = 0
    for i, raw in enumerate(bits):
        v = int(raw)
        b = 1 if v == -1 else (0 if v == 0 else int(v in model))
        x |= b << i
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("meta", type=Path)
    ap.add_argument("solver_output", type=Path)
    ap.add_argument("result", type=Path)
    args = ap.parse_args()
    meta = parse_meta(args.meta)
    model, status = parse_model(args.solver_output)
    if status != "SAT":
        args.result.write_text(f"status={status}\n")
        return 1
    p = recover(meta["P"], model)
    q = recover(meta["Q"], model)
    checks = {
        "product": p * q == N,
        "mask": (p & MASK) == PMASK,
        "p_bits": p.bit_length() == 1024,
        "q_bits": q.bit_length() == 1024,
        "gcd": math.gcd(p, q) == 1,
    }
    if not all(checks.values()):
        args.result.write_text(
            "status=INVALID_MODEL\n" + "\n".join(f"{k}={v}" for k,v in checks.items()) +
            f"\np=0x{p:x}\nq=0x{q:x}\n")
        return 2
    phi = (p - 1) * (q - 1)
    priv = pow(E, -1, phi)
    m = pow(CT, priv, N)
    plaintext = m.to_bytes((m.bit_length() + 7) // 8, "big")
    a = (p >> 150) & 15
    d = (p >> 920) & 15
    out = (
        "status=VERIFIED\n"
        f"a={a}\nd={d}\n"
        f"p={p}\nq={q}\n"
        f"p_hex=0x{p:x}\nq_hex=0x{q:x}\n"
        f"plaintext_hex={plaintext.hex()}\n"
        f"plaintext_repr={plaintext!r}\n"
        "verify_product=True\nverify_mask=True\n"
    )
    args.result.write_text(out)
    print(out)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
