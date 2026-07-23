#!/usr/bin/env python3
import argparse

N = int(
    "e505004fb5d34eb712d48ff4bbe8d27fc388133c6c0e734001061c0ee0a4edc6"
    "37c04fe8dd376185de8ba04d0ccdbabb93ab7c371b88d92e865eec42b028c61d"
    "d7004ebf2ebb5d69d0a09142be5c9de4da16e514eea318172ecda6cd192073eb"
    "afb1e02d522ec05334590ea6d75960c4937bf64f9700db177a4aa3da6aae6807"
    "e5e32c0d0e428a0db68d299f20c235d84ef459b0cf11828659c31663c9ea8204"
    "4b28152c89a9c36c3ec4303bd36664fd77fb02c58340bdae21120326d83fc017"
    "34bc90048dec9fe35f08c8fdc523abf84a91ec430f49567237c3153a2035ff62"
    "5613b6dc3e6cb14d50e18b8a79b25d678465b3ad02f5b7d818a1e2d635a0baf1", 16)
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


def bvhex(x: int, bits: int) -> str:
    return f"#x{x & ((1 << bits) - 1):0{(bits + 3)//4}x}"


def ceildiv(a, b):
    return (a + b - 1) // b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("d", type=int)
    ap.add_argument("output")
    args = ap.parse_args()
    d = args.d
    assert 0 <= d < 16
    all1024 = (1 << 1024) - 1
    dmask = 15 << 920
    remaining = (all1024 ^ MASK) & (all1024 ^ dmask)
    pmin = PMASK | (d << 920)
    pmax = pmin | remaining
    qmin = ceildiv(N, pmax)
    qmax = N // pmin
    assert qmin <= qmax
    common = 0
    while common < 1024:
        i = 1023 - common
        if ((qmin >> i) & 1) != ((qmax >> i) & 1):
            break
        common += 1

    lines = [
        "(set-logic QF_BV)",
        "(set-option :produce-models true)",
        "(declare-fun p () (_ BitVec 1024))",
        "(declare-fun q () (_ BitVec 1024))",
        f"(assert (= (bvand p {bvhex(MASK,1024)}) {bvhex(PMASK,1024)}))",
        f"(assert (= ((_ extract 923 920) p) {bvhex(d,4)}))",
        "(assert (= ((_ extract 1023 1023) p) #b1))",
        "(assert (= ((_ extract 1023 1023) q) #b1))",
        "(assert (= ((_ extract 0 0) p) #b1))",
        "(assert (= ((_ extract 0 0) q) #b1))",
        f"(assert (bvuge q {bvhex(qmin,1024)}))",
        f"(assert (bvule q {bvhex(qmax,1024)}))",
    ]
    if common:
        lo = 1024 - common
        prefix = qmin >> lo
        lines.append(f"(assert (= ((_ extract 1023 {lo}) q) {bvhex(prefix,common)}))")

    mod = 1 << 265
    nlow = N & (mod - 1)
    for a in range(16):
        plow = (PMASK | (a << 150)) & (mod - 1)
        qlow = (nlow * pow(plow, -1, mod)) & (mod - 1)
        lines.append(
            f"(assert (=> (= ((_ extract 153 150) p) {bvhex(a,4)}) "
            f"(= ((_ extract 264 0) q) {bvhex(qlow,265)})))")

    lines += [
        f"(assert (= (bvmul ((_ zero_extend 1024) p) ((_ zero_extend 1024) q)) {bvhex(N,2048)}))",
        "(check-sat)",
        "(get-value (p q))",
        "(exit)",
    ]
    with open(args.output, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"d={d} common_high={common} qmin_bits={qmin.bit_length()} qmax_bits={qmax.bit_length()}")

if __name__ == "__main__":
    main()
