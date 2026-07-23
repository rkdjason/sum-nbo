#!/usr/bin/env python3
import argparse
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


def ceildiv(a,b): return (a+b-1)//b


def parse_meta(path):
    data={}
    for line in Path(path).read_text().splitlines():
        if line.strip():
            k,*v=line.split(); data[k]=v
    return data


def lit_for_value(entry, expected):
    """Literal that is true exactly when the represented bit equals expected.
    Return True/False for constants, otherwise a signed DIMACS literal."""
    x=int(entry)
    if x == -1: return bool(expected)
    if x == 0: return not bool(expected)
    return x if expected else -x


def mismatch_lit(entry, expected):
    # True when bit differs from expected.
    v=lit_for_value(entry, expected)
    return (not v) if isinstance(v,bool) else -v


def simplify_clause(items):
    out=[]; seen=set()
    for x in items:
        if x is True: return None
        if x is False: continue
        x=int(x)
        if -x in seen: return None
        if x not in seen: seen.add(x); out.append(x)
    return out


def implication_clause(cond_entries, conclusion):
    # cond_entries are (entry, expected), conclusion is (entry, expected).
    items=[mismatch_lit(e,v) for e,v in cond_entries]
    items.append(lit_for_value(*conclusion))
    return simplify_clause(items)


def common_prefix(a,b,width=1024):
    n=0
    for i in range(width-1,-1,-1):
        if ((a>>i)&1) != ((b>>i)&1): break
        n+=1
    return n


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('meta')
    ap.add_argument('output')
    ap.add_argument('--low-depth',type=int,default=14)
    ap.add_argument('--high-depth',type=int,default=14)
    args=ap.parse_args()
    md=parse_meta(args.meta)
    d=int(md['D'][0]); pmap=md['P']; qmap=md['Q']; sels=md['SEL']
    clauses=[]

    # Low Hensel truth-table hints: selected a plus r low bits of the 84-bit block
    # determine q modulo 2^(265+r), so emit the newly determined q bit.
    for a in range(16):
        sel=sels[a]
        for r in range(1,args.low_depth+1):
            k=265+r
            mod=1<<k
            for pref in range(1<<r):
                plow=(PMASK | (a<<150) | (pref<<265)) & (mod-1)
                qlow=(N*pow(plow,-1,mod)) & (mod-1)
                expected=(qlow>>(k-1))&1
                cond=[(sel,1)]
                cond += [(pmap[265+j], (pref>>j)&1) for j in range(r)]
                c=implication_clause(cond,(qmap[k-1],expected))
                if c is not None: clauses.append(c)

    # High reciprocal interval hints: a prefix of the top unknown 46-bit block
    # fixes additional leading q bits. Emit every newly common q bit.
    all1024=(1<<1024)-1
    dmask=15<<920
    base_unknown=(all1024^MASK) & (all1024^dmask)
    pbase=PMASK | (d<<920)
    pmin0=pbase
    pmax0=pbase | base_unknown
    base_common=common_prefix(ceildiv(N,pmax0),N//pmin0)
    assigned_mask=0
    for r in range(1,args.high_depth+1):
        pos=830-r
        assigned_mask |= 1<<pos
        remaining=base_unknown & (all1024^assigned_mask)
        for pref in range(1<<r):
            add=0
            cond=[]
            for j in range(r):
                bitpos=829-j
                val=(pref>>(r-1-j))&1
                if val: add |= 1<<bitpos
                cond.append((pmap[bitpos],val))
            pmin=pbase|add
            pmax=pmin|remaining
            qmin,qmax=ceildiv(N,pmax),N//pmin
            com=common_prefix(qmin,qmax)
            for depth in range(base_common+1,com+1):
                qpos=1024-depth
                expected=(qmin>>qpos)&1
                c=implication_clause(cond,(qmap[qpos],expected))
                if c is not None: clauses.append(c)

    with open(args.output,'w') as f:
        for c in clauses: f.write(' '.join(map(str,c))+' 0\n')
    print(len(clauses))

if __name__=='__main__': main()
