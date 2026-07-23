#!/usr/bin/env python3
"""Enhanced root extraction for the centered Herrmann-May lattice.
Adds short pairwise sums/differences of reduced basis vectors before polynomial
reconstruction, then rejects wrong candidates with a finite-field check before
running numerical Newton."""
from pathlib import Path

import hm_flatter_solver as hm

_original_reconstruct = hm.reconstruct


def enhanced_reconstruct(rows, mons, xp, yp, p_lower, t, beta, max_polys=18):
    ranked = sorted(rows, key=lambda r: sum(z*z for z in r))
    seed = ranked[:16]
    extra = []
    for i in range(len(seed)):
        for j in range(i):
            a, b = seed[i], seed[j]
            extra.append([x+y for x,y in zip(a,b)])
            extra.append([x-y for x,y in zip(a,b)])
    return _original_reconstruct(ranked + extra, mons, xp, yp, p_lower, t, beta,
                                 max_polys=max_polys)


hm.reconstruct = enhanced_reconstruct

import hm_newton_solver as solver

_original_try_newton = solver.try_newton
_first_prime = hm.primes_30(1)[0]


def filtered_try_newton(polys, constant, log_path: Path):
    rr = hm.groebner_residue(polys, _first_prime, singular="Singular", timeout=180)
    with log_path.open("a") as log:
        log.write(f"PREFILTER prime={_first_prime} residue={rr}\n")
    if rr is False:
        return None
    return _original_try_newton(polys, constant, log_path)


solver.try_newton = filtered_try_newton

if __name__ == "__main__":
    raise SystemExit(solver.main())
