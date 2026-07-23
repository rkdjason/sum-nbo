#!/usr/bin/env python3
"""Run hm_newton_solver, but gate expensive numerical Newton iterations behind
one finite-field consistency check. Wrong outer-bit candidates normally produce
the unit ideal and are rejected immediately."""
import sys
from pathlib import Path

import hm_flatter_solver as hm
import hm_newton_solver as solver

_original_try_newton = solver.try_newton
_first_prime = hm.primes_30(1)[0]


def filtered_try_newton(polys, constant, log_path: Path):
    rr = hm.groebner_residue(polys, _first_prime, singular="Singular", timeout=180)
    with log_path.open("a") as log:
        log.write(f"PREFILTER prime={_first_prime} residue={rr}\n")
    if rr is False:
        return None
    # If Singular timed out or did not isolate the root, retain Newton as a
    # robust independent extraction method. A concrete residue also proceeds.
    return _original_try_newton(polys, constant, log_path)


solver.try_newton = filtered_try_newton

if __name__ == "__main__":
    raise SystemExit(solver.main())
