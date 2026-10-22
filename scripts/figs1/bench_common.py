"""Timing helpers shared by the Figure S1 benchmarks."""

import time

import numpy as np

SLOW = 30  # seconds: above this, a single repeat is used
GIVE_UP = 900  # seconds: stop increasing the tree size once one measurement takes longer than this


def timeit(fn, repeats=3):
    """Median wall time (s) of `fn`; one repeat if the first call is slow."""
    start = time.perf_counter()
    fn()
    first = time.perf_counter() - start
    if first > SLOW:
        return first
    times = [first]
    for _ in range(repeats - 1):
        start = time.perf_counter()
        fn()
        times.append(time.perf_counter() - start)
    return float(np.median(times))
