"""
Parallel utilities for running CPU-bound independent tasks using all cores.

Designed to work safely with:
- macOS (Apple Silicon, spawn mode)
- Jupyter notebooks (when imported)

Usage (from notebook):

    from parallel_utils import run_parallel

    def my_task(x):
        return x * x

    data = list(range(1000))
    results = run_parallel(my_task, data)

Notes
-----
- The worker function MUST be defined at top level (not nested).
- Avoid lambdas.
- Inputs must be pickleable.
"""

from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
from typing import Callable, Iterable, List, Any


def _get_context():
    """
    Ensures safe multiprocessing context on macOS.
    """
    return mp.get_context("spawn")


def run_parallel(
    func: Callable,
    data: Iterable,
    n_cores: int | None = None,
    chunksize: int = 1,
) -> List[Any]:
    """
    Run a function in parallel over an iterable.

    Parameters
    ----------
    func : callable
        Function to apply. Must be defined at module level.
    data : iterable
        Inputs to process.
    n_cores : int or None
        Number of CPU cores to use. If None, uses all available.
    chunksize : int
        Number of tasks sent per worker batch (important for performance).

    Returns
    -------
    list
        Results in the same order as input.
    """
    if n_cores is None:
        n_cores = mp.cpu_count()

    ctx = _get_context()

    with ProcessPoolExecutor(
        max_workers=n_cores,
        mp_context=ctx
    ) as executor:

        results = list(
            executor.map(func, data, chunksize=chunksize)
        )

    return results


# ---------------------------------------------------------------------
# Example worker (optional, for testing)
# ---------------------------------------------------------------------

def example_worker(x: int) -> int:
    """
    Simple test function.

    Example
    -------
    >>> from parallel_utils import run_parallel, example_worker
    >>> run_parallel(example_worker, [1,2,3])
    [1, 4, 9]
    """
    return x * x


# ---------------------------------------------------------------------
# Script mode (optional)
# ---------------------------------------------------------------------

if __name__ == "__main__":
    # This block allows running the file directly for testing
    data = list(range(10))
    results = run_parallel(example_worker, data)
    print(results)