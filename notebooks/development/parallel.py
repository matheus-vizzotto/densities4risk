"""
Parallel utilities with timing support.
"""

from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
from typing import Callable, Iterable, List, Any, Optional, Tuple
import time


def _get_context():
    return mp.get_context("spawn")


def _timed_worker(args: Tuple[Callable, Any]):
    """
    Internal wrapper to measure execution time of each task.
    """
    func, x = args

    start = time.perf_counter()
    result = func(x)
    end = time.perf_counter()

    return {
        "input": x,
        "result": result,
        "time": end - start
    }


def run_parallel(
    func: Callable,
    data: Iterable,
    n_cores: Optional[int] = None,
    chunksize: int = 1,
    return_task_times: bool = False,
):
    """
    Run a function in parallel and measure execution time.

    Parameters
    ----------
    func : callable
        Function to apply
    data : iterable
        Inputs
    n_cores : int or None
        Number of cores
    chunksize : int
        Chunk size for batching
    return_task_times : bool
        If True, returns per-task timing

    Returns
    -------
    dict with:
        - results
        - total_time
        - task_times (optional)
    """
    if n_cores is None:
        n_cores = mp.cpu_count()

    ctx = _get_context()

    start_total = time.perf_counter()

    with ProcessPoolExecutor(
        max_workers=n_cores,
        mp_context=ctx
    ) as executor:

        if return_task_times:
            wrapped_data = [(func, x) for x in data]
            raw_results = list(
                executor.map(_timed_worker, wrapped_data, chunksize=chunksize)
            )

            results = [r["result"] for r in raw_results]
            task_times = [r["time"] for r in raw_results]

        else:
            results = list(
                executor.map(func, data, chunksize=chunksize)
            )
            task_times = None

    end_total = time.perf_counter()

    output = {
        "results": results,
        "total_time": end_total - start_total,
    }

    if return_task_times:
        output["task_times"] = task_times

    return output


# ---------------------------------------------------------------------
# Sequential version (for benchmarking)
# ---------------------------------------------------------------------

def run_sequential(func: Callable, data: Iterable):
    """
    Run sequentially with timing (baseline comparison).
    """
    start = time.perf_counter()

    results = []
    task_times = []

    for x in data:
        t0 = time.perf_counter()
        results.append(func(x))
        t1 = time.perf_counter()
        task_times.append(t1 - t0)

    end = time.perf_counter()

    return {
        "results": results,
        "total_time": end - start,
        "task_times": task_times
    }


# ---------------------------------------------------------------------
# Example worker
# ---------------------------------------------------------------------

def worker(x):
    time.sleep(0.01)  # simulate 10ms work
    return x * x


if __name__ == "__main__":
    data = list(range(1000))

    seq = run_sequential(worker, data)
    par = run_parallel(worker, data, return_task_times=True)

    print("Sequential time:", seq["total_time"])
    print("Parallel time:", par["total_time"])