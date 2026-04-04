import numpy as np
import scipy
import pandas as pd

import numpy as np
import scipy.stats as stats

def generate_base_density(grid, kind="gaussian", **params):
    """
    Generates a normalized density f(x) evaluated on a grid using 
    flexible distribution parameters.

    Parameters
    ----------
    grid : ndarray
        The support points where the density is evaluated.
    kind : str
        The distribution type: 'gaussian', 'student_t', or 'mixture'.
    **params : dict
        Arbitrary keyword arguments passed to the scipy.stats distribution.
        Examples:
            - gaussian: mu=0, sigma=1 (mapped to loc/scale)
            - student_t: df=5, loc=0, scale=1
            - mixture: list of dicts under 'components'

    Examples
    --------
    >>> # 1. Standard Gaussian
    >>> f_norm = generate_base_density(grid, kind="gaussian", mu=0, sigma=1)
    
    >>> # 2. Heavy-Tailed Student-t 
    >>> f_heavy = generate_base_density(grid, kind="student_t", df=2.5, loc=0, scale=1)
    
    >>> # 4. Bimodal Mixture 
    >>> mix_params = {
    ...     "components": [
    ...         {"type": "norm", "weight": 0.4, "loc": -2, "scale": 0.8},
    ...         {"type": "norm", "weight": 0.6, "loc": 2, "scale": 1.2}
    ...     ]
    ... }
    >>> f_mix = generate_base_density(grid, kind="mixture", **mix_params)
    """
    if kind == "gaussian":
        # Map user-friendly names to scipy standard loc/scale
        loc = params.pop("mu", params.get("loc", 0))
        scale = params.pop("sigma", params.get("scale", 1))
        f = stats.norm.pdf(grid, loc=loc, scale=scale, **params)
        
    elif kind == "student_t":
        # Allows df, loc, scale, and even 'moments' or other scipy args
        f = stats.t.pdf(grid, **params)
        
    elif kind == "mixture":
        # Expects a list of dicts: [{'type': 'norm', 'weight': 0.5, 'mu': 0...}, ...]
        components = params.get("components", [
            {"type": "norm", "weight": 0.5, "loc": -2, "scale": 1},
            {"type": "norm", "weight": 0.5, "loc": 2, "scale": 1}
        ])
        
        f = np.zeros_like(grid)
        for comp in components:
            w = comp.pop("weight", 1/len(components))
            dist_type = comp.pop("type", "norm")
            # Dynamically call the distribution from scipy.stats
            dist_func = getattr(stats, dist_type)
            f += w * dist_func.pdf(grid, **comp)
            
    else:
        # Fallback: try to find the distribution directly in scipy.stats
        try:
            dist_func = getattr(stats, kind)
            f = dist_func.pdf(grid, **params)
        except AttributeError:
            raise ValueError(f"Distribution '{kind}' not recognized in scipy.stats")
    
    # Final normalization check
    total_area = np.trapezoid(f, grid)
    return f / total_area if total_area > 0 else f

def transform_to_l2(f, grid):
    """
    Apply modified Log Quantile Density Transform (mLQDT).
    """
    # Pseudo-steps (adapt to your implementation)
    F = np.cumsum(f) / np.sum(f)
    Q = np.interp(grid, F, grid)  # quantile function
    
    lqd = np.log(np.gradient(Q, grid))
    
    return lqd

def simulate_ar_scores(n, phis, sigma):
    """
    Simulates AR(1) score processes W_{kt}.
    
    Returns:
        scores: (K, n)
    """
    K = len(phis)
    scores = np.zeros((K, n))
    
    for k in range(K):
        shocks = np.random.normal(scale=sigma, size=n)
        for t in range(1, n):
            scores[k, t] = phis[k] * scores[k, t-1] + shocks[t]
    
    return scores

def build_basis(K, u, basis_type="cosine"):
    """
    Constructs an orthonormal trigonometric basis in the L2[0,1] Hilbert space.

    These basis functions are used to represent the stochastic components of 
    the Functional Time Series (FTS), such as the latent process X_t(u) or 
    the error term epsilon_t(u).

    Parameters
    ----------
    K : int
        Number of basis functions to generate (the dimension of the subspace).
    u : ndarray
        Discretization grid in the unit interval [0, 1]. Shape (M,).
    basis_type : {'cosine', 'sine'}, default='cosine'
        The type of trigonometric functions to use. 
        - 'cosine': {sqrt(2) * cos(j * pi * u)} for j = 1, ..., K
        - 'sine': {sqrt(2) * sin(j * pi * u)} for j = 1, ..., K

    Returns
    -------
    basis : ndarray
        An array of shape (K, M) containing the evaluated basis functions.

    Notes
    -----
    The basis is normalized such that the L2 norm of each function is 1, 
    assuming a uniform grid on [0,1]. Specifically:
    <phi_j, phi_l> = integral_{0}^{1} phi_j(u)phi_l(u)du = delta_{jl}
    """
    u = np.asarray(u)
    M = len(u)
    basis = np.zeros((K, M))

    for k in range(K):
        j = k + 1  # Standard indexing starting from 1
        if basis_type == "cosine":
            # Orthonormal Cosine Basis
            basis[k, :] = np.sqrt(2) * np.cos(j * np.pi * u)
        elif basis_type == "sine":
            # Orthonormal Sine Basis
            basis[k, :] = np.sqrt(2) * np.sin(j * np.pi * u)
        else:
            raise ValueError("basis_type must be 'sine' or 'cosine'")

    return basis

def build_dependence(scores, basis):
    """
    Computes X_t(u) = sum_k W_{kt} * phi_k(u)
    
    Returns:
        X: (nt, n)
    """
    return (scores.T @ basis).T

def simulate_brownian_bridge(n, nt, u, sigma_bb=0.05):
    """
    Simulates Brownian bridge noise ε_t(u).
    
    Returns:
        BB: (nt, n)
    """
    BB = np.zeros((nt, n))
    dt = 1 / (nt - 1)

    for t in range(n):
        dw = np.random.normal(scale=np.sqrt(dt), size=nt)
        w = np.cumsum(dw)
        bb = w - u * w[-1]
        BB[:, t] = sigma_bb * bb   # ← key change

    return BB

def simulate_white_noise(n, nt, sigma):
    """
    Simple iid functional noise.
    """
    return np.random.normal(scale=sigma, size=(nt, n))

def assemble_l2_process(base_lqd, X, noise):
    """
    Combines all components:
    
    Y_t(u) = base_lqd + X_t(u) + ε_t(u)
    """
    if base_lqd.ndim == 1:
        base_lqd = base_lqd[:, None]
        
    return base_lqd + X + noise

def simulate_l2_process(
    n,
    base_lqd,
    u,
    phis,
    sigma=0.05,
    noise_type="bridge",
    basis="cosine",
    as_dataframe=False,
    date_init=None,
    date_end=None
):
    """
    Simulate a functional time series in L2 via basis-driven dynamics.

    Generates curves Y_t(u) = base_lqd(u) + X_t(u) + ε_t(u), where X_t(u)
    is constructed from AR(1) score processes and an orthonormal cosine basis,
    and ε_t(u) is either Brownian bridge or white noise.

    Parameters
    ----------
    n : int
        Number of curves (time points).
    base_lqd : array-like, shape (nt,)
        Mean function in L2 space.
    u : array-like, shape (nt,)
        Grid over [0, 1].
    phis : array-like, shape (K,)
        AR(1) coefficients.
    sigma : float, optional
        Innovation standard deviation.
    noise_type : {"bridge", "white"}, optional
        Type of functional noise.
    as_dataframe : bool, optional
        If True, return Y as a pandas DataFrame with a daily date index.
    date_init : str or datetime-like, optional
        Start date for the index.
    date_end : str or datetime-like, optional
        End date for the index.

    Returns
    -------
    Y : ndarray or DataFrame, shape (nt, n)
        Simulated functional data.
    scores : ndarray, shape (K, n)
    X : ndarray, shape (nt, n)
    noise : ndarray, shape (nt, n)

    Notes
    -----
    If `as_dataframe=True` and no dates are provided, columns are set to the
    last `n` calendar days ending today.

    Example
    -------
    >>> Y, *_ = simulate_l2_process(
    ...     n=30,
    ...     base_lqd=base_lqd,
    ...     u=u,
    ...     phis=[0.8, 0.5],
    ...     as_dataframe=True
    ... )
    """

    nt = len(u)
    K = len(phis)

    scores = simulate_ar_scores(n, phis, sigma)
    basis = build_basis(K, u, basis_type=basis)
    X = build_dependence(scores, basis)

    if noise_type == "bridge":
        noise = simulate_brownian_bridge(n, nt, u)
    else:
        noise = simulate_white_noise(n, nt, sigma)

    Y = assemble_l2_process(base_lqd, X, noise)

    if as_dataframe:
        # --- Date logic ---
        if date_init is None and date_end is None:
            dates = pd.date_range(end=pd.Timestamp.today().normalize(), periods=n, freq="D")
        elif date_end is not None and date_init is None:
            dates = pd.date_range(end=date_end, periods=n, freq="D")
        elif date_init is not None and date_end is None:
            dates = pd.date_range(start=date_init, periods=n, freq="D")
        else:
            dates = pd.date_range(start=date_init, end=date_end, freq="D")
            if len(dates) != n:
                raise ValueError("date range length must match n")

        Y = pd.DataFrame(Y, index=u, columns=dates)

    return Y, scores, X, noise