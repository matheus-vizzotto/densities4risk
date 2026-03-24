import numpy as np
from scipy.stats import norm, t 
from scipy.integrate import cumulative_trapezoid

#---------------------------------------------
#------------- L2 curves ---------------------
#---------------------------------------------

def generate_coeffs(d):
    c = []
    for l0 in range(d):
        l = l0 + 1
        c0 = (-1)**l * (0.9 - ((0.5*l)/d))
        # print(c0, "\n")
        c.append(c0)
    return c

def simulate_curves(
        n: int,           # number of curves
        nt: int,          # grid size
        u: np.ndarray,    # nt x 1 grid over [0,1]
        phis,             # list of AR parameters for xi_k
        variances=None   # optional list of variances
    ):
    """_summary_

    Args:
        n (int): _description_

    Returns:
        _type_: _description_
    
    Example:
    n = 100          # sample size (curves)
    d = 2            # dimension parameter
    nt = 256         # number of grid points
    u = np.linspace(0.01, 0.99, nt)[:, None]  # nt x 1 grid
    phis = [-0.775, 0.65, -0.525, 0.4]
    variance = 0.01
    Y, X, mEps = simulate_curves(n,nt,u,phis, variances=np.full(len(phis), variance))
    """

    phis = np.asarray(phis)
    K = len(phis)

    # variances
    if variances is None:
        variances = np.ones(K)

    # --- simulate latent scores xi_k(t)
    Xi = np.zeros((K, n))
    for k in range(K):
        # eps = np.random.normal(scale=np.sqrt(variances[k]), size=n)
        eps = np.random.normal(scale=variances[k], size=n)
        for t in range(1, n):
            Xi[k, t] = phis[k] * Xi[k, t - 1] + eps[t]

    # --- build latent functional signal X
    X = np.zeros((nt, n))
    for k in range(K):
        X += Xi[k][None, :] * np.sqrt(2) * np.cos((k+1) * np.pi * u)

    # ADICIONA RUIDO À SÉRIE DE SENO 
    mEps = np.zeros((nt, n))
    for ii in range(n):
        for jj in range(1, 11): # c.l. de 10 componentes vem do Bathia et al (2010) e Rodney & Pinheiro (2020) 
            mEps[:, ii] += (
                norm.rvs(scale=1.0) * np.sqrt(2) *
                np.sin(np.pi * u[:, 0] * jj) / (2 ** (jj - 1))
            )

    # DADOS FUNCIONAIS OBSERVADOS: Y(t) = X(t) + epsilon(t)
    Y = X + mEps

    return Y, X, mEps

#---------------------------------------------
#------------- Densities ---------------------
#---------------------------------------------
def normal_density(x, loc=0, scale=1, normalize=True):
    """
    Generate a Gaussian (normal) density on grid x.

    Parameters
    ----------
    x : array-like
        Support grid.
    loc : float
        Mean of the distribution.
    scale : float
        Standard deviation.
    normalize : bool
        Whether to normalize numerically.

    Returns
    -------
    pdf : ndarray
        Density evaluated on x.

    Example
    -------
    >>> x = np.linspace(-10, 10, 2000)
    >>> pdf = normal_density(x, loc=0, scale=1)
    """
    pdf = norm.pdf(x, loc=loc, scale=scale)
    
    if normalize:
        pdf = pdf / np.trapezoid(pdf, x)
        
    return pdf

def t_density(x, df=5, loc=0, scale=1, normalize=True):
    """
    Generate a Student-t density on grid x.

    Parameters
    ----------
    x : array-like
        Support grid.
    df : float
        Degrees of freedom.
    loc : float
        Location parameter.
    scale : float
        Scale parameter.
    normalize : bool
        Whether to normalize numerically.

    Returns
    -------
    pdf : ndarray
        Density evaluated on x.

    Example
    -------
    >>> x = np.linspace(-10, 10, 2000)
    >>> pdf = t_density(x, df=3)  # heavy tails
    """
    pdf = t.pdf(x, df=df, loc=loc, scale=scale)
    
    if normalize:
        pdf = pdf / np.trapezoid(pdf, x)
        
    return pdf

def gaussian_mixture_density(x, weights, locs, scales, normalize=True):
    """
    Generate a Gaussian mixture density on grid x.

    Parameters
    ----------
    x : array-like
        Support grid.
    weights : list of float
        Mixture weights (must sum to 1).
    locs : list of float
        Means of each component.
    scales : list of float
        Standard deviations of each component.
    normalize : bool
        Whether to normalize numerically.

    Returns
    -------
    pdf : ndarray
        Density evaluated on x.

    Example
    -------
    >>> x = np.linspace(-10, 10, 2000)
    >>> pdf = gaussian_mixture_density(
    ...     x,
    ...     weights=[0.9, 0.1],
    ...     locs=[0, 0],
    ...     scales=[1, 5]
    ... )
    """
    
    if not (len(weights) == len(locs) == len(scales)):
        raise ValueError("weights, locs, and scales must have same length")
    
    if not np.isclose(sum(weights), 1):
        raise ValueError("weights must sum to 1")
    
    pdf = np.zeros_like(x, dtype=float)
    
    for w, mu, sigma in zip(weights, locs, scales):
        pdf += w * norm.pdf(x, loc=mu, scale=sigma)
    
    if normalize:
        pdf = pdf / np.trapezoid(pdf, x)
        
    return pdf

def t_mixture_density(x, weights, dfs, locs, scales, normalize=True):
    """
    Generate a Student-t mixture density on grid x.

    Parameters
    ----------
    x : array-like
        Support grid.
    weights : list of float
        Mixture weights (must sum to 1).
    dfs : list of float
        Degrees of freedom for each component.
    locs : list of float
        Location parameters.
    scales : list of float
        Scale parameters.
    normalize : bool
        Whether to normalize numerically.

    Returns
    -------
    pdf : ndarray
        Density evaluated on x.

    Example
    -------
    >>> x = np.linspace(-10, 10, 2000)
    >>> pdf = t_mixture_density(
    ...     x,
    ...     weights=[0.8, 0.2],
    ...     dfs=[10, 3],
    ...     locs=[0, 0],
    ...     scales=[1, 2]
    ... )
    """

    if not (len(weights) == len(dfs) == len(locs) == len(scales)):
        raise ValueError("All parameter lists must have same length")
    
    if not np.isclose(sum(weights), 1):
        raise ValueError("weights must sum to 1")
    
    pdf = np.zeros_like(x, dtype=float)
    
    for w, df_, mu, sigma in zip(weights, dfs, locs, scales):
        pdf += w * t.pdf(x, df=df_, loc=mu, scale=sigma)
    
    if normalize:
        pdf = pdf / np.trapezoid(pdf, x)
        
    return pdf

def simulate_dependent_densities(n_days, lqd_template, phis, sigma=0.1):
    """
    Simulates a sequence of LQD curves with temporal dependence.
    
    Args:
        n_days (int): Number of days to simulate.
        lqd_template (array): Your smooth T-distribution LQD (shape: nt,).
        phis (list): AR(1) coefficients for each 'flavor' of movement.
        sigma (float): Volatility of the daily changes.
    """
    nt = len(lqd_template)
    K = len(phis)
    u = np.linspace(0, 1, nt)
    
    # 1. Initialize the 'Scores' (How much of each cosine we use each day)
    # Shape: (K, n_days)
    scores = np.zeros((K, n_days))
    for k in range(K):
        # Daily shocks (Innovation)
        shocks = np.random.normal(scale=sigma, size=n_days)
        for t in range(1, n_days):
            scores[k, t] = phis[k] * scores[k, t-1] + shocks[t]
            
    # 2. Define the 'Shapes' (The basis functions)
    # We use sqrt(2)*cos to keep things orthonormal in L2
    basis = np.array([np.sqrt(2) * np.cos(np.pi * (k+1) * u) for k in range(K)])
    
    # 3. Combine: Template + (Scores * Basis)
    # We use matrix multiplication (@) for the 'perturbation' part
    perturbations = (scores.T @ basis).T # Result shape: (nt, n_days)
    
    # Add the template to every single day
    simulated_lqds = lqd_template[:, None] + perturbations
    
    return simulated_lqds, scores

def sample_from_density(x, y, n_samples=1000):
    # 1. Compute the numerical CDF using the trapezoidal rule
    # initial=0 ensures the resulting array has the same length as x and y
    cdf = cumulative_trapezoid(y, x, initial=0)
    
    # 2. Normalize to handle any tiny integration residuals 
    # (Ensures the max value is exactly 1.0)
    cdf /= cdf[-1]
    
    # 3. Generate uniform random numbers
    u = np.random.uniform(0, 1, n_samples)
    
    # 4. Inverse transform: map U back to the X support
    return np.interp(u, cdf, x)