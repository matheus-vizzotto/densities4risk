import numpy as np
from scipy.stats import norm

def simulate_curves(
        n: int,           # number of time points
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
        eps = np.random.normal(scale=np.sqrt(variances[k]), size=n)
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