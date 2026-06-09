import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pywt
from typing import Tuple



################################ NOVA CLASSE #####################

# TODO: Unificar notação e termos. Exemplo: Yhat ou fitted_values?


# ============================================================
# Auxiliary Functions
# ============================================================

def inner_product(f, g, du):
    """
    Compute the L^2 inner product ⟨f, g⟩ = ∫ f(u) g(u) du
    in its discrete approximation.

    Parameters
    ----------
    f : ndarray (m, k1)
        First function(s) evaluated on a grid.
    g : ndarray (m, k2)
        Second function(s) evaluated on the same grid.
    du : float
        Grid spacing.

    Returns
    -------
    ndarray (k1, k2)
        Matrix of inner products.
    """
    return np.dot(f.T, g) * du


def L2norm(f, du):
    """
    Compute the L^2 norm of a function vector f.

    Parameters
    ----------
    f : ndarray (m,)
        Function values on a grid.
    du : float
        Grid spacing.

    Returns
    -------
    float
        The L^2 norm.
    """
    return np.sqrt(inner_product(f, f, du))


# ============================================================
# Main Class
# ============================================================

import numpy as np

def _bootstrap_select_n_components(
    Y, Ydev, Ybar, thetahat_old, gammahat_old, 
    n, p, m, du, B=100, lag_max=10, alpha=0.05
    ):
    """
    Selects the number of components d0 using a bootstrap testing procedure.
    
    Parameters:
    -----------
    Y : np.ndarray (m, n)
        Original data matrix.
    Ydev : np.ndarray (m, n)
        Mean-centered data.
    Ybar : np.ndarray (m, 1)
        Mean vector.
    thetahat_old : list/array
        Initial eigenvalues from the original data.
    gammahat_old : np.ndarray (n-p, lag_max)
        Eigenvectors/basis from the original K* matrix.
    du : float
        Step size for the inner product (integration).
    B : int
        Number of bootstrap iterations.
    lag_max : int
        Maximum number of components to test.
    alpha : float
        Significance level for the test.
    """
    bs_pvalues = np.zeros(lag_max)

    for d0 in range(1, lag_max + 1):
        thetahatH0 = thetahat_old[d0]
        gammahat = gammahat_old[:, :d0]

        # Projection and normalization
        psihat_root = np.dot(Ydev[:, :(n - p)], gammahat)
        psihat = np.zeros((m, d0))
        for i in range(d0):
            norm_val = L2norm(psihat_root[:, i], du)
            psihat[:, i] = psihat_root[:, i] / norm_val if norm_val > 1e-15 else 0

        # etahat calculation (assuming inner_product handles (m, d0) and (m, n))
        # Replacing with explicit matrix mult if inner_product is standard
        etahat = np.dot(psihat.T, Ydev) * du
        Yhat = Ybar + np.dot(psihat, etahat)
        epsilonhat = Y - Yhat

        bs_thetahat = np.zeros(B)
        
        for b in range(B):
            # Bootstrap residuals
            bs_epsilon = epsilonhat[:, np.random.permutation(n)]
            bs_Y = Yhat + bs_epsilon

            bs_Ybar = np.mean(bs_Y, axis=1, keepdims=True)
            bs_Ydev = bs_Y - bs_Ybar

            # Core inner product matrix
            bs_core = np.dot(bs_Ydev.T, bs_Ydev) * du
            bs_Kstar_core0 = bs_core[:(n - p), :(n - p)]

            # Vectorized construction of K* core sum
            bs_Kstar_sum = np.zeros((n - p, n - p))
            for k in range(1, p + 1):
                # Slicing the precomputed core
                bs_Kstar_sum += bs_core[k:(n - (p - k)), k:(n - (p - k))]

            bs_Kstar = (1.0 / (n - p)**2) * np.dot(bs_Kstar_sum, bs_Kstar_core0)

            # Eigenvalue computation
            bs_eigs = np.linalg.eigvals(bs_Kstar)
            # Sort descending and take the (d0)-th eigenvalue
            bs_thetahat[b] = np.sort(np.real(bs_eigs))[::-1][d0]

        bs_pvalues[d0 - 1] = np.mean(bs_thetahat >= thetahatH0)

    # Component selection logic
    d0_candidates = np.where(bs_pvalues < alpha)[0]
    final_d0 = d0_candidates[0] + 1 if len(d0_candidates) > 0 else 1

    return final_d0, bs_pvalues

class K_dFPC:
    """
    Implements the estimation of a KLE-based dynamic factor model
    for functional time series Y(u, t).

    This class takes as input a matrix Y of shape (m, T), where:
        m = number of grid points
        T = number of time points (curves)

    After calling `.fit(...)`, the class computes and stores:
        - Ybar          : mean curve
        - psihat        : estimated spatial basis functions
        - etahat        : temporal score series
        - thetahat      : eigenvalues
        - gammahat      : eigenvectors of the K* operator
        - Yhat          : fitted reconstruction
        - fitted_values : same as Yhat
        - epsilonhat    : residual curves
        - d0            : selected dimension

    Methods
    -------
    fit(...)
        Fit the dynamic KLE model.

    plot_psihat()
        Plot the spatial basis functions.

    plot_etahat()
        Plot the temporal scores.

    predict(etahat_values)
        Reconstruct curves from arbitrary scores.
    
    Example
    -------
    model = K_dFPC(Y.values)
    model.fit(lag_max=6, B=1000, alpha=0.10, du=0.05, p=5, m=Y.shape[0],
                u=Y.index.values, select_ncomp=False,dimension=10
    )
    model.plot_psihat()
    model.plot_etahat()
    scores_pred = model.forecast_scores(h=1)
    predicted_curves = model.predict(scores_pred)
    """

    # ------------------------------------------------------------

    def __init__(self, Y):
        """
        Initialize the model with data.

        Parameters
        ----------
        Y : ndarray (m, T)
            Functional time series sample evaluated on an m-point grid.
        """
        self.Y = Y
        self.m, self.T = Y.shape

        # Filled after fit()
        self.Ybar = None
        self.thetahat = None
        self.gammahat = None
        self.psihat = None
        self.etahat = None
        self.Yhat = None
        self.epsilonhat = None
        self.d0 = None
        self.u = None
        self.bs_pvalues = None
        self.fitted_values = None

    # ------------------------------------------------------------
    def _compute_kstar_eigenvalues(self, Y_mat, p, du, return_vecs=False):
            """
            Unified helper to compute the eigenvalues (and optionally eigenvectors)
            of the K* matrix for any given data matrix.
            """
            n = Y_mat.shape[1]
            
            # 1. Centering
            Ybar = np.mean(Y_mat, axis=1, keepdims=True)
            Ydev = Y_mat - Ybar
            
            # 2. Core inner product matrix (m, n) -> (n, n)
            # Using np.dot directly is faster than calling a wrapper inside a loop
            core = np.dot(Ydev.T, Ydev) * du
            
            # 3. Construct K* components
            Kstar_core0 = core[:(n - p), :(n - p)]
            Kstar_sum = np.zeros((n - p, n - p))
            for k in range(1, p + 1):
                Kstar_sum += core[k:(n - (p - k)), k:(n - (p - k))]
                
            Kstar = (n - p)**(-2) * np.dot(Kstar_sum, Kstar_core0)
            
            # 4. Solve
            if return_vecs:
                vals, vecs = np.linalg.eig(Kstar)
                idx = np.argsort(np.real(vals))[::-1]

                return np.real(vals)[idx], np.real(vecs)[:, idx], Ydev, Ybar
            else:
                # eigvals is much faster if we don't need vectors (perfect for bootstrap)
                vals = np.linalg.eigvals(Kstar)

                return np.sort(np.real(vals))[::-1]

    def fit(self, p, u, du=0.05,  lag_max=5, B=1_000, alpha=0.05, 
            select_ncomp=False, dimension=None):
        """
        Fit the dynamic KLE model.

        Parameters
        ----------
        lag_max : int
            Maximum number of components to test under bootstrapping.
        
        B : int
            Number of bootstrap replications.

        alpha : float
            Significance level for component selection.

        du : float
            Grid spacing in the domain of the curves.

        p : int
            Time-lag order used to form K*.

        m : int
            Number of grid points (redundant, but kept for compatibility).

        u : ndarray (m,)
            Grid support.

        select_ncomp : bool, optional (default=False)
            Whether to run the bootstrap procedure to determine d0.

        dimension : int, optional
            Fixed number of components to extract (used only if select_ncomp=False).

        Returns
        -------
        self : K_dFPC
            The fitted model.
        """
        self.u, self.m = u, self.Y.shape[0]
        n = self.Y.shape[1]

        # 1. Initial Fit using the unified method
        thetahat_old, gammahat_old, Ydev, Ybar = self._compute_kstar_eigenvalues(
                                                                            self.Y, p, du, return_vecs=True
                                                                        )

        # 2. Component Selection
        if select_ncomp:
            d0, bs_pvalues = _bootstrap_select_n_components(
                                            self.Y, 
                                            Ydev, 
                                            Ybar, 
                                            thetahat_old, 
                                            gammahat_old, 
                                            n, 
                                            p, 
                                            self.m, 
                                            du, 
                                            B, 
                                            lag_max, 
                                            alpha
            )
            self.bs_pvalues = bs_pvalues
        else:
            d0 = dimension if dimension is not None else 1

        self.d0 = d0

        # --------------------------------------------------------
        # Step 4 — Final estimation
        # --------------------------------------------------------

        thetahat = thetahat_old[:d0]
        gammahat = gammahat_old[:, :d0]

        psihat_root = np.dot(Ydev[:, :(n - p)], gammahat)

        psihat = np.zeros((self.m, d0))
        for i in range(d0):
            psihat[:, i] = psihat_root[:, i] / L2norm(psihat_root[:, i], du)

        etahat = inner_product(psihat, Ydev, du)

        Yhat = Ybar + np.dot(psihat, etahat)
        epsilonhat = self.Y - Yhat

        # Store results
        self.Ybar = Ybar
        self.thetahat = thetahat
        self.gammahat = gammahat
        self.psihat = pd.DataFrame(psihat, columns=[f"basis_{i}" for i in range(1,d0+1)])
        self.etahat = pd.DataFrame(etahat.T, columns=[f"scores_{i}" for i in range(1,d0+1)])
        self.Yhat = Yhat
        self.epsilonhat = epsilonhat

        # Prediction storage (Ybar + ψ * η_pred_)
        self.fitted_values = Ybar + psihat @ etahat

        return self

    # ------------------------------------------------------------
    # Plotting Methods
    # ------------------------------------------------------------

    def plot_psihat(self):
        """
        Plot spatial basis functions ψ_i(u).

        Produces a line plot of estimated eigenfunctions on the grid.
        """
        if self.psihat is None:
            raise ValueError("Model not fitted yet.")

        d = self.psihat.shape[1]

        plt.figure()
        for i in range(d):
            plt.plot(self.u, self.psihat[:, i], alpha=0.6, label=f"ψ_{i+1}")

        plt.title("KLE-based decomposition: $\\hat{\\psi}$")
        plt.legend()
        plt.show()

    def plot_etahat(self):
        """
        Plot temporal score series η_i(t).

        Produces a line plot for each score over time.
        """
        if self.etahat is None:
            raise ValueError("Model not fitted yet.")

        plt.figure()
        for i in range(self.etahat.shape[0]):
            plt.plot(self.etahat[i, :], alpha=0.6, label=f"η_{i+1}")

        plt.title("KLE-based temporal scores $\\hat{\\eta}$")
        plt.legend()
        plt.show()

    # ------------------------------------------------------------
    # Prediction Method
    # ------------------------------------------------------------

    def predict(self, etahat_values):
        """
        Reconstruct curves from temporal scores. Can be applied to
        fitted or forecasted scores.

        Compute:
            Ŷ = Ȳ + ψ  η

        Parameters
        ----------
        etahat_values : ndarray (d0, T)
            Scores used to reconstruct the curves.

        Returns
        -------
        ndarray (m, T)
            Reconstructed functional observations.
        """
        if self.Ybar is None or self.psihat is None or self.etahat is None:
            raise ValueError("Model not fitted yet.")
            

        return self.Ybar + self.psihat @ etahat_values

######################## WAVELETS ########################

def wavedec_sizes(signal_length, wavelet_name, level):
    """
    Compute expected pywt.wavedec sizes:
      returns sizes = [len(A_N), len(D_N), ..., len(D_1)]
      and total = sum(sizes)
    """
    w = pywt.Wavelet(wavelet_name)
    F = w.dec_len
    L = signal_length

    lengths = []
    Lj = L
    for j in range(level):
        Lj = (Lj + F - 1) // 2
        lengths.append(Lj)

    sizes = [lengths[-1]] + lengths[::-1]

    total = sum(sizes)
    return sizes, total


# split wavelet vector into coeff lists
def vec2coeffs(vec, sizes):
    coeffs = []
    pos = 0
    for s in sizes:
        coeffs.append(vec[pos:pos+s])
        pos += s
    return coeffs

class W_dFPC:
    """
    Wavelet dynamic Functional Principal Components (W-dFPC)

    Implements the Wavelet Dimensionality Estimation procedure described in
    Fonseca & Pinheiro (2020).

    Given a sample of functional observations Y(t),
    the procedure computes:

        - wavelet decomposition coefficients,
        - mean wavelet coefficients,
        - wavelet-domain covariance operator,
        - eigenvalues and eigenfunctions,
        - scores,
        - reconstructed sample curves.

    Parameters
    ----------
    Y : np.ndarray
        Functional data matrix of shape (nt, n), where
        nt = number of time grid points in each curve,
        n  = number of observed curves.
        Columns of Y are curves Y_t(u).

    Example
    ----------
    n  = 100
    d  = 2
    nt = 256
    u  = np.linspace(0.01, 0.99, nt)[:, None]
    Y (nt X T)
    N = 3
    wavelet = 'db2'
    p = 5
    dimensions = 10
    model = W_dFPC(Y)  
    model.fit(
        nt=nt,
        N=N,
        wavelet=wavelet,
        p=p,
        d=dimensions)
    model.plot_hm(shared=True)
    model.plot_scores(shared=False)
    t = 0
    fitted_scores = model.scores[:, t]
    predicted_Y = model.predict(fitted_scores)
    """

    def __init__(self, Y):
        """
        Construct the estimator.

        Parameters
        ----------
        Y : np.ndarray
            Functional data of shape (nt, n). Each column is one curve.
        """
        self.Y = Y
        self.m, self.T = Y.shape

        #Filled after fit()
        self.Ybar      = None
        self.H         = None
        self.L         = None
        self.B         = None
        self.scores    = None
        self.Yhat      = None

    def fit(self,
        # Y :      np.ndarray,
        nt:      np.ndarray,
        N :      int         = 3,
        wavelet: str         = 'db2',
        p :      int         = 5,
        d :      int         = 3
        ):
        """
        Estimate Wavelet FPComponents.

        Parameters
        ----------
        nt : int
            Number of grid points in each curve (same as Y.shape[0]).
        N : int, optional
            Wavelet decomposition depth (default = 3).
        wavelet : str, optional
            PyWavelets wavelet name, e.g. 'db2'.
        p : int, optional
            Number of lags used to estimate the lagged covariance
            operator in wavelet domain (default = 5).
        d : int, optional
            Number of principal components to retain (default = 3).

        Returns
        -------
        self : W_dFPC
            Instance containing:
                - H eigenfunctions (nt X T)
                - scores (d X T)
                - eigenvalues L (nt X 1)
                - eigenvectors B (nt X nt)
                - reconstructed functional sample Yhat (nt X T)
                - mean function Ybar (nt X 1)
        """

        Y = self.Y
        n = self.T
        sizes, J = wavedec_sizes(nt, wavelet, N)
        
        # ----- "A" matrix -----
        A = np.zeros((J, n))

        for ii in range(n):
            coeffs = pywt.wavedec(Y[:, ii], wavelet=wavelet, level=N)
            A[:, ii] = np.concatenate(coeffs)

        # ----- mean coeffs -----
        mu_A = np.mean(A, axis=1, keepdims=True)
        C = A - mu_A

        # ----- lagged covariance D -----
        C1 = C[:, :n - p]
        D1 = np.zeros((n - p, n - p))

        for k in range(1, p + 1):
            Ct = C[:, k:(n - p + k)]
            D1 += Ct.T @ Ct

        D = C1 @ D1 @ C1.T / ((n - p) ** 2)

        # eigendecomposition (+ ordering)
        L, B = np.linalg.eig(D)
        idx = np.argsort(-L)

        L = L[idx]
        B = B[:, idx]
        
        mu_hat = pywt.waverec(vec2coeffs(mu_A[:, 0], sizes), wavelet)

        # ----- reconstruct eigenfunctions h_m -----
        H = np.zeros((nt, d))
        for m in range(d):
            coeffs_m = vec2coeffs(B[:, m], sizes)
            H[:, m] = pywt.waverec(coeffs_m, wavelet)

        # ----- compute scores -----
        scores = B[:, :d].T @ C      # shape (d, n)

        # ----- reconstruct curves -----
        Yhat = np.zeros_like(Y)
        for t in range(n):
            Yhat[:, t] = mu_hat + H @ scores[:, t]

        # Store results
        self.Ybar      = mu_hat
        self.H         = H
        self.L         = L
        self.B         = B
        self.scores    = scores
        self.Yhat      = Yhat

        return self
    
    def plot_hm(self, shared=True):
        if self.H is None:
            raise ValueError("Model not fitted yet")
        
        d = self.H.shape[1]
        alpha_ = 1  
        increment = 1/d

        if shared:
            plt.figure()
        
        for i in range(d):
            if not shared:
                plt.figure()
                plt.plot(self.H[:, i], label=f"$h_{i}$")
                plt.legend()
            else:
                plt.plot(self.H[:, i], label=f"$h_{i}$", alpha=alpha_)
            alpha_ -= increment

        plt.title("Wavelet-based decomposition: $\hat{h}$")
        plt.legend()
        plt.show()

    def plot_scores(self, shared=True):
        if self.scores is None:
            raise ValueError("Model not fitted yet")

        d = self.scores.shape[0]
        alpha_ = 1  
        increment = 1/d

        if shared:
            plt.figure()
        
        for i in range(d):
            if not shared:
                plt.figure()
                plt.plot(self.scores[i, :], label=f"$\eta_{i}$")
                plt.legend()
            else:
                plt.plot(self.scores[i, :], label=f"$\eta_{i}$", alpha=alpha_)
            alpha_ -= increment

        plt.title("Wavelet-based decomposition: $\eta{h}$")
        plt.legend()
        plt.show()

    def predict(self, new_scores: np.ndarray) -> np.ndarray:
        """
        Reconstruct new functional observations from supplied W-dFPC scores. 
        Useful for reconstructing curves after forecasting scores.

        Parameters
        ----------
        new_scores : np.ndarray
            Matrix of scores of shape (d, k), where:
            d = number of retained components
            k = number of new prediction instants

        Returns
        -------
        Ynew : np.ndarray
            Reconstructed curves of shape (nt, k).
        """
        if new_scores.ndim == 1:
            new_scores = new_scores[:, None]
        # number of new curves
        k = new_scores.shape[1]

        # allocate
        Yhat_predicted = np.zeros((self.Y.shape[0], k))

        for t in range(k):
            Yhat_predicted[:, t] = self.Ybar + self.H @ new_scores[:, t]

        return Yhat_predicted
    
    def eig_proportion(
            self, 
            i:int=2) -> float:
        # proportion of the total eigenvalue mass contained in the first 4 eigenvalues
        prop = np.linalg.norm(self.L[:i], 1) / np.linalg.norm(self.L, 1)
        print(f"Proportion of the total eigenvalue mass contained in the first {i} eigenvalues: {prop:.5%}")
        return prop 
        

########################################
############## COMPARISONS #############
########################################
def super_fun(Y, lag_max, B, alpha, du, p, m, u,
              select_ncomp=False, dimension=None):
    """
    Dynamic KLE estimator (no bootstrap version).

    Parameters
    ----------
    Y : ndarray (m x n)
        Functional observations (each column is a curve).

    lag_max, B, alpha :
        Included for compatibility (not used since bootstrap removed).

    du : float
        Grid spacing in the domain.

    p : int
        Time-lag order.

    m : int
        Number of grid points.

    u : ndarray
        Grid support.

    select_ncomp : bool
        If False, uses fixed dimension.

    dimension : int
        Number of components to extract.

    Returns
    -------
    dict
        Dictionary containing model components.
    """

    # --------------------------------------------------------
    # Step 0 — Basic quantities
    # --------------------------------------------------------

    n = N = Y.shape[1]

    # Mean function
    Ybar = np.mean(Y, axis=1, keepdims=True)

    # Deviations
    Ydev = Y - Ybar

    # --------------------------------------------------------
    # Step 1 — Build core matrices
    # --------------------------------------------------------

    core = inner_product(Ydev, Ydev, du)

    Kstar_core0 = core[:(n - p), :(n - p)]

    Kstar_core = np.zeros((n - p, n - p, p))
    for k in range(1, p + 1):
        Kstar_core[:, :, k - 1] = core[k:(n - (p - k)),
                                       k:(n - (p - k))]

    Kstar_sum = np.zeros((n - p, n - p))
    for k in range(p):
        Kstar_sum += Kstar_core[:, :, k]

    Kstar = (n - p) ** (-2) * (Kstar_sum @ Kstar_core0)

    # --------------------------------------------------------
    # Step 2 — Eigen-decomposition
    # --------------------------------------------------------

    eigvals, eigvecs = np.linalg.eig(Kstar)

    eigvals = np.real(eigvals)
    eigvecs = np.real(eigvecs)

    # Sort descending
    idx = np.argsort(eigvals)[::-1]
    thetahat_old = eigvals[idx]
    gammahat_old = eigvecs[:, idx]

    # --------------------------------------------------------
    # Step 3 — Select number of components
    # --------------------------------------------------------

    if not select_ncomp:
        d0 = dimension
    else:
        raise NotImplementedError("Bootstrap selection removed as requested.")

    # --------------------------------------------------------
    # Step 4 — Final estimation
    # --------------------------------------------------------

    thetahat = thetahat_old[:d0]
    gammahat = gammahat_old[:, :d0]

    # Eigenfunctions (root form)
    psihat_root = Ydev[:, :(n - p)] @ gammahat

    # Normalize
    psihat = np.zeros((m, d0))
    for i in range(d0):
        psihat[:, i] = psihat_root[:, i] / L2norm(psihat_root[:, i], du)

    # Scores
    etahat = inner_product(psihat, Ydev, du)

    # Reconstruction
    Yhat = Ybar + psihat @ etahat

    # --------------------------------------------------------
    # Step 5 — Density correction (exactly as in R)
    # --------------------------------------------------------

    Yhat_fix = Yhat.copy()

    # Enforce positivity
    Yhat_fix[Yhat_fix < 0] = 0

    # Renormalize to integrate to 1
    for t in range(N):
        integral = np.sum(Yhat_fix[:, t]) * du
        if integral > 0:
            Yhat_fix[:, t] /= integral

    # Residuals (same as R: using Yhat, not Yhat_fix)
    epsilonhat = Y - Yhat

    # --------------------------------------------------------
    # Return dictionary (R-style list)
    # --------------------------------------------------------

    return {
        "Y": Y,
        "Ybar": Ybar,
        "thetahat": thetahat,
        "gammahat": gammahat,
        "psihat": psihat,
        "etahat": etahat,
        "Yhat": Yhat,
        "Yhat_fix": Yhat_fix,
        "epsilonhat": epsilonhat,
        "u": u,
        "d0": d0
    }

class K_dFPC_HZ:
    """
    Implements the estimation of a KLE-based dynamic factor model
    for functional time series Y(u, t), with density normalization.

    This class takes as input a matrix Y of shape (m, T), where:
        m = number of grid points
        T = number of time points (curves)

    After calling `.fit(...)`, the class computes and stores:
        - Ybar          : mean curve
        - psihat        : estimated spatial basis functions
        - etahat        : temporal score series
        - thetahat      : eigenvalues
        - gammahat      : eigenvectors of the K* operator
        - Yhat          : fitted reconstruction (raw)
        - Yhat_fix      : fitted reconstruction (normalized/positive)
        - epsilonhat    : residual curves
        - d0            : selected dimension

    Methods
    -------
    fit(...)
        Fit the dynamic KLE model with density correction.

    predict(etahat_values)
        Reconstruct curves from arbitrary scores.
    
    Example
    -------
    model = K_dFPC_Density(Y)
    model.fit(p=5, u=u_grid, du=0.05, dimension=10)
    """

    def __init__(self, Y):
        """
        Initialize the model with data.

        Parameters
        ----------
        Y : ndarray (m, T)
            Functional time series sample evaluated on an m-point grid.
        """
        self.Y = Y
        self.m, self.T = Y.shape

        # Filled after fit()
        self.Ybar = None
        self.thetahat = None
        self.gammahat = None
        self.psihat = None
        self.etahat = None
        self.Yhat = None
        self.Yhat_fix = None
        self.epsilonhat = None
        self.d0 = None
        self.u = None

    # ------------------------------------------------------------

    # def fit(self, p, u, du=0.05, lag_max=5, B=1000, alpha=0.05, 
    #         select_ncomp=False, dimension=None):
    def fit(self, **kwargs):
        """
        Fit the dynamic KLE model by passing arguments to super_fun.
        
        Parameters
        ----------
        **kwargs : 
            Arguments like lag_max, B, alpha, du, p, u, select_ncomp, dimension.
        """
        # Check for bootstrap selection early as per super_fun logic
        if kwargs.get('select_ncomp', False):
            raise NotImplementedError("Bootstrap selection removed as requested.")

        # Call the standalone function using the unpacked dictionary
        # We include m=self.m automatically since it's stored in the class
        sp = super_fun(
            Y=self.Y,
            m=self.m,
            **kwargs
        )

        # Store results directly from the returned dictionary
        self.Ybar        = sp["Ybar"]
        self.thetahat    = sp["thetahat"]
        self.gammahat    = sp["gammahat"]
        self.d0          = sp["d0"]
        self.psihat      = pd.DataFrame(sp["psihat"], columns=[f"basis_{i}" for i in range(1,self.d0+1)])
        self.etahat      = pd.DataFrame(sp["etahat"].T, columns=[f"scores_{i}" for i in range(1,self.d0+1)])
        self.Yhat        = sp["Yhat"]
        self.Yhat_fix    = sp["Yhat_fix"]
        self.epsilonhat  = sp["epsilonhat"]
        self.u           = sp["u"]

        return self

    # ------------------------------------------------------------

    def predict(self, etahat_values):
        """
        Reconstruct curves from arbitrary scores.
        """
        return self.Ybar + self.psihat @ etahat_values