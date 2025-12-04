import numpy as np
import matplotlib.pyplot as plt
import pywt
from typing import Tuple



################################ NOVA CLASSE #####################

# TODO: K_dFPC.forecast_scores() gerando erro: IndexError: index 0 is out of bounds for axis 0 with size 0

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

    def fit(self, lag_max, B, alpha, du, p, m, u,
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
        Y = self.Y
        n = N = Y.shape[1]
        self.u = u

        # Mean and deviations
        Ybar = np.mean(Y, axis=1, keepdims=True)
        Ydev = Y - Ybar

        # --------------------------------------------------------
        # Step 1 — Build core matrices
        # --------------------------------------------------------

        core = inner_product(Ydev, Ydev, du)
        Kstar_core0 = core[:(n - p), :(n - p)]

        Kstar_core = np.zeros((n - p, n - p, p))
        for k in range(1, p + 1):
            Kstar_core[:, :, k - 1] = core[k:(n - (p - k)), k:(n - (p - k))]

        Kstar_sum = np.sum(Kstar_core, axis=2)
        Kstar = (n - p) ** (-2) * np.dot(Kstar_sum, Kstar_core0)

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
        # Step 3 — Select number of components (optional bootstrap)
        # --------------------------------------------------------

        if select_ncomp:
            bs_pvalues = np.zeros(lag_max)

            def sampleCols(A):
                return A[:, np.random.permutation(A.shape[1])]

            for d0 in range(1, lag_max + 1):
                thetahatH0 = thetahat_old[d0]
                gammahat = gammahat_old[:, :d0]

                psihat_root = np.dot(Ydev[:, :(n - p)], gammahat)

                psihat = np.zeros((m, d0))
                for i in range(d0):
                    psihat[:, i] = psihat_root[:, i] / L2norm(psihat_root[:, i], du)

                etahat = inner_product(psihat, Ydev, du)
                Yhat = Ybar + np.dot(psihat, etahat)
                epsilonhat = Y - Yhat

                # Bootstrap eigenvalues
                bs_thetahat = np.zeros(B)
                for b in range(B):
                    bs_epsilon = sampleCols(epsilonhat)
                    bs_Y = Yhat + bs_epsilon

                    bs_Ybar = np.mean(bs_Y, axis=1, keepdims=True)
                    bs_Ydev = bs_Y - bs_Ybar

                    bs_core = inner_product(bs_Ydev, bs_Ydev, du)
                    bs_Kstar_core0 = bs_core[:(n - p), :(n - p)]

                    bs_Kstar_core = np.zeros((n - p, n - p, p))
                    for k in range(1, p + 1):
                        bs_Kstar_core[:, :, k - 1] = bs_core[k:(n - (p - k)),
                                                             k:(n - (p - k))]

                    bs_Kstar_sum = np.sum(bs_Kstar_core, axis=2)
                    bs_Kstar = (n - p) ** (-2) * np.dot(bs_Kstar_sum, bs_Kstar_core0)

                    bs_eigs, _ = np.linalg.eig(bs_Kstar)
                    bs_thetahat[b] = np.sort(np.real(bs_eigs))[::-1][d0]

                bs_pvalues[d0 - 1] = np.mean(bs_thetahat >= thetahatH0)

            d0_candidates = np.where(bs_pvalues < alpha)[0]
            d0 = d0_candidates[0] + 1 if len(d0_candidates) > 0 else 1

            self.bs_pvalues = bs_pvalues

        else:
            d0 = dimension

        # --------------------------------------------------------
        # Step 4 — Final estimation
        # --------------------------------------------------------

        thetahat = thetahat_old[:d0]
        gammahat = gammahat_old[:, :d0]

        psihat_root = np.dot(Ydev[:, :(n - p)], gammahat)

        psihat = np.zeros((m, d0))
        for i in range(d0):
            psihat[:, i] = psihat_root[:, i] / L2norm(psihat_root[:, i], du)

        etahat = inner_product(psihat, Ydev, du)
        Yhat = Ybar + np.dot(psihat, etahat)
        epsilonhat = Y - Yhat

        # Store results
        self.Ybar = Ybar
        self.thetahat = thetahat
        self.gammahat = gammahat
        self.psihat = psihat
        self.etahat = etahat
        self.Yhat = Yhat
        self.epsilonhat = epsilonhat
        self.d0 = d0

        # Prediction storage (Ybar + ψ * η_pred_)
        self.fitted_values = self.Ybar + self.psihat @ self.etahat

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

    def forecast_scores(self, h=1, max_var_lag=3, var_type="c"):
        """
        Forecast the next functional observation using ARIMA (if d0=1)
        or VAR (if d0>1), following the logic of the original R code.

        Parameters
        ----------
        h : int, default=1
            Forecast horizon.
        max_var_lag : int, default=3
            Maximum lag for VAR model.
        var_type : {"const", "none", "trend", "both"}
            Equivalent to R's VAR type. Automatically translated
            to statsmodels-compatible trend codes.

        Returns
        -------
        etahat_pred_val : ndarray (d0, h)
            Forecasted score vectors.
        Yhat_future : ndarray (m, h)
            Forecasted functional curves.
        """

        import numpy as np
        from statsmodels.tsa.arima.model import ARIMA
        from statsmodels.tsa.api import VAR

        if self.etahat is None:
            raise ValueError("Model must be fitted before forecasting.")

        score_object = self.etahat
        d0, T = score_object.shape

        # -------------------------------------------------------
        # Case 1: Univariate ARIMA
        # -------------------------------------------------------
        if d0 == 1:
            series = score_object.flatten()
            model = ARIMA(series, order=(1, 0, 0))
            fitted = model.fit()
            forecast_vals = fitted.forecast(steps=h)
            etahat_pred_val = forecast_vals.reshape(1, h)

        # -------------------------------------------------------
        # Case 2: Multivariate VAR
        # -------------------------------------------------------
        else:
            score_Td = score_object.T  # shape: T × d0

            var_model = VAR(score_Td)
            sel = var_model.select_order(maxlags=max_var_lag)

            chosen_lag = sel.selected_orders["aic"]
            if chosen_lag is None:
                chosen_lag = 1
            chosen_lag = min(chosen_lag, max_var_lag)

            var_fit = var_model.fit(maxlags=chosen_lag, trend=var_type)

            # FIX: use endog, not y
            pred = var_fit.forecast(var_fit.endog, steps=h)

            etahat_pred_val = pred.T

        return etahat_pred_val
    

######################## WAVELETS

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

# DimEstWav/DimEst_boot.m
def DimEst_boot(Y, NREP, B, p, N, wname):
    """
    Python translation of MATLAB DimEst_boot.

    Parameters
    ----------
    Y : array (nt x n)
        Functional observations
    NREP : int
        Number of bootstrap replications
    B : array (J x d0)
        Eigenvector matrix (wavelet coefficients of estimated eigenfunctions)
    p : int
        Maximum lag
    N : int
        Wavelet decomposition level
    wname : str
        Wavelet name (e.g. 'db2')

    Returns
    -------
    vd0p1_boot : array (NREP,)
        Bootstrap estimates of eigenvalue_{d0+1}
    """

    J, d0 = B.shape
    nt, n = Y.shape

    vd0p1_boot = np.zeros(NREP)
    Yhat = np.zeros((nt, n))
    H = np.zeros((nt, d0))
    A = np.zeros((J, n))
    Yboot = np.zeros((nt, n))
    Aboot = np.zeros((J, n))

    # =====================================================================
    # Compute wavelet decomposition of Y
    # =====================================================================
    # Get coefficient slicing lengths using first column
    coeffs0 = pywt.wavedec(Y[:, 0], wname, level=N)
    lengths = [len(c) for c in coeffs0]

    for ii in range(n):
        coeffs = pywt.wavedec(Y[:, ii], wname, level=N)
        A[:, ii] = np.concatenate(coeffs)

    mu_A = A.mean(axis=1, keepdims=True)
    C = A - mu_A

    # =====================================================================
    # Reconstruct mean function using inverse wavelet transform
    # =====================================================================
    def waverec_from_vector(v):
        out = []
        idx = 0
        for L in lengths:
            out.append(v[idx:idx+L])
            idx += L
        return pywt.waverec(out, wname)

    mu_hat = waverec_from_vector(mu_A[:, 0])

    # =====================================================================
    # Reconstruct eigenfunctions H(:,k)
    # =====================================================================
    for ii in range(d0):
        H[:, ii] = waverec_from_vector(B[:, ii])

    # =====================================================================
    # Compute fitted curves Yhat
    # =====================================================================
    for ii in range(n):
        Yhat[:, ii] = mu_hat.copy()
        for k in range(d0):
            coef = C[:, ii].T @ B[:, k]   # inner product
            Yhat[:, ii] += coef * H[:, k]

    mEps_hat = Y - Yhat

    # =====================================================================
    # Bootstrap replications
    # =====================================================================
    for jj in range(NREP):

        # Resample residual curves with replacement
        for ii in range(n):
            idx = np.random.randint(0, n)
            Yboot[:, ii] = Yhat[:, ii] + mEps_hat[:, idx]

            coeffs = pywt.wavedec(Yboot[:, ii], wname, level=N)
            Aboot[:, ii] = np.concatenate(coeffs)

        mu_Aboot = Aboot.mean(axis=1, keepdims=True)
        Cb = Aboot - mu_Aboot

        C1 = Cb[:, :n-p]
        D1 = np.zeros((n-p, n-p))

        for k in range(1, p+1):
            Ct = Cb[:, k:(n-p+k)]
            D1 += Ct.T @ Ct

        Dboot = C1 @ D1 @ C1.T / ((n - p)**2)

        # eigen-decomposition
        Lboot, Bboot = np.linalg.eig(Dboot)

        # MATLAB takes Lboot(d0+1,d0+1), but eigenvalues are unordered.
        # We sort them descending to match MATLAB's behavior.
        L_sorted = np.sort(Lboot)[::-1]

        vd0p1_boot[jj] = L_sorted[d0]

    return vd0p1_boot

def W_dFPC(
        Y: np.ndarray,
        nt: int,
        N: int = 3,
        wavelet_: str = 'db2',
        p = 5
    ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute Wavelet-domain Functional Principal Components (WdFPC)
    for a sample of functional data.

    Parameters
    ----------
    Y : np.ndarray
        Functional observations of shape (nt, n), where:
        - nt is the number of grid points for each curve.
        - n  is the number of curves.
        Each column Y[:, i] is one function sampled on a common grid.

    nt : int
        Number of grid points of the functional domain.
        Must match Y.shape[0].

    N : int, default=3
        Wavelet decomposition level used in the discrete wavelet transform
        (DWT). Higher values give coarser approximations.

    wavelet : str, default='db2'
        Wavelet family to use for the DWT. Must be any valid PyWavelets
        wavelet name (e.g., 'db2', 'haar', 'coif1', ...).

    p : int, default=5
        Number of lags for the covariance operator.

    Returns
    -------
    L : np.ndarray
        Eigenvalues of the estimated covariance operator in the wavelet
        domain, sorted in decreasing order. Shape: (J,).

    B : np.ndarray
        Corresponding eigenvectors (principal directions) in the wavelet
        coefficient space. Shape: (J, J), where J is the total number of
        stacked wavelet coefficients.

    Notes
    -----
    - The function uses the **stacked wavelet coefficients** of each curve
      as features, following the methodology of wavelet-domain PCA.

    - The sample size `n` is inferred from `Y.shape[1]`.

    - The returned eigenvectors B span the wavelet coefficient space, not
      the original functional domain.

    - The temporal dependence structure is incorporated by forming a
      lagged autocovariance estimator using p=5 lags.

    - The function assumes Y has no missing values and is properly aligned
      on a common grid.
    """
        
    n = Y.shape[1]  # infer sample size from data

    # compute total number of stacked DWT coefficients
    sizes, total = wavedec_sizes(nt, wavelet_, N)
    J = total

    # matrix of DWT coefficient vectors
    A = np.zeros((J, n))

    for ii in range(n):
        coeffs = pywt.wavedec(Y[:, ii], wavelet=wavelet_, level=N)
        A[:, ii] = np.concatenate(coeffs)

    # mean-adjust coefficients
    mu_A = np.mean(A, axis=1, keepdims=True)
    C = A - mu_A

    # compute lagged covariance operator estimate
    C1 = C[:, :n - p]
    D1 = np.zeros((n - p, n - p))

    for k in range(1, p + 1):
        Ct = C[:, k:(n - p + k)]
        D1 += Ct.T @ Ct

    D = C1 @ D1 @ C1.T / ((n - p) ** 2)

    # eigendecomposition
    L, B = np.linalg.eig(D)

    # sort eigenvalues descending
    idx = np.argsort(-L)

    return L[idx], B[:, idx]