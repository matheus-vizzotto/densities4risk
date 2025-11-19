import numpy as np
import matplotlib.pyplot as plt

# def inner_product(f, g, du):
#     """
#     f, g = vectors
#     """
#     ip = (f.T @ g)*du
#     return ip

# def L2norm(f,du):
#     norm = np.sqrt(inner_product(f, f, du))
#     return norm

# def sampleCols(A):
#     """
#     Randomly shuffle columns of a matrix.
    
#     Parameters:
#     A : numpy.ndarray
#         Input matrix (2D array)
    
#     Returns:
#     numpy.ndarray : Matrix with columns randomly shuffled
#     """
#     # Get number of columns
#     n_cols = A.shape[1]
    
#     # Generate random permutation of column indices
#     idx = np.random.permutation(n_cols)
    
#     # Create new matrix with shuffled columns
#     M = A[:, idx].copy()
    
#     return M

# def K_dFPC(
#         Y : np.ndarray, 
#         m : int, 
#         p : int, 
#         u : list,
#         du : float,
#         dimension : int,
#         select_ncomp : bool = False,
#         lag_max : int       = 6, 
#         B : int             = 1000, 
#         alpha : float       = 0.10, 
#         ) -> dict:
#     """Function to perform the Dynamic Functional Principal Component Regression analysis from Bathia et al. (2010)

#     Args:
#         Y (np.ndarray): m x T matrix, where m is the sample size for the T transformed densities.
#         m (int): number of grid points.
#         p (int): number of backward parameters.
#         u (list): support for the transformed densities.
#         du (float): space between grid points of the support.
#         dimension (int): dimension the the curve process if select_ncomp = False.
#         select_ncomp (bool): identify the number of dimensions via bootstrap.
#         lag_max (int): maximum number of lags to test with bootstrap if select_ncomp = True.
#         B (int): number of repetitions to be performed in the bootstrap if select_ncomp = True.
#         alpha (float): confidence level for the bootstrap if select_ncomp = True.

#     Returns:
#         dict: dictionary with outputs of the dynamicFPCR:
#             ["Y", "Ybar", "thetahat", "gammahat", "psihat", "etahat", "Yhat", "epsilonhat", "u", "d0", ["bs_pvalues"]]
#     """

#     n = N = Y.shape[1]  # columns in Y

#     # Mean over time
#     Ybar = np.mean(Y, axis=1, keepdims=True)

#     # Deviation from mean
#     Ydev = Y - Ybar

#     ##############################
#     # 3. Creating the matrix core
#     ##############################

#     core = inner_product(Ydev, Ydev, du)
#     Kstar_core0 = core[:(n - p), :(n - p)]

#     Kstar_core = np.zeros((n - p, n - p, p))
#     for k in range(1, p + 1):
#         Kstar_core[:, :, k - 1] = core[k:(n - (p - k)), k:(n - (p - k))]

#     # Summing matrices
#     Kstar_sum = np.sum(Kstar_core, axis=2)

#     # Define Kstar
#     Kstar = (n - p) ** (-2) * np.dot(Kstar_sum, Kstar_core0)

#     ##########################
#     # 4. Eigen decomposition #
#     ##########################

#     eigvals, eigvecs = np.linalg.eig(Kstar)

#     # Keep real parts
#     tol = 1e-4
#     for j in range(min(10, len(eigvals))):
#         if abs(np.imag(eigvals[j])) > tol:
#             print("Complex eigenvalue found.")

#     eigvals = np.real(eigvals)
#     eigvecs = np.real(eigvecs)

#     # Sort eigenvalues and eigenvectors (descending)
#     idx = np.argsort(eigvals)[::-1]
#     thetahat = eigvals[idx]
#     gammahat = eigvecs[:, idx]

#     thetahat_old = thetahat.copy()
#     gammahat_old = gammahat.copy()

#     #############################
#     # Select number of components
#     #############################
#     if select_ncomp:
#         bs_pvalues = np.zeros(lag_max)

#         def sampleCols(A):
#             idx = np.random.permutation(A.shape[1])
#             return A[:, idx]

#         for d0 in range(1, lag_max + 1):
#             thetahatH0 = np.real(thetahat_old[d0])
#             thetahat = np.real(thetahat_old[:d0])
#             gammahat = np.real(gammahat_old[:, :d0])

#             psihat_root = np.dot(Ydev[:, :(n - p)], gammahat)

#             psihat = np.zeros((m, d0))
#             for i in range(d0):
#                 psihat[:, i] = psihat_root[:, i] / L2norm(psihat_root[:, i], du)

#             etahat = inner_product(psihat, Ydev, du)
#             Yhat = Ybar + np.dot(psihat, etahat)

#             epsilonhat = Y - Yhat

#             # Bootstrap
#             bs_thetahat = np.zeros(B)
#             for i in range(B):
#                 bs_epsilon = sampleCols(epsilonhat)
#                 # bs_Y = Yhat_fix + bs_epsilon
#                 bs_Y = Yhat + bs_epsilon
#                 bs_Ybar = np.mean(bs_Y, axis=1, keepdims=True)
#                 bs_Ydev = bs_Y - bs_Ybar

#                 bs_core = inner_product(bs_Ydev, bs_Ydev, du)
#                 bs_Kstar_core0 = bs_core[:(n - p), :(n - p)]

#                 bs_Kstar_core = np.zeros((n - p, n - p, p))
#                 for k in range(1, p + 1):
#                     bs_Kstar_core[:, :, k - 1] = bs_core[k:(n - (p - k)), k:(n - (p - k))]

#                 bs_Kstar_sum = np.sum(bs_Kstar_core, axis=2)
#                 bs_Kstar = (n - p) ** (-2) * np.dot(bs_Kstar_sum, bs_Kstar_core0)
#                 bs_eigvals, _ = np.linalg.eig(bs_Kstar)
#                 bs_eigvals = np.sort(np.real(bs_eigvals))[::-1]
#                 bs_thetahat[i] = bs_eigvals[d0]

#             # Bootstrap p-value
#             bs_pvalues[d0 - 1] = np.sum(bs_thetahat >= thetahatH0) / B

#         d0_candidates = np.where(bs_pvalues < alpha)[0]
#         d0 = d0_candidates[0] + 1 if len(d0_candidates) > 0 else 1

#     else:
#         d0 = dimension

#     ##################################
#     # 5. Estimation of Yhat and epsilons
#     ##################################

#     thetahat = np.real(thetahat_old[:d0])
#     gammahat = np.real(gammahat_old[:, :d0])

#     psihat_root = np.dot(Ydev[:, :(n - p)], gammahat)

#     psihat = np.zeros((m, d0))
#     for i in range(d0):
#         psihat[:, i] = psihat_root[:, i] / L2norm(psihat_root[:, i], du)

#     etahat = inner_product(psihat, Ydev, du)
#     Yhat = Ybar + np.dot(psihat, etahat)

#     epsilonhat = Y - Yhat

#     result = {
#         "Y": Y,
#         "Ybar": Ybar,
#         "thetahat": thetahat,
#         "gammahat": gammahat,
#         "psihat": psihat,
#         "etahat": etahat,
#         "Yhat": Yhat,
#         "epsilonhat": epsilonhat,
#         "u": u,
#         "d0": d0
#     }

#     if select_ncomp:
#         result["bs_pvalues"] = bs_pvalues

#     return result




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
            plt.plot(self.etahat[i, :], alpha=0.6)

        plt.title("KLE-based temporal scores $\\hat{\\eta}$")
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