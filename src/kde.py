import pandas as pd
import numpy as np
from sklearn.model_selection import GridSearchCV, LeaveOneOut, KFold
from typing import List, Union, Optional, Iterable, Any
from scipy.stats import norm, t
from sklearn.neighbors import KernelDensity
import warnings

def get_rot_bandwidth_ck():
    # Kernel-dependent constants for bandwidth selection
    CK_GAUSSIAN_SILVERMAN = 1.06
    CK_GAUSSIAN_SILVERMAN_ADJUSTED = 0.9
    CK_EPANECHNIKOV_SILVERMAN = 2.34

    ROT_BANDWIDTH_CK = {
        "gaussian": {
            "silverman":          CK_GAUSSIAN_SILVERMAN,
            "silverman_adjusted": CK_GAUSSIAN_SILVERMAN_ADJUSTED,

        },
        "epanechnikov": {
            "silverman":          CK_EPANECHNIKOV_SILVERMAN,
            # "silverman_adjusted": CK_GAUSSIAN_SILVERMAN_ADJUSTED * CK_EPANECHNIKOV_SILVERMAN / CK_GAUSSIAN_SILVERMAN # scaled to preserve Silverman/Gaussian ratio
        },
    }

    return ROT_BANDWIDTH_CK

# t_kernel_rot = {
#     '3': 1.192,
#     '4': 1.141,
#     '5': 1.116
#     }

def _gaussian_kernel(u):
    """
    Gaussian kernel.
    """
    return norm.pdf(u)

def _t_student_kernel(u, df):
    """
    Student-t kernel.

    Parameters
    ----------
    df : int
        Degrees of freedom.
    """
    return t.pdf(u, df=df)

def _epanechnikov_kernel(u):
    """
    Epanechnikov kernel.
    """
    return 0.75 * (1 - u**2) * (np.abs(u) <= 1)

def _pilot_density(
        X, 
        h,
        kernel: str = "gaussian",
        **kwargs
        ):
    """
    Compute pilot density estimate at sample points.
    """
    # Convert to NumPy array
    X = np.asarray(X, dtype=float)

    # Pairwise standardized distances:
    # u_ij = (X_i - X_j) / h
    u = (X[:, None] - X[None, :]) / h

    # Gaussian kernel evaluation
    if kernel=="t-student":
        pdf = _t_student_kernel(u, **kwargs)
    else:
        pdf = _gaussian_kernel(u)

    # KDE formula evaluated at sample points
    f_pilot = np.mean(pdf / h, axis=1)

    return f_pilot

def _optimize_sklearn(X, bandwidths, cv, kernel):
    """Internal helper for standard Scikit-Learn kernels."""
    grid = GridSearchCV(
        estimator=KernelDensity(kernel=kernel),
        param_grid={"bandwidth": bandwidths},
        cv=cv,
        n_jobs=-1
    )
    grid.fit(X.reshape(-1, 1))
    return grid.best_params_["bandwidth"]

def _optimize_t_student(X, bandwidths, cv, df):
    """Internal helper for manual Student-t LCV."""

    max_log_lik = -np.inf
    best_h = bandwidths[0]
    
    # If cv is LeaveOneOut, cv.split(X) works natively. 
    # If int/KFold, it also works.
    for h in bandwidths:
        fold_log_liks = []
        for train_idx, val_idx in cv.split(X):
            train, val = X[train_idx], X[val_idx]
            
            # Distance matrix (n_val, n_train)
            u = (val[:, None] - train[None, :]) / h
            # Student-t Density: f_hat = (1/nh) * sum( K(u) )
            densities = t.pdf(u, df=df) / h
            f_hat = np.mean(densities, axis=1)
            
            # Log-likelihood (with numerical floor)
            fold_log_liks.append(np.sum(np.log(np.maximum(f_hat, 1e-300))))
            
        avg_log_lik = np.mean(fold_log_liks)
        if avg_log_lik > max_log_lik:
            max_log_lik = avg_log_lik
            best_h = h
    return best_h

class BandwidthSelector:
    def __init__(
            self,
            X: pd.Series,
            kernel: str = "gaussian",
            **kwargs
        ):
        self.X = X
        self.kernel = kernel
        self.n = X.shape[0]
        self.kwargs = kwargs

    def rot_constants(self) -> List[tuple]:
        rot_elements = get_rot_bandwidth_ck()
        combinations = [
                        (kernel, rule, value) 
                        for kernel, rules in rot_elements.items() 
                        for rule, value in rules.items()
        ]
        
        return combinations

    def rot(
            self,
            rule:         str = "silverman",
            sigma_robust: bool = False,
            return_list:  bool = True,
            kernel:       str = None
        ) -> Union[float, np.array]:
        """
            Compute the Rule of Thumb (ROT) bandwidth for Kernel Density Estimation.

            This method provides a closed-form estimate of the optimal bandwidth (h) 
            based on the assumption that the underlying distribution is approximately 
            normal. It is significantly faster than Cross-Validation but may 
            over-smooth multimodal distributions.

            Parameters
            ----------
            kernel_rule : list of str, default=["gaussian", "silverman"]. See get_rot_bandwidth_ck()
            for available constants.
                A pair defining the kernel type and the specific rule:
                - Index 0: Kernel name (e.g., "gaussian", "epanechnikov").
                - Index 1: Rule variant (e.g., "silverman").
            sigma_robust : bool, default=False
                If True, uses a robust estimator of scale: min(std, IQR/1.349).
                This prevents outliers from artificially inflating the bandwidth 
                (over-smoothing).
            return_list : bool, default=True
                - If True, returns a NumPy array of length 'n' where every entry is 
                the calculated bandwidth.
                - If False, returns a single float value for the global bandwidth.

            Returns
            -------
            bandwidth : float or np.ndarray
                The calculated bandwidth(s).

            Notes
            -----
            The bandwidth is calculated using the formula:
                h = C_k * σ * n^(-1/5)
            
            where:
            - C_k is a kernel-specific constant (e.g., 1.06 for Silverman's Gaussian).
            - σ is the standard deviation (or robust alternative).
            - n is the sample size.

            References
            ----------
            Silverman, B. W. (1986). Density Estimation for Statistics and Data Analysis.
            Scott, D. W. (1992). Multivariate Density Estimation.
        """
        if kernel is None:
            kernel = self.kernel
        c_k = get_rot_bandwidth_ck()[kernel][rule]

        standard_deviation = np.std(self.X, axis=0)
        if sigma_robust:
            iqr_x = np.subtract(*np.percentile(self.X, [75, 25])) / 1.349
            standard_deviation = np.min([standard_deviation, iqr_x]) if iqr_x > 0 else standard_deviation

        bandwidth = c_k * standard_deviation * self.n **(-1/5)

        if return_list:
            bandwidth = np.full(self.n, bandwidth)

        return bandwidth
    
    def adaptive_bandwidth(
            self,
            h: Union[str,float] = "silverman",
            alpha=0.5
        ) -> np.array:
        """
        Compute adaptive local bandwidths using the Abramson (1982) variable kernel method.

        This function implements the "square-root law" for bandwidth adaptation. It 
        scales the global bandwidth 'h' based on a pilot density estimate, allowing 
        the kernel to expand in sparse regions (tails) and contract in dense regions 
        (peaks).

        Parameters
        ----------
        X : ndarray of shape (n_samples,)
            The data points at which to evaluate the pilot density and compute bandwidths.
        h : float
            The global (initial) bandwidth. Usually determined by a rule of thumb 
            (e.g., Silverman's or Scott's rule).
        alpha : float, default=0.5
            The sensitivity parameter. 
            - alpha = 0.5: Abramson’s optimal value, which reduces the bias of the 
                KDE from O(h²) to O(h⁴).
            - alpha = 0: Results in a fixed bandwidth h.
            - 0 < alpha < 1: Varying degrees of local adaptation.

        Returns
        -------
        bandwidths : ndarray of shape (n_samples,)
            The local bandwidth values λ(X_i) * h for each point in X.

        Notes
        -----
        The adaptive bandwidth is calculated as:
            h_i = h * (f_pilot(X_i) / g)^(-alpha)

        where 'g' is the geometric mean of the pilot density f_pilot(X):
            log(g) = (1/n) * Σ log(f_pilot(X_i))

        This geometric mean normalization ensures that the product of the local 
        bandwidths remains proportional to the global bandwidth h, maintaining 
        the overall smoothing scale.

        References
        ----------
        Abramson, I. S. (1982). On bandwidth variation in kernel estimates-a 
        square root law. Annals of Statistics, 10(4), 1217-1223.
        """

        if isinstance(h, str):
            h = self.rot(kernel="gaussian", rule="silverman", return_list=False)

        f_pilot = _pilot_density(self.X, h=h, kernel=self.kernel, **self.kwargs)

        # Numerical safety
        eps = np.finfo(float).eps
        f_pilot = np.maximum(f_pilot, eps)

        g = np.exp(np.mean(np.log(f_pilot)))

        bandwidths = h * (f_pilot / g) ** (-alpha)

        return bandwidths
        
    def cross_validate_bandwidth(
            self,
            bandwidths_grid: Optional[Iterable[float]] = None,
            cv: Union[int, str] = 5,
            # kernel: str = "gaussian",
            return_list: bool = True,
            **kernel_kwargs: Any
        ) -> Union[float, np.array]:
        """
        Select optimal KDE bandwidth via Likelihood Cross-Validation (LCV).

        Finds the global scaling factor (h) that maximizes the log-likelihood of 
        out-of-sample observations. This data-driven approach automates the 
        bias-variance tradeoff by optimizing the model's predictive density.

        Parameters
        ----------
        bandwidths_grid : array-like, optional
            Candidate bandwidth values to evaluate. If None, a logarithmic grid 
            is automatically constructed, centered around Silverman's Rule of 
            Thumb (ROT) and scaled by the data's robust volatility.
        cv : int or "LOO", default=5
            The cross-validation strategy:
            - int: Number of folds for K-Fold CV (Recommended for large n).
            - "LOO": Leave-One-Out CV (Classic approach, high computational cost).
        kernel : str, default="gaussian"
            Kernel function for density estimation. Supports Scikit-Learn kernels 
            (gaussian, epanechnikov, etc.) and custom "t-student" implementation.
        return_list : bool, default=True
            - If True: Returns an ndarray of length self.n filled with h_opt.
            - If False: Returns h_opt as a single float.
        **kernel_kwargs : Any
            Additional parameters passed to custom kernels. 
            Example: df=3 for 't-student'.

        Returns
        -------
        bandwidth : float or np.ndarray
            The optimal bandwidth(s) that maximize the cross-validated 
            log-likelihood.

        Notes
        -----
        Likelihood CV is asymptotically equivalent to minimizing the 
        Kullback-Leibler divergence between the estimated and true densities. 
        Unlike Least-Squares CV, LCV is particularly effective for heavy-tailed 
        distributions, though it can be sensitive to extreme outliers.

        Warnings
        --------
        Issues a UserWarning if the selected bandwidth lies on the boundary of 
        the provided or generated grid, suggesting the search space is insufficient.
        """

        # 1. Input Validation and Preparation
        X = np.asarray(self.X, dtype=float).flatten()
        if X.size < 2:
            raise ValueError("Bandwidth estimation requires at least two data points.")
        
        # 2. Automated Grid Generation (Data-Dependent)
        if bandwidths_grid is None:
            std_x = np.std(X, ddof=1)
            iqr_x = np.subtract(*np.percentile(X, [75, 25])) / 1.349
            sigma_robust = np.min([std_x, iqr_x]) if iqr_x > 0 else std_x
            h_rot = 0.9 * sigma_robust * (X.size ** -0.2)
            bandwidths_grid = np.logspace(np.log10(h_rot/10), np.log10(h_rot*10), 50)

        # 3. CV Strategy Dispatch
        cv_strategy = LeaveOneOut() if cv == "LOO" else KFold(n_splits=cv, shuffle=True)

        # 4. Optimization Path Dispatch
        if self.kernel == "t-student":
            bandwidth = _optimize_t_student(X, bandwidths_grid, cv_strategy, self.kwargs.get('df'))
        else:
            bandwidth = _optimize_sklearn(X, bandwidths_grid, cv_strategy, self.kernel)

        # 5. Boundary Warning
        if np.isclose(bandwidth, np.min(bandwidths_grid)) or np.isclose(bandwidth, np.max(bandwidths_grid)):
            
            warnings.warn(
                f"Optimal bandwidth ({bandwidth:.4f}) is at the boundary of the grid for kernel '{self.kernel}'.", 
                UserWarning
            )

        if return_list:
            bandwidth = np.full(self.n, bandwidth)
        return bandwidth
    
def df_bandwidth_selector(
        X: pd.DataFrame,
        kernel: str = "gaussian",
        method: str = "rot",
        # return_list: bool = True,
        **kwargs
    ) -> pd.DataFrame:
 
   
    results = {}
        
    # 1. Prepare Init Parameters
    # We pull 'df' out of kwargs if it exists to pass it to the constructor
    init_params = {'kernel': kernel}
    if kernel == "t-student" and "df" in kwargs:
        init_params["df"] = kwargs.pop("df") 
    
    # 2. Remaining kwargs are specific to the method (e.g., alpha, cv, rule)
    method_params = kwargs

    for col_name in X.columns:
        # Drop NaNs to perform calculation
        clean_col = X[col_name].dropna()
        
        if clean_col.empty:
            results[col_name] = np.full(len(X), np.nan)
            continue

        # Initialize the selector with kernel and potentially df
        selector = BandwidthSelector(clean_col, **init_params)
        
        # Get the specific method (rot, adaptive_bandwidth, etc.)
        calc_func = getattr(selector, method)
        
        # Execute the calculation
        res = calc_func(**method_params)
        
        # Map results back to the original index to preserve NaN positions
        if np.isscalar(res):
            results[col_name] = np.full(len(X), res)
        else:
            full_series = pd.Series(np.nan, index=X.index)
            full_series.loc[clean_col.index] = res
            results[col_name] = full_series.values
                
    return pd.DataFrame(results, index=X.index)
    
class KDE:
    """
    Description
    ---------
    Kernel Density Estimation for (functional) data.

    Supports multiple kernels and observation-specific bandwidths.

    Example
    ---------
    >>>x1 = np.random.randn(200)
    >>>x2 = np.random.standard_t(df=3, size=200)
    >>>grid = np.linspace(-4, 4, 400)

    >>>kde = FDA_KDE(kernel="gaussian", bandwidth=0.3)
    >>>density = kde.fit_transform(x1, grid)
    >>>density_grid = kde.grid

    >>>params = {"kernel": "t_student", "bandwidth": "scorr", "df": 3, "adaptive": True}
    >>>t_kde = FDA_KDE(**params)
    >>>t_density = t_kde.fit_transform(x2)
    """

    def __init__(
            self, 
            kernel="gaussian", 
            **kernel_params
            ):
        """
        Parameters
        ----------
        kernel : str
            Kernel name. Available: 'gaussian', 't_student', 'epanechnikov'.
        bandwidth : float or array-like
            Scalar or observation-specific bandwidth(s).
        **kernel_params :
            Additional parameters passed to the kernel
            (e.g., df for Student-t kernel).
        """
        self.kernel = kernel
        self.kernel_params = kernel_params
        # self.df = df

        # Kernel registry (instance-level, extensible)
        self._kernel_map = {
            "gaussian": _gaussian_kernel,
            "t-student": _t_student_kernel,
            "epanechnikov": _epanechnikov_kernel,
        }

        # Attributes set after fitting
        self.X    = None
        self.grid = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def transform(
            self, 
            X: pd.Series,
            h: np.array,
            grid: np.array = None,
            ensure_integral_constraint: bool=True
        ):
        """
        Evaluate the KDE on a given grid.

        Parameters
        ----------
        grid : array-like, shape (n_grid,)
            Evaluation points.

        Returns
        -------
        density : ndarray, shape (n_grid,)
            Estimated density values.
        """
        self.X    = X.values
        self.n    = len(X)
        if grid is None:
            grid  = np.linspace(X.min(), X.max(), self.n)
        self.grid = grid

        grid = np.asarray(grid, dtype=float)

        # Standardized distances: u_ij = (grid_i - X_j) / h_j
        u = (self.grid[:, None] - self.X[None, :]) / h[None, :]
        pdf = self._evaluate_kernel(u)
        # KDE: mean_j K(u_ij) / h_j
        density = np.mean(pdf / h[None, :], axis=1).astype(float)

        if ensure_integral_constraint:
            density /= np.trapezoid(y=density, x=self.grid, axis=0)

        return grid, density

    # ------------------------------------------------------------------
    # Kernel evaluation
    # ------------------------------------------------------------------
 
    def _evaluate_kernel(
            self, 
            u : np.array
        ) -> np.array:
        """
        Evaluate the selected kernel on standardized distances.

        This method dispatches the evaluation of the kernel function
        specified by ``self.kernel`` using the internal kernel registry.
        The kernel is evaluated pointwise on the standardized distances

            u_ij = (x_i - X_j) / h_j,

        where ``x_i`` are evaluation points, ``X_j`` are training
        observations, and ``h_j`` are observation-specific bandwidths.

        Parameters
        ----------
        u : ndarray, shape (n_grid, n_samples)
            Matrix of standardized distances between evaluation points
            and training observations.

        Returns
        -------
        pdf : ndarray, shape (n_grid, n_samples)
            Kernel values evaluated at ``u``.

        Raises
        ------
        ValueError
            If ``self.kernel`` is not a recognized kernel name.

        Notes
        -----
        This method assumes that the selected kernel integrates to one.
        """
        try:
            kernel_fn = self._kernel_map[self.kernel]
        except KeyError:
            raise ValueError(
                f"Unknown kernel '{self.kernel}'. "
                f"Available kernels: {list(self._kernel_map)}"
            )

        kde = kernel_fn(u, **self.kernel_params)
        
        return kde
    
def df_to_kde(
        X: pd.DataFrame,
        h: pd.DataFrame,
        kernel: str = "gaussian",
        m=5001,
        **kwargs
        ) -> pd.DataFrame:   

        base_grid = np.linspace(X.min().min(), X.max().max(), m)

        results = {}
        results_grids = {}
        for col_name in X.columns:
                col = X.loc[:, col_name]
                h_col = h.loc[:, col_name].values
                kde = KDE(kernel=kernel, **kwargs)
                grid, density = kde.transform(col, h=h_col, grid=base_grid)
                results_grids[col_name] = grid
                results[col_name] = density

        df_grids = pd.DataFrame(results_grids)
        df_densities = pd.DataFrame(results)

        return df_grids, df_densities