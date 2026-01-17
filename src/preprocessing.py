import pandas as pd
import numpy as np
from scipy.interpolate import interp1d, PchipInterpolator
import matplotlib.pyplot as plt
from typing import Tuple, List
from scipy.integrate import trapezoid
from scipy.stats import norm, t


def align_and_normalize_density(x_obs, f_obs, x_hat, f_hat,
                                x_common=None, n_points=256,
                                interp_kind='linear', eps=1e-12):
    """
    Interpolate two density curves onto a common grid, clip negatives,
    and renormalize so each integrates to 1.

    Returns: x_common, f_obs_common, f_hat_common
    """

    # Choose a common grid covering both supports
    if x_common is None:
        lo = min(np.nanmin(x_obs), np.nanmin(x_hat))
        hi = max(np.nanmax(x_obs), np.nanmax(x_hat))
        x_common = np.linspace(lo, hi, n_points)

    # safe interpolation: fill outside range with 0
    interp_obs = interp1d(x_obs, f_obs, kind=interp_kind,
                          bounds_error=False, fill_value=0.0)
    interp_hat = interp1d(x_hat, f_hat, kind=interp_kind,
                          bounds_error=False, fill_value=0.0)

    f_obs_c = interp_obs(x_common)
    f_hat_c = interp_hat(x_common)

    # sanitize: remove NaN/inf and negatives
    f_obs_c = np.nan_to_num(f_obs_c, nan=0.0, posinf=0.0, neginf=0.0)
    f_hat_c = np.nan_to_num(f_hat_c, nan=0.0, posinf=0.0, neginf=0.0)
    f_obs_c = np.maximum(f_obs_c, 0.0)
    f_hat_c = np.maximum(f_hat_c, 0.0)

    # renormalize each to integrate to 1
    area_obs = np.trapz(f_obs_c, x_common)
    area_hat = np.trapz(f_hat_c, x_common)

    # avoid division by zero
    if area_obs <= eps:
        raise ValueError("Observed density integrates to ~0 after interpolation.")
    if area_hat <= eps:
        raise ValueError("Forecast density integrates to ~0 after interpolation.")

    f_obs_c = f_obs_c / area_obs
    f_hat_c = f_hat_c / area_hat

    return x_common.copy(), f_obs_c.copy(), f_hat_c.copy()


def align_densities(
        densities_support : pd.DataFrame,
        densities : pd.DataFrame,
        fc_densities_support : pd.DataFrame,
        fc_densities : pd.DataFrame,
        cols : List,
        n : int = 256
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Align observed and forecast density functions onto a common grid
    at each time index and renormalize them for divergence evaluation.

    This function is designed for density forecast evaluation under
    time-varying supports, as in the modified log-quantile-density
    transformation (LQDT) framework of Kokoszka et al. (2019).
    For each selected column (typically corresponding to a forecast
    origin or time index), the observed and forecast densities are
    interpolated onto a common grid defined on the union of their
    supports. The resulting curves are truncated at zero and
    renormalized to integrate to one.

    Parameters
    ----------
    densities_support : pandas.DataFrame
        Support grids of the observed densities. Each column corresponds
        to a time index and contains the support on which the observed
        density is defined.

    densities : pandas.DataFrame
        Observed density values evaluated on `densities_support`.
        Must have the same shape and column labels as
        `densities_support`.

    fc_densities_support : pandas.DataFrame
        Support grids of the forecast densities. Each column corresponds
        to a time index and may differ from the observed support.

    fc_densities : pandas.DataFrame
        Forecast density values evaluated on `fc_densities_support`.
        Must have the same shape and column labels as
        `fc_densities_support`.

    cols : list
        List of column names (time indices) for which densities are
        aligned and normalized.

    Returns
    -------
    df_supp : pandas.DataFrame
        Common support grids obtained for each column after alignment.
        Each column contains the grid on which both observed and
        forecast densities are evaluated.

    df_f : pandas.DataFrame
        Observed densities interpolated onto the common grids and
        renormalized to integrate to one.

    df_fhat : pandas.DataFrame
        Forecast densities interpolated onto the common grids and
        renormalized to integrate to one.

    Notes
    -----
    * Alignment is performed independently for each column (time index),
      allowing for time-varying supports.
    * Interpolation, truncation at zero, and renormalization are handled
      by the auxiliary function `align_and_normalize_density`.
    * The returned densities are suitable for the computation of
      divergence-based accuracy measures such as the Kullback–Leibler
      divergence or Jensen–Shannon divergence.
    """

    base = pd.DataFrame(np.full([n,len(cols)], np.nan), columns = cols)

    df_supp = base.copy()
    df_f = base.copy()
    df_fhat = base.copy()
    for col in cols:
        supp_col, f_col, f_hat_col = align_and_normalize_density(densities_support.loc[:,col], densities.loc[:,col], fc_densities_support.loc[:,col], fc_densities.loc[:,col])
        df_supp.loc[:,col] = supp_col
        df_f.loc[:,col] = f_col
        df_fhat.loc[:,col] = f_hat_col
    
    return df_supp, df_f, df_fhat


def weigh_norm_densities(
        df_densities: pd.DataFrame, 
        support: np.array, 
        norm:int ='area'
        ) -> pd.DataFrame:
    """
    Reweight and normalize a collection of discretized density functions.

    Each density is multiplied pointwise by the cross-sectional mean density
    and then normalized either to unit integral or by a user-specified constant.

    Parameters
    ----------
    df_densities : pd.DataFrame
        DataFrame of shape (m, n) where each column represents a density
        evaluated on a common grid.
    support : np.array
        One-dimensional grid corresponding to the rows of `df_densities`,
        used for numerical integration.
    norm : {'area', int}, default='area'
        Normalization method. If 'area', densities are normalized to integrate
        to one using the trapezoidal rule. Otherwise, `norm` is treated as a
        fixed normalization constant.

    Returns
    -------
    pd.DataFrame
        DataFrame of reweighted and normalized densities with the same shape
        as `df_densities`.
    """
    df2 = df_densities.copy()
    for col in df2.columns:
        # Calculate the transformed density
        transformed = df_densities.loc[:,col] * df_densities.mean(axis=1)
        
        if norm == 'area':
            # Normalize by area to integrate to 1
            area = trapezoid(transformed, support)
            df2.loc[:,col] = transformed / area
        else:
            # Normalize by given constant
            df2.loc[:,col] = transformed / norm

    return df2


########################################################################
######################## DENSITY ESTIMATION ############################ 
########################################################################

def rule_of_thumb_bandwidth(
        x, 
        rule="silverman"
    ) -> float:
    """
    Estimate bandwidth using a rule-of-thumb method.

    Implements Scott's rule and Silverman's rule for
    univariate kernel density estimation.

    Parameters
    ----------
    x : array-like, shape (n_samples,)
        One-dimensional sample.
    rule : {"silverman", "scott"}, default="silverman"
        Bandwidth selection rule.

    Returns
    -------
    h : float
        Estimated bandwidth.

    Notes
    -----
    Scott's rule:
        h = sigma * n^{-1/5}

    Silverman's rule:
        h = 0.9 * min(sigma, IQR / 1.34) * n^{-1/5}

    where sigma is the sample standard deviation and IQR is the
    interquartile range.

    References
    ----------
    Scott, D. W. (1992). Multivariate Density Estimation.
    Silverman, B. W. (1986). Density Estimation for Statistics
    and Data Analysis.
    """
    x = np.asarray(x, dtype=float)
    if x.ndim != 1:
        raise ValueError("x must be one-dimensional")

    n = len(x)
    if n < 2:
        raise ValueError("At least two observations are required")

    sigma = np.std(x, ddof=1)

    if rule.lower() == "scott":
        return sigma * n ** (-1 / 5)

    elif rule.lower() == "silverman":
        q75, q25 = np.percentile(x, [75, 25])
        iqr = q75 - q25
        scale = min(sigma, iqr / 1.34)
        return 0.9 * scale * n ** (-1 / 5)

    else:
        raise ValueError("rule must be 'silverman' or 'scott'")
    
class KernelDensityEstimation:
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

    >>>t_kde = FDA_KDE(kernel="t_student", bandwidth='scott', df = 3)
    >>>t_density = t_kde.fit_transform(t_samples, adaptive=True)
    """

    def __init__(
            self, 
            kernel="gaussian", 
            bandwidth=1.0, 
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
        self.bandwidth = bandwidth
        self.kernel_params = kernel_params

        # Kernel registry (instance-level, extensible)
        self._kernel_map = {
            "gaussian": self._gaussian_kernel,
            "t_student": self._t_student_kernel,
            "epanechnikov": self._epanechnikov_kernel,
        }

        # Attributes set after fitting
        self.X_   = None
        self.h_   = None
        self.grid = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def fit(
            self, 
            X,
            adaptive:bool
        ):
        """
        Store training data and validate bandwidth.

        Parameters
        ----------
        X : array-like, shape (n_samples,)
            Training observations.

        Returns
        -------
        self
        """
        X = np.asarray(X, dtype=float)
        if X.ndim != 1:
            raise ValueError("X must be one-dimensional")
        
        self.X_ = X
        
        # Step 1: global bandwidth
        if isinstance(self.bandwidth, str):
            h0 = rule_of_thumb_bandwidth(X, self.bandwidth)
        else:
            h0 = float(self.bandwidth)

        # Step 2: adaptive adjustment (optional)
        if adaptive:
            self.h_ = self._adaptive_bandwidth(X, h0)
        else:
            self.h_ = np.full_like(X, h0)

        return self

    def transform(
            self, 
            grid
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
        self._check_is_fitted()

        grid = np.asarray(grid, dtype=float)

        # Standardized distances: u_ij = (grid_i - X_j) / h_j
        u = (grid[:, None] - self.X_[None, :]) / self.h_[None, :]

        pdf = self._evaluate_kernel(u)

        # KDE: mean_j K(u_ij) / h_j
        density = np.mean(pdf / self.h_[None, :], axis=1)

        return density

    def fit_transform(
            self, 
            X : np.array, 
            m : int = 256,
            adaptive : bool = False):
        """
        Fit the model and evaluate the KDE on a grid.
        """
        grid = np.linspace(X.min(), X.max(), m)

        self.grid = grid
        return self.fit(X, adaptive=adaptive).transform(grid)

    # ------------------------------------------------------------------
    # Kernel evaluation
    # ------------------------------------------------------------------


    def _gaussian_kernel(self, u):
        """
        Gaussian kernel.
        """
        return norm.pdf(u)

    def _t_student_kernel(self, u, df=5):
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
        Bandwidth scaling is applied outside this function when computing
        the KDE estimator.
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

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------

    @staticmethod
    def _validate_bandwidth(bandwidth, X):
        """
        Validate and broadcast bandwidth to match X.
        """
        h = np.asarray(bandwidth, dtype=float)

        if h.ndim == 0:
            h = np.full_like(X, h)

        if h.ndim != 1 or len(h) != len(X):
            raise ValueError(
                "bandwidth must be a scalar or an array with "
                "the same length as X"
            )

        if np.any(h <= 0):
            raise ValueError("bandwidth values must be positive")

        return h

    def _check_is_fitted(self):
        if self.X_ is None or self.h_ is None:
            raise RuntimeError("The estimator is not fitted yet.")
        
    def _pilot_density(self, X, h):
        """
        Compute pilot density estimate at sample points.
        """
        u = (X[:, None] - X[None, :]) / h
        pdf = self._gaussian_kernel(u)  # pilot is Gaussian by convention
        return np.mean(pdf / h, axis=1)
        
    def _adaptive_bandwidth(self, X, h):
        """
        Compute adaptive (Abramson) bandwidths.
        """
        f_pilot = self._pilot_density(X, h)

        # Numerical safety
        eps = np.finfo(float).eps
        f_pilot = np.maximum(f_pilot, eps)

        g = np.exp(np.mean(np.log(f_pilot)))

        return h * (f_pilot / g) ** (-0.5)