import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from typing import Dict, Union, List, Tuple, Optional
import numpy as np
from statsmodels.tsa.api import VAR
from statsmodels.tsa.vector_ar.var_model import VARResults
import pmdarima as pm
from statsmodels.tsa.stattools import adfuller, kpss#, phillips_perron
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.graphics.tsaplots import plot_acf

from src.preprocessing import align_densities
from src.transformations import (
                            mLQDT,
                            obtain_densities_from_lqd
                            )
from src.dynamicFPC import K_dFPC
from scipy.stats import norm

def adf_test(x):
    res = adfuller(x, autolag='AIC')
    # print('ADF Statistic: %f' % res[0])
    # print('p-value: %f' % res[1])
    # print('Critical Values:')
    # for k,v in res[4].items():
    #     print(f'   {k}: {v}')
    if not res[1] < 0.05:
        print('\t ADF => reject H0 (unit root) ? ', res[1] < 0.05)
        plt.figure()
        plt.plot(x)
        plt.show()
    print()

def kpss_test(x, regression='c'):
    statistic, p_value, lags, crit = kpss(x, regression=regression)
    # print('KPSS Statistic: %f' % statistic)
    # print('p-value: %f' % p_value)
    # print('Critical Values:')
    # for k,v in crit.items():
    #     print(f'   {k}: {v}')
    if p_value < 0.05:
        print('\t KPSS => reject H0 (stationary) ? ', p_value < 0.05)
        plt.figure()
        plt.plot(x)
        plt.show()
    print()

def check_autocorrelation(
        series, 
        name="Series", 
        lags=10,
        verbose=False
        ):
    # Perform the test
    results = acorr_ljungbox(series, lags=[lags], return_df=True)
    p_value = results.lb_pvalue.values[0]

    if verbose:
        print(f"--- Ljung-Box Test for {name} (Lag {lags}) ---")
        print(f"P-value: {p_value:.4f}")
        
        if p_value < 0.05:
            print(f"VERDICT: The {name} HAS significant autocorrelation.")
            print("Reason: P-value is < 0.05. We reject the Null Hypothesis (White Noise).")
        else:
            print(f"VERDICT: The {name} does NOT have significant autocorrelation.")
            print("Reason: P-value is >= 0.05. The series is consistent with White Noise.")
        print("-" * 40)
    else:
        return p_value
    
def check_stationarity(
        series: Union[pd.Series, np.ndarray], 
        name: str = "Series", 
        verbose: bool = False
    ) -> float:
    """
    Performs the Augmented Dickey-Fuller (ADF) test for stationarity.

    This test checks the Null Hypothesis (H0) that a unit root is present in 
    a univariate time series (indicating it is non-stationary). A p-value 
    below 0.05 suggests the series is stationary.

    Interpretation for your Dissertation:
    -------------------------------------
    - P-value < 0.05: The series is stationary. It reverts to a constant mean.
    - P-value >= 0.05: The series is I(1) or non-stationary. It may exhibit 
      stochastic drift, making long-term forecasting or fixed-bandwidth 
      density estimation unreliable.

    Args:
        series: The time series to test (e.g., daily returns or bandwidth values).
        name: Label for the series in the printed output.
        verbose: If True, prints a detailed diagnostic report.

    Returns:
        float: The p-value of the ADF test.
    """
    # dropna is essential for financial data which often has leading/trailing NaNs
    clean_series = np.asarray(series)
    clean_series = clean_series[~np.isnan(clean_series)]
    
    result = adfuller(clean_series)
    p_value = float(result[1])

    if verbose:
        print(f"--- ADF Stationarity Test for {name} ---")
        print(f"ADF Statistic: {result[0]:.4f}")
        print(f"P-value: {p_value:.4f}")
        
        if p_value < 0.05:
            print(f"VERDICT: The {name} is STATIONARY.")
            print("Reason: P-value < 0.05. We reject H0 (Unit Root present).")
        else:
            print(f"VERDICT: The {name} is NON-STATIONARY (Unit Root).")
            print("Reason: P-value >= 0.05. The series behaves like a Random Walk.")
        print("-" * 45)
    
    return p_value

def select_order_ic(
        data: pd.DataFrame, 
        maxlags: int = 10
        ) -> Dict[str, int]:
    """
    Select lag order using AIC and BIC from statsmodels VAR.select_order.
    Returns dict with keys 'aic', 'bic', 'hqic' (if available).
    """
    model = VAR(data)
    sel = model.select_order(maxlags)
    # statsmodels returns object with attributes aic, bic, hqic that are integers (lags)
    return {'aic': int(sel.aic), 'bic': int(sel.bic), 'hqic': int(sel.hqic)}

def fit_var(data: pd.DataFrame, nlags: int) -> VAR:
    """
    Fit a VAR model and return the fitted results object.
    """
    model = VAR(data, trend='nc')
    res = model.fit(nlags, trend='n')
    return res

def forecast_var(res, steps: int = 1) -> np.ndarray:
    """
    Forecast using a fitted statsmodels VARResults object.
    Returns a numpy array (steps x k).
    """
    return res.forecast(res.endog[-res.k_ar:], steps=steps)

class dynamics_forecaster:
    """
    A wrapper class for Vector Autoregression (VAR) modeling of FPCA scores.
    
    This class facilitates the modeling and forecasting of the dynamics underlying 
    functional data, typically used to model the temporal evolution of density 
    functions via their score components.

    Attributes:
        Y (pd.DataFrame): Time series data where columns are FPCA scores.
        fitted_model (VARResults): The statsmodels VAR results object after fitting.
    """

    def __init__(self, Y: pd.DataFrame):
        """
        Initializes the forecaster with the score matrix.

        Args:
            Y (pd.DataFrame): T x K dataframe of scores (T time points, K components).
        """
        self.Y = Y
        self.fitted_model: Optional[VARResults] = None
    
    def select_order_ic(
        self,
        maxlags: int = 10,
        criteria: str = 'bic',
        prevent_zero_lag: bool = False
    ) -> int:
        """
        Selects the optimal lag order based on Information Criteria.

        Args:
            maxlags (int): Maximum number of lags to check.
            criteria (str): 'aic', 'bic', or 'hqic'.
            prevent_zero_lag (bool): If True, forces the lag to be at least 1.

        Returns:
            int: The selected number of lags.
        """
        model = VAR(self.Y)
        sel = model.select_order(maxlags)   
        criteria_dict = {
            'aic': int(sel.aic), 
            'bic': int(sel.bic), 
            'hqic': int(sel.hqic)
        }

        selected_n_lags = criteria_dict[criteria] 

        if prevent_zero_lag:
            selected_n_lags = max(1, selected_n_lags)
        
        if selected_n_lags == 0:
            print("Warning: Selected lag is zero. Forecast will default to the mean.")    
        
        return selected_n_lags
    
    def fit_var(self, nlags: int) -> VARResults:
        """
        Fits the VAR model to the scores using the specified lag order.

        Args:
            nlags (int): Number of lags to include in the model.

        Returns:
            VARResults: The fitted statsmodels results object.
        """
        model = VAR(self.Y)
        # Fitting without trend ('n') as scores are typically centered
        res = model.fit(nlags, trend='n')
        self.fitted_model = res
        return res

    def forecast(self, h: int) -> np.ndarray:
        """
        Produces out-of-sample forecasts for the scores.

        Args:
            h (int): The forecast horizon (number of steps ahead).

        Returns:
            np.ndarray: K x h matrix of forecasted scores, transposed for 
                compatibility with dFPC reconstruction.
        """
        if self.fitted_model is None:
            raise ValueError("Model must be fitted via 'fit_var' before forecasting.")
            
        model = self.fitted_model
        # Use the last 'k_ar' observations from the endogenous data to start forecast
        fc_h = model.forecast(model.endog[-model.k_ar:], steps=h)

        return fc_h.T

    def residual_diagnostics(
            self, 
            nlags: int = 10, 
            plot: bool = False
        ) -> Dict[str, Union[float, bool]]:
        """
        Performs multivariate diagnostic tests on the VAR residuals.

        Tests for Whiteness (Portmanteau), Normality (Jarque-Bera), and 
        Mathematical Stability (Unit Root check).

        Args:
            nlags (int): Lags to use for the Portmanteau whiteness test.
            plot (bool): If True, displays ACF and Histogram/KDE plots for each score.

        Returns:
            Dict[str, Union[float, bool]]: P-values for tests and stability status.
        """
        if self.fitted_model is None:
            raise ValueError("Model must be fitted before running diagnostics.")
        
        res = self.fitted_model
        residuals = np.asarray(res.resid) 
        n_scores = residuals.shape[1]
        
        # --- 1. Statistical Tests ---
        white_test = res.test_whiteness(nlags=nlags)
        norm_test = res.test_normality()
        is_stable = res.is_stable(verbose=False)
        
        # --- 2. Descriptive Print ---
        header = "VAR RESIDUAL DIAGNOSTICS"
        print(f"\n{'='*50}\n{header.center(50)}\n{'='*50}")
        
        w_status = "PASSED" if white_test.pvalue > 0.05 else "FAILED"
        print(f"Whiteness P-value: {white_test.pvalue:.4f} -> [{w_status}]")
        print(f"Description: Residual series {'is White Noise' if w_status == 'PASSED' else 'is NOT White Noise'}.")
        
        print("-" * 50)
        
        n_status = "PASSED" if norm_test.pvalue > 0.05 else "FAILED"
        print(f"Normality P-value: {norm_test.pvalue:.4f} -> [{n_status}]")
        print(f"Description: Residual distribution {'is Gaussian' if n_status == 'PASSED' else 'is NOT Gaussian'}.")
        
        print("-" * 50)
        
        s_status = "STABLE" if is_stable else "UNSTABLE"
        print(f"Model Stability:   {s_status}")
        print(f"Description: System {'is stationary and forecastable' if is_stable else 'is non-stationary'}.")
        print("="*50)

        # --- 3. Visualization ---
        if plot:
            fig, axes = plt.subplots(n_scores, 2, figsize=(12, 3 * n_scores))
            
            # Ensure axes is 2D even if n_scores == 1
            if n_scores == 1:
                axes = np.expand_dims(axes, axis=0)

            for i in range(n_scores):
                # Column 0: ACF
                plot_acf(residuals[:, i], ax=axes[i, 0], lags=nlags)
                axes[i, 0].set_title(f"ACF: Score {i+1}")
                
                # Column 1: Distribution (Blue Histogram with Red KDE)
                sns.histplot(
                    residuals[:, i], 
                    ax=axes[i, 1], 
                    color='#8ac926', # Light green bars
                    stat='density',
                    alpha=0.6
                )
                sns.kdeplot(
                    residuals[:, i],
                    ax=axes[i, 1],
                    color='red',
                    lw=2,
                    label='KDE'
                )
                axes[i, 1].set_title(f"Distribution: Score {i+1}")
                
            plt.tight_layout()
            plt.show()
        
        return {
            'whiteness_p': float(white_test.pvalue),
            'normality_p': float(norm_test.pvalue),
            'is_stable': bool(is_stable)
        }

def run_forecaster(
        scores    : np.array,
        maxlags_  : int,
        criteria_ : str,
        h_        : int,
        selected_nlags : int = None
        ):
    
    forecaster = dynamics_forecaster(scores)

    if selected_nlags is None:
        selected_nlags = forecaster.select_order_ic(
                                        maxlags_,
                                        criteria=criteria_)
        
    if selected_nlags == 0:
        mean_vec = scores.mean(axis=0)
        fc = np.tile(mean_vec, (h_, 1)).T

    else:
        forecaster.fit_var(nlags=selected_nlags)
        fc = forecaster.forecast(h=h_)

    return fc

def train_test_split(
    Y: pd.DataFrame,
    Y_support : pd.DataFrame,
    train_size: Union[float, int]
    ):
    """
    Parameters
    ----------
    Y : pd.DataFrame
        Columns correspond to time (e.g. dates).
    train_size : float or int
        - float in (0, 1]: fraction of data used for training
        - int >= 1: number of observations in the test set
    """

    n_obs = Y.shape[1]

    if isinstance(train_size, float):
        if not (0 < train_size < 1):
            raise ValueError("train_size as float must be in (0, 1)")
        train_index = int(n_obs * train_size)

    elif isinstance(train_size, int):
        if train_size <= 0 or train_size >= n_obs:
            raise ValueError("train_size as int must be between 1 and n_obs-1")
        train_index = n_obs - train_size

    else:
        raise TypeError("train_size must be float or int")

    Y_train = Y.iloc[:, :train_index]
    Y_train_support = Y_support.iloc[:, :train_index]
    Y_test = Y.iloc[:, train_index:]
    Y_test_support = Y_support.iloc[:, train_index:]

    return Y_train_support, Y_train, Y_test_support, Y_test

# CROSS VALIDATION
def expanding_window_cv(
    T: int, 
    h: int = 1, 
    initial_window: int = 50
) -> List[Tuple[np.ndarray, np.ndarray]]:
    """
    Generate expanding window cross-validation splits for time series data.
    
    This function creates training/testing index pairs for time series cross-validation
    using an expanding window approach, where the training set grows while maintaining
    a fixed forecast horizon.
    
    Parameters
    ----------
    T : int
        Total number of time series observations (length of the dataset).
    h : int, default=1
        Forecast horizon (number of steps to predict ahead).
    initial_window : int, default=50
        Minimum number of observations required for the initial training window.
        Must be less than T - h.
    
    Returns
    -------
    List[Tuple[np.ndarray, np.ndarray]]
        List of (train_indices, test_indices) tuples, where each tuple contains:
        - train_indices : np.ndarray
            Array of indices for the training set (0 to t-1)
        - test_indices : np.ndarray
            Array of indices for the test set (t to t+h-1)
    """

    splits = []
    for t in range(initial_window, T - h + 1):
        train_idx = np.arange(0, t)
        test_idx = np.arange(t, t + h)
        splits.append((train_idx, test_idx))

    return splits

def slice_fold(Y, Y_support, train_idx, test_idx):
    return {
        "Y_train": Y.iloc[:, train_idx],
        "Y_support_train": Y_support.iloc[:, train_idx],
        "Y_test":  Y.iloc[:, test_idx],
        "Y_support_test":  Y_support.iloc[:, test_idx],
    }


############################################################################################
##################################### ACCURACY METRICS #####################################
############################################################################################
 
def KLdiv_matrix(object, eps=1e-4, overlap=True):
    """
    Kullback–Leibler divergence between columns of a matrix.
    Adapted from R package flexmix and function KLdiv.

    Parameters
    ----------
    object : np.ndarray, shape (n, p)
        Each column is a discrete density / pmf.
    eps : float
        Small positive constant to avoid log(0).
    overlap : bool
        If False, requires positive overlap between densities.

    Returns
    -------
    z : np.ndarray, shape (p, p)
        Pairwise KL divergence matrix.
    """
    if not np.issubdtype(object.dtype, np.number):
        raise ValueError("object must be a numeric matrix")

    object = object.copy()
    p = object.shape[1]

    # initialize result
    z = np.full((p, p), np.nan)

    # replace small values
    w = object < eps
    if np.any(w):
        object[w] = eps

    # normalize columns (like sweep(..., 2, colSums, "/"))
    object = object / object.sum(axis=0, keepdims=True)

    for k in range(p - 1):
        for l in range(1, p):
            ok = (object[:, k] > eps) & (object[:, l] > eps)
            if (not overlap) or np.any(ok):
                z[k, l] = np.sum(
                    object[:, k] *
                    (np.log(object[:, k]) - np.log(object[:, l]))
                )
                z[l, k] = np.sum(
                    object[:, l] *
                    (np.log(object[:, l]) - np.log(object[:, k]))
                )

    np.fill_diagonal(z, 0.0)
    return z

def KLdiv(
        p : np.array, 
        q : np.array, 
        eps=1e-4, 
        overlap=True) -> np.float64:
    """
    Kullback–Leibler divergence KL(p || q) for two discrete densities.
    Adapted from R package flexmix and function KLdiv.

    Parameters
    ----------
    p, q : array-like, shape (n,)
        Discrete densities defined on the same grid.
    eps : float
        Small positive constant to avoid log(0).
    overlap : bool
        If False, requires positive overlap between p and q
        (same logic as the original R code).

    Returns
    -------
    float
        KL(p || q)

    Example
    -------
    x = np.linspace(-3, 3, 200)
    from scipy.stats import norm, uniform, t
    y = pd.DataFrame(np.array([uniform.pdf(x), norm.pdf(x), t.pdf(x, df=10)]).T, columns=["Uniform", "Normal", "t"])
    y.plot()
    KLdiv(y.Uniform, y.Normal)
    KLdiv(y.t, y.Normal) # more similar than Uniform and Normal
    """
    p = np.asarray(p, dtype=float).copy()
    q = np.asarray(q, dtype=float).copy()

    if p.ndim != 1 or q.ndim != 1:
        raise ValueError("p and q must be one-dimensional arrays")
    if p.shape[0] != q.shape[0]:
        raise ValueError("p and q must have the same length")

    # replace small values
    p[p < eps] = eps
    q[q < eps] = eps

    # renormalize (same as sweep + colSums)
    p /= p.sum()
    q /= q.sum()

    # overlap condition (faithful to the R code)
    ok = (p > eps) & (q > eps)
    if overlap and not np.any(ok):
        return np.nan
    
    kld = np.sum(p * (np.log(p) - np.log(q)))

    return kld

def JSdiv(
        p: np.array,
        q: np.array,
        eps: float = 1e-4,
        overlap: bool = True
    ) -> np.float64:
    """
    Jensen–Shannon divergence between two discrete probability densities.
    Adapted from Kokoszka 2019.

    This function computes the (symmetric) Jensen–Shannon divergence (JSD)
    between two discrete densities `p` and `q` defined on a common grid.
    The JSD is defined as

        JSD(p, q) = 0.5 * KL(p || m) + 0.5 * KL(q || m),

    where

        m = 0.5 * (p + q)

    is the arithmetic mean mixture of the two densities, and `KL(·||·)`
    denotes the Kullback–Leibler divergence.

    The implementation follows standard practice in density forecast
    evaluation and information theory, using thresholding and
    renormalization to ensure numerical stability.

    Parameters
    ----------
    p, q : array-like of shape (n,)
        Discrete probability densities evaluated on the same support/grid.
        The inputs need not be perfectly normalized; normalization is
        enforced internally.

    eps : float, optional (default=1e-4)
        Small positive constant used to threshold the densities in order
        to avoid undefined logarithms. Values smaller than `eps` are
        replaced by `eps` prior to normalization.

    overlap : bool, optional (default=True)
        If True, the overlap condition is enforced inside the KL divergence
        computations, mimicking the behavior of the original R
        implementation from the `flexmix` package. If the overlap condition
        is violated, the function returns `np.nan`.

    Returns
    -------
    jsd : np.float64
        Jensen–Shannon divergence between `p` and `q`. The value is
        non-negative and finite. Using natural logarithms, the divergence
        is bounded above by `log(2)`.

    Notes
    -----
    - The Jensen–Shannon divergence is symmetric, i.e.,
      JSD(p, q) = JSD(q, p).
    - Unlike the Kullback–Leibler divergence, JSD is always finite when
      computed with the arithmetic mean mixture.
    - The square root of JSD defines a metric on the space of probability
      distributions.
    - This implementation relies on the `KLdiv` function and inherits its
      numerical and overlap-handling conventions.

    See Also
    --------
    KLdiv : Kullback–Leibler divergence for discrete densities.

    Examples
    --------
    x = np.linspace(-3, 3, 200)
    from scipy.stats import norm, uniform, t
    y = pd.DataFrame(np.array([uniform.pdf(x), norm.pdf(x), t.pdf(x, df=10)]).T, columns=["Uniform", "Normal", "t"])
    y.plot()
    JSdiv(y.Uniform, y.Normal) 
    JSdiv(y.t, y.Normal) # more similar than Uniform and Normal
    """
    p = np.asarray(p, dtype=float).copy()
    q = np.asarray(q, dtype=float).copy()

    # Threshold to avoid log(0)
    p[p < eps] = eps
    q[q < eps] = eps

    # Renormalize to ensure valid probability densities
    p /= p.sum()
    q /= q.sum()

    # Arithmetic mean mixture
    m = 0.5 * (p + q)

    jsd = 0.5 * KLdiv(p, m, eps=eps, overlap=overlap) \
        + 0.5 * KLdiv(q, m, eps=eps, overlap=overlap)

    return jsd

def L_norm(
        p: np.array,
        q: np.array,
        norm_type: str = "L1"
    ) -> np.float64:
    """
    Discrete Lp-type norm between two functions evaluated on a common grid.

    Parameters
    ----------
    p, q : array-like of shape (n,)
        Discrete evaluations of two functions (e.g., densities) on the
        same support/grid.

    norm_type : {"L1", "L2", "L_inf"}, optional (default="L1")
        Type of norm to compute:
        - "L1"    : sum_i |p_i - q_i|
        - "L2"    : sqrt( sum_i (p_i - q_i)^2 )
        - "LINF" : max_i |p_i - q_i|

    Returns
    -------
    norm : np.float64
        Value of the selected discrete norm.

    Notes
    -----
    - These are discrete norms; no integration weights are applied.
    - The norms treat densities as generic functions and do not exploit
      their probabilistic structure.
    """
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)

    if p.shape != q.shape:
        raise ValueError("p and q must have the same shape")

    if norm_type == "L1":
        norm = np.sum(np.abs(p - q))

    elif norm_type == "L2":
        norm = np.sqrt(np.sum((p - q) ** 2))

    elif norm_type == "LINF":
        norm = np.max(np.abs(p - q))

    else:
        raise ValueError("norm_type must be one of {'L1', 'L2', 'LINF'}")
    
    return norm
    

def get_metrics(
        p: np.array,
        q: np.array
    ) -> Dict:
    """
    Compute a collection of functional and divergence-based error measures
    between two discrete densities.

    This function evaluates multiple discrepancy measures between two
    density-valued functions `p` and `q` defined on a common grid. The
    returned metrics include both information-theoretic divergences and
    classical functional norms, allowing for a comprehensive assessment
    of density forecast accuracy.

    The following measures are computed:
    - Kullback–Leibler divergence (KLD)
    - Jensen–Shannon divergence (JSD)
    - L1 norm (integrated absolute error, discrete)
    - L2 norm (root integrated squared error, discrete)
    - L-infinity norm (supremum norm, discrete)

    Parameters
    ----------
    p, q : array-like of shape (n,)
        Discrete probability densities evaluated on the same support/grid.
        The inputs need not be perfectly normalized; normalization and
        numerical stabilization are handled internally by the divergence
        functions.

    Returns
    -------
    metrics_dict : dict
        Dictionary containing the following key–value pairs:

        - "KLD"       : float
            Kullback–Leibler divergence KL(p || q).
        - "JSD"       : float
            Jensen–Shannon divergence between p and q.
        - "L1_norm"   : float
            Discrete L1 norm ||p − q||₁.
        - "L2_norm"   : float
            Discrete L2 norm ||p − q||₂.
        - "LINF_norm" : float
            Discrete L-infinity norm ||p − q||∞.

    Notes
    -----
    - KLD and JSD explicitly exploit the probabilistic structure of
      densities and quantify distributional discrepancy and information
      loss.
    - L1, L2, and L-infinity norms treat densities as generic functions
      and measure pointwise approximation error on the grid.
    - All norms are discrete and depend on the grid resolution; no
      integration weights are applied.
    - This function is intended for forecast evaluation in functional
      time series of densities, where complementary error measures
      provide a more complete comparison than any single metric.

    See Also
    --------
    KLdiv : Kullback–Leibler divergence for discrete densities.
    JSdiv : Jensen–Shannon divergence.
    L_norm : Dispatcher for discrete Lp-type norms.

    Examples
    --------
    x = np.linspace(-3, 3, 200)
    from scipy.stats import norm, uniform, t
    y = pd.DataFrame(np.array([uniform.pdf(x), norm.pdf(x), t.pdf(x, df=10)]).T, columns=["Uniform", "Normal", "t"])
    metrics = get_metrics(y.Uniform, y.Normal)
    metrics["JSD"]
    metrics["L2_norm"]
    """
    KLD  = KLdiv(p, q)
    JSD  = JSdiv(p, q)
    L1   = L_norm(p, q, norm_type="L1")
    L2   = L_norm(p, q, norm_type="L2")
    LINF = L_norm(p, q, norm_type="LINF")

    metrics_dict = {
        "KLD": KLD,
        "JSD": JSD,
        "L1_norm":   L1,
        "L2_norm":   L2,
        "LINF_norm": LINF
    }

    return metrics_dict

def overall_KLD(
        X : pd.DataFrame,
        Y : pd.DataFrame,
        n_test : int
    ):
    """
    Compute the overall symmetric Kullback–Leibler divergence between
    two sequences of discretized density functions.

    For each evaluation index, this function computes the Kullback–Leibler
    divergence between the observed density and its forecast in both
    directions, KL(f || f̂) and KL(f̂ || f). The resulting divergences are
    then averaged over the evaluation sample and summed, yielding the
    time-averaged Jeffreys divergence as an overall measure of density
    forecast accuracy.

    Parameters
    ----------
    X : pandas.DataFrame
        Observed densities evaluated on a common grid. Each column
        corresponds to a time index (or forecast origin), and each
        column is assumed to represent a discrete probability density
        function.

    Y : pandas.DataFrame
        Forecast densities evaluated on the same grids as `X`.
        Must have the same shape and column ordering as `X`.

    n_test : int
        Number of evaluation points (columns) over which the divergence
        is computed.

    Returns
    -------
    summary : float
        Overall density forecast error measured as the average symmetric
        Kullback–Leibler divergence (Jeffreys divergence) across the
        evaluation period.

    Notes
    -----
    * Divergences are computed independently for each column (time index).
    * The function `KLdiv` is assumed to compute KL(p || q) for two
      discrete densities p and q.
    * The final summary is given by
        (1 / n_test) * sum_t [ KL(f_t || f̂_t) + KL(f̂_t || f_t) ],
      where f_t and f̂_t denote the observed and forecast densities at
      time t.
    * The resulting measure is symmetric but is not itself a proper
      scoring rule.
    """
    
    measures = np.full((n_test, 2), np.nan)
    for p in range(n_test):
        p_true = X.iloc[:, p]
        p_fc = Y.iloc[:, p]
        measures[p, 0] = KLdiv(p_true, p_fc)
        measures[p, 1] = KLdiv(p_fc, p_true)

    summary = np.round(np.nanmean(measures, axis=0).sum(), 6)

    return summary

def overall_JSD(
        X: pd.DataFrame,
        Y: pd.DataFrame,
        n_test: int
    ) -> float:
    """
    Compute the overall Jensen–Shannon divergence between two sequences
    of discretized density functions using the arithmetic mean mixture.

    For each evaluation index, the Jensen–Shannon divergence between the
    observed and forecast densities is computed. The resulting values
    are then summed over the evaluation sample to obtain an overall
    measure of density forecast accuracy.

    Parameters
    ----------
    X : pandas.DataFrame
        Observed densities evaluated on a common grid. Each column
        corresponds to a time index.

    Y : pandas.DataFrame
        Forecast densities evaluated on the same grid as `X`.
        Must have the same shape and column ordering as `X`.

    n_test : int
        Number of evaluation points (columns).

    Returns
    -------
    summary : float
        Overall Jensen–Shannon divergence summed over the evaluation
        period.

    Notes
    -----
    * The Jensen–Shannon divergence is computed as
          JSD(f_t, f̂_t) =
          0.5 * KL(f_t || m_t) + 0.5 * KL(f̂_t || m_t),
      where m_t = 0.5 * (f_t + f̂_t).
    * Divergences are computed independently for each time index and
      aggregated by summation, matching the original R implementation.
    * The resulting measure is symmetric, finite, and bounded above
      by log(2) when natural logarithms are used.

    See Also
    --------
    KLdiv : Kullback–Leibler divergence for discrete densities.
    JSdiv : Jensen–Shannon divergence for discrete densities.
    """

    jsd = np.full(n_test, np.nan)

    for t in range(n_test):
        p_true = X.iloc[:, t].to_numpy()
        p_fc   = Y.iloc[:, t].to_numpy()

        jsd[t] = JSdiv(p_true, p_fc)

    summary = np.round(np.nansum(jsd), 6)

    return summary

def overall_Lnorm(
        X: pd.DataFrame,
        Y: pd.DataFrame,
        n_test: int,
        norm: str = "L1"
    ) -> float:
    """
    Compute the overall L^p-type norm between observed and forecast
    density functions over an evaluation sample.

    Parameters
    ----------
    X : pandas.DataFrame
        Observed densities evaluated on a common grid.
        Each column corresponds to a time index.

    Y : pandas.DataFrame
        Forecast densities evaluated on the same grid as `X`.

    n_test : int
        Number of evaluation points (columns).

    norm : {"L1", "L2", "LINF"}, optional
        Type of norm to compute:
        - "L1"   : sum of absolute deviations
        - "L2"   : square root of sum of squared deviations
        - "LINF" : maximum absolute deviation

    Returns
    -------
    summary : float
        Mean norm value over the evaluation period, rounded to 4 decimals.

    Notes
    -----
    - Norms are computed pointwise for each time index and then averaged
      over the evaluation sample.
    - This function reproduces the aggregation used in the original R
      implementation exactly.
    """

    values = np.full(n_test, np.nan)

    for t in range(n_test):
        diff = X.iloc[:, t].to_numpy() - Y.iloc[:, t].to_numpy()

        if norm.upper() == "L1":
            values[t] = np.sum(np.abs(diff))

        elif norm.upper() == "L2":
            values[t] = np.sqrt(np.sum(diff ** 2))

        elif norm.upper() == "LINF":
            values[t] = np.max(np.abs(diff))

        else:
            raise ValueError("norm must be one of {'L1', 'L2', 'LINF'}")

    return np.round(np.nanmean(values), 6)

def overall_measures(
    test: pd.DataFrame,
    forecast: pd.DataFrame,
) -> Dict:

    n = forecast.shape[1]

    measures = {
        "KLD": overall_KLD(test, forecast, n),
        "JSD": overall_JSD(test, forecast, n),
        "L_1": overall_Lnorm(test, forecast, n, "L1"),
        "L_2": overall_Lnorm(test, forecast, n, "L2"),
        "L_INFTY": overall_Lnorm(test, forecast, n, "LINF"),
    }

    return measures

#------------------------------------------------
#------------ CROSS-VALIDATION ------------------
#------------------------------------------------

def cv(
        Y, 
        Y_support, 
        KdFPC_kwargs, 
        horizon=1, 
        initial_window=100,
        return_curves=False,
        var_lags: int = None
        ):
    windows = expanding_window_cv(Y.shape[1], h=horizon, initial_window=initial_window)

    measures = []
    for fold, window in enumerate(windows):
        print(f"\t\t>>> cv {fold+1}/{len(windows)}")
        idx_train = window[0]
        idx_test  = window[1]
        
        # TRAIN-TEST SPLIT FOR DENSITIES AND SUPPORTS
        Y_train_support, Y_train = Y_support.iloc[:,idx_train], Y.iloc[:,idx_train]
        Y_test_support , Y_test = Y_support.iloc[:,idx_test],  Y.iloc[:,idx_test]

        # DENSITIES TO LQDENSITIES
        bovespa_mLQDT = mLQDT(
                            Y_train,
                            Y_train_support
                        )
        bovespa_mLQDT.densities_to_lqdensities(verbose=False)

        # L2 EXPANSION
        df_lqds = bovespa_mLQDT.lqd.copy()
        df_lqds_support = bovespa_mLQDT.lqd_support.copy()
        df_lqds_support_du = df_lqds_support[1] - df_lqds_support[0]
        # KdFPC_kwargs = {
        #     "lag_max": 5,
        #     "alpha": 0.10,
        #     "du": 0.05,
        #     "B": 1000,
        #     "p": 5,
        #     "u": df_lqds_support,
        #     "select_ncomp": False,
        #     "dimension": 2
        # }
        KdFPC_kwargs.update(
            {
                "u": df_lqds_support,
                "du": df_lqds_support_du
            }
        )

        KdFPC_model = K_dFPC(df_lqds.values)
        KdFPC_model.fit(**KdFPC_kwargs)
        k_scores = KdFPC_model.etahat.real.T    

        # FORECASTING
        maxlags_  = 10
        criteria_ = 'bic'
        ## SCORES
        k_etahat_fc = run_forecaster(k_scores, maxlags_, criteria_, horizon, selected_nlags=var_lags)
        ## c=F(0)
        model = pm.auto_arima(
            bovespa_mLQDT.c,                         # univariate series
            seasonal=False,            # True if SARIMA
            error_action='ignore',     # ignore non-invertible models
            suppress_warnings=True,
            information_criterion='bic',
            trace=False
        )
        c_forecast, conf_int = model.predict(n_periods=horizon, return_conf_int=True)
        ## RECONSTRUCT FORECASTED CURVES
        k_curve_forecast = KdFPC_model.predict(k_etahat_fc)
        df_k_forecast = pd.DataFrame(k_curve_forecast, columns=Y_test.columns)

        # LQDENSITIES TO DENSITIES
        kle_bkw_supports, kle_bkw_densities = obtain_densities_from_lqd(
                                                                    df_k_forecast,
                                                                    bovespa_mLQDT.lqd_support,
                                                                    c_forecast,
                                                                    t_=bovespa_mLQDT.t,
                                                                    verbose=False
                                                                    )
        
        # PUT KDE AND FORECAST INTO THE SAME GRID FOR EVALUATION
        df_supp, df_f_kle, df_kle_fhat = align_densities(
                                        Y_test_support, 
                                        Y_test, 
                                        kle_bkw_supports, 
                                        kle_bkw_densities, 
                                        kle_bkw_densities.columns
                                        )

        # COMPUTE ACCURACY MEASURES AND STORE THEM
        oa_measures = overall_measures(test=df_f_kle, forecast=df_kle_fhat)
        d1 = {
            "fold": fold,
            "method": "KLE",
            }
        d1.update(oa_measures)
        if return_curves:
            d1.update({"df_supports": df_supp, "df_kde": df_f_kle, "df_forecast": df_kle_fhat})
        measures.append(d1)

    return measures

#------------------------------------------------
#------------ STATIONARITY --- ------------------
#------------------------------------------------

def simulate_stationary_curves(
        n_days=250, 
        n_grid=5001, 
        scenario="stationary"
    ):
    """
    Generates synthetic functional data (densities) to validate stationarity tests.

    This utility creates a time series of curves discretization over a spatial grid, 
    allowing for the simulation of three distinct data-generating processes (DGP) 
    common in financial functional data analysis.

    Scenarios:
    ----------
    1. "stationary":
       The curves oscillate around a fixed global mean density with stable variance.
       Used to test for Type I Error (ensuring the test doesn't find a break 
       where none exists).

    2. "break":
       Simulates a structural regime shift at t=125 (k*). The distribution jumps 
       from a low-volatility state (sigma=1.0) to a high-volatility state (sigma=1.5).
       Used to test for Test Power (ensuring the test correctly identifies 
       structural changes).

    3. "unit_root":
       Simulates a Functional Random Walk, or I(1) process. Each day's curve 
       is the previous day's curve plus a small functional innovation.
       Used to observe the behavior of CUSUM statistics under non-stationary 
       stochastic drift, where k* becomes a spurious artifact of the wander.

    Args:
        n_days: Number of functional observations (time steps/days). Defaults to 250.
        n_grid: Number of discretization points for each curve. Defaults to 5001.
        scenario: The DGP to simulate ("stationary", "break", or "unit_root").

    Returns:
        Tuple[pd.DataFrame, np.ndarray]:
            - A DataFrame of shape (n_grid, n_days) containing the curves.
            - A NumPy array containing the spatial grid (x-axis) values.

    Examples
    ----------
        >>>df_stat, grid = generate_test_densities(scenario="stationary")
        >>>df_break, _ = generate_test_densities(scenario="break")
        >>>df_unitroot, _ = generate_test_densities(scenario="unit_root")
    """
    grid = np.linspace(-5, 5, n_grid)
    densities = np.zeros((n_grid, n_days))
    
    np.random.seed(42)
    
    if scenario == "stationary":
        # Mean and Vol remain constant, only small random noise
        for t in range(n_days):
            mu = 0 + np.random.normal(0, 0.02)
            sigma = 1 + np.random.normal(0, 0.02)
            densities[:, t] = norm.pdf(grid, mu, sigma)
            
    elif scenario == "break":
        # A structural break occurs at Day 125 (k*)
        # Shift from low volatility to high volatility
        for t in range(n_days):
            if t < 125:
                mu, sigma = 0, 1.0
            else:
                mu, sigma = 0, 1.5 # Volatility jump
            
            # Add a tiny bit of noise so it's not a perfect step function
            mu += np.random.normal(0, 0.01)
            sigma += np.random.normal(0, 0.01)
            densities[:, t] = norm.pdf(grid, mu, sigma)
            
    elif scenario == "unit_root":
            # Functional I(1): Each day is the previous day + functional noise
            current_curve = norm.pdf(grid, 0, 1)
            for t in range(n_days):
                # Innovation: small random shifts in mean and scale
                innov_mu = np.random.normal(0, 0.05)
                innov_sigma = 1 + np.random.normal(0, 0.05)
                innovation = norm.pdf(grid, innov_mu, innov_sigma) * 0.1
                
                current_curve = current_curve + innovation
                densities[:, t] = current_curve

    return pd.DataFrame(densities), grid

class FunctionalStationarityTest:
    """
    Implements the Horváth, Kokoszka, and Rice (2014) test for stationarity 
    in functional time series.

    This class provides tools to test whether a sequence of density functions 
    (or any functional data) is stationary or contains a structural break. 
    It uses a pivotal CUSUM-based approach and functional principal component 
    analysis (FPCA) for dimension reduction.

    Attributes:
    ---------------
        sample (np.ndarray): The J x N data matrix.
        J (int): Number of grid points (discretization level).
        N (int): Number of functional observations (time points).
        grid_vals (np.ndarray): The spatial grid for the functional data.
        cumulative_var (float): Variance threshold for choosing the number of FPCs.
        results (Dict): Results from the run_test method.
    """

    def __init__(
        self, 
        sample: Union[np.ndarray, pd.DataFrame], 
        grid_vals: Optional[np.ndarray] = None, 
        cumulative_var: float = 0.90
    ) -> None:
        """
        Initializes the test with functional data.

        Args:
            sample: A J x N matrix where columns are functions and rows are grid points.
            grid_vals: The x-axis values for the functions. Defaults to [0, 1].
            cumulative_var: Fraction of variance (0.0 to 1.0) to retain in FPCA.
        Example:
        ---------------
            >>>df_break, grid = generate_test_densities(scenario="break")
            >>>tester = FunctionalStationarityTest(df_stat, grid_vals=grid)
            >>>tester.run_test(mc_rep=2000)
            >>>tester.get_summary()
            >>>tester.plot_cusum()
        """
        self.sample: np.ndarray = np.asarray(sample)
        self.J, self.N = self.sample.shape
        self.grid_vals: np.ndarray = (
            grid_vals if grid_vals is not None else np.linspace(0, 1, self.J)
        )
        self.cumulative_var: float = cumulative_var
        
        self.results: Dict[str, Union[float, int]] = {}
        self.mean_density: np.ndarray = np.mean(self.sample, axis=1, keepdims=True)
        self.centered_data: np.ndarray = self.sample - self.mean_density

    def _get_kernel_weight(self, x: float, kernel_type: int = 2) -> float:
        """Computes Bartlett-type kernel weights for long-run covariance."""
        if kernel_type == 1:
            return float(np.maximum(0, np.minimum(1, 1.1 - np.abs(x))))
        return float(np.maximum(0, np.minimum(1, 2 - 2 * np.abs(x))))

    def run_test(
            self, 
            mc_rep: int = 1000, 
            kernel_type: int = 2, 
            h: Optional[float] = None,
            d: Optional[int] = None,
            alpha: float = 0.05
        ) -> Dict[str, Union[float, int]]:
        """
        Executes the pivotal CUSUM-based stationarity test via FPCA and Monte Carlo simulation.

        This method implements the statistical framework of Horváth, Kokoszka, and Rice (2014)
        to test the null hypothesis (H0) that the functional time series is stationary. 
        The test is 'pivotal' because the limit distribution of the test statistic does 
        not depend on unknown parameters, allowing for robust p-value estimation.

        Mathematical Procedure:
        -----------------------
        1. Long-Run Covariance Estimation: 
           Estimates the operator $C$ using a Newey-West type kernel estimator to 
           account for temporal dependence between curves.
        2. Dimension Reduction (FPCA): 
           Projects the J-dimensional functional data into a d-dimensional subspace 
           defined by the eigenfunctions of the long-run covariance operator.
        3. CUSUM Process Calculation: 
           Constructs a CUSUM bridge for the scores of each principal component.
        4. Test Statistic ($\hat{M}_N$): 
           Calculates the squared $L^2$ norm of the CUSUM bridge, weighted by the 
           reciprocal of the eigenvalues.
        5. Monte Carlo P-Value: 
           Simulates the distribution of $\sum_{j=1}^d \lambda_j \int B_j^2(t)dt$ 
           where $B_j$ are independent Brownian bridges to determine significance.

        Args:
            mc_rep (int): Number of Monte Carlo replications. Higher values improve 
                p-value precision but increase computational load. Defaults to 1000.
            kernel_type (int): Type of weight function for the long-run covariance.
                - 1: Flat-top kernel (Ker1).
                - 2: Bartlett-type kernel (Ker2/Trapezoidal). Defaults to 2.
            h (float, optional): The bandwidth (lag truncation parameter) for the 
                long-run covariance estimator. If None, uses the rule of thumb $h = \sqrt{N}$. 
                Controls the balance between bias and variance in temporal smoothing.
            d (int, optional): Number of principal components. If None, uses heuristic.
            alpha (float): Significance level for the test (e.g., 0.05, 0.01). 
                Determines the critical value and rejection decision.

        Returns:
            Dict[str, Union[float, int]]: A dictionary containing:
                - 'p_value': Probability of observing the statistic under H0.
                - 'statistic': The calculated CUSUM-based test statistic.
                - 'd': Number of functional principal components (FPCs) retained 
                  to explain the specified 'cumulative_var'.
                - 'k_star': The estimated time index of the most likely structural 
                  break (location of the CUSUM maximum).

        Notes:
            A rejection (p < 0.05) suggests that the mean or covariance structure of the 
            densities has shifted, indicating either a structural break or an I(1) process.
        """
        N, J = self.N, self.J
        h_val: float = h if h is not None else N**0.5
        X: np.ndarray = self.centered_data

        # 1. Estimate Long-Run Covariance Matrix (Z)
        Z: np.ndarray = (X @ X.T) / N
        for lag in range(1, int(h_val) + 1):
            weight = self._get_kernel_weight(lag / h_val, kernel_type)
            if weight > 0:
                gamma_lag = (X[:, lag:] @ X[:, :-lag].T) / N
                Z += weight * (gamma_lag + gamma_lag.T)

        # 2. Eigen-decomposition
        eigenvalues, eigenvectors = np.linalg.eigh(Z)
        idx = eigenvalues.argsort()[::-1]
        eigenvalues = eigenvalues[idx] / J
        eigenvectors = eigenvectors[:, idx] * np.sqrt(J)

        # 3. Determine 'd'
        if d is not None:
            d_used = d
        else:
            exp_var = np.cumsum(eigenvalues) / np.sum(eigenvalues)
            d_used = int(np.searchsorted(exp_var, self.cumulative_var) + 1)
        
        # 4. Compute Test Statistic (M_N)
        scores = (X.T @ eigenvectors[:, :d_used]).T / J
        cusum_process = np.cumsum(scores, axis=1)
        bridge = (1/np.sqrt(N)) * (
            cusum_process - (np.arange(1, N+1)/N) * cusum_process[:, -1:]
        )
        stat: float = float(np.sum(np.mean(bridge**2, axis=1) / eigenvalues[:d_used]))

        # 5. Monte Carlo Simulation for Null Distribution
        # The limit distribution is the sum of squared L2 norms of d_used Brownian bridges
        k_grid = np.arange(1, 101)
        z = np.random.normal(0, 1, (mc_rep, d_used, 100))
        simulated_stats = np.sum(
            np.sum(z**2 / (np.pi**2 * k_grid**2), axis=2), axis=1
        )
        
        # Calculate P-value and Critical Value for user-defined alpha
        p_val: float = float(np.mean(simulated_stats > stat))
        critical_value: float = float(np.percentile(simulated_stats, 100 * (1 - alpha)))
        reject_h0: bool = stat > critical_value

        # Estimate Break Point (k*)
        norms = np.linalg.norm(np.cumsum(X, axis=1), axis=0)
        k_star: int = int(np.argmax(norms))

        self.results = {
            "p_value": p_val, 
            "statistic": stat, 
            "critical_value": critical_value,
            "reject_h0": reject_h0,
            "alpha": alpha,
            "d": d_used, 
            "k_star": k_star
        }
        return self.results

    def plot_cusum(self, n_lines: int = 10) -> None:
        """
        Visualizes the Functional CUSUM process to diagnose structural instability.

        The plot renders snapshots of the Brownian Bridge-like process over time. 
        Each curve represents the accumulated discrepancy between the densities 
        up to that day and the global mean density of the entire period.

        Interpretation:
        ---------------------------------------
        1. Curve Amplitude: High peaks or deep valleys indicate regions of the 
           return grid where the local distribution significantly deviates from 
           the 250-day average (e.g., higher peaks = lower local volatility).
           
        2. Temporal Convergence: 
           - If the final curves (darker colors) remain far from the zero-axis, 
             it visually confirms a 'Permanent Structural Break'.
           - If the curves 'bloom' outward and then collapse back toward zero, 
             it indicates a 'Temporary Volatility Shock'.
             
        3. Shape of the Bridge: A consistent 'bulge' in the center suggests 
           changes in the mean or variance, while 'wiggles' in the far ends 
           suggest instability in the tail-risk (heavy tails).

        Args:
            n_lines (int): Number of snapshots to plot across the time horizon.
        """
        if not self.results:
            print("Warning: Run .run_test() first for an accurate plot title.")
            
        bridge = (1 / np.sqrt(self.N)) * np.cumsum(self.centered_data, axis=1)
        plot_indices = np.linspace(0, self.N - 1, n_lines, dtype=int)

        plt.figure(figsize=(10, 6))
        colors = sns.color_palette("rocket", n_lines)
        
        for i, day in enumerate(plot_indices):
            plt.plot(self.grid_vals, bridge[:, day], color=colors[i], 
                     label=f"Day {day}", alpha=0.7)
            
        plt.axhline(0, color='black', lw=1, ls='--')
        p_val_str = f"{self.results.get('p_value', 'N/A')}"
        plt.title(f"Functional CUSUM (p-value: {p_val_str})")
        plt.xlabel("Grid (e.g., Financial Returns)")
        plt.ylabel("Cumulative Deviation")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Time Step")
        plt.grid(alpha=0.2)
        plt.tight_layout()
        plt.show()

    def plot_regime_comparison(self) -> None:
        """
        Plots the average functional signal before and after the detected break k*.
        This is the most direct way to visualize a permanent shift in shape.
        """
        if "k_star" not in self.results:
            self.analyze_persistence()
            
        k = self.results['k_star']
        mean_pre = np.mean(self.sample[:, :k], axis=1)
        mean_post = np.mean(self.sample[:, k:], axis=1)
        
        plt.figure(figsize=(10, 6))
        plt.plot(self.grid_vals, mean_pre, label=f"Mean: Day 1 to {k}", color='steelblue', lw=2)
        plt.plot(self.grid_vals, mean_post, label=f"Mean: Day {k} to {self.N}", color='crimson', lw=2, ls='--')
        
        plt.title(f"Regime Shift Visualization")
        plt.xlabel("Grid Values")
        plt.ylabel("Functional Value")
        plt.legend()
        plt.grid(alpha=0.3)
        plt.show()

    def get_summary(self) -> str:
        """
        Returns a formatted summary string of the test results.
        """
        res = self.results
        if not res:
            return "No results found. Please call .run_test() first."
            
        status = "REJECTED" if res['p_value'] < 0.05 else "NOT REJECTED"
        print (
            f"--- Functional Stationarity Analysis ---\n"
            f"Null Hypothesis (H0): Series is stationary\n"
            f"Result: {status} (p = {res['p_value']:.4f})\n"
            f"Components (d): {res['d']} FPCs\n"
            f"Detected Break (k*): Day {res['k_star']}\n"
            f"----------------------------------------"
        )


###############################################
############ Horta & Ziegelman (2018) #########
###############################################

from src.dynamicFPC import super_fun
def reconstruct_curves(etahat_pred, psihat, Ybar, du, enforce_density=False):
    """
    Reconstruct functional curves from predicted FPC scores.

    Parameters
    ----------
    etahat_pred : ndarray (d0 x T_pred)
        Predicted FPC scores.

    psihat : ndarray (m x d0)
        Estimated eigenfunctions.

    Ybar : ndarray (m x 1)
        Mean function.

    du : float
        Grid spacing.

    enforce_density : bool
        If True, enforces positivity and renormalizes to integrate to 1.

    Returns
    -------
    Yhat_pred : ndarray (m x T_pred)
        Reconstructed curves.
    """

    # Linear reconstruction
    Yhat_pred = Ybar + psihat @ etahat_pred

    if enforce_density:
        Yhat_fix = Yhat_pred.copy()

        # Enforce positivity
        Yhat_fix[Yhat_fix < 0] = 0

        # Renormalize each curve
        for t in range(Yhat_fix.shape[1]):
            integral = np.sum(Yhat_fix[:, t]) * du
            if integral > 0:
                Yhat_fix[:, t] /= integral

        return Yhat_fix

    return Yhat_pred

def cv_dfpc_horta_zieg(
        Y, 
        Y_support, 
        KdFPC_kwargs, 
        horizon=1, 
        initial_window=100,
        return_curves=False,
        var_lags: int = None
        ):
    windows = expanding_window_cv(Y.shape[1], h=horizon, initial_window=initial_window)

    measures = []
    for fold, window in enumerate(windows):
        print(f"\t\t>>> cv {fold+1}/{len(windows)}")
        idx_train = window[0]
        idx_test  = window[1]
        
        # TRAIN-TEST SPLIT FOR DENSITIES AND SUPPORTS
        Y_train_support, Y_train = Y_support.iloc[:,idx_train], Y.iloc[:,idx_train]
        Y_test_support , Y_test = Y_support.iloc[:,idx_test],   Y.iloc[:,idx_test]

        sp = super_fun(
            Y=Y_train.values,
            lag_max=KdFPC_kwargs["lag_max"],
            B=KdFPC_kwargs["B"],
            p=KdFPC_kwargs["p"],
            m=Y_train_support.shape[0],
            du=0.05,
            dimension=KdFPC_kwargs["dimension"],
            alpha=0.05,
            u=Y_train_support.iloc[:,0]
        ) 
        k_scores = pd.DataFrame(sp["etahat"].T)

        # FORECASTING
        maxlags_  = 10
        criteria_ = 'bic'
        ## SCORES
        k_etahat_fc = run_forecaster(k_scores, maxlags_, criteria_, horizon, selected_nlags=3)

        ## RECONSTRUCT FORECASTED CURVES
        # k_curve_forecast = KdFPC_model.predict(k_etahat_fc)
        k_curve_forecast = reconstruct_curves(k_etahat_fc, sp["psihat"], sp["Ybar"], du=0.05)
        df_k_forecast = pd.DataFrame(k_curve_forecast, columns=Y_test.columns)
        
        # PUT KDE AND FORECAST INTO THE SAME GRID FOR EVALUATION
        df_supp, df_f_kle, df_kle_fhat = align_densities(
                                        Y_test_support, 
                                        Y_test, 
                                        Y_test_support, 
                                        df_k_forecast, 
                                        Y_test_support.columns
                                        )

        # COMPUTE ACCURACY MEASURES AND STORE THEM
        oa_measures = overall_measures(test=df_f_kle, forecast=df_kle_fhat)
        d1 = {
            "fold": fold,
            "method": "KLE",
            }
        d1.update(oa_measures)
        d1.update({"df_supports": df_supp, "df_kde": df_f_kle, "df_forecast": df_kle_fhat})
        measures.append(d1)
    
    return pd.DataFrame(measures)