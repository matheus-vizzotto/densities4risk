import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.integrate import simpson
import matplotlib.pyplot as plt
from statsmodels.nonparametric.kernel_regression import KernelReg
from sklearn.preprocessing import PolynomialFeatures
import statsmodels.api as sm
from statsmodels.regression.quantile_regression import QuantReg


def simulate_rv_series(
        n_days=1000, 
        seed=42
        ) -> pd.Series:
    """
    Simulates a Realized Volatility (RV) time series using a latent AR(1) Log-Volatility process.
    
    The simulation mimics the stylized facts of financial volatility, including high 
    persistence, positivity, and occasional jumps. This provides a robust dataset for 
    testing Median Regression and forecasting models in heavy-tailed scenarios.

    Parameters
    ----------
    n_days : int, default 1000
        The number of business days to simulate.
    seed : int, default 42
        The random seed for reproducibility of the Monte Carlo simulation.

    Returns
    -------
    pd.Series
        A pandas Series containing the simulated Realized Volatility levels, 
        indexed by business day dates starting from 2020-01-01.

    Notes
    -----
    The process is modeled as:
        log(RV_t) = mu + phi * (log(RV_{t-1}) - mu) + epsilon_t + Jump_t
    where epsilon_t ~ N(0, sigma_v) and Jump_t is a Bernoulli-driven outlier.
    """
    np.random.seed(seed)
    
    # Parameters for the Log-RV process
    # Financial RV is highly persistent (phi close to 1)
    phi = 0.95 
    mu = -9.0    # Long-run mean of log-volatility
    sigma_v = 0.2 # Volatility of volatility (vov)
    
    log_rv = np.zeros(n_days)
    log_rv[0] = mu
    
    # Generate the AR(1) process in log-space
    for t in range(1, n_days):
        # We add some occasional "jumps" (heavy tails) using a t-distribution 
        # or Bernoulli-Gaussian to make Median Regression relevant
        jump = 2.0 if np.random.rand() > 0.98 else 0.0
        
        log_rv[t] = mu + phi * (log_rv[t-1] - mu) + np.random.normal(0, sigma_v) + jump

    # Convert back to levels (Exponential)
    rv_series = np.exp(log_rv)
    
    # Create a time index (Business days)
    dates = pd.date_range(start='2020-01-01', periods=n_days, freq='B')
    
    y = pd.Series(rv_series, index=dates, name='Realized_Volatility')

    return y


def prepare_har_data(rv_series):
    """
    Constructs the feature matrix for the Heterogeneous Autoregressive (HAR) model.
    
    This function transforms a daily Realized Volatility (RV) series into a 
    multi-component dataset representing the daily, weekly (5-day), and 
    monthly (22-day) average volatility components, as proposed by Corsi (2009).

    Parameters
    ----------
    rv_series : pd.Series
        A time series of daily realized volatility levels.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the lagged components ('RV_d', 'RV_w', 'RV_m') 
        and the lead target variable ('Target'), with NaNs removed.

    Notes
    -----
    The HAR model assumes that agents with different time horizons (speculators, 
    institutional investors, and central banks) perceive and react to 
    volatility differently, leading to the cascaded structure of the components.
    """
    df = pd.DataFrame({'RV_d': rv_series})
    # Create the weekly (5-day) and monthly (22-day) averages
    df['RV_w'] = df['RV_d'].rolling(window=5).mean()
    df['RV_m'] = df['RV_d'].rolling(window=22).mean()
    
    # Target is the next day's RV
    df['Target'] = df['RV_d'].shift(-1)
    return df.dropna()

def forecast_strategy_ben1(train_df, next_day_features):
    """
    Implements Strategy Ben1: Log-HAR model with Median Regression.
    
    train_df: DataFrame containing 'Target', 'RV_d', 'RV_w', 'RV_m'
    next_day_features: Series/Dict containing the RV components for the forecast date n
    """
    
    # 1. Fit the Log-HAR using Quantile Regression (q=0.5 is Median)
    # Using np.log directly in the formula string
    formula = 'np.log(Target) ~ np.log(RV_d) + np.log(RV_w) + np.log(RV_m)'
    model = smf.quantreg(formula, data=train_df)
    res = model.fit(q=0.5, max_iter=2000)
    
    # 2. Extract coefficients (ahat_0, ahat_1, ahat_2, ahat_3)
    intercept = res.params['Intercept']
    a1 = res.params['np.log(RV_d)']
    a2 = res.params['np.log(RV_w)']
    a3 = res.params['np.log(RV_m)']
    
    # 3. Calculate the forecast for n + 1 (Equation 22)
    # log_forecast = a0 + a1*log(RV_d_n) + a2*log(RV_w_n) + a3*log(RV_m_n)
    log_pred = (intercept + 
                a1 * np.log(next_day_features['RV_d']) + 
                a2 * np.log(next_day_features['RV_w']) + 
                a3 * np.log(next_day_features['RV_m']))
    
    # 4. Transform back to original scale
    forecast_n_plus_1 = np.exp(log_pred)
    
    return forecast_n_plus_1


def forecast_har_strategy(
        rv_series, 
        h=1, 
        train_size=0.8,
        plot = True
        ):
    """
    Constructs HAR components, splits data, and forecasts using Median Regression.
    
    Parameters:
    - rv_series: pd.Series of Realized Volatility
    - h: Forecast horizon (days ahead)
    - train_size: Fraction of data for initial training
    """
    df = rv_series.to_frame(name='Target')
    
    # 1. Create HAR components (daily, weekly, monthly averages)
    df['RV_d'] = df['Target'].shift(h)
    df['RV_w'] = df['Target'].shift(h).rolling(window=5).mean()
    df['RV_m'] = df['Target'].shift(h).rolling(window=22).mean()
    
    # Drop rows where we don't have enough history for the monthly component
    df = df.dropna()
    
    # 2. Split indices
    n = len(df)
    split_idx = int(n * train_size)
    
    actuals = []
    forecasts = []
    dates = []

    # 3. Expanding Window Walk-Forward
    # We iterate through the 'test' portion of the data
    for i in range(split_idx, n):
        train_sub = df.iloc[:i]
        test_row = df.iloc[i]
        
        # Fit Strategy Ben1 (Median Regression on Logs)
        # Note: np.log() is applied within the formula
        res = smf.quantreg('np.log(Target) ~ np.log(RV_d) + np.log(RV_w) + np.log(RV_m)', 
                           data=train_sub).fit(q=0.5)
        
        # Forecast for date n+h using Equation 22 logic
        # Predict uses the features available at the current index i
        log_pred = res.predict(test_row[['RV_d', 'RV_w', 'RV_m']])
        pred = np.exp(log_pred.values[0])
        
        forecasts.append(pred)
        actuals.append(test_row['Target'])
        dates.append(df.index[i])

    # 4. Results and Plotting
    results = pd.DataFrame({'Actual': actuals, 'Forecast': forecasts}, index=dates)
    
    if plot:
        plt.figure(figsize=(12, 6))
        plt.plot(results['Actual'], label='Actual RV', alpha=0.6, color='gray')
        plt.plot(results['Forecast'], label=f'HAR-Median Forecast (h={h})', color='red', linestyle='--')
        plt.title(f'Strategy Ben1: Realized Volatility Forecasting (h={h})')
        plt.legend()
        plt.yscale('log') # Log scale is often better for visualizing volatility
        plt.show()
    
    return results

class Ben1Model:
    """
    Implements a Log-HAR model estimated via Median Regression (q=0.5).
    
    This class handles the transformation of daily realized volatility into 
        daily, weekly, and monthly components and provides a one-step-ahead 
        forecast using a robust median-based estimator.

    Examples
    --------
    >>> import pandas as pd
    >>> import numpy as np
    >>> # Simulate some RV data
    >>> dates = pd.date_range('2026-01-01', periods=50, freq='B')
    >>> rv_data = pd.Series(np.random.lognormal(-9, 0.2, 50), index=dates)
    >>> 
    >>> # Initialize and fit the model
    >>> forecaster = HARMedianForecaster()
    >>> forecaster.fit(rv_data)
    >>> 
    >>> # Forecast the next day using the tail of the history
    >>> next_rv = forecaster.predict(rv_history=rv_data)
    >>> print(f"Forecasted RV: {next_rv:.6f}")
    """

    def __init__(self):
        self.model_fitted = None
        self.formula = 'np.log(Target) ~ np.log(RV_d) + np.log(RV_w) + np.log(RV_m)'

    def prepare_har_data(self, rv_series):
        """
        Transforms a daily RV series into the HAR component structure.
        """
        df = pd.DataFrame({'RV_d': rv_series})
        df['RV_w'] = df['RV_d'].rolling(window=5).mean()
        df['RV_m'] = df['RV_d'].rolling(window=22).mean()
        df['Target'] = df['RV_d'].shift(-1)
        return df.dropna()

    def fit(self, rv_series):
        """
        Prepares the data and estimates the HAR coefficients.
        
        Parameters
        ----------
        rv_series : pd.Series
            The historical realized volatility levels used for training.
        """
        train_df = self.prepare_har_data(rv_series)
        model = smf.quantreg(self.formula, data=train_df)
        self.model_fitted = model.fit(q=0.5, max_iter=2000)

    def predict(self, rv_history=None, next_day_features=None):
        """
        Generates a one-step-ahead forecast for RV_{n+1}.

        Parameters
        ----------
        rv_history : pd.Series, optional
            A series of past RV values. If provided, the function will 
            extract the most recent day, week, and month averages.
        next_day_features : dict or pd.Series, optional
            Pre-computed HAR components. Used if rv_history is not provided.
        """
        if self.model_fitted is None:
            raise ValueError("Model must be fitted before calling predict.")

        # Logic to determine feature source
        if next_day_features is None:
            if rv_history is None:
                raise ValueError("Must provide either 'rv_history' or 'next_day_features'.")
            
            # Use the tail of the history to compute current HAR components (at time n)
            next_day_features = {
                'RV_d': rv_history.iloc[-1],
                'RV_w': rv_history.iloc[-5:].mean(),
                'RV_m': rv_history.iloc[-22:].mean()
            }

        params = self.model_fitted.params
        log_pred = (params['Intercept'] + 
                    params['np.log(RV_d)'] * np.log(next_day_features['RV_d']) + 
                    params['np.log(RV_w)'] * np.log(next_day_features['RV_w']) + 
                    params['np.log(RV_m)'] * np.log(next_day_features['RV_m']))

        return np.exp(log_pred)
    

############################################################################################
###################################### FcstV1 ##############################################
############################################################################################

def _calculate_density_mean(
        u_grid: np.array, 
        density: np.array
        ) -> float:
    """
    Calculates the first raw moment: mu = integral(u * f(u) du)
    """
    mean = simpson(y=u_grid * density, x=u_grid)

    return mean

def _calculate_density_variance(
        u_grid: np.array, 
        density: np.array
        ) -> float:
    """
    Calculates the second central moment: Var = integral((u - mu)^2 * f(u) du)
    """
    mu = _calculate_density_mean(u_grid, density)
    
    # Calculate (u - mu)^2
    squared_deviation = (u_grid - mu)**2
    
    # Integrate squared deviation weighted by the density
    variance = simpson(y=squared_deviation * density, x=u_grid)
    
    return variance

def calculate_volatility_forecast(
        variance : float, 
        n_next : int
        ) -> float:
    """
    Implements one-step ahead variance forecast
        sigma^2_{n+1|n} = n_{n+1} * Var(f_{n+1|n})
    where n_{n+1} is the (scalar) sample size for day d+1.
    """
    variance_t_next = n_next * variance

    return variance_t_next

from scipy.integrate import simpson

# class Fcst1VModel:
#     """
#     Implements the FcstV1 volatility forecasting model based on density moments.
    
#     This model maps a forecasted functional density curve back to a scalar 
#     volatility estimate by scaling the density's variance by the 
#     intraday sample size.
#     """
    
#     def __init__(self):
#         """
#         Parameters
#         ----------

#         """
#         self.u_grid = None

#     def forecast(self, forecasted_density, n_observations):
#         """
#         Generates the scalar volatility forecast.

#         Parameters
#         ----------
#         forecasted_density : np.array
#             The 1D forecasted density curve (e.g., from FPCA or KDE).
#         n_observations : int
#             The number of intraday observations (n_{n+1}) for the target day.
#         u_grid : np.array
#             The common grid of evaluation points (e.g., log-returns) 
#             used for the functional data.

#         Returns
#         -------
#         float
#             The forecasted realized volatility (sigma^2).
#         """
#         # Extract the variance of the forecasted density
#         f_var = _calculate_density_variance(forecasted_density)
        
#         # Apply the scaling factor
#         sigma_sq_hat = n_observations * f_var
        
#         return sigma_sq_hat
    
def Fcst1VModel(
        u_grid: np.array, 
        forecasted_density: pd.Series, 
        n_observations: int
        ) -> float:
    """
    Maps a forecasted functional density curve to a scalar volatility estimate.
    """
    f_var = _calculate_density_variance(
                            u_grid=u_grid, 
                            density=forecasted_density
                            )
    rv_tp1 = n_observations * f_var # realized volatility forecast at t+1

    return rv_tp1

import numpy as np
from statsmodels.nonparametric.kernel_regression import KernelReg

class Fcst2VForecaster:
    """
    Implements Strategies Fcst2V and Fcst2H using Non-parametric Kernel Regression.
    
    This model estimates the relationship log(sigma^2) = m(eta) + epsilon 
    using a Nadaraya-Watson or Local Linear estimator.
    """

    def __init__(self, reg_type='lc', bw='cv_ls'):
        """
        Parameters
        ----------
        reg_type : str, default 'lc'
            Type of regression: 'lc' for Local Constant (Nadaraya-Watson) 
            or 'll' for Local Linear.
        bw : str or array_like, default 'cv_ls'
            Bandwidth selection method. 'cv_ls' uses least-squares cross-validation.
        """
        self.reg_type = reg_type
        self.bw = bw
        self.model = None

    def fit(self, eta_train, rv_train):
        """
        Fits the non-parametric regression m(eta).

        Parameters
        ----------
        eta_train : np.array or pd.DataFrame
            The predictors derived from functional objects (e.g., FPC scores).
        rv_train : np.array or pd.Series
            The observed realized volatility levels.
        """
        eta_train = np.asarray(eta_train)
        rv_train = np.asarray(rv_train)

        # Ensure 2D shape for exog (KernelReg requirement)
        if eta_train.ndim == 1:
            eta_train = eta_train.reshape(-1, 1)

        # We model the LOG of RV to match the dissertation equation
        # Add epsilon for numerical stability (avoid log(0))
        eps = 1e-10
        log_rv = np.log(rv_train + eps)
        
        # Initialize KernelReg. 
        # 'c' indicates the predictors are continuous. 
        # For multiple predictors, use 'ccc' for three continuous variables.
        var_types = 'c' * eta_train.shape[1]
        
        self.model = KernelReg(
            endog=log_rv,
            exog=eta_train, 
            var_type=var_types, 
            reg_type=self.reg_type, 
            bw=self.bw
        )

    def forecast(self, eta_next):
        """
        Generates the one-step-ahead forecast using the exponential of the 
        non-parametric estimate.

        Example
        -------
        >>> # Simulated training data
        >>> eta_train = np.random.randn(100, 2)  # 2-dimensional scores
        >>> rv_train = np.exp(np.random.randn(100))  # positive RV values
        >>>
        >>> model = Fcst2VForecaster(reg_type='lc', bw='cv_ls')
        >>> model.fit(eta_train, rv_train)
        >>>
        >>> # Forecast using a new eta vector
        >>> eta_next = np.array([0.1, -0.2])
        >>> forecast = model.forecast(eta_next)
        >>> print(forecast)
        """
        if self.model is None:
            raise ValueError("Model must be fitted before forecasting.")

        eta_next = np.asarray(eta_next)

        # Ensure correct shape (1, d)
        if eta_next.ndim == 1:
            eta_next = eta_next.reshape(1, -1)

        # 1. Obtain the non-parametric estimate m_hat(eta_next)
        # Returns (mean_estimate, marginal_effects)
        m_hat_log, _ = self.model.fit(eta_next)
        
        # 2. Apply the exponential transformation to return to sigma^2 levels
        return float(np.exp(m_hat_log[0]))
    
class Fcst3VForecaster:
    """
    Implements Strategies Fcst3V and Fcst3H.
    
    Estimates log(sigma^2) = m(eta) + error, where m is a second-degree 
    polynomial on eta, estimated via Median Regression (q=0.5).
    """

    def __init__(self, degree=2):
        self.degree = degree
        self.poly = PolynomialFeatures(degree=self.degree, include_bias=False)
        self.model_fitted = None
        self.feature_names_ = None

    def _prepare_features(self, eta, fit=False):
        """
        Generates polynomial terms: [eta1, eta2, eta1^2, eta1*eta2, eta2^2]
        """
        eta = np.asarray(eta)

        # Ensure 2D shape
        if eta.ndim == 1:
            eta = eta.reshape(1, -1)

        if fit:
            eta_poly = self.poly.fit_transform(eta)
            self.feature_names_ = [f'feat_{i}' for i in range(eta_poly.shape[1])]
        else:
            eta_poly = self.poly.transform(eta)

        return pd.DataFrame(eta_poly, columns=self.feature_names_)

    def fit(self, eta_train, rv_train):
        """
        Fits the polynomial median regression.
        """
        eta_train = np.asarray(eta_train)
        rv_train = np.asarray(rv_train)

        # 1. Prepare polynomial predictors
        X = self._prepare_features(eta_train, fit=True)
        
        # 2. Prepare target: log(sigma_t^2)
        # Add epsilon for numerical stability
        eps = 1e-10
        X['log_target'] = np.log(rv_train + eps)
        
        # 3. Construct formula: log_target ~ feat_0 + feat_1 + ...
        formula = "log_target ~ " + " + ".join(self.feature_names_)
        
        # 4. Fit Median Regression (Quantile Regression at q=0.5)
        model = smf.quantreg(formula, data=X)
        self.model_fitted = model.fit(q=0.5, max_iter=2000)

    def forecast(self, eta_next):
        """
        Generates the forecast using the exponential of the polynomial estimate.

        Parameters
        ----------
        eta_next : np.array
            The next-period predictor vector (e.g., forecasted FPC scores),
            of shape (d,) or (1, d).

        Returns
        -------
        float
            The one-step-ahead forecast of sigma^2.

        Example
        -------
        >>> eta_train = np.random.randn(100, 2)
        >>> rv_train = np.exp(np.random.randn(100))
        >>>
        >>> model = Fcst3VForecaster(degree=2)
        >>> model.fit(eta_train, rv_train)
        >>>
        >>> eta_next = np.array([0.1, -0.2])
        >>> forecast = model.forecast(eta_next)
        >>> print(forecast)
        """
        if self.model_fitted is None:
            raise ValueError("Model must be fitted before forecasting.")

        # 1. Transform next-day features to polynomial space
        X_next = self._prepare_features(eta_next, fit=False)
        
        # 2. Get prediction in log-space
        log_pred = self.model_fitted.predict(X_next)
        
        # 3. Final forecast: exp(m_hat)
        return float(np.exp(log_pred.iloc[0]))


def qlike_loss(rv, sigma2_hat, eps=1e-10):
    """
    Computes QLIKE loss.

    Parameters
    ----------
    rv : array-like
        Realized variance
    sigma2_hat : array-like
        Forecasted variance
    eps : float
        Small constant for numerical stability

    Returns
    -------
    np.ndarray
        QLIKE loss values
    """
    rv = np.asarray(rv)
    sigma2_hat = np.asarray(sigma2_hat)

    ratio = (rv + eps) / (sigma2_hat + eps)
    return ratio - np.log(ratio) - 1


def mse(x, y):
    return np.mean(np.square(x-y))

class Fcst2VOLS:
    """
    Implements a parametric OLS version of the Fcst2V strategy.
    
    This model estimates the relationship log(sigma^2) = beta * eta + epsilon 
    using Ordinary Least Squares.
    """

    def __init__(self, add_constant=True):
        """
        Parameters
        ----------
        add_constant : bool, default True
            Whether to include an intercept (alpha) in the regression.
        """
        self.add_constant = add_constant
        self.model = None
        self.params = None

    def fit(self, eta_train, rv_train):
        """
        Fits the OLS regression on the log-transformed RV.

        Parameters
        ----------
        eta_train : np.array or pd.DataFrame
            The predictors (e.g., FPC scores).
        rv_train : np.array or pd.Series
            The observed realized volatility levels.
        """
        eta_train = np.asarray(eta_train)
        rv_train = np.asarray(rv_train)

        if eta_train.ndim == 1:
            eta_train = eta_train.reshape(-1, 1)

        # Log transform for the target variable
        log_rv = np.log(rv_train + 1e-10)

        # Handle the constant (intercept)
        X = eta_train
        if self.add_constant:
            X = sm.add_constant(X, has_constant='add')

        # Fit the OLS model
        self.model = sm.OLS(log_rv, X).fit()
        self.params = self.model.params

    def forecast(self, eta_next):
        """
        Generates the forecast using the linear model and exponential transformation.
        """
        if self.model is None:
            raise ValueError("Model must be fitted before forecasting.")

        eta_next = np.asarray(eta_next)

        # Ensure correct shape (1, d)
        if eta_next.ndim == 1:
            eta_next = eta_next.reshape(1, -1)

        # Add constant to prediction input if necessary
        X_next = eta_next
        if self.add_constant:
            X_next = np.hstack([np.ones((X_next.shape[0], 1)), X_next])

        # 1. Predict in log-space: log(sigma^2)_hat = X * beta
        m_hat_log = self.model.predict(X_next)
        
        # 2. Return to sigma^2 levels via exponential
        # If input was a single point, return a float
        return float(np.exp(m_hat_log[0]))

def forecast_reg1(rv_train, scores_train, current_score):
    """
    Implements Strategy Reg1: log(varsigma^2_t) = m(eta_{t-1}) + error
    
    Parameters:
    rv_train: Array of realized variances (varsigma^2) for t=2...n
    scores_train: Array of lagged scores (eta) for t=1...n-1
    current_score: The score at time n (eta_n) to produce forecast for n+1
    """
    
    # 1. Prepare target: log of realized variance
    log_rv = np.log(rv_train)
    
    # 2. Fit the nonparametric regression: m(eta_{t-1})
    # 'c' indicates a continuous variable for the score
    model = KernelReg(endog=log_rv, exog=scores_train, var_type='c')
    
    # 3. Predict m(eta_n)
    # The predict method returns (mean, marginal_effects)
    m_hat_n, _ = model.fit([current_score])
    
    # 4. Apply Strategy Reg1 forecast equation (Eq 24): exp{m_hat(eta_n)}
    forecast_n_plus_1 = np.exp(m_hat_n[0])
    
    return forecast_n_plus_1

# Example Usage with your data structure:
# scores = objects["mdfpc"]["scores_fitted"]["scores_1"]
# rv = rv_train.values

# Ensure alignment: rv[t] is regressed on scores[t-1]
# y = rv[1:] (t=2...n)
# X = scores[:-1] (t=1...n-1)
# current_eta = scores[-1] (t=n)

class Reg1Forecaster(Fcst2VForecaster):
    """
    Implements Strategy Reg1 by regressing current RV on LAGGED scores.
    """
    def fit(self, eta_history, rv_history):
        # Align data: rv[t] depends on eta[t-1]
        # We drop the first RV observation and the last eta observation
        y = rv_history[1:]    # t=2 to n
        X = eta_history[:-1]  # t=1 to n-1
        
        # Use the parent fit method with aligned data
        super().fit(X, y)

    def forecast(self, eta_current):
        # Here, eta_current is eta_n (the most recent observed score)
        # This produces the forecast for n+1 directly.
        return super().forecast(eta_current)
    
class Reg2Forecaster:
    """
    Implements Strategy Reg2: log(varsigma^2_t) = m(eta_{t-1}) + error
    where m is a second-degree polynomial estimated via Median Regression.
    """
    def __init__(self, degree=2):
        self.poly = PolynomialFeatures(degree=degree, include_bias=True)
        self.model_result = None

    def _prepare_features(self, eta):
        eta = np.asarray(eta)
        if eta.ndim == 1:
            eta = eta.reshape(-1, 1)
        # Generates [1, n1, n2, n1^2, n1*n2, n2^2]
        return self.poly.fit_transform(eta)

    def fit(self, eta_history, rv_history):
        # Align: rv[t] depends on eta[t-1]
        y = np.log(rv_history[1:])
        X_raw = eta_history[:-1]
        
        # Transform scores into polynomial features
        X_poly = self._prepare_features(X_raw)
        
        # Fit Median Regression (Quantile = 0.5)
        quant_mod = QuantReg(y, X_poly)
        self.model_result = quant_mod.fit(q=0.5, max_iter=2000)

    def forecast(self, eta_current):
        if self.model_result is None:
            raise ValueError("Model must be fitted before forecasting.")
        
        # Transform current score (eta_n) to polynomial features
        X_next = self._prepare_features(eta_current)
        
        # Predict m_hat(eta_n)
        m_hat_n = self.model_result.predict(X_next)
        
        # Forecast for n+1: exp{m_hat(eta_n)}
        return float(np.exp(m_hat_n[0]))