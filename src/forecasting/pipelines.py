import numpy as np
import pandas as pd
import src.fda.dfpc as dfpc
import src.fda.transformations.lqdt as lqdt
import src.forecasting.models as fm
import pmdarima as pm

def run_multivariate_forecaster(
        scores    : np.array,
        maxlags_  : int,
        criteria_ : str,
        h_        : int,
        selected_nlags : int = None
        ):
    
    forecaster = fm.multivariate_forecaster(scores)

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

class DensityForecaster:
    """
    Time-series forecaster for Probability Density Functions (PDFs).

    This class implements a forecasting pipeline for functional data (densities) 
    by mapping them into a Hilbert space using the Log-Quantile Density (LQD) 
    transformation. It reduces dimensionality via Karhunen-Loève dynamic Functional 
    Principal Component (K-dFPC) Analysis and forecasts the resulting 
    coefficients and LQD constants using VAR and ARIMA models.

    The pipeline follows three main stages:
    1. **Transformation**: Densities are mapped to LQD space to remove 
       non-negativity and integration constraints.
    2. **Reduction**: K-dFPC extracts principal components (eigenfunctions) 
       and scores from the LQD curves.
    3. **Forecasting**: Multivariate forecasting is applied to the scores, 
       while the LQD constant 'c' is modeled via auto-ARIMA.

    Parameters
    ----------
    kdfpc_kwargs : dict, optional
        Keyword arguments passed to the `dfpc.K_dFPC` model (e.g., n_components).
    maxlags : int, default=10
        Maximum number of lags to consider for the multivariate score forecaster.
    criteria : {'bic', 'aic', 'hqic'}, default='bic'
        Information criterion used for ARIMA and VAR model selection.

    Attributes
    ----------
    model_lqd : lqdt.mLQDT
        The fitted LQD transformation object containing LQD curves and constants.
    model_kdfpc : dfpc.K_dFPC
        The fitted Functional PCA model.
    model_arima : pmdarima.ARIMA
        The fitted auto-ARIMA model for the LQD normalization constant 'c'.
    
    Notes
    -----
    The LQD transformation is defined as:
    $$L(f(t)) = \log(f(Q(t)))$$
    where $Q(t)$ is the quantile function. This allows for unconstrained 
    forecasting in $L^2$ space.

    Example
    -------
    >>> from my_module import DensityForecaster
    >>> forecaster = DensityForecaster(maxlags=3, criteria='aic')
    >>> # Y_train: list/array of densities, Y_support: corresponding grids
    >>> forecaster.fit(Y_train, Y_support)
    >>> results = forecaster.predict(horizon=5)
    >>> print(results['future_densities'].shape)
    """
    def __init__(self, 
                 kdfpc_kwargs=None, 
                 maxlags=10, 
                 criteria='bic'
                 ):
        self.kdfpc_kwargs = kdfpc_kwargs or {}
        self.maxlags = maxlags
        self.criteria = criteria
        self.model_lqd = None
        self.model_kdfpc = None
        self.model_arima = None

    def fit(self, 
            Y_train:    pd.DataFrame, 
            Y_support:  pd.DataFrame
            ):
        # 1. Transform Densities
        mlqdt = lqdt.mLQDT()
        self.model_lqd = mlqdt.transform(
            densities=Y_train,
            densities_supports=Y_support, 
            verbose=False
            )
        # self.model_lqd.densities_to_lqdensities(verbose=False)
        
        # 2. L2 Expansion (K_dFPC)
        lqd_values  = self.model_lqd.lqd.values
        lqd_support = self.model_lqd.lqd_support
        
        self.kdfpc_kwargs.update({
            "u": lqd_support,
            "du": self.model_lqd.du
        })
        
        self.model_kdfpc = dfpc.K_dFPC(lqd_values)
        self.model_kdfpc.fit(**self.kdfpc_kwargs)
        
        # 3. Fit ARIMA for the 'c' constant
        self.model_arima = pm.auto_arima(
            self.model_lqd.c, 
            seasonal=False, 
            information_criterion=self.criteria,
            suppress_warnings=True
        )

        return self

    def predict(self, 
                horizon, 
                var_lags=None,
                forecast_index=None
                ):
        # Builds forecasted dates for prediction output
        if forecast_index is not None:
            if len(forecast_index) != horizon:
                raise ValueError(f"forecast_index length ({len(forecast_index)}) "
                                f"must match horizon ({horizon})")
            future_dates = forecast_index
        elif self.index_freq is not None:
            # Automatic generation if freq was caught in .fit()
            future_dates = pd.date_range(start=self.last_date, 
                                        periods=horizon + 1, 
                                        freq=self.index_freq)[1:]
        else:
            # Fallback to integer steps if no dates are available
            future_dates = np.arange(horizon)
        # 1. Forecast Scores
        k_scores = self.model_kdfpc.etahat.values.real
        k_etahat_fc = run_multivariate_forecaster(
            k_scores, self.maxlags, self.criteria, horizon, selected_nlags=var_lags
        )
        
        # 2. Forecast 'c'
        c_forecast = self.model_arima.predict(n_periods=horizon)
        
        # 3. Reconstruct LQD
        L2_curve_forecast = self.model_kdfpc.predict(k_etahat_fc)
        L2_curve_forecast.columns = future_dates
        # print(type(L2_curve_forecast))
        # df_L2_curve_forecast = pd.DataFrame(L2_curve_forecast,future_dates)
        # print(df_L2_curve_forecast)

        fc_lqd_obj = lqdt.LQDRepresentation(
                            lqd         = L2_curve_forecast,
                            c           = c_forecast,
                            lqd_support = self.model_lqd.lqd_support,
                            t0          = self.model_lqd.t0 
                            )
        
        # 4. Back to densities
        mlqdt = lqdt.mLQDT()
        df_supports, df_densities = mlqdt.inverse_transform(fc_lqd_obj)
        
        return {
                "future_L2_curves": L2_curve_forecast,
                "future_scores":    k_etahat_fc,
                "future_cs":        c_forecast,
                "future_supports":  df_supports,
                "future_densities": df_densities
        }

class DensityForecaster_HZ:
    """
    A pipeline for forecasting functional curves using K_dFPC_HZ.
    Matches the DensityForecaster blueprint.
    """
    def __init__(self, kdfpc_kwargs=None, maxlags=10, criteria='AIC'):
        self.kdfpc_kwargs = kdfpc_kwargs or {}
        self.maxlags = maxlags
        self.criteria = criteria
        
        self.model_kdfpc = None
        self.u = None
        self.du = None

    def fit(self, Y_train, df_support):
        # 1. Prepare Support
        self.u = df_support.iloc[:, 0]
        self.du = self.u[1] - self.u[0]
        
        # 2. Update model parameters
        self.kdfpc_kwargs.update({
            "u": self.u,
            "du": self.du
        })
        
        # 3. Fit K_dFPC_HZ
        self.model_kdfpc = dfpc.K_dFPC_HZ(Y_train.values)
        self.model_kdfpc.fit(**self.kdfpc_kwargs)
        
        return self

    def predict(self, 
                horizon, 
                var_lags=2,
                forecast_index=None
                ):
        """
        Generate future curve forecasts. 
        Matches blueprint signature: predict(horizon, var_lags)
        """
        if forecast_index is not None:
            if len(forecast_index) != horizon:
                raise ValueError(f"forecast_index length ({len(forecast_index)}) "
                                f"must match horizon ({horizon})")
            future_dates = forecast_index
        elif self.index_freq is not None:
            # Automatic generation if freq was caught in .fit()
            future_dates = pd.date_range(start=self.last_date, 
                                        periods=horizon + 1, 
                                        freq=self.index_freq)[1:]
        else:
            # Fallback to integer steps if no dates are available
            future_dates = np.arange(horizon)

        # 1. Forecast Scores (etahat)
        k_etahat_fc = run_multivariate_forecaster(
            self.model_kdfpc.etahat, 
            maxlags_       = self.maxlags, 
            criteria_      = self.criteria, 
            h_             = horizon, 
            selected_nlags = var_lags
        )
        
        # 2. Reconstruct curves from forecasted scores
        k_curve_forecast = self.model_kdfpc.predict(k_etahat_fc)
        # Yhat_fc = k_curve_forecast.copy()
        # A. Enforce positivity: Replace negative "probabilities" with 0
        k_curve_forecast[k_curve_forecast < 0] = 0
        k_curve_forecast.columns = forecast_index
        # B. Renormalize: Iterate through columns by positional index
        for t in range(k_curve_forecast.shape[1]):
            # Use .iloc to access the t-th column by position
            current_col = k_curve_forecast.iloc[:, t]
            integral = current_col.sum() * self.du
            if integral > 0:
                # Update the column in-place using .iloc
                k_curve_forecast.iloc[:, t] = current_col / integral

        # 3. Return a dictionary similar to the blueprint
        # We return the raw reconstruction; alignment happens outside
        df_support = self.u.copy()
        df_support.name = forecast_index[0]
        return {
            "future_scores":    k_etahat_fc,
            "future_support":   df_support,
            "future_curves":    pd.DataFrame(k_curve_forecast)
        }