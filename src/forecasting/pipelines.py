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
    A pipeline for forecasting probability density functions using LQD 
    transformation and Functional Principal Component Analysis (K-dFPC).

    Attributes
    ----------
    model_kdfpc : object
        The fitted K_dFPC model instance.
    eigenvalues : numpy.ndarray
        Eigenvalues (λ) representing variance explained by each component.
    eigenfunctions : numpy.ndarray
        The basis functions (φ) in the LQD space.
    scores : numpy.ndarray
        The expansion coefficients (η) for the training data.

    Example
    -------
    >>> forecaster = DensityForecaster(maxlags=5)
    >>> # Fit the model
    >>> forecaster.fit(Y_train, Y_support_train)
    >>> # Access stored K_dFPC attributes
    >>> forecaster.eigenvalues
    >>> # Generate Forecast
    >>> supports, densities = forecaster.predict(horizon=10)
    """
    def __init__(self, kdfpc_kwargs=None, maxlags=10, criteria='bic'):
        self.kdfpc_kwargs = kdfpc_kwargs or {}
        self.maxlags = maxlags
        self.criteria = criteria
        self.model_lqd = None
        self.model_kdfpc = None
        self.model_arima = None

    def fit(self, Y_train, Y_support):
        # 1. Transform Densities
        mlqdt = lqdt.mLQDT()
        self.model_lqd = mlqdt.transform(
            densities=Y_train,
            densities_supports=Y_support)
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

    def predict(self, horizon, var_lags=None):
        # 1. Forecast Scores
        k_scores = self.model_kdfpc.etahat.values.real
        k_etahat_fc = run_multivariate_forecaster(
            k_scores, self.maxlags, self.criteria, horizon, selected_nlags=var_lags
        )
        
        # 2. Forecast 'c'
        c_forecast = self.model_arima.predict(n_periods=horizon)
        
        # 3. Reconstruct LQD and back to Densities
        L2_curve_forecast = self.model_kdfpc.predict(k_etahat_fc)
        df_L2_curve_forecast = pd.DataFrame(L2_curve_forecast)

        fc_lqd_obj = lqdt.LQDRepresentation(
                            lqd         = df_L2_curve_forecast,
                            c           = c_forecast,
                            lqd_support = self.model_lqd.lqd_support,
                            t0          = self.model_lqd.t0 
                            )
        
        mlqdt = lqdt.mLQDT()
        df_supports, df_densities = mlqdt.inverse_transform(fc_lqd_obj)
        
        return {
                "future_L2_curves": df_L2_curve_forecast,
                "future_scores": k_etahat_fc,
                "future_cs": c_forecast,
                "future_supports": df_supports,
                "future_densities": df_densities
        }