import numpy as np
import pandas as pd
import src.fda.dfpc as dfpc
import src.fda.transformations.lqdt as lqdt
import src.fda.kde.estimators           as kde
import src.forecasting.models as fm
import src.forecasting.simulations as fsim
import pmdarima as pm
from scipy.optimize import minimize

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
            Y_support:  pd.DataFrame,
            fpc_style:  str = "dynamic"
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
        
        if   fpc_style == "dynamic":
            self.model_kdfpc = dfpc.K_dFPC(lqd_values)
        elif fpc_style == "static":
            self.model_kdfpc = dfpc.K_sFPC(lqd_values)
        else:
            raise ValueError("FPC style not available. Choose one of the following: ['dynamic', 'static'].")
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
    

class SimulationForecaster:
    """
    Time-series forecaster for Probability Density Functions (PDFs).

    This class implements a forecasting pipeline for functional data (densities) 
    by mapping them into a Hilbert space using the Log-Quantile Density (LQD) 
    transformation. It reduces dimensionality via Karhunen-Loève dynamic Functional 
    Principal Component (K-dFPC) Analysis and forecasts the resulting 
    coefficients using a multivariate ScoreForecaster.

    Parameters
    ----------
    kdfpc_kwargs : dict, optional
        Keyword arguments passed to the `dfpc.K_dFPC` model (e.g., n_components).
    forecaster_model : str, default='ridge'
        The underlying model type used by the ScoreForecaster.
    lag : int, default=3
        The number of lags to consider for the score forecaster.

    Attributes
    ----------
    model_lqd : lqdt.mLQDT
        The fitted LQD transformation object containing LQD curves and constants.
    model_kdfpc : dfpc.K_dFPC
        The fitted Functional PCA model.
    score_forecaster : ScoreForecaster
        The fitted multivariate forecaster for PCA scores.

    Example
    ----------
    forecaster = SimulationForecaster(kdfpc_kwargs=KdFPC_kwargs, forecaster_model="elasticnet")
    forecaster.fit(Y_train, Y_supp_train, returns=returns_train)
    results = forecaster.predict(horizon=1, forecast_index=Y_test.to_frame().columns)
    """
    
    def __init__(self, kdfpc_kwargs=None, forecaster_model="ridge", lag=3):
        self.kdfpc_kwargs = kdfpc_kwargs or {}
        self.forecaster_model = forecaster_model
        self.lag = lag
        
        # Initialized during fitting
        self.returns = None
        self.model_lqd = None
        self.model_kdfpc = None
        self.score_forecaster = None

    def fit(self, Y_train: pd.DataFrame, Y_support: pd.DataFrame, returns: pd.DataFrame, fpc_style:  str = "dynamic"):
        """
        Fits the LQD transformation, Functional PCA, and the internal ScoreForecaster.
        """
        self.returns = returns
        # 1. Transform Densities to LQD Space
        mlqdt = lqdt.mLQDT()
        self.model_lqd = mlqdt.transform(
            densities=Y_train,
            densities_supports=Y_support, 
            verbose=False
        )
        
        # 2. Extract L2 Expansion (K_dFPC)
        lqd_values = self.model_lqd.lqd.values
        lqd_support = self.model_lqd.lqd_support
        
        # Use a copy to avoid mutating constructor configurations dynamically
        kdfpc_params = self.kdfpc_kwargs.copy()
        kdfpc_params.update({
            "u": lqd_support,
            "du": self.model_lqd.du
        })
        
        if   fpc_style == "dynamic":
            self.model_kdfpc = dfpc.K_dFPC(lqd_values)
        elif fpc_style == "wavelet":
            self.model_kdfpc = dfpc.W_dFPC2(lqd_values)
        elif fpc_style == "static":
            self.model_kdfpc = dfpc.K_sFPC(lqd_values)
        else:
            raise ValueError("FPC style not available. Choose one of the following: ['dynamic', 'static'].")
        self.model_kdfpc.fit(**kdfpc_params)
        
        # 3. Fit Multivariate Score Forecaster
        scores = self.model_kdfpc.etahat.values.real
        
        self.score_forecaster = fm.ScoreForecaster(
            model=self.forecaster_model,
            lag=self.lag
        )
        self.score_forecaster.fit(returns.values, scores)

        return self

    def predict(self, horizon: int, forecast_index: pd.Index = None) -> dict:
        """
        Forecasts future scores, reconstructs LQD curves, and inverts them back to densities.
        """
        if self.score_forecaster is None:
            raise RuntimeError("The forecaster must be fitted before calling predict.")

        # Build forecasted index labels
        if forecast_index is not None:
            if len(forecast_index) != horizon:
                raise ValueError(f"forecast_index length ({len(forecast_index)}) must match horizon ({horizon})")
            future_dates = forecast_index
        else:
            future_dates = np.arange(horizon)
        
        # 1. Forecast Scores
        k_etahat_fc = self.score_forecaster.predict_next(self.returns.values)
        
        # 2. Reconstruct LQD Curves
        L2_curve_forecast = pd.DataFrame(self.model_kdfpc.predict(k_etahat_fc))
        L2_curve_forecast.columns = future_dates
        
        
        return {
            "future_L2_curves": L2_curve_forecast,
            "future_scores": k_etahat_fc
        }

class ASTParameterEstimator:
    """
    Daily MLE estimator for the standardized Fernández-Steel skew-t model.

    State vector:
        theta = (mu, log_sigma, log_a)

    where
        sigma = exp(log_sigma)
        a     = exp(log_a)
    """

    def __init__(self, nu=8):
        self.nu = nu

    def _negative_loglikelihood(self, theta, z):

        mu = theta[0]
        sigma = np.exp(theta[1])
        a = np.exp(theta[2])

        dist = fsim.StandardizedSkewStudentT(
            a=a,
            nu=self.nu
        )

        x = (z - mu) / sigma

        loglik = np.sum(
            -np.log(sigma)
            + dist.density(x, log=True)
        )

        return -loglik

    def fit_day(self, z, x0=None):

        z = np.asarray(z)

        if x0 is None:
            x0 = np.array([
                np.mean(z),
                np.log(np.std(z)),
                0.0
            ])

        res = minimize(
            self._negative_loglikelihood,
            x0=x0,
            args=(z,),
            method="L-BFGS-B"
        )

        return res.x

    def fit(self, samples: pd.DataFrame):

        params = []

        x0 = None

        for col in samples.columns:

            z = samples[col].dropna().values

            theta = self.fit_day(
                z=z,
                x0=x0
            )

            params.append(theta)

            x0 = theta.copy()

        params = pd.DataFrame(
            params,
            columns=[
                "mu",
                "log_sigma",
                "log_a"
            ],
            index=samples.columns
        )

        return params


class ParametricDensityForecaster:
    """
    Parametric benchmark based on daily MLE estimation of an
    asymmetric Student-t distribution followed by VAR forecasting.

    Parameters
    ----------
    forecast : {"density", "kde"}
        density : returns the analytical AST density.
        kde     : draws Monte Carlo samples from the AST forecast and
                  estimates a KDE on the supplied support grid.
    """

    def __init__(
        self,
        nu=3,
        maxlags=10,
        criteria="bic",
        kde_bw_params=None,
        kde_params=None,
        n_kde_forecasts=1,
    ):

        self.nu = nu
        self.maxlags = maxlags
        self.criteria = criteria

        self.kde_bw_params = (
            kde_bw_params
            if kde_bw_params is not None
            else {
                "method": "rot",
                "kernel": "gaussian",
                "sigma_robust": False,
            }
        )

        self.kde_params = (
            kde_params
            if kde_params is not None
            else {
                "kernel": "gaussian",
            }
        )

        self.n_kde_forecasts = self._validate_n_kde_forecasts(
            n_kde_forecasts
        )
        self.model_estimator = None
        self.parameter_history = None

    @staticmethod
    def _validate_n_kde_forecasts(n_kde_forecasts):
        if not isinstance(n_kde_forecasts, (int, np.integer)):
            raise ValueError("n_kde_forecasts must be a positive integer.")

        if n_kde_forecasts < 1:
            raise ValueError("n_kde_forecasts must be a positive integer.")

        return int(n_kde_forecasts)

    def _ast_kde_density(
        self,
        dist,
        mu,
        sigma,
        sample_size,
        support,
        rng,
    ):
        sample = (
            mu
            + sigma
            * dist.random_numbers(
                sample_size,
                random_state=rng,
            )
        )

        returns_df = pd.DataFrame(
            sample,
            columns=["returns"],
        )

        df_h = kde.df_bandwidth_selector(
            returns_df,
            **self.kde_bw_params,
        )

        kde_params = self.kde_params.copy()
        kernel = kde_params.pop("kernel", "gaussian")
        kde_model = kde.KDE(kernel=kernel, **kde_params)

        _, density = kde_model.transform(
            X=returns_df["returns"],
            h=df_h["returns"].values,
            grid=support,
        )

        return density

    def fit(
        self,
        samples: pd.DataFrame
    ):

        self.model_estimator = ASTParameterEstimator(
            nu=self.nu
        )

        self.parameter_history = self.model_estimator.fit(samples)

        return self

    def predict(
        self,
        horizon,
        support,
        forecast="density",
        sample_size=288,
        random_state=None,
        var_lags=None,
        forecast_index=None,
        n_kde_forecasts=None,
    ):

        if forecast not in {"density", "kde"}:
            raise ValueError(
                "forecast must be either 'density' or 'kde'."
            )

        if n_kde_forecasts is None:
            n_kde_forecasts = self.n_kde_forecasts

        n_kde_forecasts = self._validate_n_kde_forecasts(
            n_kde_forecasts
        )

        rng = np.random.default_rng(random_state)
        support = np.asarray(support, dtype=float)

        # --------------------------------------------------
        # Build forecast dates
        # --------------------------------------------------

        if forecast_index is not None:

            if len(forecast_index) != horizon:
                raise ValueError(
                    f"forecast_index length ({len(forecast_index)}) "
                    f"must match horizon ({horizon})"
                )

            future_dates = forecast_index

        elif (
            hasattr(self, "index_freq")
            and self.index_freq is not None
            and self.last_date is not None
        ):

            future_dates = pd.date_range(
                start=self.last_date,
                periods=horizon + 1,
                freq=self.index_freq,
            )[1:]

        else:

            future_dates = np.arange(horizon)

        # --------------------------------------------------
        # Forecast AST parameters
        # --------------------------------------------------

        theta_fc = run_multivariate_forecaster(
            scores=self.parameter_history.values,
            maxlags_=self.maxlags,
            criteria_=self.criteria,
            h_=horizon,
            selected_nlags=var_lags,
        )

        mu_fc = theta_fc[0]
        sigma_fc = np.exp(theta_fc[1])
        a_fc = np.exp(theta_fc[2])

        density_list = []
        support_list = []
        kde_density_list = []
        kde_support_list = []
        kde_columns = []

        # --------------------------------------------------
        # Build forecast densities
        # --------------------------------------------------

        for h in range(horizon):

            dist = fsim.StandardizedSkewStudentT(
                a=a_fc[h],
                nu=self.nu,
            )

            if forecast == "density":

                x = (support - mu_fc[h]) / sigma_fc[h]

                density = dist.density(x) / sigma_fc[h]

            else:

                density_draws = []

                for a in range(n_kde_forecasts):

                    density_a = self._ast_kde_density(
                        dist=dist,
                        mu=mu_fc[h],
                        sigma=sigma_fc[h],
                        sample_size=sample_size,
                        support=support,
                        rng=rng,
                    )

                    density_draws.append(density_a)

                    if n_kde_forecasts > 1:
                        kde_density_list.append(density_a)
                        kde_support_list.append(support)
                        kde_columns.append((future_dates[h], a + 1))

                density = np.mean(density_draws, axis=0)

            density_list.append(density)
            support_list.append(support)

        # --------------------------------------------------
        # Density DataFrames
        # --------------------------------------------------

        future_densities = pd.DataFrame(
            np.asarray(density_list).T,
            index=support,
            columns=future_dates,
        )

        future_supports = pd.DataFrame(
            np.asarray(support_list).T,
            columns=future_dates,
        )

        future_parameters = pd.DataFrame(
            {
                "mu": mu_fc,
                "sigma": sigma_fc,
                "a": a_fc,
            },
            index=future_dates,
        )

        # --------------------------------------------------
        # LQD transformation
        # --------------------------------------------------

        mlqdt = lqdt.mLQDT()

        if forecast == "kde" and n_kde_forecasts > 1:

            kde_columns = pd.MultiIndex.from_tuples(
                kde_columns,
                names=["forecast_date", "kde_forecast"],
            )

            kde_densities = pd.DataFrame(
                np.asarray(kde_density_list).T,
                index=support,
                columns=kde_columns,
            )

            kde_supports = pd.DataFrame(
                np.asarray(kde_support_list).T,
                columns=kde_columns,
            )

            model_lqd_all = mlqdt.transform(
                densities=kde_densities,
                densities_supports=kde_supports,
                verbose=False,
            )

            lqd_forecast = (
                model_lqd_all.lqd
                .T
                .groupby(level="forecast_date")
                .mean()
                .T
            )
            lqd_forecast = lqd_forecast.loc[:, future_dates]

            c_forecast = (
                pd.Series(model_lqd_all.c, index=kde_columns)
                .groupby(level="forecast_date")
                .mean()
                .loc[future_dates]
                .values
            )

            self.model_lqd = lqdt.LQDRepresentation(
                lqd=lqd_forecast,
                lqd_support=model_lqd_all.lqd_support,
                c=c_forecast,
                t0=model_lqd_all.t0,
            )

        else:

            self.model_lqd = mlqdt.transform(
                densities=future_densities,
                densities_supports=future_supports,
                verbose=False,
            )

            lqd_forecast = pd.DataFrame(
                self.model_lqd.lqd.values,
                index=self.model_lqd.lqd_support,
                columns=future_dates,
            )

        return {
            "future_parameters": future_parameters,
            "future_densities": future_densities,
            "future_L2_curves": lqd_forecast,
            "supp": future_supports,
            "n_kde_forecasts": n_kde_forecasts,
        }
