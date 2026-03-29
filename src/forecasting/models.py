import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from statsmodels.tsa.api import VAR
import pmdarima as pm
from statsmodels.tsa.vector_ar.var_model import VARResults
from statsmodels.graphics.tsaplots import plot_acf

from typing import Optional, Dict, Union


class multivariate_forecaster:
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

def univariate_forecaster(
        series:  pd.Series, 
        horizon: int  = 1, 
        ic:      str  = 'bic', 
        verbose: bool = False
        ):
    """
    Ajusta um modelo Auto-ARIMA a uma série de coeficientes e gera previsões.
    
    Esta função automatiza a seleção da ordem (p, d, q) com base em critérios 
    de informação, facilitando o forecasting de objetos funcionais.

    Args:
        series (pd.Series ou np.ndarray): A série temporal do coeficiente (ex: c_it).
        horizon (int): Número de períodos à frente para prever.
        ic (str): Critério de informação ('bic' ou 'aic'). O BIC é geralmente 
                  preferido em finanças por ser mais parcimonioso.
        seasonal (bool): Se True, busca componentes sazonais (SARIMA).

    Returns:
        tuple: (previsão_pontual, intervalo_confiança)
    """
    model = pm.auto_arima(
        series,
        information_criterion   = ic,
        seasonal                = False,
        error_action            = 'ignore',
        suppress_warnings       = True,
        stepwise                = True,
        trace                   = False
    )

    # Argumento para imprimir o sumário do modelo (p, d, q) selecionado
    if verbose:
        print(f"\n--- Model Summary for {series.name if hasattr(series, 'name') else 'Coefficient'} ---")
        print(model.summary())
    
    # Gera a previsão e o intervalo de confiança
    forecast, conf_int = model.predict(
        n_periods=horizon, 
        return_conf_int=True
    )
    
    return forecast