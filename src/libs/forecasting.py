import matplotlib.pyplot as plt
import pandas as pd
from statsmodels.tsa.stattools import adfuller, kpss#, phillips_perron

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
    model = VAR(data)
    res = model.fit(nlags)
    return res

def forecast_var(res, steps: int = 1) -> np.ndarray:
    """
    Forecast using a fitted statsmodels VARResults object.
    Returns a numpy array (steps x k).
    """
    return res.forecast(res.endog[-res.k_ar:], steps=steps)