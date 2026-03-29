import pandas as pd
import numpy as np
from src.fda.dfpc import K_dFPC
from src.forecasting.models import run_forecaster
from src.fda.transformations.lqdt import mLQDT
from typing import List, Tuple, Union


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
        k_scores = KdFPC_model.etahat.values    

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
        k_curve_forecast = KdFPC_model.predict(k_etahat_fc).values
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