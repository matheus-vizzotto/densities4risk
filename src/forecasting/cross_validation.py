import pandas as pd
import numpy as np
from src.fda.dfpc import K_dFPC
from src.forecasting.pipelines import run_multivariate_forecaster
from src.fda.transformations.lqdt import mLQDT
import src.fda.utils             as fdaUtils
import src.fda.dfpc              as dfpc
import src.forecasting.pipelines as fp
import src.forecasting.accuracy  as acc 

from typing import List, Tuple, Union, Literal

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

def cv(
        Y, 
        Y_support, 
        KdFPC_kwargs, 
        horizon=1, 
        initial_window=100,
        var_lags: int = None,
        return_curves=False
        ):
    windows = expanding_window_cv(Y.shape[1], h=horizon, initial_window=initial_window)

    curves_hist = []
    measures = []
    for fold, window in enumerate(windows):
        print(f"\t\t>>> cv {fold+1}/{len(windows)}")
        idx_train = window[0]
        idx_test  = window[1]
        
        # TRAIN-TEST SPLIT FOR DENSITIES AND SUPPORTS
        Y_train_support, Y_train = Y_support.iloc[:,idx_train], Y.iloc[:,idx_train]
        Y_test_support , Y_test = Y_support.iloc[:,idx_test],  Y.iloc[:,idx_test]

        # initialize class
        forecaster = fp.DensityForecaster(kdfpc_kwargs=KdFPC_kwargs, maxlags=10)
        # fit model
        forecaster.fit(Y_train, Y_train_support)
        # forecast dict
        mdfpc_fc = forecaster.predict(horizon=1, var_lags=var_lags, forecast_index=Y_test.columns)
        
        # PUT KDE AND FORECAST INTO THE SAME GRID FOR EVALUATION
        df_supp, df_kde, df_fc = fdaUtils.align_densities(
                                        Y_test_support, 
                                        Y_test, 
                                        mdfpc_fc["future_supports"], 
                                        mdfpc_fc["future_densities"], 
                                        mdfpc_fc["future_densities"].columns
                                        )

                        
        # COMPUTE ACCURACY MEASURES AND STORE THEM
        oa_measures = acc.overall_measures(test=df_kde, forecast=df_fc)
        d1 = {
            "fold": fold,
            "method": "KLE",
            }
        d1.update(oa_measures)
        measures.append(d1)

        if return_curves:
            # d1.update({"df_supports": df_supp, "df_kde": df_f_kle, "df_forecast": df_kle_fhat})
            temp_df = pd.DataFrame({
                    "support":  df_supp.iloc[:, 0].values,
                    "actual":   df_kde.iloc[:, 0].values,
                    "forecast": df_fc.iloc[:, 0].values,
                    "date":     df_fc.columns[0],
                    "fold":     fold
                    })
            temp_df.set_index(["fold", "date", "support"], inplace=True)
            curves_hist.append(temp_df) 
        
    return measures, pd.concat(curves_hist)





def cv_dfpc(
        Y, 
        Y_support, 
        KdFPC_kwargs, 
        horizon=1, 
        initial_window=100,
        return_curves=False,
        var_lags: int = None
        ):
    windows = expanding_window_cv(Y.shape[1], h=horizon, initial_window=initial_window)

    curves_hist = []
    measures = []
    for fold, window in enumerate(windows):
        print(f"\t\t>>> cv {fold+1}/{len(windows)}")
        idx_train = window[0]
        idx_test  = window[1]
        
        # TRAIN-TEST SPLIT FOR DENSITIES AND SUPPORTS
        Y_train_support, Y_train = Y_support.iloc[:,idx_train], Y.iloc[:,idx_train]
        Y_test_support , Y_test = Y_support.iloc[:,idx_test],   Y.iloc[:,idx_test]

        # initialize class
        KdFPC_kwargs.update({
            'lag_max': 5, 
            'B': 1000,
            'alpha': 0.01
        })
        forecaster = fp.DensityForecaster_HZ(kdfpc_kwargs=KdFPC_kwargs, maxlags=10)
        # fit model
        forecaster.fit(Y_train, Y_train_support)
        # forecast dict
        dfpc_fc = forecaster.predict(horizon=1, var_lags=var_lags, forecast_index=Y_test.columns)
        dfpc_fc_dens = dfpc_fc["future_curves"]
        dfpc_fc_supp = dfpc_fc["future_support"].to_frame()
        
        # PUT KDE AND FORECAST INTO THE SAME GRID FOR EVALUATION
        df_supp, df_kde, df_fc = fdaUtils.align_densities(
                                        Y_test_support, 
                                        Y_test, 
                                        dfpc_fc_supp, 
                                        dfpc_fc_dens, 
                                        Y_test_support.columns
                                        )

        # COMPUTE ACCURACY MEASURES AND STORE THEM
        oa_measures = acc.overall_measures(test=df_kde, forecast=df_fc)
        d1 = {
            "fold": fold,
            "method": "KLE",
            }
        d1.update(oa_measures)
        d1.update({"df_supports": df_supp, "df_kde": df_kde, "df_forecast": df_fc})
        measures.append(d1)


        if return_curves:
            temp_df = pd.DataFrame({
                    "support":  df_supp.iloc[:, 0].values,
                    "actual":   df_kde.iloc[:, 0].values,
                    "forecast": df_fc.iloc[:, 0].values,
                    "date":     df_fc.columns[0],
                    "fold":     fold
                    })
            temp_df.set_index(["fold", "date", "support"], inplace=True)
            curves_hist.append(temp_df) 
    
    return measures, pd.concat(curves_hist)

def cv_sim(
        Y, 
        Y_support, 
        simulated_Y,
        simulated_Y_support,
        KdFPC_kwargs, 
        horizon=1, 
        initial_window=100,
        var_lags: int = None,
        return_curves=False
        ):

        windows = expanding_window_cv(Y.shape[1], h=horizon, initial_window=initial_window)

        curves_hist = []
        measures = []
        for fold, window in enumerate(windows):
                print(f"\t\t>>> cv {fold+1}/{len(windows)}")
                idx_train      = window[0]
                idx_test       = window[1]
                test_date      = Y.columns[idx_test]

                # TRAIN-TEST SPLIT FOR DENSITIES AND SUPPORTS
                Y_train_support, Y_train = Y_support.iloc[:,idx_train], Y.iloc[:,idx_train]
                densities_test = simulated_Y.loc[:,test_date]
                support_test   = simulated_Y_support.loc[:,test_date] 
                # Y_test  = df_densities.iloc[:,idx_test]
                # Y_test_support = densities_test.copy()
                # Y_test_support.loc[:,:] = Y_test_support.index.to_numpy()[:, None]

                # initialize class
                forecaster = fp.DensityForecaster(kdfpc_kwargs=KdFPC_kwargs, maxlags=10)
                # fit model
                forecaster.fit(Y_train, Y_train_support)
                # forecast dict
                mdfpc_fc = forecaster.predict(horizon=1, var_lags=var_lags, forecast_index=densities_test.columns)

                # PUT KDE AND FORECAST INTO THE SAME GRID FOR EVALUATION
                df_supp, df_kde, df_fc = fdaUtils.align_densities(
                                                support_test, 
                                                densities_test, 
                                                mdfpc_fc["future_supports"], 
                                                mdfpc_fc["future_densities"], 
                                                mdfpc_fc["future_densities"].columns
                                                )

                                
                # COMPUTE ACCURACY MEASURES AND STORE THEM
                oa_measures = acc.overall_measures(test=df_kde, forecast=df_fc)
                d1 = {
                        "fold": fold,
                        "method": "KLE",
                        }
                d1.update(oa_measures)
                measures.append(d1)

                if return_curves:
                        temp_df = pd.DataFrame({
                                "support":  df_supp.iloc[:, 0].values,
                                "actual":   df_kde.iloc[:, 0].values,
                                "forecast": df_fc.iloc[:, 0].values,
                                "date":     df_fc.columns[0],
                                "fold":     fold
                                })
                        temp_df.set_index(["fold", "date", "support"], inplace=True)
                        curves_hist.append(temp_df) 
                        d1.update({"curves": temp_df})
                        d1.update({"curves_fc_date": df_fc.columns[0]})
                
        return measures, pd.concat(curves_hist)


from typing import List, Tuple, Literal, Optional
import numpy as np


from typing import List, Tuple, Literal
import numpy as np


def cv_window(
    T: int,
    h: int = 1,
    window_type: Literal["expanding", "rolling"] = "expanding",
    window_size: int = 50,
    step: int = 1,
) -> List[Tuple[np.ndarray, np.ndarray]]:
    """
    Generate time series cross-validation splits using
    expanding or rolling windows.

    Parameters
    ----------
    T : int
        Total number of observations.

    h : int, default=1
        Forecast horizon.

    window_type : {"expanding", "rolling"}, default="expanding"
        Cross-validation scheme.

        - "expanding":
            Training set grows over time.

        - "rolling":
            Fixed-size moving training window.

    window_size : int, default=50
        Training window size.

        - Expanding:
            Initial training size.

        - Rolling:
            Fixed rolling training size.

    step : int, default=1
        Step size between folds.

    Returns
    -------
    List[Tuple[np.ndarray, np.ndarray]]
        List of `(train_idx, test_idx)` tuples.

    Examples
    --------
    Expanding window:

        train = [0, 1, 2]
        test  = [3]

        train = [0, 1, 2, 3]
        test  = [4]

    Rolling window:

        train = [0, 1, 2]
        test  = [3]

        train = [1, 2, 3]
        test  = [4]
    """

    # =========================
    # Validation
    # =========================

    if T <= 0:
        raise ValueError("T must be positive.")

    if h <= 0:
        raise ValueError("h must be positive.")

    if window_size < 0:
        raise ValueError("window_size must be >= 0.")

    if step <= 0:
        raise ValueError("step must be positive.")

    if window_type not in {"expanding", "rolling"}:
        raise ValueError(
            "window_type must be either "
            "'expanding' or 'rolling'."
        )

    splits = []

    # =====================================
    # Expanding window
    # =====================================

    if window_type == "expanding":

        for t in range(window_size, T, step):

            test_end = t + h

            if test_end > T:
                break

            train_idx = np.arange(0, t)
            test_idx = np.arange(t, test_end)

            splits.append((train_idx, test_idx))

    # =====================================
    # Rolling window
    # =====================================

    else:

        for test_start in range(window_size, T, step):

            test_end = test_start + h

            if test_end > T:
                break

            train_start = test_start - window_size
            train_end = test_start

            train_idx = np.arange(train_start, train_end)
            test_idx = np.arange(test_start, test_end)

            splits.append((train_idx, test_idx))

    return splits