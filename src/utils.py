import pandas as pd
from src.preprocessing import align_densities
# import pmdarima as pm

# from src.transformations import (
#                             obtain_densities, 
#                             mLQDT,
#                             obtain_densities_from_lqd
#                             )
# from src.forecasting import run_forecaster, train_test_split, overall_measures, expanding_window_cv, slice_fold
# from src.dynamicFPC import K_dFPC

# def cv(Y, Y_support, horizon=1, initial_window=100):
#     windows = expanding_window_cv(Y.shape[1], h=horizon, initial_window=initial_window)

#     measures = []
#     for fold, window in enumerate(windows):
#         print(f"cv {fold}/{len(windows)}")
#         idx_train = window[0]
#         idx_test  = window[1]
        
#         Y_train_support, Y_train = Y_support.iloc[:,idx_train], Y.iloc[:,idx_train]
#         Y_test_support , Y_test = Y_support.iloc[:,idx_test],  Y.iloc[:,idx_test]

#         bovespa_mLQDT = mLQDT(
#                     Y_train,
#                     Y_train_support
#                 )
#         bovespa_mLQDT.densities_to_lqdensities(verbose=False)
#         df_lqds = bovespa_mLQDT.lqd.iloc[1:-1,:]
#         df_lqds_support = bovespa_mLQDT.lqd_support[1:-1]

#         KdFPC_kwargs = {
#             "lag_max": 5,
#             "alpha": 0.10,
#             "du": 0.05,
#             "B": 1000,
#             "p": 5,
#             # "m": df_lqds.shape[0],
#             "u": df_lqds_support,
#             "select_ncomp": False,
#             "dimension": 2
#         }

#         KdFPC_model = K_dFPC(df_lqds.values)
#         KdFPC_model.fit(**KdFPC_kwargs)
#         k_scores = KdFPC_model.etahat.real.T

#         maxlags_  = 10
#         criteria_ = 'bic'

#         k_etahat_fc = run_forecaster(k_scores, maxlags_, criteria_, horizon)

#         model = pm.auto_arima(
#             bovespa_mLQDT.c,                         # univariate series
#             seasonal=False,            # True if SARIMA
#             error_action='ignore',     # ignore non-invertible models
#             suppress_warnings=True,
#             information_criterion='bic',
#             trace=False
#         )

#         # Forecast h steps ahead
#         c_forecast, conf_int = model.predict(n_periods=horizon, return_conf_int=True)

#         k_curve_forecast = KdFPC_model.predict(k_etahat_fc)

#         df_k_forecast = pd.DataFrame(k_curve_forecast, columns=Y_test.columns)

#         kle_bkw_supports, kle_bkw_densities = obtain_densities_from_lqd(
#                                                                     df_k_forecast,
#                                                                     bovespa_mLQDT.lqd_support,
#                                                                     c_forecast,
#                                                                     verbose=False
#                                                                     )
        
#         df_supp, df_f_kle, df_kle_fhat = align_densities(
#                                         Y_test_support, 
#                                         Y_test, 
#                                         kle_bkw_supports, 
#                                         kle_bkw_densities, 
#                                         kle_bkw_densities.columns)

        
#         d1 = {
#             "fold": fold,
#             "method": "KLE",
#         }
#         d1.update(overall_measures(test=df_f_kle, forecast=df_kle_fhat))
#         measures.append(d1)

#     return measures

def mark_repeated_columns(df, threshold=10, dropna=True):
    """
    Identify columns that contain any value repeated more than `threshold` times.
    Function to detect weird data, like returns on Ash Wednesday (started on 13pm).
    Already detected: '2025-03-05', '2025-08-01'.

    Parameters
    ----------
    df : pd.DataFrame
        Columns represent days.
    threshold : int
        Maximum allowed repetitions of a single value.
    dropna : bool
        Whether to ignore NaNs when counting repetitions.

    Returns
    -------
    bad_cols : list
        Column names with excessive repetition.
    mask : pd.Series
        Boolean mask indexed by columns.
    """

    def has_too_many_repeats(col):
        counts = col.value_counts(dropna=dropna)
        return counts.max() > threshold if not counts.empty else False

    mask = df.apply(has_too_many_repeats, axis=0)

    bad_cols = mask[mask].index.tolist()

    return bad_cols, mask

def read_hdfstore(file_path, name):
    with pd.HDFStore(file_path) as store:
        db = {
            "df_h": store[f"{name}/df_h"],
            "df_grids": store[f"{name}/df_grids"],
            "df_densities": store[f"{name}/df_densities"]
        }

        return db