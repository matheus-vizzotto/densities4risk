from datetime import datetime
from pathlib import Path
import joblib
import multiprocessing as mp

import os
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import src.forecasting.simulations      as sim
import src.fda.kde.estimators           as kde
import src.fda.transformations.lqdt     as lqdt
from src.forecasting.parallel_utils import run_parallel
import src.fda.utils                    as fdaUtils
import src.forecasting.pipelines        as fp
import src.forecasting.accuracy         as acc
import src.forecasting.cross_validation as crossVal

EXECUTION_DATE : str = datetime.now().strftime('%Y%m%d')
LOG_PATH       : str = f'logs/simulations/{EXECUTION_DATE}.log'
FILES_PATH     : str = f'data/interim/simulation/{EXECUTION_DATE}/'
DATABASES_PATH : str = f'data/interim/simulation/{EXECUTION_DATE}/databases/'
SIM_LQD_PATH       : str = f'data/interim/simulation/{EXECUTION_DATE}/sim_lqds/'
KDE_LQD_PATH   : str = f'data/interim/simulation/{EXECUTION_DATE}/kde_lqds/'
ALL_KDES_PATH  : str = f'data/interim/simulation/{EXECUTION_DATE}/all_kdes.jbl'
CV_FILES_PATH  : str = f'data/interim/simulation/{EXECUTION_DATE}/progress/'
CV_FC_PATH     : str = f'data/interim/simulation/{EXECUTION_DATE}/cv/'


Path(FILES_PATH).mkdir(parents=True, exist_ok=True)
Path(DATABASES_PATH).mkdir(parents=True, exist_ok=True)
Path(CV_FILES_PATH).mkdir(parents=True, exist_ok=True)
Path(SIM_LQD_PATH).mkdir(parents=True, exist_ok=True)
Path(KDE_LQD_PATH).mkdir(parents=True, exist_ok=True)
Path(CV_FC_PATH).mkdir(parents=True, exist_ok=True)


# Simulations
simulations_database = joblib.load(f"{FILES_PATH}simulation_database.jbl")

# KDE
kde_jbls_path = CV_FILES_PATH
kdes_paths = [''.join([kde_jbls_path, x]) for x in os.listdir(kde_jbls_path)]
kde_addresses = []
for path_name in kdes_paths:
    obj = joblib.load(path_name)
    kde_address = {
        "scenario":   obj["scenario"],
        "n_rep":      obj["n_rep"],
        "model_name": obj["model_name"],
        "address":    path_name
    }
    kde_addresses.append(kde_address)
kde_addresses = sorted(
    kde_addresses,
    key=lambda x: (x["scenario"], x["n_rep"], x["model_name"])
)
kde_index = {
    (d["scenario"], d["n_rep"], d["model_name"]): d["address"]
    for d in kde_addresses
}

# LQD(kde)
# kdes_lqd_paths = [''.join([KDE_LQD_PATH, x]) for x in os.listdir(KDE_LQD_PATH) if "scenario" in x]
kdes_lqd_paths = [''.join([KDE_LQD_PATH, x]) for x in os.listdir(KDE_LQD_PATH)]
kdes_lqd_addresses = []
for path_name in kdes_lqd_paths:
    obj = joblib.load(path_name)
    kde_address = {
        "scenario":   obj["scenario"],
        "n_rep":      obj["n_rep"],
        "model_name": obj["model_name"],
        "address":    path_name
    }
    kdes_lqd_addresses.append(kde_address)

kde_lqd_index = {
    (d["scenario"], d["n_rep"], d["model_name"]): d["address"]
    for d in kdes_lqd_addresses
}

KdFPC_kwargs = {
    "p": 5,
    # "select_ncomp": "variance",
    "dimension": 3
}

def cv_worker(kde_address):
    import joblib
    import pandas as pd

    initial_window = 200
    horizon = 1

    scenario_kde_dict = joblib.load(kde_address["address"])

    # Metadata
    id_scenario_kde = scenario_kde_dict["scenario"]
    sim_info = simulations_database[id_scenario_kde]["params"]
    n_rep = scenario_kde_dict["n_rep"]
    model_name_path = scenario_kde_dict["model_name"]
    model_name = model_name_path.replace("_", " ")
    
    print(
        f"[START] scenario={id_scenario_kde} rep={n_rep} model={model_name_path}",
        flush=True,
    )


    df_support = scenario_kde_dict["df_support"]
    df_densities = scenario_kde_dict["df_densities"]

    # Simulated densities
    sim_densities = simulations_database[id_scenario_kde]["replications"][n_rep]["densities"]
    sim_densities_supp = sim_densities.copy()
    sim_densities_supp.loc[:, :] = sim_densities_supp.index.to_numpy()[:, None]

    # LQD
    kde_lqd_address = kde_lqd_index.get((id_scenario_kde, n_rep, model_name_path))
    scenario_kde_lqd = joblib.load(kde_lqd_address)["kde_mlqdt"]
    scenario_kde_lqd_supp = scenario_kde_lqd.copy()
    scenario_kde_lqd_supp.loc[:, :] = scenario_kde_lqd_supp.index.to_numpy()[:, None]

    # CV
    windows = crossVal.expanding_window_cv(
        df_densities.shape[1],
        h=horizon,
        initial_window=initial_window
    )

    measures = []

    for fold, window in enumerate(windows, start=1):

        idx_train, idx_test = window
        test_date = df_densities.columns[idx_test]

        # Targets
        Y_t_support = scenario_kde_lqd_supp.loc[:, test_date]
        Y_t = scenario_kde_lqd.loc[:, test_date]

        f_hat_t_supp = df_support.loc[:, test_date]
        f_hat_t = df_densities.loc[:, test_date]

        f_t_support = sim_densities_supp.loc[:, test_date]
        f_t = sim_densities.loc[:, test_date]

        # Train
        kde_train_support = df_support.iloc[:, idx_train]
        kde_train = df_densities.iloc[:, idx_train]

        # Forecast
        forecaster = fp.DensityForecaster(kdfpc_kwargs=KdFPC_kwargs, maxlags=10)
        forecaster.fit(kde_train, kde_train_support)
        mdfpc_fc = forecaster.predict(horizon=1, forecast_index=test_date)

        Y_hat_t = mdfpc_fc["future_L2_curves"]
        lambda_inv_Y_hat_t_supp = mdfpc_fc["future_supports"]
        lambda_inv_Y_hat_t = mdfpc_fc["future_densities"]

        # ---------------- MEASURES ---------------- #

        # 1
        acc_measures = acc.overall_measures(forecast=Y_hat_t, test=Y_t)
        measures.append({
            "scenario": id_scenario_kde,
            "n_rep": n_rep,
            "model_name": model_name,
            "fold": fold,
            "comparison": "yHatFc_Y",
            "fc_date": test_date[0],
            **acc_measures
        })

        # 2
        df_supp, df_kde, df_fc = fdaUtils.align_densities(
            f_hat_t_supp, f_hat_t,
            lambda_inv_Y_hat_t_supp, lambda_inv_Y_hat_t,
            lambda_inv_Y_hat_t.columns
        )
        acc_measures = acc.overall_measures(forecast=df_fc, test=df_kde)
        measures.append({
            "scenario": id_scenario_kde,
            "n_rep": n_rep,
            "model_name": model_name,
            "fold": fold,
            "comparison": "fHatFc_fHat",
            "fc_date": test_date[0],
            **acc_measures
        })

        # 3
        df_supp, df_kde, df_fc = fdaUtils.align_densities(
            f_t_support, f_t,
            lambda_inv_Y_hat_t_supp, lambda_inv_Y_hat_t,
            lambda_inv_Y_hat_t.columns
        )
        acc_measures = acc.overall_measures(forecast=df_fc, test=df_kde)
        measures.append({
            "scenario": id_scenario_kde,
            "n_rep": n_rep,
            "model_name": model_name,
            "fold": fold,
            "comparison": "fHatTp1Fc_f",
            "fc_date": test_date[0],
            **acc_measures
        })

    # Save INSIDE worker (important to avoid memory blowup)
    output_path = f"{CV_FC_PATH}cvFc___scenario_{id_scenario_kde}___rep_{n_rep}___kde_{model_name_path}.jbl"
    joblib.dump(measures, output_path)

    print(
        f"[DONE] scenario={id_scenario_kde} rep={n_rep} model={model_name_path} "
        f"measures={len(measures)} output={output_path}",
        flush=True,
    )

    return {
        "scenario": id_scenario_kde,
        "n_rep": n_rep,
        "model": model_name,
        "n_measures": len(measures)
    }

def main():
    n_cores = 2
    tasks = kde_addresses

    print(f"Running {len(tasks)} KDE CV tasks with n_cores={n_cores}...", flush=True)

    if n_cores == 1:
        results = []
        for i, task in enumerate(tasks, start=1):
            print(
                f"[TASK {i}/{len(tasks)}] scenario={task['scenario']} "
                f"rep={task['n_rep']} model={task['model_name']}",
                flush=True,
            )
            results.append(cv_worker(task))
        return results

    return run_parallel(
        func=cv_worker,
        data=tasks,
        chunksize=1,      # important: tasks are heavy
        n_cores=n_cores
    )


if __name__ == "__main__":
    results = main()
