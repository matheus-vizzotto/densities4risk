import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
import joblib
from sklearn.model_selection import ParameterGrid
from tqdm.auto import tqdm

import os
import logging
import sys
# sys.path.append('')
import src.forecasting.simulations      as sim
import src.fda.kde.estimators           as kde
import src.fda.transformations.lqdt     as lqdt

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