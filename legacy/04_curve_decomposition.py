# performs function decomposition with K_dFPC and W_dFPC to obtain multivariate time series model.

################### SETUP DO AMBIENTE ###################
import sys
from pathlib import Path
root = Path.cwd()
if str(root) not in sys.path:
    sys.path.append(str(root))
#--------------------------------------------------------

import pandas as pd
import numpy as np
import json
from src.libs.dynamicFPC import K_dFPC

########################## Obtain c_1, ..., c_T from mLQDT ##########################
json_file_path = root / "data/processed/lqdensities.json"
with open(json_file_path, 'r') as j:
     densities = json.loads(j.read())
c_t = []
for t in densities.keys():
    c_t.append((t,densities[t]['lqd_c']))

########################## Data ##########################
lqdensities_path = root / "data/processed/lqdensities.xlsx"
lqdensities = pd.read_excel(lqdensities_path, index_col="x")
Y = lqdensities.reset_index(drop=True)
Y = Y.values
u = lqdensities.index.values

########################## KLE Dynamic Functional Principal Components ##########################
m=200
lag_maximum = 6
alpha_val = 0.10
no_boot = 1000
du=0.05
p=5
m=200
D_val = 10

model = K_dFPC(lqdensities.values)
model.fit(
    lag_max=lag_maximum,
    B=no_boot,
    alpha=0.10,
    du=0.05,
    p=5,
    m=lqdensities.shape[0],
    u=u,
    select_ncomp=False,
    dimension=10
)
K_dFPC_fitted_scores = model.etahat

# Saves single file with both the etahat matrix (dFPC) and c_i (mLQDT) needed for forecasting. ({object}_{type}_{dimension})
np.savez(root / "data/processed/K_dFPC_scores.npz", etahat_matrix_d0_T=model.etahat, c_array_T=c_t)
# To read the file:
# data = np.load("../data/processed/K_dFPC_scores.npz")
# matrix = data["etahat_d0_T"]
# array  = data["c_T"]