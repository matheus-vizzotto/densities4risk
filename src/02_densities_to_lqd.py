import sys
import os
import json
from pathlib import Path
import numpy as np
import pandas as pd

# Add the project root (parent of "notebooks") to sys.path
# root = Path.cwd().parent
root = Path.cwd()
if str(root) not in sys.path:
    sys.path.append(str(root))

from src.libs.transformations import T_q


# Reads json file of densities with support values
json_file_path = root / "data/processed/densities.json"
with open(json_file_path, 'r') as j:
     densities = json.loads(j.read())

# Obtains functions in L^2[0,1]
lqds = {}
lqds_values = []
col_names   = []
for density_index in densities.keys():
    dens_support = np.array(densities[density_index]["x_grid"])
    dens = np.array(densities[density_index]["y_kde"])
    lqd_support, lqd, c = T_q(dens, dens_support)
    # No need to store support, since the transformation maps all densities to the domain [0,1]
    lqds[density_index] = {
        "lqd_support": lqd_support.tolist(),
        "lqd_c": c,
        "lqd": lqd.tolist()
        }
    # For dataframe
    lqds_values.append(lqd)
    col_names.append(f"t_{density_index}")

# Writes transformed data
file_path = "data/processed/lqdensities.json"
with open(file_path, 'w') as json_file:
    json.dump(lqds, json_file, indent=4)

# Stores Excel file to facilitate the curve decomposition step
df_lqds = pd.DataFrame(lqds_values).T
df_lqds.columns = col_names
df_lqds["x"] = lqd_support # domain is equal for all lqdensities
df_lqds.set_index("x", inplace=True)
df_lqds.to_excel("data/processed/lqdensities.xlsx")