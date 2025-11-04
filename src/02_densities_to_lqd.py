import sys
import os
import json
from pathlib import Path
import numpy as np

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
for density_index in densities.keys():
    dens_support = np.array(densities[density_index]["x_grid"])
    dens = np.array(densities[density_index]["y_kde"])
    lqd_support, lqd, c = T_q(dens, dens_support)
    lqds[density_index] = {
        "lqd_support": lqd_support.tolist(),
        "lqd": lqd.tolist()
        }
    
# Writes transformed data
file_path = "data/processed/lqdensities.json"
with open(file_path, 'w') as json_file:
    json.dump(lqds, json_file, indent=4)