import numpy as np
from plot_utils import *

# Retrieve list of k-pts
filename = "tdmft_triqs_3o_DS3_EIG"

k_mesh = []

with open(filename, "r") as f:
    for line in f:
        if "kpt#" in line:
            line_split = line.split()
            k_mesh.append(float(line_split[-4]))

# Retrieve k-resolved spectral function
w_mesh = []
filename = "tdmft_triqs_3o_DS3_DFTDMFT_SpectralFunction_kres"

ikpt = 1
data = []
data_spectral = []

with open(filename, "r") as f:
    for line in f:
        line_split = line.split()
        if ("#" in line) or (len(line_split) == 0):
            continue
        ikpt_ = int(line_split[-1])
        if ikpt_ != ikpt: # New kpt
            data_spectral.append(data)
            w_mesh = []
            data = []
        ikpt = ikpt_
        w_mesh.append(float(line_split[0]))
        data.append(float(line_split[1]))

data_spectral.append(data)
data_spectral = np.array(data_spectral)
data_spectral = np.transpose(data_spectral)

fig = Figure()

params = {}
params["ylim"] = (-5,5)
params["xticks"] = [0,0.1,0.2,0.3,0.4,0.5]
params["xticklabels"] = [r"$\Gamma$","","","","",r"N"]
params["ylabel"] = r"Energy (eV)"

fig.add_data(k_mesh, w_mesh, params, z=data_spectral)
fig.plot(colorbar=True)

