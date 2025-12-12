import numpy as np
from plot_utils import *

# --- Load data and extract Fermi energy ---
filename = "tdmft_triqs_1o_DS1_DOS_AT0001"

with open(filename, "r") as f:
    for line in f:
        if "Fermi energy" in line:
            fermi = float(line.split()[-1])
            break

data = np.loadtxt(filename)
data[:, 0] -= fermi  # Shift energy axis
data[:, 0] *= HA_EV
data[:, 1:] /= HA_EV

fig = Figure()

params = {}
params["xlabel"] = r"E - E$_{\mathrm{F}}$ (eV)"
params["ylabel"] = "Number of states / eV"
params["xlim"]   = (-8,5)

fig.add_data(data[:, 0], data[:, 3], params)
fig.plot()
