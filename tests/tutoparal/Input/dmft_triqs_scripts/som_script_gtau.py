import numpy as np
import triqs.utility.mpi as mpi
from triqs.gf import *
from som import Som, fill_refreq, reconstruct
from plot_utils import *

# Retrieve data from G(tau) file
filename = "tdmft_triqs_2o_DS2_Gtau_diag_DLR_iatom0001.dat"
data = np.loadtxt(filename)

beta  = data[-1, 0]
n_tau = data.shape[0]
n_orb = data.shape[1] - 1

# Build input Green's function for SOM
g_tau = GfImTime(target_shape = [n_orb,n_orb], beta = beta, n_points = n_tau)

for i in range(n_orb):
    g_tau.data[:, i, i] = data[:, i+1]

# Assume constant error bars for every tau (only relative weights matter)
error_bars = g_tau.copy()
error_bars.data[:, :, :] = 0.01

# SOM only analytically continues the diagonal components, whose spectral functions have norm 1
norms = np.ones(n_orb)

# Parameters for the analytical solution: do not use these as reference, you need to carefully select them
acc_params = {}

energy_window = (-4.0, 4.0)

# Support of the spectral function
acc_params['energy_window'] = energy_window
# Number of particular solutions to accumulate
acc_params['l'] = 10
# Number of global updates
acc_params['f'] = 100
# Number of local updates per global update
acc_params['t'] = 50

cont = Som(g_tau, error_bars, kind="FermionGf", norms=norms)

# Accumulate solutions
cont.accumulate(**acc_params)
# Build final solution
cont.compute_final_solution(good_chi_abs=10.0, good_chi_rel=2.0)

# Build local Green's function on the real axis
n_w = 1000
g_w = GfReFreq(window=energy_window, n_points=n_w, indices=g_tau.indices)
fill_refreq(g_w, cont)

# Extract spectral function
w_mesh = np.fromiter(g_w.mesh.values(), float)
data_spectral = np.zeros((n_w,n_orb+1))
data_spectral[:, 0] = w_mesh

for i in range(n_orb):
    data_spectral[:, i+1] = - g_w.data[:, i, i].imag / np.pi

if mpi.is_master_node():
    # Save on file
    np.savetxt("spectral.dat", data_spectral)

    # Reconstruct G(tau) from spectral function and compare with input
    g_rec = g_tau.copy()
    reconstruct(g_rec, cont)

    fig = Figure()

    params = {}
    params["label"]  = r"Input $G(\tau)$"
    params["xlabel"] = r"$\tau$ (Ha$^{-1}$)"
    params["ylabel"] = r"$G(\tau)$"

    fig.add_data(data[:, 0], g_tau.data[:, 0, 0], params)

    params["label"] = r"Reconstructed $G(\tau)$"

    fig.add_data(data[:, 0], g_rec.data[:, 0, 0], params)
    fig.plot(legend=True)

