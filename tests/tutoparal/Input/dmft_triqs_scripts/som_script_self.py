import numpy as np
import os
import triqs.utility.mpi as mpi
from triqs.gf import *
from som import Som, fill_refreq, reconstruct
from plot_utils import *

# Count number of Maxent files to get the number of orbitals
listdir = os.listdir()
filtered_maxent = list(filter(lambda s: "Selfmxent0001" in s, listdir))
n_orb = len(filtered_maxent)

moments = np.zeros(4, dtype="complex")
norms   = np.zeros(n_orb)

# Retrieve frequency mesh
root_name = "tdmft_triqs_3o_DS2_Selfmxent0001_is1_iflav"
iw_mesh   = np.loadtxt(root_name+"0001")[:, 0]

beta = np.pi / iw_mesh[0]
n_iw = len(iw_mesh)

self_iw = GfImFreq(target_shape = [n_orb,n_orb], beta = beta, n_points = n_iw)

# Retrieve self-energy data from files
for i in range(n_orb):

    filename = root_name + str(i+1).rjust(4,"0")

    data = []

    with open(filename, "r") as f:

        for line in f:

            line_split = line.split()

            if "#moments" in line:
                i_mom = int(line_split[0].split("_")[-1])
                moments[i_mom-1] = float(line_split[1]) + 1j*float(line_split[2])

            if not("#" in line):
                data.append(float(line_split[1])+1j*float(line_split[2]))

    # Fill up negative frequencies by taking the conjugate of the positive frequencies
    data  = np.array(list(np.conj(data[::-1])) + data)
    data -= moments[0]  # Subtract 0th order moment

    norms[i] = moments[1] # The norm is the 1st order moment

    self_iw.data[:, i, i] = data[:]

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

# Assume constant error bars for every frequency (only relative weights matter)
error_bars = self_iw.copy()
error_bars.data[:, :, :] = 0.01

cont = Som(self_iw, error_bars, kind="FermionGf", norms=norms)

# Accumulate solutions
cont.accumulate(**acc_params)
# Build final solution
cont.compute_final_solution(good_chi_abs=10.0, good_chi_rel=2.0)

# Build self-energy on the real axis
n_w = 1000
self_w = GfReFreq(window=energy_window, n_points=n_w, indices=self_iw.indices)
fill_refreq(self_w, cont)

# Extract spectral function
w_mesh = np.fromiter(self_w.mesh.values(), float)
data_spectral = np.zeros((n_w,n_orb+1))
data_spectral[:, 0] = w_mesh

# For ABINIT, you need to provide -2*Im(Sigma)
for i in range(n_orb):
    data_spectral[:, i+1] = - 2 * self_w.data[:, i, i].imag

if mpi.is_master_node():

    # Write frequency grid on file
    filename = "tdmft_triqs_3o_DS3_spectralfunction_realgrid"

    with open(filename, "w") as f:

        f.write(str(len(w_mesh))+"\n")
        for w in w_mesh:
            f.write(str(w)+"\n")

    # Write spectral function on file
    filename = "tdmft_triqs_3i_DS3_Self_ra-omega_iatom0001_isppol1"

    with open(filename, "w") as f:

        for i in range(n_orb):
            for j in range(len(w_mesh)):
                f.write(str(w_mesh[j])+"\t"+str(data_spectral[j, i+1])+"\n")

    # Reconstruct Sigma(iw_n) from spectral function and compare with input
    self_rec = self_iw.copy()
    reconstruct(self_rec, cont)

    fig = Figure()

    params = {}
    params["label"]  = r"Input $\Sigma(i\omega_n)$"
    params["xlabel"] = r"$\omega_n$ (Ha)"
    params["ylabel"] = r"Im[$\Sigma(i\omega_n)$] (Ha)"

    mesh = np.fromiter(self_iw.mesh.values(), complex).imag

    fig.add_data(mesh, self_iw.data[:, 0, 0].imag, params)

    params["label"] = r"Reconstructed $\Sigma(i\omega_n)$"

    fig.add_data(mesh, self_rec.data[:, 0, 0].imag, params)
    fig.plot(legend=True)


