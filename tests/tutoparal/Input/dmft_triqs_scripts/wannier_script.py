import numpy as np

filename = "tdmft_triqs_4o_DS2_Wannier_functions_iatom0001_001"

data_wannier = np.loadtxt(filename)
# Average over all angular momentum
data_wannier[:, 1] = np.mean(data_wannier[:, 1:], axis=1)

filename = "Fe_Wannier_0001"
with open(filename, "w") as f:
    f.write(str(data_wannier.shape[0])+"\n") # First line is number of radial points
    for x in data_wannier[:, 1]:
        f.write(str(x)+"\n")

