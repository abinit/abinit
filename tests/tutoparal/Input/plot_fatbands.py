from plot_utils import *

filename = "tdmft_triqs_1o_DS2_FATBANDS_at0001_Fe_is1_l0002"

fig = Figure()

params = {}

params["xticks"] = [0, 18, 30, 51, 76, 94]
params["xticklabels"] = [r"$\Gamma$", "N", "P", r"$\Gamma$", "H", "N"]
params["xlabel"] = r"$k$-point"
params["ylabel"] = "Energy (eV)"

bi = 6   # Initial band
bf = 20  # Last band

b_ind = 0 # Band index


with open(filename, "r") as f:

    eigen = []
    error = []

    for line in f:

        if "BAND number" in line:

            b_ind += 1
            if bi <= b_ind <= bf:
                params["yerr"] = error
                fig.add_data(range(len(eigen)), eigen, params)
            eigen = []
            error = []
            continue

        line_split = line.split()

        if len(line_split) != 3:
            continue

        eigen.append(float(line_split[1]))
        error.append(float(line_split[2]))

fig.plot()


