import matplotlib.pyplot as plt
import numpy as np

with open("tdmft_triqs_1o_DS1_DOS_AT0001","r") as f:
    lines = f.readlines()

for line in lines:
    if "Fermi energy" in line:
        fermi = float(line.split()[-1])
        break

data = np.loadtxt("tdmft_triqs_1o_DS1_DOS_AT0001")

data[:,0]-=fermi

plt.plot(data[:,0],data[:,3],lw=5,color="red")

plt.xlim(-0.3,0.3)

plt.show()
