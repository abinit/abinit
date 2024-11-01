import numpy as np

nkpt = 168
nband = 38

eigs = np.zeros((nband, nkpt))

print(nband/8.)
nrows_eig = np.ceil(nband/8.)
print(nrows_eig)

with open("tw90_6_2o_DS4_EIG", "r") as fopen:
    for l, line in enumerate(fopen.readlines()):
        print(l, line)
        if not l:
            continue
        if (l - 1) % (nrows_eig + 1) == 0:
            continue
        else:
            row = int(np.floor((l - 1) % (nrows_eig + 1)) - 1)
            kpt = int(np.floor((l - 1) / (nrows_eig + 1)))
            print(row, kpt)
            for e, eig in enumerate(line.split()):
                eigs[row*8+e, kpt] = float(eig)

with open('band_struct.dat', 'w') as fopen:
    for b in range(nband):
        for k in range(nkpt):
            fopen.write("%.10f  %.10f\n" % (k, eigs[b, k]))
        fopen.write("\n")

for i in range(10):
    print(eigs[i])
