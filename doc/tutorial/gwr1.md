---
authors: MG
---

# First tutorial on GWR (GW in real-space and imaginary time)

## The quasi-particle band structure of Silicon in the GW approximation.

This tutorial aims at showing how to calculate self-energy corrections to the
DFT Kohn-Sham (KS) eigenvalues in the one-shot GW approximation using the GWR code

The user should be familiarized with the four basic tutorials of ABINIT,
see the [tutorial home page](/tutorial),
and is strongly encouraged to read the [introduction to the GWR code](/tutorial/gwr_intro) before running these examples.

This tutorial should take about 2 hours.

[TUTORIAL_README]

### Generation of the WFK file with empty states

*Before beginning, you might consider creating a different subdirectory to work in.
Why not create Work_gwr ?*

The file *tgwr_1.abi* is the input file for the first step
(GS run followed by direct diagonalization of the KS Hamiltonian).
Copy it to the working directory with:

```sh
mkdir Work_gwr
cd Work_gwr
cp ../tgwr_1.abi .
```

{% dialog tests/tutorial/Input/tgwr_1.abi %}

This step might be quite time-consuming so you may want to immediately start the job in background with:

```sh
abinit tgwr_1.abi > tgwr_1.log 2> err &
```


The calculation is done for silicon.
The first dataset produces the density file that used to compute the band structure in the second dataset.
We will also use this DEN file to generate the WFK with empty states in the next section.


!!! important

    gwr_task2 "HDIAGO_FULL"


The direct diagonalization is MPI-parallelized across three different levels:
collinear spin $\sigma$, k-points and scalapack distribution of the $H^\sigma_\kk(\bg,\bg')$ matrix.
Abinit will try to find an optimal distribution of the workload at runtime yet there are a couple
of things worth keeping in mind when choosing the number of MPI processes for this step.
Ideally the total number of cores should be a multiple of [[nkpt]] * [[nsppol]] to avoid load imbalance.

### QP corrections with GWR

{% dialog tests/tutorial/Input/tgwr_2.abi %}

The k-points for the QP corrections can be specified in different ways.

Explicitly via:

[[nkptgw]], [[kptgw]] and [[bdgw]]

Implicitly via [[gw_qprange]]

For the spectral function

[[nfreqsp]], [[freqspmax]]

Notes on the MPI parallelization

The GWR code employs [[gwr_np_kgts]].
Ideally the total number of MPI processes should be a multiple of [[gwr_ntau]] * [[nsppol]].
