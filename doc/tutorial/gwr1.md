---
authors: MG
---

# First tutorial on GWR (GW in real-space and imaginary time)

## The quasi-particle band structure of Silicon in the GW approximation.

This tutorial aims at showing how to calculate self-energy corrections to the
DFT Kohn-Sham (KS) eigenvalues in the on-shot GW approximation using the GWR code

The user should be familiarized with the four basic tutorials of ABINIT,
see the [tutorial home page](/tutorial).

This tutorial should take about 2 hours.

WARNING : THIS TUTORIAL IS WORK IN PROGRESS ! IT IS NOT YET COMPLETE ...

[TUTORIAL_README]

### Ground-state computation

The first dataset produces the density file that used to compute the band structure in the second dataset.
We will also use this DEN file to generate the WFK with empty states in the next section.

### Generation of the WFK file with empty states

```
###########################################
# Dataset 2: Direct diago with empty states
###########################################
optdriver2  6        # Activate GWR code
gwr_task2 "HDIAGO"   # Direct diagonalization
getden2       1
nband2     40        # Number of (occ + empty) bands
```


!!! important

    gwr_task2 "HDIAGO_FULL"   


The direct diagonalization is MPI-parallelized across three different levels: 
collinear spin $\sigma$, k-points and scalapack distribution of the $H^\sigma_\kk(\bg,\bg')$ matrix.
Abinit will try to find an optimal distribution of the workload at runtime yet there are a couple 
of things worth keeping in mind when choosing the number of MPI processes for this step.
Ideally the total number of cores should be a multiple of [[nkpt]] * [[nsppol]] to avoid load imbalance.


### QP corrections with GWR

The k-points for the QP corrections can be specified in different ways.

Explicitly via:

[[nkptgw]], [[kptgw]] and [[bdgw]]

Implicitly via [[gw_qprange]]

For the spectral function

[[nfreqsp]], [[freqspmax]]

Notes on the MPI parallelization

The GWR code employs [[gwr_np_kgts]].
Ideally the total number of MPI processes should be a multiple of [[gwr_ntau]] * [[nsppol]].
