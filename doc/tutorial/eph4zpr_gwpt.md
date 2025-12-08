---
authors: MG
---

# Zero-point renormalization of the band gap and temperature-dependent band gaps with GWPT

This tutorial is similar to [eph4zpr](/tutorial/eph4zpr).
Also in this lesson, we will compute the electron self-energy due to phonons, obtain the zero-point
renormalization (ZPR) of the band gap and temperature-dependent band energies within the harmonic approximation.
The main difference with respect to [eph4zpr](/tutorial/eph4zpr), is that, in this lesson,
the e-ph matrix elements are computed within the GWPT formalism [[cite:Li2019]].

It is assumed the user has already completed the two tutorials [RF1](/tutorial/rf1) and [RF2](/tutorial/rf2),
and that they are familiar with the calculation of ground state (GS) and response properties
in particular phonons, Born effective charges and the high-frequency dielectric tensor.

The user should have read the [introduction tutorial for the EPH code](/tutorial/eph_intro),
the description of the [gstore-based approach](/tutorial/gstore),
before running these examples.

Also, you are kindly invidited to read the [first GW tutorial](/tutorial/gw1) if you are not familiar with GW.
A brief description of the formalism and of the equations implemented in the
code can be found in the [[theory:mbt|GW_notes]].

This lesson should take about 2.0 hours.

## Formalism

In the GWPT method [[cite:Li2019]], the e-ph matrix elements are computed by replacing the
first-order change of the KS Hamiltonian due to a phonon with the variation of the $GW$ self-energy.


## Typical workflow for ZPR with GWPT

A typical workflow for ZPR-GWPT requires the same step as the ones
[needed for KS-ZPR](/tutorial/eph_intro#typical_workflow_for_zpr)
plus additional computations for the screened interaction $W$.

## Getting started

[TUTORIAL_README]

Before beginning, you might consider to work in a different subdirectory as for the other tutorials.
Why not create Work_eph4zpr in $ABI_TESTS/tutorespfn/Input?

```sh
cd $ABI_TESTS/tutorespfn/Input
mkdir Work_eph4zpr_gwpt
cd Work_eph4zpr_gwpt
```

In this tutorial, we prefer to focus on the use of the EPH code hence
we will be using **pre-computed** DDB and DFPT POT files to bypass the DFPT part.
We also provide a DEN.nc file to initialize the NSCF calculations
and a POT file with the GS KS potential required to solve the Sternheimer equation.

If *git* is installed on your machine, one can easily fetch the entire repository (23 MB) with:

```sh
git clone https://github.com/abinit/MgO_eph_zpr.git
```

Alternatively, use *wget*:

```sh
wget https://github.com/abinit/MgO_eph_zpr/archive/master.zip
```

or *curl*:

```sh
curl -L https://github.com/abinit/MgO_eph_zpr/archive/master.zip -o master.zip
```

or simply copy the tarball by clicking the "download button" available in the github web page,
unzip the file and rename the directory with:

```sh
unzip master.zip
mv MgO_eph_zpr-master MgO_eph_zpr_gwpt
```

!!! warning

    The directory with the precomputed files must be located in the same working directory
    in which you will be executing the tutorial and must be named `MgO_eph_zpr_gwpt`.


TODO: For the discussion on how to merge the DDB and DVDB files, I can use the link to the MgO ZPR lesson

## Computing the WFK files with empty states

[HDIAGO_README]

## Computing the screened interaction W

In this section, we use the WFK file generated in the previous section to compute the RPA polarizability and the
screened interaction $W$.
Note that here we generate a screening file with only two frequencies (default behaviour).
The SCR file will then be used to construct the plasmon-pole approximation when computing the self-energy using
[[ppmodel]].

## Computing QP corrections with one-shot GW

At this point, one should stop here and perform convergence studies for [[nband]], [[ecuteps]], [[ecutsigx]] and [[ngkpt]]
to make sure that our GW results are relative stable.
Then these values can be used in our next GWPT computation.

## Computing e-ph matrix elements with GWPT

To activate the computation of the GWPT matrix elements, we use the following two variables:

[[optdriver]] 7  # Enter EPH driver.
[[eph_task]] 17  # GWPT computation.

GWPT computations require several external files in input.
In what follows, we describe all the files step by step and discuss the connection with the parameters appearing in Eq.

The KS states are read from the WFK file via [[getwfk_filepath]].
This file defines the list of $\kk$-points in the e-ph matrix elements.
The value of [[ngkpt]], [[nshiftk]] and [[shiftk]] must be consistent with the ones used to generate the WFK file.
The number of bands in the $n'$ sum is given by [[nband]].
Clearly this value cannot be greater than the number of bands stored in the WFK file.

The screening is read from the SCR file specified with [[getscr_filepath]].
The cutoff energy in $W$ is given [[ecuteps]], while [[ecutsigx]] defines
the cutoff-energy for the exchange part of the self-energy.
Note that [[ecuteps]] cannot be larger than the value used in the screening calculation.
The SCR file defines the $\pp$-mesh for the integration over transferred momennta in Eq.
This $\pp$-mesh must be identical to, or a submesh of, the $\kk$-mesh associated with the WFK file.
No interpolation in $\pp$-space is possible at this level.

Also, [[ppmodel]] defines the kind of plasmon-pole approximation.
By default we use the Godby-Needs model
At the time of writing only ppmodel 1 and 2 are supported in the GWPT code.

The DFPT KS potentials are read from [[getdvdb_filepath]], while the DFPT KS densities
are taken from [[getdrhodb_filepath]].
These two files define the list of $\qq$-points in the GWPT matrix elements.
This $\qq$-mesh must be identical to, or a submesh of, the $\kk$-mesh associated with the WFK file.
Note that, in this case, it is possible to densify the $\qq$-mesh by using [[eph_ngqpt_fine]].
In a typical scenario, one generates a WFK file on a $\kk$-mesh much denser than the one used in the DFPT part
and then use the Fourier interpolation of the DFPT potentials to reach a $\qq$-mesh that is equal or half the $kk$-mesh.

Finally, the GWPT code needs to read the GS KS potential from the file specified with [[getpot_filepath]].
This file is produced at the end of the GS SCF cycle by setting [[prtpot]] to 1 (note that the default if 0).
The first order change of the KS wavefunctions due to an atomic perturbation is computed on-the-fly
by solving the NSCF Sterheimer equation.
There are two variables controlling the NSCF cycle:
[[nstep]] defines the maximum number of iterations while [[tolwfr]] gives the stopping criterion.

[[zcut]]

TODO
[[gwcomp]] 2

## Computing the ZPR with GWPT

In this section, we can finally compute the ZPR of MgO using the results produced previously.

[[optdriver] 7
[[eph_task]] 24

[[gstore_vname]] "gvals"      # Use GWPT e-ph matrix elements from GSTORE (default)
[[eph_stern]] 1               # Activate Sterheimer to compute contribution given by states above nband

It is important to understand that at this level of the calculations
there are few parameters that can be changed this run is essentially a post-processing
of the e-ph matrix elements stored in the GSTORE file.

The temperature mesh is defined by [[tmesh]]

[[gstore_vname]]

[[zcut]]
