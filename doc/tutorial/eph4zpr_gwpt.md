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

ALSO GW TUTORIALS

This lesson should take about 2.0 hours.

## Formalism

In the GWPT method [[cite:Li2019]], the e-ph matrix elements are computed by replacing the
first-order change of the KS Hamiltonian due to a phonon with the variation of the $GW$ self-energy.


## Typical workflow for ZPR with GWPT

A typical workflow for ZPR-GWPT requires the same step as the ones
[needed for KS-ZPR](/tutorial/eph_intro/typical_workflow_for_zpr)
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
Note that we rely on the plasmon-pole approximation (default behaviour)

## Computing QP corrections with one-shot GW

At this point, one should stop here and perform convergence studies for [[nband]], [[ecuteps]], [[ecutsigx]] and [[ngkpt]]
to make sure that our GW results are relative stable.
Then these values can be used in our next GWPT computation.

## Computing e-ph matrix elements with GWPT

To activate the computation of the GWPT matrix elements, we need to use:

[[optdriver]] 7  # Enter EPH driver.
[[eph_task]] 17  # GWPT computation.


[[gwcomp]] 2

## Computing the ZPR with GWPT

In this section, we can finally compute the ZPR of MgO using the results produced previously.

[[gstore_gname]] "gvals"      # Use GWPT e-ph matrix elements from GSTORE (default)
[[eph_stern]] 1               # Activate Sterheimer to compute contribution given by states above nband
