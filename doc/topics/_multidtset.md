---
description: How to set parameters for a multi dataset calculation
authors: XG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to set parameters for a multi dataset calculation with the ABINIT package.

## Introduction

The simplest use of ABINIT corresponds to one task, with one set of data: for
example, determination of the total energy for some geometry, with some set of
plane waves and some set of k-points.

It is often needed to redo the calculations for different values of some
parameter, letting all the other things equal. As typical examples, we have
convergence studies needed to determine which cut-off energy gives the needed
accuracy. In other cases, one makes chains of calculations in order to compute
the band structure: first a self-consistent calculation of the density and
potential, then the eigenenergy computation along different lines. Similarly,
DFPT, GW or BSE calculations rely on a preliminary calculation of ground-state
wavefunctions.

For such purpose, the multi-dataset mode has been implemented.

It allows the code to treat, in one run, different sets of data, and to chain
them. The number of datasets to be treated is specified by the variable
[[ndtset]], while the indices of the datasets (by default 1, 2, 3, and so on)
can be eventually provided by the arrays [[jdtset]] or [[udtset]].

A full description of the multidataset capabilities of ABINIT can be found in
[[help:abinit#multidataset|the multidataset section of the ABINIT help file]].

A very important mechanism allows to pass information obtained from some
earlier calculation, by defining `get*` input variables. Important examples
are [[getden]] for chaining a self-consistent determination of the density
with a non-self-consistent calculation of the Kohn-Sham band structure, or
[[getwfk]] for chaining a ground-state determination of wavefunctions with a
DFPT or GW computation.

!!! tip

    |AbiPy| provides a programmatic interface to generate input files from python. 
    For futher detail, please consult the |AbinitInputNb|.


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

