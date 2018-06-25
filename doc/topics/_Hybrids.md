---
description: How to use hybrid functionals
authors: FJ
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to use hybrid functionals with the ABINIT package.

## Introduction

The Fock exchange term has been implemented in ABINIT, both in the norm-
conserving pseudopotential framework and in the PAW one. Some details about
the implementation in ABINIT can be found [[pdf:hybrids-2017.pdf|here]].  
For an ABINIT user, to make a calculation of Fock exchange:  
\- do a first GGA dataset for the ground state  
\- do a second dataset for the Fock calculation choosing [[ixc]]=40-42 (HF,
PBE0, PBE0-1/3), -402 (B3LYP-Libxc), -406 (PBE0-Libxc), -456 (PBE0-1/3 Libxc),
or -427, -428 (HSE03, HSE06).

  
The energy, forces and stresses are available in the norm-conserving and PAW
frameworks.  
A one-shot G0W0 calculation can follow, only in the norm-conserving case at present.

!!! warning 

    Use [[istwfk]]=1, [[iscf]]=2, [[paral_kgb]]=0, [[paral_atom]]=0.  
    The efficiency of the calculation is not optimal. Work is in progress
    concerning this point.


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

