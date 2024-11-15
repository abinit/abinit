---
description: How to calculate NMR chemical shieldings
authors: JZ
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to calculate NMR chemical shieldings
with the ABINIT package.

## Introduction

A key observable measured in NMR spectroscopy is the chemical
shielding, which is usually thought of as the shielding of the
external magnetic field caused by the electrons in the sample, the
effect of which is to reduce slightly the Zeeman splitting of the
energy levels of the nuclear magnetic dipole from that of a bare
nucleus [[cite:Slichter1978]]. More precisely, the external magnetic
field induces an orbital electronic current, which itself generates a
small secondary magnetic field opposite to the bare field. One
approach to computing this effect in a DFT context is provided by the
Gauge Including Projector Augmented Wave formalism (GIPAW) to compute the
induced current, and from that, the effective induced field
[[cite:Pickard2001]].

From an energetic perspective, though, chemical shielding is just the
effect on the total energy of both a nuclear magnetic dipole and an
external magnetic field:

$$ \sigma_{ij} = \frac{\partial^2 E}{\partial \mathbf{m}_i\partial \mathbf{B}_j} $$

for nuclear dipole $\mathbf{m}$ and magnetic field $\mathbf{B}$. From
this perspective, the shielding is either the induced field acting on
the bare dipole, or the induced dipole acting on the bare field
[[cite:Thonhauser2009]]. The latter approach is implemented in ABINIT,
where the first order energy

$$ E^{(1)} = \frac{\partial E}{\partial\mathbf{B}_j} $$

is computed in the presence of a small nuclear magnetic dipole [[cite:Zwanziger2023]].

## Execution

While not a true response function, $E^{(1)}$ turns out to dependn on
both the ground state wavefunctions and the DDK wavefunctions,
$|\partial u_{n\mathbf{k}}/\partial k\rangle$. Thus, to compute the
effect in ABINIT, first, ground state wavefunctions are computed in
the presence of a small nuclear magnetic dipole moment. The moment is
described by the variable [[nucdipmom]], which is input as a set of
3-vectors, one for each atom in the unit cell. Most of these will be
zero, and typically just "1 0 0" or "0 1 0" or "0 0 1" will be input
for the single atom one wishes to study. The triple of numbers are the
Cartesian directions, so to compute the full spatial dependence, three
separate calculations will be carried out.

Once the ground state is computed with a dipole on the atom and in the
direction of interest, a DDK calculation is carried out, with
[[rfddk]] 1 and [[rfdir]] 1 1 1, again with a dipole imposed with
[[nucdipmom]] as in ground state. For this calculation in addition,
set [[orbmag]] 2 to initiate computation of the orbital magnetization
and hence shielding at the end of the DDK calculation.

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:nuc|The tutorial on the properties of the nuclei]] shows how to compute the NMR chemical shielding.Prerequisite: PAW1, DFPT1.

