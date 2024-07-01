---
description: How to calculate electric fields gradients and Mossbauer Fermi contact interaction
authors: JZ
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to calculate electric fields gradients 
and Mossbauer Fermi contact interaction with the ABINIT package.

## Introduction

Because the PAW formalism provides a robust way to reconstruct the
all-electron wavefunctions in the valence space, it is suitable for
computing expectation values of observables localized even very close
to the nuclei.  Obtaining equivalent accuracy within the
norm-conserving pseudopotential framework would require very small
atomic radii for the pseudization procedure, and concomitantly high
planewave cutoff energies and lengthy calculations. There remains the
question of whether even all-electron accuracy in the valence space is
sufficient for accurate representation of observables close to the
nuclei, where conventional wisdom would suggest that deep core
polarizations might be quite significant for properties such as the
electric field gradient or Fermi contact interaction. Such concerns
turn out to be unwarranted, however, as our experience and others have
shown that the PAW formalism together with a typical chemical
valence/core separation are sufficient for accurate description of
nuclear point properties such as the electric field gradient
[[cite:Petrilli1998]], [[cite:Profeta2003]], [[cite:Zwanziger2008]],
Fermi contact interaction [[cite:Zwanziger2009]] and magnetic chemical
shielding [[cite:Pickard2001]], [[cite:Zwanziger2023]].

Both the electric field gradient and Fermi contact interaction are ground-
state observables, and their computation adds negligible time to a normal
self-consistent ground state calculation. The total charge density in the PAW
formalism contains the pseudovalence density, the nuclear ionic charges, and
the all-electron and pseudo charge densities within the PAW spheres. As done
in earlier work, the electric field gradient due to the pseudovalence density
is computed in reciprocal space, and the gradient due to the (fixed) ionic
charges is computed with an Ewald sum approach. The PAW sphere charge
densities contribute matrix elements of
$(3x_\alpha x_\beta -r^2\delta_{\alpha\beta})/r^5$, weighted by the
charge densities in each channel determined by the self-consistent
minimization procedure. This treatment [[cite:Zwanziger2008]] is more flexible
than for example assuming all bands are doubly occupied, and permits the
current approach to be used with more complex electronic and magnetic states
than just insulators.

Within ABINIT, the electric field gradient computation is invoked with
the key word [[nucefg]] (for NUClear site EFG), optionally together
with the key word [[quadmom]], at the end of a normal ground state
calculation. The PAW formalism is required, and the EFG calculation
adds only a negligible amount of time to the total. The [[nucefg]] key
word takes the values 1--3. For value 1, the electric field gradient
in atomic units and SI units (V/m$^2$) is reported, along with the
eigenvectors showing its orientation in the crystal, and the
contributions of the planewave density, the PAW on-site terms, and the
ionic contributions. When [[nucefg]] is input as 2, the electric field
gradient coupling in MHz and the asymmetry are also reported, where
the conversion is made for each atom by combining the gradient with
the nuclear quadrupole moments supplied by [[quadmom]].
Finally, [[nucefg]] input as 3 allows additional computation of
a point-charge model of the gradient, for comparison purposes. The
point charges by atom are supplied through the additional variable
[[ptcharge]]. Detailed examples of the use of ABINIT to compute EFG's
can be found in [[cite:Zwanziger2008]], [[cite:Zwanziger2009a]].

The Fermi contact interaction, which is just the electron density evaluated
precisely at the nuclear location, is an important component of the isomer
shift measured in Moessbauer spectroscopy. The isomer shift is directly
proportional to
$n_{\mathrm{abs}}(\mathbf{R})-n_{\mathrm{src}}(\mathbf{R})$,
the difference in electron
density at the absorber (the sample) and the source. Evaluating the density at a
nuclear position can be accomplished by treating
$\delta(\mathbf{r}-\mathbf{R})$ as the
observable, that is, the three-dimensional Dirac delta function centered on
the nuclear position $\mathbf{R}$. Because of the short-range nature of the delta
function, in the PAW-transformed version of the observable only matrix
elements of the on-site all-electron valence functions are required
[[cite:Zwanziger2009]], and these are evaluated from a linear fit to the first
few points of the PAW radial functions.

Within ABINIT the Fermi contact interaction is invoked by setting the key word
[[nucfc]] (for NUClear site Fermi Contact) to the value 1. When
called, the electron density at each nuclear position is reported, in atomic
units (electrons per cubic Bohr). The isomer shift as measured in Moessbauer
spectroscopy is typically reported in velocity units and is obtained from the
formula

$$\delta = (2\pi cZe^2E_\gamma /3) [n_{\mathrm{abs}}(\mathbf{R})-n_{\mathrm{src}}(\mathbf{R})]\Delta\langle r^2\rangle$$

where $c$ is the speed of light, $E_\gamma$ the $\gamma$-ray energy,
$Z$ the atomic number, $e$ the
electron charge, and $\Delta\langle r^2\rangle$ the change in the
mean square nuclear radius for
the transition. The electronic densities
$n_{\mathrm{abs}}$ and $n_{\mathrm{src}}$ refer to the absorber
and source respectively. Because of the linearity of this formula in the
density at the absorber (sample) nucleus, the only unknown
($\Delta\langle r^2\rangle$) can be
obtained by comparing the calculated values in several standards to experiment
and then the computations can be used to interpret the measurements of new
materials. In [[cite:Zwanziger2009]] it is shown how to perform such studies
on a variety of compounds.


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:nuc|The tutorial on the properties of the nuclei]] shows how to compute the electric field gradient and Mossbauer Fermi contact interaction. Prerequisite: PAW1.

