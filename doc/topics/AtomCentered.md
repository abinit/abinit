---
description: How to compute atom-centered properties
authors: XG
---
<!--- This  file has been generated automatically from the corresponding _* source file. DO NOT EDIT. Edit the source file instead. -->

This page gives hints on how to compute atom-centered properties like the integral of the charge
within a sphere, or orbital magnetization within a sphere, etc,  with the ABINIT package.

## Introduction

For the purpose of understanding a material, it is interesting to be able
to integrate the density,
the spin-magnetization, or the orbital magnetization in
a sphere around an atom. Also, 
taking into account only s-, p-, or d- contributions might be desirable.

The definition of the list of atoms and the radii of the spheres is made through [[natsph]], [[iatsph]] and [[ratsph]],
and apply to all such sphere-based calculations. In the PAW case, [[ratsph]] 
is the PAW sphere cut-off from the PAW atomic dataset, but can be adjusted in the NC case.
It is even possible to define spheres that are not centered on an atom, using the
inptut variables  [[natsph_extra]], [[xredsph_extra]] and [[ratsph_extra]].

The computed integrated charges or magnetizations are printed inside the main output file of ABINIT.

For the DOS, use the intput variable [[prtdos]] with value 3 or 4, possibly m-decomposed if [[prtdosm]] is activated.
In the PAW case, use [[pawprtdos]].

To perform an atom-by-atom orbital magnetization integration inside the PAW spheres, use the prt_lorbmag input variable.

A band structure can even be represented using weights proportional to the
orbital content (so-called "Fat Bands"), in case of PAW calculation, see
[[pawfatbnd]], and related variables.

This topic is also strongly related to the two topics [[topic:ConstrainedDFT]] and [[topic:MagMom]].

## Related Input Variables

*compulsory:*

- [[abinit:iatsph]]  Index for the ATomic SPHeres of the atom-projected density-of-states
- [[abinit:natsph]]  Number of ATomic SPHeres for the atom-projected density-of-states
- [[abinit:natsph_extra]]  Number of ATomic SPHeres for the l-projected density-of-states in EXTRA set
 
*basic:*

- [[abinit:ratsph]]  Radii of the ATomic SPHere(s)
 
*useful:*

- [[abinit:pawfatbnd]]  PAW: print band structure in the FAT-BaND representation
- [[abinit:pawprtdos]]  PAW: PRinT partial DOS contributions
- prt_lorbmag  PRinT L ORBital MAGnetic moment inside PAW spheres
- [[abinit:prtdos]]  PRinT the Density Of States
- [[abinit:prtdosm]]  PRinT the Density Of States with M decomposition
- [[abinit:ratsph_extra]]  Radii of the ATomic SPHere(s) in the EXTRA set
- [[abinit:xredsph_extra]]  X(position) in REDuced coordinates of the SPHeres for dos projection in the EXTRA set
 

## Selected Input Files

*tutorial:*

- [[tests/tutorial/Input/tspin_3.abi]]
 
*v4:*

- [[tests/v4/Input/t38.abi]]
 
*v5:*

- [[tests/v5/Input/t19.abi]]
- [[tests/v5/Input/t20.abi]]
 
*v9:*

- [[tests/v9/Input/t108.abi]]
 

