
This directory contains example usages of the module ElectronPhononCoupling.

To run the examples, you must first run the Abinit calculations
in Calculations/01-LiF-dynamical/ and Calculations/02-LiF-static/.


The difference between these two calculations is that the former will split
the contributions to the second-order derivatives of the eigenvalues into
the 'active' and the 'Sternheimer' parts, while the latter sums the two
contributions in the same file. The separation into active and Sternheimer
parts is more flexible. It allows for dynamical AHC calculation, and the
computation of the spectral function. See [PRB 92, 085137 (2015)].


