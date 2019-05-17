# The tdep utility  

The Temperature Dependent Effective Potential (TDEP) method
has been developped by O. Hellman *et al.* [[cite:Hellman2011]],
[[cite:Hellman2013]], [[cite:Hellman2013a]] in 2011 and the |aTDEP| implementation
in ABINIT has been performed and used for the first time in 2015 by
J. Bouchet and F. Bottin [[cite:Bouchet2015]], [[cite:Bouchet2017]].

This manual can be found as a pdf file: [[pdf:TDEP_Guide| TDEP guide]]

## Prerequisite and theory

The approach used in this code is detailed in a publication dedicated to the development
of all formula (see [[pdf:TDEP_Paper|TDEP paper]]). We strongly encourage all the users to carefully read
this paper before beginning. All the vibrational, elastic and thermodynamic
quantities computed by |aTDEP| are
presented with the same writing conventions as the ones used in the output files of |aTDEP|.
In the same manner, a comprehensive understanding of some ABINIT basic variables is also required
in order to fill the input file and read the output file of |aTDEP|.

In addition, this paper is also useful to understand
the limitations and convergences which are inherent to the present method.
These particular points are sometimes discussed in the
article, with some references and illustrating examples.

## The ABINIT computation

To run |aTDEP|, a preliminary
ABINIT simulation is needed. This one could be a molecular dynamic trajectory
or a set of "ground state" calculations on specific configurations (representative of a given thermodynamic state).
After that, all the configurations have to be merged:
(i) in a single *NetCDF* file `HIST.nc` or (ii) in three separated *ASCII* files `fcart.dat`, `xred.dat` and
`etot.dat` (forces in cartesian coordinates, positions in reduced coordinates, total energies in Ha),
as they are written in the output file of ABINIT. In the later case, the 3 files can be built easily
by concatenating in each one all the time steps or configurations (using `agrep` shell instruction, for example).

## The |aTDEP| computation

In a same manner as performed for ABINIT, the use of |aTDEP| is quite simple. 
One has just to execute `tdep` as follows:

```sh
    tdep < input.files > log
```

with the `input.files` file containing 3 lines. The first one defines the input
file, the second one is the *NetCDF* file (if present, see above) and the third one
defines the root of all the output files:

    input.in
    HIST.nc
    output

The detection of the `HIST.nc` file is performed at the beginning; so, if this
one is absent, the code will automatically search the 3 `ASCII.dat` files.

## The input files

An example of a |aTDEP|  calculation (in the special case where the *NetCDF* file `HIST.nc` is employed)
can be found in [[test:v8_37]]. The 2 input files are
given in the `tests/v8/Input` directory.  
Let us describe briefly this [[test:v8_37]] file:

{% dialog tests/v8/Input/t37.in %}

The input file format is fixed. So:

1. This file begins with a `NormalMode` or `DebugMode` keyword and finishes with `TheEnd` (all the lines after are not read).
2. All the lines between `# Unit cell definition` and `# Optional inputs` are fixed.
3. Between `# Optional inputs` and `TheEnd`, the format is free.

More details:

* The section `# Unit cell definition` defines the bravais lattice [[brav@tdep|brav]]
  (here, a simple cubic), the number of atoms in the unit cell [[natom_unitcell@tdep|natom_unitcell]]
  (here, 5 atoms), their reduced coordinates in the unit cell [[xred_unitcell@tdep|xred_unitcell]]
  (here, a perovskite) and the type of atoms in the unit cell [[typat_unitcell@tdep|typat_unitcell]]
  (here, one atom A, one atom B and 3 atoms C).
* The section `# Supercell definition` defines the multiplicity of the
  supercell with respect to the unit cell multiplicity (here, it is a simple
  2x2x2 multiplication of the unit cell) and the temperature of the system
  temperature(here, 495.05 K).
* The section `# Computation details` defines the range [[nstep_max@tdep|nstep_max]]...[[nstep_min@tdep|nstep_min]]
  of time steps or configurations (here, 100 time steps) and the
  cutoff radius for the pair interactions [[rcut@tdep|Rcut]] (here, all the interaction pairs
  with a bond length larger than 7.426 bohr will not be considered).
* The section `# Optional inputs` can define a large number of optional
  keywords (here [[ngqpt2@tdep|ngqpt2]] defining the q-point grid for the vDOS integration
  is set to 2 2 2 in order to have a test sufficiently fast, which means that
  all the thermodynamic quantities have no sense.)
All the input variables are defined in the `tdep` section of the input variables description.
Note that some input variables, not defined in the `input.in` file, are obtained
from the `HIST.nc` file. In particular, the features of the supercell.

<sub><sup>TODO: Explain the extra input variables when the 3 ASCII files are employed.</sup></sub>

## The output files

A large number of output files are obtained after an execution of |aTDEP|.

{% dialog tests/v8/Refs/t37.out %}
{% dialog tests/v8/Refs/t37omega.dat %}
{% dialog tests/v8/Refs/t37thermo.dat %}

1. `*.out` is the main output file. It includes an echo of the input variables,
   some intermediary results, the definition of the various shells of interaction,
   the second order IFCs for all the atoms in each shell, the elastic constants and moduli,
   the energy of the model...
2. `*omega.dat` contains the dispersion of phonon frequencies (in meV) along a path in the Brillouin Zone.
3. `*thermo.dat` lists all the thermodynamic quantities obtained by considering
   the system as a quantum harmonic crystal: specific heat, vibrational
   energy, entropy and free energy. It also gives all these contributions as a
   function of temperature in the harmonic approximation.
4. `sym.dat` details all the symmetry operations of the bravais lattice,
5. `qpt.dat` defines the q-point grid used to compute the phonon frequencies
   contained in the `omega.dat` file.
6. `xredaverage.xyz` includes the ideal and average positions in the supercell.
7. `Indym*.dat` contain all the symmetry relations between one or two atoms
   in the unit cell or the supercell.
8. `vdos.dat` displays the vibrational density of states (in meV).
9. `dij.dat` lists the dynamical matrices for a particular set of q-points.
10. `etotMDvsTDEP2.dat` compares the MD trajectory with the one computed
    using the second order IFCs (these ones must be superimposed, as much
    as possible).
11. `fcartMDvsTDEP2.dat` plots the MD forces wrt the forces computed using
    the second order IFCs (the cloud of points must be closer to the first bisector).
12. `eigenvectors.dat` lists all the eigenvectors for a particular set of q-points.
13. `nbcoeff-phij.dat` shows how the number of IFC coefficients are reduced (for each shell and each symmetry).
14. ...
