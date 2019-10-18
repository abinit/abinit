---
authors: FJ, XG
---

# ABINIT features  

## An overview of ABINIT settings and features, for beginners and more experienced users.

This document gives an overview of the features implemented in the ABINIT
package, grouped in different topics, for the beginner as well as more experienced user. 
It might answer the question "How to ... with ABINIT ?", to some extent.  
It also gives a synthetic view on the needed settings.

## 1 Foreword
  
Documenting the features of a large scientific code is a complex task. The
present list of features refers to different "topics". Each topic has a
dedicated page, which should be quick to read, unlike the
[[lesson:index|lessons of the tutorial]], each of which is usually at least
one hour work. Many of the topics make the link between a physical property or
quantity (including bibliographical references for the theory, and sometimes
pointing to published work using this feature) and the way it is to be
computed with ABINIT (e.g. corresponding input variable, example input files,
and possibly tutorial(s)), and the associated post-processing tools (OPTIC, ANADDB, MULTIBINIT ...).

This "topic"-based approach might be used by the beginner to get a broad
overview of the capabilities of ABINIT applications as well as to the more expert user to
quickly find the way to compute some existing quantity, or to remember which
input variable is useful or mandatory for the calculation of a given property.
Some topics are rather "input"-oriented (e.g. how to specify the atomic
geometry, the occupation numbers, etc), other are more "property"-oriented
(e.g. how to compute the elastic constants, the temperature-dependence of the
electronic structure, etc), other are related to proper/better usage of the code.

Care is taken not to duplicate existing more complete documentation in ABINIT,
but to point to it if appropriate. Not all the ABINIT documentation is covered
by the Web-accessible documents, there are still a few unlinked documents in
the subdirectories of ~abinit/doc (work is in progress to make it all available).
Discussions on the [ABINIT forum](https://forum.abinit.org) might also allow to get information.

## 2 ABINIT specifications for static DFT calculations
  
The topic [[topic:GSintroduction|Building an input file]] briefly explains the
content of an ABINIT input file. The following topics go more into the
details, however restricting to **static DFT calculations**, without doing
anything fancy with the exchange-correlation functionals. Going beyond these
is left for other sections (Section 3 and beyond). In
particular, for any accurate electronic properties, e.g. correct band
structure, optical response, or for strongly correlated electrons, please go
beyond the present sec. 2. Also, topics related to global control parameters,
that apply, generally speaking, to all types of calculations are explained later.

### 2.1 Settings for atoms: cell, atoms, atomic positions, and symmetries
  
  1. [[topic:UnitCell|Unit cell]]
  2. [[topic:AtomTypes|Types of atoms and alchemy]]
  3. [[topic:crystal|Crystalline structure and symmetries]]
  4. [[topic:SmartSymm|Smart symmetrizer]]
  5. [[topic:AtomManipulator|Atom manipulator]] (advanced topic)

### 2.2 Physical settings for electrons: XC functionals, atomic/pseudo potentials, metal/insulator, spin, Coulomb interaction ...
  
  1. [[topic:xc|Overview of available exchange and correlation functionals]]
  2. [[topic:Hybrids|Hybrid functionals]]
  3. [[topic:vdw|Van der Waals functionals]]
  4. [[topic:RPACorrEn|Correlation energy within RPA]]
  5. [[topic:PseudosPAW|Pseudopotentials and PAW atomic datasets]]
  6. [[topic:BandOcc|Bands and occupation numbers for metals and insulators]]
  7. [[topic:spinpolarisation|Spin-polarised systems and spin-orbit coupling]]
  8. [[topic:Coulomb|Coulomb interaction and charged cells]]

### 2.3 Numerical settings for electrons: basis set, planewaves and real space sampling, Brillouin zone sampling ...
  
  1. [[topic:Planewaves|Planewaves and real space sampling]]
  2. [[topic:PAW|PAW special settings]]
  3. [[topic:Wavelets|Wavelets in ABINIT]]
  4. [[topic:k-points|Wavevector sampling (k point grid)]]

### 2.4 SCF algorithms, tuning and stopping criteria
  
  1. [[topic:SCFAlgorithms|SCF algorithms]]
  2. [[topic:SCFControl|SCF control, tolerances and stopping criteria]]
  3. [[topic:ForcesStresses|Forces and stresses]]
  4. [[topic:TuningSpeed|Tuning the speed of the calculation]]
  5. [[topic:Recursion|Recursion methods and orbital free calculations]] (not in production)

### 2.5 Added electric/magnetic field, other artificial constraints/modifications, and related properties ...
  
  1. [[topic:Berry|Electric polarization and finite electric field]]
  2. [[topic:MagField|External magnetic field]]
  3. [[topic:ConstrainedDFT|Constrained Density-Functional Theory]]
  4. [[topic:MagMom|Constrained atomic magnetic moment]]
  5. [[topic:EFG|Electric field gradients and Mossbauer Fermi contact interaction]]
  6. [[topic:Artificial|Artificial modifications of the system]]

## 3 Global control parameters: flow, parallelism, output files, output content, timing and memory control ...
  
  1. [[topic:multidtset|Multi-dataset calculations]]
  2. [[topic:parallelism|Parallelism and ABINIT]]
  3. [[topic:printing|Printing files]]
  4. [[topic:Output|Tuning the output content in different files]]
  5. [[topic:Control|Time and memory control]]

## 4 Molecular dynamics, geometry optimization, transition paths
  
  1. [[topic:GeoOpt|Geometry optimization]]
  2. [[topic:MolecularDynamics|Molecular dynamics]]
  3. [[topic:TransPath|Transition path searches: NEB and string method]]
  4. [[topic:GeoConstraints|Geometry constraints]]
  5. [[topic:PIMD|Path-integral molecular dynamics (PIMD)]]
  5. [[topic:CrossingBarriers|Crossing barrier search, linear combination of constrained DFT energies and ensemble DFT]]
  6. [[topic:LOTF|Learn-on-the-flight (LOTF)]] (not in production)

## 5 Correlated electrons
  
When correlated electrons are to be considered (in most cases, when *d* and *f*
orbitals play an active role), it is necessary to go beyond the standard DFT
framework. ABINIT enables the following possibilities:

  1. [[topic:Hybrids|Hybrid functionals]]
  2. [[topic:DFT+U|DFT+U approximation]]
  3. [[topic:DMFT|Dynamical Mean Field Theory (DMFT)]]
  4. [[topic:CRPA|Calculation of the effective Coulomb interaction]]

## 6 Adiabatic response properties (phonons, low-frequency dielectric, Raman, elasticity, temperature dependence ...)
  
Many properties can be obtained in the approximation that the electrons **stay
in their ground state** (adiabatic responses). The poweful Density-Functional
Perturbation Theory (DFPT) framework allows ABINIT to address directly all
such properties in the case that are connected to derivatives of the total
energy with respect to some perturbation. This includes all dynamical effects
due to phonons and their coupling, thus also temperature-dependent properties due to phonons.

  1. [[topic:DFPT|Generalities about DFPT]] 
  2. [[topic:q-points|Wavevectors for phonons (q-points)]] 
  3. [[topic:Phonons|Vibrational and dielectric properties 
     (phonon frequencies and modes, IR and Raman spectra, Born effective charges)]]
  4. [[topic:PhononBands|Phonon bands and DOS, interatomic force constants, sound velocity]]
  5. [[topic:Temperature|Temperature dependent properties (free energy, entropy, specific heat, 
     atomic temperature factors, thermal expansion)]]
  6. [[topic:Elastic|Elasticity and piezoelectricity]]
  7. [[topic:nonlinear|Raman intensities and electro-optic properties]]
  8. [[topic:ElPhonInt|Electron-phonon interaction]]
  9. [[topic:PhononWidth|Phonon linewidth due to the electron-phonon interaction]]
  10. [[topic:ElPhonTransport|Electronic transport properties from electron-phonon interaction 
     (resistivity, superconductivity, thermal)]]
  11. [[topic:TDepES|Temperature dependence of the electronic structure from electron-phonon interaction]]
  12. [[topic:ConstrainedPol|Constrained polarization geometry optimization]] (advanced topic)


## 7 Excited state calculations, and frequency-dependent electronic and optical properties
  
Excited-state calculations and frequency-dependent properties (for frequencies
that are non-negligible with respect to the electronic gap), can be addressed
by a variety of methodologies, usually trading accuracy for speed. At the
lowest level, one encounters the independent-particle approximation, building
upon some previously obtained band structure (e.g. Kohn-Sham band structure
from DFT). For charged excitations, allowing to obtain a quasiparticle band
structure without the well-known DFT band gap problem, one has to resort to
(costly) GW calculations. For neutral excitations (i.e. optical), the (costly)
Bethe-Salpeter approach is the most accurate formalism presently available.
TDDFT and Δ-SCF calculations are cheaper but will work well for molecules and
for isolated defects in a solid, not for e.g. correcting the band gap.

  1. [[topic:Optic|Linear and non-linear optical properties in the independent-particle approximation]]
  2. [[topic:FrequencyMeshMBPT|Definition of frequency meshes for Many-body perturbation theory]]
  3. [[topic:Susceptibility|Frequency-dependent susceptibility matrix, and related dielectric matrix and screened Coulomb interaction]]
  4. [[topic:SelfEnergy|Electronic self-energy]]
  5. [[topic:GW|GW calculations for accurate band structure, including self-consistency]]
  6. [[topic:BSE|Bethe-Salpeter calculations for accurate optical properties]]
  7. [[topic:TDDFT|TDDFT calculations]]
  8. [[topic:DeltaSCF|DeltaSCF calculations]]
  9. [[topic:RandStopPow|Random electronic stopping power]]
  10. [[topic:GWls|GW- Lanczos-Sternheimer method]] (not in production)

## 8 Second-principles calculations with MULTIBINIT: handling millions of atoms with first-principles accuracy

By constructing model Hamiltonians whose linear and selected non-linear characteristics
are identical to those from first-principles calculations, and simulating millions
of atoms with these model Hamiltonians, one can study phase transitions, polarization boundaries,  
and other properties for large-scale systems that cannot be reached from first-principles algorithms
implemented in ABINIT and most DFT codes. Even with respect to linear-scaling codes, the prefactor
is much smaller. 
This is implemented in the MULTIBINIT application.

  1. [[topic:LatticeModel|Lattice model at the harmonic level]]
  2. [[topic:FitProcess|FitProcess]]
  3. [[topic:BoundingProcess|BoundingProcess]]
  4. [[topic:DynamicsMultibinit|DynamicsMultibinit]]

## 9 Electronic properties and analysis tools (DOS, STM, Wannier, band plotting and interpolating...)

Many properties are directly deduced from the knowledge of the electronic
wavefunctions, eigenenergies, density, potential, etc. Some necessitates
additional tuning parameters or are linked to postprocessing tools and are
described in the following topics. Some others are actually activated through
a single printing parameter, such as the Electron Localization Function (ELF -
see [[prtelf]]). See the list of "printing" input variables in [[topic:printing]].

  1. [[topic:ElecBandStructure|Electronic band structure and related topics]]
  2. [[topic:ElecDOS|Electronic DOS and related topics]]
  3. [[topic:EffectiveMass|Effective mass calculations]]
  4. [[topic:Unfolding|Unfolding supercell band structures]]
  5. [[topic:DensityPotential|Manipulating the density and potential]]
  6. [[topic:Macroave|Macroscopic average of density and potential]]
  7. [[topic:STM|Scanning Tunneling Microscopy map]]
  8. [[topic:Wannier|Wannier functions]]
  9. [[topic:Bader|Bader Atom-In-Molecule analysis]]

## 10 Other physical properties (e.g. positron)
  
  1. [[topic:positron|Positron calculations]]
  2. [[topic:LDAminushalf|The LDA-1/2 approach]] 

## 11 Analysis/postprocessing tools
  
  1. [[topic:Abipy|Abipy - ABINIT swiss knife]]
  2. [[topic:APPA|Abinit Post-Processor Application (APPA), for molecular-dynamics trajectory analysis]]
  3. [[topic:Band2eps|Band2eps for phonon dispersion curves]]
  4. [[topic:Tdep|Temperature Dependent Effective Potential, for thermodynamical properties]]

## 12 Miscellaneous topics
  
  1. [[topic:Verification|Verification of the implementation]]
  2. [[topic:PortabilityNonRegression|Portability and non-regression tests]]
  3. [[topic:Git|Git, gitlab and github for the ABINIT project]]
  4. [[topic:Dev|Miscellaneous for developers]]
