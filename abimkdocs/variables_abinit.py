# coding: utf-8
from __future__ import print_function, division, unicode_literals, absolute_import

executable = "abinit"

from abimkdocs.variables import ValueWithUnit, MultipleValue, Range, ValueWithConditions
#from abipy.abio.abivar_database.variables import ValueWithUnit, MultipleValue, Range, ValueWithConditions
Variable = dict

variables = [
Variable(
    abivarname="accuracy",
    varset="basic",
    vartype="integer",
    topics=['Planewaves_basic', 'SCFControl_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="ACCURACY",
    added_in_version="before_v9",
    text=r"""
Allows to tune the accuracy of a calculation by setting automatically the
variables according to the following table:

accuracy         | 1         | 2          | 3            | 4            | 5         | 6
---              |---        |---         |---           |---           |---        |---
[[ecut]]         | E_min     | E_med      | E_med        | E_max        | E_max     | E_max
[[pawecutdg]]    | ecut      | ecut       | 1.2 * ecut   | 1.5 * ecut   | 2 * ecut  | 2 * ecut
[[fband]]        | 0.5       | 0.5        | 0.5          | 0.5          | 0.75      | 0.75
[[boxcutmin]]    | 1.5       | 1.8        | 1.8          | 2.0          | 2.0       | 2.0
[[bxctmindg]]    | 1.5       | 1.8        | 1.8          | 2.0          | 2.0       | 2.0
[[pawxcdev]]     | 1         | 1          | 1            | 1            | 2         | 2
[[pawmixdg]]     | 0         | 0          | 0            | 0            | 1         | 1
[[pawovlp]]      | 10        | 7          | 7            | 5            | 5         | 5
[[pawnhatxc]]    | 0         | 1          | 1            | 1            | 1         | 1
[[tolvrs]]       | 1.0d-3    | 1.0d-5     | 1.0d-7       | 1.0d-9       | 1.0d-10   | 1.0d-12
[[tolmxf]]       | 1.0d-3    | 5.0d-4     | 1.0d-4       | 5.0d-5       | 1.0d-6    | 1.0d-6
[[optforces]]    | 1         | 1          | 2            | 2            | 2         | 2
[[timopt]]       | 0         | 0          | 1            | 1            | 1         | 1
[[npulayit]]     | 4         | 7          | 7            | 7            | 15        | 15
[[nstep]]        | 30        | 30         | 30           | 30           | 50        | 50
[[prteig]]       | 0         | 0          | 1            | 1            | 1         | 1
[[prtden]]       | 0         | 0          | 1            | 1            | 1         | 1


*accuracy* = 4 corresponds to the default tuning of ABINIT. It is already a very accurate tuning.
For a parallel calculation, [[timopt]] is enforced to be 0.
E_min, E_med and E_max may be read from the pseudopotential file (available
only for XML PAW atomic data files). If E_min, E_med and E_max are not given
in the pseudopotential file, [[ecut]] must be given in the input file and E_max=E_med=E_max=ecut.
If the user wants to modify one of the input variable automatically tuned by *accuracy*,
they must put it in the input file. The other input variables automatically tuned
by *accuracy* will not be affected.
*accuracy* = 0 means that this input variable is deactivated.
""",
),

Variable(
    abivarname="acell",
    varset="basic",
    vartype="real",
    topics=['UnitCell_basic'],
    dimensions=[3],
    defaultval=MultipleValue(number=3, value=1),
    mnemonics="CELL lattice vector scaling",
    characteristics=['[[EVOLVING]]', '[[LENGTH]]'],
    commentdims="represented internally as acell(3,[[nimage]])",
    added_in_version="before_v9",
    text=r"""
Gives the length scales by which dimensionless primitive translations ([[rprim]])
are to be multiplied. By default, given in Bohr atomic units (1
Bohr=0.5291772108 Angstroms), although Angstrom can be specified, if
preferred, since *acell* has the [[LENGTH]] characteristics. See further
description of *acell* related to the [[rprim]] input variable, the
[[scalecart]] input variable, and the associated internal [[rprimd]] input variable.

Note that *acell* is **NOT** the length of the conventional orthogonal basis
vectors, but the scaling factors of the primitive vectors. Use [[scalecart]]
to scale the cartesian coordinates.
""",
),

Variable(
    abivarname="adpimd",
    varset="rlx",
    vartype="integer",
    topics=['PIMD_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="ADiabatic Path-Integral Molecular Dynamics",
    requires="[[imgmov]] == 9 or [[imgmov]] == 13",
    added_in_version="before_v9",
    text=r"""
Controls whether adiabatic Path-Integral Molecular Dynamics is performed or not.
The corresponding adiabaticity parameter is given by [[adpimd_gamma]].

If equal to 0, no adiabatic Path-Integral Molecular Dynamics (standard PIMD) is performed.
If equal to 1, adiabatic Path-Integral Molecular Dynamics is activated.
Only relevant with [[pitransform]] = 1 (normal mode transformation). In that case,

- the mass associated with the zero-frequency mode is the true mass [[amu]],
- the mass associated to the other higher frequency modes of the polymer
chains is equal to the normal mode mass divided by [[adpimd_gamma]] (adiabaticity parameter),
- the equation of motion on the zero-frequency mode is not thermostated.

NOT YET USABLE
""",
),

Variable(
    abivarname="adpimd_gamma",
    varset="rlx",
    vartype="real",
    topics=['PIMD_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="ADiabatic Path-Integral Molecular Dynamics: GAMMA factor",
    requires="[[adpimd]] == 1 and [[imgmov]] in [9,13]",
    added_in_version="before_v9",
    text=r"""
Adiabaticity parameter to be used in adiabatic Path-Integral Molecular Dynamics.
NOT YET USABLE
""",
),

Variable(
    abivarname="algalch",
    varset="gstate",
    vartype="integer",
    topics=['AtomTypes_expert'],
    dimensions=['[[ntypalch]]'],
    defaultval=MultipleValue(number='[[ntypalch]]', value=1),
    mnemonics="ALGorithm for generating ALCHemical pseudopotentials",
    added_in_version="before_v9",
    text=r"""
Used for the generation of alchemical pseudopotentials, that is, when
[[ntypalch]] is non-zero.

Give the algorithm to be used to generate the [[ntypalch]] alchemical
potentials from the different [[npspalch]] pseudopotentials dedicated to this use.

Presently, **algalch** can only have the value 1, that is:

  * the local potentials are mixed, thanks to the [[mixalch]] mixing coefficients
  * the form factors of the non-local projectors are all preserved,
    and all considered to generate the alchemical potential
  * the scalar coefficients of the non-local projectors are multiplied by the
    proportion of the corresponding type of atom that is present in [[mixalch]]
  * the characteristic radius for the core charge is a linear combination of the characteristic radii
    of the core charges, build with the [[mixalch]] mixing coefficients
  * the core charge function $f(r/r_c)$ is a linear combination of the core charge functions,
    build with the [[mixalch]] mixing coefficients

Later, other algorithms for the mixing might be included.

!!! important

    Note that alchemical mixing cannot be used with PAW.
""",
),

Variable(
    abivarname="amu",
    varset="rlx",
    vartype="real",
    topics=['PIMD_useful', 'Phonons_useful', 'AtomTypes_basic', 'Artificial_basic'],
    dimensions=['[[ntypat]]'],
    mnemonics="Atomic Mass Units",
    characteristics=['[[EVOLVING]]'],
    commentdefault="provided by a database of atomic masses.",
    added_in_version="before_v9",
    text=r"""
Gives the masses in atomic mass units for each kind of atom in the input cell. These
masses are used in performing molecular dynamical atomic motion if
[[ionmov]] = 1, 6, 7 or 8. They are also used in phonon calculations during the
diagonalization of the dynamical matrix. Note that one may set all masses to 1
for certain cases in which merely structural relaxation is desired and not
actual molecular dynamics.

Using the recommended values of [[cite:Martin1987]], 1 atomic mass unit = 1.6605402e-27 kg. In this
unit the mass of Carbon 12 is exactly 12.

A database of atomic masses is provided which provides the default values. Note that the
default database uses mixed isotope masses (for Carbon the natural occurrence
of Carbon 13 is taken into account). The values are those recommended by the
commission on Atomic Weights and Isotopic Abundances, Inorganic Chemistry
Division, IUPAC [[cite:Martin1987]]. For Tc, Pm, Po to
Ac, Pa and beyond U, none of the isotopes have a half-life greater than 3.0d10
years, and the values provided in the database do not come from that source.

For alchemical pseudoatoms, the masses of the constituents atoms are mixed,
according to the alchemical mixing coefficients [[mixalch]]

In most cases, the use of *amu* will be as a static (non-evolving) variable.
However, the possibility to have different values of *amu* for different
images has been coded. A population of cells with different atomic
characteristics can thus be considered, and can be made to evolve, e.g. with a
genetic algorithm (not coded in v7.0.0 though).
""",
),

Variable(
    abivarname="angdeg",
    varset="basic",
    vartype="real",
    topics=['UnitCell_useful'],
    dimensions=[3],
    mnemonics="ANGles in DEGrees",
    characteristics=['[[INPUT_ONLY]]'],
    commentdefault="deduced from '[[rprim]]'",
    added_in_version="before_v9",
    text=r"""
Gives the angles between directions of primitive vectors of the unit cell (in
degrees), as an **alternative** to the input array [[rprim]]. Will be used to set
up [[rprim]], that, together with the array [[acell]], will be used to define
the primitive vectors.

  * *angdeg(1)* is the angle between the 2nd and 3rd vectors,
  * *angdeg(2)* is the angle between the 1st and 3rd vectors,
  * *angdeg(3)* is the angle between the 1st and 2nd vectors,

If the three angles are equal within 1.0d-12 (except if they are exactly 90
degrees), the three primitive vectors are chosen so that the trigonal symmetry
that exchange them is along the z cartesian axis:

    R1 = ( a,    0,           c)
    R2 = (-a/2,  sqrt(3)/2*a, c)
    R3 = (-a/2, -sqrt(3)/2*a, c)

where a  2 + c 2 = 1.0
If the angles are not all equal (or if they are all 90 degrees), one will have
the following generic form:

  * R1 = (1, 0, 0)
  * R2 = (a, b, 0)
  * R3 = (c, d, e)

where each of the vectors is normalized, and form the desired angles with the others.
""",
),

Variable(
    abivarname="asr",
    varset="eph",
    vartype="integer",
    topics=['DFPT_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Acoustic Sum Rule",
    added_in_version="before_v9",
    text=r"""
Govern the imposition of the Acoustic Sum Rule (ASR) in phonon calculations.
Same meaning as the corresponding anaddb variable.
""",
),

Variable(
    abivarname="atvshift",
    varset="ffield",
    vartype="real",
    topics=['DFT+U_expert'],
    dimensions=['[[natvshift]]', '[[nsppol]]', '[[natom]]'],
    defaultval=MultipleValue(number=None, value=0.0),
    mnemonics="ATomic potential (V) energy SHIFTs",
    characteristics=['[[DEVELOP]]'],
    requires="[[usepawu]] /= 0 and [[natvshift]] in [5,7]",
    added_in_version="before_v9",
    text=r"""
Defines for each atom and each spin channel (at present, can only be used with
[[nsppol]] = 1 or 2, like the +U scheme), a possible potential shift, for the d
(with [[lpawu]] = 2, [[natvshift]] = 5), or f states (with [[lpawu]] = 3,
[[natvshift]] = 7). In the case of d states, and 2 spin channels, a set of 10
numbers for each atom must be defined. The first set of 5 numbers corresponds
to real spherical harmonics m=-2 to m=+2 for the spin-up channel, the second
set of 5 numbers corresponds to real spherical harmonics m=-2 to m=+2 for the
spin-down channel. In the case of f states, the same ordering applies, for
sets of 7 numbers, corresponding to m=-3 to m=+3.

!!! important

    [[usepawu]] should be non-zero, [[lpawu]] should be 2 or 3.
""",
),

Variable(
    abivarname="autoparal",
    varset="paral",
    vartype="integer",
    topics=['parallelism_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="AUTOmatisation of the PARALlelism",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
This input variable is used only when running ABINIT in parallel and for Ground-State calculations.
It controls the automatic determination of parameters related to parallel work
distribution (if not imposed in input file). Given a total number of
processors, ABINIT can find a suitable distribution that fill (when possible)
all the different levels of parallelization. ABINIT can also determine optimal
parameters for the use of parallel Linear Algebra routines (using Scalapack or Cuda, at present).
The different values are:

  * 0 --> The *autoparal* feature is deactivated. For ground-state and response function calculations,
    ABINIT can only activate automatically the parallelism over spins and k-points.

  * 1 --> The number of processors per parallelization level is determined by mean of a simple
    (but relatively efficient) heuristic method. A scaling factor is attributed to each level and an simple
    speedup factor is computed. The selected parameters are those giving the best speedup factor.
    Possibly concerned parameters: [[npimage]], [[npkpt]], [[npspinor]], [[npfft]], [[npband]], [[bandpp]].

  * 2 --> The number of processors per parallelization level is first determined by mean of a simple
    (but relatively efficient) heuristic method (see 1 above). Then the code performs a series of
    small benchmarks using the scheme applied for the LOBPCG algorithm (see [[wfoptalg]] = 4 or 14).
    The parallel distribution is then changed according to the benchmarks.
    Possibly concerned parameters: [[npimage]], [[npkpt]], [[npspinor]], [[npfft]], [[npband]], [[bandpp]].

  * 3 --> Same as *autoparal* = 1, plus automatic determination of Linear Algebra routines parameters.
    In addition, the code performs a series of small benchmarks using the Linear
    Algebra routines (ScaLapack or Cuda-GPU). The parameters used to optimize
    Linear Algebra work distribution are then changed according to the benchmarks.
    Possibly concerned parameters (in addition to those modified for
    *autoparal* = 1): [[use_slk]], [[np_slk]], [[gpu_linalg_limit]]

  * 4 --> combination of *autoparal* 2 and 3.

Note that *autoparal* = 1 can be used on every set of processors;
*autoparal* > 1 should be used on a sufficiently large number of MPI process.
Also note that *autoparal* can be used simultaneously with [[max_ncpus]]; in
this case, ABINIT performs an optimization of process distribution for each
total number of processors from 2 to [[max_ncpus]]. A weight is associated to
each distribution and the higher this weight is the better the distribution
is. After having printed out the weights, the code stops.
""",
),

Variable(
    abivarname="auxc_ixc",
    varset="gstate",
    vartype="integer",
    topics=['Hybrids_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="AUxiliary XC functional for hybrid functional, IXC number",
    added_in_version="before_v9",
    text=r"""
Specification of an auxiliary exchange-correlation functional, thanks to its
[[ixc]] value, to possibly replace the heavy evaluation of an hybrid
functional at specific occasions, e.g. when the Fock operator is frozen during
the self-consistent cycle, thanks to [[fockoptmix]] == 11, or when evaluating
the correction to forces due to the density residual. This auxiliary exchange-
correlation functional might be rescaled, thanks to [[auxc_scal]] when
[[fockoptmix]] == 11. If [[gwcalctyp]] == 5, 15 or 25, *auxc_ixc* refers to
[[ixc_sigma]] instead of [[ixc]].
""",
),

Variable(
    abivarname="auxc_scal",
    varset="gstate",
    vartype="real",
    topics=['Hybrids_useful'],
    dimensions="scalar",
    defaultval=1.0,
    mnemonics="AUxiliary XC functional for hybrid functional- SCALing factor",
    added_in_version="before_v9",
    text=r"""
Possible scaling factor for the auxiliary exchange-correlation functional
defined by [[auxc_ixc]] that has the goal to replace the Fock operator or
hybrid functional when [[fockoptmix]] == 11.

The default value 1.0 corresponds to the unmodified xc functional. When the
auxiliary functional is used to replace the hybrid functional in SCF loops, a
value of 1.5 has been observed to speed up somehow the convergence.
""",
),

Variable(
    abivarname="awtr",
    varset="gw",
    vartype="integer",
    topics=['Susceptibility_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics=r"evaluate the Adler-Wiser expression of $\chi^{0}_{KS}$ assuming Time-Reversal",
    requires="[[optdriver]] == 3",
    added_in_version="before_v9",
    text=r"""
This input variable defines whether the irreducible polarizability
$\chi^{0}_{KS}$ is evaluated taking advantage of time-reversal symmetry or not.

* 0 --> Use the "standard" Adler-Wiser expression without assuming time-reversal symmetry.
  In this case, the irreducible polarizability is calculated summing over all possible
  electronic transitions (both resonant and anti-resonant).
* 1 --> Take advantage of time-reversal symmetry to halve the number of transitions to be explicitly considered.
  This method leads to a decrease in the CPU time by a factor two with respect to the *awtr* = 0 case.

!!! important

    Note that the parallel algorithm [[gwpara]] = 2 is not compatible with the choice *awtr* = 0.
""",
),

Variable(
    abivarname="bandpp",
    varset="paral",
    vartype="integer",
    topics=['parallelism_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="BAND Per Processor",
    characteristics=['[[DEVELOP]]'],
    requires="[[paral_kgb]] == 1",
    added_in_version="before_v9",
    text=r"""
Control the size of the block in the LOBPCG algorithm. This keyword works only
with [[paral_kgb]] = 1 and has to be either 1 or a multiple of 2.

With [[npband]] = 1:

* 1 --> band-per-band algorithm
* n --> The minimization is performed using [[nband]] blocks of n bands.

!!! warning

    [[nband]] has to be an integer.

With [[npband]] $\ne$ 1:

* 1 --> The minimization is performed using [[nband]] / [[npband]] blocks of [[npband]] bands.
* n --> The minimization is performed using [[nband]] / ([[npband]] $\times$ n) blocks of [[npband]]  $\times$ n bands.

!!! warning

    [[nband]] / ([[npband]]  $\times$ n) has to be an integer.

By minimizing a larger number of bands together in LOBPCG, we increase the
convergence of the residual. The better minimization procedure (as concerns
the convergence, but not as concerns the speed) is generally performed by
using *bandpp*  $\times$ [[npband]] = [[nband]]. Put *bandpp* = 2 when [[istwfk]] = 2
(the time spent in FFTs is divided by two).
""",
),

Variable(
    abivarname="bdberry",
    varset="ffield",
    vartype="integer",
    topics=['Berry_basic'],
    dimensions=[4],
    defaultval=MultipleValue(number=4, value=0),
    mnemonics="BanD limits for BERRY phase",
    requires="[[berryopt]] in [1, 2, 3] and [[nberry]] > 0",
    added_in_version="before_v9",
    text=r"""
Give the lower band and the upper band of the set of bands for which the Berry
phase must be computed. Irrelevant if [[nberry]] is not positive. When
[[nsppol]] is 1 (no spin-polarisation), only the two first numbers, giving the
lower and highest bands, are significant. Their occupation number is assumed
to be 2. When [[nsppol]] is 2 (spin-polarized calculation), the two first
numbers give the lowest and highest bands for spin up, and the third and
fourth numbers give the lowest and highest bands for spin down. Their
occupation number is assumed to be 1.

!!! important

    Presently, *bdberry* MUST be initialized by the user in case of a Berry
    phase calculation with [[berryopt]] = 1, 2, or 3: the above-mentioned default
    will cause an early exit.
""",
),

Variable(
    abivarname="bdeigrf",
    varset="dfpt",
    vartype="integer",
    topics=['TDepES_basic'],
    dimensions="scalar",
    defaultval=-1,
    mnemonics="BanD for second-order EIGenvalues from Response-Function",
    requires="[[ieig2rf]] in [1,2,3,4,5]",
    added_in_version="before_v9",
    text=r"""

The variable *bdeigrf* is the maximum number of bands for which the second-
order eigenvalues must be calculated: the full number of bands is still used
during the computation of these corrections.

If *bdeigrf* is set to -1, the code will automatically set *bdeigrf* equal to [[nband]].
""",
),

Variable(
    abivarname="bdgw",
    varset="gw",
    vartype="integer",
    topics=['GW_basic', 'SelfEnergy_basic'],
    dimensions=[2, '[[nkptgw]]', '[[nsppol]]'],
    defaultval=MultipleValue(number=None, value=0),
    mnemonics="BanDs for GW calculation",
    requires="[[optdriver]] in [4, 7]",
    added_in_version="before_v9",
    text=r"""
For each k-point with number `ikptgw` in the range (1:[[nkptgw]]) and each spin
index `isppol`, **bdgw**(1,`ikptgw`,`isppol`) is the number of the lowest band for
which the self-energy computation must be done.
**bdgw**(2,`ikptgw`,`isppol`) gives the index of the highest band for which the self-energy computation must be done.

!!! note

    The initial values given in the input file might be changed inside
    the code so that all the degenerate states at a given k-point and spin are
    included. This might happen when [[symsigma]] = 1 is used or in the case of
    self-consistent GW calculations.
    When [[symsigma]] == 1, indeed, the diagonal matrix elements of the self-energy are
    obtained by averaging the unsymmetrized results in the subspace spanned by the degenerate states.

When [[gwcalctyp]] >= 20, the quasi-particle wavefunctions are computed and
represented as linear combination of Kohn-Sham wavefunctions. In this case
**bdgw** designates the range of KS wavefunctions used as basis set. For each
k-point, indeed, the quasi-particle wavefunctions are expanded considering only
the KS states between **bdgw**(1,`ikptgw`,`isppol`) and **bdgw**(2,`ikptgw`,`isppol`).

For self-consistent calculations, on the other hand, the basis set used to
expand the GW wavefunctions should include all the degenerate states belonging
to the same irreducible representation. Only in this case, indeed, the initial
symmetries and energy degenerations are preserved.
""",
),

Variable(
    abivarname="berryopt",
    varset="ffield",
    vartype="integer",
    topics=['Berry_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="BERRY phase OPTions",
    added_in_version="before_v9",
    text=r"""
Specifies the use of Berry phase for the computation of either the
polarization, the derivatives with respect to the wavevector, or finite
electric field calculations.

  * 0 --> no computation of expressions relying on a Berry phase (default)
  * 1 --> the computation of Berry phases is activated (berryphase routine)
  * 2 --> the computation of derivatives with respect to the wavevector, thanks to the
    Berry phase finite-difference formula, is activated (uderiv routine)
  * 3 --> same as option 1 and 2 together

!!! note
    Note that options 1 to 3 require the use of a serial build of Abinit.**

  * -1 --> alternative computation of Berry phases (berryphase_new routine)
  * -2 --> alternative computation of derivatives with respect to the wavevector,
    thanks to the Berry phase finite-difference formula (berryphase_new routine)
  * -3 --> same as option -1 and -2 together

!!! note
    Options -1 to -3 permit use of a parallel build and will be preferred by
    most users.

  * 4 --> finite electric field calculation (unreduced E-field)
  * 6 --> finite electric displacement field calculation (unreduced D-field)
  * 14 --> finite reduced electric field calculation
  * 16 --> finite electric displacement field calculation
  * 17 --> mixed electric boundary condition: finite reduced electric field in some directions,
  finite reduced electric displacement field along other directions. See variable [[jfielddir]] for more details.

Other related input variables are:

  * in case of *berryopt* = 1,2, or 3: [[bdberry]] and [[kberry]]; also, [[nberry]] must be larger than 0;
  * in case of *berryopt* = -1,-2, or -3: the variable [[rfdir]] must be used to specify the primitive vector along which the projection of the polarization or the ddk will be computed. For example if *berryopt* = -1 and [[rfdir]] = 1 0 0, the projection of the polarization along the reciprocal lattice vector $G_1$ is computed. In case [[rfdir]] = 1 1 1, ABINIT computes the projection of P along $G_1$, $G_2$ and $G_3$ and transforms the results to cartesian coordinates;
  * in cases where *berryopt* is negative, [[berrystep]] allow a computation of multiple-step Berry phase in order to accelerate the convergence.
  * [[efield]] and [[rfdir]] in case of *berryopt* = 4;

The cases *berryopt* = -1,-2,-3, 4, 6, 7, 14, 16, and 17 have to be used with [[occopt]] = 1.

The cases *berryopt* = -1 and 4, 6, 7, 14, 16, 17 are compatible with PAW,
howevever, if in these cases one uses [[kptopt]] /= 3, one must also use only
symmorphic symmetries (either because the space group is symmorphic or the
variable [[symmorphi]] is set to zero).

For a phonon calculation under a finite electric field, respect the following procedure.

  * a) Run a scf ground-state calculation at zero electric field to get wavefunctions to initialize the ground-state calculation in finite electric fields.
  * b) Run a scf ground-state calculation in finite electric field. The electric field is controlled by the input variable [[efield]]. *berryopt* should be 4. The input variable [[kptopt]] should be set to be 2.
  * c) Based on the wave functions obtained in step (2), perform phonon calculation by setting *berryopt* = 4, [[kptopt]] = 3 and The same value of [[efield]] than in step 2. [[nsym]] should be set to 1 currently but this restriction may be removed later. The other parameters are the same as phonon calculation at zero electric field.

!!! important

    The choice of k-point sampling N x N x N should be the same in the three runs and N should be an even number.

In case of finite electric and displacement field calculations
(*berryopt* = 4,6,7,14,16,17), see also the input variables [[berrysav]],
[[dfield]], [[red_dfield]], [[red_efield]], [[ddamp]]
""",
),

Variable(
    abivarname="berrysav",
    varset="ffield",
    vartype="integer",
    topics=['Berry_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="BERRY SAVe",
    added_in_version="before_v9",
    text=r"""
  * 0 --> for finite electric field calculation ([[berryopt]] = 4/14),
    the polarization branch will be chosen on each iteration from (-pi, pi).
    For finite electric displacement field calculation([[berryopt]] = 6/7/16/17),
    the polarization will be chosen to minimize the internal energy.
  * 1 --> the polarization will be kept in the same branch on each iteration.
    At the end of the run, a file "POLSAVE" will be saved containing the reduced polarization in atomic units.

    !!! note

        Make sure that "POLSAVE" is empty or it does not exist before the calculation, or else that
        it specifies the desired polarization branch.
""",
),

Variable(
    abivarname="berrystep",
    varset="ffield",
    vartype="integer",
    topics=['Berry_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="BERRY phase: multiple STEP",
    requires="0 > [[berryopt]]",
    added_in_version="before_v9",
    text=r"""
If [[berryopt]] is negative, this variable is used to compute berry phases
using multiple discrete steps, in order to accelerate convergence. The single-
step berry phase is the standard calculation using strings of k-points based
on overlap of Bloch function separated by $dk$, while the two-step berry phase
use strings use overlaps based on dk and $2*dk$, the three-step use overlaps
based on dk, $2*dk$ and $3*dk$...

The default value of this variable is 1, meaning that only the single-step
berry phase calculation is done. If a larger value is set, ABINIT will compute
all the multiple-step berry phase from the single-step to the
*berrystep*-step, and use the large-step values of berry phase to correct
the single-step berry phase. Use with care: while experience is still to be
gained with this procedure, the outlook is promising.
""",
),

Variable(
    abivarname="bfield",
    varset="ffield",
    vartype="real",
    topics=['MagField_expert'],
    dimensions=[3],
    defaultval=MultipleValue(number=3, value=0.0),
    mnemonics="finite B FIELD calculation",
    added_in_version="before_v9",
    text=r"""
Perform finite magnetic field calculation.

!!! important

    THIS CODE IS UNDER DEVELOPMENT AND IS NOT READY FOR USE.
""",
),

Variable(
    abivarname="bmass",
    varset="rlx",
    vartype="real",
    topics=['MolecularDynamics_expert'],
    dimensions="scalar",
    defaultval=10,
    mnemonics="Barostat MASS",
    added_in_version="before_v9",
    text=r"""
bmass is the mass of the barostat when [[ionmov]] = 13 (constant pressure molecular dynamics)
""",
),

Variable(
    abivarname="boxcenter",
    varset="gstate",
    vartype="real",
    topics=['TDDFT_basic'],
    dimensions=[3],
    defaultval=[0.5, 0.5, 0.5],
    mnemonics="BOX CENTER",
    added_in_version="before_v9",
    text=r"""
Defines the center of the box, in reduced coordinates. At present, this
information is only used in the case of Time-Dependent DFT computation of the
oscillator strength. One must take *boxcenter* such as to be roughly the center
of the cluster or molecule. The default is sensible when the vacuum
surrounding the cluster or molecule has xred 0 or 1. On the contrary, when the
cluster or molecule is close to the origin, it is better to take *boxcenter* = [0.0, 0.0, 0.0].
""",
),

Variable(
    abivarname="boxcutmin",
    varset="gstate",
    vartype="real",
    topics=['Planewaves_useful', 'TuningSpeed_basic'],
    dimensions="scalar",
    defaultval=2.0,
    mnemonics="BOX CUT-off MINimum",
    added_in_version="before_v9",
    text=r"""
The box cutoff ratio is the ratio between the wavefunction plane wave sphere
radius, and the radius of the sphere that can be inserted in the FFT box, in reciprocal space.

In order for the density to be exact (in the case of the plane wave part, not the PAW on-site terms),
this ratio should be at least two. If one uses a smaller ratio (e.g 1.5), one will gain speed, at the expense of accuracy.
In the case of pure ground state calculation (e.g. for the determination of geometries), this is sensible.
It should also be noticed that using a value of [[boxcutmin]] too close to one can lead to runtime errors in the FFT routines.
A value larger than to 1.1 is therefore recommended.

Prior to v8.9, the use of boxcutmin for DFPT calculations was forbidden. However, after testing, it was seen that
the deterioration in phonon band structures could be alleviated to a large extent by the imposition
of the Acoustic Sum Rule [[asr]].
""",
),

Variable(
    abivarname="brvltt",
    varset="geo",
    vartype="integer",
    topics=['UnitCell_useful', 'SmartSymm_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="BRaVais LaTTice type",
    requires="[[spgroup]] != 0",
    added_in_version="before_v9",
    text=r"""
Set the type of Bravais lattice. The cell defined by [[acell]] and [[rprim]]
or [[angdeg]] should be the CONVENTIONAL cell.

If brvltt=0, the code will assign brvltt from the space group information
[[spgroup]], and produce the symmetry operations for the conventional unit
cell. If the conventional cell is not primitive, the user should set [[chkprim]] = 0.

If brvltt=-1, the code will assign brvltt from the space group information,
then reduce the unit cell to a primitive unit cell. The echo of [[acell]] and
[[rprim]] might thus differ from those derived directly from the input
variables. Also, the input variable [[xred]] will refer to the CONVENTIONAL
unit cell, but its echo will refer to the preprocessed PRIMITIVE unit cell.
There is of course no problem with [[xangst]] and [[xcart]], as they are
independent of the unit cell.

The echo of *brvltt* in the output file will be one of the following Bravais lattices:

* 1 = Primitive with no associated translations
* 2 = Inner centered with (a/2 + b/2 + c/2) associated translation
* 3 = Face centered with (a/2 + b/2; b/2 + c/2; c/2 + a/2) associated translations
* 4 = C - centered with (a/2 + b/2) associated translation
* 5 = A - centered with (b/2 + c/2) associated translation
* 6 = B - centered with (c/2 + a/2) associated translation
* 7 = Rhombohedral lattice.

The user might also input directly these values, although they might not be
consistent with [[spgroup]].

The space groups 146, 148, 155, 160, 161, 166, 167, when used with
[[spgaxor]] = 1 (hexagonal axes) will have *brvltt* = 7 and two associated
translations: (2/3, 1/3, 1/3) and (1/3, 2/3, 2/3).
For more details see the [[help:spacegroup]].
""",
),

Variable(
    abivarname="bs_algorithm",
    varset="bse",
    vartype="integer",
    topics=['BSE_basic'],
    dimensions="scalar",
    defaultval=2,
    mnemonics="Bethe-Salpeter ALGORITHM",
    requires="[[optdriver]] == 99",
    added_in_version="before_v9",
    text=r"""
This input variable defines the algorithm employed to calculate the macroscopic dielectric function.
Possible values are in [1, 2, 3]:

* 1 --> The macroscopic dielectric is obtained by performing a direct diagonalization
  of the excitonic Hamiltonian. Advantages: It gives direct access to the excitonic eigenvalues
  as well as to the oscillator strengths. Drawbacks: It is a very CPU- and memory-consuming approach
  as the size of the Hamiltonian scales as $(n_k * n_c * n_v)^2$ where $n_k$ is the number of k-point
  in the **full** Brillouin zone, and $n_c$ and $n_v$ are the number of conduction and valence states, respectively.
  Pros: It can be used both for resonant-only and resonant + coupling calculations (non Tamm-Dancoff approximation).

* 2 --> Haydock iterative method. The macroscopic dielectric function is obtained by iterative applications
  of the Hamiltonian on a set of vectors in the electron-hole space.
  Advantages: It is less memory demanding and usually faster than the direct diagonalization provided
  that [[zcut]] is larger than the typical energy spacing of the eigenvalues. Drawbacks:
  It is an iterative method therefore the convergence with respect to [[bs_haydock_niter]] should be checked.
  It is not possible to have direct information on the exciton spectrum, oscillator strengths and excitonic wave functions.
  For the time being *bs_algorithm* = 2 cannot be used for calculations in which the coupling
  term is included (Tamm-Dancoff approximation).

* 3 --> Conjugate-gradient method. This method allows one to find the few first excitonic eigenvalues.
  Only available for resonant calculations (Tamm-Dancoff approximation).
""",
),

Variable(
    abivarname="bs_calctype",
    varset="bse",
    vartype="integer",
    topics=['BSE_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Bethe-Salpeter CALCulation TYPE",
    requires="[[optdriver]] == 99",
    added_in_version="before_v9",
    text=r"""
Possible values are in [1, 2, 3].

* 1 --> use the KS eigenvalues and wave functions stored in the WFK file to construct the transition space
* 2 --> The transition space is constructed with Kohn-Sham orbitals but the energies are read from an external GW file
* 3 --> QP amplitudes and energies will be read from the QPS file and used to construct H_ex.
  Not coded yet because <\psi|r|\psj>^QP should be calculated taking into account the non-locality
  of the self-energy in the commutator [H,r].
""",
),

Variable(
    abivarname="bs_coulomb_term",
    varset="bse",
    vartype="integer",
    topics=['BSE_expert'],
    dimensions="scalar",
    defaultval=11,
    mnemonics="Bethe-Salpeter COULOMB TERM",
    requires="[[optdriver]] == 99",
    added_in_version="before_v9",
    text=r"""
This variable governs the choice among the different options that are
available for the treatment of Coulomb term of the Bethe-Salpeter Hamiltonian.
**bs_coulomb_term** is the concatenation of two digits, labelled (A) and (B).

The first digit (A) can assume the values 0, 1, 2:

  * 0 --> The Coulomb term is not computed. This choice is equivalent to computing the RPA spectrum
    but using the representation in transition space instead of the more efficient approach based on the sum over states.

  * 1 --> The Coulomb term is computed using the screened interaction read
    from an external SCR file (standard excitonic calculation).

  * 2 --> The Coulomb term is computed using a model screening function
    (useful for convergence studies or for reproducing published results).

The second digit (B) can assume the values 0,1:

  * 0 --> Use a diagonal approximation for $W_{\GG\GG'}$ (mainly used for accelerating convergence studies).

  * 1 --> The Coulomb term is correctly evaluated using the truly non-local screening $W(\rr,\rr')$.
""",
),

Variable(
    abivarname="bs_coupling",
    varset="bse",
    vartype="integer",
    topics=['BSE_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Bethe-Salpeter COUPLING",
    requires="[[optdriver]] == 99",
    added_in_version="before_v9",
    text=r"""
The **bs_coupling** input variable defines the treatment of the coupling block
of the Bethe-Salpeter Hamiltonian. Possible values are 0, 1.

  * 0 --> The coupling block is neglected (the so-called Tamm-Dancoff approximation).
    The code runs faster and the Hamiltonian matrix requires less memory (factor 4).
    It is a good approximation for the absorption spectrum which only requires the knowledge of $\Im(\epsilon)$.
    The reliability of this approximation should be tested in the case of EELF calculations.

  * 1 --> The coupling term is included (non Tamm-Dancoff approximation).
""",
),

Variable(
    abivarname="bs_eh_cutoff",
    varset="bse",
    vartype="integer",
    topics=['BSE_expert'],
    dimensions=[2],
    defaultval=["-inf", "inf"],
    mnemonics="Bethe-Salpeter Electron-Hole CUTOFF",
    requires="[[optdriver]] == 99",
    added_in_version="before_v9",
    text=r"""
Used to define a cutoff in the e-h basis set. Only those transitions
whose energy is between bs_eh_cutoff(1) and bs_eh_cutoff(2) will be considered
in the construction of the e-h Hamiltonian.
""",
),

Variable(
    abivarname="bs_exchange_term",
    varset="bse",
    vartype="integer",
    topics=['BSE_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Bethe-Salpeter EXCHANGE TERM",
    requires="[[optdriver]] == 99",
    added_in_version="before_v9",
    text=r"""
* 0 --> The exchange term is not calculated. This is equivalent to neglecting
  local field effects in the macroscopic dielectric function.
* 1 --> The exchange term is calculated and added to the excitonic Hamiltonian.
""",
),

Variable(
    abivarname="bs_freq_mesh",
    varset="bse",
    vartype="real",
    topics=['BSE_basic'],
    dimensions=[3],
    defaultval=[0.0, 0.0, 0.01],
    mnemonics="Bethe-Salpeter FREQuency MESH",
    characteristics=['[[ENERGY]]'],
    requires="[[optdriver]] == 99",
    added_in_version="before_v9",
    text=r"""
**bs_freq_mesh(1)** defines the first frequency for the calculation of the macroscopic dielectric function.

**bs_freq_mesh(2)** gives the last frequency for the calculation of the
macroscopic dielectric function. If zero, **bs_freq_mesh(2)** is set automatically to MAX(resonant_energy) + 10%.

**bs_freq_mesh(3)** gives the step of the linear mesh used for evaluating the macroscopic dielectric function.
""",
),

Variable(
    abivarname="bs_hayd_term",
    varset="bse",
    vartype="integer",
    topics=['BSE_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Bethe-Salpeter HAYdock TERMinator",
    requires="[[optdriver]] == 99 and [[bs_algorithm]] == 2",
    added_in_version="before_v9",
    text=r"""
Defines how to terminate the continued fraction expression for the dielectric function.
The terminator reduces the number of iterations needed to converge
by smoothing the oscillation in the high energy part of the spectrum

  * 0 --> No terminator. The contribution given by the terms missing in the Lanczos chain are set to zero.
  * 1 --> Use the terminator function. The particular expression depends on the type of calculation:
    In the resonant-only case, the $a_i$ and $b_i$ coefficients for $i > \text{niter}$, are replaced
    by their values at $i = \text{niter}$.
    If the coupling block is included, the terminator function is the one described in [[cite:Rocca2008]].
""",
),

Variable(
    abivarname="bs_haydock_niter",
    varset="bse",
    vartype="integer",
    topics=['BSE_expert'],
    dimensions="scalar",
    defaultval=100,
    mnemonics="Bethe-Salpeter HAYDOCK Number of ITERations",
    requires="[[optdriver]] == 99 and [[bs_algorithm]] == 2",
    added_in_version="before_v9",
    text=r"""
**bs_haydock_niter** defines the maximum number of iterations used to
calculate the macroscopic dielectric function.
The iterative algorithm stops when the difference between two consecutive
evaluations of the optical spectra is less than [[bs_haydock_tol]].
""",
),

Variable(
    abivarname="bs_haydock_tol",
    varset="bse",
    vartype="real",
    topics=['BSE_expert'],
    dimensions=[2],
    defaultval=[0.02, 0.0],
    mnemonics="Bethe-Salpeter HAYDOCK TOLerance",
    requires="[[optdriver]] == 99 and [[bs_algorithm]] == 2",
    added_in_version="before_v9",
    text=r"""
Defines the convergence criterion for the Haydock iterative method.
The iterative algorithm stops when the difference between two consecutive
evaluations of the macroscopic dielectric function is less than **bs_haydock_tol(1)**.
The sign of **bs_haydock_tol(1)** defines how to estimate the convergence error.

A negative value signals that the converge should be reached for each frequency (strict criterion), while a positive
value indicates that the converge error is estimated by averaging over the entire frequency range (mild criterion).

**bs_haydock_tol(2)** defines the quantity that will be checked for convergence:

  * 0.0 --> both the real and the imaginary part must converge
  * 1.0 --> only the real part
  * 2.0 --> only the imaginary part

(The latter are real numbers, tolerance is 1.0d-6).
""",
),

Variable(
    abivarname="bs_interp_kmult",
    varset="bse",
    vartype="integer",
    topics=['BSE_useful'],
    dimensions=[3],
    defaultval=[0, 0, 0],
    mnemonics="Bethe-Salpeter INTERPolation K-point MULTiplication factors",
    requires="[[bs_interp_mode]] > 0 and [[bs_algorithm]] == 2 and [[bs_coupling]] == 0",
    added_in_version="before_v9",
    text=r"""
**bs_interp_kmult** defines the number of divisions used to generate the dense mesh in the interpolation.
[[ngkpt]] of the dense mesh = **bs_interp_kmult(:)** * [[ngkpt]] of the coarse mesh.
""",
),

Variable(
    abivarname="bs_interp_m3_width",
    varset="bse",
    vartype="real",
    topics=['BSE_useful'],
    dimensions="scalar",
    defaultval=1.0,
    mnemonics="Bethe-Salpeter INTERPolation Method3 WIDTH",
    requires="[[bs_interp_mode]] == 3 and [[bs_algorithm]] == 2 and [[bs_coupling]] == 0",
    added_in_version="before_v9",
    text=r"""
Defines the width of the region where divergence treatment is applied for BSE interpolation
""",
),

Variable(
    abivarname="bs_interp_method",
    varset="bse",
    vartype="integer",
    topics=['BSE_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Bethe-Salpeter INTERPolation METHOD",
    requires="[[bs_interp_mode]] > 0 and [[bs_algorithm]] == 2 and [[bs_coupling]] == 0",
    added_in_version="before_v9",
    text=r"""
*bs_interp_method* selects the interpolation method::

  * 0 --> Interpolate using Y. Gillet technique with 8 neighbours (see [[cite:Gillet2016]]).
  * 1 --> Interpolation using Rohlfing & Louie technique (see above-mentioned article and [[cite:Rohlfing2000]])
""",
),

Variable(
    abivarname="bs_interp_mode",
    varset="bse",
    vartype="integer",
    topics=['BSE_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Bethe-Salpeter INTERPolation MODE",
    requires="[[bs_interp_mode]] > 0 and [[bs_algorithm]] == 2 and [[bs_coupling]] == 0",
    added_in_version="before_v9",
    text=r"""
*bs_interp_mode* selects the mode of interpolation:

  * 0 --> No interpolation. Standard Bethe-Salpeter computation is performed
  * 1 --> Simple interpolation
  * 2 --> Treatment of the divergence on the whole set of dense k-points
  * 3 --> Treatment of the divergence along the diagonal in k-space and simple interpolation elsewhere.
""",
),

Variable(
    abivarname="bs_interp_prep",
    varset="bse",
    vartype="integer",
    topics=['BSE_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Bethe-Salpeter INTERPolation PREParation",
    requires="[[bs_interp_mode]] > 0 and [[bs_algorithm]] == 2 and [[bs_coupling]] == 0",
    added_in_version="before_v9",
    text=r"""
*bs_interp_prep* allows one to trigger the preparation of the interpolation with method 2 or method 3.
It generates the decomposition of BSR in a, b, c coefficients used for the interpolation.
""",
),

Variable(
    abivarname="bs_interp_rl_nb",
    varset="bse",
    vartype="integer",
    topics=['BSE_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Bethe-Salpeter INTERPolation Rohlfing & Louie NeighBour",
    requires="[[bs_interp_mode]] > 0 and [[bs_algorithm]] == 2 and [[bs_interp_method]] == 1 and [[bs_coupling]] == 0",
    added_in_version="before_v9",
    text=r"""
Gives the index of the neighbour that is used in the Rohlfing and Louie method ([[cite:Rohlfing2000]])
""",
),

Variable(
    abivarname="bs_loband",
    varset="bse",
    vartype="integer",
    topics=['BSE_compulsory'],
    dimensions=['[[nsppol]]'],
    defaultval=0,
    mnemonics="Bethe-Salpeter Lowest Occupied BAND",
    requires="[[optdriver]] == 99",
    added_in_version="before_v9",
    text=r"""
This variable defines the index of the lowest occupied band used for the
construction of the electron-hole basis set. For spin polarized calculations,
one must provide two separated indices for spin up and spin down.
An additional cutoff energy can be applied by means of the [[bs_eh_cutoff]] input variable.
""",
),

Variable(
    abivarname="bs_nstates",
    varset="bse",
    vartype="integer",
    topics=['BSE_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Bethe-Salpeter Number of STATES",
    requires="[[optdriver]] == 99 and [[bs_algorithm]] in [2, 3]",
    added_in_version="before_v9",
    text=r"""
**bs_nstates** defines the maximum number of excitonic states calculated in
the direct diagonalization of the excitonic matrix or in the conjugate-gradient method.
The number of states should be sufficiently large for a
correct description of the optical properties in the frequency range of interest.
""",
),

Variable(
    abivarname="builtintest",
    varset="dev",
    vartype="integer",
    topics=['Control_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="BUIT-IN TEST number",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
When [[builtintest]] is non-zero, the input file is a special one, that runs
very quickly, and that is accompanied by a specific analysis by ABINIT, at the
end of the run, against a hard-coded value of total energy (and possibly
stresses, forces ...). The echo of the analysis is done in the STATUS file. In
particular, such built-in tests can be used to check quickly whether ABINIT
fallbacks have been connected or not (bigdft, netcdf, libxc, wannier90). At
present, [[builtintest]] = 1... 7 are allowed. See more information in tests/built-in/README.
""",
),

Variable(
    abivarname="bxctmindg",
    varset="paw",
    vartype="real",
    topics=['PAW_expert'],
    dimensions="scalar",
    defaultval=2.0,
    mnemonics="BoX CuT-off MINimum for the Double Grid (PAW)",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
The box cut-off ratio is the ratio between the wavefunction plane wave sphere
radius, and the radius of the sphere that can be inserted in the FFT box, in
reciprocal space.

If the density was generated only from wavefunctions, this ratio should be at
least two in order for the density to be exact. If one uses a smaller ratio,
one will gain speed, at the expense of accuracy. In case of pure ground state
calculation (e.g. for the determination of geometries), this is sensible.
However, the wavefunctions that are obtained CANNOT be used for starting
response function calculation.

However, some augmentation charge is always added in PAW, and even with the
box cut-off ratio larger than two, the density is never exact. Sometimes, this
ratio must be much larger than two for the computation to be converged at the
required level of accuracy.
""",
),

Variable(
    abivarname="cd_customnimfrqs",
    varset="gw",
    vartype="integer",
    topics=['FrequencyMeshMBPT_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Contour Deformation CUSTOM IMaginary FReQuencieS",
    requires="([[optdriver]] ==3 or [[optdriver]] ==4) and [[gwcalctyp]] in [2,9,12,19,22,29]",
    added_in_version="before_v9",
    text=r"""
[[cd_customnimfrqs]] lets the user define the grid points along the imaginary
axis by hand. Set this to the number of frequencies you want. The frequencies
are specified with [[cd_imfrqs]].
""",
),

Variable(
    abivarname="cd_frqim_method",
    varset="gw",
    vartype="integer",
    topics=['FrequencyMeshMBPT_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Contour Deformation FReQuency integration on IMaginary axis Method",
    requires="[[optdriver]] ==4 and [[gwcalctyp]] in [2, 9, 12, 19, 22, 29]",
    added_in_version="before_v9",
    text=r"""
[[cd_frqim_method]] defines the choice of integration method along the
imaginary frequency axis for Contour Deformation calculations. The default
method is very robust, fast and optimal for the vast majority of cases.
However, for very accurate ("paranoid level") convergence studies, ABINIT
offers the possibility of a variety of methods and grids. Note that as one
starts to change the defaults, one needs to carefully consider the grid used.
Therefore we recommend that in addition to reading the information below, the
user reads the description of the input variables [[freqim_alpha]],
[[nfreqim]], [[ppmfrq]], [[gw_frqim_inzgrid]].

The integration to be performed for each matrix element of the self energy
along the imaginary axis is of the form:

$$ \langle i|\Sigma(\omega)|j \rangle  \propto
   \sum_s C_s \int_0^\infty d(i\omega^\prime)
   \frac{(\omega - \epsilon_s)}{(\omega-\epsilon_s)^2 + i{\omega^\prime}^2}
   f(i\omega^\prime),
$$

where  $\omega$ is the frequency point along the real axis, $\epsilon_s$ is an
eigenvalue, and  $i\omega^\prime$ is the variable along the imaginary axis. Thus the
function to be integrated is a Lorentzian weight function centred on the
origin (whose FWHM is decided by $|\omega-\epsilon_s|$), times a function. The
function is related to the inverse dielectric matrix. It might have a peaked
structure near the origin and is very smooth otherwise. the function decays
asymptotically as  $1/i\omega^\prime$, so the whole integral converges as this to
the third power.

  * **cd_frqim_method = 1 - Histogram:** This is the **default** method where the function  $f(i\omega^\prime)$ is approximated by a histogram, and the Lorentzian is integrated analytically in each sub-interval. See the section on grids below for a description of the default grid. This method combined with the default grid is the fastest and optimised for the use of few points along the imaginary axis.
  * **cd_frqim_method = 2 - Trapezoid:** The next step up from the histogram approximation in the previous method. The integration region is transformed  $[0, \infty] \rightarrow [0,1]$ with a proper weight depending on the width of the Lorentzian. In this space  $f(i\omega^\prime)$ is approximated by a linear function between grid points (trapezoids), and the integrand is integrated analytically in each sub-interval. This method tends to slightly overestimate contributions while the default method tends to slightly underestimate them, so the results from methods 1 and 2 should bracket the converged values. The asymptotic behaviour is explicitly taken into account by a fit using the last two grid points.
  * **cd_frqim_method = 3, 4, 5 - Natural Spline:** The function is transformed  $[0, \infty] \rightarrow [0,1]$. In this space  $f(i\omega^\prime)$ is approximated by a natural spline function whose starting and ending sections are linear. This transform is chosen so that the function should approach a linear function asymptotically as the integration interval approaches 1, so that the asymptotic behaviour is automatically taken into account. For each Lorentzian width (determined by $|\omega-\epsilon_s|$) the integrand is appropriately scaled in the interval $[0,1]$, and a nested Gauss-Kronrod (GK) numerical integration rule is performed. The integrand is evaluated at the GK nodes by means of a spline-fit. The order of the GK rule is controlled by the index of the method:
    * 3 --> Gauss 7 point, Kronrod 15 point rule.
    * 4 --> Gauss 11 point, Kronrod 23 point rule.
    * 5 --> Gauss 15 point, Kronrod 31 point rule.
There is rarely any difference to machine precision between these rules, and
the code will issue a warning if a higher-order rule is recommended.

**Grids for the integral along the imaginary axis:**

All the methods above should execute no matter what grid is used along the
imaginary axis, so this is very much under the control of the user. The only
requirement is that the grid be strictly increasing. The point at zero
frequency is assumed to lie on the real axis, so the calculation of that point
is controlled by [[nfreqre]] and corresponding variables. We highly recommend
extracting various elements of the dielectric matrix from the _SCR file using
the **Mrgscr** utility and plotting them for visual inspection.

  * **Default** - The default grid is an exponentially increasing grid given by the formula:

$$ i\omega^\prime_k = \frac{\omega_p}{\alpha-2}\left[ e^{\frac{2k}{N+1} \ln(\alpha-1)} -1 \right]. $$

Here  $\omega_p$ is the plasma frequency (by default determined by the average
density of the system, but this can be overridden by setting [[ppmfrq]]).  $N$
is the total number of grid points (set by [[nfreqim]]). $\alpha$ is a parameter
which determines how far out the final grid point will lie. The final point
will be at  $\alpha*\omega_p$ (the default is  $\alpha = 5$, and was hard-coded in
older versions of ABINIT). This grid is designed so that approximately half
the grid points are always distributed to values lower than the plasma
frequency, in order to resolve any peaked structure. If one seeks to increase
the outermost reach by increasing [[ppmfrq]] one must simultaneously take care
to increase [[nfreqim]] in order to have the appropriate resolution for the
low-frequency region. In more recent versions of ABINIT one can also simply
adjust the parameter  $\alpha$ by using [[freqim_alpha]]. This grid is optimised
for speed and accurate results with few grid points for **cd_frqim_method = 1**.

  * **Inverse z transform** - This grid is activated by the use of the variable [[gw_frqim_inzgrid]].
    This is the standard  $[0, \infty] \rightarrow [0,1]$ transform using the formula:

$$ i\omega^\prime = \omega_p \frac{z}{1-z}. $$

Here $\omega_p$ is the plasma frequency (default can be overridden by setting
[[ppmfrq]]). The grid points are then picked by an equidistant grid (number of
points set by [[nfreqim]]) in the interval  $z \subset [0,1]$. This grid can
easily be uniquely converged by just increasing [[nfreqim]]. Again the points
are distributed so that approximately half of them lie below the plasma frequency.

  * **User defined** - The user can also define their own grid using the variables [[cd_customnimfrqs]] and [[cd_imfrqs]].
    _With great power comes great responsibility!_

The **mrgscr** utility is handy in optimising the numerical effort expended in
convergence studies. By estimating the densest grid one can afford to
calculate in the SCR file, and successively removing frequencies from a single
file (using the utility), one only needs to perform the screening calculation
**once** on the dense mesh for a given convergence study. One can also use the
utility to merge independent screening calculations over q-points and frequency sections.
""",
),

Variable(
    abivarname="cd_full_grid",
    varset="gw",
    vartype="integer",
    topics=['FrequencyMeshMBPT_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Contour Deformation FULL GRID in complex plane",
    requires="[[optdriver]] == 3 and [[gwcalctyp]] in [2, 9, 12, 19, 22, 29]",
    added_in_version="before_v9",
    text=r"""
[[cd_full_grid]] enables the calculation of the screening [both chi0 and
epsilon^(-1)] on a grid in the first quadrant of the complex plane. The grid
is determined by the (tensor) product of the grid in real frequency and the
grid in imaginary frequency. In the SUS and SCR files the grid points are stored as follows:

      **Index:**
      1   . . .   nfreqre   nfrqre+1 . . . nfreqre+nfreqim   nfreqre+nfreqim+1 . . . nfreqre*nfreqim
      **Entry:**
      | purely real freq.  |     purely imaginary freq.     |      gridpoints in complex plane        |


The grid in the complex plane is stored looping over the real dimension as the
inner loop and the imaginary as the outer loop. The contents of the generated
SUS and SCR files can be extracted for visualisation and further analysis with
the **Mrgscr** utility.
""",
),

Variable(
    abivarname="cd_halfway_freq",
    varset="gw",
    vartype="real",
    topics=['FrequencyMeshMBPT_expert'],
    dimensions="scalar",
    defaultval=ValueWithUnit(units='eV', value=100.0),
    mnemonics="Contour Deformation tangent grid HALFWAY FREQuency",
    characteristics=['[[ENERGY]]'],
    requires="([[optdriver]] == 3 or [[optdriver]] == 4) and [[gwcalctyp]] in [2,9,12,19,22,29]",
    added_in_version="before_v9",
    text=r"""
[[cd_halfway_freq]] determines the frequency where half of the number of
points defined in [[nfreqre]] are used up. The tangent transformed grid is
approximately linear up to this point. To be used in conjunction with [[gw_frqre_tangrid]].
""",
),

Variable(
    abivarname="cd_imfrqs",
    varset="gw",
    vartype="real",
    topics=['FrequencyMeshMBPT_expert'],
    dimensions=['[[cd_customnimfrqs]]'],
    mnemonics="Contour Deformation IMaginary FReQuencieS",
    requires="[[optdriver]] == 3 and [[gwcalctyp]] in [2, 9, 12, 19, 22, 29] and [[cd_customnimfrqs]] != 0",
    added_in_version="before_v9",
    text=r"""
[[cd_imfrqs]] specifies the grid points for the imaginary axis. The number of
frequencies is set by the value of [[cd_customnimfrqs]]. For example,

    cd_customnimfrqs   5
    nfreqim            5
    cd_imfrqs          0.1  0.2  0.5  1.0  5.0

If [[nfreqim]] is not equal to [[cd_customnimfrqs]] a warning will be issued.

**Use at own risk!** The use of a custom grid makes it your responsibility
that the SUS and SCR files are valid in self-energy (i.e. [[optdriver]] = 4)
calculations, so caution is advised. Note that frequencies have to be strictly
increasing, and the point at zero frequency is **not** considered to be part
of the imaginary grid, but rather the grid along the real axis. The
calculation of that point should be controlled by [[nfreqre]] and related variables.
""",
),

Variable(
    abivarname="cd_max_freq",
    varset="gw",
    vartype="real",
    topics=['FrequencyMeshMBPT_expert'],
    dimensions="scalar",
    defaultval=ValueWithUnit(units='eV', value=1000.0),
    mnemonics="Contour Deformation grid MAXimum FREQuency",
    characteristics=['[[ENERGY]]'],
    requires="([[optdriver]] == 3 or [[optdriver]] == 4) and [[gwcalctyp]] in [2,9,12,19,22,29]",
    added_in_version="before_v9",
    text=r"""
[[cd_max_freq]] determines the frequency where all the points defined in
[[nfreqre]] are used up. To be used in conjunction with [[gw_frqre_tangrid]].
""",
),

Variable(
    abivarname="cd_subset_freq",
    varset="gw",
    vartype="integer",
    topics=['FrequencyMeshMBPT_expert'],
    dimensions=[2],
    defaultval=[1, '[[nfreqre]]'],
    mnemonics="Contour Deformation grid calculate SUBSET of FREQuencies",
    requires="[[optdriver]] == 3 and [[gwcalctyp]] in [2, 9, 12, 19, 22, 29] and  [[gw_frqre_tangrid]] == 0",
    added_in_version="before_v9",
    text=r"""
[[cd_subset_freq]] Specifies that only a subset of the frequencies defined by
[[nfreqre]] are to be calculated. The first index is the start and the second
the end, with index number 1 always being the origin. For example a
calculation with **[[nfreqre]] = 100** could be separated into two datasets with:

    subset_freq1   1   50
    subset_freq2   51  100

Any resulting susceptibility (_SUS) and screening (_SCR) files can then be
merged with the **mrgscr** utility.
""",
),

Variable(
    abivarname="charge",
    varset="gstate",
    vartype="real",
    topics=['Coulomb_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="CHARGE",
    added_in_version="before_v9",
    text=r"""
Used to establish charge balance between the number of electrons filling the
bands and the nominal [[charge]] associated with the atomic cores.
The code adds up the number of valence electrons provided by the
pseudopotentials of each type (call this "zval"), then add [[charge]], to get
the number of electrons per unit cell, [[nelect]].
Then, if [[iscf]] is positive, the code adds up the band occupancies (given in
array [[occ]]) for all bands at each k point, then multiplies by the k point
weight [[wtk]] at each k point. Call this sum "nelect_occ" (for the number of
electrons from occupation numbers). It is then required that: nelect_occ = [[nelect]].
To treat a neutral system, which is desired in nearly all cases, one must use
[[charge]] = 0. To treat a system missing one electron per unit cell, set [[charge]] = +1.
""",
),

Variable(
    abivarname="chempot",
    varset="geo",
    vartype="real",
    topics=['Artificial_expert'],
    dimensions=[3, '[[nzchempot]]', '[[ntypat]]'],
    defaultval=0.0,
    mnemonics="spatially varying CHEMical POTential",
    requires="[[nzchempot]] /= 0",
    added_in_version="before_v9",
    text=r"""
For each type of atoms, from 1 to [[ntypat]], specifies the spatially varying
chemical potential, through the specification of [[nzchempot]] triplets of
real numbers. They give data for [[nzchempot]] delimiting planes, all parallel
to each other, each determined by its z reduced coordinate.

The first real number is the z reduced coordinate of the delimiting plane.
The second real number is the value of the chemical potential for this type of atom on this
plane. The third real number is the derivative of the chemical potential for
this type of atom with respect to the z reduced coordinate, evaluated on this plane.

In the space between delimiting planes, a piecewise cubic polynomial
interpolation is determined: the cubic polynomial between two delimiting
planes will have the imposed chemical potentials and derivatives on the two
delimiting planes. The z reduced coordinates must be ordered in increasing
values, and cannot span more than 1.0. There is an automatic periodic
boundary condition imposed. Specifying two identical z reduced coordinates is
allowed, and means that the first one applies to the adjacent space with lower
values of z, while the second applies to the adjacent space with higher values
of z. When the spatial chemical potential is defined only for one type of atom
(and no chemical potential is present for the other atoms), simply set the
related values to *0.0 in the [[chempot]] array. In the present input array,
reduced positions, energies and derivatives of energies are mixed. Hence,
although the chemical potential is an energy, one cannot use the usual energy
definitions (i.e. the chemical potential is always to be input in Hartree atomic units).
""",
),

Variable(
    abivarname="chkdilatmx",
    varset="rlx",
    vartype="integer",
    topics=['GeoOpt_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="CHecK DILATMX",
    added_in_version="before_v9",
    text=r"""
If 0, the code will not stop execution if the threshold of [[dilatmx]] is exceeded,
it will simply issue a warning. There will be no rescaling. If 1, after tentative
rescaling as described in [[dilatmx]], the code will stop execution.
Also, the use of [[chkdilatmx]] = 0 allows one to set [[dilatmx]] to a larger value than 1.15,
otherwise forbidden as being a waste of CPU and memory.
""",
),

Variable(
    abivarname="chkexit",
    varset="gstate",
    vartype="integer",
    topics=['Control_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="CHecK whether the user want to EXIT",
    added_in_version="before_v9",
    text=r"""
If [[chkexit]] is 1 or 2, ABINIT will check whether the user wants to
interrupt the run (using the keyword "exit" on the top of the input file or
creating a file named "abinit.exit": see the end of section 3.2
of the [[help:abinit#parameters]]).

  * 0 --> the check is not performed at all
  * 1 --> the check is not performed frequently (after each SCF step)
  * 2 --> the check is performed frequently (after a few bands, at each k point)

In all cases, the check is performed at most every 2 seconds of CPU time.
""",
),

Variable(
    abivarname="chkprim",
    varset="gstate",
    vartype="integer",
    topics=['crystal_useful', 'UnitCell_useful', 'SmartSymm_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="CHecK whether the cell is PRIMitive",
    added_in_version="before_v9",
    text=r"""
If the symmetry finder is used (see [[nsym]]), a non-zero value of [[chkprim]]
will make the code stop if a non-primitive cell is used. If [[chkprim]] = 0, a
warning is issued, but the run does not stop.

If you are generating the atomic and cell geometry using [[spgroup]], you
might generate a PRIMITIVE cell using [[brvltt]] = -1.
""",
),

Variable(
    abivarname="chksymbreak",
    varset="gstate",
    vartype="integer",
    topics=['k-points_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="CHecK SYMmetry BREAKing",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
This variable governs the behaviour of the code when there are potential
source of symmetry breaking, related e.g. to the k point grid or the presence
of non-symmorphic translations which might not be coherent with the exchange-correlation grid.

When **chksymbreak** = 1, the code stops (or issue a warning) if:

  * (1) The k point grid is non-symmetric, in case [[kptopt]] =1, 2, or 4;
  * (2) The non-symmorphic translation part of the symmetry operations has components that are not zero,
    or simple fractions, with 2, 3, 4, 6, 8 or 12 as denominators.

Note that the check is disabled when the number of k-points in the BZ is greater than 40 ** 3.

When **chksymbreak** = 0, there is no such check.

When **chksymbreak** = -1, the code stops if the condition (1) is met,
but in case the condition (2) is met, there will be a trial to shift the
atomic coordinates such as to obtain symmetry operations with the adequate non-symmorphic part.

Explanation:
In the ground-state calculation, such breaking of the symmetry is usually
harmless. However, if the user is doing a calculation of phonons using DFPT
([[rfphon]] = 1), the convergence with respect to the number of k points will be
much worse with a non-symmetric grid than with a symmetric one. Also, if the
user is doing a GW calculation, the presence of non-symmorphic translations
that are not coherent with the FFT grid might cause problems. In the GW part,
indeed, one needs to reconstruct the wavefunctions in the full Brillouin zone
for calculating both the polarizability and the self-energy. The wavefunctions
in the full Brillouin zone are obtained from the irreducible wedge by applying
the symmetry operations of the space group of the crystal. In the present
implementation, the symmetrisation of the wavefunctions is done in real space
on the FFT mesh that, therefore, has to be coherent both with the rotational
part as well as with the fractional translation of each symmetry operation. If
the condition (2) is met, the GW code will not be able to find a symmetry
preserving FFT mesh.

So, it was decided to warn the user about these possible problems already at
the level of the ground state calculations, although such warning might be irrelevant.

If you encounter a problem outlined above, you have two choices: change your
atomic positions (translate them) such that the origin appears as the most
symmetric point; or ignore the problem, and set **chksymbreak** = 0.
""",
),

Variable(
    abivarname="chneut",
    varset="eph",
    vartype="integer",
    topics=['Phonons_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="CHarge NEUTrality treatment",
    added_in_version="before_v9",
    text=r"""
Set the treatment of the Charge Neutrality requirement for the effective charges.
Same meaning as the corresponding anaddb variable.
Note the different default value in abinit and anaddb
""",
),

Variable(
    abivarname="cineb_start",
    varset="rlx",
    vartype="integer",
    topics=['TransPath_expert'],
    dimensions="scalar",
    defaultval=7,
    mnemonics="Climbing-Image Nudged Elastic Band: STARTing iteration",
    requires="[[imgmov]] == 5 and [[neb_algo]] == 2",
    added_in_version="before_v9",
    text=r"""
Gives the index of the first CI-NEB iteration.
The CI-NEB method constitutes a small modification to the NEB method allowing
a rigorous convergence to the saddle point. As the image with the highest
energy has to be identified, the calculation begins with several iterations of
the standard NEB algorithm. The effective CI-NEB begins at the [[cineb_start]]
iteration. See [[cite:Henkelman2000a]] for additional details of this method.
""",
),

Variable(
    abivarname="constraint_kind",
    varset="gstate",
    vartype="integer",
    topics=['ConstrainedDFT_basic'],
    dimensions=['[[ntypat]]'],
    defaultval=0,
    mnemonics="CONSTRAINT KIND in constrained DFT",
    requires="[[iscf]] > 1 and [[iscf]] < 10 and [[ionmov]] /= 4",
    added_in_version="before_v9",
    text=r"""
If [[constraint_kind]] is non-zero for at least one type of atom,
the constrained DFT algorithm is activated.
[[constraint_kind]] defines, for each type of atom, the kind of constraint(s) imposed by constrained DFT.
When [[constraint_kind]] is zero for an atom type, there is not constraint applied to this atom type.
Otherwise, different constraints can be imposed on the total charge (ion+electronic) and/or magnetization, computed
inside a sphere of radius [[ratsph]], possibly smeared within a width [[ratsm]].
Such integrated ion+electronic charge might be imposed to be equal to [[chrgat]], while the magnetization might be compared to [[spinat]].
The first digit of [[constraint_kind]] defines the constraint on the charge, while the second digit defines the constraint on the
magnetization.

When [[constraint_kind]] is 10 or above, the charge constraint will be imposed.

When [[constraint_kind]]=1 or 11, the exact value (vector in the non-collinear case, amplitude and sign in the collinear case) of the magnetization is constrained;
When [[constraint_kind]]=2 or 12, only the direction is constrained (only meaningful in the non-collinear case);
When [[constraint_kind]]=3 or 13, only the magnitude is constrained.

For the algorithm, see [[topic:ConstrainedDFT]]. It makes important use of the potential residual,
so the algorithm works only with [[iscf]] between 2 and 9.
The balance between the potential residual, and the density/magnetization constraint is governed by [[magcon_lambda]]. The spherical integral is governed by [[ratsph]] and [[ratsm]].

Note that while a spherical integral around an atom might reasonably well capture the magnetization of an atom within a solid or within a molecule,
 so that the sum of such magnetizations might be reasonably close to the total magnetization of the solid,
such a procedure hardly gives the total charge of the solid: the space between the spheres is too large when the spheres do not overlap,
while overlapping spheres will not deliver the correct total charge of the system.

Note that [[constraint_kind]] defines constraints for types of atoms, not for specific atoms.
Atoms of the same type are supposed to incur the same constraint.
If the user wants to impose different constraints on atoms of the same type (in principle), it is possible (and easy) to pretend
that they belong to different types, even if the same pseudopotential file is used for these atoms. There is an example
in test [[test:v8_24]], the hydrogen dimer, where the charge around the first atom is constrained, and the charge around the second atom is left free.

Incidentally, [[ionmov]]==4 is not allowed in the present implementation of constrained DFT because the motion of atoms and simultaneous computation of constraints would be difficult to handle.
""",
),


Variable(
    abivarname="cpuh",
    varset="gstate",
    vartype="real",
    topics=['Control_useful'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="CPU time limit in Hours",
    characteristics=['[[NO_MULTI]]', '[[INPUT_ONLY]]'],
    excludes="[[cpum]] or [[cpus]]",
    added_in_version="before_v9",
    text=r"""
Only one of the three real parameters [[cpus]], [[cpum]] and [[cpuh]] can be
defined in the input file to set up a CPU time limit. When the job reaches
that limit, it will try to end smoothly. However, note that this might still
take some time. If the user want a firm CPU time limit, the present parameter
must be reduced sufficiently. Intuition about the actual margin to be taken
into account should come with experience.
A zero value has no action of the job.

!!! tip

    One can pass the timelimit to abinit via the command line option:

        abinit --timelimit hours:minutes:seconds

    This approach is much more powerful especially when the job must be submitted
    to the queue via a submission script e.g. a Slurm script.
    In this case, indeed, one can define a shell variable for the time limit
    and use this variable to pass the time limit to Slurm and Abinit at the same time.

    Use `abinit --help` for further information.
""",
),

Variable(
    abivarname="cpum",
    varset="gstate",
    vartype="real",
    topics=['Control_useful'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="CPU time limit in Minutes",
    characteristics=['[[NO_MULTI]]', '[[INPUT_ONLY]]'],
    excludes="[[cpuh]] or [[cpus]]",
    added_in_version="before_v9",
    text=r"""
Only one of the three real parameters [[cpus]], [[cpum]] and [[cpuh]] can be
defined in the input file to set up a CPU time limit. When the job reaches
that limit, it will try to end smoothly. However, note that this might still
take some time. If the user want a firm CPU time limit, the present parameter
must be reduced sufficiently. Intuition about the actual margin to be taken
into account should come with experience.
A zero value has no action of the job.
""",
),

Variable(
    abivarname="cpus",
    varset="gstate",
    vartype="real",
    topics=['Control_useful'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="CPU time limit in seconds",
    characteristics=['[[NO_MULTI]]', '[[INPUT_ONLY]]'],
    excludes="[[cpuh]] or [[cpum]]",
    added_in_version="before_v9",
    text=r"""
Only one of the three real parameters [[cpus]], [[cpum]] and [[cpuh]] can be
defined in the input file to set up a CPU time limit. When the job reaches
that limit, it will try to end smoothly. However, note that this might still
take some time. If the user want a firm CPU time limit, the present parameter
must be reduced sufficiently. Intuition about the actual margin to be taken
into account should come with experience.
A zero value has no action of the job.
""",
),

Variable(
    abivarname="dvdb_qcache_mb",
    varset="eph",
    vartype="real",
    topics=['ElPhonInt_useful'],
    dimensions="scalar",
    defaultval=1024,
    mnemonics="DVDB Q-CACHE size in Megabytes",
    added_in_version="before_v9",
    text=r"""
This variable activates a caching mechanism for the DFPT potentials used in the EPH part.
The code will store in memory multiple q-points up to this size in Megabytes in order
to reduce the number of IO operations required to read the potentials from the DVDB file.

This option leads to a **significant speedup** of calculations requiring integrations
in q-space ([[eph_task]] == 4) at the price of an increase of the memory requirements.
The speedup is important especially if the QP corrections are computed for several k-points.

A negative value signals to the code that all the q-points in the DVDB should be stored in memory.
Use zero value disables the cache.
""",
),

Variable(
    abivarname="d3e_pert1_atpol",
    varset="dfpt",
    vartype="integer",
    topics=['nonlinear_basic'],
    dimensions=[2],
    defaultval=[1, 1],
    mnemonics="3rd Derivative of Energy, mixed PERTurbation 1: limits of ATomic POLarisations",
    requires="[[optdriver]] == 5 (non-linear response computations)",
    added_in_version="before_v9",
    text=r"""
Controls the range of atoms for which displacements will be considered in non-linear
computations (using the 2n+1 theorem), for the 1st perturbation.
May take values from 1 to [[natom]], with **d3e_pert1_atpol** (1)<=
**d3e_pert1_atpol** (2). See [[rfatpol]] for additional details.
""",
),

Variable(
    abivarname="d3e_pert1_dir",
    varset="dfpt",
    vartype="integer",
    topics=['nonlinear_basic'],
    dimensions=[3],
    defaultval=[0, 0, 0],
    mnemonics="3rd Derivative of Energy, mixed PERTurbation 1: DIRections",
    requires="[[optdriver]] == 5 (non-linear response computations)",
    added_in_version="before_v9",
    text=r"""
Gives the directions to be considered in non-linear computations (using the
2n+1 theorem), for the 1st perturbation.
The three elements corresponds to the three primitive vectors, either in real
space (atomic displacement), or in reciprocal space (electric field perturbation).
See [[rfdir]] for additional details.
""",
),

Variable(
    abivarname="d3e_pert1_elfd",
    varset="dfpt",
    vartype="integer",
    topics=['nonlinear_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="3rd Derivative of Energy, mixed PERTurbation 1: ELectric FielD",
    requires="[[optdriver]] == 5 (non-linear response computations)",
    added_in_version="before_v9",
    text=r"""
Turns on electric field perturbation in non-linear computation, as 1st
perturbation. Actually, such calculations requires first the non-self-
consistent calculation of derivatives with respect to k, independently of the
electric field perturbation itself. See [[rfelfd]] for additional details.
""",
),

Variable(
    abivarname="d3e_pert1_phon",
    varset="dfpt",
    vartype="integer",
    topics=['nonlinear_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="3rd Derivative of Energy, mixed PERTurbation 1: PHONons",
    requires="[[optdriver]] == 5 (non-linear response computations)",
    added_in_version="before_v9",
    text=r"""
Turns on atomic displacement perturbation in non-linear computation, as 1st perturbation.
See [[rfphon]] for additional details.
""",
),

Variable(
    abivarname="d3e_pert2_atpol",
    varset="dfpt",
    vartype="integer",
    topics=['nonlinear_basic'],
    dimensions=[2],
    defaultval=[1, 1],
    mnemonics="3rd Derivative of Energy, mixed PERTurbation 2: limits of ATomic POLarisations",
    requires="[[optdriver]] == 5 (non-linear response computations)",
    added_in_version="before_v9",
    text=r"""
Controls the range of atoms for which displacements will be considered in non-
linear computations (using the 2n+1 theorem), for the 2nd perturbation.
May take values from 1 to [[natom]], with **d3e_pert2_atpol** (1) <= **d3e_pert2_atpol** (2).
See [[rfatpol]] for additional details.
""",
),

Variable(
    abivarname="d3e_pert2_dir",
    varset="dfpt",
    vartype="integer",
    topics=['nonlinear_basic'],
    dimensions=[3],
    defaultval=[0, 0, 0],
    mnemonics="3rd Derivative of Energy, mixed PERTurbation 2: DIRections",
    requires="[[optdriver]] == 5 (non-linear response computations)",
    added_in_version="before_v9",
    text=r"""
Gives the directions to be considered in non-linear computations (using the
2n+1 theorem), for the 2nd perturbation.
The three elements corresponds to the three primitive vectors, either in real
space (atomic displacement), or in reciprocal space (electric field perturbation).
See [[rfdir]] for additional details.
""",
),

Variable(
    abivarname="d3e_pert2_elfd",
    varset="dfpt",
    vartype="integer",
    topics=['nonlinear_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="3rd Derivative of Energy, mixed PERTurbation 2: ELectric FielD",
    requires="[[optdriver]] == 5 (non-linear response computations)",
    added_in_version="before_v9",
    text=r"""
Turns on electric field perturbation in non-linear computation, as 2nd
perturbation. Actually, such calculations requires first the non-self-
consistent calculation of derivatives with respect to k, independently of the
electric field perturbation itself. See [[rfelfd]] for additional details.
""",
),

Variable(
    abivarname="d3e_pert2_phon",
    varset="dfpt",
    vartype="integer",
    topics=['nonlinear_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="3rd Derivative of Energy, mixed PERTurbation 2: PHONons",
    requires="[[optdriver]] == 5 (non-linear response computations)",
    added_in_version="before_v9",
    text=r"""
Turns on atomic displacement perturbation in non-linear computation, as 2nd perturbation.
See [[rfphon]] for additional details.
""",
),

Variable(
    abivarname="d3e_pert3_atpol",
    varset="dfpt",
    vartype="integer",
    topics=['nonlinear_basic'],
    dimensions=[2],
    defaultval=[1, 1],
    mnemonics="3rd Derivative of Energy, mixed PERTurbation 3: limits of ATomic POLarisations",
    requires="[[optdriver]] == 5 (non-linear response computations)",
    added_in_version="before_v9",
    text=r"""
Controls the range of atoms for which displacements will be considered in non-
linear computations (using the 2n+1 theorem), for the 3rd perturbation.
May take values from 1 to [[natom]], with **d3e_pert3_atpol** (1)<= **d3e_pert3_atpol** (2).
See [[rfatpol]] for additional details.
""",
),

Variable(
    abivarname="d3e_pert3_dir",
    varset="dfpt",
    vartype="integer",
    topics=['nonlinear_basic'],
    dimensions=[3],
    defaultval=[0, 0, 0],
    mnemonics="3rd Derivative of Energy, mixed PERTurbation 3: DIRections",
    requires="[[optdriver]] == 5 (non-linear response computations)",
    added_in_version="before_v9",
    text=r"""
Gives the directions to be considered in non-linear computations (using the
2n+1 theorem), for the 3rd perturbation.
The three elements corresponds to the three primitive vectors, either in real
space (atomic displacement), or in reciprocal space (electric field perturbation).
See [[rfdir]] for additional details.
""",
),

Variable(
    abivarname="d3e_pert3_elfd",
    varset="dfpt",
    vartype="integer",
    topics=['nonlinear_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="3rd Derivative of Energy, mixed PERTurbation 3: ELectric FielD",
    requires="[[optdriver]] == 5 (non-linear response computations)",
    added_in_version="before_v9",
    text=r"""
Turns on electric field perturbation in non-linear computation, as 3rd
perturbation. Actually, such calculations requires first the non-self-
consistent calculation of derivatives with respect to k, independently of the
electric field perturbation itself.
See [[rfelfd]] for additional details.
""",
),

Variable(
    abivarname="d3e_pert3_phon",
    varset="dfpt",
    vartype="integer",
    topics=['nonlinear_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="3rd Derivative of Energy, mixed PERTurbation 3: PHONons",
    requires="[[optdriver]] == 5 (non-linear response computations)",
    added_in_version="before_v9",
    text=r"""
Turns on atomic displacement perturbation in non-linear computation, as 3rd perturbation.
See [[rfphon]] for additional details.
""",
),

Variable(
    abivarname="ddamp",
    varset="ffield",
    vartype="real",
    topics=['Berry_useful'],
    dimensions="scalar",
    defaultval=0.1,
    mnemonics="electric Displacement field DAMPing parameter",
    requires="[[berryopt]] in [6, 16]",
    added_in_version="before_v9",
    text=r"""
In case [[berryopt]] = 6, the electric field is updated after each SCF iteration
according to $E_{n+1}=$[[ddamp]]$(D-4 \pi P_{n})+(1-$[[ddamp]]$)E_{n}$, where
$P_{n}$ and $E_{n}$ are the polarization and electric field after $n_{th}$ SCF
iteration. [[ddamp]] is a damping parameter used to control the convergence speed.
In case [[berryopt]] = 16, the electric field is updated after each SCF
iteration according to $e_{n+1}=$[[ddamp]]$(d-p_{n})+(1-$[[ddamp]]$)e_{n}$.
If you have difficulty getting convergence, try to reduce this value or reduce
maxestep. This parameter is used in finite electric displacement field
calculations (berryopt=6,16,17).
""",
),

Variable(
    abivarname="ddb_ngqpt",
    varset="eph",
    vartype="integer",
    topics=['ElPhonInt_basic'],
    dimensions=[3],
    defaultval=[0, 0, 0],
    mnemonics="Derivative DataBase: Number of Grid points for Q-PoinTs",
    added_in_version="before_v9",
    text=r"""
This variable is mandatory when [[optdriver]] == 7. It defines the number of
divisions in the (homogeneous) q-mesh used to generate the DDB file. See also
the description of the [[getddb]], [[getddb_path]] input variables.
""",
),

Variable(
    abivarname="ddb_shiftq",
    varset="eph",
    vartype="real",
    topics=['ElPhonInt_expert'],
    dimensions=[3],
    defaultval=[0.0, 0.0, 0.0],
    mnemonics="Derivative DataBase: SHIFT of the Q-points",
    added_in_version="before_v9",
    text=r"""
Only relevant when [[optdriver]] == 7. It defines the shift in the q-mesh used
to generate the DDB file, which is defined by the [[ddb_ngqpt]] input
variable. See [[shiftk]] for more information on the definition.
""",
),

Variable(
    abivarname="delayperm",
    varset="rlx",
    vartype="integer",
    topics=['MolecularDynamics_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="DELAY between trials to PERMUTE atoms",
    added_in_version="before_v9",
    text=r"""
Delay (number of time steps) between trials to permute two atoms, in view of
accelerated search of minima. Still in development.

See the routine moldyn.F90 and [[signperm]] for additional information.

When [[delayperm]] is zero, there are no permutation trials.
""",
),

Variable(
    abivarname="chrgat",
    varset="gstate",
    vartype="real",
    topics=['ConstrainedDFT_useful'],
    dimensions=ValueWithConditions({'[[natrd]]<[[natom]]': '[ [[natrd]] ]', 'defaultval': '[ [[natom]] ]'}),
    defaultval=0.0,
    mnemonics="CHARGE of the AToms",
    added_in_version="before_v9",
    text=r"""
Gives the target integrated charge in case of constrained DFT calculations, see [[constraint_kind]].
Given in atomic unit of charge (=minus the charge of the electron).
Note that this number is the net positive charge inside the sphere: one subtract from the
nucleus charge [[ziontypat]] the integrated valence electron density in a sphere defined by [[ratsph]].
The latter has indeed a negative value. Note that if the sphere radius [[ratsph]] is not sufficiently large,
the amount of electrons will be smaller than expected based on chemical intuition. This means that there
is in this case a bias toward too positive integrated charges. By contrast, if the sphere radius is too large,
the spheres will overlap, and the electrons in the interatomic region will be double counted.
""",
),

Variable(
    abivarname="densfor_pred",
    varset="dev",
    vartype="integer",
    topics=['SCFAlgorithms_expert', 'MolecularDynamics_expert'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[paral_kgb]] == 1': '6', 'defaultval': 2}),
    mnemonics="DENSity and FORces PREDictor",
    characteristics=['[[DEVELOP]]'],
    requires="[[iscf]] > 0",
    added_in_version="before_v9",
    text=r"""
Used when [[iscf]] > 0, to define:

  - the way a change of density is derived from a change of atomic position,
  - the way forces are corrected when the SCF cycle is not converged.

Supported values:

  * 0 --> density not changed (fixed charge), forces not corrected
  * 1 --> density not changed, forces corrected with rigid ion hypothesis (atomic charge moved with atom)
  * 2 --> density changed and forces corrected with rigid ion hypothesis (atomic charge moves with atom)
  * 3 --> density changed and forces corrected with a different implementation of the rigid ion hypothesis
  * 4 --> density not changed, forces corrected with the use of Harris functional formula (see note)
  * 5 --> density changed using D. Alfe 2nd-order algorithm (see notes), forces not corrected
  * 6 --> density changed using D. Alfe 2nd-order algorithm (see notes) and forces corrected with the use of Harris functional formula

Similar negative values are also allowed (see the meaning later), for
development purposes only. No meaning for RF calculations.

For the time being,

  - [[densfor_pred]] = 3 must be used with [[ionmov]] = 4 and [[iscf]] = 5.
  - [[densfor_pred]] = 4, 5 or 6 must be used when band-FFT parallelism is selected.
  Otherwise, use [[densfor_pred]] = 2


!!! note "concerning the correction of forces (use of [[densfor_pred]] = 1, 2, 3, 4 or 6)"

    The force on the atom located at R is corrected by the addition of the following
    term: $F_{residual}=\int dr V_{residual} \frac{d \rho_{atomic}}{dR}$,
    where $\rho_{atomic}$ is an atomic (spherical) density.

    - When such an atomic density ($\rho_{atomic}$) is found in the pseudopotential or
    PAW file, it is used. If not, a gaussian density (defined by [[densty]] parameter) is used.
    - When SCF mixing is done on the density ([[iscf]] >= 10), the potential
    residual ($V_residual$) is obtained from the density residual with the first
    order formula $V_{residual}=\frac{dV}{d \rho} \rho_{residual}$
    and uses the exchange-correlation kernel
    $ \frac{dV_{xc}}{d\rho}=K_{xc}$ whose computation is time-consuming for GGA
    functionals. By default (positive values of [[densfor_pred]]), the local-
    density part of the GGA exchange-correlation kernel is used (even for GGA, for
    which it seems to give a reasonable accuracy). Using the full GGA exchange
    correlation kernel (so, including derivatives with respect to the gradient of
    the density) is always possible by giving a negative value to
    [[densfor_pred]]. In case of hybrid functionals, a similar correction term is
    added, although in the density mixing scheme, the related GGA kernel is used
    instead of the hybrid functional kernel.

!!! note "concerning the use of [[densfor_pred]] = 5 or 6 (density prediction)"

    The algorithm is described in [[cite:Alfe1999]].
    It uses an atomic (spherical) density. When such an atomic density
    is found in the pseudopotential or PAW file, it is used. If not, a gaussian
    density (defined by [[densty]] parameter) is used.
    Also note that, to be efficient, this algorithm requires a minimum convergence
    of the SCF cycle; Typically, vres2 (or nres2) has to be small enough (let's say smaller than 10e-4).
""",
),

Variable(
    abivarname="densty",
    varset="dev",
    vartype="real",
    topics=['xc_expert'],
    dimensions=['[[ntypat]]'],
    defaultval=0.0,
    mnemonics="initial DENSity for each TYpe of atom",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Gives a rough description of the initial GS density, for each type of atom.
This value is only used to create the first exchange and correlation
potential, and is not used anymore afterwards. For the time being, it
corresponds to an average radius (a.u.) of the density, and is used to
generate a gaussian density. If set to 0.0, an optimized value is used.
No meaning for RF calculations.
""",
),

Variable(
    abivarname="dfield",
    varset="ffield",
    vartype="real",
    topics=['Berry_basic'],
    dimensions=[3],
    defaultval=MultipleValue(number=3, value=0.0),
    mnemonics="Displacement FIELD",
    requires="[[berryopt]] == 6 and [[efield]]",
    added_in_version="before_v9",
    text=r"""
In case [[berryopt]] = 6, [[dfield]] specifies the (unreduced) finite electric
displacement field vector, in atomic units, that is to be imposed as a
constraint during the calculation.
""",
),

Variable(
    abivarname="dfpt_sciss",
    varset="dfpt",
    vartype="real",
    topics=['DFPT_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="DFPT SCISSor operator",
    characteristics=['[[ENERGY]]'],
    added_in_version="before_v9",
    text=r"""
It is the value of the "scissors operator", the shift of conduction band
eigenvalues, used in response function calculations.
Can be specified in Ha (the default), Ry, eV or Kelvin, since [[ecut]] has the
[[ENERGY]] characteristics (1 Ha = 27.2113845 eV).
Typical use is for response to electric field ([[rfelfd]] = 3), but NOT for d/dk
([[rfelfd]] = 2) and phonon responses.
""",
),

Variable(
    abivarname="diecut",
    varset="gstate",
    vartype="real",
    topics=['SCFAlgorithms_expert'],
    dimensions="scalar",
    defaultval=2.2,
    mnemonics="DIElectric matrix energy CUToff",
    characteristics=['[[ENERGY]]'],
    added_in_version="before_v9",
    text=r"""
Kinetic energy cutoff that controls the number of planewaves used to represent
the dielectric matrix:
$(1/2) [ 2 \pi \GG_{diel,max}]^2$ =[[diecut]] with $\GG_{diel,max}$ being the maximum
length of the reciprocal space planewave wavevectors for the dielectric matrix.
Can be specified in Ha (the default), Ry, eV or Kelvin, since [[diecut]] has
the [[ENERGY]] characteristics. (1 Ha = 27.2113845 eV)
All planewaves inside this "basis sphere" centered at $\GG$=0 are included in the
basis. This is useful only when [[iprcel]] > =21, which means that a
preconditioning scheme based on the dielectric matrix is used.

NOTE: a negative [[diecut]] will define the same dielectric basis sphere as
the corresponding positive value, but the FFT grid will be identical to the
one used for the wavefunctions. The much smaller FFT grid, used when
[[diecut]] is positive, gives exactly the same results.

No meaning for RF calculations yet.
""",
),

Variable(
    abivarname="diegap",
    varset="gstate",
    vartype="real",
    topics=['SCFAlgorithms_expert'],
    dimensions="scalar",
    defaultval=0.1,
    mnemonics="DIElectric matrix GAP",
    characteristics=['[[ENERGY]]'],
    added_in_version="before_v9",
    text=r"""
Gives a rough estimation of the dielectric gap between the highest energy
level computed in the run, and the set of bands not represented. Used to
extrapolate dielectric matrix when [[iprcel]] >= 21.
Can be specified in Ha (the default), Ry, eV or Kelvin, since [[diegap]] has
the [[ENERGY]] characteristics (1 Ha = 27.2113845 eV).

No meaning for RF calculations yet.
""",
),

Variable(
    abivarname="dielam",
    varset="gstate",
    vartype="real",
    topics=['SCFAlgorithms_expert'],
    dimensions="scalar",
    defaultval=0.5,
    mnemonics="DIElectric matrix LAMbda",
    requires="[[iprcel]] >= 21",
    added_in_version="before_v9",
    text=r"""
Gives the amount of occupied states with mean energy given by the highest
level computed in the run, included in the extrapolation of the dielectric matrix.

No meaning for RF calculations yet.
""",
),

Variable(
    abivarname="dielng",
    varset="gstate",
    vartype="real",
    topics=['SCFAlgorithms_useful'],
    dimensions="scalar",
    defaultval="1.0774841",
    mnemonics="model DIElectric screening LeNGth",
    characteristics=['[[LENGTH]]'],
    added_in_version="before_v9",
    text=r"""
Used for screening length (in Bohr) of the model dielectric function, diagonal
in reciprocal space. By default, given in Bohr atomic units (1
Bohr=0.5291772108 Angstrom), although Angstrom can be specified, if preferred,
since [[dielng]] has the [[LENGTH]] characteristics.
This model dielectric function is as follows (${\bf K}$ being a wavevector):

\begin{equation}
\epsilon({\bf K}) = \frac{ 1 + [[dielng]]^2 {\bf K}^2 }{ \left( 1/[[diemac]] + [[dielng]]^2 {\bf K}^2  \right) [[diemix]] } \nonumber
\end{equation}

The inverse of this model dielectric function will be applied to the residual,
to give the preconditioned change of potential. Right at ${\bf K}$=0, $\epsilon({\bf K})$ is imposed to be 1.

If the preconditioning were perfect, the change of potential would lead to an
exceedingly fast solution of the self-consistency problem (two or three
steps). The present model dielectric function is excellent for rather
homogeneous unit cells.
When ${\bf K}$->0, it tends to the macroscopic dielectric constant, eventually
divided by the mixing factor [[diemix]] (or [[diemixmag]]  for magnetization).
For metals, simply put [[diemac]] to a very large value ($10^6$ is OK)
The screening length [[dielng]] governs the length scale to go from the
macroscopic regime to the microscopic regime, where it is known that the
dielectric function should tend to 1. It is on the order of 1 Bohr for metals
with medium density of states at the Fermi level, like Molybdenum, and for
Silicon. For metals with a larger DOS at the Fermi level (like Iron), the
screening will be more effective, so that [[dielng]] has to be decreased by a factor of 2-4.
This works for GS and RF calculations.
""",
),

Variable(
    abivarname="diemac",
    varset="gstate",
    vartype="real",
    topics=['SCFAlgorithms_useful'],
    dimensions="scalar",
    defaultval=1000000.0,
    mnemonics="model DIElectric MACroscopic constant",
    added_in_version="before_v9",
    text=r"""
A rough knowledge of the macroscopic dielectric constant [[diemac]] of the
system is a useful help to speed-up the SCF procedure: a model dielectric
function, see the keyword [[dielng]], is used for that purpose. It is
especially useful for speeding up the treatment of rather homogeneous unit cells.

Some hint:
The value of [[diemac]] should usually be bigger than 1.0, on physical grounds.

  * For metals, simply put [[diemac]] to a very large value (the default $10^6$ is OK)
  * For silicon, use 12.0. A similar value is likely to work well for other semiconductors
  * For wider gap insulators, use 2.0 ... 4.0
  * For molecules in an otherwise empty big box, try 1.5 ... 3.0

Systems that combine a highly polarisable part and some vacuum are rather
badly treated by the model dielectric function. One has to use the
"extrapolar" technique, activated by the input variable [[iprcel]].
In sufficiently homogeneous systems, you might have to experiment a bit to
find the best [[diemac]]. If you let [[diemac]] to its default value, you
might even never obtain the self-consistent convergence!
For response function calculations, use the same values as for GS. The
improvement in speed can be considerable for small (but non-zero) values of the wavevector.
""",
),

Variable(
    abivarname="diemix",
    varset="gstate",
    vartype="real",
    topics=['SCFAlgorithms_useful'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[usepaw]] == 0 or [[iprcel]] !=0': 1.0,
 '[[usepaw]] == 1 or [[iprcel]] == 0': 0.7,
 'defaultval': None}),
    mnemonics="model DIElectric MIXing factor",
    requires="[[diemix]] >= 0.0 and [[diemix]] <=  1.0",
    added_in_version="before_v9",
    text=r"""
Gives overall factor of the preconditioned residual density/potential to be
transferred in the SCF cycle.
It should be between 0.0 and 1.0.

If the model dielectric function were perfect, [[diemix]] should be 1.0. By
contrast, if the model dielectric function does nothing (when [[diemac]] = 1.0
or [[dielng]] is larger than the size of the cell), [[diemix]] can be used to
damp the amplifying factor inherent to the SCF loop.
For molecules, a value on the order 0.5 or 0.33 is rather usual.
When mod([[iscf]],10)=3, 4,5 or 7, [[diemix]] is only important at the few
first iterations when anharmonic effects are important, since these schemes
compute their own mixing factor for self-consistency.
Also note that a different value of diemix can be used for the magnetization (see [[diemixmag]]).
""",
),

Variable(
    abivarname="diemixmag",
    varset="gstate",
    vartype="real",
    topics=['spinpolarisation_useful'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'70 < [[iprcel]] and [[iprcel]] < 80': '[[diemix]]',
 '[[iprcel]] == 0': '[[diemix]]',
 '[[iscf]]<10': '[[diemix]]',
 'defaultval': '-[[diemix]]'}),
    mnemonics="model DIElectric MIXing factor for the MAGgnetization",
    added_in_version="before_v9",
    text=r"""
Gives overall factor of the preconditioned residual magnetization/magnetic
field to be transferred in the SCF cycle (see [[diemix]] for further
information).
For the time being, apply only when the SCF mixing is done on the density ([[iscf]] > =10).

A negative value of [[diemixmag]] means that magnetization is only preconditioned
by ABS([[diemixmag]]), without the use of any preconditioner.

When SCF cycle has some difficulties to converge, changing the value of
[[diemixmag]] can have a positive effect.
In particular [[diemixmag]] = -4 is a good choice (i.e. [[diemixmag]] = 4, no other
preconditioner on magnetization).
""",
),

Variable(
    abivarname="diismemory",
    varset="rlx",
    vartype="integer",
    topics=['MolecularDynamics_expert'],
    dimensions="scalar",
    defaultval=8,
    mnemonics="Direct Inversion in the Iterative Subspace MEMORY",
    added_in_version="before_v9",
    text=r"""
Gives the maximum number of "time" steps for which the forces and stresses are
stored, and taken into account in the DIIS algorithm ([[ionmov]] = 20) to find
zero-force and stress configurations.
""",
),

Variable(
    abivarname="dilatmx",
    varset="rlx",
    vartype="real",
    topics=['GeoOpt_basic'],
    dimensions="scalar",
    defaultval=1.0,
    mnemonics="lattice DILATation: MaXimal value",
    added_in_version="before_v9",
    text=r"""
[[dilatmx]] is an auxiliary variable used to book additional memory (see detailed description later) for possible
on-the-flight variations the plane wave basis set, due to cell optimization by ABINIT.
Useful only when [[ionmov]] == 2 or 22 and [[optcell]]/=0, that is, cell optimization.

In the default mode ([[chkdilatmx]] = 1), when the [[dilatmx]] threshold is exceeded,
ABINIT will rescale uniformly the
tentative new primitive vectors to a value that leads at most to 90% of the
maximal allowed [[dilatmx]] deviation from 1. It will do this three times (to
prevent the geometry optimization algorithms to have taken a too large trial
step), but afterwards will exit. Setting [[chkdilatmx]] == 0 allows one to
book a larger planewave basis, but will not rescale the tentative new primitive vectors
nor lead to an exit when the [[dilatmx]] threshold is exceeded.
The obtained optimized primitive vectors will not be exactly the ones corresponding to the planewave basis set
determined using [[ecut]] at the latter primitive vectors. Still, as an intermediate step in a geometry search
this might be sufficiently accurate. In such case, [[dilatmx]] might even be let at its default value 1.0.

Detailed explanation: The memory space for the planewave basis set is defined
by multiplying [[ecut]] by [[dilatmx]] squared (the result is an "effective ecut", called
internally "ecut_eff". Other uses of [[ecut]] are not modified when [[dilatmx]] > 1.0.
Still, operations (like scalar products) are taking into account these fake non-used planewaves,
thus slowing down the ABINIT execution.
Using [[dilatmx]]<1.0 is equivalent to changing [[ecut]] in all its uses. This
is allowed, although its meaning is no longer related to a maximal expected scaling.

Setting [[dilatmx]] to a large value leads to waste of CPU time and memory.
By default, ABINIT will not accept that you define [[dilatmx]] bigger than 1.15.
This behaviour will be overcome by using [[chkdilatmx]] == 0.
Supposing you think that the optimized [[acell]] values might be 10% larger
than your input values, use simply [[dilatmx]] 1.1. This will already lead to
an increase of the number of planewaves by a factor (1.1)  3  =1.331, and a
corresponding increase in CPU time and memory.
It is possible to use [[dilatmx]] when [[optcell]] =0, but a value larger than
1.0 will be a waste.
""",
),

Variable(
    abivarname="dipdip",
    varset="eph",
    vartype="integer",
    topics=['Phonons_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="DIPole-DIPole interaction",
    added_in_version="before_v9",
    text=r"""
This variable defines the treatment of the dipole-dipole interaction. Same
meaning as the corresponding anaddb variable [[dipdip@anaddb]]
""",
),

Variable(
    abivarname="dmatpawu",
    varset="paw",
    vartype="real",
    topics=['DFT+U_useful', 'ElecDOS_useful'],
    dimensions=['2*max([[lpawu]])+1',
 '2*max([[lpawu]])+1',
 'max([[nsppol]], [[nspinor]])',
 '[[natpawu]]'],
    defaultval=MultipleValue(number=None, value=-10.0),
    mnemonics="initial Density MATrix for PAW+U",
    requires="[[usepaw]] == 1 and [[usepawu]] == 1 and [[usedmatpu]] != 0",
    added_in_version="before_v9",
    text=r"""
For Ground state calculations only.
Gives the value of an initial density matrix used in LDA+U and kept fixed
during the first abs([[usedmatpu]]) SCF iterations.
Only components corresponding to [[lpawu]] angular momentum are requested.
Restriction: In order to use dmatpawu, [[lpawu]] must be identical for all atom types (or -1).

The occupation matrix is in the basis of real spherical harmonics Slm (note
that this differs from the choice made when [[prtdosm]] = 1, that is in the
basis of complex spherical harmonics). They are ordered by increasing m, and
are defined e.g. in [[cite:Blancoa1997]]. For the case l=2 (d states), the five
columns corresponds respectively to (the normalisation factor has been dropped)

  * m=-2, $xy$
  * m=-1, $yz$
  * m=0, $3z^{2}-r^{2}$
  * m=1, $xz$
  * m=2, $x^{2}-y^{2}$

[[dmatpawu]] must always be given as a "spin-up" occupation matrix (and
eventually a "spin-down" matrix). Be aware that its physical meaning depends
on the magnetic properties imposed to the system (with [[nsppol]],
[[nspinor]], [[nspden]]):

  * **Non-magnetic system** ([[nsppol]] = 1, [[nspinor]] = 1, [[nspden]] = 1):
One (2lpawu+1)x(2lpawu+1) [[dmatpawu]] matrix is given for each atom on which
+U is applied.
It contains the "spin-up" occupations.

  * **Ferromagnetic spin-polarized (collinear) system** ([[nsppol]] = 2, [[nspinor]] = 1, [[nspden]] = 2):
Two (2lpawu+1)x(2lpawu+1) [[dmatpawu]] matrices are given for each atom on
which +U is applied.
They contain the "spin-up" and "spin-down" occupations.

  * **Anti-ferromagnetic spin-polarized (collinear) system** ([[nsppol]] = 1, [[nspinor]] = 1, [[nspden]] = 2):
One (2lpawu+1)x(2lpawu+1) [[dmatpawu]] matrix is given for each atom on which +U is applied.
It contains the "spin-up" occupations.

  * **Non-collinear magnetic system** ([[nsppol]] = 1, [[nspinor]] = 2, [[nspden]] = 4):
Two (2lpawu+1)x(2lpawu+1) [[dmatpawu]] matrices are given for each atom on
which +U is applied.
They contains the "spin-up" and "spin-down" occupations (defined as
n_up=(n+|m|)/2 and n_dn=(n-|m|)/2), where m is the integrated magnetization
vector).
The direction of the magnetization (which is also the direction of n_up and n_dn) is given by [[spinat]].
_Warning: unlike collinear case, atoms having the same magnetization magnitude
with different directions must be given the same occupation matrix;
the magnetization will be oriented by the value of [[spinat]] (this is the
case for antiferro-magnetism). _

  * **Non-collinear magnetic system with zero magnetization** ([[nsppol]] = 1, [[nspinor]] = 2, [[nspden]] = 1):
Two (2lpawu+1)x(2lpawu+1) [[dmatpawu]] matrices are given for each atom on
which +U is applied.
They contain the "spin-up" and "spin-down" occupations;
But, as "spin-up" and "spin-down" are constrained identical, the "spin-down" one is ignored by the code.
""",
),

Variable(
    abivarname="dmatpuopt",
    varset="paw",
    vartype="integer",
    topics=['DFT+U_expert'],
    dimensions="scalar",
    defaultval=2,
    mnemonics="Density MATrix for PAW+U OPTion",
    requires="[[usepaw]] == 1 and [[usepawu]] == 1",
    added_in_version="before_v9",
    text=r"""
This option governs the way occupations of localized atomic levels are computed:

  * [[dmatpuopt]] = 1: atomic occupations are projections on atomic orbitals (Eq. (6) of [[cite:Amadon2008a]]).

  * [[dmatpuopt]] = 2: atomic occupations are integrated values in PAW spheres of angular-momentum-decomposed charge densities (Eq. (7) of [[cite:Amadon2008a]]).

  * [[dmatpuopt]] = 3: only for tests

  * [[dmatpuopt]] = 4: Extrapolations of occupancies outside the PAW-sphere. This Definition gives normalized operator for occupation.

In the general case [[dmatpuopt]] = 2 is suitable. The use of [[dmatpuopt]] = 1 is
restricted to PAW datasets in which the first atomic wavefunction of the
correlated subspace is a normalized atomic eigenfunction.
""",
),

Variable(
    abivarname="dmatudiag",
    varset="paw",
    vartype="integer",
    topics=['DFT+U_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Density MATrix for paw+U, DIAGonalization",
    requires="[[usepaw]] == 1 and [[usepawu]] == 1 and [[nspden]] != 4",
    added_in_version="before_v9",
    text=r"""
Relevant only for Ground-State calculations.
This option can be used to diagonalize the occupation matrix Nocc_{m,m_prime}.
Relevant values are:

  * 0: deactivated.
  * 1: occupation matrix is diagonalized and printed in log file at each SCF cycle
    (eigenvectors are also given in the log file).
  * 2: for testing purpose.
""",
),

Variable(
    abivarname="dmft_charge_prec",
    varset="dmft",
    vartype="real",
    topics=['DMFT_expert'],
    dimensions="scalar",
    defaultval=1e-06,
    mnemonics="Dynamical Mean Field Theory: charge density precision",
    added_in_version="before_v9",
    text=r"""
Precision to achieve in determining the charge density in the computation of the fermi level.
Should be decreased to increase precision. However, for a large system, it can increase importantly computer time.
""",
),

Variable(
    abivarname="dmft_dc",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Dynamical Mean Field Theory: Double Counting",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""

Value of double counting used for DMFT. Only value 1 is currently activated.
It corresponds to the "Full Localized Limit" double counting.
""",
),

Variable(
    abivarname="dmft_entropy",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamical Mean Field Theory: ENTROPY",
    requires="[[usedmft]] == 1 and [[dmft_solv]] == 5",
    added_in_version="before_v9",
    text=r"""
If 1, enables the calculation of the entropy within the DMFT framework and so
allows one the calculation of the total energy (free energy). In the current
implementation, this is only possible with [[dmft_solv]] = 5 (Continuous Time
Quantum Monte Carlo). See also the input variable [[dmft_nlambda]].
""",
),

Variable(
    abivarname="dmft_kspectral_func",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamical Mean Field Theory: compute K-resolved SPECTRAL FUNCtion",
    characteristics=['[[DEVELOP]]'],
    added_in_version="9.0.0",
    text=r"""

When activated, in conjunction with [[iscf]] = -2 or -3, a calculation
of k-resolved spectral function (or density of state) is possible.
However, the calculation requires as input the self-energy computed in the real
axis using an external analytical continuation code.
The section 7 of the DFT+DMFT tutorial  details how to obtain this data
and related informations.
""",
),


Variable(
    abivarname="dmft_iter",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamical Mean Field Theory: number of ITERation",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Number of iterations for the DMFT inner loop.
""",
),

Variable(
    abivarname="dmft_mxsf",
    varset="dmft",
    vartype="real",
    topics=['DMFT_useful'],
    dimensions="scalar",
    defaultval=0.3,
    mnemonics="Dynamical Mean Field Theory: MiXing parameter for the SelF energy",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Mixing parameter for the simple mixing of the self-energy (0.3 is safe, but it can be increased most of the time to 0.6).
""",
),

Variable(
    abivarname="dmft_nlambda",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_expert'],
    dimensions="scalar",
    defaultval=6,
    mnemonics="Dynamical Mean Field Theory: Number of LAMBDA points",
    characteristics=['[[DEVELOP]]'],
    requires="[[usedmft]] == 1 and [[dmft_entropy]] == 1",
    added_in_version="before_v9",
    text=r"""
[[dmft_nlambda]] gives the number of integration points for the
thermodynamic integration in case of free energy calculation within DMFT.
Its value must be greater or equal to 3.
""",
),

Variable(
    abivarname="dmft_nwli",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamical Mean Field Theory: Number of frequency omega (W) in the LInear mesh",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Number of Matsubara frequencies (linear mesh)
""",
),

Variable(
    abivarname="dmft_nwlo",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamical Mean Field Theory: Number of frequency omega (W) in the LOg mesh",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Number of frequencies in the log mesh.
""",
),

Variable(
    abivarname="dmft_occnd_imag",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Dynamical Mean Field Theory: Occupation non-diagonal imaginary part",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
When 0 force non-diagonal occupations imaginary parts to be null. Do not use this, it is only for compatibility with old tests.
""",
),

Variable(
    abivarname="dmft_read_occnd",
    varset="dev",
    vartype="integer",
    topics=['DMFT_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamical Mean Field Theory: READ OCCupations (Non Diagonal)",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Flag to read/write Occupations as computed in DMFT. This flag is useful to
restart a DFT+DMFT calculation with self-consistency over electronic density.
The occupations are written each time a DMFT loop is finished. So if the
calculation stops because the time limit is reached, this option offers the
possibility to restart the self-consistent loop over density at the point
where it stopped (assuming a restart with the wave functions, see [[getwfk]]).

  * 0 --> Occupations are written but never read.
  * 1 --> Occupations are read from I_DMFTOCCND, where I is the root for input files.
  * 2 --> Occupations are read from O_DMFTOCCND, where O is the root for output files.

An alternative and more simple way to restart a DFT+DMFT calculation is to use
the density file (obtained with [[prtden]] = 1 or [[prtden]] = -1) and the self-energy (see [[dmft_rslf]]).
In this case, use [[dmft_read_occnd]]=0.
""",
),

Variable(
    abivarname="dmft_rslf",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamical Mean Field Theory: Read SeLF energy",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Flag to read/write Self-Energy. If put to one, self-energy is written and read at each DFT iteration.
If self-energy file is missing, the self-energy is initialized to the double counting at the first iteration.
Importantly, in order to the calculation to restart easily, the self-energy is read and write in the same file.
""",
),

Variable(
    abivarname="dmft_solv",
    varset="dmft",
    vartype="real",
    topics=['DMFT_basic'],
    dimensions="scalar",
    defaultval=5,
    mnemonics="Dynamical Mean Field Theory: choice of SOLVer",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Choice of solver for the Impurity model.

  * 0 --> No solver and U=0, J=0 (see [[upawu]] and [[jpawu]]).
  * 1 --> LDA+U self-energy is used (for testing purpose)
  * 2 --> Hubbard one solver in the density density approximation of the Coulomb interaction. The Hubbard one solver is an approximation which gives a rough description of correlated Mott insulators. It should not be used for metals.
  * 5 --> Use the Continuous Time Quantum Monte Carlo (CTQMC) solver CT-Hyb of ABINIT in the density density approximation of the Coulomb interaction. The calculation is fully parallelised over MPI processes.
  * 6 --> Continuous Time Quantum Monte Carlo (CTQMC) solver CT-Hyb of TRIQS in the density density representation.
  * 7 --> Continuous Time Quantum Monte Carlo (CTQMC) solver CT-Hyb of TRIQS with the rotationally invariant formulation.
  * 8 --> Same as 5, but off-diagonal elements of the hybridization function are taken into account (useful for low symetry systems or with spin orbit coupling).
  * 9 --> Python invocation. Give a symbolic link to your python interpreter as an input like 'input-tag'_TRIQS_python_lib and the python script as an input like 'input-tag'_TRIQS_script.py. The inputs for the script will be written in dft_for_triqs.nc and the output as triqs_for_dft.nc.

The CT Hyb algorithm is described in [[cite:Werner2006]]. For a
discussion of density-density approximation with respect with the
rotationnally invariant formulation, see e.g. [[cite:Antipov2012]].
The ABINIT/CT Hyb implementation is discussed in [[cite:Gonze2016]].
The TRIQS/CT Hyb implementation is described in [[cite:Seth2016]].
Before using it, it has to be installed following instructions available [here](https://triqs.github.io/triqs/2.1.x).
Until release 8.10 included, the
interface was valid only for TRIQS 1.4 and TRIQS/CTHYB 1.4. It has then been upgraded to TRIQS 2.1 afterwards.
An example of a config.ac file to compile ABINIT with TRIQS can be found in [[ac:higgs_gnu_7.3_triqs2.ac]].
See the useful variables for CT-QMC solver: [[dmftctqmc_basis]],
[[dmftctqmc_check]], [[dmftctqmc_correl]], [[dmftctqmc_gmove]],
[[dmftctqmc_grnns]], [[dmftctqmc_meas]], [[dmftctqmc_mrka]],
[[dmftctqmc_mov]], [[dmftctqmc_order]], [[dmftctqmc_triqs_nleg]],
[[dmftqmc_l]], [[dmftqmc_n]], [[dmftqmc_seed]], [[dmftqmc_therm]]
""",
),

Variable(
    abivarname="dmft_t2g",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamical Mean Field Theory: t2g orbitals",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""

Can be set to 1 only if in cubic symmetry. It enables one to carry a DFT+DMFT
calculations only on _t<sub>2g</sub>_ orbitals.
""",
),

Variable(
    abivarname="dmft_tolfreq",
    varset="dmft",
    vartype="real",
    topics=['DMFT_expert'],
    dimensions="scalar",
    defaultval=0.0001,
    mnemonics="Dynamical Mean Field Theory: TOLerance on DFT correlated electron occupation matrix for the definition of the FREQuency grid",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""

The LDA occupation matrix for correlated electrons can be computed directly.
It can be compared to the calculation of the same quantity using LDA Green's
function, a sum over Matsubara frequencies and a projection over correlated
orbitals. Because the Matsubara grid is finite, the two quantities differ. If
this difference is larger than dmft_tolfreq, then the code stops and an error
message is given.
""",
),

Variable(
    abivarname="dmft_tollc",
    varset="dmft",
    vartype="real",
    topics=['DMFT_useful'],
    dimensions="scalar",
    defaultval=1e-05,
    mnemonics="Dynamical Mean Field Theory: TOLerance on Local Charge for convergence of the DMFT loop",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Tolerance for the variation of Local Charge for convergence of the DMFT Loop.
Most of the time however, DFT+DMFT calculations can converge fastly using [[dmft_iter]]=1, so
that this variable is not required.
""",
),

Variable(
    abivarname="dmftbandf",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamical Mean Field Theory: BAND: Final",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
[[dmftbandf]] is the last band taken into account in the Projected Local
Orbitals scheme of DFT+DMFT. With [[dmftbandi]], they define the energy window
used to define Wannier Functions (see [[cite:Amadon2008]]).
""",
),

Variable(
    abivarname="dmftbandi",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamical Mean Field Theory: BAND: Initial",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
[[dmftbandi]] is the first band taken into account in the Projected Local
Orbitals scheme of LDA+DMFT. With [[dmftbandf]], they define the energy window
used to define Wannier Functions (see [[cite:Amadon2008]]).
""",
),

Variable(
    abivarname="dmftcheck",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamical Mean Field Theory: CHECKs",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Only for developer purposes.
""",
),

Variable(
    abivarname="dmftctqmc_basis",
    varset="dev",
    vartype="integer",
    topics=['DMFT_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Dynamical Mean Field Theory: Continuous Time Quantum Monte Carlo BASIS",
    characteristics=['[[DEVELOP]]'],
    requires="[[dmft_solv]] == 5",
    added_in_version="before_v9",
    text=r"""
Choose the basis to perform CTQMC calculation.

  * 0 --> Use the local basis in the spherical harmonics basis.
  Can be useful if the Hamiltonian has weak off diagonal terms and for this reason,
  one want to keep the original basis for simplicity and for physical insight.
  * 1 --> Default value, diagonalize the local Hamiltonian (but only if it is not diagonal).
  The best choice in general.
  * 2 --> Diagonalise the local correlated occupation matrix. Can lead to non
  diagonal Hamiltonian that cannot be handled by CTQMC. This option should be thus avoided.
""",
),

Variable(
    abivarname="dmftctqmc_check",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamical Mean Field Theory: Continuous Time Quantum Monte Carlo CHECK",
    characteristics=['[[DEVELOP]]'],
    requires="[[dmft_solv]] == 5",
    added_in_version="before_v9",
    text=r"""
Check the fast calculations during the Monte Carlo simulation with very slow
but robust methods. Should only be used for debugging.

  * 0 --> No check.
  * 1 --> Check the overlap calculations (Impurity operator).
  * 2 --> Check the update of M matrix calculation (Bath operator).
  * 3 --> Check both.
""",
),

Variable(
    abivarname="dmftctqmc_correl",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamical Mean Field Theory: Continuous Time Quantum Monte Carlo CORRELations",
    characteristics=['[[DEVELOP]]'],
    requires="[[dmft_solv]] == 5",
    added_in_version="before_v9",
    text=r"""
Flag to compute statistics about segments and anti-segments during the
simulation. Slow down the simulation.

  * 0 --> Nothing done
  * 1 --> Calculations performed and written in "Correlation.dat" file
""",
),

Variable(
    abivarname="dmftctqmc_gmove",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamical Mean Field Theory: Continuous Time Quantum Monte Carlo Global MOVEs",
    characteristics=['[[DEVELOP]]'],
    requires="[[dmft_solv]] == 5",
    added_in_version="before_v9",
    text=r"""
Default is no global moves. The value of this variable is the modulo used to
try a global move. A value of 5000 means that a global move is tried every 5000 Monte Carlo sweep.
""",
),

Variable(
    abivarname="dmftctqmc_grnns",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamical Mean Field Theory: Continuous Time Quantum Monte Carlo GReeNs NoiSe",
    characteristics=['[[DEVELOP]]'],
    requires="[[dmft_solv]] == 5",
    added_in_version="before_v9",
    text=r"""
Compute the statistical noise for each time slice of each green function. This
is a good approximation only if there is enough Monte Carlo sweeps per cpu.

  * 0 --> Nothing
  * 1 --> Do it and write the noise in the "Gtau.dat" file.
""",
),

Variable(
    abivarname="dmftctqmc_meas",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Dynamical Mean Field Theory: Continuous Time Quantum Monte Carlo MEASurements",
    characteristics=['[[DEVELOP]]'],
    requires="[[dmft_solv]] == 5",
    added_in_version="before_v9",
    text=r"""
The modulo used to measure the interaction energy and the number of electrons.
Example: 2 means the measure is perform every two sweeps.
""",
),

Variable(
    abivarname="dmftctqmc_mov",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamical Mean Field Theory: Continuous Time Quantum Monte Carlo MOVie",
    characteristics=['[[DEVELOP]]'],
    requires="[[dmft_solv]] == 5",
    added_in_version="before_v9",
    text=r"""
Print a latex file per cpu displaying the full simulation. This option should
only be use with very small number (<1000) of Monte Carlo sweeps since it
requires a lot of I/O band width.

  * 0 --> Nothing
  * 1 --> Write the "Movie_id.dat" file where id is the MPI rank of each process
""",
),

Variable(
    abivarname="dmftctqmc_mrka",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamical Mean Field Theory: Continuous Time Quantum Monte Carlo MARKov Analysis",
    characteristics=['[[DEVELOP]]'],
    requires="[[dmft_solv]] == 5",
    added_in_version="before_v9",
    text=r"""
Measure the time evolution of the number of electrons for each orbital and
perform a fourier transform. The result can be plotted using the "Markov_id.dat" file

  * 0 --> Nothing
  * 1 --> Do it and write the noise in the "Markov_id.dat" file where id is the rank of each MPI process.
""",
),

Variable(
    abivarname="dmftctqmc_order",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamical Mean Field Theory: Continuous Time Quantum Monte Carlo perturbation ORDER",
    characteristics=['[[DEVELOP]]'],
    requires="[[dmft_solv]] == 5",
    added_in_version="before_v9",
    text=r"""
Print a file containing the statistic distribution of the number of segments
per orbital. The maximal order taken into account [[dmftctqmc_order]]: 50
means that we have the statistic distribution from 0 to 50 segments. The
result is written in the "Perturbation.dat" file.
""",
),

Variable(
    abivarname="dmftctqmc_triqs_nleg",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_expert'],
    dimensions="scalar",
    defaultval=30,
    mnemonics="Dynamical Mean Field Theory: Continuous Time Quantum Monte Carlo perturbation of TRIQS, Number of LEGendre polynomials",
    characteristics=['[[DEVELOP]]'],
    requires="[[dmft_solv]] in [6, 7]",
    added_in_version="before_v9",
    text=r"""
Specify the number of Legendre polynomials used for the calculation of Green's
function in CTQMC code from the library TRIQS. Default is 30. The value of
coefficients are given in file whose name ending is
"Legendre_coefficient.dat" (see also [[cite:Boehnke2011]]).
""",
),

Variable(
    abivarname="dmftqmc_l",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamical Mean Field Theory: Quantum Monte Carlo time sLices",
    characteristics=['[[DEVELOP]]'],
    requires="[[dmft_solv]] >= 5",
    added_in_version="before_v9",
    text=r"""
Number of time slices used to represent the time green function. This value
should be carefully chosen according to Niquist frequency and the [[tsmear]] value.
""",
),

Variable(
    abivarname="dmftqmc_n",
    varset="dmft",
    vartype="real",
    topics=['DMFT_compulsory'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="Dynamical Mean Field Theory: Quantum Monte Carlo Number of sweeps",
    characteristics=['[[DEVELOP]]'],
    requires="[[dmft_solv]] >= 5",
    added_in_version="before_v9",
    text=r"""
Number of Monte Carlo sweeps. Should be at least 10<sup>6<\sup>.
""",
),

Variable(
    abivarname="dmftqmc_seed",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_useful'],
    dimensions="scalar",
    defaultval="[[jdtset]]",
    mnemonics="Dynamical Mean Field Theory: Quantum Monte Carlo SEED",
    characteristics=['[[DEVELOP]]'],
    requires="[[dmft_solv]] >= 5",
    added_in_version="before_v9",
    text=r"""
Seed to initialize the random number generator.
Should not be relevant except for testing purpose.
NOTE: If the CT-QMC ([[dmft_solv]] = 5) is used on many CPUs, each CPU
initializes its random number generator with dmftqmc_seed+rank where rank is
the rank of the cpu in the MPI communicator.
""",
),

Variable(
    abivarname="dmftqmc_therm",
    varset="dmft",
    vartype="integer",
    topics=['DMFT_compulsory'],
    dimensions="scalar",
    defaultval=1000,
    mnemonics="Dynamical Mean Field Theory: Quantum Monte Carlo THERMalization",
    characteristics=['[[DEVELOP]]'],
    requires="[[dmft_solv]] == 5",
    added_in_version="before_v9",
    text=r"""
Number of Monte Carlo sweeps for the thermalization
""",
),

Variable(
    abivarname="dosdeltae",
    varset="gstate",
    vartype="real",
    topics=['printing_prdos', 'ElecDOS_useful'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="DOS DELTA in Energy",
    characteristics=['[[ENERGY]]'],
    added_in_version="before_v9",
    text=r"""
Defines the linear grid resolution (energy increment) to be used for the
computation of the Density-Of-States, when [[prtdos]] is non-zero.
If [[dosdeltae]] is set to zero (the default value), the actual increment is
0.001 Ha if [[prtdos]] = 1, and the much smaller value 0.00005 Ha if
[[prtdos]] = 2. This different default value arises because the [[prtdos]] = 1
case, based on a smearing technique, gives a quite smooth DOS, while the DOS
from the tetrahedron method, [[prtdos]] = 2, is rapidly varying.
""",
),

Variable(
    abivarname="dtion",
    varset="rlx",
    vartype="real",
    topics=['PIMD_compulsory', 'MolecularDynamics_compulsory', "GeoOpt_compulsory"],
    dimensions="scalar",
    defaultval=100,
    mnemonics="Delta Time for IONs",
    added_in_version="before_v9",
    text=r"""
Used for controlling ion time steps. If [[ionmov]] is set to 1, 6, 7 and 15, then
molecular dynamics is  used to update atomic positions in response to forces.
The parameter [[dtion]] is a time step in atomic units of time. (One atomic
time unit is 2.418884e-17 seconds, which is the value of Planck's constant in
hartree*sec.) In this case the atomic masses, in amu (given in array " [[amu]]
"), are used in Newton's equation and the viscosity (for [[ionmov]] =1) and
number of time steps are provided to the code using input variables "[[vis]]"
and "[[ntime]]". The code actually converts from masses in amu to masses in
atomic units (in units of electron masses) but the user enters masses in
[[amu]]. (The conversion from amu to atomic units (electron masses) is
1822.88851 electron masses/amu.)

A typical good value for [[dtion]] is about 100. The user must try several
values for [[dtion]] in order to establish the stable and efficient choice for
the accompanying amu, atom types and positions, and [[vis]] (viscosity).
For quenched dynamics ([[ionmov]] = 7), a larger time step might be taken, for
example 200. No meaning for RF calculations.
It is also used in geometric relaxation calculation with the FIRE algorithm
([[ionmov]]=15), where the time is virtual. A small dtion should be set, for example 0.03.
""",
),

Variable(
    abivarname="dynimage",
    varset="rlx",
    vartype="integer",
    topics=['PIMD_expert', 'TransPath_expert'],
    dimensions=['[[nimage]]'],
    defaultval=MultipleValue(number=None, value=1),
    mnemonics="list of DYNamic IMAGEs",
    commentdefault="if [[imgmov]] in [2,5] (String Method, NEB), <b>dynimage(1)</b>=0 and <b>dynimage([[nimage]])</b>=0.",
    added_in_version="before_v9",
    text=r"""
This input variable is relevant when sets of images are activated (see
[[imgmov]]). Not all images might be required to evolve from one time step to
the other. Indeed, in the String Method or the Nudged Elastic Band, one might
impose that the extremal configurations of the string are fixed. In case
[[dynimage]](iimage)=0, the image with index "iimage" will be consider as
fixed. Thus, there is no need to compute forces and stresses for this image at
each time step. The purpose of defining extremal images is to make the input/output easier.

In order to save CPU time, the computation of properties of static images
([[dynimage]](iimage)=0) can be avoided: see [[istatimg]] keyword.
""",
),

Variable(
    abivarname="ecut",
    varset="basic",
    vartype="real",
    topics=['Planewaves_compulsory'],
    dimensions="scalar",
    mnemonics="Energy CUToff",
    characteristics=['[[ENERGY]]'],
    added_in_version="before_v9",
    text=r"""
Used to define the kinetic energy cutoff which controls the number of planewaves at given k point. The allowed
plane waves are those with kinetic energy lower than **ecut**, which translates to the following constraint
on the planewave vector $\vec{G}$ in reciprocal space

$$\frac{1}{2}(2\pi)^2 (\vec{k}+\vec{G})^2 < \text{ecut}.$$

All planewaves inside this **basis sphere** centered at k are included in the basis (except if [[dilatmx]] is defined).
The cutoff can be specified in Ha units (the default), Ry, eV or Kelvin, since **ecut** has the
[[ENERGY]] characteristics. (1 Ha = 27.2113845 eV)

This is the single parameter which can have an enormous effect on the quality
of a calculation; basically the larger **ecut** is, the better converged the
calculation is. For fixed geometry, the total energy MUST always decrease as
**ecut** is raised because of the variational nature of the problem.

_Usually one runs at least several calculations at various **ecut** to
investigate the convergence needed for reliable results._

For k-points whose coordinates are build from 0 or 1/2, the implementation of
time-reversal symmetry that links coefficients of the wavefunctions in
reciprocal space has been realized. See the input variable [[istwfk]]. If
activated (which corresponds to the default mode), this input variable
[[istwfk]] will allow to divide the number of plane wave (npw) treated
explicitly by a factor of two. Still, the final result should be identical
with the 'full' set of plane waves.

See the input variable [[ecutsm]], for the smoothing of the kinetic energy,
needed to optimize unit cell parameters.
""",
),

Variable(
    abivarname="ecuteps",
    varset="gw",
    vartype="real",
    topics=['Susceptibility_compulsory', 'RandStopPow_compulsory'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="Energy CUT-off for EPSilon (the dielectric matrix)",
    characteristics=['[[ENERGY]]'],
    requires="[[optdriver]] in [3, 4]",
    added_in_version="before_v9",
    text=r"""
[[ecuteps]] determines the cut-off energy of the planewave set used to
represent the independent-particle susceptibility $\chi^{0}_{KS}$, the
dielectric matrix $\epsilon$, and its inverse.
It is not worth to take [[ecuteps]] bigger than four times [[ecutwfn]], this
latter limit corresponding to the highest Fourier components of a wavefunction
convoluted with itself. Usually, even twice the value of [[ecutwfn]] might
overkill. A value of [[ecuteps]] between 5 and 10 Hartree often (but not
always) leads to converged results (at the level of 0.01 eV for the energy
gap). In any case, a convergence study is worth.
""",
),

Variable(
    abivarname="ecutsigx",
    varset="gw",
    vartype="real",
    topics=['SelfEnergy_compulsory'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="Energy CUT-off for SIGma eXchange",
    characteristics=['[[ENERGY]]'],
    requires="[[optdriver]] == 4",
    added_in_version="before_v9",
    text=r"""
[[ecutsigx]] determines the cut-off energy of the planewave set used to
generate the exchange part of the self-energy operator. For norm-conserving
calculations, it is pointless to have [[ecutsigx]] bigger than 4*[[ecut]],
while for PAW calculations, the maximal useful value is [[pawecutdg]]. Thus,
if you do not care about CPU time, please use these values. If you want to
spare some CPU time, you might try to use a value between [[ecut]] and these upper limits.
""",
),

Variable(
    abivarname="ecutsm",
    varset="rlx",
    vartype="real",
    topics=['Planewaves_basic'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="Energy CUToff SMearing",
    characteristics=['[[ENERGY]]'],
    added_in_version="before_v9",
    text=r"""
This input variable is important when performing relaxation of unit cell size
and shape (non-zero [[optcell]]). Using a non-zero [[ecutsm]], the total
energy curves as a function of [[ecut]], or [[acell]], can be smoothed,
keeping consistency with the stress (and automatically including the Pulay
stress). The recommended value is 0.5 Ha. Actually, when [[optcell]]/=0,
ABINIT requires [[ecutsm]] to be larger than zero. If you want to optimize
cell shape and size without smoothing the total energy curve (a dangerous
thing to do), use a very small [[ecutsm]], on the order of one microHartree.

Technical information:
See Appendix B of [[cite:Laflamme2016]].
[[ecutsm]] allows one to define an effective kinetic energy for plane waves, close
to, but lower than, the maximal kinetic energy [[ecut]]. For kinetic energies
less than [[ecut]]-[[ecutsm]], nothing is modified, while between
[[ecut]]-[[ecutsm]] and [[ecut]], the kinetic energy is multiplied by
$1.0 / ( x^2(3+x-6x^2+3x^3 ))$,
where x = ([[ecut]] - kinetic_energy)/[[ecutsm]].
Note that $x^2(3+x-6x^2+3x^3)$ is 0 at x=0, with vanishing derivative,
and that at x=1, it is 1, with also vanishing derivative.
If [[ecutsm]] is zero, the unmodified kinetic energy is used.
[[ecutsm]] can be specified in Ha (the default), Ry, eV or Kelvin, since
[[ecutsm]] has the [[ENERGY]] characteristics. (1 Ha = 27.2113845 eV).
A few test for Silicon (diamond structure, 2 k-points) have shown 0.5 Ha to be
largely enough for [[ecut]] between 2Ha and 6Ha, to get smooth curves. It is
likely that this value is OK as soon as [[ecut]] is larger than 4Ha.
""",
),

Variable(
    abivarname="ecutwfn",
    varset="gw",
    vartype="real",
    topics=['Susceptibility_compulsory', 'SelfEnergy_compulsory'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[optdriver]] in [3, 4]': '[[ecut]]', 'defaultval': 0.0}),
    mnemonics="Energy CUT-off for WaveFunctioNs",
    characteristics=['[[ENERGY]]'],
    requires=" [[optdriver]] in [3, 4]",
    added_in_version="before_v9",
    text=r"""
[[ecutwfn]] determines the cut-off energy of the planewave set used to
represent the wavefunctions in the formula that generates the independent-
particle susceptibility $\chi^{0}_{KS}$ (for [[optdriver]] = 3), or the self-
energy (for [[optdriver]] = 4).
Usually, [[ecutwfn]] is smaller than [[ecut]], so that the wavefunctions are
filtered, and some components are ignored. As a side effect, the wavefunctions
are no more normalized, and also, no more orthogonal. Also, the set of plane
waves can be much smaller for [[optdriver]] = 3, than for [[optdriver]] = 4,
although a convergence study is needed to choose correctly both values.

The size of this set of planewaves is [[npwwfn]].
""",
),

Variable(
    abivarname="effmass_free",
    varset="dev",
    vartype="real",
    topics=['Artificial_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="EFFective MASS for the FREE electron",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
This parameter allows one to change the free electron mass, with respect to its experimental value.
The electron mass is simply changed in the Schrodinger equation.
Only for testing purposes, of course.
""",
),

Variable(
    abivarname="efield",
    varset="ffield",
    vartype="real",
    topics=['Berry_basic'],
    dimensions=[3],
    defaultval=MultipleValue(number=3, value=0.0),
    mnemonics="Electric FIELD",
    requires="[[berryopt]] in [4, 6]",
    added_in_version="before_v9",
    text=r"""
In case [[berryopt]] = 4, a finite electric field calculation is performed. The
value of this electric field, and its direction is determined by [[efield]].
It must be given in atomic units (1 a.u. of electric field= 514220624373.482
V/m, see note below), in cartesian coordinates.

References for the calculation under electric field (based on multi k point Berry phase):

  * [[cite:Nunes1994]]: real-space version of the finite-field Hamiltonian
  * [[cite:Nunes2001]]: reciprocal-space version of the finite-field Hamiltonian (the one presently implemented), and extensive theoretical analysis
  * [[cite:Souza2002]]: implementation of the finite-field Hamiltonian (reciprocal-space version)
  * [[cite:Zwanziger2012]]: extension to PAW formalism

See also [[cite:Umari2003]].

!!! note

    The atomic unit of electric field strength is: $\frac{e_{Cb}}{4\pi\varepsilon_0a_0^2}$, where
    $e_{Cb}$ is the electronic charge in Coulomb (1.60217653$^{-19}$), $\varepsilon_0$ is the
    electric constant (8.854187817d-12 F/m), and $a_0$ is the Bohr radius in meter
    (0.5291772108$^{-10}$).
""",
),

Variable(
    abivarname="efmas",
    varset="dfpt",
    vartype="integer",
    topics=['EffectiveMass_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="EFfective MASs",
    added_in_version="before_v9",
    text=r"""
Turns on effective mass tensor calculations. Such calculations requires the
non-self-consistent calculation of derivatives with respect to k, in the same
dataset. It must therefore be used with [[rfelfd]] = 2 (or 1).

  * 0 --> no effective mass tensor calculation
  * 1 --> effective mass tensor calculation

!!! note
    At the present time, both norm-conserving (NC) and PAW calculations are
    supported. Also, for PAW calculations only, [[nspinor]] == 2 and
    [[pawspnorb]] == 1 (i.e. spin-orbit (SO) calculations) is supported. NC SO
    calculations are NOT currently supported. Also, for both NC and PAW,
    [[nspden]]/=1 and [[nsppol]]/=1 are NOT supported.
""",
),

Variable(
    abivarname="efmas_bands",
    varset="dfpt",
    vartype="integer",
    topics=['EffectiveMass_useful'],
    dimensions=[2, '[[nkpt]]'],
    defaultval="The full range of band available in the calculation for each k-point.",
    mnemonics="EFfective MASs, BANDS to be treated.",
    requires="[[efmas]] == 1",
    added_in_version="before_v9",
    text=r"""
This variable controls the range of bands for which the effective mass is to
be calculated. If a band is degenerate, all other bands of the degenerate
group will automatically be treated, even if they were not part of the user specified range.
""",
),

Variable(
    abivarname="efmas_calc_dirs",
    varset="dfpt",
    vartype="integer",
    topics=['EffectiveMass_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="EFfective MASs, CALCulate along DIRectionS",
    requires="[[efmas]] == 1",
    added_in_version="before_v9",
    text=r"""
Allows the user to calculate the scalar effective mass of all bands specified
by [[efmas_bands]] along specific directions in reciprocal space. This is
particularly useful when considering degenerate bands, which are usually
warped, and thus cannot have their dispersion (hessian) and effective mass
expressed as a tensor. This allows the user to see the more complex angular
behavior of effective masses in these cases, for instance.

When [[efmas_calc_dirs]] == 0, no directions are read from the input file (using
[[efmas_dirs]]) and the effective masses along the 3 cartesian directions are output by default.

When [[efmas_calc_dirs]] == 1, 2 or 3, [[efmas_n_dirs]] directions are read from
[[efmas_dirs]], assuming cartesian, reduced or angular ($\theta$,$\phi$)
coordinates, respectively. In the case [[efmas_calc_dirs]] == 3, 2 real values
per directions are read, whereas 3 real values are read in the two other cases.
""",
),

Variable(
    abivarname="efmas_deg",
    varset="dfpt",
    vartype="integer",
    topics=['EffectiveMass_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="EFfective MASs, activate DEGenerate formalism",
    requires="[[efmas]] > 0",
    added_in_version="before_v9",
    text=r"""
Activate (==1) or not (==0) the treatment of degenerate bands
(criterion [[efmas_deg_tol]] is used to determine whether bands are degenerate).
Also compute the transport equivalent effective mass (see [[cite:Mecholsky2014]]).

[[efmas]] = 0 should only be used for testing purposes.
""",
),

Variable(
    abivarname="efmas_deg_tol",
    varset="dfpt",
    vartype="real",
    topics=['EffectiveMass_useful'],
    dimensions="scalar",
    defaultval=1e-05,
    mnemonics="EFfective MASs, DEGeneracy TOLerance",
    requires="[[efmas_deg]] == 1",
    added_in_version="before_v9",
    text=r"""
Energy difference below which 2 bands are considered degenerate (and treated
using the formalism activated with [[efmas_deg]] == 1). [[efmas_deg_tol]] has
the [[ENERGY]] characteristics.
""",
),

Variable(
    abivarname="efmas_dim",
    varset="dfpt",
    vartype="integer",
    topics=['EffectiveMass_useful'],
    dimensions="scalar",
    defaultval=3,
    mnemonics="EFfective MASs, DIMension of the effective mass tensor",
    requires="[[efmas]] == 1",
    added_in_version="before_v9",
    text=r"""
For 2D or 1D systems, the band dispersion goes to 0 perpendicular to the
system, which causes the inverse effective mass to be singular, i.e. the
effective mass to be NaN. This keyword circumvents the problem by eliminating
the troublesome dimensions from the inverse effective mass.

In 2D, the Z axis is ignored and, in 1D, the Z and Y axis are ignored.

Also, note that in the 2D degenerate case, a subtlety arises: the 'transport
equivalent' effective mass does not determine the scale of the transport
tensors (conductivity and others). Therefore, for this specific case, the
factor by which these transport tensors should be scaled once determined from
the 'transport equivalent' effective mass tensor is output separately on the
line immediately after the effective mass.
""",
),

Variable(
    abivarname="efmas_dirs",
    varset="dfpt",
    vartype="real",
    topics=['EffectiveMass_basic'],
    dimensions=['3 or 2', '[[efmas_n_dirs]]'],
    defaultval=0,
    mnemonics="EFfective MASs, DIRectionS to be calculated",
    requires="[[efmas_calc_dirs]] > 0",
    added_in_version="before_v9",
    text=r"""
List of [[efmas_n_dirs]] directions to be considered according to the value of
[[efmas_calc_dirs]]. The directions are specified by 3 real values if
[[efmas_calc_dirs]] == 1 or 2 and by 2 real values if [[efmas_calc_dirs]] == 3.
""",
),

Variable(
    abivarname="efmas_n_dirs",
    varset="dfpt",
    vartype="integer",
    topics=['EffectiveMass_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="EFfective MASs, Number of DIRectionS",
    requires="[[efmas_calc_dirs]] > 0",
    added_in_version="before_v9",
    text=r"""
Number of directions in [[efmas_dirs]], to be considered according to [[efmas_calc_dirs]].
""",
),

Variable(
    abivarname="efmas_ntheta",
    varset="dfpt",
    vartype="integer",
    topics=['EffectiveMass_basic'],
    dimensions="scalar",
    defaultval=1000,
    mnemonics="EFfective MASs, Number of points for integration w/r to THETA",
    requires="[[efmas]] == 1 and [[efmas_bands]] == (degenerate band index)",
    added_in_version="before_v9",
    text=r"""
When a band is degenerate, the usual definition of effective mass becomes
invalid. However, it is still possible to define a 'transport equivalent mass
tensor' that reproduces the contribution of the band to the conductivity
tensor. To obtain this tensor, an integration over the solid sphere is
required. The angular variables are sampled using [[efmas_ntheta]] points for the theta coordinate,
and twice [[efmas_ntheta]] points for the phi coordinate.
The default value gives a tensor accurate to the 4th decimal in Ge.
""",
),

Variable(
    abivarname="einterp",
    varset="basic",
    vartype="real",
    topics=['ElecBandStructure_useful', 'SelfEnergy_expert'],
    dimensions=[4],
    defaultval=[0, 0, 0, 0],
    mnemonics="Electron bands INTERPolation",
    added_in_version="before_v9",
    text=r"""
This variable activates the interpolation of the electronic eigenvalues. It
can be used to interpolate KS eigenvalues at the end of the GS run or to
interpolate GW energies in sigma calculations ([[optdriver]] = 4). The k-path
can be specified with [[kptbounds]] and [[nkpath]].
einterp consists of 4 entries.
The first element specifies the interpolation method.

  * 0 --> No interpolation (default)
  * 1 --> Star-function interpolation (Shankland-Koelling-Wood Fourier interpolation scheme, see [[cite:Pickett1988]]

The meaning of the other entries depend on the interpolation technique selected.
In the case of star-function interpolation:

  * einterp(2): Number of star-functions per ab-initio k-point
  * einterp(3): If non-zero, activate Fourier filtering according to Eq 9 of [[cite:Uehara2000]].
    In this case, rcut is given by einterp(2) * Rmax where Rmax is the maximum length of
    the lattice vectors included in the star expansion
  * einterp(4): Used if einterp(2) /= 0. It defines rsigma in Eq 9

""",
),

Variable(
    abivarname="elph2_imagden",
    varset="dfpt",
    vartype="real",
    topics=['TDepES_basic'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="ELectron-PHonon interaction at 2nd order: IMAGinary shift of the DENominator",
    characteristics=['[[ENERGY]]'],
    requires="[[ieig2rf]] != 0",
    added_in_version="before_v9",
    text=r"""

The variable [[elph2_imagden]] determines the imaginary shift of the
denominator of the sum-over-states in the perturbation,
$(e_{nk}-e_{n'k'}+i$[[elph2_imagden]]). One should use a width comparable with
the Debye frequency or the maximum phonon frequency.
Can be specified in Ha (the default), Ry, eV or Kelvin, since [[ecut]] has the
[[ENERGY]] characteristics (1 Ha = 27.2113845 eV).
""",
),

Variable(
    abivarname="enunit",
    varset="gstate",
    vartype="integer",
    topics=['Output_useful','ElecBandStructure_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="ENergy UNITs",
    added_in_version="before_v9",
    text=r"""
Governs the units to be used for output of eigenvalues (and eventual phonon frequencies)

  * 0  --> print eigenvalues in hartree;
  * 1  --> print eigenvalues in eV;
  * 2  --> print eigenvalues in both hartree and eV.

If phonon frequencies are to be computed:

  * 0  -->  phonon frequencies in Hartree and cm-1;
  * 1  -->  phonon frequencies in eV and THz;
  * 2  -->  phonon frequencies in hartree, eV, cm-1, Thz and Kelvin.
""",
),

Variable(
    abivarname="eph_extrael",
    varset="eph",
    vartype="real",
    topics=['ElPhonInt_expert'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="Electron-PHonon: EXTRA ELectrons",
    added_in_version="before_v9",
    text=r"""
Number of electrons per unit cell to be added to the initial value computed
from the pseudopotentials and unit cell.
""",
),

Variable(
    abivarname="eph_fermie",
    varset="eph",
    vartype="real",
    topics=['ElPhonInt_expert'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="Electron-PHonon: FERMI Energy",
    characteristics=['[[ENERGY]]'],
    added_in_version="before_v9",
    text=r"""
This variable can be used to change the value of the Fermi level when
performing electron-phonon calculations with [[optdriver]] == 7. This variable
has effect only if set to a non-zero value. See also [[eph_extrael]].
""",
),

Variable(
    abivarname="eph_frohlichm",
    varset="eph",
    vartype="integer",
    topics=['ElPhonInt_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Electron-PHonon: FROHLICH Model",
    added_in_version="before_v9",
    text=r"""
Only relevant for [[optdriver]]=7 and [[eph_task]]=6.
If set to 1, use the dynamical matrix at Gamma, the Born effective charges, the dielectric tensor, as well as
the effective masses (must give a _EFMAS file as input, see [[prtefmas]] and [[getefmas]] or [[irdefmas]]),
as the parameters of a Frohlich Hamiltonian.
Then use these to compute the change of electronic eigenvalues due to electron-phonon interaction,
using second-order time-dependent perturbation theory. Can deliver (approximate) zero-point renormalisation
as well as temperature dependence.
""",
),

Variable(
    abivarname="eph_fsewin",
    varset="eph",
    vartype="real",
    topics=['ElPhonInt_expert'],
    dimensions="scalar",
    defaultval="0.01 Hartree",
    mnemonics="Electron-Phonon: Fermi Surface Energy WINdow",
    characteristics=['[[ENERGY]]'],
    added_in_version="before_v9",
    text=r"""
This variable defines the energy window around the Fermi level used for e-ph
calculations ([[optdriver]] = 7). Only states located in the energy range
[efermi - eph_fsewin, efermi + eph_fsewin] are included in the e-ph calculation.

Related input variables: [[eph_intmeth]], [[eph_fsmear]], [[eph_extrael]] and [[eph_fermie]].
""",
),

Variable(
    abivarname="eph_fsmear",
    varset="eph",
    vartype="real",
    topics=['ElPhonInt_expert'],
    dimensions="scalar",
    defaultval="0.01 Hartree",
    mnemonics="Electron-PHonon: Fermi surface SMEARing",
    characteristics=['[[ENERGY]]'],
    requires="[[eph_intmeth]] == 1",
    added_in_version="before_v9",
    text=r"""
This variable defines the gaussian broadening used for the integration over
the Fermi surface when [[eph_intmeth]] == 1.
""",
),

Variable(
    abivarname="eph_intmeth",
    varset="eph",
    vartype="integer",
    topics=['ElPhonInt_expert'],
    dimensions="scalar",
    defaultval="2 (tetra) except when [[eph_task]] = +4 where 1 is used as default.",
    mnemonics="Electron-Phonon: INTegration METHod",
    added_in_version="before_v9",
    text=r"""
This variable defines the technique for the integration over the Brillouin zone in the EPH code.

* 1 --> Gaussian technique with broadening factor
* 2 --> Tetrahedron method.

Note that the default value depends on the value of [[eph_task]] i.e. on the physical properties
we are computing.

Phonon linewidths in metals (**eph_task** = 1):

:   The default approach for the integration of the double-delta over the Fermi surface is 2 (tetrahedron).
    When the gaussian method is used, the broadening is given by [[eph_fsmear]].
    See also [[eph_fsewin]].

Electron-phonon self-energy

:   The default is gaussian method with broadening specified by [[zcut]].
""",
),

Variable(
    abivarname="eph_mustar",
    varset="eph",
    vartype="real",
    topics=['ElPhonInt_expert'],
    dimensions="scalar",
    defaultval=0.1,
    mnemonics="Electron-PHonon: MU STAR (electron-electron interaction strength)",
    added_in_version="before_v9",
    text=r"""
Average electron-electron interaction strength, for the computation of the
superconducting Tc using Mc-Millan's formula.
""",
),

Variable(
    abivarname="eph_ngqpt_fine",
    varset="eph",
    vartype="integer",
    topics=['ElPhonInt_useful'],
    dimensions=[3],
    defaultval=[0, 0, 0],
    mnemonics="Electron-PHonon: Number of Grid Q-PoinTs in FINE grid.",
    added_in_version="before_v9",
    text=r"""
This variable activates the **interpolation** of the first-order variation of the
self-consistent potential in the electron-phonon code ([[optdriver]] == 7).

If eph_nqgpt_fine differs from [0, 0, 0], the code will use the Fourier transform to interpolate
the DFPT potentials on this fine q-mesh starting from the irreducible set of
q-points read from the DVDB file. This approach is similar to the one used to
interpolate the interatomic force constants in q-space. If **eph_ngqpt_fine** is
not given, the EPH code uses the list of irreducible q-points reported in the
DDB file i.e. [[ddb_ngqpt]] (default behavior).

!!! important

    The computation of the e-ph matrix elements requires the knowledge of $\psi_{\bf k}$
    and $\psi_{\bf k + q}$. This means that the k-mesh for electrons found in the WFK must be
    compatible with the one given in *eph_ngqpt_fine*.
    The code can interpolate DFPT potentials but won't try to interpolate KS wavefunctions.
    and will stop if ${\bf k + q}$ is not found in the WFK file.
""",
),

Variable(
    abivarname="eph_task",
    varset="eph",
    vartype="integer",
    topics=['ElPhonInt_basic'],
    dimensions="scalar",
    defaultval=1,
    requires="[[optdriver]] == 7",
    mnemonics="Electron-PHonon: Task",
    added_in_version="before_v9",
    text=r"""
Select the electron-phonon task to be performed when [[optdriver]] == 7.
The choice is among:

* 0 --> No computation (mainly used to access the post-processing tools)
* 1 --> Compute phonon linewidths in metals and superconducting properties (isotropic formalism).
* 2 --> Compute e-ph matrix elements. Save results in GKK.nc file.
* -2 --> Compute e-ph matrix elements. Save results in GKQ.nc file that can be post-processed with AbiPy.
* 3 --> Compute phonon self-energy.
* 4 --> Compute electron self-energy (Fan-Migdal + Debye-Waller) and QP corrections. Generate SIGEPH.nc file.
* -4 --> Compute electron lifetimes due to e-ph interaction (imaginary part of Fan-Migdal self-energy). Generate SIGEPH.nc file.
* 5 --> Interpolate DFPT potentials to produce a new DVDB file on the [[eph_ngqpt_fine]] q-mesh that can be read with [[getdvdb]]
* -5 --> Interpolate DFPT potentials on the q-path specified by [[ph_qpath]] and [[ph_nqpath]]. Note that, in this case,
         the user has to provide the full list of q-points in the input, [[ph_ndivsm]] is not used to generate the q-path.
* 6 --> Estimate correction to the ZPR in polar materials using the Frohlich model. Requires EFMAS.nc file.
* 7 --> Compute phonon limited transport in semiconductors using lifetimes taken from SIGEPH.nc file.
* 15, -15 --> Write the average in r-space of the DFPT potentials to the V1QAVG.nc file.
              In the first case (+15) the q-points are specified via [[ph_nqpath]] and [[ph_qpath]]. The code assumes the
              input DVDB contains q-points in the IBZ and the potentials along the path are interpolated with Fourier transform.
              An array D(R) with the decay of the W(R,r) as a function of R is computed and saved to file
              In the second case (-15) the q-points are taken directly from the DVDB file.
""",
),

Variable(
    abivarname="eph_transport",
    varset="eph",
    vartype="integer",
    topics=['ElPhonInt_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Electron-PHonon: TRANSPORT flag",
    added_in_version="before_v9",
    text=r"""
NB - this does not work yet. This variable can be used to turn on the
calculation of transport quantities in the eph module of abinit.
Value of 1 corresponds to elastic LOVA as in [[cite:Savrasov1996]].
""",
),

Variable(
    abivarname="eshift",
    varset="dev",
    vartype="real",
    topics=['SCFAlgorithms_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Energy SHIFT",
    characteristics=['[[DEVELOP]]', '[[ENERGY]]'],
    requires="[[wfoptalg]] == 3",
    added_in_version="before_v9",
    text=r"""
[[eshift]] gives the shift of the energy used in the shifted Hamiltonian
squared. The algorithm will determine eigenvalues and eigenvectors centered on [[eshift]].
Can be specified in Ha (the default), Ry, eV or Kelvin, since [[eshift]] has the
[[ENERGY]] characteristics. (1 Ha = 27.2113845 eV)
""",
),

Variable(
    abivarname="esmear",
    varset="dfpt",
    vartype="real",
    topics=['TDepES_useful'],
    dimensions="scalar",
    defaultval=0.01,
    mnemonics="Eigenvalue SMEARing",
    characteristics=['[[ENERGY]]'],
    requires="[[smdelta]] != 0",
    added_in_version="before_v9",
    text=r"""
The variable [[esmear]] determines the width of the functions approximating
the delta function, $\delta(e_{nk}-e_{n'k'})$, present in the expression of the
lifetimes. One should use a width comparable with the Debye frequency or the
maximum phonon frequency.
Can be specified in Ha (the default), Ry, eV or Kelvin, since [[ecut]] has the
[[ENERGY]] characteristics (1 Ha = 27.2113845 eV).
""",
),

Variable(
    abivarname="exchmix",
    varset="dev",
    vartype="real",
    topics=['xc_useful'],
    dimensions="scalar",
    defaultval=0.25,
    mnemonics="EXCHange MIXing",
    characteristics=['[[DEVELOP]]'],
    requires="[[useexexch]] == 1",
    added_in_version="before_v9",
    text=r"""
[[exchmix]] allows one to tune the ratio of exact exchange when [[useexexch]] is
used. The default value of 0.25 corresponds to PBE0.
""",
),

Variable(
    abivarname="exchn2n3d",
    varset="dev",
    vartype="integer",
    topics=['TuningSpeed_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="EXCHange N2 and N3 Dimensions",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
If [[exchn2n3d]] is 1, the internal representation of the FFT arrays in
reciprocal space will be array(n1,n3,n2), where the second and third
dimensions have been switched. This is to allow to be coherent with the
[[exchn2n3d]] = 4xx FFT treatment.
""",
),

Variable(
    abivarname="extrapwf",
    varset="dev",
    vartype="integer",
    topics=['TuningSpeed_expert', 'MolecularDynamics_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="flag - EXTRAPolation of the Wave-Functions",
    characteristics=['[[DEVELOP]]'],
    requires="[[densfor_pred]] in [5, 6]",
    added_in_version="before_v9",
    text=r"""
This flag activates the extrapolation of wave-functions from one Molecular
Dynamics (or Structural Relaxation) step to another. The wave functions are
extrapolated using 2nd-order algorithm of [[cite:Arias1992]].
Note that, when activated, this extrapolation requires non-negligible
additional memory resources as the wave functions are stored for the two
previous time steps. Also, it can only be activated if a consistent density
extrapolation is activated (see [[densfor_pred]]).
ABINIT 7.10: this option is **under development** and might give wrong results.
""",
),

Variable(
    abivarname="f4of2_sla",
    varset="paw",
    vartype="real",
    topics=['DFT+U_expert'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'d electrons': 0.625,
 'f electrons': 0.6681,
 'defaultval': 0}),
    mnemonics="F4 Over F2 ratio of Slater integrals",
    requires="[[usepaw]] == 1 and ([[usepawu]] == 1 or [[usedmft]] == 1)",
    added_in_version="before_v9",
    text=r"""
This gives the ratio of Slater Integrals F4 and F2. It is used in DFT+U or
DFT+DMFT for the calculation of the orbital dependent screened coulomb interaction.
""",
),

Variable(
    abivarname="f6of2_sla",
    varset="paw",
    vartype="real",
    topics=['DFT+U_expert'],
    dimensions="scalar",
    defaultval=0.4943,
    mnemonics="F6 Over F2 ratio of Slater integrals",
    requires="([[usepawu]] == 1 or [[usedmft]] == 1) and [[lpawu]] == 3",
    added_in_version="before_v9",
    text=r"""
Gives the ratio of Slater Integrals F6 and F2. It is used with
[[f4of2_sla]] == 3 in DFT+U or DFT+DMFT for the calculation of the orbital
dependent screened coulomb interaction.
""",
),

Variable(
    abivarname="fband",
    varset="gstate",
    vartype="real",
    topics=['BandOcc_useful'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[occopt]] == 1': 0.125,
 '[[occopt]] > 2': 0.5,
 '[[usewvl]] == 1': 0.0,
 'defaultval': 0.0}),
    mnemonics="Factor for the number of BANDs",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Governs the number of bands to be used in the code in the case the parameter
[[nband]] is not defined in the input file (which means that [[occopt]] is not
equal to 0 or 2).

In case [[fband]] is 0.0, the code computes from the pseudopotential files
and the geometry data contained in the input file, the number of electrons
present in the system. Then, it computes the minimum number of bands that can
accommodate them, and use that value for [[nband]].

In case [[fband]] differs from zero, other bands will be added, just larger
than [[fband]] times the number of atoms. This parameter is not echoed in the
top of the main output file, but only the parameter [[nband]] that it allowed
to compute. It is also not present in the dtset array (no internal).
The default values are chosen such as to give naturally some conduction bands.
This improves the robustness of the code, since this allows one to identify lack
of convergence coming from (near-)degeneracies at the Fermi level. In the
metallic case, the number of bands generated might be too small if the
smearing factor is large. The occupation numbers of the higher bands should be
small enough such as to neglect higher bands. It is difficult to automate
this, so a fixed default value has been chosen.
""",
),

Variable(
    abivarname="fermie_nest",
    varset="dev",
    vartype="real",
    topics=['printing_prfermi'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="FERMI Energy for printing the NESTing function",
    added_in_version="before_v9",
    text=r"""
This input variable is only effective when [[prtnest]] = 1.
The energy is relative to the calculated Fermi energy.
""",
),

Variable(
    abivarname="fftalg",
    varset="dev",
    vartype="integer",
    topics=['TuningSpeed_useful'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[FFTW3]] and [[usedmft]] == 0': 312,
 '[[paral_kgb]] == 1': 401,
 'defaultval': 112}),
    mnemonics="Fast Fourier Transform ALGorithm",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
This keyword is **irrelevant** when Fast Fourier Transforms are done using
**Graphics Processing Units** (GPU), i.e. when [[use_gpu_cuda]] = 1
(in that case, it is ignored).

Allows to choose the algorithm for Fast Fourier Transforms. These have to be
used when applied to wavefunctions (routine `fourwf`), as well as when
applied to densities and potentials (routine `fourdp.F90`). Presently, it is the
concatenation of three digits, labelled (A), (B) and (C).

The first digit (A) is to be chosen among 1, 2, 3, 4 or 5:

  * 1 --> use FFT routines written by S. Goedecker.
  * 2 -->  not available anymore
  * 3 -->  use serial or multi-threaded FFTW3 Fortran routines ([http://www.fftw.org ](http://www.fftw.org)).
    Currently implemented with [[fftalg]] = 312.
  * 4 -->  use FFT routines written by S. Goedecker, 2002 version, that will be suited for MPI and OpenMP parallelism.
  * 5 -->  use serial or multi-threaded MKL routines Currently implemented with [[fftalg]] = 512.

The second digit (B) is related to `fourdp`:

  * 0 -->  only use Complex-to-complex FFT
  * 1 -->  real-to-complex is also allowed (only coded for A==1, A==3 and A==5)

The third digit (C) is related to `fourwf`:

  * 0 --> no use of zero padding
  * 1 --> use of zero padding (only coded for A==1, A==4)
  * 2 --> use of zero padding, and also combines actual FFT operations (using 2 routines from S. Goedecker)
    with important pre- and post-processing operations, in order to maximize cache data reuse.
    This is very efficient for cache architectures. (coded for A==1 and A==4, but A==4 is not yet sufficiently tested)

Internal representation as [[ngfft]](7).
""",
),

Variable(
    abivarname="fftcache",
    varset="dev",
    vartype="integer",
    topics=['TuningSpeed_expert'],
    dimensions="scalar",
    defaultval=16,
    mnemonics="Fast Fourier Transform CACHE size",
    characteristics=['[[DEVELOP]]'],
    commentdefault="todo: Not yet machine-dependent",
    added_in_version="before_v9",
    text=r"""
Gives the cache size of the current machine, in Kbytes.
Internal representation as [[ngfft]](8).
""",
),

Variable(
    abivarname="fftgw",
    varset="gw",
    vartype="integer",
    topics=['GW_expert', 'Susceptibility_expert', 'SelfEnergy_expert'],
    dimensions="scalar",
    defaultval=21,
    mnemonics="FFT for GW calculation",
    requires=" [[optdriver]] in [3, 4]",
    added_in_version="before_v9",
    text=r"""
The basic ingredients needed to perform both a screening and a sigma
calculation are the so-called oscillator matrix elements defined as

$$ \langle \mathbf{k-q},b_1 | e^{-i (\mathbf{q+G)} \mathbf{r}} | \mathbf{k}, b_2 \rangle $$

In reciprocal space, this expression is evaluated by a convolution in which
the number of reciprocal lattice vectors employed to describe the
wavefunctions is given by [[ecutwfn]]. In the case of screening calculations,
the number of **G** vectors in the above expression is defined by [[ecuteps]],
while [[ecutsigx]] defined the number of **G** used in sigma calculations. To
improve the efficiency of the code, the oscillator matrix elements are
evaluated in real space through FFT techniques, and the [[fftgw]] input
variable is used to select the FFT mesh to be used.

[[fftgw]] is the concatenation of two digits, labelled (A) and (B) whose value
is internally used to define the value of [[ngfft]](1:3) (see the setmesh.F90 routine).

The first digit (A) defines the augmentation of the FFT grid. Possible values are 1, 2 and 3.

  * 0 --> Use the FFT grid specified by the user through [[ngfft]](1:3)
  * 1 --> Use a coarse FFT grid which encloses a sphere in reciprocal space whose radius
    depends on the largest value between [[ecutwfn]] and [[ecuteps]]
  * 2 --> Use a slightly augmented FFT which is sufficient for the correct treatment of the convolution
  * 3 --> Doubled FFT grid (same mesh as that used for GS calculations).

The second digit (B) can be chosen between 0 and 1. It defines whether a FFT
grid compatible with all the symmetries of the space group must be enforced or not:

  * 0 --> Use the smallest FFT mesh which is compatible with the FFT library (faster, save memory but is less accurate)
  * 1 --> Enforce a FFT grid which is compatible with all the symmetry operations of the space group.
    This method leads to an increase both of CPU time and memory, but the matrix elements are more accurate
    since the symmetry properties of the system are preserved.

The behaviour of ABINIT before v5.5 corresponds to the default value 11.
""",
),

Variable(
    abivarname="fockdownsampling",
    varset="gstate",
    vartype="integer",
    topics=['Hybrids_useful'],
    dimensions=[3],
    defaultval="3*1",
    mnemonics="FOCK operator, k-grid DOWNSAMPLING",
    requires="[[usefock]] == 1",
    added_in_version="before_v9",
    text=r"""
The k wavevector grid used to compute the Fock operator in the full Brillouin
Zone might be a subset of the full Brillouin Zone of the k point grid used for
the wavefunctions (see [[kptopt]] for its definition). [[fockdownsampling]]
allows one to define such a reduced set of k wavevectors, yielding a
significant speed up at the expense (possibly) of accuracy. The final grid has
[[nkpthf]] k points, with coordinates [[kptns_hf]].

In the simplest case, the three values of [[fockdownsampling]] are equal, and
simply gives the factor by which the initial k point grid will be downsampled
along every direction. If [[fockdownsampling]] is 3*1, the sampling for the
Fock operator is the same as the sampling for the wavefunctions. If
[[fockdownsampling]] is 3*2 the sampling will have 8 times less k points.
Conventionally, if [[fockdownsampling]] is 3*0, then the Fock operator is
obtained solely from the Gamma point. Also, as soon as [[fockdownsampling]] is
not 3*1 or 3*0, the k point grid from which a subset will be taken is obtained
by imposing [[nshiftk]] = 1.

A more accurate description is now given, as one can achieve a better control
than described above, with differing values of [[fockdownsampling]] and also
with negative numbers. One starts from the k point grid defined by
[[kptrlatt]], with [[nshiftk]] = 1. The absolute value of each of the three
numbers of [[fockdownsampling]] is used to sample the corresponding axis (in
reduced coordinate), as described above. Moreover, the obtained k grid might
even be further downsampled by specifying negative numbers: if all three are
negative, one further downsamples each 2x2x2 subset of k points by taking its
FCC subset (so, 4 points instead of 8); if two are negative, one downsamples
each 2x2x2 subset of k points by taking its BCC subset (so, 2 points instead
of 8); is only one is negative, then the two other axes are sampled using a
face-centered sampling. Finally, some of the values might be taken as 0, in
which case the corresponding direction is sampled by only one layer of points
(if two are zero, a line of points is obtained).
""",
),

Variable(
    abivarname="fockoptmix",
    varset="gstate",
    vartype="integer",
    topics=['Hybrids_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="FOCK operator: OPTions for MIXing",
    requires="[[usefock]] == 1",
    added_in_version="before_v9",
    text=r"""
Governs the mixing algorithm at the level of the Fock operator, i.e. how to
mix it, and how the underlying SCF calculation is to be performed. It is the
most relevant when the Fock operator is not updated at each SCF step ([[nnsclohf]]/=0).

The last digit of [[fockoptmix]] governs what happens at the level of the SCF
algorithm, when the Fock operator is updated.

  1. If [[fockoptmix]] == 0: the SCF algorithm is not restarted
  (it continues to use the previous potential/density pairs without worrying).
  2. If [[fockoptmix]] == 1: the SCF algorithm is restarted (the previous potential/density pairs are discarded).

The second-to-last (dozen) digit governs the possible modification of the XC
functional inside the SCF loop to take into account the lack of update of the
Fock operator. Irrelevant when the unit digit is 0. If the value 1 is used
(so, e.g. [[fockoptmix]] == 11), an auxiliary xc functional is used inside the
SCF loop, for a frozen ACE Fock operator. This auxiliary functional is
specified thanks to [[auxc_ixc]] and [[auxc_scal]].

The third-to-last (hundreds) digit governs the mixing of the Fock operator
itself with its previous occurrences. Irrelevant when the unit digit is 0.
""",
),

Variable(
    abivarname="freqim_alpha",
    varset="gw",
    vartype="real",
    topics=['SelfEnergy_expert'],
    dimensions="scalar",
    defaultval=5.0,
    mnemonics="FREQuencies along the IMaginary axis ALPHA parameter",
    requires="[[optdriver]] == 4",
    added_in_version="before_v9",
    text=r"""
[[freqim_alpha]] is used only for numerical integration of the GW self-energy
([[gwcalctyp]] in  [2, 12, 22, 9, 19, 29]).
[[freqim_alpha]] determines the location of the maximum frequency point along
the imaginary axis if the default grid is used in Contour Deformation
(numerical integration) calculations. It is set as  $\alpha*\omega_p$, where $\omega_p$
is the plasma frequency determined by the average density of the system
(this can be set by hand by using the variable [[ppmfrq]]). See the section on
grids in the descriptive text for [[cd_frqim_method]] for a detailed
""",
),

Variable(
    abivarname="freqremax",
    varset="gw",
    vartype="real",
    topics=['FrequencyMeshMBPT_basic'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="FREQuencies along the Real axis MAXimum",
    characteristics=['[[ENERGY]]'],
    requires="[[optdriver]] == 3",
    added_in_version="before_v9",
    text=r"""
[[freqremax]] is used only for numerical integration of the GW self-energy
([[gwcalctyp]] in  [2, 12, 22, 9, 19, 29]).
[[freqremax]] sets the maximum real frequency used to calculate the dielectric
matrix in order to perform the numerical integration of the GW self-energy.
[[freqremax]], [[freqremin]] and [[nfreqre]] define the spacing of the frequency mesh along the real axis.
""",
),

Variable(
    abivarname="freqremin",
    varset="gw",
    vartype="real",
    topics=['FrequencyMeshMBPT_basic'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="FREQuencies along the Real axis MINimum",
    characteristics=['[[ENERGY]]'],
    requires="[[optdriver]] == 3",
    added_in_version="before_v9",
    text=r"""
[[freqremin]] is used only for numerical integration of the GW self-energy
([[gwcalctyp]] in  [2, 12, 22, 9, 19, 29]).
[[freqremin]] sets the minimum real frequency used to calculate the dielectric
matrix in order to perform the numerical integration of the GW self-energy.
[[freqremin]] can be used to split a wide frequency interval into smaller
subintervals that can be calculated independently. The different subintervals
can then be merged together with the **mrgscr** utility thus obtaining a
single screening file that can used for self-energy calculations. Note that
[[freqremax]], [[freqremin]] and [[nfreqre]] define the spacing of the
frequency mesh along the real axis.
""",
),

Variable(
    abivarname="freqspmax",
    varset="gw",
    vartype="real",
    topics=['SelfEnergy_basic'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="FREQuencies for the SPectral function MAXimum",
    characteristics=['[[ENERGY]]'],
    requires="[[optdriver]] == 4",
    added_in_version="before_v9",
    text=r"""
[[freqspmax]] sets the maximum real frequency used to calculate the spectral
function from the GW Green's function. [[freqspmin]], [[freqspmax]] and
[[nfreqsp]] define the spacing of an equidistant frequency mesh along the real
axis. Alternatively, the variables [[gw_customnfreqsp]] and [[gw_freqsp]] can
be used to make a user-defined grid.
""",
),

Variable(
    abivarname="freqspmin",
    varset="gw",
    vartype="real",
    topics=['SelfEnergy_basic'],
    dimensions="scalar",
    defaultval="-[[freqspmax]]",
    mnemonics="FREQuencies for the SPectral function MINimum",
    characteristics=['[[ENERGY]]'],
    requires="[[optdriver]] == 4",
    added_in_version="before_v9",
    text=r"""
[[freqspmin]] sets the minimum real frequency used to calculate the spectral
function from the GW Green's function. [[freqspmin]] is set to -[[freqspmax]]
if left undefined. [[freqspmin]], [[freqspmax]], and [[nfreqsp]] define the
spacing of an equidistant frequency mesh along the real axis. Alternatively,
the variables [[gw_customnfreqsp]] and [[gw_freqsp]] can be used to make a user-defined grid.
""",
),

Variable(
    abivarname="friction",
    varset="rlx",
    vartype="real",
    topics=['MolecularDynamics_useful'],
    dimensions="scalar",
    defaultval=0.001,
    mnemonics="internal FRICTION coefficient",
    added_in_version="before_v9",
    text=r"""
Gives the internal friction coefficient (atomic units) for Langevin dynamics
(when [[ionmov]] = 9): fixed temperature simulations with random forces.

The equation of motion is:
$M_I \frac{d^2 R_I}{dt^2}= F_I -$[[friction]]$*M_I \frac{d R_I}{dt} - F_{random,I}$,
where $F_{random,I}$ is a Gaussian random force with average zero and variance [[friction]]$*2M_IkT$.
The atomic unit of [[friction]] is Hartree*Electronic mass*(atomic unit of Time)/Bohr^2. See [[cite:Chelikowsky2000]] for additional information.
""",
),

Variable(
    abivarname="frzfermi",
    varset="dfpt",
    vartype="integer",
    topics=['DFPT_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="FReeZe FERMI energy",
    added_in_version="before_v9",
    text=r"""
Can be used to suppress artificially the first-order change of Fermi energy,
in case of Response Function calculation for metals at Q=0. If the input variable
[[frzfermi]] is set to 1, this contribution is suppressed, even though this is incorrect.
""",
),

Variable(
    abivarname="fxcartfactor",
    varset="rlx",
    vartype="real",
    topics=['TransPath_expert', 'GeoOpt_expert'],
    dimensions="scalar",
    defaultval=ValueWithUnit(units='(Bohr^2)/Hartree', value=1),
    mnemonics="Forces to (X) CARTesian coordinates FACTOR",
    added_in_version="before_v9",
    text=r"""
The forces multiplied by [[fxcartfactor]] will be treated like difference in
cartesian coordinates in the process of optimization. This is a simple preconditioner.
TO BE UPDATED See ([[ionmov]] = 2 or 22, non-zero [[optcell]]). For example, the
stopping criterion defined by [[tolmxf]] relates to these scaled stresses.
""",
),

Variable(
    abivarname="ga_algor",
    varset="rlx",
    vartype="integer",
    topics=['GeoOpt_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Genetic Algorithm - ALGOrithm selection",
    added_in_version="before_v9",
    text=r"""
Choosing method to make the structure selection. Only the enthalpy is used now
but we plan to include, energy, electronic band gap and alchemical potentials.
Right now only value of 1 (enthalpy) works.
""",
),

Variable(
    abivarname="ga_fitness",
    varset="rlx",
    vartype="integer",
    topics=['GeoOpt_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Genetic Algorithm FITNESS function selection",
    added_in_version="before_v9",
    text=r"""
Different methodologies to perform the roulette-wheel selection of parents.
Even though, the objective function is the crystalline enthalpy (H_i), the
weight of the population elements to be chosen from in a roulette-wheel
selection can be given through different functions. We consider the following cases.

1. F = H_i / Sum H_i
2. F = exp(-(H_i-H_min)) / Sum exp(-(H_i-H_min))
3. F = (1/n_i) / Sum (1/n_i). Where n_i is the position in the ordered list of enthalpies
""",
),

Variable(
    abivarname="ga_n_rules",
    varset="rlx",
    vartype="integer",
    topics=['GeoOpt_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Genetic Algorithm Number of RULES",
    added_in_version="before_v9",
    text=r"""
Different genetic rules have been implemented and the user has the change to
choose between any of them. Right now we have 4 rules. See [[ga_rules]]
""",
),

Variable(
    abivarname="ga_opt_percent",
    varset="rlx",
    vartype="real",
    topics=['GeoOpt_expert'],
    dimensions="scalar",
    defaultval=0.2,
    mnemonics="Genetic Algorithm OPTimal PERCENT",
    added_in_version="before_v9",
    text=r"""
Percentage of the population that according to the fitness function passes to the following iteration.
""",
),

Variable(
    abivarname="ga_rules",
    varset="rlx",
    vartype="integer",
    topics=['GeoOpt_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Genetic Algorithm RULES",
    added_in_version="before_v9",
    text=r"""
Different genetic rules have been implemented and the user has the change to
choose between any of them. The chosen number of rules have been defined in [[ga_n_rules]]

Implemented rules are
1) crossover. Two parents are randomly chosen and two springs are mixed from
the two by (a) choosing randomly (through Fitness function) two parents and
then randomly rotating and shifting the coordinates withing that particular
cell. (b) Slice every one of the unit cell of the parents along a random
direction and creating the spring offs from the pieces of the two parents.
2) Vector flip mutation. From the coordinates from a given parent, a piece of
it is inverted.
3) Random strain. A random anisotropic deformation is given to the unit cell.
4) Coordinates mutation of 1/4 of the whole coordinates.
""",
),

Variable(
    abivarname="genafm",
    varset="geo",
    vartype="real",
    topics=['spinpolarisation_useful', 'SmartSymm_useful'],
    dimensions=[3],
    defaultval=MultipleValue(number=3, value=0),
    mnemonics="GENerator of the translation for Anti-FerroMagnetic space group",
    added_in_version="before_v9",
    text=r"""
This input variable might be used to define a Shubnikov type IV magnetic space
group (anti-ferromagnetic space group). The user is advised to consult [[cite:Bradley1972]]
A Shubnikov type IV magnetic space group might be defined by its Fedorov space
group (set of spatial symmetries, that do not change the magnetization), and
one translation associated with a change of magnetization. [[genafm]] is
precisely this translation, in reduced coordinates (like [[xred]])
Thus, one way to specify a Shubnikov IV magnetic space group, is to define
both [[spgroup]] and [[genafm]]. Alternatively, one might define [[spgroup]]
and [[spgroupma]], or define by hand the set of symmetries, using [[symrel]], [[tnons]] and [[symafm]].
""",
),

Variable(
    abivarname="get1den",
    varset="files",
    vartype="integer",
    topics=['nonlinear_useful', 'ElPhonInt_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET the first-order density from _1DEN file",
    added_in_version="before_v9",
    text=r"""
Relevant only for non self consistent RF calculations (e.g. to get electron
phonon matrix elements) or for non linear RF calculations (to get mixed higher
order derivatives you need several perturbed densities and wave functions).
Indicate the files from which first-order densities must be obtained, in
multi-dataset mode (in single dataset mode, use [[ird1den]]).
NOTE: a negative value of a "get" variable indicates the number of datasets
to go backwards; it is not the number to be subtracted from the current
dataset to find the proper dataset. As an example:

      ndtset 3   jdtset 1 2 4  getXXX -1

refers to dataset 2 when dataset 4 is initialized.
""",
),

Variable(
    abivarname="get1wf",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET the first-order wavefunctions from _1WF file",
    added_in_version="before_v9",
    text=r"""
Eventually used when [[ndtset]] > 0 (in the multi-dataset mode), to indicate
starting wavefunctions, as an alternative to [[ird1wf]].
One should first read the explanations given for these latter variables.
This variable is  typically used to chain the calculations in the multi-dataset mode, since they
describe from which dataset the OUTPUT wavefunctions are to be taken, as INPUT
wavefunctions of the present dataset.
See also discussion in [[getwfk]].
""",
),

Variable(
    abivarname="getbscoup",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET the Bethe-Salpeter COUPling block from...",
    added_in_version="before_v9",
    text=r"""
Eventually used when [[ndtset]] > 0 (multi-dataset mode) and, in the case of a
Bethe-Salpeter calculation to indicate that the starting coupling block of the
excitonic Hamiltonian will be taken from the output of a previous dataset. It
is used to chain the calculations, since it describes from which dataset the
OUTPUT coupling block is to be taken, as INPUT of the present dataset.
If [[getbscoup]] == 0, no such use of previously computed coupling block file is done.
If [[getbscoup]] is positive, its value gives the index of the dataset to be used as input.
If [[getbscoup]] is -1, the output of the previous dataset must be taken,
which is a frequently occurring case.
If [[getbscoup]] is a negative number, it indicates the number of datasets to
go backward to find the needed file. In this case, if one refers to a non
existent data set (prior to the first), the coupling block is not initialised
from a disk file, so that it is as if [[getbscoup]] = 0 for that initialisation.
""",
),

Variable(
    abivarname="getbseig",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET the Bethe-Salpeter EIGenstates from...",
    added_in_version="before_v9",
    text=r"""
Eventually used when [[ndtset]] > 0 (multi-dataset mode) and, in the case of a
Bethe-Salpeter calculation to indicate that the starting excitonic eigenstates
are to be taken from the output of a previous dataset. It is used to chain the
calculations, since it describes from which dataset the OUTPUT eigenstates are
to be taken, as INPUT eigenstates of the present dataset.
If [[getbseig]] == 0, no such use of previously computed output eigenstates file is done.
If [[getbseig]] is positive, its value gives the index of the dataset from
which the output states is to be used as input.
If [[getbseig]] is -1, the output eigenstates of the previous dataset must be
taken, which is a frequently occurring case.
If [[getbseig]] is a negative number, it indicates the number of datasets to
go backward to find the needed file. In this case, if one refers to a non
existent data set (prior to the first), the eigenstates are not initialised
from a disk file, so that it is as if [[getbseig]] = 0 for that initialisation.
""",
),

Variable(
    abivarname="getbsreso",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET the Bethe-Salpeter RESOnant block from...",
    added_in_version="before_v9",
    text=r"""
Eventually used when [[ndtset]] > 0 (multi-dataset mode) and, in the case of a
Bethe-Salpeter calculation to indicate that the starting resonant block of the
excitonic Hamiltonian will be taken from the output of a previous dataset. It
is used to chain the calculations, since it describes from which dataset the
OUTPUT resonant block is to be taken, as INPUT of the present dataset.
If [[getbsreso]] == 0, no such use of previously computed resonant block file is done.
If [[getbsreso]] is positive, its value gives the index of the dataset to be used as input.
If [[getbsreso]] is -1, the output of the previous dataset must be taken,
which is a frequently occurring case.
If [[getbsreso]] is a negative number, it indicates the number of datasets to
go backward to find the needed file. In this case, if one refers to a non
existent data set (prior to the first), the resonant block is not initialised
from a disk file, so that it is as if [[getbsreso]] = 0 for that initialisation.
""",
),

Variable(
    abivarname="getcell",
    varset="rlx",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET CELL parameters from...",
    added_in_version="before_v9",
    text=r"""
This variable is typically used to chain the calculations, in the multi-
dataset mode ([[ndtset]] > 0), since it describes from which dataset [[acell]]
and [[rprim]] are to be taken, as input of the present dataset. The cell
parameters are [[EVOLVING]] variables, for which such a chain of calculations is useful.
If 0, no previously computed values are used.
If >0, the value must be the index of the dataset from which the
cell data is to be used as input data. It must be the index of a dataset already
computed in the SAME run.
If equal to -1, the output data of the previous dataset must be taken, which
is a frequently occurring case. However, if the first dataset is treated, -1
is equivalent to 0, since no dataset has been computed in the same run.
If another negative number, it indicates the number of datasets to go backward
to find the needed data (once again, going back beyond the first dataset is
equivalent to using a null get variable).
""",
),

Variable(
    abivarname="getddb",
    varset="files",
    vartype="integer",
    topics=['ElPhonInt_basic', 'TDepES_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET the DDB from...",
    added_in_version="before_v9",
    text=r"""
This variable should be used when performing electron-phonon or temperature-dependent calculations in semiconductors
with the legacy implementation that computes the e-ph matrix elements at the end of the DFPT run
(for the new EPH code, see [[eph_task]]).

More detailed explanation:

The Born effective charge as well as the dielectric
tensor will be read from a previous DFPT calculations of the electric field at
q=Gamma. The use of this variable will trigger the cancellation of a residual
dipole that leads to an unphysical divergence of the GKK with vanishing
q-points. The use of this variable greatly improves the k-point convergence
speed as the density of the k-point grid required to obtain the fulfillment of
the charge neutrality sum rule is usually prohibitively large.
If [[getddb]] == 0, no such use of previously computed Born effective charge and dielectric tensor is done.
If [[getddb]] is positive, its value gives the index of the dataset from which the output density is to be used as input.
If [[getddb]] is -1, the output density of the previous dataset must be taken, which is a frequently occurring case.
If [[getddb]] is a negative number, it indicates the number of datasets to go backward to find the needed file.

NOTE: a negative value of a "get" variable indicates the number of datasets
to go backwards; it is not the number to be subtracted from the current
dataset to find the proper dataset. As an example:

      ndtset 3   jdtset 1 2 4  getXXX -1

refers to dataset 2 when dataset 4 is initialized.

Note also that, starting Abinit v9, one can also use [[getddb_path]] to specify the path of the file directly.
""",
),

Variable(
    abivarname="getddk",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET the DDK wavefunctions from _1WF file",
    added_in_version="before_v9",
    text=r"""
Eventually used when [[ndtset]] > 0 (in the multi-dataset mode), to indicate
starting wavefunctions, as an alternative to
[[irdwfk]],[[irdwfq]],[[ird1wf]],[[irdddk]]. One should first read the
explanations given for these latter variables.
The **getwfk**, **getwfq**, **get1wf** and [[getddk]] variables are
typically used to chain the calculations in the multi-dataset mode, since they
describe from which dataset the OUTPUT wavefunctions are to be taken, as INPUT
wavefunctions of the present dataset.

We now focus on the **getwfk** input variable (the only one used in ground-
state calculations), but the rules for **getwfq** and **get1wf** are similar,
with _WFK replaced by _WFQ or _1WF.
If **getwfk** ==0, no use of previously computed output wavefunction file
appended with _DSx_WFK is done.
If **getwfk** is positive, its value gives the index of the dataset for which
the output wavefunction file appended with _WFK must be used.
If **getwfk** is -1, the output wf file with _WFK of the previous dataset must
be taken, which is a frequently occurring case.
If **getwfk** is a negative number, it indicates the number of datasets to go
backward to find the needed wavefunction file. In this case, if one refers to
a non existent data set (prior to the first), the wavefunctions are not
initialised from a disk file, so that it is as if **getwfk** =0 for that
initialisation. Thanks to this rule, the use of **getwfk** -1 is rather
straightforward: except for the first wavefunctions, that are not initialized
by reading a disk file, the output wavefunction of one dataset is input of the
next one.
In the case of a ddk calculation in a multi dataset run, in order to compute
correctly the localisation tensor, it is mandatory to declare give getddk the
value of the current dataset (i.e. getddk3 3 ) - this is a bit strange and
should be changed in the future.
NOTE: a negative value of a "get" variable indicates the number of datasets
to go backwards; it is not the number to be subtracted from the current
dataset to find the proper dataset. As an example:

      ndtset 3   jdtset 1 2 4  getXXX -1

refers to dataset 2 when dataset 4 is initialized.
""",
),

Variable(
    abivarname="getdelfd",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET the 1st derivative of wavefunctions with respect to ELectric FielD, from _1WF file",
    added_in_version="before_v9",
    text=r"""
Eventually used when [[ndtset]] > 0 (in the multi-dataset mode), to indicate
starting wavefunctions, as an alternative to
[[irdwfk]],[[irdwfq]],[[ird1wf]],[[irdddk]]. One should first read the
explanations given for these latter variables.
The **getwfk**, **getwfq**, **get1wf** and [[getddk]] variables are
typically used to chain the calculations in the multi-dataset mode, since they
describe from which dataset the OUTPUT wavefunctions are to be taken, as INPUT
wavefunctions of the present dataset.

We now focus on the **getwfk** input variable (the only one used in ground-
state calculations), but the rules for **getwfq** and **get1wf** are similar,
with _WFK replaced by _WFQ or _1WF.
If **getwfk** ==0, no use of previously computed output wavefunction file
appended with _DSx_WFK is done.
If **getwfk** is positive, its value gives the index of the dataset for which
the output wavefunction file appended with _WFK must be used.
If **getwfk** is -1, the output wf file with _WFK of the previous dataset must
be taken, which is a frequently occurring case.
If **getwfk** is a negative number, it indicates the number of datasets to go
backward to find the needed wavefunction file. In this case, if one refers to
a non existent data set (prior to the first), the wavefunctions are not
initialised from a disk file, so that it is as if **getwfk** =0 for that
initialisation. Thanks to this rule, the use of **getwfk** -1 is rather
straightforward: except for the first wavefunctions, that are not initialized
by reading a disk file, the output wavefunction of one dataset is input of the
next one.
In the case of a ddk calculation in a multi dataset run, in order to compute
correctly the localisation tensor, it is mandatory to declare give getddk the
value of the current dataset (i.e. getddk3 3 ) - this is a bit strange and
should be changed in the future.
NOTE: a negative value of a "get" variable indicates the number of datasets
to go backwards; it is not the number to be subtracted from the current
dataset to find the proper dataset. As an example:

      ndtset 3   jdtset 1 2 4  getXXX -1

refers to dataset 2 when dataset 4 is initialized.
""",
),

Variable(
    abivarname="getdkdk",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET the 2nd derivative of wavefunctions with respect to K, from _1WF file",
    added_in_version="before_v9",
    text=r"""
Eventually used when [[ndtset]] > 0 (in the multi-dataset mode), to indicate
starting wavefunctions, as an alternative to
[[irdwfk]],[[irdwfq]],[[ird1wf]],[[irdddk]]. One should first read the
explanations given for these latter variables.
The **getwfk**, **getwfq**, **get1wf** and [[getddk]] variables are
typically used to chain the calculations in the multi-dataset mode, since they
describe from which dataset the OUTPUT wavefunctions are to be taken, as INPUT
wavefunctions of the present dataset.

We now focus on the **getwfk** input variable (the only one used in ground-
state calculations), but the rules for **getwfq** and **get1wf** are similar,
with _WFK replaced by _WFQ or _1WF.
If **getwfk** ==0, no use of previously computed output wavefunction file
appended with _DSx_WFK is done.
If **getwfk** is positive, its value gives the index of the dataset for which
the output wavefunction file appended with _WFK must be used.
If **getwfk** is -1, the output wf file with _WFK of the previous dataset must
be taken, which is a frequently occurring case.
If **getwfk** is a negative number, it indicates the number of datasets to go
backward to find the needed wavefunction file. In this case, if one refers to
a non existent data set (prior to the first), the wavefunctions are not
initialised from a disk file, so that it is as if **getwfk** =0 for that
initialisation. Thanks to this rule, the use of **getwfk** -1 is rather
straightforward: except for the first wavefunctions, that are not initialized
by reading a disk file, the output wavefunction of one dataset is input of the
next one.
In the case of a ddk calculation in a multi dataset run, in order to compute
correctly the localisation tensor, it is mandatory to declare give getddk the
value of the current dataset (i.e. getddk3 3) - this is a bit strange and
should be changed in the future.
NOTE: a negative value of a "get" variable indicates the number of datasets
to go backwards; it is not the number to be subtracted from the current
dataset to find the proper dataset. As an example:

      ndtset 3   jdtset 1 2 4  getXXX -1

refers to dataset 2 when dataset 4 is initialized.
""",
),

Variable(
    abivarname="getdkde",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET the mixed 2nd derivative of wavefunctions with respect to K and electric field, from _1WF file",
    added_in_version="before_v9",
    text=r"""
Eventually used when [[ndtset]] > 0 (in the multi-dataset mode), to indicate
starting wavefunctions, as an alternative to
[[irdwfk]],[[irdwfq]],[[ird1wf]],[[irdddk]]. One should first read the
explanations given for these latter variables.
The **getwfk**, **getwfq**, **get1wf** and [[getddk]] variables are
typically used to chain the calculations in the multi-dataset mode, since they
describe from which dataset the OUTPUT wavefunctions are to be taken, as INPUT
wavefunctions of the present dataset.

We now focus on the **getwfk** input variable (the only one used in ground-
state calculations), but the rules for **getwfq** and **get1wf** are similar,
with _WFK replaced by _WFQ or _1WF.
If **getwfk** ==0, no use of previously computed output wavefunction file
appended with _DSx_WFK is done.
If **getwfk** is positive, its value gives the index of the dataset for which
the output wavefunction file appended with _WFK must be used.
If **getwfk** is -1, the output wf file with _WFK of the previous dataset must
be taken, which is a frequently occurring case.
If **getwfk** is a negative number, it indicates the number of datasets to go
backward to find the needed wavefunction file. In this case, if one refers to
a non existent data set (prior to the first), the wavefunctions are not
initialised from a disk file, so that it is as if **getwfk** =0 for that
initialisation. Thanks to this rule, the use of **getwfk** -1 is rather
straightforward: except for the first wavefunctions, that are not initialized
by reading a disk file, the output wavefunction of one dataset is input of the
next one.
In the case of a ddk calculation in a multi dataset run, in order to compute
correctly the localisation tensor, it is mandatory to declare give getddk the
value of the current dataset (i.e. getddk3 3 ) - this is a bit strange and
should be changed in the future.
NOTE: a negative value of a "get" variable indicates the number of datasets
to go backwards; it is not the number to be subtracted from the current
dataset to find the proper dataset. As an example:

      ndtset 3   jdtset 1 2 4  getXXX -1

refers to dataset 2 when dataset 4 is initialized.
""",
),

Variable(
    abivarname="getden",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET the DENsity from...",
    added_in_version="before_v9",
    text=r"""
Eventually used when [[ndtset]] > 0 (multi-dataset mode) and, in the case of a
ground-state calculation, if [[iscf]]<0 (non-SCF calculation), to indicate
that the starting density is to be taken from the output of a previous
dataset. It is used to chain the calculations, since it describes from which
dataset the OUTPUT density are to be taken, as INPUT density of the present dataset.

If [[getden]] == 0, no such use of previously computed output density file is done.

If [[getden]] is positive, its value gives the index of the dataset from which
the output density is to be used as input.

If [[getden]] is -1, the output density of the previous dataset must be taken,
which is a frequently occurring case.

If [[getden]] is a negative number, it indicates the number of datasets to go
backward to find the needed file. In this case, if one refers to a non
existent data set (prior to the first), the density is not initialised from a
disk file, so that it is as if [[getden]] = 0 for that initialisation. Thanks to
this rule, the use of [[getden]] -1 is rather straightforward: except for the
first density, that is not initialized by reading a disk file, the output
density of one dataset is input of the next one.
Be careful: the output density file of a run with non-zero [[ionmov]] does
not have the proper name (it has a "TIM" indication) for use as an input of an [[iscf]]<0 calculation.
One should use the output density of a [[ionmov]] == 0 run.
NOTE: a negative value of a "get" variable indicates the number of datasets
to go backwards; it is not the number to be subtracted from the current
dataset to find the proper dataset. As an example:

      ndtset 3   jdtset 1 2 4  getXXX -1

refers to dataset 2 when dataset 4 is initialized.
""",
),

Variable(
    abivarname="getdvdb",
    varset="files",
    vartype="integer",
    topics=['ElPhonInt_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET the DVDB from...",
    added_in_version="before_v9",
    text=r"""
This variable can be used when performing electron-phonon calculations with [[optdriver]] = 7
to read a DVDB file produced in a previous dataset.
For example, one can concatenate a dataset in which an initial set of DFPT potentials
on a relatively coarse q-mesh is interpolated on a denser q-mesh using [[eph_task]] = 5 and [[eph_ngqpt_fine]].

Note also that, starting Abinit v9, one can also use [[getdvdb_path]] to specify the path of the file directly.
"""
),

Variable(
    abivarname="getefmas",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET the EFfective MASses from...",
    added_in_version="before_v9",
    text=r"""
Eventually used when [[ndtset]] > 0 (multi-dataset mode).
Only relevant for [[optdriver]]=7 and [[eph_task]]=6.
If set to 1, take the data from a _EFMAS file as input. The latter must have been produced using [[prtefmas]].
""",
),

Variable(
    abivarname="getgam_eig2nkq",
    varset="dev",
    vartype="integer",
    topics=['multidtset_useful', 'TDepES_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET the GAMma phonon data EIG2NKQ from dataset",
    requires="[[ieig2rf]] != 0 and [[qpt]] != (0.0, 0.0, 0.0)",
    added_in_version="before_v9",
    text=r"""
Relevant for second-order eigenvalue calculations using response-functions
([[ieig2rf]] != 0), and only for non-zero wavevectors [[qpt]].
From the electron-phonon matrix elements at some wavevector only, it is not
possible to determine the Debye-Waller contribution: one has to know also the
q=Gamma electron-phonon matrix elements.
The variable [[getgam_eig2nkq]] allows one to transmit the information about the
second-order derivatives of the eigenvalues for q=Gamma from the dataset where
the calculation at Gamma was done, to the datasets for other wavevectors.
""",
),

Variable(
    abivarname="gethaydock",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET the HAYDOCK restart file from...",
    added_in_version="before_v9",
    text=r"""
Eventually used when [[ndtset]] > 0 (multi-dataset mode) and, in the case of a
Bethe-Salpeter calculation to indicate that the Haydock iterative technique
will be restarted from the output of a previous dataset.
If [[gethaydock]] == 0, no such use of previously computed coupling block file is done.
If [[gethaydock]] is positive, its value gives the index of the dataset to be used as input.
If [[gethaydock]] is -1, the output of the previous dataset must be taken,
which is a frequently occurring case.
If [[gethaydock]] is a negative number, it indicates the number of datasets to
go backward to find the needed file. In this case, if one refers to a non
existent data set (prior to the first), the coupling block is not initialised
from a disk file, so that it is as if [[gethaydock]] = 0 for that initialisation.
""",
),

Variable(
    abivarname="getocc",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET OCC parameters from...",
    added_in_version="before_v9",
    text=r"""
This variable is typically used to chain the calculations, in the multi-
dataset mode ([[ndtset]] > 0), since it describes from which dataset the array
[[occ]] is to be taken, as input of the present dataset. The occupation
numbers are [[EVOLVING]] variables, for which such a chain of calculations is useful.
If [[getocc]] == 0, no such use of previously computed output occupations is done.
If [[getocc]] is positive, its value gives the index of the dataset from which
the data are to be used as input data. It must be the index of a dataset
already computed in the SAME run.
If [[getocc]] is -1, the output data of the previous dataset must be taken,
which is a frequently occurring case.
If [[getocc]] is a negative number, it indicates the number of datasets to go
backward to find the needed data. In this case, if one refers to a non
existent data set (prior to the first), the date is not initialised from a
disk file, so that it is as if [[getocc]] == 0 for that initialisation.
NOTE that a non-zero [[getocc]] MUST be used with [[occopt]] == 2, so that the
number of bands has to be initialized for each k point. Of course, these
numbers of bands must be identical to the numbers of bands of the dataset from
which [[occ]] will be copied. The same is true for the number of k points.
NOTE: a negative value of a "get" variable indicates the number of datasets
to go backwards; it is not the number to be subtracted from the current
dataset to find the proper dataset. As an example:

      ndtset 3   jdtset 1 2 4  getXXX -1

refers to dataset 2 when dataset 4 is initialized.
""",
),

Variable(
    abivarname="getqps",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful', 'GW_useful', 'Susceptibility_useful', 'SelfEnergy_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET QuasiParticle Structure",
    added_in_version="before_v9",
    text=r"""
Used when [[ndtset]] > 0 (multi-dataset mode) and [[optdriver]] = 3, or 4
(screening or sigma step of a GW calculation), to indicate that the
eigenvalues and possibly the wavefunctions have to be taken from a previous
quasi-particle calculation (instead of the usual LDA starting point). This is
to achieve quasi-particle self-consistency. See also [[irdqps]]
NOTE: a negative value of a "get" variable indicates the number of datasets
to go backwards; it is not the number to be subtracted from the current
dataset to find the proper dataset. As an example:

      ndtset 3   jdtset 1 2 4  getXXX -1

refers to dataset 2 when dataset 4 is initialized.
""",
),

Variable(
    abivarname="getscr",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful', 'GW_useful', 'SelfEnergy_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET SCReening (the inverse dielectric matrix) from...",
    added_in_version="before_v9",
    text=r"""
Used when [[ndtset]] > 0 (multi-dataset mode) and [[optdriver]] = 4 (sigma step of
a GW calculation), to indicate that the dielectric matrix (_SCR file) is to be
taken from the output of a previous dataset. It is used to chain the
calculations, since it describes from which dataset the OUTPUT dielectric
matrix is to be taken, as INPUT of the present dataset.
Note also that, starting Abinit v9, one can also use [[getscr_path]] to specify the path of the file directly.

If [[getscr]] == 0, no such use of previously computed output _SCR file is done.
If [[getscr]] is positive, its value gives the index of the dataset from which
the output _SCR file is to be used as input.
If [[getscr]] is -1, the output _SCR file of the previous dataset must be
taken, which is a frequently occurring case.
If [[getscr]] is a negative number, it indicates the number of datasets to go
backward to find the needed file. In this case, if one refers to a non
existent data set (prior to the first), the _SCR file is not initialised from
a disk file, so that it is as if [[getscr]] = 0 for that initialisation.
NOTE: a negative value of a "get" variable indicates the number of datasets
to go backwards; it is not the number to be subtracted from the current
dataset to find the proper dataset. As an example:

      ndtset 3   jdtset 1 2 4  getXXX -1

refers to dataset 2 when dataset 4 is initialized.
""",
),

Variable(
    abivarname="getsuscep",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful', 'GW_useful', 'SelfEnergy_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET SUSCEPtibility (the irreducible polarizability) from...",
    added_in_version="before_v9",
    text=r"""
Used when [[ndtset]] > 0 (multi-dataset mode) and [[optdriver]] = 4 (sigma step of
a GW calculation), to indicate that the irreducible polarizability (_SUSC
file) is to be taken from the output of a previous dataset. It is used to
chain the calculations, since it describes from which dataset the OUTPUT
susceptibility is to be taken, as INPUT of the present dataset. Performing a
GW calculations starting from the _SUSC file instead of the _SCR file presents
the advantage that starting from the irreducible polarizability, one can
calculate the screened interaction using different expressions without having
to perform a screening calculation from scratch. For example, it is possible
to apply a cutoff to the Coulomb interaction in order to facilitate the
convergence of the GW correction with respect to the size of the supercell
(see [[vcutgeo]] and [[icutcoul]])
If [[getsuscep]] == 0, no such use of previously computed output _SUSC file is done.
If [[getsuscep]] is positive, its value gives the index of the dataset from
which the output _SUSC file is to be used as input.
If [[getsuscep]] is -1, the output _SUSC file of the previous dataset must be
taken, which is a frequently occurring case.
If [[getsuscep]] is a negative number, it indicates the number of datasets to
go backward to find the needed file. In this case, if one refers to a non
existent data set (prior to the first), the _SUSC file is not initialised from
a disk file, so that it is as if [[getsuscep]] = 0 for that initialisation.
NOTE: a negative value of a "get" variable indicates the number of datasets
to go backwards; it is not the number to be subtracted from the current
dataset to find the proper dataset. As an example:

      ndtset 3   jdtset 1 2 4  getXXX -1

refers to dataset 2 when dataset 4 is initialized.
""",
),

Variable(
    abivarname="getvel",
    varset="rlx",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET VEL from...",
    added_in_version="before_v9",
    text=r"""
These variables are typically used to chain the calculations, in the multi-
dataset mode ([[ndtset]] > 0) since they describe from which dataset the
corresponding output variables are to be taken, as input of the present
dataset. The atomic positions and velocities are [[EVOLVING]] variables, for
which such a chain of calculation is useful.
Note that the use of [[getxcart]] and [[getxred]] differs when [[acell]] and
[[rprim]] are different from one dataset to the other.
If 0, no previously computed values are used.
If >0, the integer should correspond to the index of the dataset from which the VEL data should be used. It must be the index of a dataset already computed in the SAME run.
If equal to -1, the output data of the previous dataset is taken, which
is a frequently occurring case. However, if the first dataset is treated, -1
is equivalent to 0, since no dataset has yet been computed in the same run.
If another negative number, it indicates the number of datasets to go backward
to find the needed data (once again, going back beyond the first dataset is
equivalent to using a null get variable).
Note: [[getxred]] and [[getxcart]] cannot be simultaneously non-zero for the
same dataset. On the other hand the use of [[getvel]] with [[getxred]] is
allowed, despite the different coordinate system.
""",
),

Variable(
    abivarname="getwfk",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET the wavefunctions from _WFK file",
    added_in_version="before_v9",
    text=r"""
Eventually used when [[ndtset]] > 0 (in the multi-dataset mode), to indicate
starting wavefunctions, as an alternative to [[irdwfk]],.
Note also that, starting Abinit v9, one can also use [[getwfk_path]] to specify the path of the file directly.

The [[getwfk]], **getwfq**, **get1wf** and **getddk** variables are typically
used to chain the calculations in the multi-dataset mode, since they describe
from which dataset the OUTPUT wavefunctions are to be taken, as INPUT wavefunctions of the present dataset.

We now focus on the [[getwfk]] input variable (the only one used in ground-state calculations),
but the rules for **getwfq** and **get1wf** are similar, with _WFK replaced by _WFQ or _1WF.
If [[getwfk]] == 0, no use of previously computed output wavefunction file appended with _DSx_WFK is done.
If [[getwfk]] is positive, its value gives the index of the dataset for which
the output wavefunction file appended with _WFK must be used.
If [[getwfk]] is -1, the output wf file with _WFK of the previous dataset must
be taken, which is a frequently occurring case.
If [[getwfk]] is a negative number, it indicates the number of datasets to go
backward to find the needed wavefunction file. In this case, if one refers to
a non existent data set (prior to the first), the wavefunctions are not
initialised from a disk file, so that it is as if [[getwfk]] = 0 for that
initialisation. Thanks to this rule, the use of [[getwfk]] -1 is rather
straightforward: except for the first wavefunctions, that are not initialized
by reading a disk file, the output wavefunction of one dataset is input of the
next one.
NOTE: a negative value of a "get" variable indicates the number of datasets
to go backwards; it is not the number to be subtracted from the current
dataset to find the proper dataset. As an example:

      ndtset 3   jdtset 1 2 4  getXXX -1

refers to dataset 2 when dataset 4 is initialized.
""",
),

Variable(
    abivarname="getwfkfine",
    varset="dev",
    vartype="integer",
    topics=['multidtset_useful', 'DFPT_useful', 'TDepES_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET the fine grid wavefunctions from _WFK file",
    added_in_version="before_v9",
    text=r"""
Eventually used when [[ndtset]] > 0 (in the multi-dataset mode), to indicate
starting wavefunctions, as an alternative to [[irdwfkfine]]. One should first
read the explanations given for these latter variables.
The [[getwfkfine]] variables is typically used to chain the calculations in
the multi-dataset mode, since they describe from which dataset the OUTPUT
wavefunctions are to be taken, as INPUT wavefunctions of the present dataset.
If [[getwfkfine]] == 0, no use of previously computed output wavefunction file
appended with _DSx_WFK is done.
If [[getwfkfine]] is positive, its value gives the index of the dataset for
which the output wavefunction file appended with _WFK must be used.
If [[getwfkfine]] is -1, the output wf file with _WFK of the previous dataset
must be taken, which is a frequently occurring case.
If [[getwfkfine]] is a negative number, it indicates the number of datasets to
go backward to find the needed wavefunction file. In this case, if one refers
to a non existent data set (prior to the first), the wavefunctions are not
initialised from a disk file, so that it is as if [[getwfkfine]] = 0 for that
initialisation. Thanks to this rule, the use of [[getwfkfine]] -1 is rather
straightforward: except for the first wavefunctions, that are not initialized
by reading a disk file, the output wavefunction of one dataset is input of the
next one.
NOTE: a negative value of a "get" variable indicates the number of datasets
to go backwards; it is not the number to be subtracted from the current
dataset to find the proper dataset. As an example:

     ndtset 3   jdtset 1 2 4  getXXX -1

refers to dataset 2 when dataset 4 is initialized. Response-function calculation:

  * one and only one of [[getwfkfine]] or [[irdwfkfine]] MUST be non-zero
  * if [[getwfkfine]] = 1: read ground state k -wavefunctions from a disk file appended with _WFK,
    produced in a previous ground state calculation.
  * Reading the fine grid wavefunction will trigger the k-points interpolation
  technique of the temperature dependent calculations.

Bethe-Salpeter calculation:

  * one and only one of [[getwfkfine]] or [[irdwfkfine]] MUST be non-zero
  * if [[getwfkfine]] = 1: read ground state k -wavefunctions from a disk file appended with _WFK,
    produced in a previous ground state calculation
  * This variable or [[irdwfkfine]] is mandatory when [[bs_interp_mode]] == 1

For further information about the *files file*, consult the [[help:abinit#files-file]].

**This variable is experimental. In development.**
""",
),

Variable(
    abivarname="getwfq",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET the wavefunctions from _WFQ file",
    added_in_version="before_v9",
    text=r"""
Eventually used when [[ndtset]] > 0 (in the multi-dataset mode), to indicate
starting wavefunctions, as an alternative to [[irdwfq]].
Note also that, starting Abinit v9, one can also use [[getwfq_path]] to specify the path of the file directly.

The **getwfk**, [[getwfq]], **get1wf** and **getddk** variables are typically
used to chain the calculations in the multi-dataset mode, since they describe
from which dataset the OUTPUT wavefunctions are to be taken, as INPUT wavefunctions of the present dataset.
See discussion in [[getwfk]]
""",
),

Variable(
    abivarname="getxcart",
    varset="rlx",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET XCART from...",
    added_in_version="before_v9",
    text=r"""
These variables are typically used to chain the calculations, in the multi-
dataset mode ([[ndtset]] > 0) since they describe from which dataset the
corresponding output variables are to be taken, as input of the present
dataset. The atomic positions and velocities are [[EVOLVING]] variables, for
which such a chain of calculation is useful.
Note that the use of [[getxcart]] and [[getxred]] differs when [[acell]] and
[[rprim]] are different from one dataset to the other.
If 0, no previously computed values are used.
If >0, the integer must correspond to the index of the dataset from which the
data are to be used as input data. It must be the index of a dataset already
computed in the SAME run.
If equal to -1, the output data of the previous dataset must be taken, which
is a frequently occurring case. However, if the first dataset is treated, -1
is equivalent to 0, since no dataset has yet been computed in the same run.
If another negative number, it indicates the number of datasets to go backward
to find the needed data (once again, going back beyond the first dataset is
equivalent to using a null get variable).
Note: [[getxred]] and [[getxcart]] cannot be simultaneously non-zero for the
same dataset. On the other hand the use of [[getvel]] with [[getxred]] is
allowed, despite the different coordinate system.
""",
),

Variable(
    abivarname="getxred",
    varset="rlx",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GET XRED from...",
    added_in_version="before_v9",
    text=r"""
These variables are typically used to chain the calculations, in the multi-
dataset mode ([[ndtset]] > 0) since they describe from which dataset the
corresponding output variables are to be taken, as input of the present
dataset. The atomic positions and velocities are [[EVOLVING]] variables, for
which such a chain of calculation is useful.
Note that the use of [[getxcart]] and [[getxred]] differs when [[acell]] and
[[rprim]] are different from one dataset to the other.
If 0, no use of previously computed values must occur.
If >0, the integer must correspond to the index of the dataset from which the
data are to be used as input data. It must be the index of a dataset already
computed in the SAME run.
If equal to -1, the output data of the previous dataset must be taken, which
is a frequently occurring case. However, if the first dataset is treated, -1
is equivalent to 0, since no dataset has yet been computed in the same run.
If another negative number, it indicates the number of datasets to go backward
to find the needed data (once again, going back beyond the first dataset is
equivalent to using a null get variable).
Note: [[getxred]] and [[getxcart]] cannot be simultaneously non-zero for the
same dataset. On the other hand the use of [[getvel]] with [[getxred]] is
allowed, despite the different coordinate system.
""",
),

Variable(
    abivarname="goprecon",
    varset="rlx",
    vartype="integer",
    topics=['GeoOpt_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Geometry Optimization PRECONditioner equations",
    added_in_version="before_v9",
    text=r"""
Set the kind of preconditioner to be used for Geometry Optimization
(Note: Under development now (2011.05.20))

  * [[goprecon]] = 0: No preconditioner
  * [[goprecon]] = [1-9]: Linear preconditioner
  * [[goprecon]] = [11-19]: Non-linear preconditioner
""",
),

Variable(
    abivarname="goprecprm",
    varset="rlx",
    vartype="real",
    topics=['GeoOpt_expert'],
    dimensions=[3],
    defaultval=0,
    mnemonics="Geometry Optimization PREconditioner PaRaMeters equations",
    added_in_version="before_v9",
    text=r"""
Set the parameters use by the preconditioner to be used for Geometry
Optimization
(Note: Under development now (2011.06.06))
""",
),

Variable(
    abivarname="gpu_devices",
    varset="paral",
    vartype="integer",
    topics=['parallelism_expert'],
    dimensions=[5],
    defaultval=[-1, -1, -1, -1, -1],
    mnemonics="GPU: choice of DEVICES on one node",
    requires="[[use_gpu_cuda]] == 1 (CUDA functionality)",
    added_in_version="before_v9",
    text=r"""
To be used when several GPU devices are present on each node, assuming the
same number of devices on all nodes.
Allows to choose in which order the GPU devices are chosen and distributed
among MPI processes (see examples below). When the default value (-1) is set,
the GPU devices are chosen by order of performance (FLOPS, memory).

Examples:

  * 2 GPU devices per node, 4 MPI processes per node, **gpu_device** =[-1,-1,-1,-1,-1] (default):
MPI processes 0 and 2 use the best GPU card, MPI processes 1 and 3 use the
slowest GPU card.

  * 3 GPU devices per node, 5 MPI processes per node, **gpu_device** =[1,0,2,-1,-1]:
MPI processes 0 and 3 use GPU card 1, MPI processes 1 and 4 use GPU card 0,
MPI process 2 uses GPU card 2.

  * 3 GPU devices per node, 5 MPI processes per node, **gpu_device** =[0,1,-1,-1,-1]:
MPI processes 0, 2 and 4 use GPU card 0, MPI processes 1 and 3 use GPU card 1;
the 3rd GPU card is not used.

GPU card are numbered starting from 0; to get the GPU devices list, type
"nvidia-smi" or "lspci | grep -i nvidia".
""",
),

Variable(
    abivarname="gpu_linalg_limit",
    varset="paral",
    vartype="integer",
    topics=['parallelism_expert'],
    dimensions="scalar",
    defaultval=2000000,
    mnemonics="GPU (Cuda): LINear ALGebra LIMIT",
    requires="[[use_gpu_cuda]] == 1 (CUDA functionality)",
    added_in_version="before_v9",
    text=r"""
Use of linear algebra and matrix algebra on GPU is only efficient if the size
of the involved matrices is large enough. The [[gpu_linalg_limit]] parameter
defines the threshold above which linear (and matrix) algebra operations are
done on the Graphics Processing Unit.
The considered matrix size is equal to:

* SIZE=([[mpw]] $\times$ [[nspinor]] / [[npspinor]]) $\times$ ([[npband]] $\times$ [[bandpp]]) $^2$

When SIZE>=[[gpu_linalg_limit]], [[wfoptalg]] parameter is automatically set
to 14 which corresponds to the use of LOBPCG algorithm for the calculation of
the eigenstates.
""",
),

Variable(
    abivarname="gw_customnfreqsp",
    varset="gw",
    vartype="integer",
    topics=['FrequencyMeshMBPT_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GW CUSTOM FREQuencies for SPectral function",
    requires="[[optdriver]] == 4 and [[gwcalctyp]] in [2, 9, 12, 19, 22, 29]",
    added_in_version="before_v9",
    text=r"""
[[gw_customnfreqsp]] lets the user define the grid points along the real
frequency axis by hand for the calculation of the self-energy along the real
axis. Set this to the number of frequencies you want. The frequencies are specified with [[gw_freqsp]].
""",
),

Variable(
    abivarname="gw_freqsp",
    varset="gw",
    vartype="real",
    topics=['FrequencyMeshMBPT_expert'],
    dimensions=['[[gw_customnfreqsp]]'],
    defaultval=Range({'start': 1, 'stop': '[[gw_customnfreqsp]]'}),
    mnemonics="GW SPectral FREQuencies",
    requires="[[optdriver]] == 4 and [[gw_customnfreqsp]] > 0",
    added_in_version="before_v9",
    text=r"""
[[gw_freqsp]] specifies the grid points for the real frequency axis when the
real and imaginary (spectral function) parts of sigma are calculated explicitly
for post-processing or plotting. Only activated if [[gw_customnfreqsp]] is not
equal to 0. The number of frequencies is set by the value of
[[gw_customnfreqsp]]. Example:

    gw_customnfreqsp   5
    nfreqsp            5
    gw_freqsp         -0.5  -0.1  0.0  1.0  10.0 eV

If [[nfreqsp]] is not equal to [[gw_customnfreqsp]] a warning will be issued.
""",
),

Variable(
    abivarname="gw_frqim_inzgrid",
    varset="gw",
    vartype="integer",
    topics=['FrequencyMeshMBPT_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GW Contour Deformation FReQuencies on IMaginary axis Inverse Z Grid",
    requires="[[optdriver]] in [3,4] and [[gwcalctyp]] in [2, 9, 12, 19, 22, 29]",
    added_in_version="before_v9",
    text=r"""
[[gw_frqim_inzgrid]] creates grid points along the **imaginary** frequency axis
by using an equidistant grid in the variable  $z \subset [0,1]$ where the transform
is:

$$ i\omega^\prime = w_p \frac{z}{1-z}. $$

Here  $\omega_p$ is the plasma frequency (default can be overridden by setting
[[ppmfrq]]). The equidistant grid in z is determined uniquely by [[nfreqim]])
and the points are distributed so that half of them lie below the plasma frequency.
""",
),

Variable(
    abivarname="gw_frqre_inzgrid",
    varset="gw",
    vartype="integer",
    topics=['FrequencyMeshMBPT_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GW Contour Deformation FReQuencies on REal axis Inverse Z Grid",
    requires="[[optdriver]] in [3, 4] and [[gwcalctyp]] in [2, 9, 12, 19, 22, 29]",
    added_in_version="before_v9",
    text=r"""
[[gw_frqre_inzgrid]] creates grid points along the **real** frequency axis by
using an equidistant grid in the variable  $z \subset [0, 1]$ where the transform is:

$$ \omega = \omega_p \frac{z}{1-z}. $$

Here $\omega_p$ is the plasma frequency (default can be overridden by setting
[[ppmfrq]]). The equidistant grid in z is determined uniquely by [[nfreqre]] )
and the points are distributed so that half of them lie below the plasma
frequency. This is useful in conjunction with [[gw_frqim_inzgrid]] if one needs
to use a grid which maps $[0, \infty] \rightarrow [0,1]$. Note that typically _many_ more
points are needed along the real axis in order to properly resolve peak
structures. In contrast, both the screening and self-energy are very smooth
along the imaginary axis. Also, please note that this is **not** an efficient
grid for **standard** Contour Deformation calculations, where typically only a
smaller range of frequencies near the origin is required. The maximum value
needed along the real frequency axis is output in the logfile during Contour
Deformation sigma calculations.
""",
),

Variable(
    abivarname="gw_frqre_tangrid",
    varset="gw",
    vartype="integer",
    topics=['FrequencyMeshMBPT_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GW Contour Deformation FReQencies on REal axis - Use Tangent Grid",
    requires="[[optdriver]] in [3,4] and [[gwcalctyp]] in [2,9,12,19,22,29]",
    added_in_version="before_v9",
    text=r"""
[[gw_frqre_tangrid]] defines a nonuniform grid to be used in frequency, with
stepsize increasing proportional to $\tan(x)$. This makes the grid approximately
linear to start with, with a rapid increase towards the end. Also, this is the
grid which gives equal importance to each point used in the integration of a
function which decays as $1/x^2$. To be used in conjunction with [[nfreqre]],
[[cd_max_freq]] and [[cd_halfway_freq]] which determine the parameters of the transformed grid.
""",
),

Variable(
    abivarname="gw_invalid_freq",
    varset="gw",
    vartype="integer",
    topics=['Susceptibility_expert', 'SelfEnergy_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GW treatment of INVALID FREQuency for Hybertsen-Louie PPM",
    requires="[[optdriver]] in [3,4] and [[ppmodel]] in [2]",
    added_in_version="before_v9",
    text=r"""
[[gw_invalid_freq]] sets the procedure to follow when a PPM frequency is
invalid (negative or imaginary).

  * [[gw_invalid_freq]] = 0: Drop them as proposed in Appendix B of [[cite:Hybertsen1986]].
  * [[gw_invalid_freq]] = 1: Set them to 1 hartree, as done for the PPM of Godby-Needs [[cite:Godby1989]].
  * [[gw_invalid_freq]] = 2: Set them to infinity.
""",
),

Variable(
    abivarname="gw_nqlwl",
    varset="gw",
    vartype="integer",
    topics=['GW_expert', 'BSE_expert', 'Susceptibility_expert', 'SelfEnergy_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="GW, Number of Q-points for the Long Wave-Length Limit",
    requires="[[optdriver]] in [3,4,99]",
    added_in_version="before_v9",
    text=r"""
Only relevant if [[optdriver]] = 3,4,99 that is, screening, sigma or Bethe-
Salpeter calculations, although the actual meaning of the variable depends on
the particular run-level (see discussion below).

[[gw_nqlwl]] defines the number of directions in reciprocal space used to
describe the non-analytical behaviour of the heads ($G = G'=0$) and the wings
($G=0$ or $G'=0$) of the dielectric matrix in the optical limit (i.e. for $q$
tending to zero). The number of directions is specified by the additional
variable [[gw_qlwl]].

When [[optdriver]] = 3, [[gw_nqlwl]] and **gw_qlwl** define the set of "small" $q$
that will be calculated and stored in the final SCR file. Therefore, the two
variables can be used to analyze how the optical spectra depend on the
direction of the incident phonon (useful especially in anisotropic systems).

When [[optdriver]] = 4, [[gw_nqlwl]] and **gw_qlwl** can be used to specify the
heads and the wings to be used to perform the quadrature of the correlated
part of the self-energy in the small region around the origin. (NB: not yet
available, at present the quadrature is performed using a single direction in
q-space)

When [[optdriver]] = 99, [[gw_nqlwl]] and **gw_qlwl** define the set of
directions in q-space along which the macroscopic dielectric function is
evaluated. By default the Bethe-Salpeter code calculates the macroscopic
dielectric function using six different directions in q-space (the three basis
vectors of the reciprocal lattice and the three Cartesian axis).
""",
),

Variable(
    abivarname="gw_nstep",
    varset="gw",
    vartype="integer",
    topics=['GW_basic'],
    dimensions="scalar",
    defaultval=30,
    mnemonics="GW Number of self-consistent STEPs",
    requires="[[optdriver]] == 8",
    added_in_version="before_v9",
    text=r"""
Gives the maximum number of self-consistent GW cycles (or "iterations") in
which G and/or W will be updated until the quasi-particle energies are
converged within [[gw_toldfeig]]. [[gwcalctyp]] and [[gw_sctype]] are used to
define the type of self-consistency.
""",
),

Variable(
    abivarname="gw_qlwl",
    varset="gw",
    vartype="real",
    topics=['Susceptibility_expert', 'SelfEnergy_expert', 'BSE_expert'],
    dimensions=[3, '[[gw_nqlwl]]'],
    defaultval=[1e-05, 2e-05, 3e-05],
    mnemonics="GW, Q-points for the Long Wave-Length limit",
    requires="[[optdriver]] in [3,4,99]",
    added_in_version="before_v9",
    text=r"""
When [[optdriver]] = 3, [[gw_qlwl]] defines the set of q-points around Gamma
that are considered during the evaluation of the non-analytical behaviour of
the dielectric matrix. Optical spectra (with and without non-local field
effects) are evaluated for each direction specified by [[gw_qlwl]].
""",
),

Variable(
    abivarname="gw_qprange",
    varset="gw",
    vartype="integer",
    topics=['SelfEnergy_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GW QuasiParticle RANGE policy",
    requires="[[optdriver]] == 4",
    added_in_version="before_v9",
    text=r"""
[[gw_qprange]] is active only when [[nkptgw]] is equal to zero (default
value). This variable simplifies the specification of the list of kpoints and
of the bands to be used for the computation of the quasi-particle corrections.
The possible values are:

  * 0 --> Compute the QP corrections only for the fundamental and the optical gap
  * +num --> Compute the QP corrections for all the k-points in the irreducible zone,
     and include `num` bands above and below the Fermi level.
  * -num --> Compute the QP corrections for all the k-points in the irreducible zone.
     Include all occupied states and `num` empty states.

The default value is 0 and is very handy for one-shot calculations. It is
important to stress, however, that the position of the optical/fundamental
gaps is deduced from the energies computed on the k-mesh used for the WFK
file. Therefore the computed gaps might differ from the correct ones that can
only be obtained with an appropriate sampling of the irreducible zone.
Positive values are useful if we do not know the position of the GW HOMO, LOMO
and we want to investigate the effect of the GW corrections on the states
close to the gap Negative values are usually used for self-consistent
calculations Note that, in the case of self-consistency or [[symsigma]] == 1, the
code might change the bands range so that all the degenerate states are
included. Note also that [[kptgw]], and [[bdgw]] are ignored when this options
is used. If you want to select manually the list of k-points and bands, you
have to provide the three variables [[nkptgw]], [[kptgw]], and [[bdgw]].
""",
),

Variable(
    abivarname="gw_sctype",
    varset="gw",
    vartype="integer",
    topics=['GW_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="GW, Self-Consistency TYPE",
    requires="[[optdriver]] in [3,4]",
    added_in_version="before_v9",
    text=r"""
This variable is used to partially define the kind of self-consistency for GW
calculations. The other piece of information is given by [[gwcalctyp]] that
defines the particular approximation for the self-energy operator as well as
whether the wavefunctions have to replaced by quasi-particle amplitudes.

If [[gw_sctype]] is specified in the input file, the code will perform an
iterative update of the quantities entering the GW equations until the quasi-
particle energies are converged within [[gw_toldfeig]]. The maximum number of
iterations is specified by [[gw_nstep]]. Possible values are:

  * 1 --> standard one-shot method (one screening calculation followed by a single sigma run)
  * 2 --> self-consistency only on W (iterative update of W followed by a sigma run in which G is approximated with the Kohn-Sham independent-particle Green's function G0)
  * 3 --> self-consistency only of G (a single screening calculation to obtain the Kohn-Sham polarizability followed by an iterative update of the Green's functions in the self-energy)
  * 4 --> fully self-consistent algorithm (iterative update of both G and W)

It is possible to initialize the self-consistent procedure by reading a
previously calculated SCR or SUSC file via the variables [[getscr]] or
[[getsuscep]], respectively. [[getqps]] can be used to read a previous QPS
file thus initializing the Green functions to be used in the first self-
consistent iteration.
""",
),

Variable(
    abivarname="gw_sigxcore",
    varset="gw",
    vartype="integer",
    topics=['SelfEnergy_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GW, SIGma (self-energy) for the CORE contribution",
    requires="[[optdriver]] == 4 and [[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
Only available for PAW and relevant if [[optdriver]] = 4 that is, sigma
calculations.

Theoretical introduction: GW calculations performed on top of electronic
calculations relying when the frozen-core approximation is used to separate
inner-core electrons from valence electrons, only the contribution to the
self-energy arising from valence electrons is explicitly accounted for. In the
standard approach based on pseudopotentials the contribution to the self-
energy due to core electrons is approximated by means of the KS exchange-
correlation potential generated by the core density. In the case of GW
calculations employing the PAW method, the core contribution to the self-
energy can be more accurately estimated in terms of the Fock operator
generated by the core wavefunctions. In the simplest approach, the only
ingredients required for this more refined treatment are the wave functions of
the core electrons in the reference atomic configuration that are calculated
during the generation of the PAW setup. This is a good approximation provided
that the core wave functions are strictly localized inside the PAW spheres.

[[gw_sigxcore]] defines the approximation used to evaluate the core
contribution to sigma.

  * [[gw_sigxcore]] = 0, standard approach, the core contribution is approximated with vxc.
  * [[gw_sigxcore]] = 1, the core term is approximated with the Fock operator inside the PAW spheres.
""",
),

Variable(
    abivarname="gw_toldfeig",
    varset="gw",
    vartype="real",
    topics=['GW_basic'],
    dimensions="scalar",
    defaultval=ValueWithUnit(units='eV', value=0.1),
    mnemonics="GW TOLerance on the DiFference of the EIGenvalues",
    characteristics=['[[ENERGY]]'],
    requires="[[optdriver]] == 8",
    added_in_version="before_v9",
    text=r"""
Sets a tolerance for absolute differences of QP energies that will cause one
self-consistent GW cycle to stop.
Can be specified in Ha (the default), Ry, eV or Kelvin, since **toldfe** has
the [[ENERGY]] characteristics (1 Ha = 27.2113845 eV)
""",
),

Variable(
    abivarname="gwcalctyp",
    varset="gw",
    vartype="integer",
    topics=['GW_basic', 'SelfEnergy_basic', 'RPACorrEn_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GW CALCulation TYPe",
    requires="[[optdriver]] in [3,4]",
    added_in_version="before_v9",
    text=r"""
**gwcalctyp** governs the choice between the different capabilities of the GW code.

  * 0 <= **gwcalctyp** <= 9: standard "1 shot" quasi-particle method.
  * 10 <= **gwcalctyp** <= 19: self-consistent quasi-particle method on energies only.
  * 20 <= **gwcalctyp** <= 29: self-consistent quasi-particle method on energies and wavefunctions.

  * **gwcalctyp** = 0, 10, or 20: standard Plasmon-Pole model GW calculation.
  * **gwcalctyp** = 1: GW calculation where the self-energy along the real axis is obtained by performing the analytic continuation from the imaginary axis to the full complex plane via the Pade approximant. Only available for standard "1 shot" quasi-particle method.
  * **gwcalctyp** = 2, 12, or 22: GW calculation using numerical integration (contour deformation method, see e.g. [[cite:Lebegue2003]]).
  * **gwcalctyp** = 5, 15, or 25: Hybrid functional or Hartree-Fock calculation, with the identifier of the functional given by [[ixc_sigma]]. See the latter for the definition of other related variables.
  * **gwcalctyp** = 6, 16, or 26: Screened Exchange calculation.
  * **gwcalctyp** = 7, 17, or 27: COHSEX calculation.
  * **gwcalctyp** = 8, 18, or 28: model GW calculation following [[cite:Faleev2004]] using a Plasmon-Pole model.
  * **gwcalctyp** = 9, 19, or 29: model GW calculation following [[cite:Faleev2004]] using numerical integration (contour deformation method).
""",
),

Variable(
    abivarname="gwcomp",
    varset="gw",
    vartype="integer",
    topics=['SelfEnergy_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GW COMPleteness",
    requires="[[optdriver]] in [3,4]",
    added_in_version="before_v9",
    text=r"""
[[gwcomp]] governs the use of an extrapolar approximation. If [[gwcomp]] == 1,
one improves the completeness in a truncated sum over states. In practice,
this permits one to reduce quite much the number of bands required in the
calculation of the screening or of the self-energy. The energy parameter
needed in the extrapolar approximation is set by [[gwencomp]]. See
[[cite:Bruneval2008]] for a description of the methodology.
""",
),

Variable(
    abivarname="gwencomp",
    varset="gw",
    vartype="real",
    topics=['SelfEnergy_useful', 'Susceptibility_useful'],
    dimensions="scalar",
    defaultval=2.0,
    mnemonics="GW ENergy for COMPleteness",
    requires="[[optdriver]] in [3,4] and [[gwcomp]] == 1",
    added_in_version="before_v9",
    text=r"""
[[gwencomp]] sets the energy parameter used in the extrapolar approximation
used to improve completeness and make the convergence against the number of
bands much faster.

See [[cite:Bruneval2008]] for a description of the methodology.
""",
),

Variable(
    abivarname="gwgamma",
    varset="gw",
    vartype="integer",
    topics=['Susceptibility_expert', 'SelfEnergy_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GW GAMMA",
    requires="[[optdriver]] = 3 or 4 (Sigma calculations)",
    added_in_version="before_v9",
    text=r"""
If [[gwgamma]] is 1, the vertex correction will be included leading to what is
known as "GW-Gamma" approximation. see R. Del Sole, L. Reining, and R. W.
Godby, Phys. Rev. B 49, 8024 (1994). Note that, in order to include the vertex
correction in W, one has to start the sigma calculation from the
susceptibility file_SUSC instead of the _SCR file (see [[getsuscep]]   and
[[irdsuscep]]  ) Not available for PAW calculations.

[[gwgamma]] = -4 activates the bootstrap kernel of Sharma et al. [[cite:Sharma2011]]
in the test-charge-test-charge dielectric function [[cite:Chen2015]].

[[gwgamma]] = -6 uses the same bootstrap kernel as with [[gwgamma]] = -4
but with only the head of the kernel. As such, the self-consistent iteration in the kernel
can be disregarded [[cite:Chen2016]].

[[gwgamma]] = -8 activates the RPA bootstrap-like kernel (one-shot) (see [[cite:Berger2015]]
and [[cite:Rigamonti2015]]).
""",
),

Variable(
    abivarname="gwls_band_index",
    varset="gw",
    vartype="integer",
    topics=['GWls_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="GWLS BAND INDEX",
    requires="[[optdriver]] == 66",
    added_in_version="before_v9",
    text=r"""
Governs the DFT eigenstate $|e\rangle$ in which the self-energy will be evaluated, as
shown in Eq. (7) of [[cite:Laflamme2015]]. That is, it is the state
to be corrected in the G0W0 scheme.
""",
),

Variable(
    abivarname="gwls_correlation",
    varset="gw",
    vartype="integer",
    topics=['GWls_expert'],
    dimensions="scalar",
    defaultval=3,
    mnemonics="GWLS CORRELATION",
    requires="[[optdriver]] == 66",
    added_in_version="before_v9",
    text=r"""
Governs the use of a dielectric model (as explained in Sec. V of
[[cite:Laflamme2015]] and the use of the Lanczos scheme to solve Eqs. (30) and
(35) of the same reference at all external [[gw_freqsp]] and integration (as
generated from [[gwls_npt_gauss_quad]]) frequencies. The different choices
are:

  * [[gwls_correlation]] == 1: GWLS calculation **with** the dielectric model and **without** the shift Lanczos technique,
  * [[gwls_correlation]] == 2: GWLS calculation **without** the dielectric model and **without** the shift Lanczos technique,
  * [[gwls_correlation]] == 3: GWLS calculation **with** the dielectric model and **with** the shift Lanczos technique,
  * [[gwls_correlation]] == 4: GWLS calculation **without** the dielectric model and **with** the shift Lanczos technique,
  * [[gwls_correlation]] == 5: Not a GWLS calculation; just calculate and print the eigenvalues of the (static) dielectric matrix (for debugging purposes).

The default, ([[gwls_correlation]] == 3), is the most performant option and
should be kept by the user. Option 1, 2 and 5 are deprecated and will be
removed.
""",
),

Variable(
    abivarname="gwls_diel_model",
    varset="gw",
    vartype="integer",
    topics=['GWls_expert'],
    dimensions="scalar",
    defaultval=2,
    mnemonics="GWLS dielectric model",
    requires="[[optdriver]] == 66",
    added_in_version="before_v9",
    text=r"""
Not used yet.
""",
),

Variable(
    abivarname="gwls_exchange",
    varset="gw",
    vartype="integer",
    topics=['GWls_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="GWLS exact EXCHANGE",
    requires="[[optdriver]] == 66",
    added_in_version="before_v9",
    text=r"""
Governs whether the exact exchange for the state to be corrected
([[gwls_band_index]]) is calculated ([[gwls_exchange]] == 1) or not
([[gwls_exchange]] = =0).
""",
),

Variable(
    abivarname="gwls_first_seed",
    varset="gw",
    vartype="integer",
    topics=['GWls_expert'],
    dimensions="scalar",
    defaultval="[[gwls_band_index]]",
    mnemonics="GWLS FIRST SEED vector",
    requires="[[optdriver]] == 66",
    added_in_version="before_v9",
    text=r"""
This variable sets the band index to be used to generate the first seed vector
to be used in the construction of the Lanczos basis for the (static)
dielectric matrix in a GWLS calculation. See Sec. IV of [[cite:Laflamme2015]].
Together with [[gwls_nseeds]], this defines the seeds for the
Lanczos procedure. That is, the states associated to band index
[[gwls_first_seed]] to [[gwls_first_seed]]+[[gwls_nseeds]]-1 are used to
generate the seed vectors.

The default [[gwls_first_seed]] == [[gwls_band_index]] and [[gwls_nseeds]] == 1
has been thoroughly tested and seems to be the most performant. Users should
therefore keep the default value.
""",
),

Variable(
    abivarname="gwls_kmax_analytic",
    varset="gw",
    vartype="integer",
    topics=['GWls_expert'],
    dimensions="scalar",
    defaultval=8,
    mnemonics="GWLS KMAX for the ANALYTIC term",
    requires="[[optdriver]] == 66",
    added_in_version="before_v9",
    text=r"""
Governs the number of iterations to be done in the shift Lanczos solution of
Eq. (35) of [[cite:Laflamme2015]] to solve it at all external
frequencies requested by the user ([[gw_freqsp]]). The default value is
converged to a few 10s of meV for all molecules studied so far.
""",
),

Variable(
    abivarname="gwls_kmax_complement",
    varset="gw",
    vartype="integer",
    topics=['GWls_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="GWLS KMAX for the COMPLEMENT space.",
    requires="[[optdriver]] == 66",
    added_in_version="before_v9",
    text=r"""
The G0W0 formalism involves the calculation of a summation conceptually linked
to the trace of the dielectric matrix [see Eq. (38) of [[cite:Laflamme2015]]\].
Since the eigenvalues spectrum of the dielectric matrix of formed by a
few large discrete eigenvalues and an integrable divergence in the density of
eigenvalues around 0, it is expensive to sample accurately this divergence
using the exact dielectric operator. It this becomes interesting to calculate
the 'trace' of the 'exact - model' dielectric matrix in a small basis and add
it to the 'trace' of the 'model' dielectric matrix obtained in a large basis.
In the context where the model dielectric matrix is used in the calculations,
[[gwls_stern_kmax]] determines the size of the 'small' basis and
[[gwls_kmax_complement]] determines the size of the 'large' basis.

For more information on the exact role of these bases and on the model
dielectric operator used, see Sec. V of [[cite:Laflamme2015]].
""",
),

Variable(
    abivarname="gwls_kmax_numeric",
    varset="gw",
    vartype="integer",
    topics=['GWls_expert'],
    dimensions="scalar",
    defaultval=16,
    mnemonics="GWLS KMAX for the NUMERIC term",
    requires="[[optdriver]] == 66",
    added_in_version="before_v9",
    text=r"""
Governs the number of iterations to be done in the shift Lanczos solution of
Eq. (30) of [[cite:Laflamme2015]] to solve it simultaneously at all
integration frequencies (generated automatically by the number of points
[[gwls_npt_gauss_quad]] to use in the gaussian quadrature) and all external
frequencies requested by the user ([[gw_freqsp]]). The default value is
converged to a few 10s of meV for all molecules studied so far.
""",
),

Variable(
    abivarname="gwls_kmax_poles",
    varset="gw",
    vartype="integer",
    topics=['GWls_expert'],
    dimensions="scalar",
    defaultval=4,
    mnemonics="GWLS KMAX for the calculation of the POLES residue",
    requires="[[optdriver]] == 66",
    added_in_version="before_v9",
    text=r"""
The contour deformation technique, in the G0W0 context, will involve the
calculation of pole residues associated to states lying between the one
corrected ([[gwls_band_index]]) and the Fermi level. These residues take the
form of a matrix element of the inverse dielectric matrix at a real frequency
[see Eq. (11) of [[cite:Laflamme2015]]\]. Therefore, the dielectric
matrix must be constructed in some basis at these frequencies and inverted to
calculate the matrix element. The present input variable sets the size of the
Lanczos basis to be constructed for this purpose. The default value has proven
to be very robust for many molecular systems and should therefore be left to
the default value by the user.

For more information on the Lanczos basis constructed for the calculation of
the residues, see Sec. IV of [[cite:Laflamme2015]].
""",
),

Variable(
    abivarname="gwls_list_proj_freq",
    varset="gw",
    vartype="real",
    topics=['GWls_expert'],
    dimensions=['[[gwls_n_proj_freq]]'],
    defaultval="*0.0",
    mnemonics="GWLS LIST of the PROJection FREQuencies",
    requires="[[optdriver]] == 66",
    added_in_version="before_v9",
    text=r"""
This variable sets the frequencies to be used to construct the basis in which
the Hamiltonian is projected to accelerate the solution of the Sternheimer
equations involved by the construction of the dielectric matrix at finite
frequencies. See Sec. VI of [[cite:Laflamme2015]]. For most cases,
since the frequencies $\infty$ and 0.0 (if [[gwls_recycle]] > 0) are used at no
computational cost, [[gwls_n_proj_freq]] == 0 (which means no ADDITIONAL
frequency is to be used) is fine and no frequencies need to be picked up.
""",
),

Variable(
    abivarname="gwls_model_parameter",
    varset="gw",
    vartype="real",
    topics=['GWls_expert'],
    dimensions="scalar",
    defaultval=1.0,
    mnemonics="GWLS MODEL PARAMETER",
    characteristics=['[[ENERGY]]'],
    requires="[[optdriver]] == 66",
    added_in_version="before_v9",
    text=r"""
This is the width of the lorentzian, in Ha, used to model the frequency
dependence of the dielectric matrix in the GWLS calculation [see Eqs. (12-16)
and (34) of [[cite:Laflamme2015]]\]. More
precisely, this parameter is the value of $\alpha$ used in Eq. (34). This model
is then used to separate the integration over frequencies into a 'model' part
[second term of Eq. (12)] and an 'exact - model' part [first term of Eq. (12)].
Since the 'model' part can be integrated analytically [see Eqs. (15), (16) and
(34)], only the 'exact - model' part needs to be integrated numerically.

The only effect of this model is therefore to alleviate the numerical cost of
the integration over frequencies in the G0W0 calculation. The value of the
associated parameter has thus an impact on the convergence rate of the GWLS
calculation with respect to the number of frequencies of integration
([[gwls_npt_gauss_quad]]), but no impact on the converged result of the GWLS
calculation. Typically, the default ([[gwls_model_parameter]] == 1.0) is optimal.
""",
),

Variable(
    abivarname="gwls_n_proj_freq",
    varset="gw",
    vartype="integer",
    topics=['GWls_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GWLS Number of PROJection FREQuencies",
    requires="[[optdriver]] == 66",
    added_in_version="before_v9",
    text=r"""
This variable sets the number of frequencies, on top of $\infty$ and 0.0 (if
[[gwls_recycle]] > 0), to be used for the construction of the basis in which
the Hamiltonian is projected to accelerate the solution of the Sternheimer
equations involved in the construction of the dielectric matrix at finite
frequencies. See Sec. VI of [[cite:Laflamme2015]]. For most cases,
the default ([[gwls_n_proj_freq]] == 0) is fine.
""",
),

Variable(
    abivarname="gwls_npt_gauss_quad",
    varset="gw",
    vartype="integer",
    topics=['GWls_expert'],
    dimensions="scalar",
    defaultval=10,
    mnemonics="GWLS Number of PoinTs to use for the GAUSSian QUADrature",
    requires="[[optdriver]] == 66",
    added_in_version="before_v9",
    text=r"""
This variable defines the number of points used for the numerical integration
of the self-energy over frequencies in GWLS computations [see Eq. (12) of
[[cite:Laflamme2015]]\]. The default is fine for most cases.
""",
),

Variable(
    abivarname="gwls_nseeds",
    varset="gw",
    vartype="integer",
    topics=['GWls_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="GWLS Number of SEED vectorS",
    requires="[[optdriver]] == 66",
    added_in_version="before_v9",
    text=r"""
This variable sets the number of seed vectors to be used in the construction
of the Lanczos basis for the (static) dielectric matrix in a GWLS calculation.
See Sec. IV of [[cite:Laflamme2015]]. Only [[gwls_nseeds]] == 1 has
been tested for now and users should keep this value.
""",
),

Variable(
    abivarname="gwls_print_debug",
    varset="gw",
    vartype="integer",
    topics=['GWls_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GWLS PRINT level for DEBUGging",
    requires="[[optdriver]] == 66",
    added_in_version="before_v9",
    text=r"""
Influences the level of verbosity for debugging purposes in a GWLS
calculation. Users should keep its value at the default.
""",
),

Variable(
    abivarname="gwls_recycle",
    varset="gw",
    vartype="integer",
    topics=['GWls_expert'],
    dimensions="scalar",
    defaultval=2,
    mnemonics="GWLS RECYCLE",
    requires="[[optdriver]] == 66",
    added_in_version="before_v9",
    text=r"""
This variable let the user choose if and how he wants to recycle the solutions
of the Sternheimer equations involved in the construction of the static dielectric matrix.

  * [[gwls_recycle]] = 0: No recycling of the solutions.
  * [[gwls_recycle]] = 1: Recycle the solutions. To do so, store them in RAM.
  * [[gwls_recycle]] = 2: Recycle the solutions. To do so, store them on disk.

If the user choose to recycle the solutions, they are used to construct the
basis in which the Hamiltonian is projected for the solution of the
Sternheimer equations involved by the calculation of the dielectric matrix at
finite frequencies. The other solutions used will be those at $\omega \to \infty$
(always used) and those at $\omega=$[[gwls_list_proj_freq]]. For more
information of the basis constructed, see Sec. IV of [[cite:Laflamme2015]].

It is important to note that the solutions rapidly take much space to store.
Therefore, it is often not possible to store them in RAM in production
calculations, yet still desirable to retain them. This is when it becomes
interesting to store them on disk. It is particularly efficient to choose the
path of the file to be on disk space local to the processor in large MPI
calculations, since each processor need only his own solutions in the
construction of the basis.
""",
),

Variable(
    abivarname="gwls_stern_kmax",
    varset="gw",
    vartype="integer",
    topics=['GWls_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="GWLS Kmax",
    requires="[[optdriver]] == 66",
    added_in_version="before_v9",
    text=r"""
This variable sets the dimension of the dielectric matrix used in a GWLS
calculation [see Sec. IV of [[cite:Laflamme2015]]\]. Typically
converged at a value of a few hundreds to a few thousands for a convergence
criterion of 50 meV on the eigenenergies.
""",
),

Variable(
    abivarname="gwmem",
    varset="gw",
    vartype="integer",
    topics=['Susceptibility_expert', 'SelfEnergy_expert'],
    dimensions="scalar",
    defaultval=11,
    mnemonics="GW MEMory",
    requires="[[optdriver]] in [3,4]",
    added_in_version="before_v9",
    text=r"""
[[gwmem]] governs the memory strategy during a screening and/or a sigma run.

  * [[gwmem]] = 1x, the screening matrix are read for all q-vectors and stored in the memory.
  * [[gwmem]] = 0x, the screening matrix are read just a q-vector after another.

  * [[gwmem]] = x1, the real-space wavefunctions are stored in the memory.
  * [[gwmem]] = x0, the real-space wavefunctions are not stored, but rather recalculated on-fly each abinit needs them using FFTs.

The default is [[gwmem]] = 11, which is the fastest, but also the most memory
consuming. When experiencing memory shortage, one should try [[gwmem]] = 0.
The first digit is only meaningful when performing sigma calculations.
""",
),

Variable(
    abivarname="gwpara",
    varset="paral",
    vartype="integer",
    topics=['parallelism_useful', 'GW_basic', 'Susceptibility_basic', 'SelfEnergy_basic'],
    dimensions="scalar",
    defaultval=2,
    mnemonics="GW PARAllelization level",
    commentdefault="The default value has been changed in v8. From 1 to 2",
    requires="[[optdriver]] in [3,4]",
    added_in_version="before_v9",
    text=r"""
gwpara is used to choose between the two different parallelization levels
available in the GW code. The available options are:

  * 1 --> parallelisation on k points.
  * 2 --> parallelisation on bands.

In the present status of the code, only the parallelization over bands
([[gwpara]] = 2) allows one to reduce the memory allocated by each processor.
Using [[gwpara]] = 1, indeed, requires the same amount of memory as a sequential
run, irrespectively of the number of CPUs used.
""",
),

Variable(
    abivarname="gwrpacorr",
    varset="gw",
    vartype="integer",
    topics=['RPACorrEn_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GW RPA CORRelation energy",
    requires="[[optdriver]] == 3 and [[gwcalctyp]] in [1,11,21]",
    added_in_version="before_v9",
    text=r"""
[[gwrpacorr]] governs the calculation of the RPA correlation energy.

  * [[gwrpacorr]] = 0, no RPA correlation energy is calculated.
  * [[gwrpacorr]] = 1, the RPA correlation energy is calculated using an exact integration
    over the coupling constant: it requires one diagonalization of the polarizability matrix.
  * [[gwrpacorr]] = $n$ > 1, the RPA correlation energy is calculated using $n$ values
    for the coupling constant: it requires $n$ inversions of the polarizability matrix.
""",
),

Variable(
    abivarname="hmctt",
    varset="rlx",
    vartype="integer",
    topics=['MolecularDynamics_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Hybrid Monte Carlo Trial Trajectory",
    requires="[[ionmov]] == 25",
    added_in_version="before_v9",
    text=r"""
Number of steps per MC trial trajectory, for the Hybrid Monte Carlo algorithm [[ionmov]]=25.
""",
),

Variable(
    abivarname="hmcsst",
    varset="rlx",
    vartype="integer",
    topics=['MolecularDynamics_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Hybrid Monte Carlo Strain Step Trajectory",
    requires="[[ionmov]] == 25",
    added_in_version="before_v9",
    text=r"""
Number of strain teps per MC trial trajectory, for the Hybrid Monte Carlo algorithm [[ionmov]]=25.
""",
),

Variable(
    abivarname="hyb_mixing",
    varset="gstate",
    vartype="real",
    topics=['Hybrids_expert'],
    dimensions="scalar",
    defaultval="-999.0",
    mnemonics="HYBrid MIXING coefficient for unscreened fock operator",
    commentdefault="With the default, [[hyb_mixing]] is initialized from [[ixc]].",
    requires="[[usefock]] > 0",
    added_in_version="before_v9",
    text=r"""
Mixing coefficient for the unscreened Fock operator in case of hybrid
functionals. Hartree-Fock corresponds to 1.0, PBE0 to 0.25.

ABINIT knows the correct value from [[ixc]]. Experts might nevertheless tune this mixing coefficient.
""",
),

Variable(
    abivarname="hyb_mixing_sr",
    varset="gstate",
    vartype="real",
    topics=['Hybrids_expert'],
    dimensions="scalar",
    defaultval="-999.0",
    mnemonics="HYBrid MIXING coefficient for Short-Range screened fock operator",
    commentdefault="With the default, [[hyb_mixing_sr]] is initialized from [[ixc]].",
    requires="[[usefock]] > 0",
    added_in_version="before_v9",
    text=r"""
Mixing coefficient for the screened Fock operator in case of hybrid
functionals. HSE has 0.25, B3LYP has 0.2.

ABINIT knows the correct value from [[ixc]]. Experts might nevertheless tune
this mixing coefficient.
""",
),

Variable(
    abivarname="hyb_range_dft",
    varset="gstate",
    vartype="real",
    topics=['Hybrids_expert'],
    dimensions="scalar",
    defaultval="-999.0 or [[hyb_range_fock]] if it is defined by the user",
    mnemonics="HYBrid RANGE for the DFT leftover from the screened fock operator",
    commentdefault="With the default=-999.0, [[hyb_range_dft]] is initialized from [[ixc]].",
    requires="[[usefock]] > 0",
    added_in_version="before_v9",
    text=r"""
Range of the DFT leftover from the screened Fock operator in case of hybrid
functionals (actually, coefficient of the distance appearing in the erf
function, thus it has the dimension of an inverse distance).

As described in the LibXC sources (and copied in the ABINIT doc, see
[[ixc]] = -428), there is a mess due to an error in the original publication.
ABINIT knows the LibXC value from [[ixc]], that might not agree with the
definitions from other codes. Usually, [[hyb_range_dft]] is the same as
[[hyb_range_fock]], see the latter for the different values. However, there is
a noticeable exception, the HSE03 from the original paper (not the HSE03 from VASP),
for which [[hyb_range_dft]] = 0.188988 while [[hyb_range_fock]] = 0.106066.
""",
),

Variable(
    abivarname="hyb_range_fock",
    varset="gstate",
    vartype="real",
    topics=['Hybrids_expert'],
    dimensions="scalar",
    defaultval="-999.0 or [[hyb_range_dft]] if it is defined by the user",
    mnemonics="HYBrid RANGE for the screened FOCK operator",
    commentdefault="With the default=-999.0, [[hyb_range_fock]] is initialized from [[ixc]].",
    requires="[[usefock]] > 0",
    added_in_version="before_v9",
    text=r"""
Range of the screened Fock operator in case of hybrid functionals (actually,
coefficient of the distance appearing in the erf function, thus it has the
dimension of an inverse distance).

As described in the LibXC sources (and copied in the ABINIT doc, see
[[ixc]] = -428), there is a mess due to an error in the original publication.
ABINIT knows the LibXC value from [[ixc]], that might not agree with the
definitions from other codes. Usually, [[hyb_range_dft]] is the same as
[[hyb_range_fock]], with one exception explained in [[hyb_range_dft]].
The HSE06 value from LibCX is 0.11, the one of Quantum Espresso is 0.106, the one of
VASP is 0.105835 (=0.2 $\AA^{-1}$).
The HSE03 value from LibCX is 0.106066 ($=0.15/\sqrt{2})$), the one of VASP is
0.1587531 (=0.3 $\AA^{-1}$).
""",
),

Variable(
    abivarname="iatcon",
    varset="rlx",
    vartype="integer",
    topics=['GeoConstraints_useful'],
    dimensions=['[[natcon]]', '[[nconeq]]'],
    defaultval=0,
    mnemonics="Indices of AToms in CONstraint equations",
    characteristics=['[[NO_MULTI]]', '[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Gives the indices of the atoms appearing in each of the [[nconeq]] independent
equations constraining the motion of atoms during structural optimization or
molecular dynamics (see [[nconeq]], [[natcon]], and [[wtatcon]]).
(Note: combined with [[wtatcon]] to give internal representation of the latter)
""",
),

Variable(
    abivarname="iatfix",
    varset="rlx",
    vartype="integer",
    topics=['GeoConstraints_basic'],
    dimensions=['[[natfix]]'],
    mnemonics="Indices of AToms that are FIXed",
    requires="[[natfix]] > 0",
    added_in_version="before_v9",
    text=r"""
Give the index (in the range 1 to [[natom]] ) of each atom which is to be held
fixed for structural optimization or molecular dynamics. The variable
[[iatfix]] lists those fixed in the three directions, while the variables
[[iatfixx]], [[iatfixy]], and [[iatfixz]], allow to fix some atoms along x, y
or z directions, or a combination of these.

WARNING: The implementation is inconsistent !! For [[ionmov]] ==1, the fixing
of directions was done in cartesian coordinates, while for the other values of
[[ionmov]], it was done in reduced coordinates. Sorry for this.

There is no harm in fixing one atom in the three directions using [[iatfix]],
then fixing it again in other directions by mentioning it in **iatfixx**,
**iatfixy** or **iatfixz**.
The internal representation of these input data is done by the mean of one
variable [[iatfix]](3,[[natom]]), defined for each direction and each atom,
being 0 if the atom is not fixed along the direction, and 1 if the atom is
fixed along the direction. When some atoms are fixed along 1 or 2 directions,
the use of symmetries is restricted to symmetry operations whose (3x3)
matrices [[symrel]] are diagonal.
If the atom manipulator is used, [[iatfix]] will be related to the
preprocessed set of atoms, generated by the atom manipulator. The user must
thus foresee the effect of this atom manipulator (see [[objarf]]).
""",
),

Variable(
    abivarname="iatfixx",
    varset="rlx",
    vartype="integer",
    topics=['GeoConstraints_basic'],
    dimensions=['[[natfixx]]'],
    mnemonics="Indices of AToms that are FIXed along the X direction",
    characteristics=['[[INPUT_ONLY]]'],
    requires="[[natfixx]] > 0",
    added_in_version="before_v9",
    text=r"""
Give the index (in the range 1 to [[natom]] ) of each atom which is to be held
fixed ALONG THE X direction for structural optimization or molecular dynamics.
The variable [[iatfix]] lists those fixed in the three directions, while the
variables [[iatfixx]], [[iatfixy]], and [[iatfixz]], allow to fix some atoms
along x, y or z directions, or a combination of these. See the variable
[[iatfix]] for more information.
""",
),

Variable(
    abivarname="iatfixy",
    varset="rlx",
    vartype="integer",
    topics=['GeoConstraints_basic'],
    dimensions=['[[natfixy]]'],
    mnemonics="Indices of AToms that are FIXed along the Y direction",
    characteristics=['[[INPUT_ONLY]]'],
    requires="[[natfixy]] > 0",
    added_in_version="before_v9",
    text=r"""
Give the index (in the range 1 to [[natom]] ) of each atom which is to be held
fixed ALONG THE Y direction for structural optimization or molecular dynamics.
The variable [[iatfix]] lists those fixed in the three directions, while the
variables [[iatfixx]], [[iatfixy]], and [[iatfixz]], allow to fix some atoms
along x, y or z directions, or a combination of these.
See the variable [[iatfix]] for more information.
""",
),

Variable(
    abivarname="iatfixz",
    varset="rlx",
    vartype="integer",
    topics=['GeoConstraints_basic'],
    dimensions=['[[natfixz]]'],
    mnemonics="Indices of AToms that are FIXed along the Z direction",
    characteristics=['[[INPUT_ONLY]]'],
    requires="[[natfixz]] > 0",
    added_in_version="before_v9",
    text=r"""
Give the index (in the range 1 to [[natom]] ) of each atom which is to be held
fixed ALONG THE Z direction for structural optimization or molecular dynamics.
The variable [[iatfix]] lists those fixed in the three directions, while the
variables [[iatfixx]], [[iatfixy]], and [[iatfixz]], allow to fix some atoms
along x, y or z directions, or a combination of these. See the variable
[[iatfix]] for more information.
""",
),

Variable(
    abivarname="iatsph",
    varset="gstate",
    vartype="integer",
    topics=['printing_prdos', 'ElecBandStructure_useful', 'ElecDOS_useful'],
    dimensions=['[[natsph]]'],
    defaultval=Range(start=1, stop='[[natsph]]'),
    mnemonics="Index for the ATomic SPHeres of the atom-projected density-of-states",
    requires="[[prtdos]] == 3 or [[pawfatbnd]] in [1,2]",
    added_in_version="before_v9",
    text=r"""
[[iatsph]] gives the number of the [[natsph]] atoms around which the sphere
for atom-projected density-of-states will be build, in the [[prtdos]] = 3 case.
The radius of these spheres is given by [[ratsph]].
If [[pawfatbnd]] = 1 or 2, it gives the number of the [[natsph]] atoms around
which atom-projected band structure will be built.
""",
),

Variable(
    abivarname="iboxcut",
    varset="paw",
    vartype="integer",
    topics=['TuningSpeed_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer governing the internal use of BOXCUT - not a very good choice of variable name",
    added_in_version="before_v9",
    text=r"""
Concern all summations in the reciprocal space and is allowed in PAW and norm-conserving.

  * if set to 0 all reciprocal space summations are done in a sphere contained in the FFT box.
  * if set to 1 all reciprocal space summations are done in the whole FFT box (useful for tests).
""",
),

Variable(
    abivarname="icoulomb",
    varset="gstate",
    vartype="integer",
    topics=['Coulomb_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Index for the COULOMB treatment",
    added_in_version="before_v9",
    text=r"""
Defines the type of computation used for Hartree potential, local part of
pseudo-potential and ion-ion interaction:

  * [[icoulomb]] = 0: usual reciprocal space computation, using $1/\GG^2$ for the Hartree potential and using Ewald correction.
  * [[icoulomb]] = 1: free boundary conditions are used when the Hartree potential is computed,
    real space expressions of pseudo-potentials are involved (restricted to GTH pseudo-potentials)
    and simple coulomb interaction gives the ion-ion energy.
""",
),

Variable(
    abivarname="icutcoul",
    varset="gw",
    vartype="integer",
    topics=['GWls_compulsory', 'Susceptibility_basic', 'Coulomb_useful', 'SelfEnergy_basic'],
    dimensions="scalar",
    defaultval=6,
    mnemonics="Integer that governs the CUT-off for COULomb interaction",
    requires="[[optdriver]] in [3,4]",
    added_in_version="before_v9",
    text=r"""
Many-body calculations for isolated systems present a slow convergence with
respect to the size of the supercell due to the long ranged Coulomb
interaction and the high degree of non-locality of the operators involved. A
similar issue also occurs in fully periodic systems due to the presence of the
integrable Coulomb singularity at $\mathbf{G}=0$ that hinders the convergence with
respect to the number of q-points used to sample the Brillouin zone. The
convergence can be accelerated by replacing the true bare Coulomb interaction
with other expressions.

[[icutcoul]] defines the particular expression to be used for the Coulomb term
in reciprocal space. The choice of [[icutcoul]] depends on the dimensionality
of the system. Possible values of [[icutcoul]] are from 0 to 6. The
corresponding influential variables are [[vcutgeo]] and [[rcut]].

  * 0 --> sphere (molecules but also 3D-crystals).
  * 1 --> cylinder (nanowires, nanotubes).
  * 2 --> surface.
  * 3 --> 3D crystal (no cut-off, integration in a spherical mini-Brillouin Zone, legacy value).
  * 4 --> ERF, long-range only Coulomb interaction.
  * 5 --> ERFC, short-range only Coulomb interaction (e.g. as used in the HSE functional).
  * 6 --> auxiliary function integration for 3D systems from [[cite:Carrier2007]].
  * 7 --> auxiliary function for 3D systems of Gygi and Baldereschi [[cite:Gygi1986]].
  * 14 --> Monte-Carlo integration in the mini-Brillouin zone for ERF, long-range only Coulomb interaction.
  * 15 --> Monte-Carlo integration in the mini-Brillouin zone for ERFC, short-range only Coulomb interaction.
  * 16 --> Monte-Carlo integration in the mini-Brillouin zone for Full Coulomb interaction.

Note that Spencer and Alavi showed that the
spherical cutoff can efficiently be used also for 3D systems [[cite:Spencer2008]].
In the latter case, use a negative value for the cutoff radius of the sphere ([[rcut]]<0),
which is automatically calculated so that the volume enclosed in the sphere is
equal to the volume of the solid.
""",
),

Variable(
    abivarname="ieig2rf",
    varset="dfpt",
    vartype="integer",
    topics=['TDepES_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer for second-order EIGenvalues from Response-Function",
    added_in_version="before_v9",
    text=r"""
If [[ieig2rf]] is greater than 0, the code will produce a file, named with the
suffix _EIGR2D, containing the second-order electronic eigenvalues for the perturbation.
These files are used in the calculation of the thermal
correction to the electronic eigenvalues.

  * If [[ieig2rf]] is set to 1, the second-order electronic eigenvalues will be
    calculated from the DFPT method (Sternheimer).

  * If [[ieig2rf]] is set to 2, the second-order electronic eigenvalues will be
    calculated from the Allen-Cardona method (sum over states).

  * If [[ieig2rf]] is set to 3, the second-order electronic eigenvalues will be
    calculated from the DFPT method (sum over states) but using a different part
    of the code. This is equivalent to [[ieig2rf]] = 1 [debuging].

  * If [[ieig2rf]] is set to 4, the second-order electronic eigenvalues will be
    calculated from the dynamical DFPT method (Sternheimer). The code will
    generate _EIGR2D.nc files that contain the electron-phonon matrix element
    squared on the space orthogonal to the active space. The code will also
    produce _FAN.nc files that contain the electron-phonon matrix elements
    squared.

  * If [[ieig2rf]] is set to 5, the second-order electronic eigenvalues will be
    calculated from the dynamical DFPT method (Sternheimer). The code will
    generate _EIGR2D.nc files that contain the electron-phonon matrix element
    square on the space orthogonal to the active space. The code will also produce
    _GKK.nc files that contain electron-phonon matrix elements. This option is
    preferable for large system to [[ieig2rf]] = 4 as the GKK files take much
    less disk space and memory (but run a little bit slower).

!!! note

    [[ieig2rf]] = 4 and 5 can only be used if Abinit is compiled with NETCDF support.


Related variables: [[bdeigrf]], [[elph2_imagden]], [[getgam_eig2nkq]], [[smdelta]]
""",
),

Variable(
    abivarname="imgmov",
    varset="rlx",
    vartype="integer",
    topics=['CrossingBarriers_useful', 'PIMD_compulsory', 'TransPath_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="IMaGe MOVEs",
    added_in_version="before_v9",
    text=r"""
Control the collective changes of images (see [[nimage]],[[npimage]],
[[dynimage]], [[ntimimage]], [[tolimg]], [[istatimg]], [[prtvolimg]]).
Similar to [[ionmov]] in spirit, although here, a population of self-consistent
calculations for possibly different (evolving) geometries is managed, while with
[[ionmov]], only self-consistent calculation for one (evolving) geometry is managed.
In this respect the maximal number of time step for image propagation is
[[ntimimage]], corresponding to the input variable [[ntime]] of the single
geometry case. Also, the stopping criterion is governed by [[tolimg]],
corresponding to the input variable [[toldfe]] of the single geometry case.
The stopping condition is crude: the image propagation is stopped when the
mean value (over dynamic images) of the absolute difference of total energy
(previous and current time step) is less than [[tolimg]].

Actually, there might be combinations of [[ionmov]] and [[imgmov]] in which
the two mechanisms are at work. Usually, however, only one mechanism will be
activated (so, usually, either [[ntimimage]] is bigger than one OR [[ntime]]
is bigger than one). In order for the user to acquire a mental representation
of the interplay between [[ionmov]] and [[imgmov]], here is a F90 pseudo-code
presenting the interplay between the different above-mentioned input
variables, as well as with the parallelism (see input variable [[npimage]]).

```fortran

    do itimimage=1,ntimimage
      do iimage=1,nimage
        (possibly, parallelisation over images)
        do itime=1,ntime
          Compute the forces and stresses for image(iimage)
          Examine whether the stopping criterion defined by tolmxf is fulfilled
          Predict the next geometry for image(iimage) using ionmov
        enddo
      enddo
      Examine whether the stopping criterion defined by tolimg is fulfilled
      Predict the next geometries for all images using imgmov
    enddo
```

  * = 0 --> simply **copy** images from previous itimimage step.
  * = 1 --> move images according to **Steepest Descent** following the (scaled) forces,
    the scaling factor being [[fxcartfactor]].
  * = 2 --> **String Method** for finding Minimal Energy Path (MEP) connecting two minima
    (see [[cite:Weinan2002]]), or even two configurations that are not local minima; the algorithm variant can be selected with the [[string_algo]] keyword
    (Simplified String Method by default). The solver for the Ordinary Differential Equation (ODE)
    can be selected with [[mep_solver]] (steepest-descent by default). See also [[mep_mxstep]] keyword.
  * = 3 --> (tentatively, not yet coded) **Metadynamics**.
  * = 4 --> (tentatively, not yet coded) **Genetic Algorithm**.
  * = 5 --> **Nudged Elastic Band (NEB)** for finding Minimal Energy Path (MEP) connecting two minima;
    the algorithm variant can be selected with the [[neb_algo]] keyword (NEB+improved tangent by default).
    The solver for the Ordinary Differential Equation (ODE) can be selected with [[mep_solver]] (steepest-descent by default).
    The spring constant connecting images along the path is defined by [[neb_spring]]. See also [[mep_mxstep]] keyword.
  * = 6 --> **Linear Combination of Constrained DFT Energies**. The images can have different electronic structure ([[occ]] can differ),
    and their total energies are combined linearly using the factors in [[mixesimgf]], giving the actual total energy of the ensemble
    of constrained DFT images. The geometry is the same for all images, forces and stresses are computed, and all usual
    algorithms for MD or geometry optimization are allowed, using [[ionmov]] (instead of [[imgmov]], this is the exception to the rule)
    and related variables.
  * = 9 or 13 --> **Path-Integral Molecular Dynamics** (see e.g. [[cite:Marx1996]]).
    Will use 9 for **Langevin thermostat** (associated friction coefficient given by [[vis]])
    and 13 for **Nose-Hoover thermostat chains** (associated input variables are the number of thermostats in the chains,
    [[nnos]], and the masses of these thermostats [[qmass]]). [[nimage]] is the Trotter number
    (no use of [[dynimage]]); possible transformations of coordinates are defined by [[pitransform]];
    Fictitious masses of the atoms (possibly different from the true masses given by [[amu]]) can be specified by [[pimass]].
    At present, it is only possible to perform calculations in the (N,V,T) ensemble ([[optcell]] = 0).

No meaning for RF calculations.
""",
),

Variable(
    abivarname="imgwfstor",
    varset="rlx",
    vartype="integer",
    topics=['CrossingBarriers_useful', 'PIMD_useful', 'TransPath_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="IMaGe WaveFunction STORage",
    requires="[[extrapwf]] == 0 and [[ntimimage]] > 0",
    added_in_version="before_v9",
    text=r"""
Govern the storage of wavefunctions at the level of the loop over images, see [[ntimimage]].
Possible values of [[imgwfstor]] are 0 or 1.
If [[imgwfstor]] is 1, the wavefunctions for each image are stored in a big array of
size [[nimage]] more than the storage needed for one set of wavefunctions.
When the specific computation (optimization/SCF cycle etc) for this image is started,
the past wavefunctions are used, to speed up the computation. If [[imgwfstor]]==0,
the wavefunctions are reinitialised, either at random or from the initial wavefunction file (so, without
any modification to take into account the computations at the previous value of itimimage.

If [[nimage]] is large, the increase of memory need can be problematic, unless the wavefunctions
are spread over many processors, which happens when [[paral_kgb]] == 1.
For some algorithms, e.g. when some geometry optimization
is performed, [[imgmov]]==2 or 5, the gain in speed of choosing [[imgwfstor]]=1 can be quite large, e.g. two to four.
For algorithms of the molecular dynamics type, [[imgmov]]==9 or 13, the expected gain is smaller.
Of course, with adequate memory resources, [[imgwfstor]]==1 should always be preferred.
""",
),

Variable(
    abivarname="inclvkb",
    varset="gw",
    vartype="integer",
    topics=['Susceptibility_expert', 'BSE_expert'],
    dimensions="scalar",
    defaultval=2,
    mnemonics="INCLude VKB",
    requires="[[optdriver]] in [3,99]",
    added_in_version="before_v9",
    text=r"""
Possible values of [[inclvkb]] are 0,1,2. If [[inclvkb]] is 1 or 2, the
commutator of the non-local part of the pseudopotential with the position
operator is correctly included in the q --> 0 contribution. This is
unfortunately time-consuming and in particular when the old algorithm
implemented by [[inclvkb]] = 1 is used ([[inclvkb]] = 2 is the recommended option). When
[[inclvkb]] is 0, this contribution is incorrectly omitted, but the computation is much faster.

The importance of this contribution depends on the number of k points. Turning
off [[inclvkb]] is to let to the choice of the user.

In general, the use of [[inclvkb]] = 0 is fine for GW calculations in
crystalline systems provided that the k-point sampling is sufficiently converged.

!!! important

    The use of [[inclvkb]] = 2 is strongly recommended for the calculation of optical properties.
""",
),

Variable(
    abivarname="intxc",
    varset="dev",
    vartype="integer",
    topics=['xc_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="INTerpolation for eXchange-Correlation",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
  * 0 --> do "usual" xc quadrature on fft grid
  * 1 --> do higher accuracy xc quadrature using fft grid and additional points at the centers of each cube
   (doubles number of grid points)--the high accuracy version is only valid for boxcut>=2. If boxcut < 2, the code stops.

For RF calculations only [[intxc]] = 0 is allowed yet. Moreover, the GS
preparation runs (giving the density file and zero-order wavefunctions) must be done with [[intxc]] = 0

Prior to ABINITv2.3, the choice [[intxc]] = 1 was favoured (it was the default),
but the continuation of the development of the code lead to prefer the default
[[intxc]] = 0. Indeed, the benefit of [[intxc]] = 1 is rather small, while making
it available for all cases is a non-negligible development effort.
Other targets are prioritary. You will notice that many automatic tests use
[[intxc]] = 1. Please, do not follow this historical choice for your production runs.
""",
),

Variable(
    abivarname="iomode",
    varset="dev",
    vartype="integer",
    topics=['parallelism_expert'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[MPI_IO]] and [[paral_kgb]] == 1': 1, 'defaultval': 0}),
    mnemonics="Input-Output MODE",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
This option selects the format used to produce "large" binary files such as the output wavefunction files,
the files with densities and potentials (DEN, POT) as well as the SCR file produced by the GW code.
Other "small" files such as the GSR.nc are always produced indipendently of the value of **iomode**.

Note that this variable mainly defines the format of the output files since Abinit is able to read
data from files independently of their format (either Fortran binary files or netcdf files).
The possible values are:

  * 0 --> Use standard Fortran IO (ok for sequential runs, not suitable for large parallel runs)
  * 1 --> Use MPI/IO routines (ok both for sequential and large parallel runs)
  * 3 --> Use NetCDF library to produce files according to the ETSF specification [[cite:Gonze2008]]
    (ok for sequential, requires netcdf4 + hdf5 + MPI-IO support for large parallel runs)

By default, Abinit produces Fortran files and uses the MPI-IO API when these operations
cannot be implemented in terms of simple Fortran write/read statements.
For example, [[paral_kgb]] = 1 uses the MPI-IO API to generate a Fortran binary file that can be read with
plain Fortran read statements.

There are cases, however, in which you would like to change the default behaviour.
For example, you may want to generate WFK or DEN files in netcdf
format because you need data in this format.
In this case, you have to use iomode == 3 in the input file to override the default behaviour.
Note, however, that you still need parallel IO capabilities enabled in the netcdf library if
you want to produce netcdf files in parallel with [[paral_kgb]] = 1
(i.e. netcdf4 + hdf5 + MPI-IO).
At present, the internal fallbacks provided by Abinit do not support netcdf4 so you have
to link against an external netcdf library that supports hdf5+MPI-IO
and is compatible with the mpif90 used to compile Abinit.
See ~abinit/doc/build/config-examples/ubu_intel_17.0_openmpi.ac for a typical configuration file.

!!! important

    The use of the ETSF_IO library [[cite:Caliste2008]] has been disabled, and replaced
    by direct NetCDF calls since the ETSF_IO library is not maintained anymore.
    The netcdf files, however, are still written following the ETSF-IO specifications [[cite:Gonze2008]]
    and extended with Abinit-specific quantities.
""",
),

Variable(
    abivarname="ionmov",
    varset="rlx",
    vartype="integer",
    topics=['MolecularDynamics_compulsory', 'GeoOpt_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="IONic MOVEs",
    added_in_version="before_v9",
    text=r"""
Choice of algorithm to control the displacements of ions, and eventually (see
[[optcell]]) changes of cell shape and size.

  * 0 --> Do not move ions;

  * 1 --> Move atoms using molecular dynamics with optional viscous damping (friction linearly proportional to velocity). The viscous damping is controlled by the parameter "[[vis]]". If actual undamped molecular dynamics is desired, set [[vis]] to 0. The implemented algorithm is the generalisation of the Numerov technique (6th order), but is NOT invariant upon time-reversal, so that the energy is not conserved. The value **ionmov** = 6 will usually be preferred, although the algorithm that is implemented is lower-order. The time step is governed by [[dtion]].
**Purpose:** Molecular dynamics (if [[vis]] = 0), Structural optimization (if
[[vis]] >0)
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:** Viscous parameter [[vis]], time step [[dtion]], index of atoms fixed [[iatfix]]

  * 2 --> Conduct structural optimization using the Broyden-Fletcher-Goldfarb-Shanno minimization (BFGS). This is much more efficient for structural optimization than viscous damping, when there are less than about 10 degrees of freedom to optimize. Another version of the BFGS is available with **ionmov** == 22, and is apparently more robust and efficient than **ionmov** == 2.
**Purpose:** Structural optimization
**Cell optimization:** Yes (if [[optcell]]/=0)
**Related variables:**

  * 3 --> Conduct structural optimization using the Broyden-Fletcher-Goldfarb-Shanno minimization (BFGS), modified to take into account the total energy as well as the gradients (as in usual BFGS).
See [[cite:Schlegel1982]]. Might be better than **ionmov** = 2 for few degrees of freedom (less than 3 or 4). Can be very
unstable - use with caution!
**Purpose:** Structural optimization
**Cell optimization:** Yes (if [[optcell]]/=0)
**Related variables:**

  * 4 --> Conjugate gradient algorithm for simultaneous optimization of potential and ionic degrees of freedom. It can be used with [[iscf]] = 2 and [[iscf]] =5 or 6 (WARNING: this is under development, and does not work very well in many cases).
**Purpose:** Structural optimization
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:**

  * 5 --> Simple relaxation of ionic positions according to (converged) forces. Equivalent to **ionmov** = 1 with zero masses, albeit the relaxation coefficient is not [[vis]], but [[iprcfc]].
**Purpose:** Structural optimization
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:**

  * 6 --> Molecular dynamics using the Verlet algorithm, see [[cite:Allen1987a]] p 81]. The only related parameter is the time step ([[dtion]]).
**Purpose:** Molecular dynamics
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:** time step [[dtion]], index of atoms fixed [[iatfix]]

  * 7 --> Quenched Molecular dynamics using the Verlet algorithm, and stopping each atom for which the scalar product of velocity and force is negative. The only related parameter is the time step ([[dtion]]). The goal is not to produce a realistic dynamics, but to go as fast as possible to the minimum. For this purpose, it is advised to set all the masses to the same value (for example, use the Carbon mass, i.e. set [[amu]] to 12 for all type of atoms).
**Purpose:** Structural optimization
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:** time step [[dtion]], index of atoms fixed [[iatfix]]

  * 8 --> Molecular dynamics with Nose-Hoover thermostat, using the Verlet algorithm.
**Purpose:** Molecular dynamics
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:** time step ([[dtion]]), Temperatures ([[mdtemp]]), and
thermostat mass ([[noseinert]]).

  * 9 --> Langevin molecular dynamics.
**Purpose:** Molecular dynamics
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:** time step ([[dtion]]), temperatures ([[mdtemp]]) and
friction coefficient ([[friction]]).

  * 10 --> Delocalized internal coordinates with BFGS simple
**Purpose:** Structural optimization
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:**

  * 11 --> Delocalized internal coordinates with BFGS using total energy
**Purpose:** Structural optimization
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:**

  * 12 --> Isokinetic ensemble molecular dynamics.
The equation of motion of the ions in contact with a thermostat are solved with the algorithm proposed in [[cite:Zhang1997]],
as worked out in [[cite:Minary2003]].
The conservation of the kinetic energy is obtained within machine precision, at each step.
As in [[cite:Evans1983]], when there is no fixing of atoms, the number of degrees of freedom in which the
microscopic kinetic energy is hosted is 3*natom-4. Indeed, the total kinetic energy is constrained, which accounts for
minus one degree of freedom (also mentioned in [[cite:Minary2003]]), but also there are three degrees of freedom
related to the total momentum in each direction, that cannot be counted as microscopic degrees of freedom, since the
total momentum is also preserved (but this is not mentioned in [[cite:Minary2003]]). When some atom is fixed in one or more direction,
e.g. using [[natfix]], [[natfixx]], [[natfixy]], or [[natfixz]], the number of degrees of freedom is decreased accordingly,
albeit taking into account that the total momentum is not preserved
anymore (e.g. fixing the position of one atom gives 3*natom-4, like in the non-fixed case).
**Purpose:** Molecular dynamics
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:** time step ([[dtion]]) and the first temperature in [[mdtemp]] in case the velocities [[vel]] are not initialized, or all initialized to zero.

  * 13 --> Isothermal/isenthalpic ensemble. The equation of motion of the ions in contact with a thermostat and a barostat are solved with the algorithm proposed in [[cite:Martyna1996]].
If optcell=1 or 2, the mass of the barostat ([[bmass]]) must be given in
addition.
**Purpose:** Molecular dynamics
**Cell optimization:** Yes (if [[optcell]]/=0)
**Related variables:** The time step ([[dtion]]), the temperatures
([[mdtemp]]), the number of thermostats ([[nnos]]), and the masses of
thermostats ([[qmass]]).

  * 14 --> Simple molecular dynamics with a symplectic algorithm proposed in [[cite:Blanes2002]]  (called SRKNa14] of the kind first published in [[cite:Yoshida1990]]This algorithm requires at least 14 evaluation of the forces (actually 15 are done within Abinit) per time step. At this cost it usually gives much better energy conservation than the verlet algorithm (**ionmov** 6) for a 30 times bigger value of [[dtion]]. Notice that the potential energy of the initial atomic configuration is never evaluated using this algorithm.
**Purpose:** Molecular dynamics
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:**

  * 15 --> Fast inertial relaxation engine (FIRE) algorithm proposed by
Erik Bitzek, Pekka Koskinen, Franz Gähler, Michael Moseler, and Peter Gumbsch in [[cite:Bitzek2006]].
According to the authors, the efficiency of this method is nearly the same as L-bfgs (**ionmov** = 22).
It is based on conventional molecular dynamics with additional velocity modifications and adaptive time steps.
The initial time step is set with [[dtion]]. Note that the physical meaning and unit of [[dtion]] are different from the default ones.
The purpose of this algorithm is relaxation, not molecular dynamics. [[dtion]] governs the ion position changes, but the cell parameter changes as well.
The positions are in reduced coordinates instead of in cartesian coordinates. The suggested first guess of dtion is 0.03.
**Purpose:** Relaxation
**Cell optimization:** Yes (if [[optcell]]/=0)
**Related variables:** The initial time step [[dtion]]

  * 20 --> Direct inversion of the iterative subspace. Given a starting point [[xred]] that is a vector of length 3*[[natom]] (reduced nuclei coordinates), and unit cell parameters ([[rprimd]]) this routine uses the DIIS (direct inversion of the iterative subspace) to minimize the gradient (forces) on atoms. The preconditioning used to compute errors from gradients is using an inverse hessian matrix obtained by a BFGS algorithm. This method is known to converge to the nearest point where gradients vanish. This is efficient to refine positions around a saddle point for instance.
**Purpose:** Structural optimization
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:** DIIS memory [[diismemory]]

  * 22 --> Conduct structural optimization using the Limited-memory Broyden-Fletcher-Goldfarb-Shanno minimization (L-BFGS) [[cite:Nocedal1980]]. The working routines were based on the original implementation of J. Nocedal available on netlib.org. This algorithm can be much better than the native implementation of BFGS in ABINIT (**ionmov** = 2) when one approaches convergence, perhaps because of better treatment of numerical details.
**Purpose:** Structural optimization
**Cell optimization:** Yes (if [[optcell]]/=0)
**Related variables:**

  * 23 --> Use of Learn on The Fly method (LOTF) for Molecular Dynamics. In the framework of isokinetic MD, the atomic forces and positions are computed by using LOTF interpolation. A SCF computation is performed only any [[lotf_nitex]] steps. The results of the SCF are used to compute the parameters of a short range classical potential (for the moment only the glue potential for gold is implemented). Then these parameters are continuously tuned to compute atomic trajectories. LOTF has to be enabled at configure time. If LOTF is not enabled and **ionmov** = 23, abinit will set automatically **ionmov** = 12.
The LOTF cycle is divided in the following steps:
a) Initialization (SFC at t=0) and computation of potential parameters.
b) Extrapolation of the atomic forces and positions for [[lotf_nitex]] time
step. To perform this extrapolation, the potential computed in a) is used
(Verlet algorithm).
c) SFC at t=[[lotf_nitex]]. Computation of the potential parameters.
d) LOTF interpolation, linear interpolation of the potential parameters and
computation of the atomic forces and positions between t=0 and t=lotf_nitex.
**Purpose:** Molecular Dynamics
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:** [[dtion]], [[lotf_classic]], [[lotf_nitex]],
[[lotf_nneigx]], [[lotf_version]].

  * 24 --> Simple constant energy molecular dynamics using the velocity Verlet symplectic algorithm (second order), see [[cite:Hairer2003]]. The only related parameter is the time step ([[dtion]]).
**Purpose:** Molecular dynamics
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:** time step [[dtion]]

  * 25 --> Hybrid Monte Carlo sampling of the ionic positions at fixed temperature and unit cell geometry (NVT ensemble). The underlying molecular dynamics corresponds to **ionmov** = 24. The related parameters are the time step ([[dtion]]) and thermostat temperature ([[mdtemp]]).
Within the HMC algorithm [[cite:Duane1987]], the trial states are generated via short $NVE$ trajectories (ten **ionmov** = 24 steps in current implementation).
 The initial momenta for each trial are randomly sampled from Boltzmann distribution, and the final trajectory state is either accepted or rejected based on the Metropolis criterion.
 Such strategy allows to simultaneously update all reduced coordinates, achieve higher acceptance ratio than classical Metropolis Monte Carlo and better sampling efficiency for shallow energy landscapes [[cite:Prokhorenko2018]].
**Purpose:** Monte Carlo sampling
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:** time step [[dtion]], thermostat temperature [[mdtemp]],

No meaning for RF calculations.
""",
),

Variable(
    abivarname="iprcel",
    varset="gstate",
    vartype="integer",
    topics=['SCFAlgorithms_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer for PReConditioning of ELectron response",
    added_in_version="before_v9",
    text=r"""
Used when [[iscf]] > 0, to define the SCF preconditioning scheme. Potential-
based preconditioning schemes for the SCF loop (electronic part) are still a
subject of active research. The present parameter (electronic part) describes
the way the change of potential is derived from the residual.
The possible values of [[iprcel]] correspond to:

  * 0 --> model dielectric function described by [[diemac]], [[dielng]] and [[diemix]].
  * larger or equal to 21 --> will compute the dielectric matrix according to [[diecut]], [[dielam]], [[diegap]]. This methodology is described in [[cite:Anglade2008]].
  * Between 21 and 29 --> for the first few steps uses the same as option 0 then compute RPA dielectric function, and use it as such.
  * Between 31 and 39 --> for the first few steps uses the same as option 0 then compute RPA dielectric function, and use it, with the mixing factor [[diemix]].
  * Between 41 and 49 --> compute the RPA dielectric matrix at the first step, and recompute it at a later step, and take into account the mixing factor [[diemix]].
  * Between 51 and 59 --> same as between 41 and 49, but compute the RPA dielectric matrix by another mean
  * Between 61 and 69 --> same as between 41 and 49, but compute the electronic dielectric matrix instead of the RPA one.
  * Between 71 and 78 --> STILL UNDER DEVELOPMENT -- NOT USABLE; Use the modified Kerker preconditioner with a real-space formulation (basic formulation is shown at [[dielng]]). The dielectric matrix is approximated thanks to [[diemac]] and [[dielng]]. Note that [[diemix]] is also used.
  * 79 --> STILL UNDER DEVELOPMENT -- NOT USABLE; same as previous but with an alternate algorithm.
  * 141 to 169 --> same as Between 41 and 69 (but, the dielectric matrix is also recomputed every iprcel modulo 10 step).

The computation of the dielectric matrix (for 0 [100]< [[iprcel]] < 70 [100])
is based on the **extrapolar** approximation, see [[cite:Anglade2008]]. This approximation can be tuned
with [[diecut]], [[dielam]], and [[diegap]]. Yet its accuracy mainly depends
on the number of conduction bands included in the system. Having 2 to 10 empty
bands in the calculation is usually enough (use [[nband]]).

NOTES:

  * The step at which the dielectric matrix is computed or recomputed is determined by modulo([[iprcel]],10). The recomputation happens just once in the calculation for [[iprcel]]  < 100.
  * For non-homogeneous relatively large cells [[iprcel]] = 45 will likely give a large improvement over [[iprcel]] = 0.
  * In case of PAW and [[iprcel]] > 0, see [[pawsushat]] input variable. By default, an approximation (which can be suppressed) is done for the computation of susceptibility matrix.
  * For extremely large inhomogeneous cells where computation of the full dielectric matrix takes too many weeks, 70 < [[iprcel]] < 80 is advised.
  * For [[nsppol]] = 2 or [[nspinor]] = 2 with metallic [[occopt]], only **mod(iprcel,100)** <50 is allowed.
  * No meaning for RF calculations yet.
  * The exchange term in the full dielectric matrix diverges for vanishing densities. Therefore the values of [[iprcel]] beyond 60 must not be used for cells containing vacuum, unless ones computes this matrix for every step ([[iprcel]] = 161).
""",
),

Variable(
    abivarname="iprcfc",
    varset="dev",
    vartype="integer",
    topics=['GeoOpt_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer for PReConditioner of Force Constants",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Used when [[iscf]] > 0, to define the SCF preconditioning scheme. Potential-
based preconditioning schemes for the SCF loop are still under development.
The present parameter (force constant part) describes the way a change of
force is derived from a change of atomic position.
Supported values:

  * 0 --> hessian is the identity matrix
  * 1 --> hessian is 0.5 times the identity matrix
  * 2 --> hessian is 0.25 times the identity matrix
  * -1 --> hessian is twice the identity matrix
  *... (simply corresponding power of 2 times the identity matrix)

No meaning for RF calculations.
""",
),

Variable(
    abivarname="iqpt",
    varset="gstate",
    vartype="integer",
    topics=['q-points_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Index for QPoinT generation",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Only used if [[nqpt]] = 1, and [[qptopt]] = 1 to 4.

Defines the index of the Q point to be selected in the list of q points
generated by [[ngqpt]], [[qptrlatt]], [[nshiftq]], and [[shiftq]].

If [[iqpt]] = 0, then the q point is Gamma (0 0 0).

The usual working mode is to define a series of values for [[iqpt]], starting
with [[iqpt]] = 0 or 1 (so through the definition of **iqpt:** ), and increasing
it by one for each dataset (thanks to **iqpt+** ).
""",
),

Variable(
    abivarname="irandom",
    varset="dev",
    vartype="integer",
    topics=['PIMD_expert'],
    dimensions="scalar",
    defaultval=3,
    mnemonics="Integer for the choice of the RANDOM number generator",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
For the time being, only used when [[imgmov]] = 9 (Langevin Path-Integral Molecular Dynamics).
[[irandom]] defines the random number generator.

Supported values:

  * 1 --> "uniformrandom", delivered with ABINIT package (initially comes from numerical recipes).
  * 2 --> intrinsic Fortran 90 random number generator.
  * 3 --> "ZBQ" non-deterministic random number generator by R. Chandler and P. Northrop.
  [Documentation](http://www.ucl.ac.uk/~ucakarc/work/software/randgen.txt) and
  [Source code](http://www.ucl.ac.uk/~ucakarc/work/software/randgen.f)

[[irandom]] = 3 is strongly advised when performing Molecular Dynamics restarts (avoids bias).
""",
),

Variable(
    abivarname="ird1den",
    varset="files",
    vartype="integer",
    topics=['nonlinear_useful'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[iscf]] < 0': 1, 'defaultval': 0}),
    mnemonics="Integer that governs the ReaDing of 1st-order DEN file",
    added_in_version="before_v9",
    text=r"""
If first order density is needed in single dataset mode (for example in
nonlinear optical response), use [[ird1den]] = 1 to read first-order densities
from _DENx files produced in other calculations. In multi-dataset mode use [[get1den]].

When [[iscf]] < 0, the reading of a DEN file is always enforced.

A non-zero value of **ird1den** is treated in the same way as other "ird" variables.
For further information about the *files file*, consult the [[help:abinit#files-file]].
""",
),

Variable(
    abivarname="ird1wf",
    varset="files",
    vartype="integer",
    topics=['DFPT_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer that governs the ReaDing of _1WF files",
    added_in_version="before_v9",
    text=r"""
Indicates eventual starting wavefunctions. As alternative, one can use the
input variables [[getwfk]], [[getwfq]], [[get1wf]] or [[getddk]].

Ground-state calculation:

  * only **irdwfk** and [[getwfk]] have a meaning
  * at most one of **irdwfk** or [[getwfk]] can be non-zero
  * if **irdwfk** and [[getwfk]] are both zero, initialize wavefunctions with random numbers for ground state calculation.
  * if **irdwfk** = 1: read ground state wavefunctions from a disk file appended with _WFK,
     produced in a previous ground state calculation.

Response-function calculation:

  * one and only one of **irdwfk** or [[getwfk]] MUST be non-zero
  * if **irdwfk** = 1: read ground state k -wavefunctions from a disk file appended with _WFK,
    produced in a previous ground state calculation.
  * only one of **irdwfq** or [[getwfq]] can be non-zero, if both of them are non-zero,
     use as k + q file the one defined by **irdwfk** and/or [[getwfk]]
  * if **irdwfq** = 1: read ground state k+q -wavefunctions from a disk file appended with _WFQ,
    produced in a previous ground state calculation.
  * at most one of [[ird1wf]] or [[get1wf]] can be non-zero
  * if both are zero, initialize first order wavefunctions to zeroes
  * if [[ird1wf]] = 1: read first-order wavefunctions from a disk file appended with _1WFx,
    produced in a previous response function calculation.
  * at most one of **irdddk** or [[getddk]] can be non-zero
  * one of them must be non-zero if an homogeneous electric field calculation is done
     (presently, a ddk calculation in the same dataset is not allowed)
  * if **irdddk** = 1: read first-order ddk wavefunctions from a disk file appended with _1WFx,
    produced in a previous response function calculation.

For further information about the *files file*, consult the [[help:abinit#files-file]].
""",
),

Variable(
    abivarname="irdbscoup",
    varset="files",
    vartype="integer",
    topics=['BSE_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer that governs the ReaDing of COUPling block",
    added_in_version="before_v9",
    text=r"""
Start the Bethe-Salpeter calculation from the BSC file containing the coupling block produced in a previous run.
""",
),

Variable(
    abivarname="irdbseig",
    varset="files",
    vartype="integer",
    topics=['BSE_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer that governs the ReaDing of BS_EIG file",
    added_in_version="before_v9",
    text=r"""
Start the Bethe-Salpeter calculation from the BS_EIG containing the exciton eigenvectors produced in a previous run.
""",
),

Variable(
    abivarname="irdbsreso",
    varset="files",
    vartype="integer",
    topics=['BSE_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer that governs the ReaDing of RESOnant block",
    added_in_version="before_v9",
    text=r"""
Start the Bethe-Salpeter calculation from the BSR file containing the resonant
block produced in a previous run.
""",
),

Variable(
    abivarname="irdddb",
    varset="files",
    vartype="integer",
    topics=['ElPhonInt_useful'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[iscf]] < 0': '1', 'defaultval': 0}),
    mnemonics="Integer that governs the ReaDing of DDB file",
    added_in_version="before_v9",
    text=r"""
This variable should be used when performing electron-phonon or temperature-
dependence calculations. The Born effective charge as well as the dielectric
tensor will be read from a previous DFPT calculations of the electric field at
q=Gamma. The use of this variable will trigger the cancellation of a residual
dipole that leads to an unphysical divergence of the GKK with vanishing
q-points. The use of this variable greatly improves the k-point convergence
speed as the density of the k-point grid required to obtain the fulfillment of
the charge neutrality sum rule is usually prohibitively large.

A non-zero value of [[irdddb]] is treated in the same way as other "ird" variables.
For further information about the *files file*, consult the [[help:abinit#files-file]].
Note also that, starting Abinit v9, one can also use [[getddb_path]] to specify the path of the DDB file directly.
""",

),

Variable(
    abivarname="irdddk",
    varset="files",
    vartype="integer",
    topics=['DFPT_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer that governs the ReaDing of DDK wavefunctions, in _1WF files",
    added_in_version="before_v9",
    text=r"""
Indicates eventual starting wavefunctions. As alternative, one can use the
input variables [[getwfk]], [[getwfq]], [[get1wf]] or [[getddk]].

Ground-state calculation:

  * only **irdwfk** and [[getwfk]] have a meaning
  * at most one of **irdwfk** or [[getwfk]] can be non-zero
  * if **irdwfk** and [[getwfk]] are both zero, initialize wavefunctions with random numbers for ground state calculation.
  * if **irdwfk** = 1: read ground state wavefunctions from a disk file appended with _WFK,
    produced in a previous ground state calculation

Response-function calculation:

  * one and only one of **irdwfk** or [[getwfk]] MUST be non-zero
  * if **irdwfk** = 1: read ground state k -wavefunctions from a disk file appended with _WFK,
    produced in a previous ground state calculation
  * only one of **irdwfq** or [[getwfq]] can be non-zero, if both of them are non-zero,
    use as k + q file the one defined by **irdwfk** and/or [[getwfk]]
  * if **irdwfq** = 1: read ground state k+q -wavefunctions from a disk file appended with _WFQ,
    produced in a previous ground state calculation
  * at most one of **ird1wf** or [[get1wf]] can be non-zero
  * if both are zero, initialize first order wavefunctions to zeroes
  * if **ird1wf** = 1: read first-order wavefunctions from a disk file appended with _1WFx,
    produced in a previous response function calculation
  * at most one of [[irdddk]] or [[getddk]] can be non-zero
  * one of them must be non-zero if an homogeneous electric field calculation is done (presently,
    a ddk calculation in the same dataset is not allowed)
  * if [[irdddk]] = 1: read first-order ddk wavefunctions from a disk file appended with _1WFx,
    produced in a previous response function calculation

For further information about the *files file*, consult the [[help:abinit#files-file]].
""",
),

Variable(
    abivarname="irdden",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[iscf]] < 0': '1', 'defaultval': 0}),
    mnemonics="Integer that governs the ReaDing of DEN file",
    added_in_version="before_v9",
    text=r"""
Start the ground-state calculation from the density file of a previous run.
When [[iscf]] < 0, the reading of a DEN file is always enforced.

A non-zero value of [[irdden]] is treated in the same way as other "ird" variables.
For further information about the *files file*, consult the [[help:abinit#files-file]].
""",
),

Variable(
    abivarname="irddvdb",
    varset="files",
    vartype="integer",
    topics=['ElPhonInt_useful'],
    dimensions="scalar",
    mnemonics="Integer that governs the ReaDing of DVDB file",
    added_in_version="before_v9",
    text=r"""
This variable can be used when performing electron-phonon calculations with [[optdriver]] = 7
to read an *input* DVDB file. See also [[getdvdb]]
""",
),

Variable(
    abivarname="irdefmas",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer to ReaD the EFfective MASses from...",
    added_in_version="before_v9",
    text=r"""
Eventually used when [[ndtset]] > 0 (multi-dataset mode).
Only relevant for [[optdriver]]=7 and [[eph_task]]=6.
If set to 1, take the data from a _EFMAS file as input. The latter must have been produced using [[prtefmas]] in another run.
""",
),

Variable(
    abivarname="irdhaydock",
    varset="files",
    vartype="integer",
    topics=['BSE_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer that governs the ReaDing of the HAYDOCK restart file",
    added_in_version="before_v9",
    text=r"""
Used to re-start the Haydock iterative technique from the HAYDR_SAVE file produced in a previous run.
""",
),

Variable(
    abivarname="irdqps",
    varset="files",
    vartype="integer",
    topics=['GW_useful', 'multidtset_useful', 'Susceptibility_useful', 'SelfEnergy_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer that governs the ReaDing of QuasiParticle Structure",
    added_in_version="before_v9",
    text=r"""
Relevant only when [[optdriver]] = 3 or 4. Indicate the file from which the
eigenvalues and possibly the wavefunctions must be obtained, in order to
achieve a self-consistent quasi-particle calculations. See also [[getqps]]
""",
),

Variable(
    abivarname="irdscr",
    varset="files",
    vartype="integer",
    topics=['GW_useful', 'multidtset_useful', 'SelfEnergy_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer that governs the ReaDing of the SCReening",
    added_in_version="before_v9",
    text=r"""
Relevant only when [[optdriver]] = 4. Indicate the file from which the
dielectric matrix must be obtained. As alternative, one can use the input variable [[getscr]].
When [[optdriver]] = 4, at least one of [[irdscr]] or [[getscr]] (alternatively,
[[irdsuscep]] or [[getsuscep]]) must be non-zero.

A non-zero value of [[irdscr]] is treated in the same way as other "ird" variables.
For further information about the *files file*, consult the [[help:abinit#files-file]].
""",
),

Variable(
    abivarname="irdsuscep",
    varset="files",
    vartype="integer",
    topics=['GW_useful', 'multidtset_useful', 'SelfEnergy_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer that governs the ReaDing of the SUSCEPtibility",
    added_in_version="before_v9",
    text=r"""
Relevant only when [[optdriver]] = 4. Indicate the file from which the
irreducible polarizability must be obtained.
As alternative, one can use the input variable [[getsuscep]].
When [[optdriver]] = 4, at least one of [[irdsuscep]] or [[getsuscep]]
(alternatively, [[irdscr]] or [[getscr]]) must be non-zero.

A non-zero value of [[irdsuscep]] is treated in the same way as other "ird" variables.
For further information about the *files file*, consult the [[help:abinit#files-file]].
""",
),

Variable(
    abivarname="irdvdw",
    varset="vdw",
    vartype="integer",
    topics=['vdw_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer that governs the ReaDing of _VDW files",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]] > 0, to read previously calculated vdW-DF variables.
Supported values:

  * 0: do not read vdW-DF variables
  * 1: read vdW-DF variables
""",
),

Variable(
    abivarname="irdwfk",
    varset="files",
    vartype="integer",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer that governs the ReaDing of _WFK files",
    added_in_version="before_v9",
    text=r"""
Indicates eventual starting wavefunctions. As alternative, one can use the
input variables [[getwfk]], [[getwfq]], [[get1wf]] or [[getddk]].

Ground-state calculation:

  * only [[irdwfk]] and [[getwfk]] have a meaning
  * at most one of [[irdwfk]] or [[getwfk]] can be non-zero
  * if [[irdwfk]] and [[getwfk]] are both zero, initialize wavefunctions with random numbers for ground state calculation.
  * if [[irdwfk]] = 1: read ground state wavefunctions from a disk file appended with _WFK,
    produced in a previous ground state calculation.

Response-function calculation:

  * one and only one of [[irdwfk]] or [[getwfk]] MUST be non-zero
  * if [[irdwfk]] = 1: read ground state k -wavefunctions from a disk file appended with _WFK,
    produced in a previous ground state calculation
  * only one of **irdwfq** or [[getwfq]] can be non-zero, if both of them are non-zero,
     use as k + q file the one defined by [[irdwfk]] and/or [[getwfk]]
  * if **irdwfq** = 1: read ground state k+q -wavefunctions from a disk file appended with _WFQ,
    produced in a previous ground state calculation
  * at most one of **ird1wf** or [[get1wf]] can be non-zero
  * if both are zero, initialize first order wavefunctions to 0's.
  * if **ird1wf** = 1: read first-order wavefunctions from a disk file appended with _1WFx,
    produced in a previous response function calculation
  * at most one of **irdddk** or [[getddk]] can be non-zero
  * one of them must be non-zero if an homogeneous electric field calculation is done
    (presently, a ddk calculation in the same dataset is not allowed)
  * if **irdddk** = 1: read first-order ddk wavefunctions from a disk file appended with _1WFx,
    produced in a previous response function calculation

For further information about the *files file*, consult the [[help:abinit#files-file]].
""",
),

Variable(
    abivarname="irdwfkfine",
    varset="dev",
    vartype="integer",
    topics=['multidtset_useful', 'TDepES_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer that governs the ReaDing of the grid _WFK file on the FINE grid",
    added_in_version="before_v9",
    text=r"""
Indicates eventual starting wavefunctions. As alternative, one can use the input variables [[getwfkfine]].

Ground-state calculation:

  * only [[irdwfkfine]] and [[getwfkfine]] have a meaning
  * at most one of [[irdwfkfine]] or [[getwfkfine]] can be non-zero
  * if [[irdwfkfine]] = 1: read ground state wavefunctions from a disk file appended with _WFK,
    produced in a previous ground state fine grid calculation

Response-function calculation:

  * one and only one of [[irdwfkfine]] or [[getwfkfine]] MUST be non-zero
  * if [[irdwfkfine]] = 1: read ground state k -wavefunctions from a disk file appended with _WFK,
    produced in a previous ground state calculation
  * Reading the fine grid wavefunction will trigger the k-points interpolation technique of the temperature dependent calculations.

Bethe-Salpeter calculation:

  * one and only one of [[irdwfkfine]] or [[getwfkfine]] MUST be non-zero
  * if [[irdwfkfine]] = 1: read ground state k -wavefunctions from a disk file appended with _WFK,
     produced in a previous ground state calculation
  * This variable or [[getwfkfine]] is mandatory when [[bs_interp_mode]] = 1

For further information about the *files file*, consult the [[help:abinit#files-file]].

**This variable is experimental. In development.**
""",
),

Variable(
    abivarname="irdwfq",
    varset="files",
    vartype="integer",
    topics=['DFPT_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer that governs the ReaDing of _WFQ files",
    added_in_version="before_v9",
    text=r"""
Indicates eventual starting wavefunctions. As alternative, one can use the
input variables [[getwfk]], [[getwfq]], [[get1wf]] or [[getddk]].

Ground-state calculation:

  * only **irdwfk** and [[getwfk]] have a meaning
  * at most one of **irdwfk** or [[getwfk]] can be non-zero
  * if **irdwfk** and [[getwfk]] are both zero, initialize wavefunctions with random numbers for ground state calculation.
  * if **irdwfk** = 1: read ground state wavefunctions from a disk file appended with _WFK,
    produced in a previous ground state calculation

Response-function calculation:

  * one and only one of **irdwfk** or [[getwfk]] MUST be non-zero
  * if **irdwfk** = 1: read ground state k -wavefunctions from a disk file appended with _WFK,
    produced in a previous ground state calculation
  * only one of [[irdwfq]] or [[getwfq]] can be non-zero, if both of them are non-zero,
    use as k + q file the one defined by **irdwfk** and/or [[getwfk]]
  * if [[irdwfq]] = 1: read ground state k+q -wavefunctions from a disk file appended with _WFQ,
    produced in a previous ground state calculation
  * at most one of **ird1wf** or [[get1wf]] can be non-zero
  * if both are zero, initialize first order wavefunctions to 0's.
  * if **ird1wf** = 1: read first-order wavefunctions from a disk file appended with _1WFx,
    produced in a previous response function calculation
  * at most one of **irdddk** or [[getddk]] can be non-zero
  * one of them must be non-zero if an homogeneous electric field calculation is done (presently,
    a ddk calculation in the same dataset is not allowed)
  * if **irdddk** = 1: read first-order ddk wavefunctions from a disk file appended with _1WFx,
    produced in a previous response function calculation

For further information about the *files file*, consult the [[help:abinit#files-file]].
""",
),

Variable(
    abivarname="iscf",
    varset="basic",
    vartype="integer",
    topics=['SCFAlgorithms_basic', 'TDDFT_compulsory', 'ElecBandStructure_basic'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[usepaw]] == 1': 17, '[[usewvl]] == 1': 0, 'defaultval': 7}),
    mnemonics="Integer for Self-Consistent-Field cycles",
    added_in_version="before_v9",
    text=r"""
Controls the self-consistency algorithm.

Positive values correspond to the usual choice for doing the usual ground state
(GS) calculations or for structural relaxations, where the potential has to be
determined self-consistently while negative values correspond to non-self-consistent calculations.

Different SCF algorithms (> 0) can be selected:

  * 0 --> SCF cycle, direct minimization scheme on the gradient of the wavefunctions. This algorithm is faster than diagonalisation and mixing but is working only for systems with a gap. It is implemented only on the wavelet basis set, when [[usewvl]] = 1.
  * 1 --> get the largest eigenvalue of the SCF cycle
([[DEVELOP]] option, used with [[irdwfk]] = 1 or [[irdwfq]] = 1)

  * 2 --> SCF cycle, simple mixing of the potential with the preconditioned potential residual (in the usual case, the latter is defined by [[diemix]], [[diemac]] and [[dielng]])
  * 3 --> SCF cycle, Anderson mixing of the potential
  * 4 --> SCF cycle, Anderson mixing of the potential based on the two previous iterations
  * 5 --> SCF cycle, CG based on the minim. of the energy with respect to the potential
  * 6 --> SCF cycle, CG based on the minim. of the energy with respect to the potential (alternate algo., [[DEVELOP]])
  * 7 --> SCF cycle, Pulay mixing of the potential based on the [[npulayit]] previous iterations
  * 12 --> SCF cycle, simple mixing of the density
  * 13 --> SCF cycle, Anderson mixing of the density
  * 14 --> SCF cycle, Anderson mixing of the density based on the two previous iterations
  * 15 --> SCF cycle, CG based on the minim. of the energy with respect to the density
  * 16 --> SCF cycle, CG based on the minim. of the energy with respect to the potential (alternate algo., [[DEVELOP]])
  * 17 --> SCF cycle, Pulay mixing of the density based on the [[npulayit]] previous iterations

!!! warning

    Other positive values, including zero ones, are not allowed.

Such algorithms for treating the "SCF iteration history" should be coupled
with accompanying algorithms for the SCF "preconditioning". See the input
variable [[iprcel]]. The default value [[iprcel]] = 0 is often a good choice,
but for inhomogeneous systems, you might gain a lot with [[iprcel]] = 45.

(Warning: if **iscf** > 10, at present (v4.6), the energy printed at each SCF
cycle is not variational - this should not affect the other properties, and at
convergence, all values are OK)

- In the norm-conserving case, the default option is **iscf** = 7, which is a
compromise between speed and reliability. The value **iscf** =  2 is safer but slower.
- In the PAW case, default option is **iscf** = 17. In PAW you have the
possibility to mix density/potential on the fine or coarse FFT grid (see [[pawmixdg]]).
- Note that a Pulay mixing (**iscf** = 7 or 17) with [[npulayit]] = 1 (resp. 2)
is equivalent to an Anderson mixing with **iscf** = 3 or 13 (resp. 4 or 14).
- Also note that:
* when mixing is done on potential (iscf < 10), total energy is computed by "direct" decomposition.
* when mixing is done on density (iscf >= 10), total energy is computed by "double counting" decomposition.
"Direct" and "double counting" decomposition of energy are equal when SCF
cycle is converged. Note that, when using GGA XC functionals, these
decompositions of energy can be slightly different due to imprecise
computation of density gradients on FFT grid (difference decreases as size of
FFT grid increases - see [[ecut]] for NC pseudopotentials, [[pawecutdg]] for PAW).

Other (**negative**) options:

  * -2 --> a non-self-consistent calculation is to be done;
    in this case an electron density rho(r) on a real space grid (produced in a previous calculation)
    will be read from a disk file (automatically if [[ndtset]] = 0, or according to the value of [[getden]] if [[ndtset]]/=0).
    The name of th density file must be given as indicated [[help:abinit#files-file|here]].
    **iscf** = -2 would be used for band structure calculations, to permit
    computation of the eigenvalues of occupied and unoccupied states at arbitrary
    k points in the fixed self consistent potential produced by some integration
    grid of k points. Due to this typical use, ABINIT insist that either
    [[prtvol]] > 2 or [[prteig]] does not vanish when there are more than 50 k points.
    To compute the eigenvalues (and wavefunctions) of unoccupied states in a
    separate (non-selfconsistent) run, the user should save the self-consistent
    rho(r) and then run **iscf** = -2 for the intended set of k-points and bands.
    To prepare a run with **iscf** = -2, a density file can be produced using the
    parameter [[prtden]] (see its description). When a self-consistent set of
    wavefunctions is already available, abinit can be used with [[nstep]] = 0 (see
    Test_v2/t47.in), and the adequate value of [[prtden]].

  * -3 --> like -2, but initialize [[occ]] and [[wtk]], directly or indirectly
    (using [[ngkpt]] or [[kptrlatt]]) depending on the value of [[occopt]].
    For GS, this option might be used to generate Density-of-states (thanks to
    [[prtdos]]), or to produce STM charge density map (thanks to [[prtstm]]).
    For RF, this option is needed to compute the response to ddk perturbation.

  * -1 --> like -2, but the non-self-consistent calculation is followed by the determination of excited states within TDDFT.
    This is only possible for [[nkpt]] = 1, with [[kpt]] = 0 0 0.
    Note that the oscillator strength needs to be defined with respect to an origin of coordinate,
    thanks to the input variable [[boxcenter]]. The maximal number of Kohn-Sham excitations to be used
    to build the excited state TDDFT matrix can be defined by [[td_mexcit]], or indirectly
    by the maximum Kohn-Sham excitation energy [[td_maxene]].
""",
),

Variable(
    abivarname="isecur",
    varset="dev",
    vartype="integer",
    topics=['SCFAlgorithms_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer for level of SECURity choice",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
In the presently used algorithms, there is a compromise between speed and
robustness, that can be tuned by using [[isecur]].
If [[isecur]] =0, an extrapolation of out-of-line data is allowed, and might
save one non-SCF calculation every two line minimisation when some stability
conditions are fulfilled (since there are 2 non-SCF calculations per line
minimisation, 1 out of 4 is saved)
Using [[isecur]] = 1 or higher integers will raise gradually the threshold to make extrapolation.
Using [[isecur]] = -2 will allow to save 2 non-SCF calculations every three line
minimisation, but this can make the algorithm unstable. Lower values of
[[isecur]] allows one for more (tentative) savings. In any case, there must be one
non-SCF computation per line minimisation.
No meaning for RF calculations yet.
""",
),

Variable(
    abivarname="istatimg",
    varset="rlx",
    vartype="integer",
    topics=['PIMD_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Integer governing the computation of STATic IMaGes",
    added_in_version="before_v9",
    text=r"""
This input variable is relevant when sets of images are activated (see [[imgmov]]).
Not all images might be required to evolve from one time step to the other
(see[[dynimage]]): these are static images.
If [[istatimg]] = 0, the total energy of static images is not computed (but
static images are used to make the dynamic images evolve).
This can be useful to save CPU time.
If [[istatimg]] = 1, the total energy of static images is computed.
""",
),

Variable(
    abivarname="istatr",
    varset="dev",
    vartype="integer",
    topics=['printing_prmisc'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integer for STATus file rate",
    characteristics=['[[DEVELOP]]', '[[NO_MULTI]]'],
    commentdefault="Values lower than 10 may not work on some machines.",
    added_in_version="before_v9",
    text=r"""
Govern the rate of output of the status file. This status file is written when
the number of the call to the status subroutine is equal to [[istatshft]]
modulo [[istatr]], so that it is written once every [[istatr]] call.
When [[istatr]] = 0, there is no writing of a status file (which is the default).
""",
),

Variable(
    abivarname="istatshft",
    varset="dev",
    vartype="integer",
    topics=['printing_prmisc'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Integer for STATus file SHiFT",
    characteristics=['[[DEVELOP]]', '[[NO_MULTI]]'],
    added_in_version="before_v9",
    text=r"""
Govern the rate of output of the status file. This status file is written when
the number of the call to the status subroutine is equal to '[[istatshft]]'
modulo ' **istatr** ', so that it is written once every ' **istatr** ' call.
There is also a writing for each of the 5 first calls, and the 10th call.
""",
),

Variable(
    abivarname="istwfk",
    varset="dev",
    vartype="integer",
    topics=['k-points_useful', 'TuningSpeed_basic'],
    dimensions=['[[nkpt]]'],
    defaultval=MultipleValue(number=None, value=0),
    mnemonics="Integer for choice of STorage of WaveFunction at each k point",
    commentdefault="For RF calculations, the Default is not used:  **istwfk** is forced to be 1 deep inside the code, for all k points. For spin-orbit calculations ([[nspinor]] = 2), **istwfk** is also forced to be 1, for all k points.",
    added_in_version="before_v9",
    text=r"""
Control the way the wavefunction for each k-point is stored inside ABINIT, in reciprocal space.
For the GS calculations, in the "cg" array containing the wavefunction
coefficients, there is for each k-point and each band, a segment
cg(1:2,1:npw). The 'full' number of plane wave is determined by [[ecut]].

However, if the k-point coordinates are build only from zeroes and halves (see
list below), the use of time-reversal symmetry (that connects coefficients)
has been implemented, in order to use real-to-complex FFTs (see [[fftalg]]),
and to treat explicitly only half of the number of plane waves (this being used as 'npw').
For the RF calculations, there is not only the "cg" array, but also the "cgq"
and "cg1" arrays. For the time-reversal symmetry to decrease the number of
plane waves of these arrays, the q vector MUST be (0 0 0). Then, for each k
point, the same rule as for the RF can be applied.
WARNING (991018): for the time being, the time-reversal symmetry cannot be
used in the RF calculations.

  * 1 --> do NOT take advantage of the time-reversal symmetry
  * 2 --> use time-reversal symmetry for k = (0, 0, 0,)
  * 3 --> use time-reversal symmetry for k = (1/2, 0, 0)
  * 4 --> use time-reversal symmetry for k = (0, 0, 1/2)
  * 5 --> use time-reversal symmetry for k = (1/2, 0, 1/2)
  * 6 --> use time-reversal symmetry for k = (0, 1/2, 0)
  * 7 --> use time-reversal symmetry for k = (1/2, 1/2, 0)
  * 8 --> use time-reversal symmetry for k = (0, 1/2, 1/2)
  * 9 --> use time-reversal symmetry for k = (1/2, 1/2, 1/2)
  * 0 --> (preprocessed) for each k point, choose automatically
          the appropriate time-reversal option when it is allowed, and chose **istwfk** = 1 for all the other k-points.
""",
),

Variable(
    abivarname="ixc",
    varset="basic",
    vartype="integer",
    topics=['xc_basic', 'Hybrids_compulsory', 'TDDFT_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Index of eXchange-Correlation functional",
    commentdefault="Default corresponds to Teter parametrization. However, if all the pseudopotentials have the same value of pspxc, the initial value of ixc will be that common value",
    added_in_version="before_v9",
    text=r"""
Controls the choice of exchange and correlation (xc). The list of XC
functionals is given below. Positive values are for ABINIT native library of
XC functionals, while negative values are for calling the much wider set of
functionals from the ETSF LibXC library (by M. Marques), available at the [
LibXC home page ](http://octopus-code.org/wiki/Libxc)
Note that the choice made here should preferably agree with the choice made in
generating the original pseudopotential, except for [[ixc]] = 0 (usually only
used for debugging). A warning is issued if this is not the case.
Unfortunately, pseudopotential (or PAW) generators for hybrid functionals and
mGGA are currently under development, so that one usually uses GGA or LDA
pseudopotentials instead. The error should be limited when GGA or LDA
pseudopotentials with semi-core states are used. Still this is a non-
controlled error. Moreover, the choices [[ixc]] = 1, 2, 3 and 7 are fits to the
same data, from Ceperley-Alder, and are rather similar, at least for spin-unpolarized systems.
The choice between the non-spin-polarized and spin-polarized case is governed
by the value of [[nsppol]] (see below).

**Native ABINIT XC functionals**

NOTE: in the implementation of the spin-dependence of these functionals, and
in order to avoid divergences in their derivatives, the interpolating function
between spin-unpolarized and fully-spin-polarized function has been slightly
modified, by including a zeta rescaled by 1.d0-1.d-6. This should affect total
energy at the level of 1.d-6Ha, and should have an even smaller effect on
differences of energies, or derivatives.
The value [[ixc]] = 10 is used internally: gives the difference between
[[ixc]] = 7 and [[ixc]] = 9, for use with an accurate RPA correlation energy.

  * 0 --> NO xc.

  * 1 --> LDA or LSD, Teter Pade parametrization (4/93, published in [[cite:Goedecker1996]], which reproduces Perdew-Wang 92 [[cite:Perdew1992a]] (which reproduces Ceperley-Alder [[cite:Ceperley1980]]!).
  * 2 --> LDA, Perdew-Zunger-Ceperley-Alder (no spin-polarization) [[cite:Perdew1981]]
  * 3 --> LDA, old Teter rational polynomial parametrization (4/91) fit to Ceperley-Alder data (no spin-polarization) [[cite:Ceperley1980]]
  * 4 --> LDA, Wigner functional (no spin-polarization)
  * 5 --> LDA, Hedin-Lundqvist functional (no spin-polarization) [[cite:Hedin1971]]
  * 6 --> LDA, "X-alpha" functional (no spin-polarization)
  * 7 --> LDA or LSD, Perdew-Wang 92 functional [[cite:Perdew1992a]]
  * 8 --> LDA or LSD, x-only part of the Perdew-Wang 92 functional [[cite:Perdew1992a]]
  * 9 --> LDA or LSD, x- and RPA correlation part of the Perdew-Wang 92 functional [[cite:Perdew1992a]]

  * 11 --> GGA, Perdew-Burke-Ernzerhof GGA functional [[cite:Perdew1996]]
  * 12 --> GGA, x-only part of Perdew-Burke-Ernzerhof GGA functional [[cite:Perdew1996]]
  * 13 --> GGA potential of van Leeuwen-Baerends [[cite:VanLeeuwen1994]], while for energy, Perdew-Wang 92 functional [[cite:Perdew1992a]]
  * 14 --> GGA, revPBE of [[cite:Zhang1998]]
  * 15 --> GGA, RPBE of [[cite:Hammer1999]]
  * 16 --> GGA, HTCH93 of [[cite:Hamprecht1998]]
  * 17 --> GGA, HTCH120 of [[cite:Boese2000]] - The usual HCTH functional.
  * 18 --> (NOT AVAILABLE: used internally for GGA BLYP pseudopotentials from [[cite:Krack2005]], available from the [ CP2K repository ](https://github.com/cp2k/cp2k/tree/master/data/POTENTIAL) \- use the LibXC instead, with [[ixc]] = -106131.
  * 19 --> (NOT AVAILABLE: used internally for GGA BP86 pseudopotentials from [[cite:Krack2005]], available from the [ CP2K repository ](https://github.com/cp2k/cp2k/tree/master/data/POTENTIAL) \- use the LibXC instead, with [[ixc]] = -106132.

  * 20 --> Fermi-Amaldi xc ( -1/N Hartree energy, where N is the number of electrons per cell; G=0 is not taken into account however),
                            for TDDFT tests. No spin-pol. Does not work for RF.
  * 21 --> same as 20, except that the xc-kernel is the LDA ([[ixc]] = 1) one, for TDDFT tests.
  * 22 --> same as 20, except that the xc-kernel is the Burke-Petersilka-Gross hybrid, for TDDFT tests.
  * 23 --> GGA of [[cite:Wu2006]].
  * 24 --> GGA, C09x exchange of [[cite:Cooper2010]].
  * 26 --> GGA, HTCH147 of [[cite:Boese2000]].
  * 27 --> GGA, HTCH407 of [[cite:Boese2001]].
  * 28 --> (NOT AVAILABLE: used internally for GGA OLYP pseudopotentials from [[cite:Krack2005]], available from the [ CP2K repository ](https://github.com/cp2k/cp2k/tree/master/data/POTENTIAL) \- use the LibXC instead, with [[ixc]] = -110131.

  * 40 --> Hartree-Fock
  * 41 --> PBE0, [[cite:Adamo1999]].
  * 42 --> PBE0-1/3, [[cite:Guido2013]].

**ETSF Lib XC functionals**

Note that you must compile ABINIT with the LibXC plug-in in order to be able
to access these functionals.
The LibXC functionals are accessed by **negative values** of [[ixc]]. The
LibXC contains functional forms for either exchange-only functionals,
correlation-only functionals, or combined exchange and correlation
functionals. Each of them is to be specified by a three-digit number. In case
of a combined exchange and correlation functional, only one such three-digit
number has to be specified as value of [[ixc]], with a minus sign (to indicate
that it comes from the LibXC). In the case of separate exchange functional
(let us represent its identifier by XXX) and correlation functional (let us
represent its identified by CCC), a six-digit number will have to be specified
for [[ixc]], by concatenation, be it XXXCCC or CCCXXX. As an example,
[[ixc]] = -020 gives access to the Teter93 LDA, while [[ixc]] = -101130 gives
access to the PBE GGA. In version 0.9 of LibXC (December 2008), there are 16
three-dimensional (S)LDA functionals (1 for X, 14 for C, 1 for combined XC),
and there are 41 three-dimensional GGA (23 for X, 8 for C, 10 for combined
XC). Note that for a meta-GGA, the kinetic energy density is needed.
This means having [[usekden]] = 1.

==(S)LDA functionals== (do not forget to add a minus sign, as discussed above)

  * 001 --> XC_LDA_X  [[cite:Dirac1930]], [[cite:Bloch1929]]
  * 002 --> XC_LDA_C_WIGNER  Wigner parametrization [[cite:Wigner1938]]
  * 003 --> XC_LDA_C_RPA  Random Phase Approximation [[cite:GellMann1957]]
  * 004 --> XC_LDA_C_HL  Hedin & Lundqvist [[cite:Hedin1971]]
  * 005 --> XC_LDA_C_GL  Gunnarson & Lundqvist [[cite:Gunnarsson1976]]
  * 006 --> XC_LDA_C_XALPHA  Slaters Xalpha
  * 007 --> XC_LDA_C_VWN  [[cite:Vosko1980]]
  * 008 --> XC_LDA_C_VWN_RPA  [[cite:Vosko1980]]
  * 009 --> XC_LDA_C_PZ  Perdew & Zunger [[cite:Perdew1981]]
  * 010 --> XC_LDA_C_PZ_MOD  Perdew & Zunger (Modified) [[cite:Perdew1981]] Modified to improve the matching between the low and high rs part
  * 011 --> XC_LDA_C_OB_PZ  Ortiz & Ballone (PZ) [[cite:Ortiz1994]] [[cite:Ortiz1997]] [[cite:Perdew1981]]
  * 012 --> XC_LDA_C_PW  Perdew & Wang [[cite:Perdew1992a]]
  * 013 --> XC_LDA_C_PW_MOD  Perdew & Wang (Modified) [[cite:Perdew1992a]]; Added extra digits to some constants as in the PBE routine.
  * 014 --> XC_LDA_C_OB_PW  Ortiz & Ballone (PW) [[cite:Ortiz1994]] [[cite:Ortiz1997]] [[cite:Perdew1992a]]
  * 017 --> XC_LDA_C_vBH  von Barth & Hedin [[cite:Barth1972]]
  * 020 --> XC_LDA_XC_TETER93  Teter 93 parametrization [[cite:Goedecker1996]]
  * 022 --> XC_LDA_C_ML1  Modified LSD (version 1) of Proynov and Salahub [[cite:Proynov1994]]
  * 023 --> XC_LDA_C_ML2  Modified LSD (version 2) of Proynov and Salahub [[cite:Proynov1994]]
  * 024 --> XC_LDA_C_GOMBAS  Gombas parametrization [[cite:Gombas1967]]
  * 025 --> XC_LDA_C_PW_RPA  Perdew & Wang fit of the RPA [[cite:Perdew1992a]]
  * 027 --> XC_LDA_C_RC04  Ragot-Cortona [[cite:Ragot2004]]
  * 028 --> XC_LDA_C_VWN_1  Vosko, Wilk, & Nussair (1) [[cite:Vosko1980]]
  * 029 --> XC_LDA_C_VWN_2  Vosko, Wilk, & Nussair (2) [[cite:Vosko1980]]
  * 030 --> XC_LDA_C_VWN_3  Vosko, Wilk, & Nussair (3) [[cite:Vosko1980]]
  * 031 --> XC_LDA_C_VWN_4  Vosko, Wilk, & Nussair (4) [[cite:Vosko1980]]

==GGA functionals== (do not forget to add a minus sign, as discussed above)

  * 084 --> XC_GGA_C_OP_XALPHA  one-parameter progressive functional (G96 version) [[cite:Tsuneda1999]]
  * 085 --> XC_GGA_C_OP_G96  one-parameter progressive functional (G96 version) [[cite:Tsuneda1999]]
  * 086 --> XC_GGA_C_OP_PBE  one-parameter progressive functional (PBE version) [[cite:Tsuneda1999]]
  * 087 --> XC_GGA_C_OP_B88  one-parameter progressive functional (B88 version) [[cite:Tsuneda1999]]
  * 088 --> XC_GGA_C_FT97  Filatov & Thiel correlation [[cite:Filatov1997a]]  [[cite:Filatov1997]]

!!! warning
    this functional is not tested. Use at your own risks.

  * 089 --> XC_GGA_C_SPBE  PBE correlation to be used with the SSB exchange [[cite:Swart2009]]
  * 090 --> XC_GGA_X_SSB_SW  Swart, Sola and Bickelhaupt correction to PBE [[cite:Swart2009a]]
  * 091 --> XC_GGA_X_SSB  [[cite:Swart2009]]

!!! warning
    This functional gives NaN on IBM (XG20130608).

  * 092 -->  XC_GGA_X_SSB_D  [[cite:Swart2009]]

!!! warning
    This functional gives NaN on IBM (XG20130608).

  * 093 -->  XC_GGA_XC_HCTH_407P  HCTH/407+ [[cite:Boese2003]]
  * 094 -->  XC_GGA_XC_HCTH_P76  HCTH p=7/6 [[cite:Menconi2001]]
  * 095 -->  XC_GGA_XC_HCTH_P14  HCTH p=1/4 [[cite:Menconi2001]]
  * 096 -->  XC_GGA_XC_B97_GGA1  Becke 97 GGA-1 [[cite:Cohen2000]]
  * 097 -->  XC_GGA_XC_HCTH_A  HCTH-A [[cite:Hamprecht1998]]
  * 098 -->  XC_GGA_X_BPCCAC  BPCCAC (GRAC for the energy) [[cite:Bremond2012]]
  * 099 -->  XC_GGA_C_REVTCA  Tognetti, Cortona, Adamo (revised) [[cite:Tognetti2008]]
  * 100 -->  XC_GGA_C_TCA  Tognetti, Cortona, Adamo [[cite:Tognetti2008a]]
  * 101 -->  XC_GGA_X_PBE  Perdew, Burke & Ernzerhof exchange [[cite:Perdew1996]]  [[cite:Perdew1997]]
  * 102 -->  XC_GGA_X_PBE_R  Perdew, Burke & Ernzerhof exchange (revised) [[cite:Zhang1998]]
  * 103 -->  XC_GGA_X_B86  Becke 86 Xalfa,beta,gamma [[cite:Becke1986]]
  * 104 -->  XC_GGA_X_HERMAN  Herman Xalphabeta GGA [[cite:Herman1969]]  [[cite:Herman2009]]
  * 105 -->  XC_GGA_X_B86_MGC  Becke 86 Xalfa,beta,gamma (with mod. grad. correction) [[cite:Becke1986]]  [[cite:Becke1986a]]
  * 106 -->  XC_GGA_X_B88  Becke 88 [[cite:Becke1988]]
  * 107 -->  XC_GGA_X_G96  Gill 96 [[cite:Gill1996]]
  * 108 -->  XC_GGA_X_PW86  Perdew & Wang 86 [[cite:Perdew1986a]]
  * 109 -->  XC_GGA_X_PW91  Perdew & Wang 91 [JP Perdew, in Proceedings of the 21st Annual International Symposium on the Electronic Structure of Solids, ed. by P Ziesche and H Eschrig (Akademie Verlag, Berlin, 1991), p. 11. ]  [[cite:Perdew1992]]  [[cite:Perdew1993]]
  * 110 -->  XC_GGA_X_OPTX  Handy & Cohen OPTX 01 [[cite:Handy2001]]
  * 111 -->  XC_GGA_X_DK87_R1  dePristo & Kress 87 (version R1) [[cite:DePristo1987]]
  * 112 -->  XC_GGA_X_DK87_R2  dePristo & Kress 87 (version R2) [[cite:DePristo1987]]
  * 113 -->  XC_GGA_X_LG93  Lacks & Gordon 93 [[cite:Lacks1993]]
  * 114 -->  XC_GGA_X_FT97_A  Filatov & Thiel 97 (version A) [[cite:Filatov1997a]]
  * 115 -->  XC_GGA_X_FT97_B  Filatov & Thiel 97 (version B) [[cite:Filatov1997a]]
  * 116 -->  XC_GGA_X_PBE_SOL  Perdew, Burke & Ernzerhof exchange (solids) [[cite:Perdew2008]]
  * 117 -->  XC_GGA_X_RPBE  Hammer, Hansen & Norskov (PBE-like) [[cite:Hammer1999]]
  * 118 -->  XC_GGA_X_WC  Wu & Cohen [[cite:Wu2006]]
  * 119 -->  XC_GGA_X_mPW91  Modified form of PW91 by Adamo & Barone [[cite:Adamo1998]]
  * 120 -->  XC_GGA_X_AM05  Armiento & Mattsson 05 exchange [[cite:Armiento2005]]  [[cite:Mattsson2008]]
  * 121 -->  XC_GGA_X_PBEA  Madsen (PBE-like) [[cite:Madsen2007]]
  * 122 -->  XC_GGA_X_MPBE  Adamo & Barone modification to PBE [[cite:Adamo2002]]
  * 123 -->  XC_GGA_X_XPBE  xPBE reparametrization by Xu & Goddard [[cite:Xu2004]]
  * 125 -->  XC_GGA_X_BAYESIAN  Bayesian best fit for the enhancement factor [[cite:Mortensen2005]]
  * 126 -->  XC_GGA_X_PBE_JSJR  PBE JSJR reparametrization by Pedroza, Silva & Capelle [[cite:Pedroza2009]]
  * 130 -->  XC_GGA_C_PBE  Perdew, Burke & Ernzerhof correlation [[cite:Perdew1996]]  [[cite:Perdew1997]]
  * 131 -->  XC_GGA_C_LYP  Lee, Yang & Parr [[cite:Lee1988]]  [[cite:Miehlich1989]]
  * 132 -->  XC_GGA_C_P86  Perdew 86 [[cite:Perdew1986]]
  * 133 -->  XC_GGA_C_PBE_SOL  Perdew, Burke & Ernzerhof correlation SOL [[cite:Perdew2008]]
  * 134 -->  XC_GGA_C_PW91  Perdew & Wang 91 [[cite:Perdew1992]]
  * 135 -->  XC_GGA_C_AM05  Armiento & Mattsson 05 correlation [[cite:Armiento2005]]  [[cite:Mattsson2008]]
  * 136 -->  XC_GGA_C_XPBE  xPBE reparametrization by Xu & Goddard [[cite:Xu2004]]
  * 137 -->  XC_GGA_C_LM  Langreth and Mehl correlation [[cite:Langreth1981]]
  * 138 -->  XC_GGA_C_PBE_JRGX  JRGX reparametrization by Pedroza, Silva & Capelle [[cite:Pedroza2009]]
  * 139 -->  XC_GGA_X_OPTB88_VDW  Becke 88 reoptimized to be used with vdW functional of Dion et al [[cite:Klimes2011]]
  * 140 -->  XC_GGA_X_PBEK1_VDW  PBE reparametrization for vdW [[cite:Klimes2011]]
  * 141 -->  XC_GGA_X_OPTPBE_VDW  PBE reparametrization for vdW [[cite:Klimes2011]]
  * 142 -->  XC_GGA_X_RGE2  Regularized PBE [[cite:Ruzsinszky2009]]
  * 143 -->  XC_GGA_C_RGE2  Regularized PBE [[cite:Ruzsinszky2009]]
  * 144 -->  XC_GGA_X_RPW86  refitted Perdew & Wang 86 [[cite:Murray2009]]
  * 145 -->  XC_GGA_X_KT1  Keal and Tozer version 1 [[cite:Keal2003]]
  * 146 -->  XC_GGA_XC_KT2 Keal and Tozer version 2 [[cite:Keal2003]]

!!! warning
    This functional gives NaN on IBM (XG20130608).

  * 147 -->  XC_GGA_C_WL  Wilson & Levy [[cite:Wilson1990]]
  * 148 -->  XC_GGA_C_WI  Wilson & Ivanov [[cite:Wilson1998]]
  * 149 -->  XC_GGA_X_MB88  Modified Becke 88 for proton transfer [[cite:Tognetti2009]]
  * 150 -->  XC_GGA_X_SOGGA  Second-order generalized gradient approximation [[cite:Zhao2008]]
  * 151 -->  XC_GGA_X_SOGGA11  Second-order generalized gradient approximation 2011 [[cite:Peverati2011]]
  * 152 -->  XC_GGA_C_SOGGA11  Second-order generalized gradient approximation 2011 [[cite:Peverati2011]]
  * 153 -->  XC_GGA_C_WI0  Wilson & Ivanov initial version [[cite:Wilson1998]]
  * 154 -->  XC_GGA_XC_TH1  Tozer and Handy v. 1 [[cite:Tozer1998]]

!!! warning
    This functional is not tested. Use at your own risks.

  * 155 -->  XC_GGA_XC_TH2  Tozer and Handy v. 2 [[cite:Tozer1998a]]
  * 156 -->  XC_GGA_XC_TH3  Tozer and Handy v. 3 [[cite:Handy1998]]
  * 157 -->  XC_GGA_XC_TH4  Tozer and Handy v. 4 [[cite:Handy1998]]
  * 158 -->  XC_GGA_X_C09X  C09x to be used with the VdW of Rutgers-Chalmers [[cite:Cooper2010]]
  * 159 -->  XC_GGA_C_SOGGA11_X  To be used with hyb_gga_x_SOGGA11-X [[cite:Peverati2011a]]
  * 161 -->  XC_GGA_XC_HCTH_93  HCTH functional fitted to 93 molecules [[cite:Hamprecht1998]]
  * 162 -->  XC_GGA_XC_HCTH_120  HCTH functional fitted to 120 molecules [[cite:Boese2000]]
  * 163 -->  XC_GGA_XC_HCTH_147  HCTH functional fitted to 147 molecules [[cite:Boese2000]]
  * 164 -->  XC_GGA_XC_HCTH_407  HCTH functional fitted to 407 molecules [[cite:Boese2001]]
  * 165 -->  XC_GGA_XC_EDF1  Empirical functionals from Adamson, Gill, and Pople [[cite:Adamson1998]]
  * 166 -->  XC_GGA_XC_XLYP  XLYP functional [[cite:Xu2004a]]
  * 167 -->  XC_GGA_XC_B97  Becke 97 [[cite:Becke1997]]
  * 168 -->  XC_GGA_XC_B97_1  Becke 97-1 [[cite:Hamprecht1998]]  [[cite:Becke1997]]
  * 169 -->  XC_GGA_XC_B97_2  Becke 97-2 [[cite:Becke1997]]
  * 170 -->  XC_GGA_XC_B97_D  Grimme functional to be used with C6 vdW term [[cite:Grimme2006]]
  * 171 -->  XC_GGA_XC_B97_K  Boese-Martin for Kinetics [[cite:Boese2004]]
  * 172 -->  XC_GGA_XC_B97_3  Becke 97-3 [[cite:Keal2005]]
  * 173 -->  XC_GGA_XC_PBE1W  Functionals fitted for water [[cite:Dahlke2005]]
  * 174 -->  XC_GGA_XC_MPWLYP1W  Functionals fitted for water [[cite:Dahlke2005]]
  * 175 -->  XC_GGA_XC_PBELYP1W  Functionals fitted for water [[cite:Dahlke2005]]
  * 176 -->  XC_GGA_XC_SB98_1a  Schmider-Becke 98 parameterization 1a [[cite:Schmider1998]]
  * 177 -->  XC_GGA_XC_SB98_1b  Schmider-Becke 98 parameterization 1b [[cite:Schmider1998]]
  * 178 -->  XC_GGA_XC_SB98_1c  Schmider-Becke 98 parameterization 1c [[cite:Schmider1998]]
  * 179 -->  XC_GGA_XC_SB98_2a  Schmider-Becke 98 parameterization 2a [[cite:Schmider1998]]
  * 180 -->  XC_GGA_XC_SB98_2b  Schmider-Becke 98 parameterization 2b [[cite:Schmider1998]]
  * 181 -->  XC_GGA_XC_SB98_2c  Schmider-Becke 98 parameterization 2c [[cite:Schmider1998]]
  * 183 -->  XC_GGA_X_OL2  Exchange form based on Ou-Yang and Levy v.2 [[cite:Fuentealba1995]]  [[cite:OuYang1991]]
  * 184 -->  XC_GGA_X_APBE  mu fixed from the semiclassical neutral atom [[cite:Constantin2011]]
  * 186 -->  XC_GGA_C_APBE  mu fixed from the semiclassical neutral atom [[cite:Constantin2011]]
  * 191 -->  XC_GGA_X_HTBS  Haas, Tran, Blaha, and Schwarz [[cite:Haas2011]]
  * 192 -->  XC_GGA_X_AIRY  Constantin et al based on the Airy gas [[cite:Constantin2009]]
  * 193 -->  XC_GGA_X_LAG  Local Airy Gas [[cite:Vitos2000]]
  * 194 -->  XC_GGA_XC_MOHLYP  Functional for organometallic chemistry [[cite:Schultz2005]]
  * 195 -->  XC_GGA_XC_MOHLYP2  Functional for barrier heights [[cite:Zheng2009]]
  * 196 -->  XC_GGA_XC_TH_FL  Tozer and Handy v. FL [[cite:Tozer1997]]
  * 197 -->  XC_GGA_XC_TH_FC  Tozer and Handy v. FC [[cite:Tozer1997]]
  * 198 -->  XC_GGA_XC_TH_FCFO  Tozer and Handy v. FCFO [[cite:Tozer1997]]
  * 199 -->  XC_GGA_XC_TH_FCO  Tozer and Handy v. FCO [[cite:Tozer1997]]
  * 200 -->  XC_GGA_C_OPTC  Optimized correlation functional of Cohen and Handy [[cite:Cohen2001]]
  * (for MetaGGA and Hybrid functionals, with indices in the 200-499 range, see the later sections)
  * 524 -->  XC_GGA_X_WPBEH  short-range version of the PBE [[cite:Heyd2003]]
  * 525 -->  XC_GGA_X_HJS_PBE  HJS screened exchange PBE version [[cite:Henderson2008]]
  * 526 -->  XC_GGA_X_HJS_PBE_SOL  HJS screened exchange PBE_SOL version [[cite:Henderson2008]]
  * 527 -->  XC_GGA_X_HJS_B88  HJS screened exchange B88 version [[cite:Henderson2008]]

!!! warning
    This functional is not tested. Use at your own risks.

  * 528 -->  XC_GGA_X_HJS_B97X  HJS screened exchange B97x version [[cite:Henderson2008]]
  * 529 -->  XC_GGA_X_ITYH  short-range recipe for exchange GGA functionals [[cite:Iikura2001]]

!!! warning
    This functional is not tested. Use at your own risks.


==MetaGGA functionals== (do not forget to add a minus sign, as discussed above).
See [[cite:Sun2011]] for the formulas.

  * 202 -->  XC_MGGA_X_TPSS  Tao, Perdew, Staroverov & Scuseria [[cite:Tao2003]]  [[cite:Perdew2004]]
  * 203 -->  XC_MGGA_X_M06L  Zhao, Truhlar exchange [[cite:Zhao2006]]  [[cite:Zhao2007]]
  * 204 -->  XC_MGGA_X_GVT4  GVT4 (X part of VSXC) from van Voorhis and Scuseria [[cite:Voorhis1998]]
  * 205 -->  XC_MGGA_X_TAU_HCTH  tau-HCTH from Boese and Handy [[cite:Boese2002]]
  * 207 -->  XC_MGGA_X_BJ06  Becke & Johnson correction to Becke-Roussel 89 [[cite:Becke2006]]

!!! warning
    This Vxc-only mGGA can only be used with a LDA correlation, typically Perdew-Wang 92 [[cite:Perdew1992a]].

  * 208 -->  XC_MGGA_X_TB09  Tran-blaha - correction to Becke & Johnson correction to Becke-Roussel 89 [[cite:Tran2009]]

!!! warning
    This Vxc-only mGGA can only be used with a LDA correlation, typically Perdew-Wang 92 [[cite:Perdew1992a]].

  * 209 -->  XC_MGGA_X_RPP09  Rasanen, Pittalis, and Proetto correction to Becke & Johnson [[cite:Rasanen2010]]

!!! warning
    This Vxc-only mGGA can only be used with a LDA correlation, typically Perdew-Wang 92 [[cite:Perdew1992a]].

  * 232 -->  XC_MGGA_C_VSXC  VSxc from Van Voorhis and Scuseria (correlation part) [[cite:Voorhis1998]]

==Hybrid functionals== (do not forget to add a minus sign, as discussed above).

  * 402 -->  XC_HYB_GGA_XC_B3LYP  The (in)famous B3LYP [[cite:Stephens1994]]
  * 406 -->  XC_HYB_GGA_XC_PBEH  PBEH (PBE0) [[cite:Adamo1999]]  [[cite:Ernzerhof1999]]
  * 427 -->  XC_HYB_GGA_XC_HSE03  The 2003 version of the screened hybrid HSE
                                  (this case corresponds to [[hyb_range_fock]]=$\omega^{HF} = 0.15/\sqrt{2}$
                                  and [[hyb_range_dft]]=$\omega^{PBE} = 0.15*(2.0)^{1/3}$ )
  * 428 -->  XC_HYB_GGA_XC_HSE06  The 2006 version of the screened hybrid HSE
                                  (this case corresponds to [[hyb_range_fock]]=[[hyb_range_dft]]=$\omega^{HF} = \omega^{PBE} = 0.11$)
                                  [[cite:Heyd2003]] [[cite:Heyd2006]] [[cite:Krukau2006]]

!!! warning
    (The following section is taken from the LibXC sources. In ABINIT, we stick to the LibXC choice.)

    Note that there is an enormous mess in the literature
    concerning the values of omega in HSE. This is due to an error in the original
    paper that stated that they had used $\omega=0.15$. This was in fact not true,
    and the real value used was $\omega^{HF} = 0.15 / \sqrt{2} \sim 0.1061$
    and $\omega^{PBE} = 0.15 * (2.0)^{1/3} \sim 0.1890$.

    In 2006 Krukau et al [[cite:Krukau2006]] tried
    to clarify the situation, called HSE03 the above choice of parameters,
    and called HSE06 to the functional where $\omega^{HF}=\omega^{PBE}$. By testing
    several properties for atoms they reached the conclusion that the best value
    for $\omega=0.11$. Of course, codes are just as messy as the papers. In Quantum Espresso
    HSE06 has the value $\omega=0.106.$ VASP, on the other hand, uses for HSE03 the
    same value $\omega^{HF} = \omega^{PBE} = 0.3 (A^{-1}) \sim 0.1587$,
    and for HSE06 $\omega^{HF} = \omega^{PBE} = 0.2 (A^{-1}) \sim 0.1058$.

  * 456 -->  XC_HYB_GGA_XC_PBE0_13  PBE0-1/3 [[cite:Cortona2012]]
""",
),

Variable(
    abivarname="ixc_sigma",
    varset="gw",
    vartype="integer",
    topics=['xc_expert', 'Hybrids_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Index of eXchange-Correlation functional used for self-energy calculations (SIGMA)",
    commentdefault="Default corresponds to Teter parametrization.",
    requires="mod([[gwcalctyp]],10)==5",
    added_in_version="before_v9",
    text=r"""
When [[gwcalctyp]] == 5, 15 or 25, [[ixc_sigma]] gives the identifier of the
advanced functional (usually a hybrid) that is used perturbatively or self-
consistently to obtain the improved electronic structure.

The meaning of the values of [[ixc_sigma]] is the same as the ones of [[ixc]],
so we refer to the latter for the list of possible values.

This input variable is introduced because in such calculation with
[[gwcalctyp]] == 5, 15 or 25, there is an underlying primary exchange-
correlation functional, that was used to obtain the starting wavefunctions and
eigenenergies, whose identified is [[ixc]]. The definition of both [[ixc]] and
[[ixc_sigma]] allows one to bypass possible sources of confusion.

Note however that in the case where [[gwcalctyp]] == 5, 15 or 25, the values of
the input variables [[auxc_ixc]], [[hyb_mixing]], [[hyb_mixing_sr]],
[[hyb_range_fock]] and [[hyb_range_dft]] refers to the advanced functional,
and not the primary one. Also, [[icutcoul]] and [[rcut]] have precedence over
[[hyb_mixing]], [[hyb_mixing_sr]], [[hyb_range_fock]] and [[hyb_range_dft]] to
define the parameters of the hybrid functional.
""",
),

Variable(
    abivarname="ixcpositron",
    varset="gstate",
    vartype="integer",
    topics=['positron_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Integer for the eXchange-Correlation applied to the electron-POSITRON interaction",
    commentdefault="(Teter parameterization). However, if all the pseudopotentials have the same value of pspxc, the initial value of ixc will be that common value",
    requires="[[positron]]/=0",
    added_in_version="before_v9",
    text=r"""
Relevant only when [[positron]]/=0.
Define the type of electron-positron correlation that is used in case of a
electron-positron two-component DFT calculation.
Define also the analytical formula of the enhancement factor used to compute
the electron-positron annihilation rate:

Electron-positron correlation functional:

  * ixcpositron=1: LDA zero positron density limit parametrized by
Arponen & Pajanne and provided by Boronski & Nieminen [[cite:Arponen1979a]],[[cite:Boronski1986]]
  * ixcpositron=11: LDA zero positron density limit parametrized by
Arponen & Pajanne and fitted by Sterne & Kaiser [[cite:Arponen1979a]],[[cite:Sterne1991]]
  * ixcpositron=2: LDA electron-positron correlation provided by
Puska, Seitsonen, and Nieminen [[cite:Arponen1979a]],[[cite:Puska1995]]
  * ixcpositron=3: GGA zero positron density limit parametrized by
Arponen & Pajanne and provided by Boronski & Nieminen [[cite:Arponen1979a]],[[cite:Boronski1986]],[[cite:Barbiellini1995]]
  * ixcpositron=31: GGA zero positron density limit parametrized
by Arponen & Pajanne and fitted by Sterne & Kaiser [[cite:Arponen1979a]],[[cite:Sterne1991]],[[cite:Barbiellini1995]]

Annihilation rate enhancement factor:

  * ixcpositron=1: Boronski and Nieminen full modelisation and RPA limit [[cite:Arponen1979a]]
  * ixcpositron=11: Sterne and Kaiser [[cite:Boronski1986]]
  * ixcpositron=2: Puska, Seitsonen and Nieminen [[cite:Sterne1991]]
  * ixcpositron=3: Boronski and Nieminen full modelisation and RPA limit [[cite:Arponen1979a]], with GGA corrections
  * ixcpositron=31: Sterne and Kaiser [[cite:Boronski1986]], with GGA corrections
""",
),

Variable(
    abivarname="ixcrot",
    varset="dfpt",
    vartype="integer",
    topics=['DFPT_expert', 'xc_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Index of the XC ROTation method used to calculate first-order exchange-correlation potential in non-collinear DFPT calculations",
    added_in_version="before_v9",
    text=r"""
Method of calculation of the 1st order XC potential in non-collinear DFPT
calculations. The possible values 1,2 and 3 correspond to the following
methods:

  * If ixcrot=1, the spinor rotation matrix U at each FFT point is not calculated explicitly. Instead the needed expressions involving U are derived based on the general properties of the U matrix.
  * If ixcrot=2, U is computed explicitly
  * If ixcrot=3, the brute force evaluation of the 1st order XC potential as a functional derivative is used. Rotation matrices are not computed.

In theory, all methods give identical results. However, due to different
implementation approaches, the round-off errors can lead to slight differences
intermediate and final results obtained using methods 1,2 and 3. The choice of
the method can also affect the convergence. For more details, see [[cite:Ricci2019]] or [[cite:Gonze2020]].
WARNING: in [[cite:Ricci2019]], the meaning of [[ixcrot]]=2 or 3 is inverted with respect to the implementation. On the contrary,
the implemention and [[cite:Gonze2020]] agree. More explicitly, the method refered to as method 1' in [[cite:Ricci2019]]
is [[ixcrot]]=2, and method refereed to as method 2 in [[cite:Ricci2019]] is [[ixcrot]]=3.

!!! note
    For non-zero perturbation wavevector ([[qpt]]/=0), only the [[ixcrot]]=3 implementation is currently available.
    The code will stop with the default [[ixcrot]] value for non-zero perturbation wavevector. The user should then set [[ixcrot]]=3 and restart.
""",
),

Variable(
    abivarname="jdtset",
    varset="basic",
    vartype="integer",
    topics=['multidtset_basic'],
    dimensions=['[[ndtset]]'],
    defaultval=Range({'start': 1, 'stop': '[[ndtset]]'}),
    mnemonics="index -J- for DaTaSETs",
    characteristics=['[[NO_MULTI]]'],
    added_in_version="before_v9",
    text=r"""
Gives the dataset index of each of the datasets. This index will be used:

  * to determine which input variables are specific to each dataset, since the variable names for this dataset will be made from the bare variable name concatenated with this index, and only if such a composite variable name does not exist, the code will consider the bare variable name, or even, the Default;
  * to characterize output variable names, if their content differs from dataset to dataset;
  * to characterize output files ( root names appended with _DSx where 'x' is the dataset index ).

The allowed index values are between 1 and 9999.
An input variable name appended with 0 is not allowed.
When [[ndtset]] == 0, this array is not used, and moreover, no input variable
name appended with a digit is allowed. This array might be initialized thanks
to the use of the input variable [[udtset]]. In this case, [[jdtset]] cannot be used.
""",
),

Variable(
    abivarname="jellslab",
    varset="gstate",
    vartype="integer",
    topics=['Artificial_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="include a JELLium SLAB in the cell",
    added_in_version="before_v9",
    text=r"""
If set to 1, a slab of uniform positive background charge density, that is, a
jellium slab, is included in the calculation cell. A portion of the unit cell
is filled with such positive charge density distribution which is equal to a
bulk-mean value $n_{bulk}$ between two edges and zero in the vacuum region if present.
For the sake of convenience the unit cell is supposed to have the third
crystal primitive lattice vector orthogonal to the other ones so that the
portion of the cell filled by the jellium slab can be defined through its edges along z.
The bulk-mean positive charge density is fixed by the input variable
[[slabwsrad]], while the position of the slab edges along z is defined through
the input variables [[slabzbeg]] and [[slabzend]].
""",
),

Variable(
    abivarname="jfielddir",
    varset="ffield",
    vartype="integer",
    topics=['Berry_basic'],
    dimensions=[3],
    defaultval=MultipleValue(number=3, value=0),
    mnemonics="electric/displacement FIELD DIRection",
    requires="[[berryopt]] = 17",
    added_in_version="before_v9",
    text=r"""
When specifying mixed electric field boundary conditions ( [[berryopt]] = 17),
[[jfielddir]] controls whether reduced electric field ([[jfielddir]] = 1) or reduced
electric displacement field ([[jfielddir]] = 2) is chosen to be fixed, in each
of the three lattice directions (i.e., in the reduced, not the Cartesian,
frame). For example, [[jfielddir]] = (1 1 2) tells the code to use fixed $\bar{e}_1$
and $\bar{e}_2$ along the first two lattice directions and fixed $d_3$ along the third.
For the case of mixed electric field boundary conditions, [[red_efieldbar]]
and [[red_dfield]] are used to control $\bar{e}$ and $d$, respectively. For example,
for electric boundary conditions corresponding to a material in a parallel-
plate capacitor, if you want to control $d_3=d_0$, while fixing $\bar{e}_1=\bar{e}_1=0$,
then the input files should have [[berryopt]] = 17, [[jfielddir]] = (1 1 2),
[[red_efieldbar]] = (0.0 0.0 a), and [[red_dfield]] = ($b\ c\ d_0$). Here a, b, and c
are the starting values. They can be chosen in this way: do a single run for
fixed d calculation ([[red_dfield]] = 0,0,$d_0$), from the final results you will
have $\bar{e}_3$, which is a good guess for a. Then do another single run for fixed
ebar calculation ([[red_efieldbar]] = (0 0 0)), from the final results you will
have $d_1$,$d_2$, these are good guesses for b, c.
""",
),

Variable(
    abivarname="jpawu",
    varset="paw",
    vartype="real",
    topics=['DFT+U_compulsory'],
    dimensions=['[[ntypat]]'],
    defaultval=MultipleValue(number=None, value=0),
    mnemonics="value of J for PAW+U",
    characteristics=['[[ENERGY]]'],
    requires="[[usepaw]] == 1 and [[usepawu]] == 1",
    added_in_version="before_v9",
    text=r"""
Gives the value of the screened exchange interaction between correlated
electrons corresponding to [[lpawu]] for each species.
In the case where [[lpawu]] =-1, the value is not used.
""",
),

Variable(
    abivarname="kberry",
    varset="ffield",
    vartype="integer",
    topics=['Berry_basic'],
    dimensions=[3, '[[nberry]]'],
    defaultval=MultipleValue(number=None, value=0),
    mnemonics="K wavevectors for BERRY phase computation",
    requires="[[berryopt]] = 1, 2, or 3",
    added_in_version="before_v9",
    text=r"""
Used for values of [[berryopt]] = 1, 2, or 3.

This array defines, for each Berry phase calculation (the number of such
calculations is defined by [[nberry]]), the difference of wavevector between k
points for which the overlap matrix must be computed. The polarisation vector
will be projected on the direction of that wavevector, and the result of the
computation will be the magnitude of this projection. Doing more than one
wavevector, with different independent direction, allows one to find the full
polarisation vector. However, note that converged results need oriented grids,
denser along the difference wavevector than usual Monkhorst-Pack grids.

The difference of wavevector is computed in the coordinate system defined by
the k-points grid (see [[ngkpt]] and [[kptrlatt]]), so that the values of
[[kberry]] are integers. Of course, such a k point grid must exist, and all
the corresponding wavefunctions must be available, so that the computation is
allowed only when [[kptopt]] is equal to 3. In order to save computing time,
it is suggested to make a preliminary calculation of the wavefunctions on the
irreducible part of the grid, with [[kptopt]] equal to 1, and then use these
converged wavefunctions in the entire Brillouin zone, by reading them to
initialize the [[kptopt]] = 3 computation.
""",
),

Variable(
    abivarname="kpt",
    varset="basic",
    vartype="real",
    topics=['k-points_useful'],
    dimensions=[3, '[[nkpt]]'],
    defaultval=[0, 0, 0],
    mnemonics="K - PoinTs",
    commentdefault="Adequate for one molecule in a supercell",
    added_in_version="before_v9",
    text=r"""
Contains the k points in terms of reciprocal space primitive translations (NOT
in cartesian coordinates!).
Needed ONLY if [[kptopt]] = 0, otherwise deduced from other input variables.

It contains dimensionless numbers in terms of which the cartesian coordinates
would be:
k_cartesian = k1*G1+k2*G2+k3*G3
where  (k1,k2,k3)  represent the dimensionless "reduced coordinates" and  G1,
G2, G3  are the cartesian coordinates of the primitive translation vectors.
G1,G2,G3 are related to the choice of direct space primitive translation
vectors made in [[rprim]]. Note that an overall norm for the k points is
supplied by [[kptnrm]]. This allows one to avoid supplying many digits for the
k points to represent such points as (1,1,1)/3.
Note: one of the algorithms used to set up the sphere of G vectors for the
basis needs components of k-points in the range [-1,1], so the remapping is
easily done by adding or subtracting 1 from each component until it is in the
range [-1,1]. That is, given the k point normalization [[kptnrm]] described
below, each component must lie in [ -[[kptnrm]], [[kptnrm]] ].
Note: a global shift can be provided by [[qptn]]
Not read if [[kptopt]]/=0.
""",
),

Variable(
    abivarname="kptbounds",
    varset="gstate",
    vartype="real",
    topics=['k-points_useful', 'ElecBandStructure_basic'],
    dimensions=[3, 'abs([[kptopt]])+1)'],
    mnemonics="K PoinT BOUNDarieS",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
It is used to generate the circuit to be followed by the band structure, when
[[kptopt]] is negative (it is not read if [[kptopt]] is zero or positive).

There are abs([[kptopt]]) segments to be defined, each of which starting from
the end point of the preceeding one. Thus, the number of points to be input is
abs([[kptopt]])+1. They form a circuit starting at
[[kptbounds]](1:3,1)/[[kptnrm]] and ending at
[[kptbounds]](1:3,abs([[kptopt]])+1)/[[kptnrm]]. The number of divisions of
each segment can be defined either using the array [[ndivk]] or the variable
[[ndivsm]] that just defines the number of divisions for the smallest segment

As for [[kpt]], [[kptbounds]] is specified using the primitive vectors in
reciprocal space. If your Bravais lattice is simple, then it should be quite
easy to find the coordinates of the end points. On the other hand, for
centered, body-centered, face-centered, hexagonal, and rhombohedral Bravais
lattice, the conversion might be more difficult. See the description of
[[kpt]] for an explanation of how to convert data from the "conventional"
cartesian coordinates to the primitive vectors in the reciprocal space. In
order to help a bit, we list below a series of typical values, for the FCC,
BCC, hexagonal and rhombohedral Bravais lattices. Note: all the data below
are given in dimensionless units; they have to be rescaled by the actual
lengths defined by the [[acell]] values. However, [[kptbounds]] values can be
used as such, if the values of [[rprim]] given below are adopted.

A. **FCC lattice**

Suppose the primitive vectors in real space are given by

      rprim   0 1 1    1 0 1    1 1 0
or

      rprim   0 1/2 1/2    1/2 0 1/2    1/2 1/2 0

(these two possibilities only differ by a scaling factor, irrelevant for the
definition of the k points in the primitive vectors in reciprocal space).
Then, the reciprocal primitive vectors (in conventional cartesian coordinates) are

      (-1/2 1/2 1/2), (1/2 -1/2 1/2), (1/2 1/2 -1/2)

or

      (-1 1 1), (1 -1 1), (1 1 -1)

and, in both cases, the coordinates of several special points with respect to
primitive vectors in reciprocal space are

      X (0   1/2 1/2)   (conventional cartesian coordinate 1/2 0 0)
      X'(1/2 1/2 1  )   (conventional cartesian coordinate 1/2 1/2 0)  (an other instance of X, in another Brillouin zone)
      L (1/2 1/2 1/2)   (conventional cartesian coordinate  1/4 1/4 1/4)
      L'(1/2 0   0  )   (conventional cartesian coordinate -1/4 1/4 1/4) (an other instance of L, on another face of the BZ)
      W (1/4 1/2 3/4)   (conventional cartesian coordinate 1/2 1/4 0)
      U (1/4 5/8 5/8)   (conventional cartesian coordinate 1/2 1/8 1/8)
      K (3/8 3/8 3/4)   (conventional cartesian coordinate 3/8 3/8 0)

Note that K is actually equivalent to U, by spatial and translational
symmetry. So, if you want to specify a typical circuit, the following might do
the work: L-Gamma-X-W-K,U-L-W-X-K,U-Gamma with

      kptbounds  1/2 0 0  0 0 0  0 1/2 1/2  1/4 1/2 3/4  3/8 3/8 3/4  1/2 1/2 1/2  1/4 1/2 3/4  1/2 1/2 1  3/8 3/8 3/4  0 0 0

The lengths of segments (this information is useful to draw the band
structure, with the correct relative scale between special points) can be
found using the conventional cartesian coordinates:

$l$(L-Gamma)=$\sqrt{3}/4=$0.433... \n
$l$(Gamma-X)=$1/2$=0.5 \n
$l$(X-W)=$1/4$=0.25 \n
$l$(W-K)=$\sqrt{2}/8$=0.177... \n
$l$(K-L)=$\sqrt{6}/8$=0.306... \n
$l$(L-W)=$\sqrt{2}/4$=0.354... \n
$l$(W-X)=$1/4$=0.25 \n
$l$(X-K)=$\sqrt{2}/8$=0.177... \n
$l$(K-Gamma)=$3\sqrt{2}/8$=0.530... \n

B. **BCC lattice**

Suppose the primitive vectors in real space are given by

      rprim  -1 1 1    1 -1 1    1 1 -1

(as for the FCC lattice, there is a scale invariance). Then, the reciprocal
primitive vectors (in conventional cartesian coordinates) are (0 1/2 1/2),
(1/2 0 1/2), and (1/2 1/2 0) and the coordinates of several special points
with respect to primitive vectors in reciprocal space are

      H (-1/2 1/2 1/2)   (conventional cartesian coordinate 1/2 0 0)
      N ( 0   0   1/2)   (conventional cartesian coordinate 1/4 1/4 0)
      P ( 1/4 1/4 1/4)   (conventional cartesian coordinate 1/4 1/4 1/4)

So, if you want to specify a typical circuit, the following might do the work: Gamma-H-N-Gamma-P-N-P-H

      kptbounds  0 0 0  -1/2 1/2 1/2  0 0 1/2  0 0 0   1/4 1/4 1/4  0 0 1/2  1/4 1/4 1/4  -1/2 1/2 1/2

The lengths of segments (this information is useful to draw the band
structure, with the correct relative scale between special points) can be
found using the conventional cartesian coordinates:

$l$(Gamma-H)=$1/2$=0.5 \n
$l$(H-N)=$\sqrt{2}/4$=0.354... \n
$l$(N-Gamma)=$\sqrt{2}/4$=0.354... \n
$l$(Gamma-P)=$\sqrt{3}/4$=0.433... \n
$l$(P-N)=$1/4$=0.25 \n
$l$(N-P)=$1/4$=0.25 \n
$l$(P-H)=$\sqrt{3}/4$=0.433... \n

C. **Hexagonal lattices**

Suppose the primitive vectors in real space are given by

      rprim  1 0 0    -1/2 sqrt(0.75) 0    0 0 1

The coordinates of several special points with respect to primitive vectors in
reciprocal space are

      M (1/2 0 0) or (0 1/2 0) or (-1/2 1/2 0)
      L (1/2 0 1/2) or (0 1/2 1/2) or (-1/2 1/2 1/2)
      K (1/3 1/3 0) or (2/3 -1/3 0) or (-1/3 2/3 0)
      H (1/3 1/3 1/2) or (2/3 -1/3 1/2) or (-1/3 2/3 1/2)
      A (0 0 1/2)

So, if you want to specify a typical circuit, the following
might do the work: K-Gamma-M-K-H-A-L-H-L-M-Gamma-A

      kptbounds  1/3 1/3 0  0 0 0  1/2 0 0  1/3 1/3 0  1/3 1/3 1/2  0 0 1/2  1/2 0 1/2  1/3 1/3 1/2  1/2 0 1/2  1/2 0 0  0 0 0  0 0 1/2

In order to find the lengths of segments (this information is useful to draw
the band structure, with the correct relative scale between special points)
one needs to know the a and c lattice parameters. Also, in what follows, we
omit the 2$\pi$ factor sometimes present in the definition of the reciprocal
space vectors. The reciprocal vectors are $(1/a\: 1/(\sqrt{3}a)\: 0)$, $(0\: 2/(\sqrt{3}a)\: 0)$,
 $(0\: 0\: 1/c)$. The lengths of the above-mentioned segments can
be computed as:

$l$(K-Gamma)=$2/(3a)$=0.666.../a \n
$l$(Gamma-M)=$1/(\sqrt{3}a)$=0.577.../a \n
$l$(M-K)=$1/(3a)$=0.333.../a \n
$l$(K-H)=$1/(2c)$=0.5.../c \n
$l$(H-A)=$2/(3a)$=0.666.../a \n
$l$(A-L)=$1/(\sqrt{3}a)$=0.577.../a \n
$l$(L-H)=$1/(3a)$=0.333.../a \n
$l$(H-L)=$1/(3a)$=0.333.../a \n
$l$(L-M)=$1/(2c)$=0.5.../c \n
$l$(M-Gamma)=$1/(\sqrt{3}a)$=0.577.../a \n
$l$(Gamma-A)=$1/(2c)$=0.5.../c \n

D. **Rhombohedral lattices**

Rhombohedral lattices are characterised by two parameters, the length of the
primitive vectors, that we will denote a0, and the angle they form, $\gamma$.
These can be directly input of ABINIT, as [[acell]] and [[angdeg]]

This will generate the primitive vectors in real space, with

      acell a0 a0 a0    and      rprim  a 0 c    -a/2 a*sqrt(0.75) c    -a/2 -a*sqrt(0.75) c

with,

 * $a^2+c^2=1$,
 * $a^2=2/3(1-\cos(\gamma))$,
 * $c^2=1/3(1+2\cos(\gamma))$,
 * $(a/c)^2=2(1-\cos(\gamma))/(1+2\cos(\gamma))$, and also
 * $\cos(\gamma)=(1-(a/c)^2/2)/(1+(a/c)^2)$.

Alternatively, these values of [[rprim]]
might directly be the input of ABINIT (then, the balance of the scaling factor
might be adjusted between [[acell]] and [[rprim]]).

Unlike for the simple cubic, FCC, BCC, hexagonal (and some other) Bravais
lattice, the topology of the Brillouin zone will depend on the $\gamma$ (or $a/c$)
value. We give below information concerning the case when $\cos(\gamma)$ is
positive, that is, $(a/c)^2$ lower than 2.

The coordinates of several special points with respect to primitive vectors in
reciprocal space will not depend on the $a/c$ ratio, but some others will depend
on it. So, some care has to be exercised. Notations for the Brillouin Zone
special points are the same as in [[cite:Gonze1990]].

      L (1/2 0 0) or (0 1/2 0) or (0 0 1/2) (or with negative signs)
      T (1/2 1/2 1/2)
      X (1/2 1/2 0) or (1/2 0 1/2) or (0 1/2 1/2) (or with separate negative signs)
      W (5/6 - (a/c)^2/6, 1/2, 1/6 + (a/c)^2/6 ) = (1 0 -1)*(1-(a/c)^2/2)/3 + (1 1 1)/2
      U ( (1+(a/c)^2)/6, (8-(a/c)^2)/12, (8-(a/c)^2)/12 ) = (-1 1/2 1/2)*(1-(a/c)^2/2)/3 + (1 1 1)/2
      K (1 0 -1)*(1+(a/c)^2/4)/3

So, if you want to specify a typical circuit, the following might do the work
(the representative points on lines of symmetry are indicated - there are
sometimes more than one way to go from one point to another): X-V-K-Sigma-
Gamma-Lambda-T-Q-W-Y-L-sigma-Gamma-sigma-X. The suggestion is to sample this
path with the following coordinates for the special points X, Gamma, T, L, Gamma, X:

      kptbounds  1/2 0 -1/2   0 0 0    1/2 1/2 1/2  1 1/2 0   1 0 0  1 1/2 1/2

In order to find the lengths of segments (this information is useful to draw
the band structure, with the correct relative scale between special points)
one needs to know the a and c lattice parameters. Also, in what follows, we
omit the $2\pi$ factor sometimes present in the definition of the reciprocal
space vectors. The reciprocal vectors are $( 2/(3a)\: 0\: 1/(3c) )$, $( -1/(3a)\:
1/(\sqrt{3}a)\: 1/(3c) )$, $( -1/(3a)\: -1/(\sqrt{3}a)\: 1/(3c) )$. The lengths of the
above-mentioned segments can be computed as:

$l$(X-Gamma)=$2/(\sqrt{3}a)$=1.155.../a \n
$l$(K-Gamma)=$4(1+(a/c)^2/4)/(3\sqrt{3}a)$ \n
$l$(Gamma-T)=$1/(2c)$ \n
$l$(T-L)=$2/(\sqrt{3}a)$=1.155.../a \n
$l$(T-W)=$4(1-(a/c)^2/2)/(3\sqrt{3}a)$ \n
$l$(L-Gamma)=$\sqrt{4/(a^2)+1/(c^2)}/3$ \n
$l$(Gamma-X)=$2\sqrt{1/(a^2)+1/(c^2)}/3$ \n
""",
),

Variable(
    abivarname="kptgw",
    varset="gw",
    vartype="real",
    topics=['SelfEnergy_basic'],
    dimensions=[3, '[[nkptgw]]'],
    defaultval=MultipleValue(number=None, value=0.0),
    mnemonics="K-PoinTs for GW calculations",
    requires="[[optdriver]] in [4, 7]",
    added_in_version="before_v9",
    text=r"""
For each k-point with number igwpt in the range (1:[[nkptgw]]),
[[kptgw]](1,igwpt) is the reduced coordinate of the k-point where the self-energy
corrections are required while [[bdgw]] (1:2,igwpt) specifies the range of
bands to be considered.

At present, not all k-points are possible. Only those corresponding to the
k-point grid defined with the same repetition parameters ( [[kptrlatt]], or
[[ngkpt]] ) than the GS one, but **without** any shift, are allowed.
""",
),

Variable(
    abivarname="kptnrm",
    varset="basic",
    vartype="real",
    topics=['k-points_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="K - PoinTs NoRMalization",
    added_in_version="before_v9",
    text=r"""
Establishes a normalizing denominator for each k point. Needed only if
[[kptopt]]<=0, otherwise deduced from other input variables.
The k point coordinates as fractions of reciprocal lattice translations are
therefore [[kpt]](mu,ikpt)/[[kptnrm]]. [[kptnrm]] defaults to 1 and can be
ignored by the user. It is introduced to avoid the need for many digits in
representing numbers such as 1/3. It cannot be smaller than 1.0
""",
),

Variable(
    abivarname="kptns",
    varset="internal",
    vartype="real",
    topics=['k-points_internal'],
    dimensions=[3, '[[nkpt]]'],
    mnemonics="K-PoinTs re-Normalized and Shifted",
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
If [[nqpt]] = 0, or if one is doing a reponse calculation, this internal
variable is derived from [[kpt]] and [[kptnrm]]: [[kptns]](1:3,:)=
[[kpt]](1:3,:)/ [[kptnrm]], so that it is [[kpt]] renormalized by [[kptnrm]].

If [[nqpt]] = 1 and one is not doing a ground-state calculation, this internal
variable is derived from [[kpt]], [[kptnrm]] and [[qptn]]: [[kptns]](1:3,:)=
[[kpt]](1:3,:)/ [[kptnrm]]+ [[qptn]](1:3), so that it is [[kpt]] renormalized
by [[kptnrm]], then shifted by [[qptn]](1:3).
""",
),

Variable(
    abivarname="kptns_hf",
    varset="internal",
    vartype="real",
    topics=['k-points_internal'],
    dimensions=[3, '[[nkpthf]]'],
    mnemonics="K-PoinTs re-Normalized and Shifted, for the Hartree-Fock operator",
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
[[kptns_hf]] is the subset of the full Brillouin Zone k point grid for
wavefunctions, used to build the Fock operator, see [[fockdownsampling]].
""",
),

Variable(
    abivarname="kptopt",
    varset="basic",
    vartype="integer",
    topics=['k-points_basic', 'ElecBandStructure_basic'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[nspden]] == 4': 4, 'defaultval': 1}),
    mnemonics="KPoinTs OPTion",
    added_in_version="before_v9",
    text=r"""
Controls the set up of the k-points list. The aim will be to initialize, by
straight reading or by a preprocessing approach based on other input
variables, the following input variables, giving the k points, their number,
and their weight: [[kpt]], [[kptnrm]], [[nkpt]], and, for [[iscf]]/=-2, [[wtk]].

Often, the k points will form a lattice in reciprocal space. In this case, one
will also aim at initializing input variables that give the reciprocal of this
k-point lattice, as well as its shift with respect to the origin: [[ngkpt]] or
[[kptrlatt]], as well as on [[nshiftk]] and [[shiftk]].

A global additional shift can be provided by [[qptn]]

The use of symmetries (spatial and/or time-reversal) is crucial to determine the action of [[kptopt]].

  * 0 --> read directly [[nkpt]], [[kpt]], [[kptnrm]] and [[wtk]].
  * 1 --> rely on [[ngkpt]] or [[kptrlatt]], as well as on [[nshiftk]] and [[shiftk]] to set up the k points.
    Take fully into account the symmetry to generate the k points in the Irreducible Brillouin Zone only,
    with the appropriate weights. (This is the usual mode for GS calculations)

  * 2 --> rely on [[ngkpt]] or [[kptrlatt]], as well as on [[nshiftk]] and [[shiftk]] to set up the k points.
    Take into account only the time-reversal symmetry: k points will be generated in half the Brillouin zone,
    with the appropriate weights.
    (This is the usual mode when preparing or executing a RF calculation at q=(0 0 0) without non-collinear magnetism)

  * 3 --> rely on [[ngkpt]] or [[kptrlatt]], as well as on [[nshiftk]] and [[shiftk]] to set up the k points.
    Do not take into account any symmetry: k points will be generated in the full Brillouin zone, with the appropriate weights.
    (This is the usual mode when preparing or executing a RF calculation at non-zero q, or with non-collinear magnetism)

  * 4 --> rely on [[ngkpt]] or [[kptrlatt]], as well as on [[nshiftk]] and [[shiftk]] to set up the k points.
   Take into account all the symmetries EXCEPT the time-reversal symmetry to generate the k points
   in the Irreducible Brillouin Zone, with the appropriate weights.
   This has to be used when performing calculations with non-collinear magnetism allowed ([[nspden]] = 4)

  * A negative value  --> rely on [[kptbounds]], and [[ndivk]] ([[ndivsm]]) to set up a band structure calculation
    along different lines (allowed only for [[iscf]] == -2). The absolute value of [[kptopt]] gives
    the number of segments of the band structure. Weights are usually irrelevant with this option,
    and will be left to their default value.

In the case of a grid of k points, the auxiliary variables [[kptrlen]],
[[ngkpt]] and [[prtkpt]] might help you to select the optimal grid.
""",
),

Variable(
    abivarname="kptrlatt",
    varset="gstate",
    vartype="integer",
    topics=['k-points_useful'],
    dimensions=[3, 3],
    defaultval=MultipleValue(number=None, value=0),
    mnemonics="K - PoinTs grid: Real space LATTice",
    excludes="[[ngkpt]]",
    added_in_version="before_v9",
    text=r"""
This input variable is used only when [[kptopt]] is positive. It partially
defines the k point grid. The other piece of information is contained in
[[shiftk]]. [[kptrlatt]] cannot be used together with [[ngkpt]].

The values [[kptrlatt]](1:3,1), [[kptrlatt]](1:3,2), [[kptrlatt]](1:3,3) are the
coordinates of three vectors in real space, expressed in the [[rprimd]]
coordinate system (reduced coordinates). They defines a super-lattice in real
space. The k point lattice is the reciprocal of this super-lattice, possibly shifted (see [[shiftk]]).

If neither [[ngkpt]] nor [[kptrlatt]] are defined, ABINIT will automatically
generate a set of k point grids, and select the best combination of
[[kptrlatt]] and [[shiftk]] that allows one to reach a sufficient value of
[[kptrlen]]. See this latter variable for a complete description of this procedure.
""",
),

Variable(
    abivarname="kptrlen",
    varset="gstate",
    vartype="real",
    topics=['k-points_useful'],
    dimensions="scalar",
    defaultval=30.0,
    mnemonics="K - PoinTs grid: Real space LENgth",
    added_in_version="before_v9",
    text=r"""
This input variable is used only when [[kptopt]] is positive and non-zero.

Preliminary explanation:
The k point lattice defined by [[ngkpt]] or [[kptrlatt]] is used to perform
integrations of periodic quantities in the Brillouin Zone, like the density or
the kinetic energy. One can relate the error made by replacing the continuous
integral by a sum over k point lattice to the Fourier transform of the
periodic quantity. Erroneous contributions will appear only for the vectors in
real space that belong to the reciprocal of the k point lattice, except the
origin. Moreover, the expected size of these contributions usually decreases
exponentially with the distance. So, the length of the smallest of these real
space vectors is a measure of the accuracy of the k point grid.

When either [[ngkpt]] or [[kptrlatt]] is defined, [[kptrlen]] is not used as
an input variable, but the length of the smallest vector will be placed in
this variable, and echoed in the output file.

On the other hand, when neither [[ngkpt]] nor [[kptrlatt]] are defined, ABINIT
will automatically generate a large set of possible k point grids, and select
among this set, the grids that give a length of smallest vector LARGER than
[[kptrlen]], and among these grids, the one that, when used with [[kptopt]] = 1,
reduces to the smallest number of k points. Note that this procedure can be
time-consuming. It is worth doing it once for a given unit cell and set of
symmetries, but not use this procedure by default. The best is then to set
[[prtkpt]] = 1, in order to get a detailed analysis of the set of grids.

If some layer of vacuum is detected in the unit cell (see the input variable
[[vacuum]]), the computation of [[kptrlen]] will ignore the dimension related
to the direction perpendicular to the vacuum layer, and generate a bi-
dimensional k point grid. If the system is confined in a tube, a one-
dimensional k point grid will be generated. For a cluster, this procedure will
only generate the Gamma point.
""",
),

Variable(
    abivarname="kssform",
    varset="files",
    vartype="integer",
    topics=['Susceptibility_expert', 'SelfEnergy_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Kohn Sham Structure file FORMat",
    added_in_version="before_v9",
    text=r"""
Governs the choice of the format for the file that contains the Kohn-Sham
electronic structure information, for use in GW calculations, see the input
variables [[optdriver]] and [[nbandkss]].

  * [[kssform]] = 1, a single file.kss (double precision) containing complete information on the Kohn Sham Structure (eigenstates and the pseudopotentials used) will be generated through full diagonalization of the complete Hamiltonian matrix. The file has at the beginning the standard abinit header.
  * [[kssform]] = 3, a single file.kss (double precision) containing complete information on the Kohn Sham Structure (eigenstates and the pseudopotentials used) will be generated through the usual conjugate gradient algorithm (so, a restricted number of states). The file has at the beginning the standard abinit header.

!!! warning

    For the time being, [[istwfk]] must be 1 for all the k-points.
""",
),

Variable(
    abivarname="ldaminushalf",
    varset="paw",
    vartype="integer",
    topics=['LDAminushalf_compulsory'],
    dimensions=['[[ntypat]]'],
    defaultval=MultipleValue(number=None, value=0),
    mnemonics="LDA minus half",
    added_in_version="before_v9",
    text=r"""
For each type of atom, gives whether a LDA-${\frac{1}{2}}$ calculation is to be performed.
**ldaminushalf** =0: the LDA-$\frac{1}{2}$ approach is not used.
**ldaminushalf** =1: the LDA-$\frac{1}{2}$ approach is used.
""",
),

Variable(
    abivarname="lexexch",
    varset="paw",
    vartype="integer",
    topics=['xc_useful'],
    dimensions=['[[ntypat]]'],
    defaultval=-1,
    mnemonics="value of angular momentum L for EXact EXCHange",
    requires="[[useexexch]] == 1",
    added_in_version="before_v9",
    text=r"""
Give for each species the value of the angular momentum (only values 2 or 3
are allowed) on which to apply the exact exchange correction.
""",
),

Variable(
    abivarname="localrdwf",
    varset="paral",
    vartype="integer",
    topics=['parallelism_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="LOCAL ReaD WaveFunctions",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
This input variable is used only when running abinit in parallel. If
[[localrdwf]] = 1, the input wavefunction disk file or the KSS/SCR file in case
of GW calculations, is read locally by each processor, while if
[[localrdwf]] = 0, only one processor reads it, and broadcast the data to the other processors.

The option [[localrdwf]] = 0 is NOT allowed when parallel I/O are activated
(MPI-IO access), i.e. when [[iomode]] == 1.

In the case of a parallel computer with a unique file system, both options are
as convenient for the user. However, if the I/O are slow compared to
communications between processors,, [[localrdwf]] = 0 should be much more
efficient; if you really need temporary disk storage, switch to localrdwf=1).

In the case of a cluster of nodes, with a different file system for each
machine, the input wavefunction file must be available on all nodes if
[[localrdwf]] = 1, while it is needed only for the master node if [[localrdwf]] = 0.
""",
),

Variable(
    abivarname="lotf_classic",
    varset="dev",
    vartype="integer",
    topics=['LOTF_expert'],
    dimensions="scalar",
    defaultval=5,
    mnemonics="LOTF CLASSIC model for glue model",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Glue model used in LOTF. For the moment it is imposed to be 5.
""",
),

Variable(
    abivarname="lotf_nitex",
    varset="dev",
    vartype="integer",
    topics=['LOTF_expert'],
    dimensions="scalar",
    defaultval=10,
    mnemonics="LOTF Number of ITerations",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Set the number of Molecular Dynamics iterations which are computed by LOTF.
""",
),

Variable(
    abivarname="lotf_nneigx",
    varset="dev",
    vartype="integer",
    topics=['LOTF_expert'],
    dimensions="scalar",
    defaultval=5,
    mnemonics="LOTF max Number of NEIGhbours",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Set the max number of Neighbours used in the LOTF method.
For the moment it is imposed to be 40.
""",
),

Variable(
    abivarname="lotf_version",
    varset="dev",
    vartype="integer",
    topics=['LOTF_expert'],
    dimensions="scalar",
    defaultval=2,
    mnemonics="LOTF VERSION of MD algorithm",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Set the MD algorithm in the LOTF method. For the moment it is imposed to be 2.
""",
),

Variable(
    abivarname="lpawu",
    varset="paw",
    vartype="integer",
    topics=['DFT+U_compulsory'],
    dimensions=['[[ntypat]]'],
    defaultval=MultipleValue(number=None, value=-1),
    mnemonics="value of angular momentum L for PAW+U",
    requires="[[usepawu]] == 1 or 2",
    added_in_version="before_v9",
    text=r"""
Give for each species the value of the angular momentum (only values 2 or 3
are allowed)  on which to apply the LDA+U correction.

  * If equal to 2 (d-orbitals)  or 3 (f-orbitals), values of [[upawu]] and  [[jpawu]] are used in the calculation.
  * If equal to -1: do not apply LDA+U correction on the species.
""",
),

Variable(
    abivarname="macro_uj",
    varset="dev",
    vartype="integer",
    topics=['DFT+U_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="MACRO variable that activates the determination of the U and J parameter (for the PAW+U calculations)",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Sets proper input values for the determination of U and J i.e. for [[pawujat]]
(first atom treated with PAW+U), [[irdwfk]] (=1), [[tolvrs]] (=10^(-8)),
[[nstep]] (=255), [[diemix]] (=0.45), [[atvshift]] ([[pawujat]]) [[pawujv]]).
Do not overwrite these variables manually unless you know what you are doing.

  * [[macro_uj]] = 1 (and [[nsppol]] = 2) Standard procedure to determine U on atom pawujat through a shift of the potential on both spin channels.
  * [[macro_uj]] = 1 (and [[nsppol]] = 1) Non standard procedure to determine U from potential shift on atom pawujat (experimental).
  * [[macro_uj]] = 2 (and [[nsppol]] = 2) Non standard procedure to determine U from potential shift on atom pawujat through a shift on spin channel 1 on this atom and the response on this channel (experimental).
  * [[macro_uj]] = 3 (and [[nsppol]] = 2) Standard procedure to determine J from potential shift on spin channel 1 on atom pawujat and response on spin channel 2 (experimental).

Determination of U and J can be done only if the symmetry of the atomic
arrangement is reduced and the atom pawujat is not connected to any other atom
by symmetry relations (either input reduced symmetries manually, define
concerned atom as a separate atomic species or shift concerned atom from ideal
position).
""",
),

Variable(
    abivarname="magcon_lambda",
    varset="gstate",
    vartype="real",
    topics=['MagMom_useful'],
    dimensions="scalar",
    defaultval=0.01,
    mnemonics="MAGnetization CONstraint LAMBDA parameter",
    added_in_version="before_v9",
    text=r"""
This variable gives the amplitude of the penalty function imposed on the
magnetization vectors on each atom (turned on with flag variable
[[magconon]]=1 to 3). Typical values for [[magcon_lambda]] are 0.001 to 0.1. The SCF convergence
will be difficult if [[magcon_lambda]] is too large. If [[magcon_lambda]] is too small, the
penalty will not be very effective and it will give magnetization not close
to the desired [[spinat]] target. In case of convergence problem, it can help
to start with a small value of [[magcon_lambda]] and to increase it by reading the
wavefunction obtained with a lower [[magcon_lambda]] value. See variable [[magconon]] for more details.
""",
),

Variable(
    abivarname="magconon",
    varset="gstate",
    vartype="integer",
    topics=['MagMom_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="turn MAGnetization CONstraint ON",
    added_in_version="before_v9",
    text=r"""
Turns on the imposition of a constraint on the magnetization, using a penalty function. For
each atom, the magnetization is calculated in a sphere (radius [[ratsph]]) and
a penalty function is applied to bring it to the input values of [[spinat]].
The constraint can be either on the direction only ([[magconon]] = 1) or on the full
vector ([[magconon]] = 2). The penalty function has an amplitude
[[magcon_lambda]] that should be neither too big (bad or impossible convergence) nor too small (no effect).
The penalty function is documented in [[cite:Ma2015]] as being a Lagrange
approach, which is a misnomer for the algorithm that they describe. It has the drawback of being unable to deliver
the exact sought value for the magnetization. So, the true Lagrange approach has to be preferred, except for testing purposes.
This is provided by the algorithm governed by the input variable [[constraint_kind]], which is actually also much more flexible
than the implementation corresponding to [[magconon]].
""",
),

Variable(
    abivarname="max_ncpus",
    varset="paral",
    vartype="integer",
    topics=['parallelism_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="MAXimum Number of CPUS",
    added_in_version="before_v9",
    text=r"""
If [[autoparal]] > 1 and [[max_ncpus]] is greater than 0, ABINIT analyzes the
efficiency of the process distribution for each possible number of processors
from 2 to [[max_ncpus]]. After having printed out the efficiency, the code stops.
""",
),

Variable(
    abivarname="maxestep",
    varset="ffield",
    vartype="real",
    topics=['Berry_useful'],
    dimensions="scalar",
    defaultval=0.005,
    mnemonics="MAXimum Electric field STEP",
    requires="[[berryopt]] = 6, 16, or 17",
    added_in_version="before_v9",
    text=r"""
This variable controls the maximum change of electric field when updating the
electric field after each SCF iteration. When the calculation is difficult to
converge, try reducing this value or reducing [[ddamp]]. This variable is used
in finite electric displacement field calculations ([[berryopt]] = 6,16,17).
""",
),

Variable(
    abivarname="maxnsym",
    varset="dev",
    vartype="integer",
    topics=['crystal_expert'],
    dimensions="scalar",
    defaultval=384,
    mnemonics="MAXimum Number of SYMetries",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Gives the maximum number of spatial symmetries allowed in the memory.
The default value is sufficient for most applications. It might have to be
increased in the case of the use of a supercell (unit cell identically
repeated).
""",
),

Variable(
    abivarname="mband",
    varset="internal",
    vartype="integer",
    topics=['BandOcc_internal'],
    dimensions="scalar",
    mnemonics="Maximum number of BANDs",
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
This internal variable derives the maximum number of bands over all k-points
and spin-polarisation from [[nband]](1:nkpt*nsppol).
""",
),

Variable(
    abivarname="mbpt_sciss",
    varset="gw",
    vartype="real",
    topics=['GW_useful', 'Susceptibility_useful', 'SelfEnergy_useful'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="Many Body Perturbation Theory SCISSor operator",
    characteristics=['[[ENERGY]]'],
    requires="[[optdriver]] in [3,4,99]",
    added_in_version="before_v9",
    text=r"""
The scissor operator energy added to the conductions states. In some cases,
it mimics a second iteration self-consistent GW calculation.
""",
),

Variable(
    abivarname="mdf_epsinf",
    varset="gw",
    vartype="real",
    topics=['BSE_useful'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="Model Dielectric Function, EPSilon INFinity",
    requires="[[optdriver]] == 99 and [[bs_coulomb_term]] in [20,21] (Bethe-Salpeter calculations with a model dielectric function",
    added_in_version="before_v9",
    text=r"""
[[mdf_epsinf]] specifies the value of the macroscopic dielectric function used
to model the screening function (see [[cite:Bechstedt1992]]).
The proper spatial symmetry of the screening $W(\mathbf{r},\mathbf{r}^\prime)$ is enforced using
Eq. (7) of [[cite:vonderLinden1988]].
""",
),

Variable(
    abivarname="mdtemp",
    varset="rlx",
    vartype="real",
    topics=['PIMD_compulsory', 'MolecularDynamics_compulsory'],
    dimensions=[2],
    defaultval=[300, 300],
    mnemonics="Molecular Dynamics TEMPeratures",
    added_in_version="before_v9",
    text=r"""
Give the initial and final temperature of the Nose-Hoover thermostat
([[ionmov]] = 8) and Langevin dynamics ([[ionmov]] = 9), in Kelvin. This
temperature will change linearly from the initial temperature **mdtemp(1)** at
itime=1 to the final temperature **mdtemp(2)** at the end of the [[ntime]] timesteps.

In the case of the isokinetic molecular dynamics ([[ionmov]] = 12), **mdtemp(1)** allows ABINIT
to generate velocities ([[vel]]) to start the run if they are not provided by the user or if they all vanish. However **mdtemp(2)** is not used (even if it must be defined to please the parser). If some velocities are non-zero, **mdtemp** is not used, the kinetic energy computed from the velocities is kept constant during the run.
""",
),

Variable(
    abivarname="mdwall",
    varset="rlx",
    vartype="real",
    topics=['MolecularDynamics_expert'],
    dimensions="scalar",
    defaultval=10000.0,
    mnemonics="Molecular Dynamics WALL location",
    commentdefault="the walls are extremely far away",
    added_in_version="before_v9",
    text=r"""
Gives the location (atomic units) of walls on which the atoms will bounce
back when [[ionmov]] = 6, 7, 8 or 9. For each cartesian direction idir=1, 2 or
3, there is a pair of walls with coordinates xcart(idir)=-wall and
xcart(idir)=rprimd(idir,idir)+wall. Supposing the particle will cross the
wall, its velocity normal to the wall is reversed, so that it bounces back.
By default, given in Bohr atomic units (1 Bohr=0.5291772108 Angstroms),
although Angstrom can be specified, if preferred, since [[mdwall]] has the [[LENGTH]] characteristics.
""",
),

Variable(
    abivarname="mem_test",
    varset="dev",
    vartype="integer",
    topics=['Control_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="MEMory TEST",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
This variable controls the memory test done in the `memana` routine. Possible values:

  * 0 no test on the available memory is performed
  * 1 the routine tries to allocate the estimated memory, for testing purposes, and if a failure occurs, the routine stops.
  * 2 like 1, but before stopping, the routine will provide an estimation of the available memory.
""",
),

Variable(
    abivarname="mep_mxstep",
    varset="rlx",
    vartype="real",
    topics=['TransPath_useful'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[imgmov]] == 5': 0.4, 'defaultval': 100.0}),
    mnemonics="Minimal Energy Path search: MaXimum allowed STEP size",
    characteristics=['[[LENGTH]]'],
    added_in_version="before_v9",
    text=r"""
Relevant only when [[imgmov]] = 1 (Steepest-Descent), 2 (String Method) or 5
(Nudged Elastic Band).
The optimizer used to solve the Ordinary Differential Equation (ODE) can be
constrained with a maximum allowed step size for each image. By default this
feature is only activated for Nudged Elastic Band (NEB) and the value is
inspired by [[cite:Sheppard2008]].
Note that the step size is defined for each image as _step = SQRT[SUM(R_i dot
R_i)]_ where the _R_i_ are the positions of the atoms in the cell.
""",
),

Variable(
    abivarname="mep_solver",
    varset="rlx",
    vartype="integer",
    topics=['TransPath_basic'],
    dimensions="scalar",
    mnemonics="Minimal Energy Path ordinary differential equation SOLVER",
    added_in_version="before_v9",
    text=r"""
Relevant only when [[imgmov]] = 2 (String Method) or 5 (Nudged Elastic Band).
Gives the algorithm used to solve the Ordinary Differential Equation (ODE)
when searching for a Minimal Energy Path (MEP).
Possible values can be:

  * 0 --> **Steepest-Descent algorithm** following the (scaled) forces, the scaling factor being [[fxcartfactor]] (forward Euler method).
Compatible with all MEP search methods.

  * 1 --> **Quick-min optimizer** following the (scaled) forces, the scaling factor being [[fxcartfactor]]. The "quick minimizer" improves upon the steepest-descent method by accelerating the system in the direction of the forces. The velocity (of the image) is projected long the force and cancelled if antiparallel to it.
Compatible only with Nudged Elastic Band ([[imgmov]] = 5).
See [[cite:Sheppard2008]].

  * 2 --> **Local Broyden-Fletcher-Goldfarb-Shanno (L-BFGS) algorithm**; each image along the band is minimized with a different instance of the BFGS optimizer.
Compatible only with Nudged Elastic Band ([[imgmov]] = 5).
See [[cite:Sheppard2008]].
IN [[DEVELOP]]PMENT - NOT RELIABLE

  * 3 --> **Global Broyden-Fletcher-Goldfarb-Shanno (GL-BFGS) algorithm**; all images along the band are minimized with a single instance of the BFGS optimizer.
Compatible only with Nudged Elastic Band ([[imgmov]] = 5).
See [[cite:Sheppard2008]].
IN [[DEVELOP]]PMENT - NOT RELIABLE

  * 4 --> **Fourth-order Runge-Kutta method**; the images along the band are moved every four steps (1 <=istep<=[[ntimimage]]) following the Runge-Kutta algorithm, the time step being [[fxcartfactor]].
Compatible only with Simplified String Method ([[imgmov]] = 2 and
[[string_algo]] = 1 or 2).
See [[cite:Weinan2007]].

All of the optimizers can be constrained with a maximum allowed step size for
each image; see [[mep_mxstep]]. This is by default the case of the Nudged
Elastic Band ([[imgmov]] = 5).
""",
),

Variable(
    abivarname="mixprec",
    varset="dev",
    vartype="integer",
    topics=['TuningSpeed_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="MIXed PRECision",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
This variable activates FFT transforms in single precision.
The code thus works in mixed-precision mode in the sense that all the high-level
operations are done in double-precision while the FFT of wavefunctions, densities and potentials
are done with workspace arrays in single precision.

This option **requires** the linkage with external FFT libraries (FFTW3 or MKL-DFTI, see also [[fftalg]])
Tests showed a speedup of ~25% in calculations in which FFTs (in particular fourwf%pot) represent the dominant part.
Typical examples are EPH calculation with [[optdriver]] = 7.

At present ( |today| ), only selected kernels support mixed-precision, in particular MPI-FFTs
in mixed precision **are not yet supported**.
""",
),

Variable(
    abivarname="mgfft",
    varset="internal",
    vartype="integer",
    topics=['Planewaves_internal'],
    dimensions="scalar",
    mnemonics="Maximum of nGFFT",
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
This internal variable contains the maximum of [[ngfft]](1:3).
""",
),

Variable(
    abivarname="mgfftdg",
    varset="internal",
    vartype="integer",
    topics=['Planewaves_internal'],
    dimensions="scalar",
    mnemonics="Maximum of nGFFT for the Double Grid",
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
This internal variable contains the maximum of [[ngfftdg]](1:3).
""",
),

Variable(
    abivarname="mixalch",
    varset="gstate",
    vartype="real",
    topics=['AtomTypes_useful'],
    dimensions=['[[npspalch]]', '[[ntypalch]]'],
    mnemonics="MIXing coefficients for ALCHemical potentials",
    characteristics=['[[EVOLVING]]'],
    added_in_version="before_v9",
    text=r"""
Used for the generation of alchemical pseudoatoms, that is, when [[ntypalch]] is non-zero.

This array gives, for each type of alchemical pseudo-atom (there are
[[ntypalch]] such pseudoatoms), the mixing coefficients of the basic
[[npspalch]] pseudopotentials for alchemical use. For each type of alchemical
pseudoatom, the sum of the mixing coefficients must equal 1.

The actual use of the mixing coefficients is defined by the input variable
[[algalch]]. Note that the masses of the atoms, [[amu]] are also mixed
according to the value of [[mixalch]], by default.

Example 1. Suppose that we want to describe Ba(0.25) Sr(0.75) Ti O$_3$.
The input variables related to the construction of the alchemical Ba(0.25)
Sr(0.75) potential will be:

      npsp   4                 ! 4 pseudopotentials should be read.
      znucl  8 40 56 38        ! The nuclear charges. Note that the two
                               ! atoms whose pseudopotentials are to be mixed
                               ! are mentioned at the end of the series.
      ntypat  3                ! There will be three types of atoms.
      ntypalch   1             ! One pseudoatom will be alchemical.
                               ! Hence, there will be ntyppure=2 pure pseudo-atoms,
                               ! with znucl 8 (O) and 40 (Ti), corresponding to
                               ! the two first pseudopotentials. Out of the
                               ! four pseudopotentials, npspalch=2 are left
                               ! for alchemical purposes, with znucl 56 (Ba)
                               ! and 38 (Sr).
      mixalch    0.25  0.75    ! For that unique pseudo-atom to be
                               ! generated, here are the mixing coefficients,
                               ! to be used to combine the Ba and Sr pseudopotentials.


Example 2. More complicated, and illustrate some minor drawback of the design
of input variables. Suppose that one wants to generate Al(0.25)Ga(0.75)
As(0.10)Sb(0.90).
The input variables will be:

      npsp  4                  ! 4 pseudopotentials should be read
      znucl  13 31 33 51       ! The atomic numbers. All pseudopotentials
                               ! will be used for some alchemical purpose
      ntypat  2                ! There will be two types of atoms.
      ntypalch   2             ! None of the atoms will be "pure".
                               ! Hence, there will be npspalch=4 pseudopotentials
                               !  to be used for alchemical purposes.
      mixalch    0.25  0.75 0.0  0.0   ! This array is a (4,2) array, arranged in the
                 0.0   0.0  0.1  0.9   ! usual Fortran order.

Minor drawback: one should not forget to fill [[mixalch]] with the needed
zero's, in this later case.

In most cases, the use of [[mixalch]] will be as a static (non-evolving)
variable. However, the possibility to have different values of [[mixalch]] for
different images has been coded. A population of cells with different atomic
characteristics can thus be considered, and can be made to evolve, e.g. with a
genetic algorithm (not coded in v7.0.0 though). There is one restriction to
this possibility: the value of [[ziontypat]] for the atoms that are mixed
should be identical.
""",
),

Variable(
    abivarname="mixesimgf",
    varset="rlx",
    vartype="real",
    topics=['CrossingBarriers_useful'],
    dimensions=['[[nimage]]'],
    mnemonics="MIXing Electronic Structure IMAGE Factors",
    added_in_version="before_v9",
    text=r"""
Used in the algorithm Linear Combination of Constrained DFT Energies, that is, when [[imgmov]]==6.

This array gives, for each one of the [[nimage]] images, the factor
by which the total energies for systems with same geometry but different electronic structures (occupation numbers) are linearly combined.
The sum of these factors must equal 1.
""",
),

Variable(
    abivarname="mpw",
    varset="internal",
    vartype="integer",
    topics=['Planewaves_internal'],
    dimensions="scalar",
    mnemonics="Maximum number of Plane Waves",
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
This internal variable gives the maximum of the number of plane waves over all
k-points. It is computed from [[ecut]] and the description of the cell,
provided by [[acell]], [[rprim]], and/or [[angdeg]].
""",
),

Variable(
    abivarname="mqgrid",
    varset="dev",
    vartype="integer",
    topics=['Planewaves_expert'],
    dimensions="scalar",
    defaultval=3001,
    mnemonics="Maximum number of Q-space GRID points for pseudopotentials",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Govern the size of the one-dimensional information related to
pseudopotentials, in reciprocal space: potentials, or projector functions.
""",
),

Variable(
    abivarname="mqgriddg",
    varset="paw",
    vartype="integer",
    topics=['PAW_expert'],
    dimensions="scalar",
    defaultval=3001,
    mnemonics="Maximum number of Q-wavevectors for the 1-dimensional GRID  for the Double Grid in PAW",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Maximum number of wavevectors used to sample the local part of the potential,
in PAW. Actually referred to as mqgrid_vl internally. Should change name to
the latter. See also [[mqgrid]].
""",
),

Variable(
    abivarname="natcon",
    varset="rlx",
    vartype="integer",
    topics=['GeoConstraints_useful'],
    dimensions=['[[nconeq]]'],
    defaultval=0,
    mnemonics="Number of AToms in CONstraint equations",
    characteristics=['[[NO_MULTI]]', '[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Gives the number of atoms appearing in each of the [[nconeq]] independent
equations constraining the motion of atoms during structural optimization or
molecular dynamics (see [[nconeq]], [[iatcon]], and [[wtatcon]]).
""",
),

Variable(
    abivarname="natfix",
    varset="rlx",
    vartype="integer",
    topics=['GeoConstraints_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of Atoms that are FIXed",
    characteristics=['[[INPUT_ONLY]]'],
    commentdefault="(no atoms held fixed)",
    added_in_version="before_v9",
    text=r"""
Gives the number of atoms (not to exceed [[natom]]) which are to be held fixed
during a structural optimization or molecular dynamics.
When [[natfix]] > 0, [[natfix]] entries should be provided in array [[iatfix]]
.
""",
),

Variable(
    abivarname="natfixx",
    varset="rlx",
    vartype="integer",
    topics=['GeoConstraints_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of Atoms that are FIXed along the X direction",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Gives the number of atoms (not to exceed [[natom]]) which are to be held fixed
along the X direction during a structural optimization or molecular dynamics.
When [[natfixx]] > 0, [[natfixx]] entries should be provided in array [[iatfixx]].
""",
),

Variable(
    abivarname="natfixy",
    varset="rlx",
    vartype="integer",
    topics=['GeoConstraints_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of Atoms that are FIXed along the Y direction",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Gives the number of atoms (not to exceed [[natom]]) which are to be held fixed
along the Y direction during a structural optimization or molecular dynamics.
When [[natfixy]] > 0, [[natfixy]] entries should be provided in array [[iatfixy]]
""",
),

Variable(
    abivarname="natfixz",
    varset="rlx",
    vartype="integer",
    topics=['GeoConstraints_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of Atoms that are FIXed along the Z direction",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Gives the number of atoms (not to exceed [[natom]]) which are to be held fixed
along the Z direction during a structural optimization or molecular dynamics.
When [[natfixz]] > 0, [[natfixz]] entries should be provided in array [[iatfixz]].
""",
),

Variable(
    abivarname="natom",
    varset="basic",
    vartype="integer",
    topics=['crystal_basic', 'SmartSymm_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Number of ATOMs",
    added_in_version="before_v9",
    text=r"""
Gives the total number of atoms in the unit cell. Default is 1 but you will
obviously want to input this value explicitly.
Note that [[natom]] refers to all atoms in the unit cell, not only to the
irreducible set of atoms in the unit cell (using symmetry operations, this set
allows one to recover all atoms). If you want to specify only the irreducible set
of atoms, use the symmetriser, see the input variable [[natrd]].
""",
),

Variable(
    abivarname="natpawu",
    varset="internal",
    vartype="integer",
    topics=['DFT+U_internal'],
    dimensions="scalar",
    mnemonics="Number of AToms on which PAW+U is applied",
    characteristics=['[[INTERNAL_ONLY]]'],
    requires="[[usepawu]] == 1",
    added_in_version="before_v9",
    text=r"""
This internal variable gives the number of atoms on which the LDA/GGA+U method
is applied. This value is determined from [[lpawu]].
""",
),

Variable(
    abivarname="natrd",
    varset="geo",
    vartype="integer",
    topics=['AtomManipulator_basic', 'SmartSymm_basic'],
    dimensions="scalar",
    defaultval="[[natom]]",
    mnemonics="Number of AToms ReaD",
    added_in_version="before_v9",
    text=r"""
Gives the number of atoms to be read from the input file, in the case the atom
manipulator or the smart symmetriser is used. In this case, [[natrd]] is also
used to dimension the array [[typat]], and the arrays [[xred]], [[xangst]] and [[xcart]].
Must take into account the vacancies (see [[vacnum]] and [[vaclst]]).
Despite possible vacancies, cannot be bigger than [[natom]].
""",
),

Variable(
    abivarname="natsph",
    varset="gstate",
    vartype="integer",
    topics=['printing_prdos', 'ElecBandStructure_useful', 'ElecDOS_useful'],
    dimensions="scalar",
    defaultval="[[natom]]",
    mnemonics="Number of ATomic SPHeres for the atom-projected density-of-states",
    requires="[[prtdos]] == 3 or [[pawfatbnd]] in [1,2]",
    added_in_version="before_v9",
    text=r"""
[[natsph]] gives the number of atoms around which the sphere for atom-projected
density-of-states will be built, in the [[prtdos]] = 3 case. The
indices of these atoms are given by [[iatsph]]. The radius of these spheres is given by [[ratsph]].
If [[pawfatbnd]] = 1 or 2, it gives the number of atoms around which atom-projected
band structure will be built (the indices of these atoms are given by [[iatsph]]).
""",
),

Variable(
    abivarname="natsph_extra",
    varset="gstate",
    vartype="integer",
    topics=['printing_prdos'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of ATomic SPHeres for the l-projected density-of-states in EXTRA set",
    requires="[[prtdos]] == 3 or [[pawfatbnd]] in [1,2]",
    added_in_version="before_v9",
    text=r"""
[[natsph_extra]] gives the number of extra spheres for which the angular-
momentum-projected density-of-states will be built, in the [[prtdos]] = 3 case.
The radius of these spheres is given by [[ratsph_extra]]. This simulates the
STS signal for an STM tip atom placed at the sphere position, according to the
chemical nature of the tip (s- p- d- wave etc...).
If [[pawfatbnd]] = 1 or 2, it gives the number of spheres in which $l$-projected
band structure will be built.
The position of the spheres is given by the [[xredsph_extra]] variable.
""",
),

Variable(
    abivarname="natvshift",
    varset="ffield",
    vartype="integer",
    topics=['DFT+U_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of ATomic potential (V) energy SHIFTs (per atom)",
    requires="[[usepawu]] /= 0, [[atvshift]]",
    added_in_version="before_v9",
    text=r"""
Number of atomic potential energy shifts (per atom), to be used to define the
array [[atvshift]]. If non-zero, only two possibilities exist: 5 for d states
(with [[lpawu]] = 2), and 7 for f states (with [[lpawu]] = 3). If non-zero, one
should define [[usepawu]], [[lpawu]] and [[atvshift]].
""",
),

Variable(
    abivarname="nband",
    varset="basic",
    vartype="integer",
    topics=['BandOcc_basic', 'GW_basic', 'RPACorrEn_basic', 'Susceptibility_basic'],
    dimensions="scalar",
    mnemonics="Number of BANDs",
    commentdefault=" the estimated number of occupied bands +1 (TODO provide the mathematical formulation)",
    added_in_version="before_v9",
    text=r"""
Gives number of bands, occupied plus possibly unoccupied, for which
wavefunctions are being computed along with eigenvalues.
Note: if the parameter [[occopt]] (see below) is not set to 2, [[nband]] is a
scalar integer, but if the parameter [[occopt]] is set to 2, then [[nband]]
must be an array [[nband]]([[nkpt]]* [[nsppol]]) giving the number of bands
explicitly for each k point. This option is provided in order to allow the
number of bands treated to vary from k point to k point.
For the values of [[occopt]] not equal to 0 or 2, [[nband]] can be omitted.
The number of bands will be set up thanks to the use of the variable
[[fband]]. The present Default will not be used.

If [[nspinor]] is 2, nband must be even for each k point.

In the case of a GW calculation ([[optdriver]] = 3 or 4), [[nband]] gives the
number of bands to be treated to generate the screening (susceptibility and
dielectric matrix), as well as the self-energy. However, to generate the _KSS
file (see [[kssform]]) the relevant number of bands is given by [[nbandkss]].
""",
),

Variable(
    abivarname="nbandhf",
    varset="basic",
    vartype="integer",
    topics=['Hybrids_useful'],
    dimensions="scalar",
    mnemonics="Number of BANDs for (Hartree)-Fock exact exchange",
    commentdefault="the estimated number of occupied bands (TODO: provide the mathematical formulation)",
    added_in_version="before_v9",
    text=r"""
Gives the maximum number of occupied bands with which Fock exact exchange is
being computed for the wavefunctions.
""",
),

Variable(
    abivarname="nbandkss",
    varset="gw",
    vartype="integer",
    topics=['GW_useful', 'Susceptibility_useful', 'SelfEnergy_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of BANDs in the KSS file",
    added_in_version="before_v9",
    text=r"""
This input variable is used for the preparation of a GW calculation: it is
used in a GS run (where [[optdriver]] = 0) to generate a _KSS file. In this run,
[[nbandkss]] should be non-zero. The generated _KSS file can be subsequently
used to calculate the irreducible polarizabilty $\chi^{(0)}_{KS}$ using
[[optdriver]] = 3 or to calculate GW corrections setting [[optdriver]] = 4.

  * If [[nbandkss]] = 0, no _KSS file is created.
  * If [[nbandkss]] = -1, all the available eigenstates (energies and eigenfunctions) are stored in the abo_KSS file at the end of the ground state calculation. The number of states is forced to be the same for all k-points: it will be the minimum of the number of plane waves over all k-points.
  * If [[nbandkss]] is greater than 0, abinit stores (about) [[nbandkss]] eigenstates in the abo_KSS file. This number of states is forced to be the same for all k-points.

See [[npwkss]] for the selection of the number of the planewave components of
the eigenstates to be stored.
The input variable [[iomode]] can be used to read and write KSS files
according to different fileformat (presently only [[iomode]] = 0 and 3 are
available in the GW part).
The precision of the KSS file can be tuned through the input variable [[kssform]].
For more details about the format of the abo_KSS file, see the routine outkss.F90.

!!! warning

    For the time being, [[istwfk]] must be 1 for all the k-points in order to generate a _KSS file.
""",
),

Variable(
    abivarname="nbdblock",
    varset="dev",
    vartype="integer",
    topics=['TuningSpeed_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Number of BanDs in a BLOCK",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
In case of non-standard, blocked algorithms for the optimization of the
wavefunctions (that is, if [[wfoptalg]] = 4):

  * if [[wfoptalg]] = 4, [[nbdblock]] defines the number of blocks (the number of bands in the block is then [[nband]]/[[nbdblock]] ).
""",
),

Variable(
    abivarname="nbdbuf",
    varset="gstate",
    vartype="integer",
    topics=['SCFControl_useful', 'BandOcc_useful'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[optdriver]] == 0 and [[iscf]]<0': '2*[[nspinor]]',
 '[[optdriver]] == 1 and 3<=[[occopt]] and [[occopt]]<= 8': '2*[[nspinor]]',
 'defaultval': 0}),
    mnemonics="Number of BanDs for the BUFfer",
    added_in_version="before_v9",
    text=r"""
[[nbdbuf]] gives the number of bands, the highest in energy, that, among the
[[nband]] bands, are to be considered as part of a buffer. This concept is
useful in three situations: in non-self-consistent calculations, for the
determination of the convergence tolerance; for response functions of metals,
to avoid instabilities, and also when finite electric fields or non-linear
responses (with electric field perturbations) are considered. For the two
first, the need of a buffer is a natural requirement of the problem, so that
the default value is changed to 2 automatically, as explained in the
following. The third case is only for implementation convenience.

In non-self-consistent GS calculations ([[iscf]]<0), the highest levels might
be difficult to converge, if they are degenerate with another level, that does
not belong to the set of bands treated. Then, it might take extremely long to
reach [[tolwfr]], although the other bands are already extremely well-
converged, and the energy of the highest bands (whose residual are not yet
good enough), is also rather well converged.
In response to this problem, for non-zero [[nbdbuf]], the largest residual
(residm), to be later compared with [[tolwfr]], will be computed only in the
set of non-buffer bands (this modification applies for non-self-consistent as
well as self-consistent calculation, for GS as well as RF calculations).
For a GS calculation, with [[iscf]]<0, supposing [[nbdbuf]] is not initialized
in the input file, then ABINIT will overcome the default [[nbdbuf]] value, and
automatically set [[nbdbuf]] to 2.

In metallic RF calculations, in the conjugate gradient optimisation of first-
order wavefunctions, there is an instability situation when the q wavevector
of the perturbation brings the eigenenergy of the highest treated band at some
k point higher than the lowest untreated eigenenergy at some k+q point. If one
accepts a buffer of frozen states, this instability can be made to disappear.
Frozen states receive automatically a residual value of -0.1
For a RF calculation, with 3<=[[occopt]]<=7, supposing [[nbdbuf]] is not
initialized in the input file, then ABINIT will overcome the default
[[nbdbuf]] value, and automatically set [[nbdbuf]] to 2. This value might be
too low in some cases.

Also, the number of active bands, in all cases, is imposed to be at least 1,
irrespective of the value of [[nbdbuf]].
""",
),

Variable(
    abivarname="nberry",
    varset="ffield",
    vartype="integer",
    topics=['Berry_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Number of BERRY phase computations",
    requires="[[berryopt]] = 1, 2, or 3",
    added_in_version="before_v9",
    text=r"""
Gives the number of Berry phase computations of polarisation, or finite-
difference estimations of the derivative of wavefunctions with respect to the
wavevector, each of which might be characterized by a different change of
wavevector [[kberry]].

When equal to 0, no Berry phase calculation of polarisation is performed. The
maximal value of [[nberry]] is 20.

Note that the computation of the polarisation for a set of bands having
different occupation numbers is meaningless (although in the case of spin-
polarized calculations, the spin up bands might have an identical occupation
number, that might differ from the identical occupation number of spin down
bands). Although meaningless, ABINIT will perform such computation, if
required by the user. The input variable [[bdberry]] governs the set of bands
for which a Berry phase is computed.

For the [[berryopt]] = 1, 2, and 3 cases, spinor wavefunctions are not
allowed, nor are parallel computations.
""",
),

Variable(
    abivarname="nc_xccc_gspace",
    varset="dev",
    vartype="integer",
    topics=['Planewaves_expert'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[usepaw]] == 0': 0, '[[usepaw]] == 1': 1, 'defaultval': 0}),
    mnemonics="Norm-Conserving pseudopotentials - use XC Core-Correction in G-SPACE",
    characteristics=['[[DEVELOP]]'],
    commentdefault="0 when [[usepaw]] = 0, 1 when [[usepaw]] = 1",
    added_in_version="before_v9",
    text=r"""
Historically, Abinit treats the model core charge used for the non-linear core
correction in real space. Alternatively, it is possible to instruct the code
to compute the core charge in G-space following the same approach used in the
PAW code. The G-space formalism is more accurate than the interpolation in
real space, especially when derivatives of the model core charge are needed,
e.g. DFPT. Preliminary tests showed that the violation of the acoustic sum
rule is reduced when [[nc_xccc_gspace]] == 1, especially for LDA. It is worth
stressing, however, that [[nc_xccc_gspace]] == 1 should be used only in
conjunction with NC pseudos with a model core charge that decays quickly in
G-space. Several NC pseudos available in the Abinit table are not optimized
for the G-space formalism and users are strongly invited to perform
convergence studies with respect to [[ecut]] before activating this option in production runs.
""",
),

Variable(
    abivarname="nconeq",
    varset="rlx",
    vartype="integer",
    topics=['GeoConstraints_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of CONstraint EQuations",
    characteristics=['[[NO_MULTI]]'],
    added_in_version="before_v9",
    text=r"""
Gives the number of independent equations constraining the motion of atoms
during structural optimization or molecular dynamics (see [[natcon]], [[iatcon]], and [[wtatcon]]).
""",
),

Variable(
    abivarname="nctime",
    varset="dev",
    vartype="integer",
    topics=['MolecularDynamics_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="NetCdf TIME between output of molecular dynamics informations",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
When [[nctime]] is non-zero, the molecular dynamics information is output in
NetCDF format, every [[nctime]] time step. Here is the content of an example file:

```
    netcdf md32.outH_moldyn1 {
    dimensions:
       time = UNLIMITED; // (11 currently)
       DimTensor = 6;
       DimCoord = 3;
       NbAtoms = 32;
       DimVector = 3;
       DimScalar = 1;
    variables:
       double E_pot(time);
          E_pot:units = "hartree";
       double E_kin(time);
          E_kin:units = "hartree";
       double Stress(time, DimTensor);
          Stress:units = "hartree/Bohr^3";
       double Position(time, DimCoord, NbAtoms);
          Position:units = "Bohr";
       double Celerity(time, DimCoord, NbAtoms);
          Celerity:units = "Bohr/(atomic time unit)";
       double PrimitiveVector1(DimVector);
       double PrimitiveVector2(DimVector);
       double PrimitiveVector3(DimVector);
       double Cell_Volume(DimScalar);
          Cell_Volume:units = "Bohr^3";
    }
```
""",
),

Variable(
    abivarname="ndivk",
    varset="gstate",
    vartype="integer",
    topics=['k-points_useful', 'ElecBandStructure_useful'],
    dimensions=['abs([[kptopt]])'],
    mnemonics="Number of DIVisions of K lines",
    characteristics=['[[INPUT_ONLY]]'],
    commentdefault="Will be generated automatically from [[ndivsm]] if the latter is defined.",
    excludes="[[ndivsm]]",
    requires="[[kptopt]] < 0",
    added_in_version="before_v9",
    text=r"""
Gives the number of divisions of each of the segments of the band structure,
whose path is determined by [[kptopt]] and [[kptbounds]]. In this case, the
absolute value of [[kptopt]] is the number of such segments.

For example, suppose that the number of segment is just one ([[kptopt]] = -1), a
value [[ndivk]] = 4 will lead to the computation of points with relative
coordinates 0.0, 0.25, 0.5, 0.75 and 1.0, along the segment in consideration.

Now, suppose that there are two segments ([[kptopt]] = -2), with [[ndivk]](1)=4
and [[ndivk]](2)=2, the computation of the eigenvalues will be done at 7
points, 5 belonging to the first segment, with relative coordinates 0.0, 0.25,
0.5, 0.75 and 1.0, the last one being also the starting point of the next
segment, for which two other points must be computed, with relative
coordinates 0.5 and 1.0.

It is easy to compute disconnected circuits (non-chained segments), by
separating the circuits with the value [[ndivk]] = 1 for the intermediate
segment connecting the end of one circuit with the beginning of the next one
(in which case no intermediate point is computed along this segment).

Alternatively it is possible to generate automatically the array [[ndivk]] by
just specifying the number of divisions for the smallest segment. See the
related input variable [[ndivsm]].
""",
),

Variable(
    abivarname="ndivsm",
    varset="gstate",
    vartype="integer",
    topics=['k-points_useful', 'ElecBandStructure_basic'],
    dimensions="scalar",
    mnemonics="Number of DIVisions for the SMallest segment",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
This variable defines the number of divisions used to sample the smallest
segment of the circuit employed in a band structure calculations (see related
input variables [[kptopt]] and [[kptbounds]]). If [[ndivsm]] is given in the
input file, there is no need to specify the number of divisions to be used for
the other segments. Indeed [[ndivk]] is automatically calculated inside the
code in order to generate a path where the number of divisions in each segment
is proportional to the length of the segment itself. This option is activated
only when [[kptopt]] is negative. In this case, the absolute value of
[[kptopt]] is the number of such segments.
""",
),

Variable(
    abivarname="ndtset",
    varset="basic",
    vartype="integer",
    topics=['multidtset_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of DaTaSETs",
    characteristics=['[[NO_MULTI]]'],
    added_in_version="before_v9",
    text=r"""
Gives the number of data sets to be treated.
If 0, means that the multi-data set treatment is not used, so that the root
filenames will not be appended with _DSx, where 'x' is the dataset index
defined by the input variable [[jdtset]], and also that input names with a
dataset index are not allowed. Otherwise, [[ndtset]] = 0 is equivalent to [[ndtset]] = 1.
""",
),

Variable(
    abivarname="ndynimage",
    varset="internal",
    vartype="integer",
    topics=['PIMD_internal'],
    dimensions="scalar",
    mnemonics="Number of DYNamical IMAGEs",
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
This internal variable gives the number of dynamical images, immediately
deduced from the number of non-zero values present in [[dynimage]]. It is used
to dimension many memory-consuming arrays (one copy for each image), e.g. the
wavefunction array (cg), the density array (rho), etc.
""",
),

Variable(
    abivarname="neb_algo",
    varset="rlx",
    vartype="integer",
    topics=['TransPath_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Nudged Elastic Band ALGOrithm",
    requires="[[imgmov]] == 5",
    added_in_version="before_v9",
    text=r"""
Gives the variant of the NEB method used. Possible values can be:

  * 0 --> **Original NEB method**.
See [[cite:Berne1998]] pp. 385-404

  * 1 --> **NEB + improved tangent**.
The Improved Tangent Method builds on the NEB with an improved estimate of the
tangent direction and a resulting change of the component of the spring force
acting on the images.
See [[cite:Henkelman2000]]

  * 2 --> **Climbing-Image NEB (CI-NEB)**.
The CI-NEB method constitutes a small modification to the NEB method.
Information about the shape of the MEP is retained, but a rigorous convergence
to a saddle point is also obtained. By default the spring constants are
variable (see [[neb_spring]]). As the image with the highest energy has to be
identified, the calculation begins with several iterations of the standard NEB
algorithm. The effective CI-NEB begins at the [[cineb_start]] iteration.
See [[cite:Henkelman2000a]]

Note that, in all cases, it is possible to define the value of the spring
constant connecting images with [[neb_spring]], keeping it constant or
allowing it to vary between 2 values (to have higher resolution close to the saddle point).
""",
),

Variable(
    abivarname="neb_spring",
    varset="rlx",
    vartype="real",
    topics=['TransPath_useful'],
    dimensions=[2],
    defaultval=ValueWithConditions({'[[neb_algo]] == 2': ValueWithUnit(units='Hartree/Bohr^2', value=[0.02, 0.15]), 'defaultval': ValueWithUnit(units='Hartree/Bohr^2', value=[0.05, 0.05])}),
    mnemonics="Nudged Elastic Band: SPRING constant",
    requires="[[imgmov]] == 5",
    added_in_version="before_v9",
    text=r"""
Gives the minimal and maximal values (in Hartree/Bohr^2) of the spring constant connecting images
for the NEB method.
In the standard "Nudged Elastic Band" method, the spring constant is constant
along the path, but, in order to have higher resolution close to the saddle
point, it can be better to have stronger springs close to it.
See [[cite:Henkelman2000a]]
""",
),

Variable(
    abivarname="nelect",
    varset="internal",
    vartype="real",
    topics=['BandOcc_internal'],
    dimensions="scalar",
    defaultval="[[AUTO_FROM_PSP]]",
    mnemonics="Number of ELECTrons",
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
This internal variable gives the number of electrons per unit cell, as
computed from the sum of the valence electrons related to each atom (given in
the pseudopotential, where it is called "zion"), and the input variable [[charge]]:
[[nelect]] = zion-[[charge]].
""",
),

Variable(
    abivarname="nfft",
    varset="internal",
    vartype="integer",
    topics=['Planewaves_internal'],
    dimensions="scalar",
    mnemonics="Number of FFT points",
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
If space parallelization is not used (that is, if [[paral_kgb]] == 0), this
internal variable gives the number of Fast Fourier Transform points in the
grid generated by [[ngfft]](1:3). It is simply the product of the three
components of [[ngfft]].

If space parallelisation is used (that is, if [[paral_kgb]] == 1), then it
becomes the number of Fast Fourier Transform points attributed to the
particular processor. It is no longer the above-mentioned simple product, but
a number usually close to this product divided by the number of processors on
which the space is shared.
""",
),

Variable(
    abivarname="nfftdg",
    varset="internal",
    vartype="integer",
    topics=['Planewaves_internal'],
    dimensions="scalar",
    mnemonics="Number of FFT points for the Double Grid",
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
If space parallelisation is not used (that is, if [[paral_kgb]] == 0), this
internal variable gives the number of Fast Fourier Transform points in the
(double) grid generated by [[ngfftdg]](1:3). It is simply the product of the
three components of [[ngfftdg]].

If space parallelisation is used (that is, if [[paral_kgb]] == 1), then it
becomes the number of Fast Fourier Transform points attributed to the
particular processor. It is no longer the above-mentioned simple product, but
a number usually close to this product divided by the number of processors on
which the space is shared.
""",
),

Variable(
    abivarname="nfreqim",
    varset="gw",
    vartype="integer",
    topics=['FrequencyMeshMBPT_basic', 'RPACorrEn_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of FREQuencies along the IMaginary axis",
    requires="[[optdriver]] == 3 and [[gwcalctyp]] in [2,12,22,9,19,29]",
    added_in_version="before_v9",
    text=r"""
[[nfreqim]] sets the number of pure imaginary frequencies used to calculate
the dielectric matrix in order to perform the numerical integration of the GW self-energy.
""",
),

Variable(
    abivarname="nfreqmidm",
    varset="gw",
    vartype="integer",
    topics=['FrequencyMeshMBPT_basic'],
    dimensions="scalar",
    mnemonics="Nth FREQuency Moment of the Imaginary part of the Dielectric Matrix",
    requires="[[optdriver]] == 4",
    added_in_version="before_v9",
    text=r"""
depending on the value of [[nfreqmidm]] will calculate the frequency moment of
the dielectric matrix or its inverse,

  * if [[nfreqmidm]] is positive: calculate (nth=[[nfreqmidm]]) frequency moment of the dielectric matrix.
  * if [[nfreqmidm]] is negative: calculate (nth=[[nfreqmidm]]) frequency moment of the inverse dielectric matrix.
  * if [[nfreqmidm]] = 0: calculate first frequency moment of the full polarizability.

See [[cite:Taut1985]].
""",
),

Variable(
    abivarname="nfreqre",
    varset="gw",
    vartype="integer",
    topics=['FrequencyMeshMBPT_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of FREQuencies along the REal axis",
    requires="[[optdriver]] == 3 and [[gwcalctyp]] in [2,12,22,9,19,29]",
    added_in_version="before_v9",
    text=r"""
[[nfreqre]] sets the number of real frequencies used to calculate the
dielectric matrix in order to perform the numerical integration of the GW self-energy.

It can be used also in case of GW calculations with plasmon-pole models, _i.e._
[[gwcalctyp]]<10, to reduce the number of frequencies used to evaluate the
dielectric matrix from the (default) two to one frequency (omega=0) by setting
[[nfreqre]] = 1. This might be a good idea in case one is planning to use
[[ppmodel]] > 1. This will force the calculation of the screening on a single
frequency ($\omega=0$) and hence reduce memory and disk space requirement. The
only draw back is that the user will not be able to perform self energy
calculation using [[ppmodel]] = 1, since in the last case the dielectric matrix
calculated on two frequencies is required. If the user is not sure which
ppmodel to use, then s/he is not advised to use this input variable. Using the
default values, one must be able to get a screening file that can be used with any [[ppmodel]].
""",
),

Variable(
    abivarname="nfreqsp",
    varset="gw",
    vartype="integer",
    topics=['SelfEnergy_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of FREQuencies for the SPectral function",
    requires="[[optdriver]] == 4",
    added_in_version="before_v9",
    text=r"""
[[nfreqsp]] defines the number of real frequencies used to calculate the
spectral function of the GW Green's function.
""",
),

Variable(
    abivarname="ngfft",
    varset="gstate",
    vartype="integer",
    topics=['Planewaves_useful'],
    dimensions=[3],
    defaultval=[0, 0, 0],
    mnemonics="Number of Grid points for Fast Fourier Transform",
    commentdefault="(automatic selection of optimal values)",
    added_in_version="before_v9",
    text=r"""
Gives the size of fast Fourier transform (FFT) grid in three dimensions. Each
number must be composed of the factors 2, 3, and 5 to be consistent with the
radices available in our FFT.

If no [[ngfft]] is provided or if [[ngfft]] is set to 0 0 0, the code will automatically
provide an optimal set of [[ngfft]] values, based on [[acell]], [[rprim]] and [[ecut]]
(see also [[boxcutmin]] for speed/accuracy concerns).
This is the recommended procedure, of course.
The total number of FFT points is the product: [[ngfft]](1) x [[ngfft]](2) x [[ngfft]](3)=[[nfft]].

When [[ngfft]] is made smaller than recommended values (e.g. by setting
[[boxcutmin]] to a value smaller than 2.0 or by setting [[ngfft]] manually),
the code runs faster and the equations in effect are approximated by a low
pass Fourier filter. The code reports to standard output (unit 06) a parameter
"boxcut" which is the smallest ratio of the FFT box side to the $\GG$ vector basis
sphere diameter. When boxcut is less than 2 the Fourier filter approximation
is being used. When boxcut gets less than about 1.5 the approximation may be
too severe for realistic results and should be tested against larger values of
[[ngfft]]. When boxcut is larger than 2, [[ngfft]] could be reduced without
loss of accuracy. In this case, the small variations that are observed are
solely due to the xc quadrature, that may be handled with [[intxc]] = 1 to even
reduce this effect.

Internally, [[ngfft]] is an array of size 18. The present components are
stored in [[ngfft]](1:3), while

  * [[ngfft]](4:6) contains slightly different (larger) values, modified for efficiency of the FFT
  * [[ngfft]](7) is [[fftalg]]
  * [[ngfft]](8) is [[fftcache]]
  * [[ngfft]](9) is set to 0 if the parallelization of the FFT is not activated, while it is set to 1 if it is activated.
  * [[ngfft]](10) is the number of processors of the FFT group
  * [[ngfft]](11) is the index of the processor in the group of processors
  * [[ngfft]](12) is n2proc, the number of x-z planes, in reciprocal space, treated by the processor
  * [[ngfft]](13) is n3proc, the number of x-y planes, in real space, treated by the processor
  * [[ngfft]](14) is mpi_comm_fft, the handle on the MPI communicator in charge of the FFT parallelisation
  * [[ngfft]](15:18) are not yet used

The number of points stored by this processor in real space is n1*n2*n3proc,
while in reciprocal space, it is n1*n2proc*n3.
""",
),

Variable(
    abivarname="ngfftdg",
    varset="paw",
    vartype="integer",
    topics=['PAW_useful'],
    dimensions=[3],
    defaultval=[0, 0, 0],
    mnemonics="Number of Grid points for Fast Fourier Transform: Double Grid",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
This variable has the same meaning as ngfft (gives the size of fast Fourier
transform (fft) grid in three dimensions) but concerns the "double grid" only
used for PAW calculations.
""",
),

Variable(
    abivarname="ngkpt",
    varset="basic",
    vartype="integer",
    topics=['k-points_basic'],
    dimensions=[3],
    defaultval=[0, 0, 0],
    mnemonics="Number of Grid points for K PoinTs generation",
    characteristics=['[[INPUT_ONLY]]'],
    excludes="[[kptrlatt]]",
    requires="[[kptopt]] >=0,",
    added_in_version="before_v9",
    text=r"""
Used when [[kptopt]] >= 0, if [[kptrlatt]] has not been defined ([[kptrlatt]]
and [[ngkpt]] are exclusive of each other).
Its three positive components give the number of k points of Monkhorst-Pack
grids (defined with respect to primitive axis in reciprocal space) in each of
the three dimensions. [[ngkpt]] will be used to generate the corresponding
[[kptrlatt]] input variable. The use of [[nshiftk]] and [[shiftk]], allows one to
generate shifted grids, or Monkhorst-Pack grids defined with respect to
conventional unit cells.

When [[nshiftk]] = 1, [[kptrlatt]] is initialized as a diagonal (3x3) matrix,
whose diagonal elements are the three values [[ngkpt]](1:3). When [[nshiftk]]
is greater than 1, ABINIT will try to generate [[kptrlatt]] on the basis of
the primitive vectors of the k-lattice: the number of shifts might be reduced,
in which case [[kptrlatt]] will not be diagonal anymore.

Monkhorst-Pack grids are usually the most efficient when their defining
integer numbers are even. For a measure of the efficiency, see the input
variable [[kptrlen]].
""",
),

Variable(
    abivarname="ngqpt",
    varset="gstate",
    vartype="integer",
    topics=['q-points_basic'],
    dimensions=[3],
    defaultval=[0, 0, 0],
    mnemonics="Number of Grid points for Q PoinTs generation",
    characteristics=['[[INPUT_ONLY]]'],
    excludes="[[qptrlatt]]",
    requires="[[nqpt]] == 1 and [[kptopt]] >= 0",
    added_in_version="before_v9",
    text=r"""
At variance with [[ngkpt]], note that only one q point is selected per dataset
(see [[iqpt]]).
Its three positive components give the number of q points of Monkhorst-Pack
grids (defined with respect to primitive axis in reciprocal space) in each of
the three dimensions. The use of [[nshiftq]] and [[shiftq]], allows one to
generate shifted grids, or Monkhorst-Pack grids defined with respect to
conventional unit cells.

For more information on Monkhorst-Pack grids, see [[ngkpt]].
""",
),

Variable(
    abivarname="nimage",
    varset="rlx",
    vartype="integer",
    topics=['PIMD_compulsory', 'TransPath_compulsory'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Number of IMAGEs",
    added_in_version="before_v9",
    text=r"""
Give the number of images (or replicas) of the system, for which the forces
and stresses might be computed independently, in the context of the string
method, the genetic algorithm, hyperdynamics or Path-Integral Molecular
Dynamics depending on the value of [[imgmov]]). Related input variables:
[[dynimage]], [[npimage]], [[ntimimage]] and [[prtvolimg]].
Images might differ by the position of atoms in the unit cell, their velocity,
as well as by their cell geometry. The following input variables might be used
to define the images:

  * [[acell]]
  * [[amu]]
  * [[angdeg]]
  * [[dmatpawu]]
  * [[jpawu]]
  * [[mixalch]]
  * [[rprim]]
  * [[upawu]]
  * [[vel]]
  * [[vel_cell]]
  * [[xangst]]
  * [[xcart]]
  * [[xred]]

These input variables, non-modified, will be used to define the image with
index 1. For the image with the last index, the input file might specify the
values of such input variables, appended with "_lastimg", e.g.:

  * acell_lastimg
  * rprim_lastimg
  * xcart_lastimg
  * ...

By default, these values will be interpolated linearly to define values for
the other images, unless there exist specific values for some images, for
which the string "last" has to be replaced by the index of the image, e.g. for
the image number 4:

  * acell_4img
  * rprim_4img
  * xcart_4img
  * ...

It is notably possible to specify the starting point and the end point of the
path (of images), while specifying intermediate points.

It usually happen that the images do not have the same symmetries and space
group. ABINIT has not been designed to use different set of symmetries for
different images. ABINIT will use the symmetry and space group of the image
number 2, that is expected to have a low number of symmetries. This might lead
to erroneous calculations, in case some image has even less symmetry. By
contrast, there is no problem if some other image has more symmetries than
those of the second image.
""",
),

Variable(
    abivarname="nkpath",
    varset="basic",
    vartype="integer",
    topics=['k-points_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of K-points defining the PATH",
    added_in_version="before_v9",
    text=r"""
This variable is used to define the number of high-symmetry k-points in the
[[kptbounds]] array when [[kptopt]] > 0. Historically, [[kptbounds]] is used
in conjuction with a negative value of [[kptopt]] when performing a NSCF band
structure calculation. In this case, the number of k-points in kptbounds is
given by abs(kptopt) + 1. There are, however, other cases in which one has to
specify a k-path in the input file in order to activate some kind of post-
processing tool. Typical examples are the interpolation of the GW corrections
at the end of the sigma run or the interpolation of the KS eigenvalues along a
path at the end of the SCF run (see also [[einterp]]) In a nutshell, nkpath
replaces [[kptopt]] when we are not performing a NSCF calculation. Note that,
unlike [[kptopt]], nkpath represents the total number of points in the
[[kptbounds]] array.
""",
),

Variable(
    abivarname="nkpt",
    varset="basic",
    vartype="integer",
    topics=['k-points_useful'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[kptopt]] == 0': 1, 'defaultval': 0}),
    mnemonics="Number of K - Points",
    added_in_version="before_v9",
    text=r"""
If non-zero, [[nkpt]] gives the number of k points in the k point array
[[kpt]]. These points are used either to sample the Brillouin zone, or to
build a band structure along specified lines.

If [[nkpt]] is zero, the code deduces from other input variables (see the list
in the description of [[kptopt]]) the number of k points, which is possible
only when [[kptopt]]/=0. If [[kptopt]]/=0 and the input value of [[nkpt]]/=0,
then ABINIT will check that the number of k points generated from the other
input variables is exactly the same than [[nkpt]].

If [[kptopt]] is positive, [[nkpt]] must be coherent with the values of
[[kptrlatt]], [[nshiftk]] and [[shiftk]].
For ground state calculations, one should select the k point in the
irreducible Brillouin Zone (obtained by taking into account point symmetries
and the time-reversal symmetry).
For response function calculations, one should select k points in the full
Brillouin zone, if the wavevector of the perturbation does not vanish, or in a
half of the Brillouin Zone if q=0. The code will automatically decrease the
number of k points to the minimal set needed for each particular perturbation.

If [[kptopt]] is negative, [[nkpt]] will be the sum of the number of points on
the different lines of the band structure. For example, if [[kptopt]] = -3, one
will have three segments; supposing [[ndivk]] is 10 12 17, the total number of
k points of the circuit will be 10+12+17+1(for the final point)=40.
""",
),

Variable(
    abivarname="nkptgw",
    varset="gw",
    vartype="integer",
    topics=['SelfEnergy_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of K-PoinTs for GW corrections",
    requires="[[optdriver]] == 4",
    added_in_version="before_v9",
    text=r"""
[[nkptgw]] gives the number of k-points for which the GW calculation must be
done. It is used to dimension [[kptgw]].
""",
),

Variable(
    abivarname="nkpthf",
    varset="basic",
    vartype="integer",
    topics=['Hybrids_expert'],
    dimensions="scalar",
    mnemonics="Number of K - Points for (Hartree) Fock exact exchange",
    added_in_version="before_v9",
    text=r"""
[[nkpthf]] gives the number of k points used to sample the full Brillouin zone
for the Fock exact exchange contribution. It is obtained from the
specification of the wavefunction k point grid (see [[kptopt]]), possibly
downfolded as specified by [[fockdownsampling]].
""",
),

Variable(
    abivarname="nline",
    varset="gstate",
    vartype="integer",
    topics=['SCFControl_expert'],
    dimensions="scalar",
    defaultval=4,
    mnemonics="Number of LINE minimisations",
    added_in_version="before_v9",
    text=r"""
Gives maximum number of line minimizations allowed in preconditioned conjugate
gradient minimization for each band. The default, 4, is fine.
Special cases, with degeneracies or near-degeneracies of levels at the Fermi
energy may require a larger value of [[nline]] (5 or 6 ?). Line minimizations
will be stopped anyway when improvement gets small (governed by [[tolrde]]).
With the input variable [[nnsclo]], governs the convergence of the
wavefunctions for fixed potential.
Note that [[nline]] = 0 can be used to diagonalize the Hamiltonian matrix in the
subspace spanned by the input wavefunctions.
""",
),

Variable(
    abivarname="nloc_alg",
    varset="dev",
    vartype="integer",
    topics=['TuningSpeed_expert'],
    dimensions="scalar",
    defaultval=4,
    mnemonics="Non LOCal ALGorithm",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Allows to choose the algorithm for non-local operator application. On super-
scalar architectures, the default [[nloc_alg]] = 4 is the best.
More detailed explanations:

- [[nloc_alg]] = 2: Should be efficient on vector machines. It is indeed the
  fastest algorithm for the NEC, but actual tests on Fujitsu machine did not
  gave better performances than the other options.
- [[nloc_alg]] = 3: same as [[nloc_alg]] == 2, but the loop order is inverted.
- [[nloc_alg]] = 4: same as [[nloc_alg]] == 3, but maximal use of registers has
   been coded. This should be especially efficient on scalar and super-scalar
   machines. This has been confirmed by tests.

Note: internally, [[nloc_alg]] is stored in `dtset%nloalg(1)`. See also
[[nloc_mem]] for the tuning of the memory used in the non-local operator application.
""",
),

Variable(
    abivarname="nloc_mem",
    varset="dev",
    vartype="integer",
    topics=['TuningSpeed_expert'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[usepaw]] == 1': 2, 'defaultval': 1}),
    mnemonics="Non LOCal MEMOry",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Controls the memory use for the application of the non-local operator.
More detailed explanations:

- [[nloc_mem]] == 1: (k+G) vectors are not precomputed, in order to save memory space.
- [[nloc_mem]] == 2: (k+G) vectors are precomputed, once per k-point.
- [[nloc_mem]] == -1 or -2: Negative values of [[nloc_mem]] correspond
  positive ones, where the phase precomputation has been suppressed, in order to
  save memory space, as an array `double precision: ph3d(2,npw,[[natom]])` is
  saved (typically half the space needed for the wavefunctions at 1 k point -
  this corresponds to the silicon case). However, the computation of phases
  inside nonlop is somehow time-consuming.

Note: internally, sign([[nloc_mem]]) is stored in `dtset%nloalg(2)` and
abs([[nloc_mem]])-1 is stored in `dtset%nloalg(3)`. See also [[nloc_alg]] for the
algorithm for the non-local operator application.
""",
),

Variable(
    abivarname="nnos",
    varset="rlx",
    vartype="integer",
    topics=['PIMD_basic', 'MolecularDynamics_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of NOSe masses",
    added_in_version="before_v9",
    text=r"""
Gives the number of thermostats in the chain of oscillators
thermostats as proposed in [[cite:Martyna1996]]. The thermostat chains can be used either to perform Molecular Dynamics (MD) ([[ionmov]] = 13) or to perform Path Integral Molecular Dynamics
(PIMD) ([[imgmov]] = 13).
The mass of these thermostats is given by [[qmass]].
""",
),

Variable(
    abivarname="nnsclo",
    varset="dev",
    vartype="integer",
    topics=['SCFControl_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of Non-Self Consistent LOops",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Gives the maximum number of non-self-consistent loops of [[nline]] line
minimisations, in the SCF case (when [[iscf]] >0). In the case [[iscf]] <=0,
the number of non-self-consistent loops is determined by [[nstep]].

    * The Default value of 0 -- for standard plane-wave calculations -- corresponds
    to make the two first fixed potential determinations of wavefunctions have 2
    non-self consistent loops, and the next ones to have only 1 non-self
    consistent loop.

    * The Default value of 0 -- for wavelets calculations ([[usewvl]] = 1) --
    corresponds to make 2 steps with 3 non-self consistent loops, 2 steps with 2
    non-self consistent loops, then the next ones with 1 non-self consistent loop.

    * A negative value corresponds to make the abs([[nnsclo]]) first fixed potential determinations
    of wavefunctions have 5 non-self consistent loops, and the next ones to have only 1 non-self
    consistent loop.

""",
),

Variable(
    abivarname="nnsclohf",
    varset="dev",
    vartype="integer",
    topics=['Hybrids_expert'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[usefock]] == 1': 1, 'defaultval': 0}),
    mnemonics="Number of Non-Self Consistent LOops for (Hartree)-Fock exact exchange",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Gives the maximum number of loops with non-self-consistent occupied states
used to calculate Fock exact exchange, in the SCF case.
The Default value is 0 when [[usefock]] = 0. Default value is 1 when
[[usefock]] = 1 and correspond to update occupied wavefunctions at each self-consistent loop.
""",
),

Variable(
    abivarname="nobj",
    varset="geo",
    vartype="integer",
    topics=['AtomManipulator_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of OBJects",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Gives the number of 'objects' to be used by the atom manipulator in order to
find the full set of atoms. At present, only one or two objects can be
defined, identified as objects 'a' and 'b'.
Related variables for object 'a' are: [[objan]], [[objaat]], [[objarf]],
[[objatr]], [[objaro]], [[objaax]]. Related variables for object 'b' are:
[[objbn]], [[objbat]], [[objbrf]], [[objbtr]], [[objbro]], [[objbax]].

More detailed explanation: when the atom manipulator is used (i.e. when
[[nobj]] == 1 or [[nobj]] == 2), the code will be given a primitive set of atoms,
from which it will have to deduce the full set of atoms.
An object will be specified by the number of atoms it includes ([[objan]] or
[[objbn]] ), and the list of these atoms ([[objaat]] or [[objbat]] ).
Examples of physical realisation of an object can be a molecule, or a group of
atom to be repeated, or a part of a molecule to be rotated. The atom
amnipulator can indeed repeat these objects ([[objarf]] or [[objbrf]] ),
rotate them ([[objaro]] or [[objbro]] ) with respect to an axis ([[objaax]] or
[[objbax]] ), and translate them ([[objatr]] or [[objbtr]] ). After having
generated a geometry thanks to rotation, translation and repetition of
objects, it is possible to remove some atoms, in order to create vacancies
([[vacnum]] and [[vaclst]]). The number of atoms in the primitive set, those
that will be read from the input file, is specified by the variable [[natrd]].
It will be always smaller than the final number of atoms, given by the
variable [[natom]]. The code checks whether the primitive number of atoms plus
those obtained by the repetition operation is coherent with the variable
[[natom]], taking into account possible vacancies.
You should look at the other variables for more information.
Go to [[objan]], for example.
""",
),

Variable(
    abivarname="nomegasf",
    varset="gw",
    vartype="integer",
    topics=['Susceptibility_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of OMEGA to evaluate the Spectral Function",
    characteristics=['[[ENERGY]]'],
    requires="[[optdriver]] == 3 and [[spmeth]]!=0",
    added_in_version="before_v9",
    text=r"""
[[nomegasf]] defines the number of real frequencies used to describe the
spectral function associated to the irreducible polarizability
$\chi^{(0)}_{KS}$. The frequency mesh will cover the interval between 0 and
the maximum (positive) transition energy between occupied and empty states.
The delta function entering the expression defining the spectral function is
approximated using two different methods according to the value of the [[spmeth]] input variable.

It is important to notice that an accurate description of the imaginary part
of $\chi^{(0)}_{KS}$ requires an extremely dense frequency mesh. It should be
kept in mind, however, that the memory required grows fast with the value of [[nomegasf]].
""",
),

Variable(
    abivarname="nomegasi",
    varset="gw",
    vartype="integer",
    topics=['SelfEnergy_basic'],
    dimensions="scalar",
    defaultval=12,
    mnemonics="Number of OMEGA(S) along the Imaginary axis",
    requires="[[optdriver]] == 4 and [[gwcalctyp]] == 1",
    added_in_version="before_v9",
    text=r"""
[[nomegasi]] defines the number of frequency points used to sample the self-
energy along the imaginary axis. The frequency mesh is linear and covers the
interval between `omegasimin`=0.01 Hartree and [[omegasimax]].
""",
),

Variable(
    abivarname="nomegasrd",
    varset="gw",
    vartype="integer",
    topics=['SelfEnergy_useful'],
    dimensions="scalar",
    defaultval=9,
    mnemonics="Number of OMEGA to evaluate the Sigma Real axis Derivative",
    requires="[[optdriver]] == 4",
    added_in_version="before_v9",
    text=r"""
The number of real frequencies around the KS energy where the self-energy
$\Sigma$ is evaluated. From these values, the derivative of $\Sigma$ at the KS
energy is numerically estimated through linear interpolation.
""",
),

Variable(
    abivarname="nonlinear_info",
    varset="dfpt",
    vartype="integer",
    topics=['nonlinear_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Output NON-LINEAR INFOrmation",
    requires="[[optdriver]] == 5 and [[usepead]] == 0, or [[rf2_dkdk]]/=0 or [[rf2_dkde]]/=0",
    added_in_version="before_v9",
    text=r"""
Control the output of the non-linear implementation (only when [[usepead]] == 0).
The default value, [[nonlinear_info]] == 0 does nothing. If [[nonlinear_info]] == 1,
different contributions of 3rd derivatives of the energy are written in the
output file (non time consuming).

Higher values activate some internal tests for
checking the implementation correctness (time consuming, not useable in parallel).
If [[nonlinear_info]] == 2, same effect than 1 and tests are done in non-linear
([[optdriver]]==5 and [[usepead]] == 0).
If [[nonlinear_info]] == 3, same effect than 1 and tests are done in rf2_init
([[rf2_dkdk]]/=0 or [[rf2_dkde]]/=0).
If [[nonlinear_info]] == 4, same effect than 1 and tests are done in both non-linear and rf2_init.
A line containining "NOT PASSED" (and other information) is added to the output file
for each test that does not pass, otherwise nothing is printed. However, more information concerning
the tests is always printed in the **standard** output file.
""",
),

Variable(
    abivarname="normpawu",
    varset="dev",
    vartype="integer",
    topics=['DFT+U_expert'],
    dimensions=['[[ntypat]]'],
    defaultval=0,
    mnemonics="NORMalize atomic PAW+U projector",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Defines whether the atomic wave function (used as projectors in PAW+U) should
be renormalized to 1 within PAW sphere.

  * [[normpawu]] = 0: leave projector
  * [[normpawu]] = 1: renormalize
""",
),

Variable(
    abivarname="noseinert",
    varset="rlx",
    vartype="real",
    topics=['MolecularDynamics_useful'],
    dimensions="scalar",
    defaultval=100000,
    mnemonics="NOSE thermostat INERTia factor",
    requires="[[ionmov]] == 8",
    added_in_version="before_v9",
    text=r"""
Give the inertia factor WT of the Nose-Hoover thermostat (when [[ionmov]] = 8),
in atomic units of weight*length2, that is (electron mass)*(Bohr)2. The
equations of motion are: MI d2RI/dt2= FI - dX/dt MI dRI/dt and WT d2X/dt2=
Sum(I) MI (dRI/dt)2 - 3NkBT where I represent each nucleus, MI is the mass of
each nucleus (see [[amu]]), RI is the coordinate of each nucleus (see
[[xcart]]), dX/dt is a dynamical friction coefficient, and T is the
temperature of the thermostat (see [[mdtemp]]).
""",
),

Variable(
    abivarname="np_slk",
    varset="paral",
    vartype="integer",
    topics=['parallelism_expert'],
    dimensions="scalar",
    defaultval=1000000,
    mnemonics="Number of mpi Processors used for ScaLapacK calls",
    characteristics=['[[DEVELOP]]'],
    requires="[[optdriver]] == 1 and [[paral_kgb]] == 1 (Ground-state calculations with LOBPCG algorithm)",
    added_in_version="before_v9",
    text=r"""
When using Scalapack (or any similar Matrix Algebra library), the efficiency
of the eigenproblem resolution saturates as the number of CPU cores increases.
It is better to use a smaller number of CPU cores for the LINALG calls.
This maximum number of cores can be set with [[np_slk]].
A large number for [[np_slk]] (i.e. 1000000) means that all cores are used for
the Linear Algebra calls.
np_slk must divide the number of processors involved in diagonalizations
([[npband]] $\times$ [[npfft]] $\times$ [[npspinor]]).
Note (bef v8.8): an optimal value for this parameter can be automatically found by using
the [[autoparal]] input keyword.
Note (since v8.8 and only for LOBPCG ([[wfoptalg]] == 114) with [[paral_kgb]] = 1):
* If set to ***0*** then scalapack is disabled.
* If set to its ***default value***, then abinit uses between 2 and [[npband]]*[[npfft]]*[[npspinor]] cpus according to the system size (default auto behaviour).
  See [[slk_rankpp]] for more customization.
* If set to a number ***>1** then forces abinit to use exactly this number of cpus.
  Due to legacy behaviour, although it is not mandatory in theory, this value *must* divide [[npband]]*[[npfft]]*[[npspinor]].
See also [[slk_rankpp]] for a better tuning instead of this variable.
""",
),

Variable(
    abivarname="npband",
    varset="paral",
    vartype="integer",
    topics=['parallelism_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Number of Processors at the BAND level",
    requires="[[paral_kgb]] == 1",
    added_in_version="before_v9",
    text=r"""
Relevant only for the band/FFT parallelisation (see the [[paral_kgb]] input variable).
[[npband]] gives the number of processors among which the work load over the
band level is shared. [[npband]], [[npfft]], [[npkpt]] and [[npspinor]] are
combined to give the total number of processors (nproc) working on the
band/FFT/k-point parallelisation.
See [[npfft]], [[npkpt]], [[npspinor]] and [[paral_kgb]] for the additional
information on the use of band/FFT/k-point parallelisation. [[npband]] has to
be a divisor or equal to [[nband]]
Note: an optimal value for this parameter can be automatically found by using
the [[autoparal]] input keyword.
""",
),

Variable(
    abivarname="npfft",
    varset="paral",
    vartype="integer",
    topics=['parallelism_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Number of Processors at the FFT level",
    requires="[[paral_kgb]] == 1",
    added_in_version="before_v9",
    text=r"""
Relevant only for the band/FFT/k-point parallelisation (see the [[paral_kgb]]
input variable).
[[npfft]] gives the number of processors among which the work load over the
FFT level is shared. [[npfft]], [[npkpt]], [[npband]] and [[npspinor]] are
combined to give the total number of processors (nproc) working on the
band/FFT/k-point parallelisation.
See [[npband]], [[npkpt]], [[npspinor]], and [[paral_kgb]] for the additional
information on the use of band/FFT/k-point parallelisation.

Note: [[ngfft]] is automatically adjusted to [[npfft]]. If the number of
processor is changed from a calculation to another one, [[npfft]] may change,
and then [[ngfft]] also.
Note: an optimal value for this parameter can be automatically found by using
the [[autoparal]] input keyword.
""",
),

Variable(
    abivarname="nphf",
    varset="paral",
    vartype="integer",
    topics=['Hybrids_useful', 'parallelism_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Number of Processors for (Hartree)-Fock exact exchange",
    added_in_version="before_v9",
    text=r"""
Relevant only for the k-point/fock parallelisation (option [[paral_kgb]] input
variable).
[[nphf]] gives the number of processors among which the work load over the
occupied states level is shared. [[nphf]] and [[npkpt]] are combined to give
the total number of processors (nproc) working on the parallelisation.

Note: [[nphf]] should be a divisor or equal to the number of k-point times
the number of bands for exact exchange ([[nkpthf]] $\times$ [[nbandhf]]) in order to
have the better load-balancing and efficiency.
""",
),

Variable(
    abivarname="npimage",
    varset="paral",
    vartype="integer",
    topics=['parallelism_useful', 'PIMD_useful', 'TransPath_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Number of Processors at the IMAGE level",
    added_in_version="before_v9",
    text=r"""
Relevant only when sets of images are activated (see [[imgmov]] and [[nimage]]
).
[[npimage]] gives the number of processors among which the work load over the
image level is shared. It is compatible with all other parallelization levels
available for ground-state calculations.
Note: an optimal value for this parameter can be automatically found by using
the [[autoparal]] input keyword.

See [[paral_kgb]], [[npkpt]], [[npband]], [[npfft]] and [[npspinor]] for the
additional information on the use of k-point/band/FFT parallelisation.
""",
),

Variable(
    abivarname="npkpt",
    varset="paral",
    vartype="integer",
    topics=['parallelism_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Number of Processors at the K-Point Level",
    requires="[[paral_kgb]] == 1",
    added_in_version="before_v9",
    text=r"""
Relevant only for the band/FFT/k-point parallelisation (see the [[paral_kgb]]
input variable).
[[npkpt]] gives the number of processors among which the work load over the
k-point/spin-component level is shared. [[npkpt]], [[npfft]], [[npband]] and
[[npspinor]] are combined to give the total number of processors (nproc)
working on the band/FFT/k-point parallelisation.
See [[npband]], [[npfft]], [[npspinor]] and [[paral_kgb]] for the additional
information on the use of band/FFT/k-point parallelisation.

[[npkpt]] should be a divisor or equal to with the number of k-point/spin-
components ([[nkpt]] $\times$ [[nsppol]]) in order to have the better load-balancing
and efficiency.
Note: an optimal value for this parameter can be automatically found by using
the [[autoparal]] input keyword.
""",
),

Variable(
    abivarname="nppert",
    varset="paral",
    vartype="integer",
    topics=['parallelism_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Number of Processors at the PERTurbation level",
    requires="[[paral_rf]] == 1",
    added_in_version="before_v9",
    text=r"""
This parameter is used in connection to the parallelization over
perturbations (see [[paral_rf]] ), for a linear response calculation.
[[nppert]] gives the number of processors among which the work load over the
perturbation level is shared. It can even be specified separately for each
dataset.
""",
),

Variable(
    abivarname="npsp",
    varset="gstate",
    vartype="integer",
    topics=['AtomTypes_useful', 'PseudosPAW_expert'],
    dimensions="scalar",
    defaultval="[[ntypat]]",
    mnemonics="Number of PSeudoPotentials",
    characteristics=['[[NO_MULTI]]'],
    added_in_version="before_v9",
    text=r"""
Usually, the number of pseudopotentials to be read is equal to the number of
type of atoms. However, in the case an alchemical mixing of pseudopotential is
to be used, often the number of pseudopotentials to be read will not equal the
number of types of atoms.

Alchemical pseudopotentials will be present when [[ntypalch]] is non-zero. See
[[ntypalch]] to understand how to use alchemical potentials in ABINIT. The
input variables ([[ntypalch]], [[algalch]], [[mixalch]]) are active, and
generate alchemical potentials from the available pseudopotentials. Also, the
inner variables ([[ntyppure]], [[npspalch]]) become active. See these input
variables, especially [[mixalch]], to understand how to use alchemical
potentials in ABINIT.
""",
),

Variable(
    abivarname="npspalch",
    varset="gstate",
    vartype="integer",
    topics=['AtomTypes_internal'],
    dimensions="scalar",
    defaultval="[[npsp]]-[[ntyppure]]",
    mnemonics='Number of PSeudoPotentials that are "ALCHemical"',
    characteristics=['[[INTERNAL_ONLY]]'],
    requires="[[ntypalch]]/=0",
    added_in_version="before_v9",
    text=r"""
Gives the number of pseudopotentials that are used for alchemical mixing (when [[ntypalch]] is non-zero):

[[npspalch]] = [[npsp]]-[[ntyppure]]
""",
),

Variable(
    abivarname="npspinor",
    varset="paral",
    vartype="integer",
    topics=['parallelism_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Number of Processors at the SPINOR level",
    requires="[[paral_kgb]] == 1",
    added_in_version="before_v9",
    text=r"""
Can be 1 or 2 (if [[nspinor]] = 2).
Relevant only for the band/FFT/k-point parallelisation (see the [[paral_kgb]]
input variable).
[[npspinor]] gives the number of processors among which the work load over the
spinorial components of wave-functions is shared. [[npspinor]], [[npfft]],
[[npband]] and [[npkpt]] are combined to give the total number of processors
(nproc) working on the band/FFT/k-point parallelisation.
Note: an optimal value for this parameter can be automatically found by using
the [[autoparal]] input keyword.

See [[npkpt]], [[npband]], [[npfft]], and [[paral_kgb]] for the additional
information on the use of band/FFT/k-point parallelisation.
""",
),

Variable(
    abivarname="npulayit",
    varset="dev",
    vartype="integer",
    topics=['SCFAlgorithms_useful'],
    dimensions="scalar",
    defaultval=7,
    mnemonics="Number of PULAY ITerations for SC mixing",
    characteristics=['[[DEVELOP]]'],
    requires="[[iscf]] in [7,17]",
    added_in_version="before_v9",
    text=r"""
Gives the number of previous iterations involved in Pulay mixing (mixing
during electronic SC iterations).
""",
),

Variable(
    abivarname="npvel",
    varset="gw",
    vartype="integer",
    topics=['RandStopPow_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of Particle VELocities",
    requires="[[optdriver]] == 3",
    added_in_version="before_v9",
    text=r"""
In the context of the electronic stopping power of impinging ion in matter,
[[npvel]] sets the number of the ion velocities to be calculated via linear response.
When [[npvel]] = 0, no stopping power calculation is performed.
The direction and the velocity maximum are set with the input variable
[[pvelmax]]. Note that the results are output for a $Z=1$ impinging ion, i.e. a proton.
""",
),

Variable(
    abivarname="npweps",
    varset="internal",
    vartype="integer",
    topics=['Susceptibility_internal'],
    dimensions="scalar",
    mnemonics="Number of PlaneWaves for EPSilon (the dielectric matrix)",
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
[[npweps]] determines the size of the planewave set used to represent the
independent-particle susceptibility $\chi^{(0)}_{KS}$, the dielectric matrix
$\epsilon$ and its inverse.
It is an internal variable, determined from [[ecuteps]].
""",
),

Variable(
    abivarname="npwkss",
    varset="gw",
    vartype="integer",
    topics=['Susceptibility_expert', 'SelfEnergy_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of PlaneWaves in the KSS file",
    added_in_version="before_v9",
    text=r"""
This input variable is used for the preparation of a GW calculation: the GS
run (where [[optdriver]] = 1 and [[nbandkss]]/=0) should be followed with a run
where [[optdriver]] = 3. Also, if [[nbandkss]] = 0, no use of [[npwkss]].

[[npwkss]] defines the number of planewave components of the Kohn-Sham states
to build the Hamiltonian, in the routine outkss.F90, and so, the size of the
matrix, the size of eigenvectors, and the number of available states, to be
stored in the abo_KSS file. If it is set to 0, then, the planewave basis set
defined by the usual Ground State input variable [[ecut]] is used to generate
the superset of all planewaves used for all k points. Note that this (large)
planewave basis is the same for all k points.

!!! warning

    For the time being, [[istwfk]] must be 1 for all the k points.
""",
),

Variable(
    abivarname="npwsigx",
    varset="internal",
    vartype="integer",
    topics=['SelfEnergy_internal'],
    dimensions="scalar",
    mnemonics="Number of PlaneWaves for SIGma eXchange",
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
[[npwsigx]] determines the cut-off energy of the planewave set used to
generate the exchange part of the self-energy operator.
It is an internal variable, determined from [[ecutsigx]].
""",
),

Variable(
    abivarname="npwwfn",
    varset="internal",
    vartype="integer",
    topics=['Susceptibility_internal', 'SelfEnergy_internal'],
    dimensions="scalar",
    mnemonics="Number of PlaneWaves for WaveFunctioNs",
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
[[npwwfn]] is the size of the planewave set used to represent the
wavefunctions in the formula that generates the independent-particle
susceptibility $\chi^{(0)}_{KS}$. It is an internal variable, determined from [[ecutwfn]].
""",
),

Variable(
    abivarname="nqpt",
    varset="gstate",
    vartype="integer",
    topics=['q-points_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of Q - POINTs",
    added_in_version="before_v9",
    text=r"""
Determines whether one q point must be read (See the variable [[qptn]]).
Can be either 0 or 1.

If 1 and used in ground-state calculation, a global shift of all the k-points
is applied, to give calculation at k+q. In this case, the output wavefunction
will be appended by _WFQ instead of _WFK
Also, if 1 and a RF calculation is done, defines the wavevector of the perturbation.

For further information about the *files file*, consult the [[help:abinit#files-file]].
""",
),

Variable(
    abivarname="nqptdm",
    varset="gw",
    vartype="integer",
    topics=['Susceptibility_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of Q-PoinTs for the Dielectric Matrix",
    requires="[[optdriver]] == 3",
    added_in_version="before_v9",
    text=r"""
If [[nqptdm]] is equal to 0, the set of q points for computing the dielectric
matrix is determined automatically considering all the possible differences
between the k-points contained in the _KSS file. When [[nqptdm]] is non-zero,
the list of q points is read from [[qptdm]]. This allows one to split the big
calculation of all the dielectric matrices into smaller calculations that can
be performed independently. The _SCR files generated in different runs can be
merged thanks to the **Mrgscr** utility. If [[nqptdm]] is equal to -1, the
code reports the list of q-points in the log file (YAML format) and then stops.
""",
),

Variable(
    abivarname="nscforder",
    varset="dev",
    vartype="integer",
    topics=['Coulomb_expert'],
    dimensions="scalar",
    defaultval=16,
    mnemonics="Nth - SCaling Function ORDER",
    added_in_version="before_v9",
    text=r"""
This variable controls the order of used scaling functions when the Hartree
potential is computed using the Poisson solver (see [[icoulomb]] input
variable). This variable is of seldom use since the default value is large
enough. Nonetheless, possible values are 8, 14, 16, 20, 24, 30, 40, 50, 60,
100. Values greater than 20 are included in ABINIT for test purposes only.
""",
),

Variable(
    abivarname="nshiftk",
    varset="basic",
    vartype="integer",
    topics=['k-points_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Number of SHIFTs for K point grids",
    added_in_version="before_v9",
    text=r"""
This parameter gives the number of shifted grids to be used concurrently to
generate the full grid of k points. It can be used with primitive grids
defined either from [[ngkpt]] or [[kptrlatt]]. The maximum allowed value of
[[nshiftk]] is 8. The values of the shifts are given by [[shiftk]].
""",
),

Variable(
    abivarname="nshiftq",
    varset="gstate",
    vartype="integer",
    topics=['q-points_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Number of SHIFTs for Q point grids",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
This parameter gives the number of shifted grids to be used concurrently to
generate the full grid of q points. It can be used with primitive grids
defined either from [[ngqpt]] or [[qptrlatt]]. The maximum allowed value of
[[nshiftq]] is 8. The values of the shifts are given by [[shiftq]].
""",
),

Variable(
    abivarname="nspden",
    varset="gstate",
    vartype="integer",
    topics=['spinpolarisation_basic'],
    dimensions="scalar",
    defaultval="[[nsppol]]",
    mnemonics="Number of SPin-DENsity components",
    added_in_version="before_v9",
    text=r"""
If [[nspden]] = 1, no spin-magnetization the density matrix is diagonal, with
same values spin-up and spin-down (compatible with [[nsppol]] = 1 only, for both
[[nspinor]] = 1 or 2)

If [[nspden]] = 2, scalar magnetization (the axis is arbitrarily fixed in the z
direction) the density matrix is diagonal, with different values for spin-up
and spin-down (compatible with [[nspinor]] = 1, either with [[nsppol]] = 2
-general collinear magnetization- or [[nsppol]] = 1 -antiferromagnetism)

If [[nspden]] = 4, vector magnetization: the density matrix is full, with
allowed x, y and z magnetization (useful only with [[nspinor]] = 2 and
[[nsppol]] = 1, either because there is spin-orbit without time-reversal
symmetry - and thus spontaneous magnetization, or with spin-orbit, if one
allows for spontaneous non-collinear magnetism). Available for
response functions [[cite:Ricci2019]]. Also note that, with [[nspden]] = 4, time-reversal symmetry
is not taken into account (at present; this has to be checked) and thus
[[kptopt]] has to be different from 1 or 2.

The default ([[nspden]] = [[nsppol]]) does not suit the case of vector magnetization.
Note that the choice of [[nspden]] has an influence on the treatment of symmetries. See [[symafm]].
""",
),

Variable(
    abivarname="nspinor",
    varset="gstate",
    vartype="integer",
    topics=['spinpolarisation_basic'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[pawspnorb]] == 1': 2, 'defaultval': 1}),
    mnemonics="Number of SPINORial components of the wavefunctions",
    added_in_version="before_v9",
    text=r"""
If [[nspinor]] = 1, usual case: scalar wavefunction (compatible with
([[nsppol]] = 1, [[nspden]] = 1) as well as ([[nsppol]] = 2, [[nspden]] = 2) )

If [[nspinor]] = 2, the wavefunction is a spinor (compatible with [[nsppol]] = 1,
with [[nspden]] = 1 or 4, but not with [[nsppol]] = 2)

When [[nspinor]] is 2, the values of [[istwfk]] are automatically set to 1.
Also, the number of bands, for each k-point, should be even.
""",
),

Variable(
    abivarname="nsppol",
    varset="basic",
    vartype="integer",
    topics=['spinpolarisation_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Number of SPin POLarization",
    added_in_version="before_v9",
    text=r"""
Give the number of INDEPENDENT spin polarisations, for which there are non-
related wavefunctions. Can take the values 1 or 2.

If [[nsppol]] = 1, one has an unpolarized calculation ([[nspinor]] = 1,
[[nspden]] = 1) or an antiferromagnetic system ([[nspinor]] = 1, [[nspden]] = 2), or
a calculation in which spin up and spin down cannot be disentangled
([[nspinor]] = 2), that is, either non-collinear magnetism or presence of spin-
orbit coupling, for which one needs spinor wavefunctions.

If [[nsppol]] = 2, one has a spin-polarized (collinear) calculation with
separate and different wavefunctions for up and down spin electrons for each
band and k point. Compatible only with [[nspinor]] == 1, [[nspden]] == 2. If
[[nsppol]] = 2, one usually uses a metallic value for [[occopt]], in order to
let ABINIT find the magnetization. On the contrary, if [[occopt]] == 1 is used,
the user has to impose the magnetization, using [[spinmagntarget]], except for
the case of a single isolated Hydrogen atom.

In the present status of development, with [[nsppol]] = 1, all values of [[ixc]]
are allowed, while with [[nsppol]] = 2, some values of [[ixc]] might not be
allowed (e.g. 2, 3, 4, 5, 6, 20, 21, 22 are not allowed).

See also the input variable [[nspden]] for the components of the density
matrix with respect to the spin-polarization.
""",
),

Variable(
    abivarname="nstep",
    varset="basic",
    vartype="integer",
    topics=['SCFControl_basic'],
    dimensions="scalar",
    defaultval=30,
    mnemonics="Number of (non-)self-consistent field STEPS",
    added_in_version="before_v9",
    text=r"""
Gives the maximum number of cycles (or "iterations") in a SCF or non-SCF run.
Full convergence from random numbers is usually achieved in 12-20 SCF
iterations. Each can take from minutes to hours. In certain difficult cases,
usually related to a small or zero band gap or magnetism, convergence
performance may be much worse. When the convergence tolerance [[tolwfr]] on
the wavefunctions is satisfied, iterations will stop, so for well converged
calculations you should set [[nstep]] to a value larger than you think will be
needed for full convergence, e.g. if using 20 steps usually converges the
system, set [[nstep]] to 30.
For non-self-consistent runs ( [[iscf]] < 0) nstep governs the number of
cycles of convergence for the wavefunctions for a fixed density and Hamiltonian.

NOTE that a choice of [[nstep]] = 0 is permitted; this will either read
wavefunctions from disk (with [[irdwfk]] = 1 or [[irdwfq]] = 1, or non-zero
[[getwfk]] or [[getwfq]] in the case of multi-dataset) and compute the
density, the total energy and stop, or else (with all of the above vanishing)
will initialize randomly the wavefunctions and compute the resulting density
and total energy. This is provided for testing purposes.
Also NOTE that [[nstep]] = 0 with [[irdwfk]] = 1 will exactly give the same result
as the previous run only if the latter is done with [[iscf]]<10 (potential mixing).
One can output the density by using [[prtden]].
The forces and stress tensor are computed with [[nstep]] = 0.
""",
),

Variable(
    abivarname="nsym",
    varset="basic",
    vartype="integer",
    topics=['crystal_useful', 'GW_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of SYMmetry operations",
    added_in_version="before_v9",
    text=r"""
Gives number of space group symmetries to be applied in this problem.
Symmetries will be input in array "[[symrel]]" and (nonsymmorphic)
translations vectors will be input in array "[[tnons]]". If there is no
symmetry in the problem then set [[nsym]] to 1, because the identity is still
a symmetry.
In case of a RF calculation, the code is able to use the symmetries of the
system to decrease the number of perturbations to be calculated, and to
decrease of the number of special k points to be used for the sampling of the
Brillouin zone. After the response to the perturbations have been calculated,
the symmetries are used to generate as many as possible elements of the 2DTE
from those already computed.

**Symmetry finder mode** (Default mode).
If [[nsym]] is 0, all the atomic coordinates must be explicitly given (one
cannot use the atom manipulator neither the smart symmetrizer): the code will
then find automatically the symmetry operations that leave the lattice and
each atomic sublattice invariant. It also checks whether the cell is primitive
(see [[chkprim]]).
Note that the tolerance on symmetric atomic positions and lattice is rather
stringent: for a symmetry operation to be admitted, the lattice and atomic
positions must map on themselves within 1.0e-8.

The user is allowed to set up systems with non-primitive unit cells (i.e.
conventional FCC or BCC cells, or supercells without any distortion). In this
case, pure translations will be identified as symmetries of the system by the
symmetry finder. Then, the combined "pure translation + usual rotation and
inversion" symmetry operations can be very numerous. For example, a
conventional FCC cell has 192 symmetry operations, instead of the 48 ones of
the primitive cell. A maximum limit of 384 symmetry operations is hard-coded.
This corresponds to the maximum number of symmetry operations of a 2x2x2
undistorted supercell. Going beyond that number will make the code stop very
rapidly. If you want nevertheless, for testing purposes, to treat a larger
number of symmetries, change [[maxnsym]].

For GW calculation, the user might want to select only the symmetry operations
whose non-symmorphic translation vector [[tnons]] is zero. This can be done
with the help of the input variable [[symmorphi]]
""",
),

Variable(
    abivarname="ntime",
    varset="rlx",
    vartype="integer",
    topics=['MolecularDynamics_basic', 'GeoOpt_basic'],
    dimensions="scalar",
    defaultval="0 if ionmvov == 0, set to 1000 if ionvmov != 0 and imgmov != 0 and the variable is not specified.",
    mnemonics="Number of TIME steps",
    added_in_version="before_v9",
    text=r"""
Gives the maximum number of molecular dynamics time steps or structural
optimization steps to be done if [[ionmov]] is non-zero.
Starting with Abinit9, ntime is automatically set to 1000 if [[ionmov]] is non-zero,
[[ntimimage]] is zero and [[ntime]] is not specified in the input file.
Users are encouraged to pass a **timelimit** to Abinit using the command line and the syntax:

        abinit --timelimit hours:minutes:seconds

so that the code will try to stop smoothly before the timelimit and produce the DEN and the WFK files
that may be used to restart the calculation.

Note that at the present the option [[ionmov]] = 1 is initialized with four
Runge-Kutta steps which costs some overhead in the startup. By contrast, the
initialisation of other [[ionmov]] values is only one SCF call.
Note that **ntime** is ignored if [[ionmov]] = 0.
""",
),

Variable(
    abivarname="ntimimage",
    varset="rlx",
    vartype="integer",
    topics=['PIMD_compulsory', 'TransPath_compulsory'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Number of TIME steps for IMAGE propagation",
    added_in_version="before_v9",
    text=r"""
Gives the maximal number of molecular dynamics time steps or structural
optimization steps to be done for the set of images, referred to as 'image-
timesteps'. At each image-timestep, all the images are propagated
simultaneously, each according to the algorithm determined by [[imgmov]] and
the usual accompanying input variables, and then the next positions and
velocities for each image are determined from the set of results obtained for all images.
""",
),

Variable(
    abivarname="ntypalch",
    varset="gstate",
    vartype="integer",
    topics=['AtomTypes_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics='Number of TYPe of atoms that are "ALCHemical"',
    added_in_version="before_v9",
    text=r"""
Used for the generation of alchemical pseudopotentials: when **ntypalch** is
non-zero, alchemical mixing will be used.

Among the [[ntypat]] types of atoms, the last **ntypalch** will be
"alchemical" pseudoatoms, while only the first [[ntyppure]] will be uniquely
associated with a pseudopotential (the [[ntyppure]] first of these, actually).
The **ntypalch** types of alchemical pseudoatoms are to be made from the
remaining [[npspalch]] pseudopotentials.

In this case, the input variables [[algalch]], [[mixalch]] are active, and
generate alchemical potentials from the available pseudopotentials. See these
input variables, especially [[mixalch]], to understand how to use alchemical
potentials in ABINIT.
""",
),

Variable(
    abivarname="ntypat",
    varset="basic",
    vartype="integer",
    topics=['AtomTypes_compulsory', 'PseudosPAW_compulsory', 'crystal_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Number of TYPes of AToms",
    characteristics=['[[NO_MULTI]]'],
    added_in_version="before_v9",
    text=r"""
Gives the number of types of atoms in the unit cell.
E.g. for a homopolar system (e.g. pure Si) **ntypat** is 1.

The code tries to read the same number of pseudopotential files. The first
pseudopotential is assigned type number 1, and so on...

There is an exception in the case of alchemical mixing of potentials, for
which there is a different number of pseudopotentials atomic types. See [[mixalch]].
""",
),

Variable(
    abivarname="ntyppure",
    varset="gstate",
    vartype="integer",
    topics=['AtomTypes_internal'],
    dimensions="scalar",
    defaultval="[[ntypat]]-[[ntypalch]]",
    mnemonics='Number of TYPe of atoms that are "PURE"',
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Gives the number of type of atoms that are "pure" when alchemical mixing is
used ([[ntypalch]] /= 0):

[[ntyppure]] = [[ntypat]] - [[ntypalch]]
""",
),

Variable(
    abivarname="nucdipmom",
    varset="gstate",
    vartype="real",
    topics=['MagField_expert'],
    dimensions=[3, '[[natom]]'],
    defaultval=0.0,
    mnemonics="NUClear DIPole MOMents",
    requires="[[usepaw]] = 1; [[pawcpxocc]] = 2; [[kptopt]] > 2",
    added_in_version="before_v9",
    text=r"""
Places an array of nuclear magnetic dipole moments on the atomic
positions, useful for computing the magnetization in the presence of
nuclear dipoles and thus the chemical shielding by the converse method
[[cite:Thonhauser2009]]. The presence of these dipoles breaks time
reversal symmetry and lowers the overall spatial symmetry.  The dipole
moment values are entered in atomic units. For reference, note that
one Bohr magneton has value $1/2$ in atomic units, while one nuclear
Bohr magneton has value $2.7321\times 10^{-4}$ in atomic units.
""",
),

Variable(
    abivarname="nwfshist",
    varset="gstate",
    vartype="integer",
    topics=['Wavelets_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of WaveFunctionS HISTory",
    added_in_version="before_v9",
    text=r"""
In the wavelet basis set, the ground state is found by direct minimisation.
The algorithm used can be either the steepest descent or the DIIS (Direct Inversion of Iteration Space).
When [[nwfshist]] = 0, the steepest descent is used ( _i.e._ there is no history storage of the previous iterations).

If [[nwfshist]] is strictly positive, a DIIS is used. A typical value is 6. Using
a DIIS increases the memory required by the program since N previous
wavefunctions are stored during the electronic minimisation.
""",
),

Variable(
    abivarname="nzchempot",
    varset="geo",
    vartype="integer",
    topics=['Artificial_expert'],
    dimensions="scalar",
    mnemonics="Number of Z reduced coordinates that define the spatial CHEMical POTential",
    added_in_version="before_v9",
    text=r"""
Defines the number of z reduced coordinates that defines the spatially varying
chemical potential. See the input variable [[chempot]], of which [[nzchempot]]
is the second dimension.
""",
),

Variable(
    abivarname="objaat",
    varset="geo",
    vartype="integer",
    topics=['AtomManipulator_basic'],
    dimensions=['[[objan]]'],
    mnemonics="OBJect A: list of AToms",
    characteristics=['[[INPUT_ONLY]]'],
    requires="'[[nobj]] == 1'",
    added_in_version="before_v9",
    text=r"""
Gives the list of atoms in object a. This list is specified by giving, for
each atom, its index in the list of coordinates ([[xred]], [[xangst]] or
[[xcart]]), that also corresponds to a type of atom (given by the array type).
These objects can be thought as molecules, or groups of atoms, or parts of
molecules, to be repeated, rotated and translated to generate the full set of atoms.
Look at [[objarf]] for further explanations.
""",
),

Variable(
    abivarname="objaax",
    varset="geo",
    vartype="real",
    topics=['AtomManipulator_useful'],
    dimensions=[6],
    mnemonics="OBJect A: AXis",
    characteristics=['[[INPUT_ONLY]]', '[[LENGTH]]'],
    commentdefault="""[[objaax]] must be provided if ([[nobj]] == 1 and one component of [[objaro]] != 0). Moreover,
[[objaax]] AND [[objbax]] must be provided if ( [[nobj]] == 2 and one component of [[objbro]] != 0 ).""",
    added_in_version="before_v9",
    text=r"""
Gives, for each object, the cartesian coordinates of two points (first point:
[[objaax]](1:3) second point: [[objaax]](4:6). By default, given in Bohr
atomic units (1 Bohr=0.5291772108 Angstroms), although Angstrom can be
specified, if preferred, since these variables have the [[LENGTH]] characteristics.
The two points define an axis of rotation of the corresponding object.
Note that the rotation of the object is done BEFORE the object is translated.
The sign of the rotation angle is positive if the object is to be rotated
clockwise when looking to it along the axis, from point 1 (coordinates 1:3)
toward point 2 (coordinates 4:6).
""",
),

Variable(
    abivarname="objan",
    varset="geo",
    vartype="integer",
    topics=['AtomManipulator_basic'],
    dimensions="scalar",
    mnemonics="OBJect A: Number of atoms",
    characteristics=['[[INPUT_ONLY]]'],
    commentdefault=""" [[objan]] MUST be provided if [[nobj]] == 1.
 [[objan]] and [[objbn]] MUST be provided if [[nobj]] == 2.""",
    added_in_version="before_v9",
    text=r"""
Gives the number of atoms in object a. The list of atoms is given by the variables [[objaat]].
""",
),

Variable(
    abivarname="objarf",
    varset="geo",
    vartype="integer",
    topics=['AtomManipulator_basic'],
    dimensions=[3],
    defaultval=[1, 1, 1],
    mnemonics="OBJect A: Repetition Factors",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Gives three repetition factors of the objects a.
This gives the opportunity to generate a three-dimensional set of repeated
objects, although a simple one-dimensional repetition will be easily obtained
through the specification of `nrep 1 1`
where `nrep` is the 1D repetition factor.
The initial rotation and translation of the object, as well as the increment
of rotation or translation from one object to the next are specified by the
variable [[objaro]].
Note that the atom manipulator will generate the full set of atoms from the
primitive set of atoms using the following order: it will process each atom
in the primitive list one by one, determine whether it belongs to either
object a or object b, and then repeat it taking into account the proper
rotation and translation, with the fastest varying repetition factor being the
first, then the second, then the third.
In the final list of atoms, one will first find the atoms generated from atom
1 in the primitive list, then those generated from atom 2 in the primitive
list, and so on.
If the atom manipulator is only used to rotate or translate an object, without
repeating it, simply use `1 1 1`, which is also the Default value.
""",
),

Variable(
    abivarname="objaro",
    varset="geo",
    vartype="real",
    topics=['AtomManipulator_useful'],
    dimensions=[4],
    defaultval=MultipleValue(number=4, value=0.0),
    mnemonics="OBJect A: ROtations",
    characteristics=['[[INPUT_ONLY]]'],
    commentdefault="(no rotation)",
    added_in_version="before_v9",
    text=r"""
Give, for each object, the angles of rotation in degrees to be applied to the
corresponding object.
The rotation is applied before the translation, and the axis is defined by the
variable [[objaax]]. See the latter variable for the
definition of the sign of the rotation.
The first component [[objaro]](1) gives the angle of
rotation to be applied to the first instance of the object. The second, third
or fourth component (resp.) gives the increment of rotation angle from one
instance to the next instance, defined by the first, second or third
repetition factor (resp.). This allows one to generate 3D arrays of molecules
with different rotation angles.
""",
),

Variable(
    abivarname="objatr",
    varset="geo",
    vartype="real",
    topics=['AtomManipulator_useful'],
    dimensions=[12],
    defaultval=MultipleValue(number=12, value=0.0),
    mnemonics="OBJect A: TRanslations",
    characteristics=['[[INPUT_ONLY]]', '[[LENGTH]]'],
    commentdefault="(no translation)",
    added_in_version="before_v9",
    text=r"""
Give, for each object, the vectors of translations, in cartesian coordinates,
to be applied to the corresponding object. By default, given in Bohr atomic
units (1 Bohr=0.5291772108 Angstroms), although Angstrom can be specified, if
preferred, since these variables have the [[LENGTH]] characteristics.
The translation is applied after the rotation.
The first vector [[objatr]](3,1) gives the translation to
be applied to the first instance of the object. The second, third or fourth
component (resp.) gives the increment of translation from one instance to the
next instance, defined by the first, second or third repetition factor (resp.).
This allows one to generate 3D arrays of molecules.
In general, when the objects are repeated, a translation vector must be given,
since otherwise, the repeated objects pack in the same region of space. As an
exception, one can have a set of molecules regularly spaced on a circle, in
which case, only rotations are needed.
Not present in the dtset array (no internal).
""",
),

Variable(
    abivarname="objbat",
    varset="geo",
    vartype="integer",
    topics=['AtomManipulator_useful'],
    dimensions=['[[objbn]]'],
    mnemonics="OBJect B: list of AToms",
    characteristics=['[[INPUT_ONLY]]'],
    requires="[[nobj]] == 2",
    added_in_version="before_v9",
    text=r"""
Gives the list of atoms in object b. This list is specified by giving, for
each atom, its index in the list of coordinates ([[xred]], [[xangst]] or
[[xcart]]), that also corresponds to a type of atom (given by the array type).
These objects can be thought as molecules, or groups of atoms, or parts of
molecules, to be repeated, rotated and translated to generate the full set of
atoms.
Look at [[objbrf]] for further explanations.
""",
),

Variable(
    abivarname="objbax",
    varset="geo",
    vartype="real",
    topics=['AtomManipulator_useful'],
    dimensions=[6],
    mnemonics="OBJect B: AXis",
    characteristics=['[[INPUT_ONLY]]', '[[LENGTH]]'],
    commentdefault="""[[objbax]] must be provided if ([[nobj]] == 1 and one component of [[objaro]] != 0). Moreover,
[[objaax]] AND [[objbax]] must be provided if ( [[nobj]] == 2 and one component of [[objbro]] != 0 ).""",
    added_in_version="before_v9",
    text=r"""
Gives, for each object, the cartesian coordinates of two points (first point:
[[objbax]](1:3) second point: [[objbax]](4:6). By default, given in Bohr
atomic units (1 Bohr=0.5291772108 Angstroms), although Angstrom can be
specified, if preferred, since these variables have the [[LENGTH]] characteristics.
The two points define an axis of rotation of the corresponding object.
Note that the rotation of the object is done BEFORE the object is translated.
The sign of the rotation angle is positive if the object is to be rotated
clockwise when looking to it along the axis, from point 1 (coordinates 1:3)
toward point 2 (coordinates 4:6).
""",
),

Variable(
    abivarname="objbn",
    varset="geo",
    vartype="integer",
    topics=['AtomManipulator_useful'],
    dimensions="scalar",
    mnemonics="OBJect B: Number of atoms",
    characteristics=['[[INPUT_ONLY]]'],
    commentdefault=" [[objan]] and [[objbn]] MUST be provided if [[nobj]] == 2.",
    added_in_version="before_v9",
    text=r"""
Gives the number of atoms in either object b. The list of atoms is given by
the variables [[objbat]].
""",
),

Variable(
    abivarname="objbrf",
    varset="geo",
    vartype="integer",
    topics=['AtomManipulator_useful'],
    dimensions=[3],
    defaultval=[1, 1, 1],
    mnemonics="OBJect B: Repetition Factors",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Gives three repetition factors of the objects a or b.
This gives the opportunity to generate a three-dimensional set of repeated
objects, although a simple one-dimensional repetition will be easily obtained
through the specification of
`nrep 1 1` where `nrep` is the 1D repetition factor.
The initial rotation and translation of the object, as well as the increment
of rotation or translation from one object to the next are specified by the
variables [[objbro]] and [[objbtr]], for object b.
Note that the atom manipulator will generate the full set of atoms from the
primitive set of atoms using the following order: it will process each atom
in the primitive list one by one, determine whether it belongs to either
object a or object b, and then repeat it taking into account the proper
rotation and translation, with the fastest varying repetition factor being the
first, then the second, then the third.
In the final list of atoms, one will first find the atoms generated from atom
1 in the primitive list, then those generated from atom 2 in the primitive
list, and so on.
If the atom manipulator is only used to rotate or translate an object, without
repeating it, simply use `1 1 1`, which is also the Default value.
""",
),

Variable(
    abivarname="objbro",
    varset="geo",
    vartype="real",
    topics=['AtomManipulator_useful'],
    dimensions=[4],
    defaultval=MultipleValue(number=4, value=0.0),
    mnemonics="OBJect B: ROtations",
    characteristics=['[[INPUT_ONLY]]'],
    commentdefault="(no rotation)",
    added_in_version="before_v9",
    text=r"""
Give, for each object, the angles of rotation in degrees to be applied to the
corresponding object.
The rotation is applied before the translation, and the axis is defined by the
variable [[objbax]]. See the latter variable for the
definition of the sign of the rotation.
The first component [[objbro]](1) gives the angle of
rotation to be applied to the first instance of the object. The second, third
or fourth component (resp.) gives the increment of rotation angle from one
instance to the next instance, defined by the first, second or third
repetition factor (resp.). This allows one to generate 3D arrays of molecules
with different rotation angles.
""",
),

Variable(
    abivarname="objbtr",
    varset="geo",
    vartype="real",
    topics=['AtomManipulator_useful'],
    dimensions=[12],
    defaultval=MultipleValue(number=12, value=0.0),
    mnemonics="OBJect B: TRanslations",
    characteristics=['[[INPUT_ONLY]]', '[[LENGTH]]'],
    commentdefault="(no translation)",
    added_in_version="before_v9",
    text=r"""
Give, for each object, the vectors of translations, in cartesian coordinates,
to be applied to the corresponding object. By default, given in Bohr atomic
units (1 Bohr=0.5291772108 Angstroms), although Angstrom can be specified, if
preferred, since these variables have the [[LENGTH]] characteristics.
The translation is applied after the rotation.
The first vector [[objatr]](3,1) and [[objbtr]](3,1) gives the translation to
be applied to the first instance of the object. The second, third or fourth
component (resp.) gives the increment of translation from one instance to the
next instance, defined by the first, second or third repetition factor (resp.).
This allows one to generate 3D arrays of molecules.
In general, when the objects are repeated, a translation vector must be given,
since otherwise, the repeated objects pack in the same region of space. As an
exception, one can have a set of molecules regularly spaced on a circle, in
which case, only rotations are needed.
""",
),

Variable(
    abivarname="occ",
    varset="gstate",
    vartype="real",
    topics=['BandOcc_basic'],
    dimensions=['[[nband]]', "[[mband]]", "[[nsppol]]"],
    defaultval=MultipleValue(number=None, value=0),
    mnemonics="OCCupation numbers",
    characteristics=['[[EVOLVING]]'],
    added_in_version="before_v9",
    text=r"""
Gives occupation numbers for all bands in the problem. Needed if [[occopt]] == 0
or [[occopt]] == 2. Ignored otherwise. Also ignored when [[iscf]] = -2.
Typical band occupancy is either 2 or 0, but can be 1 for half-occupied band
or other choices in special circumstances.

If [[occopt]] is not 2, then the occupancies must be the same for each k point.

If [[occopt]] = 2, then the band occupancies must be provided explicitly for
each band, EACH k POINT, and EACH SPIN-POLARIZATION, in an array which runs
over all bands, k points, and spin-polarizations.
The order of entries in the array would correspond to all bands at the first k
point (spin up), then all bands at the second k point (spin up), etc, then all
k-points spin down.
The total number of array elements which must be provided is
( [[nband]](1)+[[nband]](2)+...+ [[nband]]([[nkpt]]) ) * [[nsppol]].
The occupation numbers evolve only for metallic occupations, that is, [[occopt]] ≥ 3.
""",
),

Variable(
    abivarname="occopt",
    varset="basic",
    vartype="integer",
    topics=['BandOcc_basic', 'STM_compulsory'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="OCCupation OPTion",
    added_in_version="before_v9",
    text=r"""
Controls how input parameters [[nband]], [[occ]], and [[wtk]] are handled.

  * [[occopt]] = 0:
All k points and spins have the same number of bands. All k points have the same occupancies of bands for a given spin
(but these occupancies may differ for spin up and spin down - typical for ferromagnetic insulators).
[[nband]] is given as a single number, and [[occ]]([[nband]] * [[nsppol]]) is an array of
[[nband]] * [[nsppol]] elements, read in by the code.
The k point weights in array [[wtk]]([[nkpt]]) are automatically normalized by
the code to add to 1. They cannot differ for differing spins.

  * [[occopt]] = 1:
Same as [[occopt]] = 0, except that the array [[occ]] is automatically generated
by the code, to give a semiconductor.
An error occurs when filling cannot be done with occupation numbers equal to 2
or 0 in each k-point (non-spin-polarized case), or with occupation numbers
equal to 1 or 0 in each k-point (spin-polarized case). If [[nsppol]] = 2 and
[[occopt]] == 1 is used, the user has to impose the magnetization, using
[[spinmagntarget]], except for the case of a single isolated Hydrogen atom.

  * [[occopt]] = 2:
k points may optionally have different numbers of bands and different
occupancies. [[nband]]([[nkpt]] * [[nsppol]]) is given explicitly as an array of
[[nkpt]] * [[nsppol]] elements. [[occ]]() is given explicitly for all bands at
each k point, and eventually for each spin -- the total number of elements is
the sum of [[nband]](ikpt) over all k points and spins. The k point weights
[[wtk]] ([[nkpt]]) are NOT automatically normalized under this option.

  * [[occopt]] = 3, 4, 5, 6 and 7:
Metallic occupation of levels, using different occupation schemes (see below).
The corresponding thermal broadening, or cold smearing, is defined by the
input variable [[tsmear]] (see below: the variable xx is the energy in Ha,
divided by [[tsmear]])
Like for [[occopt]] = 1, the variable [[occ]] is not read.
All k points have the same number of bands, [[nband]] is given as a single
number, read by the code.
The k point weights in array [[wtk]]([[nkpt]]) are automatically normalized by
the code to add to 1.

    * [[occopt]] = 3:
Fermi-Dirac smearing (finite-temperature metal) Smeared delta function:
0.25/(cosh(xx/2.0)**2)

    * [[occopt]] = 4:
"Cold smearing" of N. Marzari (see his thesis work), with a=-.5634
(minimization of the bump)
Smeared delta function:
exp(-xx  2  )/sqrt(pi) * (1.5+xx*(-a*1.5+xx*(-1.0+a*xx)))

    * [[occopt]] = 5:
"Cold smearing" of N. Marzari (see his thesis work), with a=-.8165 (monotonic
function in the tail)
Same smeared delta function as [[occopt]] = 4, with different a.

    * [[occopt]] = 6:
Smearing of Methfessel and Paxton [[cite:Methfessel1989]] with Hermite polynomial
of degree 2, corresponding to "Cold smearing" of N. Marzari with a=0 (so, same
smeared delta function as [[occopt]] = 4, with different a).

    * [[occopt]] = 7:
Gaussian smearing, corresponding to the 0 order Hermite polynomial of
Methfessel and Paxton.
Smeared delta function: 1.0*exp(-xx**2)/sqrt(pi)

    * [[occopt]] = 8:
Uniform smearing (the delta function is replaced by a constant function of
value one over ]-1/2,1/2[ (with one-half value at the boundaries). Used for
testing purposes only.

!!! note

    One can use metallic occupation of levels in the case of a molecule,
    in order to avoid any problem with degenerate levels. However, it is advised
    NOT to use [[occopt]] = 6 (and to a lesser extent [[occopt]] = 4 and 5), since the
    associated number of electron versus the Fermi energy is NOT guaranteed to be
    a monotonic function. For true metals, AND a sufficiently dense sampling of
    the Brillouin zone, this should not happen, but be cautious ! As an indication
    of this problem, a small variation of input parameters might lead to a jump of
    total energy, because there might be two or even three possible values of the
    Fermi energy, and the bisection algorithm finds one or the other.
""",
),

Variable(
    abivarname="omegasimax",
    varset="gw",
    vartype="real",
    topics=['SelfEnergy_basic'],
    dimensions="scalar",
    defaultval=ValueWithUnit(units='eV', value=50),
    mnemonics="OMEGA to evaluate Sigma along the Imaginary axis D: MAXimal value",
    characteristics=['[[ENERGY]]'],
    requires="[[optdriver]] == 4 and [[gwcalctyp]] == 1",
    added_in_version="before_v9",
    text=r"""
[[omegasimax]] defines the maximum frequency along the imaginary the axis. In
conjunction with [[nomegasi]], [[omegasimax]] uniquely defines the linear mesh
employed to sample the self-energy along the imaginary axis.
""",
),

Variable(
    abivarname="omegasrdmax",
    varset="gw",
    vartype="real",
    topics=['SelfEnergy_basic'],
    dimensions="scalar",
    defaultval=ValueWithUnit(units='eV', value=1.0),
    mnemonics="OMEGA to evaluate the Sigma Real axis Derivative: MAXimal value",
    characteristics=['[[ENERGY]]'],
    requires="[[optdriver]] == 4",
    added_in_version="before_v9",
    text=r"""
The maximum distance from the KS energy where to evaluate Sigma. Sigma is
evaluated at [ KS_energy - [[omegasrdmax]], KS_energy + [[omegasrdmax]] ]
sampled [[nomegasrd]] times.
""",
),

Variable(
    abivarname="optcell",
    varset="rlx",
    vartype="integer",
    topics=['PIMD_basic', 'GeoOpt_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="OPTimize the CELL shape and dimensions",
    added_in_version="before_v9",
    text=r"""
Allows one to optimize the unit cell shape and dimensions, when [[ionmov]] >= 2 or
3. The configuration for which the stress almost vanishes is iteratively
determined, by using the same algorithms as for the nuclei positions. Will
eventually modify [[acell]] and/or [[rprim]]. The ionic positions are ALWAYS
updated, according to the forces. A target stress tensor might be defined, see [[strtarget]].

  * **optcell** = 0: modify nuclear positions, since [[ionmov]] = 2 or 3, but no cell shape and dimension optimisation.
  * **optcell** = 1: optimisation of volume only (do not modify [[rprim]], and allow an homogeneous dilatation of the three components of [[acell]])
  * **optcell** = 2: full optimization of cell geometry (modify [[acell]] and [[rprim]] \- normalize the vectors of [[rprim]] to generate the [[acell]]). This is the usual mode for cell shape and volume optimization. It takes into account the symmetry of the system, so that only the effectively relevant degrees of freedom are optimized.
  * **optcell** = 3: constant-volume optimization of cell geometry (modify [[acell]] and [[rprim]] under constraint \- normalize the vectors of [[rprim]] to generate the [[acell]])
  * **optcell** = 4, 5 or 6: optimize [[acell]](1), [[acell]](2), or [[acell]](3), respectively (only works if the two other vectors are orthogonal to the optimized one, the latter being along its cartesian axis).
  * **optcell** = 7, 8 or 9: optimize the cell geometry while keeping the first, second or third vector unchanged (only works if the two other vectors are orthogonal to the one left unchanged, the latter being along its cartesian axis).

A few details require attention when performing unit cell optimisation:

  * one has to get rid of the discontinuities due to discrete changes of plane wave number with cell size, by using a suitable value of [[ecutsm]];
  * one has to allow for the possibility of a larger sphere of plane waves, by using [[dilatmx]];
  * one might have to adjust the scale of stresses to the scale of forces, by using [[strfact]].
  * if all the reduced coordinates of atoms are fixed by symmetry, one cannot use [[toldff]] to stop the SCF cycle. (Suggestion: use [[toldfe]] with a small value, like 1.0d-10)

It is STRONGLY suggested first to optimize the ionic positions without cell
shape and size optimization (**optcell** = 0), then start the cell shape and
size optimization from the cell with relaxed ionic positions. Presently
(v3.1), one cannot restart ([[restartxf]]) a calculation with a non-zero
**optcell** value from the (x,f) history of another run with a different non-
zero **optcell** value. There are still a few problems at that level.
""",
),

Variable(
    abivarname="optdriver",
    varset="gstate",
    vartype="integer",
    topics=['nonlinear_compulsory',
 'GWls_compulsory',
 'ElPhonInt_compulsory',
 'GW_compulsory',
 'BSE_compulsory',
 'DFPT_compulsory',
 'Susceptibility_compulsory',
 'SelfEnergy_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="OPTions for the DRIVER",
    added_in_version="before_v9",
    text=r"""
For each dataset, choose the task to be done, at the level of the "driver" routine.

The choice is among:

  * 0 --> ground-state calculation (GS), routine *gstate*
  * 1 --> response-function calculation (RF), routine *respfn*
  * 2 --> susceptibility calculation (SUS), routine *suscep*
  * 3 --> susceptibility and dielectric matrix calculation (SCR), routine *screening*
  * 4 --> self-energy calculation (SIG), routine *sigma*.
  * 5 --> non-linear response functions (NONLINEAR), using the 2n+1 theorem, routine *nonlinear*.
  * 7 --> electron-phonon coupling (EPH)
  * 8 --> Post-processing of WFK file, routine *wfk_analyze*. See also [[wfk_task]] input variable.
  * 66 --> GW using Lanczos-Sternheimer, see input variables whose name start with `gwls_*`.
  * 99 --> Bethe-Salpeter calculation (BSE), routine *bethe_salpeter*

If one of [[rfphon]], [[rfddk]], [[rfelfd]], or [[rfstrs]] is non-zero, while
**optdriver** is not defined in the input file, ABINIT will set **optdriver**
to 1 automatically. These input variables ([[rfphon]], [[rfddk]], [[rfelfd]],
and [[rfstrs]]) must be zero if **optdriver** is not set to 1.
""",
),

Variable(
    abivarname="optforces",
    varset="dev",
    vartype="integer",
    topics=['ForcesStresses_basic'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[toldff]] or [[tolrff]] != 0': 1, 'defaultval': 2}),
    mnemonics="OPTions for the calculation of FORCES",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Allows to choose options for the calculation of forces.

  * **optforces** = 0: the forces are set to zero, and many steps of the computation of forces are skipped
  * **optforces** = 1: calculation of forces at each SCF iteration, allowing to use forces as criterion to stop the SCF cycles
  * **optforces** = 2: calculation of forces at the end of the SCF iterations (like the stresses)
""",
),

Variable(
    abivarname="optnlxccc",
    varset="dev",
    vartype="integer",
    topics=['xc_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="OPTion for the calculation of Non-Linear eXchange-Correlation Core Correction",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Allows to choose options for the calculation of non-linear XC correction. At
present, only relevant for the FHI type of pseudopotentials, with pspcod=6.

  * **optnlxccc** = 1: uses the old `psp6cc.f` routine, with inconsistent treatment of real-space derivatives of the core
  function (computed in this routine, while splined in the other parts of the code)
  * **optnlxccc** = 2: consistent calculation derivatives, in the `psp6cc_dhr.f` routine from DHamann.
""",
),

Variable(
    abivarname="optstress",
    varset="gstate",
    vartype="integer",
    topics=['ForcesStresses_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="OPTion for the computation of STRESS",
    added_in_version="before_v9",
    text=r"""
If set to 1, the computation of stresses is done, in the SCF case (under the
conditions [[iscf]] > 0, [[prtstm]] == 0, [[positron]] == 0, and either
[[nstep]] >0, or [[usepaw]] == 0 or [[irdwfk]] == 1).
Otherwise, to save CPU time, if no optimization of the cell is required, one
can skip the computation of stresses. The CPU time saving might be interesting
for some PAW calculations.
""",
),

Variable(
    abivarname="orbmag",
    varset="gstate",
    vartype="integer",
    topics=['MagField_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="ORBital MAGnetization",
    characteristics=['[[DEVELOP]]'],
    requires="""[[usepaw]] == 1;
[[usexcnhat]] == 0;
[[nspinor]] == 1;
[[paral_atom]] == 0;
[[paral_kgb]] == 0;
[[kptopt]] == 3 """,
    added_in_version="before_v9",
    text=r"""
Compute quantities related to orbital magnetization. The
    implementation assumes an insulator, so no empty or partially
    filled bands, and currently restricted to [[nspinor]] 1. Such
    insulators have orbital magnetization zero, except in the presence
    of nonzero nuclear dipole moments, see [[nucdipmom]].  [[orbmag]]
    is parallelized over k points only. The implementation follows the
    theory outlined in [[cite:Gonze2011a]] extended to the PAW case;
    see also [[cite:Ceresoli2006]]. The computed results are returned in the
    standard output file, search for "Orbital magnetization" and "Chern number".

* [[orbmag]] = 1: Compute Chern number (really, the integral of the Berry curvature
over the Brillouin zone) [[cite:Ceresoli2006]]. This computation is
faster than the full [[orbmag]] calculation, and a nonzero value indicates a circulating
electronic current.
* [[orbmag]] = 2: Compute electronic orbital magnetization.
* [[orbmag]] = 3: Compute both Chern number and electronic orbital magnetization.

The above settings use an implementation based on a discretization of the wavefunction
derivatives, as in [[cite:Ceresoli2006]]. Using [[orbmag]] -1, -2, -3 delivers the
same computations as the corresponding 1, 2, 3 values, but based on an implementation
using a discretization of the density operator itself. Both methods should converge to
the same values but in our experience the wavefunction-based method converges faster.
""",
),

Variable(
    abivarname="ortalg",
    varset="dev",
    vartype="integer",
    topics=['TuningSpeed_expert'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[wfoptalg]] >= 10 ': -2, 'defaultval': 2}),
    mnemonics="ORThogonalisation ALGorithm",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Allows to choose the algorithm for orthogonalisation.
Positive or zero values make two projections per line minimisation, one before
the preconditioning, one after. This is the clean application of the band-by-
band CG gradient for finding eigenfunctions.
Negative values make only one projection per line minimisation.
The orthogonalisation step is twice faster, but the convergence is less good.
This actually calls to a better understanding of this effect.
**ortalg** = 0, 1 or -1 is the conventional coding.
**ortalg** = 2 or -2 try to make better use of existing registers on the
particular machine one is running.
More demanding use of registers is provided by **ortalg** = 3 or -3, and so on.
The maximal value is presently 4 and -4.
Tests have shown that **ortalg** = 2 or -2 is suitable for use on the available platforms.
""",
),

Variable(
    abivarname="papiopt",
    varset="dev",
    vartype="integer",
    topics=['Control_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PAPI OPTion",
    added_in_version="before_v9",
    text=r"""
[PAPI](http://icl.cs.utk.edu/papi/index.html) aims to provide the tool
designer and application engineer with a consistent interface and methodology
for use of the performance counter hardware found in most major
microprocessors. PAPI enables software engineers to see, in near real time,
the relation between software performance and processor events.
This option can be used only when ABINIT has been compiled with the `
--enable-papi ` configure option.
If **papiopt** = 1, then PAPI counters are used instead of the usual time()
routine. All the timing output of ABINIT is then done with PAPI values. The
measurements are more accurate and give also access to the flops of the calculation.
""",
),

Variable(
    abivarname="paral_atom",
    varset="paral",
    vartype="integer",
    topics=['parallelism_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="activate PARALelization over (paw) ATOMic sites",
    added_in_version="before_v9",
    text=r"""
Relevant only for PAW calculations.
This keyword controls the parallel distribution of memory over atomic sites.
Calculations are also distributed using the "kpt-band" communicator.
Compatible with ground-state calculations and response function calculations
""",
),

Variable(
    abivarname="paral_kgb",
    varset="paral",
    vartype="integer",
    topics=['parallelism_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="activate PARALelization over K-point, G-vectors and Bands",
    added_in_version="before_v9",
    text=r"""
**If paral_kgb is not explicitely put in the input file**, ABINIT
automatically detects if the job has been sent in sequential or in parallel.
In this last case, it detects the number of processors on which the job has
been sent and calculates values of [[npkpt]], [[npfft]], [[npband]],
[[bandpp]], [[npimage]] and [[npspinor]] that are compatible with the number
of processors. It then set **paral_kgb** to 0 or 1 (see hereunder) and launches the job.

**If paral_kgb = 0**, the parallelization over k-points only is activated. In
this case, [[npkpt]], [[npspinor]], [[npfft]] and [[npband]] are ignored.
Require compilation option --enable-mpi="yes".

**If paral_kgb = 1**, the parallelization over bands, FFTs, and k-point/spin-
components is activated (see [[npkpt]], [[npfft]] [[npband]] and eventually
[[npspinor]]). With this parallelization, the work load is split over four
levels of parallelization (three level of parallelisation (kpt-band-fft )+
spin) The different communications almost occur along one dimension only.
Require compilation option --enable-mpi="yes".

HOWTO fix the number of processors along one level of parallelisation:
At first, try to parallelise over the k point and spin (see
[[npkpt]],[[npspinor]]). Otherwise, for unpolarized calculation at the gamma
point, parallelise over the two other levels: the band and FFT ones. For nproc $\leq$ 50, the best speed-up is achieved for [[npband]] = nproc and [[npfft]] = 1
(which is not yet the default). For nproc $\geq$ 50, the best speed-up is achieved
for [[npband]] $\geq$ 4 $\times$ [[npfft]].

For additional information, download F. Bottin presentation at the
[ABINIT workshop 2007](https://www.abinit.org/sites/default/files/oldsites/workshop_07/program.html)

Suggested acknowledgments:
[[cite:Bottin2008]], also available on arXiv, http://arxiv.org/abs/0707.3405.

If the total number of processors used is compatible with the four levels of
parallelization, the values for [[npkpt]], [[npspinor]], [[npfft]], [[npband]]
and [[bandpp]] will be filled automatically, although the repartition may not
be optimal. To optimize the repartition use:

**If paral_kgb = 1** and **max_ncpus = n $\ne$ 0** ABINIT will test automatically
if all the processor numbers between 2 and n are convenient for a parallel
calculation and print the possible values in the log file. A weight is
attributed to each possible processors repartition. It is adviced to select a
processor repartition for which the weight is high (as closed to the number of
processors as possible). The code will then stop after the printing. This test
can be done as well with a sequential as with a parallel version of the code.
The user can then choose the adequate number of processor on which he can run
his job. He must put again paral_kgb = 1 in the input file and put the
corresponding values for [[npkpt]], [[npfft]], [[npband]],[[bandpp]] and
eventually [[npspinor]] in the input file.
""",
),

Variable(
    abivarname="paral_rf",
    varset="paral",
    vartype="integer",
    topics=['parallelism_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Activate PARALlelization over Response Function perturbations",
    added_in_version="before_v9",
    text=r"""
This parameter activates the parallelization over perturbations which can be
used during RF-Calculation. It is possible to use this type of parallelization
in combination to the parallelization over k-points.

Currently total energies calculated by groups, where the master process is not
in, are saved in.status_LOGxxxx files.

If **paral_rf** is set to -1, the code reports the list of irreducible
perturbations for the specified q-point in the log file (YAML format) and then stops.

**paral_rf** can be specified separately for each dataset.
""",
),

Variable(
    abivarname="pawcpxocc",
    varset="paw",
    vartype="integer",
    topics=['PAW_expert'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[optdriver]] == 0 and [[ionmov]] < 6 and [[pawspnorb]] == 1 and [[iscf]] >= 10 and ([[kptopt]] !=1 or [[kptopt]]!=2) and [[usepaw]] == 1': 2,
 'defaultval': 1}),
    mnemonics="PAW - use ComPleX rhoij OCCupancies",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
The only possible values for [[pawcpxocc]] are 1 or 2.
When [[pawcpxocc]] == 1, "direct" decomposition of total energy cannot be
printed out.
When [[pawcpxocc]] == 2, PAW augmentation occupancies are treated as
COMPLEX; else they are considered as REAL.
This is needed when time-reversal symmetry is broken (typically when spin-
orbit coupling is activated).

Note for ground-state calculations ([[optdriver]] == 0):
The imaginary part of PAW augmentation occupancies is only used for the
computation of the total energy by "direct scheme"; this is only necessary
when SCF mixing on potential is chosen ([[iscf]] <10).
When SCF mixing on density is chosen ([[iscf]] >= 10), the "direct"
decomposition of energy is only printed out without being used. It is thus
possible to use [[pawcpxocc]] = 1 in the latter case.
In order to save CPU time, when molecular dynamics is selected ([[ionmov]] >= 6)
and SCF mixing done on density ([[iscf]] >= 10), [[pawcpxocc]] = 2 is (by default) set to **1**.
""",
),

Variable(
    abivarname="pawcross",
    varset="paw",
    vartype="integer",
    topics=['DFPT_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PAW - add CROSS term in oscillator strengths",
    requires="([[optdriver]] == 3 or [[optdriver]] == 4) and [[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
When **pawcross=1**, the overlap between the plane-wave part of one band and
the on-site part of an other is taken into account in the computation of the
oscillator strengths. Hence, the completeness of the on-site basis is no longer assumed.
""",
),

Variable(
    abivarname="pawecutdg",
    varset="paw",
    vartype="real",
    topics=['Planewaves_compulsory', 'PAW_compulsory'],
    dimensions="scalar",
    defaultval=-1,
    mnemonics="PAW - Energy CUToff for the Double Grid",
    characteristics=['[[ENERGY]]'],
    commentdefault="pawecutdg MUST be specified for PAW calculations.",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
Define the energy cut-off for the fine FFT grid (the "double grid", that
allows one to transfer data from the normal, coarse, FFT grid to the spherical
grid around each atom).
[[pawecutdg]] must be larger or equal to [[ecut]]. If it is equal to it, then
no fine grid is used. The results are not very accurate, but the computations
proceed quite fast.
For typical PAW computations, where [[ecut]] is on the order of 15 Ha,
[[pawecutdg]] must be tested according to what you want to do. For
calculations that do not require a high accuracy (molecular dynamics for
instance) a value of 20 Ha is enough. For calculations that require a high
accuracy (response functions for instance) it should be on the order of 30 Ha.
Choosing a larger value should not increase the accuracy, but does not slow
down the computation either, only the memory. The choice made for this
variable DOES have a bearing on the numerical accuracy of the results, and, as
such, should be the object of a convergence study. The convergence test might
be made on the total energy or derived quantities, like forces, but also on
the two values of the "Compensation charge inside spheres", a quantity written
in the log file.
""",
),

Variable(
    abivarname="pawfatbnd",
    varset="paw",
    vartype="integer",
    topics=['PAW_useful', 'ElecBandStructure_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PAW: print band structure in the FAT-BaND representation",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
For Ground-State calculations and non self-consistent calculations only.
This option can be used to plot band structure. For each atom (specified by
[[natsph]] and [[iatsph]]), each angular momentum, and each spin polarisation,
the band structure is written in files (such as e.g.
FATBANDS_at0001_Ni_is2_l2_m-1). Each file contains the eigenvalue, and the
contribution of angular momentum L, and projection of angular momentum M, (for
the corresponding wavefunction) to the PAW density inside the PAW sphere as a
function of the index of the k-point. The output can be readily plotted with
the software [xmgrace](http://plasma-gate.weizmann.ac.il/Grace/) (e.g
xmgrace FATBANDS_at0001_Ni_is2_l2_m-1). Relevant values are:

  * 0: desactivated.
  * 1: The fatbands are only resolved in L.
  * 2: The fatbands are resolved in L and M.
""",
),

Variable(
    abivarname="pawlcutd",
    varset="paw",
    vartype="integer",
    topics=['PAW_expert'],
    dimensions="scalar",
    defaultval=10,
    mnemonics="PAW - L angular momentum used to CUT the development in moments of the Densities",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
The expansion of the densities in angular momenta is performed up to
l=[[pawlcutd]].
Note that, for a given system, the maximum value of [[pawlcutd]] is
**2*l_max**, where l_max is the maximum l of the PAW partial waves basis.

The choice made for this variable DOES have a bearing on the numerical
accuracy of the results, and, as such, should be the object of a convergence
study. The convergence test might be made on the total energy or derived
quantities, like forces, but also on the two values of the "Compensation
charge inside spheres", a quantity written in the log file.
""",
),

Variable(
    abivarname="pawlmix",
    varset="paw",
    vartype="integer",
    topics=['PAW_expert'],
    dimensions="scalar",
    defaultval=10,
    mnemonics="PAW - maximum L used in the spherical part MIXing",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
The choice made for this variable determine how the spherical part of the
density is mixed during electronic iterations.

Only parts of rhoij quantities associated with l angular momenta up to
l=pawlmix are mixed. Other parts of augmentation occupancies are not included
in the mixing process.
This option is useful to save CPU time but DOES have a bearing on the
numerical accuracy of the results.
""",
),

Variable(
    abivarname="pawmixdg",
    varset="paw",
    vartype="integer",
    topics=['PAW_expert'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[npfft]] == 1': 0, 'defaultval': 1}),
    mnemonics="PAW - MIXing is done (or not) on the (fine) Double Grid",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
The choice made for this variable determines the grid on which the density (or
potential) is mixed during the SCF cycle.

If **pawmixdg=1** the density/potential is mixed in REAL space using the
fine FFT grid (defined by [[pawecutdg]] or [[ngfftdg]]).

If **pawmixdg=0** the density/potential is mixed in RECIPROCAL space using
the coarse FFT grid (defined by [[ecut]] or [[ngfft]]). Only components of the
coarse grid are mixed using the scheme defined by [[iscf]]; other components
are only precondionned by [[diemix]] and simply mixed.
This option is useful to save memory and does not affect numerical accuracy of
converged results. If **pawmixdg=1**, density and corresponding residual are
stored for previous iterations and are REAL arrays of size [[nfftdg]]. If
**pawmixdg=0**, density and corresponding residual are stored for previous
iterations and are COMPLEX arrays of size [[nfft]]. The memory saving is
particularly efficient when using the Pulay mixing ([[iscf]] = 7 or 17).

In **wavelet** calculations [[usewvl]] = 1:

    - pawmixdg is set to 1 by default.
    - A value of 0 is not allowed.
    - Density/potential is mixed in REAL space (Here only one grid is used).
""",
),

Variable(
    abivarname="pawnhatxc",
    varset="paw",
    vartype="integer",
    topics=['PAW_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="PAW - Flag for exact computation of gradients of NHAT density in eXchange-Correlation.",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
Relevant only when a GGA exchange-correlation functional is used.
When this flag is activated, the gradients of compensation charge density
(n_hat) are exactly computed (i.e. analytically); when it is deactivated, they
are computed with a numerical scheme in reciprocal space (which can produce
inaccurate results if the compensation charge density is highly localized).
As analytical treatment of compensation charge density gradients is CPU time
demanding, it is possible to bypass it with [[pawnhatxc]] = 0; but the numerical
accuracy can be affected by this choice. It is recommended to test the
validity of this approximation before use.
""",
),

Variable(
    abivarname="pawnphi",
    varset="paw",
    vartype="integer",
    topics=['PAW_expert'],
    dimensions="scalar",
    defaultval=13,
    mnemonics="PAW - Number of PHI angles used to discretize the sphere around each atom.",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
Number of phi angles (longitude) used to discretize the data on the atomic
spheres. This discretization is completely defined by [[pawnphi]] and [[pawntheta]].
""",
),

Variable(
    abivarname="pawntheta",
    varset="paw",
    vartype="integer",
    topics=['PAW_expert'],
    dimensions="scalar",
    defaultval=12,
    mnemonics="PAW - Number of THETA angles used to discretize the sphere around each atom.",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
Number of theta angles (latitude) used to discretize the data on the atomic
spheres. This discretization is completely defined by [[pawntheta]] and [[pawnphi]].
""",
),

Variable(
    abivarname="pawnzlm",
    varset="paw",
    vartype="integer",
    topics=['PAW_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="PAW - only compute Non-Zero LM-moments of the contributions to the density from the spheres",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
Concerns the computation of the contributions to the density from the spheres
(named rho_1 - rho_tild_1).
If set to 0, all lm-moments of the sphere contributions to the density are
computed at each electronic iteration.
If set to 1, only non-zero lm-moments of the sphere contributions to the
density are computed at each electronic iteration (they are all computed at
the first iteration then only those found to be non-zero will be computed;
thus the first iteration is more cpu intensive)
""",
),

Variable(
    abivarname="pawoptmix",
    varset="paw",
    vartype="integer",
    topics=['PAW_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PAW - OPTion for the MIXing of the spherical part",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
In the case of PAW computations, during the self-consistent cycle, ABINIT
mixes the density $\rho(r)= \tilde{\rho}(r) +\hat{\rho}(r)$ and the occupancy matrix $\rho_{ij}$. ($\tilde{\rho}(r)$ is
the pseudo density, $\hat{\rho}(r)$ is the compensation charge density). It can be
redundant as $\rho_{ij}$ is contained in $\hat{\rho}(r)$.

  * If **pawoptmix** =0:
ABINIT mixes $\rho(r)$ and $\rho_{ij}$ but the residual used to control the mixing
algorithm is only based on $\rho(r)$.

  * If **pawoptmix** =1:
ABINIT mixes $\rho(r)$ and $\rho_{ij}$ and the residual used to control the mixing
algorithm is based on $\rho(r)$ and $\rho_{ij}$.

This has only an influence on the efficiency of the mixing algorithm.
In case of mixing problems, the first suggestion is to increase the size of the
history (see [[npulayit]]). Then it is also possible to play with the
parameters of the Kerker mixing: [[diemix]], [[diemac]], etc...
""",
),

Variable(
    abivarname="pawoptosc",
    varset="paw",
    vartype="integer",
    topics=['Susceptibility_expert', 'BSE_expert', 'SelfEnergy_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PAW - OPTion for the computation of the OSCillator matrix elements",
    added_in_version="before_v9",
    text=r"""
Only relevant for GW or Bethe-Salpeter calculations with PAW.
This variable defines the approach used for the evaluation of the oscillator
matrix elements within the PAW formalism. Possible values are 0,1,2.
If [[pawoptosc]] = 0 the code uses its internal default value (2 for SCREENING
calculations, 1 for SIGMA calculations, 2 for Bethe-Salpeter
If [[pawoptosc]] = 1 the matrix elements are computed with the expression given
by [[cite:Arnaud2000]]. The equation is exact provided that the
set of PAW partial waves is complete.
If [[pawoptosc]] = 2 the matrix elements are computed with the approximated
expression proposed by [[cite:Shishkin2006]].
""",
),

Variable(
    abivarname="pawovlp",
    varset="paw",
    vartype="real",
    topics=['PAW_useful'],
    dimensions="scalar",
    defaultval=5.0,
    mnemonics="PAW - spheres OVerLaP allowed (in percentage)",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
When PAW is activated, a localized atomic basis is added to describe wave
functions. Spheres around atoms are defined and they are IN PRINCIPLE not
allowed to overlap. However, a small overlap can be allowed without
compromising the accuracy of results. Be aware that too high overlaps can lead
to unphysical results.
With the **pawovlp** variable, the user can control the (voluminal) overlap
percentage allowed without stopping the execution.
**pawovlp** is the value (in percentage: 0...100%) obtained by dividing the
volume of the overlap of two spheres by the volume of the smallest sphere.
The following values are permitted for **pawovlp**:

- **pawovlp** < 0 --> overlap is always allowed
- **pawovlp** = 0 --> no overlap is allowed
- **pawovlp** > 0 and < 100 --> overlap is allowed only if it is less than **pawovlp** %
""",
),

Variable(
    abivarname="pawprt_b",
    varset="dev",
    vartype="integer",
    topics=['PAW_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PAW PRinT band",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Forces the output of the all-electron wavefunction for only a single band.
To be used in conjunction with: [[pawprtwf]] = 1 and [[pawprt_k]].
The indexing of the bands start with one for the lowest occupied band and goes up from there.
""",
),

Variable(
    abivarname="pawprt_k",
    varset="dev",
    vartype="integer",
    topics=['PAW_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PAW PRinT K-point",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Forces the output of the all-electron wavefunction for only a single k-point.
To be used in conjunction with: [[pawprtwf]] = 1 and [[pawprt_b]].
The indexing follows the order in output of the internal variable **kpt** in
the beginning of the run.
""",
),

Variable(
    abivarname="pawprtden",
    varset="paw",
    vartype="integer",
    topics=['PAW_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PAW: PRinT total physical electron DENsity",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
**Deprecated:** See the [[prtden]].
""",
),

Variable(
    abivarname="pawprtdos",
    varset="paw",
    vartype="integer",
    topics=['PAW_useful', 'ElecDOS_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PAW: PRinT partial DOS contributions",
    requires="[[usepaw]] == 1 and [[prtdos]] == 3",
    added_in_version="before_v9",
    text=r"""
This input variable controls the computation and/or printing of contributions
to the PAW partial DOS in _DOS file(s):

* Plane-waves contribution

  $+$ "on-site" all-electron contribution ($\phi$)

  $-$ "on-site" pseudo contribution ($\tilde{\phi}$).

If **pawprtdos=0:**

- The 3 contributions are computed; only the total partial DOS is output in _DOS file.

If **pawprtdos=1:**

- The 3 contributions are computed and output in _DOS file.
- In that case, integrated DOS is not output.

If **pawprtdos=2:**

- Only "on-site" all-electron contribution is computed and output in _DOS file.
- This a (very) good approximation of total DOS, provided that (1) the PAW
  local basis is complete, (2) the electronic charge is mostly contained in PAW spheres.
- In that case, the [[ratsph]] variable is automatically set to the PAW radius.
""",
),

Variable(
    abivarname="pawprtvol",
    varset="paw",
    vartype="integer",
    topics=['PAW_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PAW: PRinT VOLume",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
Control print volume and debugging output for PAW in log file or standard
output. If set to 0, the print volume is at its minimum.
**pawprtvol** can have values from -3 to 3:

- **pawprtvol** = -1 or 1: matrices $\rho_{ij}$ (atomic occupancies) and $D_{ij}$ (psp
  strength) are printed at each SCF cycle with details about their contributions.
- **pawprtvol** = -2 or 2: like -1 or 1 plus additional printing: moments of
  "on-site" densities, details about local exact exchange.
- **pawprtvol** = -3 or 3: like -2 or 2 plus additional printing: details about
  PAW+U, rotation matrices of spherical harmonics.

When **pawprtvol** >= 0, up to 12 components of $\rho_{ij}$ and $D_{ij}$ matrices for the
1st and last atom are printed.
When **pawprtvol** < 0, all components of $\rho_{ij}$ and $D_{ij}$ matrices for all atoms are printed.
""",
),

Variable(
    abivarname="pawprtwf",
    varset="paw",
    vartype="integer",
    topics=['PAW_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PAW: PRinT WaveFunctions",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
This input variable controls the output of the **full** PAW wave functions
including the on-site contributions inside each PAW sphere needed to
reconstruct the correct nodal shape in the augmentation region.
**pawprtwf = 1** causes the generation of a file _PAWAVES.nc containing the full wavefunctions
in **real space** on the fine FFT grid defined by [[pawecutdg]] or [[ngfftdg]].

Limitations: At present (v8.0), **pawprtwf = 1** is not compatible with [[nspinor]] =2 and parallel executions
Therefore the output of the _AE_WFK has to be done in sequential.
Moreover, in order to use this feature, one has to enable the support for netcdf at configure-time
as the _PAWAVES file is written using the NETCDF file format following the ETSF-IO
specifications for wavefunctions in real space.

If the code is run entirely in serial, additional output is made of various contributions to the all-electron
wavefunction. By default the full available set of bands and k-points are
output, but a single band and k-point index can be requested by using the
variables [[pawprt_b]] and [[pawprt_k]].
""",
),

Variable(
    abivarname="pawspnorb",
    varset="paw",
    vartype="integer",
    topics=['PAW_useful', 'spinpolarisation_useful'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[nspinor]] == 2': 1, 'defaultval': 0}),
    mnemonics="PAW - option for SPiN-ORBit coupling",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
When PAW is activated, the **spin-orbit coupling** can be added without the
use of specific PAW datasets (pseudopotentials).
If [[pawspnorb]] = 1, spin-orbit will be added.
If the wavefunction is spinorial (that is, if [[nspinor]] = 2), there is no
reason not to include the spin-orbit interaction, so that the default value of
[[pawspnorb]] becomes 1 when [[nspinor]] = 2.
Note that only the all-electron "on-site" contribution to the Hamiltonian is
taken into account; this is a very good approximation but requires the
following conditions to be fulfilled:

1- the  $\tilde{\phi}_{i}$  basis is complete enough

2- the electronic density is mainly contained in the PAW sphere


Also note that, when spin-orbit coupling is activated and there is some
magnetization [[nspden]] = 4, the time-reversal symmetry is broken.
The use of [[kptopt]] = 1 or [[kptopt]] = 2 is thus forbidden. It is advised to
use [[kptopt]] = 3 (no symmetry used to generate k-points) or [[kptopt]] = 4 (only
spatial symmetries used to generate k-points).
Be careful if you choose to use [[kptopt]] = 0 (k-points given by hand); Time-
reversal symmetry has to be avoided.
An artificial scaling of the spin-orbit can be introduced thanks to the
[[spnorbscl]] input variable.
""",
),

Variable(
    abivarname="pawstgylm",
    varset="paw",
    vartype="integer",
    topics=['PAW_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="PAW - option for the STorage of G_l(r).YLM(r)",
    requires="[[usepaw]] = 1",
    added_in_version="before_v9",
    text=r"""
When PAW is activated, the computation of compensation charge density (so
called "hat" density) requires the computation of $g_{l}(r).Y_{lm}(r)$ factors (and
cartesian derivatives) at each point of real space contained in PAW spheres.
The number of atoms, of (l,m) quantum numbers and the sharpness of the real
FFT grid can lead to a very big {$g_{l}.Y_{lm}$} datastructure. One can save memory
by putting [[pawstgylm]] = 0; but, in that case, $g_{l}(r).Y_{lm}(r)$ factors a re-
computed each time they are needed and CPU time increases.

Possible choices:

- [[pawstgylm]] = 0: $g_{l}(r).Y_{lm}(r)$ are not stored in memory and recomputed.
- [[pawstgylm]] = 1: $g_{l}(r).Y_{lm}(r)$ are stored in memory.

Note:
$g_{l}(r)$ are shape functions (analytically known)

$Y_{lm}(r)$ are real spherical harmonics
""",
),

Variable(
    abivarname="pawsushat",
    varset="paw",
    vartype="integer",
    topics=['PAW_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PAW - SUSceptibility, inclusion of HAT (compensation charge) contribution",
    requires="[[usepaw]] == 1 and [[optdriver]] == 0",
    added_in_version="before_v9",
    text=r"""
Ground-State calculation only.
When a sophisticated preconditioning scheme is selected for the SCF cycle of a
Ground-State calculation ([[iprcel]] > 0), the computation of the susceptibility
matrix is required several times during the cycle. This computation is
computer time consuming, especially -- within PAW -- because of the inclusion
of additional terms due to the compensation charge density. As only a crude
valuation of the susceptibilty matrix is needed (to evaluate a preconditioning
matrix), the compensation charge contribution can be neglected to save CPU
time (select [[pawsushat]] = 0). This approximation could be unfavourable in
some cases; in the latter, we advise to put [[pawsushat]] = 1.

Possible choices:

- [[pawsushat]] = 0: only plane-wave contribution to suscep. matrix is computed.
- [[pawsushat]] = 1: the whole suscep. matrix (PW + PAW on-site) is computed.
""",
),

Variable(
    abivarname="pawujat",
    varset="dev",
    vartype="integer",
    topics=['DFT+U_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="PAW+macro_UJ, ATom number",
    characteristics=['[[DEVELOP]]'],
    commentdefault=" i.e. the first atom treated with PAW+U.",
    added_in_version="before_v9",
    text=r"""
Determines the atom for which U (or J) should be determined. See also [[macro_uj]].
""",
),

Variable(
    abivarname="pawujrad",
    varset="dev",
    vartype="real",
    topics=['DFT+U_expert'],
    dimensions="scalar",
    defaultval=ValueWithUnit(units='a.u.', value=20),
    mnemonics="PAW+macro_UJ, sphere RADius",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
The sphere radius serves to extrapolate the U value calculated at r_paw to a
larger sphere radius. See also [[macro_uj]]. As most projector functions are
localized within r_paw to ≈80%, 20 a.u. contains ≈100% of the wavefunction and
corresponds to r_paw -> ∞.
""",
),

Variable(
    abivarname="pawujv",
    varset="dev",
    vartype="real",
    topics=['DFT+U_expert'],
    dimensions="scalar",
    defaultval=ValueWithUnit(units='eV', value=0.1),
    mnemonics="PAW+macro_UJ, potential shift (V)",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Amplitude of the potential shift for the determination of U (or J). See also [[macro_uj]].
""",
),

Variable(
    abivarname="pawusecp",
    varset="paw",
    vartype="integer",
    topics=['PAW_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="PAW - option for the USE of CPrj in memory (cprj=WF projected with NL projector)",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
When PAW is activated, the computation of cprj arrays is memory and time consuming.
When [[pawusecp]] = 0, then the cprj are never kept in memory, they are
recomputed when needed (this is CPU-time consuming). When [[pawusecp]] = 1, then
the cprj are computed once and then kept in memory.
Change the value of the keyword only if you are an experienced user (developer).
Remember:

$$ cprj = \langle\tilde{\psi}_{m}.p_{i} \rangle $$

with $\tilde{\psi}_{n}$ is the wave function, and $p_{i}$ the non-local projector.

For the time being, only activated for RF calculations.
""",
),

Variable(
    abivarname="pawxcdev",
    varset="paw",
    vartype="integer",
    topics=['PAW_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="PAW - choice for eXchange-Correlation DEVelopment (spherical part)",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
  * If set to 0, the exchange-correlation term in the spherical part of energy is totally computed on the angular mesh (time consuming but exact!)
  * If set to 1, the exchange-correlation term in the spherical part of energy is developed onto lm-moments at order 1
  * If set to 2, the exchange-correlation term in the spherical part of energy is developed onto lm-moments at order 2

Be careful: Response function (DFPT) + PAW + GGA requires [[pawxcdev]] = 0. But if you plan to do DFPT calculations, it is better to use this option also in the preliminary ground state calculation.
""",
),

Variable(
    abivarname="ph_intmeth",
    varset="eph",
    vartype="integer",
    topics=['q-points_useful'],
    dimensions="scalar",
    defaultval=2,
    mnemonics="PHonons: INTegration METHod",
    added_in_version="before_v9",
    text=r"""
Select the integration technique for computing the phonon DOS and the
Eliashberg function $\alpha^2F(\omega)$.

 * 1 --> Gaussian scheme (see also [[ph_smear]] for the broadening).
 * 2 --> Tetrahedron method (no other input is needed but at least 4 q-points in the BZ are required).
""",
),

Variable(
    abivarname="ph_ndivsm",
    varset="eph",
    vartype="integer",
    topics=['q-points_useful'],
    dimensions="scalar",
    defaultval=20,
    mnemonics="PHonons: Number of DIVisions for sampling the SMallest segment",
    added_in_version="before_v9",
    text=r"""
This variable is used in conjunction with [[ph_nqpath]] and [[ph_qpath]] to
define the q-path used for phonon band structures and phonon linewidths. It
gives the number of points used to sample the smallest segment in the q-path
specified by [[ph_qpath]].
""",
),

Variable(
    abivarname="ph_ngqpt",
    varset="eph",
    vartype="integer",
    topics=['q-points_useful'],
    dimensions=[3],
    defaultval=[20, 20, 20],
    mnemonics="PHonons: Number of Grid points for Q-PoinT mesh.",
    added_in_version="before_v9",
    text=r"""
This variable defines the q-mesh used to compute the phonon DOS and the
Eliashberg function via Fourier interpolation.
Related input variables: [[ph_qshift]] and [[ph_nqshift]].
""",
),

Variable(
    abivarname="ph_nqpath",
    varset="eph",
    vartype="integer",
    topics=['q-points_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PHonons: Number of Q-points defining the PATH",
    added_in_version="before_v9",
    text=r"""
This integer defines the number of points in the [[ph_qpath]] array.
""",
),

Variable(
    abivarname="ph_nqshift",
    varset="eph",
    vartype="integer",
    topics=['q-points_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="PHonons: Number of Q-SHIFTs",
    added_in_version="before_v9",
    text=r"""
This variable defines the number of shifts in the q-mesh used for the phonon
DOS and for the Eliashberg functions (see [[ph_ngqpt]]). If not given, the
code assumes a Gamma-centered mesh. The shifts are specified by [[ph_qshift]].
""",
),

Variable(
    abivarname="ph_qpath",
    varset="eph",
    vartype="real",
    topics=['q-points_useful'],
    dimensions=[3, 'ph_nqpath'],
    defaultval="None",
    mnemonics="Phonons: Q-PATH",
    requires="specified([[ph_nqpath]])",
    added_in_version="before_v9",
    text=r"""
This array contains the list of special q-points used to construct the q-path
used to (Fourier) interpolate phonon band structures and phonon linewidths.
See also [[ph_nqpath]] and [[ph_ndivsm]].
""",
),

Variable(
    abivarname="ph_qshift",
    varset="eph",
    vartype="real",
    topics=['q-points_useful'],
    dimensions=[3, 'ph_nqshift'],
    defaultval=[0, 0, 0],
    mnemonics="PHonons: Q-SHIFTs for mesh.",
    requires="[[ph_nqshift]]",
    added_in_version="before_v9",
    text=r"""
This array gives the shifts to be used to construct the q-mesh for computing
the phonon DOS and the Eliashberg functions (see also [[ph_nqshift]]).
If not given, a Gamma-centered mesh is used.
""",
),

Variable(
    abivarname="ph_smear",
    varset="eph",
    vartype="real",
    topics=['q-points_useful'],
    dimensions="scalar",
    defaultval="0.00002 Hartree",
    mnemonics="PHonons: SMEARing factor",
    characteristics=['[[ENERGY]]'],
    requires="[[ph_intmeth]] == 1",
    added_in_version="before_v9",
    text=r"""
The Gaussian broadening used for the integration of the phonon DOS and the
Eliashberg function. See also [[ph_intmeth]] and [[ph_ngqpt]].
""",
),

Variable(
    abivarname="ph_wstep",
    varset="eph",
    vartype="real",
    topics=['q-points_useful'],
    dimensions="scalar",
    defaultval="0.1 meV",
    mnemonics="PHonons: frequency(W)  STEP.",
    characteristics=['[[ENERGY]]'],
    added_in_version="before_v9",
    text=r"""
The step used to generate the (linear) frequency mesh for the phonon DOS and
the Eliashberg function. The extrema of the mesh are automatically computed by the code.
""",
),

Variable(
    abivarname="pimass",
    varset="rlx",
    vartype="real",
    topics=['PIMD_useful'],
    dimensions=['[[ntypat]]'],
    defaultval="[[ntypat]]",
    mnemonics="Path Integral fictitious MASSes",
    requires="[[imgmov]] = 9 or 13",
    added_in_version="before_v9",
    text=r"""
Only relevant if [[imgmov]] = 9 or 13 (Path-Integral Molecular Dynamics).
Gives the fictitious masses ( [[cite:Marx1996]]) in atomic mass units for each kind of atom in cell. These masses are the inertial masses used in performing Path Integral Molecular
Dynamics (PIMD), they are different from the true masses ([[amu]]) used to
define the quantum spring that relates the different beads in PIMD. They can
be chosen arbitrarily, but an appropriate choice will lead the different
variables to move on the same time scale in order to optimize the sampling
efficiency of the PIMD trajectory.
If [[pitransform]] = 1 (normal mode transformation), or [[pitransform]] = 2
(staging transformation), [[pimass]] is automatically set to its optimal value.
""",
),

Variable(
    abivarname="pimd_constraint",
    varset="rlx",
    vartype="integer",
    topics=['PIMD_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Path-Integral Molecular Dynamics: CONSTRAINT to be applied on a reaction coordinate",
    requires="[[imgmov]] = 9 or 13",
    added_in_version="before_v9",
    text=r"""
Only relevant for Path-Integral Molecular Dynamics.
Selects a constraint to be applied during the PIMD trajectory. The constraint
is holonomic (it is a relation between the position variables). In practice,
the total forces applied to the atomic positions are modified so as to respect
the constraint.

To date, the available constraints are:

  * **0**: no constraint
  * **1**: _"Blue Moon Ensemble" method_.
The constraint is a linear combination of the positions of atomic centroids
(this linear combination is kept constant during the simulation).
Sum[W_i * X_i] = constant
The X_i are the coordinates of the atomic centroids. The weights W_i have to
be specified with the [[wtatcon]](3,[[natcon]],[[nconeq]]),
[[iatcon]]([[natcon]]) and [[natcon]] input parameters (where [[nconeq]] is fixed to 1).
More details on the implementation in [[cite:Komeiji2007]].
""",
),

Variable(
    abivarname="pitransform",
    varset="rlx",
    vartype="integer",
    topics=['PIMD_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Path Integral coordinate TRANSFORMation",
    added_in_version="before_v9",
    text=r"""
Only relevant if [[imgmov]] = 9 or 13 (Path-Integral Molecular Dynamics).
Coordinate transformation used in the integration of the Path Integral
Molecular Dynamics equations of motion. The transformation, with an
appropriate choice of fictitious masses ([[pimass]]), is used to force the
different modes to move on the same time scale, and thus optimize the
efficiency of the statistical sampling in the corresponding statistical
ensemble. Available with a Langevin thermostat ([[imgmov]] = 9) or with Nose-
Hoover chains ([[imgmov]] = 13). See [[cite:Tuckerman1996]].

If equal to 0, no transformation is applied (primitive coordinates).
If equal to 1, normal mode transformation (in that case, [[nimage]] must be absolutely EVEN).
If equal to 2, staging transformation.
""",
),

Variable(
    abivarname="plowan_bandf",
    varset="dev",
    vartype="integer",
    topics=['Wannier_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Projected Local Orbital WANnier functions BAND Final",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""

Gives the upper band to include in the calculation of Wannier functions
""",
),

Variable(
    abivarname="plowan_bandi",
    varset="dev",
    vartype="integer",
    topics=['Wannier_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Projected Local Orbital WANnier functions BAND Initial",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Gives the lower band to include in the calculation of Wannier functions
""",
),

Variable(
    abivarname="plowan_compute",
    varset="dev",
    vartype="integer",
    topics=['Wannier_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Projected Local Orbital WANnier functions COMPUTATION",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""

Activate computation of Projected Local Orbital Wannier functions (PLO
Wannier) and corresponding band structure. Variables [[plowan_bandi]],
[[plowan_bandf]], [[plowan_natom]], [[plowan_nbl]], [[plowan_iatom]],
[[plowan_lcalc]], [[plowan_projcalc]] are mandatory to precise the nature of
the projections.

  * 0 --> Default value: do not activate calculation of PLO Wannier.
  * 1 --> Compute PLO Wannier and band structure
  * 2 --> Compute PLO Wannier and band structure. In this case, the
  coupling in k-space between blocks of Wannier functions belonging to
  different angular momenta or atoms is removed.

Other related variables are [[plowan_realspace]], [[plowan_nt]],
[[plowan_it]]. The implementation is not symmetrized over k-point and not
parallelized. (The calculation of projections is detailed in [[cite:Amadon2008]] )
""",
),

Variable(
    abivarname="plowan_iatom",
    varset="dev",
    vartype="integer",
    topics=['Wannier_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Projected Local Orbital WANnier functions, Index of ATOM",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""

Gives the indices of the [[plowan_natom]] atoms on which the projections will be done.
""",
),

Variable(
    abivarname="plowan_it",
    varset="dev",
    vartype="integer",
    topics=['Wannier_useful'],
    dimensions=[3, '[[plowan_nt]]'],
    defaultval=0,
    mnemonics="Projected Local Orbital WANnier functions,  Index of Translation.",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Requires [[plowan_realspace]] to be greater than 0 and [[plowan_nt]] to be
greater than 0. Precise a given set of selected real space translation by
using the real space vectors basis. These atoms are used to define Wannier
functions in real space. These real space Wannier functions are used as a
basis to compute the Hamiltonian.
""",
),

Variable(
    abivarname="plowan_lcalc",
    varset="dev",
    vartype="integer",
    topics=['Wannier_compulsory'],
    dimensions=['sum([[plowan_nbl]])'],
    defaultval=-1,
    mnemonics="Projected Local Orbital WANnier functions,  L values to use for CALCulation",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Gives the [[plowan_nbl]] values of angular momenta for each atom, in the order
of the atoms as given in [[plowan_iatom]].
""",
),

Variable(
    abivarname="plowan_natom",
    varset="dev",
    vartype="integer",
    topics=['Wannier_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Projected Local Orbital WANnier functions, Number of ATOMs",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Gives the number of atoms on which the projection will be done
""",
),

Variable(
    abivarname="plowan_nbl",
    varset="dev",
    vartype="integer",
    topics=['Wannier_compulsory'],
    dimensions=['[[plowan_natom]]'],
    defaultval=0,
    mnemonics="Projected Local Orbital WANnier functions,  NumBer of L values",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Gives the total number of angular momenta (over all atoms) to compute the projections.
""",
),

Variable(
    abivarname="plowan_nt",
    varset="dev",
    vartype="integer",
    topics=['Wannier_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="""Projected Local Orbital WANnier functions,  Number of Translation on which the real space values of
energy are computed""",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Requires [[plowan_realspace]] to be greater than 0. Gives a number of selected
atoms. These atoms are used to define Wannier functions in real space. These
real space Wannier functions are used as a basis to compute the Hamiltonian.
""",
),

Variable(
    abivarname="plowan_projcalc",
    varset="dev",
    vartype="integer",
    topics=['Wannier_compulsory'],
    dimensions=['sum([[plowan_nbl]])'],
    defaultval=-1,
    mnemonics="Projected Local Orbital WANnier functions,  PROJectors values to use for CALCulation",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Gives the [[plowan_nbl]] values of projectors for each atom, in the order of
the atoms as given in [[plowan_iatom]]. The index i for the projectors refers
to the i-th number on line orbitals of the PAW atomic data file.
""",
),

Variable(
    abivarname="plowan_realspace",
    varset="dev",
    vartype="integer",
    topics=['Wannier_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Projected Local Orbital WANnier functions,  activate REAL SPACE calculation.",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Can take the following values:

  * 0 --> Default value: do not activate calculation of real space Wannier functions.
  * 1 --> Compute PLO Wannier in real space for analysis. These data can also be used
    in a following dataset to perform a Wannier interpolation.
  * 2 --> Do simple Wannier Interpolation for a given k points starting from real space Wannier function
    Hamiltonian computed in a preceding dataset.
""",
),

Variable(
    abivarname="polcen",
    varset="ffield",
    vartype="real",
    topics=['Berry_useful'],
    dimensions=[3],
    defaultval=MultipleValue(number=3, value=0),
    mnemonics="POLarization for CENtrosymmetric geometry",
    added_in_version="before_v9",
    text=r"""
When doing a finite electric displacement field calculation, if the structure
is centrosymmetric but the polarization is non-zero (such as for AlAs), this
non-zero polarization should be specified as [[polcen]] (in REDUCED
coordinates, in atomic units) in the input file. See Eq.(24) in the Suppl. of [[cite:Stengel2009]]
""",
),

Variable(
    abivarname="posdoppler",
    varset="gstate",
    vartype="integer",
    topics=['positron_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="POSitron computation of DOPPLER broadening",
    added_in_version="before_v9",
    text=r"""
Relevant only when [[positron]]/=0.
This input parameter activates the calculation of the Doppler broadening of
the electron-positron annihilation radiation.
An output file containing the momentum distributions of annihilating electron-positron pairs is created.
Such a computation needs a core wave-function file (per atom type) to be
provided. This core WF file should be named '**psp_file_name**.corewf' (where
**pspfile_name** is the name of the pseudo-potential (or PAW) file) or
'corewf.abinit**ityp**' (where **ityp** is the index of the atom type). Core WF
files can be obtained with the atompaw tool by the use of **prtcorewf** keyword.
""",
),

Variable(
    abivarname="positron",
    varset="gstate",
    vartype="integer",
    topics=['positron_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="POSITRON calculation",
    added_in_version="before_v9",
    text=r"""
This input parameter can be positive or negative.
Negative values for [[positron]] are only relevant for PAW calculations.
Electron-positron correlation functional is defined by [[ixcpositron]].
Other relevant input parameter: [[posocc]] (occupation number for the
positron).

Positive values for [[positron]]:
For [[positron]] = 1 or 2, will perform the calculation of positron
lifetime (and annihilation rate).

  * [[positron]] = 1:\n
Starting from a previous electronic GS density (with [[positron]] = 0), a
positronic ground-state calculation is performed, considering that the
electrons are not perturbed by the presence of the positron.
This is almost correct for a positron in a perfect bulk material. But this
approximation fails when defects exist in the material (for instance: the
positron might be trapped by a vacancy).
The electronic density will be automatically read from a _DEN file (with or
without [[irdden]] keyword).
At the end of the SCF cycle, the positron lifetime and annihilation rate are
printed out.

Additional information for the use of pseudopotentials:

    PAW datasets: nothing to do; simply use usual electronic PAW datasets
    Norm-conserving pseudopotentials: One has to use specific pseudopotentials for the positron calculation. They must be of the FHI type (pspcod=6), and must contain at their end, the all-electrons core density generated with FHI98PP. They must have lmax=lloc=0 (check that this works for the electronic GS !! No ghost, etc...). Otherwise, their are similar to an usual FHI pseudopotential.


  * [[positron]] = 2:\n
Starting from a previous positronic GS density (with [[positron]] = 1 ), an
electronic ground-state calculation is performed, keeping the positronic
density constant.
The positronic density will be automatically read from a _DEN file (with or
without [[getden]]/[[irdden]] keyword).
At the end of the SCF cycle, the positron lifetime and annihilation rate are
printed out.

Additional information for the use of pseudopotentials:

    PAW datasets: nothing to do; simply use usual electronic PAW datasets
    Norm-conserving pseudopotentials: One has to use specific pseudopotentials for the electron calculation. They must be of the FHI type (pspcod=6), and must contain at their end, the all-electrons core density generated with FHI98PP.


  * **Typical use**:\n
The calculation is done in several steps:
The first one is a normal GS calculation for the electrons, with [[positron]] = 0.
The only specific thing to do is to set [[prtden]] = 1 (this is the default
for ABINIT v6.x+). This will create the associated _DEN file which will be
used as input file for the positronic GS calculation.
The second step is the GS calculation of the positron and subsequently its
lifetime, with [[positron]] =1. One has to define also [[ixcpositron]].
Then, it is possible to perform an additional step, computing the GS
electronic density in presence of the positron, with [[positron]] = 2 and so on...
This procedure can be automated (for PAW only) by the use of a negative value
for [[positron]].
At the end, a converged value of the positron lifetime (decomposed in several
contributions) is printed.
See also [[posdoppler]] keyword for the calculation of Doppler broadening.


Negative values for [[positron]]:\n
For [[positron]]<0, will perform an automatic calculation of electrons and
positron densities in the two-component DFT context; then will compute
positron lifetime (and annihilation rate).

  * [[positron]] = -1:\n
Starting from scratch, will first perform a usual electronic ground-state
calculation until convergence (controlled by the use of one of the _tolerance_
keywords).
Then will perform a positronic ground state calculation in presence of the
electrons and ions; then an electronic ground state calculation in presence of
the positron and the ions and so on until the total energy is converged.
The convergence of the total energy of the ions+electrons+positron system is
controlled by the use of the [[postoldfe]], [[postoldff]] and [[posnstep]]
input keywords.
With [[positron]] = -1, at the beginning of each new electronic/positronic
step, the wave functions are unknown.

  * [[positron]] = -10:\n
Same as [[positron]] = -1 except that the electronic/positronic wave functions
are stored in memory.
Consequently, the total number of iterations to reach the convergence
($\Delta$Etotal<[[postoldfe]] or $\Delta$Forces<[[postoldff]]) is smaller.
But, this can increase the total amount of memory needed by the code.

  * [[positron]] = -2:\n
Same as [[positron]] = -1 except that the two-component DFT cycle is forced to
stop at the end of an electronic step.

  * [[positron]] = -20:\n
Same as [[positron]] = -10 except that the two-component DFT cycle is forced to
stop at the end of an electronic step.


Advice for use:
There are two typical cases which have to be differently treated:

  * **A positron in a perfect _bulk_ system**:\n
In that case, the positron is delocalized in the whole crystal. Its density is
almost zero.
Thus, the "zero density positron limit" has to be used. [[ixcpositron]] has to
be chosen accordingly.
In order to have the zero density positron limit it is advised to follow these
points:

  * 1- Put a small positronic charge (by setting a [[posocc]] to a small value) **OR** use a big supercell.
  * 2- Use only k=$\Gamma$ wave vector for the positronic calculation.
  * 3- Use the manual procedure in 2 steps: first [[positron]] = 0 and then [[positron]] = 1; avoid the [[positron]] = 2 step and the automatic procedure ([[positron]]<0).

In principle, the positron lifetime should converge with the value of
[[posocc]] or the size of the supercell.

  * **A positron trapped in a _default_ (vacancy)**:\n
In that case, the positron is localized in the default. Its density can be
localized in the simulation cell (provided that the cell is sufficiently
large) and influences the electronic density.
So, it is advised to use the automatic procedure ([[positron]]<0) or the
manual procedure with several [[positron]] = 0,1,2,1,... steps.
K-points can be used as in usual electronic calculations.
Also note that it is possible to use forces and stresses to perform structural
minimization.

References: [[cite:Arponen1979a]], [[cite:Boronski1986]], [[cite:Sterne1991]], [[cite:Puska1995]], [[cite:Barbiellini1995]]
""",
),

Variable(
    abivarname="posnstep",
    varset="gstate",
    vartype="integer",
    topics=['positron_basic'],
    dimensions="scalar",
    defaultval=50,
    mnemonics="POSitron calculation: max. Number of STEPs for the two-component DFT",
    added_in_version="before_v9",
    text=r"""
Relevant only when [[positron]]<0.
Sets the maximum number of electronic/positronic iterations that, when
reached, will cause the two-component DFT SCF cycle to stop.
The code will first compute the electronic ground-state, then the positronic
ground state in the electronic density, then the electronic ground state in
the positronic density until $\Delta$Etotal<[[postoldfe]] or $\Delta$Forces<[[postoldff]] or the number
of electronic/positronic steps is [[posnstep]].
""",
),

Variable(
    abivarname="posocc",
    varset="gstate",
    vartype="real",
    topics=['positron_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="POSitron calculation: OCCupation number for the positron",
    added_in_version="before_v9",
    text=r"""
Relevant only when [[positron]]/=0.
Sets the occupation number for the positron. Has to be <=1.
Changing [[posocc]] is only useful for bulk calculation when one wants to
perform lifetime computations using a small simulation cell (can avoid the use
of a supercell). It simulates the dispersion of the positron in the whole
crystal.
""",
),

Variable(
    abivarname="postoldfe",
    varset="gstate",
    vartype="real",
    topics=['positron_basic'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[postoldff]] = 0': 1e-06, 'defaultval': 0.0}),
    mnemonics="POSitron calculation: TOLerance on the DiFference of total Energy",
    characteristics=['[[ENERGY]]'],
    added_in_version="before_v9",
    text=r"""
Relevant only when [[positron]]<0.
Sets a tolerance for absolute difference of total energy (of
_ions+electrons+positron_ system) that, when reached, will cause the SCF cycle
to stop before the number of steps is [[nstep]] or the number of
electronic/positronic steps is [[posnstep]].

Can be specified in Ha (the default), Ry, eV or Kelvin, since [[postoldfe]] has
the [[ENERGY]] characteristics.
One and only one of [[postoldfe]] or [[postoldff]] can be set.
""",
),

Variable(
    abivarname="postoldff",
    varset="gstate",
    vartype="real",
    topics=['positron_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="POSitron calculation: TOLerance on the DiFference of Forces",
    added_in_version="before_v9",
    text=r"""
Relevant only when [[positron]] < 0.
Sets a tolerance for absolute difference of maximum force (in hartree/Bohr) acting on ions (due
to _ions+electrons+positron_ system) that, when reached, will cause the SCF
cycle to stop before the number of SCF steps is [[nstep]] or the number of
electronic/positronic steps is [[posnstep]].
One and only one of [[postoldfe]] or [[postoldff]] can be set.
""",
),

Variable(
    abivarname="ppmfrq",
    varset="gw",
    vartype="real",
    topics=['SelfEnergy_basic'],
    dimensions="scalar",
    defaultval=ValueWithUnit(units='Ha', value=0.0),
    mnemonics="Plasmon Pole Model FReQuency",
    characteristics=['[[ENERGY]]'],
    requires="[[optdriver]] in [3,4]",
    added_in_version="before_v9",
    text=r"""
**In plasmon-pole calculations**

Usually only effective if GW corrections are evaluated using the plasmon-pole
model of Godby-Needs ([[ppmodel]] == 1).

In the present status of the GW code, the convolution in frequency space
defining the self-energy operator can be evaluated using two different
approaches: numerical integration and plasmon-pole models.
Methods based on the numerical integration (contour deformation, analytic
continuation) require the knowledge of the screened interaction for several
frequencies. These approaches give the most accurate results but at the price
of an increase in the CPU time required.
Alternatively, it is possible to approximate the dynamical behaviour of the
screened interaction through simple analytical expressions, the so-called
plasmon-pole models. In the plasmon-pole model proposed by Godby-Needs
([[ppmodel]] = 1), the screening must be available at zero frequency, as well as
at another imaginary frequency, of the order of the plasmon frequency (the
peak in the EELS spectrum). This information is used to model the behaviour of
the dielectric matrix for all frequencies. During the calculation of the
screening, [[ppmfrq]] defines the imaginary frequency where the dielectric
matrix is evaluated, in addition to the zero frequency. During the self-energy
run, [[ppmfrq]] can be used to define the second frequency to be used to
calculate the plasmon-pole parameters. This is particularly useful when the
SCR file contains several frequencies along the imaginary axis. In this case
the frequency whose value is the closest one to [[ppmfrq]] will be selected.
Note that, if the plasmon-pole approximation is good, then, the choice of
[[ppmfrq]] should have no influence on the final result. One should check
whether this is the case. In general, the plasmon frequencies of bulk solids
are of the order of 0.5 Hartree.

**In Contour Deformation calculations**

[[ppmfrq]] is here used to **override** the default value calculated from the
average electronic density per unit cell. This can affect the distribution of
gridpoints along the imaginary and real frequency axes. See
[[cd_frqim_method]], [[gw_frqim_inzgrid]] and [[gw_frqre_inzgrid]] for more
details.
""",
),

Variable(
    abivarname="ppmodel",
    varset="gw",
    vartype="integer",
    topics=['SelfEnergy_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Plasmon Pole MODEL",
    requires="[[optdriver]] in [3,4]",
    added_in_version="before_v9",
    text=r"""
  * **ppmodel** = 1: PP model of Godby and Needs [[cite:Godby1989]].
  * **ppmodel** = 2: PP model of Hybertsen and Louie [[cite:Hybertsen1986]].
  * **ppmodel** = 3: PP model of W. von der Linden and P. Horsh [[cite:vonderLinden1988]].
  * **ppmodel** = 4: PP model of Farid and Engel [[cite:Engel1993]].
  * **ppmodel** = 0: no PP model, numerical integration (contour deformation method [[cite:Lebegue2003]]).

Please note the difference between **ppmodel** 1 and **ppmodel** 2,3,4. In the
first case (**ppmodel** = 1), the plasmon-pole parameters are determined in
order to reproduce the behaviour of the dielectric matrix at two calculated
frequencies: the static limit ($\omega=0$) and the imaginary frequency defined by
[[ppmfrq]]. In the last three cases, instead, the plasmon-pole parameters are
found by using the dielectric matrix calculated only at $\omega=0$ and enforcing
the so-called f-sum rule. See also [[nfreqre]].

Please note also that in the case of **ppmodel** 4, the plasmon energies are
not simple mathematical parameters, but rather have a physical meaning (at
least the lowest ones). Thus the calculated plasmon band structure (plasmon
energy vs q vector) is reported in the output file for the lowest 10 bands.
""",
),

Variable(
    abivarname="prepanl",
    varset="dfpt",
    vartype="integer",
    topics=['nonlinear_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PREPAre Non-Linear response calculation",
    added_in_version="before_v9",
    text=r"""
The computation of third-order derivatives from the 2n+1 theorem requires the
first-order wavefunctions and densities obtained from a linear response
calculation. The standard approach in a linear response calculation is:

  * compute only the irreducible perturbations;
  * use symmetries to reduce the number of k-points for the k-point integration.

This approach cannot be applied, presently (v4.1), if the first-order
wavefunctions are to be used to compute third-order derivatives. First, for
electric fields, the code needs the derivatives along the three directions.
Still, in case of phonons, only the irreducible perturbations are required.
Second, for both electric fields and phonons, the wavefunctions must be
available in half the BZ (kptopt=2), or the full BZ (kptopt=3).
During the linear response calculation, in order to prepare a non-linear
calculation, one should put [[prepanl]] to 1 in order to force ABINIT to
compute the electric field perturbation along the three directions explicitly,
and to keep the full number of k-points.

In the case of a 2nd derivative of wavefunction ([[rf2_dkdk]] or [[rf2_dkde]]),
[[prepanl]] == 1 can be used in order to skip directions of perturbations that
will not be used by the non-linear routine (see [[rf2_dkdk]] for more details).
""",
),

Variable(
    abivarname="prepgkk",
    varset="dfpt",
    vartype="integer",
    topics=['ElPhonInt_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PREPAre GKK calculation",
    added_in_version="before_v9",
    text=r"""
The calculation of electron-phonon coupling quantities requires the presence
of all the perturbations (all atoms in all directions) for the chosen set of
(irreducible) q-points. To impose this and prevent ABINIT from using symmetry
to reduce the number of perturbations, set [[prepgkk]] to 1. Use in
conjunction with [[prtgkk]].
""",
),

Variable(
    abivarname="prepscphon",
    varset="dev",
    vartype="integer",
    topics=['printing_prngs'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PREPare Self-Consistent PHONon calculation",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Print PCINFO, PHFREQ, and PHVEC files, for use with self-consistent phonon
runs, after a perturbation calculation. Only prints out files for the present
q-point, and there is presently no tool to symmetrize or merge these files, so
use anaddb instead (with prtscphon input variable). The abinit input variable
is destined to someday bypass the use of anaddb for scphon calculations.
""",
),

Variable(
    abivarname="prt1dm",
    varset="files",
    vartype="integer",
    topics=['printing_prgs'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT 1-DiMensional potential and density",
    added_in_version="before_v9",
    text=r"""
If set >= 1, provide one-dimensional projection of potential and density, for
each of the three axis. This corresponds to averaging the potential or the
density on bi-dimensional slices of the FFT grid.
""",
),

Variable(
    abivarname="prtatlist",
    varset="rlx",
    vartype="integer",
    topics=['printing_prgeo', 'Output_useful'],
    dimensions=['[[natom]]'],
    defaultval=0,
    mnemonics="PRinT by ATom LIST of ATom",
    characteristics=['[[NO_MULTI]]'],
    added_in_version="before_v9",
    text=r"""
This is an array of the numbers associated to the index atoms that the user
want to print in the output or log files, this is useful when you have a large
number of atoms and you are only interested to follow specific atoms, the
numbers associated should be consistent with the list in [[xcart]] or
[[xred]]. This input variable does not affect the contents of the "OUT.nc" or
"HIST.nc", those are NetCDF files containing the information about all the atoms.
""",
),

Variable(
    abivarname="prtbbb",
    varset="dfpt",
    vartype="integer",
    topics=['printing_prngs', 'Output_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT Band-By-Band decomposition",
    added_in_version="before_v9",
    text=r"""
If [[prtbbb]] is 1, print the band-by-band decomposition of Born effective
charges and localization tensor, in case they are computed. See [[cite:Ghosez2000]].
""",
),

Variable(
    abivarname="prtbltztrp",
    varset="dev",
    vartype="integer",
    topics=['printing_prgs'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT output for BoLTZTRaP code",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Print out geometry _BLZTRP_GEOM and eigenenergy _BLZTRP_EIGEN files for
the [BoltzTraP code](https://www.imc.tuwien.ac.at/forschungsbereich_theoretische_chemie/forschungsgruppen/prof_dr_gkh_madsen_theoretical_materials_chemistry/boltztrap/) by Georg Madsen.
""",
),

Variable(
    abivarname="prtcif",
    varset="dev",
    vartype="integer",
    topics=['printing_prgeo'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT Crystallographic Information File",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
If set to 1, a CIF file is output with the crystallographic data for the
present run (cell size shape and atomic positions).
""",
),

Variable(
    abivarname="prtden",
    varset="files",
    vartype="integer",
    topics=['printing_prden'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[nimage]] > 1': 0, 'defaultval': 1}),
    mnemonics="PRinT the DENsity",
    added_in_version="before_v9",
    text=r"""
If set to 1 or a larger value, provide output of electron density in real
space rho(r), in units of electrons/Bohr^3.
If [[ionmov]] == 0, the name of the density file will be the root output name,
followed by _DEN.
If [[ionmov]] /= 0, density files will be output at each time step, with
the name being made of

  * the root output name,
  * followed by _TIMx, where x is related to the time step (see later)
  * then followed by _DEN

The file structure of this unformatted output file is described in [[help:abinit#denfile|this section]].
If [[prtden]] is lower than 0, two files will be printed for restart every
[[prtden]] step, with the names being made of

  * the root temporary name,
  * followed by _DEN_x, where x is 0000 or 0001 alternatively.
  * The most recent of the two files should be used for restart, and copied to root input name_DS2_DEN
  * To perform a restart, in a multidataset mode, use ndtset 2 and jdtset 2 3 (that is 2 datasets, numbered 2 and 3)
  * In the dataset 2, get the density you just copied (getden2 -1), perform a non self-consistent calculation and print the wave function (prtwf2 1)
  * In the dataset 3, get the previous wf(getwfk3 -1), and continue the calculation
  * This complicated procedure is due to the fact that reading the density is only allowed for a non sc calculation, and also for a dataset different of 0 or the previous one, the option we choose here.

Please note that in the case of PAW ([[usepaw]] = 1) calculations, the _DEN
density output is not the full physical electron density. If what is wanted is
the full physical electron density, say for post-processing with [[help:aim|AIM]]
or visualization, prtden > 1 will produce physical electron density or other interesting quantities (see below).
Nevertheless, even in the PAW case, when chaining together calculations where
the density from one calculation is to be used in a subsequent calculation, it
is necessary to use the _DEN files and **not** one of the other files produced
with prtden  > 1, i.e. _PAWDEN, ATMDEN_xxx or else. Note that the usual _DEN
file is always generated as soon as prtden >= 1. Options 2 to 6 for prtden are
relevant only for [[usepaw]] = 1 and control the output of the full electron
density in the PAW case:


**prtden=2** causes generation of a file _PAWDEN that contains the bulk
**valence** charge density together with the PAW on-site contributions, and
has the same format as the other density files.
**prtden=3** causes generation of a file _PAWDEN that contains the bulk
**full** charge density (valence+core)
**prtden=4** causes generation of three files _ATMDEN_CORE, _ATMDEN_VAL and
_ATMDEN_FULL which respectively contain the core, valence and full atomic
protodensity (the density of the individual component atoms in vacuum
superposed at the bulk atomic positions). This can be used to generate various
visualizations of the bonding density.
**prtden=5** options 2 and 4 taken together.
**prtden=6** options 3 and 4 taken together.
**prtden=7** causes the generation of all the individual contributions to the
bulk **valence** charge density: n_tilde-n_hat (_N_TILDE), n_onsite (_N_ONE)
and n_tilde_onsite (_NT_ONE). This is for diagnosis purposes only.

Options 3 to 6 currently require the user to supply the atomic core and
valence density in external files in the working directory. The files must be
named properly; for example, the files for an atom of type 1 should be named:
"core_density_atom_type1.dat" and "valence_density_atom_type1.dat". The file
should be a text file, where the first line is assumed to be a comment, and
the subsequent lines contain two values each, where the first one is a radial
coordinate and the second the value of the density n(r). Please note that it
is n(r) which should be supplied, **not** n(r)/r^2. The first coordinate point
must be the origin, i.e. **_r = 0_ **. The atomic densities are spherically
averaged, so assumed to be completely spherically symmetric, even for open
shells.

NOTE: in the PAW case, **DO NOT** use _PAWDEN or _ATMDEN_xxx files produced by
prtden  > 1 to chain the density output from one calculation as the input to
another, use the _DEN file for that.
""",
),

Variable(
    abivarname="prtdensph",
    varset="gstate",
    vartype="integer",
    topics=['printing_prden'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'defaultval': 1}),
    mnemonics="PRinT integral of DENsity inside atomic SPHeres",
    added_in_version="before_v9",
    text=r"""
When this flag is activated, values of integral(s) of total density inside
sphere(s) around each atom are printed in output file (for each spin
component). Spheres around atoms are defined by a radius given by [[ratsph]] keyword.
Note: integral of density inside a sphere around an atom can be used to
determine a rough approximation of the local magnetic moment; this is
particularly useful for antiferromagnetic systems.
The algorithm to compute this integral is particularly primitive: the points
on the FFT grids, belonging to the interior of the sphere are determined, and
the value of the functions on these points are summed, taking into account a
fixed volume attributed to each point. In particular, the integral as a
function of the radius will be a constant, except when a new point enters the
sphere, in which case a sudden jump occurs. However, since the purpose of this
output is to get a rough idea of the repartition of the density, this is not a
real problem. If you are interested in a more accurate estimation of the
density within a sphere, you should use the cut3d postprocessor ([[help:cut3d]]).
""",
),

Variable(
    abivarname="prtdipole",
    varset="dev",
    vartype="integer",
    topics=['printing_prgs'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT DIPOLE",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Print out dipole of unit cell, calculated in real space for the primitive cell
only. Under development.
""",
),

Variable(
    abivarname="prtdos",
    varset="files",
    vartype="integer",
    topics=['printing_prdos', 'ElecDOS_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT the Density Of States",
    added_in_version="before_v9",
    text=r"""
Provide output of Density of States if set to 1, 2 or 3. Can either use a
smearing technique ([[prtdos]] = 1), or the tetrahedron method ([[prtdos]] = 2).
If [[prtdos]] = 3, provide output of Local Density of States inside a sphere
centered on an atom, as well as the angular-momentum projected DOS, in the
same sphere. The resolution of the linear grid of energies for which the DOS
is computed can be tuned thanks to [[dosdeltae]].

If [[prtdos]] = 1, the smeared density of states is obtained from the
eigenvalues, properly weighted at each k point using [[wtk]], and smeared
according to [[occopt]] and [[tsmear]]. All levels that are present in the
calculation are taken into account (occupied and unoccupied).
In order to compute the DOS of an insulator with [[prtdos]] = 1, compute its
density thanks to a self-consistent calculation (with a non-metallic
[[occopt]] value, 0, 1 or 2), then use [[prtdos]] = 1, together with
[[iscf]] = -3, and a metallic [[occopt]], between 3 and 7, providing the needed
smearing. If [[prtdos]] = 1, the name of the DOS file is the root name for the
output files, followed by "_DOS".

 * Note 1: [[occopt]] must be between 3 and 7.
 * Note 2: The sampling of the Brillouin Zone that is needed to get a converged DOS
 is usually much finer than the sampling needed to converge the total energy or the geometry of the
system, unless [[tsmear]] is very large (hence the DOS is not obtained
properly). A separate convergence study is needed.


If [[prtdos]] = 2, the DOS is computed using the tetrahedron method. As in the
case of [[prtdos]] = 1, all levels that are present in the calculation are taken
into account (occupied and unoccupied). In this case, the k-points must have
been defined using the input variable [[ngkpt]] or the input variable
[[kptrlatt]]. There must be at least two non-equivalent points in the
Irreducible Brillouin Zone to use [[prtdos]] = 2. It is strongly advised that you use
a non-shifted k-point grid ([[shiftk]] 0 0 0): such grids naturally contain
more extremal points (band minima and maxima at Gamma or at the zone-
boundaries) than shifted grids, and lead to more non-equivalent points than
shifted grids, for the same grid spacing. There is no need to take care of the
[[occopt]] or [[tsmear]] input variables, and there is no subtlety to be taken
into account for insulators. The computation can be done in the self-
consistent case as well as in the non-self-consistent case, using [[iscf]] = -3.
This allows one to refine the DOS at fixed starting density.
In that case, if [[ionmov]] == 0, the name of the potential file will be the
root output name, followed by _DOS (like in the [[prtdos]] = 1 case).
However, if [[ionmov]] /= 0, potential files will be output at each time
step, with the name being made of

  * the root output name,
  * followed by _TIMx, where x is related to the time step (see later)
  * then followed by _DOS.

If [[prtdos]] = 3, the same tetrahedron method as for [[prtdos]] = 2 is used, but
the DOS inside a sphere centered on some atom is delivered, as well as the
angular-momentum projected (l=0,1,2,3,4) DOS in the same sphere. The
preparation of this case, the parameters under which the computation is to be
done, and the file denomination is similar to the [[prtdos]] = 2 case. However,
three additional input variables might be provided, describing the atoms that
are the center of the sphere (input variables [[natsph]] and [[iatsph]]), as
well as the radius of this sphere (input variable [[ratsph]]).
In case of PAW, [[ratsph]] radius has to be greater or equal to the largest PAW
radius of the atom types considered (which is read from the PAW atomic data
file; see rc_sph or r_paw). Additionally, printing and/or approximations in PAW
mode can be controlled with [[pawprtdos]] keyword (in
particular,[[pawprtdos]] = 2 can be used to compute quickly a very good
approximation of the DOS).

 * Note 1: when [[prtdos]] = 3, it is possible to output m-decomposed LDOS in _DOS
file; simply use [[prtdosm]] keyword.
 * Note 2: the integrated total DOS in spheres around atoms can be obtained when
[[prtdensph]] flag is activated. It can be compared to the integrated DOS
provided in _DOS file when [[prtdos]] = 3.

[[prtdos]] = 4 delivers the sphere-projected DOS (like [[prtdos]] = 3), on the
basis of a smearing approach (like [[prtdos]] = 1)

[[prtdos]] = 5 delivers the spin-spin DOS in the [[nspinor]] == 2 case, using the
tetrahedron method (as [[prtdos]] = 2).
""",
),

Variable(
    abivarname="prtdosm",
    varset="files",
    vartype="integer",
    topics=['printing_prdos', 'ElecDOS_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT the Density Of States with M decomposition",
    added_in_version="before_v9",
    text=r"""
Relevant only when [[prtdos]] = 3.
If set to 1, the m-decomposed LDOS is delivered in DOS file.
Note that [[prtdosm]] computes the M-resolved partial dos for complex
spherical harmonics,giving e.g. DOS(L,M) == DOS(L,-M) (without spin-orbit). In
the contrary, the LDA+U occupation matrix, see [[dmatpawu]] is in the real
spherical harmonics basis.
If set to 2, the m-decomposed LDOS is delivered in DOS file.
In this case, [[prtdosm]] computes the M-resolved partial dos for real
spherical harmonics in the same basis as the LDA+U occupation matrix.
""",
),

Variable(
    abivarname="prtebands",
    varset="gstate",
    vartype="integer",
    topics=['printing_prden'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[nimage]] > 1': 0, 'defaultval': 1}),
    mnemonics="PRinT Electron BANDS",
    added_in_version="before_v9",
    text=r"""
This option activates the output of the electron eigenvalues. Possible values:

  * 0- Disable the output of the band energies.
  * 1- Write eigenvalues in xmgrace format. A file with extension `EBANDS.agr` is produced at the end of the run.
    Use `xmgrace file_EBANDS.agr` to visualize the band energies
  * 2- Write eigenvalues in gnuplot format. The code produces a `EBANDS.dat` file with the eigenvalues
    and a `file_EBANDS.gnuplot` script. Use `gnuplot file_EBANDS.gnuplot` to visualize the band energies.
""",
),

Variable(
    abivarname="prtefg",
    varset="paw",
    vartype="integer",
    topics=['printing_prngs', 'EFG_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRint Electric Field Gradient",
    requires="[[usepaw]] == 1, [[quadmom]]",
    added_in_version="before_v9",
    text=r"""
If nonzero, calculate the electric field gradient at each atomic site in the unit cell.
Using this option requires [[quadmom]] to be set as well.
Values will be written to main output file (search for Electric Field Gradient).
If prtefg=1, only the quadrupole coupling in MHz and asymmetry are reported.
If prtefg=2, the full electric field gradient tensors in atomic units are also given,
showing separate contributions from the valence electrons, the ion cores, and the PAW reconstruction.
If prtefg=3, then in addition to the prtefg=2 output, the EFGs are computed using an ionic point charge model.
This is useful for comparing the accurate PAW-based results to those of simple ion-only models.
Use of prtefg=3 requires that the variable [[ptcharge]] be set as well.
The option prtefg is compatible with spin polarized calculations (see
[[nspden]]) and also LDA+U (see [[usepawu]]).
""",
),

Variable(
    abivarname="prtefmas",
    varset="dfpt",
    vartype="integer",
    topics=['printing_prngs', 'EffectiveMass_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="PRint EFfective MASs data",
    requires="[[efmas]] == 1",
    added_in_version="before_v9",
    text=r"""
If 1, at the end of an effective mass calculation ([[efmas]] = 1), create a file *_EFMAS, that
contains the generalized second-order k-derivatives, see Eq.(66) in [[cite:Laflamme2016]],
in view of further processing.
""",
),

Variable(
    abivarname="prteig",
    varset="files",
    vartype="integer",
    topics=['printing_prden'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[nimage]] > 1': 0, 'defaultval': 1}),
    mnemonics="PRinT EIGenenergies",
    added_in_version="before_v9",
    text=r"""
If set to 1, a file *_EIG, containing the k-points and one-electron eigenvalues is printed.
""",
),

Variable(
    abivarname="prtelf",
    varset="files",
    vartype="integer",
    topics=['printing_prgs'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT Electron Localization Function (ELF)",
    added_in_version="before_v9",
    text=r"""
If set to 1 or a larger value, provide output of ELF in real space elf(r).
This is a dimensionless quantity bounded between 0 and 1.
The name of the ELF file will be the root output name, followed by _ELF.
Like a _DEN file, it can be analyzed by cut3d. However unlike densities, in
case of spin polarized calculations, the spin down component can not be
obtained by subtracting the spin up component to the total ELF. Hence when
spin polarized calculations are performed the code produces also output files
with _ELF_UP and _ELF_DOWN extensions. (For technical reasons these files
contain also two components but the second is zero. So to perform analysis of
_ELF_UP and _ELF_DOWN files with cut3d you have to answer "ispden= 0 --> Total
density" when cut3d ask you which ispden to choose. Also remember that spin
down component can not be obtained by using cut3d on the _ELF file. Sorry for
the inconvenience, this will be fixed in the next release.)
ELF is not yet implemented in non collinear spin case.
If prtelf is set to 2, in the case of spin polarized calculation, the total
ELF is computed from an alternative approach which should better take into
account the existence of spin dependent densities (see the documentation in
/doc/theory/ELF of your ABINIT repository)

Please note that ELF is **not** yet implemented in the case of PAW
([[usepaw]] = 1) calculations.
""",
),

Variable(
    abivarname="prtfc",
    varset="paw",
    vartype="integer",
    topics=['printing_prden', 'EFG_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT Fermi Contact term",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
  * If set to 1, print the Fermi contact interaction at each nuclear site, that is, the electron density at each site. The result appears in the main output file (search for FC). Note that this calculation is different than what is done by cut3d, because it also computes the PAW on-site corrections in addition to the contribution from the valence pseudo-wavefunctions.
""",
),

Variable(
    abivarname="prtfull1wf",
    varset="dfpt",
    vartype="integer",
    topics=['DFPT_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT FULL 1st-order WaveFunction",
    added_in_version="before_v9",
    text=r"""
If set to 1, the output _1WF files will contain the full 1st-order wavefunctions, for both valence and conduction bands.
Otherwise, the _1WF files are not really 1st-order perturbed wavefunctions, but merely a set of perturbed wavefunctions that yield the correct perturbed density.
This is used when one expect to perform post-processing of the 1st-order wavefunctions.
""",
),

Variable(
    abivarname="prtfsurf",
    varset="files",
    vartype="integer",
    topics=['printing_prfermi'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT Fermi SURFace file",
    added_in_version="before_v9",
    text=r"""
If set to 1, provide Fermi surface file in the BXSF format (Xcrysden) If
[[prtfsurf]] = 1, a _BXSF file readable by [XCrySDen](http://www.xcrysden.org)
will be produced at the end of the calculation. The file contains information
on the band structure of the system and can be used to visualize the Fermi
surface or any other energy isosurface. [[prtfsurf]] = 1 is compatible only with
SCF calculations ([[iscf]] > 1) or NSCF runs in which the occupation factors
and Fermi level are recalculated once convergence is achieved ([[iscf]] = -3).
The two methods should produce the same Fermi surface provided that the
k-meshes are sufficiently dense. The k-mesh used for the sampling of the Fermi
surface can be specified using the standard variables [[ngkpt]], ([[shiftk]],
and [[nshiftk]]. Note, however, that the mesh must be homogeneous and centered
on gamma (multiple shifts are not supported by Xcrysden)
""",
),

Variable(
    abivarname="prtgden",
    varset="files",
    vartype="integer",
    topics=['printing_prden'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT the Gradient of electron DENsity",
    added_in_version="before_v9",
    text=r"""
If set to 1 or a larger value, provide output of gradient of electron density
in real space grho(r), in units of Bohr^-(5/2).
The names of the gradient of electron density files will be the root output
name, followed by _GDEN1, _GDEN2, GDEN3 for each principal direction (indeed it is a vector).
Like a _DEN file, it can be analyzed by cut3d.
The file structure of this unformatted output file is described in [[help:abinit#denfile|this section]].
""",
),

Variable(
    abivarname="prtgeo",
    varset="files",
    vartype="integer",
    topics=['printing_prgeo'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT the GEOmetry analysis",
    added_in_version="before_v9",
    text=r"""
If set to 1 or a larger value, provide output of geometrical analysis (bond
lengths and bond angles). The value of [[prtgeo]] is taken by the code to be
the maximum coordination number of atoms in the system.
It will deduce a maximum number of "nearest" and "next-nearest" neighbors
accordingly, and compute corresponding bond lengths.
It will compute bond angles for the "nearest" neighbours only.
If [[ionmov]] == 0, the name of the file will be the root output name, followed by _GEO.
If [[ionmov]] /= 0, one file will be output at each time step, with the
name being made of

  * the root output name,
  * followed by _TIMx, where x is related to the time step (see later)
  * then followed by _GEO

The content of the file should be rather self-explanatory.
No output is provided by [[prtgeo]] is lower than or equal to 0.
If [[prtgeo]] > 0, the maximum number of atoms ([[natom]]) is 9999.
""",
),

Variable(
    abivarname="prtgkk",
    varset="files",
    vartype="integer",
    topics=['printing_prngs', 'ElPhonInt_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT the GKK matrix elements file",
    added_in_version="before_v9",
    text=r"""
If set to 1, provide output of electron-phonon "gkk" matrix elements, for
further treatment by mrggkk utility or anaddb utility. Note that symmetry will
be disabled for the calculation of the perturbation, forcing the inclusion of
all k-points and all perturbation directions. Additional information on
electron-phonon treatment in ABINIT is given in the tutorial [[tutorial:eph]].
""",
),

Variable(
    abivarname="prtgsr",
    varset="files",
    vartype="integer",
    topics=['printing_prgs'],
    dimensions="scalar",
    defaultval="prtgsr = 0",
    mnemonics="PRinT the GSR file",
    added_in_version="before_v9",
    text=r"""
If set to 1, ABINIT will produce a GSR file at the end of the GS calculation.
The GSR file contains the most important GS results (band structure, forces,
stresses, electronic density). The GSR file can be read by AbiPy and used for further post-processing.
""",
),

Variable(
    abivarname="prtkden",
    varset="files",
    vartype="integer",
    topics=['printing_prden'],
    dimensions="scalar",
    defaultval="1 if [[usekden]] == 1 and [[nimage]] == 1 else 0",
    mnemonics="PRinT the Kinetic energy DENsity",
    added_in_version="before_v9",
    text=r"""
If set to 1 or a larger value, provide output of kinetic energy density in
real space tau(r), in units of Bohr^-5.
The name of the kinetic energy density file will be the root output name, followed by _KDEN.
Like a _DEN file, it can be analyzed by cut3d.
The file structure of this unformatted output file is described in [[help:abinit#denfile|this section]].
Note that the computation of the kinetic energy density must be activate,
thanks to the input variable [[usekden]].
Please note that kinetic energy density is **not** yet implemented in the case
of PAW ([[usepaw]] = 1) calculations.
""",
),

Variable(
    abivarname="prtkpt",
    varset="files",
    vartype="integer",
    topics=['printing_prden', 'Output_useful', 'k-points_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT the K-PoinTs sets",
    added_in_version="before_v9",
    text=r"""
If set /= 0, proceeds to a detailed analysis of different k point grids.
Works only if [[kptopt]] is positive, and neither [[kptrlatt]] nor [[ngkpt]]
are defined. ABINIT will stop after this analysis.

Different sets of k point grids are defined, with common values of [[shiftk]].
In each set, ABINIT increases the length of vectors of the supercell (see
[[kptrlatt]]) by integer steps. The different sets are labelled by "iset". For
each k point grid, [[kptrlen]] and [[nkpt]] are computed (the latter always
invoking [[kptopt]] = 1, that is, full use of symmetries). A series is finished
when the computed [[kptrlen]] is twice larger than the input variable
[[kptrlen]]. After the examination of the different sets, ABINIT summarizes,
for each [[nkpt]], the best possible grid, that is, the one with the largest
computed [[kptrlen]].

Note that this analysis is also performed when [[prtkpt]] = 0, as soon as
neither [[kptrlatt]] nor [[ngkpt]] are defined. But, in this case, no analysis
report is given, and the code selects the grid with the smaller [[ngkpt]] for
the desired [[kptrlen]]. However, this analysis takes some times (well
sometimes, it is only a few seconds - it depends on the value of the input
[[kptrlen]]), and it is better to examine the full analysis for a given cell
and set of symmetries, [[shiftk]] for all the production runs.

If set to -2, the code stops in invars1 after the computation of the
irreducible set and a file named kpts.nc with the list of the k-points and the
corresponding weights is produced
""",
),

Variable(
    abivarname="prtlden",
    varset="files",
    vartype="integer",
    topics=['printing_prden'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT the Laplacian of electron DENsity",
    added_in_version="before_v9",
    text=r"""
If set to 1 or a larger value, provide output of Laplacian of electron density
in real space grho(r), in units of Bohr^-(7/2).
The name of the Laplacian of electron density file will be the root output name, followed by _LDEN.
Like a _DEN file, it can be analyzed by cut3d.
The file structure of this unformatted output file is described in [[help:abinit#denfile|this section]].
""",
),

Variable(
    abivarname="prtnabla",
    varset="paw",
    vartype="integer",
    topics=['printing_prgs'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRint NABLA",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
  * If set to 1, calculate the matrix elements <Psi_n|-inabla|Psi_m> and write it in file _OPT to be read by the code conducti (see [[cite:Mazevet2010]]).
""",
),

Variable(
    abivarname="prtnest",
    varset="dev",
    vartype="integer",
    topics=['printing_prfermi'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT NESTing function",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
If set to 1, the nesting function for the k-point grid is printed. For the
moment the path in q space for the nesting function is fixed, but will become
an input as well.
""",
),

Variable(
    abivarname="prtphbands",
    varset="eph",
    vartype="integer",
    topics=['printing_prngs'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="PRinT PHonon BANDS",
    added_in_version="before_v9",
    text=r"""
This option activates the output of the phonon frequencies in the EPH code.
Possible values:

  * 0 Disable the output of the phonon frequencies.
  * 1 Write frequencies in |xmgrace| format. A file with extension `PHBANDS.agr` is produced.
    Use `xmgrace file_PHBANDS.agr` to visualize the data
  * 2 Write frequencies in |gnuplot| format. The code produces a `PHBANDS.dat` file
    with the eigenvalues and a `PHBANDS.gnuplot` script.
    Use `gnuplot file_PHBANDS.gnuplot` to visualize the phonon band structure.
""",
),

Variable(
    abivarname="prtphdos",
    varset="eph",
    vartype="integer",
    topics=['printing_prngs', 'ElPhonInt_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="PRinT the PHonon Density Of States",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Print the phonon density of states. It is activated by default when [[optdriver]] == 7.

Note also that this variable activates the computation of the generalized Eliashberg function
associated to the electron-phonon self-energy when [[eph_task]] in [-4, 4].
""",
),

Variable(
    abivarname="prtphsurf",
    varset="eph",
    vartype="integer",
    topics=['printing_prngs', 'ElPhonInt_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT PHonon iso-SURFace",
    added_in_version="before_v9",
    text=r"""
Print a bxsf file (Xcrysden format) with the (interpolated) phonon frequencies
computed of the q-mesh determined by [[ph_ngqpt]]. The file can be use to
visualize iso-surfaces with Xcrysden or other similar tools supporting the bxsf
format. Note that the (dense) q-mesh must be Gamma-centered, shifted meshes are
not supported by Xcrysden. This variable requires [[optdriver]] == 7.
""",
),

Variable(
    abivarname="prtposcar",
    varset="dev",
    vartype="integer",
    topics=['printing_prgeo'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT POSCAR file",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Print out VASP-style POSCAR and FORCES files, for use with PHON or frophon
codes for frozen phonon calculations. See the associated script in
{% modal ../scripts/post_processing/phondisp2abi.py %} for further details on
interfacing with PHON, PHONOPY, etc...
""",
),

Variable(
    abivarname="prtpot",
    varset="files",
    vartype="integer",
    topics=['printing_prpot'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT total POTential",
    added_in_version="before_v9",
    text=r"""
If set >=1, provide output of the total (Kohn-Sham) potential (sum of local
pseudo-potential, Hartree potential, and xc potential).

If [[ionmov]] == 0, the name of the potential file will be the root output name,
followed by _POT.
If [[ionmov]] /= 0, potential file will be output at each time step, with
the name being made of

  * the root output name,
  * followed by _TIMx, where x is related to the time step (see later)
  * then followed by _POT.

The file structure of this unformatted output file is described in [[help:abinit#localpotfile|this section]].
No output is provided by a negative value of this variable.
""",
),

Variable(
    abivarname="prtpsps",
    varset="files",
    vartype="integer",
    topics=['printing_prgs'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRint the PSPS file",
    added_in_version="before_v9",
    text=r"""
If set to 1, the code produces a netcdf file (PSPS.nc) with the internal
tables used by Abinit to apply the pseudopotential part of the KS Hamiltonian.
The data can be visualized with AbiPy. If prtpsps is set to -1, the code will
exit after the output of the PSPS.nc file.
""",
),

Variable(
    abivarname="prtspcur",
    varset="files",
    vartype="integer",
    topics=['printing_prgs'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT the SPin CURrent density",
    added_in_version="before_v9",
    text=r"""
If set to 1 or a larger value, provide output of the current density of
different direction spins (x,y,z) in the whole unit cell. Should require
spinorial wave functions [[nspinor]] = 2. Experimental: this does not work yet.
""",
),

Variable(
    abivarname="prtstm",
    varset="files",
    vartype="integer",
    topics=['printing_prgs', 'STM_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT the STM density",
    added_in_version="before_v9",
    text=r"""
If set to 1 or a larger value, provide output of the electron density in real
space rho(r), made only from the electrons close to the Fermi energy, in a
range of energy (positive or negative), determined by the (positive or
negative, but non-zero) value of the STM bias [[stmbias]].
This is a very approximate way to obtain STM profiles: one can choose an
equidensity surface, and consider that the STM tip will follow this surface.
Such equidensity surface might be determined with the help of Cut3D, and
further post-processing of it (to be implemented). The big approximations of
this technique are: neglect of the finite size of the tip, and position-
independent transfer matrix elements between the tip and the surface.
The charge density is provided in units of electrons/Bohr^3. The name of the
STM density file will be the root output name, followed by _STM. Like a _DEN
file, it can be analyzed by cut3d.
The file structure of this unformatted output file is described in [[help:abinit#denfile|this section]].
For the STM charge density to be generated, one must give, as an input file,
the converged wavefunctions obtained from a previous run, at exactly the same
k-points and cut-off energy, self-consistently determined, using the
occupation numbers from [[occopt]] = 7.
In the run with positive [[prtstm]], one has to use:

  * positive [[iscf]]
  * [[occopt]] = 7, with specification of [[tsmear]]
  * [[nstep]] = 1
  * the [[tolwfr]] convergence criterion
  * [[ionmov]] = 0 (this is the default value)
  * [[optdriver]] = 0 (this is the default value)

Note that you might have to adjust the value of [[nband]] as well, for the
treatment of unoccupied states, because the automatic determination of
[[nband]] will often not include enough unoccupied states.
When [[prtstm]] is non-zero, the stress tensor is set to zero.
No output of _STM file is provided by [[prtstm]] lower or equal to 0.
No other printing variables for density or potentials should be activated
(e.g. [[prtden]] has to be set to zero).
""",
),

Variable(
    abivarname="prtsuscep",
    varset="files",
    vartype="integer",
    topics=['printing_prngs'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT the SUSCEPtibility file (the irreducible polarizability)",
    added_in_version="before_v9",
    text=r"""
If set to 0, no _SUSC file will be produced after the screening calculation,
only the _SCR file will be output.
""",
),

Variable(
    abivarname="prtvclmb",
    varset="files",
    vartype="integer",
    topics=['printing_prpot'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT V CouLoMB",
    added_in_version="before_v9",
    text=r"""
If set >= 0 outputs a file with the Coulomb potential, defined as Hartree +
local Pseudopotential.

If **prtvclmb=1** and in case of PAW ([[usepaw]] > 0), the full core potential
is added for the Hartree part, with the on-site corrections vh1 - vht1.

If **prtvclmb=2**, only the smooth part of the Coulomb potential is output.
""",
),

Variable(
    abivarname="prtvdw",
    varset="vdw",
    vartype="integer",
    topics=['printing_prgs'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT Van Der Waals file",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Print out a NetCDF file containing a vdW-DF kernel.
""",
),

Variable(
    abivarname="prtvha",
    varset="files",
    vartype="integer",
    topics=['printing_prpot'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT V_HArtree",
    added_in_version="before_v9",
    text=r"""
If set >=1, provide output of the Hartree potential.

If [[ionmov]] == 0, the name of the potential file will be the root output name,
followed by _VHA.
If [[ionmov]] /= 0, potential files will be output at each time step, with
the name being made of

  * the root output name,
  * followed by _TIMx, where x is related to the time step (see later)
  * then followed by _VHA.

The file structure of this unformatted output file is described in [[help:abinit#localpotfile|this section]].
No output is provided by a negative value of this variable.
""",
),

Variable(
    abivarname="prtvhxc",
    varset="files",
    vartype="integer",
    topics=['printing_prpot'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT V_HXC",
    added_in_version="before_v9",
    text=r"""
If set >=1, provide output of the sum of the Hartree potential and xc potential.

If [[ionmov]] == 0, the name of the potential file will be the root output name,
followed by _VHXC.
If [[ionmov]] /= 0, potential files will be output at each time step, with
the name being made of

  * the root output name,
  * followed by _TIMx, where x is related to the time step (see later)
  * then followed by _VHXC.

The file structure of this unformatted output file is described in [[help:abinit#localpotfile|this section]].
No output is provided by a negative value of this variable.
""",
),

Variable(
    abivarname="prtvol",
    varset="files",
    vartype="integer",
    topics=['printing_prgs', 'Output_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT VOLume",
    added_in_version="before_v9",
    text=r"""
Control the volume of printed output. In particular, this concerns the
explicit echo of eigenenergies and residuals for all bands and k points in the
main output file. Also, the analysis of the value and location of the maximal
density (and magnetization).
Standard choice is 0. Positive values (all are allowed) generally print more and more in the output and log files,
while negative values are for debugging (or preprocessing only), and cause the
code to stop at some point.

  * 0 --> The eigenenergies and residuals for all bands and k points are not echoed in the main output file. There are exceptions: the eigenvalues of the first k point are printed at the end of the SCF loop, and also, if [[iscf]] = -2 and [[kptopt]]<=0, the eigenvalues for all the k points are printed anyway, for a maximum of 50 k-points. Due to some subtlety, if for **some** dataset [[prtvol]] is non-zero, the limit for input and output echoes cannot be enforced, so it is like if [[prtvol]] = 1 for **all** the datasets for which [[prtvol]] was set to 0.
  * 1 --> the eigenvalues for the first k-point are printed in all cases, at the end of the SCF loop.
  * 2 --> all the eigenvalues and the residuals are printed at the end of the SCF loop. Also, the analysis of the value and location of the maximal density (and magnetization) is printed.
  * 3 --> Print memory information for lobpcg.
  * 4 --> Like 3 and prints information of lobpcg algorithm convergence.
  * 10 --> the eigenvalues are printed for every SCF iteration, as well as other additions.
  * 11 --> even more information ...

Debugging options:

  * = -1 --> stop in abinit (main program), before call driver. Useful to see the effect of the preprocessing of input variables (memory needed, effect of symmetries, k points...) without going further. Run very fast, on the order of the second.
  * =-2 --> same as -1, except that print only the first dataset. All the non default input variables associated to all datasets are printed in the output file, but only for the first dataset. Also all the input variables are written in the NetCDF file "OUT.nc", even if the value is the default.
  * = -3 --> stop in gstate, before call scfcv, move or brdmin. Useful to debug pseudopotentials
  * = -4 --> stop in move, after completion of all loops
  * = -5 --> stop in brdmin, after completion of all loops
  * = -6 --> stop in scfcv, after completion of all loops
  * = -7 --> stop in vtorho, after the first rho is obtained
  * = -8 --> stop in vtowfk, after the first k point is treated
  * = -9 --> stop in cgwf, after the first wf is optimized
  * = -10 --> stop in getghc, after the Hamiltonian is applied once

This debugging feature is not yet activated in the RF routines. Note that
[[fftalg]] offers another option for debugging.
""",
),

Variable(
    abivarname="prtvolimg",
    varset="files",
    vartype="integer",
    topics=['printing_prgs', 'Output_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT VOLume for IMaGes",
    added_in_version="before_v9",
    text=r"""
Control the volume of printed output when an algorithm using images of the
cell is used ([[nimage]] > 1).
When such an algorithm is activated, the printing volume (in output file) can
be large and difficult to read.
Using **prtvolimg=1**, the printing volume, for each image, is reduced to
unit cell, atomic positions, total energy, forces, stresses, velocities and
convergence residuals.
Using **prtvolimg=2**, the printing volume, for each image, is reduced to
total energy and convergence residuals only.
""",
),

Variable(
    abivarname="prtvpsp",
    varset="files",
    vartype="integer",
    topics=['printing_prpot'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT V_PSeudoPotential",
    added_in_version="before_v9",
    text=r"""
If set >=1, provide output of the local pseudo potential.

If [[ionmov]] == 0, the name of the potential file will be the root output name, followed by _VPSP.
If [[ionmov]] /= 0, potential files will be output at each time step, with the name being made of

  * the root output name,
  * followed by _TIMx, where x is related to the timestep (see later)
  * then followed by _VPSP.

The file structure of this unformatted output file is described in [[help:abinit#localpotfile|this section]].
No output is provided by a negative value of this variable.
""",
),

Variable(
    abivarname="prtvxc",
    varset="files",
    vartype="integer",
    topics=['printing_prpot'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT V_XC",
    added_in_version="before_v9",
    text=r"""
If set >=1, provide output of the exchange-correlation potential.

If [[ionmov]] == 0, the name of the potential file will be the root output name,
followed by _VXC.
If [[ionmov]] /= 0, potential files will be output at each time step, with
the name being made of

  * the root output name,
  * followed by _TIMx, where x is related to the timestep (see later)
  * then followed by _VXC.

The file structure of this unformatted output file is described in [[help:abinit#localpotfile|this section]].
No output is provided by a negative value of this variable.
""",
),

Variable(
    abivarname="prtwant",
    varset="files",
    vartype="integer",
    topics=['printing_prgs', 'Wannier_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT WANT file",
    added_in_version="before_v9",
    text=r"""
Flag used to indicate that either the Wannier90 or the WanT interfaces will be used.

  * [[prtwant]] = 1 --> Use the **ABINIT- WanT** interface.

Provide an output file that can be used by the WanT postprocessing program
(see [http://www.wannier-transport.org](http://www.wannier-transport.org)). The value of the prtwant indicates the
version of the WanT code that can read it. Currently only the value
[[prtwant]] = 1 is implemented, corresponding to WanT version 1.0.1, available
since Oct. 22, 2004.

Notes: Several requirements must be fulfilled by the wavefunction. Among
them, two are mandatory:

* An uniform grid of k-points, including the GAMMA point must be used.
* The use of time reversal symmetry is not allowed (istwfk=1)
* The list of k-points must be ordered, such that the coordinates, namely
  three-components vectors has the third index varying the most rapidly, then the second index, then the first index

If these requirement are not fulfilled, the program will stop and an error message is returned.

As an example of k-point grid in case of systems that have some 3D character
(1D systems are easy):

    nkpt 8
    kpt
    0   0   0
    0   0   1/2
    0   1/2 0
    0   1/2 1/2
    1/2 0   0
    1/2 0   1/2
    1/2 1/2 0
    1/2 1/2 1/2
    istwfk *1

Also, in order to use WanT as a post-processing program for ABINIT you might
have to recompile it with the appropriate flags (see ABINIT makefile). Up to
now only the -convert big-endian was found to be mandatory, for machines with
little-endian default choice.

  * [[prtwant]] = 2 --> Use the **ABINIT- Wannier90** interface.

ABINIT will produce the input files required by Wannier90 and it will run
Wannier90 to produce the Maximally-locallized Wannier functions (see [
http://www.wannier.org ](http://www.wannier.org) ).

!!! Notes

    * The files that are created can also be used by Wannier90 in stand-alone mode.
    * In order to use Wannier90 as a post-processing program for ABINIT you might have to recompile it with the appropriate flags (see ABINIT makefile). You might use ./configure --enable-wannier90
    * There are some other variables related to the interface of Wannier90 and ABINIT. See [[varset:w90]].

  * [[prtwant]] = 3 --> Use the **ABINIT- Wannier90** interface after converting the input wavefunctions to **quasi-particle** wavefunctions.

ABINIT will produce the input files required by Wannier90 and it will run
Wannier90 to produce the Maximally-localized Wannier functions (see [
http://www.wannier.org ](http://www.wannier.org) ).

!!! Notes

    * An input file of LDA wave functions is required which is completely consistent with the _KSS file used in the self-consistent GW calculation. This means that [[kssform]] 3 must be used to create the _KSS file and the output _WFK file from the same run must be used as input here.
    * Wannier90 requires [[nshiftk]] = 1, and [[shiftk]] =  0 0 0 is recommended. The k-point set used for the GW calculation, typically the irreducible BZ set created using [[kptopt]] = 1, and that for the Abinit- Wannier90 interface must be consistent.
    * Full-BZ wavefunctions should be generated in the run calling the interface by setting [[kptopt]] = 3, [[iscf]] = -2, and [[nstep]] = 3. This will simply use symmetry to transform the input IBZ wavefunctions to the full BZ set, still consistent with the GW _KSS input.
    * The final _QPS file created by the self-consistent GW run is required as input.
    * Any value of [[gwcalctyp]] between between 20 and 29 should be suitable, so, for example, Hartree-Fock maximally-localized Wannier functions could be generated setting [[gwcalctyp]] = 25.
""",
),

Variable(
    abivarname="prtwf",
    varset="files",
    vartype="integer",
    topics=['printing_prden', 'vdw_useful'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[nimage]] > 1': 0, 'defaultval': 1}),
    mnemonics="PRinT the WaveFunction",
    added_in_version="before_v9",
    text=r"""
If [[prtwf]] = 1, provide output of wavefunction and eigenvalue file
The file structure of this unformatted output file is described in [[help:abinit#wfkfile|this section]].
For a standard ground-state calculation, the name of the wavefunction file
will be the root output name, followed by _WFK. If [[nqpt]] = 1, the root name
will be followed by _WFQ. For response-function calculations, the root name
will be followed by _1WFx, where x is the number of the perturbation. The
dataset information will be added as well, if relevant.
No wavefunction output is provided by [[prtwf]] = 0.
If [[prtwf]] = -1, the code writes the wavefunction file only if convergence is
not achieved in the self-consistent cycle.


If [[prtwf]] = 2, a file pwfn.data is produced, to be used as input for the
CASINO QMC code. See more explanation at the end of this section.
If [[prtwf]] = 3, the file that is created is nearly the same as with
[[prtwf]] = 1, except that the records that should contain the wavefunction is
empty (so, such records exist, but store nothing). This is useful to generate
size-reduced DDK files, to perform an optic run. Indeed, in the latter case,
only matrix elements are needed [so, no wavefunction], but possibly a large
number of conduction bands, so that the DDK file might be huge if it contains
the wavefunctions.

Further explanation for the [[prtwf]] = 2 case. To produce a wave function
suitable for use as a CASINO trial wave function, certain ABINIT parameters
must be set correctly. Primarily, CASINO (and QMC methods generally) can only
take advantage of time-reversal symmetry, and not the full set of symmetries
of the crystal structure. Therefore, ABINIT must be instructed to generate
k-points not just in the Irreducible Brillouin Zone, but in a full half of the
Brillouin Zone (using time-reversal symmetry to generate the other half).
Additionally, unless instructed otherwise, Abinit avoids the need for internal
storage of many of the coefficients of its wave functions for k-points that
have the property 2k=G_latt, where G_latt is a reciprocal lattice vector, by
making use of the property that c_k(G)=c^*_k(-G-G_latt). Abinit must be
instructed not to do this in order to output the full set of coefficients for
use in CASINO. See the ABINIT theoretical background documents
ABINIT/Infos/Theory/geometry.pdf and ABINIT/Infos/Theory/1WF.pdf for more
information.
The first of these requirements is met by setting the ABINIT input variable
[[kptopt]] to 2 and the second by setting
[[istwfk]] to 1 for all the k points. Since
CASINO is typically run with relatively small numbers of k-points, this is
easily done by defining an array of "1" in the input file.
For example, for the 8 k-points generated with ngkpt 2 2 2, we add the
following lines to the input file:

    # Turn off special storage mode for time-reversal k-points
    istwfk 1 1 1 1 1 1 1 1
    # Use only time reversal symmetry, not full set of symmetries.
    kptopt 2


Other useful input variables of relevance to the plane waves ABINIT will
produce include ecut, nshiftk, shiftk, nband, occopt, occ, spinat and nsppol
(see relevant input variable documents in ABINIT/Infos/). If ABINIT is run in
multiple dataset mode, the different wave functions for the various datasets
are exported as pwfn1.data, pwfn2.data, ..., pwfnn.data where the numbers are
the contents of the contents of the input array jdtset (defaults to
1,2,...,ndtset).
Once the routine is incorporated into the ABINIT package it is anticipated
that there will be an input variable to control whether or not a CASINO
pwfn.data file is written.

Other issues related to [[prtwf]] = 2.
The exporter does not currently work when ABINIT is used in parallel mode on
multiple processors if k-point parallelism is chosen. ABINIT does not store
the full wave function on each processor but rather splits the k-points
between the processors, so no one processor could write out the whole file.
Clearly this could be fixed but we have not done it yet. The sort of plane
wave DFT calculations usually required to generate QMC trial wave functions
execute very rapidly anyway and will generally not require a parallel
machines. The outqmc routine currently bails out with an error if this
combination of modes is selected - this will hopefully be fixed later.
There has not been very extensive testing of less common situations such as
different numbers of bands for different k-points, and more complicated spin
polarized systems, so care should be taken when using the output in these circumstances.
If there is any doubt about the output of this routine, the first place to
look is the log file produced by ABINIT: if there are any warnings about
incorrectly normalized orbitals or non-integer occupation numbers there is
probably something set wrong in the input file.
""",
),

Variable(
    abivarname="prtwf_full",
    varset="files",
    vartype="integer",
    topics=['printing_prden'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT Wavefunction file on the FULL mesh",
    requires="[[prtwf]] == 1",
    added_in_version="before_v9",
    text=r"""
If set to 1 in a ground-state calculation, the code will output another WFK
file (with extension FULL_WFK) containing the wavefunctions in the full BZ as
well as a text file with the tables used for the tetrahedron method. Note that
prtwf_full requires [[prtwf]] == 1 and a ground-state calculation done on a
homogeneous k-mesh (see [[ngkpt]] and [[shiftk]]). The tetrahedron table is
produced only if the number of k-points in the irreducible zone ([[nkpt]]) is
greater than 3.
""",
),

Variable(
    abivarname="prtxml",
    varset="files",
    vartype="integer",
    topics=['printing_prgs'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT an XML output",
    added_in_version="before_v9",
    text=r"""
Create an XML output with common values. The corresponding DTD is distributed
in sources as extras/post_processing/abinitRun.dtd. All the DTD is not yet
implemented and this one is currently restricted to ground-state computations
(and derivative such as geometry optimisation).
""",
),

Variable(
    abivarname="ptcharge",
    varset="paw",
    vartype="real",
    topics=['EFG_basic'],
    dimensions=['[[ntypat]]'],
    defaultval=MultipleValue(number=None, value=0),
    mnemonics="PoinT CHARGEs",
    requires="[[usepaw]] == 1 and [[prtefg]]>=3",
    added_in_version="before_v9",
    text=r"""
  * Array of point charges, in atomic units, of the nuclei. In the normal computation of electric field gradients (see [[prtefg]]) the ionic contribution is calculated from the core charges of the atomic sites. Thus for example in a PAW data set for oxygen where the core is $1s^{2}$, the core charge is +6 (total nuclear charge minus core electron charge). In point charge models, which are much less accurate than PAW calculations, all atomic sites are treated as ions with charges determined by their valence states. In such a case oxygen almost always would have a point charge of -2. The present variable taken together with [[prtefg]] performs a full PAW computation of the electric field gradient and also a simple point charge computation. The user inputs whatever point charges he/she wishes for each atom type.
""",
),

Variable(
    abivarname="ptgroupma",
    varset="geo",
    vartype="integer",
    topics=['spinpolarisation_internal', 'SmartSymm_internal'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PoinT GROUP number for the MAgnetic space group",
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
This internal variable characterizes a Shubnikov type III magnetic space group
(anti-ferromagnetic space group). The user is advised to consult [[cite:Bradley1972]].
In a Shubnikov type III magnetic space group, the primitive cell is the same if one takes
into account the spin-flipping operations or if one does not take into account such
spin-flipping operations. Explicitely, there is no pure translation with a spin-flip.

A Shubnikov type III magnetic space group might be defined by its Fedorov
space group (set of all spatial symmetries, irrespective of their magnetic
action), and the halving space group (only the symmetries that do not change
the magnetization).
The specification of the halving space group might be done by specifying, for
each point symmetry, the magnetic action. See Table 7.1 of the above-mentioned
reference. Magnetic point groups are numbered from 1 to 58.
The halving space group is the group of the symmetry operations for which [[symafm]]=1,
namely those operations that are not accompanied with a spin flip..
Note that the definition of a spin flip is different for the [[nspden]]=2 and the [[nspden]]=4 cases,
see the description of [[symafm]].

Related input variables: [[spgroup]], [[spgroupma]], [[genafm]], [[symafm]].
""",
),

Variable(
    abivarname="pvelmax",
    varset="gw",
    vartype="real",
    topics=['RandStopPow_basic'],
    dimensions=[3],
    defaultval=MultipleValue(number=3, value=1.0),
    mnemonics="Particle VELocity MAXimum",
    requires="[[optdriver]] == 3",
    added_in_version="before_v9",
    text=r"""
When [[npvel]] is larger than 0, it performs electronic stopping power
calculations on a velocity grid along the direction determined by [[pvelmax]].
The vector [[pvelmax]] defines both the direction and the maximum velocity.
[[pvelmax]] is input in Cartesian coordinates.
""",
),

Variable(
    abivarname="pw_unbal_thresh",
    varset="paral",
    vartype="real",
    topics=['parallelism_expert'],
    dimensions="scalar",
    defaultval="40%",
    mnemonics="Plane Wave UNBALancing: THRESHold for balancing procedure",
    requires="[[paral_kgb]] == 1",
    added_in_version="before_v9",
    text=r"""
This parameter (in %) activates a load balancing procedure when the
distribution of plane wave components over MPI processes is not optimal. The
balancing procedure is activated when the ratio between the number of plane
waves treated by a processor and the ideal one is higher than
_pw_unbal_thresh_ %.
""",
),

Variable(
    abivarname="qmass",
    varset="rlx",
    vartype="real",
    topics=['PIMD_basic', 'MolecularDynamics_basic'],
    dimensions=['[[nnos]]'],
    defaultval=MultipleValue(number=None, value=10.0),
    mnemonics="Q thermostat MASS",
    added_in_version="before_v9",
    text=r"""
This are the masses of the chains of [[nnos]] thermostats to be used when
[[ionmov]] = 13 (Molecular Dynamics) or [[imgmov]] = 13 (Path Integral Molecular
Dynamics).

If [[ionmov]] = 13 (Molecular Dynamics), this temperature control can be used
with  [[optcell]] =0, 1 (homogeneous cell deformation) or 2 (full cell deformation).
If [[imgmov]] = 13 (Path Integral Molecular Dynamics), this temperature control
can be used with  [[optcell]] =0 (NVT ensemble) or 2 (fully flexible NPT
ensemble). In that case, [[optcell]] = 2 is NOT USABLE yet.
""",
),

Variable(
    abivarname="qprtrb",
    varset="ffield",
    vartype="integer",
    topics=['Artificial_useful'],
    dimensions=[3],
    defaultval=[0, 0, 0],
    mnemonics="Q-wavevector of the PERTurbation",
    characteristics=['[[DEVELOP]]'],
    requires="[[vprtrb]]",
    added_in_version="before_v9",
    text=r"""
Gives the wavevector, in units of reciprocal lattice primitive translations,
of a perturbing potential of strength [[vprtrb]]. See [[vprtrb]] for more info.
""",
),

Variable(
    abivarname="qpt",
    varset="gstate",
    vartype="real",
    topics=['q-points_useful'],
    dimensions=[3],
    defaultval=[0, 0, 0],
    mnemonics="Q PoinT",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Only used if [[nqpt]] = 1.

Combined with [[qptnrm]], define the q vector [[qptn]](1:3) in the case
[[qptopt]] = 0.

This input variable is not internal ([[qptn]](1:3) is used instead), but is
used to echo the value of [[qptn]](1:3), with renormalisation factor one.
""",
),

Variable(
    abivarname="qptdm",
    varset="gw",
    vartype="real",
    topics=['Susceptibility_useful'],
    dimensions=[3, '[[nqptdm]]'],
    defaultval=MultipleValue(number=None, value=0.0),
    mnemonics="Q-PoinTs for the Dielectric Matrix",
    requires="[[optdriver]] == 3 and [[nqptdm]]!=0",
    added_in_version="before_v9",
    text=r"""
[[qptdm]] contains the set of q points used in the screening part of ABINIT,
instead of the automatic generation of the q points when [[nqptdm]] = 0. These q
points are given in terms of reciprocal space primitive translations (**not** in
cartesian coordinates!). For further explanation, see the input variable [[nqptdm]].
""",
),

Variable(
    abivarname="qptn",
    varset="internal",
    vartype="real",
    topics=[' DFPT_internal'],
    dimensions=[3],
    defaultval=MultipleValue(number=3, value=0),
    mnemonics="Q-PoinT re-Normalized",
    characteristics=['[[INTERNAL_ONLY]]'],
    requires="[[nqpt]] == 1",
    added_in_version="before_v9",
    text=r"""
Only used if [[nqpt]] = 1.

In ground-state calculation, the vector **qptn**(1:3) is added to each
renormalized k point (whatever the value of [[kptopt]] used) to
generate the normalized, shifted, set of k-points [[kptns]](1:3,1: **nkpt** ).
In response-function calculations, **qptn**(1:3) defines the wavevector of the
phonon-type calculation.

**qptn**(1:3) can be produced on the basis of the different methods described
in [[qptopt]], like using [[qpt]](1:3) with renormalisation provided by
[[qptnrm]], or using the other possibilities defined by [[iqpt]], [[ngqpt]],
[[nshiftq]], [[qptrlatt]], [[shiftq]].

For insulators, there is no restriction on the q-points to be used for the
perturbations. By contrast, for metals, for the time being, it is advised to
take q points for which the k and k+q grids are the same (when the periodicity
in reciprocal space is taken into account). Tests remain to be done to see
whether other q points might be allowed (perhaps with some modification of the code).
""",
),

Variable(
    abivarname="qptnrm",
    varset="gstate",
    vartype="real",
    topics=['q-points_useful'],
    dimensions="scalar",
    defaultval=1.0,
    mnemonics="Q PoinTs NoRMalization",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Only used if [[nqpt]] = 1 and [[qptopt]] = 0

Provides re-normalization of [[qpt]]. Must be positive, non-zero. The actual q
vector (renormalized) is [[qptn]](1:3)= [[qpt]](1:3)/[[qptnrm]].
""",
),

Variable(
    abivarname="qptopt",
    varset="gstate",
    vartype="integer",
    topics=['q-points_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="QPoinTs OPTion",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Only used if [[nqpt]] = 1.

Controls the set up to generate the Q point [[qptn]](1:3) to be used for the
specific dataset, either as a shift of k-point grid in ground-state
calculations, or as a stand-alone phonon wavevector.

There are two basic techniques to generate the Q point: either by specifying
it directly, possibly with a renormalisation factor ([[qptopt]] = 0), or
extracting it from a grid a Q points ([[qptopt]] = 1 to 4), using the index
[[iqpt]]. At variance with the similar generation of k points, only ONE q
point can be used per dataset.

With [[qptopt]] = 1 to 4, rely on [[ngqpt]] or [[qptrlatt]], as well as on
[[nshiftq]] and [[shiftq]] to set up a q point grid, from which the q point
with number [[iqpt]] will be selected. The values [[qptopt]] = 1 to 4 differ by
the treatment of symmetries. Note that the symmetries are recomputed starting
from the values of [[rprimd]], [[xred]] and [[spinat]]. So, the explicit value
of [[symrel]] are not used. This is to allow doing calculations with
[[nsym]] = 1, sometimes needed for T-dependent electronic structure, still
decreasing the number of q points in the case [[qptopt]] = 1 or [[qptopt]] = 3.

  * 0 --> read directly [[qpt]], and its (eventual) renormalisation factor [[qptnrm]].

  * 1 --> Take fully into account the symmetry to generate the grid of q points in the Irreducible Brillouin Zone only.
    (This is the usual mode for RF calculations)

  * 2 --> Take into account only the time-reversal symmetry: q points will be generated in half the Brillouin zone.

  * 3 --> Do not take into account any symmetry: q points will be generated in the full Brillouin zone.

  * 4 --> Take into account all the symmetries EXCEPT the time-reversal symmetry to generate
    the k points in the Irreducible Brillouin Zone.

In the case of a grid of q points, the auxiliary variables [[kptrlen]],
[[ngkpt]] and [[prtkpt]] might help you to select the optimal grid, similarly
to the case of the K point grid.
""",
),

Variable(
    abivarname="qptrlatt",
    varset="gstate",
    vartype="integer",
    topics=['q-points_useful'],
    dimensions=[3, 3],
    defaultval=MultipleValue(number=None, value=0),
    mnemonics="Q - PoinTs grid: Real space LATTice",
    characteristics=['[[INPUT_ONLY]]'],
    excludes="[[ngqpt]]",
    added_in_version="before_v9",
    text=r"""
This input variable is used only when [[qptopt]] is positive. It partially
defines the q point grid. The other piece of information is contained in
[[shiftq]]. [[qptrlatt]] cannot be used together with [[ngqpt]].

The values [[qptrlatt]](1:3,1), [[qptrlatt]](1:3,2), [[qptrlatt]](1:3,3) are
the coordinates of three vectors in real space, expressed in the [[rprimd]]
coordinate system (reduced coordinates). They defines a super-lattice in real
space. The k point lattice is the reciprocal of this super-lattice, possibly
shifted (see [[shiftq]]).

If neither [[ngqpt]] nor [[qptrlatt]] are defined, ABINIT will automatically
generate a set of k point grids, and select the best combination of
[[qptrlatt]] and [[shiftq]] that allows one to reach a sufficient value of
[[kptrlen]]. See this latter variable for a complete description of this procedure.
""",
),

Variable(
    abivarname="quadmom",
    varset="paw",
    vartype="real",
    topics=['EFG_basic'],
    dimensions=['[[ntypat]]'],
    defaultval=MultipleValue(number=None, value=0),
    mnemonics="QUADrupole MOMents",
    requires="[[usepaw]] == 1 and [[prtefg]]>=1",
    added_in_version="before_v9",
    text=r"""
  * Array of quadrupole moments, in barns, of the nuclei. These values are used in conjunction with the electric field gradients computed with [[prtefg]] to calculate the quadrupole couplings in MHz, as well as the asymmetries. Note that the electric field gradient at a nuclear site is independent of the nuclear quadrupole moment, thus the quadrupole moment of a nucleus can be input as 0, and the option [[prtefg]] = 2 used to determine the electric field gradient at the site.
""",
),

Variable(
    abivarname="random_atpos",
    varset="rlx",
    vartype="integer",
    topics=['crystal_expert', 'GeoOpt_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="RANDOM ATomic POSitions",
    added_in_version="before_v9",
    text=r"""
Control the inner coordinates, which can be generated randomly by using 4
different methods depending ont its value
(0) if zero, no random generation and xred are taken as they have been introduced by the user
(1) if one, particles are generated completely random within the unit cell.
(2) if two, particles are generated randomly but the inner particle distance
is always larger than a factor of the sum of the covalent bonds between the
atoms (note: this is incompatible with the definition of alchemical mixing,
in which [[ntypat]] differs from [[npsp]])
""",
),

Variable(
    abivarname="ratsm",
    varset="gstate",
    vartype="real",
    topics=['printing_prdos', 'MagMom_useful', 'ConstrainedDFT_useful'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'any([[constraint_kind]] > 1)': 0.05, 'defaultval': 0.00}),
    mnemonics="Radii of the ATomic spheres SMearing",
    added_in_version="before_v9",
    text=r"""
Smearing width for the atomic spheres whose radius is determined by [[ratsph]].
For each spherical zone around each atom, the integrating function goes
from 1.0 to 0.0 in an interval from [[ratsph]]-[[ratsm]] to [[ratsph]].
The function is the same as the one used to smear the kinetic energy, see [[ecutsm]].
""",
),


Variable(
    abivarname="ratsph",
    varset="gstate",
    vartype="real",
    topics=['printing_prdos', 'MagMom_useful', 'ElecBandStructure_useful', 'ElecDOS_useful', 'ConstrainedDFT_basic'],
    dimensions=['[[ntypat]]'],
    defaultval=ValueWithConditions({'[[usepaw]] == 1': '[[AUTO_FROM_PSP]]', 'defaultval': 2.00}),
    mnemonics="Radii of the ATomic SPHere(s)",
    added_in_version="before_v9",
    text=r"""
Relevant only when [[prtdensph]] = 1, or [[magconon]]/=0, or any([[constraint_kind]](:)/=0) (that is, constrained DFT), or [[prtdos]] = 3.
In most cases (see later for [[prtdos]] = 3), provides the radius of the spheres around each atom in which the total
charge density or magnetization will be integrated.
The integral within the sphere is obtained by a sum over real space FFT points
inside the sphere, multiplied by a function that is one inside the sphere, except in a small boundary zone determined by [[ratsm]],
where this fonction goes smoothly from 1 to 0.
In case of PAW, [[ratsph]] radius has to be greater or equal to PAW radius of
considered atom type (which is read from the PAW dataset file; see **rc_sph** or **r_paw**).
In case of constrained DFT, note that the sphere for different atoms are not allowed to overlap.

When [[prtdos]] = 3:

Provides the radius of the spheres around the [[natsph]] atoms of indices
[[iatsph]], in which the local DOS and its angular-momentum projections will
be analysed.

The choice of this radius is quite arbitrary. In a plane-wave
basis set, there is no natural definition of an atomic sphere. However, it
might be wise to use the following well-defined and physically motivated procedure:
from the Bader analysis, one can define the radius of the sphere that contains
the same charge as the Bader volume. This "Equivalent Bader charge atomic
radius" might then be used to perform the present analysis. See the
[[help:aim]] for more explanations. Another physically motivated choice would
be to rely on another charge partitioning, like the Hirshfeld one (see the
cut3d utility [[help:cut3d]]). The advantage of using charge partitioning schemes comes from
the fact that the sum of atomic DOS, for all angular momenta and atoms,
integrated on the energy range of the occupied states, gives back the total
charge. If this is not an issue, one could rely on the half of the nearest-
neighbour distances, or any scheme that allows one to define an atomic radius.
Note that the choice of this radius is however critical for the balance
between the s, p and d components. Indeed, the integrated charge within a
given radius, behave as a different power of the radius, for the different
channels s, p, d. At the limit of very small radii, the s component dominates
the charge contained in the sphere.
""",
),

Variable(
    abivarname="ratsph_extra",
    varset="gstate",
    vartype="real",
    topics=['printing_prdos'],
    dimensions="scalar",
    defaultval=ValueWithUnit(units='Bohr', value=2.0),
    mnemonics="Radii of the ATomic SPHere(s) in the EXTRA set",
    characteristics=['[[LENGTH]]'],
    added_in_version="before_v9",
    text=r"""
Radius for extra spheres the DOS is projected into. See [[natsph_extra]] and
[[xredsph_extra]] for the number and positions of the spheres.
""",
),

Variable(
    abivarname="rcut",
    varset="gw",
    vartype="real",
    topics=['GWls_compulsory', 'Susceptibility_basic', 'SelfEnergy_basic'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="Radius of the CUT-off for coulomb interaction",
    added_in_version="before_v9",
    text=r"""
Truncation of the Coulomb interaction in real space. The meaning of [[rcut]]
is governed by the cutoff shape option [[icutcoul]].

If [[rcut]] is negative, the cutoff is automatically calculated so to enclose
the same volume inside the cutoff as the volume of the primitive cell.
""",
),

Variable(
    abivarname="recefermi",
    varset="dev",
    vartype="real",
    topics=['Recursion_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="RECursion - initial guess  of the FERMI Energy",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Used in Recursion method ([[tfkinfunc]] = 2). In the first SCF calculation it
fixes the initial guess for the Fermi energy.
""",
),

Variable(
    abivarname="recgratio",
    varset="dev",
    vartype="integer",
    topics=['Recursion_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="RECursion - Grid RATIO",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Used in Recursion method ([[tfkinfunc]] = 2). It represents the ratio of the two
grid step: [[recgratio]] = fine_step/coarse_step and it is bigger or equal than
1. It introduces a double-grid system which permits to compute the electronic
density on a coarse grid, using a fine grid (defined by [[ngfft]]) in the
discretisation of the green kernel (see [[recptrott]]). Successively the
density and the recursion coefficients are interpolated on the fine grid by
FFT interpolation. Note that ngfft/recgratio=number of points of the coarse
grid has to be compatible with the parallelization parameters.
""",
),

Variable(
    abivarname="recnpath",
    varset="dev",
    vartype="integer",
    topics=['Recursion_expert'],
    dimensions="scalar",
    defaultval=500,
    mnemonics="RECursion - Number of point for PATH integral calculations",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Used in Recursion method ([[tfkinfunc]] = 2). Determine the number of
discretisation points to compute some path integral in the recursion method;
those path integrals are used to compute the entropy and the eigenvalues
energy. during the latest SFC cycles.
""",
),

Variable(
    abivarname="recnrec",
    varset="dev",
    vartype="integer",
    topics=['Recursion_expert'],
    dimensions="scalar",
    defaultval=10,
    mnemonics="RECursion - Number of RECursions",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Used in Recursion method ([[tfkinfunc]] = 2). Determine the maximum order of
recursion, that is the dimension of the krylov space we use to compute
density. If the precision set by [[rectolden]] is reached before that order,
the recursion method automatically stops.
""",
),

Variable(
    abivarname="recptrott",
    varset="dev",
    vartype="integer",
    topics=['Recursion_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="RECursion - TROTTer parameter",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Used in Recursion method ([[tfkinfunc]] = 2). Determine the trotter parameter
used to compute the exponential of the hamiltonian in the recursion method:

$$ e^{-\beta(-\Delta + V)} \approx
\left(
    e^{-\frac{\beta}{4c} V}
    e^{-\frac{\beta}{4c} \Delta}
    e^{-\frac{\beta}{4c} V}
\right)^{2c} $$

where $c$=[[recptrott]].
If set to 0, we use [[recptrott]] = $1/2$ in the above formula.
Increasing [[recptrott]] improve the
accuracy of the trotter formula, but increase the dicretisation error: it may
be necessary to increase [[ngfft]]. The discretisation error is essentially
the discretisation error of the green kernel $e^{\frac{c}{\beta|r|^2}}$ on
the ngfft grid.
""",
),

Variable(
    abivarname="recrcut",
    varset="dev",
    vartype="integer",
    topics=['Recursion_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="RECursion - CUTing Radius",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Used in Recursion method ([[tfkinfunc]] = 2). Used to improve the computational
time in the case of the recursion method in a large cell: the density at a
point will be computed with taking account only of a sphere of radius [[recrcut]].
""",
),

Variable(
    abivarname="rectesteg",
    varset="dev",
    vartype="integer",
    topics=['Recursion_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="RECursion - TEST on Electron Gas",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Used in Recursion method ([[tfkinfunc]] = 2). It is used to test an electron gas
by putting the ion potential equal to zero.
""",
),

Variable(
    abivarname="rectolden",
    varset="dev",
    vartype="real",
    topics=['Recursion_basic'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="RECursion - TOLerance on the difference of electronic DENsity",
    characteristics=['[[DEVELOP]]'],
    commentdefault="Default value to be changed.",
    added_in_version="before_v9",
    text=r"""
Used in Recursion method ([[tfkinfunc]] = 2). Sets a tolerance for differences
of electronic density that, reached TWICE successively, will cause one SCF
cycle to stop. That electronic density difference is computed in the infinity
norm (that is, it is computed point-by-point, and then the maximum difference is computed).
""",
),

Variable(
    abivarname="red_dfield",
    varset="ffield",
    vartype="real",
    topics=['Berry_useful'],
    dimensions=[3],
    defaultval=MultipleValue(number=3, value=0.0),
    mnemonics="REDuced Displacement FIELD",
    requires="[[berryopt]] = 16, [[red_efield]]",
    added_in_version="before_v9",
    text=r"""
In case [[berryopt]] = 16, a reduced finite electric displacement field
calculation is performed. The value of this displacement field, and its
direction is determined by [[red_dfield]]. It must be given in atomic units.

[[red_dfield]] is defined via Eq.(26) in the Supplement of [[cite:Stengel2009]].
""",
),

Variable(
    abivarname="red_efield",
    varset="ffield",
    vartype="real",
    topics=['Berry_useful'],
    dimensions=[3],
    defaultval=MultipleValue(number=3, value=0.0),
    mnemonics="REDuced Electric FIELD",
    requires="[[berryopt]] = 16",
    added_in_version="before_v9",
    text=r"""
In case [[berryopt]] = 16, a reduced finite electric displacement field
calculation is performed. In this case, the parameter [[red_efield]] specifies the
initial electric field used on the first iteration, in atomic units.

[[red_efield]] is defined via Eq.(25) in the Supplement [[cite:Stengel2009]].
""",
),

Variable(
    abivarname="red_efieldbar",
    varset="ffield",
    vartype="real",
    topics=['Berry_useful'],
    dimensions=[3],
    defaultval=MultipleValue(number=3, value=0.0),
    mnemonics="REDuced Electric FIELD BAR",
    requires="[[berryopt]] = 14",
    added_in_version="before_v9",
    text=r"""
In case [[berryopt]] = 14, a reduced finite electric field calculation is
performed. The magnitude and direction of this electric field are determined
by [[red_efieldbar]]. It must be given in atomic units.

[[red_efieldbar]] is defined via Eq.(28) in the Supplement of [[cite:Stengel2009]].
""",
),

Variable(
    abivarname="restartxf",
    varset="rlx",
    vartype="integer",
    topics=['PIMD_useful', 'MolecularDynamics_useful', 'GeoOpt_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="RESTART from (X,F) history",
    added_in_version="before_v9",
    text=r"""
Control the restart of a molecular dynamics or structural optimization job.

**restartxf > 0 (Deprecated) **: The code reads from the input WFK file, the
previous history of atomic coordinates and corresponding forces, in order to
continue the work done by the job that produced this wf file. If
[[optcell]] /= 0, the history of [[acell]] and [[rprim]] variables is also taken
into account. The code will take into consideration the whole history (if
**restartxf** = 1), or discard the few first (x,f) pairs, and begin only at the
pair whose number corresponds to **restartxf**.
Works only for [[ionmov]] = 2 or 22 (Broyden) and when an input wavefunction file is
specified, thanks to the appropriate values of [[irdwfk]] or [[getwfk]].

NOTES:
* The input WFK file must have been produced by a run that exited cleanly.
  It cannot be one of the temporary wf files that exist when a job crashed.
* One cannot restart a calculation with a non-zero [[optcell]] value from the (x,f) history of another run with a different non-zero [[optcell]] value. Starting a non-zero [[optcell]] run from a zero [[optcell]] run should work.
* Deprecated, the use of the new options (-1 and -2) is preferred.

**restartxf = 0 (Default)**: No restart procedure is enabled and the code will start a
Molecular dynamics or structural optimization from scratch.

**restartxf = -1 (New)**: Use the HIST.nc file to reconstruct a partial
calculation. It will reconstruct the different configurations using the forces
and the stresses stored in the HIST.nc file, instead of calling the SCF procedure.
Using **restartxf = -1** from the beginning is harmless. The only condition is
to keep the input file the same in such a way that the same predictor is used
and it will predict the same structure recorded in the HIST.nc file.
This option will always compute extra [[ntime]] iterations independent of the
number of iterations recovered previously.

**restartxf = -2 (New)**: Read the HIST.nc file and select the atomic positions and
cell parameters with the lowest energy. Forget all the history and start the
calculation using those values. The original atomic coordinates and cell
parameters are irrelevant in this case.

**restartxf = -3 (New)**: Read **ONLY** the last required atomic positions and cell parameters
in the HIST.nc file to restart the Molecular dynamics or structural optimization.

NOTES:

* You can use **restartxf=-1, -2 or -3** for all predictors that make no use of random numbers.
* You can use **restartxf=-1, -2 or -3** to restart a calculation that was not completed. The HIST.nc file is written on each iteration. So you always have something to recover from.
* You can take advantage of the appropriate values of [[irdwfk]] or [[getwfk]] to get a good wave function to continue your job.
""",
),

Variable(
    abivarname="rf2_dkdk",
    varset="dfpt",
    vartype="integer",
    topics=['DFPT_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Response Function: 2nd Derivative of wavefunctions with respect to K",
    added_in_version="before_v9",
    text=r"""
If is equal to 1, activates computation of second derivatives of wavefunctions with respect to
wavevectors (ipert = natom+10 is activated). This is not strictly a response function but is a needed
auxiliary quantity in the calculations of 3rd-order derivatives of the energy
(non-linear response) if [[usepead]] == 0. The directions for the derivatives are determined by
[[rf2_pert1_dir]] and [[rf2_pert2_dir]] and [[prepanl]] as the following:

The computation of the 2nd derivative of wavefunction with respect to "lambda_1" and "lambda_2" is computed if
if rf2_pert1_dir[idir1] AND rf2_pert2_dir[idir2] are equal to 1, where "idir1" ("idir2") is direction of
the perturbation "lambda_1" ("lambda_2").
If ALL directions are activated (default behavior) AND [[prepanl]] == 1, then the code automatically selects
only the directions that will be used by the non-linear routine ([[optdriver]] == 5) using crystal symmetries.
""",
),

Variable(
    abivarname="rf2_dkde",
    varset="dfpt",
    vartype="integer",
    topics=['DFPT_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Response Function: mixed 2nd Derivative of wavefunctions with respect to K and electric field",
    added_in_version="before_v9",
    text=r"""
If is equal to 1, activates computation of mixed second derivatives of wavefunctions with respect to
wavevector and electric field (ipert = natom+11 is activated). This is not strictly a response function
but is a needed auxiliary quantity in the calculations of 3rd-order derivatives of the energy
(non-linear response) if [[usepead]] == 0. The directions for the derivatives are determined by
[[rf2_pert1_dir]], [[rf2_pert2_dir]] and [[prepanl]] in the same way than [[rf2_dkdk]].
""",
),

Variable(
    abivarname="rf2_pert1_dir",
    varset="dfpt",
    vartype="integer",
    topics=['DFPT_useful'],
    dimensions=[3],
    defaultval=[1, 1, 1],
    mnemonics="Response Function (2nd order Sternheimer equation): 1st PERTurbation DIRection",
    added_in_version="before_v9",
    text=r"""
Gives the directions of the 1st perturbation to be considered when solving the 2nd order Sternheimer equation.
The three elements corresponds to the three primitive vectors, either in real
space (phonon calculations), or in reciprocal space ($\,d/ \,d k$, homogeneous
electric field, homogeneous magnetic field calculations).
If equal to 1, the 2nd order wavefunctions, as defined by [[rf2_dkdk]] or [[rf2_dkde]], are computed for the
corresponding direction. If 0, this direction is not considered.
See [[rf2_dkdk]] for more details.
""",
),

Variable(
    abivarname="rf2_pert2_dir",
    varset="dfpt",
    vartype="integer",
    topics=['DFPT_useful'],
    dimensions=[3],
    defaultval=[1, 1, 1],
    mnemonics="Response Function (2nd order Sternheimer equation): 2nd PERTurbation DIRection",
    added_in_version="before_v9",
    text=r"""
Gives the directions of the 2nd perturbation to be considered when solving the 2nd order Sternheimer equation.
The three elements corresponds to the three primitive vectors, either in real
space (phonon calculations), or in reciprocal space ($\,d/ \,d k$, homogeneous
electric field, homogeneous magnetic field calculations).
If equal to 1, the 2nd order wavefunctions, as defined by [[rf2_dkdk]] or [[rf2_dkde]], are computed for the
corresponding direction. If 0, this direction is not considered.
See [[rf2_dkdk]] for more details.
""",
),

Variable(
    abivarname="rfasr",
    varset="dfpt",
    vartype="integer",
    topics=['Phonons_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Response Function: Acoustic Sum Rule",
    added_in_version="before_v9",
    text=r"""
Control the evaluation of the acoustic sum rule in effective charges and
dynamical matrix at Gamma within a response function calculation (not active
at the level of producing the DDB, but at the level of the phonon
eigenfrequencies output).

  * 0 --> no acoustic sum rule imposed
  * 1 --> acoustic sum rule imposed for dynamical matrix at Gamma, and charge neutrality
    imposed with extra charge evenly distributed among atoms
  * 2 --> acoustic sum rule imposed for dynamical matrix at Gamma, and charge neutrality
    imposed with extra charge given proportionally to those atoms with the largest effective charge.

The treatment of the acoustic sum rule and charge neutrality sum rule is finer
at the level of the ANADDB utility, with the two independent input variables
[[anaddb:asr]] and [[anaddb:chneut]].
""",
),

Variable(
    abivarname="rfatpol",
    varset="dfpt",
    vartype="integer",
    topics=['DFPT_basic', 'Elastic_compulsory', 'Phonons_compulsory'],
    dimensions=[2],
    defaultval=[1, 1],
    mnemonics="Response Function: ATomic POLarisation",
    added_in_version="before_v9",
    text=r"""
Control the range of atoms for which displacements will be considered in
phonon calculations (atomic polarizations).
These values are only relevant to phonon response function calculations.
May take values from 1 to [[natom]], with [[rfatpol]](1)<=[[rfatpol]](2).
The atoms to be moved will be defined by the
do-loop variable iatpol:

  - do iatpol=[[rfatpol]](1),[[rfatpol]](2)

For the calculation of a full dynamical matrix, use [[rfatpol]](1)=1 and
[[rfatpol]](2)=[[natom]], together with [[rfdir]] 1 1 1. For selected
elements of the dynamical matrix, use different values of [[rfatpol]] and/or
[[rfdir]]. The name 'iatpol' is used for the part of the internal variable
ipert when it runs from 1 to [[natom]]. The internal variable ipert can also
assume values larger than [[natom]], denoting perturbations of electric field
or stress type (see [the DFPT help file](../guide/respfn)).
""",
),

Variable(
    abivarname="rfddk",
    varset="dfpt",
    vartype="integer",
    topics=['DFPT_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Response Function with respect to Derivative with respect to K",
    added_in_version="before_v9",
    text=r"""
Activates computation of derivatives of ground state wavefunctions with
respect to wavevectors. This is not strictly a response function but is a
needed auxiliary quantity in the electric field calculations (see [[rfelfd]]).
The directions for the derivatives are determined by [[rfdir]].

  * 0 --> no derivative calculation
  * 1 --> calculation of first derivatives of wavefunctions with respect to k points ($\,d/ \,d k$ calculation).
    The exact same functionality is provided by [[rfelfd]] = 2.
""",
),

Variable(
    abivarname="rfdir",
    varset="dfpt",
    vartype="integer",
    topics=['DFPT_compulsory', 'Elastic_compulsory', 'Phonons_compulsory'],
    dimensions=[3],
    defaultval=[0, 0, 0],
    mnemonics="Response Function: DIRections",
    added_in_version="before_v9",
    text=r"""
Gives the directions to be considered for response function calculations (also
for the Berry phase computation of the polarization, see the [[berryopt]]
input variable).
The three elements corresponds to the three primitive vectors, either in real
space (phonon calculations), or in reciprocal space ($\,d/ \,d k$, homogeneous
electric field, homogeneous magnetic field calculations). So, they generate a
basis for the generation of the dynamical matrix or the macroscopic dielectric
tensor or magnetic susceptibility and magnetic shielding, or the effective charge tensors.
If equal to 1, response functions, as defined by [[rfddk]], [[rfelfd]],
[[rfphon]], [[rfdir]] and [[rfatpol]], are to be computed for the
corresponding direction. If 0, this direction should not be considered.
""",
),

Variable(
    abivarname="rfelfd",
    varset="dfpt",
    vartype="integer",
    topics=['EffectiveMass_compulsory', 'DFPT_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Response Function with respect to the ELectric FielD",
    added_in_version="before_v9",
    text=r"""
Turns on electric field response function calculations. Actually, such
calculations requires first the non-self-consistent calculation of derivatives
with respect to k, independently of the electric field perturbation itself.

  * 0 --> no electric field perturbation
  * 1 --> full calculation, with first the derivative of ground-state wavefunction with respect to k ($\,d / \,d k$ calculation), by a non-self-consistent calculation, then the generation of the first-order response to an homogeneous electric field
  * 2 --> only the derivative of ground-state wavefunctions with respect to k
  * 3 --> only the generation of the first-order response to the electric field, assuming that the data on derivative of ground-state wavefunction with respect to k is available on disk.

!!! note
    Because the tolerances to be used for derivatives or homogeneous
    electric field are different, one often does the calculation of derivatives in
    a separate dataset, followed by calculation of electric field response as well as phonon.
    The options 2 and 3 proves useful in that context; also, in case a scissor
    shift is to be used, it is usually not applied for the $\,d / \,d k$ response).
""",
),

Variable(
    abivarname="rfmagn",
    varset="dfpt",
    vartype="integer",
    topics=['DFPT_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Response Function with respect to MAGNetic B-field perturbation",
    added_in_version="before_v9",
    text=r"""
[[rfmagn]] allows one to run response function calculations with respect to
external magnetic field if set to 1. Currently, orbital magnetism is not taken into
account and the perturbing potential has Zeeman form. For more details, see [[cite:Ricci2019]].
""",
),

Variable(
    abivarname="rfmeth",
    varset="dfpt",
    vartype="integer",
    topics=['DFPT_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Response Function METHod",
    added_in_version="before_v9",
    text=r"""
Selects method used in response function calculations. Presently, only abs([[rfmeth]]) = 1 is
allowed. This corresponds to storing matrix elements of the 2DTE computed using non-stationary expressions,
instead of stationary ones.

The difference between positive and negative values is rather technical. Very often, the symmetries can be used in such a way that
some matrix elements can be proven to be zero even without doing any computation. Positive values of [[rfmeth]] activate
this use of symmetries, while it is denied when [[rfmeth]] is negative. There is an indirect additional outcome of this,
as a symmetrization of the whole 2DTE is sometimes rendered possible when the additional knowledge of the zero matrix elements
is available. Thus, the results obtained for positive and negative values of [[rfmeth]] might slightly differ for non-zero elements of the 2DTE,
if they are computed in both cases.
""",
),

Variable(
    abivarname="rfphon",
    varset="dfpt",
    vartype="integer",
    topics=['DFPT_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Response Function with respect to PHONons",
    added_in_version="before_v9",
    text=r"""
It must be equal to 1 to run phonon response function calculations.
""",
),

Variable(
    abivarname="rfstrs",
    varset="dfpt",
    vartype="integer",
    topics=['DFPT_basic', 'Elastic_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Response Function with respect to STRainS",
    added_in_version="before_v9",
    text=r"""
Used to run strain response-function calculations (e.g. needed to get elastic
constants). Define, with [[rfdir]], the set of perturbations.

  * 0 --> no strain perturbation
  * 1 --> only uniaxial strain(s) (ipert=natom+3 is activated)
  * 2 --> only shear strain(s) (ipert=natom+4 is activated)
  * 3 --> both uniaxial and shear strain(s) (both ipert=natom+3 and ipert=natom+4 are activated)

See the possible restrictions on the use of strain perturbations, in the [[help:respfn]].
""",
),

Variable(
    abivarname="rfuser",
    varset="dfpt",
    vartype="integer",
    topics=['DFPT_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Response Function, USER-defined",
    added_in_version="before_v9",
    text=r"""
Available to the developers, to activate the use of ipert=natom+6 and
ipert=natom+7, two sets of perturbations that the developers can define.

  * 0 --> no computations for ipert=natom+6 or ipert=natom+7
  * 1 --> response with respect to perturbation natom+6 will be computed
  * 2 --> response with respect to perturbation natom+7 will be computed
  * 3 --> responses with respect to perturbations natom+6 and natom+7 will be computed

!!! important

    In order to define and use correctly the new perturbations, the developer
    might have to include code lines or additional routines at the level of the
    following routines: dfpt_cgwf.F90, dfpt_dyout.F90, dfpt_symph.F90,
    dfpt_dyout.F90, dfpt_etot.F90, littlegroup_pert.F90, dfpt_looppert.F90,
    dfpt_mkcor.F90, dfpt_nstdy.F90, dfpt_nstwf.F90, respfn.F90, dfpt_scfcv.F90,
    irreducible_set_pert.F90, dfpt_vloca.F90, dfpt_vtorho.F90, dfpt_vtowfk.F90. In
    these routines, the developer should pay a particular attention to the rfpert
    array, defined in the routine respfn (in m_respfn_driver.F90), as well as to the ipert local variable.
""",
),

Variable(
    abivarname="rhoqpmix",
    varset="gw",
    vartype="real",
    topics=['GW_useful'],
    dimensions="scalar",
    defaultval=1.0,
    mnemonics="RHO QuasiParticle MIXing",
    added_in_version="before_v9",
    text=r"""
For self-consistent GW runs, [[rhoqpmix]] sets the mixing coefficient between
the new and the previous electronic densities. This mixing damps the spurious
oscillations in the Hartree potential when achieving self-consistency.
[[rhoqpmix]] is meaningful only when doing self-consistency on the
wavefunctions with [[gwcalctyp]] >= 20.
""",
),

Variable(
    abivarname="rprim",
    varset="basic",
    vartype="real",
    topics=['UnitCell_basic'],
    dimensions=[3, 3],
    defaultval=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    mnemonics="Real space PRIMitive translations",
    characteristics=['[[EVOLVING]]'],
    commentdims="Internally, it is represented as rprim(3,3,[[nimage]])",
    added_in_version="before_v9",
    text=r"""
Give the three dimensionless primitive translations in
real space, to be rescaled by [[acell]] and [[scalecart]].
The three first numbers are the coordinates of the first vector, the next three numbers are the coordinates
of the second, and the last three the coordinates of the third.
It is [[EVOLVING]] only if [[ionmov]] == 2 or 22 and [[optcell]]/=0, otherwise it is fixed.

If the Default is used, that is, **rprim** is the unity matrix, the three
dimensionless primitive vectors are three unit vectors in cartesian
coordinates. The coordinates (and hence the length) of each vector will be (possibly)
multiplied by the corresponding [[acell]] value, then (possibly) stretched along the
cartesian coordinates by the corresponding [[scalecart]] value, to give the dimensional primitive vectors, called [[rprimd]].

In the general case, the dimensional cartesian coordinates of the crystal
primitive translations R1p, R2p and R3p, see [[rprimd]], are

  * R1p(i) = [[scalecart]](i) x [[rprim]](i,1) x [[acell]](1)
  * R2p(i) = [[scalecart]](i) x [[rprim]](i,2) x [[acell]](2)
  * R3p(i) = [[scalecart]](i) x [[rprim]](i,3) x [[acell]](3)

where i=1,2,3 is the component of the primitive translation (i.e. x, y, and z).

The [[rprim]] variable, scaled by [[scalecart]], is thus used to define
directions of the primitive vectors, that will be multiplied (so keeping the
direction unchanged) by the appropriate length scale [[acell]](1),
[[acell]](2), or [[acell]](3), respectively to give the dimensional primitive
translations in real space in cartesian coordinates.
Presently, it is requested that the mixed product (R1xR2).R3 is positive. If
this is not the case, simply exchange a pair of vectors.

To be more specific, keeping the default value of [[scalecart]] = 1 to simplify
the matter, [[rprim]] 1 2 3 4 5 6 7 8 9 corresponds to input of the three
primitive translations R1=(1,2,3) (to be multiplied by [[acell]](1)),
R2=(4,5,6) (to be multiplied by [[acell]](2)), and R3=(7,8,9) (to be
multiplied by [[acell]](3)).
Note carefully that the first three numbers input are the first column of
[[rprim]], the next three are the second, and the final three are the third.
This corresponds with the usual Fortran order for arrays. The matrix whose
columns are the reciprocal space primitive translations is the inverse
transpose of the matrix whose columns are the direct space primitive
translations.

Alternatively to [[rprim]], directions of dimensionless primitive vectors can
be specified by using the input variable [[angdeg]]. This is especially useful
for hexagonal lattices (with 120 or 60 degrees angles). Indeed, in order for
symmetries to be recognized, rprim must be symmetric up to [[tolsym]] (10
digits by default), inducing a specification such as

      rprim  0.86602540378  0.5  0.0
            -0.86602540378  0.5  0.0
             0.0            0.0  1.0

that can be avoided thanks to [[angdeg]]:

      angdeg 90 90 120

Note that the following might work as well:


      rprim  sqrt(0.75)  0.5  0.0
            -sqrt(0.75)  0.5  0.0
             0.0         0.0  1.0

Although the use of [[scalecart]] or [[acell]] is rather equivalent when the
primitive vectors are aligned with the cartesian directions, it is not the
case for non-orthogonal primitive vectors. In particular, beginners often make
the error of trying to use [[acell]] to define primitive vectors in face-
centered tetragonal lattice, or body-centered tetragonal lattice, or similarly
in face or body-centered orthorhombic lattices. Let us take the example of a
body-centered tetragonal lattice, that might be defined using the following
("a" and "c" have to be replaced by the appropriate conventional cell vector length):

      rprim  "a"      0        0
              0      "a"       0
             "a/2"   "a/2"    "c/2"
    acell 3*1     scalecart 3*1    !  ( These are default values)


The following is a valid, alternative way to define the same primitive vectors:

      rprim   1        0       0
              0        1       0
              1/2      1/2     1/2
    scalecart  "a"  "a"  "c"
    acell 3*1    !  ( These are default values)

Indeed, the cell has been stretched along the cartesian coordinates, by "a",
"a" and "c" factors.

At variance, the following is **WRONG**:

      rprim   1       0       0
              0       1       0
              1/2     1/2     1/2
    acell  "a"  "a"  "c"    !   THIS IS WRONG
    scalecart 3*1    !  ( These are default values)

Indeed, the latter would correspond to:

      rprim  "a"      0       0
              0      "a"      0
             "c/2"   "c/2"   "c/2"
    acell 3*1     scalecart 3*1    !  ( These are default values)

Namely, the third vector has been rescaled by "c". It is not at all in the
center of the tetragonal cell whose basis vectors are defined by the scaling factor "a".
As another difference between [[scalecart]] or [[acell]], note that
[[scalecart]] is [[INPUT_ONLY]]: its content will be immediately applied to
rprim, at parsing time, and then scalecart will be set to the default values
(3*1). So, in case [[scalecart]] is used, the echo of [[rprim]] in the output
file is not the value contained in the input file, but the value rescaled by [[scalecart]].
""",
),

Variable(
    abivarname="rprimd",
    varset="basic",
    vartype="real",
    topics=['UnitCell_internal'],
    dimensions=[3, 3],
    mnemonics="Real space PRIMitive translations, Dimensional",
    characteristics=['[[INTERNAL_ONLY]]', '[[EVOLVING]]'],
    commentdims="Internally, it is represented as rprimd(3,3,[[nimage]]).",
    added_in_version="before_v9",
    text=r"""
This internal variable gives the dimensional real space primitive vectors,
computed from [[acell]], [[scalecart]], and [[rprim]].

  * R1p(i) = [[rprimd]](i,1) = [[scalecart]](i) x [[rprim]](i,1) x [[acell]](1) for i=1,2,3 (x,y,and z)
  * R2p(i) = [[rprimd]](i,2) = [[scalecart]](i) x [[rprim]](i,2) x [[acell]](2) for i=1,2,3
  * R3p(i) = [[rprimd]](i,3) = [[scalecart]](i) x [[rprim]](i,3) x [[acell]](3) for i=1,2,3

It is [[EVOLVING]] only if [[ionmov]] == 2 or 22 and [[optcell]]/=0, otherwise it is fixed.
""",
),

Variable(
    abivarname="scalecart",
    varset="basic",
    vartype="real",
    topics=['UnitCell_useful'],
    dimensions=[3],
    defaultval=MultipleValue(number=3, value=1),
    mnemonics="SCALE CARTesian coordinates",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Gives the scaling factors of cartesian coordinates by which dimensionless
primitive translations (in "[[rprim]]") are to be multiplied.
See the [[rprim]] input
variable, the [[acell]] input variable, and the associated internal [[rprimd]]
internal variable.
Especially useful for body-centered and face-centered tetragonal lattices, as
well as body-centered and face-centered orthorhombic lattices, see [[rprimd]].
Note that this input variable is [[INPUT_ONLY]]: its content will be
immediately applied to rprim, at parsing time, and then scalecart will be set
to the default values. So, it will not be echoed.
""",
),

Variable(
    abivarname="scphon_supercell",
    varset="gstate",
    vartype="integer",
    topics=['DFPT_expert'],
    dimensions=[3],
    defaultval=[1, 1, 1],
    mnemonics="Self Consistent PHONon SUPERCELL",
    added_in_version="before_v9",
    text=r"""
Give extent, in number of primitive unit cells, of the supercell being used
for a self-consistent phonon calculation. Presumes the phonon frequencies and
eigenvectors have been calculated in the original primitive unit cell, on a
grid of q-points which corresponds to the supercell in the present calculation.
Experimental.
""",
),

Variable(
    abivarname="scphon_temp",
    varset="gstate",
    vartype="real",
    topics=['DFPT_expert'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="Self Consistent PHONon TEMPerature",
    characteristics=['[[ENERGY]]'],
    added_in_version="before_v9",
    text=r"""
Temperature which is imposed on phonon distribution, in the self-consistent
scheme of [[cite:Souvatzis2008]]. Determines the extent of the
finite displacements used, and consequent anharmonic effects. Experimental.
""",
),

Variable(
    abivarname="shiftk",
    varset="basic",
    vartype="real",
    topics=['k-points_useful'],
    dimensions=[3, '[[nshiftk]]'],
    defaultval=ValueWithConditions({'[[nshiftk]]>1': None, 'defaultval': [0.5, 0.5, 0.5]}),
    mnemonics="SHIFT for K points",
    added_in_version="before_v9",
    text=r"""
It is used only when [[kptopt]] >= 0, and must be defined if [[nshiftk]] is
larger than 1.
[[shiftk]](1:3,1:[[nshiftk]]) defines [[nshiftk]] shifts of the homogeneous
grid of k points based on [[ngkpt]] or [[kptrlatt]].
The shifts induced by [[shiftk]] corresponds to the reduced coordinates in the
coordinate system defining the k-point lattice. For example, if the k point
lattice is defined using [[ngkpt]], the point whose reciprocal space reduced
coordinates are ( [[shiftk]](1,ii)/[[ngkpt]](1) [[shiftk]](2,ii)/[[ngkpt]](2)
[[shiftk]](3,ii)/[[ngkpt]](3) ) belongs to the shifted grid number ii.

The user might rely on ABINIT to suggest suitable and efficient combinations
of [[kptrlatt]] and [[shiftk]]. The procedure to be followed is described with
the input variables [[kptrlen]]. In what follows, we suggest some interesting
values of the shifts, to be used with even values of [[ngkpt]]. This list is
much less exhaustive than the above-mentioned automatic procedure.

1) When the primitive vectors of the lattice do NOT form a FCC or a BCC
lattice, the default (shifted) Monkhorst-Pack grids are formed by using
[[nshiftk]] = 1 and [[shiftk]] 0.5 0.5 0.5. This is often the preferred k point
sampling, as the shift improves the sampling efficiency. However, it can also
break symmetry, if the 111 direction is not an axis of rotation, e.g. in
tetragonal or hexagonal systems. Abinit will complain about this breaking, and
you should adapt [[shiftk]]. For a non-shifted Monkhorst-Pack grid, use
[[nshiftk]] = 1 and [[shiftk]] 0.0 0.0 0.0, which will be compatible with all
symmetries, and is necessary for some features such as k-point interpolation.

2) When the primitive vectors of the lattice form a FCC lattice, with [[rprim]]

      0.0 0.5 0.5
      0.5 0.0 0.5
      0.5 0.5 0.0

the (very efficient) usual Monkhorst-Pack sampling will be generated by using
[[nshiftk]] =  4 and [[shiftk]]

      0.5 0.5 0.5
      0.5 0.0 0.0
      0.0 0.5 0.0
      0.0 0.0 0.5

3) When the primitive vectors of the lattice form a BCC lattice, with [[rprim]]

      -0.5  0.5  0.5
       0.5 -0.5  0.5
       0.5  0.5 -0.5

the usual Monkhorst-Pack sampling will be generated by using [[nshiftk]] =  2
and [[shiftk]]

      0.25  0.25  0.25
     -0.25 -0.25 -0.25

However, the simple sampling [[nshiftk]] = 1 and [[shiftk]] 0.5 0.5 0.5 is excellent.

4) For hexagonal lattices with hexagonal axes, e.g. [[rprim]]

      1.0  0.0       0.0
     -0.5  sqrt(3)/2 0.0
      0.0  0.0       1.0

one can use [[nshiftk]] =  1 and [[shiftk]] 0.0 0.0 0.5

In rhombohedral axes, e.g. using [[angdeg]] 3*60., this corresponds to
[[shiftk]] 0.5 0.5 0.5, to keep the shift along the symmetry axis.
""",
),

Variable(
    abivarname="shiftq",
    varset="gstate",
    vartype="real",
    topics=['q-points_useful'],
    dimensions=[3, '[[nshiftq]]'],
    defaultval=ValueWithConditions({'[[nshiftq]]>1': None, 'defaultval': [0.5, 0.5, 0.5]}),
    mnemonics="SHIFT for Q points",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
It is used only when [[qptopt]] >= 0, and must be defined if [[nshiftq]] is larger than 1.
[[shiftq]](1:3,1:[[nshiftq]]) defines [[nshiftq]] shifts of the homogeneous
grid of q points based on [[ngqpt]] or [[qptrlatt]].

See [[shiftk]] for more information on the definition, use, and suitable
values for these shifts.
""",
),

Variable(
    abivarname="signperm",
    varset="rlx",
    vartype="integer",
    topics=['MolecularDynamics_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="SIGN of PERMutation potential",
    added_in_version="before_v9",
    text=r"""
+1 favors alternation of species -1 favors segregation
""",
),

Variable(
    abivarname="slabwsrad",
    varset="gstate",
    vartype="real",
    topics=['Artificial_expert'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="jellium SLAB Wigner-Seitz RADius",
    characteristics=['[[LENGTH]]'],
    added_in_version="before_v9",
    text=r"""
Fix the bulk-mean positive charge density $n_{bulk}$ of a jellium slab (if the
latter is employed, e.g. [[jellslab]]/=0). Often called $r_s$ (see for example
[[cite:Lang1970]]), [[slabwsrad]] is the radius of a
sphere which has the same volume as the average volume per particle in a
homogeneous electron gas with density $n_{bulk}$, so:

\begin{equation}
      \frac{1}{n_{bulk}} = \frac{4 \pi}{3} [[slabwsrad]]^3 \nonumber
\end{equation}

For example, the bulk aluminum fcc lattice constant is $a$=4.0495 Angstroms
[WebElements](https://www.webelements.com/), each cubic centered cell includes 4 Al atoms and each atom
has 3 valence electrons, so the average volume per electron is $a^3/12$=37.34
Bohr$^3$ which has to be equal to $\frac{4 \pi}{3} r_s^3$. Consequently Al has approximately
$r_s$=2.07 Bohr, while for example magnesium has $r_s$=2.65 Bohr, sodium 3.99 Bohr.
By default, given in Bohr atomic units (1 Bohr=0.5291772108 Angstroms).
""",
),

Variable(
    abivarname="slabzbeg",
    varset="gstate",
    vartype="real",
    topics=['Artificial_expert'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="jellium SLAB BEGinning edge along the z-direction",
    added_in_version="before_v9",
    text=r"""
Define the edges of the jellium slab (if used, so if [[jellslab]]/=0) along
z, namely the slab starts at a point along z which is expressed in Bohr by
[[slabzbeg]] and it ends at a point expressed in Bohr by [[slabzend]].
The z-direction is parallel to the third crystal primitive lattice vector which has
to be orthogonal to the other ones, so the length of the cell along z is
[[rprimd]](3,3). In addition [[slabzbeg]] and [[slabzend]] have to be such that:

      0 ≤ [[slabzbeg]]  < [[slabzend]] ≤ [[rprimd]](3,3)

Together with [[slabwsrad]] they define the jellium positive charge density
distribution $n_{+}(x,y,z)$ in this way:

\begin{eqnarray}
      n_{+}(x,y,z) &=& n_{bulk} \quad \text{if} \quad [[slabzbeg]]  \leq z \leq [[slabzend]]  \nonumber\\
                &=& 0       \quad \text{otherwise}                           \nonumber
\end{eqnarray}

so the positive charge density is invariant along the xy plane as well as the
electrostatic potential generated by it.
""",
),

Variable(
    abivarname="slabzend",
    varset="gstate",
    vartype="real",
    topics=['Artificial_expert'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="jellium SLAB ENDing edge along the z-direction",
    added_in_version="before_v9",
    text=r"""
Define the edges of the jellium slab (if used, so if [[jellslab]]/=0) along
z, namely the slab starts at a point along z which is expressed in Bohr by
[[slabzbeg]] and it ends at a point expressed in Bohr by [[slabzend]].
The z-direction is parallel to the third crystal primitive lattice vector which has
to be orthogonal to the other ones, so the length of the cell along z is
[[rprimd]](3,3). In addition [[slabzbeg]] and [[slabzend]] have to be such that:

      0 ≤ [[slabzbeg]] < [[slabzend]]  ≤ [[rprimd]](3,3)

Together with [[slabwsrad]] they define the jellium positive charge density
distribution $n_{+}(x,y,z)$ in this way:

\begin{eqnarray}
      n_{+}(x,y,z) &=& n_{bulk} \quad  \text{if} \quad [[slabzbeg]] \leq z \leq [[slabzend]] \nonumber \\
                   &=& 0        \quad  \text{otherwise}                                    \nonumber
\end{eqnarray}

so the positive charge density is invariant along the xy plane as well as the
electrostatic potential generated by it.
""",
),

Variable(
    abivarname="slk_rankpp",
    varset="gstate",
    vartype="integer",
    topics=['parallelism_expert'],
    dimensions="scalar",
    defaultval=[1000],
    mnemonics="ScaLapacK matrix RANK Per Process",
    added_in_version="before_v9",
    text=r"""
This variable controls how the number of processes to be used in Scalapack diagonalization algorithm: [[np_slk]] will be calculated according to this value.
This value is the matrix rank that each process will hold for the diagonalization.
For a 1000x1000 matrix with default value, scalapack won't be used (Lapack will be used).
For a 2000x2000 matrix with default value, scalapack will be used with 2000/1000=2 MPI processes.
For a 2000x2000 matrix with a slk_rank=500, scalapack will be used with 2000/500=4 MPI processes.
In case of hybrid MPI+OpenMP, the number of thread is also taken into account.

***WARNING*** None of the available scalapack library are thread-safe  (2018). Therefore using both scalapack *and* OpenMP is highly unpredictable.
Furthermore, using multithreaded linear algebra library (MKL ACML...) is more efficient than pure MPI scalapack.

Usually it is better to define this variable and let the code do the rest.
""",
),

Variable(
    abivarname="smdelta",
    varset="dfpt",
    vartype="integer",
    topics=['TDepES_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="SMeared DELTA function",
    added_in_version="before_v9",
    text=r"""
When [[smdelta]] in non-zero, it will trigger the calculation of the imaginary
part of the second-order electronic eigenvalues, which can be related to the
electronic lifetimes. The delta function is evaluated using:

  * when [[smdelta]] == 1, Fermi-Dirac smearing: $\frac{0.25}{(cosh(\frac{x}{2.0}))^2}$
  * when [[smdelta]] == 2, Cold smearing by Marzari using the parameter $a=-0.5634$ (minimization of the bump): $\frac{e^{-x^2}}{\sqrt{\pi}}\left(1.5+x(-a\ 1.5+x(-1.0+a\ x))\right)$
  * when [[smdelta]] == 3, Cold smearing by Marzari using the parameter $a=-0.8165$ (monotonic function in the tail): as 2 but different $a$
  * when [[smdelta]] == 4, Smearing of Methfessel and Paxton ([[cite:Methfessel1989]]) with Hermite polynomial of degree 2, corresponding to "Cold smearing" of N. Marzari with $a=0$ (so, same smeared delta function as smdelta=2, with different $a$).
  * when [[smdelta]] == 5, Gaussian smearing: $\frac{e^{-x^2}}{\sqrt{\pi}}$
""",
),

Variable(
    abivarname="so_psp",
    varset="gstate",
    vartype="integer",
    topics=['spinpolarisation_useful'],
    dimensions=['[[npsp]]'],
    defaultval=MultipleValue(number='[[npsp]]', value=1),
    mnemonics="Spin-Orbit treatment for each PSeudoPotential",
    requires="[[nspinor]] == 2 and [[usepaw]] == 0",
    added_in_version="before_v9",
    text=r"""
For each type of atom (each pseudopotential), specify the treatment of spin-orbit
interaction (if [[nspinor]] == 2 and Norm-conserving pseudopotentials i.e. [[usepaw]] == 0)
For PAW calculations with SOC, please refer to [[pawspnorb]].

  * If 0: no spin-orbit interaction, even if [[nspinor]] = 2
  * If 1: treat spin-orbit as specified in the pseudopotential file.
  * If 2: treat spin-orbit in the HGH form (usual form, although not allowed for all pseudopotentials)
  * If 3: treat spin-orbit in the HFN form (Hemstreet-Fong-Nelson) (actually, not implemented).

For typical usage, the default value is OK. If the spin-orbit needs to be
turned off for one atom, 0 might be relevant. Note however, that the code will
stop if [[nspinor]] == 2 is used and one of the pseudopotential does not contain
the information about the spin-orbit interaction (this is the case for some
old pseudopotentials). Indeed, for spinorial calculations, turning off the
spin-orbit interaction is unphysical, and also does not save CPU time.
It should only be done for test purposes

Note that if [[nspinor]] == 1, the spin-orbit cannot be treated anyhow, so the
value of [[so_psp]] is irrelevant.

Prior to v5.4, the input variable **so_typat** was used, in place of
[[so_psp]]. Because the values 0 and 1 have been switched between [[so_psp]]
and **so_typat**, it was dangerous to continue to allow the use of **so_typat**.
""",
),

Variable(
    abivarname="spbroad",
    varset="gw",
    vartype="real",
    topics=['Susceptibility_useful'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="SPectral BROADening",
    characteristics=['[[ENERGY]]'],
    requires="[[optdriver]] == 3 and [[spmeth]] == 2",
    added_in_version="before_v9",
    text=r"""
When a screening calculation ([[optdriver]] == 3) uses a spectral representation
of the irreducible polarizability in which the delta function is replaced by
the gaussian approximant ([[spmeth]] == 2), the standard deviation of the
gaussian is given by [[spbroad]].
""",
),

Variable(
    abivarname="spgaxor",
    varset="geo",
    vartype="integer",
    topics=['SmartSymm_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="SPace Group: AXes ORientation",
    added_in_version="before_v9",
    text=r"""
It is taken into account only when [[spgroup]]/=0; it allows one to define the
axes orientation for the specific space groups for which this is needed.
Trigonal groups (number 146,148,155,160,161,166,167):

  * 1 represents the hexagonal axes
  * 2 represents the rhombohedral axes

Orthorhombic space groups: there are six possibilities corresponding to the
possible axes permutations

  * 1 abc -> abc
  * 2 abc -> cab
  * 3 abc -> bca
  * 4 abc -> acb
  * 5 abc -> bac
  * 6 abc -> cba

Monoclinic: there are 3 or 9 possibilities depending on the space group.
For more details see the [[help:spacegroup]].
In the log/output file the notation used to describe the monoclinic
groups is for example: 15:c1, A2/a_c = C2/c where,

  * 15 represents the space group number,
  * c1 the orientation as it appears on the web page,
  * A is the real Bravais type lattice,
  * 2/a the existent symmetry elements,
  * _c marks the orientation of the two-fold axis or of the mirror plane,
  * C2/c represents the parent space group.

How to determine which [[spgaxor]] you need:

  1. check the reduced positions you have, for more symmetric positions, e.g. 1/2 1/4 3/4 etc...
     Let us say your symmetric positions are in the first coordinate (a axis) and you are using spgroup 62.
  2. look up the raw space group Wyckoff positions on
     [the Bilbao server](http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-wp-list) to see where they put
     the corresponding symmetric positions. For spgroup 62 Bilbao puts the 1/4 3/4 in the second coordinate, ie along the b axis.
  3. in this case you need to swap the axes from the original abc order to a new order where the
     Bilbao axis (b) is in the first position. In this case you have 2 possibilities, [[spgaxor]] 3 or 5.
     If you have more than one highly symmetric coordinate you may have only a single possibility.
""",
),

Variable(
    abivarname="spgorig",
    varset="geo",
    vartype="integer",
    topics=['SmartSymm_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="SPace Group: ORIGin",
    requires="[[spgroup]]!=0",
    added_in_version="before_v9",
    text=r"""
Gives the choice of origin for the axes system.
It is defined according to the origin choice in the International Tables of
Crystallography.
It applies only to the space groups 48, 50, 59, 70, 85, 86, 88, 125, 126, 129,
130, 133, 134, 137, 141, 142, 201, 203, 222, 224, 227, 228.
For more details see the [[help:spacegroup]].
""",
),

Variable(
    abivarname="spgroup",
    varset="geo",
    vartype="integer",
    topics=['crystal_useful', 'UnitCell_useful', 'SmartSymm_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="SPace GROUP number",
    added_in_version="before_v9",
    text=r"""
Gives the Fedorov space group number of the system.
It should be between 1 and 230, see [[help:spacegroup]].
Alternatively, if [[spgroup]] is 0, the code assumes that all the symmetries are input
through the [[symrel]] matrices and the [[tnons]] vectors, or obtained from
the symmetry finder (the default when [[nsym]] == 0). Then, ABINIT computes the value of [[spgroup]].

The list of symmetry operations that is available when [[spgroup]] is defined can be used to obtain all the
atoms in the unit cell, starting from the asymmetric unit cell, see [[natrd]].

The references for the numbering of space groups, and their list of symmetry operations is:

  * International Tables for Crystallography [[cite:Hahn1983]]
  * The mathematical theory of symmetry in solids, Representation theory for point groups and space groups [[cite:Bradley1972]]

Related input variables: [[symrel]], [[tnons]], [[symafm]], [[spgroupma]],
""",
),

Variable(
    abivarname="spgroupma",
    varset="geo",
    vartype="integer",
    topics=['spinpolarisation_useful', 'SmartSymm_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="SPace GROUP number defining a MAgnetic space group",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
This input variable might be used to define a Shubnikov magnetic space group
(anti-ferromagnetic space group). The user is advised to consult [[cite:Bradley1972]].

A Shubnikov type IV magnetic space group contains a specific type of symmetry operation,
which is a translation in real space followed by a spin flip.
Such a Shubnikov type IV magnetic space group might be defined by its Fedorov space
group (set of spatial symmetries that do not change the magnetization), and an
additional magnetic space group number [[spgroupma]].

A Shubnikov type III magnetic space group does not contain a translation in real space
wollowed by a spin flip. It might be defined by its Fedorov
space group (set of all spatial symmetries, irrespective of their magnetic
action), and an additional magnetic space group number [[spgroupma]].
For the additional number [[spgroupma]], we follow the definition of Table 7.4
of the above-mentioned [[cite:Bradley1972]].

Thus, one way to specify a Shubnikov type IV magnetic space group, is to define
both [[spgroup]] and [[spgroupma]].
For example, the group P2_1/c_prime has [[spgroup]] = 14 and [[spgroupma]] = 78.
Alternatively, for Shubnikov type IV magnetic groups, one might define [[spgroup]]
and [[genafm]]. For both the types III and IV, one might define by hand the set
of symmetries, using [[symrel]], [[tnons]] and [[symafm]].

Note that the meaning of the spin-flip operation of symmetry is different in the [[nspden]]=2 or in the [[nspden]]=4
case, see detailed explanations in the section on the [[symafm]] input variable. Thus, the same atomic
positions and [[spinat]] vectors might yield different [[symafm]] values depending on [[nspden]],
and thus different Shubnikov magnetic space groups.
""",
),

Variable(
    abivarname="spinat",
    varset="gstate",
    vartype="real",
    topics=['spinpolarisation_basic', 'crystal_useful', 'MagMom_useful', 'ConstrainedDFT_useful'],
    dimensions=ValueWithConditions({'[[natrd]]<[[natom]]': '[3, [[natrd]] ]', 'defaultval': '[3, [[natom]] ]'}),
    defaultval=0.0,
    mnemonics="SPIN for AToms",
    added_in_version="before_v9",
    text=r"""
Gives the initial electronic spin-magnetization for each atom, in unit of $\hbar/2$,
as well as, in case of fixed magnetization calculations (see [[constraint_kind]] and [[magconon]]), the target value of the magnetization.

Note that if [[nspden]] = 2, the z-component must be given for each atom, in
triplets (0 0 z-component).
For example, the electron of an hydrogen atom can be spin up (0 0 1.0) or spin
down (0 0 -1.0).

This value is only used to create the first exchange and correlation
potential.
It is not checked against the initial occupation numbers [[occ]] for each spin
channel.
It is meant to give an easy way to break the spin symmetry, and to allow to
find stable local spin fluctuations, for example: antiferromagnetism, or the
spontaneous spatial spin separation of elongated H$_2$ molecule.

  * If the atom manipulator is used, [[spinat]] will be related to the preprocessed set of atoms,
  generated by the atom manipulator. The user must thus foresee the effect of this atom manipulator (see [[objarf]]).

  * If the atom manipulator is not used, and the symmetries are not specified by the user ([[nsym]] = 0),
spinat will be used, if present, to determine the anti-ferromagnetic characteristics of the symmetry operations, see [[symafm]].
In case of collinear antiferromagnetism ([[nsppol]] = 1, [[nspinor]] = 1,
[[nspden]] = 2), these symmetries are used to symmetrize the density.
In case of non-collinear magnetism ([[nsppol]] = 1, [[nspinor]] = 2,
[[nspden]] = 4), they are also used to symmetrize the density. In the latter
case, this strongly constrains the magnetization (imposing its direction). If
the user want to let all degrees of freedom of the magnetization evolve, it is
then recommended to put [[nsym]] = 1.

* If the symmetries are specified, and the irreducible set of atoms is specified, the anti-ferromagnetic characteristics of the symmetry operations [[symafm]] will be used to generate [[spinat]] for all the non-irreducible atoms.

* In the case of PAW+U calculations using the [[dmatpawu]] initial occupation matrix, and if [[nspden]] = 4, [[spinat]] is also used to determine the direction of the integrated magnetization matrix.
""",
),

Variable(
    abivarname="spinmagntarget",
    varset="ffield",
    vartype="real",
    topics=['spinpolarisation_useful'],
    dimensions="scalar",
    defaultval=-99.99,
    mnemonics="SPIN-MAGNetization TARGET",
    added_in_version="before_v9",
    text=r"""
This input variable is active only in the [[nsppol]] = 2 case. If
[[spinmagntarget]] is not the "magic" value of -99.99, the spin-
magnetization of the primitive cell will be fixed (or optimized, if it is not
possible to impose it) to the value of [[spinmagntarget]], in Bohr magneton
units (for an Hydrogen atom, it is 1).
If [[occopt]] is a metallic one, the Fermi energies for spin up and spin down
are adjusted to give the target spin-polarisation (this is equivalent to an
exchange splitting). If [[occopt]] = 1 and [[nsppol]] = 2, the occupation numbers
for spin up and spin down will be adjusted to give the required spin-
magnetization (occupation numbers are identical for all k-points, with
[[occopt]] = 1). The definition of [[spinmagntarget]] is actually requested in
this case, except for the single isolated Hydrogen atom.
If [[spinmagntarget]] is the default one, the spin-magnetization will not be
constrained, and will be determined self-consistently, by having the same spin
up and spin down Fermi energy in the metallic case, while for the other cases,
there will be no spin-magnetization, except for an odd number of electrons if
[[occopt]] = 1 and [[nsppol]] = 2.

!!! note
    For the time being, only the spin down Fermi energy is written out in
    the main output file. In the fixed magnetic moment case, it differs from the
    spin up Fermi energy.
""",
),

Variable(
    abivarname="spmeth",
    varset="gw",
    vartype="integer",
    topics=['Susceptibility_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="SPectral METHod",
    requires="[[optdriver]] == 4",
    added_in_version="before_v9",
    text=r"""
The [[spmeth]] input variable defines the method used to calculate the
irreducible polarizability $\chi^{(0)}_{KS}$.

By default $\chi^{(0)}_{KS}$ is calculated employing the Adler-Wiser
expression ([[spmeth]] = 0) with a CPU effort that scales linearly with the
number of frequencies. This approach is convenient when few frequencies are
required, and is usually used in conjunction with plasmon-pole models in which
only one or two frequencies are calculated, according to the value of [[ppmodel]].
Unfortunately a calculation based on the Adler-Wiser expression might be quite
CPU demanding if the matrix elements of the self-energy operator are evaluated
by performing numerically the convolution defining the self-energy. The
integrand function, indeed, has poles above and below the real axis, and the
screened interaction has to be evaluated on a dense frequency mesh in order to
obtain accurate results.

In the spectral method ([[spmeth]] = 1 or 2) the irreducible polarizability is
expressed as the Hilbert transform of the imaginary part. The advantage in
using this approach consists in the fact that, once the spectral function is
known, the irreducible polarizability for an arbitrary frequency can be easily
obtained through inexpensive integrations. On the other hand, an accurate
evaluation of the imaginary part requires a dense frequency mesh due to the
presence of delta functions. Two different approaches can be used to
approximate these delta functions thus allowing the use of affordable
frequency grids.

Summarizing:

  * 0 -->  use Adler-Wiser expression to calculate $\chi^{(0)}_{KS}$
  * 1 -->  use the spectral method approximating the delta function with a triangular approximant as proposed in **REF TO BE ADDED**
  * 2 --> use spectral method but approximating the delta function with a Taylor expansion of the exponential as proposed in **REF TO BE ADDED**
""",
),

Variable(
    abivarname="spnorbscl",
    varset="paw",
    vartype="real",
    topics=['PAW_expert', 'spinpolarisation_useful'],
    dimensions="scalar",
    defaultval=1.0,
    mnemonics="SPin-ORBit SCaLing",
    requires="[[usepaw]] == 1 and [[pawspnorb]] >= 1",
    added_in_version="before_v9",
    text=r"""
Scaling of the spin-orbit interaction. The default values gives the first-
principles value, while other values are used for the analysis of the effect
of the spin-orbit interaction, but are not expected to correspond to any
physical situation.
""",
),

Variable(
    abivarname="stmbias",
    varset="gstate",
    vartype="real",
    topics=['STM_compulsory'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="Scanning Tunneling Microscopy BIAS voltage",
    characteristics=['[[ENERGY]]'],
    added_in_version="before_v9",
    text=r"""
Gives, in Hartree, the bias of the STM tip, with respect to the sample, in
order to generate the STM density map.
Used with positive [[iscf]], [[occopt]] = 7 (metallic, gaussian), [[nstep]] = 1,
and positive [[prtstm]], this value is used to generate a charge density map
from electrons close to the Fermi energy, in a (positive or negative) energy
range. Positive [[stmbias]] will lead to the inclusion of occupied (valence)
states only, while negative [[stmbias]] will lead to the inclusion of
unoccupied (conduction) states only.
Can be specified in Ha (the default), Ry, eV or Kelvin, since [[stmbias]] has
the [[ENERGY]] characteristics (0.001 Ha = 27.2113845 meV = 315.773 Kelvin).
With [[occopt]] = 7, one has also to specify an independent broadening [[tsmear]].
""",
),

Variable(
    abivarname="strfact",
    varset="rlx",
    vartype="real",
    topics=['GeoOpt_basic'],
    dimensions="scalar",
    defaultval=100,
    mnemonics="STRess FACTor",
    added_in_version="before_v9",
    text=r"""
The stresses multiplied by [[strfact]] will be treated like forces in the
process of optimization ([[ionmov]] = 2 or 22, non-zero [[optcell]]).
For example, the stopping criterion defined by [[tolmxf]] relates to these
scaled stresses.
""",
),

Variable(
    abivarname="string_algo",
    varset="rlx",
    vartype="integer",
    topics=['TransPath_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="STRING method ALGOrithm",
    added_in_version="before_v9",
    text=r"""
Relevant only when [[imgmov]] = 2 (String Method).
Gives the variant of the String Method method used.
Possible values can be:

  * 0 -->  **Original String Method**.
NOT YET IMPLEMENTED.
See [[cite:Weinan2002]]

  * 1 --> **Simplified String Method** with parametrization by **equal arc length**.
Instead of using the normal force (wrt the band), the full force is used; the
reparametrization is enforced by keeping the points of the string equally
spaced.
See [[cite:Weinan2007]]

  * 2 --> **Simplified String Method** with parametrization by **energy-weighted arc length**.
A variant of the Simplified String Method (like 1-); the reparametrization is
done by using energy-weight arc-lengths, giving a finer distribution near the saddle point.
See [[cite:Weinan2007]] and [[cite:Goodrow2009]]
""",
),

Variable(
    abivarname="strprecon",
    varset="rlx",
    vartype="real",
    topics=['ForcesStresses_useful', 'GeoOpt_useful'],
    dimensions="scalar",
    defaultval=1.0,
    mnemonics="STRess PRECONditioner",
    added_in_version="before_v9",
    text=r"""
This is a scaling factor to initialize the part of the Hessian related to the
treatment of the stresses (optimisation of the unit cell). In case there is an
instability, decrease the default value, e.g. set it to 0.1.
""",
),

Variable(
    abivarname="strtarget",
    varset="rlx",
    vartype="real",
    topics=['ForcesStresses_useful', 'GeoOpt_useful'],
    dimensions=[6],
    defaultval=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    mnemonics="STRess TARGET",
    added_in_version="before_v9",
    text=r"""
The components of the stress tensor must be stored according to: (1,1) -->1;
(2,2) --> 2; (3,3) --> 3; (2,3) --> 4; (3,1) --> 5; (1,2) -->6. The conversion factor
between Ha/Bohr**3 and GPa is: 1 Ha/Bohr**3 = 29421.033 GPa.
Not used if [[optcell]] == 0.
""",
),

Variable(
    abivarname="symafm",
    varset="gstate",
    vartype="integer",
    topics=['spinpolarisation_useful'],
    dimensions=['[[nsym]]'],
    defaultval=MultipleValue(number='[[nsym]]', value=1),
    mnemonics="SYMmetries, Anti-FerroMagnetic characteristics",
    added_in_version="before_v9",
    text=r"""
In case the material is magnetic, [[nspden]]=2 or 4, additional symmetry operations might appear,
that change the sign of the magnetization (spin-flip). They have been introduced by
Shubnikov in 1951 [[cite:Bradley1972]]. They can be used by ABINIT to decrease the CPU time, either
by decreasing the number of k-points or by suppressing the explicit treatment of one spin channel,
or by decreasing the number of perturbations in DFPT.

[[symafm]] should be set to +1 for all the usual symmetry operations, that do
not change the sign of the magnetization, while it should be set to -1 for the
magnetization-changing operations (spin-flip).
If the symmetry operations are not specified by the user in the input file,
that is, if [[nsym]] = 0, then ABINIT will use the values of [[spinat]] to
determine the content of [[symafm]].

The symmetries that can act on the magnetization can yield decreased CPU time  (and usually also memory decrease)
in the following cases:

  * antiferromagnetism ([[nsppol]] = 1, [[nspinor]] = 1, [[nspden]] = 2)
  * non-collinear magnetism ([[nsppol]] = 1, [[nspinor]] = 2, [[nspden]] = 4)

Also in the case [[nsppol]] = 2, [[nspinor]] = 1, [[nspden]] = 2  they might simply yield better accuracy (or faster convergence),
but there is no automatic gain of CPU time or memory, although it is not as clear cut as in the above cases.

IMPORTANT : The meaning of [[symafm]] is different in the [[nspden]] = 2 case (collinear magnetism),
and in the [[nspden]] = 4 case (non-collinear magnetism, with explicit treatment of magnetization as a vector).
Indeed in the first case, it is supposed that the magnetization vector is not affected by the real space symmetry operations
(so-called black and white symmetry groups).
By contrast, in the second case, the real space symmetry operations act on the magnetization vector.
The rationale for such different treatment comes from the fact that the treatment of spin-orbit coupling is incompatible with collinear magnetism [[nspden]]=2,
so there is no need to worry about it in this case. On the contrary, many calculations with [[nspden]]=2
will include spin-orbit coupling. The symmetry operations should thus act coherently on the spin-orbit coupling, which implies
that the real space operations should act also on the magnetization vector in the [[nspden]]=4 case. So, with
[[nspden]]=4, even with [[symafm]]=1,
symmetry operations might change the magnetization vector, e.g. possibly reverse it from one atom to another atom.
Still, when real space operations also act on the magnetization vector, nothing prevents to have ADDITIONAL "spin-flip" operations, which
is indeed then the meaning of [[symafm]]=-1 in the [[nspden]]=4 case.

Let's illustrate this with an example. Take an H$_2$ system, with the two H atoms quite distant from each other.
The electron on the first H atom might be 1s spin up, and the electron on the second atom might be 1s spin down.
With [[nspden]]=2, the inversion symmetry centered in the middle of the segment joining the two atoms will NOT act on the spin,
so that the actual symmetry operation that leaves the system invariant is an inversion ([[symrel]]= -1 0 0  0 -1 0  0 0 -1) accompanied
by a spin-flip with [[symafm]]=-1.
By contrast, with [[nspden]]=4, the inversion symmetry centered in the middle of the segment joining the two atoms will reverse the spin direction as well,
so that the proper symmetry operation is [[symrel]]= -1 0 0  0 -1 0  0 0 -1 but no additional spin-flip is needed to obtain a symmetry operation that leaves the system invariant, so that [[symafm]]=1.

Although this might seem confusing, ABINIT is able to recognise the correct symmetry operations from the available atomic coordinates, from [[spinat]],
and from [[nspden]], so that the user should hardly be affected by such different conventions. However, the use of [[ptgroupma]] and [[spgroupma]] to define
the antiferromagnetic operations of symmetry should be done carefully.

Related variables: [[symrel]], [[tnons]], [[ptgroupma]], [[spgroupma]].
""",
),

Variable(
    abivarname="symchi",
    varset="gw",
    vartype="integer",
    topics=['Susceptibility_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics=r"SYMmetryze $\chi_0$",
    characteristics=['[[DEVELOP]]'],
    requires="[[optdriver]] == 3",
    added_in_version="before_v9",
    text=r"""
The evaluation of the irreducible polarizability for a given q point requires
an integration over the Brillouin zone (BZ) which is approximated by a
discrete sum over k points. In principle the integrand function should be
evaluated for each k-point in the BZ, however it is possible to reduce the
number of points to be explicitly considered by taking advantage of symmetry
properties. The development input variable [[symchi]] is used to choose
between these two equivalent methods:

  * 0 -->  the summation over k points is performed considering **all** the points in the BZ (useful for testing and debugging).
  * 1 -->  the summation is restricted to the k points belonging to the irreducible wedge defined by the little group associated to the external vector q.
""",
),

Variable(
    abivarname="symdynmat",
    varset="eph",
    vartype="integer",
    topics=['Phonons_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="SYMmetrize the DYNamical MATrix",
    added_in_version="before_v9",
    text=r"""
If symdynmat is equal to 1, the dynamical matrix is symmetrized before the
diagonalization (same meaning as the corresponding anaddb variable). Note that
[[symdynmat]] == 1 will automatically enable the symmetrization of the electron-
phonon linewidths.
""",
),

Variable(
    abivarname="symmorphi",
    varset="dev",
    vartype="integer",
    topics=['crystal_expert', 'GW_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="SYMMORPHIc symmetry operation selection",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
With [[symmorphi]] = 1, symmetry operations with a non-symmorphic vector are
allowed. With [[symmorphi]] = 0, they are not allowed. In the latter case, if
the symmetry operations are specified in the input file, the code will stop
and print an error message if a non-symmorphic vector is encountered. By
contrast, if the symmetry operations are to be determined automatically (if
[[nsym]] = 0), then the set of symmetries will not include the non-symmorphic
operations.

Note: this feature exist because in a previous status of the GW calculations,
non-symmorphic symmetry operations could not be exploited. Thus, the k points
were restricted to the IBZ. In order to prepare GW calculations, and to
perform GW calculations, [[symmorphi]] = 0 was to be used, together with [[nsym]] = 0.
""",
),

Variable(
    abivarname="symrel",
    varset="basic",
    vartype="integer",
    topics=['crystal_useful'],
    dimensions=[3, 3, '[[nsym]]'],
    defaultval=ValueWithConditions({'[[nsym]] == 1': [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 'defaultval': None}),
    mnemonics="SYMmetry in REaL space",
    added_in_version="before_v9",
    text=r"""
Gives "[[nsym]]" 3x3 matrices expressing space group symmetries in terms of
their action on the direct (or real) space primitive translations.
It turns out that these can always be expressed as integers.
Always give the identity matrix even if no other symmetries hold, e.g.
[[symrel]] 1 0 0 0 1 0 0 0 1.
Also note that for this array, as for all others, the array elements are filled
in a columnwise order as is usual for Fortran.
The relation between the above symmetry matrices [[symrel]], expressed in the
basis of primitive translations, and the same symmetry matrices expressed in
cartesian coordinates, is as follows. Denote the matrix whose columns are the
primitive translations as R, and denote the cartesian symmetry matrix as S.
Then [[symrel]] = R(inverse) * S * R
where matrix multiplication is implied.
When the symmetry finder is used (see [[nsym]]), [[symrel]] will be computed automatically.
Also see the accompanying input variables [[tnons]] and [[symafm]], for the full definition of the symmetry operations.
Such variables are used to infer [[spgroup]], [[spgroupma]] and [[genafm]] if they are not user-defined.
""",
),

Variable(
    abivarname="symsigma",
    varset="gw",
    vartype="integer",
    topics=['SelfEnergy_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="SYMmetrization of SIGMA matrix elements",
    requires="[[optdriver]] in [4, 7]",
    commentdefault="The default value changed in Abinitv9 from 0 to 1",
    added_in_version="before_v9",
    text=r"""
This option activates the symmetrization of the self-energy matrix elements ([[symsigma]] = 1).
In this case the BZ integration defining the self-energy
matrix elements is reduced to an appropriate irreducible wedge defined
by the point group of the wave-vector k specified in the [[kptgw]] list.

The symmetrized expression leads to a considerable speedup of the run, especially
for high-symmetry k points e.g. $\Gamma$.
Unfortunately, this option is not yet compatible with self-consistent GW
calculations (see [[gwcalctyp]]).

The code constructs a symmetric invariant
for the diagonal matrix elements of the self-energy by averaging the self-energy matrix
elements within the degenerate subspace. Therefore, particular care has to be
taken in the presence of accidental degeneracies. Since calculations
performed with [[symsigma]] = 1 will not be able to remove the initial
accidental degeneracy. This is the reason why this option is not activated by default.
""",
),

Variable(
    abivarname="td_maxene",
    varset="dfpt",
    vartype="real",
    topics=['TDDFT_useful'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="Time-Dependent dft: MAXimal kohn-sham ENErgy difference",
    added_in_version="before_v9",
    text=r"""
The Matrix to be diagonalized in the Casida framework (see [[cite:Casida1995]])
is a NxN matrix, where, by default, N is the product of the number of occupied
states by the number of unoccupied states. The input variable [[td_maxene]]
allows one to diminish N: it selects only the pairs of occupied and unoccupied
states for which the Kohn-Sham energy difference is less than [[td_maxene]].
The default value 0.0 means that all pairs are taken into account.
See [[td_mexcit]] for an alternative way to decrease N.
""",
),

Variable(
    abivarname="td_mexcit",
    varset="dfpt",
    vartype="real",
    topics=['TDDFT_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Time-Dependent dft: Maximal number of EXCITations",
    added_in_version="before_v9",
    text=r"""
The Matrix to be diagonalized in the Casida framework (see [[cite:Casida1995]])
is a NxN matrix, where, by default, N is the product of the number of occupied
states by the number of unoccupied states. The input variable [[td_mexcit]]
allows one to diminish N: it selects the first [[td_mexcit]] pairs of occupied and
unoccupied states, ordered with respect to increasing Kohn-Sham energy difference.
However, when [[td_mexcit]] is zero, all pairs are allowed.
See [[td_maxene]] for an alternative way to decrease N.
""",
),

Variable(
    abivarname="tfkinfunc",
    varset="dev",
    vartype="integer",
    topics=['Recursion_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Thomas-Fermi KINetic energy FUNCtional",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
  * [[tfkinfunc]] = 1: Thomas-Fermi kinetic functional (explicit functional of the density) is used instead of Kohn-Sham kinetic
  energy functional (implicit functional of the density through Kohn-Sham wavefunctions). See [[cite:Perrot1979]].

  * [[tfkinfunc]] = 11: Thomas-Fermi-Weizsacker kinetic functional with Gradient Corrections is used.
The convergence of a calculation with this functional needs to be initialized
from a calculation without Gradient Correction. This is automatically done
with [[tfkinfunc]] = 11. For the initialization steps, the [[tfw_toldfe]]
criterion is used. When it is reached, then the Gradient Correction is added
and the SCF cycle continues.
Note: to obtain the convergence of a Molecular Dynamics simulation with TFW,
it is necessary to find the best set of preconditionning parameters
([[diemix]], [[diemac]], [[dielng]]) and the best value of [[npulayit]] (if
the default Pulay mixing is used).

  * [[tfkinfunc]] = 12: same as **tfkinfunc** =11, but without the initialization steps.
  Gradient correction is directly added.

  * [[tfkinfunc]] = 2: the Recursion Method is used in order to compute electronic density,
  entropy, Fermi energy and eigenvalues energy. This method computes the density
  without computing any orbital, is efficient at high temperature, with a efficient
  parallelization (almost perfect scalability).
  When that option is in use, the [[ecut]] input variable is no longer a convergence
  parameter; [[ngfft]] becomes the main convergence parameter: you should adapt ecut
  for the ngfft grid you need (it is not yet automatically computed).
  Other convergence parameter are for the energetic values: [[recnrec]], [[recptrott]], [[recnpath]].

Since the convergence of the self-consistent cycle is determined directly by
the convergence of the density: [[toldfe]], [[toldff]], [[tolrff]],
[[tolvrs]], [[tolwfr]] are not used, and are replaced by [[rectolden]]; the
energetic values, except for the fermi energy, are only computed during the
latest SFC cycle: the output file will show a jump of the total energy at the
end, but it is not because of a bad convergence behavior. Computational speed
can be improved by the use of [[recrcut]] and [[recgratio]]. The recursion
method has not been tested in the case of non cubic cell or with the use of
symmetries.
In the recursion method the following variables are set to: [[useylm]] = 1,
[[userec]] = 1.
""",
),

Variable(
    abivarname="tfw_toldfe",
    varset="dev",
    vartype="real",
    topics=['Recursion_useful'],
    dimensions="scalar",
    defaultval="1.0E-6 or [[toldfe]] is present",
    mnemonics="Thomas-Fermi-Weizsacker: TOLerance on the DiFference of total Energy, for initialization steps",
    characteristics=['[[ENERGY]]'],
    requires="[[tfkinfunc]] = 11",
    added_in_version="before_v9",
    text=r"""
This input variable has the same definition as [[toldfe]] and is only relevant
when [[tfkinfunc]] = 11.
It sets a tolerance for absolute differences of total energy that, reached
TWICE successively, will cause the initialization steps (without gradient
correction) to stop and the gradient correction to be added.
Can be specified in Ha (the default), Ry, eV or Kelvin, since it has the [[ENERGY]] characteristics.
""",
),

Variable(
    abivarname="tim1rev",
    varset="dfpt",
    vartype="integer",
    topics=['DFPT_expert', 'DFPT_internal'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="TIMe 1st order REVersal",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Allowed values are 0 or 1.

If tim1rev is equal to 1, the Sternheimer equation is solved simultaneously at
+q and -q perturbation wavevectors. The first order potential at -q is taken
to be equal to the Hermitian conjugate of the first order potential at +q.
The wavefunctions from both +q and -q are then combined to generate the first order density.
Relevant in the case of magnetic field perturbation (but will be relevant also in case of non-zero frequency DFPT, when implemented).
""",
),

Variable(
    abivarname="timopt",
    varset="gstate",
    vartype="integer",
    topics=['Control_expert'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[SEQUENTIAL]]': 1, 'defaultval': 0}),
    mnemonics="TIMing OPTion",
    characteristics=['[[NO_MULTI]]'],
    added_in_version="before_v9",
    text=r"""
This input variable allows one to modulate the use of the timing routines.

  * If 0  -->  as soon as possible, suppresses all calls to timing routines
  * If 1  -->  usual timing behaviour, with short analysis, appropriate for
sequential execution
  * If 2  -->  close to [[timopt]] = 1, except that the analysis routine does not time
the timer, appropriate for parallel execution.
  * If 3  -->  close to [[timopt]] = 1, except that the different parts of the lobpcg
routine are timed in detail.
  * If 4  -->  close to [[timopt]] = 1, except that the different parts of the lobpcg
routine are timed in detail. A different splitting of lobpcg than for
[[timopt]] = -3 is provided.
  * If -1  -->  a full analysis of timings is delivered
  * If -2  -->  a full analysis of timings is delivered, except timing the timer
  * If -3  -->  a full analysis of timings is delivered, including the detailed
timing of the different parts of the lobpcg routine. (this takes time, and is
discouraged for too small runs - the timing would take more time than the run
!). The timer is timed.
  * If -4  -->  a full analysis of timings is delivered, including the detailed
timing of the different parts of the lobpcg routine. A different splitting of
lobpcg than for [[timopt]] = -3 is provided (this takes time, and is discouraged
for too small runs - the timing would take more time than the run !). The
timer is timed. The sum of the independent parts is closer to 100% than for [[timopt]] = -3.
""",
),

Variable(
    abivarname="tl_nprccg",
    varset="gstate",
    vartype="integer",
    topics=['Wavelets_expert'],
    dimensions="scalar",
    defaultval=30,
    mnemonics="TaiL maximum Number of PReConditioner Conjugate Gradient iterations",
    added_in_version="before_v9",
    text=r"""
This variable is similar to [[wvl_nprccg]] but for the preconditioner
iterations during the tail corrections (see [[tl_radius]]).
""",
),

Variable(
    abivarname="tl_radius",
    varset="gstate",
    vartype="real",
    topics=['Wavelets_useful'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="TaiL expansion RADIUS",
    characteristics=['[[LENGTH]]'],
    added_in_version="before_v9",
    text=r"""
In the wavelet computation case, the linkage between the grid and the free
boundary conditions can be smoothed using an exponential decay. This means a
correction on the energy at the end on each wavefunction optimisation run. If
this parameter is set to zero, no tail computation is done. On the contrary,
put it to a positive value makes the tail correction available. The value
correspond to a length in atomic units being the spacial expansion with the
exponential decay around the grid.
""",
),

Variable(
    abivarname="tnons",
    varset="basic",
    vartype="real",
    topics=['crystal_useful'],
    dimensions=[3, '[[nsym]]'],
    mnemonics="Translation NON-Symmorphic vectors",
    added_in_version="before_v9",
    text=r"""
Gives the (nonsymmorphic) translation vectors associated with the symmetries
expressed in "[[symrel]]".
These may all be 0, or may be fractional (nonprimitive) translations expressed
relative to the real space primitive translations (so, using the "reduced"
system of coordinates, see "[[xred]]"). If all elements of the space group
leave 0 0 0 invariant, then these are all 0.
When the symmetry finder is used (see [[nsym]]), [[tnons]] is computed
automatically.

See also [[symafm]] for the complete description of the symmetry operation.
""",
),

Variable(
    abivarname="toldfe",
    varset="basic",
    vartype="real",
    topics=['SCFControl_basic'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="TOLerance on the DiFference of total Energy",
    characteristics=['[[ENERGY]]'],
    commentdefault="The default value implies that this stopping condition is ignored. For the SCF case, one and only one of the input tolerance criteria [[tolwfr]], [[toldff]], [[tolrff]], [[toldfe]] or [[tolvrs]] must differ from zero.",
    excludes="[[tolwfr]] or [[toldff]] or [[tolrff]] or [[tolvrs]]",
    added_in_version="before_v9",
    text=r"""
Sets a tolerance for absolute differences of total energy that, reached TWICE
successively, will cause one SCF cycle to stop (and ions to be moved).
Can be specified in Ha (the default), Ry, eV or Kelvin, since [[toldfe]] has
the [[ENERGY]] characteristics (1 Ha = 27.2113845 eV).
If set to zero, this stopping condition is ignored.
Effective only when SCF cycles are done ([[iscf]]>0).
Because of machine precision, it is not worth to try to obtain differences in
energy that are smaller than about 1.0d-12 of the total energy. To get
accurate stresses may be quite demanding.
When the geometry is optimized (relaxation of atomic positions or primitive
vectors), the use of [[toldfe]] is to be avoided. The use of [[toldff]] or
[[tolrff]] is by far preferable, in order to have a handle on the geometry
characteristics. When all forces vanish by symmetry (e.g. optimization of the
lattice parameters of a high-symmetry crystal), then place [[toldfe]] to
1.0d-12, or use (better) [[tolvrs]].
Since [[toldfe]], [[toldff]], [[tolrff]], [[tolvrs]] and [[tolwfr]] are aimed
at the same goal (causing the SCF cycle to stop), they are seen as a unique
input variable at reading. Hence, it is forbidden that two of these input
variables have non-zero values for the same dataset, or generically (for all
datasets). However, a non-zero value for one such variable for one dataset
will have precedence on the non-zero value for another input variable defined
generically.
""",
),

Variable(
    abivarname="toldff",
    varset="basic",
    vartype="real",
    topics=['SCFControl_basic', 'ForcesStresses_basic'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="TOLerance on the DiFference of Forces",
    commentdefault="The default value implies that this stopping condition is ignored. For the SCF case, one and only one of the input tolerance criteria [[tolwfr]], [[toldff]], [[tolrff]], [[toldfe]] or [[tolvrs]] must differ from zero.",
    excludes="[[tolwfr]] or [[toldfe]] or [[tolrff]] or [[tolvrs]]",
    added_in_version="before_v9",
    text=r"""
Sets a tolerance for differences of forces (in hartree/Bohr) that, reached
TWICE successively, will cause one SCF cycle to stop (and ions to be moved).
If set to zero, this stopping condition is ignored.
Effective only when SCF cycles are done ([[iscf]]>0). This tolerance applies
to any particular cartesian component of any atom, INCLUDING fixed ones. This
is to be used when trying to equilibrate a structure to its lowest energy
configuration (select [[ionmov]]), or in case of molecular dynamics ([[ionmov]] = 1)
A value ten times smaller than [[tolmxf]] is suggested (for example 5.0d-6
hartree/Bohr).
This stopping criterion is not allowed for RF calculations.
Since [[toldfe]], [[toldff]], [[tolrff]], [[tolvrs]] and [[tolwfr]] are aimed
at the same goal (causing the SCF cycle to stop), they are seen as a unique
input variable at reading. Hence, it is forbidden that two of these input
variables have non-zero values for the same dataset, or generically (for all
datasets). However, a non-zero value for one such variable for one dataset
will have precedence on the non-zero value for another input variable defined generically.
""",
),

Variable(
    abivarname="tolimg",
    varset="rlx",
    vartype="real",
    topics=['TransPath_basic'],
    dimensions="scalar",
    defaultval=5e-05,
    mnemonics="TOLerance on the mean total energy for IMaGes",
    characteristics=['[[ENERGY]]'],
    added_in_version="before_v9",
    text=r"""
Sets a maximal absolute energy tolerance (in hartree, averaged over dynamic
images) below which iterations on images (the one governed by the
[[ntimimage]] input variable) will stop.
This is to be used when trying to optimize a population of structures to their
lowest energy configuration, taking into account the particular algorithm
defined by [[imgmov]]
A value of about 5.0d-5 hartree or smaller is suggested (this corresponds to about 3.7d-7 eV).
No meaning for RF calculations.
""",
),

Variable(
    abivarname="tolmxde",
    varset="rlx",
    vartype="real",
    topics=['GeoOpt_basic'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="TOLerance on the MaXimal Difference in Energy",
    characteristics=['[[ENERGY]]'],
    added_in_version="before_v9",
    text=r"""
Sets a maximal difference in energy with respect to the two previous steps
below which BFGS structural relaxation iterations will stop.
A value of about 0.0005 eV/atom or smaller is suggested.
In order to use tolmxde, you should explicitly set tolmxf to 0.0.
No meaning for RF calculations.
""",
),

Variable(
    abivarname="tolmxf",
    varset="rlx",
    vartype="real",
    topics=['GeoOpt_basic'],
    dimensions="scalar",
    defaultval=5e-05,
    mnemonics="TOLerance on the MaXimal Force",
    added_in_version="before_v9",
    text=r"""
Sets a maximal absolute force tolerance (in hartree/Bohr) below which BFGS
structural relaxation iterations will stop.
Can also control tolerance on stresses, when [[optcell]] /=0, using the
conversion factor [[strfact]]. This tolerance applies to any particular
cartesian component of any atom, excluding fixed ones. See the parameter
[[ionmov]].
This is to be used when trying to equilibrate a structure to its lowest energy
configuration.
A value of about 5.0d-5 hartree/Bohr or smaller is suggested (this corresponds
to about 2.5d-3 eV/Angstrom).
No meaning for RF calculations.
""",
),

Variable(
    abivarname="tolrde",
    varset="dev",
    vartype="real",
    topics=['SCFControl_expert'],
    dimensions="scalar",
    defaultval=0.005,
    mnemonics="TOLerance on the Relative Difference of Eigenenergies",
    added_in_version="before_v9",
    text=r"""
Sets a tolerance for the ratio of differences of eigenenergies in the line
minimisation conjugate-gradient algorithm. It compares the decrease of the
eigenenergy due to the last line minimisation, with the one observed for the
first line minimisation. When the ratio is lower than [[tolrde]], the next
line minimisations are skipped.
The number of line minimisations is limited by [[nline]] anyhow.
This stopping criterion is present for both GS and RF calculations. In RF
calculations, [[tolrde]] is actually doubled before comparing with the above-
mentioned ratio, for historical reasons.
""",
),

Variable(
    abivarname="tolrff",
    varset="basic",
    vartype="real",
    topics=['SCFControl_basic', 'ForcesStresses_basic'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="TOLerance on the Relative diFference of Forces",
    commentdefault="The default value implies that this stopping condition is ignored. For the SCF case, one and only one of the input tolerance criteria [[tolwfr]], [[toldff]], [[tolrff]], [[toldfe]] or [[tolvrs]] must differ from zero.",
    excludes="[[tolwfr]] or [[toldfe]] or [[toldff]] or [[tolvrs]]'",
    added_in_version="before_v9",
    text=r"""
Sets a tolerance for the ratio of differences of forces (in hartree/Bohr) to
maximum force, that, reached TWICE successively, will cause one SCF cycle to
stop (and ions to be moved): diffor < tolrff * maxfor.
If set to zero, this stopping condition is ignored.
Effective only when SCF cycles are done ([[iscf]]>0). This tolerance applies
to any particular cartesian component of any atom, INCLUDING fixed ones. This
is to be used when trying to equilibrate a structure to its lowest energy
configuration (select [[ionmov]]), or in case of molecular dynamics ([[ionmov]] = 1)
A value of 0.02 is suggested.
This stopping criterion is not allowed for RF calculations.
Since [[toldfe]], [[toldff]], [[tolrff]], [[tolvrs]] and [[tolwfr]] are aimed
at the same goal (causing the SCF cycle to stop), they are seen as a unique
input variable at reading. Hence, it is forbidden that two of these input
variables have non-zero values for the same dataset, or generically (for all
datasets). However, a non-zero value for one such variable for one dataset
will have precedence on the non-zero value for another input variable defined
generically.
""",
),

Variable(
    abivarname="tolsym",
    varset="geo",
    vartype="real",
    topics=['crystal_useful'],
    dimensions="scalar",
    defaultval=1e-08,
    mnemonics="TOLERANCE for SYMmetries",
    added_in_version="before_v9",
    text=r"""
Gives the tolerance on the atomic positions (reduced coordinates), primitive
vectors, or magnetization, to be considered equivalent, thanks to symmetry
operations. This is used in the recognition of the set of symmetries of the
system, or the application of the symmetry operations to generate from a
reduced set of atoms, the full set of atoms. Note that a value larger than
0.01 is considered to be unacceptable, whatever the value of [[tolsym]]
(so, it is not worth to set [[tolsym]] bigger than 0.01).

Note: ABINIT needs the atomic positions to be symmetric to each others
within 1.e-8, irrespective of [[tolsym]].
So, if [[tolsym]] is set to a larger value than 1.e-8, then the
input atomic coordinates will be nevertheless automatically symmetrized by the symmetry
operations that will have been found.
"""
),

Variable(
    abivarname="tolvrs",
    varset="basic",
    vartype="real",
    topics=['SCFControl_basic'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="TOLerance on the potential V(r) ReSidual",
    commentdefault="The default value implies that this stopping condition is ignored. For the SCF case, one and only one of the input tolerance criteria [[tolwfr]], [[toldff]], [[tolrff]], [[toldfe]] or [[tolvrs]] must differ from zero.",
    excludes="[[tolwfr]] or [[toldfe]] or [[toldff]] or [[tolrff]]'",
    added_in_version="before_v9",
    text=r"""
Sets a tolerance for potential residual that, when reached, will cause one SCF
cycle to stop (and ions to be moved).
If set to zero, this stopping condition is ignored.
Effective only when SCF cycles are done ([[iscf]]>0).
To get accurate stresses may be quite demanding. For simple materials with
internal positions determined by symmetries, a value of [[tolvrs]] = 10^-12
empirically leads to a very approximate 10^-6 atomic unit accuracy for the
optimized lattice parameter.

Additional explanation: the residual of the potential is the difference
between the input potential and the output potential, when the latter is
obtained from the density determined from the eigenfunctions of the input
potential. When the self-consistency loop is achieved, both input and output
potentials must be equal, and the residual of the potential must be zero. The
tolerance on the potential residual is imposed by first subtracting the mean
of the residual of the potential (or the trace of the potential matrix, if the
system is spin-polarized), then summing the square of this function over all
FFT grid points. The result should be lower than [[tolvrs]].
Since [[toldfe]], [[toldff]], [[tolrff]], [[tolvrs]] and [[tolwfr]] are aimed
at the same goal (causing the SCF cycle to stop), they are seen as a unique
input variable at reading. Hence, it is forbidden that two of these input
variables have non-zero values for the same dataset, or generically (for all
datasets). However, a non-zero value for one such variable for one dataset
will have precedence on the non-zero value for another input variable defined
generically.
""",
),

Variable(
    abivarname="tolwfr",
    varset="basic",
    vartype="real",
    topics=['SCFControl_basic'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="TOLerance on WaveFunction squared Residual",
    commentdefault="The default value implies that this stopping condition is ignored. For the SCF case, one and only one of the input tolerance criteria [[tolwfr]], [[toldff]], [[tolrff]], [[toldfe]] or [[tolvrs]] must differ from zero.",
    excludes="[[toldfe]] or [[toldff]] or [[tolrff]] or [[tolvrs]]",
    added_in_version="before_v9",
    text=r"""
The signification of this tolerance depends on the basis set. In plane waves,
it gives a convergence tolerance for the largest squared "residual" (defined
below) for any given band. The squared residual is: < nk| (H-E)^2 |nk>,    E = < nk|H|nk >

$$
\langle \nk| (H - \enk)^2 |\nk \rangle, \,\text{with}\; \enk = \langle \nk|H|\nk \rangle
$$

which clearly is non-negative and goes to 0 as the iterations converge to an
eigenstate. With the squared residual expressed in Hartrees^2  (Hartrees
squared), the largest squared residual (called *residm* in the code) encountered over all
bands and k-points must be less than **tolwfr** for iterations to halt due to successful convergence.
Note that if [[iscf]] > 0, this criterion should be replaced by those based on
[[toldfe]] (preferred for [[ionmov]] == 0), [[toldff]] [[tolrff]] (preferred for
[[ionmov]] /= 0), or [[tolvrs]] (preferred for theoretical reasons!).
When **tolwfr** is 0.0, this criterion is ignored, and a finite value of
[[toldfe]], [[toldff]] or [[tolvrs]] must be specified. This also imposes a
restriction on taking an ion step; ion steps are not permitted unless the
largest squared residual is less than **tolwfr**, ensuring accurate forces.
To get accurate stresses may be quite demanding.
Note that the preparatory GS calculations before a RF calculations must be highly converged.
Typical values for these preparatory runs are **tolwfr** between 1.0d-16 and 1.0d-22.

Note that **tolwfr** is often used in the test cases, but this is **tolwfr**
purely for historical reasons: except when [[iscf]] < 0, other criteria should be used.

In the wavelet case (see [[usewvl]] = 1), this criterion is the favoured one.
It is based on the norm 2 of the gradient of the wavefunctions.
Typical values range from 5*10^-4  to 5*10^-5.

Since [[toldfe]], [[toldff]], [[tolrff]], [[tolvrs]] and **tolwfr** are aimed
at the same goal (causing the SCF cycle to stop), they are seen as a unique
input variable at reading. Hence, it is forbidden that two of these input
variables have non-zero values for the same dataset, or generically (for all
datasets). However, a non-zero value for one such variable for one dataset
will have precedence on the non-zero value for another input variable defined generically.
""",
),

Variable(
    abivarname="tphysel",
    varset="gstate",
    vartype="real",
    topics=['BandOcc_useful'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="Temperature (PHYSical) of the ELectrons",
    characteristics=['[[ENERGY]]'],
    added_in_version="before_v9",
    text=r"""
Gives, in Hartree, the physical temperature of the system, in case [[occopt]] = 4, 5, 6, or 7.

Can be specified in Ha (the default), Ry, eV or Kelvin, since [[tphysel]] has the
[[ENERGY]] characteristics (0.001 Ha = 27.2113845 meV = 315.773 Kelvin). One
has to specify an independent broadening [[tsmear]]. The combination of the
two parameters [[tphysel]] and [[tsmear]] is described in [[cite:Verstraete2002]].
Note that the signification of the entropy is modified with respect to the usual entropy.
The choice has been made to use [[tsmear]] as a prefactor of the entropy, to
define the entropy contribution to the free energy.
""",
),

Variable(
    abivarname="tsmear",
    varset="gstate",
    vartype="real",
    topics=['BandOcc_basic', 'STM_basic'],
    dimensions="scalar",
    defaultval=0.01,
    mnemonics="Temperature of SMEARing",
    characteristics=['[[ENERGY]]'],
    added_in_version="before_v9",
    text=r"""
Gives the broadening of occupation numbers [[occ]], in the metallic cases
([[occopt]] = 3, 4, 5, 6 and 7). Can be specified in Ha (the default), eV, Ry,
or Kelvin, since [[tsmear]] has the [[ENERGY]] characteristics (0.001 Ha =
27.2113845 meV = 315.773 Kelvin).
Default is 0.01 Ha. This should be OK using gaussian like smearings ([[occopt]] = 4,5,6,7)
for a free-electron metal like Al. For d-band metals, you may need to
use less.
Always check the convergence of the calculation with respect to this
parameter, and simultaneously, with respect to the sampling of k-points (see
[[nkpt]])
If [[occopt]] = 3, [[tsmear]] is the physical temperature, as the broadening is
based on Fermi-Dirac statistics. However, if [[occopt]] = 4, 5, 6, or 7, the
broadening is not based on Fermi-Dirac statistics, and [[tsmear]] is only a
convergence parameter. It is still possible to define a physical temperature,
thanks to the input variable [[tphysel]] (See also [[cite:Verstraete2002]]).
""",
),

Variable(
    abivarname="typat",
    varset="basic",
    vartype="integer",
    topics=['crystal_basic', 'AtomTypes_basic'],
    dimensions=ValueWithConditions({'[[natrd]]<[[natom]]': [3, '[[natrd]]'], 'defaultval': [3, '[[natom]]']}),
    defaultval=ValueWithConditions({'[[natom]] == 1': 1, 'defaultval': None}),
    mnemonics="TYPe of AToms",
    added_in_version="before_v9",
    text=r"""
Array giving an integer label to every atom in the unit cell to denote its
type.
The different types of atoms are constructed from the pseudopotential files.
There are at most [[ntypat]] types of atoms.
As an example, for BaTiO3, where the pseudopotential for Ba is number 1, the
one of Ti is number 2, and the one of O is number 3, the actual value of the
[[typat]] array might be:

      typat 1 2 3 3 3

The array [[typat]] has to agree with the actual locations of atoms given in
[[xred]], [[xcart]] or [[xangst]], and the input of pseudopotentials has to
be ordered to agree with the atoms identified in [[typat]].
The nuclear charge of the elements, given by the array [[znucl]], also must
agree with the type of atoms designated in "[[typat]]".
The array [[typat]] is not constrained to be increasing, so

      typat 3 3 2 3 1

is permitted. An internal
representation of the list of atoms, deep in the code (array atindx), groups
the atoms of same type together. This should be transparent to the user, while
keeping efficiency.
""",
),

Variable(
    abivarname="ucrpa",
    varset="gw",
    vartype="integer",
    topics=['CRPA_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="calculation of the screened interaction U with the Constrained RPA method",
    requires="[[nspinor]] == 1",
    added_in_version="before_v9",
    text=r"""
When equal to one or two, this variable allows one to calculate U with
the cRPA method. An explicit test is shown in automatic tests
[[test:v7_23]], [[test:v7_24]], [[test:v7_25]], [[test:v7_68]], and [[test:v7_69]].
The present implementation is parallelized (as for usual GW
calculations), use symmetry over k points only for calculations involving one
correlated atom, and can be use when correlated bands are entangled or not.
The constrained calculation of the polarisability can be done by eliminating
transition betweens correlated bands (and not orbitals) with the variable
[[ucrpa_bands]].

For [[ucrpa]] = 1, two solutions are possible. The first one is to specify
(with the variable [[ucrpa_bands]]) the bands to exclude from the
polarisability calculation. The second solution is to provide an energy window
(with the variable [[ucrpa_window]]). The electronic transitions inside this
window will not be taken into account in the polarisability calculation.

For [[ucrpa]] = 2, the ucrpa_bands should be equal to the [[dmftbandi]] and
[[dmftbandf]] values, and the polarisability of the correlated subspace is
constructed with a band and k point dependent weight.

The implementation is restricted to the case of [[nspinor]] = 1 (collinear case).

A short presentation of the method and some aspect of the implementation can
be found in Sec. II and Appendix A of [[cite:Amadon2014]].
""",
),

Variable(
    abivarname="ucrpa_bands",
    varset="gw",
    vartype="integer",
    topics=['CRPA_basic'],
    dimensions=[2],
    defaultval=[-1, -1],
    mnemonics="For the calculation of U with the Constrained RPA method, gives correlated BANDS",
    commentdefault="That is, the default includes no band.",
    added_in_version="before_v9",
    text=r"""
Gives the first and last correlated bands for the cRPA calculation of the polarisability.
""",
),

Variable(
    abivarname="ucrpa_window",
    varset="gw",
    vartype="real",
    topics=['CRPA_basic'],
    dimensions=[2],
    defaultval=[-1, -1],
    mnemonics="For the calculation of U with the Constrained RPA method, gives energy WINDOW",
    commentdefault="That is, the energy window is empty by default.",
    added_in_version="before_v9",
    text=r"""
Specify a window of energy for the cRPA calculation of the polarisability. The
transition inside this window will not be taken into account in the
constrained polarisabilty calculations.

The lower bound and the upper bound energies must be specified (two real
numbers) with respect to the position of the Fermi level.
""",
),

Variable(
    abivarname="udtset",
    varset="basic",
    vartype="integer",
    topics=['multidtset_basic'],
    dimensions=[2],
    mnemonics="Upper limit on DaTa SETs",
    commentdefault="It is not used when it is not defined",
    added_in_version="before_v9",
    text=r"""
Used to define the set of indices in the multi-data set mode, when a double
loop is needed (see later).
The values of [[udtset]](1) must be between 1 and 999, the values of
[[udtset]](2) must be between 1 and 9, and their product must be equal to
[[ndtset]].
The values of [[jdtset]] are obtained by looping on the two indices defined by
[[udtset]](1) and [[udtset]](2) as follows:

```fortran
do i1=1,intarr(1)
 do i2=1,intarr(2)
  idtset=idtset+1
  dtsets(idtset)%jdtset=i1*10+i2
 end do
end do
```

So, [[udtset]](2) sets the largest value for the unity digit, that varies
between 1 and [[udtset]](2).
If [[udtset]] is used, the input variable [[jdtset]] cannot be used.
""",
),

Variable(
    abivarname="upawu",
    varset="paw",
    vartype="real",
    topics=['DFT+U_compulsory'],
    dimensions=['[[ntypat]]'],
    defaultval=MultipleValue(number=None, value=0),
    mnemonics="value of U for PAW+U",
    characteristics=['[[ENERGY]]'],
    requires="[[usepaw]] == 1 and [[usepawu]] == 1",
    added_in_version="before_v9",
    text=r"""
Gives the value of the screened coulomb interaction between correlated
electrons corresponding to [[lpawu]] for each species.
In the case where [[lpawu]] =-1, the value is not used.
In the case of a GW calculation, the U interaction defined by [[upawu]] will
be REMOVED from the self energy. In particular, for G0 W0 calculations
(perturbative calculations), the energy eigenvalues obtained after an
underlying DFT+U calculation will be

$$E_{GW} = E_{DFT+U} + < \phi | Self-energy - U | \phi>$$

Actually, in order to perform a GW @ DFT+U calculation, one should define the
same value of U in the self-energy calculation, than the one defined in the
DFT calculation. The easiest is actually to define the value of U for the
whole set of calculations (for the different datasets), including the
screening, even if the U value does not play explicitly a role in the
computation of the latter (well, the input wavefunctions will be different
anyhow).
It is possible to perform calculations of the type GW+$U^{'}$ @ DFT+U, so
keeping a $U^{'}$ interaction (usually smaller than the initial U) in the GW
calculation.
This value will be subtracted in the GW correction calculation, as outlined
above.
Explicitly, in order to do a calculation of a material with a DFT U value of
7.5 eV, followed by a GW calculation where there is a residual U value of 2
eV, one has to define:

      upawu1   7.5 eV   ! This is for the DFT calculation
    ...
    optdriver4  4
    upawu4   5.5 eV   ! This is for the screening calculation
""",
),

Variable(
    abivarname="use_gemm_nonlop",
    varset="dev",
    vartype="integer",
    topics=['parallelism_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="USE the GEMM routine for the application of the NON-Local OPerator",
    characteristics=['[[DEVELOP]]'],
    commentdefault="because it is not usually worth using it unless bandpp is large and it requires additional memory",
    added_in_version="before_v9",
    text=r"""
This keyword tells abinit to use a BLAS routine to speed up the computation of
the non-local operator. This requires the pre-computation of a large matrix,
and has a significant memory overhead. In exchange, it provides improved
performance when used on several bands at once (Chebyshev or LOBPCG algorithm
with [[bandpp]]

The memory overhead is proportional to the number of atoms, the number of
plane waves, and the number of projectors per atom. It can be mitigated by
distributing the array with [[npfft]]

The performance depends crucially on having a good BLAS installed. Provided
the BLAS supports OpenMP, this option also yields very good scaling for the
nonlocal operator.
""",
),

Variable(
    abivarname="use_gpu_cuda",
    varset="paral",
    vartype="integer",
    topics=['parallelism_expert'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[optdriver]] == 0 and [[CUDA]]': 1, 'defaultval': 0}),
    mnemonics="activate USE of GPU accelerators with CUDA (nvidia)",
    added_in_version="before_v9",
    text=r"""
Only available if ABINIT executable has been compiled with cuda nvcc compiler.
This parameter activates the use of NVidia graphic accelerators (GPU) if present.
If [[use_gpu_cuda]] = 1, some parts of the computation are transmitted to the GPUs.
If [[use_gpu_cuda]] = 0, no computation is done on GPUs, even if present.

Note that, while running ABINIT on GPUs, it is recommended to use MAGMA
external library (i.e. Lapack on GPUs). The latter is activated during
compilation stage (see "configure" step of ABINIT compilation process). If
MAGMA is not used, ABINIT performances on GPUs can be poor.
""",
),

Variable(
    abivarname="use_nonscf_gkk",
    varset="dev",
    vartype="integer",
    topics=['ElPhonInt_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="USE NON-SCF calculation of GKK matrix elements (electron phonon)",
    characteristics=['[[DEVELOP]]'],
    commentdefault="Default is 0 for the moment. Do not use non-scf method.",
    added_in_version="before_v9",
    text=r"""
When this flag is activated during a phonon calculation with abinit, all of
the perturbations are cycled through, but only the symmetry-irreducible ones
are calculated self-consistently. For the others the perturbed density is
rotated by the appropriate symop and the gkk matrix elements are calculated
non-self-consistently. As they do not depend on the perturbed wave functions,
they are correct from the first iteration, and nstep is set to 1 for those
perturbations. Note that the resulting 1DEN files are simply the
rotate/symmetric ones and that the resulting 1WF files are garbage (completely
unconverged) except the matrix elements in the header (equivalent to GKK
files, but please use the latter much smaller files for el-ph calculations).
The new default behavior with [[use_nonscf_gkk]] = 1 should be transparent for
the user, with the same output files but a much quicker execution.

Caveat: Note that very tight convergence of ground state and phonon
calculations is necessary to get good GKK matrix elements! [[tolwfr]] = 1.e-24
or so is recommended everywhere. There may be problems using use_nonscf_gkk =
1 with non-symmorphic symmetries - please check (at least) that lifetimes for
phonons go to 0 for acoustic modes at Gamma.
""",
),

Variable(
    abivarname="use_slk",
    varset="paral",
    vartype="integer",
    topics=['parallelism_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="USE ScaLapacK",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
If set to 1, enable the use of ScaLapack within LOBPCG.
""",
),

Variable(
    abivarname="usedmatpu",
    varset="paw",
    vartype="integer",
    topics=['DFT+U_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="USE of an initial Density MATrix in Paw+U",
    requires="[[usepaw]] == 1 and [[usepawu]] == 1",
    added_in_version="before_v9",
    text=r"""
When [[usedmatpu]]/=0, an initial density matrix (given by [[dmatpawu]]
keyword) is used and kept fixed during the first ABS([[usedmatpu]]) SCF steps.
This starting value of the density matrix can be useful to find the correct
ground state. Within LDA+U formalism, finding the minimal energy of the system
is tricky; thus it is advised to test several values of the initial density
matrix.
Note also that the density matrix has to respect some symmetry rules
determined by the space group. If the symmetry is not respected in the input,
the matrix is however automatically symmetrised.

The sign of [[usedmatpu]] has influence only when [[ionmov]] /= 0 (dynamics or relaxation):

- When [[usedmatpu]]>0, the density matrix is kept constant only at first ionic step
- When [[usedmatpu]]<0, the density matrix is kept constant at each ionic step
""",
),

Variable(
    abivarname="usedmft",
    varset="dev",
    vartype="integer",
    topics=['DMFT_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="USE Dynamical Mean Field Theory",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
If set to 1, enable the use of DFT+DMFT, see in particular the important
variables [[dmft_solv]], [[dmftbandi]], [[dmftbandf]], [[dmft_nwli]],
[[dmft_nwlo]], [[dmft_tollc]], [[dmft_tolfreq]], and [[dmft_iter]].

The current implementation uses Wannier functions obtained from
[[ cite:Amadon2008 | projected local orbitals ]] as
correlated orbitals (see [[dmftbandi]] and [[dmftbandf]] input variables to define them).

The Green functions are computed on a mesh of linear Matsubara frequencies.
However, most of the code uses logarithmic Matsubara grid to lower the
computational cost. Both [[dmft_nwli]] and [[dmft_nwlo]] are thus convergence parameters.

DMFT is currently available for collinear ([[nspinor]] = 1) polarized or
unpolarized calculations ([[nspden]] = [[nsppol]] = 2 or [[nspden]] = [[nsppol]] = 1)
and for non collinear calculations ([[nspinor]] = 2,[[nspden]] = 4,[[nsppol]] = 1).
However it is not yet available for collinear antiferromagnetic calculations
([[nspden]] = 2,[[nsppol]] = 1) and non collinear non magnetic calculations
([[nspden]] = 1, [[nsppol]] = 1,[[nspinor]] = 2). CTQMC calculations
([[dmft_solv]] = 5) are not yet possible if [[nspinor]] = 2.

Only static calculations without relaxation or dynamics are possible (forces
and stress are not computed in the scheme: so the computed values should NOT
be trusted).

When correlated density matrices are diagonal, all values of [[upawu]] and
[[jpawu]] are possible. If the correlated density matrices are non diagonal,
only [[jpawu]] = 0 is implemented.

Relevant direct output quantities from converged DMFT calculations are total
energy and occupation of correlated orbitals. For Hubbard I calculation
([[dmft_solv]] = 2), total and partial spectral functions can be obtained with
prtdos=1 and can be found in files OUTSpFunc* (where OUT is the root for
output files). For CTQMC calculations ([[dmft_solv]] = 5), imaginary time
impurity Green function are output of the calculations and can be used to
produce spectral function using an external Maximum Entropy Code.

A typical DFT+DMFT calculation involves two runs. First, a DFT calculation is
fully converged (even unoccupied wavefunctions have to be converged). Then,
the DFT+DMFT calculation is started using DFT wavefunctions or density files.
As DFT+DMFT calculations (with CTQMC) are computationally expensive, it is
convenient to use prtden=-1, to write DEN file at each DFT iteration, in order
to be able to restart the calculation easily.

For details of the implementation with Wannier functions see [[cite:Amadon2008]],
for self-consistency and Hubbard I
implementation see [[cite:Amadon2012]]. If [[usedmft]] = 1 and [[nbandkss]]/=0, then, the DFT+DMFT
calculation is not done and only projections are computed at the end of the
calculation. They can be used by an external code or used to compute the
screened interaction (see variable [[ucrpa]]).
""",
),

Variable(
    abivarname="useexexch",
    varset="paw",
    vartype="integer",
    topics=['xc_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="USE of EXact EXCHange",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
When [[useexexch]] = 1, the hybrid functional PBE0 is used in PAW, inside PAW
spheres only, and only for correlated orbitals given by [[lexexch]]. To change
the ratio of exact exchange, see also [[exchmix]].
""",
),

Variable(
    abivarname="usefock",
    varset="internal",
    vartype="integer",
    topics=['Hybrids_internal'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="USE FOCK exact exchange",
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
This internal variable is automatically set to 1 when the value of [[ixc]]
refers to an Hartree-Fock calculation or hybrid functionals.

  * 0  -->  No use of exact exchange.
  * 1  -->  exact exchange is required for the calculation.
""",
),

Variable(
    abivarname="usekden",
    varset="gstate",
    vartype="integer",
    topics=['xc_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="USE Kinetic energy DENsity",
    added_in_version="before_v9",
    text=r"""
If [[usekden]] = 1 the kinetic energy density will be computed during
the self-consistent loop, in a way similar to the computation of the density.
This is needed if a meta-GGA is to be used as XC functional. By default
([[usekden]] = 0), the kinetic energy density is not computed during the self-
consistent loop.
""",
),

Variable(
    abivarname="usepaw",
    varset="internal",
    vartype="integer",
    topics=['PAW_internal'],
    dimensions="scalar",
    defaultval="[[AUTO_FROM_PSP]]",
    mnemonics="USE Projector Augmented Waves method",
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
This variable is determined by the pseudopotentials files. PAW calculations
(see [[tutorial:paw1]]) can only be performed with PAW atomic data input files,
while pseudopotential calculations are performed in ABINIT with norm-conserving
pseudopotential input files. Most functionalities in ABINIT are
available with either type of calculation.
""",
),

Variable(
    abivarname="usepawu",
    varset="paw",
    vartype="integer",
    topics=['DFT+U_compulsory', 'PAW_useful', 'GW_useful', 'SelfEnergy_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="USE PAW+U (spherical part)",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""
Must be non-zero if a DFT+U calculation is done, or if a GW calculation
following a DFT+U calculation is done (important!).

  * If set to 0, the LDA+U method is not used.

  * If set to 1, 2 or 4, the LDA+U method (cf [[cite:Anisimov1991a]]) is used.
The full rotationally invariant formulation is used (see Eq. (3) of [[cite:Liechtenstein1995]]) for the interaction term of the energy.
Three choices are allowed concerning the double counting term:

    * If [[usepawu]] = 1, the Full Localized Limit (FLL) (or Atomic limit) double counting is used (cf Eq. (4) of [[cite:Liechtenstein1995]] or Eq. (8) of [[cite:Czyzyk1994]]).

    * If [[usepawu]] = 2, the Around Mean Field (AMF) double counting is used (cf Eq. (7) of [[cite:Czyzyk1994]]). Not valid if nspinor=2.

    * If [[usepawu]] = 4, the FLL double counting is used. However, and in comparison to usepaw=1, the calculation is done without
    polarization in the exchange correlation functional (cf [[cite:Park2015]] and [[cite:Chen2016a]]). In this case, one must use [[iscf]]<10.

If LDA+U is activated ([[usepawu]] = 1 or 2), the [[lpawu]], [[upawu]] and
[[jpawu]] input variables are read.
The implementation is done inside PAW augmentation regions only (cf [[cite:Bengone2000]]).
The initial density matrix can be given in the input file (see [[usedmatpu]]).
The expression of the density matrix is chosen thanks to [[dmatpuopt]].
In the case of a GW calculation on top of a DFT+U, the absence of definition
of a U value in the self-energy will LEAVE the underlying U from the DFT
calculation. Thus, the code will actually do a GW+U @ DFT+U calculation. Note
that the screening calculation will not be affected by the presence/absence of
a U value.
Actually, in order to perform a GW @ DFT+U calculation, one should define the
same value of U in the self-energy calculation, than the one defined in the
DFT calculation. The code will know that the interaction corresponding to that
value has to be SUBTRACTED inside the self-energy. The easiest is actually to
define the presence of U for the whole set of calculations (for the different
datasets), including the screening, even if the U value does not play
explicitly a role in the computation of the latter (well, the input
wavefunctions will be different anyhow).
It is possible to perform calculations of the type GW+$U^{'}$ @ DFT+U, so
keeping a smaller U interaction in the GW calculation, by subtracting a
smaller U than the one used in the DFT calculation. See the description of the
[[upawu]] input variable.


Suggested acknowledgment:[[cite:Amadon2008a]].

""",
),

Variable(
    abivarname="usepead",
    varset="dfpt",
    vartype="integer",
    topics=['nonlinear_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="USE of PEAD formalism",
    requires="[[optdriver]] == 5 (non-linear response computations)",
    added_in_version="before_v9",
    text=r"""
Determine which non-linear implementation is used. If [[usepead]]=1, the Perturbation
Expansion After Discretization formalism is used, as in [[cite:Veithen2005]].
In that method, the electric field is treated numerically, i.e the k-space
gradient operator appearing in the expression of the electric field potential
is discretized (see Eq.7 and 10 of [[cite:Veithen2005]]).
If [[usepead]]=0, the electric field is treated analytically, leading to a better
k-points convergence. Furthermore, the current implementation is compatible with PAW
pseudopentials, while [[usepead]]=1 is not. The drawback of the analytical method
is one has to solve a second order Sternheimer equation before actually computing
third derivatives of the energy, using [[rf2_dkdk]] and [[rf2_dkde]].
This is not the most time-consumming part though.
Look at the inputs of related tests in the testsuite to see examples of the workflow.
""",
),

Variable(
    abivarname="usepotzero",
    varset="dev",
    vartype="integer",
    topics=['Coulomb_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="USE POTential ZERO",
    added_in_version="before_v9",
    text=r"""
Fix the convention for the choice of the average value of the Hartree potential, as described in [[cite:Bruneval2014]].
  * [[usepotzero]] = 0, the usual convention: the smooth potential is set to zero average value.
  * [[usepotzero]] = 1, the new convention: the all-electron physical potential is set to zero average value.
  * [[usepotzero]] = 2, the PWscf convention: the potential of equivalent point charges is set to
  zero average value (convention also valid for NC pseudopotentials).
""",
),

Variable(
    abivarname="userec",
    varset="internal",
    vartype="integer",
    topics=['Recursion_internal'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="USE RECursion",
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
This internal variable is set to 1 when the recursion method is activated (see [[tfkinfunc]]).
""",
),

Variable(
    abivarname="useria",
    varset="dev",
    vartype="integer",
    topics=['Dev_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="USER Integer variable A",
    added_in_version="before_v9",
    text=r"""
These are user-definable integers which the user may input and then utilize in
subroutines of his/her own design. They are not used in the official versions
of the ABINIT code, and should ease independent developments (hopefully
integrated in the official version afterwards).
Internally, they are available in the dtset structured datatype, e.g. `dtset%useria`.
""",
),

Variable(
    abivarname="userib",
    varset="dev",
    vartype="integer",
    topics=['Dev_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="USER Integer variable B",
    added_in_version="before_v9",
    text=r"""
These are user-definable integers which the user may input and then utilize in
subroutines of his/her own design. They are not used in the official versions
of the ABINIT code, and should ease independent developments (hopefully
integrated in the official version afterwards).
Internally, they are available in the dtset structured datatype, e.g. `dtset%useria`.
""",
),

Variable(
    abivarname="useric",
    varset="dev",
    vartype="integer",
    topics=['Dev_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="USER Integer variable C",
    added_in_version="before_v9",
    text=r"""
These are user-definable integers which the user may input and then utilize in
subroutines of his/her own design. They are not used in the official versions
of the ABINIT code, and should ease independent developments (hopefully
integrated in the official version afterwards).
Internally, they are available in the dtset structured datatype, e.g. `dtset%useria`.
""",
),

Variable(
    abivarname="userid",
    varset="dev",
    vartype="integer",
    topics=['Dev_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="USER Integer variable D",
    added_in_version="before_v9",
    text=r"""
These are user-definable integers which the user may input and then utilize in
subroutines of his/her own design. They are not used in the official versions
of the ABINIT code, and should ease independent developments (hopefully
integrated in the official version afterwards).
Internally, they are available in the dtset structured datatype, e.g. `dtset%useria`.
""",
),

Variable(
    abivarname="userie",
    varset="dev",
    vartype="integer",
    topics=['Dev_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="USER Integer variable E",
    added_in_version="before_v9",
    text=r"""
These are user-definable integers which the user may input and then utilize in
subroutines of his/her own design. They are not used in the official versions
of the ABINIT code, and should ease independent developments (hopefully
integrated in the official version afterwards).
Internally, they are available in the dtset structured datatype, e.g. `dtset%useria`.
""",
),

Variable(
    abivarname="userra",
    varset="dev",
    vartype="real",
    topics=['Dev_expert'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="USER Real variable A",
    added_in_version="before_v9",
    text=r"""
These are user-definable with the same purpose as [[useria]] and cie.
""",
),

Variable(
    abivarname="userrb",
    varset="dev",
    vartype="real",
    topics=['Dev_expert'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="USER Real variable B",
    added_in_version="before_v9",
    text=r"""
These are user-definable with the same purpose as [[useria]] and cie.
""",
),

Variable(
    abivarname="userrc",
    varset="dev",
    vartype="real",
    topics=['Dev_expert'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="USER Real variable C",
    added_in_version="before_v9",
    text=r"""
These are user-definable with the same purpose as [[useria]] and cie.
""",
),

Variable(
    abivarname="userrd",
    varset="dev",
    vartype="real",
    topics=['Dev_expert'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="USER Real variable D",
    added_in_version="before_v9",
    text=r"""
These are user-definable with the same purpose as [[useria]] and cie.
""",
),

Variable(
    abivarname="userre",
    varset="dev",
    vartype="real",
    topics=['Dev_expert'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="USER Real variable E",
    added_in_version="before_v9",
    text=r"""
These are user-definable with the same purpose as [[useria]] and cie.
""",
),

Variable(
    abivarname="usewvl",
    varset="basic",
    vartype="integer",
    topics=['Wavelets_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Use WaVeLet basis set",
    commentdefault="use plane-wave basis set",
    added_in_version="before_v9",
    text=r"""
Used to define if the calculation is done on a wavelet basis set or not.
The values of [[usewvl]] must be 0 or 1. Putting [[usewvl]] to 1, makes
[[icoulomb]] mandatory to 1. The number of band ([[nband]]) must be set
manually to the strict number need for an isolator system ( _i.e._ number of
electron over two). The cut-off is not relevant in the wavelet case, use
[[wvl_hgrid]] instead.
In wavelet case, the system must be isolated systems (molecules or clusters).
All geometry optimization are available (see [[ionmov]], especially the
geometry optimisation and the molecular dynamics).
The spin computation is not currently possible with wavelets and metallic
systems may be slow to converge.
""",
),

Variable(
    abivarname="usexcnhat",
    varset="paw",
    vartype="integer",
    topics=['PAW_useful'],
    dimensions="scalar",
    defaultval=-1,
    mnemonics="USE eXchange-Correlation with NHAT (compensation charge density)",
    requires="[[usepaw]] == 1",
    added_in_version="before_v9",
    text=r"""

This flag determines how the exchange-correlation terms are computed for the
pseudo-density.
When **usexcnhat** = 0, the exchange-correlation potential does not include the
compensation charge density, i.e. $V_{xc}=V_{xc}(\tilde{n}_{core} + \tilde{n}_{valence})$.
When **usexcnhat** = 1, the exchange-correlation potential includes the compensation
charge density, i.e. $V_{xc}=V_{xc}(\tilde{n}_{core} + \tilde{n}_{valence}+\hat{n})$.
When **usexcnhat** = -1,the value of **usexcnhat** is determined from the
reading of the PAW dataset file (pseudopotential file). When PAW datasets with
different treatment of $V_{xc}$ are used in the same run, the code stops.
""",
),

Variable(
    abivarname="useylm",
    varset="dev",
    vartype="integer",
    topics=['TuningSpeed_expert'],
    dimensions="scalar",
    defaultval=ValueWithConditions({'[[tfkinfunc]] == 1': 1, '[[usepaw]] == 1': 1, 'defaultval': 0}),
    mnemonics="USE YLM (the spherical harmonics)",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
When this flag is activated, the non-local operator is applied using an
algorithm based on spherical harmonics. Non-local projectors are used with
their usual form:

$$P_{lmn}(r) = Y_{lm}(r) p_{ln}(r)$$

When [[useylm]] = 0, the sum over $Y_{lm}$ can be reduced to a Legendre polynomial form.
""",
),

Variable(
    abivarname="vaclst",
    varset="geo",
    vartype="integer",
    topics=['AtomManipulator_useful'],
    dimensions=['[[vacnum]]'],
    mnemonics="VACancies LiST",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Gives the identification number(s) of atoms to be subtracted from the set of
atoms that are obtained after having rotated, translated and repeated the objects.
Useful to created vacancies.
""",
),

Variable(
    abivarname="vacnum",
    varset="geo",
    vartype="integer",
    topics=['AtomManipulator_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="VACancies NUMber",
    added_in_version="before_v9",
    text=r"""
Gives the number of atoms to be subtracted from the list of atoms after the
rotations, translations and repetitions have been done. The list of these
atoms is contained in [[vaclst]].
""",
),

Variable(
    abivarname="vacuum",
    varset="gstate",
    vartype="integer",
    topics=['k-points_expert'],
    dimensions=[3],
    mnemonics="VACUUM identification",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Establishes the presence (if [[vacuum]] = 1) or absence (if [[vacuum]] = 0) of a vacuum layer, along the
three possible directions normal to the primitive axes.

This information might be used to generate k-point grids, if [[kptopt]] = 0 and
neither [[ngkpt]] nor [[kptrlatt]] are defined (see explanations with the
input variable [[prtkpt]]).
It will allow to select a zero-, one-, two- or three-dimensional grid of k
points. The coordinate of the k points along vacuum directions is
automatically set to zero.

If [[vacuum]] is not defined, the input variable [[vacwidth]] will be used to
determine automatically whether the distance between atoms is sufficient to
have the presence or absence of vacuum.
""",
),

Variable(
    abivarname="vacwidth",
    varset="gstate",
    vartype="real",
    topics=['k-points_expert'],
    dimensions="scalar",
    defaultval=10.0,
    mnemonics="VACuum WIDTH",
    characteristics=['[[INPUT_ONLY]]', '[[LENGTH]]'],
    added_in_version="before_v9",
    text=r"""
Give a minimum "projected" distance between atoms to be found in order to
declare that there is some [[vacuum]] present for each of the three
directions. By default, given in Bohr atomic units (1 Bohr=0.5291772108
Angstroms), although Angstrom can be specified, if preferred, since
[[vacwidth]] has the [[LENGTH]] characteristics.
The precise requirement is that a slab of width [[vacwidth]], delimited by two
planes of constant reduced coordinates in the investigated direction, must be empty of atoms.
""",
),

Variable(
    abivarname="vcutgeo",
    varset="gw",
    vartype="real",
    topics=['GWls_compulsory', 'Susceptibility_basic', 'SelfEnergy_basic'],
    dimensions=[3],
    defaultval=MultipleValue(number=3, value=0.0),
    mnemonics="V (potential) CUT-off GEOmetry",
    requires="[[icutcoul]] in [1,2]",
    added_in_version="before_v9",
    text=r"""
[[vcutgeo]] is used in conjunction with [[icutcoul]] to specify the geometry
used to truncate the Coulomb interaction, as well as the particular approach
to be used. It has a meaning only for the cylindrical symmetry
([[icutcoul]] = 1) or in the case of surfaces ([[icutcoul]] = 2). For each
geometry, two different definitions of the cutoff region are available (see
Phys. Rev. B 73, 233103 and Phys. Rev. B 73, 205119 for a complete description
of the methods)

In the method of Ismail-Beigi [[cite:Ismail-Beigi2006]], the cutoff region is given by the
Wigner-Seitz cell centered on the axis of the cylinder. The cutoff region is
thus automatically defined by the unit cell and there is no need to specify
When [[rcut]].

To define a cylinder along the z-axis use the following lines:
```
icutcoul 1
vcutgeo  0 0 1
```

Please note that the method of Ismail-Beigi is implemented only in the case if an
orthorhombic Bravais lattice. For hexagonal lattices, one has to use the method
of Rozzi [[cite:Rozzi2006]]. In this case, the interaction is truncated
in a finite cylinder. Contrarily to the first approach, here one has to
specify both the radius of the cylinder with [[rcut]] as well as the length of
the cylinder along the periodic dimension that should always be smaller than
the extension of the Born von Karman box. The length of the cylinder is given
in terms of the fraction of the primitive vector along the periodic direction.

For example, in order to define a finite cylinder along z of radius 2.5 Bohr
and length 3*R3,
```
icutcoul 1
vcutgeo  0 0 -3.0 # note the minus sign
rcut     2.5
```

For surface calculations ([[icutcoul]] = 2), [[vcutgeo]] is used to define the
two periodic directions defining the surface. Also in this case two different
techniques are available. In the method of Ismail-Beigi, the (positive) non-zero
components of vcutgeo define the periodic directions of the infinite surface.
The interaction is truncated within a slab of width L where L is the length of
the primitive vector of the lattice along the non-periodic dimension. For
example:
```
icutcoul 2
vcutgeo  1 1 0
```

It is also possible to define a finite
surface by employing negative values. For example:
```
icutcoul 2
vcutgeo -3 -2 0
```
**Definition to be added**
""",
),

Variable(
    abivarname="vdw_df_acutmin",
    varset="vdw",
    vartype="real",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=10,
    mnemonics="vdW-DF MINimum Angular CUT-off",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, to build angular meshes for the vdW-DF kernel.
""",
),

Variable(
    abivarname="vdw_df_aratio",
    varset="vdw",
    vartype="real",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=30,
    mnemonics="""vdW-DF Angle RATIO between the highest and lowest angles.""",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, to build angular meshes for the vdW-DF kernel.
""",
),

Variable(
    abivarname="vdw_df_damax",
    varset="vdw",
    vartype="real",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=0.5,
    mnemonics="vdW-DF Delta for Angles, MAXimum",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, to build angular meshes for the vdW-DF kernel.
""",
),

Variable(
    abivarname="vdw_df_damin",
    varset="vdw",
    vartype="real",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=0.01,
    mnemonics="vdW-DF Delta for Angles, MINimum",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, to build angular meshes for the vdW-DF kernel.
""",
),

Variable(
    abivarname="vdw_df_dcut",
    varset="vdw",
    vartype="real",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=30,
    mnemonics="vdW-DF D-mesh CUT-off",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, to build the vdW-DF kernel.
""",
),

Variable(
    abivarname="vdw_df_dratio",
    varset="vdw",
    vartype="real",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=20,
    mnemonics="""vdW-DF, between the highest and
lowest D, RATIO.""",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, to build the vdW-DF kernel.
""",
),

Variable(
    abivarname="vdw_df_dsoft",
    varset="vdw",
    vartype="real",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=1.0,
    mnemonics="vdW-DF Distance for SOFTening.",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, to build the vdW-DF kernel.
""",
),

Variable(
    abivarname="vdw_df_gcut",
    varset="vdw",
    vartype="real",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=5,
    mnemonics="vdW-DF G-space CUT-off",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, to filter the vdW-DF kernel in reciprocal space.
""",
),

Variable(
    abivarname="vdw_df_ndpts",
    varset="vdw",
    vartype="integer",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=20,
    mnemonics="vdW-DF Number of D-mesh PoinTS",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, to build the vdW-DF kernel.
""",
),

Variable(
    abivarname="vdw_df_ngpts",
    varset="vdw",
    vartype="integer",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=-1,
    mnemonics="vdW-DF Number of G-mesh PoinTS",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, to build the vdW-DF kernel.
""",
),

Variable(
    abivarname="vdw_df_nqpts",
    varset="vdw",
    vartype="integer",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=30,
    mnemonics="vdW-DF Number of Q-mesh PoinTS",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, to build the vdW-DF kernel.
""",
),

Variable(
    abivarname="vdw_df_nrpts",
    varset="vdw",
    vartype="integer",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=2048,
    mnemonics="vdW-DF Number of R-PoinTS",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, to define the sampling of the vdW-DF-kernel in real-space.
""",
),

Variable(
    abivarname="vdw_df_nsmooth",
    varset="vdw",
    vartype="integer",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=12,
    mnemonics="vdW-DF Number of SMOOTHening iterations",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, to exponentially smoothen q near q0.
""",
),

Variable(
    abivarname="vdw_df_phisoft",
    varset="vdw",
    vartype="real",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=-1.0,
    mnemonics="vdW-DF PHI value SOFTening.",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, to build the vdW-DF kernel.
""",
),

Variable(
    abivarname="vdw_df_qcut",
    varset="vdw",
    vartype="real",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=5,
    mnemonics="vdW-DF Q-mesh CUT-off",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, to build the vdW-DF kernel.
""",
),

Variable(
    abivarname="vdw_df_qratio",
    varset="vdw",
    vartype="real",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=20,
    mnemonics="vdW-DF, between highest and lowest Q, RATIO.",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0,.
""",
),

Variable(
    abivarname="vdw_df_rcut",
    varset="vdw",
    vartype="real",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=100,
    mnemonics="vdW-DF Real-space CUT-off",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, to define the vdW-DF kernel cut-off radius.
""",
),

Variable(
    abivarname="vdw_df_rsoft",
    varset="vdw",
    vartype="real",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="vdW-DF radius SOFTening.",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, to build the vdW-DF kernel.
""",
),

Variable(
    abivarname="vdw_df_threshold",
    varset="vdw",
    vartype="real",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=0.01,
    mnemonics="vdW-DF energy calculation THRESHOLD",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Sets a threshold for the energy gradient that, when reached, will cause the
vdW-DF interactions to be calculated.
Adjust it to a big value (e.g. 1e12) to enable it all along the SCF
calculation. Too small values, as well as negative values, will result in the
vdW-DF energy contributions never being calculated.
""",
),

Variable(
    abivarname="vdw_df_tolerance",
    varset="vdw",
    vartype="real",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=1e-13,
    mnemonics="vdW-DF global TOLERANCE.",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, to build the vdW-DF kernel.
""",
),

Variable(
    abivarname="vdw_df_tweaks",
    varset="vdw",
    vartype="integer",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="vdW-DF TWEAKS.",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, to build the vdW-DF kernel.

!!! important

    Modifying this variable will likely transform the
    calculated energies and their gradients into garbage.
    You have been warned!
""",
),

Variable(
    abivarname="vdw_df_zab",
    varset="vdw",
    vartype="real",
    topics=['vdw_expert'],
    dimensions="scalar",
    defaultval=-0.8491,
    mnemonics="vdW-DF ZAB parameter",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]]>0",
    added_in_version="before_v9",
    text=r"""
Used when [[vdw_xc]]>0, as introduced in [[cite:Dion2004]].
""",
),

Variable(
    abivarname="vdw_nfrag",
    varset="vdw",
    vartype="integer",
    topics=['vdw_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Van Der Waals Number of interacting FRAGments",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]] in [10,11]",
    added_in_version="before_v9",
    text=r"""
The absolute value of vdw_nfrag is the number of vdW interacting fragments in
the unit cell. As wannierization takes place in reciprocal space, the MLWF
center positions could be translated by some lattice vector from the cell
where atoms are placed. If [[vdw_nfrag]] >= 1 then MLWFs are translated to the
original unit cell, otherwise the program will keep the positions obtained by
Wannier90. The later is usually correct if some atoms are located at the
corners or at limiting faces of the unit cell.
""",
),

Variable(
    abivarname="vdw_supercell",
    varset="vdw",
    vartype="integer",
    topics=['vdw_basic'],
    dimensions=[3],
    defaultval=[0, 0, 0],
    mnemonics="Van Der Waals correction from Wannier functions in SUPERCELL",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]] in [10,11]",
    added_in_version="before_v9",
    text=r"""
Set of dimensionless positive numbers which define the maximum multiples of
the primitive translations ([[rprimd]]) in the supercell construction. Each
component of vdw_supercell indicates the maximum number of cells along both
positive or negative directions of the corresponding primitive vector i.e. the
components of [[rprimd]]. In the case of layered systems for which vdW
interactions occur between layers made of tightly bound atoms, the evaluation
of vdW corrections coming from MLWFs in the same layer (fragment) must be
avoided. Both a negative or null value for one component of [[vdw_supercell]]
will indicate that the corresponding direction is normal to the layers.
""",
),

Variable(
    abivarname="vdw_tol",
    varset="vdw",
    vartype="real",
    topics=['vdw_compulsory'],
    dimensions="scalar",
    defaultval=1e-10,
    mnemonics="Van Der Waals TOLerance",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]] == 5",
    added_in_version="before_v9",
    text=r"""
The DFT-D methods [[cite:Grimme2010]] dispersion potentials, [[vdw_xc]] == 5 or
6 or 7, include a pair potential. The number of pairs of atoms contributing to
the potential is necessarily limited. To be included in the potential a pair
of atom must have contribution to the energy larger than [[vdw_tol]].
""",
),

Variable(
    abivarname="vdw_tol_3bt",
    varset="vdw",
    vartype="real",
    topics=['vdw_basic'],
    dimensions="scalar",
    defaultval=-1,
    mnemonics="Van Der Waals TOLerance for 3-Body Term",
    characteristics=['[[DEVELOP]]'],
    commentdefault="Do include the 3-body term in the correction",
    requires="[[vdw_xc]] == 6",
    added_in_version="before_v9",
    text=r"""
Control the computation of the 3-body correction inside DFT-D3 dispersion
correction (Grimme approach) to the total energy:

  * If **vdw_tol_3bt** <0, no 3-body correction.
  * If **vdw_tol_3bt** >0, the 3-body term is included with a tolerance = **vdw_tol_3bt**.

DFT-D3 as proposed by S. Grimme adds two contributions to the total energy in
order to take into account of the dispersion:

  * A pair-wise potential for which the tolerance is controlled by [[vdw_tol]]

  * A 3-body term which is obtained by summing over all triplets of atoms. Each individual contribution depends of the distances and angles between the three atoms. As it is impossible to sum over all the triplets in a periodic system, one has to define a stopping criterium which is here that an additional contribution to the energy must be higher than **vdw_tol_3bt**

The last term has been predicted to have an important effect for large
molecules [[cite:Grimme2010]]. It is however quite costly in computational
time for periodic systems and seems to lead to an overestimation of lattice
parameters for weakly bound systems [[cite:Grimme2011]]. Still, its
contribution to energy, to forces and to stress is available (not planned for
elastic constants, dynamical matrix and internal strains).
""",
),

Variable(
    abivarname="vdw_typfrag",
    varset="vdw",
    vartype="integer",
    topics=['vdw_basic'],
    dimensions=['[[natom]]'],
    defaultval=MultipleValue(number=1, value='[[natom]]'),
    mnemonics="Van Der Waals TYPe of FRAGment",
    characteristics=['[[DEVELOP]]'],
    requires="[[vdw_xc]] in [10,11]",
    added_in_version="before_v9",
    text=r"""
This array defines the interacting fragments by assigning to each atom an
integer index from 1 to **vdw_nfrag**. The ordering of [[vdw_typfrag]] is the
same as [[typat]] or [[xcart]]. Internally each MLWF is assigned to a given
fragment by computing the distance to the atoms. MLWFs belong to the same
fragment as their nearest atom. The resulting set of MLWFs in each interacting
fragment can be found in the output file in xyz format for easy visualization.
""",
),

Variable(
    abivarname="vdw_xc",
    varset="vdw",
    vartype="integer",
    topics=['vdw_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Van Der Waals eXchange-Correlation functional",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
Selects a van-der-Waals density functional to apply the corresponding
correction to the exchange-correlation energy. If set to zero, no correction
will be applied.
Possible values are:

  * 0: no correction.
  * 1: apply vdW-DF1 (DRSLL) from [[cite:Dion2004]].
  * 2: apply vdw-DF2 (LMKLL) from [[cite:Lee2010]].
  * 5: apply vdw-DFT-D2 as proposed by S. Grimme [[cite:Grimme2006]] (adding a semi-empirical dispersion potential). Available only for ground-state calculations and response functions; see [[vdw_tol]] variable to control convergence.
  * 6: apply vdw-DFT-D3 as proposed by S. Grimme [[cite:Grimme2010]] (refined version of DFT-D2). Available only for ground-state calculations and response functions; see [[vdw_tol]] variable to control convergence and [[vdw_tol_3bt]] variable to include 3-body corrections.
  * 7: apply vdw-DFT-D3(BJ) as proposed by Grimme (based on Becke-Jonhson method from [[cite:Becke2006]]). Available only for ground-state calculations and response functions; see [[vdw_tol]] variable to control convergence.
  * 10: evaluate the vdW correlation energy from maximally localized Wannier functions, as proposed by P. L. Silvestrelli, also known as vdW-WF1 method [[cite:Silvestrelli2008]]. For details on this implementation please check [[cite:Espejo2012]]. The improvements introduced by Andrinopoulos _et al._ [[cite:Andrinopoulos2011]], namely the amalgamation procedure, splitting of p-like MLWFs into two s-like Wannier functions and fractional occupation of MLWFs are performed automatically.
  * 11: evaluate the vdW correlation energy from maximally localized Wannier functions, as proposed by A. Ambrosetti and P. L. Silvestrelli, also known as vdW-WF2 method [[cite:Ambrosetti2012]].
  * 14: apply DFT/vdW-QHO-WF method as proposed by Silvestrelli, which combines the quantum harmonic oscillator-model with localized Wannier functions [[cite:Silvestrelli2013]]. For periodic systems a supercell approach has to be used since **vdw_supercell** is not enabled in this case.

For [[vdw_xc]] = 1 and [[vdw_xc]] = 2, the implementation follows the strategy
devised in the article of Roman-Perez and Soler [[cite:Romanperez2009]].
""",
),

Variable(
    abivarname="vel",
    varset="rlx",
    vartype="real",
    topics=['PIMD_useful', 'MolecularDynamics_basic'],
    dimensions=[3, '[[natom]]'],
    defaultval=MultipleValue(number=None, value=0),
    mnemonics="VELocity",
    characteristics=['[[EVOLVING]]'],
    commentdims="It is represented internally as [[vel]](3,[[natom]],[[nimage]])",
    requires="[[ionmov]] > 0",
    added_in_version="before_v9",
    text=r"""
Gives the starting velocities of atoms, in cartesian coordinates, in
Bohr/atomic time units (atomic time units given where [[dtion]] is described).
For [[ionmov]] = 8 (Nose thermostat), if [[vel]] is not initialized, a random
initial velocity giving the right kinetic energy will be generated.
If the atom manipulator is used, [[vel]] will be related to the preprocessed
set of atoms, generated by the atom manipulator. The user must thus foresee
the effect of this atom manipulator (see [[objarf]]).
Velocities evolve is [[ionmov]] == 1.
""",
),

Variable(
    abivarname="vel_cell",
    varset="rlx",
    vartype="real",
    topics=['PIMD_expert'],
    dimensions=[3, 3],
    defaultval=MultipleValue(number=None, value=3),
    mnemonics="VELocity of the CELL parameters",
    characteristics=['[[EVOLVING]]'],
    commentdims="It is represented internally as [[vel_cell]](3,3,[[nimage]])",
    requires="""[[imgmov]] in [9,13] and [[optcell]] > 0
(Path-Integral Molecular Dynamics with NPT algorithm)""",
    added_in_version="before_v9",
    text=r"""
Irrelevant unless [[imgmov]] = 9 or 13 and [[optcell]]>0 (Path-Integral
Molecular Dynamics with NPT algorithm).
Gives the starting velocities of the dimensional cell parameters in
Bohr/atomic time units (atomic time units given where [[dtion]] is described).
""",
),

Variable(
    abivarname="vis",
    varset="rlx",
    vartype="real",
    topics=['PIMD_basic', 'MolecularDynamics_basic'],
    dimensions="scalar",
    defaultval=100,
    mnemonics="VIScosity",
    added_in_version="before_v9",
    text=r"""
The equation of motion is:

M  I  d  2  R  I  /dt  2  = F  I  - [[vis]] dR  I  /dt

The atomic unit of viscosity is hartree * (atomic time units)/Bohr 2. Units
are not critical as this is a fictitious damping used to relax structures. A
typical value for silicon is 400 with [[dtion]] of 350 and atomic mass 28
[[amu]]. Critical damping is most desirable and is found only by optimizing
[[vis]] for a given situation.

In the case of Path-Integral Molecular Dynamics using the Langevin Thermostat
([[imgmov]] = 9), [[vis]] defines the friction coefficient, in atomic units.
Typical value range is 0.00001-0.001.
""",
),

Variable(
    abivarname="vprtrb",
    varset="ffield",
    vartype="real",
    topics=['Artificial_useful'],
    dimensions=[2],
    defaultval=[0.0, 0.0],
    mnemonics="potential -V- for the PeRTuRBation",
    characteristics=['[[DEVELOP]]', '[[ENERGY]]'],
    requires="[[qprtrb]]",
    added_in_version="before_v9",
    text=r"""
Gives the real and imaginary parts of a scalar potential perturbation. Can be
specified in Ha (the default), Ry, eV or Kelvin, since [[vprtrb]] has the
[[ENERGY]] characteristics.
This is made available for testing responses to such perturbations. The form
of the perturbation, which is added to the local potential, is:

  * ([[vprtrb]](1)$+I$[[vprtrb]](2)$)/2$ at $G=$[[qprtrb]] and
  * ([[vprtrb]](1)$-I$[[vprtrb]](2)$)/2$ at $G=-$[[qprtrb]] (see [[qprtrb]] also).
""",
),

Variable(
    abivarname="w90iniprj",
    varset="w90",
    vartype="integer",
    topics=['Wannier_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Wannier90- INItial PROJections",
    requires="[[prtwant]] == 2 or [[prtwant]] == 3",
    added_in_version="before_v9",
    text=r"""
In order to find the Maximally Localized Wannier Functions, the user has to
provide an initial guess. A set of localized trial orbitals is chosen
corresponding to some rough initial guess at the Wannier Functions, and these
are projected onto the Bloch eigenstates. See [[cite:Souza2002a]].
These initial projections are stored in a file **.amn** and the variable
**w90iniprj** is used to construct them:

  * **w90iniprj** =1: Random projections.

  * **w90iniprj** =2: The initial projections will be a linear combination of hydrogenic atomic orbitals.
The user has to define the projections in the secondary input file
wannier90.win.
Information about how to define them can be found in the manual of Wannier90.
See [www.wannier.org](http://www.wannier.org)
""",
),

Variable(
    abivarname="w90prtunk",
    varset="w90",
    vartype="integer",
    topics=['Wannier_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Wannier90- PRINT UNKp.s file",
    commentdefault="""The default is set to zero because UNKp.s files occupy a lot of memory.""",
    requires="[[prtwant]] == 2 or [[prtwant]] == 3",
    added_in_version="before_v9",
    text=r"""
Defines whether or not the UNKp.s file will be printed.

  * [[w90prtunk]] = 0: Do not print the UNKp.s files

  * [[w90prtunk]] = 1: Print the UNKp.s files on a fine grid

  * [[w90prtunk]]>1: Print the UNKp.s files on a coarse grid

Instead of printing every record we will print every [[w90prtunk]] records. This
is useful to reduce the size of the UNKp.s files, but, the quality is also reduced.

These files contain the periodic part of the bloch states represented on a
regular real space grid. They are indexed by k-point **p** (from 1 to nkpt)
and spin **s** ('1' for 'up','2' for 'down').

The name of the wavefunction file is assumed to have the form:

write(wfnname,200) **p**, **spin**
200 format ('UNK',i5.5,'.',i1)

These file are unformatted. The first line of each file contains 5 integers:
the number of grid points in each direction ( **n1**, **n2** and **n3** ),
the k-point number **ikpt** and the total number of bands mband in the file.
The following rows contain the wavefunctions in real space.

These files are written in the following way for the coarse grid:

```fortran
     write(iun_plot) n1/w90prtunk,n2/w90prtunk,n3/w90prtunk,ikpt,nband
     write(iun_plot) (((fofr(1,jj1,jj2,jj3),fofr(2,jj1,jj2,jj3),&
    &      jj1=1,n1,w90prtunk),jj2=1,n2,w90prtunk),jj3=1,n3,w90prtunk)
```

Where **fofr** is a double precision variable which contains the wavefunctions
in real space. Note that in order to reduce the size of the UNK files we are
just including records in the wavefunctions for 1/(w90prtunk$^3$) of the grid
points. That is why we divide **n1**, **n2** and **n3** by [[w90prtunk]]. The output .xsf files
for plotting with XCrysDen will also be on the coarse grid. When this does not
produce an acceptable plot, [[w90prtunk]] can be set to 1 to output every grid point.
(You should try spline interpolation in XCrysDen first.)
""",
),

Variable(
    abivarname="wfmix",
    varset="gstate",
    vartype="real",
    topics=['Hybrids_useful'],
    dimensions="scalar",
    defaultval=1.0,
    mnemonics="WaveFunctions MIXing factor",
    requires="[[usefock]] > 0 and [[nnsclohf]] >0 and [[fockoptmix]]/100 > 0",
    added_in_version="before_v9",
    text=r"""
When the wavefunctions are determined using a SCF double loop (hybrid
functionals), [[wfmix]] provides the mixing factor to obtain the new input
wavefunctions by the combination of the earlier input wavefunctions and
corresponding (DFT-preconditioned) output wavefunctions at the level of the
outer loop, according to the algorithm specified by [[fockoptmix]]/100. If
[[wfmix]] is 1.0, the output wavefunctions only will determine the new input
wavefunctions. This might possibly lead to instabilities. If [[wfmix]] is
smaller than 1.0, the whole iteration procedure is damped, which might allow
better stability, but might be slower. If it is larger than 1.0, perhaps less
iterations will be needed (if there is no instability).
""",
),

Variable(
    abivarname="wfoptalg",
    varset="dev",
    vartype="integer",
    topics=['SCFAlgorithms_basic'],
    dimensions="scalar",
    defaultval="[[AUTO_FROM_PSP]]",
    mnemonics="WaveFunction OPTimisation ALGorithm",
    characteristics=['[[DEVELOP]]'],
    commentdefault="0 when [[usepaw]] = 0 (norm-conserving pseudopotentials), 10 when [[usepaw]] = 1 (PAW); 114 if [[paral_kgb]] = 1.",
    added_in_version="before_v9",
    text=r"""
Allows one to choose the algorithm for the optimisation of the wavefunctions.
The different possibilities are:

  * [[wfoptalg]] = 0: standard state-by-state conjugate gradient algorithm, with no possibility to parallelize over the states;

  * [[wfoptalg]] = 2: minimisation of the residual with respect to different shifts, in order to cover the whole set of occupied bands,
  with possibility to parallelize over blocks of states (or bands). The number of states in a block is defined in [[nbdblock]]. THIS IS STILL IN DEVELOPMENT.

  * [[wfoptalg]] = 3: minimisation of the residual with respect to a shift. Available only in the non-self-consistent case [[iscf]] = -2,
  in order to find eigenvalues and wavefunctions close to a prescribed value.

  * [[wfoptalg]] = 4: (see also [[wfoptalg]] = 14), a parallel code based on the Locally Optimal Block Preconditioned Conjugate Gradient (LOBPCG)
  method of [[cite:Knyazev2001 ]].
  The implementation rests on the [matlab program by Knyazev](http://www.mathworks.com/matlabcentral/fileexchange/48-lobpcg-m) [[cite:Knyazev2007]].
  For more information see [[cite:Bottin2008]]

  * [[wfoptalg]] = 10: (for PAW) standard state-by-state conjugate gradient algorithm, with no possibility to parallelize over the states,
  but modified scheme described in [[cite:Kresse1996]] (modified kinetic energy, modified preconditionning, minimal orthogonalization, ...);

  * [[wfoptalg]] = 14: the recommended for parallel code, the same as [[wfoptalg]] = 4 except that the preconditioning of the block vectors does not
  depend on the kinetic energy of each band, and the orthogonalization after the LOBPCG algorithm is no longer performed. The first modification increases the convergence and the second one the efficiency.

  * [[wfoptalg]] = 114: A new version of [[wfoptalg]] = 14 which is more efficient for few blocks and can take advantage of OpenMP if abinit is compiled with a multithreaded linear algebra library.
  With more than 1 thread [[npfft]] shoud NOT be used for the time being.

  * [[wfoptalg]] = 1: new algorithm based on Chebyshev filtering, designed for very large number of processors, in the regime
  where LOBPCG does not scale anymore. It is not able to use preconditionning and therefore might converge slower than other algorithms.
  By design, it will **not** converge the last bands: it is recommended to use slightly more bands than necessary.
  For usage with [[tolwfr]], it is imperative to use [[nbdbuf]]. For more performance, try [[use_gemm_nonlop]].
  For more information, see the [performance guide](../../theory/howto_chebfi.pdf) and the [[cite:Levitt2015]]. Status: experimental but usable.
  Questions and bug reports should be sent to antoine (dot) levitt (at) gmail.com.
""",
),

Variable(
    abivarname="wtatcon",
    varset="rlx",
    vartype="real",
    topics=['GeoConstraints_useful'],
    dimensions=[3, '[[natcon]]', '[[nconeq]]'],
    defaultval=0,
    mnemonics="WeighTs for AToms in CONstraint equations",
    characteristics=['[[NO_MULTI]]'],
    added_in_version="before_v9",
    text=r"""
Gives the weights determining how the motion of atoms is constrained during
structural optimization or molecular dynamics (see [[nconeq]], [[natcon]],
and [[iatcon]]). For each of the [[nconeq]] independent constraint equations,
wtatcon is a 3*[[natcon]] array giving weights, W  I, for the x, y, and z
components of each of the atoms (labeled by I) in the list of indices
[[iatcon]]. Prior to taking an atomic step, the calculated forces, F  I, are
replaced by projected forces, F'  I, which satisfy the set of constraint equations

Sum  mu=x,y,z; I=1,natcon: W  mu,I  * F'  mu,I  = 0 for each of the [[nconeq]] arrays W  I.

Different types of motion constraints can be implemented this way. For example,

    nconeq 1 natcon 2 iatcon 1 2 wtatcon 0 0 +1 0 0 -1

could be used to constrain the relative height difference of two adsorbate
atoms on a surface (assuming their masses are equal), since F'  z,1  - F'
z,2  = 0 implies z  1  - z  2  = constant.
""",
),

Variable(
    abivarname="wtk",
    varset="basic",
    vartype="real",
    topics=['k-points_useful'],
    dimensions=['[[nkpt]]'],
    defaultval=MultipleValue(number='[[nkpt]]', value=1.0),
    mnemonics="WeighTs for K points",
    commentdefault="Except when [[kptopt]]/=0",
    added_in_version="before_v9",
    text=r"""
Gives the k point weights.
The k point weights will have their sum (re)normalized to 1 (unless
[[occopt]] = 2 and [[kptopt]] = 0; see description of [[occopt]]) within the
program and therefore may be input with any arbitrary normalization. This
feature helps avoid the need for many digits in representing fractional weights such as 1/3.
[[wtk]] is ignored if [[iscf]] is not positive, except if [[iscf]] = -3.
""",
),

Variable(
    abivarname="wtq",
    varset="gstate",
    vartype="real",
    topics=['q-points_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="WeighTs for the current Q-points",
    commentdefault="Except when [[qptopt]]/=0",
    added_in_version="before_v9",
    text=r"""
Gives the current q-point weight.
""",
),

Variable(
    abivarname="wvl_bigdft_comp",
    varset="gstate",
    vartype="integer",
    topics=['Wavelets_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="WaVeLet BIGDFT Comparison",
    added_in_version="before_v9",
    text=r"""
This variable is used for the wavelets capabilities of ABINIT (see [[usewvl]]).
It is used to compare the results obtained with ABINIT with those obtained
with BigDFT stand-alone. When it is set to 1, ABINIT will follow the workflow
as in BigDFT stand-alone. Therefore, the results must be exactly the same with the two codes.
""",
),

Variable(
    abivarname="wvl_crmult",
    varset="gstate",
    vartype="real",
    topics=['Wavelets_basic'],
    dimensions="scalar",
    defaultval=6.0,
    mnemonics="WaVeLet Coarse grid Radius MULTiplier",
    added_in_version="before_v9",
    text=r"""
This factor is used to define the expansion of the coarse resolution grid in
the case of wavelets (see [[usewvl]]). The grid is made of points inside
spheres centered on atoms. The radius of these spheres are the product between
this factor and the covalent radius of element (read from the pseudo-potential file).
This factor is responsible for the amount of used memory (see also [[wvl_hgrid]]).
""",
),

Variable(
    abivarname="wvl_frmult",
    varset="gstate",
    vartype="real",
    topics=['Wavelets_basic'],
    dimensions="scalar",
    defaultval=10.0,
    mnemonics="WaVeLet Fine grid Radius MULTiplier",
    added_in_version="before_v9",
    text=r"""
This factor is used to define the expansion of the fine resolution grid in
the case of wavelets (see [[usewvl]]). This fine resolution grid has the same
grid step than the coarse one (see [[wvl_crmult]] ), but on each point, 8
coefficients are stored instead of one, increasing the precision of the
calculation in this area. The grid is made of points inside spheres centered
on atoms. The radius of these spheres are the product between this factor and
a value read from the pseudo-potential file.
This factor is responsible for the amount of used memory (see also [[wvl_hgrid]]).
""",
),

Variable(
    abivarname="wvl_hgrid",
    varset="basic",
    vartype="real",
    topics=['Wavelets_basic'],
    dimensions="scalar",
    defaultval=0.5,
    mnemonics="WaVeLet H step GRID",
    characteristics=['[[LENGTH]]'],
    added_in_version="before_v9",
    text=r"""
It gives the step size in real space for the grid resolution in the wavelet
basis set. This value is highly responsible for the memory occupation in the
wavelet computation. The value is a length in atomic units.
""",
),

Variable(
    abivarname="wvl_ngauss",
    varset="gstate",
    vartype="integer",
    topics=['Wavelets_expert'],
    dimensions=[2],
    defaultval=[1, 100],
    mnemonics="WaVeLet Number of GAUSSians",
    added_in_version="before_v9",
    text=r"""
In the wavelet-PAW computation case, projectors may be fitted to a sum of
complex Gaussians. The fit is done for [[wvl_ngauss]](1), [[wvl_ngauss]](1)+1... up
to [[wvl_ngauss]](2) Gaussians.
""",
),

Variable(
    abivarname="wvl_nprccg",
    varset="gstate",
    vartype="integer",
    topics=['Wavelets_expert'],
    dimensions="scalar",
    defaultval=5,
    mnemonics="WaVeLet maximum Number of PReConditioner Conjugate Gradient iterations",
    added_in_version="before_v9",
    text=r"""
In the wavelet computation case, the wavefunctions are directly minimised
using a real-space preconditioner. This preconditioner has internally some
conjugate gradient iterations. This value defines a boundary for the number of
conjugate gradient iterations on each wavefunction convergence step.
""",
),

Variable(
    abivarname="xangst",
    varset="basic",
    vartype="real",
    topics=['crystal_compulsory'],
    dimensions=[3, 'min([[natom]],[[natrd]])'],
    mnemonics="vectors (X) of atom positions in cartesian coordinates -length in ANGSTrom-",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Gives the cartesian coordinates of atoms within unit cell, in angstrom. This
information is redundant with that supplied by array [[xred]] or [[xcart]].
If [[xred]] and [[xangst]] are ABSENT from the input file and [[xcart]] is
provided, then the values of [[xred]] will be computed from the provided
[[xcart]] (i.e. the user may use xangst instead of [[xred]] or [[xcart]] to
provide starting coordinates).
One and only one of [[xred]], [[xcart]] and [[xangst]] must be provided.
The conversion factor between Bohr and Angstrom is 1 Bohr=0.5291772108
Angstrom, see the [NIST site](http://physics.nist.gov/cuu/Constants/index.html).
Atomic positions evolve if [[ionmov]]/=0. In constrast with [[xred]] and
[[xcart]], [[xangst]] is not internal.
""",
),

Variable(
    abivarname="xc_denpos",
    varset="dev",
    vartype="real",
    topics=['xc_expert'],
    dimensions="scalar",
    defaultval=1e-14,
    mnemonics="eXchange-Correlation - DENsity POSitivity value",
    characteristics=['[[DEVELOP]]'],
    added_in_version="before_v9",
    text=r"""
For the evaluation of the exchange-correlation functionals, the density cannot
be negative, or even too small (e.g. the LDA exchange kernel behaves like the
density at power -(2/3), and the density is used at the denominator of
different factors in GGAs and metaGGAs. [[xc_denpos]] is the smallest value
that the density can assume at the time of the evaluation of a XC functional,
in ABINIT. When then computed density drops below [[xc_denpos]] before
attacking the evaluation of the XC functional, then it will be (only for that
purpose) replaced by [[xc_denpos]]. Note that the evaluation of the gradients
or other quantities that are density-dependent is performed before this replacement.

It has been observed that the SCF cycle of the Tran-Blaha mGGA can be quite
hard to make converge, for systems for which there is some vacuum. In this
case, setting [[xc_denpos]] to 1.0e-7 ... 1.0e-6 has been seen to allow good
convergence. Of course, this will affect the numerical results somehow, and
one should play a bit with this value to avoid incorrect calculations.
""",
),

Variable(
    abivarname="xc_tb09_c",
    varset="dev",
    vartype="real",
    topics=['xc_expert'],
    dimensions="scalar",
    defaultval=99.99,
    mnemonics="Value of the c parameter in the eXchange-Correlation TB09 functional",
    added_in_version="before_v9",
    text=r"""
The modified Becke-Johnson exchange-correlation functional by
[[cite:Tran2009 | Tran and Blaha]] reads:

$$ V_x(r) =
c V_x^{BR}(r) +
(3c - 2) \frac{1}{\pi} \sqrt{\frac{5}{12}}
\sqrt{2 \frac{t(r)}{\rho(r)}} $$

where $\rho(r)$ is the electron density,
$t(r)$ is the kinetic-energy density, and
$ V_x^{BR}(r)$ is the Becke-Roussel potential.

In this equation the parameter $c$ can be evaluated at each SCF step according
to the following equation:

$$ c = \alpha + \beta
\left( \frac{1}{V_{cell}}
\int_{V_{cell}} \frac{|\nabla \rho(r)|}{\rho(r)}
d^3r \right)^{1/2} $$

The $c$ parameter is evaluated thanks to the previous equation when xc_tb09_c is
equal to the "magic" default value 99.99. The $c$ parameter can also be fixed to
some (property-optimized or material-optimized) value by using this variable.
""",
),

Variable(
    abivarname="xcart",
    varset="basic",
    vartype="real",
    topics=['crystal_compulsory'],
    dimensions=[3, 'min([[natom]],[[natrd]])'],
    mnemonics="vectors (X) of atom positions in CARTesian coordinates",
    characteristics=['[[EVOLVING]]', '[[LENGTH]]'],
    added_in_version="before_v9",
    text=r"""
Gives the cartesian coordinates of atoms within unit cell. This information is
redundant with that supplied by array [[xred]] or [[xangst]]. By default,
[[xcart]] is given in Bohr atomic units (1 Bohr=0.5291772108 Angstroms),
although Angstrom can be specified, if preferred, since [[xcart]] has the [[LENGTH]] characteristics.
If [[xred]] and [[xangst]] are ABSENT from the input file and [[xcart]] is
provided, then the values of [[xred]] will be computed from the provided
[[xcart]] (i.e. the user may use [[xcart]] instead of [[xred]] or [[xangst]]
to provide starting coordinates).
One and only one of [[xred]], [[xcart]] and **xangst** must be provided.
Atomic positions evolve if [[ionmov]]/=0.
""",
),

Variable(
    abivarname="xclevel",
    varset="internal",
    vartype="integer",
    topics=['xc_internal', 'TDDFT_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="eXchange Correlation functional LEVEL",
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Automatically determined from the value of [[ixc]].

  * 0  -->  No XC contribution.
  * 1  -->  LDA functional.
  * 2  -->  GGA functional or hybrid functional based on GGA.
  * 3  -->  Functional for TDDFT.
""",
),

Variable(
    abivarname="xred",
    varset="basic",
    vartype="real",
    topics=['crystal_compulsory'],
    dimensions=[3, 'min([[natom]],[[natrd]])'],
    defaultval=MultipleValue(number=None, value=0.0),
    mnemonics="vectors (X) of atom positions in REDuced coordinates",
    characteristics=['[[EVOLVING]]'],
    commentdims="represented internally as xred(3,[[natom]],[[nimage]])",
    added_in_version="before_v9",
    text=r"""
Gives the atomic locations within unit cell in coordinates relative to real
space primitive translations (**NOT in cartesian coordinates**). Thus these are
fractional numbers typically between 0 and 1 and are dimensionless. The
cartesian coordinates of atoms (in Bohr) are given by:
R_cartesian = xred1 * rprimd1 + xred2 * rprimd2 + xred3 * rprimd3
where (xred1,xred2,xred3) are the "reduced coordinates" given in columns of
"[[xred]]", (rprimd1,rprimd2,rprimd3) are the columns of primitive vectors
array "[[rprimd]]" in Bohr.

If you prefer to work only with cartesian coordinates, you may work entirely
with "[[xcart]]" or "[[xangst]]" and ignore [[xred]], in which case [[xred]]
must be absent from the input file.
One and only one of [[xred]], [[xcart]] and [[xangst]] must be provided.
Atomic positions evolve if [[ionmov]]/=0.
""",
),

Variable(
    abivarname="xredsph_extra",
    varset="gstate",
    vartype="real",
    topics=['printing_prdos'],
    dimensions=[3, '[[natsph_extra]]'],
    defaultval=MultipleValue(number=None, value=0.0),
    mnemonics="X(position) in REDuced coordinates of the SPHeres for dos projection in the EXTRA set",
    requires="[[natsph_extra]] > 0",
    added_in_version="before_v9",
    text=r"""
The positions in reduced coordinates of extra spheres used in the DOS
projection, simulating an STS signal. See [[natsph_extra]] for a more complete description.
""",
),

Variable(
    abivarname="xyzfile",
    varset="geo",
    vartype="string",
    topics=['crystal_useful'],
    dimensions="scalar",
    mnemonics="XYZ FILE input for geometry",
    characteristics=['[[INPUT_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Gives the name of a xyz format file, to take [[natom]], [[ntypat]], [[typat]],
[[znucl]], and [[xangst]] from. This input can not be mixed with normal atom
specifications for other datasets.

Notes: do not quote the file name in the abinit input file, simply leave a
space after xyzfile. The xyz format is the number of atoms on the first line,
a comment line, then one line per atom, with the element as a 2 letter symbol
("As" "O" or "Pu") and the three cartesian coordinates in Angstrom.
""",
),

Variable(
    abivarname="zcut",
    varset="gw",
    vartype="real",
    topics=['Susceptibility_expert', 'BSE_expert', 'SelfEnergy_expert'],
    dimensions="scalar",
    defaultval=0.0036749326,
    mnemonics="Z-CUT",
    characteristics=['[[ENERGY]]'],
    commentdefault="0.1 eV (0.0036749326 Ha)",
    requires="[[optdriver]] in [3, 4, 7, 99]",
    added_in_version="before_v9",
    text=r"""
It is meant to avoid some divergences that might occur during the evaluation
of the Adler-Wiser expression of the irreducible polarizability
([[optdriver]] = 3) or during the numerical treatment of the integrals defining
the contribution to the self-energy matrix elements ([[optdriver]] = 4). If the
denominator becomes smaller than **zcut**, a small imaginary part (depending
on **zcut**) is added, in order to avoid the divergence.

When [[optdriver]] = 99, **zcut** defines the small complex shift used to avoid
divergences in the expression for the macroscopic dielectric function. It
simulates the experimental uncertainty and the finite lifetime of the
quasi-particles (although the true lifetime should be k- and band-dependent).
The value of **zcut** affects the number of iteration needed to achieve
convergence in the Haydock iterative method. In this case, **zcut** should be
larger than the typical distance between the eigenvalues of the exciton Hamiltonian.
Ideally, one should make a convergence study decreasing the value of **zcut**
for increasing number of k points.

When [[optdriver]] = 7, **zcut** defines the small complex shift used to avoid
divergences in the expression for the Fan-Migdal e-ph self-energy.
Note that the default value is to large for e-ph calculations, smaller values of the order
of 0.001 or 0.001 eV should be used (and carefully tested).
""",
),

Variable(
    abivarname="zeemanfield",
    varset="ffield",
    vartype="real",
    topics=['MagField_basic'],
    dimensions=[3],
    defaultval=0,
    mnemonics="ZEEMAN FIELD",
    characteristics=['[[MAGNETIC_FIELD]]'],
    added_in_version="before_v9",
    text=r"""
Give the value of the Zeeman field, $H$, acting on the spinorial wavefunctions.
Note that Tesla are admitted. This sets the magnitude of $\mu_0H$, in Tesla,
with H in Amperes/metre.
""",
),

Variable(
    abivarname="ziontypat",
    varset="internal",
    vartype="real",
    topics=['AtomTypes_internal', 'PseudosPAW_internal'],
    dimensions=['[[ntypat]]'],
    defaultval="[[AUTO_FROM_PSP]]",
    mnemonics="Z (charge) of the IONs for the different TYPes of AToms",
    characteristics=['[[INTERNAL_ONLY]]'],
    added_in_version="before_v9",
    text=r"""
Charge of the pseudo-ion (defined as the number of valence electrons that are needed to
screen exactly the pseudopotential).
""",
),

Variable(
    abivarname="znucl",
    varset="basic",
    vartype="real",
    topics=['AtomTypes_compulsory', 'PseudosPAW_compulsory'],
    dimensions=['[[npsp]]'],
    mnemonics="charge -Z- of the NUCLeus",
    characteristics=['[[NO_MULTI]]'],
    added_in_version="before_v9",
    text=r"""
Gives nuclear charge for each type of pseudopotential, in order.
If [[znucl]] does not agree with nuclear charge, as given in pseudopotential
files, the program writes an error message and stops.

!!! note

    In the pseudopotential files, [[znucl]] is called "zatom".

For a "dummy" atom, with [[znucl]] = 0, as used in the case of calculations
with only a jellium surface, ABINIT sets arbitrarily the covalent radius to one.
""",
),
#{"abinit_version": "8.7.3"}
Variable(
    abivarname="tmesh",
    varset="eph",
    topics=['ElPhonInt_basic'],
    vartype="real",
    defaultval=[5.0, 59.0, 6.0],
    dimensions=[3],
    mnemonics="Temperature MESH",
    added_in_version="8.7.3",
    text=r"""
This variable defines the linear mesh of temperatures used in the EPH code ([[optdriver]] = 7).
The first entry gives the initial temperature in Kelvin, the second entry the linear step in Kelvin,
the third entry is the number of points in the mesh. The default value corresponds to 6 points between 5 K and 300 K.
""",
),

Variable(
    abivarname="prtkbff",
    varset="files",
    topics=['printing_prden'],
    vartype="integer",
    defaultval=0,
    dimensions="scalar",
    requires="[[iomode]] == 3",
    mnemonics="PRinT Kleynman-Bylander Form Factors",
    added_in_version="8.7.3",
    text=r"""
This input variable activates the output of the Kleynman-Bylander form factors in the **netcdf** WFK file
produced at the end of the ground-state calculation. Remember to set [[iomode]] to 3.

The form factors are needed to compute the matrix elements of the commutator [Vnl, r]
of the non-local part of the (NC) pseudopotentials.
This WFK file can therefore be used to perform optical and/or many-body calculations with external codes such as DP/EXC and Yambo.
The option is ignored if PAW.

!!! important

    At the time of writing (|today|, [[istwfk]] must be set to 1 for all k-points in the IBZ
    since external codes do not support wavefunctions given on the reduced G-sphere.
    Moreover [[useylm]] must be 0 (default if NC pseudos).
""",
),

#{"abinit_version": "9.0.0"}
Variable(
    abivarname="sigma_ngkpt",
    varset="gw",
    topics=['SelfEnergy_useful'],
    vartype="integer",
    defaultval=0,
    dimensions=[3],
    requires="[[optdriver]] in [4, 7]",
    mnemonics="SIGMA: Number of Grid points for K PoinTs generation",
    added_in_version="9.0.0",
    text=r"""
This variable allows the user to specify the list of k-points in the self-energy $\Sigma_{n\kk}$
in terms of a homogeneous mesh in the IBZ instead of the traditional approach based
on [[nkptgw]], [[kptgw]], [[bdgw]].

The specification in terms of sigma_ngkpt is easier to use in particular when
the self-energy is needed on a sub-mesh.
The use of this variables requires a range of bands specified via [[gw_qprange]].

!!! important

    sigma_ngkpt and [[nkptgw]] and [[sigma_erange]] are mutually exclusive.
""",
),

Variable(
    abivarname="sigma_nshiftk",
    varset="gw",
    topics=['SelfEnergy_useful'],
    vartype="integer",
    defaultval=0,
    dimensions="scalar",
    requires="[[optdriver]] in [4, 7] and [[sigma_shiftk]]",
    mnemonics="SIGMA: Number of SHIFTs for K point grids",
    added_in_version="9.0.0",
    text=r"""
The number of shifts in [[sigma_shiftk]].
""",
),

Variable(
    abivarname="sigma_shiftk",
    varset="gw",
    topics=['SelfEnergy_useful'],
    vartype="integer",
    defaultval=[0, 0, 0],
    dimensions=[3, '[[sigma_nshiftk]]'],
    requires="[[optdriver]] in [4, 7] and [[sigma_nshiftk]]",
    mnemonics="SHIFT for K points",
    excludes="[[sigma_erange]] or [[nkptgw]]",
    added_in_version="9.0.0",
    text=r"""
The shifts of the k-mesh used to define the list of k-points for the computation of the
electron self-energy $\Sigma_{n\kk}$.
See also [[sigma_nshiftk]].

!!! important

   This variable is not compatible with [[nkptgw]] and [[sigma_erange]].
""",
),

Variable(
    abivarname="wfk_task",
    varset="gstate",
    topics=['ElecBandStructure_useful'],
    vartype="string",
    defaultval=0,
    dimensions="scalar",
    requires="[[optdriver]] == 8",
    mnemonics="WFK TASK",
    added_in_version="9.0.0",
    text=r"""

This variable defines the quantity to compute starting from a previously generated WFK file.
Possible values are:

  * "wfk_full" --> Read WFK file and produce new WFK file with k-points in the full BZ.
        Wavefunctions with [[istwfk]] > 2 are automatically converted into the full G-sphere representation.
        This option can be used to interface Abinit with external tools requiring k-points in the full BZ.

  * "wfk_einterp" --> Read energies from WFK file and interpolate band structure using the parameters specified by [[einterp]].

  * "wfk_ddk" --> Compute DDK matrix elements for all bands and k-points in the WFK file.
     The contribution due to the non-local part of the pseudopotential can be ignored
     by setting [[inclvkb]] = 0 (not recommended unless you know what you are doing).

  * "wfk_kpts_erange" --> Read WFK file,
        use star-function and [[einterp]] parameters to interpolate electron energies onto fine k-mesh
        defined by [[sigma_ngkpt]] and [[sigma_shiftk]].
        Find k-points inside (electron/hole) pockets according to the values specified in [[sigma_erange]].
        Write KERANGE.nc file with the tables required by the code to automate NSCF band structure calculations
        inside the pocket(s) and electron lifetime computation in the EPH code when [[eph_task]] = -4.
""",
),

Variable(
    abivarname="sigma_bsum_range",
    varset="gw",
    topics=['SelfEnergy_expert'],
    vartype="integer",
    defaultval=[0, 0],
    dimensions=[2],
    requires="[[optdriver]] in [7]",
    mnemonics="SIGMA: Band SUM RANGE",
    added_in_version="9.0.0",
    text=r"""
This variable allows the user to specify the range of bands in the sum over states for the e-ph self-energy $\Sigma_{n\kk}$.
If not specified, the code includes all the states from 1 up to [[nband]].
Note that this option can be used only when computing both the real and imaginary part of the self-energy.
In the calculation of electron linewidths, indeed, the states are automatically selected using an energy window
that takes into account the maximum phonon frequency.
""",
),

# TODO: Remove
#Variable(
#    abivarname="frohl_params",
#    varset="gw",
#    topics=['SelfEnergy_expert'],
#    vartype="real",
#    defaultval=[0, 0, 0, 0],
#    dimensions=[4],
#    requires="[[optdriver]] in [7]",
#    mnemonics="FROHLich PARAMeterS",
#    text=r"""
#This variable is still under development. User interface will change
#""",
#),

Variable(
    abivarname="prteliash",
    varset="eph",
    topics=['SelfEnergy_expert'],
    vartype="integer",
    defaultval=0,
    dimensions="scalar",
    requires="[[optdriver]] in [7]",
    mnemonics="PRINT ELIASHberg function.",
    added_in_version="9.0.0",
    text=r"""
This variable controls the output of the generalized Eliashberg function when [[eph_task]] is +4 or -4.
If set 1, the EPH code will compute the generalized Eliashberg function and will save the results in the SIGEPH.nc file.
""",
),

Variable(
    abivarname="sigma_erange",
    varset="eph",
    topics=['SelfEnergy_expert'],
    vartype="real",
    defaultval=[-1.0, -1.0],
    dimensions=[2],
    mnemonics="SIGMA Energy-range.",
    characteristics=['[[ENERGY]]'],
    added_in_version="9.0.0",
    text=r"""
This variable selects the k-points and the bands in the self-energy matrix elements on the basis
of their position with respect to the band edges (energy differences are **always positive**, even for holes).

Only the k-points and the bands whose energy difference if less than this value will be included in the calculation.
The first entry refers to holes, the second one to electrons.
A negative entry can be used to exclude either holes or electrons from the calculation.
This variable is not compatible with [[nkptgw]] and [[sigma_ngkpt]].

!!! important

    By default, this variable is given in Hartree. Use

        sigma_erange 1 1 eV

    to specify the energy intervals in eV units.
""",
),

Variable(
    abivarname="eph_tols_idelta",
    varset="eph",
    topics=['SelfEnergy_expert'],
    vartype="real",
    defaultval=[1e-12, 1e-12],
    dimensions=[2],
    mnemonics="EPH TOLeranceS on Integral of DELTA.",
    added_in_version="9.0.0",
    text=r"""
This variable can be used to introduce a cutoff on the q-points when computing the imaginary
part of the electron-phonon self-energy ([[eph_task]] = -4) with the tetrahedron method ([[eph_intmeth]] = 2).
The first entry refers to phonon absorption while the second one is associated to phonon emission.
A q-point is included in the sum of the tetrahedron weights for phonon absorption/emission are larger that these values.
""",
),

Variable(
    abivarname="eph_phrange",
    varset="eph",
    topics=['SelfEnergy_expert'],
    vartype="real",
    defaultval=[0, 0],
    dimensions=[2],
    mnemonics="EPH PHonon mode RANGE.",
    added_in_version="9.0.0",
    text=r"""
This variable is used to select the range of phonon modes included in the computation of the electron-phonon self-energy.
By default all phonon modes are included ([0, 0]), otherwise only the phonon modes with index between the first
and second entry are included.
""",
),

Variable(
    abivarname="transport_ngkpt",
    varset="eph",
    topics=['SelfEnergy_expert'],
    vartype="integer",
    defaultval=[0, 0, 0],
    dimensions=[3],
    requires="[[optdriver]] == 7 and [[eph_task]] in [-4,7]",
    mnemonics="TRANSPORT: Number of Grid points for K PoinTs integration in transport computations",
    added_in_version="9.0.0",
    text=r"""
This in an advanced option that is mainly used to downsample the k-points used
to compute the carrier mobility obtained with ([[eph_task]] = -4 or [[eph_task]] = 7).

If this variable is not specified, the code uses the k-mesh specified by [[ngkpt]]
(i.e. the k-mesh corresponding to the WFK file) for computing the mobility integral in the BZ.
In some cases, however, one may want to employ a submesh of k-points to analyze the convergence behaviour.
For instance one may have performed a calculation with a 100x100x100 k-mesh and may be interested in the values
obtained with a 50x50x50 without having to perform a full lifetime calculation on a 50x50x50 from scratch.
""",
),

Variable(
    abivarname="eph_restart",
    varset="eph",
    topics=['ElPhonInt_basic'],
    vartype="integer",
    defaultval=0,
    dimensions="scalar",
    mnemonics="EPH RESTART.",
    added_in_version="9.0.0",
    text=r"""
This variable can be used to restart an EPH calculation.
At present, this feature is supported only when computing the electron-phonon self-energy ([[eph_task]] = 4, -4).
In this case, the code will look for a **pre-existing** SIGEPH.nc file and will compute the remaining k-points
provided that the metadata found in the netcdf file is compatible with the input variables specified in the input file.
The code aborts if the metadata reported in the SIGEPH.nc file is not compatible with the input file.
Note that the restart in done **in-place** that is the output SIGEPH.nc is used as input of the calculation so there is no
need to specify getsigeph or irdsigeph input variables.
""",
),

Variable(
    abivarname="eph_stern",
    varset="eph",
    vartype="integer",
    topics=['ElPhonInt_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Electron-PHonon: use STERNheimer approach to replace sum over empty states.",
    requires="[[tolwfr]] > 0",
    added_in_version="9.0.0",
    text=r"""
This variable activates the Sternheimer method in the calculation of the e-ph self-energy ([[eph_task]] == 4)
This technique replaces the explicit sum over empty states **above** [[nband]]
with the NSCF computation of the first order derivative of the KS wavefunctions (actually
the projection in the subspace orthogonal to the nband states).

The Sternheimer approach requires an external file with the KS potential produced by setting [[prtpot]] = 1
during the GS run and the specification of [[tolwfr]] in the EPH input file.
The path to the POT file used in the EPH calculation is specified via [[getpot_path]].
The number of line minimisations for the Sternheimer solver is defined by [[nline]].

!!! important

    The Sternheimer approach approximates the e-ph self-energy with the adiabatic expression
    in which phonon frequencies are neglected and the frequency dependence of $\Sigma_{n\kk}(\omega)$ is
    replaced by $\Sigma_{n\kk}(\ee_{n\kk})$.
    This approximation is valid provided that **enough** bands above the states of interest are explicitly included.
    The calculation should therefore be converged with respect to the value of [[nband]].
    Note however that the memory requirements and the computational cost of the Sternheimer solver increases with **nband**
    as this part is not yet parallelized.
""",
),

Variable(
    abivarname="getkerange_path",
    varset="eph",
    vartype="string",
    topics=['ElPhonInt_expert'],
    dimensions="scalar",
    defaultval=None,
    mnemonics="KERANGE PATH",
    added_in_version="9.0.0",
    text=r"""
This variable defines the path of the external KERANGE.nc file with the list of k-points in the electron/hole pockets.
The tables stored in the file are used for the calculation of the imaginary part of the e-ph self-energy ([[eph_task]] == -4)
This file is generated by running a preliminary step with [[wfk_task]] = "wfk_einterp".
""",
),

Variable(
    abivarname="symv1scf",
    varset="eph",
    vartype="integer",
    topics=['ElPhonInt_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="SYMmetrize V1 DFPT SCF potentials",
    added_in_version="9.0.0",
    text=r"""
If *symv1scf* is equal to 1, the spatial-symmetry on the first-order DFPT potentials
is enforced every time a set of potentials in the BZ is reconstructed by symmetry
starting from the initial values in the IBZ.
This option is similar to [[symdynmat]] but it acts on the DFPT potentials instead of
the dynamical matrix.
""",
),

Variable(
    abivarname="dvdb_add_lr",
    varset="eph",
    vartype="integer",
    topics=['ElPhonInt_expert'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="DVDB ADD Long-Range part when interpolating DFPT potentials.",
    added_in_version="9.0.0",
    text=r"""
This flag is used in the Fourier interpolation in q-space of the DFPT potentials.
In polar materials there is a long range (LR) component in the first-order variation
of the KS potential that can be modeled in terms of the Born effective charges and
the macroscopic dielectric tensor [[cite:Verdi2015]], [[cite:Giustino2017]].
Possible values are [0, -1, 1].

Setting this flag to 0 deactivates the treatment of the LR contribution (not recommended in polar materials).

If *dvdb_add_lr* is set to 1, the LR part is removed when computing the real-space representation
of the DFPT potentials so that the potential in real space is short-ranged and ameneable to Fourier interpolation.
The long-range contribution is then added back when interpolating the DFPT potentials at arbitrary q-points

If *dvdb_add_lr* is set to -1, the LR part is removed before computing the real-space representation
but the LR term is **not** reintroduced during the interpolation in $\qq$-space.
This option is mainly used for debugging purposes.

By default, the code will always treat the LR term if the DDB file contains the Born effective charges
and the macroscopic dielectric tensor.
This option is similar to [[dipdip]] but it acts on the DFPT potentials instead of the dynamical matrix.
""",
),

Variable(
    abivarname="eph_np_pqbks",
    varset="eph",
    vartype="integer",
    topics=['ElPhonInt_expert'],
    dimensions=[5],
    defaultval=0,
    mnemonics="EPH Number of Processors for Perturbations, Q-points, Bands, K-points, Spin.",
    added_in_version="9.0.0",
    text=r"""
This variable defines the Cartesian grid of MPI processors used for EPH calculations.
If not specified in the input, the code will generate this grid automatically using the total number of processors
and the basic dimensions of the job computed at runtime.
At present (|today|), this variable is supported only in the calculation of the phonon linewidths ([[eph_task]] 1)
and in the computation of the e-ph self-energy ([[eph_task]] 4 or -4).
In all the other tasks, this variable is ignored.

Preliminary considerations:

EPH calculations require very dense samplings of the BZ to converge and the memory requirements
increase quickly with the number of k-points, q-points and [[natom]].
The EPH code can MPI-distribute the most important datastructures but non all the MPI-levels
present the same scalability and the same parallel efficiency.
Besides the maximum number of MPI processes that can be used for the different MPI-levels is related
to the basic dimensions of the calculation.

In what follows, we briefly explain the pros and cons of the different MPI-levels, then we specialize
the discussion to the different calculations activated by [[eph_task]].

The parallelization over perturbations (**np**) is network intensive but it allows one to decrease the memory
needed for the DFPT potentials especially when computing the e-ph self-energy.
The maximum value for **np** is 3 * [[natom]] and the workload is equally distributed provided **np**
divides 3 * [[natom]] equally.
Using **np** == [[natom]] usually gives good parallel efficiency.

The parallelization over bands (**nb**) has limited scalability that depends on the number of bands included
in the self-energy but it allows one to reduce the memory
allocated for the wavefunctions, especially when we have to sum over empty states in the e-ph self-energy.

[[eph_task]] = +1
    By default, the code uses all the processes for the (k-point, spin) parallelism.
    Since the number of k-points around the FS is usually large, this parallelization scheme is OK in most of the cases.
    When the number of processes becomes comparable to the number of k-points around the FS,
    it makes sense to activate the q-point parallelism.
    The parallelism over perturbations should be used to reduce the memory allocated for the interpolation of the DFPT potentials.
    The band parallelism is not supported in this part.

[[eph_task]] = +4
    Parallelization over bands allows one to reduce the memory needed for the wavefunctions but
    this level is less efficient than the parallelization over q-points and perturbations.
    To avoid load and memory imbalance, **nb** should divide [[nband]].
    We suggest to increase the number of procs for bands until the memory allocated for the wavefunctions
    decreases to a reasonable level and then use the remaining procs for **nq** and **np** in this order
    until these levels start to saturate.
    The MPI parallelism over k-points and spins is efficient at the level of the wall-time
    but it requires HDF5 + MPI-IO support and memory does not scale. Use these additional levels if the memory requirements
    are under control and you need to boost the calculation. Note also that in this case the output results are written to
    different text files, only the SIGEPH.nc file will contains all the k-points and spins.

[[eph_task]] = -4
    The number of bands in the self-energy sum is usually small so it does not make sense to
    parallelize along this dimension. The parallelization over q-points seem to be more efficient than
    the one over perturbations although it introduces some load imbalance because, due to memory reasons,
    the code distributes the q-points in the IBZ (nqibz) instead of the q-points in the full BZ (nqbz).
    Moreover non all the q-points in the IBZ contribute to the imaginary part of $$\Sigma_nk$$.
    The MPI parallelism over k-points and spins is supported with similar behaviour as in **eph_task** +4.


!!! important

    The total number of MPI processes must be equal to the product of the different entries.
""",
),

Variable(
    abivarname="eph_use_ftinterp",
    varset="eph",
    vartype="integer",
    topics=['ElPhonInt_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="EPH FORCE Fourier Transform Interpolation of DFPT potentials.",
    added_in_version="9.0.0",
    text=r"""
This is an *advanced option* used for testing/debugging the interpolation of the DFPT potentials when [[eph_task]] in (2, -2).
By default, the code seeks for the q-point in the input DVDB file when *eph_use_ftinterp* is set to zero (default)
and stops is the q-point in not found in the file.
When *eph_use_ftinterp* is set to 1, the input DVDB file (assumed to contain the [[ddb_ngqpt]] q-mesh)
will be used to generate the real-space representation of the DFPT potentials and interpolate the potential
at the input [[qpt]].
""",
),

#Variable(
#    abivarname="eph_alpha_gmin",
#    varset="eph",
#    vartype="real",
#    topics=['ElPhonInt_expert'],
#    dimensions="scalar",
#    defaultval=0,
#    mnemonics="EPH ALPHA times norm of GMIN.",
#    added_in_version="9.0.0",
#    text=r"""
#This is an *advanced option* used to compute the long-range part of the DFTP potential.
#TO BE DESCRIBED WHEN WE ENTER PRODUCTION
#""",
#),

Variable(
    abivarname="getpot_path",
    varset="files",
    vartype="string",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=None,
    mnemonics="GET the KS POTential from PATH",
    added_in_version="9.0.0",
    text=r"""
This variable defines the path of the POT file containing the KS ground-state potential
that should be used in input.
At present, it is mainly used in the EPH code when performing calculation with the Sternheimer equation.
Note that the path must be inserted between quotation marks.
Note also that relative paths are interpreted according to the working directory in which Abinit is executed!
""",
),

Variable(
    abivarname="getwfk_path",
    varset="files",
    vartype="string",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=None,
    mnemonics="GET the wavefunctions from WFK PATH",
    added_in_version="9.0.0",
    text=r"""
Specify the path of the WFK file using a string instead of the dataset index.
Alternative to [[getwfk]] and [[irdwfk]]. The string must be enclosed between quotation marks:

    getwfk_path "../outdata/out_WFK"
"""
),


Variable(
    abivarname="getwfkfine_path",
    varset="files",
    vartype="string",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=None,
    mnemonics="GET the fine wavefunctions from PATH",
    added_in_version="9.0.0",
    text=r"""
Specify the path of the fine WFK file using a string instead of the dataset index.
Alternative to [[getwfkfine]] and [[irdwfkfine]]. The string must be enclosed between quotation marks:

    getwfkfine_path "../outdata/out_WFK"
"""
),


Variable(
    abivarname="getwfq_path",
    varset="files",
    vartype="string",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=None,
    mnemonics="GET the k+q wavefunctions from WFQ PATH",
    added_in_version="9.0.0",
    text=r"""
Specify the path of the WFQ file using a string instead of the dataset index.
Alternative to [[getwfq]] and [[irdwfq]]. The string must be enclosed between quotation marks:

    getwfq_path "../outdata/out_WFQ"
"""
),

Variable(
    abivarname="getddb_path",
    varset="files",
    vartype="string",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval="None",
    mnemonics="GET the DDB from PATH",
    added_in_version="9.0.0",
    text=r"""
Specify the path of the DDB file using a string instead of the dataset index.
Alternative to [[getddb]] and [[irdddb]]. The string must be enclosed between quotation marks:

    getddb_path "../outdata/out_DDB"
"""
),

Variable(
    abivarname="getdvdb_path",
    varset="files",
    vartype="string",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=None,
    mnemonics="GET the DVDB file from PATH",
    added_in_version="9.0.0",
    text=r"""
Specify the path of the DVDB file using a string instead of the dataset index.
Alternative to [[getdvdb]] and [[irddvdb]]. The string must be enclosed between quotation marks:

    getdvdb_path "../outdata/out_DVDB"
"""
),

Variable(
    abivarname="getden_path",
    varset="files",
    vartype="string",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=None,
    mnemonics="GET the DEN file from PATH",
    added_in_version="9.0.0",
    text=r"""
Specify the path of the DEN file using a string instead of the dataset index.
Alternative to [[getden]] and [[irdden]]. The string must be enclosed between quotation marks:

    getden_path "../outdata/out_DEN"
"""
),

Variable(
    abivarname="getscr_path",
    varset="files",
    vartype="string",
    topics=['multidtset_useful'],
    dimensions="scalar",
    defaultval=None,
    mnemonics="GET the SCR file from PATH",
    added_in_version="9.0.0",
    text=r"""
Specify the path of the SCR file using a string instead of the dataset index.
Alternative to [[getscr]] and [[irdscr]]. The string must be enclosed between quotation marks:

    getscr_path "../outdata/out_SCR"
"""
),

Variable(
    abivarname="eph_ecutosc",
    varset="eph",
    vartype="real",
    topics=['ElPhonInt_expert'],
    dimensions="scalar",
    defaultval="0.0 Hartree",
    mnemonics="Electron-Phonon: Energy CUToff for OSCillator matrix elements",
    characteristics=['[[ENERGY]]'],
    added_in_version="9.0.0",
    text=r"""
This variable defines the energy cutoff defining the number of G-vectors in the oscillator matrix elements:

$$ \langle \mathbf{k+q},b_1 | e^{+i (\mathbf{q+G)} \mathbf{r}} | \mathbf{k}, b_2 \rangle $$

These quantities are used to compute the long-range part of the e-ph matrix elements that are then used
to integrate the Frohlich divergence.

Possible values:

    - = 0 --> Approximate oscillators with $ \delta_{b_1 b_2} $
    - > 0 --> Use full expression with G-dependence
    - < 0 --> Deactivate computation of oscillators.

!!! important

    eph_ecutosc cannot be greater than [[ecut]]

This variable is still under development!
""",
),

Variable(
    abivarname="output_file",
    varset="files",
    vartype="string",
    topics=['Control_useful'],
    dimensions="scalar",
    defaultval=None,
    mnemonics="OUTPUT FILE",
    added_in_version="9.0.0",
    text=r"""
String specifying the name of the main output file
when Abinit is executed with the new syntax:

    abinit run.abi > run.log 2> run.err &

If not specified, the name of the output file is automatically generated by replacing
the file extension of the input file with `.abo`.
To specify the filename in the input use the syntax:

    output_file "t01.out"

with the string enclosed between double quotation marks.
"""
),

Variable(
    abivarname="indata_prefix",
    varset="files",
    vartype="string",
    topics=['Control_useful'],
    dimensions="scalar",
    defaultval=None,
    mnemonics="INput DATA PREFIX",
    added_in_version="9.0.0",
    text=r"""
Prefix for input files. Replaces the analogous entry in the obsolete *files_file* .
This variable is used when Abinit is executed with the new syntax:

    abinit run.abi > run.log 2> run.err &

If this option is not specified, a prefix is automatically constructed from the input file name
provided the file ends with e.g. `.ext`. (`.abi` is recommended)
If the input file does not have a file extension, a default is provided.
"""
),

Variable(
    abivarname="outdata_prefix",
    varset="files",
    vartype="string",
    topics=['Control_useful'],
    dimensions="scalar",
    defaultval=None,
    mnemonics="OUTput DATA PREFIX",
    added_in_version="9.0.0",
    text=r"""
Prefix for output files. Replaces the analogous entry in the obsolete *files_file* .
This variable is used when Abinit is executed with the new syntax:

    abinit run.abi > run.log 2> run.err &

If this option is not specified, a prefix is automatically constructed from the input file name
provided the file ends with e.g. `.ext`. (`.abi` is recommended)
If the input file does not have a file extension, a default is provided.
"""
),

Variable(
    abivarname="tmpdata_prefix",
    varset="files",
    vartype="string",
    topics=['Control_useful'],
    dimensions="scalar",
    defaultval=None,
    mnemonics="TeMPorary DATA PREFIX",
    added_in_version="9.0.0",
    text=r"""
Prefix for temporary files. Replaces the analogous entry in the obsolete *files_file* .
This variable is used when Abinit is executed with the new syntax:

    abinit run.abi > run.log 2> run.err &

If this option is not specified, a prefix is automatically constructed from the input file name
provided the file ends with `.ext`.

If the input file does not have a file extension, a default is provided.
"""
),

Variable(
    abivarname="pp_dirpath",
    varset="files",
    vartype="string",
    topics=['Control_useful'],
    dimensions="scalar",
    defaultval="",
    mnemonics="PseudoPotential DIRectory PATH",
    added_in_version="9.0.0",
    text=r"""
This variable specifies the directory that will prependeded to the names of the pseudopotentials
specified in [[pseudos]].
This option is useful when all your pseudos are gathered in a single directory in your file system
and you don not want to type the absolute path for each pseudopotential file.

This variable is used when Abinit is executed with the new syntax:

    abinit run.abi > run.log 2> run.err &

The string must be quoted in double quotation marks:

    pp_dirpath "/home/user/my_pseudos/"
    pseudos "al.psp8, as.psp8"

If not present, the filenames specified in [[pseudos]] are used directly.
"""
),

Variable(
    abivarname="pseudos",
    varset="files",
    vartype="string",
    topics=['Control_useful'],
    dimensions="scalar",
    defaultval="",
    mnemonics="PSEUDOpotentialS",
    added_in_version="9.0.0",
    text=r"""
String defining the list of pseudopotential files when Abinit is executed with the new syntax:

    abinit run.abi > run.log 2> run.err &

The string must be quoted in double quotation marks and multiple files should be separated by a comma, e.g.

    pseudos = "al.psp8, as.psp8"

The **mandatory** list must contain [[ntypat]] pseudos ordered according to the [[znucl]] array.
The directory where all pseudos are located can be specified with [[pp_dirpath]].

!!! important

    Shell variables e.g. $HOME or tilde syntax `~` for user home are not supported.
"""
),

Variable(
    abivarname="structure",
    varset="basic",
    vartype="string",
    topics=['crystal_useful'],
    dimensions="scalar",
    defaultval="",
    mnemonics="initialize the crystalline STRUCTURE from ...",
    added_in_version="9.0.0",
    text=r"""
This variable provides a simplified interface to build the crystalline structure from an external file.
The idea is to keep the **geometrical information** separated from the input file so that one
can perform multiple calculations in different input files sharing the same structure
without having to copy & paste the description of the unit cell inside the input.
The single source of truth is now given by an external file that can be easily shared.
As a side effect one can easily restart structure relaxations in place by reading
the structure from the output file of a previous run.

The [[structure]] variable is a string in the format **filetype:filepath** where:

* filetype specifies the format of the external file
* filepath gives the path to the file *relative* to the directory where the input file is located.

Variables such as [[natom]], [[ntypat]], [[typat]] and [[znucl]] are automatically initialized from
the external file and need not to be specified in the ABINIT input.

At present ( |today| ), the allowed values for **filetype** are:

* abifile --> An output file produced by Abinit (only netcdf files are supported for the time being)
* abivars --> A txt input file with Abinit variables
* poscar  --> POSCAR files in VASP-5 format (element symbol after the atomic position is required).

Some examples will help clarify.

To read the structure from an external netcdf file produced by Abinit (e.g. *out_GSR.nc*)
use the **abifile** prefix and the syntax:

    structure "abifile:out_GSR.nc"

Other Abinit output files such as the `WFK.nc`, the `DEN.nc` and the `HIST.nc` file
are supported as well e.g.

    structure "abifile:out_HIST.nc"

In the case of structural relaxations, these files contain the final geometry (not necessarily relaxed
within the given tolerance) hence [[structure]] can be used to perform an in-place restart by reading
the output of a previous run.

To read the structure from an external file with the structure in Abinit format, use:

    structure "abivars:my_text_file"

where *my_text_file* specifies the lattice in terms [[acell], ([[rprimd]] or [[angdeg]])
while the atomic positions are specified with [[natom]] and [[xred_symbols]].

```
# MgB2 lattice structure

acell   2*3.086  3.523 Angstrom

rprim   0.866025403784439  0.5  0.0
       -0.866025403784439  0.5  0.0
        0.0                0.0  1.0

natom   3

# Reduced positions followed by element symbol.

xred_symbols
 0.0  0.0  0.0  Mg
 1/3  2/3  0.5  B
 2/3  1/3  0.5  B
```


To read the structure from an external POSCAR file, use:

    structure "poscar:t04_POSCAR"

where `t04_POSCAR` is the name of external file.

!!! critical

    Note that the Abinit parser does not support all the possible variants of the POSCAR format.
    More specifically, we assume a subset of the POSCAR specifications in VASP-5 format.
    This means that the list of element symbols must be present and **written on a single line**.
    Only the "direct" and "cartesian" keywords before the atomic positions are allowed (case insensitive).

A typical POSCAR file for hexagonal MgB2 looks like (ignore comments):

```
Mg1 B2                           # Title string (ignored by Abinit)
1.0                              # Scaling factor for lattice vectors.
2.672554  1.543000 0.000000      # Lattice vectors in Angstrom (NOTE: ABINIT uses BOHR by default)
-2.672554 1.543000 0.000000
0.000000  0.000000 3.523000
Mg B                             # List of element symbols
1 2                              # Number of atoms for each symbol
direct                           # "direct" for reduced coordinates or "cartesian".
0.000000 0.000000 0.0 Mg         # Coordinates followed by the chemical symbol of the atom.
0.333333 0.666667 0.5 B
0.666667 0.333333 0.5 B
```

A positive scaling factor can be used to rescale the lattice vectors
whereas a negative value is interpreted as the volume of the unit cell in Angstrom**3.
Obviously, zero is not allowed!

The [[typat]] variable is automatically initialized from the list of chemical symbols
according to their position in the list.
In this example, Mg is of type 1 while B is of type 2.

The ABINIT variables associated to this POSCAR are therefore:

    ntypat 2
    typat 1 2 2
    znucl  12.0   5.0

These are the **ONLY QUANTITIES** that are initialized from the external POSCAR so please make sure that
your POSCAR resembles the example given above and do not expect ABINIT to understand other entries
such as `Selective dynamics` or velocities.

!!! important

    Several POSCAR files available on the internet give atomic positions and lattice vectors with ~6 digits.
    The ABINIT routines use tight tolerances to detect the space group thus it may happen that ABINIT does not
    detect all the symmetry operations with a consequent **INCREASE** of the number of k-points in the IBZ
    and the associated computational cost. This is especially true for hexagonal or rhombohedral lattices.
    A possible solution is to increase the value of [[tolsym]] in the input file to e.g. 1e-4
    so that ABINIT will automatically refine the atomic positions.

Note the following important remarks:

- The structure is initialized by the parser at the very beginning of the calculation
  hence the external files **must exist** when Abinit starts to analyze the input file.
  In a nutshell, [[structure]] variables cannot be used to pass the output geometry from one dataset
  to the next one.

- Multidatasets are supported but mind that some variables such as
  [[ntypat]], [[typat]] and [[znucl]] are tagged as [[NO_MULTI]].
  In other words, one can read different files via [[structure]] and the multidatset syntax
  provided these quantities do not change.
  ABINIT syntax such as `xred+` are, obviously, not supported.

- The value of [[typat]] and [[znucl]] given in the input file (if any) is ignored by the parser.
  The value of [[natom]], [[ntypat]] is checked for consistency.

- As a rule of thumb, do not try to mix the two approaches: either use [[structure]] or the standard
  (more verbose) approach based on [[ntypat]], [[typat]] and [[znucl]] to define the unit cell.


Limitations:

- The specification of structures for calculations with images is not supported.
- Alchemical mixing is not supported.
- Reading structure from Fortran file is not yet implemented. It is just a technical problem
  that will be hopefully solved in the next releases.

In all these cases in which the [[structure]] variable is not supported,
one has to resort to the **standard approach** to define the list of atoms and their type.
"""
),

]
