# coding: utf-8
from __future__ import print_function, division, unicode_literals, absolute_import

executable = "atdep"

from abimkdocs.variables import ValueWithUnit, MultipleValue, Range
ValueWithConditions = dict

Variable=dict
variables = [

Variable(
    abivarname="alloy@atdep",
    varset="atdep",
    vartype="integer",
    topics=['aTDEP_expert'],
    dimensions="3",
    defaultval="0",
    mnemonics="ALLOY treatment",
    added_in_version="9.5.1",
    text=r"""
OPTIONAL: Defines the treatment of the alloy. The first value defines the approximation used (0=nothing and 1=Virtual Crystal Approximation). The second and the third values define the [[typat_unitcell] of the atoms which have to be alloyed.
""",
),

Variable(
    abivarname="amu@atdep",
    varset="atdep",
    vartype="real",
    topics=['aTDEP_basic'],
    dimensions=['[[atdep:ntypat]]'],
    defaultval="[[atdep:ntypat]]*0.d0",
    mnemonics="Atomic masses in Mass Units",
    added_in_version="before_v9",
    text=r"""
Defines the masses in atomic mass units for each kind of atom. See the ABINIT variable [[amu]] for more details. (Only required when the NetCDF file is absent).
""",
),

Variable(
    abivarname="angle@atdep",
    varset="atdep",
    vartype="real",
    topics=['aTDEP_basic'],
    dimensions="scalar",
    defaultval="90.d0",
    mnemonics="ANGLE alpha",
    added_in_version="before_v9",
    text=r"""
This angle has to be defined if the bravais lattice is monoclinic. That is to say if [[atdep:brav]](1)=2. 
""",
),

Variable(
    abivarname="born_charge@atdep",
    varset="atdep",
    vartype="real",
    topics=['aTDEP_expert'],
    dimensions=['[[atdep:ntypat]]'],
    defaultval="[[atdep:ntypat]]*0.d0",
    mnemonics="BORN effective CHARGE",
    added_in_version="before_v9",
    text=r"""
OPTIONAL : Defines the Born effective charge (for each kind of atom) used to compute the dipole-dipole interaction.
""",
),

Variable(
    abivarname="brav@atdep",
    varset="atdep",
    vartype="integer",
    topics=['aTDEP_basic'],
    dimensions=[2],
    defaultval="2*0",
    mnemonics="BRAVais",
    added_in_version="before_v9",
    text="""
These two parameters define the Bravais lattice (as defined in the ABINIT code) and the primitive vectors [[rprim]] in the aTDEP code.

- For bravais(1): The holohedral groups are numbered as follows (see international tables for crystallography (1983), p. 13):
iholohedry=1   triclinic      1bar
iholohedry=2   monoclinic     2/m
iholohedry=3   orthorhombic   mmm
iholohedry=4   tetragonal     4/mmm
iholohedry=5   trigonal       3bar m
iholohedry=6   hexagonal      6/mmm
iholohedry=7   cubic          m3bar m

- For bravais(2): The centering is defined as follows:
center=0        no centering
center=-1       body-centered
center=-3       face-centered
center=1        A-face centered
center=2        B-face centered
center=3        C-face centered

""",
),

Variable(
    abivarname="bzlength@atdep",
    varset="atdep",
    vartype="integer+real",
    topics=['aTDEP_expert'],
    dimensions="'[[atdep:bzlength]](1)'+1",
    defaultval="0",
    mnemonics="Brillouin Zone LENGTH",
    added_in_version="before_v9",
    text=r"""
OPTIONAL: Defines the length of the Brillouin Zone for the phonon spectrum calculation. The first value defines the number of segments used in the path. The other values define the size of each segment.
""",
),

Variable(
    abivarname="bzpath@atdep",
    varset="atdep",
    vartype="integer+letter",
    topics=['aTDEP_expert'],
    dimensions="'[[atdep:bzpath]](1)'+1",
    defaultval="0",
    mnemonics="Brillouin Zone PATH",
    added_in_version="before_v9",
    text=r"""
OPTIONAL: Defines the path in the Brillouin Zone for the phonon spectrum calculation. The first value defines the number of special points used in the path. The other values define the special points of the BZ (only the letters fixed by convention for the present lattice are allowed: L, X, M... and G for $\Gamma$). 
""",
),

Variable(
    abivarname="dielec_constant@atdep",
    varset="atdep",
    vartype="real",
    topics=['aTDEP_expert'],
    dimensions="scalar",
    defaultval="0",
    mnemonics="DIELECtric CONSTANT",
    added_in_version="before_v9",
    text="""
OPTIONAL: Defines the dielectric constant used to compute the dipole-dipole interaction.
""",
),

Variable(
    abivarname="dosdeltae@atdep",
    varset="atdep",
    vartype="integer",
    topics=['aTDEP_expert'],
    dimensions="scalar",
    defaultval="0.2 cm$^{-1}$",
    mnemonics="DOS delta Energy",
    added_in_version="before_v9",
    text="""
OPTIONAL: Defines the smearing used for the phonon Density Of State calculation.

Prior to v9.10, the default was 4.5d-6.
""",
),

Variable(
    abivarname="enunit@atdep",
    varset="atdep",
    vartype="integer",
    topics=['aTDEP_expert'],
    dimensions="scalar",
    defaultval="0",
    mnemonics="ENergy UNIT",
    added_in_version="before_v9",
    text="""
OPTIONAL: Defines the energy unit used for the phonon spectrum (0 for meV, 1 for cm-1, 2 for Ha and 3 for THz).
""",
),

Variable(
    abivarname="multiplicity@atdep",
    varset="atdep",
    vartype="real",
    topics=['aTDEP_basic'],
    dimensions=[3,3],
    defaultval="9*0.d0",
    mnemonics="MULTIPLICITY",
    added_in_version="before_v9",
    text=r"""
Defines the multiplicity of the SUPERCELL with respect to the primitive UNICELL. See the ABINIT variables [[rprimd]], [[acell]] and [[rprim]] for more details. The multiplicity [[atdep:multiplicity]] and the SUPERCELL lattice parameters [[atdep:rprimd]] are used to find the UNITCELL lattice parameters acell_unitcell such as:

$$ \text{rprimd}_{i,j}=\sum_{k=1}^3 \text{acell_unitcell}_i * \text{multiplicity}_{i,k}*\text{rprim_tmp}_{k,j} $$

For example:

- for a fcc lattice: rprim = ( 0 1/2 1/2 ; 1/2 0 1/2 ; 1/2 1/2 0) and acell = (a a a). If the SUPERCELL is rprimd = (3a 0 0 ; 0 3a 0 ; 0 0 3a), the multiplicity = ( -3 3 3 ; 3 -3 3 ; 3 3 -3) 
- for a bcc lattice: rprim = ( -1/2 1/2 1/2 ; 1/2 -1/2 1/2 ; 1/2 1/2 -1/2) and acell = (a a a). If the SUPERCELL is rprimd = (3a 0 0 ; 0 3a 0 ; 0 0 3a), the multiplicity = ( 0 3 3 ; 3 0 3 ; 3 3 0) 
""",
),

Variable(
    abivarname="natom@atdep",
    varset="atdep",
    vartype="integer",
    topics=['aTDEP_basic'],
    dimensions="scalar",
    defaultval="0",
    mnemonics="NATOM",
    added_in_version="before_v9",
    text="""
Defines the number of atoms in the SUPERCELL. See the ABINIT variable [[natom]] for more details. (Only required when the NetCDF file is absent).
""",
),

Variable(
    abivarname="natom_unitcell@atdep",
    varset="atdep",
    vartype="integer",
    topics=['aTDEP_basic'],
    dimensions="scalar",
    defaultval="0",
    mnemonics="NATOM in the UNITCELL",
    added_in_version="before_v9",
    text="""
Defines the number of atoms in the UNITCELL.
""",
),

Variable(
    abivarname="ngqpt1@atdep",
    varset="atdep",
    vartype="integer",
    topics=['aTDEP_expert'],
    dimensions=[3],
    defaultval=[8, 8, 8],
    mnemonics="Number of Grid points for Q PoinTs generation (coarse)",
    added_in_version="before_v9",
    text=r"""
OPTIONAL: Defines the COARSE grid of q-points for the dynamical matrix output (in DDB).
""",
),

Variable(
    abivarname="ngqpt2@atdep",
    varset="atdep",
    vartype="integer",
    topics=['aTDEP_expert'],
    dimensions=[3],
    defaultval=[32, 32, 32],
    mnemonics="Number of Grid points for Q PoinTs generation (fine)",
    added_in_version="before_v9",
    text="""
OPTIONAL: Defines the FINE grid of q-points for the DOS and thermodynamic quantity calculations.
""",
),

Variable(
    abivarname="nstep_max@atdep",
    varset="atdep",
    vartype="integer",
    topics=['aTDEP_basic'],
    dimensions="scalar",
    defaultval="0",
    mnemonics="NSTEP at MAX",
    added_in_version="before_v9",
    text="""
Defines the upper limit in the range of configurations that one wants to use. This number has to be lower than the maximum number of configurations present in the NetCDF or ASCII file.
""",
),

Variable(
    abivarname="nstep_min@atdep",
    varset="atdep",
    vartype="integer",
    topics=['aTDEP_basic'],
    dimensions="scalar",
    defaultval="0",
    mnemonics="NSTEP at MIN",
    added_in_version="before_v9",
    text="""
Defines the lower limit in the range of configurations that one wants to use. This number has to be larger than the minimum number of configurations present in the NetCDF or ASCII file.
""",
),

Variable(
    abivarname="ntypat@atdep",
    varset="atdep",
    vartype="integer",
    topics=['aTDEP_basic'],
    dimensions="scalar",
    defaultval="0",
    mnemonics="NTYPAT",
    added_in_version="before_v9",
    text="""
Defines the number of atom types. See the ABINIT variable [[ntypat]] for more details. (Only required when the NetCDF file is absent).
""",
),

Variable(
    abivarname="order@atdep",
    varset="atdep",
    vartype="integer+real",
    topics=['aTDEP_expert'],
    dimensions="2",
    defaultval="2",
    mnemonics="ORDER for the IFC",
    added_in_version="before_v9",
    text="""
OPTIONAL: Defines at which order the calculation of the IFCs is performed. If the first value [[atdep:order]](1)=3, that turns on a third order calculation and the second value [[atdep:order]](2) defines the cutoff radius. 
""",
),

Variable(
    abivarname="rcut@atdep",
    varset="atdep",
    vartype="real",
    topics=['aTDEP_basic'],
    dimensions="scalar",
    defaultval="0.d0",
    mnemonics="Radius CUToff",
    added_in_version="before_v9",
    text="""
Defines the cutoff radius used when the second order IFCs are computed. This ones has to be lower than half the smallest SUPERCELL lattice parameter.
""",
),

Variable(
    abivarname="readifc@atdep",
    varset="atdep",
    vartype="integer",
    topics=['aTDEP_expert'],
    dimensions="scalar",
    defaultval="0",
    mnemonics="READ the Interatomic Force Constants",
    added_in_version="before_v9",
    text="""
OPTIONAL: Defines the IO strategy used for the IFC. If :

- [[atdep:readifc]] = 1 : Read the IFC coming from an input file.
- [[atdep:readifc]] = 2 : Write and read the IFC coming from the calculation (for tests).
""",
),

Variable(
    abivarname="rprimd@atdep",
    varset="atdep",
    vartype="real",
    topics=['aTDEP_basic'],
    dimensions=[3,3],
    defaultval="9*0.d0",
    mnemonics="RPRIMD",
    added_in_version="before_v9",
    text="""
Defines the dimensional real space primitive vectors of the SUPERCELL. See [[rprimd]] for more details. (Only required when the NetCDF file is absent).
""",
),

Variable(
    abivarname="slice@atdep",
    varset="atdep",
    vartype="integer",
    topics=['aTDEP_expert'],
    dimensions="1",
    defaultval="1",
    mnemonics="SLICE",
    added_in_version="before_v9",
    text="""
OPTIONAL: Defines the slice used to include some configurations in the calculations. Only the ([[atdep:nstep_max]]-[[atdep:nstep_min]])/[[atdep:slice]] configurations will be considered in the calculations of the IFCs.
""",
),


Variable(
    abivarname="temperature@atdep",
    varset="atdep",
    vartype="real",
    topics=['aTDEP_basic'],
    dimensions="scalar",
    defaultval="0",
    mnemonics="TEMPERATURE",
    added_in_version="before_v9",
    text="""
Defines the temperature of the system.
""",
),

Variable(
    abivarname="together@atdep",
    varset="atdep",
    vartype="integer",
    topics=['aTDEP_expert'],
    dimensions="scalar",
    defaultval="1",
    mnemonics="Are the different orders solved TOGETHER?",
    added_in_version="9.5.1",
    text="""
OPTIONAL: Defines if the different [[atdep:order]] are solved together or not :

- simultaneously : [[atdep:together]] = 1.
- successively : [[atdep:together]] = 0.
""",
),

Variable(
    abivarname="typat@atdep",
    varset="atdep",
    vartype="integer",
    topics=['aTDEP_basic'],
    dimensions=['[[atdep:natom]]'],
    defaultval="[[atdep:natom]]*0",
    mnemonics="TYPAT",
    added_in_version="before_v9",
    text="""
Defines the type of atoms in the SUPERCELL. See [[typat]] for more details. (Only required when the NetCDF file is absent).
""",
),

Variable(
    abivarname="typat_unitcell@atdep",
    varset="atdep",
    vartype="integer",
    topics=['aTDEP_basic'],
    dimensions=['[[atdep:natom_unitcell]]'],
    defaultval="[[atdep:natom_unitcell]]*0",
    mnemonics="TYPAT in the UNITCELL",
    added_in_version="before_v9",
    text="""
Defines the type of atoms in the UNITCELL.
""",
),

Variable(
    abivarname="use_ideal_positions@atdep",
    varset="atdep",
    vartype="integer",
    topics=['aTDEP_expert'],
    dimensions="scalar",
    defaultval="1",
    mnemonics="USE IDEAL POSITIONS",
    added_in_version="before_v9",
    text="""
OPTIONAL: Defines if the ideal ([[atdep:use_ideal_positions]]=1) or averaged ([[atdep:use_ideal_positions]]=0) positions are used during the calculations. It can affect strongly the phonon spectrum (and other quantities) if the system is close to an instability (soft mode,...).
""",
),

Variable(
    abivarname="use_weights@atdep",
    varset="atdep",
    vartype="integer",
    topics=['aTDEP_expert'],
    dimensions="scalar",
    defaultval="0",
    mnemonics="USE the WEIGHTS",
    added_in_version="9.5.1",
    text="""
OPTIONAL: Defines if a specific weight has to be used for each configuration ([[atdep:use_weights]]=1).
""",
),

Variable(
    abivarname="xred_unitcell@atdep",
    varset="atdep",
    vartype="real",
    topics=['aTDEP_basic'],
    dimensions=[3, '[[atdep:natom_unitcell]]'],
    defaultval="(3*[[atdep:natom_unitcell]])*0.d0",
    mnemonics="XRED in the UNITCELL",
    added_in_version="before_v9",
    text="""
Defines the reduced coordinates of atoms in the UNITCELL.
""",
),

]
