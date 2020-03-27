# coding: utf-8
from __future__ import print_function, division, unicode_literals, absolute_import

executable = "tdep"

from abimkdocs.variables import ValueWithUnit, MultipleValue, Range
ValueWithConditions = dict

Variable=dict
variables = [

Variable(
    abivarname="amu@tdep",
    varset="tdep",
    vartype="real",
    topics=['Tdep_basic'],
    dimensions=['[[tdep:ntypat]]'],
    defaultval="[[tdep:ntypat]]*0.d0",
    mnemonics="Atomic masses in Mass Units",
    added_in_version="before_v9",
    text=r"""
Defines the masses in atomic mass units for each kind of atom. See the ABINIT variable [[amu]] for more details. (Only required when the NetCDF file is absent).
""",
),

Variable(
    abivarname="angle@tdep",
    varset="tdep",
    vartype="real",
    topics=['Tdep_basic'],
    dimensions="scalar",
    defaultval="90.d0",
    mnemonics="ANGLE alpha",
    added_in_version="before_v9",
    text=r"""
This angle has to be defined if the bravais lattice is monoclinic. That is to say if [[tdep:brav]](1)=2. 
""",
),

Variable(
    abivarname="brav@tdep",
    varset="tdep",
    vartype="integer",
    topics=['Tdep_basic'],
    dimensions=[2],
    defaultval="2*0",
    mnemonics="BRAVais",
    added_in_version="before_v9",
    text="""
These two parameters define the Bravais lattice (as defined in the ABINIT code) and the primitive vectors [[rprim]] in the TDEP code.

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
    abivarname="bzpath@tdep",
    varset="tdep",
    vartype="integer+letter",
    topics=['Tdep_expert'],
    dimensions="'[[tdep:bzpath]](1)'+1",
    defaultval="0",
    mnemonics="Brillouin Zone PATH",
    added_in_version="before_v9",
    text=r"""
OPTIONAL: Defines the path in the Brillouin Zone for the phonon spectrum calculation. The first value defines the number of special points used in the path. The other values define the special points of the BZ (only the letters fixed by convention for the present lattice are allowed: L, X, M... and G for $\Gamma$). 
""",
),

Variable(
    abivarname="dosdeltae@tdep",
    varset="tdep",
    vartype="integer",
    topics=['Tdep_expert'],
    dimensions="scalar",
    defaultval="4.5d-6",
    mnemonics="DOS delta Energy",
    added_in_version="before_v9",
    text="""
OPTIONAL: Defines the smearing used for the phonon Density Of State calculation.
""",
),

Variable(
    abivarname="enunit@tdep",
    varset="tdep",
    vartype="integer",
    topics=['Tdep_expert'],
    dimensions="scalar",
    defaultval="0",
    mnemonics="ENergy UNIT",
    added_in_version="before_v9",
    text="""
OPTIONAL: Defines the energy unit used for the phonon spectrum (0 for meV, 1 for cm-1 and 2 for Ha).
""",
),

Variable(
    abivarname="multiplicity@tdep",
    varset="tdep",
    vartype="real",
    topics=['Tdep_basic'],
    dimensions=[3,3],
    defaultval="9*0.d0",
    mnemonics="MULTIPLICITY",
    added_in_version="before_v9",
    text=r"""
Defines the multiplicity of the SUPERCELL with respect to the primitive UNICELL. See the ABINIT variables [[rprimd]], [[acell]] and [[rprim]] for more details. The multiplicity [[tdep:multiplicity]] and the SUPERCELL lattice parameters [[tdep:rprimd]] are used to find the UNITCELL lattice parameters acell_unitcell such as:

$$ \text{rprimd}_{i,j}=\sum_{k=1}^3 \text{acell_unitcell}_i * \text{multiplicity}_{i,k}*\text{rprim_tmp}_{k,j} $$

For example:

- for a fcc lattice: rprim = ( 0 1/2 1/2 ; 1/2 0 1/2 ; 1/2 1/2 0) and acell = (a a a). If the SUPERCELL is rprimd = (3a 0 0 ; 0 3a 0 ; 0 0 3a), the multiplicity = ( -3 3 3 ; 3 -3 3 ; 3 3 -3) 
- for a bcc lattice: rprim = ( -1/2 1/2 1/2 ; 1/2 -1/2 1/2 ; 1/2 1/2 -1/2) and acell = (a a a). If the SUPERCELL is rprimd = (3a 0 0 ; 0 3a 0 ; 0 0 3a), the multiplicity = ( 0 3 3 ; 3 0 3 ; 3 3 0) 
""",
),

Variable(
    abivarname="natom@tdep",
    varset="tdep",
    vartype="integer",
    topics=['Tdep_basic'],
    dimensions="scalar",
    defaultval="0",
    mnemonics="NATOM",
    added_in_version="before_v9",
    text="""
Defines the number of atoms in the SUPERCELL. See the ABINIT variable [[natom]] for more details. (Only required when the NetCDF file is absent).
""",
),

Variable(
    abivarname="natom_unitcell@tdep",
    varset="tdep",
    vartype="integer",
    topics=['Tdep_basic'],
    dimensions="scalar",
    defaultval="0",
    mnemonics="NATOM in the UNITCELL",
    added_in_version="before_v9",
    text="""
Defines the number of atoms in the UNITCELL.
""",
),

Variable(
    abivarname="ngqpt1@tdep",
    varset="tdep",
    vartype="integer",
    topics=['Tdep_expert'],
    dimensions=[3],
    defaultval=[8, 8, 8],
    mnemonics="Number of Grid points for Q PoinTs generation (coarse)",
    added_in_version="before_v9",
    text=r"""
OPTIONAL: Defines the COARSE grid of q-points for the dynamical matrix output (in DDB).
""",
),

Variable(
    abivarname="ngqpt2@tdep",
    varset="tdep",
    vartype="integer",
    topics=['Tdep_expert'],
    dimensions=[3],
    defaultval=[32, 32, 32],
    mnemonics="Number of Grid points for Q PoinTs generation (fine)",
    added_in_version="before_v9",
    text="""
OPTIONAL: Defines the FINE grid of q-points for the DOS and thermodynamic quantity calculations.
""",
),

Variable(
    abivarname="nstep_max@tdep",
    varset="tdep",
    vartype="integer",
    topics=['Tdep_basic'],
    dimensions="scalar",
    defaultval="0",
    mnemonics="NSTEP at MAX",
    added_in_version="before_v9",
    text="""
Defines the upper limit in the range of configurations that one wants to use. This number has to be lower than the maximum number of configurations present in the NetCDF or ASCII file.
""",
),

Variable(
    abivarname="nstep_min@tdep",
    varset="tdep",
    vartype="integer",
    topics=['Tdep_basic'],
    dimensions="scalar",
    defaultval="0",
    mnemonics="NSTEP at MIN",
    added_in_version="before_v9",
    text="""
Defines the lower limit in the range of configurations that one wants to use. This number has to be larger than the minimum number of configurations present in the NetCDF or ASCII file.
""",
),

Variable(
    abivarname="ntypat@tdep",
    varset="tdep",
    vartype="integer",
    topics=['Tdep_basic'],
    dimensions="scalar",
    defaultval="0",
    mnemonics="NTYPAT",
    added_in_version="before_v9",
    text="""
Defines the number of atom types. See the ABINIT variable [[ntypat]] for more details. (Only required when the NetCDF file is absent).
""",
),

Variable(
    abivarname="order@tdep",
    varset="tdep",
    vartype="integer+real",
    topics=['Tdep_expert'],
    dimensions="2",
    defaultval="2",
    mnemonics="ORDER for the IFC",
    added_in_version="before_v9",
    text="""
OPTIONAL: Defines at which order the calculation of the IFCs is performed. If the first value [[tdep:order]](1)=3, that turns on a third order calculation and the second value [[tdep:order]](2) defines the cutoff radius. 
""",
),

Variable(
    abivarname="rcut@tdep",
    varset="tdep",
    vartype="real",
    topics=['Tdep_basic'],
    dimensions="scalar",
    defaultval="0.d0",
    mnemonics="Radius CUToff",
    added_in_version="before_v9",
    text="""
Defines the cutoff radius used when the second order IFCs are computed. This ones has to be lower than half the smallest SUPERCELL lattice parameter.
""",
),

Variable(
    abivarname="rprimd@tdep",
    varset="tdep",
    vartype="real",
    topics=['Tdep_basic'],
    dimensions=[3,3],
    defaultval="9*0.d0",
    mnemonics="RPRIMD",
    added_in_version="before_v9",
    text="""
Defines the dimensional real space primitive vectors of the SUPERCELL. See [[rprimd]] for more details. (Only required when the NetCDF file is absent).
""",
),

Variable(
    abivarname="slice@tdep",
    varset="tdep",
    vartype="integer",
    topics=['Tdep_expert'],
    dimensions="1",
    defaultval="1",
    mnemonics="SLICE",
    added_in_version="before_v9",
    text="""
OPTIONAL: Defines the slice used to include some configurations in the calculations. Only the ([[tdep:nstep_max]]-[[tdep:nstep_min]])/[[tdep:slice]] configurations will be considered in the calculations of the IFCs.
""",
),


Variable(
    abivarname="temperature@tdep",
    varset="tdep",
    vartype="real",
    topics=['Tdep_basic'],
    dimensions="scalar",
    defaultval="0",
    mnemonics="TEMPERATURE",
    added_in_version="before_v9",
    text="""
Defines the temperature of the system.
""",
),

Variable(
    abivarname="typat@tdep",
    varset="tdep",
    vartype="integer",
    topics=['Tdep_basic'],
    dimensions=['[[tdep:natom]]'],
    defaultval="[[tdep:natom]]*0",
    mnemonics="TYPAT",
    added_in_version="before_v9",
    text="""
Defines the type of atoms in the SUPERCELL. See [[typat]] for more details. (Only required when the NetCDF file is absent).
""",
),

Variable(
    abivarname="typat_unitcell@tdep",
    varset="tdep",
    vartype="integer",
    topics=['Tdep_basic'],
    dimensions=['[[tdep:natom_unitcell]]'],
    defaultval="[[tdep:natom_unitcell]]*0",
    mnemonics="TYPAT in the UNITCELL",
    added_in_version="before_v9",
    text="""
Defines the type of atoms in the UNITCELL.
""",
),

Variable(
    abivarname="use_ideal_positions@tdep",
    varset="tdep",
    vartype="integer",
    topics=['Tdep_expert'],
    dimensions="scalar",
    defaultval="1",
    mnemonics="USE IDEAL POSITIONS",
    added_in_version="before_v9",
    text="""
OPTIONAL: Defines if the ideal ([[tdep:use_ideal_positions]]=1) or averaged ([[tdep:use_ideal_positions]]=0) positions are used during the calculations. It can affect strongly the phonon spectrum (and other quantities) if the system is close to an instability (soft mode,...).
""",
),

Variable(
    abivarname="xred_unitcell@tdep",
    varset="tdep",
    vartype="real",
    topics=['Tdep_basic'],
    dimensions=[3, '[[tdep:natom_unitcell]]'],
    defaultval="(3*[[tdep:natom_unitcell]])*0.d0",
    mnemonics="XRED in the UNITCELL",
    added_in_version="before_v9",
    text="""
Defines the reduced coordinates of atoms in the UNITCELL.
""",
),

]
