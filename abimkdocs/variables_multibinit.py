# coding: utf-8
from __future__ import print_function, division, unicode_literals, absolute_import

executable = "multibinit"

from abimkdocs.variables import ValueWithUnit, MultipleValue, Range
ValueWithConditions = dict

Variable=dict
variables = [
Variable(
    abivarname="dipdip@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['PhononBands_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="DIPole-DIPole interaction",    
    text="""
* 1 --> Recompute the dipole-dipole interaction.
* 0 --> Do not recompute the dipole-dipole interaction.
""",
),

Variable(
    abivarname="dipdip_prt@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['PhononBands_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="DIPole-DIPole PRinT",    
    text="""
* 1 --> Print the dipole-dipole interaction into the XML.
* 0 --> Do not print the dipole-dipole interaction into the XML.
""",
),

Variable(
    abivarname="dipdip_range@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['PhononBands_useful'],
    dimensions=[3],
    defaultval=0,
    mnemonics="Dipole-Dipole interaction",    
    text="""
Depending of the cases, the range of the dipole-dipole interaction will be parameted by:

* dipdip_range if superior to n_cell and superior to short-range interaction
* n_cell if dipdip_range inferior to n_cell
* short-range if dipdip_range inferior to short-range interaction

For example:

    * if dipdip_range = 2 2 2 and the short range interaction if 3 3 3, the dipdip interaction will be set on 3 3 3
    
    * if n_cell = 15 15 15 and the dipdip_range is 6 6 6, the dipdip interaction will be set on 15 15 15
""",
),


Variable(
    abivarname="energy_reference@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['PhononBands_useful'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="Energy of the refences structure",
    characteristics=['[[ENERGY]]'],
    text="""
Set the energy of the reference structure (from the DFT calculation)
if the energy of the reference is not specified in the DDB,
(for example if the DDB file of the ground states is not merged),
or not specified in the XML file), (by default in Hartree).
""",
),

Variable(
    abivarname="ngqpt@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['PhononBands_useful'],
    dimensions=[3],
    defaultval="3*1",
    mnemonics="Number of Grids points for Q PoinTs",
    text="""
The Monkhorst-Pack grid linear dimensions, for the DDB (coarse grid).
""",
),

Variable(
    abivarname="nqshft@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['PhononBands_useful'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Number of Q SHiFTs",
    text="""
The number of vector shifts of the simple Monkhorst and Pack grid, needed to
generate the coarse grid of q points (for the series of fine grids, the number
of shifts it is always taken to be 1). Usually, put it to 1. Use 2 if BCC
sampling (Warning: not BCC lattice, BCC *sampling*), and 4 for FCC sampling
(Warning: not FCC lattice, FCC *sampling*).
""",
),
    
Variable(
    abivarname="prt_model@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['PhononBands_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Effective potential XML output",
    text="""
* 0  -->  do nothing (Default).
* 1  --> Generate the XML file with:

         * The system definition and the model (Harmonic + Anharmonic) in _model.xml

* 2  --> Generate two XML files with:

         * The system definition and the model (Harmonic) in _sys.XML
         * The model (Anharmonic) in _coeffs.xml

* 3  --> Generate only one XML file with:

         * The system definition and the model (Harmonic) in _sys.XML

* 4  --> Generate only one XML file with:

         * The model (Anharmonic) in _coeffs.xml
""",
),


Variable(
    abivarname="fit_coeff@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="FIT anharmonic COEFFficients",
    text="""
* 0  --> do not active the fit process

* 1  --> active the fit process

* -1 --> only for developers, print the files for the scripts
""",
),

Variable(
    abivarname="fit_ncoeff@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="FIT Number of COEFFicients",
    text="""
Give the number of anharmonic coefficients to add in the model during the fit process
""",
),

Variable(
    abivarname="fit_generateCoeff@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="FIT GENERATE anharmonic COEFFicient ",
    text="""
Flag to activate the generation of the anharmonic coefficient for the fit process
""",
),

Variable(
    abivarname="fit_initializeData@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="FIT INITIALIZE DATA for the fit",
    text="""
Flag to de/activate  the precomputing and  storage of all the data for the fit, it will reduce the computation time but increase a lot the memory...
""",
),

Variable(
    abivarname="fit_rangePower@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_basic'],
    dimensions=[2],
    defaultval="3 4",
    mnemonics="FIT RANGE POWER for the coefficients",
    text="""
Set the range of the powers for the anharmonic coefficients
""",
),

Variable(
    abivarname="fit_cutoff@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['FitProcess_basic'],
    dimensions="scalar",
    defaultval="Unit cell",
    mnemonics="FIT CUT-OFF of the anharmonic phonon interaction",
    text="""
Cut-off for the  anharmonic phonon interaction (in Bohr)
""",
),
    
Variable(
    abivarname="fit_anhaStrain@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="FIT ANHARmonic STRAIN coefficients",
    text="""
Flag to activate the anharmonic strain. This option will add coefficients like  (eta^4)
""",
),

Variable(
    abivarname="fit_SPCoupling@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="FIT anharmonic Strain-Phonon COUPLING coefficients",
    text="""
Flag to activate the strain  phonon coupling. This option will add coefficients like  (Sr-Ti)^1(eta^4)
""",
),


Variable(
    abivarname="fit_tolMSDE@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['FitProcess_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="FIT TOLerance on Mean Standard Deviation of the Energy",
    text="""
Tolerance  of the fit based on the Mean Standard Deviation of the Energy in (meV/atm) 
""",
),

Variable(
    abivarname="fit_tolMSDS@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['FitProcess_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="FIT TOLerance on Mean Standard Deviation of the Stresses",
    text="""
Tolerance of the fit based on the Mean Standard Deviation of the Stresses in (eV^2/A^2)  
""",
),

Variable(
    abivarname="fit_tolMSDF@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['FitProcess_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="FIT TOLerance on Mean Standard Deviation of the Forces",
    text="""
Tolerance of the fit based on the Mean Standard Deviation of the  Forces (eV^2/A^2)
""",
),

Variable(
    abivarname="fit_tolMSDFS@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['FitProcess_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="FIT TOLerance on Mean Standard Deviation of the Forces and Stresses",
    text="""
Tolerance of the fit based on the Mean Standard Deviation of the  Forces and Sressses (eV^2/A^2)
""",
),

Variable(
    abivarname="fit_nfixcoeff@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="FIT Number of FIXed COEFFicients",
    text="""
Number of imposed coefficients during the fit process for the model: 

* -1 -->  fix all the coefficients

* 0  -->  do not fix coefficients 

* n  -->  fix n coefficients (requires [[multibinit:fit_fixcoeff]] input variable) 
""",
),

Variable(
    abivarname="fit_fixcoeff@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_expert'],
    dimensions="fit_nfixcoeff",
    defaultval=0,
    mnemonics="FIT FIXed COEFFicients",
    text="""
Indexes of the imposed coefficients during the fit process for the model: 
""",
),


Variable(
    abivarname="fit_nbancoeff@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="FIT Number of BANed COEFFicients",
    text="""
Number of imposed coefficients during the fit process of the model: 

* 0  -->  do not ban coefficients 

* n  -->  ban n coefficients (requires [[multibinit:fit_bancoeff]] input variable) 
""",
),

Variable(
    abivarname="fit_bancoeff@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_expert'],
    dimensions="fit_nbancoeff",
    defaultval=0,
    mnemonics="FIT BANed COEFFicients",
    text="""
Indexes of the banned coefficients during the fit process of the model: 
""",
),

Variable(
    abivarname="fit_ts_option@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="FIT Training Set OPTION",
    text="""
* 0 --> the Training is hist from ABINIT 

* 1 --> the Training contains -1 * stress  (usualy output from VASP)
""",
),

Variable(
    abivarname="bound_coeff@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['BoundProcess_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="BOUND COEFFicient",
    text="""
flag to activate the bound process:

* 0 --> no bound process 

* 1 --> bound process 

* 2 -->  new version of the bound process 
""",
),

Variable(
    abivarname="bound_maxCoeff@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['BoundProcess_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="BOUND MAX COEFFicient",
    text="""
Number of maximum additional coefficients for the bound process
""",
),

Variable(
    abivarname="bound_rangePower@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['BoundProcess_basic'],
    dimensions=[2],
    defaultval="6,6",
    mnemonics="BOUND RANGE POWER",
    text="""
range of the power for the additional coefficients in the bound process
""",
),

Variable(
    abivarname="bound_cutoff@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['BoundProcess_useful'],
    dimensions="scalar",
    defaultval="1 unit cell",
    mnemonics="BOUND CUT OFF",
    text="""
Cut-off for the anharmonic phonon interaction in the bound process (in Bohr)
""",
),


Variable(
    abivarname="bound_anhaStrain@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['BoundProcess_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="BOUND ANHArmonic STRAIN coefficients",
    text="""
Flag to activate the anharmonic strain. This option will add terms like (eta^4)
""",
),

Variable(
    abivarname="bound_SPCoupling@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['BoundProcess_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="BOUND Strain Phonon COUPLING coefficients",
    text="""
Flag to activate the strain phonon coupling. This option will add term like (Sr-Ti)^1eta^2
""",
),

Variable(
    abivarname="bound_cell@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['BoundProcess_expert'],
    dimensions=[3],
    defaultval="6,6,6",
    mnemonics="BOUND superCELL size",
    text="""
Set the size of the supercell during the bound process 
""",
),

Variable(
    abivarname="bound_temp@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['BoundProcess_expert'],
    dimensions="scalar",
    defaultval=500,
    mnemonics="BOUND",
    text="""
Set the temperature of the MD during the bound process 
""",
),

Variable(
    abivarname="bound_step@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['BoundProcess_expert'],
    dimensions="scalar",
    defaultval=1000,
    mnemonics="BOUND",
    text="""
Set the maximum number of MD step during the bound process
""",
),

]
