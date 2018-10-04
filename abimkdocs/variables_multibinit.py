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
    topics=['LatticeModel_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="DIPole-DIPole interaction",    
    text="""
* 0 --> Do not recompute the dipole-dipole interaction.
* 1 --> Recompute the dipole-dipole interaction based on ewald summation .
""",
),

Variable(
    abivarname="dipdip_prt@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['LatticeModel_expert'],
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
    topics=['LatticeModel_expert'],
    dimensions=[3],
    defaultval=0,
    mnemonics="Dipole-Dipole interaction",    
    text="""
Depending of the cases, the range of the dipole-dipole interaction will be parameted by:

* dipdip_range if superior to ncell and superior to short-range interaction
* ncell if dipdip_range inferior to ncell
* short-range if dipdip_range inferior to short-range interaction

For example:

    * if dipdip_range = 2 2 2 and the short range interaction if 3 3 3, the dipdip interaction will be set on 3 3 3
    
    * if ncell = 15 15 15 and the dipdip_range is 6 6 6, the dipdip interaction will be set on 15 15 15
""",
),


Variable(
    abivarname="energy_reference@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['LatticeModel_useful'],
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
    abivarname="ncoeff@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['LatticeModel_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of anharmonic COEFFicients",
    text="""
Set the number of anharmonic coefficients in the model. This number have to be in agreement with the number of coefficients present in the XML file.

If ncoeff /= 0, [[multibinit:coefficients]] have to be present in the input files
""",
),

Variable(
    abivarname="coefficients@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['LatticeModel_useful'],
    dimensions=['[[multibinit:ncoeff]]'],
    defaultval=0.0,
    mnemonics="values of the COEFFICIENTS",
    text="""
Set the values of the coefficients present in the XML file
""",
),

    
Variable(
    abivarname="ngqpt@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['LatticeModel_useful'],
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
    topics=['LatticeModel_useful'],
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
    topics=['LatticeModel_basic'],
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

* 1 -->  Activate the fit process. This option will first generate a set of coefficients if [[multibinit:fit_generateCoeff]] is set to one. This generation is mainly parametrized by [[multibinit:fit_rangePower]] and [[multibinit:fit_cutoff]]. You can also provided a list of coefficients with the  model_anharmonic.MXL (see [[help:multibinit]]). Then the fit process will select the coefficients one by one up to [[multibinit:fit_ncoeff]] (see this [[cite:Escorihuela-Sayalero2017|paper]] for the details of the procedure). 

* -1 --> **only for developers**, print the files for the scripts
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

**Related variables:** The  power range of the coefficients ([[multibinit:fit_rangePower]]), the cut off of the interactions ([[multibinit:fit_cutoff]]), the flag to add ahnarmonic strain ([[multibinit:fit_anhaStrain]]), the flag to add phonon strain coupling ([[multibinit:fit_SPCoupling]])
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
Flag to activate the strain  phonon coupling. This option will add coefficients like  (Sr-Ti)^1 (eta^4)
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
    dimensions=['[[multibinit:fit_nfixcoeff]]'],
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
    dimensions=['[[multibinit:fit_nbancoeff]]'],
    defaultval=0,
    mnemonics="FIT BANed COEFFicients",
    text="""
Indexes of the banned coefficients during the fit process of the model
""",
),

Variable(
    abivarname="ts_option@multibinit",
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
    abivarname="bound_model@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['BoundingProcess_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="BOUND COEFFicient",
    text="""
Flag to activate the bound process:

* 0 --> Do not activate the bound process 

* 1 --> This option will generate all the possible combinaisons of coefficients from 1 to [[multibinit:bound_maxCoeff]]. Some constrains are imposed during the generation and the fit of the coefficients, they have to be positive and with even power. Finaly, the code will try all the possible combinaisons and try to find a bounded model.

* 2 -->  **new version** This option will generate a set of coefficients with a power range defined by [[multibinit:bound_rangePower]] and keep only the coefficients with even power. Then the procedure is similar to the fit process with the constrains to only keep positive coefficients. The bound process will select the coefficients one by one up to [[multibinit:bound_maxCoeff]] and try if the model is bound at each step of the process. 

**Related variables:** The number of maximum additional coefficient in the polynome ([[multibinit:bound_maxCoeff]]), the  power range for the additional coefficients ([[multibinit:bound_rangePower]]), the cut off of the additional interactions ([[multibinit:bound_cutoff]])
""",
),
     
Variable(
    abivarname="bound_maxCoeff@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['BoundingProcess_basic'],
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
    topics=['BoundingProcess_basic'],
    dimensions=[2],
    defaultval="6,6",
    mnemonics="BOUND RANGE POWER",
    text="""
Range of the power for the additional coefficients in the bound process
""",
),

Variable(
    abivarname="bound_cutoff@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['BoundingProcess_basic'],
    dimensions="scalar",
    defaultval="1 unit cell",
    mnemonics="BOUND CUT OFF",
    text="""
Cut-off for the anharmonic phonon interaction during the bound process (in Bohr)
""",
),


Variable(
    abivarname="bound_anhaStrain@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['BoundingProcess_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="BOUND ANHArmonic STRAIN coefficients",
    text="""
Flag to activate the anharmonic strain. When the bound process will generate the possible coefficients for the fit, if this input variable is set to 1, the generator will consider coefficients like eta^4
""",
),

Variable(
    abivarname="bound_SPCoupling@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['BoundingProcess_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="BOUND Strain Phonon COUPLING coefficients",
    text="""
Flag to activate the strain phonon coupling. When the bound process will generate the possible coefficients for the fit, if this input variable is set to 1, the generator will consider coefficients like (Sr-Ti)^2 eta^2
""",
),

Variable(
    abivarname="bound_cell@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['BoundingProcess_expert'],
    dimensions=[3],
    defaultval="6,6,6",
    mnemonics="BOUND superCELL size for the molecular dynamics",
    text="""
When the process will try a given model, this input variable is used to  set the size of the supercell for the molecular dynamics
""",
),

Variable(
    abivarname="bound_temp@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['BoundingProcess_expert'],
    dimensions="scalar",
    defaultval=500,
    mnemonics="BOUND TEMPerature for the molecular dynamics (in Kelvin)",
    text="""
When the process will try a given model, this input variable is used to set the temperature for the molecular dynamics
""",
),

Variable(
    abivarname="bound_step@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['BoundingProcess_expert'],
    dimensions="scalar",
    defaultval=1000,
    mnemonics="BOUND number of STEP for the molecular dynamics",
    text="""
When the process will try a given model, this input variable is used to set the maximum number of molecular dynamics steps
""",
),

Variable(
    abivarname="dynamics@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Dynamics option for Multibinit",
    text="""
Set the Dynamics option for Multibinit. This option is equivalent to [[abinit:ionmov]]:

* 0 --> do nothing

* 12 --> Isokinetic ensemble molecular dynamics. The equation of motion of the ions in contact with a thermostat are solved with the algorithm proposed by Zhang [J. Chem. Phys. 106, 6102 (1997)], as worked out by Minary et al [J. Chem. Phys. 188, 2510 (2003)]. The conservation of the kinetic energy is obtained within machine precision, at each step.
**Purpose:** Molecular dynamics
**Cell optimization:** No (Use [[optcell]]=0 only)

* 13 --> Isothermal/isenthalpic ensemble. The equation of motion of the ions in contact with a thermostat and a barostat are solved with the algorithm proposed by Martyna, Tuckermann Tobias and Klein [Mol. Phys., 1996, p. 1117].
If optcell=1 or 2, the mass of the barostat ([[bmass]]) must be given in
addition.
**Purpose:** Molecular dynamics
**Cell optimization:** Yes (if [[optcell]]/=0)
**Related variables:** The time step ([[dtion]]), the temperatures
([[mdtemp]]), the number of thermostats ([[nnos]]), and the masses of
thermostats ([[qmass]]).
""",
),

Variable(
    abivarname="dtion@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=100,
    mnemonics="Delta Time for IONs",
    text="""
See [[abinit:dtion]]
""",
),

Variable(
    abivarname="ntime@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=200,
    mnemonics="Number of TIME step",
    text="""
Number of step for the dynamics
""",
),

Variable(
    abivarname="nnos@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Number of NOSe masses",
    text="""
See [[abinit:nnos]]
""",
),

Variable(
    abivarname="qmass@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['DynamicsMultibinit_basic'],
    dimensions=['[[abinit:nnos]]'],
    defaultval=0,
    mnemonics="Q thermostat MASS",
    text="""
See [[abinit:qmass]]
""",
),

Variable(
    abivarname="nctime@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="NetCdf TIME between output of molecular dynamics informations ",
    text="""
Set the number of step between output the molecular dynamics informations in the NetCDF file
""",
),

Variable(
    abivarname="temperature@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=325,
    mnemonics="molecular dynamics TEMPERATURE (in Kelvin)",
    text="""
Give the temperature of the dynamics in Kelvin
""",
),

Variable(
    abivarname="ncell@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['DynamicsMultibinit_basic'],
    dimensions=[3],
    defaultval=[6,6,6],
    mnemonics="Number of Cell",
    text="""
Give the size of the supercell for the dynamics
""",
),

Variable(
    abivarname="strtarget@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['DynamicsMultibinit_basic'],
    dimensions=[6],
    defaultval=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    mnemonics="STRess TARGET",
    text="""
See [[abinit:strtarget]]
""",
),

Variable(
    abivarname="bmass@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=10,
    mnemonics="Barostat MASS",
    text="""
See [[abinit:bmass]]
""",
),

Variable(
    abivarname="optcell@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="OPTimize the CELL shape and dimensions",
    text="""
See [[abinit:optcell]]
""",
),

Variable(
    abivarname="restartxf@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="RESTART from (X,F) history",
    text="""
See [[abinit:restartxf]]
""",
),

]
