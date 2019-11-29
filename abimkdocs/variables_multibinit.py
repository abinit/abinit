# coding: utf-8
from __future__ import print_function, division, unicode_literals, absolute_import

executable = "multibinit"

from abimkdocs.variables import ValueWithUnit, MultipleValue, Range
#from abipy.abio.abivar_database.variables import ValueWithUnit, MultipleValue, Range, ValueWithConditions
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
Flag to activate the bound process:

* 0 --> Do not activate the bound process

* 1 --> This option will generate all the possible combinaisons of coefficients from 1 to [[multibinit:bound_maxCoeff]]. Some constrains are imposed during the generation and the fit of the coefficients, they have to be positive and with even power. Finaly, the code will try all the possible combinaisons and try to find a bounded model.

* 2 -->  **new version** This option will generate a set of coefficients with a power range defined by [[multibinit:bound_rangePower]] and keep only the coefficients with even power. Then the procedure is similar to the fit process with the constrains to only keep positive coefficients. The bound process will select the coefficients one by one up to [[multibinit:bound_maxCoeff]] and try if the model is bound at each step of the process.

**Related variables:1 and 2** The number of maximum additional coefficient in the polynome ([[multibinit:bound_maxCoeff]]), the  power range for the additional coefficients ([[multibinit:bound_rangePower]]), the cut off of the additional interactions ([[multibinit:bound_cutoff]])

*3 --> Check each anharmonic term in the effective potential. If the term contains has a negative coefficient and is even in its displacement or contains odd powers in the displacement generate high order bounding terms of the same combination of displacement within the range of powers defined by the user ([[multibinit:bound_rangePower]]). The coefficients of the added high-order terms are optimized until the precision of the original effective potential is retained. 

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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
Set the Dynamics option for Multibinit. This option is equivalent to [[abinit:ionmov]] for numbers < 100. For numbers >100, it uses algorithms implemented inside Multibinit:

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

* 101 --> NVE ensemble with velocity Verlet algorithm  [[cite:Swope1982]] . 
**Purpose:** Molecular dynamics
**Cell optimization:** No (Use [[optcell]]=0 only)
**Related variables:** The time step ([[dtion]]), the temperatures
([[multibinit:temperature]]).


* 102 --> NVT ensemble with Langevin algorithm. [[cite:Vanden2006]] . 
**Purpose:** Molecular dynamics
**Cell optimization:** No (Use [[optcell]]=0 only)
**Related variables:** The time step ([[dtion]]), the temperatures
([[multibinit:temperature]]), the friction [[multibinit:latt_friction]].


* 103 --> NVT ensemble. The temperature is approached by scaling the velocity of atoms. The method is proposed by Berendsen et al. in  J. Chem. Phys., 81 3684–3690 (1984) [[cite:Berendsen1984]]. Note that this method does NOT generate properly the thermostated ensemble. It does not have the correct distribution of the kinetic energy but have the correct average.  However, it approches the target temperature exponentially without oscillation, for which the steps can be easily controlled.
**Purpose:** Molecular dynamics
**Cell optimization:** No (Use [[optcell]]=0 only)
**Related variables:** The time step ([[dtion]]), the temperatures
([[multibinit:temperature]]), the ion relaxation time [[multibinit:latt_taut]].

* 104 --> NPT ensemble with method. Similar to option 103, except the pressure is also scaled. 
**Purpose:** Molecular dynamics
**Cell optimization:** No (Use [[optcell]]=0 only)
**Related variables:** The time step ([[dtion]]), the temperatures
([[multibinit:temperature]]), the ion relaxation time [[multibinit:latt_taut]], the pressure relaxation time [[multibinit:latt_taup]].

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
    text=r"""
See [[abinit:dtion]]
""",
),

Variable(
    abivarname="latt_friction@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=1e-4,
    mnemonics="LATTice dynamics FRICTION parameter",
    text=r"""
    Parameter of the friction used in Langevin dynamcis [[multibinit:dynamics]] =101.
""",
),


Variable(
    abivarname="latt_taut@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=1000,
    mnemonics="LATTice dynamics relaxation time TAUT",
    text=r"""
    Parameter used in Berendsen lattice dynamcis [[multibinit:dynamics]] =102 and 103, in which the temperature is relaxed exponentially to the target temperature, with the characteristic time of latt_taut.
    The unit is atomic unit, same as [[dtion]].
""",
),

Variable(
    abivarname="latt_taup@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=1000,
    mnemonics="LATTice dynamics relaxation time TAUT",
    text=r"""
    Parameter used in Berendsen lattice dynamcis [[multibinit:dynamics]] =103, in which the pressure is relaxed exponentially to the target temperature, with the characteristic time of latt_taup.
    The unit is atomic unit, same as [[dtion]].
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
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
    text=r"""
See [[abinit:restartxf]]
""",
),

Variable(
    abivarname="spin_calc_correlation_obs@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="SPIN CALCulate CORRELATION OBServables",
    text=r"""
Flag to calculate spin correlation function based observables.

* 0 --> do not calculate.

* 1 --> calculate.
""",
),


Variable(
    abivarname="spin_calc_traj_obs@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="SPIN CALCulate TRAJectory based OBServables",
    text=r"""
Flag to calculate spin trajectory based observables. (Nothing included yet.)

* 0 --> do not calculate.

* 1 --> calculate.
""",
),


Variable(
    abivarname="spin_calc_thermo_obs@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="SPIN CALCulate THERMO dynamics OBServables",
    text=r"""
Flag to calculate spin thermo dynamics observables,
including the specific heat, magnetic susceptibility, Binder U4 value.
It's recommend to always calculate these observables.

* 0 --> do not calculate.

* 1 --> calculate.
""",
),


Variable(
    abivarname="spin_damping@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=-1.0,
    mnemonics="SPIN gilbert DAMPING factor",
    text=r"""
Gilbert damping factor in LLG equation for spin dynamics.

* negative value --> use damping factor from spin xml file.

* positive value --> use as damping factor. The value should be between 0.0 and 1.0 (both included).
""",
),

Variable(
    abivarname="spin_dipdip@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="SPIN DIPole DIPole interaction",
    text=r"""
* 0 --> Switch off spin dipole-dipole interaction.

* 1 --> Switch on spin dipole-dipole interaction.
    (Not yet implemented.)
""",
),

Variable(
    abivarname="spin_dt@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=100,
    mnemonics="SPIN Delta Time",
    text=r"""
Time step for spin dynamics. Default value is 100.
Default unit is atomic unit (2.419e-17 s).
S, Sec or Second can be appended as unit. (e.g. 1e-16 Sec).
""",
),


Variable(
    abivarname="spin_dynamics@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="SPIN DYNAMICS",
    text=r"""
Flag to run spin dynamics.

* 0 --> Do not run spin dynamics.

* 1 --> Run spin dynamics with HeunP integration method.

* 2 --> Run spin dynamics with Depondt-Mertens method [[cite:Depondt2009]].

* 3 --> Run Monte Carlo.

The HeunP method does less computation for each step,
whereas the Depondt-Mertens method allow larger time step.
For system with very simple interaction terms, HeunP could be faster.
Otherwise, use Depondt-Mertens method can be more efficient.
""",
),

Variable(
    abivarname="spin_init_state@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="SPIN INITial STATE",
    text=r"""
Flag to initialize spin state. (only option 1 and 2 are implemented.)

* 0 --> Read from spinhist netcdf file.

* 1 --> Random spin state using uniform random numbers.

* 2 --> Ferromagnetic state.

* 3 --> State with q-vector using [[multibinit:spin_qpoint]]

* 4 --> Random spin state with temperature of [[multibinit:spin_temperature]]
""",
),


Variable(
    abivarname="spin_mag_field@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions=[3],
    defaultval=[0,0,0],
    mnemonics="SPIN Magnetic Field",
    text=r"""
External magnetic field. Unit: Tesla.
""",
),


Variable(
    abivarname="spin_nctime@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="SPIN NetCdf write per number of TIME steps",
    text=r"""
Write spin into netcdf file in every spin_nctime of spin dynamics time steps.
""",
),

Variable(
    abivarname="spin_ntime@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="SPIN dynamics total Number of TIME steps",
    text=r"""
Total number of spin dynamics time  steps.
""",
),


Variable(
    abivarname="spin_ntime_pre@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="SPIN dynamics total Number of TIME steps for PREparing",
    text=r"""
Total number of spin dynamics time  steps for preparing the system.
The results of these time step are not written to trajectory spinhist.nc file,
And they are not used for calculating the observables.
""",
),

Variable(
    abivarname="spin_qpoint@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions=[3],
    defaultval=[0,0,0],
    mnemonics="SPIN QPOINT",
    text=r"""
Spin wave vector. It is used for getting the total spin. $M_{tot}=\sum_i M_i exp(i q \cdot R_i)$. The unit is the reciprocal lattice vectors of the unitcell.
""",
),


Variable(
    abivarname="spin_sia_add@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="SPIN Single Ion Anistropy ADD",
    text=r"""
Add single ion anistropy term to the spin model hamiltonian.
with user defined values (see [[multibinit:spin_sia_k1amp]] and [[multibinit:spin_sia_k1dir]].

* 0 --> Do not add, use the term defined in the spin model xml file.

* 1 --> Override the term in spin model xml file.

* 2 --> Add to the value defined in spin model xml file.

Default value: 0.
""",
),


Variable(
    abivarname="spin_sia_k1amp@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="SPIN Single Ion Anistropy K1 AMPtitude",
    text=r"""
User defined amplitude of single ion anistropy. Only used when [[multibinit:spin_sia_add]] is not 0.
The direction is defined with [[multibinit:spin_sia_k1dir]]. The unit is Ha. To use eV or Ry as unit, 
put eV or Ry at the end.
""",
),


Variable(
    abivarname="spin_sia_k1dir@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions=[3],
    defaultval=[0.0,0.0,1.0],
    mnemonics="SPIN Single Ion Anistropy K1 DIRection",
    text=r"""
User defined direction of single ion anistropy. Only used when [[multibinit:spin_sia_add]] is not 0.
It will be automatically normalized to 1.0.  The amplitude is defined with [[multibinit:spin_sia_k1amp]].
Default value: [0.0, 0.0,1.0].
""",
),


Variable(
    abivarname="spin_temperature@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=325,
    mnemonics="SPIN TEMPERATURE",
    text=r"""
Temperature of spin for spin dynamics. Unit: Kelvin.
Default value: 325.
""",
),


Variable(
    abivarname="spin_var_temperature@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="SPIN  VARiable TEMPERATURE",
    text=r"""
Switch for variable temperature calculation. 0: off. 1: on.
If switched on, a series of spin dynamics calculation with temperatures from
[[multibinit:spin_temperature_start]] to [[multibinit:spin_temperature_end]],
with number of steps [[multibinit:spin_temperature_nstep]] will be done.
The corresponding _spinhist.nc  file has the corresponding temperature in the filename.
""",
),


Variable(
    abivarname="spin_write_traj@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="SPIN WRITE TRAJectory to spinhist.nc file",
    text="""
Switch for writting of spin trajectory file. 0: off. 1 on.
The trajectory is needed for postprocessing of correlation functions.
""",
),


Variable(
    abivarname="spin_temperature_start@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="SPIN TEMPERATURE START",
    text=r"""
Start point of variable temperature spin dynamcis calculation (see [[multibinit:spin_var_temperature]]) in spin dynamics calculation.
""",
),

Variable(
    abivarname="spin_temperature_end@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="SPIN TEMPERATURE END",
    text=r"""
End point of variable temperature spin dynamics calculation (see [[multibinit:spin_var_temperature]]) in spin dynamics calculation.
""",
),

Variable(
    abivarname="spin_temperature_nstep@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="SPIN TEMPERATURE Number of STEPs",
    text=r"""
Number of steps in the variable temperature spin dynamics calculation (see [[multibinit:spin_var_temperature]]) in spin dynamics calculation.
""",
),

Variable(
    abivarname="test_effpot@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['LatticeModel_Basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="TEST EFFective POTential",
    text=r"""
* 0 --> nothing.
* 1 --> Evaluate the effective potential with respect to given test-set and calculate mean square differences between ab-initio energy/forces and model energy/forces""",
),

Variable(
    abivarname="analyze_anh_pot@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['LatticeModel_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="ANALYZE ANHarmonic POTential",
    text=r"""
* 0 --> nothing.
* 1 --> Print energy contribution of each anharmonic term in the effective Potential. 
        If it is a Molecular Dynamics (MD) run the contribution of each term is printed for each MD-step into MD_anharmonic_terms_energy.dat
        If the effective potential is tested against a test set the contribution of each term for each configuration in the test is set is printed in TES_anharmonic_terms_energy.dat 
        If the a effective potential is fitted the contribution of each selected term for each configuration in the training set is printed in TRS_anharmonic_terms_energy.dat""",
),


Variable(
    abivarname="opt_effpot@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="OPTimize EFFective POTential",
    text=r"""
* 0 --> nothing.
* 1 --> Turn on reading of optimization of effective potential keywords (opt_)
        The optimization process gives the user the ability to refit the coefficients of specified terms with respect to the training set while keeping the rest fixed.

**Related variables:** The number of coefficients to refit ([[multibinit:opt_ncoeff]]), the  indexes of the coefficients to optimize ([[multibinit:opt_coeff]])""" 
),

Variable(
  abivarname="opt_ncoeff@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="OPTimize NUMBER of COEFFicients",
    text=r"""
* Number of anharmonic terms to refit in the effective potential""" 
),

Variable(
    abivarname="opt_coeff@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_expert'],
    dimensions=['[[multibinit:opt_ncoeff]]'],
    defaultval=0,
    mnemonics="OPTimize Cofficients",
    text=r"""
Indexes of the terms to refit in the effective potential. """,
),


]


