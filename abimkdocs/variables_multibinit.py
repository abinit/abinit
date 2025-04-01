# coding: utf-8
from __future__ import print_function, division, unicode_literals, absolute_import

executable = "multibinit"

try:
    from abimkdocs.variables import ValueWithUnit, MultipleValue, Range
except ImportError:
    # This is needed for importing this module within Abipy
    from abipy.abio.abivar_database.variables import ValueWithUnit, MultipleValue, Range

ValueWithConditions = dict
Variable = dict

variables = [
Variable(
    abivarname="dipdip@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['LatticeModel_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="DIPole-DIPole interaction",
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    mnemonics="Dipole-Dipole range",
    added_in_version="before_v9",
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
    added_in_version="before_v9",
    text=r"""
Set the energy of the reference structure (from the DFT calculation)
if the energy of the reference is not specified in the DDB,
(for example if the DDB file of the ground states is not merged),
or not specified in the XML file), (by default in Hartree).
""",
),

Variable(
    abivarname="lwf_dynamics@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=["LWFModel_basic"],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Lattice Wannier Function DYNAMICS",
    added_in_version="9.8",
    text=r"""
Kind of LWF dynamics to run. Currently there is only the option 3.

* 0: Do not run LWF dynamics.

* 3: Run NVT LWF dynamics with the Berendsen thermalstat.

""",
),


Variable(
    abivarname="lwf_init_hist_fname@multibinit",
    varset="multibinit",
    vartype="string",
    topics=['LWFModel_basic'],
    dimensions="scalar",
    defaultval="",
    mnemonics="LWF INITIAL state HISTory file name",
    added_in_version="9.8",
    text=r"""
Specify the initial state of the multibinit LWF dynamics calculation, which can be a lwf_hist netcdf file, usually generated from previous LWF dynamics calculations. It is used when [[multibinit:lwf_init_state]]=4 The string must be enclosed between quotation marks, for example:

    lwf_init_hist_fname "last_step_lwf_hist.nc"
""",
),


Variable(
    abivarname="lwf_init_state@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=["LWFModel_basic"],
    dimensions="scalar",
    defaultval=1,
    mnemonics="Lattice Wannier Function INITial STATE",
    added_in_version="before_v9",
    text=r"""
Flag to initialize spin state.

* 1 --> The LWF amplitudes are homogenous random numbers between -0.1 to 0.1 Bohr.

* 2 --> The LWF amplitudes are 0.

* 4 --> Restart from an input spin hist file, as specified in [[multibinit:lwf_init_hist_fname]].
""",
),

Variable(
    abivarname="lwf_constraint@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=["LWFModel_expert"],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Lattice Wannier Function use CONSTRAINT",
    added_in_version="9.8",
    text=r"""
    Whether to use constraint in Lattice Wannier function dynamics. The constraints are defined in a LWF initial state file by three parameters:
    n_fixed_lwf: number of fixed LWF.
    fixed_lwf_ids: indices of fixed LWFs.
    fixed_lwf_values: values of fixed LWFs.
""",
),

Variable(
    abivarname="lwf_dt@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['LWFModel_basic'],
    dimensions="scalar",
    defaultval=100,
    mnemonics="Lattice Wannier Function Delta Time",
    added_in_version="9.8",
    text=r"""
Time step for lwf dynamics. Default value is 100.
Default unit is atomic unit (2.419e-17 s).
S, Sec or Second can be appended as unit. (e.g. 1e-16 Sec).
""",
),



Variable(
    abivarname="lwf_nctime@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=["LWFModel_basic"],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Lattice Wannier function dynamics NetCdf write per number of TIME steps",
    added_in_version="before_v9",
    text=r"""
Write LWF amplitude into netcdf file in every lwf_nctime of spin dynamics time steps.
""",
),



Variable(
    abivarname="lwf_ntime@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=["LWFModel_basic"],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Lattice Wannier function dynamics total Number of TIME steps",
    added_in_version="9.8",
    text=r"""
Total number of lattice Wannier function dynamics time steps.
""",
),

Variable(
    abivarname="lwf_pot_fname@multibinit",
    varset="multibinit",
    vartype="string",
    topics=["LWFModel_basic"],
    dimensions="scalar",
    defaultval="",
    mnemonics="Lattice Wannier function POTential File NAME",
    added_in_version="9.8",
    text=r"""
Specify the LWF potential file name in the multibinit lwf dynamics calculation, which is a netcdf file. The string must be enclosed between quotation marks:

    lwf_pot_fname "BaTiO3_lwf_pot.nc"
"""

),


Variable(
    abivarname="lwf_taut@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['LWFModel_basic'],
    dimensions="scalar",
    defaultval=1000,
    mnemonics="Lattice Wannier function dynamics relaxation time TAUT",
    added_in_version="9.8",
    text=r"""
    Parameter used in Berendsen lattice dynamics [[multibinit:lwf_dynamics]] = 3, in which the temperature is relaxed exponentially to the target temperature, with the characteristic time of lwf_taut.
    The default unit is atomic unit. But it is possible to use the second as the unit by adding Second or S at the end of the line, for example:

    lwf_taut 1d-15 S
""",
),



Variable(
    abivarname="lwf_temperature_start@multibinit",
    varset="multibinit",
    vartype="real",
    topics=["LWFModel_basic"],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="Lattice Wannier function TEMPERATURE START",
    added_in_version="9.8",
    text=r"""
Start point of variable temperature LWF dynamcis calculation (see [[multibinit:lwf_var_temperature]]) in lwf dynamics calculation.
""",
),

Variable(
    abivarname="lwf_temperature_end@multibinit",
    varset="multibinit",
    vartype="real",
    topics=["LWFModel_basic"],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="Lattice Wannier function TEMPERATURE END",
    added_in_version="9.8",
    text=r"""
End point of variable temperature LWF dynamics calculation (see [[multibinit:LWF_var_temperature]]) in LWF dynamics calculation.
""",
),

Variable(
    abivarname="lwf_temperature_nstep@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=["LWFModel_basic"],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Lattice Wannier function TEMPERATURE Number of STEPs",
    added_in_version="9.8",
    text=r"""
Number of steps in the variable temperature LWF dynamics calculation (see [[multibinit:lwf_var_temperature]]) in lwf dynamics calculation.
""",
),

Variable(
    abivarname="lwf_var_temperature@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=["LWFModel_basic"],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Lattice Wannier function VARiable TEMPERATURE",
    added_in_version="9.8",
    text=r"""
Switch for variable temperature calculation in LWF dynamics. 0: off. 1: on.
If switched on, a series of LWF dynamics calculation with temperatures from
[[multibinit:lwf_temperature_start]] to [[multibinit:lwf_temperature_end]],
with number of steps [[multibinit:lwf_temperature_nstep]] will be done.
The corresponding _lwf_hist.nc  file has the corresponding temperature in the filename.
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
    text=r"""
The number of vector shifts of the simple Monkhorst and Pack grid, needed to
generate the coarse grid of q points (for the series of fine grids, the number
of shifts it is always taken to be 1). Usually, put it to 1. Use 2 if BCC
sampling (Warning: not BCC lattice, BCC *sampling*), and 4 for FCC sampling
(Warning: not FCC lattice, FCC *sampling*).
""",
),

Variable(
    abivarname="prt_GF_csv@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['LatticeModel_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Print the Goal-Function values in a CSV file",
    added_in_version="v9",
    text=r"""
* 0 --> do nothing (Default)
* 1 --> Print the Goal-Function Values (GF) for all coefficients on a given processor
        at a given fit iteration into a csv file. Each iteration each processor
        prints a csv file. The colums are the GF on Energy, Force+Stresses, Forces, Stresses.
""",
),

Variable(
    abivarname="prt_model@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['LatticeModel_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="PRinT the MODEL",
    added_in_version="before_v9",
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
    abivarname="test_prt_ph@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['LatticeModel_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Prt test-set evaluation into file ph_test.nc",
    added_in_version="before_v9",
    text=r"""
Flag to activate the printing of the evaluation of the effective potential on to a test set into  a seperate netcdf file called ph_test.nc.

Forces, Energies, Stresses and Atomic Positions are written in ph_test.nc.
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
    added_in_version="before_v9",
    text=r"""
* 0  --> do not active the fit process

* 1 -->  Activate the fit process. This option will first generate a set of coefficients if [[multibinit:fit_generateCoeff]] is set to one. This generation is mainly parametrized by [[multibinit:fit_rangePower]] and [[multibinit:fit_cutoff]]. You can also provided a list of coefficients with the  model_anharmonic.MXL (see [[help:multibinit]]). Then the fit process will select the coefficients one by one up to [[multibinit:fit_ncoeff]] (see this [[cite:Escorihuela-Sayalero2017|paper]] for the details of the procedure).

* -1 --> **only for developers**, print the files for the scripts
""",
),

Variable(
    abivarname="fit_EFS@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_basic'],
    dimensions=[3],
    defaultval=[0,1,1],
    mnemonics="FIT on Energy, Forces, and or, Stresses",
    added_in_version="v9",
    text=r"""
Specifies on which first-principles quantities the anharmonic coefficients will be fitted.
The first number flags the fitting on the energies, the second the fitting on the forces, and the third on the stressses.

Default value is 0 1 1, so anharmonic coefficients get fitted on Forces and Stresses but not on energies
""",
),

Variable(
    abivarname="fit_factors@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['FitProcess_basic'],
    dimensions=[3],
    defaultval=[1,1,1],
    mnemonics="FIT FACTORS for Goal Function of Energy, Forces, and Stresses",
    added_in_version="v9",
    text=r"""
Specifies three factors for Energy, Forces and Stresses in the calcluation of the Goal Function which is to be minimized during the
fit process allowing to change the relative weight of the three quantities.

Default value is 1 1 1, equally balancing energy, forces and stresses.
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
    added_in_version="before_v9",
    text=r"""
Give the number of anharmonic coefficients to add in the model during the fit process
""",
),

Variable(
    abivarname="fit_ncoeff_per_iatom@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="FIT Number of COEFFicients per Irreducible ATOM",
    added_in_version="before_v9",
    text=r"""
Give the number of anharmonic coefficients per symmetric irreducible atom to add during fit process.
[[multibinit:fit_ncoeff]]/(nirred_atoms*fit_ncoeff_per_iatom) gives the number of fitting loops performed during the fit process, where in each loop fit_ncoeff_per_iatom coefficients for each irreducible atom will be added to the anharmonic potential.
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
    added_in_version="before_v9",
    text=r"""
Flag to activate the generation of the anharmonic coefficient for the fit process

**Related variables:** The  power range of the coefficients ([[multibinit:fit_rangePower]]), the cut off of the interactions ([[multibinit:fit_cutoff]]), the flag to add anharmonic strain ([[multibinit:fit_anhaStrain]]), the flag to add phonon strain coupling ([[multibinit:fit_SPCoupling]])
""",
),

Variable(
    abivarname="fit_iatom@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="FIT anharmonic terms around ATOM I",
    added_in_version="before_v9",
    text=r"""Gives the index of the atom in the reference structure around which the anharmonic terms will be generated.
If 0 (default) a loop over all atoms in the reference structure will be perforemed and fit_ncoeff coefficienst will be fitted and selected per atom.
If -1 all possible cross terms will be generated (e.G. (A_x-B_x)^2*(C_y-D_y)^1. This options generates much more terms.
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
    text=r"""
Flag to activate the strain  phonon coupling. This option will add coefficients like  (Sr-Ti)^1 (eta^4)
""",
),

Variable(
    abivarname="fit_dispterms@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="FIT anharmonic Strain-Phonon COUPLING coefficients",
    added_in_version="before_v9",
    text=r"""
Flag to activate the generation of pure displacement coefficients. This option will generate coefficients like (Sr-Ti)^2*(Sr-O), where only atomic displacements occur.

Default value: 1 -> displacement terms are generated.
""",
),

Variable(
    abivarname="fit_SPC_maxS@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="FIT Strain Phonon Coupling maximum Strain",
    added_in_version="v9",
    text=r"""
Set maximum power of strain body in strain-phonon coupling terms.
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
    text=r"""
Indices of the imposed coefficients during the fit process for the model:
""",
),


Variable(
    abivarname="fit_nimposecoeff@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="FIT Number of IMPOSEd COEFFicients",
    added_in_version="before_v9",
    text=r"""
Number of coefficients imposed with fixed value as in the input xml during the fit process for the model:

* -1 -->  fix all the coefficients

* 0  -->  do not fix coefficients

* n  -->  fix n coefficients (requires [[multibinit:fit_imposecoeff]] input variable)
""",
),

Variable(
    abivarname="fit_imposecoeff@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_expert'],
    dimensions=['[[multibinit:fit_nimposecoeff]]'],
    defaultval=0,
    mnemonics="FIT Number of IMPOSEd COEFFicients",
    added_in_version="before_v9",
    text=r"""
Indices of the imposed coefficients with fixed coefficient value during the fit process for the model.
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
    text=r"""
Indices of the banned coefficients during the fit process of the model
""",
),

Variable(
    abivarname="ts_option@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_useful'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="fit Training Set OPTION",
    added_in_version="before_v9",
    text=r"""
* 0 --> the Training is hist from ABINIT

* 1 --> the Training contains -1 * stress  (usualy output from VASP)
""",
),

Variable(
    abivarname="bound_factors@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['FitProcess_basic'],
    dimensions=[3],
    defaultval=[1,1,1],
    mnemonics="FACTORS for Goal Function of Energy, Forces, and Stresses during bounding process",
    added_in_version="v9",
    text=r"""
Specify three factors for Energy, Forces and Stresses in the calculation of the Goal Function which is to be minimized during the
bounding process, allowing to change the relative weights of the three quantities.

Default value is 1 1 1, equally balancing energy, forces and stresses.
""",
),

Variable(
    abivarname="bound_model@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['BoundingProcess_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="BOUND MODEL",
    added_in_version="before_v9",
    text=r"""
Flag to activate the bound process:

* 0 --> Do not activate the bound process

* 1 --> This option will generate all the possible combinations of coefficients from 1 to [[multibinit:bound_maxCoeff]]. Some constrains are imposed during the generation and the fit of the coefficients, they have to be positive and with even power. Finaly, the code will try all the possible combinations and try to find a bounded model.

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
    added_in_version="before_v9",
    text=r"""
Number of maximum additional coefficients for the bound process
""",
),

Variable(
    abivarname="bound_penalty@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['FitProcess_basic'],
    dimensions="scalar",
    defaultval=1.001,
    mnemonics="Goal Function penalty for determination of bounding coefficients",
    added_in_version="v9",
    text=r"""
Relative penalty for the determination of bounding coefficient values. The penalty defines the ration of the goal function before and after adding the coefficient. If the optimum value of the coefficient (-one that decreases the value of the goal function-) is negative a positive value that .
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
    text=r"""
Set the Dynamics option for Multibinit. This option is equivalent to [[abinit:ionmov]] for numbers < 100. For numbers >100, it uses algorithms implemented inside Multibinit:

* 0 --> do nothing

  * 1 --> Move atoms using molecular dynamics with optional viscous damping (friction linearly proportional to velocity). The viscous damping is controlled by the parameter "[[vis]]". If actual undamped molecular dynamics is desired, set [[vis]] to 0. The implemented algorithm is the generalisation of the Numerov technique (6th order), but is NOT invariant upon time-reversal, so that the energy is not conserved. The value [[ionmov]] = 6 will usually be preferred, although the algorithm that is implemented is lower-order. The time step is governed by [[dtion]].
**Purpose:** Molecular dynamics (if [[vis]] = 0), Structural optimization (if
[[vis]] >0)
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:** Viscous parameter [[vis]], time step [[dtion]], index
of atoms fixed [[iatfix]]

  * 2 --> Conduct structural optimization using the Broyden-Fletcher-Goldfarb-Shanno minimization (BFGS). This is much more efficient for structural optimization than viscous damping, when there are less than about 10 degrees of freedom to optimize. Another version of the BFGS is available with [[ionmov]]==22, and is apparently more robust and efficient than [[ionmov]]==2.
**Purpose:** Structural optimization
**Cell optimization:** Yes (if [[optcell]]/=0)
**Related variables:**

  * 6 --> Molecular dynamics using the Verlet algorithm, see [[cite:Allen1987a]] p 81]. The only related parameter is the time step ([[dtion]]).
**Purpose:** Molecular dynamics
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:** time step [[dtion]], index of atoms fixed [[iatfix]]

  * 7 --> Quenched Molecular dynamics using the Verlet algorithm, and stopping each atom for which the scalar product of velocity and force is negative. The only related parameter is the time step ([[dtion]]). The goal is not to produce a realistic dynamics, but to go as fast as possible to the minimum. For this purpose, it is advised to set all the masses to the same value (for example, use the Carbon mass, i.e. set [[amu]] to 12 for all type of atoms).
**Purpose:** Structural optimization
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:** time step [[dtion]], index of atoms fixed [[iatfix]]

  * 9 --> Langevin molecular dynamics.
**Purpose:** Molecular dynamics
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:** time step ([[dtion]]), temperatures ([[mdtemp]]) and
friction coefficient ([[friction]]).

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

  * 22 --> Conduct structural optimization using the Limited-memory Broyden-Fletcher-Goldfarb-Shanno minimization (L-BFGS) [[cite:Nocedal1980]]. The working routines were based on the original implementation of J. Nocedal available on netlib.org. This algorithm can be much better than the native implementation of BFGS in ABINIT ([[ionmov]] = 2) when one approaches convergence, perhaps because of better treatment of numerical details.
**Purpose:** Structural optimization
**Cell optimization:** Yes (if [[optcell]]/=0)
**Related variables:**

  * 24 --> Simple constant energy molecular dynamics using the velocity Verlet symplectic algorithm (second order), see [[cite:Hairer2003]]. The only related parameter is the time step ([[dtion]]).
**Purpose:** Molecular dynamics
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:** time step [[dtion]]

  * 25 --> Hybrid Monte Carlo sampling of the ionic positions at fixed temperature and unit cell geometry (NVT ensemble). The underlying molecular dynamics corresponds to [[ionmov]]=24. The related parameters are the time step ([[dtion]]) and thermostat temperature ([[mdtemp]]).
Within the HMC algorithm [[cite:Duane1987]], the trial states are generated via short $NVE$ trajectories (ten [[ionmov]]=24 steps in current implementation).
 The initial momenta for each trial are randomly sampled from Boltzmann distribution, and the final trajectory state is either accepted or rejected based on the Metropolis criterion.
 Such strategy allows to simultaneously update all reduced coordinates, achieve higher acceptance ratio than classical Metropolis Monte Carlo and better sampling efficiency for shallow energy landscapes [[cite:Prokhorenko2018]].
**Purpose:** Monte Carlo sampling
**Cell optimization:** No (Use [[optcell]] = 0 only)
**Related variables:** time step [[dtion]], thermostat temperature [[mdtemp]],

* 101 --> NVE ensemble with velocity Verlet algorithm  [[cite:Swope1982]] .
**Purpose:** Molecular dynamics
**Cell optimization:** No (Use [[optcell]]=0 only)
**Related variables:** The time step ([[dtion]]), the temperatures
([[multibinit:temperature]]). The time step should be small enough to make the energy conserved. The temperature is set to intialize the velocities of the atoms, which is in principle not preserved during the NVE run.

* 102 --> NVT ensemble with Langevin algorithm. [[cite:Vanden2006]] .
**Purpose:** Molecular dynamics
**Cell optimization:** No (Use [[optcell]]=0 only)
**Related variables:** The time step ([[dtion]]), the temperatures
([[multibinit:temperature]]), the friction [[multibinit:latt_friction]].
The atoms are coupled to the heat bath, which is represented by a gauss noise  in the forces, whose amplitude is defined by the temperature, and a friction term.


* 103 --> NVT ensemble. The temperature is approached by scaling the velocity of atoms. The method is proposed by Berendsen et al. in  J. Chem. Phys., 81 3684–3690 (1984) [[cite:Berendsen1984]]. Note that this method does NOT generate properly the thermostated ensemble. It does not have the correct distribution of the kinetic energy but have the correct average.  However, it approches the target temperature exponentially without oscillation, for which the steps can be easily controlled.
**Purpose:** Molecular dynamics
**Cell optimization:** No (Use [[optcell]]=0 only)
**Related variables:** The time step ([[dtion]]), the temperatures
([[multibinit:temperature]]), the ion relaxation time [[multibinit:latt_taut]].

* 120 --> Dummy mover. Atoms does not move. For testing only.
""",

# Not yet fully implemented. Need to be properly documented and tested. Disactivated temporarily.
#* 104 --> NPT ensemble with method. Similar to option 103, except the pressure is also scaled.
#**Purpose:** Molecular dynamics
#**Cell optimization:** No (Use [[optcell]]=0 only)
#**Related variables:** The time step ([[dtion]]), the temperatures
#([[multibinit:temperature]]), the ion relaxation time [[multibinit:latt_taut]], the pressure relaxation time [[multibinit:latt_taup]].




),

Variable(
    abivarname="dyn_chksym@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="DYNamics CHeK SYMmetry",
    added_in_version="v9",
    text=r"""
Flag to activate symmetry finder and imposition of symmetry of the restart structure before dynamics run, when restartxf is negativ.
Useful to do symmetry constrained relaxation with structural realxations algorithms.
Be cautious to use it with large number of atoms, symmetry detection might take a long time.

**Related variables:** Restart flag for multibinit dynamcis ([[multibinit:restartxf]]), symmetry on symmetry finder ([[multibinit:dyn_tolsym]]))
""",
),

Variable(
    abivarname="dyn_tolsym@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=1e-10,
    mnemonics="DYNamics TOLerance on SYMmetries",
    added_in_version="v9",
    text=r"""
Tolerance on symmetry finder.
**Related variables:** Activation flag for symmetry finder ([[multibinit:dyn_chksym]])
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
    text=r"""
    Parameter of the friction coefficient used in Langevin dynamics [[multibinit:dynamics]] =102. Typical value is 1e-4 to 1e-2.
""",
),

Variable(
    abivarname="latt_anharm_pot_fname@multibinit",
    varset="multibinit",
    vartype="string",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval="",
    mnemonics="LATTice HARMornic POTential File NAME",
    added_in_version="9.8",
    text=r"""
Specify the input coefficients from fitted polynomial in multibinit lattice calculation, which can be a xml file. The string must be enclosed between quotation marks:

    lat_anharm_pot_fname "BaTiO3_coeff.xml"

"""
),


Variable(
    abivarname="latt_harm_pot_fname@multibinit",
    varset="multibinit",
    vartype="string",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval="",
    mnemonics="LATTice HARMonic POTential File NAME",
    added_in_version="9.8",
    text=r"""
Specify the input derivative database of reference structure in multibinit lattice calculation, which can be a DDB file or a xml file. The string must be enclosed between quotation marks:

    latt_harm_pot_fname "BaTiO3.ddb"

"""
),

Variable(
    abivarname="latt_training_set_fname@multibinit",
    varset="multibinit",
    vartype="string",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval="",
    mnemonics="LATTice potential TRAINING SET File NAME",
    added_in_version="9.3.3",
    text=r"""
Specify the training set file name for building multibinit lattice potential, which is usually a abinit molecular dynamics history  netcdf file. The string must be enclosed between quotation marks:

    latt_training_set_fname "BaTiO3_md_hist.nc"

"""
),

Variable(
    abivarname="latt_test_set_fname@multibinit",
    varset="multibinit",
    vartype="string",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval="",
    mnemonics="LATTice potential TEST SET File NAME",
    added_in_version="9.8",
    text=r"""
Specify the test set file name for building multibinit lattice potential, which is usually a abinit molecular dynamics history  netcdf file. The string must be enclosed between quotation marks:

    latt_test_set_fname "BaTiO3_md_hist.nc"

"""
),


Variable(
    abivarname="latt_pot_fname@multibinit",
    varset="multibinit",
    vartype="string",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval="",
    mnemonics="LATTice POTential FileNAME",
    added_in_version="9.8",
    text=r"""
Specify the lattice potential file name in the multibinit lattice dynamics calculation, which can be a netcdf file. This variable is only used in the harmonic-only lattice potential for testing only. The string must be enclosed between quotation marks:

    latt_pot_fname "BaTiO3.nc"
"""
),


Variable(
    abivarname="latt_taut@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=1000,
    mnemonics="LATTice dynamics relaxation time TAUT",
    added_in_version="before_v9",
    text=r"""
    Parameter used in Berendsen lattice dynamics [[multibinit:dynamics]] =102 and 103, in which the temperature is relaxed exponentially to the target temperature, with the characteristic time of latt_taut.
    The default unit is atomic unit. But it is possible to use the second as the unit by adding Second or S at the end of the line, for example:

    latt_taut 1d-15 S

""",
),

Variable(
    abivarname="latt_taup@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=1000,
    mnemonics="LATTice dynamics relaxation time TAUP",
    added_in_version="before_v9",
    text=r"""
    Parameter used in Berendsen lattice dynamics [[multibinit:dynamics]] =103, in which the pressure is relaxed exponentially to the target temperature, with the characteristic time of latt_taup.
    The default unit is atomic unit. But it is possible to use the second as the unit by adding Second or S at the end of the line, for example:

    latt_taup 1d-15 S


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
    added_in_version="before_v9",
    text=r"""
Number of step for the dynamics.
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
    text=r"""
Give the size of the supercell for the dynamics
""",
),

Variable(
    abivarname="ncellmat@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['DynamicsMultibinit_basic'],
    dimensions=[3, 3],
    defaultval=[[1,0,0],[0,1,0], [0,0,1]],
    mnemonics="Number of superCELL MATtrix",
    added_in_version="9.5.0",
    text=r"""
Give the size of the supercell for the dynamics in the format of a matrix.
Currently allowed in spin dynamics and spin and LWF dynamics.
It will override the [[multibinit:ncell]] if specified.
"""
),

Variable(
    abivarname="outdata_prefix",
    varset="multibinit",
    vartype="string",
    topics=['Control_useful'],
    dimensions="scalar",
    defaultval=None,
    mnemonics="OUTput DATA PREFIX",
    added_in_version="9.3.3",
    text=r"""
Prefix for output files. Replaces the analogous entry in the obsolete *files_file* .
This variable is used when Abinit is executed with the new syntax:

    multibinit run.abi > run.log 2> run.err &

If this option is not specified, a prefix is automatically constructed from the input file name
provided the filename ends with an extension, e.g. `.ext`. (`.abi` is recommended)
If the input filename does not have a file extension, a default is provided.
"""
),


Variable(
    abivarname="strfact@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=100.0,
    mnemonics="STRess FACTor",
    added_in_version="v9.1",
    text=r"""
See [[abinit:strfact]]
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
    text=r"""
See [[abinit:restartxf]]
""",
),

Variable(
    abivarname="sel_EFS@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_basic'],
    dimensions=[3],
    defaultval=[0,1,1],
    added_in_version="v9",
    mnemonics="Select on Energy, Forces, and or, Stresses",
    text=r"""
Specifies on which goal function quantities the anharmonic coefficients will be selected.
The first number flags the selecting on the energies, the second the fitting on the forces, and the third on the stressses.

Default value is 0 1 1, so anharmonic coefficients get selected on Forces and Stresses but not on energies
""",
),

# The below are not yet functioning, comment out temporarily.
#Variable(
#    abivarname="spin_calc_correlation_obs@multibinit",
#    varset="multibinit",
#    vartype="integer",
#    topics=['SpinDynamicsMultibinit_basic'],
#    dimensions="scalar",
#    defaultval=0,
#    mnemonics="SPIN CALCulate CORRELATION OBServables",
#    added_in_version="before_v9",
#    text=r"""
#Flag to calculate spin correlation function based observables.
#
#* 0 --> do not calculate.
#
#* 1 --> calculate.
#""",
#),
#
#
#Variable(
#    abivarname="spin_calc_traj_obs@multibinit",
#    varset="multibinit",
#    vartype="integer",
#    topics=['SpinDynamicsMultibinit_basic'],
#    dimensions="scalar",
#    defaultval=0,
#    mnemonics="SPIN CALCulate TRAJectory based OBServables",
#    added_in_version="before_v9",
#    text=r"""
#Flag to calculate spin trajectory based observables. (Nothing included yet.)
#
#* 0 --> do not calculate.
#
#* 1 --> calculate.
#""",
#),

Variable(
    abivarname="spin_calc_thermo_obs@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="SPIN CALCulate THERMO dynamics OBServables",
    added_in_version="before_v9",
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
    added_in_version="before_v9",
    text=r"""
Gilbert damping factor in LLG equation for spin dynamics.

* negative value --> use damping factor from spin xml file.

* positive value --> use as damping factor. The value should be between 0.0 and 1.0 (both included).
""",
),

#Variable(
#    abivarname="spin_dipdip@multibinit",
#    varset="multibinit",
#    vartype="integer",
#    topics=['SpinDynamicsMultibinit_basic'],
#    dimensions="scalar",
#    defaultval=0,
#    mnemonics="SPIN DIPole DIPole interaction",
#    added_in_version="before_v9",
#    text=r"""
#* 0 --> Switch off spin dipole-dipole interaction.
#
#* 1 --> Switch on spin dipole-dipole interaction.
#    (Not yet implemented.)
#""",
#),

Variable(
    abivarname="spin_dt@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=100,
    mnemonics="SPIN Delta Time",
    added_in_version="before_v9",
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
    added_in_version="before_v9",
    text=r"""
Flag to run spin dynamics.

* 0 --> Do not run spin dynamics.

* 1 --> Run spin dynamics with HeunP integration method.

* 2 --> Run spin dynamics with Depondt-Mertens method [[cite:Depondt2009]].

* 3 --> Run Monte Carlo.

* 20 --> Dummy mover. Spin will not rotate. For test only.

The HeunP method does less computation for each step,
whereas the Depondt-Mertens method allow larger time step.
For system with very simple interaction terms, HeunP could be faster.
Otherwise, use Depondt-Mertens method can be more efficient.
""",
),

Variable(
    abivarname="spin_init_hist_fname@multibinit",
    varset="multibinit",
    vartype="string",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval="",
    mnemonics="SPIN INITIAL state HISTory file name",
    added_in_version="9.3.3",
    text=r"""
Specify the initial state of the multibinit spin dynamics calculation, which can be a spin_hist netcdf file, usually generated from previous spin dynamics calculations. It is used when [[multibinit:spin_init_state]]=4 The string must be enclosed between quotation marks, for example:

    spin_init_hist_fname "last_step_spin_hist.nc"
"""

),



Variable(
    abivarname="spin_init_state@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="SPIN INITial STATE",
    added_in_version="before_v9",
    text=r"""
Flag to initialize spin state.

* 1 --> Random spin state using uniform random numbers.

* 2 --> Reference spin state from potential file if present.

* 3 --> State with q-vector using [[multibinit:spin_init_qpoint]], [[multibinit:spin_init_rotate_axis]], and [[multibinit:spin_init_orientation]]. Please check default values for those variables.

* 4 --> Restart from last step of input spin hist file, as specified in [[multibinit:spin_init_hist_fname]].
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
    text=r"""
Total number of spin dynamics time  steps for preparing the system.
The results of these time step are not written to trajectory spinhist.nc file,
And they are not used for calculating the observables.
""",
),


Variable(
    abivarname="spin_init_orientation@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions=[3],
    defaultval=[0,0,1],
    mnemonics="SPIN INITial ORIENTATION",
    added_in_version="before_v9",
    text=r"""
Spin initial orientation. It is used for setting the initial spin in a supercell.
    For a spin in a cell labeled with R, the rotation angle is $2\pi Q\cdot R$
    from the initial orientation along the rotate axis.
    Default is along z(0,0,1) direction.
""",
),




Variable(
    abivarname="spin_init_qpoint@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions=[3],
    defaultval=[0,0,0],
    mnemonics="SPIN INITial QPOINT",
    added_in_version="before_v9",
    text=r"""
Spin wave vector. It is used for setting the initial spin in a supercell.
    For a spin in a cell labeled with R, the rotation angle is $2\pi Q\cdot R$
    from the initial orientation along the rotate axis.
    Default is Gamma (0, 0, 0).
""",
),



Variable(
    abivarname="spin_init_rotate_axis@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions=[3],
    defaultval=[1,0,0],
    mnemonics="SPIN INITial ROTATE AXIS",
    added_in_version="before_v9",
    text=r"""
Spin initial rotate axis. It is used for setting the initial spin in a supercell.
    For a spin in a cell labeled with R, the rotation angle is $2\pi Q\cdot R$
    from the initial orientation along the rotate axis.
    Default is along x axis (1,0,0).
""",
),

Variable(
    abivarname="spin_pot_fname@multibinit",
    varset="multibinit",
    vartype="string",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval="",
    mnemonics="SPIN POTential File NAME",
    added_in_version="9.3.3",
    text=r"""
Specify the spin potential file name in the multibinit spin dynamics calculation, which can be either a xml or a netcdf file. The string must be enclosed between quotation marks:

    spin_pot_fname "Fe.xml"
"""

),

Variable(
    abivarname="spin_projection_qpoint@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions=[3],
    defaultval=[0,0,0],
    mnemonics="SPIN PROJECTION QPOINT",
    added_in_version="before_v9",
    text=r"""
Spin wave vector. It is used for getting the total spin. $M_{tot}=\sum_i M_i exp(i q \cdot R_i)$. The unit is the reciprocal lattice vectors of the unitcell.
    Default is Gamma. (0, 0, 0)
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
    added_in_version="before_v9",
    text=r"""
Add single ion anisotropy term to the spin model hamiltonian.
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
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
    added_in_version="before_v9",
    text=r"""
Number of steps in the variable temperature spin dynamics calculation (see [[multibinit:spin_var_temperature]]) in spin dynamics calculation.
""",
),

Variable(
    abivarname="test_effpot@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['LatticeModel_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="TEST EFFective POTential",
    added_in_version="before_v9",
    text=r"""
* 0 --> nothing.
* 1 --> Evaluate the effective potential with respect to given test-set and calculate mean square differences between ab-initio energy/forces and model energy/forces
""",
),

Variable(
    abivarname="tolmxf@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['DynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=2e-5,
    mnemonics="TOLerance on the MaXimal Force",
    added_in_version="v9",
    text=r"""
Sets a maximal absolute force tolerance (in hartree/Bohr) below which BFGS structural relaxation iterations will stop. Corresponds to [[tolmxf]] of Abinit.
""",
),


Variable(
abivarname="analyze_anh_pot@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['LatticeModel_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="ANALYZE ANHarmonic POTential",
    added_in_version="before_v9",
    text=r"""
* 0 --> Nothing.
* 1 --> Print energy contribution of each anharmonic term in the effective Potential.
        If it is a Molecular Dynamics (MD) run, the contribution of each term is printed for each MD-step into MD_anharmonic_terms_energy.dat .
        If the effective potential is tested against a test set the contribution of each term for each configuration in the test set is printed in TES_anharmonic_terms_energy.dat .
        If the effective potential is fitted, the contribution of each selected term for each configuration in the training set is printed in TRS_anharmonic_terms_energy.dat
""",
),


Variable(
    abivarname="opt_effpot@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="OPTimize EFFective POTential",
    added_in_version="before_v9",
    text=r"""
* 0 --> Nothing.
* 1 --> Turn on reading of optimization of effective potential keywords (opt_).
        The optimization process gives the user the ability to refit the coefficients of specified terms with respect to the training set while keeping the rest fixed.

**Related variables:** The number of coefficients to refit ([[multibinit:opt_ncoeff]]), the  indices of the coefficients to optimize ([[multibinit:opt_coeff]]).
""",
),

Variable(
    abivarname="opt_factors@multibinit",
    varset="multibinit",
    vartype="real",
    topics=['FitProcess_basic'],
    dimensions=[3],
    defaultval=[1,1,1],
    mnemonics="FACTORS for Goal Function of Energy, Forces, and Stresses during optimization of coefficients",
    added_in_version="v9",
    text=r"""
Specifies three factors for Energy, Forces and Stresses in the calcluation of the Goal Function which is to be minimized during the
optimization process allowing to change the relative weight of the three quantities.

Default value is 1 1 1, equally balancing energy, forces and stresses.
""",
),

Variable(
  abivarname="opt_ncoeff@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_expert'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="OPTimize NUMBER of COEFFicients",
    added_in_version="before_v9",
    text=r"""
* Number of anharmonic terms to refit in the effective potential.
**Related variables:** [[multibinit:opt_coeff]]
""",
),

Variable(
    abivarname="opt_coeff@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['FitProcess_expert'],
    dimensions=['[[multibinit:opt_ncoeff]]'],
    defaultval=0,
    mnemonics="OPTimize Cofficients",
    added_in_version="before_v9",
    text=r"""
Indices of the terms to refit in the effective potential.
""",
),

Variable(
    abivarname="randomseed@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['SpinDynamicsMultibinit_expert'],
    dimensions='scalar',
    defaultval=0,
    mnemonics="RANDOM SEED",
    added_in_version="9.8",
    text=r"""
Random seed to be used in Multibinit spin/LWF dynamics.
It should be 0, or a large positive integer.
The default value 0 means it will use the current clock time.
DO NOT set this number unless you want to repeat the previous result. If a series
of dynamics is done with the same seed, the results could be wrong due to the
artificial periodicity of the random number that is generated. Even [[randomseed@multibinit]] is set, it is not guranteed
that the previous result can be recovered, as the generation of numbers is also affected by the number of
processors, type of type of CPU, compiler,  and version of MULTIBINIT.
""",
),



Variable(
    abivarname="slc_coupling@multibinit",
    varset="multibinit",
    vartype="integer",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="SpinLatticeCoupling_Coupling",
    added_in_version="9.0.0",
    text=r"""
Which spin-lattice coupling terms are used in the calculation, different terms can be combined in a binary fashion, i.e. 1010 turns on all terms quadratic in spin.

* 0    --> No coupling.

* 0001 --> Coupling term linear in spin and lattice coordinate

* 0010 --> Coupling term quadratic in spin and linear in lattice coordinate

* 0100 --> Coupling term linear in spin and quadratic in lattice coordinate

* 1000 --> Coupling term quadratic in spin and lattice coordinate
""",
),

Variable(
    abivarname="slc_pot_fname@multibinit",
    varset="multibinit",
    vartype="string",
    topics=['SpinDynamicsMultibinit_basic'],
    dimensions="scalar",
    defaultval="",
    mnemonics="SLC POTential FileNAME",
    added_in_version="9.3.3",
    text=r"""
Specify the spin-lattice-coupling potential file name in the coupled spin-lattice multibinit dynamics calculation, which can be a netcdf file. The string must be enclosed between quotation marks:

    slc_pot_fname "BaTiO3_slc.nc"
"""

),

Variable(
    abivarname="outdata_prefix@multibinit",
    varset="multibinit",
    vartype="string",
    topics=['Control_useful'],
    dimensions="scalar",
    defaultval=None,
    mnemonics="OUTput DATA PREFIX",
    added_in_version="9.8.0",
    text=r"""
Prefix for output files. Replaces the analogous entry in the obsolete *files_file* .
This variable is used when MULTIBINIT is executed with the new syntax:

    multibinit run.abi > run.log 2> run.err &

If this option is not specified, a prefix is automatically constructed from the input file name
provided the filename ends with an extension, e.g. `.ext`. (`.abi` is recommended)
If the input filename does not have a file extension, a default is provided.

Example:

    outdata_prefix = "t01_o"

See also [[outdata_prefix@abinit]] and [[outdata_prefix@anaddb]]

"""
),





]
