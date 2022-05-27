<!--
---
authors: AM, FR, LB
---

# Second tutorial on MULTIBINIT

## Build a harmonic lattice model and run a dynamics

This lesson aims at showing how to build a harmonic model by using a second-principles approach
for lattice dynamics simulations based on atomic potentials fitted on first-principles calculations.

**Before beginning, it is very important to read the reference [[cite:Wojdel2013]].**

With this lesson, you will learn to:

  * Compute the model from the DDB
  * Generate the XML for the model
  * Run a dynamics within MULTIBINIT

In this tutorial, all the knowledge about the Density Functional Theory (DFT) and Density Functional Perturbation Theory (DFPT) should have been already acquired.
In addition, the DFPT is a key feature of ABINIT for MULTIBINIT thus, in order to learn how to use the DFPT (producing the related DDB) and the associated code to merge different DDB files,
please have a look at the tutorials on [[lesson:rf1| Phonon response]], [[lesson:elastic|strain response]] and [[lesson:rf2| mrgddb]] before to continue.
After these tutorials you should be able to compute a full DFPT calculation and a complete DDB file.
Thereby this tutorial will not provide the inputs for ABINIT, that you can find, however, in the previously cited tutorials.

The AGATE software, used to make the analysis of the results, is also required for this tutorial. You can install it on debian with:

    sudo add-apt-repository ppa:piti-diablotin/abiout
    sudo apt-get update && sudo apt-get install abiout

[TUTORIAL_README]

## 1 The Harmonic part of the lattice model

As described in [[cite:Wojdel2013]], the construction of a model starts by defining the reference structure (RS) of the system under study:

![Schema 1](lattice_model_assets/reference.png)
: Fig. 1: Example of cubic RS made by two different (black and red) atomic species.

The choice of this RS is fundamental since the model will be based on perturbations acting on it.
Once this choice is done, make sure than your system is fully relaxed and perform a single DFT calculation on this system in order to generate a DDB file.
In this tutorial, we will take as an example of a material without instabilities, the perovskite CaTiO3 in the Pnma phase.

From the RS, now we consider the perturbations: the set of atomic displacements $\boldsymbol{u}$ and lattice strain $\boldsymbol{\eta}$:

![Schema 2](lattice_model_assets/deformation.png)
: Fig. 2: Example of lattice perturbations: atomic displacements and strain.

At **the harmonic level**, we can express the potential energy as a sum of three contributions as a function of the set of perturbations ($\boldsymbol{u}$, $\boldsymbol{\eta}$):

$$\displaystyle  E^{harm}(\boldsymbol{u},\boldsymbol{\eta}) =  E^{harm}(\boldsymbol{u}) + E^{harm}(\boldsymbol{u},\boldsymbol{\eta}) + E^{harm}(\boldsymbol{\eta})$$

This calculation requires:

  * the computation of the phonon response (including short range and dipole-dipole interactions):

$$\displaystyle  E^{harm}(\boldsymbol{u}) \Longrightarrow \underbrace{\frac{\partial^2 E}{\partial
          \boldsymbol{u}^2}}_{\substack{\text{Inter-atomic}\\\text{force constants}}} ,
     \underbrace{\frac{\partial^2 E}{\partial
          {\boldsymbol{\cal{E}}}^2}}_{\text{Dielectric tensor}} ,
     \underbrace{\frac{\partial^2 E}{\partial{\boldsymbol{\cal{E}}} \partial \boldsymbol{u}}}_{\text{Effective charges}} $$

  * the computation of the strain response:
$$\displaystyle  E^{harm}(\boldsymbol{\eta}) \Longrightarrow \underbrace{\frac{\partial^2 E}{\partial
            \boldsymbol{\eta}^2}}_{\text{Elastic constants}} $$

  * the computation of the strain-phonon coupling:
$$\displaystyle E^{harm}(\boldsymbol{u},\boldsymbol{\eta}) \Longrightarrow \underbrace{\frac{\partial^2
            E}{\partial\boldsymbol{\eta}\partial\boldsymbol{u}}}_{\substack{\text{Internal strain}}} $$

We note that all the needed quantities related to the harmonic part can be computed by using the [[topic:DFPT|DFPT]] method and, once obtained, they are stored in a database (DDB) file.
In the case of an *insulator*, the dipole-dipole interaction will be recomputed directly by MULTIBINIT.
Thereby you need to provide into the DDB file the clamped-ion dielectric tensor and the Born effective charges.

!!! note

    Do not forget to include in the final DDB file the DDB obtained from the RS single DFT calculation (you can still use [[help:mrgddb | mrgddb]]). Indeed, the DDB file is also an output of a DFT calculation and, in order to build of a model, it is important to include it in the DDB file to be used.
    
-->
---
authors: LB, PhG
---

# First tutorial on MULTIBINIT

## Build a second-principles effective atomistic model and run finite-temperature lattice dynamics simulations

This lesson aims at learning how to build an effective atomistic model from a set of first-principles data and then use it for simulations at finite temperature.

**Before beginning, it is very important to read the reference [[cite:Wojdel2013]].**

Within this lesson, we will describe :

  * the complete set of first-principles data to be provided.
  * the steps for constructing a model for a prototypical compound (BaHfO$_3$).
  * the way to perform a finite temperature simulation from the previous model.

In this tutorial, we make the hypothesis that you have been already acquired a practical knowledge regarding Density Functional Theory (DFT) and Density Functional Perturbation Theory (DFPT) .
In particular, DFPT is a key feature of ABINIT directly exploited by MULTIBINIT. In order to learn how to use the DFPT (producing the related DDB) and the associated code to merge different DDB files, please have a look at the tutorials on [[lesson:rf1| Phonon response]], [[lesson:elastic|strain response]] and [[lesson:rf2| Mrgddb]]. After these tutorials, you should be able to compute a full DFPT calculation and a complete DDB file.
In this tutorial will not provide the inputs for ABINIT DFPT calculations (that you can be found in the previously cited tutorials) but instead the final DDB resulting from them.


The [AGATE](https://github.com/piti-diablotin/agate) software is also required for this tutorial, as a tool for the analysis of the results. You can install it on debian with:

    sudo add-apt-repository ppa:piti-diablotin/abiout
    sudo apt-get update && sudo apt-get install abiout

[TUTORIAL_README]

## 1 Method and first-principles inputs

As described in [[cite:Wojdel2013]], the construction of a lattice model with MULTIBINIT consists in determining an explicit form of the Born-Oppenheimer (BO) energy surface around a reference structure (RS), in terms of individual atomic displacements $\boldsymbol{u}$ and macroscopic strains $\boldsymbol{\eta}$ :

$$\displaystyle  E^{tot}(\boldsymbol{u},\boldsymbol{\eta}) =  E^{0} + E^{phonon}(\boldsymbol{u}) + E^{elastic}(\boldsymbol{\eta}) + E^{coupling}(\boldsymbol{u},\boldsymbol{\eta}) $$

The methodology followed in MULTIBINIT consists in making a Taylor expansion around the RF, which is assumed to be a stationary point of the BO energy surface. As such, the energy expression can but further decomposed as follows :

$$\displaystyle  E^{tot}(\boldsymbol{u},\boldsymbol{\eta}) =  E^{0} + [ E^{phonon}_{harm}(\boldsymbol{u})  + E^{elastic}_{harm}(\boldsymbol{\eta})  + E^{coupling}_{harm}(\boldsymbol{u},\boldsymbol{\eta})] + [E^{phonon}_{anharm}(\boldsymbol{u}) + E^{elastic}_{anharm}(\boldsymbol{\eta}) + E^{coupling}_{anharm}(\boldsymbol{u},\boldsymbol{\eta})] $$

The first term $E_O$ is the energy of the RS, which has be be fully relaxed (e.g. [[ionmov]]=2) with very strict tolerance criterium ([[tolmxf]] < 1E-7) since we assume all first energy derivatives are zero.  This $E_0$ energy has to be included in the global DDB file, by including the ground-state DDB when merging all partial DDBs with [[lesson:rf2| mrgddb]] 
The coefficients of the set of harmonic terms correspond to various second derivatives of the energy respect to atomic displacements and macroscopic strains. They can be directly computed with ABINIT using DFPT ([[lesson:rf1| Phonon response]], [[lesson:elastic|strain response]], [[lesson:ffield|electric field response]]) and used as parameters of our model. As such, our second-principles model reproduces exactly the first-principles results at the harmonic level (i.e. full phonon dispersion curves, elastic and piezoelectric constants of the RS). In pratcice, the global DDB file produced by ABINIT is so used as an input file for MULTIBINIT containing all the harmonic coedfficients.  This file must contain second energy derivatives respect to (i) all atomic displacements ([[rfphon]] 1; [[rfatpol]] 1 natom; [[rfdir]] 1 1 1) on a converged grid of q-points (defining the range of interactions in real space), (ii) macroscopic strains ([[rfstrs]] 3; [[rfdir]] 1 1 1) and also for insulators (iii) electric fields ([[rfelfd]] 1; [[rfdir]] 1 1 1) in order to provide the Born effective charges and dielectric constant used for the description of long-range dipole-dipole interactions.

The coefficients of the set of anharmonic terms correspond to higher-order derivatives of the energy respect to atomic displacements and macroscopic strains. They are numerous and not computed individually at the first-principles level. Instead, the most important terms will be selected by MULTIBINIT and related coefficients fitted in order to reproduce the BO energy surface. To that end, a training set (TS) of ABINIT data needs to be provided on which the fit will be realized. This TS consists in a set of of atomistic configurations realized on a suitable supercell depending on the range of anharmonic interactions (typically 2x2x2 supercell) and for which energy, forces and stresses are provided. This takes the form of an ABINIT netcdf "_HIST.nc" file. Providing an appropriate TS, properly sampling the BO surface, is crucial to obtain an appropriate model. How to built it depends on the kind of system (stable or with instabilities) and will not be further discussed here. 

In summary, constructing a second-principles lattice model with MULTIBINIT requires two input files which are direct output of ABINIT : (i) a full "DDB" file containing the reference energy and second energy derivatives which correspond to harmonic coefficients of the model and (ii) a "_HIST.nc" file containing the energy, forces and stresses of an appropriate training set of configurations from which the anharmonic terms will be selected and fitted. 

For this tutoral both these files will be provided.  


## 2 Fitting procedure: creating anharmonicities

In this tutorial, we will take as an exemple of a material without lattice instabilities: the perovskite $\mathrm{BaHfO_3}$ in its cubic phase.

*Optional exercise $\Longrightarrow$ Compute the phonon band structure with [[help:anaddb|Anaddb]]. You can download the complete DDB file (resulting from the previous calculations) here:*
{% dialog tests/tutomultibinit/Input/tmulti6_DDB %}

**Before starting, you might to consider working in a different subdirectory than for the other lessons. Why not create "Work_fitFatticeModel"?**

The file ~abinit/tests/tutomultibinit/Input/tmulti1_1.files lists the file names and root names.
You can copy it in the **Work_fitLatticeModel** directory and look at this file content, you should see:

      tutomulti6_1.abi
      tutomulti6_1.abo
      tmulti6_DDB
      no
      tmulti6_HIST.nc
      no

As mentioned in the guide of [[help:multibinit | MULTIBINIT]]:

   * "tutomulti6_1.abi" is the main input file
   * "tutomulti6_1.abo" is the main output
   * "tmulti6_DDB" is the DDB which contains the system definition and the list of the energy derivatives
   * "tmulti6_HIST.nc" is the set of DFT configurations to fit

It is now time to copy the file ~abinit/tests/tutomultibinit/Input/tutomulti6_1.abi, ~abinit/tests/tutomultibinit/Input/tmulti6_DDB and  tmulti6_HIST.nc in your **Work_fitLatticeModel** directory.
You should read carefully the input file:

{% dialog tests/tutomultibinit/Input/tmulti6_1.abi %}

You should now run (it would take less than 15 minutes):

    mpirun -np 10 multibinit < multi6_1.files > tmulti6_1_stdout& 
    
The resulting output file, tmulti1_1.abo, should be similar to the one below.
{% dialog tests/tutomultibinit/Refs/tmulti1_1.abo %}

The fitted anharmonocites are stored in tmulti6_1_coeffs.xml and informations about the differences between the the DFT data and the model are stored in "TRS\_fit\_diff\_energy.dat" and "TRS\_fit\_diff\_stress.dat". The global information about the reproduction of the DFT data is written in the output file.
Before the fit, the goal function is equal to:

```
Goal function values at the begining of the fit process (eV^2/A^2):
   Energy          :   2.2780100662032291E-03
   Forces+Stresses :   3.9733621903060741E-02
   Forces          :   2.4711894644538119E-02
   Stresses        :   1.5021727258522627E-02
```

After adding the anharmonicities, the goal function value is equal to:


```
Goal function values at the end of the fit process (eV^2/A^2):
   Energy          :   5.3930653437262267E-05
   Forces+Stresses :   6.3907163407037606E-03
   Forces          :   3.9644715403426073E-03
   Stresses        :   2.4262448003611538E-03
```

*Optional exercise $\Longrightarrow$ Try to fit on all irreducible atoms with [fit\_iatom](fit_iatom) = 0. This procedure is time consumming (around 15 min). You can also play with [fit\_cutoff](fit_cutoff) to see if there is other terms selected. *

## 3 Bounding of the model

Since the approach of the procedure is based on a polynomial expansion of the energy, it is common that the produced model is diverging at high temperature. In order to avoid this divergence, we will produce additional terms (order 6 and 8 terms) that assure the bounding of the model.

**Before starting, you might to consider working in a different subdirectory than for the other lessons. Why not create "Work_boundingLatticeModel"?**

The file ~abinit/tests/tutomultibinit/Input/tmulti7\_1.files lists the file names and root names.
You can copy it in the **Work_boundingLatticeModel** directory and look at this file content, you should see:

      tutomulti7_1.abi
      tutomulti7_1.abo
      tmulti6_DDB
      tmulti7_1_coeffs.xml
      tmulti6_HIST.nc
      no

"tmulti7\_1\_coeffs.xml" is the model that you produced in the previous section and we will bound.

It is now time to copy the file ~abinit/tests/tutomultibinit/Input/tutomulti7\_1.abi, ~abinit/tests/tutomultibinit/Input/tmulti6\_DDB, tmulti7\_1\_coeffs.xml and tmulti6\_HIST.nc in your **Work_boundingLatticeModel** directory.
You should read carefully the input file:

{% dialog tests/tutomultibinit/Input/tmulti7_1.abi %}

You should now run (it would take less than 1 minute):

    multibinit < multi7_1.files > tmulti7_1_stdout&
    

    
## 4 Running an molecular dynamics with the effective model

The aim of the construction of effective model is to be able to run realistic molecular dynamics. 

The file ~abinit/tests/tutomultibinit/Input/tmulti8\_1.files lists the file names and root names.
You can copy it in the **Work_MDLatticeModel** directory and look at this file content, you should see:


      tutomulti8_1.abi
      tutomulti8_1.abo
      tmulti6_DDB
      tmulti8_1.xml
      no
      no

"tmulti8\_1\_coeffs.xml" is the model that have been bounded in the previous step.

It is now time to copy the file ~abinit/tests/tutomultibinit/Input/tutomulti7\_1.abi, ~abinit/tests/tutomultibinit/Input/tmulti6\_DDB, tmulti7\_1\_coeffs.xml and tmulti6\_HIST.nc in your **Work_MDLatticeModel** directory.
You should read carefully the input file:

{% dialog tests/tutomultibinit/Input/tmulti8_1.abi %}

You should now run (it would take less than 2 minutes):

    multibinit -np 10 < multi8_1.files > tmulti8_1_stdout&
    
You can visualize your dynamics with the agate software:

    agate tmulti1_3.abo_HIST.nc
    
*Optional exercise $\Longrightarrow$ Try to recover the phase transition of $\mathrm{SrTiO_3}$ (DDB is located in ~abinit/tests/tutomultibinit/Input/tutomulti9_1.ddb and the model in ~abinit/tests/tutomultibinit/Input/tmulti9_1.xml).*


![Schema 1](multibinit_assets/HeatingRot.pdf) 


<!--In this tutorial, we will take as an example of a material without lattice instabilities: the perovskite CaTiO$_3$ in its $Pnma$ phase.

*Optional exercise $\Longrightarrow$ Compute the phonon band structure with [[help:anaddb|Anaddb]]. You can download the complete DDB file (resulting from the previous calculations) here:*
{% dialog tests/tutomultibinit/Input/tmulti1_DDB %}


**Before starting, you might to consider working in a different subdirectory than for the other lessons. Why not create "Work_latticeModel"?**

The file ~abinit/tests/tutomultibinit/Input/tmulti1_1.files lists the file names and root names.
You can copy it in the **Work_latticeModel** directory and look at this file content, you should see:

      tutomulti1_1.abi
      tutomulti1_1.abo
      tmulti1_DDB

As mentioned in the guide of [[help:multibinit | MULTIBINIT]]:

   * "tutomulti1_1.abi" is the main input file
   * "tutomulti1_1.abo" is the main output
   * "tmulti1_DDB" is the DDB which contains the system definition and the list of the energy derivatives

It is now time to copy the file ~abinit/tests/tutomultibinit/Input/tmulti1_DDB and ~abinit/tests/tutomultibinit/Input/tmulti1_DDB in your **Work_latticeModel** directory.
You should read carefully the input file:

{% dialog tests/tutomultibinit/Input/tmulti1_1.abi %}

You should now run (it would take less than a second):

    multibinit < tmulti1_1.files > tmulti1_1_stdout

The resulting output file, tmulti1_1.abo, should be similar to the one below.
{% dialog tests/tutomultibinit/Refs/tmulti1_1.abo %}


The run you performed was aimed at reading the DDB file, generating the short range interatomic force constants and extract all the other informations related to the harmonic part of the model.
You can find inside the output file, the Born effective charges, the clamped-ion elastic tensor and the internal strain coupling parameters. Take some time to open and read the tmulti1_1.abo file.
If the DDB file is complete, the generation of the XML file requires only few input variables:

   * [[multibinit:prt_model]] = 1 $\Longrightarrow$ the generation of the XML file is activated, takes the time to read the possible options for [[multibinit:prt_model]].
   * [[multibinit:ngqpt]]    = 2 2 2 $\Longrightarrow$ specified the q-point mesh included in the tmulti1_DDB file
   * [[multibinit:dipdip]]   = 0 $\Longrightarrow$  it disables the computation of the dipole-dipole interaction: we don''t need to compute it to generate the XML file, MULTIBINIT will be able to regenerate the dipole-dipole interaction in a next run. We remind you that the default for this option is 1 and in the most part of your runs you will not use this option.

After this run, you should see in your directory tmulti1_1_model.xml, you can take some time to open and read this file. As you can see, it contains all the informations about the system definition (lattice parameters, atomic positions) and the data for the harmonic part of the potential.

Your XML file is now ready and you can use it as input for MULTIBINIT. To do that, copy now in your work directory the file ~abinit/tests/tutomultibinit/Input/tmulti1_2.files; you should see inside it:

      tutomulti1_2.abi
      tutomulti1_2.abo
      tmulti1_1_model.xml

Here, the DDB file is replaced by the XML file. Do not forget to copy the ~abinit/tests/tutomultibinit/Input/tutomulti1_2.abi in your directory before you run:

    multibinit < tmulti1_2.files > tmulti1_2_stdout

In tutomulti1_2.abi, [[multibinit:prt_model]] is still set to one so multibinit will write again the model XML file, which is useless at this stage, being a copy of the one read as input. Set this input variable to 0 and, in this case, MULTIBINIT will just read (and not write) the XML file.

With the two last examples, we have shown that MULTIBINIT is able to read either a DDB file or a XML as inputs for the system definition and the harmonic part of the potential.

We can now run our *first dynamics*: you can copy the files ~abinit/tests/tutomultibinit/Input/tutomulti1_3.* in your work directly and have a look them.

{% dialog tests/tutomultibinit/Input/tmulti1_3.abi %}

The simulation starts from the DDB to correctly account for the dipole-dipole interactions. You can visualize your dynamics with the agate software:

    agate tmulti1_3.abo_HIST.nc

Also try to use the effective potential from the xml file instead, in which the dipole-dipole interactions were not corrected. What do you see when you visualize the track?

* * *

This MULTIBINIT tutorial is now finished.-->
