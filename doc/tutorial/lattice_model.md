---
authors: AM
---

# First lesson on Multibinit

## How build a lattice model at the harmonic level and run a dynamics.  

This lesson aims at showing how to build a harmonic model by using a second-principles approach for lattice dynamics simulations based on atomic potentials fitted on first-principles calculations.

**Before beginning, it is very important to read the reference [[cite:Wojdel2013]].**

With this lesson, you will learn how to :

  * Compute the model from the DDB
  * Generate the XML for the model 
  * Run a dynamics within Multibinit

In this tutorial, all the knowledge about the density-functional theory (DFT) and  density-functional perturbation theory (DFPT) should be already acquired.
Morever, the DFPT is the key feature of ABINIT for multibinit. In order to learn  the DFPT and the associated codes Mrgddb to produce the DDB file,
please see the tutorials on [[lesson:rf1| Phonon response]] and [[lesson:elastic|strain response]] and [[lesson:rf2| Mrgddb]] before to continue.
After theses tutorials you should be able to compute a full DFPT calculation and a complete DDB file.
Thereby this tutorial will not provide the inputs for abinit.
  
The AGATE software is also required for this tutorial, you can install it on debian with:

    sudo add-apt-repository ppa:piti-diablotin/abiout
    sudo apt-get update && sudo apt-get install abiout

## 1 The Harmonic part of the lattice model

As mentionned in [[cite:Wojdel2013]], define the reference of your model is the starting point in the construction of the model:

![Schema 1](lattice_model_assets/reference.png)

The choice of this reference is important because we will construct the model based on an expansion around this reference structure.
Once this choice is done, make sure than your system is fully relaxed and compute a single DFT calculation on this system in order to generate a DDB file.
In this tutorial, we will take as an example of a material without unstabilities, the perovskite CaTiO3 in the Pnma phase.

Now from this reference, we consider the set of perturbation, the atomic displacements $u$ and the strain $\eta$:

![Schema 2](lattice_model_assets/deformation.png)

At **the harmonic part level**, we can express the potential as a sum of three contributions (with respect to the set of pertubations ($u$,$\eta$)):

  $$\displaystyle E^{harm}(u,\eta) =  E^{harm}(u) + E^{harm}(u,\eta) + E^{harm}(\eta)$$

  * A Computation of the phonon response (short range + dipole-dipole interaction):
      $$\displaystyle E^{harm}(u) => \underbrace{\frac{\partial^2 E}{\partial
          u^2}}_{\substack{\text{Inter-atomic}\\\text{force constants}}} +
     \underbrace{\frac{\partial^2 E}{\partial
          {\cal{E}}^2}}_{\text{Dielectric tensor}} +
      \underbrace{\frac{\partial^2 E}{\partial{\cal{E}} \partial u}}_{\text{Effective charges}} $$
  
  * A Computation of the strain response:
$$\displaystyle  E^{harm}(\eta) =>\underbrace{\frac{\partial^2 E}{\partial
            \eta^2}}_{\text{Elastic constants}} $$

  * A Computation of the strain-phonon coupling:
$$\displaystyle E^{harm}(u,\eta) =>
      \underbrace{\frac{\partial^2
            E}{\partial\eta \partial u}}_{\substack{\text{Internal strain}}} $$

We note that for the harmonic part, All the needed quantities can be computed with the DFPT features of abinit and the resulting file is called the DDB file.
In the case of insulator, the dipole-dipole interaction will be recomputed directly within Multibinit.
Thereby you need to provide into the DDB file the clamped ion dielectric tensor and the Born effective charges.
Don't forget to add the in the final DDB file the DDB from the single DFT calculation done on the reference system (you can still use mrgddb).

Moreover, the ddb file is also an output of a DFT calculation. In the case of the generation of a model with multibinit, it is important to merge into the final DDB file,

In this tutorial, we will take as an example of a material without unstabilities, the perovskite CaTiO3 in the Pnma phase.

*Optional exercice => Compute the phonon band structure with Anaddb. you can download the complete DDB file (resulting from the previous calculations) here:*
{% dialog tests/tutomultibinit/Input/tmulti1_DDB %}


*Before beginning, you might consider to work in a different subdirectory as
for the other lessons. Why not create "Work_latticeModel"*


The file ~abinit/tests/tutomultibinit/Input/tmulti1.files lists the file names and root names.
You can copy it in the Work_latticeModel directory and look at this file, you should see:

      tutomulti1_1.in
      tutomulti1_1.out
      tmulti1_DDB

As mentionned in the [[help:multibinit | guide of multibinit]]:

   * tutomulti1_1.in is the main input file
   * tutomulti1_1.out is the main output
   * tmulti1_DDB is the DDB which contains the system definition and a list of derivatives

It is now time to copy the file ~abinit/tests/tutomultibinit/Input/tmulti1_DDB and ~abinit/tests/tutomultibinit/Input/tmulti1_DDB in your work_latticeModel directory. 
You should read carefully the input file:

{% dialog tests/tutomultibinit/Input/tmulti1_1.in %}

You now should make the run (less than a seconds):

  ./multibinit < tmulti1.files > tmulti1_1_stdout

The resulting main output file, trf1_1.out, should be similar to the one below.
{% dialog tests/tutomultibinit/Refs/tmulti1_1.out %}


The aim of the run that you just made was to read the DDB file, generate the short range interatomic force constants and extract all the other informations for the harmonic part of the model.
You can see into the output file, the Born effective charges, The clamped ion elastic tensor and the internal strain coupling parameters. Take some time to open and read the tmulti1_1.out file.
Once the DDB file is complete, the generation of the XML file within multibinit requires only few input variables:

   * [[multibinit:prt_model]] = 1     => active the generation of the XML file, takes the time to read the possible options for [[multibinit:prt_model]].
   * [[multibinit:ngqpt]]    = 2 2 2 => specified the q-point mesh included in the tmulti1_DDB file
   * [[multibinit:dipdip]]   = 0     => here the computation of the dipole-dipole interaction is disable because we don't need to compute it to generate the XML, indeed in multibinit will be able to regenerate the dipole-dipole interaction in a next run. Let's remind you that the default for this option is 1, so in many of you multibinit run you will never use this input.

After this run, you should see in your directory tmulti1_1_model.xml, you can take some time to open and read this file. You will find all the informations about the system definition (lattice parameters, atomic positions) and the data for the harmonic part of the potential.

You XML file is now generated you can now use it as input for multibinit. To do that, open copy now in your work directory the files ~abinit/tests/tutomultibinit/Input/tmulti2.files; you should see inside:

      tutomulti1_2.in
      tutomulti1_2.out
      tmulti1_1_model.xml

Here the DDB file is remplace by the XML file, do not forget to copy the ~abinit/tests/tutomultibinit/Input/tutomulti1_2.in in your directory and run:

   ./multibinit < tmulti2.files > tmulti1_2_stdout
  
In tutomulti1_2.in, [[multibinit:prt_model]] is still set to one so multibinit will generate again a new XML which is not useful. You can set this input variable to 0 but in this case, multibinit will just read the XML.

With the two last examples, we showed that multibinit is able to read either the DDB file or the XML as input for the system definition and the harmonic part of the potential.

We can now run our first dynamic within multibinit, you can copy into you directory the files ~abinit/tests/tutomultibinit/Input/tutomulti1_3*. and have a look the

{% dialog tests/tutomultibinit/Input/tmulti1_3.in %}
  
* * *

This MULTIBINIT tutorial is now finished. You are advised to read now the [[lesson:fit_process|second tutorial on Multibinit]]
