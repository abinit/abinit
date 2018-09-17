---
authors: AM,ACGC,FR
---

# The multibinit software

The MULTIBINIT software is using a second-principles approach for lattice dynamics simulations based on atomic potentials fitted on first-principles calculations [[cite:Wojdel2013]]. The second-principles effective potential accounts for harmonic (short-range and long-range dipole-dipole interactions) and anharmonic contributions, as well as the explicit treatment of homogeneous strain and its coupling with the lattice. On the one hand, parameters associated to second energy derivatives (harmonic interatomic force constants -IFC-, elastic constants, strain-phonon coupling, etc.) are determined exactly and provided by the Density Functional Perturbation Theory ([[topic:DFPT|DFPT]]). On the other hand, the anharmonic lattice contribution is restricted to a limited number of terms (short-range interaction and low order) and treated in a more effective way: it is fitted [[cite:Escorihuela-Sayalero2017]] to reproduce stresses and forces on a training set which should include representative configurations properly sampling the phase space.


## 0 Installation  

The MULTIBINIT software is included in the [[help:abinit|ABINIT]] package, thus to install this code, you can follow the instructions furnished by the [[help:../installation|installation guide]] of the package.

However, to be able to use all the MULTIBINIT features, you need to recompile ABINIT adding some flags in the configure (config.ac) file:
  
 * MULTIBINIT generates the output of the molecular dynamics in the NetCDF format: the NetCDF library is linked by using the flag:

        with_trio_flavor="netcdf"

* MULTIBINIT uses Libxml2 to parse the XML files (not mandatory, but more efficient for heavy XML files): this library is linked by using the flags:
  
        enable_xml=yes
        CFLAGS_EXTRA="-I/usr/include/libxml2"


## 1 How to run the code

### 1.1 Introducing the "files" file

Once prepared the input file (see below for the parameters description) and the required files for the generation of the model, the user musts to create a ".files" file which lists (one for each line) the filenames the job will require: the main input file, the main output file, the file for the model (model.DDB or model.XML), the XML file for the anharmonic part of the model and the NetCDF file for the training set.
The files file (called for example multibinit.files) looks like:
 
      multibinit.in
      multibinit.out
      model.DDB or model.XML
      model_anharmonic.XML
      training_set_HIST.nc

In this example:

  * The main input file is called "multibinit.in".
  * The main output will be written as "multibinit.out".
  * model.DDB or model.XML is the Database from ABINIT or XML.
  * model_anharmonic.XML (_optional_) is the XML with the coefficients from fitted polynomial.
  * training_set_HIST.nc is the history file in NetCDF format containing the training set.

The model.DDB (or model.XML) file contains the system definition and the list of the total energy derivatives with respect to three kind of perturbations: phonons, electric field and strain. The _optional_ XML file contains the list of coefficients obtained by the energy fitting by a polynome. The last file is mandatory to obtain the "model_anharmonic.XML" file.

### 1.2 Running the code

The main executable file is called "multibinit". Supposing that the executable is located in your working
directory, you can run it interactively (in Unix) with the command:

    multibinit < multibinit.files > log
  
or, in the background, with the command:

    multibinit < multibinit.files > log &

Here, the standard output and standard error are piped to the file called "log".

The user has the full freedom to change the filenames. Moreover, modifications of the above commands could be needed, depending on the UNIX flavour that is used on the platform supposed to execute the code.

The syntax of the input file is completely similar to the one of the main ABINIT: this file is parsed, keywords are identified, comments are also identified. However, the [[topic:multidtset|multidataset]] mode is not available.

## 2 Input variables
 
This Software is able to perform many different tasks. For example, it is possible to generate a second principle model by extracting the harmonic part from a DFPT calculation and fitting the related anharmonic part with a polynomial expression. Such a model can contain phonons instabilities and can be bound with an automatic procedure. In addition, for a given model, it is also possible to perform a molecular dynamics using MULTIBINIT. Each feature can be selected by a set of input variables for which the corresponding 'flag' variables activate the different tasks, e.g. [[multibinit:prt_model]], [[multibinit:fit_coeff]], [[multibinit:bound_model]] and [[multibinit:dynamics]]. The full list of the input variables are presented in the [[varset:multibinit|variable set]].

## 3 How to generate a model

Before to learn how to generate the model, we encourage you to read the papers [[cite:Wojdel2013]] and [[cite:Escorihuela-Sayalero2017]].

### 3.1 How to compute the harmonic contribution using [[help:abinit|ABINIT]]

MULTIBINIT requires a DDB file from ABINIT to build the harmonic part of the potential. To compute this file, you are first supposed to learn how to perform [[topic:DFPT|DFPT]] calculations with ABINIT. For this purpose, it is useful to follow the different tutorials: start with the second lesson about response functions ([[lesson:rf2]]) to compute the phonon dispersion and build the related DDB file, then follow the lesson on elasticity ([[lesson:elastic]]) to compute elastic contribution and build the related DDB. At the end, use the [[help:mrgddb|DDB merge tool]] to generate the "model.DDB" file seen before.

To learn the procedure to compute the **harmonic** part of the potential, you can follow this [[lesson:lattice_model | tutorial]].

### 3.2 How to compute the anharmonic contribution using a *training set* 

To compute the anharmonic part of the model, MULTIBINIT requires a training set containing several DFT calculations. Such training set will contain an ensemble of atomic configurations sampling the potential energy surface. This file will be used to fit a polynomial following the method developed by [[cite:Escorihuela-Sayalero2017]]. To generate the training set, the user requieres to get some experience in [[topic:MolecularDynamics|molecular dynamics]] simulation within ABINIT. Alternatively, the user can generate the training set configurations by means of phonon population procedure at fixed temperature (DDB file mandatory) activated by imposing [[abinit:ionmov]]=27.

### 3.3 How to bind a model

MULTIBINIT includes an automatic binding process.
In this [[lesson:fit_process | tutorial]], you can get the basic knowlegde to bind a model.

### 3.4 How to run a dynamics

MULTIBINIT uses routines of ABINIT which perform atomic relaxation, molecular dynamics, and Monte Carlo simulations.
In this [[lesson:lattice_model | tutorial]], you can get the basic knowlegde for these processes.
