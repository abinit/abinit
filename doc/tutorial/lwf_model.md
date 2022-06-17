---
authors: xuhe
---

# Lattice Wannier function (LWF) model in MULTIBINIT


## Build a LWF model and run LWF dynamics in MULTIBINIT

This lesson aims at showing how to build a LWF model and run a spin dynamics calculation.

**Before beginning, we recommend the reading of the theory on Lattice Wannier functions in the literature [[cite:rabe1995]]. The construction of LWF from the phonon band structure is covered in another [tutorial](/tutorial/lattice_wannier). 

With this lesson, you will learn to:

  * Run the dynamics of LWF.

[TUTORIAL_README]

*Before beginning, you might consider to work in a subdirectory for this tutorial. Why not Work_lwfdyn?


### Basic: how to use MULTIBINIT to run LWF dynamics

With the LWF model built and saved into a netcf file, we can run LWF dynamics calculation with MULTIBINIT.  The input files for an example can be found at  ~abinit/tests/tutomultibinit/Input/
The netcdf file tmulti6.nc contains the LWF model,  and tmulti6.abi  is the input file. The input file is shown below. Note that the first line "lwf_pot_fname"




```
lwf_pot_fname = "tmulti6.nc"
ncellmat = 0 2 0  2 0 0 0 0 2
lwf_dynamics=3
lwf_dt =1e-16 s
lwf_ntime = 3000
lwf_nctime = 100
lwf_taut = 1e-14 s
lwf_init_state=1

lwf_var_temperature=1
lwf_temperature_start=600.1
lwf_temperature_end=0.1
lwf_temperature_nstep=7
```


### Basic: how to use MULTIBINIT to run spin dynamics

Once we have the spin model xml file, we can run a spin dynamics calculation with MULTIBINIT. Example input files can be found at ~abinit/tests/tutomultibinit/Input/tmulti5_1.* .  There are three files:

* "tmulti5_1.files" is the "files" file, which gives the names of the input and output files for  MULTIBINIT.
* "tmulti5_1.abi" is the main input file containing the parameters for the spin dynamics simulation.
* "tmulti5_1.xml" is the file containing the Heisenberg model parameters. 

You can copy these three files into a directory (e.g. Work_spindyn).

In tmulti5_1.files, three file names are given:

```
tmulti5_1.abi
tmulti5_1.abo
tmulti5_1.xml
```

which gives the input, output and xml file names. The file tmulti5_1.xml contains the $J_{ij}$ values for a simple toy system which has a cubic lattice and one atom per unit cell. Its critical temperature is around 600K.

In tmulti5_1.abi, the variables for running a spin dynamics calculation are given:

```
prt_model = 0
ncell =   16 16 16              ! number of unit cells in supercell

spin_dynamics = 1               ! switch on spin dynamics
spin_init_state = 2             ! ferromagnetic initial state

spin_temperature = 600          ! temperature of spin (Kelvin)
spin_ntime_pre = 10000          ! time steps for thermolization
spin_ntime = 20000              ! time steps for measurement
spin_nctime = 100               ! Number of time steps between two writes
                                ! into netcdf
spin_dt = 1e-16 s               ! Time step (seconds)
spin_projection_qpoint= 0.0 0.0 0.0       ! Wave vector for summation of spin in each
                                ! sublattice.

spin_write_traj = 0             ! do not write spin trajectory to netcdf file
```

To run spin dynamics with MULTIBINIT

```
cd Work_spindyn
multibinit --F03 < tmulti5_1.files > tmulti5_1.txt 
```

Note that the .files file will be deprecated in the next version of ABINIT and MULTIBINIT. Then only two files are required. The following variables in the input file can be used to specify the spin potential file and the prefix of the output files. 

```
spin_pot_fname = "tmulti5_1.xml"
outdata_prefix = "tmulti5_1.abo"
```

To run the spin dynamics without using the files file, 

```
multibinit tmulti5_1.abi --F03 > tmulti5_1.txt 
```

After the calculation is done, you will find an output file named tmulti5_1.abo and a netcdf file tmulti5_1.abo_spinhist.nc.

With the LWF model built and saved into a 

