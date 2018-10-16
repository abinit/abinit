---
authors: xuhe
---

# Spin model in Multibinit

## Build a spin model and run spin dynamics in Multibinit

This lesson aims at showing how to build a spin model and run a spin dynamics calculation.

**Before beginning, it is important to know the theory of spin dynamics, which can be found in the literature (e.g.[[cite: Evans2014]] ).**

With this lesson, you will learn to:

  * Generate the XML for the model 
  * Run a spin dynamics calculation with MULTIBINIT

The TB2J python package, which can be used to generate a spin model, can be found on the abinit gitlab website at https://gitlab.abinit.org/xuhe/TB2J. This package will be included in the Abinit package in the future. 



[TUTORIAL_README]

## 1 Heisenberg Model formalism

The spin model, as implemented in Multibinit, is defined as a classical Heisenberg model. In the current version of Multibinit, we consider the following interactions: exchange interaction, single ion anisotropy (SIA), Dzyaloshinski-Moriya (DM) interaction, external magnetic fields. The total energy then reads as 

$$E = E^{exc}+E^{SIA} + E^{DM}+E^{ext}$$

The exchange energy $E^{exc}$ can be written as

$$E^{exc} =- \sum_{i\neq j} J_{ij} \vec{S_i}\cdot{\vec{S_j}}$$ 

A few conventions used in Multibinit should be noted:

- all the $\vec{S}$ are normalized to 1. 
- Both $J_{ij}$ and $J_{ji}$ are in the Hamiltonian.
- There is a minus sign, which means that the interaction is ferromagnetic for $J_{ij} >0$.

As the sites $i$ and $j$ defined in the model are in a finite cell (often the primitive cell), there are interactions between sites in the same cell and between sites in two different cells. The position vector of a site $j'$ in a different cell than site $i$ is denoted as $\vec{r}_{j'}=\vec{r}_j+\vec{R}$ with $\vec{R}$ being a combination of lattice vectors. For a site $j$ in the same cell as site $i$ the lattice vector $\vec{R}$ is $(0,0,0)$. Due to translation symmetry we can choose the lattice vector for site $i$ to be $\vec{R}=\vec{0}$. Hence, we denote the Heisenberg coefficients as $J_{ij}(\vec{R})$ and drop the prime for the sites in different cells.

The SIA term can be written as 

$$E^{SIA}=-k_u \sum_i (\vec{S_i}\cdot \vec{e})^2,$$

where $k_u$ and $\vec{e}$ are the amplitude and direction of the single ion anisotropy.

The DM term can be written as

$$ E^{DM} = \sum_{i\neq j} \vec{D}_{ij}\cdot \vec{S_i}\times{\vec{S_j}},$$

where $\vec{D_{ij}}$ is the amplitude of the DM interaction.

The external magnetic field term can be written as

$$E^{ext}=- \sum_i   m_{i}   \vec{S_i}\cdot \vec{H},$$

where $m_i$ denotes the magnetic moment of site $i$, and $\vec{H}$ is the magnetic field.

## 2. Format of xml file for spin model (TODO: remove?)

## 3. Build spin model file

One way to calculate the Heisenberg model parameters is to use the spin force theorem (see [[cite:Liechtenstein1983]], [[cite:Katsnelson2000]]), which take rotations of localized spins as perturbations. In Abinit, the Hamiltonian uses plane waves as a basis set, thus the localized spin is not directly accessible. We can construct localized Wannier functions rewrite the Hamiltonian in the Wannier basis. Then, the exchange parameters can be calculated from this Hamiltonian ( [[cite:Korotin2015]] ). 

For building the Wannier function Hamiltonian from Abinit, see the tutorial [wannier90](wannier90). Other DFT codes interfaced with [Wannier90](http://www.wannier.org) can also be used. Then, the  [TB2J](https://gitlab.abinit.org/xuhe/TB2J) package can be used to get the Heisenberg model parameters and generate the input model for Multibinit. 

## 4. Run spin dynamics

### Basic

Now that we have the spin model xml file, we can run a spin dynamics calculation with multibinit.  Example input files can be found at ~abinit/tests/tutomultibinit/Input/tmulti5_1.* .  There are three files: 

* "tmulti5.files" is the "files" file, which gives the names of the input and output files for  Multibinit.

* "tmulti5.xml" is the file containing the Heisenberg model parameters.

* "tmulti5.in" is the main input file, where you can put the parameters for the spin dynamics simulation.

You can copy these three files into a directory (e.g. tutor_spindyn). 

In tmulti5.files, three file names are given: 

```
tmulti5_1.in
tmulti5_1.out
tmulti5_1.xml
```

which gives the input, output and xml file names. 

In tmulti5.in, the variables for running spin dynamics are given:

```
prt_model = 0
ncell =   15 15 15             ! number of unitcell in supercell

spin_mag_field= 0.0 0.0 0.0    ! external magnetic field (Tesla)
spin_dynamics=1	               ! switch on spin dynamics
spin_temperature = 600         ! temperature of spin. (Kelvin)
spin_ntime =10000              ! Total number of time steps.
spin_nctime=100                ! Number of time steps between two writes into netcdf file
spin_dt=1e-16 s                ! Time step, unit s.  
spin_qpoint = 0.0 0.0 0.0      ! Wave vector for summation of spin in each sublattice.
```

To run spin dynamics with multibinit,

```
cd tutor_spindyn
multibinit < tmulti5.files > tmulti5.txt
```

After that, an output file named tmulti5.out and a netcdf file tmulti5.out_spinhist.nc will be found.  

In the .out file, lines below can be found, which gives a overview of the evolution of the system with time. 

```
   Begining spin dynamic steps :
 =================================================================
     Iteration          time(s)       Avg_Mst/Ms      Energy (Ha)
 -----------------------------------------------------------------
           100      9.90000E-15      6.38643E-01     -1.45634E+01
           200      1.99000E-14      5.49196E-01     -1.27619E+01
           300      2.99000E-14      5.20451E-01     -1.24310E+01
           400      3.99000E-14      4.94333E-01     -1.20624E+01
           500      4.99000E-14      4.93464E-01     -1.20125E+01
......
```

Here, the Avg_mst ($||<m_i e^{2\pi \vec{q}\cdot\vec{R_i}}>||$) means the average staggered magnetic moment, Ms is the saturate magnetic moment . If all the spins align to the wavevector of ($\vec{q}$) [[multibinit:spin_qpoint]], this value would be 1.0. And it would deviate from 1.0 with thermo fluctuation.

At the end of the run, there is a summary of the 

In the netcdf file, the trajectories of the spins can be found. They can be further analyzed using postprocessing tools. 

It is essential to choose the parameters so that the calculation result is meaningful. 

* time step ([[multibinit: spin_dt]]):

    Typical time step is about $10^{-15}  $ to $10^{-17}$ s. A optimal time step can be got by tring several values and comparing the results (equilibrium magnetic order, moments, etc) with a calculation with small time step. At this stage, a small box and temperature close to zero can be used.   (<TODO: there must be a better way.>)

* supercell size ([[multibinit:ncell]])

  Due to the periodic boundary condition, the spins between periods could be correlated with each other, which could lead to artificially higher phase transition temperature, etc. A small box is also not enough for sampling of some quantities and the To avoid this, it is required to 

  For antiferromagnetic structure, or more generally, structure with non-zero wave vector, the box size should allow the spins to align with the q-vector, i.e. ($\vec{q}\cdot \vec{n}$) should be integers. For some structure, it is not easy or sometimes impossible to find such $\vec{n}$. In these cases, a large box is usually required.  

* Total time ([[multibinit: spin_nstep]])

  The total should allow the spins to relax to the equilibrium state, and in order to calculate some observables, longer time are required so enough ensembles can be generated.  To see how much time is needed for the system to get to the equilibrium state, we can plot the magnetic moment. 

### Temperature dependent calculation

A most common usage of spin dynamics is to calculate the quantities (e.g. magnetic moments, susceptibility, specific heat ) as functions of temperature and find the phase transition temperature. By setting [[multibinit:spin_var_temperature]] to 1 and the starting, ending temperature, and number of steps as follows, a series of calculation will be carried out. (See ~abinit/tests/tutomultibinit/Input/tmulti5_2.* )

```
spin_var_temperature = 1       ! switch on variable temperature calculation
spin_temperature_start = 0.0   ! starting temperature
spin_temperature_end = 1500.0  ! ending temperature  
spin_temperature_nstep = 31    ! number of steps, both starting and ending T included.
```

After the run, the trajectories for each temperature will be written into the T0001_spin_hist.nc to T0031_spin_hist.nc files. 

There are several ways to find the critical temperature. The most natural way is to use the M-T curve. However, there are some difficulties due to that the change of mangetic moment is not abrupt, and it's sensitive to box size. There is divegences of  specific heat and magnetic susceptibility, which could be better for finding $T_c$. Another quantity is Binder cumulant, defined as $U_4= 1.0- \frac{<m^4>}{3 <m^2> }$. It is less sensitive to the box size and the change is abrupt. 



##### To add

- observables: magnetization versus temperature, determine Curie temperature
- spectral function
- magnetic susceptibility
- comparison with experiments?
- alternating spins for antiferromagnets

## 5. Postprocessing






