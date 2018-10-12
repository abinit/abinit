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

Now that we have the spin model xml file, we can run a spin dynamics calculation with multibinit.  Example input files can be found at ~abinit/tests/tutomultibinit/Input/tmulti5.* .  There are three files: 

* "tmulti5.files" is the "files" file, which gives the names of the input and output files for  Multibinit.

* "tmulti5.xml" is the file containing the Heisenberg model parameters.

* "tmulti5.in" is the main input file, where you can put the parameters for the spin dynamics simulation.

You can copy these three files into a directory (e.g. tutor_spindyn). 

In tmulti5.files, three file names are given: 

```
tmulti5.in
tmulti5.out
tmulti5.xml
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
spin_dt=1e-16 s                ! Time step
spin_qpoint = 0.0 0.0 0.0      ! Wave vector for summation of spin in each sublattice.
```

To run spin dynamics with multibinit,

```
cd tutor_spindyn
multibinit < tmulti5.files > tmulti5.txt
```

After that, an output file named tmulti5.out and a netcdf file tmulti5.out_spinhist.nc will be found.  

(<!--TODO: Impove the output. Energy should be in Ha instead of eV. average M should be for each sublattice (or not)-->: )

In the .out file, lines below can be found, which gives a overview of the evolution of the system with time. 

```
 Begining spin dynamic steps :
  ==================================================================================
     Iteration          time(s)        average M           Energy
  -----------------------------------------------------------------------------------
           100      9.90000E-15      6.38643E-01     -6.34928E-17
           200      1.99000E-14      5.49196E-01     -5.56385E-17
           300      2.99000E-14      5.20451E-01     -5.41959E-17
           400      3.99000E-14      4.94333E-01     -5.25889E-17
           500      4.99000E-14      4.93464E-01     -5.23715E-17
......
```

In the netcdf file, the trajectories of the spins can be found. They can be further analyzed using postprocessing tools. 

##### To add

- convergence wrt. time step (energy conservation), box size
- observables: magnetization versus temperature, determine Curie temperature
- spectral function
- magnetic susceptibility
- comparison with experiments?
- alternating spins for antiferromagnets

## 5. Postprocessing






