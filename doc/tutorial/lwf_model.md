---
authors: xuhe
---

# Lattice Wannier function (LWF) model in MULTIBINIT


## Build a LWF model and run LWF dynamics in MULTIBINIT

This lesson aims at showing how to build a LWF model and run a spin dynamics calculation.

**Before beginning, we recommend the reading of the theory on Lattice Wannier functions in the literature [[cite:rabe1995]]. The construction of LWF from the phonon band structure is covered in another [tutorial](/tutorial/lattice_wannier). **

With this lesson, you will learn to:

  * Run the dynamics of LWF.

[TUTORIAL_README]

*Before trying to run the LWF dynamics, you might consider to work in a subdirectory for this tutorial. Why not Work_lwfdyn?*


### Basic: how to use MULTIBINIT to run LWF dynamics

With the LWF model built and saved into a netcf file, we can run LWF dynamics calculation with MULTIBINIT.  The input files for an example can be found at  ~abinit/tests/tutomultibinit/Input/ .
The input file can also be downloaded from:
{% dialog tests/tutomultibinit/Input/tmulti6_1.abi %}

The netcdf file tmulti6_1.nc contains the LWF model,  and tmulti6_1.abi  is the input file. The input file is shown below. Note that the first line "lwf_pot_fname"


```
ncellmat = 8 0 0  0 8 0 0 0 8 # size of supercell, 3x3 matrix.
lwf_dynamics=3               # type of dynamics. 3: Berendsen thermostat
lwf_dt =1e-16 s              # time step.
lwf_ntime = 3000             # number of time steps
lwf_nctime = 100             # write to netcdf history file per 100 steps.
lwf_taut = 1e-14 s           # characteristic time for berendsen dynamics.
lwf_init_state=1             # initial state. 1: random amplitude.

lwf_var_temperature=1        # whether to do variable temperature calculation.
lwf_temperature_start=600.1  # Starting temperature.
lwf_temperature_end=0.1      # Ending temperature.
lwf_temperature_nstep=7      # number of temperature steps.
```

For a more realistic simulation, the size of the supercell, and the number of temperature steps should be increased. 


To run the LWF dynamics, use the command below:

```
multibinit tmulti6_1.abi --F03 > tmulti6_1.log
```

For each temperature point, the temperature, kinetic energy, potential energy, and
total energy will be outputed for every lwf_nctime steps. 

```
 Iteration     temperature(K)        Ekin(Ha/uc)        Epot(Ha/uc)        ETOT(Ha/uc)
      100          601.99469        1.90641E-03        3.70420E-04        2.27683E-03
      200          419.24864        1.32768E-03        2.73706E-03        4.06475E-03
      300          421.02202        1.33330E-03        5.28296E-03        6.61626E-03
      ......
```

A file ends with "_lwfhist.nc", which contains the trajectory of the LWF amplitudes, will be generated. By analyzing it, we 
can extract the information like heat capacity and susceptibility. If there is structural phase transition, the critical temperature can also 
be found. The tools for doing these analysis are still under development. 

