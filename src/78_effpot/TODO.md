---
This file list the things to do for multibinit
---
# Todo list for Multibinit

## Genenal structure

### [TODO] Clean lattice part in 98_main/multibinit.F90 to manager or lattice code

State: 


### [TODO] Separate the 78_effpot directory into layers 

State: 

### [TODO] Make a example harmonic lattice potential. (hexu) [time: July, 2019]


#### [DONE] implement primitive potentail (read from netcdf file). Tested.
#### [DONE] implement supercell potential, fill_supercell
#### [DONE] add the data structure to the manager.

###[TODO] move random number generator to manager and use pointer in each mover.

---

## Documenation
### [TODO] Document general data structure  (hexu) [time: July, 2019]

State: ?

#### [TODO] Make a developer's guide on how to add new things. (hexu) [time: July, 2019]

State: First draft in multibinit_note.md
 
### [TODO?] Make a full technical note (including all the methods used) (everyone) 

State: One file exist for spin/spin-lattice coupling.

### [TODO] 

---

## Lattice

### [TODO] Adapt to general data structure

State: 

### [TODO] Mover

#### [TODO] set intial state:
       * [TODO] implement get_temperature and get_kinetic_energy
       * [TODO] More options other than Maxwell-Boltzmann distribution
       * [TODO] make center of mass fixed.
       * [TODO] make sure there is no rotation.

#### [TODO] implement movers. Langevin (NVT) and Verlet velocity (NVE) implemented. Not test yet.
      
      * [TODO] add friction parameter for langevin

---

## Spin

### [TODO] Documentaion of codes (hexu) [time: July, 2019] 

### [TODO] Implement dipole-dipole 

### [TODO?] Implement IMP/SIB integration mover.

### [TODO?] biqudratic terms (as computed from TB2J)

### [TODO] 

---

## Spin-Lattice coupling

### [TODO] Test two dynamics without coupling (hexu) [time: July,2019]

### [TODO] Implement primitive Oiju potential

### [TODO] Implement supercell Oiju potential

### [TODO??] Implement advanced integration algorithm for spin/lattice dynamics

---
## Lattice wannier function

---

## Electron

---
