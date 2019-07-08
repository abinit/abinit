---
This file list the things to do for multibinit
---
# Todo list for Multibinit

## Genenal structure

### [TODO] Clean lattice part in 98_main/multibinit.F90 to manager or lattice code

State: 


### [TODO] Separate the 78_effpot directory into layers 

State: 

### [DONE] Make a example harmonic lattice potential. (hexu) [time: July, 2019]


#### [DONE] implement primitive potentail (read from netcdf file). Tested.
#### [DONE] implement supercell potential, fill_supercell
#### [DONE] add the data structure to the manager.
###  [DONE] move random number generator to manager and use pointer in each mover.

### [ISSUE] Where to store the total energy.
	When there are multiple potential, each mover get part of the energies but not
       all of them.
	e.g. For spin+lattice. Spin mover get E(spin) + E(spin-lattice-coupling).
            Lattice_mover get E(lattice) + E(spin-lattice-couping) 
        Another problem is from the kinetic energy, which is not calculated by the potential
	but the lattic mover.
        Also note that the energy is related to specific heat, susceptibility, etc. 

---

## Documenation
### [DONE] Document general data structure  (hexu) [time: July, 2019]

State: ?

#### [TODO] Make a developer's guide on how to add new things. (hexu) [time: July, 2019]

State: First draft in multibinit_note.md
 
### [TODO?] Make a full technical note (including all the methods used) (everyone) 

State: One file exist for spin/spin-lattice coupling.


---

## Lattice

### [TODO] Adapt to general data structure

State: 

### [TODO] Mover

#### [TODO] set intial state:  (hexu) [time: July 2019]
       * [DONE] implement get_temperature and get_kinetic_energy
       * [TODO] More options other than Maxwell-Boltzmann distribution
       * [DONE] make center of mass fixed.
       * [TODO] make sure there is no global rotation. (useful for non-periodic structure only.)

#### [TODO] implement movers. (hexu) [time: July 2019]
      * [DONE] implement Velocity Verlet mover
      * [DONE] implement Langevin mover
      * [DONE] implement Berendsen NVT
      * [Pending] implement Berendsen NPT  (Pending due to no testcase available)
      * [DONE] add friction parameter for langevin, taut for Berendsen, (taup) for Berendsen NPT
      * [DONE] Add documentation to input variables


---

## Spin

### [TODO] Documentaion of codes (hexu) [time: July, 2019] 

### [TODO] Implement dipole-dipole 

### [TODO?] Implement IMP/SIB integration mover.

### [TODO?] biqudratic terms (as computed from TB2J)


---

## Spin-Lattice coupling

### [DONE] Test two dynamics without coupling (hexu) [time: July,2019]

### [TODO] Implement primitive Oiju potential

### [TODO] Implement supercell Oiju potential

### [TODO??] Implement advanced integration algorithm for spin/lattice dynamics

---
## Lattice wannier function

---

## Electron

---
