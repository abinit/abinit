---
This file list the things to do for multibinit
---
# Todo list for Multibinit

## Temporary: for short term todos. Will be moved to other sections or removed once it is done.
### [TODO] add energy table in manager.
### [TODO] Check acoustic sum rule of Oiju and Tijuv term. 
   \sum_u Oiju = 0, \sum_u Tijuv=0, \sum_v Tijuv=0, \sum_uv Tijuv=0
### [TODO] check if Oiju parameter fitting has correct sign for ij permutation.
### [TODO] implement sorting so that  (@Nicole, Can you confirm which ones?)  
 - Oiju could be grouped by ju?
 - Tijuv could be grouped by uv? ju? 


## Preparing for Abinit9

### [TODO] Update the automatic tests (hexu) [time: November, 2019]

### [TODO] Add more autotests. (Nicole) [time: November, 2019??]

      - list of SLC tests?
[TODO] Autotests with Multibinit --F03.

### [Ongoing] Finalize the format of potential netcdf file and spin hist file. (hexu & Nicole) [time: December]

      - Oiju format updated.
### [TODO] documentation for the file formats. (Nicole & hexu)

#### [TODO] use reduced coordinates instead of cartersian? Do we have the final decision? (hexu)  [time: November]
### [Ongoing] Full workflow of fitting parameters. (hexu) [time: November]
### [TODO] Update the tutorials and other documentations. (Nicole & hexu) [time: November]
### [TODO] Tutorial for spin-lattice coupling???


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


#### [DONE] Make a developer's guide on how to add new things. (hexu) [time: July, 2019]

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

### [TODO] Implement primitive Oiju potential (Nicole) 
 - Other terms are also implemented. Test to be done.

### [TODO] Implement supercell Oiju potential (Nicole)
 - Other terms are also implemented. Test to be done.

### [TODO??] Implement advanced integration algorithm for spin/lattice dynamics

### [TODO] Optimize Tijuv terms so the time is acceptable. (Nicole & hexu) [time: Now]
    - The firt time consuming part is the sorting of the ijuv indices. Merge sort (O(NlogN)) is now (07/11/2019) instead of insertion sort (O(N^2)). More test needed to see if calculating the force and bfield is fast enough, especially on larger supercell. 


### [TODO] Implement full procedure of coupled SLC dynamics. (Nicole)

---
## Lattice wannier function

---

## Electron

---


## Long term issues
### [ISSUE] damping parameter for spin dynamics. (1) how to get the damping parameter for spin dynamics. (2) Counting spin-phonon coupling


