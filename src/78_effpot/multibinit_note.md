---
Copyright (C) 2001-2020 ABINIT Group (hexu)
This file is distributed under the terms of the
GNU General Public License, see ~abinit/COPYING or http://www.gnu.org/copyleft/gpl.txt.
For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
---

# Developer's notes on Multibinit

## Manager
### The whole procedure
* Read files (outside of manager) file and input file.

* Allocation of data structure (using information from input params). 

* Read primitive potential and put them to primitive_potential_list. This includes the reading of primitive cell. So it works like this: 1. allocate a object for the primitive potential. 2. set the primitive pointer. 3. call load_from_files method of primitive potential. 4. append the potential to the primitive potential list.

* Build supercell and supercell potential. 

* Initialize movers, which use information from the input file and the supercell.

* Call the movers to run the dynamics. The movers will do: 1. call supercell potentials to calculate energies and the derivatives. 2. update system states. 3. calculate observables if needed. 4. write to hist file.

* Finalize. One thing to note is that not all the objects in the manager are allocated or initialized because they are not needed in that calculation. So for pointers, they should be initially null() and the manager should check assciated(pointer) before finalising it.  And for objects, there should be a tag to tell if it is initialized, which should be taken care of inside the objects themselves. 

## Potential
A potential usually has a primitive cell and a supercell presentation.
The primitive cell potential is the one saved in a file (xml, netcdf, etc), and of which the parameters are fitted.
The supercell potential can be used as a "calculator" to get the energies and their derivatives.
There are abstract type for both of them. To add a new potential one should add both of them. And the primitive potential should 
have a fill_supercell.

### Add has_varibles in the abstract class
    * in abstract_primitive_potential 
    * in abstract_potential
    
    TODO: since there are a few of this, and they are universal in the whole program. 
    Would it be better if we define a data structure to save it and use pointers.

### Add a section to a cell
    * The mb_cell_t type has sections (lattice, spin, lwf, ...).

### Adding primitive potential

    * make a derived type from abstract_primitive_potential_t

    * modify the function in the manager  (read_potentials). It usually has 4 steps.
    1. allocate a (abstract_primitive_t) pointer to the primitive potential needed.
    2. set the primitive_cell pointer to the potential. This should be done before the loading from file function, because some of the load function also load information for the primitive cell type.
    3. load the primitive potential from file. 
    4. append the primitive potential to the primitive potential list, which means there will be a pointer pointed to the allocated type.
    
   The finalize function will be called by the primitive potential list when itself is finalized. 

### Adding supercell potential

    * make a derived type from abstract_potential_t 

    * override the methods if necessary

    * There is NO need to do anything in the manager for building the supercell potential. ALL supercell potentials are built by the fill_supercell method of the primitve potetial list automatically.

    * Also there is no need to worry about calling the finalize function: it will be called by the potential_list finalize function. 
    

##  Movers
The movers
### Adding a new type of mover
    * make a derived type from abstract_mover_t
    * override the methods if necessary
    * add a polymorphic class pointer of the mover to the manager, so that any mover derived from this type can be used. It should be null by default, so one can know if it gets initialized.
    * add the initialization to the set_movers method in the manager. 
    * implement a run function to make the mover work. 
    * The finalization function need to be called inside the manager finalize function. One thing has to be noted: the mover might not be initialize. Therefore, before finalize it, always check using the associated(pointer). 

### Adding a new lattice mover
    * make a derived type from lattice_mover_t
    * override (at least) the initialize, finalize, and run_one_step function
    * modify the set_movers function in the manager to activate the mover (if the parameters ask for it.)
    
### Write hist file
    * A hist file 
    NOTE: remember to put the netcdf related code between if have netcdf

## How to use supercell_maker
The supercell maker use the translation symmetry to 
### translation

## 

