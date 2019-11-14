# Multibinit Enhancement Proposals

## MEP0: keep a list of proposed changes as MEP.
If some changes to be basic structure of multibinit, put them here.

## MEP1: add label property to potentials (for debugging)
For each potential, it has a default label, which can be modified.
After it is finalized, it becomes "Destroyed potential". 

## MEP2: add get_delta_energy method to each potential. (For Monte Carlo)
- get_delta_energy( displacement, strain, spin, lwf, ispin, dS, ilatt, dtau, ilwf, dlwf, energy)
all params are optional, but should be given in pairs (ispin, dS).
This is limited to one move per step (simple metropolis-Hastings)
possible extension: arrays of changes for other MC methods.

## MEP3: Do not use mover list
- We need very fine control of each mover. Putting them into a list and loop over them is not useful. 
And FORTRAN polymorphous  pointer is quite verbose to use. 

## MEP4: save global state in a type (extends Multibinit_dtset_type).
- The potentials and movers need to share some states, mostly the input parameters and things derived from them, e.g. current temperature from temperature range.  . We can extend multibinit_dtset_type and add some attributes. And a pointer to it can be added to the movers and potentials if necessary. 

## MEP5: specify which data structures are mpi aware.
- If yes, it has to be outside if(iam_master) 
     else it has to be inside if(iam_master)
     
## MEP6: Unit test mode
 Unit tests are useful but not so easy in Abinit/Multibinit. We can make a module m_multibinit unit test and call unit test functions from there.
 Unit test functions can be implemented in the modules they are for.
 By running "multibinit --unittest", the unit tests will run instead of the main Multibinit.
 
## MEP7: use/not use sub type in primitive/supercell?
- The lattice part use crystal_t from abinit. we can either extend this type, or make one type for each nature. Which one is better?

     


- 
