---
authors: XG, DCA
---

# The mrgddb utility  

The user is advised to be familiar with the main [[help:abinit]] help file
before reading the present file. It is important to read also the [[help:anaddb]]
to complement the present reading.

## 1 Introduction
  
The mrgddb code has the purpose to merge "transfer" DDBs (that were generated
from the ABINIT code) to make a complete DDB that can be exploited by the Anaddb code.

The input is very simple, and could be given directly at the screen, or more
conveniently, piped from a file. The user should provide first the name of the
new (output) DDB. He/she should then give a short description (one line) of
this new DDB that will be created. This line will be printed at the beginning of the DDB.

The user should then give the number of DDBs that will be merged, then the
whole set of filenames for the DDBs to be merged, one on each line.

## 2 What is the usefulness of a merging code ?
  
The ABINIT code in its present version is only able to produce results for one
q wavevector for each dataset. A database for more than one q point can thus
be created using MRGDDB. Also, it is useful to be able to merge different DDBs
if they are produced independently on different machines.

## 3 How does the merging code work ?
  
The DDBs are made of two parts. The first part is a list of the parameters
that were used to make the DDB, and the second part lists the 2DTEs and 3DTEs.
The merging code will check if the following variables are exactly the same in
the different input DDB: natom, [[ntypat]], [[nband]], [[acell]], [[amu]],
[[ecut]], [[ixc]], lloc, [[ngfft]], [[occ]], [[rprim]], [[typat]], [[xred]],
zion. For [[nband]] and [[occ]], the value of [[occopt]] is taken into account (see abinit help file). 
If possible, MRGDDB will produce a DDB with occopt=0.

In case of two different data sets, the code will print an error message and
stop. The code cannot merge two DDB that have been generated using two
different geometries or convergence ([[ecut]], ...) parameters. The only
exception is connected to the possibility to use Time-reversal symmetry to
decrease the number of special k points when the wavevector of the
perturbation is Gamma. In that case, the code will merge the DDBs and put the
largest set of k-points inside the new DDB. MRGDDB will copy the latest date
of the transfer or current DDB and copy it in the new DDB. It will also take
the less accurate tolwfr and copy it in the new DDB. This ends the first part
of the action of MRGDDB, namely to compare the information of the two
different DDB. 

When the checking is done, MRGDDB will check the content of the
different data blocks and constitute the new DDB by copying sequentially the
non-identical blocks and merging the identical blocks. In case two elements
are identical, MRGDDB copies the value of the transfer DDB. (This latter
property makes it easy to get rid of old, erroneous data and put new, correct
data in its place) Finally, the summary of the block content of the DDB is
provided at the end of the DDB file.
