---
authors: XG
---

# File format for ABINIT : Norm-Conserving Pseudopotentials 

The pseudopotential files for ABINIT consist of a number of header lines, followed by the data on a radial grid.  
The three first lines of the header have the same format and meaning for all norm-conserving pseudopotential files that can be read by ABINIT. They are :

     title  (single 80 character line)
     zatom, zion, pspdat
     pspcod, pspxc, lmax, lloc, mmax, r2well
     
The data may be located anywhere on the line as long as it is provided
in the order indicated (it is read with free format).
In the case of Si with lmax=2, the header may look like the following three lines:

     Si  Fri Oct 08 11:18:59 1993
     14.00000   4.00000    930920                zatom, zion, pspdat
     1    1    2    2      2001    .00050        pspcod,pspxc,lmax,lloc,mmax,r2well

The first line is free, and might contain information useful for the pseudopotential generator only.
In the present case, it means that this silicon pseudopotential was created in 1993.
The next two lines are important for ABINIT, with the following meaning:

     zatom   : atomic number of the atom (14 for Si)
     zion    : number of valence electrons (4 for Si)
     pspdat  : code revision date (930920 for this case)
     pspcod  : identifier for the pseudopotential format (crucial for further reading, see below !)
     pspxc   : the choice of exchange-correlation, coherent with the ABINIT nomenclature
     lmax    : highest angular momentum for which a pseudopotential
      is defined, which is also used for the local potential 
     lloc    : angular momentum used for the local potential
     mmax    : number of grid points 
     r2well  : prefactor of a harmonic well sometimes used to bind
      electrons which would otherwise be unbound in lda (usually 0.000)

At this point, differences might appear in the different file formats, that are described in the following pages:

  * pspcod=1, [psp1 format](./psp1_info.md)
  * pspcod=3, [psp3 format](./psp3_info.md)
  * pspcod=4 or 5, [psp45 format](./psp45_info.md)
  * pspcod=5, [psp5 spinorbit format](./psp5spinorbit_info.md)
  * pspcod=6, [psp6 format](./psp6_info.md)
  * pspcod=8, [psp8 format](./psp8_info.md)

By far, the most used format for NC pseudopotentials for ABINIT is psp8.
