Information on the format 1 for pseudopotentials.

The format 1 for ABINIT pseudopotentials allows to use pseudopotentials
from the set of LDA pseudotentials for the whole periodic table,
build by DC Allan and A.Khein. 
They have been generated according to the Troullier-Martins technique.
See ~abinit/doc/users/bibliography.html for the corresponding references.

The pspcod=1 psp files are formatted data files
which give potentials and projector functions on a real space
radial grid.

Firstly, the radial grid runs from index 0 to 2000 (2001 points),
with grid points given by the following equation (given as fortran):

      nmax=2000
      do j=0,nmax
       x=dble(j)/dble(nmax)
       r(j)=100.d0*(x+.01d0)**5-1.d-8
      enddo

The psp file consists of a number of header lines followed by the
data on the radial grid.  The header section is as follows:

     title  (single 80 character line)
     zatom, zion, pspdat
     pspcod, pspxc, lmax, lloc, mmax, r2well
     ...then, for l=0 to lmax, the following 2 lines:
      l,e99.0,e99.9,nproj,rcpsp
      rms,ekb1,ekb2,epsatm
     ...finally one more line:
      rchrg,fchrg,qchrg

The data may be located anywhere on the line as long as it is provided
in the order indicated (it is read with free format).
In the case of Si with lmax=2, the header may look like the following 10 lines:

    Si  Fri Oct 08 11:18:59 1993
    14.00000   4.00000    930920                zatom, zion, pspdat
      1    1    2    2      2001    .00050      pspcod,pspxc,lmax,lloc,mmax,r2well
      0  19.464  25.000    2   1.8971118        l,e99.0,e99.9,nproj,rcpsp
      .00112760   6.1457108933   4.4765165955  29.74712295   rms,ekb1,ekb2,epsatm
      1  21.459  28.812    2   1.8971118        l,e99.0,e99.9,nproj,rcpsp
      .00119946   3.2090654032   2.0935248528  19.11150542   rms,ekb1,ekb2,epsatm
      2   8.223  21.459    0   1.8971118        l,e99.0,e99.9,nproj,rcpsp
      .00098688    .0000000000    .0000000000  -3.97301006   rms,ekb1,ekb2,epsatm
      1.70000000000000     .22513330685109     .96523597101781   rchrg,fchrg,qchrg

zatom is the atomic number of the atom (14 for Si)
zion is the number of valence electrons (4 for Si)
pspdat is a code revision date (930920 for this case)
pspcod is another index describing the code (1 for this case)
pspxc is an index showing the choice of exchange-correlation (1)
lmax is the highest angular momentum for which a pseudopotential
 is defined, which is also used for the local potential (2)
lloc is the angular momentum used for the local potential (2)
mmax is the number of grid points (2001)
r2well is the prefactor of a harmonic well sometimes used to bind
 electrons which would otherwise be unbound in lda (.00050)
l is the angular momentum (0, 1, or 2 for Si for this case)
e99.0 is the planewave cutoff needed to converge to within 99.0% of the
 kinetic energy of the atom (various numbers for various l)
e99.9 is the planewave cutoff needed to converge to within 99.9% of the
 kinetic energy of the atom (various numbers for various l)
nproj is the number of projection functions used for each l (2)
rcpsp is the pseudopotential core radius
rms is a measure of pseudopotential quality reflecting the value of the
 penalty function in designing the potential
ekb1, ekb2 are the Kleinman-Bylander energies for each projection
 function for each l
epsatm is the integral Int[0 to Inf] (4*Pi*r*(r*V(r)+Zion))
rchrg is the core charge radius for additional core charge used to
 match the xc contribution to the hardness matrix
fchrg is the prefactor of the core charge expression
qchrg is the total (integrated) core charge

Following the header are, for l=0 to lmax, the pseudopotentials
in the form of a title line followed by 667 lines of data, each line
containing three numbers so that the radial grid values from 0 to 2000
are given.  The title line gives the value of l first followed by some
text.  For the case of Si, e.g., for l=0, this line is

      0 =l for Teter pseudopotential

followed by the 667 lines, 3 numbers each, giving the l=0
potential on the radial grid described above.

Following the pseudopotentials are the first projection functions,
again given for each l with a title line followed by 667 data lines.

Following the first projection functions for each l are the second
projection functions, if any (determined by nproj).

Following the second projection functions, if any, are several lines
of additional data which is not read by plane_wave but is provided
to describe more about the details of the construction of the
pseudopotential.  Omission of these lines does not affect the running
of plane_wave (at this time).

If you do not wish to use core charges, simply set fchrg to 0 and
use rchrg=1, qchrg=0.

If you wish to make a local potential, use lmax=lloc=0, nproj=0, and
you need not provide any projection function(s) (at this time).
The values of rms, ekb1, ekb2, epsatm, e99.0, e99.9 are used only
for information (at this time) so may be set to 0 when creating
a pseudopotential file.  rcpsp is still used as the definition of the
pseudopotential core radius so it must be provided.


To best understand this data structure, it is recommended to study
several pseudopotential files and compare their contents with the
description given here.


Inside ABINIT, a pseudopotential with format 1 will be treated by
the routine psp1in.f, that calls psp1lo.f (local part),
psp1nl.f (non-local part), and psp1cc.f (XC core correction).

As a matter of numerical accuracy, note that the integral
of (V(r)+Zion/r) r^2 in psp1lo.f is performed from 0 to the highest
allowed radius (usually about 100 a.u.), without cut-off.
V(r)+Zion/r should tend rapidly to zero
for large radii (beyond 5 a.u.), but this correct behaviour will
not be enforced by the routine. If the tail of V(r) is inaccurate
(i.e. if the pseudopotential is in single precision), there
will be large inaccuracies in the integral, because of the r^2 factor.
