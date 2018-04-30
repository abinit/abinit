## This file contains information about formats of pspcod=4 and pspcod=5 pseudopotentials.

Pspcod=4 corresponds to the case of Teter pseudopotentials
generated in Louvain-La-Neuve using the code ATOM. Pspcod=5 corresponds
to "Phoney" pseudopotentials built on a Hamman grid. At this stage (22 July
1998) it is possible to treat these pseudopotentials (in the version 1.5 and
later).  However, one has to be careful about the formats of these
pseudopotentials. Indeed, to be read by the ABINIT code, they have to be slightly modified. 
Let's take examples to explain these modifications:

# pspcod = 4 case (Teter pseudopotentials). 
 
The example corresponds to Lead.

OLD HEADER

    (Xe+4f14)+6s1.8 5d10 6p0.2 5f0.05;rcs=rcd=2.0(exnc11),rcp=2.0(26),rcf=1.3(11) no chem-hard; ecut 19/25
     82.00000  14.00000  0.000000000000000000E+00      z,zion,etot(wrong)
            4         3         2      2001            iexc,ipsp,lmax-1,ngrid
     .000   .000000000000000000E+00   .000000000000000000E+00 rchrg,fchrg,qchrg
    2.00426660272461010 2.00426660272461010 2.00426660272461010 1.29915156996312131
    2 2 2 0

NEW HEADER (READABLE BY ABINIT)

    (Xe+4f14)+6s1.8 5d10 6p0.2 5f0.05;rcs=rcd=2.0(exnc11),rcp=2.0(26),rcf=1.3(11) no chem-hard; ecut 19/25
      82.00000  14.00000  960808                 zatom,zion,pspdat
      4         3         3   3      2001    0  pspcod,pspxc,lmax,lloc,mmax,r2well
      0  0  0   2    2.00426660272461010          l,e99.0,e99.9,nproj,rcpsp
      .000 .000 .000 .000                         rms,ekb1,ekb2,epsatm
      1  0  0   2    2.00426660272461010          l,e99.0,e99.9,nproj,rcpsp
      .000 .000 .000 .000                         rms,ekb1,ekb2,epsatm
      2  0  0   2    2.00426660272461010          l,e99.0,e99.9,nproj,rcpsp
      .000 .000 .000 .000                         rms,ekb1,ekb2,epsatm
      3  0  0   0    1.29915156996312131          l,e99.0,e99.9,nproj,rcpsp
      .000 .000 .000 .000                         rms,ekb1,ekb2,epsatm
      .000   .000   .000                          rchrg,fchrg,qchrg


This header corresponds exactly to these of the TM pseudoptentials (pspcod=1).
In addition, in the new form of the Teter pseudopotentials, one indicates the
number of angular momentum before each set of data corresponding to the
considered angular momentum. This allows a clearer pseudopotential file.

# pspcod = 5 case ("Phoney" pseudopotentials).

The example corresponds to Oxygen.

OLD HEADER

    Compromise psp for oxygen with rc=1.5 ec=25 double reference
       8.00000000000000000       6.00000000000000000      -15.6648128583845843
               4           2           0         600
      0.999999999999999955E-06  0.307523885541775704E-01
       1.49157651759552068       1.49157651759552068
               2           0

NEW HEADER  (READABLE BY ABINIT)

    Compromise psp for oxygen with rc=1.5 ec=25 double reference
       8.000       6.000      980710     zatom,zion,pspdat
       5   3   1   1   600     0         pspcod,pspxc,lmax,lloc,mmax,r2well
       0.999999999999999955E-06  0.307523885541775704E-01  r1,al
       0  0  0   2     1.49157651759552068               l,e99.0,e99.9,nproj,rcpsp
       0  0  0   0                                       rms,ekb1,ekb2,epsatm
       1  0  0   0      1.49157651759552068              l,e99.0,e99.9,nproj,rcpsp
       0  0  0   0                                       rms,ekb1,ekb2,epsatm
       0  0  0                                           rchrg,fchrg,qchrg

Note that this header is almost the same as for the pspcod=1 and pspcod=4 case.
The only difference lies in the presence of "r1" and "al" parameters (defining
the Hamman grid). Here also, we indicate in addition the number of each angular
momentum before the set of data corresponding to the considered angular
momentum.

Inside ABINIT, a pseudopotential with format 4 will be treated by
the routine psp1in.f, that calls psp1lo.f (local part),
psp1nl.f (non-local part), and psp1cc (core correction).

As a matter of numerical accuracy, note that the integral
of (V(r)+Zion/r) r^2 in psp1lo.f is performed from 0 to the highest
allowed radius (usually about 100 a.u.), without cut-off.
V(r)+Zion/r should tend rapidly to zero
for large radii (beyond 5 a.u.), but this correct behaviour will
not be enforced by psp1lo . If the tail of V(r) is inaccurate
(i.e. if the pseudopotential is in single precision), there
will be large inaccuracies in the integral, because of the r^2
factor.

By contrast, a pseudopotential with format 5 will be treated by
the routine psp5in.f, that calls psp5lo.f (local part),
psp5nl.f (non-local part), and psp5cc (core correction).

The integral of (V(r)+Zion/r) r^2 in psp5lo.f is performed
from 0 to the highest allowed radius (usually about 100 a.u.),
except that beyond 20 a.u. no value of abs(V(r)) larger
than 2.0d-8 is tolerated. This will allow to cut off spurious
behaviour of pseudopotential whose data file is written in
single precision.
