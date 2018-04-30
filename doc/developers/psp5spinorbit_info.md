This file contains information about the use of spin-orbit
for format pspcod=5 of pseudopotentials.

The first line containing Haman grid parameters must be completed by
an information about the spin-orbit format of the pseudopotential.
So, one replaces:

    2.508991628593723E-4  0.0125           r1,al

by

    2.508991628593723E-4  0.0125  2        r1,al,pspso

pspso is 1 for non spin-orbit (optional),
         2 for spin-orbit


The next lines describing the number of projectors and some of their
characteristics must be duplicated, as soon as l/=0.
The number of projectors must be 2.

So, one replaces

    0    0    0    1   2.76        l,e99.0,e99.9,nproj,rcpsp
    .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
    1    0    0    1   3.91        l,e99.0,e99.9,nproj,rcpsp
    .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
    2    0    0    1   1.57        l,e99.0,e99.9,nproj,rcpsp
    .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm

by

    0    0    0    1   2.76        l,e99.0,e99.9,nproj,rcpsp
    .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
    1    0    0    2   3.91        l,e99.0,e99.9,nproj,rcpsp
    .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
    1    0    0    2   3.91        l,e99.0,e99.9,nproj,rcpsp
    .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
    2    0    0    2   1.57        l,e99.0,e99.9,nproj,rcpsp
    .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
    2    0    0    2   1.57        l,e99.0,e99.9,nproj,rcpsp
    .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm


The pseudopotentials and pseudowavefunctions are then entered,
in the order  s, p1/2, p3/2, d3/2, d5/2, f5/2, f7/2.


The local potential is specified by lloc with the following values:

    lloc=0  refers to the s potential
    lloc=1  refers to the p3/2 potential
    lloc=2  refers to the d5/2 potential
    lloc=3  refers to the f7/2 potential
    lloc=-1 refers to the p1/2 potential
    lloc=-2 refers to the d3/2 potential
    lloc=-3 refers to the f5/2 potential
