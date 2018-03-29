Information on the format 3 for pseudopotentials.

(Note: the implementation of format 3 was done by Fr. Detraux).

The format 3 for ABINIT pseudopotentials allows to use pseudopotentials
from the table I in the Phys. Rev. B 58, 3641 (1998) paper
by C. Hartwigsen, S. Goedecker and J. Hutter (HGH). This paper
presents LDA pseudopotentials for all elements from H to Rn. Some of them
are even presented twice, because of the possibility to include semi-core
states. Their accuracy has been demonstrated in the HGH paper, but note
that the energy cut-off needed to get this high accuracy might be
larger than the one usually needed for Troullier-Martins or other
pseudopotentials. So, convergence studies are very important !

A few lines must be added to the data mentioned in that table,
and are described in the present file. ABINITv1.5 is able to
read format 3 pseudopotential files, as well as later versions of ABINIT.

We will suppose that the user has the HGH table.
For the description that follows, we will focus on the Sn pseudopotential.
In the table (p. 3651), one finds:

    Sn  4  0.605000  4.610912
           0.663544  1.648791 -.141974  -.576546
           0.745865  0.769355 -.445070
                     0.103931 0.005057
           0.944459  0.225115
                     0.007066

These data must be put in the following format to be read by ABINIT
(this file is also found in the `~abinit/tests/Psps_for_tests directory`,
with the name `50sn.psphgh`):

    Hartwigsen-Goedecker-Hutter psp for Tin,  from PRB58, 3641 (1998) paper
    50  4   980509                     zatom,zion,pspdat
    3   1   2    0    2001    0        pspcod,pspxc,lmax,lloc,mmax,r2well
    0.605000  4.610912 0         0         0  rloc, c1, c2, c3, c4
    0.663544  1.648791 -.141974  -.576546     rs, h11s, h22s, h33s
    0.745865  0.769355 -.445070  0            rp, h11p, h22p, h33p
              0.103931 0.005057  0                k11p, k22p, k33p
    0.944459  0.225115 0         0            rd, h11d, h23d, h33d
              0.007066 0         0                k11d, k22d, k33d


First, the following three lines have been added at the beginning
of the data from the table:

    Hartwigsen-Goedecker-Hutter psp for Tin,  from PRB58, 3641 (1998) paper
    50  4   980509                     zatom,zion,pspdat
    3   1   2    0    2001    0        pspcod,pspxc,lmax,lloc,mmax,r2well

Similar lines must be added to the data for other elements, in order
to make them readable by ABINIT.

Line 1 is simply a header, that might include any information, and that will
be printed without modification in the output of ABINIT.

Line 2 describes:

- the atomic number (zatom);
- the ionic charge (zion, number of valence electrons);
- the date of pseudopotential generation.

The two first information are crucial, the third one is not
really important. The atomic number is 1 for Hydrogen, 8 for Oxygen, and so on.
The ionic charge should be the same number as the second number mentioned
in the data table (described as 'Zion' in section VI of the paper)

Line 3 describes:

- the format of the pseudopotential (pspcod ; must be 3 for this format);
- the XC functional used to generate the pseudopotential (pspxc; here, must be 1)
- the maximal angular momentum of the wavefunctions described
   in the pseudopotential file
   (lmax=0 if purely local pseudopotential: H, He, Li^sc, Be^sc ;
    lmax=1 for s and p projectors : Li, B ... Ar, K^sc, Au ;
    lmax=2 for s, p and d projectors : K, Ca, ... Cs, Ba, Hf^sc ... Pt^sc,
     Au^sc ... Rn ;
    lmax=3 for s, p, d and f projectors : Cs^sc, Ba^sc, La^sc ... Lu^sc )
- lloc has no meaning, and is set to 0
- mmax has no meaning, and is set to 2001
- the last number can be set to 0 .

The lines that follow these three lines
are generated from the data in the table, though the name of the
element and Zion are not reproduced. Note also that
each line has to be completed with zeroes, to give the format presented
at the beginning of the section VI of the HGH paper :
five columns for the first of these lines,
then 4 columns for each projector, and 3 columns for the spin-orbit splitting.
For readability, the meaning of each number has also been added
in our example pseudopotential file
(for example, 'rloc, c1, c2, c3, c4' ) but these are NOT read by ABINIT.
Thus, unlike the zeroes, it is not important to add them to the
data from the table.

Inside ABINIT, the pseudopotential with format 1 will be treated by
the routine psp2in.f, that calls psp2lo.f (local part) and
psp2nl.f (non-local part).

As a matter of numerical accuracy, note that the integral
of (V(r)+Zion/r) r^2 in psp2lo.f is performed analytically without cut-off.
