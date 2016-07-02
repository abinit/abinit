
Known problems

This file contains the list of known problems
at the time a version of ABINIT is released.
Note that bugs discovered later than the
release are NOT described here...
(In particular, at the time of release, a particular
 version is not tested already on ALL platforms,
 only on a subset of 5, usually).

For version 4.7.x, list closed on 20060131
after installation and tests for: 	
seq and paral        PC/Linux with PGI compiler (~ABINIT/Mach*/*P6/makefile_macros.PGI4.0.2_sleepy)
no                   PC/Linux with PGI compiler (~ABINIT/Mach*/*P6/makefile _macros.PGIstatic_dummy)
seq and paral        PC/Linux with IFC compiler (~ABINIT/Mach*/*P6/makefile_macros.IFC8.1-32_sleepy)
no                   PC/Linux with IFC compiler (~ABINIT/Mach*/*P6/makefile_macros.IFCstatic_dummy)
no                   PC/Linux with g95 compiler (~ABINIT/Mach*/*P6/makefile_macros.g95_sleepy)
no                   PC/Linux with pathscale compiler (~ABINIT/Mach*/*P6/makefile_macros.Pathscale2.1-64_sleepy)
no                   Itanium2 machine with IFC compiler (~ABINIT/Mach*/It*/makefile_macros_chpit_ifc8.1-ia64)
seq and paral        Compaq/Dec EV67 OSF (~ABINIT/Mach*/DEC/makefile_macros.decci)
seq and paral*       Compaq/Dec EV56 Linux (~ABINIT/Mach*/DEC/makefile_macros.DecLinux.tux)
seq                  IBM P44 AIX (~ABINIT/Mach*/DEC/makefile_macros.dirac_nc)
no                   SGI Octane IRIX (~ABINIT/Mach*/DEC/makefile_macros.spinoza)

Compilation of Netcdf was not tested
Compilation of LibXC  was not tested

Warning : paral with Compaq/Dec EV56 did not work, see P47.5

All the problems are listed cumulatively.
Those specific to version 4.6 are found at the end
of this file

===============================================

The type of the problem is identified as follows :

NUMERICAL ERROR : the physical result is incorrect,
 but the job does not crash. Moreover,
 the use of the input variables is legitimate.
 This is a very dangerous type of problem.
 There are two subclasses :
 - INTERNALLY INCOHERENT
  (the bug can be seen because two sets of
   input variables should give the same results
   and do not)
 - INTRINSIC
  (the bug does not cause incoherency, it can
   only be seen when comparing results with
   those of another code or with experimental results)
 The latter one is the most dangerous, of course ...

CRASH : the code does not produce the result, because
 it stops unduly.

SLOWING DOWN : the physical result is correct,
 but the code does not go as fast as it should.
 For example, the reading of a restart file might
 not be correct, or the SCF convergence
 might not be as fast as it should.

OTHER : some information (not numerical results)
might be erroneous ...

In each class of problem, one can identify
MACHINE-DEPENDENT problems (on some platform, this
problem does not occur, not on others) as well as
GENERAL problems.

Some of the problems listed below
have been reported by user, but have not yet been
examined for confirmation, in which case they are
labelled "TO BE CONFIRMED".

=================================================

P0. (Many reports)
SLOWING DOWN
The iterative procedures (SCF and determination
of geometry) still present instabilities,
in particular for metallic slabs

P2. (XG, 001223 + later)
OTHER  (TO BE DONE)
GENERAL
The memory needs are not yet evaluated correctly
for the spin-orbit case and non-collinear magnetism
case.

P8. (GMR, November 2000)
SLOWING DOWN
GENERAL
restartxf for optcell/=0 might have a bug,
despite Test_v1#80
Identified in March 2001 by GMR : only happens when
restarting from a different, non-zero, optcell.
Should be corrected, but documented, now.

P9. (Massimilio Stengel, October 2000)
CRASH
MACHINE-DEPENDENT : Compaq/DEC
Problem with non-linear XC correction, it seems

P15. (XG, 010103)
OTHER
GENERAL
Should enforce coding rule related to length of file

P117. (Laurent Pizzagalli, 010622)
CRASH
Observed on PC AMD, atomic relaxation (ionmov==7),
no cell optimisation, stops in xfpack complaining
about ndim. See mail 010622.

P118. (XG, 010710)
OTHER
GENERAL (d/dk by finite differences)
Test_v3 #5 . Sensitivity to phase choice.

P121. (GMR, 010801)
NUMERICAL RESULT
MACHINE_DEPENDENT (Compaq ES40, Turing)
The use of fourdp , with fftalg 200 or 1xx
gives complex conjugate results. This was seen
in the process of merging the GW code with ABINIT.
To be examined in detail.

P204. (XG, 011020)
OTHER (TO BE DONE)
MACHINE_DEPENDENT
Test_v3#14 : the output file is machine-dependent,
because the spin-phase of spin-degenerate wavefunctions
has not been fixed...
Or, is there something more fundamental ??

P205. (XG, 011020)
SLOWING DOWN
MACHINE_DEPENDENT ??
On a Intel/Linux PIII 900 MHz, with PGI
compiler, slowing down in getghc,
due to the use of datastructures...
For example, 18% of the total time
in Test_v3/t10.out ... Argh ...
even more in Test_v2/t30.out ...

P220. (XG, 020114)
OTHER
GENERAL
Refs files in Test_physics should be updated.

P230. (MTijssens, 020306)
TIMING PROBLEM
MACHINE-DEPENDENT (P6, with IFC compiler v5.0)
The wall clock timing of one routine can be completely
wrong. Seems to come from the date_and_time
intrinsics, called by timein.f .
Previously, it was even causing
crash, in the headwr.f routine.

P250. (LSi, 020521)
CRASH
MACHINE-DEPENDENT (DOS/Windows)
Access violation Test_v2#7,26,30,90,92,96
and Tutorial#55,56

P302. (XG, 020305; DDKlug, 020313)
NUMERICAL ERROR
MACHINE_DEPENDENT (IBM 44P; SGI)
RF GGA : Test_v3#8,16,18,60,61,62,66
Problem with Response Function GGA tests.
The accuracy of the RF GGA Test_v3#16
is very bad on these machines : the acoustic
modes have large frequencies...

P331. (XG, 020528)
SLOWING DOWN
MACHINE_DEPENDENT (MACOSX)
Test_v2 #53 and 54
The computation of the dielectric matrix
is numerically incorrect. Fortunately, this
dielectric matrix is only used to precondition the
SCF cycle.

P332. (XG,020528)
OTHER
MACHINE_DEPENDENT (MACOSX)
Test_v3 #31 and 32
The Fermi energy is not printed out correctly.

P349. (XG021016)
OTHER
LIKELY GENERAL
Test_v3 #77. This test is not yet portable : accuracy
of fldiff was already tuned, but still problems
of portability.

P350. (XG030225)
OTHER
GENERAL
psp*cc.f the fifth derivative is not computed.

P402. (XG021227)
NUMERICAL ERROR
LIKELY GENERAL
Lattice Wannier Function only
Test_v3#88 and 89 are still undergoing changes

P404. (XG030103)
NUMERICAL ERROR
LIKELY GENERAL
fftalg 401 works, but gives a slightly
different answer than usual for Tfast#10.
To test this, replace the default value in indefo.f

P405. (XG030106)
CONVERGENCE PROBLEM
LIKELY GENERAL
Test_v3#15 : with very well non-SCF converged wavefunction (nnsclo=2),
the total energy at intermediate steps is lower than the
final energy !

P406. (XG030106)
CONVERGENCE PROBLEM
LIKELY GENERAL
Test_v3#85 : set iscf to 5 in the RF, one sees a very bad convergence...

P407. (XG030108)
NUMERICAL ACCURACY
INTEL with IFC compiler, EV67 under OSF, EV56 under Linux
Test_v3#14,77,85 with diffs compared to refs, that are a bit too large.

P411. (XG030515)
PORTABILITY PROBLEM
INTEL with IFC compiler, EV67 under OSF, EV56 under Linux
Test_v4#40 differences compared to refs are a bit large
Test_v4#42 sign of some eigenvectors might not be portable

P416. (XG030522)
CRASH
Seen on Intel/PC with PGI compiler, in parallel
Test_v1#75,80  do not work
Problem to read XFhistory part of a wavefunction file.

P41.9. (XG030719)
PORTABILITY PROBLEM
GENERAL
Test_v4#37 is not portable : the spin decomposition of the
spinor wavefunctions is not gauge invariant. This
test has been disabled.

P41.12. (XG031127)
CRASH
HP 8000
Test_v3#88 and 89 : LWF CRASH

P42.7. (DH030929)
OTHER
Should examine the d2sym3.f routine, giving extra
zeros in the output file. See mail
from DHamman, on 2003-9-29.

P43.3. (XG040218)
PORTABILITY of AUTOMATIC TESTS
On P6 with INTEL compiler, DEC EV67, DEC EV56
Test_v4#54,67,68,71,75,77,78 and
Tutorial#61-69 :
Accuracy criterion to be adjusted,
or portability improved.

P43.12. (DRH040227)
NUMERICAL ERRORS
On a PC with IFC compiler version 5.0.1
Test_v4#77
Erroneous polarization

P43.14. (XG040309)
PORTABILITY
On SGI "Spinoza"
Test_v4#54 , anaddb analysis of eigenvectors
gives a different sign.

P43.17. (XG040123)
CRASH
On Itanium2 (Chpit)
Tutorial#69
Error message "Unaligned access to ..."

P44.2. (DHamann20040527)
NUMERICAL INACCURACY
For the treatment of pseudopotentials with pspcod=6 (FHI),
the non-linear XC core charge is not treated consistently
in psp6cc.f and in the rest of the code.
See Don's mail to the developers of 20040527
See also input variable optnlxccc .

P44.12 (MM20040901)
PORTABILITY / NUMERICAL ERROR
IFC v7.0
From Mikami-san
Test_v4#62 : completely erroneous ?!

P44.18 (XG041102)
OTHER : PARALLEL CASE
GW is not parallelized.
However, all the automatic tests run "in parallel" now
(so, in a few cases without speed-up).

P44.19 (XG041216 - updated XG050614 - updated XG050723)
Problem with the script abirules
"make abirules"
causes unwanted changes in qmcout.F90 .

P45.5. (XG050220)
OTHER : PARALLEL TESTING
Test_v4#90 (CASINO output) crash in parallel (test_paral on grumpy).

P45.6. (XG050225)
NUMERICAL PROBLEM (with ionmov=12)
Test_v4#97 has one strange output

P45.8. (XG050227)
SCRIPT PROBLEM
Too much errors observed when using the script parents : should change the rules
or the script

P45.14. (XG050619)
Libxc from Miguel is still to be tested.

P45.16. (N. Choudhury 050107)
In case of a very high symmetry unit cell, with atoms
placed such that the space group is of low symmetry, there
is a problem to recognize the space group.
This is illustrated by Test_v4#98.

P45.17 (B. Andriyevsky 050405)
Problem of normalisation of the epsilon'' from
the "conducti" code ?!
A factor of pi (3.1416) might be missing.
See his mail of 5 April 2005 to the forum.
See also other mails exchanged on the forum.

P45.18 (D Hamann 050510)
NUMERICAL ERROR
With finite-difference ddk calculation
See his mail of 10 May, 2005.
(Excerpt : ...

The problem is that 4.4.4 and 4.5.2 WILL run drf.in and produce incorrect
results, all related to the dataset 3 finite-differend ddk calculation.
The results are only somewhat incorrect, not obvious garbage. Here is a
small sample, Born charges from "2nd-order matrix" in the .out files:

Correct (from erf.out):
   1    1   1    4   -4.1367712586    0.0000000000
   1    1   2    4   -0.0863319602    0.0000000000
   1    1   3    4    0.0706059161    0.0000000000

Incorrect (from drf.out):
   1    1   1    4   -3.4782061716    0.0000000000
   1    1   2    4   -0.0832639192    0.0000000000
   1    1   3    4    0.0560776887    0.0000000000


Putting in some istwfk printing, I found that setting berryopt/=0 now forces
istwfk=1 for dataset 3.  The problem apparently is that dataset 2 has some
istwfk>1, and apparently dataset 3 is getting incorrect WFK's from dataset 2
and using them.  Why it doesn't re-converge them first isn't clear.

Setting istwfk=1 in dataset 2 (erf.in) fixes the problem and gets correct
results.  (My hack also had the effect of setting istwfk=1 in dataset 2 as
well as 3.)

Curiously, removing the dataset 2 calculation altogether and feeding the
kptopt1=1 WFK's from dataset 1 directly to the finite-difference ddk
calculation (nrf.in) also fixes the problem, despite the fact that dataset 1
has some istwfk>1.

Clearly something bad is going on in WFK readin, but only in certain cases.
I've never looked at this part of the code, and I'm afraid I can't divert
any more time now towards tracking this down.



P46.1. (XG050701)
PORTABILITY PROBLEM
Test tutorial Finite electric field
alas_relaxdion.in (the last one)
does not work on my PC : it stops before the end.

P46.4. (XG050714)
NUMERICAL INACCURACY / PORTABILITY
Test_v2#13,15,22
Relaxed ion dielectric tensor varies quite a lot between PC/IFC8.1 and PC/PGI4.0.2

P46.6. (XG050716)
MISC.
Variables nqptdm and qptdm should be substituted with nqpt and qpt.

P46.7. (XG050716)
PORTABILITY
DEC EV67 (Deccint)
afterscfloop.F90 only compiles with -O1 .

P46.8. (XG050716 - modified XG050723)
PORTABILITY
Itanium 2 (Chpit)
Test_v4#30, initialisation differs widely from other machines ?!

P46.13. (XG050724)
PORTABILITY
Tutorial#elast_6, on Tux and Decci
The number of SCF steps differs from the others....
Should try to set   nstep33  20

P46.15. (XG051002)
(PORTABILITY)
Tutorial#optic_1 and optic_3, on PC/IFC, SGI (Spinoza)
and IBM 44P (Dirac).
Some of the printed residuals are 10-310 .
Other results seem fine.

P46.16. (XG051002)
(PORTABILITY)
Tutorial#nlo_4, on DEC EV67 (Deccint) and SGI (Spinoza)
and P6 under g95
DDB contains numbers like 0.12370492835605D-55
instead of zeroes.

P46.19. (XG051003)
(PORTABILITY)
Test_v4#66 on IBM 44P (Dirac)
Strain response not accurate in the DDB

P46.21. (XG051017)
Memory leaks
Using the g95 compiler with debugging option,
makefile_macros.g95_sleepy_debug  (-g)
Identified remaining memory leaks in
Test_v2#01 : hdr not cleaned (init in hdr_init routine)
Test_v3#30 : hdr not cleaned (init in hdr_io routine)
Test_v4#50 : in sortph.F90

P47.1. (XG060130)
In the density mixing scheme (see input variable iscf),
the total energy as printed during the SCF procedure,
is not variational. Hence, it does not decrease monotonically.
Note however, that it converges to the correct value,
and that all the other quantities are not affected by this
lack of variational behaviour of the total energy.

P47.2. (XG060131)
Test_paral Cases K and L do not work yet.

P47.3. (XG060131)
(PORTABILITY)
IBM 44P
Test_v4#01 and 03 crash
Use Ylm for non-local norm-conserving psps

P47.4. (XG060131)
(PORTABILITY)
IBM 44P
Test_v4#55 and 78 crash
Tutorial#ffield_6 crash
Finite electric field

P47.5. (XG060131)
(PORTABILITY)
Test_paral#A2 and following
Alpha/compaq EV56 under Linux
Parallelism does not work ?!
Hence, Test_paral/Refs/t*10.out have not been updated
