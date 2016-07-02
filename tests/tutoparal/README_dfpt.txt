
To run the dfpt test case.

-----------------------------------------

First test : FCC aluminum, computation of the phonon frequencies at q=0.25 -0.125 0.125 .

This test, with only one atom, is very fast, but does not scale very well,
as the sequential part will quickly dominate. It is actually about 1.5% of the
code on one processor. So the maximum speed-up is about 60.
One relies on a k-point grid of 4x4x4 x 4 shifts (=256 k points), and 5 bands.
There are three perturbations. For the two first perturbations, no symmetry can be used,
while for the the third, two symmetries can be used to reduce the number of k points to 128
Hence, for the perfectly scalable sections of the code, the maximum speed up is 128*5=640 (on 640 cores).
However, the sequential parts of the code dominate at a much lower value...
The scaling limitation is mostly due to the reading of ground state wavefunctions in the inwffil.F90
routine. This reading is not parallelized over k points/bands. In the present
case, 10k points are present in the irreducible Brillouin zone (256 in the
full BZ) and must be read.

There is one preparatory steps, before running the DFPT calculation.

Preparatory step 1
(mpirun ...)  abinit < tdfpt_01.files > tdfpt_01.log
cp tdfpt_01.o_WFK tdfpt_02.i_WFK
cp tdfpt_01.o_WFK tdfpt_02.i_WFQ

Test case, step 2 (DFPT calculation)
(mpirun ...)  abinit < tdfpt_02.files > tdfpt_02.log

-------------------------------------------

Second test : BaTiO3 slab (29 atoms), 
computation of the phonon frequencies at qpt 0.0 0.375 0.0

This test, with 29 atom, is quite slow, but scales very well.

There is one preparatory step, before running the DFPT calculation.
The preparatory step can be run on 16 processors at most with the current
input file. It might use more processors as well, with the kgb parallelism
(but the input file has to be modified).
On 8 processors, the preparatory step is about three hours.
It generates well-converged wavefunctions. For a quick trial,
simply set nstep 1   instead of nstep 50 ,
this will run in about 6 minutes.

The test case itself is an underconverged calculation of the response with
respect to one perturbation (atomic displacement). It is underconverged
because nstep has been set to 10, while more than 30 are needed.
Moreover, obtaining the interatomic force constants would need computing
many more perturbations than the present one.
In any case, the present test case run in about 45 minutes on a 8 core
machine.
Since the number of k points to be kept for the present perturbation is is 8x8x1 with 4 symmetries,
that is 16, and the number of bands is 120, the perfectly scalable part of the
test case should have a maximum speed up of 1920.

From tests for the 8 core case, on a total of 20200 secs, there
were 305 secs for vtorho3:synchro (sequential) and
260.460 for inwffil (sequential).
The latter will not increase with a bigger value of nstep, and for more
perturbations, while the former will increase proportionally.

Hence, in the present status, for 8 cores, the sequential part is about 3%,
leading to a maximum speed-up with respect to sequential, of about 240.
For a larger test case (bigger nstep, more perturbations), the maximum speed up might
be twice bigger.



Preparatory step 1
(mpirun ...)  abinit < tdfpt_03.files > tdfpt_03.log
cp tdfpt_03.o_WFK tdfpt_04.i_WFK
cp tdfpt_03.o_WFK tdfpt_04.i_WFQ

Test case, step 2 (DFPT calculation)
(mpirun ...)  abinit < tdfpt_04.files > tdfpt_04.log

-------------------------------------------

