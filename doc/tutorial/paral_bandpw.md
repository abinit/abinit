---
title: Parallelization of ground state
authors: FB, JB, MT
---

# Parallelization of ground state calculations

## Explore the *k-points/plane waves/bands* parallelization

This tutorial discusses how to perform ground-state calculations on hundreds/thousands
of computing unit (CPUs) cores using ABINIT.

You will learn how to use some keywords related to the **KGB** parallelization scheme where
**K** stands for *k-point*, **G** refers to the wavevector of a *planewave*, and **B** stands for *Band*.
It is possible to use ABINIT with other levels of parallelism but this is not the focus of this tutorial.
You will learn how to speedup your calculations and how to improve their convergence rate.

This tutorial should take about 1.5 hour and requires access to 64 CPU core parallel computer except for the last section which requires 256 CPU cores.

You are supposed to know already some basics of ABINIT.
Some useful references: [[cite:Levitt2015]], [[cite:Bottin2008]], [[cite:Knyazev2001]]

[TUTORIAL_README]

## Introduction

*Before continuing you might work in a different subdirectory as for the other
tutorials. Why not work_paral_kgb?*

!!! important 

    In what follows, the names of files are mentioned as if you were in this subdirectory.
    All the input files can be found in the `\$ABI_TESTS/tutoparal/Input` directory.
    You can compare your results with reference output files located in `\$ABI_TESTS/tutoparal/Refs`.

    In the following, when "run ABINIT over _nn_ CPU cores" appears, you have to use
    a specific command line according to the operating system and architecture of
    the computer you are using. This can be for instance: `mpirun -n nn abinit input.abi`
    or the use of a specific submission file.


When the size of the system increases up to 100 or 1000 atoms, it is usually
impossible to perform *ab initio* calculations with a single computing core.
This is because the basis sets used to solve the problem (plane waves, bands, ...) increase
&mdash; linearly, as the square, or even exponentially &mdash;.
The computational resources are limited by two factors:

* The memory, *i.e*. the amount of data stored in RAM,
* The computing efficiency, with specific bottlenecks.

Therefore, it is mandatory to adopt a parallelization strategy:

1. Distribute or share the data across a large number of computing nodes,
2. Parallelize the time consuming routines.

In this tutorial, we will discuss:

* How to improve performance by using a large number of computing units (CPU cores),
* How to decrease the computational time for a given number of CPU cores by:

    1. Reducing the time needed to perform one electronic iteration (improve efficiency)
    2. Reducing the number of electronic iterations (improve convergence)

The tests are performed on a 107 gold atom system.
In this tutorial the plane-wave
cutoff energy is strongly reduced, for practical reasons.

## A simple way to begin: automatic distributed parallelism

The easiest way to activate the KGB parallelization in ABINIT is to
add just one input variable in the input file, [[paral_kgb]], which controls
everything concerning the KGB parallelization, including the choice of the iterative
eigensolver ([[wfoptalg]] **= 1, 4, 14, 114**) and the use of a parallel 3dim-FFT.
Then you have to choose between 2 strategies:

* Activate the &mdash; not so good &mdash; flag [[autoparal]]**=1** (automatic parallelization)
and use the associated [[max_ncpus]] variable (maximum number of CPU cores you want),

or

* Manually define the number of processes associated to each level of parallelism:
  [[np_spkpt]] (number of processes for spin and k points),
  [[npband]] (number of processes for bands),
  [[npfft]] (number of processes for plane-waves/FFT).

OK, let's start!

Copy the `tparal_bandpw_01.abi` file the tutorial directory into your working directory.

{% dialog tests/tutoparal/Input/tparal_bandpw_01.abi %}

Then run ABINIT on 1 CPU core (using 1 MPI process and 1 *OpenMP* thread).

ABINIT should stop without starting a calculation (don't pay attention to the error message).
At the end of the output file `tparal_bandpw_01.abo`, you will see:

```md
 Searching for all possible proc distributions for this input with #CPUs<=64:

 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 |    np_spkpt|       npfft|      npband|      bandpp|  #MPI(proc)|    WEIGHT|
 |    1<<    1|    1<<   22|    1<<   64|    1<<  640|    1<<   64|  <=    64|
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 |           1|           4|          16|          10|          64|    42.328|
 |           1|           8|           8|          20|          64|    41.844|
 |           1|          16|           4|          40|          64|    41.730|
 |           1|           6|          10|          16|          60|    40.093|
 |           1|          15|           4|          40|          60|    39.054|
 |           1|           4|          16|           8|          64|    39.043|
 |           1|          12|           5|          32|          60|    39.026|
 |           1|           8|           8|          16|          64|    38.598|
 |           1|          16|           4|          32|          64|    38.493|
 |           1|           3|          20|           8|          60|    38.319|
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Only the best possible choices for nproc are printed...

```

A weight is assigned to each distribution of processors. As indicated, you are advised
to select a processor distribution with a high *weight*. If we just focus on [[npband]]
and [[npfft]], we see that, for 64 MPI processes, the recommended distribution [[npband]]x[[npfft]] is 16x4.

In a second step you can launch ABINIT in parallel on 64 cores by
changing your input file as follows:

```diff
- autoparal     1
- max_ncpus     64
+ npband        16
+ npfft         4
```

You can now perform your calculations using the KGB parallelization in ABINIT.
But somehow, you did it without understanding how you got the result...

## A more sophisticated method

In this part we will try to recover the previous process distribution, but with a better
understanding. As shown above, the pair ([[npband]] x [[npfft]]) of input variables
can have various values: (64x1), (32x2), (16x4), (8x8), (4x16).
In order to perform these 5 calculations you can use the `tparal_bandpw_02.abi`.

{% dialog tests/tutoparal/Input/tparal_bandpw_02.abi %}

Change the line corresponding to the processor distribution.
A first calculation with:
```diff
- npband         16
- npfft          4
+ npband         64
+ npfft          1
```

A second one with:
```diff
- npband         64
- npfft          1
+ npband         32
+ npfft          2
```

And so on, until you have run all 5 calculations ([[npband]]x[[npfft]]) : (64x1), (32x2), (16x4), (8x8), (4x16).

Store all the output files by renaming them as follows:
The timing of each calculation can be retrieved using the shell command:
`tparal_bandpw_02.64.01.abo`, `tparal_bandpw_02.32.02.abo`,
`tparal_bandpw_02.16.04.abo`, `tparal_bandpw_02.08.08.abo` and
`tparal_bandpw_02.04.16.abo`

The timing of each calculation can be retrieved using the shell command:
```bash
grep Proc *.abo
```

```bash
tparal_bandpw_02.04.16.abo:- Proc.   0 individual time (sec): cpu=         62.9  wall=         69.4
tparal_bandpw_02.08.08.abo:- Proc.   0 individual time (sec): cpu=         57.9  wall=         64.0
tparal_bandpw_02.16.04.abo:- Proc.   0 individual time (sec): cpu=         56.0  wall=         61.8
tparal_bandpw_02.32.02.abo:- Proc.   0 individual time (sec): cpu=         57.1  wall=         62.8
tparal_bandpw_02.64.01.abo:- Proc.   0 individual time (sec): cpu=         60.7  wall=         66.4
```

As far as the timing is concerned, the best distributions are the ones
proposed in previous section: (16x4) seems to be the best one.
The prediction using [[autoparal]]**=1** was pretty good.

Up to now, we have not learned more than before. We have so far only
considered the timing of 10 electronic steps.

However the *Locally Optimal Block Preconditioned Conjugate Gradient* algorithm (LOBPCG)
&mdash; used in ABINIT by default when [[paral_kgb]] == 1 &mdash; operates a diagonalization by block of eigenvectors.
Each block of eigenvectors is concurrently diagonalized by the [[npband]] processes,
one block after the other. When the [[npband]] value is modified, the size of the block
changes accordingly (it is exactly equal to `npband`), and the solutions of the eigensolver are modified.
One calculation can be the quickest if we look at the time needed by one iteration
but the slowest at the end because many more steps are performed.
In order to see this, we can have a look at the convergence rate at the end
of the calculations.
The last iterations of the SCF loops are:

```bash
grep "ETOT  5" *.abo
tparal_bandpw_02.04.16.abo: ETOT  5  -3654.9080524851    -1.683E-04 6.622E-05 2.509E-04
tparal_bandpw_02.08.08.abo: ETOT  5  -3654.9081359132    -2.710E-04 6.111E-05 2.318E-04
tparal_bandpw_02.16.04.abo: ETOT  5  -3654.9082768015    -6.727E-05 5.277E-05 1.490E-04
tparal_bandpw_02.32.02.abo: ETOT  5  -3654.9081759902    -2.737E-04 2.495E-05 1.968E-04
tparal_bandpw_02.64.01.abo: ETOT  5  -3654.9083410155    -1.580E-04 6.181E-05 1.440E-04
```

The last column indicates the convergence of the **residual of the density** (because we are using PAW, otherwise it would be the residual of the potential in norm-conserving).
You can see that this quantity is the smallest when [[npband]] is
the highest. This result is expected: the convergence is better when the size
of the block is the largest. But this **best convergence** is obtained for
the (64x1) distribution... when the **worst timing** is measured.

Already at the third iteration this behaviour is observed !
```bash
grep "ETOT  3" *.abo
tparal_bandpw_02.04.16.abo: ETOT  3  -3654.8846449612    -2.877E-01 5.690E-04 7.055E-02
tparal_bandpw_02.08.08.abo: ETOT  3  -3654.8848503583    -2.884E-01 8.889E-04 7.093E-02
tparal_bandpw_02.16.04.abo: ETOT  3  -3654.8858758622    -2.798E-01 5.805E-04 6.792E-02
tparal_bandpw_02.32.02.abo: ETOT  3  -3654.8866790037    -2.689E-01 6.794E-05 6.472E-02
tparal_bandpw_02.64.01.abo: ETOT  3  -3654.8885890816    -2.528E-01 4.112E-05 5.918E-02
```

So, you face a dilemma. The calculation with the smallest number of iterations
(the best convergence) is not the best concerning the timing of one iteration
(the best efficiency), and vice versa... The best choice is a compromise.

In the following we will choose the (16x4) pair, because it
definitively offers more guarantees concerning the convergence and the timing.

!!! note

    You could check that the convergence is not changed when the [[npfft]] value is modified.


## Even more sophisticated: BANDs Per Process (bandpp)

We have seen in the previous section that the best convergence is obtained
when the size of the block is the largest. This size was exactly equal to the
[[npband]] value. It was only possible to increase the block size by increasing
the number of MPI processes.

* *Is it possible to do better?*

    The answer is **yes**! The input variable named [[bandpp]] (BANDs Per Process)
    enables to increase the block size without changing the number of processes
    dedicated to bands.

* *How does this work?*

    As previoulsy, each block of bands is diagonalized by [[npband]] MPI processes in parallel.
    But, if [[bandpp]] is activated, each process handles [[bandpp]] bands.
    The block size is now equal to `npband x bandpp`. 
    Accordingly the block size can be
    modified (usually increased) by playing with the value of [[bandpp]], without changing
    the number of MPI processes.
    Note that FFTs are still performed by the [[npfft]] MPI processes.

In the following we use the same settings as previously, just performing more electronic step:

```diff
+ nstep         10
+ bandpp        1 # This is the default value
```

{% dialog tests/tutoparal/Input/tparal_bandpw_03.abi %}

Copy the input file `tparal_bandpw_03.abi` then run ABINIT over 64 CPUs,
setting [[bandpp]]**=1** and then [[bandpp]]**=2**. The output files will be
named `tparal_bandpw_03.bandpp1.abo` and `tparal_bandpw_03.bandpp2.abo`, respectively. A
comparison of these two files shows that the convergence is better in the
second case.
Conclusion: for a given number of processors, it is possible to improve
the convergence by increasing bandpp.

However, as you can see, the CPU time per iteration
increases when [[bandpp]] increases: 

```bash
grep Proc tparal_bandpw_03.bandpp1.abo tparal_bandpw_03.bandpp2.abo
tparal_bandpw_03.bandpp1.abo:- Proc.   0 individual time (sec): cpu=         90.0  wall=         95.6
tparal_bandpw_03.bandpp2.abo:- Proc.   0 individual time (sec): cpu=         90.6  wall=         96.4
```

Now look at the last iteration
```bash
grep "ETOT 10" tparal_bandpw_03.bandpp1.abo tparal_bandpw_03.bandpp2.abo
tparal_bandpw_03.bandpp1.abo: ETOT 10  -3654.9085401615    -2.882E-08 5.100E-05 1.247E-08
tparal_bandpw_03.bandpp2.abo: ETOT 10  -3654.9085401624    -3.744E-08 2.632E-05 8.003E-09
```

With [[bandpp]]=2, the calculation is more converged ! 

We can also compare the (16x4)+[[bandpp]]**=2** distribution (16x2=32) with
the (32x2)+[[bandpp]]**=1** (32x1=32) one.
These 2 calculation have [[npband]]x[[bandpp]]=32 (and [[nband]]x[[npfft]]=64).
Thus the convergence is the same !
Use the same input file and change it according to:

```diff
- npband        16 
- bandpp        2  
- npfft         4  
+ npband        32 
+ bandpp        1  
+ npfft         2  
```

Then run ABINIT over 64 MPI processes and name the output file `tparal_bandpw_03.32.02.bandpp1.abo` 
(This is one of the calculation you already did in the previous section but with [[nstep]]=10 instead of 5).
Perform a `diff` between the two output files `tparal_bandpw_03.bandpp2.abo` 
and `tparal_bandpw_03.32.02.bandpp1.abo`. 
As you can see, the two calculations give exactly the same 
convergence rate. This was expected since, in both cases, the block sizes are equal (to 32)
and the number of FFT processors [[npfft]] does not affect the convergence.

!!! tip

    It is possible to adjust the distribution of processes, without changing the
    convergence, by reducing [[npband]] and increasing [[bandpp]] proportionally.


* *Where does this CPU time consumption come from?*

    As previously explained, each MPI processes handles `bandpp` bands **sequentially**.
    Thus the *sequential part* of the code increases when [[bandpp]] increases.

!!! note

    Be carefull ! Depending on the number of plane waves and the the number of bands, increasing [[bandpp]]
    can increase or decrease the cpu time ! Usually starting from 1, increasing [[bandpp]] first decreases or maintains the cpu time,
    and continuing increasing [[bandpp]] will then increase the cpu time. 
    Experience it by yourself (also depends on the hardware).

We will see in the next section how the use of *hybrid parallelism* can again improve this...

!!! important

    Using only MPI parallelism, the timing of a single electronic step should theoretically increase
    when [[bandpp]] increases but the convergence rate is better.


## Hybrid parallelism: MPI+*OpenMP*

In modern supercomputers, the computing units (CPU cores) are no more equally distributed.
They are grouped by **nodes** in which they share the same memory access.
In so-called *many-core* architecture CPU cores can be numerous on the same node.
You could continue to use them as if they were not sharing the memory (using MPI only)
but this is not the most efficient way to take benefit from the computer.
The best practice is to use *hybrid* parallelism, mixing distributed memory parallelism
(MPI, between nodes) and shared memory parallelism ( *OpenMP*, inside a node). As you will
see, this will also have consequences on the performance of the iterative diagonalization
algorithm (LOBPCG).
Let's try!

!!! important
    
    To use ABINIT with both MPI and OpenMP, make sure to compile abinit with
    `--enable-openmp` to activate OpenMP support
    **and** to link against a linear algebra library that support multithreading like [mkl](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl/link-line-advisor.html)

!!! note

    When we say "set `OMP_NUM_THREADS=XX`", set your environment with `export OMP_NUM_THREADS=XX` for `bash` or `setenv OMP_NUM_THREADS XX` for `csh`

!!! important

    When using threads, we have to impose [[npfft]] **= 1**.
    The best is to suppress it from the input file.

The `tparal_bandpw_04.abi` input file has the parallelization set to [[npband]]=32 and [[npfft]]=2.
Thus it uses 32 MPI processes. We have set [[bandpp]]=2.

1. Run ABINIT using 64 MPI processes and 1 *OpenMP* threads (set `OMP_NUM_THREADS=1`). Copy the output file to `tparal_bandpw_04.bandpp2.1thread.abo`
2. Set [[npfft]]=1 and run ABINIT using 32 MPI processes and 2 *OpenMP* threads (set `OMP_NUM_THREADS=2`). Copy the output file to `tparal_bandpw_04.bandpp2.2thread.abo`

!!! note

    32 MPI x 2 threads = 64 cores.


{% dialog tests/tutoparal/Input/tparal_bandpw_04.abi %}

Let's have a look at the timings and compare them :

```bash
grep Proc tparal_bandpw_04.bandpp2.1thread.abo tparal_bandpw_04.bandpp2.2thread.abo
tparal_bandpw_04.bandpp2.1thread.abo:- Proc.   0 individual time (sec): cpu=         97.0  wall=        102.5
tparal_bandpw_04.bandpp2.2thread.abo:- Proc.   0 individual time (sec): cpu=        148.3  wall=         86.5
```

As you can wee, the new output file show a larger computing time for process 0: (cpu=....) disappointing?
Not really: you have to keep in mind that this timing is for one MPI process, adding the timings
of all the *OpenMP* tasks for this process. In the pure MPI case, we thus have `97 sec.` per task;
but in the *hybrid* case, we have `148/2=74 sec.` per task. 
As you can see, the "wall= ...." time is closer to thos value and is more or less the real user time.
Therefore with 2 threads, the run last 86s whereas with only 1 thread, the run last 102s !!

For the total 64 cores, the total cpu time used by ABINIT is `97x64=6208 sec.` in the MPI case, `74*32=2368 sec.` in the hybrid case.
This is much better!
The best way to confirm that is to look at the *Wall Time* (cumulated on all tasks) at the end
of the output file:

```bash
grep Overall tparal_bandpw_04.bandpp2.1thread.abo tparal_bandpw_04.bandpp2.2thread.abo
tparal_bandpw_04.bandpp2.1thread.abo:+Overall time at end (sec) : cpu=       6419.3  wall=       6555.8
tparal_bandpw_04.bandpp2.2thread.abo:+Overall time at end (sec) : cpu=       4853.7  wall=       2766.3
```

*How does ABINIT distribute the workload?*

Each block of bands is diagonalized by [[npband]] MPI processes.
As previously, each process handles `bandpp` bands but now using the *OpenMP* tasks.
This means that `bandpp x npband` bands are computed in parallel using
`nthreads x npband` tasks (`bandpp` is supposed to be a multiple of `nthreads` to get the best performances).
This is in principle more efficient than in the pure MPI case.
Scalar products, matrix multiplications and other algebra operations are done in parallel by the *OpenMP* tasks ([[npfft]] not used).
Doing so, the communications between MPI processes are also reduced !

!!! important

    When using threads, [[bandpp]] needs not to be a multiple of threads but it is **highly recommanded**!

* *How do we choose the number of threads?*

  Well, it strongly depends on the computer architecture and the case!
  
  A computer is made of *nodes*. On each node, there are *sockets* containing a given number
  of CPU cores. All the cores of one *node* can access the RAM of all the *sockets* but this access
  is faster on their own *socket*. This is the origin of the famous *Non Uniform Memory
  Access* effect (NUMA). The number of *threads* has thus to be a divisor of the total
  number of cores in the node, but it is better to choose a divisor of the number of
  cores in a *socket*. Indeed ABINIT performance is very sensitive to NUMA effect.

In the following, let's learn how to use ABINIT on 4 *OpenMP* threads...
First we change the input file as follows:

```diff
- npband        32
- bandpp        2
+ npband        16
+ bandpp        4
```

Then we run ABINIT using 16 MPI processes and 4 threads (still 64 cores) and rename the output file `tparal_bandpw_04.bandpp4.4thread.abo`
And we compare the timing of this "4 threads" case with the previous "2 threads" case:

```bash
grep Overall tparal_bandpw_04.bandpp2.2thread.abo tparal_bandpw_04.bandpp4.4thread.abo
tparal_bandpw_04.bandpp2.2thread.abo:+Overall time at end (sec) : cpu=       4853.7  wall=       2766.3
tparal_bandpw_04.bandpp4.4thread.abo:+Overall time at end (sec) : cpu=       4423.2  wall=       1316.8
```

We again have improved ABINIT performances!

* *Can we do better?*

  In principle, yes. As previously explained, we have to increase the block size for
  the LOBPCG diagonalization algorithm. 

Let's try it, just changing the [[bandpp]] value in input file:

```diff
- bandpp        4
+ bandpp        8
```

We don't change here the number of threads (keeping 4).
And we obtain the following timings:

```diff
grep Overall tparal_bandpw_04.bandpp4.4thread.abo tparal_bandpw_04.bandpp8.4thread.abo
tparal_bandpw_04.bandpp4.4thread.abo:+Overall time at end (sec) : cpu=       4423.2  wall=       1316.8
tparal_bandpw_04.bandpp8.4thread.abo:+Overall time at end (sec) : cpu=       4704.6  wall=       1417.3
```

The new settings do not give a better result...**but** for a few second you have a better convergence
```diff
grep "ETOT 10" tparal_bandpw_04.bandpp4.4thread.abo tparal_bandpw_04.bandpp8.4thread.abo
tparal_bandpw_04.bandpp4.4thread.abo: ETOT 10  -3654.9085401627    -5.105E-09 2.274E-05 1.673E-09
tparal_bandpw_04.bandpp8.4thread.abo: ETOT 10  -3654.9085401627    -1.376E-09 1.912E-05 6.053E-10
```
Convergence has a cost...

To help you in choosing the distribution of processes/tasks, you can launch ABINIT with
the [[autoparal]]**=1** and [[max_ncpus]] keywords. [[max_ncpus]] should be equal the
total number of targeted CPU cores, i.e. `nthreads x nMPI` and you should launch ABINIT on 1 MPI
process with `OMP_NUM_THREADS=nn`.
You can try this with `max_ncpus=64` and `OMP_NUM_THREADS=4`...

!!! tip

    The rules to distribute the workload are:

    * `npband x bandpp` (size of a block) should be maximized.
       It has to divide the number of bands (`nband`)
    * `bandpp` should be a multiple of the number of *OpenMP* tasks
    * `nband` has to be a multiple of `npband x bandpp`.
    * In any case, the ideal distribution is system dependent!

## 6 The KGB parallelization

!!! note
    For this part, a cluster was used and therefore the timing for `tparal_bandpw_03.abi` is different than the previous run in the previous section.
    **Use the same cluster to run `tparal_bandpw_03.abi` and the following run.**

Up to now, we only performed a "GB" parallelization, using 2 levels (bands and/or FFT).
If the system has more than 1 *k-point*, one can add a third level of parallelization
and perform a full "KBG" parallelization. There is no difficulty in
adding processes to this level.

To test the full parallelism, we restart
with the same input file as in `tparal_bandpw_03.abi` and add a denser *k-point* grid.
In this case, the system has 4 *k-points* in
the *irreducible Brillouin zone* (IBZ) so the calculation can be parallelized over (at most) 4 *k-points*
MPI processes. This is done using the [[np_spkpt]] input variable:

```diff
-     nkpt          1
+    ngkpt         4 4 4
+ np_spkpt         4
```

We need 4 times more processes than before, so run ABINIT over
256 CPU cores (only MPI) with the `tparal_bandpw_05.abi` file.
The timing obtained in the output file `tparal_bandpw_05.abo` and `tparal_bandpw_03.abo` are:

```bash
grep Proc tparal_bandpw_03.abo tparal_bandpw_05.abo 
tparal_bandpw_03.abo:- Proc.   0 individual time (sec): cpu=         44.2  wall=         45.5
tparal_bandpw_05.abo:- Proc.   0 individual time (sec): cpu=         49.3  wall=         50.4
```
They are quasi-identical !
This means that the scalability of ABINIT is quasi-linear on the *k-point* level.

!!!important

    When you want to parallelize a calculation, begin by the k-point level, then
    follow with the band level; then activate the FFT level or OpenMP threads.

Here, the timing obtained for the output `tparal_bandpw_05.abo` leads to a hypothetical
speedup of $45.5/50.4 \times 256\approx 231$, which is good, but not 256 as expected if the scaling was
linear. Indeed, the time needed here is slightly longer (5 sec. more) than
the one obtained in `tparal_bandpw_03.abo`.
To go further, let's compare the time spent in all the routines.
All the input files you have used contain the input variable [[timopt]]=-3 which activates the timing of different subparts of the run. 
You can find the results a the end of the output files before the references.

A first clue to undestand this not perfect speedup comes from the sequential code which is not parallelized. 
This sequential part is mainly (99%) done outside the "vtowfk level".

Let's have a look at the time spend in the well parallelized subroutine `vtowfk`:

```bash
grep -e '- vtowfk *[[:digit:]]' tparal_bandpw_03.abo tparal_bandpw_05.abo
tparal_bandpw_03.abo:- vtowfk                      2574.132  89.4   2589.639  89.3            640                   0.99       0.99
tparal_bandpw_03.abo:- vtowfk                      2574.132  89.4   2589.639  89.3            640                   0.99       0.99
tparal_bandpw_05.abo:- vtowfk                     10521.057  82.3  10595.093  82.3           2560                   0.99       0.99
tparal_bandpw_05.abo:- vtowfk                     10521.057  82.3  10595.093  82.3           2560                   0.99       0.99
```

We see that the KGB parallelization performs really well, since the wall time
spent within `vtowfk` is approximatively equal: $10521.057/256\approx 2574.132/64\approx 40$.
So, the speedup is quasi-linear in `vtowfk`. The problem comes from parts outside
`vtowfk` which are not parallelized and are responsible for the negligible
(100-82.3)% of time spent in sequential. These parts are no longer negligible
when you parallelize over hundreds of processors. The time spent in `vtowfk`
corresponds to 89.4% of the overall time when you don't parallelize over k-points,
and only 82.3% when you parallelize.

This behaviour is related to the [Amdhal's law](https://en.wikipedia.org/wiki/Amdahl%27s_law):

> The speedup of a program using multiple processes in parallel computing is
> limited by the time needed for the sequential fraction of the program.

For example, if a program needs 20 hours using a single processor core,
and a particular portion of 1 hour cannot be parallelized, while the remaining
portion of 19 hours (95%) can be parallelized, then regardless of how
many processor cores we have, the minimum execution time cannot be
less than that critical 1 hour. Hence the speedup is limited to 20.

In our case, the part above the loop over *k-points* in not parallelized by the
KGB parallelization. Even if this part is very small &mdash; less than 1% &mdash;
it determines an upper bound for the speedup.

<!--
<sub><sup>To do in the future: discuss convergence, wfoptalg, nline</sup></sub>
-->
