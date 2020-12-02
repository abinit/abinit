---
title: Parallelization of ground state
authors: FB, JB, MT
---

# Parallelization of ground state calculations

## Explore the *k-points/plane waves/bands* parallelization

This tutorial discusses how to perform ground-state calculations on hundreds/thousands
of computing units (CPUs) using ABINIT.

You will learn how to use some keywords related to the **KGB** parallelization scheme where
**K** stands for "k-point", **G** refers to the wavevector of a planewave, and **B** stands for "Band".
It is possible to use ABINIT with other levels of parallelism but this is not the focus of this tutorial.
You will learn how to speedup your calculations and how to improve their convergence rate.

This tutorial should take about 1.5 hour and requires access to at least a 200
CPU core parallel computer.

You are supposed to know already some basics of ABINIT.
Some useful references: [[cite:Levitt2015]], [[cite:Bottin2008]], [[cite:Knyazev2001]]

[TUTORIAL_README]

## 1 Introduction

*Before continuing you might work in a different subdirectory, as for the
other tutorials. Why not `work_paral`?*

All the input files can be found in the `Input` directory in the directory dedicated to this tutorial.
You might have to adapt them to the path of the working directory.
You can compare your results with reference output files located in `Refs`.

!!! note

    In the following, when "run ABINIT over nn CPUs" appears, you have to
    use a specific command line or submission file, according to the operating system/architecture
    of your computer.

<!---
Some scripts are given as examples in the directory `\$ABI_HOME/doc/tutorial/paral_gspw_assets/`.
You can adapt them to your own calculations.
--->

When the size of the system increases up to 100 or 1000 atoms, it is usually
impossible to perform *ab initio* calculations with a single computing core.
This is because the basis sets used to solve the problem (plane waves, bands, ...) increase
&mdash; linearly, as the square, or even exponentially &mdash;.
The computational resources are limited by two factors:

* The memory, i.e. the amount of data stored in RAM,
* The computing efficiency, with specific bottlenecks.

Therefore, it is mandatory to adopt a parallelization strategy:

1. Distribute or share the data across a large number of computing nodes,
2. Parallelize the time consuming routines.

In this tutorial, we will discuss:

* How to improve performance by using a large number of computing units (CPU cores),
* How to decrease the computational time for a given number of CPU cores by:

    1. Reducing the time needed to perform one electronic iteration (improve efficiency)
    2. Reducing the number of electronic iterations (improve convergence)

The tests are performed on a 108 gold atom system.
In this tutorial the plane-wave
cutoff energy is strongly reduced, for practical reasons.

## 2 A simple way to begin: automatic distributed parallelism

The easiest way to activate the KGB parallelization in ABINIT is to
add just one input variable in the input file, [[paral_kgb]], which controls
everything concerning the KGB parallelization, including the choice of the iterative
eigensolver ([[wfoptalg]] **= 1, 4, 14, 114**) and the use of a parallel 3dim-FFT.
Then you have to choose between 2 strategies:

* Activate the &mdash; not so good &mdash; flag [[autoparal]]**=1** (automatic parallelization)
and use the associated [[max_ncpus]] variable (maximum number of CPU cores you want),

or

* Manually define the number of processes associated to each level of parallelism:
  [[npkpt]] (number of processes for k points),
  [[npband]] (number of processes for bands),
  [[npfft]] (number of processes for plane-waves/FFT).

OK, let's start!
Copy the `tgspw_01.in` file and the related `tgspw_01.files` from the tutorial directory
into your working directory.

{% dialog tests/tutoparal/Input/tgspw_01.files tests/tutoparal/Input/tgspw_01.in %}

Then run ABINIT on 1 CPU core (using 1 MPI process and 1 *openMP* thread).
<!--- ADDON1 --->

ABINIT should stop without starting a calculation (don't pay attention to the error message).
At the end of the log file `*log`, you will see:

```md
 Computing all possible proc distributions for this input with #CPUs<=108:

 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 |       npkpt|       npfft|      npband|      bandpp|  #MPI(proc)|    WEIGHT|
 |    1<<    1|    1<<   22|    1<<  108|    1<<  648|    1<<  108|  <=   108|
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 |           1|           6|          18|          12|         108|    74.762|
 |           1|           9|          12|          18|         108|    74.530|
 |           1|          12|           9|          24|         108|    73.687|
 |           1|          18|           6|          36|         108|    73.560|
 |           1|           4|          27|           8|         108|    73.037|
 |           1|           3|          36|           6|         108|    70.188|
 |           1|          12|           9|          18|         108|    70.065|
 |           1|           4|          27|           6|         108|    69.448|
 |           1|           8|          12|          18|          96|    67.032|
 |           1|          16|           6|          36|          96|    65.931|
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Only the best possible choices for nproc are printed...
......................
```

A weight is assigned to each distribution of processors. As indicated, you are advised
to select a processor distribution with a high "weight"". If we just focus on [[npband]]
and [[npfft]], we see that, for 108 processes, the recommended distribution is (18x6).

In a second step you can launch ABINIT in parallel on 108 processors by
changing your input file as follows:

```diff
- paral_kgb 1 autoparal 1 max_ncpus 108
+ paral_kgb 1 npband 18 npfft 6
```

<!--- ADDON2 --->

You can now perform your calculations using the KGB parallelization in ABINIT.
But somehow, you did it without understanding how you got the result...

## 3 A more sophisticated method

In this part we will try to recover the previous process distribution, but with a better
understanding. As shown above, the pair ([[npband]] x [[npfft]]) of input variables
can have various values: (108x1), (54x2), (36x3), (27x4), (18x6), (12x9) or (9x12).
In order to perform these 7 calculations you can use the `tgspw_02.in`
and `tgspw_02.files` files.

{% dialog tests/tutoparal/Input/tgspw_02.files tests/tutoparal/Input/tgspw_02.in %}

Change the line corresponding to the processor distribution.
A first calculation with:
```diff
+ npband 108 npfft 1
```

A second one with:
```diff
- npband 108 npfft 1
+ npband  54 npfft 2
```

And so on, using (36x3), (27x4), (18x6), (12x9) and (9x12).

Alternatively, this can be performed using a shell script including:
```bash
cp tgspw_02.in tmp.file
echo "npband 108 npfft 1" >> tgspw_02.in
mpirun -n 108 abinit < tgspw_02.files > log
cp tgspw_02.out tgspw_02.108-01.out
cp tmp.file tgspw_02.in
echo "npband 54 npfft 2" >> tgspw_02.in
...
```

Store all the output files by renaming them as follows:
`tgspw_02.108-01.out`, `tgspw_02.054-02.out`, `tgspw_02.036-03.out`,
`tgspw_02.027-04.out`, `tgspw_02.018-06.out`, `tgspw_02.012-09.out`
and `tgspw_02.009-12.out`.
The timing of each calculation can be retrieved using the shell command:
```
grep Proc *02*out
```

```bash
tgspw_02.009-12.out:- Proc.   0 individual time (sec): cpu=         66.1  wall=         67.8
tgspw_02.012-09.out:- Proc.   0 individual time (sec): cpu=         63.9  wall=         65.6
tgspw_02.018-06.out:- Proc.   0 individual time (sec): cpu=         61.1  wall=         62.6
tgspw_02.027-04.out:- Proc.   0 individual time (sec): cpu=         59.0  wall=         60.8
tgspw_02.036-03.out:- Proc.   0 individual time (sec): cpu=         60.6  wall=         62.3
tgspw_02.054-02.out:- Proc.   0 individual time (sec): cpu=         63.1  wall=         64.8
tgspw_02.108-01.out:- Proc.   0 individual time (sec): cpu=         74.4  wall=         76.0
```

As far as the timing is concerned, the best distributions are the ones
proposed in section 2; (27x4) seems to be the best one.
The prediction using [[autoparal]]**=1** was pretty good.

Up to now, we have not learned more than before. We have so far only
considered the timing of 10 electronic steps.
However the "Locally Optimal Block Preconditioned Conjugate Gradient" algorithm (LOBPCG)
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
grep "ETOT 10" *02*out
tgspw_02.009-12.out: ETOT 10  -4191.7097103563     3.739E-04 1.171E-05 6.693E+00
tgspw_02.012-09.out: ETOT 10  -4191.7096733078     2.111E-04 1.702E-05 6.999E+00
tgspw_02.018-06.out: ETOT 10  -4191.7102597852     1.384E-04 2.298E-05 4.033E+00
tgspw_02.027-04.out: ETOT 10  -4191.7107673046    -6.372E-06 1.772E-05 1.947E+00
tgspw_02.036-03.out: ETOT 10  -4191.7111202227    -1.791E-04 1.875E-05 1.097E+00
tgspw_02.054-02.out: ETOT 10  -4191.7113973688    -2.561E-04 1.181E-05 3.156E-01
tgspw_02.108-01.out: ETOT 10  -4191.7114913729    -4.742E-05 1.646E-05 3.480E-02
```

The last column indicates the convergence of the **residual of the potential**.
You can see that this quantity is the smallest when [[npband]] is
the highest. This result is expected: the convergence is better when the size
of the block is the largest. But this (best) convergence is obtained for
the (108x1) distribution... when the worst timing is measured.

So, you face a dilemma. The calculation with the smallest number of iterations
(the best convergence) is not the best concerning the timing of one iteration
(the best efficiency), and vice versa... The best choice is a compromise.

In the following we will choose the (27x4) pair, because it
definitively offers more guarantees concerning the convergence and the timing.

!!! note

    You could check that the convergence is not changed when the [[npfft]] value is modified.

## 4 Even more sophisticated: BANDs Per Process (bandpp)

We have seen in the previous section that the best convergence is obtained
when the size of the block is the largest. This size was exactly equal to the
[[npband]] value. It was only possible to increase the block size by increasing
the number of MPI processes.

*Is it possible to do better?*
The answer is *yes*! The input variable named [[bandpp]] (BANDs Per Process)
enables an increasing of the block size without changing the number of processes
dedicated to bands.

*How does this work?*
As previoulsy, each block of bands is diagonalized by [[npband]] MPI processes in parallel.
But, if [[bandpp]] is activated, each process handles `bandpp` bands (sequentially).
The block size &mdash; exactly equal to the number of bands handled by the *band
processes* &mdash; is now equal to `npband x bandpp`. Accordingly the block size can be
modified (usually increased) by playing with the value of [[bandpp]], without changing
the number of MPI processes.
Note also that scalar products and FFTs are done by the [[npfft]] MPI processes.

In the following we use the same settings as previously, just changing:

```diff
+ nstep 20
+ paral_kgb 1 npband 27 npfft 4 bandpp 1
```

Copy the input files `tgspw_03.files` and `tgspw_03.in` then run ABINIT over 108 CPUs,
setting [[bandpp]]**=1** and then [[bandpp]]**=2**. The output files will be
named `tgspw_03.bandpp1.out` and `tgspw_03.bandpp2.out`, respectively. A
comparison of these two files shows that the convergence is better in the
second case.
Conclusion: for a given number of processors, it is possible to improve
the convergence by increasing bandpp.

We can also compare the (27x4)+[[bandpp]]**=2** distribution with
the (54x2)+[[bandpp]]**=1** one.
Use the same input file and change it according to:

```diff
- paral_kgb 1 npband 27 npfft 4 bandpp 2
+ paral_kgb 1 npband 54 npfft 2 bandpp 1
```

Then run ABINIT over 108 CPUs and name the output file `tgspw_03.054-02.out`.
Perform a `diff` between the two output files `tgspw_03.bandpp1.out`
and `tgspw_03.054-02.out`. As you can see, the two calculations give exactly the same
convergence rate. This was expected since, in both cases, the block sizes are equal (to 54)
and the number of FFT processors [[npfft]] does not affect the convergence.

!!!tip

    It is possible to adjust the distribution of processes, without changing the
    convergence, by reducing [[npband]] and increasing [[bandpp]] proportionally.

However, as you can see in the previous calculations, the CPU time per iteration
increases when [[bandpp]] increases (note that the 2<sup>nd</sup> run performed less iterations
than the first one):

```bash
grep Proc tgspw_03.bandpp1.out tgspw_03.bandpp2.out
tgspw_03.bandpp1.out:- Proc.   0 individual time (sec): cpu=         97.4  wall=         99.3
tgspw_03.bandpp2.out:- Proc.   0 individual time (sec): cpu=         96.3  wall=         98.2
```

*Where does this CPU time consumption come from?*
As previously explained, each MPI processes handles `bandpp` bands **sequentially**.
Thus the *sequential part* of the code increases when [[bandpp]] increases.
So, [[bandpp]]**>1** is usually mandatory but do not increase it too much;
do it only if you want to improve
the convergence rate whatever the cost in total timing.

We will see in the next section how the use of *hybrid parallelism* can improve this...

!!! important

    Using only MPI parallelism, the timing of a single electronic step increases
    when [[bandpp]] increases but the convergence rate is better.

    The only exception is when [[istwfk]] = 2, i.e. when the wavefunctions are real.
    This occurs when only the &Gamma; point is used in the Brillouin Zone.
    For even values of [[bandpp]], the real wavefunctions are associated in pairs
    in the complex FFTs, leading to a reduction by a factor of their cost.
    When calculations are performed at &Gamma; point you are strongly
    encouraged to use `bandpp=2, 4,...` (even).


## 5 Hybrid parallelism: MPI+*openMP*

In modern supercomputers, the computing units (CPU cores) are no more equally distributed.
They are grouped by **nodes** in which they share the same memory access.
In so-called *many-core* architecture CPU cores can be numerous on the same node.
You could continue to use them as if they were not sharing the memory (using MPI only)
but this is not the most efficient way to take benefit from the computer.
The best practice is to used *hybrid* parallelism, mixing distributed memory parallelism
(MPI, between nodes) and shared memory parallelism ( *openMP*, inside a node). As you will
see, this will also have consequences on the performance of the iterative diagonalization
algorithm (LOBPCG).
Let's try!

We are going to run ABINIT using 2 threads. Copy the `tgspw_04.files` and `tgspw_04.in`
files and run ABINIT using 54 MPI processes and 2 *openMP* threads (by setting `OMP_NUM_THREADS=2` in your environment).
Note: `54MPI x 2threads = 108 CPU cores`.
<!--- ADDON3 --->

!!! important

    When using threads, we have to impose [[npfft]] **= 1**.
    The best is to suppress it from the input file.


Let's have a look at the timings and compare them to the same run without *openMP* threads:

```bash
grep Proc tgspw_03.bandpp2.out tgspw_04.bandpp2.2threads.out
tgspw_03.bandpp2.out:         - Proc.   0 individual time (sec): cpu=         96.3  wall=         98.2
tgspw_04.bandpp2.2threads.out:- Proc.   0 individual time (sec): cpu=        147.9  wall=        148.8
```

As you can wee, the new output file show a larger computing time for process 0: disappointing?
Not really: you have to keep in mind that this timing is for one MPI process, adding the timings
of all the *openMP* tasks for this process. In the pure MPI case, we thus have `96 sec.` per task;
but in the *hybrid* case, we have `148/2=74 sec.` per task. For the total 108 CPUS, the time
used by ABINIT is `96x108=10368 sec.` in the MPI case, `74*108=7992 sec.` in the hybrid case.
This is better!
The best way to confirm that is to look at the *Wall Time* (cumulated on all tasks) at the end
of the output file:

```bash
grep Overall tgspw_03.bandpp2.out tgspw_04.bandpp2.2threads.out
tgspw_03.bandpp2.out:         +Overall time at end (sec) : cpu=      10429.3  wall=      10607.3
tgspw_04.bandpp2.2threads.out:+Overall time at end (sec) : cpu=       8003.9  wall=       8033.4
```

*How does ABINIT distribute the workload?*
Each block of bands is diagonalized by [[npband]] MPI processes in parallel.
As previously, each process handles `bandpp` bands but now using the *openMP* tasks.
This means that `bandpp x npband` bands are computed in parallel using
`nthreads x npband` tasks (`bandpp` has thus to be a multiple of `nthreads`).
This is in principle more efficient than in the pure MPI case.
Scalar products and FFTs are done in parallel by the *openMP* tasks ([[npfft]] not used).
So, note that there are subtle differences with the pure MPI case.

!!! important

    When using threads, `bandpp` has to be a multiple of the number of threads.

    *How do we choose the number of threads?*
    Well, it strongly depends on the computer architecture!

    A computer is made of `nodes`. On each node, there are `sockets` containing a given number
    of CPU cores. All the cores of the node can access the RAM of all the `sockets` but this access
    is faster on their own `socket`. This is the origin of the famous *Non Uniform Memory
    Access* effect (NUMA). The number of `threads` has thus to be a divisor of the total
    number of CPU cores in the node, but it is better to choose a divisor of the number of
    cores in a `socket`. Indeed ABINIT performance is very sensitive to NUMA effect.

<!--- ADDON5 --->

In the following, let's learn how to use ABINIT on 6 *openMP* threads...
First we change the input file as follows:

```diff
- paral_kgb 1 npband 54 bandpp 2
+ paral_kgb 1 npband 18 bandpp 6
```

<!--- ADDON6 --->

Then we run ABINIT using 18 MPI processes and 6 threads (still 108 CPUs).
And we compare the timing of this "6 threads" case with the "2 threads" case:

```bash
grep Overall tgspw_04.bandpp2.2threads.out tgspw_04.bandpp6.6threads.out
tgspw_04.bandpp2.2threads.out:+Overall time at end (sec) : cpu=       8003.9  wall=       8033.4
tgspw_04.bandpp6.6threads.out:+Overall time at end (sec) : cpu=       5452.4  wall=       5442.8
```

We again have improved ABINIT performances!

*Can we do better?*
In principle, yes. As previously explained, we have to increase the block size for
the LOBPCG diagonalization algorithm. Let's try it, just changing the [[bandpp]]
value in input file:

```diff
- paral_kgb 1 npband 54 bandpp 2
+ paral_kgb 1 npband 18 bandpp 12
```

We don't change here the number of threads (keeping 6).
And we obtain the following timings:

```diff
grep Overall tgspw_04.bandpp6.6threads.out tgspw_04.bandpp12.6threads.out
tgspw_04.bandpp6.6threads.out: +Overall time at end (sec) : cpu=       5452.4  wall=       5442.8
tgspw_04.bandpp12.6threads.out:+Overall time at end (sec) : cpu=       6426.2  wall=       6414.9
```

The new settings do not give a better result...

To help you in choosing the distribution of processes/tasks, you can launch ABINIT with
the [[autoparal]]**=1** and [[max_ncpus]] keywords. [[max_ncpus]] should be equal the
total number of targeted CPU cores, i.e. `nthreads x nMPI` and you should launch ABINIT on 1 MPI
process with `OMP_NUM_THREADS=nn`.
You can try this with `max_ncpus=108` and `OMP_NUM_THREADS=6`...

!!! tip

    The rules to distribute the workload are:

    * `npband x bandpp` (size of a block) should be maximized.
       It has to divide the number of bands (`nband`)
    * `bandpp` has to be a multiple of the number of *openMP* tasks
    * `nband` has to be a multiple of `npband x bandpp`.
      Using [[autoparal]], the code gives you some good values for `nband`.
    * In any case, the ideal distribution is system dependent!

## 6 The KGB parallelization

Up to now, we only performed a "GB" parallelization, using 2 levels (bands and/or FFT).
If the system has more than 1 *k-point*, one can add a third level of parallelization
and perform a full "KBG" parallelization. There is no difficulty in
adding processes to this level.

To test the full parallelism, we restart
with the same input file as in section 3 and add a denser *k-point* grid.
In this case, the system has 4 *k-points* in
the *irreducible Brillouin zone* (IBZ) so the calculation can be parallelized over (at most) 4 *k-points*
MPI processes. This is done using the [[npkpt]] input variable:

```diff
- paral_kgb 1 npband 27 npfft 4 bandpp 1
+ paral_kgb 1 npkpt 4 npband 27 npfft 4 bandpp 1
```

We need 4 times more processes than before, so run ABINIT over
432 CPU cores (only MPI) with the `tgspw_05.in` and `tgspw_05.files` files.
The timing obtained in the output file *tgspw_05.out*...

```bash
grep Proc tgspw_05.out
- Proc.   0 individual time (sec): cpu=         69.8  wall=         71.9
```

is quasi-identical to the one obtained for 1 *k-point* (see `tgspw_02.027-04.out`).
This means that the scalability of ABINIT is quasi-linear on the *k-point* level.

!!!important

    When you want to parallelize a calculation, begin by the k-point level, then
    follow with the band level; then activate the FFT level or openMP threads.

Here, the timing obtained for the output `tgspw_05.out` leads to a hypothetical
speedup of 346, which is good, but not 432 as expected if the scaling was
linear. Indeed, the time needed here is slightly longer (10 sec. more) than
the one obtained in `tgspw_02.027-04.out`.
To go further, let's compare the time spent in all the routines.
A first clue comes from the timing done outside the "vtowfk level",
which contains 99% of sequential processor time:

```bash
grep "vtowfk   " tgspw_05.out tgspw_02.027-04.out
tgspw_05.out:       - vtowfk   21028.382  68.3  21042.797  67.9  4320
tgspw_05.out:       - vtowfk   21028.382  68.3  21042.797  67.9  4320
tgspw_02.027-04.out:- vtowfk    5154.681  80.8   5155.008  78.6  1080
tgspw_02.027-04.out:- vtowfk    5154.681  80.8   5155.008  78.6  1080
```

We see that the KGB parallelization performs really well, since the wall time
spent within `vtowfk` is approximatively equal: $21028/432\approx 5154/108$.
So, the speedup is quasi-linear in `vtowfk`. The problem comes from parts outside
`vtowfk` which are not parallelized and are responsible for the negligible
(1-99.xyz)% of time spent in sequential. These parts are no longer negligible
when you parallelize over hundreds of processors. The time spent in `vtowfk`
corresponds to 84.8% of the overall time when you don't parallelize over k-points,
and only 70.7% when you parallelize.

This behaviour is related to the [Amdhal's law](https://en.wikipedia.org/wiki/Amdahl%27s_law):

> The speedup of a program using multiple processes in parallel computing is
  limited by the time needed for the sequential fraction of the program.

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
