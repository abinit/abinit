---
title: Parallelization of ground state
authors: FB
---

# Parallelism for the ground state  

## Explore the various features of the KGB parallelization scheme  

This tutorial discusses how to perform GS calculations on hundreds of processors using ABINIT.

You will learn how to use some keywords related to the "KGB" parallelization scheme where
"K" stands for "k-point", "G" refers to the wavevector of a planewave, and "B" stands for a "Band". 
On the one hand, you will see how to improve the speedup of your calculations and, on the other hand, 
how to increase the convergence rate of a self consistent field calculation.  
  
This tutorial should take about 1.5 hour and requires access to at least a 200
CPU core parallel computer.

You are supposed to know already some basics of parallelism in ABINIT,
explained in the tutorial [A first introduction to ABINIT in parallel](basepar).
On the contrary, the present tutorial is not a prerequisite for other
tutorials about parallelism in ABINIT.

[TUTORIAL_README]

## 1 Introduction

*Before continuing you might work in a different subdirectory, as for the
other tutorials. Why not Work_paral_gspw?* 

All the input files can be found in the *\$ABI_TUTOPARAL/Input* directory. 
You might have to adapt them to the path of the directory in which you have decided to perform your runs. 
You can compare your results with reference output files located in *\$ABI_TUTOPARAL/Refs*.

In the following, when "run ABINIT over nn CPU cores" appears, you have to
use a specific command line, according to the operating system and
architecture of the computer you are using. This can be, for instance:
    
```bash
mpirun -n nn abinit < abinit.files > log 2> err
``` 

or the use of a specific submission file. 
Some scripts are given as examples in the directory *\$ABI_HOME/doc/tutorial/paral_gspw_assets/*. 
You can adapt them to your own calculations.  

When the size of the system increases up to 100 or 1000 atoms, it is usually
impossible to perform ab initio calculations with a single computing core. 
This is because the basis sets used to solve the problem (PWs, bands, ...) increase
proportionally (linearly, as the square, or even exponentially...). 
The trouble has two origins: 

* the memory i.e. the amount of data stored in RAM
* the computation time, with specific bottlenecks which are, for
  example, the eigensolver, the 3dim-FFTs...

Therefore, it is generally mandatory to adopt a parallelization strategy. 
That is to say: (i) to distribute the data in a MPI-sense on a large number of
processors, and (ii) to parallelize the routines responsible for the increase
of the computation time.

In this tutorial, we will go beyond the simple k-point and spin parallelism explained 
in the tutorial [A first introduction to ABINIT in parallel](basepar) and show:  

* how to improve performance by using a large number of processors, even when the number of k-points is not large,

* how to decrease the computation time for a given number of processors. The aim is twofold:  

    1. reduce the time needed to perform one electronic iteration (to improve the efficiency)  
    2. reduce the number of electronic iterations (to improve the convergence) 

* how to use other keywords or features related to the KGB parallelization scheme.  

The tests are performed on a system with 108 atoms of gold; this is a
benchmark used for a long time during the development and the implementation of the KGB parallelization. 
With respect to the input file used for the benchmark, the cutoff energy is strongly reduced in this tutorial, for
practical reasons. For real tests, you can see the results (in particular the scaling) in:

  * the publication concerning the KGB parallelisation: [[cite:Bottin2008]]
  ([available on Arxiv.org](https://arxiv.org/abs/0707.3405)). 

  * the Abinit paper [[cite:Gonze2009]] 

  * the presentation of F. Bottin at the [ABINIT workshop 2007](https://www.abinit.org/sites/default/files/oldsites/workshop_07/program.html) (Monday 29, session 2).

You are strongly suggested to read these documents before beginning this tutorial. 
  
However, even if you read them with attention, you won't learn the answer to
the most frequently asked question. Why this parallelization is named KGB? We
don't know. Some people say that the reason comes from the K-point, the plane
waves G and the Bands, but you can imagine everything you want.  

## 2 A simple way to begin...
  
One of the most simple way to launch the KGB parallelization in ABINIT is to
add just one input variable to the sequential input file. 
This is [[paral_kgb]] and controls everything concerning the KGB parallelization: the
use of the LOBPCG eigensolver ([[wfoptalg]] = 4 or 14 or 114) of A. Knyazev, the
parallel 3dim-FFT ([[fftalg]] = 401) written by S. Goedecker, and some other tricks... 
Then you have to chose between tuning you parallelism with the 3 variables
to define the number or processors needed on each level of the KGB parallelization: [[npband]], [[npfft]] and [[npkpt]],
**or** use the not so good input variable [[autoparal]] = 1 with 
[[max_ncpus]]="the maximum number of processors you want" which will do the work for you.
This tutorial is mainly aimed at learning how to chose the [[npband]], [[npfft]] and [[npkpt]] variables.

But first things first, copy the file *tgspw_01.in* and the
related *tgspw_01.files* file from *\$ABI_TUTOPARAL/Input* to your working directory. 

```sh
cd $ABI_TUTOPARAL/Input
mkdir Work_paral_gspw
cd Work_paral_gspw
cp ../tgspw_01.files . 
cp ../tgspw_01.in .
```

{% dialog tests/tutoparal/Input/tgspw_01.files tests/tutoparal/Input/tgspw_01.in %}

Then run ABINIT on 1 CPU core.
At the end of the log file *tgspw_01.log*, you will see:  
    
```
   npimage|       npkpt|    npspinor|       npfft|      npband|     bandpp |       nproc|      weight|  
1 - >    1|   1 ->    1|   1 ->    1|   1 ->   22|   1 ->  108|   1 ->   65|   8 ->  108|   1 ->  108|  
         1|           1|           1|          12|           9|           1|         108|      55.13 |  
         1|           1|           1|          10|           9|           1|          90|      54.98 |  
         1|           1|           1|           9|           9|           1|          81|      54.85 |  
         1|           1|           1|           9|          12|           1|         108|      54.37 |  
         1|           1|           1|          12|           9|           2|         108|      52.86 |  
......................  
```

A weight is assigned to each distribution of processors. As indicated in the
documentation, you are advised to select a processor distribution with a
higher weight. If we just focus on [[npband]] and [[npfft]], we see that for
108 processors the recommended distribution is (12x9).  
  
In a second step you can launch ABINIT in parallel on 108 processors by
changing your input file as follows:  
  
```diff
- paral_kgb 1 autoparal 1 max_ncpus 108  
+ paral_kgb 1 npband 12 npfft 9  
```
  
You can now perform your calculations using the KGB parallelization in ABINIT.
But somehow, you did it without understanding how you got the result...  

## 3 ... which is however coherent with a more sophisticated method
  
In this part we will try to recover the previous result, but with a more
comprehensive understanding of what is happening. As shown above, the couple
([[npband]] x [[npfft]]) of input variables can have various values: (108x1),
(54x2), (36x3), (27x4), (18x6) and (12x9). But also (9x12) ... which is not
indicated.
In order to perform these seven calculations you can use the input file 
*tgspw_02.in* and *tgspw_02.files*

{% dialog tests/tutoparal/Input/tgspw_02.files tests/tutoparal/Input/tgspw_02.in %}

and change the line corresponding to the processor distribution.
A first calculation with:  

```diff
+ npband 108 npfft 1  
```
  
A second one with another distribution:  
  
```diff
- npband 108 npfft 1  
+ npband  54 npfft 2  
```
    
And so on... Alternatively, this can be performed using a shell script including:  

```bash
cp tgspw_02.in tmp.file  
echo "npband 108 npfft 1" >> tgspw_02.in  
mpirun -n 108 abinit < tgspw_02.files  
cp tgspw_02.out tgspw_02.108-01.out  
cp tmp.file tgspw_02.in  
echo "npband 54 npfft 2" >> tgspw_02.in  
...   
```

By reference to the couple ([[npband]] x [[npfft]]), all these results are named: 
*tgspw_02.108-01.out*, *tgspw_02.054-02.out*, *tgspw_02.036-03.out*,
*tgspw_02.027-04.out*, *tgspw_02.018-06.out*, *tgspw_02.012-09.out* and *tgspw_02.009-12.out*. 
The timing of each calculation can be retrieved using the shell command:  
    
```bash
grep Proc *out  
tgspw_02.009-12.out:- Proc.   0 individual time (sec): cpu=         88.3  wall=         88.3  
tgspw_02.012-09.out:- Proc.   0 individual time (sec): cpu=         75.2  wall=         75.2  
tgspw_02.018-06.out:- Proc.   0 individual time (sec): cpu=         63.7  wall=         63.7  
tgspw_02.027-04.out:- Proc.   0 individual time (sec): cpu=         69.9  wall=         69.9  
tgspw_02.036-03.out:- Proc.   0 individual time (sec): cpu=        116.0  wall=        116.0  
tgspw_02.054-02.out:- Proc.   0 individual time (sec): cpu=        104.7  wall=        104.7  
tgspw_02.108-01.out:- Proc.   0 individual time (sec): cpu=        141.5  wall=        141.5  
```
    
As far as the timing is concerned, the best distributions are then the ones
proposed above in section 2.: that is to say the couples (18x6) and (27x4). 
So the prediction was surprisingly pretty good.  

Up to now, we have not learned more than before. We have so far only
considered the timing (the efficiency) of one electronic step, or 10
electronic steps as this is limited in the input file. However, when the
[[npband]] value is modified, the size of the block in LOBPCG changes, and
finally the solutions of this blocked eigensolver are also affected. In other
words, we never had in mind that the convergence of these calculations is also
very important. One calculation can be the quickest if we look at time required
for one iteration but the slowest at the end of the SCF cycle because it takes many more steps.
In order to see this without performing any additional calculations, we can have
a look at the degree of convergence at the end of the calculations we already
have. The last iterations of the SCF loop give:  
    
```bash
grep "ETOT 10" *.out  
tgspw_02.009-12.out: ETOT 10  -3754.4454784191    -1.549E-03 7.222E-05 1.394E+00 
tgspw_02.012-09.out: ETOT 10  -3754.4458434046    -7.875E-04 6.680E-05 2.596E-01 
tgspw_02.018-06.out: ETOT 10  -3754.4457793663    -1.319E-03 1.230E-04 6.962E-01 
tgspw_02.027-04.out: ETOT 10  -3754.4459048995    -1.127E-03 1.191E-04 5.701E-01 
tgspw_02.036-03.out: ETOT 10  -3754.4460493339    -1.529E-03 7.121E-05 3.144E-01 
tgspw_02.054-02.out: ETOT 10  -3754.4460393029    -1.646E-03 7.096E-05 7.284E-01 
tgspw_02.108-01.out: ETOT 10  -3754.4464631635    -6.162E-05 2.151E-05 7.457E-02
```

The last column indicates the convergence of the density (or potential)
residual. You can see that this quantity is the smallest when [[npband]] is
the highest. This result is expected: the convergence is better when the size
of the block is as large as possible. But the best convergence is obtained for
the (108x1) distribution... when the worst timing is measured.  
  
So, you face a dilemma. The calculation with the smallest number of iterations
(the best convergence) is not the best concerning the timing of one iteration
(the best efficiency), and vice versa... So you have to check both of these
features for all the processor distributions. On one hand, the timing of one
iteration and, on the other hand, the number of iterations needed to converge.
The best choice is a compromise between them, not necessarily the independent optima.  
  
In the following we will choose the (27x4) couple, because this one
definitively offers more guarantees concerning the convergence and the
timing even if the (18x6) one is slightly quicker per electronic step.  
  
Note: you can verify that the convergence is not changed when the [[npfft]]
value is modified. The same results will be obtained, step by step.  

## 4 Meaning of bandpp: part 1 (convergence)
  
We have seen in the previous section that the best convergence is obtained
when the size of the block is the largest. This size is related to the
[[npband]] input variable. But not only. It is possible to increase the size
of the block without increasing drastically the number of band processors.
This means that it's possible to decrease the number of electronic steps
without increasing strongly the timing of one electronic step. For systems
with peculiar convergence, when some trouble leads the calculation to diverge,
this is not just useful but indispensable to converge at all.  
  
The input variable enabling an increasing of the size block without increasing
the number of band processors is named [[bandpp]]. The size block is then
defined as: [[bandpp]] x [[npband]]. In the following, we keep the same input
file as previously and add:  
    
```diff
- nstep 10  
+ nstep 20  
+ npband 27 npfft 4  
```

You can copy the input file *tgspw_03.in* then run ABINIT over 108 cores with
*tgspw_02.027-04.out*, [[bandpp]] = 1 and [[bandpp]] = 2. The output files will be
named *tgspw_03.bandpp1.out* and *tgspw_03.bandpp2.out*, respectively. A
comparison between these two files shows that the convergence is better in the
second case. The convergence is even achieved before the input file limit of
20 electronic iterations. Quod Erat Demonstrandum:  
  
For a given number of processors, it is possible to improve the convergence by increasing bandpp.  
  
We can also compare the result obtained for the (27x4) distribution and
[[bandpp]] = 2 with the (54x2) one and [[bandpp]] = 1.
Use the same input file and add:
    
```diff
- npband 27 npfft 4  
+ npband 54 npfft 2 bandpp 1  
```

Then run ABINIT over 108 cores and copy the output file to
*tgspw_03.054-02.out*. Perform a diff (with vimdiff for example) between the two
output files *tgspw_03.bandpp1.out* and *tgspw_03.054-02.out*. You can convince
yourself that the two calculations (54x2) with [[bandpp]] = 1 and (27x4) with
[[bandpp]] = 2, give exactly the same convergence. This result is expected,
since the sizes of the block are equal (to 54) and the number of FFT
processors [[npfft]] does not affect the convergence.  

It is possible to modify the distribution of processors, without changing the
convergence, by reducing [[npband]] and increasing [[bandpp]] proportionally.  

## 5 Meaning of [[bandpp]]: part 2 (efficiency)
  
In the previous section, we showed that the convergence doesn't change if
[[bandpp]] and [[npband]] change in inverse proportions. What about the
influence of [[bandpp]] if you fix the distribution?
Two cases have to be treated separately.

You can see, in the previous calculations of section 4, that the timing
increases when [[bandpp]] increases:
    
```bash
grep Proc tgspw_03.bandpp1.out tgspw_03.bandpp2.out  
tgspw_03.bandpp1.out:- Proc.   0 individual time (sec): cpu=        121.4  wall=        121.4  
tgspw_03.bandpp2.out:- Proc.   0 individual time (sec): cpu=        150.7  wall=        150.7  
```

while there are fewer electronic iterations for [[bandpp]] = 2 (19) than for
[[bandpp]] = 1 (20). If you perform a diff between these two files, you will see
that the increase in time is essentially due to the section "zheegv-dsyegv".
    
```bash
grep "zheegv-dsyegv" tgspw_03.bandpp1.out tgspw_03.bandpp2.out  
tgspw_03.bandpp1.out:- zheegv-dsyegv               1321.797  10.1   1323.215  10.1         513216  
tgspw_03.bandpp2.out:- zheegv-dsyegv               5166.002  31.8   5164.574  31.8         244944  
```

The "zheegv-dsyegv" is a part of the LOBPCG algorithm which is performed in
sequential, so the same calculation is done on each processor. In the second
calculation, the size of the block being larger (27x2=54) than in the first
(27), the computational time of this diagonalization is more expensive. To sum
up, the timing of a single electronic step increases by increasing [[bandpp]],
but the convergence improves.  

Do not increase too much the bandpp value, unless you decrease proportionally
[[npband]] or if you want to improve the convergence whatever the cost in total timing.  
  
The only exception is when [[istwfk]] = 2, i.e. when real wavefunctions are
employed. This occurs when the Gamma point alone is used to sample the
Brillouin Zone. You can use the input file *tgspw_04.in* in order to check that.
The input is modified with respect to the previous input files in order to be
more realistic and use only the Gamma point. Add [[bandpp]] = 1,2,4 or 6 in the
input file *tgspw_04.in* and run ABINIT in each case over 108 cores. You will
obtain four output files named *tgspw_04.bandpp1.out*, *tgspw_04.bandpp2.out*,
*tgspw_04.bandpp4.out* and *tgspw_04.bandpp6.out* in reference to [[bandpp]]. If
you compare the outputs of these calculations:  
    
```bash
grep Proc *out  
tgspw_04.bandpp1.out:- Proc.   0 individual time (sec): cpu= 61.4  wall= 61.4  
tgspw_04.bandpp2.out:- Proc.   0 individual time (sec): cpu= 49.0  wall= 49.0  
tgspw_04.bandpp4.out:- Proc.   0 individual time (sec): cpu= 62.5  wall= 62.5  
tgspw_04.bandpp6.out:- Proc.   0 individual time (sec): cpu= 75.3  wall= 75.3
```

you can see that the timing decreases for [[bandpp]] = 2 and increases
thereafter. This behaviour comes from the FFTs. For [[bandpp]] = 2, the real
wavefunctions are associated in pairs in the complex FFTs, leading to a
reduction by a factor of 2 of the timing in this part (you can see this
reduction by a diff of the output files). Above [[bandpp]] = 2, there is longer
any gain in the FFTs, whereas some significant losses in computational time
appear in the "zheegv-dsyegv" section.  
  
When calculations are performed at the Gamma point, you are strongly
encouraged to use bandpp = 2... or more if you need to improve the convergence whatever the timing.  

## 6 The KGB parallelization
  
Up to now, we only performed a GB parallelization. This implies
parallelization over 2 levels of PWs or over 2 levels of bands and FFTs, for
different sections of the code (see the [paper](https://arxiv.org/abs/0707.3405) or
[presentation](https://www.abinit.org/sites/default/files/oldsites/workshop_07/home.html)). 
If the system has more than 1 k-point, one can add a third level of parallelization
and perform a real KBG parallelization. There is no additional difficulty in
adding processors on this level. In order to explain the procedure, we restart
with the same input file that was used in section 3 and add a denser M-P grid
(see the input file *tgspw_05.in*). In this case, the system has 4 k-points in
the IBZ so the calculation can be parallelized over (at most) 4 k-point
processors. This is done using the [[npkpt]] input variable:  
    
```diff
+ npkpt 4  
```

This implies we use four times more processors than before, so run ABINIT over
432 CPU cores. The timing obtained in the output file *tgspw_05.out*:
    
```bash
grep Proc tgspw_05.out  
- Proc.   0 individual time (sec): cpu=         87.3  wall=         87.3
```

is quasi-identical to the one obtained for 1 k-point (69.9 sec, see the output
file *tgspw_02.027-04.out*. This means that a calculation 4 times larger (due to
an increase of the number of k-points) gives approximatively the same human
time if you parallelize over all the k-points. You have just re-derived a well
established result: the scaling (the speedup) is quasi-linear on the k-point level.  

!!! importan

    When you want to parallelize a calculation, begin by the k-point level, then
    follow with the band level (up to [[npband]] = 50 typically) then activate the FFT level.  

Here, the timing obtained for the output *tgspw_05.out* leads to a hypothetical
speedup of 346, which is good, but not 432 as expected if the scaling was
linear as a function of the number of the k-point processors. Indeed, in order
to be comprehensive, we have to mention that the timing obtained in this
output is slightly longer (17 sec. more) than the one obtained in
*tgspw_02.027-04.out*.  Compare the time spent in all the routines. A first clue
comes from the timing done below the "vtowfk level", which contains 99% of
sequential processor time:  
    
```bash
grep "vtowfk   " tgspw_05.out tgspw_02.027-04.out  
tgspw_05.out:- vtowfk          26409.565  70.7  26409.553  70.7  4320  
tgspw_05.out:- vtowfk          26409.565  70.7  26409.553  70.7  4320  
tgspw_02.027-04.out:- vtowfk    6372.940  84.8   6372.958  84.8  1080  
tgspw_02.027-04.out:- vtowfk    6372.940  84.8   6372.958  84.8  1080  
```

We see that the KGB parallelization performs really well, since the wall time
spent within vtowfk is approximatively equal: $26409/432\approx 6372/108$. So, the
speedup is quasi-linear below vtowfk. The problem comes from parts above
vtowfk which are not parallelized and are responsible for the negligible
(1-99.xyz)% of time spent in sequential. These parts are no longer negligible
when you parallelize over hundreds of processors.  
  
You can also prove this using the percentages rather than the values of the
overall times (shown above): the time spent in vtowfk corresponds to 84.8% of
the overall time when you don't parallelize over k-points, and only 70.7% when
you parallelize. This means you lose time above vtowfk in this case.

This behaviour is related to the [Amdhal's law](https://en.wikipedia.org/wiki/Amdahl%27s_law): 
> The speedup of a program using multiple processors in parallel computing is 
> limited by the time needed for the sequential fraction of the program. 

For example, if a program needs 20 hours using a single processor core, 
and a particular portion of 1 hour cannot
be parallelized, while the remaining promising portion of 19 hours (95%) can
be parallelized, then regardless of how many processors we devote to a
parallelized execution of this program, the minimum execution time cannot be
less than that critical 1 hour. Hence the speedup is limited up to 20."

In our case, the part above the loop over k-point in not parallelized by the
KGB parallelization. Even if this part is very small, less than 1%, when the
biggest part below is strongly (perfectly) parallelized, this small part
determines an upper bound for the speedup.  
  
To do in the future: Discuss convergence and [[wfoptalg]], [[nline]]
