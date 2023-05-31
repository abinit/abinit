---
authors: T. Rangel, O. Gingras
---

# Wannier90 postprocessor

## The Wannier90 interface tutorial

This tutorial aims at showing how to use the Wannier90 interface to compute
Maximally Localized Wannier Functions (MLWFs).

You will learn how to get MLWFs with ABINIT and Wannier90 and what are the
basic variables to govern the numerical efficiency.

This tutorial should take about 2 hours and it is important to note that the examples in this tutorial
are not converged, they are just examples to show how to use the code.

[TUTORIAL_README]

## 1. Summary of Wannier90 in ABINIT

Wannier90 is a code that computes MLWFs (see [www.wannier.org](http://www.wannier.org) ).
Wannier90 uses the methodology introduced by N. Marzari and D. Vanderbilt in 1997 and it is
highly recommended to read the following papers to understand its basics:
[[cite:Marzari1997]] and [[cite:Souza2002a]].

Wannier functions (WFs) can be obtained from Bloch states by means of the formulas 1-3 of [[cite:Souza2002a]].
As you may note there is a freedom of choice in the Bloch orbital's phase which is reflected
in the shape and the spatial extent of the WF.
This means that for different phases there will be WFs with different spatial localizations.

To obtain the MLWFs we minimize the spread of the WF with respect to the choice of phase.
This is done by using a steepest-descent algorithm, see section D of [[cite:Marzari1997]].
After a ground state calculation the Wannier90 code will obtain the MLWFs
requiring just two ingredients:

* The overlaps $M_{mn} = \langle u_{mk} | u_{nk+b} \rangle$ between the cell periodic part
  of the Bloch states $|u_{nk}\rangle$. See Eq. 25 of [[cite:Marzari1997]]).

* As a starting guess the projection, $A_{mn} = \langle \psi_{mk} | g_{n} \rangle$,
  of the Bloch states $|\psi_{nk} \rangle$ onto trial
  localized orbitals $|g_{n}\rangle$ (See section D of [[cite:Souza2002a]])

What ABINIT does is to take the Bloch functions from a ground state calculation
and compute these two ingredients. Then, Wannier90 is executed. Wannier90 is
included as a library in ABINIT and the process is automatic, so that in a
single run you can do both the ground state calculation and the computation of MLWFs.

## 2. A first example: silicon

Before starting make sure that you compiled abinit enabling Wannier90.
You may have to recompile the code with

```
configure --with-config-file=myconf.ac9
```

where *myconf.ac9* defines:

```sh
#Install prefix of the PSML I/O library (e.g.
#                          /usr/local).
with_wannier90="/usr/local"

```

Now we will compute a set of MLWFs for silicon.
We are going to extract the Wannier functions corresponding to the four valence states of silicon.

*Before beginning, you might consider to work in a different sub-directory as
for the other tutorials. Why not Work_w90?*
Then copy the files *tw90_1.abi* and *wannier90.win* from
the *$ABI_TESTS/tutoplugs/Input* directory to *Work_w90*:

    cd $ABI_TESTS/tutoplugs/Input
    mkdir Work_w90
    cd Work_w90
    cp ../tw90_1.abi .

Wannier90 also uses a secondary input file called *wannier90.win*.
Therefore, you must include this file in the folder:

    cp ../wannier90.win .

Now you are ready to run abinit. Please type in:

    abinit  tw90_1.abi > log 2> err &

Let's examine the input file *tw90_1.abi*, while the calculation is running.

{% dialog tests/tutoplugs/Input/tw90_1.abi %}

The input file should look familiar to you. It is indeed the primitive cell of silicon.
It has two data sets: first a SCF calculation and then a NSCF calculation which
will call the Wannier90 library. The only new input variable is [[prtwant]]
which has to be set to 2 in order to use the Wannier90 utility.

Now lets look at the second input file *wannier90.win*.

{% dialog tests/tutoplugs/Input/wannier90.win %}

This is a mandatory input file required by the Wannier90 library.
There are many variables that can be defined inside this file.
In our case we used **num_wann** and **num_iter**.
These variables are used in the minimization of the spread
to obtain the MLWF. In particular, **num_wann** defines the number of Wannier
functions to extract while **num_iter** sets the maximum number of iterations. There
are also variables to govern the disentanglement procedure outlined in [[cite:Souza2002a]]
which are not used in this simple case. The complete list of input variables
can be found in the Wannier90 user guide (see [www.wannier.org](http://www.wannier.org)).

We can now examine the log file. After the convergence of the SCF cycle is
reached. We can see that the Wannier90 library is called. You will find the following lines:

      Calculation of overlap and call to Wannier90 library
      to obtain Maximally Localized Wannier functions
      - wannier90.win is a mandatory secondary input
      - wannier90.wout is the output for the library
      - wannier90.amn contains projections
      - wannier90random.amn contains random projections
      - wannier90.mmn contains the overlap
      - wannier90.eig contains the eigenvalues

This is an explanation of the input and output files for the Wannier90
library. As you can see many new files were created. The input files for
Wannier90 which were created by ABINIT are:

**wannier90random.amn**
:   Contains a list of projections to be used as a starting guess of the WF.
    This is the $A_{mn}$ matrix which was mentioned before in this tutorial.

**wannier90.eig**
:   Contains a list of eigenvalues for each k-point and band.

**wannier90.mmn**
:   Contains the overlaps between the cell periodic part of the Bloch states.
    This is the M_mn matrix mentioned before in this tutorial.

**UNK**
:   Files containing the wavefunction in real space for every k-point.
    Once these files were computed by ABINIT the Wannier90 library was used.
    The output files of Wannier90 are:

**wannier90.wout**
:   This is the main output file of the library.
    You should read it carefully to see the details of the calculation.

**wannier90.chk**
:   This file is required to restart a calculation is case you use Wannier90 in standalone mode.
    In our case it is not used.

To obtain information about the steepest-descent minimization just issue:

    grep CONV wannier90.wout

You will obtain a table of the following form:

     +--------------------------------------------------------------------+<-- CONV
     | Iter  Delta Spread     RMS Gradient      Spread (Ang^2)      Time  |<-- CONV
     +--------------------------------------------------------------------+<-- CONV
          0     0.438E+02     0.0000000000       43.7939618280       0.08  <-- CONV
          1    -0.946E+01    10.5484513508       34.3387915333       0.09  <-- CONV


You can verify that the final spread you get is around 4.0 Å$^2$.

Similarly to obtain information about the disentanglement procedure (not used in this example)
just type:

    grep DIS wannier90.wout

You will obtain a table of the following form:

     +---------------------------------------------------------------------+<-- DIS
     |  Iter     Omega_I(i-1)      Omega_I(i)      Delta (frac.)    Time   |<-- DIS
     +---------------------------------------------------------------------+<-- DIS

!!! tip

    If |AbiPy| is installed on your machine, you can use the |abiopen| script
    with the `wout` command and the `--expose` option to visualize the iterations

        abiopen.py wannier90.wout --expose -sns

    ![](wannier90_assets/abiopen_tw90_1o_DS2_w90.png)


### Visualize the Wannier functions

You can see the Wannier functions in |xcrysden| format. Just type:

    xcrysden --xsf wannier90_00001.xsf

To see the isosurface click on: Tools->Data Grid -> OK
And modify the isovalue, put, for instance, 0.3 and check the option "Render +/- isovalue" then click on OK
Alternatively, one can read the xsf file with |vesta|.
MacOsx users can use the command line:

    open wannier90_00003.xsf -a vesta

to invoke vesta directly from the terminal:

![](wannier90_assets/tw90_1o_DS2_w90_00001.png)


!!! important

    * It is important to set [[istwfk]] equal to 1 for every k-point avoiding using symmetries.
      The reason is that the formalism used to extract the MLWF assumes that you have a uniform grid of k-points.
      See section IV of [[cite:Marzari1997]].

    * The number of Wannier functions to extract should be minor or equal to the number of bands.
      If _nband > nwan_ then the disentanglement routines will be called.

    * The number of k-points should be equal to ngkpt(1)*ngkpt(2)*ngkpt(3).
      This is achieved by using the input variables [[kptopt]]= 3, [[ngkpt]] and [[nshiftk]]= 1.

    * The prefix of all wannier90 files in this sample case is _wannier90_.
      Other possible prefixes are w90_ and **abo** __w90_ , where **abo** is the fourth line on your .file file.
      To setup the prefix, ABINIT will first look for a file named **abo**
      __w90.win_ if it is not found then it will look for _w90.win_ and finally for _wannier90.win_.

## 3. The PAW case

Before starting it is assumed that you have already completed the tutorials [PAW1](/tutorial/paw1) and [PAW2](/tutorial/paw2).

For silicon, we just have to add the variable [[pawecutdg]] and the PAW Atomic Data is included in the pseudopotential file.
An example has already been prepared.

Just copy the file  *tw90_2.abi* into *Work_w90*:

    cp ../tw90_2.abi .

We are going to reuse the wannier90.win of the previous example.
Now, just run abinit again

    abinit tw90_2.abi > log 2> err &

As it is expected, the results should be similar than those of the PW case.

!!! important

    For the PAW case the UNK files are not those of the real
    wavefunctions. The contribution inside the spheres is not computed, however,
    they can be used to plot the Wannier functions.

## 4. Defining the initial projections

Up to now we have obtained the MLWF for the four valence bands of silicon. It
is important to note that for valence states the MLWF can be obtained starting
from a random initial position. However, for conduction states we have to give
a very accurate starting guess to get the MLWF.

We are going to extract the $sp^3$ hybrid orbitals of Silane SiH$_4$. You can start
by copying from the tests/tutoplugs directory the following files:

    cp ../tw90_3.abi .
    cp ../tw90_3o_DS2_w90.win .

Now run abinit

    abinit  tw90_3.abi > log 2> err &

While it is running, we can start to examine the input files.
Open the main input file *tw90_3.abi*. The file is divided into three datasets.

{% dialog tests/tutoplugs/Input/tw90_3.abi %}

First a SCF calculation is done. What follows is a NSCF calculation including
more bands. Finally, in the third dataset we just read the wavefunction from
the previous one and the Wannier90 library is called. **w90iniprj** is a keyword
used to indicate that the initial projections will be given in the **.win** file.

!!! note

    You may notice that the **.win** file is now called *tw90_3o_DS2_w90.win*.
    It has the following form: *prefix_dataset_w90.win*, where the prefix is taken from
    the third line of your .file file. and dataset is the dataset number
    at which you call Wannier90 (dataset 2 in this example).

Now open the **.win** file. The initial projections will be the $sp^3$ hybrid
orbitals centered in the position of the silicon atom.
This is written explicitly as:

    begin projections
    Si:sp3
    end projections

There is an enormous freedom of choice for the initial projections. For
instance, you can define the centers in Cartesian coordinates or rotate the
axis. Please refer to the Wannier90 user guide and see the part related to
projections (see [www.wannier90.org](http://www.wannier.org)).

## 5. More on Wannier90 + ABINIT

Now we will redo the silicon case but defining different initial projections.

This calculation will be more time consuming, so you can start by running the
calculation while reading:

    cp ../tw90_4.abi .
    cp ../tw90_4o_DS3_w90.win .
    abinit  tw90_4.abi > log 2> err &

**Initial projections:**

In this example we extract sp$^3$ orbitals centered on the silicon atoms. But you
could also extract bonding and anti-bonding orbitals by uncommenting and
commenting the required lines as it is indicated in *tw90_4o_DS3_w90.win*.

You can see that we are using **r=4** in the initial projections block. This
is to indicate that the radial part will be a Gaussian function whose width
can be controlled by the value of the variable **zona**. The main advantage over
radial functions in the form of hydrogenic orbitals is that the time to write
the .amn file will be reduced.

**Interpolated band structure**

We are going to run Wannier90 in standalone mode. Just comment out the following
lines of the **.win** file:

    postproc_setup = .true.   !used to write .nnkp file at first run
    num_iter = 100

And uncomment the following two lines:

    !restart = plot
    !bands_plot = .true.

Now run Wannier90:

    $ABI_HOME/fallbacks/exports/bin/wannier90.x-abinit tw90_4o_DS3_w90

The interpolated band structure is in *tw90_4o_DS3_w90_band.dat*

To plot the bands, open |gnuplot| in the terminal and type:

```
gnuplot "tw90_4o_DS3_w90_band.gnu"
```

![](wannier90_assets/load_tw90_4o_DS3_w90_band.png)


## 6. A second example: La2CuO4

As a next example, we will compute a single MLWF for the cuprate La2CuO4.
This material is the parent compound of a set of non-conventional superconductors.
DFT predicts it is a metal, yet experimentally it is a well known antiferromagnetic insulator.
The reason is that it is a strongly correlated material which undergoes a Mott
transition.

Because DFT with the LDA functional is not able to account for strong electronic
correlations, one needs to use many-body methods to correctly describe the insulating
phase of this material.
For example, one can apply DFT+DMFT to describe the Mott transition correctly.

In this tutorial, we will show how one can construct an effective Hamiltonian in the
correlated orbital subset using Wannier90.
Once this Hamiltonian is generated, DMFT can be applied to study the Mott transition.



### 6.1 Orbital characters of the band structure

Let's start by running the first calculations. Copy this file:

    cp ../tw90_5.abi .

and run abinit using

    mpirun -n 4 abinit tw90_5.abi > tw90_5.log &

This calculation should take a couple of minutes.

The most important physics at low energy derives from the orbitals near the Fermi level.
In order to correctly select the subset of correlated orbitals that we want to model
in our effective Hamiltonian and how to project it, we need to study the orbital character
of the band structure of the material.

Let's look at the input file:
{% dialog tests/tutoplugs/Input/tw90_5.abi %}

The first dataset is simply useful to compute the ground state.
DFT will predict it to be metallic, so we use an occupation option [[occopt]] for metals.

In the second dataset, we compute a band structure along with the orbital character
of the bands. We use [[pawfatbnd]] to specify that we want the L and M resolved orbitals,
meaning we want to separate the contributions for the angular momentum L and its projection M.
We select a high symmetry path using [[kptbounds]].
A few parameters are useful for the postprocessing step following, that is 
[[pawprtdos]], [[prtdos]] and [[prtdosm]].

Now before looking at the output from the calculation, let's discuss what we expect to find.
The correlated atom here is the cupper, which has a partially filled 3d shell.
The d means the orbitals have angular momentum l=2.
There are 5 such orbitals, with two spins.

In vaccuum, these 10 states are degenerate.
However, La2CuO4 here has a layered perosvkite structure, which means the cupper atom
feels a crystal field produced by the surrounding atoms.
In particular, an octaedron of oxygen with lift the degeneracy of the 3d shell into
two subshells: the t2g and eg shells. Furthermore, the elongation of this 
octaedron splits even more these subshells.
The resulting splits are illustrated in the following figure.

![](wannier90_assets/crystal_field_LCO.png)

Additionally, by counting the stoechiometry of this material, we find that there
should be 9 electrons occupying this shell, thus we expect only the x2-y2 orbital
to be partially filled, at the Fermi level.

!!! tip
    Let's now look at the resulting fatbands now.
    We use abipy with the following commands in a python script:
    
        from abipy.abilab import abiopen                                                
                                                                                    
        with abiopen("tw90_5o_DS2_FATBANDS.nc") as fb:                                 
            fb.plot_fatbands_mview(iatom=2, fact=1.5, lmax=2, ylims=[-8, 2])
    
    These are the fatbands for the atom 2 (third atom as python starts at 0), thus the
    cupper atom. First column is the s character, second are the 3 p characters and
    last column are the 5 d characters. The last one is the dx2-y2, the only one
    that is partially filled, as expected.
    
    ![](wannier90_assets/LCO_fatbands.png) 

    This block contains 17 bands, mostly formed by the hybridization of the 4 times 3 O-2p orbitals
    with the 5 Cu-3d orbitals. It includes band numbers 19 to 33.


## 6.2 Wannier90 projection

We will now use the ABINIT + Wannier90 interface to compute an effective Hamiltonian
for the Cu 3d x^2 -y^2 orbital.
First copy this file:

    cp ../tw90_6.abi .
    cp ../tw90_6o_DS3_w90.win .

and run abinit using

    mpirun -n 4 abinit tw90_6.abi > tw90_6.log &

This calculation should take a couple of minutes. It uses two files.
The first one is

{% dialog tests/tutoplugs/Input/tw90_6.abi %}

which contains four datasets. In the first one, we simply compute the density and wave function
of the ground state on a coarse grid.
In the second dataset, we use this coarse-grid ground state to do a non-self-consistent
calculation on a finer grid.
In the third dataset, we call the Wannier90 library to perform the projection,
which uses the data computed in the second dataset.
The fourth dataset computes a band structure, that we use to compare with the projected band structure
computed by Wannier90.

The important part here happens in the third set. The Wannier90 library uses the following input file:

{% dialog tests/tutoplugs/Input/tw90_6o_DS3_w90.win %}

The keywords can be understood using Wannier90's user guide found on their web page ([https://raw.githubusercontent.com/wannier-developers/wannier90/v3.1.0/doc/compiled_docs/user_guide.pdf](https://raw.githubusercontent.com/wannier-developers/wannier90/v3.1.0/doc/compiled_docs/user_guide.pdf)). 

We want to make a one band model of the dx^2 -y^2 cupper orbital, which is defined by the projection keyword.
As seen from the fatbands calculation, this orbital has overlap with many bands, across a specific energy window.
The non-important bands which can give a non-accurate representation of the band structure are excluded using the exclude_bands keyword.
The energy window is defined using the dis_win_min and dis_win_max keywords.

Moreover, this orbital has to be disentangled with other orbital, so we need the disentanglement procedure.
The number of iterations in specified by dis_num_iter.

Finally, we want to compare the resulting band structure with the one coming from the Kohn-Sham basis.
This is enforced using the bands_plot keyword, along with a defined path which corresponds to the same path
as the one calculated in ABINIT.

To plot the resulting band structures, you will need a little bit of work. Run this python script:

    python ../EIG_to_col.py

Then you can use gnuplot to plot both the band structures together to compare them, with something like:

    echo 'plot [][0:9] "tw90_6o_DS3_w90_band.dat" u ($1/3.0470170):2 tit "Wannier90" w l lw 4, "band_struct.dat" u ($1/167):2 tit "DFT" w l
    pause -1' > plot.gnu
    gnuplot plot.gnu

Performing the calculation with more converged parameters, we find the following band strcture:

![](wannier90_assets/Wannier90_vs_ABINIT.png) 

The resulting Hamiltonian is written as the tw90_6o_DS3_w90_hr.dat file.
It can be used for example using the TRIQS package to apply DMFT, thus correcting the
DFT calculation to compute correctly the strong electronic correlations.
We would then observe a Mott transition with respect to temperature.
