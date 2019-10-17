---
authors: GZ, MT, EB, MJV
---

# Tutorial about the spin  

## Properties related to spin (spin polarized calculations, ferro- ferri- magnetic materials, and spin-orbit coupling).  

This tutorial aims at showing how to get the following physical properties:

* the total magnetization of a ferromagnetic material  
* the magnetization of an antiferromagnetic material  
* analyse the total density of states per spin direction    
* analyse the density of states per atom and per spin direction  
* look at the effect of spin-orbit coupling for a non magnetic system  
* non-collinear magnetism (not yet)
* spin-orbit coupling and magnetocristalline anisotropy (not yet)

You will learn to use features of ABINIT which deal with spin.  

This tutorial should take about 1.5 hour.

[TUTORIAL_README]

## 1 A ferromagnetic material: *bcc* Fe
  
*Before beginning, you might consider to work in a different subdirectory, as
for the other tutorials. Why not Work_spin?*

The file *tspin_x.files* in *\$ABI_TESTS/tutorial/Input* lists the file names and root names. 
while *tspin_1.in* is our input file.
You can copy these two files in the *Work_spin* directory with:

```sh
cd $ABI_TESTS/tutorial/Input
mkdir Work_spin
cd Work_spin
cp ../tspin_x.files .  # Change it, when needed, as usual.
cp ../tspin_1.in .
```

{% dialog tests/tutorial/Input/tspin_x.files tests/tutorial/Input/tspin_1.in %}

You can now run the calculation with:

```sh
abinit < tspin_x.files > log 2> err &
```

then you should edit the input file, and read it carefully.  
Because we are going to perform magnetic calculations, there a two new types of variables:  

* [[nsppol]]
* [[spinat]]

You can read their description in the help file. You will work at fixed
[[ecut]] (=18Ha) and k-point grid, defined by [[ngkpt]] (the 4x4x4
Monkhorst-Pack grid). It is implicit that in *real life*, you should do a
convergence test with respect to both convergence parameters (NB: one needs a
minimal cut-off to exhibit magnetic effects).
This run takes about 12 seconds on a modern PC.

We will compare the output with and without magnetization. (Hence, there are two datasets in the run)  
We now look at the output file: 
In the magnetic case, the electronic density is split into two parts, 
the "Spin-up" and the "Spin-down" parts to which correspond different Kohn-Sham
potentials and different sets of eigenvalues whose occupations are given by the
Fermi-Dirac function (without the ubiquitous factor 2)  
  
For the first k-point, for instance, we get:  

```
(no magnetization)  
occ   2.00000   1.99989   1.99989   1.22915   1.22915   0.28676   0.00000   0.00000  
(magnetic case)  
occ   1.00000   0.99999   0.99999   0.98396   0.98396   0.69467   0.00000   0.00000 (spin-up)  
      1.00000   0.99730   0.99730   0.00898   0.00898   0.00224   0.00000   0.00000 (spin-down)  
```
  
We note that the occupations are very different for up and down spins, which
means that the eigenvalues are shifted, which is in turn due to a shift of the
exchange-correlation potential, and therefore of the effective potential. 
You can indeed have a look at the output file to compare spin-up and down eigenvalues:  

    -0.48411  -0.38615  -0.38615  -0.30587  -0.30587  -0.27293   0.33747   0.33747 (up, kpt#1)  
    -0.46638  -0.32383  -0.32383  -0.21767  -0.21767  -0.20371   0.36261   0.36261 (dn, kpt#1)  
  
The magnetization density (in unit of $\mu_B$ - Bohr's magneton) is the
difference between the up and down densities. The magnetization density,
divided by the total density, is denoted "zeta". 
This quantity "zeta" can vary between -1 and 1. It is zero everywhere in the non-magnetic case. 
In the magnetic case, we can read for instance its minimal and maximal values in the output file:
  
    Min spin pol zeta= -4.8326E-02 at reduced coord.  0.7222  0.5000  0.2222  
         next min= -4.8326E-02 at reduced coord.  0.5000  0.7222  0.2222  
    Max spin pol zeta=  5.7306E-01 at reduced coord.  0.0000  0.8889  0.8889  
         next max=  5.7306E-01 at reduced coord.  0.8889  0.0000  0.8889  
  
The total magnetization, i.e. the integrated in the unit cell, is now:  
  
    Magnetization (Bohr magneton)=  1.96743463E+00  
    Total spin up =  4.98371731E+00   Total spin down =  3.01628269E+00  
  
We observe that the total density (up + down) yields 8.000 as expected.  
  
The magnetization density is not the only changed quantity. 
The energy is changed too, and we get:  
  
    etotal1  -2.4661707268E+01  (no magnetization)  
    etotal2  -2.4670792868E+01  (with magnetization)  
  
The energy of the magnetized system is the lowest and therefore energetically
favoured, as expected since bcc iron is a ferromagnet.  
Finally, one also notes that the stress tensor is affected by the
magnetization. This would also be true for the forces, for a less symmetric material.  
  
It is interesting to consider in more detail the distribution of eigenvalues
for each direction of magnetization, which is best done by looking at the
respective densities of state.  
To this end we have set [[prtdos]] = 1 in the input file, in order to obtain
the density of states corresponding to spin-up and spin-down electrons (as soon as [[nsppol]] = 2).  
The values of the DOS are in the files *tspin_1o_DS1_DOS* and *tspin_1o_DS2_DOS*
for the magnetic and non-magnetic cases respectively. We can extract the values
for use in a plotting software.  
Traditionally, in order to enhance visibility, one plots 
the DOS of minority spin electrons using negative values. 
If we compare the DOS of the magnetized system 

![](spin_assets/bccfe_mag_dos2.jpg)

and the non-magnetized system

![](spin_assets/bccfe_nonmag_dos2.jpg)

we observe that the up and down DOS have been "shifted" with respect each other.  
The integrated density of states yields the number of electrons for each spin
direction, and we see the magnetization which arises from the fact that there
are more up than down electrons at the Fermi level.  
  
That the magnetization points upwards is fortuitous, and we can get it
pointing downwards by changing the sign of the initial [[spinat]].  
Indeed, in the absence of spin-orbit coupling, there is no relation between
the direction of magnetization and the crystal axes.  
If we start with a [[spinat]] of 0, the magnetization remains 0. [[spinat]]
serves two purposes: it is a way to initially break the spin symmetry (up/down), and also
to start with a reasonable magnetic moment, close enough to the final one (in spin DFT,
as opposed to the original flavor, there can be several local minima for the total energy).  

## 2 An antiferromagnetic example: *fcc* Fe

Well sort of....  
  
Actually, fcc Fe, displays many complicated structures, in particular spin spirals. 
A spiral is characterized by a direction along an axis, an angle of
the magnetization with respect to this axis and a step after which the magnetization comes full circle.  
A very simple particular case is when the angle is 90°, the axis is <100> and
the step is the unit cell side: spin directions alternate between
planes perpendicular to the <100> axis yielding a "spiral stairway":

![](spin_assets/fcc_fe_conv.jpg)

For instance, if the atom at [x,y,0] possesses an "up" magnetization, the atom
at [x+1/2,y,1/2] would possess a down magnetization etc...  
To describe such a structure, a unit cell with two atoms is sufficient, [0,0,0] and
[1/2,0,1/2].  
The atoms will be given opposite magnetization with the help of the variable [[spinat]].  
  
Copy the file *$ABI_TESTS/tutorial/Input/tspin_2.in* in *Work_spin*.

{% dialog tests/tutorial/Input/tspin_2.in %}

This is your input file. Modify the *tspin_x.files* file accordingly. 

You can run the calculation, then you should edit the *tspin_2.in* file, and briefly
look at the two changes with respect to the file *tspin_1.in*: the 
unit cell basis vectors [[rprim]], and the new [[spinat]].  
  
Note also we use now [[nsppol]] = 1 and [[nspden]] = 2: this combination of values
is only valid when performing a strictly antiferromagnetic calculation: nspden = 2 means
that we have 2 independent components for the charge density while nsppol = 1
means that we have 1 independent component for the wave-functions.  
In that case, ABINIT uses the so-called Shubnikov symmetries, to perform
calculations twice faster than with [[nsppol]] = 2 and [[nspden]] = 2. The
symmetry of the crystal is not the full fcc symmetry anymore, since the
symmetry must now preserve the magnetization of each atom.  ABINIT is
nevertheless able to detect such symmetry belonging to the Shubnikov groups
and correctly finds that the cell is primitive, which would not be the case
if we had the same vector [[spinat]] on each atom.  
  
If we now run the calculation again, this total computation time is
approximately 30 seconds on a recent CPU.  
If we look at the eigenvalues and occupations, they are again filled with a
factor 2, which comes from the symmetry considerations alluded to above, and
not from the "usual" spin degeneracy: the potential for spin-up is equal to
the potential for spin-down, shifted by the antiferromagnetic translation
vector. Eigenenergies are identical for spin-up and spin-down, but
wavefunctions are shifted one with respect to the other.  
  
```
kpt#   1, nband= 16, wtk=  0.05556, kpt=  0.0833  0.0833  0.1250 (reduced coord)  
-0.60539  -0.47491  -0.42613  -0.39022  -0.35974  -0.34377  -0.28895  -0.28828  
-0.25314  -0.24042  -0.22943  -0.14218   0.20264   0.26203   0.26641   0.62158  
    occupation numbers for kpt#   1  
 2.00000   2.00000   2.00000   1.99997   1.99945   1.99728   1.50632 1.48106  
 0.15660   0.04652   0.01574   0.00000   0.00000   0.00000   0.00000 0.00000  
```
  
How do we know we have magnetic order?  
The density of states used for bcc Fe will not be useful since the net
magnetization is zero and we have as many up and down electrons.  
The magnetization is reflected in the existence of distinct up and down
electronic densities, whose sum is the total density and whose difference yields
the net magnetization density at each point in real space.  
  
In particular, the integral of the magnetization around each atom will give an
indication of the magnetic moment carried by this particular atom. A first
estimation is printed out by ABINIT. You can read:  

```
Integrated total density in atomic spheres:  
-------------------------------------------  
Atom  Sphere radius  Integrated_up_density  Integrated_dn_density
Total(up+dn)   Diff(up-dn)  
   1        2.00000             3.32789225             2.99093607    6.31882832    0.33695617  
   2        2.00000             2.98670724             3.32364305    6.31035029   -0.33693581  
Note: Diff(up-dn) can be considered as a rough approximation of a local magnetic moment.  
```
  
and obtain a rough estimation of the magnetic moment of each atom 
(strongly dependent on the radius used to project the charge density):  

    magnetization of atom 1= 0.33696  
    magnetization of atom 2=-0.33693  
  
But here we want more precise results...  
To perform the integration, we will use the utility *cut3d* which yields an
interpolation of the magnetization at any point in space. *cut3d* is one of the
executables of the ABINIT package and is installed together with abinit.  
For the moment cut3d is interactive, and we will use it through a very primitive script
(written in Python) to perform a rough estimate of the magnetization on each atom.  
You can have a look at the [magnetization.py program](spin_assets/magnetization.py), and note
(or believe) that it does perform an integration of the magnetization in a cube of
side acell/2 around each atom; if applicable, you might consider adjusting the
value of the "CUT3D" string in the Python script.  
  
Copy it in your *Work_spin* directory. If you run the program, by typing 

```
python magnetization.py
```

you will see the result:  

    magnetization of atom 1= 0.38212  
    magnetization of atom 2=-0.38203  
  
which shows that the magnetizations of the two atoms are really opposite.  
With the next input file *tspin_3.in*, we will consider this same problem, but
in a different way. We note, for future reference, that the total energy is:
Etotal=-4.92489592898935E+01  

## 3 Another look at *fcc* Fe
  
Instead of treating fcc Fe directly as an antiferromagnetic material, we will
not make any hypotheses on its magnetic structure, and run the calculation
like the one for fcc Fe, anticipating only that the two spin directions are going to be different.
We will not even assume that the initial spins are of the same magnitude.  
  
You can copy the file *$ABI_TESTS/tutorial/Input/tspin_3.in* to *Work_spin*.

{% dialog tests/tutorial/Input/tspin_3.in %}

This is your input file. You can modify the file *tspin_x.files* and immediately
start running the calculation. Then, you should edit it to understand its contents. 

Note the values of [[spinat]]. In this job, we wish again to characterize the magnetic structure.  
We are not going to use zeta as in the preceding calculation, but we will here
use another feature of abinit: atom and angular momentum projected densities of state.  
These are densities of states weighted by the projection of the wave functions
on angular momentum channels (that is spherical harmonics) centered on each atom of the system.  
Note that these DOS are computed with the tetrahedron method, which is rather
time consuming and produces more accurate but less smooth DOS than the smearing method. The time
is strongly dependent on the number of k-points, and we use here only a reduced set.  
(This will take about 1.5 minutes on a modern computer)  
  
To specify this calculation we need new variables, in addition to [[prtdos]] set now to 3:  

* [[natsph]]
* [[iatsph]]
* [[ratsph]]

This will specify the atoms around which the calculation will be performed, and the radius of the sphere.  
We specifically select a new dataset for each atom, a non self-consistent
calculation being run to generate the projected density of states.  
First, we note that the value of the energy is: Etotal=-4.92489557316370E+01,
which shows that we have attained essentially the same state as above.  
  
The density of states will be in the files *tspin_3o_DS2_DOS_AT0001* for the
first atom, and *tspin_3o_DS3_DOS_AT0002* for the second atom.  
We can extract the density of d states, which carries most of the magnetic
moment and whose integral up to the Fermi level will yield an estimate of the
magnetization on each atom. We note the Fermi level (echoed the file *tspin_3o_DS1_DOS*):  
  
    Fermi energy :      -0.28270392  
  
If we have a look at the integrated site-projected density of states, we can
compute the total moment on each atom. To this end, one can open the file
*tspin_3o_DS3_DOS_AT0002*, which contains information pertaining to atom 2. This
file is self-documented, and describes the line content, for spin up and spin down:  

```
# energy(Ha)  l=0   l=1   l=2   l=3   l=4    (integral=>)  l=0   l=1   l=2 l=3   l=4  
```
  
If we look for the lines containing  an energy of "-0.28250", we find  
 
up  -0.28250  0.8026  2.6082  23.3966  0.7727  0.1687  0.30  0.34  <font color="red">3.42</font>  0.04  0.01  
dn  -0.28250  0.3381  1.8716  24.0456  0.3104  0.1116  0.30  0.33  <font color="red">2.74</font>  0.04  0.01  
  
There are apparently changes in the densities of states for all the channels,
but besides the d-channels, these are indeed fluctuations. This is confirmed
by looking at the integrated density of states which is different only for the
d-channel. The difference between up and down is 0.68, in rough agreement
(regarding our very crude methods of integration) with the previous
calculation. Using a calculation with the same number of k-points for the
projected DOS, we can plot the up-down integrated dos difference for the d-channel. 

![](spin_assets/energy_diff_fccfe.jpg)

Note that there is some scatter in this graph, due to the finite number of digits (2 decimal
places) of the integrated dos given in the file *tspin_3o_DS3_DOS_AT0002*.  
  
If we now look at the up and down DOS for each atom, we can see that the
corner atom and the face atom possess opposite magnetizations, which roughly
cancel each other. The density of states computed with the tetrahedron method
is not as smooth as by the smearing method, and a running average allows for a better view.  

## 4 Ferrimagnetic (not yet)
  
Some materials can display a particular form of ferromagnetism, which also can
be viewed as non compensated antiferromagnetism, called ferrimagnetism.  
Some atoms possess up spin and other possess down spin, but the total spin magnetization is non zero.  
This happens generally for system with different type of atoms, and sometimes
in rather complicated structures such as magnetite.  

## 5 The spin-orbit coupling
  
For heavy atoms a relativistic description of the electronic structure becomes
necessary, and this can be accomplished through the relativistic LDA approximation.

### 5.1 Norm-conserving pseudo-potentials  

For atoms, the Dirac equation is solved and the 2(2l+1) l-channel
degeneracy is lifted according to the eigenvalues of the $L+S$ operator 
(l+1/2 and l-1/2 of degeneracy 2l+2 and 2l). 
After pseudization, the associated wave functions can be recovered by adding to usual pseudo-potential projectors a
spin-orbit term of the generic form $v(r).|l,s\rangle L.S \langle l,s|$. 
Not all potentials include this additional term, but the HGH type pseudopotentials do systematically.

In a plane wave calculation, the wavefunctions will be two-component
spinors, that is they will have a spin-up and a spin-down component, and these
components will be coupled. This means the size of the Hamiltonian matrix is quadrupled.

We will consider here a heavier atom than Iron: *Tantalum*.
You will have to change the "files" file accordingly, as we want to use the
potential: *73ta.hghsc*. It is a HGH pseudopotential, with semicore states.
Replace the last line of the tspin_x.files by
    
    ../../../Psps_for_tests/73ta.hghsc  

You can copy the file *$ABI_TESTS/tutorial/Input/tspin_5.in* in *Work_spin*.

{% dialog tests/tutorial/Input/tspin_5.in %}

Change accordingly the file names in *tspin_x.files*, then run the calculation.
It takes about 20 secs on a recent computer.  

The input file contains one new variable:  

  * [[nspinor]]

Have a look at it. You should also look at [[so_psp]]; it is not set explicitly here,
because the SO information is directly read from the pseudopotential file.
One could force a non-SO calculation by setting [[so_psp]] to 0.
 
In this run, we check that we recover the splitting of the atomic levels by
performing a calculation in a big box. Two calculations are launched with and
without spin-orbit.  
  
We can easily follow the symmetry of the different levels of the non spin orbit calculation:  

```
  kpt#   1, nband= 26, wtk=  1.00000, kpt=  0.0000  0.0000  0.0000 (reduced coord)  
-2.44760    
-1.46437  -1.46437  -1.46437    
-0.17045    
-0.10852  -0.10852  -0.10852  -0.10740  -0.10740  
```
  
That is, the symmetry: s, p, s, d  
After application of the spin-orbit coupling, we now have to consider twice as many levels:  

```
 kpt#   1, nband= 26, wtk=  1.00000, kpt=  0.0000  0.0000  0.0000 (reduced coord)  
-2.43258  -2.43258    
-1.67294  -1.67294  -1.35468  -1.35468  -1.35468  -1.35468   
-0.16788  -0.16788    
-0.11629  -0.11629  -0.11629  -0.11629  -0.09221  -0.09221 -0.09120  -0.09120  -0.09120  -0.09120   
```
  
The levels are not perfectly degenerate, due to the finite size of the simulation box,
and in particular the cubic shape, which gives a small crystal field splitting of the d orbitals
between $e_g$ and $t_{2g}$ states. 
We can nevetheless compute the splitting of the levels, and we obtain, for e.g. the p-channel: 1.67294-1.35468=0.31826 Ha  
  
If we now consider the 
[NIST table](https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations-tantalum) 
of atomic data, we obtain:  

    5p splitting, table: 1.681344-1.359740=0.321604 Ha  
    5d splitting, table:   .153395-.131684=0.021711 Ha  
  
We obtain a reasonable agreement.  
A more converged (and more expensive calculation) would yield:  

    5p splitting, abinit: 1.64582-1.32141=0.32441 Ha  
    5d splitting, abinit:   .09084-.11180=0.02096 Ha  

### 5.2 Projector Augmented-Wave

Within the Projector Augmented-Wave method, the usual (pseudo-)Hamiltonian can be expressed as:  

$$
H  =  K + V_{eff} + \Sigma_{ij} D_{ij}  |p_i \rangle \langle p_j|
$$
  
If the two following conditions are satisfied:  

(1) the local PAW basis is complete enough;  
(2) the electronic density is mainly contained in the PAW augmentation regions,  
  
it can be shown that a very good approximation of the PAW Hamiltonian --
including spin-orbit coupling -- is:  

$$
H  \simeq  K + V_{eff} + \Sigma (D_{ij}+D^{SO}_{ij})  |p_i \rangle \langle p_j|
$$
  
where $D^{SO}_{ij}$ is the projection of the ($L.S$) operator into the PAW augmentation regions.  
As an immediate consequence , we thus have the possibility to use the standard $p_i$ PAW projectors; 
in other words, it is possible to use the standard PAW datasets (pseudopotentials) to perform 
calculations including spin-orbit coupling.
But, of course, it is still necessary to express the wave-functions as two
components spinors (spin-up and a spin-down components).  
Let's have a look at the following keyword:  

  * [[pawspnorb]]

This activates the spin-orbit coupling within PAW (forcing [[nspinor]]=2).  
  
Now the practice:  
We consider Bismuth.  
You will have to change the "files" file accordingly, to use the new
potential *83bi.paw*. This is a PAW dataset with 5d, 6s and 6p electrons in the valence.  
Replace the last line of the tspin_x.files by:

    ../../../Psps_for_tests/83bi.paw  

You can copy the file *$ABI_TESTS/tutorial/Input/tspin_6.in* in *Work_spin*
(one Bismuth atom in a large cell). Change the file names in
*tspin_x.files* accordingly, then run the calculation. It takes about 10 seconds on a recent computer. 

{% dialog tests/tutorial/Input/tspin_6.in %}

Two datasets are executed: the first without spin-orbit coupling, the second one using  [[pawspnorb]]=1.  

The resulting eigenvalues are:  

```
 Eigenvalues (hartree) for nkpt=   1  k points:  
 kpt#   1, nband= 24, wtk=  1.00000, kpt=  0.0000  0.0000  0.0000 (reduced
coord)  
 5d   -0.93353  -0.93353  -0.93353  -0.93353  -0.82304  -0.82304  -0.82304
-0.82304 -0.82291  -0.82291  
 6s   -0.42972  -0.42972  
 6p   -0.11089  -0.11089  -0.03810  -0.03810  -0.03810  -0.03810  
```
  
Again, the levels are not perfectly degenerate, due to the finite size and non spherical 
shape of the simulation box.
We can compute the splitting of the levels, and we obtain:  

    5d-channel: 0.93353-0.82304=0.11048 Ha  
    6p-channel: 0.11089-0.03810=0.07289 Ha  
  
If we now consider the 
[NIST table](https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations-bismuth) 
of atomic data, we obtain:  

    5d-channel: 1.063136-0.952668=0.11047 Ha  
    6p-channel: 0.228107-0.156444=0.07166 Ha  
  
A perfect agreement even with a small simulation cell and very small values of plane-wave cut-offs. 
This comes from the generation of the PAW dataset, where the SOC is calculated very accurately
and for an atomic reference. The exchange correlation functional has little impact on large SOC
splittings, which are mainly a kinetic energy effect.

## 6 Rotation of the magnetization and spin-orbit coupling
  
The most spectacular manifestation of the spin-orbit coupling is the energy
associated with a rotation of the magnetisation with respect with the crystal axis. 
It is at the origin of the magneto crystalline anisotropy of paramount technological importance.  

* * *
GZ would like to thank B. Siberchicot for useful comments.
