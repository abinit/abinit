---
authors: PhG,  MVeithen,  XG, NAP
---

# Tutorial on static non-linear properties  

## Electronic non-linear susceptibility, non-resonant Raman tensor, electro-optic effect.  

This tutorial aims at showing how to get the following non-linear physical properties, for an insulator:

  * The non-linear optical susceptibilities 
  * The Raman tensor of TO and LO modes 
  * The electro-optic coefficients 

We will work on AlAs. During the preliminary steps needed to compute 
non-linear properties, one will also obtain several linear response properties:

  * The Born effective charges 
  * The dielectric constant 
  * The proper piezoelectric tensor (clamped and relaxed ions) 

Finally, we will also compute the derivative of the susceptibility tensor with
respect to atomic positions (Raman tensor) thanks to finite differences.

The user should have already passed through several advanced tutorials of the
tutorial: the [tutorial Response-Function 1](rf1), the [tutorial Response-Function 2](rf2), 
the [tutorial on Polarization and finite electric field](ffield), and the 
[tutorial on Elastic properties](elastic).

This tutorial should take about 1 hour and 30 minutes.

[TUTORIAL_README]

## 1 Ground-state properties of AlAs and general parameters
  
*Before beginning, you might consider to work in a different subdirectory as for the other tutorials. 
Why not create Work-NLO in \$ABI_TESTS/tutorespfn/Input?*

In order to save some time, you might immediately start running a calculation.
Copy the file *tnlo_2.in* from *\$ABI_TESTS/tutorespfn/Input* to *Work-NLO*. Copy also
*tnlo_x.files* in *Work-NLO*, and modify it so that
all occurrences of tnlo_x are replaced by tnlo_2, then run abinit with these
data. This calculation might be one or two minutes on a PC 3GHz.

In this tutorial we will assume that the ground-state properties of AlAs have
been previously obtained, and that the corresponding convergence studies have
been done. Let us emphasize that, from our experience, accurate and trustable
results for third order energy derivatives usually require a extremely high
level of convergence. So, even if it is not always very apparent here, careful
convergence studies must be explicitly performed in any practical study you
could start after this tutorial. As usual, convergence tests must be done on
the property of interest (i.e. it is wrong to determine parameters giving
proper convergence on the total energy, and to use them blindly for non-linear
properties)

We will adopt the following set of generic parameters (the same than in the
[tutorial on Polarization and finite electric field](ffield)):
    
       acell               10.53
       ixc                 3
       ecut                2.8     (results with ecut = 5 are also reported
                                   in the discussion)
       ecutsm              0.5
       dilatmx             1.05
       nband               4 (=number of occupied bands)
       nbdbuf              0
       ngkpt               6 6 6
       nshiftk             4
       shiftk   0.5 0.5 0.5
                0.5 0.0 0.0
                0.0 0.5 0.0
                0.0 0.0 0.5
    
       pseudopotentials    13al.pspnc
                           33as.pspnc
    

In principle, the [[acell]] to be used should be the one corresponding to the
optimized structure at the [[ecut]], and [[ngkpt]] combined with [[nshiftk]]
and [[shiftk]], chosen for the calculations. Unfortunately, for the purpose of
this tutorial, in order to limit the duration of the runs, we have to work at
an unusually low cutoff of 2.8 Ha for which the optimized lattice constant is
unrealistic and equal to 7.45 Bohr (instead of the converged value of 10.64).
In what follows, the lattice constant has been arbitrarily fixed to 10.53
Bohr. For comparison, results with [[ecut]] = 5 are also reported and, in that
case, were obtained at the optimized lattice constant of 10.64 Bohr. For those
who would like to try later, convergence tests and structural optimizations
can be done using the file *\$ABI_TESTS/tutorespfn/Input/tnlo_1.in*. Before
going further, you might refresh your mind concerning the other variables:
[[ixc]], [[ecutsm]], [[dilatmx]], [[nbdbuf]].

## 2 Linear and non-linear responses from density functional perturbation theory (DFPT)
  
As a theoretical support to this section of the tutorial, you might consider
reading [[cite:Veithen2005]]:

In the first part of this tutorial, we will describe how to compute various
linear and non-linear responses directly connected to second-order and third-
order derivatives of the energy, using DFPT. From the (2n+1) theorem,
computation of energy derivatives up to third order only requires the
knowledge of the ground-state and first-order wavefunctions, see 
[[cite:Gonze1989]] and [[cite:Gonze1995]].
Our study will therefore include the following steps : (i) resolution
of the ground-state problem, (ii) determination of the first-order
wavefunctions and construction of the related databases for second and third-
order energy derivatives, (iii) combination of the different databases and
analysis to get the physical properties of interest.

This is closely related to what was done for the dielectric and dynamical
properties, except that an additional database for third-order energy
derivatives will now be set up during the run. Only selected third-order
derivatives are implemented at this stage and concern responses to electric
field and atomic displacements:

  * non-linear optical susceptibilities 
    (related to a third-order derivative of the energy with respect to electric fields at clamped nuclei positions)

  * Raman susceptibilities (mixed third-order derivative of the energy, twice with respect 
    to electric fields at clamped nuclei positions, and once with respect to atomic displacement)

  * Electro-optic coefficients (related to a third-order derivative of the energy with respect to electric fields, 
    two of them being optical fields -clamped nuclei positions- and one of them being a static field -the ions are allowed to move-)

Many different steps can be combined in one input file. For the sake of
clarity we have decomposed the calculation into individual inputs that are now described.

**Responses to electric fields and atomic displacements.**

Let us examine the file *tnlo_2.in*. Its purpose is to build databases for
second and third energy derivatives with respect to electric fields and atomic
displacements. You can edit it. It is made of 5 datasets. The first four data
sets are nearly the same as for a usual linear response calculation: (1)
self-consistent calculation in the IBZ; (2) non self-consistent calculations
to get the wave-functions over the full BZ; (3) ddk calculation, (4)
derivatives with respect to electric field and atomic displacements. Some
specific features must however be explicitly specified in order to prepare the
non-linear response step (dataset 5). First, from dataset 2 it is mandatory to specify:
    
            nbdbuf  0
            nband   4 (= number of valence bands)
    
Also, in dataset 4, it is required to impose [[prtden]], and [[prepanl]]
    
            prtden4    1
            prepanl4   1

The purpose for this is (i) to constrain [[kptopt]] = 2 even in the
computation of phonons where ABINIT usually take advantages of symmetry
irrespective of kptopt and (ii) compute the electric field derivatives in the
3 directions of space, irrespective of the symmetry.

If it was not done at the beginning of this tutorial, you can now make the
run. You can have a quick look to the output file to verify that everything is
OK. It contains the values of second and third energy derivatives. It has
however no immediate interest since the information is not presented in a very
convenient format. The relevant information is in fact also stored in the
database files (DDB) generated independently for second and third energy
derivatives at the end of run steps 4 and 5. Keep these databases that will be
used later for a global and convenient analysis of the results using ANADDB.

**Responses to strain.** We combine the above-mentioned computation of the
response to electric field and atomic displacements with the response to
strain. This is not at all mandatory for the computation of the presently
accessible non-linear response coefficients. However, this was used in 
[[cite:Veithen2005]], already mentioned,
to add corrections corresponding to
free boundary conditions, thanks to a further finite difference calculation on
top of linear response calculations. The DFPT implementation of the
computation of this correction is not available at present.

You can now copy the file *\$ABI_TESTS/tutorespfn/Input/tnlo_3.in* in *Work-NLO*,
and modify the *tnlo_x.files* accordingly (or create a file *tnlo_3.files* -
in any case, this new file should contain tnlo_3 instead of tnlo_x or tnlo_2).
You can launch the calculation, it might last about 1 minute on a PC 3 GHz.
The purpose of this run is to build databases for second energy derivatives
with respect to strains. You can edit tnlo_3.in . It is made of 4 datasets:
(1) self-consistent calculation in the IBZ; (2) non self-consistent
calculations to get the wave-functions over the full BZ; (3) ddk calculation;
(4) strain perturbation. The ddk calculations has been included in order to
access to the piezoelectric tensor.

You can have a quick look to the output *tnlo_3.out*, when it is ready. It
contains rigid ions elastic and piezoelectric constants as well as the
internal strain coupling parameters. This information is also stored in a
database file (DDB) for further convenient analysis with ANADDB.

**Merge of the DDB.**

At this stage, all the relevant energy derivatives have been obtained and are
stored in individual databases. These must be combined with the MRGDDB merge
utility in order to get a unique database *tnlo_4.ddb.out*. Explicitely, you
should merge the files *tnlo_2o_DS4_DDB*, *tnlo_3o_DS4_DDB*, and *tnlo_2o_DS5_DDB*.
You might have a look at the input file for MRGDDB named *tnlo_4.in*, and use
it to perform the merge. You already used MRGDDB previously. It might be
located in *\$ABI_HOME/src/98_main* or another (build) directory. 
You might copy it, or make an alias.

**Analysis of the DDB.**

We are now ready for the analysis of the results using ANADDB. You can copy
the files *\$ABI_TESTS/tutorespfn/Input/tnlo_5.in* and
*\$ABI_TESTS/tutorespfn/Input/tnlo_5.files* in *Work-NLO*. You already used
ANADDB previously. It is located in the same directory as *abinit*.
You might copy it, or make an alias. The present input is in
principle very similar to the one you have used for the analysis of dynamical
and dielectric responses except that some new flags need to be activated.

For the strain perturbation you need to specify [[anaddb:elaflag]], [[anaddb:piezoflag]], and [[anaddb:instrflag]]:
    
            elaflag 3
            piezoflag  3
            instrflag  1

For the non-linear responses you need
    
            nlflag  1
            ramansr 1
            alphon  1
            prtmbm  1

[[anaddb:nlflag]]=1 activates the non-linear response.

[[anaddb:ramansr]] = 1 will impose the sum rule on the first-order change of the
electronic dielectric susceptibility under atomic displacement, hereafter
referred to as $\frac{d \chi}{d \tau}$. It is a condition of invariance of chi under
translation of the whole crystal, similar to the acoustic sum rules for
phonons at Gamma or the charge neutrality sum rule on Z*.

[[anaddb:prtmbm]]=1 will allow to get the mode by mode phonon contributions of
the ions to the electro-optic coefficients.

[[anaddb:alphon]]=1 will allow to get the previous mode by mode contribution
when aligning the eigenvectors with the cartesian axis of coordinates (in the
input, the principal axis should always be aligned with z for a convenient
analysis of the results).

Finally, the second list of phonon, specified with [[anaddb:nph2l]] and
[[anaddb:qph2l]], must also be explicitely considered to obtain the Raman
efficiencies of longitudinal modes (in a way similar to the computation of
frequencies of longitudinal mode frequencies at Gamma):
    
            # Wave vector list no. 2
            #***********************
                   nph2l  1
                   qph2l  1.0 0.0 0.0 0.0
    

You can now run the code ANADDB. The results are in the file tnlo_5.out.
Various interesting physical properties are now directly accessible in this
output in meaningful units. You can go through the file and look in order to
identify the results mention hereafter. Note that the order in which they are
given below is not the same than the order in which they appear in the
tnlo_5.out. You will have to jump between different sections of tnlo_5.out to find them.

For comparison, we report in parenthesis (...) the values obtained with ecut =
5, and for nonlinear responses in brackets [...] the fully converged result as
reported in [[cite:Veithen2005]].

  * Born effective charge of Al:
    
        Z*_Al = 2.043399E+00 (2.105999E+00)

  * Optical phonon frequencies at Gamma :
    
        w_TO (cm^-1) = 3.694366E+02 (3.602635E+02)
    w_LO (cm^-1) = 4.031189E+02 (3.931598E+02) 

  * Linear optical dielectric constant :
    
        Electronic dielectric tensor = 9.20199931 (9.94846084) 

  * Static dielectric constant :
    
        relaxed ion dielectric tensor = 10.95642097 (11.84823634) 

Some other quantities, as the piezoelectric coefficients, are related to the
strain response as it is more extensively discussed in the tutorial on the strain perturbation.

  * Proper piezoelectric coefficients :
    
        clamped ion (Unit:c/m^2) = -0.65029623 (-0.69401363)
    relaxed ion (Unit:c/m^2) =  0.03754602 (-0.04228777) 

Finally, different quantities are related to non-linear responses.

  * Nonlinear optical susceptibility :   
They are directly provided in the output in pm/V. As you can see the value
computed here is far from the well converged result as reported in 
[[cite:Veithen2005]].

    
        d_36 (pm/V)  = 21.175523 (32.772254) [fully converged :35] 

  * Electro-optic coefficients:   
As we asked for mode by mode decomposition the output provides individual
contributions. We report below a summary of the results. It concern the
clamped r_63 coefficient.

    
                Electronic EO constant (pm/V): -1.000298285 (-1.324507791) [-1.69]
            Full Ionic EO constant (pm/V):  0.543837671 (0.533097548)  [0.64]
            Total EO constant (pm/V):      -0.456460614 (-0.791410242) [-1.05] 

  * Raman properties   
The code directly report the Raman susceptibilities for both transverse (TO)
and longitudinal (LO) optic modes at Gamma:

    
                alpha(TO) = -0.008489212 (-0.009114814)
            alpha(LO) = -0.011466211 (-0.013439375) 
  
The basic quantity to get the Raman susceptibilities are the $\frac{d \chi}{d \tau}$ that
are also reported separately:
    
        dchi_23/dtau_1 (Bohr^-1) of Al = -0.094488281 (-0.099889084) 
  
In cubic semiconductors, it is usual to report the Raman polarizability of
optical phonon modes at Gamma which is defined as
    
        a = Omega_0 * dchi/dtau = Sqrt[mu * Omega_0] alpha 
  
where Omega_0 is the primitive unit cell volume (i.e. one quarter of the cubic
unit cell volume, to be expressed here in Ang) and mu is the reduced mass of
the system (1/mu = 1/m_Al + 1/m_As). From the previous data, we get :
    
                a(TO) (Unit: Ang^2)= -7.7233 (-8.4222112)  [-8.48]
            a(LO) (Unit: Ang^2)= -10.4317 (-12.418168) [-12.48] 

## 3 Finite difference calculation of the Raman tensor
  
For comparison with the DPFT calculation, we can compute $\frac{d \chi}{d \tau}$ for the Al
nucleus from finite differences. In practice, this is achieved by computing
the linear optical susceptibility for 3 different positions of the Al nucleus.
This is done with the file *\$ABI_TESTS/tutorespfn/Input/tnlo_6.in*, however
with the unrealistic cutoff of 2.8 Ha. The calculation is about 2 or 3 minutes
on a PC 3 GHz). For those who want to do it you anyway, you can copy
*\$ABI_TESTS/tutorespfn/Input/tnlo_6.in* in your working directory. If you
have time, you should modify the cutoff to [[ecut]] = 5 Ha, in order to obtain
realistic results. So, you might as well start the run after this modification
(the run is about two times more time-consuming than with 2.8 Ha).

You can have a look at this input file. It contains 8 datasets. We need to
compute the linear optical susceptibility (4 datasets for SC calculation in
the IBZ, NSC calculation in the full BZ, ddk, ddE) for different atomic
positions. We will do this for 2 sets of atomic positions, the reference
symmetric structure (referred to as tau=0), and a distorted structure
(referred to as tau= +0.01), for which the Al atom has been displaced to the
right by 0.01 Bohr (look at xcart to identify the differences). In the first
case, the dielectric tensor must be diagonal, isotropic, while in the second
case, a off-diagonal yz component will appear, that is an odd function of the
Al atomic displacement.

Supposing you are running the calculation, you have now time for a Belgian
Beer, why not a Gouyasse ?! ... Or you can look at the results as summarized below.

For tau = 0:
    
    Dielectric tensor, in cartesian coordinates,
         j1       j2             matrix element
      dir pert dir pert     real part    imaginary part
    
       1    4   1    4         9.2020015668        -0.0000000000
       1    4   2    4         0.0000000000        -0.0000000000
       1    4   3    4         0.0000000000        -0.0000000000
      
       2    4   1    4         0.0000000000        -0.0000000000
       2    4   2    4         9.2020015668        -0.0000000000
       2    4   3    4         0.0000000000        -0.0000000000
      
       3    4   1    4         0.0000000000        -0.0000000000
       3    4   2    4         0.0000000000        -0.0000000000
       3    4   3    4         9.2020015668        -0.0000000000
  
For tau = +0.01 :
    
    Dielectric tensor, in cartesian coordinates,
         j1       j2             matrix element
      dir pert dir pert     real part    imaginary part
    
       1    4   1    4         9.2023220436        -0.0000000000
       1    4   2    4        -0.0000000000        -0.0000000000
       1    4   3    4        -0.0000000000        -0.0000000000
      
       2    4   1    4        -0.0000000000        -0.0000000000
       2    4   2    4         9.2021443491        -0.0000000000
       2    4   3    4        -0.0123700617        -0.0000000000
      
       3    4   1    4        -0.0000000000        -0.0000000000
       3    4   2    4        -0.0123700617        -0.0000000000
       3    4   3    4         9.2021443491        -0.0000000000
  
Note that the following results would have been obtained for tau = -0.01 (with
obvious even / odd behaviour with respect to tau of the different components,
and some very small numerical noise):
    
    Dielectric tensor, in cartesian coordinates,
         j1       j2             matrix element
      dir pert dir pert     real part    imaginary part
    
       1    4   1    4         9.2023220610        -0.0000000000
       1    4   2    4        -0.0000000000        -0.0000000000
       1    4   3    4         0.0000000000        -0.0000000000
      
       2    4   1    4         0.0000000000        -0.0000000000
       2    4   2    4         9.2021443663        -0.0000000000
       2    4   3    4         0.0123700529        -0.0000000000
      
       3    4   1    4         0.0000000000        -0.0000000000
       3    4   2    4         0.0123700529        -0.0000000000
       3    4   3    4         9.2021443663        -0.0000000000
  
You can extract the value of dchi_23/dtau_1 for Al from the dielectric tensor
(hereafter called eps) above using the following finite-difference formula [unit of bohr^-1] :
    
     dchi_23/dtau_1= (1/4 pi) (eps_23[tau=+0.01] +eps_23[tau=0.00])/tau
                   = (1/4 pi) (-0.0123700 -0.0)/(0.01)
                   = -0.098437
    
This value is close to that obtained at [[ecut]]=5 from DFPT (-0.099889084).
When convergence is reached (beware, the k point convergence is extremely
slow, much slower than for other properties), both approaches allow to get the
right answer. You might therefore ask which approach is the most convenient
and should be used in practice.

As a guide, we can mention that the finite-difference approach give results
very similar to the DFPT ones for a similar cutoff and k-point grid. It is
however more tedious because, individual atomic displacement must be
successively considered (heavy for complex crystals) and the results must then
be converted into appropriate units with risk of error of manipulations.

The DFPT approach is the most convenient and avoid a lot of human work.
Everything is reported together (not only $\frac{d \chi}{d \tau}$ but also the full Raman
polarizability tensors) and in appropriate units. It should therefore be
considered as the best choice (when available, as in ABINIT).

## Calculation of the Raman Spectra
AFter an ANADDB calculation, one can visualize the Raman spectra using the post-processing script Raman_spec.py (The script can be found in the post-processing scripts repository ( ~/scripts/post_processing/)). Take a moment to explore the help menu (try Raman_spec.py --help) and maybe look at a typical input file (Raman_spec.py --input).  When you are done with that, execute the calculation:

   python Raman_spec.py "input file name"

On a normal computer, this calculation may take several minutes. 

This python script reads the output file generated by your ANADDB calculation, extracts the Raman tensor and phonon frequencies, and calculates the polarization dependent and powder-averaged Raman spectra.  All the calculated intensities (the 6 polarization dependent spectra and the powder-average spectra) are printed to a file.

We can view the Raman spectra with 

    xmgrace *_spec.out

The resulting powder-average spectra, plotted here with Gnuplot, is shown below. For the cubic structure calculated here, the resulting spectra contains a single Raman mode corresponding to an XY polarization. 

!!! tip

    ![](nlo_assets/ramanspec_tnlo5_spec.pdf)

A typical input file for the Raman_spec script contains the following variables
  * filename   - name of the ANADDB output file
  * outname    - uses specified output file name
  * temp       - temperature 
  * laser_freq - laser frequency
  * freq_unit  - output frequency (default cm$^{-1}$
  * spread     - spread of the Lorentzian (same for all modes)
  * n_freq     - number of output frequencies (default: 1000)
  * min_freq   - minimum frequency (default: 0.95 times the lowest active mode)
  * max_freq   - maximum frequency (default) 1.05 times the highest active mode)
  * relative_intensity - plots the relative intensity (if present)
  * keep_file  - keep the output file, any existing file is removed (if present)

Finally, if one includes a calculation of the frequency dependent dielectric tensor during the ANADDB calculation, then this program extracts that dielectric tensor and prints it to its own file. 

