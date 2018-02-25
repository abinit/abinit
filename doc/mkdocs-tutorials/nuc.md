---
authors: JWZ, XG
---

# Lesson on properties at the nuclei  

## Observables near the atomic nuclei.  

The purpose of this lesson is to show how to compute several observables of
interest in Moessbauer, NMR, and NQR spectroscopy, namely:

  * the electric field gradient, 
  * the isomer shift (not yet available), 
  * and the electronic density itself (not yet available). 

This lesson should take about 1 hour.


## 1 Electric field gradient

  
Various spectroscopies, including nuclear magnetic resonance and nuclear
quadrupole resonance (NMR and NQR), as well as Moessbauer spectroscopy, show
spectral features arising from the electric field gradient at the nuclear
sites. Note that the electric field gradient (EFG) considered here arises from
the distribution of charge within the solid, not due to any external electric
fields.

The way that the EFG is observed in spectroscopic experiments is through its
coupling to the nuclear electric quadrupole moment. The physics of this
coupling is described in various texts, for example Principles of Magnetic
Resonance, 3rd ed., C. P. Slichter (Springer, New York, 1989). ABINIT computes
the field gradient at each site, and then reports the gradient and its
coupling based on input values of the nuclear quadrupole moments.

The electric field and its gradient at each nuclear site arises from the
distribution of charge, both electronic and ionic, in the solid. The gradient
especially is quite sensitive to the details of the distribution at short
range, and so it is necessary to use the PAW formalism to compute the gradient
accurately. The various sources of charge in the PAW decomposition are
summarized in the following equation:  
![](../documents/lesson_nuc/charge.gif)  

  
Here the "v" subscript indicates valence, "c" indicates core, and "Z"
indicates the ions. Essentially the gradient must be computed for each source
of charge, which is done in the code as follows:

  * Valence space described by planewaves: expression for gradient is Fourier-transformed at each nuclear site. 
  * Ion cores: gradient is computed by an Ewald sum method 
  * On-site PAW contributions: moments of densities are integrated in real space around each atom, weighted by the gradient operator 
The code reports each contribution separately if requested.

The electric field gradient computation is performed at the end of a ground-
state calculation, and takes almost no additional time. The tutorial file is
for stishovite, a polymorph of SiO2. In addition to typical ground state
variables, only two additional variables are added:

    
    
            prtefg  2
            quadmom 0.0 -0.02558
    

The first variable instructs Abinit to compute and print the electric field
gradient, and the second gives the quadrupole moments of the nuclei, one for
each type of atom. Here we are considering silicon and oxygen, and in
particular Si-29, which as zero quadrupole moment, and O-17, the only stable
isotope of oxygen with a non-zero quadrupole moment.

After running the file tnuc_1.in through abinit, you can find the following
near the end of the output file:

    
    
     Electric Field Gradient Calculation 
    
     Atom   1, typat   1: Cq =      0.000000 MHz     eta =      0.000000
    
          efg eigval :     -0.165960
    -         eigvec :     -0.000001    -0.000001    -1.000000
          efg eigval :     -0.042510
    -         eigvec :      0.707107    -0.707107     0.000000
          efg eigval :      0.208470
    -         eigvec :      0.707107     0.707107    -0.000002
    
          total efg :      0.082980     0.125490    -0.000000
          total efg :      0.125490     0.082980    -0.000000
          total efg :     -0.000000    -0.000000    -0.165960
    

This fragment gives the gradient at the first atom, which was silicon. Note
that the gradient is not zero, but the coupling is---that's because the
quadrupole moment of Si-29 is zero, so although there's a gradient there's
nothing in the nucleus for it to couple to.

Atom 2 is an oxygen atom, and its entry in the output is:

    
    
     Atom   2, typat   2: Cq =      6.603688 MHz     eta =      0.140953
    
          efg eigval :     -1.098710
    -         eigvec :     -0.707107     0.707107     0.000000
          efg eigval :      0.471922
    -         eigvec :     -0.000270    -0.000270     1.000000
          efg eigval :      0.626789
    -         eigvec :      0.707107     0.707107     0.000382
    
          total efg :     -0.235961     0.862750     0.000042
          total efg :      0.862750    -0.235961     0.000042
          total efg :      0.000042     0.000042     0.471922
    
    
          efg_el :     -0.044260    -0.065290     0.000042
          efg_el :     -0.065290    -0.044260     0.000042
          efg_el :      0.000042     0.000042     0.088520
    
          efg_ion :     -0.017255     0.306132    -0.000000
          efg_ion :      0.306132    -0.017255    -0.000000
          efg_ion :     -0.000000    -0.000000     0.034509
    
          efg_paw :     -0.174446     0.621908     0.000000
          efg_paw :      0.621908    -0.174446     0.000000
          efg_paw :      0.000000     0.000000     0.348892
    

Now we see the electric field gradient coupling, in frequency units, along
with the asymmetry of the coupling tensor, and, finally, the three
contributions to the total. Note that the valence part, efg_el, is quite
small, while the ionic part and the on-site PAW part are larger. In fact, the
PAW part is largest--this is why these calculations give very poor results
with norm-conserving pseudopotentials, and need the full accuracy of PAW.



