---
authors: JWZ, XG
---

# Tutorial on properties at the nuclei

## Observables near the atomic nuclei.

The purpose of this tutorial is to show how to compute several observables of
interest in M&ouml;ssbauer, NMR, and NQR spectroscopy, namely:

  * the electric field gradient,
  * the isomer shift

This tutorial should take about 1 hour.

[TUTORIAL_README]

## Electric field gradient

Various spectroscopies, including nuclear magnetic resonance and nuclear
quadrupole resonance (NMR and NQR), as well as M&ouml;ssbauer spectroscopy, show
spectral features arising from the electric field gradient at the nuclear
sites. Note that the electric field gradient (EFG) considered here arises from
the distribution of charge within the solid, not due to any external electric fields.

The way that the EFG is observed in spectroscopic experiments is through its
coupling to the nuclear electric quadrupole moment. The physics of this
coupling is described in various texts, for example [[cite:Slichter1978]].
Abinit computes the field gradient at each site, and then reports the gradient and its
coupling based on input values of the nuclear quadrupole moments.

The electric field and its gradient at each nuclear site arises from the
distribution of charge, both electronic and ionic, in the solid. The gradient
especially is quite sensitive to the details of the distribution at short
range, and so it is necessary to use the PAW formalism to compute the gradient
accurately. For charge density $n({\mathbf r})$, the potential $V$ is given
by

$$ V({\mathbf r})=\int \frac{n({\mathbf r'})}{ |\mathbf{r}-\mathbf{r'}| } d{\mathbf r'} $$

and the electric field gradient is

$$V_{ij} = -\frac{\partial^2}{\partial x_i\partial x_j}V({\mathbf r}).$$

The gradient is computed at each nuclear site, for each source of charge arising
from the PAW decomposition (see [the tutorial PAW1](paw1) ).
This is done in the code as follows:

  * Valence space described by planewaves: expression for gradient is Fourier-transformed at each nuclear site.
  * Ion cores: gradient is computed by an Ewald sum method
  * On-site PAW contributions: moments of densities are integrated in real space around each atom, weighted by the gradient operator

The code reports each contribution separately if requested.

The electric field gradient computation is performed at the end of a ground-state calculation,
and takes almost no additional time.
The tutorial file is for stishovite, a polymorph of SiO$_2$. In addition to typical ground state
variables, only two additional variables are added:

    prtefg  2
    quadmom 0.0 -0.02558

{% dialog tests/tutorial/Input/tnuc_1.in %}

The first variable instructs Abinit to compute and print the electric field
gradient, and the second gives the quadrupole moments of the nuclei, in
barns, one for each type of atom.
A standard source for quadrupole moments is [[cite:Pyykko2008]].
Here we are considering silicon and oxygen, and in
particular Si-29, which has zero quadrupole moment, and O-17, the only stable
isotope of oxygen with a non-zero quadrupole moment.

After running the file *tnuc_1.in* through Abinit, you can find the following
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
PAW part is largest -- this is why these calculations give very poor results
with norm-conserving pseudopotentials, and need the full accuracy of PAW.
Experimentally, the nuclear quadrupole coupling for O-17 in stishovite is
reported as $6.5\pm 0.1$ MHz, with asymmetry $0.125\pm 0.05$ [[cite:Xianyuxue1994]].

## Fermi contact interaction

The Fermi contact interaction arises from overlap of the electronic wavefunctions
with the atomic nucleus, and is an observable for example in
M&ouml;ssbauer spectroscopy [[cite:Greenwood1971]]. In M&ouml;ssbauer spectra,
the isomer shift $\delta$ is expressed in velocity units as

$$
\begin{equation}
\label{eq:mossbauershift}
\delta = \frac{c}{E_\gamma}\frac{2\pi Z e^2}{3}(|\Psi(R)_A|^2-|\Psi(R)_S|^2)\Delta\langle r^2\rangle ,
\end{equation}
$$

where $\Psi(R)$ is the electronic
wavefunction at nuclear site $R$, for the absorber (A) and source (S);
$c$ is the speed of light, $E_\gamma$ is the nuclear transition energy, and $Z$ the atomic number;
and $\Delta\langle r^2\rangle$ the change in the nuclear size squared. All these quantities
are assumed known in the M&ouml;ssbauer spectrum of interest, except $|\Psi(R)|^2$, the
Fermi contact term.

Abinit computes the Fermi contact term in the PAW formalism by using as observable
$\delta(R)$, that is, the Dirac delta function at the nuclear site [[cite:Zwanziger2009]].
Like the electric field gradient computation, the Fermi contact calculation
is performed at the end of a ground-
state calculation, and takes almost no time. There is a tutorial file for
SnO$_2$, which, like stishovite studied above, has the rutile structure.
In addition to typical ground state
variables, only one additional variable is needed:

    prtfc  1

{% dialog tests/tutorial/Input/tnuc_2.in %}

After running this file, inspect the output and look for the phrase
"Fermi-contact Term Calculation". There you'll find the FC output for
each atom; in this case, the Sn atoms, [[typat]] 1, yield a contact term
of 72.2969 atomic units (charge per volume $e/a^3_0$).

To interpret M&ouml;ssbauer spectra you need really both a source and
an absorber; in the tutorial we provide also a file for $\alpha$-Sn (grey
tin, which is non-metallic).

{% dialog tests/tutorial/Input/tnuc_3.in %}

If you run this file, you should find a contact term of 102.3008.

To check your results, you can use experimental data for the isomer shift $\delta$
for known compounds to compute $\Delta\langle r^2\rangle$ in Eq.\ref{eq:mossbauershift}
(see [[cite:Zwanziger2009]]). Using our results above together with standard
tin M&ouml;ssbauer parameters of $E_\gamma = 23.875$ keV and an experimental shift
of 2.2 mm/sec for $\alpha$-Sn relative to SnO$_2$, we find
$\Delta\langle r^2\rangle = 5.74\times 10^{-3}\mathrm{fm}^2$, in decent agreement
with other calculations of 6--7$\times 10^{-3}\mathrm{fm}^2$ [[cite:Svane1987]], [[cite:Svane1997]].
