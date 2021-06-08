---
authors: MG
---

# Ground-state FAQs

This page collects FAQs related to the GS part of Abinit.

## What are the dimensions defining the workload of a typical GS run?

[[nkpt]]
[[nband]]
[[nspinor]] and [[nsppol]]
The cutoff energy [[ecut]] defines the maximum number of planewaves [[mpw]] 
[[dilatmx]]
[[istwfk]]

The FFT mesh [[ngfft]]] is comuted from [[boxcutmin]] and 

[[natom]]

These variables are reported at the beginning of the main output file.

Please keep in mind that there is a strong interplay between these dimensions and the MPI algorithm

## Which GS algorithm do you recommend for relatively large systems?

[[paral_kgb]] == 1 is the most efficient and scalable eigenvalue solver

## Is it a good idea to set the number of bands to the highest occupied state to accelerate a GS calculation?

The answer is **No** and the explanation requires some technical discussion about the KS solver.
In principle, the workload and the memory footprint is proportional to [[nband]] so, at first sight,
reducing the number of states seems the most efficient choice.

In the case of non-magnetic semiconductors ([[nsppol]] == 1), one can in principle set [[nband]]
to the highest occupied band index i.e. the total number of valence electrons divided by 2 (4 for Silicon).
On the other hand, the GS solver may require less iteration including empty states in the GS run
usually helps to reach convergence


!!! important

    In the case of metals, one needs partially occupied states to account for the tail of the occupation
    function [[occoopt]] and the broadening parameter [[tsmear]] so this is a **very bad idea**

## The SCF cycle does not converge

The default SCF algorithm for SCF (see [[iscf]]) is an excellent compromise between speed and robustness.
It is referred to as Pulay II in Fig. 10 of
N. D. Woods and M. C. Payne and P. J. Hasnip, J. Phys.: Condens. Matter 31, 453001 (2019),
and appears on the Pareto frontier.
Some parameters in the Kerker preconditioning have also an effect on the balance between speed and robustness,
and the default values are good for metallic systems.
As thouroughly discussed in this reference, there is no fool-proof fast algorithm.

If your SCF cycle do not converge (see nres2 or vres2 column in the output file), the reasons can be:
(1) Insufficient underlying accuracy in the solution of the Schrödinger at fixed potential
(2) Transient non-linear behaviour of the SCF, due either to
    (2a) sudden change of occupation numbers (usually only for metals), or
    (2b) long-wavelength fluctuations of potential, bigger than the gap.

Start first to address (2), by some tuning which can come without significantly slowing ABINIT.
Test different values of [[diemac]].
The default is very large, and suitable fo metals.
Try running with a small value (perhaps 5) if your system is a bulk semiconductor or a molecule.
If it is a doped solid, some intermediate value (perhaps 50) might be appropriate.

If this does not work, try to address (1), set tolrde to 0.001, increase the value of [[nline]]
(e.g. 6 or 8 instead of the default 4), as well as of [[nnsclo]]  (e.g. 2 instead of the default).
Your residm should be significantly lower than without such modifications, and perhaps your problem will be solved.

If this still does not work, but your residm did not look bad after all before addressing (1),
then revert back to the default values of [[tolrde]], [[nline]] and [[nnsclo]], as this indeed slows down ABINIT.

Then, try to use [[iscf]] 2.
This is a very slow but inconditionally convergent algorithm provided [[diemix]] is small enough.
At some stage of convergence, you might restart with the obtained DEN a better algorithm,
as non-linear effects should have been eliminated at some stage.


## How can I visualize the evolution of the cell during a structural relaxation run.

```
abiopen.py HIST_FILE -e

abiview.py hist out_HIST.nc -a ovito     ==> Visualize relaxation/molecular dynamics results with ovito.
abiview.py hist out_HIST.nc --xdatcar    ==> Convert HIST.nc into XDATCAR format (caveat: assume fixed unit
```

```
abiopen.py HIST_FILE -p
```
