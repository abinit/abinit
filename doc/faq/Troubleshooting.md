## Troubleshooting

Something going wrong ? Better to realize this than to think everything is OK ...

In this page, the user will find tips to fix several problems with ABINIT.   

Many other tips might be found on the [[https://discourse.abinit.org|forum]], or in the [[https://docs.abinit.org|standard ABINIT documentation]].


### Configure problems

    config.status: error: cannot find input file: `tests/Makefile.in'

This problem might come from the unpacking of the tar.gz. It was reported that using the Krusader file manager does not unpack correctly the tar.gz, so some files are missing. A manual unpack (tar -xvf ...) should solve the problem.


    configure script won't accept hostnames containing a hyphen

Also apparently an error due to incorrect unpacking of the tar.gz . 


### ABINIT does not run 

(to be filled ...)
### ABINIT stops without finishing its tasks 

Make sure you **get and read the error messages**. They are usually written at the end of the log file, and look like:

<code>
--- !ERROR
src_file: m_kg.F90
src_line: 179
mpi_rank: 1
message:   (here follows the error message)
</code>

If you do not see such error message at the end of the log file, the error message might be found alternatively in an "error" file, depending on the command that you used to launch ABINIT.
Also, if you are running in batch, the behaviour might be different, and error messages might end up in another file than the log file.
If you run in parallel, the error messages might have been written in still another, separate file, entitled 
<code>
__ABI_MPIABORTFILE__
</code>

If ABINIT stops quickly, consider to **run in sequential, interactively** (not in parallel neither in batch) to ease the understanding of what is happening.

Sometimes, the error message is prepared by another message giving preliminary information. Please, be sure to have identified the relevant information in the log file.

If you think to have obtained all information that ABINIT was supposed to give you, try to identify whether ABINIT stops because of geometry optimization problem or SCF problem. In the first case, your input geometry is perhaps crazy. Anyhow, see the next items in this troubleshooting page.

By the way, there might also be memory problems ... See <https://docs.abinit.org/topics/TuningSpeedMem>.


### Incorrect initial geometry 

Many mistakes done by beginners are related to incorrect starting geometry.

Here is a check list :

  * Check that the **units** are correct for your cell parameters and atomic positions. __ABINIT uses atomic unit by default__, but can use several other units if specified by the user, see [[https://docs.abinit.org/guide/abinit/#parameters|ABINIT input parameters]].
  * Check that the [[https://docs.abinit.org/variables/basic/#typat|typat]] atom types are correct with respect to [[https://docs.abinit.org/variables/basic/#xred|xred]] or [[https://docs.abinit.org/variables/basic/#xcart|xcart]]
  * Check that the **number of atoms** [[https://docs.abinit.org/variables/basic/#natom|natom]] is coherent with your list of coordinates, [[https://docs.abinit.org/variables/basic/#xred|xred]] or [[https://docs.abinit.org/variables/basic/#xcart|xcart]]. __ABINIT reads only the coordinates of [[https://docs.abinit.org/variables/basic/#natom|natom]] nuclei, and **ignore others**__.
  * Relax first the atomic positions **at fixed primitive vectors** before optimizing the cell. Explicitly, use a first datadet with [[https://docs.abinit.org/variables/rlx/#optcell|optcell]] = 0, then a second dataset with non-zero [[https://docs.abinit.org/variables/rlx/#optcell|optcell]], in which you tell ABINIT to read optimized atomic positions using [[https://docs.abinit.org/variables/rlx/#getxred|getxred]] or [[https://docs.abinit.org/variables/rlx/#getxcart|getxcart]]. In this second dataset, do not forget to use[[https://docs.abinit.org/variables/rlx/#dilatmx|dilatmx]] bigger than 1 if you expect the volume of the cell to increase during the optimization. Possibly after the atomic position relaxation, make a run with [[https://docs.abinit.org/variables/rlx/#chkdilatmx|chkdilatmx]]=0, then a third run with [[https://docs.abinit.org/variables/rlx/#dilatmx|dilatmx]]=1. See the additional suggestions in the documentation of [[https://docs.abinit.org/variables/rlx/#optcell|optcell]] . 
  * ABINIT can read VASP POSCAR external files containing unit cell parameters and atomic positions. See the input variable [[https://docs.abinit.org/variables/basic/#structure|structure]]. This might help in setting the geometry correctly.
  * Try to visualize your primitive cell. There are numerous possibilities, using Abipy, VESTA, XCrysDen (to be documented).


### SCF cycle does not converge

The default ABINIT algorithm for SCF (see input variable iscf) is an excellent compromise between speed and robustness.
It is referred to as Pulay II in Fig. 10 of
N. D. Woods and M. C. Payne and P. J. Hasnip, J. Phys.: Condens. Matter 31, 453001 (2019), and appear on the Pareto frontier.
Some parameters in the Kerker preconditioning have also an effect on the balance between speed and robustness, and the ABINIT default values
are good for metallic systems. As thouroughly discussed in this reference, there is no fool-proof fast algorithm.

If your SCF cycle do not converge (see nres2 or vres2 column in the output file), the reasons can be:
(1) Insufficient underlying accuracy in the solution of the Schrödinger at fixed potential
(2) Transient non-linear behaviour of the SCF, due either to (2a) sudden change of occupation numbers (usually only for metals), or (2b) long-wavelength fluctuations of potential, bigger than the gap.

Start first to address (2), by some tuning which can come without significantly slowing ABINIT. Test different values of diemac. The default is very large, and suitable fo metals. Try running with a small value (perhaps 5) if your system is a bulk semiconductor or a molecule. If it is a doped solid, some intermediate value (perhaps 50) might be appropriate.

If this does not work, try to address (1), set tolrde to 0.001, increase the value of nline (e.g. 6 or 8 instead of the default 4), as well as of nnsclo (e.g. 2 instead of the default). Your residm should be significantly lower than without such modifications, and perhaps your problem will be solved.

If this still does not work, but your residm did not look bad after all before addressing (1), then revert back to the default values of tolrde, nline an nnsclo, as this indeed slows down ABINIT.

Then, try to use iscf 2. This is a very slow but inconditionally convergent algorithm provided diemix is small enough. At some stage of convergence, you might restart with the obtained _DEN a better algorithm, as non-linear effects should have been eliminated at some stage.


### Geometry optimization does not converge 

(to be filled)

### ABINIT delivers strange results

(to be filled)

