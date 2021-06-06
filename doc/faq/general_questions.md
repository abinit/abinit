---
authors: MG
---

# General questions

This page gives a beginner's introduction to the ABINIT resources,
the package, and the main ABINIT applications.

## How can I print the Abinit version?

```
abinit --version
```

See ``abinit --help`` for additional command line options.

## Is there any option to check the validity of the input without running a full calculation?

```
abinit --dry-run run.abi
```

## What if ABINIT stops without finishing its tasks?

Make sure you get and read the error messages.
They are usually written at the end of the log file, and look like:

```yaml
--- !ERROR
src_file: m_kg.F90
src_line: 179
mpi_rank: 1
message:   (here follows the error message)
...
```

If you do not see such error message at the end of the log file, the error message
might be found alternatively in an “error” file, depending on the command that you used to launch ABINIT.
Also, if you are running in batch, the behaviour might be different, and error messages
might end up in another file than the log file.
If you run in parallel, the error messages might have been written in still another, separate file, entitled

__ABI_MPIABORTFILE__

If ABINIT stops quickly, consider to run in sequential, interactively (not in parallel neither in batch)
to ease the understanding of what is happening.
Sometimes, the error message is prepared by another message giving preliminary information.
Please, be sure to have identified the relevant information in the log file.

If you think to have obtained all information that ABINIT was supposed to give you, try to identify whether
ABINIT stops because of geometry optimization problem or SCF problem.
In the first case, your input geometry is perhaps crazy.
Anyhow, see the next items in this troubleshooting page.

## Is the memory estimate provided by Abinit at the beginning of the run reliable?

Well, **it depends**.
The value reported by the code in the case of ground-state or DFPT calculations should
represent a reasonble estimate
Very likely, it is smaller than the real value because not all workspace arrays are taken into account.

Unfortunately, it is not easy to estimate the memory requirements for GW, BSE, EPH calculations
so do not rely on the values printed by the code at the beginning.
To have a better idea of the memory footprint, we suggest to use the log file:

```
grep "<<< MEM" run.log -A 3
```

## Is the paral_kgb = 1 algorithm supported everywhere?

No, [[paral_kgb]] = 1 refers to a particular MPI-distribution of the wavefunctions that is only available
for GS calculations (total energy calculations, relaxations, molecular dynamics as well
as NSCF band structure calculations).
Variables such as [[npkpt]], [[npband]], [[npfft]], [[npspinor]] are relevant only if [[paral_kgb]] = 1.
All the other Abinit drivers (e.g. DFPT, GW, BSE, EPH) use a completely different approach to parallelize
the calculation and distribute memory.

## Where can I find "old pseudos"?

Legacy pseudopotential files are available
[here](https://www.abinit.org/sites/default/files/PrevAtomicData/psp-links/pseudopotentials.html)

## Can I mix NC-pseudos with PAW?

No, the two implementations are mutually exclusive.

## Can I use UPF pseudos with Abinit?

Pseudos in UPF2 format are not supported.
Note that it is not just a matter of format.
There are indeed fundamental differences between the formalism implemented in
Abinit and the one used by QE, especially for PAW and NC-pseudos with SOC.
The pseudodojo project provides NC-pseudos

## Can I run calculations with an XC functional different from the one used to generate the pseudos?

Yes, use [[ixc]] to change the XC functional in the input file.

## Is there any tool to print a pseudo in a more human-readable format?

```
abiopen.py PSEUDO_FILE -p
```

## chkprim and unit cell multiplicity

In highly symmetric crystals you may end up with the following error message from abinit:

```
chkprimit : ERROR -
According to the symmetry finder, the unit cell is
NOT primitive. The multiplicity is 4 .
The use of non-primitive unit cells is allowed
only when the input variable chkprim is 0.
Action : either change your unit cell (rprim or angdeg),
or set chkprim to 0.

leave_new : decision taken to exit ...
```

By default abinit checks that the unit cell is primitive and contains the smallest possible number of atoms. 
For example, instead of the conventional cubic FCC unit cell with 4 atoms

```
rprim
1 0 0
0 1 0
0 0 1
natom 4
xred
0 0 0
0 0.5 0.5
0.5 0 0.5
0.5 0.5 0
```

you should be using the following unit cell vectors:

```
rprim
0.0 0.5 0.5
0.5 0.0 0.5
0.5 0.5 0.0
natom 1
xred
0 0 0
```

to get a primitive unit cell with 3 axes separated by angles of 60 degrees.
Another possibility is:

```
angdeg 3*60.
natom 1
xred
0 0 0
```

Using the primitive cell ensures the fastest calculation and the best use of symmetry operations, 
so in general you should listen to abinit and reduce your unit cell.

If you have a good reason to do so, you can override abinit and force it to use a non-primitive unit cell, 
by setting [[chkprim]] 1 in the input file. 
This will, however, disable symmetry operations in certain places.

One possible reason to do this is, for example, making a large supercell of a crystal (say 3x3x3 primitive unit cells) 
in which you want to introduce a defect. Doing the pristine crystal calculation in the 3x3x3 supercell is possible, 
but not useful (you will just get 27 times the energy). Once you have introduced the defect, of course, 
you will lower the supercell's symmetry and abinit will no longer complain that the cell is not primitive. 
chkprim 1 is no longer to be used for the defected cell, in which you want all possible symmetries to be exploited once again.

## pspxc from pseudopotential not equal to ixc

Abinit complains if the exchange correlation functional (variable [[ixc]]) 
used in the input is not the same as that specified in the pseudopotential / atomic data files used.

```
 pspatm: WARNING -
  Pseudopotential file pspxc=       7,
  not equal to input ixc=       1.
  These parameters must agree to get the same xc 
  in ABINIT code as in psp construction.
  Action : check psp design or input file.
  Assume experienced user. Execution will continue.
```

This is usually dangerous, as you are making uncontrollable errors in compensating for the Vxc of the core 
in the pseudopotential with a different functional.

However, there is an important exception: all the LDA functionals (ixc 1-7) are basically identical, 
but with different functional forms to fit the same data. 
As a result, mixing these ixc is mostly harmless (as above). 
The warning appears often, as the many pseudopotentials on the web site are created with a variety of ixc, 
not necessarily the default value 1.
