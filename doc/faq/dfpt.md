---
authors: MG, MJV
---

# DFPT FAQs

This page collects FAQs related to the DFPT part of Abinit.

## Can I use paral_kgb 1 in the DFPT part?

No, the DFPT code uses a completely different parallel algorithm in which
the workload and memory is automatically distributed at runtime depending on the
number of MPI processes and the perburbation to be computed. The parallelism is over the
k points and the bands, so whatever the size of the system, one can use several dozen processors efficiently.

In principle it is also possible to parallelize over the perturbations
(see [[paral_rf]] and [[nppert]]) but keep in mind that this approach is not optimal
from the point of view of the worload distribution as each perturbation has its own list of irreducible k-points
(see discussion below).
Moreover if one perturbation does not converge withing [[nstep]] iterations,
ABINIT will abort and we may end up with zero DDB files produced!
Finally, it is also possible to run different d points concurrently (with different ABINIT runs), and to merge their DDB files.
The workflows with Abipy can do the launch of many concurrent jobs for you.

## How to get the irreducible q-points for phonon calculations

When doing a phonon calculation, only the irreducible perturbations need to be calculated,
and abinit and anaddb will complete them by symmetry.
In particular, the dynamical matrix will be completed and symmetrized at each q-point.

To get the irreducible q-points for your system, you can do a small auxiliary calculation
(with low [[ecut]] and minimal [[nband]]).
Make a new directory, and copy your input file for the ground state there.
If you set the k-point grid [[ngkpt]] to be equal to the q-point grid you want for the phonons, and set [[kptopt]] to 1.
Also, don't forget to set [[shiftk]] to 0 0 0, in order to get a non-shifted grid with the Gamma point.
From the header of the log or output file, the [[kpt]] points can be copied directly to the phonon input file, as q-points.

One can also use this one, but one should change the output format

    abistruct.py abikmesh si.cif --ngkpt 2 2 2 --shiftk 0 0 0

## warning The dynamical matrix was incomplete

> I am getting the following warning in phonon calculations:
  The dynamical matrix was incomplete : phonon frequencies may be wrong
  Does anyone knows what is wrong and how can I fix it?

This warning is given by abinit but is not a problem: when you calculate certain perturbations
at a given q-point abinit tries to immediately complete and diagonalize the dynamical matrix.
In general this is not possible, as you have not done a full set of irreducible perturbations
all at once (hence the warning).
The frequencies are output, but many may be zero because certain dynamical matrices have not been calculated yet.

You need to continue the calculations and use mrgddb to merge the DDB files, and then anaddb
to do the reconstruction of the dynamical matrix.
Follow out the whole of the rf tutorials and you will see the phonon frequencies will come out when anaddb is run.

## Why some perturbations are faster to compute than others?

Each perturbation breaks the initial symmetry of the crystal thus only a subset of the
crystalline symmetries can be used to define the irreducible set of wavevectors used in the DFPT equations.
For instance, phonon calculations at $\Gamma$ are much faster that calculations done at non-zero $\qq$-points
since more symmetries can be exploited.
Note, however, that even perturbations with the same $\qq$-point may have a different workload.
For best performance, each perturbation should be in principle computed in a different run and with a different
number of MPI procs.

## What to do if a perturbation does not converge?

First of all, try to restart from the first-order WFK file using [[ird1wf]], [[get1wf]].
In principle, it is also possible to restart from the first-order density via [[get1den]]
but use this approach only if the first-order WFK file is not available.
You may want to use [[prtwf]] = -1 in the DFPT part to produce the first-order WFK file only
when the DFPT SCF cycle does not converge in order to reduce the amount of data written to disk.

You may also try to increase [[nline]].

## Is it a good idea to compute all the perturbations in parallel with a single input file?

If you care about perfomance, the answer is **definitely NO!**
As already mentioned, each perturbation has a different workload so it is much better
the compute the different perturbations with separate input file and use more MPI processes
for the heaviest perturbations.

## Is there an easy way to compute the phonon band structure from the DDB file?

```
abiview.py ddb DDB_FILE
```

## How can I reduce the breaking of the ASR?

## Is there an easy way to visualize the breaking of the ASR?

```
abiview.py ddb_asr DDB_FILE
```

## Can I convert the DDB into phonopy format?

Yes, but you will need to use AbiPy to convert from the two formats.

