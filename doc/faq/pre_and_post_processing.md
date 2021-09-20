---
authors: MG
---

# Pre-processing and post-processing FAQs

This page collects FAQs related to pre-processing and post-processing tools as well
as tricks to be more productive when preparing ABINIT input files.

## Can I use a xyx file to specify the unit cell?

Yes, see [[xyzfile]]

## Can I use a cif file to specify the unit cell?

ABINIT does not accept cif files in input but it is possible to convert cif to a format that ABINIT understands.
Perhaps, the easiest solution is to use the AbiPy |abistruct| script:

```
abistruct.py convert CIF_FILE
```

It is also possible to generate a cif file at the end of the GS calculation with the [[prtcif]] variable,
although it is much easier to use AbiPy to perform such task with the command line interface.

```
abistruct.py convert FILE -f cif
```

where FILE is **any file** with a structural info such as the *GSR.nc*, the **DDB** file, 
an ABINIT input file, etc.

!!! important

    In the cif format, the lattice is usually specified in terms of the three angles 
    formed by the direct lattice vectors ($\alpha, \beta, \gamma$ and their lengths).
    Note that these six parameters are not enough to specify the nine entries of 
    the lattice matrix [[rprimd]] expected by Abinit since any rigid rotation of the lattice
    will give the same angles and the lengths reported in the CIF file.
    The CIF parser implemented in pymatgen fixed this arbitrariness using some conventions.

    This can lead to unexpected behaviour if one tries to convert an Abinit structure to CIF format
    and then reuse this CIF file to reconstruct the structure since you are not guranteed to 
    get that same orientation of the lattice.
    Vectors or tensors components expressed in the reduced coordinates obviously depend 
    on the value of [[rprimd]] used in the calculation.
    Obviously one can always transform from one lattice to the other one but, in order to be on the safe side,
    we suggest using CIF files only when specifying the initial lattice and then rely on the ABINIT usual variables.

## Can I use POSCAR files with ABINIT?

Yes, please consult the documentation of the [[structure]] variable.
Note that it is also possible to convert a POSCAR to the Abinit format with:

```
abistruct.py convert POSCAR
```

## How can I convert a structure to primitive?

Use abipy :

```
abistruct.py primitive FILE
```

## Wyckoff positions and site symmetries

Use abipy :

```
abistruct.py wyckoff FILE
```

## How can I compare my structure with the materials project database?

Use abipy :

```
abistruct.py mp_match FILE
```

## How can I visualize a structure?

With abitk, one can read a DEN/WFK file and produce an xsf file

```
abistruct.py visualize FILE -a vesta
```

```
abistruct.py visualize --help

  ...

  -a APPNAME, --appname APPNAME
                        Application name. Possible options: avogadro, ovito, v_sim, vesta, xcrysden, mpl
                        (matplotlib), mayavi, vtk
```


## Can I include external files in the input?

Sure you can.
Abinit supports the `include` statement that can be used to include external files. In your input file,
add a line with the syntax:

```
include "name_of_file_to_be_included_in_the_input_file"
```

where you must obviously change the argument name_of_file... but keep the quotation marks.

!!! note

    **include** is not a standard variable thus it not listed in the variables page.


## Is there an easy way to print the warnings/comments in the log file?

Warnings, Comments and Errors are written in Yaml format.

```
abiopen.py LOG_FILE --print   # or -p
```

## Is there an easy way to obtain a high-symmetry k-path?

Use abipy :

```
abistruct.py kpath FILE
```

## Is there an easy way to plot an electronic band structure?

The most important results of a GS run, including the KS energies, 
are stored in the GSR.nc output file.
To plot the band structure with matplotlib, use the |abiopen| script and the syntax:

```
abiopen.py GSR_FILE --expose  # or -e
```

Other options are available. See `abiopen.py --help`.

To generate an xmgrace file with electron bands, use |abiview| and the command:

```
abiview.py ebands out_GSR.nc --xmgrace
```

## How can I visualize the evolution of the cell during a structural relaxation run.

The lattice paramenters, atomic positions, forces at each relaxation step are stored in the 
HIST.nc output file.

To print

```
abiopen.py HIST_FILE -p
```

To visualize the results with matplotlib, use:

```
abiopen.py HIST_FILE -e
```

To visualize the results with matplotlib, use:

```
abiview.py hist out_HIST.nc -a ovito     ==> Visualize relaxation/molecular dynamics results with ovito.
abiview.py hist out_HIST.nc --xdatcar    ==> Convert HIST.nc into XDATCAR format (caveat: assume fixed unit
```
