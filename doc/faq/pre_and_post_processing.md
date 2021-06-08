---
authors: MG
---

# Pre-processing and post-processing FAQs

This page collects FAQs related to pre-processing, post-processing tools as well
as tricks to be more productive.

## Can I use a xyx file to specify the unit cell?

Yes, see [[xyzfile]]

## Can I use a cif file to specify the unit cell?

ABINIT does not accept cif files in input but it is possible to convert cif
to a format that ABINIT can understand.
Perhaps, the easiest solution is to use the AbiPy |abistruct| script and the syntax:

```
abistruct.py convert CIF_FILE
```

[[structure]]

It is also possible to generate a cif file at the end of the GS calculation with the [[prtcif]] variable,
although it is much easier to use AbiPy to perform such task with:

```
abistruct.py convert FILE -f cif
```

where FILE is **any file** with a structural info such as the *GSR.nc*, the **DDB** file etc.

## Can I use POSCAR files with ABINIT?

Yes, see the documentaion of the [[structure]] variable.
|abistruct|

## How can I convert a structure to primitive ?

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
Abinit supports the `include` statement that can be used to include external files with the
syntax:

```
include
```

!!! note

    **include** is not a stadard variable thus it not listed in the variables page.


## Is there an easy way to print the warnings/comments in the log file?

Warnings, Comments and Errors are written in Yaml format.

```
abiopen.py LOG_FILE --print   # or -p 
```

## Is there an easy way to plot the band structure?

```
abiopen.py GSR_FILE --expose  # or -e
```
