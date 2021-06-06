---
authors: MG
---

# Pre-processing and post-processing FAQs

This page collects FAQs related to pre-processing, post-processing tools as well
as tricks to be more productive.

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

## Can I include external files in the input?

Sure you can. 
Abinit supports the `include` statement that can be used to

```
include
```

!!! note

    **include** is not a stadard variable thus it not listed in the variables page.


## Is there an easy way to print the warnings/comments in the log file?

```
abiopen.py LOG_FILE -p
```

## Is there an easy way to plot the band structure?

```
abiopen.py GSR_FILE -e
```
