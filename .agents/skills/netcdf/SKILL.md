---
name: netcdf
description: Diagnose ABINIT NetCDF content against its Fortran writer, including dimensions, slices, units, and Fortran/C index ordering. Use for .nc schema or value discrepancies.
---

Purpose

This skill teaches an AI agent how to inspect NetCDF files produced by scientific codes (e.g. ABINIT)
and compare their contents with the corresponding Fortran implementation.

The goal is to determine whether the values stored in the file agree with the in-memory Fortran arrays,
identify indexing mistakes, and understand how multidimensional arrays are serialized.

When to use

Use this skill whenever you need to

* debug NetCDF output
* compare Fortran arrays with data stored in a .nc file
* verify I/O routines
* understand how a variable is written
* investigate apparent transposition or permutation of indices

⸻

The most common pitfall

The NetCDF Fortran interface maps Fortran array order to the file's declared dimension order.
Tools using the C data model commonly display those dimensions in the reverse order from the
corresponding Fortran declaration.

Suppose the Fortran code writes

```fortran
real(dp) :: x(n1,n2,n3)
call nf90_put_var(ncid, varid, x)
```

The NetCDF variable may therefore appear as

(n3, n2, n1)

when inspected with

* Python (netCDF4)
* xarray
* ncdump
* h5dump

This does not mean the data were transposed incorrectly.
Confirm this from the actual dimension definition and writing call.
Do not reverse axes mechanically when the call uses `start`, `count`, a slice, or an intermediate buffer.

⸻

Never compare dimensions only

A frequent mistake is to compare

Fortran:
A(i,j,k)

Python:
A[k,j,i]

and conclude that the data are wrong.

Instead:

1. identify the Fortran declaration
2. identify the NetCDF dimensions
3. understand the ordering convention
4. compare corresponding physical indices

⸻

Always locate the writing routine
Before comparing numbers, find where the variable is written.

Typical workflow:

1. locate nf90_put_var
2. identify the Fortran array passed
3. inspect its declaration
4. inspect how it is filled
5. only then inspect the NetCDF file

Also check for unit conversion, packing, symmetry operations, and parallel decomposition before the
write.

Many apparent NetCDF bugs are actually caused by incorrect array construction before writing.

⸻

Beware of slicing

Many writes are of the form

call nf90_put_var(..., x(:,:,iq))

or

call nf90_put_var(..., x(:,iband))

The resulting variable may have fewer dimensions than the original array.

Always compare against the exact slice being written.

⸻

Check dimension names

Do not rely only on dimension sizes.

Instead, inspect

* dimension names
* variable attributes
* metadata
* units

For example

number_of_kpoints
number_of_spins
number_of_bands

is much more informative than

64
2
128

⸻

Compare values, not memory

When validating data,

* choose a few representative indices
* compute the expected value from the Fortran code
* compare with the NetCDF value

Avoid comparing raw memory layouts.

⸻

Recommended tools

Python

```python
from netCDF4 import Dataset
nc = Dataset("file.nc")
print(nc.variables.keys())
print(nc.variables["eigenvalues"].shape)
```

Command line

```bash
ncdump -h file.nc
```

or

```bash
ncdump -v variable file.nc
```
