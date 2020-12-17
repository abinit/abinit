
This directory contains git submodules with precomputed netcdf or text files.
These files are used in selected ABINIT tutorials to bypass the preliminary steps (e.g. expensive GS/DFPT calculation).
We use netcdf or text files as these formats are portable across different architectures 
contrarily to Fortran binary files that are not portable due to hardware-endianess.

The advantage of such git-module-based approach is that, in the tutorials, we can focus 
on the most important part of the calculation without having to waste time to generate the required files.
Moreover this improves the stability of the automatic tests as well as user experience.
The disadvantage is that including all these files in the official distribution would 
increase the size of tarball file.
This is the reason, why we decided to host these external files on github repos 
that are included here as submodules. These submodules are not included in the tarball file.
(see ~abinit/config/specs/junk.conf).

Users interested in running tutorials based on git submodules can easily download the data from the internet
or use git directly to fetch the repository from the internet.
For a quick intro to submodules, please consult the [official documentation](https://git-scm.com/book/en/v2/Git-Tools-Submodules).

Also, note that the ABINIT test suite is submodule-aware ...
A test that requires a submodule, must declared it in the files section with the syntax.

```
#%% [files]
#%% use_git_submodule = MgO_eph_zpr
```

and files can be passed to ABINIT via strings, e.g.:

```
getden_filepath "MgO_eph_zpr/flow_zpr_mgo/w0/t0/outdata/out_DEN.nc"

structure "abifile:MgO_eph_zpr/flow_zpr_mgo/w0/t0/outdata/out_DEN.nc"
```

See e.g. tutorespfn/Input/teph4zpr_3.in


    git submodule update --remote --init

[submodule "tests/modules_with_data/MgO_eph_zpr"]
	path = tests/modules_with_data/MgO_eph_zpr
	url = https://github.com/abinit/MgO_eph_zpr.git
	branch = master

To add a new submodule, use:

    git submodule add https://github.com/abinit/MgO_eph_zpr.git

[submodule "tests/modules_with_data/MgO_eph_zpr"]
	path = tests/modules_with_data/MgO_eph_zpr
	url = git://github.com/abinit/MgO_eph_zpr.git
	branch = master
