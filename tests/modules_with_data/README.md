
This directory contains git submodules with precomputed netcdf or text files.
These files are used in selected ABINIT tutorials to bypass the preliminary steps (e.g. expensive GS/DFPT calculation).
Note that we use netcdf or text files as these formats are portable across different architectures 
contrarily to Fortran binary files that are not portable due to hardware-endianess.

The advantage of such git-module-based approach is that in the tutorials we can focus 
on the most important part of the calculation without having to waste time to generate the required files.
This improves the stability of the automatic tests as well as the user experience.
The disadvantage is that including all these files in the official distribution would 
increase the size of tarball file.
This is the reason, why we decide to host these external files on github repos.

Users interested in running tutorials based on git submodules can easily download the data from the internet
or use git directly to fetch the repository from the internet.

The ABINIT test suite is submodule-aware ...

    git submodule update --remote --init

[submodule "tests/modules_with_data/MgO_eph_zpr"]
	path = tests/modules_with_data/MgO_eph_zpr
	url = https://github.com/abinit/MgO_eph_zpr.git
	branch = master

[submodule "tests/modules_with_data/MgO_eph_zpr"]
	path = tests/modules_with_data/MgO_eph_zpr
	url = git://github.com/abinit/MgO_eph_zpr.git
	branch = master

To add a new submodule, use:

    git submodule add https://github.com/abinit/MgO_eph_zpr.git
    git submodule add https://github.com/abinit/MgO_eph_zpr.git

[submodule "tests/modules_with_data/MgO_eph_zpr"]
	path = tests/modules_with_data/MgO_eph_zpr
	url = git://github.com/abinit/MgO_eph_zpr.git
	branch = master
