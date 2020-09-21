
This directory contains git submodules with precomputed (netcdf of txt) files that can be used
in the ABINIT tutorials to bypass the preliminary steps (e.g. expensive GS/DFPT calculation).
We use netcdf or txt as these formats are portable across different architectures contrarily to Fortran binary files.

The advantage of such approach is that in the tutorial we can focus on the most important part of the calculation
without having to waste time to generate the required files.
This improves the stability of the test as well as user experience.

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


[submodule "tests/modules_with_data/MgO_eph_zpr"]
	path = tests/modules_with_data/MgO_eph_zpr
	url = git://github.com/abinit/MgO_eph_zpr.git
	branch = master
