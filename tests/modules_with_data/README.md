
This directory contains git submodules providing precomputed files that are then used 
in the ABINIT tutorials to bypass the preliminary steps of the computation (e.g. GS/DFPT calculation).

The advantage of such approach is that in the tutorials we can focus on the most important part of the calculation
without having to waste time to prepare the calculation. Moreover, it is possible to provide input ingredients that 
much closer to convergence and the preliminary steps need not to be executed.

The disadvantage is that including these files in the official package would increase the size of tarball file.
This is the reason, why these external files are hosted in git repos on github in order not to increase the size of the package.

Users interested in running tutorials based on git submodules can easily download the data from the internet or use git 
directly to fetch the repository.

The ABINIT test suite is submodule-aware ...

git submodule update --remote --init
