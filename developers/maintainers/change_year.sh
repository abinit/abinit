#!/usr/bin/env bash
# Copyright (C) 1998-2020 ABINIT group (XG)
# 
# The purpose of this script is to change the copyright year
# in nearly all files in the ABINIT package. 
# First you should reapply the present script without changing it, because it has been seen that some people bring
# routines with erroneous date during a few months after the beginning of a new one ... !
#
# Then one should update the present script to put the current year (at present, valid for going from 2019 to 2020 !) !
#
# Then should be called from the top directory  (here a list of generic filenames, ordered on the basis of the alphanumeric string)
# developers/maintainers/change_year.sh *.ac */*.ac */*/*.am */*/*.bat */*/*/*.bat */*/*.c */*/*/*.c */*/*.cf */*/*.cnf */*/*.com */*/*.conf */*/*.cu */*/*.csh 
# developers/maintainers/change_year.sh */*/*.dep */*/*.dir */*.env */*/*.finc */*/*.f90 */*/*.F90 */*/*/*.F90 *.in */*.in */*/*.in 
# developers/maintainers/change_year.sh */*/*.h */*/*.help */*/*/*.help */*/*.html */*/*/*/html */*/*/*.log */*/*.m */*/*/*.m */*/make* 
# developers/maintainers/change_year.sh */*/*.mk */*/*.m4 */*/*/*.m4 */*/Makefile */*/*/*.out */*/*.pl */*/*/*.pl bindings/README fallbacks/README tests/README */*/README 
# developers/maintainers/change_year.sh */*/*.sav */*.sh */*/*.src */*/*/*.stdout */*/*.tex */*/*/*.tex */*/*.txt */*/*_ 

# Please do not change the permission of py files. Not all py modules must be executable! 
# Simply issue the in-place command
# sed -i 's!2019 ABINIT!2020 ABINIT!' */*.py */*/*.py */*/*/*.py
#
# In the previous list, files without an extension are not treated (except the Makefile and README files - warning some README are only links ...), 
# and */*/*.sh are not treated (except tests/*/*.sh), because of conflict with the present script file extension !!
# Also config/scripts/abilint cannot be treated automatically ...
#
# So, also issue, one after the other (cut and paste the following):
# developers/maintainers/change_year.sh config/scripts/a* config/scripts/clean* config/scripts/u* 
# developers/maintainers/change_year.sh developers/maintainers/change2.sh developers/maintainers/change.sh developers/various/fixed_to_free tests/cpu/Refs/changeref 
# developers/maintainers/change_year.sh developers/maintainers/suppress.sh scripts/*/*.sh
# developers/maintainers/change_year.sh developers/various/*.sh developers/various/fixed_to_free doc/config/scripts/make* fallbacks/config/scripts/make* INSTALL 
# developers/maintainers/change_year.sh tests/config/scripts/make-makefiles-tests tests/cpu/Refs/changeref scripts/configure/upgrade-build-config packages/debian/copyright 
# 
# Moreover, one should complement the present script with a search 
# grep 'past_year ABINIT' * */* */*/* */*/*/* */*/*/*/*
# and treat by hand the remaining files ...
#
#XG 100118 Still other problems with copyrights might be detected by using the following command (replace 2018 by the present year !):
# grep -i opyright * */* */*/* */*/*/* */*/*/*/* | grep -v 2018 | grep -v '!! COPYRIGHT' | grep -v 'Oldenburg' | grep -v 'Stefan Goedecker' | grep -v 'doc/rel' | grep -v 'Remove' | grep -v 'tests/' | grep -v 'EXC group' | grep -v 'PWSCF group' | grep -v 'Makefile' | grep -v 'abinit.d' | grep -v 'fallbacks' | grep -v 'doc/features/features' | grep -v 'doc/install_notes/install' | grep -v 'COPYING' | grep -v 'gui' | grep -v 'default' | grep -v js_files

for file in "$@"
do
 echo "working on $file"
 rm -f tmp.yr*  
#The primitive command was
#sed -i 's&(C) 2018-2019 ABINIT&(C) 2018-2020 ABINIT&' $file
#for different values replacing 2018 
#The following uses a regexp. This change the two last digits of the year.
 sed -i 's&\(([cC]) ....-20\)19\( ABINIT\)&\120\2&' $file
#The next line is also needed, as some developers decide to use this syntax, and some the other ...
 sed -i 's&\(([cC]) 20\)19\( ABINIT\)&\119-2020\2&' $file
 echo "file $file treated "
done
#chmod 755 */*/*.sh */*/*.py */*/*.pl */*/*.com config/*/make* developers/*/make*  */config/scripts/* */*.sh
# Please do not change the permission of py files. Not all py modules must be executable!
chmod 755 */*/*.sh */*/*.pl config/*/make* */config/scripts/* */*.sh
chmod 755 config/scripts/* developers/*/* tests/cpu/Refs/changeref
