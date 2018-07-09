#!/usr/bin/env bash
# Copyright (C) 1998-2017 ABINIT group (XG)
# 
# The purpose of this script is to change the copyright year
# in nearly all files in the ABINIT package. 
# First you should reapply the present script without changing it, because it has been seen that some people bring
# routines with erroneous date during a few months after the beginning of a new one ... !
#
# Then one should update the present script to put the current year (at present, valid for going from 2017 to 2018 !) !
#
# Then should be called from the top directory  (here a list of generic filenames, ordered on the basis of the alphanumeric string)
# developers/maintainers/change_year.sh *.ac */*.ac */*/*.am */*/*.bat */*/*/*.bat */*/*.c */*/*/*.c */*/*.cf */*/*.cnf */*/*.com */*/*.conf */*/*.cu */*/*.csh 
# developers/maintainers/change_year.sh */*/*.dep */*/*.dir */*.env */*/*.finc */*/*.f90 */*/*.F90 */*/*/*.F90 *.in */*.in */*/*.in 
# developers/maintainers/change_year.sh */*/*.h */*/*.help */*/*/*.help */*/*.html */*/*/*/html */*/*/*.log */*/*.m */*/*/*.m */*/make* 
# developers/maintainers/change_year.sh */*/*.mk */*/*.m4 */*/*/*.m4 */*/Makefile */*/*/*.out */*/*.pl */*/*/*.pl */README */*/README 
# developers/maintainers/change_year.sh */*/*.sav */*.sh */*/*.src */*/*/*.stdout */*/*.tex */*/*/*.tex */*/*.txt */*/*_ */*/*/*.yml

# Please do not change the permission of py files. Not all py modules must be executable! So, the following command should be used simply to see whether the copyright date has to be changed ... and if not,
# please, restart from the previous version, and if yes, do it by hand !
# developers/maintainers/change_year.sh */*/*.py */*/*/*.py 
#
# In the previous list, files without an extension are not treated (except the Makefile and README files - warning some README are only links ...), 
# and */*/*.sh are not treated (except tests/*/*.sh), because of conflict with the present script file extension !!
# Also config/scripts/abilint cannot be treated automatically ...
#
# So, also issue, one after the other (cut and paste the following):
# developers/maintainers/change_year.sh autom4te.cache/tr* config/scripts/a* config/scripts/clean* config/scripts/u* config/wrappers/wrap-fc
# developers/maintainers/change_year.sh developers/bzr_helpers/abinit-forge-branch developers/bzr_helpers/bzr-make-patch 
# developers/maintainers/change_year.sh developers/maintainers/change2.sh developers/maintainers/change.sh developers/various/change_perl.sh developers/various/fixed_to_free tests/cpu/Refs/changeref 
# developers/maintainers/change_year.sh developers/various/*.sh developers/various/fixed_to_free doc/config/scripts/make* doc/manpages/abinit.1 fallbacks/config/scripts/make* INSTALL 
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
 sed -e 's&Copyright (c)&Copyright (C)&' $file > tmp.yrup
 sed -e 's&(C) 1987-2017 ABINIT&(C) 1987-2018 ABINIT&' tmp.yrup > tmp.yr87
 sed -e 's&(C) 1991-2017 ABINIT&(C) 1991-2018 ABINIT&' tmp.yr87 > tmp.yr91
 sed -e 's&(C) 1992-2017 ABINIT&(C) 1992-2018 ABINIT&' tmp.yr91 > tmp.yr92
 sed -e 's&(C) 1993-2017 ABINIT&(C) 1993-2018 ABINIT&' tmp.yr92 > tmp.yr93
 sed -e 's&(C) 1996-2017 ABINIT&(C) 1996-2018 ABINIT&' tmp.yr93 > tmp.yr96
 sed -e 's&(C) 1997-2017 ABINIT&(C) 1997-2018 ABINIT&' tmp.yr96 > tmp.yr97
 sed -e 's&(C) 1998-2017 ABINIT&(C) 1998-2018 ABINIT&' tmp.yr97 > tmp.yr98
 sed -e 's&(C) 1999-2017 ABINIT&(C) 1999-2018 ABINIT&' tmp.yr98 > tmp.yr99
 sed -e 's&(C) 2000-2017 ABINIT&(C) 2000-2018 ABINIT&' tmp.yr99 > tmp.yr00
 sed -e 's&(C) 2001-2017 ABINIT&(C) 2001-2018 ABINIT&' tmp.yr00 > tmp.yr01
 sed -e 's&(C) 2002-2017 ABINIT&(C) 2002-2018 ABINIT&' tmp.yr01 > tmp.yr02
 sed -e 's&(C) 2003-2017 ABINIT&(C) 2003-2018 ABINIT&' tmp.yr02 > tmp.yr03
 sed -e 's&(C) 2004-2017 ABINIT&(C) 2004-2018 ABINIT&' tmp.yr03 > tmp.yr04
 sed -e 's&(C) 2005-2017 ABINIT&(C) 2005-2018 ABINIT&' tmp.yr04 > tmp.yr05
 sed -e 's&(C) 2006-2017 ABINIT&(C) 2006-2018 ABINIT&' tmp.yr05 > tmp.yr06
 sed -e 's&(C) 2007-2017 ABINIT&(C) 2007-2018 ABINIT&' tmp.yr06 > tmp.yr07
 sed -e 's&(C) 2008-2017 ABINIT&(C) 2008-2018 ABINIT&' tmp.yr07 > tmp.yr08
 sed -e 's&(C) 2009-2017 ABINIT&(C) 2009-2018 ABINIT&' tmp.yr08 > tmp.yr09
 sed -e 's&(C) 2010-2017 ABINIT&(C) 2010-2018 ABINIT&' tmp.yr09 > tmp.yr10
 sed -e 's&(C) 2011-2017 ABINIT&(C) 2011-2018 ABINIT&' tmp.yr10 > tmp.yr11
 sed -e 's&(C) 2012-2017 ABINIT&(C) 2012-2018 ABINIT&' tmp.yr11 > tmp.yr12
 sed -e 's&(C) 2013-2017 ABINIT&(C) 2013-2018 ABINIT&' tmp.yr12 > tmp.yr13
 sed -e 's&(C) 2014-2017 ABINIT&(C) 2014-2018 ABINIT&' tmp.yr13 > tmp.yr14
 sed -e 's&(C) 2015-2017 ABINIT&(C) 2015-2018 ABINIT&' tmp.yr14 > tmp.yr15
 sed -e 's&(C) 2016-2017 ABINIT&(C) 2016-2018 ABINIT&' tmp.yr15 > tmp.yr16
#The next lines are both needed, as some developers decide to use one, and some the other ...
 sed -e 's&(C) 2017-2017 ABINIT&(C) 2017-2018 ABINIT&' tmp.yr16 > tmp.yr17
 sed -e 's&(C) 2017 ABINIT&(C) 2017-2018 ABINIT&' tmp.yr17 > tmp.yr
 echo "changes done "
 # put the modified file at the correct place
 mv tmp.yr $file
 echo "file $file written "
done
rm -f tmp.yr*  
#chmod 755 */*/*.sh */*/*.py */*/*.pl */*/*.com config/*/make* developers/*/make*  */config/scripts/* */*.sh
# Please do not change the permission of py files. Not all py modules must be executable!
chmod 755 */*/*.sh */*/*.pl config/*/make* */config/scripts/* */*.sh
chmod 755 config/scripts/* developers/bzr_helpers/* developers/*/* tests/cpu/Refs/changeref
