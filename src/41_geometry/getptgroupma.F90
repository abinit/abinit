!{\src2tex{textfont=tt}}
!!****f* ABINIT/getptgroupma
!! NAME
!! getptgroupma
!!
!! FUNCTION
!! Return magnetic point group number from the full point group number
!! and the point group number of the non-magnetic symmetry operations.
!! The (normal) point group numbers are taken from
!! The International Tables for Crystallography
!! Volume A, 1983 Ed. Theo Hahn, D. Reidel Publishing Company
!! The magnetic point group number are taken from
!! The mathematical theory of symmetry in solids, Representation theory for point
!! groups and space groups, 1972, C.J. Bradley and A.P.
!! Cracknell, Clarendon Press, Oxford.
!! In particular, see table 7.1 of the latter reference
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! ptgroup = character(len=5) point group of all the symmetry operation
!! ptgroupha = character(len=5) point group of the non-magnetic symmetry operation (halved point group)
!!
!! OUTPUT
!! ptgroupma = magnetic point group number
!!
!! NOTES
!!
!! PARENTS
!!      symanal
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine getptgroupma(ptgroup,ptgroupha,ptgroupma)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getptgroupma'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: ptgroupma
 character(len=5),intent(in) :: ptgroup,ptgroupha

! *************************************************************************

!DEBUG
!write(std_out,*)' getptgroupma : enter '
!write(std_out,*)' ptgroup="',ptgroup,'"'
!write(std_out,*)' ptgroupha="',ptgroupha,'"'
!ENDDEBUG

 ptgroupma=0
 select case (ptgroup)
 case("   -1")
   ptgroupma=1
 case("    2")
   ptgroupma=2
 case("   -2")
   ptgroupma=3
 case("  2/m")
   if(ptgroupha=="    2")ptgroupma=4
   if(ptgroupha=="   -2")ptgroupma=5
   if(ptgroupha=="   -1")ptgroupma=6
 case("  222")
   ptgroupma=7
 case("  mm2")
   if(ptgroupha=="    2")ptgroupma=8
   if(ptgroupha=="   -2")ptgroupma=9
 case("  mmm")
   if(ptgroupha=="  222")ptgroupma=10
   if(ptgroupha=="  mm2")ptgroupma=11
   if(ptgroupha=="  2/m")ptgroupma=12
 case("    4")
   ptgroupma=13
 case("   -4")
   ptgroupma=14
 case("  422")
   if(ptgroupha=="    4")ptgroupma=15
   if(ptgroupha=="  222")ptgroupma=16
 case("  4/m")
   if(ptgroupha=="    4")ptgroupma=17
   if(ptgroupha=="   -4")ptgroupma=18
   if(ptgroupha=="  2/m")ptgroupma=19
 case("  4mm")
   if(ptgroupha=="    4")ptgroupma=20
   if(ptgroupha=="  mm2")ptgroupma=21
 case(" -42m")
   if(ptgroupha=="   -4")ptgroupma=22
   if(ptgroupha=="  222")ptgroupma=23
   if(ptgroupha=="  mm2")ptgroupma=24
 case("4/mmm")
   if(ptgroupha=="  422")ptgroupma=25
   if(ptgroupha=="  4mm")ptgroupma=26
   if(ptgroupha=="  mmm")ptgroupma=27
   if(ptgroupha==" -42m")ptgroupma=28
   if(ptgroupha=="  4/m")ptgroupma=29
 case("   32")
   ptgroupma=30
 case("   3m")
   ptgroupma=31
 case("   -6")
   ptgroupma=32
 case(" -62m")
   if(ptgroupha=="   -6")ptgroupma=33
   if(ptgroupha=="   3m")ptgroupma=34
   if(ptgroupha=="   32")ptgroupma=35
 case("    6")
   ptgroupma=36
 case("   -3")
   ptgroupma=37
 case("  -3m")
   if(ptgroupha=="   -3")ptgroupma=38
   if(ptgroupha=="   3m")ptgroupma=39
   if(ptgroupha=="   32")ptgroupma=40
 case("  622")
   if(ptgroupha=="    6")ptgroupma=41
   if(ptgroupha=="   32")ptgroupma=42
 case("  6/m")
   if(ptgroupha=="    6")ptgroupma=43
   if(ptgroupha=="   -3")ptgroupma=44
   if(ptgroupha=="   -6")ptgroupma=45
 case("  6mm")
   if(ptgroupha=="    6")ptgroupma=46
   if(ptgroupha=="   3m")ptgroupma=47
 case("6/mmm")
   if(ptgroupha==" -62m")ptgroupma=48
   if(ptgroupha=="  -3m")ptgroupma=49
   if(ptgroupha=="  622")ptgroupma=50
   if(ptgroupha=="  6mm")ptgroupma=51
   if(ptgroupha=="  6/m")ptgroupma=52
 case("  m-3")
   ptgroupma=53
 case(" -43m")
   ptgroupma=54
 case("  432")
   ptgroupma=55
 case(" m-3m")
   if(ptgroupha=="  432")ptgroupma=56
   if(ptgroupha==" -43m")ptgroupma=57
   if(ptgroupha=="  m-3")ptgroupma=58
 end select

!DEBUG
!write(std_out,*)' getptgroupma : exit '
!write(std_out,*)' ptgroupma="',ptgroupma,'"'
!ENDDEBUG

end subroutine getptgroupma
!!***
