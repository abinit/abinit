!{\src2tex{textfont=tt}}
!!****f* ABINIT/ptgmadata
!! NAME
!! ptgmadata
!!
!! FUNCTION
!! Return magnetic point group symbol from the magnetic point group number
!! The symbols and numbers are taken from  The Internationl Tables for Crystallography
!! Volume A, 1983 Ed. Theo Hahn, D. Reidel Publishing Company and
!! The mathematical theory of symmetry in solids, Representation theory for point
!! groups and space groups, 1972, C.J. Bradley and A.P.
!! Cracknell, Clarendon Press, Oxford.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! ptgroupma = space group number
!!
!! OUTPUT
!! ptgrpmasb= symbol
!!
!! NOTES
!!
!! PARENTS
!!      prtspgroup
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ptgmadata(ptgroupma,ptgrpmasb)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ptgmadata'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ptgroupma
 character(len=10),intent(out) :: ptgrpmasb

!Local variables-------------------------------
!character(len=500) :: message

! *************************************************************************

 select case (ptgroupma)
 case(1)
   ptgrpmasb="-1'"
 case(2)
   ptgrpmasb="2'"
 case(3)
   ptgrpmasb="m'"
 case(4)
   ptgrpmasb="2/m'"
 case(5)
   ptgrpmasb="2'/m"
 case(6)
   ptgrpmasb="2'/m'"
 case(7)
   ptgrpmasb="2'2'2"
 case(8)
   ptgrpmasb="m'm'2"
 case(9)
   ptgrpmasb="m'm2'"
 case(10)
   ptgrpmasb="m'm'm'"
 case(11)
   ptgrpmasb="mmm'"
 case(12)
   ptgrpmasb="m'm'm"
 case(13)
   ptgrpmasb="4'"
 case(14)
   ptgrpmasb="-4'"
 case(15)
   ptgrpmasb="42'2'"
 case(16)
   ptgrpmasb="4'22'"
 case(17)
   ptgrpmasb="4/m'"
 case(18)
   ptgrpmasb="4'/m'"
 case(19)
   ptgrpmasb="4'/m"
 case(20)
   ptgrpmasb="4m'm'"
 case(21)
   ptgrpmasb="4'mm'"
 case(22)
   ptgrpmasb="-42'm'"
 case(23)
   ptgrpmasb="-4'2m'"
 case(24)
   ptgrpmasb="-4'm2'"
 case(25)
   ptgrpmasb="4/m'm'm'"
 case(26)
   ptgrpmasb="4/m'mm"
 case(27)
   ptgrpmasb="4'/mmm'"
 case(28)
   ptgrpmasb="4'/m'm'm"
 case(29)
   ptgrpmasb="4/mm'm'"
 case(30)
   ptgrpmasb="32'"
 case(31)
   ptgrpmasb="3m'"
 case(32)
   ptgrpmasb="-6'"
 case(33)
   ptgrpmasb="-6m'2'"
 case(34)
   ptgrpmasb="-6'm2'"
 case(35)
   ptgrpmasb="-6'm'2"
 case(36)
   ptgrpmasb="6'"
 case(37)
   ptgrpmasb="-3'"
 case(38)
   ptgrpmasb="-3m'"
 case(39)
   ptgrpmasb="-3'm"
 case(40)
   ptgrpmasb="-3'm'"
 case(41)
   ptgrpmasb="62'2'"
 case(42)
   ptgrpmasb="6'2'2"
 case(43)
   ptgrpmasb="6/m'"
 case(44)
   ptgrpmasb="6'/m'"
 case(45)
   ptgrpmasb="6'/m"
 case(46)
   ptgrpmasb="6m'm'"
 case(47)
   ptgrpmasb="6'm'm"
 case(48)
   ptgrpmasb="6'/mmm'"
 case(49)
   ptgrpmasb="6'/m'm'm"
 case(50)
   ptgrpmasb="6/m'm'm'"
 case(51)
   ptgrpmasb="6/m'mm"
 case(52)
   ptgrpmasb="6/mm'm'"
 case(53)
   ptgrpmasb="m'3"
 case(54)
   ptgrpmasb="-4'3m'"
 case(55)
   ptgrpmasb="4'32'"
 case(56)
   ptgrpmasb="m'3m'"
 case(57)
   ptgrpmasb="m'3m"
 case(58)
   ptgrpmasb="m3m'"
 end select

end subroutine ptgmadata
!!***
