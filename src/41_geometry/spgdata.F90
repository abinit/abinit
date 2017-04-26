!{\src2tex{textfont=tt}}
!!****f* ABINIT/spgdata
!! NAME
!! spgdata
!!
!! FUNCTION
!! Return point and space group data : Bravais lattice symbol,
!! international symbol, Schonflies symbol, multiplicity
!! The symbols are taken from  The International Tables for Crystallography
!! Volume A, 1983 Ed. Theo Hahn, D. Reidel Publishing Company and
!! The mathematical theory of symmetry in solids, Representation theory for point
!! groups and space groups, 1972, C.J. Bradley and A.P.
!! Cracknell, Clarendon Press, Oxford.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2017 ABINIT group (RC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! spgroup = space group number
!! spgorig = space group origin
!! spgaxor = space group axis orientation
!!
!! OUTPUT
!! brvsb=Bravais lattice symbol (P, I, F, A, B, C, R)
!! intsb=international symbol (like m3m, 222, 2_12_12_1)
!! intsbl=international symbol in long format like P2_b = P121)
!! ptintsb=International point group symbol
!! ptschsb=Schoenflies point group symbol
!! sporder=multiplicity of the space group
!! schsb=Schoenflies symbol
!!
!! NOTES
!! brvsb, intsb, and schsb have been extensively checked, while
!! more checking should be done for the others
!! XG20160612 : in particular, at present it might be that spgaxor and spgorig are indetermined (e.g. spgaxor=-1;spgorig=-1) at input.
!! When this has a bearing on some of the output variables (even brvsb or intsb !), these are mentioned as being X, unknown, or to be determined.
!!
!! PARENTS
!!      m_ab7_symmetry,prt_cif,prtspgroup,symsgcube,symsghexa,symsgmono
!!      symsgortho,symsgtetra,symspgr
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine spgdata(brvsb,intsb,intsbl,ptintsb,ptschsb,&
& schsb,spgaxor,spgroup,sporder,spgorig)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spgdata'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spgaxor,spgorig,spgroup
 integer,intent(out) :: sporder
 character(len=1),intent(out) :: brvsb
 character(len=15),intent(out) :: intsb,ptintsb,ptschsb,schsb
 character(len=35),intent(out) :: intsbl

!Local variables-------------------------------
!character(len=500) :: message

! *************************************************************************

 intsbl="same"
!defaults for case spgroup is not well defined (eg chkprim 0)
 brvsb="P"
 intsb="1"
 schsb="C1^1"
 sporder=1

 select case (spgroup)
 case(1)
   brvsb="P"; intsb="1"; schsb="C1^1"; sporder=1
 case(2)
   brvsb="P"; intsb="-1"; schsb="Ci^1"; sporder=2
 case(3)
   brvsb="P"; intsb="2"; schsb="C2^1"; sporder=2
   select case (spgaxor)
   case(1)
     intsbl="P 2 _b = P 1 2 1"
   case(2)
     intsbl="P 2_a = P 2 1 1"
   case(3)
     intsbl="P 2 _c = P 1 1 2"
   case default
     intsbl="intsbl to be determined"
   end select
 case(4)
   brvsb="P"; intsb="2_1"; schsb="C2^2"; sporder=2
   select case (spgaxor)
   case(1)
     intsbl="P 2 1 _b = P 1 2_1 1"
   case(2)
     intsbl="P 2 1 _a = P 2_1 1 1"
   case(3)
     intsbl="P 2 1 _c = P 1 1 2_1"
   case default
     intsbl="intsbl to be determined"
   end select
 case(5)
   brvsb="C"; intsb="2"; schsb="C2^3"; sporder=2
   select case (spgaxor)
   case(1)
     intsbl="C 2 _b1 =  C 1 2 1"
   case(2)
     intsbl="C 2 _a1 =  B 2 1 1"
   case(3)
     intsbl="C 2 _a2 =  C 2 1 1"
   case(4)
     intsbl="C 2 _a3 =  I 2 1 1"
   case(5)
     intsbl="C 2 _b2 =  A 1 2 1"
   case(6)
     intsbl="C 2 _b3 =  I 1 2 1"
   case(7)
     intsbl="C 2 _c1 =  A 1 1 2"
   case(8)
     intsbl="C 2 _c2 =  B 1 1 2 = B 2"
   case(9)
     intsbl="C 2 _c3 =  I 1 1 2"
   case default
     intsbl="intsbl to be determined"
   end select
 case(6)
   brvsb="P"; intsb="m"; schsb="Cs^1"; sporder=2
   select case (spgaxor)
   case(1)
     intsbl="P m _b = P 1 m 1"
   case(2)
     intsbl="P m _a = P m 1 1"
   case(3)
     intsbl="P m _c = P 1 1 m"
   case default
     intsbl="intsbl to be determined"
   end select
 case(7)
   brvsb="P"; intsb="c"; schsb="Cs^2"; sporder=2
   select case (spgaxor)
   case(1)
     intsbl="P c _b1 = P 1 c 1"
   case(2)
     intsbl="P c _a1 = P b 1 1"
   case(3)
     intsbl="P c _a2 = P n 1 1"
   case(4)
     intsbl="P c _a3 = P c 1 1"
   case(5)
     intsbl="P c _b2 = P 1 n 1"
   case(6)
     intsbl="P c _b3 = P 1 a 1"
   case(7)
     intsbl="P c _c1 = P 1 1 a"
   case(8)
     intsbl="P c _c2 = P 1 1 n"
   case(9)
     intsbl="P c _c3 = P 1 1 b = P b"
   case default
     intsbl="intsbl to be determined"
   end select
 case(8)
   brvsb="C"; intsb="m"; schsb="Cs^3"; sporder=4
   select case (spgaxor)
   case(1)
     intsbl="C m _b1 = C 1 m 1"
   case(2)
     intsbl="C m _a1 = B m 1 1"
   case(3)
     intsbl="C m _a2 = C m 1 1"
   case(4)
     intsbl="C m _a3 = I m 1 1"
   case(5)
     intsbl="C m _b2 = A 1 m 1"
   case(6)
     intsbl="C m _b3 = I 1 m 1"
   case(7)
     intsbl="C m _c1 = A 1 1 m"
   case(8)
     intsbl="C m _c2 = B 1 1 m = B m"
   case(9)
     intsbl="C m _c3 = I 1 1 m"
   case default
     intsbl="intsbl to be determined"
   end select
 case(9)
   brvsb="C"; intsb="c"; schsb="Cs^4"; sporder=4
   select case (spgaxor)
   case(1)
     intsbl="C c _b1 = C 1 c 1"
   case(2)
     intsbl="C c _a1 = B b 1 1"
   case(3)
     intsbl="C c _a2 = C n 1 1"
   case(4)
     intsbl="C c _a3 = I c 1 1"
   case(5)
     intsbl="C c _b2 = A 1 n 1"
   case(6)
     intsbl="C c _b3 = I 1 a 1"
   case(7)
     intsbl="C c _c1 = A 1 1 a"
   case(8)
     intsbl="C c _c2 = B 1 1 n"
   case(9)
     intsbl="C c _c3 = I 1 1 b"
   case default
     intsbl="intsbl to be determined"
   end select
 case(10)
   brvsb="P"; intsb="2/m"; schsb="C2h^1"; sporder=4
   select case (spgaxor)
   case(1)
     intsbl="P 2/m _b = P 1 2/m 1"
   case(2)
     intsbl="P 2/m _a = P 2/m 1 1"
   case(3)
     intsbl="P 2/m _c = P 1 1 2/m"
   case default
     intsbl="intsbl to be determined"
   end select
 case(11)
   brvsb="P"
   intsb="2_1/m"
   schsb="C2h^2"
   sporder=4
   select case (spgaxor)
   case(1)
     intsbl="P 2_1/m _b = P 1 2_1/m 1"
   case(2)
     intsbl="P 2_1/m _a = P 2_1/m 1 1"
   case(3)
     intsbl="P 2_1/m _c = P 1 1 2_1/m"
   case default
     intsbl="intsbl to be determined"
   end select
 case(12)
   brvsb="C"; intsb="2/m"; schsb="C2h^3"; sporder=8
   select case (spgaxor)
   case(1)
     intsbl="C 2/m _b1 = C 1 2/m 1"
   case(2)
     intsbl="C 2/m _a1 = B 2/m 1 1"
   case(3)
     intsbl="C 2/m _a2 = C 2/m 1 1"
   case(4)
     intsbl="C 2/m _a3 = I 2/m 1 1"
   case(5)
     intsbl="C 2/m _b2 = A 1 2/m 1"
   case(6)
     intsbl="C 2/m _b3 = I 1 2/m 1"
   case(7)
     intsbl="C 2/m _c1 = A 1 1 2/m"
   case(8)
     intsbl="C 2/m _c2 = B 1 1 2/m = B 2/m"
   case(9)
     intsbl="C 2/m _c3 = I 1 1 2/m"
   case default
     intsbl="intsbl to be determined"
   end select
 case(13)
   brvsb="P"; intsb="2/c"; schsb="C2h^4"; sporder=4
   select case (spgaxor)
   case(1)
     intsbl="P 2/c _b1 = P 1 2/c 1"
   case(2)
     intsbl="P 2/c _a1 = P 2/b 1 1"
   case(3)
     intsbl="P 2/c _a2 = P 2/n 1 1"
   case(4)
     intsbl="P 2/c _a3 = P 2/c 1 1"
   case(5)
     intsbl="P 2/c _b2 = P 1 2/n 1"
   case(6)
     intsbl="P 2/c _b3 = P 1 2/a 1"
   case(7)
     intsbl="P 2/c _c1 = P 1 1 2/a"
   case(8)
     intsbl="P 2/c _c2 = P 1 1 2/n"
   case(9)
     intsbl="P 2/c _c3 = P 1 1 2/b = P 2/b"
   case default
     intsbl="intsbl to be determined"
   end select
 case(14)
   brvsb="P"; intsb="2_1/c"; schsb="C2h^5"; sporder=4
   select case (spgaxor)
   case(1)
     intsbl="P 2_1/c _b1 = P 1 2_1/c 1"
   case(2)
     intsbl="P 2_1/c _a1 = P 2_1/b 1 1"
   case(3)
     intsbl="P 2_1/c _a2 = P 2_1/n 1 1"
   case(4)
     intsbl="P 2_1/c _a3 = P 2_1/c 1 1"
   case(5)
     intsbl="P 2_1/c _b2 = P 1 2_1/n 1"
   case(6)
     intsbl="P 2_1/c _b3 = P 1 2_1/a 1"
   case(7)
     intsbl="P 2_1/c _c1 = P 1 1 2_1/a"
   case(8)
     intsbl="P 2_1/c _c2 = P 1 1 2_1/n"
   case(9)
     intsbl="P 2_1/c _c3 = P 1 1 2_1/b = P 2_1/b"
   case default
     intsbl="intsbl to be determined"
   end select
 case(15)
   brvsb="C"; intsb="2/c"; schsb="C2h^6"; sporder=8
   select case (spgaxor)
   case(1)
     intsbl="C 2/c _b1 = C 1 2/c 1"
   case(2)
     intsbl="C 2/c _a1 = B 2/b 1 1"
   case(3)
     intsbl="C 2/c _a2 = C 2/n 1 1"
   case(4)
     intsbl="C 2/c _a3 = I 2/c 1 1"
   case(5)
     intsbl="C 2/c _b2 = A 1 2/n 1"
   case(6)
     intsbl="C 2/c _b3 = I 1 2/a 1"
   case(7)
     intsbl="C 2/c _c1 = A 1 1 2/a"
   case(8)
     intsbl="C 2/c _c2 = B 1 1 2/n"
   case(9)
     intsbl="C 2/c _c3 = I 1 1 2/b"
   case default
     intsbl="intsbl to be determined"
   end select
 case(16)
   brvsb="P"; intsb="2 2 2"; schsb="D2^1"; sporder=4
 case(17)
   brvsb="P"; intsb="2 2 2_1"; schsb="D2^2"; sporder=4
   select case (spgaxor)
   case(1)
     intsbl="P 2 2 2_1"
   case(2)
     intsbl="P 2_1 2 2"
   case(3)
     intsbl="P 2 2_1 2"
   case default
     intsbl="intsbl to be determined"
   end select
 case(18)
   brvsb="P"; intsb="2_1 2_1 2"; schsb="D2^3"; sporder=4
   select case (spgaxor)
   case(1)
     intsbl="P 2_1 2_1 2"
   case(2)
     intsbl="P 2 2_1 2_1"
   case(3)
     intsbl="P 2_1 2 2_1"
   case default
     intsbl="intsbl to be determined"
   end select
 case(19)
   brvsb="P"; intsb="2_1 2_1 2_1"; schsb="D2^4"; sporder=4
 case(20)
   schsb="D2^5"; sporder=8
   select case (spgaxor)
   case(1)
     brvsb="C"; intsb="2 2 2_1"
   case(2)
     brvsb="A"; intsb="2_1 2 2"
   case(3)
     brvsb="B"; intsb="2 2_1 2"
   case default
     brvsb="X"
     intsbl="intsbl to be determined"
   end select
 case(21)
   schsb="D2^6"; sporder=8
   select case (spgaxor)
   case(1)
     brvsb="C"; intsb="2 2 2"
   case(2)
     brvsb="A"; intsb="2 2 2"
   case(3)
     brvsb="B"; intsb="2 2 2"
   case default
     brvsb="X"
     intsbl="intsbl to be determined"
   end select
 case(22)
   brvsb="F"; intsb="2 2 2"; schsb="D2^7"; sporder=16
 case(23)
   brvsb="I"; intsb="2 2 2"; schsb="D2^8"; sporder=8
 case(24)
   brvsb="I"; intsb="2_1 2_1 2_1"; schsb="D2^9"; sporder=8
 case(25)
   brvsb="P"; schsb="C2v^1"; sporder=4
   select case (spgaxor)
   case(1)
     intsb="m m 2"
   case(2)
     intsb="2 m m"
   case(3)
     intsb="m 2 m"
   case default
     intsb="intsb unknown"
   end select
 case(26)
   brvsb="P"; schsb="C2v^2"; sporder=4
   select case (spgaxor)
   case(1)
     intsb="m c 2_1"
   case(2)
     intsb="2_1 m a"
   case(3)
     intsb="b 2_1 m"
   case(4)
     intsb="m 2_1 b"
   case(5)
     intsb="c m 2_1"
   case(6)
     intsb="2_1 a m"
   case default
     intsb="intsb unknown"
   end select
 case(27)
   brvsb="P"; schsb="C2v^3"; sporder=4
   select case (spgaxor)
   case(1)
     intsb="c c 2"
   case(2)
     intsb="2 a a"
   case(3)
     intsb="b 2 b"
   case default
     intsb="intsb unknown"
   end select
 case(28)
   brvsb="P"; schsb="C2v^4"; sporder=4
   select case (spgaxor)
   case(1)
     intsb="m a 2"
   case(2)
     intsb="2 m b"
   case(3)
     intsb="c 2 m"
   case(4)
     intsb="m 2 a"
   case(5)
     intsb="b m 2"
   case(6)
     intsb="2 c m"
   case default
     intsb="intsb unknown"
   end select
 case(29)
   brvsb="P"; schsb="C2v^5"; sporder=4
   select case (spgaxor)
   case(1)
     intsb="c a 2_1"
   case(2)
     intsb="2_1 a b"
   case(3)
     intsb="c 2_1 b"
   case(4)
     intsb="b 2_1 a"
   case(5)
     intsb="b c 2_1"
   case(6)
     intsb="2_1 c a"
   case default
     intsb="intsb unknown"
   end select
 case(30)
   brvsb="P"; schsb="C2v^6"; sporder=4
   select case (spgaxor)
   case(1)
     intsb="n c 2"
   case(2)
     intsb="2 n a"
   case(3)
     intsb="b 2 n"
   case(4)
     intsb="n 2 b"
   case(5)
     intsb="c n 2"
   case(6)
     intsb="2 a n"
   case default
     intsb="intsb unknown"
   end select
 case(31)
   brvsb="P"; schsb="C2v^7"; sporder=4
   select case (spgaxor)
   case(1)
     intsb="m n 2_1"
   case(2)
     intsb="2_1 m n"
   case(3)
     intsb="n 2_1 m"
   case(4)
     intsb="m 2_1 n"
   case(5)
     intsb="n m 2_1"
   case(6)
     intsb="2_1 n m"
   case default
     intsb="intsb unknown"
   end select
 case(32)
   brvsb="P"; schsb="C2v^8"; sporder=4
   select case (spgaxor)
   case(1)
     intsb="b a 2"
   case(2)
     intsb="2 c b"
   case(3)
     intsb="c 2 a"
   case default
     intsb="intsb unknown"
   end select
 case(33)
   brvsb="P"; schsb="C2v^9"; sporder=4
   select case (spgaxor)
   case(1)
     intsb="n a 2_1"
   case(2)
     intsb="2_1 n b"
   case(3)
     intsb="c 2_1 n"
   case(4)
     intsb="n 2_1 a"
   case(5)
     intsb="b n 2_1"
   case(6)
     intsb="2_1 c n"
   case default
     intsb="intsb unknown"
   end select
 case(34)
   brvsb="P"; schsb="C2v^10"; sporder=4
   select case (spgaxor)
   case(1)
     intsb="n n 2"
   case(2)
     intsb="2 n n"
   case(3)
     intsb="n 2 n"
   case default
     intsb="intsb unknown"
   end select
 case(35)
   schsb="C2v^11"; sporder=8
   select case (spgaxor)
   case(1)
     brvsb="C"; intsb="m m 2"
   case(2)
     brvsb="A"; intsb="2 m m"
   case(3)
     brvsb="B"; intsb="m 2 m"
   case default
     brvsb="X"
     intsb="intsb unknown"
   end select
 case(36)
   schsb="C2v^12"; sporder=8
   select case (spgaxor)
   case(1)
     brvsb="C"; intsb="m c 2_1"
   case(2)
     brvsb="A"; intsb="2_1 m a"
   case(3)
     brvsb="B"; intsb="b 2_1 m"
   case(4)
     brvsb="B"; intsb="m 2_1 b"
   case(5)
     brvsb="C"; intsb="c m 2_1"
   case(6)
     brvsb="A"; intsb="2_1 a m"
   case default
     brvsb="X"
     intsb="intsb unknown"
   end select
 case(37)
   schsb="C2v^13"; sporder=8
   select case (spgaxor)
   case(1)
     brvsb="C"; intsb="c c 2"
   case(2)
     brvsb="A"; intsb="2 a a"
   case(3)
     brvsb="B"; intsb="b 2 b"
   case default
     brvsb="X"
     intsb="intsb unknown"
   end select
 case(38)
   schsb="C2v^14"; sporder=8
   select case (spgaxor)
   case(1)
     brvsb="A"; intsb="m m 2"
   case(2)
     brvsb="B"; intsb="2 m m"
   case(3)
     brvsb="C"; intsb="m 2 m"
   case(4)
     brvsb="A"; intsb="m 2 m"
   case(5)
     brvsb="B"; intsb="m m 2"
   case(6)
     brvsb="C"; intsb="2 m m"
   case default
     brvsb="X"
     intsb="intsb unknown"
   end select
 case(39)
   schsb="C2v^15"; sporder=8
   select case (spgaxor)
   case(1)
     brvsb="A"; intsb="b m 2"
   case(2)
     brvsb="B"; intsb="2 c m"
   case(3)
     brvsb="C"; intsb="m 2 a"
   case(4)
     brvsb="A"; intsb="c 2 m"
   case(5)
     brvsb="B"; intsb="m a 2"
   case(6)
     brvsb="C"; intsb="2 m b"
   case default
     brvsb="X"
     intsb="intsb unknown"
   end select
 case(40)
   schsb="C2v^16"; sporder=8
   select case (spgaxor)
   case(1)
     brvsb="A"; intsb="m a 2"
   case(2)
     brvsb="B"; intsb="2 m b"
   case(3)
     brvsb="C"; intsb="c 2 m"
   case(4)
     brvsb="A"; intsb="m 2 a"
   case(5)
     brvsb="B"; intsb="b m 2"
   case(6)
     brvsb="C"; intsb="2 c m"
   case default
     brvsb="X"
     intsb="intsb unknown"
   end select
 case(41)
   schsb="C2v^17"; sporder=8
   select case (spgaxor)
   case(1)
     brvsb="A"; intsb="b a 2"
   case(2)
     brvsb="B"; intsb="2 c b"
   case(3)
     brvsb="C"; intsb="c 2 a"
   case(4)
     brvsb="A"; intsb="c 2 a"
   case(5)
     brvsb="B"; intsb="b a 2"
   case(6)
     brvsb="C"; intsb="2 c b"
   case default
     brvsb="X"
     intsb="intsb unknown"
   end select
 case(42)
   brvsb="F"; schsb="C2v^18"; sporder=16
   select case (spgaxor)
   case(1)
     intsb="m m 2"
   case(2)
     intsb="2 m m"
   case(3)
     intsb="m 2 m"
   case default
     intsb="intsb unknown"
   end select
 case(43)
   brvsb="F"; schsb="C2v^19"; sporder=16
   select case (spgaxor)
   case(1)
     intsb="d d 2"
   case(2)
     intsb="2 d d"
   case(3)
     intsb="d 2 d"
   case default
     intsb="intsb unknown"
   end select
 case(44)
   brvsb="I"; schsb="C2v^20"; sporder=8
   select case (spgaxor)
   case(1)
     intsb="m m 2"
   case(2)
     intsb="2 m m"
   case(3)
     intsb="m 2 m"
   case default
     intsb="intsb unknown"
   end select
 case(45)
   brvsb="I"; schsb="C2v^21"; sporder=8
   select case (spgaxor)
   case(1)
     intsb="b a 2"
   case(2)
     intsb="2 c b"
   case(3)
     intsb="c 2 a"
   case default
     intsb="intsb unknown"
   end select
 case(46)
   brvsb="I"; schsb="C2v^22"; sporder=8
   select case (spgaxor)
   case(1)
     intsb="m a 2"
   case(2)
     intsb="2 m b"
   case(3)
     intsb="c 2 m"
   case(4)
     intsb="m 2 a"
   case(5)
     intsb="b m 2"
   case(6)
     intsb="2 c m"
   case default
     intsb="intsb unknown"
   end select
 case(47)
   brvsb="P"; intsb="m m m"; schsb="D2h^1"; sporder=8
 case(48)
   brvsb="P"; intsb="n n n"; schsb="D2h^2"; sporder=8
   select case (spgorig)
   case(1)
     intsbl="n n n _1"
   case(2)
     intsbl="n n n _2"
   case default
     intsbl="intsbl to be determined"
   end select
 case(49)
   brvsb="P"; schsb="D2h^3"; sporder=8
   select case (spgaxor)
   case(1)
     intsb="c c m"
   case(2)
     intsb="m a a"
   case(3)
     intsb="b m b"
   case default
     intsb="intsb unknown"
   end select
 case(50)
   brvsb="P"; schsb="D2h^4"; sporder=8
   select case(spgorig)
   case(1)
     select case(spgaxor)
     case(1)
       intsb="b a n"; intsbl="b a n _1"
     case(2)
       intsb="n c b"; intsbl="n c b _1"
     case(3)
       intsb="c n a"; intsbl="c n a _1"
     case default
       intsb="intsb unknown"
       intsbl="intsbl to be determined"
     end select
   case(2)
     select case(spgaxor)
     case(5)
       intsb="b a n"; intsbl="b a n _2"
     case(6)
       intsb="n c b"; intsbl="n c b _2"
     case(4)
       intsb="c n a"; intsbl="c n a _2"
     case default
       intsb="intsb unknown"
       intsbl="intsbl to be determined"
     end select
   case default
     intsb="intsb unknown"
     intsbl="intsbl to be determined"
   end select
 case(51)
   brvsb="P"; schsb="D2h^5"; sporder=8
   select case (spgaxor)
   case(1)
     intsb="m m a"
   case(2)
     intsb="b m m"
   case(3)
     intsb="m c m"
   case(4)
     intsb="m a m"
   case(5)
     intsb="m m b"
   case(6)
     intsb="c m m"
   case default
     intsb="intsb unknown"
   end select
 case(52)
   brvsb="P"; schsb="D2h^6"; sporder=8
   select case (spgaxor)
   case(1)
     intsb="n n a"
   case(2)
     intsb="b n n"
   case(3)
     intsb="n c n"
   case(4)
     intsb="n a n"
   case(5)
     intsb="n n b"
   case(6)
     intsb="c n n"
   case default
     intsb="intsb unknown"
   end select
 case(53)
   brvsb="P"
   schsb="D2h^7"
   sporder=8
   select case (spgaxor)
   case(1)
     intsb="m n a"
   case(2)
     intsb="b m n"
   case(3)
     intsb="n c m"
   case(4)
     intsb="m a n"
   case(5)
     intsb="n m b"
   case(6)
     intsb="c n m"
   case default
     intsb="intsb unknown"
   end select
 case(54)
   brvsb="P"; schsb="D2h^8"; sporder=8
   select case (spgaxor)
   case(1)
     intsb="c c a"
   case(2)
     intsb="b a a"
   case(3)
     intsb="b c b"
   case(4)
     intsb="b a b"
   case(5)
     intsb="c c b"
   case(6)
     intsb="c a a"
   case default
     intsb="intsb unknown"
   end select
 case(55)
   brvsb="P"; schsb="D2h^9"; sporder=8
   select case (spgaxor)
   case(1)
     intsb="b a m"
   case(2)
     intsb="m c b"
   case(3)
     intsb="c m a"
   case default
     intsb="intsb unknown"
   end select
 case(56)
   brvsb="P"; schsb="D2h^10"; sporder=8
   select case (spgaxor)
   case(1)
     intsb="c c n"
   case(2)
     intsb="n a a"
   case(3)
     intsb="b n b"
   case default
     intsb="intsb unknown"
   end select
 case(57)
   brvsb="P"; schsb="D2h^11"; sporder=8
   select case (spgaxor)
   case(1)
     intsb="b c m"
   case(2)
     intsb="m c a"
   case(3)
     intsb="b m a"
   case(4)
     intsb="c m b"
   case(5)
     intsb="c a m"
   case(6)
     intsb="m a b"
   case default
     intsb="intsb unknown"
   end select
 case(58)
   brvsb="P"; schsb="D2h^12"; sporder=8
   select case (spgaxor)
   case(1)
     intsb="n n m"
   case(2)
     intsb="m n n"
   case(3)
     intsb="n m n"
   case default
     intsb="intsb unknown"
   end select
 case(59)
   brvsb="P"; schsb="D2h^13"; sporder=8
   if (spgorig==1) then
     select case (spgaxor)
     case(1)
       intsb="m m n"; intsbl="m m n _1"
     case(2)
       intsb="m m n"; intsbl="n m m _1"
     case(3)
       intsb="m m n"; intsbl="m n m _1"
     case default
       intsb="intsb unknown"
       intsbl="intsbl to be determined"
     end select
   else if(spgorig==2) then
     select case (spgaxor)
     case(5)
       intsb="m m n"; intsbl="m m n _2"
     case(6)
       intsb="m m n"; intsbl="n m m _2"
     case(4)
       intsb="m m n"; intsbl="m n m _2"
     case default
       intsb="intsb unknown"
       intsbl="intsbl to be determined"
     end select
   else
     intsb="intsb unknown"
     intsbl="intsbl to be determined"
   end if
 case(60)
   brvsb="P"; schsb="D2h^14"; sporder=8
   select case (spgaxor)
   case(1)
     intsb="b c n"
   case(2)
     intsb="n c a"
   case(3)
     intsb="b n a"
   case(4)
     intsb="c n b"
   case(5)
     intsb="c a n"
   case(6)
     intsb="n a b"
   case default
     intsb="intsb unknown"
   end select
 case(61)
   brvsb="P"; schsb="D2h^15"; sporder=8
   if (spgaxor==1)then
     intsb="b c a"
   else if (spgaxor==2)then
     intsb="c a b"
   else
     intsb="intsb unknown"
   end if
 case(62)
   brvsb="P"; schsb="D2h^16"; sporder=8
   select case (spgaxor)
   case(1)
     intsb="n m a"
   case(2)
     intsb="b n m"
   case(3)
     intsb="m c n"
   case(4)
     intsb="n a m"
   case(5)
     intsb="m n b"
   case(6)
     intsb="c m n"
   case default
     intsb="intsb unknown"
   end select
 case(63)
   schsb="D2h^17"; sporder=16
   select case (spgaxor)
   case(1)
     brvsb="C"; intsb="m c m"
   case(2)
     brvsb="A"; intsb="m m a"
   case(3)
     brvsb="B"; intsb="b m m"
   case(4)
     brvsb="B"; intsb="m m b"
   case(5)
     brvsb="C"; intsb="c m m"
   case(6)
     brvsb="A"; intsb="m a m"
   case default
     brvsb="X"
     intsbl="intsbl unknown"
   end select
 case(64)
   schsb="D2h^18"; sporder=16
   select case (spgaxor)
   case(1)
     brvsb="C"; intsb="m c a"
   case(2)
     brvsb="A"; intsb="b m a"
   case(3)
     brvsb="B"; intsb="b c m"
   case(4)
     brvsb="B"; intsb="m a b"
   case(5)
     brvsb="C"; intsb="c m b"
   case(6)
     brvsb="A"; intsb="c a m"
   case default
     brvsb="X"
     intsb="intsb unknown"
   end select
 case(65)
   schsb="D2h^19"; sporder=16
   select case (spgaxor)
   case(1)
     brvsb="C"; intsb="m m m"
   case(2)
     brvsb="A"; intsb="m m m"
   case(3)
     brvsb="B"; intsb="m m m"
   case default
     brvsb="X"
     intsb="intsb unknown"
   end select
 case(66)
   schsb="D2h^20"; sporder=16
   select case (spgaxor)
   case(1)
     brvsb="C"; intsb="c c m"
   case(2)
     brvsb="A"; intsb="m a a"
   case(3)
     brvsb="B"; intsb="b m b"
   case default
     brvsb="X"
     intsb="intsb unknown"
   end select
 case(67)
   schsb="D2h^21"; sporder=16
   select case (spgaxor)
   case(1)
     brvsb="C"; intsb="m m a"
   case(2)
     brvsb="A"; intsb="b m m"
   case(3)
     brvsb="B"; intsb="m c m"
   case(4)
     brvsb="B"; intsb="m a m"
   case(5)
     brvsb="C"; intsb="m m b"
   case(6)
     brvsb="A"; intsb="c m m"
   case default
     brvsb="X"
     intsb="intsb unknown"
   end select
 case(68)
   schsb="D2h^22"; sporder=16
   if (spgorig==1) then
     select case (spgaxor)
     case(1)
       brvsb="C"; intsb="c c a"; intsbl="c c a _1"
     case(2)
       brvsb="A"; intsb="b a a"; intsbl="b a a _1"
     case(3)
       brvsb="B"; intsb="b c b"; intsbl="b c b _1"
     case(4)
       brvsb="B"; intsb="b a b"; intsbl="b a b _1"
     case(5)
       brvsb="C"; intsb="c c b"; intsbl="c c b _1"
     case(6)
       brvsb="A"; intsb="c a a"; intsbl="c a a _1"
     case default
       brvsb="X"
       intsb="intsb unknown"
       intsbl="intsbl to be determined"
     end select
   else if(spgorig==2)then
     select case (spgaxor)
     case(1)
       brvsb="C"; intsb="c c a"; intsbl="c c a _2"
     case(2)
       brvsb="A"; intsb="b a a"; intsbl="b a a _2"
     case(3)
       brvsb="B"; intsb="b c b"; intsbl="b c b _2"
     case(4)
       brvsb="B"; intsb="b a b"; intsbl="b a b _2"
     case(5)
       brvsb="C"; intsb="c c b"; intsbl="c c b _2"
     case(6)
       brvsb="A"; intsb="c a a"; intsbl="c a a _2"
     case default
       brvsb="X"
       intsb="intsb unknown"
       intsbl="intsbl to be determined"
     end select
   else
     brvsb="X"
     intsb="intsb unknown"
     intsbl="intsbl to be determined"
   end if
 case(69)
   brvsb="F"; intsb="m m m"; schsb="D2h^23"; sporder=32
 case(70)
   brvsb="F"; intsb="d d d"; schsb="D2h^24"; sporder=32
   if (spgorig==1)then
     intsbl="d d d _1"
   else if (spgorig==2)then
     intsbl="d d d _2"
   else
     intsbl="intsbl to be determined"
   end if
 case(71)
   brvsb="I"; intsb="m m m"; schsb="D2h^25"; sporder=16
 case(72)
   brvsb="I"; schsb="D2h^26"; sporder=16
   select case (spgaxor)
   case(1)
     intsb="b a m"
   case(2)
     intsb="m c b"
   case(3)
     intsb="c m a"
   case default
     intsb="intsb unknown"
   end select
 case(73)
   brvsb="I"; schsb="D2h^27"; sporder=16
   if (spgorig==1)then
     intsb="b c a"
   else if (spgorig==2)then
     intsb="c a b"
   else
     intsb="intsb unknown"
   end if
 case(74)
   brvsb="I"; schsb="D2h^28"; sporder=16
   select case (spgaxor)
   case(1)
     intsb="m m a"
   case(2)
     intsb="b m m"
   case(3)
     intsb="m c m"
   case(4)
     intsb="m a m"
   case(5)
     intsb="m m b"
   case(6)
     intsb="c m m"
   case default
     intsb="intsb unknown"
   end select
 case(75)
   brvsb="P"; intsb="4"; schsb="C4^1"; sporder=4
 case(76)
   brvsb="P"; intsb="4_1"; schsb="C4^2"; sporder=4
 case(77)
   brvsb="P"; intsb="4_2"; schsb="C4^3"; sporder=4
 case(78)
   brvsb="P"; intsb="4_3"; schsb="C4^4"; sporder=4
 case(79)
   brvsb="I"; intsb="4"; schsb="C4^5"; sporder=8
 case(80)
   brvsb="I"; intsb="4_1"; schsb="C4^6"; sporder=8
 case(81)
   brvsb="P"; intsb="-4"; schsb="S4^1"; sporder=4
 case(82)
   brvsb="I"; intsb="-4"; schsb="S4^2"; sporder=8
 case(83)
   brvsb="P"; intsb="4/m"; schsb="C4h^1"; sporder=8
 case(84)
   brvsb="P"; intsb="4_2/m"; schsb="C4h^2"; sporder=8
 case(85)
   brvsb="P"; intsb="4/n"; schsb="C4h^3"; sporder=8
   if (spgorig==1)then
     intsbl="4/n _1"
   else if (spgorig==2)then
     intsbl="4/n _2"
   else
     intsbl="intsbl to be determined"
   end if
 case(86)
   brvsb="P"; intsb="4_2/n"; schsb="C4h^4"; sporder=8
   if (spgorig==1)then
     intsbl="4_2/n _1"
   else if (spgorig==2)then
     intsbl="4_2/n _2"
   else
     intsbl="intsbl to be determined"
   end if
 case(87)
   brvsb="I"; intsb="4/m"; schsb="C4h^5"; sporder=16
 case(88)
   brvsb="I"; intsb="4_1/a"; schsb="C4h^6"; sporder=16
   if (spgorig==1)then
     intsbl="4_1/a _1"
   else if (spgorig==2)then
     intsbl="4_1/a _2"
   else
     intsbl="intsbl to be determined"
   end if
 case(89)
   brvsb="P"; intsb="4 2 2"; schsb="D4^1"; sporder=8
 case(90)
   brvsb="P"; intsb="4 2_1 2"; schsb="D4^2"; sporder=8
 case(91)
   brvsb="P"; intsb="4_1 2 2"; schsb="D4^3"; sporder=8
 case(92)
   brvsb="P"; intsb="4_1 2_1 2"; schsb="D4^4"; sporder=8
 case(93)
   brvsb="P"; intsb="4_2 2 2"; schsb="D4^5"; sporder=8
 case(94)
   brvsb="P"; intsb="4_2 2_1 2"; schsb="D4^6"; sporder=8
 case(95)
   brvsb="P"; intsb="4_3 2 2"; schsb="D4^7"; sporder=8
 case(96)
   brvsb="P"; intsb="4_3 2_1 2"; schsb="D4^8"; sporder=8
 case(97)
   brvsb="I"; intsb="4 2 2"; schsb="D4^9"; sporder=16
 case(98)
   brvsb="I"; intsb="4_1 2 2"; schsb="D4^10"; sporder=16
 case(99)
   brvsb="P"; intsb="4 m m"; schsb="C4v^1"; sporder=8
 case(100)
   brvsb="P"; intsb="4 b m"; schsb="C4v^2"; sporder=8
 case(101)
   brvsb="P"; intsb="4_2 c m"; schsb="C4v^3"; sporder=8
 case(102)
   brvsb="P"; intsb="4_2 n m"; schsb="C4v^4"; sporder=8
 case(103)
   brvsb="P"; intsb="4 c c"; schsb="C4v^5"; sporder=8
 case(104)
   brvsb="P"; intsb="4 n c"; schsb="C4v^6"; sporder=8
 case(105)
   brvsb="P"; intsb="4_2 m c"; schsb="C4v^7"; sporder=8
 case(106)
   brvsb="P"; intsb="4_2 b c"; schsb="C4v^8"; sporder=8
 case(107)
   brvsb="I"; intsb="4 m m"; schsb="C4v^9"; sporder=16
 case(108)
   brvsb="I"; intsb="4 c m"; schsb="C4v^10"; sporder=16
 case(109)
   brvsb="I"; intsb="4_1 m d"; schsb="C4v^11"; sporder=16
 case(110)
   brvsb="I"; intsb="4_1 c d"; schsb="C4v^12"; sporder=16
 case(111)
   brvsb="P"; intsb="-4 2 m"; schsb="D2d^1"; sporder=8
 case(112)
   brvsb="P"; intsb="-4 2 c"; schsb="D2d^2"; sporder=8
 case(113)
   brvsb="P"; intsb="-4 2_1 m"; schsb="D2d^3"; sporder=8
 case(114)
   brvsb="P"; intsb="-4 2_1 c"; schsb="D2d^4"; sporder=8
 case(115)
   brvsb="P"; intsb="-4 m 2"; schsb="D2d^5"; sporder=8
 case(116)
   brvsb="P"; intsb="-4 c 2"; schsb="D2d^6"; sporder=8
 case(117)
   brvsb="P"; intsb="-4 b 2"; schsb="D2d^7"; sporder=8
 case(118)
   brvsb="P"; intsb="-4 n 2"; schsb="D2d^8"; sporder=8
 case(119)
   brvsb="I"; intsb="-4 m 2"; schsb="D2d^9"; sporder=16
 case(120)
   brvsb="I"; intsb="-4 c 2"; schsb="D2d^10"; sporder=16
 case(121)
   brvsb="I"; intsb="-4 2 m"; schsb="D2d^11"; sporder=16
 case(122)
   brvsb="I"; intsb="-4 2 d"; schsb="D2d^12"; sporder=16
 case(123)
   brvsb="P"; intsb="4/m m m"; schsb="D4h^1"; sporder=16
 case(124)
   brvsb="P"; intsb="4/m c c"; schsb="D4h^2"; sporder=16
 case(125)
   brvsb="P"; intsb="4/n b m"; schsb="D4h^3"; sporder=16
   if (spgorig==1)then
     intsbl="4/n b m _1"
   else if (spgorig==2)then
     intsbl="4/n b m _2"
   else
     intsbl="intsbl to be determined"
   end if
 case(126)
   brvsb="P"; intsb="4/n n c"; schsb="D4h^4"; sporder=16
   if (spgorig==1)then 
     intsbl="4/n n c _1"
   else if (spgorig==2)then
     intsbl="4/n n c _2"
   else
     intsbl="intsbl to be determined"
   end if
 case(127)
   brvsb="P"; intsb="4/m b m"; schsb="D4h^5"; sporder=16
 case(128)
   brvsb="P"; intsb="4/m n c"; schsb="D4h^6"; sporder=16
 case(129)
   brvsb="P"; intsb="4/n m m"; schsb="D4h^7"; sporder=16
   if (spgorig==1)then
     intsbl="4/n m m _1"
   else if (spgorig==2)then
     intsbl="4/n m m _2"
   else
     intsbl="intsbl to be determined"
   end if
 case(130)
   brvsb="P"; intsb="4/n c c"; schsb="D4h^8"; sporder=16
   if (spgorig==1)then
     intsbl="4/n c c _1"
   else if (spgorig==2) then
     intsbl="4/n c c _2"
   else   
     intsbl="intsbl to be determined"
   end if
 case(131)
   brvsb="P"; intsb="4_2/m m c"; schsb="D4h^9"; sporder=16
 case(132)
   brvsb="P"; intsb="4_2/m c m"; schsb="D4h^10"; sporder=16
 case(133)
   brvsb="P"; intsb="4_2/n b c"; schsb="D4h^11"; sporder=16
   if (spgorig==1)then
     intsbl="4_2/n b c _1"
   else if (spgorig==2)then 
     intsbl="4_2/n b c _2"
   else
     intsbl="intsbl to be determined"
   end if
 case(134)
   brvsb="P"; intsb="4_2/n n m"; schsb="D4h^12"; sporder=16
   if (spgorig==1)then 
     intsbl="4_2/n n m _1"
   else if (spgorig==2)then 
     intsbl="4_2/n n m _2"
   else
     intsbl="intsbl to be determined"
   end if
 case(135)
   brvsb="P"; intsb="4_2/m b c"; schsb="D4h^13"; sporder=16
 case(136)
   brvsb="P"; intsb="4_2/m n m"; schsb="D4h^14"; sporder=16
 case(137)
   brvsb="P"; intsb="4_2/n m c"; schsb="D4h^15"; sporder=16
   if (spgorig==1)then  
     intsbl="4_2/n m c _1"
   else if (spgorig==2)then 
     intsbl="4_2/n m c _2"
   else
     intsbl="intsbl to be determined" 
   end if
 case(138)
   brvsb="P"; intsb="4_2/n c m"; schsb="D4h^16"; sporder=16
   if (spgorig==1)then 
     intsbl="4_2/n c m _1"
   else if (spgorig==2)then 
     intsbl="4_2/n c m _2"
   else
     intsbl="intsbl to be determined"
   end if
 case(139)
   brvsb="I"; intsb="4/m m m"; schsb="D4h^17"; sporder=32
 case(140)
   brvsb="I"; intsb="4/m c m"; schsb="D4h^18"; sporder=32
 case(141)
   brvsb="I"; intsb="4_1/a m d"; schsb="D4h^19"; sporder=32
   if (spgorig==1)then
     intsbl="4_1/a m d _1"
   else if (spgorig==2)then 
     intsbl="4_1/a m d _2"
   else
     intsbl="intsbl to be determined"
   end if
 case(142)
   brvsb="I"; intsb="4_1/a c d"; schsb="D4h^20"; sporder=32
   if (spgorig==1)then
     intsbl="4_1/a c d _1"
   else if (spgorig==2)then
     intsbl="4_1/a c d _2"
   else
     intsbl="intsbl to be determined"
   end if
 case(143)
   brvsb="P"; intsb="3"; schsb="C3^1"; sporder=3
 case(144)
   brvsb="P"; intsb="3_1"; schsb="C3^2"; sporder=3
 case(145)
   brvsb="P"; intsb="3_2"; schsb="C3^3"; sporder=3
 case(146)
   brvsb="R"; intsb="3"; schsb="C3^4"
   if (spgorig==1)then
     intsbl="3 _H" ; sporder=9
   else if (spgorig==2)then
     intsbl="3 _R" ; sporder=3
   else
     intsbl="intsbl to be determined"
   end if
 case(147)
   brvsb="P"; intsb="-3"; schsb="C3i^1"; sporder=6
 case(148)
   brvsb="R"; intsb="-3"; schsb="C3i^2"
   if (spgorig==1) then
     intsbl="-3 _H" ; sporder=9
   else if (spgorig==2) then 
     intsbl="-3 _R" ; sporder=3
   else
     intsbl="intsbl to be determined"
   end if
 case(149)
   brvsb="P"; intsb="3 1 2"; schsb="D3^1"; sporder=6
 case(150)
   brvsb="P"; intsb="3 2 1"; schsb="D3^2"; sporder=6
 case(151)
   brvsb="P"; intsb="3_1 1 2"; schsb="D3^3"; sporder=6
 case(152)
   brvsb="P"; intsb="3_1 2 1"; schsb="D3^4"; sporder=6
 case(153)
   brvsb="P"; intsb="3_2 1 2"; schsb="D3^5"; sporder=6
 case(154)
   brvsb="P"; intsb="3_2 2 1"; schsb="D3^6"; sporder=6
 case(155)
   brvsb="R"; intsb="3 2"; schsb="D3^7"
   if (spgorig==1) then
     intsbl="3 2 _H" ; sporder=18
   else if (spgorig==2) then
     intsbl="3 2 _R" ; sporder=6
   else
     intsbl="intsbl to be determined"
   end if
 case(156)
   brvsb="P"; intsb="3 m 1"; schsb="C3v^1"; sporder=6
 case(157)
   brvsb="P"; intsb="3 1 m"; schsb="C3v^2"; sporder=6
 case(158)
   brvsb="P"; intsb="3 c 1"; schsb="C3v^3"; sporder=6
 case(159)
   brvsb="P"; intsb="3 1 c"; schsb="C3v^4"; sporder=6
 case(160)
   brvsb="R"; intsb="3 m"; schsb="C3v^5"
   if (spgorig==1) then
     intsbl="3 m _H" ; sporder=18
   else if (spgorig==2) then
     intsbl="3 m _R" ; sporder=6
   else
     intsbl="intsbl to be determined"
   end if
 case(161)
   brvsb="R"; intsb="3 c"; schsb="C3v^6"
   if (spgorig==1) then
     intsbl="3 m _H" ; sporder=18
   else if (spgorig==2)then
     intsbl="3 m _R" ; sporder=6
   else
     intsbl="intsbl to be determined"
   end if
 case(162)
   brvsb="P"; intsb="-3 1 m"; schsb="D3d^1"; sporder=12
 case(163)
   brvsb="P"; intsb="-3 1 c"; schsb="D3d^2"; sporder=12
 case(164)
   brvsb="P"; intsb="-3 m 1"; schsb="D3d^3"; sporder=12
 case(165)
   brvsb="P"; intsb="-3 c 1"; schsb="D3d^4"; sporder=12
 case(166)
   brvsb="R"; intsb="-3 m"; schsb="D3d^5"
   if (spgorig==1) then
     intsbl="3 m _H"; sporder=18
   else if (spgorig==2) then
     intsbl="3 m _R"; sporder=6 
   else
     intsbl="intsbl to be determined"
   end if
 case(167)
   brvsb="R"; intsb="-3 c"; schsb="D3d^6"
   if (spgorig==1) then
     intsbl="-3 c _H"; sporder=36
   else if (spgorig==2) then 
     intsbl="-3 c _R"; sporder=12
   else
     intsbl="intsbl to be determined"
     sporder=-1
   end if
 case(168)
   brvsb="P"; intsb="6"; schsb="C6^1"; sporder=6
 case(169)
   brvsb="P"; intsb="6_1"; schsb="C6^2"; sporder=6
 case(170)
   brvsb="P"; intsb="6_5"; schsb="C6^3"; sporder=6
 case(171)
   brvsb="P"; intsb="6_2"; schsb="C6^4"; sporder=6
 case(172)
   brvsb="P"; intsb="6_4"; schsb="C6^5"; sporder=6
 case(173)
   brvsb="P"; intsb="6_3"; schsb="C6^6"; sporder=6
 case(174)
   brvsb="P"; intsb="-6"; schsb="C3h^1"; sporder=6
 case(175)
   brvsb="P"; intsb="6/m"; schsb="C6h^1"; sporder=12
 case(176)
   brvsb="P"; intsb="6_3/m"; schsb="C6h^2"; sporder=12
 case(177)
   brvsb="P"; intsb="6 2 2"; schsb="D6^1"; sporder=12
 case(178)
   brvsb="P"; intsb="6_1 2 2"; schsb="D6^2"; sporder=12
 case(179)
   brvsb="P"; intsb="6_5 2 2"; schsb="D6^3"; sporder=12
 case(180)
   brvsb="P"; intsb="6_2 2 2"; schsb="D6^4"; sporder=12
 case(181)
   brvsb="P"; intsb="6_4 2 2"; schsb="D6^5"; sporder=12
 case(182)
   brvsb="P"; intsb="6_3 2 2"; schsb="D6^6"; sporder=12
 case(183)
   brvsb="P"; intsb="6 m m"; schsb="C6v^1"; sporder=12
 case(184)
   brvsb="P"; intsb="6 c c"; schsb="C6v^2"; sporder=12
 case(185)
   brvsb="P"; intsb="6_3 c m"; schsb="C6v^3"; sporder=12
 case(186)
   brvsb="P"; intsb="6_3 m c"; schsb="C6v^4"; sporder=12
 case(187)
   brvsb="P"; intsb="-6 m 2"; schsb="D3h^1"; sporder=12
 case(188)
   brvsb="P"; intsb="-6 c 2"; schsb="D3h^2"; sporder=12
 case(189)
   brvsb="P"; intsb="-6 2 m"; schsb="D3h^3"; sporder=12
 case(190)
   brvsb="P"; intsb="-6 2 c"; schsb="D3h^4"; sporder=12
 case(191)
   brvsb="P"; intsb="6/m m m"; schsb="D6h^1"; sporder=24
 case(192)
   brvsb="P"; intsb="6/m c c"; schsb="D6h^2"; sporder=24
 case(193)
   brvsb="P"; intsb="6_3/m c m"; schsb="D6h^3"; sporder=24
 case(194)
   brvsb="P"; intsb="6_3/m m c"; schsb="D6h^4"; sporder=24
 case(195)
   brvsb="P"; intsb="2 3"; schsb="T^1"; sporder=12
 case(196)
   brvsb="F"; intsb="2 3"; schsb="T^2"; sporder=48
 case(197)
   brvsb="I"; intsb="2 3"; schsb="T^3"; sporder=24
 case(198)
   brvsb="P"; intsb="2_1 3"; schsb="T^4"; sporder=12
 case(199)
   brvsb="I"; intsb="2_1 3"; schsb="T^5"; sporder=24
 case(200)
   brvsb="P"; intsb="m -3"; schsb="Th^1"; sporder=24
 case(201)
   brvsb="P"; intsb="n -3"; schsb="Th^2"; sporder=24
   if (spgorig==1) then
     intsbl="n -3 _1"
   else if (spgorig==2)then 
     intsbl="n -3 _2"
   else
     intsbl="intsbl to be determined"
   end if
 case(202)
   brvsb="F"; intsb="m -3"; schsb="Th^3"; sporder=96
 case(203)
   brvsb="F"; intsb="d -3"; schsb="Th^4"; sporder=96
   if (spgorig==1) then
     intsbl="d -3 _1"
   else if (spgorig==2) then
     intsbl="d -3 _2"
   else
     intsbl="intsbl to be determined"
   end if
 case(204)
   brvsb="I"; intsb="m -3"; schsb="Th^5"; sporder=48
 case(205)
   brvsb="P"; intsb="a -3"; schsb="Th^6"; sporder=24
 case(206)
   brvsb="I"; intsb="a -3"; schsb="Th^7"; sporder=48
 case(207)
   brvsb="P"; intsb="4 3 2"; schsb="O^1"; sporder=24
 case(208)
   brvsb="P"; intsb="4_2 3 2"; schsb="O^2"; sporder=24
 case(209)
   brvsb="F"; intsb="4 3 2"; schsb="O^3"; sporder=96
 case(210)
   brvsb="F"; intsb="4_1 3 2"; schsb="O^4"; sporder=96
 case(211)
   brvsb="I"; intsb="4 3 2"; schsb="O^5"; sporder=48
 case(212)
   brvsb="P"; intsb="4_3 3 2"; schsb="O^6"; sporder=24
 case(213)
   brvsb="P"; intsb="4_1 3 2"; schsb="O^7"; sporder=24
 case(214)
   brvsb="I"; intsb="4_1 3 2"; schsb="O^8"; sporder=48
 case(215)
   brvsb="P"; intsb="-4 3 m"; schsb="Td^1"; sporder=24
 case(216)
   brvsb="F"; intsb="-4 3 m"; schsb="Td^2"; sporder=96
 case(217)
   brvsb="I"; intsb="-4 3 m"; schsb="Td^3"; sporder=48
 case(218)
   brvsb="P"; intsb="-4 3 n"; schsb="Td^4"; sporder=24
 case(219)
   brvsb="F"; intsb="-4 3 c"; schsb="Td^5"; sporder=96
 case(220)
   brvsb="I"; intsb="-4 3 d"; schsb="Td^6"; sporder=48
 case(221)
   brvsb="P"; intsb="m -3 m"; schsb="Oh^1"; sporder=48
 case(222)
   brvsb="P"; intsb="n -3 n"; schsb="Oh^2"; sporder=48
   if (spgorig==1) then
     intsbl="n -3 n _1"
   else if (spgorig==2) then 
     intsbl="n -3 n _2"
   else
     intsbl="intsbl to be determined"
   end if
 case(223)
   brvsb="P"; intsb="m -3 n"; schsb="Oh^3"; sporder=48
 case(224)
   brvsb="P"; intsb="n -3 m"; schsb="Oh^4"; sporder=48
   if (spgorig==1) then
     intsbl="n -3 m _1"
   else if (spgorig==2)then
     intsbl="n -3 m _2"
   else
     intsbl="intsbl to be determined"
   end if
 case(225)
   brvsb="F"; intsb="m -3 m"; schsb="Oh^5"; sporder=192
 case(226)
   brvsb="F"; intsb="m -3 c"; schsb="Oh^6"; sporder=192
 case(227)
   brvsb="F"; intsb="d -3 m"; schsb="Oh^7"; sporder=192
   if (spgorig==1) then
     intsbl="d -3 m _1"
   else if (spgorig==2) then
     intsbl="d -3 m _2"
   else
     intsbl="intsbl to be determined"
   end if
 case(228)
   brvsb="F"; intsb="d -3 c"; schsb="Oh^8"; sporder=192
   if (spgorig==1) then
     intsbl="d -3 c _1"
   else if (spgorig==2) then 
     intsbl="d -3 c _2"
   else
     intsbl="intsbl to be determined"
   end if
 case(229)
   brvsb="I"; intsb="m -3 m"; schsb="Oh^9"; sporder=96
 case(230)
   brvsb="I"; intsb="a -3 d"; schsb="Oh^10"; sporder=96
 end select

 if(trim(intsbl)=="same")intsbl=intsb

!Assignment of the point group number
 if(spgroup<=2)then  ! Triclinic system
   select case(spgroup)
   case (1)
     ptintsb="1"; ptschsb="C1"
   case (2)
     ptintsb="-1"; ptschsb="Ci"
   end select
 else if(spgroup<=15)then  ! Monoclinic system
   select case(spgroup)
   case (3:5)
     ptintsb="2"; ptschsb="C2"
   case (6:9)
     ptintsb="m"; ptschsb="Cs = C1h "
   case (10:15)
     ptintsb="2/m"; ptschsb="C2h"
   end select
 else if(spgroup<=74)then  ! Orthorhombic system
   select case(spgroup)
   case (16:24)
     ptintsb="2 2 2"; ptschsb="D2"
   case (25:46)
     ptintsb="m m 2"; ptschsb="C2v"
   case (47:74)
     ptintsb="m m m"; ptschsb="D2h"
   end select
 else if(spgroup<=142)then  ! Tetragonal system
   select case(spgroup)
   case (75:80)
     ptintsb="4"; ptschsb="C4"
   case (81,82)
     ptintsb="-4"; ptschsb="S4"
   case (83:88)
     ptintsb="4/m"; ptschsb="C4h"
   case (89:98)
     ptintsb="4 2 2"; ptschsb="D4"
   case (99:110)
     ptintsb="4 m m"; ptschsb="C4v"
   case (111:114,121,122)
     ptintsb="-4 2 m"; ptschsb="D2d^1"
   case (115:120)
     ptintsb="-4 m 2"; ptschsb="D2h^2"
   case (123:142)
     ptintsb="4/m m m"
     ptschsb="D4h"
   end select
 else if(spgroup<=167)then  ! Trigonal system
   select case(spgroup)
   case (143:146)
     ptintsb="3"; ptschsb="C3"
   case (147,148)
     ptintsb="-3"; ptschsb="C3i"
   case (149,151,153)
     ptintsb="3 1 2"; ptschsb="D3^1"
   case (150,152,154,155)
     ptintsb="3 2 1"; ptschsb="D3^2"
   case (156,158,160,161)
     ptintsb="3 m 1"; ptschsb="C3v^1"
   case (157,159)
     ptintsb="3 1 m"; ptschsb="C3v^2"
   case (162,163)
     ptintsb="-3 1 m"; ptschsb="D3d^1"
   case (164:167)
     ptintsb="-3 m 1"; ptschsb="D3d^2"
   end select
 else if(spgroup<=194)then  ! Hexagonal system
   select case(spgroup)
   case (168:173)
     ptintsb="6"; ptschsb="C6"
   case (174)
     ptintsb="-6"; ptschsb="C3h"
   case (175,176)
     ptintsb="6/m"; ptschsb="C6h"
   case (177:182)
     ptintsb="6 2 2"; ptschsb="D6"
   case (183:186)
     ptintsb="6 m m"; ptschsb="C6v"
   case (187,188)
     ptintsb="-6 m 2"; ptschsb="D3h^1"
   case (189,190)
     ptintsb="-6 2 m"; ptschsb="D3h^2"
   case (191:194)
     ptintsb="6/m m m"; ptschsb="D6h"
   end select
 else                        ! Cubic system
   select case(spgroup)
   case (195:199)
     ptintsb="2 3"; ptschsb="T"
   case (200:206)
     ptintsb="m 3"; ptschsb="Th"
   case (207:214)
     ptintsb="4 3 2"; ptschsb="O"
   case (215:220)
     ptintsb="4 3 m"; ptschsb="Td"
   case (221:230)
     ptintsb="m -3 m"; ptschsb="Oh"
   end select
 end if

end subroutine spgdata
!!***
