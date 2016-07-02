!{\src2tex{textfont=tt}}
!!****f* ABINIT/gensymshub
!! NAME
!! gensymshub
!!
!! FUNCTION
!! Analyse the Shubnikov space group, from the input of the
!! the Fedorov space group number, and the Shubnikov
!! space group number :
!! 1) determine the type (III or IV)
!! 2) for type (IV), determine the translation generating
!!   the anti-ferromagnetic operations
!! At present, follow strictly the specifications of Bradley
!! and Cracknell. However, one should take into account the
!! orientation of the symmetry group (spgaxor).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (RC,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! spgroup = number of space group
!! spgroupma = number of magnetic space group
!!
!! OUTPUT
!! genafm(3) = in case of shubnikov type IV, translation, generator of the
!!  anti-ferromagnetic symmetry operations
!! shubnikov = type of the shubnikov group
!!
!! PARENTS
!!      ingeo
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine gensymshub(genafm,spgroup,spgroupma,shubnikov)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gensymshub'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spgroup,spgroupma
 integer,intent(out) :: shubnikov
!arrays
 real(dp),intent(out) :: genafm(3)

!Local variables ------------------------------
!scalars
 integer :: brvlttbw=0,spgrmatch=1
 character(len=500) :: message

! *************************************************************************

!List of the input parameters
!DEBUG
!write(std_out,*)
!write(std_out,*)' gensymshub : enter with:'
!write(std_out,*)' brvlttbw  = ',brvlttbw
!write(std_out,*)' spgroup = ',spgroup
!write(std_out,*)' spgroupma = ',spgroupma
!ENDDEBUG

!Test for consistency the magnetic and non-magnetic space group
 select case (spgroup)
 case(1)
   if (.not.(spgroupma>=3  .and.  3 >=spgroupma) ) spgrmatch=0
 case(2)
   if (.not.(spgroupma>=6  .and.  7 >=spgroupma) ) spgrmatch=0
 case(3)
   if (.not.(spgroupma>=3   .and.  6 >=spgroupma) ) spgrmatch=0
 case(4)
   if (.not.(spgroupma>=9   .and.  12 >=spgroupma) ) spgrmatch=0
 case(5)
   if (.not.(spgroupma>=15   .and.  17 >=spgroupma) ) spgrmatch=0
 case(6)
   if (.not.(spgroupma>=20   .and.  23 >=spgroupma) ) spgrmatch=0
 case(7)
   if (.not.(spgroupma>=26   .and.  31 >=spgroupma) ) spgrmatch=0
 case(8)
   if (.not.(spgroupma>=34   .and.  36 >=spgroupma) ) spgrmatch=0
 case(9)
   if (.not.(spgroupma>=39   .and.  41 >=spgroupma) ) spgrmatch=0
 case(10)
   if (.not.(spgroupma>=44   .and.  49 >=spgroupma) ) spgrmatch=0
 case(11)
   if (.not.(spgroupma>=52   .and.  57 >=spgroupma) ) spgrmatch=0
 case(12)
   if (.not.(spgroupma>=60   .and.  64 >=spgroupma) ) spgrmatch=0
 case(13)
   if (.not.(spgroupma>=67   .and.  74 >=spgroupma) ) spgrmatch=0
 case(14)
   if (.not.(spgroupma>=77   .and.  84 >=spgroupma) ) spgrmatch=0
 case(15)
   if (.not.(spgroupma>=87   .and.  91 >=spgroupma) ) spgrmatch=0
 case(16)
   if (.not.(spgroupma>= 3  .and.   6>=spgroupma) ) spgrmatch=0
 case(17)
   if (.not.(spgroupma>= 9  .and.  15 >=spgroupma) ) spgrmatch=0
 case(18)
   if (.not.(spgroupma>=18   .and.  24 >=spgroupma) ) spgrmatch=0
 case(19)
   if (.not.(spgroupma>=27   .and.  30 >=spgroupma) ) spgrmatch=0
 case(20)
   if (.not.(spgroupma>=33   .and.  37 >=spgroupma) ) spgrmatch=0
 case(21)
   if (.not.(spgroupma>=40   .and.  44 >=spgroupma) ) spgrmatch=0
 case(22)
   if (.not.(spgroupma>=47   .and.  48 >=spgroupma) ) spgrmatch=0
 case(23)
   if (.not.(spgroupma>=51   .and.  52 >=spgroupma) ) spgrmatch=0
 case(24)
   if (.not.(spgroupma>=55   .and.  56 >=spgroupma) ) spgrmatch=0
 case(25)
   if (.not.(spgroupma>=59   .and.  65 >=spgroupma) ) spgrmatch=0
 case(26)
   if (.not.(spgroupma>=68   .and.  77 >=spgroupma) ) spgrmatch=0
 case(27)
   if (.not.(spgroupma>=80   .and.  86 >=spgroupma) ) spgrmatch=0
 case(28)
   if (.not.(spgroupma>=89   .and.  98 >=spgroupma) ) spgrmatch=0
 case(29)
   if (.not.(spgroupma>=101  .and.  110 >=spgroupma) ) spgrmatch=0
 case(30)
   if (.not.(spgroupma>=113   .and. 122  >=spgroupma) ) spgrmatch=0
 case(31)
   if (.not.(spgroupma>=125   .and. 134  >=spgroupma) ) spgrmatch=0
 case(32)
   if (.not.(spgroupma>=137   .and. 143  >=spgroupma) ) spgrmatch=0
 case(33)
   if (.not.(spgroupma>=146   .and.  155 >=spgroupma) ) spgrmatch=0
 case(34)
   if (.not.(spgroupma>=158   .and.  164 >=spgroupma) ) spgrmatch=0
 case(35)
   if (.not.(spgroupma>=167   .and.   171>=spgroupma) ) spgrmatch=0
 case(36)
   if (.not.(spgroupma>=174   .and.   179>=spgroupma) ) spgrmatch=0
 case(37)
   if (.not.(spgroupma>=182   .and.  186 >=spgroupma) ) spgrmatch=0
 case(38)
   if (.not.(spgroupma>=189   .and.  194 >=spgroupma) ) spgrmatch=0
 case(39)
   if (.not.(spgroupma>=197   .and.  202 >=spgroupma) ) spgrmatch=0
 case(40)
   if (.not.(spgroupma>=205   .and.  210 >=spgroupma) ) spgrmatch=0
 case(41)
   if (.not.(spgroupma>=213   .and.  218 >=spgroupma) ) spgrmatch=0
 case(42)
   if (.not.(spgroupma>=221   .and.   223>=spgroupma) ) spgrmatch=0
 case(43)
   if (.not.(spgroupma>=226   .and.   228>=spgroupma) ) spgrmatch=0
 case(44)
   if (.not.(spgroupma>=231   .and.  234 >=spgroupma) ) spgrmatch=0
 case(45)
   if (.not.(spgroupma>=237   .and.  240 >=spgroupma) ) spgrmatch=0
 case(46)
   if (.not.(spgroupma>=243   .and.  248 >=spgroupma) ) spgrmatch=0
 case(47)
   if (.not.(spgroupma>=251   .and.  256 >=spgroupma) ) spgrmatch=0
 case(48)
   if (.not.(spgroupma>=259   .and.  264 >=spgroupma) ) spgrmatch=0
 case(49)
   if (.not.(spgroupma>=267   .and.  276 >=spgroupma) ) spgrmatch=0
 case(50)
   if (.not.(spgroupma>=279   .and.  288 >=spgroupma) ) spgrmatch=0
 case(51)
   if (.not.(spgroupma>=291   .and.  304 >=spgroupma) ) spgrmatch=0
 case(52)
   if (.not.(spgroupma>=307   .and.  320 >=spgroupma) ) spgrmatch=0
 case(53)
   if (.not.(spgroupma>=323   .and.  336 >=spgroupma) ) spgrmatch=0
 case(54)
   if (.not.(spgroupma>=339   .and.  352 >=spgroupma) ) spgrmatch=0
 case(55)
   if (.not.(spgroupma>=355   .and.  364 >=spgroupma) ) spgrmatch=0
 case(56)
   if (.not.(spgroupma>=367   .and.  376 >=spgroupma) ) spgrmatch=0
 case(57)
   if (.not.(spgroupma>=379   .and.  392 >=spgroupma) ) spgrmatch=0
 case(58)
   if (.not.(spgroupma>=395   .and.  404 >=spgroupma) ) spgrmatch=0
 case(59)
   if (.not.(spgroupma>=407   .and.  416 >=spgroupma) ) spgrmatch=0
 case(60)
   if (.not.(spgroupma>=419   .and.  432 >=spgroupma) ) spgrmatch=0
 case(61)
   if (.not.(spgroupma>=435   .and.  440 >=spgroupma) ) spgrmatch=0
 case(62)
   if (.not.(spgroupma>=443   .and.  456 >=spgroupma) ) spgrmatch=0
 case(63)
   if (.not.(spgroupma>=459   .and.  468 >=spgroupma) ) spgrmatch=0
 case(64)
   if (.not.(spgroupma>=471   .and.  480 >=spgroupma) ) spgrmatch=0
 case(65)
   if (.not.(spgroupma>=483   .and.  490 >=spgroupma) ) spgrmatch=0
 case(66)
   if (.not.(spgroupma>=493   .and.  500 >=spgroupma) ) spgrmatch=0
 case(67)
   if (.not.(spgroupma>=503   .and.  510 >=spgroupma) ) spgrmatch=0
 case(68)
   if (.not.(spgroupma>=513   .and.  520 >=spgroupma) ) spgrmatch=0
 case(69)
   if (.not.(spgroupma>=523   .and.  526 >=spgroupma) ) spgrmatch=0
 case(70)
   if (.not.(spgroupma>=529   .and.  532 >=spgroupma) ) spgrmatch=0
 case(71)
   if (.not.(spgroupma>=535   .and.  538 >=spgroupma) ) spgrmatch=0
 case(72)
   if (.not.(spgroupma>=541   .and.  547 >=spgroupma) ) spgrmatch=0
 case(73)
   if (.not.(spgroupma>=550   .and.  553 >=spgroupma) ) spgrmatch=0
 case(74)
   if (.not.(spgroupma>=556   .and.  562 >=spgroupma) ) spgrmatch=0
 case(75)
   if (.not.(spgroupma>= 3  .and.   6>=spgroupma) ) spgrmatch=0
 case(76)
   if (.not.(spgroupma>= 9  .and.  12 >=spgroupma) ) spgrmatch=0
 case(77)
   if (.not.(spgroupma>= 15  .and.   18>=spgroupma) ) spgrmatch=0
 case(78)
   if (.not.(spgroupma>= 21  .and.   24>=spgroupma) ) spgrmatch=0
 case(79)
   if (.not.(spgroupma>= 27  .and.  28 >=spgroupma) ) spgrmatch=0
 case(80)
   if (.not.(spgroupma>= 31  .and.  32 >=spgroupma) ) spgrmatch=0
 case(81)
   if (.not.(spgroupma>= 35  .and.  38 >=spgroupma) ) spgrmatch=0
 case(82)
   if (.not.(spgroupma>= 41  .and.  42 >=spgroupma) ) spgrmatch=0
 case(83)
   if (.not.(spgroupma>=45   .and.   50>=spgroupma) ) spgrmatch=0
 case(84)
   if (.not.(spgroupma>=53   .and.   58>=spgroupma) ) spgrmatch=0
 case(85)
   if (.not.(spgroupma>=61   .and.   66>=spgroupma) ) spgrmatch=0
 case(86)
   if (.not.(spgroupma>=69   .and.   74>=spgroupma) ) spgrmatch=0
 case(87)
   if (.not.(spgroupma>=77   .and.   80>=spgroupma) ) spgrmatch=0
 case(88)
   if (.not.(spgroupma>=83   .and.   86>=spgroupma) ) spgrmatch=0
 case(89)
   if (.not.(spgroupma>=89   .and.   94>=spgroupma) ) spgrmatch=0
 case(90)
   if (.not.(spgroupma>=97   .and.   102>=spgroupma) ) spgrmatch=0
 case(91)
   if (.not.(spgroupma>=105   .and.   110>=spgroupma) ) spgrmatch=0
 case(92)
   if (.not.(spgroupma>=113   .and.   118>=spgroupma) ) spgrmatch=0
 case(93)
   if (.not.(spgroupma>=121   .and.   126>=spgroupma) ) spgrmatch=0
 case(94)
   if (.not.(spgroupma>=129   .and.   134>=spgroupma) ) spgrmatch=0
 case(95)
   if (.not.(spgroupma>=137   .and.   142>=spgroupma) ) spgrmatch=0
 case(96)
   if (.not.(spgroupma>=145   .and.   150>=spgroupma) ) spgrmatch=0
 case(97)
   if (.not.(spgroupma>=153   .and.   156>=spgroupma) ) spgrmatch=0
 case(98)
   if (.not.(spgroupma>=159   .and.   162>=spgroupma) ) spgrmatch=0
 case(99)
   if (.not.(spgroupma>=165   .and.   170>=spgroupma) ) spgrmatch=0
 case(100)
   if (.not.(spgroupma>=173   .and.   178>=spgroupma) ) spgrmatch=0
 case(101)
   if (.not.(spgroupma>=181   .and.   186>=spgroupma) ) spgrmatch=0
 case(102)
   if (.not.(spgroupma>=189   .and.   194>=spgroupma) ) spgrmatch=0
 case(103)
   if (.not.(spgroupma>=197   .and.   202>=spgroupma) ) spgrmatch=0
 case(104)
   if (.not.(spgroupma>=205   .and.   210>=spgroupma) ) spgrmatch=0
 case(105)
   if (.not.(spgroupma>=213   .and.   218>=spgroupma) ) spgrmatch=0
 case(106)
   if (.not.(spgroupma>=221   .and.   226>=spgroupma) ) spgrmatch=0
 case(107)
   if (.not.(spgroupma>=229  .and.   232>=spgroupma) ) spgrmatch=0
 case(108)
   if (.not.(spgroupma>=235   .and.   238>=spgroupma) ) spgrmatch=0
 case(109)
   if (.not.(spgroupma>=241   .and.   244>=spgroupma) ) spgrmatch=0
 case(110)
   if (.not.(spgroupma>=247   .and.   250>=spgroupma) ) spgrmatch=0
 case(111)
   if (.not.(spgroupma>=253   .and.   258>=spgroupma) ) spgrmatch=0
 case(112)
   if (.not.(spgroupma>=261   .and.   266>=spgroupma) ) spgrmatch=0
 case(113)
   if (.not.(spgroupma>=269   .and.   274>=spgroupma) ) spgrmatch=0
 case(114)
   if (.not.(spgroupma>=277   .and.   282>=spgroupma) ) spgrmatch=0
 case(115)
   if (.not.(spgroupma>=285   .and.   290>=spgroupma) ) spgrmatch=0
 case(116)
   if (.not.(spgroupma>=293   .and.   298>=spgroupma) ) spgrmatch=0
 case(117)
   if (.not.(spgroupma>=301   .and.   306>=spgroupma) ) spgrmatch=0
 case(118)
   if (.not.(spgroupma>=309   .and.   314>=spgroupma) ) spgrmatch=0
 case(119)
   if (.not.(spgroupma>=317   .and.   320>=spgroupma) ) spgrmatch=0
 case(120)
   if (.not.(spgroupma>=323   .and.   326>=spgroupma) ) spgrmatch=0
 case(121)
   if (.not.(spgroupma>=329   .and.   332>=spgroupma) ) spgrmatch=0
 case(122)
   if (.not.(spgroupma>=335   .and.   338>=spgroupma) ) spgrmatch=0
 case(123)
   if (.not.(spgroupma>=341   .and.   350>=spgroupma) ) spgrmatch=0
 case(124)
   if (.not.(spgroupma>=353   .and.   362>=spgroupma) ) spgrmatch=0
 case(125)
   if (.not.(spgroupma>=365   .and.   374>=spgroupma) ) spgrmatch=0
 case(126)
   if (.not.(spgroupma>=377   .and.   386>=spgroupma) ) spgrmatch=0
 case(127)
   if (.not.(spgroupma>=389   .and.   398>=spgroupma) ) spgrmatch=0
 case(128)
   if (.not.(spgroupma>=401   .and.   410>=spgroupma) ) spgrmatch=0
 case(129)
   if (.not.(spgroupma>=413   .and.   422>=spgroupma) ) spgrmatch=0
 case(130)
   if (.not.(spgroupma>=425   .and.   434>=spgroupma) ) spgrmatch=0
 case(131)
   if (.not.(spgroupma>=437   .and.   446>=spgroupma) ) spgrmatch=0
 case(132)
   if (.not.(spgroupma>=449   .and.   458>=spgroupma) ) spgrmatch=0
 case(133)
   if (.not.(spgroupma>=461   .and.   470>=spgroupma) ) spgrmatch=0
 case(134)
   if (.not.(spgroupma>=473   .and.   482>=spgroupma) ) spgrmatch=0
 case(135)
   if (.not.(spgroupma>=485   .and.   494>=spgroupma) ) spgrmatch=0
 case(136)
   if (.not.(spgroupma>=497  .and.   506>=spgroupma) ) spgrmatch=0
 case(137)
   if (.not.(spgroupma>=509   .and.   518>=spgroupma) ) spgrmatch=0
 case(138)
   if (.not.(spgroupma>=521   .and.   530>=spgroupma) ) spgrmatch=0
 case(139)
   if (.not.(spgroupma>=533   .and.   540>=spgroupma) ) spgrmatch=0
 case(140)
   if (.not.(spgroupma>=543   .and.   550>=spgroupma) ) spgrmatch=0
 case(141)
   if (.not.(spgroupma>=553   .and.   560>=spgroupma) ) spgrmatch=0
 case(142)
   if (.not.(spgroupma>=563   .and.   570>=spgroupma) ) spgrmatch=0
 case(143)
   if (.not.(spgroupma>=3   .and.   3>=spgroupma) ) spgrmatch=0
 case(144)
   if (.not.(spgroupma>=6   .and.   6>=spgroupma) ) spgrmatch=0
 case(145)
   if (.not.(spgroupma>=9   .and.   9>=spgroupma) ) spgrmatch=0
 case(146)
   if (.not.(spgroupma>=12   .and.   12>=spgroupma) ) spgrmatch=0
 case(147)
   if (.not.(spgroupma>=15   .and.   16>=spgroupma) ) spgrmatch=0
 case(148)
   if (.not.(spgroupma>=19   .and.   20>=spgroupma) ) spgrmatch=0
 case(149)
   if (.not.(spgroupma>=23   .and.   24>=spgroupma) ) spgrmatch=0
 case(150)
   if (.not.(spgroupma>=27   .and.   28>=spgroupma) ) spgrmatch=0
 case(151)
   if (.not.(spgroupma>=31   .and.   32>=spgroupma) ) spgrmatch=0
 case(152)
   if (.not.(spgroupma>=35   .and.   36>=spgroupma) ) spgrmatch=0
 case(153)
   if (.not.(spgroupma>=39   .and.   40>=spgroupma) ) spgrmatch=0
 case(154)
   if (.not.(spgroupma>=43   .and.   44>=spgroupma) ) spgrmatch=0
 case(155)
   if (.not.(spgroupma>=47   .and.   48>=spgroupma) ) spgrmatch=0
 case(156)
   if (.not.(spgroupma>=51   .and.   52>=spgroupma) ) spgrmatch=0
 case(157)
   if (.not.(spgroupma>=55   .and.   56>=spgroupma) ) spgrmatch=0
 case(158)
   if (.not.(spgroupma>=59   .and.   60>=spgroupma) ) spgrmatch=0
 case(159)
   if (.not.(spgroupma>=63   .and.   64>=spgroupma) ) spgrmatch=0
 case(160)
   if (.not.(spgroupma>=67   .and.   68>=spgroupma) ) spgrmatch=0
 case(161)
   if (.not.(spgroupma>=71   .and.   72>=spgroupma) ) spgrmatch=0
 case(162)
   if (.not.(spgroupma>=75   .and.   78>=spgroupma) ) spgrmatch=0
 case(163)
   if (.not.(spgroupma>=81   .and.   84>=spgroupma) ) spgrmatch=0
 case(164)
   if (.not.(spgroupma>=87   .and.   90>=spgroupma) ) spgrmatch=0
 case(165)
   if (.not.(spgroupma>=93   .and.   96>=spgroupma) ) spgrmatch=0
 case(166)
   if (.not.(spgroupma>=99   .and.   102>=spgroupma) ) spgrmatch=0
 case(167)
   if (.not.(spgroupma>=105   .and.   108>=spgroupma) ) spgrmatch=0
 case(168)
   if (.not.(spgroupma>=111   .and.   112>=spgroupma) ) spgrmatch=0
 case(169)
   if (.not.(spgroupma>=115   .and.   116>=spgroupma) ) spgrmatch=0
 case(170)
   if (.not.(spgroupma>=119   .and.   120>=spgroupma) ) spgrmatch=0
 case(171)
   if (.not.(spgroupma>=123   .and.   124>=spgroupma) ) spgrmatch=0
 case(172)
   if (.not.(spgroupma>=127   .and.   128>=spgroupma) ) spgrmatch=0
 case(173)
   if (.not.(spgroupma>=131   .and.   132>=spgroupma) ) spgrmatch=0
 case(174)
   if (.not.(spgroupma>=135   .and.   136>=spgroupma) ) spgrmatch=0
 case(175)
   if (.not.(spgroupma>=139   .and.   142>=spgroupma) ) spgrmatch=0
 case(176)
   if (.not.(spgroupma>=145   .and.   148>=spgroupma) ) spgrmatch=0
 case(177)
   if (.not.(spgroupma>=151   .and.   154>=spgroupma) ) spgrmatch=0
 case(178)
   if (.not.(spgroupma>=157   .and.   160>=spgroupma) ) spgrmatch=0
 case(179)
   if (.not.(spgroupma>=163   .and.   166>=spgroupma) ) spgrmatch=0
 case(180)
   if (.not.(spgroupma>=169   .and.   172>=spgroupma) ) spgrmatch=0
 case(181)
   if (.not.(spgroupma>=175   .and.   178>=spgroupma) ) spgrmatch=0
 case(182)
   if (.not.(spgroupma>=181   .and.   184>=spgroupma) ) spgrmatch=0
 case(183)
   if (.not.(spgroupma>=187   .and.   190>=spgroupma) ) spgrmatch=0
 case(184)
   if (.not.(spgroupma>=193   .and.   196>=spgroupma) ) spgrmatch=0
 case(185)
   if (.not.(spgroupma>=199   .and.   202>=spgroupma) ) spgrmatch=0
 case(186)
   if (.not.(spgroupma>=205   .and.   208>=spgroupma) ) spgrmatch=0
 case(187)
   if (.not.(spgroupma>=211   .and.   214>=spgroupma) ) spgrmatch=0
 case(188)
   if (.not.(spgroupma>=217   .and.   220>=spgroupma) ) spgrmatch=0
 case(189)
   if (.not.(spgroupma>=223   .and.   226>=spgroupma) ) spgrmatch=0
 case(190)
   if (.not.(spgroupma>=229   .and.   232>=spgroupma) ) spgrmatch=0
 case(191)
   if (.not.(spgroupma>=235   .and.   242>=spgroupma) ) spgrmatch=0
 case(192)
   if (.not.(spgroupma>=245   .and.   252>=spgroupma) ) spgrmatch=0
 case(193)
   if (.not.(spgroupma>=255   .and.   262>=spgroupma) ) spgrmatch=0
 case(194)
   if (.not.(spgroupma>=265   .and.   272>=spgroupma) ) spgrmatch=0
 case(195)
   if (.not.(spgroupma>=3   .and.   3>=spgroupma) ) spgrmatch=0
 case(196)
   if (.not.(spgroupma>=6   .and.   6>=spgroupma) ) spgrmatch=0
 case(197)
   spgrmatch=0
 case(198)
   if (.not.(spgroupma>=11  .and.   11>=spgroupma) ) spgrmatch=0
 case(199)
   spgrmatch=0
 case(200)
   if (.not.(spgroupma>=16   .and.   17>=spgroupma) ) spgrmatch=0
 case(201)
   if (.not.(spgroupma>=20   .and.   21>=spgroupma) ) spgrmatch=0
 case(202)
   if (.not.(spgroupma>=24   .and.   25>=spgroupma) ) spgrmatch=0
 case(203)
   if (.not.(spgroupma>=28   .and.   29>=spgroupma) ) spgrmatch=0
 case(204)
   if (.not.(spgroupma>=32   .and.   32>=spgroupma) ) spgrmatch=0
 case(205)
   if (.not.(spgroupma>=35   .and.   36>=spgroupma) ) spgrmatch=0
 case(206)
   if (.not.(spgroupma>=39   .and.   39>=spgroupma) ) spgrmatch=0
 case(207)
   if (.not.(spgroupma>=42   .and.   43>=spgroupma) ) spgrmatch=0
 case(208)
   if (.not.(spgroupma>=46   .and.   47>=spgroupma) ) spgrmatch=0
 case(209)
   if (.not.(spgroupma>=50   .and.   51>=spgroupma) ) spgrmatch=0
 case(210)
   if (.not.(spgroupma>=54   .and.   55>=spgroupma) ) spgrmatch=0
 case(211)
   if (.not.(spgroupma>=58   .and.   58>=spgroupma) ) spgrmatch=0
 case(212)
   if (.not.(spgroupma>=61  .and.   62>=spgroupma) ) spgrmatch=0
 case(213)
   if (.not.(spgroupma>=65   .and.   66>=spgroupma) ) spgrmatch=0
 case(214)
   if (.not.(spgroupma>=69   .and.   69>=spgroupma) ) spgrmatch=0
 case(215)
   if (.not.(spgroupma>=72   .and.   73>=spgroupma) ) spgrmatch=0
 case(216)
   if (.not.(spgroupma>=76   .and.   77>=spgroupma) ) spgrmatch=0
 case(217)
   if (.not.(spgroupma>=80   .and.   80>=spgroupma) ) spgrmatch=0
 case(218)
   if (.not.(spgroupma>=83   .and.   84>=spgroupma) ) spgrmatch=0
 case(219)
   if (.not.(spgroupma>=87   .and.   88>=spgroupma) ) spgrmatch=0
 case(220)
   if (.not.(spgroupma>=91   .and.   91>=spgroupma) ) spgrmatch=0
 case(221)
   if (.not.(spgroupma>=94   .and.   97>=spgroupma) ) spgrmatch=0
 case(222)
   if (.not.(spgroupma>=100  .and.   103>=spgroupma) ) spgrmatch=0
 case(223)
   if (.not.(spgroupma>=106   .and.   109>=spgroupma) ) spgrmatch=0
 case(224)
   if (.not.(spgroupma>=112   .and.   115>=spgroupma) ) spgrmatch=0
 case(225)
   if (.not.(spgroupma>=118   .and.   121>=spgroupma) ) spgrmatch=0
 case(226)
   if (.not.(spgroupma>=124   .and.   127>=spgroupma) ) spgrmatch=0
 case(227)
   if (.not.(spgroupma>=130   .and.   133>=spgroupma) ) spgrmatch=0
 case(228)
   if (.not.(spgroupma>=136   .and.   139>=spgroupma) ) spgrmatch=0
 case(229)
   if (.not.(spgroupma>=142   .and.   144>=spgroupma) ) spgrmatch=0
 case(230)
   if (.not.(spgroupma>=147   .and.   149>=spgroupma) ) spgrmatch=0
 case default
   write(message, '(3a,i8,4a)' )&
&   'The non-magnetic spacegroup is not specified ',ch10,&
&   'while the magnetic space group is specified, spgroupma= ',spgroupma,ch10,&
&   'This is not allowed.  ',ch10,&
&   'Action: specify spgroup in the input file.'
   MSG_ERROR(message)
 end select

 if (spgrmatch==0) then
   write(message, '(a,i8,a,a,i8,4a)' )&
&   'mismatch between the non-magnetic spacegroup ',spgroup,ch10,&
&   'and the magnetic space group ',spgroupma,ch10,&
&   'This is not allowed.  ',ch10,&
&   'Action: modify spgroup or spgroupma in the input file.'
   MSG_ERROR(message)
 end if

!DEBUG
!write(std_out,*) ' gensymshub, after check the spgroup ... '
!ENDDEBUG

!Assign the magnetic Bravais lattice type from the magnetic space group number
!As the magnetic space group number begins with 1 for EACH crystal system
!we must first make our choice as a function of spgroup

!Convention :
!brvlttbw = input variable giving Bravais black-and-white translation
!(from 1 to 8 : Shubnikov type IV space group)
!1,2,3 = 1/2 translation along a, b, or c respectively
!4,5,6 = translation corresponding to A,B,C centering, respectively
!7 = I centering
!8 = (1/2 0 0) centering corresponding to normal F lattice
!9 -> Shubnikov type III space group
!Note that the use of the table 7.3 (p585) of Bradney and Cracknell
!for the definition of the translation vectors
!is extremely confusing, due to the strange choice of basis
!vectors of table 3.1.
!See table 7.4 (p588) for the spgroupma interpretation

 select case(spgroup)
 case(1,2)     !Triclinic
   select case(spgroupma)
   case(6)
     brvlttbw=9      !ShubIII
   case(3,7)
     brvlttbw=1      !Ps (note that it is not body centered, according to Table
!        7.3 of Bradley and Cracknell)
   end select
 case(3:15)     !Monoclinic
   select case(spgroupma)
   case(3,9,15,20,26,34,39,44:46,52:54,60:62,67:69,77:79,87:89)
     brvlttbw=9      !ShubIII
   case(4,10,17,21,27,36,41,47,55,64,70,80,91)
     brvlttbw=1     !a
   case(5,11,22,29,48,56,71,81)
     brvlttbw=2     !b
   case(16,28,35,40,63,72,82,90)
     brvlttbw=3     !c
   case(31,73,83)
     brvlttbw=4     !A
   case(6,12,23,30,49,57,74,84)
     brvlttbw=6     !C
   end select
 case(16:74)     !Orthorhombic
   select case(spgroupma)
   case(3,9,10,18,19,27,33,34,40,41,47,51,55,59,60,68:70,&
&     80,81,89:91,101:103,113:115,125:127,137,138,146:148,158,159,&
&     167,168,174:176,182,183,189:191,197:199,205:207,213:215,221,&
&     222,226,227,231,232,237,238,243:245,251:253,259:261,267:271,&
&     279:283,291:297,307:313,323:329,339:345,355:359,367:371,&
&     379:385,395:399,407:411,419:425,435:437,443:449,459:465,&
&     471:477,483:487,493:497,503:507,513:517,523:525,529:531,&
&     535:537,541:545,550:552,556:560)
     brvlttbw=9      !ShubIII
   case(4,11,20,28,36,43,62,71,83,92,104,116,128,140,&
&     149,160,170,178,185,192,200,208,216,234,240,247,254,262,272,&
&     284,298,314,330,346,360,372,386,400,412,426,438,450,467,&
&     479,489,499,509,519,547,562)
     brvlttbw=1     !a
   case(72,93,105,117,129,150,248,299,315,331,347,387,427,451)
     brvlttbw=2     !b
   case(12,21,35,42,52,56,61,73,82,94,106,118,130,139,&
&     151,161,169,177,184,193,201,209,217,233,239,246,273,285,300,&
&     316,332,348,361,373,388,401,413,428,452,466,478,488,498,&
&     508,518,538,546,553,561)
     brvlttbw=3     !c
   case(13,22,37,44,64,74,85,95,107,119,131,142,152,162,&
&     171,179,186,274,286,301,317,333,349,362,374,389,402,414,429,&
&     453,468,480,490,500,510,520)
     brvlttbw=4     !A
   case(75,96,108,120,132,153,302,318,334,350,390,430,454)
     brvlttbw=5     !B
   case(5,14,23,29,63,76,84,97,109,121,133,141,154,163,&
&     194,202,210,218,255,263,275,287,303,319,335,351,363,375,&
&     391,403,415,431,439,455)
     brvlttbw=6     !C
   case(6,15,24,30,65,77,86,98,110,122,134,143,155,164,256,&
&     264,276,288,304,320,336,352,364,376,392,404,416,432,440,456)
     brvlttbw=7     !I
   case(48,223,228,526,532)
     brvlttbw=8     !s
   end select
 case(75:142)       !Tetragonal
   select case(spgroupma)
   case(3,9,15,21,27,31,35,41,45:47,53:55,61:63,69:71,&
&     77:79,83:85,89:91,97:99,105:107,113:115,121:123,129:131,&
&     137:139,145:147,153:155,159:161,165:167,173:175,181:183,&
&     189:191,197:199,205:207,213:215,221:223,229:231,235:237,&
&     241:243,247:249,253:255,261:263,269:271,277:279,285:287,&
&     293:295,301:303,309:311,317:319,323:325,329:331,335:337,&
&     341:347,353:359,365:371,377:383,389:395,401:407,413:419,&
&     425:431,437:443,449:455,461:467,473:479,485:491,497:503,&
&     509:515,521:527,533:539,543:549,553:559,563:569)
     brvlttbw=9      !ShubIII
   case(4,10,16,22,28,32,36,42,48,56,64,72,80,86,92,100,&
&     108,116,124,132,140,148,156,162,168,176,184,192,200,208,&
&     216,224,232,238,244,250,256,264,272,280,288,296,304,312,&
&     320,326,332,338,348,360,372,384,396,408,420,432,444,456,&
&     468,480,492,504,516,528,540,550,560,570)
     brvlttbw=3     !c
   case(5,11,17,23,37,49,57,65,73,93,101,109,117,125,133,141,&
&     149,169,177,185,193,201,209,217,225,257,265,273,281,289,297,&
&     305,313,349,361,373,385,397,409,421,433,445,457,469,481,493,&
&     505,517,529)
     brvlttbw=6     !C
   case(6,12,18,24,38,50,58,66,74,94,102,110,118,126,134,142,&
&     150,170,178,186,194,202,210,218,226,258,266,274,282,290,298,&
&     306,314,350,362,374,386,398,410,422,434,446,458,470,482,494,&
&     506,518,530)
     brvlttbw=7     !I
   end select
 case(143:194)      !Hexagonal
   select case(spgroupma)
   case(15,19,23,27,31,35,39,43,47,51,55,59,63,67,71,&
&     75:77,81:83,87:89,93:95,99:101,105:107,111,115,119,123,127,&
&     131,135,139:141,145:147,151:153,157:159,163:165,169:171,&
&     175:177,181:183,187:189,193:195,199:201,205:207,211:213,&
&     217:219,223:225,229:231,235:241,245:251,255:261,265:271)
     brvlttbw=9      !ShubIII
   case(3,6,9,16,24,28,32,36,40,44,52,56,60,64,78,84,&
&     90,96,112,116,120,124,128,132,136,142,148,154,160,166,172,&
&     178,184,190,196,202,208,214,220,226,232,242,252,262,272)
     brvlttbw=6     !C
   case(12,20,48,68,72,102,108)
     brvlttbw=7     !I
   end select
 case(195:230)       !Cubic
   select case(spgroupma)
   case(16,20,24,28,32,35,39,42,46,50,54,58,61,65,69,&
&     72,76,80,83,87,91,94:96,100:102,106:108,112:114,118:120,&
&     124:126,130:132,136:138,142:144,147:149)
     brvlttbw=9      !ShubIII
   case(3,11,17,21,36,43,47,62,66,73,84,97,103,109,115)
     brvlttbw=7     !I
   case(6,25,29,51,55,77,88,121,127,133,139)
     brvlttbw=8     !s
   end select
 end select

 if(brvlttbw==9)shubnikov=3 !Shubnikov type III
 if(brvlttbw>=1 .and. brvlttbw<=8) shubnikov=4 !Shubnikov type IV

 genafm(:)=zero
 if(shubnikov==4)then
   if(brvlttbw==1)genafm(1)=half
   if(brvlttbw==2)genafm(2)=half
   if(brvlttbw==3)genafm(3)=half
   if(brvlttbw>=4 .and. brvlttbw<=8)then
     genafm(:)=half
     if(brvlttbw==4)genafm(1)=zero
     if(brvlttbw==5)genafm(2)=zero
     if(brvlttbw==6)genafm(3)=zero
   end if
 end if

!DEBUG
!write(std_out,*) 'gensymshub : end '
!write(std_out,*) 'gensymshub, shubnikov =',shubnikov
!ENDDEBUG
end subroutine gensymshub


!!***
