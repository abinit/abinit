!{\src2tex{textfont=tt}}
!!****f* ABINIT/int_ang
!! NAME
!! int_ang
!!
!! FUNCTION
!! Evaluate angular part for <phi_i|nabla|phi_j> and <tphi_i|nabla|tphi_j>
!!
!! COPYRIGHT
!! Copyright (C) 2005-2018 ABINIT group (VR,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  mpsang=1+ max. angular momentum
!!
!! OUTPUT
!!  ang_phipphj :: angular part for <phi_i|nabla|phi_j> and <tphi_i|nabla|tphi_j>
!!  ang_phipphj(i,j,1)=\int sin\theta cos\phi Si Sj d\omega
!!  ang_phipphj(i,j,2)=\int cos\theta cos\phi Si \frac{d}{d\theta}Sj d\Omega
!!  ang_phipphj(i,j,3)=\int -sin\phi  Si \frac{d}{d\phi}Sj d\Omega
!!  ang_phipphj(i,j,4)=\int sin\theta sin\phi Si Sj d\Omega
!!  ang_phipphj(i,j,5)=\int cos\theta sin\phi Si \frac{d}{d\theta}Sj d\Omega
!!  ang_phipphj(i,j,6)=\int cos\phi Si \frac{d}{d\phi}Sj d\Omega
!!  ang_phipphj(i,j,7)=\int cos\theta  Si Sj d\Omega
!!  ang_phipphj(i,j,8)=\int -sin\theta Si \frac{d}{d\theta}Sj d\Omega
!!
!! PARENTS
!!      optics_paw,optics_paw_core,pawnabla_init
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine int_ang(ang_phipphj,mpsang)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'int_ang'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpsang
!arrays
 real(dp),intent(out) :: ang_phipphj(mpsang**2,mpsang**2,8)

!Local variables-------------------------------
 character(len=500) :: message
 real(dp) :: ang_phipphj_tmp(16,16,8)

! ************************************************************************

 if (mpsang>4) then
   message = '  Not designed for angular momentum greater than 3 !'
   MSG_ERROR(message)
 end if

!8 angular integrals for l=0..3, m=-l..+l
!ang_phipphj(1,4,1)=\frac{1}{\sqrt{3}}
!ang_phipphj(2,5,1)=\frac{1}{\sqrt{5}}
!ang_phipphj(3,8,1)=\frac{1}{\sqrt{5}}
!ang_phipphj(4,1,1)=\frac{1}{\sqrt{3}}
!ang_phipphj(4,7,1)=-\frac{1}{\sqrt{15}}
!ang_phipphj(4,9,1)=\frac{1}{\sqrt{5}}
!ang_phipphj(5,2,1)=\frac{1}{\sqrt{5}}
!ang_phipphj(5,10,1)=\sqrt{\frac{3}{14}}
!ang_phipphj(5,12,1)=-\frac{1}{\sqrt{70}}
!ang_phipphj(6,11,1)=\frac{1}{\sqrt{7}}
!ang_phipphj(7,4,1)=-\frac{1}{\sqrt{15}}
!ang_phipphj(7,14,1)=\sqrt{\frac{6}{35}}
!ang_phipphj(8,3,1)=\frac{1}{\sqrt{5}}
!ang_phipphj(8,13,1)=-\sqrt{\frac{3}{35}}
!ang_phipphj(8,15,1)=\frac{1}{\sqrt{7}}
!ang_phipphj(9,4,1)=\frac{1}{\sqrt{5}}
!ang_phipphj(9,14,1)=-\frac{1}{\sqrt{70}}
!ang_phipphj(9,16,1)=\sqrt{\frac{3}{14}}
!ang_phipphj(10,5,1)=\sqrt{\frac{3}{14}}
!ang_phipphj(11,6,1)=\frac{1}{\sqrt{7}}
!ang_phipphj(12,5,1)=-\frac{1}{\sqrt{70}}
!ang_phipphj(13,8,1)=-\sqrt{\frac{3}{35}}
!ang_phipphj(14,7,1)=\sqrt{\frac{6}{35}}
!ang_phipphj(14,9,1)=-\frac{1}{\sqrt{70}}
!ang_phipphj(15,8,1)=\frac{1}{\sqrt{7}}
!ang_phipphj(16,9,1)=\sqrt{\frac{3}{14}}
!ang_phipphj(1,4,2)=\frac{1}{2 \sqrt{3}}
!ang_phipphj(1,14,2)=-\frac{\sqrt{\frac{7}{6}}}{2}
!ang_phipphj(2,5,2)=\frac{1}{2 \sqrt{5}}
!ang_phipphj(3,8,2)=\frac{1}{2 \sqrt{5}}
!ang_phipphj(4,7,2)=-\sqrt{\frac{3}{5}}
!ang_phipphj(4,9,2)=\frac{1}{2 \sqrt{5}}
!ang_phipphj(5,2,2)=\frac{1}{4 \sqrt{5}}
!ang_phipphj(5,10,2)=\frac{\sqrt{\frac{3}{14}}}{2}
!ang_phipphj(5,12,2)=-2 \sqrt{\frac{2}{35}}
!ang_phipphj(6,11,2)=\frac{1}{2 \sqrt{7}}
!ang_phipphj(7,4,2)=\frac{1}{\sqrt{15}}
!ang_phipphj(7,14,2)=\frac{13}{2 \sqrt{210}}
!ang_phipphj(8,3,2)=-\frac{1}{\sqrt{5}}
!ang_phipphj(8,13,2)=-4 \sqrt{\frac{3}{35}}
!ang_phipphj(8,15,2)=\frac{1}{2 \sqrt{7}}
!ang_phipphj(9,4,2)=\frac{1}{4 \sqrt{5}}
!ang_phipphj(9,14,2)=-2 \sqrt{\frac{2}{35}}
!ang_phipphj(9,16,2)=\frac{\sqrt{\frac{3}{14}}}{2}
!ang_phipphj(10,5,2)=\frac{1}{\sqrt{42}}
!ang_phipphj(11,6,2)=-\frac{1}{4 \sqrt{7}}
!ang_phipphj(12,5,2)=\sqrt{\frac{2}{35}}
!ang_phipphj(13,8,2)=2 \sqrt{\frac{3}{35}}
!ang_phipphj(14,7,2)=-2 \sqrt{\frac{6}{35}}
!ang_phipphj(14,9,2)=\sqrt{\frac{2}{35}}
!ang_phipphj(15,8,2)=-\frac{1}{4 \sqrt{7}}
!ang_phipphj(16,9,2)=\frac{1}{\sqrt{42}}
!ang_phipphj(1,4,3)=\frac{\sqrt{3}}{2}
!ang_phipphj(1,14,3)=\frac{\sqrt{\frac{7}{6}}}{2}
!ang_phipphj(2,5,3)=\frac{\sqrt{5}}{2}
!ang_phipphj(3,8,3)=\frac{\sqrt{5}}{2}
!ang_phipphj(4,9,3)=\frac{\sqrt{5}}{2}
!ang_phipphj(5,2,3)=-\frac{\sqrt{5}}{4}
!ang_phipphj(5,10,3)=\frac{\sqrt{\frac{21}{2}}}{2}
!ang_phipphj(6,11,3)=\frac{\sqrt{7}}{2}
!ang_phipphj(7,14,3)=\frac{\sqrt{\frac{35}{6}}}{2}
!ang_phipphj(8,15,3)=\frac{\sqrt{7}}{2}
!ang_phipphj(9,4,3)=-\frac{\sqrt{5}}{4}
!ang_phipphj(9,16,3)=\frac{\sqrt{\frac{21}{2}}}{2}
!ang_phipphj(10,5,3)=-\sqrt{\frac{7}{6}}
!ang_phipphj(11,6,3)=-\frac{\sqrt{7}}{4}
!ang_phipphj(15,8,3)=-\frac{\sqrt{7}}{4}
!ang_phipphj(16,9,3)=-\sqrt{\frac{7}{6}}
!ang_phipphj(1,2,4)=\frac{1}{\sqrt{3}}
!ang_phipphj(2,1,4)=\frac{1}{\sqrt{3}}
!ang_phipphj(2,7,4)=-\frac{1}{\sqrt{15}}
!ang_phipphj(2,9,4)=-\frac{1}{\sqrt{5}}
!ang_phipphj(3,6,4)=\frac{1}{\sqrt{5}}
!ang_phipphj(4,5,4)=\frac{1}{\sqrt{5}}
!ang_phipphj(5,4,4)=\frac{1}{\sqrt{5}}
!ang_phipphj(5,14,4)=-\frac{1}{\sqrt{70}}
!ang_phipphj(5,16,4)=-\sqrt{\frac{3}{14}}
!ang_phipphj(6,3,4)=\frac{1}{\sqrt{5}}
!ang_phipphj(6,13,4)=-\sqrt{\frac{3}{35}}
!ang_phipphj(6,15,4)=-\frac{1}{\sqrt{7}}
!ang_phipphj(7,2,4)=-\frac{1}{\sqrt{15}}
!ang_phipphj(7,12,4)=\sqrt{\frac{6}{35}}
!ang_phipphj(8,11,4)=\frac{1}{\sqrt{7}}
!ang_phipphj(9,2,4)=-\frac{1}{\sqrt{5}}
!ang_phipphj(9,10,4)=\sqrt{\frac{3}{14}}
!ang_phipphj(9,12,4)=\frac{1}{\sqrt{70}}
!ang_phipphj(10,9,4)=\sqrt{\frac{3}{14}}
!ang_phipphj(11,8,4)=\frac{1}{\sqrt{7}}
!ang_phipphj(12,7,4)=\sqrt{\frac{6}{35}}
!ang_phipphj(12,9,4)=\frac{1}{\sqrt{70}}
!ang_phipphj(13,6,4)=-\sqrt{\frac{3}{35}}
!ang_phipphj(14,5,4)=-\frac{1}{\sqrt{70}}
!ang_phipphj(15,6,4)=-\frac{1}{\sqrt{7}}
!ang_phipphj(16,5,4)=-\sqrt{\frac{3}{14}}
!ang_phipphj(1,2,5)=\frac{1}{2 \sqrt{3}}
!ang_phipphj(1,12,5)=-\frac{\sqrt{\frac{7}{6}}}{2}
!ang_phipphj(2,7,5)=-\sqrt{\frac{3}{5}}
!ang_phipphj(2,9,5)=-\frac{1}{2 \sqrt{5}}
!ang_phipphj(3,6,5)=\frac{1}{2 \sqrt{5}}
!ang_phipphj(4,5,5)=\frac{1}{2 \sqrt{5}}
!ang_phipphj(5,4,5)=\frac{1}{4 \sqrt{5}}
!ang_phipphj(5,14,5)=-2 \sqrt{\frac{2}{35}}
!ang_phipphj(5,16,5)=-\frac{\sqrt{\frac{3}{14}}}{2}
!ang_phipphj(6,3,5)=-\frac{1}{\sqrt{5}}
!ang_phipphj(6,13,5)=-4 \sqrt{\frac{3}{35}}
!ang_phipphj(6,15,5)=-\frac{1}{2 \sqrt{7}}
!ang_phipphj(7,2,5)=\frac{1}{\sqrt{15}}
!ang_phipphj(7,12,5)=\frac{13}{2 \sqrt{210}}
!ang_phipphj(8,11,5)=\frac{1}{2 \sqrt{7}}
!ang_phipphj(9,2,5)=-\frac{1}{4 \sqrt{5}}
!ang_phipphj(9,10,5)=\frac{\sqrt{\frac{3}{14}}}{2}
!ang_phipphj(9,12,5)=2 \sqrt{\frac{2}{35}}
!ang_phipphj(10,9,5)=\frac{1}{\sqrt{42}}
!ang_phipphj(11,8,5)=-\frac{1}{4 \sqrt{7}}
!ang_phipphj(12,7,5)=-2 \sqrt{\frac{6}{35}}
!ang_phipphj(12,9,5)=-\sqrt{\frac{2}{35}}
!ang_phipphj(13,6,5)=2 \sqrt{\frac{3}{35}}
!ang_phipphj(14,5,5)=\sqrt{\frac{2}{35}}
!ang_phipphj(15,6,5)=\frac{1}{4 \sqrt{7}}
!ang_phipphj(16,5,5)=-\frac{1}{\sqrt{42}}
!ang_phipphj(1,2,6)=\frac{\sqrt{3}}{2}
!ang_phipphj(1,12,6)=\frac{\sqrt{\frac{7}{6}}}{2}
!ang_phipphj(2,9,6)=-\frac{\sqrt{5}}{2}
!ang_phipphj(3,6,6)=\frac{\sqrt{5}}{2}
!ang_phipphj(4,5,6)=\frac{\sqrt{5}}{2}
!ang_phipphj(5,4,6)=-\frac{\sqrt{5}}{4}
!ang_phipphj(5,16,6)=-\frac{\sqrt{\frac{21}{2}}}{2}
!ang_phipphj(6,15,6)=-\frac{\sqrt{7}}{2}
!ang_phipphj(7,12,6)=\frac{\sqrt{\frac{35}{6}}}{2}
!ang_phipphj(8,11,6)=\frac{\sqrt{7}}{2}
!ang_phipphj(9,2,6)=\frac{\sqrt{5}}{4}
!ang_phipphj(9,10,6)=\frac{\sqrt{\frac{21}{2}}}{2}
!ang_phipphj(10,9,6)=-\sqrt{\frac{7}{6}}
!ang_phipphj(11,8,6)=-\frac{\sqrt{7}}{4}
!ang_phipphj(15,6,6)=\frac{\sqrt{7}}{4}
!ang_phipphj(16,5,6)=\sqrt{\frac{7}{6}}
!ang_phipphj(1,3,7)=\frac{1}{\sqrt{3}}
!ang_phipphj(2,6,7)=\frac{1}{\sqrt{5}}
!ang_phipphj(3,1,7)=\frac{1}{\sqrt{3}}
!ang_phipphj(3,7,7)=\frac{2}{\sqrt{15}}
!ang_phipphj(4,8,7)=\frac{1}{\sqrt{5}}
!ang_phipphj(5,11,7)=\frac{1}{\sqrt{7}}
!ang_phipphj(6,2,7)=\frac{1}{\sqrt{5}}
!ang_phipphj(6,12,7)=2 \sqrt{\frac{2}{35}}
!ang_phipphj(7,3,7)=\frac{2}{\sqrt{15}}
!ang_phipphj(7,13,7)=\frac{3}{\sqrt{35}}
!ang_phipphj(8,4,7)=\frac{1}{\sqrt{5}}
!ang_phipphj(8,14,7)=2 \sqrt{\frac{2}{35}}
!ang_phipphj(9,15,7)=\frac{1}{\sqrt{7}}
!ang_phipphj(11,5,7)=\frac{1}{\sqrt{7}}
!ang_phipphj(12,6,7)=2 \sqrt{\frac{2}{35}}
!ang_phipphj(13,7,7)=\frac{3}{\sqrt{35}}
!ang_phipphj(14,8,7)=2 \sqrt{\frac{2}{35}}
!ang_phipphj(15,9,7)=\frac{1}{\sqrt{7}}
!ang_phipphj(1,3,8)=\frac{2}{\sqrt{3}}
!ang_phipphj(2,6,8)=\frac{3}{\sqrt{5}}
!ang_phipphj(3,7,8)=2 \sqrt{\frac{3}{5}}
!ang_phipphj(4,8,8)=\frac{3}{\sqrt{5}}
!ang_phipphj(5,11,8)=\frac{4}{\sqrt{7}}
!ang_phipphj(6,2,8)=-\frac{1}{\sqrt{5}}
!ang_phipphj(6,12,8)=8 \sqrt{\frac{2}{35}}
!ang_phipphj(7,3,8)=-\frac{2}{\sqrt{15}}
!ang_phipphj(7,13,8)=\frac{12}{\sqrt{35}}
!ang_phipphj(8,4,8)=-\frac{1}{\sqrt{5}}
!ang_phipphj(8,14,8)=8 \sqrt{\frac{2}{35}}
!ang_phipphj(9,15,8)=\frac{4}{\sqrt{7}}
!ang_phipphj(11,5,8)=-\frac{2}{\sqrt{7}}
!ang_phipphj(12,6,8)=-4 \sqrt{\frac{2}{35}}
!ang_phipphj(13,7,8)=-\frac{6}{\sqrt{35}}
!ang_phipphj(14,8,8)=-4 \sqrt{\frac{2}{35}}
!ang_phipphj(15,9,8)=-\frac{2}{\sqrt{7}}


 ang_phipphj_tmp=zero
!
 ang_phipphj_tmp(1,4,1)=0.57735026918962576451_dp
 ang_phipphj_tmp(2,5,1)=0.44721359549995793928_dp
 ang_phipphj_tmp(3,8,1)=0.44721359549995793928_dp
 ang_phipphj_tmp(4,1,1)=0.57735026918962576451_dp
 ang_phipphj_tmp(4,7,1)=-0.25819888974716112568_dp
 ang_phipphj_tmp(4,9,1)=0.44721359549995793928_dp
 ang_phipphj_tmp(5,2,1)=0.44721359549995793928_dp
 ang_phipphj_tmp(5,10,1)=0.46291004988627573078_dp
 ang_phipphj_tmp(5,12,1)=-0.11952286093343936400_dp
 ang_phipphj_tmp(6,11,1)=0.37796447300922722721_dp
 ang_phipphj_tmp(7,4,1)=-0.25819888974716112568_dp
 ang_phipphj_tmp(7,14,1)=0.41403933560541253068_dp
 ang_phipphj_tmp(8,3,1)=0.44721359549995793928_dp
 ang_phipphj_tmp(8,13,1)=-0.29277002188455995381_dp
 ang_phipphj_tmp(8,15,1)=0.37796447300922722721_dp
 ang_phipphj_tmp(9,4,1)=0.44721359549995793928_dp
 ang_phipphj_tmp(9,14,1)=-0.11952286093343936400_dp
 ang_phipphj_tmp(9,16,1)=0.46291004988627573078_dp
 ang_phipphj_tmp(10,5,1)=0.46291004988627573078_dp
 ang_phipphj_tmp(11,6,1)=0.37796447300922722721_dp
 ang_phipphj_tmp(12,5,1)=-0.11952286093343936400_dp
 ang_phipphj_tmp(13,8,1)=-0.29277002188455995381_dp
 ang_phipphj_tmp(14,7,1)=0.41403933560541253068_dp
 ang_phipphj_tmp(14,9,1)=-0.11952286093343936400_dp
 ang_phipphj_tmp(15,8,1)=0.37796447300922722721_dp
 ang_phipphj_tmp(16,9,1)=0.46291004988627573078_dp
!
 ang_phipphj_tmp(1,4,2)=0.28867513459481288225_dp
 ang_phipphj_tmp(1,14,2)=-0.54006172486732168591_dp
 ang_phipphj_tmp(2,5,2)=0.22360679774997896964_dp
 ang_phipphj_tmp(3,8,2)=0.22360679774997896964_dp
 ang_phipphj_tmp(4,7,2)=-0.77459666924148337704_dp
 ang_phipphj_tmp(4,9,2)=0.22360679774997896964_dp
 ang_phipphj_tmp(5,2,2)=0.11180339887498948482_dp
 ang_phipphj_tmp(5,10,2)=0.23145502494313786539_dp
 ang_phipphj_tmp(5,12,2)=-0.47809144373375745599_dp
 ang_phipphj_tmp(6,11,2)=0.18898223650461361361_dp
 ang_phipphj_tmp(7,4,2)=0.25819888974716112568_dp
 ang_phipphj_tmp(7,14,2)=0.44854261357253024157_dp
 ang_phipphj_tmp(8,3,2)=-0.44721359549995793928_dp
 ang_phipphj_tmp(8,13,2)=-1.1710800875382398152_dp
 ang_phipphj_tmp(8,15,2)=0.18898223650461361361_dp
 ang_phipphj_tmp(9,4,2)=0.11180339887498948482_dp
 ang_phipphj_tmp(9,14,2)=-0.47809144373375745599_dp
 ang_phipphj_tmp(9,16,2)=0.23145502494313786539_dp
 ang_phipphj_tmp(10,5,2)=0.15430334996209191026_dp
 ang_phipphj_tmp(11,6,2)=-0.094491118252306806804_dp
 ang_phipphj_tmp(12,5,2)=0.23904572186687872799_dp
 ang_phipphj_tmp(13,8,2)=0.58554004376911990761_dp
 ang_phipphj_tmp(14,7,2)=-0.82807867121082506136_dp
 ang_phipphj_tmp(14,9,2)=0.23904572186687872799_dp
 ang_phipphj_tmp(15,8,2)=-0.094491118252306806804_dp
 ang_phipphj_tmp(16,9,2)=0.15430334996209191026_dp
!
 ang_phipphj_tmp(1,4,3)=0.86602540378443864676_dp
 ang_phipphj_tmp(1,14,3)=0.54006172486732168591_dp
 ang_phipphj_tmp(2,5,3)=1.1180339887498948482_dp
 ang_phipphj_tmp(3,8,3)=1.1180339887498948482_dp
 ang_phipphj_tmp(4,9,3)=1.1180339887498948482_dp
 ang_phipphj_tmp(5,2,3)=-0.55901699437494742410_dp
 ang_phipphj_tmp(5,10,3)=1.6201851746019650577_dp
 ang_phipphj_tmp(6,11,3)=1.3228756555322952953_dp
 ang_phipphj_tmp(7,14,3)=1.2076147288491198811_dp
 ang_phipphj_tmp(8,15,3)=1.3228756555322952953_dp
 ang_phipphj_tmp(9,4,3)=-0.55901699437494742410_dp
 ang_phipphj_tmp(9,16,3)=1.6201851746019650577_dp
 ang_phipphj_tmp(10,5,3)=-1.0801234497346433718_dp
 ang_phipphj_tmp(11,6,3)=-0.66143782776614764763_dp
 ang_phipphj_tmp(15,8,3)=-0.66143782776614764763_dp
 ang_phipphj_tmp(16,9,3)=-1.0801234497346433718_dp
!
 ang_phipphj_tmp(1,2,4)=0.57735026918962576451_dp
 ang_phipphj_tmp(2,1,4)=0.57735026918962576451_dp
 ang_phipphj_tmp(2,7,4)=-0.25819888974716112568_dp
 ang_phipphj_tmp(2,9,4)=-0.44721359549995793928_dp
 ang_phipphj_tmp(3,6,4)=0.44721359549995793928_dp
 ang_phipphj_tmp(4,5,4)=0.44721359549995793928_dp
 ang_phipphj_tmp(5,4,4)=0.44721359549995793928_dp
 ang_phipphj_tmp(5,14,4)=-0.11952286093343936400_dp
 ang_phipphj_tmp(5,16,4)=-0.46291004988627573078_dp
 ang_phipphj_tmp(6,3,4)=0.44721359549995793928_dp
 ang_phipphj_tmp(6,13,4)=-0.29277002188455995381_dp
 ang_phipphj_tmp(6,15,4)=-0.37796447300922722721_dp
 ang_phipphj_tmp(7,2,4)=-0.25819888974716112568_dp
 ang_phipphj_tmp(7,12,4)=0.41403933560541253068_dp
 ang_phipphj_tmp(8,11,4)=0.37796447300922722721_dp
 ang_phipphj_tmp(9,2,4)=-0.44721359549995793928_dp
 ang_phipphj_tmp(9,10,4)=0.46291004988627573078_dp
 ang_phipphj_tmp(9,12,4)=0.11952286093343936400_dp
 ang_phipphj_tmp(10,9,4)=0.46291004988627573078_dp
 ang_phipphj_tmp(11,8,4)=0.37796447300922722721_dp
 ang_phipphj_tmp(12,7,4)=0.41403933560541253068_dp
 ang_phipphj_tmp(12,9,4)=0.11952286093343936400_dp
 ang_phipphj_tmp(13,6,4)=-0.29277002188455995381_dp
 ang_phipphj_tmp(14,5,4)=-0.11952286093343936400_dp
 ang_phipphj_tmp(15,6,4)=-0.37796447300922722721_dp
 ang_phipphj_tmp(16,5,4)=-0.46291004988627573078_dp
!
 ang_phipphj_tmp(1,2,5)=0.28867513459481288225_dp
 ang_phipphj_tmp(1,12,5)=-0.54006172486732168591_dp
 ang_phipphj_tmp(2,7,5)=-0.77459666924148337704_dp
 ang_phipphj_tmp(2,9,5)=-0.22360679774997896964_dp
 ang_phipphj_tmp(3,6,5)=0.22360679774997896964_dp
 ang_phipphj_tmp(4,5,5)=0.22360679774997896964_dp
 ang_phipphj_tmp(5,4,5)=0.11180339887498948482_dp
 ang_phipphj_tmp(5,14,5)=-0.47809144373375745599_dp
 ang_phipphj_tmp(5,16,5)=-0.23145502494313786539_dp
 ang_phipphj_tmp(6,3,5)=-0.44721359549995793928_dp
 ang_phipphj_tmp(6,13,5)=-1.1710800875382398152_dp
 ang_phipphj_tmp(6,15,5)=-0.18898223650461361361_dp
 ang_phipphj_tmp(7,2,5)=0.25819888974716112568_dp
 ang_phipphj_tmp(7,12,5)=0.44854261357253024157_dp
 ang_phipphj_tmp(8,11,5)=0.18898223650461361361_dp
 ang_phipphj_tmp(9,2,5)=-0.11180339887498948482_dp
 ang_phipphj_tmp(9,10,5)=0.23145502494313786539_dp
 ang_phipphj_tmp(9,12,5)=0.47809144373375745599_dp
 ang_phipphj_tmp(10,9,5)=0.15430334996209191026_dp
 ang_phipphj_tmp(11,8,5)=-0.094491118252306806804_dp
 ang_phipphj_tmp(12,7,5)=-0.82807867121082506136_dp
 ang_phipphj_tmp(12,9,5)=-0.23904572186687872799_dp
 ang_phipphj_tmp(13,6,5)=0.58554004376911990761_dp
 ang_phipphj_tmp(14,5,5)=0.23904572186687872799_dp
 ang_phipphj_tmp(15,6,5)=0.094491118252306806804_dp
 ang_phipphj_tmp(16,5,5)=-0.15430334996209191026_dp
!
 ang_phipphj_tmp(1,2,6)=0.86602540378443864676_dp
 ang_phipphj_tmp(1,12,6)=0.54006172486732168591_dp
 ang_phipphj_tmp(2,9,6)=-1.1180339887498948482_dp
 ang_phipphj_tmp(3,6,6)=1.1180339887498948482_dp
 ang_phipphj_tmp(4,5,6)=1.1180339887498948482_dp
 ang_phipphj_tmp(5,4,6)=-0.55901699437494742410_dp
 ang_phipphj_tmp(5,16,6)=-1.6201851746019650577_dp
 ang_phipphj_tmp(6,15,6)=-1.3228756555322952953_dp
 ang_phipphj_tmp(7,12,6)=1.2076147288491198811_dp
 ang_phipphj_tmp(8,11,6)=1.3228756555322952953_dp
 ang_phipphj_tmp(9,2,6)=0.55901699437494742410_dp
 ang_phipphj_tmp(9,10,6)=1.6201851746019650577_dp
 ang_phipphj_tmp(10,9,6)=-1.0801234497346433718_dp
 ang_phipphj_tmp(11,8,6)=-0.66143782776614764763_dp
 ang_phipphj_tmp(15,6,6)=0.66143782776614764763_dp
 ang_phipphj_tmp(16,5,6)=1.0801234497346433718_dp
!
 ang_phipphj_tmp(1,3,7)=0.57735026918962576451_dp
 ang_phipphj_tmp(2,6,7)=0.44721359549995793928_dp
 ang_phipphj_tmp(3,1,7)=0.57735026918962576451_dp
 ang_phipphj_tmp(3,7,7)=0.51639777949432225136_dp
 ang_phipphj_tmp(4,8,7)=0.44721359549995793928_dp
 ang_phipphj_tmp(5,11,7)=0.37796447300922722721_dp
 ang_phipphj_tmp(6,2,7)=0.44721359549995793928_dp
 ang_phipphj_tmp(6,12,7)=0.47809144373375745599_dp
 ang_phipphj_tmp(7,3,7)=0.51639777949432225136_dp
 ang_phipphj_tmp(7,13,7)=0.50709255283710994651_dp
 ang_phipphj_tmp(8,4,7)=0.44721359549995793928_dp
 ang_phipphj_tmp(8,14,7)=0.47809144373375745599_dp
 ang_phipphj_tmp(9,15,7)=0.37796447300922722721_dp
 ang_phipphj_tmp(11,5,7)=0.37796447300922722721_dp
 ang_phipphj_tmp(12,6,7)=0.47809144373375745599_dp
 ang_phipphj_tmp(13,7,7)=0.50709255283710994651_dp
 ang_phipphj_tmp(14,8,7)=0.47809144373375745599_dp
 ang_phipphj_tmp(15,9,7)=0.37796447300922722721_dp
!
 ang_phipphj_tmp(1,3,8)=1.1547005383792515290_dp
 ang_phipphj_tmp(2,6,8)=1.3416407864998738178_dp
 ang_phipphj_tmp(3,7,8)=1.5491933384829667541_dp
 ang_phipphj_tmp(4,8,8)=1.3416407864998738178_dp
 ang_phipphj_tmp(5,11,8)=1.5118578920369089089_dp
 ang_phipphj_tmp(6,2,8)=-0.44721359549995793928_dp
 ang_phipphj_tmp(6,12,8)=1.9123657749350298240_dp
 ang_phipphj_tmp(7,3,8)=-0.51639777949432225136_dp
 ang_phipphj_tmp(7,13,8)=2.0283702113484397860_dp
 ang_phipphj_tmp(8,4,8)=-0.44721359549995793928_dp
 ang_phipphj_tmp(8,14,8)=1.9123657749350298240_dp
 ang_phipphj_tmp(9,15,8)=1.5118578920369089089_dp
 ang_phipphj_tmp(11,5,8)=-0.75592894601845445443_dp
 ang_phipphj_tmp(12,6,8)=-0.95618288746751491198_dp
 ang_phipphj_tmp(13,7,8)=-1.0141851056742198930_dp
 ang_phipphj_tmp(14,8,8)=-0.95618288746751491198_dp
 ang_phipphj_tmp(15,9,8)=-0.75592894601845445443_dp

 ang_phipphj(:,:,:)=ang_phipphj_tmp(1:mpsang**2,1:mpsang**2,:)


 end subroutine int_ang
!!***
