!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_b_matrix
!! NAME
!! calc_b_matrix
!!
!! FUNCTION
!!  calculate values of derivatives of internal coordinates as a function of
!!  cartesian ones =  B matrix
!!
!! COPYRIGHT
!! Copyright (C) 2003-2017 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! angs= number of angles
!! bonds(2,2,nbond)=for a bond between iatom and jatom
!!              bonds(1,1,nbond) = iatom
!!              bonds(2,1,nbond) = icenter
!!              bonds(1,2,nbond) = jatom
!!              bonds(2,2,nbond) = irshift
!! carts(2,ncart)= index of total primitive internal, and atom (carts(2,:))
!! dihedrals(2,4,ndihed)=indexes to characterize dihedrals
!! nang(2,3,nang)=indexes to characterize angles
!! nbond=number of bonds
!! ncart=number of auxiliary cartesian atom coordinates (used for constraints)
!! ndihed= number of dihedrals
!! ninternal=nbond+nang+ndihed+ncart: number of internal coordinates
!! nrshift= dimension of rshift
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! rshift(3,nrshift)=shift in xred that must be done to find all neighbors of
!!                   a given atom within a given number of neighboring shells
!! xcart(3,natom)=cartesian coordinates of atoms (bohr)
!!
!! OUTPUT
!! b_matrix(ninternal,3*natom)=matrix of derivatives of internal coordinates
!!   wrt cartesians
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      xcart2deloc
!!
!! CHILDREN
!!      acrossb
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine calc_b_matrix(deloc,natom,rprimd,xcart,b_matrix)

 use defs_basis
 use m_abimover
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_b_matrix'
 use interfaces_45_geomoptim, except_this_one => calc_b_matrix
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
 type(delocint),intent(inout) :: deloc

!arrays
 real(dp),intent(in) :: rprimd(3,3),xcart(3,natom)
 real(dp),intent(out) :: b_matrix(deloc%ninternal,3*natom)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,i4,iang,ibond,icart,idihed,iprim,s1,s2,s3,s4
!arrays
 real(dp) :: bb(3),r1(3),r2(3),r3(3),r4(3)

! *************************************************************************

 iprim=0
 b_matrix(:,:) = zero

 do ibond=1,deloc%nbond
   i1 = deloc%bonds(1,1,ibond)
   s1 = deloc%bonds(2,1,ibond)
   r1(:) = xcart(:,i1)+deloc%rshift(1,s1)*rprimd(:,1)&
&   +deloc%rshift(2,s1)*rprimd(:,2)&
&   +deloc%rshift(3,s1)*rprimd(:,3)
   i2 = deloc%bonds(1,2,ibond)
   s2 = deloc%bonds(2,2,ibond)
   r2(:) = xcart(:,i2)+deloc%rshift(1,s2)*rprimd(:,1)&
&   +deloc%rshift(2,s2)*rprimd(:,2)&
&   +deloc%rshift(3,s2)*rprimd(:,3)
   iprim=iprim+1
   call dbond_length_d1(r1,r2,bb)
   b_matrix(iprim,3*(i1-1)+1:3*i1) = b_matrix(iprim,3*(i1-1)+1:3*i1) + bb(:)
   call dbond_length_d1(r2,r1,bb)
   b_matrix(iprim,3*(i2-1)+1:3*i2) = b_matrix(iprim,3*(i2-1)+1:3*i2) + bb(:)
 end do

!second: angle values (ang)
 do iang=1,deloc%nang
   i1 = deloc%angs(1,1,iang)
   s1 = deloc%angs(2,1,iang)
   r1(:) = xcart(:,i1)+deloc%rshift(1,s1)*rprimd(:,1)&
&   +deloc%rshift(2,s1)*rprimd(:,2)&
&   +deloc%rshift(3,s1)*rprimd(:,3)
   i2 = deloc%angs(1,2,iang)
   s2 = deloc%angs(2,2,iang)
   r2(:) = xcart(:,i2)+deloc%rshift(1,s2)*rprimd(:,1)&
&   +deloc%rshift(2,s2)*rprimd(:,2)&
&   +deloc%rshift(3,s2)*rprimd(:,3)
   i3 = deloc%angs(1,3,iang)
   s3 = deloc%angs(2,3,iang)
   r3(:) = xcart(:,i3)+deloc%rshift(1,s3)*rprimd(:,1)&
&   +deloc%rshift(2,s3)*rprimd(:,2)&
&   +deloc%rshift(3,s3)*rprimd(:,3)
   iprim=iprim+1
   call dang_d1(r1,r2,r3,bb)
   b_matrix(iprim,3*(i1-1)+1:3*i1) = b_matrix(iprim,3*(i1-1)+1:3*i1) + bb(:)
   call dang_d2(r1,r2,r3,bb)
   b_matrix(iprim,3*(i2-1)+1:3*i2) = b_matrix(iprim,3*(i2-1)+1:3*i2) + bb(:)
   call dang_d1(r3,r2,r1,bb)
   b_matrix(iprim,3*(i3-1)+1:3*i3) = b_matrix(iprim,3*(i3-1)+1:3*i3) + bb(:)
 end do

!third: dihedral values
 do idihed=1,deloc%ndihed
   i1 = deloc%dihedrals(1,1,idihed)
   s1 = deloc%dihedrals(2,1,idihed)
   r1(:) = xcart(:,i1)+deloc%rshift(1,s1)*rprimd(:,1)&
&   +deloc%rshift(2,s1)*rprimd(:,2)&
&   +deloc%rshift(3,s1)*rprimd(:,3)
   i2 = deloc%dihedrals(1,2,idihed)
   s2 = deloc%dihedrals(2,2,idihed)
   r2(:) = xcart(:,i2)+deloc%rshift(1,s2)*rprimd(:,1)&
&   +deloc%rshift(2,s2)*rprimd(:,2)&
&   +deloc%rshift(3,s2)*rprimd(:,3)
   i3 = deloc%dihedrals(1,3,idihed)
   s3 = deloc%dihedrals(2,3,idihed)
   r3(:) = xcart(:,i3)+deloc%rshift(1,s3)*rprimd(:,1)&
&   +deloc%rshift(2,s3)*rprimd(:,2)&
&   +deloc%rshift(3,s3)*rprimd(:,3)
   i4 = deloc%dihedrals(1,4,idihed)
   s4 = deloc%dihedrals(2,4,idihed)
   r4(:) = xcart(:,i4)+deloc%rshift(1,s4)*rprimd(:,1)&
&   +deloc%rshift(2,s4)*rprimd(:,2)&
&   +deloc%rshift(3,s4)*rprimd(:,3)
!  write(std_out,*) 'dihed ',idihed
!  write(std_out,*) r1
!  write(std_out,*) r2
!  write(std_out,*) r3
!  write(std_out,*) r4

   iprim=iprim+1
   call ddihedral_d1(r1,r2,r3,r4,bb)
   b_matrix(iprim,3*(i1-1)+1:3*i1) = b_matrix(iprim,3*(i1-1)+1:3*i1) + bb(:)
   call ddihedral_d2(r1,r2,r3,r4,bb)
   b_matrix(iprim,3*(i2-1)+1:3*i2) = b_matrix(iprim,3*(i2-1)+1:3*i2) + bb(:)
   call ddihedral_d2(r4,r3,r2,r1,bb)
   b_matrix(iprim,3*(i3-1)+1:3*i3) = b_matrix(iprim,3*(i3-1)+1:3*i3) + bb(:)
   call ddihedral_d1(r4,r3,r2,r1,bb)
   b_matrix(iprim,3*(i4-1)+1:3*i4) = b_matrix(iprim,3*(i4-1)+1:3*i4) + bb(:)
 end do

 do icart=1,deloc%ncart
   iprim=iprim+1
   b_matrix(iprim,3*(deloc%carts(2,icart)-1)+deloc%carts(1,icart)) = &
&   b_matrix(iprim,3*(deloc%carts(2,icart)-1)+deloc%carts(1,icart)) + one
 end do

!DEBUG
! write (200,*) 'calc_b_matrix : b_matrix = '
! do iprim=1,deloc%ninternal
!   do i1=1, 3*natom
!     write (200,'(E16.6,2x)',ADVANCE='NO') b_matrix(iprim,i1)
!   end do
!   write (200,*)
! end do
!ENDDEBUG


end subroutine calc_b_matrix
!!***


!!****f* ABINIT/dbond_length_d1
!! NAME
!! dbond_length_d1
!!
!! FUNCTION
!!
!! PARENTS
!!      calc_b_matrix
!!
!! CHILDREN
!!      acrossb
!!
!! SOURCE
!!

subroutine dbond_length_d1(r1,r2,bb)

 use defs_basis
 use m_profiling_abi
 use m_abimover

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dbond_length_d1'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: r1(3),r2(3)
 real(dp),intent(out) :: bb(3)

!Local variables ------------------------------
!arrays
 real(dp) :: rpt(3)

!************************************************************************
 rpt(:) = r1(:)-r2(:)
 bb(:) = rpt(:)/bond_length(r1,r2)

end subroutine dbond_length_d1
!!***


!!****f* ABINIT/dang_d1
!! NAME
!! dang_d1
!!
!! FUNCTION
!!
!! PARENTS
!!      calc_b_matrix
!!
!! CHILDREN
!!      acrossb
!!
!! SOURCE
!!

subroutine dang_d1(r1,r2,r3,bb)

 use defs_basis
 use m_profiling_abi
 use m_abimover

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dang_d1'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: r1(3),r2(3),r3(3)
 real(dp),intent(out) :: bb(3)

!Local variables ------------------------------
!scalars
 real(dp) :: cos_ang,n1,n1232,n2,tmp
!arrays
 real(dp) :: cp1232(3),rpt(3),rpt12(3),rpt32(3)

!************************************************************************
 n1=bond_length(r1,r2)
 n2=bond_length(r3,r2)

 rpt12(:) = r1(:)-r2(:)
 rpt32(:) = r3(:)-r2(:)

 cos_ang = (rpt12(1)*rpt32(1)+rpt12(2)*rpt32(2)+rpt12(3)*rpt32(3))/n1/n2
 if (cos_ang > one - epsilon(one)*two) then
   cos_ang = one
 else if(cos_ang < -one + epsilon(one)*two) then
   cos_ang = -one
 end if

 rpt(:) = rpt32(:)/n1/n2 - rpt12(:)*cos_ang/n1/n1

 tmp = sqrt(one-cos_ang**2)
 bb(:) = zero
 if (tmp > epsilon(one)) then
   bb(:) = rpt(:) * (-one)/tmp
 end if

!TEST: version from MOLECULAR VIBRATIONS EB Wilson
 call acrossb(rpt12,rpt32,cp1232)
 n1232 = sqrt(cp1232(1)**2+cp1232(2)**2+cp1232(3)**2)
 rpt(:) = (cos_ang*rpt12(:)*n2/n1 - rpt32(:))/n1232
 if (abs(bb(1)-rpt(1))+abs(bb(2)-rpt(2))+abs(bb(3)-rpt(3)) > tol10) then
   write(std_out,*) 'Compare bb ang 1 : '
   write(std_out,*) bb(:), rpt(:), bb(:)-rpt(:)
 end if
 bb(:) = rpt(:)

end subroutine dang_d1
!!***


!!****f* ABINIT/dang_d2
!! NAME
!! dang_d2
!!
!! FUNCTION
!!
!! PARENTS
!!      calc_b_matrix
!!
!! CHILDREN
!!      acrossb
!!
!! SOURCE
!!

subroutine dang_d2(r1,r2,r3,bb)

 use defs_basis
 use m_profiling_abi
 use m_abimover

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dang_d2'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: r1(3),r2(3),r3(3)
 real(dp),intent(out) :: bb(3)

!Local variables ------------------------------
!scalars
 real(dp) :: cos_ang,n1,n1232,n2,tmp
!arrays
 real(dp) :: cp1232(3),rpt(3),rpt12(3),rpt32(3)

!************************************************************************
 n1=bond_length(r1,r2)
 n2=bond_length(r3,r2)

 rpt12(:) = r1(:)-r2(:)
 rpt32(:) = r3(:)-r2(:)

 cos_ang = (rpt12(1)*rpt32(1)+rpt12(2)*rpt32(2)+rpt12(3)*rpt32(3))/n1/n2
 if (cos_ang > one - epsilon(one)*two) then
   cos_ang = one
 else if(cos_ang < -one + epsilon(one)*two) then
   cos_ang = -one
 end if

 rpt(:) = -rpt32(:)/n1/n2 - rpt12(:)/n1/n2 &
& + rpt12(:)*cos_ang/n1/n1 + rpt32(:)*cos_ang/n2/n2

 tmp = sqrt(one-cos_ang**2)
 bb(:) = zero
 if (tmp > tol12) then
   bb(:) = rpt(:) * (-one)/tmp
 end if

!TEST: version from MOLECULAR VIBRATIONS EB Wilson
 call acrossb(rpt12,rpt32,cp1232)
 n1232 = sqrt(cp1232(1)**2+cp1232(2)**2+cp1232(3)**2)
 rpt(:) = ((n1-n2*cos_ang)*rpt12(:)/n1 + (n2-n1*cos_ang)*rpt32(:)/n2) / n1232
 if (abs(bb(1)-rpt(1))+abs(bb(2)-rpt(2))+abs(bb(3)-rpt(3))  > tol10) then
   write(std_out,*) 'Compare bb ang 2 : '
   write(std_out,*) bb(:), rpt(:), bb(:)-rpt(:)
 end if
 bb(:) = rpt(:)


end subroutine dang_d2
!!***


!!****f* ABINIT/ddihedral_d1
!! NAME
!! ddihedral_d1
!!
!! FUNCTION
!!
!! PARENTS
!!      calc_b_matrix
!!
!! CHILDREN
!!      acrossb
!!
!! SOURCE
!!

subroutine ddihedral_d1(r1,r2,r3,r4,bb)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddihedral_d1'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: r1(3),r2(3),r3(3),r4(3)
 real(dp),intent(out) :: bb(3)

!Local variables ------------------------------------
!scalars
 real(dp) :: cos_dihedral,dih_sign,n1,n2,n23,sin_dihedral,tmp
!arrays
 real(dp) :: cp1232(3),cp32_1232(3),cp32_3432(3),cp3432(3),cpcp(3),rpt(3)
 real(dp) :: rpt12(3),rpt32(3),rpt34(3)

!******************************************************************
 rpt12(:) = r1(:)-r2(:)
 rpt32(:) = r3(:)-r2(:)
 rpt34(:) = r3(:)-r4(:)

 call acrossb(rpt12,rpt32,cp1232)
 call acrossb(rpt34,rpt32,cp3432)

!DEBUG
!write(std_out,*) ' cos_dihedral : cp1232 = ', cp1232
!write(std_out,*) ' cos_dihedral : cp3432 = ', cp3432
!ENDDEBUG

 n1 = sqrt(cp1232(1)**2+cp1232(2)**2+cp1232(3)**2)
 n2 = sqrt(cp3432(1)**2+cp3432(2)**2+cp3432(3)**2)

 cos_dihedral = (cp1232(1)*cp3432(1)+cp1232(2)*cp3432(2)+cp1232(3)*cp3432(3))/n1/n2
 if (cos_dihedral > one - epsilon(one)*two) then
   cos_dihedral = one
 else if(cos_dihedral < -one + epsilon(one)*two) then
   cos_dihedral = -one
 end if
!we use complementary of standard angle, so
!cos_dihedral = -cos_dihedral

 call acrossb(cp1232,cp3432,cpcp)
 cpcp(:) = cpcp(:)/n1/n2
!we use complementary of standard angle, but sin is invariant
 sin_dihedral = -(cpcp(1)*rpt32(1)+cpcp(2)*rpt32(2)+cpcp(3)*rpt32(3))&
& /sqrt(rpt32(1)**2+rpt32(2)**2+rpt32(3)**2)
 dih_sign = one
 if (sin_dihedral < -epsilon(one)) then
   dih_sign = -one
 end if

!DEBUG
!write(std_out,'(a,3E16.6)') 'ddihedral_d1 : cos abs(sin) dih_sign= ',&
!&    cos_dihedral,sin_dihedral,dih_sign
!ENDDEBUG

!ddihedral_d1 = dih_sign* acos(cos_dihedral)
 call acrossb(rpt32,cp1232,cp32_1232)
 call acrossb(rpt32,cp3432,cp32_3432)

 rpt(:) = cp32_3432(:)/n1/n2 - cp32_1232(:)/n1/n1 * cos_dihedral
 bb(:) = zero

!DEBUG
!write(std_out,*) 'ddihedral_d1 cp1232 cp3432 = ',cp1232,cp3432,rpt32
!write(std_out,*) 'ddihedral_d1 cp32_1232 cp32_3432 = ',cp32_1232,cp32_3432,cos_dihedral,n1,n2
!write(std_out,*) 'ddihedral_d1 rpt = ',rpt
!ENDDEBUG

 tmp = sqrt(one-cos_dihedral**2)
 if (tmp > tol12) then
!  we use complementary of standard angle, so cosine in acos has - sign,
!  and it appears for the derivative
   bb(:) = -dih_sign * rpt(:) * (-one) / tmp
 else
   bb(:) = dih_sign * cp32_3432(:) / n1 / n2 / &
&   sqrt(cp32_3432(1)**2+cp32_3432(2)**2+cp32_3432(3)**2)
 end if

!TEST: version from MOLECULAR VIBRATIONS EB Wilson

 n23 = sqrt(rpt32(1)*rpt32(1)+rpt32(2)*rpt32(2)+rpt32(3)*rpt32(3))
 rpt(:) = cp1232(:)*n23/n1/n1
!if (abs(bb(1)-rpt(1))+abs(bb(2)-rpt(2))+abs(bb(3)-rpt(3))  > tol10) then
!write(std_out,*) 'Compare bb1 : '
!write(std_out,*) bb(:), rpt(:), bb(:)-rpt(:)
!end if
 bb(:) = rpt(:)

end subroutine ddihedral_d1
!!***


!!****f* ABINIT/ddihedral_d2
!! NAME
!! ddihedral_d2
!!
!! FUNCTION
!!
!! PARENTS
!!      calc_b_matrix
!!
!! CHILDREN
!!      acrossb
!!
!! SOURCE
!!

subroutine ddihedral_d2(r1,r2,r3,r4,bb)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddihedral_d2'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: r1(3),r2(3),r3(3),r4(3)
 real(dp),intent(out) :: bb(3)

!Local variables
!scalars
 real(dp) :: cos_dihedral,dih_sign,n1,n2,n23,sin_dihedral,sp1232,sp3432,tmp
!arrays
 real(dp) :: cp1232(3),cp1232_12(3),cp1232_34(3),cp32_1232(3),cp32_3432(3)
 real(dp) :: cp3432(3),cp3432_12(3),cp3432_34(3),cpcp(3),rpt(3),rpt12(3)
 real(dp) :: rpt32(3),rpt34(3)

! *************************************************************************
 rpt12(:) = r1(:)-r2(:)
 rpt32(:) = r3(:)-r2(:)
 rpt34(:) = r3(:)-r4(:)

 call acrossb(rpt12,rpt32,cp1232)
 call acrossb(rpt34,rpt32,cp3432)

!DEBUG
!write(std_out,*) ' cos_dihedral : cp1232 = ', cp1232
!write(std_out,*) ' cos_dihedral : cp3432 = ', cp3432
!ENDDEBUG

 n1 = sqrt(cp1232(1)**2+cp1232(2)**2+cp1232(3)**2)
 n2 = sqrt(cp3432(1)**2+cp3432(2)**2+cp3432(3)**2)

 cos_dihedral = (cp1232(1)*cp3432(1)+cp1232(2)*cp3432(2)+cp1232(3)*cp3432(3))/n1/n2
 if (cos_dihedral > one - epsilon(one)*two) then
   cos_dihedral = one
 else if(cos_dihedral < -one + epsilon(one)*two) then
   cos_dihedral = -one
 end if
!we use complementary of standard angle, so
!cos_dihedral = -cos_dihedral

 call acrossb(cp1232,cp3432,cpcp)
 cpcp(:) = cpcp(:)/n1/n2
!we use complementary of standard angle, but sin is invariant
 sin_dihedral = -(cpcp(1)*rpt32(1)+cpcp(2)*rpt32(2)+cpcp(3)*rpt32(3))&
& /sqrt(rpt32(1)**2+rpt32(2)**2+rpt32(3)**2)
 dih_sign = one
 if (sin_dihedral <  -tol12) then
   dih_sign = -one
 end if

!DEBUG
!write(std_out,'(a,3E16.6)') 'ddihedral_d2 : cos abs(sin) dih_sign= ',&
!&    cos_dihedral,sin_dihedral,dih_sign
!ENDDEBUG

!ddihedral_d2 = dih_sign* acos(cos_dihedral)
 call acrossb(rpt32,cp3432,cp32_3432)
 call acrossb(cp3432,rpt12,cp3432_12)
 call acrossb(cp1232,rpt34,cp1232_34)

 call acrossb(rpt32,cp1232,cp32_1232)
 call acrossb(cp1232,rpt12,cp1232_12)
 call acrossb(cp3432,rpt34,cp3432_34)

 rpt(:) = -(cp32_3432(:) + cp3432_12(:) + cp1232_34(:))/n1/n2 &
& +cos_dihedral*(cp32_1232(:)/n1/n1 + cp1232_12(:)/n1/n1 + cp3432_34(:)/n2/n2)
 bb(:) = zero
 tmp = sqrt(one-cos_dihedral**2)
 if (tmp > tol12) then
!  we use complementary of standard angle, so cosine in acos has - sign,
!  and it appears for derivative
   bb(:) = -dih_sign * rpt(:) * (-one) / tmp
 else
   bb(:) = dih_sign * cos_dihedral * &
&   ( cp32_1232(:)/n1/n1/sqrt(cp32_1232(1)**2+cp32_1232(2)**2+cp32_1232(3)**2) &
&   +cp1232_12(:)/n1/n1/sqrt(cp1232_12(1)**2+cp1232_12(2)**2+cp1232_12(3)**2) &
&   +cp3432_34(:)/n2/n2/sqrt(cp3432_34(1)**2+cp3432_34(2)**2+cp3432_34(3)**2) )
 end if

!TEST: version from MOLECULAR VIBRATIONS EB Wilson p. 61
 n23 = sqrt(rpt32(1)*rpt32(1)+rpt32(2)*rpt32(2)+rpt32(3)*rpt32(3))
 sp1232 = rpt12(1)*rpt32(1)+rpt12(2)*rpt32(2)+rpt12(3)*rpt32(3)
 sp3432 = rpt34(1)*rpt32(1)+rpt34(2)*rpt32(2)+rpt34(3)*rpt32(3)

 rpt(:) = -cp1232(:)*(n23-sp1232/n23)/n1/n1 - cp3432(:)*sp3432/n23/n2/n2
!if (abs(bb(1)-rpt(1))+abs(bb(2)-rpt(2))+abs(bb(3)-rpt(3))  > tol10) then
!write(std_out,*) 'Compare bb2 : '
!write(std_out,*) bb(:), rpt(:), bb(:)-rpt(:)
!write(std_out,*) -cp1232(:)*(n23-sp1232/n23)/n1/n1, -cp3432(:)*sp3432/n23/n2/n2
!end if
 bb(:) = rpt(:)


end subroutine ddihedral_d2
!!***
