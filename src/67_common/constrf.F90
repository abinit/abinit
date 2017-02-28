!{\src2tex{textfont=tt}}
!!****f* ABINIT/constrf
!! NAME
!! constrf
!!
!! FUNCTION
!! Computes projected forces, fpcart, which satisfy a set of
!! constraint equations of the form
!!  Sum[mu,iatom]: wtatcon(mu,iatom,iconeq)*fpcart(mu,iatom) = 0 (iconeq=1,nconeq).
!! These projected forces are returned in fcart and thus replace
!! the original forces.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (SCE, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  iatfix(3,natom)=1 for frozen atom along each direction, 0 for unfrozen
!!  natom=number of atoms in cell
!!  nconeq=number of atomic constraint equations
!!  prtvol=control print volume and debugging
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  wtatcon(3,natom,nconeq)=weights for atomic constraints
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  diffor=maximum absolute change in component of projected fcart between present
!!          and previous SCF cycle
!!  fred(3,natom)=grads of Etot wrt reduced coordinates (hartree)
!!  maxfor=maximum absolute value of fcart
!!
!! SIDE EFFECTS
!!  fcart(3,natom)=cartesian forces (hartree/bohr) on input, projected forces on output
!!  forold(3,natom)=cartesian forces of previous SCF cycle (hartree/bohr)
!!
!! TODO
!!
!! PARENTS
!!      forces
!!
!! CHILDREN
!!      dposv,prtxvf,wrtout,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine constrf(diffor,fcart,forold,fred,iatfix,ionmov,maxfor,natom,&
& nconeq,prtvol,rprimd,wtatcon,xred)

 use defs_basis
 use m_linalg_interfaces
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'constrf'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_67_common, except_this_one => constrf
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ionmov,natom,nconeq,prtvol
 real(dp),intent(out) :: diffor,maxfor
!arrays
 integer,intent(in) :: iatfix(3,natom)
 real(dp),intent(in) :: rprimd(3,3),wtatcon(3,natom,nconeq)
 real(dp),intent(inout) :: fcart(3,natom),forold(3,natom),xred(3,natom)
 real(dp),intent(inout) :: fred(3,natom) !vz_i

!Local variables -------------------------
!scalars
 integer :: iatom,iconeq,iconeq1,iconeq2,index,info,mu,prtvel
 character(len=500) :: message
!arrays
 real(dp),allocatable :: fcartvec(:),fpcart(:,:),fvector(:),vel_dummy(:,:)
 real(dp),allocatable :: wmatrix(:,:),wtatconvec(:,:),wtcoeffs(:),xcart(:,:)

!************************************************************************

!Allocate temporary variables
 ABI_ALLOCATE(fpcart,(3,natom))
 ABI_ALLOCATE(fcartvec,(3*natom))
 ABI_ALLOCATE(fvector,(nconeq))
 ABI_ALLOCATE(vel_dummy,(3,natom))
 ABI_ALLOCATE(wmatrix,(nconeq,nconeq))
 ABI_ALLOCATE(wtatconvec,(3*natom,nconeq))
 ABI_ALLOCATE(wtcoeffs,(nconeq))
 ABI_ALLOCATE(xcart,(3,natom))

!If prtvol>10, output coordinates and forces prior to projecting
 if(prtvol>=10)then
   write(message,'(a)')' constrf - coordinates and forces prior to constraint projections:'
   call wrtout(std_out,message,'COLL')
   call xred2xcart(natom,rprimd,xcart,xred)
   prtvel=0
   call prtxvf(fcart,fred,iatfix,06,natom,prtvel,vel_dummy,xcart,xred)
 end if

!Transfer fcart and wtatcon to flat vectors
 index=0
 do iatom=1,natom
   do mu=1,3
     index=index+1
     fcartvec(index)=fcart(mu,iatom)
     wtatconvec(index,:)=wtatcon(mu,iatom,:)
   end do
 end do

!Compute a matrix (wmatrix) and vector (fvector) such that solving
!the linear equations wmatrix*wcoeffs=fvector gives the coefficients
!of wtatcon (wcoeffs) needed to compute the projected forces
 do iconeq2=1,nconeq
   fvector(iconeq2)=ddot(3*natom,fcartvec,1,wtatconvec(1,iconeq2),1)
   do iconeq1=1,nconeq
     wmatrix(iconeq1,iconeq2)=ddot(3*natom,wtatconvec(1,iconeq1),1,wtatconvec(1,iconeq2),1)
   end do
 end do

!Solve the system of linear equations, wmatrix*wcoeffs=fvector
 call dposv('U',nconeq,1,wmatrix,nconeq,fvector,nconeq,info)

 if (info/=0) then
   write(message, '(a,a,a,a,a)' )&
&   '  Constraint matrix is not positive definite,',ch10,&
&   '  probably because constraints are linearly dependent.',ch10,&
&   '  Action : Check for linear dependence of constraints.'
   MSG_ERROR(message)
 end if

!The solution vector is returned in fvector, so copy it to a more sensible location
 wtcoeffs(:)=fvector(:)

!Compute the projected forces, which now satisfy all the constraint equations
 fpcart(:,:)=fcart(:,:)
 do iconeq=1,nconeq
   fpcart(:,:)=fpcart(:,:)-wtcoeffs(iconeq)*wtatcon(:,:,iconeq)
 end do

!Reconvert constrained forces back from fpcart to fred
 do iatom=1,natom
   do mu=1,3
     fred(mu,iatom)= - (rprimd(1,mu)*fpcart(1,iatom)+&
&     rprimd(2,mu)*fpcart(2,iatom)+&
&     rprimd(3,mu)*fpcart(3,iatom))
   end do
 end do

!If prtvol>=10, output coordinates and forces after projecting
 if(prtvol>=10)then
   write(message,'(a)')' constrf - coordinates and forces after constraint projections:'
   call wrtout(std_out,message,'COLL')
   prtvel=0
   call prtxvf(fpcart,fred,iatfix,06,natom,prtvel,vel_dummy,xcart,xred)
 end if

!Copy the constrained forces, fpcart, back to fcart
 fcart(:,:)=fpcart(:,:)

!Compute maximal force and maximal difference of the projected forces,
!overriding the values already computed in forces
 maxfor=0.0_dp
 diffor=0.0_dp
 do iatom=1,natom
   do mu=1,3
     if (iatfix(mu,iatom) /= 1) then
       maxfor=max(maxfor,abs(fcart(mu,iatom)))
       diffor=max(diffor,abs(fcart(mu,iatom)-forold(mu,iatom)))
     else if (ionmov==4 .or. ionmov==5) then
!      Make the force vanish on fixed atoms when ionmov=4 or 5
!      This is because fixing of atom cannot be imposed at the
!      level of a routine similar to brdmin or moldyn for these options.
       fcart(mu,iatom)=0.0_dp
     end if
   end do
 end do
 forold(:,:)=fcart(:,:)

!Dellocate temporary variables
 ABI_DEALLOCATE(fpcart)
 ABI_DEALLOCATE(fcartvec)
 ABI_DEALLOCATE(fvector)
 ABI_DEALLOCATE(vel_dummy)
 ABI_DEALLOCATE(wmatrix)
 ABI_DEALLOCATE(wtatconvec)
 ABI_DEALLOCATE(wtcoeffs)
 ABI_DEALLOCATE(xcart)

end subroutine constrf
!!***
