!{\src2tex{textfont=tt}}
!!****f* ABINIT/mag_constr_e
!! NAME
!! mag_constr_e
!!
!! FUNCTION
!! This routine is called to compute the energy corresponding to constrained magnetic moments.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (ILuk)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  magconon=constraining option (on/off); 1=fix only the direction, 2=fix the direction and size
!!  spinat=fixed magnetic moments vectors
!!  magcon_lambda=the size of the penalty terms
!!
!! OUTPUT
!!  Epen=penalty contribution to the total energy corresponding to the constrained potential
!!  Econstr=???
!!  Eexp=???
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      calcdensph,metric,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mag_constr_e(magconon,magcon_lambda,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,ratsph,rhor,rprimd,spinat,typat,xred)

!use functions
use defs_basis
use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mag_constr_e'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_54_abiutil
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,magconon,nspden,nfft,ntypat
 real(dp),intent(in) :: magcon_lambda
!arrays
 integer, intent(in) :: ngfft(18),typat(natom)
 real(dp),intent(in) :: spinat(3,natom), rprimd(3,3)
 real(dp),intent(in) :: ratsph(ntypat),rhor(nfft,nspden),xred(3,natom)
 type(MPI_type),intent(in) :: mpi_enreg

!Local variables-------------------------------
!scalars
 integer :: iatom,ii
 real(dp) :: intgden_proj, Epen,Econstr,lVp, norm
!arrays
 real(dp) :: intmm(3), mag_1atom(3)
 real(dp), allocatable :: intgden(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),ucvol
 real(dp) :: spinat_norm(3,natom)
 character(len=500) :: message

! *********************************************************************

!We need the metric because it is needed in calcdensph.F90
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 ABI_ALLOCATE (intgden, (nspden,natom))

!We need the integrated magnetic moments
 call calcdensph(gmet,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,std_out,ratsph,rhor,rprimd,typat,ucvol,xred,&
& intgden)

 Epen=0
 Econstr=0
 lVp=0

!Loop over atoms
!-------------------------------------------
 do iatom=1,natom

   norm = sqrt(sum(spinat(:,iatom)**2))
   spinat_norm(:,iatom) = zero
   if (norm > tol10) then
     spinat_norm(:,iatom) = spinat(:,iatom) / norm
   else if (magconon == 1) then
!    if spinat = 0 and we are imposing the direction only, skip this atom
     cycle
   end if
!  Calculate the scalar product of the fixed mag. mom. vector and calculated mag. mom. vector
!  This is actually the size of the projection of the calc. mag. mom. vector on the fixed mag. mom. vector

! for the collinear spin case, set up a fictitious 3D vector along z
   if (nspden == 4) then
     mag_1atom(1:3) = intgden(2:4,iatom)
   else if (nspden == 2) then
     mag_1atom = zero
     mag_1atom(3) = intgden(1,iatom)-intgden(2,iatom)
   end if

   intgden_proj = zero
   intmm = zero
!  Calculate the square bracket term
   if (magconon==1) then
     intgden_proj=spinat_norm(1,iatom)*mag_1atom(1)+ &
&     spinat_norm(2,iatom)*mag_1atom(2)+ &
&     spinat_norm(3,iatom)*mag_1atom(3)

     do ii=1,3
       intmm(ii)=mag_1atom(ii)-spinat_norm(ii,iatom)*intgden_proj
     end do

!    Calculate the energy Epen corresponding to the constraining potential
!    Econstr and lVp do not have a clear meaning (yet)
     Epen=Epen+magcon_lambda*(intmm(1)*intmm(1)+intmm(2)*intmm(2)+intmm(3)*intmm(3))
     Econstr=Econstr-magcon_lambda*(intmm(1)*mag_1atom(1)+intmm(2)*mag_1atom(2)+intmm(3)*mag_1atom(3))
     lVp=lVp+2*magcon_lambda*(intmm(1)*mag_1atom(1)+intmm(2)*mag_1atom(2)+intmm(3)*mag_1atom(3))

   else if (magconon==2) then
     do ii=1,3
       intmm(ii)=mag_1atom(ii)-spinat(ii,iatom)
     end do

!    Calculate the energy Epen corresponding to the constraining potential
!    Epen = -Econstr - lVp
!    Econstr = -M**2 + spinat**2 
!    lVp = +2 M \cdot spinat
     Epen=Epen+magcon_lambda*(intmm(1)*intmm(1)+intmm(2)*intmm(2)+intmm(3)*intmm(3))
     Econstr=Econstr-magcon_lambda*(mag_1atom(1)*mag_1atom(1)+&
&     mag_1atom(2)*mag_1atom(2)+&
&     mag_1atom(3)*mag_1atom(3)) &
&     +magcon_lambda*(spinat(1,iatom)*spinat(1,iatom)+&
&     spinat(2,iatom)*spinat(2,iatom)+&
&     spinat(3,iatom)*spinat(3,iatom))
     lVp=lVp+2*magcon_lambda*(intmm(1)*mag_1atom(1)+intmm(2)*mag_1atom(2)+intmm(3)*mag_1atom(3))
   end if

   write(message, *) 'atom             constraining magnetic field'
   call wrtout(std_out,message,'COLL')
   write(message, '(I3,A2,E12.5,A2,E12.5,A2,E12.5)') &
   iatom,'  ',magcon_lambda*intmm(1),'  ',magcon_lambda*intmm(2),'  ',magcon_lambda*intmm(3)
   call wrtout(std_out,message,'COLL')

!  End loop over atoms
!  -------------------------------------------
 end do

!Printing
 write(message, '(A17,E10.3)' ) ' magcon_lambda    = ',magcon_lambda
 call wrtout(std_out,message,'COLL')
 write(message, '(A17,E12.5)' ) ' Lagrange penalty = ',Epen
 call wrtout(std_out,message,'COLL')
 write(message, '(A17,E12.5)' ) ' E_constraint     = ',Econstr
 call wrtout(std_out,message,'COLL')
 write(message, '(A17,E12.5)' ) ' lVp = ',lVp
 call wrtout(std_out,message,'COLL')

 ABI_DEALLOCATE (intgden)

end subroutine mag_constr_e

!!***
