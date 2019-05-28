!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paw_onsite
!! NAME
!!  m_paw_onsite
!!
!! FUNCTION
!!  This module contains a set of routines to compute various PAW on-site quantities
!!  i.e. quantities expressed with <Phi_i|...|Phi_j> and/or <tild_Phi_i|...|tild_Phi_j>.
!!
!! COPYRIGHT
!! Copyright (C) 2013-2019 ABINIT group (MT,FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

MODULE m_paw_onsite

 USE_DEFS
 USE_MSG_HANDLING
 USE_MEMORY_PROFILING

 use m_pawrad,      only : pawrad_type, pawrad_deducer0, simp_gen, nderiv_gen
 use m_pawtab,      only : pawtab_type
 use m_paw_sphharm, only : setnabla_ylm

 implicit none

 private

!public procedures.
 public ::  pawnabla_init      ! Evaluate valence-valence on-site contribs of the nabla operator in cart. coord.
 public ::  pawnabla_core_init ! Evaluate core-valence on-site contribs of the nabla operator in cart. coord.

!!***

CONTAINS

!===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_onsite/pawnabla_init
!! NAME
!! pawnabla_init
!!
!! FUNCTION
!! Evaluate all valence-valence onsite contributions of the nabla operator in cartesian coordinates.
!! Store values in the pawtab% data structure.
!!
!! INPUTS
!!  mpsang=1+maximum angular momentum
!!  ntypat=Number of types of atoms in cell
!!  Pawrad(ntypat)<Pawrad_type>=PAW radial mesh and related data:
!!    %rad(mesh_size)=The coordinates of all the points of the radial mesh
!!  Pawtab(ntypat) <type(pawtab_type>=PAW tabulated starting data:
!!    %mesh_size=Dimension of radial mesh
!!    %lmn_size=Number of (l,m,n) elements for the PAW basis
!!
!! OUTPUT
!!  See side effects
!!
!! SIDE EFFECTS
!!  Pawtab(ntypat) <type(pawtab_type>=PAW tabulated starting data:
!!    %has_nabla=set to 1 in matrix elements are calculated and stored
!!    %nabla_ij(3,lmn_size,lmn_size)= <phi_i|nabla|phi_j>-<tphi_i|nabla|tphi_j>
!!
!! NOTES
!!  MG extracted this piece of code from optics_paw.F90 in order to have something more 
!!  reusable! Note however the storage mode of nabla_ij differs from optics_paw 
!!  (here Cartesian coordinates run faster). Besides nabla_ij contains the matrix 
!!  elements of \nabla instead of the elements of the momentum operator p.
!!
!! PARENTS
!!
!! CHILDREN
!!      setnabla_ylm,nderiv_gen,pawrad_deducer0,simp_gen
!!
!! SOURCE

subroutine pawnabla_init(mpsang,ntypat,pawrad,pawtab)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpsang,ntypat
!arrays
 type(pawtab_type),target,intent(inout) :: pawtab(ntypat)
 type(pawrad_type),intent(in) :: pawrad(ntypat)

!Local variables-------------------------------
!scalars
 integer :: nln,il,ilm,ilmn,iln,itypat
 integer :: jl,jlm,jlmn,jln,lmn_size,mesh_size 
 real(dp) :: intg
 character(len=500) :: msg      
!arrays
 integer, LIBPAW_CONTIGUOUS pointer :: indlmn(:,:)
 real(dp) :: ang_phipphj(mpsang**2,mpsang**2,8)
 real(dp),allocatable :: dphi(:),dtphi(:),ff(:),int1(:,:),int2(:,:),rad(:)
 
! *************************************************************************

 if (mpsang>4)then
   write(msg,'(3a)')&
&   'Not designed for angular momentum greater than 3 ',ch10,&
&   'Modification in the table defined in routine setnabla_ylm is required.'
   MSG_BUG(msg)
 end if

!Integration of the angular part: all angular integrals have been computed 
!outside Abinit and tabulated for each (l,m) value up to l=3
 call setnabla_ylm(ang_phipphj,mpsang)

 do itypat=1,ntypat

!  COMPUTE nabla_ij := <phi_i|nabla|phi_j>-<tphi_i|nabla|tphi_j> for this type
   mesh_size=pawtab(itypat)%mesh_size
   lmn_size=pawtab(itypat)%lmn_size
   nln=pawtab(itypat)%basis_size

   if (allocated(pawtab(itypat)%nabla_ij)) then
     LIBPAW_DEALLOCATE(pawtab(itypat)%nabla_ij)
   end if
   LIBPAW_ALLOCATE(pawtab(itypat)%nabla_ij,(3,lmn_size,lmn_size))
   pawtab(itypat)%has_nabla=1

   LIBPAW_ALLOCATE(ff,(mesh_size))
   LIBPAW_ALLOCATE(rad,(mesh_size))
   LIBPAW_ALLOCATE(int2,(lmn_size,lmn_size))
   LIBPAW_ALLOCATE(int1,(lmn_size,lmn_size))
   LIBPAW_ALLOCATE(dphi,(mesh_size))
   LIBPAW_ALLOCATE(dtphi,(mesh_size))

   indlmn => pawtab(itypat)%indlmn
   rad(1:mesh_size)=pawrad(itypat)%rad(1:mesh_size)

!  int1=\int phi phj/r dr - \int tphi tphj /r dr
   do jln=1,nln
     do iln=1,nln
       ff(2:mesh_size)= ( &
&       pawtab(itypat)%phi (2:mesh_size,iln)*pawtab(itypat)%phi (2:mesh_size,jln) &
&       -pawtab(itypat)%tphi(2:mesh_size,iln)*pawtab(itypat)%tphi(2:mesh_size,jln) ) /rad(2:mesh_size)
       call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intg,ff,pawrad(itypat))
       int1(iln,jln)=intg
     end do
   end do

!  int2=\int phi/r d/dr(phj/r) r^2dr - \int tphi/r d/dr(tphj/r)r^2 dr
   do jln=1,nln
     ff(1:mesh_size)=pawtab(itypat)%phi(1:mesh_size,jln)
     call nderiv_gen(dphi,ff,pawrad(itypat))
     ff(1:mesh_size)=pawtab(itypat)%tphi(1:mesh_size,jln)
     call nderiv_gen(dtphi,ff,pawrad(itypat))

     do iln=1,nln
       ff(2:mesh_size)= &
&       pawtab(itypat)%phi (2:mesh_size,iln)*dphi (2:mesh_size) &
&       -pawtab(itypat)%phi (2:mesh_size,iln)*pawtab(itypat)%phi (2:mesh_size,jln)/rad(2:mesh_size) &
&       -( pawtab(itypat)%tphi(2:mesh_size,iln)*dtphi(2:mesh_size) &
&       -pawtab(itypat)%tphi(2:mesh_size,iln)*pawtab(itypat)%tphi(2:mesh_size,jln)/rad(2:mesh_size) )
       call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intg,ff,pawrad(itypat))
       int2(iln,jln)=intg
     end do
   end do

!  Integration of the radial part, Note unpacked loop
   do jlmn=1,lmn_size
     jlm=indlmn(4,jlmn)
     jl =indlmn(5,jlmn)
     do ilmn=1,lmn_size
       ilm=indlmn(4,ilmn)
       il =indlmn(5,ilmn)

       pawtab(itypat)%nabla_ij(1,ilmn,jlmn)= &
&       int2(il,jl)* ang_phipphj(ilm,jlm,1) &
&       +int1(il,jl)*(ang_phipphj(ilm,jlm,2)+ang_phipphj(ilm,jlm,3))

       pawtab(itypat)%nabla_ij(2,ilmn,jlmn)= &
&       int2(il,jl)* ang_phipphj(ilm,jlm,4) &
&       +int1(il,jl)*(ang_phipphj(ilm,jlm,5)+ang_phipphj(ilm,jlm,6))

       pawtab(itypat)%nabla_ij(3,ilmn,jlmn)= &
&       int2(il,jl)* ang_phipphj(ilm,jlm,7) &
&       +int1(il,jl)* ang_phipphj(ilm,jlm,8)

     end do !ilmn
   end do !jlmn

   pawtab(itypat)%has_nabla=2
   LIBPAW_DEALLOCATE(ff)
   LIBPAW_DEALLOCATE(rad)
   LIBPAW_DEALLOCATE(int2)
   LIBPAW_DEALLOCATE(int1)
   LIBPAW_DEALLOCATE(dphi)
   LIBPAW_DEALLOCATE(dtphi)

 end do !itypat

end subroutine pawnabla_init
!!***

!----------------------------------------------------------------------

!!****f* m_paw_onsite/pawnabla_core_init
!! NAME
!! pawnabla_core_init
!!
!! FUNCTION
!! Evaluate core-valence onsite contributions of the nabla operator in cartesian coordinates.
!! Core wave-functions are only given for one atom type.
!! Store values in the pawtab% data structure.
!!
!! INPUTS
!!  mpsang=1+maximum angular momentum
!!  ntypat=Number of types of atoms in cell
!!  Pawrad(ntypat)<Pawrad_type>=PAW radial mesh and related data:
!!    %rad(mesh_size)=The coordinates of all the points of the radial mesh
!!  Pawtab(ntypat) <type(pawtab_type>=PAW tabulated starting data:
!!    %mesh_size=Dimension of radial mesh
!!    %lmn_size=Number of (l,m,n) elements for the PAW basis
!!  phi_cor(mesh_size,nphicor)=core wave-functions for the current type of atoms.
!!  indlmn_cor(6,nlmn_core)= array giving l,m,n,lm,ln,s for i=lmn, for the core wave functions.
!!
!! OUTPUT
!!  See side effects
!!
!! SIDE EFFECTS
!!  Pawtab(ntypat) <type(pawtab_type>=PAW tabulated starting data:
!!    %has_nabla=set to 1 in matrix elements are calculated and stored
!!    %nabla_ij(3,lmn_size,lmn_size)= <phi_i|nabla|phi_j>-<tphi_i|nabla|tphi_j>
!!
!! NOTES
!!  MG extracted this piece of code from optics_paw.F90 in order to have something more
!!  reusable! Note however the storage mode of nabla_ij differs from optics_paw
!!  (here Cartesian coordinates run faster). Besides nabla_ij contains the matrix
!!  elements of \nabla instead of the elements of the momentum operator p.
!!
!! PARENTS
!!
!! CHILDREN
!!      setnabla_ylm,nderiv_gen,pawrad_deducer0,simp_gen
!!
!! SOURCE

subroutine pawnabla_core_init(mpsang,ntypat,pawrad,pawtab,phi_cor,indlmn_cor)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpsang,ntypat
!arrays
 integer,intent(in) :: indlmn_cor(:,:)
 real(dp) :: phi_cor(:,:)
 type(pawtab_type),target,intent(inout) :: pawtab(ntypat)
 type(pawrad_type),intent(in) :: pawrad(ntypat)

!Local variables-------------------------------
!scalars
 integer :: nln,nln_cor,il,ilm,ilmn,iln,itypat
 integer :: jl,jlm,jlmn,jln,lmn_size,lmcmax,lmncmax,mesh_size,mesh_size_cor
 real(dp) :: intg
 character(len=500) :: msg
!arrays
 integer, LIBPAW_CONTIGUOUS pointer :: indlmn(:,:)
 real(dp) , allocatable:: ang_phipphj(:,:,:)
 real(dp),allocatable :: dphi(:),ff(:),int1(:,:),int2(:,:),rad(:)

! *************************************************************************

 mesh_size_cor=size(phi_cor,1)
 nln_cor=size(phi_cor,2)
 lmncmax=size(indlmn_cor,2)
 lmcmax=indlmn_cor(4,lmncmax)
 LIBPAW_ALLOCATE(ang_phipphj,(lmcmax,lmcmax,8))

 if (mpsang>4)then
   write(msg,'(3a)')&
&   'Not designed for angular momentum greater than 3 ',ch10,&
&   'Modification in the table defined in routine setnabla_ylm is required.'
   MSG_BUG(msg)
 end if

 if (mesh_size_cor/=pawrad(1)%mesh_size) then
   write(msg,'(a)') 'Wrong mesh_size_cor value (1)!'
   MSG_BUG(msg)
 end if

 if (any(mesh_size_cor/=pawtab(:)%mesh_size)) then
   write(msg,'(3a)') 'Wrong mesh_size_cor value (2)!',ch10,&
&                    'Should have only one type of atom.'
   MSG_ERROR(msg)
 end if

!Integration of the angular part: all angular integrals have been computed
!outside Abinit and tabulated for each (l,m) value up to l=3
 call setnabla_ylm(ang_phipphj,mpsang)

 do itypat=1,ntypat

!  COMPUTE nabla_ij := <phi_i|nabla|phi_cor> for this type
   mesh_size=pawtab(itypat)%mesh_size
   lmn_size=pawtab(itypat)%lmn_size
   nln=pawtab(itypat)%basis_size

   if (allocated(pawtab(itypat)%nabla_ij)) then
     LIBPAW_DEALLOCATE(pawtab(itypat)%nabla_ij)
   end if
   LIBPAW_ALLOCATE(pawtab(itypat)%nabla_ij,(3,lmn_size,lmncmax))
   pawtab(itypat)%has_nabla=1

   LIBPAW_ALLOCATE(ff,(mesh_size))
   LIBPAW_ALLOCATE(rad,(mesh_size))
   LIBPAW_ALLOCATE(int2,(lmn_size,lmncmax))
   LIBPAW_ALLOCATE(int1,(lmn_size,lmncmax))
   LIBPAW_ALLOCATE(dphi,(mesh_size))

   indlmn => pawtab(itypat)%indlmn
   rad(1:mesh_size)=pawrad(itypat)%rad(1:mesh_size)

!  int1=\int phi phi_core/r dr
   do jln=1,nln_cor
     do iln=1,nln
       ff(2:mesh_size)=(pawtab(itypat)%phi(2:mesh_size,iln)*phi_cor(2:mesh_size,jln))/rad(2:mesh_size)
       call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intg,ff,pawrad(itypat))
       int1(iln,jln)=intg
     end do
   end do

!  int2=\int phi/r d/dr(phi_core/r) r^2dr - phi phi_core/r dr
   do jln=1,nln_cor
     ff(1:mesh_size)=phi_cor(1:mesh_size,jln)
     call nderiv_gen(dphi,ff,pawrad(itypat))
     do iln=1,nln
       ff(2:mesh_size)=pawtab(itypat)%phi(2:mesh_size,iln)*dphi(2:mesh_size) &
&       -pawtab(itypat)%phi(2:mesh_size,iln)*phi_cor(2:mesh_size,jln)/rad(2:mesh_size)
       call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intg,ff,pawrad(itypat))
       int2(iln,jln)=intg
     end do
   end do

!  Integration of the radial part, Note unpacked loop
   do jlmn=1,lmncmax
     jlm=indlmn_cor(4,jlmn)
     jl =indlmn_cor(5,jlmn)
     do ilmn=1,lmn_size
       ilm=indlmn(4,ilmn)
       il =indlmn(5,ilmn)

       pawtab(itypat)%nabla_ij(1,ilmn,jlmn)= &
&       int2(il,jl)* ang_phipphj(ilm,jlm,1) &
&       +int1(il,jl)*(ang_phipphj(ilm,jlm,2)+ang_phipphj(ilm,jlm,3))

       pawtab(itypat)%nabla_ij(2,ilmn,jlmn)= &
&       int2(il,jl)* ang_phipphj(ilm,jlm,4) &
&       +int1(il,jl)*(ang_phipphj(ilm,jlm,5)+ang_phipphj(ilm,jlm,6))

       pawtab(itypat)%nabla_ij(3,ilmn,jlmn)= &
&       int2(il,jl)* ang_phipphj(ilm,jlm,7) &
&       +int1(il,jl)* ang_phipphj(ilm,jlm,8)

     end do !ilmn
   end do !jlmn

   pawtab(itypat)%has_nabla=3
   LIBPAW_DEALLOCATE(ff)
   LIBPAW_DEALLOCATE(rad)
   LIBPAW_DEALLOCATE(int2)
   LIBPAW_DEALLOCATE(int1)
   LIBPAW_DEALLOCATE(dphi)

 end do !itypat
 LIBPAW_DEALLOCATE(ang_phipphj)

end subroutine pawnabla_core_init
!!***

!----------------------------------------------------------------------

end module m_paw_onsite
!!***

