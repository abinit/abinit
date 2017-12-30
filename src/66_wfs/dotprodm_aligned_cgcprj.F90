!{\src2tex{textfont=tt}}
!!****f* ABINIT/dotprodm_sumdiag_cgcprj
!!
!! NAME
!! dotprodm_sumdiag_cgcprj
!!
!! FUNCTION
!! For one k point and spinpol, compute the matrix of sum of diagonal scalar products 
!! between sets of nband wavefunctions that are known in the cg+cprj representation
!! The sets of wavefunctions must be contained in a big array of sets of wavefunctions.
!! Sij=Sum_m <wfm(set i)|wfm(set j)>
!!
!! Can also treat the case of the computation of scalar products within one set of wavefunctions
!!
!! This implementation is NOT band-parallelized
!! Also, it is far of being optimal at the level of linear algebra, and involves extra copying
!! that are detrimental for performance...
!!
!! COPYRIGHT
!! Copyright (C) 2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cg_set(2,mcg,mset)= plane wave wavefunction coefficients for several sets of wavefunctions (all k points and spinpol)
!!  cprj_set(natom,mcprj,mset) <type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk> with NL projectors in the different sets
!!  dimcprj(natom)=number of lmn components in the <p_{lmn}^i|\psi> for the i-th atom
!!  ibg=shift in cprj_set array to locate current k-point 
!!  icg=shift in cg_set array to locate current k-point
!!  ikpt=current k point index
!!  isppol=current spin polarization index
!!  istwf=input option parameter that describes the storage of wfs
!!  mband=maximum number of bands (used in the dimensioning of cprj_set)
!!  mcg=second dimension of cg array (mpw*nspinor*mband*mkmem*nsppol)
!!  mcprj=second dimension of cprj array 
!!  mkmem=number of k points which can fit in memory
!!  mpi_enreg=information about MPI parallelization
!!  mset=third dimension of cg_set and cprj_set, maximum number of sets 
!!  natom=number of atoms
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nbd=number of bands for each set of wavefunctions
!!  npw=number of planewaves in basis at this k point
!!  nset1=number of sets of wavefunctions to be considered in the left side of the scalar products
!!  nset2=number of sets of wavefunctions to be considered in the right side of the scalar products
!!  nspinor=number of spinor components
!!  nsppol=number of spin polarizations
!!  ntypat=number of types of atoms
!!  pawtab(dtset%ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  shift_set1=shift that defines the first set of wavefunctions to be considered in the left side of the scalar products
!!  shift_set2=shift that defines the first set of wavefunctions to be considered in the right side of the scalar products
!!  usepaw=1 if PAW is activated
!!
!! OUTPUT
!!  smn(2,nset1,nset2)=matrix of sum of diagonal scalar products between the first set of wavefunctions and the second set of wavefunctions
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dotprod_sumdiag_cgcprj(atindx1,cg_set,cprj_set,dimcprj,&
& ibg,icg,ikpt,isppol,istwf,mband,mcg,mcprj,mkmem,&
& mpi_enreg,mset,natom,nattyp,nbd,npw,nset1,nset2,nspinor,nsppol,ntypat,&
& shift_set1,shift_set2,pawtab,smn,usepaw)

 use defs_basis
 use defs_abitypes
 use m_cgtools
 use m_errors
 use m_xmpi
 use m_pawtab, only : pawtab_type
 use m_pawcprj, only : pawcprj_type, pawcprj_alloc, pawcprj_get, pawcprj_free

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dotprod_sumdiag_cgcprj'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ibg,icg,ikpt,isppol,istwf
 integer, intent(in) :: mband,mcg,mcprj,mkmem,mset,
 integer, intent(in) :: natom,nbd,npw,nset1,nset2,nspinor,nsppol,ntypat
 integer, intent(in) :: shift_set1,shift_set2,usepaw
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer, intent(in) :: atindx1(natom),dimcprj(natom),nattyp(ntypat)
 real(dp), intent(in) :: cg_set(2,mcg,mset)
 real(dp), intent(out) :: smn(2,nset1,nset2)
 type(pawcprj_type),intent(in) :: cprj_set(natom,mcprj,mset)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: ia,iat,itypat,ibd,icgb,ig
 integer :: ilmn1,ilmn2,ind_set1,ind_set2,iset1,iset2,klmn
 real(dp) :: dotr,doti
!arrays
 real(dp),allocatable :: cwavef1(:,:),cwavef2(:,:)
 type(pawcprj_type),allocatable :: cprj1_k(:,:),cprj2_k(:,:)

! *************************************************************************

!DEBUG
 write(std_out,*)' dotprod_sumdiag_cgcprj : enter '
!ENDDEBUG

 ABI_ALLOCATE(cwavef1,(2,npw*nspinor))
 ABI_ALLOCATE(cwavef2,(2,npw*nspinor))
 if(usepaw==1) then
   ABI_DATATYPE_ALLOCATE(cprj1_k,(natom,nspinor*nbd))
   ABI_DATATYPE_ALLOCATE(cprj2_k,(natom,nspinor*nbd))
 endif
 if(usepaw==1) then
   call pawcprj_alloc(cprj1_k,cprj_set(1,1,1)%ncpgr,dimcprj)
 endif

 smn(:,:,:)=zero

 icgb=icg
 do ibd=1,nbd

   do iset1=1,nset1

     ind_set1=iset1+shift_set1

!    Extract wavefunction information
     do ig=1,npw*nspinor
       cwavef1(1,ig)=cg_set(1,ig+icgb,ind_set1)
       cwavef1(2,ig)=cg_set(2,ig+icgb,ind_set1)
     end do
     if(usepaw==1) then
       call pawcprj_get(atindx1,cprj1_k,cprj_set(:,:,ind_set1),natom,1,ibg,ikpt,1,isppol,mband,&
&         mkmem,natom,nbd,nbd,nspinor,nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
     endif

     do iset2=1,nset2

       ind_set2=iset2+shift_set2
       if(ind_set2<ind_set1 .and. ind_set2>shift_iset1)then
         continue ! These matrix elements have already been computed, the smn matrix will be completed later.

       else if(ind_set1/=ind_set2)then

!        Extract wavefunction information
         do ig=1,npw*nspinor
           cwavef2(1,ig)=cg_set(1,ig+icgb,ind_set2)
           cwavef2(2,ig)=cg_set(2,ig+icgb,ind_set2)
         end do

         if(usepaw==1) then
           call pawcprj_get(atindx1,cprj2_k,cprj_set(:,:,ind_set2),natom,1,ibg,ikpt,1,isppol,mband,&
&           mkmem,natom,nbd,nbd,nspinor,nsppol,0,&
&           mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         endif

!        Calculate Smn=<cg1|cg2>
         call dotprod_g(dotr,doti,istwf,npw*nspinor,2,cwavef1,cwavef2,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

         if(usepaw==1) then
           ia =0
           do itypat=1,ntypat
             do iat=1+ia,nattyp(itypat)+ia
               do ilmn1=1,pawtab(itypat)%lmn_size
                 do ilmn2=1,ilmn1
                   klmn=((ilmn1-1)*ilmn1)/2+ilmn2
                   dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj1_k(iat,ibd)%cp(1,ilmn1)*cprj2_k(iat,ibd)%cp(1,ilmn2)+&
&                   cprj1_k(iat,ibd)%cp(2,ilmn1)*cprj2_k(iat,ibd)%cp(2,ilmn2))
                   doti=doti+pawtab(itypat)%sij(klmn)*(cprj1_k(iat,ibd)%cp(1,ilmn1)*cprj2_k(iat,ibd)%cp(2,ilmn2)-&
&                   cprj1_k(iat,ibd)%cp(2,ilmn1)*cprj2_k(iat,ibd)%cp(1,ilmn2))
                 end do
                 do ilmn2=ilmn1+1,pawtab(itypat)%lmn_size
                   klmn=((ilmn2-1)*ilmn2)/2+ilmn1
                   dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj1_k(iat,ibd)%cp(1,ilmn1)*cprj2_k(iat,ibd)%cp(1,ilmn2)+&
&                   cprj1_k(iat,ibd)%cp(2,ilmn1)*cprj2_k(iat,ibd)%cp(2,ilmn2))
                   doti=doti+pawtab(itypat)%sij(klmn)*(cprj1_k(iat,ibd)%cp(1,ilmn1)*cprj2_k(iat,ibd)%cp(2,ilmn2)-&
&                   cprj1_k(iat,ibd)%cp(2,ilmn1)*cprj2_k(iat,ibd)%cp(1,ilmn2))
                 end do
               end do
             end do
             ia=ia+nattyp(itypat)
           end do

       else 
!        Diagonal part : no need to extract another wavefunction, and the scalar product must be real

!        Calculate Smn=<cg1|cg1>
         call dotprod_g(dotr,doti,istwf,npw*nspinor,1,cwavef1,cwavef1,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

         if(usepaw==1) then
           ia =0
           do itypat=1,ntypat
             do iat=1+ia,nattyp(itypat)+ia
               do ilmn1=1,pawtab(itypat)%lmn_size
                 do ilmn2=1,ilmn1
                   klmn=((ilmn1-1)*ilmn1)/2+ilmn2
                   dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj1_k(iat,ibd)%cp(1,ilmn1)*cprj1_k(iat,ibd)%cp(1,ilmn2)+&
&                   cprj1_k(iat,ibd)%cp(2,ilmn1)*cprj1_k(iat,ibd)%cp(2,ilmn2))
                 end do
                 do ilmn2=ilmn1+1,pawtab(itypat)%lmn_size
                   klmn=((ilmn2-1)*ilmn2)/2+ilmn1
                   dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj1_k(iat,ibd)%cp(1,ilmn1)*cprj1_k(iat,ibd)%cp(1,ilmn2)+&
&                   cprj1_k(iat,ibd)%cp(2,ilmn1)*cprj1_k(iat,ibd)%cp(2,ilmn2))
                 end do
               end do
             end do
             ia=ia+nattyp(itypat)
           end do

       endif
 
       smn(1,iset1,iset2)=smn(1,iset1,iset2)+dotr
       smn(2,iset1,iset2)=smn(2,iset1,iset2)+doti

     enddo ! iset2
   enddo ! iset1

!  End loop over bands ibd 
   icgb=icgb+npw*nspinor
 end do

!Complete the matrix, using its hermitian property.
 do iset1=1,nset1
   ind_set1=iset1+shift_iset1
   do iset2=1,nset2
     ind_set2=iset2+shift_iset2
     if(ind_set2<ind_set1 .and. ind_set2>shift_iset1)then
       smn(1,iset1,iset2)= smn(1,iset2,iset1)
       smn(2,iset1,iset2)=-smn(2,iset2,iset1)
     endif
   enddo
 enddo

!DEBUG
!write(std_out,*)' smn=',smn
!ENDDEBUG

 ABI_DEALLOCATE(cwavef1)
 ABI_DEALLOCATE(cwavef2)
 if(usepaw==1) then
   call pawcprj_free(cprj1_k)
   call pawcprj_free(cprj2_k)
   ABI_DATATYPE_DEALLOCATE(cprj1_k)
   ABI_DATATYPE_DEALLOCATE(cprj2_k)
 end if

end subroutine dotprod_sumdiag_cgcprj
!!***
