!!****m* ABINIT/m_cgcprj
!! NAME
!!  m_cgcprj
!!
!! FUNCTION
!!  Functions operating on wavefunctions in the cg+cprj representation.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_cgcprj

 use defs_basis
 use defs_abitypes
 use m_cgtools
 use m_errors
 use m_xmpi

 use m_pawtab, only : pawtab_type
 use m_pawcprj, only : pawcprj_type, pawcprj_alloc, pawcprj_get, pawcprj_free, pawcprj_lincom

 implicit none

 private
!!***

 public :: dotprod_set_cgcprj
 public :: dotprodm_sumdiag_cgcprj
 public :: lincom_cgcprj
 public :: cgcprj_cholesky
!!***

contains
!!***

!!****f* ABINIT/dotprod_set_cgcprj
!! NAME
!! dotprod_set_cgcprj
!!
!! FUNCTION
!! For one k point and spinpol, compute the matrix of scalar products between two sets of nband wavefunctions
!! that are known in the cg+cprj representation
!! Smn=<wf1m|wf2n>
!!
!! Can also treat the case of the computation of scalar products within one set of nband wavefunctions
!! without recomputing already computed matrix elements (define hermitian=1)
!!
!! This implementation is NOT band-parallelized
!! Also, it is far of being optimal at the level of linear algebra, and involves extra copying
!! that are detrimental for performance...
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cg1(2,mcg1)= plane wave wavefunction coefficients for the first set of wavefunctions (all k points and spinpol)
!!  cg2(2,mcg2)= plane wave wavefunction coefficients for the second set of wavefunctions (all k points and spinpol)
!!  cprj1(natom,mcprj1) <type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk> with NL projectors in the first set,
!!  cprj2(natom,mcprj2) <type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk> with NL projectors in the second set,
!!  dimcprj(natom)=number of lmn components in the <p_{lmn}^i|\psi> for the i-th atom
!!  hermitian= if 1, consider that the Smn matrix is hermitian, and do not recompute already computed matrix elements.
!!  ibg1=shift in cprj1 array to locate current k-point.
!!    Might be 0, in which case cprj1 is not copied internally, which saves some time/space
!!  ibg2=shift in cprj2 array to locate current k-point.
!!    Might be 0, in which case cprj2 is not copied internally, which saves some time/space
!!  icg1=shift in cg1 array to locate current k-point
!!  icg2=shift in cg2 array to locate current k-point
!!  ikpt=current k point index
!!  isppol=current spin polarization index
!!  istwf=input option parameter that describes the storage of wfs
!!  mcg1=second dimension of cg1 array (mpw*nspinor*mband1*mkmem*nsppol)
!!  mcg2=second dimension of cg2 array (mpw*nspinor*mband2*mkmem*nsppol)
!!  mcprj1=second dimension of cprj1 array
!!  mcprj2=second dimension of cprj2 array
!!  mkmem=number of k points which can fit in memory
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nbd1=number of bands for the first set of wavefunctions
!!  nbd2=number of bands for the second set of wavefunctions
!!  npw=number of planewaves in basis at this k point
!!  nspinor=number of spinor components
!!  nsppol=number of spin polarizations
!!  ntypat=number of types of atoms
!!  pawtab(dtset%ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  usepaw=1 if PAW is activated
!!
!! OUTPUT
!!  smn(2,nbd1,nbd2)=matrix of scalar products between the first set of wavefunctions and the second set of wavefunctions
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_cgcprj,m_extraprho,m_scfcv_core
!!
!! CHILDREN
!!      dotprod_set_cgcprj,lincom_cgcprj,zpotrf,ztrsm
!!
!! SOURCE

subroutine dotprod_set_cgcprj(atindx1,cg1,cg2,cprj1,cprj2,dimcprj,hermitian,&
& ibg1,ibg2,icg1,icg2,ikpt,isppol,istwf,mband,mcg1,mcg2,mcprj1,mcprj2,mkmem,&
& mpi_enreg,natom,nattyp,nbd1,nbd2,npw,nspinor,nsppol,ntypat,pawtab,smn,usepaw)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: hermitian,ibg1,ibg2,icg1,icg2,ikpt,isppol,istwf
 integer, intent(in) :: mkmem,mband,mcg1,mcg2,mcprj1,mcprj2
 integer, intent(in) :: natom,nbd1,nbd2,npw,nspinor,nsppol,ntypat,usepaw
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer, intent(in) :: atindx1(natom),dimcprj(natom),nattyp(ntypat)
 real(dp), intent(in) :: cg1(2,mcg1),cg2(2,mcg2)
 real(dp), intent(out) :: smn(2,nbd1,nbd2)
 type(pawcprj_type),intent(in) :: cprj1(natom,mcprj1),cprj2(natom,mcprj2)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: ia,iat,itypat,ibd1,ibd2,icgb1,icgb2,ier,ig,ii,i1,i2,iorder
 integer :: ilmn1,ilmn2,klmn,max_nbd2,nbd
 real(dp) :: dotr,doti
!arrays
 real(dp),allocatable :: cwavef1(:,:),cwavef2(:,:),proj(:,:,:)
 real(dp),allocatable :: eigval(:),eigvec(:,:,:),matrx(:,:),zhpev1(:,:),zhpev2(:)
 type(pawcprj_type),allocatable :: cprj1_k(:,:),cprj2_k(:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' dotprod_set_cgcprj : enter '
!write(std_out,*)' dotprod_set_cgcprj : npw, nspinor=',npw,nspinor
!write(std_out,*)' dotprod_set_cgcprj : usepaw,nbd1,nbd2=',usepaw,nbd1,nbd2
!call flush(std_out)
!ENDDEBUG

 if(hermitian==1)then
   if(nbd1/=nbd2)then
     ABI_ERROR(' With hermitian==1, nb1 and nb2 must be equal ')
   end if
 end if

 ABI_ALLOCATE(cwavef1,(2,npw*nspinor))
 ABI_ALLOCATE(cwavef2,(2,npw*nspinor))
 if(usepaw==1) then
   ABI_DATATYPE_ALLOCATE(cprj1_k,(natom,nspinor*nbd1))
   ABI_DATATYPE_ALLOCATE(cprj2_k,(natom,nspinor*nbd2))
   iorder=0 ! There is no change of ordering of cprj when copying wavefunctions
 end if
 if(usepaw==1 .and. ibg1/=0) then
   call pawcprj_alloc(cprj1_k,cprj1(1,1)%ncpgr,dimcprj)
 end if
 if(usepaw==1 .and. ibg2/=0) then
   call pawcprj_alloc(cprj2_k,cprj1(1,1)%ncpgr,dimcprj)
 end if

 icgb1=icg1
 do ibd1=1,nbd1

!  Extract wavefunction information
   do ig=1,npw*nspinor
     cwavef1(1,ig)=cg1(1,ig+icgb1)
     cwavef1(2,ig)=cg1(2,ig+icgb1)
   end do
   if(usepaw==1 .and. ibg1/=0) then
     call pawcprj_get(atindx1,cprj1_k,cprj1,natom,1,ibg1,ikpt,iorder,isppol,mband,&
&     mkmem,natom,nbd1,nbd1,nspinor,nsppol,0,&
&     mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
   end if

   icgb2=icg2
   max_nbd2=nbd2
   if(hermitian==1)max_nbd2=ibd1
   do ibd2=1,max_nbd2

!    XG171222 Note that this copy step, being inside the ibd1 loop, is quite detrimental.
!    It might be reduced by copying several cwavef2, and use a ZGEMM type of approach.

!    Extract wavefunction information
     do ig=1,npw*nspinor
       cwavef2(1,ig)=cg2(1,ig+icgb2)
       cwavef2(2,ig)=cg2(2,ig+icgb2)
     end do

     if(usepaw==1 .and. ibg2/=0) then
       call pawcprj_get(atindx1,cprj2_k,cprj2,natom,1,ibg2,ikpt,iorder,isppol,mband,&
&       mkmem,natom,nbd2,nbd2,nspinor,nsppol,0,&
&       mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
     end if

!    Calculate Smn=<cg1|cg2>
     call dotprod_g(dotr,doti,istwf,npw*nspinor,2,cwavef1,cwavef2,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

     if(usepaw==1) then
       ia =0
       do itypat=1,ntypat
         if(ibg1/=0 .and. ibg2/=0)then
           do iat=1+ia,nattyp(itypat)+ia
             do ilmn1=1,pawtab(itypat)%lmn_size
               do ilmn2=1,ilmn1
                 klmn=((ilmn1-1)*ilmn1)/2+ilmn2
                 dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj1_k(iat,ibd1)%cp(1,ilmn1)*cprj2_k(iat,ibd2)%cp(1,ilmn2)+&
&                 cprj1_k(iat,ibd1)%cp(2,ilmn1)*cprj2_k(iat,ibd2)%cp(2,ilmn2))
                 doti=doti+pawtab(itypat)%sij(klmn)*(cprj1_k(iat,ibd1)%cp(1,ilmn1)*cprj2_k(iat,ibd2)%cp(2,ilmn2)-&
&                 cprj1_k(iat,ibd1)%cp(2,ilmn1)*cprj2_k(iat,ibd2)%cp(1,ilmn2))
               end do
               do ilmn2=ilmn1+1,pawtab(itypat)%lmn_size
                 klmn=((ilmn2-1)*ilmn2)/2+ilmn1
                 dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj1_k(iat,ibd1)%cp(1,ilmn1)*cprj2_k(iat,ibd2)%cp(1,ilmn2)+&
&                 cprj1_k(iat,ibd1)%cp(2,ilmn1)*cprj2_k(iat,ibd2)%cp(2,ilmn2))
                 doti=doti+pawtab(itypat)%sij(klmn)*(cprj1_k(iat,ibd1)%cp(1,ilmn1)*cprj2_k(iat,ibd2)%cp(2,ilmn2)-&
&                 cprj1_k(iat,ibd1)%cp(2,ilmn1)*cprj2_k(iat,ibd2)%cp(1,ilmn2))
               end do
             end do
           end do
         else if(ibg1/=0 .and. ibg2==0)then
           do iat=1+ia,nattyp(itypat)+ia
             do ilmn1=1,pawtab(itypat)%lmn_size
               do ilmn2=1,ilmn1
                 klmn=((ilmn1-1)*ilmn1)/2+ilmn2
                 dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj1_k(iat,ibd1)%cp(1,ilmn1)*cprj2(iat,ibd2)%cp(1,ilmn2)+&
&                 cprj1_k(iat,ibd1)%cp(2,ilmn1)*cprj2(iat,ibd2)%cp(2,ilmn2))
                 doti=doti+pawtab(itypat)%sij(klmn)*(cprj1_k(iat,ibd1)%cp(1,ilmn1)*cprj2(iat,ibd2)%cp(2,ilmn2)-&
&                 cprj1_k(iat,ibd1)%cp(2,ilmn1)*cprj2(iat,ibd2)%cp(1,ilmn2))
               end do
               do ilmn2=ilmn1+1,pawtab(itypat)%lmn_size
                 klmn=((ilmn2-1)*ilmn2)/2+ilmn1
                 dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj1_k(iat,ibd1)%cp(1,ilmn1)*cprj2(iat,ibd2)%cp(1,ilmn2)+&
&                 cprj1_k(iat,ibd1)%cp(2,ilmn1)*cprj2(iat,ibd2)%cp(2,ilmn2))
                 doti=doti+pawtab(itypat)%sij(klmn)*(cprj1_k(iat,ibd1)%cp(1,ilmn1)*cprj2(iat,ibd2)%cp(2,ilmn2)-&
&                 cprj1_k(iat,ibd1)%cp(2,ilmn1)*cprj2(iat,ibd2)%cp(1,ilmn2))
               end do
             end do
           end do
         else if(ibg1==0 .and. ibg2/=0)then
           do iat=1+ia,nattyp(itypat)+ia
             do ilmn1=1,pawtab(itypat)%lmn_size
               do ilmn2=1,ilmn1
                 klmn=((ilmn1-1)*ilmn1)/2+ilmn2
                 dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj1(iat,ibd1)%cp(1,ilmn1)*cprj2_k(iat,ibd2)%cp(1,ilmn2)+&
&                 cprj1(iat,ibd1)%cp(2,ilmn1)*cprj2_k(iat,ibd2)%cp(2,ilmn2))
                 doti=doti+pawtab(itypat)%sij(klmn)*(cprj1(iat,ibd1)%cp(1,ilmn1)*cprj2_k(iat,ibd2)%cp(2,ilmn2)-&
&                 cprj1(iat,ibd1)%cp(2,ilmn1)*cprj2_k(iat,ibd2)%cp(1,ilmn2))
               end do
               do ilmn2=ilmn1+1,pawtab(itypat)%lmn_size
                 klmn=((ilmn2-1)*ilmn2)/2+ilmn1
                 dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj1(iat,ibd1)%cp(1,ilmn1)*cprj2_k(iat,ibd2)%cp(1,ilmn2)+&
&                 cprj1(iat,ibd1)%cp(2,ilmn1)*cprj2_k(iat,ibd2)%cp(2,ilmn2))
                 doti=doti+pawtab(itypat)%sij(klmn)*(cprj1(iat,ibd1)%cp(1,ilmn1)*cprj2_k(iat,ibd2)%cp(2,ilmn2)-&
&                 cprj1(iat,ibd1)%cp(2,ilmn1)*cprj2_k(iat,ibd2)%cp(1,ilmn2))
               end do
             end do
           end do
         else if(ibg1==0 .and. ibg2==0)then
           do iat=1+ia,nattyp(itypat)+ia
             do ilmn1=1,pawtab(itypat)%lmn_size
               do ilmn2=1,ilmn1
                 klmn=((ilmn1-1)*ilmn1)/2+ilmn2
                 dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj1(iat,ibd1)%cp(1,ilmn1)*cprj2(iat,ibd2)%cp(1,ilmn2)+&
&                 cprj1(iat,ibd1)%cp(2,ilmn1)*cprj2(iat,ibd2)%cp(2,ilmn2))
                 doti=doti+pawtab(itypat)%sij(klmn)*(cprj1(iat,ibd1)%cp(1,ilmn1)*cprj2(iat,ibd2)%cp(2,ilmn2)-&
&                 cprj1(iat,ibd1)%cp(2,ilmn1)*cprj2(iat,ibd2)%cp(1,ilmn2))
               end do
               do ilmn2=ilmn1+1,pawtab(itypat)%lmn_size
                 klmn=((ilmn2-1)*ilmn2)/2+ilmn1
                 dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj1(iat,ibd1)%cp(1,ilmn1)*cprj2(iat,ibd2)%cp(1,ilmn2)+&
&                 cprj1(iat,ibd1)%cp(2,ilmn1)*cprj2(iat,ibd2)%cp(2,ilmn2))
                 doti=doti+pawtab(itypat)%sij(klmn)*(cprj1(iat,ibd1)%cp(1,ilmn1)*cprj2(iat,ibd2)%cp(2,ilmn2)-&
&                 cprj1(iat,ibd1)%cp(2,ilmn1)*cprj2(iat,ibd2)%cp(1,ilmn2))
               end do
             end do
           end do
         end if
         ia=ia+nattyp(itypat)
       end do
     end if
     smn(1,ibd1,ibd2)=dotr
     smn(2,ibd1,ibd2)=doti
!    End loop over bands ibd2
     icgb2=icgb2+npw*nspinor

   end do

!  End loop over bands ibd1
   icgb1=icgb1+npw*nspinor
 end do

!Complete the matrix if hermitian
 if(hermitian==1)then
   do ibd1=1,nbd1-1
     do ibd2=ibd1+1,nbd2
       smn(1,ibd1,ibd2)= smn(1,ibd2,ibd1)
       smn(2,ibd1,ibd2)=-smn(2,ibd2,ibd1)
     end do
   end do
 end if

!DEBUG
!write(std_out,*)' smn=',smn
!ENDDEBUG

!====== Debugging section ==========
 if(.false.)then
!DEBUG
!Compute the eigenvalues of the projector S herm(S) or herm(S) S, depending on which has lowest dimension.
!write(std_out,*)' dotprod_set_cgcprj : compute the projector matrix '
   nbd=min(nbd1,nbd2)
   ABI_ALLOCATE(proj,(2,nbd,nbd))
   proj(:,:,:)=zero
   if(nbd1<=nbd2)then
     do ibd1=1,nbd1
       do ibd2=1,nbd2
         proj(1,:,ibd1)=proj(1,:,ibd1)+smn(1,:,ibd2)*smn(1,ibd1,ibd2)+smn(2,:,ibd2)*smn(2,ibd1,ibd2)
         proj(2,:,ibd1)=proj(2,:,ibd1)-smn(1,:,ibd2)*smn(2,ibd1,ibd2)+smn(2,:,ibd2)*smn(1,ibd1,ibd2)
       end do
     end do
   else
     do ibd2=1,nbd2
       do ibd1=1,nbd1
         proj(1,:,ibd2)=proj(1,:,ibd2)+smn(1,ibd1,:)*smn(1,ibd1,ibd2)+smn(2,ibd1,:)*smn(2,ibd1,ibd2)
         proj(2,:,ibd2)=proj(2,:,ibd2)+smn(1,ibd1,:)*smn(2,ibd1,ibd2)-smn(2,ibd1,:)*smn(1,ibd1,ibd2)
       end do
     end do
   end if

!write(std_out,*)' proj=',proj

!write(std_out,*)' dotprod_set_cgcprj : compute the eigenvalues of the projector '
   ABI_ALLOCATE(matrx,(2,(nbd*(nbd+1))/2))
   ii=1
   do i2=1,nbd
     do i1=1,i2
       matrx(1,ii)=proj(1,i1,i2)
       matrx(2,ii)=proj(2,i1,i2)
       ii=ii+1
     end do
   end do

   ABI_ALLOCATE(zhpev1,(2,2*nbd-1))
   ABI_ALLOCATE(zhpev2,(3*nbd-2))
   ABI_ALLOCATE(eigval,(nbd))
   ABI_ALLOCATE(eigvec,(2,nbd,nbd))

   call ZHPEV ('V','U',nbd,matrx,eigval,eigvec,nbd,zhpev1,zhpev2,ier)

!write(std_out,*)' eigval=',eigval

   ABI_DEALLOCATE(matrx)
   ABI_DEALLOCATE(zhpev1)
   ABI_DEALLOCATE(zhpev2)
   ABI_DEALLOCATE(eigval)
   ABI_DEALLOCATE(eigvec)

   ABI_DEALLOCATE(proj)
!stop
!ENDDEBUG
 end if
!====== End of debugging section ==========

 ABI_DEALLOCATE(cwavef1)
 ABI_DEALLOCATE(cwavef2)
 if(usepaw==1) then
   if(ibg1/=0)then
     call pawcprj_free(cprj1_k)
   end if
   if(ibg2/=0)then
     call pawcprj_free(cprj2_k)
   end if
   ABI_DATATYPE_DEALLOCATE(cprj1_k)
   ABI_DATATYPE_DEALLOCATE(cprj2_k)
 end if

end subroutine dotprod_set_cgcprj
!!***

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
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cg_set(2,mcg,mset)= plane wave wavefunction coefficients for several sets of wavefunctions (all k points and spins)
!!  cprj_set(natom,mcprj,mset) <type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk>
!!    with NL projectors in the different sets
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
!!  smn(2,nset1,nset2)=matrix of sum of diagonal scalar products between the first set
!!    of wavefunctions and the second set of wavefunctions
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_scfcv_core
!!
!! CHILDREN
!!      dotprod_set_cgcprj,lincom_cgcprj,zpotrf,ztrsm
!!
!! SOURCE

subroutine dotprodm_sumdiag_cgcprj(atindx1,cg_set,cprj_set,dimcprj,&
& ibg,icg,ikpt,isppol,istwf,mband,mcg,mcprj,mkmem,&
& mpi_enreg,mset,natom,nattyp,nbd,npw,nset1,nset2,nspinor,nsppol,ntypat,&
& shift_set1,shift_set2,pawtab,smn,usepaw)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ibg,icg,ikpt,isppol,istwf
 integer, intent(in) :: mband,mcg,mcprj,mkmem,mset
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
 integer :: ia,iat,itypat,ibd,icgb,ig,iorder
 integer :: ilmn1,ilmn2,ind_set1,ind_set2,iset1,iset2,klmn
 real(dp) :: dotr,doti
!arrays
 real(dp),allocatable :: cwavef1(:,:),cwavef2(:,:)
 type(pawcprj_type),allocatable :: cprj1_k(:,:),cprj2_k(:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' dotprodm_sumdiag_cgcprj : enter '
!call flush(std_out)
!ENDDEBUG

 ABI_ALLOCATE(cwavef1,(2,npw*nspinor))
 ABI_ALLOCATE(cwavef2,(2,npw*nspinor))
 if(usepaw==1) then
   ABI_DATATYPE_ALLOCATE(cprj1_k,(natom,nspinor*nbd))
   ABI_DATATYPE_ALLOCATE(cprj2_k,(natom,nspinor*nbd))
   iorder=0 ! There is no change of ordering in the copy of wavefunctions
   call pawcprj_alloc(cprj1_k,cprj_set(1,1,1)%ncpgr,dimcprj)
 end if

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
       call pawcprj_get(atindx1,cprj1_k,cprj_set(:,:,ind_set1),natom,1,ibg,ikpt,iorder,isppol,mband,&
&       mkmem,natom,nbd,nbd,nspinor,nsppol,0,&
&       mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
     end if

     do iset2=1,nset2

       ind_set2=iset2+shift_set2
       if(ind_set2<ind_set1 .and. ind_set2>shift_set1)then
         continue ! These matrix elements have already been computed, the smn matrix will be completed later.

       else if(ind_set1/=ind_set2)then

!        Extract wavefunction information
         do ig=1,npw*nspinor
           cwavef2(1,ig)=cg_set(1,ig+icgb,ind_set2)
           cwavef2(2,ig)=cg_set(2,ig+icgb,ind_set2)
         end do

         if(usepaw==1) then
           call pawcprj_get(atindx1,cprj2_k,cprj_set(:,:,ind_set2),natom,1,ibg,ikpt,iorder,isppol,mband,&
&           mkmem,natom,nbd,nbd,nspinor,nsppol,0,&
&           mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         end if

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
         end if ! usepaw

!      if(.false.)then
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
         end if ! usepaw
         doti=zero

       end if ! Compare ind_set1 and ind_set2

       smn(1,iset1,iset2)=smn(1,iset1,iset2)+dotr
       smn(2,iset1,iset2)=smn(2,iset1,iset2)+doti

     end do ! iset2
   end do ! iset1

!  End loop over bands ibd
   icgb=icgb+npw*nspinor
 end do

!Complete the matrix, using its hermitian property.
 do iset1=1,nset1
   ind_set1=iset1+shift_set1
   do iset2=1,nset2
     ind_set2=iset2+shift_set2
     if(ind_set2<ind_set1 .and. ind_set2>shift_set1)then
       smn(1,iset1,iset2)= smn(1,iset2,iset1)
       smn(2,iset1,iset2)=-smn(2,iset2,iset1)
     end if
   end do
 end do

 ABI_DEALLOCATE(cwavef1)
 ABI_DEALLOCATE(cwavef2)
 if(usepaw==1) then
   call pawcprj_free(cprj1_k)
   call pawcprj_free(cprj2_k)
   ABI_DATATYPE_DEALLOCATE(cprj1_k)
   ABI_DATATYPE_DEALLOCATE(cprj2_k)
 end if

end subroutine dotprodm_sumdiag_cgcprj
!!***

!!****f* ABINIT/lincom_cgcprj
!!
!! NAME
!! lincom_cgcprj
!!
!! FUNCTION
!! For one k point and spin, compute a set (size nband_out) of linear combinations of nband_in wavefunctions,
!! that are known in the cg+cprj representation :
!! cgout_n(:,:) <--- Sum_m [ cg_m(:,:) . alpha_mn ]
!! cprjout_n(:,:) <--- Sum_m [ cprj_m(:,:) . alpha_mn ]
!! If nband_out is smaller or equal to nband_in, the result might be in-place
!! output in cg instead of cgout, and in cprj instead of cprjout).
!! Otherwise, it is contained in the optional cgout+cprjout pair.
!!
!! In the present status, the cg and cgout relates to all the k points and spins, and rely on the icg index,
!! while it is assumed that cprj and cprjout refer to the specific k point and spin.
!! This is not coherent.
!! THIS MIGHT BE CHANGED IN THE FUTURE !
!!
!! This implementation is NOT band-parallelized
!! Also, it is far of being optimal at the level of linear algebra, and involves extra copying
!! that are detrimental for performance...
!!
!! INPUTS
!!  alpha_mn(2,nband_in,nband_out)=complex matrix of coefficients of the linear combinations to be computed
!!  dimcprj(natom)=number of lmn components in the <p_{lmn}^i|\psi> for the i-th atom
!!  icg=shift in cg array to locate current k-point and spinpol (for input, and possibly for in-place output)
!!  inplace= if 0, output in cgout and cprjout ; if 1, output in cg and cprj
!!  mcg=second dimension of cg array (mpw*nspinor*mband*mkmem*nsppol)
!!  mcprj=second dimension of cprj array
!!  natom=number of atoms
!!  nband_in=number of bands, size of the input set of wavefunctions
!!  nband_out=number of bands, size of the output set of wavefunctions (should be equal to nband_in if inplace==1)
!!  npw=number of planewaves in basis at this k point
!!  nspinor=number of spinor components
!!  usepaw=1 if PAW is activated
!!  [icgout= shift in cgout array to locate current k-point and spinpol (for output)]
!!  [mcgout=second dimension of cgout array (mpw*nspinor*mband*mkmem*nsppol)]
!!  [mcprjout=second dimension of cprjout array]
!!
!! OUTPUT
!!  [cgout(2,mcgout)= plane wave wavefunction coefficients for the set of output wavefunctions]
!!  [cprjout(natom,mcprjout) <type(pawcprj_type)>= projected output wave functions <Proj_i|Cnk> with NL projectors]
!!
!! SIDE EFFECTS
!!  (this quantities are input, and possibly updated output when inplace==1)
!!  cg(2,mcg)= plane wave wavefunction coefficients for the set of input wavefunctions (all k points and spinpol)
!!  cprj(natom,mcprj) <type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk> with NL projectors
!!
!! PARENTS
!!      m_cgcprj,m_extraprho,m_scfcv_core
!!
!! CHILDREN
!!      dotprod_set_cgcprj,lincom_cgcprj,zpotrf,ztrsm
!!
!! SOURCE

 subroutine lincom_cgcprj(alpha_mn,cg,cprj,dimcprj,&
& icg,inplace,mcg,mcprj,natom,nband_in,nband_out,npw,nspinor,usepaw, &
& cgout,cprjout,icgout) ! optional args

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: icg,inplace,mcg,mcprj
 integer, intent(in) :: natom,nband_in,nband_out,npw,nspinor,usepaw
 integer, intent(in),optional :: icgout
!arrays
 integer, intent(in) :: dimcprj(natom)
 real(dp), intent(inout) :: cg(2,mcg)
 real(dp), intent(in) :: alpha_mn(2,nband_in,nband_out)
 real(dp), intent(out),optional :: cgout(:,:)
 type(pawcprj_type),intent(inout) :: cprj(natom,mcprj)
 type(pawcprj_type),intent(out),optional :: cprjout(:,:)

!Local variables-------------------------------
!scalars
 integer :: iband_in,iband_out,ii
!arrays
 real(dp),allocatable :: al(:,:),cgout_(:,:)
 type(pawcprj_type),allocatable :: cprjout_(:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' lincom_cgcprj : enter '
!write(std_out,*)' lincom_cgcprj : npw, nspinor=',npw,nspinor
!write(std_out,*)' lincom_cgcprj : icgout=',icgout
!ENDDEBUG

 if(inplace==0)then
   if(.not.present(cgout))then
     ABI_ERROR(' inplace==0 while .not.present(cgout) is not permitted ')
   end if
   if(usepaw==1) then
     if(.not.present(cprjout))then
       ABI_ERROR(' inplace==0 and usepaw==1 while .not.present(cprjout) is not permitted ')
     end if
   end if
 end if

!Take care of the plane wave part
 ABI_ALLOCATE(cgout_,(2,npw*nspinor*nband_out))

 call zgemm('N','N',npw*nspinor,nband_out,nband_in,dcmplx(1._dp), &
& cg(:,icg+1:icg+npw*nspinor*nband_in),npw*nspinor, &
& alpha_mn,nband_in,dcmplx(0._dp),cgout_,npw*nspinor)

 if(inplace==1)then
   cg(:,icg+1:icg+npw*nspinor*nband_out)=cgout_
 else
   cgout(:,icgout+1:icgout+npw*nspinor*nband_out)=cgout_
 end if
 ABI_DEALLOCATE(cgout_)

!Take care of the cprj part
 if(usepaw==1) then

   ABI_DATATYPE_ALLOCATE(cprjout_,(natom,nspinor*nband_out))
   call pawcprj_alloc(cprjout_,cprj(1,1)%ncpgr,dimcprj)
   ABI_ALLOCATE(al,(2,nband_in))
   do iband_out=1,nband_out
     ii=(iband_out-1)*nspinor
     do iband_in=1,nband_in
       al(1,iband_in)=alpha_mn(1,iband_in,iband_out)
       al(2,iband_in)=alpha_mn(2,iband_in,iband_out)
     end do
     call pawcprj_lincom(al,cprj,cprjout_(:,ii+1:ii+nspinor),nband_in)
   end do
   ABI_DEALLOCATE(al)

   if(inplace==1)then
     cprj=cprjout_
   else
     cprjout=cprjout_
   end if
   call pawcprj_free(cprjout_)
   ABI_DATATYPE_DEALLOCATE(cprjout_)

 end if

end subroutine lincom_cgcprj
!!***

!!****m* ABINIT/cgcprj_cholesky
!! NAME
!!  cgcprj_cholesky
!!
!! FUNCTION
!! Cholesky orthonormalization of the vectors stored in cg+cprj mode.
!!
!! This implementation is NOT band-parallelized
!! Also, it is far of being optimal at the level of linear algebra
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dimcprj(natom)=number of lmn components in the <p_{lmn}^i|\psi> for the i-th atom
!!  icg=shift in cg array to locate current k-point and spinpol
!!  ikpt=current k point index
!!  isppol=current spin polarization index
!!  istwf=input option parameter that describes the storage of wfs
!!  mcg=second dimension of cg array (mpw*nspinor*mband*mkmem*nsppol)
!!  mcprj=second dimension of cprj_k array
!!  mkmem=number of k points which can fit in memory
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nband=number of bands
!!  npw=number of planewaves in basis at this k point
!!  nspinor=number of spinor components
!!  nsppol=number of spin polarizations
!!  ntypat=number of types of atoms
!!  pawtab(dtset%ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  usepaw=1 if PAW is activated
!!
!! SIDE EFFECTS
!!  cg(2,mcg)= plane wave wavefunction coefficients for the set of input wavefunctions (all k points and spinpol)
!!  cprj_k(natom,mcprj) <type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk> with NL projectors for the specific k point and spinpol
!!
!! PARENTS
!!      m_extraprho,m_scfcv_core
!!
!! CHILDREN
!!      dotprod_set_cgcprj,lincom_cgcprj,zpotrf,ztrsm
!!
!! SOURCE

 subroutine cgcprj_cholesky(atindx1,cg,cprj_k,dimcprj,icg,ikpt,isppol,istwf,mcg,mcprj,mkmem,&
&  mpi_enreg,natom,nattyp,nband,npw,nspinor,nsppol,ntypat,pawtab,usepaw)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,ikpt,isppol,istwf,mcg,mcprj,mkmem
 integer,intent(in) :: natom,nband,npw,nspinor,nsppol,ntypat,usepaw
!arrays
 integer, intent(in) :: atindx1(natom),dimcprj(natom),nattyp(ntypat)
 real(dp), intent(inout) :: cg(2,mcg)
 type(pawcprj_type),intent(inout) :: cprj_k(natom,mcprj)
 type(MPI_type),intent(in) :: mpi_enreg
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables ------------------------------
!scalars
 integer :: hermitian,ierr,ii,inplace
!arrays
 real(dp), allocatable :: dmn(:,:,:),smn(:,:,:)

! *************************************************************************

 ABI_ALLOCATE(smn,(2,nband,nband))
 ABI_ALLOCATE(dmn,(2,nband,nband))

 hermitian=1
 call dotprod_set_cgcprj(atindx1,cg,cg,cprj_k,cprj_k,dimcprj,hermitian,&
& 0,0,icg,icg,ikpt,isppol,istwf,nband,mcg,mcg,mcprj,mcprj,mkmem,&
& mpi_enreg,natom,nattyp,nband,nband,npw,nspinor,nsppol,ntypat,pawtab,smn,usepaw)

!Cholesky factorization: O = U^H U with U upper triangle matrix.
 call ZPOTRF('U',nband,smn,nband,ierr)

!Solve X U = 1.
 dmn=zero
 do ii=1,nband
   dmn(1,ii,ii)=one
 end do
 call ZTRSM('Right','Upper','Normal','Normal',nband,nband,cone,smn,nband,dmn,nband)

 inplace=1
!This call does not take into account the fact that X=dmn is an upper triangular matrix...
!The number of operations might be divided by two.
 call lincom_cgcprj(dmn,cg,cprj_k,dimcprj,&
& icg,inplace,mcg,mcprj,natom,nband,nband,npw,nspinor,usepaw)

 ABI_DEALLOCATE(smn)
 ABI_DEALLOCATE(dmn)

end subroutine cgcprj_cholesky
!!***

end module m_cgcprj
!!***
