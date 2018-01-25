!{\src2tex{textfont=tt}}
!!****f* ABINIT/dotprod_set_cgcprj
!!
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
!! COPYRIGHT
!! Copyright (C) 2017-2018 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cg1(2,mcg1)= plane wave wavefunction coefficients for the first set of wavefunctions (all k points and spinpol)
!!  cg2(2,mcg2)= plane wave wavefunction coefficients for the second set of wavefunctions (all k points and spinpol)
!!  cprj1(natom,mcprj1) <type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk> with NL projectors in the first set,
!!  cprj2(natom,mcprj2) <type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk> with NL projectors in the second set,
!!  dimcprj(natom)=number of lmn components in the <p_{lmn}^i|\psi> for the i-th atom
!!  hermitian= if 1, consider that the Smn matrix is hermitian, and do not recompute already computed matrix elements.
!!  ibg1=shift in cprj1 array to locate current k-point ! Might be 0, in which case cprj1 is not copied internally, which saves some time/space
!!  ibg2=shift in cprj2 array to locate current k-point ! Might be 0, in which case cprj2 is not copied internally, which saves some time/space
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
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dotprod_set_cgcprj(atindx1,cg1,cg2,cprj1,cprj2,dimcprj,hermitian,&
& ibg1,ibg2,icg1,icg2,ikpt,isppol,istwf,mband,mcg1,mcg2,mcprj1,mcprj2,mkmem,&
& mpi_enreg,natom,nattyp,nbd1,nbd2,npw,nspinor,nsppol,ntypat,pawtab,smn,usepaw)

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
#define ABI_FUNC 'dotprod_set_cgcprj'
!End of the abilint section

 implicit none

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
 integer :: ia,iat,itypat,ibd1,ibd2,icgb1,icgb2,ier,ig,ii,i1,i2
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
     MSG_ERROR(' With hermitian==1, nb1 and nb2 must be equal ')
   endif
 endif

 ABI_ALLOCATE(cwavef1,(2,npw*nspinor))
 ABI_ALLOCATE(cwavef2,(2,npw*nspinor))
 if(usepaw==1) then
   ABI_DATATYPE_ALLOCATE(cprj1_k,(natom,nspinor*nbd1))
   ABI_DATATYPE_ALLOCATE(cprj2_k,(natom,nspinor*nbd2))
 endif
 if(usepaw==1 .and. ibg1/=0) then
   call pawcprj_alloc(cprj1_k,cprj1(1,1)%ncpgr,dimcprj)
 endif
 if(usepaw==1 .and. ibg2/=0) then
   call pawcprj_alloc(cprj2_k,cprj1(1,1)%ncpgr,dimcprj)
 endif

 icgb1=icg1
 do ibd1=1,nbd1

!  Extract wavefunction information
   do ig=1,npw*nspinor
     cwavef1(1,ig)=cg1(1,ig+icgb1)
     cwavef1(2,ig)=cg1(2,ig+icgb1)
   end do
   if(usepaw==1 .and. ibg1/=0) then
     call pawcprj_get(atindx1,cprj1_k,cprj1,natom,1,ibg1,ikpt,1,isppol,mband,&
&         mkmem,natom,nbd1,nbd1,nspinor,nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
   endif

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
       call pawcprj_get(atindx1,cprj2_k,cprj2,natom,1,ibg2,ikpt,1,isppol,mband,&
&         mkmem,natom,nbd2,nbd2,nspinor,nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
     endif

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
         endif
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
     enddo
   enddo
 endif

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
     enddo
   enddo
 else
   do ibd2=1,nbd2
     do ibd1=1,nbd1
       proj(1,:,ibd2)=proj(1,:,ibd2)+smn(1,ibd1,:)*smn(1,ibd1,ibd2)+smn(2,ibd1,:)*smn(2,ibd1,ibd2)
       proj(2,:,ibd2)=proj(2,:,ibd2)+smn(1,ibd1,:)*smn(2,ibd1,ibd2)-smn(2,ibd1,:)*smn(1,ibd1,ibd2)
     enddo
   enddo
 endif

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
 endif
!====== End of debugging section ==========

 ABI_DEALLOCATE(cwavef1)
 ABI_DEALLOCATE(cwavef2)
 if(usepaw==1) then
   if(ibg1/=0)then
     call pawcprj_free(cprj1_k)
   endif
   if(ibg2/=0)then
     call pawcprj_free(cprj2_k)
   endif
   ABI_DATATYPE_DEALLOCATE(cprj1_k)
   ABI_DATATYPE_DEALLOCATE(cprj2_k)
 end if

end subroutine dotprod_set_cgcprj
!!***
