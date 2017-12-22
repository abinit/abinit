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
!!  cg1(2,mcg1)= plane wave wavefunction coefficients for the first set of wavefunctions (all k points and spinpol)
!!  cg2(2,mcg2)= plane wave wavefunction coefficients for the second set of wavefunctions (all k points and spinpol)
!!  cprj1(natom,mcprj1) <type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk> with NL projectors in the first set,
!!  cprj2(natom,mcprj2) <type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk> with NL projectors in the second set,
!!  dimcprj(natom)=number of lmn components in the <p_{lmn}^i|\psi> for the i-th atom
!!  ibg1=shift in cprj1 array to locate current k-point
!!  ibg2=shift in cprj2 array to locate current k-point
!!  icg1=shift in cg1 array to locate current k-point
!!  icg2=shift in cg2 array to locate current k-point
!!  ikpt=current k point index
!!  isppol=current spin polarization index
!!  istwf=input option parameter that describes the storage of wfs
!!  mcg1=second dimension of cg1 array (mpw*nspinor*mband1*mkmem*nsppol)
!!  mcg2=second dimension of cg2 array (mpw*nspinor*mband2*mkmem*nsppol)
!!  mcprj1=second dimension of cprj1 array 
!!  mcprj2=second dimension of cprj2 array 
!!  mpi_enreg=information about MPI parallelization
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nbd1=number of bands for the first set of wavefunctions
!!  nbd2=number of bands for the second set of wavefunctions
!!  npw=number of planewaves in basis at this k point
!!  nspinor=number of spinor components
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

subroutine dotprod_set_cgcprj(atindx1,cg1,cg2,cprj1,cprj2,dimcprj,&
& ibg1,ibg2,icg1,icg2,ikpt,isppol,istwf,mband,mcg1,mcg2,mcprj1,mcprj2,mkmem,&
& mpi_enreg,natom,nattyp,nbd1,nbd2,npw,nspinor,nsppol,ntypat,pawtab,smn,usepaw)

 use defs_basis
 use defs_abitypes
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
 integer, intent(in) :: ibg1,ibg2,icg1,icg2,ikpt,isppol,istwf
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
 integer :: ia,iat,itypat,ibd1,ibd2,icgb1,icgb2,ig
 integer :: ilmn1,ilmn2,klmn
 real(dp) :: dotr,doti
!arrays
 real(dp),allocatable :: cwavef1(:,:),cwavef2(:,:)
 type(pawcprj_type),allocatable :: cprj1_k(:,:),cprj2_k(:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' dotprod_set_cgcprj : enter '
!ENDDEBUG

 ABI_ALLOCATE(cwavef1,(2,npw*nspinor))
 ABI_ALLOCATE(cwavef2,(2,npw*nspinor))
 if(usepaw==1) then
   ABI_DATATYPE_ALLOCATE(cprj1_k,(natom,nspinor*nbd1))
   ABI_DATATYPE_ALLOCATE(cprj2_k,(natom,nspinor*nbd2))
 endif

 icgb1=icg1
 do ibd1=1,nbd1

!  Extract wavefunction information
   do ig=1,npw*nspinor
     cwavef1(1,ig)=cg1(1,ig+icgb1)
     cwavef1(2,ig)=cg1(2,ig+icgb1)
   end do
   if(usepaw==1) then
     call pawcprj_alloc(cprj1_k,cprj1(1,1)%ncpgr,dimcprj)
     call pawcprj_get(atindx1,cprj1_k,cprj1,natom,1,ibg1,ikpt,1,isppol,mband,&
&         mkmem,natom,nbd1,nbd1,nspinor,nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
   endif

   icgb2=icg2
   do ibd2=1,nbd2

!    XG171222 Note that this copy step, being inside the ibd1 loop, is quite detrimental.
!    It might be reduced by copying several cwavef2, and use a ZGEMM type of approach.

!    Extract wavefunction information
     do ig=1,npw*nspinor
       cwavef2(1,ig)=cg2(1,ig+icgb2)
       cwavef2(2,ig)=cg2(2,ig+icgb2)
     end do

     if(usepaw==1) then
       call pawcprj_alloc(cprj2_k,cprj2(1,1)%ncpgr,dimcprj)
       call pawcprj_get(atindx1,cprj2_k,cprj2,natom,1,ibg2,ikpt,1,isppol,mband,&
&         mkmem,natom,nbd2,nbd2,nspinor,nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
     endif

!    Calculate Smn=<cg1|cg2>
     call dotprod_g(dotr,doti,istwf,npw*nspinor,2,cwavef1,cwavef2,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

     if(usepaw==1) then
       ia =0
       do itypat=1,ntypat
         do iat=1+ia,nattyp(itypat)+ia
           do ilmn1=1,pawtab(itypat)%lmn_size
             do ilmn2=1,ilmn1
               klmn=((ilmn1-1)*ilmn1)/2+ilmn2
               dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj1(iat,ibd1)%cp(1,ilmn1)*cprj2(iat,ibd2)%cp(1,ilmn2)+&
&               cprj1(iat,ibd1)%cp(2,ilmn1)*cprj2(iat,ibd2)%cp(2,ilmn2))
               doti=doti+pawtab(itypat)%sij(klmn)*(cprj1(iat,ibd1)%cp(1,ilmn1)*cprj2(iat,ibd2)%cp(2,ilmn2)-&
&               cprj1(iat,ibd1)%cp(2,ilmn1)*cprj2(iat,ibd2)%cp(1,ilmn2))
             end do
             do ilmn2=ilmn1+1,pawtab(itypat)%lmn_size
               klmn=((ilmn2-1)*ilmn2)/2+ilmn1
               dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj1(iat,ibd1)%cp(1,ilmn1)*cprj2(iat,ibd2)%cp(1,ilmn2)+&
&               cprj1(iat,ibd1)%cp(2,ilmn1)*cprj2(iat,ibd2)%cp(2,ilmn2))
               doti=doti+pawtab(itypat)%sij(klmn)*(cprj1(iat,ibd1)%cp(1,ilmn1)*cprj2(iat,ibd2)%cp(2,ilmn2)-&
&               cprj1(iat,ibd1)%cp(2,ilmn1)*cprj2(iat,ibd2)%cp(1,ilmn2))
             end do
           end do
         end do
         ia=ia+nattyp(itypat)
       end do
     end if
     smn(1,ibd1,ibd2)=dotr
     smn(2,ibd1,ibd2)=doti
!    End loop over bands ibd1
     icgb1=icgb1+npw*nspinor

   end do

!  End loop over bands ibd2
   icgb2=icgb2+npw*nspinor
 end do

 ABI_DEALLOCATE(cwavef1)
 ABI_DEALLOCATE(cwavef2)
 if(usepaw==1) then
   call pawcprj_free(cprj1_k)
   ABI_DATATYPE_DEALLOCATE(cprj1_k)
   call pawcprj_free(cprj2_k)
   ABI_DATATYPE_DEALLOCATE(cprj2_k)
 end if

end subroutine dotprod_set_cgcprj
!!***
