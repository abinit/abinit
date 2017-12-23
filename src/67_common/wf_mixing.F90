!{\src2tex{textfont=tt}}
!!****f* ABINIT/wf_mixing
!!
!! NAME
!! wf_mixing
!!
!! FUNCTION
!! Mixing of wavefunctions in the outer loop of a double loop SCF approach.
!! Different algorithms are implemented, depending on the value of wfmixalg.
!!
!! COPYRIGHT
!! Copyright (C) 2017 ABINIT group (XG,MT,FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx1(dtset%natom)=index table for atoms, inverse of atindx
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  istep=number of call the routine (usually the outer loop in the SCF double loop)
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of cprj array 
!!  mpi_enreg=information about MPI parallelization
!!  nattyp(dtset%ntypat)=number of atoms of each type in cell.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  pawtab(dtset%ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!
!! SIDE EFFECTS
!!  cg(2,mcg)= plane wave wavefunction coefficient
!!                          Value from previous SCF cycle is input and stored in some form
!!                          Extrapolated value is output
!!  cprj(natom,mcprj) <type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk> with NL projectors
!!                          Value from previous SCF cycle is input and stored in some form
!!                          Extrapolated value is output
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
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

subroutine wf_mixing(atindx1,cg,cprj,dtset,istep,mcg,mcprj,mpi_enreg,&
& nattyp,npwarr,pawtab,scf_history)

 use defs_basis
 use defs_abitypes
 use m_scf_history
 use m_xmpi
 use m_profiling_abi
 use m_errors
 use m_cgtools

 use m_pawtab, only : pawtab_type
 use m_pawcprj, only : pawcprj_type, pawcprj_alloc, pawcprj_copy, pawcprj_get, pawcprj_lincom, &
&                      pawcprj_free, pawcprj_zaxpby, pawcprj_put, pawcprj_getdim

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wf_mixing'
 use interfaces_32_util
 use interfaces_66_wfs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep,mcg,mcprj
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(scf_history_type),intent(inout) :: scf_history
!arrays
 integer,intent(in) :: atindx1(dtset%natom),nattyp(dtset%ntypat)
 integer,intent(in) :: npwarr(dtset%nkpt)
 real(dp), intent(inout) :: cg(2,mcg)
 type(pawcprj_type),intent(inout) :: cprj(dtset%natom,mcprj)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer :: ia,iat,iatom,iband_max,iband_max1,iband_min,iband_min1,ibd,ibg,iblockbd,iblockbd1,icg,icgb,icgb1
 integer :: ierr,ig,ii,ikpt,ilmn1,ilmn2,inc,indh,ind2
 integer :: isize,isppol,istwf_k,itypat,kk,klmn,me_distrb,my_nspinor
 integer :: nband_k,nblockbd,nprocband,npw_k,npw_nk,ntypat,ortalgo,spaceComm_band,usepaw,wfmixalg
 real(dp) :: dotr,dotr1,doti,doti1
 !character(len=500) :: message
!arrays
 real(dp) :: alpha(2),beta(2)
 integer,allocatable :: bufsize(:),bufsize_wf(:),bufdisp(:),bufdisp_wf(:)
 integer,allocatable :: ipiv(:),dimcprj(:),npw_block(:),npw_disp(:)
 real(dp),allocatable :: al(:,:),cwavef(:,:),cwavefh(:,:),cwavef_tmp(:,:)
 real(dp),allocatable :: dum(:,:)
 real(dp),allocatable :: dmn(:,:,:),dmn_debug(:,:,:),mmn(:,:,:)
 real(dp),allocatable :: smn(:,:,:),smn_(:,:,:),smn_debug(:,:,:)
 real(dp),allocatable :: work(:,:),work1(:,:)
 type(pawcprj_type),allocatable :: cprj_k(:,:),cprj_kh(:,:),cprj_k3(:,:)
!DEBUG
 real(dp),allocatable :: cg_ref(:,:)
 type(pawcprj_type),allocatable :: cprj_ref(:,:)
!ENDDEBUG

! *************************************************************************

!DEBUG
 write(std_out,*)' wf_mixing : enter '
 write(std_out,*)' istep,scf_history%alpha=',istep,scf_history%alpha
 write(std_out,*)' cg(1:2,1:2)=',cg(1:2,1:2)
 write(std_out,*)' scf_history%cg(1:2,1:2,1)=',scf_history%cg(1:2,1:2,1)
 ABI_ALLOCATE(cg_ref,(2,mcg))
 cg_ref(:,:)=cg(:,:)
 ABI_DATATYPE_ALLOCATE(cprj_ref,(dtset%natom,mcprj))
 cprj_ref(:,:)=cprj(:,:)

!ENDDEBUG

 if (istep==0) return

 ntypat=dtset%ntypat
 usepaw=dtset%usepaw
 wfmixalg=scf_history%wfmixalg

!Useful array
 ABI_ALLOCATE(dimcprj,(dtset%natom))
 if (usepaw==1) then
   call pawcprj_getdim(dimcprj,dtset%natom,nattyp,ntypat,dtset%typat,pawtab,'O')
 end if

!Index for the wavefunction stored in scf_history
 indh=1

!DEBUG
! if(istep==1)then
!   scf_history%cg(:,:,1)=cg(:,:)
! else
!   cg(:,:)=scf_history%cg(:,:,1)
! endif
! return
!ENDDEBUG

!First step
 if (istep==1 .or. (wfmixalg==2 .and. abs(scf_history%alpha-one)<tol8) ) then
   scf_history%cg(:,:,1)=cg(:,:)
   if(usepaw==1) then
     scf_history%cprj(:,:,1)=cprj(:,:)
   end if

 else
!From 2nd step

!  Init parallelism
   me_distrb=mpi_enreg%me_kpt
   if (mpi_enreg%paral_kgb==1.or.mpi_enreg%paralbd==1) then
     spaceComm_band=mpi_enreg%comm_band
     nprocband=mpi_enreg%nproc_band
   else
     spaceComm_band=xmpi_comm_self
     nprocband=1
   end if

!  For the moment no band-fft parallelism
   nprocband=1

!  Additional statements if band-fft parallelism
   if (nprocband>1) then
     ABI_ALLOCATE(npw_block,(nprocband))
     ABI_ALLOCATE(npw_disp,(nprocband))
     ABI_ALLOCATE(bufsize,(nprocband))
     ABI_ALLOCATE(bufdisp,(nprocband))
     ABI_ALLOCATE(bufsize_wf,(nprocband))
     ABI_ALLOCATE(bufdisp_wf,(nprocband))
   end if

   icg=0
   ibg=0

!  LOOP OVER SPINS
   my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
   do isppol=1,dtset%nsppol

!    BIG FAT k POINT LOOP
     do ikpt=1,dtset%nkpt

!      Select k point to be treated by this proc
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
       if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) cycle

       istwf_k=dtset%istwfk(ikpt)

!      Retrieve number of plane waves
       npw_k=npwarr(ikpt)
       if (nprocband>1) then
!        Special treatment for band-fft //
         call xmpi_allgather(npw_k,npw_block,spaceComm_band,ierr)
         npw_nk=sum(npw_block);npw_disp(1)=0
         do ii=2,nprocband
           npw_disp(ii)=npw_disp(ii-1)+npw_block(ii-1)
         end do
       else
         npw_nk=npw_k
       end if

!      Allocate arrays for a wave-function (or a block of WFs)
       ABI_ALLOCATE(cwavef,(2,npw_nk*my_nspinor))
       ABI_ALLOCATE(cwavefh,(2,npw_nk*my_nspinor))


       if (nprocband>1) then
         isize=2*my_nspinor;bufsize(:)=isize*npw_block(:);bufdisp(:)=isize*npw_disp(:)
         isize=2*my_nspinor*npw_k;bufsize_wf(:)=isize
         do ii=1,nprocband
           bufdisp_wf(ii)=(ii-1)*isize
         end do
       end if

!      Space biorthogonalization

!      Loop over bands or blocks of bands
       nblockbd=nband_k/nprocband
       icgb=icg

       if(usepaw==1) then
         ABI_DATATYPE_ALLOCATE( cprj_k,(dtset%natom,my_nspinor*nblockbd))
         call pawcprj_alloc(cprj_k,cprj(1,1)%ncpgr,dimcprj)
         call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,ibg,ikpt,1,isppol,dtset%mband,&
&         dtset%mkmem,dtset%natom,nblockbd,nblockbd,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         ABI_DATATYPE_ALLOCATE( cprj_kh,(dtset%natom,my_nspinor*nblockbd))
         call pawcprj_alloc(cprj_kh,scf_history%cprj(1,1,indh)%ncpgr,dimcprj)
         call pawcprj_get(atindx1,cprj_kh,scf_history%cprj(:,:,indh),dtset%natom,1,ibg,ikpt,1,isppol,&
&         dtset%mband,dtset%mkmem,dtset%natom,nblockbd,nblockbd,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
       end if  !end usepaw=1

       ABI_ALLOCATE(smn,(2,nblockbd,nblockbd))

!DEBUG
       write(std_out,*)' Compute the S matrix, whose matrix elements are scalar products.'
!ENDDEBUG

       call dotprod_set_cgcprj(atindx1,cg,scf_history%cg(:,:,indh),cprj,scf_history%cprj(:,:,indh),dimcprj,&
&        ibg,ibg,icg,icg,ikpt,isppol,istwf_k,dtset%mband,mcg,mcg,mcprj,mcprj,dtset%mkmem,&
&        mpi_enreg,dtset%natom,nattyp,nband_k,nband_k,npw_nk,my_nspinor,dtset%nsppol,ntypat,pawtab,smn,usepaw)

!DEBUG
       do iblockbd=1,nband_k
         write(std_out, '(a,i4)')' iblockbd=',iblockbd
         write(std_out, '(a,8f12.4)')' Real:',smn(1,1:nblockbd,iblockbd)
         write(std_out, '(a,8f12.4)')' Imag:',smn(2,1:nblockbd,iblockbd)
       enddo
!ENDDEBUG

       ABI_ALLOCATE(smn_debug,(2,nblockbd,nblockbd))
       smn_debug=zero

       do iblockbd=1,nblockbd
         iband_min=1+(iblockbd-1)*nprocband
         iband_max=iblockbd*nprocband

         if(xmpi_paral==1.and.mpi_enreg%paral_kgb/=1) then
           if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband_min,iband_max,isppol,me_distrb)) cycle
         end if

!        Extract wavefunction information
         if (nprocband>1) then
!          Special treatment for band-fft //
           ABI_ALLOCATE(cwavef_tmp,(2,npw_k*my_nspinor*nprocband))
           do ig=1,npw_k*my_nspinor*nprocband
             cwavef_tmp(1,ig)=cg(1,ig+icgb)
             cwavef_tmp(2,ig)=cg(2,ig+icgb)
           end do
           call xmpi_alltoallv(cwavef_tmp,bufsize_wf,bufdisp_wf,cwavef,bufsize,bufdisp,spaceComm_band,ierr)
           ABI_DEALLOCATE(cwavef_tmp)
         else
           do ig=1,npw_k*my_nspinor
             cwavef(1,ig)=cg(1,ig+icgb)
             cwavef(2,ig)=cg(2,ig+icgb)
           end do
         end if

         icgb1=icg

         do iblockbd1=1,nblockbd
           iband_min1=1+(iblockbd1-1)*nprocband
           iband_max1=iblockbd1*nprocband

           if(xmpi_paral==1.and.mpi_enreg%paral_kgb/=1) then
             if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband_min1,iband_max1,isppol,me_distrb)) cycle
           end if

!          Extract wavefunction information

           if (nprocband>1) then
!            Special treatment for band-fft //
             ABI_ALLOCATE(cwavef_tmp,(2,npw_k*my_nspinor*nprocband))
             do ig=1,npw_k*my_nspinor*nprocband
               cwavef_tmp(1,ig)=scf_history%cg(1,ig+icgb1,indh)
               cwavef_tmp(2,ig)=scf_history%cg(2,ig+icgb1,indh)
             end do
             call xmpi_alltoallv(cwavef_tmp,bufsize_wf,bufdisp_wf,cwavefh,bufsize,bufdisp,spaceComm_band,ierr)
             ABI_DEALLOCATE(cwavef_tmp)
           else
             do ig=1,npw_k*my_nspinor
               cwavefh(1,ig)=scf_history%cg(1,ig+icgb1,indh)
               cwavefh(2,ig)=scf_history%cg(2,ig+icgb1,indh)
             end do
           end if

!          Calculate Smn=<cg|S|cg_hist>
           call dotprod_g(dotr,doti,istwf_k,npw_k*my_nspinor,2,cwavef,cwavefh,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
           if(usepaw==1) then
             ia =0
             do itypat=1,ntypat
               do iat=1+ia,nattyp(itypat)+ia
                 do ilmn1=1,pawtab(itypat)%lmn_size
                   do ilmn2=1,ilmn1
                     klmn=((ilmn1-1)*ilmn1)/2+ilmn2
                     dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_kh(iat,iblockbd1)%cp(1,ilmn2)+&
&                     cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_kh(iat,iblockbd1)%cp(2,ilmn2))
                     doti=doti+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_kh(iat,iblockbd1)%cp(2,ilmn2)-&
&                     cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_kh(iat,iblockbd1)%cp(1,ilmn2))
                   end do
                   do ilmn2=ilmn1+1,pawtab(itypat)%lmn_size
                     klmn=((ilmn2-1)*ilmn2)/2+ilmn1
                     dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_kh(iat,iblockbd1)%cp(1,ilmn2)+&
&                     cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_kh(iat,iblockbd1)%cp(2,ilmn2))
                     doti=doti+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_kh(iat,iblockbd1)%cp(2,ilmn2)-&
&                     cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_kh(iat,iblockbd1)%cp(1,ilmn2))
                   end do
                 end do
               end do
               ia=ia+nattyp(itypat)
             end do
           end if
           smn_debug(1,iblockbd,iblockbd1)=dotr
           smn_debug(2,iblockbd,iblockbd1)=doti
!          End loop over bands iblockbd1
           icgb1=icgb1+npw_k*my_nspinor*nprocband

         end do

!        End loop over bands iblockbd
         icgb=icgb+npw_k*my_nspinor*nprocband
       end do

!DEBUG
       do iblockbd=1,nband_k
         write(std_out, '(a,i4)')' iblockbd=',iblockbd
         write(std_out, '(a,8f12.4)')' Real:',smn_debug(1,1:nblockbd,iblockbd)
         write(std_out, '(a,8f12.4)')' Imag:',smn_debug(2,1:nblockbd,iblockbd)
       enddo
       if(maxval(abs(smn-smn_debug))>tol8)then
         write(std_out,*)' wf_mixing : smn and smn_debug do not agree '
         stop
       endif
!      smn=smn_debug
!ENDDEBUG

!      Invert S matrix, which is NOT hermitian. 

!      Calculate M=S^-1
       ABI_ALLOCATE(mmn,(2,nband_k,nband_k))
       mmn=zero
       do kk=1,nband_k
         mmn(1,kk,kk)=one
       end do

!      Cholesky factorisation of smn=Lx(trans(L)*. On output mkl=L being a lower triangular matrix.
!      call zpotrf("L",nband_k,smn,nband_k,ierr)
!      call ztrtrs("L","N","N",nband_k,nband_k,smn,nband_k,mmn,nband_k,ierr)

       ABI_ALLOCATE(smn_,(2,nband_k,nband_k))
       ABI_ALLOCATE(ipiv,(nband_k))
!      The smn_ arrays stores a copy of the smn array, and will be destroyed by the following inverse call
       smn_=smn
       call zgesv(nband_k,nband_k,smn_,nband_k,ipiv,mmn,nband_k,ierr)
       ABI_DEALLOCATE(ipiv)
       ABI_DEALLOCATE(smn_)
!DEBUG
       if(ierr/=0)then
         write(std_out,*)' wf_mixing : the call to cgesv general inversion routine returned an error code ierr=',ierr
         stop
       endif
!ENDDEBUG

!DEBUG
!Print the M matrix
       write(std_out,*)' Print the M matrix.'
       do iblockbd=1,nblockbd
         write(std_out, '(a,i4)')' iblockbd=',iblockbd
         write(std_out, '(a,8f12.4)')' Real:',mmn(1,1:nblockbd,iblockbd)
         write(std_out, '(a,8f12.4)')' Imag:',mmn(2,1:nblockbd,iblockbd)
       end do
       write(std_out,*)' Check M * S = 1.'
       ABI_ALLOCATE(dmn,(2,nband_k,nband_k))
       dmn=zero
       do iblockbd=1,nblockbd
         do iblockbd1=1,nblockbd
           dmn(1,:,iblockbd)=dmn(1,:,iblockbd)&
&           +mmn(1,:,iblockbd1)*smn(1,iblockbd1,iblockbd)-mmn(2,:,iblockbd1)*smn(2,iblockbd1,iblockbd)
           dmn(2,iblockbd,:)=dmn(2,iblockbd,:)&
&           +mmn(1,:,iblockbd1)*smn(2,iblockbd1,iblockbd)+mmn(2,:,iblockbd1)*smn(1,iblockbd1,iblockbd)
         enddo
       enddo
       write(std_out,*)' Print the M*S matrix.'
       do iblockbd=1,nblockbd
         write(std_out, '(a,i4)')' iblockbd=',iblockbd
         write(std_out, '(a,8f12.4)')' Real:',dmn(1,1:nblockbd,iblockbd)
         write(std_out, '(a,8f12.4)')' Imag:',dmn(2,1:nblockbd,iblockbd)
       end do
       ABI_ALLOCATE(dmn_debug,(2,nband_k,nband_k))
       dmn_debug=zero
       do kk=1,nband_k
         dmn_debug(1,kk,kk)=one
       end do
       if(maxval(abs(dmn-dmn_debug))>tol8)then
         write(std_out,*)' wf_mixing : dmn and dmn_debug do not agree '
         stop
       endif
       ABI_DEALLOCATE(dmn)
       ABI_DEALLOCATE(dmn_debug)
!ENDDEBUG


!      This is the simple mixing case : the wavefunction from scf_history is biorthogonalized to cg, taken as reference
!      Wavefunction alignment (istwfk=1 ?)
       ABI_ALLOCATE(work,(2,npw_nk*my_nspinor*nblockbd))
       ABI_ALLOCATE(work1,(2,npw_nk*my_nspinor*nblockbd))
       work1(:,:)=scf_history%cg(:,icg+1:icg+my_nspinor*nblockbd*npw_nk,indh)
       call zgemm('N','N',npw_nk*my_nspinor,nband_k,nband_k,dcmplx(1._dp), &
&       work1,npw_nk*my_nspinor, &
&       mmn,nblockbd,dcmplx(0._dp),work,npw_nk*my_nspinor)
       scf_history%cg(:,1+icg:npw_nk*my_nspinor*nblockbd+icg,indh)=work(:,:)

!      If paw, must also align cprj from history
       if (usepaw==1) then
!        New version (MT):
         ABI_DATATYPE_ALLOCATE(cprj_k3,(dtset%natom,my_nspinor))
         call pawcprj_alloc(cprj_k3,cprj_kh(1,1)%ncpgr,dimcprj)
         ABI_ALLOCATE(al,(2,nblockbd))
         do iblockbd=1,nblockbd
           ii=(iblockbd-1)*my_nspinor
           do iblockbd1=1,nblockbd
             al(1,iblockbd1)=mmn(1,iblockbd,iblockbd1)
             al(2,iblockbd1)=mmn(2,iblockbd,iblockbd1)
           end do
           call pawcprj_lincom(al,cprj_kh,cprj_k3,nblockbd)
           call pawcprj_copy(cprj_k3,cprj_kh(:,ii+1:ii+my_nspinor))
         end do
         ABI_DEALLOCATE(al)
         call pawcprj_free(cprj_k3)
         ABI_DATATYPE_DEALLOCATE(cprj_k3)
       end if
       ABI_DEALLOCATE(mmn)
       ABI_DEALLOCATE(work)
       ABI_DEALLOCATE(work1)

!DEBUG
!      This is a check that now the scf_history%cg(:,:,indh) is biorthogonal to the cg
!      Calculate Smn=<cg|S|cg_hist>
       write(std_out,*)' Check that the biorthogonalized scf_history%cg is indeed biorthogonal to cg'
       call dotprod_set_cgcprj(atindx1,cg,scf_history%cg(:,:,indh),cprj,scf_history%cprj(:,:,indh),dimcprj,&
&        ibg,ibg,icg,icg,ikpt,isppol,istwf_k,dtset%mband,mcg,mcg,mcprj,mcprj,dtset%mkmem,&
&        mpi_enreg,dtset%natom,nattyp,nband_k,nband_k,npw_nk,my_nspinor,dtset%nsppol,ntypat,pawtab,smn,usepaw)
       do iblockbd=1,nband_k
         write(std_out, '(a,i4)')' iblockbd=',iblockbd
         write(std_out, '(a,8f12.4)')' Real:',smn(1,1:nblockbd,iblockbd)
         write(std_out, '(a,8f12.4)')' Imag:',smn(2,1:nblockbd,iblockbd)
       enddo
!ENDDEBUG

!      Wavefunction extrapolation, simple mixing case
       ibd=0  
       inc=npw_nk*my_nspinor
       do iblockbd=1,nblockbd
!        scf_history%alpha contains dtset%wfmix in the simple mixing case.
         cg(:,icg+1+ibd:icg+inc+ibd)=scf_history%cg(:,1+icg+ibd:icg+ibd+inc,indh)&
&          +scf_history%alpha*(cg(:,icg+1+ibd:ibd+icg+inc)-scf_history%cg(:,1+icg+ibd:icg+ibd+inc,indh))
         if(usepaw==1) then
           alpha(1)=one-scf_history%alpha;alpha(2)=zero
           beta(1)=scf_history%alpha;beta(2)=zero
           call pawcprj_zaxpby(alpha,beta,cprj_kh(:,iblockbd:iblockbd),cprj_k(:,iblockbd:iblockbd))
         endif
         ibd=ibd+inc
       end do ! end loop on iblockbd

!DEBUG
!      This is a check that the new cg is biorthogonal to the cg_ref
!      Calculate Smn=<cg|S|cg_hist>
       write(std_out,*)' Check that the new extrapolated cg is biorthogonal to cg_ref'
       call dotprod_set_cgcprj(atindx1,cg,cg_ref,cprj,cprj_ref,dimcprj,&
&        ibg,ibg,icg,icg,ikpt,isppol,istwf_k,dtset%mband,mcg,mcg,mcprj,mcprj,dtset%mkmem,&
&        mpi_enreg,dtset%natom,nattyp,nband_k,nband_k,npw_nk,my_nspinor,dtset%nsppol,ntypat,pawtab,smn,usepaw)
       do iblockbd=1,nband_k
         write(std_out, '(a,i4)')' iblockbd=',iblockbd
         write(std_out, '(a,8f12.4)')' Real:',smn(1,1:nblockbd,iblockbd)
         write(std_out, '(a,8f12.4)')' Imag:',smn(2,1:nblockbd,iblockbd)
       enddo
!ENDDEBUG


!      Back to usual orthonormalization (borrowed from vtowfk and wfconv - which is problematic for PAW).
       ortalgo=mpi_enreg%paral_kgb
!      There are dummy arguments as PAW is not implemented with gsc in the present status !! 
       ABI_ALLOCATE(dum,(2,0))
!      Warning :  for the band-fft parallelism, perhaps npw_k has to be used instead of npw_nk
       call pw_orthon(icg,0,istwf_k,mcg,0,npw_nk*my_nspinor,nband_k,ortalgo,dum,usepaw,cg,&
&        mpi_enreg%me_g0,mpi_enreg%comm_bandspinorfft)
       ABI_DEALLOCATE(dum)

!DEBUG
!      This is a check that the new cg is orthonotmalized
!      Calculate Smn=<cg|S|cg>
       write(std_out,*)' Check that the final extrapolated cg is orthonormalized '
       call dotprod_set_cgcprj(atindx1,cg,cg,cprj,cprj,dimcprj,&
&        ibg,ibg,icg,icg,ikpt,isppol,istwf_k,dtset%mband,mcg,mcg,mcprj,mcprj,dtset%mkmem,&
&        mpi_enreg,dtset%natom,nattyp,nband_k,nband_k,npw_nk,my_nspinor,dtset%nsppol,ntypat,pawtab,smn,usepaw)
       do iblockbd=1,nband_k
         write(std_out, '(a,i4)')' iblockbd=',iblockbd
         write(std_out, '(a,8f12.4)')' Real:',smn(1,1:nblockbd,iblockbd)
         write(std_out, '(a,8f12.4)')' Imag:',smn(2,1:nblockbd,iblockbd)
       enddo
!      stop
!ENDDEBUG


!      Store the newly extrapolated wavefunctions, orthonormalized, in scf_history
       ibd=0
       inc=npw_nk*my_nspinor
       do iblockbd=1,nblockbd
         scf_history%cg(:,1+icg+ibd:icg+ibd+inc,indh)=cg(:,icg+1+ibd:ibd+icg+inc)
         if(usepaw==1) then
           call pawcprj_put(atindx1,cprj_k,scf_history%cprj(:,:,indh),dtset%natom,1,ibg,ikpt,1,isppol,&
&           dtset%mband,dtset%mkmem,dtset%natom,nblockbd,nblockbd,dimcprj,my_nspinor,dtset%nsppol,0,&
&           mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
         end if
         ibd=ibd+inc
       end do ! end loop on iblockbd

       ABI_DEALLOCATE(cwavef)
       ABI_DEALLOCATE(cwavefh)
       ABI_DEALLOCATE(smn)
       ABI_DEALLOCATE(smn_debug)
       if(usepaw==1) then
         call pawcprj_free(cprj_k)
         ABI_DATATYPE_DEALLOCATE(cprj_k)
         call pawcprj_free(cprj_kh)
         ABI_DATATYPE_DEALLOCATE(cprj_kh)
       end if

       ibg=ibg+my_nspinor*nband_k
       icg=icg+my_nspinor*nband_k*npw_k

!      End big k point loop
     end do
!    End loop over spins
   end do

   if (nprocband>1) then
     ABI_DEALLOCATE(npw_block)
     ABI_DEALLOCATE(npw_disp)
     ABI_DEALLOCATE(bufsize)
     ABI_DEALLOCATE(bufdisp)
     ABI_DEALLOCATE(bufsize_wf)
     ABI_DEALLOCATE(bufdisp_wf)
   end if

 end if ! istep>=2

!DEBUG
 write(std_out,*)' wf_mixing : exit '
 write(std_out,*)' cg(1:2,1:2)=',cg(1:2,1:2)
 write(std_out,*)' scf_history%cg(1:2,1:2,1)=',scf_history%cg(1:2,1:2,1)
 ABI_DEALLOCATE(cg_ref)
 ABI_DATATYPE_DEALLOCATE(cprj_ref)
!ENDDEBUG

 ABI_DEALLOCATE(dimcprj)

end subroutine wf_mixing
!!***
