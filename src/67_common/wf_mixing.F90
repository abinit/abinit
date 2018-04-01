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
!! Copyright (C) 2017-2018 ABINIT group (XG,MT,FJ)
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
!!  scf_history_wf <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      cgcprj_cholesky,dotprod_set_cgcprj,dotprodm_sumdiag_cgcprj
!!      lincom_cgcprj,pawcprj_alloc,pawcprj_axpby,pawcprj_free,pawcprj_get
!!      pawcprj_getdim,pawcprj_lincom,pawcprj_put,timab,xmpi_sum,zgesv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wf_mixing(atindx1,cg,cprj,dtset,istep,mcg,mcprj,mpi_enreg,&
& nattyp,npwarr,pawtab,scf_history_wf)

 use defs_basis
 use defs_abitypes
 use m_scf_history
 use m_xmpi
 use m_profiling_abi
 use m_errors
 use m_cgtools

 use m_time,    only : timab
 use m_pawtab,  only : pawtab_type
 use m_pawcprj, only : pawcprj_type, pawcprj_alloc, pawcprj_copy, pawcprj_get, pawcprj_lincom, &
&                      pawcprj_free, pawcprj_axpby, pawcprj_put, pawcprj_getdim

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wf_mixing'
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_66_wfs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep,mcg,mcprj
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(scf_history_type),intent(inout) :: scf_history_wf
!arrays
 integer,intent(in) :: atindx1(dtset%natom),nattyp(dtset%ntypat)
 integer,intent(in) :: npwarr(dtset%nkpt)
 real(dp), intent(inout) :: cg(2,mcg)
 type(pawcprj_type),intent(inout) :: cprj(dtset%natom,mcprj)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer :: hermitian
 integer :: ibdmix,ibdsp,ibg,ibg_hist,icg,icg_hist
 integer :: ierr,ikpt,indh,ind_biorthog,ind_biorthog_eff,ind_newwf,ind_residual,inplace
 integer :: iset2,isppol,istep_cycle,istep_new,istwf_k,kk,me_distrb,my_nspinor
 integer :: nband_k,nbdmix,npw_k,nset1,nset2,ntypat
 integer :: shift_set1,shift_set2,spaceComm_band,spare_mem,usepaw,wfmixalg
 real(dp) :: alpha,beta
 complex(dpc) :: sum_coeffs
!arrays
 integer,allocatable :: ipiv(:),dimcprj(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: al(:,:),mmn(:,:,:)
 real(dp),allocatable :: dotprod_res(:,:,:),dotprod_res_k(:,:,:),res_mn(:,:,:),smn(:,:,:)
 complex(dpc),allocatable :: coeffs(:)
 type(pawcprj_type),allocatable :: cprj_k(:,:),cprj_kh(:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)
!write(std_out,*)' wf_mixing : enter, istep= ',istep
!call flush(std_out)
!write(std_out,*)' istep,scf_history_wf%alpha=',istep,scf_history_wf%alpha
!write(std_out,*)' cg(1,1)=',cg(1,1)
!write(std_out,*)' scf_history_wf%cg(1,1,1:5)=',scf_history_wf%cg(1,1,1:5)
!ABI_ALLOCATE(cg_ref,(2,mcg))
!cg_ref(:,:)=cg(:,:)
!ABI_DATATYPE_ALLOCATE(cprj_ref,(dtset%natom,mcprj))
!cprj_ref(:,:)=cprj(:,:)
!      write(std_out,*)' scf_history_wf%dotprod_sumdiag_cgcprj_ij(:,2,2)=',&
!&       scf_history_wf%dotprod_sumdiag_cgcprj_ij(:,2,2)
!       call flush(std_out)
!ENDDEBUG

 if (istep==0) return

 ntypat=dtset%ntypat
 usepaw=dtset%usepaw
 wfmixalg=scf_history_wf%wfmixalg
 nbdmix=dtset%nbandhf
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 me_distrb=mpi_enreg%me_kpt
 spaceComm_band=xmpi_comm_self

 spare_mem=0
 if(scf_history_wf%history_size==wfmixalg-1)spare_mem=1

!scf_history_wf%alpha contains dtset%wfmix
 alpha=scf_history_wf%alpha
 beta=one-scf_history_wf%alpha

 icg=0
 icg_hist=0
 ibg=0
 ibg_hist=0

!Useful array
 ABI_ALLOCATE(dimcprj,(dtset%natom))
 if (usepaw==1) then
   call pawcprj_getdim(dimcprj,dtset%natom,nattyp,ntypat,dtset%typat,pawtab,'O')
 end if

 if(istep==1)then
   do indh=1,scf_history_wf%history_size
     call pawcprj_alloc(scf_history_wf%cprj(:,:,indh),0,dimcprj)
   end do
 end if

 ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,my_nspinor*nbdmix))
 ABI_DATATYPE_ALLOCATE(cprj_kh,(dtset%natom,my_nspinor*nbdmix))
 if(usepaw==1) then
   call pawcprj_alloc(cprj_k,0,dimcprj)
   call pawcprj_alloc(cprj_kh,0,dimcprj)
 end if
 ABI_ALLOCATE(smn,(2,nbdmix,nbdmix))
 ABI_ALLOCATE(mmn,(2,nbdmix,nbdmix))

 if(wfmixalg>2)then
   nset1=1 
   nset2=min(istep-1,wfmixalg-1) 
   ABI_ALLOCATE(dotprod_res_k,(2,1,nset2))
   ABI_ALLOCATE(dotprod_res,(2,1,nset2))
   ABI_ALLOCATE(res_mn,(2,wfmixalg-1,wfmixalg-1)) 
   dotprod_res=zero
   if(istep==1)then
     scf_history_wf%dotprod_sumdiag_cgcprj_ij=zero
   end if
 end if

!Explanation for the index for the wavefunction stored in scf_history_wf
!The reference is the cg+cprj output after the wf optimization at istep 1. 
!It comes as input to the present routine as cgcprj input at step 2, and is usually found at indh=1.

!In the simple mixing case (wfmixalg==1), the reference is never stored, because it is used "on-the-fly" to biothogonalize the
!previous input (that was stored in indh=1), then generate the next input, which is stored again in indh=1

!When the storage is not spared: 
!- the values of indh from 2 to wfmixalg store the (computed here) biorthogonalized input cgcprj, then the residual
!- the values of indh from wfmixalg+1 to 2*wfmixalg-1 store the biorthogonalized output cgcprj (coming as argument)

!First step
 if (istep==1 .or. (wfmixalg==2 .and. abs(scf_history_wf%alpha-one)<tol8) ) then

   indh=2   ! This input wavefunction is NOT the reference
   if(wfmixalg==2)indh=1 ! But this does not matter in the simple mixing case that has history_size=1

!  Simply store the wavefunctions and cprj. However, nband_k might be different from nbandhf...
!  LOOP OVER SPINS
   do isppol=1,dtset%nsppol

!    BIG FAT k POINT LOOP
     do ikpt=1,dtset%nkpt

!      Select k point to be treated by this proc
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
       if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) cycle

       npw_k=npwarr(ikpt)

       scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,indh)=cg(:,icg+1:icg+my_nspinor*npw_k*nbdmix)
       if(usepaw==1) then
!        scf_history_wf%cprj(:,ibg_hist+1:ibg_hist+my_nspinor*nbdmix,1)=cprj(:,ibg+1:ibg+my_nspinor*nbdmix)
         call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,ibg,ikpt,1,isppol,dtset%mband,&
&         dtset%mkmem,dtset%natom,nbdmix,nband_k,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         call pawcprj_put(atindx1,cprj_k,scf_history_wf%cprj(:,:,indh),dtset%natom,1,ibg_hist,ikpt,1,isppol,&
&         nbdmix,dtset%mkmem,dtset%natom,nbdmix,nbdmix,dimcprj,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
       end if

!      Update the counters
       ibg=ibg+my_nspinor*nband_k
       ibg_hist=ibg_hist+my_nspinor*nbdmix
       icg=icg+my_nspinor*nband_k*npw_k
       icg_hist=icg_hist+my_nspinor*nbdmix*npw_k

     end do
   end do

 else
!  From istep==2

!  First part of the computation : biorthogonalization, and computation of the residual (possibly, prediction of the next input in the case of simple mixing)
!  Index for the wavefunctions stored in scf_history_wf whose scalar products with the argument cgcprj will have to be computed.
   indh=1   ! This input wavefunction is the reference
   if(wfmixalg/=2 .and. istep==2)indh=2 ! except for istep=2 in the rmm-diis

   if(wfmixalg>2)then
!    istep inside the cycle defined by wfmixalg, and next index. Then, indices of the wavefunction sets.
     istep_cycle=mod((istep-2),wfmixalg-1)
     istep_new=mod((istep-1),wfmixalg-1)
     ind_biorthog=1+wfmixalg+istep_cycle
     ind_residual=2+istep_cycle
     ind_newwf=2+istep_new
     shift_set1=ind_residual-1
     shift_set2=1
   end if

!  LOOP OVER SPINS
   do isppol=1,dtset%nsppol

!    BIG FAT k POINT LOOP
     do ikpt=1,dtset%nkpt

!      Select k point to be treated by this proc
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
       if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) cycle

       istwf_k=dtset%istwfk(ikpt)
       npw_k=npwarr(ikpt)

!      Biorthogonalization

       if(usepaw==1) then
         call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,ibg,ikpt,1,isppol,dtset%mband,&
&         dtset%mkmem,dtset%natom,nbdmix,nband_k,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         call pawcprj_get(atindx1,cprj_kh,scf_history_wf%cprj(:,:,indh),dtset%natom,1,ibg_hist,ikpt,1,isppol,&
&         nbdmix,dtset%mkmem,dtset%natom,nbdmix,nbdmix,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
       end if  !end usepaw=1

       hermitian=0
       if(wfmixalg==2 .or. istep==2)then
         call dotprod_set_cgcprj(atindx1,cg,scf_history_wf%cg(:,:,indh),cprj_k,cprj_kh,dimcprj,hermitian,&
&         0,0,icg,icg_hist,ikpt,isppol,istwf_k,nbdmix,mcg,mcg,mcprj,mcprj,dtset%mkmem,&
&         mpi_enreg,dtset%natom,nattyp,nbdmix,nbdmix,npw_k,my_nspinor,dtset%nsppol,ntypat,pawtab,smn,usepaw)
       else
         call dotprod_set_cgcprj(atindx1,scf_history_wf%cg(:,:,indh),cg,cprj_kh,cprj_k,dimcprj,hermitian,&
&         0,0,icg,icg_hist,ikpt,isppol,istwf_k,nbdmix,mcg,mcg,mcprj,mcprj,dtset%mkmem,&
&         mpi_enreg,dtset%natom,nattyp,nbdmix,nbdmix,npw_k,my_nspinor,dtset%nsppol,ntypat,pawtab,smn,usepaw)
       end if

!      Invert S matrix, that is NOT hermitian. 
!      Calculate M=S^-1
       mmn=zero
       do kk=1,nbdmix
         mmn(1,kk,kk)=one
       end do

       ABI_ALLOCATE(ipiv,(nbdmix))
!      The smn is destroyed by the following inverse call
       call zgesv(nbdmix,nbdmix,smn,nbdmix,ipiv,mmn,nbdmix,ierr)
       ABI_DEALLOCATE(ipiv)
!DEBUG
       if(ierr/=0)then
         MSG_ERROR(' The call to cgesv general inversion routine failed')
       end if
!ENDDEBUG

!      The M matrix is used to compute the biorthogonalized set of wavefunctions, and to store it at the proper place
       if(wfmixalg==2 .or. istep==2)then
         inplace=1
         call lincom_cgcprj(mmn,scf_history_wf%cg(:,:,indh),cprj_kh,dimcprj,&
&         icg_hist,inplace,mcg,my_nspinor*nbdmix,dtset%natom,nbdmix,nbdmix,npw_k,my_nspinor,usepaw)
       else
         inplace=0
         call lincom_cgcprj(mmn,cg,cprj_k,dimcprj,&
&         icg,inplace,mcg,my_nspinor*nbdmix,dtset%natom,nbdmix,nbdmix,npw_k,my_nspinor,usepaw,&
&         cgout=scf_history_wf%cg(:,:,ind_biorthog),cprjout=scf_history_wf%cprj(:,:,ind_biorthog),icgout=icg_hist)
       end if

!      The biorthogonalised set of wavefunctions is now stored at the proper place

!      Finalize this first part of the computation, depending on the algorithm and the step.

       if(wfmixalg==2)then

!        Wavefunction extrapolation, simple mixing case
!        alpha contains dtset%wfmix, beta contains one-alpha
         cg(:,icg+1:icg+my_nspinor*npw_k*nbdmix)=&
&         alpha*cg(:,icg+1:icg+my_nspinor*npw_k*nbdmix)&
&         +beta*scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,indh) 
         if(usepaw==1) then
           do ibdmix=1,nbdmix
             call pawcprj_axpby(beta,alpha,cprj_kh(:,ibdmix:ibdmix),cprj_k(:,ibdmix:ibdmix))
           end do ! end loop on ibdmix
         end if

!        Back to usual orthonormalization 
         call cgcprj_cholesky(atindx1,cg,cprj_k,dimcprj,icg,ikpt,isppol,istwf_k,mcg,my_nspinor*nband_k,dtset%mkmem,&
&         mpi_enreg,dtset%natom,nattyp,nbdmix,npw_k,my_nspinor,dtset%nsppol,ntypat,pawtab,usepaw)

!        Store the newly extrapolated wavefunctions, orthonormalized, in scf_history_wf
         scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,indh)=cg(:,icg+1:icg+my_nspinor*npw_k*nbdmix)
         if(usepaw==1) then
           do ibdmix=1,nbdmix
             call pawcprj_put(atindx1,cprj_k,scf_history_wf%cprj(:,:,indh),dtset%natom,1,ibg_hist,ikpt,1,isppol,&
&             nbdmix,dtset%mkmem,dtset%natom,nbdmix,nbdmix,dimcprj,my_nspinor,dtset%nsppol,0,&
&             mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
           end do ! end loop on ibdmix
         end if

       else  !  wfmixalg/=2
!        RMM-DIIS

         if (istep==2)then
!          Store the argument wf as the reference for all future steps, in scf_history_wf with index 1. 
           scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,1)=cg(:,icg+1:icg+my_nspinor*npw_k*nbdmix)
           if(usepaw==1) then
             do ibdmix=1,nbdmix
               call pawcprj_put(atindx1,cprj_k,scf_history_wf%cprj(:,:,1),dtset%natom,1,ibg_hist,ikpt,1,isppol,&
&               nbdmix,dtset%mkmem,dtset%natom,nbdmix,nbdmix,dimcprj,my_nspinor,dtset%nsppol,0,&
&               mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
             end do ! end loop on ibdmix
           end if
         end if

         ind_biorthog_eff=ind_biorthog
         if(istep==2)ind_biorthog_eff=1 ! The argument wf has not been stored in ind_biorthog

!        Compute the residual of the wavefunctions for this istep, 
!        that replaces the previously stored set of (biorthogonalized) input wavefunctions
         scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,ind_residual)=&
&         scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,ind_biorthog_eff)&
&         -scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,ind_residual)
         if(usepaw==1) then
           do ibdmix=1,nbdmix
             call pawcprj_axpby(one,-one,scf_history_wf%cprj(:,ibdmix:ibdmix,ind_biorthog_eff),&
&             scf_history_wf%cprj(:,ibdmix:ibdmix,ind_residual))
           end do ! end loop on ibdmix
         end if
         
!        Compute the new scalar products to fill the res_mn matrix
         call dotprodm_sumdiag_cgcprj(atindx1,scf_history_wf%cg,scf_history_wf%cprj,dimcprj,&
&         ibg,icg,ikpt,isppol,istwf_k,nbdmix,mcg,mcprj,dtset%mkmem,&
&         mpi_enreg,scf_history_wf%history_size,dtset%natom,nattyp,nbdmix,npw_k,nset1,nset2,my_nspinor,dtset%nsppol,ntypat,&
&         shift_set1,shift_set2,pawtab,dotprod_res_k,usepaw)

         dotprod_res=dotprod_res+dotprod_res_k

!        scf_history_wf for index ind_biorthog will contain the extrapolated wavefunctions (and no more the output of the SCF loop).
         scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,ind_biorthog)=&
&         scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,ind_biorthog_eff)+&
&         (alpha-one)*scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,ind_residual)
         if(usepaw==1) then
           do ibdmix=1,nbdmix
             if(ind_biorthog/=ind_biorthog_eff)then
               scf_history_wf%cprj(:,ibdmix:ibdmix,ind_biorthog)=scf_history_wf%cprj(:,ibdmix:ibdmix,ind_biorthog_eff) 
             end if
             call pawcprj_axpby((alpha-one),one,scf_history_wf%cprj(:,ibdmix:ibdmix,ind_residual),&
&             scf_history_wf%cprj(:,ibdmix:ibdmix,ind_biorthog))
           end do ! end loop on ibdmix
         end if

       end if

       ibg=ibg+my_nspinor*nband_k
       ibg_hist=ibg_hist+my_nspinor*nbdmix
       icg=icg+my_nspinor*nband_k*npw_k
       icg_hist=icg_hist+my_nspinor*nbdmix*npw_k

!      End big k point loop
     end do
!    End loop over spins
   end do

 end if ! istep>=2

 if(wfmixalg>2 .and. istep>1)then

!DEBUG
!      write(std_out,*)' '
!      write(std_out,*)' Entering the residual minimisation part '
!      write(std_out,*)' '
!      call flush(std_out)
!ENDDEBUG

   call timab(48,1,tsec)
   call xmpi_sum(dotprod_res,mpi_enreg%comm_kpt,ierr)
   call timab(48,2,tsec)

   scf_history_wf%dotprod_sumdiag_cgcprj_ij(:,1+shift_set1,1+shift_set2:nset2+shift_set2)=dotprod_res(:,1,1:nset2)
   scf_history_wf%dotprod_sumdiag_cgcprj_ij(1,1+shift_set2:nset2+shift_set2,1+shift_set1)=dotprod_res(1,1,1:nset2)
   scf_history_wf%dotprod_sumdiag_cgcprj_ij(2,1+shift_set2:nset2+shift_set2,1+shift_set1)=-dotprod_res(2,1,1:nset2)

 end if ! wfmixalg>2 and istep>1

 if(wfmixalg>2 .and. istep>2)then

!  Extract the relevant matrix R_mn
   res_mn(:,1:nset2,1:nset2)=&
&   scf_history_wf%dotprod_sumdiag_cgcprj_ij(:,1+shift_set2:nset2+shift_set2,1+shift_set2:nset2+shift_set2)

!DEBUG
!      write(std_out,*)' The matrix res_mn(:,1:nset2,1:nset2) is :'
!      write(std_out,*)res_mn(:,1:nset2,1:nset2)
!      call flush(std_out)
!ENDDEBUG

!  Solve R_mn \alpha_n = 1_m
   ABI_ALLOCATE(ipiv,(nset2))
   ABI_ALLOCATE(coeffs,(nset2))
   coeffs(:)=cone
!  The res_mn is destroyed by the following inverse call
   call zgesv(nset2,1,res_mn,wfmixalg-1,ipiv,coeffs,nset2,ierr)
   ABI_DEALLOCATE(ipiv)
!  The coefficients must sum to one
   sum_coeffs=sum(coeffs)
   coeffs=coeffs/sum_coeffs

!DEBUG
!      write(std_out,*)' The coefficients that minimize the residual have been found'
!      write(std_out,*)' coeffs =',coeffs
!      call flush(std_out)
!ENDDEBUG
 end if ! wfmixalg>2 and istep>2

 if(wfmixalg>2 .and. istep>1)then

!  Find the new "input" wavefunction, bi-orthogonalized, and store it replacing the adequate "old" input wavefunction.

   icg=0
   icg_hist=0
   ibg=0
   ibg_hist=0
   ABI_ALLOCATE(al,(2,nset2))
   if(istep>2)then
     do iset2=1,nset2
       al(1,iset2)=real(coeffs(iset2)) ; al(2,iset2)=aimag(coeffs(iset2))
     end do
   else
     al(1,1)=one ; al(2,1)=zero
   end if

!DEBUG
!      write(std_out,*)' Overload the coefficients, in order to simulate a simple mixing with wfmix '
!      write(std_out,*)' Set al(1,ind_biorthog-3)=one, for ind_biorthog=',ind_biorthog
!      write(std_out,*)' This will feed scf_history for set ind_biorthog-3+wfmixalg=',ind_biorthog-3+wfmixalg
!      al(:,:)=zero
!      al(1,ind_biorthog-3)=one
!      call flush(std_out)
!ENDDEBUG

!  LOOP OVER SPINS
   do isppol=1,dtset%nsppol

!    BIG FAT k POINT LOOP
     do ikpt=1,dtset%nkpt

!      Select k point to be treated by this proc
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
       if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) cycle

       istwf_k=dtset%istwfk(ikpt)
       npw_k=npwarr(ikpt)

       if(istep>2)then
!        Make the appropriate linear combination (from the extrapolated wfs)
         cg(:,icg+1:icg+my_nspinor*npw_k*nband_k)=zero
         do iset2=1,nset2
           cg(1,icg+1:icg+my_nspinor*npw_k*nband_k)=cg(1,icg+1:icg+my_nspinor*npw_k*nband_k)&
&           +al(1,iset2)*scf_history_wf%cg(1,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,iset2+wfmixalg)&
&           -al(2,iset2)*scf_history_wf%cg(2,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,iset2+wfmixalg)
           cg(2,icg+1:icg+my_nspinor*npw_k*nband_k)=cg(2,icg+1:icg+my_nspinor*npw_k*nband_k)&
&           +al(1,iset2)*scf_history_wf%cg(2,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,iset2+wfmixalg)&        
&           +al(2,iset2)*scf_history_wf%cg(1,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,iset2+wfmixalg)
         end do
       else ! One needs a simple copy from the extrapolated wavefunctions
         cg(:,icg+1:icg+my_nspinor*npw_k*nband_k)=scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,1+wfmixalg)
       end if
!      Note the storage in cprj_k. By the way, a simple copy might also be used in case istep=2.
       if(usepaw==1) then
         do ibdsp=1,my_nspinor*nbdmix
           call pawcprj_lincom(al,scf_history_wf%cprj(:,ibdsp,1+wfmixalg:nset2+wfmixalg),cprj_k(:,ibdsp:ibdsp),nset2)
         end do
       end if

!      Store the newly extrapolated wavefunctions for this k point, still bi-orthonormalized, in scf_history_wf
       scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,ind_newwf)=cg(:,icg+1:icg+my_nspinor*npw_k*nbdmix)
       if(usepaw==1) then
         call pawcprj_put(atindx1,cprj_k,scf_history_wf%cprj(:,:,ind_newwf),dtset%natom,1,ibg_hist,ikpt,1,isppol,&
&         nbdmix,dtset%mkmem,dtset%natom,nbdmix,nbdmix,dimcprj,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
       end if

!      Back to usual orthonormalization for the cg and cprj_k
       call cgcprj_cholesky(atindx1,cg,cprj_k,dimcprj,icg,ikpt,isppol,istwf_k,mcg,my_nspinor*nband_k,dtset%mkmem,&
&       mpi_enreg,dtset%natom,nattyp,nbdmix,npw_k,my_nspinor,dtset%nsppol,ntypat,pawtab,usepaw)

!      Need to transfer cprj_k to cprj
       if(usepaw==1) then
         call pawcprj_put(atindx1,cprj_k,cprj,dtset%natom,1,ibg_hist,ikpt,1,isppol,&
&         nbdmix,dtset%mkmem,dtset%natom,nbdmix,nbdmix,dimcprj,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
       end if

       ibg=ibg+my_nspinor*nband_k
       ibg_hist=ibg_hist+my_nspinor*nbdmix
       icg=icg+my_nspinor*nband_k*npw_k
       icg_hist=icg_hist+my_nspinor*nbdmix*npw_k

!      End big k point loop
     end do
!    End loop over spins
   end do

   if(istep>2)then
     ABI_DEALLOCATE(coeffs)
   end if
   ABI_DEALLOCATE(al)

 end if ! wfmixalg>2 and istep>1

!DEBUG
! write(std_out,*)' wf_mixing : exit '
!      write(std_out,*)' scf_history_wf%dotprod_sumdiag_cgcprj_ij(:,2,2)=',&
!&       scf_history_wf%dotprod_sumdiag_cgcprj_ij(:,2,2)
! write(std_out,*)' cg(1:2,1:2)=',cg(1:2,1:2)
! write(std_out,*)' scf_history_wf%cg(1:2,1:2,1)=',scf_history_wf%cg(1:2,1:2,1)
! ABI_DEALLOCATE(cg_ref)
! ABI_DATATYPE_DEALLOCATE(cprj_ref)
!ENDDEBUG

 if(usepaw==1) then
   call pawcprj_free(cprj_k)
   call pawcprj_free(cprj_kh)
 end if
 ABI_DATATYPE_DEALLOCATE(cprj_k)
 ABI_DATATYPE_DEALLOCATE(cprj_kh)
 ABI_DEALLOCATE(dimcprj)
 ABI_DEALLOCATE(mmn)
 ABI_DEALLOCATE(smn)
 if(wfmixalg>2)then
   ABI_DEALLOCATE(dotprod_res_k)
   ABI_DEALLOCATE(dotprod_res)
   ABI_DEALLOCATE(res_mn)
 end if
end subroutine wf_mixing
!!***
