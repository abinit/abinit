!{\src2tex{textfont=tt}}
!!****f* ABINIT/extrapwf
!!
!! NAME
!! extrapwf
!!
!! FUNCTION
!! Extrapolate wavefunctions for new ionic positions
!! from values of wavefunctions of previous SCF cycle.
!! Use algorithm proposed by T. A.  Arias et al. in PRB 45, 1538 (1992)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (MT,FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  istep=number of call the routine
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=information about MPI parallelization
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  ngfft(18)=contain all needed information about 3D FFT
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  ntypat=number of types of atoms in cell
!!  pawtab(ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps<type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  xred_old(3,natom)=old reduced coordinates for atoms in unit cell
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!
!! SIDE EFFECTS
!!  cg(2,mcg)= plane wave wavefunction coefficient
!!                          Value from previous SCF cycle is input
!!                          Extrapolated value is output
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!
!! NOTES
!!  THIS ROUTINE IS NOT USEABLE AT PRESENT.
!!  SHOULD BE CAREFULY TESTED AND DEBUGGED (ESPECIALLY WITHIN PAW).
!!
!! PARENTS
!!      extraprho
!!
!! CHILDREN
!!      ctocprj,dotprod_g,getph,hermit,metric,pawcprj_alloc,pawcprj_copy
!!      pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_lincom,pawcprj_put
!!      pawcprj_zaxpby,xmpi_allgather,xmpi_alltoallv,zgemm,zhpev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine extrapwf(atindx,atindx1,cg,dtset,istep,kg,mcg,mgfft,mpi_enreg,&
& nattyp,ngfft,npwarr,ntypat,pawtab,psps,rprimd,scf_history,usepaw,xred_old,ylm)

 use defs_basis
 use defs_abitypes
 use m_scf_history
 use m_xmpi
 use m_profiling_abi
 use m_errors
 use m_cgtools

 use defs_datatypes, only : pseudopotential_type
 use m_pawtab, only : pawtab_type
 use m_pawcprj, only : pawcprj_type, pawcprj_alloc, pawcprj_copy, pawcprj_get, pawcprj_lincom, &
&                      pawcprj_free, pawcprj_zaxpby, pawcprj_put, pawcprj_getdim

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'extrapwf'
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_56_recipspace
 use interfaces_66_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep,mcg,mgfft,ntypat,usepaw
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(scf_history_type),intent(inout) :: scf_history
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom),kg(3,dtset%mpw*dtset%mkmem),nattyp(ntypat),ngfft(18)
 integer,intent(in) :: npwarr(dtset%nkpt)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(inout) :: cg(2,mcg)
 real(dp),intent(in) :: xred_old(3,dtset%natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: ia,iat,iatom,iband_max,iband_max1,iband_min,iband_min1,ibd,ibg,iblockbd,iblockbd1,icg,icgb,icgb1,icgb2
 integer :: ierr,ig,ii,ikpt,ilmn1,ilmn2,inc,ind1,ind2,iorder_cprj
 integer :: isize,isppol,istep1,istwf_k,itypat,klmn,me_distrb,my_nspinor
 integer :: nband_k,nblockbd,nprocband,npw_k,npw_nk,spaceComm_band
 real(dp) :: dotr,dotr1,doti,doti1,eigval
 !character(len=500) :: message
!arrays
 real(dp) :: alpha(2),beta(2),gmet(3,3),gprimd(3,3),rmet(3,3),ph1d(2,3*(2*mgfft+1)*dtset%natom),ucvol
 integer,allocatable :: bufsize(:),bufsize_wf(:),bufdisp(:),bufdisp_wf(:),dimcprj(:),npw_block(:),npw_disp(:)
 real(dp),allocatable :: al(:,:),anm(:),cwavef(:,:),cwavef1(:,:),cwavef_tmp(:,:),deltawf1(:,:),deltawf2(:,:)
 real(dp),allocatable :: eig(:),evec(:,:)
 real(dp),allocatable :: unm(:,:,:)
 real(dp),allocatable :: work(:,:),work1(:,:),wf1(:,:),ylmgr_k(:,:,:),zhpev1(:,:),zhpev2(:)
 complex(dpc),allocatable :: unm_tmp(:,:),anm_tmp(:,:)
 type(pawcprj_type),allocatable :: cprj(:,:),cprj_k(:,:),cprj_k1(:,:),cprj_k2(:,:),cprj_k3(:,:),cprj_k4(:,:)
!complex(dpc) :: aa

! *************************************************************************

 if (istep==0) return

!Useful array
 if (usepaw==1) then
   ABI_ALLOCATE(dimcprj,(dtset%natom))
   call pawcprj_getdim(dimcprj,dtset%natom,nattyp,ntypat,dtset%typat,pawtab,'O')
 end if

!Metric
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!History indexes
 ind1=1;ind2=2

!First step
 if (istep==1) then
   scf_history%cg(:,:,ind1)=cg(:,:)
!  scf_history%cg(:,:,ind2)=zero
   scf_history%cg(:,:,ind2)= cg(:,:)
   if(usepaw==1) then
!    WARNING: THIS SECTION IS USELESS; NOW crpj CAN BE READ FROM SCFCV
     call getph(atindx,dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,xred_old)
     iatom=0 ; iorder_cprj=0
     call pawcprj_alloc(scf_history%cprj(:,:,ind1),0,dimcprj)
     call pawcprj_alloc(scf_history%cprj(:,:,ind2),0,dimcprj)
     ABI_ALLOCATE(ylmgr_k,(dtset%mpw,3,0))
     call ctocprj(atindx,cg,1,scf_history%cprj(:,:,ind1),gmet,gprimd,&
&     iatom,0,iorder_cprj,dtset%istwfk,kg,dtset%kptns,dtset%mband,mcg,scf_history%mcprj,&
&     dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,&
&     dtset%natom,nattyp,dtset%nband,dtset%natom,ngfft,dtset%nkpt,&
&     dtset%nloalg,npwarr,dtset%nspinor,dtset%nsppol,dtset%ntypat,&
&     dtset%paral_kgb,ph1d,psps,rmet,dtset%typat,ucvol,0,&
&     xred_old,ylm,ylmgr_k)
     ABI_DEALLOCATE(ylmgr_k)
!    call pawcprj_set_zero(scf_history%cprj(:,:,ind2))
     call pawcprj_copy(scf_history%cprj(:,:,ind1),scf_history%cprj(:,:,ind2))
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

!  For the moment sequential part only
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

   if(usepaw==1) then
!    WARNING: THIS SECTION IS USELESS; NOW cprj CAN BE READ FROM SCFCV
     call getph(atindx,dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,xred_old)
     ABI_DATATYPE_ALLOCATE(cprj,(dtset%natom,scf_history%mcprj))
     call pawcprj_alloc(cprj,0,dimcprj)
     iatom=0 ; iorder_cprj=0
     ABI_ALLOCATE(ylmgr_k,(dtset%mpw,3,0))
     call ctocprj(atindx,cg,1,cprj,gmet,gprimd,iatom,0,iorder_cprj,&
&     dtset%istwfk,kg,dtset%kptns,dtset%mband,mcg,scf_history%mcprj,dtset%mgfft,&
&     dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%natom,&
&     nattyp,dtset%nband,dtset%natom,ngfft,dtset%nkpt,dtset%nloalg,&
&     npwarr,dtset%nspinor,dtset%nsppol,dtset%ntypat,dtset%paral_kgb,&
&     ph1d,psps,rmet,dtset%typat,ucvol,0,xred_old,&
&     ylm,ylmgr_k)
     ABI_DEALLOCATE(ylmgr_k)
   end if  ! end usepaw=1

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
       ABI_ALLOCATE(cwavef1,(2,npw_nk*my_nspinor))
       if (nprocband>1) then
         isize=2*my_nspinor;bufsize(:)=isize*npw_block(:);bufdisp(:)=isize*npw_disp(:)
         isize=2*my_nspinor*npw_k;bufsize_wf(:)=isize
         do ii=1,nprocband
           bufdisp_wf(ii)=(ii-1)*isize
         end do
       end if

!      Subspace alignment

!      Loop over bands or blocks of bands
       nblockbd=nband_k/nprocband
       icgb=icg

       if(usepaw==1) then
         ABI_DATATYPE_ALLOCATE( cprj_k,(dtset%natom,my_nspinor*nblockbd))
         call pawcprj_alloc(cprj_k,cprj(1,1)%ncpgr,dimcprj)
         call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,ibg,ikpt,1,isppol,dtset%mband,&
&         dtset%mkmem,dtset%natom,nblockbd,nblockbd,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         ABI_DATATYPE_ALLOCATE( cprj_k1,(dtset%natom,my_nspinor*nblockbd))
         call pawcprj_alloc(cprj_k1,scf_history%cprj(1,1,ind1)%ncpgr,dimcprj)
         call pawcprj_get(atindx1,cprj_k1,scf_history%cprj(:,:,ind1),dtset%natom,1,ibg,ikpt,1,isppol,&
&         dtset%mband,dtset%mkmem,dtset%natom,nblockbd,nblockbd,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         ABI_DATATYPE_ALLOCATE( cprj_k2,(dtset%natom,my_nspinor*nblockbd))
         call pawcprj_alloc(cprj_k2,scf_history%cprj(1,1,ind2)%ncpgr,dimcprj)
         call pawcprj_get(atindx1,cprj_k2,scf_history%cprj(:,:,ind2),dtset%natom,1,ibg,ikpt,1,isppol,&
&         dtset%mband,dtset%mkmem,dtset%natom,nblockbd,nblockbd,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
       end if  !end usepaw=1

       ABI_ALLOCATE(unm,(2,nblockbd,nblockbd))
       unm=zero
       icgb2=0

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
               cwavef_tmp(1,ig)=scf_history%cg(1,ig+icgb1,ind1)
               cwavef_tmp(2,ig)=scf_history%cg(2,ig+icgb1,ind1)
             end do
             call xmpi_alltoallv(cwavef_tmp,bufsize_wf,bufdisp_wf,cwavef1,bufsize,bufdisp,spaceComm_band,ierr)
             ABI_DEALLOCATE(cwavef_tmp)
           else
             do ig=1,npw_k*my_nspinor
               cwavef1(1,ig)=scf_history%cg(1,ig+icgb1,ind1)
               cwavef1(2,ig)=scf_history%cg(2,ig+icgb1,ind1)
             end do
           end if

!          Calculate Unm=<psi_nk(t)|S|psi_mk(t-dt)>
           call dotprod_g(dotr,doti,istwf_k,npw_k*my_nspinor,2,cwavef,cwavef1,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
           if(usepaw==1) then
             ia =0
             do itypat=1,ntypat
               do iat=1+ia,nattyp(itypat)+ia
                 do ilmn1=1,pawtab(itypat)%lmn_size
                   do ilmn2=1,ilmn1
                     klmn=((ilmn1-1)*ilmn1)/2+ilmn2
                     dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k1(iat,iblockbd1)%cp(1,ilmn2)+&
&                     cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k1(iat,iblockbd1)%cp(2,ilmn2))
                     doti=doti+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k1(iat,iblockbd1)%cp(2,ilmn2)-&
&                     cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k1(iat,iblockbd1)%cp(1,ilmn2))
                   end do
                   do ilmn2=ilmn1+1,pawtab(itypat)%lmn_size
                     klmn=((ilmn2-1)*ilmn2)/2+ilmn1
                     dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k1(iat,iblockbd1)%cp(1,ilmn2)+&
&                     cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k1(iat,iblockbd1)%cp(2,ilmn2))
                     doti=doti+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k1(iat,iblockbd1)%cp(2,ilmn2)-&
&                     cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k1(iat,iblockbd1)%cp(1,ilmn2))
                   end do
                 end do
               end do
               ia=ia+nattyp(itypat)
             end do
           end if
!          unm(1,iblockbd,iblockbd1)=dotr
!          unm(2,iblockbd,iblockbd1)=doti
           unm(1,iblockbd1,iblockbd)=dotr
           unm(2,iblockbd1,iblockbd)=doti
!          End loop over bands iblockbd1
           icgb1=icgb1+npw_k*my_nspinor*nprocband

         end do

!        End loop over bands iblockbd
         icgb2=icgb2+npw_k*my_nspinor*nprocband
         icgb=icgb+npw_k*my_nspinor*nprocband
       end do

!      write(std_out,*) 'UNM'
!      do iblockbd=1,nblockbd
!      write(std_out,11) (unm(1,iblockbd,iblockbd1),unm(2,iblockbd,iblockbd1),iblockbd1=1,nblockbd)
!      end do
!      11 format(12(1x,f9.5),a)
!      Compute A=tU^*U
       ABI_ALLOCATE(unm_tmp,(nblockbd,nblockbd))
       ABI_ALLOCATE(anm_tmp,(nblockbd,nblockbd))
       ABI_ALLOCATE(anm,(nblockbd*(nblockbd+1)))
       unm_tmp(:,:)=cmplx(unm(1,:,:),unm(2,:,:),kind=dp)
       call zgemm('C','N',nblockbd,nblockbd,nblockbd,dcmplx(1._dp), unm_tmp,nblockbd, &
&       unm_tmp,nblockbd,dcmplx(0._dp),anm_tmp,nblockbd)
       do iblockbd=1,nblockbd
         do iblockbd1=iblockbd,nblockbd
           ii=iblockbd1*(iblockbd1-1)+2*(iblockbd-1)+1
           anm(ii)=real(anm_tmp(iblockbd,iblockbd1))
           anm(ii+1)=aimag(anm_tmp(iblockbd,iblockbd1))
         end do
       end do
       call hermit(anm,anm,ierr,nblockbd)
!      aa=dcmplx(0._dp)
!      do iblockbd=1,nblockbd
!      aa=aa+conjg(unm_tmp(iblockbd,1))*unm_tmp(iblockbd,1)
!      end do
!      write(std_out,*) 'tU*U', aa
!      write(std_out,*) 'ANM_tmp'
!      do iblockbd=1,nblockbd
!      write(std_out,11) (anm_tmp(iblockbd,iblockbd1),iblockbd1=1,nblockbd)
!      end do
!      write(std_out,*) 'ANM'
!      do iblockbd=1,nblockbd*(nblockbd+1)
!      write(std_out,11) anm(iblockbd)
!      end do

!      Diagonalize A
       ABI_ALLOCATE(eig,(nblockbd))
       ABI_ALLOCATE(evec,(2*nblockbd,nblockbd))
       ABI_ALLOCATE(zhpev1,(2,2*nblockbd-1))
       ABI_ALLOCATE(zhpev2,(3*nblockbd-2))
       call zhpev('V','U',nblockbd,anm,eig,evec,nblockbd,zhpev1,&
&       zhpev2,ierr)
       ABI_DEALLOCATE(anm)
       ABI_DEALLOCATE(zhpev1)
       ABI_DEALLOCATE(zhpev2)
!      aa=dcmplx(0._dp)
!      do iblockbd=1,nblockbd
!      aa=aa+anm_tmp(1,iblockbd)*cmplx(evec((2*iblockbd-1),1),evec(2*iblockbd,1),kind=dp)
!      end do
!      write(std_out,*) 'EIG', aa, eig(1)*evec(1,1),eig(1)*evec(2,1)

!      Compute A'=evec*tU^/sqrt(eig)
       call zgemm('C','C',nblockbd,nblockbd,nblockbd,dcmplx(1._dp),evec,nblockbd, &
&       unm_tmp,nblockbd,dcmplx(0._dp),anm_tmp,nblockbd)
       do iblockbd=1,nblockbd
         eigval=dsqrt(eig(iblockbd))
         do iblockbd1=1,nblockbd
           anm_tmp(iblockbd,iblockbd1)=anm_tmp(iblockbd,iblockbd1)/eigval
         end do
       end do

!      Compute tA^A'to come back to the initial subspace for the cg's

       call zgemm('N','N',nblockbd,nblockbd,nblockbd,dcmplx(1._dp),evec,nblockbd, &
&       anm_tmp,nblockbd,dcmplx(0._dp),unm_tmp,nblockbd)
       anm_tmp=unm_tmp
!      write(std_out,*) 'ANM_tmp'
!      do iblockbd=1,nblockbd
!      write(std_out,11) (anm_tmp(iblockbd,iblockbd1),iblockbd1=1,nblockbd)
!      end do

!      Wavefunction alignment (istwfk=1 ?)
       ABI_ALLOCATE(work,(2,npw_nk*my_nspinor*nblockbd))
       ABI_ALLOCATE(work1,(2,my_nspinor*nblockbd*npw_nk))
       work1(:,:)=scf_history%cg(:,icg+1:icg+my_nspinor*nblockbd*npw_nk,ind1)
       call zgemm('N','N',npw_nk*my_nspinor,nblockbd,nblockbd,dcmplx(1._dp), &
&       work1,npw_nk*my_nspinor, &
&       anm_tmp,nblockbd,dcmplx(0._dp),work,npw_nk*my_nspinor)
       scf_history%cg(:,1+icg:npw_nk*my_nspinor*nblockbd+icg,ind1)=work(:,:)

       work1(:,:)=scf_history%cg(:,icg+1:icg+my_nspinor*nblockbd*npw_nk,ind2)
       call zgemm('N','N',npw_nk*my_nspinor,nblockbd,nblockbd,dcmplx(1._dp), &
&       work1,npw_nk*my_nspinor, &
&       anm_tmp,nblockbd,dcmplx(0._dp),work,npw_nk*my_nspinor)
       scf_history%cg(:,1+icg:npw_nk*my_nspinor*nblockbd+icg,ind2)=work(:,:)
       ABI_DEALLOCATE(work1)
!      If paw, must also align cprj:
       if (usepaw==1) then
!        New version (MT):
         ABI_DATATYPE_ALLOCATE(cprj_k3,(dtset%natom,my_nspinor))
         ABI_DATATYPE_ALLOCATE(cprj_k4,(dtset%natom,my_nspinor))
         call pawcprj_alloc(cprj_k3,cprj_k1(1,1)%ncpgr,dimcprj)
         call pawcprj_alloc(cprj_k4,cprj_k2(1,1)%ncpgr,dimcprj)
         ABI_ALLOCATE(al,(2,nblockbd))
         do iblockbd=1,nblockbd
           ii=(iblockbd-1)*my_nspinor
           do iblockbd1=1,nblockbd
             al(1,iblockbd1)=real (anm_tmp(iblockbd,iblockbd1))
             al(2,iblockbd1)=aimag(anm_tmp(iblockbd,iblockbd1))
           end do
           call pawcprj_lincom(al,cprj_k1,cprj_k3,nblockbd)
           call pawcprj_lincom(al,cprj_k2,cprj_k4,nblockbd)
           call pawcprj_copy(cprj_k3,cprj_k1(:,ii+1:ii+my_nspinor))
           call pawcprj_copy(cprj_k4,cprj_k2(:,ii+1:ii+my_nspinor))
         end do
         ABI_DEALLOCATE(al)
!        Old version (FJ):
!        allocate( cprj_k3(dtset%natom,my_nspinor*nblockbd))
!        call pawcprj_alloc(cprj_k3,cprj_k1(1,1)%ncpgr,dimcprj)
!        allocate( cprj_k4(dtset%natom,my_nspinor*nblockbd))
!        call pawcprj_alloc(cprj_k4,cprj_k2(1,1)%ncpgr,dimcprj)
!        beta(1)=one;beta(2)=zero
!        do iblockbd=1,nblockbd*my_nspinor
!        do iblockbd1=1,nblockbd*my_nspinor
!        alpha(1)=real(anm_tmp(iblockbd,iblockbd1));alpha(2)=aimag(anm_tmp(iblockbd,iblockbd1))
!        call pawcprj_zaxpby(alpha,beta,cprj_k1(:,iblockbd1:iblockbd1),cprj_k3(:,iblockbd:iblockbd))
!        call pawcprj_zaxpby(alpha,beta,cprj_k2(:,iblockbd1:iblockbd1),cprj_k4(:,iblockbd:iblockbd))
!        end do
!        end do
!        call pawcprj_copy(cprj_k3,cprj_k1)
!        call pawcprj_copy(cprj_k4,cprj_k2)

         call pawcprj_free(cprj_k3)
         call pawcprj_free(cprj_k4)
         ABI_DATATYPE_DEALLOCATE(cprj_k3)
         ABI_DATATYPE_DEALLOCATE(cprj_k4)
       end if
       ABI_DEALLOCATE(anm_tmp)
       ABI_DEALLOCATE(unm_tmp)
       ABI_DEALLOCATE(work)

!      Wavefunction extrapolation
       ibd=0
       inc=npw_nk*my_nspinor
       ABI_ALLOCATE(deltawf2,(2,npw_nk*my_nspinor))
       ABI_ALLOCATE(wf1,(2,npw_nk*my_nspinor))
       ABI_ALLOCATE(deltawf1,(2,npw_nk*my_nspinor))
       do iblockbd=1,nblockbd
         deltawf2(:,:)=scf_history%cg(:,1+icg+ibd:icg+ibd+inc,ind2)
         wf1(:,:)=scf_history%cg(:,1+icg+ibd:icg+ibd+inc,ind1)
!        wf1(2,1)=zero;deltawf2(2,1)=zero

         call dotprod_g(dotr,doti,istwf_k,npw_nk*my_nspinor,2,cg(:,icg+1+ibd:ibd+icg+inc),cg(:,icg+1+ibd:ibd+icg+inc),&
&         mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
         call dotprod_g(dotr1,doti1,istwf_k,npw_nk*my_nspinor,2,cg(:,icg+1+ibd:ibd+icg+inc),wf1,&
&         mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
         if(usepaw==1) then
           ia =0
           do itypat=1,ntypat
             do iat=1+ia,nattyp(itypat)+ia
               do ilmn1=1,pawtab(itypat)%lmn_size
                 do ilmn2=1,ilmn1
                   klmn=((ilmn1-1)*ilmn1)/2+ilmn2
                   dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k(iat,iblockbd)%cp(1,ilmn2)+&
&                   cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k(iat,iblockbd)%cp(2,ilmn2))
                   doti=doti+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k(iat,iblockbd)%cp(2,ilmn2)-&
&                   cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k(iat,iblockbd)%cp(1,ilmn2))
                   dotr1=dotr1+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k1(iat,iblockbd)%cp(1,ilmn2)+&
&                   cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k1(iat,iblockbd)%cp(2,ilmn2))
                   doti1=doti1+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k1(iat,iblockbd)%cp(2,ilmn2)-&
&                   cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k1(iat,iblockbd)%cp(1,ilmn2))
                 end do
                 do ilmn2=ilmn1+1,pawtab(itypat)%lmn_size
                   klmn=((ilmn2-1)*ilmn2)/2+ilmn1
                   dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k(iat,iblockbd)%cp(1,ilmn2)+&
&                   cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k(iat,iblockbd)%cp(2,ilmn2))
                   doti=doti+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k(iat,iblockbd)%cp(2,ilmn2)-&
&                   cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k(iat,iblockbd)%cp(1,ilmn2))
                   dotr1=dotr1+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k1(iat,iblockbd)%cp(1,ilmn2)+&
&                   cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k1(iat,iblockbd)%cp(2,ilmn2))
                   doti1=doti1+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k1(iat,iblockbd)%cp(2,ilmn2)-&
&                   cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k1(iat,iblockbd)%cp(1,ilmn2))
                 end do
               end do
             end do
             ia=ia+nattyp(itypat)
           end do
         end if
         dotr=sqrt(dotr**2+doti**2)
         dotr1=sqrt(dotr1**2+doti1**2)
         write(std_out,*)'DOTR, DOTR1',dotr,dotr1
         dotr=dotr1/dotr
         write(std_out,*)'DOTR',dotr
         deltawf1=zero
         if(dotr>=0.9d0) then
           deltawf1(:,:)=cg(:,icg+1+ibd:ibd+icg+inc)-wf1(:,:)
           if(usepaw==1) then
             alpha(1)=one;alpha(2)=zero
             beta(1)=-one;beta(2)=zero
             ia =0
             call pawcprj_zaxpby(alpha,beta,cprj_k(:,iblockbd:iblockbd),cprj_k1(:,iblockbd:iblockbd))
           end if
           istep1=istep
         else
           istep1=1
         end if
         scf_history%cg(:,1+icg+ibd:icg+ibd+inc,ind1)=cg(:,icg+1+ibd:ibd+icg+inc)
         scf_history%cg(:,1+icg+ibd:icg+ibd+inc,ind2)=deltawf1(:,:)
         if(usepaw==1) then
           call pawcprj_put(atindx1,cprj_k,scf_history%cprj(:,:,ind1),dtset%natom,1,ibg,ikpt,1,isppol,&
&           dtset%mband,dtset%mkmem,dtset%natom,nblockbd,nblockbd,dimcprj,my_nspinor,dtset%nsppol,0,&
&           mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
           call pawcprj_put(atindx1,cprj_k1,scf_history%cprj(:,:,ind2),dtset%natom,1,ibg,ikpt,1,isppol,&
&           dtset%mband,dtset%mkmem,dtset%natom,nblockbd,nblockbd,dimcprj,my_nspinor,dtset%nsppol,0,&
&           mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
         end if

!        if(istep1>=3) then
         cg(:,icg+1+ibd:ibd+icg+inc)=cg(:,icg+1+ibd:ibd+icg+inc)+scf_history%alpha*deltawf1(:,:) &
&         +scf_history%beta *deltawf2(:,:)

!        to be used later
!        if(usepaw==1) then
!        alpha(2)=zero
!        beta(1)=one;beta(2)=zero
!        alpha(1)=scf_history%alpha
!        call pawcprj_zaxpby(alpha,beta,cprj_k1(:,iblockbd:iblockbd),cprj_k(:,iblockbd:iblockbd))
!        alpha(1)=scf_history%beta
!        call pawcprj_zaxpby(alpha,beta,cprj_k2(:,iblockbd:iblockbd),cprj_k(:,iblockbd:iblockbd))
!        call pawcprj_put(atindx1,cprj_k,cprj,dtset%natom,1,ibg,ikpt,1,isppol,&
!        &    dtset%mband,dtset%mkmem,dtset%natom,nblockbd,nblockbd,dimcprj,my_nspinor,dtset%nsppol,0,&
!        &    mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
!        end if
!        else if (istep1==2) then
!          cg(:,icg+1+ibd:ibd+icg+inc)=cg(:,icg+1+ibd:ibd+icg+inc)+scf_history%alpha*deltawf1(:,:)+scf_history%beta*wf1(:,:)
!       !     cg(:,icg+1+ibd:ibd+icg+inc)=cg(:,icg+1+ibd:ibd+icg+inc)+deltawf1(:,:)
!        if(usepaw==1) then
!        alpha(2)=zero
!        beta(1)=one;beta(2)=zero
!        alpha(1)=scf_history%alpha
!        call pawcprj_zaxpby(alpha,beta,cprj_k1(:,iblockbd:iblockbd),cprj_k(:,iblockbd:iblockbd))
!        alpha(1)=scf_history%beta
!        call pawcprj_zaxpby(alpha,beta,cprj_k2(:,iblockbd:iblockbd),cprj_k(:,iblockbd:iblockbd))
!        call pawcprj_put(atindx1,cprj_k,cprj,dtset%natom,1,ibg,ikpt,1,isppol,&
!        &    dtset%mband,dtset%mkmem,dtset%natom,nblockbd,nblockbd,dimcprj,my_nspinor,dtset%nsppol,0,&
!        &    mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
!        end if
!        end if
         ibd=ibd+inc
       end do ! end loop on iblockbd

       ABI_DEALLOCATE(deltawf1)
       ABI_DEALLOCATE(deltawf2)
       ABI_DEALLOCATE(wf1)
       ABI_DEALLOCATE(cwavef)
       ABI_DEALLOCATE(cwavef1)
       ABI_DEALLOCATE(eig)
       ABI_DEALLOCATE(evec)
       ABI_DEALLOCATE(unm)
       if(usepaw==1) then
         call pawcprj_free(cprj_k)
         ABI_DATATYPE_DEALLOCATE(cprj_k)
         call pawcprj_free(cprj_k1)
         ABI_DATATYPE_DEALLOCATE(cprj_k1)
         call pawcprj_free(cprj_k2)
         ABI_DATATYPE_DEALLOCATE(cprj_k2)
       end if

       ibg=ibg+my_nspinor*nband_k
       icg=icg+my_nspinor*nband_k*npw_k

!      End big k point loop
     end do
!    End loop over spins
   end do

   if(usepaw==1) then
     call pawcprj_free(cprj)
     ABI_DATATYPE_DEALLOCATE(cprj)
   end if
   if (nprocband>1) then
     ABI_DEALLOCATE(npw_block)
     ABI_DEALLOCATE(npw_disp)
     ABI_DEALLOCATE(bufsize)
     ABI_DEALLOCATE(bufdisp)
     ABI_DEALLOCATE(bufsize_wf)
     ABI_DEALLOCATE(bufdisp_wf)
   end if

 end if ! istep>=2

 if (usepaw==1) then
   ABI_DEALLOCATE(dimcprj)
 end if

end subroutine extrapwf
!!***
