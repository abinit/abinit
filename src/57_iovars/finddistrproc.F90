!{\src2tex{textfont=tt}}
!!****f* ABINIT/finddistrproc
!! NAME
!! finddistrproc
!!
!! FUNCTION
!!   Given a total number of processors, find a suitable distribution
!!   that fill all the different levels of parallelization
!!   (npimage, nppert, npkpt, npspinor, npband, npfft, bandpp)
!!   Also determine parameters of parallel Linear Algebra routines
!!   (use_slk, np_slk, gpu_linalg_limit)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (FJ,MT,FD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables,
!!   for all datasets; at this stage only datasets with index lower than
!!   idtset are already initalized
!!  filnam(5)=character strings giving file names
!!  idtset=number of the current dataset
!!  mpi_enreg=information about MPI parallelization
!!  mband=maximum number of bands.
!!  ndtset_alloc=number of datasets, corrected for allocation of at least one data set
!!  tread(11)=flags indicating wether parallel input parameters were read from input file
!!            tread(1)  : paral_kgb      tread(6) : npfft
!!            tread(2)  : npimage        tread(7) : npband
!!            tread(3)  : nppert         tread(8) : bandpp
!!            tread(4)  : npkpt          tread(9) : use_slk
!!            tread(5)  : nspinor        tread(10): np_slk
!!            tread(11) : gpu_linalg_limit
!!
!! SIDE EFFECTS
!!  iexit= if incremented, an exit is required
!!  dtset%paral_kgb= flag for band-fft parallelism
!!  dtset%npimage  = number of processors for parallelisation over image
!!  dtset%nppert   = number of processors for parallelisation over perturbations
!!  dtset%npspinor = number of processors for parallelisation on k points
!!  dtset%npkpt    = number of processors for parallelisation on k points
!!  dtset%npfft    = number of processors for parallelisation on fft grid
!!  dtset%npband   = number of processors for parallelisation on bands
!!  dtset%nphf     = number of processors for parallelisation on occupied states for fock exchange
!!  dtset%bandpp   = internal parameter for lobpcg parallelisation algorithm
!!  dtset%use_slk  = flag for ScalaPAck use
!!  dtset%np_slk   = number of processors used in ScaLapack routines
!!  dtset%gpu_linalg_limit=threshold activating Linear Algebra on GPU
!!
!! PARENTS
!!      mpi_setup
!!
!! CHILDREN
!!      compute_kgb_indicator,get_npert_rbz,hdr_free,hdr_read_from_fname
!!      initmpi_world,kpgcount,metric,mkfilename,mkrdim,sort_dp,wrtout
!!      xmpi_bcast
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine finddistrproc(dtsets,filnam,idtset,iexit,mband,mpi_enreg,ndtset_alloc,tread)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_xomp
 use m_hdr
 use m_sort

 use m_fftcore, only : kpgcount

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'finddistrproc'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_51_manage_mpi
 use interfaces_54_abiutil
 use interfaces_57_iovars, except_this_one => finddistrproc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idtset,mband,ndtset_alloc
 integer,intent(inout) :: iexit
 type(dataset_type),intent(inout),target :: dtsets(0:ndtset_alloc)
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: tread(11)
 character(len=fnlen),intent(in) :: filnam(5)

!Local variables-------------------------------
!scalars
!128 should be a reasonable maximum for npfft (scaling is very poor for npfft>20)
 integer,parameter :: NPFMAX=128
 integer,parameter :: MAXCOUNT=250,MAXBENCH=25,NPF_CUTOFF=20
 integer :: bpp,bpp_max,bpp_min,optdriver,autoparal
 integer :: npi_max,npi_min,npc,npc_max,npc_min
 integer :: npk,npk_max,npk_min,npp_max,npp_min
 integer :: nps,nps_max,nps_min,npf,npf_max,npf_min
 integer :: npb,npb_max,npb_min,max_ncpus,ount
 integer :: work_size,nks_per_proc,tot_ncpus
 integer :: icount,ii,imin,jj,mcount,mcount_eff,mpw
 integer :: n2,n3,ncell_eff,ncount,nimage_eff,nkpt_eff,npert_eff
 integer :: nproc,nproc1,nprocmin,np_slk,use_linalg_gpu,omp_ncpus
 logical,parameter :: new_version=.true.
 logical :: dtset_found,file_found,first_bpp,iam_master
 real(dp):: acc_c,acc_k,acc_kgb,acc_kgb_0,acc_s,ecut_eff,ucvol,weight0
 real(dp):: eff
 character(len=9) :: suffix
 character(len=500) :: message
 character(len=fnlen) :: filden
 type(hdr_type) :: hdr0
!arrays
 integer :: idum(1),idum3(3),ngmax(3),ngmin(3)
 integer,allocatable :: isort(:),jdtset_(:),my_distp(:,:)
 integer,pointer :: nkpt_rbz(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: weight(:)
 real(dp),pointer :: nband_rbz(:,:)
 type(dataset_type),pointer :: dtset
!Cut-off function for npfft
! cutoff(nn)= &
!&    0.2_dp+(one-0.2_dp)*(sin((pi*(nn-NPF_CUTOFF))/(one*(NPFMAX-NPF_CUTOFF))) &
!&                           /((pi*(nn-NPF_CUTOFF))/(one*(NPFMAX-NPF_CUTOFF))))**2

!******************************************************************

 DBG_ENTER("COLL")

!Select current dataset
 dtset => dtsets(idtset)

!Is automatic parallelization activated?
 autoparal = dtset%autoparal
 if (autoparal==0) return


 ! Handy local variables
 iam_master = (mpi_enreg%me==0)
 optdriver = dtset%optdriver
 max_ncpus = dtset%max_ncpus

 if (max_ncpus > 0 .and. autoparal/=0) then
   iexit = iexit + 1 ! will stop in the parent.
 end if

 ! Unit number used for outputting the autoparal sections
 ount = std_out
 ount = ab_out

 ! Small hack: set paral_kgb to -max_ncpus so that I don't have to change the previous implementation.
 !if (dtset%paral_kgb == 1 .and. max_ncpus > 0) then
 !  dtset%paral_kgb = -max_ncpus
 !end if

 if (optdriver==RUNL_GSTATE .and. dtset%paral_kgb==0 .and. &
& max_ncpus>0 .and. autoparal/=0) then
   if (iam_master) then
     ! This corresponds to the simplest algorithm for GS (band-by-band CG)
     ! with distribution of k-points and spin.
     work_size = dtset%nkpt * dtset%nsppol
     write(ount,"(2a)")ch10,"--- !Autoparal"
     write(ount,"(a)")"# Autoparal section for GS run (band-by-band CG method)"
     write(ount,"(a)")   "info:"
     write(ount,"(a,i0)")"    autoparal: ",autoparal
     write(ount,"(a,i0)")"    paral_kgb: ",dtset%paral_kgb
     write(ount,"(a,i0)")"    max_ncpus: ",max_ncpus
     write(ount,"(a,i0)")"    nspinor: ",dtset%nspinor
     write(ount,"(a,i0)")"    nsppol: ",dtset%nsppol
     write(ount,"(a,i0)")"    nkpt: ",dtset%nkpt
     write(ount,"(a,i0)")"    mband: ",mband

     ! List of configurations.
     ! Assuming an OpenMP implementation with perfect speedup!
     write(ount,"(a)")"configurations:"

     do ii=1,max_ncpus
       if (ii > work_size) cycle
       do omp_ncpus=1,xomp_get_max_threads()
         nks_per_proc = work_size / ii
         nks_per_proc = nks_per_proc + MOD(work_size, ii)
         eff = (one * work_size) / (ii * nks_per_proc)

         write(ount,"(a,i0)")"    - tot_ncpus: ",ii * omp_ncpus
         write(ount,"(a,i0)")"      mpi_ncpus: ",ii
         write(ount,"(a,i0)")"      omp_ncpus: ",omp_ncpus
         write(ount,"(a,f12.9)")"      efficiency: ",eff
         !write(ount,"(a,f12.2)")"      mem_per_cpu: ",mempercpu_mb
       end do
     end do
     write(ount,'(a)')"..."
   end if
   ! Return immediately, will stop in the parent.
   iexit = iexit + 1
   RETURN
 end if


 nproc=mpi_enreg%nproc
 !if (xmpi_paral==1.and.dtset%paral_kgb <0) nproc=-dtset%paral_kgb
 if (max_ncpus > 0) nproc = dtset%max_ncpus
 !if (xmpi_paral==1.and.dtset%paral_kgb <0) nproc=dtset%max_ncpus
 if (xmpi_paral==0.and.dtset%paral_kgb>=0) nproc=1

 if (dtset%paral_kgb>=0) then
   if (nproc==1) then
     if (tread(1)==0.or.xmpi_paral==0) dtset%paral_kgb= 0
     if (tread(2)==0.or.xmpi_paral==0) dtset%npimage  = 1
     if (tread(3)==0.or.xmpi_paral==0) dtset%nppert   = 1
     if (tread(4)==0.or.xmpi_paral==0) dtset%npspinor = 1
     if (tread(5)==0.or.xmpi_paral==0) dtset%npkpt    = 1
     if (tread(6)==0.or.xmpi_paral==0) dtset%npfft    = 1
     if (tread(7)==0.or.xmpi_paral==0) dtset%npband   = 1
     if (tread(8)==0.or.xmpi_paral==0) dtset%bandpp   = 1
     if (tread(9)==0.or.xmpi_paral==0) dtset%use_slk  = 0
     if (tread(10)==0.or.xmpi_paral==0) dtset%np_slk  = 1000000
     return
   end if
   if ((dtset%optdriver/=RUNL_GSTATE.and.dtset%optdriver/=RUNL_RESPFN.and.dtset%optdriver/=RUNL_GWLS).or. &
&   (dtset%optdriver==RUNL_GSTATE.and.dtset%usewvl==1)) then
     dtset%paral_kgb= 0
     dtset%npimage  = max(1,dtset%npimage)
     dtset%nppert   = max(1,dtset%nppert)
     dtset%npspinor = max(1,dtset%npspinor)
     dtset%npkpt    = max(1,dtset%npkpt)
     dtset%npfft    = max(1,dtset%npfft)
     dtset%npband   = max(1,dtset%npband)
     dtset%bandpp   = max(1,dtset%bandpp)
     return
   end if
 end if

 nprocmin=2
 if (xmpi_paral==1.and.dtset%paral_kgb>=0) nprocmin=max(2,nproc-100)
 if (max_ncpus > 0 .and. autoparal/=0) nprocmin = 1

!Need the metric tensor
 call mkrdim(dtset%acell_orig(1:3,1),dtset%rprim_orig(1:3,1:3,1),rprimd)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Determine some quantities related to plane waves
!  - Crude estimation of the number of PW
!  - Number of G vectors in each direction
 mpw=0;ngmin=0;ngmax=0
 if (optdriver==RUNL_GSTATE) then
   ecut_eff = dtset%ecut*dtset%dilatmx**2
   mpw = nint(ucvol*((two*ecut_eff)**1.5_dp)/(six*pi**2)) ! Crude estimation
   if (all(dtset%istwfk(1:dtset%nkpt)>1)) mpw=mpw/2+1
   call kpgcount(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,dtset%kpt,ngmax,ngmin,dtset%nkpt)
   write(message,'(a,i8)') ' getmpw sequential formula gave: ',mpw
   call wrtout(std_out,message,'COLL')
 end if

 write(message,'(2a,i0)')  ch10,&
& ' Computing all possible proc distrib for this input with nproc less than ',nproc
 call wrtout(std_out,message,'COLL')

!Parallelization over images
 npi_min=1;npi_max=1;nimage_eff=1
 if (optdriver==RUNL_GSTATE) then
   nimage_eff=dtset%ndynimage
   if (dtset%ntimimage<=1) nimage_eff=dtset%nimage
   npi_min=max(1,dtset%npimage)
   npi_max=min(nproc,nimage_eff)
   if (tread(2)==1) npi_max=dtset%npimage
 end if

!Parallelization over k-points and spin components (GS)
 npk_min=1;npk_max=1;nkpt_eff=0
 if (optdriver==RUNL_GSTATE) then
   nkpt_eff=dtset%nkpt*dtset%nsppol
   npk_min=max(1,dtset%npkpt)
   npk_max=min(nproc,nkpt_eff)
   if (tread(4)==1) npk_max=dtset%npkpt
 end if

!Parallelization over perturbations, k-points and spin components (DFPT)
 npp_min=1;npp_max=1;npert_eff=1
 if (optdriver==RUNL_RESPFN) then
   if (dtset%paral_rf==1) then
     call get_npert_rbz(dtset,nband_rbz,nkpt_rbz,npert_eff)
     do jj=1,npert_eff
       ii=dtset%nsppol*nkpt_rbz(jj)*maxval(nband_rbz(:,jj))
       nkpt_eff=max(nkpt_eff,ii)
     end do
     npp_min=max(1,dtset%nppert)
     npp_max=min(nproc,npert_eff)
     if (tread(3)==1) then
       npp_max=dtset%nppert
       if (npp_max>npert_eff) then
         npp_min=npert_eff;npp_max=npert_eff
         message='nppert is bigger than npert; we set nppert=npert'
         MSG_WARNING(message)
       end if
     end if
     npk_min=1
     npk_max=min(nproc,nkpt_eff)
     ABI_DEALLOCATE(nkpt_rbz)
     ABI_DEALLOCATE(nband_rbz)
   else
     nkpt_eff=nproc
     npk_min=nproc-5
     npk_max=nproc
   end if
 end if

!Parallelization over spinorial components
 nps_min=1;nps_max=1
 if (optdriver==RUNL_GSTATE) then
   nps_min=max(1,dtset%npspinor)
   nps_max=min(nproc,dtset%nspinor)
   if (tread(5)==1) nps_max=dtset%npspinor
 end if

!KGB Parallelization

!>> FFT level
 npf_min=1;npf_max=1
 npb_min=1;npb_max=1
 bpp_min=1;bpp_max=1
 n2=0;n3=0
 if (dtset%optdriver==RUNL_GSTATE) then
   npf_min=max(1,dtset%npfft)
   npf_min=min(npf_min,ngmin(2))
   npf_max=min(nproc,NPFMAX)
   if (tread(6)==1) then
     npf_max=dtset%npfft
     if (npf_max>ngmin(2)) then
       write(message,'(3a)') &
&       "Value of npfft given in input file is too high for the FFT grid!",ch10,&
&       "Action: decrease npfft or increase FFT grid (ecut, ngfft, ...)."
       MSG_ERROR(message)
     end if
   end if
   npf_max=min(npf_max,ngmin(2))
   if (dtset%use_gpu_cuda==1) then
     npf_min=1;npf_max=1
   end if

!  Number of FFT procs has to be a multiple of FFT grid sizes
!  In case of a restart from a density file, it has to be
!  compatible with the FFT grid used for the density
   n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
   if (n2==0.and.n3==0.and.(dtset%getden/=0.or.dtset%irdden/=0.or.dtset%iscf<0)) then
     dtset_found=.false.;file_found=.false.
     !1-Try to find ngfft from previous dataset
     if (dtset%getden/=0) then
       do ii=1,ndtset_alloc
         jj=dtset%getden;if (jj<0) jj=dtset%jdtset+jj
         if (dtsets(ii)%jdtset==jj) then
           dtset_found=.true.
!          n2=dtsets(ii)%nfftdg;n3=0
           n2=dtsets(ii)%ngfftdg(2);n3=dtsets(ii)%ngfftdg(3)
         end if
       end do
     end if
     !2-If not found, try to extract ngfft from density file
     if (.not.dtset_found) then
       !Retrieve file name
       suffix='_DEN';if (dtset%nimage>1) suffix='_IMG1_DEN'
       ABI_ALLOCATE(jdtset_,(0:ndtset_alloc))
       jdtset_=0;if(ndtset_alloc/=0) jdtset_(0:ndtset_alloc)=dtsets(0:ndtset_alloc)%jdtset
       call mkfilename(filnam,filden,dtset%getden,idtset,dtset%irdden,jdtset_,ndtset_alloc,suffix,'den',ii)
       ABI_DEALLOCATE(jdtset_)
       !Retrieve ngfft from file header
       idum3=0
       if (mpi_enreg%me==0) then
         inquire(file=trim(filden),exist=file_found)
         if (file_found) then
           call hdr_read_from_fname(hdr0,filden,ii,xmpi_comm_self)
           idum3(1:2)=hdr0%ngfft(2:3);if (file_found) idum3(3)=1
           call hdr_free(hdr0)
           MSG_WARNING("Cannot find filden"//filden)
         end if
       end if
       call xmpi_bcast(idum3,0,mpi_enreg%comm_world,ii)
       n2=idum3(1);n3=idum3(2);file_found=(idum3(3)/=0)
     end if
   end if

!  >> Band level
   npb_min=max(1,dtset%npband)
   npb_max=min(nproc,mband)
   if (tread(7)==1) npb_max=dtset%npband

!  >> banddp level
   bpp_min=max(1,dtset%bandpp)
   bpp_max=max(4,nint(mband/10.)) ! reasonnable bandpp max
   if (tread(8)==1) bpp_max=dtset%bandpp
 end if

!Disable KGB parallelisation in some cases:
!  - no GS
!  - paral_kgb=0 present in input file
!  - nstep=0
!  - Self-consistent DMFT
!  - Hartree-Fock or hybrid calculation (for now on)
 if ( (optdriver/=RUNL_GSTATE) .or. (dtset%paral_kgb==0.and.tread(1)==1) .or. &
& (dtset%nstep==0).or. (dtset%usedmft==1.and.dtset%nstep>1) .or. &
& (dtset%usefock==1) ) then
   nps_min=1; nps_max=1
   npf_min=1; npf_max=1
   npb_min=1; npb_max=1
   bpp_min=1; bpp_max=1
 end if

!Print title
 if (iam_master) then
   if (optdriver==RUNL_GSTATE) then
     write(message, '(8(a12,a1),a,8(i4,a4,i4,a1))' )  &
     'npimage','|','npkpt','|','npspinor','|','npfft','|','npband','|',' bandpp ' ,'|','nproc','|','weight','|', ch10, &
     npi_min,' -> ',npi_max,'|',npk_min,' -> ',npk_max,'|',nps_min,' -> ',nps_max,'|', &
     npf_min,' -> ',npf_max,'|',npb_min,' -> ',npb_max,'|',bpp_min,' -> ',bpp_max,'|', &
     nprocmin,' -> ',nproc,'|', 1 ,' -> ',nproc,'|'
   end if
   if (optdriver==RUNL_RESPFN) then
     write(message, '(4(a12,a1),a,4(i4,a4,i4,a1))' )  &
     'nppert','|','npkpt','|','nproc','|','weight','|', ch10, &
     npp_min,' -> ',npp_max,'|',      npk_min,' -> ',npk_max,'|', &
     nprocmin,' -> ',nproc,'|', 1 ,' -> ',nproc,'|'
   end if
   call wrtout(std_out,message,'COLL')
   if(max_ncpus>0) then
     call wrtout(ab_out,message,'COLL')
   end if
 end if

!Allocate lists
 ABI_ALLOCATE(my_distp,(10,MAXCOUNT))
 ABI_ALLOCATE(weight,(MAXCOUNT))
 my_distp(1:7,:)=0;weight(:)=zero
 my_distp(8,:)=dtset%use_slk
 my_distp(9,:)=dtset%np_slk
 my_distp(10,:)=dtset%gpu_linalg_limit
 icount=0;imin=1

 npc_min=1;npc_max=1;ncell_eff=1
 if (optdriver==RUNL_GSTATE) then
   ncell_eff=nimage_eff;npc_min=npi_min;npc_max=npi_max
 end if
 if (optdriver==RUNL_RESPFN) then
   ncell_eff=npert_eff;npc_min=npp_min;npc_max=npp_max
 end if

!Loop over all possibilities
!Computation of weight~"estimated acceleration"
 if (new_version) then

!  ======= NEW VERSION ========
   do npc=npc_min,npc_max
     acc_c=one;if (npc>1) acc_c=0.99_dp*speedup_fdp(ncell_eff,npc)

     do npk=npk_min,npk_max
!      -> for DFPT runs, impose that nsppol divide npk
       if (optdriver==RUNL_RESPFN.and.modulo(npk,dtset%nsppol)>0.and.npk>1) cycle
       acc_k=one;if (npk>1) acc_k=0.96_dp*speedup_fdp(nkpt_eff,npk)

       do nps=nps_min,nps_max
         acc_s=one;if (nps>1) acc_s=0.85_dp*speedup_fdp(dtset%nspinor,nps)

         do npf=npf_min,npf_max
!          -> npf should divide ngfft if set (if unset, ngfft=0 so the modulo test is ok)
           if((modulo(n2,npf)>0).or.(modulo(n3,npf)>0)) cycle
!          -> npf should be only divisible by 2, 3 or 5
           ii=npf
           do while (modulo(ii,2)==0)
             ii=ii/2
           end do
           do while (modulo(ii,3)==0)
             ii=ii/3
           end do
           do while (modulo(ii,5)==0)
             ii=ii/5
           end do
           if(ii/=1) cycle

           do npb=npb_min,npb_max
             nproc1=npc*npk*nps*npf*npb
             if (nproc1<nprocmin)     cycle
             if (nproc1>nproc)        cycle
             if (modulo(mband,npb)>0) cycle

!            Base speedup
             acc_kgb_0=one;if (npb*npf>1) acc_kgb_0=0.7_dp*speedup_fdp(mpw,(npb*npf))

             if (npb*npf>4) then
!              Promote npb=npf
               acc_kgb_0=acc_kgb_0*min((one*npf)/(one*npb),(one*npb)/(one*npf))
!              Promote npf<=20
               if (npf>20)then
                 acc_kgb_0=acc_kgb_0* &
&                 0.2_dp+(one-0.2_dp)*(sin((pi*(npf-NPF_CUTOFF))/(one*(NPFMAX-NPF_CUTOFF))) &
&                 /((pi*(npf-NPF_CUTOFF))/(one*(NPFMAX-NPF_CUTOFF))))**2
               end if
             end if

             first_bpp=.true.
             do bpp=bpp_min,bpp_max
               if (modulo(mband/npb,bpp)>0) cycle
               if ((bpp>1).and.(modulo(bpp,2)>0)) cycle
               if (one*npb*bpp >max(1.,mband/3.).and.(mband>30)) cycle
               if (npb*npf<=4.and.(.not.first_bpp)) cycle
               first_bpp=.false.

               acc_kgb=acc_kgb_0
!              Promote bpp*npb>mband/3
               if (npb*npf>4.and.mband>30) acc_kgb=acc_kgb*(one-(three*bpp*npb)/(one*mband))

!              Resulting speedup
!              weight0=acc_c*acc_k*acc_s*acc_kgb
               weight0=nproc1*(acc_c+acc_k+acc_s+acc_kgb)/(npc+npk+nps+(npf*npb))

!              Store data
               icount=icount+1
               if (icount<=MAXCOUNT) then
                 my_distp(1:7,icount)=(/npc,npk,nps,npf,npb,bpp,nproc1/)
                 weight(icount)=weight0
                 if (weight0<weight(imin)) imin=icount
               else
                 if (weight0>weight(imin)) then
                   my_distp(1:7,imin)=(/npc,npk,nps,npf,npb,bpp,nproc1/)
                   weight(imin)=weight0
                   idum=minloc(weight);imin=idum(1)
                 end if
               end if

             end do ! bpp
           end do ! npb
         end do ! npf
       end do ! nps
     end do ! npk
   end do ! npc
 else

!  ======= OLD VERSION ========
   do npc=npc_min,npc_max
     acc_c=one;if (npc>1) acc_c = 0.99_dp*ncell_eff/((ncell_eff+npc-1)/npc)

     do npk=npk_min,npk_max
       acc_k=one;if (npk>1) acc_k = 0.96_dp*nkpt_eff/((nkpt_eff+npk-1)/npk)

       do nps=nps_min,nps_max
         acc_s=one;if (nps>1) acc_s = 0.85_dp*dtset%nspinor/ ((dtset%nspinor+nps-1)/nps)

         do npf=npf_min,npf_max
!          -> npf should divide ngfft if set (if unset, ngfft=0 so the modulo test is ok)
           if((modulo(n2,npf)>0).or.(modulo(n3,npf)>0)) cycle
!          -> npf should be only divisible by 2, 3, 5, 7 or 11
           npb=npf ! Note that here, npb is used as a temp var
           do while (modulo(npb,2)==0)
             npb=npb/2
           end do
           do while (modulo(npb,3)==0)
             npb=npb/3
           end do
           do while (modulo(npb,5)==0)
             npb=npb/5
           end do
           do while (modulo(npb,7)==0)
             npb=npb/7
           end do
           do while (modulo(npb,11)==0)
             npb=npb/11
           end do
           if(npb/=1) cycle

           do npb=npb_min,npb_max
             nproc1=npc*npk*nps*npf*npb
             if (nproc1<nprocmin) cycle
             if (nproc1>nproc) cycle
             if(modulo(mband,npb)>0) cycle

             do bpp=bpp_max,bpp_min,-1
               if(modulo(mband/npb,bpp)>0) cycle
               if((bpp>1).and.(modulo(bpp,2)>0)) cycle
               if (1.*npb*bpp >max(1.,mband/3.)) cycle

               acc_kgb=one
               if (npb*npf>4) then
                 acc_kgb=min((one*npf)/(one*npb),(one*npb)/(one*npf))  * &
                 (mpw/(mpw/(npb*npf)))*(one-(three*bpp*npb)/mband)
               else if (npb*npf >1) then
                 acc_kgb=(mpw*mband/(mband*mpw/(npb*npf)))*0.7_dp
               end if

!              Weight average for efficiency and estimated acceleration
               weight0=(acc_c+acc_k+acc_s+acc_kgb)/(npc+npk+nps+(npf*npb))
               weight0=weight0*nproc1

!              Store data
               icount=icount+1
               if (icount<=MAXCOUNT) then
                 my_distp(1:7,icount)=(/npc,npk,nps,npf,npb,bpp,nproc1/)
                 weight(icount)=weight0
                 if (weight0<weight(imin)) imin=icount
               else
                 if (weight0>weight(imin)) then
                   my_distp(1:7,imin)=(/npc,npk,nps,npf,npb,bpp,nproc1/)
                   weight(imin)=weight0
                   idum=minloc(weight);imin=idum(1)
                 end if
               end if

             end do ! bpp
           end do ! npb
         end do ! npf
       end do ! nps
     end do ! npk
   end do ! npc

!  New or old version
 end if

 mcount_eff=icount
 mcount=min(mcount_eff,MAXCOUNT)

 if (mcount==0) then
   write(message,'(a,i0,2a,i0,a)')  &
   'Your input dataset does not let Abinit find an appropriate process distribution with nproc=',nproc,ch10, &
   'Try to comment all the np* vars and set paral_kgb=',-nproc,' to have advices on process distribution.'
   MSG_WARNING(message)
!  Override here the 0 default value changed in indefo1
   dtset%npimage  = max(1,dtset%npimage)
   dtset%nppert   = max(1,dtset%nppert)
   dtset%npkpt    = max(1,dtset%npkpt)
   dtset%npspinor = max(1,dtset%npspinor)
   dtset%npfft    = max(1,dtset%npfft)
   dtset%npband   = max(1,dtset%npband)
   dtset%bandpp   = max(1,dtset%bandpp)
   ABI_DEALLOCATE(my_distp)
   ABI_DEALLOCATE(weight)
   return
 end if

!* HF or hybrid calculation: no use of the fonction "autoparal"
 if ((dtset%usefock==1).AND.(dtset%nphf/=1)) then
   write(message,'(a,i5,2a,i6,a)')  &
   'Hartree-Fock or hybrid calculation : Your input dataset does not let Abinit find an appropriate process distribution.'
   MSG_WARNING(message)
!  Override here the 0 default value changed in indefo1
   dtset%npimage  = max(1,dtset%npimage)
   dtset%npkpt    = max(1,dtset%npkpt)
   dtset%npspinor = max(1,dtset%npspinor)
   dtset%npfft    = max(1,dtset%npfft)
   dtset%npband   = max(1,dtset%npband)
   dtset%bandpp   = max(1,dtset%bandpp)
   ABI_DEALLOCATE(my_distp)
   ABI_DEALLOCATE(weight)
   return
 end if

!Sort data by increasing weight
 ABI_ALLOCATE(isort,(mcount))
 isort=(/(ii,ii=1,mcount)/)
 call sort_dp(mcount,weight,isort,tol6)

 ncount=mcount;if (dtset%paral_kgb>=0) ncount=min(mcount,5)
 if (iam_master) then
   do jj=mcount,mcount-ncount+1,-1
     ii=isort(jj)
     if (optdriver==RUNL_GSTATE) then
       write(message, '(7(i12,a1),f11.2,a2)') &
&       my_distp(1,ii),'|',my_distp(2,ii),'|',my_distp(3,ii),'|',my_distp(4,ii),'|', &
&       my_distp(5,ii),'|',my_distp(6,ii),'|',my_distp(7,ii),'|',weight(jj),' |'
     end if
     if (optdriver==RUNL_RESPFN) then
       write(message, '(3(i12,a1),f11.2,a2)') &
&       my_distp(1,ii),'|',my_distp(2,ii),'|',my_distp(7,ii),'|',weight(jj),' |'
     end if
     call wrtout(std_out,message,'COLL')
     if(max_ncpus>0) then
       call wrtout(ab_out,message,'COLL')
     end if
   end do
 end if

 if (max_ncpus>0.and.(mcount_eff>MAXCOUNT)) then
   write(message,'(a,i0,a,i0,a)') &
&   ' Received max_ncpus ',max_ncpus,' possible choices for nproc; only the first ',MAXCOUNT,' ones are printed...'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

 !if (iam_master .and. dtset%paral_kgb<0) then
 if (iam_master .and. max_ncpus>0) then
   write(ount,'(2a)')ch10,"--- !Autoparal"

   if (optdriver==RUNL_GSTATE) then
     write(ount,"(a)")"#Autoparal section for GS calculations with paral_kgb"
   else if (optdriver==RUNL_RESPFN) then
     write(ount,"(a)")'#Autoparal section for DFPT calculations'
   else
     MSG_ERROR("Unsupported optdriver")
   end if

   write(ount,"(a)")   "info:"
   write(ount,"(a,i0)")"    autoparal: ",autoparal
   write(ount,"(a,i0)")"    paral_kgb: ",dtset%paral_kgb
   write(ount,"(a,i0)")"    max_ncpus: ",max_ncpus
   write(ount,"(a,i0)")"    nspinor: ",dtset%nspinor
   write(ount,"(a,i0)")"    nsppol: ",dtset%nsppol
   write(ount,"(a,i0)")"    nkpt: ",dtset%nkpt
   write(ount,"(a,i0)")"    mband: ",mband

   ! List of configurations.
   write(ount,"(a)")"configurations:"

   if (optdriver==RUNL_GSTATE) then

     do jj=mcount,mcount-ncount+1,-1
       ii=isort(jj)
       tot_ncpus = my_distp(7,ii)
       eff = weight(jj) / tot_ncpus

       write(ount,"(a,i0)")"    - tot_ncpus: ",tot_ncpus
       write(ount,"(a,i0)")"      mpi_ncpus: ",tot_ncpus
       !write(ount,"(a,i0)")"      omp_ncpus: ",omp_ncpus !OMP not supported  (yet)
       write(ount,"(a,f12.9)")"      efficiency: ",eff
       !write(ount,"(a,f12.2)")"      mem_per_cpu: ",mempercpu_mb

       ! list of variables to use.
       !'npimage','|','npkpt','|','npspinor','|','npfft','|','npband','|',' bandpp ' ,'|','nproc','|','weight','|'
       write(ount,"(a)"   )"      vars: {"
       write(ount,"(a,i0,a)")"            npimage: ",my_distp(1,ii),","
       write(ount,"(a,i0,a)")"            npkpt: ", my_distp(2,ii),","
       write(ount,"(a,i0,a)")"            npspinor: ",my_distp(3,ii),","
       write(ount,"(a,i0,a)")"            npfft: ", my_distp(4,ii),","
       write(ount,"(a,i0,a)")"            npband: ",my_distp(5,ii),","
       write(ount,"(a,i0,a)")"            bandpp: ",my_distp(6,ii),","
       write(ount,"(a)")   "            }"
     end do

   else if (optdriver==RUNL_RESPFN) then

     do jj=mcount,mcount-ncount+1,-1
       ii=isort(jj)
       tot_ncpus = my_distp(7,ii)
       eff = weight(jj) / tot_ncpus

       write(ount,"(a,i0)")"    - tot_ncpus: ",tot_ncpus
       write(ount,"(a,i0)")"      mpi_ncpus: ",tot_ncpus
       !write(ount,"(a,i0)")"      omp_ncpus: ",omp_ncpus !OMP not supported  (yet)
       write(ount,"(a,f12.9)")"      efficiency: ",eff
       !write(ount,"(a,f12.2)")"      mem_per_cpu: ",mempercpu_mb
       ! list of variables to use.
       !'nppert','|','npkpt','|','nproc','|','weight','|',
       write(ount,"(a)"   )"      vars: {"
       write(ount,"(a,i0,a)")"             nppert: ", my_distp(1,ii),","
       write(ount,"(a,i0,a)")"             npkpt: ", my_distp(2,ii),","
       write(ount,"(a)")   "            }"
     end do

   end if
   write(ount,'(a)')"..."
   iexit = iexit + 1
 end if

 icount=isort(mcount)

!Refinement of the process distribution by mean of a LinAlg routines benchmarking
 if (optdriver==RUNL_GSTATE.and.autoparal/=1) then
   if (autoparal/=3) then
     if (autoparal==2) then
       write(message,'(5a,9(a10,a1))') ch10, &
&       ' Values below have been tested with respect to Linear Algebra performance;',ch10,&
&       ' Weights below are corrected according:',ch10,&
&       'npimage','|','npkpt' ,'|','npspinor'  ,'|','npfft'     ,'|','npband','|',' bandpp ' ,'|',&
&       'nproc'  ,'|','weight','|','new weight','|'
     else
       write(message,'(5a,11(a10,a1))') ch10, &
&       ' Values below have been tested with respect to Linear Algebra performance;',ch10,&
&       ' Weights below are corrected according:',ch10,&
&       'npimage','|','npkpt' ,'|','npspinor'  ,'|','npfft'     ,'|','npband','|',' bandpp ' ,'|',&
&       'nproc'  ,'|','weight','|','new weight','|','best npslk','|','linalggpu' ,'|'
     end if
     call wrtout(std_out,message,'COLL')
     if (max_ncpus > 0) then
       call wrtout(ab_out,message,'COLL')
     end if
   end if
   acc_k=zero
   ncount=min(MAXBENCH,mcount);if (autoparal==3) ncount=1
   do jj=mcount,mcount-ncount+1,-1
     ii=isort(jj)
     npf=my_distp(4,ii);npb=my_distp(5,ii);bpp=my_distp(6,ii)
     if ((npb*npf*bpp>1).and.(npf*npb<=mpi_enreg%nproc)) then
       use_linalg_gpu=dtset%use_gpu_cuda
       call compute_kgb_indicator(acc_kgb,bpp,xmpi_world,mband,mpw,npb,npf,np_slk,use_linalg_gpu)
       if (autoparal/=2) then
         my_distp(9,ii)=np_slk
         if (np_slk>0) my_distp(8,ii)=1
!        * gpu_linalg_limit:
!        No use of GPU: huge value ~2  *vectsize*blocksize**2 tested
!        Use of GPU:    tiny value ~0.5*vectsize*blocksize**2 tested
         my_distp(10,ii)=2*dtset%mpw*(npb*bpp)**2/npf
         if (use_linalg_gpu==1) my_distp(10,ii)=my_distp(10,ii)/4
       end if
       if (abs(acc_k)<=tol12) acc_k=acc_kgb ! Ref value : the first one computed
!      * Weight (corrected by 10% of the computed ratio)
       weight0=weight(jj)*(one + 0.1_dp*acc_k/acc_kgb)
       if (autoparal==2) then
         write(message, '(7(i10,a1),f9.2,a2,f9.5,a2)') &
&         my_distp(1,ii),'|',my_distp(2,ii),'|',my_distp(3,ii),'|',my_distp(4,ii),'|',&
&         my_distp(5,ii),'|',my_distp(6,ii),'|',my_distp(7,ii),'|',weight(jj),'=>', weight0,' |'
       else if (autoparal==3) then
         write(message,'(a,5(a,i3))') ch10,' For npband=',npb,', npfft=',npf,' and bandpp=',bpp, &
&         ', compute_kgb_indicator recommends you to set np_slk=',my_distp(9,ii),&
&         ' and use_linalg_gpu=',use_linalg_gpu
       else
         write(message, '(7(i10,a1),f9.2,a2,f9.5,a2,2(i10,a1))') &
&         my_distp(1,ii),'|',my_distp(2,ii),'|',my_distp(3,ii),'|',my_distp(4,ii),'|',&
&         my_distp(5,ii),'|',my_distp(6,ii),'|',my_distp(7,ii),'|',weight(jj),'=>', weight0,' |',&
&         my_distp(9,ii),'|',use_linalg_gpu,'|'
       end if
       call wrtout(std_out,message,'COLL')
       if (max_ncpus>0) then
         call wrtout(ab_out,message,'COLL')
       end if
!      We store the best value in weight(mcount) and keep icount
       if (weight0 > weight(mcount)) then
         icount=ii;weight(mcount)=weight0
       end if
     end if
   end do
 end if

!Final advice in case max_ncpus > 0
 if (max_ncpus>0) then
   write(message,'(6a)') ch10,&
&   ' Launch a parallel version of ABINIT with a number of processors among the above list,',ch10,&
&   ' and the associated input variables npkpt, npband, npfft and bandpp. ',ch10,&
&   ' The optimal weight is close to nproc and the higher should be better.'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   iexit=iexit+1
   GOTO 100
 end if

!Store new process distribution
 if (dtset%paral_kgb>=0) then
   nproc1=my_distp(7,icount)
!  Work load distribution
   if (optdriver==RUNL_GSTATE) then
     dtset%npimage= my_distp(1,icount)
     dtset%nppert = 1
     dtset%npkpt  = my_distp(2,icount)
   end if
   if (optdriver==RUNL_RESPFN) then
     dtset%npimage= 1
     dtset%nppert = my_distp(1,icount)
     dtset%npkpt  = 1
   end if
   dtset%npspinor = my_distp(3,icount)
   dtset%npfft    = my_distp(4,icount)
   dtset%npband   = my_distp(5,icount)
   dtset%bandpp   = my_distp(6,icount)
!  The following lines are mandatory : the DFT+DMFT must use ALL the
!  available procs specified by the user. So nproc1=nproc.
!  Works only if paral_kgb is not activated.
   if (dtset%usedmft/=0.and.optdriver==RUNL_GSTATE.and.dtset%paral_kgb==0) then
     dtset%npspinor = 1
     dtset%npfft    = 1
     dtset%npband   = 1
     dtset%bandpp   = 1
     dtset%npimage  = 1
     nproc1         = nproc
   end if
   if (dtset%npband*dtset%npfft*dtset%bandpp>1) dtset%paral_kgb=1
!  LinAlg parameters: we change values only if they are not present in input file
   if (dtset%paral_kgb==1) then
     if (tread(9)==0) dtset%use_slk=my_distp(8,icount)
     if (tread(10)==0) dtset%np_slk=my_distp(9,icount)
     if (tread(11)==0) dtset%gpu_linalg_limit=my_distp(10,icount)
   end if
!  New definition of "world" MPI communicator
   if (optdriver==RUNL_RESPFN.and.dtset%paral_rf==1) then
     nproc1=max(nproc1,dtset%nsppol*dtset%nkpt) ! Take into account the code in respfn.F90
     nproc1=min(nproc1,nproc)
     nproc1=(nproc1/dtset%nppert)*dtset%nppert
   end if
   call initmpi_world(mpi_enreg,nproc1)
!   call initmpi_world(mpi_enreg,nproc1)
 end if

 100 continue

 ABI_DEALLOCATE(isort)
 ABI_DEALLOCATE(my_distp)
 ABI_DEALLOCATE(weight)

 DBG_EXIT("COLL")

 contains

   function speedup_fdp(nn,mm)
   !Expected linear speedup for a nn-sized problem and mm processes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'speedup_fdp'
!End of the abilint section

   real(dp) :: speedup_fdp
   integer,intent(in) :: nn,mm
   speedup_fdp=(one*nn)/(one*((nn/mm)+merge(0,1,mod(nn,mm)==0)))
 end function speedup_fdp

end subroutine finddistrproc
!!***
