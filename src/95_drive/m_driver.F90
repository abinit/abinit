!!****m* ABINIT/m_driver
!! NAME
!!  m_driver
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group ()
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

module m_driver

 use defs_basis
 use defs_wvltypes
 use m_errors
 use m_dtset
 use m_results_out
 use m_results_respfn
 use m_yaml
 use m_xmpi
 use m_xomp
 use m_abi_linalg
 use m_abicore
 use m_exit
 use m_dtfil
 use m_fftcore
 use libxc_functionals
#if defined DEV_YP_VDWXC
 use m_xc_vdw
#endif
 use m_xgScalapack

 use defs_datatypes, only : pseudopotential_type, pspheader_type
 use defs_abitypes,  only : MPI_type
 use m_fstrings,     only : sjoin, itoa
 use m_time,         only : timab
 use m_xg,           only : xg_finalize
 use m_libpaw_tools, only : libpaw_write_comm_set
 use m_geometry,     only : mkrdim, xcart2xred, xred2xcart, chkdilatmx
 use m_pawang,       only : pawang_type, pawang_free
 use m_pawrad,       only : pawrad_type, pawrad_free
 use m_pawtab,       only : pawtab_type, pawtab_nullify, pawtab_free
 use m_fftw3,        only : fftw3_init_threads, fftw3_cleanup
 use m_psps,         only : psps_init_global, psps_init_from_dtset, psps_free
 use m_mpinfo,       only : mpi_distrib_is_ok
 use m_respfn_driver,    only : respfn
 use m_screening_driver, only : screening
 use m_sigma_driver,     only : sigma
 use m_bethe_salpeter,   only : bethe_salpeter
 use m_eph_driver,       only : eph
 use m_wfk_analyze,      only : wfk_analyze
 use m_gstateimg,        only : gstateimg
 use m_gwls_sternheimer, only : gwls_sternheimer
 use m_nonlinear,        only : nonlinear
 use m_drivexc,          only : echo_xc_name

#if defined HAVE_BIGDFT
 use BigDFT_API,   only: xc_init, xc_end, XC_MIXED, XC_ABINIT,&
&                        mpi_environment_set,bigdft_mpi, f_malloc_set_status
#endif

 use m_longwave

 implicit none

 private
!!***

 public :: driver
!!***

contains
!!***

!!****f* ABINIT/driver
!! NAME
!! driver
!!
!! FUNCTION
!! Driver for ground state, response function, screening
!! and sigma calculations. The present routine drives the following operations.
!! An outer loop allows computation related to different data sets.
!! For each data set, either a GS calculation, a RF calculation,
!! a SUS calculation, a SCR calculation or a SIGMA calculation is made.
!! In both cases, the input variables are transferred in the proper variables,
!! selected big arrays are allocated, then the gstate, respfn, ...  subroutines are called.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2020 ABINIT group (XG,MKV,MM,MT,FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! codvsn= code version
!! cpui=initial CPU time
!! filnam(5)=character strings giving file names
!! filstat=character strings giving name of status file
!! mpi_enregs=information about MPI parallelization
!! ndtset=number of datasets
!! ndtset_alloc=number of datasets, corrected for allocation of at least one data set.
!! npsp=number of pseudopotentials
!! pspheads(npsp)=<type pspheader_type>all the important informatio from the
!!   pseudopotential file header, as well as the psp file name
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  Input/Output
!! results_out(0:ndtset_alloc)=<type results_out_type>contains the results
!!   needed for outvars, including evolving variables
!!   Default values are set up in the calling routine
!! dtsets(0:ndtset_alloc)=<type datasets_type>
!!   input: all input variables initialized from the input file.
!!   output: the effective set of variables used in the different datasets.
!!           Some variables, indeed, might have been redefined in one of the children.
!!           of this routine.
!!
!! NOTES
!! The array filnam is used for the name of input and output files,
!! and roots for generic input, output or temporary files.
!! Pseudopotential file names are set in pspini and pspatm,
!! using another name. The name filstat will be needed beyond gstate to check
!! the appearance of the "exit" flag, to make a hasty exit, as well as
!! in order to output the status of the computation.
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      abi_io_redirect,abi_linalg_finalize,abi_linalg_init,bethe_salpeter
!!      chkdilatmx,destroy_results_out,destroy_results_respfn,dtfil_init
!!      dtfil_init_img,dtset%free,dtsets,echo_xc_name,eph,exit_check
!!      f_malloc_set_status,fftw3_cleanup,fftw3_init_threads,find_getdtset
!!      gather_results_out,gstateimg,gwls_sternheimer,init_results_respfn
!!      libpaw_write_comm_set,libxc_functionals_end,libxc_functionals_init
!!      longwave,mkrdim,mpi_environment_set,nonlinear,nullify_wvl_data
!!      pawang_free,pawrad_free,pawtab_free,pawtab_nullify,psps_free
!!      psps_init_from_dtset,psps_init_global,respfn,screening,sigma,timab
!!      wfk_analyze,wrtout,xc_vdw_done,xc_vdw_init,xc_vdw_libxc_init
!!      xc_vdw_memcheck,xc_vdw_read,xc_vdw_show,xc_vdw_trigger,xc_vdw_write
!!      xcart2xred,xg_finalize,xgscalapack_config,xmpi_bcast,xred2xcart
!!      yaml_iterstart,ydoc%add_ints,ydoc%add_reals,ydoc%write_and_free
!!
!! SOURCE

subroutine driver(codvsn,cpui,dtsets,filnam,filstat,&
&                 mpi_enregs,ndtset,ndtset_alloc,npsp,pspheads,results_out)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: ndtset,ndtset_alloc,npsp
 real(dp),intent(in) :: cpui
 character(len=8),intent(in) :: codvsn
 character(len=fnlen),intent(in) :: filstat
 type(MPI_type),intent(inout) :: mpi_enregs(0:ndtset_alloc)
 !arrays
 character(len=fnlen),intent(in) :: filnam(5)
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
 type(pspheader_type),intent(in) :: pspheads(npsp)
 type(results_out_type),target,intent(inout) :: results_out(0:ndtset_alloc)

 !Local variables-------------------------------
 !scalars
#if defined DEV_YP_VDWXC
 character(len=fnlen) :: vdw_filnam
#endif
 integer,parameter :: level=2,mdtset=9999
 integer,save :: paw_size_old=-1
 integer :: idtset,ierr,iexit,iget_cell,iget_occ,iget_vel,iget_xcart,iget_xred
 integer :: ii,iimage,iimage_get,jdtset,jdtset_status,jj,kk,linalg_max_size
 integer :: mtypalch,mu,mxnimage,nimage,openexit,paw_size,prtvol, omp_nthreads
 real(dp) :: etotal
 character(len=500) :: msg, dilatmx_errmsg
 logical :: converged,results_gathered,test_img,use_results_all
 type(dataset_type) :: dtset
 type(datafiles_type) :: dtfil
 type(pawang_type) :: pawang
 type(pseudopotential_type) :: psps
 type(results_respfn_type) :: results_respfn
 type(wvl_data) :: wvl
 type(yamldoc_t) :: ydoc
#if defined DEV_YP_VDWXC
 type(xc_vdw_type) :: vdw_params
#endif
 !arrays
 integer :: mkmems(3)
 integer,allocatable :: jdtset_(:),npwtot(:)
 real(dp) :: acell(3),rprim(3,3),rprimd(3,3),tsec(2)
 real(dp),allocatable :: acell_img(:,:),amu_img(:,:),rprim_img(:,:,:)
 real(dp),allocatable :: fcart_img(:,:,:),fred_img(:,:,:),intgres_img(:,:,:)
 real(dp),allocatable :: etotal_img(:),mixalch_img(:,:,:),strten_img(:,:),miximage(:,:)
 real(dp),allocatable :: occ(:),xcart(:,:),xred(:,:),xredget(:,:)
 real(dp),allocatable :: occ_img(:,:),vel_cell_img(:,:,:),vel_img(:,:,:),xred_img(:,:,:)
 type(pawrad_type),allocatable :: pawrad(:)
 type(pawtab_type),allocatable :: pawtab(:)
 type(results_out_type),pointer :: results_out_all(:)

!******************************************************************

 DBG_ENTER("COLL")

 call timab(640,1,tsec)
 call timab(641,3,tsec)

!Structured debugging if prtvol==-level
 prtvol=dtsets(1)%prtvol
 if(prtvol==-level)then
   write(msg,'(80a,a,a)')  ('=',ii=1,80),ch10,' driver : enter , debug mode '
   call wrtout(std_out, msg)
 end if

 if(ndtset>mdtset)then
   write(msg,'(a,i0,a,i0,a)')'  The maximal allowed ndtset is ',mdtset,' while the input value is ',ndtset,'.'
   MSG_BUG(msg)
 end if

 mtypalch=dtsets(1)%ntypalch
 do ii=1,ndtset_alloc
   mtypalch=max(dtsets(ii)%ntypalch,mtypalch)
 end do
 call psps_init_global(mtypalch, npsp, psps, pspheads)

 ABI_ALLOCATE(jdtset_,(0:ndtset))
 if(ndtset/=0)then
   jdtset_(:)=dtsets(0:ndtset)%jdtset
 else
   jdtset_(0)=0
 end if

 do idtset=1,ndtset_alloc
   call init_results_respfn(dtsets,ndtset_alloc,results_respfn)
 end do

 call timab(641,2,tsec)

!*********************************************************************
!Big loop on datasets

!Do loop on idtset (allocate statements are present)
 do idtset=1,ndtset_alloc

   if(mpi_enregs(idtset)%me<0) cycle

   call timab(642,1,tsec)
   call abi_io_redirect(new_io_comm=mpi_enregs(idtset)%comm_world)
   call libpaw_write_comm_set(mpi_enregs(idtset)%comm_world)

   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=1

   if(ndtset>=2)then
     jdtset_status=jdtset
   else
     jdtset_status=0
   end if

!  Copy input variables into a local dtset.
   dtset = dtsets(idtset)%copy()

!  Set other values
   dtset%jdtset = jdtset
   dtset%ndtset = ndtset

!  Print DATASET number
   write(msg,'(83a)') ch10,('=',mu=1,80),ch10,'== DATASET'
   if (jdtset>=100) then
     write(msg,'(2a,i4,65a)') trim(msg),' ',jdtset,' ',('=',mu=1,64)
   else
     write(msg,'(2a,i2,67a)') trim(msg),' ',jdtset,' ',('=',mu=1,66)
   end if

#ifdef HAVE_OPENMP
   omp_nthreads = xomp_get_num_threads(open_parallel=.True.)
#else
   omp_nthreads = -1
#endif
   write(msg,'(2a,2(a, i0), a)') trim(msg),ch10,&
     '-   mpi_nproc: ',mpi_enregs(idtset)%nproc, ", omp_nthreads: ", omp_nthreads, " (-1 if OMP is not activated)"

   if (dtset%optdriver == RUNL_GSTATE) then
     if (.not. mpi_distrib_is_ok(mpi_enregs(idtset),dtset%mband,dtset%nkpt,dtset%mkmem,dtset%nsppol)) then
       write(msg,'(3a)') trim(msg),ch10,'-    --> not optimal distribution: autoparal keyword recommended in input file <--'
     end if
   end if

   write(msg,'(3a)') trim(msg),ch10,' '
   call wrtout(ab_out, msg)
   call wrtout(std_out, msg, 'PERS')     ! PERS is choosen to make debugging easier

   call yaml_iterstart('dtset', jdtset, ab_out, dtset%use_yaml)

   if (mpi_enregs(idtset)%me_cell == 0) then
     ydoc = yamldoc_open('DatasetInfo')

     ! Write basic dimensions.
     call ydoc%add_ints( &
       "natom, nkpt, mband, nsppol, nspinor, nspden, mpw", &
       [dtset%natom, dtset%nkpt, dtset%mband, dtset%nsppol, dtset%nspinor, dtset%nspden, dtset%mpw], &
       dict_key="dimensions")

     ! Write cutoff energies.
     call ydoc%add_reals("ecut, pawecutdg", [dtset%ecut, dtset%pawecutdg], &
       real_fmt="(f5.1)", dict_key="cutoff_energies")

     ! Write info on electrons.
     call ydoc%add_reals("nelect, charge, occopt, tsmear", &
       [dtset%nelect, dtset%charge, one * dtset%occopt, dtset%tsmear], &
       dict_key="electrons")

     ! This part depends on optdriver.
     ! Third-party Yaml parsers may use optdriver to specialize the logic used to interpret the next data.
     select case (dtset%optdriver)
     case (RUNL_GSTATE)
       call ydoc%add_ints("optdriver, ionmov, optcell, iscf, paral_kgb", &
         [dtset%optdriver, dtset%ionmov, dtset%optcell, dtset%iscf, dtset%paral_kgb], &
         dict_key="meta")

     case(RUNL_RESPFN)
       call ydoc%add_ints("optdriver, rfddk, rfelfd, rfmagn, rfphon, rfstrs", &
         [dtset%optdriver, dtset%rfddk, dtset%rfelfd, dtset%rfmagn, dtset%rfphon, dtset%rfstrs], &
         ignore=0, dict_key="meta")
         ! dtset%rfdir ??

     case(RUNL_NONLINEAR)
       call ydoc%add_ints("optdriver", [dtset%optdriver], &
         dict_key="meta")

     case(RUNL_GWLS)
       call ydoc%add_ints("optdriver", [dtset%optdriver], &
         dict_key="meta")

     case(RUNL_WFK)
       call ydoc%add_ints("optdriver, wfk_task", &
         [dtset%optdriver, dtset%wfk_task] , dict_key="meta")

     case (RUNL_SCREENING, RUNL_SIGMA)
       call ydoc%add_ints("optdriver, gwcalctyp", &
         [dtset%optdriver, dtset%gwcalctyp] , dict_key="meta")

     case (RUNL_BSE)
       call ydoc%add_ints("optdriver, bs_calctype, bs_algorithm", &
         [dtset%optdriver, dtset%bs_calctype, dtset%bs_algorithm] , dict_key="meta")

     case (RUNL_EPH)
       call ydoc%add_ints("optdriver, eph_task", &
         [dtset%optdriver, dtset%eph_task] , dict_key="meta")

     case(RUNL_LONGWAVE)
       call ydoc%add_ints("optdriver", [dtset%optdriver], &
         dict_key="meta")

     case default
       MSG_ERROR(sjoin('Add a meta section for optdriver: ', itoa(dtset%optdriver)))
     end select

     !if (dtset%use_yaml == 1) then
     call ydoc%write_and_free(ab_out)
     !else
     !  call ydoc%write_and_free(std_out)
     !end if
   end if

   if ( dtset%np_slk == 0 ) then
     call xgScalapack_config(SLK_DISABLED,dtset%slk_rankpp)
   else if ( dtset%np_slk == 1000000 ) then
     call xgScalapack_config(SLK_AUTO,dtset%slk_rankpp)
   else if ( dtset%np_slk > 1 ) then
     call xgScalapack_config(dtset%np_slk,dtset%slk_rankpp)
   else
     call xgScalapack_config(SLK_AUTO,dtset%slk_rankpp)
   end if

!  Copy input values
   mkmems(1) = dtset%mkmem
   mkmems(2) = dtset%mkqmem
   mkmems(3) = dtset%mk1mem

!  Initialize MPI data for the parallelism over images
   nimage=mpi_enregs(idtset)%my_nimage

!  Retrieve evolving arrays (most of them depend on images)
   ABI_ALLOCATE(acell_img,(3,nimage))
   ABI_ALLOCATE(amu_img,(dtset%ntypat,nimage))
   ABI_ALLOCATE(mixalch_img,(dtset%npspalch,dtset%ntypalch,nimage))
   ABI_ALLOCATE(occ_img,(dtset%mband*dtset%nkpt*dtset%nsppol,nimage))
   ABI_ALLOCATE(rprim_img,(3,3,nimage))
   ABI_ALLOCATE(vel_img,(3,dtset%natom,nimage))
   ABI_ALLOCATE(vel_cell_img,(3,3,nimage))
   ABI_ALLOCATE(xred_img,(3,dtset%natom,nimage))
   do iimage=1,nimage
     ii=mpi_enregs(idtset)%my_imgtab(iimage)
     acell_img   (:  ,iimage) = dtset%acell_orig(:,ii)
     amu_img   (:  ,iimage)   = dtset%amu_orig(1:dtset%ntypat,ii)
     mixalch_img   (:,:,iimage) =dtset%mixalch_orig(1:dtset%npspalch,1:dtset%ntypalch,ii)
     rprim_img   (:,:,iimage) = dtset%rprim_orig(:,:,ii)
     vel_img     (:,:,iimage) = dtset%vel_orig(:,1:dtset%natom,ii)
     vel_cell_img(:,:,iimage) = dtset%vel_cell_orig(:,:,ii)
     xred_img (:,:,iimage) = dtset%xred_orig(:,1:dtset%natom,ii)
     occ_img  (:  ,iimage) = dtset%occ_orig(1:dtset%mband*dtset%nkpt*dtset%nsppol,ii)
   end do

!  ****************************************************************************
!  Treat the file names (get variables)

!  In the case of multiple images, the file names will be overwritten later (for each image)
   call dtfil_init(dtfil,dtset,filnam,filstat,idtset,jdtset_,mpi_enregs(idtset),ndtset)
   if (dtset%optdriver==RUNL_GSTATE.and.dtset%nimage>1) then
     call dtfil_init_img(dtfil,dtset,dtsets,idtset,jdtset_,ndtset,ndtset_alloc)
   end if

!  ****************************************************************************
!  Treat other get variables

!  If multi dataset mode, and already the second dataset,
!  treatment of other get variables.
   if( ndtset>1 .and. idtset>1 )then

!    Check if parallelization over images is activated
     mxnimage=maxval(dtsets(1:ndtset_alloc)%nimage)
     ABI_ALLOCATE(miximage,(mxnimage,mxnimage))
     test_img=(mxnimage/=1.and.maxval(dtsets(:)%npimage)>1)
     use_results_all=.false.
     if (test_img.and.mpi_enregs(idtset)%me_cell==0) then
       use_results_all=.true.
       ABI_MALLOC(results_out_all,(0:ndtset_alloc))
     else
       results_out_all => results_out
     end if
     results_gathered=.false.

     call find_getdtset(dtsets,dtset%getocc,'getocc',idtset,iget_occ,miximage,mxnimage,ndtset_alloc)
     if(iget_occ/=0)then
!      Gather contributions to results_out from images, if needed
       if (test_img.and.mpi_enregs(idtset)%me_cell==0.and.(.not.results_gathered)) then
         call gather_results_out(dtsets,mpi_enregs,results_out,results_out_all,use_results_all, &
           allgather=.true.,only_one_per_img=.true.)
         results_gathered=.true.
       end if
       if ((.not.test_img).or.mpi_enregs(idtset)%me_cell==0) then
         do iimage=1,nimage
           ii=mpi_enregs(idtset)%my_imgtab(iimage)
           occ_img(:,iimage)=zero
           do iimage_get=1,dtsets(iget_occ)%nimage
             do jj=1,dtset%mband*dtset%nkpt*dtset%nsppol
               occ_img(jj,iimage)=occ_img(jj,iimage)+miximage(ii,iimage_get) &
&               *results_out_all(iget_occ)%occ(jj,iimage_get)
             end do
           end do
         end do
       end if
     end if

!    Getcell has to be treated BEFORE getxcart since acell and rprim will be used
     call find_getdtset(dtsets,dtset%getcell,'getcell',idtset,iget_cell,miximage,mxnimage,ndtset_alloc)
     if(iget_cell/=0)then
!      Gather contributions to results_out from images, if needed
       if (test_img.and.mpi_enregs(idtset)%me_cell==0.and.(.not.results_gathered)) then
         call gather_results_out(dtsets,mpi_enregs,results_out,results_out_all,use_results_all, &
&         allgather=.true.,only_one_per_img=.true.)
         results_gathered=.true.
       end if
       if ((.not.test_img).or.mpi_enregs(idtset)%me_cell==0) then
         do iimage=1,nimage
           ii=mpi_enregs(idtset)%my_imgtab(iimage)
           acell_img(:,iimage)=zero
           rprim_img(:,:,iimage)=zero
           do iimage_get=1,dtsets(iget_cell)%nimage
             do jj=1,3
               acell_img(jj,iimage)=acell_img(jj,iimage)+&
&               miximage(ii,iimage_get)*results_out_all(iget_cell)%acell(jj,iimage_get)
               do kk=1,3
                 rprim_img(kk,jj,iimage)=rprim_img(kk,jj,iimage)+&
&                 miximage(ii,iimage_get)*results_out_all(iget_cell)%rprim(kk,jj,iimage_get)
               end do
             end do
             dtset%rprimd_orig(:,:,iimage)=dtsets(iget_cell)%rprimd_orig(:,:,iimage)
!            Check that the new acell and rprim are consistent with the input dilatmx
             call mkrdim(acell_img(:,iimage),rprim_img(:,:,iimage),rprimd)
             call chkdilatmx(dtset%chkdilatmx,dtset%dilatmx,rprimd,&
&             dtset%rprimd_orig(1:3,1:3,iimage), dilatmx_errmsg)
             _IBM6("dilatxm_errmsg: "//TRIM(dilatmx_errmsg))
             if (LEN_TRIM(dilatmx_errmsg) /= 0) then
               acell_img(1,iimage) = sqrt(sum(rprimd(:,1)**2))
               acell_img(2,iimage) = sqrt(sum(rprimd(:,2)**2))
               acell_img(3,iimage) = sqrt(sum(rprimd(:,3)**2))
               rprim_img(:,1,iimage) = rprimd(:,1) / acell_img(1,iimage)
               rprim_img(:,2,iimage) = rprimd(:,2) / acell_img(1,iimage)
               rprim_img(:,3,iimage) = rprimd(:,3) / acell_img(1,iimage)
               MSG_WARNING(dilatmx_errmsg)
             end if
           end do
         end do
       end if
     end if

     call find_getdtset(dtsets,dtset%getxred,'getxred',idtset,iget_xred,miximage,mxnimage,ndtset_alloc)
     if(iget_xred/=0)then
!      Gather contributions to results_out from images, if needed
       if (test_img.and.mpi_enregs(idtset)%me_cell==0.and.(.not.results_gathered)) then
         call gather_results_out(dtsets,mpi_enregs,results_out,results_out_all,use_results_all, &
&         allgather=.true.,only_one_per_img=.true.)
         results_gathered=.true.
       end if
       if ((.not.test_img).or.mpi_enregs(idtset)%me_cell==0) then
         do iimage=1,nimage
           ii=mpi_enregs(idtset)%my_imgtab(iimage)
           xred_img(:,:,iimage)=zero
           do iimage_get=1,dtsets(iget_xred)%nimage
             do jj=1,dtset%natom
               do kk=1,3
                 xred_img(kk,jj,iimage)=xred_img(kk,jj,iimage)+&
&                 miximage(ii,iimage_get)*results_out_all(iget_xred)%xred(kk,jj,iimage_get)
               end do
             end do
           end do
         end do
       end if
     end if

     call find_getdtset(dtsets,dtset%getxcart,'getxcart',idtset,iget_xcart,miximage,mxnimage,ndtset_alloc)
     if(iget_xcart/=0)then
!      Gather contributions to results_out from images, if needed
       if (test_img.and.mpi_enregs(idtset)%me_cell==0.and.(.not.results_gathered)) then
         call gather_results_out(dtsets,mpi_enregs,results_out,results_out_all,use_results_all, &
&         allgather=.true.,only_one_per_img=.true.)
         results_gathered=.true.
       end if
       if ((.not.test_img).or.mpi_enregs(idtset)%me_cell==0) then
         ABI_ALLOCATE(xcart,(3,dtset%natom))
         ABI_ALLOCATE(xredget,(3,dtset%natom))
         do iimage=1,nimage
           ii=mpi_enregs(idtset)%my_imgtab(iimage)
           xred_img(:,:,iimage)=zero
           do iimage_get=1,dtsets(iget_xcart)%nimage
!            Compute xcart of the previous dataset
             call  mkrdim(results_out_all(iget_xcart)%acell(:,iimage_get),&
&             results_out_all(iget_xcart)%rprim(:,:,iimage_get),rprimd)
             do jj=1,dtset%natom
               do kk=1,3
                 xredget (kk,jj)=results_out_all(iget_xcart)%xred(kk,jj,iimage_get)
               end do
             end do
             call xred2xcart(dtset%natom,rprimd,xcart,xredget)
!            xcart from previous dataset is computed. Now, produce xred for the new dataset,
!            with the new acell and rprim ...
             call mkrdim(acell_img(:,iimage),rprim_img(:,:,iimage),rprimd)
             call xcart2xred(dtset%natom,rprimd,xcart,xredget(:,:))
             do jj=1,dtset%natom
               do kk=1,3
                 xred_img(kk,jj,iimage)=xred_img(kk,jj,iimage)+miximage(ii,iimage_get)*xredget(kk,jj)
               end do
             end do
           end do
         end do
         ABI_FREE(xcart)
         ABI_FREE(xredget)
       end if
     end if

     call find_getdtset(dtsets,dtset%getvel,'getvel',idtset,iget_vel,miximage,mxnimage,ndtset_alloc)
     if(iget_vel/=0)then
!      Gather contributions to results_out from images, if needed
       if (test_img.and.mpi_enregs(idtset)%me_cell==0.and.(.not.results_gathered)) then
         call gather_results_out(dtsets,mpi_enregs,results_out,results_out_all,use_results_all, &
&         allgather=.true.,only_one_per_img=.true.)
         results_gathered=.true.
       end if
       if ((.not.test_img).or.mpi_enregs(idtset)%me_cell==0) then
         do iimage=1,nimage
           ii=mpi_enregs(idtset)%my_imgtab(iimage)
           vel_img(:,:,iimage)=zero
           vel_cell_img(:,:,iimage)=zero
           do iimage_get=1,dtsets(iget_vel)%nimage
             do jj=1,dtset%natom
               do kk=1,3
                 vel_img(kk,jj,iimage)=vel_img(kk,jj,iimage)+&
&                 miximage(ii,iimage_get)*results_out_all(iget_vel)%vel(kk,jj,iimage_get)
                 vel_cell_img(kk,jj,iimage)=vel_cell_img(kk,jj,iimage)+&
&                 miximage(ii,iimage_get)*results_out_all(iget_vel)%vel_cell(kk,jj,iimage_get)
               end do
             end do
           end do
         end do
       end if
     end if

!    In the case of parallelization over images, has to distribute data
     if (test_img) then
       if (iget_occ/=0) then
         call xmpi_bcast(occ_img,0,mpi_enregs(idtset)%comm_cell,ierr)
       end if
       if (iget_cell/=0) then
         call xmpi_bcast(acell_img,0,mpi_enregs(idtset)%comm_cell,ierr)
         call xmpi_bcast(rprim_img,0,mpi_enregs(idtset)%comm_cell,ierr)
       end if
       if (iget_vel/=0) then
         call xmpi_bcast(vel_img,0,mpi_enregs(idtset)%comm_cell,ierr)
         call xmpi_bcast(vel_cell_img,0,mpi_enregs(idtset)%comm_cell,ierr)
       end if
       if (iget_xred/=0.or.iget_xcart/=0) then
         call xmpi_bcast(xred_img,0,mpi_enregs(idtset)%comm_cell,ierr)
       end if
     end if

!    Clean memory
     ABI_FREE(miximage)
     if (test_img.and.mpi_enregs(idtset)%me_cell==0) then
       if (results_gathered) then
         call destroy_results_out(results_out_all)
       end if
       ABI_FREE(results_out_all)
     else
       nullify(results_out_all)
     end if

   end if

!  ****************************************************************************
!  Treat the pseudopotentials: initialize the psps/PAW variable

   call psps_init_from_dtset(dtset, idtset, psps, pspheads)

!  The correct dimension of pawrad/tab is ntypat. In case of alchemical psps
!  pawrad/tab(ipsp) is invoked with ipsp<=npsp. So, in order to avoid any problem,
!  declare pawrad/tab at paw_size=max(ntypat,npsp).
   paw_size=0;if (psps%usepaw==1) paw_size=max(dtset%ntypat,dtset%npsp)
   if (paw_size/=paw_size_old) then
     if (paw_size_old/=-1) then
       call pawrad_free(pawrad)
       call pawtab_free(pawtab)
       ABI_FREE(pawrad)
       ABI_FREE(pawtab)
     end if
     ABI_MALLOC(pawrad,(paw_size))
     ABI_MALLOC(pawtab,(paw_size))
     call pawtab_nullify(pawtab)
     paw_size_old=paw_size
   end if

!  ****************************************************************************
!  WVL allocations.

!  Nullify wvl_data. It is important to do so irregardless of the value of usewvl:
   call nullify_wvl_data(wvl)

!  Set up mpi information from the dataset
   if (dtset%usewvl == 1) then
#if defined HAVE_BIGDFT
     call f_malloc_set_status(iproc=mpi_enregs(idtset)%me_wvl)
     call mpi_environment_set(bigdft_mpi,mpi_enregs(idtset)%me_wvl,&
&     mpi_enregs(idtset)%nproc_wvl,mpi_enregs(idtset)%comm_wvl,&
&     mpi_enregs(idtset)%nproc_wvl)
#endif
   end if

!  ****************************************************************************
!  At this stage, all the data needed for the treatment of one dataset
!  have been transferred from multi-dataset arrays.

   iexit=0

!  Smaller integer arrays
   etotal=zero
   ABI_ALLOCATE(npwtot,(dtset%nkpt))
   npwtot = 0
   ABI_ALLOCATE(xred,(3,dtset%natom))
   ABI_ALLOCATE(occ,(dtset%mband*dtset%nkpt*dtset%nsppol))
   if(dtset%optdriver/=RUNL_GSTATE)then
     occ(:)=occ_img(:,1)
     acell(:)=acell_img(:,1)
     rprim(:,:)=rprim_img(:,:,1)
     xred(:,:)=xred_img(:,:,1)
   end if

!  ****************************************************************************
!  Exchange-correlation

   call echo_xc_name(dtset%ixc)

   if (dtset%ixc<0) then
     call libxc_functionals_init(dtset%ixc,dtset%nspden,xc_tb09_c=dtset%xc_tb09_c)

#if defined DEV_YP_VDWXC
     if ( (dtset%vdw_xc > 0) .and. (dtset%vdw_xc < 3) ) then
       vdw_params%functional = dtset%vdw_xc
       vdw_params%acutmin = dtset%vdw_df_acutmin
       vdw_params%aratio = dtset%vdw_df_aratio
       vdw_params%damax = dtset%vdw_df_damax
       vdw_params%damin = dtset%vdw_df_damin
       vdw_params%dcut = dtset%vdw_df_dcut
       vdw_params%dratio = dtset%vdw_df_dratio
       vdw_params%dsoft = dtset%vdw_df_dsoft
       vdw_params%gcut = dtset%vdw_df_gcut
       vdw_params%ndpts = dtset%vdw_df_ndpts
       vdw_params%ngpts = dtset%vdw_df_ngpts
       vdw_params%nqpts = dtset%vdw_df_nqpts
       vdw_params%nrpts = dtset%vdw_df_nrpts
       vdw_params%nsmooth = dtset%vdw_df_nsmooth
       vdw_params%phisoft = dtset%vdw_df_phisoft
       vdw_params%qcut = dtset%vdw_df_qcut
       vdw_params%qratio = dtset%vdw_df_qratio
       vdw_params%rcut = dtset%vdw_df_rcut
       vdw_params%rsoft = dtset%vdw_df_rsoft
       vdw_params%tolerance = dtset%vdw_df_tolerance
       vdw_params%tweaks = dtset%vdw_df_tweaks
       vdw_params%zab = dtset%vdw_df_zab
       write(msg,'(a,1x,a)') ch10,'[vdW-DF] *** Before init ***'
       call wrtout(std_out,msg)
       call xc_vdw_show(std_out,vdw_params)
       if ( dtset%irdvdw == 1 ) then
         write(vdw_filnam,'(a,a)') trim(filnam(3)),'_VDW.nc'
         call xc_vdw_read(vdw_filnam)
       else
         call xc_vdw_init(vdw_params)
       end if
       call xc_vdw_libxc_init(vdw_params%functional)
       write(msg,'(a,1x,a)') ch10,'[vdW-DF] *** After init ***'
       call wrtout(std_out,msg)
!      call xc_vdw_get_params(vdw_params)
       call xc_vdw_memcheck(std_out)
       call xc_vdw_show(std_out)
       call xc_vdw_show(ab_out)

       write (msg,'(a,1x,a,e10.3,a)')ch10,&
&       '[vdW-DF] activation threshold: vdw_df_threshold=',dtset%vdw_df_threshold,ch10
       call xc_vdw_trigger(.false.)
       call wrtout(std_out,msg)
     end if
#else
     if ( (dtset%vdw_xc > 0) .and. (dtset%vdw_xc < 3) ) then
       write(msg,'(3a)')&
&       'vdW-DF functionals are not fully operational yet.',ch10,&
&       'Action: modify vdw_xc'
       MSG_ERROR(msg)
     end if
#endif
   end if

   ! FFTW3 threads initialization
   if (dtset%ngfft(7) / 100 == FFT_FFTW3) call fftw3_init_threads()

   ! Set precision for FFT libs.
   ii = fftcore_set_mixprec(dtset%mixprec)

!  linalg initialisation
   linalg_max_size=maxval(dtset%nband(:))
   call abi_linalg_init(linalg_max_size,dtset%optdriver,dtset%wfoptalg,dtset%paral_kgb,&
        dtset%use_gpu_cuda,dtset%use_slk,dtset%np_slk,mpi_enregs(idtset)%comm_bandspinorfft)

   call timab(642,2,tsec)

!  ****************************************************************************
!  Main case selection in driver

   select case(dtset%optdriver)

   case(RUNL_GSTATE)

     ABI_ALLOCATE(fcart_img,(3,dtset%natom,nimage))
     ABI_ALLOCATE(fred_img,(3,dtset%natom,nimage))
     ABI_ALLOCATE(intgres_img,(dtset%nspden,dtset%natom,nimage))
     ABI_ALLOCATE(etotal_img,(nimage))
     ABI_ALLOCATE(strten_img,(6,nimage))

     call gstateimg(acell_img,amu_img,codvsn,cpui,dtfil,dtset,etotal_img,fcart_img,&
       fred_img,iexit,intgres_img,mixalch_img,mpi_enregs(idtset),nimage,npwtot,occ_img,&
       pawang,pawrad,pawtab,psps,rprim_img,strten_img,vel_cell_img,vel_img,wvl,xred_img,&
       filnam,filstat,idtset,jdtset_,ndtset)

   case(RUNL_RESPFN)

     call respfn(codvsn,cpui,dtfil,dtset,etotal,iexit,mkmems,mpi_enregs(idtset),&
       npwtot,occ,pawang,pawrad,pawtab,psps,results_respfn,xred)

   case(RUNL_SCREENING)
     call screening(acell,codvsn,dtfil,dtset,pawang,pawrad,pawtab,psps,rprim)

   case(RUNL_SIGMA)
     call sigma(acell,codvsn,dtfil,dtset,pawang,pawrad,pawtab,psps,rprim,converged)

   case(RUNL_NONLINEAR)
     call nonlinear(codvsn,dtfil,dtset,etotal,mpi_enregs(idtset),npwtot,occ,pawang,pawrad,pawtab,psps,xred)

   case (RUNL_BSE)
     call bethe_salpeter(acell,codvsn,dtfil,dtset,pawang,pawrad,pawtab,psps,rprim,xred)

   case(RUNL_GWLS)
     ! For running G0W0 calculations with Lanczos basis for dielectric operator
     ! and Sternheimer equation for avoiding the use of conduction states (MC+JJL)
     ABI_ALLOCATE(etotal_img,(nimage))
     ABI_ALLOCATE(fcart_img,(3,dtset%natom,nimage))
     ABI_ALLOCATE(fred_img,(3,dtset%natom,nimage))
     ABI_ALLOCATE(intgres_img,(dtset%nspden,dtset%natom,nimage))
     ABI_ALLOCATE(strten_img,(6,nimage))

     call gwls_sternheimer(acell_img,amu_img,codvsn,cpui,dtfil,dtset,etotal_img,fcart_img,&
       fred_img,iexit,intgres_img,mixalch_img,mpi_enregs(idtset),nimage,npwtot,occ_img,&
       pawang,pawrad,pawtab,psps,rprim_img,strten_img,vel_cell_img,vel_img,xred_img,&
       filnam,filstat,idtset,jdtset_,ndtset)

     ABI_FREE(etotal_img)
     ABI_FREE(fcart_img)
     ABI_FREE(fred_img)
     ABI_FREE(intgres_img)
     ABI_FREE(strten_img)

   case (RUNL_WFK)
     call wfk_analyze(acell,codvsn,dtfil,dtset,pawang,pawrad,pawtab,psps,rprim,xred)

   case (RUNL_EPH)
     call dtsets(0)%free_nkpt_arrays()
     call dtsets(idtset)%free_nkpt_arrays()
     call eph(acell,codvsn,dtfil,dtset,pawang,pawrad,pawtab,psps,rprim,xred)

   case(RUNL_LONGWAVE)

     call longwave(codvsn,dtfil,dtset,etotal,mpi_enregs(idtset),npwtot,occ,&
&     pawrad,pawtab,psps,xred)

   case default
     ! Bad value for optdriver
     write(msg,'(a,i0,4a)')&
      'Unknown value for the variable optdriver: ',dtset%optdriver,ch10,&
      'This is not allowed.',ch10, 'Action: modify optdriver in the input file.'
     MSG_ERROR(msg)
   end select

   call timab(643,1,tsec)
!  ****************************************************************************

!  Transfer of multi dataset outputs from temporaries:
!  acell, xred, occ rprim, and vel might be modified from their input values
!  etotal, fcart, fred, intgres, and strten have been computed
!  npwtot was already computed before, but is stored only now

   if(dtset%optdriver==RUNL_GSTATE)then
     do iimage=1,nimage
       results_out(idtset)%etotal(iimage)                 =etotal_img(iimage)
       results_out(idtset)%acell(:,iimage)                =acell_img(:,iimage)
       results_out(idtset)%amu(1:dtset%ntypat,iimage)     =amu_img(:,iimage)
       results_out(idtset)%rprim(:,:,iimage)              =rprim_img(:,:,iimage)
       results_out(idtset)%strten(:,iimage)                =strten_img(:,iimage)
       results_out(idtset)%fcart(1:3,1:dtset%natom,iimage)=fcart_img(:,:,iimage)
       results_out(idtset)%fred(1:3,1:dtset%natom,iimage) =fred_img(:,:,iimage)
       if(dtset%nspden/=2)then
         results_out(idtset)%intgres(1:dtset%nspden,1:dtset%natom,iimage) =intgres_img(:,:,iimage)
       else
         results_out(idtset)%intgres(1,1:dtset%natom,iimage) =intgres_img(1,:,iimage)
         results_out(idtset)%intgres(4,1:dtset%natom,iimage) =intgres_img(2,:,iimage)
       endif
       results_out(idtset)%mixalch(1:dtset%npspalch,1:dtset%ntypalch,iimage) &
&       =mixalch_img(1:dtset%npspalch,1:dtset%ntypalch,iimage)
       results_out(idtset)%npwtot(1:dtset%nkpt,iimage)    =npwtot(1:dtset%nkpt)
       results_out(idtset)%occ(1:dtset%mband*dtset%nkpt*dtset%nsppol,iimage)=&
&       occ_img(1:dtset%mband*dtset%nkpt*dtset%nsppol,iimage)
       results_out(idtset)%vel(:,1:dtset%natom,iimage)    =vel_img(:,:,iimage)
       results_out(idtset)%vel_cell(:,:,iimage)           =vel_cell_img(:,:,iimage)
       results_out(idtset)%xred(:,1:dtset%natom,iimage)   =xred_img(:,:,iimage)
     end do
     ABI_FREE(etotal_img)
     ABI_FREE(fcart_img)
     ABI_FREE(fred_img)
     ABI_FREE(intgres_img)
     ABI_FREE(strten_img)
   else
     results_out(idtset)%acell(:,1)                =acell(:)
     results_out(idtset)%etotal(1)                 =etotal
     results_out(idtset)%rprim(:,:,1)              =rprim(:,:)
     results_out(idtset)%npwtot(1:dtset%nkpt,1)    =npwtot(1:dtset%nkpt)
     results_out(idtset)%occ(1:dtset%mband*dtset%nkpt*dtset%nsppol,1)=&
&     occ(1:dtset%mband*dtset%nkpt*dtset%nsppol)
     results_out(idtset)%xred(:,1:dtset%natom,1)   =xred(:,:)
   end if
   ABI_FREE(acell_img)
   ABI_FREE(amu_img)
   ABI_FREE(mixalch_img)
   ABI_FREE(occ_img)
   ABI_FREE(rprim_img)
   ABI_FREE(vel_img)
   ABI_FREE(vel_cell_img)
   ABI_FREE(xred_img)

   if (dtset%ngfft(7) / 100 == FFT_FFTW3) call fftw3_cleanup()

   if (dtset%ixc<0) then
     call libxc_functionals_end()

#if defined DEV_YP_VDWXC
     if ( (dtset%vdw_xc > 0) .and. (dtset%vdw_xc < 10) ) then
       if ( dtset%prtvdw /= 0 ) then
         write(vdw_filnam,'(a,a)') trim(filnam(4)),'_VDW.nc'
         call xc_vdw_write(vdw_filnam)
       end if
       call xc_vdw_show(std_out,vdw_params)
       call xc_vdw_done(vdw_params)
     end if
#endif
   end if

   ! MG: There are routines such as GW and Berry phase that can change/compute
   ! entries in the datatype at run-time. These changes won't be visibile
   ! in the main output file since we are passing a copy of dtsets.
   ! I tried to update the results with the call below but this creates
   ! several problems in outvars since one should take into account
   ! the new dimensions (e.g. nkptgw) and their maximum value.
   ! For the time being, we continue to pass a copy of dtsets(idtset).
   !call dtsets(idtset)%free()
   !call dtsets(idtset) = dtset%copy()

   call dtset%free()

   ABI_FREE(occ)
   ABI_FREE(xred)
   ABI_FREE(npwtot)

   call abi_linalg_finalize()
   call xg_finalize()

   ! Check whether exiting was required by the user.
   ! If found then beat a hasty exit from time steps
   openexit=1; if(dtset%chkexit==0) openexit=0
   call exit_check(zero,dtfil%filnam_ds(1),iexit,ab_out,mpi_enregs(idtset)%comm_cell,openexit)

   call timab(643,2,tsec)

   if (iexit/=0) exit
 end do ! idtset (allocate statements are present - an exit statement is present)

!*********************************************************************

 call timab(644,1,tsec)

 !PSP deallocation
 call psps_free(psps)

 !XG 121126 : One should not use dtset or idtset in this section, as these might not be defined for all processors.

 !PAW deallocation
 if (allocated(pawrad)) then
   call pawrad_free(pawrad)
   ABI_FREE(pawrad)
 end if
 if (allocated(pawtab)) then
   call pawtab_free(pawtab)
   ABI_FREE(pawtab)
 end if
 call pawang_free(pawang)

 ABI_FREE(jdtset_)

!Results_respfn deallocation
 call destroy_results_respfn(results_respfn)

#if defined HAVE_BIGDFT
!XG 121126 : NOTE that the next debugging section was quite problematic : indeed we are
!outside the loop over datasets, so the value of dtset%usewvl, that is referring to the last
!dtset, might not be initialized, if the dataset is NOT treated by the processor.
!See line 253 for the initialisation of dtset, while this might not happen if
!if(mpi_enregs(idtset)%me<0) cycle  , as defined by line 216
!if (dtset%usewvl == 1) then
!!  WVL - debugging
!call memocc_abi(0,mpi_enregs(1)%me_wvl,'count','stop')
!call wvl_timing(mpi_enregs(1)%me_wvl,'WFN_OPT','PR')
!call wvl_timing(mpi_enregs(1)%me_wvl,'              ','RE')
!end if
#endif

 call timab(644,2,tsec)
 call timab(640,2,tsec)

 DBG_EXIT("COLL")

end subroutine driver
!!***

end module m_driver
!!***
