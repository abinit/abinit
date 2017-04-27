!{\src2tex{textfont=tt}}
!!****f* ABINIT/gstateimg
!! NAME
!! gstateimg
!!
!! FUNCTION
!! Routine for conducting DFT calculations for a set of (dynamical) images
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (XG, AR, GG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  codvsn=code version
!!  cpui=initial CPU time
!!  nimage=number of images of the cell (treated by current proc)
!!  === Optional arguments (needed when nimage>1) ===
!!    filnam(5)=character strings giving file names
!!    filstat=character strings giving name of status file
!!    idtset=index of the dataset
!!    jdtset(0:ndtset)=actual index of the datasets
!!    ndtset=number of datasets
!!
!! OUTPUT
!!  etotal_img=total energy, for each image
!!  fcart_img(3,natom,nimage)=forces, in cartesian coordinates, for each image
!!  fred_img(3,natom,nimage)=forces, in reduced coordinates, for each image
!!  npwtot(nkpt) = total number of plane waves at each k point
!!  strten_img(6,nimage)=stress tensor, for each image
!!
!! SIDE EFFECTS
!!  acell_img(3,nimage)=unit cell length scales (bohr), for each image
!!  amu_img(ntypat,nimage)=value of mass for each atomic type, for each image
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | mband =maximum number of bands (IN)
!!   | mgfft =maximum single fft dimension (IN)
!!   | mkmem =maximum number of k points which can fit in core memory (IN)
!!   | mpw   =maximum number of planewaves in basis sphere (large number) (IN)
!!   | natom =number of atoms in unit cell (IN)
!!   | nfft  =(effective) number of FFT grid points (for this processor) (IN)
!!   | nkpt  =number of k points (IN)
!!   | nspden=number of spin-density components (IN)
!!   | nsppol=number of channels for spin-polarization (1 or 2) (IN)
!!   | nsym  =number of symmetry elements in space group
!!  iexit= exit flag
!!  mixalch_img(npspalch,ntypalch,nimage)=value of alchemical mixing factors,for each image
!!  mpi_enreg=MPI-parallelisation information (some already initialized,
!!            some others to be initialized here)
!!  occ_img(mband*nkpt*nsppol,nimage) = occupation number for each band and k, for each image
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   Before entering the first time in gstateimg, a significant part of
!!   psps has been initialized :
!!   the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid,
!!     ntypat,n1xccc,usepaw,useylm, and the arrays dimensioned to npsp
!!   All the remaining components of psps are to be initialized in the call
!!   to pspini .
!!   The next time the code enters gstateimg, psps might be identical to the
!!   one of the previous dtset, in which case, no reinitialisation is scheduled
!!   in pspini.f .
!!  rprim_img(3,3,nimage)=dimensionless real space primitive translations, for each image
!!  vel_cell_img(3,3,nimage)=value of cell parameters velocities, for each image
!!  vel_img(3,natom,nimage)=value of atomic velocities,for each image
!!  xred_img(3,natom,nimage) = reduced atomic coordinates, for each image
!!
!! NOTES
!! USE OF FFT GRIDS:
!! =================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut)
!!      for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ...
!!      are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg)
!!      for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ...
!!      Total density, potentials, ...
!!      are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!! In case of wavelets:
!! --------------------
!!    - Only the usual FFT grid (defined by wvl_crmult) is used.
!!      It is defined by nfft, ngfft, mgfft, ... This is strictly not
!!      an FFT grid since its dimensions are not suited for FFTs. They are
!!      defined by wvl_setngfft().
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! TODO
!! Not yet possible to use restartxf in parallel when localrdwf==0
!!
!! PARENTS
!!      driver,gwls_sternheimer
!!
!! CHILDREN
!!      args_gs_free,args_gs_init,copy_results_img,destroy_results_img
!!      dtfil_init,ga_destroy,ga_init,gstate,init_results_img
!!      libpaw_spmsg_mpisum,localfilnam,localrdfile,localredirect,localwrfile
!!      mep_destroy,mep_init,pimd_destroy,pimd_init,predictimg,prtimg
!!      scf_history_free,scf_history_nullify,specialmsg_mpisum,status,timab
!!      wrtout,xmpi_barrier,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine gstateimg(acell_img,amu_img,codvsn,cpui,dtfil,dtset,etotal_img,fcart_img,&
&                    fred_img,iexit,mixalch_img,mpi_enreg,nimage,npwtot,occ_img,&
&                    pawang,pawrad,pawtab,psps,&
&                    rprim_img,strten_img,vel_cell_img,vel_img,wvl,xred_img,&
&                    filnam,filstat,idtset,jdtset,ndtset) ! optional arguments


 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use defs_parameters
 use defs_rectypes
 use m_profiling_abi
 use m_abihist
 use m_mep
 use m_ga
 use m_pimd
 use m_xmpi
 use m_errors
 use m_rec
 use m_args_gs
 use m_results_img
 use m_scf_history
 use m_io_redirect
 use m_specialmsg, only : specialmsg_mpisum

 use m_libpaw_tools, only : libpaw_spmsg_mpisum
 use m_pawang,       only : pawang_type
 use m_pawrad,       only : pawrad_type
 use m_pawtab,       only : pawtab_type

#if defined  HAVE_BIGDFT
 use BigDFT_API, only: mpi_environment_set
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gstateimg'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_45_geomoptim
 use interfaces_67_common
 use interfaces_95_drive, except_this_one => gstateimg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nimage
 integer,optional,intent(in) :: idtset,ndtset
 integer,intent(inout) :: iexit
 real(dp),intent(in) :: cpui
 character(len=6),intent(in) :: codvsn
 character(len=fnlen),optional,intent(in) :: filstat
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),target,intent(inout) :: dtfil
 type(dataset_type),target,intent(inout) :: dtset
 type(pawang_type),intent(inout) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
 type(wvl_data),intent(inout) :: wvl
!arrays
 integer,optional,intent(in) :: jdtset(:)
 integer,intent(out) :: npwtot(dtset%nkpt)
 character(len=fnlen),optional,intent(in) :: filnam(:)
 real(dp), intent(out) :: etotal_img(nimage),fcart_img(3,dtset%natom,nimage)
 real(dp), intent(out) :: fred_img(3,dtset%natom,nimage),strten_img(6,nimage)
 real(dp),intent(inout) :: acell_img(3,nimage),amu_img(dtset%ntypat,nimage)
 real(dp),intent(inout) :: mixalch_img(dtset%npspalch,dtset%ntypalch,nimage)
 real(dp),intent(inout) :: occ_img(dtset%mband*dtset%nkpt*dtset%nsppol,nimage)
 real(dp),intent(inout) :: rprim_img(3,3,nimage),vel_cell_img(3,3,nimage),vel_img(3,dtset%natom,nimage)
 real(dp),intent(inout) :: xred_img(3,dtset%natom,nimage)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!Define file format for different type of files. Presently,
!only one file format is supported for each type of files, but this might
!change soon ...
!2   for wavefunction file, new format (version 2.0 and after)    (fform)   NOT USED
!52  for density rho(r)       (fformr)
!102 for potential V(r) file. (fformv)  NOT USED
!scalars
 integer,parameter :: formeig=0,level=100,ndtpawuj=0,response=0
 integer :: history_size,idelta,idynimage,ierr,ifirst
 integer :: ii,iimage,ih,itimimage,itimimage_eff,itimimage_prev,ndynimage,nocc
 integer :: ntimimage,ntimimage_stored,ntimimage_max,similar
 logical :: check_conv,compute_all_images,compute_static_images
 logical :: isVused,isARused,is_master,is_mep,is_pimd
 logical :: call_predictor,use_hist,use_hist_prev
 real(dp) :: delta_energy,dtion
 character(len=500) :: hist_filename,msg
 type(args_gs_type) :: args_gs
 type(mep_type) :: mep_param
 type(ga_type) :: ga_param
 type(pimd_type) :: pimd_param
!arrays
 integer,allocatable :: list_dynimage(:),scf_initialized(:)
 character(len=60),parameter :: imagealgo_str(0:13)=(/ &
&   'IMAGE COPY                                                  ',& ! 0
&   'IMAGE STEEPEST DESCENT                                      ',& ! 1
&   'STRING METHOD                                               ',& ! 2
&   'METADYNAMICS                                                ',& ! 3
&   'GENETIC ALGORITHM                                           ',& ! 4
&   'NUDGED ELASTIC BAND                                         ',& ! 5
&   '                                                            ',& ! 6
&   '                                                            ',& ! 7
&   '                                                            ',& ! 8
&   'PATH-INTEGRAL MOLECULAR DYNAMICS (LANGEVIN)                 ',& ! 9
&   'PATH-INTEGRAL MOLECULAR DYNAMICS (QUANTUM THERMAL BATH)     ',& ! 10
&   '                                                            ',& ! 11
&   '                                                            ',& ! 12
&   'PATH-INTEGRAL MOLECULAR DYNAMICS (CHAIN OF THERMOSTATS)     '/) ! 13
 real(dp) :: acell(3),rprim(3,3),rprimd(3,3),tsec(2)
 real(dp),allocatable :: amass(:,:),occ(:),vel(:,:),vel_cell(:,:),xred(:,:)
 type(abihist),allocatable :: hist(:),hist_prev(:)
 type(results_img_type),pointer :: results_img_timimage(:,:),res_img(:)
 type(scf_history_type),allocatable :: scf_history(:)

! ***********************************************************************

 DBG_ENTER("COLL")

 call timab(700,1,tsec)
 call timab(703,3,tsec)

 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Arguments check
 if (dtset%nimage>1) then
   if ((.not.present(filnam)).or.(.not.present(filnam)).or.(.not.present(idtset)).or.&
&   (.not.present(ndtset)).or.(.not.present(jdtset))) then
     write(msg,'(3a)') &
&     'When nimage>1, all the following argument should be present:',ch10,&
&     'filnam, filstat, idtset, ndtset, jdtset  !'
     MSG_BUG(msg)
   end if
 end if

!Set flag for the effective computation (once) of static images
!For the time being only set when parallelization is activated
!Note: if you modify this flag, do not forget to change it in outvars and outvar1
 compute_static_images=(dtset%istatimg>0)

!Prepare the calculation, by computing flags and dimensions
 is_pimd=(dtset%imgmov==9.or.dtset%imgmov==10.or.dtset%imgmov==13)
 is_mep =(dtset%imgmov==1.or.dtset%imgmov== 2.or.dtset%imgmov== 5)
 ntimimage=dtset%ntimimage
 ntimimage_stored=ntimimage;if(is_pimd)ntimimage_stored=2
 nocc=dtset%mband*dtset%nkpt*dtset%nsppol
 is_master=(mpi_enreg%me_cell==0.and.mpi_enreg%me_img==0)
 delta_energy=zero

!Management of dynamics/relaxation history (positions, forces, stresses, ...)
 use_hist=(dtset%imgmov/=0.and.nimage>0) ; use_hist_prev=.false.
 isVused=is_pimd;isARused=(dtset%optcell/=0)
 if (use_hist) then
   !Read history from file (and broadcast if MPI)
#if defined HAVE_NETCDF
   use_hist_prev=(dtset%restartxf==-1.and.nimage>0)
#endif
   hist_filename=trim(dtfil%filnam_ds(4))//'_HIST.nc'
   if (use_hist_prev)then
     ABI_DATATYPE_ALLOCATE(hist_prev,(nimage))
     if (mpi_enreg%me_cell==0) then
       call read_md_hist_img(hist_filename,hist_prev,isVused,isARused,&
&                            imgtab=mpi_enreg%my_imgtab)
     end if
     call abihist_bcast(hist_prev,0,mpi_enreg%comm_cell)
     if (nimage>0) then
       if (any(hist_prev(:)%mxhist/=hist_prev(1)%mxhist)) then
         msg='History problem: all images should have the same number of time steps!'
         MSG_ERROR(msg)
       end if
       use_hist_prev=(hist_prev(1)%mxhist>0)
       if (use_hist_prev) ntimimage=ntimimage+hist_prev(1)%mxhist
     end if
     if (.not.use_hist_prev) then
       call abihist_free(hist_prev)
       ABI_DATATYPE_DEALLOCATE(hist_prev)
     end if
   end if
   !Initialize a variable to write the history
   ABI_DATATYPE_ALLOCATE(hist,(nimage))
   call abihist_init(hist,dtset%natom,ntimimage,isVused,isARused)
 end if ! imgmov/=0

!Allocations
 ABI_ALLOCATE(occ,(nocc))
 ABI_ALLOCATE(vel,(3,dtset%natom))
 ABI_ALLOCATE(vel_cell,(3,3))
 ABI_ALLOCATE(xred,(3,dtset%natom))
 ABI_DATATYPE_ALLOCATE(results_img_timimage,(nimage,ntimimage_stored))
 ABI_ALLOCATE(list_dynimage,(dtset%ndynimage))
 do itimimage=1,ntimimage_stored
   res_img => results_img_timimage(:,itimimage)
   call init_results_img(dtset%natom,dtset%npspalch,dtset%nsppol,dtset%ntypalch,&
&                        dtset%ntypat,res_img)
   do iimage=1,nimage
     res_img(iimage)%acell(:)     =acell_img(:,iimage)
     res_img(iimage)%amu(:)       =amu_img(:,iimage)
     res_img(iimage)%mixalch(:,:) =mixalch_img(:,:,iimage)
     res_img(iimage)%rprim(:,:)   =rprim_img(:,:,iimage)
     res_img(iimage)%xred(:,:)    =xred_img(:,:,iimage)
     res_img(iimage)%vel(:,:)     =vel_img(:,:,iimage)
     res_img(iimage)%vel_cell(:,:)=vel_cell_img(:,:,iimage)
   end do
 end do
 ndynimage=0
 do iimage=1,nimage
   ii=mpi_enreg%my_imgtab(iimage)
   if (dtset%dynimage(ii)==1) then
     ndynimage=ndynimage+1
     list_dynimage(ndynimage)=iimage
   end if
 end do

!Management of SCF history (density/WF predictions from one time step to another)
 ABI_DATATYPE_ALLOCATE(scf_history,(nimage))
 ABI_ALLOCATE(scf_initialized,(nimage))
 scf_initialized=0
 history_size=-1
 if (dtset%ntimimage<=1) then
   if (dtset%usewvl==0.and.dtset%ionmov>0.and. &
&   (abs(dtset%densfor_pred)==5.or.abs(dtset%densfor_pred)==6)) history_size=2
 else
   if (abs(dtset%densfor_pred)==2.or.abs(dtset%densfor_pred)==3) history_size=0
   if (dtset%usewvl==0.and.(abs(dtset%densfor_pred)==5.or.abs(dtset%densfor_pred)==6)) history_size=2
 end if
 do iimage=1,nimage
   call scf_history_nullify(scf_history(iimage))
   scf_history(iimage)%history_size=history_size
 end do

!MEP search: fill in eventually the data structure mep_param
 call mep_init(dtset,mep_param)

!GA search: fill in eventually the data structure ga_param
 call ga_init(dtset,ga_param)

!PIMD: fill in eventually the data structure pimd_param
 call pimd_init(dtset,pimd_param,is_master)
 dtion=one;if (is_pimd) dtion=pimd_param%dtion

!In some cases, need amass variable
 if (use_hist.and.is_pimd) then
   ABI_ALLOCATE(amass,(dtset%natom,nimage))
   do iimage=1,nimage
     if (any(amu_img(:,iimage)/=amu_img(:,1))) then
       msg='HIST file is not compatible with variable masses!'
       MSG_ERROR(msg)
     end if
     amass(:,iimage)=amu_emass*amu_img(dtset%typat(:),iimage)
   end do
 end if

!In the case of the 4th-order Runge-Kutta solver,
!one must have a number of step multiple of 4.
 ntimimage_max=ntimimage;idelta=1
 if (dtset%imgmov==2.and.mep_param%mep_solver==4) then
   ntimimage_max=4*(ntimimage_max/4)
   idelta=4
 end if

 call timab(703,2,tsec)

!-----------------------------------------------------------------------------------------
!Big loop on the propagation of all images
 itimimage_eff=1
 do itimimage=1,ntimimage

   res_img => results_img_timimage(:,itimimage_eff)
   call_predictor=(ntimimage>1)

!  If history is activated and if current image is inside it: do not compute anything
   if (use_hist_prev) then
     if (all(hist_prev(:)%ihist<=hist_prev(:)%mxhist)) then
       do iimage=1,nimage
         ih=hist_prev(iimage)%ihist
         call abihist_copy(hist_prev(iimage),hist(iimage))
         call mkradim(hist_prev(iimage)%acell(:,ih),rprim,hist_prev(iimage)%rprimd(:,:,ih))
         res_img(iimage)%acell(:)=hist_prev(iimage)%acell(:,ih)
         res_img(iimage)%rprim(:,:)=rprim
         res_img(iimage)%xred(:,:)=hist_prev(iimage)%xred(:,:,ih)
         res_img(iimage)%vel(:,:)=hist_prev(iimage)%vel(:,:,ih)
         res_img(iimage)%vel_cell(:,:)=hist_prev(iimage)%vel_cell(:,:,ih)
         res_img(iimage)%results_gs%fcart(:,:)=hist_prev(iimage)%fcart(:,:,ih)
         res_img(iimage)%results_gs%strten(:)=hist_prev(iimage)%strten(:,ih)
         res_img(iimage)%results_gs%etotal=hist_prev(iimage)%etot(ih)
         res_img(iimage)%results_gs%energies%entropy=hist_prev(iimage)%entropy(ih)
         call fcart2fred(res_img(iimage)%results_gs%fcart,res_img(iimage)%results_gs%fred,&
&                        hist_prev(iimage)%rprimd(:,:,ih),dtset%natom)
         hist_prev(iimage)%ihist=hist_prev(iimage)%ihist+1
       end do
       !PI-QTB: skip a record in random force file
       if (pimd_param%use_qtb==1) call pimd_skip_qtb(pimd_param)
       !call_predictor=.false.
       goto 110 ! This is temporary
     end if
   end if

   call timab(704,1,tsec)
   call localfilnam(mpi_enreg%comm_img,mpi_enreg%comm_cell,mpi_enreg%comm_world,filnam,'_IMG',dtset%nimage)
   compute_all_images=(compute_static_images.and.itimimage==1)

!  Print title for time step
   if (dtset%nimage>1.or.dtset%ntimimage>1) then
     if (dtset%prtvolimg<2) then
       msg=ch10;if (itimimage >1) write(msg,'(2a)') ch10,ch10
       write(msg,'(5a)') trim(msg),&
&       '================================================================================',&
&       ch10,' ',trim(imagealgo_str(dtset%imgmov))
     else
       msg='';if (itimimage >1) msg=ch10
       write(msg,'(5a)') trim(msg),&
&       '--------------------------------------------------------------------------------',&
&       ch10,' ',trim(imagealgo_str(dtset%imgmov))
     end if
     if (dtset%ntimimage==1) write(msg,'(2a)')    trim(msg),' FOR 1 TIME STEP'
     if (dtset%ntimimage >1) write(msg,'(2a,i5)') trim(msg),' - TIME STEP ',itimimage
     if (dtset%prtvolimg<2) then
       write(msg,'(3a)') trim(msg),ch10,&
&       '================================================================================'
     end if
     call wrtout(ab_out ,msg,'COLL')
     call wrtout(std_out,msg,'PERS')
   end if

   call timab(704,2,tsec)

!  Loop on the dynamical images
   idynimage=0
   do iimage=1,nimage

     call timab(705,1,tsec)

     ii=mpi_enreg%my_imgtab(iimage)
     if (dtset%dynimage(ii)==1) idynimage=idynimage+1

!    Compute static image only at first time step
     if (dtset%dynimage(ii)==1.or.compute_all_images) then

!      Change file names according to image index (if nimage>1)
       if (dtset%nimage>1) then
         call dtfil_init(dtfil,dtset,filnam,filstat,idtset,jdtset,mpi_enreg,ndtset,&
&         image_index=ii)
         if (itimimage>1) then
           dtfil%ireadwf=0;dtfil%ireadden=0;dtfil%ireadkden=0
         end if
       end if

!      Redefine output units
       call localwrfile(mpi_enreg%comm_cell,ii,dtset%nimage,mpi_enreg%paral_img,dtset%prtvolimg)

!      Print title for image
       if (dtset%nimage>1.and.(dtset%prtvolimg==0.or.do_write_log)) then
         if (ii==1) write(msg,'(a)' ) ch10
         if (ii >1) write(msg,'(2a)') ch10,ch10
         write(msg,'(6a,i4,a,i4,3a)') trim(msg),&
&         '--------------------------------------------------------------------------------',ch10,&
&         ' ',trim(imagealgo_str(dtset%imgmov)),' - CELL # ',ii,'/',dtset%nimage,ch10,&
&         '--------------------------------------------------------------------------------',ch10
         if (dtset%prtvolimg==0) then
           call wrtout(ab_out ,msg,'COLL')
         end if
         if (do_write_log) then
           call wrtout(std_out,msg,'PERS')
         end if
       end if

       acell(:)     =res_img(iimage)%acell(:)
       rprim(:,:)   =res_img(iimage)%rprim(:,:)
       vel(:,:)     =res_img(iimage)%vel(:,:)
       vel_cell(:,:)=res_img(iimage)%vel_cell(:,:)
       xred(:,:)    =res_img(iimage)%xred(:,:)
       occ(:)       =occ_img(:,iimage)

       call args_gs_init(args_gs, &
&       res_img(iimage)%amu(:),&
&       res_img(iimage)%mixalch(:,:),&
&       dtset%dmatpawu(:,:,:,:,ii),dtset%upawu(:,ii),dtset%jpawu(:,ii))

       call timab(705,2,tsec)

       call status(idynimage+100*itimimage,dtfil%filstat,iexit,level,'call gstate   ')
       call gstate(args_gs,acell,codvsn,cpui,dtfil,dtset,iexit,scf_initialized(iimage),&
&       mpi_enreg,npwtot,occ,pawang,pawrad,pawtab,psps,&
&       res_img(iimage)%results_gs,&
&       rprim,scf_history(iimage),vel,vel_cell,wvl,xred)

       call timab(706,1,tsec)

       call args_gs_free(args_gs)

       if (dtset%dynimage(ii)==1) then
         res_img(iimage)%acell(:)     =acell(:)
         res_img(iimage)%rprim(:,:)   =rprim(:,:)
         res_img(iimage)%vel(:,:)     =vel(:,:)
         res_img(iimage)%vel_cell(:,:)=vel_cell(:,:)
         res_img(iimage)%xred(:,:)    =xred(:,:)
         occ_img(:,iimage)            =occ(:)
       end if

!    check change de rprim et reecriture dans hist
!    check change de xred et reecriture dans hist

!      Close output units ; restore defaults
       call localredirect(mpi_enreg%comm_cell,mpi_enreg%comm_world,dtset%nimage,mpi_enreg%paral_img,dtset%prtvolimg)
       call timab(706,2,tsec)

     else if (itimimage>1) then ! For static images, simply copy one time step to the other
       itimimage_prev=itimimage_eff-1
       if (itimimage_prev<1) itimimage_prev=ntimimage_stored
       call copy_results_img(results_img_timimage(iimage,itimimage_prev), &
&                            results_img_timimage(iimage,itimimage_eff ))
     end if

!    Store results in hist datastructure
     if (use_hist) then
       ih=hist(iimage)%ihist
       call mkrdim(res_img(iimage)%acell(:),res_img(iimage)%rprim(:,:),rprimd)
       call var2hist(res_img(iimage)%acell(:),hist(iimage),dtset%natom,&
&                    rprimd,res_img(iimage)%xred(:,:),.FALSE.)
       call vel2hist(amass(:,iimage),hist(iimage),res_img(iimage)%vel(:,:),&
&                    res_img(iimage)%vel_cell(:,:))
       hist(iimage)%fcart(:,:,ih)=res_img(iimage)%results_gs%fcart(:,:)
       hist(iimage)%strten(:,ih)=res_img(iimage)%results_gs%strten(:)
       hist(iimage)%etot(ih)=res_img(iimage)%results_gs%etotal
       hist(iimage)%entropy(ih)=res_img(iimage)%results_gs%energies%entropy
       hist(iimage)%time(ih)=real(itimimage,kind=dp)*dtion
     end if

   end do ! iimage

   if(mpi_enreg%paral_img==1)then
     call timab(702,1,tsec)
     call xmpi_barrier(mpi_enreg%comm_img)
     call timab(702,2,tsec)
   end if

   call timab(707,1,tsec)

!  Output when images are used
   if (dtset%nimage>1) then
!    === 1st option: reduced outputs ===
     if (dtset%prtvolimg>0) then
       call prtimg(dtset%dynimage,imagealgo_str(dtset%imgmov),dtset%imgmov,ab_out,&
&       mpi_enreg,nimage,dtset%nimage,compute_all_images,dtset%prtvolimg,res_img)
     end if
   end if

!  Manage log files when images are used
   call localrdfile(mpi_enreg%comm_img,mpi_enreg%comm_world,compute_all_images,&
&   dtset%nimage,mpi_enreg%paral_img,dtset%prtvolimg,dyn=dtset%dynimage)

!  Write hist datastructure in HIST file
#if defined HAVE_NETCDF
   if (use_hist.and.mpi_enreg%me_cell==0) then
     ifirst=merge(0,1,itimimage>1)
     call write_md_hist_img(hist,hist_filename,ifirst,dtset%natom,dtset%ntypat,&
&                           dtset%typat,amu_img(:,1),dtset%znucl,dtion,&
&                           nimage=dtset%nimage,comm_img=mpi_enreg%comm_img,&
&                           imgtab=mpi_enreg%my_imgtab)
   end if
#endif

!  TESTS WHETHER ONE CONTINUES THE LOOP
!  Here we calculate the change in energy, and exit if delta_energy < tolimg
   delta_energy=zero
!  Doesn't check convergence in case of PIMD
   check_conv=((.not.is_pimd).and.itimimage>1)
!  In case of 4th-order Runge-Kutta, does check convergence every 4 steps
   if (dtset%imgmov==2.and.mep_param%mep_solver==4) then
     check_conv=(mod(itimimage,4)==0.and.itimimage>4)
   end if
   if (check_conv) then
     do idynimage=1,ndynimage
       iimage=list_dynimage(idynimage)
       delta_energy=delta_energy &
&       +abs(results_img_timimage(iimage,itimimage)%results_gs%etotal &
&           -results_img_timimage(iimage,itimimage-idelta)%results_gs%etotal)
     end do
     if (mpi_enreg%paral_img==1) then
       call xmpi_sum(delta_energy,mpi_enreg%comm_img,ierr)
     end if
     delta_energy=delta_energy/dtset%ndynimage
     if (delta_energy<dtset%tolimg) then
       if (dtset%prtvolimg<2) then
         write(msg,'(5a,i5,6a,es11.3,a,es11.3,2a)') ch10,ch10,&
&         '================================================================================',ch10,&
&         ' At time step ',itimimage,ch10,&
&         ' ',trim(imagealgo_str(dtset%imgmov)),' has reached energy convergence',ch10,&
&         ' with Average[Abs(Etotal(t)-Etotal(t-dt))]=',delta_energy,'<tolimg=',dtset%tolimg,ch10,&
&         '================================================================================'
       else
         write(msg,'(4a,i5,6a,es11.3,a,es11.3)') ch10,&
&         '--------------------------------------------------------------------------------',ch10,&
&         ' At time step ',itimimage,ch10,&
&         ' ',trim(imagealgo_str(dtset%imgmov)),' has reached energy convergence',ch10,&
&         ' with Average[Abs(Etotal(t)-Etotal(t-dt))]=',delta_energy,'<tolimg=',dtset%tolimg
       end if
       call wrtout(ab_out ,msg,'COLL')
       call wrtout(std_out,msg,'COLL')
       call timab(707,2,tsec)
       exit   ! exit itimimage
     end if
   end if

!Temporary statement
110 continue

!  Dont call the predictor at last time step
   if (itimimage>=ntimimage_max) call_predictor=(call_predictor.and.is_pimd)

!  Predict the next value of the images
   if (call_predictor) then
     call predictimg(delta_energy,imagealgo_str(dtset%imgmov),dtset%imgmov,itimimage,&
&     itimimage_eff,list_dynimage,ga_param,mep_param,mpi_enreg,dtset%natom,ndynimage,&
&     nimage,dtset%nimage,ntimimage_stored,pimd_param,dtset%prtvolimg,results_img_timimage)
   end if

!  Increment indexes
   if (itimimage>=ntimimage_max) exit
   itimimage_eff=itimimage_eff+1;if (itimimage_eff>ntimimage_stored) itimimage_eff=1
   if (use_hist) then
     do iimage=1,nimage
       hist(iimage)%ihist=hist(iimage)%ihist+1
     end do
   end if

   call timab(707,2,tsec)

 end do ! itimimage
!-----------------------------------------------------------------------------------------

 call timab(708,1,tsec)

!Copy the results of the computation in the appropriate arguments of the routine
 do iimage=1,nimage
   ii=mpi_enreg%my_imgtab(iimage)
   if (dtset%dynimage(ii)==1) then
     acell_img(:,iimage)     =results_img_timimage(iimage,itimimage_eff)%acell(:)
     amu_img(:,iimage)       =results_img_timimage(iimage,itimimage_eff)%amu(:)
     mixalch_img(:,:,iimage) =results_img_timimage(iimage,itimimage_eff)%mixalch(:,:)
     rprim_img(:,:,iimage)   =results_img_timimage(iimage,itimimage_eff)%rprim(:,:)
     vel_img(:,:,iimage)     =results_img_timimage(iimage,itimimage_eff)%vel(:,:)
     vel_cell_img(:,:,iimage)=results_img_timimage(iimage,itimimage_eff)%vel_cell(:,:)
     xred_img(:,:,iimage)    =results_img_timimage(iimage,itimimage_eff)%xred(:,:)
     etotal_img(iimage)      =results_img_timimage(iimage,itimimage_eff)%results_gs%etotal
     fcart_img(:,:,iimage)   =results_img_timimage(iimage,itimimage_eff)%results_gs%fcart(:,:)
     fred_img(:,:,iimage)    =results_img_timimage(iimage,itimimage_eff)%results_gs%fred(:,:)
     strten_img(:,iimage)    =results_img_timimage(iimage,itimimage_eff)%results_gs%strten(:)
   else if (compute_static_images) then
     etotal_img(iimage)    =results_img_timimage(iimage,1)%results_gs%etotal
     fcart_img(:,:,iimage) =results_img_timimage(iimage,1)%results_gs%fcart(:,:)
     fred_img(:,:,iimage)  =results_img_timimage(iimage,1)%results_gs%fred(:,:)
     strten_img(:,iimage)  =results_img_timimage(iimage,1)%results_gs%strten(:)
   end if
 end do

!When parallelizattion over images is activated, has to sum number of warnings
!and comments written in log file
 if (mpi_enreg%paral_img==1.and.mpi_enreg%comm_cell==0) then
   call specialmsg_mpisum(mpi_enreg%comm_img)
   call libpaw_spmsg_mpisum(mpi_enreg%comm_img)
 end if
!Final deallocations
 ABI_DEALLOCATE(occ)
 ABI_DEALLOCATE(vel)
 ABI_DEALLOCATE(vel_cell)
 ABI_DEALLOCATE(xred)
 ABI_DEALLOCATE(list_dynimage)
 if (allocated(amass)) then
   ABI_DEALLOCATE(amass)
 end if

 do itimimage=1,ntimimage_stored
   call destroy_results_img(results_img_timimage(:,itimimage))
 end do
 ABI_DATATYPE_DEALLOCATE(results_img_timimage)
 do iimage=1,nimage
   call scf_history_free(scf_history(iimage))
 end do
 ABI_DATATYPE_DEALLOCATE(scf_history)
 ABI_DEALLOCATE(scf_initialized)
 if (allocated(hist_prev)) then
   call abihist_free(hist_prev)
   ABI_DATATYPE_DEALLOCATE(hist_prev)
 end if
 if (allocated(hist)) then
   call abihist_free(hist)
   ABI_DATATYPE_DEALLOCATE(hist)
 end if

 call mep_destroy(mep_param)
 call ga_destroy(ga_param)
 call pimd_destroy(pimd_param)

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(708,2,tsec)
 call timab(700,2,tsec)

 DBG_EXIT("COLL")

end subroutine gstateimg
!!***
