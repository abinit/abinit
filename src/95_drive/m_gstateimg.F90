!!****m* ABINIT/m_gstateimg
!! NAME
!!  m_gstateimg
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (XG, AR, GG, MT)
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

module m_gstateimg

 use defs_basis
 use defs_wvltypes
 use defs_rectypes
 use m_abicore
 use m_abihist
 use m_mep
 use m_ga
 use m_use_ga
 use m_pimd
 use m_xmpi
 use m_errors
 use m_rec
 use m_args_gs
 use m_results_img
 use m_scf_history
 use m_io_redirect
 use m_m1geo
 use m_abimover
 use m_yaml
 use m_dtfil

 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes, only : MPI_type
 use m_time,         only : timab
 use m_geometry,     only : mkradim, mkrdim, fcart2fred, xred2xcart, metric
 use m_specialmsg,   only : specialmsg_mpisum
 use m_libpaw_tools, only : libpaw_spmsg_mpisum
 use m_pawang,       only : pawang_type
 use m_pawrad,       only : pawrad_type
 use m_pawtab,       only : pawtab_type
 use m_dtfil,        only : dtfil_init
 use m_gstate,       only : gstate
 use m_predtk,       only : prtxvf
 use m_precpred_1geo, only : precpred_1geo
 use m_pred_simple,  only : prec_simple

#if defined  HAVE_BIGDFT
 use BigDFT_API, only: mpi_environment_set
#endif

 implicit none

 private
!!***

 public :: gstateimg
!!***

contains
!!***

!!****f* ABINIT/gstateimg
!! NAME
!! gstateimg
!!
!! FUNCTION
!! Routine for conducting DFT calculations for a set of (dynamical) images
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
!!  intgres_img(nspden,natom,nimage)=gradient wrt constraints, for each image
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
!!      abihist_bcast,abihist_copy,abihist_free,abihist_init,args_gs_free
!!      args_gs_init,copy_results_img,destroy_results_img,dtfil_init,fcart2fred
!!      ga_destroy,ga_init,gstate,init_results_img,libpaw_spmsg_mpisum
!!      localfilnam,localrdfile,localredirect,localwrfile,mep_destroy,mep_init
!!      mkradim,mkrdim,pimd_destroy,pimd_init,pimd_skip_qtb,predictimg,prtimg
!!      read_md_hist_img,scf_history_free,scf_history_nullify,specialmsg_mpisum
!!      status,timab,var2hist,vel2hist,write_md_hist_img,wrtout,xmpi_barrier
!!      xmpi_sum
!!
!! SOURCE

subroutine gstateimg(acell_img,amu_img,codvsn,cpui,dtfil,dtset,etotal_img,fcart_img,&
&                    fred_img,iexit,intgres_img,mixalch_img,mpi_enreg,nimage,npwtot,occ_img,&
&                    pawang,pawrad,pawtab,psps,&
&                    rprim_img,strten_img,vel_cell_img,vel_img,wvl,xred_img,&
&                    filnam,filstat,idtset,jdtset,ndtset) ! optional arguments

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
 real(dp), intent(out) :: fred_img(3,dtset%natom,nimage)
 real(dp), intent(out) :: intgres_img(dtset%nspden,dtset%natom,nimage)
 real(dp), intent(out) :: strten_img(6,nimage)
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
 integer :: ntimimage,ntimimage_stored,ntimimage_max
 logical :: check_conv,compute_all_images,compute_static_images
 logical :: isVused,isARused,is_master,is_mep,is_pimd
 logical :: call_predictor,use_hist,use_hist_prev
 real(dp) :: delta_energy,dtion
 character(len=500) :: hist_filename,msg
 type(args_gs_type) :: args_gs
 type(mep_type) :: mep_param
 type(ga_type) :: ga_param
 type(m1geo_type) :: m1geo_param
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
&   'LINEAR COMBINATION OF CONSTRAINED DFT ENERGIES              ',& ! 6
&   '                                                            ',& ! 7
&   '                                                            ',& ! 8
&   'PATH-INTEGRAL MOLECULAR DYNAMICS (LANGEVIN)                 ',& ! 9
&   'PATH-INTEGRAL MOLECULAR DYNAMICS (QUANTUM THERMAL BATH)     ',& ! 10
&   '                                                            ',& ! 11
&   '                                                            ',& ! 12
&   'PATH-INTEGRAL MOLECULAR DYNAMICS (CHAIN OF THERMOSTATS)     '/) ! 13
 character(len=24),parameter :: stgalgo_str(0:2)=(/ &
&   'ORIGINAL ALGO.          ',& ! 0
&   'SIMPLIFIED + EQUAL ARC  ',& ! 1
&   'SIMPLIFIED + ENERGY-WGTH'/) ! 2
 character(len=20),parameter :: nebalgo_str(0:2)=(/ &
&   'ORIGINAL ALGO.      ',& ! 0
&   'IMPROVED TANGENT    ',& ! 1
&   'CLIMBING IMAGE      '/) ! 2
 character(len=20),parameter :: mepsolver_str(0:4)=(/ &
&   'STEEPEST-DESCENT    ',& ! 0
&   'QUICK-MIN OPT.      ',& ! 1
&   'L-BFGS              ',& ! 2
&   'GL-BFGS             ',& ! 3
&   'ORDER 4 RUNGE-KUTTA '/) ! 4
 real(dp) :: acell(3),rprim(3,3),rprimd(3,3),tsec(2),vel_cell(3,3)
 real(dp),allocatable :: amass(:,:),occ(:),vel(:,:),xred(:,:)
 type(abihist),allocatable :: hist(:),hist_prev(:)
 type(results_img_type),pointer :: results_img(:,:),res_img(:)
 type(scf_history_type),allocatable :: scf_history(:)
 !type(abiforstr) :: preconforstr ! Preconditioned forces and stress ... Only needed to deallocate an internal matrix in prec_simple

! ***********************************************************************

 DBG_ENTER("COLL")

 call timab(700,1,tsec)
 call timab(703,3,tsec)

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
&       imgtab=mpi_enreg%my_imgtab)
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
 ABI_ALLOCATE(xred,(3,dtset%natom))
 ABI_DATATYPE_ALLOCATE(results_img,(nimage,ntimimage_stored))
 ABI_ALLOCATE(list_dynimage,(dtset%ndynimage))
 do itimimage=1,ntimimage_stored
   res_img => results_img(:,itimimage)
   call init_results_img(dtset%natom,dtset%npspalch,dtset%nspden,dtset%nsppol,dtset%ntypalch,&
&   dtset%ntypat,res_img)
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
&   (abs(dtset%densfor_pred)==5.or.abs(dtset%densfor_pred)==6)) then
      history_size=2
      if(dtset%extrapwf==2) history_size=3
    end if
 else
   if (abs(dtset%densfor_pred)==2.or.abs(dtset%densfor_pred)==3) history_size=0
   if (dtset%imgwfstor==1) history_size=1
   if (dtset%usewvl==0.and.(abs(dtset%densfor_pred)==5.or.abs(dtset%densfor_pred)==6)) history_size=2
 end if
 do iimage=1,nimage
   call scf_history_nullify(scf_history(iimage))
   scf_history(iimage)%history_size=history_size
 end do

!In some cases, need amass variable
 if (use_hist) then
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

!MEP search: fill in eventually the data structure mep_param
 call mep_init(dtset,mep_param)

!GA search: fill in eventually the data structure ga_param
 call ga_init(dtset,ga_param)

!Move 1GEO approach: fill the data structure m1geo_param
 call m1geo_init(dtfil,dtset,m1geo_param)

!PIMD: fill in eventually the data structure pimd_param
 call pimd_init(dtset,pimd_param,is_master)
 dtion=one;if (is_pimd) dtion=pimd_param%dtion

 call timab(703,2,tsec)

!-----------------------------------------------------------------------------------------
!Big loop on the propagation of all images
 itimimage_eff=1
 do itimimage=1,ntimimage

   res_img => results_img(:,itimimage_eff)
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
&         hist_prev(iimage)%rprimd(:,:,ih),dtset%natom)
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
     if (dtset%imgmov==2) then
       write(msg,'(6a)') trim(msg),' (',trim(stgalgo_str(mep_param%string_algo)),' + ',&
&                        trim(mepsolver_str(mep_param%mep_solver)),')'
     end if
     if (dtset%imgmov==5) then
       ii=merge(mep_param%neb_algo,1,mep_param%neb_algo/=2.or.itimimage>=mep_param%cineb_start)
       write(msg,'(6a)') trim(msg),' (',trim(nebalgo_str(ii)),' + ',&
&                        trim(mepsolver_str(mep_param%mep_solver)),')'
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
   call yaml_iterstart('timimage', itimimage, ab_out, dtset%use_yaml)

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
           '--------------------------------------------------------------------------------',ch10,&
           ' ',trim(imagealgo_str(dtset%imgmov)),' - CELL # ',ii,'/',dtset%nimage,ch10,&
           '--------------------------------------------------------------------------------',ch10
         if (dtset%prtvolimg==0) call wrtout(ab_out ,msg,'COLL')
         if (do_write_log) call wrtout(std_out,msg,'PERS')
       end if
       call yaml_iterstart('image', iimage, ab_out, dtset%use_yaml)

       acell(:)     =res_img(iimage)%acell(:)
       rprim(:,:)   =res_img(iimage)%rprim(:,:)
       vel(:,:)     =res_img(iimage)%vel(:,:)
       vel_cell(:,:)=res_img(iimage)%vel_cell(:,:)
       xred(:,:)    =res_img(iimage)%xred(:,:)
       occ(:)       =occ_img(:,iimage)

       call args_gs_init(args_gs, &
&       res_img(iimage)%amu(:),res_img(iimage)%mixalch(:,:),&
&       dtset%dmatpawu(:,:,:,:,ii),dtset%upawu(:,ii),dtset%jpawu(:,ii),&
&       dtset%rprimd_orig(:,:,ii))

       call timab(705,2,tsec)

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

!    check change of rprim and rewriting in hist
!    check change of xred and rewriting in hist

!      Close output units ; restore defaults
       call localredirect(mpi_enreg%comm_cell,mpi_enreg%comm_world,dtset%nimage,mpi_enreg%paral_img,dtset%prtvolimg)
       call timab(706,2,tsec)

     else if (itimimage>1) then ! For static images, simply copy one time step to the other
       itimimage_prev=itimimage_eff-1
       if (itimimage_prev<1) itimimage_prev=ntimimage_stored
       call copy_results_img(results_img(iimage,itimimage_prev), &
&       results_img(iimage,itimimage_eff ))
     end if

!    Store results in hist datastructure
     if (use_hist) then
       ih=hist(iimage)%ihist
       call mkrdim(res_img(iimage)%acell(:),res_img(iimage)%rprim(:,:),rprimd)
       call var2hist(res_img(iimage)%acell(:),hist(iimage),dtset%natom,&
&       rprimd,res_img(iimage)%xred(:,:),.FALSE.)
       call vel2hist(amass(:,iimage),hist(iimage),res_img(iimage)%vel(:,:),&
&       res_img(iimage)%vel_cell(:,:))
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
     call write_md_hist_img(hist,hist_filename,ifirst,itimimage,dtset%natom,dtset%ntypat,&
&     dtset%typat,amu_img(:,1),dtset%znucl,dtion,&
&     nimage=dtset%nimage,imgmov=dtset%imgmov,mdtemp=dtset%mdtemp,comm_img=mpi_enreg%comm_img,&
&     imgtab=mpi_enreg%my_imgtab)
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
       _IBM6("hello world")
       iimage=list_dynimage(idynimage)
       delta_energy=delta_energy &
&       +abs(results_img(iimage,itimimage)%results_gs%etotal &
&       -results_img(iimage,itimimage-idelta)%results_gs%etotal)
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
       call wrtout([std_out, ab_out] ,msg,'COLL')
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
&     itimimage_eff,list_dynimage,ga_param,mep_param,mpi_enreg,m1geo_param,dtset%natom,ndynimage,&
&     nimage,dtset%nimage,ntimimage_stored,pimd_param,dtset%prtvolimg,results_img)
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
     acell_img(:,iimage)     =results_img(iimage,itimimage_eff)%acell(:)
     amu_img(:,iimage)       =results_img(iimage,itimimage_eff)%amu(:)
     mixalch_img(:,:,iimage) =results_img(iimage,itimimage_eff)%mixalch(:,:)
     rprim_img(:,:,iimage)   =results_img(iimage,itimimage_eff)%rprim(:,:)
     vel_img(:,:,iimage)     =results_img(iimage,itimimage_eff)%vel(:,:)
     vel_cell_img(:,:,iimage)=results_img(iimage,itimimage_eff)%vel_cell(:,:)
     xred_img(:,:,iimage)    =results_img(iimage,itimimage_eff)%xred(:,:)
     etotal_img(iimage)      =results_img(iimage,itimimage_eff)%results_gs%etotal
     fcart_img(:,:,iimage)   =results_img(iimage,itimimage_eff)%results_gs%fcart(:,:)
     fred_img(:,:,iimage)    =results_img(iimage,itimimage_eff)%results_gs%fred(:,:)
     intgres_img(:,:,iimage) =results_img(iimage,itimimage_eff)%results_gs%intgres(:,:)
     strten_img(:,iimage)    =results_img(iimage,itimimage_eff)%results_gs%strten(:)
   else if (compute_static_images) then
     etotal_img(iimage)    =results_img(iimage,1)%results_gs%etotal
     fcart_img(:,:,iimage) =results_img(iimage,1)%results_gs%fcart(:,:)
     fred_img(:,:,iimage)  =results_img(iimage,1)%results_gs%fred(:,:)
     intgres_img(:,:,iimage)=results_img(iimage,1)%results_gs%intgres(:,:)
     strten_img(:,iimage)  =results_img(iimage,1)%results_gs%strten(:)
   end if
 end do

!When parallelizattion over images is activated, has to sum number of warnings
!and comments written in log file
 if (mpi_enreg%paral_img==1.and.mpi_enreg%comm_cell==0) then
   call specialmsg_mpisum(mpi_enreg%comm_img)
   call libpaw_spmsg_mpisum(mpi_enreg%comm_img)
 end if


!Final deallocations

!This call is needed to free internal storages in different routines (prec_simple, pred_bfgs ...)
 if(dtset%imgmov==6)then
   m1geo_param%iexit=1
   call predictimg(delta_energy,imagealgo_str(dtset%imgmov),dtset%imgmov,itimimage,&
&   itimimage_eff,list_dynimage,ga_param,mep_param,mpi_enreg,m1geo_param,dtset%natom,ndynimage,&
&   nimage,dtset%nimage,ntimimage_stored,pimd_param,dtset%prtvolimg,results_img)
 endif

 ABI_DEALLOCATE(occ)
 ABI_DEALLOCATE(vel)
 ABI_DEALLOCATE(xred)
 ABI_DEALLOCATE(list_dynimage)

 if (allocated(amass)) then
   ABI_DEALLOCATE(amass)
 end if

 do itimimage=1,ntimimage_stored
   call destroy_results_img(results_img(:,itimimage))
 end do
 ABI_DATATYPE_DEALLOCATE(results_img)
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
 call m1geo_destroy(m1geo_param)
 call pimd_destroy(pimd_param)

 call timab(708,2,tsec)
 call timab(700,2,tsec)

 DBG_EXIT("COLL")

end subroutine gstateimg
!!***

!!****f* ABINIT/prtimg
!! NAME
!! prtimg
!!
!! FUNCTION
!! Print out results obtained by as ground-state calculation of
!! an image of the cell. The printing format is condensed in order
!! to facilitate the reading.
!!
!! INPUTS
!!  dynimage(nimagetot)=flags defining static/dynamic state of images
!!  imagealgo_str=name of the algorithm (with images) used
!!  imgmov=index of algorithm (with images) used
!!  iout=unit number for output
!!  mpi_enreg=MPI-parallelisation information
!!  nimage=number of images stored on current proc
!!  nimage_tot=total number of images (should be dtset%nimage)
!!  prt_all_images=true if all images have to be printed out (ignoring dynimage)
!!  prtvolimg=printing volume for each image
!!           <0 : nothing
!!            0 : only a title
!!            1 : energy, residuals, forces, stresses, velocities, atomic positions
!!            2 : energy, residuals
!!  resimg(nimage) <type(results_img_type)>=results of the ground-state computations
!!                                          for all images treated by current proc
!!
!! OUTPUT
!!  (data written to unit iout)
!!
!! PARENTS
!!      gstateimg
!!
!! CHILDREN
!!      destroy_results_img,gather_results_img,metric,prtxvf,wrtout,xred2xcart
!!
!! SOURCE

subroutine prtimg(dynimage,imagealgo_str,imgmov,iout,mpi_enreg,nimage,nimage_tot,&
&                 prt_all_images,prtvolimg,resimg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nimage_tot,dynimage(nimage_tot),imgmov,iout,nimage,prtvolimg !vz_d
 logical,intent(in) :: prt_all_images
 character(len=60),intent(in) :: imagealgo_str
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 type(results_img_type),target,intent(inout) :: resimg(nimage)

!Local variables-------------------------------
!scalars
 integer :: ii,prtvel
 logical :: test_img
 real(dp) :: ucvol_img
 character(len=500) :: msg
!arrays
 integer,allocatable :: iatfix_img(:,:)
 real(dp),allocatable :: gmet_img(:,:),gprimd_img(:,:),rmet_img(:,:),xcart_img(:,:)
 type(results_img_type),pointer :: resimg_all(:)

! ****************************************************************

 DBG_ENTER('COLL')

 if (prtvolimg<=0) return
 if (mpi_enreg%me_cell/=0) return

!Gather data
 if (prtvolimg==1.or.prtvolimg==2) then
   test_img=(nimage_tot/=1.and.mpi_enreg%paral_img==1)
   if (test_img) then
     if (mpi_enreg%me==0)  then
       ABI_DATATYPE_ALLOCATE(resimg_all,(nimage_tot))
     end if
     call gather_results_img(mpi_enreg,resimg,resimg_all,master=0,&
&     allgather=.false.,only_one_per_img=.true.)
   else
     resimg_all => resimg
   end if
 end if

!===== First option for the printing volume ===
 if (prtvolimg==1.and.mpi_enreg%me==0) then

   prtvel=0;if (imgmov==0.or.imgmov==9.or.imgmov==10.or.imgmov==13) prtvel=1

   do ii=1,nimage_tot
     if (dynimage(ii)==1.or.prt_all_images) then

!      Title
       write(msg,'(6a,i4,a,i4,2a)') ch10,&
         '----------------------------------------------------------------------',ch10,&
         ' ',trim(imagealgo_str),' - CELL # ',ii,'/',nimage_tot,ch10,&
         '----------------------------------------------------------------------'
       call wrtout(iout,msg,'COLL')

!      Total energy
       write(msg,'(2a,es20.12)') ch10,' Total energy for the cell [Ha]: ',resimg_all(ii)%results_gs%etotal
       call wrtout(iout,msg,'COLL')

!      Residuals of the SCF cycle
       write(msg,'(3a,4(a,es16.8,a))') ch10,&
         ' Residuals from SCF cycle: ',ch10,&
         '    Total energy difference        =',resimg_all(ii)%results_gs%deltae,ch10,&
         '    Maximal forces difference      =',resimg_all(ii)%results_gs%diffor,ch10,&
         '    Max. residual of wave-functions=',resimg_all(ii)%results_gs%residm,ch10,&
         '    Density/potential residual (^2)=',resimg_all(ii)%results_gs%res2,ch10
       call wrtout(iout,msg,'COLL')

!      Cell parameters
       ABI_ALLOCATE(rmet_img,(3,3))
       ABI_ALLOCATE(gmet_img,(3,3))
       ABI_ALLOCATE(gprimd_img,(3,3))
       call metric(gmet_img,gprimd_img,iout,rmet_img,resimg_all(ii)%rprim,ucvol_img)
       ABI_DEALLOCATE(rmet_img)
       ABI_DEALLOCATE(gmet_img)
       ABI_DEALLOCATE(gprimd_img)

!      Positions, forces and velocities
       ABI_ALLOCATE(iatfix_img,(3,resimg_all(ii)%natom))
       ABI_ALLOCATE(xcart_img,(3,resimg_all(ii)%natom))
       iatfix_img=0
       call xred2xcart(resimg_all(ii)%natom,resimg_all(ii)%rprim,xcart_img,resimg_all(ii)%xred)
       call prtxvf(resimg_all(ii)%results_gs%fcart,resimg_all(ii)%results_gs%fred,&
&       iatfix_img,iout,resimg_all(ii)%natom,prtvel,&
&       resimg_all(ii)%vel,xcart_img,resimg_all(ii)%xred)
       ABI_DEALLOCATE(iatfix_img)
       ABI_DEALLOCATE(xcart_img)

!      Stress tensor
       write(msg, '(a,es12.4,a)' ) &
&       '-Cartesian components of stress tensor (GPa)         [Pressure=',&
&       -(resimg_all(ii)%results_gs%strten(1)+resimg_all(ii)%results_gs%strten(2) &
&       +resimg_all(ii)%results_gs%strten(3))*HaBohr3_GPa/three,' GPa]'
       call wrtout(iout,msg,'COLL')
       write(msg, '(2(a,1p,e16.8))' ) '- sigma(1 1)=',resimg_all(ii)%results_gs%strten(1)*HaBohr3_GPa,&
&       '  sigma(3 2)=',resimg_all(ii)%results_gs%strten(4)*HaBohr3_GPa
       call wrtout(iout,msg,'COLL')
       write(msg, '(2(a,1p,e16.8))' ) '- sigma(2 2)=',resimg_all(ii)%results_gs%strten(2)*HaBohr3_GPa,&
&       '  sigma(3 1)=',resimg_all(ii)%results_gs%strten(5)*HaBohr3_GPa
       call wrtout(iout,msg,'COLL')
       write(msg, '(2(a,1p,e16.8))' ) '- sigma(3 3)=',resimg_all(ii)%results_gs%strten(3)*HaBohr3_GPa,&
&       '  sigma(2 1)=',resimg_all(ii)%results_gs%strten(6)*HaBohr3_GPa
       call wrtout(iout,msg,'COLL')
     end if
   end do
 end if


!===== 2nd option for the printing volume ===
 if (prtvolimg==2.and.mpi_enreg%me==0) then
   write(msg,'(a,1x,a)') ch10,'Cell   Total_energy[Ha]     deltae       diffor       residm         res2'
   call wrtout(iout,msg,'COLL')
   do ii=1,nimage_tot
     if (dynimage(ii)==1.or.prt_all_images) then
       write(msg,'(1x,i4,2x,es16.8,4(1x,es13.5))') &
&       ii,resimg_all(ii)%results_gs%etotal,resimg_all(ii)%results_gs%deltae,&
&       resimg_all(ii)%results_gs%diffor,resimg_all(ii)%results_gs%residm,&
&       resimg_all(ii)%results_gs%res2
       call wrtout(iout,msg,'COLL')
     end if
   end do
 end if

!=====
 if (prtvolimg==1.or.prtvolimg==2) then
   if (test_img.and.mpi_enreg%me==0) then
     call destroy_results_img(resimg_all)
     ABI_DATATYPE_DEALLOCATE(resimg_all)
   end if
   nullify(resimg_all)
 end if

 DBG_EXIT('COLL')

end subroutine prtimg
!!***

!!****f* ABINIT/predictimg
!! NAME
!! predictimg
!!
!! FUNCTION
!! Given the past history of images, predict the new set of images
!!
!! INPUTS
!! deltae=averaged energy difference used to control convergence over images
!! imagealgo_str=name of the algorithm (with images) used
!! imgmov=gives the algorithm to be used for prediction of new set of images
!! itimimage=time index for image propagation (itimimage+1 is to be predicted here)
!! itimimage_eff=time index in the history
!! list_dynimage(nimage)=list of dynamical images. The non-dynamical ones will not change.
!!       Example : in the NEB method, or in the string method, one expect the two end images to be fixed.
!! mep_param=several parameters for Minimal Energy Path (MEP) search
!! mpi_enreg=MPI-parallelisation information
!! natom= number of atoms
!! ndynimage=number of dynamical images
!! nimage=number of images (treated by current proc)
!! nimage_tot=total number of images
!! ntimimage_stored=number of time steps stored in the history
!! pimd_param=several parameters for Path-Integral MD
!! prtvolimg=printing volume
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! results_img(ntimimage_stored,nimage)=datastructure that holds the history of previous computations.
!!   results_img(:,:)%acell(3)
!!    at input, history of the values of acell for all images
!!    at output, the predicted values of acell for all images
!!   results_img(:,:)%results_gs
!!    at input, history of the values of energies and forces for all images
!!   results_img(:,:)%rprim(3,3)
!!    at input, history of the values of rprim for all images
!!    at output, the predicted values of rprim for all images
!!   results_img(:,:)%vel(3,natom)
!!    at input, history of the values of vel for all images
!!    at output, the predicted values of vel for all images
!!   results_img(:,:)%vel_cell(3,3)
!!    at input, history of the values of vel_cell for all images
!!    at output, the predicted values of vel_cell for all images
!!   results_img(:,:)%xred(3,natom)
!!    at input, history of the values of xred for all images
!!    at output, the predicted values of xred for all images
!!
!! PARENTS
!!      gstateimg
!!
!! CHILDREN
!!      predict_copy,predict_ga,predict_neb,predict_pimd,predict_steepest
!!      predict_string,wrtout
!!
!! SOURCE

subroutine predictimg(deltae,imagealgo_str,imgmov,itimimage,itimimage_eff,list_dynimage,&
&                     ga_param,mep_param,mpi_enreg,m1geo_param,natom,ndynimage,nimage,nimage_tot,&
&                     ntimimage_stored,pimd_param,prtvolimg,results_img)

 use m_results_gs , only : results_gs_type
 use m_predict_neb, only : predict_neb
 use m_predict_steepest, only : predict_steepest
 use m_predict_pimd,    only : predict_pimd
 use m_predict_string, only : predict_string

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: imgmov,itimimage,itimimage_eff,natom,ndynimage
 integer,intent(in) :: nimage,nimage_tot,ntimimage_stored,prtvolimg
 character(len=60),intent(in) :: imagealgo_str
 real(dp),intent(in) :: deltae
 type(mep_type),intent(inout) :: mep_param
 type(m1geo_type),intent(inout) :: m1geo_param
 type(ga_type),intent(inout) :: ga_param
 type(pimd_type),intent(inout) :: pimd_param
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: list_dynimage(ndynimage)
 type(results_img_type) :: results_img(nimage,ntimimage_stored)

!Local variables-------------------------------
!scalars
 integer,save :: idum=5
 logical :: is_pimd
 character(len=500) :: msg

! *************************************************************************

 is_pimd=(imgmov==9.or.imgmov==10.or.imgmov==13)

!Write convergence info
 write(msg,'(3a)') ch10,&
& '------------------------------------------------------------',ch10
 if (prtvolimg<2) write(msg,'(5a)') trim(msg),' ',trim(imagealgo_str),':',ch10

!Specific case of 4th-order RK algorithm
 if (mep_param%mep_solver==4) then
   if (mod(itimimage,4)==0) then
     write(msg,'(4a)') trim(msg),&
&     ' Fourth-order Runge-Kutta algorithm - final step',ch10
     if (itimimage>4) write(msg,'(2a,es11.3,2a)') trim(msg),&
&     ' Average[Abs(Etotal(t)-Etotal(t-dt))]=',deltae,' Hartree',ch10
     write(msg,'(2a)') trim(msg),' Moving images of the cell...'
   else
     write(msg,'(2a,i1,2a)') trim(msg),&
&     ' Fourth-order Runge-Kutta algorithm - intermediate step ',mod(itimimage,4),ch10
     write(msg,'(2a)') trim(msg),' Computing new intermediate positions...'
   end if
 else if (is_pimd) then

!  PIMD
   write(msg,'(2a)') trim(msg),' Moving images of the cell...'
 else

!  Other cases
   if (itimimage>1) write(msg,'(2a,es11.3,2a)') trim(msg),&
&   ' Average[Abs(Etotal(t)-Etotal(t-dt))]=',deltae,' Hartree',ch10
   write(msg,'(2a)') trim(msg),' Moving images of the cell...'

 end if

!Prevent writing if iexit==1, which at present only happens for imgmov==6 algo
 if(imgmov/=6 .or. m1geo_param%iexit==0) call wrtout([std_out, ab_out] ,msg)

 select case(imgmov)

 case(0)
   call predict_copy(itimimage_eff,list_dynimage,ndynimage,nimage,&
&   ntimimage_stored,results_img)

 case(1)
   call predict_steepest(itimimage,itimimage_eff,list_dynimage,mep_param,natom,ndynimage,nimage,&
&   ntimimage_stored,results_img)

 case(2)
   call predict_string(itimimage,itimimage_eff,list_dynimage,mep_param,mpi_enreg,natom,&
&   ndynimage,nimage,nimage_tot,ntimimage_stored,results_img)

 case(4)
   call predict_ga(itimimage_eff,idum,ga_param,natom,nimage,&
&   ntimimage_stored,results_img)

 case(5)
   call predict_neb(itimimage,itimimage_eff,list_dynimage,mep_param,mpi_enreg,natom,&
&   ndynimage,nimage,nimage_tot,ntimimage_stored,results_img)

 case(6)
   call move_1geo(itimimage_eff,m1geo_param,mpi_enreg,nimage,ntimimage_stored,results_img)

 case(9, 10, 13)
!    Path Integral Molecular Dynamics
   call predict_pimd(imgmov,itimimage,itimimage_eff,mpi_enreg,natom,nimage,nimage_tot,&
&   ntimimage_stored,pimd_param,prtvolimg,results_img)

 case default

 end select

end subroutine predictimg
!!***

!!****f* ABINIT/predict_copy
!! NAME
!! predict_copy
!!
!! FUNCTION
!! Given the past history of images, predict the new set of images.
!! Here, simple copy of the previous image.
!!
!! INPUTS
!! itimimage_eff=time index in the history
!! list_dynimage(nimage)=list of dynamical images. The non-dynamical ones will not change.
!!       Example : in the NEB of string method, one expect the two end images to be fixed.
!! ndynimage=number of dynamical images
!! nimage=number of images
!! ntimimage_stored=number of time steps stored in the history
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! results_img(ntimimage_stored,nimage)=datastructure that holds the history of previous computations.
!!   results_img(:,:)%acell(3)
!!    at input, history of the values of acell for all images
!!    at output, the predicted values of acell for all images
!!   results_img(:,:)%results_gs
!!    at input, history of the values of energies and forces for all images
!!   results_img(:,:)%rprim(3,3)
!!    at input, history of the values of rprim for all images
!!    at output, the predicted values of rprim for all images
!!   results_img(:,:)%vel(3,natom)
!!    at input, history of the values of vel for all images
!!    at output, the predicted values of vel for all images
!!   results_img(:,:)%vel_cell(3,3)
!!    at input, history of the values of vel_cell for all images
!!    at output, the predicted values of vel_cell for all images
!!   results_img(:,:)%xred(3,natom)
!!    at input, history of the values of xred for all images
!!    at output, the predicted values of xred for all images
!!
!! PARENTS
!!      predictimg
!!
!! CHILDREN
!!
!! SOURCE

subroutine predict_copy(itimimage_eff,list_dynimage,ndynimage,nimage,&
&                       ntimimage_stored,results_img)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itimimage_eff,ndynimage,nimage,ntimimage_stored
!arrays
 integer,intent(in) :: list_dynimage(ndynimage)
 type(results_img_type),intent(inout) :: results_img(nimage,ntimimage_stored)

!Local variables-------------------------------
!scalars
 integer :: idynimage,iimage,next_itimimage

! *************************************************************************

 next_itimimage=itimimage_eff+1
 if (next_itimimage>ntimimage_stored) next_itimimage=1

 do idynimage=1,ndynimage

   iimage=list_dynimage(idynimage)

   results_img(iimage,next_itimimage)%acell(:)     =results_img(iimage,itimimage_eff)%acell(:)
   results_img(iimage,next_itimimage)%rprim(:,:)   =results_img(iimage,itimimage_eff)%rprim(:,:)
   results_img(iimage,next_itimimage)%vel(:,:)     =results_img(iimage,itimimage_eff)%vel(:,:)
   results_img(iimage,next_itimimage)%vel_cell(:,:)=results_img(iimage,itimimage_eff)%vel_cell(:,:)
   results_img(iimage,next_itimimage)%xred(:,:)    =results_img(iimage,itimimage_eff)%xred(:,:)

 end do  ! idynimage

end subroutine predict_copy
!!***

!!****f* ABINIT/move_1geo
!! NAME
!! move_1geo
!!
!! FUNCTION
!! This subroutine uses the forces, stresses and other results obtained for several images with one, common, geometry,
!! weight them to deliver averaged forces, stresses, etc, and uses these to predict the next common geometry.
!! All images must be dynamical.
!! WARNING : at present, only forces are used, to change atomic positions. No change of cell geometry.
!! Since this is not the PIMD algorithm, suppose ntimimage_stored=ntimimage, and itimimage=itimimage_eff.
!!
!! INPUTS
!! itimimage_eff=time index in the history
!! nimage=number of images
!! ntimimage_stored=number of time steps stored in the history
!!  mpi_enreg=MPI-parallelisation information
!! m1geo_param=parameters for the 1geo algorithms
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! results_img(ntimimage_stored,nimage)=datastructure that holds the history of previous computations.
!!   results_img(:,:)%acell(3)
!!    at input, history of the values of acell for all images
!!    at output, the predicted values of acell for all images
!!   results_img(:,:)%results_gs
!!    at input, history of the values of energies and forces for all images
!!   results_img(:,:)%rprim(3,3)
!!    at input, history of the values of rprim for all images
!!    at output, the predicted values of rprim for all images
!!   results_img(:,:)%vel(3,natom)
!!    at input, history of the values of vel for all images
!!    at output, the predicted values of vel for all images
!!   results_img(:,:)%vel_cell(3,3)
!!    at input, history of the values of vel_cell for all images
!!    at output, the predicted values of vel_cell for all images
!!   results_img(:,:)%xred(3,natom)
!!    at input, history of the values of xred for all images
!!    at output, the predicted values of xred for all images
!!
!! PARENTS
!!      predictimg
!!
!! CHILDREN
!!
!! SOURCE

subroutine move_1geo(itimimage_eff,m1geo_param,mpi_enreg,nimage,ntimimage_stored,results_img)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itimimage_eff,nimage,ntimimage_stored
 type(MPI_type),intent(in) :: mpi_enreg
 type(m1geo_type),intent(inout) :: m1geo_param
!arrays
 type(results_img_type),intent(inout) :: results_img(nimage,ntimimage_stored)

!Local variables-------------------------------
!scalars
 integer :: ihist,iimage,natom,next_itimimage
 real(dp) :: acell(3),rprim(3,3),rprimd(3,3),strten(6),vel_cell(3,3)
 real(dp),allocatable :: fcart(:,:),vel(:,:),xred(:,:)
 logical :: DEBUG=.FALSE.

! *************************************************************************

 natom=m1geo_param%ab_mover%natom
 ihist=m1geo_param%hist_1geo%ihist

 ABI_ALLOCATE(fcart,(3,natom))
 ABI_ALLOCATE(vel,(3,natom))
 ABI_ALLOCATE(xred,(3,natom))

!Of course, assume that the geometry parameters are the same for all images, so take them from the first one.
 xred(:,:)    =results_img(1,itimimage_eff)%xred(:,:)
 acell(:)     =results_img(1,itimimage_eff)%acell(:)
 rprim(:,:)   =results_img(1,itimimage_eff)%rprim(:,:)
 vel(:,:)     =results_img(1,itimimage_eff)%vel(:,:)
 vel_cell(:,:)=results_img(1,itimimage_eff)%vel_cell(:,:)

 call mkrdim(acell,rprim,rprimd)

!Fill history with the values of xred, acell and rprimd
 call var2hist(acell,m1geo_param%hist_1geo,natom,rprimd,xred,DEBUG)

!Fill history with velocities and ionic kinetic energy
 call vel2hist(m1geo_param%ab_mover%amass,m1geo_param%hist_1geo,vel,vel_cell)
 m1geo_param%hist_1geo%time(ihist)=zero

!Compute forces and stresses for the 1geo : take the weighted average.
 fcart(:,:)=zero
 strten(:)=zero
 do iimage=1,nimage
   fcart(:,:)=fcart(:,:)+results_img(iimage,itimimage_eff)%results_gs%fcart(:,:)*m1geo_param%mixesimgf(iimage)
   strten(:) =strten(:) +results_img(iimage,itimimage_eff)%results_gs%strten(:)*m1geo_param%mixesimgf(iimage)
 enddo

!Store them in hist_1geo
 m1geo_param%hist_1geo%fcart(:,:,ihist)=fcart(:,:)
 m1geo_param%hist_1geo%strten(:,ihist) =strten(:)

!Store them in ab_xfh
!THIS IS TO BE DONE !

!Compute new atomic positions and cell characteristics in the single geometry
 call precpred_1geo(m1geo_param%ab_mover,&
& m1geo_param%ab_xfh_1geo,&
& m1geo_param%ab_mover%amu_curr,&
& m1geo_param%deloc,&
& m1geo_param%dt_chkdilatmx,&
& mpi_enreg%comm_cell,&
& m1geo_param%dilatmx,&
& m1geo_param%filnam_ds4,&
& m1geo_param%hist_1geo,&
& m1geo_param%hmctt,&
& m1geo_param%icycle,&
& m1geo_param%iexit,&
!& m1geo_param%itime,&
  itimimage_eff,&       ! m1geo_param%itime should be eliminated, no need for it
& m1geo_param%mttk_vars,&
& m1geo_param%nctime,&
& m1geo_param%ncycle,&
& m1geo_param%nerr_dilatmx,&
& m1geo_param%npsp,&
& m1geo_param%ntime,&
& m1geo_param%rprimd_orig,&
& m1geo_param%skipcycle,&
& m1geo_param%usewvl)

!Retrieve the new positions, cell parameters [and velocities ?!]
 call hist2var(acell,m1geo_param%hist_1geo,natom,rprimd,xred,DEBUG)

!Store acell, rprim, xred and vel for the new iteration
 next_itimimage=itimimage_eff+1
 if (next_itimimage>ntimimage_stored)then
   MSG_ERROR('next_itimimage>ntimimage_stored')
 endif

 do iimage=1,nimage
   results_img(iimage,next_itimimage)%xred(:,:)    =xred(:,:)
   results_img(iimage,next_itimimage)%acell(:)     =acell(:)
   results_img(iimage,next_itimimage)%rprim(:,:)   =rprim(:,:)
!  WARNING : Should also store vel and vel_cell of course ...
!  results_img(iimage,next_itimimage)%vel(:,:)     =vel(:,:)
!  results_img(iimage,next_itimimage)%vel_cell(:,:)=vel_cell(:,:)
 end do

 ABI_DEALLOCATE(fcart)
 ABI_DEALLOCATE(vel)
 ABI_DEALLOCATE(xred)

end subroutine move_1geo
!!***

end module m_gstateimg
!!***
