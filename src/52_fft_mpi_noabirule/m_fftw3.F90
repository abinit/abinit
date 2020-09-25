!!****m* ABINIT/m_fftw3
!! NAME
!! m_fftw3
!!
!! FUNCTION
!!  This module provides wrappers for the FFTW3 routines: in-place and out-of-place version.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2020 ABINIT group (MG, FD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  1) MPI parallelism is in testing stage
!!  2) For better performance the FFT divisions should contain small factors  [2, 3, 5, 7, 11]
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

! It seems that MKL wrappers do not like the advanced interfaces for
! r2c and c2r transforms although they work fine if the true FFTW3 library is used.
!#define DEV_RC_BUG
#undef DEV_RC_BUG

#define FFTLIB "FFTW3"
#define FFT_PREF(name) CONCAT(fftw3_,name)
#define SPAWN_THREADS_HERE(ndat, nthreads) fftw3_spawn_threads_here(ndat, nthreads)
#define FFT_DOUBLE 1
#define FFT_SINGLE 2
#define FFT_MIXPREC 3

MODULE m_fftw3

 use defs_basis
 use m_abicore
 use m_errors
 use m_xomp
 use m_xmpi
 use m_hide_blas
 use m_cgtools
 use m_cplxtools
 use m_distribfft
 use m_fftcore
 use iso_c_binding

 use m_time,           only : timab
 use m_numeric_tools,  only : imax_loc
 use defs_abitypes,    only : MPI_type
 use m_mpinfo,         only : ptabs_fourwf
 use m_fstrings,       only : strcat
 use m_fft_mesh,       only : zpad_t, zpad_init, zpad_free

 implicit none

#ifdef HAVE_FFTW3_MPI
 include 'fftw3-mpi.f03'
#endif

!This should be done but MKL fftw hasn't always this include file
!#ifdef HAVE_FFT_FFTW3
! include 'fftw3.f03'
!#endif

 private

! Entry points for client code
 public :: fftw3_seqfourdp      ! 3D FFT of lengths nx, ny, nz. Mainly used for densities or potentials.
 public :: fftw3_seqfourwf      ! FFT transform of wavefunctions (high-level interface).
 public :: fftw3_fftrisc
 public :: fftw3_fftug          ! G-->R. 3D zero-padded FFT of lengths nx, ny, nz. Mainly used for wavefunctions
 public :: fftw3_fftur          ! R-->G, 3D zero-padded FFT of lengths nx, ny, nz. Mainly used for wavefunctions
 public :: fftw3_use_lib_threads

 public :: fftw3_mpifourdp

! Low-level routines.
 public :: fftw3_cleanup        ! Reset FFTW to the pristine state it was in when you started your program,
 public :: fftw3_init_threads   ! one-time initialization required to use FFTW3 threads.
 public :: fftw3_set_nthreads   ! Set the number of threads you want FFTW3 to use when HAVE_FFT_FFTW3_THREADS is defined.
 public :: fftw3_r2c_op         ! Real to complex transform (out-of-place version).
 public :: fftw3_c2r_op         ! Complex to real transform (out-of-place version).
 public :: fftw3_c2c_op         ! complex to complex transform (out-of-place version).
 public :: fftw3_c2c_ip         ! complex to complex transform (in-place version).
 public :: fftw3_many_dft_op    ! Driver routine for many out-of-place 3D complex-to-complex FFTs.
 public :: fftw3_many_dft_ip    ! Driver routine for many in-place 3D complex-to-complex FFTs.
 public :: fftw3_fftpad         ! Driver routines for zero-padded FFT of wavefunctions.
 public :: fftw3_fftpad_dp      ! Driver routines for zero-padded FFT of wavefunctions.
 public :: fftw3_fftug_dp       ! Driver routines for zero-padded FFT of wavefunctions.
 public :: fftw3_poisson        ! Solve the poisson equation in G-space starting from n(r).

 ! MPI version
 public :: fftw3_mpiback_wf
 public :: fftw3_mpiback_manywf
 public :: fftw3_mpiforw_wf
 public :: fftw3_mpiforw_manywf
 public :: fftw3_mpiback
 public :: fftw3_mpiforw
 public :: fftw3_applypot
 public :: fftw3_applypot_many
 public :: fftw3_accrho

#ifdef HAVE_FFTW3_MPI
! flags copied from fftw3.f
 integer,public,parameter :: ABI_FFTW_FORWARD = FFTW_FORWARD
 integer,public,parameter :: ABI_FFTW_BACKWARD = FFTW_BACKWARD
 integer,public,parameter :: ABI_FFTW_ESTIMATE = FFTW_ESTIMATE
 integer,public,parameter :: ABI_FFTW_MEASURE = FFTW_MEASURE
 ! end flags copied from fftw3.f
 integer,public,parameter :: ABI_FFTW_MPI_TRANSPOSED_IN = FFTW_MPI_TRANSPOSED_IN
 integer,public,parameter :: ABI_FFTW_MPI_TRANSPOSED_OUT = FFTW_MPI_TRANSPOSED_OUT
 ! end flags copies from fftw3-mpi.f03
#else
 integer,public,parameter :: ABI_FFTW_FORWARD = -1
 integer,public,parameter :: ABI_FFTW_BACKWARD = +1
 integer,public,parameter :: ABI_FFTW_ESTIMATE = 64
 integer,public,parameter :: ABI_FFTW_MEASURE = 0
! end flags copied from fftw3.f
 integer,public,parameter :: ABI_FFTW_MPI_TRANSPOSED_IN = 536870912
 integer,public,parameter :: ABI_FFTW_MPI_TRANSPOSED_OUT = 1073741824
! end flags copies from fftw3-mpi.f03
#endif

! ==========================================================================================
! ==== Variables introduced for the FFTW3 interface in abinit. Not belonging to fftw3.f ====
! ==========================================================================================

 integer,public,parameter :: NULL_PLAN = 0
 ! MKL wrappers might return NULL_PLAN if a particular FFTW3 feature is not available

 integer,public,parameter :: KIND_FFTW_PLAN = 8
 ! It should be at least integer*@SIZEOF_INT_P@
 ! MKL wrappers requires it to be integer*8, so do _not_ use C_INTPTR_T.

#ifdef HAVE_FFTW3_THREADS
 integer,private,save :: THREADS_INITED = 0
 ! 1 if treads have been initialized. 0 otherwise.
#endif

 logical,private,save :: USE_LIB_THREADS = .FALSE.
!!***

!----------------------------------------------------------------------

!!****t* m_fftw3/fftw3_plan3_t
!! NAME
!! fftw3_plan3_t
!!
!! FUNCTION
!!  Structure storing the pointer to the FFTW plan as well as the options used to generate it.
!!
!! SOURCE

 type,private :: fftw3_plan3_t
   integer :: isign=0                           ! Sign of the exponential in the FFT
   integer :: ndat=-1                           ! Number of FFTs associated to the plan
   integer :: flags=-HUGE(0)                    ! FFTW3 flags used to construct the plan.
   integer(KIND_FFTW_PLAN) :: plan=NULL_PLAN    ! FFTW3 plan.
   integer :: nthreads=1                        ! The number of threads associated to the plan.
   integer :: idist=-1
   integer :: odist=-1
   integer :: istride=-1
   integer :: ostride=-1
   integer :: n(3)=-1                           ! The number of FFT divisions.
   integer :: inembed(3)=-1
   integer :: onembed(3)=-1
   !integer(C_INT) :: alignment(2)              ! The alignment of the arrays used to construct the plan.
 end type fftw3_plan3_t
!!***

 interface fftw3_fftrisc
   module procedure fftw3_fftrisc_sp
   module procedure fftw3_fftrisc_dp
 end interface fftw3_fftrisc

 interface fftw3_fftug
   module procedure fftw3_fftug_dp
   module procedure fftw3_fftug_spc
   module procedure fftw3_fftug_dpc
 end interface fftw3_fftug

 interface fftw3_fftur
   module procedure fftw3_fftur_dp
   module procedure fftw3_fftur_spc
   module procedure fftw3_fftur_dpc
 end interface fftw3_fftur

 interface fftw3_c2c_op
   module procedure fftw3_c2c_op_spc
   module procedure fftw3_c2c_op_dpc
 end interface fftw3_c2c_op

 interface fftw3_c2c_ip
   module procedure fftw3_c2c_ip_spc
   module procedure fftw3_c2c_ip_dpc
 end interface fftw3_c2c_ip

 interface fftw3_fftpad
   module procedure fftw3_fftpad_dp
   module procedure fftw3_fftpad_spc
   module procedure fftw3_fftpad_dpc
 end interface fftw3_fftpad

#ifdef HAVE_FFTW3
  ! Overloaded planner.
 interface fftw3_plan_many_dft
   module procedure dplan_many_dft_1D
   !module procedure dplan_many_dft_2D
   module procedure cplan_many_dft
   module procedure zplan_many_dft
 end interface fftw3_plan_many_dft

 interface fftw3_execute_dft
   module procedure fftw3_execute_dft_dp
   module procedure fftw3_execute_dft_spc
   module procedure fftw3_execute_dft_dpc
 end interface fftw3_execute_dft

 interface fftw3_alloc_real
   module procedure fftw3_alloc_real1d_dp
   module procedure fftw3_alloc_real2d_dp
   !module procedure fftw3_alloc_real3d_dp
 end interface fftw3_alloc_real

 interface fftw3_alloc_complex
   module procedure fftw3_alloc_complex1d_spc
   module procedure fftw3_alloc_complex1d_dpc
 end interface fftw3_alloc_complex

!!   FDahm :In case of ffw3-mpi flavor, one must include fftw3-mpi.F03 so
!! the next few lines cause compiler errors because of name redefinition
!! In deed, i think this should be avoid and replace by 'include fftw3.f03'
!! MG: These bindings are needed when we use the FFTW3 wrappers provided by
!!     the MKL library. I agree that we should use the include file fftw3.F03
!!     but this implies that
!!     1) we cannot use the procedures defined in this module to call the MKL wrappers
!!     2) we drop support for FFTW3 versions < 3.2 since these version do not provide
!!        the F2003 interface
!!     I don't have any problem in dropping support for old versions of FFTW3
!!     We just have to handle the migration on the different slaves of the test farm
!!     builders with MKL should use m_dfti.F90
!!     builders with FFTW3 should provide a recent version of the library
!!
#ifndef HAVE_FFTW3_MPI
 ! Fortran binding for fftw_malloc
 interface fftw_malloc
   type(C_PTR) function fftw_malloc(alloc_size) bind(C, name='fftw_malloc')
     import
     integer(C_SIZE_T), value :: alloc_size
   end function fftw_malloc
 end interface fftw_malloc
 ! Fortran binding for fftw_free
 interface fftw_free
    subroutine fftw_free(cptr) bind(C, name='fftw_free')
     import
     type(C_PTR), value :: cptr
   end subroutine fftw_free
end interface
#endif
#endif

CONTAINS  !===========================================================

!!****f* m_fftw3/fftw3_seqfourdp
!! NAME
!!  fftw3_seqfourdp
!!
!! FUNCTION
!! Driver routine for 3D FFT of lengths nx, ny, nz. Mainly used for densities or potentials.
!! FFT Transform is out-of-place
!!
!! INPUTS
!! cplex=1 if fofr is real, 2 if fofr is complex
!! nx,ny,nz=Number of point along the three directions.
!! ldx,ldy,ldz=Leading dimension of the array.
!! ndat = Number of FFTS
!! isign= +1 : fofg(G) => fofr(R);
!!        -1 : fofr(R) => fofg(G)
!! fofg(2,ldx*ldy*ldz*ndat)=The array to be transformed.
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator. Defaults to FFTW_ESTIMATE.
!!
!! OUTPUT
!! fofr(cplex,ldx*ldy*ldz*ndat)=The FFT of fofg
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_seqfourdp(cplex,nx,ny,nz,ldx,ldy,ldz,ndat,isign,fofg,fofr,fftw_flags)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nx,ny,nz,ldx,ldy,ldz,ndat,isign
 integer,optional,intent(in) :: fftw_flags
!arrays
 real(dp),intent(inout) :: fofg(2*ldx*ldy*ldz*ndat)
 real(dp),intent(inout) :: fofr(cplex*ldx*ldy*ldz*ndat)

!Local variables-------------------------------
!scalars
 integer :: my_flags,ii,jj
 complex(spc), allocatable :: work_sp(:)

! *************************************************************************

 my_flags = ABI_FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags= fftw_flags

 select case (cplex)
 case (2)
   ! Complex to Complex.
   if (fftcore_mixprec == 1) then
     ! Mixed precision: copyin + in-place + copyout
     ABI_MALLOC(work_sp, (ldx*ldy*ldz*ndat))
     if (isign == ABI_FFTW_BACKWARD) then ! +1
       work_sp(:) = cmplx(fofg(1::2), fofg(2::2), kind=spc)
     else if (isign == ABI_FFTW_FORWARD) then ! -1
       work_sp(:) = cmplx(fofr(1::2), fofr(2::2), kind=spc)
     else
       MSG_BUG("Wrong isign")
     end if

     call fftw3_c2c_ip_spc(nx,ny,nz,ldx,ldy,ldz,ndat,isign,work_sp,fftw_flags=my_flags)

     if (isign == ABI_FFTW_BACKWARD) then ! +1
       jj = 1
       do ii=1,ldx*ldy*ldz*ndat
         fofr(jj) = real(work_sp(ii), kind=dp)
         fofr(jj+1) = aimag(work_sp(ii))
         jj = jj + 2
       end do
     else if (isign == ABI_FFTW_FORWARD) then  ! -1
       jj = 1
       do ii=1,ldx*ldy*ldz*ndat
         fofg(jj) = real(work_sp(ii), kind=dp)
         fofg(jj+1) = aimag(work_sp(ii))
         jj = jj + 2
       end do
     end if
     ABI_FREE(work_sp)

   else
     ! double precision version.
     select case (isign)
     case (ABI_FFTW_BACKWARD) ! +1
       call fftw3_many_dft_op(nx,ny,nz,ldx,ldy,ldz,ndat,isign,fofg,fofr,fftw_flags=my_flags)
     case (ABI_FFTW_FORWARD)  ! -1
       call fftw3_many_dft_op(nx,ny,nz,ldx,ldy,ldz,ndat,isign,fofr,fofg,fftw_flags=my_flags)
     case default
       MSG_BUG("Wrong isign")
     end select
   end if

 case (1)
   ! Real case.
   select case (isign)
   case (ABI_FFTW_FORWARD) ! -1; R --> G
     call fftw3_r2c_op(nx,ny,nz,ldx,ldy,ldz,ndat,fofr,fofg,fftw_flags=my_flags)
   case (ABI_FFTW_BACKWARD) ! +1; G --> R
     call fftw3_c2r_op(nx,ny,nz,ldx,ldy,ldz,ndat,fofg,fofr,fftw_flags=my_flags)
   case default
     MSG_BUG("Wrong isign")
   end select

 case default
   MSG_BUG(" Wrong value for cplex")
 end select

end subroutine fftw3_seqfourdp
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_seqfourwf
!! NAME
!! fftw3_seqfourwf
!!
!! FUNCTION
!! Carry out composite Fourier transforms between real and reciprocal (G) space.
!! Wavefunctions, contained in a sphere in reciprocal space,
!! can be FFT to real space. They can also be FFT from real space
!! to a sphere. Also, the density maybe accumulated, and a local potential can be applied.
!!
!! The different options are :
!! - option=0 --> reciprocal to real space and output the result.
!! - option=1 --> reciprocal to real space and accumulate the density.
!! - option=2 --> reciprocal to real space, apply the local potential to the wavefunction
!!                in real space and produce the result in reciprocal space.
!! - option=3 --> real space to reciprocal space.
!!                NOTE that in this case, fftalg=1x1 MUST be used. This may be changed in the future.
!!
!! INPUTS
!! cplex= if 1 , denpot is real, if 2 , denpot is complex
!!    (cplex=2 only allowed for option=2, and istwf_k=1)
!!    not relevant if option=0 or option=3, so cplex=0 can be used to minimize memory
!! fofgin(2,npwin)=holds input wavefunction in G vector basis sphere.
!!                 (intent(in) but the routine sphere can modify it for another iflag)
!! gboundin(2*mgfft+8,2)=sphere boundary info for reciprocal to real space
!! gboundout(2*mgfft+8,2)=sphere boundary info for real to reciprocal space
!! istwf_k=option parameter that describes the storage of wfs
!! kg_kin(3,npwin)=reduced planewave coordinates, input
!! kg_kout(3,npwout)=reduced planewave coordinates, output
!! mgfft=maximum size of 1D FFTs
!! ndat=number of FFT to do in //
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! npwin=number of elements in fofgin array (for option 0, 1 and 2)
!! npwout=number of elements in fofgout array (for option 2 and 3)
!! ldx,ldy,ldz=ngfft(4),ngfft(5),ngfft(6), dimensions of fofr.
!! option= if 0: do direct FFT
!!         if 1: do direct FFT, then sum the density
!!         if 2: do direct FFT, multiply by the potential, then do reverse FFT
!!         if 3: do reverse FFT only
!! weight_r=weight to be used for the accumulation of the density in real space
!!         (needed only when option=1)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! for option==0, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
!!                fofr(2,ldx*ldy*ldz) contains the output Fourier Transform of fofgin;
!!                no use of denpot, fofgout and npwout.
!! for option==1, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
!!                denpot(cplex*ldx,ldy,ldz) contains the input density at input,
!!                and the updated density at output (accumulated);
!!                no use of fofgout and npwout.
!! for option==2, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
!!                denpot(cplex*ldx,ldy,ldz) contains the input local potential;
!!                fofgout(2,npwout*ndat) contains the output function;
!! for option==3, fofr(2,ldx*ldy*ldz*ndat) contains the input real space wavefunction;
!!                fofgout(2,npwout*ndat) contains its output Fourier transform;
!!                no use of fofgin and npwin.
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_seqfourwf(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&
&  kg_kin,kg_kout,mgfft,ndat,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,istwf_k,ldx,ldy,ldz,ndat,npwin,npwout,option,mgfft
 real(dp),intent(in) :: weight_i,weight_r
!arrays
 integer,intent(in) :: gboundin(2*mgfft+8,2),gboundout(2*mgfft+8,2)
 integer,intent(in) :: kg_kin(3,npwin),kg_kout(3,npwout),ngfft(18)
 real(dp),intent(inout) :: denpot(cplex*ldx,ldy,ldz),fofgin(2,npwin*ndat)
 real(dp),intent(inout) :: fofr(2,ldx*ldy*ldz*ndat)
 real(dp),intent(out) :: fofgout(2,npwout*ndat)

!Local variables-------------------------------
!scalars
 integer,parameter :: me_g0=1,ndat1=1
 integer :: nx,ny,nz,fftalg,fftalga,fftalgc,fftcache,dat,ptg,ptr,ptgin,ptgout,nthreads
 character(len=500) :: msg
 logical :: use_fftrisc
!arrays
 !real(dp),allocatable :: saveden(:,:,:)
#if 0
 logical :: use_fftbox
 integer,parameter :: shiftg(3)=(/0,0,0/)
 integer :: symm(3,3)
#endif

! *************************************************************************

 nx=ngfft(1); ny=ngfft(2); nz=ngfft(3)
 fftalg=ngfft(7); fftalga=fftalg/100; fftalgc=mod(fftalg,10)
 fftcache=ngfft(8)

 if (ALL(option /= (/0,1,2,3/))) then
   write(msg,'(a,i0,a)')' The option number',option,' is not allowed. Only option=0, 1, 2 or 3 are allowed presently.'
   MSG_ERROR(msg)
 end if

 if (option==1 .and. cplex/=1) then
   write(msg,'(a,i0)')' With the option number 1, cplex must be 1 but it is cplex=',cplex
   MSG_ERROR(msg)
 end if

 if (option==2 .and. (cplex/=1 .and. cplex/=2)) then
   write(msg,'(a,i0)')' With the option number 2, cplex must be 1 or 2, but it is cplex=',cplex
   MSG_ERROR(msg)
 end if

 use_fftrisc = (fftalgc==2)
 if (istwf_k==2.and.option==3) use_fftrisc = .FALSE.
 if (istwf_k>2.and.ANY(option==(/0,3/))) use_fftrisc = .FALSE.

 nthreads = xomp_get_num_threads(open_parallel=.TRUE.)

 if (use_fftrisc) then
   !call wrtout(std_out, calling fftw3_fftrisc","COLL")

   if (ndat == 1) then
     if (fftcore_mixprec == 0) then
       call fftw3_fftrisc_dp(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
         mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
     else
       call fftw3_fftrisc_mixprec(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
         mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
     end if

   else
     ! All this boilerplate code is needed because the caller might pass zero-sized arrays
     ! for the arguments that are not referenced and we don't want to have problems at run-time.
     ! Moreover option 1 requires a special treatment when threads are started at this level.

     SELECT CASE (option)
     CASE (0)
       !
       ! fofgin -> fofr, no use of denpot, fofgout and npwout.
       if (.not.fftw3_spawn_threads_here(ndat,nthreads)) then
         do dat=1,ndat
           ptg = 1 + (dat-1)*npwin
           ptr = 1 + (dat-1)*ldx*ldy*ldz
           call fftw3_fftrisc_dp(cplex,denpot,fofgin(1,ptg),fofgout,fofr(1,ptr),gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
&            mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
         end do
       else
!$OMP PARALLEL DO PRIVATE(ptg,ptr)
         do dat=1,ndat
           ptg = 1 + (dat-1)*npwin
           ptr = 1 + (dat-1)*ldx*ldy*ldz
           call fftw3_fftrisc_dp(cplex,denpot,fofgin(1,ptg),fofgout,fofr(1,ptr),gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
&            mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
         end do
       end if

     CASE (1)
       !fofgin -> local ur and accumulate density in denpot
       ! TODO this is delicate part to do in parallel, as one should OMP reduce denpot.
       ! but this causes problems with the stack.

       do dat=1,ndat
         ptg = 1 + (dat-1)*npwin
         ptr = 1 + (dat-1)*ldx*ldy*ldz
         call fftw3_fftrisc_dp(cplex,denpot,fofgin(1,ptg),fofgout,fofr,gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
&          mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
       end do

       ! This version doesn't seem efficient
       !!!  !$OMP PARALLEL PRIVATE(ptg,ptr,saveden)
       !!!         ABI_MALLOC(saveden, (ldx,ldy,ldz))
       !!!         saveden = zero
       !!!  !$OMP DO
       !!!         do dat=1,ndat
       !!!           ptg = 1 + (dat-1)*npwin
       !!!           ptr = 1 + (dat-1)*ldx*ldy*ldz
       !!!           call fftw3_fftrisc_dp(cplex,saveden,fofgin(1,ptg),fofgout,fofr,gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
       !!!  &          mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r)
       !!!         end do
       !!!  !$OMP END DO NOWAIT
       !!!  !$OMP CRITICAL (OMPC_addrho)
       !!!         denpot = denpot + saveden
       !!!  !$OMP END CRITICAL (OMPC_addrho)
       !!!         ABI_FREE(saveden)
       !!!  !$OMP END PARALLEL

     CASE (2)
       ! <G|vloc(r)|fofgin(r)> in fofgout
       if (.not.fftw3_spawn_threads_here(ndat,nthreads)) then
         do dat=1,ndat
           ptgin  = 1 + (dat-1)*npwin
           ptgout = 1 + (dat-1)*npwout
           if (fftcore_mixprec == 0) then
             call fftw3_fftrisc_dp(cplex,denpot,fofgin(1,ptgin),fofgout(1,ptgout),fofr,gboundin,gboundout,&
                 istwf_k,kg_kin,kg_kout,mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
           else
             call fftw3_fftrisc_mixprec(cplex,denpot,fofgin(1,ptgin),fofgout(1,ptgout),fofr,gboundin,gboundout,&
                 istwf_k,kg_kin,kg_kout,mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
           end if
         end do
       else
!$OMP PARALLEL DO PRIVATE(ptgin,ptgout)
         do dat=1,ndat
           ptgin  = 1 + (dat-1)*npwin
           ptgout = 1 + (dat-1)*npwout
           call fftw3_fftrisc_dp(cplex,denpot,fofgin(1,ptgin),fofgout(1,ptgout),fofr,gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
&            mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
         end do
       end if

     CASE (3)
       ! fofr -> fofgout
       if (.not.fftw3_spawn_threads_here(ndat,nthreads)) then
         do dat=1,ndat
           ptr    = 1 + (dat-1)*ldx*ldy*ldz
           ptgout = 1 + (dat-1)*npwout
           call fftw3_fftrisc_dp(cplex,denpot,fofgin,fofgout(1,ptgout),fofr(1,ptr),gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
&            mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
         end do
       else
!$OMP PARALLEL DO PRIVATE(ptr,ptgout)
         do dat=1,ndat
           ptr    = 1 + (dat-1)*ldx*ldy*ldz
           ptgout = 1 + (dat-1)*npwout
           call fftw3_fftrisc_dp(cplex,denpot,fofgin,fofgout(1,ptgout),fofr(1,ptr),gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
&            mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
         end do
       end if

     CASE DEFAULT
       write(msg,'(a,i0,a)')'Option',option,' is not allowed. Only option=0, 1, 2 or 3 are allowed presently.'
       MSG_ERROR(msg)
     END SELECT

   end if

 else

#if 1
   SELECT CASE (option)
   CASE (0)
     !
     ! FFT u(g) --> u(r)
     if (.not.fftw3_spawn_threads_here(ndat,nthreads)) then
       call fftw3_fftug_dp(fftalg,fftcache,npwin,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_kin,gboundin,fofgin,fofr)
     else
!$OMP PARALLEL DO PRIVATE(ptg, ptr)
       do dat=1,ndat
         ptg = 1 + (dat-1)*npwin
         ptr = 1 + (dat-1)*ldx*ldy*ldz
         call fftw3_fftug_dp(fftalg,fftcache,npwin,nx,ny,nz,ldx,ldy,ldz,ndat1,&
&          istwf_k,mgfft,kg_kin,gboundin,fofgin(1,ptg),fofr(1,ptr))
       end do
     end if

   CASE (1)
     ! TODO this is delicate part to do in parallel, as one should OMP reduce denpot.
     call fftw3_fftug_dp(fftalg,fftcache,npwin,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_kin,gboundin,fofgin,fofr)
     call cg_addtorho(nx,ny,nz,ldx,ldy,ldz,ndat,weight_r,weight_i,fofr,denpot)

   CASE (2)

     if (.not.fftw3_spawn_threads_here(ndat,nthreads)) then
       call fftw3_fftug_dp(fftalg,fftcache,npwin,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_kin,gboundin,fofgin,fofr)
       call cg_vlocpsi(nx,ny,nz,ldx,ldy,ldz,ndat,cplex,denpot,fofr)

       !  The data for option==2 is now in fofr.
       call fftw3_fftpad_dp(fofr,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,-1,gboundout)

       call cg_box2gsph(nx,ny,nz,ldx,ldy,ldz,ndat,npwout,kg_kout,fofr,fofgout)
     else

!$OMP PARALLEL DO PRIVATE(ptg, ptr)
       do dat=1,ndat
         ptg = 1 + (dat-1)*npwin
         ptr = 1 + (dat-1)*ldx*ldy*ldz
         call fftw3_fftug_dp(fftalg,fftcache,npwin,nx,ny,nz,ldx,ldy,ldz,ndat1,&
&          istwf_k,mgfft,kg_kin,gboundin,fofgin(1,ptg),fofr(1,ptr))

         call cg_vlocpsi(nx,ny,nz,ldx,ldy,ldz,ndat1,cplex,denpot,fofr(1,ptr))

         !  The data for option==2 is now in fofr.
         call fftw3_fftpad_dp(fofr(1,ptr),nx,ny,nz,ldx,ldy,ldz,ndat1,mgfft,-1,gboundout)

         ptg = 1 + (dat-1)*npwout
         call cg_box2gsph(nx,ny,nz,ldx,ldy,ldz,ndat1,npwout,kg_kout,fofr(1,ptr),fofgout(1,ptg))
       end do
     end if

   CASE (3)
     !  The data for option==3 is already in fofr.
     if (.not.fftw3_spawn_threads_here(ndat,nthreads)) then
       call fftw3_fftpad_dp(fofr,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,-1,gboundout)
       call cg_box2gsph(nx,ny,nz,ldx,ldy,ldz,ndat,npwout,kg_kout,fofr,fofgout)
     else
!$OMP PARALLEL DO PRIVATE(ptg, ptr)
       do dat=1,ndat
         ptg = 1 + (dat-1)*npwout
         ptr = 1 + (dat-1)*ldx*ldy*ldz
         call fftw3_fftpad_dp(fofr(1,ptr),nx,ny,nz,ldx,ldy,ldz,ndat1,mgfft,-1,gboundout)
         call cg_box2gsph(nx,ny,nz,ldx,ldy,ldz,ndat1,npwout,kg_kout,fofr(1,ptr),fofgout(1,ptg))
       end do
     end if

   CASE DEFAULT
     write(msg,'(a,i0,a)')'Option',option,' is not allowed. Only option=0, 1, 2 or 3 are allowed presently.'
     MSG_ERROR(msg)
   END SELECT


#else
   symm=0; symm(1,1)=1; symm(2,2)=1; symm(3,3)=1
   use_fftbox = .FALSE.
#ifdef HAVE_OPENMP
   use_fftbox = (ndat>1)
#endif
   !use_fftbox = .TRUE.

   SELECT CASE (option)
   CASE (0)
     !
     ! FFT u(g) --> u(r)
     call sphere(fofgin,ndat,npwin,fofr,nx,ny,nz,ldx,ldy,ldz,kg_kin,istwf_k,1,me_g0,shiftg,symm,one)

     if (use_fftbox) then
       call fftw3_many_dft_ip(nx,ny,nz,ldx,ldy,ldz,ndat,ABI_FFTW_BACKWARD,fofr)
     else
       call fftw3_fftpad_dp(fofr,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,ABI_FFTW_BACKWARD,gboundin)
     end if

   CASE (1)
     ! TODO this is delicate part to do in parallel, as one should OMP reduce denpot.

     call sphere(fofgin,ndat,npwin,fofr,nx,ny,nz,ldx,ldy,ldz,kg_kin,istwf_k,1,me_g0,shiftg,symm,one)

     if (use_fftbox) then
       call fftw3_many_dft_ip(nx,ny,nz,ldx,ldy,ldz,ndat,ABI_FFTW_BACKWARD,fofr)
     else
       call fftw3_fftpad_dp(fofr,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,ABI_FFTW_BACKWARD,gboundin)
     end if

     call cg_addtorho(nx,ny,nz,ldx,ldy,ldz,ndat,weight_r,weight_i,fofr,denpot)

   CASE (2)

     call sphere(fofgin,ndat,npwin,fofr,nx,ny,nz,ldx,ldy,ldz,kg_kin,istwf_k,1,me_g0,shiftg,symm,one)

     if (use_fftbox) then
       call fftw3_many_dft_ip(nx,ny,nz,ldx,ldy,ldz,ndat,ABI_FFTW_BACKWARD,fofr)
     else
       call fftw3_fftpad_dp(fofr,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,ABI_FFTW_BACKWARD,gboundin)
     end if

     call cg_vlocpsi(nx,ny,nz,ldx,ldy,ldz,ndat,cplex,denpot,fofr)

     ! The data for option==2 is now in fofr.
     if (use_fftbox) then
       call fftw3_many_dft_ip(nx,ny,nz,ldx,ldy,ldz,ndat,ABI_FFTW_FORWARD,fofr)
     else
       call fftw3_fftpad_dp(fofr,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,ABI_FFTW_FORWARD,gboundout)
     end if

     call cg_box2gsph(nx,ny,nz,ldx,ldy,ldz,ndat,npwout,kg_kout,fofr,fofgout)

   CASE (3)
     !  The data for option==3 is already in fofr.
     call fftw3_fftpad_dp(fofr,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,ABI_FFTW_FORWARD,gboundout)

     call cg_box2gsph(nx,ny,nz,ldx,ldy,ldz,ndat,npwout,kg_kout,fofr,fofgout)

   CASE DEFAULT
     write(msg,'(a,i0,a)')'Option',option,' is not allowed. Only option=0, 1, 2 or 3 are allowed presently.'
     MSG_ERROR(msg)
   END SELECT
#endif

 end if

end subroutine fftw3_seqfourwf
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_fftrisc_sp
!! NAME
!! fftw3_fftrisc_sp
!!
!! FUNCTION
!! Carry out Fourier transforms between real and reciprocal (G) space,
!! for wavefunctions, contained in a sphere in reciprocal space,
!! in both directions. Also accomplish some post-processing.
!! See fftw3_fftrisc_dp for API doc.
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_fftrisc_sp(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
& mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,istwf_k,mgfft,ldx,ldy,ldz,npwin,npwout,option
 real(dp),intent(in) :: weight_i,weight_r
!arrays
 integer,intent(in) :: gboundin(2*mgfft+8,2),gboundout(2*mgfft+8,2)
 integer,intent(in) :: kg_kin(3,npwin),kg_kout(3,npwout),ngfft(18)
 real(sp),intent(in) :: fofgin(2,npwin)
 real(dp),intent(inout) :: denpot(cplex*ldx,ldy,ldz)
 real(sp),intent(inout) :: fofr(2,ldx*ldy*ldz)
 real(sp),intent(inout) :: fofgout(2,npwout)

! *************************************************************************

#ifdef HAVE_FFTW3

#undef  FFT_PRECISION
#undef  MYKIND
#undef  MYCZERO
#undef  MYCMPLX
#undef  MYCONJG

#define FFT_PRECISION FFT_SINGLE
#define MYKIND SPC
#define MYCZERO (0._sp,0._sp)
#define MYCMPLX  CMPLX
#define MYCONJG  CONJG

#include "fftw3_fftrisc.finc"

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/cplex,gboundin(1,1),gboundout(1,1),istwf_k,kg_kin(1,1),kg_kout(1,1)/))
 ABI_UNUSED((/mgfft,ngfft(1),npwin,npwout,ldx,ldy,ldz,option/))
 ABI_UNUSED((/denpot(1,1,1),weight_r,weight_i/))
 ABI_UNUSED((/fofgin(1,1),fofgout(1,1),fofr(1,1)/))
#endif

end subroutine fftw3_fftrisc_sp
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_fftrisc_dp
!! NAME
!! fftw3_fftrisc_dp
!!
!! FUNCTION
!! Carry out Fourier transforms between real and reciprocal (G) space,
!! for wavefunctions, contained in a sphere in reciprocal space,
!! in both directions. Also accomplish some post-processing.
!!
!! NOTES
!! Specifically uses rather sophisticated algorithms, based on S Goedecker
!! routines, specialized for superscalar RISC architecture.
!! Zero padding : saves 7/12 execution time
!! Bi-dimensional data locality in most of the routine : cache reuse
!! For k-point (0 0 0) : takes advantage of symmetry of data.
!! Note however that no blocking is used, in both 1D z-transform
!! or subsequent 2D transform. This should be improved.
!!
!! INPUTS
!!  cplex= if 1 , denpot is real, if 2 , denpot is complex
!!     (cplex=2 only allowed for option=2 when istwf_k=1)
!!     one can also use cplex=0 if option=0 or option=3
!!  fofgin(2,npwin)=holds input wavefunction in G vector basis sphere.
!!  gboundin(2*mgfft+8,2)=sphere boundary info for reciprocal to real space
!!  gboundout(2*mgfft+8,2)=sphere boundary info for real to reciprocal space
!!  istwf_k=option parameter that describes the storage of wfs
!!  kg_kin(3,npwin)=reduced planewave coordinates, input
!!  kg_kout(3,npwout)=reduced planewave coordinates, output
!!  mgfft=maximum size of 1D FFTs
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  npwin=number of elements in fofgin array (for option 0, 1 and 2)
!!  npwout=number of elements in fofgout array (for option 2 and 3)
!!  ldx,ldy,ldz=ngfft(4),ngfft(5),ngfft(6), dimensions of fofr.
!!  option= if 0: do direct FFT
!!          if 1: do direct FFT, then sum the density
!!          if 2: do direct FFT, multiply by the potential, then do reverse FFT
!!          if 3: do reverse FFT only
!!  weight=weight to be used for the accumulation of the density in real space
!!          (needed only when option=1)
!!
!! OUTPUT
!!  (see side effects)
!!
!! OPTIONS
!!  The different options are:
!!  - reciprocal to real space and output the result (when option=0),
!!  - reciprocal to real space and accumulate the density (when option=1) or
!!  - reciprocal to real space, apply the local potential to the wavefunction
!!    in real space and produce the result in reciprocal space (when option=2)
!!  - real space to reciprocal space (when option=3).
!!  option=0 IS NOT ALLOWED when istwf_k>2
!!  option=3 IS NOT ALLOWED when istwf_k>=2
!!
!! SIDE EFFECTS
!!  for option==0, fofgin(2,npwin)=holds input wavefunction in G sphere;
!!                 fofr(2,ldx,ldy,ldz) contains the Fourier Transform of fofgin;
!!                 no use of denpot, fofgout and npwout.
!!  for option==1, fofgin(2,npwin)=holds input wavefunction in G sphere;
!!                 denpot(cplex*ldx,ldy,ldz) contains the input density at input,
!!                 and the updated density at output;
!!                 no use of fofgout and npwout.
!!  for option==2, fofgin(2,npwin)=holds input wavefunction in G sphere;
!!                 denpot(cplex*ldx,ldy,ldz) contains the input local potential;
!!                 fofgout(2,npwout) contains the output function;
!!  for option==3, fofr(2,ldx,ldy,ldz) contains the real space wavefunction;
!!                 fofgout(2,npwout) contains its Fourier transform;
!!                 no use of fofgin and npwin.
!!
!! PARENTS
!!      m_fftw3
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_fftrisc_dp(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
& mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,istwf_k,mgfft,ldx,ldy,ldz,npwin,npwout,option
 real(dp),intent(in) :: weight_r,weight_i
!arrays
 integer,intent(in) :: gboundin(2*mgfft+8,2),gboundout(2*mgfft+8,2)
 integer,intent(in) :: kg_kin(3,npwin),kg_kout(3,npwout),ngfft(18)
 real(dp),intent(in) :: fofgin(2,npwin)
 real(dp),intent(inout) :: denpot(cplex*ldx,ldy,ldz),fofr(2,ldx*ldy*ldz)
 real(dp),intent(inout) :: fofgout(2,npwout)

! *************************************************************************

#ifdef HAVE_FFTW3

#undef  FFT_PRECISION
#undef  MYKIND
#undef  MYCZERO
#undef  MYCMPLX
#undef  MYCONJG

#define FFT_PRECISION FFT_DOUBLE
#define MYKIND DPC
#define MYCZERO (0._dp,0._dp)
#define MYCMPLX  DCMPLX
#define MYCONJG  DCONJG

#include "fftw3_fftrisc.finc"

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/cplex,gboundin(1,1),gboundout(1,1),istwf_k,kg_kin(1,1),kg_kout(1,1)/))
 ABI_UNUSED((/mgfft,ngfft(1),npwin,npwout,ldx,ldy,ldz,option/))
 ABI_UNUSED((/denpot(1,1,1),fofgin(1,1),fofgout(1,1),fofr(1,1),weight_r,weight_i/))
#endif

end subroutine fftw3_fftrisc_dp
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_fftrisc_mixprec
!! NAME
!! fftw3_fftrisc_mixprec
!!
!! FUNCTION
!!  Mixed precision version of fftrisc: input/output in dp, computation done in sp.
!!  See fftw3_fftrisc_dp for API docs.
!!
!! PARENTS
!!      m_fftw3
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_fftrisc_mixprec(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
& mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,istwf_k,mgfft,ldx,ldy,ldz,npwin,npwout,option
 real(dp),intent(in) :: weight_r,weight_i
!arrays
 integer,intent(in) :: gboundin(2*mgfft+8,2),gboundout(2*mgfft+8,2)
 integer,intent(in) :: kg_kin(3,npwin),kg_kout(3,npwout),ngfft(18)
 real(dp),intent(in) :: fofgin(2,npwin)
 real(dp),intent(inout) :: denpot(cplex*ldx,ldy,ldz),fofr(2,ldx*ldy*ldz)
 real(dp),intent(inout) :: fofgout(2,npwout)

! *************************************************************************

#ifdef HAVE_FFTW3

#undef  FFT_PRECISION
#undef  MYKIND
#undef  MYCZERO
#undef  MYCMPLX
#undef  MYCONJG

#define FFT_PRECISION FFT_MIXPREC
#define MYKIND SPC
#define MYCZERO (0._sp,0._sp)
#define MYCMPLX  CMPLX
#define MYCONJG  CONJG

#include "fftw3_fftrisc.finc"

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/cplex,gboundin(1,1),gboundout(1,1),istwf_k,kg_kin(1,1),kg_kout(1,1)/))
 ABI_UNUSED((/mgfft,ngfft(1),npwin,npwout,ldx,ldy,ldz,option/))
 ABI_UNUSED((/denpot(1,1,1),fofgin(1,1),fofgout(1,1),fofr(1,1),weight_r,weight_i/))
#endif

end subroutine fftw3_fftrisc_mixprec
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_fftug_dp
!! NAME
!! fftw3_fftug_dp
!!
!! FUNCTION
!! Compute ndat zero-padded FFTs from G to R space.
!! Mainly used for the transform of wavefunctions.
!! TARGET: dp arrays with real and imaginary part
!!
!! INPUTS
!! fftalg=FFT algorith (see input variable)
!! fftcache=size of the cache (kB)
!! npw_k=number of plane waves for this k-point.
!! nx,ny,nz=Number of point along the three directions.
!! ldx,ldy,ldz=Leading dimensions of the array.
!! ndat=Number of transforms
!! istwf_k=Option describing the storage of the wavefunction.
!! mgfft=Max number of FFT divisions (used to dimension gbound)
!! kg_k(3,npw_k)=G-vectors in reduced coordinates
!! gbound(2*mgfft+8,2)=Table for zero-padded FFT. See sphereboundary.
!!  ug(npw_k*ndat)=wavefunctions in reciprocal space.
!!
!! OUTPUT
!!  ur(ldx*ldy*ldz*ndat)=wavefunctions in real space.
!!
!! PARENTS
!!      m_fftw3
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_fftug_dp(fftalg,fftcache,npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k,gbound,ug,ur)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalg,fftcache
 integer,intent(in) :: npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2),kg_k(3,npw_k)
 real(dp),target,intent(in) :: ug(2*npw_k*ndat)
 real(dp),target,intent(inout) :: ur(2*ldx*ldy*ldz*ndat)

#ifdef HAVE_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: dist=2
 real(dp) :: fofgout(2,0)
 real(dp),ABI_CONTIGUOUS pointer :: real_ug(:,:),real_ur(:,:)

! *************************************************************************

#undef TK_PREF
#define TK_PREF(name) CONCAT(cg_,name)

#include "fftug.finc"

#else
 ! Silence compiler warning
 MSG_ERROR("FFT_FFTW3 support not activated")
 ABI_UNUSED((/fftalg,fftcache,npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k(1,1),gbound(1,1)/))
 ABI_UNUSED((/ug(1),ur(1)/))
#endif

end subroutine fftw3_fftug_dp
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_fftug_spc
!! NAME
!! fftw3_fftug_spc
!!
!! FUNCTION
!! Compute ndat zero-padded FFTs from G-->R.
!! Mainly used for the transform of wavefunctions.
!! TARGET: spc arrays
!!
!! INPUTS
!! fftalg=FFT algorith (see input variable)
!! fftcache=size of the cache (kB)
!! npw_k=number of plane waves for this k-point.
!! nx,ny,nz=Number of point along the three directions.
!! ldx,ldy,ldz=Leading dimensions of the array.
!! ndat=Number of transforms
!! istwf_k=Option describing the storage of the wavefunction.
!! mgfft=Max number of FFT divisions (used to dimension gbound)
!! kg_k(3,npw_k)=G-vectors in reduced coordinates
!! gbound(2*mgfft+8,2)=Table for padded-FFT. See sphereboundary.
!!  ug(npw_k*ndat)=wavefunctions in reciprocal space.
!!
!! OUTPUT
!!  ur(ldx*ldy*ldz*ndat)=wavefunctions in real space
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_fftug_spc(fftalg,fftcache,npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k,gbound,ug,ur)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalg,fftcache
 integer,intent(in) :: npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2),kg_k(3,npw_k)
 complex(spc),target,intent(in) :: ug(npw_k*ndat)
 complex(spc),target,intent(inout) :: ur(ldx*ldy*ldz*ndat)

#ifdef HAVE_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: dist=1
!arrays
 real(sp) :: fofgout(2,0)
 real(sp),ABI_CONTIGUOUS pointer :: real_ug(:,:),real_ur(:,:)

! *************************************************************************

#undef TK_PREF
#define TK_PREF(name) CONCAT(cplx_,name)

#include "fftug.finc"

#else
 ! Silence compiler warning
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/fftalg,fftcache,npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k(1,1),gbound(1,1)/))
 ABI_UNUSED((/ug(1),ur(1)/))
#endif

end subroutine fftw3_fftug_spc
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_fftug_dpc
!! NAME
!! fftw3_fftug_dpc
!!
!! FUNCTION
!! Compute ndat zero-padded FFTs.
!! Mainly used for the transform of wavefunctions.
!! TARGET: DPC arrays
!!
!! INPUTS
!! fftalg=FFT algorith (see input variable)
!! fftcache=size of the cache (kB)
!! npw_k=number of plane waves for this k-point.
!! nx,ny,nz=Number of point along the three directions.
!! ldx,ldy,ldz=Leading dimensions of the array.
!! ndat=Number of transforms
!! istwf_k=Option describing the storage of the wavefunction.
!! mgfft=Max number of FFT divisions (used to dimension gbound)
!! kg_k(3,npw_k)=G-vectors in reduced coordinates
!! gbound(2*mgfft+8,2)=Table for padded-FFT. See sphereboundary.
!!  ug(npw_k*ndat)=wavefunctions in reciprocal space
!!
!! OUTPUT
!!  ur(ldx*ldy*ldz*ndat)=wavefunctions in real space.
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_fftug_dpc(fftalg,fftcache,npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k,gbound,ug,ur)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalg,fftcache
 integer,intent(in) :: npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2),kg_k(3,npw_k)
 complex(dpc),target,intent(in) :: ug(npw_k*ndat)
 complex(dpc),target,intent(inout) :: ur(ldx*ldy*ldz*ndat)

#ifdef HAVE_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: dist=1
!arrays
 real(dp) :: fofgout(2,0)
 real(dp),ABI_CONTIGUOUS pointer :: real_ug(:,:),real_ur(:,:)

! *************************************************************************

#undef TK_PREF
#define TK_PREF(name) CONCAT(cplx_,name)

#include "fftug.finc"

#else
 ! Silence compiler warning
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/fftalg,fftcache,npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k(1,1),gbound(1,1)/))
 ABI_UNUSED((/ug(1),ur(1)/))
#endif

end subroutine fftw3_fftug_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_fftur_dp
!! NAME
!! fftw3_fftur_dp
!!
!! FUNCTION
!! Compute ndat zero-padded FFTs from R- to G-space .
!! Mainly used for the transform of wavefunctions.
!! TARGET: dp arrays
!!
!! INPUTS
!! fftalg=FFT algorith (see input variable)
!! fftcache=size of the cache (kB)
!! npw_k=number of plane waves for this k-point.
!! nx,ny,nz=Number of point along the three directions.
!! ldx,ldy,ldz=Leading dimensions of the array.
!! ndat=Number of transforms
!! istwf_k=Option describing the storage of the wavefunction.
!! mgfft=Max number of FFT divisions (used to dimension gbound)
!! kg_k(3,npw_k)=G-vectors in reduced coordinates
!! gbound(2*mgfft+8,2)=Table for padded-FFT. See sphereboundary.
!!
!! SIDE EFFECT
!! ur(ldx*ldy*ldz*ndat)= In input: wavefunctions in real space.
!!                       Destroyed in output. Do not use ur anymore!
!! OUTPUT
!! ug(npw_k*ndat)=wavefunctions in reciprocal space.
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_fftur_dp(fftalg,fftcache,npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k,gbound,ur,ug)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalg,fftcache
 integer,intent(in) :: npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2),kg_k(3,npw_k)
 real(dp),target,intent(inout) :: ur(2*ldx*ldy*ldz*ndat)
 real(dp),target,intent(inout) :: ug(2*npw_k*ndat)

#ifdef HAVE_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: dist=2
!arrays
 real(dp) :: dum_ugin(2,0)
 real(dp),ABI_CONTIGUOUS pointer :: real_ug(:,:),real_ur(:,:)

! *************************************************************************

#undef TK_PREF
#define TK_PREF(name) CONCAT(cg_,name)

#include "fftur.finc"

#else
 ! Silence compiler warning
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/fftalg,fftcache/))
 ABI_UNUSED((/npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k(1,1),gbound(1,1)/))
 ABI_UNUSED((/ug(1),ur(1)/))
#endif

end subroutine fftw3_fftur_dp
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_fftur_spc
!! NAME
!! fftw3_fftur_spc
!!
!! FUNCTION
!! Compute ndat zero-padded FFTs from R- to G-space .
!! Mainly used for the transform of wavefunctions.
!! TARGET: spc arrays
!!
!! INPUTS
!! fftalg=FFT algorith (see input variable)
!! fftcache=size of the cache (kB)
!! npw_k=number of plane waves for this k-point.
!! nx,ny,nz=Number of point along the three directions.
!! ldx,ldy,ldz=Leading dimensions of the array.
!! ndat=Number of transforms
!! istwf_k=Option describing the storage of the wavefunction.
!! mgfft=Max number of FFT divisions (used to dimension gbound)
!! kg_k(3,npw_k)=G-vectors in reduced coordinates
!! gbound(2*mgfft+8,2)=Table for padded-FFT. See sphereboundary.
!!
!! SIDE EFFECT
!! ur(ldx*ldy*ldz*ndat)= In input: wavefunctions in real space.
!!                       Destroyed in output. Do not use ur anymore!
!! OUTPUT
!! ug(npw_k*ndat)=wavefunctions in reciprocal space.
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_fftur_spc(fftalg,fftcache,npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k,gbound,ur,ug)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalg,fftcache
 integer,intent(in) :: npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2),kg_k(3,npw_k)
 complex(spc),target,intent(inout) :: ur(ldx*ldy*ldz*ndat)
 complex(spc),target,intent(inout) :: ug(npw_k*ndat)

#ifdef HAVE_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: dist=1
!arrays
 real(sp) :: dum_ugin(2,0)
 real(sp),ABI_CONTIGUOUS pointer :: real_ug(:,:),real_ur(:,:)

! *************************************************************************

#undef TK_PREF
#define TK_PREF(name) CONCAT(cplx_,name)

#include "fftur.finc"

#else
 ! Silence compiler warning
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/fftalg,fftcache/))
 ABI_UNUSED((/npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k(1,1),gbound(1,1)/))
 ABI_UNUSED((/ug(1),ur(1)/))
#endif

end subroutine fftw3_fftur_spc
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_fftur_dpc
!! NAME
!! fftw3_fftur_dpc
!!
!! FUNCTION
!! Compute ndat zero-padded FFTs from R ro G.
!! Mainly used for the transform of wavefunctions.
!! TARGET: DPC arrays
!!
!! INPUTS
!! fftalg=FFT algorith (see input variable)
!! fftcache=size of the cache (kB)
!! npw_k=number of plane waves for this k-point.
!! nx,ny,nz=Number of point along the three directions.
!! ldx,ldy,ldz=Leading dimensions of the array.
!! ndat=Number of transforms
!! istwf_k=Option describing the storage of the wavefunction.
!! mgfft=Max number of FFT divisions (used to dimension gbound)
!! kg_k(3,npw_k)=G-vectors in reduced coordinates
!! gbound(2*mgfft+8,2)=Table for padded-FFT. See sphereboundary.
!!
!! SIDE EFFECT
!! ur(ldx*ldy*ldz*ndat)= In input: wavefunctions in real space.
!!                       Destroyed in output. Do not use ur anymore!
!! OUTPUT
!! ug(npw_k*ndat)=wavefunctions in reciprocal space
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_fftur_dpc(fftalg,fftcache,npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k,gbound,ur,ug)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalg,fftcache
 integer,intent(in) :: npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2),kg_k(3,npw_k)
 complex(dpc),target,intent(inout) :: ur(ldx*ldy*ldz*ndat)
 complex(dpc),target,intent(inout) :: ug(npw_k*ndat)

#ifdef HAVE_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: dist=1
!arrays
 real(dp) :: dum_ugin(2,0)
 real(dp),ABI_CONTIGUOUS pointer :: real_ug(:,:),real_ur(:,:)

! *************************************************************************

#undef TK_PREF
#define TK_PREF(name) CONCAT(cplx_,name)

#include "fftur.finc"

#else
 ! Silence compiler warning
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/fftalg,fftcache/))
 ABI_UNUSED((/npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k(1,1),gbound(1,1)/))
 ABI_UNUSED((/ug(1),ur(1)/))
#endif

end subroutine fftw3_fftur_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_c2c_ip_spc
!! NAME
!!  fftw3_c2c_ip_spc
!!
!! FUNCTION
!! Driver routine for in-place 3D complex-complex FFT.
!! TARGET: Simple precision complex arrays.
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimensions of the array.
!! ndat=Number of FFTs to be done.
!! isign= +1 : ff(G) => ff(R); -1 : ff(R) => ff(G)
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator.
!!   Defaults to ABI_FFTW_ESTIMATE.
!!
!! SIDE EFFECTS
!!  ff(ldx*ldy*ldz*ndat)=
!!    In input: the complex array to be transformed.
!!    In output: the Fourier transformed in the space specified by isign.
!!
!! PARENTS
!!      m_fftw3
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_c2c_ip_spc(nx,ny,nz,ldx,ldy,ldz,ndat,isign,ff,fftw_flags)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,isign
 integer,optional,intent(in) :: fftw_flags
!arrays
 complex(spc),intent(inout) :: ff(ldx*ldy*ldz*ndat)

#ifdef HAVE_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: rank3=3,nt_all=-1
 integer :: my_flags,dist,stride
 integer(KIND_FFTW_PLAN) :: my_plan
!arrays
 integer :: embed(rank3),n(rank3)

! *************************************************************************

 my_flags=ABI_FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags=fftw_flags

 stride = 1
 dist   = ldx*ldy*ldz
 embed  = (/ldx,ldy,ldz/)
 n      = (/nx ,ny ,nz /)

 my_plan = fftw3_plan_many_dft(rank3, n, ndat, ff, embed, stride, dist, ff, embed, stride, dist, isign, my_flags, nt_all)

 ! Now perform the 3D FFT via FFTW.
 call sfftw_execute_dft(my_plan, ff, ff)

 call fftw3_destroy_plan(my_plan)

 if (isign==ABI_FFTW_FORWARD) then ! -1, FFTW returns not normalized FTs
   call xscal(ldx*ldy*ldz*ndat, REAL(one/(nx*ny*nz),KIND=sp), ff, 1)
 end if

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,isign/))
 ABI_UNUSED(ff)
 if (PRESENT(fftw_flags)) then
   ABI_UNUSED(fftw_flags)
 end if
#endif

end subroutine fftw3_c2c_ip_spc
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_fftpad_spc
!! NAME
!!  fftw3_fftpad_spc
!!
!! FUNCTION
!!  This routine transforms wavefunctions using 3D zero-padded FFTs with FFTW3.
!!  The 3D ffts are computed only on lines and planes which have non zero elements.
!!  These lines and planes are defined by the two vectors do_fft_x(ldy*nz) and do_fft_y(nz)
!!  FFT transform is in-place. Target: complex arrays.
!!
!! INPUTS
!!   nx,ny,nz=Logical dimensions of the FFT mesh.
!!   ldx,ldy,ldz=Physical dimension of the f array (to avoid cache conflicts).
!!   ndat=Number of FFT transforms.
!!   mgfft=MAX(nx,ny,nz), only used to dimension gbound.
!!   isign=The sign of the transform.
!!   gbound(2*mgfft+8,2)= The boundaries of the basis sphere of G vectors at a given k-point.
!!     See sphereboundary for more info.
!!
!! SIDE EFFECTS
!!  ff(ldx*ldy*ldz*ndat)=
!!    input: The array with the data to be transformed.
!!    output: The results of the FFT.
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_fftpad_spc(ff,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign,gbound)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2)
 complex(spc),intent(inout) :: ff(ldx*ldy*ldz*ndat)

#ifdef HAVE_FFTW3
!Local variables-------------------------------
 integer,parameter :: dst=1
 real(sp) :: fact

! *************************************************************************

#include "fftw3_fftpad.finc"

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign/))
 ABI_UNUSED(gbound(1,1))
 ABI_UNUSED(ff(1))
#endif

end subroutine fftw3_fftpad_spc
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_c2c_ip_dpc
!! NAME
!!  fftw3_c2c_ip_dpc
!!
!! FUNCTION
!! Driver routine for in-place 3D complex-complex FFT.
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimensions of the array.
!! ndat=Number of FFTs to be done.
!! isign= +1 : ff(G) => ff(R); -1 : ff(R) => ff(G)
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator.
!!   Defaults to ABI_FFTW_ESTIMATE.
!!
!! SIDE EFFECTS
!!  ff(ldx*ldy*ldz*ndat)=
!!    In input: the complex array to be transformed.
!!    In output: the Fourier transformed in the space specified by isign.
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_c2c_ip_dpc(nx,ny,nz,ldx,ldy,ldz,ndat,isign,ff,fftw_flags)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,isign
 integer,optional,intent(in) :: fftw_flags
!arrays
 complex(dpc),intent(inout) :: ff(ldx*ldy*ldz*ndat)

#ifdef HAVE_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: rank3=3,nt_all=-1
 integer :: my_flags,dist,stride
 integer(KIND_FFTW_PLAN) :: my_plan
!arrays
 integer :: embed(rank3),n(rank3)

! *************************************************************************

 my_flags=ABI_FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags=fftw_flags

 stride = 1
 dist   = ldx*ldy*ldz
 embed  = (/ldx,ldy,ldz/)
 n      = (/nx ,ny ,nz /)

 my_plan = fftw3_plan_many_dft(rank3, n, ndat, ff, embed, stride, dist, ff, embed, stride, dist, isign, my_flags, nt_all)

 ! Now perform the 3D FFT via FFTW.
 call dfftw_execute_dft(my_plan, ff, ff)

 call fftw3_destroy_plan(my_plan)

 if (isign==ABI_FFTW_FORWARD) then ! -1, FFTW returns not normalized FTs
  call ZDSCAL(ldx*ldy*ldz*ndat, one/(nx*ny*nz), ff, 1)
 end if

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,isign/))
 ABI_UNUSED(ff)
 if (PRESENT(fftw_flags)) then
   ABI_UNUSED(fftw_flags)
 end if
#endif

end subroutine fftw3_c2c_ip_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_c2c_op_spc
!! NAME
!!  fftw3_c2c_op_spc
!!
!! FUNCTION
!! Driver routine for out-of-place 3D complex-complex FFT of lengths nx, ny, nz.
!! TARGET: single precision complex arrays
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimensions of the array.
!! ndat=Number of FFTs to be done.
!! isign= +1 : ff(G) => gg(R); -1 : ff(R) => gg(G)
!! ff(ldx*ldy*ldz*ndat)=The array to be transformed.
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator.
!!   Defaults to ABI_FFTW_ESTIMATE.
!!
!! OUTPUT
!! gg(ldx*ldy*ldz*ndat)=The FFT of ff.
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_c2c_op_spc(nx,ny,nz,ldx,ldy,ldz,ndat,isign,ff,gg,fftw_flags)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,isign,ndat
 integer,optional,intent(in) :: fftw_flags
!arrays
 complex(spc),intent(in) :: ff(ldx*ldy*ldz*ndat)
 complex(spc),intent(out) :: gg(ldx*ldy*ldz*ndat)

#ifdef HAVE_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: rank3=3,nt_all=-1
 integer :: my_flags,dist,stride
 integer(KIND_FFTW_PLAN) :: my_plan
!arrays
 integer :: embed(rank3),n(rank3)

! *************************************************************************

 my_flags=ABI_FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags= fftw_flags

 stride = 1
 dist   = ldx*ldy*ldz
 embed  = (/ldx,ldy,ldz/)
 n      = (/nx ,ny ,nz/)

 my_plan = fftw3_plan_many_dft(rank3, n, ndat, ff, embed, stride, dist, gg, embed, stride, dist, isign, my_flags, nt_all)

 ! Now perform the 3D FFT via FFTW.
 call sfftw_execute_dft(my_plan, ff, gg)

 call fftw3_destroy_plan(my_plan)

 if (isign==ABI_FFTW_FORWARD) then ! -1, FFTW returns not normalized FTs
   call xscal(ldx*ldy*ldz*ndat, REAL(one/(nx*ny*nz), KIND=sp), gg, 1)
 end if

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,isign/))
 ABI_UNUSED(ff)
 ABI_UNUSED(gg)
 if (PRESENT(fftw_flags)) then
   ABI_UNUSED(fftw_flags)
 end if
#endif

end subroutine fftw3_c2c_op_spc
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_c2c_op_dpc
!! NAME
!!  fftw3_c2c_op_dpc
!!
!! FUNCTION
!! Driver routine for out-of-place 3D complex-complex FFT of lengths nx, ny, nz.
!! TARGET: single precision complex arrays
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimensions of the array.
!! ndat=Number of FFTs to be done.
!! isign= +1 : ff(G) => gg(R); -1 : ff(R) => gg(G)
!! ff(ldx*ldy*ldz*ndat)=The array to be transformed.
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator.
!!   Defaults to ABI_FFTW_ESTIMATE.
!!
!! OUTPUT
!! gg(ldx*ldy*ldz*ndat)=The FFT of ff.
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_c2c_op_dpc(nx,ny,nz,ldx,ldy,ldz,ndat,isign,ff,gg,fftw_flags)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,isign,ndat
 integer,optional,intent(in) :: fftw_flags
!arrays
 complex(dpc),intent(in) :: ff(ldx*ldy*ldz*ndat)
 complex(dpc),intent(out) :: gg(ldx*ldy*ldz*ndat)

#ifdef HAVE_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: rank3=3,nt_all=-1
 integer :: my_flags,dist,stride
 integer(KIND_FFTW_PLAN) :: my_plan
!arrays
 integer :: embed(rank3),n(rank3)

! *************************************************************************

 my_flags=ABI_FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags= fftw_flags

 stride = 1
 dist   = ldx*ldy*ldz
 embed  = (/ldx,ldy,ldz/)
 n      = (/nx ,ny ,nz/)

 my_plan = fftw3_plan_many_dft(rank3, n, ndat, ff, embed, stride, dist, gg, embed, stride, dist, isign, my_flags, nt_all)

 ! Now perform the 3D FFT via FFTW.
 call dfftw_execute_dft(my_plan, ff, gg)

 call fftw3_destroy_plan(my_plan)

 if (isign==ABI_FFTW_FORWARD) then ! -1, FFTW returns not normalized FTs
   call xscal(ldx*ldy*ldz*ndat, one/(nx*ny*nz), gg, 1)
 end if

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,isign/))
 ABI_UNUSED(ff)
 ABI_UNUSED(gg)
 if (PRESENT(fftw_flags)) then
   ABI_UNUSED(fftw_flags)
 end if
#endif

end subroutine fftw3_c2c_op_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_r2c_op
!! NAME
!!  fftw3_r2c_op
!!
!! FUNCTION
!! Driver routine for out-of-place 3D real-to-complex FFT of lengths nx, ny, nz.
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimensions of the f array (to avoid cache conflicts).
!! ff(ldx*ldy*ldz*ndat)=The real array to be transformed.
!! ndat=Number of FFTs to be done.
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator.
!!   Defaults to ABI_FFTW_ESTIMATE.
!!
!! OUTPUT
!! gg(2,nx*ny*nz*ndat)=The forward FFT of ff.
!!
!! NOTES
!!  FIXME For the time-being. No augmentation of the mesh to reduce memory conflicts, as MKL crashes
!!  if the advanced interface is used.
!!
!! PARENTS
!!      m_fftw3
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_r2c_op(nx,ny,nz,ldx,ldy,ldz,ndat,ff,gg,fftw_flags)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat
 integer,optional,intent(in) :: fftw_flags
!arrays
 real(dp),intent(in) :: ff(ldx*ldy*ldz*ndat)
 real(dp),intent(out) :: gg(2,ldx*ldy*ldz*ndat)

#ifdef HAVE_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: rank3=3,nt_all=-1
 integer :: nhp,my_flags,idist,odist,padx,i1,i2,i3,igp,igf,imgf,stride
 integer :: i1inv,i2inv,i3inv,idat,padatf
 integer(KIND_FFTW_PLAN) :: my_plan
!arrays
 integer :: inembed(rank3),onembed(rank3),n(rank3)
 integer,allocatable :: i1inver(:),i2inver(:),i3inver(:)
 real(dp),allocatable :: gg_hp(:,:)

! *************************************************************************

 my_flags=ABI_FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags= fftw_flags

 idist = ldx*ldy*ldz
 nhp = (nx/2+1)*ny*nz
 odist = nhp

 stride = 1
 n      = (/nx,ny,nz/)
 inembed= (/ldx,ldy,ldz/)
 onembed= (/(nx/2+1),ny,nz/)

 ABI_MALLOC(gg_hp,(2,nhp*ndat))

#ifdef DEV_RC_BUG
 if (ndat/=1) MSG_ERROR("ndat/=1 + MKL not coded")

 if (ANY( n /= inembed )) then
   MSG_ERROR("Augmentation not supported")
 end if

 call dfftw_plan_dft_r2c_3d(my_plan, nx, ny, nz, ff, gg_hp, my_flags)
 if (my_plan==NULL_PLAN) then
   MSG_ERROR("dfftw_plan_dft_r2c_3d returned NULL_PLAN")
 end if

 !fftw_plan fftw_plan_many_dft_r2c(int rank3, const int *n, int howmany,
 !  double *in, const int *inembed, int istride, int idist,
 !  fftw_complex *out, const int *onembed, int ostride, int odist, unsigned flags);
#else
 my_plan = dplan_many_dft_r2c(rank3, n, ndat, ff, inembed, stride, idist, gg_hp, onembed, stride, odist, my_flags, nt_all)
#endif

 ! Now perform the 3D FFT via FFTW. r2c are always ABI_FFTW_FORWARD
 call dfftw_execute_dft_r2c(my_plan, ff, gg_hp)

 call fftw3_destroy_plan(my_plan)

 call ZDSCAL(nhp*ndat, one/(nx*ny*nz), gg_hp, 1)  ! FFTW returns not normalized FTs
 ! Reconstruct full FFT: Hermitian redundancy: out[i] is the conjugate of out[n-i]
 padx = (nx/2+1)

 ABI_MALLOC(i1inver,(padx))
 ABI_MALLOC(i2inver,(ny))
 ABI_MALLOC(i3inver,(nz))

 i1inver(1)=1
 do i1=2,padx
   i1inver(i1)=nx+2-i1
 end do

 i2inver(1)=1
 do i2=2,ny
   i2inver(i2)=ny+2-i2
 end do

 i3inver(1)=1
 do i3=2,nz
   i3inver(i3)=nz+2-i3
 end do

 igp=0
 do idat=1,ndat
   padatf = (idat-1)*ldx*ldy*ldz
   do i3=1,nz
     i3inv = i3inver(i3)
     do i2=1,ny
       i2inv = i2inver(i2)
       do i1=1,padx
         igp = igp+1
         igf = i1 + (i3-1)*ldx*ldy + (i2-1)*ldx + padatf
         gg(:,igf) =  gg_hp(:,igp)
         i1inv = i1inver(i1)
         if (i1inv/=i1) then
           imgf = i1inv + (i3inv-1)*ldx*ldy + (i2inv-1)*ldx + padatf
           gg(1,imgf) =  gg_hp(1,igp)
           gg(2,imgf) = -gg_hp(2,igp)
         end if
       end do
     end do
   end do
 end do

 ABI_FREE(i1inver)
 ABI_FREE(i2inver)
 ABI_FREE(i3inver)

 ABI_FREE(gg_hp)

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz/))
 ABI_UNUSED(ff)
 ABI_UNUSED(gg(1,1))
 if (PRESENT(fftw_flags)) then
   ABI_UNUSED(fftw_flags)
 end if
#endif

end subroutine fftw3_r2c_op
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_c2r_op
!! NAME
!!  fftw3_c2r_op
!!
!! FUNCTION
!! Driver routine for out-of-place 3D complex-to-real FFT of lengths nx, ny, nz.
!!
!! INPUTS
!! nx,ny,nz=Number of point along the three directions.
!! ldx,ldy,ldz=Physical dimension of the f array (to avoid cache conflicts).
!! ndat=Number of FFTs to be done.
!! ff(2*ldx*ldy*ldz*ndat)=The complex array to be transformed.
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator.
!!   Defaults to ABI_FFTW_ESTIMATE.
!!
!! OUTPUT
!! gg(ldx*ldy*ldz*ndat)=The backwards real FFT of ff.
!!
!! NOTES
!!  FIXME For the time-being. No augmentation of the mesh to reduce memory conflicts, as MKL crashes
!!  if the advanced interface is used.
!!
!! PARENTS
!!      m_fftw3
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_c2r_op(nx,ny,nz,ldx,ldy,ldz,ndat,ff,gg,fftw_flags)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat
 integer,optional,intent(in) :: fftw_flags
!arrays
 real(dp),intent(in) :: ff(2,ldx*ldy*ldz*ndat)
 real(dp),intent(out) :: gg(ldx*ldy*ldz*ndat)

#ifdef HAVE_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: rank3=3,nt_all=-1
 integer :: nhp,my_flags,padx,i2,i3,igp,igf,idat,padatf,padatp,idist,odist,stride
 integer(KIND_FFTW_PLAN) :: my_plan
!arrays
 integer :: inembed(rank3),onembed(rank3),n(rank3)
 real(dp),allocatable :: ff_hp(:,:)

! *************************************************************************

#ifdef DEV_RC_BUG
 if (ANY( (/nx,ny,nz/) /= (/ldx,ldy,ldz/) )) then
   MSG_ERROR("Augmentation not supported")
 end if
#endif

 my_flags=ABI_FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags= fftw_flags

 stride  = 1
 nhp     = (nx/2+1)*ny*nz
 idist   = nhp
 odist   = ldx*ldy*ldz
 n       = (/nx,ny,nz/)
 inembed = (/(nx/2+1),ny,nz/)
 onembed = (/ldx,ldy,ldz/)

 ! Fill the Hermitian part: Hermitian redundancy: out[i] is the conjugate of out[n-i]
 ABI_MALLOC(ff_hp,(2,nhp*ndat))

 padx = (nx/2+1)
 do idat=1,ndat
   padatf=(idat-1)*ldx*ldy*ldz
   padatp=(idat-1)*padx*ny*nz
!$OMP PARALLEL DO PRIVATE(igf,igp)
   do i3=1,nz
     do i2=1,ny
       igf = (i3-1)*ldx*ldy + (i2-1)*ldx   + padatf
       igp = (i3-1)*padx*ny + (i2-1)*padx  + padatp
       ff_hp(:,igp+1:igp+padx) = ff(:,igf+1:igf+padx)
     end do
   end do
 end do

 ! NOTE: The c2r transform destroys its input array even for out-of-place transforms.
#ifdef DEV_RC_BUG
 if (ndat/=1) MSG_ERROR("ndat/=1 + MKL not coded")
 call dfftw_plan_dft_c2r_3d(my_plan, nx, ny, nz, ff_hp, gg, my_flags)
 if (my_plan==NULL_PLAN) then
   MSG_ERROR("dfftw_plan_dft_c2r_3d returned NULL_PLAN")
 end if
#else
 my_plan = dplan_many_dft_c2r(rank3, n, ndat, ff_hp, inembed, stride, idist, gg, onembed, stride, odist, my_flags, nt_all)
#endif

 ! Now perform the 3D FFT via FFTW. c2r are always ABI_FFTW_BACKWARD
 call dfftw_execute_dft_c2r(my_plan, ff_hp, gg)

 call fftw3_destroy_plan(my_plan)

 ABI_FREE(ff_hp)

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz/))
 ABI_UNUSED(ff(1,1))
 ABI_UNUSED(gg(1))
 if (PRESENT(fftw_flags)) then
   ABI_UNUSED(fftw_flags)
 end if
#endif

end subroutine fftw3_c2r_op
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_many_dft_op
!! NAME
!!  fftw3_many_dft_op
!!
!! FUNCTION
!! Driver routine for many out-of-place 3D complex-to-complex FFTs of lengths nx, ny, nz.
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimension of the fin and fout arrays (to avoid cache conflicts).
!! ndat=Number of FFTs to be done.
!! fin(2*ldx*ldy*ldz*ndat)=The complex array to be transformed.
!! isign=sign of Fourier transform exponent: current convention uses
!!   +1 for transforming from G to r,
!!   -1 for transforming from r to G.
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator.
!!   Defaults to ABI_FFTW_ESTIMATE.
!!
!! OUTPUT
!! fout(2,ldx*ldy*ldz*ndat)=The Fourier transform of fin.
!!
!! PARENTS
!!      m_fftw3
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_many_dft_op(nx,ny,nz,ldx,ldy,ldz,ndat,isign,fin,fout,fftw_flags)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,isign
 integer,optional,intent(in) :: fftw_flags
!arrays
 real(dp),intent(in) :: fin(2*ldx*ldy*ldz*ndat)
 real(dp),intent(out) :: fout(2*ldx*ldy*ldz*ndat)

#ifdef HAVE_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: rank3=3,nt_all=-1
 integer :: my_flags,dist,stride
 integer(KIND_FFTW_PLAN) :: my_plan
!arrays
 integer :: embed(rank3),n(rank3)

! *************************************************************************

 my_flags=ABI_FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags= fftw_flags

 stride = 1
 dist   = ldx*ldy*ldz
 embed  = (/ldx,ldy,ldz/)
 n      = (/nx ,ny ,nz /)

 my_plan = fftw3_plan_many_dft(rank3, n, ndat, fin, embed, stride, dist, fout, embed, stride, dist, isign, my_flags, nt_all)

 ! Now perform the 3D FFT via FFTW.
 call dfftw_execute_dft(my_plan, fin, fout)

 call fftw3_destroy_plan(my_plan)

 if (isign==ABI_FFTW_FORWARD) then ! -1, FFTW returns not normalized FTs
  call ZDSCAL(ldx*ldy*ldz*ndat, one/(nx*ny*nz), fout, 1)
  !call cg_zscal(ldx*ldy*ldz*ndat, (/one/(nx*ny*nz), zero/), fout)
 end if

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,ndat,isign/))
 if (PRESENT(fftw_flags)) then
   ABI_UNUSED(fftw_flags)
 end if
 ABI_UNUSED(fin(1))
 ABI_UNUSED(fout(1))
#endif

end subroutine fftw3_many_dft_op
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_many_dft_ip
!! NAME
!!  fftw3_many_dft_ip
!!
!! FUNCTION
!! Driver routine for many in-place 3D complex-to-complex FFTs of lengths nx, ny, nz.
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimension of the finout array (to avoid cache conflicts).
!! ndat=Number of FFTs to be done.
!! isign=sign of Fourier transform exponent: current convention uses
!!   +1 for transforming from G to r,
!!   -1 for transforming from r to G.
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator.
!!   Defaults to ABI_FFTW_ESTIMATE.
!!
!! OUTPUT
!! finout(2,ldx*ldy*ldz*ndat)=
!!   In input: The complex array to be transformed.
!!   In output: The FFT results.
!!
!! PARENTS
!!      m_fftw3
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_many_dft_ip(nx,ny,nz,ldx,ldy,ldz,ndat,isign,finout,fftw_flags)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,isign
 integer,optional,intent(in) :: fftw_flags
!arrays
 real(dp),intent(inout) :: finout(2*ldx*ldy*ldz*ndat)

#ifdef HAVE_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: rank3=3,nt_all=-1
 integer :: my_flags,dist,stride
 integer(KIND_FFTW_PLAN) :: my_plan
!arrays
 integer :: embed(rank3),n(rank3)

! *************************************************************************

 my_flags=ABI_FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags= fftw_flags

 stride = 1
 dist   = ldx*ldy*ldz
 embed  = (/ldx,ldy,ldz/)
 n      = (/nx ,ny ,nz /)

 my_plan = fftw3_plan_many_dft(rank3, n, ndat, finout, embed, stride, dist, finout,embed, stride, dist, isign, my_flags, nt_all)

 ! Now perform the 3D FFT via FFTW.
 call dfftw_execute_dft(my_plan, finout, finout)

 call fftw3_destroy_plan(my_plan)

 if (isign==ABI_FFTW_FORWARD) then ! -1, FFTW returns not normalized FTs
  call ZDSCAL(ldx*ldy*ldz*ndat, one/(nx*ny*nz), finout, 1)
  !call cg_zscal(ldx*ldy*ldz*ndat, (/one/(nx*ny*nz),zero/), finout)
 end if

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,ndat,isign/))
 if (PRESENT(fftw_flags)) then
   ABI_UNUSED(fftw_flags)
 end if
 ABI_UNUSED(finout(1))
#endif

end subroutine fftw3_many_dft_ip
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_cleanup
!! NAME
!!  fftw3_cleanup
!!
!! FUNCTION
!!  Reset FFTW to the pristine state it was in when you started your program,
!!  All existing plans become undefined.
!!
!! NOTES
!!  FFTW planner saves some other persistent data, such as the accumulated wisdom and a list of
!!  algorithms available in the current configuration. If you want to deallocate all of that and reset
!!  FFTW to the pristine state it was in when you started your program, you can call fftw3_cleanup();
!!  After calling fftw3_cleanup, all existing plans become undefined, and you should not attempt to
!!  execute them nor to destroy them. You can however create and execute/destroy new plans, in which case
!!  FFTW starts accumulating wisdom information again.
!!  fftw3_cleanup does not deallocate your plans, however. To prevent memory leaks, you must still call
!!  fftw_destroy_plan before executing fftw3_cleanup
!!
!! PARENTS
!!      m_driver
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_cleanup()

! *************************************************************************

#ifdef HAVE_FFTW3_MPI
 call fftw_mpi_cleanup()
#endif
#ifdef HAVE_FFTW3_THREADS
 if (THREADS_INITED==1) then
   call dfftw_cleanup_threads()
   THREADS_INITED = 0
 end if
#elif defined HAVE_FFTW3
 call dfftw_cleanup()
#else
 MSG_ERROR("FFTW3 support not activated")
#endif

end subroutine fftw3_cleanup
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_destroy_plan
!! NAME
!!  fftw3_destroy_plan
!!
!! FUNCTION
!!  Release the memory allocate for the plan.
!!
!! INPUTS
!!
!! PARENTS
!!      m_fftw3
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_destroy_plan(plan)

!Arguments ------------------------------------
!scalars
 integer(KIND_FFTW_PLAN),intent(in) :: plan

! *************************************************************************

#ifdef HAVE_FFTW3
!$OMP CRITICAL (OMPC_fftw3_destroy_plan)
 call dfftw_destroy_plan(plan)
!$OMP END CRITICAL (OMPC_fftw3_destroy_plan)

#else
 if (.FALSE.) write(std_out,*)plan
#endif

end subroutine fftw3_destroy_plan
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_init_threads
!! NAME
!!  fftw3_init_threads
!!
!! FUNCTION
!!  This function performs the one-time initialization required to use FFTW3 threads.
!!  It does nothing if HAVE_FFT_FFTW3_THREADS is not defined.
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!  The one-time initialization required to use FFTW3 threads is performed when the routine
!!  is called for the first time.
!!
!! PARENTS
!!      fftprof,m_driver
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_init_threads()

!Local variables ------------------------------
!scalars
#ifdef HAVE_FFTW3_THREADS
 integer :: iret
#endif

! *************************************************************************

#ifdef HAVE_FFTW3_THREADS
 if (THREADS_INITED==0) then
   !call wrtout(std_out,"Calling dfftw_init_threads()","COLL")
   call dfftw_init_threads(iret)

   if (iret==0) then
     MSG_WARNING(" dfftw_init_threads returned 0; threaded FFTW3 is not being used!")
   else
     THREADS_INITED=1
   end if
   call fftw3_set_nthreads()
 end if

#ifndef HAVE_OPENMP
  MSG_WARNING("Using FFTW3 with threads but HAVE_OPENMP is not defined!")
#endif
#endif

#ifdef HAVE_FFTW3_MPI
  !call wrtout(std_out,"Calling fftw_mpi_init()","COLL")
  call fftw_mpi_init()
#endif

end subroutine fftw3_init_threads
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_set_nthreads
!! NAME
!!  fftw3_set_nthreads
!!
!! FUNCTION
!!  This function sets the number of threads you want FFTW3 to use (or actually, the maximum number).
!!  It also performs any one-time initialization required to use FFTW3 threads.
!!  All plans subsequently created with any planner routine will use nthreads threads.
!!  If you pass an nthreads argument of 1 (the default), threads are disabled for subsequent plans.
!!  It does nothing if HAVE_FFT_FFTW3_THREADS is not defined.
!!
!! INPUTS
!!  [nthreads]=The number of threads you want FFTW3 to use.  Default xomp_get_max_threads()
!!
!! PARENTS
!!      m_fft_prof,m_fftw3
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_set_nthreads(nthreads)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: nthreads

!Local variables ------------------------------
!scalars
#ifdef HAVE_FFTW3_THREADS
 integer :: istat,nt
 integer,parameter :: enough=1
 integer,save :: nwarns=0
#endif

! *************************************************************************

#ifdef HAVE_FFTW3_THREADS
 if (THREADS_INITED==0) then
   MSG_WARNING("Threads are not initialized")
 end if

 if (PRESENT(nthreads)) then
   if (nthreads<=0) then
     nt = xomp_get_max_threads()
   else
     nt = nthreads
   end if
 else
   nt = xomp_get_max_threads()
 end if

 call dfftw_plan_with_nthreads(nt)

#ifndef HAVE_OPENMP
  if (nwarns <= enough) then
    nwarns = nwarns + 1
    MSG_WARNING("Using FFTW3 with threads but HAVE_OPENMP is not defined!")
  end if
#endif

#else
 if (PRESENT(nthreads)) then
   ABI_UNUSED(nthreads)
 end if
#endif

end subroutine fftw3_set_nthreads
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_fftpad_dp
!! NAME
!!  fftw3_fftpad_dp
!!
!! FUNCTION
!!  This routine transforms wavefunctions using 3D zero-padded FFTs with FFTW3.
!!  The 3D ffts are computed only on lines and planes which have non zero elements.
!!  These lines and planes are defined by the two vectors do_fft_x(ldy*nz) and do_fft_y(nz)
!!  FFT transform is in-place.
!!
!! INPUTS
!!   nx,ny,nz=Logical dimensions of the FFT mesh.
!!   ldx,ldy,ldz=Physical dimension of the f array (to avoid cache conflicts).
!!   ndat=Number of FFT transforms.
!!   mgfft=MAX(nx,ny,nz), only used to dimension gbound
!!   isign=The sign of the transform.
!!   gbound(2*mgfft+8,2)= The boundaries of the basis sphere of G vectors at a given k-point.
!!     See sphereboundary for more info.
!!
!! SIDE EFFECTS
!!   ff(2*ldx*ldy*ldz*ndat)=
!!     input: The array with the data to be transformed.
!!     output: The results of the FFT.
!!
!! PARENTS
!!      m_fftw3
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_fftpad_dp(ff,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign,gbound)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2)
 real(dp),intent(inout) :: ff(2*ldx*ldy*ldz*ndat)

!Local variables-------------------------------
!scalars
#ifdef HAVE_FFTW3
 integer,parameter :: dst=2
 real(dp) :: fact

! *************************************************************************

#include "fftw3_fftpad.finc"

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,mgfft,isign/))
 ABI_UNUSED(gbound(1,1))
 ABI_UNUSED(ff(1))
#endif

end subroutine fftw3_fftpad_dp
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_fftpad_dpc
!! NAME
!!  fftw3_fftpad_dpc
!!
!! FUNCTION
!!  This routine transforms wavefunctions using 3D zero-padded FFTs with FFTW3.
!!  The 3D ffts are computed only on lines and planes which have non zero elements.
!!  These lines and planes are defined by the two vectors do_fft_x(ldy*nz) and do_fft_y(nz)
!!  FFT transform is in-place. Target: complex arrays.
!!
!! INPUTS
!!   nx,ny,nz=Logical dimensions of the FFT mesh.
!!   ldx,ldy,ldz=Physical dimension of the f array (to avoid cache conflicts).
!!   ndat=Number of FFT transforms.
!!   mgfft=MAX(nx,ny,nz), only used to dimension gbound.
!!   isign=The sign of the transform.
!!   gbound(2*mgfft+8,2)= The boundaries of the basis sphere of G vectors at a given k-point.
!!     See sphereboundary for more info.
!!
!! SIDE EFFECTS
!!  ff(ldx*ldy*ldz*ndat)=
!!    input: The array with the data to be transformed.
!!    output: The results of the FFT.
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_fftpad_dpc(ff,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign,gbound)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2)
 complex(dpc),intent(inout) :: ff(ldx*ldy*ldz*ndat)

#ifdef HAVE_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: dst=1
 real(dp) :: fact

! *************************************************************************

#include "fftw3_fftpad.finc"

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign/))
 ABI_UNUSED(gbound(1,1))
 ABI_UNUSED(ff(1))
#endif

end subroutine fftw3_fftpad_dpc
!!***

#ifdef HAVE_FFTW3

!----------------------------------------------------------------------

!!****f* m_fftw3/dplan_many_dft_1D
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! SOURCE

function dplan_many_dft_1D(rank,n,howmany,fin,inembed,istride,idist,fout,onembed,ostride,odist,sign,flags,nthreads) result(plan)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: rank,howmany,istride,ostride, sign,flags,idist,odist,nthreads
 integer,intent(in) :: n(rank),inembed(rank),onembed(rank)
 integer(KIND_FFTW_PLAN) :: plan
!arrays
 real(dp) :: fin(*),fout(*)

!Local variables-------------------------------
 character(len=500) :: msg,frmt

! *************************************************************************

!$OMP CRITICAL (OMPC_dfftw_plan_many_dft_1D)
 call fftw3_set_nthreads(nthreads)

 call dfftw_plan_many_dft(plan, rank, n, howmany, &
&  fin, inembed, istride, idist, fout, onembed, ostride, odist, sign, flags)
!$OMP END CRITICAL (OMPC_dfftw_plan_many_dft_1D)

 if (plan==NULL_PLAN) then
   call wrtout(std_out,"dfftw_plan_many_dft returned NULL_PLAN!","COLL")
   write(frmt,*)"(a,",rank,"(1x,i0),3(a,i0),a,2(a,",rank,"(1x,i0),2(a,i0),a))"
   write(msg,frmt)&
&    " n= ",n," howmany= ",howmany," sign= ",sign," flags= ",flags,ch10,&
&    " inembed= ",inembed," istride= ",istride," idist=",idist,ch10,    &
&    " onembed= ",onembed," ostride= ",ostride," odist=",idist,ch10
   call wrtout(std_out,msg,"COLL")
   MSG_ERROR("Check FFTW library and/or abinit code")
 end if

end function dplan_many_dft_1D
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/dplan_many_dft_2D
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! SOURCE

function dplan_many_dft_2D(rank,n,howmany,fin,inembed,istride,idist,fout,onembed,ostride,odist,sign,flags,nthreads) result(plan)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: rank,howmany,istride,ostride, sign,flags,idist,odist,nthreads
 integer,intent(in) :: n(rank),inembed(rank),onembed(rank)
 integer(KIND_FFTW_PLAN) :: plan
!arrays
 real(dp) :: fin(2,*),fout(2,*)

!Local variables-------------------------------
 character(len=500) :: msg,frmt

! *************************************************************************

!$OMP CRITICAL (OMPC_dfftw_plan_many_dft_2D)
 call fftw3_set_nthreads(nthreads)

 call dfftw_plan_many_dft(plan, rank, n, howmany, &
&  fin, inembed, istride, idist, fout, onembed, ostride, odist, sign, flags)
!$OMP END CRITICAL (OMPC_dfftw_plan_many_dft_2D)

 if (plan==NULL_PLAN) then
   call wrtout(std_out,"dfftw_plan_many_dft returned NULL_PLAN!","COLL")
   write(frmt,*)"(a,",rank,"(1x,i0),3(a,i0),a,2(a,",rank,"(1x,i0),2(a,i0),a))"
   write(msg,frmt)&
&    " n= ",n," howmany= ",howmany," sign= ",sign," flags= ",flags,ch10,&
&    " inembed= ",inembed," istride= ",istride," idist=",idist,ch10,    &
&    " onembed= ",onembed," ostride= ",ostride," odist=",idist,ch10
   call wrtout(std_out,msg,"COLL")
   MSG_ERROR("Check FFTW library and/or abinit code")
 end if

end function dplan_many_dft_2D
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/cplan_many_dft
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! SOURCE
!! FIXME  technically it should be intent(inout) since FFTW3 can destroy the input for particular flags.

function cplan_many_dft(rank,n,howmany,fin,inembed,istride,idist,fout,onembed,ostride,odist,sign,flags,nthreads) result(plan)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: rank,howmany,istride,ostride, sign,flags,idist,odist,nthreads
 integer,intent(in) :: n(rank),inembed(rank),onembed(rank)
 integer(KIND_FFTW_PLAN) :: plan
!arrays
 complex(spc) :: fin(*),fout(*)

!Local variables-------------------------------
 character(len=500) :: msg,frmt

! *************************************************************************

!$OMP CRITICAL (OMPC_cplan_many_dft)
 call fftw3_set_nthreads(nthreads)

 call sfftw_plan_many_dft(plan, rank, n, howmany, &
&  fin, inembed, istride, idist, fout, onembed, ostride, odist, sign, flags)
!$OMP END CRITICAL (OMPC_cplan_many_dft)

 if (plan==NULL_PLAN) then ! handle the error
   call wrtout(std_out,"sfftw_plan_many_dft returned NULL_PLAN (complex version)","COLL")
   write(frmt,*)"(a,",rank,"(1x,i0),3(a,i0),a,2(a,",rank,"(1x,i0),2(a,i0),a))"
   write(msg,frmt)&
&    " n = ",n," howmany = ",howmany," sign = ",sign," flags = ",flags,ch10,&
&    " inembed = ",inembed," istride = ",istride," idist =",idist,ch10,     &
&    " onembed = ",onembed," ostride = ",ostride," odist =",idist,ch10
   call wrtout(std_out,msg,"COLL")
   MSG_ERROR("Check FFTW library and/or abinit code")
 end if

end function cplan_many_dft
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/zplan_many_dft
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! SOURCE
!! FIXME  technically it should be intent(inout) since FFTW3 can destroy the input for particular flags.

function zplan_many_dft(rank,n,howmany,fin,inembed,istride,idist,fout,onembed,ostride,odist,sign,flags,nthreads) result(plan)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: rank,howmany,istride,ostride, sign,flags,idist,odist,nthreads
 integer,intent(in) :: n(rank),inembed(rank),onembed(rank)
 integer(KIND_FFTW_PLAN) :: plan
!arrays
 complex(dpc) :: fin(*),fout(*)

!Local variables-------------------------------
 character(len=500) :: msg,frmt

! *************************************************************************

!$OMP CRITICAL (OMPC_zplan_many_dft)
 call fftw3_set_nthreads(nthreads)

 call dfftw_plan_many_dft(plan, rank, n, howmany, &
&  fin, inembed, istride, idist, fout, onembed, ostride, odist, sign, flags)
!$OMP END CRITICAL (OMPC_zplan_many_dft)

 if (plan==NULL_PLAN) then ! handle the error
   call wrtout(std_out,"dfftw_plan_many_dft returned NULL_PLAN (complex version)","COLL")
   write(frmt,*)"(a,",rank,"(1x,i0),3(a,i0),a,2(a,",rank,"(1x,i0),2(a,i0),a))"
   write(msg,frmt)&
&    " n = ",n," howmany = ",howmany," sign = ",sign," flags = ",flags,ch10,&
&    " inembed = ",inembed," istride = ",istride," idist =",idist,ch10,     &
&    " onembed = ",onembed," ostride = ",ostride," odist =",idist,ch10
   call wrtout(std_out,msg,"COLL")
   MSG_ERROR("Check FFTW library and/or abinit code")
 end if

end function zplan_many_dft
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/dplan_many_dft_r2c
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! SOURCE
!! FIXME  technically it should be intent(inout) since FFTW3 can destroy the input
!! for particular flags.

function dplan_many_dft_r2c(rank,n,howmany,fin,inembed,istride,idist,fout,onembed,ostride,odist,flags,nthreads) result(plan)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: rank,howmany,istride,ostride,flags,idist,odist,nthreads
 integer,intent(in) :: n(rank),inembed(rank),onembed(rank)
 integer(KIND_FFTW_PLAN) :: plan
!arrays
 real(dp) :: fin(*),fout(*)

!Local variables-------------------------------
 character(len=500) :: msg,frmt

! *************************************************************************

!$OMP CRITICAL (OMPC_dplan_many_dft_r2c)
 call fftw3_set_nthreads(nthreads)

 call dfftw_plan_many_dft_r2c(plan, rank, n, howmany, &
&  fin, inembed, istride, idist, fout, onembed, ostride, odist, flags)
!$OMP END CRITICAL (OMPC_dplan_many_dft_r2c)

 if (plan==NULL_PLAN) then ! handle the error.
   call wrtout(std_out,"dfftw_plan_many_dft_r2c returned NULL_PLAN","COLL")
   write(frmt,*)"(a,",rank,"(1x,i0),2(a,i0),a,2(a,",rank,"(1x,i0),2(a,i0),a))"
   write(msg,frmt)&
&    " n = ",n," howmany = ",howmany," flags = ",flags,ch10,&
&    " inembed = ",inembed," istride = ",istride," idist = ",idist,ch10,&
&    " onembed = ",onembed," ostride = ",ostride," odist = ",idist,ch10
   call wrtout(std_out,msg,"COLL")
   MSG_ERROR("Check FFTW library and/or abinit code")
 end if

end function dplan_many_dft_r2c
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/dplan_many_dft_c2r
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! SOURCE

function dplan_many_dft_c2r(rank,n,howmany,fin,inembed,istride,idist,fout,onembed,ostride,odist,flags, nthreads) result(plan)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: rank,howmany,istride,ostride,flags,idist,odist,nthreads
 integer,intent(in) :: n(rank),inembed(rank),onembed(rank)
 integer(KIND_FFTW_PLAN) :: plan
!arrays
 real(dp) :: fin(*),fout(*)

!Local variables-------------------------------
 character(len=500) :: msg,frmt

! *************************************************************************

!$OMP CRITICAL (OMPC_dplan_many_dft_c2r)
 call fftw3_set_nthreads(nthreads)

 call dfftw_plan_many_dft_c2r(plan, rank, n, howmany, &
&  fin, inembed, istride, idist, fout, onembed, ostride, odist, flags)
!$OMP END CRITICAL (OMPC_dplan_many_dft_c2r)

 if (plan==NULL_PLAN) then ! handle the error.
   call wrtout(std_out,"dfftw_plan_many_dft_c2r returned NULL_PLAN","COLL")
   write(frmt,*)"(a,",rank,"(1x,i0),2(a,i0),a,2(a,",rank,"(1x,i0),2(a,i0),a))"
   write(msg,frmt)&
&    " n = ",n," howmany = ",howmany," flags = ",flags,ch10,&
&    " inembed = ",inembed," istride = ",istride," idist = ",idist,ch10,&
&    " onembed = ",onembed," ostride = ",ostride," odist = ",idist,ch10
   call wrtout(std_out,msg,"COLL")
   MSG_ERROR("Check FFTW library and/or abinit code")
 end if

end function dplan_many_dft_c2r
!!***

#endif

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_execute_dft_dp
!! NAME
!! fftw3_execute_dft_dp
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!  This interface is used to perform complex to complex FFT with real arrays
!!  containing the real and imaginary part. I have to admit that this interface
!!  is a bit ambiguous since FFTW3 provides routines for real-to-real transforms.
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

#ifdef HAVE_FFTW3

subroutine fftw3_execute_dft_dp(plan, in, out)

!Arguments ------------------------------------
!scalars
 integer(KIND_FFTW_PLAN),intent(in) :: plan
 real(C_DOUBLE),intent(inout) :: in(*)
 real(C_DOUBLE),intent(out) :: out(*)

! *************************************************************************

 call dfftw_execute_dft(plan, in, out)

end subroutine fftw3_execute_dft_dp
!!***

#endif

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_execute_dft_spc
!! NAME
!! fftw3_execute_dft_spc
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

#ifdef HAVE_FFTW3

subroutine fftw3_execute_dft_spc(plan, in, out)

!Arguments ------------------------------------
!scalars
 integer(KIND_FFTW_PLAN),intent(in) :: plan
 complex(C_FLOAT_COMPLEX),intent(inout) :: in(*)
 complex(C_FLOAT_COMPLEX),intent(out) :: out(*)

! *************************************************************************

 call sfftw_execute_dft(plan, in, out)

end subroutine fftw3_execute_dft_spc
!!***

#endif

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_execute_dft_dpc
!! NAME
!! fftw3_execute_dft_dpc
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

#ifdef HAVE_FFTW3

subroutine fftw3_execute_dft_dpc(plan, in, out)

!Arguments ------------------------------------
!scalars
 integer(KIND_FFTW_PLAN),intent(in) :: plan
 complex(C_DOUBLE_COMPLEX),intent(inout) :: in(*)
 complex(C_DOUBLE_COMPLEX),intent(out) :: out(*)

! *************************************************************************

 call dfftw_execute_dft(plan, in, out)

end subroutine fftw3_execute_dft_dpc
!!***

#endif

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_alloc_real1d_dp
!! NAME
!! fftw3_alloc_real1d_dp
!!
!! FUNCTION
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

#ifdef HAVE_FFTW3

subroutine fftw3_alloc_real1d_dp(size,cptr,fptr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: size
 real(dp),ABI_CONTIGUOUS pointer :: fptr(:)
 type(C_PTR),intent(out) :: cptr

! *************************************************************************

 cptr = fftw_malloc( INT(size*C_DOUBLE, KIND=C_SIZE_T))
 if (.not. C_ASSOCIATED(cptr)) then
   MSG_ERROR("fftw_malloc returned NULL!")
 end if

 call c_f_pointer(cptr, fptr, [size])

end subroutine fftw3_alloc_real1d_dp
!!***

#endif

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_alloc_real2d_dp
!! NAME
!! fftw3_alloc_real2d_dp
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

#ifdef HAVE_FFTW3

subroutine fftw3_alloc_real2d_dp(shape,cptr,fptr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: shape(2)
 real(dp),ABI_CONTIGUOUS pointer :: fptr(:,:)
 type(C_PTR),intent(out) :: cptr

! *************************************************************************

 cptr = fftw_malloc( INT(product(shape)*C_DOUBLE, KIND=C_SIZE_T))
 if (.not. C_ASSOCIATED(cptr)) then
   MSG_ERROR("fftw_malloc returned NULL!")
 end if

 call c_f_pointer(cptr, fptr, shape)

end subroutine fftw3_alloc_real2d_dp
!!***

#endif

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_alloc_complex1d_spc
!! NAME
!! fftw3_alloc_complex1d_spc
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

#ifdef HAVE_FFTW3

subroutine fftw3_alloc_complex1d_spc(size,cptr,fptr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: size
 complex(spc),ABI_CONTIGUOUS pointer :: fptr(:)
 type(C_PTR),intent(out) :: cptr

! *************************************************************************

 cptr = fftw_malloc( INT(2*size*C_FLOAT, KIND=C_SIZE_T))
 if (.not. C_ASSOCIATED(cptr)) then
   MSG_ERROR("fftw_malloc returned NULL!")
 end if

 call c_f_pointer(cptr, fptr, [size])

end subroutine fftw3_alloc_complex1d_spc
!!***

#endif

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_alloc_complex1d_dpc
!! NAME
!! fftw3_alloc_complex1d_dpc
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

#ifdef HAVE_FFTW3

subroutine fftw3_alloc_complex1d_dpc(size,cptr,fptr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: size
 complex(dpc),ABI_CONTIGUOUS pointer :: fptr(:)
 type(C_PTR),intent(out) :: cptr

! *************************************************************************

 cptr = fftw_malloc( INT(2*size*C_DOUBLE, KIND=C_SIZE_T))
 if (.not. C_ASSOCIATED(cptr)) then
   MSG_ERROR("fftw_malloc returned NULL!")
 end if

 call c_f_pointer(cptr, fptr, [size])

end subroutine fftw3_alloc_complex1d_dpc
!!***

#endif

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_spawn_threads_here
!! NAME
!! fftw3_spawn_threads_here
!!
!! FUNCTION
!!  Helper function that returns true if FFT calls should be OMP
!!  parallelized in the client code.
!!
!! INPUTS
!!  ndat=Number of FFT transforms to do
!!  nthreads = Number of threads available
!!
!! PARENTS
!!
!! SOURCE

function fftw3_spawn_threads_here(ndat,nthreads) result(ans)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndat,nthreads
 logical :: ans

! *************************************************************************

 ans = .FALSE.
#ifdef HAVE_OPENMP
 ans = (nthreads > 1 .and. MOD(ndat,nthreads) == 0 .and. .not. USE_LIB_THREADS)
#else
 ABI_UNUSED((/ndat,nthreads/))
#endif

end function fftw3_spawn_threads_here
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_use_lib_threads
!! NAME
!! fftw3_use_lib_threads
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_use_lib_threads(logvar)

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: logvar

! *************************************************************************

 USE_LIB_THREADS = logvar

end subroutine fftw3_use_lib_threads
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftwmpi_get_work_array
!! NAME
!!  fftwmpi_get_work_array
!!
!! FUNCTION
!! Driver routine for allocate fftw work arrray for 3D complex-to-complex FFTs of lengths nx, ny, nz.
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ndat=Number of FFTs to be done.
!! comm_fft=MPI communicator.
!!
!! OUTPUT
!! cdata_f,cdata_r: C pointers to use for fourier andreal data
!! n0,n0_tr : local size on the shared dimension (nz or ny if transposed mode is used)
!! offset,offset_tr : offset per process in continuous tabx
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftwmpi_get_work_array(cdata_f,cdata_r,rank,nx,ny,nz,ndat,comm_fft,n0,offset,n0_tr,offset_tr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ndat,rank,comm_fft
 integer(C_INTPTR_T), intent(out) :: n0, offset, n0_tr, offset_tr
 type(C_PTR), intent(out) :: cdata_f,cdata_r

!Local variables-------------------------------
#ifdef HAVE_FFTW3_MPI
!scalars
 integer(C_INTPTR_T) :: alloc_local
!arrays
 integer(C_INTPTR_T) :: fft_sizes(4)

! *************************************************************************

 ! Dimensions are inverted here (C interface).
 fft_sizes(1)=nz
 fft_sizes(2)=ny
 fft_sizes(3)=nx
 fft_sizes(4)=ndat

 alloc_local = fftw_mpi_local_size_many_transposed(rank,fft_sizes(1:3),fft_sizes(4), &
&      FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, comm_fft, &
&      n0,offset, &
&      n0_tr,offset_tr)

 cdata_f = fftw_alloc_complex(alloc_local)
 cdata_r = fftw_alloc_complex(alloc_local)

#else
  MSG_ERROR("FFTW3_MPI support not activated")
  ABI_UNUSED((/nx,ny,nz,ndat,rank,comm_fft/))
  cdata_f = C_NULL_PTR; cdata_r = C_NULL_PTR
  n0 = 0; offset = 0; n0_tr = 0; offset_tr = 0
#endif

end subroutine fftwmpi_get_work_array
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftwmpi_free_work_array
!! NAME
!!  fftwmpi_free_work_array
!!
!! FUNCTION
!!  routine for freeing fftw work arrray
!!
!! INPUTS
!!
!! OUTPUT
!! cdata_f,cdata_r: C pointers to free for fourier andreal data
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftwmpi_free_work_array(cdata_f,cdata_r)

!Arguments ------------------------------------
!scalars
 type(C_PTR), intent(inout) :: cdata_f,cdata_r

! *************************************************************************

#ifdef HAVE_FFTW3_MPI
 call fftw_free(cdata_r)
 call fftw_free(cdata_f)
#else
 MSG_ERROR("FFTW3_MPI support not activated")
 if(.false.) then
   cdata_r = C_NULL_PTR; cdata_f = C_NULL_PTR
 end if
#endif

end subroutine fftwmpi_free_work_array
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3mpi_many_dft_ip
!! NAME
!!  fftw3mpi_many_dft_ip
!!
!! FUNCTION
!! Driver routine for many out-of-place 3D complex-to-complex FFTs of lengths nx, ny, nz.
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimension of the fin and fout arrays (to avoid cache conflicts).
!! ndat=Number of FFTs to be done.
!! fin(2*ldx*ldy*ldz*ndat)=The complex array to be transformed.
!! isign=sign of Fourier transform exponent: current convention uses
!!   +1 for transforming from G to r,
!!   -1 for transforming from r to G.
!! comm_fft=MPI communicator for the FFT
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator.
!!   Defaults to ABI_FFTW_ESTIMATE.
!!
!! TODO
!!  Remove me
!!
!! OUTPUT
!! fout(2,ldx*ldy*ldz*ndat)=The Fourier transform of fin.
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3mpi_many_dft_ip(nx,ny,nz,ldx,ldy,ldz,ndat,isign,fin,fout,comm_fft,fftw_flags)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,isign,comm_fft
 integer,optional,intent(in) :: fftw_flags
!arrays
 real(dp),intent(in) :: fin(2,ldx,ldy,ldz*ndat)
 real(dp),intent(out) :: fout(2,ldx,ldy,ldz*ndat)

#ifdef HAVE_FFTW3_MPI
!Local variables-------------------------------
!scalars
 integer,parameter :: rank3=3
 integer :: my_flags
 real(dp):: factor_fft
!arrays
 type(C_PTR) :: plan, cdata
 complex(C_DOUBLE_COMPLEX), ABI_CONTIGUOUS pointer :: data(:,:,:)
 integer(C_INTPTR_T) :: i, j, k, alloc_local, local_n0, local_0_start,fft_sizes(4)

!*************************************************************************

 my_flags=ABI_FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags= fftw_flags

 ! get local data size and allocate (note dimension reversal)
 fft_sizes = [nz,ny,nx,ndat]

 alloc_local = fftw_mpi_local_size_many( &
&      rank3,fft_sizes(1:3),fft_sizes(4),&
&      FFTW_MPI_DEFAULT_BLOCK, comm_fft, &
&      local_n0,local_0_start)

 ! Allocate cdata, build the plane and copy data: fin --> data
 cdata = fftw_alloc_complex(alloc_local)
 call c_f_pointer(cdata, data, [fft_sizes(3),fft_sizes(2), local_n0])

 plan = fftw_mpi_plan_many_dft(rank3,fft_sizes(1:3),fft_sizes(4), &
&                               FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, &
&                               data,data,comm_fft,isign,my_flags)

 do k=1, local_n0*ndat
    do j=1, ny
       do i=1, nx
          data(i,j,k) = CMPLX( fin(1,i,j,k),fin(2,i,j,k),C_DOUBLE_COMPLEX)
       end do
    end do
 end do

 ! Compute transform.
 call fftw_mpi_execute_dft(plan, data, data)

 if(isign==ABI_FFTW_FORWARD) then
    ! Scale results.
    factor_fft = one / (nx*ny*nz)
    do k=1, local_n0*ndat
       do j=1, ny
          do i=1, nx
             fout(1,i,j,k) =  real(data(i,j,k)) * factor_fft
             fout(2,i,j,k) = aimag(data(i,j,k)) * factor_fft
          end do
       end do
    end do
 end if

 call fftw_destroy_plan(plan)
 call fftw_free(cdata)

#else
 MSG_ERROR("FFTW3_MPI support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,ndat,isign/))
 ABI_UNUSED(comm_fft)
 if (PRESENT(fftw_flags)) then
    ABI_UNUSED(fftw_flags)
 end if
 ABI_UNUSED(fin(1,1,1,1))
 ABI_UNUSED(fout(1,1,1,1))
#endif

end subroutine fftw3mpi_many_dft_ip
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3mpi_many_dft_tr
!! NAME
!!  fftw3mpi_many_dft_tr
!!
!! FUNCTION
!! Driver routine for many out-of-place 3D complex-to-complex FFTs of lengths nx, ny, nz.
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimension of the fin and fout arrays (to avoid cache conflicts).
!! ndat=Number of FFTs to be done.
!! fin(2*ldx*ldy*ldz*ndat)=The complex array to be transformed.
!! isign=sign of Fourier transform exponent: current convention uses
!!   +1 for transforming from G to r,
!!   -1 for transforming from r to G.
!! comm_fft=MPI communicator for the FFT.
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator.
!!   Defaults to ABI_FFTW_ESTIMATE.
!!
!! TODO
!!  Remove me
!!
!! OUTPUT
!! fout(2,ldx*ldy*ldz*ndat)=The Fourier transform of fin.
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3mpi_many_dft_tr(nx,ny,nz,ndat,isign,fin,fout,comm_fft,fftw_flags)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ndat,isign,comm_fft
 integer,optional,intent(in) :: fftw_flags
!arrays
 complex(C_DOUBLE_COMPLEX),ABI_CONTIGUOUS pointer  :: fin(:,:,:)
 complex(C_DOUBLE_COMPLEX),ABI_CONTIGUOUS pointer :: fout(:,:,:)

!Local variables-------------------------------
#ifdef HAVE_FFTW3_MPI
!scalars
 integer :: my_flags
 !FFTWMPI stuff
 type(C_PTR) :: plan
 integer(C_INTPTR_T) :: fft_sizes(4)

!*************************************************************************

 my_flags=ABI_FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags= fftw_flags
 my_flags = ior(my_flags,FFTW_DESTROY_INPUT)

 fft_sizes(1)=nz
 fft_sizes(2)=ny
 fft_sizes(3)=nx
 fft_sizes(4)=ndat

 plan = fftw_mpi_plan_many_dft(3,fft_sizes(1:3),fft_sizes(4), &
&                              FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, &
&                              fin,fout,comm_fft,isign,my_flags)

!Compute transform (as many times as desired)
 call fftw_mpi_execute_dft(plan, fin, fout)
 call fftw_destroy_plan(plan)

#else
 MSG_ERROR("FFTW3_MPI support not activated")
 ABI_UNUSED((/nx,ny,nz,ndat,isign,comm_fft/))
 if (PRESENT(fftw_flags)) then
    ABI_UNUSED(fftw_flags)
 end if
 ABI_UNUSED(fin(1,1,1))
 ABI_UNUSED(fout(1,1,1))
#endif

end subroutine fftw3mpi_many_dft_tr
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_mpifourdp_c2r
!! NAME
!!  fftw3_mpifourdp_c2r
!!
!! FUNCTION
!! Driver routine for transposed out-of-place 3D complex-to-real FFT of lengths nx, ny, nz.
!!
!! INPUTS
!! nx,ny,nz=Number of point along the three directions.
!! ndat=Number of FFTs to be done.
!! fofg(2,nx*ny*nz*ndat)=The complex array to be transformed.
!! comm_fft=MPI communicator.
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator.
!!   Defaults to ABI_FFTW_ESTIMATE.
!!
!! OUTPUT
!! fofr(2,nx*ny*nz*ndat)=The backwards real FFT of ff.
!!
!! NOTES
!! LOCAL DATA IN FOURIER SPACE : TRANSPOSED ORDER
!! real space     --> dim = [  nx  | ny | nz/np_fft]
!! fourier  space --> dim = [ nx/2 | nz | ny/np_ff ]
!!
!! PARENTS
!!      m_fftw3
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_mpifourdp_c2r(nfft,ngfft,ndat,&
  fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft,fftw_flags)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,ndat,comm_fft
 integer,optional,intent(in) :: fftw_flags
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(in) :: fftn2_distrib(ngfft(2)),ffti2_local(ngfft(2))
 integer,intent(in) :: fftn3_distrib(ngfft(3)),ffti3_local(ngfft(3))
 real(dp),intent(in) :: fofg(2,nfft*ndat)
 real(dp),intent(out) :: fofr(nfft*ndat)

!Local variables-------------------------------
#ifdef HAVE_FFTW3_MPI
!scalars
 integer,parameter :: rank3=3
 integer :: nx,ny,nz,nproc_fft
 type(C_PTR) :: plan_bw, cdata_cplx,cdata_real
 integer(C_INTPTR_T) :: i,j,jdat,k,alloc_local,fft_sizes(4),demi_nx,base,idat,kdat
 integer(C_INTPTR_T) :: local_n0, local_0_start, local_n1, local_1_start
!arrays
 complex(C_DOUBLE_COMPLEX), ABI_CONTIGUOUS pointer :: data_cplx(:,:,:)
 real(C_DOUBLE), ABI_CONTIGUOUS pointer :: data_real(:,:,:)

! *************************************************************************

 !ABI_CHECK(ndat==1, "ndat > 1 not implemented yet")

 nx=ngfft(1); ny=ngfft(2); nz=ngfft(3)
 nproc_fft = xmpi_comm_size(comm_fft)

 demi_nx = nx/2 + 1
 fft_sizes(1)=nz
 fft_sizes(2)=ny
 fft_sizes(3)=demi_nx
 fft_sizes(4)=ndat

 alloc_local = fftw_mpi_local_size_many_transposed(&
&      rank3,fft_sizes(1:3),fft_sizes(4), &
&      FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, comm_fft, &
&      local_n0,local_0_start, &
&      local_n1,local_1_start)

 cdata_cplx = fftw_alloc_complex(alloc_local)
 cdata_real = fftw_alloc_real(alloc_local*2)

! OLD BY FDHAM
 ! dimensions are  (x/2,z,y) in Fourier's Space
 call c_f_pointer(cdata_cplx, data_cplx, [demi_nx  ,fft_sizes(1),local_n1])
 ! dimensions in real space : (nx,ny,nz/nproc)
 call c_f_pointer(cdata_real, data_real, [2*demi_nx,fft_sizes(2),local_n0])

 ! dimensions are  (x/2,z,y) in Fourier's Space
 !call c_f_pointer(cdata_cplx, data_cplx, [demi_nx  ,fft_sizes(1),local_n0])

 !! dimensions in real space : (nx,ny,nz/nproc)
 !call c_f_pointer(cdata_real, data_real, [2*demi_nx,fft_sizes(2),local_n1])

 fft_sizes(3)=nx
 plan_bw =  fftw_mpi_plan_many_dft_c2r(&
&      rank3,fft_sizes(1:3),fft_sizes(4), &
&      FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, &
&      data_cplx, data_real , &
&      comm_fft,ior(ABI_FFTW_ESTIMATE,ABI_FFTW_MPI_TRANSPOSED_IN))

 do idat=1,ndat
   do k=1, nz
     do j=1, ny / nproc_fft
       jdat = j + (idat-1) * ny / nproc_fft
       base = nx*((j-1) + (ny/nproc_fft)*(k-1)) + (idat-1) * nfft
       do i=1, demi_nx
         data_cplx(i,k,jdat) = CMPLX(fofg(1, i + base), fofg(2, i + base), kind=C_DOUBLE_COMPLEX)
       end do
     end do
   end do
 end do

 ! compute transform (as many times as desired)
 call fftw_mpi_execute_dft_c2r(plan_bw, data_cplx, data_real)

 do idat=1,ndat
   do k=1,local_n0
     kdat = k + (idat - 1) * local_n0
     do j=1,ny
       base = nx*((j-1) + ny*(k-1)) + (idat - 1) * nfft
       do i=1,nx
         fofr(i+base) = data_real(i,j,kdat)
       end do
     end do
   end do
 end do

 call fftw_destroy_plan(plan_bw)
 call fftw_free(cdata_cplx)
 call fftw_free(cdata_real)

#else
 MSG_ERROR("FFTW3_MPI support not activated")
 ABI_UNUSED((/nfft,ngfft(1),ndat,comm_fft/))
 ABI_UNUSED((/fftn2_distrib(1),ffti2_local(1)/))
 ABI_UNUSED((/fftn3_distrib(1),ffti3_local(1)/))
 if (PRESENT(fftw_flags)) then
    ABI_UNUSED(fftw_flags)
 end if
 ABI_UNUSED(fofg(1,1))
 ABI_UNUSED(fofr(1))
#endif

end subroutine fftw3_mpifourdp_c2r
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_mpifourdp_r2c
!! NAME
!!  fftw3_mpifourdp_r2c
!!
!! FUNCTION
!! Driver routine for out-of-place 3D real-to-complex FFT of lengths nx, ny, nz.
!!
!! INPUTS
!! fofr(nx*ny*nz*ndat)=The real array to be transformed.
!! ndat=Number of FFTs to be done.
!! comm_fft=MPI communicator for the FFT.
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator.
!!   Defaults to ABI_FFTW_ESTIMATE.
!!
!! OUTPUT
!! fofg(2,nx*ny*nz*ndat)=The forward FFT of ff.
!!
!! NOTES
!! LOCAL DATA FOR FOURIER TRANSFORMS : TRANSPOSED ORDER AND DISTRIBUTED
!! real space     --> dim = [  nx  | ny | nz/np_fft ]
!! fourier  space --> dim = [  nx | nz | ny/np_fft ]
!! we can't take in account the symetric of the real case because after
!! fft have been computed, the symetric data needed are dispatched over
!! other process in parallel
!!
!! PARENTS
!!      m_fftw3
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE


subroutine fftw3_mpifourdp_r2c(nfft,ngfft,ndat,&
  fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft,fftw_flags)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,ndat,comm_fft
 integer,optional,intent(in) :: fftw_flags
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(in) :: fftn2_distrib(ngfft(2)),ffti2_local(ngfft(2))
 integer,intent(in) :: fftn3_distrib(ngfft(3)),ffti3_local(ngfft(3))
 real(dp),intent(in) :: fofr(nfft*ndat)
 real(dp),intent(out) :: fofg(2,nfft*ndat)

!Local variables-------------------------------
#ifdef HAVE_FFTW3_MPI
 !scalars
 integer,parameter :: rank3=3
 integer :: my_flags,nproc_fft,nx,ny,nz
 integer(C_INTPTR_T) :: i,j,k,base,alloc_local,i1,i2,i3,igf,idat,kdat,i2dat,padatf
 integer(C_INTPTR_T) :: local_n0,local_0_start,local_n1,local_1_start
 real(dp) :: factor_fft
 type(C_PTR) :: plan_fw,cdata_cplx,cdata_real
!arrays
 complex(C_DOUBLE_COMPLEX), ABI_CONTIGUOUS pointer :: data_cplx(:,:,:),data_real(:,:,:)
 integer(C_INTPTR_T) :: fft_sizes(4)

! *************************************************************************

 nproc_fft = xmpi_comm_size(comm_fft)

 nx=ngfft(1); ny=ngfft(2); nz=ngfft(3)

 my_flags=ABI_FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags= fftw_flags

 fft_sizes(1)=nz
 fft_sizes(2)=ny
 fft_sizes(3)=nx
 fft_sizes(4)=ndat

 ! Get parallel sizes
 alloc_local = fftw_mpi_local_size_many_transposed(&
&      rank3,fft_sizes(1:3),fft_sizes(4), &
&      FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, comm_fft, &
&      local_n0,local_0_start, &
&      local_n1,local_1_start)

 ! Allocate data and reference it

 ! local data in real space     --> dim = [nx | ny | nz/nproc_fft]
 cdata_real = fftw_alloc_complex(alloc_local)
 call c_f_pointer(cdata_real, data_real, [fft_sizes(3),fft_sizes(2),local_n0])

 ! local data in Fourier space --> dim = [nx | nz | ny/nproc_fft]
 cdata_cplx = fftw_alloc_complex(alloc_local)
 call c_f_pointer(cdata_cplx, data_cplx, [fft_sizes(3),fft_sizes(1),local_n1])

 ! TODO: Use true real to complex API!
 ! Create Plan C2C (nx,ny,nz)
 plan_fw =  fftw_mpi_plan_many_dft(&
&      rank3,fft_sizes(1:3),fft_sizes(4), &
&      FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, &
&      data_real, data_cplx , &
&      comm_fft,ABI_FFTW_FORWARD,ior(ABI_FFTW_ESTIMATE,ABI_FFTW_MPI_TRANSPOSED_OUT))

 ! Copy input data in correct format
 do idat=1,ndat
   do k=1,local_n0
     kdat = k + (idat-1) * local_n0
     do j=1, ny
       base = nx*((j-1) + ny*(k-1)) + (idat-1) * nfft
       do i=1, nx
         data_real(i,j,kdat) = CMPLX(fofr(i+base),zero, kind=C_DOUBLE_COMPLEX)
       end do
     end do
   end do
 end do

 ! Compute transform
 call fftw_mpi_execute_dft(plan_fw, data_real, data_cplx)

 factor_fft = one / (nx*ny*nz)

 do idat=1,ndat
    padatf=(idat-1)*nfft
    do i3=1,nz
       do i2=1,ny/nproc_fft ! equivalent a local_n1
          i2dat = i2 + (idat-1) * ny/nproc_fft
          do i1=1,nx
             igf = i1 + nx*( (i2-1) + (i3-1)*ny/nproc_fft  ) + padatf
             fofg(1,igf) = real(data_cplx(i1,i3,i2dat)) * factor_fft
             fofg(2,igf) =aimag(data_cplx(i1,i3,i2dat)) * factor_fft
          end do
       end do
    end do
 end do

 call fftw_destroy_plan(plan_fw)
 call fftw_free(cdata_cplx)
 call fftw_free(cdata_real)

#else
 MSG_ERROR("FFTW3_MPI support not activated")
 ABI_UNUSED((/nfft,ngfft(1),ndat,comm_fft/))
 ABI_UNUSED((/fftn2_distrib(1),ffti2_local(1)/))
 ABI_UNUSED((/fftn3_distrib(1),ffti3_local(1)/))
 if (PRESENT(fftw_flags)) then
    ABI_UNUSED(fftw_flags)
 end if
 ABI_UNUSED(fofg(1,1))
 ABI_UNUSED(fofr(1))
#endif

end subroutine fftw3_mpifourdp_r2c
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/old_fftw3_mpifourdp
!! NAME
!!  old_fftw3_mpifourdp
!!
!! FUNCTION
!!
!! INPUTS
!! cplex=1 if fofr is real, 2 if fofr is complex
!! nfft=(effective) number of FFT grid points (for this processor)
!! ndat=Number of FFTs to be done.
!! isign= +1 : fofg(G) => fofr(R);
!!        -1 : fofr(R) => fofg(G)
!! fftn2_distrib(n2)=  rank of the processor which own fft planes in 2nd dimension for fourdp
!! ffti2_local(n2) = local i2 indices in fourdp
!! fftn3_distrib(n3) = rank of the processor which own fft planes in 3rd dimension for fourdp
!! ffti3_local(n3) = local i3 indices in fourdp
!! comm_fft=MPI communicator for the FFT
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator.
!!
!! SIDE EFFECTS
!! Input/Output
!! fofg(2,nfft*ndat)=f(G), complex.
!! fofr(cplex*nfft*ndat)=input function f(r) (real or complex)
!!
!! PARENTS
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine old_fftw3_mpifourdp(cplex,nfft,ngfft,ndat,isign,&
  fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft,fftw_flags)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,ndat,isign,comm_fft
 integer,optional,intent(in) :: fftw_flags
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(in) :: fftn2_distrib(ngfft(2)),ffti2_local(ngfft(2))
 integer,intent(in) :: fftn3_distrib(ngfft(3)),ffti3_local(ngfft(3))
 real(dp),intent(inout) :: fofg(2,nfft*ndat),fofr(cplex*nfft*ndat)

#ifdef HAVE_FFTW3_MPI
!Local variables-------------------------------
!scalars
 integer :: nx,ny,nz,my_flags

! *************************************************************************

 my_flags=ABI_FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags= fftw_flags

 nx=ngfft(1); ny=ngfft(2); nz=ngfft(3)
 !me_fft=ngfft(11); nproc_fft=ngfft(10)

 select case (cplex)

 case (1)

   ! Complex to Complex.
   ! This one is ok when ndat > 1
   !call fftw3_mpifourdp_c2c(cplex,nfft,ngfft,ndat,isign,&
   !& fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft,fftw_flags=my_flags)
   !return

   ! r2c or c2r case.
   ! FIXME this one is buggy when ndat > 1
   select case (isign)
   case (ABI_FFTW_FORWARD)
     ! +1; R --> G
    call fftw3_mpifourdp_r2c(nfft,ngfft,ndat,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,&
&     fofg,fofr,comm_fft,fftw_flags=my_flags)

   case (ABI_FFTW_BACKWARD)
     ! -1; G --> R
    call fftw3_mpifourdp_c2r(nfft,ngfft,ndat,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,&
&     fofg,fofr,comm_fft,fftw_flags=my_flags)

   case default
     MSG_BUG("Wrong isign")
   end select

 case (2)
   ! Complex to Complex.
   call fftw3_mpifourdp_c2c(cplex,nfft,ngfft,ndat,isign,&
&    fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft,fftw_flags=my_flags)

 case default
   MSG_BUG(" Wrong value for cplex")
 end select

#else
 MSG_ERROR("FFTW3_MPI support not activated")
 ABI_UNUSED((/cplex,nfft,ngfft(1),ndat,isign,comm_fft/))
 ABI_UNUSED((/fftn2_distrib(1),ffti2_local(1)/))
 ABI_UNUSED((/fftn3_distrib(1),ffti3_local(1)/))
 if (PRESENT(fftw_flags)) then
    ABI_UNUSED(fftw_flags)
 end if
 ABI_UNUSED(fofg(1,1))
 ABI_UNUSED(fofr(1))
#endif

end subroutine old_fftw3_mpifourdp
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_mpifourdp_c2c
!! NAME
!! fftw3_mpifourdp_c2c
!!
!! FUNCTION
!! Driver routine for many out-of-place 3D complex-to-complex FFTs of lengths n1, n2, n3.
!!
!! INPUTS
!! cplex=1 if fofr is real, 2 if fofr is complex
!! nfft=(effective) number of FFT grid points (for this processor)
!! ndat=Number of FFTs to be done.
!! isign=sign of Fourier transform exponent: current convention uses
!!   +1 for transforming from G to r,
!!   -1 for transforming from r to G.
!! fftn2_distrib(n2)=  rank of the processor which own fft planes in 2nd dimension for fourdp
!! ffti2_local(n2) = local i2 indices in fourdp
!! fftn3_distrib(n3) = rank of the processor which own fft planes in 3rd dimension for fourdp
!! ffti3_local(n3) = local i3 indices in fourdp
!! comm_fft=MPI communicator for the FFT
!! [fftw_flags]=Flags used to create the plan. They can be combined with the "+" operator.
!! fin(2*ldx*ldy*ldz*ndat)=The complex array to be transformed.
!!
!! TODO
!!   Add c2r and r2c version.
!!
!! SIDE EFFECTS
!! Input/Output
!! fofg(2,nfft*ndat)=f(G), complex.
!! fofr(cplex*nfft*ndat)=input function f(r) (real or complex)
!!
!! PARENTS
!!      m_fftw3
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_mpifourdp_c2c(cplex,nfft,ngfft,ndat,isign,&
&  fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft,fftw_flags)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,isign,nfft,ndat,comm_fft
 integer,optional,intent(in) :: fftw_flags
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(in) :: fftn2_distrib(ngfft(2)),ffti2_local(ngfft(2))
 integer,intent(in) :: fftn3_distrib(ngfft(3)),ffti3_local(ngfft(3))
 real(dp),intent(inout) :: fofg(2,nfft*ndat),fofr(cplex*nfft*ndat)

#ifdef HAVE_FFTW3_MPI
!Local variables-------------------------------
!scalars
 integer,parameter :: rank3=3
 integer :: n1,n2,n3,n4,n5,n6,nd2proc,nd3proc,my_flags,me_fft,nproc_fft
 integer(C_INTPTR_T) :: alloc_local,local_n0,local_0_start,local_n1,local_1_start
 type(C_PTR) :: plan,cptr_cdata
!arrays
 integer(C_INTPTR_T) :: fft_sizes(4)
 complex(C_DOUBLE_COMPLEX), ABI_CONTIGUOUS pointer :: f03_cdata(:)

!*************************************************************************

 my_flags=ABI_FFTW_ESTIMATE; if (PRESENT(fftw_flags)) my_flags= fftw_flags

 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 ! No augmentation as FFTW3 does not support it
 n4=n1; n5=n2; n6=n3
 me_fft=ngfft(11); nproc_fft=ngfft(10)

 nd2proc=((n2-1)/nproc_fft) +1
 nd3proc=((n6-1)/nproc_fft) +1

 ! Get local data size and allocate (note dimension reversal, we call the C interface directly!)
 fft_sizes = [n3,n2,n1,ndat]

 ! Use TRANSPOSED_OUT
 my_flags = ior(ABI_FFTW_ESTIMATE, ABI_FFTW_MPI_TRANSPOSED_OUT)

 if (isign == ABI_FFTW_BACKWARD) then
   ! G --> R, Exchange n2 and n3
   fft_sizes = [n2,n3,n1,ndat]
   !my_flags = ior(ABI_FFTW_ESTIMATE, ABI_FFTW_MPI_TRANSPOSED_IN)
 end if

 alloc_local = fftw_mpi_local_size_many_transposed(&
&      rank3,fft_sizes(1:3),fft_sizes(4), &
&      FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, comm_fft, &
&      local_n0,local_0_start, &
&      local_n1,local_1_start)

 ! C to F
 !local_0_start = local_0_start + 1
 !local_1_start = local_1_start + 1
 !write(std_out,*)"local_n0,local_0_start,alloc_local",local_n0,local_0_start,alloc_local
 !write(std_out,*)"local_n1,local_1_start,alloc_local",local_n1,local_1_start,alloc_local

 ! Allocate cptr_cdata, associate to F pointer and build the plane.
 cptr_cdata = fftw_alloc_complex(alloc_local)

 call c_f_pointer(cptr_cdata, f03_cdata, [alloc_local])

 plan = fftw_mpi_plan_many_dft(rank3,fft_sizes(1:3),fft_sizes(4), &
&                              FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, &
&                              f03_cdata,f03_cdata,comm_fft,isign,my_flags)

 select case (isign)
 case (ABI_FFTW_BACKWARD)
     ! G --> R
     ABI_CHECK(local_n0 == nd2proc, "local_n0 != nd2proc")

     call mpifft_fg2dbox_dpc(nfft,ndat,fofg,n1,n2,n3,n4,nd2proc,n6,fftn2_distrib,ffti2_local,me_fft,f03_cdata)

     ! Compute transform.
     call fftw_mpi_execute_dft(plan, f03_cdata, f03_cdata)

     call mpifft_dbox2fr_dpc(n1,n2,n3,n4,n5,nd3proc,ndat,fftn3_distrib,ffti3_local,me_fft,f03_cdata,cplex,nfft,fofr)

 case (ABI_FFTW_FORWARD)
     ! R --> G
     ABI_CHECK(local_n0 == nd3proc, "local_n0 != nd3proc")

     call mpifft_fr2dbox_dpc(cplex,nfft,ndat,fofr,n1,n2,n3,n4,n5,nd3proc,fftn3_distrib,ffti3_local,me_fft,f03_cdata)

     ! Compute transform.
     call fftw_mpi_execute_dft(plan, f03_cdata, f03_cdata)

     ! Scale results.
     call mpifft_dbox2fg_dpc(n1,n2,n3,n4,nd2proc,n6,ndat,fftn2_distrib,ffti2_local,me_fft,f03_cdata,nfft,fofg)

 case default
   MSG_ERROR("Wrong sign")
 end select

 call fftw_destroy_plan(plan)
 call fftw_free(cptr_cdata)

#else
 MSG_ERROR("FFTW3_MPI support not activated")
 ABI_UNUSED((/cplex,nfft,ngfft(1),ndat,isign,comm_fft/))
 ABI_UNUSED((/fftn2_distrib(1),ffti2_local(1)/))
 ABI_UNUSED((/fftn3_distrib(1),ffti3_local(1)/))
 if (PRESENT(fftw_flags)) then
    ABI_UNUSED(fftw_flags)
 end if
 ABI_UNUSED(fofg(1,1))
 ABI_UNUSED(fofr(1))
#endif

end subroutine fftw3_mpifourdp_c2c
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_mpiback_wf
!! NAME
!!  fftw3_mpiback_wf
!!
!! FUNCTION
!!   Does multiple 3-dim backward FFTs from Fourier into real space
!!   Adopt standard convention that isign=1 for backward transform
!!
!!   CALCULATES THE DISCRETE FOURIER TRANSFORM ZF(I1,I2,I3)=
!!
!!   S_(j1,j2,j3) EXP(isign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) ZF(j1,j3,j2)
!!
!!   in parallel using MPI/OpenMP.
!!
!! INPUTS:
!!    cplexwf=1 if wavefunction is real, 2 if complex
!!    ndat=Number of wavefunctions to transform.
!!    n1,n2,n3: logical dimension of the transform. As transform lengths
!!              most products of the prime factors 2,3,5 are allowed.
!!              The detailed table with allowed transform lengths can be found in subroutine CTRIG
!!    nd1,nd2,nd3: Leading Dimension of ZR
!!    nd3proc=((nd3-1)/nproc_fft)+1 maximal number of big box 3rd dim slices for one proc
!!    max1 is positive or zero; m1 >=max1+1
!!      i1= 1... max1+1 corresponds to positive and zero wavevectors 0 ... max1
!!      then, if m1 > max1+1, one has min1=max1-m1+1 and
!!      i1= max1+2 ... m1 corresponds to negative wavevectors min1 ... -1
!!    max2 and max3 have a similar definition of range
!!    m1,m2,m3=Size of the box enclosing the G-sphere.
!!    md1,md2,md3: Dimension of ZF given on the **small** FFT box.
!!    md2proc=((md2-1)/nproc_fft)+1 maximal number of small box 2nd dim slices for one proc
!!    nproc_fft: number of processors used as returned by MPI_COMM_SIZE
!!    comm_fft=MPI communicator for the FFT.
!!    ZF: input array (note the switch of i2 and i3)
!!          real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!!          imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!!
!! OUTPUTS
!!    ZR: output array
!!          ZR(1,i1,i2,i3,idat)=real(R(i1,i2,i3,idat))
!!          ZR(2,i1,i2,i3,idat)=imag(R(i1,i2,i3,idat))
!!        i1=1,n1 , i2=1,n2 , i3=1,n3 , idat=1,ndat
!!
!! NOTES
!!   The maximum number of processors that can reasonably be used is max(n2/2,n3/2)
!!
!!   It is very important to find the optimal
!!   value of NCACHE. NCACHE determines the size of the work array ZW, that
!!   has to fit into cache. It has therefore to be chosen to equal roughly
!!   half the size of the physical cache in units of real*8 numbers.
!!   The optimal value of ncache can easily be determined by numerical
!!   experimentation. A too large value of ncache leads to a dramatic
!!   and sudden decrease of performance, a too small value to a to a
!!   slow and less dramatic decrease of performance. If NCACHE is set
!!   to a value so small, that not even a single one dimensional transform
!!   can be done in the workarray zw, the program stops with an error message.
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_mpiback_wf(cplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc,&
&  max1,max2,max3,m1,m2,m3,md1,md2proc,md3,zf,zr,comm_fft)

!Arguments ------------------------------------
 integer,intent(in) :: cplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc
 integer,intent(in) :: max1,max2,max3,m1,m2,m3,md1,md2proc,md3,comm_fft
 real(dp),intent(in) :: zf(2,md1,md3,md2proc,ndat)
 real(dp),intent(out) :: zr(2,nd1,nd2,nd3proc,ndat)

#ifdef HAVE_FFTW3
!Local variables-------------------------------
 integer,parameter :: nt1=1
 integer :: j,i1,i2,idat,ierr,includelast
 integer :: ioption,j2,j3,j2st,jp2st,jeff,lzt,m1zt,ma,mb,n1dfft,nnd3
 integer :: lot1,lot2,lot3
 integer :: m2eff,ncache,n1eff,n1half,nproc_fft,me_fft,nthreads
 integer(KIND_FFTW_PLAN) :: bw_plan1_lot,bw_plan1_rest
 integer(KIND_FFTW_PLAN) :: bw_plan2_lot,bw_plan2_rest
 integer(KIND_FFTW_PLAN) :: bw_plan3_lot,bw_plan3_rest
 !type(C_PTR) :: zw_cptr,zt_cptr
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: zmpi1(:,:,:,:),zmpi2(:,:,:,:)  ! work arrays for MPI
 real(dp),allocatable :: zw(:,:),zt(:,:,:) ! cache work array and array for transpositions
 !real(dp),ABI_CONTIGUOUS pointer :: zw(:,:),zt(:,:,:)
! FFT work arrays
 real(dp) :: tsec(2)

! *************************************************************************

 !call wrtout(std_out,"mpiback standard ALLTOALL + FFTW3","COLL")

 ! FIXME must provide a default value but which one?
 ! ioption = 0
 ioption = 1
 !if (paral_kgb==1) ioption=1

 nproc_fft = xmpi_comm_size(comm_fft); me_fft = xmpi_comm_rank(comm_fft)

 ! Find cache size that gives optimal performance on machine
 ncache=2*max(n1,n2,n3,1024)
 if (ncache/(2*max(n1,n2,n3))<1) then
   write(msg,"(5a)") &
&    'ncache has to be enlarged to be able to hold at',ch10, &
&    'least one 1-d FFT of each size even though this will',ch10,&
&    'reduce the performance for shorter transform lengths'
    MSG_ERROR(msg)
 end if

 ! Effective m1 and m2 (complex-to-complex or real-to-complex)
 n1eff=n1; m2eff=m2; m1zt=n1
 if (cplexwf==1) then
   n1eff=(n1+1)/2; m2eff=m2/2+1; m1zt=2*(n1/2+1)
 end if

 lzt=m2eff
 if (mod(m2eff,2)==0) lzt=lzt+1
 if (mod(m2eff,4)==0) lzt=lzt+1

 ! maximal number of big box 3rd dim slices for all procs
 nnd3=nd3proc*nproc_fft

 ! Allocate cache work array and work arrays for MPI transpositions.
 ABI_MALLOC(zw,(2,ncache/2))
 ABI_MALLOC(zt,(2,lzt,m1zt))

 !call fftw3_alloc_real([2,ncache/2],zw_cptr,zw)
 !call fftw3_alloc_real([2,lzt,m1zt],zt_cptr,zt)

 ABI_MALLOC(zmpi2,(2,md1,md2proc,nnd3))
 if (nproc_fft>1)  then
   ABI_MALLOC(zmpi1,(2,md1,md2proc,nnd3))
 end if

!DEBUG
! write(std_out,'(2a,3i4)' )itoa(me_fft),': fftw3_mpiback_wf,zf n1,n2,n3',n1,n2,n3
! write(std_out,'(2a,3i4)' )itoa(me_fft),': nd1,nd2,nd3proc',nd1,nd2,nd3proc
! write(std_out,'(2a,3i4)' )itoa(me_fft),': m1,m2,m3',m1,m2,m3
! write(std_out,'(2a,3i4)' )itoa(me_fft),': max1,max2,max3',max1,max2,max3
! write(std_out,'(2a,3i4)' )itoa(me_fft),': md1,md2proc,md3',md1,md2proc,md3
! write(std_out,'(2a,3i4)' )itoa(me_fft),'n1eff,m2eff,m1zt',n1eff,m2eff,m1zt
!ENDDEBUG

 ! Create plans.
 ! The prototype for sfftw_plan_many_dft is:
 ! sfftw_plan_many_dft(rank, n, howmany,
 !   fin,  iembed, istride, idist,
 !   fout, oembed, ostride, odist, isign, my_flags)

 lot3=ncache/(2*n3)
 lot1=ncache/(2*n1)
 lot2=ncache/(2*n2)

 nthreads = xomp_get_num_threads(open_parallel=.TRUE.)
 !nthreads = 1

 bw_plan3_lot = dplan_many_dft_2D(1, [n3], lot3, &
&    zw, [ncache/2], lot3, 1,                          &
&    zw, [ncache/2], lot3, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(m1, lot3) /= 0) then
   bw_plan3_rest = dplan_many_dft_2D(1, [n3], mod(m1, lot3), &
&      zw, [ncache/2], lot3, 1,                                    &
&      zw, [ncache/2], lot3, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 bw_plan1_lot = dplan_many_dft_2D(1, [n1], lot1, &
&    zw, [ncache/2],  lot1, 1,                         &
&    zt, [lzt, m1zt], lzt, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(m2eff, lot1) /= 0) then
   bw_plan1_rest = dplan_many_dft_2D(1, [n1], mod(m2eff, lot1), &
&      zw, [ncache/2],  lot1, 1,                                      &
&      zt, [lzt, m1zt], lzt, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 bw_plan2_lot = dplan_many_dft_2D(1, [n2], lot2, &
&    zw, [ncache/2], lot2, 1,                          &
&    zr, [nd1,nd2,nd3proc,ndat], nd1, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(n1eff, lot2) /= 0) then
   bw_plan2_rest = dplan_many_dft_2D(1, [n2], mod(n1eff,lot2), &
&      zw, [ncache/2], lot2, 1,                                      &
&      zr, [nd1,nd2,nd3proc,ndat], nd1, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 do idat=1,ndat
    ! transform along z axis
    ! input: G1,G3,G2,(Gp2)

    ! Loop over the y planes treated by this node and trasform n1ddft G_z lines.
    do j2=1,md2proc
      ! if (me_fft*md2proc+j2<=m2eff) then !a faire plus tard
      do i1=1,m1,lot3
        ma=i1
        mb=min(i1+(lot3-1),m1)
        n1dfft=mb-ma+1

        ! zero-pad n1dfft G_z lines
        ! input:  G1,G3,G2,(Gp2)
        ! output: G1,R3,G2,(Gp2)
        call fill_cent(md1,md3,lot3,n1dfft,max3,m3,n3,zf(1,i1,1,j2,idat),zw)

        ! Transform along z.
        if (n1dfft == lot3) then
          call dfftw_execute_dft(bw_plan3_lot, zw, zw)
        else
          call dfftw_execute_dft(bw_plan3_rest, zw, zw)
        end if

        ! Local rotation.
        ! input:  G1,R3,G2,(Gp2)
        ! output: G1,G2,R3,(Gp2)
        call scramble(i1,j2,lot3,n1dfft,md1,n3,md2proc,nnd3,zw,zmpi2)
      end do
    end do ! j2

    ! Interprocessor data transposition
    ! input:  G1,G2,R3,Rp3,(Gp2)
    ! output: G1,G2,R3,Gp2,(Rp3)
    if (nproc_fft>1) then
      call timab(543,1,tsec)
      call xmpi_alltoall(zmpi2,2*md1*md2proc*nd3proc, &
&                        zmpi1,2*md1*md2proc*nd3proc,comm_fft,ierr)
      call timab(543,2,tsec)
    end if

    ! Loop over the z treated by this node.
    do j3=1,nd3proc
      if (me_fft*nd3proc+j3 <= n3) then
        Jp2st=1; J2st=1

        ! Loop over G_y in the small box.
        do j=1,m2eff,lot1
          ma=j
          mb=min(j+(lot1-1),m2eff)
          n1dfft=mb-ma+1

          ! Zero-pad input.
          ! input:  G1,G2,R3,JG2,(Rp3)
          ! output: G2,G1,R3,JG2,(Rp3)
          if (nproc_fft==1) then
            call mpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot1,max1,md1,m1,n1,&
&             md2proc,nd3proc,nproc_fft,ioption,zmpi2,zw,max2,m2,n2)
          else
            call mpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot1,max1,md1,m1,n1,&
&             md2proc,nd3proc,nproc_fft,ioption,zmpi1,zw,max2,m2,n2)
          end if

          ! Transform along x
          ! input:  G2,G1,R3,(Rp3)
          ! output: G2,R1,R3,(Rp3)
          if (n1dfft == lot1) then
            call dfftw_execute_dft(bw_plan1_lot, zw, zt(1,j,1))
          else
            call dfftw_execute_dft(bw_plan1_rest, zw, zt(1,j,1))
          end if

        end do ! j

        ! Transform along y axis (take into account c2c or c2r case).
        ! Must loop over the full box.
        do j=1,n1eff,lot2
          ma=j
          mb=min(j+(lot2-1),n1eff)
          n1dfft=mb-ma+1
          includelast=1

          if (cplexwf==1) then
            jeff=2*j-1
            if (mb==n1eff .and. n1eff*2/=n1) includelast=0
          end if

          ! Zero-pad the input.
          ! input:  G2,R1,R3,(Rp3)
          ! output: R1,G2,R3,(Rp3)
          if (cplexwf==2) then
            call switch_cent(n1dfft,max2,m2,n2,lot2,n1,lzt,zt(1,1,j),zw)
          else
            call switchreal_cent(includelast,n1dfft,max2,n2,lot2,m1zt,lzt,zt(1,1,jeff),zw)
          end if

          ! input:  R1,G2,R3,(Rp3)
          ! output: R1,R2,R3,(Rp3)
          if (n1dfft == lot2) then
            call dfftw_execute_dft(bw_plan2_lot, zw, zr(1,j,1,j3,idat))
          else
            call dfftw_execute_dft(bw_plan2_rest, zw, zr(1,j,1,j3,idat))
          end if

        end do

        ! Treat real wavefunctions.
        if (cplexwf==1) then
          n1half=n1/2
          ! If odd
          if (n1half*2/=n1) then
            do i2=1,n2
              zr(1,n1,i2,j3,idat)=zr(1,n1eff,i2,j3,idat)
              zr(2,n1,i2,j3,idat)=zero
            end do
          end if
          do i2=1,n2
            do i1=n1half,1,-1
              zr(1,2*i1-1,i2,j3,idat)=zr(1,i1,i2,j3,idat)
              zr(1,2*i1  ,i2,j3,idat)=zr(2,i1,i2,j3,idat)
              zr(2,2*i1-1,i2,j3,idat)=zero
              zr(2,2*i1  ,i2,j3,idat)=zero
            end do
          end do
        end if

      end if
   end do ! j3
 end do ! idat

 call dfftw_destroy_plan(bw_plan3_lot)
 if (mod(m1, lot3) /= 0) then
   call dfftw_destroy_plan(bw_plan3_rest)
 end if

 call dfftw_destroy_plan(bw_plan1_lot)
 if (mod(m2eff, lot1) /= 0) then
   call dfftw_destroy_plan(bw_plan1_rest)
 end if

 call dfftw_destroy_plan(bw_plan2_lot)
 if (mod(n1eff, lot2) /= 0) then
   call dfftw_destroy_plan(bw_plan2_rest)
 end if

 ABI_FREE(zmpi2)
 ABI_FREE(zw)
 ABI_FREE(zt)
 if (nproc_fft>1)  then
   ABI_FREE(zmpi1)
 end if

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/cplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc/))
 ABI_UNUSED((/ max1,max2,max3,m1,m2,m3,md1,md2proc,md3,comm_fft/))
 ABI_UNUSED((/zf(1,1,1,1,1),zr(1,1,1,1,1)/))
#endif

end subroutine fftw3_mpiback_wf
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_mpiforw_wf
!! NAME
!!  fftw3_mpiforw_wf
!!
!! FUNCTION
!!   Does multiple 3-dim backward FFTs from real into Fourier space
!!   Adopt standard convention that isign=-1 for forward transform
!!   CALCULATES THE DISCRETE FOURIERTRANSFORM
!!
!!   ZF(I1,I3,I2)=S_(j1,j2,j3) EXP(isign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) ZR(j1,j2,j3)
!!
!!   in parallel using MPI/OpenMP.
!!
!! INPUT:
!!   ZR: input array
!!        ZR(1,i1,i2,i3,idat)=real(R(i1,i2,i3,idat))
!!        ZR(2,i1,i2,i3,idat)=imag(R(i1,i2,i3,idat))
!!        i1=1,n1 , i2=1,n2 , i3=1,n3 , idat=1,ndat
!!   NOTE that ZR is changed by the routine
!!
!!   n1,n2,n3: logical dimension of the transform. As transform lengths
!!             most products of the prime factors 2,3,5 are allowed.
!!             The detailed table with allowed transform lengths can
!!             be found in subroutine CTRIG
!!   nd1,nd2,nd3: Dimension of ZR
!!   nd3proc=((nd3-1)/nproc_fft)+1  maximal number of big box 3rd dim slices for one proc
!!
!! OUTPUT:
!!   ZF: output array (note the switch of i2 and i3)
!!        real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!!        imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!!   max1 is positive or zero ; m1 >=max1+1
!!     i1= 1... max1+1 corresponds to positive and zero wavevectors 0 ... max1
!!     then, if m1 > max1+1, one has min1=max1-m1+1 and
!!     i1= max1+2 ... m1 corresponds to negative wavevectors min1 ... -1
!!     i2 and i3 have a similar definition of range
!!   idat=1,ndat
!!   md1,md2,md3: Dimension of ZF
!!   md2proc=((md2-1)/nproc_fft)+1  maximal number of small box 2nd dim slices for one proc
!!   nproc_fft: number of processors used as returned by MPI_COMM_SIZE
!!   me_fft: [0:nproc-1] rank of the processor in the FFT communicator.
!!   comm_fft=MPI communicator for parallel FFT.
!!
!! NOTES
!!  The maximum number of processors that can reasonably be used is max(n2/2,n3/2)
!!
!!  It is very important to find the optimal
!!  value of NCACHE. NCACHE determines the size of the work array ZW, that
!!  has to fit into cache. It has therefore to be chosen to equal roughly
!!   half the size of the physical cache in units of real*8 numbers.
!!  The optimal value of ncache can easily be determined by numerical
!!  experimentation. A too large value of ncache leads to a dramatic
!!  and sudden decrease of performance, a too small value to a to a
!!  slow and less dramatic decrease of performance. If NCACHE is set
!!  to a value so small, that not even a single one dimensional transform
!!  can be done in the workarray zw, the program stops with an error message.
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_mpiforw_wf(cplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc,&
&        max1,max2,max3,m1,m2,m3,md1,md2proc,md3,zr,zf,comm_fft)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc
 integer,intent(in) :: max1,max2,max3,m1,m2,m3,md1,md2proc,md3,comm_fft
!arrays
 real(dp),intent(inout) :: zr(2,nd1,nd2,nd3proc,ndat)
 real(dp),intent(out) :: zf(2,md1,md3,md2proc,ndat)

!Local variables-------------------------------
!scalars
#ifdef HAVE_FFTW3
 integer :: j,i1,i2,i3,idat,ierr,nproc_fft,me_fft,nthreads
 integer :: ioption,j2,j3,j2st,jp2st,lot1,lot2,lot3,lzt,m1zt,ma,mb,n1dfft,nnd3
 integer :: m2eff,ncache,n1eff,n1half,i1inv,i2inv,i3inv
 integer(KIND_FFTW_PLAN) :: fw_plan1_lot,fw_plan1_rest
 integer(KIND_FFTW_PLAN) :: fw_plan2_lot,fw_plan2_rest
 integer(KIND_FFTW_PLAN) :: fw_plan3_lot,fw_plan3_rest
 character(len=500) :: msg
!arrays
 real(dp), allocatable :: zmpi1(:,:,:,:),zmpi2(:,:,:,:) ! work arrays for MPI
 real(dp),allocatable :: zw(:,:),zt(:,:,:) ! cache work array and array for transpositions
! FFT work arrays
 real(dp) :: tsec(2)

! *************************************************************************

 ! FIXME must provide a default value but which one?
 !ioption = 0
 ioption = 1
 !if (paral_kgb==1) ioption=1

 nproc_fft = xmpi_comm_size(comm_fft); me_fft = xmpi_comm_rank(comm_fft)

 ! find cache size that gives optimal performance on machine
 ncache=2*max(n1,n2,n3,1024)
 !ncache=2*max(n1,n2,n3,16*1024)

 if (ncache/(2*max(n1,n2,n3))<1) then
   write(msg,'(5a)') &
&    'ncache has to be enlarged to be able to hold at',ch10, &
&    'least one 1-d FFT of each size even though this will',ch10,&
&    'reduce the performance for shorter transform lengths'
   MSG_ERROR(msg)
 end if

 ! Effective m1 and m2 (complex-to-complex or real-to-complex)
 n1eff=n1; m2eff=m2; m1zt=n1
 if (cplexwf==1) then
   n1eff=(n1+1)/2; m2eff=m2/2+1; m1zt=2*(n1/2+1)
 end if

 lzt=m2eff
 if (mod(m2eff,2)==0) lzt=lzt+1
 if (mod(m2eff,4)==0) lzt=lzt+1

 ! maximal number of big box 3rd dim slices for all procs
 nnd3=nd3proc*nproc_fft

 ABI_MALLOC(zw,(2,ncache/2))
 ABI_MALLOC(zt,(2,lzt,m1zt))
 ABI_MALLOC(zmpi2,(2,md1,md2proc,nnd3))
 if (nproc_fft>1)  then
   ABI_MALLOC(zmpi1,(2,md1,md2proc,nnd3))
 end if

!DEBUG
! write(std_out,'(2a,3i4)' )itoa(me_fft),'fftw3_mpiforw_wf, enter, i1,i2,i3,zr,n1,n2,n3',n1,n2,n3
! write(std_out,'(2a,3i4)' )itoa(me_fft),'nd1,nd2,nd3proc',nd1,nd2,nd3proc
! write(std_out,'(2a,3i4)' )itoa(me_fft),'m1,m2,m3',m1,m2,m3
! write(std_out,'(2a,3i4)' )itoa(me_fft),'max1,max2,max3',max1,max2,max3
! write(std_out,'(2a,3i4)' )itoa(me_fft),'md1,md2proc,md3',md1,md2proc,md3
! write(std_out,'(2a,3i4)' )itoa(me_fft),'n1eff,m2eff,m1zt',n1eff,m2eff,m1zt
!ENDDEBUG

 ! Create plans.
 ! The prototype for sfftw_plan_many_dft is:
 ! sfftw_plan_many_dft(rank, n, howmany,
 !   fin,  iembed, istride, idist,
 !   fout, oembed, ostride, odist, isign, my_flags)

 lot2=ncache/(2*n2)
 lot1=ncache/(2*n1)
 lot3=ncache/(2*n3)

 nthreads = xomp_get_num_threads(open_parallel=.TRUE.)
 !nthreads = 1

 fw_plan3_lot = dplan_many_dft_2D(1, [n3], lot3, &
&    zw, [ncache/2], lot3, 1,                          &
&    zw, [ncache/2], lot3, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE,nthreads)

 if (mod(m1, lot3) /= 0) then
   fw_plan3_rest = dplan_many_dft_2D(1, [n3], mod(m1, lot3), &
&    zw, [ncache/2], lot3, 1,                                      &
&    zw, [ncache/2], lot3, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 fw_plan1_lot = dplan_many_dft_2D(1, [n1], lot1, &
&    zt, [lzt, m1zt],   lzt,  1,                       &
&    zw, [ncache/2], lot1, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(m2eff, lot1) /= 0) then
   fw_plan1_rest = dplan_many_dft_2D(1, [n1], mod(m2eff, lot1), &
&    zt, [lzt, m1zt],   lzt, 1,                                       &
&    zw, [ncache/2], lot1, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 fw_plan2_lot = dplan_many_dft_2D(1, [n2], lot2, &
&    zr, [nd1,nd2,nd3proc,ndat], nd1, 1,               &
&    zw, [ncache/2], lot2, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(n1eff, lot2) /= 0) then
   fw_plan2_rest = dplan_many_dft_2D(1, [n2], mod(n1eff,lot2), &
&    zr, [nd1,nd2,nd3proc,ndat], nd1, 1,                             &
&    zw, [ncache/2], lot2, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 do idat=1,ndat
   ! Loop over the z-planes treated by this node
   do j3=1,nd3proc

     if (me_fft*nd3proc+j3 <= n3) then
       Jp2st=1
       J2st=1

       ! Treat real wavefunctions.
       if (cplexwf==1) then
         n1half=n1/2
         do i2=1,n2
           do i1=1,n1half
             zr(1,i1,i2,j3,idat)=zr(1,2*i1-1,i2,j3,idat)
             zr(2,i1,i2,j3,idat)=zr(1,2*i1  ,i2,j3,idat)
           end do
         end do
         ! If odd
         if(n1half*2/=n1)then
           do i2=1,n2
             zr(1,n1eff,i2,j3,idat)=zr(1,n1,i2,j3,idat)
             zr(2,n1eff,i2,j3,idat)=zero
           end do
         end if
       end if

       ! transform along y axis
       ! input: R1,R2,R3,(Rp3)
       ! input: R1,G2,R3,(Rp3)
       do j=1,n1eff,lot2
         ma=j
         mb=min(j+(lot2-1),n1eff)
         n1dfft=mb-ma+1

         if (n1dfft == lot2) then
           call dfftw_execute_dft(fw_plan2_lot,  zr(1,j,1,j3,idat), zw)
         else
           call dfftw_execute_dft(fw_plan2_rest, zr(1,j,1,j3,idat), zw)
         end if

         ! input:  R1,G2,R3,(Rp3)
         ! output: G2,R1,R3,(Rp3)
         if (cplexwf==2) then
           call unswitch_cent(n1dfft,max2,m2,n2,lot2,n1,lzt,zw,zt(1,1,j))
         else
           call unswitchreal_cent(n1dfft,max2,n2,lot2,n1,lzt,zw,zt(1,1,2*j-1))
         end if
       end do

       ! transform along x axis
       ! input: G2,R1,R3,(Rp3)
       do j=1,m2eff,lot1
         ma=j
         mb=min(j+(lot1-1),m2eff)
         n1dfft=mb-ma+1

         if (n1dfft == lot1) then
           call dfftw_execute_dft(fw_plan1_lot,  zt(1,j,1), zw)
         else
           call dfftw_execute_dft(fw_plan1_rest, zt(1,j,1), zw)
         end if
         ! output: G2,G1,R3,(Rp3)

         ! input:  G2,G1,R3,Gp2,(Rp3)
         ! output: G1,G2,R3,Gp2,(Rp3)
         if (nproc_fft==1) then
           call unmpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot1,max1,md1,m1,n1,&
&            md2proc,nd3proc,nproc_fft,ioption,zw,zmpi2)
         else
           call unmpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot1,max1,md1,m1,n1,&
&            md2proc,nd3proc,nproc_fft,ioption,zw,zmpi1)
         end if
       end do

      end if
    end do ! j3

    ! Interprocessor data transposition
    ! input:  G1,G2,R3,Gp2,(Rp3)
    ! output: G1,G2,R3,Rp3,(Gp2)
    if (nproc_fft>1) then
      call timab(544,1,tsec)
      call xmpi_alltoall(zmpi1,2*md1*md2proc*nd3proc, &
&                        zmpi2,2*md1*md2proc*nd3proc,comm_fft,ierr)
      call timab(544,2,tsec)
    end if

    ! transform along z axis
    ! input: G1,G2,R3,(Gp2)

    do j2=1,md2proc
      if (me_fft*md2proc+j2 <= m2eff) then
        ! write(std_out,*)' forwf_wf : before unscramble, j2,md2proc,me_fft,m2=',j2,md2proc,me_fft,m2
        do i1=1,m1,lot3
          ma=i1
          mb=min(i1+(lot3-1),m1)
          n1dfft=mb-ma+1

          ! input:  G1,G2,R3,(Gp2)
          ! output: G1,R3,G2,(Gp2)
          call unscramble(i1,j2,lot3,n1dfft,md1,n3,md2proc,nnd3,zmpi2,zw)

          if (n1dfft == lot3) then
            call dfftw_execute_dft(fw_plan3_lot, zw, zw)
          else
            call dfftw_execute_dft(fw_plan3_rest, zw, zw)
          end if

          call unfill_cent(md1,md3,lot3,n1dfft,max3,m3,n3,zw,zf(1,i1,1,j2,idat))
          ! output: G1,G3,G2,(Gp2)
        end do
      end if
    end do

    if (cplexwf==1) then
      ! Complete missing values with complex conjugate
      ! Inverse of ix is located at nx+2-ix , except for ix=1, for which it is 1.
      do i3=1,m3
        i3inv=m3+2-i3
        if(i3==1)i3inv=1

        if (m2eff>1) then
          do i2=2,m2eff
            i2inv=m2+2-i2
            zf(1,1,i3inv,i2inv,idat)= zf(1,1,i3,i2,idat)
            zf(2,1,i3inv,i2inv,idat)=-zf(2,1,i3,i2,idat)
            do i1=2,m1
              i1inv=m1+2-i1
              zf(1,i1inv,i3inv,i2inv,idat)= zf(1,i1,i3,i2,idat)
              zf(2,i1inv,i3inv,i2inv,idat)=-zf(2,i1,i3,i2,idat)
            end do
          end do
        end if
      end do
    end if

 end do ! idat

 call dfftw_destroy_plan(fw_plan3_lot)
 if (mod(m1, lot3) /= 0) then
   call dfftw_destroy_plan(fw_plan3_rest)
 end if

 call dfftw_destroy_plan(fw_plan1_lot)
 if (mod(m2eff, lot1) /= 0) then
   call dfftw_destroy_plan(fw_plan1_rest)
 end if

 call dfftw_destroy_plan(fw_plan2_lot)
 if (mod(n1eff, lot2) /= 0) then
   call dfftw_destroy_plan(fw_plan2_rest)
 end if

 ABI_FREE(zmpi2)
 ABI_FREE(zw)
 ABI_FREE(zt)
 if (nproc_fft>1)  then
   ABI_FREE(zmpi1)
 end if

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/cplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc/))
 ABI_UNUSED((/max1,max2,max3,m1,m2,m3,md1,md2proc,md3,comm_fft/))
 ABI_UNUSED((/zf(1,1,1,1,1),zr(1,1,1,1,1)/))
#endif

end subroutine fftw3_mpiforw_wf
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_mpiback
!! NAME
!!  fftw3_mpiback
!!
!! FUNCTION
!!   CALCULATES THE DISCRETE FOURIER TRANSFORM  in parallel using MPI/OpenMP
!!
!!   ZR(I1,I2,I3)= \sum_(j1,j2,j3) EXP(isign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) ZF(j1,j3,j2)
!!
!! Adopt standard convention that isign=1 for backward transform
!!
!! INPUTS:
!!    option= 1 if call from fourwf, 2 if call from other routine
!!    cplex=1 for real --> complex, 2 for complex --> complex
!!    ZF: input array in G-space (note the switch of i2 and i3)
!!
!!         real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!!         imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!!
!!         i1=1,n1 , i2=1,n2 , i3=1,n3 , idat=1,ndat
!! OUTPUTS:
!!    ZR: output array in R space.
!!
!!         ZR(1,i1,i2,i3,idat)=real(R(i1,i2,i3,idat))
!!         ZR(2,i1,i2,i3,idat)=imag(R(i1,i2,i3,idat))
!!
!!         i1=1,n1 , i2=1,n2 , i3=1,n3 , idat=1,ndat
!!
!!    nproc_fft: number of processors used as returned by MPI_COMM_SIZE
!!    me_fft: [0:nproc_fft-1] number of processor as returned by MPI_COMM_RANK
!!    n1,n2,n3: logical dimension of the transform. As transform lengths
!!              most products of the prime factors 2,3,5 are allowed.
!!              The detailed table with allowed transform lengths can
!!              be found in subroutine CTRIG
!!    nd1,nd2,nd3: Dimension of ZF and ZR
!!    nd2proc=((nd2-1)/nproc_fft)+1 maximal number of 2nd dim slices
!!    nd3proc=((nd3-1)/nproc_fft)+1 maximal number of 3rd dim slices
!!
!! NOTES:
!!   The maximum number of processors that can reasonably be used is max(n2,n3)
!!   It is very important to find the optimal
!!   value of NCACHE. NCACHE determines the size of the work array ZW, that
!!   has to fit into cache. It has therefore to be chosen to equal roughly
!!    half the size of the physical cache in units of real*8 numbers.
!!   The optimal value of ncache can easily be determined by numerical
!!   experimentation. A too large value of ncache leads to a dramatic
!!   and sudden decrease of performance, a too small value to a to a
!!   slow and less dramatic decrease of performance. If NCACHE is set
!!   to a value so small, that not even a single one dimensional transform
!!   can be done in the workarray zw, the program stops with an error message.
!!
!! PARENTS
!!      m_fftw3
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_mpiback(cplex,ndat,n1,n2,n3,nd1,nd2,nd3,nd1eff,nd2proc,nd3proc,option,zf,zr,comm_fft)

!Arguments ------------------------------------
! real space input
 integer,intent(in) :: cplex,ndat,n1,n2,n3,nd1,nd2,nd3,nd1eff,nd2proc,nd3proc,option,comm_fft
 real(dp),intent(in) :: zf(2,nd1,nd3,nd2proc,ndat)
 real(dp),intent(out) :: zr(2,nd1eff,nd2,nd3proc,ndat)

!Local variables-------------------------------
!scalaras
#ifdef HAVE_FFTW3
 integer :: j,i1,idat,ierr,includelast,j2,j2st,j3,jeff,jp2st,lzt,nthreads
 integer :: ma,mb,n1dfft,n1eff,n2eff,n1zt,ncache,nnd3,nproc_fft,me_fft,lot1,lot2,lot3
 integer(KIND_FFTW_PLAN) :: bw_plan1_lot,bw_plan1_rest
 integer(KIND_FFTW_PLAN) :: bw_plan2_lot,bw_plan2_rest
 integer(KIND_FFTW_PLAN) :: bw_plan3_lot,bw_plan3_rest
 character(len=500) :: msg
!arrays
 real(dp), allocatable :: zmpi1(:,:,:,:),zmpi2(:,:,:,:) ! work arrays for MPI
 real(dp),allocatable :: zw(:,:),zt(:,:,:) ! cache work array and array for transpositions

! *************************************************************************

 nproc_fft = xmpi_comm_size(comm_fft); me_fft = xmpi_comm_rank(comm_fft)

 ! find cache size that gives optimal performance on machine
 ncache=2*max(n1,n2,n3,1024)

 if (ncache/(2*max(n1,n2,n3))<1) then
   write(msg,'(5a)') &
&    'ncache has to be enlarged to be able to hold at',ch10, &
&    'least one 1-d FFT of each size even though this will',ch10,&
&    'reduce the performance for shorter transform lengths'
   MSG_ERROR(msg)
 end if

! check input
 if (nd1<n1 .or. nd2<n2 .or. nd3<n3) then
   MSG_ERROR("nd1<n1 .or. nd2<n2 .or. nd3<n3")
 end if

 ! Effective n1 and n2 (complex-to-complex or real-to-complex)
 n1eff=n1; n2eff=n2; n1zt=n1
 if (cplex==1) then
   n1eff=(n1+1)/2; n2eff=n2/2+1 ; n1zt=2*(n1/2+1)
 end if

 lzt=n2eff
 if (mod(n2eff,2) == 0) lzt=lzt+1
 if (mod(n2eff,4) == 0) lzt=lzt+1

! maximal number of big box 3rd dim slices for all procs
 nnd3=nd3proc*nproc_fft

 ABI_MALLOC(zw,(2,ncache/2))
 ABI_MALLOC(zt,(2,lzt,n1zt))
 ABI_MALLOC(zmpi2,(2,n1,nd2proc,nnd3))
 if (nproc_fft>1)  then
   ABI_MALLOC(zmpi1,(2,n1,nd2proc,nnd3))
 end if

!DEBUG
! write(std_out,'(a,3i4)' )'back,zf n1,n2,n3',n1,n2,n3
! write(std_out,'(a,3i4)' )'nd1,nd2,nd3proc',nd1,nd2,nd3proc
! write(std_out,'(a,3i4)' )'m1,m2,m3',m1,m2,m3
! write(std_out,'(a,3i4)' )'max1,max2,max3',max1,max2,max3
! write(std_out,'(a,3i4)' )'md1,md2proc,md3',md1,md2proc,md3
! write(std_out,'(a,3i4)' )'n1eff,m2eff,m1zt',n1eff,m2eff,m1zt
!ENDDEBUG

 ! Create plans.
 ! The prototype for sfftw_plan_many_dft is:
 ! sfftw_plan_many_dft(rank, n, howmany,
 !   fin,  iembed, istride, idist,
 !   fout, oembed, ostride, odist, isign, my_flags)

 lot3=ncache/(2*n3)
 lot1=ncache/(2*n1)
 lot2=ncache/(2*n2)

 nthreads = xomp_get_num_threads(open_parallel=.TRUE.)
 !nthreads = 1

 bw_plan3_lot = dplan_many_dft_2D(1, [n3], lot3, &
&    zw, [ncache/2], lot3, 1,                          &
&    zw, [ncache/2], lot3, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(n1, lot3) /= 0) then
   bw_plan3_rest = dplan_many_dft_2D(1, [n3], mod(n1, lot3), &
&      zw, [ncache/2], lot3, 1,                                    &
&      zw, [ncache/2], lot3, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 bw_plan1_lot = dplan_many_dft_2D(1, [n1], lot1, &
&    zw, [ncache/2],  lot1, 1,                         &
&    zt, [lzt, n1zt], lzt, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(n2eff, lot1) /= 0) then
   bw_plan1_rest = dplan_many_dft_2D(1, [n1], mod(n2eff, lot1), &
&      zw, [ncache/2], lot1, 1,                                       &
&      zt, [lzt, n1zt],   lzt, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 bw_plan2_lot = dplan_many_dft_2D(1, [n2], lot2, &
&    zw, [ncache/2], lot2, 1,                          &
&    zr, [nd1eff,nd2,nd3proc,ndat], nd1eff, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(n1eff, lot2) /= 0) then
   bw_plan2_rest = dplan_many_dft_2D(1, [n2], mod(n1eff,lot2), &
&      zw, [ncache/2], lot2, 1,                                      &
&      zr, [nd1eff,nd2,nd3proc,ndat], nd1eff, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 do idat=1,ndat
   ! transform along z axis
   ! input: I1,I3,J2,(Jp2)

   do j2=1,nd2proc
     if (me_fft*nd2proc+j2 <= n2eff) then

       do i1=1,n1,lot3
         ma=i1
         mb=min(i1+(lot3-1),n1)
         n1dfft=mb-ma+1

         ! input:  G1,G3,G2,(Gp2)
         ! output: G1,R3,G2,(Gp2)
         call fill(nd1,nd3,lot3,n1dfft,n3,zf(1,i1,1,j2,idat),zw)

         if (n1dfft == lot3) then
           call dfftw_execute_dft(bw_plan3_lot, zw, zw)
         else
           call dfftw_execute_dft(bw_plan3_rest, zw, zw)
         end if

         ! input:  G1,R3,G2,(Gp2)
         ! output: G1,G2,R3,(Gp2)
         call scramble(i1,j2,lot3,n1dfft,n1,n3,nd2proc,nd3,zw,zmpi2)
       end do
     end if
   end do

   ! Interprocessor data transposition
   ! input:  G1,G2,R3,Rp3,(Gp2)
   ! output: G1,G2,G3,Gp2,(Rp3)
   if (nproc_fft>1) then
     call xmpi_alltoall(zmpi2,2*n1*nd2proc*nd3proc, &
&                       zmpi1,2*n1*nd2proc*nd3proc,comm_fft,ierr)
   end if

   do j3=1,nd3proc
     if (me_fft*nd3proc+j3 <= n3) then
       Jp2st=1; J2st=1

       ! transform along x axis
       do j=1,n2eff,lot1
         ma=j
         mb=min(j+(lot1-1),n2eff)
         n1dfft=mb-ma+1

         ! input:  G1,G2,R3,Gp2,(Rp3)
         ! output: G2,G1,R3,Jp2,(Rp3)
         if (nproc_fft == 1) then
           call mpiswitch(j3,n1dfft,Jp2st,J2st,lot1,n1,nd2proc,nd3proc,nproc_fft,option,zmpi2,zw)
         else
           call mpiswitch(j3,n1dfft,Jp2st,J2st,lot1,n1,nd2proc,nd3proc,nproc_fft,option,zmpi1,zw)
         end if

         ! input:  G2,G1,R3,(Rp3)
         ! output: G2,R1,R3,(Rp3)
         if (n1dfft == lot1) then
           call dfftw_execute_dft(bw_plan1_lot, zw, zt(1,j,1))
         else
           call dfftw_execute_dft(bw_plan1_rest, zw, zt(1,j,1))
         end if

       end do

       ! transform along y axis
       do j=1,n1eff,lot2
         ma=j
         mb=min(j+(lot2-1),n1eff)
         n1dfft=mb-ma+1
         includelast=1
         if (cplex==1) then
          jeff=2*j-1
          includelast=1
          if (mb==n1eff .and. n1eff*2/=n1) includelast=0
         end if

         ! input:  G2,R1,R3,(Rp3)
         ! output: R1,G2,R3,(Rp3)
         if (cplex==2) then
           call switch(n1dfft,n2,lot2,n1,lzt,zt(1,1,j),zw)
         else
           call switchreal(includelast,n1dfft,n2,n2eff,lot2,n1zt,lzt,zt(1,1,jeff),zw)
         end if

         if (n1dfft == lot2) then
           call dfftw_execute_dft(bw_plan2_lot, zw, zr(1,j,1,j3,idat))
         else
           call dfftw_execute_dft(bw_plan2_rest, zw, zr(1,j,1,j3,idat))
         end if
       end do
       ! output: R1,R2,R3,(Rp3)

     end if
   end do
 end do ! idat

 call dfftw_destroy_plan(bw_plan3_lot)
 if (mod(n1, lot3) /= 0) then
   call dfftw_destroy_plan(bw_plan3_rest)
 end if

 call dfftw_destroy_plan(bw_plan1_lot)
 if (mod(n2eff, lot1) /= 0) then
   call dfftw_destroy_plan(bw_plan1_rest)
 end if

 call dfftw_destroy_plan(bw_plan2_lot)
 if (mod(n1eff, lot2) /= 0) then
   call dfftw_destroy_plan(bw_plan2_rest)
 end if

 ABI_FREE(zmpi2)
 ABI_FREE(zw)
 ABI_FREE(zt)
 if (nproc_fft>1)  then
   ABI_FREE(zmpi1)
 end if

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/cplex,ndat,n1,n2,n3,nd1,nd2,nd1eff,nd2proc,nd3proc,option,comm_fft/))
 ABI_UNUSED((/zf(1,1,1,1,1),zr(1,1,1,1,1)/))
#endif

end subroutine fftw3_mpiback
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_mpiforw
!! NAME
!!  fftw3_mpiforw
!!
!! FUNCTION
!!   Adopt standard convention that isign=-1 for forward transform
!!   CALCULATES THE DISCRETE FOURIERTRANSFORM ZF(I1,I3,I2)=
!!   S_(j1,j2,j3) EXP(isign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) ZR(j1,j2,j3)
!!   in parallel using MPI/OpenMP and BLAS library calls.
!!
!! INPUTS
!!    ZR: input array
!!         ZR(1,i1,i2,i3,idat)=real(R(i1,i2,i3,idat))
!!         ZR(2,i1,i2,i3,idat)=imag(R(i1,i2,i3,idat))
!!         i1=1,n1 , i2=1,n2 , i3=1,n3 , idat=1,ndat
!! OUTPUTS
!!    ZF: output array (note the switch of i2 and i3)
!!         real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!!         imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!!         i1=1,n1 , i2=1,n2 , i3=1,n3 , idat=1,ndat
!!    nproc_fft: number of processors used as returned by MPI_COMM_SIZE
!!    me_fft: [0:nproc_fft-1] number of processor as returned by MPI_COMM_RANK
!!     n1,n2,n3: logical dimension of the transform. As transform lengths
!!               most products of the prime factors 2,3,5 are allowed.
!!              The detailed table with allowed transform lengths can
!!              be found in subroutine CTRIG
!!     nd1,nd2,nd3: Dimension of ZR and ZF
!!    nd2proc=((nd2-1)/nproc_fft)+1 maximal number of 2nd dim slices
!!    nd3proc=((nd3-1)/nproc_fft)+1 maximal number of 3rd dim slices
!!
!! NOTES
!!  SHOULD describe nd1eff
!!  SHOULD put cplex and nd1eff in OMP declarations
!!  SHOULD describe the change of value of nd2prod
!!
!!  The maximum number of processors that can reasonably be used is max(n2,n3)
!!
!!  It is very important to find the optimal
!!  value of NCACHE. NCACHE determines the size of the work array ZW, that
!!  has to fit into cache. It has therefore to be chosen to equal roughly
!!   half the size of the physical cache in units of real*8 numbers.
!!  The optimal value of ncache can easily be determined by numerical
!!  experimentation. A too large value of ncache leads to a dramatic
!!  and sudden decrease of performance, a too small value to a to a
!!  slow and less dramatic decrease of performance. If NCACHE is set
!!  to a value so small, that not even a single one dimensional transform
!!  can be done in the workarray zw, the program stops with an error message.
!!
!! PARENTS
!!      m_fftw3
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_mpiforw(cplex,ndat,n1,n2,n3,nd1,nd2,nd3,nd1eff,nd2proc,nd3proc,option,zr,zf,comm_fft)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,comm_fft
 integer,intent(in) :: ndat,n1,n2,n3,nd1,nd2,nd3,nd1eff,nd2proc,nd3proc,option
!arrays
 real(dp),intent(in) :: zr(2,nd1eff,nd2,nd3proc,ndat)
 real(dp),intent(out) :: zf(2,nd1,nd3,nd2proc,ndat)

!Local variables-------------------------------
!scalars
#ifdef HAVE_FFTW3
 integer :: j,i1,idat,ierr,j2,j2st,j3,jp2st,lzt,nthreads
 integer :: ma,mb,n1dfft,n1eff,n2eff,n1zt,ncache,nnd3,nproc_fft,me_fft,lot1,lot2,lot3
 integer(KIND_FFTW_PLAN) :: fw_plan1_lot,fw_plan1_rest
 integer(KIND_FFTW_PLAN) :: fw_plan2_lot,fw_plan2_rest
 integer(KIND_FFTW_PLAN) :: fw_plan3_lot,fw_plan3_rest
 character(len=500) :: msg
!arrays
 real(dp), allocatable :: zmpi1(:,:,:,:),zmpi2(:,:,:,:) ! work arrays for MPI
 real(dp),allocatable :: zw(:,:),zt(:,:,:) ! cache work array and array for transpositions

! *************************************************************************

 nproc_fft = xmpi_comm_size(comm_fft); me_fft = xmpi_comm_rank(comm_fft)

 ! find cache size that gives optimal performance on machine
 ncache=2*max(n1,n2,n3,1024)
 if (ncache/(2*max(n1,n2,n3))<1) then
   write(msg,'(5a)')&
&     'ncache has to be enlarged to be able to hold at',ch10, &
&     'least one 1-d FFT of each size even though this will',ch10,&
&     'reduce the performance for shorter transform lengths'
   MSG_ERROR(msg)
 end if

 ! check input
 if (nd1<n1 .or. nd2<n2 .or. nd3<n3) then
   MSG_ERROR("forw: assertion error nd1<n1 .or. nd2<n2 .or. nd3<n3")
 end if

!Effective n1 and n2 (complex-to-complex or real-to-complex)
 n1eff=n1; n2eff=n2; n1zt=n1
 if (cplex==1) then
   n1eff=(n1+1)/2; n2eff=n2/2+1; n1zt=2*(n1/2+1)
 end if

 lzt=n2eff
 if (mod(n2eff,2) == 0) lzt=lzt+1
 if (mod(n2eff,4) == 0) lzt=lzt+1

 ! maximal number of big box 3rd dim slices for all procs
 nnd3=nd3proc*nproc_fft

 ABI_MALLOC(zw,(2,ncache/2))
 ABI_MALLOC(zt,(2,lzt,n1zt))
 ABI_MALLOC(zmpi2,(2,n1,nd2proc,nnd3))
 if (nproc_fft>1)  then
   ABI_MALLOC(zmpi1,(2,n1,nd2proc,nnd3))
 end if

 ! Create plans.
 ! The prototype for sfftw_plan_many_dft is:
 ! sfftw_plan_many_dft(rank, n, howmany,
 !   fin,  iembed, istride, idist,
 !   fout, oembed, ostride, odist, isign, my_flags)

 lot1=ncache/(2*n1)
 lot2=ncache/(2*n2)
 lot3=ncache/(2*n3)

 nthreads = xomp_get_num_threads(open_parallel=.TRUE.)
 !nthreads = 1

 fw_plan3_lot = dplan_many_dft_2D(1, [n3], lot3, &
&    zw, [ncache/2], lot3, 1,                          &
&    zw, [ncache/2], lot3, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(n1, lot3) /= 0) then
   fw_plan3_rest = dplan_many_dft_2D(1, [n3], mod(n1, lot3), &
&    zw, [ncache/2], lot3, 1,                                      &
&    zw, [ncache/2], lot3, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 fw_plan1_lot = dplan_many_dft_2D(1, [n1], lot1, &
&    zt, [lzt, n1zt],   lzt,  1,                       &
&    zw, [ncache/2], lot1, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(n2eff, lot1) /= 0) then
   fw_plan1_rest = dplan_many_dft_2D(1, [n1], mod(n2eff, lot1), &
&    zt, [lzt, n1zt],   lzt, 1,                                       &
&    zw, [ncache/2], lot1, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 fw_plan2_lot = dplan_many_dft_2D(1, [n2], lot2, &
&    zr, [nd1eff,nd2,nd3proc,ndat], nd1eff, 1,         &
&    zw, [ncache/2], lot2, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(n1eff, lot2) /= 0) then
   fw_plan2_rest = dplan_many_dft_2D(1, [n2], mod(n1eff,lot2), &
&    zr, [nd1eff,nd2,nd3proc,ndat], nd1eff, 1,                       &
&    zw, [ncache/2], lot2, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 do idat=1,ndat

   do j3=1,nd3proc
     if (me_fft*(nd3proc)+j3 <= n3) then
       Jp2st=1; J2st=1

       ! transform along y axis
       ! input: R1,R2,R3,(Rp3)
       do j=1,n1eff,lot2
         ma=j
         mb=min(j+(lot2-1),n1eff)
         n1dfft=mb-ma+1

         if (n1dfft == lot2) then
           call dfftw_execute_dft(fw_plan2_lot,  zr(1,j,1,j3,idat), zw)
         else
           call dfftw_execute_dft(fw_plan2_rest, zr(1,j,1,j3,idat), zw)
         end if

         !  input: R1,G2,R3,(Rp3)
         ! output: G2,R1,R3,(Rp3)
         if (cplex==2) then
           call unswitch(n1dfft,n2,lot2,n1zt,lzt,zw,zt(1,1,j))
         else
           call unswitchreal(n1dfft,n2,n2eff,lot2,n1zt,lzt,zw,zt(1,1,2*j-1))
         end if
       end do

       ! transform along x axis
       ! input: G2,R1,R3,(Rp3)
       do j=1,n2eff,lot1
         ma=j
         mb=min(j+(lot1-1),n2eff)
         n1dfft=mb-ma+1

         if (n1dfft == lot1) then
           call dfftw_execute_dft(fw_plan1_lot,  zt(1,j,1), zw)
         else
           call dfftw_execute_dft(fw_plan1_rest, zt(1,j,1), zw)
         end if

         ! input:  G2,G1,R3,Gp2,(Rp3)
         ! output: G1,G2,R3,Gp2,(Rp3)
         ! write(std_out,*) 'J2st,Jp2st',J2st,Jp2st
         if (nproc_fft == 1) then
           call unmpiswitch(j3,n1dfft,Jp2st,J2st,lot1,n1,nd2proc,nd3proc,nproc_fft,option,zw,zmpi2)
         else
           call unmpiswitch(j3,n1dfft,Jp2st,J2st,lot1,n1,nd2proc,nd3proc,nproc_fft,option,zw,zmpi1)
         end if
       end do

     end if
   end do ! j3

   ! Interprocessor data transposition
   ! input:  G1,G2,R3,Gp2,(Rp3)
   ! output: G1,G2,R3,Rp3,(Gp2)
   if (nproc_fft>1) then
     call xmpi_alltoall(zmpi1,2*n1*nd2proc*nd3proc, &
&                       zmpi2,2*n1*nd2proc*nd3proc,comm_fft,ierr)
   end if

   ! transform along z axis
   ! input: G1,G2,R3,(Gp2)

   do j2=1,nd2proc
     if (me_fft*(nd2proc)+j2 <= n2eff) then
       do i1=1,n1,lot3
         ma=i1
         mb=min(i1+(lot3-1),n1)
         n1dfft=mb-ma+1

         ! input:  G1,G2,R3,(Gp2)
         ! output: G1,R3,G2,(Gp2)
         call unscramble(i1,j2,lot3,n1dfft,n1,n3,nd2proc,nd3,zmpi2,zw)

         if (n1dfft == lot3) then
           call dfftw_execute_dft(fw_plan3_lot, zw, zw)
         else
           call dfftw_execute_dft(fw_plan3_rest, zw, zw)
         end if

         call unfill(nd1,nd3,lot3,n1dfft,n3,zw,zf(1,i1,1,j2,idat))
         ! output: G1,G3,G2,(Gp2)
       end do
     end if
   end do

 end do ! idat

 call dfftw_destroy_plan(fw_plan3_lot)
 if (mod(n1, lot3) /= 0) then
   call dfftw_destroy_plan(fw_plan3_rest)
 end if

 call dfftw_destroy_plan(fw_plan1_lot)
 if (mod(n2eff, lot1) /= 0) then
   call dfftw_destroy_plan(fw_plan1_rest)
 end if

 call dfftw_destroy_plan(fw_plan2_lot)
 if (mod(n1eff, lot2) /= 0) then
   call dfftw_destroy_plan(fw_plan2_rest)
 end if

 ABI_FREE(zmpi2)
 ABI_FREE(zw)
 ABI_FREE(zt)
 if (nproc_fft>1)  then
   ABI_FREE(zmpi1)
 end if

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/cplex,ndat,n1,n2,n3,nd1,nd2,nd1eff,nd2proc,nd3proc,option,comm_fft/))
 ABI_UNUSED((/zf(1,1,1,1,1),zr(1,1,1,1,1)/))
#endif

end subroutine fftw3_mpiforw
!!***

!----------------------------------------------------------------------

!!****f* m_m_fftw3/fftw3_mpifourdp
!! NAME
!! fftw3_mpifourdp
!!
!! FUNCTION
!! Conduct Fourier transform of REAL or COMPLEX function f(r)=fofr defined on
!! fft grid in real space, to create complex f(G)=fofg defined on full fft grid
!! in reciprocal space, in full storage mode, or the reverse operation.
!! For the reverse operation, the final data is divided by nfftot.
!! REAL case when cplex=1, COMPLEX case when cplex=2
!! Usually used for density and potentials.
!!
!! INPUTS
!! cplex=1 if fofr is real, 2 if fofr is complex
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! ndat=Numbre of FFT transforms
!! isign=sign of Fourier transform exponent: current convention uses
!!    +1 for transforming from G to r
!!    -1 for transforming from r to G.
!! fftn2_distrib(2),ffti2_local(2)
!! fftn3_distrib(3),ffti3_local(3)
!! comm_fft=MPI communicator
!!
!! SIDE EFFECTS
!! Input/Output
!! fofg(2,nfft)=f(G), complex.
!! fofr(cplex*nfft)=input function f(r) (real or complex)
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_mpifourdp(cplex,nfft,ngfft,ndat,isign,&
&  fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,isign,nfft,ndat,comm_fft
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(in) :: fftn2_distrib(ngfft(2)),ffti2_local(ngfft(2))
 integer,intent(in) :: fftn3_distrib(ngfft(3)),ffti3_local(ngfft(3))
 real(dp),intent(inout) :: fofg(2,nfft*ndat),fofr(cplex*nfft*ndat)

!Local variables-------------------------------
!scalars
 integer :: n1,n2,n3,n4,n5,n6,nd2proc,nd3proc,nproc_fft,me_fft
!arrays
 real(dp),allocatable :: workf(:,:,:,:,:),workr(:,:,:,:,:)

! *************************************************************************

 ! Note the only c2c is supported in parallel.
 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 n4=ngfft(4); n5=ngfft(5); n6=ngfft(6)
 me_fft=ngfft(11); nproc_fft=ngfft(10)

 nd2proc=((n2-1)/nproc_fft) +1
 nd3proc=((n6-1)/nproc_fft) +1
 ABI_ALLOCATE(workr,(2,n4,n5,nd3proc,ndat))
 ABI_ALLOCATE(workf,(2,n4,n6,nd2proc,ndat))

 ! Complex to Complex
 ! TODO: Complex to Real
 select case (isign)
 case (1)
   ! G --> R
   call mpifft_fg2dbox(nfft,ndat,fofg,n1,n2,n3,n4,nd2proc,n6,fftn2_distrib,ffti2_local,me_fft,workf)

   call fftw3_mpiback(2,ndat,n1,n2,n3,n4,n5,n6,n4,nd2proc,nd3proc,2,workf,workr,comm_fft)

   call mpifft_dbox2fr(n1,n2,n3,n4,n5,nd3proc,ndat,fftn3_distrib,ffti3_local,me_fft,workr,cplex,nfft,fofr)

 case (-1)
   ! R --> G
   call mpifft_fr2dbox(cplex,nfft,ndat,fofr,n1,n2,n3,n4,n5,nd3proc,fftn3_distrib,ffti3_local,me_fft,workr)

   call fftw3_mpiforw(2,ndat,n1,n2,n3,n4,n5,n6,n4,nd2proc,nd3proc,2,workr,workf,comm_fft)

   ! Transfer FFT output to the original fft box.
   call mpifft_dbox2fg(n1,n2,n3,n4,nd2proc,n6,ndat,fftn2_distrib,ffti2_local,me_fft,workf,nfft,fofg)

 case default
   MSG_BUG("Wrong isign")
 end select

 ABI_DEALLOCATE(workr)
 ABI_DEALLOCATE(workf)

end subroutine fftw3_mpifourdp
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_applypot
!! NAME
!!  fftw3_applypot
!!
!! FUNCTION
!! Applies the local real space potential to multiple wavefunctions in Fourier space
!!
!! INPUTS
!!   ZF: Wavefunction (input/output) (note the switch of i2 and i3)
!!        real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!!        imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!!   max1 is positive or zero ; m1 >=max1+1
!!   i1= 1... max1+1 corresponds to positive and zero wavevectors 0 ... max1
!!   then, if m1 > max1+1, one has min1=max1-m1+1 and
!!   i1= max1+2 ... m1 corresponds to negative wavevectors min1 ... -1
!!   i2 and i3 have a similar definition of range
!!   idat=1,ndat
!!   md1,md2,md3: Dimension of ZF (input as well as output), distributed on different procs
!!   md2proc=((md2-1)/nproc_fft)+1  maximal number of small box 2nd dim slices for one proc
!!
!!   POT: Potential
!!        POT(cplex*i1,i2,i3)
!!        cplex=1 or 2 ,  i1=1,n1 , i2=1,n2 , i3=1,n3
!!   nd1,nd2,nd3: dimension of pot
!!   comm_fft: MPI communicator
!!   nproc_fft: number of processors used as returned by MPI_COMM_SIZE
!!   me_fft: [0:nproc_fft-1] number of processor as returned by MPI_COMM_RANK
!!    n1,n2,n3: logical dimension of the transform. As transform lengths
!!              most products of the prime factors 2,3,5 are allowed.
!!             The detailed table with allowed transform lengths can
!!             be found in subroutine CTRIG
!!
!! NOTES:
!!   PERFORMANCE CONSIDERATIONS:
!!   The maximum number of processors that can reasonably be used is max(n2/2,n3/2)
!!
!!   It is very important to find the optimal
!!   value of NCACHE. NCACHE determines the size of the work array ZW, that
!!   has to fit into cache. It has therefore to be chosen to equal roughly
!!    half the size of the physical cache in units of real*8 numbers.
!!   The optimal value of ncache can easily be determined by numerical
!!   experimentation. A too large value of ncache leads to a dramatic
!!   and sudden decrease of performance, a too small value to a to a
!!   slow and less dramatic decrease of performance. If NCACHE is set
!!   to a value so small, that not even a single one dimensional transform
!!   can be done in the workarray zw, the program stops with an error message.
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE


subroutine fftw3_applypot(cplexwf,cplex,ndat,n1,n2,n3,nd1,nd2,nd3,nd3proc,&
&  max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3,&
&  max1o,max2o,max3o,m1o,m2o,m3o,comm_fft,nproc_fft,me_fft,pot,zf)

!Arguments ------------------------------------
 integer,intent(in) :: cplexwf,cplex,ndat,n1,n2,n3,nd1,nd2,nd3,nd3proc
 integer,intent(in) :: max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3
 integer,intent(in) :: max1o,max2o,max3o,m1o,m2o,m3o,comm_fft,nproc_fft,me_fft
 real(dp),intent(in) :: pot(cplex*nd1,nd2,nd3)
 real(dp),intent(inout) :: zf(2,md1,md3,md2proc,ndat)

!Local variables-------------------------------
!scalars
#ifdef HAVE_FFTW3
 integer,parameter :: unused0=0
 integer :: j,i1,i2,i3,idat,ierr,j3glob,nthreads
 integer :: ioption,j2,j3,lzt,m1zt,ma,mb,n1dfft,nnd3,lot1,lot2,lot3
 integer :: m2eff,ncache,n1eff,i1inv,i2inv,i3inv,jeff,includelast,j2stb
 integer :: jx,j2stf,Jp2stb,Jp2stf,m2ieff,m2oeff
 integer(KIND_FFTW_PLAN) :: bw_plan1_lot,bw_plan1_rest
 integer(KIND_FFTW_PLAN) :: bw_plan2_lot,bw_plan2_rest
 integer(KIND_FFTW_PLAN) :: bw_plan3_lot,bw_plan3_rest
 integer(KIND_FFTW_PLAN) :: fw_plan1_lot,fw_plan1_rest
 integer(KIND_FFTW_PLAN) :: fw_plan2_lot,fw_plan2_rest
 integer(KIND_FFTW_PLAN) :: fw_plan3_lot,fw_plan3_rest
 character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 real(dp), allocatable :: zmpi1(:,:,:,:),zmpi2(:,:,:,:) ! work arrays for MPI
 real(dp),allocatable :: zw(:,:),zt(:,:,:) ! cache work array and array for transpositions
! FFT work arrays

! *************************************************************************

 !ioption=0 ! This was in the old version.
 ioption=1 ! This one is needed to be compatible with paral_kgb

 ncache=2*max(n1,n2,n3,1024)
 if (ncache/(2*max(n1,n2,n3)) < 1) then
   write(msg,"(5a)") &
&    'ncache has to be enlarged to be able to hold at',ch10,&
&    'least one 1-d FFT of each size even though this will',ch10,&
&    'reduce the performance for shorter transform lengths'
   MSG_ERROR(msg)
 end if

 !call wrtout(std_out,"applypot standard ALLTOALL + FFTW3","COLL")

 ! Effective m1 and m2 (complex-to-complex or real-to-complex)
 n1eff=n1; m2ieff=m2i; m2oeff=m2o; m1zt=n1
 if (cplexwf==1) then
   n1eff=(n1+1)/2; m2ieff=m2i/2+1; m2oeff=m2o/2+1; m1zt=2*(n1/2+1)
 end if

 m2eff=max(m2ieff,m2oeff)
 lzt=m2eff
 if (mod(m2eff,2) == 0) lzt=lzt+1
 if (mod(m2eff,4) == 0) lzt=lzt+1

 ! maximal number of big box 3rd dim slices for all procs
 nnd3=nd3proc*nproc_fft

 ABI_ALLOCATE(zw,(2,ncache/2))
 ABI_ALLOCATE(zt,(2,lzt,m1zt))
 ABI_ALLOCATE(zmpi2,(2,md1,md2proc,nnd3))
 if (nproc_fft > 1)  then
   ABI_ALLOCATE(zmpi1,(2,md1,md2proc,nnd3))
 end if

 lot3=ncache/(2*n3)
 lot1=ncache/(2*n1)
 lot2=ncache/(2*n2)

 ! The prototype for sfftw_plan_many_dft is:
 ! sfftw_plan_many_dft(rank, n, howmany,
 !   fin,  iembed, istride, idist,
 !   fout, oembed, ostride, odist, isign, my_flags)

 nthreads = xomp_get_num_threads(open_parallel=.TRUE.)
 !nthreads = 1

 ! Create plans for G --> R (see back_wf)
 bw_plan3_lot = dplan_many_dft_2D(1, [n3], lot3, &
&    zw, [ncache/2], lot3, 1,                          &
&    zw, [ncache/2], lot3, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(m1i, lot3) /= 0) then
   bw_plan3_rest = dplan_many_dft_2D(1, [n3], mod(m1i, lot3),&
&      zw, [ncache/2], lot3, 1,                                    &
&      zw, [ncache/2], lot3, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 bw_plan1_lot = dplan_many_dft_2D(1, [n1], lot1, &
&    zw, [ncache/2],  lot1, 1,                         &
&    zt, [lzt, m1zt], lzt, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(m2ieff, lot1) /= 0) then
   bw_plan1_rest = dplan_many_dft_2D(1, [n1], mod(m2ieff, lot1), &
&      zw, [ncache/2],  lot1, 1,                                       &
&      zt, [lzt, m1zt], lzt, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 !TODO this won't work if iclexwf==1
 ! Recheck this
 bw_plan2_lot = dplan_many_dft_2D(1, [n2], lot2, &
&    zw, [ncache/2], lot2, 1,                          &
&    zw, [ncache/2], lot2, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(n1eff, lot2) /= 0) then
   bw_plan2_rest = dplan_many_dft_2D(1, [n2], mod(n1eff,lot2), &
&      zw, [ncache/2], lot2, 1,                                      &
&      zw, [ncache/2], lot2, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 ! Create plans for G --> R (see forw_wf)
 fw_plan3_lot = dplan_many_dft_2D(1, [n3], lot3, &
&    zw, [ncache/2], lot3, 1,                          &
&    zw, [ncache/2], lot3, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(m1o, lot3) /= 0) then
   fw_plan3_rest = dplan_many_dft_2D(1, [n3], mod(m1o, lot3),&
&    zw, [ncache/2], lot3, 1,                                      &
&    zw, [ncache/2], lot3, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 fw_plan1_lot = dplan_many_dft_2D(1, [n1], lot1,&
&    zt, [lzt, m1zt], lzt,  1,                        &
&    zw, [ncache/2],  lot1, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(m2oeff, lot1) /= 0) then
   fw_plan1_rest = dplan_many_dft_2D(1, [n1], mod(m2oeff, lot1),&
&    zt, [lzt, m1zt], lzt,  1,                                        &
&    zw, [ncache/2],  lot1, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 fw_plan2_lot = dplan_many_dft_2D(1, [n2], lot2,&
&    zw, [ncache/2], lot2, 1,                         &
&    zw, [ncache/2], lot2, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(n1eff, lot2) /= 0) then
   fw_plan2_rest = dplan_many_dft_2D(1, [n2], mod(n1eff,lot2),&
&    zw, [ncache/2], lot2, 1,                                       &
&    zw, [ncache/2], lot2, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 do idat=1,ndat
   !
   ! transform along z axis
   ! input: G1,G3,G2,(Gp2)
   do j2=1,md2proc
     if (me_fft*md2proc+j2 <= m2ieff) then
       do i1=1,m1i,lot3
         ma=i1
         mb=min(i1+(lot3-1),m1i)
         n1dfft=mb-ma+1

         ! zero-pad n1dfft G_z lines
         ! input: G1,G3,G2,(Gp2)
         call fill_cent(md1,md3,lot3,n1dfft,max3i,m3i,n3,zf(1,i1,1,j2,idat),zw)

         if (n1dfft == lot3) then
           call dfftw_execute_dft(bw_plan3_lot, zw, zw)
         else
           call dfftw_execute_dft(bw_plan3_rest, zw, zw)
         end if

         ! Local rotation.
         ! input:  G1,R3,G2,(Gp2)
         ! output: G1,G2,R3,(Gp2)
         call scramble(i1,j2,lot3,n1dfft,md1,n3,md2proc,nnd3,zw,zmpi2)
       end do
     end if
   end do

   ! Interprocessor data transposition
   ! input:  G1,G2,R3,Rp3,(Gp2)
   ! output: G1,G2,R3,Gp2,(Rp3)
   if (nproc_fft > 1) then
     call timab(543,1,tsec)
     call xmpi_alltoall(zmpi2,2*md1*md2proc*nd3proc,&
&                       zmpi1,2*md1*md2proc*nd3proc,comm_fft,ierr)
     call timab(543,2,tsec)
   end if

   do j3=1,nd3proc
     j3glob = j3 + me_fft*nd3proc
     if (me_fft*nd3proc+j3 <= n3) then
       Jp2stb=1; J2stb=1
       Jp2stf=1; J2stf=1

       ! transform along x axis
       do j=1,m2ieff,lot1
         ma=j
         mb=min(j+(lot1-1),m2ieff)
         n1dfft=mb-ma+1

         ! Zero-pad input.
         ! input:  G1,G2,R3,G2,(Rp3)
         ! output: G2,G1,R3,G2,(Rp3)
         if (nproc_fft == 1) then
           call mpiswitch_cent(j3,n1dfft,Jp2stb,J2stb,lot1,max1i,md1,m1i,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zmpi2,zw, unused0, unused0, unused0)
         else
           call mpiswitch_cent(j3,n1dfft,Jp2stb,J2stb,lot1,max1i,md1,m1i,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zmpi1,zw, unused0, unused0, unused0)
         end if

         ! Transform along x
         ! input:  G2,G1,R3,(Rp3)
         ! output: G2,R1,R3,(Rp3)
         if (n1dfft == lot1) then
           call dfftw_execute_dft(bw_plan1_lot, zw, zt(1,j,1))
         else
           call dfftw_execute_dft(bw_plan1_rest, zw, zt(1,j,1))
         end if
       end do

       ! Transform along y axis (take into account c2c or c2r case).
       ! Must loop over the full box.
       !TODO this won't work
       if (cplexwf==1) then
         if(mod(lot2,2).ne.0) lot2=lot2-1 ! needed to introduce jeff
       end if

       do j=1,n1eff,lot2
         ma=j
         mb=min(j+(lot2-1),n1eff)
         n1dfft=mb-ma+1
         jeff=j
         includelast=1

         if (cplexwf==1) then
           jeff=2*j-1
           includelast=1
           if (mb==n1eff .and. n1eff*2/=n1) includelast=0
         end if

         ! Zero-pad the input.
         !  input: G2,R1,R3,(Rp3)
         ! output: R1,G2,R3,(Rp3)
         if (cplexwf==2) then
           call switch_cent(n1dfft,max2i,m2i,n2,lot2,n1,lzt,zt(1,1,jeff),zw)
         else
           call switchreal_cent(includelast,n1dfft,max2i,n2,lot2,m1zt,lzt,zt(1,1,jeff),zw)
         end if

         ! input:  R1,G2,R3,(Rp3)
         ! output: R1,R2,R3,(Rp3)
         ! Be careful here
         if (n1dfft == lot2) then
           call dfftw_execute_dft(bw_plan2_lot, zw, zw)
         else
           call dfftw_execute_dft(bw_plan2_rest, zw, zw)
         end if

         ! Multiply with potential in real space
         jx=cplex*(jeff-1)+1
         call multpot(cplexwf,cplex,includelast,nd1,nd2,n2,lot2,n1dfft,pot(jx,1,j3glob),zw)

         ! TRANSFORM BACK IN FOURIER SPACE
         ! transform along y axis
         ! input: R1,R2,R3,(Rp3)
         if (n1dfft == lot2) then
           call dfftw_execute_dft(fw_plan2_lot,  zw, zw)
         else
           call dfftw_execute_dft(fw_plan2_rest, zw, zw)
         end if

         ! input: R1,G2,R3,(Rp3)
         ! output: G2,R1,R3,(Rp3)
         if (cplexwf==2) then
           call unswitch_cent(n1dfft,max2o,m2o,n2,lot2,n1,lzt,zw,zt(1,1,jeff))
         else
           call unswitchreal_cent(n1dfft,max2o,n2,lot2,n1,lzt,zw,zt(1,1,jeff))
         end if
       end do ! j

       ! transform along x axis
       ! input:  R2,R1,R3,(Rp3)
       ! output: R2,G1,R3,(Rp3)
       do j=1,m2oeff,lot1
         ma=j
         mb=min(j+(lot1-1),m2oeff)
         n1dfft=mb-ma+1

         if (n1dfft == lot1) then
           call dfftw_execute_dft(fw_plan1_lot,  zt(1,j,1), zw)
         else
           call dfftw_execute_dft(fw_plan1_rest, zt(1,j,1), zw)
         end if

         ! input:  G2,G1,R3,Gp2,(Rp3)
         ! output: G1,G2,R3,Gp2,(Rp3)
         if (nproc_fft == 1) then
           call unmpiswitch_cent(j3,n1dfft,Jp2stf,J2stf,lot1,max1o,md1,m1o,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zw,zmpi2)
         else
           call unmpiswitch_cent(j3,n1dfft,Jp2stf,J2stf,lot1,max1o,md1,m1o,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zw,zmpi1)
         end if
       end do ! j
     end if
   end do

   ! Interprocessor data transposition
   ! input:  G1,G2,R3,Gp2,(Rp3)
   ! output: G1,G2,R3,Rp3,(Gp2)
   if (nproc_fft > 1) then
     call timab(544,1,tsec)
     call xmpi_alltoall(zmpi1,2*md1*md2proc*nd3proc, &
&                       zmpi2,2*md1*md2proc*nd3proc,comm_fft,ierr)
     call timab(544,2,tsec)
   end if

   ! transform along z axis
   ! input: G1,G2,R3,(Gp2)
   !lot=ncache/(4*n3)
   do j2=1,md2proc
     if (me_fft*md2proc+j2 <= m2oeff) then
       do i1=1,m1o,lot3
         ma=i1
         mb=min(i1+(lot3-1),m1o)
         n1dfft=mb-ma+1

         ! input:  G1,G2,R3,(Gp2)
         ! output: G1,R3,G2,(Gp2)
         call unscramble(i1,j2,lot3,n1dfft,md1,n3,md2proc,nnd3,zmpi2,zw)

          if (n1dfft == lot3) then
            call dfftw_execute_dft(fw_plan3_lot, zw, zw)
          else
            call dfftw_execute_dft(fw_plan3_rest, zw, zw)
          end if

         call unfill_cent(md1,md3,lot3,n1dfft,max3o,m3o,n3,zw,zf(1,i1,1,j2,idat))
         ! output: G1,G3,G2,(Gp2)
       end do
     end if
   end do

   ! Complete missing values with complex conjugate
   ! Inverse of ix is located at nx+2-ix , except for ix=1, for which it is 1.
   if (cplexwf==1) then
     do i3=1,m3o
       i3inv=m3o+2-i3
       if (i3==1) i3inv=1
       if (m2oeff>1)then
         do i2=2,m2oeff
           i2inv=m2o+2-i2
           zf(1,1,i3inv,i2inv,idat)= zf(1,1,i3,i2,idat)
           zf(2,1,i3inv,i2inv,idat)=-zf(2,1,i3,i2,idat)
           do i1=2,m1o
             i1inv=m1o+2-i1
             zf(1,i1inv,i3inv,i2inv,idat)= zf(1,i1,i3,i2,idat)
             zf(2,i1inv,i3inv,i2inv,idat)=-zf(2,i1,i3,i2,idat)
           end do
         end do
       end if
     end do
   end if

 end do ! idat

 call dfftw_destroy_plan(bw_plan3_lot)
 if (mod(m1i, lot3) /= 0) then
   call dfftw_destroy_plan(bw_plan3_rest)
 end if

 call dfftw_destroy_plan(bw_plan1_lot)
 if (mod(m2ieff, lot1) /= 0) then
   call dfftw_destroy_plan(bw_plan1_rest)
 end if

 call dfftw_destroy_plan(bw_plan2_lot)
 if (mod(n1eff, lot2) /= 0) then
   call dfftw_destroy_plan(bw_plan2_rest)
 end if

 call dfftw_destroy_plan(fw_plan3_lot)
 if (mod(m1o, lot3) /= 0) then
   call dfftw_destroy_plan(fw_plan3_rest)
 end if

 call dfftw_destroy_plan(fw_plan1_lot)
 if (mod(m2oeff, lot1) /= 0) then
   call dfftw_destroy_plan(fw_plan1_rest)
 end if

 call dfftw_destroy_plan(fw_plan2_lot)
 if (mod(n1eff, lot2) /= 0) then
   call dfftw_destroy_plan(fw_plan2_rest)
 end if

 ABI_DEALLOCATE(zmpi2)
 ABI_DEALLOCATE(zw)
 ABI_DEALLOCATE(zt)
 if (nproc_fft > 1)  then
   ABI_DEALLOCATE(zmpi1)
 end if

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/cplexwf,cplex,ndat,n1,n2,n3,nd1,nd2,nd3,nd3proc/))
 ABI_UNUSED((/max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3/))
 ABI_UNUSED((/max1o,max2o,max3o,m1o,m2o,m3o,comm_fft,nproc_fft,me_fft/))
 ABI_UNUSED((/pot(1,1,1),zf(1,1,1,1,1)/))
#endif

end subroutine fftw3_applypot
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_accrho
!! NAME
!! fftw3_accrho
!!
!! FUNCTION
!! Accumulates the real space density rho from the ndat wavefunctions zf
!! by transforming zf into real space and adding all the amplitudes squared
!!
!! INPUTS:
!!   ZF: input array (note the switch of i2 and i3)
!!         real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!!         imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!!   max1 is positive or zero ; m1 >=max1+1
!!   i1= 1... max1+1 corresponds to positive and zero wavevectors 0 ... max1
!!   then, if m1 > max1+1, one has min1=max1-m1+1 and
!!   i1= max1+2 ... m1 corresponds to negative wavevectors min1 ... -1
!!   i2 and i3 have a similar definition of range
!!   idat=1,ndat
!!   md1,md2,md3: Dimension of ZF
!!   md2proc=((md2-1)/nproc_fft)+1 ! maximal number of small box 2nd dim slices for one proc
!!   weight(ndat)= weight for the density accumulation
!!
!! OUTPUTS:
!!    RHOoutput(i1,i2,i3) = RHOinput(i1,i2,i3) + sum on idat of (FFT(ZF))**2 *weight
!!        i1=1,n1 , i2=1,n2 , i3=1,n3
!!   comm_fft: MPI communicator
!!   nproc_fft: number of processors used as returned by MPI_COMM_SIZE
!!   me_fft: [0:nproc_fft-1] number of processor as returned by MPI_COMM_RANK
!!    n1,n2,n3: logical dimension of the transform. As transform lengths
!!              most products of the prime factors 2,3,5 are allowed.
!!             The detailed table with allowed transform lengths can
!!             be found in subroutine CTRIG
!!    nd1,nd2,nd3: Dimension of RHO
!!   nd3proc=((nd3-1)/nproc_fft)+1 ! maximal number of big box 3rd dim slices for one proc
!!
!! NOTES:
!!   PERFORMANCE CONSIDERATIONS:
!!   The maximum number of processors that can reasonably be used is max(n2/2,n3/2)
!!
!!   It is very important to find the optimal
!!   value of NCACHE. NCACHE determines the size of the work array ZW, that
!!   has to fit into cache. It has therefore to be chosen to equal roughly
!!    half the size of the physical cache in units of real*8 numbers.
!!   The optimal value of ncache can easily be determined by numerical
!!   experimentation. A too large value of ncache leads to a dramatic
!!   and sudden decrease of performance, a too small value to a to a
!!   slow and less dramatic decrease of performance. If NCACHE is set
!!   to a value so small, that not even a single one dimensional transform
!!   can be done in the workarray zw, the program stops with an error message.
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_accrho(cplexwf,ndat,n1,n2,n3,nd1,nd2,nd3,nd3proc,&
&  max1,max2,max3,m1,m2,m3,md1,md2proc,md3,comm_fft,nproc_fft,me_fft,zf,rho,weight_r,weight_i)

!Arguments ------------------------------------
 integer,intent(in) :: cplexwf,ndat,n1,n2,n3,nd1,nd2,nd3,nd3proc
 integer,intent(in) :: max1,max2,max3,m1,m2,m3,md1,md2proc,md3,comm_fft,nproc_fft,me_fft
 real(dp),intent(in) :: zf(2,md1,md3,md2proc,ndat)
 real(dp),intent(in) :: weight_r(ndat) , weight_i(ndat)
 real(dp),intent(inout) :: rho(nd1,nd2,nd3)

!Local variables-------------------------------
!scalars
#ifdef HAVE_FFTW3
 integer,parameter :: unused0=0
 integer :: j,i1,idat,ierr,j3glob
 integer :: ioption,j2,j3,j2st,jp2st,lzt,m1zt,ma,mb,n1dfft,nnd3
 integer :: m2eff,ncache,n1eff,jeff,includelast,lot1,lot2,lot3,nthreads
 integer(KIND_FFTW_PLAN) :: bw_plan1_lot,bw_plan1_rest
 integer(KIND_FFTW_PLAN) :: bw_plan2_lot,bw_plan2_rest
 integer(KIND_FFTW_PLAN) :: bw_plan3_lot,bw_plan3_rest
 character(len=500) :: msg
!arrays
 real(dp), allocatable :: zmpi1(:,:,:,:),zmpi2(:,:,:,:) ! work arrays for MPI
 real(dp),allocatable :: zw(:,:),zt(:,:,:) ! cache work array and array for transpositions
 real(dp) :: tsec(2)

! *************************************************************************

 !ioption=0 ! This was in the old version.
 ioption=1 ! This one is needed to be compatible with paral_kgb

 !nproc_fft = xmpi_comm_size(comm_fft); me_fft = xmpi_comm_rank(comm_fft)

! find cache size that gives optimal performance on machine
 ncache=2*max(n1,n2,n3,1024)
 if (ncache/(2*max(n1,n2,n3)) < 1) then
    write(msg,"(5a)") &
&     'ncache has to be enlarged to be able to hold at',ch10,&
&     'least one 1-d FFT of each size even though this will',ch10,&
&     'reduce the performance for shorter transform lengths'
    MSG_ERROR(msg)
 end if

!Effective m1 and m2 (complex-to-complex or real-to-complex)
 n1eff=n1; m2eff=m2 ; m1zt=n1
 if (cplexwf==1) then
   n1eff=(n1+1)/2; m2eff=m2/2+1; m1zt=2*(n1/2+1)
 end if

 lzt=m2eff
 if (mod(m2eff,2) == 0) lzt=lzt+1
 if (mod(m2eff,4) == 0) lzt=lzt+1

 ! maximal number of big box 3rd dim slices for all procs
 nnd3=nd3proc*nproc_fft

 ABI_ALLOCATE(zw,(2,ncache/2))
 ABI_ALLOCATE(zt,(2,lzt,m1zt))
 ABI_ALLOCATE(zmpi2,(2,md1,md2proc,nnd3))
 if (nproc_fft > 1)  then
   ABI_ALLOCATE(zmpi1,(2,md1,md2proc,nnd3))
 end if

 ! Create plans.
 ! The prototype for sfftw_plan_many_dft is:
 ! sfftw_plan_many_dft(rank, n, howmany,
 !   fin,  iembed, istride, idist,
 !   fout, oembed, ostride, odist, isign, my_flags)

 lot3=ncache/(2*n3)
 lot1=ncache/(2*n1)
 lot2=ncache/(2*n2)

 nthreads = xomp_get_num_threads(open_parallel=.TRUE.)
 !nthreads = 1

 bw_plan3_lot = dplan_many_dft_2D(1, [n3], lot3, &
&    zw, [ncache/2], lot3, 1,                          &
&    zw, [ncache/2], lot3, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(m1, lot3) /= 0) then
   bw_plan3_rest = dplan_many_dft_2D(1, [n3], mod(m1, lot3), &
&      zw, [ncache/2], lot3, 1,                                    &
&      zw, [ncache/2], lot3, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 bw_plan1_lot = dplan_many_dft_2D(1, [n1], lot1, &
&    zw, [ncache/2],  lot1, 1,                         &
&    zt, [lzt, m1zt], lzt, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(m2eff, lot1) /= 0) then
   bw_plan1_rest = dplan_many_dft_2D(1, [n1], mod(m2eff, lot1), &
&      zw, [ncache/2],  lot1, 1,                                      &
&      zt, [lzt, m1zt], lzt, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 ! FIXME THis won't work if ixplexwf == 1
 bw_plan2_lot = dplan_many_dft_2D(1, [n2], lot2, &
&    zw, [ncache/2], lot2, 1,                          &
&    zw, [ncache/2], lot2, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(n1eff, lot2) /= 0) then
   bw_plan2_rest = dplan_many_dft_2D(1, [n2], mod(n1eff,lot2), &
&      zw, [ncache/2], lot2, 1,                                      &
&      zw, [ncache/2], lot2, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 do idat=1,ndat
   ! transform along z axis
   ! input: I1,I3,J2,(Jp2)
   !lot=ncache/(4*n3)

   ! Loop over the y planes treated by this node and trasform n1ddft G_z lines.
   do j2=1,md2proc
     if (me_fft*md2proc+j2 <= m2eff) then ! MG REMOVED TO BE COSISTENT WITH BACK_WF
       do i1=1,m1,lot3
         ma=i1
         mb=min(i1+(lot3-1),m1)
         n1dfft=mb-ma+1

         ! zero-pad n1dfft G_z lines
         !  input: G1,G3,G2,(Gp2)
         ! output: G1,R3,G2,(Gp2)
         call fill_cent(md1,md3,lot3,n1dfft,max3,m3,n3,zf(1,i1,1,j2,idat),zw)

         ! Transform along z.
         if (n1dfft == lot3) then
           call dfftw_execute_dft(bw_plan3_lot, zw, zw)
         else
           call dfftw_execute_dft(bw_plan3_rest, zw, zw)
         end if

         ! Local rotation.
         ! input:  G1,R3,G2,(Gp2)
         ! output: G1,G2,R3,(Gp2)
         call scramble(i1,j2,lot3,n1dfft,md1,n3,md2proc,nnd3,zw,zmpi2)
       end do
     end if
   end do

   ! Interprocessor data transposition
   ! input:  G1,G2,R3,Rp3,(Gp2)
   ! output: G1,G2,R3,Gp2,(Rp3)
   if (nproc_fft > 1) then
     call timab(543,1,tsec)
     call xmpi_alltoall(zmpi2,2*md1*md2proc*nd3proc, &
&                       zmpi1,2*md1*md2proc*nd3proc,comm_fft,ierr)
     call timab(543,2,tsec)
   end if

   ! Loop over the z treated by this node.
   do j3=1,nd3proc
     j3glob = j3 + me_fft*nd3proc

     if (me_fft*nd3proc+j3 <= n3) then
       Jp2st=1; J2st=1

       ! Loop over G_y in the small box.
       do j=1,m2eff,lot1
         ma=j
         mb=min(j+(lot1-1),m2eff)
         n1dfft=mb-ma+1

         ! Zero-pad input.
         ! input:  G1,G2,R3,JG2,(Rp3)
         ! output: G2,G1,R3,JG2,(Rp3)
         if (nproc_fft == 1) then
          call mpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot1,max1,md1,m1,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zmpi2,zw,unused0, unused0,unused0)
         else
          call mpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot1,max1,md1,m1,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zmpi1,zw, unused0,unused0,unused0)
         end if

         ! Transform along x
         ! input:  G2,G1,R3,(Rp3)
         ! output: G2,R1,R3,(Rp3)
         if (n1dfft == lot1) then
           call dfftw_execute_dft(bw_plan1_lot, zw, zt(1,j,1))
         else
           call dfftw_execute_dft(bw_plan1_rest, zw, zt(1,j,1))
         end if

       end do

       ! Transform along y axis (take into account c2c or c2r case).
       ! Must loop over the full box.
       !lot=ncache/(4*n2)
       ! FIXME THis won't work
       if (cplexwf==1) then
         if (mod(lot2,2) /=0) lot2=lot2-1 ! needed to introduce jeff
       end if

       do j=1,n1eff,lot2
         ma=j
         mb=min(j+(lot2-1),n1eff)
         n1dfft=mb-ma+1
         jeff=j
         includelast=1

         if (cplexwf==1) then
           jeff=2*j-1
           includelast=1
           if (mb==n1eff .and. n1eff*2/=n1) includelast=0
         end if

         ! Zero-pad the input.
         ! input:  G2,R1,R3,(Rp3)
         ! output: R1,G2,R3,(Rp3)
         if (cplexwf==2) then
           call switch_cent(n1dfft,max2,m2,n2,lot2,n1,lzt,zt(1,1,j),zw)
         else
           call switchreal_cent(includelast,n1dfft,max2,n2,lot2,m1zt,lzt,zt(1,1,jeff),zw)
         end if

         if (n1dfft == lot2) then
           call dfftw_execute_dft(bw_plan2_lot, zw, zw)
         else
           call dfftw_execute_dft(bw_plan2_rest, zw, zw)
         end if

         ! Accumulate
         call addrho(cplexwf,includelast,nd1,nd2,n2,lot2,n1dfft,&
&          zw,rho(jeff,1,j3glob),weight_r(idat),weight_i(idat))
       end do
       ! output: i1,i2,j3,(jp3)

      end if
    end do ! j3
 end do ! idat

 call dfftw_destroy_plan(bw_plan3_lot)
 if (mod(m1, lot3) /= 0) then
   call dfftw_destroy_plan(bw_plan3_rest)
 end if

 call dfftw_destroy_plan(bw_plan1_lot)
 if (mod(m2eff, lot1) /= 0) then
   call dfftw_destroy_plan(bw_plan1_rest)
 end if

 call dfftw_destroy_plan(bw_plan2_lot)
 if (mod(n1eff, lot2) /= 0) then
   call dfftw_destroy_plan(bw_plan2_rest)
 end if

 ABI_DEALLOCATE(zmpi2)
 ABI_DEALLOCATE(zw)
 ABI_DEALLOCATE(zt)
 if (nproc_fft > 1)  then
   ABI_DEALLOCATE(zmpi1)
 end if

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/cplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc/))
 ABI_UNUSED((/ max1,max2,max3,m1,m2,m3,md1,md2proc,md3,comm_fft,nproc_fft,me_fft/))
 ABI_UNUSED((/zf(1,1,1,1,1),rho(1,1,1),weight_r(1),weight_i(1)/))
#endif

end subroutine fftw3_accrho
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_mpiback_manywf
!! NAME
!!  fftw3_mpiback_manywf
!!
!! FUNCTION
!!   Does multiple 3-dim backward FFTs from Fourier into real space
!!   Adopt standard convention that isign=1 for backward transform
!!
!!   CALCULATES THE DISCRETE FOURIER TRANSFORM ZF(I1,I2,I3)=
!!
!!   S_(j1,j2,j3) EXP(isign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) ZF(j1,j3,j2)
!!
!!   in parallel using MPI/OpenMP.
!!
!! INPUTS:
!!    cplexwf=1 if wavefunction is real, 2 if complex
!!    ndat=Number of wavefunctions to transform.
!!    n1,n2,n3: logical dimension of the transform. As transform lengths
!!              most products of the prime factors 2,3,5 are allowed.
!!              The detailed table with allowed transform lengths can be found in subroutine CTRIG
!!    nd1,nd2,nd3: Leading Dimension of ZR
!!    nd3proc=((nd3-1)/nproc_fft)+1 maximal number of big box 3rd dim slices for one proc
!!    max1 is positive or zero; m1 >=max1+1
!!      i1= 1... max1+1 corresponds to positive and zero wavevectors 0 ... max1
!!      then, if m1 > max1+1, one has min1=max1-m1+1 and
!!      i1= max1+2 ... m1 corresponds to negative wavevectors min1 ... -1
!!    max2 and max3 have a similar definition of range
!!    m1,m2,m3=Size of the box enclosing the G-sphere.
!!    md1,md2,md3: Dimension of ZF given on the **small** FFT box.
!!    md2proc=((md2-1)/nproc_fft)+1 maximal number of small box 2nd dim slices for one proc
!!    nproc_fft: number of processors used as returned by MPI_COMM_SIZE
!!    comm_fft=MPI communicator for the FFT.
!!    ZF: input array (note the switch of i2 and i3)
!!          real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!!          imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!!
!! OUTPUTS
!!    ZR: output array
!!          ZR(1,i1,i2,i3,idat)=real(R(i1,i2,i3,idat))
!!          ZR(2,i1,i2,i3,idat)=imag(R(i1,i2,i3,idat))
!!        i1=1,n1 , i2=1,n2 , i3=1,n3 , idat=1,ndat
!!
!! NOTES
!!   The maximum number of processors that can reasonably be used is max(n2/2,n3/2)
!!
!!   It is very important to find the optimal
!!   value of NCACHE. NCACHE determines the size of the work array ZW, that
!!   has to fit into cache. It has therefore to be chosen to equal roughly
!!   half the size of the physical cache in units of real*8 numbers.
!!   The optimal value of ncache can easily be determined by numerical
!!   experimentation. A too large value of ncache leads to a dramatic
!!   and sudden decrease of performance, a too small value to a to a
!!   slow and less dramatic decrease of performance. If NCACHE is set
!!   to a value so small, that not even a single one dimensional transform
!!   can be done in the workarray zw, the program stops with an error message.
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_mpiback_manywf(cplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc,&
&  max1,max2,max3,m1,m2,m3,md1,md2proc,md3,zf,zr,comm_fft)

!Arguments ------------------------------------
 integer,intent(in) :: cplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc
 integer,intent(in) :: max1,max2,max3,m1,m2,m3,md1,md2proc,md3,comm_fft
 real(dp),intent(in) :: zf(2,md1,md3,md2proc,ndat)
 real(dp),intent(out) :: zr(2,nd1,nd2,nd3proc,ndat)

#ifdef HAVE_FFTW3
!Local variables-------------------------------
 integer,parameter :: nt1=1
 integer :: j,i1,i2,idat,ierr,includelast,nthreads
 integer :: ioption,j2,j3,j2st,jp2st,jeff,lzt,m1zt,ma,mb,n1dfft,nnd3
 integer :: lot1,lot2,lot3
 integer :: m2eff,ncache,n1eff,n1half,nproc_fft,me_fft
 integer(KIND_FFTW_PLAN) :: bw_plan1_lot,bw_plan1_rest
 integer(KIND_FFTW_PLAN) :: bw_plan2_lot,bw_plan2_rest
 integer(KIND_FFTW_PLAN) :: bw_plan3_lot,bw_plan3_rest
 !type(C_PTR) :: zw_cptr,zt_cptr
 character(len=500) :: msg
!arrays
 integer :: requests(ndat)
 real(dp) ABI_ASYNC, allocatable :: zmpi1(:,:,:,:,:),zmpi2(:,:,:,:,:)  ! work arrays for MPI
 real(dp),allocatable :: zw(:,:),zt(:,:,:) ! cache work array and array for transpositions
 !real(dp),ABI_CONTIGUOUS pointer :: zw(:,:),zt(:,:,:)
! FFT work arrays
 real(dp) :: tsec(2)

! *************************************************************************

 !call wrtout(std_out,"mpiback with non-blocking IALLTOALL + FFTW3","COLL")


 ! FIXME must provide a default value but which one?
 ! ioption = 0
 ioption = 1
 !if (paral_kgb==1) ioption=1

 nproc_fft = xmpi_comm_size(comm_fft); me_fft = xmpi_comm_rank(comm_fft)

 ! Find cache size that gives optimal performance on machine
 ncache=2*max(n1,n2,n3,1024)
 if (ncache/(2*max(n1,n2,n3))<1) then
   write(msg,"(5a)") &
&    'ncache has to be enlarged to be able to hold at',ch10, &
&    'least one 1-d FFT of each size even though this will',ch10,&
&    'reduce the performance for shorter transform lengths'
    MSG_ERROR(msg)
 end if

 ! Effective m1 and m2 (complex-to-complex or real-to-complex)
 n1eff=n1; m2eff=m2; m1zt=n1
 if (cplexwf==1) then
   n1eff=(n1+1)/2; m2eff=m2/2+1; m1zt=2*(n1/2+1)
 end if

 lzt=m2eff
 if (mod(m2eff,2)==0) lzt=lzt+1
 if (mod(m2eff,4)==0) lzt=lzt+1

 ! maximal number of big box 3rd dim slices for all procs
 nnd3=nd3proc*nproc_fft

 ! Allocate cache work array and work arrays for MPI transpositions.
 ABI_MALLOC(zw,(2,ncache/2))
 ABI_MALLOC(zt,(2,lzt,m1zt))

 !call fftw3_alloc_real([2,ncache/2],zw_cptr,zw)
 !call fftw3_alloc_real([2,lzt,m1zt],zt_cptr,zt)

 ABI_MALLOC(zmpi2,(2,md1,md2proc,nnd3,ndat))
 if (nproc_fft>1)  then
   ABI_MALLOC(zmpi1,(2,md1,md2proc,nnd3,ndat))
 end if

 ! Create plans.
 ! The prototype for sfftw_plan_many_dft is:
 ! sfftw_plan_many_dft(rank, n, howmany,
 !   fin,  iembed, istride, idist,
 !   fout, oembed, ostride, odist, isign, my_flags)

 lot3=ncache/(2*n3)
 lot1=ncache/(2*n1)
 lot2=ncache/(2*n2)

 nthreads = xomp_get_num_threads(open_parallel=.TRUE.)
 !nthreads = 1

 bw_plan3_lot = dplan_many_dft_2D(1, [n3], lot3, &
&    zw, [ncache/2], lot3, 1,                          &
&    zw, [ncache/2], lot3, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(m1, lot3) /= 0) then
   bw_plan3_rest = dplan_many_dft_2D(1, [n3], mod(m1, lot3), &
&      zw, [ncache/2], lot3, 1,                                    &
&      zw, [ncache/2], lot3, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 bw_plan1_lot = dplan_many_dft_2D(1, [n1], lot1, &
&    zw, [ncache/2],  lot1, 1,                         &
&    zt, [lzt, m1zt], lzt, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(m2eff, lot1) /= 0) then
   bw_plan1_rest = dplan_many_dft_2D(1, [n1], mod(m2eff, lot1), &
&      zw, [ncache/2],  lot1, 1,                                      &
&      zt, [lzt, m1zt], lzt, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 bw_plan2_lot = dplan_many_dft_2D(1, [n2], lot2, &
&    zw, [ncache/2], lot2, 1,                          &
&    zr, [nd1,nd2,nd3proc,ndat], nd1, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(n1eff, lot2) /= 0) then
   bw_plan2_rest = dplan_many_dft_2D(1, [n2], mod(n1eff,lot2), &
&      zw, [ncache/2], lot2, 1,                                      &
&      zr, [nd1,nd2,nd3proc,ndat], nd1, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 do idat=1,ndat
    ! transform along z axis
    ! input: G1,G3,G2,(Gp2)

    ! Loop over the y planes treated by this node and trasform n1ddft G_z lines.
    do j2=1,md2proc
      ! if (me_fft*md2proc+j2<=m2eff) then !a faire plus tard
      do i1=1,m1,lot3
        ma=i1
        mb=min(i1+(lot3-1),m1)
        n1dfft=mb-ma+1

        ! zero-pad n1dfft G_z lines
        ! input:  G1,G3,G2,(Gp2)
        ! output: G1,R3,G2,(Gp2)
        call fill_cent(md1,md3,lot3,n1dfft,max3,m3,n3,zf(1,i1,1,j2,idat),zw)

        ! Transform along z.
        if (n1dfft == lot3) then
          call dfftw_execute_dft(bw_plan3_lot, zw, zw)
        else
          call dfftw_execute_dft(bw_plan3_rest, zw, zw)
        end if

        ! Local rotation.
        ! input:  G1,R3,G2,(Gp2)
        ! output: G1,G2,R3,(Gp2)
        call scramble(i1,j2,lot3,n1dfft,md1,n3,md2proc,nnd3,zw,zmpi2(:,:,:,:,idat))
      end do
    end do ! j2

    ! Interprocessor data transposition
    ! input:  G1,G2,R3,Rp3,(Gp2)
    ! output: G1,G2,R3,Gp2,(Rp3)
    if (nproc_fft>1) then
      call timab(543,1,tsec)
      call xmpi_ialltoall(zmpi2(:,:,:,:,idat),2*md1*md2proc*nd3proc, &
&                         zmpi1(:,:,:,:,idat),2*md1*md2proc*nd3proc,comm_fft,requests(idat))
      call timab(543,2,tsec)
    end if
 end do

 do idat=1,ndat
    if (nproc_fft>1) call xmpi_wait(requests(idat),ierr)
    ! Loop over the z treated by this node.
    do j3=1,nd3proc
      if (me_fft*nd3proc+j3 <= n3) then
        Jp2st=1; J2st=1

        ! Loop over G_y in the small box.
        do j=1,m2eff,lot1
          ma=j
          mb=min(j+(lot1-1),m2eff)
          n1dfft=mb-ma+1

          ! Zero-pad input.
          ! input:  G1,G2,R3,JG2,(Rp3)
          ! output: G2,G1,R3,JG2,(Rp3)
          if (nproc_fft==1) then
            call mpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot1,max1,md1,m1,n1,&
&             md2proc,nd3proc,nproc_fft,ioption,zmpi2(:,:,:,:,idat),zw,max2,m2,n2)
          else
            call mpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot1,max1,md1,m1,n1,&
&             md2proc,nd3proc,nproc_fft,ioption,zmpi1(:,:,:,:,idat),zw,max2,m2,n2)
          end if

          ! Transform along x
          ! input:  G2,G1,R3,(Rp3)
          ! output: G2,R1,R3,(Rp3)
          if (n1dfft == lot1) then
            call dfftw_execute_dft(bw_plan1_lot, zw, zt(1,j,1))
          else
            call dfftw_execute_dft(bw_plan1_rest, zw, zt(1,j,1))
          end if

        end do ! j

        ! Transform along y axis (take into account c2c or c2r case).
        ! Must loop over the full box.
        do j=1,n1eff,lot2
          ma=j
          mb=min(j+(lot2-1),n1eff)
          n1dfft=mb-ma+1
          includelast=1

          if (cplexwf==1) then
            jeff=2*j-1
            if (mb==n1eff .and. n1eff*2/=n1) includelast=0
          end if

          ! Zero-pad the input.
          ! input:  G2,R1,R3,(Rp3)
          ! output: R1,G2,R3,(Rp3)
          if (cplexwf==2) then
            call switch_cent(n1dfft,max2,m2,n2,lot2,n1,lzt,zt(1,1,j),zw)
          else
            call switchreal_cent(includelast,n1dfft,max2,n2,lot2,m1zt,lzt,zt(1,1,jeff),zw)
          end if

          ! input:  R1,G2,R3,(Rp3)
          ! output: R1,R2,R3,(Rp3)
          if (n1dfft == lot2) then
            call dfftw_execute_dft(bw_plan2_lot, zw, zr(1,j,1,j3,idat))
          else
            call dfftw_execute_dft(bw_plan2_rest, zw, zr(1,j,1,j3,idat))
          end if

        end do

        ! Treat real wavefunctions.
        if (cplexwf==1) then
          n1half=n1/2
          ! If odd
          if (n1half*2/=n1) then
            do i2=1,n2
              zr(1,n1,i2,j3,idat)=zr(1,n1eff,i2,j3,idat)
              zr(2,n1,i2,j3,idat)=zero
            end do
          end if
          do i2=1,n2
            do i1=n1half,1,-1
              zr(1,2*i1-1,i2,j3,idat)=zr(1,i1,i2,j3,idat)
              zr(1,2*i1  ,i2,j3,idat)=zr(2,i1,i2,j3,idat)
              zr(2,2*i1-1,i2,j3,idat)=zero
              zr(2,2*i1  ,i2,j3,idat)=zero
            end do
          end do
        end if

      end if

   end do ! j3
 end do ! idat

 call dfftw_destroy_plan(bw_plan3_lot)
 if (mod(m1, lot3) /= 0) then
   call dfftw_destroy_plan(bw_plan3_rest)
 end if

 call dfftw_destroy_plan(bw_plan1_lot)
 if (mod(m2eff, lot1) /= 0) then
   call dfftw_destroy_plan(bw_plan1_rest)
 end if

 call dfftw_destroy_plan(bw_plan2_lot)
 if (mod(n1eff, lot2) /= 0) then
   call dfftw_destroy_plan(bw_plan2_rest)
 end if

 ABI_FREE(zmpi2)
 ABI_FREE(zw)
 ABI_FREE(zt)
 if (nproc_fft>1)  then
   ABI_FREE(zmpi1)
 end if

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/cplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc/))
 ABI_UNUSED((/ max1,max2,max3,m1,m2,m3,md1,md2proc,md3,comm_fft/))
 ABI_UNUSED((/zf(1,1,1,1,1),zr(1,1,1,1,1)/))
#endif

end subroutine fftw3_mpiback_manywf
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_mpiforw_manywf
!! NAME
!!  fftw3_mpiforw_manywf
!!
!! FUNCTION
!!   Does multiple 3-dim backward FFTs from real into Fourier space
!!   Adopt standard convention that isign=-1 for forward transform
!!   CALCULATES THE DISCRETE FOURIERTRANSFORM
!!
!!   ZF(I1,I3,I2)=S_(j1,j2,j3) EXP(isign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) ZR(j1,j2,j3)
!!
!!   in parallel using MPI/OpenMP.
!!
!! INPUT:
!!   ZR: input array
!!        ZR(1,i1,i2,i3,idat)=real(R(i1,i2,i3,idat))
!!        ZR(2,i1,i2,i3,idat)=imag(R(i1,i2,i3,idat))
!!        i1=1,n1 , i2=1,n2 , i3=1,n3 , idat=1,ndat
!!   NOTE that ZR is changed by the routine
!!
!!   n1,n2,n3: logical dimension of the transform. As transform lengths
!!             most products of the prime factors 2,3,5 are allowed.
!!             The detailed table with allowed transform lengths can
!!             be found in subroutine CTRIG
!!   nd1,nd2,nd3: Dimension of ZR
!!   nd3proc=((nd3-1)/nproc_fft)+1  maximal number of big box 3rd dim slices for one proc
!!
!! OUTPUT:
!!   ZF: output array (note the switch of i2 and i3)
!!        real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!!        imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!!   max1 is positive or zero ; m1 >=max1+1
!!     i1= 1... max1+1 corresponds to positive and zero wavevectors 0 ... max1
!!     then, if m1 > max1+1, one has min1=max1-m1+1 and
!!     i1= max1+2 ... m1 corresponds to negative wavevectors min1 ... -1
!!     i2 and i3 have a similar definition of range
!!   idat=1,ndat
!!   md1,md2,md3: Dimension of ZF
!!   md2proc=((md2-1)/nproc_fft)+1  maximal number of small box 2nd dim slices for one proc
!!   nproc_fft: number of processors used as returned by MPI_COMM_SIZE
!!   me_fft: [0:nproc-1] rank of the processor in the FFT communicator.
!!   comm_fft=MPI communicator for parallel FFT.
!!
!! NOTES
!!  The maximum number of processors that can reasonably be used is max(n2/2,n3/2)
!!
!!  It is very important to find the optimal
!!  value of NCACHE. NCACHE determines the size of the work array ZW, that
!!  has to fit into cache. It has therefore to be chosen to equal roughly
!!   half the size of the physical cache in units of real*8 numbers.
!!  The optimal value of ncache can easily be determined by numerical
!!  experimentation. A too large value of ncache leads to a dramatic
!!  and sudden decrease of performance, a too small value to a to a
!!  slow and less dramatic decrease of performance. If NCACHE is set
!!  to a value so small, that not even a single one dimensional transform
!!  can be done in the workarray zw, the program stops with an error message.
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_mpiforw_manywf(cplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc,&
&        max1,max2,max3,m1,m2,m3,md1,md2proc,md3,zr,zf,comm_fft)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc
 integer,intent(in) :: max1,max2,max3,m1,m2,m3,md1,md2proc,md3,comm_fft
!arrays
 real(dp),intent(inout) :: zr(2,nd1,nd2,nd3proc,ndat)
 real(dp),intent(out) :: zf(2,md1,md3,md2proc,ndat)

!Local variables-------------------------------
!scalars
#ifdef HAVE_FFTW3
 integer :: j,i1,i2,i3,idat,ierr,nproc_fft,me_fft
 integer :: ioption,j2,j3,j2st,jp2st,lot1,lot2,lot3,lzt,m1zt,ma,mb,n1dfft,nnd3
 integer :: m2eff,ncache,n1eff,n1half,i1inv,i2inv,i3inv,nthreads
 integer(KIND_FFTW_PLAN) :: fw_plan1_lot,fw_plan1_rest
 integer(KIND_FFTW_PLAN) :: fw_plan2_lot,fw_plan2_rest
 integer(KIND_FFTW_PLAN) :: fw_plan3_lot,fw_plan3_rest
 character(len=500) :: msg
!arrays
 integer :: requests(ndat)
 real(dp) ABI_ASYNC, allocatable :: zmpi1(:,:,:,:,:),zmpi2(:,:,:,:,:) ! work arrays for MPI
 real(dp),allocatable :: zw(:,:),zt(:,:,:) ! cache work array and array for transpositions
! FFT work arrays
 real(dp) :: tsec(2)

! *************************************************************************

 ! FIXME must provide a default value but which one?
 !ioption = 0
 ioption = 1
 !if (paral_kgb==1) ioption=1

 nproc_fft = xmpi_comm_size(comm_fft); me_fft = xmpi_comm_rank(comm_fft)

 ! find cache size that gives optimal performance on machine
 ncache=2*max(n1,n2,n3,1024)
 !ncache=2*max(n1,n2,n3,16*1024)

 if (ncache/(2*max(n1,n2,n3))<1) then
   write(msg,'(5a)') &
&    'ncache has to be enlarged to be able to hold at',ch10, &
&    'least one 1-d FFT of each size even though this will',ch10,&
&    'reduce the performance for shorter transform lengths'
   MSG_ERROR(msg)
 end if

 ! Effective m1 and m2 (complex-to-complex or real-to-complex)
 n1eff=n1; m2eff=m2; m1zt=n1
 if (cplexwf==1) then
   n1eff=(n1+1)/2; m2eff=m2/2+1; m1zt=2*(n1/2+1)
 end if

 lzt=m2eff
 if (mod(m2eff,2)==0) lzt=lzt+1
 if (mod(m2eff,4)==0) lzt=lzt+1

 ! maximal number of big box 3rd dim slices for all procs
 nnd3=nd3proc*nproc_fft

 ABI_MALLOC(zw,(2,ncache/2))
 ABI_MALLOC(zt,(2,lzt,m1zt))
 ABI_MALLOC(zmpi2,(2,md1,md2proc,nnd3,ndat))
 if (nproc_fft>1)  then
   ABI_MALLOC(zmpi1,(2,md1,md2proc,nnd3,ndat))
 end if

 ! Create plans.
 ! The prototype for sfftw_plan_many_dft is:
 ! sfftw_plan_many_dft(rank, n, howmany,
 !   fin,  iembed, istride, idist,
 !   fout, oembed, ostride, odist, isign, my_flags)

 lot2=ncache/(2*n2)
 lot1=ncache/(2*n1)
 lot3=ncache/(2*n3)

 nthreads = xomp_get_num_threads(open_parallel=.TRUE.)
 !nthreads = 1

 fw_plan3_lot = dplan_many_dft_2D(1, [n3], lot3, &
&    zw, [ncache/2], lot3, 1,                          &
&    zw, [ncache/2], lot3, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(m1, lot3) /= 0) then
   fw_plan3_rest = dplan_many_dft_2D(1, [n3], mod(m1, lot3), &
&    zw, [ncache/2], lot3, 1,                                      &
&    zw, [ncache/2], lot3, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 fw_plan1_lot = dplan_many_dft_2D(1, [n1], lot1, &
&    zt, [lzt, m1zt],   lzt,  1,                       &
&    zw, [ncache/2], lot1, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(m2eff, lot1) /= 0) then
   fw_plan1_rest = dplan_many_dft_2D(1, [n1], mod(m2eff, lot1), &
&    zt, [lzt, m1zt],   lzt, 1,                                       &
&    zw, [ncache/2], lot1, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 fw_plan2_lot = dplan_many_dft_2D(1, [n2], lot2, &
&    zr, [nd1,nd2,nd3proc,ndat], nd1, 1,               &
&    zw, [ncache/2], lot2, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(n1eff, lot2) /= 0) then
   fw_plan2_rest = dplan_many_dft_2D(1, [n2], mod(n1eff,lot2), &
&    zr, [nd1,nd2,nd3proc,ndat], nd1, 1,                             &
&    zw, [ncache/2], lot2, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 do idat=1,ndat
   ! Loop over the z-planes treated by this node
   do j3=1,nd3proc

     if (me_fft*nd3proc+j3 <= n3) then
       Jp2st=1
       J2st=1

       ! Treat real wavefunctions.
       if (cplexwf==1) then
         n1half=n1/2
         do i2=1,n2
           do i1=1,n1half
             zr(1,i1,i2,j3,idat)=zr(1,2*i1-1,i2,j3,idat)
             zr(2,i1,i2,j3,idat)=zr(1,2*i1  ,i2,j3,idat)
           end do
         end do
         ! If odd
         if(n1half*2/=n1)then
           do i2=1,n2
             zr(1,n1eff,i2,j3,idat)=zr(1,n1,i2,j3,idat)
             zr(2,n1eff,i2,j3,idat)=zero
           end do
         end if
       end if

       ! transform along y axis
       ! input: R1,R2,R3,(Rp3)
       ! input: R1,G2,R3,(Rp3)
       do j=1,n1eff,lot2
         ma=j
         mb=min(j+(lot2-1),n1eff)
         n1dfft=mb-ma+1

         if (n1dfft == lot2) then
           call dfftw_execute_dft(fw_plan2_lot,  zr(1,j,1,j3,idat), zw)
         else
           call dfftw_execute_dft(fw_plan2_rest, zr(1,j,1,j3,idat), zw)
         end if

         ! input:  R1,G2,R3,(Rp3)
         ! output: G2,R1,R3,(Rp3)
         if (cplexwf==2) then
           call unswitch_cent(n1dfft,max2,m2,n2,lot2,n1,lzt,zw,zt(1,1,j))
         else
           call unswitchreal_cent(n1dfft,max2,n2,lot2,n1,lzt,zw,zt(1,1,2*j-1))
         end if
       end do

       ! transform along x axis
       ! input: G2,R1,R3,(Rp3)
       do j=1,m2eff,lot1
         ma=j
         mb=min(j+(lot1-1),m2eff)
         n1dfft=mb-ma+1

         if (n1dfft == lot1) then
           call dfftw_execute_dft(fw_plan1_lot,  zt(1,j,1), zw)
         else
           call dfftw_execute_dft(fw_plan1_rest, zt(1,j,1), zw)
         end if
         ! output: G2,G1,R3,(Rp3)

         ! input:  G2,G1,R3,Gp2,(Rp3)
         ! output: G1,G2,R3,Gp2,(Rp3)
         if (nproc_fft==1) then
           call unmpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot1,max1,md1,m1,n1,&
&            md2proc,nd3proc,nproc_fft,ioption,zw,zmpi2(:,:,:,:,idat))
         else
           call unmpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot1,max1,md1,m1,n1,&
&            md2proc,nd3proc,nproc_fft,ioption,zw,zmpi1(:,:,:,:,idat))
         end if
       end do
     end if
   end do ! j3

   ! Interprocessor data transposition
   ! input:  G1,G2,R3,Gp2,(Rp3)
   ! output: G1,G2,R3,Rp3,(Gp2)
   if (nproc_fft>1) then
     call timab(544,1,tsec)
     call xmpi_ialltoall(zmpi1(:,:,:,:,idat),2*md1*md2proc*nd3proc, &
&                        zmpi2(:,:,:,:,idat),2*md1*md2proc*nd3proc,comm_fft,requests(idat))
     call timab(544,2,tsec)
   end if
 end do

 do idat=1,ndat
    if (nproc_fft>1) call xmpi_wait(requests(idat),ierr)
   ! transform along z axis
   ! input: G1,G2,R3,(Gp2)

   do j2=1,md2proc
     if (me_fft*md2proc+j2 <= m2eff) then
       ! write(std_out,*)' forwf_wf : before unscramble, j2,md2proc,me_fft,m2=',j2,md2proc,me_fft,m2
       do i1=1,m1,lot3
         ma=i1
         mb=min(i1+(lot3-1),m1)
         n1dfft=mb-ma+1

         ! input:  G1,G2,R3,(Gp2)
         ! output: G1,R3,G2,(Gp2)
         call unscramble(i1,j2,lot3,n1dfft,md1,n3,md2proc,nnd3,zmpi2(:,:,:,:,idat),zw)

         if (n1dfft == lot3) then
           call dfftw_execute_dft(fw_plan3_lot, zw, zw)
         else
           call dfftw_execute_dft(fw_plan3_rest, zw, zw)
         end if

         call unfill_cent(md1,md3,lot3,n1dfft,max3,m3,n3,zw,zf(1,i1,1,j2,idat))
         ! output: G1,G3,G2,(Gp2)
       end do
     end if
   end do

   if (cplexwf==1) then
     ! Complete missing values with complex conjugate
     ! Inverse of ix is located at nx+2-ix , except for ix=1, for which it is 1.
     do i3=1,m3
       i3inv=m3+2-i3
       if(i3==1)i3inv=1

       if (m2eff>1) then
         do i2=2,m2eff
           i2inv=m2+2-i2
           zf(1,1,i3inv,i2inv,idat)= zf(1,1,i3,i2,idat)
           zf(2,1,i3inv,i2inv,idat)=-zf(2,1,i3,i2,idat)
           do i1=2,m1
             i1inv=m1+2-i1
             zf(1,i1inv,i3inv,i2inv,idat)= zf(1,i1,i3,i2,idat)
             zf(2,i1inv,i3inv,i2inv,idat)=-zf(2,i1,i3,i2,idat)
           end do
         end do
       end if
     end do
   end if

 end do ! idat

 call dfftw_destroy_plan(fw_plan3_lot)
 if (mod(m1, lot3) /= 0) then
   call dfftw_destroy_plan(fw_plan3_rest)
 end if

 call dfftw_destroy_plan(fw_plan1_lot)
 if (mod(m2eff, lot1) /= 0) then
   call dfftw_destroy_plan(fw_plan1_rest)
 end if

 call dfftw_destroy_plan(fw_plan2_lot)
 if (mod(n1eff, lot2) /= 0) then
   call dfftw_destroy_plan(fw_plan2_rest)
 end if

 ABI_FREE(zmpi2)
 ABI_FREE(zw)
 ABI_FREE(zt)
 if (nproc_fft>1)  then
   ABI_FREE(zmpi1)
 end if

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/cplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc/))
 ABI_UNUSED((/max1,max2,max3,m1,m2,m3,md1,md2proc,md3,comm_fft/))
 ABI_UNUSED((/zf(1,1,1,1,1),zr(1,1,1,1,1)/))
#endif

end subroutine fftw3_mpiforw_manywf
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_applypot_many
!! NAME
!!  fftw3_applypot_many
!!
!! FUNCTION
!! Applies the local real space potential to multiple wavefunctions in Fourier space
!!
!! INPUTS
!!   ZF: Wavefunction (input/output) (note the switch of i2 and i3)
!!        real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!!        imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!!   max1 is positive or zero ; m1 >=max1+1
!!   i1= 1... max1+1 corresponds to positive and zero wavevectors 0 ... max1
!!   then, if m1 > max1+1, one has min1=max1-m1+1 and
!!   i1= max1+2 ... m1 corresponds to negative wavevectors min1 ... -1
!!   i2 and i3 have a similar definition of range
!!   idat=1,ndat
!!   md1,md2,md3: Dimension of ZF (input as well as output), distributed on different procs
!!   md2proc=((md2-1)/nproc_fft)+1  maximal number of small box 2nd dim slices for one proc
!!
!!   POT: Potential
!!        POT(cplex*i1,i2,i3)
!!        cplex=1 or 2 ,  i1=1,n1 , i2=1,n2 , i3=1,n3
!!   nd1,nd2,nd3: dimension of pot
!!   comm_fft: MPI communicator
!!   nproc_fft: number of processors used as returned by MPI_COMM_SIZE
!!   me_fft: [0:nproc_fft-1] number of processor as returned by MPI_COMM_RANK
!!    n1,n2,n3: logical dimension of the transform. As transform lengths
!!              most products of the prime factors 2,3,5 are allowed.
!!             The detailed table with allowed transform lengths can
!!             be found in subroutine CTRIG
!!
!! NOTES:
!!   PERFORMANCE CONSIDERATIONS:
!!   The maximum number of processors that can reasonably be used is max(n2/2,n3/2)
!!
!!   It is very important to find the optimal
!!   value of NCACHE. NCACHE determines the size of the work array ZW, that
!!   has to fit into cache. It has therefore to be chosen to equal roughly
!!    half the size of the physical cache in units of real*8 numbers.
!!   The optimal value of ncache can easily be determined by numerical
!!   experimentation. A too large value of ncache leads to a dramatic
!!   and sudden decrease of performance, a too small value to a to a
!!   slow and less dramatic decrease of performance. If NCACHE is set
!!   to a value so small, that not even a single one dimensional transform
!!   can be done in the workarray zw, the program stops with an error message.
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE


subroutine fftw3_applypot_many(cplexwf,cplex,ndat,n1,n2,n3,nd1,nd2,nd3,nd3proc,&
&  max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3,&
&  max1o,max2o,max3o,m1o,m2o,m3o,comm_fft,nproc_fft,me_fft,pot,zf)

!Arguments ------------------------------------
 integer,intent(in) :: cplexwf,cplex,ndat,n1,n2,n3,nd1,nd2,nd3,nd3proc
 integer,intent(in) :: max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3
 integer,intent(in) :: max1o,max2o,max3o,m1o,m2o,m3o,comm_fft,nproc_fft,me_fft
 real(dp),intent(in) :: pot(cplex*nd1,nd2,nd3)
 real(dp),intent(inout) :: zf(2,md1,md3,md2proc,ndat)

!Local variables-------------------------------
!scalars
#ifdef HAVE_FFTW3
 integer,parameter :: unused0=0
 integer :: j,i1,i2,i3,idat,ierr,j3glob,nthreads
 integer :: ioption,j2,j3,lzt,m1zt,ma,mb,n1dfft,nnd3,lot1,lot2,lot3
 integer :: m2eff,ncache,n1eff,i1inv,i2inv,i3inv,jeff,includelast,j2stb
 integer :: jx,j2stf,Jp2stb,Jp2stf,m2ieff,m2oeff
 integer(KIND_FFTW_PLAN) :: bw_plan1_lot,bw_plan1_rest
 integer(KIND_FFTW_PLAN) :: bw_plan2_lot,bw_plan2_rest
 integer(KIND_FFTW_PLAN) :: bw_plan3_lot,bw_plan3_rest
 integer(KIND_FFTW_PLAN) :: fw_plan1_lot,fw_plan1_rest
 integer(KIND_FFTW_PLAN) :: fw_plan2_lot,fw_plan2_rest
 integer(KIND_FFTW_PLAN) :: fw_plan3_lot,fw_plan3_rest
 character(len=500) :: msg
!arrays
 integer :: requests(ndat)
 real(dp) :: tsec(2)
 real(dp) ABI_ASYNC, allocatable :: zmpi1(:,:,:,:,:),zmpi2(:,:,:,:,:) ! work arrays for MPI
 real(dp),allocatable :: zw(:,:),zt(:,:,:) ! cache work array and array for transpositions
! FFT work arrays

! *************************************************************************

 !ioption=0 ! This was in the old version.
 ioption=1 ! This one is needed to be compatible with paral_kgb

 ncache=2*max(n1,n2,n3,1024)
 if (ncache/(2*max(n1,n2,n3)) < 1) then
   write(msg,"(5a)") &
&    'ncache has to be enlarged to be able to hold at',ch10,&
&    'least one 1-d FFT of each size even though this will',ch10,&
&    'reduce the performance for shorter transform lengths'
   MSG_ERROR(msg)
 end if

 !call wrtout(std_out,"applypot with non-blocking IALLTOALL + FFTW3","COLL")
 !write(std_out,"(a,i0)")"in applypot_many with ndat: ",ndat

 ! Effective m1 and m2 (complex-to-complex or real-to-complex)
 n1eff=n1; m2ieff=m2i; m2oeff=m2o; m1zt=n1
 if (cplexwf==1) then
   n1eff=(n1+1)/2; m2ieff=m2i/2+1; m2oeff=m2o/2+1; m1zt=2*(n1/2+1)
 end if

 m2eff=max(m2ieff,m2oeff)
 lzt=m2eff
 if (mod(m2eff,2) == 0) lzt=lzt+1
 if (mod(m2eff,4) == 0) lzt=lzt+1

 ! maximal number of big box 3rd dim slices for all procs
 nnd3=nd3proc*nproc_fft

 ABI_ALLOCATE(zw,(2,ncache/2))
 ABI_ALLOCATE(zt,(2,lzt,m1zt))
 ABI_ALLOCATE(zmpi2,(2,md1,md2proc,nnd3,ndat))
 if (nproc_fft > 1)  then
   ABI_ALLOCATE(zmpi1,(2,md1,md2proc,nnd3,ndat))
 end if

 lot3=ncache/(2*n3)
 lot1=ncache/(2*n1)
 lot2=ncache/(2*n2)

 nthreads = xomp_get_num_threads(open_parallel=.TRUE.)
 !nthreads = 1

 ! The prototype for sfftw_plan_many_dft is:
 ! sfftw_plan_many_dft(rank, n, howmany,
 !   fin,  iembed, istride, idist,
 !   fout, oembed, ostride, odist, isign, my_flags)

 ! Create plans for G --> R (see back_wf)
 bw_plan3_lot = dplan_many_dft_2D(1, [n3], lot3, &
&    zw, [ncache/2], lot3, 1,                          &
&    zw, [ncache/2], lot3, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(m1i, lot3) /= 0) then
   bw_plan3_rest = dplan_many_dft_2D(1, [n3], mod(m1i, lot3),&
&      zw, [ncache/2], lot3, 1,                                    &
&      zw, [ncache/2], lot3, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 bw_plan1_lot = dplan_many_dft_2D(1, [n1], lot1, &
&    zw, [ncache/2],  lot1, 1,                         &
&    zt, [lzt, m1zt], lzt, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(m2ieff, lot1) /= 0) then
   bw_plan1_rest = dplan_many_dft_2D(1, [n1], mod(m2ieff, lot1), &
&      zw, [ncache/2],  lot1, 1,                                       &
&      zt, [lzt, m1zt], lzt, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 !TODO this won't work if iclexwf==1
 ! Recheck this
 bw_plan2_lot = dplan_many_dft_2D(1, [n2], lot2, &
&    zw, [ncache/2], lot2, 1,                          &
&    zw, [ncache/2], lot2, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(n1eff, lot2) /= 0) then
   bw_plan2_rest = dplan_many_dft_2D(1, [n2], mod(n1eff,lot2), &
&      zw, [ncache/2], lot2, 1,                                      &
&      zw, [ncache/2], lot2, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 ! Create plans for G --> R (see forw_wf)
 fw_plan3_lot = dplan_many_dft_2D(1, [n3], lot3, &
&    zw, [ncache/2], lot3, 1,                          &
&    zw, [ncache/2], lot3, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(m1o, lot3) /= 0) then
   fw_plan3_rest = dplan_many_dft_2D(1, [n3], mod(m1o, lot3),&
&    zw, [ncache/2], lot3, 1,                                      &
&    zw, [ncache/2], lot3, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 fw_plan1_lot = dplan_many_dft_2D(1, [n1], lot1,&
&    zt, [lzt, m1zt], lzt,  1,                        &
&    zw, [ncache/2],  lot1, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(m2oeff, lot1) /= 0) then
   fw_plan1_rest = dplan_many_dft_2D(1, [n1], mod(m2oeff, lot1),&
&    zt, [lzt, m1zt], lzt,  1,                                        &
&    zw, [ncache/2],  lot1, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 fw_plan2_lot = dplan_many_dft_2D(1, [n2], lot2,&
&    zw, [ncache/2], lot2, 1,                         &
&    zw, [ncache/2], lot2, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)

 if (mod(n1eff, lot2) /= 0) then
   fw_plan2_rest = dplan_many_dft_2D(1, [n2], mod(n1eff,lot2),&
&    zw, [ncache/2], lot2, 1,                                       &
&    zw, [ncache/2], lot2, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)
 end if

 ! Here we take advantage of non-blocking IALLTOALL:
 ! Perform the first step of MPI-FFT for ndat wavefunctions.
 do idat=1,ndat
   !
   ! transform along z axis
   ! input: G1,G3,G2,(Gp2)
   do j2=1,md2proc
     if (me_fft*md2proc+j2 <= m2ieff) then
       do i1=1,m1i,lot3
         ma=i1
         mb=min(i1+(lot3-1),m1i)
         n1dfft=mb-ma+1

         ! zero-pad n1dfft G_z lines
         ! input: G1,G3,G2,(Gp2)
         call fill_cent(md1,md3,lot3,n1dfft,max3i,m3i,n3,zf(1,i1,1,j2,idat),zw)

         if (n1dfft == lot3) then
           call dfftw_execute_dft(bw_plan3_lot, zw, zw)
         else
           call dfftw_execute_dft(bw_plan3_rest, zw, zw)
         end if

         ! Local rotation.
         ! input:  G1,R3,G2,(Gp2)
         ! output: G1,G2,R3,(Gp2)
         call scramble(i1,j2,lot3,n1dfft,md1,n3,md2proc,nnd3,zw,zmpi2(:,:,:,:,idat))
       end do
     end if
   end do

   ! Interprocessor data transposition
   ! input:  G1,G2,R3,Rp3,(Gp2)
   ! output: G1,G2,R3,Gp2,(Rp3)
   if (nproc_fft > 1) then
     call timab(543,1,tsec)
     call xmpi_ialltoall(zmpi2(:,:,:,:,idat),2*md1*md2proc*nd3proc,&
&                        zmpi1(:,:,:,:,idat),2*md1*md2proc*nd3proc,comm_fft,requests(idat))
     call timab(543,2,tsec)
   end if
 end do ! idat

 ! The second step of MPI-FFT
 do idat=1,ndat
    ! Make sure communication is completed.
    if (nproc_fft>1) call xmpi_wait(requests(idat),ierr)

   do j3=1,nd3proc
     j3glob = j3 + me_fft*nd3proc
     if (me_fft*nd3proc+j3 <= n3) then
       Jp2stb=1; J2stb=1
       Jp2stf=1; J2stf=1

       ! transform along x axis
       do j=1,m2ieff,lot1
         ma=j
         mb=min(j+(lot1-1),m2ieff)
         n1dfft=mb-ma+1

         ! Zero-pad input.
         ! input:  G1,G2,R3,G2,(Rp3)
         ! output: G2,G1,R3,G2,(Rp3)
         if (nproc_fft == 1) then
           call mpiswitch_cent(j3,n1dfft,Jp2stb,J2stb,lot1,max1i,md1,m1i,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zmpi2(:,:,:,:,idat),zw, unused0, unused0, unused0)
         else
           call mpiswitch_cent(j3,n1dfft,Jp2stb,J2stb,lot1,max1i,md1,m1i,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zmpi1(:,:,:,:,idat),zw, unused0, unused0, unused0)
         end if

         ! Transform along x
         ! input:  G2,G1,R3,(Rp3)
         ! output: G2,R1,R3,(Rp3)
         if (n1dfft == lot1) then
           call dfftw_execute_dft(bw_plan1_lot, zw, zt(1,j,1))
         else
           call dfftw_execute_dft(bw_plan1_rest, zw, zt(1,j,1))
         end if
       end do

       ! Transform along y axis (take into account c2c or c2r case).
       ! Must loop over the full box.
       !TODO this won't work
       if (cplexwf==1) then
         if(mod(lot2,2).ne.0) lot2=lot2-1 ! needed to introduce jeff
       end if

       do j=1,n1eff,lot2
         ma=j
         mb=min(j+(lot2-1),n1eff)
         n1dfft=mb-ma+1
         jeff=j
         includelast=1

         if (cplexwf==1) then
           jeff=2*j-1
           includelast=1
           if (mb==n1eff .and. n1eff*2/=n1) includelast=0
         end if

         ! Zero-pad the input.
         !  input: G2,R1,R3,(Rp3)
         ! output: R1,G2,R3,(Rp3)
         if (cplexwf==2) then
           call switch_cent(n1dfft,max2i,m2i,n2,lot2,n1,lzt,zt(1,1,jeff),zw)
         else
           call switchreal_cent(includelast,n1dfft,max2i,n2,lot2,m1zt,lzt,zt(1,1,jeff),zw)
         end if

         ! input:  R1,G2,R3,(Rp3)
         ! output: R1,R2,R3,(Rp3)
         ! Be careful here
         if (n1dfft == lot2) then
           call dfftw_execute_dft(bw_plan2_lot, zw, zw)
         else
           call dfftw_execute_dft(bw_plan2_rest, zw, zw)
         end if

         ! Multiply with potential in real space
         jx=cplex*(jeff-1)+1
         call multpot(cplexwf,cplex,includelast,nd1,nd2,n2,lot2,n1dfft,pot(jx,1,j3glob),zw)

         ! TRANSFORM BACK IN FOURIER SPACE
         ! transform along y axis
         ! input: R1,R2,R3,(Rp3)
         if (n1dfft == lot2) then
           call dfftw_execute_dft(fw_plan2_lot,  zw, zw)
         else
           call dfftw_execute_dft(fw_plan2_rest, zw, zw)
         end if

         !  input: R1,G2,R3,(Rp3)
         ! output: G2,R1,R3,(Rp3)
         if (cplexwf==2) then
           call unswitch_cent(n1dfft,max2o,m2o,n2,lot2,n1,lzt,zw,zt(1,1,jeff))
         else
           call unswitchreal_cent(n1dfft,max2o,n2,lot2,n1,lzt,zw,zt(1,1,jeff))
         end if
       end do ! j

       ! transform along x axis
       ! input:  R2,R1,R3,(Rp3)
       ! output: R2,G1,R3,(Rp3)
       do j=1,m2oeff,lot1
         ma=j
         mb=min(j+(lot1-1),m2oeff)
         n1dfft=mb-ma+1

         if (n1dfft == lot1) then
           call dfftw_execute_dft(fw_plan1_lot,  zt(1,j,1), zw)
         else
           call dfftw_execute_dft(fw_plan1_rest, zt(1,j,1), zw)
         end if

         ! input:  G2,G1,R3,Gp2,(Rp3)
         ! output: G1,G2,R3,Gp2,(Rp3)
         if (nproc_fft == 1) then
           call unmpiswitch_cent(j3,n1dfft,Jp2stf,J2stf,lot1,max1o,md1,m1o,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zw,zmpi2(:,:,:,:,idat))
         else
           call unmpiswitch_cent(j3,n1dfft,Jp2stf,J2stf,lot1,max1o,md1,m1o,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zw,zmpi1(:,:,:,:,idat))
         end if
       end do ! j
     end if
   end do

   ! Interprocessor data transposition
   ! input:  G1,G2,R3,Gp2,(Rp3)
   ! output: G1,G2,R3,Rp3,(Gp2)
   if (nproc_fft > 1) then
     call timab(544,1,tsec)
     call xmpi_ialltoall(zmpi1(:,:,:,:,idat),2*md1*md2proc*nd3proc, &
&                        zmpi2(:,:,:,:,idat),2*md1*md2proc*nd3proc,comm_fft,requests(idat))
     call timab(544,2,tsec)
   end if
 end do

 do idat=1,ndat
   if (nproc_fft>1) call xmpi_wait(requests(idat),ierr)
   ! transform along z axis
   ! input: G1,G2,R3,(Gp2)
   !lot=ncache/(4*n3)
   do j2=1,md2proc
     if (me_fft*md2proc+j2 <= m2oeff) then
       do i1=1,m1o,lot3
         ma=i1
         mb=min(i1+(lot3-1),m1o)
         n1dfft=mb-ma+1

         ! input:  G1,G2,R3,(Gp2)
         ! output: G1,R3,G2,(Gp2)
         call unscramble(i1,j2,lot3,n1dfft,md1,n3,md2proc,nnd3,zmpi2(:,:,:,:,idat),zw)

          if (n1dfft == lot3) then
            call dfftw_execute_dft(fw_plan3_lot, zw, zw)
          else
            call dfftw_execute_dft(fw_plan3_rest, zw, zw)
          end if

         call unfill_cent(md1,md3,lot3,n1dfft,max3o,m3o,n3,zw,zf(1,i1,1,j2,idat))
         ! output: G1,G3,G2,(Gp2)
       end do
     end if
   end do

   ! Complete missing values with complex conjugate
   ! Inverse of ix is located at nx+2-ix , except for ix=1, for which it is 1.
   if (cplexwf==1) then
     do i3=1,m3o
       i3inv=m3o+2-i3
       if (i3==1) i3inv=1
       if (m2oeff>1)then
         do i2=2,m2oeff
           i2inv=m2o+2-i2
           zf(1,1,i3inv,i2inv,idat)= zf(1,1,i3,i2,idat)
           zf(2,1,i3inv,i2inv,idat)=-zf(2,1,i3,i2,idat)
           do i1=2,m1o
             i1inv=m1o+2-i1
             zf(1,i1inv,i3inv,i2inv,idat)= zf(1,i1,i3,i2,idat)
             zf(2,i1inv,i3inv,i2inv,idat)=-zf(2,i1,i3,i2,idat)
           end do
         end do
       end if
     end do
   end if

 end do ! idat

 call dfftw_destroy_plan(bw_plan3_lot)
 if (mod(m1i, lot3) /= 0) then
   call dfftw_destroy_plan(bw_plan3_rest)
 end if

 call dfftw_destroy_plan(bw_plan1_lot)
 if (mod(m2ieff, lot1) /= 0) then
   call dfftw_destroy_plan(bw_plan1_rest)
 end if

 call dfftw_destroy_plan(bw_plan2_lot)
 if (mod(n1eff, lot2) /= 0) then
   call dfftw_destroy_plan(bw_plan2_rest)
 end if

 call dfftw_destroy_plan(fw_plan3_lot)
 if (mod(m1o, lot3) /= 0) then
   call dfftw_destroy_plan(fw_plan3_rest)
 end if

 call dfftw_destroy_plan(fw_plan1_lot)
 if (mod(m2oeff, lot1) /= 0) then
   call dfftw_destroy_plan(fw_plan1_rest)
 end if

 call dfftw_destroy_plan(fw_plan2_lot)
 if (mod(n1eff, lot2) /= 0) then
   call dfftw_destroy_plan(fw_plan2_rest)
 end if

 ABI_DEALLOCATE(zmpi2)
 ABI_DEALLOCATE(zw)
 ABI_DEALLOCATE(zt)
 if (nproc_fft > 1)  then
   ABI_DEALLOCATE(zmpi1)
 end if

#else
 MSG_ERROR("FFTW3 support not activated")
 ABI_UNUSED((/cplexwf,cplex,ndat,n1,n2,n3,nd1,nd2,nd3,nd3proc/))
 ABI_UNUSED((/max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3/))
 ABI_UNUSED((/max1o,max2o,max3o,m1o,m2o,m3o,comm_fft,nproc_fft,me_fft/))
 ABI_UNUSED((/pot(1,1,1),zf(1,1,1,1,1)/))
#endif

end subroutine fftw3_applypot_many
!!***

!----------------------------------------------------------------------

!!****f* m_fftw3/fftw3_poisson
!! NAME
!! fftw3_poisson
!!
!! FUNCTION
!!  Solve the Poisson equation in G-space given the density, n(r),
!!  in real space of the FFT box.
!!
!! INPUTS
!! cplex=1 if fofr is real, 2 if fofr is complex
!! nx,ny,nz=Number of FFT points along the three directions.
!! ldx,ldy,ldz=Leading dimension of the array nr and vg.
!! ndat = Number of densities
!! vg(nx*ny*nz)=Potential in reciprocal space.
!!
!! SIDE EFFECTS
!! nr(cplex*ldx*ldy*ldz*ndat)
!!    input: n(r) (real or complex)
!!    output: the hartree potential in real space
!!
!! NOTES
!!   vg is given on the FFT mesh instead of the augmented mesh [ldx,ldy,ldz]
!!   in order to simplify the interface with the other routines operating of vg
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      fftw3_destroy_plan,fftw3_execute_dft
!!
!! SOURCE

subroutine fftw3_poisson(cplex,nx,ny,nz,ldx,ldy,ldz,ndat,vg,nr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nx,ny,nz,ldx,ldy,ldz,ndat
!arrays
 real(dp),intent(inout) :: nr(cplex*ldx*ldy*ldz*ndat)
 real(dp),intent(in) :: vg(nx*ny*nz)

#ifdef HAVE_FFTW3
!Local variables-------------------------------
!scalars
 integer,parameter :: rank1=1,rank2=2
 integer :: ii,jj,kk,sidx,ig,ir,vgbase,ypad
 integer, parameter :: nthreads=1
 integer(KIND_FFTW_PLAN) :: bw_plan_xy,bw_plan3
 integer(KIND_FFTW_PLAN) :: fw_plan_xy,fw_plan3
 real(dp) :: fft_fact,vg_fftfact

! *************************************************************************

 !write(std_out,*)"in poisson"
 ABI_CHECK(cplex==2,"cplex!=2 not coded")
 ABI_CHECK(ndat==1,"ndat!=1 not coded")

 fft_fact = one/(nx*ny*nz)

 ! The prototype for sfftw_plan_many_dft is:
 ! sfftw_plan_many_dft(n, howmany,
 !   fin,  iembed, istride, idist,
 !   fout, oembed, ostride, odist, isign, my_flags)

 ! 1) ldx*ldy transforms along Rz.
 fw_plan3 = fftw3_plan_many_dft(rank1, (/nz/), ldx*ldy, & ! We have to visit the entire augmented x-y plane!
&   nr, (/ldx, ldy, ldz/), ldx*ldy, 1,                  &
&   nr, (/ldx, ldy, ldz/), ldx*ldy, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)

 call fftw3_execute_dft(fw_plan3, nr, nr) ! Now we have nr(x,y,Gz)
 call fftw3_destroy_plan(fw_plan3)

 ! R --> G Transforms in x-y plane
 fw_plan_xy = fftw3_plan_many_dft(rank2, [nx,ny], 1, &
&     nr, (/ldx, ldy, ldz/), 1, 1,                   &
&     nr, (/ldx, ldy, ldz/), 1, 1, ABI_FFTW_FORWARD, ABI_FFTW_ESTIMATE, nthreads)

 ! G --> R Transforms in x-y plane
 bw_plan_xy = fftw3_plan_many_dft(rank2, [nx, ny], 1, &
&     nr, (/ldx, ldy, ldz/), 1, 1,                    &
&     nr, (/ldx, ldy, ldz/), 1, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 ! Loop on z-planes.
 do kk=1,nz
   sidx = 1 + cplex*(kk-1)*ldx*ldy  !+ cplex*(dat-1) * ldx*ldy*ldz

   call fftw3_execute_dft(fw_plan_xy, nr(sidx:), nr(sidx:))

   ! At this point we have nr(Gx,Gy,Gz) on the current plane.
   ! Multiply by vc(Gx,Gy,Gz) and then back transform immediately to get vc(x,y,Gz)
   ! Note that nr is complex whereas vg is real.
   ! Besides, FFTW returns not normalized FTs if sign=-1 so we have to scale by fft_fact
   vgbase = (kk-1)*nx*ny !;vgbase = (kk-1)*ldx*ldy

   ig = 0
   do jj=1,ny
     ypad = cplex*(jj-1)*ldx + sidx
     do ii=1,nx
       ig = ig + 1
       vg_fftfact = vg(vgbase+ig) * fft_fact

       ir = cplex*(ii-1) + ypad
       nr(ir:ir+1) = nr(ir:ir+1) * vg_fftfact
     end do
   end do

   call fftw3_execute_dft(bw_plan_xy, nr(sidx:), nr(sidx:))
 end do

 ! Free plans
 call fftw3_destroy_plan(fw_plan_xy)
 call fftw3_destroy_plan(bw_plan_xy)

 ! Final transforms of vc(x,y,Gz) along Gz to get vc(x,y,z)
 bw_plan3 = fftw3_plan_many_dft(rank1, (/nz/), ldx*ldy, & ! We have to visit the entire augmented x-y plane!
&   nr, (/ldx, ldy, ldz/), ldx*ldy, 1,                  &
&   nr, (/ldx, ldy, ldz/), ldx*ldy, 1, ABI_FFTW_BACKWARD, ABI_FFTW_ESTIMATE, nthreads)

 call fftw3_execute_dft(bw_plan3, nr, nr)
 call fftw3_destroy_plan(bw_plan3)

#else
 ABI_UNUSED((/cplex,nx,ny,nz,ldx,ldy,ldz,ndat/))
 ABI_UNUSED((/nr(1),vg(1)/))
#endif

end subroutine fftw3_poisson
!!**

!----------------------------------------------------------------------

END MODULE m_fftw3
!!***
