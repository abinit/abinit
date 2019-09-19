!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dfti
!! NAME
!! m_dfti
!!
!! FUNCTION
!!  This module provides wrappers for the MKL DFTI routines: in-place and out-of-place version.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2019 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  1) MPI parallelism is not supported
!!  2) For better performance the FFT divisions should contain small factors  (/2, 3, 5, 7, 11, 13/)
!!     see http://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_userguide_lnx/index.htm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

#ifdef HAVE_DFTI

! Include and generate MKL_DFTI module
#include "mkl_dfti.f90"

! Macros for template files.
#define FFTLIB "DFTI"
#define FFT_PREF(name) CONCAT(dfti_,name)
#define SPAWN_THREADS_HERE(ndat, nthreads) dfti_spawn_threads_here(ndat, nthreads)

#define FFT_DOUBLE 1
#define FFT_SINGLE 2
#define FFT_MIXPREC 3

#endif

MODULE m_dfti

 use defs_basis
 use m_abicore
 use m_errors
 use m_xomp
 use m_cgtools
 use m_cplxtools
 use m_fftcore
 use iso_c_binding
#ifdef HAVE_DFTI
 use MKL_DFTI
#endif

 use m_fstrings,  only : basename, strcat, int2char10
 use m_hide_blas, only : xcopy
 use m_fft_mesh,  only : zpad_t, zpad_init, zpad_free

 implicit none

 private

! Entry points for client code
 public :: dfti_seqfourdp      ! 3D FFT of lengths nx, ny, nz. Mainly used for densities or potentials.
 public :: dfti_seqfourwf      ! FFT transform of wavefunctions (high-level interface).
 public :: dfti_fftrisc
 public :: dfti_fftug          ! G-->R, 3D zero-padded FFT of lengths nx, ny, nz. Mainly used for wavefunctions
 public :: dfti_fftur          ! R-->G, 3D zero-padded FFT of lengths nx, ny, nz. Mainly used for wavefunctions

! Low-level routines.
 public :: dfti_r2c_op         ! Real to complex transform (out-of-place version).
 public :: dfti_c2r_op         ! Complex to real transform (out-of-place version).
 public :: dfti_c2c_op         ! complex to complex transform (out-of-place version).
 public :: dfti_c2c_ip         ! complex to complex transform (in-place version).
 public :: dfti_many_dft_op    ! Driver routine for many out-of-place 3D complex-to-complex FFTs.
 public :: dfti_many_dft_ip    ! Driver routine for many in-place 3D complex-to-complex FFTs.
 public :: dfti_fftpad         ! Driver routines for zero-padded FFT of wavefunctions.

 !FIXME I don't know why gcc does not recognize this one
 ! Perhaps I have to provide an interfaces for fofr(:,:,:,:)
 public :: dfti_fftpad_dp      ! Driver routines for zero-padded FFT of wavefunctions.
 public :: dfti_fftug_dp       ! Driver routines for zero-padded FFT of wavefunctions.
 public :: dfti_use_lib_threads
!!***

 interface dfti_fftrisc
   module procedure dfti_fftrisc_sp
   module procedure dfti_fftrisc_dp
 end interface dfti_fftrisc

 interface dfti_fftug
   module procedure dfti_fftug_dp
   module procedure dfti_fftug_spc
   module procedure dfti_fftug_dpc
 end interface dfti_fftug

 interface dfti_fftur
   module procedure dfti_fftur_dp
   module procedure dfti_fftur_spc
   module procedure dfti_fftur_dpc
 end interface dfti_fftur

 interface dfti_r2c_op
   module procedure dfti_r2c_op_dp
   module procedure dfti_r2c_op_dpc
 end interface dfti_r2c_op

 interface dfti_c2r_op
   module procedure dfti_c2r_op_dp
   module procedure dfti_c2r_op_dpc
 end interface dfti_c2r_op

 interface dfti_c2c_op
   module procedure dfti_c2c_op_spc
   module procedure dfti_c2c_op_dpc
 end interface dfti_c2c_op

 interface dfti_c2c_ip
   module procedure dfti_c2c_ip_spc
   module procedure dfti_c2c_ip_dpc
 end interface dfti_c2c_ip

 !interface dfti_many_dft_op
 !  module procedure dfti_many_dft_op
 !  module procedure dfti_many_dft_op
 !end interface dfti_many_dft_op

 !interface dfti_many_dft_ip
 !  module procedure dfti_many_dft_ip
 !  module procedure dfti_many_dft_ip
 !end interface dfti_many_dft_ip

 interface dfti_fftpad
   module procedure dfti_fftpad_dp
   module procedure dfti_fftpad_spc
   module procedure dfti_fftpad_dpc
 end interface dfti_fftpad

 logical,private,save :: USE_LIB_THREADS = .FALSE.

#ifdef HAVE_DFTI
 ! dfti_alloc_* allocates arrays aligned on DFTI_DEFAULT_ALIGNMENT boundaries.
 integer(C_INT),private,parameter :: DFTI_DEFAULT_ALIGNMENT_SP = 64
 integer(C_INT),private,parameter :: DFTI_DEFAULT_ALIGNMENT_DP = 64

 interface dfti_alloc_real
   !module procedure dfti_alloc_real_sp
   module procedure dfti_alloc_real_dp
 end interface dfti_alloc_real

 interface dfti_alloc_complex
   module procedure dfti_alloc_complex_spc
   module procedure dfti_alloc_complex_dpc
 end interface dfti_alloc_complex

 ! Fortran binding for MKL_malloc
 interface mkl_malloc
   type(C_PTR) function mkl_malloc(alloc_size, alignment) bind(C, name='MKL_malloc')
     import
     integer(C_SIZE_T), value :: alloc_size
     integer(C_INT), value :: alignment
   end function mkl_malloc
 end interface mkl_malloc

 ! Fortran binding for MKL_free
 interface dfti_free
   subroutine mkl_free(cptr) bind(C, name='MKL_free')
     import
     type(C_PTR), value :: cptr
   end subroutine mkl_free
 end interface dfti_free
#endif

!----------------------------------------------------------------------

CONTAINS  !===========================================================
!!***

!!****f* m_dfti/dfti_seqfourdp
!! NAME
!!  dfti_seqfourdp
!!
!! FUNCTION
!! Driver routine for 3D FFT of lengths nx, ny, nz. Mainly used for densities or potentials.
!! FFT Transform is out-of-place
!!
!! INPUTS
!! cplex=1 if fofr is real, 2 if fofr is complex
!! nx,ny,nz=Number of point along the three directions.
!! ldx,ldy,ldz=Leading dimensions of the array.
!! ndat = Number of FFTS
!! isign= +1 : fofg(G) => fofr(R);
!!        -1 : fofr(R) => fofg(G)
!! fofg(2,ldx*ldy*ldz*ndat)=The array to be transformed.
!!
!! OUTPUT
!! fofr(cplex,ldx*ldy*ldz*ndat)=The FFT of fofg
!!
!! PARENTS
!!      fourdp,m_fft
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_seqfourdp(cplex,nx,ny,nz,ldx,ldy,ldz,ndat,isign,fofg,fofr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nx,ny,nz,ldx,ldy,ldz,ndat,isign
!arrays
 real(dp),intent(inout) :: fofg(2*ldx*ldy*ldz*ndat)
 real(dp),intent(inout) :: fofr(cplex*ldx*ldy*ldz*ndat)

!Local variables-------------------------------
!scalars
 integer :: ii,jj
 complex(spc), allocatable :: work_sp(:)

! *************************************************************************

 select case (cplex)
 case (2)
   ! Complex to Complex.
   if (fftcore_mixprec == 1) then
     ! Mixed precision: copyin + in-place + copyout
     ABI_MALLOC(work_sp, (ldx*ldy*ldz*ndat))
     if (isign == +1) then
       work_sp(:) = cmplx(fofg(1::2), fofg(2::2), kind=spc)
     else if (isign == -1) then
       work_sp(:) = cmplx(fofr(1::2), fofr(2::2), kind=spc)
     else
       MSG_BUG("Wrong isign")
     end if

     call dfti_c2c_ip_spc(nx,ny,nz,ldx,ldy,ldz,ndat,isign,work_sp)

     if (isign == +1) then
       jj = 1
       do ii=1,ldx*ldy*ldz*ndat
         fofr(jj) = real(work_sp(ii), kind=dp)
         fofr(jj+1) = aimag(work_sp(ii))
         jj = jj + 2
       end do
     else if (isign == -1) then
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
     case (+1)
       call dfti_many_dft_op(nx,ny,nz,ldx,ldy,ldz,ndat,isign,fofg,fofr)
     case (-1) ! -1
       call dfti_many_dft_op(nx,ny,nz,ldx,ldy,ldz,ndat,isign,fofr,fofg)
     case default
       MSG_BUG("Wrong isign")
     end select
   end if

 case (1) ! Real case.

   select case (isign)
   case (+1) ! G --> R
     call dfti_c2r_op(nx,ny,nz,ldx,ldy,ldz,ndat,fofg,fofr)
   case (-1) ! R --> G
     call dfti_r2c_op(nx,ny,nz,ldx,ldy,ldz,ndat,fofr,fofg)
   case default
     MSG_BUG("Wrong isign")
   end select

 case default
   MSG_BUG("Wrong value for cplex")
 end select

end subroutine dfti_seqfourdp
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_seqfourwf
!! NAME
!! dfti_seqfourwf
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

!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! for option==0, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
!!                fofr(2,ldx,ldy,ldz) contains the output Fourier Transform of fofgin;
!!                no use of denpot, fofgout and npwout.
!! for option==1, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
!!                denpot(cplex*ldx,ldy,ldz) contains the input density at input,
!!                and the updated density at output (accumulated);
!!                no use of fofgout and npwout.
!! for option==2, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
!!                denpot(cplex*ldx,ldy,ldz) contains the input local potential;
!!                fofgout(2,npwout*ndat) contains the output function;
!! for option==3, fofr(2,ldx,ldy,ldz*ndat) contains the input real space wavefunction;
!!                fofgout(2,npwout*ndat) contains its output Fourier transform;
!!                no use of fofgin and npwin.
!!
!! PARENTS
!!      fourwf
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_seqfourwf(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&
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
 integer,parameter :: ndat1=1
 integer :: nx,ny,nz,fftalg,fftalga,fftalgc,fftcache,dat,ptg,ptr,ptgin,ptgout,nthreads
 logical :: use_fftrisc
 character(len=500) :: msg

! *************************************************************************

 if (option==1 .and. cplex/=1) then
   write(msg,'(a,i0)')' With the option number 1, cplex must be 1 but it is cplex=',cplex
   MSG_ERROR(msg)
 end if

 if (option==2 .and. (cplex/=1 .and. cplex/=2)) then
   write(msg,'(a,i0)')' With the option number 2, cplex must be 1 or 2, but it is cplex=',cplex
   MSG_ERROR(msg)
 end if

 nx=ngfft(1); ny=ngfft(2); nz=ngfft(3)
 fftalg=ngfft(7); fftalga=fftalg/100; fftalgc=MOD(fftalg,10)
 fftcache=ngfft(8)

 use_fftrisc = (fftalgc==2)
 if (istwf_k==2.and.option==3) use_fftrisc = .FALSE.
 if (istwf_k>2.and.ANY(option==(/0,3/))) use_fftrisc = .FALSE.

 nthreads = xomp_get_num_threads(open_parallel=.TRUE.)

 if (use_fftrisc) then
   !call wrtout(std_out, " calls dfti_fftrisc")
   if (ndat == 1) then
     if (fftcore_mixprec == 0) then
       call dfti_fftrisc_dp(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
         mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
     else
       call dfti_fftrisc_mixprec(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
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
       if (.not.dfti_spawn_threads_here(ndat,nthreads)) then
         do dat=1,ndat
           ptg = 1 + (dat-1)*npwin
           ptr = 1 + (dat-1)*ldx*ldy*ldz
           call dfti_fftrisc_dp(cplex,denpot,fofgin(1,ptg),fofgout,fofr(1,ptr),gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
&            mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
         end do
       else
!$OMP PARALLEL DO PRIVATE(ptg,ptr)
         do dat=1,ndat
           ptg = 1 + (dat-1)*npwin
           ptr = 1 + (dat-1)*ldx*ldy*ldz
           call dfti_fftrisc_dp(cplex,denpot,fofgin(1,ptg),fofgout,fofr(1,ptr),gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
&            mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
         end do
       end if

     CASE (1)
       !fofgin -> local ur and accumulate density in denpot
       ! TODO this is delicate part to do in parallel, as one should OMP reduce denpot.
       do dat=1,ndat
         ptg = 1 + (dat-1)*npwin
         ptr = 1 + (dat-1)*ldx*ldy*ldz
         call dfti_fftrisc_dp(cplex,denpot,fofgin(1,ptg),fofgout,fofr,gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
&          mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
       end do

     CASE (2)
       ! <G|vloc(r)|fofgin(r)> in fofgout
       if (.not.dfti_spawn_threads_here(ndat,nthreads)) then
         do dat=1,ndat
           ptgin  = 1 + (dat-1)*npwin
           ptgout = 1 + (dat-1)*npwout
           if (fftcore_mixprec == 0) then
             call dfti_fftrisc_dp(cplex,denpot,fofgin(1,ptgin),fofgout(1,ptgout),fofr,gboundin,gboundout,istwf_k,&
               kg_kin,kg_kout,mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
           else
             call dfti_fftrisc_mixprec(cplex,denpot,fofgin(1,ptgin),fofgout(1,ptgout),fofr,gboundin,gboundout,istwf_k,&
               kg_kin,kg_kout,mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
           end if
         end do
       else
!$OMP PARALLEL DO PRIVATE(ptgin,ptgout)
         do dat=1,ndat
           ptgin  = 1 + (dat-1)*npwin
           ptgout = 1 + (dat-1)*npwout
           call dfti_fftrisc_dp(cplex,denpot,fofgin(1,ptgin),fofgout(1,ptgout),fofr,gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
&            mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
         end do
       end if

     CASE (3)
       !fofr -> fofgout
       if (.not.dfti_spawn_threads_here(ndat,nthreads)) then
         do dat=1,ndat
           ptr    = 1 + (dat-1)*ldx*ldy*ldz
           ptgout = 1 + (dat-1)*npwout
           call dfti_fftrisc_dp(cplex,denpot,fofgin,fofgout(1,ptgout),fofr(1,ptr),gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
&            mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
         end do
       else
!$OMP PARALLEL DO PRIVATE(ptr,ptgout)
         do dat=1,ndat
           ptr    = 1 + (dat-1)*ldx*ldy*ldz
           ptgout = 1 + (dat-1)*npwout
           call dfti_fftrisc_dp(cplex,denpot,fofgin,fofgout(1,ptgout),fofr(1,ptr),gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
&            mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)
         end do
       end if

     CASE DEFAULT
       write(msg,'(a,i0,a)')'Option',option,' is not allowed. Only option=0, 1, 2 or 3 are allowed presently.'
       MSG_ERROR(msg)
     END SELECT

   end if

 else
   SELECT CASE (option)
   CASE (0)
     !
     ! FFT u(g) --> u(r)
     if (.not.dfti_spawn_threads_here(ndat,nthreads)) then
       call dfti_fftug_dp(fftalg,fftcache,npwin,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_kin,gboundin,fofgin,fofr)
     else
!$OMP PARALLEL DO PRIVATE(ptg, ptr)
       do dat=1,ndat
         ptg = 1 + (dat-1)*npwin
         ptr = 1 + (dat-1)*ldx*ldy*ldz
         call dfti_fftug_dp(fftalg,fftcache,npwin,nx,ny,nz,ldx,ldy,ldz,ndat1,&
&          istwf_k,mgfft,kg_kin,gboundin,fofgin(1,ptg),fofr(1,ptr))
       end do
     end if

   CASE (1)
     ! TODO this is delicate part to do in parallel, as one should OMP reduce denpot.
     call dfti_fftug_dp(fftalg,fftcache,npwin,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_kin,gboundin,fofgin,fofr)
     call cg_addtorho(nx,ny,nz,ldx,ldy,ldz,ndat,weight_r,weight_i,fofr,denpot)

   CASE (2)

     if (.not.dfti_spawn_threads_here(ndat,nthreads)) then
       call dfti_fftug_dp(fftalg,fftcache,npwin,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_kin,gboundin,fofgin,fofr)
       call cg_vlocpsi(nx,ny,nz,ldx,ldy,ldz,ndat,cplex,denpot,fofr)

       !  The data for option==2 is now in fofr.
       call dfti_fftpad_dp(fofr,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,-1,gboundout)

       call cg_box2gsph(nx,ny,nz,ldx,ldy,ldz,ndat,npwout,kg_kout,fofr,fofgout)
     else

!$OMP PARALLEL DO PRIVATE(ptg, ptr)
       do dat=1,ndat
         ptg = 1 + (dat-1)*npwin
         ptr = 1 + (dat-1)*ldx*ldy*ldz
         call dfti_fftug_dp(fftalg,fftcache,npwin,nx,ny,nz,ldx,ldy,ldz,ndat1,&
&          istwf_k,mgfft,kg_kin,gboundin,fofgin(1,ptg),fofr(1,ptr))

         call cg_vlocpsi(nx,ny,nz,ldx,ldy,ldz,ndat1,cplex,denpot,fofr(1,ptr))

         !  The data for option==2 is now in fofr.
         call dfti_fftpad_dp(fofr(1,ptr),nx,ny,nz,ldx,ldy,ldz,ndat1,mgfft,-1,gboundout)

         ptg = 1 + (dat-1)*npwout
         call cg_box2gsph(nx,ny,nz,ldx,ldy,ldz,ndat1,npwout,kg_kout,fofr(1,ptr),fofgout(1,ptg))
       end do
     end if

   CASE (3)
     !  The data for option==3 is already in fofr.
     if (.not.dfti_spawn_threads_here(ndat,nthreads)) then
       call dfti_fftpad_dp(fofr,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,-1,gboundout)
       call cg_box2gsph(nx,ny,nz,ldx,ldy,ldz,ndat,npwout,kg_kout,fofr,fofgout)
     else
!$OMP PARALLEL DO PRIVATE(ptg, ptr)
       do dat=1,ndat
         ptg = 1 + (dat-1)*npwout
         ptr = 1 + (dat-1)*ldx*ldy*ldz
         call dfti_fftpad_dp(fofr(1,ptr),nx,ny,nz,ldx,ldy,ldz,ndat1,mgfft,-1,gboundout)
         call cg_box2gsph(nx,ny,nz,ldx,ldy,ldz,ndat1,npwout,kg_kout,fofr(1,ptr),fofgout(1,ptg))
       end do
     end if

   CASE DEFAULT
     write(msg,'(a,i0,a)')'Option',option,' is not allowed. Only option=0, 1, 2 or 3 are allowed presently.'
     MSG_ERROR(msg)
   END SELECT
 end if

end subroutine dfti_seqfourwf
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_fftrisc_sp
!! NAME
!! dfti_fftrisc_sp
!!
!! FUNCTION
!! Carry out Fourier transforms between real and reciprocal (G) space,
!! for wavefunctions, contained in a sphere in reciprocal space,
!! in both directions. Also accomplish some post-processing.
!! See dfti_fftrisc_dp for API doc.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_fftrisc_sp(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
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
 real(sp),intent(inout) :: fofgout(2,npwout)    !vz_i

! *************************************************************************

#ifdef HAVE_DFTI

#undef FFT_PRECISION
#undef MYKIND
#undef MYCZERO
#undef MYCMPLX
#undef MYCONJG

#define FFT_PRECISION DFTI_SINGLE
#define MYKIND SPC
#define MYCZERO (0._sp,0._sp)
#define MYCMPLX  CMPLX
#define MYCONJG  CONJG

#include "dfti_fftrisc.finc"

#else
 MSG_ERROR("DFTI support not activated")
 ABI_UNUSED((/cplex,gboundin(1,1),gboundout(1,1),istwf_k,kg_kin(1,1),kg_kout(1,1)/))
 ABI_UNUSED((/mgfft,ngfft(1),npwin,npwout,ldx,ldy,ldz,option/))
 ABI_UNUSED((/denpot(1,1,1),weight_r,weight_i/))
 ABI_UNUSED((/fofgin(1,1),fofgout(1,1),fofr(1,1)/))
#endif

end subroutine dfti_fftrisc_sp
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_fftrisc_dp
!! NAME
!! dfti_fftrisc_dp
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
!!      m_dfti
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_fftrisc_dp(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
& mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,istwf_k,mgfft,ldx,ldy,ldz,npwin,npwout,option
 real(dp),intent(in) :: weight_i,weight_r
!arrays
 integer,intent(in) :: gboundin(2*mgfft+8,2),gboundout(2*mgfft+8,2)
 integer,intent(in) :: kg_kin(3,npwin),kg_kout(3,npwout),ngfft(18)
 real(dp),intent(in) :: fofgin(2,npwin)
 real(dp),intent(inout) :: denpot(cplex*ldx,ldy,ldz)
 real(dp),intent(inout) :: fofr(2,ldx*ldy*ldz)
 real(dp),intent(inout) :: fofgout(2,npwout)    !vz_i

! *************************************************************************

#ifdef HAVE_DFTI

#undef  FFT_PRECISION
#undef  MYKIND
#undef  MYCZERO
#undef  MYCMPLX
#undef  MYCONJG

#define FFT_PRECISION DFTI_DOUBLE
#define MYKIND DPC
#define MYCZERO (0._dp,0._dp)
#define MYCMPLX  DCMPLX
#define MYCONJG  DCONJG

#include "dfti_fftrisc.finc"

#else
 MSG_ERROR("DFTI support not activated")
 ABI_UNUSED((/cplex,gboundin(1,1),gboundout(1,1),istwf_k,kg_kin(1,1),kg_kout(1,1)/))
 ABI_UNUSED((/mgfft,ngfft(1),npwin,npwout,ldx,ldy,ldz,option/))
 ABI_UNUSED((/denpot(1,1,1),fofgin(1,1),fofgout(1,1),fofr(1,1),weight_r,weight_i/))
#endif

end subroutine dfti_fftrisc_dp
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_fftrisc_mixprec
!! NAME
!! dfti_fftrisc_mixprec
!!
!! FUNCTION
!! Carry out Fourier transforms between real and reciprocal (G) space,
!! for wavefunctions, contained in a sphere in reciprocal space,
!! in both directions. Also accomplish some post-processing.
!! This is the Mixed Precision version (dp in input, FFT done with sp data, output is dp)
!! See dfti_fftrisc_dp for API doc.
!!
!! PARENTS
!!      m_dfti
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_fftrisc_mixprec(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,kg_kin,kg_kout,&
& mgfft,ngfft,npwin,npwout,ldx,ldy,ldz,option,weight_r,weight_i)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,istwf_k,mgfft,ldx,ldy,ldz,npwin,npwout,option
 real(dp),intent(in) :: weight_i,weight_r
!arrays
 integer,intent(in) :: gboundin(2*mgfft+8,2),gboundout(2*mgfft+8,2)
 integer,intent(in) :: kg_kin(3,npwin),kg_kout(3,npwout),ngfft(18)
 real(dp),intent(in) :: fofgin(2,npwin)
 real(dp),intent(inout) :: denpot(cplex*ldx,ldy,ldz)
 real(dp),intent(inout) :: fofr(2,ldx*ldy*ldz)
 real(dp),intent(inout) :: fofgout(2,npwout)    !vz_i

! *************************************************************************

#ifdef HAVE_DFTI

#undef  FFT_PRECISION
#undef  MYKIND
#undef  MYCZERO
#undef  MYCMPLX
#undef  MYCONJG

#define FFT_PRECISION DFTI_SINGLE
#define MYKIND SPC
#define MYCZERO (0._sp,0._sp)
#define MYCMPLX  CMPLX
#define MYCONJG  CONJG

#define HAVE_DFTI_MIXED_PRECISION 1

#include "dfti_fftrisc.finc"

#undef HAVE_DFTI_MIXED_PRECISION

#else
 MSG_ERROR("DFTI support not activated")
 ABI_UNUSED((/cplex,gboundin(1,1),gboundout(1,1),istwf_k,kg_kin(1,1),kg_kout(1,1)/))
 ABI_UNUSED((/mgfft,ngfft(1),npwin,npwout,ldx,ldy,ldz,option/))
 ABI_UNUSED((/denpot(1,1,1),fofgin(1,1),fofgout(1,1),fofr(1,1),weight_r,weight_i/))
#endif

end subroutine dfti_fftrisc_mixprec
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_fftug_dp
!! NAME
!! dfti_fftug_dp
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
!! gbound(2*mgfft+8,2)=Table for padded-FFT. See sphereboundary.
!!  ug(npw_k*ndat)=wavefunctions in reciprocal space.
!!
!! OUTPUT
!!  ur(ldx*ldy*ldz*ndat)=wavefunctions in real space.
!!
!! PARENTS
!!      m_dfti
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_fftug_dp(fftalg,fftcache,npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k,gbound,ug,ur)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalg,fftcache
 integer,intent(in) :: npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2),kg_k(3,npw_k)
 real(dp),target,intent(in) :: ug(2*npw_k*ndat)
 real(dp),target,intent(inout) :: ur(2*ldx*ldy*ldz*ndat)

#ifdef HAVE_DFTI
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
 MSG_ERROR("FFT_DFTI support not activated")
 ABI_UNUSED((/fftalg,fftcache,npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k(1,1),gbound(1,1)/))
 ABI_UNUSED((/ug(1),ur(1)/))
#endif

end subroutine dfti_fftug_dp
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_fftug_spc
!! NAME
!! dfti_fftug_spc
!!
!! FUNCTION
!! Compute ndat zero-padded FFTs from G- to R-space .
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
!!  ur(ldx*ldy*ldz*ndat)=wavefunctions in real space.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_fftug_spc(fftalg,fftcache,npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k,gbound,ug,ur)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalg,fftcache
 integer,intent(in) :: npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2),kg_k(3,npw_k)
 complex(spc),target,intent(in) :: ug(npw_k*ndat)
 complex(spc),target,intent(inout) :: ur(ldx*ldy*ldz*ndat)    !vz_i

#ifdef HAVE_DFTI
!Local variables-------------------------------
!scalars
 integer,parameter :: dist=1
 real(sp) :: fofgout(2,0)
 real(sp),ABI_CONTIGUOUS pointer :: real_ug(:,:),real_ur(:,:)

! *************************************************************************

#undef TK_PREF
#define TK_PREF(name) CONCAT(cplx_,name)

#include "fftug.finc"

#else
 ! Silence compiler warning
 MSG_ERROR("FFT_DFTI support not activated")
 ABI_UNUSED((/fftalg,fftcache,npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k(1,1),gbound(1,1)/))
 ABI_UNUSED((/ug(1),ur(1)/))
#endif

end subroutine dfti_fftug_spc
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_fftug_dpc
!! NAME
!! dfti_fftug_dpc
!!
!! FUNCTION
!! Compute ndat zero-padded FFTs from G ro R.
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
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_fftug_dpc(fftalg,fftcache,npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k,gbound,ug,ur)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalg,fftcache
 integer,intent(in) :: npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2),kg_k(3,npw_k)
 complex(dpc),target,intent(in) :: ug(npw_k*ndat)
 complex(dpc),target,intent(inout) :: ur(ldx*ldy*ldz*ndat)    !vz_i

#ifdef HAVE_DFTI
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
 MSG_ERROR("FFT_DFTI support not activated")
 ABI_UNUSED((/fftalg,fftcache,npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k(1,1),gbound(1,1)/))
 ABI_UNUSED((/ug(1),ur(1)/))
#endif

end subroutine dfti_fftug_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_fftur_dp
!! NAME
!! dfti_fftur_dp
!!
!! FUNCTION
!! Compute ndat zero-padded FFTs from R- to G-space .
!! Mainly used for the transform of wavefunctions.
!! TARGET: dp arrays with real and imaginary part.
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
!! ur(2,ldx*ldy*ldz*ndat)= In input: wavefunctions in real space.
!!                         Destroyed in output. Do not use ur anymore!
!! OUTPUT
!! ug(2,npw_k*ndat)=wavefunctions in reciprocal space.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_fftur_dp(fftalg,fftcache,npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k,gbound,ur,ug)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalg,fftcache
 integer,intent(in) :: npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2),kg_k(3,npw_k)
 real(dp),target,intent(inout) :: ur(2*ldx*ldy*ldz*ndat)
 real(dp),target,intent(inout) :: ug(2*npw_k*ndat)    !vz_i

#ifdef HAVE_DFTI
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
 MSG_ERROR("FFT_DFTI support not activated")
 ABI_UNUSED((/fftalg,fftcache/))
 ABI_UNUSED((/npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k(1,1),gbound(1,1)/))
 ABI_UNUSED((/ug(1),ur(1)/))
#endif

end subroutine dfti_fftur_dp
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_fftur_spc
!! NAME
!! dfti_fftur_spc
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
!!
!! OUTPUT
!! ug(npw_k*ndat)=wavefunctions in reciprocal space.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_fftur_spc(fftalg,fftcache,npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k,gbound,ur,ug)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalg,fftcache
 integer,intent(in) :: npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2),kg_k(3,npw_k)
 complex(spc),target,intent(inout) :: ur(ldx*ldy*ldz*ndat)
 complex(spc),target,intent(inout) :: ug(npw_k*ndat)    !vz_i

#ifdef HAVE_DFTI
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
 MSG_ERROR("FFT_DFTI support not activated")
 ABI_UNUSED((/fftalg,fftcache/))
 ABI_UNUSED((/npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k(1,1),gbound(1,1)/))
 ABI_UNUSED((/ug(1),ur(1)/))
#endif

end subroutine dfti_fftur_spc
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_fftur_dpc
!! NAME
!! dfti_fftur_dpc
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
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_fftur_dpc(fftalg,fftcache,npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k,gbound,ur,ug)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalg,fftcache
 integer,intent(in) :: npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2),kg_k(3,npw_k)
 complex(dpc),target,intent(inout) :: ur(ldx*ldy*ldz*ndat)
 complex(dpc),target,intent(inout) :: ug(npw_k*ndat)    !vz_i

#ifdef HAVE_DFTI
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
 MSG_ERROR("FFT_DFTI support not activated")
 ABI_UNUSED((/fftalg,fftcache/))
 ABI_UNUSED((/npw_k,nx,ny,nz,ldx,ldy,ldz,ndat,istwf_k,mgfft,kg_k(1,1),gbound(1,1)/))
 ABI_UNUSED((/ug(1),ur(1)/))
#endif

end subroutine dfti_fftur_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_c2c_ip_spc
!! NAME
!!  dfti_c2c_ip_spc
!!
!! FUNCTION
!! Driver routine for in-place 3D complex-complex FFT. TARGET: SPC arrays
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimensions of the array.
!! ndat=Number of FFTs to be done.
!! isign= +1 : ff(G) => ff(R); -1 : ff(R) => ff(G)
!!
!! SIDE EFFECTS
!!  ff(ldx*ldy*ldz*ndat)=
!!    In input: the complex array to be transformed.
!!    In output: the Fourier transform in the space specified by isign.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_c2c_ip_spc(nx,ny,nz,ldx,ldy,ldz,ndat,isign,ff)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,isign
!arrays
 complex(spc),intent(inout) :: ff(ldx*ldy*ldz*ndat)

! *************************************************************************

! Include Fortran template
#undef DEV_DFTI_PRECISION
#define DEV_DFTI_PRECISION DFTI_SINGLE

#include "dfti_c2c_ip.finc"

end subroutine dfti_c2c_ip_spc
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_c2c_ip_dpc
!! NAME
!!  dfti_c2c_ip_dpc
!!
!! FUNCTION
!! Driver routine for in-place 3D complex-complex FFT. TARGET: DPC arrays
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimensions of the array.
!! ndat=Number of FFTs to be done.
!! isign= +1 : ff(G) => ff(R); -1 : ff(R) => ff(G)
!!
!! SIDE EFFECTS
!!  ff(ldx*ldy*ldz*ndat)=
!!    In input: the complex array to be transformed.
!!    In output: the Fourier transformed in the space specified by isign.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_c2c_ip_dpc(nx,ny,nz,ldx,ldy,ldz,ndat,isign,ff)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,isign
!arrays
 complex(dpc),intent(inout) :: ff(ldx*ldy*ldz*ndat)

! *************************************************************************

! Include Fortran template
#undef DEV_DFTI_PRECISION
#define DEV_DFTI_PRECISION DFTI_DOUBLE

#include "dfti_c2c_ip.finc"

end subroutine dfti_c2c_ip_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_c2c_op_spc
!! NAME
!!  dfti_c2c_op_spc
!!
!! FUNCTION
!! Driver routine for out-of-place 3D complex-complex FFT of lengths nx, ny, nz.
!! TARGET: spc arrays
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimensions of the array.
!! ndat=Number of FFTs to be done.
!! isign= +1 : ff(G) => gg(R); -1 : ff(R) => gg(G)
!! ff(ldx*ldy*ldz*ndat)=The array to be transformed.
!!
!! OUTPUT
!!   gg(ldx*ldy*ldz*ndat)=The FFT of ff.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_c2c_op_spc(nx,ny,nz,ldx,ldy,ldz,ndat,isign,ff,gg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,isign,ndat
!arrays
 complex(spc),intent(in) :: ff(ldx*ldy*ldz*ndat)
 complex(spc),intent(out) :: gg(ldx*ldy*ldz*ndat)

! *************************************************************************

! Include Fortran template
#undef DEV_DFTI_PRECISION
#define DEV_DFTI_PRECISION DFTI_SINGLE

#include "dfti_c2c_op.finc"

end subroutine dfti_c2c_op_spc
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_c2c_op_dpc
!! NAME
!!  dfti_c2c_op_dpc
!!
!! FUNCTION
!! Driver routine for out-of-place 3D complex-complex FFT of lengths nx, ny, nz.
!! TARGET: DPC arrays
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimensions of the array.
!! ndat=Number of FFTs to be done.
!! isign= +1 : ff(G) => gg(R); -1 : ff(R) => gg(G)
!! ff(ldx*ldy*ldz*ndat)=The array to be transformed.
!!
!! OUTPUT
!!   gg(ldx*ldy*ldz*ndat)=The FFT of ff.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_c2c_op_dpc(nx,ny,nz,ldx,ldy,ldz,ndat,isign,ff,gg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,isign,ndat
!arrays
 complex(dpc),intent(in) :: ff(ldx*ldy*ldz*ndat)
 complex(dpc),intent(out) :: gg(ldx*ldy*ldz*ndat)

! *************************************************************************

! Include Fortran template
#undef DEV_DFTI_PRECISION
#define DEV_DFTI_PRECISION DFTI_DOUBLE

#include "dfti_c2c_op.finc"

end subroutine dfti_c2c_op_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_many_dft_op
!! NAME
!!  dfti_many_dft_op
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
!!
!! OUTPUT
!! fout(2,ldx*ldy*ldz*ndat)=The Fourier transform of fin.
!!
!! PARENTS
!!      m_dfti
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_many_dft_op(nx,ny,nz,ldx,ldy,ldz,ndat,isign,fin,fout)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,isign
!arrays
 real(dp),target,intent(in) :: fin(2*ldx*ldy*ldz*ndat)
 real(dp),target,intent(out) :: fout(2*ldx*ldy*ldz*ndat)

#ifdef HAVE_DFTI
!Local variables-------------------------------
!scalars
 type(C_ptr) :: fin_cptr, fout_cptr
!arrays
 complex(dpc),ABI_CONTIGUOUS pointer :: fin_fptr(:),fout_fptr(:)

! *************************************************************************

 ! Associate complex pointers with real inputs via the C pointers
 fin_cptr = C_loc(fin)
 call C_F_pointer(fin_cptr,fin_fptr, shape=(/ldx*ldy*ldz*ndat/))

 fout_cptr = C_loc(fout)
 call C_F_pointer(fout_cptr,fout_fptr, shape=(/ldx*ldy*ldz*ndat/))

 ! Call complex version --> a lot of boilerplate code avoided
 call dfti_c2c_op(nx,ny,nz,ldx,ldy,ldz,ndat,isign,fin_fptr,fout_fptr)

#else
 MSG_ERROR("FFT_DFTI support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,ndat,isign/))
 ABI_UNUSED(fin(1))
 ABI_UNUSED(fout(1))
#endif

end subroutine dfti_many_dft_op
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_many_dft_ip
!! NAME
!!  dfti_many_dft_ip
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
!!
!! OUTPUT
!! finout(2,ldx*ldy*ldz*ndat)=
!!   In input: The complex array to be transformed.
!!   In output: The FFT results.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_many_dft_ip(nx,ny,nz,ldx,ldy,ldz,ndat,isign,finout)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,isign
!arrays
 real(dp),target,intent(inout) :: finout(2*ldx*ldy*ldz*ndat)

#ifdef HAVE_DFTI
!Local variables-------------------------------
!scalars
 type(C_ptr) :: finout_cptr
!arrays
 complex(dpc),ABI_CONTIGUOUS pointer :: finout_fptr(:)

! *************************************************************************

 ! Associate complex finout_fptr with real ffinout via the C pointer
 finout_cptr = C_loc(finout)
 call C_F_pointer(finout_cptr,finout_fptr, shape=(/ldx*ldy*ldz*ndat/))

 ! Call complex version --> a lot of boilerplate code avoided
 call dfti_c2c_ip(nx,ny,nz,ldx,ldy,ldz,ndat,isign,finout_fptr)

#else
 MSG_ERROR("FFT_DFTI support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,ndat,isign/))
 ABI_UNUSED(finout(1))
#endif

end subroutine dfti_many_dft_ip
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_fftpad_dp
!! NAME
!!  dfti_fftpad_dp
!!
!! FUNCTION
!!  This routine transforms wavefunctions using 3D zero-padded FFTs with DFTI.
!!  The 3D ffts are computed only on lines and planes which have non zero elements (see zpad_init)
!!  FFT transform is in-place.
!!
!! INPUTS
!!   nx,ny,nz=Logical dimensions of the FFT mesh.
!!   ldx,ldy,ldz=Physical dimension of the f array (to avoid cache conflicts).
!!   ndat=Numer of FFTs
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
!!      m_dfti
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_fftpad_dp(ff,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign,gbound)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2)
 real(dp),target,intent(inout) :: ff(2*ldx*ldy*ldz*ndat)

!Local variables-------------------------------
#ifdef HAVE_DFTI
!scalars
 type(C_ptr) :: cptr
!arrays
 complex(dpc),ABI_CONTIGUOUS pointer :: fptr(:)

! *************************************************************************

 ! Associate complex fptr with real ff via the C pointer
 cptr = C_loc(ff)
 call C_F_pointer(cptr,fptr, SHAPE=(/ldx*ldy*ldz*ndat/))

 ! Call complex version --> a lot of boilerplate code avoided
 call dfti_fftpad_dpc(fptr,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign,gbound)

#else
 MSG_ERROR("FFT_DFTI support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign/))
 ABI_UNUSED(gbound(1,1))
 ABI_UNUSED(ff(1))
#endif

end subroutine dfti_fftpad_dp
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_fftpad_dpc
!! NAME
!!  dfti_fftpad_dpc
!!
!! FUNCTION
!!  This routine transforms wavefunctions using 3D zero-padded FFTs with DFTI.
!!  The 3D ffts are computed only on lines and planes which have non zero elements (see zpad_init)
!!  FFT transform is in-place. Target: complex DPC arrays.
!!
!! INPUTS
!!   nx,ny,nz=Logical dimensions of the FFT mesh.
!!   ldx,ldy,ldz=Physical dimension of the f array (to avoid cache conflicts).
!!   ndat=Number of FFTs.
!!   mgfft=MAX(nx,ny,nz), only used to dimension gbound
!!   isign=The sign of the transform.
!!   gbound(2*mgfft+8,2)= The boundaries of the basis sphere of G vectors at a given k-point.
!!     See sphereboundary for more info.
!!
!! SIDE EFFECTS
!!   ff(ldx*ldy*ldz*ndat)=
!!     input: The array with the data to be transformed.
!!     output: The results of the FFT.
!!
!! PARENTS
!!      m_dfti
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_fftpad_dpc(ff,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign,gbound)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2)
 complex(dpc),intent(inout) :: ff(ldx*ldy*ldz*ndat)

!Local variables-------------------------------
!scalars
#ifdef HAVE_DFTI

! *************************************************************************

! Include Fortran template
#undef DEV_DFTI_PRECISION
#define DEV_DFTI_PRECISION DFTI_DOUBLE

#include "dfti_fftpad.finc"

#else
 MSG_ERROR("FFT_DFTI support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign/))
 ABI_UNUSED(gbound(1,1))
 ABI_UNUSED(ff(1))
#endif

end subroutine dfti_fftpad_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_fftpad_spc
!! NAME
!!  dfti_fftpad_spc
!!
!! FUNCTION
!!  This routine transforms wavefunctions using 3D zero-padded FFTs with DFTI.
!!  The 3D ffts are computed only on lines and planes which have non zero elements (see zpad_init)
!!  FFT transform is in-place. Target: complex SPC arrays.
!!
!! INPUTS
!!   nx,ny,nz=Logical dimensions of the FFT mesh.
!!   ldx,ldy,ldz=Physical dimension of the f array (to avoid cache conflicts).
!!   ndat=Number of FFTs.
!!   mgfft=MAX(nx,ny,nz), only used to dimension gbound
!!   isign=The sign of the transform.
!!   gbound(2*mgfft+8,2)= The boundaries of the basis sphere of G vectors at a given k-point.
!!     See sphereboundary for more info.
!!
!! SIDE EFFECTS
!!   ff(ldx*ldy*ldz*ndat)=
!!     input: The array with the data to be transformed.
!!     output: The results of the FFT.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_fftpad_spc(ff,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign,gbound)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2)
 complex(spc),intent(inout) :: ff(ldx*ldy*ldz*ndat)
#ifdef HAVE_DFTI

! *************************************************************************

! Include Fortran template
#undef DEV_DFTI_PRECISION
#define DEV_DFTI_PRECISION DFTI_SINGLE

#include "dfti_fftpad.finc"

#else
 MSG_ERROR("FFT_DFTI support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign/))
 ABI_UNUSED(gbound(1,1))
 ABI_UNUSED(ff(1))
#endif

end subroutine dfti_fftpad_spc
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_r2c_op_dpc
!! NAME
!!  dfti_r2c_op_dpc
!!
!! FUNCTION
!! Driver routine for out-of-place 3D real-to-complex FFT of lengths nx, ny, nz.
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimensions of the f array (to avoid cache conflicts).
!! ff(ldx*ldy*ldz*ndat)=The real array to be transformed.
!! ndat=Number of FFTs to be done.
!!
!! OUTPUT
!! gg(ldx*ldy*ldz*ndat)=The forward FFT of ff (complex valued)
!!
!! PARENTS
!!      m_dfti
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_r2c_op_dpc(nx,ny,nz,ldx,ldy,ldz,ndat,ff,gg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat
!arrays
 real(dp),intent(in) :: ff(ldx*ldy*ldz*ndat)
 complex(dpc),intent(out) :: gg(ldx*ldy*ldz*ndat)

#ifdef HAVE_DFTI
!Local variables-------------------------------
!scalars
 integer :: status,nhp,padx,i1,i2,i3,igp,igf,imgf,ii
 integer :: i1inv,i2inv,i3inv,idat,padatf
 type(DFTI_DESCRIPTOR),pointer :: Desc
 type(C_PTR) :: cptr
!arrays
 integer,allocatable :: i1inver(:),i2inver(:),i3inver(:)
 complex(dpc),ABI_CONTIGUOUS pointer :: gg_hp(:)

! *************************************************************************

 padx = (nx/2+1)
 nhp = (nx/2+1)*ny*nz

 call dfti_alloc_complex(nhp*ndat,cptr,gg_hp)

 status = DftiCreateDescriptor(Desc, DFTI_DOUBLE, DFTI_REAL, 3, (/nx,ny,nz/) )
 DFTI_CHECK(status)

 status = DftiSetValue(Desc, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
 status = DftiSetValue(Desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
 status = DftiSetValue(Desc, DFTI_INPUT_STRIDES,  (/0, 1, ldx,  ldx*ldy/))
 status = DftiSetValue(Desc, DFTI_INPUT_DISTANCE, ldx*ldy*ldz)
 status = DftiSetValue(Desc, DFTI_OUTPUT_STRIDES, (/0, 1, padx, padx*ny/))
 status = DftiSetValue(Desc, DFTI_OUTPUT_DISTANCE, nhp)
 status = DftiSetValue(Desc, DFTI_NUMBER_OF_TRANSFORMS, ndat)
 status = DftiSetValue(Desc, DFTI_FORWARD_SCALE, one / DBLE(nx*ny*nz) )
 DFTI_CHECK(status)

 DFTI_CHECK( DftiCommitDescriptor(Desc) )

 DFTI_CHECK( DftiComputeForward(Desc, ff, gg_hp) )

 DFTI_CHECK( DftiFreeDescriptor(Desc) )

 ! Reconstruct full FFT: Hermitian redundancy: out[i] is the conjugate of out[n-i]
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
   padatf=(idat-1)*ldx*ldy*ldz
   do i3=1,nz
     i3inv = i3inver(i3)
     do i2=1,ny
       i2inv = i2inver(i2)
       do i1=1,padx
         igp=igp+1
         igf = i1 + (i3-1)*ldx*ldy + (i2-1)*ldx + padatf
         gg(igf) =  gg_hp(igp)
         i1inv = i1inver(i1)
         if (i1inv/=i1) then
           imgf = i1inv + (i3inv-1)*ldx*ldy + (i2inv-1)*ldx + padatf
           gg(imgf) = DCONJG(gg_hp(igp))
         end if
       end do
     end do
   end do
 end do

 ABI_FREE(i1inver)
 ABI_FREE(i2inver)
 ABI_FREE(i3inver)

 call dfti_free(cptr)

#else
 MSG_ERROR("FFT_DFTI support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz/))
 ABI_UNUSED(ff)
 ABI_UNUSED(gg(1))
#endif

end subroutine dfti_r2c_op_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_r2c_op_dp
!! NAME
!!  dfti_r2c_op_dp
!!
!! FUNCTION
!! Driver routine for out-of-place 3D real-to-complex FFT of lengths nx, ny, nz.
!!
!! INPUTS
!! nx,ny,nz=Number of points along the three directions.
!! ldx,ldy,ldz=Physical dimensions of the f array (to avoid cache conflicts).
!! ndat=Number of FFTs to be done.
!! ff(ldx*ldy*ldz*ndat)=The real array to be transformed.
!!
!! OUTPUT
!! gg(2*ldx*ldy*ldz*ndat)=The forward FFT of ff (real valued)
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_r2c_op_dp(nx,ny,nz,ldx,ldy,ldz,ndat,ff,gg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat
!arrays
 real(dp),intent(in) :: ff(ldx*ldy*ldz*ndat)
 real(dp),target,intent(out) :: gg(2*ldx*ldy*ldz*ndat)

#ifdef HAVE_DFTI
!Local variables-------------------------------
!scalars
 type(C_ptr) :: gg_cptr
!arrays
 complex(dpc),ABI_CONTIGUOUS pointer :: gg_fptr(:)

! *************************************************************************

 gg_cptr = C_loc(gg)
 call C_F_pointer(gg_cptr,gg_fptr, shape=(/ldx*ldy*ldz*ndat/))

 call dfti_r2c_op_dpc(nx,ny,nz,ldx,ldy,ldz,ndat,ff,gg_fptr)

#else
 MSG_ERROR("FFT_DFTI support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz/))
 ABI_UNUSED(ff)
 ABI_UNUSED(gg(1))
#endif

end subroutine dfti_r2c_op_dp
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_c2r_op_dpc
!! NAME
!!  dfti_c2r_op_dpc
!!
!! FUNCTION
!! Driver routine for out-of-place 3D complex-to-real FFT of lengths nx, ny, nz.
!!
!! INPUTS
!! nx,ny,nz=Number of point along the three directions.
!! ldx,ldy,ldz=Physical dimension of the f array (to avoid cache conflicts).
!! ndat=Number of FFTs to be done.
!! ff(2,ldx*ldy*ldz*ndat)=The complex array to be transformed.
!!
!! OUTPUT
!! gg(2,ldx*ldy*ldz*ndat)=The backwards real FFT of ff.
!!
!! PARENTS
!!      m_dfti
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_c2r_op_dpc(nx,ny,nz,ldx,ldy,ldz,ndat,ff,gg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat
!arrays
 complex(dpc),intent(in) :: ff(ldx*ldy*ldz*ndat)
 real(dp),intent(out) :: gg(ldx*ldy*ldz*ndat)

#ifdef HAVE_DFTI
!Local variables-------------------------------
!scalars
 integer :: status,nhp,padx,i2,i3,igp,igf,idat,padatf,padatp,ii
 type(DFTI_DESCRIPTOR),pointer :: Desc
 type(C_PTR) :: cptr
!arrays
 complex(dpc),ABI_CONTIGUOUS pointer :: ff_hp(:)

! *************************************************************************

 !stride  = 1
 !idist   = nhp
 !odist   = nx*ny*nz
 !n       = (/nx,ny,nz/)
 !inembed = (/(nx/2+1),ny,nz/)
 !onembed = (/nx,ny,nz/) ! check this
 !my_plan = retrieve_plan3(n,ndat,inembed,stride,idist,onembed,stride,odist,FFTW_BACKWARD,my_flags,Saved_plans)

 ! Fill the Hermitian part: Hermitian redundancy: out[i] is the conjugate of out[n-i]
 padx    = (nx/2+1)
 nhp     = padx*ny*nz

 call dfti_alloc_complex(nhp*ndat,cptr,ff_hp)

!!$OMP PARALLEL DO PRIVATE(padatf,padatp,igf,igp)
 do idat=1,ndat
   padatf=(idat-1)*ldx*ldy*ldz
   padatp=(idat-1)*padx*ny*nz
   do i3=1,nz
     do i2=1,ny
       igf = (i3-1)*ldx*ldy + (i2-1)*ldx   + padatf
       igp = (i3-1)*padx*ny + (i2-1)*padx  + padatp
       ff_hp(igp+1:igp+padx) = ff(igf+1:igf+padx)
     end do
   end do
 end do

 status = DftiCreateDescriptor(Desc, DFTI_DOUBLE, DFTI_REAL, 3, (/nx,ny,nz/) )
 DFTI_CHECK(status)

 status = DftiSetValue(Desc, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
 status = DftiSetValue(Desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
 status = DftiSetValue(Desc, DFTI_INPUT_STRIDES, (/0, 1, padx, padx*ny/))
 status = DftiSetValue(Desc, DFTI_INPUT_DISTANCE, nhp)
 status = DftiSetValue(Desc, DFTI_OUTPUT_STRIDES, (/0, 1, ldx, ldx*ldy/))
 status = DftiSetValue(Desc, DFTI_OUTPUT_DISTANCE, ldx*ldy*ldz)
 status = DftiSetValue(Desc, DFTI_NUMBER_OF_TRANSFORMS, ndat)
 DFTI_CHECK(status)

 DFTI_CHECK( DftiCommitDescriptor(Desc) )

 DFTI_CHECK( DftiComputeBackward(Desc, ff_hp, gg) )

 DFTI_CHECK( DftiFreeDescriptor(Desc) )

 call dfti_free(cptr)

#else
 MSG_ERROR("FFT_DFTI support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz/))
 ABI_UNUSED(ff(1))
 ABI_UNUSED(gg(1))
#endif

end subroutine dfti_c2r_op_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_c2r_op_dp
!! NAME
!!  dfti_c2r_op_dp
!!
!! FUNCTION
!! Driver routine for out-of-place 3D complex-to-real FFT of lengths nx, ny, nz.
!!
!! INPUTS
!! nx,ny,nz=Number of point along the three directions.
!! ldx,ldy,ldz=Physical dimension of the f array (to avoid cache conflicts).
!! ndat=Number of FFTs to be done.
!! ff(2,ldx*ldy*ldz*ndat)=The complex array to be transformed.
!!
!! OUTPUT
!! gg(ldx*ldy*ldz*ndat)=The backwards real FFT of ff.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_c2r_op_dp(nx,ny,nz,ldx,ldy,ldz,ndat,ff,gg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat
!arrays
 real(dp),target,intent(in) :: ff(2*ldx*ldy*ldz*ndat)
 real(dp),intent(inout) :: gg(ldx*ldy*ldz*ndat)    !vz_i

#ifdef HAVE_DFTI
!Local variables-------------------------------
!scalars
 type(C_ptr) :: ff_cptr
!arrays
 complex(dpc),ABI_CONTIGUOUS pointer :: ff_fptr(:)

! *************************************************************************

 ff_cptr = C_loc(ff)
 call C_F_pointer(ff_cptr,ff_fptr, shape=(/ldx*ldy*ldz*ndat/))

 call dfti_c2r_op_dpc(nx,ny,nz,ldx,ldy,ldz,ndat,ff_fptr,gg)

#else
 MSG_ERROR("FFT_DFTI support not activated")
 ABI_UNUSED((/nx,ny,nz,ldx,ldy,ldz/))
 ABI_UNUSED((/ff(1),gg(1)/))
#endif

end subroutine dfti_c2r_op_dp
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_check_status
!! NAME
!!  dfti_check_status
!!
!! FUNCTION
!!  Error handler for DFTI wrappers. Print error message and abort.
!!
!! INPUTS
!!  status = status error reported by dfti.
!!  [file] = file name
!!  [line] = line number
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

#ifdef HAVE_DFTI

subroutine dfti_check_status(status,file,line)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: status
 integer,optional,intent(in) :: line
 character(len=*),optional,intent(in) :: file

!Local variables-------------------------------
 integer :: f90line
 character(len=10) :: lnum
 character(len=500) :: f90name
 character(len=500) :: my_msg
 character(len=DFTI_MAX_MESSAGE_LENGTH+500) :: err_msg

! *************************************************************************

 if (PRESENT(line)) then
   f90line=line
 else
   f90line=0
 end if
 call int2char10(f90line,lnum)

 if (PRESENT(file)) then
   f90name = basename(file)
 else
   f90name='Subroutine Unknown'
 end if

 my_msg=strcat(f90name,":",lnum,":")

 if (status /= 0) then
   if (.not. DftiErrorClass(status, DFTI_NO_ERROR)) then
     err_msg = strcat(my_msg," Error: ",DftiErrorMessage(status))
     MSG_ERROR(err_msg)
   end if
 end if

end subroutine dfti_check_status
!!***

#endif

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_spawn_threads_here
!! NAME
!! dfti_spawn_threads_here
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

function dfti_spawn_threads_here(ndat,nthreads) result(ans)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndat,nthreads
 logical :: ans
!arrays

! *************************************************************************

 ans = .FALSE.
#ifdef HAVE_OPENMP
 ans = (nthreads > 1 .and. MOD(ndat,nthreads) == 0 .and. .not. USE_LIB_THREADS)
#else
 ABI_UNUSED((/ndat,nthreads/))
#endif

end function dfti_spawn_threads_here
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_use_lib_threads
!! NAME
!! dfti_use_lib_threads
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine dfti_use_lib_threads(logvar)

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: logvar

! *************************************************************************

 USE_LIB_THREADS = logvar

end subroutine dfti_use_lib_threads
!!***

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_alloc_real_dp
!! NAME
!! dfti_alloc_real_dp
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

#ifdef HAVE_DFTI

subroutine dfti_alloc_real_dp(size,cptr,fptr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: size
 real(dp),ABI_CONTIGUOUS pointer :: fptr(:)
 type(C_PTR),intent(out) :: cptr
!arrays

! *************************************************************************

 cptr = mkl_malloc( INT(size*C_DOUBLE, KIND=C_SIZE_T), DFTI_DEFAULT_ALIGNMENT_DP)
 if (.not. C_ASSOCIATED(cptr)) then
   MSG_ERROR("mkl_malloc returned NULL!")
 end if

 call c_f_pointer(cptr, fptr, [size])

end subroutine dfti_alloc_real_dp
!!***

#endif

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_alloc_complex_spc
!! NAME
!! dfti_alloc_complex_spc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

#ifdef HAVE_DFTI

subroutine dfti_alloc_complex_spc(size,cptr,fptr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: size
 complex(spc),ABI_CONTIGUOUS pointer :: fptr(:)
 type(C_PTR),intent(out) :: cptr
!arrays

! *************************************************************************

 cptr = mkl_malloc( INT(2*size*C_FLOAT, KIND=C_SIZE_T), DFTI_DEFAULT_ALIGNMENT_SP)
 if (.not. C_ASSOCIATED(cptr)) then
   MSG_ERROR("mkl_malloc returned NULL!")
 end if

 call c_f_pointer(cptr, fptr, [size])

end subroutine dfti_alloc_complex_spc
!!***

#endif

!----------------------------------------------------------------------

!!****f* m_dfti/dfti_alloc_complex_dpc
!! NAME
!! dfti_alloc_complex_dpc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

#ifdef HAVE_DFTI

subroutine dfti_alloc_complex_dpc(size,cptr,fptr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: size
 complex(dpc),ABI_CONTIGUOUS pointer :: fptr(:)
 type(C_PTR),intent(out) :: cptr
!arrays

! *************************************************************************

 cptr = mkl_malloc( INT(2*size*C_DOUBLE, KIND=C_SIZE_T), DFTI_DEFAULT_ALIGNMENT_DP)
 if (.not. C_ASSOCIATED(cptr)) then
   MSG_ERROR("mkl_malloc returned NULL!")
 end if

 call c_f_pointer(cptr, fptr, [size])

end subroutine dfti_alloc_complex_dpc
!!***

#endif

!----------------------------------------------------------------------

END MODULE m_dfti
!!***
