!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_FFT_prof
!! NAME
!! m_FFT_prof
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_FFT_prof

 use defs_basis
 use m_xomp
 use m_errors
 use m_abicore
 use m_fftw3
 use m_fft
 use m_distribfft

 use defs_abitypes,    only : MPI_type
 use m_numeric_tools,  only : arth
 use m_time,           only : cwtime
 use m_io_tools,       only : open_file
 use m_geometry,       only : metric
 use m_hide_blas,      only : xcopy
 use m_cgtools,        only : set_istwfk
 use m_fftcore,        only : get_kg, print_ngfft, fftalg_info, kgindex, getng, sphereboundary
 use m_fft_mesh,       only : calc_eigr, calc_ceigr
 use m_mpinfo,         only : nullify_mpi_enreg, destroy_mpi_enreg, copy_mpi_enreg, initmpi_seq
 use m_oscillators,    only : rho_tw_g

 implicit none

 private

 integer,private,parameter :: TNAME_LEN=100
!!***
!----------------------------------------------------------------------

!!****t* m_fft_mesh/FFT_test_t
!! NAME
!! FFT_test_t
!!
!! FUNCTION
!! Structure storing the set of paramenters passed to the FFT routines used in
!! abinit (fourdp|fourwf).
!!
!! SOURCE

 type,public :: FFT_test_t
   integer :: available=0
   integer :: istwf_k=-1
   integer :: mgfft=-1
   integer :: ndat=-1
   integer :: nfft=-1
   integer :: nthreads=1
   integer :: npw_k=-1
   integer :: npw_kout=-1
   integer :: paral_kgb=-1

   real(dp) :: ecut=zero

   integer :: ngfft(18)=-1

   real(dp) :: kpoint(3) = [zero,zero,zero]
   real(dp) :: rprimd(3,3),rmet(3,3)
   real(dp) :: gprimd(3,3),gmet(3,3)

   integer,allocatable :: kg_k(:,:)
   integer,allocatable :: kg_kout(:,:)
   integer,allocatable :: indpw_k(:)

   type(MPI_type) :: MPI_enreg
 end type FFT_test_t

 public :: fft_test_init
 public :: fft_test_nullify
 public :: fft_test_free
 public :: fft_test_print
!!***

 interface fft_test_free
   module procedure fft_test_free_0D
   module procedure fft_test_free_1D
 end interface fft_test_free

!----------------------------------------------------------------------

!!****t* m_fft_mesh/FFT_prof_t
!! NAME
!! FFT_prof_t
!!
!! FUNCTION
!!  The results of the tests
!!
!! SOURCE

 type,public :: FFT_prof_t
   integer :: ncalls
   integer :: ndat
   integer :: nthreads
   real(dp) :: cpu_time
   real(dp) :: wall_time
   real(dp) :: gflops
   character(len=TNAME_LEN) :: test_name
   complex(dpc),allocatable :: results(:)
 end type FFT_prof_t

 public :: fftprof_init
 public :: fftprof_free
 public :: fftprof_print

 interface fftprof_free
   module procedure fftprof_free_0D
   module procedure fftprof_free_1D
 end interface fftprof_free
!!***

!----------------------------------------------------------------------

 public :: fftprof_ncalls_per_test
 !
 ! Timing routines.
 public :: time_fourdp
 public :: time_fftbox
 public :: time_fourwf
 public :: time_rhotwg
 public :: time_fftu

 ! Routines for benchmarks.
 public :: prof_fourdp
 public :: prof_fourwf
 public :: prof_rhotwg

!----------------------------------------------------------------------
 ! Number of calls of each FFT algo, used to have a betters statistics for timing.
 integer,save,private :: NCALLS_FOR_TEST=10

 integer,private,parameter :: CACHE_KBSIZE=0
! Argument of empty_cache. Set it to zero if the cache should not be emptied.

CONTAINS  !====================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/fft_test_init
!! NAME
!!  fft_test_init
!!
!! FUNCTION
!!  Creation method for the FFT_test_t structured datatype.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof,m_fft_prof
!!
!! CHILDREN
!!      fft_test_free,fft_test_init,fft_test_nullify,fftprof_free,get_kg
!!      random_number,time_rhotwg
!!
!! SOURCE

subroutine fft_test_init(Ftest,fft_setup,kpoint,ecut,boxcutmin,rprimd,nsym,symrel,MPI_enreg_in)

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: nsym
 real(dp),intent(in) :: ecut,boxcutmin
 type(FFT_test_t),intent(inout) :: Ftest
 type(MPI_type),intent(in) :: MPI_enreg_in
!arrays
 integer,intent(in) :: fft_setup(5),symrel(3,3,nsym)
 real(dp),intent(in) :: kpoint(3),rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: fftalg,fftcache,ndat
 real(dp) :: ucvol
!arrays
 logical,allocatable :: mask(:)
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)

! *************************************************************************

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 !@FFT_test_t
 Ftest%rprimd = rprimd
 Ftest%rmet   = rmet
 Ftest%gprimd = gprimd
 Ftest%gmet   = gmet
 Ftest%ecut   = ecut

 fftalg   = fft_setup(1)
 fftcache = fft_setup(2)
 ndat     = fft_setup(3)

 Ftest%nthreads  = fft_setup(4)
 Ftest%available = fft_setup(5)

 Ftest%paral_kgb = 0
 Ftest%kpoint    = kpoint
 Ftest%ndat      = ndat

 Ftest%istwf_k = set_istwfk(kpoint)

 call get_kg(Ftest%kpoint,Ftest%istwf_k,ecut,gmet,Ftest%npw_k,Ftest%kg_k)
 call get_kg(Ftest%kpoint,Ftest%istwf_k,ecut,gmet,Ftest%npw_kout,Ftest%kg_kout)

 call copy_mpi_enreg(MPI_enreg_in,Ftest%MPI_enreg)

 Ftest%ngfft(7) = fftalg
 Ftest%ngfft(8) = fftcache

 ! Fill part of ngfft
 call getng(boxcutmin,ecut,gmet,k0,Ftest%MPI_enreg%me_fft,Ftest%mgfft,Ftest%nfft,Ftest%ngfft,Ftest%MPI_enreg%nproc_fft,nsym,&
&  Ftest%MPI_enreg%paral_kgb,symrel, unit=dev_null)

 call init_distribfft(Ftest%MPI_enreg%distribfft,'c',Ftest%MPI_enreg%nproc_fft,Ftest%ngfft(2),Ftest%ngfft(3))

 ! Compute the index of each plane wave in the FFT grid.
 ABI_MALLOC(Ftest%indpw_k,(Ftest%npw_k))

 ABI_MALLOC(mask,(Ftest%npw_k))
 call kgindex(Ftest%indpw_k,Ftest%kg_k,mask,Ftest%MPI_enreg,Ftest%ngfft,Ftest%npw_k)
 ABI_CHECK(ALL(mask),"FFT parallelism not supported in fftprof")
 ABI_FREE(mask)

end subroutine fft_test_init
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/fft_test_nullify
!! NAME
!!  fft_test_nullify
!!
!! FUNCTION
!!  Nullify all pointers.
!!
!! INPUTS
!!
!! PARENTS
!!      fftprof,m_fft_prof
!!
!! CHILDREN
!!      fft_test_free,fft_test_init,fft_test_nullify,fftprof_free,get_kg
!!      random_number,time_rhotwg
!!
!! SOURCE

subroutine fft_test_nullify(Ftest)

!Arguments -----------------------------------
!scalars
 type(FFT_test_t),intent(inout) :: Ftest

! *************************************************************************

 ! @FFT_test_t
 call nullify_mpi_enreg(Ftest%MPI_enreg)

end subroutine fft_test_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/fft_test_free_0D
!! NAME
!!  fft_test_free_0D
!!
!! FUNCTION
!!  Destruction method for the FFT_test_t structured datatype.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_fft_prof
!!
!! CHILDREN
!!      fft_test_free,fft_test_init,fft_test_nullify,fftprof_free,get_kg
!!      random_number,time_rhotwg
!!
!! SOURCE

subroutine fft_test_free_0D(Ftest)

!Arguments -----------------------------------
!scalars
 type(FFT_test_t),intent(inout) :: Ftest

! *********************************************************************

 !@FFT_test_t
 Ftest%available=0
 Ftest%istwf_k=-1
 Ftest%mgfft=-1
 Ftest%ndat=-1
 Ftest%nfft=-1
 Ftest%nthreads=1
 Ftest%npw_k=-1
 Ftest%npw_kout=-1
 Ftest%paral_kgb=-1

 Ftest%ecut=zero

 Ftest%ngfft=-1

 Ftest%kpoint =zero
 Ftest%rprimd =zero
 Ftest%rmet   =zero
 Ftest%gprimd =zero
 Ftest%gmet   =zero

 if (allocated(Ftest%indpw_k)) then
   ABI_FREE(Ftest%indpw_k)
 end if
 if (allocated(Ftest%kg_k)) then
   ABI_FREE(Ftest%kg_k)
 end if
 if (allocated(Ftest%kg_kout)) then
   ABI_FREE(Ftest%kg_kout)
 end if

 call destroy_mpi_enreg(Ftest%MPI_enreg)

end subroutine fft_test_free_0D
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/fft_test_free_1D
!! NAME
!!  fft_test_free_1D
!!
!! FUNCTION
!!  Destruction method for the FFT_test_t structured datatype.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      fft_test_free,fft_test_init,fft_test_nullify,fftprof_free,get_kg
!!      random_number,time_rhotwg
!!
!! SOURCE

subroutine fft_test_free_1D(Ftest)

!Arguments -----------------------------------
!scalars
 type(FFT_test_t),intent(inout) :: Ftest(:)

!Local variables-------------------------------
!scalars
 integer :: ii

! *********************************************************************

 do ii=LBOUND(Ftest,DIM=1),UBOUND(Ftest,DIM=1)
   call fft_test_free_0D(Ftest(ii))
 end do

end subroutine fft_test_free_1D
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/fft_test_print
!! NAME
!!  fft_test_print
!!
!! FUNCTION
!!  Printout of the FFT_test_t structured datatype.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof
!!
!! CHILDREN
!!      fft_test_free,fft_test_init,fft_test_nullify,fftprof_free,get_kg
!!      random_number,time_rhotwg
!!
!! SOURCE

subroutine fft_test_print(Ftest,header,unit,mode_paral,prtvol)

!Arguments -----------------------------------
!scalars
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header
 type(FFT_test_t),intent(in) :: Ftest

!Local variables-------------------------------
!scalars
 integer :: my_unt,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg
! *********************************************************************

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 !msg=' ==== Info on the FFT test object ==== '
 if (PRESENT(header)) then
   msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
   call wrtout(my_unt,msg,my_mode)
 end if

 !TODO add additional info
 write(msg,'(a,i3)')"FFT setup for fftalg ",Ftest%ngfft(7)
 call print_ngfft(Ftest%ngfft,header=msg,unit=my_unt,mode_paral="COLL")

end subroutine fft_test_print
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/name_of
!! NAME
!!  name_of
!!
!! FUNCTION
!!  Returns a string with info on the test.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

function name_of(Ftest)

!Arguments -----------------------------------
!scalars
 character(len=TNAME_LEN) :: name_of
 type(FFT_test_t),intent(in) :: Ftest

!Local variables-------------------------------
!scalars
 character(len=TNAME_LEN) :: library_name,cplex_mode,padding_mode

! *********************************************************************

 call fftalg_info(Ftest%ngfft(7),library_name,cplex_mode,padding_mode)
 !name_of = TRIM(library_name)//"; "//TRIM(cplex_mode)//"; "//TRIM(padding_mode)

 write(name_of,'(i3)')Ftest%ngfft(7)
 name_of = TRIM(library_name)//" ("//TRIM(name_of)//")"

end function name_of
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/fftprof_init
!! NAME
!!  fftprof_init
!!
!! FUNCTION
!!  Creation method for the FFT_prof_t structured datatype.
!!
!! INPUTS
!!  gflops = Gigaflops
!!
!! OUTPUT
!!
!! PARENTS
!!      m_fft_prof
!!
!! CHILDREN
!!      fft_test_free,fft_test_init,fft_test_nullify,fftprof_free,get_kg
!!      random_number,time_rhotwg
!!
!! SOURCE

subroutine fftprof_init(Ftprof,test_name,nthreads,ncalls,ndat,cpu_time,wall_time,gflops,results)

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: ncalls,nthreads,ndat
 real(dp),intent(in) :: cpu_time,wall_time,gflops
 character(len=*),intent(in) :: test_name
 type(FFT_prof_t),intent(out) :: Ftprof
!arrays
 complex(dpc),optional,intent(in) :: results(:)

! *************************************************************************

 !@FFT_prof_t
 Ftprof%ncalls    =  ncalls
 Ftprof%nthreads  =  nthreads
 Ftprof%ndat      = ndat
 Ftprof%cpu_time  = cpu_time
 Ftprof%wall_time = wall_time
 Ftprof%gflops    = gflops
 Ftprof%test_name = test_name

 if (PRESENT(results)) then
   if (allocated(Ftprof%results)) then
     ABI_FREE(Ftprof%results)
   end if
   ABI_MALLOC(Ftprof%results,(SIZE(results)))
   Ftprof%results = results
 end if

end subroutine fftprof_init
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/fftprof_free_0D
!! NAME
!!  fftprof_free_0D
!!
!! FUNCTION
!!  Destruction method for the FFT_prof_t structured datatype.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_fft_prof
!!
!! CHILDREN
!!      fft_test_free,fft_test_init,fft_test_nullify,fftprof_free,get_kg
!!      random_number,time_rhotwg
!!
!! SOURCE

subroutine fftprof_free_0D(Ftprof)

!Arguments -----------------------------------
!scalars
 type(FFT_prof_t),intent(inout) :: Ftprof

! *********************************************************************

 !@FFT_prof_t
 Ftprof%ncalls=0
 Ftprof%nthreads=0
 Ftprof%cpu_time=zero
 Ftprof%wall_time=zero
 Ftprof%test_name = "None"

 if (allocated(Ftprof%results))  then
   ABI_FREE(Ftprof%results)
 end if

end subroutine fftprof_free_0D
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/fftprof_free_1D
!! NAME
!!  fftprof_free_1D
!!
!! FUNCTION
!!  Destruction method for the FFT_prof_t structured datatype.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      fft_test_free,fft_test_init,fft_test_nullify,fftprof_free,get_kg
!!      random_number,time_rhotwg
!!
!! SOURCE

subroutine fftprof_free_1D(Ftprof)

!Arguments -----------------------------------
!scalars
 type(FFT_prof_t),intent(inout) :: Ftprof(:)

!Local variables-------------------------------
!scalars
 integer :: ii
! *********************************************************************

 !@FFT_prof_t
 do ii=LBOUND(Ftprof,DIM=1),UBOUND(Ftprof,DIM=1)
   call fftprof_free_0D(Ftprof(ii))
 end do

end subroutine fftprof_free_1D
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/fftprof_print
!! NAME
!!  fftprof_print
!!
!! FUNCTION
!!  Printout of the FFT_prof_t structured datatype.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof
!!
!! CHILDREN
!!      fft_test_free,fft_test_init,fft_test_nullify,fftprof_free,get_kg
!!      random_number,time_rhotwg
!!
!! SOURCE

subroutine fftprof_print(Fprof,header,unit,mode_paral,prtvol)

!Arguments -----------------------------------
!scalars
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header
 type(FFT_prof_t),intent(in) :: Fprof(:)

!Local variables-------------------------------
!scalars
 integer :: my_unt,my_prtvol,ncalls
 integer :: field1_w,ii,ref_lib !ifft
 real(dp) :: mabs_err,mean_err,check_mabs_err,check_mean_err
 real(dp) :: ref_wtime,para_eff
 character(len=4) :: my_mode
 character(len=500) :: ofmt,hfmt,nafmt,msg

! *********************************************************************

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg='==== Info on the FFT_prof_t object ===='
 if (PRESENT(header)) msg='==== '//TRIM(ADJUSTL(header))//' ===='

 call wrtout(my_unt,ch10//REPEAT("=",LEN_TRIM(msg)))
 call wrtout(my_unt,msg,my_mode)
 call wrtout(my_unt,REPEAT("=",LEN_TRIM(msg)))

 field1_w=0 ! Width of the field used to print key names.
 do ii=1,SIZE(Fprof)
   field1_w = MAX(field1_w, LEN_TRIM(Fprof(ii)%test_name))
 end do

 ! TODO Add gflops
 if (field1_w==0) RETURN
 field1_w = field1_w + 2 ! To account for ". "
 write(ofmt,*)"(a",field1_w,",2x,2(f7.4,4x),1x,i2,1x,a,i3,a,1x,i0,4x,2(es9.2,3x))"

 write(hfmt,*)"(a",field1_w,",2x,a)"
 write(std_out,hfmt)" Library      ","CPU-time   WALL-time   nthreads  ncalls  Max_|Err|   <|Err|>"
 !
 ! Find reference library.
 ref_lib=0
 do ii=1,SIZE(Fprof)
   if (Fprof(ii)%ncalls>0) then
     ref_lib = ii
     EXIT
   end if
 end do
 !ref_lib=3
 !
 ! Write timing analysis and error wrt reference library if available.
 check_mabs_err=zero; check_mean_err=zero
 do ii=1,SIZE(Fprof)
   ncalls = Fprof(ii)%ncalls
   if (ncalls>0) then
     mabs_err = zero; mean_err=zero
     if (ref_lib>0) then
       mabs_err = MAXVAL( ABS(Fprof(ii)%results - Fprof(ref_lib)%results) )
       mean_err = SUM( ABS(Fprof(ii)%results - Fprof(ref_lib)%results) ) / SIZE(Fprof(ref_lib)%results)
       ! Relative error is not a good estimator because some components are close to zero within machine accuracy.
       !mean_err = 100 * MAXVAL( ABS(Fprof(ii)%results - Fprof(1)%results)/ ABS(Fprof(1)%results ))
       !ifft = imax_loc( ABS(Fprof(ii)%results - Fprof(1)%results)/ ABS(Fprof(1)%results) )
       !write(std_out,*) Fprof(ii)%results(ifft),Fprof(1)%results(ifft)
     end if
     if (Fprof(ii)%nthreads==1) ref_wtime = Fprof(ii)%wall_time
     para_eff = 100 * ref_wtime / ( Fprof(ii)%wall_time * Fprof(ii)%nthreads)
     write(std_out,ofmt)&
&      "- "//Fprof(ii)%test_name,Fprof(ii)%cpu_time/ncalls,Fprof(ii)%wall_time/ncalls,&
&       Fprof(ii)%nthreads,"(",NINT(para_eff),"%)",ncalls,mabs_err,mean_err
     check_mabs_err = MAX(check_mabs_err, mabs_err)
     check_mean_err = MAX(check_mean_err, mean_err)
   else
     write(nafmt,*)"(a",field1_w,",2x,a)"
     write(std_out,nafmt)"- "//Fprof(ii)%test_name,"   N/A        N/A        N/A     N/A       N/A        N/A"
   end if
 end do

 if (ref_lib>0) then
   write(std_out,'(/,2(a,es9.2),2a)')&
&    " Consistency check: MAX(Max_|Err|) = ",check_mabs_err,&
&    ", Max(<|Err|>) = ",check_mean_err,", reference_lib: ",TRIM(Fprof(ref_lib)%test_name)
 end if
 write(std_out,*)

end subroutine fftprof_print
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/time_fourdp
!! NAME
!!  time_fourdp
!!
!! FUNCTION
!!  Profiling of the fourdp routine.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof,m_fft_prof
!!
!! CHILDREN
!!      fft_test_free,fft_test_init,fft_test_nullify,fftprof_free,get_kg
!!      random_number,time_rhotwg
!!
!! SOURCE

subroutine time_fourdp(Ftest,isign,cplex,header,Ftprof)

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: cplex,isign
 character(len=500),intent(out) :: header
 type(FFT_test_t),intent(inout) :: Ftest
 type(FFT_prof_t),intent(out) :: Ftprof

!Local variables-------------------------------
!scalars
 integer,parameter :: nspinor1=1
 integer :: icall,i1,i2,i3,n1,n2,n3,ifft
 real(dp) :: cpu_time,wall_time,gflops
 real(dp) :: gsq
 character(len=500) :: msg
 character(len=TNAME_LEN) :: test_name
!arrays
 integer :: gg(3)
 integer,parameter :: g0(3)=(/1,2,-1/)
 real(dp),allocatable :: fofg(:,:),fofr(:)
 complex(dpc),allocatable :: results(:),ctmp(:)

! *********************************************************************

 test_name = name_of(Ftest)
 n1=Ftest%ngfft(1)
 n2=Ftest%ngfft(2)
 n3=Ftest%ngfft(3)

 write(header,'(2(a,i2),a)')" fourdp with cplex ",cplex,", isign ",isign,", ndat 1"

 if (Ftest%available==0) then
   call fftprof_init(Ftprof,test_name,0,0,0,zero,zero,zero)
   RETURN
 end if

 call xomp_set_num_threads(Ftest%nthreads)

 if (Ftest%ngfft(7)/100 == FFT_FFTW3) then
   call fftw3_set_nthreads(Ftest%nthreads)
 end if

 ABI_MALLOC(fofg,(2,Ftest%nfft))
 ABI_MALLOC(fofr,(cplex*Ftest%nfft))
 !
 ! Initialize input data.
 if (isign==1) then ! initialize fofg
   fofg = zero
   ifft=0
   do i3=1,n3
     gg(3)=i3-1; if (i3>1+n3/2) gg(3)=i3-n3-1 ! TODO recheck this
     do i2=1,n2
       gg(2)=i2-1; if (i2>1+n2/2) gg(2)=i2-n2-1
       do i1=1,n1
         gg(1)=i1-1; if (i1>1+n1/2) gg(1)=i1-n1-1
         gsq = two_pi**2 * DOT_PRODUCT(gg,MATMUL(Ftest%gmet,gg))
         ifft=ifft+1
         fofg(1,ifft) = EXP(-gsq)
         fofg(2,ifft) = zero
       end do
     end do
   end do

 else ! init fofr
   if (cplex==2) then
     call calc_eigr(g0,Ftest%nfft,Ftest%ngfft,fofr)
   else if (cplex==1) then
     ABI_MALLOC(ctmp,(Ftest%nfft))
     call calc_ceigr(g0,Ftest%nfft,nspinor1,Ftest%ngfft,ctmp)
     fofr = REAL(ctmp)
     ABI_FREE(ctmp)
   else
     write(msg,'(a,i0)')" Wrong cplex: ",cplex
     MSG_ERROR(msg)
   end if
 end if

 ABI_MALLOC(results,(Ftest%nfft))
 results=czero

 call cwtime(cpu_time,wall_time,gflops,"start")

 do icall=1,NCALLS_FOR_TEST
   !
   ifft = empty_cache(CACHE_KBSIZE)
   call fourdp(cplex,fofg,fofr,isign,Ftest%MPI_enreg,Ftest%nfft,1,Ftest%ngfft,0)
   !
   ! Store results at the first call.
   if (icall==1) then
     if (isign==1) then
       if (cplex==1) then
         results = DCMPLX(fofr, zero)
       else if (cplex==2) then
         results = DCMPLX( fofr(1:2*Ftest%nfft:2), fofr(2:2*Ftest%nfft:2) )
       end if
     else if (isign==-1) then
       results = DCMPLX(fofg(1,:),fofg(2,:))
     end if
   end if
 end do

 call cwtime(cpu_time,wall_time,gflops,"stop")

 ABI_FREE(fofg)
 ABI_FREE(fofr)

 call fftprof_init(Ftprof,test_name,Ftest%nthreads,NCALLS_FOR_TEST,&
&  Ftest%ndat,cpu_time,wall_time,gflops,results=results)

 ABI_FREE(results)

end subroutine time_fourdp
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/time_fftbox
!! NAME
!!  time_fftbox
!!
!! FUNCTION
!!  Profiling of the fftbox_[io]p routines.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof
!!
!! CHILDREN
!!      fft_test_free,fft_test_init,fft_test_nullify,fftprof_free,get_kg
!!      random_number,time_rhotwg
!!
!! SOURCE

subroutine time_fftbox(Ftest,isign,inplace,header,Ftprof)

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: isign,inplace
 character(len=500),intent(inout) :: header !vz_i
 type(FFT_test_t),intent(inout) :: Ftest
 type(FFT_prof_t),intent(out) :: Ftprof

!Local variables-------------------------------
!scalars
 integer,parameter :: nspinor1=1
 integer :: icall,i1,i2,i3,n1,n2,n3,ifft,ndat,nfft,dat,padat
 real(dp) :: cpu_time,wall_time,gflops,gsq
 character(len=500) :: msg
 character(len=TNAME_LEN) :: test_name
 type(fftbox_plan3_t) :: plan
!arrays
 integer,parameter :: g0(3) = (/1,-2,1/)
 integer :: gg(3)
 complex(dpc),allocatable :: ffc(:),ggc(:),results(:)
! *********************************************************************

 test_name = name_of(Ftest)

 if (Ftest%available==0) then
   call fftprof_init(Ftprof,test_name,0,0,0,zero,zero,zero)
   RETURN
 end if

 call xomp_set_num_threads(Ftest%nthreads)

 if (Ftest%ngfft(7)/100 == FFT_FFTW3) then
   call fftw3_set_nthreads(Ftest%nthreads)
 end if

 n1=Ftest%ngfft(1); n2=Ftest%ngfft(2); n3=Ftest%ngfft(3)
 nfft = Ftest%nfft
 ndat = Ftest%ndat

 write(header,'(3(a,i2))')" fftbox with isign ",isign,", in-place ",inplace,", ndat ",ndat

 ABI_MALLOC(ffc,(nfft*ndat))
 ABI_MALLOC(ggc,(nfft*ndat))
 ABI_MALLOC(results,(nfft*ndat))
 ffc=czero; ggc=czero; results=czero

 if (isign==-1) then
   call calc_ceigr(g0,nfft,nspinor1,Ftest%ngfft,ffc)
 else if (isign==1) then
   ifft=0
   do i3=1,n3
     gg(3)=i3-1; if (i3>1+n3/2) gg(3)=i3-n3-1 ! TODO recheck this
     do i2=1,n2
       gg(2)=i2-1; if (i2>1+n2/2) gg(2)=i2-n2-1
       do i1=1,n1
         gg(1)=i1-1; if (i1>1+n1/2) gg(1)=i1-n1-1
         gsq = two_pi**2 * DOT_PRODUCT(gg,MATMUL(Ftest%gmet,gg))
         ifft=ifft+1
         ffc(ifft) = EXP(-gsq)
       end do
     end do
   end do
 else
   MSG_ERROR("Wrong isign")
 end if

 ! Replicate and scale input data
 do dat=2,ndat
   padat = (dat-1) * nfft
   do ifft=1,nfft
     ffc(ifft+padat) = DBLE(dat) * ffc(ifft)
   end do
 end do

 call cwtime(cpu_time,wall_time,gflops,"start")

 ! No augmentation here.
 call fftbox_plan3_many(plan,ndat,Ftest%ngfft(1:3),Ftest%ngfft(1:3),Ftest%ngfft(7),isign)

 select case (inplace)
 case (0)
   do icall=1,NCALLS_FOR_TEST
     ifft = empty_cache(CACHE_KBSIZE)
     call fftbox_execute(plan,ffc,ggc)
     ! Store results at the first call.
     if (icall==1) results = ggc
   end do
 case (1)
   do icall=1,NCALLS_FOR_TEST
     ifft = empty_cache(CACHE_KBSIZE)
     call fftbox_execute(plan,ffc)
     ! Store results at the first call.
     if (icall==1) results = ffc
   end do
 case default
   write(msg,'(a,i0)')" Wrong value for inplace= ",inplace
   MSG_ERROR(msg)
 end select

 call cwtime(cpu_time,wall_time,gflops,"stop")

 ABI_FREE(ffc)
 ABI_FREE(ggc)

 call fftprof_init(Ftprof,test_name,Ftest%nthreads,NCALLS_FOR_TEST,Ftest%ndat,&
&  cpu_time,wall_time,gflops,results=results)

 ABI_FREE(results)

end subroutine time_fftbox
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/time_fourwf
!! NAME
!!  time_fourwf
!!
!! FUNCTION
!!  Profiling of the fourwf routine.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof,m_fft_prof
!!
!! CHILDREN
!!      fft_test_free,fft_test_init,fft_test_nullify,fftprof_free,get_kg
!!      random_number,time_rhotwg
!!
!! SOURCE

subroutine time_fourwf(Ftest,cplex,option_fourwf,header,Ftprof)

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: cplex,option_fourwf
 character(len=500),intent(out) :: header
 type(FFT_test_t),intent(inout) :: Ftest
 type(FFT_prof_t),intent(out) :: Ftprof

!Local variables-------------------------------
!scalars
 integer,parameter :: tim0=0
 integer :: n1,n2,n3,n4,n5,n6,npw_out,icall,i1,i2,i3,idx,ipw,ndat,cnt,dat,padat
 integer :: fftalg,fftalga,fftalgc
 real(dp),parameter :: weight_i=one,weight_r=one
 real(dp) :: cpu_time,wall_time,gflops,gsq,g0dotr
 logical :: isbuggy,not_supported
 character(len=500) :: msg
 character(len=TNAME_LEN) :: test_name
!arrays
 integer,parameter :: g0(3)=(/1,-1,2/) !g0(3)=(/1,0,0/)
 integer :: gg(3)
 integer,allocatable :: gbound_in(:,:),gbound_out(:,:)
 real(dp),allocatable :: denpot(:,:,:),fofg_in(:,:)
 real(dp),allocatable :: fofr_4(:,:,:,:),fofg_out(:,:)
 complex(dpc),allocatable :: results(:)

! *********************************************************************

 test_name = name_of(Ftest)

 fftalg  = Ftest%ngfft(7); fftalga = fftalg/100; fftalgc = MOD(fftalg,10)

 isbuggy = &
&   (option_fourwf==3 .and. fftalga==FFT_SG2002 .and. fftalgc /= 0 .and. any(Ftest%istwf_k == [3,4,5,6,7,8,9]))  ! see sg_fourwf
 !isbuggy = .False.

 !FIXME problems with the unitary tests reference files!
 not_supported = ( &
&  (fftalgc == 2 .and. option_fourwf==0 .and. Ftest%istwf_k >  2) &
&                )
 !not_supported = .FALSE.

 npw_out= Ftest%npw_kout
 ndat   = Ftest%ndat
 n1=Ftest%ngfft(1); n2=Ftest%ngfft(2); n3=Ftest%ngfft(3)
 n4=Ftest%ngfft(4); n5=Ftest%ngfft(5); n6=Ftest%ngfft(6)

 write(header,'(4(a,i2))')" fourwf with option ",option_fourwf,", cplex ",cplex,", ndat ",ndat,", istwf_k ",Ftest%istwf_k

 if (isbuggy .or. not_supported .or. Ftest%available==0) then
   call fftprof_init(Ftprof,test_name,0,0,0,zero,zero,zero)
   RETURN
 end if

 call xomp_set_num_threads(Ftest%nthreads)

 if (Ftest%ngfft(7)/100 == FFT_FFTW3) then
   call fftw3_set_nthreads(Ftest%nthreads)
 end if

 ! FFTW3 does not need gbound_in but oh well
 ABI_MALLOC(gbound_in,(2*Ftest%mgfft+8,2))
 call sphereboundary(gbound_in,Ftest%istwf_k,Ftest%kg_k,Ftest%mgfft,Ftest%npw_k)

 ABI_MALLOC(gbound_out,(2*Ftest%mgfft+8,2))
 call sphereboundary(gbound_out,Ftest%istwf_k,Ftest%kg_kout,Ftest%mgfft,Ftest%npw_kout)

 ABI_MALLOC(denpot,(cplex*n4,n5,n6))
 ABI_MALLOC(fofg_in,(2,Ftest%npw_k*ndat))
 ABI_MALLOC(fofg_out,(2,npw_out*ndat))
 ABI_MALLOC(fofr_4,(2,n4,n5,n6*ndat))

 denpot   = zero
 fofg_in  = zero
 fofg_out = zero
 fofr_4   = zero

 ABI_MALLOC(results,(Ftest%nfft*ndat))
 results=czero

 select case (option_fourwf)
 case (0,1,2)
   !! for option==0, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
   !!                fofr(2,n4,n5,n6) contains the output Fourier Transform of fofgin;
   !!                no use of denpot, fofgout and npwout.
   !! for option==1, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
   !!                denpot(cplex*n4,n5,n6) contains the input density at input,
   !!                and the updated density at output (accumulated);
   !!                no use of fofgout and npwout.
   !! for option==2, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
   !!                denpot(cplex*n4,n5,n6) contains the input local potential;
   !!                fofgout(2,npwout*ndat) contains the output function;
   !!
   do cnt=0,(Ftest%npw_k * ndat) - 1
     ipw = 1 + MOD(cnt, Ftest%npw_k)
     gg = Ftest%kg_k(:,ipw)
     gsq = two_pi**2 * DOT_PRODUCT(gg,MATMUL(Ftest%gmet,gg))
     fofg_in(1,cnt+1) = EXP(-gsq)
     fofg_in(2,cnt+1) = zero
   end do

   if (option_fourwf==1) then ! Init denpot
     denpot = one
   end if
   !
   if (option_fourwf==2) then ! Init denpot
     !
     if (cplex==1) then
       do i3=0,n3-1
         do i2=0,n2-1
           do i1=0,n1-1
             g0dotr= two_pi*( g0(1)*(i1/DBLE(n1)) &
&                            +g0(2)*(i2/DBLE(n2)) &
&                            +g0(3)*(i3/DBLE(n3)) )
             denpot(i1+1,i2+1,i3+1)=COS(g0dotr)
           end do
         end do
       end do
     else if (cplex==2) then
       do i3=0,n3-1
         do i2=0,n2-1
           idx=1
           do i1=0,n1-1
             g0dotr= two_pi*( g0(1)*(i1/DBLE(n1)) &
&                            +g0(2)*(i2/DBLE(n2)) &
&                            +g0(3)*(i3/DBLE(n3)) )

             denpot(idx,  i2+1,i3+1)= COS(g0dotr)
             denpot(idx+1,i2+1,i3+1)= SIN(g0dotr)
             idx=idx+2
           end do
         end do
       end do
     end if
   end if

 case (3)
   !! for option==3, fofr(2,n4,n5,n6*ndat) contains the input real space wavefunction;
   !!                fofgout(2,npwout*ndat) contains its output Fourier transform;
   !!                no use of fofgin and npwin.
   idx=0
   do dat=1,ndat
     padat = (dat-1)*n6
     do i3=0,n3-1
       do i2=0,n2-1
         do i1=0,n1-1
           g0dotr= two_pi*( g0(1)*(i1/DBLE(n1)) &
                           +g0(2)*(i2/DBLE(n2)) &
                           +g0(3)*(i3/DBLE(n3)) )
           !fofr_4(1,i1+1,i2+1,i3+1+padat)=COS(g0dotr)
           !fofr_4(2,i1+1,i2+1,i3+1+padat)=SIN(g0dotr)
           fofr_4(1,i1+1,i2+1,i3+1+padat)=10 * EXP(-g0dotr**2)
           fofr_4(2,i1+1,i2+1,i3+1+padat)=zero
           idx=idx+1
           results(idx) = DCMPLX(fofr_4(1,i1+1,i2+1,i3+1+padat), fofr_4(2,i1+1,i2+1,i3+1+padat))
         end do
       end do
     end do
   end do

 case default
   write(msg,'(a,i0)')" Wrong value for option_fourwf: ",option_fourwf
   MSG_ERROR(msg)
 end select

 call cwtime(cpu_time,wall_time,gflops,"start")

 do icall=1,NCALLS_FOR_TEST

   i1 = empty_cache(CACHE_KBSIZE)

   call fourwf(cplex,denpot,fofg_in,fofg_out,fofr_4,gbound_in,gbound_out,Ftest%istwf_k,&
&    Ftest%kg_k,Ftest%kg_kout,Ftest%mgfft,Ftest%MPI_enreg,ndat,Ftest%ngfft,Ftest%npw_k,npw_out,n4,n5,n6,option_fourwf,&
&    tim0,weight_r,weight_i)
   !
   ! Store results at the first call.
   if (icall==1) then
     select case (option_fourwf)
     case (0)
       !! for option==0, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
       !!                fofr(2,n4,n5,n6) contains the output Fourier Transform of fofgin;
       !!                no use of denpot, fofgout and npwout.
       idx=0
       do dat=1,ndat
         padat = (dat-1)*n6
         do i3=1,n3
           do i2=1,n2
             do i1=1,n1
               idx=idx+1
               results(idx) = DCMPLX(fofr_4(1,i1,i2,i3+padat),fofr_4(2,i1,i2,i3+padat))
             end do
           end do
         end do
       end do
     case (1)
       !! for option==1, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
       !!                denpot(cplex*n4,n5,n6) contains the input density at input,
       !!                and the updated density at output (accumulated);
       !!                no use of fofgout and npwout.
       if (cplex==1) then
        idx=0
        do i3=1,n3
          do i2=1,n2
            do i1=1,n1
              idx=idx+1
              results(idx) = DCMPLX(denpot(i1,i2,i3),zero)
            end do
          end do
        end do

       else if (cplex==2) then
         idx=0
         do i3=1,n3
           do i2=1,n2
             do i1=1,2*n1,2
               idx=idx+1
               results(idx) = DCMPLX(denpot(i1,i2,i3),denpot(i1+1,i2,i3))
             end do
           end do
         end do
       else
         MSG_ERROR("Wrong cplex")
       end if

     case (2)
       !! for option==2, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
       !!                denpot(cplex*n4,n5,n6) contains the input local potential;
       !!                fofgout(2,npwout*ndat) contains the output function;
       do ipw=1,npw_out*ndat
         results(ipw) = DCMPLX(fofg_out(1,ipw),fofg_out(2,ipw))
       end do
     case (3)
       !! for option==3, fofr(2,n4,n5,n6*ndat) contains the input real space wavefunction;
       !!                fofgout(2,npwout*ndat) contains its output Fourier transform;
       !!                no use of fofgin and npwin.
       do ipw=1,npw_out*ndat
         results(ipw) = DCMPLX(fofg_out(1,ipw),fofg_out(2,ipw))
       end do
       !write(Ftest%ngfft(7),*)"results opt 3 Ftest%ngfft(7)",Ftest%ngfft(7)
       !do dat=1,ndat
       !  do ipw=1,npw_out
       !    idx = ipw + (dat-1)*npw_out
       !    write(Ftest%ngfft(7),*)ipw,dat,results(idx)
       !  end do
       !end do
     end select
   end if
 end do

 call cwtime(cpu_time,wall_time,gflops,"stop")

 ABI_FREE(denpot)
 ABI_FREE(fofg_in)
 ABI_FREE(fofg_out)
 ABI_FREE(fofr_4)
 ABI_FREE(gbound_in)
 ABI_FREE(gbound_out)

 call fftprof_init(Ftprof,test_name,Ftest%nthreads,NCALLS_FOR_TEST,&
&   Ftest%ndat,cpu_time,wall_time,gflops,results=results)

 ABI_FREE(results)

end subroutine time_fourwf
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/fftprof_ncalls_per_test
!! NAME
!!  fftprof_ncalls_per_test
!!
!! FUNCTION
!!  Helper function used to set the number of calls to  be used in each time_* routine.
!!
!! INPUTS
!!   ncalls=Number of calls to be used.
!!
!! SIDE EFFECTS
!!  NCALLS_FOR_TEST = ncalls
!!
!! PARENTS
!!      fftprof
!!
!! CHILDREN
!!      fft_test_free,fft_test_init,fft_test_nullify,fftprof_free,get_kg
!!      random_number,time_rhotwg
!!
!! SOURCE

subroutine fftprof_ncalls_per_test(ncalls)

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: ncalls

! *********************************************************************

 NCALLS_FOR_TEST = ncalls

end subroutine fftprof_ncalls_per_test
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/time_rhotwg
!! NAME
!!  time_rhotwg
!!
!! FUNCTION
!!  Profiling of the rho_tw_g routine.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof,m_fft_prof
!!
!! CHILDREN
!!      fft_test_free,fft_test_init,fft_test_nullify,fftprof_free,get_kg
!!      random_number,time_rhotwg
!!
!! SOURCE

subroutine time_rhotwg(Ftest,map2sphere,use_padfft,osc_npw,osc_gvec,header,Ftprof)

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: map2sphere,use_padfft,osc_npw
 character(len=500),intent(out) :: header
 type(FFT_test_t),intent(inout) :: Ftest
 type(FFT_prof_t),intent(out) :: Ftprof
!arrays
 integer,intent(in) :: osc_gvec(3,osc_npw)

!Local variables-------------------------------
!scalars
 integer,parameter :: nspinor1=1,dim_rtwg=1,istwfk1=1
 integer :: icall,ifft,itim1,itim2,nfft,dat,sprc,ptr,ndat
 integer :: n1,n2,n3,n4,n5,n6
 real(dp) :: cpu_time,wall_time,gflops
 complex(dpc) :: ktabp1=cone,ktabp2=cone
 character(len=TNAME_LEN) :: test_name
 logical :: not_implemented
 type(MPI_type) :: MPI_enreg_seq
!arrays
 !integer,parameter :: g1(3)=(/1,1,2/),g2(3)=(/-1,-2,-1/)
 integer,parameter :: g1(3)=(/-1,0,0/),g2(3)=(/1,0,0/)
 integer,allocatable :: gbound(:,:)
 integer,allocatable :: ktabr1(:),ktabr2(:)
 integer,allocatable :: igfftg0(:)
 real(dp),parameter :: spinrot1(4)=(/one,zero,zero,one/),spinrot2(4)=(/one,zero,zero,one/)
 logical,allocatable :: mask(:)
 complex(dpc),allocatable :: results(:)
 complex(gwpc),allocatable :: rhotwg(:)
 complex(gwpc),allocatable :: wfn1(:),wfn2(:)

! *********************************************************************

 test_name = name_of(Ftest)

 nfft = Ftest%nfft
 ndat = Ftest%ndat
 n1=Ftest%ngfft(1); n2=Ftest%ngfft(2); n3=Ftest%ngfft(3)
 n4=Ftest%ngfft(4); n5=Ftest%ngfft(5); n6=Ftest%ngfft(6)

 write(header,'(3(a,i2))')"rho_tw_g with use_padfft ",use_padfft,", map2sphere ",map2sphere,", ndat ",ndat

 ! TODO: zero-pad not available with SG2001 routines.
 not_implemented = (use_padfft==1.and.Ftest%ngfft(7) == 412)

 if (Ftest%available==0.or.not_implemented) then
   call fftprof_init(Ftprof,test_name,0,0,0,zero,zero,zero)
   RETURN
 end if

 call xomp_set_num_threads(Ftest%nthreads)

 if (Ftest%ngfft(7)/100 == FFT_FFTW3) then
   call fftw3_set_nthreads(Ftest%nthreads)
 end if

 call initmpi_seq(MPI_enreg_seq)
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',n2,n3,'all')

 itim1=1; itim2=1
 ABI_MALLOC(ktabr1,(nfft))
 ABI_MALLOC(ktabr2,(nfft))

 do ifft=1,nfft
   ktabr1(ifft)= ifft
   ktabr2(ifft)= ifft
 end do

 ABI_MALLOC(igfftg0,(osc_npw*map2sphere))

 if (map2sphere>0) then
   ABI_MALLOC(mask,(osc_npw))
   call kgindex(igfftg0,osc_gvec,mask,MPI_enreg_seq,Ftest%ngfft,osc_npw)
   ABI_CHECK(ALL(mask)," FFT parallelism not supported")
   ABI_FREE(mask)
 end if

 ABI_MALLOC(gbound,(2*Ftest%mgfft+8,2*use_padfft))
 if (use_padfft==1) then
   call sphereboundary(gbound,istwfk1,osc_gvec,Ftest%mgfft,osc_npw)
 end if

 ABI_MALLOC(wfn1,(nfft*nspinor1*ndat))
 ABI_MALLOC(wfn2,(nfft*nspinor1*ndat))

 do dat=1,ndat
   do sprc=1,nspinor1
     ptr = 1 + (sprc-1)*nfft + (dat-1)*nfft*nspinor1
     call calc_ceigr(g1,nfft,nspinor1,Ftest%ngfft,wfn1(ptr:))
     call calc_ceigr(g2,nfft,nspinor1,Ftest%ngfft,wfn2(ptr:))
   end do
 end do

 ABI_MALLOC(rhotwg,(osc_npw*dim_rtwg*ndat))
 ABI_MALLOC(results,(osc_npw*dim_rtwg*ndat))

 call cwtime(cpu_time,wall_time,gflops,"start")

 do icall=1,NCALLS_FOR_TEST

     ifft = empty_cache(CACHE_KBSIZE)
     call rho_tw_g(nspinor1,osc_npw,nfft,ndat,Ftest%ngfft,map2sphere,use_padfft,igfftg0,gbound,&
&      wfn1,itim1,ktabr1,ktabp1,spinrot1,&
&      wfn2,itim2,ktabr2,ktabp2,spinrot2,&
&      dim_rtwg,rhotwg)
     ! Store results at the first call.
     if (icall==1) results = rhotwg
 end do

 call cwtime(cpu_time,wall_time,gflops,"stop")

 ABI_FREE(ktabr1)
 ABI_FREE(ktabr2)
 ABI_FREE(igfftg0)
 ABI_FREE(gbound)
 ABI_FREE(rhotwg)
 ABI_FREE(wfn1)
 ABI_FREE(wfn2)

 call fftprof_init(Ftprof,test_name,Ftest%nthreads,NCALLS_FOR_TEST,&
&  Ftest%ndat,cpu_time,wall_time,gflops,results)

 call destroy_mpi_enreg(MPI_enreg_seq)
 ABI_FREE(results)

end subroutine time_rhotwg
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/time_fftu
!! NAME
!!  time_fftu
!!
!! FUNCTION
!!  Profiling of the fftu routines.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof
!!
!! CHILDREN
!!      fft_test_free,fft_test_init,fft_test_nullify,fftprof_free,get_kg
!!      random_number,time_rhotwg
!!
!! SOURCE

subroutine time_fftu(Ftest,isign,header,Ftprof)

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: isign
 character(len=500),intent(out) :: header
 type(FFT_test_t),intent(inout) :: Ftest
 type(FFT_prof_t),intent(out) :: Ftprof

!Local variables-------------------------------
!scalars
 integer,parameter :: nspinor1=1
 integer :: icall,i1,i2,i3,n1,n2,n3,ifft,npw_k,ndat,dat,nfft,cnt,ipw,padat,istwf_k
 !integer :: n4,n5,n6,n456,idx
 real(dp) :: cpu_time,wall_time,gflops,gsq,g0dotr
 logical :: not_implemented
 character(len=500) :: msg
 character(len=TNAME_LEN) :: test_name
!arrays
 integer,parameter :: g0(3) = (/1,-2,1/)
 integer :: gg(3)
 integer,allocatable :: kg_k(:,:),gbound(:,:)
 complex(dpc),allocatable :: ug(:),results(:),ur(:)
! *********************************************************************

 test_name = name_of(Ftest)

 n1=Ftest%ngfft(1); n2=Ftest%ngfft(2); n3=Ftest%ngfft(3)
 !n4=Ftest%ngfft(4); n5=Ftest%ngfft(5); n6=Ftest%ngfft(6)

 nfft = n1*n2*n3
 !n456 = n4*n5*n6
 ndat = Ftest%ndat

 write(header,'(2(a,i2))')" fftu with isign ",isign,", ndat ",ndat

 ! TODO: zero-pad not available with SG2001 routines.
 not_implemented = (Ftest%ngfft(7) == 412)

 if (Ftest%available==0.or.not_implemented) then
   call fftprof_init(Ftprof,test_name,0,0,0,zero,zero,zero)
   RETURN
 end if

 call xomp_set_num_threads(Ftest%nthreads)

 if (Ftest%ngfft(7)/100 == FFT_FFTW3) then
   call fftw3_set_nthreads(Ftest%nthreads)
 end if
 !
 istwf_k = Ftest%istwf_k
 !istwf_k = 1 ! Do not use istwf_k trick here
 call get_kg(Ftest%kpoint,istwf_k,Ftest%ecut,Ftest%gmet,npw_k,kg_k)

 ABI_MALLOC(gbound,(2*Ftest%mgfft+8,2))
 call sphereboundary(gbound,istwf_k,kg_k,Ftest%mgfft,npw_k)

 ABI_MALLOC(ug,(npw_k*ndat))
 ABI_MALLOC(ur, (nfft*ndat))
 ABI_MALLOC(results,(nfft*ndat))
 results=czero ! needed for R --> G since npw_k < nfft

 if (isign==1) then
   do cnt=0,(npw_k * ndat) - 1
     ipw = 1 + MOD(cnt, npw_k)
     gg = kg_k(:,ipw)
     gsq = two_pi**2 * DOT_PRODUCT(gg,MATMUL(Ftest%gmet,gg))
     ug(cnt+1) = DCMPLX(EXP(-gsq), zero)
   end do
   !
   ! Replicate input data
   do dat=2,ndat
     padat = (dat-1) * npw_k
     do ipw=1,npw_k
       ug(ipw+padat) = DBLE(dat) * ug(ipw)
     end do
   end do

 else if (isign==-1) then

   do i3=0,n3-1
     do i2=0,n2-1
       do i1=0,n1-1
         g0dotr= two_pi*( g0(1)*(i1/DBLE(n1)) &
&                        +g0(2)*(i2/DBLE(n2)) &
&                        +g0(3)*(i3/DBLE(n3)) )
         ifft = 1 + i1 + i2*n1 + i3*n2*n3
         ur(ifft)=DCMPLX(DCOS(g0dotr),DSIN(g0dotr))
       end do
     end do
   end do
   !
   ! Replicate input data
   do dat=2,ndat
     padat = (dat-1)*nfft
     do ifft=1,nfft
       ur(ifft+padat) = DBLE(dat) * ur(ifft)
     end do
   end do

 else
   write(msg,'(a,i0)')" Wrong isign= ",isign
   MSG_ERROR(msg)
 end if

 call cwtime(cpu_time,wall_time,gflops,"start")

 do icall=1,NCALLS_FOR_TEST
   ifft = empty_cache(CACHE_KBSIZE)

   !call fftpad(ug,Ftest%ngfft,n1,n2,n3,n4,n5,n6,ndat,Ftest%mgfft,isign,gbound)
   if (isign==+1) then
     call fft_ug(npw_k,nfft,nspinor1,ndat,Ftest%mgfft,Ftest%ngfft,istwf_k,kg_k,gbound,ug,ur)
     ! Store results at the first call.
     if (icall==1) then
       call xcopy(nfft*ndat,ur,1,results,1)
     end if

   else
     call fft_ur(npw_k,nfft,nspinor1,ndat,Ftest%mgfft,Ftest%ngfft,istwf_k,kg_k,gbound,ur,ug)
     ! Store results at the first call.
     if (icall==1) then
       call xcopy(npw_k*ndat,ug,1,results,1)
     end if
   end if
   !if (isign==-1) then
   !  write(Ftest%ngfft(7),*)"ngfft, isign ",Ftest%ngfft(7),isign
   !  do idx =1,n4*n5*n6*ndat; write(Ftest%ngfft(7),*)idx,results(idx); end do
   !end if
 end do

 call cwtime(cpu_time,wall_time,gflops,"stop")

 ABI_FREE(kg_k)
 ABI_FREE(gbound)
 ABI_FREE(ug)
 ABI_FREE(ur)

 call fftprof_init(Ftprof,test_name,Ftest%nthreads,NCALLS_FOR_TEST,ndat,&
&   cpu_time,wall_time,gflops,results=results)

 ABI_FREE(results)

end subroutine time_fftu
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/prof_fourdp
!! NAME
!!  prof_fourdp
!!
!! FUNCTION
!!  Profile fourdp
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof
!!
!! CHILDREN
!!      fft_test_free,fft_test_init,fft_test_nullify,fftprof_free,get_kg
!!      random_number,time_rhotwg
!!
!! SOURCE

subroutine prof_fourdp(fft_setups,isign,cplex,necut,ecut_arth,boxcutmin,rprimd,nsym,symrel,MPI_enreg_in)

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: nsym,isign,cplex,necut
 real(dp),intent(in) :: boxcutmin
 type(MPI_type),intent(in) :: MPI_enreg_in
!arrays
 integer,intent(in) :: fft_setups(:,:),symrel(3,3,nsym)
 real(dp),intent(in) :: ecut_arth(2)
 real(dp),intent(in) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: iec,nsetups,set,funt
 type(FFT_test_t) :: Ftest
 type(FFT_prof_t) :: Ftprof
 character(len=500) :: msg,frm,header
 character(len=fnlen) :: fname
!arrays
 integer :: ngfft_ecut(18,necut)
 real(dp),parameter :: k_gamma(3)=zero
 real(dp) :: ecut_list(necut)
 real(dp),allocatable :: prof_res(:,:,:)

! *********************************************************************

 nsetups = SIZE(fft_setups,DIM=2)
 ecut_list = arth(ecut_arth(1),ecut_arth(2),necut)
 !
 ! Open file and write header with info.
 write(fname,'(2(a,i1))')"PROF_fourdp_cplex",cplex,"_isign",isign
 if (open_file(fname,msg,newunit=funt) /= 0) then
   MSG_ERROR(msg)
 end if

 write(msg,'(2(a,i0))')"Benchmark: routine = fourdp, cplex =",cplex,", isign=",isign
 write(std_out,'(a)')" Running "//TRIM(msg)

 write(funt,'(a)')"# "//TRIM(msg)
 do set=1,nsetups
   write(funt,'(a,5(a,i0))') "#",&
&    "  fftalg = "   ,fft_setups(1,set),&
&    ", fftcache = " ,fft_setups(2,set),&
&    ", ndat = "     ,fft_setups(3,set),&
&    ", nthreads = " ,fft_setups(4,set),&
&    ", available = ",fft_setups(5,set)
 end do

 ABI_MALLOC(prof_res,(2,necut,nsetups))

 do set=1,nsetups
   !
   do iec=1,necut
     call fft_test_nullify(Ftest)
     call fft_test_init(Ftest,fft_setups(:,set),k_gamma,ecut_list(iec),boxcutmin,rprimd,nsym,symrel,MPI_enreg_in)

     call time_fourdp(Ftest,isign,cplex,header,Ftprof)

     prof_res(1,iec,set) = Ftprof%cpu_time /Ftprof%ncalls
     prof_res(2,iec,set) = Ftprof%wall_time/Ftprof%ncalls
     !
     ! Save FFT divisions.
     if (set==1) ngfft_ecut(:,iec) = Ftest%ngfft

     call fftprof_free(Ftprof)
     call fft_test_free(Ftest)
   end do
   !
 end do
 !
 ! Write the wall-time as a function of ecut.
 write(frm,*)"(f7.1,3i4,",nsetups,"(f7.4))"
 do iec=1,necut
   write(funt,frm) ecut_list(iec),ngfft_ecut(4:6,iec),prof_res(2,iec,:)
 end do

 close(funt)
 ABI_FREE(prof_res)

end subroutine prof_fourdp
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/prof_fourwf
!! NAME
!!  prof_fourwf
!!
!! FUNCTION
!!  profile fourwf
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof
!!
!! CHILDREN
!!      fft_test_free,fft_test_init,fft_test_nullify,fftprof_free,get_kg
!!      random_number,time_rhotwg
!!
!! SOURCE

subroutine prof_fourwf(fft_setups,cplex,option,kpoint,necut,ecut_arth,boxcutmin,rprimd,nsym,symrel,MPI_enreg_in)

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: nsym,cplex,necut,option
 real(dp),intent(in) :: boxcutmin
 type(MPI_type),intent(in) :: MPI_enreg_in
!arrays
 integer,intent(in) :: fft_setups(:,:),symrel(3,3,nsym)
 real(dp),intent(in) :: ecut_arth(2)
 real(dp),intent(in) :: kpoint(3),rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: iec,nsetups,set,funt,istwf_k
 type(FFT_test_t) :: Ftest
 type(FFT_prof_t) :: Ftprof
 character(len=500) :: msg,frm,header
 character(len=fnlen) :: fname
!arrays
 integer :: ngfft_ecut(18,necut)
 real(dp) :: ecut_list(necut)
 real(dp),allocatable :: prof_res(:,:,:)

! *********************************************************************

 nsetups = SIZE(fft_setups,DIM=2)
 ecut_list = arth(ecut_arth(1),ecut_arth(2),necut)
 istwf_k = set_istwfk(kpoint)
 !
 ! Open file and write header with info.
 write(fname,'(3(a,i1))')"PROF_fourwf_cplex",cplex,"_option",option,"_istwfk",istwf_k
 if (open_file(fname,msg,newunit=funt) /= 0) then
   MSG_ERROR(msg)
 end if

 write(msg,'(3(a,i1))')"Benchmark: routine = fourwf, cplex = ",cplex,", option= ",option,", istwfk= ",istwf_k
 write(std_out,'(a)')" Running "//TRIM(msg)

 write(funt,'(a)')"# "//TRIM(msg)
 do set=1,nsetups
   write(funt,'(a,5(a,i0))') "#",&
&    "  fftalg = "   ,fft_setups(1,set),&
&    ", fftcache = " ,fft_setups(2,set),&
&    ", ndat = "     ,fft_setups(3,set),&
&    ", nthreads = " ,fft_setups(4,set),&
&    ", available = ",fft_setups(5,set)
 end do

 ABI_MALLOC(prof_res,(2,necut,nsetups))

 do set=1,nsetups
   !
   do iec=1,necut
     call fft_test_nullify(Ftest)
     call fft_test_init(Ftest,fft_setups(:,set),kpoint,ecut_list(iec),boxcutmin,rprimd,nsym,symrel,MPI_enreg_in)

     call time_fourwf(Ftest,cplex,option,header,Ftprof)

     prof_res(1,iec,set) = Ftprof%cpu_time /Ftprof%ncalls
     prof_res(2,iec,set) = Ftprof%wall_time/Ftprof%ncalls
     !
     ! Save FFT divisions.
     if (set==1) ngfft_ecut(:,iec) = Ftest%ngfft

     call fftprof_free(Ftprof)
     call fft_test_free(Ftest)
   end do
   !
 end do
 !
 ! Write the wall-time as a function of ecut.
 write(frm,*)"(f7.1,3i4,",nsetups,"(f7.4))"
 do iec=1,necut
   write(funt,frm) ecut_list(iec),ngfft_ecut(4:6,iec),prof_res(2,iec,:)
 end do

 close(funt)
 ABI_FREE(prof_res)

end subroutine prof_fourwf
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/prof_rhotwg
!! NAME
!!  prof_rhotwg
!!
!! FUNCTION
!!  profile rhotwg
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof
!!
!! CHILDREN
!!      fft_test_free,fft_test_init,fft_test_nullify,fftprof_free,get_kg
!!      random_number,time_rhotwg
!!
!! SOURCE

subroutine prof_rhotwg(fft_setups,map2sphere,use_padfft,necut,ecut_arth,osc_ecut,boxcutmin,&
&  rprimd,nsym,symrel,gmet,MPI_enreg_in)

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: nsym,necut,map2sphere,use_padfft
 real(dp),intent(in) :: boxcutmin,osc_ecut
 type(MPI_type),intent(in) :: MPI_enreg_in
!arrays
 integer,intent(in) :: fft_setups(:,:),symrel(3,3,nsym)
 real(dp),intent(in) :: ecut_arth(2)
 real(dp),intent(in) :: rprimd(3,3),gmet(3,3)

!Local variables-------------------------------
!scalars
 integer :: iec,nsetups,set,funt,osc_npw
 type(FFT_test_t) :: Ftest
 type(FFT_prof_t) :: Ftprof
 character(len=500) :: msg,frm,header
 character(len=fnlen) :: fname
!arrays
 integer,allocatable :: osc_gvec(:,:)
 integer :: ngfft_ecut(18,necut)
 real(dp),parameter :: k_gamma(3)=zero
 real(dp) :: ecut_list(necut)
 real(dp),allocatable :: prof_res(:,:,:)

! *********************************************************************

 nsetups = SIZE(fft_setups,DIM=2)
 ecut_list = arth(ecut_arth(1),ecut_arth(2),necut)
 !
 ! Open file and write header with info.
 write(fname,'(2(a,i1))')"PROF_rhotwg_map2sphere",map2sphere,"_use_padfft",use_padfft

 if (open_file(fname,msg,newunit=funt) /= 0) then
   MSG_ERROR(msg)
 end if

 write(msg,'(2(a,i0),a,f5.1)')&
&  "Benchmark: routine = rho_tw_g, map2sphere = ",map2sphere,", use_padfft = ",use_padfft,", osc_ecut = ",osc_ecut
 write(std_out,'(a)')" Running "//TRIM(msg)

 write(funt,'(a)')"# "//TRIM(msg)
 do set=1,nsetups
   write(funt,'(a,5(a,i0))') "#",&
&    "  fftalg = "   ,fft_setups(1,set),&
&    ", fftcache = " ,fft_setups(2,set),&
&    ", ndat = "     ,fft_setups(3,set),&
&    ", nthreads = " ,fft_setups(4,set),&
&    ", available = ",fft_setups(5,set)
 end do

 call get_kg((/zero,zero,zero/),1,osc_ecut,gmet,osc_npw,osc_gvec)
!  TODO should reorder by shells to be consistent with the GW part!
!  Moreover I guess this ordering is more efficient when we have
!  to map the box to the G-sphere!

 ABI_MALLOC(prof_res,(2,necut,nsetups))

 do set=1,nsetups
   !
   do iec=1,necut
     call fft_test_nullify(Ftest)
     call fft_test_init(Ftest,fft_setups(:,set),k_gamma,ecut_list(iec),boxcutmin,rprimd,nsym,symrel,MPI_enreg_in)

     call time_rhotwg(Ftest,map2sphere,use_padfft,osc_npw,osc_gvec,header,Ftprof)

     prof_res(1,iec,set) = Ftprof%cpu_time /Ftprof%ncalls
     prof_res(2,iec,set) = Ftprof%wall_time/Ftprof%ncalls
     !
     ! Save FFT divisions.
     if (set==1) ngfft_ecut(:,iec) = Ftest%ngfft

     call fftprof_free(Ftprof)
     call fft_test_free(Ftest)
   end do
   !
 end do
 !
 ! Write the wall-time as a function of ecut.
 write(frm,*)"(f7.1,3i4,",nsetups,"(f7.4))"
 do iec=1,necut
   write(funt,frm) ecut_list(iec),ngfft_ecut(4:6,iec),prof_res(2,iec,:)
 end do

 close(funt)
 ABI_FREE(prof_res)
 ABI_FREE(osc_gvec)

end subroutine prof_rhotwg
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/empty_cache
!! NAME
!!  empty_cache
!!
!! FUNCTION
!!  Empty the memory cache
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

function empty_cache(kbsize) result(fake)

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: kbsize
 integer :: fake

!Local variables-------------------------------
!scalars
 integer :: sz
!arrays
 real(dp),allocatable :: chunk(:)

! *********************************************************************

 if (kbsize <= 0) RETURN

 sz = int((100._dp * kbsize) / dp)

 ABI_MALLOC(chunk,(sz))
 call random_number(chunk)
 fake = int(SUM(chunk)) ! Need a result, otherwise some smart compiler could skip the call.
 ABI_FREE(chunk)

!----------------------------------------------------------------------

end function empty_cache
!!***

END MODULE m_FFT_prof
!!***
