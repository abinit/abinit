!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_fftcore
!! NAME
!!  m_fftcore
!!
!! FUNCTION
!!  Low-level tools for FFT (sequential and MPI parallel version)
!!  It also provides helper functions to set up the list of G vectors
!!  inside a sphere or to count them.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2019 ABINIT group (SG, XG, AR, MG, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! TODO
!!  1) Pass distribfft instead of MPI_enreg to simplify the API and facilitate code-reuse.
!!
!!  2) Merge this module with m_distribfft
!!
!!  3) Get rid of paral_kgb and MPI_type! This is a low-level module that may be called by other
!!     code in which paral_kgb is meaningless! FFT tables and a MPI communicator are sufficient.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_fftcore

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_sort

 use m_time,         only : timab
 use m_fstrings,     only : itoa, sjoin
 use defs_abitypes,  only : MPI_type
 use m_mpinfo,       only : destroy_mpi_enreg, initmpi_seq

 implicit none

 private

 public :: fftalg_isavailable  ! True if the FFT library specified by fftalg is available.
 public :: fftalg_has_mpi      ! True if fftalg provides MPI-FFTs.
 public :: fftalg_for_npfft    ! Returns the default value for fftalg given the number of processors for the FFT.
 public :: fftalg_info         ! Returns strings with info on the FFT library specified by fftalg.
 public :: get_cache_kb        ! Returns the cache size in Kbs (based on CPP variables).
 public :: ngfft_seq           ! initialize ngfft(18) from the FFT divisions (assume sequential FFT)
 public :: print_ngfft         ! Print the content of ngfft(18) in explicative format.
 public :: bound               ! Find distance**2 to boundary point of fft box nearest to kpt
 public :: getng               ! From ecut and metric tensor in reciprocal space, computes recommended ngfft(1:3)
 public :: sphereboundary      ! Finds the boundary of the basis sphere of G vectors
 public :: sphere
 public :: sphere_fft      ! Insert cg inside box.
 public :: sphere_fft1     ! TODO: This one should be replaced by sphere_fft.
 public :: change_istwfk   ! Change the istwfk mode of a set of wavefunctions (sequential version, same k-point)
 public :: kpgsph          ! Set up the G vector list
 public :: kpgcount        ! Give the number of G vector in each direction
 public :: get_kg          ! Helper function to calculate the set of G-vectors at a given kpoint (no MPI FFT)
 public :: kgindex         ! Compute the index of each plane wave on a FFT grid.

 ! Low-level tools for MPI FFT
 public :: switch
 public :: switch_cent
 public :: switchreal
 public :: switchreal_cent
 public :: scramble
 public :: fill
 public :: fill_cent
 public :: unfill
 public :: unfill_cent
 public :: unmpiswitch
 public :: unswitch
 public :: unswitchreal_cent
 public :: unmpiswitch_cent
 public :: unscramble
 public :: unswitch_cent
 public :: unswitchreal
 public :: mpiswitch
 public :: mpiswitch_cent

 public :: mpifft_fg2dbox
 public :: mpifft_fg2dbox_dpc
 public :: mpifft_dbox2fg
 public :: mpifft_dbox2fg_dpc
 public :: mpifft_dbox2fr
 public :: mpifft_dbox2fr_dpc
 public :: mpifft_fr2dbox
 public :: mpifft_fr2dbox_dpc
 public :: mpifft_collect_datar           ! Collect a real-space MPI-FFT distributed array on each proc.

 public :: indfftrisc
 public :: addrho
 public :: multpot

 ! 0 for double precision version (default), 1 for mixed precision FFTs
 integer, public, protected :: fftcore_mixprec = 0
 public :: fftcore_set_mixprec

! *************************************************************************

!----------------------------------------------------------------------
! Private variables

#define FFTALGA_SIZE 5
 character(len=*),private,parameter :: fftalga2name(1:FFTALGA_SIZE)= &
& (/"Goedecker     ", &
&   "Vendor FFT    ", &
&   "FFTW3         ", &
&   "Goedecker2002 ", &
&   "DFTI          " /)

#define FFTALGB_SIZE 1
 character(len=*),private,parameter :: fftalgb2name(0:FFTALGB_SIZE)= &
& (/"C2C",&
&   "R2C"/)

#define FFTALGC_SIZE 2
 character(len=*),private,parameter :: fftalgc2name(0:FFTALGC_SIZE)= &
& (/"No pad         ",&
&   "zero-pad       ",&
&   "zero-pad+cache "/)

contains
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/fftcore_set_mixprec
!! NAME
!! fftalg_set_precision
!!
!! FUNCTION
!!  Set the precision to be used in the FFT routines: 0 for standard double precision,
!!  1 for mixed precision (dp input, sp for intermediate arrays passed to FFT libs)
!!  Return old value.
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

integer function fftcore_set_mixprec(wp) result(old_wp)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: wp

! *************************************************************************

 old_wp = fftcore_mixprec
 fftcore_mixprec = abs(wp)

 select case (fftcore_mixprec)
 case (0)
   if (old_wp /= fftcore_mixprec) call wrtout(std_out, " Activating FFT in double-precision.")
 case (1)
   if (old_wp /= fftcore_mixprec) call wrtout(std_out, " Activating FFT in mixed precision.")
 case default
   MSG_ERROR(sjoin("Wrong value for input wp:", itoa(fftcore_mixprec)))
 end select

end function fftcore_set_mixprec
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/fftalg_isavailable
!! NAME
!! fftalg_isavailable
!!
!! FUNCTION
!!  Returns TRUE if the FFT library specified by fftalg (ngfft(7)) is available
!!
!! INPUTS
!!  fftalg=Input variable.
!!
!! PARENTS
!!
!! SOURCE

pure function fftalg_isavailable(fftalg) result(ans)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalg
 logical :: ans

!Local variables-------------------------------
!scalars
 integer :: fftalga,fftalgb,fftalgc

! *************************************************************************

 ans = .TRUE.
 fftalga = fftalg/100
 fftalgb = mod(fftalg,100)/10
 fftalgc = mod(fftalg,10)

 ! Optional FFT libraries.
#ifndef HAVE_FFT_FFTW3
 if (fftalga == FFT_FFTW3) ans = .FALSE.
#endif

#ifndef HAVE_FFT_DFTI
 if (fftalga == FFT_DFTI) ans = .FALSE.
#endif

end function fftalg_isavailable
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/fftalg_has_mpi
!! NAME
!! fftalg_has_mpi
!!
!! FUNCTION
!!  True if the FFT library specified by fftalg is available.
!!
!! INPUTS
!!  fftalg=Input variable.
!!
!! PARENTS
!!
!! SOURCE

pure function fftalg_has_mpi(fftalg) result(ans)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalg
 logical :: ans

!Local variables-------------------------------
!scalars
 integer :: fftalga,fftalgb,fftalgc

! *************************************************************************

 ans = .False.
 fftalga = fftalg/100; fftalgb = mod(fftalg,100)/10; fftalgc = mod(fftalg,10)

 if (fftalga == FFT_FFTW3) ans = .True.
 !if (fftalga == FFT_DFTI)  ans = .True.
 if (fftalga == FFT_SG2002) ans = .True.

end function fftalg_has_mpi
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/fftalg_for_npfft
!! NAME
!! fftalg_for_npfft
!!
!! FUNCTION
!!  Returns the default value of fftalg given the number of MPI nodes
!!  to be used in the FFTs.
!!
!! INPUTS
!!  nproc_fft=Number of processors used for MPI FFT
!!
!! OUTPUT
!!  fftalg=Integer used to select the FFT library.
!!
!! PARENTS
!!
!! SOURCE

pure function fftalg_for_npfft(nproc_fft) result(fftalg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nproc_fft
 integer :: fftalg

! *************************************************************************

 ! Default  for the sequential case.
 fftalg = 112

 ! Use Goedecker2002 if fftalg does not support MPI (e.g 112)
 if (nproc_fft > 1) fftalg = 401
 !if (nproc_fft > 1) fftalg = 402

#ifdef HAVE_FFT_FFTW3
 fftalg = 312
#elif defined HAVE_FFT_DFTI
 fftalg = 512
 if (nproc_fft > 1) fftalg = 401  ! MPI-FFT with DFTI is not implemented yet
#endif

 !if (nproc_fft > 1) fftalg = 401 ! This is to revert to the old behavior.

end function fftalg_for_npfft
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/fftalg_info
!! NAME
!! fftalg_info
!!
!! FUNCTION
!!  Returns info on the FFT library specified by fftalg (ngfft(7))
!!
!! INPUTS
!!  fftalg=Input variable.
!!
!! OUTPUT
!!  library=String with the name of FFT library
!!  cplex_mode= String defining whether the FFT library supports real<-->complex transforms.
!!  padding_mode=Padding mode.
!!
!! PARENTS
!!      m_fft,m_fft_prof
!!
!! CHILDREN
!!      xmpi_sum,xmpi_sum_master
!!
!! SOURCE

subroutine fftalg_info(fftalg,library,cplex_mode,padding_mode)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalg
 character(len=*),intent(out) :: library,cplex_mode,padding_mode

!Local variables-------------------------------
!scalars
 integer :: fftalga,fftalgb,fftalgc

! *************************************************************************

 library = "Unknown"; cplex_mode = "Unknown"; padding_mode  = "Unknown"

 fftalga=fftalg/100
 if (fftalga>0 .and. fftalga<=FFTALGA_SIZE) library = fftalga2name(fftalga)

 fftalgb=mod(fftalg,100)/10
 if (fftalgb>=0 .and. fftalgb<=FFTALGB_SIZE) cplex_mode = fftalgb2name(fftalgb)

 fftalgc=mod(fftalg,10)
 if (fftalgc>=0 .and. fftalgc<=FFTALGC_SIZE) padding_mode = fftalgc2name(fftalgc)

end subroutine fftalg_info
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/get_cache_kb
!! NAME
!! get_cache_kb
!!
!! FUNCTION
!!  Returns the cache size in KB to be used for cache blocking algorithms in the FFT routines.
!!  The value is derived from the values of the CPP options defined in config.h
!!
!! TODO
!!   Use C to get the real cache size.
!!   See http://stackoverflow.com/questions/12594208/c-program-to-determine-levels-size-of-cache
!!
!! PARENTS
!!
!! SOURCE

pure function get_cache_kb()

!Local variables-------------------------------
!scalars
 integer :: get_cache_kb

! *************************************************************************

 ! Default value
 get_cache_kb = 16
 !get_cache_kb = 32
 !get_cache_kb = 256

end function get_cache_kb
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/ngfft_seq
!! NAME
!! ngfft_seq
!!
!! FUNCTION
!! Helper function used to initialize ngfft(18) from the FFT divisions
!! in the case of sequential execution.
!!
!! INPUTS
!!  n123(3)=FFT divisions.
!!
!! OUTPUT
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft.
!!
!! PARENTS
!!      m_dvdb,m_phgamma,m_wfk,mrgdv
!!
!! CHILDREN
!!      xmpi_sum,xmpi_sum_master
!!
!! SOURCE

pure subroutine ngfft_seq(ngfft, n123)

!Arguments ------------------------------------
!scalars
!arrays
 integer,intent(in) :: n123(3)
 integer,intent(out) :: ngfft(18)

!Local variables-------------------------------
!scalars
 integer :: fftalg

! *************************************************************************

 ! Default  for the sequential case.
 fftalg = 112
#ifdef HAVE_FFT_FFTW3
 fftalg = 312
#elif defined HAVE_FFT_DFTI
 fftalg = 512
#endif

 ngfft(1:3) = n123
 ngfft(4) = 2*(ngfft(1)/2)+1
 ngfft(5) = 2*(ngfft(2)/2)+1
 ngfft(6) = ngfft(3)
 ngfft(7)= fftalg           ! fftalg
 ngfft(8)= get_cache_kb()   ! cache_kb
 ngfft(9)= 0                ! paral_fft_
 ngfft(10)=1                ! nproc_fft
 ngfft(11)=0                ! me_fft
 ngfft(12)=0                ! n2proc
 ngfft(13)=0                ! n3proc
 ngfft(14:18)=0             ! not used

end subroutine ngfft_seq
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/print_ngfft
!! NAME
!! print_ngfft
!!
!! FUNCTION
!!  Print the content of ngfft(18) in explicative format.
!!
!! INPUTS
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft.
!!  [unit]=unit number for output (defaults to std_out).
!!  [prtvol]=verbosity level (defaults to 0).
!!  [mode_paral]=either "COLL" or "PERS" ("COLL" is default).
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      bethe_salpeter,eph,getng,m_fft,m_fft_prof,m_wfd,screening,setup_bse
!!      sigma,wfk_analyze
!!
!! CHILDREN
!!      xmpi_sum,xmpi_sum_master
!!
!! SOURCE

subroutine print_ngfft(ngfft,header,unit,mode_paral,prtvol)

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: prtvol,unit
 character(len=*),intent(in),optional :: header
 character(len=4),intent(in),optional :: mode_paral
!arrays
 integer,intent(in) :: ngfft(18)

!Local variables-------------------------------
!scalars
 integer :: my_unt,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_prtvol=0;       if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_mode  ='COLL';  if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=ch10//' ==== FFT mesh description (ngfft) ==== '
 if (PRESENT(header)) msg=ch10//' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)
 write(msg,'(2(a,3i5,a),a,i5,2a,i5)')&
&  '  FFT mesh divisions ........................ ',ngfft(1),ngfft(2),ngfft(3),ch10,&
&  '  Augmented FFT divisions ................... ',ngfft(4),ngfft(5),ngfft(6),ch10,&
&  '  FFT algorithm ............................. ',ngfft(7),ch10,&
&  '  FFT cache size ............................ ',ngfft(8)
 call wrtout(my_unt,msg,my_mode)

 if (my_prtvol>0) then
   write(msg,'(6(a,i5,a),a,4i5)')&
&    '  FFT parallelization level ................. ',ngfft(9),ch10,&
&    '  Number of processors in my FFT group ...... ',ngfft(10),ch10,&
&    '  Index of me in my FFT group ............... ',ngfft(11),ch10,&
&    '  No of xy planes in R space treated by me .. ',ngfft(12),ch10,&
&    '  No of xy planes in G space treated by me .. ',ngfft(13),ch10,&
&    '  MPI communicator for FFT .................. ',ngfft(14),ch10,&
&    '  Value of ngfft(15:18) ..................... ',ngfft(15:18)
   call wrtout(my_unt,msg,my_mode)
 end if

end subroutine print_ngfft
!!***

!!****f* m_fftcore/bound
!! NAME
!! bound
!!
!! FUNCTION
!! For given kpt, ngfft, and gmet,
!!  Find distance**2 to boundary point of fft box nearest to kpt
!!  Find distance**2 to boundary point of fft box farthest to kpt
!!
!! INPUTS
!!  kpt(3)=real input k vector (reduced coordinates)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  gmet(3,3)=reciprocal space metric (currently in Bohr**-2)
!!
!! OUTPUT
!!  dsqmax=maximum distance**2 from k to boundary in Bohr**-2.
!!  dsqmin=minimum distance**2 from k to boundary in Bohr**-2.
!!  gbound(3)=coords of G on boundary (correspnding to gsqmin)
!!  plane=which plane min occurs in (1,2, or 3 for G1,etc).
!!
!! NOTES
!! Potential trouble: this routine was written assuming kpt lies inside
!! first Brillouin zone.  No measure is taken to fold input kpt back
!! into first zone.  Given arbitrary kpt, this will cause trouble.
!!
!! PARENTS
!!      getcut,getng
!!
!! CHILDREN
!!
!! SOURCE

subroutine bound(dsqmax,dsqmin,gbound,gmet,kpt,ngfft,plane)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: plane
 real(dp),intent(out) :: dsqmax,dsqmin
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(out) :: gbound(3)
 real(dp),intent(in) :: gmet(3,3),kpt(3)

!Local variables-------------------------------
!scalars
 integer :: i1,i1min,i2,i2min,i3,i3min
 real(dp) :: dsm,dsp
 character(len=500) :: message

! *************************************************************************

! dsq(i1,i2,i3)=gmet(1,1)*(kpt(1)+dble(i1))**2&
!& +gmet(2,2)*(kpt(2)+dble(i2))**2&
!& +gmet(3,3)*(kpt(3)+dble(i3))**2&
!& +2._dp*(gmet(1,2)*(kpt(1)+dble(i1))*(kpt(2)+dble(i2))&
!& +gmet(2,3)*(kpt(2)+dble(i2))*(kpt(3)+dble(i3))&
!& +gmet(3,1)*(kpt(3)+dble(i3))*(kpt(1)+dble(i1)))

!Set plane to impossible value
 plane=0

!look at +/- g1 planes:
 dsqmax=zero
 dsqmin=dsq(ngfft(1)/2,-ngfft(2)/2,-ngfft(3)/2,gmet,kpt)+0.01_dp
 do i2=-ngfft(2)/2,ngfft(2)/2
   do i3=-ngfft(3)/2,ngfft(3)/2
     dsp = dsq(ngfft(1)/2, i2, i3,gmet,kpt)
     dsm = dsq( - ngfft(1)/2, i2, i3,gmet,kpt)
     if (dsp>dsqmax) dsqmax = dsp
     if (dsm>dsqmax) dsqmax = dsm
     if (dsp<dsqmin) then
       dsqmin = dsp
       i1min = ngfft(1)/2
       i2min = i2
       i3min = i3
       plane=1
     end if
     if (dsm<dsqmin) then
       dsqmin = dsm
       i1min =  - ngfft(1)/2
       i2min = i2
       i3min = i3
       plane=1
     end if
   end do
 end do
!
!+/- g2 planes:
 do i1=-ngfft(1)/2,ngfft(1)/2
   do i3=-ngfft(3)/2,ngfft(3)/2
     dsp = dsq(i1,ngfft(2)/2,i3,gmet,kpt)
     dsm = dsq(i1,-ngfft(2)/2,i3,gmet,kpt)
     if (dsp>dsqmax) dsqmax = dsp
     if (dsm>dsqmax) dsqmax = dsm
     if (dsp<dsqmin) then
       dsqmin = dsp
       i1min = i1
       i2min = ngfft(2)/2
       i3min = i3
       plane=2
     end if
     if (dsm<dsqmin) then
       dsqmin = dsm
       i1min = i1
       i2min =  - ngfft(2)/2
       i3min = i3
       plane=2
     end if
   end do
 end do
!
!+/- g3 planes:
 do i1=-ngfft(1)/2,ngfft(1)/2
   do i2=-ngfft(2)/2,ngfft(2)/2
     dsp = dsq(i1,i2,ngfft(3)/2,gmet,kpt)
     dsm = dsq(i1,i2,-ngfft(3)/2,gmet,kpt)
     if (dsp>dsqmax) dsqmax = dsp
     if (dsm>dsqmax) dsqmax = dsm
     if (dsp<dsqmin) then
       dsqmin = dsp
       i1min = i1
       i2min = i2
       i3min = ngfft(3)/2
       plane=3
     end if
     if (dsm<dsqmin) then
       dsqmin = dsm
       i1min = i1
       i2min = i2
       i3min =  - ngfft(3)/2
       plane=3
     end if
   end do
 end do

 if (plane==0) then
!  Trouble: missed boundary somehow
   write(message, '(a,a,a,3f9.4,a,3i5,a,a,a,a,a)' )&
&   'Trouble finding boundary of G sphere for',ch10,&
&   'kpt=',kpt(:),' and ng=',ngfft(1:3),ch10,&
&   'Action : check that kpt lies',&
&   'reasonably within first Brillouin zone; ',ch10,&
&   'else code bug, contact ABINIT group.'
   MSG_BUG(message)
 end if

 gbound(1)=i1min
 gbound(2)=i2min
 gbound(3)=i3min

 contains

   function dsq(i1,i2,i3,gmet,kpt)

     integer :: i1,i2,i3
     real(dp) :: dsq
     real(dp) :: kpt(3),gmet(3,3)

     dsq=gmet(1,1)*(kpt(1)+dble(i1))**2&
&      +gmet(2,2)*(kpt(2)+dble(i2))**2&
&      +gmet(3,3)*(kpt(3)+dble(i3))**2&
&      +2._dp*(gmet(1,2)*(kpt(1)+dble(i1))*(kpt(2)+dble(i2))&
&      +gmet(2,3)*(kpt(2)+dble(i2))*(kpt(3)+dble(i3))&
&      +gmet(3,1)*(kpt(3)+dble(i3))*(kpt(1)+dble(i1)))
   end function dsq

end subroutine bound
!!***

!!****f* m_fftcore/getng
!! NAME
!! getng
!!
!! FUNCTION
!! From ecut and metric tensor in reciprocal space, computes recommended ngfft(1:3)
!! Also computes the recommended value of nfft and mgfft
!! Pay attention that the FFT grid must be compatible with the symmetry operations (see irrzg.f).
!!
!! INPUTS
!! boxcutmin=minimum value of boxcut admitted (boxcut is the ratio
!!  between the radius of the sphere contained in the FFT box, and the
!!  radius of the planewave sphere) : usually 2.0 .
!! ecut=energy cutoff in Hartrees
!! gmet(3,3)=reciprocal space metric (bohr**-2).
!! kpt(3)=input k vector in terms of reciprocal lattice primitive translations
!! me_fft=index of the processor in the FFT set (use 0 if sequential)
!! nproc_fft=number of processors in the FFT set (use 1 if sequential)
!! nsym=number of symmetry elements in group
!! paral_fft=0 if no FFT parallelisation ; 1 if FFT parallelisation
!! symrel(3,3,nsym)=symmetry matrices in real space (integers)
!!
!! OUTPUT
!! mgfft= max(ngfft(1),ngfft(2),ngfft(3))
!! nfft=number of points in the FFT box=ngfft(1)*ngfft(2)*ngfft(3)/nproc_fft
!!
!! SIDE EFFECTS
!! Input/Output
!! ngfft(1:18)=integer array with FFT box dimensions and other
!!   information on FFTs. On input ngfft(1:3) contains
!!   optional trial values. If ngfft(1:3)/minbox is greater than value
!!   calculated to avoid wrap-around error and ngfft obeys constraint
!!   placed by the FFT routine that is used
!!   then ngfft(1:3) is left unchanged. Otherwise set to value computed in now.
!!
!! Note that there is the possibility of an undetected error if we
!! are dealing with a cubic unit cell and ngfft(1), ngfft(2) and ngfft(3)
!! are different. In the future we should handle this case.
!!
!! ngfft(4),ngfft(5),ngfft(6)= modified values to avoid cache trashing,
!!        presently : ngfft(4)=ngfft(1)+1 if ngfft(1) is even ;
!!                    ngfft(5)=ngfft(2)+1 if ngfft(2) is even.
!!           in the other cases, ngfft(4:6)=ngfft(1:3).
!!   Other choices may better, but this is left for the future.
!! ngfft(7)=choice for FFT algorithm, see the input variable fftalg
!! ngfft(8)=size of the cache, in bytes (not used here presently).!!
!!   other ngfft slots are used for parallelism see ~abinit/doc/variables/vargs.htm#ngfft
!! [unit] = Output Unit number (DEFAULT std_out)
!!
!! PARENTS
!!      fftprof,m_fft,m_fft_prof,memory_eval,mpi_setup,mrgscr,scfcv
!!
!! CHILDREN
!!      bound,print_ngfft,sort_int,wrtout
!!
!! SOURCE

subroutine getng(boxcutmin,ecut,gmet,kpt,me_fft,mgfft,nfft,ngfft,nproc_fft,nsym,paral_fft,symrel,&
&                ngfftc,use_gpu_cuda,unit) ! optional

 use defs_fftdata,  only : mg

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: me_fft,nproc_fft,nsym,paral_fft
 integer,intent(out) :: mgfft,nfft
 integer,optional,intent(in) :: unit,use_gpu_cuda
 real(dp),intent(in) :: boxcutmin,ecut
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 integer,intent(in),optional :: ngfftc(3)
 integer,intent(inout) :: ngfft(18)
 real(dp),intent(in) :: gmet(3,3),kpt(3)

!Local variables-------------------------------
!scalars
 integer,save :: first=1,msrch(3),previous_paral_mode=0
 integer :: element,ii,index,isrch,isrch1,isrch2,isrch3,isym,jj,mu,paral_fft_
 integer :: plane,testok,tobechecked,ount,fftalga
 real(dp),parameter :: minbox=0.75_dp
 real(dp) :: dsqmax,dsqmin,ecutmx,prodcurrent,prodtrial,xx,yy
 logical :: testdiv
 character(len=500) :: message
 integer,parameter :: largest_ngfft=mg ! Goedecker FFT: any powers of 2, 3, and 5 - must be coherent with defs_fftdata.F90
 integer,parameter :: maxpow2 =16      ! int(log(largest_ngfft+half)/log(two))
 integer,parameter :: maxpow3 =6       ! int(log(largest_ngfft+half)/log(three))
 integer,parameter :: maxpow5 =6       ! int(log(largest_ngfft+half)/log(five))
 integer,parameter :: maxpow7 =0
 integer,parameter :: maxpow11=0
 integer,parameter :: mmsrch=(maxpow2+1)*(maxpow3+1)*(maxpow5+1)*(maxpow7+1)*(maxpow11+1)
 integer,save :: iperm(mmsrch),srch(mmsrch,3)
 integer(i8b) :: li_srch(mmsrch)
 integer :: divisor(3,3),gbound(3),imax(3),imin(3),ngcurrent(3)
 integer :: ngmax(3),ngsav(3),ngtrial(3)

! *************************************************************************

 ount = std_out; if (PRESENT(unit)) ount = unit

 fftalga = ngfft(7)/100

!If not yet done, compute recommended (boxcut>=2) fft grid dimensions
!In case we switch for paral to sequential mode, recompute srch.
!This is the case e.g. when computing ngfftdiel in sequential mode
!after an initial computation of ngfft in parallel

 paral_fft_=paral_fft;if (nproc_fft==0) paral_fft_=0

 if(first==1.or.paral_fft_ /= previous_paral_mode) then
   first=0; previous_paral_mode=paral_fft_
   srch(:,:)=0

   ! Factors of 2
   srch(1,1)=1
   do ii=1,maxpow2
     srch(ii+1,1)=srch(ii,1)*2
   end do

   ! Factors of 3
   index=maxpow2+1
   if(maxpow3>0)then
     do ii=1,max(1,maxpow3)
       srch(1+ii*index:(ii+1)*index,1)=3*srch(1+(ii-1)*index:ii*index,1)
     end do
   end if

   ! Factors of 5
   index=(maxpow3+1)*index
   if(maxpow5>0)then
     do ii=1,max(1,maxpow5)
       li_srch = 0
       li_srch(1+ii*index:(ii+1)*index)=5_i8b*srch(1+(ii-1)*index:ii*index,1)
       where (li_srch > huge(maxpow3)) li_srch = huge(maxpow3)
       srch(1+ii*index:(ii+1)*index,1)=li_srch(1+ii*index:(ii+1)*index)
     end do
   end if

   ! Factors of 7
   index=(maxpow5+1)*index
   if(maxpow7>0)then
     do ii=1,max(1,maxpow7)
       srch(1+ii*index:(ii+1)*index,1)=7*srch(1+(ii-1)*index:ii*index,1)
     end do
   end if

   ! Factors of 11
   index=(maxpow7+1)*index
   if(maxpow11>0)then
     do ii=1,max(1,maxpow11)
       srch(1+ii*index:(ii+1)*index,1)=11*srch(1+(ii-1)*index:ii*index,1)
     end do
   end if

   call sort_int(mmsrch,srch(:,1),iperm)

   do ii=1,mmsrch
     if(srch(ii,1)>largest_ngfft)exit
   end do
   msrch(1)=ii-1

   ! In case of FFT parallelism, one need ngfft(2) and ngfft(3) to be multiple of nproc_fft
   if(paral_fft_==1)then
     msrch(2)=0
     do ii=1,msrch(1)
       if(modulo(srch(ii,1),nproc_fft)==0) then
         msrch(2)=msrch(2)+1
         srch(msrch(2),2)=srch(ii,1)
       end if
     end do
     !write(message,'(a,i0,a,i0,2a,i0)')&
     ! 'The second and third dimension of the FFT grid: ',ngfft(2),", ",ngfft(3),ch10,&
     ! 'were imposed to be multiple of the number of processors for the FFT: ', nproc_fft
     !if (ount /= dev_null) MSG_COMMENT(message)
   else
     msrch(2)=msrch(1)
     srch(:,2)=srch(:,1)
   end if

   ! The second and third search list have the same constraint
   msrch(3)=msrch(2)
   srch(:,3)=srch(:,2)

!  The set of allowed ngfft values has been found
 end if ! first==1

!Save input values of ngfft
 ngsav(1:3) = ngfft(1:3)

!As an initial guess for ngfft, use the provided coarse mesh grid
 if (PRESENT(ngfftc)) then
   ngfft(1:3)=ngfftc(1:3)
   call wrtout(ount,' Using supplied coarse mesh as initial guess.','COLL')
 else
   ngfft(1:3)=2
 end if

!Enlarge the initial guess until the set of ngfft entirely comprises the sphere
 do

   call bound(dsqmax,dsqmin,gbound,gmet,kpt,ngfft,plane)

   ! Exit the infinite do-loop if the sphere is inside the FFT box
   if (dsqmin>=(half*boxcutmin**2*ecut/pi**2)) exit

   ! Fix nearest boundary
   do ii=1,msrch(plane)-1
     if (srch(ii,plane)>=ngfft(plane)) then
!      redefine ngfft(plane) to next higher choice
       ngfft(plane)=srch(ii+1,plane)
       exit ! Exit the loop over ii
     end if

     if (ii==msrch(plane)-1)then
       ! Here, we are in trouble
       write(message, '(a,i12,5a)' ) &
&       'ngfft is bigger than allowed value =',ngfft(plane),'.',ch10,&
&       'This indicates that desired ngfft is larger than getng',ch10,&
&       'can handle. The code has to be changed and compiled.'
       MSG_BUG(message)
     end if
   end do

 end do ! End of the infinite do-loop : will either "exit", or abort

!ecutmx=maximum ecut consistent with chosen ngfft
 ecutmx=0.5_dp*pi**2*dsqmin

!Print results
 write(message, '(a,1p,e14.6,a,3i8,a,a,e14.6)' ) &
& ' For input ecut=',ecut,' best grid ngfft=',ngfft(1:3),ch10,&
& '       max ecut=',ecutmx
 call wrtout(ount,message,'COLL')

! The FFT grid is compatible with the symmetries if for each
! symmetry isym, each ii and each jj, the quantity
! (ngfft(jj)*symrel(jj,ii,isym))/ngfft(ii) is an integer.
! This relation is immediately verified for diagonal elements, since
! symrel is an integer. It is also verified if symrel(ii,jj,isym) is zero.

!Compute the biggest (positive) common divisor of each off-diagonal element of the symmetry matrices
 divisor(:,:)=0; tobechecked=0

 do ii=1,3
   do jj=1,3
     if(jj==ii)cycle
     do isym=1,nsym
       if(symrel(jj,ii,isym)==0 .or. divisor(jj,ii)==1 )cycle
       tobechecked=1
       element=abs(symrel(jj,ii,isym))
       testdiv= ( divisor(jj,ii)==0 .or. divisor(jj,ii)==element .or. element==1)
       if(testdiv)then
         divisor(jj,ii)=element
       else
!        Must evaluate common divisor between non-trivial numbers
         do
           if(divisor(jj,ii)<element)element=element-divisor(jj,ii)
           if(divisor(jj,ii)>element)divisor(jj,ii)=divisor(jj,ii)-element
           if(divisor(jj,ii)==element)exit
         end do
       end if
     end do
   end do
 end do

!Check whether there is a problem
 testok=1
 if(tobechecked==1)then
   do ii=1,3
     do jj=1,3
       xx=divisor(jj,ii)*ngfft(jj)
       yy=xx/ngfft(ii)
       if(abs(yy-nint(yy))>tol8)testok=0
     end do
   end do
 end if

!There is definitely a problem
 if(testok==0)then
!  Use a brute force algorithm
!  1) Because one knows that three identical numbers will satisfy
!  the constraint, use the maximal ngfft value to define current triplet
!  and associate total number of grid points
   ngcurrent(1:3)=maxval(ngfft(1:3))
!  Takes into account the fact that ngfft(2) and ngfft(3) must
!  be multiple of nproc_fft
   if(mod(ngcurrent(1),nproc_fft)/=0)ngcurrent(1:3)=ngcurrent(1:3)*max(1,nproc_fft)
   prodcurrent=ngcurrent(1)**3+1.0d-3
!  2) Define maximal values for each component, limited
!  by the maximal value of the list
   ngmax(1)=min(int(prodcurrent/(ngfft(2)*ngfft(3))),srch(msrch(1),1))
   ngmax(2)=min(int(prodcurrent/(ngfft(1)*ngfft(3))),srch(msrch(2),2))
   ngmax(3)=min(int(prodcurrent/(ngfft(1)*ngfft(2))),srch(msrch(3),3))
!  3) Get minimal and maximal search indices
   do ii=1,3
     do isrch=1,msrch(ii)
       index=srch(isrch,ii)
       if(index==ngfft(ii))imin(ii)=isrch
!      One cannot suppose that imax belongs to the allowed list,
!      so must use <= instead of == , to determine largest index
       if(index<=ngmax(ii))imax(ii)=isrch
     end do
   end do
!  4) Compute product of trial ngffts
!  DEBUG
!  write(ount,*)' getng : enter triple loop '
!  write(ount,*)'imin',imin(1:3)
!  write(ount,*)'imax',imax(1:3)
!  write(ount,*)'ngcurrent',ngcurrent(1:3)
!  ENDDEBUG
   do isrch1=imin(1),imax(1)
     ngtrial(1)=srch(isrch1,1)
     do isrch2=imin(2),imax(2)
       ngtrial(2)=srch(isrch2,2)
       do isrch3=imin(3),imax(3)
         ngtrial(3)=srch(isrch3,3)
         prodtrial=real(ngtrial(1))*real(ngtrial(2))*real(ngtrial(3))+1.0d-3
         if(prodtrial>prodcurrent-1.0d-4)exit
!        The trial product is lower or equal to the current product,
!        so now, checks whether the symmetry constraints are OK
         testok=1
         do ii=1,3
           do jj=1,3
             xx=divisor(jj,ii)*ngtrial(jj)
             yy=xx/ngtrial(ii)
             if(abs(yy-nint(yy))>tol8)testok=0
           end do
         end do
!        DEBUG
!        write(ount,'(a,3i6,a,i3,a,es16.6)' )' getng : current trial triplet',ngtrial(1:3),&
!        &     ' testok=',testok,' prodtrial=',prodtrial
!        ENDDEBUG
         if(testok==0)cycle
!        The symmetry constraints are fulfilled, so update current values
         ngcurrent(1:3)=ngtrial(1:3)
         prodcurrent=prodtrial
       end do
     end do
   end do

   ngfft(1:3)=ngcurrent(1:3)
   call bound(dsqmax,dsqmin,gbound,gmet,kpt,ngfft,plane)
!  ecutmx=maximum ecut consistent with chosen ngfft
   ecutmx=0.5_dp*pi**2*dsqmin
!  Give results
   write(message, '(a,3i8,a,a,e14.6)' ) &
&   ' However, must be changed due to symmetry =>',ngfft(1:3),ch10,&
&   '       with max ecut=',ecutmx
   call wrtout(ount,message,'COLL')

   if (prodcurrent>huge(ii)) then
     write(message, '(5a)' )&
&     'The best FFT grid will lead to indices larger',ch10,&
&     'than the largest representable integer on this machine.',ch10,&
&     'Action: try to deal with smaller problems. Also contact ABINIT group.'
     MSG_ERROR(message)
   end if

 end if ! testok==0

!Possibly use the input values of ngfft
 if ( int( dble(ngsav(1)) / minbox ) >= ngfft(1) .and.&
&     int( dble(ngsav(2)) / minbox ) >= ngfft(2) .and.&
&     int( dble(ngsav(3)) / minbox ) >= ngfft(3) ) then

!  Must check whether the values are in the allowed list
   testok=0
   do mu=1,3
     do ii=1,msrch(mu)
       if(srch(ii,mu)==ngsav(mu))then
         testok=testok+1
         exit
       end if
     end do
   end do
   if(testok==3)then
     write(message,'(a,3(a,i1,a,i3),a)') ' input values of',&
&     (' ngfft(',mu,') =',ngsav(mu),mu=1,3),' are alright and will be used'
     call wrtout(ount,message,'COLL')
     do mu = 1,3
       ngfft(mu) = ngsav(mu)
     end do
   end if

 end if

!mgfft needs to be set to the maximum of ngfft(1),ngfft(2),ngfft(3)
 mgfft = maxval(ngfft(1:3))

 if (paral_fft_==1) then
   ! For the time being, one need ngfft(2) and ngfft(3) to be multiple of nproc_fft
   if(modulo(ngfft(2),nproc_fft)/=0)then
     write(message,'(3a,i5,a,i5)')&
&     'The second dimension of the FFT grid, ngfft(2), should be',&
&     'a multiple of the number of processors for the FFT, nproc_fft.',&
&     'However, ngfft(2)=',ngfft(2),' and nproc_fft=',nproc_fft
     MSG_BUG(message)
   end if
   if(modulo(ngfft(3),nproc_fft)/=0)then
     write(message,'(3a,i5,a,i5)')&
&     'The third dimension of the FFT grid, ngfft(3), should be',&
&     'a multiple of the number of processors for the FFT, nproc_fft.',&
&     'However, ngfft(3)=',ngfft(3),' and nproc_fft=',nproc_fft
     MSG_BUG(message)
   end if

 else if (paral_fft_/=0) then
   write(message,'(a,i0)')'paral_fft_ should be 0 or 1, but its value is ',paral_fft_
   MSG_BUG(message)
 end if

! Compute effective number of FFT points (for this MPI node if parall FFT)
 nfft=product(ngfft(1:3))/max(1,nproc_fft)

!Set up fft array dimensions ngfft(4,5,6) to avoid cache conflicts
 ngfft(4)=2*(ngfft(1)/2)+1
 ngfft(5)=2*(ngfft(2)/2)+1
 ngfft(6)=ngfft(3)
 if (any(fftalga == [FFT_FFTW3, FFT_DFTI])) then
   ! FFTW3 supports leading dimensions but at the price of a larger number of FFTs
   ! to be executed along z when the zero-padded version is used.
   ! One should check whether the augmentation is beneficial for FFTW3.
   ngfft(4)=2*(ngfft(1)/2)+1
   ngfft(5)=2*(ngfft(2)/2)+1
   !ngfft(4)=ngfft(1)
   !ngfft(5)=ngfft(2)
   ngfft(6)=ngfft(3)
 end if

 if (present(use_gpu_cuda)) then
   if (use_gpu_cuda==1) then
     ngfft(4)=ngfft(1)
     ngfft(5)=ngfft(2)
     ngfft(6)=ngfft(3)
   end if
 end if

 ngfft(14:18)=0 ! ngfft(14) to be filled outside of getng

 if (paral_fft_==0) then
   ngfft(9)=0     ! paral_fft_
   ngfft(10)=1    ! nproc_fft
   ngfft(11)=0    ! me_fft
   ngfft(12)=0    ! n2proc
   ngfft(13)=0    ! n3proc
 else
   ngfft(9)=1     ! paral_fft_
   ngfft(10)=nproc_fft
   ngfft(11)=me_fft
   ngfft(12)=ngfft(2)/nproc_fft
   ngfft(13)=ngfft(3)/nproc_fft
 end if


 call print_ngfft(ngfft,header="FFT mesh",unit=ount,mode_paral="COLL")

end subroutine getng
!!***

!!****f* m_fftcore/sphereboundary
!! NAME
!! sphereboundary
!!
!! FUNCTION
!! Finds the boundary of the basis sphere of G vectors (at a given
!! k point) for use in improved zero padding of ffts in 3 dimensions.
!! Provides data to be used by subroutine fourwf, in the form of
!! an array gbound(2*mgfft+8,2).
!!
!! The first component (for use when mod(fftalg,10)==2))
!! provides integer values g1min,g1max,g2min,g2max
!! and then for g2 in the
!! sequence g2=0,1,2,...,g2max,g2min,g2min+1,...,-1, provides g1min, g1max.
!!
!! The second component (for use when mod(fftalg,10)==1))
!! provides integer values g1min,g1max,g3min,g3max,
!! where g3min and g3max have been corrected in case of time-reversal
!! and then for g3 in the sequence
!! g3=0,1,2,...,g3max,g3min,g3min+1,...,-1, provides g2min, g2max.
!! (also corrected in case of time-reversal)
!!
!! These values are stored in the above order in array gbound.
!! Debug mode, if fftalg is between 000 and 099
!!
!! INPUTS
!!  istwf_k=option parameter that describes the storage of wfs
!!  kg_k(3,npw)=integer coordinates of G vectors in basis sphere
!!  mgfft=maximum size of 1D FFTs (only for dimensioning purposes)
!!  npw=number of G vectors in basis at this k point
!!
!! OUTPUT
!!  gbound(2*mgfft+8,2)=defined above
!!
!! PARENTS
!!      dfpt_eltfrkin,dfpt_mkrho,fock_getghc,m_bandfft_kpt,m_cut3d,m_epjdos
!!      m_fft,m_fft_prof,m_fock,m_gsphere,m_hamiltonian,m_wfd,mkrho,mlwfovlp
!!      pawmkaewf,posdoppler,scfcv,spin_current,suscep_stat,susk,tddft,wfconv
!!
!! CHILDREN
!!
!! SOURCE

subroutine sphereboundary(gbound, istwf_k, kg_k, mgfft, npw)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k,mgfft,npw
!arrays
 integer,intent(in) :: kg_k(3,npw)
 integer,intent(out) :: gbound(2*mgfft+8,2)

!Local variables-------------------------------
!scalars
 integer :: dim_a,dim_b,fftalgc,g_a,gmax_a,gmax_b,gmax_b1,gmax_b2,gmin_a,gmin_b
 integer :: gmin_b1,gmin_b2,igb,ii,iloop,ipw,testm,testp,kgk
 character(len=500) :: message
!arrays
 integer :: gmax(2),gmin(2)

! *************************************************************************
!
!DEBUG
!write(std_out,*)' sphereboundary : enter'
!write(std_out, '(a)' )' sphereboundary : list of plane waves coordinates for k point '
!write(std_out, '(a)' )'       ipw       kg_k(1:3,ipw) '
!do ipw=1,npw
!write(std_out, '(i10,a,3i6)' )ipw,'  ',kg_k(1:3,ipw)
!end do
!gbound=-999
!ENDDEBUG

!Determine cube boundaries
 gbound(1,1)=minval(kg_k(1,:))
 gbound(2,1)=maxval(kg_k(1,:))
 gbound(1:2,2)=gbound(1:2,1)

!Treat differently the fftalgc cases
 do ii=1,2

   fftalgc=3-ii

   if(fftalgc/=2)then
     dim_a=3
     dim_b=2
   else
     dim_a=2
     dim_b=1
   end if

!  Relevant boundaries
   gbound(3,ii)=minval(kg_k(dim_a,:))
   gbound(4,ii)=maxval(kg_k(dim_a,:))
   gmin_a=gbound(3,ii)
   gmax_a=gbound(4,ii)

!  Must complete the sphere for fftalgc==1 and special storage modes.
!  Explanation : sg_fftpad is not able to take into account
!  the time-reversal symmetry, so that the boundaries will not be delimited
!  by the kg_k set, but by their symmetric also.
   if(istwf_k>=2 .and. fftalgc==1)then
     if( istwf_k==2 .or. istwf_k==3 .or. istwf_k==6 .or. istwf_k==7 )then
       gbound(4,2)=max(gmax_a,-gmin_a)
       gbound(3,2)=-gbound(4,2)
     else if( istwf_k==4 .or. istwf_k==5 .or. istwf_k==8 .or. istwf_k==9 )then
       gbound(4,2)=max(gmax_a,-gmin_a-1)
       gbound(3,2)=-gbound(4,2)-1
     end if
     gmax_a=gbound(4,2) ; gmin_a=gbound(3,2)
   end if

   igb=5

!  Consider first every g_a in range 0 ... gmax_a, then gmin_a ... -1
   gmin(1)=0         ; gmax(1)=gmax_a
   gmin(2)=gmin_a    ; gmax(2)=-1

   do iloop=1,2

     if( gmin(iloop) <= gmax(iloop) )then

       do g_a=gmin(iloop),gmax(iloop)

         if(istwf_k==1 .or. fftalgc/=1)then
           ! Select the minimal and maximal values, in the selected plane
           gmin_b=mgfft+1 ! Initialized with a value larger than all possible ones
           gmax_b=-mgfft-1 ! Initialized with a value smaller than all possible ones
           do ipw=1,npw
             if(kg_k(dim_a,ipw)==g_a)then
               kgk=kg_k(dim_b,ipw)
               if(kgk<=gmin_b)gmin_b=kgk
               if(kgk>=gmax_b)gmax_b=kgk
             end if
           end do

         else if(istwf_k>=2 .and. fftalgc==1)then

           ! Here, must take into account time-reversal symmetry explicitely

           ! Determine the boundaries for the plane g_a
           testp=0
           if(g_a<=gmax_a)then
             ! Select the minimal and maximal values, in the selected plane
             gmin_b1=mgfft+1 ! Initialized with a value larger than all possible ones
             gmax_b1=-mgfft-1 ! Initialized with a value smaller than all possible ones
             do ipw=1,npw
               if(kg_k(dim_a,ipw)==g_a)then
                 kgk=kg_k(dim_b,ipw)
                 if(kgk<=gmin_b1)gmin_b1=kgk
                 if(kgk>=gmax_b1)gmax_b1=kgk
               end if
             end do


             testp=1
           end if

           ! Determine the boundaries for the plane -g_a or -g_a-1
           testm=0
           if( istwf_k==2 .or. istwf_k==3 .or. istwf_k==6 .or. istwf_k==7 )then

             if(-g_a>=gmin_a)then
               ! Select the minimal and maximal values, in the selected plane
               ! Warning : there is an inversion of search (might be confusing)
               gmax_b2=mgfft+1 ! Initialized with a value larger than all possible ones
               gmin_b2=-mgfft-1 ! Initialized with a value smaller than all possible ones
               do ipw=1,npw
                 if(kg_k(dim_a,ipw)==-g_a)then
                   kgk=kg_k(dim_b,ipw)
                   if(kgk<=gmax_b2)gmax_b2=kgk
                   if(kgk>=gmin_b2)gmin_b2=kgk
                 end if
               end do
               testm=1
             end if

           else if( istwf_k==4 .or. istwf_k==5 .or. istwf_k==8 .or. istwf_k==9 )then

             if(-g_a-1>=gmin_a)then
               ! Select the minimal and maximal values, in the selected plane
               ! Warning : there is an inversion of search (might be confusing)
               gmax_b2=mgfft+1 ! Initialized with a value larger than all possible ones
               gmin_b2=-mgfft-1 ! Initialized with a value smaller than all possible ones
               do ipw=1,npw
                 if(kg_k(dim_a,ipw)==-g_a-1)then
                   kgk=kg_k(dim_b,ipw)
                   if(kgk<=gmax_b2)gmax_b2=kgk
                   if(kgk>=gmin_b2)gmin_b2=kgk
                 end if
               end do
               testm=1
             end if

           end if

           !  Must invert the boundaries, to use them for plane g_a
           if(testm==1)then
             ! This is needed to avoid border effect
             ! if the search did not lead to any element
             gmin_b2=max(gmin_b2,-mgfft) ; gmax_b2=min(gmax_b2,mgfft)
             if(istwf_k<=5)then
               gmax_b2=-gmax_b2 ; gmin_b2=-gmin_b2
             else
               gmax_b2=-gmax_b2-1 ; gmin_b2=-gmin_b2-1
             end if
           end if

           if( testp==1 .and. testm==1)then
             gmin_b=min(gmin_b1,gmin_b2) ; gmax_b=max(gmax_b1,gmax_b2)
           else if( testp==1 )then
             gmin_b=gmin_b1 ; gmax_b=gmax_b1
           else if( testm==1 )then
             gmin_b=gmin_b2 ; gmax_b=gmax_b2
           end if

         end if ! Endif take into account time-reversal symmetry

         if (igb+1>2*mgfft+4) then
           write(message, '(2a, 4(a,3(i0,1x),a))' )&
             "About to overwrite gbound array (FFT mesh too small) ",ch10, &
             "   iloop, igb, mgb = ",iloop,igb,2*mgfft+4, ch10, &
             "   istwfk, mgfft, npw = ",istwf_k, mgfft, npw, ch10, &
             "   minval(kg_k) = ",minval(kg_k, dim=2), ch10, &
             "   maxval(kg_k) = ",maxval(kg_k, dim=2), ch10
           MSG_BUG(message)
         end if

         gbound(igb,ii)=gmin_b
         gbound(igb+1,ii)=gmax_b

         if( iloop==1 .and. istwf_k>=2 .and. istwf_k<=5 .and. fftalgc==2 .and. g_a==0)then
!          If k_y=0 , for fftalgc==2, the g_a==0 plane must be completed
           if(istwf_k==2 .or. istwf_k==4)then
             gbound(igb+1,ii)=max(gmax_b,-gmin_b)
             gbound(igb,ii)=-gbound(igb+1,ii)
           else if(istwf_k==3 .or. istwf_k==5)then
             gbound(igb+1,ii)=max(gmax_b,-gmin_b-1)
             gbound(igb,ii)=-gbound(igb+1,ii)-1
           end if

         end if

         igb=igb+2

       end do ! g_a
     end if
   end do  ! iloop
 end do ! ii (fftalgc)

!DEBUG
!write(std_out,'(a)')' sphereoundary : list of plane waves coordinates for 1st k point '
!write(std_out,'(a)')'       ipw       kg_k(1:3,ipw) '
!do ipw=1,npw
!write(std_out, '(i10,a,3i6)' )ipw,'  ',kg_k(1:3,ipw)
!end do
!write(std_out, '(a)' )' sphereboundary : list of boundaries '
!do igb=1,2*mgfft+8
!write(std_out, '(i10,a,2i6)' )igb,'  ',gbound(igb,1),gbound(igb,2)
!end do
!write(std_out,*)' sphereboundary : exit '
!ENDDEBUG

end subroutine sphereboundary
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/sphere
!! NAME
!! sphere
!!
!! FUNCTION
!! Array cg is defined in sphere with npw points. Insert cg inside box
!! of n1*n2*n3 points to define array cfft for fft box.
!! corresponds to given element in cg.  rest of cfft is filled with 0 s.
!!
!! iflag=1==>insert cg into cfft.
!! iflag=2==>insert cg into cfft, where the second and third dimension
!! have been switched (needed for new 2002 SGoedecker FFT)
!! iflag=-1==> extract cg from cfft.
!! iflag=-2==> extract cg from cfft, where the second and third dimension
!! have been switched (needed for new 2002 SGoedecker FFT)
!!  (WARNING : iflag=-2 cannot use symmetry operations)
!!
!! There is also the possibility to apply a symmetry operation,
!! as well as to make a shift in reciprocal space, or to multiply
!! by a constant factor, in the case iflag=-1.
!! Multiplication by a constant factor is also possible in the case iflag=-2.
!!
!! INPUTS
!! cg(2,npw*ndat)= contains values for npw G vectors in basis sphere
!! ndat=number of FFT to do in //
!! npw=number of G vectors in basis at this k point
!! cfft(2,n4,n5,n6) = fft box
!! n1,n2,n3=physical dimension of the box (cfft)
!! n4,n5,n6=memory dimension of cfft
!! kg_k(3,npw)=integer coordinates of G vectors in basis sphere
!! istwf_k=option parameter that describes the storage of wfs
!! iflag=option parameter. Possible values: -1, -2, 1, 2
!! me_g0=1 if this node has G=0.
!! shiftg(3)=The shift in reciprocal space.
!! symm(3,3)=symmetry operation in reciprocal space to be applied.
!! xnorm=Normalization factor.
!!
!! SIDE EFFECTS
!! Input/Output
!! iflag=1 and 2, insert cg(input) into cfft(output)
!! iflag=-1 and -2, extract cg(output) from cfft(input)
!!
!! NOTES
!! cg and cfft are assumed to be of type COMPLEX, although this routine treats
!! them as real of twice the length to avoid nonstandard complex*16.
!! If istwf_k differs from 1, then special storage modes must be taken
!! into account, for symmetric wavefunctions coming from k=(0 0 0) or other
!! special k points.
!!
!! TODO
!! 1) Order arguments
!! 2) Split the two cases to avoid breaking intent: from and to sphere (merge with cg_box2gpsh and cg_gsph2box?)
!! 3) If symmetries are used with or without shiftg, it might happen that the FFT mesh
!!    is not large enough to accomodate the rotated G, in this case one should return ierr /= 0
!!
!! PARENTS
!!      cg_rotate,fourwf,m_fft,m_fftcore,m_fftw3,m_io_kss,wfconv
!!
!! CHILDREN
!!      xmpi_sum,xmpi_sum_master
!!
!! SOURCE

subroutine sphere(cg,ndat,npw,cfft,n1,n2,n3,n4,n5,n6,kg_k,istwf_k,iflag,me_g0,shiftg,symm,xnorm)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iflag,istwf_k,n1,n2,n3,n4,n5,n6,ndat,npw,me_g0
 real(dp),intent(in) :: xnorm
!arrays
 integer,intent(in) :: kg_k(3,npw),shiftg(3),symm(3,3)
 real(dp),intent(inout) :: cfft(2,n4,n5,n6*ndat),cg(2,npw*ndat)

!Local variables-------------------------------
!scalars
 integer :: i1,i1inv,i2,i2inv,i3,i3inv,id1,id2,id3,idat,ipw
 integer :: j1,j2,j3,l1,l2,l3,npwmin,use_symmetry,i3dat,i3invdat,i2invdat,ipwdat,i2dat
 character(len=500) :: msg
!arrays
 integer :: identity(3,3)
 integer :: i1inver(n1),i2inver(n2),i3inver(n3)

! *************************************************************************

 DBG_ENTER("COLL")

 !In the case of special k-points, invariant under time-reversal,
 !but not Gamma, initialize the inverse coordinates. !Remember indeed that
 !
 !  u_k(G) = u_{k+G0}(G-G0); u_{-k}(G) = u_k(G)^* and therefore:
 !  u_{G0/2}(G) = u_{G0/2}(-G-G0)^*.

 if (istwf_k>=2) then
   if(istwf_k==2 .or. istwf_k==4 .or. istwf_k==6 .or. istwf_k==8)then
     i1inver(1)=1
     do i1=2,n1
       i1inver(i1)=n1+2-i1
     end do
   else
     do i1=1,n1
       i1inver(i1)=n1+1-i1
     end do
   end if
   if(istwf_k>=2 .and. istwf_k<=5)then
     i2inver(1)=1
     do i2=2,n2
       i2inver(i2)=n2+2-i2
     end do
   else
     do i2=1,n2
       i2inver(i2)=n2+1-i2
     end do
   end if
   if(istwf_k==2 .or. istwf_k==3 .or. istwf_k==6 .or. istwf_k==7)then
     i3inver(1)=1
     do i3=2,n3
       i3inver(i3)=n3+2-i3
     end do
   else
     do i3=1,n3
       i3inver(i3)=n3+1-i3
     end do
   end if
 end if

 if (iflag==1 .or. iflag==2) then
   ! Insert cg into cfft with extra 0 s around outside:
   cfft = zero

   ! Take care of each plane wave, and complete cfft if needed
   if (istwf_k==1) then

     if (iflag==1) then
!$OMP PARALLEL DO PRIVATE(i1,i2,i3) IF (ndat>1)
       do idat=1,ndat
         do ipw=1,npw
           i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
           i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
           i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

           cfft(1,i1,i2,i3+n6*(idat-1))=cg(1,ipw+npw*(idat-1))
           cfft(2,i1,i2,i3+n6*(idat-1))=cg(2,ipw+npw*(idat-1))
         end do
       end do
     end if

     if (iflag==2) then
!$OMP PARALLEL DO PRIVATE(i1,i2,i3) IF (ndat>1)
       do idat=1,ndat
         do ipw=1,npw
           i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
           i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
           i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

           cfft(1,i1,i3,i2+n6*(idat-1))=cg(1,ipw+npw*(idat-1))
           cfft(2,i1,i3,i2+n6*(idat-1))=cg(2,ipw+npw*(idat-1))
         end do
       end do
     end if

   else if (istwf_k>=2) then

     npwmin=1
     if (istwf_k==2 .and. me_g0==1) then
       ! If gamma point, then cfft must be completed
       do idat=1,ndat
         cfft(1,1,1,1+n6*(idat-1))=cg(1,1+npw*(idat-1))
         cfft(2,1,1,1+n6*(idat-1))=zero
       end do
       npwmin=2
     end if

     if (iflag==1) then
!$OMP PARALLEL DO PRIVATE(i1,i1inv,i2,i2inv,i3,i3inv) IF (ndat>1)
       do idat=1,ndat
         do ipw=npwmin,npw
           i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
           i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
           i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1
           ! Construct the coordinates of -k-G
           i1inv=i1inver(i1) ; i2inv=i2inver(i2) ; i3inv=i3inver(i3)

           cfft(1,i1,i2,i3+n6*(idat-1))=cg(1,ipw+npw*(idat-1))
           cfft(2,i1,i2,i3+n6*(idat-1))=cg(2,ipw+npw*(idat-1))
           cfft(1,i1inv,i2inv,i3inv+n6*(idat-1))= cg(1,ipw+npw*(idat-1))
           cfft(2,i1inv,i2inv,i3inv+n6*(idat-1))=-cg(2,ipw+npw*(idat-1))
         end do
       end do
     end if

     if (iflag==2) then
!$OMP PARALLEL DO PRIVATE(i1,i1inv,i2,i2inv,i3,i3inv) IF (ndat>1)
       do idat=1,ndat
         do ipw=npwmin,npw
           i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
           i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
           i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

           ! Construct the coordinates of -k-G
           i1inv=i1inver(i1) ; i2inv=i2inver(i2) ; i3inv=i3inver(i3)

           cfft(1,i1,i3,i2+n6*(idat-1))=cg(1,ipw+npw*(idat-1))
           cfft(2,i1,i3,i2+n6*(idat-1))=cg(2,ipw+npw*(idat-1))
           cfft(1,i1inv,i3inv,i2inv+n6*(idat-1))= cg(1,ipw+npw*(idat-1))
           cfft(2,i1inv,i3inv,i2inv+n6*(idat-1))=-cg(2,ipw+npw*(idat-1))
         end do
       end do
     end if

   end if

 else if (iflag==-1 .or. iflag==-2) then

   use_symmetry=0
   identity(:,:)=0; identity(1,1)=1 ; identity(2,2)=1 ; identity(3,3)=1
   if(sum((symm(:,:)-identity(:,:))**2)/=0)use_symmetry=1
   if(sum(shiftg(:)**2)/=0)use_symmetry=1

   ! Extract cg from cfft, ignoring components outside range of cg:
   if (istwf_k==1) then

     if (use_symmetry==0) then
       if (iflag==-1) then
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,ipwdat,i3dat) IF (ndat>1)
         do idat=1,ndat
           do ipw=1,npw
             i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
             i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
             i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1
             ipwdat = ipw + (idat-1) * npw
             i3dat = i3 + (idat-1) * n6

             cg(1,ipwdat)=cfft(1,i1,i2,i3dat)*xnorm
             cg(2,ipwdat)=cfft(2,i1,i2,i3dat)*xnorm
           end do
         end do
       else
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,ipwdat,i2dat) IF (ndat>1)
         do idat=1,ndat
           do ipw=1,npw
             i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
             i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
             i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

             ipwdat = ipw + (idat-1) * npw
             i2dat = i2 + (idat-1) * n6

             cg(1,ipwdat)=cfft(1,i1,i3,i2dat)*xnorm
             cg(2,ipwdat)=cfft(2,i1,i3,i2dat)*xnorm
           end do
         end do
       end if
     else
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,j1,j2,j3,l1,l2,l3,ipwdat,i3dat) IF (ndat>1)
       do idat=1,ndat
         do ipw=1,npw
           l1=kg_k(1,ipw)+shiftg(1)
           l2=kg_k(2,ipw)+shiftg(2)
           l3=kg_k(3,ipw)+shiftg(3)
           j1=symm(1,1)*l1+symm(1,2)*l2+symm(1,3)*l3
           j2=symm(2,1)*l1+symm(2,2)*l2+symm(2,3)*l3
           j3=symm(3,1)*l1+symm(3,2)*l2+symm(3,3)*l3
           if(j1<0)j1=j1+n1; i1=j1+1
           if(j2<0)j2=j2+n2; i2=j2+1
           if(j3<0)j3=j3+n3; i3=j3+1

           ipwdat = ipw + (idat-1) * npw
           i3dat = i3 + (idat-1)*n6

           cg(1,ipwdat)=cfft(1,i1,i2,i3dat)*xnorm
           cg(2,ipwdat)=cfft(2,i1,i2,i3dat)*xnorm
         end do
       end do
     end if

   else if (istwf_k>=2) then

     npwmin=1
     if (istwf_k==2 .and. me_g0==1) then
       ! Extract cg from cfft, in a way that projects on a
       ! wavefunction with time-reversal symmetry
       do idat=1,ndat
         ipwdat = 1 + (idat-1) * npw
         i3dat = 1 + (idat-1)*n6
         cg(1,ipwdat)=cfft(1,1,1,i3dat)*xnorm
         cg(2,ipwdat)=zero
       end do
       npwmin=2
     end if

     if (use_symmetry==0) then

       if (iflag==-1) then
!$OMP PARALLEL DO PRIVATE(i1,i1inv,i2,i2inv,i3,i3inv,ipwdat,i3dat,i3invdat) IF (ndat>1)
         do idat=1,ndat
           do ipw=npwmin,npw
             i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
             i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
             i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

             ! Construct the coordinates of -k-G
             i1inv=i1inver(i1); i2inv=i2inver(i2); i3inv=i3inver(i3)

             ipwdat = ipw + (idat-1) * npw
             i3dat = i3 + (idat-1) * n6
             i3invdat = i3inv + (idat-1) * n6

             ! Here the time-reversal symmetry is used to project from cfft
             cg(1,ipwdat)=(cfft(1,i1,i2,i3dat) + cfft(1,i1inv,i2inv,i3invdat))*0.5d0*xnorm
             cg(2,ipwdat)=(cfft(2,i1,i2,i3dat) - cfft(2,i1inv,i2inv,i3invdat))*0.5d0*xnorm
           end do
         end do

       else
!$OMP PARALLEL DO PRIVATE(i1,i1inv,i2,i2inv,i3,i3inv,ipwdat,i2dat,i2invdat) IF (ndat>1)
         do idat=1,ndat
           do ipw=npwmin,npw
             i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
             i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
             i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

             ! Construct the coordinates of -k-G
             i1inv=i1inver(i1) ; i2inv=i2inver(i2) ; i3inv=i3inver(i3)

             ipwdat = ipw + (idat-1) * npw
             i2dat = i2 + (idat-1) * n6
             i2invdat = i2inv + (idat-1) * n6

             ! Here the time-reversal symmetry is used to project from cfft
             cg(1,ipwdat)=(cfft(1,i1,i3,i2dat) + cfft(1,i1inv,i3inv,i2invdat))*0.5d0*xnorm
             cg(2,ipwdat)=(cfft(2,i1,i3,i2dat) - cfft(2,i1inv,i3inv,i2invdat))*0.5d0*xnorm
           end do
         end do
       end if

     else ! Use symmetry
       id1=n1/2+2
       id2=n2/2+2
       id3=n3/2+2

!$OMP PARALLEL DO PRIVATE(i1,i1inv,i2,i2inv,i3,i3inv,j1,j2,j3,l1,l2,l3,ipwdat,i3dat,i3invdat) IF (ndat>1)
       do idat=1,ndat
         do ipw=npwmin,npw

           i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
           i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
           i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

           i1inv=i1inver(i1) ; i2inv=i2inver(i2) ; i3inv=i3inver(i3)

           l1=kg_k(1,ipw)+shiftg(1)
           l2=kg_k(2,ipw)+shiftg(2)
           l3=kg_k(3,ipw)+shiftg(3)
           j1=symm(1,1)*l1+symm(1,2)*l2+symm(1,3)*l3
           j2=symm(2,1)*l1+symm(2,2)*l2+symm(2,3)*l3
           j3=symm(3,1)*l1+symm(3,2)*l2+symm(3,3)*l3
           if(j1<0)j1=j1+n1 ; i1=j1+1
           if(j2<0)j2=j2+n2 ; i2=j2+1
           if(j3<0)j3=j3+n3 ; i3=j3+1

           ! Construct the coordinates of -k-G
           l1=i1inv-(i1inv/id1)*n1-1+shiftg(1)
           l2=i2inv-(i2inv/id2)*n2-1+shiftg(2)
           l3=i3inv-(i3inv/id3)*n3-1+shiftg(3)
           j1=symm(1,1)*l1+symm(1,2)*l2+symm(1,3)*l3
           j2=symm(2,1)*l1+symm(2,2)*l2+symm(2,3)*l3
           j3=symm(3,1)*l1+symm(3,2)*l2+symm(3,3)*l3
           if(j1<0)j1=j1+n1 ; i1inv=j1+1
           if(j2<0)j2=j2+n2 ; i2inv=j2+1
           if(j3<0)j3=j3+n3 ; i3inv=j3+1

           ipwdat = ipw + (idat-1) * npw
           i3dat = i3 + (idat-1) * n6
           i3invdat = i3inv + (idat-1) * n6

           ! Here the time-reversal symmetry is used to project from cfft
           cg(1,ipwdat)=(cfft(1,i1,i2,i3dat) + cfft(1,i1inv,i2inv,i3invdat))*0.5d0*xnorm
           cg(2,ipwdat)=(cfft(2,i1,i2,i3dat) - cfft(2,i1inv,i2inv,i3invdat))*0.5d0*xnorm
         end do
       end do
     end if

   end if

 else
   write(msg,'(a,i0,a)')'  iflag=',iflag,' not acceptable.'
   MSG_BUG(msg)
 end if

 DBG_EXIT("COLL")

end subroutine sphere
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/sphere_fft
!! NAME
!! sphere_fft
!!
!! FUNCTION
!! Array cg is defined in sphere with npw points. Insert cg inside box
!! of n1*n2*n3 points to define array cfft for fft box.
!! corresponds to given element in cg.  rest of cfft is filled with 0 s.
!!
!! iflag=1==>insert cg into cfft.
!! iflag=2==>insert cg into cfft, where the second and third dimension
!! have been switched (needed for new 2002 SGoedecker FFT)
!! iflag=-1==> extract cg from cfft.
!! iflag=-2==> extract cg from cfft, where the second and third dimension
!! have been switched (needed for new 2002 SGoedecker FFT)
!!  (WARNING : iflag=-2 cannot use symmetry operations)
!!
!! There is also the possibility to apply a symmetry operation,
!! as well as to make a shift in reciprocal space, or to multiply
!! by a constant factor, in the case iflag=-1.
!! Multiplication by a constant factor is also possible in the case iflag=-2.
!!
!! INPUTS
!! cg(2,npw)= contains values for npw G vectors in basis sphere
!! ndat=number of FFT to do in //
!! npw=number of G vectors in basis at this k point
!! cfft(2,n4,n5,n6) = fft box
!! n1,n2,n3=physical dimension of the box (cfft)
!! n4,n5,n6=memory dimension of cfft
!! kg_k(3,npw)=integer coordinates of G vectors in basis sphere
!! mpi_enreg=information about MPI parallelization
!! tab_fftwf2_local(n2)=local i2 indices in fourwf
!! nd2proc TO BE DESCRIBED SB 090831
!! iflag=option parameter. Possible values: -1, -2, 1, 2 ; this is used only in debug option
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! iflag=1 and 2, insert cg(input) into cfft(output)
!! iflag=-1 and -2, extract cg(output) from cfft(input)
!!
!! NOTES
!! cg and cfft are assumed to be of type COMPLEX, although this routine treats
!! them as real of twice the length to avoid nonstandard complex*16.
!!
!! WARNING
!! NO CHECK is DONE over iflag.
!!
!! TODO
!! Order arguments
!!
!! PARENTS
!!      fourwf
!!
!! CHILDREN
!!      xmpi_sum,xmpi_sum_master
!!
!! SOURCE

subroutine sphere_fft(cg,ndat,npw,cfft,n1,n2,n3,n4,n5,kg_k,tab_fftwf2_local,nd2proc)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1,n2,n3,n4,n5,nd2proc,ndat,npw
 integer,intent(in) :: tab_fftwf2_local(n2)
!arrays
 integer,intent(in) :: kg_k(3,npw)
 real(dp),intent(in) :: cg(2,npw*ndat)
 real(dp),intent(out) :: cfft(2,n4,n5,nd2proc*ndat)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i2_local,i3,idat,ipw

! *************************************************************************

!Insert cg into cfft with extra 0 s around outside:
 cfft = zero

!$OMP PARALLEL DO PRIVATE(i1,i2,i2_local,i3)
 do ipw=1,npw
   i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
   i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
   i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1
   i2_local = tab_fftwf2_local(i2)
   do idat=1,ndat
     cfft(1,i1,i3,i2_local + nd2proc*(idat-1))=cg(1,ipw+npw*(idat-1))
     cfft(2,i1,i3,i2_local + nd2proc*(idat-1))=cg(2,ipw+npw*(idat-1))
   end do
 end do

end subroutine sphere_fft
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/sphere_fft1
!! NAME
!! sphere_fft1
!!
!! FUNCTION
!! Array cg is defined in sphere with npw points. Insert cg inside box
!! of n1*n2*n3 points to define array cfft for fft box.
!! corresponds to given element in cg.  rest of cfft is filled with 0 s.
!!
!! iflag=1==>insert cg into cfft.
!! iflag=2==>insert cg into cfft, where the second and third dimension
!! have been switched (needed for new 2002 SGoedecker FFT)
!! iflag=-1==> extract cg from cfft.
!! iflag=-2==> extract cg from cfft, where the second and third dimension
!! have been switched (needed for new 2002 SGoedecker FFT)
!!  (WARNING : iflag=-2 cannot use symmetry operations)
!!
!! There is also the possibility to apply a symmetry operation,
!! as well as to make a shift in reciprocal space, or to multiply
!! by a constant factor, in the case iflag=-1.
!! Multiplication by a constant factor is also possible in the case iflag=-2.
!!
!! INPUTS
!! cg(2,npw)= contains values for npw G vectors in basis sphere
!! ndat=number of FFT to do in //
!! npw=number of G vectors in basis at this k point
!! cfft(2,n4,n5,n6) = fft box
!! n1,n2,n3=physical dimension of the box (cfft)
!! n4,n5,n6=memory dimension of cfft
!! kg_k(3,npw)=integer coordinates of G vectors in basis sphere
!! nd2proc TO BE DESCRIBED SB 090831
!! iflag=option parameter. Possible values: -1, -2, 1, 2
!! tab_fftwf2_local(n2)=local i2 indices in fourwf
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! iflag=1 and 2, insert cg(input) into cfft(output)
!! iflag=-1 and -2, extract cg(output) from cfft(input)
!!
!! NOTES
!! cg and cfft are assumed to be of type COMPLEX, although this routine treats
!! them as real of twice the length to avoid nonstandard complex*16.
!!
!! WARNING
!! NO CHECK is DONE over iflag.
!!
!! TODO
!!   Order arguments
!!   sphere_fft1 is similar to sphere_fft, the only difference being that ndat > 1 is not supported.
!!   Why? Should merge the two APIs.
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      xmpi_sum,xmpi_sum_master
!!
!! SOURCE

subroutine sphere_fft1(cg,ndat,npw,cfft,n1,n2,n3,n4,n5,n6,kg_k,tab_fftwf2_local)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1,n2,n3,n4,n5,n6,ndat,npw
!arrays
 integer,intent(in) :: kg_k(3,npw)
 integer,intent(in) :: tab_fftwf2_local(n2)
 real(dp),intent(in) :: cg(2,npw*ndat)
 real(dp),intent(inout) :: cfft(2,n4,n5,n6*ndat)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i2_local,i3,idat,ipw

! *************************************************************************

!Insert cg into cfft with extra 0 s around outside:

 cfft = zero
!$OMP PARALLEL DO PRIVATE(i1,i2,i2_local,i3)
 do idat=1,ndat
   do ipw=1,npw
     i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
     i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
     i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1
     i2_local = tab_fftwf2_local(i2) + n6*(idat-1)
     cfft(1,i1,i3,i2_local)=cg(1,ipw+npw*(idat-1))
     cfft(2,i1,i3,i2_local)=cg(2,ipw+npw*(idat-1))
   end do
 end do

end subroutine sphere_fft1
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/change_istwfk
!! NAME
!! change_istwfk
!!
!! FUNCTION
!! This function allows one to change the time-reversal storage mode (istwfk)
!! of a *full* set of u(G). It does not support MPI-FFT!
!!
!! INPUTS
!! from_npw=number of G vectors in input from_cg
!! from_kg_k(3,npw)=integer coordinates of the G vectors of from_cg
!! from_istwfk=option parameter that describes the storage in from_cg
!! to_npw=number of G vectors in output to_cg
!! to_kg_k(3,npw)=integer coordinates of the G vectors in to_cg
!! to_istwfk=option parameter that describes the storage in to_cg
!! n1,n2,n3=physical dimension of the box (must be large enough to contain the sphere, no check is done)
!! ndat=number of wavefunctions
!! from_cg(2,from_npw*ndat)= Input u(g) values
!!
!! OUTPUTS
!! to_cg(2,to_npw*ndat)= Output u(g) defined on the list of vectors to_kg_k with time-reversal mode to_istwfk
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      xmpi_sum,xmpi_sum_master
!!
!! SOURCE

subroutine change_istwfk(from_npw,from_kg,from_istwfk,to_npw,to_kg,to_istwfk,n1,n2,n3,ndat,from_cg,to_cg)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: from_npw,from_istwfk,to_npw,to_istwfk,n1,n2,n3,ndat
!arrays
 integer,intent(in) :: from_kg(3,from_npw),to_kg(3,to_npw)
 real(dp),intent(inout) :: from_cg(2,from_npw*ndat) ! out due to sphere!
 real(dp),intent(inout) :: to_cg(2,to_npw*ndat)

!Local variables-------------------------------
!scalars
 integer :: n4,n5,n6
 real(dp),parameter :: xnorm1=one
!arrays
 integer,parameter :: shiftg0(3)=0,me_g0=1
 integer,parameter :: symmE(3,3)=reshape([1,0,0,0,1,0,0,0,1],[3,3])
 real(dp),allocatable :: cfft(:,:,:,:)

! *************************************************************************

 n4=2*(n1/2)+1
 n5=2*(n2/2)+1
 n6=2*(n3/2)+1

 ABI_MALLOC(cfft, (2,n4,n5,n6*ndat))

 ! iflag=1 ==> insert from_cg into cfft.
 call sphere(from_cg,ndat,from_npw,cfft,n1,n2,n3,n4,n5,n6,from_kg,from_istwfk,+1,me_g0,shiftg0,symmE,xnorm1)

 ! iflag=-1 ==> extract to_cg from cfft.
 call sphere(to_cg,ndat,to_npw,cfft,n1,n2,n3,n4,n5,n6,to_kg,to_istwfk,-1,me_g0,shiftg0,symmE,xnorm1)

 ABI_FREE(cfft)

end subroutine change_istwfk
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/switch
!! NAME
!!  switch
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

pure subroutine switch(n1dfft,n2,lot,n1,lzt,zt,zw)


!Arguments ------------------------------------
 integer,intent(in) :: n1dfft,n2,lot,n1,lzt
 real(dp),intent(in) :: zt(2,lzt,n1)
 real(dp),intent(inout) :: zw(2,lot,n2)

!Local variables-------------------------------
 integer :: i,j
! *************************************************************************

 do j=1,n1dfft
   do i=1,n2
     zw(1,j,i)=zt(1,i,j)
     zw(2,j,i)=zt(2,i,j)
   end do
 end do

end subroutine switch
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/switch_cent
!! NAME
!! switch_cent
!!
!! FUNCTION
!!   Perform the rotation:
!!
!!     input: I2,i1,j3,(jp3)
!!     output: i1,I2,j3,(jp3)
!!
!!   and padd the signal with zeros.
!!
!! INPUTS
!!  n1dfft=Number of 1D FFTs to perform
!!  max2=Max G_y in the small box enclosing the G-sphere.
!!  m2=Size of the small box enclosing the G-sphere along y
!!  n2=Dimension of the transform along y
!!  lot=Cache blocking factor.
!!  n1=Dimension of the transform along x
!!  lzt=Second dimension of z
!!  zt(2,lzt,n1)
!!
!! OUTPUT
!!  zw(2,lot,n2)=Cache working array
!!
!! PARENTS
!!
!! SOURCE

pure subroutine switch_cent(n1dfft,max2,m2,n2,lot,n1,lzt,zt,zw)


!Arguments ------------------------------------
 integer,intent(in) :: n1dfft,max2,m2,n2,lot,n1,lzt
 real(dp),intent(in) :: zt(2,lzt,n1)
 real(dp),intent(inout) :: zw(2,lot,n2)

!Local variables-------------------------------
 integer :: i,j

! *************************************************************************

 ! Here, zero and positive frequencies
 do j=1,n1dfft
   do i=1,max2+1
     zw(1,j,i)=zt(1,i,j)
     zw(2,j,i)=zt(2,i,j)
   end do
 end do

 ! Fill the center region with zeros
 do i=max2+2,n2-m2+max2+1
   do j=1,n1dfft
     zw(1,j,i)=zero
     zw(2,j,i)=zero
   end do
 end do

 ! Here, negative frequencies
 if (m2>=max2+2) then
   do j=1,n1dfft
     do i=max2+2,m2
       zw(1,j,i+n2-m2)=zt(1,i,j)
       zw(2,j,i+n2-m2)=zt(2,i,j)
     end do
   end do
 end if

end subroutine switch_cent
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/switchreal
!! NAME
!!  switchreal
!!
!! FUNCTION
!!   Perform the rotation:
!!
!!     input: I2,i1,j3,(jp3)
!!     output: i1,I2,j3,(jp3)
!!
!!   and padd the signal with zeros.
!!   Used for real wavefunctions.
!!
!! INPUTS
!!  includelast
!!  n1dfft=Number of 1D FFTs to perform
!!  n2=Dimension of the transform along y
!!  n2eff
!!  lot=Cache blocking factor.
!!  n1zt
!!  lzt
!!  zt(2,lzt,n1zt)
!!
!! OUTPUT
!!  zw(2,lot,n2)
!!
!! PARENTS
!!
!! SOURCE

pure subroutine switchreal(includelast,n1dfft,n2,n2eff,lot,n1zt,lzt,zt,zw)


!Arguments ------------------------------------
 integer,intent(in) :: includelast,n1dfft,n2,n2eff,lot,n1zt,lzt
 real(dp),intent(in) :: zt(2,lzt,n1zt)
 real(dp),intent(inout) :: zw(2,lot,n2)

!Local variables-------------------------------
 integer :: i,j
! *************************************************************************

 if (includelast==1) then

   ! Compute symmetric and antisymmetric combinations
   do j=1,n1dfft
     zw(1,j,1)=zt(1,1,2*j-1)
     zw(2,j,1)=zt(1,1,2*j  )
   end do
   do i=2,n2eff
     do j=1,n1dfft
       zw(1,j,i)=      zt(1,i,2*j-1)-zt(2,i,2*j)
       zw(2,j,i)=      zt(2,i,2*j-1)+zt(1,i,2*j)
       zw(1,j,n2+2-i)= zt(1,i,2*j-1)+zt(2,i,2*j)
       zw(2,j,n2+2-i)=-zt(2,i,2*j-1)+zt(1,i,2*j)
     end do
   end do

 else

   ! An odd number of FFTs
   ! Compute symmetric and antisymmetric combinations
   do j=1,n1dfft-1
     zw(1,j,1)=zt(1,1,2*j-1)
     zw(2,j,1)=zt(1,1,2*j  )
   end do
   zw(1,n1dfft,1)=zt(1,1,2*n1dfft-1)
   zw(2,n1dfft,1)=zero

   do i=2,n2eff
     do j=1,n1dfft-1
       zw(1,j,i)=      zt(1,i,2*j-1)-zt(2,i,2*j)
       zw(2,j,i)=      zt(2,i,2*j-1)+zt(1,i,2*j)
       zw(1,j,n2+2-i)= zt(1,i,2*j-1)+zt(2,i,2*j)
       zw(2,j,n2+2-i)=-zt(2,i,2*j-1)+zt(1,i,2*j)
     end do
     zw(1,n1dfft,i)=      zt(1,i,2*n1dfft-1)
     zw(2,n1dfft,i)=      zt(2,i,2*n1dfft-1)
     zw(1,n1dfft,n2+2-i)= zt(1,i,2*n1dfft-1)
     zw(2,n1dfft,n2+2-i)=-zt(2,i,2*n1dfft-1)
   end do
 end if

end subroutine switchreal
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/switchreal_cent
!! NAME
!!  switchreal_cent
!!
!! FUNCTION
!!   Perform the rotation:
!!
!!     input: I2,i1,j3,(jp3)
!!     output: i1,I2,j3,(jp3)
!!
!!   and padd the signal with zeros.
!!   Used for the Fourier transform of real wavefunctions.
!!
!! INPUTS
!!  includelast
!!  n1dfft=Number of 1D FFTs to perform
!!  max2
!!  n2=Dimension of the transform along y
!!  lot=Cache blocking factor.
!!  n1zt, lzt=Dimensions of zt.
!!  zt(2,lzt,n1zt)
!!
!! OUTPUT
!!  zw(2,lot,n2)
!!
!! PARENTS
!!
!! SOURCE

pure subroutine switchreal_cent(includelast,n1dfft,max2,n2,lot,n1zt,lzt,zt,zw)


!Arguments ------------------------------------
 integer,intent(in) :: includelast,n1dfft,max2,n2,lot,n1zt,lzt
 real(dp),intent(in) :: zt(2,lzt,n1zt)
 real(dp),intent(inout) :: zw(2,lot,n2)

!Local variables-------------------------------
 integer :: i,j
! *************************************************************************

 if (includelast==1) then

   ! Compute symmetric and antisymmetric combinations
   do j=1,n1dfft
    zw(1,j,1)=zt(1,1,2*j-1)
    zw(2,j,1)=zt(1,1,2*j  )
   end do

   do i=2,max2+1
    do j=1,n1dfft
     zw(1,j,i)=      zt(1,i,2*j-1)-zt(2,i,2*j)
     zw(2,j,i)=      zt(2,i,2*j-1)+zt(1,i,2*j)
     zw(1,j,n2+2-i)= zt(1,i,2*j-1)+zt(2,i,2*j)
     zw(2,j,n2+2-i)=-zt(2,i,2*j-1)+zt(1,i,2*j)
    end do
   end do

   if(max2+1<n2-max2)then
    do i=max2+2,n2-max2
     do j=1,n1dfft
      zw(1,j,i)=zero
      zw(2,j,i)=zero
     end do
    end do
   end if

 else
   ! Compute symmetric and antisymmetric combinations
   do j=1,n1dfft-1
     zw(1,j,1)=zt(1,1,2*j-1)
     zw(2,j,1)=zt(1,1,2*j  )
   end do

   zw(1,n1dfft,1)=zt(1,1,2*n1dfft-1)
   zw(2,n1dfft,1)=zero
   do i=2,max2+1
     do j=1,n1dfft-1
       zw(1,j,i)=      zt(1,i,2*j-1)-zt(2,i,2*j)
       zw(2,j,i)=      zt(2,i,2*j-1)+zt(1,i,2*j)
       zw(1,j,n2+2-i)= zt(1,i,2*j-1)+zt(2,i,2*j)
       zw(2,j,n2+2-i)=-zt(2,i,2*j-1)+zt(1,i,2*j)
     end do
     zw(1,n1dfft,i)=      zt(1,i,2*n1dfft-1)
     zw(2,n1dfft,i)=      zt(2,i,2*n1dfft-1)
     zw(1,n1dfft,n2+2-i)= zt(1,i,2*n1dfft-1)
     zw(2,n1dfft,n2+2-i)=-zt(2,i,2*n1dfft-1)
   end do

   if(max2+1<n2-max2)then
     do i=max2+2,n2-max2
       do j=1,n1dfft
        zw(1,j,i)=zero
        zw(2,j,i)=zero
       end do
     end do
   end if
 end if

end subroutine switchreal_cent
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/scramble
!! NAME
!!  scramble
!!
!! FUNCTION
!!  This routine performs the local rotation
!!
!!     input:  G1,R3,G2,(Gp2)
!!     output: G1,G2,R3,(Gp2)
!!
!! INPUTS
!!  i1=Index of x in the small box enclosing the G-sphere.
!!  j2
!!  lot=Cache blocking factor
!!  n1dfft=Number of 1D FFTs performed.
!!  md1,md2proc,nnd3=Used to dimension zmpi2
!!  n3=Dimension of the transform along z.
!!  zw(2,lot,n3): zw(:,1:n1dfft,n3) contains the lines transformed along z
!!
!! OUTPTU
!! zmpi2(2,md1,md2proc,nnd3)
!!
!! PARENTS
!!
!! SOURCE

pure subroutine scramble(i1,j2,lot,n1dfft,md1,n3,md2proc,nnd3,zw,zmpi2)


!Arguments ------------------------------------
 integer,intent(in) :: i1,j2,lot,n1dfft,md1,n3,md2proc,nnd3
 real(dp),intent(in) :: zw(2,lot,n3)
 real(dp),intent(inout) :: zmpi2(2,md1,md2proc,nnd3)

!Local variables-------------------------------
!scalars
 integer :: i3,i

! *************************************************************************

 do i3=1,n3
   do i=0,n1dfft-1
     zmpi2(1,i1+i,j2,i3)=zw(1,i+1,i3)
     zmpi2(2,i1+i,j2,i3)=zw(2,i+1,i3)
   end do
 end do

end subroutine scramble
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/fill
!! NAME
!!  fill
!!
!! FUNCTION
!!   Receives a set of z-lines in reciprocal space,
!!   insert the values in the cache work array zw (no padding)
!!
!! INPUTS
!!  nd1,nd3=Dimensions of the input array zf.
!!  lot=Cache blocking factor.
!!  n1dfft=Number of 1D FFTs to perform
!!  n3=Dimension of the transform along z
!!  zf(2,nd1,nd3)=Input array
!!
!! OUTPUT
!!  zw(2,lot,n3)=Cache work array with the z-lines.
!!
!! PARENTS
!!
!! SOURCE

pure subroutine fill(nd1,nd3,lot,n1dfft,n3,zf,zw)


!Arguments ------------------------------------
 integer,intent(in) :: nd1,nd3,lot,n1dfft,n3
 real(dp),intent(in) :: zf(2,nd1,nd3)
 real(dp),intent(inout) :: zw(2,lot,n3)

! local variables
 integer :: i1,i3

! *************************************************************************

 do i3=1,n3
   do i1=1,n1dfft
     zw(1,i1,i3)=zf(1,i1,i3)
     zw(2,i1,i3)=zf(2,i1,i3)
   end do
 end do

end subroutine fill
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/fill_cent
!! NAME
!! fill_cent
!!
!! FUNCTION
!!   Receives a set of z-lines in reciprocal space,
!!   insert the values in cache work array defined on the FFT box
!!   and pads the central of the frequency region with zeros.
!!
!! INPUTS
!!  md1,md3=Leading dimension of zf along x and z
!!  lot=second dimension of zw (cache blocking factor)
!!  n1dfft=Number of 1d transforms to be performed along z.
!!  max3=Max G_z in the small box enclosing the G-sphere.
!!  m3=Number of points in the *small* box enclosing the G-sphere
!!  n3=Dimension of the FFT transform along z
!!  zf(2,md1,md3)=x-z planes in reciprocal space
!!
!! OUTPUT
!!   zw(2,lot,n3)= Filled cache work array.
!!     zw(:,1:n1dfft,n3) contains the lines to be transformed along.
!!
!! PARENTS
!!
!! SOURCE

pure subroutine fill_cent(md1,md3,lot,n1dfft,max3,m3,n3,zf,zw)


!Arguments ------------------------------------
 integer,intent(in) :: md1,md3,lot,n1dfft,max3,m3,n3
 real(dp),intent(in) :: zf(2,md1,md3)
 real(dp),intent(inout) :: zw(2,lot,n3)

!Local variables-------------------------------
!scalars
 integer :: i1,i3

! *************************************************************************

 ! Here, zero and positive frequencies
 do i3=1,max3+1
   do i1=1,n1dfft
     zw(1,i1,i3)=zf(1,i1,i3)
     zw(2,i1,i3)=zf(2,i1,i3)
   end do
 end do

 ! Fill the center region with zeros
 do i3=max3+2,n3-m3+max3+1
   do i1=1,n1dfft
     zw(1,i1,i3)=zero
     zw(2,i1,i3)=zero
   end do
 end do

 ! Here, negative frequencies
 do i3=max3+2,m3
   do i1=1,n1dfft
     zw(1,i1,i3+n3-m3)=zf(1,i1,i3)
     zw(2,i1,i3+n3-m3)=zf(2,i1,i3)
   end do
 end do

end subroutine fill_cent
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/unfill
!! NAME
!!  unfill
!!
!! FUNCTION
!!  Move data from the cache work array to zf
!!
!! INPUTS
!!  nd1,nd3=Dimensions of the input array zf.
!!  lot=Cache blocking factor.
!!  n1dfft=Number of 1D FFTs to perform
!!  n3=Dimension of the transform along z
!!  zw(2,lot,n3)=Cache work array with the z-lines.
!!
!! OUTPUT
!!  zf(2,nd1,nd3)= zf(:,1:n1dfft,:1:n3) is filled with the results stored in zw
!!
!! PARENTS
!!
!! SOURCE

pure subroutine unfill(nd1,nd3,lot,n1dfft,n3,zw,zf)


!Arguments ------------------------------------
 integer,intent(in) :: nd1,nd3,lot,n1dfft,n3
 real(dp),intent(in) :: zw(2,lot,n3)
 real(dp),intent(inout) :: zf(2,nd1,nd3)

!Local variables-------------------------------
 integer :: i1,i3
! *************************************************************************

 do i3=1,n3
   do i1=1,n1dfft
     zf(1,i1,i3)=zw(1,i1,i3)
     zf(2,i1,i3)=zw(2,i1,i3)
   end do
 end do

end subroutine unfill
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/unfill_cent
!! NAME
!!  unfill_cent
!!
!! FUNCTION
!!  Transfer data from the cache work array to zf.
!!  Takes into account zero padding (only the non-zero entries are moved)
!!
!! INPUTS
!!  md1,md3=Leading dimension of zf along x and z
!!  lot=Cache blocking factor.
!!  n1dfft=Number of 1d transforms performed along z.
!!  max3=Max index of G_z the small box enclosing the G-sphere.
!!  m3=Number of points in the *small* box enclosing the G-sphere
!!  n3=Dimension of the FFT transform along z
!!  zw(2,lot,n3)=Cache work array
!!
!! OUTPUT
!!  zf(2,md1,md3)= zf(:,1:n1dfft,1:m3) is filled with the non-zero components.
!!
!! PARENTS
!!
!! SOURCE

pure subroutine unfill_cent(md1,md3,lot,n1dfft,max3,m3,n3,zw,zf)


!Arguments ------------------------------------
 integer,intent(in) :: md1,md3,lot,n1dfft,max3,m3,n3
 real(dp),intent(in) :: zw(2,lot,n3)
 real(dp),intent(inout) :: zf(2,md1,md3)

!Local variables-------------------------------
 integer :: i1,i3
! *************************************************************************

 ! Here, zero and positive frequencies
 do i3=1,max3+1
   do i1=1,n1dfft
     zf(1,i1,i3)=zw(1,i1,i3)
     zf(2,i1,i3)=zw(2,i1,i3)
   end do
 end do

 ! Here, negative frequencies
 do i3=max3+2,m3
   do i1=1,n1dfft
     zf(1,i1,i3)=zw(1,i1,i3+n3-m3)
     zf(2,i1,i3)=zw(2,i1,i3+n3-m3)
   end do
 end do

end subroutine unfill_cent
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/unmpiswitch
!! NAME
!!  unmpiswitch
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_fftw3,m_sg2002
!!
!! CHILDREN
!!      xmpi_sum,xmpi_sum_master
!!
!! SOURCE

pure subroutine unmpiswitch(j3,n1dfft,Jp2st,J2st,lot,n1,nd2proc,nd3proc,nproc,ioption,zw,zmpi1)


!Arguments ------------------------------------
 integer,intent(in) :: j3,n1dfft,lot,n1,nd2proc,nd3proc,nproc,ioption
 integer,intent(inout) :: Jp2st,J2st
 real(dp),intent(in) :: zw(2,lot,n1)
 real(dp),intent(inout) :: zmpi1(2,n1,nd2proc,nd3proc,nproc)

!Local variables-------------------------------
 integer :: i1,jp2,j2,ind,jjp2,mfft,jj2

! *************************************************************************

 mfft=0
 if (ioption == 2) then
   do Jp2=Jp2st,nproc
     do J2=J2st,nd2proc
       mfft=mfft+1
       if (mfft.gt.n1dfft) then
         Jp2st=Jp2
         J2st=J2
         return
       end if
       do I1=1,n1
         zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1)
         zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1)
       end do
     end do
     J2st=1
   end do

 else
   do Jp2=Jp2st,nproc
     do J2=J2st,nd2proc
       mfft=mfft+1
       if (mfft.gt.n1dfft) then
         Jp2st=Jp2
         J2st=J2
         return
       end if
       ind=(Jp2-1) * nd2proc + J2
       jj2=(ind-1)/nproc +1

       !jjp2=modulo(ind,nproc) +1
       jjp2=modulo(ind-1,nproc)+1

       do I1=1,n1
         zmpi1(1,I1,jj2,j3,jjp2)=zw(1,mfft,I1)
         zmpi1(2,I1,jj2,j3,jjp2)=zw(2,mfft,I1)
       end do
     end do
     J2st=1
   end do
 end if

end subroutine unmpiswitch
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/unmpiswitch_cent
!! NAME
!!  unmpiswitch_cent
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_fftw3,m_sg2002
!!
!! CHILDREN
!!      xmpi_sum,xmpi_sum_master
!!
!! SOURCE

pure subroutine unmpiswitch_cent(j3,n1dfft,Jp2stf,J2stf,lot,max1,md1,m1,n1,md2proc,nd3proc,nproc,ioption,zw,zmpi1)


!Arguments ------------------------------------
 integer,intent(in) :: j3,n1dfft,lot,max1,md1,m1,n1,md2proc,nd3proc,nproc,ioption
 integer,intent(inout) :: Jp2stf,J2stf
 real(dp),intent(inout) :: zmpi1(2,md1,md2proc,nd3proc,nproc)
 real(dp),intent(in) :: zw(2,lot,n1)

!Local variables-------------------------------
 integer :: mfft,Jp2,J2,I1,ind,jj2,jjp2

! *************************************************************************

 mfft=0

 if (ioption == 2) then
   do Jp2=Jp2stf,nproc
     do J2=J2stf,md2proc
       mfft=mfft+1

       if (mfft.gt.n1dfft) then
         Jp2stf=Jp2
         J2stf=J2
         return
       end if

       ! Here, zero and positive frequencies
       do I1=1,max1+1
         zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1)
         zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1)
       end do

       ! Here, negative frequencies
       do I1=max1+2,m1
         zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1+n1-m1)
         zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1+n1-m1)
       end do

     end do
     J2stf=1
   end do

 else
   do Jp2=Jp2stf,nproc
     do J2=J2stf,md2proc
       mfft=mfft+1
       if (mfft.gt.n1dfft) then
         Jp2stf=Jp2
         J2stf=J2
         return
       end if
       ind=(Jp2-1) * md2proc + J2
       jj2=(ind-1)/nproc +1

       !jjp2=modulo(ind,nproc) +1
       jjp2=modulo(ind-1,nproc)+1

       ! Here, zero and positive frequencies
       do I1=1,max1+1
         zmpi1(1,I1,Jj2,j3,Jjp2)=zw(1,mfft,I1)
         zmpi1(2,I1,Jj2,j3,Jjp2)=zw(2,mfft,I1)
       end do

       ! Here, negative frequencies
       do I1=max1+2,m1
         zmpi1(1,I1,Jj2,j3,Jjp2)=zw(1,mfft,I1+n1-m1)
         zmpi1(2,I1,Jj2,j3,Jjp2)=zw(2,mfft,I1+n1-m1)
       end do
     end do
     J2stf=1
   end do
 end if

end subroutine unmpiswitch_cent
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/unscramble
!! NAME
!!  unscramble
!!
!! FUNCTION
!!
!! INPUTS
!!  i1
!!  j2
!!  lot=Cache blocking factor.
!!  n1dfft=Number of 1D FFTs to perform
!!  md1,n3,md2proc,nnd3
!!  zmpi2(2,md1,md2proc,nnd3)
!!
!! OUTPUT
!!  zw(2,lot,n3)= cache work array
!!
!! PARENTS
!!
!! SOURCE

pure subroutine unscramble(i1,j2,lot,n1dfft,md1,n3,md2proc,nnd3,zmpi2,zw)


!Arguments ------------------------------------
 integer,intent(in) :: i1,j2,lot,n1dfft,md1,n3,md2proc,nnd3
 real(dp),intent(in) :: zmpi2(2,md1,md2proc,nnd3)
 real(dp),intent(inout) :: zw(2,lot,n3)

!Local variables-------------------------------
!scalars
 integer :: i,i3

! *************************************************************************

 do i3=1,n3
   do i=0,n1dfft-1
     zw(1,i+1,i3)=zmpi2(1,i1+i,j2,i3)
     zw(2,i+1,i3)=zmpi2(2,i1+i,j2,i3)
   end do
 end do

end subroutine unscramble
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/unswitch
!! NAME
!!  unswitch
!!
!! FUNCTION
!!
!! INPUTS
!!  n1dfft=Number of 1D FFTs
!!  n2=Dimension of the transform along y
!!  lot=Cache blocking factor.
!!  n1=Dimension of the transform along x.
!!  lzt
!!  zw(2,lot,n2)=Cache work array
!!
!! OUTPUT
!!  zt(2,lzt,n1)
!!
!! PARENTS
!!
!! SOURCE

pure subroutine unswitch(n1dfft,n2,lot,n1,lzt,zw,zt)


!Arguments ------------------------------------
 integer,intent(in) :: n1dfft,n2,lot,n1,lzt
 real(dp),intent(in) :: zw(2,lot,n2)
 real(dp),intent(inout) :: zt(2,lzt,n1)

!Local variables-------------------------------
 integer :: i,j
! *************************************************************************

 do j=1,n1dfft
   do i=1,n2
     zt(1,i,j)=zw(1,j,i)
     zt(2,i,j)=zw(2,j,i)
   end do
 end do

end subroutine unswitch
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/unswitch_cent
!! NAME
!!  unswitch_cent
!!
!! FUNCTION
!!
!! INPUTS
!!  n1dfft=Number of 1D FFTs to perform
!!  max2=Max G_y in the small box enclosing the G-sphere.
!!  m2=Size of the small box enclosing the G-sphere along y
!!  n2=Dimension of the transform along y
!!  lot=Cache blocking factor.
!!  n1=Dimension of the transform along x
!!  lzt
!!  zw(2,lot,n2)=Cache working array
!!
!! OUTPUT
!!  zt(2,lzt,n1)
!!
!! PARENTS
!!
!! SOURCE

pure subroutine unswitch_cent(n1dfft,max2,m2,n2,lot,n1,lzt,zw,zt)


!Arguments ------------------------------------
 integer,intent(in) :: n1dfft,max2,m2,n2,lot,n1,lzt
 real(dp),intent(in) :: zw(2,lot,n2)
 real(dp),intent(inout) :: zt(2,lzt,n1)

!Local variables-------------------------------
 integer :: i,j
! *************************************************************************

! Here, zero and positive frequencies
 do j=1,n1dfft
   do i=1,max2+1
     zt(1,i,j)=zw(1,j,i)
     zt(2,i,j)=zw(2,j,i)
   end do
 end do

! Here, negative frequencies
 if(m2>=max2+2)then
   do j=1,n1dfft
     do i=max2+2,m2
       zt(1,i,j)=zw(1,j,i+n2-m2)
       zt(2,i,j)=zw(2,j,i+n2-m2)
     end do
   end do
 end if

end subroutine unswitch_cent
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/unswitchreal
!! NAME
!!  unswitchreal
!!
!! FUNCTION
!!
!! INPUTS
!!  n1dfft=Number of 1D FFTs to perform
!!  n2=Dimension of the transform along y
!!  n2eff=
!!  lot=Cache blocking factor.
!!  n1zt
!!  lzt
!!  zw(2,lot,n2)=Cache working array
!!
!! OUTPUT
!!  zt(2,lzt,n1)
!!
!! PARENTS
!!
!! SOURCE

pure subroutine unswitchreal(n1dfft,n2,n2eff,lot,n1zt,lzt,zw,zt)


!Arguments ------------------------------------
 integer,intent(in) :: n1dfft,n2,n2eff,lot,n1zt,lzt
 real(dp),intent(in) :: zw(2,lot,n2)
 real(dp),intent(inout) :: zt(2,lzt,n1zt)

!Local variables-------------------------------
 integer :: i,j
! *************************************************************************

! Decompose symmetric and antisymmetric parts
 do j=1,n1dfft
   zt(1,1,2*j-1)=zw(1,j,1)
   zt(2,1,2*j-1)=zero
   zt(1,1,2*j)  =zw(2,j,1)
   zt(2,1,2*j)  =zero
 end do

 do i=2,n2eff
   do j=1,n1dfft
     zt(1,i,2*j-1)= (zw(1,j,i)+zw(1,j,n2+2-i))*half
     zt(2,i,2*j-1)= (zw(2,j,i)-zw(2,j,n2+2-i))*half
     zt(1,i,2*j)  = (zw(2,j,i)+zw(2,j,n2+2-i))*half
     zt(2,i,2*j)  =-(zw(1,j,i)-zw(1,j,n2+2-i))*half
   end do
 end do

end subroutine unswitchreal
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/unswitchreal_cent
!! NAME
!! unswitchreal_cent
!!
!! FUNCTION
!!
!! INPUTS
!!  n1dfft=Number of 1D FFTs to perform
!!  max2=Max G_y in the small box enclosing the G-sphere.
!!  n2=Dimension of the transform along y
!!  lot=Cache blocking factor.
!!  n1zt
!!  lzt
!!  zw(2,lot,n2)=Cache working array
!!
!! OUTPUT
!!  zt(2,lzt,n1)
!!
!! PARENTS
!!
!! SOURCE

pure subroutine unswitchreal_cent(n1dfft,max2,n2,lot,n1zt,lzt,zw,zt)


!Arguments ------------------------------------
 integer,intent(in) :: n1dfft,max2,n2,lot,n1zt,lzt
 real(dp),intent(in) :: zw(2,lot,n2)
 real(dp),intent(inout) :: zt(2,lzt,n1zt)

!Local variables-------------------------------
 integer :: i,j
! *************************************************************************

 do j=1,n1dfft
   zt(1,1,2*j-1)=zw(1,j,1)
   zt(2,1,2*j-1)=zero
   zt(1,1,2*j)  =zw(2,j,1)
   zt(2,1,2*j)  =zero
 end do

 do i=2,max2+1
   do j=1,n1dfft
     zt(1,i,2*j-1)= (zw(1,j,i)+zw(1,j,n2+2-i))*half
     zt(2,i,2*j-1)= (zw(2,j,i)-zw(2,j,n2+2-i))*half
     zt(1,i,2*j)  = (zw(2,j,i)+zw(2,j,n2+2-i))*half
     zt(2,i,2*j)  =-(zw(1,j,i)-zw(1,j,n2+2-i))*half
   end do
 end do

!       Here, zero and positive frequencies
!        do 90,j=1,n1dfft
!        do 90,i=1,max2+1
!        zt(1,i,j)=zw(1,j,i)
!        zt(2,i,j)=zw(2,j,i)
!90      continue

!       Here, negative frequencies
!        if(m2>=max2+2)then
!         do 110,j=1,n1dfft
!         do 110,i=max2+2,m2
!         zt(1,i,j)=zw(1,j,i+n2-m2)
!         zt(2,i,j)=zw(2,j,i+n2-m2)
!110      continue
!        end if

end subroutine unswitchreal_cent
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/mpiswitch
!! NAME
!! mpiswitch
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_fftw3,m_sg2002
!!
!! CHILDREN
!!      xmpi_sum,xmpi_sum_master
!!
!! SOURCE

pure subroutine mpiswitch(j3,n1dfft,Jp2st,J2st,lot,n1,nd2proc,nd3proc,nproc,ioption,zmpi1,zw)


!Arguments ------------------------------------
 integer,intent(in) :: j3,n1dfft,lot,n1,nd2proc,nd3proc,nproc,ioption
 integer,intent(inout) :: Jp2st,J2st
 real(dp),intent(in) :: zmpi1(2,n1,nd2proc,nd3proc,nproc)
 real(dp),intent(inout) :: zw(2,lot,n1)

!Local variables-------------------------------
 integer :: Jp2,J2,I1,ind,jj2,mfft,jjp2

! *************************************************************************
 mfft=0

 if (ioption /= 1) then
   do Jp2=Jp2st,nproc
     do J2=J2st,nd2proc
       mfft=mfft+1
       if (mfft.gt.n1dfft) then
         Jp2st=Jp2
         J2st=J2
         return
       end if
       do I1=1,n1
         zw(1,mfft,I1)=zmpi1(1,I1,J2,j3,Jp2)
         zw(2,mfft,I1)=zmpi1(2,I1,J2,j3,Jp2)
       end do
     end do
     J2st=1
   end do

 else
   do Jp2=Jp2st,nproc
     do J2=J2st,nd2proc
       mfft=mfft+1
       if (mfft.gt.n1dfft) then
         Jp2st=Jp2
         J2st=J2
         return
       end if
       ind=(Jp2-1) * nd2proc + J2
       jj2=(ind-1)/nproc +1

       !jjp2=modulo(ind,nproc) +1
       jjp2=modulo(ind-1,nproc)+1

       !in other words: mfft=(jj2-1)*nproc+jjp2 (modulo case)
       !istead of mfft=(Jjp2-1) * nd2proc + Jj2 (slice case)
       !with 1<=jjp2<=nproc, jj2=1,nd2proc
       do I1=1,n1
         ! zw(1,mfft,I1)=zmpi1(1,I1,J2,j3,Jp2)
         ! zw(2,mfft,I1)=zmpi1(2,I1,J2,j3,Jp2)
         zw(1,mfft,I1)=zmpi1(1,I1,jj2,j3,jjp2)
         zw(2,mfft,I1)=zmpi1(2,I1,jj2,j3,jjp2)
       end do
     end do
     J2st=1
   end do
 end if

end subroutine mpiswitch
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/mpiswitch_cent
!! NAME
!! mpiswitch_cent
!!
!! FUNCTION
!!   Perform the local rotation
!!
!!     input: I1,J2,j3,Jp2,(jp3)
!!     output: J2,Jp2,I1,j3,(jp3)
!!
!!   and fill the central region of the frequency spectrum with zeros
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_fftw3,m_sg2002
!!
!! CHILDREN
!!      xmpi_sum,xmpi_sum_master
!!
!! SOURCE


pure subroutine mpiswitch_cent(j3,n1dfft,Jp2stb,J2stb,lot,max1,md1,m1,n1,md2proc,&
&  nd3proc,nproc,ioption,zmpi1,zw,max2,m2,n2)


!Arguments ------------------------------------
 integer,intent(in) :: j3,n1dfft,lot,max1,md1,m1,n1,md2proc,nd3proc,nproc,ioption
 integer,intent(in) :: m2,max2,n2
 integer,intent(inout) :: Jp2stb,J2stb
 real(dp),intent(in) :: zmpi1(2,md1,md2proc,nd3proc,nproc)
 real(dp),intent(inout) :: zw(2,lot,n1)

!Local variables-------------------------------
 integer :: mfft,jp2,j2,jjp2,jj2,i1,ind

! *************************************************************************

 ABI_UNUSED((/m2,max2,n2/))

 mfft=0

 if (ioption /= 1) then
   do Jp2=Jp2stb,nproc
     do J2=J2stb,md2proc

       mfft=mfft+1
       if (mfft.gt.n1dfft) then
         Jp2stb=Jp2
         J2stb=J2
         !MSG_WARNING("Returning from mpiswithc_cent")
         return
       end if

       ! Here, zero and positive frequencies
       ! In zmpi1, they are stored from 1 to max1+1
       do I1=1,max1+1
         zw(1,mfft,I1)=zmpi1(1,I1,J2,j3,Jp2)
         zw(2,mfft,I1)=zmpi1(2,I1,J2,j3,Jp2)
       end do

       ! Fill the center region with zeros
       do I1=max1+2,n1-m1+max1+1
         zw(1,mfft,I1)=zero
         zw(2,mfft,I1)=zero
       end do

       ! Here, negative frequencies
       ! In zmpi1, they are stored from 1 to m1half
       do I1=max1+2,m1
         zw(1,mfft,I1+n1-m1)=zmpi1(1,I1,J2,j3,Jp2)
         zw(2,mfft,I1+n1-m1)=zmpi1(2,I1,J2,j3,Jp2)
       end do
     end do
     J2stb=1
   end do

 else
   do Jp2=Jp2stb,nproc
     do J2=J2stb,md2proc

       mfft=mfft+1
       if (mfft.gt.n1dfft) then
         Jp2stb=Jp2
         J2stb=J2
         !MSG_WARNING("Returning from mpiswithc_cent")
         return
       end if

       ind=(Jp2-1) * md2proc + J2
       jj2=(ind-1)/nproc +1

       !jjp2=modulo(ind,nproc) +1
       jjp2=modulo(ind-1,nproc)+1

       ! I gather consecutive I2 indexes in mfft in the modulo case
       ! Here, zero and positive frequencies
       ! In zmpi1, they are stored from 1 to max1+1
       do I1=1,max1+1
         zw(1,mfft,I1)=zmpi1(1,I1,Jj2,j3,Jjp2)
         zw(2,mfft,I1)=zmpi1(2,I1,Jj2,j3,Jjp2)
       end do

       ! Fill the center region with zeros
       do I1=max1+2,n1-m1+max1+1
         zw(1,mfft,I1)=zero
         zw(2,mfft,I1)=zero
       end do

       ! Here, negative frequencies
       ! In zmpi1, they are stored from 1 to m1half
       do I1=max1+2,m1
         zw(1,mfft,I1+n1-m1)=zmpi1(1,I1,Jj2,j3,Jjp2)
         zw(2,mfft,I1+n1-m1)=zmpi1(2,I1,Jj2,j3,Jjp2)
       end do

     end do
     J2stb=1
   end do
 end if

end subroutine mpiswitch_cent
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/mpifft_fg2dbox
!! NAME
!!  mpifft_fg2dbox
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

pure subroutine mpifft_fg2dbox(nfft,ndat,fofg,n1,n2,n3,n4,nd2proc,n6,fftn2_distrib,ffti2_local,me_fft,workf)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,ndat,n1,n2,n3,n4,nd2proc,n6,me_fft
!arrays
 integer,intent(in) :: fftn2_distrib(n2),ffti2_local(n2)
 real(dp),intent(in) :: fofg(2,nfft*ndat)
 real(dp),intent(inout) :: workf(2,n4,n6,nd2proc*ndat)

!Local variables-------------------------------
 integer :: idat,i1,i2,i3,i2_local,i2_ldat,fgbase

! *************************************************************************

 do idat=1,ndat
   do i3=1,n3
     do i2=1,n2
       if (fftn2_distrib(i2) == me_fft) then
         i2_local = ffti2_local(i2)
         i2_ldat = i2_local + (idat-1) * nd2proc
         fgbase= n1*(i2_local-1 + nd2proc*(i3-1)) + (idat-1) * nfft
         do i1=1,n1
           workf(1,i1,i3,i2_ldat)=fofg(1,i1+fgbase)
           workf(2,i1,i3,i2_ldat)=fofg(2,i1+fgbase)
         end do
       end if
     end do
   end do
 end do

end subroutine mpifft_fg2dbox
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/mpifft_fg2dbox_dpc
!! NAME
!!  mpifft_fg2dbox_dpc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

pure subroutine mpifft_fg2dbox_dpc(nfft,ndat,fofg,n1,n2,n3,n4,nd2proc,n6,fftn2_distrib,ffti2_local,me_fft,workf)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,ndat,n1,n2,n3,n4,nd2proc,n6,me_fft
!arrays
 integer,intent(in) :: fftn2_distrib(n2),ffti2_local(n2)
 real(dp),intent(in) :: fofg(2,nfft*ndat)
 complex(dpc),intent(inout) :: workf(n4,n6,nd2proc*ndat)

!Local variables-------------------------------
 integer :: idat,i1,i2,i3,i2_local,i2_ldat,fgbase

! *************************************************************************

 do idat=1,ndat
   do i3=1,n3
     do i2=1,n2
       if (fftn2_distrib(i2) == me_fft) then
         i2_local = ffti2_local(i2)
         i2_ldat = i2_local + (idat-1) * nd2proc
         fgbase= n1*(i2_local-1 + nd2proc*(i3-1)) + (idat-1) * nfft
         do i1=1,n1
           workf(i1,i3,i2_ldat)=CMPLX(fofg(1,i1+fgbase), fofg(2,i1+fgbase), kind=dpc)
         end do
       end if
     end do
   end do
 end do

end subroutine mpifft_fg2dbox_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/mpifft_dbox2fg
!! NAME
!!  mpifft_dbox2fg
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

pure subroutine mpifft_dbox2fg(n1,n2,n3,n4,nd2proc,n6,ndat,fftn2_distrib,ffti2_local,me_fft,workf,nfft,fofg)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1,n2,n3,n4,nd2proc,n6,ndat,me_fft,nfft
!arrays
 integer,intent(in) :: fftn2_distrib(n2),ffti2_local(n2)
 real(dp),intent(in) :: workf(2,n4,n6,nd2proc*ndat)
 real(dp),intent(out) :: fofg(2,nfft*ndat)

!Local variables-------------------------------
 integer :: idat,i1,i2,i3,i2_local,i2_ldat,fgbase
 real(dp) :: xnorm

! *************************************************************************

 xnorm=one/dble(n1*n2*n3)

 ! Transfer fft output to the original fft box
 do idat=1,ndat
   do i2=1,n2
     if( fftn2_distrib(i2) == me_fft) then
       i2_local = ffti2_local(i2)
       i2_ldat = i2_local + (idat-1) * nd2proc
       do i3=1,n3
         fgbase = n1*(i2_local - 1 + nd2proc*(i3-1)) + (idat - 1) * nfft
         do i1=1,n1
           fofg(1,i1+fgbase)=workf(1,i1,i3,i2_ldat)*xnorm
           fofg(2,i1+fgbase)=workf(2,i1,i3,i2_ldat)*xnorm
         end do
       end do
     end if
   end do
 end do

end subroutine mpifft_dbox2fg
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/mpifft_dbox2fg_dpc
!! NAME
!!  mpifft_dbox2fg_dpc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

pure subroutine mpifft_dbox2fg_dpc(n1,n2,n3,n4,nd2proc,n6,ndat,fftn2_distrib,ffti2_local,me_fft,workf,nfft,fofg)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1,n2,n3,n4,nd2proc,n6,ndat,me_fft,nfft
!arrays
 integer,intent(in) :: fftn2_distrib(n2),ffti2_local(n2)
 complex(dpc),intent(in) :: workf(n4,n6,nd2proc*ndat)
 real(dp),intent(out) :: fofg(2,nfft*ndat)

!Local variables-------------------------------
 integer :: idat,i1,i2,i3,i2_local,i2_ldat,fgbase
 real(dp) :: xnorm

! *************************************************************************

 xnorm=one/dble(n1*n2*n3)

 ! Transfer fft output to the original fft box
 do idat=1,ndat
   do i2=1,n2
     if( fftn2_distrib(i2) == me_fft) then
       i2_local = ffti2_local(i2)
       i2_ldat = i2_local + (idat-1) * nd2proc
       do i3=1,n3
         fgbase = n1*(i2_local - 1 + nd2proc*(i3-1)) + (idat - 1) * nfft
         do i1=1,n1
           fofg(1,i1+fgbase)=REAL (workf(i1,i3,i2_ldat))*xnorm
           fofg(2,i1+fgbase)=AIMAG(workf(i1,i3,i2_ldat))*xnorm
         end do
       end do
     end if
   end do
 end do

end subroutine mpifft_dbox2fg_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/mpifft_dbox2fr
!! NAME
!!  mpifft_dbox2fr
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

pure subroutine mpifft_dbox2fr(n1,n2,n3,n4,n5,nd3proc,ndat,fftn3_distrib,ffti3_local,me_fft,workr,cplex,nfft,fofr)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1,n2,n3,n4,n5,nd3proc,ndat,me_fft,nfft,cplex
!!arrays
 integer,intent(in) :: fftn3_distrib(n3),ffti3_local(n3)
 real(dp),intent(in) :: workr(2,n4,n5,nd3proc*ndat)
 real(dp),intent(out) :: fofr(cplex*nfft*ndat)

!Local variables-------------------------------
 integer :: idat,i1,i2,i3,i3_local,i3_ldat,frbase

! *************************************************************************

 select case (cplex)
 case (1)

   do idat=1,ndat
     do i3=1,n3
       if( fftn3_distrib(i3) == me_fft) then
         i3_local = ffti3_local(i3)
         i3_ldat = i3_local + (idat - 1) * nd3proc
         do i2=1,n2
           frbase=n1*(i2-1+n2*(i3_local-1)) + (idat - 1) * nfft
           do i1=1,n1
             fofr(i1+frbase)=workr(1,i1,i2,i3_ldat)
           end do
         end do
       end if
     end do
   end do

 case (2)

   do idat=1,ndat
     do i3=1,n3
       if (fftn3_distrib(i3) == me_fft) then
         i3_local = ffti3_local(i3)
         i3_ldat = i3_local + (idat - 1) * nd3proc
         do i2=1,n2
           frbase=2*n1*(i2-1+n2*(i3_local-1)) + (idat - 1) * cplex * nfft
           !if (frbase > cplex*nfft*ndat - 2*n1) then
           !   write(std_out,*)i2,i3_local,frbase,cplex*nfft*ndat
           !   MSG_ERROR("frbase")
           !end if
           do i1=1,n1
             fofr(2*i1-1+frbase)=workr(1,i1,i2,i3_ldat)
             fofr(2*i1  +frbase)=workr(2,i1,i2,i3_ldat)
           end do
         end do
       end if
     end do
   end do

 case default
   !MSG_BUG("Wrong cplex")
   fofr = huge(one)
 end select

end subroutine mpifft_dbox2fr
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/mpifft_dbox2fr_dpc
!! NAME
!!  mpifft_dbox2fr_dpc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

pure subroutine mpifft_dbox2fr_dpc(n1,n2,n3,n4,n5,nd3proc,ndat,fftn3_distrib,ffti3_local,me_fft,workr,cplex,nfft,fofr)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1,n2,n3,n4,n5,nd3proc,ndat,me_fft,nfft,cplex
!!arrays
 integer,intent(in) :: fftn3_distrib(n3),ffti3_local(n3)
 complex(dpc),intent(in) :: workr(n4,n5,nd3proc*ndat)
 real(dp),intent(out) :: fofr(cplex*nfft*ndat)

!Local variables-------------------------------
 integer :: idat,i1,i2,i3,i3_local,i3_ldat,frbase

! *************************************************************************

 select case (cplex)
 case (1)

   do idat=1,ndat
     do i3=1,n3
       if( fftn3_distrib(i3) == me_fft) then
         i3_local = ffti3_local(i3)
         i3_ldat = i3_local + (idat - 1) * nd3proc
         do i2=1,n2
           frbase=n1*(i2-1+n2*(i3_local-1)) + (idat - 1) * nfft
           do i1=1,n1
             fofr(i1+frbase)=REAL(workr(i1,i2,i3_ldat))
           end do
         end do
       end if
     end do
   end do

 case (2)

   do idat=1,ndat
     do i3=1,n3
       if (fftn3_distrib(i3) == me_fft) then
         i3_local = ffti3_local(i3)
         i3_ldat = i3_local + (idat - 1) * nd3proc
         do i2=1,n2
           frbase=2*n1*(i2-1+n2*(i3_local-1)) + (idat - 1) * cplex * nfft
           do i1=1,n1
             fofr(2*i1-1+frbase)=REAL (workr(i1,i2,i3_ldat))
             fofr(2*i1  +frbase)=AIMAG(workr(i1,i2,i3_ldat))
           end do
         end do
       end if
     end do
   end do

 case default
   !MSG_BUG("Wrong cplex")
   fofr = huge(one)
 end select

end subroutine mpifft_dbox2fr_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/mpifft_fr2dbox
!! NAME
!!  mpifft_fr2dbox
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

pure subroutine mpifft_fr2dbox(cplex,nfft,ndat,fofr,n1,n2,n3,n4,n5,nd3proc,fftn3_distrib,ffti3_local,me_fft,workr)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,ndat,n1,n2,n3,n4,n5,nd3proc,me_fft
!!arrays
 integer,intent(in) :: fftn3_distrib(n3),ffti3_local(n3)
 real(dp),intent(in) :: fofr(cplex*nfft*ndat)
 real(dp),intent(inout) :: workr(2,n4,n5,nd3proc*ndat)

!Local variables-------------------------------
 integer :: idat,i1,i2,i3,i3_local,i3_ldat,frbase

! *************************************************************************

 select case (cplex)
 case (1)

   do idat=1,ndat
     do i3=1,n3
       if( me_fft == fftn3_distrib(i3) ) then
         i3_local = ffti3_local(i3)
         i3_ldat = i3_local + (idat-1) * nd3proc
         do i2=1,n2
           frbase=n1*(i2-1+n2*(i3_local-1)) + (idat-1) * nfft
           do i1=1,n1
             workr(1,i1,i2,i3_ldat)=fofr(i1+frbase)
             workr(2,i1,i2,i3_ldat)=zero
           end do
         end do
       end if
     end do
   end do

 case (2)

   do idat=1,ndat
     do i3=1,n3
       if( me_fft == fftn3_distrib(i3) ) then
         i3_local = ffti3_local(i3)
         i3_ldat = i3_local + (idat-1) * nd3proc
         do i2=1,n2
           frbase=2*n1*(i2-1+n2*(i3_local-1)) + (idat-1) * cplex * nfft
           do i1=1,n1
             workr(1,i1,i2,i3_ldat)=fofr(2*i1-1+frbase)
             workr(2,i1,i2,i3_ldat)=fofr(2*i1  +frbase)
           end do
         end do
       end if
     end do
   end do

 case default
   !MSG_BUG("Wrong cplex")
   workr = huge(one)
 end select

end subroutine mpifft_fr2dbox
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/mpifft_fr2dbox_dpc
!! NAME
!!  mpifft_fr2dbox_dpc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

pure subroutine mpifft_fr2dbox_dpc(cplex,nfft,ndat,fofr,n1,n2,n3,n4,n5,nd3proc,fftn3_distrib,ffti3_local,me_fft,workr)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,ndat,n1,n2,n3,n4,n5,nd3proc,me_fft
!!arrays
 integer,intent(in) :: fftn3_distrib(n3),ffti3_local(n3)
 real(dp),intent(in) :: fofr(cplex*nfft*ndat)
 complex(dpc),intent(inout) :: workr(n4,n5,nd3proc*ndat)

!Local variables-------------------------------
 integer :: idat,i1,i2,i3,i3_local,i3_ldat,frbase

! *************************************************************************

 select case (cplex)
 case (1)

   do idat=1,ndat
     do i3=1,n3
       if( me_fft == fftn3_distrib(i3) ) then
         i3_local = ffti3_local(i3)
         i3_ldat = i3_local + (idat-1) * nd3proc
         do i2=1,n2
           frbase=n1*(i2-1+n2*(i3_local-1)) + (idat-1) * nfft
           do i1=1,n1
             workr(i1,i2,i3_ldat)=CMPLX(fofr(i1+frbase), zero, kind=dpc)
           end do
         end do
       end if
     end do
   end do

 case (2)

   do idat=1,ndat
     do i3=1,n3
       if( me_fft == fftn3_distrib(i3) ) then
         i3_local = ffti3_local(i3)
         i3_ldat = i3_local + (idat-1) * nd3proc
         do i2=1,n2
           frbase=2*n1*(i2-1+n2*(i3_local-1)) + (idat-1) * cplex * nfft
           do i1=1,n1
             workr(i1,i2,i3_ldat)=CMPLX(fofr(2*i1-1+frbase), fofr(2*i1  +frbase), kind=dpc)
           end do
         end do
       end if
     end do
   end do

 case default
   !MSG_BUG("Wrong cplex")
   workr = huge(one)
 end select

end subroutine mpifft_fr2dbox_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/indfftrisc
!!
!! NAME
!! indfftrisc
!!
!! FUNCTION
!! Take the data for sphere boundary and list of planewave in sphere (kg_k), manipulate them
!! for convenient use in fourwf, and output them in indpw_k
!!
!! INPUTS
!! gbound(2*mgfft+4)=sphere boundary data
!! kg_k(3,npw_k)=reduced planewave coordinates
!! mgfft=maximum size of 1D FFTs
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! npw_k=number of G vectors in basis at this k point
!!
!! OUTPUT
!! indpw_k(4,npw_k)=array which gives fft box index for given basis sphere
!!   in a representation that is directly usable by sg_fftrisc.f
!! ngb=number of FFTs along z
!!
!! PARENTS
!!      m_sgfft
!!
!! CHILDREN
!!      xmpi_sum,xmpi_sum_master
!!
!! SOURCE

subroutine indfftrisc(gbound,indpw_k,kg_k,mgfft,ngb,ngfft,npw_k)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,npw_k
 integer,intent(out) :: ngb
!arrays
 integer,intent(in) :: gbound(2*mgfft+4),kg_k(3,npw_k),ngfft(18)
 integer,intent(out) :: indpw_k(4,npw_k)

!Local variables-------------------------------
!scalars
 integer :: g1,g2,i1,i2,i3,igb,index,ipw,n1,n2,n3
!arrays
 integer,allocatable :: index2d(:,:)

! *************************************************************************

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)

!First, generate a 2d index for each column of data
 ABI_ALLOCATE(index2d,(n1,n2))
 index2d(:,:)=0
 index=1
 igb=3
 do g2=0,gbound(2) ! g2max
   do g1=0,gbound(igb+1)  ! g1max
     index2d(g1+1,g2+1)=index
     index=index+1
   end do
   if(gbound(igb)<=-1)then ! g1min
     do g1=gbound(igb)+n1,n1-1
       index2d(g1+1,g2+1)=index
       index=index+1
     end do
   end if
   igb=igb+2
 end do

 if(gbound(1)<=-1)then ! g2min
   do g2=gbound(1)+n2,n2-1
     do g1=0,gbound(igb+1)
       index2d(g1+1,g2+1)=index
       index=index+1
     end do
     if(gbound(igb)<=-1)then
       do g1=gbound(igb)+n1,n1-1
         index2d(g1+1,g2+1)=index
         index=index+1
       end do
     end if
     igb=igb+2
   end do
 end if

 ngb=index-1


!The 2d index has been generated
!Now, contract indpw_k(1,ipw) and indpw_k(2,ipw) into indpw_k(4,ipw)
!indpw_k(1,ipw) and indpw_k(2,ipw) are used to hold inverse of index2d,
!and for them, the second index does not fill 1:npw . It is only
!the number of z-transform FFTs.

!$OMP PARALLEL DO PRIVATE(i1,i2,i3)
 do ipw=1,npw_k
   i1=kg_k(1,ipw); if(i1<0)i1=i1+n1 ; i1=i1+1
   i2=kg_k(2,ipw); if(i2<0)i2=i2+n2 ; i2=i2+1
   i3=kg_k(3,ipw); if(i3<0)i3=i3+n3 ; i3=i3+1
   indpw_k(4,ipw)=index2d(i1,i2)
   indpw_k(3,ipw)=i3
 end do

 do i1=1,n1
   do i2=1,n2
     index=index2d(i1,i2)
     if(index/=0)then
       indpw_k(1,index)=i1
       indpw_k(2,index)=i2
     end if
   end do
 end do

 ABI_DEALLOCATE(index2d)

end subroutine indfftrisc
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/kpgsph
!! NAME
!! kpgsph
!!
!! FUNCTION
!! Use reciprocal space metric gmet(3,3) to set up the list
!! of G vectors inside a sphere out to $ (1/2)*(2*\pi*(k+G))^2=ecut $.
!! If mkmem=0 and mpw=0, then only count the number of planewaves
!!
!! INPUTS
!!  ecut=planewave kinetic energy cutoff (hartrees)
!!  exchn2n3d=if 1, n2 and n3 are exchanged
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!  ikg=shift to be given to the location of the output data in the array kg
!!  ikpt=number of the k-point
!!  istwf_k=option parameter that describes the storage of wfs
!!  kpt(3)=reduced coords of k point (in terms of recip latt vecs)
!!  mkmem =maximum number of k points which can fit in core memory
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum number of planewaves as dimensioned in calling routine
!!
!! OUTPUT
!!  kg(3,mpw*mkmem)=dimensionless coords of resulting G vecs (integer)
!!  npw=resulting number of planewaves inside ecut centered at kpt
!!
!! SIDE EFFECTS
!!  mpi_enreg
!!    %me_g0=if 1, the plane wave G(0 0 0) is in the set of plane waves (and is the first)
!!    TODO: other SIDE EFFECTS on mpi_enreg should be described !!!
!!
!! NOTES
!!  Must take into account the time-reversal symmetry when istwf_k is not 1.
!!
!! PARENTS
!!      getmpw,initberry,initorbmag,kpgio,ks_ddiago,m_fft,m_fftcore,m_gsphere
!!      m_hamiltonian,m_wfd,mkpwind_k,wfconv
!!
!! CHILDREN
!!      xmpi_sum,xmpi_sum_master
!!
!! SOURCE


subroutine kpgsph(ecut,exchn2n3d,gmet,ikg,ikpt,istwf_k,kg,kpt,mkmem,mpi_enreg,mpw,npw)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: exchn2n3d,ikg,ikpt,istwf_k,mkmem,mpw
 integer,intent(out) :: npw
 real(dp),intent(in) :: ecut
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(inout) :: kg(3,mpw*mkmem)
 real(dp),intent(in) :: gmet(3,3),kpt(3)

!Local variables-------------------------------
!scalars
 integer :: i1,ig,ig1p,ig1pmax,ig2,ig2p,ig2pmax,ig2pmin,ig3,ig3p,ig3pmax
 integer :: ig3pmin,igtot,ii,ikpt_this_proc,in,ind,np_band,np_fft,npw_before,npw_remain,npw_split
 integer, save :: alloc_size=0
 real(dp) :: gap_pw,gmet11,gmet_trace,gmin,gs_fact,gs_part,gscut,v1,v2,v3,xx
 logical :: ipw_ok
 character(len=500) :: message
!arrays
 integer :: ngrid(3),nmax(3),nmin(3),n2,ierr
 integer,allocatable :: array_ipw(:),ig1arr(:),ig2arr(:)
 integer,allocatable :: ig3arr(:),kg_ind(:),kg_small(:,:)
 integer, allocatable :: npw_gather(:),npw_disp(:) ,kg_ind_gather(:),kg_small_gather(:,:)
 real(dp) :: kmax(3),minor(3),numer(3),tsec(2)
 real(dp),allocatable :: kg1arr(:),kg2arr(:),kg3arr(:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(23,1,tsec)
 if(istwf_k<1 .or. istwf_k>9)then
   write(message,'(3a,i0,a)' )&
&   'The variable istwf_k must be between 1 and 9, while',ch10,&
&   'the argument of the routine istwf_k =',istwf_k,'.'
   MSG_BUG(message)
 end if

 if(ikg+mpw>mkmem*mpw)then
   write(message,'(5a,i0,a,i0,a,i0,4a)')&
&   'The variables ikg, mkmem, and mpw  must satisfy ikg<=(mkmem-1)*mpw,',ch10,&
&   'while the arguments of the routine are',ch10,&
&   'ikg =',ikg,', mkmem =',mkmem,', mpw =',mpw,ch10,&
&   'Probable cause: Known error in invars1 for parallel spin-polarized case.',ch10,&
&   'Temporary solution: Change the number of parallel processes.'
   MSG_BUG(message)
 end if

 np_band=0
 if (mpw>0) then
   np_band=1; if(mpi_enreg%paral_kgb==1) np_band=max(1,mpi_enreg%nproc_band)
   alloc_size=max(alloc_size,(mpw+1)*np_band)
   ABI_ALLOCATE(kg_small,(3, alloc_size))
   ABI_ALLOCATE(kg_ind,(alloc_size))
   kg_ind(:)=0
 end if

!A larger array, that will be split on the correct processor
!G**2 cutoff, gscut=Ecut/2 /Pi^2

 gscut=0.5_dp*ecut*piinv**2

!In reduced coordinates, determine maximal value of k+G and G for each direction

 minor(1)=gmet(2,2)*gmet(3,3)-gmet(2,3)**2
 numer(1)=gmet(1,2)**2*gmet(3,3)-2.0_dp*gmet(1,2)*gmet(1,3)*gmet(2,3) +gmet(1,3)**2*gmet(2,2)
 minor(2)=gmet(1,1)*gmet(3,3)-gmet(1,3)**2
 numer(2)=gmet(2,3)**2*gmet(1,1)-2.0_dp*gmet(1,2)*gmet(1,3)*gmet(2,3) +gmet(2,1)**2*gmet(3,3)
 minor(3)=gmet(2,2)*gmet(1,1)-gmet(1,2)**2
 numer(3)=gmet(3,2)**2*gmet(1,1)-2.0_dp*gmet(1,2)*gmet(1,3)*gmet(2,3) +gmet(1,3)**2*gmet(2,2)

!Take the trace of the gmet tensor as dimensional reference
 gmet_trace=gmet(1,1)+gmet(2,2)+gmet(3,3)

 do ii=1,3
   xx=gmet(ii,ii)*minor(ii)-numer(ii)
   if(xx<tol10*gmet_trace**3 .or. minor(ii)<tol10*gmet_trace**2)then
     MSG_BUG('The metric tensor seem incorrect')
   end if
   kmax(ii)=sqrt(gscut*minor(ii)/xx)
   nmax(ii)=floor(kmax(ii)-kpt(ii)+tol10)
   nmin(ii)=ceiling(-kmax(ii)-kpt(ii)-tol10)
   ngrid(ii)=nmax(ii)-nmin(ii)+1
 end do
!perform looping over fft box grid of size ngfft(1)*ngfft(2)*ngfft(3):
 ig=0;ind=0
 in=0
 gmet11=gmet(1,1)

!Set up standard search sequence for grid points, in standard storage mode :
!0 1 2 3 ... nmax nmin ... -1
!If the mode is not standard, then some part of the FFT grid must be selected
!
 ABI_ALLOCATE(ig1arr,(ngrid(1)))
 ABI_ALLOCATE(ig2arr,(ngrid(2)))
 ABI_ALLOCATE(ig3arr,(ngrid(3)))
 ABI_ALLOCATE(kg1arr,(ngrid(1)))
 ABI_ALLOCATE(kg2arr,(ngrid(2)))
 ABI_ALLOCATE(kg3arr,(ngrid(3)))

 do ig1p=1,ngrid(1)
   ig1arr(ig1p)=ig1p-1
   if (ig1p-1>nmax(1)) ig1arr(ig1p)=ig1p-ngrid(1)-1
   kg1arr(ig1p)=kpt(1)+dble(ig1arr(ig1p))
 end do

!For the second direction, the number of points might depend on istwf_k
!---------------------------------------------------------------------
 ig2pmax=ngrid(2)
 if(istwf_k>=2 .and. exchn2n3d==0) ig2pmax=nmax(2)+1
 ABI_ALLOCATE(array_ipw,(-ig2pmax:ig2pmax))
 array_ipw(:)=0
 do ig2p=1,ig2pmax
   ig2arr(ig2p)=ig2p-1
   if (ig2p-1>nmax(2)) ig2arr(ig2p)=ig2p-ngrid(2)-1
   kg2arr(ig2p)=kpt(2)+dble(ig2arr(ig2p))
 end do

!For the third direction, the number of points might depend on istwf_k
!---------------------------------------------------------------------
 ig3pmax=ngrid(3)
 if (istwf_k>=2 .and. exchn2n3d==1) ig3pmax=nmax(3)+1

 do ig3p=1,ig3pmax
   ig3arr(ig3p)=ig3p-1
   if(ig3p-1>nmax(3)) ig3arr(ig3p)=ig3p-ngrid(3)-1
   kg3arr(ig3p)=kpt(3)+dble(ig3arr(ig3p))
 end do

!Performs loop on all grid points.
!---------------------------------------------------------------------
 igtot = 0
 if(exchn2n3d==0)then
   mpi_enreg%me_g0=0
   do ig3p=1,ngrid(3)
     ig3=ig3arr(ig3p)
     v3=kg3arr(ig3p)
     ig2pmin=1
     if( istwf_k>=2 .and. istwf_k<=5 .and. ig3<0)then
       ig2pmin=2
     end if
!    ig2pmax was initialized previously
     do ig2p=ig2pmin,ig2pmax
       ig2=ig2arr(ig2p)
!      PAY ATTENTION : the proc 0 must have me_g0=1
       ipw_ok = .true.
       if(mpi_enreg%paral_kgb==1 ) then
          n2 =mpi_enreg%distribfft%n2_coarse
          ipw_ok = ipw_ok.and.(mpi_enreg%me_fft == mpi_enreg%distribfft%tab_fftwf2_distrib( modulo(ig2,n2) + 1))
       end if
       if (ig2==0 .and. ipw_ok) mpi_enreg%me_g0=1
       v2=kg2arr(ig2p)
       gs_part=gmet(2,2)*v2*v2+gmet(3,3)*v3*v3+2.0_dp*gmet(2,3)*v2*v3
       gs_fact=2.0_dp*(gmet(1,2)*v2+gmet(3,1)*v3)
       ig1pmax=ngrid(1)
       if( (istwf_k==2.or.istwf_k==3) .and. ig3p==1 .and. ig2p==1)ig1pmax=nmax(1)+1
       do ig1p=1,ig1pmax
         v1=kg1arr(ig1p)
         gmin=gs_part+v1*(gs_fact+v1*gmet11)
!        If inside sphere:
         if (gmin<=gscut) then
           if (ipw_ok) then
             ig=ig+1  ! inside sphere
             igtot=igtot+1
             if (mpw>0.and.ig<=alloc_size) then
!              Keep coords of pw:
               kg_small(1,ig)=ig1arr(ig1p)
               kg_small(2,ig)=ig2
               kg_small(3,ig)=ig3
               kg_ind(ig)=igtot
             end if
             array_ipw(ig2)=array_ipw(ig2)+1
           else
             igtot=igtot+1
           end if
         end if
       end do !ig1p
     end do !ig2p
   end do !ig3p

 else ! if (exchn2n3d/=0)

!  ig2pmax was initialized previously
   mpi_enreg%me_g0=0
   do ig2p=1,ngrid(2)
     ig2=ig2arr(ig2p)
!    PAY ATTENTION : the proc 0 must have me_g0=1
     ipw_ok = .true.
     if(mpi_enreg%paral_kgb==1 ) then
        n2 =mpi_enreg%distribfft%n2_coarse
        ipw_ok = ipw_ok.and.(mpi_enreg%me_fft == mpi_enreg%distribfft%tab_fftwf2_distrib( modulo(ig2,n2) + 1))
     end if
     if(ig2==0 .and. istwf_k>=2 .and. ipw_ok) mpi_enreg%me_g0=1
     v2     =kg2arr(ig2p)
     ig3pmin=1
     if( (istwf_k==2 .or. istwf_k==3 .or. istwf_k==6 .or. istwf_k==7) .and. ig2<0)then
       ig3pmin=2
     end if
     do ig3p=ig3pmin,ig3pmax
       ig3=ig3arr(ig3p)
       v3=kg3arr(ig3p)
       gs_part=gmet(2,2)*v2*v2+gmet(3,3)*v3*v3+2.0_dp*gmet(2,3)*v2*v3
       gs_fact=2.0_dp*(gmet(1,2)*v2+gmet(3,1)*v3)
       ig1pmax=ngrid(1)
       if( (istwf_k==2.or.istwf_k==3) .and. ig3p==1 .and. ig2p==1)ig1pmax=nmax(1)+1
       do ig1p=1,ig1pmax
         v1=kg1arr(ig1p)
         gmin=gs_part+v1*(gs_fact+v1*gmet11)
!        If inside sphere:
         if (gmin<=gscut) then
           if (ipw_ok) then
             ig=ig+1  ! inside sphere
             igtot=igtot+1
!            Make sure not to overrun array, or simply do not store if mpw=0
             if (mpw>0.and.ig<=alloc_size) then
!              Keep coords of pw:
               kg_small(1,ig)=ig1arr(ig1p)
               kg_small(2,ig)=ig2
               kg_small(3,ig)=ig3
               kg_ind(ig)=igtot
             end if
           else
             igtot=igtot+1
           end if
         end if

       end do ! ig1p
     end do ! ig3p
!    end if ! if the ig2 plane is to be treated by this processor
   end do ! ig2p

 end if ! exchn2n3d==0 or ==1

!Total number of G vectors at this k point is assigned: npw
!when getcell = 1 it can be that ig exceeds mpw, the bound on kp_small
!here is a workaround:
 if (mpw*np_band > 0 .and. ig > alloc_size) then
   npw = mpw*np_band
 else
   npw=ig
 end if
 alloc_size = max(alloc_size,npw)

 ABI_DEALLOCATE(ig1arr)
 ABI_DEALLOCATE(ig2arr)
 ABI_DEALLOCATE(ig3arr)
 ABI_DEALLOCATE(kg1arr)
 ABI_DEALLOCATE(kg2arr)
 ABI_DEALLOCATE(kg3arr)

!BandFFT: plane-wave load balancing
 if (mpi_enreg%paral_kgb==1.and.mpi_enreg%nproc_fft>1.and. &
&    mpi_enreg%pw_unbal_thresh>zero.and. &
&    istwf_k==1) then
!  Check for reequilibration
   np_fft=max(1,mpi_enreg%nproc_fft)
   ABI_ALLOCATE(npw_gather,(np_fft)) ! Count pw before balancing
   call xmpi_allgather(npw,npw_gather,mpi_enreg%comm_fft,ierr)
   gap_pw = 100._dp*(maxval(npw_gather(:))-minval(npw_gather))/(1.*sum(npw_gather(:))/np_fft)
   write(message,'(a,f5.2)' ) ' Relative gap for number of plane waves between process (%): ',gap_pw
   call wrtout(std_out,message,'COLL')
   if(gap_pw > mpi_enreg%pw_unbal_thresh) then ! Effective reequilibration
     write(message,'(a,f5.2,a,i4,a,f5.2,a)') &
&        'Plane-wave unbalancing (',gap_pw,'%) for kpt ',ikpt,' is higher than threshold (',&
&        mpi_enreg%pw_unbal_thresh,'%); a plane-wave balancing procedure is activated!'
     call wrtout(std_out,message,'COLL')
     !Get optimal number
     npw_split=sum(npw_gather(:))
     npw=npw_split/np_fft
     npw_remain=modulo(npw_split,np_fft)
     if(mpi_enreg%me_fft < npw_remain) npw=npw+1
     ig=npw
#ifdef DEBUG_MODE
     write(message,*) 'New npw_fft = ', npw
     call wrtout(std_out,message,'COLL')
#endif
     alloc_size = max(alloc_size,npw)
     if(mpw>0 ) then  !Step for requilibration between fft process
       ABI_ALLOCATE(npw_disp,(np_fft))
       npw_disp=0
       do i1=2,np_fft
         npw_disp(i1) = npw_disp(i1-1) + npw_gather(i1-1)
       end do
       ABI_ALLOCATE(kg_ind_gather,(npw_split))
       ABI_ALLOCATE(kg_small_gather,(3,npw_split))
       call xmpi_allgatherv(kg_ind, npw_gather(mpi_enreg%me_fft+1) , &
&            kg_ind_gather,npw_gather,npw_disp,mpi_enreg%comm_fft,ierr)
       call xmpi_allgatherv(kg_small,3*npw_gather(mpi_enreg%me_fft+1), &
&            kg_small_gather,3*npw_gather, 3*npw_disp,mpi_enreg%comm_fft,ierr)
       npw_before=mpi_enreg%me_fft*(npw_split/np_fft)+min(npw_remain,mpi_enreg%me_fft)
       kg_small(:,1:npw)=kg_small_gather(:,npw_before+1:npw_before+npw)
       kg_ind(  1:npw)=kg_ind_gather(npw_before+1:npw_before+npw)
#ifdef DEBUG_MODE
       call wrtout(std_out,"keeping values done",'COLL')
#endif
       ABI_DEALLOCATE(npw_disp)
       ABI_DEALLOCATE(kg_ind_gather)
       ABI_DEALLOCATE(kg_small_gather)
     end if
   end if!End of reequilibration step for paral KGB
   ABI_DEALLOCATE(npw_gather)
 end if

!BandFFT: band load balancing
 if(mpi_enreg%paral_kgb==1.and.mpi_enreg%nproc_band>0) then
   np_band=max(1,mpi_enreg%nproc_band)
   npw_split=ig;npw=npw_split/np_band
   npw_remain=modulo(npw_split,np_band)
   if(mpi_enreg%me_band < npw_remain) npw=npw+1
   if(mpw > 0) then ! This is for the case when we only compute npw and put mpw=0
     npw_before=mpi_enreg%me_band*(npw_split/np_band)+min(npw_remain,mpi_enreg%me_band)
     kg_small(:,1:npw)=kg_small(:,npw_before+1:npw_before+npw)
     kg_ind  (  1:npw)=kg_ind  (  npw_before+1:npw_before+npw)
   end if
 end if
 if(mpw > 0) then
   do i1=1,npw
     kg(:,i1+ikg)=kg_small(:,i1)
     if (allocated(mpi_enreg%my_kgtab)) then
       ikpt_this_proc=mpi_enreg%my_kpttab(ikpt)
       mpi_enreg%my_kgtab(i1,ikpt_this_proc) = kg_ind(i1)
     end if
   end do
   ABI_DEALLOCATE(kg_small)
   ABI_DEALLOCATE(kg_ind)
 end if

 ABI_DEALLOCATE(array_ipw)

!Take care of the me_g0 flag
 if(mpi_enreg%paral_kgb==1.and.mpi_enreg%nproc_band>0) then
   if(mpi_enreg%me_band==0.and.mpi_enreg%me_g0==1) then
!    In this case, the processors had the 0 G vector before the new distribution, and still keeps it
     mpi_enreg%me_g0=1
   else
!    All other cases
     mpi_enreg%me_g0=0
   end if
 end if

!Check that npw is not zero
 if(mpi_enreg%paral_kgb==1.and.npw==0) then
   write(message,'(5a)' )&
&   'Please decrease the number of npband*npfft MPI processes!',ch10,&
&   'One of the MPI process has no plane-wave to handle.',ch10,&
&   'Action: decrease npband and/or npfft.'
   MSG_ERROR(message)
 endif

 call timab(23,2,tsec)

 DBG_EXIT("COLL")

end subroutine kpgsph
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/kpgcount
!! NAME
!! kpgcount
!!
!! FUNCTION
!!  Give the minimum and maximum number of G vectors in each direction:
!!  for each k-point, compute the number of G vectors in each direction,
!!  then store the min and max value over the set of k-points.
!!
!! INPUTS
!!  ecut=planewave kinetic energy cutoff (hartrees)
!!  exchn2n3d=if 1, n2 and n3 are exchanged
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!  istwfk(nkpt)=option parameter that describes the storage of wfs
!!  kpt(3,nkpt)=reduced coords of k point (in terms of recip latt vecs)
!!  nkpt=number of k points
!!
!! OUTPUT
!!  ngmax(3)=maximum number of G vectors in each direction (x,y,z)
!!  ngmin(3)=minimum number of G vectors in each direction (x,y,z)
!!
!! NOTES
!!  This routine has been extracted from kpgsph...
!!
!! PARENTS
!!      finddistrproc
!!
!! CHILDREN
!!      xmpi_sum,xmpi_sum_master
!!
!! SOURCE

subroutine kpgcount(ecut,exchn2n3d,gmet,istwfk,kpt,ngmax,ngmin,nkpt)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: exchn2n3d,nkpt
 integer,intent(out) :: ngmax(3),ngmin(3)
 real(dp),intent(in) :: ecut
!arrays
 integer,intent(in) :: istwfk(nkpt)
 real(dp),intent(in) :: gmet(3,3),kpt(3,nkpt)

!Local variables-------------------------------
!scalars
 integer :: ii,ikpt,istwf_k,kmax,ng1,ng2,ng3,nmin
 real(dp) :: gmet_trace,gscut,xx
!arrays
 integer :: ngrid(3),nmax(3)
 real(dp) :: minor(3),numer(3)

! *************************************************************************

 DBG_ENTER("COLL")

 gscut=0.5_dp*ecut*piinv**2
 gmet_trace=gmet(1,1)+gmet(2,2)+gmet(3,3)
 minor(1)=gmet(2,2)*gmet(3,3)-gmet(2,3)**2
 minor(2)=gmet(1,1)*gmet(3,3)-gmet(1,3)**2
 minor(3)=gmet(2,2)*gmet(1,1)-gmet(1,2)**2
 numer(1)=gmet(1,2)**2*gmet(3,3)-2.0_dp*gmet(1,2)*gmet(1,3)*gmet(2,3)+gmet(1,3)**2*gmet(2,2)
 numer(2)=gmet(2,3)**2*gmet(1,1)-2.0_dp*gmet(1,2)*gmet(1,3)*gmet(2,3)+gmet(2,1)**2*gmet(3,3)
 numer(3)=gmet(3,2)**2*gmet(1,1)-2.0_dp*gmet(1,2)*gmet(1,3)*gmet(2,3)+gmet(1,3)**2*gmet(2,2)

 ngmin(:)=1000000;ngmax(:)=0
 do ikpt=1,nkpt
   istwf_k=istwfk(ikpt)

   do ii=1,3
     xx=gmet(ii,ii)*minor(ii)-numer(ii)
     if(xx<tol10*gmet_trace**3.or.minor(ii)<tol10*gmet_trace**2)then
       MSG_BUG('The metric tensor seem incorrect')
     end if
     kmax=sqrt(gscut*minor(ii)/xx)
     nmax(ii)=floor(kmax-kpt(ii,ikpt)+tol10)
     nmin=ceiling(-kmax-kpt(ii,ikpt)-tol10)
     ngrid(ii)=nmax(ii)-nmin+1
   end do

   ng1=ngrid(1);if(istwf_k==2.or.istwf_k==3) ng1=nmax(1)+1
   if(exchn2n3d==0)then
     ng3=ngrid(3)
     ng2=ngrid(2);if(istwf_k>=2) ng2=nmax(2)+1
     if(istwf_k>=2.and.istwf_k<=5) ng2=ng2-1
   else
     ng2=ngrid(2)
     ng3=ngrid(3);if(istwf_k>=2) ng3=nmax(3)+1
     if(istwf_k==2.or.istwf_k==3.or.istwf_k==6.or.istwf_k==7) ng3=ng3-1
   end if

   if(ng1<ngmin(1)) ngmin(1)=ng1
   if(ng2<ngmin(2)) ngmin(2)=ng2
   if(ng3<ngmin(3)) ngmin(3)=ng3
   if(ng1>ngmax(1)) ngmax(1)=ng1
   if(ng2>ngmax(2)) ngmax(2)=ng2
   if(ng3>ngmax(3)) ngmax(3)=ng3

 end do

 DBG_EXIT("COLL")

end subroutine kpgcount
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/get_kg
!! NAME
!!  get_kg
!!
!! FUNCTION
!!  Helper function to calculate the set of G-vectors at a given kpoint.
!!  without taking advantage of FFT parallelism and G-vector distributions.
!!
!! INPUTS
!!  kpoint(3)=The k-point in reduced coordinates.
!!  ecut=Cutoff energy for planewave basis set.
!!  gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{-2}$).
!!  istwfk=Options defining if time-reversal is used to decrease the number of G"s.
!!
!! OUTPUT
!!  npw_k=Total number of G-vectors in the full G-sphere.
!!
!! SIDE EFFECTS
!!  kg_k(:,:):
!!   allocated inside the routine.
!!   output: kg_k(3,npw_k) contains the list of G-vectors.
!!
!! PARENTS
!!      fftprof,m_ebands,m_fft,m_fft_prof,m_gkk,m_io_kss,m_phgamma,m_phpi
!!      m_shirley,m_sigmaph,m_wfd,m_wfk,outkss
!!
!! CHILDREN
!!      xmpi_sum,xmpi_sum_master
!!
!! SOURCE

subroutine get_kg(kpoint, istwf_k, ecut, gmet, npw_k, kg_k)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k
 integer,intent(out) :: npw_k
 real(dp),intent(in) :: ecut
!arrays
 integer,allocatable,intent(out) :: kg_k(:,:)
 real(dp),intent(in) :: gmet(3,3),kpoint(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: mkmem_=1,exchn2n3d0=0,ikg0=0
 integer :: npw_k_test
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer :: kg_dum(3,0)

! *********************************************************************

 call initmpi_seq(MPI_enreg_seq)

 ! Calculate the number of G-vectors for this k-point.
 call kpgsph(ecut,exchn2n3d0,gmet,ikg0,0,istwf_k,kg_dum,kpoint,0,MPI_enreg_seq,0,npw_k)

 ! Allocate and calculate the set of G-vectors.
 ABI_MALLOC(kg_k,(3,npw_k))
 call kpgsph(ecut,exchn2n3d0,gmet,ikg0,0,istwf_k,kg_k,kpoint,mkmem_,MPI_enreg_seq,npw_k,npw_k_test)

 call destroy_mpi_enreg(MPI_enreg_seq)

end subroutine get_kg
!!***

!!****f* m_fftcore/kgindex
!! NAME
!! kgindex
!!
!! FUNCTION
!! Compute the index of each plane wave on a FFT grid.
!!
!! INPUTS
!!  kg_k(3,npw_k)=dimensionless coords of G vecs (integer)
!!  mpi_enreg=information about MPI parallelization
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  npw_k=number of planewaves
!!
!! OUTPUT
!!  indpw_k(npw_k)=linear list number (in fft box) of given G vector for the current processor (local adress)
!!                =0 if kg_k(ipw) is not treated by this procesor
!!  mask(npw_k)=True if  kg_k(ipw) belongs to this procesor, false otherwise.
!!
!! NOTES
!!   mpi_enreg is not necessary in this case (the info is also in ngfft), but much more easy to read...
!!
!! PARENTS
!!      m_fft_prof,m_gsphere,m_screening,m_shirley,m_wfd,prcref,prcref_PMA
!!
!! CHILDREN
!!
!! SOURCE

subroutine kgindex(indpw_k, kg_k, mask, mpi_enreg, ngfft, npw_k)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: kg_k(3,npw_k),ngfft(18)
 integer,intent(out) :: indpw_k(npw_k)
 logical,intent(out) :: mask(npw_k)
!Local variables-------------------------------
!scalars
 integer :: ig,ig1,ig2,ig3,me_fft,n1,n2,n3,nd2
 character(len=500) :: msg
 !arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 !integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)

! *************************************************************************

 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)

 ! Use the following indexing (N means ngfft of the adequate direction)
 ! 0 1 2 3 ... N/2    -(N-1)/2 ... -1    <= kg
 ! 1 2 3 4 ....N/2+1  N/2+2    ...  N    <= index

 me_fft=mpi_enreg%me_fft
 nd2=(n2-1)/mpi_enreg%nproc_fft+1

 if (n2== mpi_enreg%distribfft%n2_coarse) then
   fftn2_distrib => mpi_enreg%distribfft%tab_fftdp2_distrib
   ffti2_local => mpi_enreg%distribfft%tab_fftdp2_local
 else if (n2 == mpi_enreg%distribfft%n2_fine) then
   fftn2_distrib => mpi_enreg%distribfft%tab_fftdp2dg_distrib
   ffti2_local => mpi_enreg%distribfft%tab_fftdp2dg_local
 else
   MSG_BUG("Unable to find an allocated distrib for this fft grid")
 end if

 !call ptabs_fourwf(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

 do ig=1,npw_k
   ig1=modulo(kg_k(1,ig),n1)
   ig2=modulo(kg_k(2,ig),n2)
   ig3=modulo(kg_k(3,ig),n3)
   if(me_fft==fftn2_distrib(ig2+1)) then
     ig2=ffti2_local(ig2+1) - 1
     indpw_k(ig)=ig1+1+n1*(ig2+nd2*ig3)
     mask(ig)=.true.
   else
     indpw_k(ig)=0
     mask(ig)=.false.
   end if
   if (any(kg_k(:,ig) > ngfft(1:3)/2) .or. any(kg_k(:,ig) < -(ngfft(1:3)-1)/2) ) then
     write(msg,'(a,3(i0,1x),a)')" The G-vector: ",kg_k(:, ig)," falls outside the FFT box. Increase boxcutmin (?)"
     MSG_ERROR(msg)
   end if
 end do

end subroutine kgindex
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/addrho
!! NAME
!!  addrho
!!
!! FUNCTION
!!   Add the contribution to the density generated by n1dfft x-y planes
!!
!! INPUTS
!!  icplexwf=1 if u(r) is real, 2 otherwise.
!!  includelast
!!  nd1,nd2=Leading dimensions of rhopart.
!!  n2=FFT dimension along y.
!!  lot=2nd Leading dimension of zw (cache blocking factor).
!!  n1dfft=Number of 1D FFTs along y performed.
!!  zw(2,lot,n2)=Array with the x-y planes (wavefunction in real space).
!!  weigth=Weight factor for the density.
!!
!! SIDE EFFECTS
!!   rhopart(nd1,nd2)=density in the x-y plane, accumulated in output.
!!
!! PARENTS
!!
!! SOURCE

pure subroutine addrho(icplexwf,includelast,nd1,nd2,n2,lot,n1dfft,zw,rhopart,weight_r,weight_i)


!Arguments ------------------------------------
 integer,intent(in) :: icplexwf,includelast,nd1,nd2,n2,lot,n1dfft
 real(dp),intent(in) :: zw(2,lot,n2)
 real(dp),intent(inout) :: rhopart(nd1,nd2)
 real(dp),intent(in) :: weight_i,weight_r

!Local variables-------------------------------
 integer :: i2,j

! *************************************************************************

 if (icplexwf==2) then
   ! Complex wavefunction in real space.
   do i2=1,n2-1,2
     do j=1,n1dfft
       rhopart(j,i2) =   rhopart(j,i2) +   weight_r*zw(1,j,i2)**2+weight_i*zw(2,j,i2)**2
       rhopart(j,i2+1) = rhopart(j,i2+1) + weight_r*zw(1,j,i2+1)**2+weight_i*zw(2,j,i2+1)**2
     end do
   end do

   if (2*(n2/2)/=n2) then
     do j=1,n1dfft
       rhopart(j,n2  )=rhopart(j,n2  )+weight_r*zw(1,j,n2  )**2+weight_i*zw(2,j,n2  )**2
     end do
   end if
 else
   ! The wavefunction is real, in real space
   if (includelast==1) then
     do i2=1,n2
       do j=1,n1dfft
         rhopart(2*j-1,i2)=rhopart(2*j-1,i2)+zw(1,j,i2)**2*weight_r
         rhopart(2*j  ,i2)=rhopart(2*j  ,i2)+zw(2,j,i2)**2*weight_i
       end do
     end do
   else
     do i2=1,n2
       do j=1,n1dfft-1
         rhopart(2*j-1,i2)=rhopart(2*j-1,i2)+zw(1,j,i2)**2*weight_r
         rhopart(2*j  ,i2)=rhopart(2*j  ,i2)+zw(2,j,i2)**2*weight_i
       end do
       rhopart(2*n1dfft-1,i2)=rhopart(2*n1dfft-1,i2)+zw(1,n1dfft,i2)**2*weight_r
     end do
   end if

 end if

end subroutine addrho
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/multpot
!! NAME
!!  multpot
!!
!! FUNCTION
!!
!! INPUTS
!!  icplexwf=1 if u(r) is real, 2 otherwise.
!!  icplex=1 if v(r) is real, 2 otherwise.
!!  includelast
!!  nd1,nd2=Leading dimensions of pot(icplex*nd1,nd2)
!!  n2
!!  lot
!!  n1dfft
!!
!! OUTPUT
!!
!! PARENTS
!!      m_fftw3,m_sg2002
!!
!! CHILDREN
!!      xmpi_sum,xmpi_sum_master
!!
!! SOURCE

subroutine multpot(icplexwf,icplex,includelast,nd1,nd2,n2,lot,n1dfft,pot,zw)


 !Arguments ------------------------------------
 integer,intent(in) :: icplexwf,icplex,includelast,nd1,nd2,n2,lot,n1dfft
 real(dp),intent(in) :: pot(icplex*nd1,nd2)
 real(dp),intent(inout) :: zw(2,lot,n2)

!Local variables-------------------------------
 integer :: i2,j
 real(dp) :: rew,imw

! *************************************************************************

 if (icplexwf==1) then
   ! Real u(r)

   if (icplex==2) then
     MSG_BUG('multpot: icplexwf=1 and icplex=2')
   else
     ! TO BE SPEEDED UP : should use the same trick as Stefan
     if(includelast==1)then
       do i2=1,n2
         do j=1,n1dfft
           zw(1,j,i2)=zw(1,j,i2)*pot(2*j-1,i2)
           zw(2,j,i2)=zw(2,j,i2)*pot(2*j  ,i2)
         end do
       end do
     else
       do i2=1,n2
         do j=1,n1dfft-1
           zw(1,j,i2)=zw(1,j,i2)*pot(2*j-1,i2)
           zw(2,j,i2)=zw(2,j,i2)*pot(2*j  ,i2)
         end do
         zw(1,n1dfft,i2)=zw(1,n1dfft,i2)*pot(2*n1dfft-1,i2)
       end do
     end if
   end if

 else if (icplexwf==2) then
   ! Complex u(r)

   if (icplex==1) then

     do i2=1,n2-1,2
       do j=1,n1dfft
         zw(1,j,i2)=zw(1,j,i2)*pot(j,i2)
         zw(2,j,i2)=zw(2,j,i2)*pot(j,i2)
         zw(1,j,i2+1)=zw(1,j,i2+1)*pot(j,i2+1)
         zw(2,j,i2+1)=zw(2,j,i2+1)*pot(j,i2+1)
       end do
     end do

     if (2*(n2/2)/=n2) then
       do j=1,n1dfft
         zw(1,j,n2)=zw(1,j,n2)*pot(j,n2)
         zw(2,j,n2)=zw(2,j,n2)*pot(j,n2)
       end do
     end if

   else

     do i2=1,n2-1,2
       do j=1,n1dfft
         rew = zw(1,j,i2); imw = zw(2,j,i2)
         zw(1,j,i2) = rew*pot(2*j-1,i2) - imw*pot(2*j,i2)
         zw(2,j,i2) = imw*pot(2*j-1,i2) + rew*pot(2*j,i2)

         rew = zw(1,j,i2+1); imw = zw(2,j,i2+1)
         zw(1,j,i2+1) = rew*pot(2*j-1,i2+1) - imw*pot(2*j,i2+1)
         zw(2,j,i2+1) = imw*pot(2*j-1,i2+1) + rew*pot(2*j,i2+1)
       end do
     end do

     if (2*(n2/2)/=n2) then
       do j=1,n1dfft
         rew = zw(1,j,n2); imw = zw(2,j,n2)
         zw(1,j,n2) = rew*pot(2*j-1,n2) - imw*pot(2*j,n2)
         zw(2,j,n2) = imw*pot(2*j-1,n2) + rew*pot(2*j,n2)
       end do
     end if

   end if
 end if

end subroutine multpot
!!***

!----------------------------------------------------------------------

!!****f* m_fftcore/mpifft_collect_datar
!! NAME
!!  mpifft_collect_datar
!!
!! FUNCTION
!! Collect a real-space MPI-FFT distributed array on each proc.
!!
!! INPUTS
!!  ngfft(18)=contain all needed information about 3D FFT (see NOTES at beginning of scfcv)
!!  cplex=1 if real array, 2 for complex
!!  nfft=Number of FFT points treated by this MPI proc
!!  nspden=Second dimension of rhor
!!  rhor(cplex*nfft,nspden)=Array in real space (MPI-FFT distributed)
!!  fftn3_distrib(n3)=rank of the processors which own fft planes in 3rd dimension.
!!  fftn3_local(n3)=local i3 indices
!!  comm_fft=MPI-FFT communicator
!!  [master]=MPI rank, Optional. If present, the global array is available only on master node.
!!
!! OUTPUT
!!   rhor_glob(cplex*nfft_tot,nspden)=Global array
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_sum,xmpi_sum_master
!!
!! SOURCE

subroutine mpifft_collect_datar(ngfft,cplex,nfft,nspden,rhor,comm_fft,fftn3_distrib,ffti3_local,rhor_glob,master)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nspden,comm_fft
 integer,optional,intent(in) :: master
!arrays
 integer,intent(in) :: ngfft(18),fftn3_distrib(ngfft(3)),ffti3_local(ngfft(3))
 real(dp),intent(in) :: rhor(cplex*nfft,nspden)
 real(dp),intent(out) :: rhor_glob(cplex*product(ngfft(1:3)),nspden)

!Local variables-------------------------------
 integer :: ispden,i1,i2,i3,me_fft,i3_local,my_fftbase,glob_fftbase
 integer :: n1,n2,n3,ierr,nfft_tot

! *************************************************************************

 nfft_tot = product(ngfft(1:3)); me_fft = xmpi_comm_rank(comm_fft)
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)

 if (nfft_tot == nfft) then
   ! full rhor on each node, just do a copy
   rhor_glob = rhor
 else
   ! if MPI-FFT we have to gather the full rhor on each node.
   rhor_glob = zero
   do ispden=1,nspden
     do i3=1,n3
       if (me_fft == fftn3_distrib(i3)) then
         i3_local = ffti3_local(i3)
         do i2=1,n2
           my_fftbase =   cplex * ( (i2-1)*n1 + (i3_local-1)*n1*n2 )
           glob_fftbase = cplex * ( (i2-1)*n1 + (i3-1)*n1*n2 )
           do i1=1,cplex * n1
             rhor_glob(i1+glob_fftbase,ispden) = rhor(i1+my_fftbase,ispden)
           end do
         end do
       end if
     end do
   end do
   if (present(master)) then
     call xmpi_sum_master(rhor_glob,master,comm_fft,ierr)
   else
     call xmpi_sum(rhor_glob,comm_fft,ierr)
   end if
 end if

end subroutine mpifft_collect_datar
!!***

END MODULE m_fftcore
!!***
