!!****m* ABINIT/m_green
!! NAME
!!  m_green
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2025 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

! nvtx related macro definition
#include "nvtx_macros.h"

 MODULE m_green

 use defs_basis
 use m_abicore
 use m_errors

 use m_crystal, only : crystal_t
 use m_fstrings, only : int2char4
 use m_hide_lapack, only : xginv
 use m_io_tools, only : open_file
 use m_lib_four
 use m_matlu, only : add_matlu,copy_matlu,destroy_matlu,diff_matlu,init_matlu,matlu_type,print_matlu, &
                   & prod_matlu,shift_matlu,sym_matlu,trace_matlu,trace_prod_matlu,xmpi_matlu,zero_matlu
 use m_oper, only : copy_oper,copy_oper_from_ndat,copy_oper_to_ndat,destroy_oper,downfold_oper,gather_oper, &
                  & gather_oper_ks,init_oper,init_oper_ndat,inverse_oper,oper_type,print_oper,prod_oper, &
                  & trace_oper,trace_prod_oper,upfold_oper
 use m_paw_dmft, only : construct_nwli_dmft,mpi_distrib_dmft_type,paw_dmft_type
 use m_self, only : self_type
 use m_splines
 use m_time, only : timab
 use m_xmpi, only : xmpi_barrier,xmpi_sum
 use m_abi_linalg
 use, intrinsic :: iso_c_binding, only: c_size_t

#ifdef HAVE_GPU_MARKERS
 use m_nvtx_data
#endif

 implicit none

 private

 public :: init_green
 public :: destroy_green
 public :: init_green_tau
 public :: destroy_green_tau
 public :: print_green
 public :: printocc_green
 public :: compute_green
 public :: integrate_green
 public :: icip_green
 public :: fourier_green
 public :: check_fourier_green
 public :: compa_occup_ks
 public :: copy_green
 public :: occup_green_tau
 public :: add_int_fct
 public :: int_fct
 public :: fourier_fct
 public :: spline_fct
 public :: distrib_paral
 public :: greendftcompute_green
 public :: fermi_green
 public :: newton
 public :: local_ks_green
 public :: compute_moments_ks
 public :: compute_trace_moments_ks
 public :: compute_moments_loc
 public :: occup_fd
!!***

!!****t* m_green/green_type
!! NAME
!!  green_type
!!
!! FUNCTION
!!  This structured datatype contains the necessary data
!!
!! SOURCE

 type, public :: green_type ! for each atom

  integer :: dmft_nwli
  ! Linear index of the last imaginary frequency

  integer :: dmft_nwlo
  ! Number of imaginary frequencies

  integer :: dmftqmc_l
  ! Number of time slices for QMC

  integer :: fileprt_tau
  ! 1 if file created

  integer :: fileprt_w
  ! 1 if file created

  integer :: has_charge_matlu
  ! =2 if calculation of LOCAL CORRELATED occupations is done

  integer :: has_charge_matlu_prev
  ! =0 charge_matlu_prev not allocated
  ! =1 charge_matlu_prev is allocated
  ! =2 charge_matlu_prev is calculated (ie calculation of LOCAL CORRELATED occupations is done  from
  ! solver green function)

  integer :: has_charge_matlu_solver
  ! =0 charge_matlu_solver not allocated
  ! =1 charge_matlu_solver is allocated
  ! =2 charge_matlu_solver is calculated (ie calculation of LOCAL CORRELATED occupations is done from
  ! solver green function)

  integer :: has_greenmatlu_xsum
  ! =1 green%oper%matlu xsum in compute_green
  ! =0 green%oper%matlu non xsumed in compute_green
  ! used in integrate_green to checked that green function was computed in
  ! integrate_green.

  integer :: has_moments
  ! =1 if the high-frequency moments are computed

  integer :: ichargeloc_cv

  integer :: ifermie_cv

  integer :: nmoments
  ! Number of high-frequency moments which will be computed

  integer :: nw
  ! Number of frequencies

  integer :: use_oper_tau_ks
  ! 0 do not use oper_tau_ks
  ! 1 use oper_tau_ks

  character(len=4) :: w_type
  ! type of frequencies used

  !character(len=12) :: whichgreen
  ! describe the type of green function computed (DFT, DMFT, ..)

  real(dp) :: charge_ks
  ! Total charge computed from ks orbitals

  real(dp) :: ekin_imp
  ! Kinetic energy of the impurity

  real(dp) :: integral
  ! Integral of the interaction energy divided by U

  real(dp) :: trace_log
  ! Tr(log(G)) in KS space

  !integer, allocatable :: procb(:,:)

  !integer, allocatable :: proct(:,:)

  type(oper_type) :: occup
  ! Occupation in different basis

  type(oper_type) :: occup_tau
  ! Occupation in different basis

  complex(dpc) :: trace_fermie(12)
  ! Container to store useful quantities for a quick computation
  ! of the moments during the Fermi level search

  real(dp), allocatable :: charge_matlu(:,:)
  ! Total charge on correlated orbitals
! todo_ba name of charge_matlu is misleading: should be changed

  real(dp), allocatable :: charge_matlu_prev(:,:)
  ! Total charge on correlated orbitals from previous iteration

  real(dp), allocatable :: charge_matlu_solver(:,:)
  ! Total charge on correlated orbitals obtained from solver by
  ! integration over frequencies.

  real(dp), allocatable :: ecorr_qmc(:)
  ! Correlation energy for a given atom in qmc

  real(dp), allocatable :: tau(:)
  ! Value of time in imaginary space

  complex(dpc), allocatable :: trace_moments_log_ks(:)
  ! Trace of the moments of log(G)+log(iw*Id) in KS space

  complex(dpc), allocatable :: trace_moments_log_loc(:)
  ! Trace of the moments of log(G)+log(iw*Id) in local space

  type(oper_type), allocatable :: moments(:)
  ! High-frequency moments

  type(oper_type), allocatable :: oper(:)
  ! Green's function in different basis

  type(oper_type), allocatable :: oper_tau(:)
  ! Green's function in different basis

  real(dp), ABI_CONTIGUOUS pointer :: omega(:) => null()
  ! Value of frequencies

  type(mpi_distrib_dmft_type), pointer :: distrib => null()
  ! Datastructure for MPI parallelization

 end type green_type

!----------------------------------------------------------------------


CONTAINS
!!***

!!****f* m_green/init_green
!! NAME
!! init_green
!!
!! FUNCTION
!!  Allocate variables used in type green_type.
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!  energies_dmft  = datastructure for dmft energy
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  opt_oper_ksloc (optional) = option for init_oper
!!  wtype = "real" Green function will be computed for real frequencies
!!        = "imag" Green function will be computed for imaginary frequencies
!!  opt_moments = 1 : to allocate the high frequency moments
!!  opt_moments_ksloc = option to init green%moments
!!  opt_occup_ksloc = option to init green%occup
!!
!! OUTPUTS
!! green  = variable of type green_type
!!
!! SOURCE

subroutine init_green(green,paw_dmft,opt_oper_ksloc,wtype,opt_moments,opt_moments_ksloc,opt_occup_ksloc)

!Arguments ------------------------------------
 type(green_type), intent(inout) :: green
 type(paw_dmft_type), target, intent(in) :: paw_dmft
 integer, optional, intent(in) :: opt_moments,opt_moments_ksloc,opt_occup_ksloc,opt_oper_ksloc
 character(len=4), optional, intent(in) :: wtype
!Local variables ------------------------------------
 integer :: i,ifreq,mkmem,natom,nsppol,nw,optmoments
 integer :: optmoments_ksloc,optoccup_ksloc,optoper_ksloc,shift
!************************************************************************

 optoper_ksloc = 2
 if (present(opt_oper_ksloc)) optoper_ksloc = opt_oper_ksloc

 green%w_type = "imag"
 if (present(wtype)) green%w_type = wtype

 optmoments = 0
 if (present(opt_moments)) optmoments = opt_moments

 optoccup_ksloc = 3
 if (present(opt_occup_ksloc)) optoccup_ksloc = opt_occup_ksloc

 optmoments_ksloc = 3
 if (present(opt_moments_ksloc)) optmoments_ksloc = opt_moments_ksloc

 if (green%w_type == "imag") then
   nw = paw_dmft%dmft_nwlo
   green%omega => paw_dmft%omega_lo(:)
   green%distrib => paw_dmft%distrib
 else if (green%w_type == "real") then
   nw = size(paw_dmft%omega_r(:))
   green%omega => paw_dmft%omega_r(:)
   green%distrib => paw_dmft%distrib_r
 end if ! w_type

 natom  = paw_dmft%natom
 nsppol = paw_dmft%nsppol

 green%dmft_nwlo = paw_dmft%dmft_nwlo
 green%dmft_nwli = paw_dmft%dmft_nwli
 green%charge_ks = zero
 ABI_MALLOC(green%charge_matlu,(nsppol+1,natom))
 green%charge_matlu(:,:) = zero
 green%has_charge_matlu = 1
 green%has_greenmatlu_xsum = 0

 ABI_MALLOC(green%charge_matlu_solver,(nsppol+1,natom))
 green%charge_matlu_solver(:,:) = zero
 green%has_charge_matlu_solver = 1

 ABI_MALLOC(green%charge_matlu_prev,(nsppol+1,natom))
 green%charge_matlu_prev(:,:) = zero
 green%has_charge_matlu_prev = 1

 call init_oper(paw_dmft,green%occup,opt_ksloc=optoccup_ksloc)

!  build simple arrays to distribute the tasks in compute_green.
 !ABI_MALLOC(green%procb,(nw,paw_dmft%nkpt))
 !ABI_MALLOC(green%proct,(nw,0:paw_dmft%nproc-1))

 !call distrib_paral(paw_dmft%nkpt,paw_dmft%nproc,nw,nw_perproc,green%procb,green%proct)
 green%nw = nw

!  need to distribute memory over frequencies

!!  begin of temporary modificatios
! ABI_MALLOC(green%oper,(green%nw_perproc))
! do ifreq=1,green%nw_perproc
!  call init_oper(paw_dmft,green%oper(ifreq),opt_ksloc=optoper_ksloc)
! enddo
!
! do ifreq=1,green%nw
!   if(green%proct(ifreq,myproc)==1) then
!     do ikpt = 1 , paw_dmft%nkpt
!       if (green%procb(ifreq,ikpt)==myproc) then
!       endif
!     enddo ! ikpt
!   endif ! parallelisation
! enddo ! ifreq
!
!
!!  end of temporary modificatios
 ABI_MALLOC(green%oper,(nw))
 do ifreq=1,nw
   call init_oper(paw_dmft,green%oper(ifreq),opt_ksloc=optoper_ksloc)
 end do ! ifreq
 green%ichargeloc_cv = 0
 green%ifermie_cv    = 0

 if (paw_dmft%dmft_solv >= 5) then
   ABI_MALLOC(green%ecorr_qmc,(natom))
   green%ecorr_qmc(:) = zero
 end if

 green%fileprt_tau = 0
 green%fileprt_w = 0

 green%has_moments = optmoments

 green%nmoments = 0

 if (green%has_moments == 1) then
   green%nmoments = 5
   shift = green%distrib%shiftk
   mkmem = green%distrib%nkpt_mem(green%distrib%me_kpt+1)
   ABI_MALLOC(green%moments,(green%nmoments))
   call init_oper(paw_dmft,green%moments(1),nkpt=mkmem,shiftk=shift,opt_ksloc=2)
   do i=2,green%nmoments
     call init_oper(paw_dmft,green%moments(i),nkpt=mkmem,shiftk=shift,opt_ksloc=optmoments_ksloc)
   end do ! i
   if (paw_dmft%dmft_triqs_entropy == 1) then
     ABI_MALLOC(green%trace_moments_log_ks,(green%nmoments-1))
     ABI_MALLOC(green%trace_moments_log_loc,(green%nmoments-1))
   end if ! entropy
 end if ! moments

end subroutine init_green
!!***

!!****f* m_green/init_green_tau
!! NAME
!! init_green_tau
!!
!! FUNCTION
!!  Allocate variables used in type green_type.
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! OUTPUTS
!!  green  <type(green_type)>= green function data
!!
!! SOURCE

subroutine init_green_tau(green,paw_dmft)

!Arguments ------------------------------------
 type(green_type), intent(inout) :: green
 type(paw_dmft_type), intent(in) :: paw_dmft
!Local variables ------------------------------------
 integer :: itau,optksloc
!************************************************************************

 green%use_oper_tau_ks = 0
 optksloc = 2

 green%dmftqmc_l = paw_dmft%dmftqmc_l
 ABI_MALLOC(green%tau,(green%dmftqmc_l))
 do itau=1,green%dmftqmc_l
   green%tau(itau) = dble(itau-1) / dble(green%dmftqmc_l) / paw_dmft%temp
 end do

 call init_oper(paw_dmft,green%occup_tau,opt_ksloc=optksloc)

 ABI_MALLOC(green%oper_tau,(paw_dmft%dmftqmc_l))
 do itau=1,green%dmftqmc_l
   call init_oper(paw_dmft,green%oper_tau(itau),opt_ksloc=optksloc)
 end do

end subroutine init_green_tau
!!***

!!****f* m_green/destroy_green
!! NAME
!! destroy_green
!!
!! FUNCTION
!!  Deallocate green
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_green(green)

!Arguments ------------------------------------
 type(green_type), intent(inout) :: green
!Local variables-------------------------------
 integer :: i,ifreq
! *********************************************************************

 call destroy_oper(green%occup)
 if (allocated(green%oper)) then
   do ifreq=1,green%nw
     call destroy_oper(green%oper(ifreq))
   end do ! ifreq
   ABI_FREE(green%oper)
 end if ! green%oper

 ABI_SFREE(green%charge_matlu)
 green%has_charge_matlu = 0

 ABI_SFREE(green%charge_matlu_prev)
 green%has_charge_matlu_prev = 0

 ABI_SFREE(green%charge_matlu_solver)
 green%has_charge_matlu_solver = 0

 ABI_SFREE(green%trace_moments_log_ks)
 ABI_SFREE(green%trace_moments_log_loc)
 ABI_SFREE(green%ecorr_qmc)

 if (allocated(green%moments)) then
   do i=1,green%nmoments
     call destroy_oper(green%moments(i))
   end do ! i
   ABI_FREE(green%moments)
 end if  ! green%moments

 !if (allocated(green%procb))   then
 !   ABI_FREE(green%procb)
 !end if
 !if ( allocated(green%proct))   then
 !   ABI_FREE(green%proct)
 !end if
 green%distrib => null()
 green%omega => null()

end subroutine destroy_green
!!***

!!****f* m_green/destroy_green_tau
!! NAME
!! destroy_green_tau
!!
!! FUNCTION
!!  Deallocate green
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_green_tau(green)

!Arguments ------------------------------------
 type(green_type), intent(inout) :: green
! integer, optional, intent(in) :: opt_ksloc
!Local variables-------------------------------
 integer :: itau
! integer :: optksloc
! *********************************************************************
! if(present(opt_ksloc)) then
!   optksloc=opt_ksloc
! else
!   optksloc=3
! endif

 call destroy_oper(green%occup_tau)
 if (allocated(green%oper_tau)) then
   do itau=1,green%dmftqmc_l
     call destroy_oper(green%oper_tau(itau))
   end do ! itau
   ABI_FREE(green%oper_tau)
 end if ! green%oper_tau
 ABI_SFREE(green%tau)

end subroutine destroy_green_tau
!!***

!!****f* m_green/copy_green
!! NAME
!! copy_green
!!
!! FUNCTION
!!  Copy one data structure green1 into green2
!!
!! INPUTS
!!  green1  <type(green_type)>= green function data
!!  green2  <type(green_type)>= green function data
!!  opt_tw = option to specify which data to copy
!!          1: copy only green%occup_tau and green%oper_tau data
!!          2: copy only green%occup and green%oper data (frequency)
!!
!! OUTPUT
!!
!! SOURCE

subroutine copy_green(green1,green2,opt_tw)

!Arguments ------------------------------------
 type(green_type), intent(in) :: green1
 type(green_type), intent(inout) :: green2
 integer, intent(in) :: opt_tw
!Local variables-------------------------------
 integer :: i,ifreq,itau
! *********************************************************************

 if (opt_tw == 2) then
   call copy_oper(green1%occup,green2%occup)
   do ifreq=1,green1%nw
     call copy_oper(green1%oper(ifreq),green2%oper(ifreq))
   end do
   if (green1%has_moments == 1) then
     do i=1,green1%nmoments
       call copy_oper(green1%moments(i),green2%moments(i))
     end do
   end if
   ! Indicate to integrate_green that xsum has been done
   ! for matlu in compute_green.
   if (green1%has_greenmatlu_xsum == 1) green2%has_greenmatlu_xsum = 1
 else if (opt_tw == 1) then
   call copy_oper(green1%occup_tau,green2%occup_tau)
   do itau=1,green1%dmftqmc_l
     call copy_oper(green1%oper_tau(itau),green2%oper_tau(itau))
   end do ! itau
 end if ! opt_tw

end subroutine copy_green
!!***

!!****f* m_green/printocc_green
!! NAME
!! printocc_green
!!
!! FUNCTION
!!  Print occupations
!!
!! INPUTS
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  green  <type(green_type)>= green function data
!!  option= 1 :for G(w)
!!          2 :for G(tau)
!!          3 :for G(tau) and check % G(w)
!!          4
!!         <5: write diagonal part of KS occupation matrix
!!          5: for G(w)
!!          6: for G(tau)
!!          7 :for G(tau) and check % G(w)
!!         >8: write all elements of KS occup. matrix.
!!          9: for G(w)
!!  pawprtvol: flag for print
!!  opt_weissgreen = 1 for Weiss field
!!                 = 2 (default) for regular Green's function
!!  chtype: to specify the type of occupations being printed
!!
!! OUTPUT
!!
!! SOURCE

subroutine printocc_green(green,option,paw_dmft,pawprtvol,opt_weissgreen,chtype)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(green_type), intent(in) :: green
 integer, intent(in) :: option,pawprtvol
 integer, optional, intent(in) :: opt_weissgreen
 character(len=*), optional, intent(in) :: chtype
!Local variables-------------------------------
 character(len=500) :: message
 integer :: i_tau,optweissgreen
! *********************************************************************

 optweissgreen = 2
 if (present(opt_weissgreen)) optweissgreen = opt_weissgreen

 if (mod(option,4) == 1) then
   if (optweissgreen == 2) then
     if (present(chtype)) then
       write(message,'(4a)') ch10,"  == The ",trim(chtype)," occupations are  == "
     else
       write(message,'(2a)') ch10,"  == The occupations (integral of the Green's function) are  == "
     end if ! present(chtype)
   else if (optweissgreen == 1) then
     write(message,'(2a)') ch10,"  == The integrals of the Weiss function are  == "
   end if ! optweissgreen
   call wrtout(std_out,message,'COLL')
   call print_oper(green%occup,option,paw_dmft,pawprtvol)
 end if ! mod(option,4)=1

 if (mod(option,4) >= 2) then
   if (optweissgreen == 2) then
     write(message,'(2a)') ch10,"  == The occupations (value of G(tau) for tau=0-) are  == "
   else if (optweissgreen == 1) then
     write(message,'(2a)') ch10,"  == Values of G_0(tau) for tau=0- are  == "
   end if ! optweissgreen
   call wrtout(std_out,message,'COLL')
   call print_oper(green%occup_tau,option,paw_dmft,pawprtvol)
!   write(message,'(2a)') ch10," == check: occupations from Green functions are  == "
!   call wrtout(std_out,message,'COLL')
!   call print_oper(green%occup,1,paw_dmft,pawprtvol)
   if (mod(option,4) >= 3) then
     call diff_matlu("Local occup from integral of G(iw) ","Local occup from G(tau=0-) ",&
       & green%occup%matlu(:),green%occup_tau%matlu(:),paw_dmft%natom,1,tol4)
     write(message,'(2a)') ch10,&
       & '  *****  => Calculations of occupations in omega and tau spaces are coherent ****'
     call wrtout(std_out,message,'COLL')
   end if ! mod(option,4)>=3
 end if ! mod(option,4)>=2

 if (present(chtype)) then
   if (paw_dmft%prtvol >= 4 .and. &
      & (chtype == "DFT+DMFT (end of DMFT loop)" .or. chtype == "converged DMFT") &
      & .and. green%occup%has_opermatlu == 1) then
     write(message,'(4a)') ch10,"  == The DFT+DMFT occupation matrix for correlated electrons is == "
     call wrtout(ab_out,message,'COLL')
     call print_matlu(green%occup%matlu(:),paw_dmft%natom,pawprtvol,opt_ab_out=1)
     write(message,'(a)') "  "
     call wrtout(ab_out,message,'COLL')
   end if
 end if ! present(chtype)

 if (mod(option,4) >= 2) then
   i_tau = 1
   if (optweissgreen == 1) i_tau = -1
   call trace_matlu(green%occup_tau%matlu(:),paw_dmft%natom,itau=i_tau)
 end if ! mod(option,4)>=2

end subroutine printocc_green
!!***

!!****f* m_green/print_green
!! NAME
!! print_green
!!
!! FUNCTION
!!  print green function
!!
!! INPUTS
!!  char1 = character which describes the type of green function
!!  green  <type(green_type)>= green function data
!!  option=1 print local green function
!!         2 print KS green function
!!         3 print both local and KS green function
!!         4 print spectral function is green%w_type="real"
!!         5 print k-resolved spectral function is green%w_type="real"
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawprtvol = printing option
!!  opt_wt=1 print green function as a function of frequency
!!         2 print green function as a function of imaginary time
!!  opt_decim= if present, write more decimals
!!
!! OUTPUT
!!
!! SOURCE

subroutine print_green(char1,green,option,paw_dmft,opt_wt,opt_decim)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(green_type), intent(inout) :: green
 integer, intent(in) :: option
 integer, optional, intent(in) :: opt_wt,opt_decim
 character(len=*), intent(in) :: char1
!Local variables-------------------------------
 integer :: iall,iatom,ib,ifreq,ikpt,im,ispinor
 integer :: isppol,itau,lpawu,lsub,mbandc,natom
 integer :: ndim,nkpt,nspinor,nsppol,optwt,spf_unt
 integer :: spfkresolved_unt,spcorb_unt
 !real(dp) :: ima,re
 character(len=2000) :: message
 character(len=fnlen) :: tmpfil
 character(len=1) :: tag_is,tag_is2
 character(len=10) :: tag_at
 character(len=3) :: tag_ik
 integer, allocatable :: unitgreenfunc_arr(:),unitgreenloc_arr(:)
 complex(dpc), allocatable :: sf(:,:),sf_corr(:),sf2(:)
! *********************************************************************

 optwt = 1
 if (present(opt_wt)) optwt = opt_wt

 mbandc  = paw_dmft%mbandc
 natom   = paw_dmft%natom
 nkpt    = paw_dmft%nkpt
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol

!  == Print local Green Function
 if (option == 1 .or. option == 3) then
   ABI_MALLOC(unitgreenfunc_arr,(natom*nsppol*nspinor))
   iall = 0
   do iatom=1,natom
     lpawu = paw_dmft%lpawu(iatom)
     if (lpawu == -1) cycle
     ndim = 2*lpawu + 1
     call int2char4(iatom,tag_at)
     ABI_CHECK((tag_at(1:1)/='#'),'Bug: string length too short!')
     do isppol=1,nsppol
       write(tag_is,'(i1)') isppol
       do ispinor=1,nspinor
         iall = iall + 1
         write(tag_is2,'(i1)') ispinor
!       == Create names
         if (optwt == 1) then
           tmpfil = trim(paw_dmft%filapp)//'Green-'//trim(char1)//'-omega_iatom'// &
             & trim(tag_at)//'_isppol'//tag_is//'_ispinor'//tag_is2
         else
           tmpfil = trim(paw_dmft%filapp)//'Green-'//trim(char1)//'-tau_iatom'// &
             & trim(tag_at)//'_isppol'//tag_is//'_ispinor'//tag_is2
         end if ! optwt
         if (iall <= 4) then
           write(message,'(3a)') ch10,"  == Print green function on file ",tmpfil
           call wrtout(std_out,message,'COLL')
         else if (iall == 5) then
           write(message,'(3a)') ch10,"  == following values are printed in files"
           call wrtout(std_out,message,'COLL')
         end if ! iall
         unitgreenfunc_arr(iall) = 300 + iall - 1
         if ((optwt == 1 .or. green%fileprt_tau == 0) .or. (optwt == 2 .and. green%fileprt_tau == 1)) &
          & open(unit=unitgreenfunc_arr(iall),file=trim(tmpfil),status='unknown',form='formatted',position='append')

!         == Write in files
!           rewind(unitgreenfunc_arr(iall))
!           write(message,'(a,a,a,i4)') 'opened file : ', trim(tmpfil), ' unit', unitgreenfunc_arr(iall)
!           call wrtout(std_out,message,'COLL')
         write(message,'(2a)') ch10,"# New record :"
         call wrtout(unitgreenfunc_arr(iall),message,'COLL')
         if (optwt == 1) then
           do ifreq=1,green%nw
             if (present(opt_decim)) then
               write(message,'(2x,30(e23.16,2x))') green%omega(ifreq), &
                 & (green%oper(ifreq)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im+(ispinor-1)*ndim,isppol),im=1,ndim)
             else
               write(message,'(2x,30(e10.3,2x))') green%omega(ifreq), &
                 & (green%oper(ifreq)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im+(ispinor-1)*ndim,isppol),im=1,ndim)
             end if ! opt_decim
             call wrtout(unitgreenfunc_arr(iall),message,'COLL')
               !re=real(green%oper(ifreq)%matlu(iatom)%mat(1,1,isppol,ispinor,ispinor))
               !ima=aimag(green%oper(ifreq)%matlu(iatom)%mat(1,1,isppol,ispinor,ispinor))
!               write(228,*) green%omega(ifreq),re/(re**2+ima**2),ima/(re**2+ima**2)+green%omega(ifreq)
           end do ! ifreq
         else
           do itau=1,green%dmftqmc_l
             write(message,'(2x,30(e10.3,2x))') green%tau(itau), &
                & (green%oper_tau(itau)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im+(ispinor-1)*ndim,isppol),im=1,ndim)
               call wrtout(unitgreenfunc_arr(iall),message,'COLL')
           end do ! itau
           write(message,'(2x,30(e10.3,2x))') one/paw_dmft%temp,&
              & (-green%oper_tau(1)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im+(ispinor-1)*ndim,isppol)-one,im=1,ndim)
           call wrtout(unitgreenfunc_arr(iall),message,'COLL')
         end if ! optwt
         close(unitgreenfunc_arr(iall))
       end do ! ispinor
     end do ! isppol
   end do ! iatom
   ABI_FREE(unitgreenfunc_arr)
 end if ! option=1 or 3

!  == Print ks green function
 if ((option == 2 .or. option == 3) .and. green%oper(1)%has_operks == 1) then
   ABI_MALLOC(unitgreenloc_arr,(nsppol*nkpt))
   iall = 0
   do isppol=1,nsppol
     write(tag_is,'(i1)') isppol
     do ikpt=1,nkpt
       write(tag_ik,'(i3)') ikpt
!      do ib1 = 1, mbandc
       iall = iall + 1
!         == Create names
       if (optwt == 1) then
         tmpfil = trim(paw_dmft%filapp)//'Green-'//trim(char1)//'-omega_isppol'//tag_is//'_ikpt'//trim(adjustl(tag_ik))
       else
         tmpfil = trim(paw_dmft%filapp)//'Green-'//trim(char1)//'-tau_isppol'//tag_is//'_ikpt'//trim(adjustl(tag_ik))
       end if ! optwt
       if (iall <= 4) then
         write(message,'(3a)') ch10,"  == Print green function on file ",tmpfil
         call wrtout(std_out,message,'COLL')
       else if (iall == 5)  then
         write(message,'(3a)') ch10,"  == following values are printed in files"
         call wrtout(std_out,message,'COLL')
       end if ! iall
       unitgreenloc_arr(iall) = 400 + iall - 1
       open(unit=unitgreenloc_arr(iall),file=trim(tmpfil),status='unknown',form='formatted')
!      rewind(unitgreenloc_arr(iall))
!       write(message,'(a,a,a,i4)') 'opened file : ', trim(tmpfil), ' unit', unitgreenloc_arr(iall)
!       call wrtout(std_out,message,'COLL')

!      == Write in files
       write(message,'(2a)') ch10,"# New record : First 20 bands"
       call wrtout(unitgreenloc_arr(iall),message,'COLL')
!       call flush(std_out)
       do lsub=1,mbandc/20+1
         if (optwt == 1) then
           do ifreq=1,green%nw
!             call flush(std_out)
             write(message,'(2x,50(e10.3,2x))') green%omega(ifreq), &
               & (green%oper(ifreq)%ks(ib,ib,ikpt,isppol),ib=20*(lsub-1)+1,min(20*lsub,mbandc))
             call wrtout(unitgreenloc_arr(iall),message,'COLL')
             !re  = dble(green%oper(ifreq)%ks(1,1,ikpt,isppol))
             !ima = aimag(green%oper(ifreq)%ks(1,1,ikpt,isppol))
             !if (ikpt == 1) write(229,*) green%omega(ifreq),re/(re**2+ima**2),ima/(re**2+ima**2)+green%omega(ifreq)
           end do ! ifreq
         else if (green%use_oper_tau_ks == 1) then
           do itau=1,green%dmftqmc_l
!               call flush(std_out)
             write(message,'(2x,50(e10.3,2x))') green%tau(itau), &
                 & (green%oper_tau(itau)%ks(ib,ib,ikpt,isppol),ib=20*(lsub-1)+1,min(20*lsub,mbandc))
             call wrtout(unitgreenloc_arr(iall),message,'COLL')
           end do ! itau
         end if ! optwt
         if (20*lsub < mbandc) write(message,'(2a,i5,a,i5)') ch10,"# Same record, Following bands : From ", &
           & 20*(lsub),"  to ",min(20*(lsub+1),mbandc)
         call wrtout(unitgreenloc_arr(iall),message,'COLL')
       end do ! lsub
       close(unitgreenloc_arr(iall))
     end do ! ikpt
   end do ! isppol
   ABI_FREE(unitgreenloc_arr)
 end if ! option=2 or 3

 if ((green%w_type == "real" .and. option >= 4) .and. green%oper(1)%has_operks == 1) then
   write(message,'(2a)') ch10,"  == About to print spectral function"
   call wrtout(std_out,message,'COLL')
   if (option == 4) then
     tmpfil = trim(paw_dmft%filapp)//'_SpFunc-'//trim(char1)
     if (open_file(tmpfil,message,newunit=spf_unt,status='unknown',form='formatted') /= 0) &
       & ABI_ERROR(message)
     write(spf_unt,'(3a)') "# This is the total spectral function (DOS).",ch10, &
                         & "#    Real frequency (Ha)        Spectral function"

   end if ! option=4
   if (option == 5) then
     tmpfil = trim(paw_dmft%filapp)//'_DFTDMFT_SpectralFunction_kres' !//trim(char1)
     if (open_file(tmpfil,message,newunit=spfkresolved_unt,status='unknown',form='formatted') /= 0) &
       & ABI_ERROR(message)
     write(spfkresolved_unt,'(3a)') "# This is the k-resolved spectral function.",ch10, &
                                & "#  Real frequency (eV)       Spectral function (eV^-1)   ikpt"
     ABI_MALLOC(sf,(nkpt,green%nw))
     sf(:,:) = czero
     do ifreq=1,green%nw
       do isppol=1,nsppol
         do ikpt=1,nkpt
           do ib=1,mbandc
             sf(ikpt,ifreq) = sf(ikpt,ifreq) + green%oper(ifreq)%ks(ib,ib,ikpt,isppol)
           end do ! ib
         end do ! ikpt
       end do ! isppol
     end do ! ifreq
     do ikpt=1,nkpt
       do ifreq=1,green%nw
         write(message,'(2x,2(es24.16e3,2x),i5)') green%omega(ifreq)*Ha_eV,(-aimag(sf(ikpt,ifreq)))/pi/Ha_eV,ikpt
         call wrtout(spfkresolved_unt,message,'COLL')
       end do ! ifreq
       write(message,*)
       call wrtout(spfkresolved_unt,message,'COLL')
     end do ! ikpt
     write(message,*) ch10
     call wrtout(spfkresolved_unt,message,'COLL')
     ABI_FREE(sf)
     close(spfkresolved_unt)
!
!     do isppol = 1 , nsppol
!       do ikpt = 1, nkpt
!         do ib=1,mbandc
!           sf=czero
!           write(71,*)
!           write(71,*)  "#", ikpt, ib
!           do ifreq=1,green%nw
!             sf(ifreq)=sf(ifreq)+green%oper(ifreq)%ks(isppol,ikpt,ib,ib)
!             write(71,*) green%omega(ifreq)*Ha_eV,(-aimag(sf(ifreq)))/pi/Ha_eV,ikpt
!           enddo
!         enddo
!       enddo
!     enddo
   end if ! option=5

   if (option == 4) then
     ABI_MALLOC(sf2,(green%nw))
     sf2(:) = czero
     do ifreq=1,green%nw
       do isppol=1,nsppol
         do ikpt=1,nkpt
           do ib=1,mbandc
             sf2(ifreq) = sf2(ifreq) + green%oper(ifreq)%ks(ib,ib,ikpt,isppol)*green%oper(1)%wtk(ikpt)
           end do ! ib
         end do ! ikpt
       end do ! isppol
     end do ! ifreq
     do ifreq=1,green%nw
       write(message,'(2x,2(es24.16e3,2x))') green%omega(ifreq),(-aimag(sf2(ifreq)))/pi
       call wrtout(spf_unt,message,'COLL')
     end do ! ifreq
     ABI_FREE(sf2)
     close(spf_unt)
   end if ! option=4

   if (paw_dmft%dmft_kspectralfunc == 1) then
     ABI_MALLOC(sf_corr,(green%nw))
     do iatom=1,natom
       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       sf_corr(:) = czero
       call int2char4(iatom,tag_at)
       ABI_CHECK((tag_at(1:1)/='#'),'Bug: string length too short!')
       tmpfil = trim(paw_dmft%filapp)//'_DFTDMFT_SpFunloc_iatom'//trim(tag_at)
       if (open_file(tmpfil,message,newunit=spcorb_unt,status='unknown',form='formatted') /= 0) &
         & ABI_ERROR(message)
       ndim = 2*lpawu + 1
       write(message,'(3a,3(i3),i6,2a,i2,2a)') "# nspinor,nsppol,ndim,nw",ch10,"#",nspinor,nsppol,ndim,green%nw,ch10, &
                                          & "# lpawu",lpawu,ch10,"#       Frequency (Ha.)        Spectral function"
       call wrtout(spcorb_unt,message,'COLL')
       ndim = nspinor * ndim
       do ifreq=1,green%nw
         do isppol=1,nsppol
           do im=1,ndim
             sf_corr(ifreq) = sf_corr(ifreq) + green%oper(ifreq)%matlu(iatom)%mat(im,im,isppol)
           end do ! im
         end do ! isppol
       end do ! ifreq
       do ifreq=1,green%nw
         write(message,'(2x,2(es24.16e3,2x))') green%omega(ifreq),(-aimag(sf_corr(ifreq)))/pi
         call wrtout(spcorb_unt,message,'COLL')
       end do ! ifreq
       close(spcorb_unt)
     end do ! iatom
     ABI_FREE(sf_corr)
   end if ! dmft_kspectralfunc=1
 end if ! (green%w_type == "real" .and. option >= 4) .and. green%oper(1)%has_operks == 1

 if (optwt == 2 .and. (option == 1 .or. option == 3)) green%fileprt_tau = 1  ! file for G(tau) has been created here

end subroutine print_green
!!***

!!****f* m_green/compute_green_batched_core
!! NAME
!! compute_green_batched_core
!!
!! FUNCTION
!! Variant for core loop of compute green function, useful for computation on GPU
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  prtopt : option for printing
!!  self <type(self_type)>= variables related to self-energy
!!  opt_self = optional argument, if =1, upfold self-energy
!!  opt_log = if 1, compute Tr(log(G))
!!
!! OUTPUT
!!
!! SOURCE
subroutine compute_green_batched_core(green,paw_dmft,self,optself,optlog)

 use m_abi_linalg, only : abi_xgemm
 use m_matlu, only : add_matlu,sym_matlu
 use m_oper, only : downfold_oper,inverse_oper,upfold_oper
 use m_time, only : timab

!Arguments ------------------------------------
 type(green_type), target, intent(inout) :: green
 type(paw_dmft_type), target, intent(in) :: paw_dmft
 type(self_type), intent(inout) :: self
 integer, intent(in) :: optlog,optself
!Local variables-------------------------------
 logical :: oper_ndat_allocated
 integer :: diag,ib,ib1,ifreq,ikpt,info,isppol,lwork,mbandc,ifreq_beg,ifreq_end,idat
 integer :: me_kpt,mkmem,myproc,natom,nkpt,nmoments,nspinor,nsppol,gpu_option,ndat
 integer :: option,shift,shift_green,spacecomm,optoper_ksloc
 real(dp) :: fermilevel,wtk,temp
 complex(dpc) :: green_tmp,trace_tmp
 real(dp), allocatable :: eig(:),rwork(:),fac(:)
 complex(dpc), allocatable :: mat_tmp(:,:),work(:),omega_current(:)
 type(oper_type), target :: green_oper_ndat
 real(dp), ABI_CONTIGUOUS pointer :: eigen_dft(:,:,:)
 complex(dpc), ABI_CONTIGUOUS pointer :: ks(:,:,:,:),occup_ks(:,:,:,:)
! *********************************************************************

 ABI_NVTX_START_RANGE(NVTX_DMFT_COMPUTE_GREEN_BATCHED)

 diag = 1 - optself

! Initialise spaceComm, myproc, and nproc
 me_kpt = green%distrib%me_kpt
 myproc = paw_dmft%myproc
 spacecomm = paw_dmft%spacecomm
 oper_ndat_allocated = .false.

! Initialise integers
 mbandc  = paw_dmft%mbandc
 natom   = paw_dmft%natom
 nkpt    = paw_dmft%nkpt
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol
 temp    = paw_dmft%temp
 gpu_option = paw_dmft%gpu_option

 if (optlog == 1) then
   ABI_MALLOC(eig,(mbandc))
   ABI_MALLOC(work,(2*mbandc-1))
   ABI_MALLOC(rwork,(3*mbandc-2))
   ABI_MALLOC(mat_tmp,(mbandc,mbandc))
   call zheev("n","u",mbandc,mat_tmp(:,:),mbandc,eig(:),work(:),-1,rwork(:),info)
   lwork = int(work(1))
   ABI_FREE(work)
   ABI_MALLOC(work,(lwork))
   green%trace_log = zero
 end if ! optlog=1

 fermilevel = paw_dmft%fermie

 ! option for downfold_oper
 option = 1
 if (diag == 1) option = 3

! ================================
! == Compute Green function G(k)
! ================================
 green%occup%ks(:,:,:,:) = czero

 nmoments = 1
 if (green%has_moments == 1) then
   nmoments = green%nmoments
 end if ! moments

 shift = green%distrib%shiftk
 mkmem = green%distrib%nkpt_mem(green%distrib%me_kpt+1)
 shift_green = shift
 ndat = green%distrib%nw_mem_kptparal(green%distrib%me_freq+1)
 if (green%oper(1)%has_operks == 0) shift_green = 0

 if(green%oper(1)%has_opermatlu==0) optoper_ksloc=1
 if(green%oper(1)%has_opermatlu==1) optoper_ksloc=3

 do ifreq=1,green%nw
   if (green%distrib%proct(ifreq) == green%distrib%me_freq) then
     ifreq_beg = ifreq; ifreq_end = ifreq_beg + ndat - 1
     exit
   end if
 end do

 ABI_MALLOC(omega_current,(green%nw))
 ABI_MALLOC(fac,(green%nw))
! Initialise for compiler
 omega_current = czero
 do ifreq=1,green%nw
   !if(present(iii)) write(6,*) ch10,'ifreq  self', ifreq,self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)
!   ====================================================
!   First Upfold self-energy and double counting  Self_imp -> self(k)
!   ====================================================
!   if(mod(ifreq-1,nproc)==myproc) then
!    write(6,*) "compute_green ifreq",ifreq, mod(ifreq-1,nproc)==myproc,proct(ifreq,myproc)==1
   if (green%distrib%proct(ifreq) /= green%distrib%me_freq) cycle
   if (green%w_type == "imag") then
     omega_current(ifreq) = cmplx(zero,green%omega(ifreq),kind=dp)
     fac(ifreq) = two * paw_dmft%wgt_wlo(ifreq) * paw_dmft%temp
   else if (green%w_type == "real") then
     omega_current(ifreq) = cmplx(green%omega(ifreq),paw_dmft%temp,kind=dp)
   end if ! green%w_type

   if (green%oper(ifreq)%has_operks == 0) then
     ABI_MALLOC(green%oper(ifreq)%ks,(mbandc,mbandc,mkmem,nsppol))
     green%oper(ifreq)%nkpt   = mkmem
     green%oper(ifreq)%paral  = 1
     green%oper(ifreq)%shiftk = shift
   end if
   if(.not. oper_ndat_allocated) then
     call init_oper_ndat(paw_dmft,green_oper_ndat,ndat,nkpt=green%oper(ifreq)%nkpt,opt_ksloc=optoper_ksloc,gpu_option=gpu_option)
     if (green%oper(ifreq)%has_operks == 0) then
       green_oper_ndat%paral  = 1
       green_oper_ndat%shiftk = shift
     end if
     oper_ndat_allocated=.true.
   end if
 end do ! ifreq

 if (optself == 0) then
   if(green_oper_ndat%gpu_option==ABI_GPU_DISABLED) then
     green_oper_ndat%ks(:,:,:,:) = czero
   else if(green_oper_ndat%gpu_option==ABI_GPU_OPENMP) then
     call gpu_set_to_zero_complex(green_oper_ndat%ks, int(nsppol,c_size_t)*ndat*mbandc*mbandc*mkmem)
   end if
   call copy_oper_to_ndat(green%oper,green_oper_ndat,ndat,green%nw,green%distrib%proct,green%distrib%me_freq,.false.)
 else
   do ifreq=ifreq_beg,ifreq_end
     call add_matlu(self%hdc%matlu(:),self%oper(ifreq)%matlu(:),green%oper(ifreq)%matlu(:),natom,-1)
   end do
   call copy_oper_to_ndat(green%oper,green_oper_ndat,ndat,green%nw,green%distrib%proct,green%distrib%me_freq,.false.)
   call upfold_oper(green_oper_ndat,paw_dmft,procb=green%distrib%procb(:),iproc=me_kpt,gpu_option=gpu_option)
 end if ! optself

 ks => green_oper_ndat%ks
 eigen_dft => paw_dmft%eigen_dft
 if(gpu_option==ABI_GPU_DISABLED) then
   do ifreq=ifreq_beg,ifreq_end
     do isppol=1,nsppol
       do ikpt=1,mkmem
         do ib=1,mbandc
           idat=ifreq-(ifreq_end-ndat)
           green_tmp = omega_current(ifreq) + fermilevel - eigen_dft(ib,ikpt+shift,isppol)
           if (optself == 0) then
             ks(ib,ib+(idat-1)*mbandc,ikpt+shift_green,isppol) = cone / green_tmp
           else
             ks(ib,ib+(idat-1)*mbandc,ikpt+shift_green,isppol) = &
                 & ks(ib,ib+(idat-1)*mbandc,ikpt+shift_green,isppol) + green_tmp
           end if
         end do ! ib
       end do ! ikpt
     end do ! isppol
   end do ! ifreq
 else if(gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
   if (optself == 0) then
     !$OMP TARGET TEAMS DISTRIBUTE MAP(to:ks,eigen_dft) PRIVATE(ifreq,idat)
     do ifreq=ifreq_beg,ifreq_end
       idat=ifreq-(ifreq_end-ndat)
       !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(isppol,ikpt,ib,green_tmp)
       do isppol=1,nsppol
         do ikpt=1,mkmem
           do ib=1,mbandc
             green_tmp = omega_current(ifreq) + fermilevel - eigen_dft(ib,ikpt+shift,isppol)
             ks(ib,ib+(idat-1)*mbandc,ikpt+shift_green,isppol) = cone / green_tmp
           end do ! ib
         end do ! ikpt
       end do ! isppol
     end do ! ifreq
   else
     !$OMP TARGET TEAMS DISTRIBUTE MAP(to:ks,eigen_dft) PRIVATE(ifreq,idat)
     do ifreq=ifreq_beg,ifreq_end
       idat=ifreq-(ifreq_end-ndat)
       !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(isppol,ikpt,ib,green_tmp)
       do isppol=1,nsppol
         do ikpt=1,mkmem
           do ib=1,mbandc
             green_tmp = omega_current(ifreq) + fermilevel - eigen_dft(ib,ikpt+shift,isppol)
             ks(ib,ib+(idat-1)*mbandc,ikpt+shift_green,isppol) = &
                 & ks(ib,ib+(idat-1)*mbandc,ikpt+shift_green,isppol) + green_tmp
           end do ! ib
         end do ! ikpt
       end do ! isppol
     end do ! ifreq
   end if
#endif
 end if

 if (optself /= 0) then
   call inverse_oper(green_oper_ndat,1,procb=green%distrib%procb(:),iproc=me_kpt,gpu_option=gpu_option)
 end if

 if (optlog == 1) then
   call copy_oper_from_ndat(green_oper_ndat,green%oper,ndat,green%nw,green%distrib%proct,green%distrib%me_freq,.true.)
   do ifreq=1,green%nw
     if (green%distrib%proct(ifreq) /= green%distrib%me_freq) cycle
     trace_tmp = czero
     do isppol=1,nsppol
       do ikpt=1,mkmem
         wtk = green%oper(ifreq)%wtk(ikpt+shift)
         if (optself == 0) then
           do ib=1,mbandc
             trace_tmp = trace_tmp + two*temp*wtk*log(green%oper(ifreq)%ks(ib,ib,ikpt+shift_green,isppol)*omega_current(ifreq))
           end do ! ib
         else
           ! Use Tr(log(G(iw))) + Tr(log(G(-iw))) = Tr(log(G(iw)G(iw)^H)) so that we can use the much faster zheev instead of
           ! zgeev. Even if G(iw) and G(iw)^H have no reason to commute, this equality still holds since we only care about the trace.
           call abi_xgemm("n","c",mbandc,mbandc,mbandc,cone,green%oper(ifreq)%ks(:,:,ikpt+shift_green,isppol), &
                        & mbandc,green%oper(ifreq)%ks(:,:,ikpt+shift_green,isppol),mbandc,czero,mat_tmp(:,:),mbandc)
           call zheev("n","u",mbandc,mat_tmp(:,:),mbandc,eig(:),work(:),lwork,rwork(:),info)
           ! Do not use DOT_PRODUCT
           trace_tmp = trace_tmp + sum(log(eig(:)*omega_current(ifreq)))*wtk*temp
         end if ! optself
       end do ! ikpt
     end do ! isppol
     if (nsppol == 1 .and. nspinor == 1) trace_tmp = trace_tmp * two
     green%trace_log = green%trace_log + dble(trace_tmp)
   end do ! ifreq
   call copy_oper_to_ndat(green%oper,green_oper_ndat,ndat,green%nw,green%distrib%proct,green%distrib%me_freq,.true.)
 end if ! optlog=1

 if (green%w_type /= "real") then
   ABI_NVTX_START_RANGE(NVTX_DMFT_ADD_INT_FCT)
   occup_ks    => green%occup%ks
   ks          => green_oper_ndat%ks
   if(gpu_option==ABI_GPU_DISABLED) then
     do ifreq=ifreq_beg,ifreq_end

       do isppol=1,nsppol
         do ikpt=1,mkmem
           do ib1=1,mbandc
             idat=ifreq-(ifreq_end-ndat)
             if (diag == 1) then
               green%occup%ks(ib1,ib1,ikpt+shift,isppol) = green%occup%ks(ib1,ib1,ikpt+shift,isppol) + &
                  & fac(ifreq)*green_oper_ndat%ks(ib1,ib1+(idat-1)*mbandc,ikpt+shift_green,isppol)
             else
               do ib=1,mbandc
                 green%occup%ks(ib,ib1,ikpt+shift,isppol) = green%occup%ks(ib,ib1,ikpt+shift,isppol) + &
                   & fac(ifreq)*green_oper_ndat%ks(ib,ib1+(idat-1)*mbandc,ikpt+shift_green,isppol)
               end do ! ib
             end if ! diag
           end do ! ib1
         end do ! ikpt
       end do ! isppol

     end do ! ifreq
   else if(gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
     if (diag == 1) then
       !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) &
       !$OMP& MAP(tofrom:occup_ks) MAP(to:ks,fac) PRIVATE(isppol,ikpt)
       do isppol=1,nsppol
         do ikpt=1,mkmem
           !$OMP PARALLEL DO PRIVATE(ib1,ifreq,idat)
           do ib1=1,mbandc
             do ifreq=ifreq_beg,ifreq_end
               idat=ifreq-(ifreq_end-ndat)
               occup_ks(ib1,ib1,ikpt+shift,isppol) = occup_ks(ib1,ib1,ikpt+shift,isppol) + &
                 & fac(ifreq)*ks(ib1,ib1+(idat-1)*mbandc,ikpt+shift_green,isppol)
             end do ! ifreq
           end do ! ib1
         end do ! ikpt
       end do ! isppol

     else

       !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(3) &
       !$OMP& MAP(tofrom:occup_ks) MAP(to:ks,fac) PRIVATE(isppol,ikpt,ib1)
       do isppol=1,nsppol
         do ikpt=1,mkmem
           do ib1=1,mbandc
             !$OMP PARALLEL DO PRIVATE(ib,ifreq,idat)
             do ib=1,mbandc
               do ifreq=ifreq_beg,ifreq_end
                 idat=ifreq-(ifreq_end-ndat)
                 occup_ks(ib,ib1,ikpt+shift,isppol) = occup_ks(ib,ib1,ikpt+shift,isppol) + &
                   & fac(ifreq)*ks(ib,ib1+(idat-1)*mbandc,ikpt+shift_green,isppol)
               end do ! ib
             end do ! ib1
           end do ! ikpt
         end do ! isppol
       end do ! ifreq

     end if ! diag
#endif
   end if
   ABI_NVTX_END_RANGE()
 end if ! w_type/="real"

 if (paw_dmft%lchipsiortho == 1 .or. optself == 1) then
   call downfold_oper(green_oper_ndat,paw_dmft,&
   &    procb=green%distrib%procb(:),iproc=me_kpt,option=option,gpu_option=gpu_option)
 end if ! lchipsiortho=1

 call copy_oper_from_ndat(green_oper_ndat,green%oper,ndat,green%nw,green%distrib%proct,&
 &    green%distrib%me_freq,green%oper(ifreq_beg)%has_operks /= 0)

 do ifreq=1,green%nw
   if (green%distrib%proct(ifreq) /= green%distrib%me_freq) cycle
   if (green%oper(ifreq)%has_operks == 0) then
     ABI_FREE(green%oper(ifreq)%ks)
   end if
 end do ! ifreq

 call destroy_oper(green_oper_ndat)
 ABI_FREE(omega_current)
 ABI_FREE(fac)

 ABI_NVTX_END_RANGE()

end subroutine compute_green_batched_core
!!***

!!****f* m_green/compute_green
!! NAME
!! compute_green
!!
!! FUNCTION
!! compute green function from DFT and self-energy
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  prtopt : option for printing
!!  self <type(self_type)>= variables related to self-energy
!!  opt_self = optional argument, if =1, upfold self-energy
!!  opt_nonxsum = 0 : do usual xsum after calculation of green(freq)%ks
!!              = 1 : do not do xsum after calculation of green(freq)%ks: each proc as only a part
!!                    of the data: this is useful where only the total number of electron will be computed.
!!  opt_nonxsum2 0 : green(ifreq)%matlu will be broadcasted
!!               1 : green(ifreq)%matlu will not be broadcasted in compute_green: calc
!!                   if occupations will not possible.
!!                   (a keyword: compute_local_green would be in fact equivalent and more clear)
!!  opt_log = if 1, compute Tr(log(G))
!!  opt_restart_moments = if 1, activate quick restart for calculation of the high frequency moments
!!                        (useful when compute_green is called after fermi_green)
!!
!! OUTPUT
!!
!! SOURCE

subroutine compute_green(green,paw_dmft,prtopt,self,opt_self,opt_nonxsum,opt_nonxsum2,opt_log,opt_restart_moments)

!Arguments ------------------------------------
 type(green_type), intent(inout) :: green
 type(paw_dmft_type), intent(in) :: paw_dmft
 !type(MPI_type), intent(in) :: mpi_enreg
 type(self_type), intent(inout) :: self
 integer, intent(in) :: prtopt
 integer, optional, intent(in) :: opt_log,opt_nonxsum,opt_nonxsum2,opt_restart_moments,opt_self
!Local variables-------------------------------
 !logical :: lintegrate
 integer :: band_index,diag,i,ib,ierr,ifreq,ikpt,info,isppol,lwork,mbandc
 integer :: me_kpt,mkmem,myproc,natom,nband_k,nkpt,nmoments,nspinor,nsppol
 integer :: opt_quick_restart,option,optlog,optnonxsum,optnonxsum2,optself
 integer :: shift,shift_green,spacecomm,gpu_option
 real(dp) :: beta,correction,eigen,fac,fermilevel,freq2,temp,wtk
 complex(dpc) :: green_tmp,omega_current,trace_tmp
 character(len=500) :: message
 real(dp) :: tsec(2)
 real(dp), allocatable :: eig(:),rwork(:)
 complex(dpc), allocatable :: mat_tmp(:,:),omega_fac(:),work(:)
#ifdef HAVE_OPENMP_OFFLOAD
 integer :: ndat
 type(oper_type), target :: green_oper_ndat
#endif
! integer, allocatable :: procb(:,:),proct(:,:)
! *********************************************************************

 ABI_NVTX_START_RANGE(NVTX_DMFT_COMPUTE_GREEN)
 !lintegrate=.true.
 !if(lintegrate.and.green%w_type=="real") then
 !if(green%w_type=="real") then
 !  message = 'integrate_green not implemented for real frequency'
 !  ABI_BUG(message)
 !endif
 call timab(624,1,tsec(:))
 optself = 0
 if (present(opt_self)) optself = opt_self
 optnonxsum = 0
 if (present(opt_nonxsum)) optnonxsum = opt_nonxsum
 optnonxsum2 = 0
 if (present(opt_nonxsum2)) optnonxsum2 = opt_nonxsum2
 optlog = 0
 if (present(opt_log)) optlog = opt_log
 opt_quick_restart = 0
 if (present(opt_restart_moments)) opt_quick_restart = opt_restart_moments

 diag = 1 - optself

 if (prtopt > 0) then
   write(message,'(2a)') ch10," ===  Compute Green's function "
   call wrtout(std_out,message,'COLL')
 end if ! prtopt>0

 if (self%nw /= green%nw) then
   message = ' BUG: frequencies for green and self not coherent'
   ABI_BUG(message)
 end if

! Initialise spaceComm, myproc, and nproc
 me_kpt = green%distrib%me_kpt
 myproc = paw_dmft%myproc
 !nproc     = paw_dmft%nproc
 spacecomm = paw_dmft%spacecomm

! Initialise integers
 !mband   = paw_dmft%mband
 mbandc  = paw_dmft%mbandc
 natom   = paw_dmft%natom
 nkpt    = paw_dmft%nkpt
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol
 temp    = paw_dmft%temp
 gpu_option = paw_dmft%gpu_option

 if (optlog == 1) then
   ABI_MALLOC(eig,(mbandc))
   ABI_MALLOC(work,(2*mbandc-1))
   ABI_MALLOC(rwork,(3*mbandc-2))
   ABI_MALLOC(mat_tmp,(mbandc,mbandc))
   call zheev("n","u",mbandc,mat_tmp(:,:),mbandc,eig(:),work(:),-1,rwork(:),info)
   lwork = int(work(1))
   ABI_FREE(work)
   ABI_MALLOC(work,(lwork))
   green%trace_log = zero
 end if ! optlog=1

 !icomp_chloc = 0

 option = 1
 fermilevel = paw_dmft%fermie
 if (option == 123) then
   fermilevel = two
   write(message,'(2a,e14.3,a)') ch10,'  Warning (special case for check: fermi level=',fermilevel,')'
   call wrtout(std_out,message,'COLL')
 end if ! option=123

 option = merge(3,1,diag==1) ! option for downfold_oper

! ====================================================
! Upfold self-energy and double counting  Self_imp -> self(k)
! ====================================================
! if(optself==1) then
!   do ifreq=1,green%nw
!     call upfold_oper(self%oper(ifreq),paw_dmft,1)
!   enddo ! ifreq
!   call upfold_oper(self%hdc,paw_dmft,1)
! endif

! ================================
! == Compute Green function G(k)
! ================================
 green%occup%ks(:,:,:,:) = czero

 nmoments = 1
 if (green%has_moments == 1) then
   nmoments = green%nmoments
   call compute_moments_ks(green,self,paw_dmft,opt_self=optself,opt_log=optlog, &
                         & opt_quick_restart=opt_quick_restart)
   if (paw_dmft%lchipsiortho == 1 .or. optself == 1) then
     call downfold_oper(green%moments(1),paw_dmft,option=2)
     do i=2,green%nmoments
       call downfold_oper(green%moments(i),paw_dmft,option=option)
     end do ! i
     do i=1,green%nmoments
       call xmpi_matlu(green%moments(i)%matlu(:),natom,green%distrib%comm_kpt)
       call sym_matlu(green%moments(i)%matlu(:),paw_dmft)
     end do ! i
   end if ! lchipsiortho=1
 end if ! moments

 if (green%w_type /= "real") then
   ABI_MALLOC(omega_fac,(nmoments))
   do i=1,nmoments
     omega_fac(i) = czero
     do ifreq=green%nw,1,-1 ! NEVER change the summation order and DON'T use the intrinsic SUM
       omega_fac(i) = omega_fac(i) + paw_dmft%wgt_wlo(ifreq) / (paw_dmft%omega_lo(ifreq))**i
     end do
     omega_fac(i) = - two * temp * omega_fac(i) / (j_dpc)**i
     if (i == 1) omega_fac(i) = omega_fac(i) + half
     if (i == 2) omega_fac(i) = omega_fac(i) - cone/(four*temp)
     if (i == 4) omega_fac(i) = omega_fac(i) + cone/(dble(48)*(temp**3))
   end do ! i
 end if ! w_type

 shift = green%distrib%shiftk
 mkmem = green%distrib%nkpt_mem(green%distrib%me_kpt+1)
 shift_green = merge(0,shift,green%oper(1)%has_operks==0)


 if(mkmem>0 .and. gpu_option==ABI_GPU_OPENMP) then
   call compute_green_batched_core(green,paw_dmft,self,optself,optlog)
 else if(mkmem>0) then
   ABI_NVTX_START_RANGE(NVTX_DMFT_COMPUTE_GREEN_LOOP)
! Initialise for compiler
 omega_current = czero
 do ifreq=1,green%nw
   !if(present(iii)) write(6,*) ch10,'ifreq  self', ifreq,self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)
!   ====================================================
!   First Upfold self-energy and double counting  Self_imp -> self(k)
!   ====================================================
!   if(mod(ifreq-1,nproc)==myproc) then
!    write(6,*) "compute_green ifreq",ifreq, mod(ifreq-1,nproc)==myproc,proct(ifreq,myproc)==1
   if (green%distrib%proct(ifreq) /= green%distrib%me_freq) cycle
   if (green%w_type == "imag") then
     omega_current = cmplx(zero,green%omega(ifreq),kind=dp)
     fac = two * paw_dmft%wgt_wlo(ifreq) * temp
   else if (green%w_type == "real") then
     omega_current = cmplx(green%omega(ifreq),paw_dmft%temp,kind=dp)
   end if ! green%w_type

   if (green%oper(ifreq)%has_operks == 0) then
     ABI_MALLOC(green%oper(ifreq)%ks,(mbandc,mbandc,mkmem,nsppol))
     green%oper(ifreq)%nkpt   = mkmem
     green%oper(ifreq)%paral  = 1
     green%oper(ifreq)%shiftk = shift
   end if

   if (optself == 0) then
     green%oper(ifreq)%ks(:,:,:,:) = czero
   else
     call add_matlu(self%hdc%matlu(:),self%oper(ifreq)%matlu(:),green%oper(ifreq)%matlu(:),natom,-1)
! do iatom = 1 , natom
   !write(6,*) 'self matlu', ifreq, self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)
   !write(6,*) 'self hdc  ', ifreq, self%hdc%matlu(1)%mat(1,1,1,1,1)
   !write(6,*) 'self_minus_hdc_oper  ', ifreq, self_minus_hdc_oper%matlu(1)%mat(1,1,1,1,1)
! enddo ! natom
     !if(paw_dmft%dmft_solv==4)  then
     !  call shift_matlu(self_minus_hdc_oper%matlu,paw_dmft%natom,cmplx(self%qmc_shift,0.d0,kind=dp),-1)
     !  call shift_matlu(self_minus_hdc_oper%matlu,paw_dmft%natom,cmplx(self%qmc_xmu,0.d0,kind=dp),-1)
     !endif
     call upfold_oper(green%oper(ifreq),paw_dmft,procb=green%distrib%procb(:),iproc=me_kpt)
   end if ! optself

   do isppol=1,nsppol
     do ikpt=1,mkmem
       do ib=1,mbandc
         green_tmp = omega_current + fermilevel - paw_dmft%eigen_dft(ib,ikpt+shift,isppol)
         if (optself == 0) then
           green%oper(ifreq)%ks(ib,ib,ikpt+shift_green,isppol) = cone / green_tmp
         else
           green%oper(ifreq)%ks(ib,ib,ikpt+shift_green,isppol) = &
             & green%oper(ifreq)%ks(ib,ib,ikpt+shift_green,isppol) + green_tmp
         end if
       end do ! ib
     end do ! ikpt
   end do ! isppol

    ! do ib1 = 1 , paw_dmft%mbandc
    !   do ib = 1 , paw_dmft%mbandc
    !     do ikpt = 1 , paw_dmft%nkpt
    !       do is = 1 , paw_dmft%nsppol
    !         if (green%procb(ifreq,ikpt)==myproc) then
!               green%oper(ifreq)%ks(is,ikpt,ib,ib1)=       &
    !           green_temp%ks(is,ikpt,ib,ib1)=       &
!&               ( omega_current     &
!&               + fermilevel                               &
!&               - paw_dmft%eigen_dft(is,ikpt,ib)) * Id(ib,ib1) &
!&               - self_minus_hdc_oper%ks(is,ikpt,ib,ib1)
          !if(ikpt==2.and.ib==ib1) then
          !  write(6,*)
          !  "self",ib1,ib,ikpt,is,ifreq,self_minus_hdc_oper%ks(is,ikpt,ib,ib1)
          !endif
!&               -
!(self%oper(ifreq)%ks(is,ikpt,ib,ib1)-self%hdc%ks(is,ikpt,ib,ib1))
!               if(prtopt>5) then
!               if(ikpt==1.and.(ifreq==1.or.ifreq==3).and.ib==16.and.ib1==16)
!               then
!                write(std_out,*) 'omega_current
!                ',omega_current
!                write(std_out,*) 'fermilevel
!                ',fermilevel
!                write(std_out,*) ' paw_dmft%eigen_dft(is,ikpt,ib)       ',
!                paw_dmft%eigen_dft(is,ikpt,ib),Id(ib,ib1)
!                write(std_out,*)
!                'self_minus_hdc_oper%ks(is,ikpt,ib,ib1)',self_minus_hdc_oper%ks(is,ikpt,ib,ib1)
!                write(std_out,*) 'green
!                ',green%oper(ifreq)%ks(is,ikpt,ib,ib1)
!               endif
!               if(ib==1.and.ib1==3) then
!                 write(std_out,*) "ff compute",ikpt,ifreq,is,ikpt,ib,ib1
!                 write(std_out,*) "ff compute",ikpt,ifreq,
!                 green_temp%ks(is,ikpt,ib,ib1)
!                 write(std_out,*) "ff details",paw_dmft%eigen_dft(is,ikpt,ib)
!                 write(std_out,*) "ff details2",fermilevel
!                 write(std_out,*) "ff details3",Id(ib,ib1)
!        !        write(std_out,*) "ff
!        details4",self_minus_hdc_oper%ks(is,ikpt,ib,ib1)
!               endif
        !     endif
        !   enddo ! is
        ! enddo ! ikpt
       !enddo ! ib
     !enddo ! ib1

!     call print_oper(green%oper(ifreq),9,paw_dmft,3)
!       write(std_out,*) 'after print_oper'
!     if(ifreq==1.or.ifreq==3) then
!         write(std_out,*) 'after print_oper', ifreq
!       write(std_out,*) 'green1  ifreq         %ks(1,1,16,16)',ifreq,green%oper(ifreq)%ks(1,1,16,16)
!     endif
!       write(std_out,*) 'before inverse_oper'
   if (optself /= 0) then
     call inverse_oper(green%oper(ifreq),1,procb=green%distrib%procb(:),iproc=me_kpt)
   end if
    !if(ifreq==1) then
    !  write(std_out,*) "1188",green_temp%ks(1,1,8,8)
    !  write(std_out,*) "1189",green_temp%ks(1,1,8,9)
    !  write(std_out,*) "1198",green_temp%ks(1,1,9,8)
    !  write(std_out,*) "1199",green_temp%ks(1,1,9,9)
    !endif
   if (optlog == 1) then
     freq2 = paw_dmft%omega_lo(ifreq)**2
     trace_tmp = czero
     do isppol=1,nsppol
       do ikpt=1,mkmem
         wtk = green%oper(ifreq)%wtk(ikpt+shift)
         if (optself == 0) then
           do ib=1,mbandc
             trace_tmp = trace_tmp + two*temp*wtk*log(green%oper(ifreq)%ks(ib,ib,ikpt+shift_green,isppol)*omega_current)
           end do ! ib
         else
           ! Use Tr(log(G(iw))) + Tr(log(G(-iw))) = Tr(log(G(iw)G(iw)^H)) so that we can use the much faster zheev instead of
           ! zgeev. Even if G(iw) and G(iw)^H have no reason to commute, this equality still holds since we only care about the trace.
           call abi_xgemm("n","c",mbandc,mbandc,mbandc,cone,green%oper(ifreq)%ks(:,:,ikpt+shift_green,isppol), &
                        & mbandc,green%oper(ifreq)%ks(:,:,ikpt+shift_green,isppol),mbandc,czero,mat_tmp(:,:),mbandc)
           call zheev("n","u",mbandc,mat_tmp(:,:),mbandc,eig(:),work(:),lwork,rwork(:),info)
           ! Do not use DOT_PRODUCT
           trace_tmp = trace_tmp + sum(log(eig(:)*freq2))*wtk*temp
         end if ! optself
       end do ! ikpt
     end do ! isppol
     if (nsppol == 1 .and. nspinor == 1) trace_tmp = trace_tmp * two
     green%trace_log = green%trace_log + dble(trace_tmp)
   end if ! optlog=1

     !if(lintegrate) then
!    accumulate integration
   if (green%w_type /= "real") then

     !do isppol=1,nsppol
     !  do ikpt=1,nkpt
     !    if (green%distrib%procb(ikpt+(isppol-1)*nkpt) /= me_spkpt) cycle

         !do ib1=1,mbandc
         !  do ib=1,mbandc
              ! if (green%procb(ifreq,ikpt)==myproc) then
          !   call
          !   add_int_fct(ifreq,green%oper(ifreq)%ks(ib,ib1,ikpt,isppol),ib==ib1,
          !   &
          !         & omega_current,2,green%occup%ks(ib,ib1,ikpt,isppol), &
          !         & paw_dmft%temp,paw_dmft%wgt_wlo(ifreq),paw_dmft%dmft_nwlo)

               !endif
         !  end do ! ib
         !end do ! ib1
     !  end do ! ikpt
     !end do ! isppol

     if (diag == 1) then
       do ib=1,mbandc
         green%occup%ks(ib,ib,1+shift:mkmem+shift,:) = green%occup%ks(ib,ib,1+shift:mkmem+shift,:) + &
                 & fac*green%oper(ifreq)%ks(ib,ib,1+shift_green:mkmem+shift_green,:)
       end do ! ib
     else
       green%occup%ks(:,:,1+shift:mkmem+shift,:) = green%occup%ks(:,:,1+shift:mkmem+shift,:) + &
          & fac*green%oper(ifreq)%ks(:,:,1+shift_green:mkmem+shift_green,:)
     end if ! diag

   end if ! green%wtype

!        write(std_out,*) 'after inverse_oper'
!      if(ikpt==1.and.is==1.and.ib==1.and.ib1==1) then
!        write(6,*) 'occup(is,ikpt,ib,ib1)',ifreq,green%occup%ks(1,1,1,1),green_temp%ks(1,1,1,1)
!      endif
!        write(std_out,*) 'green1afterinversion  %ks(1,1,16,16)',ifreq,green%oper(ifreq)%ks(1,1,16,16)
!      endif
!        write(std_out,*) 'before flush'
!      call flush(std_out)
     !endif
! ================================
! == Compute Local Green function
! ================================
   !write(message,'(2a)') ch10,' loc'
   !call wrtout(std_out,message,'COLL')
   !call flush(std_out)
   if (paw_dmft%lchipsiortho == 1 .or. optself == 1) then
     call downfold_oper(green%oper(ifreq),paw_dmft,procb=green%distrib%procb(:),iproc=me_kpt,option=option)
    !if(ifreq==1) then
    !  write(std_out,*) "4411",green_temp%matlu(1)%mat(4,4,1,1,1)
    !  write(std_out,*) "4512",green_temp%matlu(1)%mat(4,5,1,1,2)
    !  write(std_out,*) "5421",green_temp%matlu(1)%mat(5,4,1,2,1)
    !  write(std_out,*) "5522",green_temp%matlu(1)%mat(5,5,1,2,2)
    !  write(std_out,*) "(5512)",green_temp%matlu(1)%mat(5,5,1,1,2)
    !endif
!     write(std_out,*) ifreq,nproc,'before if after loc_oper'
!     if(ifreq==1.or.ifreq==11) then
!       write(std_out,*) ifreq,nproc,'before sym'
!       call print_matlu(green%oper(ifreq)%matlu,green%oper(ifreq)%natom,2,-1,0)
!       ! ok
!     endif
! call flush(std_out)
   end if ! lchipsiortho=1
     !call copy_matlu(green_temp%matlu,green%oper(ifreq)%matlu,natom)
!     if(ifreq==1.and.ifreq==11) then
!       write(std_out,*) ifreq,nproc,'after sym'
!       call print_matlu(green%oper(ifreq)%matlu,green%oper(ifreq)%natom,2,-1,0)
!       ! ok
!     endif
   if (green%oper(ifreq)%has_operks == 0) then
     ABI_FREE(green%oper(ifreq)%ks)
   end if
! call flush(std_out)
 end do ! ifreq
   ABI_NVTX_END_RANGE()
 end if


 if (optlog == 1) then

   ABI_FREE(eig)
   ABI_FREE(work)
   ABI_FREE(rwork)
   ABI_FREE(mat_tmp)

   call xmpi_sum(green%trace_log,spacecomm,ierr)

   correction = - log(two) * mbandc * nsppol * temp

   band_index = 0
   beta = one / paw_dmft%temp
   do isppol=1,nsppol
     do ikpt=1,nkpt
       nband_k = paw_dmft%nband(ikpt+(isppol-1)*nkpt)
       do ib=1,nband_k
         if (paw_dmft%band_in(ib)) cycle
         eigen = paw_dmft%eigen(ib+band_index)
         if (eigen-fermilevel >= zero) then
           correction = correction - temp*paw_dmft%wtk(ikpt)*log(one+exp(-beta*(eigen-fermilevel)))
         else
           correction = correction - temp*paw_dmft%wtk(ikpt)*(log(one+exp(beta*(eigen-fermilevel)))-beta*(eigen-fermilevel))
         end if ! eigen-fermilevel>=0
       end do ! ib
       band_index = band_index + nband_k
     end do ! ikpt
   end do ! isppol
   if (nsppol == 1 .and. nspinor == 1) correction = correction * two

   ! Do not use DOT_PRODUCT
   green%trace_log = green%trace_log + correction + fermilevel*paw_dmft%nelectval + &
      & dble(sum(green%trace_moments_log_ks(1:nmoments-1)*omega_fac(1:nmoments-1)))

 end if ! optlog=1

! =============================================
! Build total green function (sum over procs).
! =============================================
 !call xmpi_barrier(spacecomm)
!  call xmpi_sum(green%oper(ifreq)%ks,spacecomm,ierr)
!       print *, "myproc, proct, ifreq ------------------- ", myproc, ifreq
!   do ikpt=1,paw_dmft%nkpt
!     call xmpi_bcast(green%oper(ifreq)%ks(:,ikpt,:,:),procb(ifreq,ikpt),spacecomm,ierr)
!   enddo
!! or
 if (optnonxsum == 0 .and. green%oper(1)%has_operks == 1) then
   call gather_oper(green%oper(:),green%distrib,paw_dmft,opt_ksloc=1,opt_diag=diag)
 else if (optnonxsum == 0 .and. green%oper(1)%has_operks == 0) then
   message = 'optnonxsum=0 and green%oper(1)%has_operks=0: not compatible'
   ABI_BUG(message)
 end if ! optnonxsum
!   endif
! enddo ! ifreq
!  print *,"myproc", myproc
 if (optnonxsum2 == 0) then
   if (paw_dmft%lchipsiortho == 1 .or. optself == 1) then
     call gather_oper(green%oper(:),green%distrib,paw_dmft,opt_ksloc=2,opt_commkpt=1)

     if(gpu_option==ABI_GPU_DISABLED) then
       do ifreq=1,green%nw
         if (green%distrib%procf(ifreq) /= myproc) cycle
         call sym_matlu(green%oper(ifreq)%matlu(:),paw_dmft)
       end do
     else if(gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
       !FIXME: Remove those CPU-GPU transfers once gather_oper has been ported to GPU
       ! 1) Init green_oper_ndat
       ndat = green%distrib%nw_mem_kptparal(green%distrib%me_freq+1)
       do ifreq=1,green%nw
         if (green%distrib%procf(ifreq) /= myproc) cycle
         call init_oper_ndat(paw_dmft,green_oper_ndat,ndat,nkpt=green%oper(ifreq)%nkpt,opt_ksloc=2,gpu_option=gpu_option)
         if (green%oper(ifreq)%has_operks == 0) then
           green_oper_ndat%paral  = 1
           green_oper_ndat%shiftk = shift
         end if
         exit
       end do
       ! 2) Copy green%oper(:)%matlu into green_oper_ndat (CPU->GPU transfer)
       call copy_oper_to_ndat(green%oper,green_oper_ndat,ndat,green%nw,green%distrib%proct,green%distrib%me_freq,.false.)

       ! 3) Perform sym_matlu on green_oper_ndat (GPU enabled)
       call sym_matlu(green_oper_ndat%matlu(:),paw_dmft)

       ! 4) Copy back green%oper(:)%matlu from green_oper_ndat (GPU->CPU transfer)
       call copy_oper_from_ndat(green_oper_ndat,green%oper,ndat,green%nw,green%distrib%proct,&
       &    green%distrib%me_freq,.false.)
       ! 5) Destroy green_oper_ndat
       call destroy_oper(green_oper_ndat)
#endif
     end if

     call gather_oper(green%oper(:),green%distrib,paw_dmft,opt_ksloc=2)
   end if
   green%has_greenmatlu_xsum = 1
 else if (optnonxsum2 == 1) then
   green%has_greenmatlu_xsum = 0
 end if ! optnonxsum2
!     if(ifreq==1.or.ifreq==11) then
!       write(std_out,*) ifreq,nproc,'after xsum'
!       call print_matlu(green%oper(ifreq)%matlu,green%oper(ifreq)%natom,2,-1,0)
!       ! ok
!     endif
 !end do ! ifreq
 !call xmpi_barrier(spacecomm)
! write(std_out,*) 'afterxsum sym     %matlu(1)%mat(2,5,1,1,1) 1',green%oper(1)%matlu(1)%mat(2,5,1,1,1)

 if (green%w_type /= "real") then

   if (green%distrib%me_freq == 0) then

     do ib=1,mbandc
       green%occup%ks(ib,ib,1+shift:mkmem+shift,:) = green%occup%ks(ib,ib,1+shift:mkmem+shift,:) + omega_fac(1)
     end do ! ib

     do i=2,nmoments
       if (diag == 1) then
         do ib=1,mbandc
           green%occup%ks(ib,ib,1+shift:mkmem+shift,:) = green%occup%ks(ib,ib,1+shift:mkmem+shift,:) + &
                & omega_fac(i)*green%moments(i)%ks(ib,ib,:,:)
         end do ! ib
       else
         green%occup%ks(:,:,1+shift:mkmem+shift,:) = green%occup%ks(:,:,1+shift:mkmem+shift,:) + &
           & omega_fac(i)*green%moments(i)%ks(:,:,:,:)
       end if ! diag
     end do ! i

   end if ! me_freq = 0

   call gather_oper_ks(green%occup,green%distrib,paw_dmft,opt_diag=diag)

 end if ! w_type/="real"

 ABI_SFREE(omega_fac)

 if (prtopt /= 0 .and. prtopt > -100) then
   write(message,'(2a)') ch10," ===  Green's function is computed"
   call wrtout(std_out,message,'COLL')
 end if

 if (prtopt /= 0 .and. prtopt > -100 .and. (paw_dmft%lchipsiortho == 1 .or. optself == 1)) then
   write(message,'(2a)') ch10,&
     & " ===  Local Green's function has been computed and projected on local orbitals"
   call wrtout(std_out,message,'COLL')
 end if
! useless test
 if (abs(prtopt) >= 4 .and. prtopt > -100) then
   write(message,'(2a)') ch10," == Green's function is now printed for first frequency"
   call wrtout(std_out,message,'COLL')
   call print_oper(green%oper(1),9,paw_dmft,3)
   write(message,'(2a)') ch10," == Green's function is now printed for second frequency"
   call wrtout(std_out,message,'COLL')
   call print_oper(green%oper(2),9,paw_dmft,3)
   if (paw_dmft%dmft_nwlo >= 11) then
     write(message,'(2a)') ch10," == Green's function is now printed for 11th frequency"
     call wrtout(std_out,message,'COLL')
     call print_oper(green%oper(11),9,paw_dmft,3)
   end if
 end if
! call flush(std_out)

 ABI_NVTX_END_RANGE()
 call timab(624,2,tsec(:))

end subroutine compute_green
!!***

!!****f* m_green/integrate_green
!! NAME
!! integrate_green
!!
!! FUNCTION
!!  Integrate green function
!!
!! INPUTS
!!  green <type(green_type)>=green function  (green%oper(:))
!!  paw_dmft <type(m_paw_dmft)>= paw+dmft data
!!  prtopt : flag for printing
!!  opt_ksloc=   1: only integrate on the KS basis
!!               2: do the integration on the local basis
!!                 (This can work only if chipsi are renormalized!!!)
!!               3: do both calculations and test the consistency of it
!!              -1: do the integration on the KS basis, but only
!!                      compute diagonal part of the band-band density matrix
!!                      in order to compute the total charge for fermi_green
!!  opt_after_solver= to be activated after a call to impurity_solve
!!  opt_diff= check the change in the local charge compared to the previous iteration
!!  opt_fill_occnd= 0 (default) : do not fill paw_dmft%occnd with the Green's function occupations
!!                = 1 : fill paw_dmft%occnd with the Green's function occupations
!!  opt_self= 0 : assume the green%ks is diagonal
!!          = 1 (default) : assume there is a self-energy and green%ks is non-diagonal
!!
!! OUTPUT
!!   green%occup = occupations
!!
!! SOURCE

subroutine integrate_green(green,paw_dmft,prtopt,opt_ksloc,opt_after_solver,opt_diff,opt_fill_occnd,opt_self)

!Arguments ------------------------------------
 type(green_type), intent(inout) :: green
 !type(MPI_type), intent(in) :: mpi_enreg
 type(paw_dmft_type), intent(inout) :: paw_dmft
 integer, intent(in) :: prtopt
 integer, optional, intent(in) :: opt_after_solver,opt_diff,opt_fill_occnd,opt_ksloc,opt_self
!Local variables-------------------------------
 integer :: band_index,i,iatom,ib,ib1,icomp_chloc,ifreq,ikpt,isppol
 integer :: lpawu,mbandc,myproc,natom,nband_k,nkpt,nmoments,nspinor
 integer :: nsppol,optaftsolv,optdiff,optfilloccnd,option,optksloc,optself
 real(dp) :: correction,diff_chloc,fac,temp
 character(len=12) :: tag
 character(len=500) :: message
 real(dp) :: tsec(2)
 complex(dpc), allocatable :: omega_fac(:),shift(:)
 type(matlu_type), allocatable :: matlu_temp(:)
! real(dp), allocatable :: charge_loc_old(:,:)
! type(oper_type)  :: oper_c
! *********************************************************************

 DBG_ENTER("COLL")
 call timab(625,1,tsec(:))
 ABI_NVTX_START_RANGE(NVTX_DMFT_INTEGRATE_GREEN)

 if (prtopt > 0) then
   write(message,'(2a,i3,13x,a)') ch10," ===  Integrate Green's function"
   call wrtout(std_out,message,'COLL')
 end if
 if (green%w_type == "real") then
   message = 'integrate_green not implemented for real frequency'
   ABI_BUG(message)
 end if

 optfilloccnd = 0
 if (present(opt_fill_occnd)) optfilloccnd = opt_fill_occnd

 optself = 1
 if (present(opt_self)) optself = opt_self
 option = merge(3,1,optself==0) ! option for downfold_oper

! Initialize spaceComm, myproc, and master
 !spacecomm=paw_dmft%spacecomm
 myproc = paw_dmft%myproc
 !nproc=paw_dmft%nproc

! Initialise integers
 !mband   = paw_dmft%mband
 mbandc  = paw_dmft%mbandc
 natom   = paw_dmft%natom
 nkpt    = paw_dmft%nkpt
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol
 temp    = paw_dmft%temp

 nmoments = merge(green%nmoments,1,green%has_moments==1)

! Initialize green%oper before calculation (important for xmpi_sum)
! allocate(charge_loc_old(paw_dmft%natom,paw_dmft%nsppol+1))
! if(.not.present(opt_diff)) then  ! if integrate_green is called in m_dmft after calculation of self
!   charge_loc_old=green%charge_matlu
! endif
 icomp_chloc = 0
 optdiff = 0
 if (present(opt_diff)) optdiff = opt_diff

! Choose what to compute
 optksloc = 3
 if (present(opt_ksloc)) optksloc = opt_ksloc
 optaftsolv = 0
 if (present(opt_after_solver)) optaftsolv = opt_after_solver

 if (optaftsolv == 1 .and. abs(optksloc) /= 2) then
    message = "integration of ks green function should not be done after call to solver : it has not been computed"
    ABI_BUG(message)
 end if

 if (abs(optksloc) >= 2 .and. green%has_greenmatlu_xsum == 0) then
   write(message,'(4a)') ch10,&
      & "BUG: integrate_green is asked to integrate local green function",ch10,&
      & " and local green function was non broadcasted in compute_green"
    ABI_BUG(message)
 end if

! Allocations
 ABI_MALLOC(matlu_temp,(natom))
 call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),matlu_temp(:))

! =================================================
! == Integrate Local Green function ===============
 if (abs(optksloc)/2 == 1) then ! optksloc=2 or 3
! =================================================

! ==  Calculation of \int{G_{LL'}{\sigma\sigma',s}(R)(i\omega_n)}
   if (paw_dmft%lchipsiortho == 1) then
!  - Calculation of frequency sum over positive frequency
     !if (nspinor==1) option=1
     !if (nspinor==2) option=2

    !do atom=1, natom
    !   ndim=2*green%oper(1)%matlu(iatom)%lpawu+1
    !   if(green%oper(1)%matlu(iatom)%lpawu.ne.-1) then
    !         do ispinor1 = 1, nspinor
    !       do ispinor = 1, nspinor
    !     do is = 1 , nsppol
    !           do im=1,ndim
    !             do im1=1,ndim
    !               do ifreq=1,green%nw
    !                 ff(ifreq)= &
!& green%oper(ifreq)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
!                     write(std_out,*) green%omega(ifreq),ff(ifreq)," integrate
!                     green fw_lo"
   !if(present(iii).and.im==1.and.im1==1) write(std_out_default,*)
   !ch10,ifreq,ff(ifreq),"#ff"
 !                  enddo
!                   call int_fct(ff,(im==im1).and.(ispinor==ispinor1),&
!&                   option,paw_dmft,integral)
  !                 call int_fct(ff,(im==im1).and.(ispinor==ispinor1),&
!&                   2,paw_dmft,integral)  ! test_1
 !                  green%occup%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)=integral
   !if(present(iii).and.im==1.and.im1==1) write(std_out_default,*)
   !ch10,'integral',im,im1,ifreq,integral
!                   if(im==2.and.im1==5.and.is==1.and.iatom==1) then
!                     write(std_out,*) " occup
!                     %matlu(1)%mat(2,5,1,1,1)",green%occup%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
!                   endif
           !      enddo
           !    enddo
           !  enddo ! ispinor1
        !   enddo ! ispinor
        ! enddo ! is
        ! matlu_temp(iatom)%mat=green%occup%matlu(iatom)%mat
    !   endif ! lpawu=/-1
    ! enddo ! iatom

     call zero_matlu(green%occup%matlu(:),natom)

     ABI_MALLOC(omega_fac,(nmoments))

     do i=1,nmoments
       omega_fac(i) = czero
       do ifreq=green%nw,1,-1 ! NEVER change the summation order and DON'T use the intrinsic SUM
         omega_fac(i) = omega_fac(i) + paw_dmft%wgt_wlo(ifreq) / (paw_dmft%omega_lo(ifreq))**i
       end do
       omega_fac(i) = - two * temp * omega_fac(i) / (j_dpc)**i
       if (i == 1) omega_fac(i) = omega_fac(i) + half
       if (i == 2) omega_fac(i) = omega_fac(i) - cone/(four*temp)
       if (i == 4) omega_fac(i) = omega_fac(i) + cone/(dble(48)*(temp**3))
     end do ! i

     do ifreq=1,green%nw

       if (green%distrib%procf(ifreq) /= myproc) cycle
       fac = two * temp * paw_dmft%wgt_wlo(ifreq)
       do iatom=1,natom
         lpawu = paw_dmft%lpawu(iatom)
         if (lpawu == -1) cycle
         green%occup%matlu(iatom)%mat(:,:,:) = green%occup%matlu(iatom)%mat(:,:,:) + &
            & fac*green%oper(ifreq)%matlu(iatom)%mat(:,:,:)
       end do ! iatom

     end do ! ifreq

     call xmpi_matlu(green%occup%matlu(:),natom,paw_dmft%spacecomm)

     if (green%has_moments == 1) then
       do i=1,nmoments
         do iatom=1,natom
           lpawu = paw_dmft%lpawu(iatom)
           if (lpawu == -1) cycle
           green%occup%matlu(iatom)%mat(:,:,:) = green%occup%matlu(iatom)%mat(:,:,:) + &
             omega_fac(i)*green%moments(i)%matlu(iatom)%mat(:,:,:)
         end do ! iatom
       end do ! i
     else
       ABI_MALLOC(shift,(natom))
       shift(:) = omega_fac(1)
       call shift_matlu(green%occup%matlu(:),natom,shift(:),signe=-1)
       ABI_FREE(shift)
     end if

     ABI_FREE(omega_fac)

!    Print density matrix if prtopt high
     if (abs(prtopt) > 2) then
       write(message,'(2a,i10,a)') ch10,"  = green%occup%matlu from int(gloc(w))"
       call wrtout(std_out,message,'COLL')
       call print_matlu(green%occup%matlu(:),natom,prtopt=3,opt_diag=-1)
     end if ! abs(prtopt)>2

!  - Symmetrize: continue sum over k-point: Full BZ
     call sym_matlu(green%occup%matlu(:),paw_dmft)
     if (abs(prtopt) > 2) then
       write(message,'(2a,i10,a)') ch10, &
         & "  = green%occup%matlu from int(gloc(w)) with symmetrization"
       call wrtout(std_out,message,'COLL')
       call print_matlu(green%occup%matlu(:),natom,prtopt=3,opt_diag=-1)
     end if ! abs(prtopt)>2

!  - Post-treatment for summation over negative and positive frequencies:
!    necessary in the case of nspinor==2 AND nspinor==1, but valid anywhere
!    N(ll'sigmasigma')= (N(ll'sigmasigma')+ N*(l'lsigma'sigma))/2
!    because [G_{LL'}^{sigma,sigma'}(iomega_n)]*= G_{L'L}^{sigma',sigma}(-iomega_n)
     do iatom=1,natom
       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       do isppol=1,nsppol
         green%occup%matlu(iatom)%mat(:,:,isppol) = (green%occup%matlu(iatom)%mat(:,:,isppol) + &
            & conjg(transpose(green%occup%matlu(iatom)%mat(:,:,isppol)))) * half
       end do ! isppol
       matlu_temp(iatom)%mat(:,:,:) = green%occup%matlu(iatom)%mat(:,:,:)
     end do ! iatom
     if (abs(prtopt) > 2) then
       write(message,'(2a,i10,a)') ch10, &
         & "  = green%occup%matlu from int(gloc(w)) symmetrized with post-treatment"
       call wrtout(std_out,message,'COLL')
       call print_matlu(green%occup%matlu(:),natom,prtopt=3,opt_diag=-1)
     end if ! abs(prtopt)>2

     if (optaftsolv == 0) then
       call trace_oper(green%occup,green%charge_ks,green%charge_matlu(:,:),2)
       green%has_charge_matlu = 2
       green%has_charge_matlu_solver = 0
       icomp_chloc = 1
     else if (optaftsolv == 1) then ! .and. paw_dmft%dmft_solv /= 4) then
!     else if(optaftsolv==1) then
       !if(paw_dmft%dmft_solv==4) then
       !  write(message,'(a,a,a)') ch10,&
!&        "  = WARNING: Double Counting will be computed with solver charge.." ,&
!&        "might be problematic with Hirsch Fye QMC"
    !     call wrtout(std_out,message,'COLL')
    !   endif
!      This is done only when called from impurity_solver with solver
!      green function. (And if QMC is NOT used).
       if (paw_dmft%dmft_solv >= 5 .and. green%has_charge_matlu_solver /= 2) then
         write(message,'(2a,i3)') ch10,&
           & "  = BUG : has_charge_matlu_solver should be 2 and is",green%has_charge_matlu_solver
         ABI_BUG(message)
       end if
       if (paw_dmft%dmft_solv <= 4) then
         call trace_oper(green%occup,green%charge_ks,green%charge_matlu_solver(:,:),2)
         green%has_charge_matlu_solver = 2
       end if
     end if ! optaftsolv
   else
     write(message,'(a,4x,3a,4x,a)') ch10,&
       & "  Local basis is not (yet) orthonormal:",&
       & " local Green's function is thus not integrated",ch10,&
       & "  Local occupations are computed by downfolding the Kohn-Sham occupations instead"
     call wrtout(std_out,message,'COLL')
   end if ! lchipsiortho=1

 end if ! optksloc
! =================================================

! =================================================
! == Integrate Kohn Sham Green function ===========
 if (mod(abs(optksloc),2) == 1) then ! optksloc=1 or 3 or -1
!   green%occup%ks=czero ! important for xmpi_sum
!! =================================================
!! ==  Calculation of \int{G_{\nu\nu'}{k,s}(i\omega_n)}
!   ff=czero
!   do is = 1 , nsppol
!     do ikpt = 1, nkpt
!!020212       if(mod(ikpt-1,nproc)==myproc) then
!!!         print'(A,3I4)', "P ",myproc,is,ikpt
!         do ib = 1, mbandc
!           do ib1 = 1, mbandc
!             if(optksloc==-1.and.(ib/=ib1)) cycle
!             ff(:)=czero
!             do ifreq=1,green%nw
!             !!  the following line is a quick and dirty tricks and should be removed when the data will
!             !!  be correctly distributed.
!!!               if((optnonxsum==1).and.(green%procb(ifreq,ikpt)==myproc)).or.(optnonxsum==0)) then
!               if(green%procb(ifreq,ikpt)==myproc) then
!                 ff(ifreq)=green%oper(ifreq)%ks(is,ikpt,ib,ib1)
!!!                 print'(A,5I4,3(2E15.5,3x))', "P1",myproc,is,ikpt,ib,ib1,ff(:)
!               endif
!               !endif
!             enddo
!!             call int_fct(ff,ib==ib1,nspinor,paw_dmft,integral) ! here, option==1 even if nspinor==2
!!             green%occup%ks(is,ikpt,ib,ib1)=integral
!!!                 print'(A,5I4,3(2E15.5,3x))', "PP",myproc,is,ikpt,ib,ib1,ff(:)
!             call int_fct(ff,ib==ib1,2,paw_dmft,integral,green%procb(:,ikpt),myproc) ! here, option==1 even if nspinor==2
!!020212             call int_fct(ff,ib==ib1,2,paw_dmft,integral) ! here, option==1 even if nspinor==2
!             green%occup%ks(is,ikpt,ib,ib1)=integral
!!!                 print'(A,5I4,2E15.8)', "P2",myproc,is,ikpt,ib,ib1,integral
!            !   write(std_out,*) "integral",ikpt,green%occup%ks(is,ikpt,ib,ib1)
!            ! endif
!!        write(std_out,'(a,4i6,e14.5,e14.5)') "ks",is,ikpt,ib,ib1,integral
!           enddo ! ib1
!         enddo ! ib
!!020212       endif
!     enddo ! ikpt
!   enddo ! isppol


!   call xmpi_barrier(spacecomm)
!   call xmpi_sum(green%occup%ks,spacecomm,ierr)


!!   do is = 1 , nsppol
!!     do ikpt = 1, nkpt
!!         do ib = 1, mbandc
!!           do ib1 = 1, mbandc
!!             write(6,'(A,5I4,2E15.8)') "AAAFTERXSUM",myproc,is,ikpt,ib,ib1,green%occup%ks(is,ikpt,ib,ib1)
!!             print '(A,5I4,2E15.8)', "AFTERXSUM P",myproc,is,ikpt,ib,ib1,green%occup%ks(is,ikpt,ib,ib1)
!!           enddo ! ib1
!!         enddo ! ib
!!     enddo ! ikpt
!!   enddo ! isppol
               !write(std_out,*) "occup%ks ik1",green%occup%ks(1,1,1,3)
               !write(std_out,*) "occup%ks ik2",green%occup%ks(1,2,1,3)
!  - Post-treatment for summation over negative and positive frequencies:
!    necessary in the case of nspinor==2, but valid everywhere
!    N(k,n_1,n_2)= (N(k,n_1,n_2)+ N*(k,n_2,n_1))/2
!    because [G_{k}^{n_1,n_2}(iomega_n)]*= G_{k}^{n_2,n_1}(-iomega_n)
   do isppol=1,nsppol
     do ikpt=1,nkpt
       green%occup%ks(:,:,ikpt,isppol) = (green%occup%ks(:,:,ikpt,isppol)+ &
          & conjg(transpose(green%occup%ks(:,:,ikpt,isppol)))) * half
     end do ! ikpt
   end do ! isppol
               !write(std_out,*) "occup%ks ik1 BB",green%occup%ks(1,1,1,3)
               !write(std_out,*) "occup%ks ik2 AA",green%occup%ks(1,2,1,3)
   if (optfilloccnd == 1) then
     band_index = 0
     fac = merge(two,one,nsppol==1.and.nspinor==1)
     do isppol=1,nsppol
       do ikpt=1,nkpt
         do ib1=1,mbandc
           do ib=1,mbandc
             paw_dmft%occnd(1,paw_dmft%include_bands(ib),paw_dmft%include_bands(ib1),&
                & ikpt,isppol) = fac * dble(green%occup%ks(ib,ib1,ikpt,isppol))
             paw_dmft%occnd(2,paw_dmft%include_bands(ib),paw_dmft%include_bands(ib1),&
                & ikpt,isppol) = fac * aimag(green%occup%ks(ib,ib1,ikpt,isppol))
           end do  ! ib
         end do ! ib1
         if (paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7) then
           nband_k = paw_dmft%nband(ikpt+(isppol-1)*nkpt)
           do ib=1,nband_k
             if (paw_dmft%band_in(ib)) cycle
             paw_dmft%occnd(:,:,ib,ikpt,isppol)  = zero
             paw_dmft%occnd(1,ib,ib,ikpt,isppol) = fac * &
               & occup_fd(paw_dmft%eigen(ib+band_index),paw_dmft%fermie,paw_dmft%temp)
           end do ! ib
           band_index = band_index + nband_k
         end if ! use_all_bands
       end do ! ikpt
     end do ! isppol
   end if ! optfilloccnd

   if (optksloc > 0) then
!  - Compute local occupations
     call downfold_oper(green%occup,paw_dmft,procb=green%distrib%procb(:), &
                      & iproc=green%distrib%me_kpt,option=option)
     call xmpi_matlu(green%occup%matlu(:),natom,green%distrib%comm_kpt)
     if (abs(prtopt) > 2) then
       write(message,'(2a,i10,a)') ch10,&
          & "  = green%occup%matlu from projection of int(gks(w)) without symmetrization"
       call wrtout(std_out,message,'COLL')
       call print_matlu(green%occup%matlu(:),natom,prtopt=3,opt_diag=-1)
     end if ! abs(prtopt)>2

!  - Symmetrize: continue sum over k-point: Full BZ
     call sym_matlu(green%occup%matlu(:),paw_dmft)
     if (abs(prtopt) >= 2) then
!       write(message,'(a,a,i10,a)') ch10,&
!&        "  = green%occup%matlu from projection of int(gks(w)) with symetrization"
       write(message,'(2a,i10,a)') ch10,"  ==  Occupation matrix from downfolded Kohn-Sham occupations"
       call wrtout(std_out,message,'COLL')
       call print_matlu(green%occup%matlu(:),natom,prtopt=3,opt_diag=0)
     end if ! abs(prtopt)>=2

!    - If the trace of occup matrix in the LOCAL basis was not done
!    before because lchipsiortho/=1 , do it now
     if (paw_dmft%lchipsiortho /= 1 .and. abs(prtopt) > 0) then
       call trace_oper(green%occup,green%charge_ks,green%charge_matlu(:,:),2)
       green%has_charge_matlu = 2
       icomp_chloc = 1
     end if ! lchipsiortho/=1 and abs(prtopt)>0
   end if ! optksloc>0: only diagonal elements of G\nunu' are computed

!  - Compute trace over ks density matrix
   call trace_oper(green%occup,green%charge_ks,green%charge_matlu(:,:),1)
   if (paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7) then

     band_index = 0
     correction = zero
     do isppol=1,nsppol
       do ikpt=1,nkpt
         nband_k = paw_dmft%nband(ikpt+(isppol-1)*nkpt)
         do ib=1,nband_k
           if (paw_dmft%band_in(ib)) cycle
           correction = correction + &
             & occup_fd(paw_dmft%eigen(ib+band_index),paw_dmft%fermie,paw_dmft%temp)*green%occup%wtk(ikpt)
         end do ! ib
         band_index = band_index + nband_k
       end do ! ikpt
     end do ! isppol
     if (nsppol == 1 .and. nspinor == 1) correction = correction * two
     green%charge_ks = green%charge_ks + correction
   end if ! use_all_bands
   if (abs(prtopt) > 0) then
     write(tag,'(f12.6)') green%charge_ks
     write(message,'(3a)') ch10,&
       & "  ==  Total number of electrons from integration of Kohn-Sham Green's function is : ",adjustl(tag)
     call wrtout(std_out,message,'COLL')
     write(tag,'(f12.6)') paw_dmft%nelectval
     write(message,'(8x,3a)') " (should be ",trim(adjustl(tag)),")"
     call wrtout(std_out,message,'COLL')
   end if ! abs(prtopt)>0
 end if ! optksloc
! =================================================

! =================================================
! Tests and compute precision on local charge
! =================================================
!  - Check that if, renormalized psichi are used, occupations matrices
!    obtained directly from local green function or, through kohn sham
!    occupations are the same.
 if ((abs(optksloc) == 3) .and. (paw_dmft%lchipsiortho == 1)) then ! optksloc= 3
   call diff_matlu("Local projection of Kohn-Sham occupations ",&
        & "Integration of local Green's function ",&
        & green%occup%matlu(:),matlu_temp(:),natom,1,tol4)
   write(message,'(2a)') ch10,&
       & "  ***** => Calculations of Green's function in Kohn-Sham and local spaces are coherent ****"
   call wrtout(std_out,message,'COLL')
 end if ! abs(optksloc)=3

! == Precision on charge_matlu (done only if local charge was computed ie not for optksloc=-1)
 if (icomp_chloc == 1 .and. paw_dmft%idmftloop >= 1 .and. optdiff == 1) then ! if the computation was done here.
   if (green%has_charge_matlu_prev == 2) then
     diff_chloc = zero
     do iatom=1,natom
       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       do isppol=1,nsppol
         diff_chloc = diff_chloc + &
           & (green%charge_matlu_prev(isppol,iatom)-green%charge_matlu(isppol,iatom))**2
       end do ! isppol
     end do ! iatom
     if (sqrt(diff_chloc) < paw_dmft%dmft_lcpr) then
       green%ichargeloc_cv = 1
       write(message,'(a,8x,a,e9.2,a,8x,a)') ch10,"Change of correlated number of electron =<",&
          & paw_dmft%dmft_lcpr,ch10,"DMFT Loop: Local charge is converged"
       call wrtout(std_out,message,'COLL')
     else
       green%ichargeloc_cv = 0
       write(message,'(a,8x,a,e9.2,a,8x,a)') ch10,"Change of correlated number of electron  >",&
         & paw_dmft%dmft_lcpr,ch10,"DMFT Loop: Local charge is not converged"
       call wrtout(std_out,message,'COLL')
     end if ! sqrt(diff_chloc)<dmft_lcpr
   end if ! green%has_charge_matlu_prev=2
   green%charge_matlu_prev(:,:) = green%charge_matlu(:,:)
   green%has_charge_matlu_prev = 2
 end if ! icomp_chloc=1 and present(opt_diff) and idmftloop>=1

 call destroy_matlu(matlu_temp(:),natom)
 ABI_FREE(matlu_temp)
! deallocate(charge_loc_old)
 ABI_NVTX_END_RANGE()
 call timab(625,2,tsec(:))
 DBG_EXIT("COLL")

end subroutine integrate_green
!!***

!!****f* m_green/icip_green
!! NAME
!!  icip_green
!!
!! FUNCTION
!!  init, compute, integrate and print lda green function
!!
!! INPUTS
!!  char1 = character which precises the type of green function computed.
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  green  <type(green_type)>= green function data
!!  pawprtvol  = option for printing
!!  self <type(self_type)>= variables related to self-energy
!!  opt_self = optional argument, if =1, upfold self-energy
!!  opt_moments = 1 to activate the computation of high-frequency moments
!!  opt_log = 1 to compute Tr(log(G))
!!
!! OUTPUT
!!
!! SOURCE

subroutine icip_green(char1,green,paw_dmft,pawprtvol,self,opt_self,opt_moments,opt_log)

!Arguments ------------------------------------
 !type(MPI_type), intent(in) :: mpi_enreg
 type(green_type), intent(inout) :: green
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(self_type), intent(inout) :: self
 integer, intent(in) :: pawprtvol
 character(len=*), intent(in) :: char1
 integer, optional, intent(in) :: opt_log,opt_moments,opt_self
!Local variables ------------------------------------
 integer :: option,optlog,optmoments,optself,optlocks,prtopt_for_integrate_green,opt_nonxsum_icip
 character(len=500) :: message
! *********************************************************************

 !green%whichgreen="DFT"
 prtopt_for_integrate_green = 2

 optself = 0
 if (present(opt_self)) optself = opt_self
 opt_nonxsum_icip = 1
 !if(paw_dmft%dmftcheck==1) then ! for fourier_green
 !  opt_nonxsum_icip=0
 !endif

 optlog = 0
 if (present(opt_log)) optlog = opt_log

 optmoments = 0
 if (present(opt_moments)) optmoments = opt_moments

 call init_green(green,paw_dmft,opt_moments=optmoments)

! useless test ok
! call printocc_green(green,1,paw_dmft,2)
! write(std_out,*)" printocc_green zero finished "

!== Compute green%oper(:)%ks
!== Deduce  green%oper(:)%matlu(:)%mat
 call compute_green(green,paw_dmft,pawprtvol,self,opt_self=optself,opt_nonxsum=opt_nonxsum_icip,opt_log=optlog)
 if (paw_dmft%dmft_prgn >= 1 .and. paw_dmft%dmft_prgn <= 2) then
   optlocks = paw_dmft%dmft_prgn*2 + 1 ! if dmft_prgn==2 => do not print
   if (paw_dmft%lchipsiortho == 1 .and. pawprtvol > -100) then
      call print_green(char1,green,optlocks,paw_dmft)
!     call print_green('inicip',green,1,paw_dmft,pawprtvol=1,opt_wt=1)
   end if
 end if

!== Integrate green%oper(:)%ks
!== Integrate green%oper(:)%matlu(:)%mat
 call integrate_green(green,paw_dmft,prtopt_for_integrate_green,opt_ksloc=3,opt_self=optself)!,opt_nonxsum=opt_nonxsum_icip)
!== Print green%oper(:)%ks
!== Print green%oper(:)%matlu(:)%mat
 if (char1 == "DFT") then
   option = 1
   if (self%oper(1)%matlu(1)%lpawu /= -1) then
     if (abs(dble(self%oper(1)%matlu(1)%mat(1,1,1))) > tol7) then
! todo_ab: generalise this
       write(message,'(2a)') ch10,&
          & "Warning: a DFT calculation is carried out and self is not zero"
       ABI_WARNING(message)
!       call abi_abort('COLL')
     end if
   end if
 else
   option = 5
 end if ! char1

 if (pawprtvol > -100) then
   call printocc_green(green,option,paw_dmft,3,chtype=char1)
 end if

end subroutine icip_green
!!***

!!****f* m_green/fourier_green
!! NAME
!! fourier_green
!!
!! FUNCTION
!!  integrate green function in the band index basis
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>= crystal structure data.
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! OUTPUT
!!
!! SOURCE

subroutine fourier_green(cryst_struc,green,paw_dmft,opt_ksloc,opt_tw)

!Arguments ------------------------------------
!type
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type),intent(inout) :: green
 !type(MPI_type), intent(in) :: mpi_enreg
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer,intent(in) :: opt_ksloc,opt_tw ! fourier on ks or local

!local variables-------------------------------
 integer :: iatom,ib,ib1,ierr,ifreq,ikpt,im,im1,iparal,is,ispinor,ispinor1,itau
 integer :: mband,mbandc,myproc,natom,ndim,nkpt,nproc,nspinor,nsppol,spacecomm!,opt_four
 character(len=500) :: message
! complex(dpc):: ybcbeg,ybcend
! arrays
 complex(dpc), allocatable :: fw(:)
 complex(dpc), allocatable :: ft(:)
 type(green_type) :: green_temp
! *********************************************************************
! ybcbeg=czero
! ybcend=czero

! define spaceComm, myproc, and nproc from world communicator
! and mpi_enreg
 spacecomm=paw_dmft%spacecomm
 myproc=paw_dmft%myproc
 nproc=paw_dmft%nproc

! Initialise integers
 mband   = paw_dmft%mband
 mbandc  = paw_dmft%mbandc
 nkpt    = paw_dmft%nkpt
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol
 natom   = cryst_struc%natom

! Only imaginary frequencies here
 if(green%w_type=="real") then
   message = 'fourier_green not implemented for real frequency'
   ABI_BUG(message)
 endif

! Initialise temporary green function
 call init_green(green_temp,paw_dmft)
 call init_green_tau(green_temp,paw_dmft)

 !green%oper(:)%matlu(1)%mat(1,1,1,1,1)
 ABI_MALLOC(fw,(green%nw))
 ABI_MALLOC(ft,(green%dmftqmc_l))

!==============================================
! == Inverse fourier transformation ===========
!==============================================

 if(opt_tw==-1) then

!  = For Kohn Sham green function
!==============================================
   if(opt_ksloc==1.and.green%use_oper_tau_ks==1) then

     do is = 1 , nsppol
       do ikpt = 1, nkpt
         do ib = 1, mbandc
           do ib1 = 1, mbandc
             do ifreq=1,green%nw
               fw(ifreq)=green%oper(ifreq)%ks(ib,ib1,ikpt,is)
             enddo
             call fourier_fct(fw,ft,ib==ib1,green%dmftqmc_l,-1,paw_dmft) ! inverse fourier
             do itau=1,green%dmftqmc_l
               green_temp%oper_tau(itau)%ks(ib,ib1,ikpt,is)=ft(itau)
             enddo
             if(ib==ib1) then
               green%occup_tau%ks(ib,ib1,ikpt,is)=ft(1)+one
             else
               green%occup_tau%ks(ib,ib1,ikpt,is)=ft(1)
             endif
           enddo ! ib1
         enddo ! ib
       enddo ! ikpt
     enddo ! isppol

!  = Post-treatment: necessary in the case of nspinor==2, but valid anywhere
!    because G(tau)=G_{nu,nu'}(tau)+[G_{nu,nu'}(tau)]*

     do is = 1 , nsppol
       do ikpt = 1, nkpt
         do ib = 1, mbandc
           do ib1 = 1, mbandc
             do itau=1,green%dmftqmc_l
               green%oper_tau(itau)%ks(ib,ib1,ikpt,is)=&
&                (green_temp%oper_tau(itau)%ks(ib,ib1,ikpt,is)+ &
&                 conjg(green_temp%oper_tau(itau)%ks(ib1,ib,ikpt,is)))/two
               if(ib==ib1) then
                 green%occup_tau%ks(ib,ib1,ikpt,is)=green%oper_tau(1)%ks(ib,ib1,ikpt,is)+one
               else
                 green%occup_tau%ks(ib,ib1,ikpt,is)=green%oper_tau(1)%ks(ib,ib1,ikpt,is)
               endif
             enddo
           enddo ! ib1
         enddo ! ib
       enddo ! ikpt
     enddo ! isppol
     call downfold_oper(green%occup_tau,paw_dmft)
     call sym_matlu(green%occup_tau%matlu,paw_dmft)
     write(message,'(a,a,i10,a)') ch10,"  green%occup_tau%matlu from green_occup_tau%ks"
     call wrtout(std_out,message,'COLL')
     call print_matlu(green%occup_tau%matlu,natom,prtopt=3)
   endif

!  = For local green function
!==============================================
   if(opt_ksloc ==2) then

     iparal=0
     do iatom=1, natom
       if(green%oper(1)%matlu(iatom)%lpawu.ne.-1) then
         ndim=2*green%oper(1)%matlu(iatom)%lpawu+1
         do is = 1 , nsppol
           do ispinor = 1, nspinor
             do ispinor1 = 1, nspinor
               do im=1,ndim
                 do im1=1,ndim
                   iparal=iparal+1
                   if(mod(iparal-1,nproc)==myproc) then
                     do ifreq=1,green%nw
                       fw(ifreq)=green%oper(ifreq)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,is)
                     enddo
                     ! inverse fourier
!                     write(std_out,'(a)') "fourierbeforeposttreatement,ispinor,ispinor1,is,im,im1"
!                     write(std_out,'(a,5i4,f12.5,f12.5)') "fourier",ispinor,ispinor1,is,im,im1
!                     write(std_out,'(a,e12.5,e12.5)')&
!                     &"green%oper(4)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)"&
!&                     ,green%oper(4)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
                     call fourier_fct(fw,ft,(im==im1).and.(ispinor==ispinor1),green%dmftqmc_l,-1,paw_dmft)
                     do itau=1,green%dmftqmc_l
                       green_temp%oper_tau(itau)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,is)=ft(itau)
!                     write(std_out,*) itau,green_temp%oper_tau(itau)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
                     enddo
!                    if((im==im1).and.(ispinor==ispinor1)) then
!                      green%occup_tau%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)=ft(1)+one
!                    else
!                      green%occup_tau%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)=ft(1)
!                    endif
                   endif ! iparal
                 enddo
               enddo
             enddo ! ispinor1
           enddo ! ispinor
         enddo ! isppol
       endif ! lpawu.ne.-1
     enddo ! iatom

!    Parallelisation must be finished here, because in the post
!    treatment, value from different proc will be mixed.
     call xmpi_barrier(spacecomm)
     do iatom=1,natom
       if (green%oper(1)%matlu(iatom)%lpawu==-1) cycle
       do itau=1,green%dmftqmc_l
        call xmpi_sum(green_temp%oper_tau(itau)%matlu(iatom)%mat,spacecomm,ierr)
       enddo
     enddo
     call xmpi_barrier(spacecomm)

!  = Post-treatment: necessary in the case of nspinor==2, but valid anywhere
!    because G(tau)=G_{LL'}^{sigma,sigma'}(tau)+[G_{L'L}^{sigma',sigma}(tau)]*

     if(nspinor>0) Then
       do iatom=1, natom
         if(green%oper(1)%matlu(iatom)%lpawu.ne.-1) then
           ndim=2*green%oper(1)%matlu(iatom)%lpawu+1
           do is = 1 , nsppol
             do ispinor = 1, nspinor
               do ispinor1 = 1, nspinor
                 do im=1,ndim
                   do im1=1,ndim
!                   write(std_out,'(a,5i4,f12.5,f12.5)') "fourier -1",ispinor,ispinor1,is,im,im1
                     do itau=1,green%dmftqmc_l
                       green%oper_tau(itau)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,is)=&
&                      ((green_temp%oper_tau(itau)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,is)+ &
&                      conjg(green_temp%oper_tau(itau)%matlu(iatom)%mat(im1+(ispinor1-1)*ndim,im+(ispinor-1)*ndim,is))))/two
!                       write(std_out,*) itau,green%oper_tau(itau)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
                       if((im==im1).and.(ispinor==ispinor1)) then
                         green%occup_tau%matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,is)=&
&                         green%oper_tau(1)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,is)+one
                       else
                         green%occup_tau%matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,is)=&
&                         green%oper_tau(1)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,is)
                       endif ! test diag
                     enddo  ! itau
                   enddo  ! im1
                 enddo ! im
               enddo ! ispinor1
             enddo ! ispinor
           enddo ! isppol
         endif
       enddo ! iatom
     endif ! nspinor>0
   endif ! opt_ksloc=2
 endif ! opt_tw==-1


!==============================================
! == Direct fourier transformation
!==============================================
! todo_ba ft useful only for diagonal elements ...

 if(opt_tw==1) then

!  = For local green function
!==============================================
   if(opt_ksloc ==2) then

     iparal=0
     do iatom=1, natom
!  put to zero (usefull at least for parallelism)
       if(green%oper(1)%matlu(iatom)%lpawu.ne.-1) then
         do ifreq=1,green%nw
           green%oper(ifreq)%matlu(iatom)%mat=czero
         enddo
         ndim=2*green%oper(1)%matlu(iatom)%lpawu+1
         do is = 1 , nsppol
           do ispinor = 1, nspinor
             do ispinor1 = 1, nspinor
               do im=1,ndim
                 do im1=1,ndim
                   iparal=iparal+1
                   if(mod(iparal-1,nproc)==myproc) then
                     do itau=1,green%dmftqmc_l
                       ft(itau)=green%oper_tau(itau)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,is)
                     enddo
!                     if(im==1.and.im1==1) write(std_out,*) "ft(itau=1)",is,ft(1) ! debug
                     call fourier_fct(fw,ft,(im==im1).and.(ispinor==ispinor1),green%dmftqmc_l,1,paw_dmft)
                     do ifreq=1,green%nw
                       green%oper(ifreq)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,is)=fw(ifreq)
                     enddo
!                     if(im==1.and.im1==1) write(std_out,*) "fw(ifreq=1)",is,fw(1) ! debug
                   endif
                 enddo
               enddo
             enddo ! ispinor1
           enddo ! ispinor
         enddo ! isppol
       endif ! lpawu=/-1
     enddo ! iatom
     call xmpi_barrier(spacecomm)
     do iatom=1,natom
       if (green%oper(1)%matlu(iatom)%lpawu==-1) cycle
       do ifreq=1,green%nw
       call xmpi_sum(green%oper(ifreq)%matlu(iatom)%mat,spacecomm,ierr)
       enddo
     enddo
   endif ! opt_ksloc=2
 endif ! opt_tw==-1
 ABI_FREE(fw)
 ABI_FREE(ft)
 call destroy_green_tau(green_temp)
 call destroy_green(green_temp)

end subroutine fourier_green
!!***


!!****f* m_green/check_fourier_green
!! NAME
!! check_fourier_green
!!
!! FUNCTION
!!  Check fourier transformations
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>= crystal structure data.
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawang <type(pawang)>=paw angular mesh and related data
!!
!! OUTPUT
!!
!! SOURCE

subroutine check_fourier_green(cryst_struc,green,paw_dmft)

!Arguments ------------------------------------
!type
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type),intent(inout) :: green
 !type(MPI_type), intent(in) :: mpi_enreg
 type(paw_dmft_type), intent(inout) :: paw_dmft

!local variables-------------------------------
 type(green_type) :: green_check
 character(len=500) :: message
! *********************************************************************
! Only imaginary frequencies here
 if(green%w_type=="real") then
   message = 'check_fourier_green not implemented for real frequency'
   ABI_BUG(message)
 endif

 call init_green(green_check,paw_dmft)
 call init_green_tau(green_check,paw_dmft)
 call copy_green(green,green_check,opt_tw=2)

 write(message,'(2a,i3,13x,a)') ch10,'   ===  Inverse Fourier Transform w->t of Weiss Field'
 call wrtout(std_out,message,'COLL')
 call fourier_green(cryst_struc,green_check,paw_dmft,&
& opt_ksloc=2,opt_tw=-1)

 write(message,'(3a)') ch10,' === Print (for check by user) of occupation matrix'&
&   ,' after  fourier transform with respect to initial one'
 call wrtout(std_out,message,'COLL')
 call printocc_green(green_check,6,paw_dmft,3)

 write(message,'(2a,i3,13x,a)') ch10,'   ===  Direct Fourier Transform t->w of Weiss Field'
 call wrtout(std_out,message,'COLL')
 call fourier_green(cryst_struc,green_check,paw_dmft,&
& opt_ksloc=2,opt_tw=1)
! call print_matlu(green%oper(1)%matlu,paw_dmft%natom,1) ! debug

 call integrate_green(green_check,paw_dmft,&
& prtopt=2,opt_ksloc=2)

 write(message,'(3a)') ch10,' === Print (for check by user) of occupation matrix'&
&   ,' after double fourier transform with respect to initial one'
 call wrtout(std_out,message,'COLL')
 call printocc_green(green_check,5,paw_dmft,3)

 call destroy_green_tau(green_check)
 call destroy_green(green_check)

end subroutine check_fourier_green
!!***

!!****f* m_green/compa_occup_ks
!! NAME
!! compa_occup_ks
!!
!! FUNCTION
!!  Compare occupations to Fermi Dirac Occupations
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! OUTPUT
!!
!! SOURCE

subroutine compa_occup_ks(green,paw_dmft)

!Arguments ------------------------------------
!type
 type(green_type),intent(inout) :: green
 type(paw_dmft_type), intent(inout) :: paw_dmft

!local variables-------------------------------
 character(len=500) :: message
 integer :: ib,ikpt,isppol
 integer :: ib_,ikpt_,isppol_
 integer :: ib__,ikpt__,isppol__
 real(dp) :: diffmax,occ1,occ2,occ_a,occ_b,diffrel,occ_aa,occ_bb
! *********************************************************************
 diffmax=zero
 diffrel=zero
 do isppol=1,paw_dmft%nsppol
   do ikpt=1,paw_dmft%nkpt
     do ib=1,paw_dmft%mbandc
       occ1=occup_fd(paw_dmft%eigen_dft(ib,ikpt,isppol),paw_dmft%fermie,paw_dmft%temp)
       occ2=real(green%occup%ks(ib,ib,ikpt,isppol))
       if(abs(occ1-occ2)>diffmax) then
         diffmax=abs(occ1-occ2)
         occ_a=occ1
         occ_b=occ2
         ib_=ib;isppol_=isppol;ikpt_=ikpt
       endif
       if(abs(two*(occ1-occ2)/(occ1+occ2))>diffrel) then
         diffrel=abs(two*(occ1-occ2)/(occ1+occ2))
         occ_aa=occ1
         occ_bb=occ2
         ib__=ib;isppol__=isppol;ikpt__=ikpt
       endif
     enddo
   enddo
 enddo


 write(message,'(2a)') ch10,'   ===  Compare green function occupations and Fermi Dirac occupations'
 call wrtout(std_out,message,'COLL')
 write(message,'(2a,f12.5,2a,f12.5,2a,f12.5,2a,3i5)') ch10,'     =  Max difference is',diffmax,&
&                           ch10,'        Corresponding Occupation from green function is',occ_b,&
&                           ch10,'        Corresponding Occupation Fermi Dirac weight  is',occ_a,&
&                           ch10,'        (For polarization, k-point and band index)     ',isppol_,ikpt_,ib_
 call wrtout(std_out,message,'COLL')
 write(message,'(2a,f12.5,2a,f12.5,2a,f12.5,2a,3i5)') ch10,'     =  Max relative difference is',diffrel,&
&                           ch10,'        Corresponding Occupation from green function is',occ_bb,&
&                           ch10,'        Corresponding Occupation Fermi Dirac weight  is',occ_aa,&
&                           ch10,'        (For polarization, k-point and band index)     ',isppol__,ikpt__,ib__
 call wrtout(std_out,message,'COLL')

end subroutine compa_occup_ks
!!***

!!****f* m_green/add_int_fct
!! NAME
!! add_int_fct
!!
!! FUNCTION
!! Do integration in matsubara space
!!
!! COPYRIGHT
!! Copyright (C) 2006-2025 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ff= function is frequency space
!!  ldiag    = option according to diagonal or non-diagonal elements
!!  option = nspinor
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  proct= for parallelism
!!
!! SIDE EFFECTS
!!  * integral = integral of ff over matsubara frequencies (there is an accumulation in the present routine, so intent(inout))
!!  * ft= function is time space
!!
!! SOURCE
subroutine add_int_fct(ifreq,ff,ldiag,omega_current,option,integral,temp,wgt_wlo,dmft_nwlo)

!Arguments ------------------------------------
!type
 integer,intent(in) :: ifreq
 logical,intent(in) :: ldiag
 integer,intent(in) :: option,dmft_nwlo
 complex(dpc),intent(inout) :: integral
 complex(dpc), intent(in) :: ff
 complex(dpc), intent(in) :: omega_current
 real(dp), intent(in) :: temp, wgt_wlo

!local variables-------------------------------
 real(dp)  :: omega
! *********************************************************************
 omega=aimag(omega_current)
 if(ldiag) then

  if(option==1) then ! nspinor==1
!   write(500,*) paw_dmft%omega_lo(ifreq),real(ff(ifreq)),imag(ff(ifreq))
    integral=integral+2.d0*temp *                         &
&      real( ff-one / ( j_dpc*omega ) ) *  &
&      wgt_wlo
    if(ifreq==dmft_nwlo) integral=integral+half
!    integral=integral+half
    ! the if is here, to count only one time this correction
  endif

  if(option==2) then ! nspinor==2
    integral=integral+2.d0*temp *                         &
&       ( ff-one / ( j_dpc*omega ) ) *  &
&   wgt_wlo
    if(ifreq==dmft_nwlo) integral=integral+half
!    integral=integral+half
    ! the if is here, to count only one time this correction
  endif


 else   ! ldiag

! write(std_out,*) "nondiag"
  if(option==1) then
    integral=integral+2.d0*temp *   &
&    real( ff ) *                     &
&   wgt_wlo
  endif
  if(option==2) then
    integral=integral+2.d0*temp *   &
&        ff   *                     &
&    wgt_wlo
  endif
 endif  ! ldiag

end subroutine add_int_fct
!!***

!!****m* m_green/int_fct
!! NAME
!! int_fct
!!
!! FUNCTION
!! Do integration in matsubara space
!!
!! COPYRIGHT
!! Copyright (C) 2006-2025 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ff= function is frequency space
!!  ldiag    = option according to diagonal or non-diagonal elements
!!  option = nspinor
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  proct= for parallelism
!!
!! OUTPUT
!!  integral = integral of ff over matsubara frequencies
!!
!! SIDE EFFECTS
!!  ft= function is time space
!!
!! SOURCE
subroutine int_fct(ff,ldiag,option,paw_dmft,integral,procb,myproc)

!Arguments ------------------------------------
!type
 logical,intent(in) :: ldiag
 integer,intent(in) :: option
 complex(dpc),intent(out) :: integral
 type(paw_dmft_type), intent(in) :: paw_dmft
 complex(dpc), intent(in) :: ff(paw_dmft%dmft_nwlo)
 integer, optional, intent(in) :: procb(paw_dmft%dmft_nwlo)
 integer, optional, intent(in) :: myproc

!local variables-------------------------------
 logical, allocatable :: procb2(:)
 character(len=500) :: message
 integer :: ifreq
! *********************************************************************
 ABI_MALLOC(procb2,(paw_dmft%dmft_nwlo))
 if(present(procb).and.present(myproc)) then
   do ifreq=1,paw_dmft%dmft_nwlo
    procb2(ifreq)=(procb(ifreq)==myproc)
   enddo
 else if(present(procb).and..not.present(myproc)) then
   write(message,'(a,a,2(e15.4))') ch10,&
&    "BUG: procb is present and not myproc in int_fct"
   ABI_BUG(message)
 else if(.not.present(procb).and.present(myproc)) then
   write(message,'(a,a,2(e15.4))') ch10,&
&    "BUG: procb is not present and myproc is in int_fct"
   ABI_BUG(message)
 else
   do ifreq=1,paw_dmft%dmft_nwlo
     procb2(ifreq)=(1==1)
   enddo
 endif

 integral=czero
 if(ldiag) then

  if(option==1) then ! nspinor==1
    do ifreq=1,paw_dmft%dmft_nwlo
      if(procb2(ifreq)) then
!     write(500,*) paw_dmft%omega_lo(ifreq),real(ff(ifreq)),imag(ff(ifreq))
       integral=integral+2.d0*paw_dmft%temp *                         &
&        real( ff(ifreq)-one / ( j_dpc*paw_dmft%omega_lo(ifreq) ) ) *  &
&        paw_dmft%wgt_wlo(ifreq)
      endif
    enddo
    if(procb2(paw_dmft%dmft_nwlo)) integral=integral+half
!    integral=integral+half
    ! the if is here, to count only one time this correction
  endif

  if(option==2) then ! nspinor==2
    do ifreq=1,paw_dmft%dmft_nwlo
      if(procb2(ifreq)) then
        integral=integral+2.d0*paw_dmft%temp *                         &
&           ( ff(ifreq)-one / ( j_dpc*paw_dmft%omega_lo(ifreq) ) ) *  &
&       paw_dmft%wgt_wlo(ifreq)
      endif
    enddo
    if(procb2(paw_dmft%dmft_nwlo)) integral=integral+half
!    integral=integral+half
    ! the if is here, to count only one time this correction
  endif


 else   ! ldiag

! write(std_out,*) "nondiag"
  if(option==1) then
    do ifreq=1,paw_dmft%dmft_nwlo
      if(procb2(ifreq)) then
        integral=integral+2.d0*paw_dmft%temp *   &
&        real( ff(ifreq) ) *                     &
&       paw_dmft%wgt_wlo(ifreq)
      endif
    enddo
  endif
  if(option==2) then
    do ifreq=1,paw_dmft%dmft_nwlo
      if(procb2(ifreq)) then
        integral=integral+2.d0*paw_dmft%temp *   &
&            ff(ifreq)   *                     &
&        paw_dmft%wgt_wlo(ifreq)
      endif
    enddo
  endif
 endif  ! ldiag
 ABI_FREE(procb2)

end subroutine int_fct
!!***

!!****f* m_green/fourier_fct
!! NAME
!! fourier_fct
!!
!! FUNCTION
!! Do fourier transformation from matsubara space to imaginary time
!! (A spline is performed )
!!
!! COPYRIGHT
!! Copyright (C) 2006-2025 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ldiag    = option according to diagonal or non-diagonal elements
!!  opt_four = option for direct (1) or inverse (-1) fourier transform
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  fw= function is frequency space
!!  ft= function is time space
!!
!! SOURCE
subroutine fourier_fct(fw,ft,ldiag,ltau,opt_four,paw_dmft)

!Arguments ------------------------------------
!type
 logical,intent(in) :: ldiag
 integer,intent(in) :: ltau,opt_four
 type(paw_dmft_type), intent(in) :: paw_dmft
 complex(dpc), intent(inout) :: fw(paw_dmft%dmft_nwlo)
 complex(dpc), intent(inout) :: ft(ltau)

!local variables-------------------------------
 complex(dpc), allocatable ::  splined_li(:)
 complex(dpc), allocatable ::  tospline_li(:)
! complex(dpc), allocatable :: fw1(:)
 real(dp), allocatable :: ftr(:)
 real(dp) :: beta
 complex(dpc) :: xsto
 integer :: iflag,ifreq,itau,iwarn,log_direct
 character(len=500) :: message
 real(dp), allocatable :: omega_li(:)
! *********************************************************************
 beta=one/paw_dmft%temp
 iflag=0
 log_direct=1

 if(ldiag) iflag=1

! == inverse fourier transform
 if(opt_four==-1) then

   ABI_MALLOC(splined_li,(paw_dmft%dmft_nwli))
!   allocate(fw1(0:paw_dmft%dmft_nwlo-1))
   if(paw_dmft%dmft_log_freq==1) then
     ABI_MALLOC(omega_li,(1:paw_dmft%dmft_nwli))
     call construct_nwli_dmft(paw_dmft,paw_dmft%dmft_nwli,omega_li)
     call spline_c(paw_dmft%dmft_nwlo,paw_dmft%dmft_nwli,paw_dmft%omega_lo,&
&                   omega_li,splined_li,fw)
     ABI_FREE(omega_li)
   else
     splined_li=fw
   endif
   call invfourier(splined_li,ft,paw_dmft%dmft_nwli,ltau,iflag,beta)
!   deallocate(fw1)
   ABI_FREE(splined_li)

! == direct fourier transform
 else if(opt_four==1) then

   ABI_MALLOC(ftr,(ltau))

   iwarn=0
   do itau=1,ltau
     if(abs(aimag(ft(itau)))>tol12) then
       if(ldiag) then
         write(message,'(a,a,2(e15.4))') ch10,&
&          "green function is not real in imaginary time space",ft(itau)
         ABI_ERROR(message)
       else
         iwarn=iwarn+1
         ftr(itau)=real(ft(itau))
         xsto=ft(itau)
       endif
     else
       ftr(itau)=real(ft(itau))
     endif
!       write(std_out,*) itau,ftr(itau)
   enddo

   ABI_MALLOC(tospline_li,(paw_dmft%dmft_nwli))
   if (log_direct==1) then
!   do not have a physical meaning..because the log frequency is not
!   one of the linear frequency.
!   if log frequencies are also one of linear frequency, it should be good
     do ifreq=1, paw_dmft%dmft_nwlo
       call nfourier2(ftr,fw(ifreq),iflag,paw_dmft%omega_lo(ifreq),ltau,beta)
     enddo
   else
     call nfourier(ftr,tospline_li,iflag,paw_dmft%dmft_nwli-1,ltau,beta)
     do ifreq=1,paw_dmft%dmft_nwli
        write(1112,*) paw_dmft%temp*pi*real(2*ifreq-1,kind=dp),real(tospline_li(ifreq),kind=dp),aimag(tospline_li(ifreq))
        !write(1112,*) paw_dmft%omega_li(ifreq),real(tospline_li(ifreq)),aimag(tospline_li(ifreq))
     enddo
     if(paw_dmft%dmft_log_freq==1) then
       ABI_MALLOC(omega_li,(1:paw_dmft%dmft_nwli))
       call construct_nwli_dmft(paw_dmft,paw_dmft%dmft_nwli,omega_li)
       call spline_c(paw_dmft%dmft_nwli,paw_dmft%dmft_nwlo,omega_li,&
&                 paw_dmft%omega_lo,fw,tospline_li)
       ABI_FREE(omega_li)
     else
       fw=tospline_li
     endif
   endif

   ABI_FREE(tospline_li)

   ABI_FREE(ftr)
   if(iwarn>0) then
     write(message,'(a,a,2(e15.4))') ch10,&
&     "WARNING: off-diag green function is not real in imaginary time space",xsto
     call wrtout(std_out,message,'COLL')
   endif

 endif

end subroutine fourier_fct
!!***

!!****f* m_green/spline_fct
!! NAME
!! spline_fct
!!
!! FUNCTION
!! Do fourier transformation from matsubara space to imaginary time
!! (A spline is performed )
!!
!! COPYRIGHT
!! Copyright (C) 2006-2025 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  opt_spline = option for the spline
!!              -1 log to linear.
!!               1 linear to log.
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  fw= function is frequency space
!!  ft= function is time space
!!
!! SOURCE

subroutine spline_fct(fw1,fw2,opt_spline,paw_dmft)

!Arguments ------------------------------------
!type
 integer,intent(in) :: opt_spline
 type(paw_dmft_type), intent(in) :: paw_dmft
 complex(dpc), intent(inout) :: fw1(:)
 complex(dpc), intent(inout) :: fw2(:)
 integer :: size_fw1
 integer :: size_fw2
 real(dp), allocatable :: omega_li(:)

! *********************************************************************

 size_fw1 = size(fw1)
 size_fw2 = size(fw2)

! == inverse fourier transform
 if(opt_spline==-1) then

!   allocate(fw1(0:paw_dmft%dmft_nwlo-1))
   if(paw_dmft%dmft_log_freq==1) then
     ABI_MALLOC(omega_li,(1:size_fw2))
     call construct_nwli_dmft(paw_dmft,size_fw2,omega_li)
     call spline_c(size_fw1,size_fw2,paw_dmft%omega_lo(1:size_fw1),&
&                   omega_li(1:size_fw2),fw2,fw1)
     ABI_FREE(omega_li)
   else
     fw2=fw1
   endif

! == direct fourier transform
 else if(opt_spline==1) then


   if(paw_dmft%dmft_log_freq==1) then
     ABI_MALLOC(omega_li,(1:size_fw2))
     call construct_nwli_dmft(paw_dmft,size_fw2,omega_li)
     call spline_c(size_fw2,size_fw1,omega_li(1:size_fw2),&
&                 paw_dmft%omega_lo(1:size_fw1),fw1,fw2)
     ABI_FREE(omega_li)
   else
     fw1=fw2
   endif

 endif

end subroutine spline_fct
!!***

!!****f* m_green/occup_green_tau
!! NAME
!! occup_green_tau
!!
!! FUNCTION
!!  Compute occup_tau from green%oper_tau
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!
!! OUTPUT
!!
!! SOURCE

subroutine occup_green_tau(green)

!Arguments ------------------------------------
 type(green_type), intent(inout) :: green
!Local variables-------------------------------
 integer :: natom
 complex(dpc), allocatable :: shift(:)
! *********************************************************************

 natom = green%oper_tau(1)%natom

 ABI_MALLOC(shift,(natom))

 shift(:) = - cone

 call copy_matlu(green%oper_tau(1)%matlu(:),green%occup_tau%matlu(:),natom)
 call shift_matlu(green%occup_tau%matlu(:),natom,shift(:))

 ABI_FREE(shift)

end subroutine occup_green_tau
!!***

!!****f* m_green/occupfd
!! NAME
!! occupfd
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!!
!! SOURCE

 !function occupfd(eig,fermie,temp)


!Arguments ------------------------------------
!type
! Integrate analytic tail 1/(iw-mu)
 !real(dp),intent(in) :: eig,fermie,temp
 !real(dp) :: occupfd
!Local variables-------------------------------
! *********************************************************************

 !if ((eig-fermie) > zero) then
 !  occupfd = exp(-(eig-fermie)/temp)/(one+exp(-(eig-fermie)/temp))
 !else
 !  occupfd = one / (one+exp((eig-fermie)/temp))
 !end if

 !end function occupfd
!!***

!!****f* m_green/distrib_paral
!! NAME
!! distrib_paral
!!
!! FUNCTION
!!
!! INPUTS
!!  nw : number of frequencies
!!  nkpt : number of k-points
!!  nproc : number of procs
!!
!! OUTPUT
!!  procb(iw,ikpt) :  number of the proc  that is in charge of the combination {iw,ikpt}.
!!  proct(iw,iproc) : 1 if the frequency "iw" should be computed by the proc "iproc"
!!                    0 if iw should not   "      "
!!                    careful: if proct=1, it does not mean that each combination {iw,ikpt} is treated by this proc.
!! SOURCE

 subroutine distrib_paral(nkpt,nproc,nw,nw_perproc,procb,proct)

!Arguments ------------------------------------
!type
! Integrate analytic tail 1/(iw-mu)
 integer, intent(in) :: nw,nkpt,nproc
 integer, intent(out) :: nw_perproc
 integer, intent(out):: procb(nw,nkpt),proct(nw,0:nproc-1)
!Local variables-------------------------------
 integer :: ikpt,iw,ir,ratio,m,iproc
 integer, allocatable :: proca(:,:)
! *********************************************************************


 proct(:,:)=0
 procb(:,:)=-1

 if(nproc.ge.2*nw) then
 !write(6,*) "AA"
   ratio=nproc/nw
   ABI_MALLOC(proca,(nw,nproc/nw))
   do ir=1,ratio
     do iw=1,nw
       proca(iw,ir)=iw+(ir-1)*nw-1
       proct(iw,iw+(ir-1)*nw-1)=1
     enddo
   enddo
!     do iw=1,nw
!        write(6,*) " freq, procs"
!        write(6,*) iw,proca(iw,:)
!     enddo
   do iw=1,nw
     do ikpt=1,nkpt
       if(ikpt.le.ratio) procb(iw,ikpt)=proca(iw,ikpt)
       if(ikpt.ge.ratio) then
         m=mod(ikpt,ratio)
         if(m==0) then
            procb(iw,ikpt)=proca(iw,ratio)
         else
            procb(iw,ikpt)=proca(iw,m)
         endif
       endif
!       write(6,*) " freq, k-point, proc"
!       write(6,*) iw,ikpt,procb(iw,ikpt)
     enddo
   enddo
   nw_perproc=1
   ABI_FREE(proca)

 else if (nproc.ge.nw) then
 !write(6,*) "BB"
   do iw=1,nw
     procb(iw,:)= iw
     proct(iw,iw)=1
   enddo
!  do iw=1,nw
!        write(6,*) " freq, proc"
!        write(6,*) iw, procb(iw,1)
!  enddo
   nw_perproc=1

 else if (nproc.le.nw) then
 !write(6,*) "CC"
   ratio=nw/nproc
 !write(6,*) "ratio", ratio
   do iproc=0,nproc-1
     do iw=1,nw
       if (mod(iw-1,nproc)==iproc) then
         procb(iw,:)=iproc
         proct(iw,iproc)=1
       endif
     enddo
   enddo
    ! do iw=1,nw
    !   write(6,*) "iw, iproc", iw, procb(iw,1)
    ! enddo
   nw_perproc=ratio+1
!  some procs will compute a number of frequency which is ratio and some
!  other will compute a number of frequency which is ratio+1.

!  do iw=1,nw
!        write(6,*) " freq, proc"
!        write(6,*) iw, procb(iw,1)
!  enddo
  endif

!     do iw=1,nw
!      do iproc=0,nproc-1
!        write(6,*) " freq, procs,?"
!        write(6,*) iw,iproc,proct(iw,iproc)
!      enddo
!     enddo



 end subroutine distrib_paral
!!***

!!****f* ABINIT/greendftcompute_green
!! NAME
!! greendftcompute_green
!!
!! FUNCTION
!! Compute levels for ctqmc
!!
!! COPYRIGHT
!! Copyright (C) 1999-2025 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

 subroutine greendftcompute_green(green,paw_dmft)

!Arguments ------------------------------------
!scalars
 type(paw_dmft_type), intent(in)  :: paw_dmft
 type(green_type),intent(inout) :: green

!Local variables ------------------------------
! scalars
 integer :: ifreq,iband,ikpt,isppol
 integer :: mbandc,natom,nspinor,nsppol,nkpt
 character(len=500) :: message
! arrays
!************************************************************************

 mbandc=paw_dmft%mbandc
 nkpt=paw_dmft%nkpt
 nsppol=paw_dmft%nsppol
 nspinor=paw_dmft%nspinor
 natom=paw_dmft%natom

     !  write(6,*) green%oper(1)%ks(1,1,1,1)
 if(green%oper(1)%has_operks==0) then
  ABI_ERROR("greendft%oper(1)%ks not allocated")
 endif

!======================================
!Get Green's function G=1/(iw_n+mu-e_nk)
!======================================
 do ifreq=1,green%nw
   do iband=1,mbandc
     do ikpt=1,nkpt
       do isppol=1,nsppol
     !  write(6,*) green%oper(ifreq)%ks(isppol,ikpt,iband,iband)
     !  write(6,*) ifreq,iband,ikpt,isppol
     !  write(6,*) cmplx(0.d0,paw_dmft%omega_lo(ifreq),kind=dp)
     !  write(6,*) paw_dmft%fermie
     !  write(6,*) paw_dmft%eigen_dft(isppol,ikpt,iband)
         green%oper(ifreq)%ks(iband,iband,ikpt,isppol)=&
         cone/(cmplx(0.d0,paw_dmft%omega_lo(ifreq),kind=dp)+paw_dmft%fermie-paw_dmft%eigen_dft(iband,ikpt,isppol))
       end do
     end do
   end do


!======================================================================
! Compute local Green's function
!======================================================================
   call downfold_oper(green%oper(ifreq),paw_dmft)

!======================================================================
! Symetrize
!======================================================================
   call sym_matlu(green%oper(ifreq)%matlu,paw_dmft)
 enddo
 write(message,'(a,2x,a,f13.5)') ch10," == Print DFT Green's function for last frequency"
 call wrtout(std_out,message,'COLL')
 call print_matlu(green%oper(paw_dmft%dmft_nwlo)%matlu,natom,1)

 end subroutine greendftcompute_green
!!***

!!****f* m_green/fermi_green
!! NAME
!! fermi_green
!!
!! FUNCTION
!!  Compute Fermi level for DMFT or DFT.
!!
!! COPYRIGHT
!! Copyright (C) 2006-2025 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  self <type(self_type)>= variables related to self-energy
!!
!! OUTPUT
!!
!! SOURCE

subroutine fermi_green(green,paw_dmft,self)

!Arguments ------------------------------------
 type(green_type), intent(inout) :: green
 type(paw_dmft_type), intent(inout) :: paw_dmft
 !type(MPI_type), intent(in) :: mpi_enreg
 type(self_type), intent(inout) :: self
!Local variables-------------------------------
 integer :: ierr_hh,max_iter
 real(dp) :: f_precision,fermi_old,x_precision
! real(dp) :: hx
 character(len=13) :: tag
 character(len=500) :: message
!************************************************************************
!
 ABI_NVTX_START_RANGE(NVTX_DMFT_FERMI_GREEN)
 write(message,'(a,8x,a)') ch10,"  == Compute Fermi level"
 call wrtout(std_out,message,'COLL')

!=============
!headers
!=============
 write(message,'(2a)') ch10,"  |--- Newton method to search for Fermi level ------------|"
 call wrtout(std_out,message,'COLL')
 write(tag,'(f13.6)') paw_dmft%fermie
 write(message,'(3a)') ch10,"  |--- Initial value for Fermi level ",adjustl(tag)
 call wrtout(std_out,message,'COLL')

!========================================
!Define precision and nb of iterations
!=========================================
 fermi_old = paw_dmft%fermie
 ierr_hh = 0
 f_precision = paw_dmft%dmft_charge_prec
 !f_precision=0.01
 x_precision = tol5
!if(option==1) then
!f_precision=(erreursurlacharge)/100d0
!else
!f_precision=tol11
!endif
 !max_iter=1 ! for tests only
 !write(6,*) "for tests max_iter=1"
 max_iter = 50
 !opt_noninter = 4

 ! Precompute some useful quantities so that the update of the moments can be done in constant time at each iteration
 green%trace_fermie(1) = cmplx(paw_dmft%nsppol*paw_dmft%mbandc,zero,kind=dp)
 if (paw_dmft%nsppol == 1 .and. paw_dmft%nspinor == 1) green%trace_fermie(1) = &
               & green%trace_fermie(1) * two
 if (green%has_moments == 1) then
   call compute_trace_moments_ks(green,self,paw_dmft)
 end if

!=====================
!Call newton method
!=====================
 write(message,'(a,4x,a,e13.6)') ch10," Precision required :",f_precision
 call wrtout(std_out,message,'COLL')
 if (f_precision < ten) then
   call newton(green,self,paw_dmft,paw_dmft%fermie,x_precision,max_iter,f_precision,ierr_hh)
 end if

!===========================
!Deals with errors signals
!===========================
 if (ierr_hh == -314) then
   write(message,'(a)') "Warning, check Fermi level"
   call wrtout(std_out,message,'COLL')
!  call abi_abort('COLL')
   write(message,'(2a,f13.6)') ch10,"  |--- Final value for Fermi level (check)",paw_dmft%fermie
   call wrtout(std_out,message,'COLL')
 else if (ierr_hh == -123) then
   write(message,'(a,f13.6)') " Fermi level is set to",fermi_old
   paw_dmft%fermie = fermi_old
   call wrtout(std_out,message,'COLL')

!  =====================================
!  If fermi level search was successful
!  =====================================
 else
   write(tag,'(e13.6)') x_precision
   write(message,'(a,4x,2a)') ch10," Precision achieved on Fermi Level : ",adjustl(tag)
   call wrtout(std_out,message,'COLL')
   write(tag,'(e13.6)') f_precision
   write(message,'(4x,2a)') " Precision achieved on number of electrons : ",adjustl(tag)
   call wrtout(std_out,message,'COLL')
   write(tag,'(f13.6)') paw_dmft%fermie
   write(message,'(3a)') ch10,"  |--- Final value for Fermi level ",adjustl(tag)
   call wrtout(std_out,message,'COLL')
 end if ! ierr_hh

!========================================================
! Check convergence of fermi level during DMFT iterations
!========================================================
 if (paw_dmft%idmftloop >= 2) then
   if (abs(paw_dmft%fermie-fermi_old) <= paw_dmft%dmft_fermi_prec) then
!    write(message,'(a,8x,a,e9.2,a,8x,a,e12.5)') ch10,"|fermie(n)-fermie(n-1)|=<",paw_dmft%dmft_fermi_prec,ch10,&
     write(message,'(a,8x,a,e9.2,a,e9.2,a,8x,a,e12.5)') ch10,"|fermie(n)-fermie(n-1)|=",&
       & abs(paw_dmft%fermie-fermi_old),"<",paw_dmft%dmft_fermi_prec,ch10,&
       & "=> DMFT Loop: Fermi level is converged to:",paw_dmft%fermie
     call wrtout(std_out,message,'COLL')
     green%ifermie_cv = 1
   else
     write(tag,'(f12.5)') paw_dmft%fermie
     write(message,'(a,8x,2a)') ch10,"DMFT Loop: Fermi level is not converged: ",adjustl(tag)
     call wrtout(std_out,message,'COLL')
     green%ifermie_cv = 0
   end if ! convergence of Fermi level
 end if ! idmftloop>=2
 write(message,'(2a)') ch10, "  |---------------------------------------------------|"
 call wrtout(std_out,message,'COLL')
!
!==========================================================
!Recompute full green function including non diag elements
!==========================================================
 !call compute_green(green,paw_dmft,0,self,opt_self=1,opt_nonxsum=1)
 !call integrate_green(green,paw_dmft,prtopt=0,opt_ksloc=3) !,opt_nonxsum=1)

 ABI_NVTX_END_RANGE()
 !return
end subroutine fermi_green
!!***

!!****f* m_green/newton
!! NAME
!! newton
!!
!! FUNCTION
!!  Compute root of a function with newton methods (newton/halley)
!!
!! COPYRIGHT
!! Copyright (C) 2006-2025 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!  self <type(self_type)>= variables related to self-energy
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  x_input      : input value for x
!!  max_iter     : maximum number of iterations
!!  f_precision  : required precision on function F
!!  opt_algo     : 1 for Newton, 2 for Halley
!!
!! OUTPUT
!!  x_precision  : output precision on x
!!  ierr_hh      : different from zero if an error occurs
!!
!! SOURCE

subroutine newton(green,self,paw_dmft,x_input,x_precision,max_iter,&
                & f_precision,ierr_hh,opt_algo)

!Arguments ------------------------------------
 type(green_type), intent(inout) :: green
 type(self_type), intent(inout) :: self
 !type(MPI_type), intent(in) :: mpi_enreg
 type(paw_dmft_type), intent(inout) :: paw_dmft
 integer, intent(in) :: max_iter
 integer, intent(out) :: ierr_hh
 real(dp), intent(in) :: f_precision
 real(dp), intent(inout) :: x_input,x_precision
 integer, optional, intent(in) :: opt_algo
!Local variables-------------------------------
 integer :: iter,option
 logical :: dmft_optim,l_minus,l_plus
 real(dp) :: Fx,Fxdouble,Fxoptimum,Fxprime,nb_elec_x,step
 real(dp) :: x_minus,x_optimum,x_plus,xold
 character(len=500) :: message
! *********************************************************************

 x_minus = - dble(10)
 x_plus  = - dble(11)
 xold    = - dble(12)

 ierr_hh = 0
 option = 1
 if (present(opt_algo)) option = opt_algo
 step = paw_dmft%dmft_fermi_step

!write(std_out,*) "ldaprint",opt_noninter

!--- Start of iterations
 write(message,'(a,3a)') "     Fermi level   Charge    Difference"
 call wrtout(std_out,message,'COLL')
!do iter=1, 40
!x_input=float(iter)/100_dp
!call function_and_deriv(cryst_struc,f_precision,green,iter,mpi_enreg,paw_dmft,pawang,self,&
!&  x_input,x_before,x_precision,Fx,Fxprime,Fxdouble,opt_noninter,option)
!write(std_out,*) x_input,Fx
!enddo
!call abi_abort('COLL')

 l_minus = .false.
 l_plus  = .false.
 Fxoptimum = one
 x_optimum = zero

 dmft_optim = (paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7)

 if (dmft_optim) xold = x_input

!========================================
! Start iteration to find fermi level
!========================================
 do iter=1,max_iter

!  ========================================
!  If zero is located between two values: apply newton method or dichotomy
!  ========================================
   if ((l_minus .and. l_plus) .or. dmft_optim) then

!    ==============================================
!    Compute the function and derivatives for newton
!    ==============================================
     call function_and_deriv(green,self,paw_dmft,x_input,x_precision, &
                           & f_precision,Fx,Fxprime,Fxdouble,option)

!    Apply stop criterion on Fx
     if (abs(Fx) < f_precision) then
!      write(message,'(a,2f12.6)') "Fx,f_precision",Fx,f_precision
!      call wrtout(std_out,message,'COLL')
       x_precision = x_input - xold
       return
     end if ! abs(Fx)<f_precision
     if (iter == max_iter) then
       write(message,'(a,2f12.6)') "   Fermi level could not be found"
       call wrtout(std_out,message,'COLL')
       x_input = x_optimum
       ierr_hh = -123
       return
     end if ! iter=max_iter

!    Cannot divide by Fxprime if too small
     if (abs(Fxprime) <= tol15) then
       ierr_hh = -314
       write(message,'(a,f12.7)') "Fxprime=",Fxprime
       call wrtout(std_out,message,'COLL')
       return
     end if ! abs(Fxprime)<=tol15

     x_precision = x_input - xold

!    ==============================================
!    Newton/Halley's formula for next iteration
!    ==============================================
     xold = x_input
     if (option == 1) then
       if (dmft_optim) then
         x_input = x_input - sign(one,Fx)*merge(step,min(step,abs(Fx/Fxprime)),Fxprime<0)
       else
         x_input = x_input - Fx/Fxprime
       end if ! dmft_optim
     end if ! option=1
     if (option == 2) x_input = x_input - two*Fx*Fxprime/(two*(Fxprime**2)-Fx*Fxdouble)

!    ==============================================
!    If newton does not work well, use dichotomy.
!    ==============================================

     if (dmft_optim) then
       if (Fx < 0) then
         l_minus = .true.
         x_minus = xold
       end if
       if (Fx > 0) then
         l_plus = .true.
         x_plus = xold
       end if
     end if ! dmft_optim

     if ((x_input < x_minus .or. x_input > x_plus) .and. (l_minus .and. l_plus)) then

       if (.not. dmft_optim) then
         call compute_nb_elec(green,self,paw_dmft,Fx,nb_elec_x,xold)
       end if

       write(message,'(a,3f12.6)') " ---",x_input,Fx+paw_dmft%nelectval,Fx
       call wrtout(std_out,message,'COLL')
       if (.not. dmft_optim) then
         if (Fx > 0) then
           x_plus = xold
         else if (Fx < 0) then
           x_minus = xold
         end if ! Fx>0
       end if
       x_input = (x_plus+x_minus) * half

     end if ! x_input<x_minus or x_input>x_plus
!    write(std_out,'(a,2f12.6)') " Q(xold) and dQ/dx=",Fx,Fxprime
!    write(std_out,'(a,f12.6)') " =>  new Fermi level",x_input
!    ========================================
!    Locate zero between two values
!    ========================================
   else
     call compute_nb_elec(green,self,paw_dmft,Fx,nb_elec_x,x_input)
     write(message,'(a,3f12.6)') "  --",x_input,nb_elec_x,Fx
!    Possible improvement for large systems, removed temporarely for
!    automatic tests: more study is necessary: might worsen the convergency
!    if(iter==1) then
!    f_precision=max(abs(Fx/50),f_precision)
!    write(message,'(a,4x,a,e12.6)') ch10," Precision required changed to:",f_precision
!    call wrtout(std_out,message,'COLL')
!    endif
     call wrtout(std_out,message,'COLL')
     if (Fx > 0) then
       l_plus  = .true.
       x_plus  = x_input
       x_input = x_input - step
     else if (Fx < 0) then
       l_minus = .true.
       x_minus = x_input
       x_input = x_input + step
     end if ! Fx>0

   end if ! l_minus and l_plus

   if (abs(Fx) < abs(Fxoptimum) .or. (iter == 1 .and. dmft_optim)) then
     Fxoptimum = Fx
     x_optimum = merge(xold,x_input,dmft_optim)
   end if ! abs(Fx)<abs(Fxoptimum)



!  if(myid==master) then
!  write(std_out,'(a,i4,3f12.6)') "i,xnew,F,Fprime",i,x_input,Fx,Fxprime
!  endif


!  ! Apply stop criterion on x
!  if(abs(x_input-xold) .le. x_input*x_precision) then
!  !    write(std_out,'(a,4f12.6)') "x_input-xold, x_precision*x_input   "&
!  !&    ,x_input-xold,x_precision*x_input,x_precision
!  f_precision=Fx
!  return
!  endif

 end do ! iter
!--- End of iterations


 ierr_hh = 1
 !return

 CONTAINS  !========================================================================================
!-----------------------------------------------------------------------
!!***

!!****f* newton/function_and_deriv
!! NAME
!!  function_and_deriv
!!
!! FUNCTION
!!  Compute value of a function and its numerical derivatives
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!  self <type(self_type)>= variables related to self-energy
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  x_input = input value for x
!!  option = 1 to only compute the first derivative
!!         = 2 to compute the two first derivatives.
!!
!! OUTPUTS
!!  Fx           : Value of F(x)
!!  Fxprime      : Value of F'(x)
!!  Fxdouble     : Value of F''(x)
!!
!! SOURCE

subroutine function_and_deriv(green,self,paw_dmft,x_input,x_precision, &
                            & f_precision,Fx,Fxprime,Fxdouble,option)

!Arguments ------------------------------------
 type(green_type), intent(inout) :: green
 type(self_type), intent(inout) :: self
 !type(MPI_type), intent(in) :: mpi_enreg
 type(paw_dmft_type), intent(inout) :: paw_dmft
 integer, intent(in) :: option
 real(dp), intent(in) :: f_precision,x_input,x_precision
 real(dp), intent(out) :: Fx,Fxprime,Fxdouble
!Local variables-------------------------------
 logical :: dmft_optim
 real(dp) :: deltax,Fxminus,Fxplus,nb_elec_x,xminus,x0,xplus
 character(len=500) :: message
! *********************************************************************

   dmft_optim = (paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7)

   if (.not. dmft_optim) then

!  Choose deltax: for numeric evaluation of derivative
   !if(iter==1) then
!    deltax=0.02
   !end if
!  deltax=max((x_input-x_old)/10.d0,min(0.00001_dp,x_precision/100_dp))
     deltax = min(tol5,x_precision/dble(100))  ! small but efficient
!  endif
!  write(std_out,*) "iter,x_input,deltax",iter,x_input,deltax
     x0 = x_input
     xminus = x0 - deltax
     xplus  = x0 + deltax

     call compute_nb_elec(green,self,paw_dmft,Fx,nb_elec_x,x0)

     write(message,'(a,3f12.6)') "  - ",x0,nb_elec_x,Fx
     call wrtout(std_out,message,'COLL')

!  write(std_out,*) "Fx", Fx
     if (abs(Fx) < f_precision) return

     call compute_nb_elec(green,self,paw_dmft,Fxplus,nb_elec_x,xplus)

     write(message,'(a,3f12.6)') "  - ",xplus,nb_elec_x,Fxplus
     call wrtout(std_out,message,'COLL')

     if (option == 2) then
       call compute_nb_elec(green,self,paw_dmft,Fxminus,nb_elec_x,xminus)

       write(message,'(a,3f12.6)') "  - ",xminus,nb_elec_x,Fxminus
       call wrtout(std_out,message,'COLL')
     end if ! option=2

     if (option == 1) then
       Fxprime = (Fxplus-Fx) / deltax
     else if (option == 2) then
       Fxprime  = (Fxplus-Fxminus) / (two*deltax)
       Fxdouble = (Fxplus+Fxminus-two*Fx) / (deltax**2)
     end if ! option
!  write(std_out,*) "after computation of Fxprime",myid

   else

     call compute_nb_elec(green,self,paw_dmft,Fx,nb_elec_x,x_input,Fxprime=Fxprime)

     write(message,'(a,3f12.6)') "  - ",x_input,nb_elec_x,Fx
     call wrtout(std_out,message,'COLL')

   end if ! dmft_optim

   if (Fxprime < zero) then
     write(message,'(a,f12.6)') "  Warning: slope of charge versus fermi level is negative !",Fxprime
     call wrtout(std_out,message,'COLL')
   end if

 end subroutine function_and_deriv
!!***

!!****f* newton/compute_nb_elec
!! NAME
!! compute_nb_elec
!!
!! FUNCTION
!!  Compute nb of electrons as a function of Fermi level
!!
!! INPUTS
!! fermie       : input of energy
!! opt_noninter
!!
!! OUTPUTS
!! Fx           : Value of F(x).
!! nb_elec_x    : Number of electrons for the value of x
!!
!! SOURCE

subroutine compute_nb_elec(green,self,paw_dmft,Fx,nb_elec_x,fermie,Fxprime)

!Arguments ------------------------------------
 type(green_type), intent(inout) :: green
 type(self_type), intent(inout) :: self
 !type(MPI_type), intent(in) :: mpi_enreg
 type(paw_dmft_type), intent(inout) :: paw_dmft
 !integer,intent(in) :: opt_noninter
 real(dp), intent(in) :: fermie
 real(dp), intent(out) :: Fx,nb_elec_x
 real(dp), optional, intent(out) :: Fxprime
!Local variables-------------------------------
 integer :: band_index,i,ib,ierr,ifreq,ikpt,isppol,mbandc
 integer :: mkmem,nband_k,nkpt,nmoments,nspinor,nsppol,shift
 logical :: dmft_optim
 real(dp) :: correction,correction_prime,eig
 real(dp) :: fac,occ_prime,temp,wtk
 complex(dpc) :: omega
 type(oper_type) :: oper_tmp
 complex(dpc), allocatable :: omega_fac(:),trace_moments(:),trace_moments_prime(:)
! *********************************************************************

   ABI_NVTX_START_RANGE(NVTX_DMFT_COMPUTE_NB_ELEC)
   dmft_optim = (paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7)

   if (.not. dmft_optim) then

     paw_dmft%fermie = fermie
     call compute_green(green,paw_dmft,0,self,opt_self=1,opt_nonxsum=1,opt_nonxsum2=1)
     call integrate_green(green,paw_dmft,0,opt_ksloc=-1) !,opt_nonxsum=1)
        !  opt_ksloc=-1, compute total charge

   else

     green%charge_ks = zero

     mkmem = green%distrib%nkpt_mem(green%distrib%me_kpt+1)
     shift = green%distrib%shiftk
     nmoments = merge(green%nmoments,1,green%has_moments==1)

     mbandc  = paw_dmft%mbandc
     nkpt    = paw_dmft%nkpt
     nspinor = paw_dmft%nspinor
     nsppol  = paw_dmft%nsppol
     temp    = paw_dmft%temp

     ABI_MALLOC(trace_moments,(nmoments))
     ABI_MALLOC(trace_moments_prime,(nmoments))

     trace_moments(1) = green%trace_fermie(1)
     trace_moments_prime(:) = czero
     Fxprime = zero

     if (green%has_moments == 1) then
       call compute_trace_moments(fermie,green%trace_fermie(:),trace_moments(:),trace_moments_prime(:))
     end if

     call init_oper(paw_dmft,oper_tmp,nkpt=mkmem,shiftk=shift)

     do ifreq=1,green%nw
       if (green%distrib%proct(ifreq) /= green%distrib%me_freq) cycle
       omega = cmplx(zero,green%omega(ifreq),kind=dp)
       fac = two * paw_dmft%wgt_wlo(ifreq) * merge(two*temp,temp,nsppol==1.and.nspinor==1)

       call add_matlu(self%hdc%matlu(:),self%oper(ifreq)%matlu(:),oper_tmp%matlu(:),paw_dmft%natom,-1)
       call upfold_oper(oper_tmp,paw_dmft)

       do isppol=1,nsppol
         do ikpt=1,mkmem
           wtk = paw_dmft%wtk(ikpt+shift)
           do ib=1,mbandc
             oper_tmp%ks(ib,ib,ikpt,isppol) = oper_tmp%ks(ib,ib,ikpt,isppol) + omega + &
                & fermie - paw_dmft%eigen_dft(ib,ikpt+shift,isppol)
           end do ! ib
           call xginv(oper_tmp%ks(:,:,ikpt,isppol),mbandc)

           Fxprime = Fxprime - dble(sum(oper_tmp%ks(:,:,ikpt,isppol)*transpose(oper_tmp%ks(:,:,ikpt,isppol))))*wtk*fac
           do ib=1,mbandc
             green%charge_ks = green%charge_ks + dble(oper_tmp%ks(ib,ib,ikpt,isppol))*wtk*fac
           end do ! ib
         end do ! ikpt
       end do ! isppol

     end do ! ifreq

     call xmpi_sum(green%charge_ks,paw_dmft%spacecomm,ierr)
     call xmpi_sum(Fxprime,paw_dmft%spacecomm,ierr)

     call destroy_oper(oper_tmp)

     ABI_MALLOC(omega_fac,(nmoments))

     do i=1,nmoments
       omega_fac(i) = czero
       do ifreq=green%nw,1,-1 ! NEVER change the summation order and DON'T use the intrinsic SUM
         omega_fac(i) = omega_fac(i) + paw_dmft%wgt_wlo(ifreq) / (paw_dmft%omega_lo(ifreq))**i
       end do
       omega_fac(i) = - two * temp * omega_fac(i) / (j_dpc)**i
       if (i == 1) omega_fac(i) = omega_fac(i) + half
       if (i == 2) omega_fac(i) = omega_fac(i) - cone/(four*temp)
       if (i == 4) omega_fac(i) = omega_fac(i) + cone/(dble(48)*(temp**3))
     end do ! i

     ! Do not use DOT_PRODUCT
     green%charge_ks = green%charge_ks + dble(sum(trace_moments(1:nmoments)*omega_fac(1:nmoments)))
     ! Do not use DOT_PRODUCT
     Fxprime = Fxprime + dble(sum(trace_moments_prime(1:nmoments)*omega_fac(1:nmoments)))

     ABI_FREE(trace_moments)
     ABI_FREE(trace_moments_prime)
     ABI_FREE(omega_fac)

     if (paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7) then
       band_index = 0
       correction = zero
       correction_prime = zero
       do isppol=1,nsppol
         do ikpt=1,nkpt
           nband_k = paw_dmft%nband(ikpt+(isppol-1)*nkpt)
           wtk = paw_dmft%wtk(ikpt)
           do ib=1,nband_k
             if (paw_dmft%band_in(ib)) cycle
             eig = paw_dmft%eigen(ib+band_index)
             correction = correction + occup_fd(eig,fermie,temp)*wtk
             if ((eig-fermie) > zero) then
               occ_prime = exp(-(eig-fermie)/temp) / (one+exp(-(eig-fermie)/temp))**2 / temp
             else
               occ_prime = exp((eig-fermie)/temp) / (one+exp((eig-fermie)/temp))**2 / temp
             end if
             correction_prime = correction_prime + occ_prime*wtk
           end do ! ib
           band_index = band_index + nband_k
         end do ! ikpt
       end do ! isppol
       if (nsppol == 1 .and. nspinor == 1) correction = correction * two
       green%charge_ks = green%charge_ks + correction
       if (nsppol == 1 .and. nspinor == 1) correction_prime = correction_prime * two
       Fxprime = Fxprime + correction_prime
     end if ! use_all_bands

   end if ! dmft_optim

   nb_elec_x = green%charge_ks
   Fx = green%charge_ks - paw_dmft%nelectval

   ABI_NVTX_END_RANGE()
 end subroutine compute_nb_elec
!!***

subroutine compute_trace_moments(fermie,trace_fermie,trace_moments,trace_moments_prime)

!Arguments ------------------------------------
 real(dp), intent(in) :: fermie
 complex(dpc), intent(in) :: trace_fermie(:)
 complex(dpc), intent(inout) :: trace_moments(:),trace_moments_prime(:)
!Local variables-------------------------------
 real(dp) :: fermie2,fermie3,fermie4
! *********************************************

  fermie2 = fermie * fermie
  fermie3 = fermie2 * fermie
  fermie4 = fermie3 * fermie

  trace_moments(2) = trace_fermie(2) - fermie*trace_fermie(1)
  trace_moments(3) = trace_fermie(3) + trace_fermie(4) - two*fermie*trace_fermie(2) + &
         & fermie2*trace_fermie(1)
  trace_moments(4) = trace_fermie(5) + two*(trace_fermie(6)-fermie*trace_fermie(3)) + &
         & trace_fermie(7) - three*fermie*trace_fermie(4) + &
         & three*fermie2*trace_fermie(2) - fermie3*trace_fermie(1)
  trace_moments(5) = trace_fermie(8) + two*(trace_fermie(9)-fermie*trace_fermie(5)) + &
         & trace_fermie(10) + three*(trace_fermie(11)-two*fermie*trace_fermie(6)+ &
         & fermie2*trace_fermie(3)) + trace_fermie(12) - four*fermie*trace_fermie(7) + &
         & six*fermie2*trace_fermie(4) - four*fermie3*trace_fermie(2) + &
         & fermie4*trace_fermie(1)

  trace_moments_prime(2) = - trace_fermie(1)
  trace_moments_prime(3) = two * (fermie*trace_fermie(1)-trace_fermie(2))
  trace_moments_prime(4) = - two*trace_fermie(3) - three*trace_fermie(4) + &
         & six*fermie*trace_fermie(2) - three*fermie2*trace_fermie(1)
  trace_moments_prime(5) = - two*trace_fermie(5) + three*(-two*trace_fermie(6)+ &
         & two*fermie*trace_fermie(3)) - four*trace_fermie(7) + &
         & dble(12)*fermie*trace_fermie(4) - dble(12)*fermie2*trace_fermie(2) + &
         & four*fermie3*trace_fermie(1)

 end subroutine compute_trace_moments
!!***

end subroutine newton
!!***

!!****f* m_green/local_ks_green
!! NAME
!! local_ks_green
!!
!! FUNCTION
!! Compute the sum over k-point of ks green function.
!! do the fourier transformation and print it
!!
!! COPYRIGHT
!! Copyright (C) 1999-2025 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cryst_struc
!!  istep    =  step of iteration for DFT.
!!  dft_occup
!!  mpi_enreg=information about MPI parallelization
!!  paw_dmft =  data for self-consistent DFT+DMFT calculations.
!!
!! OUTPUT
!!  paw_dmft =  data for self-consistent DFT+DMFT calculations.
!!
!! SOURCE

subroutine local_ks_green(green,paw_dmft,prtopt)

!Arguments ------------------------------------
!scalars
 type(green_type), intent(in) :: green
 type(paw_dmft_type), intent(in)  :: paw_dmft
 integer, intent(in) :: prtopt

!Local variables ------------------------------
 character(len=500) :: message
 integer :: iband,ifreq,ikpt,isppol,itau,lsub,ltau,mbandc,nkpt,nsppol
 character(len=1) :: tag_is
 character(len=fnlen) :: tmpfil
 integer,allocatable :: unitgreenlocks_arr(:)
 real(dp) :: beta
 real(dp), allocatable :: tau(:)
 complex(dpc), allocatable :: loc_ks(:,:,:)
 complex(dpc), allocatable :: loc_ks_tau(:,:,:),fw(:),ft(:)
!scalars
!************************************************************************
 mbandc=paw_dmft%mbandc
 nkpt=paw_dmft%nkpt
 nsppol=paw_dmft%nsppol
 ltau=128
 ABI_MALLOC(tau,(ltau))
 do itau=1,ltau
   tau(itau)=float(itau-1)/float(ltau)/paw_dmft%temp
 end do
 beta=one/paw_dmft%temp

!Only imaginary frequencies here
 if(green%w_type=="real") then
   message = ' compute_energy not implemented for real frequency'
   ABI_BUG(message)
 end if

!=========================================
!Compute local band ks green function
! should be computed in compute_green: it would be less costly in memory.
!=========================================
 ABI_MALLOC(loc_ks,(nsppol,mbandc,paw_dmft%dmft_nwlo))
 if(green%oper(1)%has_operks==1) then
   loc_ks(:,:,:)=czero
   do isppol=1,nsppol
     do iband=1,mbandc
       do ifreq=1,paw_dmft%dmft_nwlo
         do ikpt=1,nkpt
           loc_ks(isppol,iband,ifreq)=loc_ks(isppol,iband,ifreq)+  &
&           green%oper(ifreq)%ks(iband,iband,ikpt,isppol)*paw_dmft%wtk(ikpt)
         end do
       end do
     end do
   end do
 else
   message = ' green fct is not computed in ks space'
   ABI_BUG(message)
 end if

!=========================================
!Compute fourier transformation
!=========================================

 ABI_MALLOC(loc_ks_tau,(nsppol,mbandc,ltau))
 ABI_MALLOC(fw,(paw_dmft%dmft_nwlo))
 ABI_MALLOC(ft,(ltau))
 loc_ks_tau(:,:,:)=czero
 do isppol=1,nsppol
   do iband=1,mbandc
     do ifreq=1,paw_dmft%dmft_nwlo
       fw(ifreq)=loc_ks(isppol,iband,ifreq)
     end do
     call fourier_fct(fw,ft,.true.,ltau,-1,paw_dmft) ! inverse fourier
     do itau=1,ltau
       loc_ks_tau(isppol,iband,itau)=ft(itau)
     end do
   end do
 end do
 ABI_FREE(fw)
 ABI_FREE(ft)
 do isppol=1,nsppol
   do iband=1,mbandc
     do itau=1,ltau
       loc_ks_tau(isppol,iband,itau)=(loc_ks_tau(isppol,iband,itau)+conjg(loc_ks_tau(isppol,iband,itau)))/two
     end do
   end do
 end do

!=========================================
!Print out ksloc green function
!=========================================
 if(abs(prtopt)==1) then
   ABI_MALLOC(unitgreenlocks_arr,(nsppol))
   do isppol=1,nsppol
     write(tag_is,'(i1)')isppol
     tmpfil = trim(paw_dmft%filapp)//'Gtau_locks_isppol'//tag_is
     write(message,'(3a)') ch10," == Print green function on file ",tmpfil
     call wrtout(std_out,message,'COLL')
     unitgreenlocks_arr(isppol)=500+isppol-1
     open (unit=unitgreenlocks_arr(isppol),file=trim(tmpfil),status='unknown',form='formatted')
     rewind(unitgreenlocks_arr(isppol))
     write(message,'(a,a,a,i4)') 'opened file : ', trim(tmpfil), ' unit', unitgreenlocks_arr(isppol)
     write(message,'(a,a)') ch10,"# New record : First 40 bands"
     call wrtout(unitgreenlocks_arr(isppol),message,'COLL')
     do lsub=1,mbandc/40+1
       do itau=1, ltau
         write(message,'(2x,50(e10.3,2x))') tau(itau), &
&         (real(loc_ks_tau(isppol,iband,itau)),iband=40*(lsub-1)+1,min(40*lsub,mbandc))
         call wrtout(unitgreenlocks_arr(isppol),message,'COLL')
       end do
       write(message,'(2x,50(e10.3,2x))') beta, &
&       ((-one-real(loc_ks_tau(isppol,iband,1))),iband=40*(lsub-1)+1,min(40*lsub,mbandc))
       call wrtout(unitgreenlocks_arr(isppol),message,'COLL')
       if(40*lsub<mbandc) then
         write(message,'(a,a,i5,a,i5)')    &
&         ch10,"# Same record, Following bands : From ",    &
&         40*(lsub),"  to ",min(40*(lsub+1),mbandc)
         call wrtout(unitgreenlocks_arr(isppol),message,'COLL')
       end if
     end do
!    call flush(unitgreenlocks_arr(isppol))
   end do
   ABI_FREE(unitgreenlocks_arr)
 end if

!Deallocations
 ABI_FREE(loc_ks)
 ABI_FREE(loc_ks_tau)
 ABI_FREE(tau)

end subroutine local_ks_green
!!***

!!****f* m_green/compute_moments_ks
!! NAME
!! compute_moments_ks
!!
!! FUNCTION
!! Compute the high-frequency moments of the Green's function in KS space
!!
!! COPYRIGHT
!! Copyright (C) 1999-2024 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!  self <type(self_type)>= variables related to self-energy
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  opt_self = 0 : do not use the self-energy
!!           = 1 (default) : use the self-energy, by upfolding its moments
!!  opt_log = 0 (default): the moments of G are computed in green%moments using self%moments
!!          = 1 : the trace of the moments of log(G)+log(iw*Id)=-log(Id+(mu-sigma)/iw) is computed in green%trace_moments_log_ks
!!  opt_quick_restart = 0 (default) : default behavior, recompute everything from scratch
!!                    = 1 : use precomputed quantities from fermi_green : assume that self%moments(i>=2) have
!!                    already been upfolded, that self%moments(1)-self%hdc has already been upfolded in
!!                    green%moments(2), and that (self%moments(1)-self%hdc+eigen_dft)**2 has been computed
!!                    in green%moments(3)
!!  The conventions are: green%moments(i) -> ith order moment of the Green's function in the usual case
!!                                            (the 1st one is not computed, since this is simply the identity)
!!                       self%moments(i)  -> (i-1)th order moment of the self-energy (the 1st one is not upfolded, since we directly
!!                                  upfold self%moments(1) - self%hdc
!!  where the ith order moment corresponds to the factor in front of 1/(jw_n)^i in the asymptotic expansion
!!
!! OUTPUT
!!
!! SOURCE

subroutine compute_moments_ks(green,self,paw_dmft,opt_self,opt_log,opt_quick_restart)

!Arguments ------------------------------------
!scalars
 type(green_type), intent(inout) :: green
 type(self_type), intent(inout) :: self
 type(paw_dmft_type), intent(in)  :: paw_dmft
 integer, optional, intent(in) :: opt_log,opt_quick_restart,opt_self
!Local variables ------------------------------
 integer :: diag,i,ib,ierr,mkmem,natom,nsppol
 integer :: optlog,optquickrestart,optself,shift
 real(dp) :: dum,mu,mu2,mu3,mu4
 complex(dpc) :: trace_tmp
 type(oper_type) :: oper(2)
 real(dp), allocatable :: trace_loc(:,:)
 character(len=500) :: message
!************************************************************************

 optlog = 0
 if (present(opt_log)) optlog = opt_log

 optself = 1
 if (present(opt_self)) optself = opt_self

 optquickrestart = 0
 if (present(opt_quick_restart)) optquickrestart = opt_quick_restart

 if (optself == 0 .and. optquickrestart == 1) then
   message = "Case optself=0 and optquickrestart=1 is not supposed to happen"
   ABI_ERROR(message)
 end if

 diag   = 1 - optself
 mkmem  = green%moments(1)%nkpt
 mu     = paw_dmft%fermie
 mu2    = mu * mu
 mu3    = mu2 * mu
 mu4    = mu3 * mu
 natom  = paw_dmft%natom
 nsppol = paw_dmft%nsppol
 shift  = green%moments(1)%shiftk

 if (optlog == 1 .and. optquickrestart == 0) then
   ABI_MALLOC(trace_loc,(nsppol+1,natom))
 end if

 if (optself == 1 .and. optquickrestart == 0) then
   call add_matlu(self%moments(1)%matlu(:),self%hdc%matlu(:),green%moments(2)%matlu(:),natom,-1)
   call upfold_oper(green%moments(2),paw_dmft)
   do i=2,self%nmoments
     call upfold_oper(self%moments(i),paw_dmft)
   end do ! i
 end if ! optself=1

 do ib=1,paw_dmft%mbandc
   if (optself == 1) then
     green%moments(2)%ks(ib,ib,:,:) = green%moments(2)%ks(ib,ib,:,:) + &
            & paw_dmft%eigen_dft(ib,1+shift:mkmem+shift,:) - mu
   else
     green%moments(2)%ks(ib,ib,:,:) = paw_dmft%eigen_dft(ib,1+shift:mkmem+shift,:) - mu
   end if
   ! Careful, we need to substract mu**2 below, not add it (think about it)
   if (optquickrestart == 1) green%moments(3)%ks(ib,ib,:,:) = green%moments(3)%ks(ib,ib,:,:) - mu2
 end do ! ib

 if (optquickrestart == 1) then
   green%moments(3)%ks(:,:,:,:) = green%moments(3)%ks(:,:,:,:) - two*mu*green%moments(2)%ks(:,:,:,:)
 else
   call prod_oper(green%moments(2),green%moments(2),green%moments(3),1,opt_diag=diag)
 end if ! optquickrestart
 do i=3,green%nmoments-1
   call prod_oper(green%moments(2),green%moments(i),green%moments(i+1),1,opt_diag=diag)
 end do ! i

 if (optlog == 1 .and. optquickrestart == 0) then
   do i=1,green%nmoments-1
     call trace_oper(green%moments(i+1),dum,trace_loc(:,:),1,trace_ks_cmplx=trace_tmp)
     green%trace_moments_log_ks(i) = trace_tmp / dble(i)
   end do ! i
   if (optself == 1) then
     do i=2,self%nmoments
       call trace_oper(self%moments(i),dum,trace_loc(:,:),1,trace_ks_cmplx=trace_tmp)
       green%trace_moments_log_ks(i) = green%trace_moments_log_ks(i) + trace_tmp
     end do ! i
   end if ! optself=1
 end if ! optlog=1 and optquickrestart=0

 if (optlog == 1 .and. optquickrestart == 1) then
   green%trace_moments_log_ks(1) = green%trace_fermie(2) - mu*green%trace_fermie(1)
   green%trace_moments_log_ks(2) = green%trace_fermie(3) + &
       & half*(green%trace_fermie(4)-two*mu*green%trace_fermie(2)+mu2*green%trace_fermie(1))
   green%trace_moments_log_ks(3) = green%trace_fermie(5) + green%trace_fermie(6) - &
      & mu*green%trace_fermie(3) + third*(green%trace_fermie(7)+ &
      & three*mu2*green%trace_fermie(2)-three*mu*green%trace_fermie(4)- &
      & mu3*green%trace_fermie(1))
   green%trace_moments_log_ks(4) = green%trace_fermie(8) + green%trace_fermie(9) - &
      & mu*green%trace_fermie(5) + half*green%trace_fermie(10) + &
      & green%trace_fermie(11) - two*mu*green%trace_fermie(6) + mu2*green%trace_fermie(3) + &
      & quarter*(green%trace_fermie(12)-four*mu3*green%trace_fermie(2)+ &
      & six*mu2*green%trace_fermie(4)-four*mu*green%trace_fermie(7)+ &
      & mu4*green%trace_fermie(1))
 end if ! optlog=1 and optquickrestart=1

 if (optself == 1) then
   do i=1,2
     call init_oper(paw_dmft,oper(i),nkpt=mkmem,shiftk=shift,opt_ksloc=1)
   end do
   call prod_oper(green%moments(2),self%moments(2),oper(1),1)
   if (optlog == 1 .and. optquickrestart == 0) then
     call trace_oper(oper(1),dum,trace_loc(:,:),1,trace_ks_cmplx=trace_tmp)
     green%trace_moments_log_ks(3) = green%trace_moments_log_ks(3) + trace_tmp
   end if
   call prod_oper(self%moments(2),green%moments(2),oper(2),1)
   oper(2)%ks(:,:,:,:) = oper(2)%ks(:,:,:,:) + oper(1)%ks(:,:,:,:)
   green%moments(4)%ks(:,:,:,:) = oper(2)%ks(:,:,:,:) + green%moments(4)%ks(:,:,:,:)
   call prod_oper(oper(2),green%moments(2),oper(1),1)
   call prod_oper(green%moments(3),self%moments(2),oper(2),1)
   if (optlog == 1 .and. optquickrestart == 0) then
     call trace_oper(oper(2),dum,trace_loc(:,:),1,trace_ks_cmplx=trace_tmp)
     green%trace_moments_log_ks(4) = green%trace_moments_log_ks(4) + trace_tmp
   end if
   oper(2)%ks(:,:,:,:) = oper(2)%ks(:,:,:,:) + oper(1)%ks(:,:,:,:)
   call prod_oper(green%moments(2),self%moments(3),oper(1),1)
   if (optlog == 1 .and. optquickrestart == 0) then
     call trace_oper(oper(1),dum,trace_loc(:,:),1,trace_ks_cmplx=trace_tmp)
     green%trace_moments_log_ks(4) = green%trace_moments_log_ks(4) + trace_tmp
   end if
   oper(2)%ks(:,:,:,:) = oper(2)%ks(:,:,:,:) + oper(1)%ks(:,:,:,:)
   call prod_oper(self%moments(3),green%moments(2),oper(1),1)
   oper(2)%ks(:,:,:,:) = oper(2)%ks(:,:,:,:) + oper(1)%ks(:,:,:,:)
   call prod_oper(self%moments(2),self%moments(2),oper(1),1)
   if (optlog == 1 .and. optquickrestart == 0) then
     call trace_oper(oper(1),dum,trace_loc(:,:),1,trace_ks_cmplx=trace_tmp)
     green%trace_moments_log_ks(4) = green%trace_moments_log_ks(4) + trace_tmp*half
   end if
   oper(2)%ks(:,:,:,:) = oper(2)%ks(:,:,:,:) + oper(1)%ks(:,:,:,:)
   green%moments(5)%ks(:,:,:,:) = oper(2)%ks(:,:,:,:) + green%moments(5)%ks(:,:,:,:)
   do i=2,self%nmoments
     green%moments(i+1)%ks(:,:,:,:) = green%moments(i+1)%ks(:,:,:,:) + self%moments(i)%ks(:,:,:,:)
   end do ! i
   do i=1,2
     call destroy_oper(oper(i))
   end do
 end if ! optself>=1

 if (optlog == 1 .and. optquickrestart == 0) then
   call xmpi_sum(green%trace_moments_log_ks(1:green%nmoments-1),green%distrib%comm_kpt,ierr)
 end if

 ABI_SFREE(trace_loc)

end subroutine compute_moments_ks
!!***

!!****f* m_green/compute_trace_moments_ks
!! NAME
!! compute_trace_moments_ks
!!
!! FUNCTION
!! Precompute some useful quantities for the calculation of the trace of the moments (for
!! Fermi level search)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2024 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!  self <type(self_type)>= variables related to self-energy
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  The conventions are:
!!    trace_fermie(1)  = nsppol*mbandc (=tr(Id))
!!    trace_fermie(2)  = tr(m0) (with m0=self%moments(1)-self%hdc+eigen_dft*Id
!!    trace_fermie(3)  = tr(self%moments(2))
!!    trace_fermie(4)  = tr(m0**2)
!!    trace_fermie(5)  = tr(self%moments(3))
!!    trace_fermie(6)  = tr(m0*self%moments(2))
!!    trace_fermie(7)  = tr(m0**3)
!!    trace_fermie(8)  = tr(self%moments(4))
!!    trace_fermie(9)  = tr(m0*self%moments(3))
!!    trace_fermie(10) = tr(self%moments(2)**2)
!!    trace_fermie(11) = tr(m0**2*self%moments(2))
!!    trace_fermie(12) = tr(m0**4)
!!
!! OUTPUT
!!
!! SOURCE

subroutine compute_trace_moments_ks(green,self,paw_dmft)

!Arguments ------------------------------------
!scalars
 type(green_type), intent(inout) :: green
 type(self_type), intent(inout) :: self
 type(paw_dmft_type), intent(in) :: paw_dmft
!Local variables ------------------------------
 integer :: i,ib,ierr,mkmem,natom,nsppol,shift
 real(dp) :: dum
 type(oper_type) :: oper_tmp
 real(dp), allocatable :: trace_loc(:,:)
!************************************************************************

 mkmem  = green%distrib%nkpt_mem(green%distrib%me_kpt+1)
 natom  = paw_dmft%natom
 nsppol = paw_dmft%nsppol
 shift  = green%distrib%shiftk

 ABI_MALLOC(trace_loc,(nsppol+1,natom))

 call init_oper(paw_dmft,oper_tmp,nkpt=mkmem,shiftk=shift,opt_ksloc=1)

 call add_matlu(self%moments(1)%matlu(:),self%hdc%matlu(:),green%moments(2)%matlu(:),natom,-1)
 call upfold_oper(green%moments(2),paw_dmft)
 do i=2,self%nmoments
   call upfold_oper(self%moments(i),paw_dmft)
 end do ! i

 oper_tmp%ks(:,:,:,:) = green%moments(2)%ks(:,:,:,:)

 do ib=1,paw_dmft%mbandc
   oper_tmp%ks(ib,ib,:,:) = oper_tmp%ks(ib,ib,:,:) + paw_dmft%eigen_dft(ib,1+shift:mkmem+shift,:)
 end do ! ib

 call prod_oper(oper_tmp,oper_tmp,green%moments(3),1)
 call trace_oper(oper_tmp,dum,trace_loc(:,:),1,trace_ks_cmplx=green%trace_fermie(2))
 call trace_oper(self%moments(2),dum,trace_loc(:,:),1,trace_ks_cmplx=green%trace_fermie(3))
 call trace_oper(green%moments(3),dum,trace_loc(:,:),1,trace_ks_cmplx=green%trace_fermie(4))
 call trace_oper(self%moments(3),dum,trace_loc(:,:),1,trace_ks_cmplx=green%trace_fermie(5))
 call trace_prod_oper(oper_tmp,self%moments(2),green%trace_fermie(6))
 call trace_prod_oper(oper_tmp,green%moments(3),green%trace_fermie(7))
 call trace_oper(self%moments(4),dum,trace_loc(:,:),1,trace_ks_cmplx=green%trace_fermie(8))
 call trace_prod_oper(oper_tmp,self%moments(3),green%trace_fermie(9))
 call trace_prod_oper(self%moments(2),self%moments(2),green%trace_fermie(10))
 call trace_prod_oper(green%moments(3),self%moments(2),green%trace_fermie(11))
 call trace_prod_oper(green%moments(3),green%moments(3),green%trace_fermie(12))

 call destroy_oper(oper_tmp)
 ABI_FREE(trace_loc)

 call xmpi_sum(green%trace_fermie(2:12),green%distrib%comm_kpt,ierr)

end subroutine compute_trace_moments_ks
!!***

!!****f* m_green/compute_moments_loc
!! NAME
!! compute_moments_loc
!!
!! FUNCTION
!!  Computes the high-frequency moments in local space.
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!  self <type(self_type)>= variables related to self-energy
!!  energy_level = local energy levels
!!  weiss <type(green_type)>= weiss field
!!  option = 0 : computes the hybridization moments in weiss%moments up to 4th order (including
!!               the spurious 0th order moment, which is removed from the weiss field later),
!!               using the moments of the self-energy and Green's function as inputs
!!         = 1 : compute the self-energy moments in self%moments using the moments of the hybridization
!!               and the Green's function (assuming the spurious 0th order moment of the hybridization
!!               has been removed)
!!  opt_log = if set to 1, also computes the trace of the moments of log(G)+log(iw*Id) and
!!            log(G0)+log(iw*Id) in green%trace_moments_log_loc and weiss%trace_moments_log_loc
!!
!! OUTPUTS
!!
!! SOURCE

subroutine compute_moments_loc(green,self,energy_level,weiss,option,opt_log)

!Arguments ------------------------------------
 type(green_type), intent(inout) :: green,weiss
 type(self_type), intent(inout) :: self
 type(oper_type), intent(in) :: energy_level
 integer, intent(in) :: option
 integer, optional, intent(in) :: opt_log
!Local variables ------------------------------
 integer :: i,natom,nspinor,nsppol,optlog
 complex(dpc) :: trace
 integer, allocatable :: lpawu(:)
 complex(dpc), allocatable :: trace_loc(:)
 type(matlu_type), allocatable :: matlu(:,:)
!************************************************************************

 optlog = 0
 if (present(opt_log)) optlog = opt_log

 natom   = energy_level%natom
 nspinor = energy_level%nspinor
 nsppol  = energy_level%nsppol

 ABI_MALLOC(trace_loc,(natom))

 ABI_MALLOC(lpawu,(natom))
 ABI_MALLOC(matlu,(natom,6))

 lpawu(:) = energy_level%matlu(:)%lpawu ! to prevent the creation of a temporary array when calling init_matlu

 do i=1,6
   call init_matlu(natom,nspinor,nsppol,lpawu(:),matlu(:,i))
 end do

 if (option == 0) then
   call add_matlu(green%moments(2)%matlu(:),energy_level%matlu(:),matlu(:,1),natom,-1)
   call add_matlu(matlu(:,1),self%moments(1)%matlu(:),weiss%moments(1)%matlu(:),natom,-1)
 else
   call add_matlu(green%moments(2)%matlu(:),energy_level%matlu(:),self%moments(1)%matlu(:),natom,-1)
 end if ! option

 if (optlog == 1) then
   call trace_matlu(energy_level%matlu(:),natom,itau=0,trace=weiss%trace_moments_log_loc(1))
   do i=2,green%nmoments-1
     call trace_matlu(weiss%moments(i)%matlu(:),natom,itau=0,trace=weiss%trace_moments_log_loc(i))
   end do ! i
   call prod_matlu(energy_level%matlu(:),energy_level%matlu(:),matlu(:,1),natom)
   call trace_matlu(matlu(:,1),natom,itau=0,trace=trace)
   weiss%trace_moments_log_loc(2) = weiss%trace_moments_log_loc(2) + trace*half
   call trace_prod_matlu(energy_level%matlu(:),weiss%moments(2)%matlu(:),natom,trace_loc(:),trace_tot=trace)
   weiss%trace_moments_log_loc(3) = weiss%trace_moments_log_loc(3) + trace
   call trace_prod_matlu(energy_level%matlu(:),matlu(:,1),natom,trace_loc(:),trace_tot=trace)
   weiss%trace_moments_log_loc(3) = weiss%trace_moments_log_loc(3) + trace*third
   call trace_prod_matlu(energy_level%matlu(:),weiss%moments(3)%matlu(:),natom,trace_loc(:),trace_tot=trace)
   weiss%trace_moments_log_loc(4) = weiss%trace_moments_log_loc(4) + trace
   call trace_prod_matlu(weiss%moments(2)%matlu(:),weiss%moments(2)%matlu(:),natom,trace_loc(:),trace_tot=trace)
   weiss%trace_moments_log_loc(4) = weiss%trace_moments_log_loc(4) + trace*half
   call trace_prod_matlu(weiss%moments(2)%matlu(:),matlu(:,1),natom,trace_loc(:),trace_tot=trace)
   weiss%trace_moments_log_loc(4) = weiss%trace_moments_log_loc(4) + trace
   call trace_prod_matlu(matlu(:,1),matlu(:,1),natom,trace_loc(:),trace_tot=trace)
   weiss%trace_moments_log_loc(4) = weiss%trace_moments_log_loc(4) + trace*quarter

   call trace_matlu(green%moments(2)%matlu(:),natom,itau=0,trace=green%trace_moments_log_loc(1))
 end if ! optlog=1

 call prod_matlu(green%moments(2)%matlu(:),green%moments(2)%matlu(:),matlu(:,1),natom) ! matlu=m0m0

 call add_matlu(green%moments(3)%matlu(:),matlu(:,1),matlu(:,2),natom,-1) ! matlu2=m1

 if (optlog == 1) then
   call trace_matlu(matlu(:,2),natom,itau=0,trace=green%trace_moments_log_loc(2))
   call trace_matlu(matlu(:,1),natom,itau=0,trace=trace)
   green%trace_moments_log_loc(2) = green%trace_moments_log_loc(2) + half*trace
 end if ! optlog

 if (option == 0) then
   call add_matlu(matlu(:,2),self%moments(2)%matlu(:),weiss%moments(2)%matlu(:),natom,-1)
 else if (option == 1) then
   call add_matlu(matlu(:,2),weiss%moments(2)%matlu(:),self%moments(2)%matlu(:),natom,-1)
 end if ! option

 call prod_matlu(green%moments(2)%matlu(:),matlu(:,1),matlu(:,3),natom) ! matlu3=m0m0m0
 call add_matlu(green%moments(4)%matlu(:),matlu(:,3),matlu(:,1),natom,-1) ! matlu=green%moments(4)-m0m0m0

 if (optlog == 1) then
   call trace_matlu(matlu(:,3),natom,itau=0,trace=trace)
   green%trace_moments_log_loc(3) = trace * third
 end if ! optlog

 call prod_matlu(green%moments(2)%matlu(:),matlu(:,2),matlu(:,4),natom) ! matlu4=m0m1
 call prod_matlu(matlu(:,2),green%moments(2)%matlu(:),matlu(:,5),natom) ! matlu5=m1m0

 call add_matlu(matlu(:,4),matlu(:,5),matlu(:,6),natom,1) ! matlu6=m0m1+m1m0

 if (optlog == 1) then
   call trace_matlu(matlu(:,4),natom,itau=0,trace=trace)
   green%trace_moments_log_loc(3) = trace + green%trace_moments_log_loc(3)
 end if ! optlog

 call add_matlu(matlu(:,1),matlu(:,6),matlu(:,5),natom,-1) ! matlu5=m2

 if (optlog == 1) then
   call trace_matlu(matlu(:,5),natom,itau=0,trace=trace)
   green%trace_moments_log_loc(3) = trace + green%trace_moments_log_loc(3)
 end if ! optlog

 if (option == 0) then
   call add_matlu(matlu(:,5),self%moments(3)%matlu(:),weiss%moments(3)%matlu(:),natom,-1)
 else if (option == 1) then
   call add_matlu(matlu(:,5),weiss%moments(3)%matlu(:),self%moments(3)%matlu(:),natom,-1)
 end if ! option

 call prod_matlu(green%moments(2)%matlu(:),matlu(:,3),matlu(:,1),natom) ! matlu=m0m0m0m0
 call prod_matlu(matlu(:,6),green%moments(2)%matlu(:),matlu(:,3),natom) ! matlu3 = m0m1m0+m1m0m0
 call prod_matlu(green%moments(2)%matlu(:),matlu(:,4),matlu(:,6),natom) ! matlu6 = m0m0m1
 call add_matlu(matlu(:,3),matlu(:,6),matlu(:,4),natom,1) ! matlu4 = m0m1m0+m1m0m0+m0m0m1

 if (optlog == 1) then
   call trace_matlu(matlu(:,4),natom,itau=0,trace=trace)
   green%trace_moments_log_loc(4) = trace * third
   call trace_matlu(matlu(:,1),natom,itau=0,trace=trace)
   green%trace_moments_log_loc(4) = green%trace_moments_log_loc(4) + trace*quarter
 end if ! optlog

 call add_matlu(matlu(:,1),matlu(:,4),matlu(:,3),natom,1) ! matlu3 = m0m1m0+m1m0m0+m0m0m1+m0m0m0m0

 call prod_matlu(green%moments(2)%matlu(:),matlu(:,5),matlu(:,1),natom) ! matlu=m0m2
 call prod_matlu(matlu(:,5),green%moments(2)%matlu(:),matlu(:,4),natom) ! matlu4=m2m0
 call add_matlu(matlu(:,1),matlu(:,4),matlu(:,5),natom,1) ! matlu5=m0m2+m2m0
 call prod_matlu(matlu(:,2),matlu(:,2),matlu(:,1),natom) ! matlu=m1m1
 call add_matlu(matlu(:,1),matlu(:,5),matlu(:,2),natom,1) ! matlu2=m1m1+m0m2+m2m0

 if (optlog == 1) then
   call trace_matlu(matlu(:,2),natom,itau=0,trace=trace)
   green%trace_moments_log_loc(4) = green%trace_moments_log_loc(4) + trace*half
 end if ! optlog

 call add_matlu(green%moments(5)%matlu(:),matlu(:,2),matlu(:,1),natom,-1)
 call add_matlu(matlu(:,1),matlu(:,3),matlu(:,2),natom,-1) ! matlu2=m3

 if (optlog == 1) then
   call trace_matlu(matlu(:,2),natom,itau=0,trace=trace)
   green%trace_moments_log_loc(4) = green%trace_moments_log_loc(4) + trace
 end if ! optlog

 if (option == 0) then
   call add_matlu(matlu(:,2),self%moments(4)%matlu(:),weiss%moments(4)%matlu(:),natom,-1)
 else if (option == 1) then
   call add_matlu(matlu(:,2),weiss%moments(4)%matlu(:),self%moments(4)%matlu(:),natom,-1)
 end if ! option

 do i=1,6
   call destroy_matlu(matlu(:,i),natom)
 end do

 ABI_FREE(lpawu)
 ABI_FREE(matlu)
 ABI_FREE(trace_loc)

end subroutine compute_moments_loc
!!***

!!****f* m_green/occup_fd
!! NAME
!! occup_fd
!!
!! FUNCTION
!!  Computes the Fermi-Dirac occupation
!!
!! INPUTS
!!  eig = eigenvalues
!!  fermie = Fermi level
!!  temp = temperature
!!
!! OUTPUT
!!
!!
!! SOURCE

 function occup_fd(eig,fermie,temp)

!Arguments ------------------------------------
!type
! Integrate analytic tail 1/(iw-mu)
 real(dp), intent(in) :: eig,fermie,temp
 real(dp) :: occup_fd
!Local variables-------------------------------
! *********************************************************************

 if ((eig-fermie) > zero) then
   occup_fd = exp(-(eig-fermie)/temp) / (one+exp(-(eig-fermie)/temp))
 else
   occup_fd = one / (one+exp((eig-fermie)/temp))
 end if ! eig-fermie>0

 end function occup_fd

END MODULE m_green
!!***
