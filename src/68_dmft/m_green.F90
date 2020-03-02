!!****m* ABINIT/m_green
!! NAME
!!  m_green
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2020 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

 MODULE m_green

 use defs_basis
 use m_xmpi
 use m_abicore
 use m_errors
 use m_lib_four
 use m_splines

 !use defs_abitypes, only : MPI_type
 use m_crystal, only : crystal_t
 use m_oper,     only : init_oper, destroy_oper, copy_oper, print_oper, inverse_oper, upfold_oper, &
                        loc_oper, trace_oper, oper_type
 use m_paw_dmft, only: paw_dmft_type, construct_nwli_dmft
 use m_matlu,    only : diff_matlu,print_matlu,trace_matlu, matlu_type, &
                        sym_matlu,zero_matlu,add_matlu,shift_matlu,init_matlu,destroy_matlu,copy_matlu
 use m_fstrings, only : int2char4
 use m_pawang,   only : pawang_type
 use m_self,     only : self_type
 use m_io_tools, only : open_file
 use m_time,     only : timab

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
 private :: occupfd
 private :: distrib_paral
 public :: greenldacompute_green
 public :: fermi_green
 public :: newton
 public :: local_ks_green
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

  integer :: dmft_nwlo
  ! dmft frequencies

  character(len=4) :: w_type
  ! type of frequencies used

  character(len=12) :: whichgreen
  ! describe the type of green function computed (LDA, DMFT, KS..)

  integer :: nw
  ! dmft frequencies

  integer :: fileprt_w
  ! 1 if file created

  integer :: fileprt_tau
  ! 1 if file created

  integer :: dmft_nwli
  ! dmft frequencies

  integer :: dmftqmc_l
  ! number of time slices for QMC

  integer :: ifermie_cv

  integer :: ichargeloc_cv

  integer :: use_oper_tau_ks
  ! 0 do not use oper_tau_ks
  ! 1 use oper_tau_ks

  real(dp) :: charge_ks
  ! Total charge computed from ks orbitals

  integer :: has_charge_matlu_solver
  ! =0 charge_matlu_solver not allocated
  ! =1 charge_matlu_solver is allocated
  ! =2 charge_matlu_solver is calculated (ie calculation of LOCAL CORRELATED occupations is done  from
  ! solver green function)

  integer :: has_greenmatlu_xsum
  ! =1 green%oper%matlu xsum in compute_green
  ! =0 green%oper%matlu non xsumed in compute_green
  ! used in integrate_green to checked that green function was computed in
  ! integrate_green.

  integer :: has_charge_matlu
  ! =2 if calculation of LOCAL CORRELATED occupations is done  from

  integer :: has_charge_matlu_prev
  ! =0 charge_matlu_prev not allocated
  ! =1 charge_matlu_prev is allocated
  ! =2 charge_matlu_prev is calculated (ie calculation of LOCAL CORRELATED occupations is done  from
  ! solver green function)

  integer, allocatable :: procb(:,:)

  integer, allocatable :: proct(:,:)

  real(dp), allocatable :: charge_matlu_prev(:,:)
  ! Total charge on correlated orbitals from previous iteration

  real(dp), allocatable :: charge_matlu(:,:)
  ! Total charge on correlated orbitals
! todo_ba name of charge_matlu is misleading: should be changed

  real(dp), allocatable :: charge_matlu_solver(:,:)
  ! Total charge on correlated orbitals obtained from solver by
  ! integration over frequencies.

  real(dp), allocatable :: tau(:)
  ! value of time in imaginary space

  real(dp), pointer :: omega(:) => null()
  ! value of frequencies

  real(dp), allocatable :: ecorr_qmc(:)
  ! Correlation energy for a given atom in qmc

  type(oper_type), allocatable :: oper(:)
  ! green function  in different basis

  type(oper_type), allocatable :: oper_tau(:)
  ! green function  in different basis

  type(oper_type) :: occup
  ! occupation in different basis

  type(oper_type) :: occup_tau
  ! occupation in different basis

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
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  opt_oper_ksloc (optional) = option for init_oper
!!  wtype = "real" Green function will be computed for real frequencies
!!        = "imag" Green function will be computed for imaginary frequencies
!!
!! OUTPUTS
!! green  = variable of type green_type
!!
!! PARENTS
!!      m_dmft,dyson,m_hubbard_one,m_green,qmc_prep_ctqmc,spectral_function
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

subroutine init_green(green,paw_dmft,opt_oper_ksloc,wtype)

!Arguments ------------------------------------
!scalars
!type
 type(green_type), intent(out) :: green
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer, optional, intent(in) :: opt_oper_ksloc
 character(len=4), optional :: wtype
!local variables ------------------------------------
 integer :: ifreq,nw,optoper_ksloc,nw_perproc
!************************************************************************
 if(present(opt_oper_ksloc)) then
   optoper_ksloc=opt_oper_ksloc
 else
   optoper_ksloc=2
 endif
 if(present(wtype)) then
   green%w_type=wtype
 else
   green%w_type="imag"
 endif


 if(green%w_type=="imag") then
   nw=paw_dmft%dmft_nwlo
   green%omega=>paw_dmft%omega_lo
 else if(green%w_type=="real") then
   nw=size(paw_dmft%omega_r)
   green%omega=>paw_dmft%omega_r
 endif

 green%dmft_nwlo=paw_dmft%dmft_nwlo
 green%dmft_nwli=paw_dmft%dmft_nwli
 green%charge_ks=zero
 green%has_charge_matlu=0
 ABI_ALLOCATE(green%charge_matlu,(paw_dmft%natom,paw_dmft%nsppol+1))
 green%charge_matlu=zero
 green%has_charge_matlu=1
 green%has_greenmatlu_xsum=0

 green%has_charge_matlu_solver=0
 ABI_ALLOCATE(green%charge_matlu_solver,(paw_dmft%natom,paw_dmft%nsppol+1))
 green%charge_matlu_solver=zero
 green%has_charge_matlu_solver=1

 green%has_charge_matlu_prev=0
 ABI_ALLOCATE(green%charge_matlu_prev,(paw_dmft%natom,paw_dmft%nsppol+1))
 green%charge_matlu_prev=zero
 green%has_charge_matlu_prev=1

 call init_oper(paw_dmft,green%occup,opt_ksloc=3)

!  built simple arrays to distribute the tasks in compute_green.
 ABI_ALLOCATE(green%procb,(nw,paw_dmft%nkpt))
 ABI_ALLOCATE(green%proct,(nw,0:paw_dmft%nproc-1))

 call distrib_paral(paw_dmft%nkpt,paw_dmft%nproc,nw,nw_perproc,green%procb,green%proct)
 green%nw=nw

!  need to distribute memory over frequencies

!!  begin of temporary modificatios
! ABI_DATATYPE_ALLOCATE(green%oper,(green%nw_perproc))
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
 ABI_DATATYPE_ALLOCATE(green%oper,(green%nw))
 do ifreq=1,green%nw
  call init_oper(paw_dmft,green%oper(ifreq),opt_ksloc=optoper_ksloc)
 enddo
 green%ifermie_cv=0
 green%ichargeloc_cv=0

 if(paw_dmft%dmft_solv>=4) then
   ABI_ALLOCATE(green%ecorr_qmc,(paw_dmft%natom))
 end if
 if(paw_dmft%dmft_solv>=4) green%ecorr_qmc(:)=zero


 green%fileprt_w=0
 green%fileprt_tau=0


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
!! PARENTS
!!      impurity_solve,m_green
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

subroutine init_green_tau(green,paw_dmft,opt_ksloc)

!Arguments ------------------------------------
!scalars
!type
 type(green_type), intent(inout) :: green
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer, optional, intent(in) :: opt_ksloc
!local variables ------------------------------------
 integer :: optksloc
 integer :: itau
!************************************************************************
 if(present(opt_ksloc)) then
   optksloc=opt_ksloc
 else
   optksloc=3
 endif
 green%use_oper_tau_ks=0
 if(green%use_oper_tau_ks==0) then
   optksloc=2
 endif

 green%dmftqmc_l=paw_dmft%dmftqmc_l
 ABI_ALLOCATE(green%tau,(green%dmftqmc_l))
 do itau=1,green%dmftqmc_l
  green%tau(itau)=float(itau-1)/float(green%dmftqmc_l)/paw_dmft%temp
 enddo

 call init_oper(paw_dmft,green%occup_tau,optksloc)

 ABI_DATATYPE_ALLOCATE(green%oper_tau,(paw_dmft%dmftqmc_l))
 do itau=1,green%dmftqmc_l
  call init_oper(paw_dmft,green%oper_tau(itau),optksloc)
 enddo

end subroutine init_green_tau
!!***

!!****f* m_green/destroy_green
!! NAME
!! destroy_green
!!
!! FUNCTION
!!  deallocate green
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!
!! OUTPUT
!!
!! PARENTS
!!      m_dmft,dyson,m_hubbard_one,m_green,qmc_prep_ctqmc,spectral_function
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

subroutine destroy_green(green)

!Arguments ------------------------------------
!scalars
 type(green_type),intent(inout) :: green

!local variables-------------------------------
 integer :: ifreq

! *********************************************************************
 call destroy_oper(green%occup)
 if ( allocated(green%oper))       then
   do ifreq=1,green%nw
    call destroy_oper(green%oper(ifreq))
   enddo
   ABI_DATATYPE_DEALLOCATE(green%oper)
 end if
 if ( allocated(green%charge_matlu))       then
   ABI_DEALLOCATE(green%charge_matlu)
 end if
 green%has_charge_matlu=0

 if ( allocated(green%charge_matlu_prev))       then
   ABI_DEALLOCATE(green%charge_matlu_prev)
 end if
 green%has_charge_matlu_prev=0

 if ( allocated(green%charge_matlu_solver))       then
   ABI_DEALLOCATE(green%charge_matlu_solver)
 end if
 green%has_charge_matlu_solver=0

 if ( allocated(green%ecorr_qmc))   then
    ABI_DEALLOCATE(green%ecorr_qmc)
 end if
 if ( allocated(green%procb))   then
    ABI_DEALLOCATE(green%procb)
 end if
 if ( allocated(green%proct))   then
    ABI_DEALLOCATE(green%proct)
 end if
 green%omega => null()

end subroutine destroy_green
!!***

!!****f* m_green/destroy_green_tau
!! NAME
!! destroy_green_tau
!!
!! FUNCTION
!!  deallocate green
!!
!! INPUTS
!!  green  <type(green_type)>= green function data
!!
!! OUTPUT
!!
!! PARENTS
!!      impurity_solve,m_green
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

subroutine destroy_green_tau(green)

!Arguments ------------------------------------
!scalars
 type(green_type),intent(inout) :: green
! integer, optional, intent(in) :: opt_ksloc
!local variables-------------------------------
 integer :: itau
! integer :: optksloc
! *********************************************************************
! if(present(opt_ksloc)) then
!   optksloc=opt_ksloc
! else
!   optksloc=3
! endif

 call destroy_oper(green%occup_tau)
 if ( allocated(green%oper_tau))       then
   do itau=1,green%dmftqmc_l
    call destroy_oper(green%oper_tau(itau))
   enddo
   ABI_DATATYPE_DEALLOCATE(green%oper_tau)
 end if
 if ( allocated(green%tau))            then
   ABI_DEALLOCATE(green%tau)
 end if

end subroutine destroy_green_tau
!!***

!!****f* m_green/copy_green
!! NAME
!! copy_green
!!
!! FUNCTION
!!  copy one data structure green1 into green2
!!
!! INPUTS
!!  green1  <type(green_type)>= green function data
!!  green2  <type(green_type)>= green function data
!!  opt_tw = option to precise which data to copy
!!          1: copy only green%occup_tau and green%oper_tau data
!!          2: copy only green%occup and green%oper  data (frequency)
!!
!! OUTPUT
!!
!! PARENTS
!!      dyson,impurity_solve,m_green,qmc_prep_ctqmc,spectral_function
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

subroutine copy_green(green1,green2,opt_tw)

!Arguments ------------------------------------
!type
 type(green_type),intent(in) :: green1
 type(green_type),intent(inout) :: green2
 integer, intent(in) :: opt_tw

!local variables-------------------------------
 integer :: ifreq,itau
! *********************************************************************

 if(opt_tw==2) then
   call copy_oper(green1%occup,green2%occup)
   do ifreq=1, green1%nw
     call copy_oper(green1%oper(ifreq),green2%oper(ifreq))
     if(green1%has_greenmatlu_xsum==1) then ! indicate to integrate_green  that xsum has been done
!  for matlu in compute_green.
        green2%has_greenmatlu_xsum=1
     endif
   enddo
 else if (opt_tw==1) then
   call copy_oper(green1%occup_tau,green2%occup_tau)
   do itau=1,green1%dmftqmc_l
     call copy_oper(green1%oper_tau(itau),green2%oper_tau(itau))
   enddo
 endif

end subroutine copy_green
!!***

!!****f* m_green/printocc_green
!! NAME
!! printocc_green
!!
!! FUNCTION
!!  print occupations
!!
!! INPUTS
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  green  <type(green_type)>= green function data
!! option= 1 :for G(w)
!!         2 :for G(tau)
!!         3 :for G(tau) and check % G(w)
!!         4
!!         <5: write diagonal part of KS occupation matrix
!!         5: for G(w)
!!         6: for G(tau)
!!         7 :for G(tau) and check % G(w)
!!         >8: write all elements of KS occup. matrix.
!!         9: for G(w)
!!
!! OUTPUT
!!
!! PARENTS
!!      m_dmft,impurity_solve,m_green,qmc_prep_ctqmc
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

subroutine printocc_green(green,option,paw_dmft,pawprtvol,opt_weissgreen,chtype)

!Arguments ------------------------------------
!type
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(green_type),intent(in) :: green
 integer,intent(in) :: option,pawprtvol
 integer,optional,intent(in) :: opt_weissgreen
 character(len=*), optional, intent(in) :: chtype

!local variables-------------------------------
 character(len=500) :: message
 integer :: optweissgreen,i_tau
! *********************************************************************
 if(present(opt_weissgreen)) then
   optweissgreen=opt_weissgreen
 else
   optweissgreen=2
 endif


 if(mod(option,4)==1) then
   if(optweissgreen==2) then
     if(present(chtype)) then
       write(message,'(4a)') ch10,"  == The ",trim(chtype)," occupations are  == "
     else
       write(message,'(2a)') ch10,"  == The occupations (integral of the Green function) are  == "
     endif
   else if(optweissgreen==1) then
     write(message,'(2a)') ch10,"  == The integrals of the Weiss function are  == "
   endif
   call wrtout(std_out,message,'COLL')
   call print_oper(green%occup,option,paw_dmft,pawprtvol)
 endif
 if(mod(option,4)>=2) then
   if(optweissgreen==2) then
     write(message,'(2a)') ch10,"  == The occupations (value of G(tau) for tau=0-) are  == "
   else if(optweissgreen==1) then
     write(message,'(2a)') ch10,"  == Values of G_0(tau) for tau=0- are  == "
   endif
   call wrtout(std_out,message,'COLL')
   call print_oper(green%occup_tau,option,paw_dmft,pawprtvol)
!   write(message,'(2a)') ch10," == check: occupations from Green functions are  == "
!   call wrtout(std_out,message,'COLL')
!   call print_oper(green%occup,1,paw_dmft,pawprtvol)
   if(mod(option,4)>=3) then
     call diff_matlu("Local occup from integral of G(w) ",&
&       "Local occup from G(tau=0-) ",&
&       green%occup%matlu,green%occup_tau%matlu,paw_dmft%natom,1,tol4)
     write(message,'(2a)') ch10,&
&     '  *****  => Calculations of occupations in omega and tau spaces are coherent ****'
     call wrtout(std_out,message,'COLL')
   endif
 endif

 if(present(chtype)) then
   if(paw_dmft%prtvol>=4.and.&
&      (chtype=="DMFT (end of DMFT loop)".or.chtype=="converged DMFT")&
&      .and.green%occup%has_opermatlu==1) then
     write(message,'(4a)') ch10,"  == The DFT+DMFT occupation matrix for correlated electrons is == "
     call wrtout(ab_out,message,'COLL')
     call print_matlu(green%occup%matlu,paw_dmft%natom,pawprtvol,opt_ab_out=1)
     write(message,'(a)') "  "
     call wrtout(ab_out,message,'COLL')
   endif
 endif

 if(mod(option,4)>=2) then
   i_tau=1
   if(optweissgreen==1) i_tau=-1
   call trace_matlu(green%occup_tau%matlu,paw_dmft%natom,itau=i_tau)
 endif

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
!!
!! OUTPUT
!!
!! PARENTS
!!      impurity_solve,m_green,qmc_prep_ctqmc,spectral_function
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

subroutine print_green(char1,green,option,paw_dmft,pawprtvol,opt_wt,opt_decim)

!Arguments ------------------------------------
!type
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(green_type),intent(inout) :: green
 integer,intent(in) :: option,pawprtvol
 integer,intent(in),optional :: opt_wt,opt_decim
 character(len=*), intent(in) :: char1

!local variables-------------------------------
 integer :: iall,iatom,ib,ifreq,ikpt,im,ispinor,isppol,itau
 integer :: lsub,mbandc,natom,ndim,nkpt,nspinor,nsppol,optwt,spf_unt,spfkresolved_unt,spcorb_unt
 character(len=2000) :: message
 integer,allocatable :: unitgreenfunc_arr(:)
 integer,allocatable :: unitgreenloc_arr(:)
 character(len=fnlen) :: tmpfil
 character(len=1) :: tag_is,tag_is2
 character(len=10) :: tag_at
 character(len=3) :: tag_ik
 real(dp) :: re, ima
 complex(dpc), allocatable :: sf(:),sf_corr(:)
! *********************************************************************
 if(present(opt_wt)) then
   optwt=opt_wt
 else
   optwt=1
 endif


 if(pawprtvol>200) then
 endif
 natom=green%oper(1)%natom
 nsppol=paw_dmft%nsppol
 nspinor=paw_dmft%nspinor
 mbandc=paw_dmft%mbandc
 nkpt=paw_dmft%nkpt

!  == Print local Green Function
 if(option==1.or.option==3) then
   ABI_ALLOCATE(unitgreenfunc_arr,(natom*nsppol*nspinor))
   iall=0
   do iatom=1,natom
     if(green%oper(1)%matlu(iatom)%lpawu.ne.-1) then
       ndim=2*green%oper(1)%matlu(iatom)%lpawu+1
       call int2char4(iatom,tag_at)
       ABI_CHECK((tag_at(1:1)/='#'),'Bug: string length too short!')
       do isppol=1,nsppol
         write(tag_is,'(i1)')isppol
         do ispinor=1,nspinor
           iall=iall+1
           write(tag_is2,'(i1)')ispinor

!         == Create names
           if(optwt==1) then
             tmpfil = &
&             trim(paw_dmft%filapp)//'Green-'//trim(char1)//'-omega_iatom'// &
&             trim(tag_at)//'_isppol'//tag_is//'_ispinor'//tag_is2
           else
             tmpfil = &
&             trim(paw_dmft%filapp)//'Green-'//trim(char1)//'-tau_iatom'// &
&             trim(tag_at)//'_isppol'//tag_is//'_ispinor'//tag_is2
           endif
           if(iall<=4) then
             write(message,'(3a)') ch10,"  == Print green function on file ",tmpfil
             call wrtout(std_out,message,'COLL')
           elseif(iall==5)  then
             write(message,'(3a)') ch10,"  == following values are printed in files"
             call wrtout(std_out,message,'COLL')
           endif
           unitgreenfunc_arr(iall)=300+iall-1
           if(optwt==1.or.green%fileprt_tau==0) then
             open(unit=unitgreenfunc_arr(iall),file=trim(tmpfil),status='unknown',form='formatted',position='append')
           else if(optwt==2.and.green%fileprt_tau==1) then
             open(unit=unitgreenfunc_arr(iall),file=trim(tmpfil),status='unknown',form='formatted',position='append')
           endif

!         == Write in files
!           rewind(unitgreenfunc_arr(iall))
!           write(message,'(a,a,a,i4)') 'opened file : ', trim(tmpfil), ' unit', unitgreenfunc_arr(iall)
!           call wrtout(std_out,message,'COLL')
           write(message,'(a,a)') ch10,"# New record :"
           call wrtout(unitgreenfunc_arr(iall),message,'COLL')
           if(optwt==1) then
             do ifreq=1,green%nw
               if(present(opt_decim)) then
                 write(message,'(2x,30(e23.16,2x))') &
&                green%omega(ifreq),&
&                (green%oper(ifreq)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor),im=1,ndim)
               else
                 write(message,'(2x,30(e10.3,2x))') &
&                green%omega(ifreq),&
&                (green%oper(ifreq)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor),im=1,ndim)
               endif
               call wrtout(unitgreenfunc_arr(iall),message,'COLL')
               re=real(green%oper(ifreq)%matlu(iatom)%mat(1,1,isppol,ispinor,ispinor))
               ima=aimag(green%oper(ifreq)%matlu(iatom)%mat(1,1,isppol,ispinor,ispinor))
!               write(228,*) green%omega(ifreq),re/(re**2+ima**2),ima/(re**2+ima**2)+green%omega(ifreq)
             enddo
           else
             do itau=1,green%dmftqmc_l
                 write(message,'(2x,30(e10.3,2x))') &
&                green%tau(itau),&
&               (green%oper_tau(itau)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor),im=1,ndim)
               call wrtout(unitgreenfunc_arr(iall),message,'COLL')
             enddo
             write(message,'(2x,30(e10.3,2x))') &
&            one/paw_dmft%temp,&
&            (-green%oper_tau(1)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)-one,im=1,ndim)
             call wrtout(unitgreenfunc_arr(iall),message,'COLL')
           endif
         enddo
       enddo ! isppol
     endif ! lpawu=/-1
   enddo ! iatom
   ABI_DEALLOCATE(unitgreenfunc_arr)
 endif

!  == Print ks green function
 if((option==2.or.option==3).and.green%oper(1)%has_operks==1) then
   ABI_ALLOCATE(unitgreenloc_arr,(nsppol*nkpt))
   iall=0
   do isppol = 1 , nsppol
     write(tag_is,'(i1)')isppol
     do ikpt = 1, nkpt
       write(tag_ik,'(i3)')ikpt
!      do ib1 = 1, mbandc
       iall=iall+1

!         == Create names
       if(optwt==1) then
         tmpfil = trim(paw_dmft%filapp)//'Green-'//trim(char1)//'-omega_isppol'//tag_is//'_ikpt'//trim(adjustl(tag_ik))
       else
         tmpfil = trim(paw_dmft%filapp)//'Green-'//trim(char1)//'-tau_isppol'//tag_is//'_ikpt'//trim(adjustl(tag_ik))
       endif
       if(iall<=4)  then
         write(message,'(3a)') ch10,"  == Print green function on file ",tmpfil
         call wrtout(std_out,message,'COLL')
       elseif(iall==5)  then
         write(message,'(3a)') ch10,"  == following values are printed in files"
         call wrtout(std_out,message,'COLL')
       endif
       unitgreenloc_arr(iall)=400+iall-1
       open (unit=unitgreenloc_arr(iall),file=trim(tmpfil),status='unknown',form='formatted')
!      rewind(unitgreenloc_arr(iall))
!       write(message,'(a,a,a,i4)') 'opened file : ', trim(tmpfil), ' unit', unitgreenloc_arr(iall)
!       call wrtout(std_out,message,'COLL')

!         == Write in files
       write(message,'(a,a)') ch10,"# New record : First 20 bands"
       call wrtout(unitgreenloc_arr(iall),message,'COLL')
!       call flush(std_out)
       do lsub=1,mbandc/20+1
         if(optwt==1) then
           do ifreq=1,green%nw
!             call flush(std_out)
             write(message,'(2x,50(e10.3,2x))') &
&            green%omega(ifreq), &
&            (green%oper(ifreq)%ks(isppol,ikpt,ib,ib),ib=20*(lsub-1)+1,min(20*lsub,mbandc))
             call wrtout(unitgreenloc_arr(iall),message,'COLL')
               re=real(green%oper(ifreq)%ks(isppol,ikpt,1,1))
               ima=aimag(green%oper(ifreq)%ks(isppol,ikpt,1,1))
               if(ikpt==1) write(229,*) green%omega(ifreq),re/(re**2+ima**2),ima/(re**2+ima**2)+green%omega(ifreq)
           enddo
         else
           if(green%use_oper_tau_ks==1) then
             do itau=1,green%dmftqmc_l
!               call flush(std_out)
               write(message,'(2x,50(e10.3,2x))') &
&              green%tau(itau), &
&              (green%oper_tau(itau)%ks(isppol,ikpt,ib,ib),ib=20*(lsub-1)+1,min(20*lsub,mbandc))
               call wrtout(unitgreenloc_arr(iall),message,'COLL')
             enddo
           endif
         endif
         if(20*lsub<mbandc) write(message,'(a,a,i5,a,i5)')    &
&           ch10,"# Same record, Following bands : From ",    &
&           20*(lsub),"  to ",min(20*(lsub+1),mbandc)
         call wrtout(unitgreenloc_arr(iall),message,'COLL')
       enddo
     enddo ! ikpt
   enddo ! isppol
   ABI_DEALLOCATE(unitgreenloc_arr)
 endif

 if((green%w_type=="real".and.option>=4).and.green%oper(1)%has_operks==1) then
   write(message,'(a,a)') ch10,"  == About to print spectral function"
   call wrtout(std_out,message,'COLL')
   if (option==4) then
     tmpfil = trim(paw_dmft%filapp)//'SpFunc-'//trim(char1)
     if (open_file(tmpfil, message, newunit=spf_unt, status='unknown', form='formatted') /= 0) then
       MSG_ERROR(message)
     end if
   endif
   if (option==5) then
     tmpfil = trim(paw_dmft%filapp)//'_DFTDMFT_SpectralFunction_kresolved_'//trim(char1)
     if (open_file(tmpfil, message, newunit=spfkresolved_unt, status='unknown', form='formatted') /= 0) then
       MSG_ERROR(message)
     end if
   endif
   ABI_ALLOCATE(sf,(green%nw))
   ABI_ALLOCATE(sf_corr,(green%nw))
   iall=0
   if (option==5) then
       do ikpt = 1, nkpt
         sf=czero
         do ifreq=1,green%nw
           do isppol = 1 , nsppol
             do ib=1,mbandc
               sf(ifreq)=sf(ifreq)+green%oper(ifreq)%ks(isppol,ikpt,ib,ib)
             enddo
           enddo
           write(message,*) green%omega(ifreq)*Ha_eV,(-aimag(sf(ifreq)))/pi/Ha_eV,ikpt
           call wrtout(spfkresolved_unt,message,'COLL')
         enddo
           write(message,*)
           call wrtout(spfkresolved_unt,message,'COLL')
       enddo
       write(message,*) ch10
       call wrtout(spfkresolved_unt,message,'COLL')
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
   endif
   if (option==4) then
         sf=czero
     do isppol = 1 , nsppol
       do ikpt = 1, nkpt
         do ib=1,mbandc
           do ifreq=1,green%nw
             sf(ifreq)=sf(ifreq)+green%oper(ifreq)%ks(isppol,ikpt,ib,ib)*green%oper(1)%wtk(ikpt)
           enddo
         enddo
       enddo
     enddo
   endif
   if(paw_dmft%dmft_kspectralfunc==1) then
     do iatom=1,natom
       if(green%oper(1)%matlu(iatom)%lpawu.ne.-1) then
         sf_corr=czero
         call int2char4(iatom,tag_at)
         ABI_CHECK((tag_at(1:1)/='#'),'Bug: string length too short!')
         tmpfil = trim(paw_dmft%filapp)//'_DFTDMFT_spectralfunction_orb_'//trim(char1)//'_iatom'//trim(tag_at)
         if (open_file(tmpfil, message, newunit=spcorb_unt, status='unknown', form='formatted') /= 0) then
           MSG_ERROR(message)
         end if
         write(message,*) "#", nspinor,nsppol,ndim,green%nw
         call wrtout(spcorb_unt,message,'COLL')
         write(message,*) "#", green%oper(1)%matlu(iatom)%lpawu
         call wrtout(spcorb_unt,message,'COLL')
         ndim=2*green%oper(1)%matlu(iatom)%lpawu+1
         do isppol = 1 , nsppol
           do ispinor=1, nspinor
             do im=1,ndim
               do ifreq=1,green%nw
                 sf_corr(ifreq)=sf_corr(ifreq)+ green%oper(ifreq)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)
               enddo
             enddo
           enddo
         enddo
         do ifreq=1,green%nw
           write(message,*) green%omega(ifreq),(-aimag(sf_corr(ifreq)))/3.141592653589793238_dp
           call wrtout(spcorb_unt,message,'COLL')
         enddo
         close(spcorb_unt)
       endif
     enddo
   endif
   if (option==4) then
     do ifreq=1,green%nw
       write(message,*) green%omega(ifreq),(-aimag(sf(ifreq)))/3.141592653589793238_dp
       call wrtout(spf_unt,message,'COLL')
     enddo
   endif
   close(spf_unt)
   ABI_DEALLOCATE(sf)
   ABI_DEALLOCATE(sf_corr)
 endif

 if(optwt==2.and.(option==1.or.option==3)) then
    green%fileprt_tau=1  ! file for G(tau) has been created here
 endif

end subroutine print_green
!!***

!!****f* m_green/compute_green
!! NAME
!! compute_green
!!
!! FUNCTION
!! compute green function from LDA and self-energy
!!
!! INPUTS
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  opt_self = optional argument, if =1, upfold self-energy
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  green  <type(green_type)>= green function data
!!  opt_nonxsum = 0 : do usual xsum after calculation of green(freq)%ks
!!              = 1 : do not do xsum after calculation of green(freq)%ks: each proc as only a part
!!                    of the data: this is useful where only the total number of electron will be computed.
!!  opt_nonxsum2 0 : green(ifreq)%matlu will be broadcasted
!!               1 : green(ifreq)%matlu will not be broadcasted in compute_green: calc
!!                   if occupations will not possible.
!!                   (a keyword: compute_local_green would be in fact equivalent and more clear)
!!
!! OUTPUT
!!
!! PARENTS
!!      m_dmft,fermi_green,m_green,newton,spectral_function
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

subroutine compute_green(cryst_struc,green,paw_dmft,pawang,prtopt,self,opt_self,&
&           opt_nonxsum,opt_nonxsum2)

!Arguments ------------------------------------
!type
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type),intent(inout) :: green
 type(paw_dmft_type), intent(in) :: paw_dmft
 !type(MPI_type), intent(in) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(self_type), intent(inout) :: self
 integer, intent(in) :: prtopt
 integer, optional, intent(in) :: opt_self
 integer, optional, intent(in) :: opt_nonxsum
 integer, optional, intent(in) :: opt_nonxsum2

!local variables-------------------------------
 !logical :: lintegrate
 integer :: iatom,ib,ib1,ierr,ifreq,ikpt,is
 integer :: icomp_chloc
 integer :: natom,mband,mbandc,nkpt
 integer :: myproc,nproc,spacecomm,nspinor,nsppol,option,optself,optnonxsum
 integer :: optnonxsum2
 integer, allocatable :: procb_ifreq(:)
 character(len=500) :: message
 real(dp) :: fermilevel
 real(dp), allocatable :: Id(:,:)
! integer, allocatable :: procb(:,:),proct(:,:)
 real(dp) :: tsec(2)
 type(oper_type) :: self_minus_hdc_oper
 complex(dpc) :: omega_current
 type(oper_type) :: green_temp
! *********************************************************************

 !lintegrate=.true.
 !if(lintegrate.and.green%w_type=="real") then
 !if(green%w_type=="real") then
 !  message = 'integrate_green not implemented for real frequency'
 !  MSG_BUG(message)
 !endif
 call timab(624,1,tsec)
 if(present(opt_self)) then
   optself=opt_self
 else
   optself=0
 endif
 if(present(opt_nonxsum)) then
   optnonxsum=opt_nonxsum
 else
   optnonxsum=0
 endif
 if(present(opt_nonxsum2)) then
   optnonxsum2=opt_nonxsum2
 else
   optnonxsum2=0
 endif

 if(prtopt>0)  then
   write(message,'(2a,i3,13x,a)') ch10,'  ===  Compute green function '
   call wrtout(std_out,message,'COLL')
 endif

 if(self%nw/=green%nw)  then
   message = ' BUG: frequencies for green and self not coherent'
   MSG_BUG(message)
 endif

! Initialise spaceComm, myproc, and nproc
 spacecomm=paw_dmft%spacecomm
 myproc=paw_dmft%myproc
 nproc=paw_dmft%nproc

! Initialise integers
 mband   = paw_dmft%mband
 mbandc  = paw_dmft%mbandc
 nkpt    = paw_dmft%nkpt
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol
 natom   = paw_dmft%natom

! Initialise for compiler
 omega_current=czero
 icomp_chloc=0

! ==============================
! Initialise Identity
! ==============================
 ABI_ALLOCATE(Id,(paw_dmft%mbandc,paw_dmft%mbandc))
 Id=zero
 do ib=1, paw_dmft%mbandc
   Id(ib,ib)=one
 enddo ! ib

 if(prtopt/=0.and.prtopt>-100)  then
   write(message,'(2a)') ch10,'  == Green function is computed:'
   call wrtout(std_out,message,'COLL')
 endif
 option=1
 fermilevel=paw_dmft%fermie
 if(option==123) then
   fermilevel=2.d0
   write(message,'(2a,e14.3,a)') ch10,&
&  '  Warning (special case for check: fermi level=',fermilevel,')'
   call wrtout(std_out,message,'COLL')
 endif

! ====================================================
! Upfold self-energy and double counting  Self_imp -> self(k)
! ====================================================
! if(optself==1) then
!   do ifreq=1,green%nw
!     call upfold_oper(self%oper(ifreq),paw_dmft,1)
!   enddo ! ifreq
!   call upfold_oper(self%hdc,paw_dmft,1)
! endif

! =================================================================
! Initialize green%oper before calculation (important for xmpi_sum)
! =================================================================
 do ifreq=1,green%nw
   call zero_matlu(green%oper(ifreq)%matlu,green%oper(ifreq)%natom)
 enddo ! ifreq
 !call xmpi_barrier(spacecomm)

! ================================
! == Compute Green function G(k)
! ================================
 green%occup%ks=czero
 ABI_ALLOCATE(procb_ifreq,(paw_dmft%nkpt))
 do ifreq=1,green%nw
   !if(present(iii)) write(6,*) ch10,'ifreq  self', ifreq,self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)
!   ====================================================
!   First Upfold self-energy and double counting  Self_imp -> self(k)
!   ====================================================
!   if(mod(ifreq-1,nproc)==myproc) then
!    write(6,*) "compute_green ifreq",ifreq, mod(ifreq-1,nproc)==myproc,proct(ifreq,myproc)==1
   if(green%proct(ifreq,myproc)==1) then
     if(green%w_type=="imag") then
       omega_current=cmplx(zero,green%omega(ifreq),kind=dp)
     else if(green%w_type=="real") then
       omega_current=cmplx(green%omega(ifreq),paw_dmft%temp,kind=dp)
     endif
     call init_oper(paw_dmft,self_minus_hdc_oper)
     call init_oper(paw_dmft,green_temp)

     call add_matlu(self%oper(ifreq)%matlu,self%hdc%matlu,&
&                   self_minus_hdc_oper%matlu,green%oper(ifreq)%natom,-1)
! do iatom = 1 , natom
   !write(6,*) 'self matlu', ifreq, self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)
   !write(6,*) 'self hdc  ', ifreq, self%hdc%matlu(1)%mat(1,1,1,1,1)
   !write(6,*) 'self_minus_hdc_oper  ', ifreq, self_minus_hdc_oper%matlu(1)%mat(1,1,1,1,1)
! enddo ! natom
     if(paw_dmft%dmft_solv==4)  then
       call shift_matlu(self_minus_hdc_oper%matlu,paw_dmft%natom,cmplx(self%qmc_shift,0.d0,kind=dp),-1)
       call shift_matlu(self_minus_hdc_oper%matlu,paw_dmft%natom,cmplx(self%qmc_xmu,0.d0,kind=dp),-1)
     endif
     if(optself==0) then
       call zero_matlu(self_minus_hdc_oper%matlu,paw_dmft%natom)
     end if

     procb_ifreq=green%procb(ifreq,:)
     call upfold_oper(self_minus_hdc_oper,paw_dmft,1,procb=procb_ifreq,iproc=myproc,prt=1)
     do ib1 = 1 , paw_dmft%mbandc
       do ib = 1 , paw_dmft%mbandc
         do ikpt = 1 , paw_dmft%nkpt
           do is = 1 , paw_dmft%nsppol
             if (green%procb(ifreq,ikpt)==myproc) then
!               green%oper(ifreq)%ks(is,ikpt,ib,ib1)=       &
               green_temp%ks(is,ikpt,ib,ib1)=       &
&               ( omega_current     &
&               + fermilevel                               &
&               - paw_dmft%eigen_lda(is,ikpt,ib)) * Id(ib,ib1) &
&               - self_minus_hdc_oper%ks(is,ikpt,ib,ib1)
          !if(ikpt==2.and.ib==ib1) then
          !  write(6,*) "self",ib1,ib,ikpt,is,ifreq,self_minus_hdc_oper%ks(is,ikpt,ib,ib1)
          !endif
!&               - (self%oper(ifreq)%ks(is,ikpt,ib,ib1)-self%hdc%ks(is,ikpt,ib,ib1))
!               if(prtopt>5) then
!               if(ikpt==1.and.(ifreq==1.or.ifreq==3).and.ib==16.and.ib1==16) then
!                write(std_out,*) 'omega_current                         ',omega_current
!                write(std_out,*) 'fermilevel                            ',fermilevel
!                write(std_out,*) ' paw_dmft%eigen_lda(is,ikpt,ib)       ', paw_dmft%eigen_lda(is,ikpt,ib),Id(ib,ib1)
!                write(std_out,*) 'self_minus_hdc_oper%ks(is,ikpt,ib,ib1)',self_minus_hdc_oper%ks(is,ikpt,ib,ib1)
!                write(std_out,*) 'green                                 ',green%oper(ifreq)%ks(is,ikpt,ib,ib1)
!               endif
!               if(ib==1.and.ib1==3) then
!                 write(std_out,*) "ff compute",ikpt,ifreq,is,ikpt,ib,ib1
!                 write(std_out,*) "ff compute",ikpt,ifreq, green_temp%ks(is,ikpt,ib,ib1)
!                 write(std_out,*) "ff details",paw_dmft%eigen_lda(is,ikpt,ib)
!                 write(std_out,*) "ff details2",fermilevel
!                 write(std_out,*) "ff details3",Id(ib,ib1)
!        !        write(std_out,*) "ff details4",self_minus_hdc_oper%ks(is,ikpt,ib,ib1)
!               endif
             endif
           enddo ! is
         enddo ! ikpt
       enddo ! ib
     enddo ! ib1
!     call print_oper(green%oper(ifreq),9,paw_dmft,3)
!       write(std_out,*) 'after print_oper'
!     if(ifreq==1.or.ifreq==3) then
!         write(std_out,*) 'after print_oper', ifreq
!       write(std_out,*) 'green1  ifreq         %ks(1,1,16,16)',ifreq,green%oper(ifreq)%ks(1,1,16,16)
!     endif
!       write(std_out,*) 'before inverse_oper'
    call inverse_oper(green_temp,2,prtopt,procb=procb_ifreq,iproc=myproc)
     if (green%oper(1)%has_operks==1) green%oper(ifreq)%ks=green_temp%ks
    !if(ifreq==1) then
    !  write(std_out,*) "1188",green_temp%ks(1,1,8,8)
    !  write(std_out,*) "1189",green_temp%ks(1,1,8,9)
    !  write(std_out,*) "1198",green_temp%ks(1,1,9,8)
    !  write(std_out,*) "1199",green_temp%ks(1,1,9,9)
    !endif

     !if(lintegrate) then
!    accumulate integration
     if(green%w_type/="real") then
       do is = 1 , paw_dmft%nsppol
         do ib = 1 , paw_dmft%mbandc
           do ib1 = 1 , paw_dmft%mbandc
             do ikpt = 1 , paw_dmft%nkpt
               if (green%procb(ifreq,ikpt)==myproc) then
                 call add_int_fct(ifreq,green_temp%ks(is,ikpt,ib,ib1),ib==ib1,    &
&                      omega_current,2,green%occup%ks(is,ikpt,ib,ib1),            &
&                      paw_dmft%temp,paw_dmft%wgt_wlo(ifreq),paw_dmft%dmft_nwlo)
               endif
             enddo ! ikpt
           enddo ! ib1
         enddo ! ib
       enddo ! is
     endif
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
     call loc_oper(green_temp,paw_dmft,1,procb=procb_ifreq,iproc=myproc)
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
     call sym_matlu(cryst_struc,green_temp%matlu,pawang,paw_dmft)
     call copy_matlu(green_temp%matlu,green%oper(ifreq)%matlu,natom)
!     if(ifreq==1.and.ifreq==11) then
!       write(std_out,*) ifreq,nproc,'after sym'
!       call print_matlu(green%oper(ifreq)%matlu,green%oper(ifreq)%natom,2,-1,0)
!       ! ok
!     endif
     call destroy_oper(self_minus_hdc_oper)
     call destroy_oper(green_temp)
! call flush(std_out)
   endif ! parallelisation
 enddo ! ifreq
 ABI_DEALLOCATE(procb_ifreq)

! =============================================
! built total green function (sum over procs).
! =============================================
 !call xmpi_barrier(spacecomm)
 do ifreq=1,green%nw
!  call xmpi_sum(green%oper(ifreq)%ks,spacecomm,ierr)
!       print *, "myproc, proct, ifreq ------------------- ", myproc, ifreq
!   do ikpt=1,paw_dmft%nkpt
!     call xmpi_bcast(green%oper(ifreq)%ks(:,ikpt,:,:),procb(ifreq,ikpt),spacecomm,ierr)
!   enddo
!! or
   if(optnonxsum==0.and.green%oper(ifreq)%has_operks==1) then
     call xmpi_sum(green%oper(ifreq)%ks,spacecomm,ierr)
   else if(optnonxsum==0.and.green%oper(ifreq)%has_operks==0) then
     message = 'optnonxsum==0 and green%oper(ifreq)%has_operks==0: not compatible'
     MSG_BUG(message)
   endif
!   endif
! enddo ! ifreq
!  print *,"myproc", myproc
   if(optnonxsum2==0) then
      do iatom=1,green%oper(ifreq)%natom
       if(green%oper(ifreq)%matlu(iatom)%lpawu.ne.-1) then
         call xmpi_sum(green%oper(ifreq)%matlu(iatom)%mat,spacecomm,ierr)
       endif
     enddo ! iatom
     green%has_greenmatlu_xsum=1
   else if(optnonxsum2==1) then
     green%has_greenmatlu_xsum=0
   endif
!     if(ifreq==1.or.ifreq==11) then
!       write(std_out,*) ifreq,nproc,'after xsum'
!       call print_matlu(green%oper(ifreq)%matlu,green%oper(ifreq)%natom,2,-1,0)
!       ! ok
!     endif
 enddo ! ifreq
 !call xmpi_barrier(spacecomm)
 call xmpi_sum(green%occup%ks,spacecomm,ierr)
! write(std_out,*) 'afterxsum sym     %matlu(1)%mat(2,5,1,1,1) 1',green%oper(1)%matlu(1)%mat(2,5,1,1,1)

 if(prtopt/=0.and.prtopt>-100)  then
   write(message,'(2a)') ch10,&
&   '  == Local Green function has been computed and projected on local orbitals'
   call wrtout(std_out,message,'COLL')
 endif
! useless test
 if(abs(prtopt)>=4.and.prtopt>-100) then
   write(message,'(2a)') ch10,' == Green function is now printed for first frequency'
   call wrtout(std_out,message,'COLL')
   call print_oper(green%oper(1),9,paw_dmft,3)
   write(message,'(2a)') ch10,' == Green function is now printed for second frequency'
   call wrtout(std_out,message,'COLL')
   call print_oper(green%oper(2),9,paw_dmft,3)
   if(paw_dmft%dmft_nwlo>=11) then
     write(message,'(2a)') ch10,' == Green function is now printed for 11th frequency'
     call wrtout(std_out,message,'COLL')
     call print_oper(green%oper(11),9,paw_dmft,3)
   endif
 endif
! call flush(std_out)

 ABI_DEALLOCATE(Id)
 call timab(624,2,tsec)

end subroutine compute_green
!!***

!!****f* m_green/integrate_green
!! NAME
!! integrate_green
!!
!! FUNCTION
!!  integrate green function in the band index basis
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  green <type(green_type)>=green function  (green%oper(:))
!!  opt_ksloc=   1: do only integrate on the KS basis
!!               2: do the integration on the local basis
!!                 (This can work only if psichi are renormalized!!!)
!!               3: do both calculations and test the consistency of it
!!              -1: do the integration on the KS basis, but only
!!                      compute diagonal part of the band-band density matrix
!!                      in order to compute the total charge for fermi_green
!!  paw_dmft <type(m_paw_dmft)>= paw+dmft data
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  opt_nonxsum 0 : green(ifreq)%ks was broadcast in compute_green.
!!              1 : green(ifreq)%ks was not broadcast in compute_green, so a xsum will be
!!                  necessary here to compute the total number of electron.
!!
!! OUTPUT
!!   green%occup = occupations
!!
!! PARENTS
!!      m_dmft,fermi_green,impurity_solve,m_green,newton
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

subroutine integrate_green(cryst_struc,green,paw_dmft&
&  ,pawang,prtopt,opt_ksloc,opt_after_solver,opt_diff,opt_nonxsum)

!Arguments ------------------------------------
!type
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type),intent(inout) :: green
 !type(MPI_type), intent(in) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(paw_dmft_type), intent(inout) :: paw_dmft
 integer, intent(in) :: prtopt
 integer, optional, intent(in) :: opt_ksloc
 integer, optional, intent(in) :: opt_after_solver
 integer, optional, intent(in) :: opt_diff
 integer, optional, intent(in) :: opt_nonxsum

!local variables-------------------------------
 real(dp) :: tsec(2)
 integer :: iatom,ib,ib1,icomp_chloc,ifreq,ikpt,im,im1,ispinor,ispinor1,is
 integer :: mband,mbandc,myproc,natom,ndim,nkpt,nproc,nspinor
 integer :: nsppol,option
 integer :: optksloc,spacecomm,optaftsolv,optnonxsum
 complex(dpc) :: integral
 character(len=500) :: message
 complex(dpc), allocatable :: ff(:)
 type(matlu_type), allocatable :: matlu_temp(:)
 type(oper_type) :: occup_temp
 real(dp) :: diff_chloc
! real(dp), allocatable :: charge_loc_old(:,:)
! type(oper_type)  :: oper_c
! *********************************************************************

 DBG_ENTER("COLL")
 call timab(625,1,tsec)

 if(prtopt>0) then
   write(message,'(2a,i3,13x,a)') ch10,'   ===  Integrate green function'
   call wrtout(std_out,message,'COLL')
 endif
 if(green%w_type=="real") then
   message = 'integrate_green not implemented for real frequency'
   MSG_BUG(message)
 endif
 if(present(opt_nonxsum)) then
   optnonxsum=opt_nonxsum
 else
   optnonxsum=0
 endif

! Initialise spaceComm, myproc, and master
 spacecomm=paw_dmft%spacecomm
 myproc=paw_dmft%myproc
 nproc=paw_dmft%nproc


! Initialise integers
 mband   = paw_dmft%mband
 mbandc  = paw_dmft%mbandc
 nkpt    = paw_dmft%nkpt
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol
 natom   = paw_dmft%natom

! Initialize green%oper before calculation (important for xmpi_sum)
! allocate(charge_loc_old(paw_dmft%natom,paw_dmft%nsppol+1))
! if(.not.present(opt_diff)) then  ! if integrate_green is called in m_dmft after calculation of self
!   charge_loc_old=green%charge_matlu
! endif
 icomp_chloc=0

! Choose what to compute
 if(present(opt_ksloc)) then
   optksloc=opt_ksloc
 else
   optksloc=3
 endif
 if(present(opt_after_solver)) then
   optaftsolv=opt_after_solver
 else
   optaftsolv=0
 endif
 if(optaftsolv==1.and.abs(optksloc)/=2) then
    message = "integration of ks green function should not be done after call to solver : it has not been computed"
    MSG_BUG(message)
 endif
 if(abs(optksloc)>=2.and.green%has_greenmatlu_xsum==0) then
    write(message,'(4a)') ch10,&
&     "BUG: integrate_green is asked to integrate local green function",ch10,&
&    " and local green function was non broadcasted in compute_green"
    MSG_BUG(message)
 endif

! Allocations
 ABI_DATATYPE_ALLOCATE(matlu_temp,(natom))
 call init_matlu(natom,nspinor,nsppol,green%oper(1)%matlu(:)%lpawu,matlu_temp)

 ABI_ALLOCATE(ff,(green%nw))

! =================================================
! == Integrate Local Green function ===============
 if(abs(optksloc)/2==1) then ! optksloc=2 or 3
! =================================================
   call zero_matlu(green%occup%matlu,green%occup%natom)

! ==  Calculation of \int{G_{LL'}{\sigma\sigma',s}(R)(i\omega_n)}
   if(paw_dmft%lpsichiortho==1) then
!  - Calculation of frequency sum over positive frequency
     if (nspinor==1) option=1
     if (nspinor==2) option=2
     do iatom=1, natom
       ndim=2*green%oper(1)%matlu(iatom)%lpawu+1
       if(green%oper(1)%matlu(iatom)%lpawu.ne.-1) then
             do ispinor1 = 1, nspinor
           do ispinor = 1, nspinor
         do is = 1 , nsppol
               do im=1,ndim
                 do im1=1,ndim
                   do ifreq=1,green%nw
                     ff(ifreq)= &
&                     green%oper(ifreq)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
!                     write(std_out,*) green%omega(ifreq),ff(ifreq)," integrate green fw_lo"
   !if(present(iii).and.im==1.and.im1==1) write(std_out_default,*) ch10,ifreq,ff(ifreq),"#ff"
                   enddo
!                   call int_fct(ff,(im==im1).and.(ispinor==ispinor1),&
!&                   option,paw_dmft,integral)
                   call int_fct(ff,(im==im1).and.(ispinor==ispinor1),&
&                   2,paw_dmft,integral)  ! test_1
                   green%occup%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)=integral
   !if(present(iii).and.im==1.and.im1==1) write(std_out_default,*) ch10,'integral',im,im1,ifreq,integral
!                   if(im==2.and.im1==5.and.is==1.and.iatom==1) then
!                     write(std_out,*) " occup        %matlu(1)%mat(2,5,1,1,1)",green%occup%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
!                   endif
                 enddo
               enddo
             enddo ! ispinor1
           enddo ! ispinor
         enddo ! is
         matlu_temp(iatom)%mat=green%occup%matlu(iatom)%mat
       endif ! lpawu=/-1
     enddo ! iatom

!   Print density matrix if prtopt high
     if(abs(prtopt)>2) then
       write(message,'(a,a,i10,a)') ch10,"  = green%occup%matlu from int(gloc(w))"
       call wrtout(std_out,message,'COLL')
       call print_matlu(green%occup%matlu,natom,prtopt=3,opt_diag=-1)
     endif

!  - Symetrise: continue sum over k-point: Full BZ
     call sym_matlu(cryst_struc,green%occup%matlu,pawang,paw_dmft)
     if(abs(prtopt)>2) then
       write(message,'(a,a,i10,a)') ch10,&
&       "  = green%occup%matlu from int(gloc(w)) with symetrization"
       call wrtout(std_out,message,'COLL')
       call print_matlu(green%occup%matlu,natom,prtopt=3,opt_diag=-1)
     endif
     call sym_matlu(cryst_struc,matlu_temp,pawang,paw_dmft)

!  - Post-treatment for summation over negative and positive frequencies:
!    necessary in the case of nspinor==2 AND nspinor==1, but valid anywhere
!    N(ll'sigmasigma')= (N(ll'sigmasigma')+ N*(l'lsigma'sigma))/2
!    because [G_{LL'}^{sigma,sigma'}(iomega_n)]*= G_{L'L}^{sigma',sigma}(-iomega_n)
     if(nspinor>=1) Then
       do iatom=1, natom
         ndim=2*green%oper(1)%matlu(iatom)%lpawu+1
         if(green%oper(1)%matlu(iatom)%lpawu.ne.-1) then
               do ispinor1 = 1, nspinor
             do ispinor = 1, nspinor
           do is = 1 , nsppol
                 do im=1,ndim
                   do im1=1,ndim
                     green%occup%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)= &
&                          (matlu_temp(iatom)%mat(im,im1,is,ispinor,ispinor1)+ &
&                           conjg(matlu_temp(iatom)%mat(im1,im,is,ispinor1,ispinor)))/two
                   enddo
                 enddo
               enddo ! ispinor1
             enddo ! ispinor
           enddo ! isppol
           matlu_temp(iatom)%mat=green%occup%matlu(iatom)%mat
         endif
       enddo ! iatom
     endif
     if(abs(prtopt)>2) then
       write(message,'(a,a,i10,a)') ch10,&
&       "  = green%occup%matlu from int(gloc(w)) symetrized with post-treatment"
       call wrtout(std_out,message,'COLL')
       call print_matlu(green%occup%matlu,natom,prtopt=3,opt_diag=-1)
     endif
     if(optaftsolv==0) then
       call trace_oper(green%occup,green%charge_ks,green%charge_matlu,2)
       green%has_charge_matlu=2
       green%has_charge_matlu_solver=0
       icomp_chloc=1
     else if(optaftsolv==1.and.paw_dmft%dmft_solv/=4) then
!     else if(optaftsolv==1) then
       if(paw_dmft%dmft_solv==4) then
         write(message,'(a,a,a)') ch10,&
&        "  = WARNING: Double Counting will be computed with solver charge.." ,&
&        "might be problematic with Hirsch Fye QMC"
         call wrtout(std_out,message,'COLL')
       endif
!      This is done only when called from impurity_solver with solver
!      green function. (And if QMC is NOT used).
       if(paw_dmft%dmft_solv>=4.and.green%has_charge_matlu_solver/=2) then
         write(message,'(a,a,i3)') ch10,&
&        "  = BUG : has_charge_matlu_solver should be 2 and is",&
&        green%has_charge_matlu_solver
         MSG_BUG(message)
       endif
       if(paw_dmft%dmft_solv<=4) then
         call trace_oper(green%occup,green%charge_ks,green%charge_matlu_solver,2)
         green%has_charge_matlu_solver=2
       endif
     endif
   else
     write(message,'(a,4x,a,a,a,4x,a)') ch10,&
&     " Local basis is not (yet) orthonormal:",&
&     " local green function is thus not integrated",ch10,&
&     " Local occupations are computed from KS occupations"
     call wrtout(std_out,message,'COLL')
   endif

 endif ! optksloc
! =================================================



! =================================================
! == Integrate Kohn Sham Green function ===========
 if(mod(abs(optksloc),2)==1) then ! optksloc=1 or 3 or -1
!   green%occup%ks=czero ! important for xmpi_sum
!! =================================================
!! ==  Calculation of \int{G_{\nu\nu'}{k,s}(i\omega_n)}
   call init_oper(paw_dmft,occup_temp)
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
   occup_temp%ks=green%occup%ks
               !write(std_out,*) "occup%ks ik1",green%occup%ks(1,1,1,3)
               !write(std_out,*) "occup%ks ik2",green%occup%ks(1,2,1,3)
!  - Post-treatment for summation over negative and positive frequencies:
!    necessary in the case of nspinor==2, but valid everywhere
!    N(k,n_1,n_2)= (N(k,n_1,n_2)+ N*(k,n_2,n_1))/2
!    because [G_{k}^{n_1,n_2}(iomega_n)]*= G_{k}^{n_2,n_1}(-iomega_n)
         do ib1 = 1, mbandc
       do ib = 1, mbandc
     do ikpt = 1, nkpt
   do is = 1 , nsppol
           green%occup%ks(is,ikpt,ib,ib1)=&
&            (       occup_temp%ks(is,ikpt,ib,ib1)+ &
&              conjg(occup_temp%ks(is,ikpt,ib1,ib))   )/two
         enddo ! ib1
       enddo ! ib
     enddo ! ikpt
   enddo ! isppol
               !write(std_out,*) "occup%ks ik1 BB",green%occup%ks(1,1,1,3)
               !write(std_out,*) "occup%ks ik2 AA",green%occup%ks(1,2,1,3)
   call destroy_oper(occup_temp)
         do ib1 = 1, mbandc
       do ib = 1, mbandc
     do ikpt = 1, nkpt
   do is = 1 , nsppol
           paw_dmft%occnd(1,paw_dmft%include_bands(ib),&
&           paw_dmft%include_bands(ib1),ikpt,is)=dreal(green%occup%ks(is,ikpt,ib,ib1))
           paw_dmft%occnd(2,paw_dmft%include_bands(ib),&
&           paw_dmft%include_bands(ib1),ikpt,is)=dimag(green%occup%ks(is,ikpt,ib,ib1))
           if(nspinor==1 .and. nsppol==1) then
             paw_dmft%occnd(1,paw_dmft%include_bands(ib),&
&             paw_dmft%include_bands(ib1),ikpt,is)=two*dreal(green%occup%ks(is,ikpt,ib,ib1))
             paw_dmft%occnd(2,paw_dmft%include_bands(ib),&
&             paw_dmft%include_bands(ib1),ikpt,is)=two*dimag(green%occup%ks(is,ikpt,ib,ib1))
           endif
         enddo
       enddo
     enddo
   enddo

   if(optksloc>0) then
!  - Compute local occupations
     call loc_oper(green%occup,paw_dmft,1)
     if(abs(prtopt)>2) then
       write(message,'(a,a,i10,a)') ch10,&
&        "  = green%occup%matlu from projection of int(gks(w)) without symetrization"
       call wrtout(std_out,message,'COLL')
       call print_matlu(green%occup%matlu,natom,prtopt=3,opt_diag=-1)
     endif

!  - Symetrise: continue sum over k-point: Full BZ
     call sym_matlu(cryst_struc,green%occup%matlu,pawang,paw_dmft)
     if(abs(prtopt)>=2) then
!       write(message,'(a,a,i10,a)') ch10,&
!&        "  = green%occup%matlu from projection of int(gks(w)) with symetrization"
       write(message,'(a,a,i10,a)') ch10,&
&        "  = Occupation matrix from KS occupations"
       call wrtout(std_out,message,'COLL')
       call print_matlu(green%occup%matlu,natom,prtopt=3,opt_diag=0)
     endif

!  - If the trace of occup matrix in the LOCAL basis was not done
!  before because lpsichiortho/=1 , do it now
     if(paw_dmft%lpsichiortho/=1) then
       if(abs(prtopt)>0) then
         call trace_oper(green%occup,green%charge_ks,green%charge_matlu,2)
         green%has_charge_matlu=2
         icomp_chloc=1
       endif
     endif
   endif ! optksloc>0: only diagonal elements of G\nunu' are computed

!  - Compute trace over ks density matrix
   call trace_oper(green%occup,green%charge_ks,green%charge_matlu,1)
   if(abs(prtopt)>0) then
     write(message,'(a,a,f12.6)') ch10,&
&    "  ==  Total number of electrons from KS green function is :", green%charge_ks
     call wrtout(std_out,message,'COLL')
     write(message,'(8x,a,f12.6,a)') " (should be",paw_dmft%nelectval,")"
     call wrtout(std_out,message,'COLL')
   endif
 endif ! optksloc
! =================================================


! =================================================
! Tests and compute precision on local charge
! =================================================
!  - Check that if, renormalized psichi are used, occupations matrices
!    obtained directly from local green function or, through kohn sham
!    occupations are the same.
 if(abs(optksloc)==3) then ! optksloc= 3
   if(paw_dmft%lpsichiortho==1) then
     call diff_matlu("Local_projection_of_kohnsham_occupations ",&
&     "Integration_of_local_green_function ",&
&       green%occup%matlu,matlu_temp,natom,1,tol4)
     write(message,'(2a)') ch10,&
&     '  ***** => Calculations of Green function in KS and local spaces are coherent ****'
     call wrtout(std_out,message,'COLL')
   endif
 endif

!!***

! == Precision on charge_matlu (done only if local charge was computed ie not for optksloc=-1)
 if(icomp_chloc==1.and.paw_dmft%idmftloop>=1.and.present(opt_diff)) then ! if the computation was done here.
   if(green%has_charge_matlu_prev==2) then
     diff_chloc=zero
     do iatom=1, natom
       if(green%oper(1)%matlu(iatom)%lpawu.ne.-1) then
         do is=1,nsppol
           diff_chloc=diff_chloc+&
&           (green%charge_matlu_prev(iatom,is)-green%charge_matlu(iatom,is))**2
         enddo
       endif
     enddo
     if(sqrt(diff_chloc)<paw_dmft%dmft_lcpr) then
       green%ichargeloc_cv=1
       write(message,'(a,8x,a,e9.2,a,8x,a)') ch10,"Change of correlated number of electron =<",&
&        paw_dmft%dmft_lcpr,&
&       ch10,"DMFT Loop: Local charge is converged"
       call wrtout(std_out,message,'COLL')
     else
       green%ichargeloc_cv=0
       write(message,'(a,8x,a,e9.2,a,8x,a)') ch10,"Change of correlated number of electron  >",&
&        paw_dmft%dmft_lcpr,&
&       ch10,"DMFT Loop: Local charge is not converged"
       call wrtout(std_out,message,'COLL')
     endif
   endif
   green%charge_matlu_prev=green%charge_matlu
   green%has_charge_matlu_prev=2
 endif


 ABI_DEALLOCATE(ff)
 call destroy_matlu(matlu_temp,natom)
 ABI_DATATYPE_DEALLOCATE(matlu_temp)
! deallocate(charge_loc_old)
 call timab(625,2,tsec)
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
!!  cryst_struc <type(crystal_t)>= crystal structure data.
!!  mpi_enreg=informations about MPI parallelization
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  green  <type(green_type)>= green function data
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawprtvol  = option for printing
!!  self <type(self_type)>= variables related to self-energy
!!  opt_self = optional argument, if =1, upfold self-energy
!!
!! OUTPUT
!!
!! PARENTS
!!      m_dmft
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

subroutine icip_green(char1,cryst_struc,green,paw_dmft,pawang,pawprtvol,self,opt_self)

!Arguments ------------------------------------
!type
 type(crystal_t),intent(in) :: cryst_struc
 !type(MPI_type), intent(in) :: mpi_enreg
 type(green_type),intent(inout) :: green
 type(pawang_type),intent(in) :: pawang
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(self_type), intent(inout) :: self
 integer, intent(in) :: pawprtvol
 character(len=*), intent(in) :: char1
 integer, optional, intent(in) :: opt_self
!Local variables ------------------------------------
 integer :: option,optself,optlocks,prtopt_for_integrate_green,opt_nonxsum_icip
 character(len=500) :: message
! *********************************************************************
 green%whichgreen="LDA"
 prtopt_for_integrate_green=2

 if(present(opt_self)) then
   optself=opt_self
 else
   optself=0
 endif
 opt_nonxsum_icip=1
 if(paw_dmft%dmftcheck==1) then ! for fourier_green
 !  opt_nonxsum_icip=0
 endif

 call init_green(green,paw_dmft)

! useless test ok
! call printocc_green(green,1,paw_dmft,2)
! write(std_out,*)" printocc_green zero finished "

!== Compute green%oper(:)%ks
!== Deduce  green%oper(:)%matlu(:)%mat
 call compute_green(cryst_struc,green,paw_dmft,&
& pawang,pawprtvol,self,optself,opt_nonxsum=opt_nonxsum_icip)
 if(paw_dmft%dmft_prgn>=1.and.paw_dmft%dmft_prgn<=2) then
   optlocks=paw_dmft%dmft_prgn*2+1 ! if dmft_prgn==2 => do not print
   if(paw_dmft%lpsichiortho==1.and.pawprtvol>-100)  then
     call print_green(char1,green,optlocks,paw_dmft,pawprtvol)
!     call print_green('inicip',green,1,paw_dmft,pawprtvol=1,opt_wt=1)
   endif
 endif

!== Integrate green%oper(:)%ks
!== Integrate green%oper(:)%matlu(:)%mat
 call integrate_green(cryst_struc,green,paw_dmft,pawang,prtopt_for_integrate_green,opt_ksloc=3)!,opt_nonxsum=opt_nonxsum_icip)
!== Print green%oper(:)%ks
!== Print green%oper(:)%matlu(:)%mat
 if(char1=="LDA") then
   option=1
   if(self%oper(1)%matlu(1)%lpawu/=-1) then
     if(abs(real(self%oper(1)%matlu(1)%mat(1,1,1,1,1)))>tol7) then
! todo_ab: generalise this
       write(message,'(a,a,2(e15.4))') ch10,&
&        "Warning:  a LDA calculation is carried out and self is not zero"
       call wrtout(std_out,message,'COLL')
!       call abi_abort('COLL')
     endif
   endif
 else
   option=5
 endif
 if(pawprtvol>-100) then
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
!!  pawang <type(pawang)>=paw angular mesh and related data
!!
!! OUTPUT
!!
!! PARENTS
!!      impurity_solve,m_green,qmc_prep_ctqmc
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

subroutine fourier_green(cryst_struc,green,paw_dmft,pawang,opt_ksloc,opt_tw)

!Arguments ------------------------------------
!type
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type),intent(inout) :: green
 !type(MPI_type), intent(in) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
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
   MSG_BUG(message)
 endif

! Initialise temporary green function
 call init_green(green_temp,paw_dmft)
 call init_green_tau(green_temp,paw_dmft)

 !green%oper(:)%matlu(1)%mat(1,1,1,1,1)
 ABI_ALLOCATE(fw,(green%nw))
 ABI_ALLOCATE(ft,(green%dmftqmc_l))

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
               fw(ifreq)=green%oper(ifreq)%ks(is,ikpt,ib,ib1)
             enddo
             call fourier_fct(fw,ft,ib==ib1,green%dmftqmc_l,-1,paw_dmft) ! inverse fourier
             do itau=1,green%dmftqmc_l
               green_temp%oper_tau(itau)%ks(is,ikpt,ib,ib1)=ft(itau)
             enddo
             if(ib==ib1) then
               green%occup_tau%ks(is,ikpt,ib,ib1)=ft(1)+one
             else
               green%occup_tau%ks(is,ikpt,ib,ib1)=ft(1)
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
               green%oper_tau(itau)%ks(is,ikpt,ib,ib1)=&
&                (green_temp%oper_tau(itau)%ks(is,ikpt,ib,ib1)+ &
&                 conjg(green_temp%oper_tau(itau)%ks(is,ikpt,ib1,ib)))/two
               if(ib==ib1) then
                 green%occup_tau%ks(is,ikpt,ib,ib1)=green%oper_tau(1)%ks(is,ikpt,ib,ib1)+one
               else
                 green%occup_tau%ks(is,ikpt,ib,ib1)=green%oper_tau(1)%ks(is,ikpt,ib,ib1)
               endif
             enddo
           enddo ! ib1
         enddo ! ib
       enddo ! ikpt
     enddo ! isppol
     call loc_oper(green%occup_tau,paw_dmft,1)
     call sym_matlu(cryst_struc,green%occup_tau%matlu,pawang,paw_dmft)
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
                       fw(ifreq)=green%oper(ifreq)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
                     enddo
                     ! inverse fourier
!                     write(std_out,'(a)') "fourierbeforeposttreatement,ispinor,ispinor1,is,im,im1"
!                     write(std_out,'(a,5i4,f12.5,f12.5)') "fourier",ispinor,ispinor1,is,im,im1
!                     write(std_out,'(a,e12.5,e12.5)')&
!                     &"green%oper(4)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)"&
!&                     ,green%oper(4)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
                     call fourier_fct(fw,ft,(im==im1).and.(ispinor==ispinor1),green%dmftqmc_l,-1,paw_dmft)
                     do itau=1,green%dmftqmc_l
                       green_temp%oper_tau(itau)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)=ft(itau)
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
                       green%oper_tau(itau)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)=&
&                      ((green_temp%oper_tau(itau)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)+ &
&                      conjg(green_temp%oper_tau(itau)%matlu(iatom)%mat(im1,im,is,ispinor1,ispinor))))/two
!                       write(std_out,*) itau,green%oper_tau(itau)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
                       if((im==im1).and.(ispinor==ispinor1)) then
                         green%occup_tau%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)=&
&                         green%oper_tau(1)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)+one
                       else
                         green%occup_tau%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)=&
&                         green%oper_tau(1)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
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
       do ifreq=1,green%nw
         green%oper(ifreq)%matlu(iatom)%mat=czero
       enddo
       if(green%oper(1)%matlu(iatom)%lpawu.ne.-1) then
         ndim=2*green%oper(1)%matlu(iatom)%lpawu+1
         do is = 1 , nsppol
           do ispinor = 1, nspinor
             do ispinor1 = 1, nspinor
               do im=1,ndim
                 do im1=1,ndim
                   iparal=iparal+1
                   if(mod(iparal-1,nproc)==myproc) then
                     do itau=1,green%dmftqmc_l
                       ft(itau)=green%oper_tau(itau)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)
                     enddo
!                     if(im==1.and.im1==1) write(std_out,*) "ft(itau=1)",is,ft(1) ! debug
                     call fourier_fct(fw,ft,(im==im1).and.(ispinor==ispinor1),green%dmftqmc_l,1,paw_dmft)
                     do ifreq=1,green%nw
                       green%oper(ifreq)%matlu(iatom)%mat(im,im1,is,ispinor,ispinor1)=fw(ifreq)
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
       do ifreq=1,green%nw
       call xmpi_sum(green%oper(ifreq)%matlu(iatom)%mat,spacecomm,ierr)
       enddo
     enddo
   endif ! opt_ksloc=2
 endif ! opt_tw==-1
 ABI_DEALLOCATE(fw)
 ABI_DEALLOCATE(ft)
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
!! PARENTS
!!      m_dmft
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

subroutine check_fourier_green(cryst_struc,green,paw_dmft,pawang)

!Arguments ------------------------------------
!type
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type),intent(inout) :: green
 !type(MPI_type), intent(in) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(paw_dmft_type), intent(inout) :: paw_dmft

!local variables-------------------------------
 type(green_type) :: green_check
 character(len=500) :: message
! *********************************************************************
! Only imaginary frequencies here
 if(green%w_type=="real") then
   message = 'check_fourier_green not implemented for real frequency'
   MSG_BUG(message)
 endif

 call init_green(green_check,paw_dmft)
 call init_green_tau(green_check,paw_dmft)
 call copy_green(green,green_check,opt_tw=2)

 write(message,'(2a,i3,13x,a)') ch10,'   ===  Inverse Fourier Transform w->t of Weiss Field'
 call wrtout(std_out,message,'COLL')
 call fourier_green(cryst_struc,green_check,paw_dmft,&
& pawang,opt_ksloc=2,opt_tw=-1)

 write(message,'(3a)') ch10,' === Print (for check by user) of occupation matrix'&
&   ,' after  fourier transform with respect to initial one'
 call wrtout(std_out,message,'COLL')
 call printocc_green(green_check,6,paw_dmft,3)

 write(message,'(2a,i3,13x,a)') ch10,'   ===  Direct Fourier Transform t->w of Weiss Field'
 call wrtout(std_out,message,'COLL')
 call fourier_green(cryst_struc,green_check,paw_dmft,&
& pawang,opt_ksloc=2,opt_tw=1)
! call print_matlu(green%oper(1)%matlu,paw_dmft%natom,1) ! debug

 call integrate_green(cryst_struc,green_check,paw_dmft,&
& pawang,prtopt=2,opt_ksloc=2)

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
!! PARENTS
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
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
       occ1=occupfd(paw_dmft%eigen_lda(isppol,ikpt,ib),paw_dmft%fermie,paw_dmft%temp)
       occ2=real(green%occup%ks(isppol,ikpt,ib,ib))
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
!! Copyright (C) 2006-2020 ABINIT group (BAmadon)
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
!! PARENTS
!!      m_green
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
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
!! Copyright (C) 2006-2020 ABINIT group (BAmadon)
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
!! PARENTS
!!      m_green,qmc_prep_ctqmc
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
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
 ABI_ALLOCATE(procb2,(paw_dmft%dmft_nwlo))
 if(present(procb).and.present(myproc)) then
   do ifreq=1,paw_dmft%dmft_nwlo
    procb2(ifreq)=(procb(ifreq)==myproc)
   enddo
 else if(present(procb).and..not.present(myproc)) then
   write(message,'(a,a,2(e15.4))') ch10,&
&    "BUG: procb is present and not myproc in int_fct"
   MSG_BUG(message)
 else if(.not.present(procb).and.present(myproc)) then
   write(message,'(a,a,2(e15.4))') ch10,&
&    "BUG: procb is not present and myproc is in int_fct"
   MSG_BUG(message)
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
 ABI_DEALLOCATE(procb2)

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
!! Copyright (C) 2006-2020 ABINIT group (BAmadon)
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
!! PARENTS
!!      local_ks_green,m_green
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
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

   ABI_ALLOCATE(splined_li,(paw_dmft%dmft_nwli))
!   allocate(fw1(0:paw_dmft%dmft_nwlo-1))
   if(paw_dmft%dmft_log_freq==1) then
     ABI_ALLOCATE(omega_li,(1:paw_dmft%dmft_nwli))
     call construct_nwli_dmft(paw_dmft,paw_dmft%dmft_nwli,omega_li)
     call spline_c(paw_dmft%dmft_nwlo,paw_dmft%dmft_nwli,paw_dmft%omega_lo,&
&                   omega_li,splined_li,fw)
     ABI_DEALLOCATE(omega_li)
   else
     splined_li=fw
   endif
   call invfourier(splined_li,ft,paw_dmft%dmft_nwli,ltau,iflag,beta)
!   deallocate(fw1)
   ABI_DEALLOCATE(splined_li)

! == direct fourier transform
 else if(opt_four==1) then

   ABI_ALLOCATE(ftr,(ltau))

   iwarn=0
   do itau=1,ltau
     if(abs(aimag(ft(itau)))>tol12) then
       if(ldiag) then
         write(message,'(a,a,2(e15.4))') ch10,&
&          "green function is not real in imaginary time space",ft(itau)
         MSG_ERROR(message)
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

   ABI_ALLOCATE(tospline_li,(paw_dmft%dmft_nwli))
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
       ABI_ALLOCATE(omega_li,(1:paw_dmft%dmft_nwli))
       call construct_nwli_dmft(paw_dmft,paw_dmft%dmft_nwli,omega_li)
       call spline_c(paw_dmft%dmft_nwli,paw_dmft%dmft_nwlo,omega_li,&
&                 paw_dmft%omega_lo,fw,tospline_li)
       ABI_DEALLOCATE(omega_li)
     else
       fw=tospline_li
     endif
   endif

   ABI_DEALLOCATE(tospline_li)

   ABI_DEALLOCATE(ftr)
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
!! Copyright (C) 2006-2020 ABINIT group (BAmadon)
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
!! PARENTS
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
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
     ABI_ALLOCATE(omega_li,(1:size_fw2))
     call construct_nwli_dmft(paw_dmft,size_fw2,omega_li)
     call spline_c(size_fw1,size_fw2,paw_dmft%omega_lo(1:size_fw1),&
&                   omega_li(1:size_fw2),fw2,fw1)
     ABI_DEALLOCATE(omega_li)
   else
     fw2=fw1
   endif

! == direct fourier transform
 else if(opt_spline==1) then


   if(paw_dmft%dmft_log_freq==1) then
     ABI_ALLOCATE(omega_li,(1:size_fw2))
     call construct_nwli_dmft(paw_dmft,size_fw2,omega_li)
     call spline_c(size_fw2,size_fw1,omega_li(1:size_fw2),&
&                 paw_dmft%omega_lo(1:size_fw1),fw1,fw2)
     ABI_DEALLOCATE(omega_li)
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
!! PARENTS
!!      qmc_prep_ctqmc
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

subroutine occup_green_tau(green)

!Arguments ------------------------------------
!type
 type(green_type),intent(inout) :: green

!local variables-------------------------------
 integer :: iatom,lpawu
 complex(dpc), allocatable :: shift(:)
! *********************************************************************

 ABI_ALLOCATE(shift,(green%oper_tau(1)%natom))
 do iatom=1,green%oper_tau(1)%natom
   shift(iatom)=-cone
   lpawu=green%oper_tau(1)%matlu(iatom)%lpawu
   if(lpawu/=-1) then
!     do isppol=1,green%oper_tau(1)%nsppol
     green%occup_tau%matlu(iatom)%mat(:,:,:,:,:)= &
&      green%oper_tau(1)%matlu(iatom)%mat(:,:,:,:,:)
   endif
 enddo
 call shift_matlu(green%occup_tau%matlu,green%occup_tau%natom,shift)
 ABI_DEALLOCATE(shift)


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
!! PARENTS
!!
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

 function occupfd(eig,fermie,temp)


!Arguments ------------------------------------
!type
! Integrate analytic tail 1/(iw-mu)
 real(dp),intent(in) :: eig,fermie,temp
 real(dp) :: occupfd
!Local variables-------------------------------
! *********************************************************************

 if((eig-fermie) > zero) then
   occupfd=exp(-(eig-fermie)/temp)/(one+exp(-(eig-fermie)/temp))
 else
   occupfd=one/(one+exp((eig-fermie)/temp))
 endif

 end function occupfd
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
!! PARENTS
!!      m_green
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
!!
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
   ABI_ALLOCATE(proca,(nw,nproc/nw))
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
   ABI_DEALLOCATE(proca)

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

!!****f* ABINIT/greenldacompute_green
!! NAME
!! greenldacompute_green
!!
!! FUNCTION
!! Compute levels for ctqmc
!!
!! COPYRIGHT
!! Copyright (C) 1999-2020 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

 subroutine greenldacompute_green(cryst_struc,green,pawang,paw_dmft)

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(pawang_type), intent(in) :: pawang
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
  MSG_ERROR("greenlda%oper(1)%ks not allocated")
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
     !  write(6,*) paw_dmft%eigen_lda(isppol,ikpt,iband)
         green%oper(ifreq)%ks(isppol,ikpt,iband,iband)=&
         cone/(cmplx(0.d0,paw_dmft%omega_lo(ifreq),kind=dp)+paw_dmft%fermie-paw_dmft%eigen_lda(isppol,ikpt,iband))
       end do
     end do
   end do


!======================================================================
! Compute local Green's function
!======================================================================
   call loc_oper(green%oper(ifreq),paw_dmft,1)

!======================================================================
! Symetrize
!======================================================================
   call sym_matlu(cryst_struc,green%oper(ifreq)%matlu,pawang,paw_dmft)
 enddo
 write(message,'(a,2x,a,f13.5)') ch10," == Print LDA Green's function for last frequency"
 call wrtout(std_out,message,'COLL')
 call print_matlu(green%oper(paw_dmft%dmft_nwlo)%matlu,natom,1)

 end subroutine greenldacompute_green
!!***


!!****f* m_green/fermi_green
!! NAME
!! fermi_green
!!
!! FUNCTION
!!  Compute Fermi level for DMFT or LDA.
!!
!! COPYRIGHT
!! Copyright (C) 2006-2020 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! fermie: current value
!! f_precision: precision of f required
!! ITRLST: =1 if last iteration of DMFT
!! opt_noninter   : if one wants the LDA fermi level
!! max_iter : max number of iterations.
!!
!! OUTPUT
!! fermie: output value
!!
!! PARENTS
!!      m_dmft
!!
!! CHILDREN
!!      compute_green,integrate_green,newton,wrtout
!!
!! SOURCE

subroutine fermi_green(cryst_struc,green,paw_dmft,pawang,self)

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type),intent(inout) :: green
 type(paw_dmft_type), intent(inout) :: paw_dmft
 !type(MPI_type), intent(in) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(self_type), intent(inout) :: self

!Local variables-------------------------------
 integer :: ierr_hh,opt_noninter,max_iter
 real(dp):: x_precision,f_precision,fermi_old
! real(dp) :: hx
 character(len=500) :: message
!************************************************************************
!
 write(message,'(a,8x,a)') ch10,"  == Compute Fermi level"
 call wrtout(std_out,message,'COLL')

!=============
!headers
!=============
 write(message,'(2a)') ch10, "  |---Newton method to search Fermi level ------------|"
 call wrtout(std_out,message,'COLL')
 write(message,'(2a,f13.6)') ch10, "  |--- Initial value for Fermi level",paw_dmft%fermie
 call wrtout(std_out,message,'COLL')

!========================================
!Define precision and nb of iterations
!=========================================
 fermi_old=paw_dmft%fermie
 ierr_hh=0
 f_precision=paw_dmft%dmft_charge_prec
 !f_precision=0.01
 x_precision=tol5
!if(option==1) then
!f_precision=(erreursurlacharge)/100d0
!else
!f_precision=tol11
!endif
 max_iter=1 ! for tests only
 !write(6,*) "for tests max_iter=1"
 max_iter=50
 opt_noninter=4

!=====================
!Call newton method
!=====================
 write(message,'(a,4x,a,e13.6)') ch10," Precision required :",f_precision
 call wrtout(std_out,message,'COLL')
 if (f_precision<10_dp)  then
   call newton(cryst_struc,green,paw_dmft,pawang,self,&
&   paw_dmft%fermie,x_precision,max_iter,f_precision,ierr_hh,opt_noninter)
 end if

!===========================
!Deals with errors signals
!===========================
 if(ierr_hh==-314) then
   write(message,'(a)') "Warning, check Fermi level"
   call wrtout(std_out,message,'COLL')
!  call abi_abort('COLL')
   write(message,'(2a,f13.6)') ch10, "  |---  Final  value for Fermi level (check)",paw_dmft%fermie
   call wrtout(std_out,message,'COLL')
 else if (ierr_hh==-123) then
   write(message,'(a,f13.6)') " Fermi level is put to",fermi_old
   paw_dmft%fermie=fermi_old
   call wrtout(std_out,message,'COLL')

!  =====================================
!  If fermi level search was successful
!  =====================================
 else
   write(message,'(a,4x,a,e13.6)') ch10," Precision achieved on Fermi Level :",x_precision
   call wrtout(std_out,message,'COLL')
   write(message,'(4x,a,e13.6)') " Precision achieved on number of electrons :",f_precision
   call wrtout(std_out,message,'COLL')
   write(message,'(2a,f13.6)') ch10, "  |---  Final  value for Fermi level",paw_dmft%fermie
   call wrtout(std_out,message,'COLL')
 end if

!========================================================
!Check convergence of fermi level during DMFT iterations
!========================================================
 if(paw_dmft%idmftloop>=2) then
   if(abs(paw_dmft%fermie-fermi_old).le.paw_dmft%dmft_fermi_prec) then
!    write(message,'(a,8x,a,e9.2,a,8x,a,e12.5)') ch10,"|fermie(n)-fermie(n-1)|=<",paw_dmft%dmft_fermi_prec,ch10,&
     write(message,'(a,8x,a,e9.2,a,e9.2,a,8x,a,e12.5)') ch10,"|fermie(n)-fermie(n-1)|=",&
&     abs(paw_dmft%fermie-fermi_old),"<",paw_dmft%dmft_fermi_prec,ch10,&
&     "=> DMFT Loop: Fermi level is converged to:",paw_dmft%fermie
     call wrtout(std_out,message,'COLL')
     green%ifermie_cv=1
   else
     write(message,'(a,8x,a,2f12.5)') ch10,"DMFT Loop: Fermi level is not converged:",&
&     paw_dmft%fermie
     call wrtout(std_out,message,'COLL')
     green%ifermie_cv=0
   end if
 end if
 write(message,'(2a)') ch10, "  |---------------------------------------------------|"
 call wrtout(std_out,message,'COLL')
!

!==========================================================
!Recompute full green function including non diag elements
!==========================================================
 call compute_green(cryst_struc,green,paw_dmft,pawang,0,self,opt_self=1,opt_nonxsum=1)
 call integrate_green(cryst_struc,green,paw_dmft,pawang,prtopt=0,opt_ksloc=3) !,opt_nonxsum=1)

 return
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
!! Copyright (C) 2006-2020 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  x_input      : input of x
!!  x_precision  : required precision on x
!!  max_iter     : maximum number of iterations
!!  opt_noninter
!!
!! OUTPUT
!!  f_precision  : output precision on function F
!!  ierr_hh      : different from zero if an error occurs
!!
!! PARENTS
!!      fermi_green
!!
!! CHILDREN
!!      compute_green,integrate_green
!!
!! SOURCE

subroutine newton(cryst_struc,green,paw_dmft,pawang,self,&
& x_input,x_precision,max_iter,f_precision,ierr_hh,opt_noninter,opt_algo)

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type),intent(inout) :: green
 !type(MPI_type), intent(in) :: mpi_enreg
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(pawang_type),intent(in) :: pawang
 type(self_type), intent(inout) :: self
 integer,intent(in) :: opt_noninter,max_iter
 integer,intent(out) :: ierr_hh
 real(dp),intent(inout) :: x_input,x_precision
 real(dp),intent(inout) :: f_precision
 real(dp),intent(in), optional :: opt_algo
!Local variables-------------------------------
 integer iter
 real(dp) Fx,Fxprime,Fxdouble,xold,x_before,Fxoptimum
 real(dp) :: nb_elec_x
 integer option,optalgo
 logical l_minus,l_plus
 real(dp) :: x_minus,x_plus,x_optimum
 character(len=500) :: message
! *********************************************************************
 x_minus=-10_dp
 x_plus=-11_dp
 xold=-12_dp

 if(present(opt_algo)) then
   optalgo=opt_algo
 else
   optalgo=1 ! newton
 end if

 x_input=paw_dmft%fermie
 ierr_hh=0
 option =2  ! Halley method
 option =1  ! Newton method

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

 l_minus=.false.
 l_plus=.false.
 Fxoptimum=1_dp
!========================================
!start iteration to find fermi level
!========================================
 do iter=1, max_iter

!  ========================================
!  If zero is located between two values: apply newton method or dichotomy
!  ========================================
   if(l_minus.and.l_plus) then

!    ==============================================
!    Compute the function and derivatives for newton
!    ==============================================
     call function_and_deriv(cryst_struc,f_precision,green,iter,paw_dmft,pawang,self,&
&     x_input,x_before,x_precision,Fx,Fxprime,Fxdouble,opt_noninter,option)

!    Apply stop criterion on Fx
     if(abs(Fx) < f_precision) then
!      write(message,'(a,2f12.6)') "Fx,f_precision",Fx,f_precision
!      call wrtout(std_out,message,'COLL')
       x_precision=x_input-xold
       return
     end if
     if(iter==max_iter) then
       write(message,'(a,2f12.6)') "   Fermi level could not be found"
       call wrtout(std_out,message,'COLL')
       x_input=x_optimum
       ierr_hh=-123
       return
     end if

!    Cannot divide by Fxprime if too small
     if(abs(Fxprime) .le. 1.e-15)then
       ierr_hh=-314
       write(message,'(a,f12.7)') "Fxprime=",Fxprime
       call wrtout(std_out,message,'COLL')
       return
     end if

     x_precision=x_input-xold

!    ==============================================
!    Newton/Halley's  formula for next iteration
!    ==============================================
     xold=x_input
     if(option==1) x_input=x_input-Fx/Fxprime
     if(option==2) x_input=x_input-2*Fx*Fxprime/(2*Fxprime**2-Fx*Fxdouble)

!    ==============================================
!    If newton does not work well, use dichotomy.
!    ==============================================
     if(x_input<x_minus.or.x_input>x_plus) then
       call compute_nb_elec(cryst_struc,green,paw_dmft,pawang,self,&
&       Fx,opt_noninter,nb_elec_x,xold)
       write(message,'(a,3f12.6)') " ---",x_input,nb_elec_x,Fx
       call wrtout(std_out,message,'COLL')
       if(Fx>0) then
         x_plus=xold
       else if(Fx<0) then
         x_minus=xold
       end if
       x_input=(x_plus+x_minus)/2.d0

     end if
!    write(std_out,'(a,2f12.6)') " Q(xold) and dQ/dx=",Fx,Fxprime
!    write(std_out,'(a,f12.6)') " =>  new Fermi level",x_input
!    ========================================
!    Locate zero between two values
!    ========================================
   else
     call compute_nb_elec(cryst_struc,green,paw_dmft,pawang,self,&
&     Fx,opt_noninter,nb_elec_x,x_input)
     write(message,'(a,3f12.6)') "  --",x_input,nb_elec_x,Fx
!    Possible improvement for large systems, removed temporarely for
!    automatic tests: more study is necessary: might worsen the convergency
!    if(iter==1) then
!    f_precision=max(abs(Fx/50),f_precision)
!    write(message,'(a,4x,a,e12.6)') ch10," Precision required changed to:",f_precision
!    call wrtout(std_out,message,'COLL')
!    endif
     call wrtout(std_out,message,'COLL')
     if(Fx>0) then
       l_plus=.true.
       x_plus=x_input
       x_input=x_input-0.02
     else if(Fx<0) then
       l_minus=.true.
       x_minus=x_input
       x_input=x_input+0.02
     end if

   end if

   if(abs(Fx)<abs(Fxoptimum)) then
     Fxoptimum=Fx
     x_optimum=x_input
   end if



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

 end do
!--- End of iterations


 ierr_hh=1
 return

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
!!  x_input      : input of x
!!  option       : if 1 compute only first derivative
!!                 if 2 compute the two first derivatives.
!!  opt_noninter
!!
!! OUTPUTS
!!  Fx           : Value of F(x)
!!  Fxprime      : Value of F'(x)
!!  Fxdouble     : Value of F''(x)
!!
!! PARENTS
!!      newton
!!
!! CHILDREN
!!      compute_green,integrate_green
!!
!! SOURCE

subroutine function_and_deriv(cryst_struc,f_precision,green,iter,paw_dmft,pawang,self&
& ,x_input,x_old,x_precision,Fx,Fxprime,Fxdouble,opt_noninter,option)



!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type),intent(inout) :: green
 !type(MPI_type), intent(in) :: mpi_enreg
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(pawang_type),intent(in) :: pawang
 type(self_type), intent(inout) :: self
 integer,intent(in) :: iter,opt_noninter,option
 real(dp),intent(inout)  :: f_precision,x_input,x_precision
 real(dp),intent(out) :: Fx,Fxprime,Fxdouble
 real(dp),intent(inout) :: x_old
!Local variables-------------------------------
 real(dp) :: deltax,nb_elec_x,Fxmoins,Fxplus,xmoins,xplus,x0
 character(len=500) :: message
! *********************************************************************

!  choose deltax: for numeric evaluation of derivative
   if(iter==1) then
!    deltax=0.02
   end if
!  deltax=max((x_input-x_old)/10.d0,min(0.00001_dp,x_precision/100_dp))
   deltax=min(0.00001_dp,x_precision/100_dp)  ! small but efficient
!  endif
!  write(std_out,*) "iter,x_input,deltax",iter,x_input,deltax
   x0=x_input
   xmoins=x0-deltax
   xplus=x0+deltax

   call compute_nb_elec(cryst_struc,green,paw_dmft,pawang,self,&
&   Fx,opt_noninter,nb_elec_x,x0)

   write(message,'(a,3f12.6)') "  - ",x0,nb_elec_x,Fx
   call wrtout(std_out,message,'COLL')
!  write(std_out,*) "Fx", Fx
   if(abs(Fx)<f_precision) return

   call compute_nb_elec(cryst_struc,green,paw_dmft,pawang,self,&
&   Fxplus,opt_noninter,nb_elec_x,xplus)

   write(message,'(a,3f12.6)') "  - ",xplus,nb_elec_x,Fxplus
   call wrtout(std_out,message,'COLL')

   if(option==2) then
     call compute_nb_elec(cryst_struc,green,paw_dmft,pawang,self,&
&     Fxmoins,opt_noninter,nb_elec_x,xmoins)

     write(message,'(a,3f12.6)') "  - ",xmoins,nb_elec_x,Fxmoins
     call wrtout(std_out,message,'COLL')
   end if

   if(option==1) then
     Fxprime=(Fxplus-Fx)/deltax
   else if (option==2) then
     Fxprime=(Fxplus-Fxmoins)/(2*deltax)
     Fxdouble=(Fxplus+Fxmoins-2*Fx)/(deltax**2)
   end if
!  write(std_out,*) "after computation of Fxprime",myid
   if(Fxprime<zero) then
     write(message,'(a,f12.6)') "  Warning: slope of charge versus fermi level is negative !", Fxprime
     call wrtout(std_out,message,'COLL')
   end if
   x_old=x_input

   return
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
!! PARENTS
!!      newton
!!
!! CHILDREN
!!      compute_green,integrate_green
!!
!! SOURCE

subroutine compute_nb_elec(cryst_struc,green,paw_dmft,pawang,self,&
&  Fx,opt_noninter,nb_elec_x,fermie)



!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type),intent(inout) :: green
 !type(MPI_type), intent(in) :: mpi_enreg
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(pawang_type),intent(in) :: pawang
 type(self_type), intent(inout) :: self
 integer,intent(in) :: opt_noninter
 real(dp),intent(in)  :: fermie
 real(dp),intent(out) :: Fx,nb_elec_x
! *********************************************************************
   paw_dmft%fermie=fermie
   call compute_green(cryst_struc,green,paw_dmft,pawang,0,self,opt_self=1,&
&   opt_nonxsum=1,opt_nonxsum2=1)
   call integrate_green(cryst_struc,green,paw_dmft,pawang,prtopt=0,&
&   opt_ksloc=-1) !,opt_nonxsum=1)
!  opt_ksloc=-1, compute total charge
   nb_elec_x=green%charge_ks
   Fx=nb_elec_x-paw_dmft%nelectval

   if(opt_noninter==1) then
   end if
 end subroutine compute_nb_elec
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
!! Copyright (C) 1999-2020 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cryst_struc
!!  istep    =  step of iteration for LDA.
!!  lda_occup
!!  mpi_enreg=informations about MPI parallelization
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!
!! OUTPUT
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!
!! PARENTS
!!      m_dmft
!!
!! CHILDREN
!!      fourier_fct,wrtout
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
 ABI_ALLOCATE(tau,(ltau))
 do itau=1,ltau
   tau(itau)=float(itau-1)/float(ltau)/paw_dmft%temp
 end do
 beta=one/paw_dmft%temp

!Only imaginary frequencies here
 if(green%w_type=="real") then
   message = ' compute_energy not implemented for real frequency'
   MSG_BUG(message)
 end if

!=========================================
!Compute local band ks green function
! should be computed in compute_green: it would be less costly in memory.
!=========================================
 ABI_ALLOCATE(loc_ks,(nsppol,mbandc,paw_dmft%dmft_nwlo))
 if(green%oper(1)%has_operks==1) then
   loc_ks(:,:,:)=czero
   do isppol=1,nsppol
     do iband=1,mbandc
       do ifreq=1,paw_dmft%dmft_nwlo
         do ikpt=1,nkpt
           loc_ks(isppol,iband,ifreq)=loc_ks(isppol,iband,ifreq)+  &
&           green%oper(ifreq)%ks(isppol,ikpt,iband,iband)*paw_dmft%wtk(ikpt)
         end do
       end do
     end do
   end do
 else
   message = ' green fct is not computed in ks space'
   MSG_BUG(message)
 end if

!=========================================
!Compute fourier transformation
!=========================================

 ABI_ALLOCATE(loc_ks_tau,(nsppol,mbandc,ltau))
 ABI_ALLOCATE(fw,(paw_dmft%dmft_nwlo))
 ABI_ALLOCATE(ft,(ltau))
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
 ABI_DEALLOCATE(fw)
 ABI_DEALLOCATE(ft)
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
   ABI_ALLOCATE(unitgreenlocks_arr,(nsppol))
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
   ABI_DEALLOCATE(unitgreenlocks_arr)
 end if

!Deallocations
 ABI_DEALLOCATE(loc_ks)
 ABI_DEALLOCATE(loc_ks_tau)
 ABI_DEALLOCATE(tau)

end subroutine local_ks_green
!!***


END MODULE m_green
!!***
