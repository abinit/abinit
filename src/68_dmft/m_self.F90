!!****m* ABINIT/m_self
!! NAME
!!  m_self
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

MODULE m_self

 use defs_basis
 use m_errors
 use m_abicore

 use m_datafordmft, only : compute_levels
 use m_fstrings, only : int2char4
 use m_hu, only : hu_type
 use m_io_tools, only : get_unit
 use m_matlu, only : add_matlu,copy_matlu,destroy_matlu,diag_matlu,fac_matlu, &
                   & init_matlu,matlu_type,print_matlu,rotate_matlu,xmpi_matlu,zero_matlu
 use m_oper, only : destroy_oper,gather_oper,init_oper,oper_type,print_oper
 use m_paw_dmft, only : mpi_distrib_dmft_type,paw_dmft_type
 use m_paw_exactDC, only : compute_exactDC
 use m_pawtab, only : pawtab_type
 use m_xmpi, only : xmpi_bcast,xmpi_sum

#ifdef HAVE_GPU_MARKERS
 use m_nvtx_data
#endif

 implicit none

 private

 public :: alloc_self
 public :: initialize_self
 public :: destroy_self
 public :: print_self
 public :: rw_self
 public :: dc_self
 public :: new_self
 !public :: make_qmcshift_self
 public :: selfreal2imag_self

!!***

!!****t* m_self/self_type
!! NAME
!!  self_type
!!
!! FUNCTION
!!  This structured datatype contains the necessary data
!!
!! SOURCE

 type, public :: self_type ! for each atom

  integer :: dmft_nwli
  ! Linear index of the last imaginary frequency

  integer :: dmft_nwlo
  ! Number of imaginary frequencies

  integer :: has_moments
  ! =1 if the high-frequency moments are computed

  integer :: iself_cv
  ! Integer for convergence of self-energy

  integer :: nmoments
  ! Number of high-frequency moments which will be computed

  integer :: nw
  ! Number of frequencies (equal to dmft_nwlo only if w_type="imag")

  character(len=4) :: w_type
  ! Type of frequencies used ("real" or "imag")

  !real(dp), allocatable :: qmc_shift(:)
  ! value of frequencies

  !real(dp), allocatable :: qmc_xmu(:)
  ! value of frequencies

  type(oper_type) :: hdc
  ! Operator for double counting

  type(oper_type), allocatable :: moments(:)
  ! High-frequency moments

  type(oper_type), allocatable :: oper(:)
  ! Operator for self-energy, for each frequency

  real(dp), ABI_CONTIGUOUS pointer :: omega(:) => null()
  ! Value of frequencies

  type(mpi_distrib_dmft_type), pointer :: distrib => null()
  ! Datastructure for MPI parallelization

 end type self_type
!!***

!----------------------------------------------------------------------


CONTAINS  !========================================================================================
!!***

!!****f* m_self/alloc_self
!! NAME
!! alloc_self
!!
!! FUNCTION
!!  Allocate variables used in type self_type.
!!
!! INPUTS
!!  self <type(self_type)>= variables related to self-energy
!!  paw_dmft <type(paw_dmft_type)> =  variables related to self-consistent DFT+DMFT calculations.
!!  opt_oper = 1  Allocate only quantities in the KS basis.
!!             2  Allocate only quantities in the local basis.
!!             3  Allocate quantities in both the KS and local basis.
!!  wtype = "real" Self energy will be computed for real frequencies
!!        = "imag" (default) Self energy will be computed for imaginary frequencies
!!  opt_moments = 1 to allocate the high-frequency moments
!!
!! OUTPUTS
!!  self <type(self_type)>= variables related to self-energy
!!
!! SOURCE

subroutine alloc_self(self,paw_dmft,opt_oper,wtype,opt_moments)

!Arguments ------------------------------------
 type(self_type), intent(inout) :: self
 type(paw_dmft_type), target, intent(in) :: paw_dmft
 integer, optional, intent(in) :: opt_moments,opt_oper
 character(len=4), optional :: wtype
!Local variables ------------------------------------
 integer :: i,ifreq,mkmem,optmoments,optoper,shift
!************************************************************************

 optoper = 2
 self%w_type = "imag"
 self%nmoments = 0
 optmoments = 0

 if (present(opt_oper)) optoper  = opt_oper
 if (present(wtype)) self%w_type = wtype
 if (present(opt_moments)) optmoments = opt_moments

 self%has_moments = optmoments

 if (self%w_type == "imag") then
   self%nw = paw_dmft%dmft_nwlo
   self%omega => paw_dmft%omega_lo(:)
   self%distrib => paw_dmft%distrib
 else if (self%w_type == "real") then
   self%nw = size(paw_dmft%omega_r(:))
   self%omega => paw_dmft%omega_r(:)
   self%distrib => paw_dmft%distrib_r
 end if ! w_type

 self%dmft_nwlo = paw_dmft%dmft_nwlo
 self%dmft_nwli = paw_dmft%dmft_nwli
 self%iself_cv  = 0

 call init_oper(paw_dmft,self%hdc,opt_ksloc=optoper)

 ABI_MALLOC(self%oper,(self%nw))
 do ifreq=1,self%nw
   call init_oper(paw_dmft,self%oper(ifreq),opt_ksloc=optoper)
 end do ! ifreq

 if (optmoments == 1) then
   self%nmoments = 4
   shift = self%distrib%shiftk
   mkmem = self%distrib%nkpt_mem(self%distrib%me_kpt+1)
   ABI_MALLOC(self%moments,(self%nmoments))
   call init_oper(paw_dmft,self%moments(1),nkpt=mkmem,shiftk=shift,opt_ksloc=2)
   do i=2,self%nmoments
     call init_oper(paw_dmft,self%moments(i),nkpt=mkmem,shiftk=shift,opt_ksloc=3)
   end do ! i
 end if ! optmoments

 !if (paw_dmft%dmft_solv == 4) then
 !  ABI_MALLOC(self%qmc_shift,(paw_dmft%natom))
 !  ABI_MALLOC(self%qmc_xmu,(paw_dmft%natom))
 !  self%qmc_shift(:) = zero
 !  self%qmc_xmu(:) = zero
 !end if ! dmft_solv=4

end subroutine alloc_self
!!***

!!****f* m_self/initialize_self
!! NAME
!! initialize_self
!!
!! FUNCTION
!!  Initialize self-energy.
!!
!! INPUTS
!!  self <type(self_type)>= variables related to self-energy
!!  paw_dmft <type(paw_dmft_type)> =  variables related to self-consistent DFT+DMFT calculations.
!!  wtype = "real" Self energy will be computed for real frequencies
!!        = "imag" (default) Self energy will be computed for imaginary frequencies
!!  opt_moments = 1 to compute the high frequency moments
!!
!! OUTPUTS
!!  self <type(self_type)>= variables related to self-energy
!!
!!
!! SOURCE

subroutine initialize_self(self,paw_dmft,wtype,opt_moments)

!Arguments ------------------------------------
 type(self_type), intent(inout) :: self
 type(paw_dmft_type), intent(in) :: paw_dmft
 character(len=4), optional, intent(in) :: wtype
 integer, optional, intent(in) :: opt_moments
!Local variables ------------------------------------
! character(len=500) :: message
 integer :: optmoments
 character(len=4) :: wtype2
!************************************************************************

 optmoments = 0
 wtype2 = "imag"
 if (present(wtype)) wtype2 = wtype
 if (present(opt_moments)) optmoments = opt_moments

 call alloc_self(self,paw_dmft,opt_oper=2,wtype=wtype2,opt_moments=optmoments) !  opt_oper=1 is not useful and not implemented
! if(paw_dmft%dmft_rslf==1.and.opt_read==1) then
!   call rw_self(cryst_struc,self,mpi_enreg,paw_dmft,pawtab,pawprtvol=2,opt_rw=1)
! endif
! write(message,'(a,a)') ch10,"   Self-energy for large frequency is"
! call wrtout(std_out,message,'COLL')
! call print_matlu(self%oper(paw_dmft%dmft_nwlo)%matlu,  &
!&                 paw_dmft%natom,3)

end subroutine initialize_self
!!***

!!****f* m_self/destroy_self
!! NAME
!! destroy_self
!!
!! FUNCTION
!!  Deallocate self
!!
!! INPUTS
!!  self <type(self_type)>= variables related to self-energy
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_self(self)

!Arguments ------------------------------------
 type(self_type), intent(inout) :: self
!Local variables-------------------------------
 integer :: i,ifreq
! *********************************************************************

 if (allocated(self%oper)) then
   do ifreq=1,self%nw
     call destroy_oper(self%oper(ifreq))
   end do
   ABI_FREE(self%oper)
 end if

 call destroy_oper(self%hdc)

 if (allocated(self%moments)) then
   do i=1,self%nmoments
     call destroy_oper(self%moments(i))
   end do ! i
   ABI_FREE(self%moments)
 end if
 !if (allocated(self%qmc_shift)) ABI_FREE(self%qmc_shift)
 !if (allocated(self%qmc_xmu)) ABI_FREE(self%qmc_xmu)
 self%distrib => null()
 self%omega => null()

end subroutine destroy_self
!!***

!!****f* m_self/print_self
!! NAME
!! print_self
!!
!! FUNCTION
!!  Print self-energy
!!
!! INPUTS
!!  self <type(self_type)>= variables related to self-energy
!!  prtdc = print double counting if equal to "print_dc"
!!  paw_dmft <type(paw_dmft_type)> =  variables related to self-consistent DFT+DMFT calculations.
!!  prtopt = integer which specifies the amount of printing in the subroutine called
!!
!! OUTPUT
!!  self <type(self_type)>= variables related to self-energy
!!
!! SOURCE

subroutine print_self(self,prtdc,paw_dmft,prtopt)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(self_type), intent(in)  :: self
 character(len=*), intent(in) :: prtdc
 integer, intent(in) :: prtopt
!Local variables-------------------------------
 character(len=500) :: message
! *********************************************************************

 write(message,'(2a)') ch10,"  == The self-energy for smallest frequency is   == "
 call wrtout(std_out,message,'COLL')
 call print_oper(self%oper(1),1,paw_dmft,prtopt)
! write(message,'(2a)') ch10,"  == The self-energy for small (3) frequency is   == "
! call wrtout(std_out,message,'COLL')
! call print_oper(self%oper(3),1,paw_dmft,prtopt)
 write(message,'(2a)') ch10,"  == The self-energy for largest frequency is   == "
 call wrtout(std_out,message,'COLL')
 call print_oper(self%oper(self%nw),1,paw_dmft,prtopt)
 if (prtdc == "print_dc") then
   write(message,'(2a)') ch10,"  == The double counting potential is  == "
   call wrtout(std_out,message,'COLL')
   call print_matlu(self%hdc%matlu(:),paw_dmft%natom,prtopt)
 end if ! prtdc

end subroutine print_self
!!***

!!****f* m_self/dc_self
!! NAME
!! dc_self
!!
!! FUNCTION
!!  Computes the double counting
!!
!! INPUTS
!!  charge_loc : local charge for each polarization and each atom
!!  hu <type(hu_type)>= U interaction
!!  paw_dmft <type(paw_dmft_type)> =  variables related to self-consistent DFT+DMFT calculations.
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!  occ_matlu : local occupation matrix
!!
!! OUTPUT
!!  hdc : double counting
!!
!! SOURCE

subroutine dc_self(charge_loc,hdc,hu,paw_dmft,pawtab,occ_matlu)

!Arguments ------------------------------------
!type
 !type(crystal_t),intent(in) :: cryst_struc
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(matlu_type), intent(inout) :: hdc(paw_dmft%natom)
 real(dp), intent(in) :: charge_loc(paw_dmft%nsppol+1,paw_dmft%natom)
 type(hu_type), intent(in) :: hu(paw_dmft%ntypat)
 type(pawtab_type), intent(inout) :: pawtab(paw_dmft%ntypat)
 type(matlu_type), intent(in) :: occ_matlu(paw_dmft%natom)
 !type(hu_type),intent(inout) :: hu(cryst_struc%ntypat)
 !integer, intent(in) :: dmft_dc
!Local variables-------------------------------
 integer :: dmft_dc,iatom,iatomc,ierr,im,ispinor,isppol
 integer :: itypat,lpawu,natom,ndim,nspinor,nsppol
 real(dp) :: dc,jpawu,ntot,upawu
 logical :: amf,fll,nmdc
 character(len=500) :: message
 complex(dpc), allocatable :: occ(:,:),vdc(:,:)
! *********************************************************************

 dmft_dc = paw_dmft%dmft_dc
 natom   = paw_dmft%natom
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol

 amf  = (dmft_dc == 2 .or. dmft_dc == 6) ! AMF double counting
 fll  = (dmft_dc == 1 .or. dmft_dc == 4 .or. dmft_dc == 5) ! FLL double counting
 nmdc = (nspinor == 2 .or. (dmft_dc >= 5 .and. dmft_dc <= 8)) ! non-magnetic double counting

 if ((.not. fll) .and. (.not. amf) .and. dmft_dc /= 7 .and. dmft_dc /= 8) &
   & ABI_ERROR("not implemented")
 if (amf .and. (nspinor == 2)) then
   write(message,'(a,i4,i4,2x,e20.10)') " AMF Double counting is under test for SOC"
   ABI_WARNING(message)
 end if

 if (dmft_dc == 8) then
   paw_dmft%edc(:)   = zero
   paw_dmft%edcdc(:) = zero
   call zero_matlu(hdc(:),natom)
 end if

 iatomc = -1

 do iatom=1,natom
   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   iatomc = iatomc + 1
   hdc(iatom)%mat(:,:,:) = czero
   ntot   = charge_loc(nsppol+1,iatom)
   ndim   = 2*lpawu + 1
   itypat = paw_dmft%typat(iatom)
   upawu  = hu(itypat)%upawu
   jpawu = merge(zero,hu(itypat)%jpawu,dmft_dc==4)

   if (dmft_dc == 8) then
     if (mod(iatomc,paw_dmft%nproc) /= paw_dmft%myproc) cycle
     ABI_MALLOC(occ,(ndim,ndim))
     ABI_MALLOC(vdc,(ndim,ndim))
     occ(:,:) = czero
     do isppol=1,nsppol
       do ispinor=1,nspinor
         occ(:,:) = occ(:,:) + occ_matlu(iatom)%mat(1+(ispinor-1)*ndim:ispinor*ndim,1+(ispinor-1)*ndim:ispinor*ndim,isppol)
       end do ! ispinor
     end do ! isppol
     if (nsppol == 1 .and. nspinor == 1) occ(:,:) = occ(:,:) * two
     call compute_exactDC(lpawu,pawtab(itypat),paw_dmft%radgrid(itypat),occ(:,:), &
                        & vdc(:,:),paw_dmft%edc(iatom),paw_dmft%edcdc(iatom),paw_dmft%ixc)
     do isppol=1,nsppol
       do ispinor=1,nspinor
         hdc(iatom)%mat(1+(ispinor-1)*ndim:ispinor*ndim,1+(ispinor-1)*ndim:ispinor*ndim,isppol) = vdc(:,:)
       end do ! ispinor
     end do ! isppol
     ABI_FREE(occ)
     ABI_FREE(vdc)
   else
     if (nmdc) then ! Non-magnetic DC
       if (fll) dc = upawu*(ntot-half) - half*jpawu*(ntot-one)
       if (amf) dc = upawu*ntot*half + (upawu-jpawu)*ntot*half*dble(2*lpawu)/dble(2*lpawu+1)
       if (dmft_dc == 7) dc = upawu*(dble(paw_dmft%dmft_nominal(iatom))-half) - &
         & half*jpawu*(dble(paw_dmft%dmft_nominal(iatom))-one)
     end if ! nmdc
     ndim = nspinor * ndim
     do isppol=1,nsppol
       if (.not. nmdc) then ! Magnetic DC
         if (fll) dc = upawu*(ntot-half) - jpawu*(charge_loc(isppol,iatom)-half)
         if (amf) dc = upawu*charge_loc(min(3-isppol,nsppol),iatom) + &
              & (upawu-jpawu)*charge_loc(isppol,iatom)*dble(2*lpawu)/dble(2*lpawu+1)
       end if ! not nmdc
       do im=1,ndim
         hdc(iatom)%mat(im,im,isppol) = dc
       end do ! im
     end do ! isppol
   end if ! dc=8
 end do ! iatom

 if (dmft_dc == 8) then
   call xmpi_sum(paw_dmft%edc(:),paw_dmft%spacecomm,ierr)
   call xmpi_sum(paw_dmft%edcdc(:),paw_dmft%spacecomm,ierr)
   call xmpi_matlu(hdc(:),natom,paw_dmft%spacecomm)
 end if ! dmft_dc=8

end subroutine dc_self
!!***

!!****f* m_self/rw_self
!! NAME
!! rw_self
!!
!! FUNCTION
!!  Read/write self-energy on file.
!!
!! INPUTS
!!  self <type(self_type)>= variables related to self-energy
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  prtopt = flag for printing
!!  opt_rw = 0 (default) Set self-energy either to double counting or 0 depending on dmft_rslf
!!           1  Read Self-Energy.
!!           2  Write Self-Energy.
!!           3  Impose Self-Energy.
!!  istep_iter = iteration step
!!  opt_char = char to add at the end of the filename
!!  opt_imagonly = only reads the imaginary part (useful when reading from Maxent)
!!  opt_selflimit = 0th order moment of the self-energy (useful when reading from Maxent)
!!  opt_hdc = double counting (useful when reading from Maxent)
!!  opt_stop = stop when encountering errors
!!  opt_maxent = if > 0, read/write Maxent files
!!
!! OUTPUT
!!
!! SOURCE

subroutine rw_self(self,paw_dmft,prtopt,opt_rw,istep_iter,opt_char,opt_imagonly,&
                 & opt_selflimit,opt_hdc,opt_stop,opt_maxent)

!Arguments ------------------------------------
!type
 type(self_type), intent(inout) :: self
 type(paw_dmft_type), intent(inout) :: paw_dmft
 integer, intent(in) :: prtopt
 integer, optional, intent(in) :: istep_iter,opt_imagonly,opt_maxent,opt_rw,opt_stop
 character(len=4), optional, intent(in) :: opt_char
 type(matlu_type), optional, intent(inout) :: opt_selflimit(paw_dmft%natom)
 type(matlu_type), optional, intent(in) :: opt_hdc(paw_dmft%natom)
!local variables-------------------------------
 integer :: i,iall,iatom,iatu,icount,ier,iexist2,iexit,iflavor,ifreq,im,im1,ioerr
 integer :: ispinor,ispinor1,isppol,istep,istep_imp,istepiter,iter,iter_imp,lpawu,master
 integer :: myproc,natom,natom_read,ncount,ndim,ndim_read,nrecl,nspinor,nspinor_read
 integer :: nsppol,nsppol_read,nw_read,optmaxent,optrw,readimagonly,spacecomm,unitrot
 real(dp) :: fermie_read,x_r,x_i,xtemp
 logical :: lexist,lexist_rot,nondiaglevels,prtself
 character(len=30000) :: message ! Big buffer to avoid buffer overflow.
 character(len=fnlen) :: stringfile,tmpfil,tmpfil2,tmpfilrot,tmpmatrot
 character(len=1) :: tag_is
 character(len=3) :: self_iter
 character(len=4) :: chtemp
 character(len=5) :: tag_freq
 character(len=10) :: tag_at,tag_iflavor
 character(len=13) :: tag
 character(len=50) :: string_format
 type(oper_type) :: energy_level
 integer, allocatable :: unitselffunc_arr(:),unitselffunc_arr2(:),unitselfrot(:,:,:,:)
 real(dp), allocatable :: s_i(:,:),s_r(:,:) !,fermie_read2(:)
 complex(dpc), allocatable :: buffer(:)
 type(matlu_type), allocatable :: eigvectmatlu(:),level_diag(:),selfmomrot(:,:)
 type(oper_type), allocatable :: selfrotmatlu(:)
! *********************************************************************

 ABI_NVTX_START_RANGE(NVTX_DMFT_RW_SELF)
 natom   = paw_dmft%natom
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol
 !mbandc = paw_dmft%mbandc
 !nkpt = paw_dmft%nkpt

! Initialise spaceComm, myproc, and nproc
 istep = 0
 iter  = 0
 istep_imp = 0
 istepiter = 0
 iter_imp  = 0
 optmaxent = 0
 optrw = 0
 prtself = (paw_dmft%dmft_prtself == 1)
 readimagonly = 0

 if (present(opt_rw)) optrw = opt_rw
 if (present(opt_maxent)) optmaxent = opt_maxent

 if (present(opt_imagonly)) then
   if (opt_imagonly == 1 .and. paw_dmft%dmft_solv >= 5) then
     readimagonly = opt_imagonly
     write(message,*)
     write(message,'(a,4x,2a)') ch10,"About to read imaginary part of Self energy"
     call wrtout(std_out,message,'COLL')
   else
     write(message,'(4x,2a)') "About to read both real and imaginary part of Self energy"
     call wrtout(std_out,message,'COLL')
   end if ! opt_imagonly
 end if ! present(opt_imagonly)

 if (present(istep_iter)) istepiter = istep_iter

 if (paw_dmft%use_fixed_self > 0) then
   istep = istepiter / 1000
   iter  = istepiter - (istepiter/1000)*1000
   istep_imp = paw_dmft%use_fixed_self / 1000
   iter_imp  = paw_dmft%use_fixed_self - (paw_dmft%use_fixed_self/1000)*1000
 end if !use_fixed_self

 if (paw_dmft%dmft_rslf <= 0 .and. optrw == 1) optrw = 0

 iexist2   = 1
 iexit     = 0
 ioerr     = 0
 lexist    = .true.
 master    = 0
 myproc    = paw_dmft%myproc
 spacecomm = paw_dmft%spacecomm
 !nproc = paw_dmft%nproc

! write(std_out,*) "myproc,master",myproc,master
 !if(prtopt>200) then
 !endif

!   - For the Tentative rotation of the self-energy file (begin init)
 if (optmaxent > 0) then
   if (optrw == 2) then
     write(message,'(a,2x,a)') ch10," == About to print self-energy for MAXENT code in basis that diagonalizes the atomic levels"
   else if (optrw == 1)  then
     write(message,'(a,2x,a)') ch10," == About to read self-energy from MAXENT code"
   end if
   call wrtout(std_out,message,'COLL')

   ABI_MALLOC(eigvectmatlu,(natom))
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),eigvectmatlu(:))
 end if ! optmaxent > 0
!   - For the Tentative rotation of the self-energy file (end init)

!   - For the Tentative rotation of the self-energy file (begin diag)
 if (optrw == 2 .and. optmaxent > 0) then
   ABI_MALLOC(unitselfrot,(2*paw_dmft%maxlpawu+1,nspinor,nsppol,natom)) ! 7 is the max ndim possible
   ABI_MALLOC(level_diag,(natom))
   ABI_MALLOC(selfrotmatlu,(self%nw))
   call init_oper(paw_dmft,energy_level,opt_ksloc=2)
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),level_diag(:))
   call compute_levels(energy_level,self%hdc,paw_dmft,nondiag=nondiaglevels)
   write(tag,'(f13.5)') paw_dmft%fermie
   if (self%nw >= 2) then
     write(message,'(a,2x,2a)') ch10," == Print non Diagonalized Self Energy for Fermi Level= ",adjustl(tag)
     call wrtout(std_out,message,'COLL')
     call print_matlu(self%oper(2)%matlu(:),natom,1,compl=1,opt_exp=1)
   end if
   call diag_matlu(energy_level%matlu(:),level_diag(:),natom,prtopt,eigvectmatlu(:),test=paw_dmft%dmft_solv)
   write(message,'(a,2x,2a)') ch10," == Print Diagonalized levels for Fermi Level= ",adjustl(tag)
   call wrtout(std_out,message,'COLL')
   call print_matlu(level_diag(:),natom,1,compl=1,opt_exp=1)
   ! Rotate self
   do ifreq=1,self%nw
     call init_oper(paw_dmft,selfrotmatlu(ifreq),opt_ksloc=2)
     if (self%distrib%procf(ifreq) /= myproc) cycle
     call copy_matlu(self%oper(ifreq)%matlu(:),selfrotmatlu(ifreq)%matlu(:),natom)
     call rotate_matlu(selfrotmatlu(ifreq)%matlu(:),eigvectmatlu(:),natom,1)
   end do ! ifreq
   call gather_oper(selfrotmatlu(:),self%distrib,paw_dmft,2,master=master)
   if (self%has_moments == 1) then
     ABI_MALLOC(selfmomrot,(natom,self%nmoments))
     do i=1,self%nmoments
       call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),selfmomrot(:,i))
       call copy_matlu(self%moments(i)%matlu(:),selfmomrot(:,i),natom)
       call rotate_matlu(selfmomrot(:,i),eigvectmatlu(:),natom,1)
     end do ! i
   end if
   do ifreq=1,min(3,self%nw)
     write(tag_freq,'(i5)') ifreq
     if (ifreq < 3 .and. myproc == master) then ! very important to call print_matlu only on master node
       write(message,'(a,2x,2a)') ch10," == Print non Rotated Self Energy for frequency ",adjustl(tag_freq)
       call wrtout(std_out,message,'COLL')
       call print_matlu(self%oper(ifreq)%matlu(:),natom,1,compl=1)
       write(message,'(a,2x,2a)') ch10," == Print Rotated Self Energy for frequency ",adjustl(tag_freq)
       call wrtout(std_out,message,'COLL')
       call print_matlu(selfrotmatlu(ifreq)%matlu(:),natom,1,compl=1)
     else if (ifreq == 3) then
       write(message,'(a,2x,a,i4)') ch10,"  (Other frequencies not printed)"
       call wrtout(std_out,message,'COLL')
     end if ! ifreq<3
   end do ! ifreq
 end if ! optrw=2 and optmaxent>0
 !  Create file for rotation
 if (optmaxent > 0 .and. myproc == master .and. (optrw == 1 .or. optrw == 2)) then
   do iatom=1,natom
     lpawu = paw_dmft%lpawu(iatom)
     if (lpawu == -1) cycle
     call int2char4(iatom,tag_at)
     ABI_CHECK((tag_at(1:1)/='#'),'Bug: string length too short!')
     if (optrw == 2) then
       tmpmatrot = trim(paw_dmft%filapp)//'_UnitaryMatrix_iatom'//trim(tag_at)
     else if (optrw == 1) then
       tmpmatrot = trim(paw_dmft%filnamei)//'_UnitaryMatrix_iatom'//trim(tag_at)
       inquire(file=trim(tmpmatrot),exist=lexist_rot)
       if (.not. lexist_rot) ABI_ERROR("File "//trim(tmpmatrot)//" does not exist !")
     end if ! optrw
     unitrot = 3189 + iatom
#ifdef FC_NAG
     open(unit=unitrot,file=trim(tmpmatrot),status='unknown',form='formatted',recl=ABI_RECL)
#else
     open(unit=unitrot,file=trim(tmpmatrot),status='unknown',form='formatted')
#endif
     write(std_out,'(2a)') "      Open file ",trim(tmpmatrot)
     rewind(unitrot)
     ndim = (2*lpawu+1) * nspinor
     do isppol=1,nsppol
       do im=1,ndim
         do im1=1,ndim
           if (optrw == 2) then
             write(message,*) dble(eigvectmatlu(iatom)%mat(im,im1,isppol)),aimag(eigvectmatlu(iatom)%mat(im,im1,isppol))
             call wrtout(unitrot,message,'COLL')
           else if (optrw == 1) then
             read(unitrot,*) x_r,x_i
             eigvectmatlu(iatom)%mat(im,im1,isppol) = cmplx(x_r,x_i,kind=dp)
           end if ! optrw
         end do ! im1
       end do ! im
     end do ! isppol
     close(unitrot)
   end do ! iatom

   if (optrw == 1) then
     write(message,'(2a)') ch10," == Print non-rotated high-frequency limit of the self-energy"
     call wrtout(std_out,message,'COLL')
     call print_matlu(opt_selflimit(:),natom,1,compl=1)
     call rotate_matlu(opt_selflimit(:),eigvectmatlu(:),natom,1)
     write(message,'(2a)') ch10," == Print rotated high-frequency limit of the self-energy"
     call wrtout(std_out,message,'COLL')
     call print_matlu(opt_selflimit(:),natom,1,compl=1)
   end if ! optrw=1

 end if ! optmaxent > 0
!   - For the Tentative rotation of the self-energy file (end diag)

 if ((optrw == 2 .or. optrw == 1) .and. myproc == master) then
   ABI_MALLOC(unitselffunc_arr,(natom*nsppol))
   if (optrw == 2) then
     ABI_MALLOC(unitselffunc_arr2,(natom*nsppol))
     if (paw_dmft%idmftloop < 10) then
       write(self_iter,'("00",i1)') paw_dmft%idmftloop
     else if (paw_dmft%idmftloop >= 10 .and. paw_dmft%idmftloop < 100) then
       write(self_iter,'("0",i2)') paw_dmft%idmftloop
     else if (paw_dmft%idmftloop >= 100 .and. paw_dmft%idmftloop < 1000) then
       write(self_iter,'(i3)') paw_dmft%idmftloop
     else
       self_iter="xxx"
     end if ! idmftloop
   end if ! optrw=2
   iall = 0
   do iatom=1,natom
     lpawu = paw_dmft%lpawu(iatom)
     if (lpawu == -1) cycle
     ndim = 2*lpawu + 1
     ABI_MALLOC(s_r,(ndim*nspinor,ndim*nspinor))
     ABI_MALLOC(s_i,(ndim*nspinor,ndim*nspinor))
!       write(std_out,*) "print_self",ndim
     call int2char4(iatom,tag_at)
     ABI_CHECK((tag_at(1:1)/='#'),'Bug: string length too short!')
     do isppol=1,nsppol
       write(tag_is,'(i1)') isppol
!         do ispinor=1,nspinor
       iall = iall + 1
!           write(tag_is2,'(i1)')ispinor

!      ===========================
!      == Create name for file
!      ===========================

       if (self%w_type == "real") then
         tmpfil = trim(merge(paw_dmft%filnamei,paw_dmft%filapp,optrw==1))//'_Self_ra-omega_iatom'//trim(tag_at)//'_isppol'//tag_is
       else
         if (present(opt_char)) then
           tmpfil = trim(paw_dmft%filapp)//'Self_ra-omega_iatom'//trim(tag_at)//'_isppol'//tag_is//opt_char
         else
           stringfile = "_iatom" // trim(tag_at) // '_isppol' // tag_is
           if (optrw == 1 .and. paw_dmft%idmftloop == 0) then
             tmpfil = trim(paw_dmft%filselfin) // stringfile
             iexist2 = paw_dmft%ireadself
           else
             tmpfil = trim(paw_dmft%filapp) // '_Self-omega' // stringfile
           end if
         end if ! opt_char
       end if ! w_type
       if (optrw == 2) tmpfil2 = trim(tmpfil)//"_"//self_iter

       if (optrw == 1 .and. iexist2 == 1) then
         write(message,'(3a)') ch10,"  == Read self-energy and Fermi Level from file ",trim(tmpfil)
         call wrtout(std_out,message,'COLL')
       else if (optrw == 2) then
         write(message,'(3a)') ch10,"  == Write self-energy and Fermi Level on file ",trim(tmpfil)
         call wrtout(std_out,message,'COLL')
       end if
           !unitselffunc_arr(iall)=300+iall-1
       unitselffunc_arr(iall) = get_unit()
       ABI_CHECK(unitselffunc_arr(iall) > 0, "Cannot find free IO unit!")

       !- For the Tentative rotation of the self-energy file (create file)
       if (optrw == 2 .and. optmaxent > 0) then
         iflavor = 0
         do ispinor=1,nspinor
           do im=1,ndim
             iflavor = iflavor + 1
             call int2char4(iflavor,tag_iflavor)
             unitselfrot(im,ispinor,isppol,iatom) = 3000 + iflavor
             !ABI_CHECK(unitselfrot(im,ispinor,isppol,iatom) > 0, "Cannot find free IO unit for unitselfrot!")
             tmpfilrot = trim(paw_dmft%filapp)//'_Selfmxent'//&
                & trim(tag_at)//'_is'//tag_is//'_iflav'//trim(tag_iflavor)
             write(std_out,*) "Create file  ",trim(tmpfilrot)," unit ",unitselfrot(im,ispinor,isppol,iatom)," for flavor",iflavor
#ifdef FC_NAG
             open(unit=unitselfrot(im,ispinor,isppol,iatom),file=trim(tmpfilrot),status='unknown',form='formatted',recl=ABI_RECL)
#else
             open(unit=unitselfrot(im,ispinor,isppol,iatom),file=trim(tmpfilrot),status='unknown',form='formatted')
#endif
             rewind(unitselfrot(im,ispinor,isppol,iatom))
             write(unitselfrot(im,ispinor,isppol,iatom),'(3a)') "#  Diagonal component of the self-energy, in the basis that diagonalizes the electronic levels.", &
                                                         & ch10,"#       Frequency (Ha)              Real part              Imaginary part"

           end do ! im
         end do ! ispinor
       end if ! optrw=2 and optmaxent>0
       !- For the Tentative rotation of the self-energy file (create file)

!           write(std_out,*) "1"

!      ===========================
!      == Read: check that the file exists
!      ===========================
       if (optrw == 1) then
!           write(std_out,*) "3"
         inquire(file=trim(tmpfil),exist=lexist,recl=nrecl)
         if (.not. lexist .and. (paw_dmft%ireadself == 1 .or. paw_dmft%idmftloop > 0)) then
!           write(std_out,*) "4"
           iexist2 = 0
           write(message,'(4x,a,i5,3a)') "File number",unitselffunc_arr(iall),&
             & " called ",trim(tmpfil)," does not exist"
!               write(std_out,*) lexist,nrecl
           call wrtout(std_out,message,'COLL')
         end if ! not lexist
       end if ! optrw=1
           !write(std_out,*) "2"

!      ===========================
!      == Open file
!      ===========================
       if (optrw == 2 .or. (optrw == 1 .and. iexist2 == 1)) then
             !write(std_out,*) "5"
#ifdef FC_NAG
         open(unit=unitselffunc_arr(iall),file=trim(tmpfil),status='unknown',form='formatted',recl=ABI_RECL)
#else
         open(unit=unitselffunc_arr(iall),file=trim(tmpfil),status='unknown',form='formatted')
#endif
         rewind(unitselffunc_arr(iall))

         if (optrw == 2 .and. prtself) then
           unitselffunc_arr2(iall) = get_unit()
           ABI_CHECK(unitselffunc_arr2(iall) > 0,"Cannot find free IO unit!")
#ifdef FC_NAG
           open(unit=unitselffunc_arr2(iall),file=trim(tmpfil2),status="unknown",form="formatted",recl=ABI_RECL)
#else
           open(unit=unitselffunc_arr2(iall),file=trim(tmpfil2),status="unknown",form="formatted")
#endif
           rewind(unitselffunc_arr2(iall))
         end if ! optrw=2

             !write(std_out,*) "61",nrecl
         if (prtopt >= 3) then
           write(message,'(3a,i4)') '     Opened file : ',trim(tmpfil),' on unit ',unitselffunc_arr(iall)
           call wrtout(std_out,message,'COLL')
         end if ! prtopt>=3
       end if
           !write(std_out,*) "6",nrecl

!      ===========================
!      == Check Header
!      ===========================

       if (optrw == 2) then

         write(message,'(11a,4i5,i6,2x,e25.17e3)') "# DFT+DMFT self-energy for each frequency",ch10, &
                           & "# Columns are ordered this way:",ch10, &
                           & "# Frequency (Ha)  ((((Re(Sigma(im,im1,ispinor,ispinor1)) ", &
                           & "Im(Sigma(im,im1,ispinor,ispinor1)),im=1,2*l+1),im1=1,2*l+1),", &
                           & "ispinor=1,nspinor),ispinor1=1,nspinor) where the leftmost index varies first", &
                           & ch10,"# natom,nsppol,nspinor,ndim,nw,fermilevel",ch10,&
                           & "####",natom,nsppol,nspinor,ndim,self%nw,paw_dmft%fermie
         call wrtout(unitselffunc_arr(iall),message,'COLL')
         if (prtself) then
           call wrtout(unitselffunc_arr2(iall),message,'COLL')
         end if
       else if (optrw == 1 .and. iexist2 == 1 .and. readimagonly == 0) then
         read(unitselffunc_arr(iall),*)
         read(unitselffunc_arr(iall),*)
         read(unitselffunc_arr(iall),*)
         read(unitselffunc_arr(iall),*)
         read(unitselffunc_arr(iall),*,iostat=ioerr) &
           & chtemp,natom_read,nsppol_read,nspinor_read,ndim_read,nw_read,fermie_read
             !if(ioerr<0) then
!              write(std_out,*)" HEADER IOERR"
!              write(std_out,'(a4,2x,31(e15.8,2x))') chtemp,natom_read,nsppol_read,nspinor_read,ndim_read,nw_read,fermie_read
             !endif
         if (ioerr == 0) then
           write(message,'(a,3x,3a,i12,2a,i11,2a,i10,2a,i13,2a,i15,2a,e25.8)') ch10,"Data in Self Energy file corresponds to",&
              & ch10,"     natom",natom_read,&
              & ch10,"     nsppol",nsppol_read,&
              & ch10,"     nspinor",nspinor_read,&
              & ch10,"     ndim",ndim_read, &
              & ch10,"     nw",nw_read, &
              & ch10,"     Fermi level",fermie_read
           call wrtout(std_out,message,'COLL')
           if ((natom /= natom_read) .or. (nsppol_read /= nsppol) .or. &
             & (nspinor /= nspinor_read) .or. (nw_read /= self%nw)) then
             write(message,'(a,3x,3a,i12,2a,i11,2a,i10,2a,i13,2a,i15,2a,e25.8)') ch10,"Data required is ",&
                & ch10,"     natom",natom,&
                & ch10,"     nsppol",nsppol,&
                & ch10,"     nspinor",nspinor,&
                & ch10,"     ndim",ndim, &
                & ch10,"     nw",self%nw, &
                & ch10,"     Fermi level",paw_dmft%fermie
             call wrtout(std_out,message,'COLL')
             message = "Dimensions in self-energy file are not correct"
             if (present(opt_stop)) then
               ABI_ERROR(message)
             else
               ABI_WARNING(message)
             end if
             iexist2 = 2
           end if
         else
           ABI_WARNING("Self-energy file is empty")
         end if ! ioerr
       end if ! optrw
           !write(std_out,*) "7"

!      ===========================
!      == Write/Read self in the file
!      ===========================

       !rewind(111)
       do ifreq=1,self%nw
         if (optrw == 2) then
!               write(std_out,'(a,2x,31(e15.8,2x))') &
!&              "SETEST",self%omega(ifreq),&
!&              (self%oper(ifreq)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)&
!&               ,im=1,ndim)
!               write(std_out,*) self%omega(ifreq),&
!&              ((self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor)&
!&               ,im=1,ndim),im1=1,ndim)
!               write(message,'(2x,393(e25.17,2x))')  self%omega(ifreq),&
!&              ((((self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)&
!&              ,im=1,ndim),im1=1,ndim),ispinor=1,nspinor),ispinor1=1,nspinor)

               !MGNAG Runtime Error: wrtout_cpp.f90, line 896: Buffer overflow on output
               !Is it possible to rewrite the code below to avoid such a long message
               !What about Netcdf binary files ?
           string_format = merge('(2x,393(e25.17e3,2x))','(2x,393(e18.10e3,2x))',nspinor==1)
           write(message,string_format) self%omega(ifreq),&
               & ((((dble(self%oper(ifreq)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,isppol)),&
               & aimag(self%oper(ifreq)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,isppol)),&
               & im=1,ndim),im1=1,ndim),ispinor=1,nspinor),ispinor1=1,nspinor)
           call wrtout(unitselffunc_arr(iall),message,'COLL')
           if (prtself) then
             call wrtout(unitselffunc_arr2(iall),message,'COLL')
           end if

           !- For the Tentative rotation of the self-energy file (begin rot)
           !----------------------------------------------------------------
           if (optmaxent > 0) then
             !iflavor = 0
             do ispinor=1,nspinor
               do im=1,ndim
                 !iflavor = iflavor + 1
                    ! if(ifreq<5) then
                    !   write(std_out,*) "Write in file unit",unitselfrot(iatom,isppol,ispinor,im),"for flavor",iflavor
                    ! endif
                 write(message,'(2x,393(es24.16e3,2x))') self%omega(ifreq),&
                  & dble(selfrotmatlu(ifreq)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im+(ispinor-1)*ndim,isppol)),&
                  & aimag(selfrotmatlu(ifreq)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im+(ispinor-1)*ndim,isppol))
!                     write(6,'(2x,393(e18.10,2x))')  self%omega(ifreq),&
!&                      real(selfrotmatlu(iatom)%mat(im,im,isppol,ispinor,ispinor)),&
!&                      aimag(selfrotmatlu(iatom)%mat(im,im,isppol,ispinor,ispinor))
                    ! if(iflavor==1) then
                    ! write(1024,*) iatom,isppol,ispinor,im,unitselfrot(iatom,isppol,ispinor,im)
                    ! write(1024,'(2x,393(e18.10,2x))')  self%omega(ifreq),&
                    ! &  real(selfrotmatlu(iatom)%mat(im,im,isppol,ispinor,ispinor)),&
                    ! &  aimag(selfrotmatlu(iatom)%mat(im,im,isppol,ispinor,ispinor))
                    ! endif
                 call wrtout(unitselfrot(im,ispinor,isppol,iatom),message,'COLL')
               end do ! im
             end do ! ispinor
           end if ! optrw=2 and optmaxent>0
               !- For the Tentative rotation of the self-energy file (end rot)

!               write(std_out,*) unitselffunc_arr(iall)
         else if (optrw == 1 .and. iexist2 == 1 .and. ioerr == 0 .and. readimagonly == 0) then
           !write(std_out,*) "8"
!               read(unitselffunc_arr(iall),'(2x,31(e15.8,2x))',iostat=ioerr) &
!&              xtemp,(s_r(im),s_i(im),im=1,ndim)
           !if (readimagonly == 0) then
           read(unitselffunc_arr(iall),*,iostat=ioerr) xtemp,&
               & ((((s_r(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim),s_i(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim), &
               & im=1,ndim),im1=1,ndim),ispinor=1,nspinor),ispinor1=1,nspinor)
!               if(ioerr<0) then
!                write(std_out,*)" SELF IOERR<"
!               else if(ioerr>0) then
!                write(std_out,*)" SELF IOERR>"
!                write(std_out,'(a4,2x,31(e15.8,2x))') xtemp,(s_r(im),s_i(im),im=1,ndim)
!               endif
           self%oper(ifreq)%matlu(iatom)%mat(:,:,isppol) = cmplx(s_r(:,:),s_i(:,:),kind=dp)
         end if ! optrw
       end do ! ifreq

       !- For the Tentative rotation of the self-energy file (begin close file)
       if (optrw == 2 .and. optmaxent > 0) then
         do ispinor=1,nspinor
           do im=1,ndim
             if (self%has_moments == 1) then
               do i=1,self%nmoments
                 write(tag_is,'(i1)') i
                 write(message,'(a,2x,2(es24.16e3,2x))') "#moments_"//tag_is, &
                   & dble(selfmomrot(iatom,i)%mat(im+(ispinor-1)*ndim,im+(ispinor-1)*ndim,isppol)), &
                   & aimag(selfmomrot(iatom,i)%mat(im+(ispinor-1)*ndim,im+(ispinor-1)*ndim,isppol))
                 call wrtout(unitselfrot(im,ispinor,isppol,iatom),message,'COLL')
               end do ! i
             end if ! nmoments
             close(unitselfrot(im,ispinor,isppol,iatom))
             write(std_out,*) "Close file unit",unitselfrot(im,ispinor,isppol,iatom)
           end do ! im
         end do ! ispinor
       end if ! optrw=2 and optmaxent>0

       !- For the Tentative rotation of the self-energy file (end close file)

       if (optrw == 1 .and. iexist2 == 1 .and. ioerr == 0 .and. readimagonly == 1) then

         ! Read self energy from Maxent (imag part) on the real axis
         !----------------------------------------------------------
         do ifreq=1,self%nw
           self%oper(ifreq)%matlu(iatom)%mat(:,:,:) = czero
         end do ! ifreq
         do ispinor=1,nspinor
           do im=1,ndim
             do ifreq=1,self%nw
               read(unitselffunc_arr(iall),*,iostat=ioerr) xtemp,s_i(im+(ispinor-1)*ndim,im+(ispinor-1)*ndim)
               ! minus sign because - Im Sigma is the output of OmegaMaxent
               self%oper(ifreq)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im+(ispinor-1)*ndim,isppol) &
                 & = cmplx(zero,-half*s_i(im+(ispinor-1)*ndim,im+(ispinor-1)*ndim),kind=dp)
             end do ! ifreq
           end do ! im
         end do ! ispinor

         write(message,'(4x,2a)') " Read only diagonal self energy from Maxent"
         call wrtout(std_out,message,'COLL')
               !write(6,*) "opt_hdc",opt_hdc(1)%mat(1,1,1,1,1)

       end if ! optrw=1

!      ===========================
!      == Write/Read hdc in the file
!      ===========================
       if (optrw == 2) then
!             write(std_out,'(a,2x,31(e15.8,2x))') &
!&            "SETEST #dc ",(self%hdc%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor),im=1,ndim)
         write(message,'(a,2x,500(e25.17e3,2x))') &
          & "#dc ",((((self%hdc%matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,isppol),&
            & im=1,ndim),im1=1,ndim),ispinor=1,nspinor),ispinor1=1,nspinor)
         call wrtout(unitselffunc_arr(iall),message,'COLL')
         if (prtself) then
           call wrtout(unitselffunc_arr2(iall),message,'COLL')
         end if
         if (self%has_moments == 1) then
           do i=1,self%nmoments
             write(tag_is,'(i1)') i
             write(message,'(a,2x,500(e25.17e3,2x))') "#moments_"//trim(tag_is),&
                & ((self%moments(i)%matlu(iatom)%mat(im,im1,isppol),im=1,nspinor*ndim),im1=1,nspinor*ndim)
             call wrtout(unitselffunc_arr(iall),message,'COLL')
             if (prtself) then
               call wrtout(unitselffunc_arr2(iall),message,'COLL')
             end if
           end do ! i
         end if ! moments
       else if (optrw == 1 .and. iexist2 == 1 .and. ioerr == 0 .and. readimagonly == 0) then
         !write(std_out,*) "8"
         read(unitselffunc_arr(iall),*,iostat=ioerr) &
           & chtemp,((((s_r(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim),s_i(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim),&
              & im=1,ndim),im1=1,ndim),ispinor=1,nspinor),ispinor1=1,nspinor)
            !if(ioerr<0) then
!              write(std_out,*)" HDC IOERR<",ioerr
             !else if(ioerr>0) then
!              write(std_out,*)" HDC IOERR>",ioerr
             !endif
         self%hdc%matlu(iatom)%mat(:,:,isppol) = cmplx(s_r(:,:),s_i(:,:),kind=dp)

         if (self%has_moments == 1) then
           do i=1,self%nmoments
             read(unitselffunc_arr(iall),*,iostat=ioerr) &
               & chtemp,((s_r(im,im1),s_i(im,im1),im=1,nspinor*ndim),im1=1,nspinor*ndim)
             self%moments(i)%matlu(iatom)%mat(:,:,isppol) = cmplx(s_r(:,:),s_i(:,:),kind=dp)
           end do ! i
         end if ! moments
             !write(6,*) "read selfhdc",self%hdc%matlu(1)%mat(1,1,1,1,1)
       else if (readimagonly == 1 .and. (.not. present(opt_hdc))) then
         self%hdc%matlu(iatom)%mat(:,:,isppol) = czero
       !else
       !  write(std_out,*) "     self%hdc fixed in kramerskronig_self"
       end if ! optrw
       close(unitselffunc_arr(iall))
       if (optrw == 2 .and. prtself) close(unitselffunc_arr2(iall))
!         enddo ! ispinor
     end do ! isppol
     ABI_FREE(s_r)
     ABI_FREE(s_i)
   end do ! iatom

   ABI_FREE(unitselffunc_arr)
   ABI_SFREE(unitselffunc_arr2)
 end if ! optrw==1 or 2 and myproc==master

!  ===========================
!  == Error messages
!  ===========================
 if (optrw == 1) then
!   call xmpi_barrier(spacecomm)
   !write(std_out,*) ncount,maxval(pawtab(:)%lpawu)*2+1
   call xmpi_bcast(iexist2,master,spacecomm,ier)
   call xmpi_bcast(ioerr,master,spacecomm,ier)
   if (iexist2 == 0 .or. ioerr /= 0) then
     message = "Self-energy file does not exist or is incomplete"
     if (readimagonly == 1 .or. present(opt_stop)) then
       if (readimagonly == 1) message = "Self-energy file does not exist or is incomplete: check the number of self-energy data in file"
       ABI_ERROR(message)
     else if (paw_dmft%ireadself == 1 .or. paw_dmft%idmftloop > 0) then
       ABI_WARNING(message)
     end if ! readimagonly=1 or present(opt_stop)
     if (iexist2 == 0 .and. (paw_dmft%ireadself == 1 .or. paw_dmft%idmftloop > 0)) then
       write(message,'(4x,2a)') "File does not exist"
       call wrtout(std_out,message,'COLL')
     end if ! iexist2=0
     if (ioerr < 0) then
       write(message,'(4x,2a)') "End of file reached"
       call wrtout(std_out,message,'COLL')
     end if ! ioerr<0
     if (ioerr > 0) then
       write(message,'(4x,2a)') "Error during read statement"
       call wrtout(std_out,message,'COLL')
     end if ! ioerr>0
     !if (paw_dmft%dmft_solv /= 4) then
     write(message,'(4x,2a,5i5,2x,e14.7)') "-> Set Self-Energy Equal to double counting term"
     !else if (paw_dmft%dmft_solv == 4) then
     !  write(message,'(4x,a,a,5i5,2x,e14.7)') "-> Put Self-Energy Equal to dc term - shift"
     !  call wrtout(std_out,message,'COLL')
     !  write(message,'(4x,a,a,5i5,2x,e14.7)') " No self energy is given, change dmft_rslf"
     !  ABI_ERROR(message)
     !end if
     call wrtout(std_out,message,'COLL')
     do ifreq=1,self%nw
!       write(std_out,*) "before",self%oper(1)%matlu(1)%mat(1,1,1,1,1)
!       write(std_out,*) "before",self%hdc%matlu(1)%mat(1,1,1,1,1)
       call copy_matlu(self%hdc%matlu(:),self%oper(ifreq)%matlu(:),natom)
       if (nspinor == 1 .and. nsppol == 2 .and. paw_dmft%dmft_dc >= 5) then
         do iatom=1,natom
           lpawu = paw_dmft%lpawu(iatom)
           if (lpawu == -1) cycle
           ndim = 2*lpawu + 1
           do im=1,ndim
             self%oper(ifreq)%matlu(iatom)%mat(im,im,1) = self%oper(ifreq)%matlu(iatom)%mat(im,im,2) + &
                 & paw_dmft%dmft_shiftself(iatom)
           end do ! im
         end do ! iatom
       end if ! nspinor=1 and nsppol=2
!       write(std_out,*) "after",self%oper(1)%matlu(1)%mat(1,1,1,1,1)
!       write(std_out,*) "before",self%hdc%matlu(1)%mat(1,1,1,1,1)
       !if (paw_dmft%dmft_solv == 4) then
!         if(ifreq==1) write(std_out,*) "shift",self%qmc_shift(1)
       !  call shift_matlu(self%oper(ifreq)%matlu(:),natom,cmplx(self%qmc_shift(:),zero,kind=dp),1)
!         if(ifreq==1) write(std_out,*) "self after dc and shift",self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)
!         if(ifreq==1) write(std_out,*) "shift",self%qmc_shift(1)
       !end if ! dmft_solv=4
     end do ! ifreq
     if (self%has_moments == 1) then
       call copy_matlu(self%oper(1)%matlu(:),self%moments(1)%matlu(:),natom)
     end if
   else ! test read successfull
!   call xmpi_barrier(spacecomm)
!! lignes 924-928 semblent inutiles puisque la valeur de paw_dmft%fermie creee
!! en ligne 927 est ecrasee en ligne 992. BA+jmb
!!     ABI_MALLOC(fermie_read2,(1))
!!     fermie_read2(1)=fermie_read
!!     call xmpi_sum(fermie_read2,spacecomm ,ier)
!!     paw_dmft%fermie=fermie_read2(1)
!!     ABI_FREE(fermie_read2)
     !ncount = natom * nsppol * (nspinor**2) * (self%nw+1) *(maxval(paw_dmft%lpawu(:))*2+1)**2
     ncount = 0
     do iatom=1,natom
       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       ndim = 2*lpawu + 1
       ncount = ncount + ndim**2
     end do ! iatom
     ncount = ncount * (nspinor**2) * (self%nw+self%nmoments+1) * nsppol

!    ===========================
!     bcast to other proc
!    ===========================
     ABI_MALLOC(buffer,(ncount))
!! BA+jmb
     !fermie_read2 = zero
   !write(std_out,*) self%nw
     if (myproc == master) then

!      == Send read data to all process
       if (readimagonly == 0) paw_dmft%fermie = fermie_read
       icount = 0
!      Self energy-----------
       do ifreq=1,self%nw
         do iatom=1,natom
           lpawu = paw_dmft%lpawu(iatom)
           if (lpawu == -1) cycle
           ndim = (2*lpawu+1) * nspinor
           do isppol=1,nsppol
             do im1=1,ndim
               buffer(icount+1:icount+ndim) = self%oper(ifreq)%matlu(iatom)%mat(:,im1,isppol)
               icount = icount + ndim
                 !if (icount > ncount) then
                 !  write(message,'(2a,2i5)') ch10,"Error buffer",icount,ncount
                 !  iexit = 1
                 !  ABI_ERROR(message)
                 !end if ! icount > ncount
             end do ! im1
           end do ! isppol
         end do ! iatom
       end do ! ifreq
!      Double counting-------
       do iatom=1,natom
         lpawu = paw_dmft%lpawu(iatom)
         if (lpawu == -1) cycle
         ndim = (2*lpawu+1) * nspinor
         do isppol=1,nsppol
           do im1=1,ndim
             buffer(icount+1:icount+ndim) = self%hdc%matlu(iatom)%mat(:,im1,isppol)
             icount = icount + ndim
             !if (icount > ncount) then
             !  write(message,'(2a,2i5)') ch10,"Error buffer",icount,ncount
             !  iexit = 1
             !  ABI_ERROR(message)
             !end if ! icount > ncount
           end do ! im1
         end do ! isppol
       end do ! iatom
       if (self%has_moments == 1) then
         do i=1,self%nmoments
           do iatom=1,natom
             lpawu = paw_dmft%lpawu(iatom)
             if (lpawu == -1) cycle
             ndim = (2*lpawu+1) * nspinor
             do isppol=1,nsppol
               do im1=1,ndim
                 buffer(icount+1:icount+ndim) = self%moments(i)%matlu(iatom)%mat(:,im1,isppol)
                 icount = icount + ndim
               end do ! im1
             end do ! isppol
           end do ! iatom
         end do ! i
       end if ! moments
     end if ! proc=master
     call xmpi_bcast(buffer(:),master,spacecomm,ier)
!    call xmpi_sum(iexit,spacecomm ,ier)
!!JB call xmpi_barrier(spacecomm)
     !call xmpi_sum(buffer,spacecomm,ier)
!!JB call xmpi_barrier(spacecomm)

! bcast fermi level
     !call xmpi_sum(fermie_read2,spacecomm,ier)
     if (readimagonly == 0) then
       call xmpi_bcast(paw_dmft%fermie,master,spacecomm,ier)
     end if

     if (ier /= 0) then
       message =  "error in xmpi_sum in rw_self"
       ABI_ERROR(message)
     end if
     !paw_dmft%fermie = fermie_read2(1)
!     write(std_out,*) "Fermi level",paw_dmft%fermie
     icount = 0
!    Self ---------------
     do ifreq=1,self%nw
       do iatom=1,natom
         lpawu = paw_dmft%lpawu(iatom)
         if (lpawu == -1) cycle
         ndim = (2*lpawu+1) * nspinor
         do isppol=1,nsppol
           do im1=1,ndim
             self%oper(ifreq)%matlu(iatom)%mat(:,im1,isppol) = buffer(icount+1:icount+ndim)
             icount = icount + ndim
                     !write(6,*)'self procs', ifreq, self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)
           end do ! im1
         end do ! isppol
       end do ! iatom
     end do ! ifreq
!    hdc  ---------------
     do iatom=1,natom
       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       ndim = (2*lpawu+1) * nspinor
       do isppol=1,nsppol
         do im1=1,ndim
           self%hdc%matlu(iatom)%mat(:,im1,isppol) = buffer(icount+1:icount+ndim)
           icount = icount + ndim
         end do ! im1
       end do ! isppol
     end do ! iatom
     if (self%has_moments == 1) then
       do i=1,self%nmoments
         do iatom=1,natom
           lpawu = paw_dmft%lpawu(iatom)
           if (lpawu == -1) cycle
           ndim = (2*lpawu+1) * nspinor
           do isppol=1,nsppol
             do im1=1,ndim
               self%moments(i)%matlu(iatom)%mat(:,im1,isppol) = buffer(icount+1:icount+ndim)
               icount = icount + ndim
             end do ! im1
           end do ! isppol
         end do ! iatom
       end do ! i
     end if ! moments
     !ABI_FREE(fermie_read2)
     ABI_FREE(buffer)
   end if  ! test read successful
 end if  ! optrw==1

 if (optmaxent > 0 .and. optrw == 1 .and. iexist2 == 1 .and. ioerr == 0 .and. readimagonly == 1) then

   write(message,'(4x,2a)') " Rotate Back self-energy in the cubic basis"
   call wrtout(std_out,message,'COLL')

   ! Kramers Kronig
   !-----------------------------

   call kramerskronig_self(self,opt_selflimit(:),opt_hdc(:),paw_dmft,paw_dmft%filapp)

   call xmpi_matlu(eigvectmatlu(:),natom,spacecomm,master=master,option=2)

   call copy_matlu(opt_hdc(:),self%hdc%matlu(:),natom)

   ! Rotate back self
   !-----------------------------

   do ifreq=1,self%nw
     if (ifreq < 20) then
       write(tag_freq,'(i5)') ifreq
       write(message,'(a,2x,2a)') ch10," == Print Rotated Self Energy on real axis for frequency ",adjustl(tag_freq)
       call wrtout(std_out,message,'COLL')
       call print_matlu(self%oper(ifreq)%matlu(:),natom,1,compl=1)
     end if ! ifreq<20
     if (self%distrib%procf(ifreq) /= myproc) cycle
     call rotate_matlu(self%oper(ifreq)%matlu(:),eigvectmatlu(:),natom,-1)
   end do ! ifreq
   call gather_oper(self%oper(:),self%distrib,paw_dmft,opt_ksloc=2)
   do ifreq=1,min(19,self%nw)
     write(tag_freq,'(i5)') ifreq
     write(message,'(a,2x,2a)') ch10," == Print Self Energy rotated back in cubic basis on real axis for frequency ",adjustl(tag_freq)
     call wrtout(std_out,message,'COLL')
     call print_matlu(self%oper(ifreq)%matlu(:),natom,1,compl=1)
   end do ! ifreq
 end if ! rotate back self-energy

! call xmpi_barrier(spacecomm)
           !write(std_out,*) "9"
!   - For the Tentative rotation of the self-energy file (begin destroy)
 if (optmaxent > 0) then

   if (optrw == 2) then
     call destroy_oper(energy_level)
     call destroy_matlu(level_diag(:),natom)
     if (self%has_moments == 1) then
       do i=1,self%nmoments
         call destroy_matlu(selfmomrot(:,i),natom)
       end do
       ABI_FREE(selfmomrot)
     end if
     do ifreq=1,self%nw
       call destroy_oper(selfrotmatlu(ifreq))
     end do ! ifreq
     ABI_FREE(level_diag)
     ABI_FREE(selfrotmatlu)
     ABI_FREE(unitselfrot) ! 7 is the max ndim possible
   end if ! optrw=2

   call destroy_matlu(eigvectmatlu(:),natom)
   ABI_FREE(eigvectmatlu)
 end if ! optmaxent > 0
!   - For the Tentative rotation of the self-energy file (end destroy)


!   call flush_unit(std_out)
!   ABI_ERROR("Aboring now")
 if (optrw == 0) then
   if (paw_dmft%dmft_rslf == 0) then
     !if (paw_dmft%dmft_solv /= 4) then
     write(message,'(4x,a)') "-> Set Self-Energy Equal to double counting term"
     !else if (paw_dmft%dmft_solv == 4) then
     !  write(message,'(4x,a,a,5i5,2x,e14.7)') "-> Put Self-Energy Equal to dc term - shift"
     !end if ! dmft_solv=4
   else if (paw_dmft%dmft_rslf == -1) then
     write(message,'(4x,a)') "-> Set Self-Energy Equal to zero"
   end if ! dmft_rslf=0
   call wrtout(std_out,message,'COLL')
   do ifreq=1,self%nw
     if (paw_dmft%dmft_rslf == 0) then
       call copy_matlu(self%hdc%matlu(:),self%oper(ifreq)%matlu(:),natom)
      ! if(ifreq==1) write(std_out,*) "self after dc",self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)
       !if (paw_dmft%dmft_solv == 4) then
      !   if(ifreq==1) write(std_out,*) "shift",self%qmc_shift(1)
       !  call shift_matlu(self%oper(ifreq)%matlu(:),natom,cmplx(self%qmc_shift(:),zero,kind=dp),1)
      !   if(ifreq==1) write(std_out,*) "self after dc and shift",self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)
      !   if(ifreq==1) write(std_out,*) "shift",self%qmc_shift(1)
       !end if ! dmft_solv=4
     else if (paw_dmft%dmft_rslf == -1) then
       call zero_matlu(self%oper(ifreq)%matlu(:),natom)
     end if ! dmft_rslf
   end do ! ifreq
   if (paw_dmft%dmft_rslf == 0 .and. self%has_moments == 1) then
     call copy_matlu(self%hdc%matlu(:),self%moments(1)%matlu(:),natom)
   end if
 end if ! optrw=0

! write(std_out,*) "optrw,use_fixed_self,istep,iter,istep_imp,iter_imp"
! write(std_out,*) optrw,paw_dmft%use_fixed_self,istep,iter,istep_imp,iter_imp
 if ((optrw == 1 .or. optrw == 3) .and. paw_dmft%use_fixed_self > 0 .and. istep <= istep_imp .and. iter <= iter_imp) then
   write(message,'(4x,a)') "-> Set Self-Energy Equal to imposed self-energy"
   call wrtout(std_out,message,'COLL')
   do ifreq=1,self%nw
     iatu = 0
     do iatom=1,natom
       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       self%oper(ifreq)%matlu(iatom)%mat(:,:,:) = czero
       iatu = iatu + 1
       ndim = 2*lpawu + 1
       do isppol=1,nsppol
         do ispinor=1,nspinor
           do im=1,ndim
             if (nspinor == 1) then
               self%oper(ifreq)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im+(ispinor-1)*ndim,isppol) = paw_dmft%fixed_self(im,im,isppol,iatu)
!            write(std_out,*) paw_dmft%fixed_self(im,im,isppol,iatu)
             else
               self%oper(ifreq)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im+(ispinor-1)*ndim,isppol) = paw_dmft%fixed_self(im,im,ispinor,iatu)
!                 write(message,'(a,i4,i4,2x,e20.10)') " Fixed self not implemented for nspinor==2"
!                 call wrtout(std_out,  message,'COLL')
!                 ABI_ERROR("Aboring now")
             end if ! nspinor
           end do ! im
         end do ! ispinor
       end do ! isppol
     end do ! iatom
   end do ! ifreq
   if (self%has_moments == 1) then
     call copy_matlu(self%oper(1)%matlu(:),self%moments(1)%matlu(:),natom)
     do i=2,self%nmoments
       call zero_matlu(self%moments(i)%matlu(:),natom)
     end do ! i
   end if

 end if ! use_fixed_self
 ABI_NVTX_END_RANGE()

end subroutine rw_self
!!***

!!****f* m_self/new_self
!! NAME
!! new_self
!!
!! FUNCTION
!!
!!  Mix Old and New self_energy with the mixing coefficient dmft_mxsf
!!
!! INPUTS
!!  self <type(self_type)>= variables related to self-energy
!!  self_new <type(self_type)>= variables related to the new self-energy
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! OUTPUT
!!  self <type(self_type)>= variables related to mixed self-energy
!!
!! SOURCE

subroutine new_self(self,self_new,paw_dmft)

!Arguments ------------------------------------
!type
! type(crystal_t),intent(in) :: cryst_struc
 type(self_type), intent(inout) :: self
 type(self_type), intent(in) :: self_new
 type(paw_dmft_type), intent(in) :: paw_dmft
!Local variables-------------------------------
 integer :: i,iatom,icount,ifreq,im,im1,isppol,lpawu,natom,ndim,nspinor,nsppol
 real(dp) :: alpha,diff_self,sum_self
 character(len=500) :: message
! *********************************************************************

 alpha   = paw_dmft%dmft_mxsf
 natom   = paw_dmft%natom
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol

 do ifreq=1,self%nw
   do iatom=1,natom
     lpawu = paw_dmft%lpawu(iatom)
     if (lpawu == -1) cycle
     self%oper(ifreq)%matlu(iatom)%mat(:,:,:) = (one-alpha)*self%oper(ifreq)%matlu(iatom)%mat(:,:,:) + &
         & alpha*self_new%oper(ifreq)%matlu(iatom)%mat(:,:,:)
       !  warning: self_new is the recent self-energy, which is mixed with self
       !  to give self= mixed self energy. self_new is deallocated just after.
   end do ! iatom
 end do ! ifreq

 diff_self = zero
 sum_self = zero
 icount = 0

 do iatom=1,natom
   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   ndim = nspinor * (2*lpawu+1)
   icount = icount + nsppol*ndim
   do isppol=1,nsppol
     do im1=1,ndim
       do im=1,ndim
         diff_self = diff_self + abs(self%hdc%matlu(iatom)%mat(im,im1,isppol)-self_new%hdc%matlu(iatom)%mat(im,im1,isppol))
         sum_self = sum_self + abs(self%hdc%matlu(iatom)%mat(im,im1,isppol))
         self%hdc%matlu(iatom)%mat(im,im1,isppol) = (one-alpha)*self%hdc%matlu(iatom)%mat(im,im1,isppol) + &
           & alpha*self_new%hdc%matlu(iatom)%mat(im,im1,isppol)
       end do ! im
     end do ! im1
   end do ! isppol
 end do ! iatom

 if (self%has_moments == 1) then
   do i=1,self%nmoments
     do iatom=1,natom
       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       self%moments(i)%matlu(iatom)%mat(:,:,:) = (one-alpha)*self%moments(i)%matlu(iatom)%mat(:,:,:) + &
           & alpha*self_new%moments(i)%matlu(iatom)%mat(:,:,:)
     end do ! iatom
   end do ! i
 end if ! moments

 !if(opt_mix==1) then
 !endif

 diff_self = diff_self / dble(icount)

 write(message,'(8x,a,e12.5)') "DMFT Loop: Precision on self-energy is",diff_self
 call wrtout(std_out,message,'COLL')
 if (diff_self < paw_dmft%dmft_fermi_prec .and. sum_self > tol6 .and. paw_dmft%idmftloop >= 2) then
    write(message,'(a,8x,a,e9.2,a,8x,a)') ch10, "Change of self =<", paw_dmft%dmft_fermi_prec,&
      & ch10,"DMFT Loop: Self Energy is converged"
    call wrtout(std_out,message,'COLL')
    self%iself_cv = 1
 else
    write(message,'(a,8x,a)') ch10,"DMFT Loop: Self Energy is not converged"
    call wrtout(std_out,message,'COLL')
    self%iself_cv = 0
 end if ! diff_self

end subroutine new_self
!!***

!!****f* m_self/make_qmcshift_self
!! NAME
!! make_qmcshift_hu
!!
!! FUNCTION
!!
!! INPUTS
!!  hu <type(hu_type)> = U interaction
!!  paw_dmft  <type(paw_dmft_type)> = paw+dmft related data
!!
!! OUTPUT
!!  self%qmc_shift in self <type(self_type)> = Self-energy
!!
!! SOURCE

!subroutine make_qmcshift_self(cryst_struc,hu,self,apply)

!Arguments ------------------------------------
!type
! type(crystal_t),intent(in) :: cryst_struc
! type(hu_type),intent(in) :: hu(cryst_struc%ntypat)
! type(self_type),intent(inout) :: self
! logical, optional :: apply

!Local variables-------------------------------
! integer :: im,iatom,ifreq,itypat,lpawu,tndim
! real(dp) :: hu_shift2
! character(len=500) :: message
! *********************************************************************

! do iatom = 1 , cryst_struc%natom
!   lpawu=self%hdc%matlu(iatom)%lpawu
!   tndim=2*lpawu+1
!   itypat=cryst_struc%typat(iatom)
!   if(lpawu/=-1) then
!     self%qmc_shift(iatom) = zero
!     do im =1, 2*tndim-1
!!       write(std_out,*)"make before",self%qmc_shift(iatom)
!       self%qmc_shift(iatom) = self%qmc_shift(iatom) + hu(itypat)%uqmc(im)
!!       write(std_out,*)"make after",self%qmc_shift(iatom)
!     enddo
!     self%qmc_shift(iatom) = self%qmc_shift(iatom) / two
!     hu_shift2 = hu(itypat)%uqmc(1)
!
!     do im = 2*tndim, 2*tndim + 2*tndim -3
!       hu_shift2 = hu_shift2 + hu(itypat)%uqmc(im)
!     enddo

!     hu_shift2 = hu_shift2 / two
!     write(message,'(2a,i4)')  ch10,'  -------> For Correlated atom',iatom
!     call wrtout(std_out,  message,'COLL')

!     if(abs(self%qmc_shift(iatom)-hu_shift2)>tol6) then
!       write(message,'(2a,2f16.7)')  "  Shift for QMC is not correctly"&
!&      ," computed",self%qmc_shift(iatom),hu_shift2
!       ABI_ERROR(message)
!     endif ! shifts not equals
!
!     write(message,'(4x,a,f16.7)')  &
!&     "  Shift for QMC (used to compute G(w)) is (in Ha) :",&
!&     self%qmc_shift(iatom)
!     call wrtout(std_out,  message,'COLL')

!     self%qmc_xmu(iatom)=-self%qmc_shift(iatom)
!     self%qmc_xmu(iatom)=zero
!     write(message,'(4x,a,f16.7)')  &
!&     "Artificial Shift used in QMC AND to compute G is (in Ha) :",self%qmc_xmu(iatom)
!     ABI_WARNING(message)

!   endif ! lpawu/=1
! enddo ! natom

! if(present(apply)) then
!   if (apply) then
!     write(message,'(5x,a,f16.7,a)')  " Shifts applied to self"
!     call wrtout(std_out,  message,'COLL')
!     do ifreq=1,self%nw
!       call shift_matlu(self%oper(ifreq)%matlu,cryst_struc%natom,cmplx(self%qmc_shift,0.d0,kind=dp),1)
!       call shift_matlu(self%oper(ifreq)%matlu,cryst_struc%natom,cmplx(self%qmc_xmu,0.d0,kind=dp),1)
!     enddo
!   endif
! endif


!end subroutine make_qmcshift_self
!!***

!!****f* m_self/kramerskronig_self
!! NAME
!! kramerskronig_self
!!
!! FUNCTION
!!  Compute the real part of the self-energy using Kramers-Kronig formula.
!!
!! INPUTS
!!  self <type(self_type)>= variables related to self-energy
!!  selflimit = 0th order moment of the self-energy for each atom
!!  selfhdc = double counting
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  filapp = name of the file
!!
!! OUTPUT
!!
!! SOURCE

subroutine kramerskronig_self(self,selflimit,selfhdc,paw_dmft,filapp)

!Arguments ------------------------------------
 type(self_type), intent(inout) :: self
 type(matlu_type), intent(in) :: selfhdc(self%hdc%natom),selflimit(self%hdc%natom)
 type(paw_dmft_type), intent(in) :: paw_dmft
 character(len=fnlen), intent(in) :: filapp
!Local variables-------------------------------
 integer :: iatom,ifreq,im,im1,ispinor,ispinor1,isppol,jfreq
 integer :: lpawu,myproc,natom,ndim,nspinor,nsppol,unt
 character(len=2) :: tag_im,tag_im1
 character(len=4) :: tag_at
 character(len=50) :: tag_is
 character(len=500) :: message
! *********************************************************************

 !delta = 0.0000000

 myproc  = paw_dmft%myproc
 natom   = self%hdc%natom
 nspinor = self%hdc%nspinor
 nsppol  = self%hdc%nsppol

 !ABI_MALLOC(selftemp_imag,(self%nw))

 write(message,'(2a,i4)') ch10,'  ------ High-frequency limit of the Self-Energy'
 call wrtout(std_out,message,'COLL')

 call print_matlu(selflimit(:),natom,3)

!print norms
 write(message,'(2a,i4)')  ch10,'  ------ Double counting'
 call wrtout(std_out,message,'COLL')

 call print_matlu(selfhdc(:),natom,3)

!  Compute limit of Real Part and put in double counting energy.
! call copy_matlu(selfhdc,self%hdc%matlu,natom)
     !!write(6,*) "selfhdc   kramerskronig",selfhdc(1)%mat(1,1,1,1,1)
     !write(6,*) "selfr%hdc kramerskronig",self%hdc%matlu(1)%mat(1,1,1,1,1)

 do ifreq=1,self%nw
   if (self%distrib%procf(ifreq) /= myproc) cycle
   do jfreq=1,self%nw-1
     if (jfreq == ifreq) cycle
     do iatom=1,natom
       lpawu = selfhdc(iatom)%lpawu
       if (lpawu == -1) cycle
       self%oper(ifreq)%matlu(iatom)%mat(:,:,:) = self%oper(ifreq)%matlu(iatom)%mat(:,:,:) - &
                     & cmplx(aimag(self%oper(jfreq)%matlu(iatom)%mat(:,:,:)) / (self%omega(ifreq)-self%omega(jfreq)) &
                     & * (self%omega(jfreq+1)-self%omega(jfreq)),zero,kind=dp)
     end do  ! iatom
   end do  ! jfreq
   do iatom=1,natom
     lpawu = selfhdc(iatom)%lpawu
     if (lpawu == -1) cycle
     self%oper(ifreq)%matlu(iatom)%mat(:,:,:) = cmplx(dble(self%oper(ifreq)%matlu(iatom)%mat(:,:,:))/pi,&
         & aimag(self%oper(ifreq)%matlu(iatom)%mat(:,:,:)),kind=dp) + selflimit(iatom)%mat(:,:,:)
     ! write(6,*) "TWO FACTOR IS PUT BECAUSE OF MAXENT CODE ??"
     !self%oper(ifreq)%matlu(iatom)%mat(:,:,:) = self%oper(ifreq)%matlu(iatom)%mat(:,:,:) * half
      !                       & aimag(self%oper(ifreq)%matlu(iatom)%mat(:,:,:)),kind=dp) * half
   end do ! iatom
 end do ! ifreq

 call gather_oper(self%oper(:),self%distrib,paw_dmft,opt_ksloc=2)

 if (myproc == 0) then

   unt = get_unit()
   open(unit=unt,file=trim(filapp)//"_DFTDMFT_Self_realaxis.dat",status='unknown',form='formatted')
   rewind(unt)
   write(unt,'(4a)') "# Self-energy on the real axis, with the real part computed using Kramers-Kronig relations.",ch10, &
                     "# Real Frequency (Ha.)         Real part                Imaginary part",ch10
   tag_is = ""
   do iatom=1,natom
     lpawu = selfhdc(iatom)%lpawu
     if (lpawu == -1) cycle
     write(tag_at,'(i4)') iatom
     ndim = 2*lpawu + 1
     do isppol=1,nsppol
       if (nsppol == 2) tag_is = trim(adjustl(" for spin " // merge("up  ","down",isppol==1))) // "and"
       do ispinor=1,nspinor
         do ispinor1=1,nspinor
           do im=1,ndim
             write(tag_im,'(i2)') im + (ispinor-1)*ndim
             do im1=1,ndim
               write(tag_im1,'(i2)') im1 + (ispinor1-1)*ndim
               write(unt,'(9a)') "## Sigma_{",trim(adjustl(tag_im)),",",trim(adjustl(tag_im1)),"}",trim(adjustl(tag_is))," for atom ",trim(adjustl(tag_at)),ch10
               do ifreq=1,self%nw
                 write(unt,'(2x,393(es24.16e3,2x))') self%omega(ifreq),&
                   & dble(self%oper(ifreq)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,isppol)), &
                   & aimag(self%oper(ifreq)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,isppol))
               end do ! ifreq
               write(unt,*)
             end do ! im1
           end do ! im
         end do ! ispinor1
       end do ! ispinor
     end do ! isppol
   end do ! iatom

   close(unt)

 end if ! master node

 !do iatom=1,natom
 !  lpawu = self%oper(1)%matlu(iatom)%lpawu
 !  if (lpawu == -1) cycle
 !  ndim = 2*lpawu + 1
 !  do isppol=1,nsppol
 !      do ispinor=1,nspinor
 !        do ispinor1=1,nspinor
 !          do im=1,ndim
 !            do im1=1,ndim
               !write(6,*)
               !"realpart",real(selflimit(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
 !              do ifreq=1,self%nw
 !                selftemp_re(ifreq)=zero
 !                selftemp_imag(ifreq)=aimag(self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
  !               do jfreq=1,self%nw-1
  !                 if(jfreq==ifreq) cycle
!                   selftemp_re(ifreq)=selftemp_re(ifreq) -   &
! &
! aimag(self%oper(jfreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))  &
! &                   *(self%omega(ifreq)-self%omega(jfreq)) &
! &                   /((self%omega(ifreq)-self%omega(jfreq))**2+delta**2)&
! &                   *(self%omega(jfreq+1)-self%omega(jfreq))
   !                selftemp_re(ifreq)=selftemp_re(ifreq) -   &
 !&
 !aimag(self%oper(jfreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))  &
 !&                   /(self%omega(ifreq)-self%omega(jfreq)) *
 !(self%omega(jfreq+1)-self%omega(jfreq))
 !                enddo
 !                selftemp_re(ifreq)=selftemp_re(ifreq)/pi
                 !write(671,*)
                 !self%omega(ifreq),selftemp_re(ifreq),selftemp_imag(ifreq)
  !             enddo
!                 TEST*************************
!               do ifreq=1,self%nw
!                 selftemp_imag(ifreq)=zero
!                 do jfreq=1,self%nw-1
!                   if(jfreq==ifreq) cycle
!!                   selftemp_re(ifreq)=selftemp_re(ifreq) -   &
!! &
!aimag(self%oper(jfreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))  &
!! &                   *(self%omega(ifreq)-self%omega(jfreq)) &
!! &                   /((self%omega(ifreq)-self%omega(jfreq))**2+delta**2)&
!! &                   *(self%omega(jfreq+1)-self%omega(jfreq))
!                   selftemp_imag(ifreq)=selftemp_imag(ifreq) +    &
! &                   selftemp_re(jfreq)  &
! &                   /(self%omega(ifreq)-self%omega(jfreq)) *
! (self%omega(jfreq+1)-self%omega(jfreq))
!                 enddo
!                 selftemp_imag(ifreq)=selftemp_imag(ifreq)/pi
!                 write(672,*)
!                 self%omega(ifreq),selftemp_re(ifreq),selftemp_imag(ifreq)
!               enddo
!                 TEST*************************
              ! write(6,*) "TWO FACTOR IS PUT BECAUSE OF MAXENT CODE ??"
  !             do ifreq=1,self%nw
!                 write(68,*)
!                 self%omega(ifreq),selftemp_re(ifreq),selftemp_imag(ifreq)
 !                selftemp_re(ifreq)=selftemp_re(ifreq)+ &
 !&                 real(selflimit(iatom)%mat(im,im1,isppol,ispinor,ispinor1)- &
 !&                 selfhdc(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
 !                self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)&
 ! &
 ! =cmplx(selftemp_re(ifreq),selftemp_imag(ifreq),kind=dp)/two
  !&                       =cmplx(0.d0,selftemp_imag(ifreq),kind=dp)/two
!  &                       =cmplx(selftemp_re(ifreq),0.d0,kind=dp)/two
  !               self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)&
  !&                       =cmplx(selftemp_re(ifreq),0.d0,kind=dp)/two
  !&                       =cmplx(0.d0,0.d0,kind=dp)/two
!    The factor two is here to compensate for the factor two in OmegaMaxent..
!  &                       =cmplx(selftemp_re(ifreq),0.0,kind=dp)
  !               write(67,*)
  !               self%omega(ifreq),real(self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))&
  !               ,aimag(self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
  !             enddo
   !            write(67,*)
               !write(68,*)
               !!!!!!!!!! Z renormalization
!               i0=389
!               slope=(selftemp_re(i0+1)-selftemp_re(i0))/&
!                     (self%omega(i0+1)-self%omega(i0))
!               y0= selftemp_re(i0)
!               do ifreq=1,self%nw
!                 selftemp_re(ifreq)=slope * (self%omega(ifreq)-self%omega(i0))
!                 + y0
!                 selftemp_imag(ifreq)=zero
!                 self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)&
!  &
!  =cmplx(selftemp_re(ifreq),selftemp_imag(ifreq),kind=dp)/two
!                 write(6777,*)
!                 self%omega(ifreq),real(self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)),aimag(self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
!               enddo
               !!!!!!!!!!
!             enddo
!           enddo
!         enddo
!       enddo
!     enddo ! isppol
! end do ! iatom

end subroutine kramerskronig_self
!!***

!!****f* m_self/selfreal2imag_self
!! NAME
!! selfreal2imag_self
!!
!! FUNCTION
!!  Print self on imaginary axis from real axis
!!
!! INPUTS
!!  selfr = self on real axis
!!  self  = self on imaginary axis
!!  filapp = name of the filename
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! OUTPUT
!!  self%qmc_shift in self <type(self_type)> = Self-energy
!!
!! SOURCE

subroutine selfreal2imag_self(selfr,self,filapp,paw_dmft)

!Arguments ------------------------------------
!type
 type(self_type), intent(in) :: selfr,self
 character(len=fnlen), intent(in) :: filapp
 type(paw_dmft_type), intent(in) :: paw_dmft
!Local variables-------------------------------
 integer :: iatom,ifreq,im,im1,ispinor,ispinor1,isppol,jfreq
 integer :: lpawu,myproc,natom,ndim,nspinor,nsppol,unt
 logical :: triqs
 complex(dpc) :: omega
 !real(dp) :: delta
 type(self_type) :: selftempmatsub
 character(len=2) :: tag_im,tag_im1
 character(len=4) :: tag_at
 character(len=50) :: tag_is
 type(matlu_type), allocatable :: matlu_tmp(:)
! *********************************************************************

 !delta=0.0000000
 call initialize_self(selftempmatsub,paw_dmft)

 myproc  = paw_dmft%myproc
 natom   = self%hdc%natom
 nspinor = self%hdc%nspinor
 nsppol  = self%hdc%nsppol
!  Compute limit of Real Part and put in double counting energy.
! call copy_matlu(selfhdc,self%hdc%matlu,natom)
     !write(6,*) "self3",aimag(selfr%oper(489)%matlu(1)%mat(1,1,1,1,1))

 triqs = (paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7)

 if (triqs) then
   ABI_MALLOC(matlu_tmp,(natom))
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),matlu_tmp(:))
 end if ! triqs

 do ifreq=1,self%nw
   if (self%distrib%procf(ifreq) /= myproc) cycle
   omega = cmplx(zero,self%omega(ifreq),kind=dp)
   do jfreq=1,selfr%nw-1
     do iatom=1,natom
       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       selftempmatsub%oper(ifreq)%matlu(iatom)%mat(:,:,:) = selftempmatsub%oper(ifreq)%matlu(iatom)%mat(:,:,:) - &
           & aimag(selfr%oper(jfreq)%matlu(iatom)%mat(:,:,:)) / &
           & (omega-selfr%omega(jfreq)) * (selfr%omega(jfreq+1)-selfr%omega(jfreq))
     end do ! iatom
   end do ! jfreq
   call fac_matlu(selftempmatsub%oper(ifreq)%matlu(:),natom,cone/pi)
   if (triqs) then ! add the missing high-frequency moment when using TRIQS
     call add_matlu(selftempmatsub%oper(ifreq)%matlu(:),self%moments(1)%matlu(:),matlu_tmp(:),natom,1)
     call copy_matlu(matlu_tmp(:),selftempmatsub%oper(ifreq)%matlu(:),natom)
   end if ! triqs
 end do ! ifreq

 if (triqs) then
   call destroy_matlu(matlu_tmp(:),natom)
   ABI_FREE(matlu_tmp)
 end if ! triqs

 call gather_oper(selftempmatsub%oper(:),self%distrib,paw_dmft,opt_ksloc=2)

 if (myproc == 0) then
   unt = get_unit()
   open(unit=unt,file=trim(filapp)//"_DFTDMFT_Self_backtransform.dat",status='unknown',form='formatted')
   write(unt,'(6a)') "# Hilbert transform of the analytically continued self-energy. To be compared with the actual",ch10, &
             & "# self-energy on the imaginary axis in the cubic basis.",ch10,"# Matsubara Frequency (Ha.)       Real part                Imaginary part",ch10
   tag_is = ""
   do iatom=1,natom
     lpawu = self%hdc%matlu(iatom)%lpawu
     if (lpawu == -1) cycle
     write(tag_at,'(i4)') iatom
     ndim = 2*lpawu + 1
     do isppol=1,nsppol
       if (nsppol == 2) tag_is = trim(adjustl(" for spin " // merge("up  ","down",isppol==1))) // "and"
       do ispinor=1,nspinor
         do ispinor1=1,nspinor
           do im=1,ndim
             write(tag_im,'(i2)') im + (ispinor-1)*ndim
             do im1=1,ndim
               write(tag_im1,'(i2)') im1 + (ispinor1-1)*ndim
               !do jfreq=1,selfr%nw-1
               !  write(6700,*)  selfr%omega(jfreq),aimag(selfr%oper(jfreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
               !enddo
               !  write(6700,*)
               write(unt,'(9a)') "## Sigma_{",trim(adjustl(tag_im)),",",trim(adjustl(tag_im1)),"}",trim(adjustl(tag_is))," for atom ",trim(adjustl(tag_at)),ch10
               do ifreq=1,self%nw
                 write(unt,'(2x,393(es24.16e3,2x))') self%omega(ifreq),&
                   & dble(selftempmatsub%oper(ifreq)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,isppol)),&
                   & aimag(selftempmatsub%oper(ifreq)%matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,isppol))
               end do ! ifreq
               write(unt,*)
             end do ! im1
           end do ! im
         end do ! ispinor1
       end do ! ispinor
     end do ! isppol
   end do ! iatom
   close(unt)
 end if ! master node

 call destroy_self(selftempmatsub)

end subroutine selfreal2imag_self
!!***

END MODULE m_self
!!***
