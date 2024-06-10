!!****m* ABINIT/m_xg_nonlop
!! NAME
!! m_xg_nonlop
!!
!! FUNCTION
!!  This module provides functions to compute the nonlocal operator by means of the BLAS GEMM
!!  routine. By treating ndat simultaneous wavefunctions, it is able to exploit BLAS3 routines,
!!  which leads to excellent CPU efficiency and OpenMP scalability.
!!
!! COPYRIGHT
!! Copyright (C) 2022-2022 ABINIT group (LB)
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

module m_xg_nonlop

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
! use m_abi_linalg
 use m_xg
 use m_xomp
#ifdef HAVE_OPENMP
 use omp_lib
#endif

 use m_time, only : timab
 use defs_abitypes, only : MPI_type

 use m_paw_ij, only : paw_ij_type
 use m_pawtab, only : pawtab_type

 implicit none

 private

 double precision, parameter :: inv_sqrt2 = 1/sqrt2

 integer, parameter :: tim_getcprj    = 1091
 integer, parameter :: tim_apply_prj  = 1092
 integer, parameter :: tim_apply_Aij  = 1093
 integer, parameter :: tim_mult_cprj  = 1094
 integer, parameter :: tim_getXSX     = 1095
 integer, parameter :: tim_getXHX     = 1096
 integer, parameter :: tim_make_k     = 1097
 integer, parameter :: tim_make_Dij   = 1098
 integer, parameter :: tim_make_Sij   = 1099
 integer, parameter :: tim_getHmeSX   = 1089
 integer, parameter :: tim_iter_refinement = 1088

 integer, parameter :: tim_getcprj_gemm     = 1030
 integer, parameter :: tim_getcprj_copy     = 1031
 integer, parameter :: tim_getcprj_mpi      = 1032

 integer, parameter :: tim_apply_prj_gemm   = 1035
 integer, parameter :: tim_apply_prj_copy   = 1036
 integer, parameter :: tim_apply_prj_mpi    = 1037

! integer, parameter :: tim_apply_Aij_gemm   = 1040
! integer, parameter :: tim_apply_Aij_copy   = 1041

 integer, parameter :: tim_mult_cprj_gemm   = 1045
 integer, parameter :: tim_mult_cprj_copy   = 1046
 integer, parameter :: tim_mult_cprj_mpi    = 1047

 type,public :: xg_nonlop_t

   integer :: cplex
   integer :: cprjdim
   integer :: comm_atom
   integer :: comm_band
   integer :: npw_k
   integer :: total_npw_k
   integer :: max_npw_k
   integer :: me_band
   integer :: my_natom
   integer :: natom
   integer :: mkmem
   integer :: nlmn_max
   integer :: ntypat
   integer :: nspinor
   integer :: space_pw
   integer :: space_Dij
   logical :: paw
   real(dp) :: weight

   integer, pointer :: mpi_atmtab(:)
   integer, pointer :: indlmn(:,:,:)
   integer, pointer :: nattyp(:)

   integer,allocatable :: l_npw_k(:)
   integer,allocatable :: l_shift_npw_k(:)

   real(dp), pointer :: sij_triangular_mat(:,:)

   integer, allocatable :: nlmn_natom(:)
   integer, allocatable :: nlmn_ntypat(:)

   type(xg_t),pointer :: projectors(:)
   type(xg_t),pointer :: projectors_k

   ! paw only:
   type(xg_t),pointer :: gram_proj(:)
   type(xg_t),pointer :: gram_proj_k

   type(xg_t) :: Dij_all
   type(xg_t) :: Sij_all
   type(xgBlock_t),allocatable :: Dij(:)
   type(xgBlock_t),allocatable :: Sij(:)

   type(xg_t) :: Sijm1_all
   type(xg_t) :: invSij_approx_all
   type(xgBlock_t),allocatable :: Sijm1(:)
   type(xgBlock_t),allocatable :: invSij_approx(:)
   ! end paw only

 end type xg_nonlop_t
!!***

  public :: xg_nonlop_getcprj
  public :: xg_nonlop_init
  public :: xg_nonlop_make_k
  public :: xg_nonlop_destroy
  public :: xg_nonlop_make_Dij    ! paw only
  public :: xg_nonlop_make_Sij    ! paw only
  public :: xg_nonlop_destroy_Dij ! paw only
  public :: xg_nonlop_destroy_Sij ! paw only
  public :: xg_nonlop_apply_Aij
  public :: xg_nonlop_precond_iterative_refinement
  public :: xg_nonlop_mult_cprj
  public :: xg_nonlop_apply_prj
  public :: xg_nonlop_colwiseXAX
  public :: xg_nonlop_getXAX
  public :: xg_nonlop_getXSX ! paw only
  public :: xg_nonlop_getXHX
  public :: xg_nonlop_getAX
  public :: xg_nonlop_getHX
  public :: xg_nonlop_getSX   ! paw only
  public :: xg_nonlop_getSm1X ! paw only
  public :: xg_nonlop_getHmeSX ! paw only

contains
!!***

!!****f* m_xg_nonlop/init_xg_nonlop
!! NAME
!! init_xg_nonlop
!!
!! FUNCTION
!! Initalization of the xg_nonlop_kpt array
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
 subroutine xg_nonlop_init(xg_nonlop,indlmn,mpi_atmtab,my_natom,nattyp,mkmem,ntypat,nspinor,ucvol,usepaw,&
     me_band,comm_band,comm_atom)

   integer ,intent(in) :: me_band,comm_band,comm_atom
   integer ,intent(in) :: my_natom
   integer ,intent(in) :: mkmem
   integer ,intent(in) :: ntypat
   integer ,intent(in) :: nspinor
   integer ,intent(in) :: usepaw
   real(dp),intent(in) :: ucvol
   type(xg_nonlop_t),intent(inout) :: xg_nonlop

   integer,intent(in),target :: mpi_atmtab(:)
   integer,intent(in),target :: indlmn(:,:,:)
   integer,intent(in),target :: nattyp(:)

   integer :: itypat,cprjdim,nlmn,nlmn_max,natom,shift

   xg_nonlop%mkmem=mkmem
   xg_nonlop%nspinor=nspinor
   xg_nonlop%comm_atom=comm_atom
   xg_nonlop%my_natom=my_natom
   xg_nonlop%me_band=me_band
   xg_nonlop%comm_band=comm_band

   xg_nonlop%paw=usepaw==1

   xg_nonlop%space_pw=0
   xg_nonlop%space_Dij=0

   xg_nonlop%weight=four_pi/sqrt(ucvol)
   xg_nonlop%ntypat=ntypat

   xg_nonlop%mpi_atmtab=>mpi_atmtab
   xg_nonlop%nattyp=>nattyp
   xg_nonlop%indlmn=>indlmn

   natom=0
   do itypat=1,ntypat
     natom = natom + nattyp(itypat)
   end do

   ABI_MALLOC(xg_nonlop%nlmn_ntypat,(ntypat))
   ABI_MALLOC(xg_nonlop%nlmn_natom,(natom))

   cprjdim=0
   shift=0
   nlmn_max=0
   do itypat=1,ntypat
     nlmn = count(indlmn(3,:,itypat)>0)
     if (nlmn>nlmn_max) nlmn_max=nlmn
     cprjdim = cprjdim + nattyp(itypat)*nlmn
     xg_nonlop%nlmn_ntypat(itypat) = nlmn
     xg_nonlop%nlmn_natom(1+shift:nattyp(itypat)+shift) = nlmn
     shift = shift + nattyp(itypat)
   end do
   xg_nonlop%nlmn_max = nlmn_max
   xg_nonlop%natom = natom
   xg_nonlop%cprjdim = cprjdim

   ABI_MALLOC(xg_nonlop%projectors,(mkmem))
   if (xg_nonlop%paw) ABI_MALLOC(xg_nonlop%gram_proj,(mkmem))

 end subroutine xg_nonlop_init
!!***

 subroutine xg_nonlop_destroy(xg_nonlop)

  type(xg_nonlop_t),intent(inout) :: xg_nonlop

  integer :: ikpt

  if (xg_nonlop%paw) then
    call xg_nonlop_destroy_Sij(xg_nonlop) ! Can be destroyed before
    call xg_nonlop_destroy_Dij(xg_nonlop) ! Can be destroyed before
  end if

  ABI_FREE(xg_nonlop%nlmn_ntypat)
  ABI_FREE(xg_nonlop%nlmn_natom)

  if (allocated(xg_nonlop%l_npw_k)) then
    ABI_FREE(xg_nonlop%l_npw_k)
  end if
  if (allocated(xg_nonlop%l_shift_npw_k)) then
    ABI_FREE(xg_nonlop%l_shift_npw_k)
  end if

  do ikpt=1,size(xg_nonlop%projectors)
    call xg_free(xg_nonlop%projectors(ikpt))
    if (xg_nonlop%paw) call xg_free(xg_nonlop%gram_proj(ikpt))
  end do
  ABI_FREE(xg_nonlop%projectors)

  if (xg_nonlop%paw) then
    ABI_FREE(xg_nonlop%gram_proj)
    if (allocated(xg_nonlop%invSij_approx)) then
      ABI_FREE(xg_nonlop%invSij_approx)
    end if
    call xg_free(xg_nonlop%invSij_approx_all)
  end if

 end subroutine xg_nonlop_destroy
!!***

 subroutine xg_nonlop_make_Dij(xg_nonlop,paw_ij,isppol,atindx)

  type(xg_nonlop_t),intent(inout) :: xg_nonlop
  type(paw_ij_type),intent(in)    :: paw_ij(:)
  integer,intent(in)              :: isppol
  integer,intent(in)              :: atindx(:)

  logical :: paral_atom
  integer :: iatom, iatom_tot, iatom_type, ierr, nlmn, nlmn_max, natom, nspinor, shift
  integer :: ilmn, jlmn, j0lmn, jjlmn, ijlmn
  integer :: cplex_alldij,cplex_dij,isp,isps,jsp,jsps,ijsp
  integer,allocatable :: l_cplex(:)
  real(dp),pointer :: Dij_(:,:),Dij_all_(:,:)
  real(dp) :: tsec(2)

! *************************************************************************

  call timab(tim_make_Dij,1,tsec)

  if (.not.xg_nonlop%paw) then
    ABI_ERROR('Not implemented with paw=False.')
  end if

  if (isppol/=1) then ! isppol must be 1 or 2 if nspinor==1, and must be 1 of nspinor==2
    if (isppol/=2.or.xg_nonlop%nspinor/=1) then
      ABI_ERROR('wrong isppol')
    end if
  end if

  nspinor  = xg_nonlop%nspinor
  natom    = xg_nonlop%natom
  nlmn_max = xg_nonlop%nlmn_max

  paral_atom=(xmpi_comm_size(xg_nonlop%comm_atom)>1)

  ABI_MALLOC(l_cplex,(natom))
  l_cplex=0
  do iatom=1,xg_nonlop%my_natom ! loop over atoms treated by this proc
    iatom_tot=iatom;if (paral_atom) iatom_tot=xg_nonlop%mpi_atmtab(iatom)
    l_cplex(iatom_tot)=paw_ij(iatom)%cplex_dij
  end do
  call xmpi_sum(l_cplex,xg_nonlop%comm_atom,ierr)
  cplex_alldij = 0
  do iatom=1,natom ! loop over all atoms
    if (cplex_alldij<l_cplex(iatom)) cplex_alldij=l_cplex(iatom)
  end do
  ABI_FREE(l_cplex)

  if (cplex_alldij==1) then
    xg_nonlop%space_Dij=SPACE_R
  else if (cplex_alldij==2) then
    xg_nonlop%space_Dij=SPACE_C
  else
    ABI_ERROR('Bad cplex_alldij')
  end if

  call xg_init(xg_nonlop%Dij_all,xg_nonlop%space_Dij,nspinor*nlmn_max,nspinor*nlmn_max*natom,xmpi_comm_null)
  ABI_MALLOC(xg_nonlop%Dij,(natom))

  do iatom=1,natom ! loop over all atoms, even if paral atoms
    nlmn=xg_nonlop%nlmn_natom(iatom)
    shift=1+(iatom-1)*nspinor*nlmn_max
    call xg_setBlock(xg_nonlop%Dij_all,xg_nonlop%Dij(iatom),nspinor*nlmn,nspinor*nlmn,fcol=shift)
  end do

  do iatom=1,xg_nonlop%my_natom ! loop over atoms treated by this proc
    iatom_tot=iatom;if (paral_atom) iatom_tot=xg_nonlop%mpi_atmtab(iatom)
    iatom_type = atindx(iatom_tot) ! convert iatom from input file to iatom ordered by type
    nlmn=xg_nonlop%nlmn_natom(iatom_type)
    cplex_dij=paw_ij(iatom)%cplex_dij
    call xgBlock_reverseMap(xg_nonlop%Dij(iatom_type),Dij_)
    do jsp=1,nspinor
      jsps = (jsp-1)*nlmn
      do jlmn=1,nlmn
        j0lmn=jlmn*(jlmn-1)/2
        jjlmn=j0lmn+jlmn
        do isp=1,nspinor
          isps = (isp-1)*nlmn
          if (nspinor==1) then
            ijsp = isppol
          else
            if (isp==jsp) then
              ijsp = isp
            else if (isp==1) then
              ijsp = 3 ! up/down
            else
              ijsp = 4 ! down/up
            end if
          end if
          ! see m_hamiltonian:pawdij2ekb
          if (cplex_dij==1) then
            Dij_(cplex_alldij*(jlmn+isps-1)+1,jlmn+jsps) = paw_ij(iatom)%dij(jjlmn,ijsp)
            if (cplex_alldij==2) Dij_(2*(jlmn+isps),jlmn+jsps) = zero
          else
            Dij_(2*(jlmn+isps-1)+1,jlmn+jsps) = paw_ij(iatom)%dij(2*jjlmn-1,ijsp)
            Dij_(2*(jlmn+isps)    ,jlmn+jsps) = paw_ij(iatom)%dij(2*jjlmn  ,ijsp)
          end if
          do ilmn=1,jlmn-1
            ! see m_hamiltonian:pawdij2ekb and opernlc_ylm
            ijlmn=j0lmn+ilmn
            if (cplex_dij==1) then
              Dij_(cplex_alldij*(ilmn+isps-1)+1,jlmn+jsps) = paw_ij(iatom)%dij(ijlmn,ijsp)
              if (cplex_alldij==2) Dij_(2*(ilmn+isps),jlmn+jsps) = zero
              Dij_(cplex_alldij*(jlmn+jsps-1)+1,ilmn+isps) = paw_ij(iatom)%dij(ijlmn,ijsp)
              if (cplex_alldij==2) Dij_(2*(jlmn+jsps),ilmn+isps) = zero
            else
              Dij_(2*(ilmn+isps-1)+1,jlmn+jsps) = paw_ij(iatom)%dij(2*ijlmn-1,ijsp)
              Dij_(2*(ilmn+isps)    ,jlmn+jsps) = paw_ij(iatom)%dij(2*ijlmn  ,ijsp)
              Dij_(2*(jlmn+jsps-1)+1,ilmn+isps) = paw_ij(iatom)%dij(2*ijlmn-1,ijsp)
              Dij_(2*(jlmn+jsps)    ,ilmn+isps) =-paw_ij(iatom)%dij(2*ijlmn  ,ijsp)
            end if
          end do
        end do
      end do
    end do
  end do

! Communication in case of distribution over atomic sites
  if (paral_atom) then
    call xgBlock_reverseMap(xg_nonlop%Dij_all%self,Dij_all_)
    call xmpi_sum(Dij_all_,xg_nonlop%comm_atom,ierr)
  end if

  call timab(tim_make_Dij,2,tsec)

 end subroutine xg_nonlop_make_Dij
!!***

 subroutine xg_nonlop_make_Sij(xg_nonlop,pawtab,inv_sij)

  type(xg_nonlop_t),intent(inout) :: xg_nonlop
  type(pawtab_type),intent(in) :: pawtab(:)
  logical,optional,intent(in) :: inv_sij

  logical :: inv_sij_
  integer :: itypat, nlmn, nlmn_max, ntypat, shift
  integer :: ilmn, jlmn, j0lmn, jjlmn, ijlmn
  real(dp),pointer :: Sij_(:,:)
  real(dp) :: tsec(2)
  type(xg_t) :: work

! *************************************************************************

  call timab(tim_make_Sij,1,tsec)

  if (.not.xg_nonlop%paw) then
    ABI_ERROR('Not implemented with paw=False.')
  end if

  ntypat   = xg_nonlop%ntypat
  nlmn_max = xg_nonlop%nlmn_max

  call xg_init(xg_nonlop%Sij_all,SPACE_R,nlmn_max,nlmn_max*ntypat,xmpi_comm_self)
  ABI_MALLOC(xg_nonlop%Sij,(ntypat))

  inv_sij_ = .false.
  if (present(inv_sij)) inv_sij_ = inv_sij

  if (inv_sij_) then
    call xg_init(xg_nonlop%Sijm1_all,SPACE_R,nlmn_max,nlmn_max*ntypat,xmpi_comm_self)
    ABI_MALLOC(xg_nonlop%Sijm1,(ntypat))
  end if

  do itypat=1,ntypat

    nlmn=xg_nonlop%nlmn_ntypat(itypat)
    shift=1+(itypat-1)*nlmn_max
    call xg_setBlock(xg_nonlop%Sij_all,xg_nonlop%Sij(itypat),nlmn,nlmn,fcol=shift)
    call xgBlock_reverseMap(xg_nonlop%Sij(itypat),Sij_)
    do jlmn=1,nlmn
      j0lmn=jlmn*(jlmn-1)/2
      jjlmn=j0lmn+jlmn
      Sij_(jlmn,jlmn) = pawtab(itypat)%sij(jjlmn)
      do ilmn=1,jlmn-1
        ijlmn=j0lmn+ilmn
        Sij_(ilmn,jlmn) = pawtab(itypat)%sij(ijlmn)
        Sij_(jlmn,ilmn) = pawtab(itypat)%sij(ijlmn)
      end do
    end do
    if (inv_sij_) then
      call xg_init(work,SPACE_R,nlmn,nlmn,xmpi_comm_self)
      call xg_setBlock(xg_nonlop%Sijm1_all,xg_nonlop%Sijm1(itypat),nlmn,nlmn,fcol=shift)
      call xgBlock_invert_sy(xg_nonlop%Sijm1(itypat),work%self,xg_input=xg_nonlop%Sij(itypat))
      call xg_free(work)
    end if

  end do

  call timab(tim_make_Sij,2,tsec)

 end subroutine xg_nonlop_make_Sij
!!***

 subroutine xg_nonlop_destroy_Dij(xg_nonlop)

  type(xg_nonlop_t),intent(inout) :: xg_nonlop

! *************************************************************************

  if (.not.xg_nonlop%paw) then
    ABI_ERROR('Not implemented with paw=False.')
  end if

  if (allocated(xg_nonlop%Dij)) then
    ABI_FREE(xg_nonlop%Dij)
  end if
  call xg_free(xg_nonlop%Dij_all)

 end subroutine xg_nonlop_destroy_Dij
!!***

 subroutine xg_nonlop_destroy_Sij(xg_nonlop)

  type(xg_nonlop_t),intent(inout) :: xg_nonlop

! *************************************************************************

  if (.not.xg_nonlop%paw) then
    ABI_ERROR('Not implemented with paw=False.')
  end if

  if (allocated(xg_nonlop%Sij)) then
    ABI_FREE(xg_nonlop%Sij)
  end if
  call xg_free(xg_nonlop%Sij_all)
  if (allocated(xg_nonlop%Sijm1)) then
    ABI_FREE(xg_nonlop%Sijm1)
  end if
  call xg_free(xg_nonlop%Sijm1_all)

 end subroutine xg_nonlop_destroy_Sij
!!***

 subroutine xg_nonlop_make_k(xg_nonlop,ikpt,istwf_k,me_g0,npw_k,ffnl_k,ph3d_k,compute_proj,compute_invS_approx,compute_gram)

  type(xg_nonlop_t),intent(inout) :: xg_nonlop

  logical ,intent(in) :: compute_proj
  integer ,intent(in) :: ikpt
  integer ,intent(in) :: istwf_k
  integer, intent(in) :: me_g0
  integer, intent(in) :: npw_k
  real(dp), intent(in) :: ffnl_k(:,:,:,:)
  real(dp), intent(in) :: ph3d_k(:,:,:)
  logical ,optional,intent(in) :: compute_invS_approx
  logical ,optional,intent(in) :: compute_gram

  logical :: compute_gram_,compute_invS_approx_
  real(dp),pointer :: projectors_k_(:,:),gram_proj_k_(:,:),Sijm1_(:,:)
  real(dp) :: cphase(2)
  complex(dp) :: cil(4),ctmp
  integer :: ierr,iatom, iblock, shift, shiftc, shift_itypat, itypat, ilmn, jlmn, nlmn, ia, icol, il, ipw
  integer :: cplex,cols,rows,nlmn_max,ntypat,nmpi,me_g0_loc,space_work
  real(dp) :: tsec(2)
  type(xg_t) :: work
  type(xgBlock_t) :: projs

! *************************************************************************

  call timab(tim_make_k,1,tsec)

  me_g0_loc = -1
  if (istwf_k==1) then
    xg_nonlop%cplex=2
    xg_nonlop%space_pw=SPACE_C
  else ! istwf_k>1
    xg_nonlop%cplex=1
    xg_nonlop%space_pw=SPACE_CR
    me_g0_loc = 0
    if (istwf_k==2.and.me_g0==1) me_g0_loc = me_g0
  end if

  xg_nonlop%npw_k = npw_k
  if (allocated(xg_nonlop%l_npw_k)) then
    ABI_FREE(xg_nonlop%l_npw_k)
  end if
  nmpi = xmpi_comm_size(xg_nonlop%comm_band)
  ABI_MALLOC(xg_nonlop%l_npw_k,(nmpi))
  xg_nonlop%l_npw_k(:) = 0
  xg_nonlop%l_npw_k(xg_nonlop%me_band+1) = npw_k
  call xmpi_sum(xg_nonlop%l_npw_k,xg_nonlop%comm_band,ierr)
  xg_nonlop%total_npw_k = sum(xg_nonlop%l_npw_k)
  xg_nonlop%max_npw_k = maxval(xg_nonlop%l_npw_k)
  if (allocated(xg_nonlop%l_shift_npw_k)) then
    ABI_FREE(xg_nonlop%l_shift_npw_k)
  end if
  ABI_MALLOC(xg_nonlop%l_shift_npw_k,(nmpi))
  xg_nonlop%l_shift_npw_k(1) = 0
  do iblock=2,nmpi
    xg_nonlop%l_shift_npw_k(iblock) = xg_nonlop%l_shift_npw_k(iblock-1) + xg_nonlop%l_npw_k(iblock-1)
  end do

  xg_nonlop%projectors_k => xg_nonlop%projectors(ikpt)
  if (xg_nonlop%paw) xg_nonlop%gram_proj_k => xg_nonlop%gram_proj(ikpt)

  if (compute_proj) then

    ntypat = xg_nonlop%ntypat

!   4pi/sqrt(ucvol) * (-i)^l
    cil(1) = ( 1.0_DP, 0.0_DP) * xg_nonlop%weight
    cil(2) = ( 0.0_DP,-1.0_DP) * xg_nonlop%weight
    cil(3) = (-1.0_DP, 0.0_DP) * xg_nonlop%weight
    cil(4) = ( 0.0_DP, 1.0_DP) * xg_nonlop%weight

    rows = npw_k
    cols = xg_nonlop%cprjdim

    call xg_init(xg_nonlop%projectors_k,xg_nonlop%space_pw,rows,cols,xg_nonlop%comm_band,me_g0=me_g0_loc)
    call xgBlock_reverseMap(xg_nonlop%projectors_k%self,projectors_k_)

    !TODO use openMP
    shift_itypat=0
    icol=1
    do itypat = 1, ntypat
      nlmn = xg_nonlop%nlmn_ntypat(itypat)
      do ia = 1, xg_nonlop%nattyp(itypat)
        iatom = ia + shift_itypat
        !! projectors = 4pi/sqrt(ucvol)* conj(ph3d) * ffnl * (-i)^l
        do ilmn=1,nlmn
          il=mod(xg_nonlop%indlmn(1,ilmn,itypat),4)+1
          do ipw=1, npw_k
            ctmp = cil(il) * cmplx( ph3d_k(1,ipw,iatom), -ph3d_k(2,ipw,iatom), kind=DP )
            cphase(1) =  dble(ctmp)
            cphase(2) = dimag(ctmp)
            projectors_k_(2*(ipw-1)+1:2*ipw,icol) = cphase(:) * ffnl_k(ipw, 1, ilmn, itypat)
          end do
          icol = icol + 1
        end do
      end do
      shift_itypat = shift_itypat + xg_nonlop%nattyp(itypat)
    end do

    compute_invS_approx_=.false.
    if (present(compute_invS_approx)) compute_invS_approx_ = compute_invS_approx
    if (compute_invS_approx_) then

      if (.not.xg_nonlop%paw) then
        ABI_ERROR('Not implemented with paw=False.')
      end if
      if (.not.allocated(xg_nonlop%invSij_approx)) then
        nlmn_max=xg_nonlop%nlmn_max
        if (xg_nonlop%space_pw==SPACE_CR) then
          call xg_init(xg_nonlop%invSij_approx_all,SPACE_R,nlmn_max,nlmn_max*ntypat,xmpi_comm_self)
        else
          call xg_init(xg_nonlop%invSij_approx_all,SPACE_C,nlmn_max,nlmn_max*ntypat,xmpi_comm_self)
        end if
        ABI_MALLOC(xg_nonlop%invSij_approx,(ntypat))

        shift_itypat=1
        do itypat = 1, ntypat
          nlmn = xg_nonlop%nlmn_ntypat(itypat)
          call xgBlock_setBlock(xg_nonlop%projectors_k%self,projs,rows,nlmn,fcol=shift_itypat)
          shift=1+(itypat-1)*nlmn_max
          call xg_setBlock(xg_nonlop%invSij_approx_all,xg_nonlop%invSij_approx(itypat),nlmn,nlmn,fcol=shift)
          if (xg_nonlop%space_pw==SPACE_CR) then
            call xgBlock_copy(xg_nonlop%Sijm1(itypat),xg_nonlop%invSij_approx(itypat))
          else
            call xgBlock_r2c(xg_nonlop%Sijm1(itypat),xg_nonlop%invSij_approx(itypat),1)
          end if
          space_work = xg_nonlop%space_pw
          if (xg_nonlop%space_pw==SPACE_CR) then
            space_work = SPACE_R
          end if
          call xg_init(work,space_work,nlmn,nlmn,xmpi_comm_self)
          call xgBlock_gemm('t','n',1.0d0,projs,projs,0.0d0,work%self,comm=xg_nonlop%comm_band)
          call xgBlock_add(xg_nonlop%invSij_approx(itypat),work%self)
          call xgBlock_invert_sy(xg_nonlop%invSij_approx(itypat),work%self)
          call xg_free(work)
          shift_itypat = shift_itypat + nlmn*xg_nonlop%nattyp(itypat)
        end do
      end if

    end if

    compute_gram_=.false.
    if (present(compute_gram)) compute_gram_ = compute_gram
    if (compute_gram_) then

      if (.not.xg_nonlop%paw) then
        ABI_ERROR('Not implemented with paw=False.')
      end if
      nlmn_max=xg_nonlop%nlmn_max
      if (xg_nonlop%space_pw==SPACE_CR) then
        cplex=1
        space_work=SPACE_R
      else
        cplex=2
        space_work=SPACE_C
      end if
      call xg_init(xg_nonlop%gram_proj_k,space_work,cols,cols,xmpi_comm_self)
      projs = xg_nonlop%projectors_k%self
      call xgBlock_gemm('t','n',1.0d0,projs,projs,0.0d0,xg_nonlop%gram_proj_k%self,comm=xg_nonlop%comm_band)
      call xgBlock_reverseMap(xg_nonlop%gram_proj_k%self,gram_proj_k_)
      shift=0
      shiftc=0
      do itypat = 1, ntypat
        nlmn = xg_nonlop%nlmn_ntypat(itypat)
        call xgBlock_reverseMap(xg_nonlop%Sijm1(itypat),Sijm1_)
        do ia = 1, xg_nonlop%nattyp(itypat)
          do jlmn=1,nlmn
            do ilmn=1,nlmn
              gram_proj_k_(shiftc+cplex*(ilmn-1)+1,shift+jlmn) = gram_proj_k_(shiftc+cplex*(ilmn-1)+1,shift+jlmn) &
                & + Sijm1_(ilmn,jlmn)
            end do
          end do
          shift  = shift  + nlmn
          shiftc = shiftc + cplex*nlmn
        end do
      end do

    end if

  end if ! compute_proj

  call timab(tim_make_k,2,tsec)

 end subroutine xg_nonlop_make_k
!!***

subroutine xg_nonlop_set_nmpi(xg_nonlop,X,cprjX,nmpi,blocksize,fft_representation)

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: X
   type(xgBlock_t), intent(in) :: cprjX
   integer, intent(out) :: blocksize,nmpi
   logical, intent(out) :: fft_representation

   integer :: nspinor,npw
   type(xgBlock_t) :: projs

   if (.not.associated(xg_nonlop%projectors_k)) then
     ABI_ERROR('projectors_k should be associated')
   end if

   projs = xg_nonlop%projectors_k%self

   npw = xg_nonlop%npw_k
   nspinor = xg_nonlop%nspinor
   blocksize = cols(cprjX)

   ! Check projector sizes
   if (npw/=rows(projs)) then
     ABI_ERROR('npw/=rows(projs)')
   end if
   if (xg_nonlop%cprjdim/=cols(projs)) then
     ABI_ERROR('cols(projs)/=cprjdim')
   end if
   ! Check cprj sizes
   if (rows(cprjX)/=xg_nonlop%cprjdim) then
     ABI_ERROR('rows(cprjX)/=cprjdim')
   end if
   ! rows(projs), cols(projs) and rows(cprj) are checked.
   ! Now we check : rows(X),cols(X) and cols(cprj) depending on mpi
   nmpi = xmpi_comm_size(comm(cprjX))
   if (nmpi==1) then ! sequential
     if (xmpi_comm_size(comm(X))/=1) then
       ABI_ERROR('size(comm(X))/=1')
     end if
     if (rows(X)/=npw*nspinor) then
       ABI_ERROR('rows(X)/=npw*nspinor')
     end if
     if (cols(X)*nspinor/=blocksize) then
       ABI_ERROR('cols(cprjX)/=cols(X)*nspinor')
     end if
   else ! MPI
     ! FFT representation (X have all rows, cols are distributed)
     if (rows(X)==xg_nonlop%total_npw_k*nspinor) then
       fft_representation = .True.
       if (xmpi_comm_size(comm(X))/=1) then
         ABI_ERROR('size(comm(X))/=1')
       end if
       if (cols(X)*nspinor/=blocksize) then
         ABI_ERROR('cols(cprjX)/=cols(X)*nspinor')
       end if
     ! Linalg representation (X have all cols, rows are distributed)
     else
       fft_representation = .False.
       if (comm(X)/=comm(cprjX)) then
         ABI_ERROR('comm(X)/=comm(cprjX)')
       end if
       if (rows(X)/=npw*nspinor) then
         ABI_ERROR('rows(X)/=npw*nspinor')
       end if
       if (cols(X)*nspinor/=blocksize*nmpi) then
         ABI_ERROR('cols(cprjX)*nmpi/=cols(X)*nspinor')
       end if
     end if
   end if

end subroutine xg_nonlop_set_nmpi

!!****f* m_xg_nonlop/xg_nonlop_getcprj
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
subroutine xg_nonlop_getcprj(xg_nonlop,X,cprjX,work_mpi)

   type(xg_nonlop_t), intent(in)    :: xg_nonlop
   type(xgBlock_t)  , intent(in   ) :: X
   type(xgBlock_t)  , intent(inout) :: cprjX,work_mpi

   real(dp) :: tsec(2)
   integer :: ierr,iblock,nmpi,npw,npw_max,blocksize,nspinor,shift,shift_row
   integer :: source,dest,tag,request,me_band,me_g0_loc
   logical :: fft_representation
   type(xgBlock_t) :: X_block,X_spinor,projs,work_mpi_npw
   type(xg_t) :: work_npw

   call timab(tim_getcprj,1,tsec)

   call xg_nonlop_set_nmpi(xg_nonlop,X,cprjX,nmpi,blocksize,fft_representation)

   nspinor = xg_nonlop%nspinor
   projs = xg_nonlop%projectors_k%self

   call xgBlock_reshape_spinor(X,X_spinor,nspinor,ROWS2COLS)

   if (nmpi==1) then

     call timab(tim_getcprj_gemm,1,tsec)
     call xgBlock_gemm('t','n',1.0d0,projs,X_spinor,0.d0,cprjX)
     call timab(tim_getcprj_gemm,2,tsec)

   else

     me_band=xg_nonlop%me_band

     if (fft_representation) then ! FFT representation (X have all rows, cols are distributed)

       npw_max = xg_nonlop%max_npw_k

       if (rows(work_mpi)/=npw_max) then
         ABI_ERROR('rows(work)/=npw_max')
       end if
       if (cols(work_mpi)/=xg_nonlop%cprjdim) then
         ABI_ERROR('rows(work)/=cprjdim')
       end if

       ! npw_max=npw_tot/Nmpi and blocksize=nband/Nmpi so size(work_npw) ~ 1/Nmpi^2
       call xg_init(work_npw,xg_nonlop%space_pw,npw_max,blocksize,xmpi_comm_null,me_g0=me_g0(X))

       do iblock=1,nmpi

         if (iblock==1) then
           npw = xg_nonlop%l_npw_k(me_band+1)
           shift_row = xg_nonlop%l_shift_npw_k(me_band+1)
           call xg_setBlock(work_npw,X_block,npw,blocksize)

           call timab(tim_getcprj_copy,1,tsec)
           call xgBlock_partialcopy(X_spinor,X_block,shift_row,0,BIG2SMALL)
           call timab(tim_getcprj_copy,2,tsec)

           call timab(tim_getcprj_gemm,1,tsec)
           call xgBlock_gemm('t','n',1.0d0,projs,X_block,0.d0,cprjX)
           call timab(tim_getcprj_gemm,2,tsec)

         else
           tag  = iblock
           dest = mod(me_band-(iblock-1),nmpi)
           if (dest<0) dest=dest+nmpi

           call timab(tim_getcprj_mpi,1,tsec)
           call xgBlock_mpi_isend(projs,dest,tag,request,comm=xg_nonlop%comm_band)
           call timab(tim_getcprj_mpi,2,tsec)

           source = mod(me_band+(iblock-1),nmpi)
           npw = xg_nonlop%l_npw_k(source+1)
           shift_row = xg_nonlop%l_shift_npw_k(source+1)
           me_g0_loc = -1
           if (xg_nonlop%space_pw==SPACE_CR) then
             me_g0_loc = me_g0(X)
             if (shift_row>0) me_g0_loc = 0
           end if
           call xg_setBlock(work_npw,X_block,npw,blocksize)

           call timab(tim_getcprj_copy,1,tsec)
           call xgBlock_partialcopy(X_spinor,X_block,shift_row,0,BIG2SMALL)
           call timab(tim_getcprj_copy,2,tsec)

           call xgBlock_setBlock(work_mpi,work_mpi_npw,npw_max,xg_nonlop%cprjdim)
           call xgBlock_free_reshape(work_mpi_npw,npw,xg_nonlop%cprjdim,new_me_g0=me_g0_loc)

           call timab(tim_getcprj_mpi,1,tsec)
           call xgBlock_mpi_recv(work_mpi_npw,source,tag,comm=xg_nonlop%comm_band)
           call timab(tim_getcprj_mpi,2,tsec)

           call timab(tim_getcprj_gemm,1,tsec)
           call xgBlock_gemm('t','n',1.0d0,work_mpi_npw,X_block,1.d0,cprjX)
           call timab(tim_getcprj_gemm,2,tsec)

           call timab(tim_getcprj_mpi,1,tsec)
           call xmpi_wait(request,ierr)
           call timab(tim_getcprj_mpi,2,tsec)

         end if

       end do

       call xg_free(work_npw)

     else ! Linalg representation (X have all cols, rows are distributed)

       call xgBlock_check(cprjX,work_mpi)

       do iblock=1,nmpi
         shift=1+(iblock-1)*blocksize
         call xgBlock_setBlock(X_spinor,X_block,rows(X_spinor),blocksize,fcol=shift)

         call timab(tim_getcprj_gemm,1,tsec)
         call xgBlock_gemm('t','n',1.0d0,projs,X_block,0.d0,work_mpi)
         call timab(tim_getcprj_gemm,2,tsec)
         ! We do the mpi sum outside xgBlock_gemm just to include the timing in tim_getcprj_mpi,
         ! (instead of tim_gemm_mpi).
         call timab(tim_getcprj_mpi,1,tsec)
         call xgBlock_mpi_sum(work_mpi,comm=xg_nonlop%comm_band)
         call timab(tim_getcprj_mpi,2,tsec)

         if (me_band==iblock-1) then
           call timab(tim_getcprj_copy,1,tsec)
           call xgBlock_copy(work_mpi,cprjX)
           call timab(tim_getcprj_copy,2,tsec)
         end if
       end do

     end if

   end if

   call timab(tim_getcprj,2,tsec)

 end subroutine xg_nonlop_getcprj
!!***

 subroutine xg_nonlop_apply_prj(xg_nonlop,cprjX,X,work_mpi)

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: cprjX
   type(xgBlock_t), intent(inout) :: X,work_mpi

   logical :: fft_representation
   integer :: nmpi,npw,blocksize,iblock,shift_col,nspinor!,nband,npw
   integer :: me_band,shift_row,npw_max
   integer :: ierr,source,dest,tag,request
   real(dp) :: tsec(2)
   type(xgBlock_t) :: projs
   type(xgBlock_t) :: X_spinor,X_block,work_mpi_npw
   type(xg_t) :: work_npw

   call timab(tim_apply_prj,1,tsec)

   call xg_nonlop_set_nmpi(xg_nonlop,X,cprjX,nmpi,blocksize,fft_representation)

   projs = xg_nonlop%projectors_k%self

   nspinor = xg_nonlop%nspinor
   call xgBlock_reshape_spinor(X,X_spinor,nspinor,ROWS2COLS)

   if (nmpi==1) then

     call timab(tim_apply_prj_gemm,1,tsec)
     call xgBlock_gemm('n','n',1.0d0,projs,cprjX,1.d0,X_spinor)
     call timab(tim_apply_prj_gemm,2,tsec)

   else

     me_band = xg_nonlop%me_band

     if (fft_representation) then ! FFT representation (X have all rows, cols are distributed)

       npw_max = xg_nonlop%max_npw_k

       if (rows(work_mpi)/=npw_max) then
         ABI_ERROR('rows(work)/=npw_max')
       end if
       if (cols(work_mpi)/=xg_nonlop%cprjdim) then
         ABI_ERROR('rows(work)/=cprjdim')
       end if

       ! npw_max=npw_tot/Nmpi and blocksize=nband/Nmpi so size(work_npw) ~ 1/Nmpi^2
       call xg_init(work_npw,xg_nonlop%space_pw,npw_max,blocksize,xmpi_comm_null)

       do iblock=1,nmpi

         if (iblock==1) then
           npw = xg_nonlop%l_npw_k(me_band+1)
           call xg_setBlock(work_npw,X_block,npw,blocksize)
           shift_row = xg_nonlop%l_shift_npw_k(me_band+1)

           call timab(tim_apply_prj_copy,1,tsec)
           call xgBlock_partialcopy(X_spinor,X_block,shift_row,0,BIG2SMALL)
           call timab(tim_apply_prj_copy,2,tsec)

           call timab(tim_apply_prj_gemm,1,tsec)
           call xgBlock_gemm('n','n',1.0d0,projs,cprjX,1.d0,X_block)
           call timab(tim_apply_prj_gemm,2,tsec)

           call timab(tim_apply_prj_copy,1,tsec)
           call xgBlock_partialcopy(X_block,X_spinor,shift_row,0,SMALL2BIG)
           call timab(tim_apply_prj_copy,2,tsec)
         else
           tag  = iblock
           dest = mod(me_band-(iblock-1),nmpi)
           if (dest<0) dest=dest+nmpi

           call timab(tim_apply_prj_mpi,1,tsec)
           call xgBlock_mpi_isend(projs,dest,tag,request,comm=xg_nonlop%comm_band)
           call timab(tim_apply_prj_mpi,2,tsec)

           source = mod(me_band+(iblock-1),nmpi)
           npw = xg_nonlop%l_npw_k(source+1)
           call xg_setBlock(work_npw,X_block,npw,blocksize)
           shift_row = xg_nonlop%l_shift_npw_k(source+1)

           call timab(tim_apply_prj_copy,1,tsec)
           call xgBlock_partialcopy(X_spinor,X_block,shift_row,0,BIG2SMALL)
           call timab(tim_apply_prj_copy,2,tsec)

           call xgBlock_setBlock(work_mpi,work_mpi_npw,npw_max,xg_nonlop%cprjdim)
           call xgBlock_free_reshape(work_mpi_npw,npw,xg_nonlop%cprjdim)

           call timab(tim_apply_prj_mpi,1,tsec)
           call xgBlock_mpi_recv(work_mpi_npw,source,tag,comm=xg_nonlop%comm_band)
           call timab(tim_apply_prj_mpi,2,tsec)

           call timab(tim_apply_prj_gemm,1,tsec)
           call xgBlock_gemm('n','n',1.0d0,work_mpi_npw,cprjX,1.d0,X_block)
           call timab(tim_apply_prj_gemm,2,tsec)

           call timab(tim_apply_prj_copy,1,tsec)
           call xgBlock_partialcopy(X_block,X_spinor,shift_row,0,SMALL2BIG)
           call timab(tim_apply_prj_copy,2,tsec)

           call timab(tim_apply_prj_mpi,1,tsec)
           call xmpi_wait(request,ierr)
           call timab(tim_apply_prj_mpi,2,tsec)
         end if

       end do

       call xg_free(work_npw)

     else ! Linalg representation (X have all cols, rows are distributed)

       call xgBlock_check(cprjX,work_mpi)

       do iblock=1,nmpi

         shift_col = 1 + mod(me_band+iblock-1,nmpi) * cols(cprjX)

         call xgBlock_setBlock(X_spinor,X_block,rows(X_spinor),cols(cprjX),fcol=shift_col)

         if (iblock==1) then

           call timab(tim_apply_prj_gemm,1,tsec)
           call xgBlock_gemm('n','n',1.0d0,projs,cprjX,1.d0,X_block)
           call timab(tim_apply_prj_gemm,2,tsec)

         else

           call timab(tim_apply_prj_mpi,1,tsec)
           tag = iblock
           dest = mod(me_band-(iblock-1),nmpi)
           if (dest<0) dest=dest+nmpi
           call xgBlock_mpi_isend(cprjX,dest,tag,request)
           source = mod(me_band+(iblock-1),nmpi)
           call xgBlock_mpi_recv(work_mpi,source,tag)
           call timab(tim_apply_prj_mpi,2,tsec)

           call timab(tim_apply_prj_gemm,1,tsec)
           call xgBlock_gemm('n','n',1.0d0,projs,work_mpi,1.d0,X_block)
           call timab(tim_apply_prj_gemm,2,tsec)

           call timab(tim_apply_prj_mpi,1,tsec)
           call xmpi_wait(request,ierr)
           call timab(tim_apply_prj_mpi,2,tsec)

         end if

       end do

     end if

   end if

   call timab(tim_apply_prj,2,tsec)

 end subroutine xg_nonlop_apply_prj
!!***

 subroutine xg_nonlop_apply_Aij(xg_nonlop,Aij,cprjin,cprjout,A_with_spin)

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in), target :: Aij(:)
   type(xgBlock_t), intent(in) :: cprjin
   type(xgBlock_t), intent(inout) :: cprjout
   logical,optional,intent(in) :: A_with_spin

   integer :: ia, iband, cprjdim, shift_itypat, iatom, itypat, nattyp, nlmn, shift
   integer :: space_aij, space_cprj, cplex, nlmn_max
   integer :: nspinor, nrows, ncols, nrows_A, ncols_A
   integer :: nlmn_1atom,nlmn_max_1atom,ncols_1atom

   type(xg_t)        :: cprjin_nlmn_max,cprjout_nlmn_max
   type(xg_t),target :: Aij_complex
   type(xgBlock_t) :: cprjin_nlmn,cprjout_nlmn
   type(xgBlock_t),pointer :: Aij_
   real(dp),pointer :: cprjin_nlmn_(:,:),cprjout_nlmn_(:,:)
   real(dp),pointer :: cprjin_(:,:),cprjout_(:,:)
   real(dp) :: tsec(2)
   logical :: aij_r2c,A_with_spin_

   call timab(tim_apply_Aij,1,tsec)

   call xgBlock_getsize(cprjin,nrows,ncols)

   cprjdim = xg_nonlop%cprjdim
   if (nrows/=cprjdim) then
     ABI_ERROR('nrows/=cprjdim')
   end if
   call xgBlock_check(cprjin,cprjout)

   space_cprj = space(cprjin)
   if (space_cprj==SPACE_C) then
     cplex=2
   else
     cplex=1
   end if
   space_aij = space(Aij(1))
   aij_r2c = .false.
   if (space_cprj/=space_aij) then
     if (space_aij==SPACE_R.and.space_cprj==SPACE_C) then
       aij_r2c = .true.
     else
       ABI_ERROR('space_aij and space_cprj are not compatible')
     end if
   end if

   nspinor=xg_nonlop%nspinor
   nlmn_max=xg_nonlop%nlmn_max

   A_with_spin_=.true.
   if (present(A_with_spin)) then
     A_with_spin_=A_with_spin
   end if
   if (A_with_spin_) then
     nlmn_max_1atom = nlmn_max * nspinor
     ncols_1atom = ncols / nspinor
   else
     nlmn_max_1atom = nlmn_max
     ncols_1atom = ncols
   end if

   ! Create work spaces for cprj of ONE atom for ALL bands
   call xg_init(cprjin_nlmn_max ,space_cprj,nlmn_max_1atom,ncols_1atom)
   call xg_init(cprjout_nlmn_max,space_cprj,nlmn_max_1atom,ncols_1atom)

   call xgBlock_reverseMap(cprjin ,cprjin_ )
   call xgBlock_reverseMap(cprjout,cprjout_)

   shift = 0
   shift_itypat = 0
   do itypat=1,xg_nonlop%ntypat
     nlmn=xg_nonlop%nlmn_ntypat(itypat)
     if (A_with_spin_) then
       nlmn_1atom = nlmn*nspinor
     else
       nlmn_1atom = nlmn
     end if
     nattyp=xg_nonlop%nattyp(itypat)
     call xg_setBlock(cprjin_nlmn_max ,cprjin_nlmn ,nlmn_1atom,ncols_1atom)
     call xg_setBlock(cprjout_nlmn_max,cprjout_nlmn,nlmn_1atom,ncols_1atom)
     call xgBlock_reverseMap(cprjin_nlmn ,cprjin_nlmn_ )
     call xgBlock_reverseMap(cprjout_nlmn,cprjout_nlmn_)
     if (aij_r2c) call xg_init(Aij_complex,SPACE_C,nlmn_1atom,nlmn_1atom)
     do ia=1,nattyp
       if (size(Aij)==xg_nonlop%ntypat) then
         Aij_ => Aij(itypat)
       else if (size(Aij)==xg_nonlop%natom) then
         iatom = ia + shift_itypat
         Aij_ => Aij(iatom)
       else
         ABI_ERROR('wrong size of Aij!')
       end if
       call xgBlock_getsize(Aij_,nrows_A,ncols_A)
       if (.not.aij_r2c) then
         if (nrows_A/=nlmn_1atom) then
           ABI_ERROR('nrows_A/=nlmn_1atom')
         end if
         if (ncols_A/=nlmn_1atom) then
           ABI_ERROR('ncols_A/=nlmn_1atom')
         end if
       else
         if (nrows_A/=nlmn) then
           ABI_ERROR('nrows_A/=nlmn')
         end if
         if (ncols_A/=nlmn) then
           ABI_ERROR('ncols_A/=nlmn')
         end if
       end if
       ! if needed, transfer real matrix to a complex one
       if (aij_r2c) then
         call xgBlock_r2c(Aij_,Aij_complex%self,nspinor)
         Aij_ => Aij_complex%self
       end if
       ! Copy cprj of ONE atom for ALL bands from cprjin to cprin_nlmn
       !call timab(tim_apply_Aij_copy,1,tsec)
       if (A_with_spin_) then
         do iband=1,ncols_1atom
           cprjin_nlmn_(1:cplex*nlmn,iband) = cprjin_(1+shift:cplex*nlmn+shift,1+nspinor*(iband-1))
         end do
         if (nspinor==2) then
           do iband=1,ncols_1atom
             cprjin_nlmn_(1+cplex*nlmn:2*cplex*nlmn,iband) = cprjin_(1+shift:cplex*nlmn+shift,nspinor*iband)
           end do
         end if
       else
         do iband=1,ncols_1atom
           cprjin_nlmn_(1:cplex*nlmn,iband) = cprjin_(1+shift:cplex*nlmn+shift,iband)
         end do
       end if
       !call timab(tim_apply_Aij_copy,2,tsec)

       !call timab(tim_apply_Aij_gemm,1,tsec)
!       call xgBlock_gemm('n','n',1.0d0,Aij_,cprjin_nlmn,0.d0,cprjout_nlmn,timing=.false.)
       call xgBlock_gemm('n','n',1.0d0,Aij_,cprjin_nlmn,0.d0,cprjout_nlmn)
       !call timab(tim_apply_Aij_gemm,2,tsec)

       !call timab(tim_apply_Aij_copy,1,tsec)
       if (A_with_spin_) then
         do iband=1,ncols_1atom
           cprjout_(1+shift:cplex*nlmn+shift,1+nspinor*(iband-1)) = cprjout_(1+shift:cplex*nlmn+shift,1+nspinor*(iband-1)) &
           & + cprjout_nlmn_(1:cplex*nlmn,iband)
         end do
         if (nspinor==2) then
           do iband=1,ncols_1atom
             cprjout_(1+shift:cplex*nlmn+shift,nspinor*iband) = cprjout_(1+shift:cplex*nlmn+shift,nspinor*iband) &
             & + cprjout_nlmn_(1+cplex*nlmn:2*cplex*nlmn,iband)
           end do
         end if
       else
         do iband=1,ncols_1atom
           cprjout_(1+shift:cplex*nlmn+shift,iband) = cprjout_(1+shift:cplex*nlmn+shift,iband) &
           & + cprjout_nlmn_(1:cplex*nlmn,iband)
         end do
       end if
       !call timab(tim_apply_Aij_copy,2,tsec)
       shift=shift+cplex*nlmn
     end do

     if (aij_r2c) call xg_free(Aij_complex)
     shift_itypat = shift_itypat + nattyp

   end do

   call xg_free(cprjin_nlmn_max)
   call xg_free(cprjout_nlmn_max)

   call timab(tim_apply_Aij,2,tsec)

 end subroutine xg_nonlop_apply_Aij
!!***

 subroutine xg_nonlop_precond_iterative_refinement(xg_nonlop,A,precond,cprj_in,cprj_out,cprj_work)

   type(xg_nonlop_t), intent(in)  :: xg_nonlop
   type(xgBlock_t), intent(in   ) :: A,precond(:)
   type(xgBlock_t), intent(in   ) :: cprj_in
   type(xgBlock_t), intent(inout) :: cprj_out,cprj_work

   integer :: iter,cprjdim,ncols,additional_steps_to_take
   real(dp), parameter :: tolerance = 1e-14 ! maximum relative error. TODO: use tolwfr ?
   type(xg_t) :: err
   real(dp) :: norm,max_err,previous_max_err,convergence_rate,tsec(2)

   call timab(tim_iter_refinement,1,tsec)

   ! Note that precond is block-diagonal whereas A is not

   cprjdim = xg_nonlop%cprjdim
   ncols   = cols(cprj_in)
   if (cprjdim/=rows(cprj_in)) then
     ABI_ERROR('Wrong size for cprj_in')
   end if
   call xgBlock_check(cprj_in,cprj_out)
   call xgBlock_check(cprj_in,cprj_work)
   if (rows(A)/=cprjdim.or.cols(A)/=cprjdim) then
     ABI_ERROR('Wrong size for A')
   end if

   call xg_init(err,SPACE_R,ncols,1,xmpi_comm_self)

   call xgBlock_colwiseNorm2(cprj_in,err%self,max_val=norm)

   ! Y_0 = PX (with P block diagonal and "close" to A^-1)
   call xgBlock_zero(cprj_out)
   call xg_nonlop_apply_Aij(xg_nonlop,precond,cprj_in,cprj_out,A_with_spin=.false.)

   additional_steps_to_take = -1
   do iter=1,30
     ! compute AY_i
     call xgBlock_gemm('n','n',1.0d0,A,cprj_out,0.0d0,cprj_work)
     ! RES = AY_i - X
     call xgBlock_saxpy(cprj_work,-1.0d0,cprj_in)
     call xgBlock_colwiseNorm2(cprj_work,err%self,max_val=max_err)
     max_err = sqrt(max_err / norm)
     if(max_err < tolerance .or. additional_steps_to_take == 1) then
       exit
       ! We might stall and never get to the specified precision because of machine errors.
       ! If we got to 1e-10, extrapolate convergence rate and determine the number of additional
       ! steps to take to reach precision
     else if(max_err < 1e-10 .and. additional_steps_to_take == -1) then
       convergence_rate = -LOG(1e-10) / iter
       additional_steps_to_take = CEILING(-LOG(tolerance/1e-10)/convergence_rate) + 1
     else if(additional_steps_to_take > 0) then
       if(previous_max_err<max_err)exit
       additional_steps_to_take = additional_steps_to_take - 1
     end if
     previous_max_err=max_err
     ! RES = X - AY_i
     call xgBlock_scale(cprj_work,-1.0d0,1)
     ! Y_(i+1) = Y_i + P RES
     call xg_nonlop_apply_Aij(xg_nonlop,precond,cprj_work,cprj_out,A_with_spin=.false.)
   end do

   call xg_free(err)

   call timab(tim_iter_refinement,2,tsec)

 end subroutine xg_nonlop_precond_iterative_refinement
!!***

 subroutine xg_nonlop_mult_cprj(xg_nonlop,cprj_left,cprj_right,res,blocksize)

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   integer, intent(in) :: blocksize
   type(xgBlock_t), intent(in) :: cprj_left,cprj_right
   type(xgBlock_t), intent(inout) :: res

   integer :: space_res
   integer :: blocksize_spinor,iblock_mpi,nblocks_mpi,shift_row,shift_col,shift_col_mpi
   integer :: iblock_left,iblock_right,nblocks_left,nblocks_right
   integer :: res_nrows,res_ncols,cprjdim
   integer :: nrows_r,nrows_l,ncols_r,ncols_l,nspinor
   integer :: comm_cprj,tag,request,ierr,me_band,source,dest
   logical :: multiblock
   real(dp) :: tsec(2)
   type(xg_t) :: res_work_mpi,res_block,res_block_mpi,cprj_work_mpi
   type(xgBlock_t) :: cprj_left_spinor,cprj_right_spinor

   call timab(tim_mult_cprj,1,tsec)

   comm_cprj = comm(cprj_right)
   nblocks_mpi = xmpi_comm_size(comm_cprj)

   cprjdim = xg_nonlop%cprjdim
   nspinor = xg_nonlop%nspinor

   call xgBlock_reshape_spinor(cprj_right,cprj_right_spinor,nspinor,COLS2ROWS)
   call xgBlock_reshape_spinor(cprj_left ,cprj_left_spinor ,nspinor,COLS2ROWS)

   nrows_r = rows(cprj_right_spinor)
   nrows_l = rows(cprj_left_spinor)

   ncols_r = cols(cprj_right_spinor)
   ncols_l = cols(cprj_left_spinor)

   me_band = xg_nonlop%me_band

   if (nrows_r/=nspinor*cprjdim) then
     ABI_ERROR("rows(cprj_right)/=nspinor*cprjdim")
   end if
   if (nrows_l/=nspinor*cprjdim) then
     ABI_ERROR("rows(cprj_left)/=nspinor*cprjdim")
   end if
   if (rows(res)/=nblocks_mpi*ncols_l) then
     ABI_ERROR("rows(res)/=nblocks_mpi*cols(cprj_left)")
   end if
   if (cols(res)/=nblocks_mpi*ncols_r) then
     ABI_ERROR("cols(res)/=nblocks_mpi*cols(cprj_right)")
   end if

   if (nblocks_mpi==1) then

     call timab(tim_mult_cprj_gemm,1,tsec)
     call xgBlock_gemm('t','n',1.0d0,cprj_left_spinor,cprj_right_spinor,1.d0,res)
     call timab(tim_mult_cprj_gemm,2,tsec)

   else

     blocksize_spinor = blocksize / nspinor
     nblocks_right = ncols_r / blocksize_spinor
     nblocks_left  = ncols_l / blocksize_spinor

     multiblock = .false.
     if (nblocks_right>1.or.nblocks_left>1) then
       multiblock = .true.
     end if

     res_nrows        = rows(res)
     res_ncols        = cols(res)

     space_res = space(cprj_right)

     call xg_init(res_block_mpi,space_res,ncols_l,ncols_r,xmpi_comm_null)
     if (multiblock) then
       call xg_init(res_block,space_res,blocksize_spinor,blocksize_spinor,xmpi_comm_null)
     end if
     call xg_init(res_work_mpi,space_res,res_nrows,res_ncols,xmpi_comm_null)
     call xgBlock_zero(res_work_mpi%self)

     call xg_init(cprj_work_mpi,space_res,nspinor*cprjdim,ncols_r,comm_cprj)

     do iblock_mpi=1,nblocks_mpi

       if (iblock_mpi==1) then

         call timab(tim_mult_cprj_gemm,1,tsec)
         call xgBlock_gemm('t','n',1.0d0,cprj_left_spinor,cprj_right_spinor,0.d0,res_block_mpi%self)
         call timab(tim_mult_cprj_gemm,2,tsec)

       else

         call timab(tim_mult_cprj_mpi,1,tsec)
         tag = iblock_mpi
         dest = mod(me_band-(iblock_mpi-1),nblocks_mpi)
         if (dest<0) dest=dest+nblocks_mpi
         call xgBlock_mpi_isend(cprj_right_spinor,dest,tag,request)
         source = mod(me_band+(iblock_mpi-1),nblocks_mpi)
         call xgBlock_mpi_recv(cprj_work_mpi%self,source,tag)
         call timab(tim_mult_cprj_mpi,2,tsec)

         call timab(tim_mult_cprj_gemm,1,tsec)
         call xgBlock_gemm('t','n',1.0d0,cprj_left_spinor,cprj_work_mpi%self,0.d0,res_block_mpi%self)
         call timab(tim_mult_cprj_gemm,2,tsec)

       end if

       if (.not.multiblock) then
         shift_row = me_band*ncols_l
         shift_col = mod(me_band+iblock_mpi-1,res_ncols)*ncols_r
         if (shift_col>=res_ncols) shift_col = shift_col - res_ncols

         call timab(tim_mult_cprj_copy,1,tsec)
         call xgBlock_partialcopy(res_block_mpi%self,res_work_mpi%self,shift_row,shift_col,SMALL2BIG)
         call timab(tim_mult_cprj_copy,2,tsec)

       else
         do iblock_right=1,nblocks_right
           do iblock_left=1,nblocks_left
             shift_row = (iblock_left-1)*blocksize_spinor
             shift_col = (iblock_right-1)*blocksize_spinor

             call timab(tim_mult_cprj_copy,1,tsec)
             call xgBlock_partialcopy(res_block_mpi%self,res_block%self,shift_row,shift_col,BIG2SMALL)
             call timab(tim_mult_cprj_copy,2,tsec)

             shift_row = (me_band*blocksize_spinor) + (iblock_left-1)*res_nrows/nblocks_left
             shift_col_mpi = mod(me_band+iblock_mpi-1,res_ncols/nblocks_right)*blocksize_spinor
             if (shift_col_mpi>=res_ncols/nblocks_right) shift_col_mpi = shift_col_mpi - res_ncols/nblocks_right
             shift_col = shift_col_mpi + (iblock_right-1)*res_ncols/nblocks_right

             call timab(tim_mult_cprj_copy,1,tsec)
             call xgBlock_partialcopy(res_block%self,res_work_mpi%self,shift_row,shift_col,SMALL2BIG)
             call timab(tim_mult_cprj_copy,2,tsec)

           end do
         end do
       end if

       if (iblock_mpi>1) then
         call timab(tim_mult_cprj_mpi,1,tsec)
         call xmpi_wait(request,ierr)
         call timab(tim_mult_cprj_mpi,2,tsec)
       end if

     end do

     call timab(tim_mult_cprj_mpi,1,tsec)
     call xgBlock_mpi_sum(res_work_mpi%self,comm=comm_cprj)
     call timab(tim_mult_cprj_mpi,2,tsec)

     call xgBlock_add(res,res_work_mpi%self)

     if (multiblock) then
       call xg_free(res_block)
     end if
     call xg_free(res_block_mpi)
     call xg_free(res_work_mpi)
     call xg_free(cprj_work_mpi)

   end if

   call timab(tim_mult_cprj,2,tsec)

 end subroutine xg_nonlop_mult_cprj
!!***

subroutine xg_nonlop_colwiseXAX(xg_nonlop,Aij,cprj,cprj_work,res)

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: cprj,Aij(:)
   type(xgBlock_t), intent(inout) :: cprj_work,res

   integer :: ncols
   type(xgBlock_t) :: cprj_spinor,cprj_work_spinor

   call xgBlock_check(cprj,cprj_work)
   ncols = cols(cprj)
   if (ncols/=xg_nonlop%nspinor*rows(res)) then
     ABI_ERROR('Wrong cols for cprj or res.')
   end if

   call xgBlock_zero(cprj_work)
   call xg_nonlop_apply_Aij(xg_nonlop,Aij,cprj,cprj_work)

   call xgBlock_reshape_spinor(cprj     ,cprj_spinor     ,xg_nonlop%nspinor,COLS2ROWS)
   call xgBlock_reshape_spinor(cprj_work,cprj_work_spinor,xg_nonlop%nspinor,COLS2ROWS)
   call xgBlock_colwiseDotProduct(cprj_spinor,cprj_work_spinor,res)

end subroutine xg_nonlop_colwiseXAX
!!***

subroutine xg_nonlop_getXAX(xg_nonlop,Aij,cprj_left,cprj_right,cprj_work,res,blocksize)

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: cprj_left,cprj_right,Aij(:)
   type(xgBlock_t), intent(inout) :: cprj_work,res
   integer,intent(in) :: blocksize

   call xgBlock_zero(cprj_work)
   call xg_nonlop_apply_Aij(xg_nonlop,Aij,cprj_right,cprj_work)

   call xg_nonlop_mult_cprj(xg_nonlop,cprj_left,cprj_work,res,blocksize)

 end subroutine xg_nonlop_getXAX
!!***

subroutine xg_nonlop_getXSX(xg_nonlop,cprj_left,cprj_right,cprj_work,res,blocksize)

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(inout) :: cprj_left,cprj_right,cprj_work,res
   integer,intent(in) :: blocksize

   real(dp) :: tsec(2)

   call timab(tim_getXSX,1,tsec)

   if (.not.xg_nonlop%paw) then
     ABI_ERROR('Not implemented with paw=False.')
   end if

   call xg_nonlop_getXAX(xg_nonlop,xg_nonlop%Sij,cprj_left,cprj_right,cprj_work,res,blocksize)

   call timab(tim_getXSX,2,tsec)

 end subroutine xg_nonlop_getXSX
!!***

subroutine xg_nonlop_getXHX(xg_nonlop,cprj_left,cprj_right,cprj_work,res,blocksize)

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(inout) :: cprj_left,cprj_right,cprj_work,res
   integer,intent(in) :: blocksize

   real(dp) :: tsec(2)

   call timab(tim_getXHX,1,tsec)

   call xg_nonlop_getXAX(xg_nonlop,xg_nonlop%Dij,cprj_left,cprj_right,cprj_work,res,blocksize)

   call timab(tim_getXHX,2,tsec)

 end subroutine xg_nonlop_getXHX
!!***

 subroutine xg_nonlop_getAX(xg_nonlop,Aij,Xin,cprjin,cprj_work,work_mpi,Xout)

   use iso_c_binding

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: cprjin,Aij(:)
   type(xgBlock_t), intent(inout) :: Xin,cprj_work,work_mpi
   type(xgBlock_t), optional, intent(inout) :: Xout

   integer :: nblocks,nspinor
   integer :: nrows,nrows_cprj
   integer :: ncols,ncols_cprj

   nblocks = xmpi_comm_size(comm(cprjin))

   nspinor = xg_nonlop%nspinor

   ! check sizes
   if (present(Xout)) then
     call xgBlock_check(Xin,Xout)
   end if
   call xgBlock_check(cprjin,cprj_work)
   call xgBlock_getsize(Xin,nrows,ncols)
   call xgBlock_getsize(cprjin,nrows_cprj,ncols_cprj)

   ! cprj_work = sum_j Saij cprjin
   call xgBlock_zero(cprj_work)
   call xg_nonlop_apply_Aij(xg_nonlop,Aij,cprjin,cprj_work)

   ! Xout = Xin + sum_ai pai cprj_work
   if (present(Xout)) then
     call xgBlock_copy(Xin,Xout)
     call xg_nonlop_apply_prj(xg_nonlop,cprj_work,Xout,work_mpi)
   else ! in place version
     call xg_nonlop_apply_prj(xg_nonlop,cprj_work,Xin,work_mpi)
   end if

 end subroutine xg_nonlop_getAX

 subroutine xg_nonlop_getHX(xg_nonlop,Xin,cprjin,cprj_work,work_mpi,Xout)

   use iso_c_binding

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: cprjin
   type(xgBlock_t), intent(inout) :: Xin,cprj_work,work_mpi
   type(xgBlock_t), optional, intent(inout) :: Xout

   call xg_nonlop_getAX(xg_nonlop,xg_nonlop%Dij,Xin,cprjin,cprj_work,work_mpi,Xout)

 end subroutine xg_nonlop_getHX

 subroutine xg_nonlop_getSX(xg_nonlop,Xin,cprjin,cprj_work,work_mpi,Xout)

   use iso_c_binding

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: cprjin
   type(xgBlock_t), intent(inout) :: Xin,cprj_work,work_mpi
   type(xgBlock_t), optional, intent(inout) :: Xout

   if (.not.xg_nonlop%paw) then
     ABI_ERROR('Not implemented with paw=False.')
   end if

   call xg_nonlop_getAX(xg_nonlop,xg_nonlop%Sij,Xin,cprjin,cprj_work,work_mpi,Xout)

 end subroutine xg_nonlop_getSX

 subroutine xg_nonlop_getSm1X(xg_nonlop,Xin,cprjin,cprj_work1,cprj_work2,work_mpi,Xout)

   use iso_c_binding

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: cprjin
   type(xgBlock_t), intent(inout) :: Xin,cprj_work1,cprj_work2,work_mpi
   type(xgBlock_t), optional, intent(inout) :: Xout

   integer :: nblocks,nspinor
   integer :: nrows,nrows_cprj
   integer :: ncols,ncols_cprj

   if (.not.xg_nonlop%paw) then
     ABI_ERROR('Not implemented with paw=False.')
   end if

   nblocks = xmpi_comm_size(comm(cprjin))

   nspinor = xg_nonlop%nspinor

   if (present(Xout)) then
     call xgBlock_check(Xin,Xout)
   end if
   call xgBlock_check(cprjin,cprj_work1)
   call xgBlock_check(cprjin,cprj_work2)
   call xgBlock_getsize(Xin,nrows,ncols)
   call xgBlock_getsize(cprjin,nrows_cprj,ncols_cprj)
   if (.not.allocated(xg_nonlop%invSij_approx)) then
     ABI_ERROR('invSij_approx not allocated')
   end if
   if (.not.associated(xg_nonlop%gram_proj_k)) then
     ABI_ERROR('gram_proj_k should be associated')
   end if

   call xg_nonlop_precond_iterative_refinement(xg_nonlop,xg_nonlop%gram_proj_k%self,xg_nonlop%invSij_approx,&
     & cprjin,cprj_work1,cprj_work2)
   call xgBlock_scale(cprj_work1,-1.0d0,1)

   ! Xout = Xin + sum_ai pai cprj_work
   if (present(Xout)) then
     call xgBlock_copy(Xin,Xout)
     call xg_nonlop_apply_prj(xg_nonlop,cprj_work1,Xout,work_mpi)
   else ! in place version
     call xg_nonlop_apply_prj(xg_nonlop,cprj_work1,Xin,work_mpi)
   end if

 end subroutine xg_nonlop_getSm1X

 subroutine xg_nonlop_getHmeSX(xg_nonlop,Xin,cprjin,Xout,eigen,cprj_work,work_mpi,no_H)

   use iso_c_binding

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: Xin,eigen,cprjin
   type(xgBlock_t), intent(inout) :: Xout,cprj_work,work_mpi
   logical,optional,intent(in) :: no_H

   real(dp) :: tsec(2)
   integer :: nblocks,shift
   integer :: nrows,ncols,nrows_out,ncols_out
   integer :: nrows_diag,ncols_diag
   integer :: nrows_cprj,ncols_cprj
   integer :: nrows_cprj_work,ncols_cprj_work
   integer :: nspinor
   logical :: no_H_

   call timab(tim_getHmeSX,1,tsec)

   if (.not.xg_nonlop%paw) then
     ABI_ERROR('Not implemented with paw=False.')
   end if

   nblocks = xmpi_comm_size(comm(cprjin))

   nspinor = xg_nonlop%nspinor

   call xgBlock_getsize(Xin,nrows,ncols)
   call xgBlock_getsize(Xout,nrows_out,ncols_out)
   call xgBlock_getsize(eigen,nrows_diag,ncols_diag)
   call xgBlock_getsize(cprjin,nrows_cprj,ncols_cprj)
   call xgBlock_getsize(cprj_work,nrows_cprj_work,ncols_cprj_work)
   if (ncols/=nrows_diag.or.ncols*nspinor/=nblocks*ncols_cprj.or.ncols_cprj/=ncols_cprj_work.or.ncols/=ncols_out) then
     ABI_ERROR('wrong ncols')
   end if
   if (nrows/=nrows_out) then
     ABI_ERROR('nrows/=nrows_out')
   end if
   if (nrows_cprj/=nrows_cprj_work) then
     ABI_ERROR('nrows_cprj/=nrows_cprj_work')
   end if
   if (ncols_diag/=1) then
     ABI_ERROR('ncols_diag should be one')
   end if

   ! cprj_work = sum_j Saij cprjin
   call xgBlock_zero(cprj_work)
   call xg_nonlop_apply_Aij(xg_nonlop,xg_nonlop%Sij,cprjin,cprj_work)

   ! cprj_work = - e cprj_work = -e sum_j Saij cprjin
   shift = xg_nonlop%me_band*ncols_cprj/xg_nonlop%nspinor
   call xgBlock_ymax(cprj_work,eigen,shift,nblocks,xg_nonlop%nspinor)

   no_H_=.False.
   if (present(no_H)) then
     no_H_ = no_H
   end if
   if (.not.no_H_) then
     ! cprj_work = sum_j Daij cprjin + cprj_work
     call xg_nonlop_apply_Aij(xg_nonlop,xg_nonlop%Dij,cprjin,cprj_work)
   end if

   ! Xout = Xout + sum_ai pai cprj_work
   call xg_nonlop_apply_prj(xg_nonlop,cprj_work,Xout,work_mpi)

   ! Xout = Xout - e Xin
   call xgBlock_yxmax(Xout,eigen,Xin)

   call timab(tim_getHmeSX,2,tsec)

 end subroutine xg_nonlop_getHmeSX
!!***
end module m_xg_nonlop
!!***
