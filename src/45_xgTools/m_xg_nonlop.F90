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
!! Copyright (C) 2022-2025 ABINIT group (LB)
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

 ! Independent timers of xg_nonlop :
 integer, parameter :: tim_getcprj    = 2101
 integer, parameter :: tim_apply_prj  = 2102
 integer, parameter :: tim_apply_Aij  = 2103
 integer, parameter :: tim_mult_cprj  = 2104
 integer, parameter :: tim_make_k     = 2105
 integer, parameter :: tim_make_Dij   = 2106
 integer, parameter :: tim_make_Sij   = 2107
 integer, parameter :: tim_make_ekb   = 2108
 integer, parameter :: tim_apply_diag = 2109
 integer, parameter :: tim_init       = 2110

 ! Timers that depend on other xg_nonlop timers :
 integer, parameter :: tim_getXSY     = 2120
 integer, parameter :: tim_getXHY     = 2121
 integer, parameter :: tim_getHmeSX   = 2122
 integer, parameter :: tim_iter_refinement = 2123

 integer, parameter :: tim_getcprj_gemm     = 2130
 integer, parameter :: tim_getcprj_copy     = 2131
 integer, parameter :: tim_getcprj_mpi      = 2132
 integer, parameter :: tim_getcprj_otf      = 2133

 integer, parameter :: tim_apply_prj_gemm   = 2135
 integer, parameter :: tim_apply_prj_copy   = 2136
 integer, parameter :: tim_apply_prj_mpi    = 2137
 integer, parameter :: tim_apply_prj_otf    = 2138

 integer, parameter :: tim_mult_cprj_gemm   = 2140
 integer, parameter :: tim_mult_cprj_copy   = 2141
 integer, parameter :: tim_mult_cprj_mpi    = 2142

 integer, parameter :: tim_forces_stress      = 2150
 integer, parameter :: tim_fst_start          = 2151
 integer, parameter :: tim_fst_cprj_deriv_f   = 2152
 integer, parameter :: tim_fst_cprj_deriv_str = 2153
 integer, parameter :: tim_fst_mult_cprj_f    = 2154
 integer, parameter :: tim_fst_mult_cprj_str  = 2155
 integer, parameter :: tim_fst_work_str       = 2156

 integer, parameter, public :: DERIV_ATOM   = 1
 integer, parameter, public :: DERIV_STRESS = 2

 type,public :: xg_nonlop_t

   integer :: cplex
   integer :: cplex_alldij
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
   integer :: space_cprj
   integer :: space_Dij
   logical :: paw
   integer :: option
   real(dp) :: weight

   integer, pointer :: mpi_atmtab(:)
   integer, pointer :: indlmn(:,:,:)
   integer, pointer :: nattyp(:)

   integer,allocatable :: l_npw_k(:)
   integer,allocatable :: l_shift_npw_k(:)

   real(dp), pointer :: sij_triangular_mat(:,:)

   integer, allocatable :: nlmn_natom(:)
   integer, allocatable :: nlmn_ntypat(:)

   real(dp), pointer :: ffnl_k(:,:,:,:)
   real(dp), pointer :: ph3d_k(:,:,:)

   real(dp), pointer :: kpg_k(:,:)

   type(xg_t),pointer :: projectors(:)
   type(xg_t),pointer :: projectors_k

   type(xg_t),pointer :: ffnl_gather(:)
   type(xg_t),pointer :: ffnl_gather_k

   type(xg_t),pointer :: ph3d_gather(:)
   type(xg_t),pointer :: ph3d_gather_k

   ! non paw only:
   type(xg_t) :: ekb
   ! end non paw only

   ! paw only:
   type(xg_t),pointer :: gram_proj(:)
   type(xg_t),pointer :: gram_proj_k

   type(xg_t) :: Dij
   type(xgBlock_t) :: Dij_spin
   type(xg_t) :: Sij

   type(xg_t) :: Sijm1
   type(xg_t), pointer :: invSij_approx(:)
   type(xg_t), pointer :: invSij_approx_k
   ! end paw only

 end type xg_nonlop_t
!!***

  ! To initialize/make/destroy xg_nonlop object
  public :: xg_nonlop_init
  public :: xg_nonlop_make_k
  public :: xg_nonlop_destroy
  public :: xg_nonlop_make_ekb    ! non paw only
  public :: xg_nonlop_destroy_ekb ! non paw only
  public :: xg_nonlop_update_weight
  public :: xg_nonlop_init_cplex_alldij ! paw only
  public :: xg_nonlop_make_Dij    ! paw only
  public :: xg_nonlop_set_Dij_spin! paw only
  public :: xg_nonlop_make_Sij    ! paw only
  public :: xg_nonlop_destroy_Dij ! paw only
  public :: xg_nonlop_destroy_Sij ! paw only
  ! Generic operations
  public :: xg_nonlop_getcprj
  public :: xg_nonlop_apply_Aij
  public :: xg_nonlop_precond_iterative_refinement
  public :: xg_nonlop_mult_cprj
  public :: xg_nonlop_apply_prj
  public :: xg_nonlop_colwiseXAX
  public :: xg_nonlop_colwiseXDX
  public :: xg_nonlop_getXAY
  public :: xg_nonlop_getXDY
  public :: xg_nonlop_getAX
  public :: xg_nonlop_getDX
  ! Specific operations (using Sij/Dij or ekb arrays)
  public :: xg_nonlop_getXHY
  public :: xg_nonlop_getHX
  public :: xg_nonlop_colwiseXHX
  public :: xg_nonlop_getXSY   ! paw only
  public :: xg_nonlop_getSX    ! paw only
  public :: xg_nonlop_getSm1X  ! paw only
  public :: xg_nonlop_getHmeSX ! paw only
  public :: xg_nonlop_forces_stress

contains
!!***

!!****f* m_xg_nonlop/xg_nonlop_init
!! NAME
!! xg_nonlop_init
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
 subroutine xg_nonlop_init(xg_nonlop,indlmn,my_natom,nattyp,mkmem,ntypat,nspinor,ucvol,usepaw,&
     xg_nonlop_option,me_band,comm_band,comm_atom,mpi_atmtab)

   integer ,intent(in) :: me_band,comm_band,comm_atom
   integer ,intent(in) :: my_natom
   integer ,intent(in) :: mkmem
   integer ,intent(in) :: ntypat
   integer ,intent(in) :: nspinor
   integer ,intent(in) :: usepaw
   integer ,intent(in) :: xg_nonlop_option
   real(dp),intent(in) :: ucvol
   type(xg_nonlop_t),intent(inout) :: xg_nonlop

   integer,optional,intent(in),target :: mpi_atmtab(:)
   integer,intent(in),target :: indlmn(:,:,:)
   integer,intent(in),target :: nattyp(:)
   real(dp) :: tsec(2)

   integer :: itypat,cprjdim,nlmn,nlmn_max,natom,shift,nmpi

   call timab(tim_init,1,tsec)

   xg_nonlop%mkmem=mkmem
   xg_nonlop%nspinor=nspinor
   xg_nonlop%comm_atom=comm_atom
   if (xmpi_comm_size(comm_atom)>1) then
     if (present(mpi_atmtab)) then
       xg_nonlop%mpi_atmtab=>mpi_atmtab
     else
       ABI_ERROR("mpi_atmtab must be present")
     end if
   else
     xg_nonlop%mpi_atmtab=>null()
   end if
   xg_nonlop%my_natom=my_natom
   xg_nonlop%me_band=me_band
   xg_nonlop%comm_band=comm_band

   if (xg_nonlop_option==0.or.xg_nonlop_option==1) then
     xg_nonlop%option = xg_nonlop_option
   else
     ABI_ERROR('Wrong value of xg_nonlop_option')
   end if

   xg_nonlop%paw=usepaw==1

   xg_nonlop%space_pw=0
   xg_nonlop%space_cprj=0
   xg_nonlop%space_Dij=0

   xg_nonlop%weight=four_pi/sqrt(ucvol)
   xg_nonlop%ntypat=ntypat

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
   ABI_MALLOC(xg_nonlop%ffnl_gather,(mkmem))
   ABI_MALLOC(xg_nonlop%ph3d_gather,(mkmem))
   if (xg_nonlop%paw) then
     ABI_MALLOC(xg_nonlop%gram_proj,(mkmem))
     ABI_MALLOC(xg_nonlop%invSij_approx,(mkmem))
   end if

  nmpi = xmpi_comm_size(xg_nonlop%comm_band)
  ABI_MALLOC(xg_nonlop%l_npw_k,(nmpi))
  ABI_MALLOC(xg_nonlop%l_shift_npw_k,(nmpi))

  call timab(tim_init,2,tsec)

 end subroutine xg_nonlop_init
!!***

!!****f* m_xg_nonlop/xg_nonlop_update_weight
!! NAME
!! xg_nonlop_init
!!
!! FUNCTION
!! Compute new weights with respect to ucvol
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
 subroutine xg_nonlop_update_weight(xg_nonlop,ucvol)

   type(xg_nonlop_t),intent(inout) :: xg_nonlop
   real(dp),intent(in) :: ucvol

   xg_nonlop%weight=four_pi/sqrt(ucvol)

 end subroutine xg_nonlop_update_weight
!!***

!!****f* m_xg_nonlop/xg_nonlop_init_cplex_alldij
!! NAME
!! xg_nonlop_init_cplex_alldij
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
 subroutine xg_nonlop_init_cplex_alldij(xg_nonlop,paw_ij)

   type(paw_ij_type),intent(in)    :: paw_ij(:)
   type(xg_nonlop_t),intent(inout) :: xg_nonlop

   logical :: paral_atom
   integer :: iatom,iatom_tot,ierr
   integer,allocatable :: l_cplex(:)
   real(dp) :: tsec(2)

   call timab(tim_init,1,tsec)

   paral_atom=(xmpi_comm_size(xg_nonlop%comm_atom)>1)

   ABI_MALLOC(l_cplex,(xg_nonlop%natom))
   l_cplex=0
   do iatom=1,xg_nonlop%my_natom ! loop over atoms treated by this proc
     iatom_tot=iatom;if (paral_atom) iatom_tot=xg_nonlop%mpi_atmtab(iatom)
     l_cplex(iatom_tot)=paw_ij(iatom)%cplex_dij
   end do
   call xmpi_sum(l_cplex,xg_nonlop%comm_atom,ierr)
   xg_nonlop%cplex_alldij = 0
   do iatom=1,xg_nonlop%natom ! loop over all atoms
     if (xg_nonlop%cplex_alldij<l_cplex(iatom)) xg_nonlop%cplex_alldij=l_cplex(iatom)
   end do
   ABI_FREE(l_cplex)

   call timab(tim_init,2,tsec)

 end subroutine xg_nonlop_init_cplex_alldij
!!***

 subroutine xg_nonlop_destroy(xg_nonlop)

  type(xg_nonlop_t),intent(inout) :: xg_nonlop

  integer :: ikpt

  if (xg_nonlop%paw) then
    call xg_nonlop_destroy_Sij(xg_nonlop) ! Can be destroyed before
    call xg_nonlop_destroy_Dij(xg_nonlop) ! Can be destroyed before
  else
    call xg_nonlop_destroy_ekb(xg_nonlop) ! Can be destroyed before
  end if

  ABI_FREE(xg_nonlop%nlmn_ntypat)
  ABI_FREE(xg_nonlop%nlmn_natom)

  ABI_FREE(xg_nonlop%l_npw_k)
  ABI_FREE(xg_nonlop%l_shift_npw_k)

  do ikpt=1,xg_nonlop%mkmem
    call xg_free(xg_nonlop%projectors(ikpt))
    call xg_free(xg_nonlop%ffnl_gather(ikpt))
    call xg_free(xg_nonlop%ph3d_gather(ikpt))
    if (xg_nonlop%paw) then
      call xg_free(xg_nonlop%gram_proj(ikpt))
      call xg_free(xg_nonlop%invSij_approx(ikpt))
    end if
  end do
  ABI_FREE(xg_nonlop%projectors)
  ABI_FREE(xg_nonlop%ffnl_gather)
  ABI_FREE(xg_nonlop%ph3d_gather)
  if (xg_nonlop%paw) then
    ABI_FREE(xg_nonlop%gram_proj)
    ABI_FREE(xg_nonlop%invSij_approx)
  end if

 end subroutine xg_nonlop_destroy
!!***

 subroutine xg_nonlop_make_ekb(xg_nonlop,ekb)

  type(xg_nonlop_t),intent(inout) :: xg_nonlop
  real(dp), intent(in) :: ekb(:,:)

  integer :: itypat, ilmn, iln, nlmn, nlmn_max, ntypat
  real(dp),pointer :: ekb_(:)
  real(dp) :: tsec(2)
  type(xgBlock_t) :: ekb_itypat

! *************************************************************************

  call timab(tim_make_ekb,1,tsec)

  if (xg_nonlop%paw) then
    ABI_ERROR('Not implemented with paw=True.')
  end if

  ntypat   = xg_nonlop%ntypat
  nlmn_max = xg_nonlop%nlmn_max

  call xg_init(xg_nonlop%ekb,SPACE_R,nlmn_max,ntypat,xmpi_comm_self)

  do itypat=1,ntypat

    nlmn=xg_nonlop%nlmn_ntypat(itypat)
    call xg_setBlock(xg_nonlop%ekb,ekb_itypat,nlmn,1,fcol=itypat)
    call xgBlock_reverseMap_1D(ekb_itypat,ekb_)
    do ilmn=1,nlmn
      iln=xg_nonlop%indlmn(5,ilmn,itypat)
      ekb_(ilmn) = ekb(iln,itypat)
    end do

  end do

  call timab(tim_make_ekb,2,tsec)

 end subroutine xg_nonlop_make_ekb
!!***

 subroutine xg_nonlop_make_Dij(xg_nonlop,paw_ij,nsppol,atindx)

  type(xg_nonlop_t),intent(inout) :: xg_nonlop
  type(paw_ij_type),intent(in)    :: paw_ij(:)
  integer,intent(in)              :: nsppol
  integer,intent(in)              :: atindx(:)

  logical :: paral_atom
  integer :: isppol, iatom, iatom_input, iatom_type, nlmn, nlmn_max, natom, nspinor, shift
  integer :: ilmn, jlmn, j0lmn, jjlmn, ijlmn
  integer :: cplex_alldij,cplex_dij,isp,isps,jsp,jsps,ijsp
  real(dp),pointer :: Dij_iatom_(:,:)
  real(dp) :: tsec(2)
  type(xgBlock_t) :: Dij_iatom

! *************************************************************************

  call timab(tim_make_Dij,1,tsec)

  if (.not.xg_nonlop%paw) then
    ABI_ERROR('Not implemented with paw=False.')
  end if

  nspinor  = xg_nonlop%nspinor
  natom    = xg_nonlop%natom
  nlmn_max = xg_nonlop%nlmn_max

  paral_atom=(xmpi_comm_size(xg_nonlop%comm_atom)>1)

  cplex_alldij = xg_nonlop%cplex_alldij
  if (cplex_alldij==1) then
    xg_nonlop%space_Dij=SPACE_R
  else if (cplex_alldij==2) then
    xg_nonlop%space_Dij=SPACE_C
  else
    ABI_ERROR('Bad cplex_alldij')
  end if

  call xg_init(xg_nonlop%Dij,xg_nonlop%space_Dij,nspinor*nlmn_max,nspinor*nlmn_max*natom*nsppol,xmpi_comm_null)

  do isppol=1,nsppol

    do iatom=1,xg_nonlop%my_natom ! loop over atoms treated by this proc
      iatom_input = iatom
      if (paral_atom) iatom_input=xg_nonlop%mpi_atmtab(iatom) ! mpi_atmtab(iatom) has the ordering of the input file
      iatom_type = atindx(iatom_input) ! convert iatom from input file to iatom ordered by type
      nlmn=xg_nonlop%nlmn_natom(iatom_type)
      shift=1+(iatom_type-1)*nspinor*nlmn_max+(isppol-1)*natom*nspinor*nlmn_max
      call xg_setBlock(xg_nonlop%Dij,Dij_iatom,nspinor*nlmn,nspinor*nlmn,fcol=shift)
      cplex_dij=paw_ij(iatom)%cplex_dij
      call xgBlock_reverseMap(Dij_iatom,Dij_iatom_)
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
              Dij_iatom_(cplex_alldij*(jlmn+isps-1)+1,jlmn+jsps) = paw_ij(iatom)%dij(jjlmn,ijsp)
              if (cplex_alldij==2) Dij_iatom_(2*(jlmn+isps),jlmn+jsps) = zero
            else
              Dij_iatom_(2*(jlmn+isps-1)+1,jlmn+jsps) = paw_ij(iatom)%dij(2*jjlmn-1,ijsp)
              Dij_iatom_(2*(jlmn+isps)    ,jlmn+jsps) = paw_ij(iatom)%dij(2*jjlmn  ,ijsp)
            end if
            do ilmn=1,jlmn-1
              ! see m_hamiltonian:pawdij2ekb and opernlc_ylm
              ijlmn=j0lmn+ilmn
              if (cplex_dij==1) then
                Dij_iatom_(cplex_alldij*(ilmn+isps-1)+1,jlmn+jsps) = paw_ij(iatom)%dij(ijlmn,ijsp)
                if (cplex_alldij==2) Dij_iatom_(2*(ilmn+isps),jlmn+jsps) = zero
                Dij_iatom_(cplex_alldij*(jlmn+jsps-1)+1,ilmn+isps) = paw_ij(iatom)%dij(ijlmn,ijsp)
                if (cplex_alldij==2) Dij_iatom_(2*(jlmn+jsps),ilmn+isps) = zero
              else
                Dij_iatom_(2*(ilmn+isps-1)+1,jlmn+jsps) = paw_ij(iatom)%dij(2*ijlmn-1,ijsp)
                Dij_iatom_(2*(ilmn+isps)    ,jlmn+jsps) = paw_ij(iatom)%dij(2*ijlmn  ,ijsp)
                Dij_iatom_(2*(jlmn+jsps-1)+1,ilmn+isps) = paw_ij(iatom)%dij(2*ijlmn-1,ijsp)
                Dij_iatom_(2*(jlmn+jsps)    ,ilmn+isps) =-paw_ij(iatom)%dij(2*ijlmn  ,ijsp)
              end if
            end do
          end do
        end do
      end do
    end do

  end do

! Communication in case of distribution over atomic sites
  if (paral_atom) then
    !call xgBlock_reverseMap(xg_nonlop%Dij%self,Dij_all_)
    !call xmpi_sum(Dij_all_,xg_nonlop%comm_atom,ierr)
    call xgBlock_mpi_sum(xg_nonlop%Dij%self,comm=xg_nonlop%comm_atom)
  end if

  call timab(tim_make_Dij,2,tsec)

 end subroutine xg_nonlop_make_Dij
!!***

 subroutine xg_nonlop_set_Dij_spin(xg_nonlop,isppol)

  integer,intent(in) :: isppol
  type(xg_nonlop_t),intent(inout) :: xg_nonlop

  integer :: shift,nlmn_max,nspinor,natom

  nspinor  = xg_nonlop%nspinor
  natom    = xg_nonlop%natom
  nlmn_max = xg_nonlop%nlmn_max

  if (isppol/=1) then ! isppol must be 1 or 2 if nspinor==1, and must be 1 of nspinor==2
    if (isppol/=2.or.nspinor/=1) then
      ABI_ERROR('wrong isppol')
    end if
  end if

  shift=1+(isppol-1)*natom*nspinor*nlmn_max
  call xg_setBlock(xg_nonlop%Dij,xg_nonlop%Dij_spin,nspinor*nlmn_max,nspinor*nlmn_max*natom,fcol=shift)

 end subroutine xg_nonlop_set_Dij_spin
!!***

 subroutine xg_nonlop_make_Sij(xg_nonlop,pawtab,inv_sij)

  type(xg_nonlop_t),intent(inout) :: xg_nonlop
  type(pawtab_type),intent(in) :: pawtab(:)
  logical,optional,intent(in) :: inv_sij

  logical :: inv_sij_
  integer :: itypat, nlmn, nlmn_max, ntypat, shift
  integer :: ilmn, jlmn, j0lmn, jjlmn, ijlmn
  real(dp),pointer :: Sij_itypat_(:,:)
  real(dp) :: tsec(2)
  type(xg_t) :: work
  type(xgBlock_t) :: Sij_itypat,Sijm1_itypat

! *************************************************************************

  call timab(tim_make_Sij,1,tsec)

  if (.not.xg_nonlop%paw) then
    ABI_ERROR('Not implemented with paw=False.')
  end if

  ntypat   = xg_nonlop%ntypat
  nlmn_max = xg_nonlop%nlmn_max

  call xg_init(xg_nonlop%Sij,SPACE_R,nlmn_max,nlmn_max*ntypat,xmpi_comm_self)

  inv_sij_ = .false.
  if (present(inv_sij)) inv_sij_ = inv_sij

  if (inv_sij_) then
    call xg_init(xg_nonlop%Sijm1,SPACE_R,nlmn_max,nlmn_max*ntypat,xmpi_comm_self)
  end if

  do itypat=1,ntypat

    nlmn=xg_nonlop%nlmn_ntypat(itypat)
    shift=1+(itypat-1)*nlmn_max
    call xg_setBlock(xg_nonlop%Sij,Sij_itypat,nlmn,nlmn,fcol=shift)
    call xgBlock_reverseMap(Sij_itypat,Sij_itypat_)
    do jlmn=1,nlmn
      j0lmn=jlmn*(jlmn-1)/2
      jjlmn=j0lmn+jlmn
      Sij_itypat_(jlmn,jlmn) = pawtab(itypat)%sij(jjlmn)
      do ilmn=1,jlmn-1
        ijlmn=j0lmn+ilmn
        Sij_itypat_(ilmn,jlmn) = pawtab(itypat)%sij(ijlmn)
        Sij_itypat_(jlmn,ilmn) = pawtab(itypat)%sij(ijlmn)
      end do
    end do
    if (inv_sij_) then
      call xg_init(work,SPACE_R,nlmn,nlmn,xmpi_comm_self)
      call xg_setBlock(xg_nonlop%Sijm1,Sijm1_itypat,nlmn,nlmn,fcol=shift)
      call xgBlock_invert_sy(Sijm1_itypat,work%self,xg_input=Sij_itypat)
      call xg_free(work)
    end if

  end do

  call timab(tim_make_Sij,2,tsec)

 end subroutine xg_nonlop_make_Sij
!!***

 subroutine xg_nonlop_destroy_ekb(xg_nonlop)

  type(xg_nonlop_t),intent(inout) :: xg_nonlop

! *************************************************************************

  if (xg_nonlop%paw) then
    ABI_ERROR('Not implemented with paw=True.')
  end if

  call xg_free(xg_nonlop%ekb)

 end subroutine xg_nonlop_destroy_ekb
!!***

 subroutine xg_nonlop_destroy_Dij(xg_nonlop)

  type(xg_nonlop_t),intent(inout) :: xg_nonlop

! *************************************************************************

  if (.not.xg_nonlop%paw) then
    ABI_ERROR('Not implemented with paw=False.')
  end if

  call xg_free(xg_nonlop%Dij)

 end subroutine xg_nonlop_destroy_Dij
!!***

 subroutine xg_nonlop_destroy_Sij(xg_nonlop)

  type(xg_nonlop_t),intent(inout) :: xg_nonlop

! *************************************************************************

  if (.not.xg_nonlop%paw) then
    ABI_ERROR('Not implemented with paw=False.')
  end if

  call xg_free(xg_nonlop%Sij)
  call xg_free(xg_nonlop%Sijm1)

 end subroutine xg_nonlop_destroy_Sij
!!***

 subroutine xg_nonlop_compute_projs(xg_nonlop)

  type(xg_nonlop_t),intent(inout) :: xg_nonlop

  complex(dp),pointer :: projectors_k_(:,:)
  real(dp),pointer :: projectors_k_real(:,:)
  real(dp),pointer :: ffnl_gather_k_(:,:)
  complex(dp),pointer :: ph3d_gather_k_(:,:)
  real(dp),pointer :: ph3d_gather_k_real(:,:)

  integer :: shift_itypat,shift_itypat_nlmn,ntypat,nattyp,shift_ipw
  integer :: icol,ilmn,nlmn,il,ipw,ia,iatom,itypat
  real(dp) :: ffnl_ipw
  complex(dp) :: cil(4),ph3d_ipw,ctmp
  logical :: compute_gather

  ntypat = xg_nonlop%ntypat

! 4pi/sqrt(ucvol) * (-i)^l
  cil(1) = ( 1.0_DP, 0.0_DP) * xg_nonlop%weight
  cil(2) = ( 0.0_DP,-1.0_DP) * xg_nonlop%weight
  cil(3) = (-1.0_DP, 0.0_DP) * xg_nonlop%weight
  cil(4) = ( 0.0_DP, 1.0_DP) * xg_nonlop%weight

  compute_gather = xg_nonlop%option==0

  if (compute_gather) then
    call xgBlock_zero(xg_nonlop%ffnl_gather_k%self)
    call xgBlock_zero(xg_nonlop%ph3d_gather_k%self)
    call xgBlock_reverseMap(xg_nonlop%ffnl_gather_k%self,ffnl_gather_k_)
  end if

  shift_ipw = xg_nonlop%l_shift_npw_k(xg_nonlop%me_band+1)

  select case(xg_nonlop%space_pw)

    case (SPACE_C)

      call xgBlock_reverseMap(xg_nonlop%projectors_k%self,projectors_k_)
      if (compute_gather) then
        call xgBlock_reverseMap(xg_nonlop%ph3d_gather_k%self,ph3d_gather_k_)
      end if
      shift_itypat=0
      shift_itypat_nlmn=0
      !$omp parallel default (none) &
      !$omp& shared(xg_nonlop,ph3d_gather_k_,ffnl_gather_k_,projectors_k_), &
      !$omp& firstprivate(compute_gather,ntypat,shift_ipw,shift_itypat,shift_itypat_nlmn,cil), &
      !$omp& private(itypat,nattyp,nlmn,ia,ilmn,ipw,iatom,il,ph3d_ipw,ffnl_ipw,icol)
      do itypat = 1, ntypat
        nlmn = xg_nonlop%nlmn_ntypat(itypat)
        nattyp = xg_nonlop%nattyp(itypat)
        !! projectors = 4pi/sqrt(ucvol) * (-i)^l * conj(ph3d) * ffnl
        !$omp do collapse(3)
        do ia = 1, nattyp
          do ilmn=1,nlmn
            do ipw=1,xg_nonlop%npw_k
              iatom = ia + shift_itypat
              il=mod(xg_nonlop%indlmn(1,ilmn,itypat),4)+1
              ph3d_ipw = cmplx( xg_nonlop%ph3d_k(1,ipw,iatom), xg_nonlop%ph3d_k(2,ipw,iatom), kind=DP)
              ffnl_ipw = xg_nonlop%ffnl_k(ipw, 1, ilmn, itypat)
              !
              if (compute_gather) then
                ph3d_gather_k_(ipw+shift_ipw,iatom) = ph3d_ipw
                ffnl_gather_k_(ipw+shift_ipw,ilmn+(itypat-1)*xg_nonlop%nlmn_max) = ffnl_ipw
              end if
              !
              icol = ilmn + (ia-1)*nlmn + shift_itypat_nlmn
              projectors_k_(ipw,icol) = cil(il) * conjg(ph3d_ipw) * ffnl_ipw
            end do
          end do
        end do
        !$omp end do
        shift_itypat      = shift_itypat      + nattyp
        shift_itypat_nlmn = shift_itypat_nlmn + nattyp*nlmn
      end do
      !$omp end parallel

    case (SPACE_CR)

      call xgBlock_reverseMap(xg_nonlop%projectors_k%self,projectors_k_real)
      if (compute_gather) then
        call xgBlock_reverseMap(xg_nonlop%ph3d_gather_k%self,ph3d_gather_k_real)
      end if
      shift_itypat=0
      shift_itypat_nlmn=0
      !$omp parallel default (none) &
      !$omp& shared(xg_nonlop,ph3d_gather_k_real,ffnl_gather_k_,projectors_k_real), &
      !$omp& firstprivate(compute_gather,ntypat,shift_ipw,shift_itypat,shift_itypat_nlmn,cil), &
      !$omp& private(itypat,nattyp,nlmn,ia,ilmn,ipw,iatom,il,ph3d_ipw,ffnl_ipw,icol,ctmp)
      do itypat = 1, ntypat
        nlmn = xg_nonlop%nlmn_ntypat(itypat)
        nattyp = xg_nonlop%nattyp(itypat)
        !! projectors = 4pi/sqrt(ucvol) * (-i)^l * conj(ph3d) * ffnl
        !$omp do collapse(3)
        do ia = 1, nattyp
          do ilmn=1,nlmn
            do ipw=1,xg_nonlop%npw_k
              iatom = ia + shift_itypat
              il=mod(xg_nonlop%indlmn(1,ilmn,itypat),4)+1
              ph3d_ipw = cmplx( xg_nonlop%ph3d_k(1,ipw,iatom), xg_nonlop%ph3d_k(2,ipw,iatom), kind=DP)
              ffnl_ipw = xg_nonlop%ffnl_k(ipw, 1, ilmn, itypat)
              !
              if (compute_gather) then
                ph3d_gather_k_real(2*(ipw+shift_ipw)-1,iatom) = dble(ph3d_ipw)
                ph3d_gather_k_real(2*(ipw+shift_ipw)  ,iatom) = dimag(ph3d_ipw)
                ffnl_gather_k_(ipw+shift_ipw,ilmn+(itypat-1)*xg_nonlop%nlmn_max) = ffnl_ipw
              end if
              !
              ctmp = cil(il) * conjg(ph3d_ipw) * ffnl_ipw
              icol = ilmn + (ia-1)*nlmn + shift_itypat_nlmn
              projectors_k_real(2*ipw-1,icol) = dble(ctmp)
              projectors_k_real(2*ipw  ,icol) = dimag(ctmp)
            end do
          end do
        end do
        !$omp end do
        shift_itypat      = shift_itypat      + nattyp
        shift_itypat_nlmn = shift_itypat_nlmn + nattyp*nlmn
      end do
      !$omp end parallel

    case default
      ABI_ERROR("Wrong space")

  end select

  if (compute_gather) then
    call xgBlock_mpi_sum(xg_nonlop%ffnl_gather_k%self,comm=xg_nonlop%comm_band)
    call xgBlock_mpi_sum(xg_nonlop%ph3d_gather_k%self,comm=xg_nonlop%comm_band)
  end if

 end subroutine xg_nonlop_compute_projs
!!***

 subroutine xg_nonlop_compute_projs_deriv_atom(xg_nonlop,projs_deriv_atom)

  type(xgBlock_t),intent(inout) :: projs_deriv_atom
  type(xg_nonlop_t),intent(in) :: xg_nonlop

  complex(dp),pointer :: projectors_k_(:,:)
  complex(dp),pointer :: projectors_deriv_atom_k_(:,:)
  real(dp),pointer :: projectors_k_real(:,:)
  real(dp),pointer :: projectors_deriv_atom_k_real(:,:)

  integer :: shift_itypat_nlmn,shift_itypat_3nlmn,ntypat,nattyp
  integer :: icol,icol_deriv,ilmn,nlmn,ipw,ia,itypat,idir
  complex(dp) :: ctmp
  real(dp) :: tmp,proj_deriv_ipw_re,proj_deriv_ipw_im

  if (.not.associated(xg_nonlop%projectors_k)) then
    ABI_ERROR('projectors_k should be associated')
  end if

  if (rows(projs_deriv_atom)/=xg_nonlop%npw_k) then
    ABI_ERROR('rows(projs_deriv_atom)/=npw_k')
  end if
  if (cols(projs_deriv_atom)/=3*xg_nonlop%cprjdim) then
    ABI_ERROR('cols(projs_deriv_atom)/=3*cprjdim')
  end if

  ntypat = xg_nonlop%ntypat

  select case(xg_nonlop%space_pw)

    case (SPACE_C)

      call xgBlock_reverseMap(xg_nonlop%projectors_k%self,projectors_k_)
      call xgBlock_reverseMap(projs_deriv_atom,projectors_deriv_atom_k_)
      shift_itypat_nlmn=0
      shift_itypat_3nlmn=0
      !$omp parallel default (none) &
      !$omp& shared(xg_nonlop,projectors_k_,projectors_deriv_atom_k_), &
      !$omp& firstprivate(ntypat,shift_itypat_nlmn,shift_itypat_3nlmn), &
      !$omp& private(itypat,nattyp,nlmn,ia,ilmn,ipw,idir,icol,icol_deriv,ctmp)
      do itypat = 1, ntypat
        nlmn = xg_nonlop%nlmn_ntypat(itypat)
        nattyp = xg_nonlop%nattyp(itypat)
        !! projectors_deriv_atom(k+G) = -i * 2pi * (k+G)_idir * projectors(k+G)
        !$omp do collapse(4)
        do ia = 1, nattyp
          do ilmn=1,nlmn
            do ipw=1,xg_nonlop%npw_k
              do idir=1,3
                icol = ilmn + (ia-1)*nlmn + shift_itypat_nlmn
                ctmp = ( 0.0_DP, -1.0_DP) * two_pi * xg_nonlop%kpg_k(ipw,idir)
                icol_deriv = ilmn + (idir-1)*nlmn + (ia-1)*3*nlmn + shift_itypat_3nlmn
                projectors_deriv_atom_k_(ipw,icol_deriv) = ctmp * projectors_k_(ipw,icol)
              end do
            end do
          end do
        end do
        !$omp end do
        shift_itypat_nlmn  = shift_itypat_nlmn  + nattyp*nlmn
        shift_itypat_3nlmn = shift_itypat_3nlmn + nattyp*3*nlmn
      end do
      !$omp end parallel

    case (SPACE_CR)

      call xgBlock_reverseMap(xg_nonlop%projectors_k%self,projectors_k_real)
      call xgBlock_reverseMap(projs_deriv_atom,projectors_deriv_atom_k_real)
      shift_itypat_nlmn=0
      shift_itypat_3nlmn=0
      !$omp parallel default (none) &
      !$omp& shared(xg_nonlop,projectors_k_real,projectors_deriv_atom_k_real), &
      !$omp& firstprivate(ntypat,shift_itypat_nlmn,shift_itypat_3nlmn), &
      !$omp& private(itypat,nattyp,nlmn,ia,ilmn,ipw,idir,icol,icol_deriv,tmp), &
      !$omp& private(proj_deriv_ipw_re,proj_deriv_ipw_im)
      do itypat = 1, ntypat
        nlmn = xg_nonlop%nlmn_ntypat(itypat)
        nattyp = xg_nonlop%nattyp(itypat)
        !! projectors_deriv_atom(k+G) = -i * 2pi * (k+G)_idir * projectors(k+G)
        !$omp do collapse(4)
        do ia = 1, nattyp
          do ilmn=1,nlmn
            do ipw=1,xg_nonlop%npw_k
              do idir=1,3
                icol = ilmn + (ia-1)*nlmn + shift_itypat_nlmn
                tmp  = - two_pi * xg_nonlop%kpg_k(ipw,idir)
                icol_deriv = ilmn + (idir-1)*nlmn + (ia-1)*3*nlmn + shift_itypat_3nlmn
                !! Re(projectors_deriv_atom) =  2pi * (k+G)_idir * Im(projectors)
                !! Im(projectors_deriv_atom) = -2pi * (k+G)_idir * Re(projectors)
                proj_deriv_ipw_re = - tmp * projectors_k_real(2*ipw  ,icol)
                proj_deriv_ipw_im =   tmp * projectors_k_real(2*ipw-1,icol)
                projectors_deriv_atom_k_real(2*ipw-1,icol_deriv) = proj_deriv_ipw_re
                projectors_deriv_atom_k_real(2*ipw  ,icol_deriv) = proj_deriv_ipw_im
              end do
            end do
          end do
        end do
        !$omp end do
        shift_itypat_nlmn  = shift_itypat_nlmn  + nattyp*nlmn
        shift_itypat_3nlmn = shift_itypat_3nlmn + nattyp*3*nlmn
      end do
      !$omp end parallel

    case default
      ABI_ERROR("Wrong space")

  end select

 end subroutine xg_nonlop_compute_projs_deriv_atom
!!***

 subroutine xg_nonlop_compute_projs_deriv_stress(xg_nonlop,projs_deriv_stress)

  type(xgBlock_t),intent(inout) :: projs_deriv_stress
  type(xg_nonlop_t),intent(in) :: xg_nonlop

  complex(dp),pointer :: projectors_k_(:,:)
  complex(dp),pointer :: projectors_deriv_stress_k_(:,:)
  real(dp),pointer :: projectors_k_real(:,:)
  real(dp),pointer :: projectors_deriv_stress_k_real(:,:)

  integer :: shift_itypat,shift_itypat_nlmn,shift_itypat_6nlmn,ntypat,nattyp
  integer :: iatom,icol_shift,icol_deriv,ilmn,nlmn,ipw,ia,itypat,idir,il
  complex(dp) :: ctmp(3),cil(4),ph3d_ipw,cipw
  real(dp) :: ffnl_ipw(3)

  if (rows(projs_deriv_stress)/=xg_nonlop%npw_k) then
    ABI_ERROR('rows(projs_deriv_atom)/=npw_k')
  end if
  if (cols(projs_deriv_stress)/=6*xg_nonlop%cprjdim) then
    ABI_ERROR('cols(projs_deriv_stress)/=6*cprjdim')
  end if

  ntypat = xg_nonlop%ntypat

! 4pi/sqrt(ucvol) * (-i)^l
  cil(1) = ( 1.0_DP, 0.0_DP) * xg_nonlop%weight
  cil(2) = ( 0.0_DP,-1.0_DP) * xg_nonlop%weight
  cil(3) = (-1.0_DP, 0.0_DP) * xg_nonlop%weight
  cil(4) = ( 0.0_DP, 1.0_DP) * xg_nonlop%weight

  select case(xg_nonlop%space_pw)

    case (SPACE_C)

      call xgBlock_reverseMap(xg_nonlop%projectors_k%self,projectors_k_)
      call xgBlock_reverseMap(projs_deriv_stress,projectors_deriv_stress_k_)
      shift_itypat=0
      shift_itypat_nlmn=0
      shift_itypat_6nlmn=0
      !$omp parallel default (none) &
      !$omp& shared(xg_nonlop,projectors_k_,projectors_deriv_stress_k_), &
      !$omp& firstprivate(cil,ntypat,shift_itypat,shift_itypat_nlmn,shift_itypat_6nlmn), &
      !$omp& private(il,iatom,itypat,nattyp,nlmn,ia,ilmn,ipw,idir,icol_shift,icol_deriv), &
      !$omp& private(ctmp,ph3d_ipw,ffnl_ipw)
      do itypat = 1, ntypat
        nlmn = xg_nonlop%nlmn_ntypat(itypat)
        nattyp = xg_nonlop%nattyp(itypat)
        !! projectors_deriv_stress(k+G)_ab = 4pi/sqrt(ucvol) * (-i)^l * conj(ph3d) * (-(k+G)_b) * d/d(K_a)[ffnl_deriv(k+G)]
        !$omp do collapse(3)
        do ia = 1, nattyp
          do ilmn=1,nlmn
            do ipw=1,xg_nonlop%npw_k
              iatom = ia + shift_itypat
              il=mod(xg_nonlop%indlmn(1,ilmn,itypat),4)+1
              ph3d_ipw = cmplx( xg_nonlop%ph3d_k(1,ipw,iatom), xg_nonlop%ph3d_k(2,ipw,iatom), kind=DP)
              ffnl_ipw(:) = xg_nonlop%ffnl_k(ipw, 2:4, ilmn, itypat)
              ctmp(:) = - cil(il) * conjg(ph3d_ipw) * xg_nonlop%kpg_k(ipw,:)
              icol_shift = ilmn + (ia-1)*6*nlmn + shift_itypat_6nlmn
              ! diagonal part
              do idir=1,3
                icol_deriv = icol_shift + (idir-1)*nlmn
                projectors_deriv_stress_k_(ipw,icol_deriv) = ctmp(idir) * ffnl_ipw(idir)
              end do
              ! off-diagonal part (which is symmetric)
              ctmp(:) = half*ctmp(:)
              icol_deriv = icol_shift + (4-1)*nlmn
              projectors_deriv_stress_k_(ipw,icol_deriv) = ctmp(2) * ffnl_ipw(3) + ctmp(3) * ffnl_ipw(2)
              icol_deriv = icol_shift + (5-1)*nlmn
              projectors_deriv_stress_k_(ipw,icol_deriv) = ctmp(1) * ffnl_ipw(3) + ctmp(3) * ffnl_ipw(1)
              icol_deriv = icol_shift + (6-1)*nlmn
              projectors_deriv_stress_k_(ipw,icol_deriv) = ctmp(1) * ffnl_ipw(2) + ctmp(2) * ffnl_ipw(1)
            end do
          end do
        end do
        !$omp end do
        shift_itypat       = shift_itypat       + nattyp
        shift_itypat_nlmn  = shift_itypat_nlmn  + nattyp*nlmn
        shift_itypat_6nlmn = shift_itypat_6nlmn + nattyp*6*nlmn
      end do
      !$omp end parallel

    case (SPACE_CR)

      call xgBlock_reverseMap(xg_nonlop%projectors_k%self,projectors_k_real)
      call xgBlock_reverseMap(projs_deriv_stress,projectors_deriv_stress_k_real)
      shift_itypat=0
      shift_itypat_nlmn=0
      shift_itypat_6nlmn=0
      !$omp parallel default (none) &
      !$omp& shared(xg_nonlop,projectors_k_real,projectors_deriv_stress_k_real), &
      !$omp& firstprivate(cil,ntypat,shift_itypat,shift_itypat_nlmn,shift_itypat_6nlmn), &
      !$omp& private(il,iatom,itypat,nattyp,nlmn,ia,ilmn,ipw,idir,icol_shift,icol_deriv), &
      !$omp& private(ctmp,cipw,ph3d_ipw,ffnl_ipw)
      do itypat = 1, ntypat
        nlmn = xg_nonlop%nlmn_ntypat(itypat)
        nattyp = xg_nonlop%nattyp(itypat)
        !! projectors_deriv_stress(k+G)_ab = 4pi/sqrt(ucvol) * (-i)^l * conj(ph3d) * (-(k+G)_b) * d/d(K_a)[ffnl_deriv(k+G)]
        !$omp do collapse(3)
        do ia = 1, nattyp
          do ilmn=1,nlmn
            do ipw=1,xg_nonlop%npw_k
              iatom = ia + shift_itypat
              il=mod(xg_nonlop%indlmn(1,ilmn,itypat),4)+1
              ph3d_ipw = cmplx( xg_nonlop%ph3d_k(1,ipw,iatom), xg_nonlop%ph3d_k(2,ipw,iatom), kind=DP)
              ffnl_ipw(:) = xg_nonlop%ffnl_k(ipw, 2:4, ilmn, itypat)
              ctmp(:) = - cil(il) * conjg(ph3d_ipw) * xg_nonlop%kpg_k(ipw,:)
              icol_shift = ilmn + (ia-1)*6*nlmn + shift_itypat_6nlmn
              ! diagonal part
              do idir=1,3
                icol_deriv = icol_shift + (idir-1)*nlmn
                projectors_deriv_stress_k_real(2*ipw-1,icol_deriv) =  dble(ctmp(idir)) * ffnl_ipw(idir)
                projectors_deriv_stress_k_real(2*ipw  ,icol_deriv) = dimag(ctmp(idir)) * ffnl_ipw(idir)
              end do
              ! off-diagonal part (which is symmetric)
              ctmp(:) = half*ctmp(:)

              icol_deriv = icol_shift + (4-1)*nlmn
              cipw = ctmp(2) * ffnl_ipw(3) + ctmp(3) * ffnl_ipw(2)
              projectors_deriv_stress_k_real(2*ipw-1,icol_deriv) =  dble(cipw)
              projectors_deriv_stress_k_real(2*ipw  ,icol_deriv) = dimag(cipw)

              icol_deriv = icol_shift + (5-1)*nlmn
              cipw = ctmp(1) * ffnl_ipw(3) + ctmp(3) * ffnl_ipw(1)
              projectors_deriv_stress_k_real(2*ipw-1,icol_deriv) =  dble(cipw)
              projectors_deriv_stress_k_real(2*ipw  ,icol_deriv) = dimag(cipw)

              icol_deriv = icol_shift + (6-1)*nlmn
              cipw = ctmp(1) * ffnl_ipw(2) + ctmp(2) * ffnl_ipw(1)
              projectors_deriv_stress_k_real(2*ipw-1,icol_deriv) =  dble(cipw)
              projectors_deriv_stress_k_real(2*ipw  ,icol_deriv) = dimag(cipw)

            end do
          end do
        end do
        !$omp end do
        shift_itypat       = shift_itypat       + nattyp
        shift_itypat_nlmn  = shift_itypat_nlmn  + nattyp*nlmn
        shift_itypat_6nlmn = shift_itypat_6nlmn + nattyp*6*nlmn
      end do
      !$omp end parallel

    case default
      ABI_ERROR("Wrong space")

  end select

 end subroutine xg_nonlop_compute_projs_deriv_stress
!!***

 subroutine xg_nonlop_compute_projs_otf(xg_nonlop,projs_otf,index_mpi)

  integer,intent(in) :: index_mpi
  type(xg_nonlop_t),intent(in) :: xg_nonlop
  type(xgBlock_t),intent(inout) :: projs_otf

  complex(dp),pointer :: projectors_k_(:,:)
  real(dp),pointer :: projectors_k_real(:,:)
  real(dp),pointer :: ffnl_gather_k_(:,:)
  complex(dp),pointer :: ph3d_gather_k_(:,:)
  real(dp),pointer :: ph3d_gather_k_real(:,:)

  integer :: nmpi,npw_k
  integer :: shift_itypat,shift_itypat_nlmn,ntypat,nattyp,shift_ipw
  integer :: icol,ilmn,nlmn,il,ipw,ia,iatom,itypat
  real(dp) :: ffnl_ipw,ph3d_ipw_r(2)
  complex(dp) :: cil(4),ph3d_ipw,ctmp

  if (xg_nonlop%option/=0) then
    ABI_ERROR('xg_nonlop%option/=0')
  end if

  nmpi = xmpi_comm_size(xg_nonlop%comm_band)
  if (index_mpi<0.or.index_mpi>nmpi-1) then
    ABI_ERROR('index_mpi should be between 0 and size(comm_band)-1')
  end if

  npw_k = xg_nonlop%l_npw_k(index_mpi+1)
  shift_ipw = xg_nonlop%l_shift_npw_k(index_mpi+1)

  if (rows(projs_otf)/=npw_k) then
    ABI_ERROR('rows(projs_otf)/=npw_k!')
  end if
  if (cols(projs_otf)/=xg_nonlop%cprjdim) then
    ABI_ERROR('cols(projs_otf)/=cprjdim!')
  end if

  ntypat = xg_nonlop%ntypat

! 4pi/sqrt(ucvol) * (-i)^l
  cil(1) = ( 1.0_DP, 0.0_DP) * xg_nonlop%weight
  cil(2) = ( 0.0_DP,-1.0_DP) * xg_nonlop%weight
  cil(3) = (-1.0_DP, 0.0_DP) * xg_nonlop%weight
  cil(4) = ( 0.0_DP, 1.0_DP) * xg_nonlop%weight

  call xgBlock_reverseMap(xg_nonlop%ffnl_gather_k%self,ffnl_gather_k_)

  select case(xg_nonlop%space_pw)

    case (SPACE_C)

      call xgBlock_reverseMap(projs_otf,projectors_k_)
      call xgBlock_reverseMap(xg_nonlop%ph3d_gather_k%self,ph3d_gather_k_)
      shift_itypat=0
      shift_itypat_nlmn=0
      !$omp parallel default (none) &
      !$omp& shared(xg_nonlop,ph3d_gather_k_,ffnl_gather_k_,projectors_k_), &
      !$omp& firstprivate(ntypat,npw_k,shift_ipw,shift_itypat,shift_itypat_nlmn,cil), &
      !$omp& private(itypat,nattyp,nlmn,ia,ilmn,ipw,iatom,il,ph3d_ipw,ffnl_ipw,icol)
      do itypat = 1, ntypat
        nlmn = xg_nonlop%nlmn_ntypat(itypat)
        nattyp = xg_nonlop%nattyp(itypat)
        !! projectors = 4pi/sqrt(ucvol)* conj(ph3d) * ffnl * (-i)^l
        !$omp do collapse(3)
        do ia = 1, xg_nonlop%nattyp(itypat)
          do ilmn=1,nlmn
            do ipw=1,npw_k
              iatom = ia + shift_itypat
              il=mod(xg_nonlop%indlmn(1,ilmn,itypat),4)+1
              ph3d_ipw = ph3d_gather_k_(ipw+shift_ipw,iatom)
              ffnl_ipw = ffnl_gather_k_(ipw+shift_ipw,ilmn+(itypat-1)*xg_nonlop%nlmn_max)
              !
              icol = ilmn + (ia-1)*nlmn + shift_itypat_nlmn
              projectors_k_(ipw,icol) = cil(il) * conjg(ph3d_ipw) * ffnl_ipw
            end do
          end do
        end do
        !$omp end do
        shift_itypat      = shift_itypat      + nattyp
        shift_itypat_nlmn = shift_itypat_nlmn + nattyp*nlmn
      end do
      !$omp end parallel

    case(SPACE_CR)

      call xgBlock_reverseMap(projs_otf,projectors_k_real)
      call xgBlock_reverseMap(xg_nonlop%ph3d_gather_k%self,ph3d_gather_k_real)
      shift_itypat=0
      shift_itypat_nlmn=0
      !$omp parallel default (none) &
      !$omp& shared(xg_nonlop,ph3d_gather_k_real,ffnl_gather_k_,projectors_k_real), &
      !$omp& firstprivate(ntypat,npw_k,shift_ipw,shift_itypat,shift_itypat_nlmn,cil), &
      !$omp& private(itypat,nattyp,nlmn,ia,ilmn,ipw,iatom,il,ph3d_ipw_r,ffnl_ipw,icol,ctmp)
      do itypat = 1, ntypat
        nlmn = xg_nonlop%nlmn_ntypat(itypat)
        nattyp = xg_nonlop%nattyp(itypat)
        !! projectors = 4pi/sqrt(ucvol)* conj(ph3d) * ffnl * (-i)^l
        !$omp do collapse(3)
        do ia = 1, nattyp
          do ilmn=1,nlmn
            do ipw=1,npw_k
              iatom = ia + shift_itypat
              il=mod(xg_nonlop%indlmn(1,ilmn,itypat),4)+1
              ph3d_ipw_r(1) = ph3d_gather_k_real(2*(ipw+shift_ipw)-1,iatom)
              ph3d_ipw_r(2) = ph3d_gather_k_real(2*(ipw+shift_ipw)  ,iatom)
              ffnl_ipw = ffnl_gather_k_(ipw+shift_ipw,ilmn+(itypat-1)*xg_nonlop%nlmn_max)
              !
              ctmp = cmplx( ph3d_ipw_r(1), ph3d_ipw_r(2), kind=DP)
              ctmp = cil(il) * conjg(ctmp) * ffnl_ipw
              icol = ilmn + (ia-1)*nlmn + shift_itypat_nlmn
              projectors_k_real(2*ipw-1,icol) = dble(ctmp)
              projectors_k_real(2*ipw  ,icol) = dimag(ctmp)
            end do
          end do
        end do
        !$omp end do
        shift_itypat      = shift_itypat      + nattyp
        shift_itypat_nlmn = shift_itypat_nlmn + nattyp*nlmn
      end do
      !$omp end parallel

    case default
      ABI_ERROR("Wrong space")

  end select

 end subroutine xg_nonlop_compute_projs_otf
!!***

 subroutine xg_nonlop_make_k(xg_nonlop,ikpt,istwf_k,me_g0,me_g0_fft,npw_k,ffnl_k,ph3d_k,kpg_k,compute_proj,&
     compute_invS_approx,compute_gram)

  type(xg_nonlop_t),intent(inout) :: xg_nonlop

  logical ,intent(in) :: compute_proj
  integer ,intent(in) :: ikpt
  integer ,intent(in) :: istwf_k
  integer, intent(in) :: me_g0
  integer, intent(in) :: me_g0_fft
  integer, intent(in) :: npw_k
  real(dp), intent(in), target :: ffnl_k(:,:,:,:)
  real(dp), intent(in), target :: ph3d_k(:,:,:)
  real(dp), intent(in), target :: kpg_k(:,:)
  logical ,optional,intent(in) :: compute_invS_approx
  logical ,optional,intent(in) :: compute_gram

  logical :: compute_gram_,compute_invS_approx_
  real(dp),pointer :: gram_proj_k_(:,:),Sijm1_(:,:)
  integer :: ierr, iblock, shift, shiftc, shift_sij, shift_itypat, itypat, ilmn, jlmn, nlmn, nlmn_max, ia
  integer :: cplex,cols,ntypat,nmpi,me_g0_loc,me_g0_fft_loc,space_cprj
  real(dp) :: tsec(2)
  type(xg_t) :: work
  type(xgBlock_t) :: projs,invSij_approx_k_itypat,Sijm1_itypat

! *************************************************************************

  call timab(tim_make_k,1,tsec)

  me_g0_loc = -1
  me_g0_fft_loc = -1
  if (istwf_k==1) then
    xg_nonlop%cplex=2
    xg_nonlop%space_pw=SPACE_C
    xg_nonlop%space_cprj=SPACE_C
  else ! istwf_k>1
    xg_nonlop%cplex=1
    xg_nonlop%space_pw=SPACE_CR
    xg_nonlop%space_cprj=SPACE_R
    me_g0_loc = 0
    me_g0_fft_loc = 0
    if (istwf_k==2.and.me_g0==1) then
      me_g0_loc = me_g0
      me_g0_fft_loc = me_g0_fft
    end if
  end if

  cplex = xg_nonlop%cplex
  space_cprj = xg_nonlop%space_cprj

  xg_nonlop%npw_k = npw_k

  xg_nonlop%l_npw_k(:) = 0
  xg_nonlop%l_npw_k(xg_nonlop%me_band+1) = npw_k
  call xmpi_sum(xg_nonlop%l_npw_k,xg_nonlop%comm_band,ierr)

  xg_nonlop%total_npw_k = sum(xg_nonlop%l_npw_k)
  xg_nonlop%max_npw_k = maxval(xg_nonlop%l_npw_k)

  xg_nonlop%l_shift_npw_k(1) = 0
  nmpi = xmpi_comm_size(xg_nonlop%comm_band)
  do iblock=2,nmpi
    xg_nonlop%l_shift_npw_k(iblock) = xg_nonlop%l_shift_npw_k(iblock-1) + xg_nonlop%l_npw_k(iblock-1)
  end do

  ntypat = xg_nonlop%ntypat
  nlmn_max = xg_nonlop%nlmn_max

  xg_nonlop%ph3d_k => ph3d_k
  xg_nonlop%ffnl_k => ffnl_k

  xg_nonlop%kpg_k => kpg_k

  xg_nonlop%projectors_k  => xg_nonlop%projectors(ikpt)
  if (xg_nonlop%option==0) then
    xg_nonlop%ffnl_gather_k => xg_nonlop%ffnl_gather(ikpt)
    xg_nonlop%ph3d_gather_k => xg_nonlop%ph3d_gather(ikpt)
  end if

  if (xg_nonlop%paw) then
    xg_nonlop%gram_proj_k => xg_nonlop%gram_proj(ikpt)
    xg_nonlop%invSij_approx_k => xg_nonlop%invSij_approx(ikpt)
  end if

  if (compute_proj) then

    call xg_init(xg_nonlop%projectors_k,xg_nonlop%space_pw,npw_k,xg_nonlop%cprjdim,&
      xg_nonlop%comm_band,me_g0=me_g0_loc)

    if (xg_nonlop%option==0) then
      call xg_init(xg_nonlop%ffnl_gather_k,SPACE_R,xg_nonlop%total_npw_k,nlmn_max*ntypat,xmpi_comm_null)
      call xg_init(xg_nonlop%ph3d_gather_k,xg_nonlop%space_pw,xg_nonlop%total_npw_k,xg_nonlop%natom,&
        xmpi_comm_null,me_g0=me_g0_fft_loc)
    end if

    call xg_nonlop_compute_projs(xg_nonlop)

    compute_invS_approx_=.false.
    if (present(compute_invS_approx)) compute_invS_approx_ = compute_invS_approx
    if (compute_invS_approx_) then

      if (.not.xg_nonlop%paw) then
        ABI_ERROR('Not implemented with paw=False.')
      end if

      !invSij_approx_k is allocated here as space_pw depends on k-point
      call xg_init(xg_nonlop%invSij_approx_k,space_cprj,nlmn_max,nlmn_max*ntypat,xmpi_comm_self)

      shift_itypat=1
      shift_sij=1
      do itypat = 1, ntypat
        nlmn = xg_nonlop%nlmn_ntypat(itypat)

        call xgBlock_setBlock(xg_nonlop%projectors_k%self,projs,npw_k,nlmn,fcol=shift_itypat)
        call xg_setBlock(xg_nonlop%invSij_approx_k,invSij_approx_k_itypat,nlmn,nlmn,fcol=shift_sij)
        call xg_setBlock(xg_nonlop%Sijm1,Sijm1_itypat,nlmn,nlmn,fcol=shift_sij)
        if (space_cprj==SPACE_R) then
          call xgBlock_copy(Sijm1_itypat,invSij_approx_k_itypat)
        else
          call xgBlock_r2c(Sijm1_itypat,invSij_approx_k_itypat,1)
        end if
        call xg_init(work,space_cprj,nlmn,nlmn,xmpi_comm_self)
        call xgBlock_gemm('t','n',1.0d0,projs,projs,0.0d0,work%self,comm=xg_nonlop%comm_band)
        call xgBlock_add(invSij_approx_k_itypat,work%self)
        call xgBlock_invert_sy(invSij_approx_k_itypat,work%self)
        call xg_free(work)
        shift_itypat = shift_itypat + nlmn*xg_nonlop%nattyp(itypat)
        shift_sij    = shift_sij + nlmn_max
      end do

    end if

    compute_gram_=.false.
    if (present(compute_gram)) compute_gram_ = compute_gram
    if (compute_gram_) then

      if (.not.xg_nonlop%paw) then
        ABI_ERROR('Not implemented with paw=False.')
      end if
      cplex=xg_nonlop%cplex
      cols = xg_nonlop%cprjdim
      call xg_init(xg_nonlop%gram_proj_k,space_cprj,cols,cols,xmpi_comm_self)
      projs = xg_nonlop%projectors_k%self
      call xgBlock_gemm('t','n',1.0d0,projs,projs,0.0d0,xg_nonlop%gram_proj_k%self,comm=xg_nonlop%comm_band)
      call xgBlock_reverseMap(xg_nonlop%gram_proj_k%self,gram_proj_k_)
      shift=0
      shiftc=0
      shift_sij=1
      do itypat = 1, ntypat
        nlmn = xg_nonlop%nlmn_ntypat(itypat)
        call xg_setBlock(xg_nonlop%Sijm1,Sijm1_itypat,nlmn,nlmn,fcol=shift_sij)
        call xgBlock_reverseMap(Sijm1_itypat,Sijm1_)
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
        shift_sij    = shift_sij + nlmn_max
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

           if (xg_nonlop%option==1) then
             call timab(tim_getcprj_mpi,1,tsec)
             call xgBlock_mpi_isend(projs,dest,tag,request,comm=xg_nonlop%comm_band)
             call timab(tim_getcprj_mpi,2,tsec)
           end if

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

           if (xg_nonlop%option==1) then
             call timab(tim_getcprj_mpi,1,tsec)
             call xgBlock_mpi_recv(work_mpi_npw,source,tag,comm=xg_nonlop%comm_band)
             call timab(tim_getcprj_mpi,2,tsec)
           else if (xg_nonlop%option==0) then
             call timab(tim_getcprj_otf,1,tsec)
             call xg_nonlop_compute_projs_otf(xg_nonlop,work_mpi_npw,source)
             call timab(tim_getcprj_otf,2,tsec)
           else
             ABI_ERROR("Wrong xg_nonlop%option")
           end if

           call timab(tim_getcprj_gemm,1,tsec)
           call xgBlock_gemm('t','n',1.0d0,work_mpi_npw,X_block,1.d0,cprjX)
           call timab(tim_getcprj_gemm,2,tsec)

           if (xg_nonlop%option==1) then
             call timab(tim_getcprj_mpi,1,tsec)
             call xmpi_wait(request,ierr)
             call timab(tim_getcprj_mpi,2,tsec)
           end if

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

!!****f* m_xg_nonlop/xg_nonlop_getcprj_deriv
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
subroutine xg_nonlop_getcprj_deriv(xg_nonlop,X,cprjX,work_mpi,option)

   type(xg_nonlop_t), intent(in)    :: xg_nonlop
   type(xgBlock_t)  , intent(in)    :: X
   integer          , intent(in)    :: option
   type(xgBlock_t)  , intent(inout) :: cprjX,work_mpi

!   real(dp) :: tsec(2)
   integer :: iblock,nmpi,npw,blocksize,nspinor,shift
   integer :: me_band,proj_size
   type(xgBlock_t) :: X_block,X_spinor
   type(xg_t) :: projs_deriv

!   call timab(tim_getcprj,1,tsec)

   npw = xg_nonlop%npw_k
   nspinor = xg_nonlop%nspinor
   blocksize = cols(cprjX)

   if (option==DERIV_ATOM) then
     proj_size = 3*xg_nonlop%cprjdim
   else if (option==DERIV_STRESS) then
     proj_size = 6*xg_nonlop%cprjdim
   else
     ABI_ERROR('Bad option')
   end if
   ! Check cprj sizes
   if (rows(cprjX)/=proj_size) then
     ABI_ERROR('rows(cprjX)/=proj_size')
   end if
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
     ! Linalg representation (X have all cols, rows are distributed)
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

   call xg_init(projs_deriv,xg_nonlop%space_pw,npw,proj_size,comm=xg_nonlop%comm_band,me_g0=me_g0(X))

   select case (option)
     case (DERIV_ATOM)
       call xg_nonlop_compute_projs_deriv_atom(xg_nonlop,projs_deriv%self)
     case (DERIV_STRESS)
       call xg_nonlop_compute_projs_deriv_stress(xg_nonlop,projs_deriv%self)
   end select

   call xgBlock_reshape_spinor(X,X_spinor,nspinor,ROWS2COLS)

   if (nmpi==1) then

!     call timab(tim_getcprj_gemm,1,tsec)
     call xgBlock_gemm('t','n',1.0d0,projs_deriv%self,X_spinor,0.d0,cprjX)
!     call timab(tim_getcprj_gemm,2,tsec)

   else

     me_band=xg_nonlop%me_band

     call xgBlock_check(cprjX,work_mpi)

     do iblock=1,nmpi
       shift=1+(iblock-1)*blocksize
       call xgBlock_setBlock(X_spinor,X_block,rows(X_spinor),blocksize,fcol=shift)

!       call timab(tim_getcprj_gemm,1,tsec)
       call xgBlock_gemm('t','n',1.0d0,projs_deriv%self,X_block,0.d0,work_mpi)
!       call timab(tim_getcprj_gemm,2,tsec)
       ! We do the mpi sum outside xgBlock_gemm just to include the timing in tim_getcprj_mpi,
       ! (instead of tim_gemm_mpi).
!       call timab(tim_getcprj_mpi,1,tsec)
       call xgBlock_mpi_sum(work_mpi,comm=xg_nonlop%comm_band)
!       call timab(tim_getcprj_mpi,2,tsec)

       if (me_band==iblock-1) then
!         call timab(tim_getcprj_copy,1,tsec)
         call xgBlock_copy(work_mpi,cprjX)
!         call timab(tim_getcprj_copy,2,tsec)
       end if
     end do

   end if

   call xg_free(projs_deriv)
!   call timab(tim_getcprj,2,tsec)

 end subroutine xg_nonlop_getcprj_deriv
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

           if (xg_nonlop%option==1) then
             call timab(tim_apply_prj_mpi,1,tsec)
             call xgBlock_mpi_isend(projs,dest,tag,request,comm=xg_nonlop%comm_band)
             call timab(tim_apply_prj_mpi,2,tsec)
           end if

           source = mod(me_band+(iblock-1),nmpi)
           npw = xg_nonlop%l_npw_k(source+1)
           call xg_setBlock(work_npw,X_block,npw,blocksize)
           shift_row = xg_nonlop%l_shift_npw_k(source+1)

           call timab(tim_apply_prj_copy,1,tsec)
           call xgBlock_partialcopy(X_spinor,X_block,shift_row,0,BIG2SMALL)
           call timab(tim_apply_prj_copy,2,tsec)

           call xgBlock_setBlock(work_mpi,work_mpi_npw,npw_max,xg_nonlop%cprjdim)
           call xgBlock_free_reshape(work_mpi_npw,npw,xg_nonlop%cprjdim)

           if (xg_nonlop%option==1) then
             call timab(tim_apply_prj_mpi,1,tsec)
             call xgBlock_mpi_recv(work_mpi_npw,source,tag,comm=xg_nonlop%comm_band)
             call timab(tim_apply_prj_mpi,2,tsec)
           else if (xg_nonlop%option==0) then
             call timab(tim_apply_prj_otf,1,tsec)
             call xg_nonlop_compute_projs_otf(xg_nonlop,work_mpi_npw,source)
             call timab(tim_apply_prj_otf,2,tsec)
           else
             ABI_ERROR("Wrong xg_nonlop%option")
           end if

           call timab(tim_apply_prj_gemm,1,tsec)
           call xgBlock_gemm('n','n',1.0d0,work_mpi_npw,cprjX,1.d0,X_block)
           call timab(tim_apply_prj_gemm,2,tsec)

           call timab(tim_apply_prj_copy,1,tsec)
           call xgBlock_partialcopy(X_block,X_spinor,shift_row,0,SMALL2BIG)
           call timab(tim_apply_prj_copy,2,tsec)

           if (xg_nonlop%option==1) then
             call timab(tim_apply_prj_mpi,1,tsec)
             call xmpi_wait(request,ierr)
             call timab(tim_apply_prj_mpi,2,tsec)
           end if
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

 subroutine xg_nonlop_apply_diag(xg_nonlop,diag_op,cprjin,cprjout)

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: diag_op
   type(xgBlock_t), intent(in) :: cprjin
   type(xgBlock_t), intent(inout) :: cprjout

   logical :: loop_over_atoms
   integer :: ia, iband, cprjdim, shift_itypat, iatom, itypat, nattyp, nlmn, shift
   integer :: space_cprj, cplex, nlmn_max
   integer :: nspinor, nrows, ncols

   type(xg_t)        :: cprjin_nlmn_max,cprjout_nlmn_max
   type(xgBlock_t) :: cprjin_nlmn,cprjout_nlmn,diag_op_iatom
   real(dp),pointer :: cprjin_nlmn_(:,:),cprjout_nlmn_(:,:)
   real(dp),pointer :: cprjin_(:,:),cprjout_(:,:)
   real(dp) :: tsec(2)

   call timab(tim_apply_diag,1,tsec)

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

   nspinor=xg_nonlop%nspinor
   nlmn_max=xg_nonlop%nlmn_max

   ! Create work spaces for cprj of ONE atom for ALL bands
   call xg_init(cprjin_nlmn_max ,space_cprj,nlmn_max,ncols)
   call xg_init(cprjout_nlmn_max,space_cprj,nlmn_max,ncols)

   call xgBlock_reverseMap(cprjin ,cprjin_ )
   call xgBlock_reverseMap(cprjout,cprjout_)

   if (cols(diag_op) == xg_nonlop%ntypat) then
     loop_over_atoms = .false.
   else if (cols(diag_op) == xg_nonlop%natom) then
     loop_over_atoms = .true.
   else
     ABI_ERROR('wrong cols for diag_op!')
   end if
   shift = 0
   shift_itypat = 0
   do itypat=1,xg_nonlop%ntypat
     nlmn=xg_nonlop%nlmn_ntypat(itypat)
     nattyp=xg_nonlop%nattyp(itypat)
     call xg_setBlock(cprjin_nlmn_max ,cprjin_nlmn ,nlmn,ncols)
     call xg_setBlock(cprjout_nlmn_max,cprjout_nlmn,nlmn,ncols)
     call xgBlock_reverseMap(cprjin_nlmn ,cprjin_nlmn_ )
     call xgBlock_reverseMap(cprjout_nlmn,cprjout_nlmn_)
     if (.not.loop_over_atoms) then
       call xgBlock_setBlock(diag_op,diag_op_iatom,nlmn,1,fcol=itypat)
     end if
     do ia=1,nattyp
       if (loop_over_atoms) then
         iatom = ia + shift_itypat
         call xgBlock_setBlock(diag_op,diag_op_iatom,nlmn,1,fcol=iatom)
       end if
       ! Copy cprj of ONE atom for ALL bands from cprjin to cprin_nlmn
       do iband=1,ncols
         cprjin_nlmn_(1:cplex*nlmn,iband) = cprjin_(1+shift:cplex*nlmn+shift,iband)
       end do
       call xgBlock_apply_diag(cprjin_nlmn,diag_op_iatom,1,Y=cprjout_nlmn)
       do iband=1,ncols
         cprjout_(1+shift:cplex*nlmn+shift,iband) = cprjout_(1+shift:cplex*nlmn+shift,iband) &
         & + cprjout_nlmn_(1:cplex*nlmn,iband)
       end do
       shift=shift+cplex*nlmn
     end do

     shift_itypat = shift_itypat + nattyp

   end do

   call xg_free(cprjin_nlmn_max)
   call xg_free(cprjout_nlmn_max)

   call timab(tim_apply_diag,2,tsec)

 end subroutine xg_nonlop_apply_diag
!!***

 subroutine xg_nonlop_apply_Aij(xg_nonlop,Aij,cprjin,cprjout,A_with_spin)

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: Aij
   type(xgBlock_t), intent(in) :: cprjin
   type(xgBlock_t), intent(inout) :: cprjout
   logical,optional,intent(in) :: A_with_spin

   logical :: loop_over_atoms
   integer :: ia, iband, cprjdim, shift_itypat, iatom, itypat, nattyp, nlmn, shift
   integer :: space_aij, space_cprj, cplex, nlmn_max
   integer :: nspinor, nrows, ncols, nrows_A, ncols_A
   integer :: nlmn_1atom,nlmn_max_1atom,ncols_1atom

   type(xg_t)        :: cprjin_nlmn_max,cprjout_nlmn_max
   type(xg_t),target :: Aij_complex
   type(xgBlock_t) :: cprjin_nlmn,cprjout_nlmn
   type(xgBlock_t) :: Aij_iatom,Aij_iatom_
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
   space_aij = space(Aij)
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

   if (aij_r2c) then
     if (rows(Aij) /= nlmn_max) then
       ABI_ERROR('wrong rows for Aij')
     end if
     if (cols(Aij) == xg_nonlop%ntypat*nlmn_max) then
       loop_over_atoms = .false.
     else if (cols(Aij) == xg_nonlop%natom*nlmn_max) then
       loop_over_atoms = .true.
     else
       ABI_ERROR('wrong cols for Aij (aij_r2c)')
     end if
   else
     if (rows(Aij) /= nlmn_max_1atom) then
       ABI_ERROR('wrong rows for Aij')
     end if
     if (cols(Aij) == xg_nonlop%ntypat*nlmn_max_1atom) then
       loop_over_atoms = .false.
     else if (cols(Aij) == xg_nonlop%natom*nlmn_max_1atom) then
       loop_over_atoms = .true.
     else
       ABI_ERROR('wrong cols for Aij')
     end if
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
     if (.not.loop_over_atoms) then
       if (aij_r2c) then
         call xgBlock_setBlock(Aij,Aij_iatom,nlmn,nlmn,fcol=1+(itypat-1)*nlmn_max)
       else
         call xgBlock_setBlock(Aij,Aij_iatom,nlmn_1atom,nlmn_1atom,fcol=1+(itypat-1)*nlmn_max_1atom)
       end if
     end if
     do ia=1,nattyp
       if (loop_over_atoms) then
         iatom = ia + shift_itypat
         if (aij_r2c) then
           call xgBlock_setBlock(Aij,Aij_iatom,nlmn,nlmn,fcol=1+(iatom-1)*nlmn_max)
         else
           call xgBlock_setBlock(Aij,Aij_iatom,nlmn_1atom,nlmn_1atom,fcol=1+(iatom-1)*nlmn_max_1atom)
         end if
       end if
       call xgBlock_getsize(Aij_iatom,nrows_A,ncols_A)
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
         call xgBlock_r2c(Aij_iatom,Aij_complex%self,nspinor)
         Aij_iatom_ = Aij_complex%self
       else
         Aij_iatom_ = Aij_iatom
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
       !call xgBlock_gemm('n','n',1.0d0,Aij_iatom_,cprjin_nlmn,0.d0,cprjout_nlmn,timing=.false.)
       call xgBlock_gemm('n','n',1.0d0,Aij_iatom_,cprjin_nlmn,0.d0,cprjout_nlmn)
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
   type(xgBlock_t), intent(in   ) :: A,precond
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
   type(xgBlock_t), intent(in) :: cprj_left,cprj_right
   type(xgBlock_t), intent(inout) :: res
   integer, intent(in),optional :: blocksize

   integer :: space_res
   integer :: blocksize_,blocksize_spinor,iblock_mpi,nblocks_mpi,shift_row,shift_col,shift_col_mpi
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

   blocksize_ = ncols_r
   if (present(blocksize)) then
     if (mod(blocksize,nspinor)/=0) then
       ABI_ERROR("wrong blocksize (nspinor)")
     end if
     if (mod(ncols_l,blocksize / nspinor)/=0) then
       ABI_ERROR("wrong blocksize")
     end if
     blocksize_ = blocksize
   end if

   if (nblocks_mpi==1) then

     call timab(tim_mult_cprj_gemm,1,tsec)
     call xgBlock_gemm('t','n',1.0d0,cprj_left_spinor,cprj_right_spinor,1.d0,res)
     call timab(tim_mult_cprj_gemm,2,tsec)

   else

     blocksize_spinor = blocksize_ / nspinor
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
   type(xgBlock_t), intent(in) :: cprj,Aij
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
   call xgBlock_colwiseDotProduct(cprj_spinor,cprj_work_spinor,res,comm_loc=xmpi_comm_null)

end subroutine xg_nonlop_colwiseXAX
!!***

subroutine xg_nonlop_colwiseXDX(xg_nonlop,diag,cprj,cprj_work,res)

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: cprj,diag
   type(xgBlock_t), intent(inout) :: cprj_work,res

   integer :: ncols,space_diag,space_res
   type(xgBlock_t) :: cprj_spinor,cprj_work_spinor
   type(xg_t) :: res_complex

   call xgBlock_check(cprj,cprj_work)
   ncols = cols(cprj)
   if (ncols/=xg_nonlop%nspinor*rows(res)) then
     ABI_ERROR('Wrong cols for cprj or res.')
   end if
   space_diag = space(diag)
   space_res = space(res)
   if (space_diag==SPACE_C) then
     if (space_res/=SPACE_C) then
       ABI_ERROR('space(res) should be SPACE_C.')
     end if
   else if (space_diag/=SPACE_R) then
     ABI_ERROR('space(diag) should be SPACE_C or SPACE_R.')
   end if

   call xgBlock_zero(cprj_work)
   call xg_nonlop_apply_diag(xg_nonlop,diag,cprj,cprj_work)

   call xgBlock_reshape_spinor(cprj     ,cprj_spinor     ,xg_nonlop%nspinor,COLS2ROWS)
   call xgBlock_reshape_spinor(cprj_work,cprj_work_spinor,xg_nonlop%nspinor,COLS2ROWS)

   ! If space_diag==SPACE_R, the result is actually real and can be stored in a xgBlock with SPACE_R
   if ( space_diag==SPACE_R .and. space_res==SPACE_R .and. space(cprj)==SPACE_C) then
     call xg_init(res_complex,SPACE_C,rows(res),cols(res))
     call xgBlock_colwiseDotProduct(cprj_spinor,cprj_work_spinor,res_complex%self,comm_loc=xmpi_comm_null)
     call xgBlock_c2r(res_complex%self,res)
     call xg_free(res_complex)
   else
     call xgBlock_colwiseDotProduct(cprj_spinor,cprj_work_spinor,res,comm_loc=xmpi_comm_null)
   end if

end subroutine xg_nonlop_colwiseXDX
!!***

subroutine xg_nonlop_colwiseXHX(xg_nonlop,cprj,cprj_work,res)

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: cprj
   type(xgBlock_t), intent(inout) :: cprj_work,res

   integer :: nmpi,ncols,nres,shift
!   type(xg_t) :: res_mpi
   type(xgBlock_t) :: res_mpi,res_tmp

   if (cols(res)/=1) then
     ABI_ERROR('cols(res)/=1')
   end if
   nmpi = xmpi_comm_size(comm(cprj))
   ncols = cols(cprj)
   nres = rows(res)
   if (nres*xg_nonlop%nspinor/=nmpi*ncols) then
     ABI_ERROR('rows(res)*nspinor/=nmpi*cols(cprj))')
   end if

   call xgBlock_zero(res)

   if (nmpi==1) then
     res_mpi = res
   else
     call xgBlock_setBlock(res,res_tmp,nres,1)
     call xgBlock_reshape(res_tmp,1,nres)
     shift = xg_nonlop%me_band*ncols
     call xgBlock_setBlock(res_tmp,res_mpi,1,ncols,fcol=1+shift)
     call xgBlock_reshape(res_mpi,ncols,1)
   end if

   if (xg_nonlop%paw) then
     call xg_nonlop_colwiseXAX(xg_nonlop,xg_nonlop%Dij_spin,cprj,cprj_work,res_mpi)
   else
     call xg_nonlop_colwiseXDX(xg_nonlop,xg_nonlop%ekb%self,cprj,cprj_work,res_mpi)
   end if

   call xgBlock_mpi_sum(res,comm=xg_nonlop%comm_band)

end subroutine xg_nonlop_colwiseXHX
!!***

subroutine xg_nonlop_getXAY(xg_nonlop,Aij,cprj_left,cprj_right,cprj_work,res,blocksize)

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: cprj_left,cprj_right,Aij
   type(xgBlock_t), intent(inout) :: cprj_work,res
   integer,intent(in),optional :: blocksize

   integer :: blocksize_

   call xgBlock_zero(cprj_work)
   call xg_nonlop_apply_Aij(xg_nonlop,Aij,cprj_right,cprj_work)

   blocksize_ = cols(cprj_right)
   if (present(blocksize)) then
     blocksize_ = blocksize
   end if
   call xg_nonlop_mult_cprj(xg_nonlop,cprj_left,cprj_work,res,blocksize=blocksize_)

 end subroutine xg_nonlop_getXAY
!!***

subroutine xg_nonlop_getXDY(xg_nonlop,diag,cprj_left,cprj_right,cprj_work,res,blocksize)

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: cprj_left,cprj_right,diag
   type(xgBlock_t), intent(inout) :: cprj_work,res
   integer,intent(in),optional :: blocksize

   integer :: blocksize_

   call xgBlock_zero(cprj_work)
   call xg_nonlop_apply_diag(xg_nonlop,diag,cprj_right,cprj_work)

   blocksize_ = cols(cprj_right)
   if (present(blocksize)) then
     blocksize_ = blocksize
   end if
   call xg_nonlop_mult_cprj(xg_nonlop,cprj_left,cprj_work,res,blocksize=blocksize_)

 end subroutine xg_nonlop_getXDY
!!***

subroutine xg_nonlop_getXSY(xg_nonlop,cprj_left,cprj_right,cprj_work,res,blocksize)

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(inout) :: cprj_left,cprj_right,cprj_work,res
   integer,intent(in),optional :: blocksize

   integer :: blocksize_

   real(dp) :: tsec(2)

   call timab(tim_getXSY,1,tsec)

   if (.not.xg_nonlop%paw) then
     ABI_ERROR('Not implemented with paw=False.')
   end if

   blocksize_ = cols(cprj_right)
   if (present(blocksize)) then
     blocksize_ = blocksize
   end if
   call xg_nonlop_getXAY(xg_nonlop,xg_nonlop%Sij%self,cprj_left,cprj_right,cprj_work,res,blocksize=blocksize_)

   call timab(tim_getXSY,2,tsec)

 end subroutine xg_nonlop_getXSY
!!***

subroutine xg_nonlop_getXHY(xg_nonlop,cprj_left,cprj_right,cprj_work,res,blocksize)

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(inout) :: cprj_left,cprj_right,cprj_work,res
   integer,intent(in),optional :: blocksize

   integer :: blocksize_

   real(dp) :: tsec(2)

   call timab(tim_getXHY,1,tsec)

   blocksize_ = cols(cprj_right)
   if (present(blocksize)) then
     blocksize_ = blocksize
   end if
   if (xg_nonlop%paw) then
     call xg_nonlop_getXAY(xg_nonlop,xg_nonlop%Dij_spin,cprj_left,cprj_right,cprj_work,res,blocksize=blocksize_)
   else
     call xg_nonlop_getXDY(xg_nonlop,xg_nonlop%ekb%self,cprj_left,cprj_right,cprj_work,res,blocksize=blocksize_)
   end if

   call timab(tim_getXHY,2,tsec)

 end subroutine xg_nonlop_getXHY
!!***

 subroutine xg_nonlop_getAX(xg_nonlop,Aij,Xin,cprjin,cprj_work,work_mpi,Xout)

   use iso_c_binding

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: cprjin,Aij
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

 subroutine xg_nonlop_getDX(xg_nonlop,diag,Xin,cprjin,cprj_work,work_mpi,Xout)

   use iso_c_binding

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: cprjin,diag
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
   call xg_nonlop_apply_diag(xg_nonlop,diag,cprjin,cprj_work)

   ! Xout = Xin + sum_ai pai cprj_work
   if (present(Xout)) then
     call xgBlock_copy(Xin,Xout)
     call xg_nonlop_apply_prj(xg_nonlop,cprj_work,Xout,work_mpi)
   else ! in place version
     call xg_nonlop_apply_prj(xg_nonlop,cprj_work,Xin,work_mpi)
   end if

 end subroutine xg_nonlop_getDX

 subroutine xg_nonlop_getHX(xg_nonlop,Xin,cprjin,cprj_work,work_mpi,Xout)

   use iso_c_binding

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: cprjin
   type(xgBlock_t), intent(inout) :: Xin,cprj_work,work_mpi
   type(xgBlock_t), optional, intent(inout) :: Xout

   if (present(Xout)) then
     if (xg_nonlop%paw) then
       call xg_nonlop_getAX(xg_nonlop,xg_nonlop%Dij_spin,Xin,cprjin,cprj_work,work_mpi,Xout=Xout)
     else
       call xg_nonlop_getDX(xg_nonlop,xg_nonlop%ekb%self,Xin,cprjin,cprj_work,work_mpi,Xout=Xout)
     end if
   else
     if (xg_nonlop%paw) then
       call xg_nonlop_getAX(xg_nonlop,xg_nonlop%Dij_spin,Xin,cprjin,cprj_work,work_mpi)
     else
       call xg_nonlop_getDX(xg_nonlop,xg_nonlop%ekb%self,Xin,cprjin,cprj_work,work_mpi)
     end if
   end if

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

   call xg_nonlop_getAX(xg_nonlop,xg_nonlop%Sij%self,Xin,cprjin,cprj_work,work_mpi,Xout)

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
   if (.not.associated(xg_nonlop%invSij_approx_k)) then
     ABI_ERROR('invSij_approx not associated')
   end if
   if (.not.associated(xg_nonlop%gram_proj_k)) then
     ABI_ERROR('gram_proj_k should be associated')
   end if

   call xg_nonlop_precond_iterative_refinement(xg_nonlop,xg_nonlop%gram_proj_k%self,xg_nonlop%invSij_approx_k%self,&
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
   type(xgBlock_t) :: cprj_work_spinor

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
   call xg_nonlop_apply_Aij(xg_nonlop,xg_nonlop%Sij%self,cprjin,cprj_work)

   ! cprj_work = - e cprj_work = -e sum_j Saij cprjin
   shift = xg_nonlop%me_band*ncols_cprj/xg_nonlop%nspinor
   call xgBlock_reshape_spinor(cprj_work,cprj_work_spinor,nspinor,COLS2ROWS)
   call xgBlock_ymax(cprj_work_spinor,eigen,shift,nblocks)

   no_H_=.False.
   if (present(no_H)) then
     no_H_ = no_H
   end if
   if (.not.no_H_) then
     ! cprj_work = sum_j Daij cprjin + cprj_work
     call xg_nonlop_apply_Aij(xg_nonlop,xg_nonlop%Dij_spin,cprjin,cprj_work)
   end if

   ! Xout = Xout + sum_ai pai cprj_work
   call xg_nonlop_apply_prj(xg_nonlop,cprj_work,Xout,work_mpi)

   ! Xout = Xout - e Xin
   call xgBlock_yxmax(Xout,eigen,Xin)

   call timab(tim_getHmeSX,2,tsec)

 end subroutine xg_nonlop_getHmeSX
!!***

!!****f* m_xg_nonlop/xg_nonlop_forces
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
subroutine xg_nonlop_forces_stress(xg_nonlop,Xin,cprjin,cprj_work,eigen,forces,stress,gprimd)

   use m_geometry, only : strconv

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: Xin,eigen
   type(xgBlock_t), intent(in) :: cprjin
   type(xgBlock_t), intent(inout) :: cprj_work
   type(xgBlock_t), optional, intent(inout) :: forces
   type(xgBlock_t), optional, intent(inout) :: stress
   real(dp), optional, intent(in) :: gprimd(:,:)

   real(dp) :: tsec(2)
   integer :: cplex,iband
   integer :: nmpi,shift
   integer :: nrows,ncols,ncols_cprj,ncols_cprj_nospin
   integer :: nrows_diag,ncols_diag
   integer :: nspinor,space_cprj

   type(xg_t) :: cprj_deriv,work_mpi_deriv,dot,dot_all
   real(dp), pointer :: dot_(:,:),dot_all_(:,:)
   logical :: do_forces,do_stress
   real(dp) :: work(6)
   real(dp), pointer :: stress_(:,:)
   type(xgBlock_t) :: cprjin_spinor,cprj_work_spinor

   call timab(tim_forces_stress,1,tsec)

   call timab(tim_fst_start,1,tsec)

   if (comm(Xin)/=xg_nonlop%comm_band) then
     ABI_ERROR('wrong communicator')
   end if

   nmpi = xmpi_comm_size(comm(Xin))
   nspinor = xg_nonlop%nspinor

   call xgBlock_getsize(Xin,nrows,ncols)
   call xgBlock_getsize(eigen,nrows_diag,ncols_diag)
   if (ncols/=nrows_diag) then
     ABI_ERROR('wrong ncols')
   end if
   if (ncols_diag/=1) then
     ABI_ERROR('ncols_diag should be one')
   end if

   ncols_cprj = ncols*nspinor/nmpi
   ncols_cprj_nospin = ncols/nmpi

   space_cprj = xg_nonlop%space_cprj

   do_forces = .false.
   if (present(forces)) then
     do_forces = .true.
     if (rows(forces)/=3*xg_nonlop%natom) then
       ABI_ERROR('rows(forces)/=3*natom')
     end if
     if (cols(forces)/=ncols) then
       ABI_ERROR('cols(forces)/=ncols')
     end if
   end if

   do_stress = .false.
   if (present(stress)) then
     do_stress = .true.
     if (.not.present(gprimd)) then
       ABI_ERROR('If stress is present, gprimd must be present too')
     end if
     if (rows(stress)/=6) then
       ABI_ERROR('rows(stress)/=6')
     end if
     if (cols(stress)/=ncols) then
       ABI_ERROR('cols(stress)/=ncols')
     end if
   end if

   call xgBlock_zero(cprj_work)
   if (xg_nonlop%paw) then
     ! cprj_work = sum_j Saij cprjin
     call xg_nonlop_apply_Aij(xg_nonlop,xg_nonlop%Sij%self,cprjin,cprj_work)
     ! cprj_work = - e cprj_work = -e sum_j Saij cprjin
     shift = xg_nonlop%me_band*ncols_cprj_nospin
     call xgBlock_reshape_spinor(cprj_work,cprj_work_spinor,nspinor,COLS2ROWS)
     call xgBlock_ymax(cprj_work_spinor,eigen,shift,nmpi)
   end if

   ! cprj_work = sum_j Daij cprjin + cprj_work
   if (xg_nonlop%paw) then
     call xg_nonlop_apply_Aij(xg_nonlop,xg_nonlop%Dij_spin,cprjin,cprj_work)
   else
     call xg_nonlop_apply_diag(xg_nonlop,xg_nonlop%ekb%self,cprjin,cprj_work)
   end if

   call timab(tim_fst_start,2,tsec)

   if (do_forces) then

     call xg_init(cprj_deriv    ,space_cprj,3*xg_nonlop%cprjdim,ncols_cprj,comm=xg_nonlop%comm_band)
     call xg_init(work_mpi_deriv,space_cprj,3*xg_nonlop%cprjdim,ncols_cprj,comm=xg_nonlop%comm_band)

     call timab(tim_fst_cprj_deriv_f,1,tsec)
     call xg_nonlop_getcprj_deriv(xg_nonlop,Xin,cprj_deriv%self,work_mpi_deriv%self,DERIV_ATOM)
     call timab(tim_fst_cprj_deriv_f,2,tsec)

     call xg_free(work_mpi_deriv)

     call xgBlock_zero(forces)

     call timab(tim_fst_mult_cprj_f,1,tsec)
     call xg_nonlop_mult_cprj_forces(xg_nonlop,cprj_work,cprj_deriv%self,forces)
     call timab(tim_fst_mult_cprj_f,2,tsec)

     call xg_free(cprj_deriv)

   end if

   if (do_stress) then

     call xg_init(cprj_deriv    ,space_cprj,6*xg_nonlop%cprjdim,ncols_cprj,comm=xg_nonlop%comm_band)
     call xg_init(work_mpi_deriv,space_cprj,6*xg_nonlop%cprjdim,ncols_cprj,comm=xg_nonlop%comm_band)

     call timab(tim_fst_cprj_deriv_str,1,tsec)
     call xg_nonlop_getcprj_deriv(xg_nonlop,Xin,cprj_deriv%self,work_mpi_deriv%self,DERIV_STRESS)
     call timab(tim_fst_cprj_deriv_str,2,tsec)

     call xg_free(work_mpi_deriv)

     call xgBlock_zero(stress)

     call timab(tim_fst_mult_cprj_str,1,tsec)
     call xg_nonlop_mult_cprj_stress(xg_nonlop,cprj_work,cprj_deriv%self,stress)
     call timab(tim_fst_mult_cprj_str,2,tsec)

     call xg_free(cprj_deriv)

     call timab(tim_fst_work_str,1,tsec)
     call xg_init(dot,space_cprj,ncols_cprj_nospin,1)

     call xgBlock_reshape_spinor(cprjin,cprjin_spinor,nspinor,COLS2ROWS)
     call xgBlock_reshape_spinor(cprj_work,cprj_work_spinor,nspinor,COLS2ROWS)
     call xgBlock_colwiseDotProduct(cprjin_spinor,cprj_work_spinor,dot%self,comm_loc=xmpi_comm_null)

     cplex=1
     if (space_cprj==SPACE_C) cplex=2

     call xg_init(dot_all,space_cprj,ncols,1)
     if (xmpi_comm_size(xg_nonlop%comm_band)>1) then
       call xgBlock_zero(dot_all%self)
       call xgBlock_reverseMap(dot%self,dot_)
       call xgBlock_reverseMap(dot_all%self,dot_all_)
       shift = cplex*xg_nonlop%me_band*ncols_cprj_nospin
       do iband=1,ncols_cprj_nospin
         dot_all_(1+cplex*(iband-1)+shift,1) = dot_(1+cplex*(iband-1),1)
       end do
       call xgBlock_mpi_sum(dot_all%self,comm=xg_nonlop%comm_band)
     else
       call xgBlock_reverseMap(dot%self,dot_all_)
     end if

     call xgBlock_reverseMap(stress,stress_)

     do iband=1,ncols
       work = stress_(:,iband)
       call strconv(work,gprimd,work)
       stress_(1:3,iband) = work(1:3) - dot_all_(1+cplex*(iband-1),1)
       stress_(4:6,iband) = work(4:6)
     end do

     call xg_free(dot)
     call xg_free(dot_all)

     call timab(tim_fst_work_str,2,tsec)

   end if

   call timab(tim_forces_stress,2,tsec)

 end subroutine xg_nonlop_forces_stress
!!***

 subroutine xg_nonlop_mult_cprj_forces(xg_nonlop,cprj,cprj_deriv,forces)

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: cprj,cprj_deriv
   type(xgBlock_t), intent(inout) :: forces

   !real(dp) :: tsec(2)
   integer :: ia,idir,ilmn,iband,iband_spinor,my_iband,itypat,nlmn,nattyp
   integer :: ispinor,iforces,icprj,icprj_deriv
   integer :: ncols_cprj,ncols_cprj_nospin
   integer :: nspinor
   integer :: shift_itypat,shift_itypat_nlmn,shift_itypat_3nlmn

   real(dp), pointer :: forces_(:,:)
   complex(dp), pointer :: cprj_(:,:),cprj_deriv_(:,:)
   real(dp), pointer :: cprj_real(:,:),cprj_deriv_real(:,:)
   real(dp) :: forces_tmp

   nspinor = xg_nonlop%nspinor

   ncols_cprj = cols(cprj)
   ncols_cprj_nospin = ncols_cprj/nspinor

   call xgBlock_reverseMap(forces,forces_)

   select case(xg_nonlop%space_cprj)

     case (SPACE_C)

       call xgBlock_reverseMap(cprj,cprj_)
       call xgBlock_reverseMap(cprj_deriv,cprj_deriv_)

       shift_itypat=0
       shift_itypat_nlmn=0
       shift_itypat_3nlmn=0
       !$omp parallel default (none) &
       !$omp& shared(xg_nonlop,forces_,cprj_deriv_,cprj_), &
       !$omp& firstprivate(shift_itypat,shift_itypat_nlmn,shift_itypat_3nlmn,ncols_cprj_nospin,nspinor), &
       !$omp& private(itypat,nattyp,nlmn,ia,ilmn,idir,iforces,my_iband), &
       !$omp& private(iband_spinor,icprj,icprj_deriv,forces_tmp)
       do itypat = 1, xg_nonlop%ntypat
         nlmn = xg_nonlop%nlmn_ntypat(itypat)
         nattyp = xg_nonlop%nattyp(itypat)
         !$omp do collapse(3)
         do iband=1,ncols_cprj_nospin
           do ia = 1, nattyp
             do idir=1,3
               forces_tmp = zero
               do ispinor=1,nspinor
                 do ilmn=1,nlmn
                   iband_spinor = ispinor + nspinor*(iband-1)
                   icprj       = ilmn + nlmn*(ia-1) + shift_itypat_nlmn
                   icprj_deriv = ilmn + nlmn*(idir-1) + 3*nlmn*(ia-1) + shift_itypat_3nlmn
                   forces_tmp = forces_tmp &
                     & + 2 * dble(conjg(cprj_deriv_(icprj_deriv,iband_spinor))*cprj_(icprj,iband_spinor))
                 end do
               end do
               iforces  = idir + 3*(ia-1) + shift_itypat
               my_iband = iband + xg_nonlop%me_band*ncols_cprj_nospin
               forces_(iforces,my_iband) = forces_(iforces,my_iband) + forces_tmp
             end do
           end do
         end do
         !$omp end do
         shift_itypat       = shift_itypat       + 3*nattyp
         shift_itypat_nlmn  = shift_itypat_nlmn  + nattyp*nlmn
         shift_itypat_3nlmn = shift_itypat_3nlmn + nattyp*3*nlmn
       end do
       !$omp end parallel

     case (SPACE_R)

       call xgBlock_reverseMap(cprj,cprj_real)
       call xgBlock_reverseMap(cprj_deriv,cprj_deriv_real)

       shift_itypat=0
       shift_itypat_nlmn=0
       shift_itypat_3nlmn=0
       !$omp parallel default (none) &
       !$omp& shared(xg_nonlop,forces_,cprj_deriv_real,cprj_real), &
       !$omp& firstprivate(shift_itypat,shift_itypat_nlmn,shift_itypat_3nlmn,ncols_cprj_nospin,nspinor), &
       !$omp& private(itypat,nattyp,nlmn,ia,ilmn,idir,iforces,my_iband), &
       !$omp& private(iband_spinor,icprj,icprj_deriv,forces_tmp)
       do itypat = 1, xg_nonlop%ntypat
         nlmn = xg_nonlop%nlmn_ntypat(itypat)
         nattyp = xg_nonlop%nattyp(itypat)
         !$omp do collapse(3)
         do iband=1,ncols_cprj_nospin
           do ia = 1, nattyp
             do idir=1,3
               forces_tmp = zero
               do ispinor=1,nspinor
                 do ilmn=1,nlmn
                   iband_spinor = ispinor + nspinor*(iband-1)
                   icprj       = ilmn + nlmn*(ia-1) + shift_itypat_nlmn
                   icprj_deriv = ilmn + nlmn*(idir-1) + 3*nlmn*(ia-1) + shift_itypat_3nlmn
                   forces_tmp = forces_tmp &
                     & + 2 * cprj_deriv_real(icprj_deriv,iband_spinor)*cprj_real(icprj,iband_spinor)
                 end do
               end do
               iforces  = idir + 3*(ia-1) + shift_itypat
               my_iband = iband + xg_nonlop%me_band*ncols_cprj_nospin
               forces_(iforces,my_iband) = forces_(iforces,my_iband) + forces_tmp
             end do
           end do
         end do
         !$omp end do
         shift_itypat       = shift_itypat       + 3*nattyp
         shift_itypat_nlmn  = shift_itypat_nlmn  + nattyp*nlmn
         shift_itypat_3nlmn = shift_itypat_3nlmn + nattyp*3*nlmn
       end do
       !$omp end parallel

     case default
       ABI_ERROR("Wrong space")

   end select

   call xgBlock_mpi_sum(forces,comm=xg_nonlop%comm_band)

 end subroutine xg_nonlop_mult_cprj_forces

!!****f* m_xg_nonlop/xg_nonlop_stress
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
subroutine xg_nonlop_mult_cprj_stress(xg_nonlop,cprj,cprj_deriv,stress)

   type(xg_nonlop_t), intent(in) :: xg_nonlop
   type(xgBlock_t), intent(in) :: cprj
   type(xgBlock_t), intent(in) :: cprj_deriv
   type(xgBlock_t), intent(inout) :: stress

   !real(dp) :: tsec(2)
   integer :: ia,idir,ilmn,iband,iband_spinor,my_iband,itypat,nlmn,nattyp
   integer :: ispinor,icprj,icprj_deriv
   integer :: ncols_cprj,ncols_cprj_nospin,nspinor
   integer :: shift_itypat,shift_itypat_nlmn,shift_itypat_6nlmn

   complex(dp), pointer :: cprj_(:,:),cprj_deriv_(:,:)
   real(dp), pointer :: cprj_real(:,:),cprj_deriv_real(:,:)
   real(dp), pointer :: stress_(:,:)
   real(dp) :: stress_tmp

!   call timab(tim_getHmeSX,1,tsec)

   nspinor = xg_nonlop%nspinor

   ncols_cprj = cols(cprj)
   ncols_cprj_nospin = ncols_cprj/nspinor

   call xgBlock_reverseMap(stress,stress_)

   select case(xg_nonlop%space_cprj)

     case (SPACE_C)

       call xgBlock_reverseMap(cprj,cprj_)
       call xgBlock_reverseMap(cprj_deriv,cprj_deriv_)

       shift_itypat=0
       shift_itypat_nlmn=0
       shift_itypat_6nlmn=0
       !$omp parallel default (none) &
       !$omp& shared(xg_nonlop,stress_,cprj_deriv_,cprj_), &
       !$omp& firstprivate(shift_itypat,shift_itypat_nlmn,shift_itypat_6nlmn,ncols_cprj_nospin,nspinor), &
       !$omp& private(itypat,nattyp,nlmn,ia,ilmn,idir,my_iband), &
       !$omp& private(iband_spinor,icprj,icprj_deriv,stress_tmp)
       do itypat = 1, xg_nonlop%ntypat
         nlmn = xg_nonlop%nlmn_ntypat(itypat)
         nattyp = xg_nonlop%nattyp(itypat)
         !$omp do collapse(2)
         do iband=1,ncols_cprj_nospin
           do idir=1,6
             stress_tmp = zero
             do ispinor=1,nspinor
               do ia = 1, nattyp
                 do ilmn=1,nlmn
                   iband_spinor = ispinor + nspinor*(iband-1)
                   icprj       = ilmn + nlmn*(ia-1) + shift_itypat_nlmn
                   icprj_deriv = ilmn + nlmn*(idir-1) + 6*nlmn*(ia-1) + shift_itypat_6nlmn
                   stress_tmp = stress_tmp &
                     & + 2 * dble(conjg(cprj_deriv_(icprj_deriv,iband_spinor))*cprj_(icprj,iband_spinor))
                 end do
               end do
             end do
             my_iband = iband + xg_nonlop%me_band*ncols_cprj_nospin
             stress_(idir,my_iband) = stress_(idir,my_iband) + stress_tmp
           end do
         end do
         !$omp end do
         shift_itypat       = shift_itypat       + 6*nattyp
         shift_itypat_nlmn  = shift_itypat_nlmn  + nattyp*nlmn
         shift_itypat_6nlmn = shift_itypat_6nlmn + nattyp*6*nlmn
       end do
       !$omp end parallel

     case (SPACE_R)

       call xgBlock_reverseMap(cprj,cprj_real)
       call xgBlock_reverseMap(cprj_deriv,cprj_deriv_real)

       shift_itypat=0
       shift_itypat_nlmn=0
       shift_itypat_6nlmn=0
       !$omp parallel default (none) &
       !$omp& shared(xg_nonlop,stress_,cprj_deriv_real,cprj_real), &
       !$omp& firstprivate(shift_itypat,shift_itypat_nlmn,shift_itypat_6nlmn,ncols_cprj_nospin,nspinor), &
       !$omp& private(itypat,nattyp,nlmn,ia,ilmn,idir,my_iband), &
       !$omp& private(iband_spinor,icprj,icprj_deriv,stress_tmp)
       do itypat = 1, xg_nonlop%ntypat
         nlmn = xg_nonlop%nlmn_ntypat(itypat)
         nattyp = xg_nonlop%nattyp(itypat)
         !$omp do collapse(2)
         do iband=1,ncols_cprj_nospin
           do idir=1,6
             stress_tmp = zero
             do ispinor=1,nspinor
               do ia = 1, nattyp
                 do ilmn=1,nlmn
                   iband_spinor = ispinor + nspinor*(iband-1)
                   icprj       = ilmn + nlmn*(ia-1) + shift_itypat_nlmn
                   icprj_deriv = ilmn + nlmn*(idir-1) + 6*nlmn*(ia-1) + shift_itypat_6nlmn
                   stress_tmp = stress_tmp &
                     & + 2 * cprj_deriv_real(icprj_deriv,iband_spinor)*cprj_real(icprj,iband_spinor)
                 end do
               end do
             end do
             my_iband = iband + xg_nonlop%me_band*ncols_cprj_nospin
             stress_(idir,my_iband) = stress_(idir,my_iband) + stress_tmp
           end do
         end do
         !$omp end do
         shift_itypat       = shift_itypat       + 6*nattyp
         shift_itypat_nlmn  = shift_itypat_nlmn  + nattyp*nlmn
         shift_itypat_6nlmn = shift_itypat_6nlmn + nattyp*6*nlmn
       end do
       !$omp end parallel

     case default
       ABI_ERROR("Wrong space")

   end select

   call xgBlock_mpi_sum(stress,comm=xg_nonlop%comm_band)


!   call timab(tim_getHmeSX,2,tsec)

 end subroutine xg_nonlop_mult_cprj_stress
!!***

end module m_xg_nonlop
!!***
