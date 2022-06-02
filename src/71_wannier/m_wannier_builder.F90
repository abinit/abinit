!!****m*ABINIT/m_wannier_builder
!! NAME
!!  m_wannier_builder
!!
!! FUNCTION
!!  Algorithms for building Wannier functions
!!  Methods:
!!  SCDM (select columns of density matrix method) and
!!  projected wannier function (PWF)
!! COPYRIGHT
!!  Copyright (C) 2005-2022 ABINIT group (hexu)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
!!
!! Todolist:
!! - output of Wannier90 format Amn, Mmn, eigenvalue files
!! - output header in netcdf file
!! - allow units in Hamiltonian
!! - projected wannier functions

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"


!===============================================================
! SCDM-k
!> @description: Select Column Density matrix method for
!   generating Wannier functions.
!===============================================================
module m_wannier_builder
  use defs_basis
  use m_abicore
  use m_errors
  use m_io_tools,        only : open_file
  use m_fstrings,        only : ltoa
  use m_scdm_math, only: complex_QRCP_piv_only, complex_svd, tpi_im, &
       & gaussian, fermi, insertion_sort_double, eigensolver
  use m_wann_netcdf, only: IOWannNC
  implicit none
  public:: WannierBuilder_t
  public:: Amn_to_H
  private


  !===============================================================
  ! scdmk type:
  !> @ description: the class for scdmk method.
  !===============================================================
  type::  WannierBuilder_t

     ! Inputs 
     real(dp), allocatable:: kpts(:, :) !(idim, ikpt)
     real(dp), allocatable:: kweights(:) !(ikpt)
     !real(dp), allocatable:: weight(:, :) !(iband, ikpt)
     integer, allocatable:: exclude_bands(:)

     integer:: method = 0  ! 1: SCDM-k 2: projected wf
     ! Disentanglement
     integer:: disentangle_func_type  ! 1:unity  2: erfc function 3: gauss function
     real(dp):: mu, sigma

     ! Projected Wannier function related
     complex(dp), allocatable:: projectors(:, :)

     ! SCDM related
     integer, allocatable:: cols(:)
     real(dp), allocatable:: anchor_kpt(:)
     integer:: anchor_ikpt
     integer, allocatable:: anchor_ibands(:)
     logical:: project_to_anchor = .False.

     ! Output Wannier function and Hamiltonian
     integer:: nwann, nkpt, nband, nbasis, nkdim
     integer:: dim  ! dimension of position
     integer:: nR
     integer, allocatable:: Rlist(:,:) !(idim, iRpt)
     complex(dp), allocatable:: Amnk(:, :, :) !(nband, nwann, nkpt)
     complex(dp), allocatable:: Hwannk(:, :, :) !(nwann, nwann, nkpt)

     complex(dp), allocatable:: psi_wann_k(:, :, :)  !(nbasis, nwann, nkpt)

     complex(dp), allocatable:: HwannR(:, :, :) !(nwann, nwann, nR)
     complex(dp), allocatable:: wannR(:,:,:) !(nbasis, nwann, nR)

   contains
     procedure:: initialize
     procedure:: finalize
     procedure:: get_psi_k
     procedure:: get_evals_k
     procedure:: run_all
     procedure:: find_kpoint
     !procedure:: remove_phase
     procedure:: auto_find_anchors
     procedure:: get_columns
     procedure:: get_scdm_Amnk  ! Amnk for all kpoints
     procedure:: get_projected_Amnk  ! Amnk for all kpoints
     procedure:: get_Amnk
     procedure:: get_weight
     procedure:: set_anchor
     procedure:: set_disp_projector
     procedure:: set_mode_projector
     procedure:: get_wannR_and_HwannR
     procedure:: select_columns
     procedure:: construct_wannier
     procedure:: write_Amnk_w90
     procedure:: write_Hwann_w90
     procedure:: create_ncfile
     procedure:: close_ncfile
     procedure:: write_wann_netcdf
     procedure:: get_wannier_eigen
     procedure:: get_wannier_eigen_klist
  end type WannierBuilder_t


  ! Extends WannierBuilder_t with pointer to eigens. 
  type, extends(WannierBuilder_t):: WannierBuilder_witheigen_t
     real(dp),  pointer:: evals(:, :) => null()   !(iband, ikpt)
     complex(dp),  pointer:: psi(:, :, :) => null() ! (ibasis, iband, ikpt)
   contains
     procedure :: set_eigen => set_eigen
     procedure :: get_psi_k => get_psi_k_from_eigen
     procedure :: get_evals_k => get_evals_k_from_eigen
  end type WannierBuilder_witheigen_t

contains

  !===============================================================
  !
  !> @
  !> kpts: kpoints. indices: (idim, ikpt)
  !> kweights: kweights. indices: (ikpt)
  !> nwann: number of Wannier functions to be calcualted.
  !> nbasis: number of basis in the original wavefunction.
  !> disentangle_func_type: the type of the disentanglement function. 1: unity function. 2. Fermi. 3. Gauss
  !> project_to_anchor: whether to multiply the weight function by the projection to the anchor states. 
  !===============================================================
  subroutine initialize(self,kpts, kweights, Rlist, nwann, nbasis, nband &
       &  disentangle_func_type, mu, sigma, exclude_bands, project_to_anchor, method)
    class(WannierBuilder_t), intent(inout):: self
    integer, intent(in):: nwann, nbasis, nband
    real(dp), intent(in):: kpts(:, :)  !(idim, ikpt)
    real(dp), intent(in):: kweights(:)  !(ikpt)
    integer, intent(in):: Rlist(:, :)  !(idim, iRpt)
    integer, optional, intent(in):: disentangle_func_type
    real(dp), optional, intent(in):: mu, sigma
    integer, optional, intent(in):: exclude_bands(:)
    logical, optional, intent(in):: project_to_anchor
    integer, intent(in):: method

    self%nkdim = size(kpts, 1)
    self%nkpt = size(kpts, 2)
    self%nwann = nwann
    self%nbasis = nbasis
    self%nband = nband
    self%method = method
    !if (present(psi_phase)) then
    !    if (psi_phase) then
    !        ABI_MALLOC(self%psi, (self%nbasis, self%nband, self%nkpt))
    !        call self%remove_phase(psi)
    !    else
    !        self%psi => psi
    !    end if
    !else
    !    self%psi => psi
    !end if

    ABI_MALLOC(self%cols, (self%nwann))
    self%cols(:) = 0
    ABI_MALLOC(self%kpts, (self%nkdim, self%nkpt))
    self%kpts = kpts
    ABI_MALLOC(self%kweights, (self%nkpt))
    self%kweights = kweights

    self%nR = size(Rlist, 2)
    ABI_MALLOC(self%Rlist, (self%nkdim, self%nR))
    self%Rlist = Rlist

    !ABI_MALLOC(self%weight, (self%nband, self%nkpt))
    !self%weight(:, :) = 0.0_dp

    if (present(disentangle_func_type)) then
       self%disentangle_func_type = disentangle_func_type
    else
       self%disentangle_func_type = 0
    end if

    if (present(mu)) then
       self%mu = mu
    else
       self%mu = 0
    end if

    if (present(sigma)) then
       self%sigma = sigma
    else
       self%sigma = sigma
    end if

    ABI_MALLOC(self%Amnk, (self%nband, self%nwann, self%nkpt))
    ABI_MALLOC(self%psi_wann_k, (self%nbasis, self%nwann, self%nkpt))
    ABI_MALLOC(self%Hwannk, (self%nwann, self%nwann, self%nkpt))

    if (present(exclude_bands)) then
       ABI_MALLOC(self%exclude_bands, (size(exclude_bands, 1)))
    end if

    if (present(project_to_anchor)) self%project_to_anchor = project_to_anchor
  end subroutine initialize



  subroutine finalize(self)
    class(WannierBuilder_t), intent(inout):: self
    ABI_SFREE(self%cols)
    ABI_SFREE(self%kpts)
    ABI_SFREE(self%kweights)
    ABI_SFREE(self%Amnk)
    ABI_SFREE(self%psi_wann_k)
    ABI_SFREE(self%Hwannk)
    ABI_SFREE(self%Rlist)
    ABI_SFREE(self%wannR)
    ABI_SFREE(self%HwannR)
    ABI_SFREE(self%exclude_bands)
    select case (self%method)
    case(1)
       ABI_SFREE(self%anchor_kpt)
       ABI_SFREE(self%anchor_ibands)
       ABI_SFREE(self%projectors)
    case(2)
       ABI_SFREE(self%projectors)
    end select
  end subroutine finalize

  function get_psi_k(self, ikpt) result(psik)
    class(WannierBuilder_t), intent(inout):: self
    integer, intent(in):: ikpt
    complex(dp),  pointer:: psik(:, :)
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(ikpt)
    ABI_UNUSED_A(psik)
    ABI_ERROR("WannierBuilder_t%get_psi_k should be overrided!")
  end function get_psi_k


  function get_evals_k(self, ikpt) result(ek)
    class(WannierBuilder_t), intent(inout):: self
    integer, intent(in):: ikpt
    real(dp),  pointer:: ek(:)
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(ikpt)
    ABI_UNUSED_A(ek)
    ABI_ERROR("WannierBuilder_t%get_evals_k should be overrided")
  end function get_evals_k


  ! automatically set the anchor points using the weight functions.
  ! The bands with the largest weights are selected as the anchor points. 
  subroutine auto_find_anchors(self, ianchors)
    class(WannierBuilder_t), intent(inout):: self
    integer, intent(out):: ianchors(self%nwann)
    integer:: i
    real(dp):: weights(self%nband)
    integer:: order(self%nband)

    call self%get_weight(self%anchor_ikpt, self%disentangle_func_type, &
         & self%mu, self%sigma, weights, project_to_anchor=.False.)
    call insertion_sort_double(weights, order)

    do i = 1, self%nwann
       ianchors(i)= order(self%nband-i+1)
    end  do
  end subroutine auto_find_anchors

  subroutine set_anchor(self, anchor_kpt, anchor_ibands)
    !> anchor_kpt: anchor kpoint (optional).
    !> anchor_ibands: the indices of mode used as anchor points (optional).
    class(WannierBuilder_t), intent(inout):: self
    real(dp), intent(in) ::  anchor_kpt(:)
    integer, optional, intent(in):: anchor_ibands(:)
    complex(dp), pointer :: psik(:,:)
    character(len = 500):: msg
    integer :: i
    ABI_MALLOC(self%anchor_ibands, (self%nwann))
    ABI_MALLOC(self%anchor_kpt, (size(anchor_kpt)))
    ABI_MALLOC(self%projectors, (self%nbasis, self%nwann))
    self%anchor_kpt = anchor_kpt
    self%anchor_ikpt = self%find_kpoint(anchor_kpt)

    if (.not. present(anchor_ibands)) then
       call wrtout( std_out, "Anchor points not specified, finding atomatically")
       call self%auto_find_anchors( self%anchor_ibands)
    else if (.not. size(anchor_ibands) == self%nwann) then
       ABI_ERROR("The number of anchor points should be equal to the number of Wannier functions.")
    else
       self%anchor_ibands = anchor_ibands
    end if
    write(msg, "(2a)") "Anchor point band indices set to ", trim(ltoa(self%anchor_ibands))
    call wrtout([ab_out, std_out], msg )
    psik=> self%get_psi_k(self%anchor_ikpt)
    do i = 1, self%nwann
        self%projectors(:, i)=psik(:, self%anchor_ibands(i))
  end do
  end subroutine set_anchor


  subroutine set_disp_projector(self, id_projectors)
    class(WannierBuilder_t), intent(inout):: self
    integer, intent(in):: id_projectors(:)
    integer:: i
    ABI_MALLOC(self%projectors, (self%nbasis, self%nwann))
    self%projectors(:,:)=zero
    do i = 1, self%nwann
       self%projectors(id_projectors(i), i)= cmplx(1.0, 0.0, dp)
    end do
  end subroutine  set_disp_projector

  subroutine set_mode_projector(self, kpoint, iband)
    class(WannierBuilder_t), intent(inout):: self
    real(dp), intent(in):: kpoint(:)
    integer:: iband(:)
    ABI_UNUSED_A(self)
    ABI_UNUSED(kpoint)
    ABI_UNUSED(iband)
  end subroutine  set_mode_projector


  subroutine select_columns(self)
    class(WannierBuilder_t), intent(inout):: self
    integer::  iband
    complex(dp):: psi_dagger(self%nband, self%nbasis)
    complex(dp):: psi_dagger_copy(self%nband, self%nbasis)
    real(dp):: weight(self%nband)
    !real(dp), pointer:: evals_anchor(:)
    !type(eigensolver):: esolver
    character(len = 500):: msg
    ! find anchor points, by default gamma
    self%anchor_ikpt = self%find_kpoint(self%anchor_kpt)
    ! TODO: add qpt if anchor_ikpt is not found
    !
    !if (size(self%anchor_ibands) /= 0) then
    !
    !end if
    ! calculate weight matrix for each kpoint
    call self%get_weight(self%anchor_ikpt, self%disentangle_func_type, self%mu, self%sigma, weight, &
            &project_to_anchor = .True.)
           ! &project_to_anchor = self%project_to_anchor)

    ! at anchor-kpoint, find cols
    ! psi: (ibasis, iband)
    psi_dagger = transpose(conjg(self%get_psi_k(self%anchor_ikpt)))
    do iband = 1, self%nband
       psi_dagger(iband, :) = psi_dagger(iband, :)*weight(iband)
    end do
    psi_dagger_copy(:,:) = psi_dagger(:,:)
    call self%get_columns(psi_dagger_copy, self%cols)
    write(msg, '(2a) ') 'Columns selected: ', trim(ltoa(self%cols))
    call wrtout([ab_out, std_out], msg )

    !psi_dagger_copy(:,:) = psi_dagger(:,:)
    !! check the eigen values:
    !call self%get_Amnk(self%anchor_ikpt, Amn)
    !! print the anchor point eigen values:
    !evals_anchor => self%get_evals_k(self%anchor_ikpt)
    !call Amn_to_H_from_evals(Amn, evals_anchor, &
    !    & self%nwann, self%nband, Hwann)
    !evals = evals_anchor(self%anchor_ibands)
    !write(msg, '(2a)') "The eigen values of the anchor points: ", &
    !     & trim(ltoa(evals))
    !call wrtout([ab_out, std_out], msg )
    !! calculate the eigen values of the Hwannk at anchor point
    !call esolver%run(evals, Hwann)
    !write(msg, '(2a)') "The eigen values of Hwann(k_anchor):   ", &
    !     & trim(ltoa(evals))
    !call wrtout([ab_out, std_out], msg )
    !call esolver%finalize()
  end subroutine select_columns

  subroutine construct_wannier(self)
    class(WannierBuilder_t), intent(inout):: self
    integer:: ikpt
    complex(dp), pointer:: p(:, :)
    !real(dp):: weight(self%nband)

    !complex(dp):: tmp(self%nwann, self%nwann)
    !real(dp):: evals(self%nwann)
    !type(eigensolver):: esolver

    if(self%method == 1)then
       call self%select_columns()
    end if

    ! For each kpoint, calculate Amn matrix, wannier function, and Hk at k
    do ikpt = 1, self%nkpt
       !Amnk (nband, nwann, nkpt)
       p => self%get_psi_k(ikpt)
       call self%get_Amnk(ikpt, self%Amnk(:, :, ikpt))
       ! psik*Amnk
       self%psi_wann_k(:, :, ikpt) = matmul(p, self%Amnk(:, :, ikpt))
       call Amn_to_H_from_evals(self%Amnk(:, :, ikpt), self%get_evals_k(ikpt), &
             self%nwann, self%nband, self%Hwannk(:, :, ikpt))
         !tmp = self%Hwannk(:,:,ikpt)
         !call esolver%run(evals, tmp)
         !call esolver%finalize()
       !print *, "ev1:", self%get_evals_k(ikpt)
      ! print *, "ev2:",evals

    end do
    ! Fourier transform of wannier function to real space
    call self%get_wannR_and_HwannR(self%Rlist)
  end subroutine construct_wannier



  subroutine run_all(self, ncfilename, Amnkfilename)
    class(WannierBuilder_t), intent(inout):: self
    character(*), intent(in):: ncfilename
    character(*), intent(in):: Amnkfilename
    type(IOWannNC):: ncfile
    call self%construct_wannier()
    call self%create_ncfile(ncfilename, ncfile)
    call self%write_wann_netcdf( ncfile,   &
        &wannR_unit='dimensionless', HwannR_unit='eV')
    call self%close_ncfile(ncfile)
    call self%write_Amnk_w90(trim(Amnkfilename))
  end subroutine run_all

  !subroutine remove_phase(self, psip)
  !    class(WannierBuilder_t), intent(inout):: self
  !    complex(dp), intent(in):: psip(:, :, :) ! (ibasis, iband, ikpt)
  !    !complex(dp), intent(out):: psi(:,:,:) ! (ibasis, iband, ikpt)
  !    integer:: ikpt, ibasis
  !    complex(dp):: phase
  !    do ikpt = 1, self%nkpt
  !        do ibasis = 1, self%nbasis
  !            phase = exp(-tpi_im*dot_product(self%kpts(:, ikpt), self%positions_red(:, ibasis)))
  !            self%psi(ibasis, :, ikpt) = psip(ibasis, :, ikpt)*phase
  !        end do
  !    end do
  !end subroutine remove_phase

  !===============================================================
  ! Find one kpoint in a list of kpoints.
  !> @
  !===============================================================
  function find_kpoint(self, kpoint) result(ik)
    class(WannierBuilder_t), intent(inout):: self
    real(dp), intent(in):: kpoint(:)
    integer:: ik, nk
    integer:: i
    real(dp):: a(size(self%kpts, 2))
    nk = size(self%kpts, 2)
    ! should transfer back to 1st BZ?
    do i = 1, nk
       a(i) = sum((self%kpts(:, i) - kpoint)**2)
    end do

    ik = minloc(a, dim = 1)
    if (a(ik) > 0.001) then
       ABI_ERROR("Error in finding kpoint from kpoint list. ")
    end if
  end function find_kpoint

  !===============================================================
  ! Calculate the weight function for each mode described by iband and ikpt
  ! The
  !> @
  !===============================================================
  subroutine get_weight(self, ikpt, disentanglement, mu, sigma, weight, project_to_anchor)
    class(WannierBuilder_t), intent(inout):: self
    integer, intent(in):: ikpt, disentanglement
    real(dp), intent(in):: mu, sigma
    real(dp), intent(inout):: weight(self%nband)
    logical,  intent(in):: project_to_anchor

    integer:: iband, ianchor
    real(dp):: proj
    complex(dp):: p
    real(dp), pointer:: ek(:)
    complex(dp), pointer:: psik(:,:)

    ek => self%get_evals_k(ikpt)
    select case (disentanglement)
    case (1)
       weight(:) = 1.0
    case (2)
       do iband = 1, self%nband
          weight(iband) = fermi(ek(iband), mu, sigma)
       end do
    case (3)
       do iband = 1, self%nband
          weight(iband) = gaussian(ek(iband), mu, sigma)
       end do
    case default
       ABI_ERROR("The disentanglement function type can only be 1:unity, 2: fermi, 3: gaussian")
    end select

    ! weight_mk *=\sum_anchor < psi_anchor | psi mk>

    if( project_to_anchor) then
      if (size(self%anchor_ibands) /= 0) then
      psik=>self%get_psi_k(ikpt)
        do iband = 1, self%nband
             proj = 0.0_dp
             do ianchor = 1, size(self%anchor_ibands)
               p = dot_product(self%projectors(:, ianchor), psik(:, iband))
               proj = proj+real(conjg(p)*p)
             end do
            weight(iband) = weight(iband)*proj
          end do
       end if
   end if
  end subroutine get_weight

  subroutine get_columns(self, psi_dagger, cols)
    class(WannierBuilder_t), intent(inout):: self
    complex(dp), intent(in):: psi_dagger(:, :)
    integer:: piv(size(psi_dagger, 2))
    integer, intent(inout):: cols(self%nwann)
    call complex_QRCP_piv_only(psi_dagger, piv)
    cols = piv(:self%nwann)
  end subroutine get_columns

  subroutine get_Amnk(self, ikpt, Amnk)
    class(WannierBuilder_t), intent(inout):: self
    integer, intent(in):: ikpt
    complex(dp), intent(inout):: Amnk(self%nband, self%nwann)
    if(self%method == 1) then
       call self%get_scdm_Amnk(ikpt, Amnk)
    else if(self%method == 2) then
       call self%get_projected_Amnk(ikpt, Amnk)
    end if
  end subroutine get_Amnk

  subroutine get_scdm_Amnk(self, ikpt, Amnk)
    class(WannierBuilder_t), intent(inout):: self
    integer, intent(in):: ikpt
    complex(dp):: psi_dagger(self%nband, self%nbasis)
    complex(dp), intent(inout):: Amnk(self%nband, self%nwann)
    real(dp):: weight(self%nband)
    complex(dp):: U(self%nband, self%nwann), VT(self%nwann, self%nwann)
    real(dp):: S(self%nband)
    !real(dp):: weights(self%nband)
    integer:: iband
    complex(dp),  pointer:: p(:, :)
 
    call self%get_weight(ikpt, self%disentangle_func_type, self%mu, self%sigma, weight, &
         &project_to_anchor = self%project_to_anchor)

    p => self%get_psi_k(ikpt)
    do iband = 1, self%nband
       psi_dagger(iband, :) = conjg(p(:, iband)) *weight(iband)
    end do

    ! orthogonalize selected columns
    call complex_svd(psi_dagger(:, self%cols), U, S, VT, 'S')
    Amnk(:, :) = matmul(U, VT)
  end subroutine get_scdm_Amnk

  
  subroutine get_projected_Amnk(self, ikpt, Amnk)
    class(WannierBuilder_t), intent(inout):: self
    integer, intent(in):: ikpt
    complex(dp), intent(inout):: Amnk(self%nband, self%nwann)
    complex(dp):: U(self%nband, self%nwann), VT(self%nwann, self%nwann)
    real(dp):: S(self%nband), weights(self%nband)
    integer:: iband, iwann
    complex(dp), pointer:: psi(:, :)

    psi => self%get_psi_k(ikpt)
    call self%get_weight(ikpt, self%disentangle_func_type, &
         & self%mu, self%sigma, weights, project_to_anchor=.False.)

    ! <proj||psi> * weight
    do iband = 1, self%nband
       do iwann = 1, self%nwann
          !Amnk(iband, iwann)= dot_product(conjg(self%projectors(:, iwann)), psi(:, iband )) * weights(iband)
          !Amnk(iband, iwann)= dot_product(self%projectors(:, iwann), conjg(psi(:, iband ))) * weights(iband)
          Amnk(iband, iwann)= dot_product(self%projectors(:, iwann), psi(:, iband )) * weights(iband)
       end do
    end do
    call complex_svd(Amnk, U, S, VT, 'S')
    Amnk(:, :) = matmul(U, VT)
  end subroutine get_projected_Amnk


  subroutine get_wannR_and_HwannR(self, Rlist)
    class(WannierBuilder_t), intent(inout):: self
    integer, intent(in):: Rlist(:, :)
    !complex(dp), intent(out):: HR(self%nwann, self%nwann, size(Rlist, 2))
    !-- H(R)= \sum_k H(k) * exp(i2pi k.R)
    integer:: ik, iR, nR
    complex(dp):: factor
    nR = size(Rlist, 2)
    ABI_MALLOC(self%HwannR, (self%nwann, self%nwann, nR))
    ABI_MALLOC(self%WannR, (self%nbasis, self%nwann, nR))

    self%HwannR(:, :, :) = cmplx(0.0, 0.0, dp)
    self%wannR(:, :, :) = cmplx(0.0, 0.0, dp)
    do ik = 1, self%nkpt
       do iR = 1, nR
          factor = exp(-tpi_Im*dot_product(self%kpts(:, ik), Rlist(:, iR))) * self%kweights(ik)
          self%HwannR(:, :, iR) = self%HwannR(:,:, iR) + self%Hwannk(:, :, ik)*factor
          self%wannR(:, :, iR) = self%wannR(:,:, iR) + self%psi_wann_k(:, :, ik)*factor
       end do
    end do
  end subroutine get_wannR_and_HwannR

  subroutine Amn_to_H_from_evals(Amn, evals, nwann, nband, Hwann)
    integer, intent(in)::  nwann, nband
    complex(dp), intent(in):: Amn(nband, nwann)
    real(dp), intent(in):: evals(nband)
    complex(dp), intent(inout):: Hwann(nwann, nwann)
    integer:: i
    complex(dp):: tmp(nwann, nband)

    ! A\dagger E
    do i = 1, nband
       tmp(:, i) = conjg(Amn(i, :))*evals(i)
    end do
    ! Hwann = A\dagger @ E @ A
    Hwann = matmul(tmp, Amn)
  end subroutine Amn_to_H_from_evals

  subroutine Amn_to_H(Amn, psi, H0, nbasis, nwann, Hwann)
    ! Hwann = psi\dagger A
    complex(dp), intent(in):: Amn(:, :), psi(:, :), H0(:, :)
    complex(dp), intent(inout):: Hwann(:, :)
    integer, intent(in):: nbasis, nwann
    complex(dp):: tmp(nbasis, nwann)
    !Hwann = A_dagger@psi_dagger@H0@psi@A
    tmp(:, :) = matmul(psi, Amn)
    Hwann = matmul(matmul(transpose(conjg(tmp)), H0), tmp)
  end subroutine Amn_to_H

  subroutine write_Amnk_w90(self, fname)
    ! write to Amnk file
    class(WannierBuilder_t), intent(inout):: self
    character(len=*), intent(in):: fname
    integer:: iwann, iband, ikpt, locibnd
    integer:: iun_amn
    character(len=500):: msg

    if (open_file(trim(fname)//".amn", msg, newunit=iun_amn, &
         & form="formatted", status="unknown", action="write") /= 0) then
       ABI_ERROR(msg)
    end if


    !IF (wan_mode=='standalone') THEN
    !   iun_amn = find_free_unit()
    !   IF (ionode) OPEN (unit = iun_amn, file = trim(seedname)//".amn",form='formatted')
    !ENDIF

    !TODO: re-enable this.
    !WRITE(stdout, '(a, i8)') '  AMN: iknum = ',iknum
    !
    !IF (wan_mode=='standalone') THEN
    !   CALL date_and_tim( cdate, ctime )
    !   header='Created on '//cdate//' at '//ctime//' with SCDM '
    !   IF (ionode) THEN
    !      WRITE (iun_amn, *) header
    !      WRITE (iun_amn, '(3i8, xxx, 2f10.6)') numbands,  iknum, n_wannier, scdm_mu, scdm_sigma
    !   ENDIF
    !ENDIF

    do ikpt = 1, self%nkpt
       do iwann = 1, self%nwann
          locibnd = 0
          do iband = 1, self%nband
             !IF (excluded_band(iband)) CYCLE
             locibnd = locibnd+1
             WRITE (iun_amn, '(3i5, 2f18.12)') locibnd, iwann, ikpt, &
                  & REAL(self%Amnk(locibnd, iwann, ikpt)), &
                  & AIMAG(self%Amnk(locibnd, iwann, ikpt))
          end do
       end do
    end do
    close (iun_amn)
  end subroutine write_Amnk_w90

  subroutine write_Hwann_w90(self, HR, Rlist, fname)
    class(WannierBuilder_t), intent(inout):: self
    complex(dp), intent(in):: HR(:, :, :)
    integer, intent(in):: Rlist(:, :)
    character(len=*), intent(in):: fname
    integer:: iR, ifile, iwann1, iwann2
    character(len=500):: msg

    if (open_file(trim(fname)//".hr", msg, newunit=ifile, &
         & form="formatted", status="unknown", action="write") /= 0) then
       ABI_ERROR(msg)
    end if


    do iR = 1, size(Rlist, 2)
       WRITE (ifile, '(3i5)') Rlist(:, iR)
       do iwann1 = 1, self%nwann
          do iwann2 = 1, self%nwann
             WRITE (ifile, '(f18.12)') HR(iwann1, iwann2, iR)
          end do
       end do
    end do
    close (ifile)
  end subroutine write_Hwann_w90

  subroutine create_ncfile(self, fname, ncfile)
    class(WannierBuilder_t), intent(inout):: self
    type(IOWannNC):: ncfile
    character(len=*), intent(in):: fname
    ABI_UNUSED_A(self)
    call ncfile%initialize(filename = fname)
  end subroutine create_ncfile

  subroutine close_ncfile(self, ncfile)
    class(WannierBuilder_t), intent(inout):: self
    type(IOWannNC):: ncfile
    ABI_UNUSED_A(self)
    call ncfile%close_file()
  end subroutine close_ncfile


  subroutine write_wann_netcdf(self, ncfile, wannR_unit, HwannR_unit)
    class(WannierBuilder_t), intent(inout):: self
    character(*), intent(in):: HwannR_unit, wannR_unit
    type(IOWannNC), intent(inout):: ncfile
    call ncfile%write_wann( nR = self%nR, ndim = self%nkdim, &
         & nwann = self%nwann, nbasis = self%nbasis, Rlist = self%Rlist, &
         & wannR = self%wannR, HwannR = self%HwannR, &
         & wannR_unit = wannR_unit, HwannR_unit = HwannR_unit)
    call ncfile%write_Amnk(nkpt = self%nkpt, nband = self%nband, nwann = self%nwann, &
         & kpoints = self%kpts, eigvals = self%evals, Amnk = self%Amnk)
  end subroutine write_wann_netcdf

  subroutine get_wannier_eigen(self, kpoint, evals, evecs)
    class(WannierBuilder_t), intent(inout):: self
    real(dp), intent(in) :: kpoint(3)
    real(dp), intent(inout) :: evals(self%nwann)
    complex(dp), optional, intent(inout) :: evecs(self%nwann, self%nwann)
    complex(dp) :: Hk(self%nwann, self%nwann), phase
    type(eigensolver):: esolver
    integer :: iR
    Hk(:,:)=0.0_dp
    do iR=1, self%nR
       phase = exp(tpi_im * dot_product(kpoint, self%Rlist(:, iR)))
       Hk = Hk + self%HwannR(:, :, iR) * phase
    end do
    call esolver%run(evals, Hk)
    ! Hk is overwritten as evecs
    if (present(evecs)) then
       evecs=Hk
    end if
    call esolver%finalize()
  end subroutine get_wannier_eigen

  subroutine get_wannier_eigen_klist(self, kpoints, nk, evals_nk, evecs_nk)
    class(WannierBuilder_t), intent(inout):: self
    integer, intent(in) :: nk
    real(dp), intent(in) :: kpoints(3, nk)
    real(dp), intent(inout) :: evals_nk(self%nwann, nk)
    complex(dp), optional, intent(inout) :: evecs_nk(self%nwann, self%nwann, nk)
    integer :: ik
    do ik =1, nk
       if (present(evecs_nk)) then
          call self%get_wannier_eigen(kpoints(:, ik), &
               & evals_nk(:, ik), evecs_nk(:,:, ik))
       else
          call self%get_wannier_eigen(kpoints(:, ik), evals_nk(:, ik))
       end if
    end do
  end subroutine get_wannier_eigen_klist

end module m_wannier_builder

!============================   WannierBuilder_witheigen_t   =========================

subroutine set_eigen(self,  evals, psi)
  class(WannierBuilder_witheigen_t) :: self
  real(dp), intent(in), target:: evals(:, :)   !(iband, ikpt)
  complex(dp), intent(in), target:: psi(:, :, :) ! (ibasis, iband, ikpt)
  self%psi => psi
  self%evals => evals
end subroutine set_eigen

function get_evals_k_from_eigen(self, ikpt) result(ek)
  class(WannierBuilder_witheigen_t), intent(inout):: self
  integer, intent(in):: ikpt
  real(dp), pointer:: ek(:)
  ek => self%evals( :, ikpt)
end function get_evals_k_from_eigen


function get_psi_k_from_eigen(self, ikpt) result(psik)
  class(WannierBuilder_witheigen_t), intent(inout):: self
  integer, intent(in):: ikpt
  complex(dp),  pointer:: psik(:, :)
  psik => self%psi(:, :, ikpt)
end function get_psi_k_from_eigen


!!***
