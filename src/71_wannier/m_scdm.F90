!!****m* ABINIT/m_scdm
!! NAME
!!  m_scdm
!!
!! FUNCTION
!!  SCDM (select columns of density matrix method) for getting wannier functions.
!!
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"


!===============================================================
! SCDM-k
!> @description: Select Column Density matrix method for
!   generating Wannier functions.
!===============================================================
module m_scdm
    use defs_basis, only: dp, PI
    use m_errors
    use m_scdm_math, only: complex_QRCP_piv_only, complex_svd, tpi_im, &
        & gaussian, fermi, insertion_sort_double, eigensolver
    use m_wann_netcdf, only: IOWannNC
    implicit none
    public :: scdmk_t
    public :: Amn_to_H
    private


    !===============================================================
    ! scdmk type:
    !> @ description: the class for scdmk method.
    !===============================================================
    type::  scdmk_t
        real(dp), pointer :: evals(:, :) => null()   !(iband, ikpt)
        complex(dp), pointer :: psi(:, :, :) => null() ! (ibasis, iband, ikpt)
        real(dp), allocatable :: kpts(:, :) !(idim, ikpt)
        integer, allocatable :: Rlist(:,:) !(idim, iRpt)
        real(dp), allocatable :: kweights(:) !(ikpt)
        real(dp), allocatable :: weight(:, :) !(iband, ikpt)
        integer, allocatable :: cols(:)
        real(dp), allocatable :: anchor_kpt(:)
        integer :: anchor_ikpt
        integer, allocatable :: anchor_ibands(:)
        integer :: nwann, nkpt, nband, nbasis, nkdim
        integer :: dim ! dimension of position
        integer :: nR
        integer :: disentangle_func_type
        real(dp) :: mu, sigma
        complex(dp), allocatable :: Amnk(:, :, :) !(nband, nwann, nkpt)
        complex(dp), allocatable :: Hwannk(:, :, :) !(nwann, nwann, nkpt)

        complex(dp), allocatable :: psi_wann_k(:, :, :)  !(nbasis, nwann, nkpt)

        complex(dp), allocatable :: HwannR(:, :, :) !(nwann, nwann, nR)
        complex(dp), allocatable :: wannR(:,:,:) !(nbasis, nwann, nR)
        integer, allocatable :: exclude_bands(:)
        logical :: project_to_anchor = .False.
    contains
        procedure :: initialize
        procedure :: finalize
        procedure :: run_all
        procedure :: find_kpoint
        !procedure :: remove_phase
        procedure :: auto_find_anchors
        procedure :: get_columns
        procedure :: get_Amnk ! Amnk for all kpoints
        !procedure :: get_anchor_projections
        procedure :: get_weight
        procedure :: set_anchor
        !procedure :: select_column
        procedure :: get_wannR_and_HwannR
        procedure :: write_Amnk_w90
        procedure :: write_Hwann_w90
        procedure :: write_wann_netcdf
    end type scdmk_t

contains

    !===============================================================
    !
    !> @
    !> evals: pointer to eigen values. 2D real matrix. The indices are (iband, ikpt).
    !> psi: pointer to wavefunctions. complex matrix. The indices are (ibasis, iband, ikpt)
    !> kpts: kpoints. indices: (idim, ikpt)
    !> kweights: kweights. indices: (ikpt)
    !> nwann: number of Wannier functions to be calcualted.
    !> disentangle_func_type: the type of the disentanglement function. 1: unity function. 2. Fermi. 3. Gauss
    !> project_to_anchor: whether to multiply the weight function by the projection to the anchor states. 
    !===============================================================
    subroutine initialize(self, evals, psi, kpts, kweights, Rlist, nwann, &
        &  disentangle_func_type, mu, sigma, exclude_bands, project_to_anchor)
        class(scdmk_t), intent(inout) :: self
        integer, intent(in) :: nwann
        real(dp), intent(in), target :: evals(:, :)   !(iband, ikpt)
        complex(dp), intent(in), target :: psi(:, :, :) ! (ibasis, iband, ikpt)
        real(dp), intent(in) :: kpts(:, :)  !(idim, ikpt)
        real(dp), intent(in) :: kweights(:)  !(ikpt)
        integer, intent(in) :: Rlist(:, :)  !(idim, iRpt)
        integer, optional, intent(in) :: disentangle_func_type
        real(dp), optional, intent(in) :: mu, sigma
        integer, optional, intent(in) :: exclude_bands(:)
        logical, optional, intent(in) :: project_to_anchor
        integer :: nR

        self%nkdim = size(kpts, 1)
        self%nkpt = size(kpts, 2)
        self%nbasis = size(psi, 1)
        self%nband = size(psi, 2)
        self%nwann = nwann
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

        self%psi => psi
        self%evals => evals
        ABI_MALLOC(self%cols, (self%nwann))
        self%cols(:) = 0
        ABI_MALLOC(self%kpts, (self%nkdim, self%nkpt))
        self%kpts = kpts
        ABI_MALLOC(self%kweights, (self%nkpt))
        self%kweights = kweights

        self%nR=size(Rlist, 2)
            ABI_MALLOC(self%Rlist, (self%nkdim, self%nR))
            self%Rlist = Rlist

            ABI_MALLOC(self%weight, (self%nband, self%nkpt))
            self%weight(:, :) = 0.0_dp

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

            if (present(project_to_anchor)) self%project_to_anchor=project_to_anchor
       end subroutine initialize

        subroutine finalize(self)
            class(scdmk_t), intent(inout) :: self
            nullify(self%psi)
            nullify(self%evals)
            ABI_SFREE(self%cols)
            ABI_SFREE(self%kpts)
            ABI_SFREE(self%kweights)
            ABI_SFREE(self%anchor_kpt)
            ABI_SFREE(self%anchor_ibands)
            ABI_SFREE(self%Amnk)
            ABI_SFREE(self%psi_wann_k)
            ABI_SFREE(self%Hwannk)
            ABI_SFREE(self%wannR)
            ABI_SFREE(self%HwannR)
        end subroutine finalize

        ! automatically set the anchor points using the weight functions.
        ! The bands with the largest weights are selected as the anchor points. 
        subroutine auto_find_anchors(self, anchor_kpt, ianchors)
            class(scdmk_t), intent(inout) :: self
            real(dp), intent(in) ::  anchor_kpt(:)
            integer, intent(out) :: ianchors(self%nwann)
            integer :: i, ikpt
            integer :: anchor_ibands(self%nwann)
            real(dp) :: weights(self%nband)
            integer :: order(self%nband)



            ikpt = self%find_kpoint(anchor_kpt)
            call self%get_weight(ikpt, self%disentangle_func_type, &
                & self%mu, self%sigma, weights, project_to_anchor=.False.)
           call insertion_sort_double(weights, order)

           do i = 1, self%nwann
            ianchors(i)= order(self%nband-i+1)
           end  do
        end subroutine auto_find_anchors

        subroutine set_anchor(self, anchor_kpt, anchor_ibands)
            !> anchor_kpt: anchor kpoint (optional).
            !> anchor_ibands: the indices of mode used as anchor points (optional).
            class(scdmk_t), intent(inout) :: self
            real(dp), optional, intent(in) ::  anchor_kpt(:)
            integer, optional, intent(in) :: anchor_ibands(:)
            ABI_MALLOC(self%anchor_ibands, (self%nwann))

            if (.not. present(anchor_ibands)) then
               print *, "Anchor points not specified, finding atomatically"
                call self%auto_find_anchors(anchor_kpt, self%anchor_ibands)
            else if (.not. size(anchor_ibands) == self%nwann) then
                ABI_ERROR("The number of anchor points should be equal to the number of Wannier functions.")
            else
                ABI_MALLOC(self%anchor_kpt, (size(anchor_kpt)))
                self%anchor_kpt = anchor_kpt
                self%anchor_ikpt = self%find_kpoint(anchor_kpt)
                self%anchor_ibands = anchor_ibands
            end if
            print *, "Anchor points set to ", self%anchor_ibands
        end subroutine set_anchor

        subroutine run_all(self)
            class(scdmk_t), intent(inout) :: self
            integer :: ikpt, iband, iwann
            complex(dp) :: psi_dagger(self%nband, self%nbasis)
            complex(dp) :: tmp(self%nband, self%nwann)

            type(eigensolver) :: esolver
            real(dp) :: evals(self%nwann)
            complex(dp) :: evecs(self%nwann, self%nwann)
            ! find anchor points, by default gamma
            print *, "fining anchor k"
            self%anchor_ikpt = self%find_kpoint(self%anchor_kpt)
            print *, "Found", self%anchor_ikpt
            !if (size(self%anchor_ibands) /= 0) then
            !
            !end if
            ! calculate weight matrix for each kpoint
            print *, "weight function"
            do ikpt = 1, self%nkpt
                call self%get_weight(ikpt, self%disentangle_func_type, self%mu, self%sigma, self%weight(:, ikpt), &
                  &project_to_anchor=self%project_to_anchor)
            end do

            ! at anchor-kpoint, find cols
            ! psi: (ibasis, iband)

            print *, "mixing weight function"
            psi_dagger = transpose(conjg(self%psi(:, :, self%anchor_ikpt)))
            do iband = 1, self%nband
                psi_dagger(iband, :) = psi_dagger(iband, :)*self%weight(iband, self%anchor_ikpt)
            end do

            print *, "select columns"
            call self%get_columns(psi_dagger, self%cols)
            print *, "columns selected: ", self%cols


            ! For each kpoint, calculate Amn matrix, wannier function, and Hk at k
            do ikpt = 1, self%nkpt
                do iband = 1, self%nband
                   psi_dagger(iband, :) = conjg(self%psi(:, iband, ikpt))*self%weight(iband, ikpt)
                end do
                !Amnk (nband, nwann, nkpt)
                call self%get_Amnk(psi_dagger, self%cols, self%Amnk(:, :, ikpt))
                ! psik * Amnk
                self%psi_wann_k(:, :, ikpt) = matmul(self%psi(:, :, ikpt), self%Amnk(:, :, ikpt))

                print *, "Getting Amnk"
                call Amn_to_H_from_evals(self%Amnk(:, :, ikpt), self%evals(:, ikpt), &
                        & self%nbasis, self%nwann, self%nband, self%Hwannk(:, :, ikpt))

                ! solve the eigens for the Hwannk
                evecs=self%Hwannk(:,:, ikpt)
                call esolver%run(evals, evecs)

            end do
            ! Fourier transform of wannier function to real space

            print *, "Hwannk -> HwannR"
            call self%get_wannR_and_HwannR(self%Rlist)
            print *, "Writting netcdf"
            call self%write_wann_netcdf("wann.nc")
            print *, "Writting Amnk"
            call self%write_Amnk_w90("Amnk.dat")
            print *, "Amnk written.  "
        end subroutine run_all

        !subroutine remove_phase(self, psip)
        !    class(scdmk_t), intent(inout) :: self
        !    complex(dp), intent(in) :: psip(:, :, :) ! (ibasis, iband, ikpt)
        !    !complex(dp), intent(out) :: psi(:,:,:) ! (ibasis, iband, ikpt)
        !    integer :: ikpt, ibasis
        !    complex(dp) :: phase
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
            class(scdmk_t), intent(inout) :: self
            real(dp), intent(in) :: kpoint(:)
            integer :: ik
            integer :: i
            real(dp) :: a(size(self%kpts, 2))
            ! should transfer back to 1st BZ?
            print *, size(self%kpts, 2)
            print *, size(a)
            print *, kpoint
            do i = 1, size(self%kpts, 2)
                a(i) = sum((self%kpts(:, i) - kpoint)**2)
            end do
            print *, a(i)
            ik = minloc(a, dim=1)
            print *, ik
            print *, a(ik)
            if (a(ik) > 0.001) then
                ABI_ERROR("Error in finding kpoint from kpoint list. ")
            end if
        end function find_kpoint

        !===============================================================
        ! Calculate the weight function for each mode described by iband and ikpt
        ! The
        !> @
        !===============================================================
        subroutine get_weight(self, ikpt, type, mu, sigma, weight, project_to_anchor)
            class(scdmk_t), intent(inout) :: self
            integer, intent(in) :: ikpt, type
            real(dp), intent(in) :: mu, sigma
            real(dp), intent(inout) :: weight(self%nband)
            logical, optional, intent(in) :: project_to_anchor

            integer :: iband, ianchor
            select case (type)
            case (1)
                weight(:) = 1.0
            case (2)
                do iband = 1, self%nband
                    weight(iband) = fermi(self%evals(iband, ikpt), mu, sigma)
                end do
            case (3)
                do iband = 1, self%nband
                    weight(iband) = gaussian(self%evals(iband, ikpt), mu, sigma)
                end do
            case default
                ABI_ERROR("The disentanglement function type can only be 1:unity, 2: fermi, 3: gaussian")
            end select

            ! weight_mk *=\sum_anchor< psi_anchor | psi mk>

            if (present(project_to_anchor)) then
            if (size(self%anchor_ibands) /= 0 .and. project_to_anchor) then
                do iband = 1, self%nband
                    do ianchor = 1, size(self%anchor_ibands)
                        weight(iband) = weight(iband)*real(dot_product(  &
                             & conjg(self%psi(:, ianchor, self%anchor_ikpt)), &
                             & self%psi(:, iband, ikpt)))
                    end do
                end do
            end if
        end if
        end subroutine get_weight

        subroutine get_columns(self, psi_dagger, cols)
            class(scdmk_t), intent(inout) :: self
            complex(dp), intent(in) :: psi_dagger(:, :)
            integer :: piv(size(psi_dagger, 2))
            integer, intent(inout) :: cols(self%nwann)
            call complex_QRCP_piv_only(psi_dagger, piv)
            cols = piv(:self%nwann)
        end subroutine get_columns

        subroutine get_Amnk(self, psi_dagger, cols, Amnk)
            class(scdmk_t), intent(inout) :: self
            complex(dp), intent(in) :: psi_dagger(self%nband, self%nbasis)
            integer :: cols(:)
            complex(dp), intent(inout) :: Amnk(self%nband, self%nwann)

            complex(dp) :: U(self%nband, self%nwann), VT(self%nwann, self%nwann)
            real(dp) :: S(self%nband)
            ! orthogonalize selected columns
            call complex_svd(psi_dagger(:, cols), U, S, VT, 'S')
            Amnk(:, :) = matmul(U, VT)
          end subroutine get_Amnk

    subroutine get_wannR_and_HwannR(self, Rlist)
        class(scdmk_t), intent(inout) :: self
        integer, intent(in) :: Rlist(:, :)
        !complex(dp), intent(out) :: HR(self%nwann, self%nwann, size(Rlist, 2))
        !-- H(R)= \sum_k H(k) * exp(i2pi k.R)
        integer :: ik, iR, nR
        complex(dp) :: factor
        nR = size(Rlist, 2)
        ABI_MALLOC(self%HwannR, (self%nwann, self%nwann,nR))
        ABI_MALLOC(self%WannR, (self%nbasis, self%nwann, nR))

        self%HwannR(:, :, :) = cmplx(0.0, 0.0, dp)
        self%wannR(:, :, :) = cmplx(0.0, 0.0, dp)
        do ik = 1, self%nkpt
            do iR = 1, nR
               factor=exp(tpi_Im*dot_product(self%kpts(:, ik), Rlist(:, iR))) * self%kweights(ik)
               self%HwannR(:, :, iR) = self%HwannR(:,:, iR) + self%Hwannk(:, :, ik)*factor
               self%wannR(:, :, iR) = self%wannR(:,:, iR) + self%psi_wann_k(:, :, ik)*factor
            end do
        end do
      end subroutine get_wannR_and_HwannR

    subroutine Amn_to_H_from_evals(Amn, evals, nbasis, nwann, nband, Hwann)
        complex(dp), intent(in) :: Amn(nband, nwann)
        real(dp), intent(in) :: evals(nband)
        integer, intent(in) :: nbasis, nwann, nband
        complex(dp), intent(inout) :: Hwann(nwann, nwann)

        integer :: i
        complex(dp) :: tmp(nwann, nband)

        ! A\dagger E
        do i = 1, nband
            tmp(:, i) = conjg(Amn(i, :))*evals(i)
        end do
        ! Hwann=A\dagger @ E @ A
        Hwann = matmul(tmp, Amn)
    end subroutine Amn_to_H_from_evals

    subroutine Amn_to_H(Amn, psi, H0, nbasis, nwann, nband, Hwann)
        ! Hwann = psi\dagger A
        complex(dp), intent(in) :: Amn(:, :), psi(:, :), H0(:, :)
        complex(dp), intent(inout) :: Hwann(:, :)
        integer, intent(in) :: nbasis, nwann, nband
        complex(dp) :: tmp(nbasis, nwann)
        !Hwann = A_dagger@psi_dagger@H0@psi@A
        tmp(:, :) = matmul(psi, Amn)
        Hwann = matmul(matmul(transpose(conjg(tmp)), H0), tmp)
    end subroutine Amn_to_H

    subroutine write_Amnk_w90(self, fname)
        ! write to Amnk file
        class(scdmk_t), intent(inout) :: self
        character(len=*), intent(in) :: fname
        integer :: iwann, iband, ikpt, locibnd
        integer :: iun_amn

        iun_amn = 101
        OPEN (unit=iun_amn, file=trim(fname)//".amn", form='formatted')

        !IF (wan_mode=='standalone') THEN
        !   iun_amn = find_free_unit()
        !   IF (ionode) OPEN (unit=iun_amn, file=trim(seedname)//".amn",form='formatted')
        !ENDIF

        !WRITE(stdout,'(a,i8)') '  AMN: iknum = ',iknum
        !
        !IF (wan_mode=='standalone') THEN
        !   CALL date_and_tim( cdate, ctime )
        !   header='Created on '//cdate//' at '//ctime//' with SCDM '
        !   IF (ionode) THEN
        !      WRITE (iun_amn,*) header
        !      WRITE (iun_amn,'(3i8,xxx,2f10.6)') numbands,  iknum, n_wannier, scdm_mu, scdm_sigma
        !   ENDIF
        !ENDIF

        do ikpt = 1, self%nkpt
            do iwann = 1, self%nwann
                locibnd = 0
                do iband = 1, self%nband
                    !IF (excluded_band(iband)) CYCLE
                    locibnd = locibnd + 1
                    WRITE (iun_amn, '(3i5,2f18.12)') locibnd, iwann, ikpt, &
                         & REAL(self%Amnk(locibnd, iwann, ikpt)), &
                         & AIMAG(self%Amnk(locibnd, iwann, ikpt))
                end do
            end do
        end do

        close (iun_amn)
    end subroutine write_Amnk_w90

    subroutine write_Hwann_w90(self, HR, Rlist, fname)
        class(scdmk_t), intent(inout) :: self
        complex(dp), intent(in) :: HR(:, :, :)
        integer, intent(in) :: Rlist(:, :)
        character(len=*), intent(in) :: fname
        integer :: iR, ifile, iwann1, iwann2
        ifile = 103
        OPEN (unit=ifile, file=trim(fname)//".hr", form='formatted')
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

    subroutine write_wann_netcdf(self, fname)
      class(scdmk_t), intent(inout) :: self
      character(len=*), intent(in) :: fname
      type(IOWannNC) :: myfile

      call myfile%write_wann(filename=fname, nR=self%nR, ndim=self%nkdim, &
           & nwann=self%nwann, nbasis=self%nbasis, Rlist=self%Rlist, &
           & wannR=self%wannR, HwannR=self%HwannR)
      call myfile%write_Amnk(nkpt=self%nkpt, nband=self%nband, nwann=self%nwann, &
           & kpoints=self%kpts, eigvals=self%evals, Amnk=self%Amnk)
      call myfile%close_file()
    end subroutine write_wann_netcdf


end module m_scdm


!!***
