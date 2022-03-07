#include "abi_common.h"

module m_wann_netcdf
  use defs_basis
  use m_nctk, only: ab_define_var, netcdf_check
  use netcdf
  implicit none
  private

  type, public :: IOWannNC
     ! id of netcdf file
     integer :: ncid
     ! d_label is the dimension id in netcdf file
     integer :: d_nR, d_ndim, d_nwann, d_nbasis
     integer :: d_nkpt, d_nband
     integer :: d_three, d_natom
     ! i_label is the variable id in netcdf file
     integer :: i_Rlist, i_wannR_real, i_wannR_imag, i_HwannR_real, i_HwannR_imag
     integer :: i_cell, i_numbers, i_masses, i_xred, i_xcart
     integer :: i_kpts, i_eigvals, i_Amnk_real, i_Amnk_imag
   contains
     procedure :: write_wann
     procedure :: write_atoms
     procedure :: close_file
     procedure :: write_Amnk
  end type IOWannNC

  contains

  subroutine write_wann( self,  filename, nR, ndim, nwann, nbasis, Rlist, WannR, HwannR)
    class(IOwannNC), intent(inout) :: self
    character(len=*),intent(in) :: filename
    integer, intent(in) :: nR, ndim, nwann, nbasis, Rlist(:,:)
    complex(dp) , intent(in) :: WannR(:,:,:), HwannR(:,:,:)
    ! id of variables
    integer :: ncerr

    ncerr = nf90_create(path=trim(filename), cmode=NF90_CLOBBER, ncid=self%ncid)
    NCF_CHECK_MSG(ncerr, "Error when creating wannier netcdf  file")

    ! define dimensions
    ncerr=nf90_def_dim(self%ncid, "ndim", ndim , self%d_ndim)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension ndim in wannier nc file.")

    ncerr=nf90_def_dim(self%ncid, "nwann", nwann,  self%d_nwann)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension nwann in wannier nc file.")

    ncerr=nf90_def_dim(self%ncid, "nbasis", nbasis,  self%d_nbasis)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension nbasis in wannier nc file.")


    ncerr=nf90_def_dim(self%ncid, "nR", nR,  self%d_nR)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension nR in wannier nc file.")

    ! define variables
    call ab_define_var(self%ncid, [self%d_ndim, self%d_nR], self%i_Rlist, &
         & NF90_INT, "Rlist", "List of R vectors", "dimensionless")

    call ab_define_var(self%ncid, [self%d_nbasis, self%d_nwann, self%d_nR], &
         & self%i_wannR_real, NF90_DOUBLE, "wannier_function_real", "The real part of the Wannier function", "Angstrom")

    call ab_define_var(self%ncid, [self%d_nbasis, self%d_nwann, self%d_nR], &
         & self%i_wannR_imag, NF90_DOUBLE, "wannier_function_imag", "The imaginary part of the Wannier function", "Angstrom")

    call ab_define_var(self%ncid, [self%d_nwann, self%d_nwann, self%d_nR], &
         & self%i_HwannR_real, NF90_DOUBLE, "HamR_real", "The real part of the Wannier Hamiltonian", "eV")

    call ab_define_var(self%ncid, [self%d_nwann, self%d_nwann, self%d_nR], &
         & self%i_HwannR_imag, NF90_DOUBLE, "HamR_imag", "The imaginary part of the Wannier Hamiltonian", "eV")
    ncerr =nf90_enddef(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error when ending def mode in wannier netcdf file")

    ! set variables
    
    ncerr = nf90_put_var(self%ncid, self%i_Rlist, Rlist, start=[1,1], count=[3, nR])
    NCF_CHECK_MSG(ncerr, "Error when writting Rlist in wannier netcdf file.")

    ncerr = nf90_put_var(self%ncid, self%i_wannR_real, real(real(wannR)), &
         & start=[1,1], count=[nbasis, nwann, nR])
    NCF_CHECK_MSG(ncerr, "Error when writting WannR_real in wannier netcdf file.")

    ncerr = nf90_put_var(self%ncid, self%i_wannR_imag, real(aimag(wannR)), &
         & start=[1,1], count=[nbasis, nwann, nR])
    NCF_CHECK_MSG(ncerr, "Error when writting WannR_imag in wannier netcdf file.")



    ncerr = nf90_put_var(self%ncid, self%i_HwannR_real, real(real(HwannR)), &
         & start=[1,1], count=[nwann, nwann, nR])
    NCF_CHECK_MSG(ncerr, "Error when writting HwannR_real in wannier netcdf file.")

    ncerr = nf90_put_var(self%ncid, self%i_HwannR_imag, real(aimag(HwannR)), &
         & start=[1,1], count=[nwann, nwann, nR])
    NCF_CHECK_MSG(ncerr, "Error when writting HwannR_imag in wannier netcdf file.")

  end subroutine write_wann

  subroutine close_file(self)
    class(IOwannNC), intent(inout) :: self
    integer :: ncerr
    ncerr = nf90_close(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error close wannier netcdf file.")
  end subroutine close_file

  subroutine write_Amnk(self, nkpt, nband, nwann, kpoints, eigvals, Amnk)
    class(IOwannNC), intent(inout) :: self
    real(dp), intent(in) :: kpoints(:, :)
    integer, intent(in) :: nkpt, nband, nwann
    real(dp),  intent(in) :: eigvals(:,:)
    complex(dp),  intent(in) :: Amnk(:,:, :)
    integer :: ncerr

    ncerr = nf90_redef(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error starting redef in wannier netcdf file")

    ncerr=nf90_def_dim(self%ncid, "nkpt", nkpt, self%d_nkpt)
    NCF_CHECK_MSG(ncerr, "Error defining nkpt in wannier netcdf file")

    ncerr=nf90_def_dim(self%ncid, "nband", nband, self%d_nband)
    NCF_CHECK_MSG(ncerr, "Error defining nband in wannier netcdf file")

    call ab_define_var(self%ncid, [self%d_ndim, self%d_nkpt], &
         & self%i_kpts, NF90_DOUBLE, "kpoints", &
         &"KPOINTS" , "dimensionless")

    call ab_define_var(self%ncid, [self%d_nband, self%d_nkpt], &
         & self%i_eigvals, NF90_DOUBLE, "eigvals", &
         &"EIGen VALueS" , "eV")

    call ab_define_var(self%ncid, [self%d_nband, self%d_nwann,  self%d_nkpt], &
         & self%i_Amnk_real, NF90_DOUBLE, "Amnk_real", &
         &"The Amnk matrix: the REAL part" , "unitless")

    call ab_define_var(self%ncid, [self%d_nband, self%d_nwann,  self%d_nkpt], &
         & self%i_Amnk_imag, NF90_DOUBLE, "Amnk_imag", &
         &"The Amnk matrix: the imaginary part" , "unitless")

    ncerr =nf90_enddef(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error ending redef in wannier netcdf file")

    ncerr = nf90_put_var(self%ncid, self%i_kpts, kpoints, start=[1,1], count=[3, nkpt])
    NCF_CHECK_MSG(ncerr, "Error when writting kpoints in wannier netcdf file.")

    ncerr = nf90_put_var(self%ncid, self%i_eigvals, eigvals, start=[1,1], count=[nband, nkpt])
    NCF_CHECK_MSG(ncerr, "Error when writting eigvals in wannier netcdf file.")

    ncerr = nf90_put_var(self%ncid, self%i_Amnk_real, real(real(Amnk)), &
         & start=[1,1], count=[nband, nwann, nkpt])
    NCF_CHECK_MSG(ncerr, "Error when writting Amnk_real in wannier netcdf file.")

    ncerr = nf90_put_var(self%ncid, self%i_Amnk_imag, real(aimag(Amnk)), &
         & start=[1,1], count=[nband, nwann, nkpt])
    NCF_CHECK_MSG(ncerr, "Error when writting Amnk_imag in wannier netcdf file.")

  end subroutine write_Amnk

  subroutine write_atoms(self, natom, cell, numbers, masses, xred, xcart)
    class(IOwannNC), intent(inout) :: self
    integer , intent(in) :: natom
    integer, intent(in) :: numbers(:)
    real(dp), intent(in):: cell(:,:), masses(:), xred(:, :), xcart(:,:)
    integer :: ncerr

    ncerr = nf90_redef(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error starting redef in wannier netcdf file")
    ncerr=nf90_def_dim(self%ncid, "three", 3, self%d_three)
    NCF_CHECK_MSG(ncerr, "Error defining three in wannier netcdf file")

    ncerr=nf90_def_dim(self%ncid, "natom", natom, self%d_natom)
    NCF_CHECK_MSG(ncerr, "Error defining natom in wannier netcdf file")

    !integer :: i_cell, i_numbers, i_masses, i_xred, i_xcart

    call ab_define_var(self%ncid, [self%d_three, self%d_three], &
         & self%i_cell, NF90_DOUBLE, "atomic_cell", &
         &"ATOMic CELL" , "Angstrom")

    call ab_define_var(self%ncid, [ self%d_natom], &
         & self%i_numbers, NF90_DOUBLE, "atomic_numbers", &
         &"ATOMic NUMBERS" , "atomic unit")

    call ab_define_var(self%ncid, [ self%d_natom], &
         & self%i_masses, NF90_DOUBLE, "atomic_masses", &
         &"ATOMic MASSES" , "atomic unit")

    call ab_define_var(self%ncid, [self%d_three, self%d_natom], &
         & self%i_xred, NF90_DOUBLE, "atomic_xred", &
         &"ATOMic structure: Reduced positions" , "dimensionless")

    call ab_define_var(self%ncid, [self%d_three, self%d_natom], &
         & self%i_xcart, NF90_DOUBLE, "atomic_xcart", &
         &"ATOMic structure: CARTesian positions" , "Angstrom")


    ncerr =nf90_enddef(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error ending redef in wannier netcdf file")

    ncerr = nf90_put_var(self%ncid, self%i_cell, cell, start=[1,1], count=[3, 3])
    NCF_CHECK_MSG(ncerr, "Error when writting atomic_cell in wannier netcdf file.")

    ncerr = nf90_put_var(self%ncid, self%i_numbers, numbers, start=[1], count=[natom])
    NCF_CHECK_MSG(ncerr, "Error when writting atomic_numbers in wannier netcdf file.")

    ncerr = nf90_put_var(self%ncid, self%i_masses, masses, start=[1], count=[natom])
    NCF_CHECK_MSG(ncerr, "Error when writting atomic_masses in wannier netcdf file.")

    ncerr = nf90_put_var(self%ncid, self%i_xred, xred, start=[1, 1], count=[3, natom])
    NCF_CHECK_MSG(ncerr, "Error when writting atomic_xred in wannier netcdf file.")

    ncerr = nf90_put_var(self%ncid, self%i_xcart, xcart, start=[1, 1], count=[3, natom])
    NCF_CHECK_MSG(ncerr, "Error when writting atomic_xred in wannier netcdf file.")

  end subroutine write_atoms

end module m_wann_netcdf
