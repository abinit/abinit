!!****m*ABINIT/m_lwf
!! NAME
!! m_lwf
!!
!! FUNCTION
!! Module for the lattice Wannier function
!! Container type is defined, and destruction, print subroutines
!! as well as the central mkphdos
!!
!! COPYRIGHT
!! Copyright (C) 1999-2021 ABINIT group (HeXu)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE


!! Todolist:
!! - LWF output filename change
!! - LWF specific outputs: lwf_masses, original crystal structure
!! - Units of Hamiltonian need to be changed to eV-Angstrom units
!! - dipole-dipole
!! - MPI
!! - exclude_bands

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_lwf

 use defs_basis
 use m_errors
 use m_xmpi
 use m_abicore
 use m_htetra
 use m_numeric_tools
 use m_cgtools
 use m_crystal
 use m_nctk
 use iso_c_binding
 use m_atprj
 use m_sortph
 use m_ddb
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_supercell
 use m_dtset
 use m_krank

 use m_fstrings,        only : itoa, ftoa, sjoin, ltoa, ktoa, strcat, basename, replace
 use m_symtk,           only : matr3inv
 use m_time,            only : cwtime, cwtime_report
 use m_io_tools,        only : open_file
 use m_geometry,        only : mkrdim, symredcart, normv, phdispl_cart2red
 use m_dynmat,          only : gtdyn9, dfpt_phfrq, dfpt_prtph, pheigvec_normalize, massmult_and_breaksym, phdispl_from_eigvec
 use m_bz_mesh,         only : isamek, make_path, kpath_t, kpath_new
 use m_ifc,             only : ifc_type
 use m_anaddb_dataset,  only : anaddb_dataset_type
! use m_kpts,            only : get_kpt_full_bz, 
 use m_kpts,            only : kpts_ibz_from_kptrlatt, get_full_kgrid

 use m_special_funcs,   only : bose_einstein
 use m_sort,            only : sort_dp
 use m_symfind,         only : symanal
 use m_scdm_math,       only : build_Rgrid
 use m_wannier_builder,            only : WannierBuilder_witheigen_t
 use m_wann_netcdf,     only : IOWannNC


 implicit none

 type LatticeWannier
    type(WannierBuilder_witheigen_t):: scdm
    type(ifc_type), pointer:: ifc
    type(crystal_t), pointer:: crystal
    integer:: qptrlatt(3, 3)
    real(dp):: shiftq(3)
    integer:: nqibz
    real(dp), allocatable:: qibz(:, :), qweights(:)
    integer:: nR, natom
    integer, allocatable:: Rlist(:, :)
    real(dp), allocatable:: eigenvalues(:, :)         ! iband, iqpt
    complex(dp), allocatable:: eigenvectors(:,:, :)  ! ibasis, iband, iqpt
    integer:: comm, nprocs, my_rank
  contains
    procedure:: initialize
    procedure:: finalize
    procedure:: init_mpi
    procedure:: sanity_check
    procedure:: prepare_qpoints
    procedure:: prepare_Rlist
    procedure:: get_ifc_eigens
    procedure:: run_all 
    procedure:: write_lwf_nc
    procedure:: print_Rlist
    procedure :: write_bands
 end type LatticeWannier

 private

 public:: run_lattice_wannier

contains
  subroutine initialize(self, ifc, crystal, dtset, comm)
    class(LatticeWannier), intent(inout):: self
    integer, intent(in):: comm
    type(ifc_type), target, intent(in):: ifc
    type(crystal_t), target, intent(in):: crystal
    type(anaddb_dataset_type), intent(in):: dtset
    real(dp):: symcart(3, 3, crystal%nsym)
    integer:: isym
    ! TODO: add exclude_bands
    integer:: exclude_bands(0)
    
    character(len = 500+dtset%lwf_nwann*10):: msg
    real(dp):: mu, sigma

    self%ifc => ifc
    self%crystal => crystal

    ! mpi vars
    call self%init_mpi(comm)

    ! sanity check
    call self%sanity_check(dtset)

    do isym = 1, crystal%nsym
       call symredcart(crystal%rprimd, crystal%gprimd, symcart(:,:,isym), crystal%symrel(:,:,isym))
    end do

    ! prepare qpoints
    call self%prepare_qpoints(crystal, dtset)
    call self%prepare_Rlist()
    ! prepare eigen values and eigen vectors

    call self%get_ifc_eigens(ifc, crystal)
    
    ! set mu and sigma from cm^-1 to eigenvalue
    if(dtset%lwfflag == 1 .or. dtset%lwfflag == 2) then
      mu = freq_to_eigenval(dtset%lwf_mu/Ha_cmm1)
      sigma = freq_to_eigenval(dtset%lwf_sigma/Ha_cmm1)
    end if

    ! set up scdm 

    call self%scdm%initialize(  &
         & kpts = self%qibz, kweights = self%qweights, Rlist = self%Rlist, &
         & nwann = dtset%lwf_nwann, nbasis=self%natom*3, nband=self%natom*3, disentangle_func_type = dtset%lwf_disentangle, &
         & mu = mu, sigma = sigma, exclude_bands = exclude_bands, &
         & project_to_anchor = (dtset%lwf_anchor_proj > 0 ), method = dtset%lwfflag)
    call self%scdm%set_eigen(evals = self%eigenvalues,  psi = self%eigenvectors)

     if(dtset%lwfflag == 1) then
        write(msg, '(a)')  ' Constructing LWF with SCDM-k method.'
        call wrtout([ab_out, std_out], msg ) 
        write(msg, '(a, i0)')  ' Number of LWF: ', dtset%lwf_nwann
        call wrtout([ab_out, std_out], msg ) 
        write(msg, '(a, 3f8.5)')  ' Anchor Points q-point: ', &
             & dtset%lwf_anchor_qpt(1), dtset%lwf_anchor_qpt(2), dtset%lwf_anchor_qpt(3)
        !write(msg, '(2a) ') 'Anchor points band indices: ', trim(ltoa(self%scdm%anchor_ibands))
        call wrtout([ab_out, std_out], msg ) 

       if(dtset%lwf_anchor_iband(1) > 0) then
         call self%scdm%set_anchor( dtset%lwf_anchor_qpt, dtset%lwf_anchor_iband)
       else
         call self%scdm%set_anchor(anchor_kpt = dtset%lwf_anchor_qpt)
       end if
    ! output information
    else if(dtset%lwfflag == 2) then
    ! TODO: projected lattice wannier function
       write(msg, '(a)')  ' Constructing LWF with projected wannier function method.'
       call wrtout([ab_out, std_out], msg ) 
       call self%scdm%set_disp_projector(dtset%lwf_projector)
       write(msg, '(2a)')  ' The projectors: ', trim(ltoa(dtset%lwf_projector))
       call wrtout([ab_out, std_out], msg ) 
    end if
   
    end subroutine initialize

  subroutine finalize(self)
    class(LatticeWannier), intent(inout):: self
    call self%scdm%finalize()
    ABI_SFREE(self%qibz)
    ABI_SFREE(self%qweights)
    ABI_SFREE(self%eigenvalues)
    ABI_SFREE(self%eigenvectors)
    ABI_SFREE(self%Rlist)
  end subroutine finalize

  subroutine init_mpi(self, comm)
    class(LatticeWannier), intent(inout):: self
    integer, intent(in):: comm
    self%comm = comm
    self%nprocs = xmpi_comm_size(comm)
    self%my_rank = xmpi_comm_rank(comm)
  end subroutine init_mpi

  subroutine sanity_check(self, dtset)
    class(LatticeWannier), intent(inout):: self
    type(anaddb_dataset_type), intent(in):: dtset
    if(self%nprocs /= 1) then
       ABI_ERROR(" MPI is not yet implemented for Lattice Wannier function.")
    end if

    if(dtset%dipdip > 0) then
       ABI_ERROR(" dipdip is not yet implemented for Lattice Wannier function.")
    end if
  end subroutine sanity_check


  subroutine prepare_qpoints(self, crystal, dtset)
    class(LatticeWannier), intent(inout):: self
    type(crystal_t), intent(in):: crystal
    type(anaddb_dataset_type), intent(in):: dtset
    integer:: nqshft = 1
    real(dp):: lwf_qshift(3, 1)
    integer:: nqbz
    real(dp), allocatable:: qbz(:, :)
    integer:: in_qptrlatt(3, 3), new_qptrlatt(3, 3)
    integer, allocatable:: bz2ibz_smap(:,:)!, bz2ibz(:)
    real(dp), allocatable::  new_shiftq(:,:)
    integer, parameter:: bcorr0 = 0, master = 0
    integer:: my_qptopt
    character(len = 500):: msg
    integer:: iqpt, nkpout
    ! Copied from m_phonons/mkphdos
    in_qptrlatt = 0
    in_qptrlatt(1, 1) = dtset%lwf_ngqpt(1)
    in_qptrlatt(2, 2) = dtset%lwf_ngqpt(2)
    in_qptrlatt(3, 3) = dtset%lwf_ngqpt(3)

    ! TODO: shift is now not supported.
    nqshft = 1
    lwf_qshift(:, :) = 0.0_dp

    ! TODO: HeXu: symmetry is not used. Check if it is available.
    ! In  m_phonons, mkphdos, there is the comment:
    ! Rotate e(q) to get e(Sq) to account for symmetrical q-points in BZ.
    ! eigenvectors indeed are not invariant under rotation. See e.g. Eq 39-40 of PhysRevB.76.165108 [[cite:Giustino2007]].
    ! In principle there's a phase due to nonsymmorphic translations but we here need |e(Sq)_iatom|**2
    my_qptopt = 3

    call kpts_ibz_from_kptrlatt(crystal, in_qptrlatt, my_qptopt, nqshft, lwf_qshift, &
         & self%nqibz, self%qibz, self%qweights, nqbz, qbz, new_kptrlatt = new_qptrlatt, &
         & new_shiftk = new_shiftq, bz2ibz = bz2ibz_smap)

    ABI_FREE(bz2ibz_smap)
    !ABI_FREE(bz2ibz)

    self%qptrlatt = new_qptrlatt
    self%shiftq(:) = new_shiftq(:, 1)  ! only one shift in output
    if (self%my_rank == master) then
       write(msg, "(3a, i0)")" LWF ngqpt: ", trim(ltoa(dtset%lwf_ngqpt)), ", qptopt: ", my_qptopt
       call wrtout([ab_out, std_out], msg)
       write(msg, "(2(a, i0))")" Number of q-points in the IBZ: ", self%nqibz, ", number of MPI processes: ", self%nprocs
       call wrtout([ab_out, std_out], msg)

       write(msg, "(a)") " List of q-points: "
       call wrtout([ab_out, std_out], msg)

       write(msg,'(a,i8)')' Grid q points for sampling in the reciprocal space : ',self%nqibz
       call wrtout([ab_out, std_out],msg,'COLL')
       nkpout=self%nqibz
       if(self%nqibz>80)then
         call wrtout([ab_out, std_out],' greater than 80, so only write 20 of them ','COLL')
         nkpout=20
       end if
       do iqpt=1,nkpout
         write(msg, '(1x,i2,a2,3es16.8)' )iqpt,') ',self%qibz(1,iqpt),self%qibz(2,iqpt),self%qibz(3,iqpt)
         call wrtout([ab_out, std_out], msg, 'COLL')
       end do
    end if
    ABI_SFREE(new_shiftq)
    ABI_SFREE(qbz)
  end subroutine prepare_qpoints

  subroutine print_Rlist(self, dtset)
    class(LatticeWannier), intent(inout):: self
    type(anaddb_dataset_type) :: dtset
    integer :: master=0
    integer :: nRout, ii
    character(len = 500):: msg

      if (self%my_rank == master) then
       write(msg, "(2a)")" LWF is transformed to the real space cells: ", trim(ltoa(dtset%lwf_ngqpt))
       call wrtout([ab_out, std_out], msg)
       write(msg, "((a, i0))")" Number of cells: ", self%nR
       call wrtout([ab_out, std_out], msg)

       write(msg,'(a,i8)')' R-vectors for cells in the real space : ',self%nR
       call wrtout([ab_out, std_out],msg,'COLL')
       nRout=self%nR
       if(self%nR>80)then
         call wrtout([ab_out, std_out],' greater than 80, so only write 20 of them ','COLL')
         nRout=20
       end if
       do ii=1,nRout
         write(msg, '(1x,i2,a2,3i8)' )ii,') ',self%Rlist(1,ii),self%Rlist(2,ii),self%Rlist(3,ii)
         call wrtout([ab_out, std_out], msg, 'COLL')
       end do
    end if
  end subroutine print_Rlist




  subroutine prepare_Rlist(self)
    class(LatticeWannier), intent(inout):: self
    integer:: qptrlatt(3), i
    self%nR=1
    do i = 1, 3
       qptrlatt(i) = self%qptrlatt(i, i)
       self%nR=self%nR*qptrlatt(i)
    end do
    call build_Rgrid(qptrlatt, self%Rlist)
  end subroutine prepare_Rlist

  elemental function freq_to_eigenval(f) result (evalue)
    real(dp), intent(in):: f
    real(dp):: evalue
    if(f < -1d-16) then
       evalue = - f*f
    else if (f > 1d-16) then
       evalue = f*f
    else
       evalue = 0.0_dp
    end if
  end function freq_to_eigenval

  subroutine get_ifc_eigens(self, ifc, crystal)
    class(LatticeWannier), intent(inout):: self
    type(ifc_type), intent(in):: ifc
    type(crystal_t), intent(in):: crystal
    real(dp):: eigvec(2, 3, Crystal%natom, 3*Crystal%natom), phfrq(3*Crystal%natom)
    real(dp):: displ(2*3*Crystal%natom*3*Crystal%natom)
    integer:: iq_ibz
    integer:: natom, natom3
    integer:: iatom, iband, i3
    complex(dp):: phase
    natom = crystal%natom
    self%natom = natom
    natom3 = natom*3
    ABI_MALLOC(self%eigenvalues, (natom3, self%nqibz))
    ABI_MALLOC(self%eigenvectors, (natom3, natom3, self%nqibz))
    !print *, "trans:", ifc%trans
    !print *, "xred:", crystal%xred
    do iq_ibz = 1, self%nqibz
       call ifc%fourq(crystal, self%qibz(:,iq_ibz), phfrq, displ, out_eigvec = eigvec)
       ! freqency to eigenvalues
       self%eigenvalues(:, iq_ibz) = freq_to_eigenval(phfrq)
       ! remove phases from eigenvector
       do iatom = 1, natom
           ! to remove the phase factor exp(iqr)
           !phase = exp(-cmplx(0.0_dp, two_pi)*dot_product(crystal%xred(:, iatom), self%qibz(:, iq_ibz) ))
           phase = exp(-cmplx(0.0_dp, two_pi)*dot_product(ifc%trans(:, iatom), self%qibz(:, iq_ibz)))
           do iband=1, natom3
             do i3 = 1, 3
              self%eigenvectors((iatom-1)*3+i3, iband, iq_ibz ) = &
                    &  CMPLX(eigvec(1, i3, iatom, iband), eigvec(2, i3, iatom, iband)) * phase
               end do
           end do
        end do
    end do
  end subroutine get_ifc_eigens

  subroutine write_lwf_nc(self, prefix)
    class(LatticeWannier), intent(inout):: self
    character(len=*), intent(in) ::  prefix
    character(len=500) :: msg 
#ifdef HAVE_NETCDF
    type(IOWannNC):: ncfile
    call self%scdm%create_ncfile(trim(prefix)//"_lwf.nc", ncfile)
    call self%scdm%write_wann_netcdf(ncfile, wannR_unit='dimensionless', HwannR_unit='Ha')
    !NCF_CHECK(self%crystal%ncwrite(ncfile%ncid))
    call self%scdm%close_ncfile(ncfile)
    write(msg, '(a)')  ' LWF construction finished.' 
    call wrtout([ab_out, std_out], msg) 
    write(msg, '(a)')  ' LWF coefficients and Hamiltonian writen to file: '//trim(prefix)//"_lwf.nc ."
    call wrtout([ab_out, std_out], msg) 

#else
    ABI_UNUSED_A(self)
    ABI_UNUSED(prefix)
    NETCDF_NOTENABLED_ERROR()
#endif
  end subroutine write_lwf_nc
  

  subroutine run_all(self, prefix, dtset)
    class(LatticeWannier), intent(inout):: self
    character(len=*), intent(in) ::  prefix
    type(anaddb_dataset_type), intent(in):: dtset
    call self%scdm%construct_wannier()
    call self%print_Rlist(dtset)
    call self%write_lwf_nc(prefix = prefix)
    if(dtset%nqpath>0) then
       call self%write_bands(prefix, dtset)
    end if
  end subroutine run_all

  subroutine run_lattice_wannier(ifc, crystal, dtset, prefix, comm)
    integer, intent(in):: comm
    character(len=*), intent(in):: prefix
    type(ifc_type), intent(in):: ifc
    type(crystal_t), intent(in):: crystal
    type(anaddb_dataset_type), intent(in):: dtset
    type(LatticeWannier):: lwf
    call lwf%initialize(ifc, crystal, dtset, comm)
    call lwf%run_all(prefix, dtset)
    call lwf%finalize()
  end subroutine run_lattice_wannier


  subroutine write_bands(self,  prefix, inp)
    class(LatticeWannier), intent(inout):: self
    character(len=*), intent(in):: prefix
    type(anaddb_dataset_type), intent(in):: inp
    integer,allocatable :: ndiv(:)
    integer :: nfineqpath
    real(dp),allocatable :: fineqpath(:,:)
    real(dp),allocatable :: phfrq(:, :)
    complex(dp), allocatable :: eigvec(:,:, :)
    !real(dp),allocatable :: weights(:)
    integer :: iq, iband
    real(dp) :: f
    ABI_MALLOC(ndiv,(inp%nqpath-1))
    call make_path(inp%nqpath,inp%qpath,self%Crystal%gmet,'G',inp%ndivsm,ndiv,nfineqpath,fineqpath,std_out)
    ABI_FREE(ndiv)

    ABI_MALLOC(phfrq, (inp%lwf_nwann, nfineqpath))
    ABI_MALLOC(eigvec, (inp%lwf_nwann, inp%lwf_nwann, nfineqpath))

    call self%scdm%get_wannier_eigen_klist(fineqpath, nfineqpath, phfrq, eigvec)
    do iq=1, nfineqpath
       do iband =1, inp%lwf_nwann
          f=phfrq(iband, iq)
          if (f>1.0d-9) then
             phfrq(iband, iq)= sqrt(f)
          else if (f<-1.0d-9) then
             phfrq(iband, iq)= -sqrt(-f)
          else
             phfrq(iband, iq)= 0.0_dp
          end if

       end do
    end do

    call write_phfrq(trim(prefix)//"_lwf_PHFRQ", inp%lwf_nwann, nfineqpath, phfrq)

    ABI_FREE(phfrq)
    ABI_FREE(eigvec)
    ABI_FREE(fineqpath)
  end subroutine write_bands

  subroutine write_phfrq(path,nlwf,nqpts,phfreq)
    
    !Arguments ------------------------------------
    !scalars
    integer,intent(in) :: nqpts, nlwf
    character(len=*),intent(in) :: path
    !arrays
    !real(dp),intent(in) :: qpoints(3,nqpts)
    real(dp),intent(in) :: phfreq(nlwf,nqpts)
    
    !Local variables-------------------------------
    !scalars
    integer :: nphmodes, iq, iunit
    !real(dp) :: dummy
    character(len=300) :: formt
    character(len=500) :: msg
    
    ! *************************************************************************
   
    nphmodes = nlwf
    !dummy = qpoints(1,1); dummy = weights(1)
    if (open_file(path, msg, newunit=iunit, form="formatted", status="unknown", action="write") /= 0) then
      ABI_ERROR(msg)
    end if
   
    write (iunit, '(a)')  '# ABINIT generated LWF phonon band structure file. All in Ha atomic units'
    write (iunit, '(a)')  '# '
    write (iunit, '(a,i0)')  '# number_of_qpoints ', nqpts
    write (iunit, '(a,i0)')  '# number_of_phonon_modes ', nphmodes
    write (iunit, '(a)')  '# '
    write (formt,'(a,i0,a)') "(I5, ", nphmodes, "E20.10)"
    do iq= 1, nqpts
      write (iunit, formt)  iq, phfreq(:,iq)
    end do
    close(iunit)

end subroutine write_phfrq



end module m_lwf
!!***
