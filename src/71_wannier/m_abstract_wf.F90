

!!****m* ABINIT/m_abstract_wf
!! NAME
!!  m_abstract_wf
!!
!! FUNCTION
!!  Interface with Wannier90
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2022 ABINIT group (BAmadon, CEspejo, FJollet, TRangel, DRH)
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

module m_abstract_wf

 use defs_basis
 use defs_wannier90
 use m_abicore
 use m_errors
 use m_atomdata
 use m_xmpi
 use m_sort
#ifdef FC_NAG
 use f90_unix_dir
#endif
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_nctk
 use m_hdr
 use m_dtset
 use m_dtfil

 use defs_datatypes, only : pseudopotential_type, ebands_t
 use defs_abitypes, only : MPI_type
 use m_io_tools, only : delete_file, get_unit, open_file
 use m_hide_lapack,     only : matrginv
 use m_fstrings,      only : strcat, sjoin
 use m_numeric_tools, only : uniformrandom, simpson_int, c2r, l2int
 use m_special_funcs,   only : besjm
 use m_geometry,  only : xred2xcart, rotmat, wigner_seitz
 use m_fftcore,  only : sphereboundary
 use m_crystal,  only : crystal_t
 use m_ebands,   only : ebands_ncwrite
 use m_pawang,   only : pawang_type
 use m_pawrad,   only : pawrad_type, simp_gen
 use m_pawtab,   only : pawtab_type
 use m_pawcprj,  only : pawcprj_type
 use m_paw_sphharm, only : ylm_cmplx, initylmr
 use m_paw_overlap, only : smatrix_pawinit
 use m_evdw_wannier, only : evdw_wannier
 use m_fft,            only : fourwf
 use m_wfd, only: wfd_t, wfd_init, wave_t, WFD_STORED

 implicit none

 private


!!***
 type, public :: abstract_wf
    logical :: has_paw = .False.
    type(crystal_t), pointer :: cryst => null()
    type(datafiles_type),pointer :: dtfil => null()
    type(dataset_type),pointer :: dtset => null()
    type(hdr_type), pointer :: hdr => null()
    type(mpi_type), pointer :: MPI_enreg => null()
    type(pseudopotential_type), pointer :: psps => null()
    type(pawtab_type),pointer :: pawtab(:)
    integer :: natom=0, nspinor=0, nsppol=0, mband=0, &
         & mkmem=0, nkpt=0, rank=0, nprocs=0
  contains
    !procedure :: init => abstract_wf_init
    procedure :: abstract_init
    procedure :: free => abstract_wf_free
    procedure :: cg_elem => abstract_wf_cg_elem
    procedure :: cg_elem_complex => abstract_wf_cg_elem_complex
    procedure :: cprj_elem =>abstract_wf_cprj_elem
    procedure :: load_cg => abstract_wf_load_cg
    procedure :: show_info
 end type abstract_wf


 type, public,extends(abstract_wf) ::  cg_cprj
    real(dp), pointer :: cg(:, :)=>null()
    type(pawcprj_type), pointer :: cprj(:,:)=>null()
    integer, pointer :: iwav(:,:,:, :)=>null()
    integer, allocatable :: icprj(:, :, :)
  contains
    procedure :: init => cg_cprj_init
    procedure :: free => cg_cprj_free
    procedure :: compute_index_cprj
    procedure :: cg_elem
    procedure :: cg_elem_complex
    procedure :: cprj_elem
    procedure :: write_cg_and_cprj_tmpfile
    procedure :: remove_tmpfile
    procedure :: load_cg
 end type cg_cprj


 type, public, extends(abstract_wf) :: wfd_wf
    logical :: is_fullbz = .False.
    type(wfd_t), pointer :: wfd => null()
    type(wave_t), pointer :: waveprt => null()
    !real(dp), allocatable :: cg_cache(:, :)
    integer :: iband_c=-1, ikpt_c=-1, isppol_c=-1 !index for the cache
  contains
    procedure :: init => wfd_wf_init
    procedure :: cg_elem => wfd_cg_elem
    procedure :: cg_elem_complex => wfd_cg_elem_complex
    procedure :: cprj_elem =>wfd_cprj_elem
    procedure :: load_cg => wfd_load_cg
    !procedure :: build_cache => wfd_build_cache
 end type wfd_wf


 public :: init_mywfc, compute_iwav, write_cg_and_cprj

! public :: mlwfovlp_pw, mlwfovlp_proj, mlwfovlp_projpaw, mlwfovlp_setup, mlwfovlp_seedname
!!***

contains
!!***

  subroutine init_mywfc(mywfc, wfd,  cg, cprj, cryst, dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank)
    class(abstract_wf), pointer, intent(inout) :: mywfc
    type(crystal_t), target, intent(in) :: cryst
    type(wfd_t), target, optional, intent(in) :: wfd
    real(dp), target, optional, intent(in):: cg(:, :)
    type(pawcprj_type), target, optional, intent(in):: cprj(:,:)
    type(dataset_type),target, intent(in) :: dtset
    type(datafiles_type),target, intent(in) :: dtfil
    type(mpi_type), target, intent(in) :: MPI_enreg
    type(pseudopotential_type), target, intent(in) :: psps
    type(pawtab_type), target, optional, intent(in) :: pawtab(:)
    type(hdr_type), target, intent(in) :: hdr
    integer, intent(in) :: nprocs, rank

    if(present(cg)) then
       ABI_MALLOC_TYPE_SCALAR(cg_cprj, mywfc)
    else if (present(wfd)) then
       ABI_MALLOC_TYPE_SCALAR(wfd_wf, mywfc)
    end if
    select type(mywfc)
    type is(cg_cprj)
       call mywfc%init( cg, cprj, cryst, dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank)
    type is(wfd_wf)
       call mywfc%init( wfd, cryst, dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank)
    end select
  end subroutine init_mywfc

  subroutine abstract_init(self,cryst, dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank)
    class(abstract_wf), intent(inout) :: self
    type(crystal_t), target, intent(in) :: cryst
    type(dataset_type),target, intent(in) :: dtset
    type(datafiles_type),target, intent(in) :: dtfil
    type(mpi_type), target, intent(in) :: MPI_enreg
    type(pseudopotential_type), target, intent(in) :: psps
    type(pawtab_type), target, optional, intent(in) :: pawtab(:)
    type(hdr_type), target, intent(in) :: hdr
    integer, intent(in) :: nprocs, rank
    if(present(pawtab)) self%pawtab => pawtab
    self%cryst => cryst
    self%dtset => dtset
    self%dtfil => dtfil
    self%hdr => hdr
    self%MPI_enreg => MPI_enreg
    self%psps => psps
    self%natom = hdr%natom
    self%nspinor = hdr%nspinor
    self%nsppol = hdr%nsppol
    self%mband= hdr%mband
    self%mkmem = dtset%mkmem
    self%nkpt = hdr%nkpt
    self%rank = rank
    self%nprocs = nprocs
  end subroutine abstract_init

  function abstract_wf_cg_elem(self, icplx, ig, ispinor, iband, ikpt, isppol ) result(ret)
    class(abstract_wf), intent(inout) :: self
    integer, intent(in) :: icplx, ig, ispinor, iband, ikpt, isppol
    real(dp) :: ret
    ABI_UNUSED_A(self) 
    ABI_UNUSED(icplx)
    ABI_UNUSED(ig)
    ABI_UNUSED(ispinor)
    ABI_UNUSED(iband)
    ABI_UNUSED(ikpt)
    ABI_UNUSED(isppol)
    ABI_UNUSED(ret)
    ABI_ERROR("Function should be overrided:")
  end function abstract_wf_cg_elem

  function abstract_wf_cg_elem_complex(self,  ig, ispinor, iband, ikpt, isppol ) result(ret)
    class(abstract_wf), intent(inout) :: self
    integer, intent(in) ::  ig, ispinor, iband, ikpt, isppol
    complex(dp) :: ret
    ABI_UNUSED_A(self)
    ABI_UNUSED(ig)
    ABI_UNUSED(ispinor)
    ABI_UNUSED(iband)
    ABI_UNUSED(ikpt)
    ABI_UNUSED(isppol)
    ABI_UNUSED(ret)
    ABI_ERROR("Function should be overrided:")
  end function abstract_wf_cg_elem_complex

  subroutine abstract_wf_load_cg(self, ikpt2, isppol, cg_read)
    class(abstract_wf), intent(inout) :: self
    integer, intent(in) :: ikpt2, isppol
    real(dp), intent(inout) :: cg_read(:, :) ! (2, nspinor*mpw*mband )
    ABI_UNUSED_A(self)
    ABI_UNUSED(ikpt2)
    ABI_UNUSED(isppol)
    ABI_UNUSED(cg_read)
    ABI_ERROR("This function abstarct_wf_load_cg should be overrided.")
  end subroutine abstract_wf_load_cg




  function abstract_wf_cprj_elem(self,icplx,ispinor, iband, ikpt, isppol, iatom, ilmn) result(ret)
    class(abstract_wf), intent(inout) :: self
    integer, intent(in) :: icplx, ispinor, iband, ikpt, isppol, ilmn, iatom
    real(dp) :: ret
    ABI_UNUSED_A(self)
    ABI_UNUSED(icplx)
    ABI_UNUSED(ispinor)
    ABI_UNUSED(iband)
    ABI_UNUSED(ikpt)
    ABI_UNUSED(isppol)
    ABI_UNUSED(iatom)
    ABI_UNUSED(ilmn)
    ABI_UNUSED(ret)
    ABI_ERROR("Function should be overrided:")
  end function abstract_wf_cprj_elem

  subroutine abstract_wf_free(self)
    class(abstract_wf), intent(inout) :: self

    nullify(self%cryst)
    nullify(self%dtset)
    nullify(self%dtfil)
    nullify(self%hdr)
    nullify(self%MPI_enreg)
    nullify(self%psps)
    nullify(self%pawtab)
  end subroutine abstract_wf_free

  subroutine show_info(self)
    class(abstract_wf), intent(inout) :: self
    ! cg_elem
    integer :: ipw=3, ispinor=1, ikpt=2, iband=4, isppol=1
    print *, "========showing wf info=============="
    print *, "ipw:", ipw, "  ispinor:", ispinor, "  ikpt:", ikpt, "  iband:", iband, "  isppol:", isppol
    print *, "cg_elem:", self%cg_elem_complex(ipw, ispinor, ikpt, iband, isppol)
    print *, "====end showing wf info=============="
  end subroutine show_info


  subroutine wfd_wf_init(self, wfd,cryst, dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank)
    class(wfd_wf), intent(inout) :: self
    type(crystal_t), target, intent(in) :: cryst
    type(wfd_t), target, intent(in) :: wfd
    type(dataset_type),target, intent(in) :: dtset
    type(datafiles_type),target, intent(in) :: dtfil
    type(mpi_type), target, intent(in) :: MPI_enreg
    type(pseudopotential_type), target, intent(in) :: psps
    type(pawtab_type), target, optional, intent(in) :: pawtab(:)
    type(hdr_type), target, intent(in) :: hdr
    integer, intent(in) :: nprocs, rank
    self%wfd=> wfd
    call self%abstract_wf%abstract_init(cryst, dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank)
  end subroutine wfd_wf_init




  subroutine wfd_free(self)
    class(wfd_wf), intent(inout) :: self
    nullify(self%wfd)
    !ABI_SFREE(self%cg_cache)
    self%iband_c=-1
    self%ikpt_c=-1
    self%isppol_c=-1
    call self%abstract_wf%free()
  end subroutine wfd_free

!
!  subroutine wfd_build_cache(self, iband, ik, isppol)
!    class(wfd_wf), intent(inout) :: self
!    integer, intent(in) :: iband, ik, isppol
!    integer :: ik_ibz, size
!    ik_ibz=ik
!    size = self%wfd%npwarr(ik_ibz) * self%wfd%nspinor
!    if(allocated(self%cg_cache)) then
!       ABI_FREE(self%cg_cache)
!    end if
!    ABI_MALLOC(self%cg_cache, (2, size))
!    call self%wfd%copy_cg(iband, ik_ibz, isppol, self%cg_cache)
!    self%iband_c = iband
!    self%ikpt_c = ik
!    self%isppol_c = isppol
!  end subroutine wfd_build_cache

  function wfd_cg_elem(self, icplx, ig, ispinor, iband, ikpt, isppol ) result(ret)
    class(wfd_Wf), intent(inout) :: self
    integer, intent(in) :: icplx, ig, ispinor, iband, ikpt, isppol
    integer :: ik_ibz
    real(dp) :: ret
    complex(dp) :: t
    integer :: npw_k
    type(wave_t),pointer :: wave
    character(len=500) :: msg
    ik_ibz=ikpt
    !if(.not. (iband==self%iband_c .and. ikpt==self%ikpt_c .and. isppol==self%isppol_c)) then
       !print *, "Building cache for : ", iband, ikpt, isppol
       !call self%build_cache(iband, ikpt, isppol)
    !end if
    ABI_CHECK(self%wfd%get_wave_ptr(iband, ik_ibz, isppol, wave, msg) == 0, msg)
    if (.not. wave%has_ug == WFD_STORED) then
       write(msg,'(a,i0,a,3i0)')" Node ",self%wfd%my_rank," doesn't have (band,ik_ibz,spin): ",iband,ik_ibz,isppol
       ABI_BUG(msg)
    end if
    npw_k = self%Wfd%npwarr(ik_ibz)
    !call xcopy(npw_k*Wfd%nspinor, wave%ug, 1, ug, 1)
    !ret = self%cg_cache(icplx, ig+self%hdr%npwarr(ik_ibz)*(ispinor-1))
    t = wave%ug(ig+self%hdr%npwarr(ik_ibz)*(ispinor-1))
    select case(icplx)
       case(1)
          ret = real(t)
       case(2)
          ret = aimag(t)
       case default
          ret=-999999.99_dp
       end select

  end function wfd_cg_elem

  function wfd_cg_elem_complex(self,  ig, ispinor, iband, ikpt, isppol ) result(ret)
    class(wfd_wf), intent(inout) :: self
    integer, intent(in) ::  ig, ispinor, iband, ikpt, isppol
    complex(dp) :: ret
    integer :: ik_ibz
    integer :: npw_k
    type(wave_t),pointer :: wave
    character(len=500) :: msg
    ik_ibz=ikpt
    ABI_CHECK(self%wfd%get_wave_ptr(iband, ik_ibz, isppol, wave, msg) == 0, msg)
    if (.not. wave%has_ug == WFD_STORED) then
       write(msg,'(a,i0,a,3i0)')" Node ",self%wfd%my_rank," doesn't have (band,ik_ibz,spin): ",iband,ik_ibz,isppol
       ABI_BUG(msg)
    end if
    npw_k = self%Wfd%npwarr(ik_ibz)
    ret = wave%ug(ig+self%hdr%npwarr(ik_ibz)*(ispinor-1))
  end function wfd_cg_elem_complex

  subroutine wfd_load_cg(self, ikpt2, isppol, cg_read)
    class(wfd_wf), intent(inout) :: self
    integer, intent(in) :: ikpt2, isppol
    real(dp), intent(inout) :: cg_read(:, :)
    integer :: iband, iblk, size
    iblk=0
    size=self%hdr%npwarr(ikpt2) * self%nspinor
    do iband =1, self%mband
        call self%wfd%copy_cg(iband, ikpt2, isppol, cg_read(:, iblk+1:iblk+size))
        iblk = iblk + size
    end do
  end subroutine wfd_load_cg



  function wfd_cprj_elem(self,icplx,ispinor, iband, ikpt, isppol, iatom, ilmn) result(ret)
    class(wfd_wf), intent(inout) :: self
    integer, intent(in) :: icplx, ispinor, iband, ikpt, isppol, ilmn, iatom
    real(dp) :: ret
    type(pawcprj_type) :: cprj_out(self%natom,self%nspinor)

    integer :: ik_ibz
    !TODO: 
    ik_ibz = ikpt
    !call self%wfd%ug2cprj(band=iband,ik_ibz=ik_ibz,spin=ispin,choice=1,idir=0,natom=self%natom,Cryst=self%Cryst ,cwaveprj,sorted=.False.)
    call self%wfd%get_cprj( band=iband, ik_ibz=ik_ibz, spin=isppol, &
         & Cryst=self%cryst, Cprj_out=cprj_out, sorted=.False.)
    ret=cprj_out(iatom, ispinor)%cp(icplx, ilmn)
  end function wfd_cprj_elem



  subroutine cg_cprj_init(self, cg, cprj, cryst, dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank)
    class(cg_cprj), intent(inout) :: self
    type(crystal_t), target, intent(in) :: cryst
    real(dp), target, optional, intent(in):: cg(:, :)
    type(pawcprj_type), target, optional, intent(in):: cprj(:,:)
    type(dataset_type),target, intent(in) :: dtset
    type(datafiles_type),target, intent(in) :: dtfil
    type(mpi_type), target, intent(in) :: MPI_enreg
    type(pseudopotential_type), target, intent(in) :: psps
    type(pawtab_type), target, optional, intent(in) :: pawtab(:)
    type(hdr_type), target, intent(in) :: hdr
    integer, intent(in) :: nprocs, rank
    call self%abstract_wf%abstract_init(cryst, dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank)
    if(present(cg))    self%cg=> cg
    if(present(cprj) .and. present(pawtab)) then
       self%has_paw=.True.
    end if
    if(present(cprj))  self%cprj => cprj
    ABI_MALLOC(self%iwav,(self%nspinor, self%mband,self%nkpt,self%nsppol))
    call compute_iwav(MPI_enreg, dtset, hdr, self%iwav, nprocs, rank)
    call self%compute_index_cprj()
    if(nprocs>1) then
       call self%write_cg_and_cprj_tmpfile()
    end if
    call self%show_info()
  end subroutine cg_cprj_init


  subroutine cg_cprj_free(self)
    class(cg_cprj), intent(inout) :: self
    call self%abstract_wf%free()
    if(self%nprocs>1) then
       call self%remove_tmpfile()
    end if
    nullify(self%cg)
    nullify(self%cprj)
    call self%abstract_wf%free()
    ABI_FREE(self%iwav)
    nullify(self%iwav)
    ABI_FREE(self%icprj)
  end subroutine cg_cprj_free

  function cg_elem(self, icplx, ig, ispinor, iband, ikpt, isppol ) result(ret)
    class(cg_cprj), intent(inout) :: self
    integer, intent(in) :: icplx, ig, ispinor, iband, ikpt, isppol
    integer :: ind
    real(dp) :: ret
    ind=ig+self%iwav(ispinor, iband,ikpt,isppol)
    ret=self%cg(icplx,ind)
  end function cg_elem



  function cg_elem_complex(self, ig,ispinor, iband, ikpt, isppol) result(ret)
    class(cg_cprj), intent(inout) :: self
    integer, intent(in) ::  ig, ispinor, iband, ikpt, isppol
    integer :: ind
    complex(dp) :: ret
    ind=ig+self%iwav(ispinor, iband,ikpt,isppol)
    ret=CMPLX(self%cg(1,ind),  self%cg(2,ind))
  end function cg_elem_complex


  subroutine compute_index_cprj(self)
    ! FIXME:hexu: this is modified from the m_mlwfovlp,
    !     but I think it should be carefully checked.
    ! mcprj=nspinor*mband*mkmem*nsppol
    ! 1. nspinor=2 case seems to be wrong.
    class(cg_cprj), intent(inout) :: self
    integer :: ii, ikpt, isppol, iband
    ABI_MALLOC(self%icprj, (self%mband,self%nkpt,self%nsppol))
    ii=0
    do isppol=1,self%nsppol
       ! FIXME: check if it should be mkmem or nkpt.
       do ikpt=1,self%nkpt
          ! FIXME: nband has the shape of (nsppol*nkpt).
          !nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
          do iband=1,self%dtset%nband(ikpt)
             ! FIXME: hexu: should cycle if the kpt is not in this node??
             ii=ii+1
             self%icprj(iband,ikpt,isppol)=ii
          end do
       end do
    end do
  end subroutine compute_index_cprj


  ! get one element of cprj
  function cprj_elem(self,icplx,ispinor, iband, ikpt, isppol, iatom, ilmn) result(ret)
    class(cg_cprj), intent(inout) :: self
    integer, intent(in) :: icplx, ispinor, iband, ikpt, isppol, ilmn, iatom
    real(dp) :: ret
    integer :: ig
    ! TODO: this seems to be better than compute_index_cprj,
    ! But should it be mband or nband(ikpt)

    ! mcprj=nspinor*mband*mkmem*nsppol
    ! this is the original version in m_mlwfovlp
    !ig=iband+(ikpt-1)*self%mband*self%nspinor + &
    !     &(isppol-1)*self%mkmem*self%mband*self%nspinor

    ig=ispinor+(iband-1)*self%nspinor+(ikpt-1)*self%mband*self%nspinor + &
         &(isppol-1)*self%mkmem*self%mband*self%nspinor
    ret= self%cprj(iatom, ig)%cp(icplx, ilmn)
  end function cprj_elem

  subroutine write_cg_and_cprj_tmpfile(self)
    class(cg_cprj), intent(inout) :: self
    call write_cg_and_cprj(self%dtset, self%cg, self%cprj, self%dtfil, self%iwav, &
         & self%hdr%npwarr, self%mband, self%natom, &
         & self%nsppol, self%nkpt,  self%MPI_enreg, &
         & self%rank, self%psps, self%pawtab)
  end subroutine write_cg_and_cprj_tmpfile

  subroutine remove_tmpfile(self)
    class(cg_cprj), intent(inout) :: self
    integer :: isppol, ikpt, ierr
    integer :: master=1
    character(len=fnlen) :: wfnname
    character(len=500) :: message
    if(self%dtset%prtvol>0) then
       write(message, '(3a)' ) ch10,&
            &       '   mlwfovlp :  Removing temporary files with cg and cprj (PAW)',ch10
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
    end if
    !
    !    Just master  node will remove the files
    !
    if(self%rank==master) then
       do isppol=1,self%nsppol
          do ikpt=1,self%nkpt
             write(wfnname,'(a,I5.5,".",I1)') trim(self%dtfil%fnametmp_cg),ikpt,isppol
             call delete_file(wfnname,ierr)
             if(self%psps%usepaw==1) then
                write(wfnname,'(a,I5.5,".",I1)') trim(self%dtfil%fnametmp_cprj),ikpt,isppol
                call delete_file(wfnname,ierr)
             end if
          end do !ikpt
       end do !isppol
    end if
  end subroutine remove_tmpfile



  subroutine load_cg(self, ikpt2, isppol, cg_read)
    class(cg_cprj), intent(inout) :: self
    integer, intent(in) :: ikpt2, isppol
    real(dp), intent(inout) :: cg_read(:, :)
    character(len=fnlen) :: cg_file
    integer :: npw_k2, ios, ii
    character(len=500) :: message

    integer :: iunit , iband2, ipw, index
    write(cg_file,'(a,I5.5,".",I1)') trim(self%dtfil%fnametmp_cg),ikpt2,isppol
    iunit=1000+ikpt2+ikpt2*(isppol-1)
    npw_k2=self%hdr%npwarr(ikpt2)

    open (unit=iunit, file=cg_file,form='unformatted',status='old',iostat=ios)
    if(ios /= 0) then
       write(message,*) " mlwfovlp_pw: file",trim(cg_file), "not found"
       ABI_ERROR(message)
    end if
    !
    do iband2=1,self%mband
       do ipw=1,npw_k2*self%nspinor
          index=ipw+(iband2-1)*npw_k2*self%nspinor
          read(iunit) (cg_read(ii,index),ii=1,2)
          !            if(me==0 .and. ikpt2==4)write(300,*)'ipw,iband2,index',ipw,iband2,index,cg_read(:,index)
          !            if(me==1 .and. ikpt2==4)write(301,*)'ipw,iband2,index',ipw,iband2,index,cg_read(:,index)
       end do
    end do
    close(iunit)
  end subroutine load_cg


  subroutine compute_iwav(MPI_enreg, dtset, hdr, iwav, nprocs, rank)
    type(mpi_type), intent(in) :: MPI_enreg
    type(dataset_type), intent(in) :: dtset
    type(hdr_type), intent(in) :: hdr
    integer, intent(inout) :: iwav(:, :, :, :)
    ! dimension: (nspinor, mband,nkpt,nsppol))
    integer, intent(in) :: nprocs, rank
    character(len=500) :: message
    integer :: icg(hdr%nsppol, hdr%nkpt)
    integer :: icgtemp, isppol, ikpt, nband_k, npw_k, iband, ispinor

    write(message, '(a,a)' ) ch10,&
         & '   mlwfovlp : compute shifts for g-points '
    call wrtout(std_out,  message,'COLL')
    !----------------------------------------------------------------------
    !Compute shifts for g points (icg,iwav)
    !(here mband is not used, because shifts are internal variables of abinit)
    !----------------------------------------------------------------------
    !write(std_out,*) mpw*dtset%nspinor*mband*mkmem*nsppol
    !ABI_MALLOC(icg,(nsppol,nkpt))
    icg=0
    icgtemp=0
    !ABI_MALLOC(iwav,(dtset%mband,nkpt,nsppol))
    iwav(:,:,:, :)=0
    do isppol=1,hdr%nsppol
       do ikpt=1,hdr%nkpt
          !    MPI:cycle over k-points not treated by this node
          if (nprocs>1 ) then !sometimes we can have just one processor
             if ( ABS(MPI_enreg%proc_distrb(ikpt,1,isppol)-rank)  /=0) CYCLE
          end if

          !    write(std_out,*)'rank',rank,'ikpt',ikpt,'isppol',isppol
          nband_k=dtset%nband(ikpt+(isppol-1)*hdr%nkpt)
          !    write(std_out,*) ikpt+(isppol-1)*nkpt,nkpt
          npw_k=hdr%npwarr(ikpt)
          do iband=1,nband_k
             if(iband.gt. dtset%mband) then
                write(message,'(a,3i0)')" mband",iband,dtset%mband,nband_k
                ABI_ERROR(message)
             end if
             do ispinor =1, dtset%nspinor
                iwav(ispinor, iband,ikpt,isppol)= &
                     !&       (iband-1)*npw_k*dtset%nspinor+icgtemp &
                     & (ispinor-1)*npw_k +icgtemp
             end do
             !icgtemp=icgtemp+ npw_k*dtset%nspinor*nband_k
             icgtemp = icgtemp + npw_k
          end do ! iband
          !    icg(isppol,ikpt)=icgtemp
          !    write(std_out,*) "icg", isppol,ikpt,icg(isppol,ikpt)
       end do  ! ikpt
    end do   ! isppol
    !write(std_out,*) "shift for cg computed"
    !
    !Shifts computed.
  end subroutine compute_iwav


 subroutine write_cg_and_cprj(dtset, cg, cprj, dtfil, iwav, npwarr, mband, natom, &
      &nsppol, nkpt,  MPI_enreg, rank, psps, pawtab)

   type(dataset_type),intent(in) :: dtset
   type(MPI_type),intent(in) :: mpi_enreg
   type(datafiles_type),intent(in) :: dtfil
   integer,intent(in):: iwav(:,:,:,:)
   real(dp),intent(in) :: cg(:, :)
   type(pawcprj_type), intent(in) :: cprj(:, :)
   integer, intent(in) :: rank, nsppol, nkpt, mband, natom
   integer, intent(in) :: npwarr(nkpt)
   type(pseudopotential_type),intent(in) :: psps
   !type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
   type(pawtab_type),intent(in) :: pawtab(:)

   integer :: ikpt, ikpt2, isppol, iun_plot, npw_k, iband,  i, ipw, ispinor, ig
   integer :: iatom,  itypat, lmn_size, ilmn
   character(len=fnlen) :: wfnname
   character(len=1000) :: message

     if(dtset%prtvol>0) then
       write(message, '(3a)' ) ch10,&
&       '   mlwfovlp :  Creating temporary files with cg and cprj (PAW)',ch10
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
     end if
!
     do isppol=1,nsppol
       do ikpt=1,nkpt
!
!        MPI:cycle over k-points not treated by this node
!
          if ( ABS(MPI_enreg%proc_distrb(ikpt,1,isppol)-rank)  /=0) CYCLE

!        write(std_out,*)'writing kpt ',ikpt,'isppol',isppol,' by node ', rank
         write(wfnname,'(a,I5.5,".",I1)') trim(dtfil%fnametmp_cg),ikpt,isppol
         iun_plot=1000+ikpt+ikpt*(isppol-1)

         open (unit=iun_plot, file=wfnname,form='unformatted')
         npw_k=npwarr(ikpt)
         do iband=1,mband
            do ispinor = 1, dtset%nspinor
               do ipw=1,npw_k
                  write(iun_plot) (cg(i,ipw+iwav(ispinor, iband,ikpt,isppol)),i=1,2)
               end do
           end do
         end do
         close(iun_plot)
       end do !ikpt
     end do !isppol
!
!    In the PAW case we also need to write out cprj into files
!
     if(psps%usepaw==1) then
!
!      big loop on atoms, kpts, bands and lmn
        print *, "nsppol", nsppol
        print *, "nkpt", nkpt
        print *, "mkmem", dtset%mkmem
!
       ikpt2=0
       do isppol=1,nsppol
         do ikpt=1,nkpt
!
!          MPI:cycle over k-points not treated by this node
!
            if ( ABS(MPI_enreg%proc_distrb(ikpt,1,isppol)-MPI_enreg%me)  /=0) CYCLE

           ikpt2=ikpt2+1 !sums just on the k-points treated by this node
!
           write(wfnname,'(a,I5.5,".",I1)') trim(dtfil%fnametmp_cprj),ikpt,isppol
           iun_plot=1000+ikpt
           open (unit=iun_plot, file=wfnname,form='unformatted')
!
           do iband=1,mband*dtset%nspinor
             !ig=iband+(ikpt2-1)*mband*dtset%nspinor +(isppol-1)*nkpt*mband*dtset%nspinor !index for cprj(:,ig)
              ! cprj: only mkmem k-points are stored in this node.
             ig=iband+(ikpt2-1)*mband*dtset%nspinor +(isppol-1)*dtset%mkmem*mband*dtset%nspinor !index for cprj(:,ig)
              !
             do iatom=1,natom
               itypat=dtset%typat(iatom)
               lmn_size=pawtab(itypat)%lmn_size
!
               do ilmn=1,lmn_size
                 write(iun_plot) (( cprj(iatom,ig)%cp(i,ilmn)),i=1,2)
               end do !ilmn
             end do !iatom
           end do !iband

           close(iun_plot)
         end do !ikpt
       end do !isppol
     end if !usepaw==1
 end subroutine write_cg_and_cprj


end module m_abstract_wf
