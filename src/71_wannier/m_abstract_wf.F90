

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

 use m_build_info,      only :  abinit_version
 use defs_wvltypes,  only : wvl_internal_type
 use defs_datatypes, only : pseudopotential_type, ebands_t
 use defs_abitypes, only : MPI_type
 use m_io_tools, only : delete_file, get_unit, open_file
 use m_hide_lapack,     only : matrginv
 use m_fstrings,      only : strcat, sjoin, itoa
 use m_numeric_tools, only : uniformrandom, simpson_int, c2r, l2int
 use m_special_funcs,   only : besjm
 use m_geometry,  only : xred2xcart, rotmat, wigner_seitz
 use m_fftcore,  only : sphereboundary, ngfft_seq, get_kg
 use m_crystal,  only : crystal_t
 use m_ebands,   only : ebands_ncwrite, ebands_expandk, ebands_free
 use m_pawang,   only : pawang_type
 use m_pawrad,   only : pawrad_type, simp_gen
 use m_pawtab,   only : pawtab_type
 use m_pawcprj,  only : pawcprj_type
 use m_paw_sphharm, only : ylm_cmplx, initylmr
 use m_paw_overlap, only : smatrix_pawinit
 use m_pawrhoij, only: pawrhoij_copy
 use m_evdw_wannier, only : evdw_wannier
 use m_fft,            only : fourwf
 use m_wfd, only: wfd_t, wfd_init, wave_t, WFD_STORED

 implicit none

 private


!!***
 type, public:: wann_ksetting_t
   logical :: has_ovikp =  .False.
   type(crystal_t), pointer :: cryst => null()
   integer :: nkpt=0, num_nnmax=0, nntot=0, nsppol=0,rank, comm, nprocs
   integer, allocatable :: ovikp(:, :)
   integer :: my_nspin,  my_nkpt, my_nkpt_pnn
   integer, allocatable :: my_spins(:), my_ikpts(:), my_ikpts_pnn(:)
   real(dp), pointer :: kkpts(:, :)=> null()
   real(dp), allocatable :: my_kkpts(:, :), my_kkpts_pnn(:, :)
 contains
   procedure :: init => wann_ksetting_init
   procedure :: set_ovikp => wann_ksetting_set_ovikp
   procedure :: free => wann_ksetting_free
   procedure :: get_mpw_gmax => wann_ksetting_get_mpw_gmax
   procedure :: distribute_mpi => wann_ksetting_distribute_mpi
   procedure :: get_bks_mask => wann_ksetting_get_bks_mask
 end type wann_ksetting_t



 type, public :: abstract_wf
   logical :: has_paw = .True.
    type(crystal_t), pointer :: cryst => null()
    type(datafiles_type),pointer :: dtfil => null()
    type(dataset_type),pointer :: dtset => null()
    type(hdr_type), pointer :: hdr => null()
    type(mpi_type), pointer :: MPI_enreg => null()
    type(pseudopotential_type), pointer :: psps => null()
    type(pawtab_type),pointer :: pawtab(:) => null()
    type(ebands_t), pointer :: ebands => null()
    type(wann_ksetting_t) :: kset
    integer :: natom=0, nspinor=0, nsppol=0, mband=0, &
         & mkmem=0, nkpt=0, rank=0, nprocs=0, comm=0
  contains
    !procedure :: init => abstract_wf_init
    procedure :: abstract_init
    procedure :: free => abstract_wf_free
    procedure :: cg_elem => abstract_wf_cg_elem
    procedure :: cg_elem_complex => abstract_wf_cg_elem_complex
    procedure :: cprj_elem =>abstract_wf_cprj_elem
    procedure :: load_cg => abstract_wf_load_cg
    procedure :: show_info
    procedure :: get_kgs=>  abstract_wf_get_kgs

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
   ! The working wfd, ebands, etc
    type(wfd_t), pointer :: wfd => null()
    type(wave_t), pointer :: waveprt => null()

    ! if wfd is in IBZ, expand to fullBZ
    logical :: expanded = .False.
    integer, allocatable :: bz2ibz(:, :)
    type(wfd_t), pointer :: wfd_ibz => null()
    type(ebands_t), pointer :: ebands_ibz => null()
    type(hdr_type), pointer :: hdr_ibz => null()

    logical, allocatable :: bks_mask(:, :, :)
    logical, allocatable :: keep_ur(:, :, :)

    type(wfd_t) :: wfd_bz
    type(ebands_t) :: ebands_bz
    type(hdr_type) :: hdr_bz
    type(dataset_type) :: dtset_bz
  contains
    procedure :: init => wfd_wf_init
    procedure :: free => wfd_wf_free
    procedure :: cg_elem => wfd_cg_elem
    procedure :: cg_elem_complex => wfd_cg_elem_complex
    procedure :: cprj_elem =>wfd_cprj_elem
    procedure :: load_cg => wfd_load_cg
 end type wfd_wf

 public :: init_mywfc, compute_iwav, write_cg_and_cprj

contains
!!***


  subroutine wann_ksetting_init(self, cryst, nkpt, &
    & nsppol, kkpts, comm, nprocs, rank)
    class(wann_ksetting_t), intent(inout) :: self
    type(crystal_t), target, intent(in):: cryst
    integer, intent(in) :: nkpt, nsppol, comm, nprocs, rank
    real(dp), target, intent(in) :: kkpts(:, :)
    self%comm=comm
    self%nprocs=nprocs
    self%rank=rank
    self%cryst=>cryst
    self%nkpt=nkpt
    self%rank=rank
    self%nprocs=nprocs
    self%nsppol=nsppol
    self%my_nspin = nsppol
    self%kkpts=> kkpts
  end subroutine wann_ksetting_init



  subroutine wann_ksetting_set_ovikp(self,  ovikp, nntot, num_nnmax)
    class(wann_ksetting_t), intent(inout) :: self
    integer, intent(in) ::  nntot, num_nnmax
    integer, intent(in) :: ovikp(:, :)
    if (self%has_ovikp) then
      ABI_ERROR("ovikp already set!")
    end if
    self%has_ovikp = .True.
    self%nntot = nntot
    self%num_nnmax=num_nnmax
    ABI_MALLOC(self%ovikp, (self%nkpt, num_nnmax))
    self%ovikp(:,:) = ovikp(:, :)
    call self%distribute_mpi()
  end subroutine wann_ksetting_set_ovikp


  subroutine wann_ksetting_free(self)
    class(wann_ksetting_t), intent(inout) :: self
    nullify(self%cryst)
    nullify(self%kkpts)
    if (self%has_ovikp) then
      ABI_SFREE(self%ovikp)
      ABI_FREE(self%my_spins)
      ABI_FREE(self%my_ikpts)
      ABI_FREE(self%my_kkpts)
      ABI_FREE(self%my_ikpts_pnn)
      ABI_FREE(self%my_kkpts_pnn)
    end if


  end subroutine wann_ksetting_free


  subroutine wann_ksetting_distribute_mpi(self)
    class(wann_ksetting_t), intent(inout) :: self
    integer :: ikpt, inn, ik_me, ik_nn, ispin
    logical :: belongs(self%nkpt)
    integer :: counter
    self%my_nkpt=0
    self%my_nkpt_pnn=0

    ABI_MALLOC(self%my_spins, (self%nsppol))
    self%my_nspin=self%nsppol
    do ispin =1, self%nsppol
      self%my_spins(ispin)=ispin
    end do
    ! split kpoints
    call xmpi_split_block(self%nkpt, self%comm, self%my_nkpt, self%my_ikpts)

    ABI_MALLOC(self%my_kkpts, (3, self%my_nkpt))

    belongs(:)=.False.
    do ikpt=1, self%my_nkpt
      ik_me=self%my_ikpts(ikpt)
      self%my_kkpts(:, ikpt) = self%kkpts(:, ik_me)
      belongs(ik_me) = .True.
      do inn=1, self%nntot
        ik_nn=self%ovikp(ikpt, inn)
        belongs(ik_nn) = .True.
      end do
    end do

    self%my_nkpt_pnn=0
    do ikpt =1, self%nkpt
      if (belongs(ikpt)) self%my_nkpt_pnn=self%my_nkpt_pnn+1
    end do

    ABI_MALLOC(self%my_ikpts_pnn, (self%my_nkpt_pnn))
    ABI_MALLOC(self%my_kkpts_pnn, (3, self%my_nkpt_pnn))

    counter=0
    do ikpt =1, self%nkpt
      if (belongs(ikpt)) then
        counter =counter+1
        self%my_ikpts_pnn(counter) = ikpt
        self%my_kkpts_pnn(:, counter) = self%kkpts(:, ikpt)
      end if
    end do

    print *, "my_nspin", self%my_nspin
    print *, "my_spin", self%my_spins
    print *, "my_nkpt", self%my_nkpt
    print *, "my_ikpt", self%my_ikpts
  end subroutine wann_ksetting_distribute_mpi


  subroutine wann_ksetting_get_bks_mask(self, bks_mask, keep_ur, nband, nsppol, keep_ur_value)
    class(wann_ksetting_t), intent(inout) :: self
    integer, intent(in) :: nband, nsppol
    logical, allocatable, intent(inout) :: bks_mask(:, :, :) ! bank, kpt, spin
    logical, allocatable, intent(inout) :: keep_ur(:, :, :)
    logical, intent(in) :: keep_ur_value
    integer :: ikpt, ispin, inn, ik_nn, ik_me

    ABI_MALLOC(bks_mask, (nband, self%nkpt, nsppol))
    ABI_MALLOC(keep_ur, (nband, self%nkpt, nsppol))
    bks_mask(:, :, :) = .False.
    keep_ur(:, :, :) = .False.
    do ispin =1, nsppol
      do ikpt =1 , self%my_nkpt
        ik_me = self%my_ikpts(ikpt)
        bks_mask(:, ik_me, :) =.True.
        keep_ur(:, ik_me, :) = keep_ur_value
        do inn =1, self%nntot
          ik_nn=self%ovikp(ikpt, inn)
          bks_mask(:, ik_nn, ispin)=.True.
          keep_ur(:, ik_nn, ispin) = keep_ur_value
        end do
      end do
    end do
  end subroutine wann_ksetting_get_bks_mask


 subroutine wann_ksetting_get_mpw_gmax(self, ecut, mpw, gmax)
   class(wann_ksetting_t), intent(inout) :: self
   real(dp),intent(in) :: ecut
   integer,intent(out) :: mpw, gmax(3)
   integer,parameter :: istwfk1 = 1
   real(dp) :: cpu, wall, gflops !weight_k,
   !type(gqk_t),pointer :: gqk
   !arrays
   integer :: my_gmax(3), onpw, ipw, ii, my_mpw, ierr
   integer,allocatable :: gtmp(:,:)
   real(dp) :: kk(3), qpt(3)

   integer :: spin, my_ik

!----------------------------------------------------------------------

 mpw = 0; gmax = 0

 ! TODO: This is an hotspot due to the double loop over k and q.
 ! Should use a geometrical approach to compute mpw and gmax.

 call wrtout(std_out, " Computing mpw. This may take some time for dense k/q meshes...")
 !call cwtime(cpu, wall, gflops, "start")

 !do my_is=1,gstore%my_nspins
 !  gqk => gstore%gqk(my_is)
 !  spin = gstore%my_spins(my_is)
 do spin =1, self%nsppol

   do my_ik=1,self%my_nkpt_pnn
     kk = self%my_kkpts_pnn(:, my_ik)

     ! Compute G sphere, returning npw. Note istwfk == 1.
     call get_kg(kk, istwfk1, ecut, self%cryst%gmet, onpw, gtmp)
     mpw = max(mpw, onpw)
     do ipw=1,onpw
       do ii=1,3
         gmax(ii) = max(gmax(ii), abs(gtmp(ii,ipw)))
       end do
     end do
     ABI_FREE(gtmp)
   end do ! my_ik
 end do ! my_is

 my_mpw = mpw; call xmpi_max(my_mpw, mpw, self%comm, ierr)
 my_gmax = gmax; call xmpi_max(my_gmax, gmax, self%comm, ierr)

 call wrtout(std_out, sjoin(' Optimal value of mpw: ', itoa(mpw)))
 !call cwtime_report(" gmax and mpw", cpu, wall, gflops)

end subroutine wann_ksetting_get_mpw_gmax



  subroutine init_mywfc(mywfc, ebands, wfd , cg, cprj, cryst, dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank, comm)
    class(abstract_wf), pointer, intent(inout) :: mywfc
    type(crystal_t), target, intent(in) :: cryst
    type(ebands_t), target, optional, intent(inout) :: ebands
    type(wfd_t), target, optional, intent(inout) :: wfd
    real(dp), target, optional, intent(in):: cg(:, :)
    type(pawcprj_type), target, optional, intent(in):: cprj(:,:)
    type(dataset_type),target, intent(in) :: dtset
    type(datafiles_type),target, intent(in) :: dtfil
    type(mpi_type), target, intent(in) :: MPI_enreg
    type(pseudopotential_type), target, intent(in) :: psps
    type(pawtab_type), target, optional, intent(in) :: pawtab(:)
    type(hdr_type), target, intent(in) :: hdr
    integer, intent(in) :: nprocs, rank, comm

    if(present(cg)) then
       ABI_MALLOC_TYPE_SCALAR(cg_cprj, mywfc)
    else if (present(wfd)) then
       ABI_MALLOC_TYPE_SCALAR(wfd_wf, mywfc)
    end if
    select type(mywfc)
    type is(cg_cprj)
       call mywfc%init( ebands, cg, cprj, cryst, dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank, comm)
    type is(wfd_wf)
       call mywfc%init( ebands, wfd, cryst, dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank, comm)
    end select
  end subroutine init_mywfc

  subroutine abstract_init(self, ebands, cryst, dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank, comm)
    class(abstract_wf), intent(inout) :: self
    type(crystal_t), target, intent(in) :: cryst
    type(ebands_t), target, intent(in) :: ebands
    type(dataset_type),target, intent(in) :: dtset
    type(datafiles_type),target, intent(in) :: dtfil
    type(mpi_type), target, intent(in) :: MPI_enreg
    type(pseudopotential_type), target, intent(in) :: psps
    type(pawtab_type), target, optional, intent(in) :: pawtab(:)
    type(hdr_type), target, intent(in) :: hdr
    integer, intent(in) :: nprocs, rank, comm
    if(present(pawtab)) self%pawtab => pawtab
    self%ebands => ebands
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
    self%comm=comm
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


  subroutine abstract_wf_get_kgs(self, ptr_kg)
    class(abstract_wf), intent(inout) :: self
    integer,  intent(inout) :: ptr_kg(:, :)
    integer :: npw_k, ik, ikg
    real(dp) :: ecut_eff
    integer, allocatable :: kg_k(:,:)
    integer, parameter :: istwfk_1=1
    ptr_kg(:,:)=zero
    !ecut_eff = dtset%ecut * dtset%dilatmx **2
    ecut_eff=self%hdr%ecut_eff
    ikg=0
    do ik=1, self%kset%my_nkpt
      npw_k = self%hdr%npwarr(ik)
      call get_kg(self%ebands%kptns(:,ik),istwfk_1,ecut_eff, &
        & self%cryst%gmet,npw_k,kg_k)
      ptr_kg(:,1+ikg:npw_k+ikg)=kg_k !wfd%Kdata(ik)%kg_k(:,:)
      ikg =ikg+npw_k
      ABI_FREE(kg_k)
    end do
  end subroutine abstract_wf_get_kgs



  subroutine abstract_wf_free(self)
    class(abstract_wf), intent(inout) :: self
    call self%kset%free()
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

  subroutine wfd_wf_free(self)
    class(wfd_wf), intent(inout) :: self
    call self%hdr_bz%free()
    call ebands_free(self%ebands_bz)
    call self%dtset_bz%free()
    ! TODO reenable this
    call self%wfd_bz%free()
    ABI_FREE(self%bz2ibz)
    call self%abstract_wf%free()

    ABI_FREE(self%bks_mask)
    ABI_FREE(self%keep_ur)
  end subroutine wfd_wf_free


  subroutine wfd_wf_init(self, ebands, wfd,cryst, dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank, comm)
    class(wfd_wf), target, intent(inout) :: self
    type(ebands_t), target, intent(in) :: ebands
    type(crystal_t), target, intent(in) :: cryst
    type(wfd_t), target, intent(inout) :: wfd
    type(dataset_type),target, intent(in) :: dtset
    type(datafiles_type),target, intent(in) :: dtfil
    type(mpi_type), target, intent(in) :: MPI_enreg
    type(pseudopotential_type), target, intent(in) :: psps
    type(pawtab_type), target, optional, intent(in) :: pawtab(:)
    type(hdr_type), target, intent(in) :: hdr
    integer, intent(in) :: nprocs, rank, comm


    self%expanded=(dtset%kptopt==1 .or. dtset%kptopt==2)

    if (self%expanded) then
      print *, "Expanding ebands and hdr"
      self%expanded=.True.
      self%ebands_ibz => ebands
      self%wfd_ibz => wfd
      self%hdr_ibz => hdr
      print *, "expanding ebands and hdr"
      call ebands_and_hdr_expandk()
      print *, "expanding dtset"
      call dtset_expandk()
      print *, "expand wfd"
      call self%kset%init(cryst, self%hdr_bz%nkpt, &
        & self%hdr_bz%nsppol, self%hdr_bz%kptns, &
        & self%comm, self%nprocs, self%rank)
      call wfd_expandk()
      self%wfd=> self%wfd_bz
      !self%wfd=> wfd
      call self%abstract_wf%abstract_init(self%ebands_bz, cryst, self%dtset_bz, dtfil, &
        & self%hdr_bz, MPI_enreg, nprocs, psps, pawtab, rank, comm)
    else
      self%expanded=.False.
      self%wfd => wfd
      self%ebands => ebands
      call self%abstract_wf%abstract_init(ebands,cryst, dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank, comm)
    end if
  contains

    subroutine dtset_expandk()
      integer :: nkpt
      self%dtset_bz = dtset%copy()
      !ABI_FREE(self%dtset_bz%kpt)
      ABI_FREE(self%dtset_bz%kptns)
      ABI_FREE(self%dtset_bz%istwfk)
      ABI_FREE(self%dtset_bz%nband)
      nkpt=self%ebands_bz%nkpt
      self%dtset_bz%kptopt=3
      self%dtset_bz%nkpt = nkpt
      !ABI_MALLOC(self%dtset_bz%kpt, (3,nkpt ))
      ABI_MALLOC(self%dtset_bz%kptns, (3,nkpt ))
      ABI_MALLOC(self%dtset_bz%istwfk, (nkpt ))
      ABI_MALLOC(self%dtset_bz%nband, (nkpt ))
      self%dtset_bz%kptns(:,:) = self%ebands_bz%kptns(:,:)
      self%dtset_bz%istwfk(:) = 1.0_dp
      self%dtset_bz%nband(:) = self%hdr_bz%nband(:)
    end subroutine dtset_expandk

    subroutine ebands_and_hdr_expandk()
      real(dp) :: ecut_eff, dksqmax
      type(wvl_internal_type) :: dummy_wvl
      integer:: kptopt3=3
      character(len=200) :: msg
      ! NOTE: is this OK to assume so?
      ecut_eff = dtset%ecut * dtset%dilatmx **2
      call ebands_expandk(inb=ebands, cryst=cryst, ecut_eff=ecut_eff, &
        & force_istwfk1=.True., dksqmax=dksqmax, &
        & bz2ibz=self%bz2ibz, outb=self%ebands_bz)
      if (dksqmax > tol12) then
        write(msg, '(3a,es16.6,4a)' )&
          'At least one of the k points could not be generated from a symmetrical one.',ch10,&
          'dksqmax=',dksqmax,ch10,&
          'Action: check your WFK file and k-point input variables',ch10,&
          '        (e.g. kptopt or shiftk might be wrong in the present dataset or the preparatory one.'
        ABI_ERROR(msg)
      end if

      call hdr_init_lowlvl(self%hdr_bz,self%ebands_bz,psps,pawtab,dummy_wvl,abinit_version,&
        hdr%pertcase,hdr%natom,hdr%nsym,hdr%nspden,hdr%ecut,dtset%pawecutdg,hdr%ecutsm,dtset%dilatmx,&
        hdr%intxc,hdr%ixc,hdr%stmbias,hdr%usewvl,dtset%pawcpxocc,dtset%pawspnorb,dtset%ngfft,dtset%ngfftdg,hdr%so_psp,&
        hdr%qptn,cryst%rprimd,cryst%xred,hdr%symrel,hdr%tnons,hdr%symafm,hdr%typat,hdr%amu,hdr%icoulomb,&
        kptopt3,dtset%nelect,dtset%ne_qFD,dtset%nh_qFD,dtset%ivalence,dtset%cellcharge(1),&
        dtset%kptrlatt_orig,dtset%kptrlatt,&
        dtset%nshiftk_orig,dtset%nshiftk,dtset%shiftk_orig,dtset%shiftk)
      print *, "hdr_bz%npwarr:",self%hdr_bz%npwarr
      ! End CP modified

      if (psps%usepaw == 1) call pawrhoij_copy(hdr%pawrhoij, self%hdr_bz%pawrhoij)
    end subroutine ebands_and_hdr_expandk


    subroutine wfd_expandk()
      integer :: ngfft
      integer, allocatable :: istwfk(:)
      integer :: ik, spin, band
      complex(gwpc), allocatable :: ug(:)
      integer :: g0_k(3), g0_kq(3), g0_q(3), work_ngfft(18),gmax(3),indkk(6,1)
      real(dp),allocatable ::  work(:,:,:,:), cg_kbz(:, :, :)
      integer ::mpw, mband, npw_kbz, size
      integer,allocatable :: kg_kbz(:,:)
      real(dp):: kk_bz(3), kk_ibz(3)

      print *, "dtset%mband:", dtset%mband
      print *, "ebands%nband:", ebands%nband
      print *, "ebands_bz%nband:", self%ebands_bz%nband
      print *, "hdr%mband:", hdr%mband
      print *, "hdr_bz%mband:", self%hdr_bz%mband

      print *, "Here1"
      mband= dtset%mband
      ABI_MALLOC(istwfk, (self%ebands_bz%nkpt))
      istwfk(:) = 1

      block
        integer :: ovikp(self%hdr_bz%nkpt, 1), i
        do i=1, self%hdr_bz%nkpt
          ovikp(i, 1)=i
        end do
        call self%kset%set_ovikp(ovikp, 1 , 1)
      end block

      print *, "Here2"
      !TODO: get mpw, gmax
      call self%kset%get_mpw_gmax(dtset%ecut, mpw, gmax)

      print *, "Here3"
      !mpw = maxval(self%hdr_bz%npwarr)
      mband = self%hdr_bz%mband
      !call gstore%get_mpw_gmax(ecut, mpw, gmax)
      gmax = gmax + 4 ! FIXME: this is to account for umklapp
      gmax = 2*gmax + 1
      call ngfft_seq(work_ngfft, gmax)
      !write(std_out,*)"work_ngfft(1:3): ",work_ngfft(1:3)

      print *, "Here4"
      ABI_MALLOC(work, (2, work_ngfft(4), work_ngfft(5), work_ngfft(6)))
      print *, "alloc kg_kbz"
      ABI_MALLOC(kg_kbz, (3, mpw))
      print *, "alloc cg_kbz"
      ABI_MALLOC(cg_kbz, (2, mpw*self%hdr_bz%nspinor, self%hdr_bz%mband))

      print *, "Here5"
      ! TODO kibz=? ebands%kptns or ebands_bz%kptns

      print *, "nband:", ebands%nband
      print *, "mband", dtset%mband
      print *, "mband", self%hdr_bz%nspinor

      call self%kset%get_bks_mask( bks_mask=self%bks_mask, keep_ur=self%keep_ur, &
        & nband=self%hdr_bz%mband, nsppol=hdr%nsppol, keep_ur_value=.False.)
      self%bks_mask(:,:,:)=.True.
      call wfd_init(wfd=self%wfd_bz,Cryst=cryst,Pawtab=pawtab,Psps=psps, &
        & keep_ur=self%keep_ur,mband=self%hdr_bz%mband,nband=self%ebands_bz%nband, &
        &nkibz=self%ebands_bz%nkpt,nsppol=dtset%nsppol,bks_mask=self%bks_mask,&
        &nspden=dtset%nspden,nspinor=hdr%nspinor,ecut=dtset%ecut, &
        &ecutsm=dtset%ecutsm,dilatmx=dtset%dilatmx, &
        &istwfk=istwfk,kibz=self%ebands_bz%kptns,ngfft=wfd%ngfft, &
        &nloalg=wfd%nloalg,prtvol=dtset%prtvol,pawprtvol=dtset%pawprtvol,comm=comm,&
        &  use_fnl_dir0der0 = .False.) ! optional

      print *, "hdr_bz", self%hdr_bz%npwarr
      print *, "ebands_bz", self%ebands_bz%npwarr
      print *, "wfd_bz", self%wfd_bz%npwarr

      print *, "Here6"
      do spin =1, dtset%nsppol
        do ik=1, self%ebands_bz%nkpt
          npw_kbz=self%hdr_bz%npwarr(ik)
          size =self%hdr_bz%nspinor*npw_kbz
          print *, "size", npw_kbz
          print *, "size", size
          ABI_MALLOC(ug, (size))
          indkk(:, 1) = self%bz2ibz(ik,: )
          !TODO: ecut or ecut_eff?
          print *, "indkk: ", indkk(:, 1)

          kk_bz=self%ebands_bz%kptns(:, ik)
          kk_ibz=self%ebands_ibz%kptns(:, self%bz2ibz(ik, 1))

          print *, "kk_bz:", kk_bz
          print *, "kk_ibz:", kk_ibz

          call wfd%sym_ug_kg(ecut=dtset%ecut, &
            & kk_bz=kk_bz, kk_ibz=kk_ibz, bstart=1, nband=mband, &
            & spin=spin, mpw=mpw, indkk=indkk, cryst=cryst, &
            & work_ngfft=work_ngfft, work=work, istwf_kbz=istwfk(ik), &
            & npw_kbz=npw_kbz, kg_kbz=kg_kbz, cgs_kbz=cg_kbz)
            do band = 1, self%ebands_bz%mband
              ! TODO push_ug
              print *, "npwarr:", self%wfd_bz%npwarr(ik)
              print *, "nspinor:", self%wfd_bz%nspinor
              print *, "size:", size

              ug(:)= cmplx(cg_kbz(1, 1:size,band), cg_kbz(2, 1:size, band), kind=gwpc)
              call self%wfd_bz%push_ug(band, ik, spin, Cryst,ug, &
                & update_ur=.True., update_cprj=.False.)
            end do
            ABI_FREE(ug)
          end do
      end do

      print *, "free cg_kbz, kg_kbz"
      ABI_FREE(cg_kbz)
      ABI_FREE(kg_kbz)
      ABI_FREE(work)
      ABI_FREE(istwfk)
    end subroutine wfd_expandk


  end subroutine wfd_wf_init



#ifdef DONTCOMPILE
  subroutine ibz_to_fullbz(wfd, ebands, dtset, hdr, &
    & cryst,  psps, pawtab, &
    & wfd_bz, ebands_bz, hdr_bz, &
    & bks_mask,  keep_ur, comm, kset)

    type(wfd_t), optional, intent(inout) :: wfd
    type(ebands_t), target, intent(inout) :: ebands
    type(dataset_type),intent(in) :: dtset
    type(hdr_type),  intent(in) :: hdr
    type(crystal_t), intent(in) :: cryst
    type(pseudopotential_type), intent(in) :: psps
    type(pawtab_type), optional,intent(in) :: pawtab(:)
    type(wann_ksetting_t), intent(inout) :: kset
    logical,intent(in) :: bks_mask(:, :, :)
    logical,intent(in) :: keep_ur(:, :, :)
    integer , intent(in)::  comm

    type(wfd_t),   intent(inout) :: wfd_bz
    type(ebands_t), target,  intent(inout) :: ebands_bz
    type(hdr_type),  intent(inout) :: hdr_bz
    logical :: is_fullbz = .True.
    integer, allocatable :: bz2ibz(:,:), indkk(:)
    integer :: nk_ibz, nk_bz
    !real(dp), pointer :: kk_ibz(:, :)=>null(), kk_bz(:, :)=> null()

    integer :: mband,prtvol,pawprtvol
    integer :: nkibz,nsppol,nspden,nspinor
    real(dp) :: ecut,ecutsm,dilatmx

    if(hdr%kptopt /=3 )then
      !call ebands_and_kpt_expandk()

      call hdr_expandk()
      call wfd_expandk()
    end if
  contains



    subroutine hdr_expandk()
      ! see wfk_tofullbz in m_wfk
      type(wvl_internal_type) :: dummy_wvl
      integer:: kptopt3=3
      call hdr_init_lowlvl(hdr_bz,ebands_bz,psps,pawtab,dummy_wvl,abinit_version,&
        hdr%pertcase,hdr%natom,hdr%nsym,hdr%nspden,hdr%ecut,dtset%pawecutdg,hdr%ecutsm,dtset%dilatmx,&
        hdr%intxc,hdr%ixc,hdr%stmbias,hdr%usewvl,dtset%pawcpxocc,dtset%pawspnorb,dtset%ngfft,dtset%ngfftdg,hdr%so_psp,&
        hdr%qptn,cryst%rprimd,cryst%xred,hdr%symrel,hdr%tnons,hdr%symafm,hdr%typat,hdr%amu,hdr%icoulomb,&
        kptopt3,dtset%nelect,dtset%ne_qFD,dtset%nh_qFD,dtset%ivalence,dtset%cellcharge(1),&
        dtset%kptrlatt_orig,dtset%kptrlatt,&
        dtset%nshiftk_orig,dtset%nshiftk,dtset%shiftk_orig,dtset%shiftk)
      ! End CP modified
      if (psps%usepaw == 1) call pawrhoij_copy(hdr%pawrhoij, hdr_bz%pawrhoij)
    end subroutine hdr_expandk



    subroutine wfd_expandk(dtset, ebands,  )
      integer :: ngfft
      integer, allocatable :: istwfk(:)
      integer :: ik, spin, band
      complex, allocatable :: ug(:, :)
      integer :: g0_k(3), g0_kq(3), g0_q(3), work_ngfft(18),gmax(3),indkk(6,1)
      real(dp),allocatable ::  work(:,:,:,:), cg_kbz(:, :, :)
      integer ::mpw, mband, npw_kbz
      integer,allocatable :: kg_kbz(:,:)

      mband= dtset%mband
      nsppol = ebands%nsppol
      nspinor = ebands%nspinor
      ABI_MALLOC(istwfk, (ebands_bz%nkpt))
      istwfk(:) = 1


      !TODO: get mpw, gmax
      call kset%get_mpw_gmax(dtset%ecut, mpw, gmax)

      mpw = maxval(hdr_bz%npwarr)
      mband = hdr_bz%mband
      !call gstore%get_mpw_gmax(ecut, mpw, gmax)
      ! TODO: Init work_ngfft
      !gmax = gmax + 4 ! FIXME: this is to account for umklapp
      !gmax = 2*gmax + 1
      call ngfft_seq(work_ngfft, gmax)
      !write(std_out,*)"work_ngfft(1:3): ",work_ngfft(1:3)

      ABI_MALLOC(work, (2, work_ngfft(4), work_ngfft(5), work_ngfft(6)))
      ABI_MALLOC(kg_kbz, (3, mpw))
      !ABI_MALLOC(cg_kbz, (2, mpw*nspinor, nb))

      ! TODO kibz=? ebands%kptns or ebands_bz%kptns
      call wfd_init(wfd=wfd_bz,Cryst=cryst,Pawtab=pawtab,Psps=psps, &
        & keep_ur=keep_ur,mband=dtset%mband,nband=ebands_bz%nband, &
        &nkibz=ebands_bz%nkpt,nsppol=dtset%nsppol,bks_mask=bks_mask,&
        &nspden=dtset%nspden,nspinor=hdr%nspinor,ecut=dtset%ecut, &
        &ecutsm=dtset%ecutsm,dilatmx=dtset%dilatmx, &
        &istwfk=istwfk,kibz=ebands%kptns,ngfft=wfd%ngfft, &
        &nloalg=wfd%nloalg,prtvol=dtset%prtvol,pawprtvol=dtset%pawprtvol,comm=comm,&
        &  use_fnl_dir0der0 = .False.) ! optional
      do spin =1, dtset%nsppol
        do ik=1, ebands_bz%nkpt
          npw_kbz=hdr_bz%npwarr(ik)
          indkk(:, 1) = bz2ibz(:, ik)
          !TODO: ecut or ecut_eff?
          call wfd%sym_ug_kg(ecut=dtset%ecut, &
            & kk_bz=ebands_bz%kptns, kk_ibz=ebands%kptns, bstart=1, nband=mband, &
            & spin=spin, mpw=mpw, indkk=indkk, cryst=cryst, &
            & work_ngfft=work_ngfft, work=work, istwf_kbz=istwfk(ik), &
            & npw_kbz=npw_kbz, kg_kbz=kg_kbz, cgs_kbz=cg_kbz)
            do band = 1, ebands_bz%mband
              ! TODO push_ug
              !call wfd_bz%push_ug(band, ik, spin, Cryst, ug, &
              !  & update_ur=.True., update_cprj=.False.)
            end do
          end do
      end do

      ABI_FREE(kg_kbz)
      ABI_FREE(work)

    end subroutine wfd_expandk

  end subroutine ibz_to_fullbz
#endif


  subroutine wfd_ibz_to_full_bz(wfc_ibz, dtset,hdr_ibz, psps, pawtab, ebands_ibz)
    type(wfd_t), target, intent(in) :: wfc_ibz
    type(pseudopotential_type),intent(in) :: psps
    type(dataset_type),intent(in) :: dtset
    type(hdr_type), target, intent(in) :: hdr_ibz
    type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)
    type(ebands_t), intent(inout) :: ebands_ibz
    type(crystal_t) :: cryst
    type(hdr_type) :: hdr_kfull
    type(hdr_type),pointer :: ihdr
    type(ebands_t),target :: ebands_full
    type(wvl_internal_type) :: dummy_wvl
    logical,parameter :: force_istwfk1=.True.

    integer,parameter :: formeig0=0,kptopt3=3
    integer :: spin,ikf,ik_ibz,nband_k,mpw_ki,mpw_kf,mband,nspinor,nkfull

    integer :: in_iomode,nsppol,nkibz,out_iomode,isym,itimrev
    integer :: npw_ki,npw_kf,istwf_ki,istwf_kf,ii,jj,iqst,nqst
    real(dp) :: ecut_eff,dksqmax,cpu,wall,gflops

    integer :: g0(3),work_ngfft(18),gmax_ki(3),gmax_kf(3),gmax(3)
    integer,allocatable :: bz2ibz(:,:),kg_ki(:,:),kg_kf(:,:),iperm(:),bz2ibz_sort(:)
    real(dp) :: kf(3),kibz(3)
    real(dp),allocatable :: cg_ki(:,:),cg_kf(:,:),eig_ki(:),occ_ki(:),work(:,:,:,:)
    real(dp), ABI_CONTIGUOUS pointer :: kfull(:,:)
    character(len=500) :: msg
    if (all(dtset%kptrlatt == 0)) then
      write(msg,"(5a)")&
        "Cannot produce full WFK file because kptrlatt == 0",ch10,&
        "Please use nkgpt and shiftk to define a homogeneous k-mesh.",ch10,&
        "Returning to caller"
      ABI_WARNING(msg)
      return
    end if

    mband = hdr_ibz%mband; mpw_ki = maxval(hdr_ibz%npwarr); nkibz = hdr_ibz%nkpt
    nsppol = hdr_ibz%nsppol; nspinor = hdr_ibz%nspinor
    ecut_eff = hdr_ibz%ecut_eff ! ecut * dilatmx**2

    ABI_MALLOC(kg_ki, (3, mpw_ki))
    ABI_MALLOC(cg_ki, (2, mpw_ki*nspinor*mband))
    !ABI_MALLOC(eig_ki, ((2*mband)**iwfk%formeig*mband) )
    ! Here we assume it is GS wavefunction and formeig=0
    ABI_MALLOC(eig_ki, ((2*mband)**0*mband) )
    ABI_MALLOC(occ_ki, (mband))

    cryst = hdr_ibz%get_crystal()

 ! Build new header for owfk. This is the most delicate part since all the arrays in hdr_full
 ! that depend on k-points must be consistent with kfull and nkfull.
 call ebands_expandk(ebands_ibz, cryst, ecut_eff, force_istwfk1, dksqmax, bz2ibz, ebands_full)

 if (dksqmax > tol12) then
   write(msg, '(3a,es16.6,4a)' )&
   'At least one of the k points could not be generated from a symmetrical one.',ch10,&
   'dksqmax=',dksqmax,ch10,&
   'Action: check your WFK file and k-point input variables',ch10,&
   '        (e.g. kptopt or shiftk might be wrong in the present dataset or the preparatory one.'
   ABI_ERROR(msg)
 end if

 nkfull = ebands_full%nkpt
 kfull => ebands_full%kptns

 ! Build new header and update pawrhoij.
 call hdr_init_lowlvl(hdr_kfull,ebands_full,psps,pawtab,dummy_wvl,abinit_version,&
   ihdr%pertcase,ihdr%natom,ihdr%nsym,ihdr%nspden,ihdr%ecut,dtset%pawecutdg,ihdr%ecutsm,dtset%dilatmx,&
   ihdr%intxc,ihdr%ixc,ihdr%stmbias,ihdr%usewvl,dtset%pawcpxocc,dtset%pawspnorb,dtset%ngfft,dtset%ngfftdg,ihdr%so_psp,&
   ihdr%qptn,cryst%rprimd,cryst%xred,ihdr%symrel,ihdr%tnons,ihdr%symafm,ihdr%typat,ihdr%amu,ihdr%icoulomb,&
   kptopt3,dtset%nelect,dtset%ne_qFD,dtset%nh_qFD,dtset%ivalence,dtset%cellcharge(1),&
   dtset%kptrlatt_orig,dtset%kptrlatt,&
   dtset%nshiftk_orig,dtset%nshiftk,dtset%shiftk_orig,dtset%shiftk)

 if (psps%usepaw == 1) call pawrhoij_copy(hdr_ibz%pawrhoij, hdr_kfull%pawrhoij)
  end subroutine wfd_ibz_to_full_bz


  subroutine wfd_free(self)
    class(wfd_wf), intent(inout) :: self
    nullify(self%wfd)
    !ABI_SFREE(self%cg_cache)
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


  subroutine wfd_ug(self, iband, ikpt, isppol, ug )
    class(wfd_Wf), intent(inout) :: self
    integer, intent(in) ::  iband, ikpt, isppol
    real(dp), intent(inout) :: ug(:, :)
    integer :: ik_ibz
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
    call xcopy(npw_k*self%Wfd%nspinor, wave%ug, 1, ug, 1)
  end subroutine wfd_ug


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
    print *, "ik_ibz", ik_ibz
    call self%wfd%get_cprj( band=iband, ik_ibz=ik_ibz, spin=isppol, &
         & Cryst=self%cryst, Cprj_out=cprj_out, sorted=.False.)
    ret=cprj_out(iatom, ispinor)%cp(icplx, ilmn)
  end function wfd_cprj_elem



  subroutine cg_cprj_init(self, ebands, cg, cprj, cryst, dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank, comm)
    class(cg_cprj), intent(inout) :: self
    type(crystal_t), target, intent(in) :: cryst
    type(ebands_t), target, intent(in) :: ebands
    real(dp), target, optional, intent(in):: cg(:, :)
    type(pawcprj_type), target, optional, intent(in):: cprj(:,:)
    type(dataset_type),target, intent(in) :: dtset
    type(datafiles_type),target, intent(in) :: dtfil
    type(mpi_type), target, intent(in) :: MPI_enreg
    type(pseudopotential_type), target, intent(in) :: psps
    type(pawtab_type), target, optional, intent(in) :: pawtab(:)
    type(hdr_type), target, intent(in) :: hdr
    integer, intent(in) :: nprocs, rank, comm
    call self%abstract_wf%abstract_init(ebands, cryst, dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank, comm)
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
    if(self%nprocs>1) then
       call self%remove_tmpfile()
    end if
    nullify(self%cg)
    nullify(self%cprj)
    call self%abstract_wf%free()
    ABI_FREE(self%iwav)
    nullify(self%iwav)
    ABI_FREE(self%icprj)
    call self%abstract_wf%free()
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
