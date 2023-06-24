

!!****m* ABINIT/m_abstract_wf
!! NAME
!!  m_abstract_wf
!!
!! FUNCTION
!!  Interface with Wannier90. 
!!  This module contains the abstract type abstract_wf and its children.
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
   integer :: nkpt=0, mband=0, num_nnmax=0, nntot=0, nsppol=0
   integer :: rank=-999, comm=-999, nprocs=-999
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
         & mkmem=0, nkpt=0, rank=-999, nprocs=-999, comm=-999
  contains
    !procedure :: init => abstract_wf_init
    procedure :: abstract_init
    procedure :: free => abstract_wf_free
    procedure :: cg_elem => abstract_wf_cg_elem
    procedure :: cg_elem_complex => abstract_wf_cg_elem_complex
    procedure :: cprj_elem =>abstract_wf_cprj_elem
    procedure :: get_cg_ptr => abstract_wf_get_cprj_ptr
    procedure :: get_cprj_ptr => abstract_wf_get_cprj_ptr
    procedure :: load_cg => abstract_wf_load_cg
    !procedure :: show_info
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
    procedure :: get_cg_ptr => cg_cprj_get_cprj_ptr
    procedure :: get_cprj_ptr => cg_cprj_get_cprj_ptr
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
    type(MPI_type) :: mpi_enreg_bz
  contains
    procedure :: init => wfd_wf_init
    procedure :: free => wfd_wf_free
    !procedure :: get_ug => wfd_wf_ug
    procedure :: cg_elem => wfd_cg_elem
    procedure :: cg_elem_complex => wfd_cg_elem_complex
    procedure :: cprj_elem =>wfd_cprj_elem
    procedure :: load_cg => wfd_load_cg
 end type wfd_wf

 public :: init_mywfc, compute_iwav, write_cg_and_cprj

contains


  subroutine wann_ksetting_init(self, cryst, nkpt, mband, &
    & nsppol, kkpts, comm, nprocs, rank)
    class(wann_ksetting_t), intent(inout) :: self
    type(crystal_t), target, intent(in):: cryst
    integer, intent(in) :: nkpt, mband, nsppol, comm, nprocs, rank
    real(dp), target, intent(in) :: kkpts(:, :)
    self%comm=comm
    self%nprocs=nprocs
    self%rank=rank
    self%cryst=>cryst
    self%nkpt=nkpt
    self%mband=mband
    self%nsppol=nsppol
    self%my_nspin = nsppol
    self%kkpts=> kkpts
  end subroutine wann_ksetting_init



  subroutine wann_ksetting_set_ovikp(self,  ovikp, nntot, num_nnmax, mpi_enreg)
    class(wann_ksetting_t), intent(inout) :: self
    integer, intent(in) ::  nntot, num_nnmax
    integer, intent(in) :: ovikp(:, :)
    type(mpi_type), intent(inout) :: mpi_enreg
    if (self%has_ovikp) then
      ABI_ERROR("ovikp already set!")
    end if
    !print *, "set_ovikp"
    self%has_ovikp = .True.
    self%nntot = nntot
    self%num_nnmax=num_nnmax
    ABI_MALLOC(self%ovikp, (self%nkpt, num_nnmax))
    self%ovikp(:,:) = ovikp(:, :)
    call self%distribute_mpi(mpi_enreg)
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


  subroutine wann_ksetting_distribute_mpi(self, mpi_enreg)
    class(wann_ksetting_t), intent(inout) :: self
    type(MPI_type), intent(inout) :: mpi_enreg
    integer :: ikpt, inn, ik_me, ik_nn, ispin
    logical :: belongs(self%nkpt)
    integer :: counter
    MPI_enreg%comm_cell=self%comm
    mpi_enreg%me=self%rank
    mpi_enreg%me_kpt=self%rank
    mpi_enreg%nproc = self%nprocs
    MPI_enreg%paral_spinor=0
    !write(std_out,*) "Distributed mpi:", "rank:",  self%rank, "comm", self%comm
    !write(std_out,*) "nprocs:", self%nprocs, "me:", mpi_enreg%me, "me_kpt:", mpi_enreg%me_kpt
    if (.not. allocated(mpi_enreg%proc_distrb))then
      ABI_MALLOC(mpi_enreg%proc_distrb, (self%nkpt, self%mband, self%nsppol) )
    end if
    mpi_enreg%proc_distrb(:, :, :) =-999


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
      MPI_enreg%proc_distrb(ik_me,:,:)= self%rank
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
   !real(dp) :: cpu, wall, gflops !weight_k,
   !type(gqk_t),pointer :: gqk
   !arrays
   integer :: my_gmax(3), onpw, ipw, ii, my_mpw, ierr
   integer,allocatable :: gtmp(:,:)
   real(dp) :: kk(3)

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




subroutine init_mywfc(mywfc, ebands, wfd , cg, cprj, cryst, &
  & dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank, comm)
    class(abstract_wf), pointer, intent(inout) :: mywfc
    type(crystal_t), target, intent(in) :: cryst
    type(ebands_t), target, optional, intent(inout) :: ebands
    type(wfd_t), target, optional, intent(inout) :: wfd
    real(dp), target, optional, intent(in):: cg(:, :)
    type(pawcprj_type), target, optional, intent(in):: cprj(:,:)
    type(dataset_type),target, intent(in) :: dtset
    type(datafiles_type),target, intent(in) :: dtfil
    type(mpi_type), target, intent(inout) :: MPI_enreg
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
      call mywfc%init( ebands=ebands, cg=cg, cprj=cprj, cryst=cryst, dtset=dtset, &
        & dtfil=dtfil, hdr=hdr, MPI_enreg=mpi_enreg, nprocs=nprocs, &
        & psps=psps, pawtab=pawtab, rank=rank, comm=comm)
    type is(wfd_wf)
       call mywfc%init( ebands, wfd, cryst, dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank, comm)
       !call wfd_print_norm(mywfc%wfd, mywfc%hdr)
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

  function abstract_wf_cg_elem(self, icplx, ig, ispinor, iband, ikpt, isppol ) result(res)
    class(abstract_wf), intent(inout) :: self
    integer, intent(in) :: icplx, ig, ispinor, iband, ikpt, isppol
    real(dp) :: res
    ABI_UNUSED_A(self)
    ABI_UNUSED(icplx)
    ABI_UNUSED(ig)
    ABI_UNUSED(ispinor)
    ABI_UNUSED(iband)
    ABI_UNUSED(ikpt)
    ABI_UNUSED(isppol)
    ABI_UNUSED(res)
    ABI_ERROR("Function should be overrided:")
  end function abstract_wf_cg_elem

  function abstract_wf_cg_elem_complex(self,  ig, ispinor, iband, ikpt, isppol ) result(res)
    class(abstract_wf), intent(inout) :: self
    integer, intent(in) ::  ig, ispinor, iband, ikpt, isppol
    complex(dp) :: res
    ABI_UNUSED_A(self)
    ABI_UNUSED(ig)
    ABI_UNUSED(ispinor)
    ABI_UNUSED(iband)
    ABI_UNUSED(ikpt)
    ABI_UNUSED(isppol)
    ABI_UNUSED(res)
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




  function abstract_wf_cprj_elem(self,icplx,ispinor, iband, ikpt, isppol, iatom, ilmn) result(res)
    class(abstract_wf), intent(inout) :: self
    integer, intent(in) :: icplx, ispinor, iband, ikpt, isppol, ilmn, iatom
    real(dp) :: res
    ABI_UNUSED_A(self)
    ABI_UNUSED(icplx)
    ABI_UNUSED(ispinor)
    ABI_UNUSED(iband)
    ABI_UNUSED(ikpt)
    ABI_UNUSED(isppol)
    ABI_UNUSED(iatom)
    ABI_UNUSED(ilmn)
    ABI_UNUSED(res)
    ABI_ERROR("Function should be overrided:")
  end function abstract_wf_cprj_elem

  function abstract_wf_get_cg_ptr(self) result(cg)
    class(abstract_wf), target, intent(inout) :: self
    real(dp), pointer :: cg(:,:)
    ABI_UNUSED_A(self)
    ABI_UNUSED(cg)
    ABI_ERROR("The function abstract_wf%get_cg_ptr is not implemented")
  end function abstract_wf_get_cg_ptr

  function abstract_wf_get_cprj_ptr(self) result(cprj)
    class(abstract_wf), target, intent(inout) :: self
    type(pawcprj_type), pointer :: cprj(:, :)
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(cprj)
    ABI_ERROR("The function abstract_wf%get_cg_ptr is not implemented")
  end function abstract_wf_get_cprj_ptr

  subroutine abstract_wf_get_kgs(self, ptr_kg)
    class(abstract_wf), intent(inout) :: self
    integer,  intent(inout) :: ptr_kg(:, :)
    integer :: npw_k, ik, ikg, ik_me
    real(dp) :: ecut_eff
    integer, allocatable :: kg_k(:,:)
    integer, parameter :: istwfk_1=1
    ptr_kg(:,:)=zero
    !ecut_eff = dtset%ecut * dtset%dilatmx **2
    ecut_eff=self%hdr%ecut_eff
    ikg=0
    do ik=1, self%kset%my_nkpt
      ik_me = self%kset%my_ikpts(ik)
      npw_k = self%hdr%npwarr(ik_me)
      call get_kg(self%ebands%kptns(:,ik_me),istwfk_1,ecut_eff, &
        & self%cryst%gmet,npw_k,kg_k)
      ptr_kg(:,1+ikg:npw_k+ikg)=kg_k(:, :) !wfd%Kdata(ik)%kg_k(:,:)
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

  ! subroutine show_info(self)
  !   class(abstract_wf), intent(inout) :: self
  !   ! cg_elem
  !   integer :: ipw=3, ispinor=1, ikpt=2, iband=4, isppol=1
  !   print *, "========showing wf info=============="
  !   print *, "ipw:", ipw, "  ispinor:", ispinor, "  ikpt:", ikpt, "  iband:", iband, "  isppol:", isppol
  !   print *, "cg_elem:", self%cg_elem_complex(ipw, ispinor, ikpt, iband, isppol)
  !   print *, "====end showing wf info=============="
  ! end subroutine show_info

  subroutine wfd_wf_free(self)
    class(wfd_wf), intent(inout) :: self
    ! TODO reenable this
    if (self%expanded) then
      call self%wfd_bz%free()
      call self%hdr_bz%free()
      call ebands_free(self%ebands_bz)
      call self%dtset_bz%free()
      !call self%mpi_enreg_bz%free()
      ABI_FREE(self%bz2ibz)
      ABI_FREE(self%bks_mask)
      ABI_FREE(self%keep_ur)
    end if
    call self%abstract_wf%free()

  end subroutine wfd_wf_free


  subroutine wfd_wf_init(self, ebands, wfd,cryst, dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank, comm)
    class(wfd_wf), target, intent(inout) :: self
    type(ebands_t), target, intent(in) :: ebands
    type(crystal_t), target, intent(in) :: cryst
    type(wfd_t), target, intent(inout) :: wfd
    type(dataset_type),target, intent(in) :: dtset
    type(datafiles_type),target, intent(in) :: dtfil
    type(mpi_type), target, intent(inout) :: MPI_enreg
    type(pseudopotential_type), target, intent(in) :: psps
    type(pawtab_type), target, optional, intent(in) :: pawtab(:)
    type(hdr_type), target, intent(in) :: hdr
    integer, intent(in) :: nprocs, rank, comm


    self%expanded=(dtset%kptopt==1 .or. dtset%kptopt==2)

    self%comm=comm
    self%rank=rank
    self%nprocs=nprocs
    !print *, "set mpi info to wfd_wf", "self%comm", "self%rank", "self%nprocs"

    if (self%expanded) then
      self%expanded=.True.
      self%ebands_ibz => ebands
      self%wfd_ibz => wfd
      self%hdr_ibz => hdr
      call ebands_and_hdr_expandk()
      call self%kset%init(cryst=cryst,nkpt=self%hdr_bz%nkpt, mband=self%hdr_bz%mband,  &
        & nsppol=self%hdr_bz%nsppol,kkpts= self%hdr_bz%kptns, &
        & comm=self%comm, nprocs=self%nprocs, rank=self%rank)

      self%hdr=> self%hdr_bz
      call set_fake_ovikp()

      call wfd_expandk()
      self%wfd=> self%wfd_bz
      call dtset_expandk()
      call self%abstract_wf%abstract_init(self%ebands_bz, cryst, self%dtset_bz, dtfil, &
        & self%hdr_bz, MPI_enreg, nprocs, psps, pawtab, rank, comm)
    else
      self%expanded=.False.
      self%wfd => wfd
      self%ebands => ebands
      self%hdr=> hdr
      call self%kset%init(cryst=cryst,nkpt=hdr%nkpt, mband=hdr%mband,  &
        & nsppol=hdr%nsppol,kkpts= hdr%kptns, &
        & comm=self%comm, nprocs=self%nprocs, rank=self%rank)
      call set_fake_ovikp()
      call self%abstract_wf%abstract_init(ebands,cryst, dtset, dtfil, hdr, MPI_enreg, nprocs, psps, pawtab, rank, comm)
    end if
  contains

    subroutine set_fake_ovikp()
      integer :: ovikp(self%hdr%nkpt, 1), i
      do i=1, self%hdr%nkpt
        ovikp(i, 1)=i
      end do
      call self%kset%set_ovikp(ovikp, 1 , 1, mpi_enreg)
    end subroutine set_fake_ovikp



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
      self%dtset_bz%mkmem = self%kset%my_nkpt
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
! TODO: test if force_istwfk1 is not set to True, force rotate
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
      ! End CP modified

      if (psps%usepaw == 1) call pawrhoij_copy(hdr%pawrhoij, self%hdr_bz%pawrhoij)
    end subroutine ebands_and_hdr_expandk

    logical  function isirr(ik)
      integer, intent(in) :: ik
      integer :: isym, itimrev, g0(3)
      !ik_ibz = bz2ibz(ikf,1)
      isym = self%bz2ibz(ik,2)
      itimrev = self%bz2ibz(ik,6)
      g0 = self%bz2ibz(ik,3:5)        ! IS(k_ibz) + g0 = k_bz
      isirr = (isym == 1 .and. itimrev == 0 .and. all(g0 == 0))
    end function isirr


    subroutine wfd_expandk()
      integer, allocatable :: istwfk(:)
      integer :: ik, spin, band
      complex(gwpc), allocatable :: ug(:)
      integer ::  work_ngfft(18),gmax(3),indkk(6,1)
      real(dp),allocatable ::  work(:,:,:,:), cg_kbz(:, :, :)
      integer ::mpw, mband, npw_kbz, size, ik_ibz
      integer,allocatable :: kg_kbz(:,:)
      real(dp):: kk_bz(3), kk_ibz(3)

      mband= dtset%mband
      ABI_MALLOC(istwfk, (self%ebands_bz%nkpt))
      istwfk(:) = 1

      call self%kset%get_mpw_gmax(dtset%ecut, mpw, gmax)

      !mpw = maxval(self%hdr_bz%npwarr)
      mband = self%hdr_bz%mband
      !call gstore%get_mpw_gmax(ecut, mpw, gmax)
      gmax = gmax + 4 ! FIXME: this is to account for umklapp
      gmax = 2*gmax + 1
      call ngfft_seq(work_ngfft, gmax)


      ABI_MALLOC(work, (2, work_ngfft(4), work_ngfft(5), work_ngfft(6)))
      ABI_MALLOC(kg_kbz, (3, mpw))
      ABI_MALLOC(cg_kbz, (2, mpw*self%hdr_bz%nspinor, self%hdr_bz%mband))


      !print *, "nband:", ebands%nband
      !print *, "mband", dtset%mband
      !print *, "mband", self%hdr_bz%nspinor

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

      do spin =1, dtset%nsppol
        do ik=1, self%ebands_bz%nkpt
          work(:, :, :, :) =0.0_dp
          kg_kbz(:, :)=0
          cg_kbz(:, :, :)=0.0_dp
          npw_kbz=self%hdr_bz%npwarr(ik)
          size =self%hdr_bz%nspinor*npw_kbz
          !print *, "size", npw_kbz
          !print *, "size", size
          indkk(:, 1) = self%bz2ibz(ik,: )
          !TODO: ecut or ecut_eff?
          !print *, "indkk: ", indkk(:, 1)

          ik_ibz =self%bz2ibz(ik, 1)
          kk_bz=self%ebands_bz%kptns(:, ik)
          kk_ibz=self%ebands_ibz%kptns(:,ik_ibz )


          !Note that we use force_rotate here.
          !Otherwise the sym_ug_kg gives different npw_kbz as in ebands_bz or hdr_bz.
          !if the kpoint is in the IBZ.
          call wfd%sym_ug_kg(ecut=dtset%ecut, &
            & kk_bz=kk_bz, kk_ibz=kk_ibz, bstart=1, nband=mband, &
            & spin=spin, mpw=mpw, indkk=indkk, cryst=cryst, &
            & work_ngfft=work_ngfft, work=work, istwf_kbz=istwfk(ik), &
            & npw_kbz=npw_kbz, kg_kbz=kg_kbz, cgs_kbz=cg_kbz, &
            & force_rotate=.True.)

          self%hdr_bz%npwarr(ik)=npw_kbz
          self%ebands_bz%npwarr(ik)=npw_kbz
          self%wfd_bz%npwarr(ik)=npw_kbz
          size =self%hdr_bz%nspinor*npw_kbz
          ABI_MALLOC(ug, (size))
          do band = 1, self%ebands_bz%nband(ik)
            ug(:)= cmplx(cg_kbz(1, 1:size,band), cg_kbz(2, 1:size, band), kind=gwpc)
            !ug(:) = ug(:) / sqrt(sum(cg_kbz(:, 1:size, band)**2))
            call self%wfd_bz%push_ug(band, ik, spin, Cryst,ug, &
              & update_ur=.True., update_cprj=.False.)
          end do
            ABI_FREE(ug)
          end do
      end do

      !print *, "Norm of wfd"
      !call wfd_print_norm(wfd)
      !print *, "Norm of wfd_bz"
      !call wfd_print_norm(self%wfd_bz)


      !print *, "free cg_kbz, kg_kbz"
      ABI_FREE(cg_kbz)
      ABI_FREE(kg_kbz)
      ABI_FREE(work)
      ABI_FREE(istwfk)
    end subroutine wfd_expandk


  end subroutine wfd_wf_init



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


  ! subroutine wfd_wf_ug(self, iband, ikpt, isppol, ug )
  !   class(wfd_Wf), intent(inout) :: self
  !   integer, intent(in) ::  iband, ikpt, isppol
  !   real(dp), intent(inout) :: ug(:, :)
  !   integer :: ik_ibz
  !   integer :: npw_k
  !   type(wave_t),pointer :: wave
  !   character(len=500) :: msg
  !   ik_ibz=ikpt
  !   !if(.not. (iband==self%iband_c .and. ikpt==self%ikpt_c .and. isppol==self%isppol_c)) then
  !      !print *, "Building cache for : ", iband, ikpt, isppol
  !      !call self%build_cache(iband, ikpt, isppol)
  !   !end if
  !   ABI_CHECK(self%wfd%get_wave_ptr(iband, ik_ibz, isppol, wave, msg) == 0, msg)
  !   if (.not. wave%has_ug == WFD_STORED) then
  !      write(msg,'(a,i0,a,3i0)')" Node ",self%wfd%my_rank," doesn't have (band,ik_ibz,spin): ",iband,ik_ibz,isppol
  !      ABI_BUG(msg)
  !   end if
  !   npw_k = self%Wfd%npwarr(ik_ibz)
  !   call xcopy(npw_k*self%Wfd%nspinor, wave%ug, 1, ug, 1)
  ! end subroutine wfd_wf_ug


  function wfd_cg_elem(self, icplx, ig, ispinor, iband, ikpt, isppol ) result(res)
    class(wfd_Wf), intent(inout) :: self
    integer, intent(in) :: icplx, ig, ispinor, iband, ikpt, isppol
    integer :: ik_ibz
    real(dp) :: res
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
    !res = self%cg_cache(icplx, ig+self%hdr%npwarr(ik_ibz)*(ispinor-1))
    t = wave%ug(ig+self%hdr%npwarr(ik_ibz)*(ispinor-1))
    select case(icplx)
       case(1)
          res = real(t)
       case(2)
          res = aimag(t)
       case default
          res=-999999.99_dp
       end select
  end function wfd_cg_elem

  function wfd_cg_elem_complex(self,  ig, ispinor, iband, ikpt, isppol ) result(res)
    class(wfd_wf), intent(inout) :: self
    integer, intent(in) ::  ig, ispinor, iband, ikpt, isppol
    complex(dp) :: res
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
    res = wave%ug(ig+self%hdr%npwarr(ik_ibz)*(ispinor-1))
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



  function wfd_cprj_elem(self,icplx,ispinor, iband, ikpt, isppol, iatom, ilmn) result(res)
    class(wfd_wf), intent(inout) :: self
    integer, intent(in) :: icplx, ispinor, iband, ikpt, isppol, ilmn, iatom
    real(dp) :: res
    type(pawcprj_type) :: cprj_out(self%natom,self%nspinor)

    integer :: ik_ibz
    !TODO:
    ik_ibz = ikpt
    !call self%wfd%ug2cprj(band=iband,ik_ibz=ik_ibz,spin=ispin,choice=1,idir=0,natom=self%natom,Cryst=self%Cryst ,cwaveprj,sorted=.False.)
    call self%wfd%get_cprj( band=iband, ik_ibz=ik_ibz, spin=isppol, &
         & Cryst=self%cryst, Cprj_out=cprj_out, sorted=.False.)
    res=cprj_out(iatom, ispinor)%cp(icplx, ilmn)
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
    if(present(cprj)) then
      !print *, "before linked", cprj(1,1)%cp(:,:)
      self%cprj => cprj
      !print *, "after linked", self%cprj(1,1)%cp(:,:)
    end if
    ABI_MALLOC(self%iwav,(self%nspinor, self%mband,self%nkpt,self%nsppol))
    call compute_iwav(MPI_enreg, dtset, hdr, self%iwav, nprocs, rank)
    call self%compute_index_cprj()
    if(nprocs>1) then
       call self%write_cg_and_cprj_tmpfile()
    end if
    !call self%show_info()
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

  ! return one entry of cg. 
  ! parameters:
  ! icplx: 1 for real part, 2 for imaginary part
  ! ig: index of G vector
  ! ispinor: index of spinor
  ! iband: index of band
  ! ikpt: index of k point
  ! isppol: index of spin
  function cg_elem(self, icplx, ig, ispinor, iband, ikpt, isppol ) result(res)
    class(cg_cprj), intent(inout) :: self
    integer, intent(in) :: icplx, ig, ispinor, iband, ikpt, isppol
    integer :: ind
    real(dp) :: res
    ind=ig+self%iwav(ispinor, iband,ikpt,isppol)
    res=self%cg(icplx,ind)
  end function cg_elem



  ! return one entry of cg in complex form.
  ! Parameters:
  ! same as cg_elem, except that icplx is not needed.
  function cg_elem_complex(self, ig,ispinor, iband, ikpt, isppol) result(res)
    class(cg_cprj), intent(inout) :: self
    integer, intent(in) ::  ig, ispinor, iband, ikpt, isppol
    integer :: ind
    complex(dp) :: res
    ind=ig+self%iwav(ispinor, iband,ikpt,isppol)
    res=CMPLX(self%cg(1,ind),  self%cg(2,ind), kind=dp)
  end function cg_elem_complex

  !return a pointer to the cg array in the cg_cprj object.
  function cg_cprj_get_cg_ptr(self) result(cg)
    class(cg_cprj), target, intent(inout) :: self
    real(dp), pointer :: cg(:,:)
    cg=> self%cg
  end function cg_cprj_get_cg_ptr

  !return a pointer to the cprj array in the cg_cprj object.
  function cg_cprj_get_cprj_ptr(self) result(cprj)
    class(cg_cprj), target, intent(inout) :: self
    type(pawcprj_type), pointer :: cprj(:, :)
    cprj=>self%cprj
  end function cg_cprj_get_cprj_ptr

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
  function cprj_elem(self,icplx,ispinor, iband, ikpt, isppol, iatom, ilmn) result(res)
    class(cg_cprj), intent(inout) :: self
    integer, intent(in) :: icplx, ispinor, iband, ikpt, isppol, ilmn, iatom
    real(dp) :: res
    integer :: ig
    ! TODO: this seems to be better than compute_index_cprj,
    ! But should it be mband or nband(ikpt)

    ! mcprj=nspinor*mband*mkmem*nsppol
    ! this is the original version in m_mlwfovlp
    !ig=iband+(ikpt-1)*self%mband*self%nspinor + &
    !     &(isppol-1)*self%mkmem*self%mband*self%nspinor

    ig=ispinor+(iband-1)*self%nspinor+(ikpt-1)*self%mband*self%nspinor + &
         &(isppol-1)*self%mkmem*self%mband*self%nspinor
    res= self%cprj(iatom, ig)%cp(icplx, ilmn)
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
       !print *, "nsppol", nsppol
       !print *, "nkpt", nkpt
       !print *, "mkmem", dtset%mkmem
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


 ! subroutine wfd_print_norm(wfd, hdr)
 !   type(wfd_t), intent(in) :: wfd
 !   type(hdr_type), intent(in) :: hdr
 !   integer :: spin, band, ikpt, size
 !   real(dp), allocatable :: cgtemp(:, :)
 !   do spin=1, wfd%nsppol
 !     do ikpt=1, wfd%nkibz
 !       size=wfd%nspinor * hdr%npwarr(ikpt)
 !       ABI_MALLOC(cgtemp, (2,  size))
 !       do band=1, wfd%mband
 !         cgtemp(:, :)=0
 !         !if(isirr(ik)) then
 !         call wfd%copy_cg(band,ikpt, spin, cgtemp)
 !         !print *, "spin:", spin, "band:", band, "ik:", ikpt, "delta:", sum(cgtemp**2)
 !       end do
 !       ABI_FREE(cgtemp)
 !     end do
 !   end do
 ! end subroutine wfd_print_norm


end module m_abstract_wf
