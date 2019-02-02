!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_transport
!! NAME
!!  m_transport
!!
!! FUNCTION
!!  Module to compute transport properties using the Boltzmann transport equation (BTE).
!!  Initially for electron mobility limited by electron-phonon scattering.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (HM)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_transport

 use defs_basis
 use defs_abitypes
 use iso_c_binding
 use m_abicore
 use m_xmpi
 use m_errors
 use m_hide_blas
 use m_copy
 use m_ebands
 use m_wfk
 use m_ddk
 use m_wfd
 use m_nctk
 use m_sigmaph
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_crystal,        only : crystal_t
 use m_fstrings,       only : strcat
 use m_io_tools,       only : iomode_from_fname, file_exists
 use m_pawang,         only : pawang_type, gauleg
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawfgr,         only : pawfgr_type

 implicit none

 private
!!****

 public :: transport !! main entry point for transport calculations
!!****

!----------------------------------------------------------------------

!!****t* m_transport/transport_rta_t
!! NAME
!! transport_rta_t
!!
!! FUNCTION
!! Container for transport quantities in the relaxation time approximation
!!
!! SOURCE

type,public :: transport_rta_t


   integer :: nsppol
   ! number of spin polarizations

   integer :: ntemp
   ! number of temperatures

   real(dp),allocatable :: kTmesh(:)
   ! a list of temperatures at which to compute the transport

   real(dp),allocatable :: eminmax_spin(:,:)
   ! min max energy of the of the original ebands object

   real(dp),allocatable :: linewidth_serta(:,:,:,:)
   ! Linewidth computed in the self-energy relaxation time aproximation

   real(dp),allocatable :: linewidth_mrta(:,:,:,:)
   ! Linewidth computed in the momentum relaxation time approximation

   real(dp),allocatable :: velocity(:,:,:,:)
   ! band velocity

   type(gaps_t) :: gaps
   ! get gaps of original ebands object

   !integer :: nmu
   ! number of dopings

   !real(dp) :: nmesh(:)
   ! a list of carrier concentrations at which to compute transport

   type(ebands_t) :: ebands
   ! container for the bandstructure used to compute the transport properties

   type(edos_t) :: edos
   ! electronic density of states

   real(dp),allocatable :: vvdos_mesh(:)
   ! velocity density of states mesh
   ! (nw)

   real(dp),allocatable :: vvdos(:,:,:,:,:,:)
   ! velocity density of states
   ! (nw, 2, 0:nsppol, 3, 3, 1+ntemps)

 end type transport_rta_t
!!***

!----------------------------------------------------------------------

contains  !=====================================================
!!***

!----------------------------------------------------------------------

!!****f* m_transport/transport
!! NAME
!! transport
!!
!! FUNCTION
!! General driver to compute transport properties
!! wk0_path=String with the path to the GS unperturbed WFK file.
!! ngfft(18),ngfftf(18)=Coarse and Fine FFT meshes.
!! dtset<dataset_type>=All input variables for this dataset.
!! ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!! pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!! pawang<pawang_type)>=PAW angular mesh and related data.
!! pawrad(ntypat*usepaw)<pawrad_type>=Paw radial mesh and related data.
!! pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! psps<pseudopotential_type>=Variables related to pseudopotentials.
!! comm=MPI communicator.
!!
!! INPUTS
!!
!! SOURCE

subroutine transport(wfk0_path, ngfft, ngfftf, dtfil, dtset, ebands, cryst, pawfgr, pawang, pawrad, pawtab, psps, comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: wfk0_path
 integer, intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(pseudopotential_type),intent(in) :: psps
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(pawfgr_type),intent(in) :: pawfgr
!arrays
 integer,intent(in) :: ngfft(18),ngfftf(18)

!Local variables ------------------------------
 type(sigmaph_t) :: sigmaph
 type(transport_rta_t) :: transport_rta
 type(wfd_t) :: wfd
 real(dp) :: ecut
 integer :: ierr, my_rank
 integer :: nsppol, nspinor, nkpt, nspden
#ifdef HAVE_NETCDF
 integer :: ncid
#endif
 character(len=fnlen) :: path
 character(len=500) :: msg
 integer,allocatable :: nband(:,:), wfd_istwfk(:)
 logical,allocatable :: bks_mask(:,:,:), keep_ur(:,:,:)

 my_rank = xmpi_comm_rank(comm)
 write(*,*) 'Transport computation driver'

 sigmaph = sigmaph_read(dtset,dtfil,xmpi_comm_self,msg,ierr,keep_open=.true.)
 if (ierr/=0) MSG_ERROR(msg)

! intialize transport
 transport_rta = transport_rta_new(sigmaph,cryst,ebands)
 sigmaph%ncid = nctk_noid
 call sigmaph_free(sigmaph)

 ! Compute transport
 call transport_rta_compute(transport_rta,cryst,dtset,comm)

 ! Master creates the netcdf file used to store the results of the calculation.
#ifdef HAVE_NETCDF
 if (my_rank == 0) then
   path = strcat(dtfil%filnam_ds(4), "_TRANSPORT.nc")
   NCF_CHECK(nctk_open_create(ncid, path, xmpi_comm_self))
   call transport_rta_ncwrite(transport_rta, cryst, ncid)

   ! Close the netcdf file
   NCF_CHECK(nf90_close(ncid))
 end if
#endif

 ! Free memory
 call transport_rta_free(transport_rta)

end subroutine transport
!***

!----------------------------------------------------------------------

!!****f* m_transport/transport
!! NAME
!! sigmaph_t
!!
!! FUNCTION
!! Compute transport quantities in the relaxation time approximation
!!
!! INPUTS
!!
!! SOURCE

type(transport_rta_t) function transport_rta_new(sigmaph,cryst,ebands) result (new)

!Arguments -------------------------------------
 type(sigmaph_t) :: sigmaph
 type(crystal_t) :: cryst
 type(ebands_t),intent(in) :: ebands

!Local variables ------------------------------
 integer :: ierr
 integer,allocatable :: indq2ebands(:)

 ! Allocate important arrays
 ABI_MALLOC(new%kTmesh,(sigmaph%ntemp))
 new%ntemp = sigmaph%ntemp
 new%kTmesh = sigmaph%kTmesh

 new%nsppol = ebands%nsppol
 ABI_MALLOC(new%eminmax_spin,(2,ebands%nsppol))
 new%eminmax_spin = get_minmax(ebands, "eig")

 ierr = get_gaps(ebands,new%gaps)

! read lifetimes to ebands object
 new%ebands = sigmaph_ebands(sigmaph,cryst,ebands,new%linewidth_serta,new%linewidth_mrta,new%velocity,xmpi_comm_self,ierr)

end function transport_rta_new
!!***

!----------------------------------------------------------------------

!!****f* m_transport/transport_rta_compute
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE

subroutine transport_rta_compute(self, cryst, dtset, comm)

!Arguments ------------------------------------
 integer,intent(in) :: comm
 type(transport_rta_t),intent(inout) :: self
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst

!Local variables ------------------------------
 integer :: nsppol, nkpt, mband, ib, ik, spin, ii, jj, itemp
 integer :: ntens, nvecs, nvals, edos_intmeth, ierr
 real(dp) :: vr(3)
 real(dp) :: emin, emax, edos_broad, edos_step
 real(dp) :: linewidth
 real(dp) :: dummy_vals(1,1,1,1), dummy_vecs(1,1,1,1,1)
 real(dp),allocatable :: vvdos_mesh(:), vvdos_tens(:,:,:,:,:,:), vv_tens(:,:,:,:,:,:)
 real(dp),allocatable :: dummy_dosvals(:,:,:,:), dummy_dosvecs(:,:,:,:,:)

 ! create alias for dimensions
 nsppol = self%ebands%nsppol
 nkpt   = self%ebands%nkpt
 mband  = self%ebands%mband
 nvals = 0; nvecs = 0

 ! Allocate vv tensors with and without the lifetimes
 ntens = 1+self%ntemp
 ABI_MALLOC(vv_tens, (3, 3, ntens, mband, nkpt, nsppol))
 do spin=1,nsppol
   do ik=1,nkpt
     do ib=1,mband
       ! Go to cartesian coordinates (same as pmat2cart routine).
       !vr = cryst%rprimd(:,1)*self%ebands%velocity(1,ib,ik,spin) &
       !    +cryst%rprimd(:,2)*self%ebands%velocity(2,ib,ik,spin) &
       !    +cryst%rprimd(:,3)*self%ebands%velocity(3,ib,ik,spin)
       vr(:) = self%velocity(:,ib,ik,spin)
       vr = vr / two_pi
       ! Store in vv_tens
       do ii=1,3
         do jj=1,3
           vv_tens(ii, jj, 1, ib, ik, spin) = vr(ii) * vr(jj)
         end do
       end do
       ! Multiply by the lifetime
       do itemp=1,self%ntemp
         linewidth = abs(self%linewidth_serta(itemp, ib, ik, spin))
         vv_tens(:, :, 1+itemp, ib, ik, spin) = 0
         if (linewidth > tol12) vv_tens(:, :, 1+itemp, ib, ik, spin) = vv_tens(:, :, 1, ib, ik, spin) / linewidth
       end do
     end do
   end do
 end do

 ! Compute DOS and VVDOS
 edos_intmeth = 2
 if (dtset%prtdos == 1) edos_intmeth = 1
 if (dtset%prtdos == -2) edos_intmeth = 3
 edos_step = dtset%dosdeltae; edos_broad = dtset%tsmear
 if (edos_step == 0) edos_step = 0.001

 ! Set default erange
 emin = minval(self%eminmax_spin(1,:)); emin = emin - 0.1_dp * abs(emin)
 emax = maxval(self%eminmax_spin(2,:)); emax = emax + 0.1_dp * abs(emax)

 ! If sigma_erange is set, get emin and emax
 do spin=1,self%ebands%nsppol
   if (dtset%sigma_erange(1) >= zero) emin = self%gaps%vb_max(spin) + tol2 * eV_Ha - dtset%sigma_erange(1)
   if (dtset%sigma_erange(2) >= zero) emax = self%gaps%cb_min(spin) - tol2 * eV_Ha + dtset%sigma_erange(2)
 end do

 ! Compute dos and vvdos multiplied by lifetimes
 self%edos = ebands_get_dos_matrix_elements(self%ebands, cryst, &
                                            dummy_vals, nvals, dummy_vecs, nvecs, vv_tens, ntens, &
                                            edos_intmeth, edos_step, edos_broad, comm, &
                                            self%vvdos_mesh, &
                                            dummy_dosvals, dummy_dosvecs, self%vvdos, emin, emax)

 ! Free memory
 ABI_FREE(vv_tens)

end subroutine transport_rta_compute
!!***

!----------------------------------------------------------------------

!!****f* m_transport/transport_rta_ncwrite
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE

subroutine transport_rta_ncwrite(self, cryst, ncid)

!Arguments --------------------------------------
 type(transport_rta_t),intent(in) :: self
 type(crystal_t),intent(in) :: cryst
 integer,intent(in) :: ncid

!Local variables --------------------------------
 integer :: ncerr

! write to netcdf file
 ncerr = nctk_def_dims(ncid, [ nctkdim_t("ntemp", self%ntemp) ], defmode=.True.)
 NCF_CHECK(self%edos%ncwrite(ncid))
 NCF_CHECK(ebands_ncwrite(self%ebands,ncid))
 NCF_CHECK(cryst%ncwrite(ncid))
 ncerr = nctk_def_arrays(ncid, [nctkarr_t('vvdos_mesh', "dp", "edos_nw")], defmode=.True.)
 ncerr = nctk_def_arrays(ncid, [nctkarr_t('kTmesh', "dp", "ntemp")])
 ncerr = nctk_def_arrays(ncid, [nctkarr_t('vvdos_vals', "dp", "edos_nw, nsppol_plus1, three, three")])
 ncerr = nctk_def_arrays(ncid, [nctkarr_t('vvdos_tau', "dp", "edos_nw, nsppol_plus1, three, three, ntemp")], defmode=.True.)
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kTmesh"), self%kTmesh))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vvdos_mesh"), self%vvdos_mesh))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vvdos_vals"), self%vvdos(:,1,:,:,:,1)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vvdos_tau"),  self%vvdos(:,1,:,:,:,2:)))

end subroutine transport_rta_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_transport/transport_rta_free
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE

subroutine transport_rta_free(transport_rta)

!Arguments --------------------------------------
 type(transport_rta_t),intent(inout) :: transport_rta

 call ebands_free(transport_rta%ebands)
 call transport_rta%gaps%free()
 call edos_free(transport_rta%edos)

 ! free the allocated arrays and datastructure
 ABI_SFREE(transport_rta%vvdos)
 ABI_SFREE(transport_rta%vvdos_mesh)
 ABI_SFREE(transport_rta%kTmesh)
 ABI_SFREE(transport_rta%eminmax_spin)

 ABI_SFREE(transport_rta%velocity)
 ABI_SFREE(transport_rta%linewidth_mrta)
 ABI_SFREE(transport_rta%linewidth_serta)

end subroutine transport_rta_free
!!***

end module m_transport
!!***
