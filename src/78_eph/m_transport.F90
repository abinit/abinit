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
 use m_io_tools,       only : iomode_from_fname, file_exists
 use m_pawang,         only : pawang_type, gauleg
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawfgr,         only : pawfgr_type

 implicit none

 private
!!****

 public transport !! main entry point for transport calculations
!!****

!!****t* m_sigmaph/sigmaph_t
!! NAME
!! transport_rta_t
!!
!! FUNCTION
!! Container for transport quantities in the relaxation time approximation
!!
!! SOURCE

 type,private :: transport_rta_t

   integer :: nsppol
   ! number of spin polarizations

   integer :: ntemp
   ! number of temperatures

   real(dp),allocatable :: tmesh(:)
   ! a list of temperatures at which to compute the transport

   !integer :: nmu
   ! number of dopings

   !real(dp) :: nmesh(:)
   ! a list of carrier concentrations at which to compute transport

   type(ebands_t) :: ebands
   ! container for the bandstructure used to compute the transport properties

   type(edos_t) :: edos
   ! electronic density of states

   real(dp),allocatable :: vvdos(:,:,:,:,:)
   ! velocity density of states
   ! (nw, 2, 0:nsppol, 3, 3)

   real(dp),allocatable :: vvdos_tau(:,:,:,:,:,:)
   ! velocity density of states times the lifetimes for different temperatures
   ! (nw, 2, 0:nsppol, 3, 3, ntemps)

 end type transport_rta_t

!----------------------------------------------------------------------

contains  !=====================================================
!!***

!----------------------------------------------------------------------

!!****f* m_transport/transport
!! NAME
!! sigmaph_t
!!
!! FUNCTION
!! General driver to compute transport properties
!!
!! INPUTS
!!
!! SOURCE

subroutine transport(wfk0_path, ngfft, ngfftf, dtfil, dtset, cryst, pawfgr, pawang, pawrad, pawtab, psps, ebands, comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: wfk0_path
 integer, intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(inout) :: ebands
 type(pseudopotential_type),intent(in) :: psps
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(pawfgr_type),intent(in) :: pawfgr
!arrays
 integer,intent(in) :: ngfft(18),ngfftf(18)

!Local variables ------------------------------
 type(sigmaph_t) :: sigma
 type(wfd_t) :: wfd
 real(dp) :: ecut
 integer, parameter :: opt1 = 1
 integer :: ierr
 integer :: nsppol, nspinor, nkpt, nspden
 integer,allocatable :: nband(:,:), wfd_istwfk(:), indq2ebands(:)
 logical,allocatable :: bks_mask(:,:,:), keep_ur(:,:,:)

 write(*,*) 'Transport computation driver'

 ! compute velocities if not present
 ! initialize important dimensions
 nsppol = ebands%nsppol; nspinor = ebands%nspinor
 nspden = dtset%nspden; nkpt = ebands%nkpt

 ! Initialize the wave function descriptor.
 ! Each node has all k-points and spins and bands between my_bstart and my_bstop
 ABI_MALLOC(nband, (nkpt, nsppol))
 ABI_MALLOC(bks_mask, (dtset%mband, nkpt, nsppol))
 ABI_MALLOC(keep_ur, (dtset%mband, nkpt ,nsppol))

 nband = dtset%mband; bks_mask = .False.; keep_ur = .False.

 write(*,*) 'mband', dtset%mband
 write(*,*) 'nkpt', nkpt
 write(*,*) 'nspinor', nspinor
 write(*,*) 'nspden', nspden
 write(*,*) 'nsppol', nsppol

 bks_mask = .True.

 ! Impose istwfk=1 for all k points. This is also done in respfn (see inkpts)
 ! wfd_read_wfk will handle a possible conversion if WFK contains istwfk /= 1.
 ABI_MALLOC(wfd_istwfk, (nkpt))
 wfd_istwfk = 1

 call wfd_init(wfd,cryst,pawtab,psps,keep_ur,dtset%mband,nband,nkpt,nsppol,bks_mask,&
   nspden,nspinor,ecut,dtset%ecutsm,dtset%dilatmx,wfd_istwfk,ebands%kptns,ngfft,&
   dtset%nloalg,dtset%prtvol,dtset%pawprtvol,comm)

 call wfd%print(header="Wavefunctions for self-energy calculation.",mode_paral='PERS')

 ABI_FREE(nband)
 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)

 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))

! get mapping from sigma to ebands
 sigma = sigmaph_read(dtset,dtfil,xmpi_comm_self,ierr)
 call sigmaph_ebands(sigma,cryst,ebands,opt1,comm,ierr,indq2ebands=indq2ebands)
 write(*,*) indq2ebands


! read lifetimes to array
! compute dos
! compute vvdos
! compute vvdos_tau

end subroutine transport

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

type(transport_rta_t) function transport_rta_new(ebands) result (transport)

!Arguments -------------------------------------
 type(ebands_t),intent(in) :: ebands

 !Allocate important arrays

end function transport_rta_new

!----------------------------------------------------------------------

!!****f* m_transport/transport_rta_compute
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE

subroutine transport_rta_compute(transport_rta, dtset)

!Arguments ------------------------------------
 type(transport_rta_t) :: transport_rta
 type(dataset_type),intent(in) :: dtset

 ! compute dos
 ! compute vvdos
 ! compute tau vvdos
 !edos = ebands_get_dos_matrix_elements(ebands, cryst, &
 !                                      bks_vals, nvals, bks_vecs, nvecs, bks_trans, ntens, )

end subroutine transport_rta_compute

!----------------------------------------------------------------------

!!****f* m_transport/transport_rta_ncwrite
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE

subroutine transport_rta_ncwrite(transport_rta, ncid)

!Arguments --------------------------------------
 type(transport_rta_t),intent(in) :: transport_rta
 integer,intent(in) :: ncid

!Local variables --------------------------------
 integer :: ncerr

! write to netcdf file
 NCF_CHECK(transport_rta%edos%ncwrite(ncid))
 ncerr = nctk_def_arrays(ncid, [ nctkarr_t('vvdos_mesh', "dp", "edos_nw")], defmode=.True.)
 ncerr = nctk_def_arrays(ncid, [ nctkarr_t('vvdos_vals', "dp", "edos_nw, nsppol_plus1, three, three")], defmode=.True.)
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vvdos_vals"), transport_rta%vvdos(:,1,:,:,:)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vvdos_tau"), transport_rta%vvdos_tau(:,1,:,:,:,:)))

end subroutine transport_rta_ncwrite

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

 ! free the allocated arrays and datastructure
 ABI_SFREE(transport_rta%vvdos)
 ABI_SFREE(transport_rta%vvdos_tau)

end subroutine transport_rta_free


end module m_transport
!!***
