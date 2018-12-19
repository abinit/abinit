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
 type(ebands_t),intent(in) :: ebands
 type(pseudopotential_type),intent(in) :: psps
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(pawfgr_type),intent(in) :: pawfgr
!arrays
 integer,intent(in) :: ngfft(18),ngfftf(18)

!Local variables ------------------------------
 type(sigmaph_t) :: sigma
 type(transport_rta_t) :: transport_rta
 type(wfd_t) :: wfd
 real(dp) :: ecut
 integer :: ierr
 integer :: nsppol, nspinor, nkpt, nspden
 integer,allocatable :: nband(:,:), wfd_istwfk(:)
 logical,allocatable :: bks_mask(:,:,:), keep_ur(:,:,:)

 write(*,*) 'Transport computation driver'

! intialize transport
 transport_rta = transport_rta_new()

! read lifetimes to ebands object
 sigma = sigmaph_read(dtset,dtfil,xmpi_comm_self,ierr)
 transport_rta%ebands = sigmaph_ebands(sigma,cryst,ebands,[1,1],comm,ierr)

 if (ierr == 99) then
   ! compute velocities if not present in the SIGEPH.nc file
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

   !TODO: actually compute velocities
 end if

 call transport_rta_compute(transport_rta,cryst,dtset,comm)

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

type(transport_rta_t) function transport_rta_new() result (transport)

!Arguments -------------------------------------

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

subroutine transport_rta_compute(transport_rta, cryst, dtset, comm)

!Arguments ------------------------------------
 integer,intent(in) :: comm
 type(transport_rta_t),intent(inout) :: transport_rta
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst

!Local variables ------------------------------
 integer :: nsppol, nkpt, mband, ib, ik, spin, ii, jj, itemp
 integer :: ntens, nvecs, nvals, edos_intmeth
 real(dp) :: eminmax_spin(2,transport_rta%ebands%nsppol), vr(3)
 real(dp) :: emin, emax, edos_broad, edos_step
 real(dp),allocatable :: dummy_vals(:,:,:,:), dummy_vecs(:,:,:,:,:), vv_tens(:,:,:,:,:,:)
 real(dp),allocatable :: vvdos_mesh(:), vvdos_tens(:,:,:,:,:,:)
 real(dp),allocatable :: dummy_dosvals(:,:,:,:), dummy_dosvecs(:,:,:,:,:)

 ! create alias for dimensions
 nsppol = transport_rta%ebands%nsppol
 nkpt   = transport_rta%ebands%nkpt
 mband  = transport_rta%ebands%mband
 nvals = 0; nvecs = 0

 ! Allocate vv tensors with and without the lifetimes
 ntens = 1+transport_rta%ntemp
 ABI_MALLOC(vv_tens, (3, 3, ntens, mband, nkpt, nsppol))
 do spin=1,nsppol
   do ik=1,nkpt
     do ib=1,mband
       ! Go to cartesian coordinates (same as pmat2cart routine).
       vr = cryst%rprimd(:,1)*transport_rta%ebands%velocity(1,ib,ik,spin) &
           +cryst%rprimd(:,2)*transport_rta%ebands%velocity(2,ib,ik,spin) &
           +cryst%rprimd(:,3)*transport_rta%ebands%velocity(3,ib,ik,spin)
       vr = vr / two_pi
       ! Store in vv_tens
       do ii=1,3
         do jj=1,3
           vv_tens(ii, jj, 1, ib, ik, spin) = vr(ii) * vr(jj)
         end do
       end do
       ! Multiply by the lifetime
       do itemp=1,transport_rta%ntemp
         vv_tens(:, :, 1+itemp, ib, ik, spin) = vv_tens(:, :, 1, ib, ik, spin) / &
                                                transport_rta%ebands%linewidth(itemp, ib, ik, spin)
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

 !set default erange
 eminmax_spin = get_minmax(transport_rta%ebands, "eig")
 emin = minval(eminmax_spin(1,:)); emin = emin - 0.1_dp * abs(emin)
 emax = maxval(eminmax_spin(2,:)); emax = emax + 0.1_dp * abs(emax)

 !compute dos and vvdos multiplied by lifetimes
 transport_rta%edos = ebands_get_dos_matrix_elements(transport_rta%ebands, cryst, &
                                       dummy_vals, nvals, dummy_vecs, nvecs, vv_tens, ntens, &
                                       edos_intmeth, edos_step, edos_broad, comm, vvdos_mesh, &
                                       dummy_dosvals, dummy_dosvecs, vvdos_tens, emin, emax )

 !unpack the data

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
