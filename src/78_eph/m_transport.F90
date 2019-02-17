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
 use m_nctk
 use m_sigmaph
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_crystal,        only : crystal_t
 use m_numeric_tools,  only : arth, simpson_int, polyn_interp
 use m_fstrings,       only : strcat
 use m_occ,            only : occ_fd, occ_dfd
 use m_pawang,         only : pawang_type
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

   integer :: nspinor
   ! number of spin components

   integer :: ntemp
   ! number of temperatures

   integer :: ndop
   ! number of carrier concentrarions at which to evaluate Fermi energy

   integer :: nw
   ! number of frequencies at which trasport quantities are computed

   real(dp),allocatable :: kTmesh(:)
   ! a list of temperatures at which to compute the transport

   real(dp),allocatable :: carriers(:)
   ! a list of different carrier concentrarions (number of electrons per unit cell)

   real(dp),allocatable :: mu_carriers(:,:)
   ! Chemical potential at this carrier concentrarion and temperature
   ! (ntemp, ndop)

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
   ! (nw, 2, 0:nsppol, 3, 3, 1+ntemp)

   real(dp),allocatable :: n(:,:)
   ! (nw,ntemp) carrier density (n/cm^3)

   real(dp),allocatable :: l0(:,:,:,:,:)
   real(dp),allocatable :: l1(:,:,:,:,:)
   real(dp),allocatable :: l2(:,:,:,:,:)
   ! onsager coeficients
   ! (nw, nsppol, 3, 3, ntemp)

   real(dp),allocatable :: sigma(:,:,:,:,:)
   real(dp),allocatable :: mobility(:,:,:,:,:)
   real(dp),allocatable :: seebeck(:,:,:,:,:)
   real(dp),allocatable :: kappa(:,:,:,:,:)
   real(dp),allocatable :: pi(:,:,:,:,:)
   ! transport coefficients
   ! (nw, nsppol, 3, 3, ntemp)

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
 integer :: ii, ierr, my_rank
 integer :: nsppol, nspinor, nspden
#ifdef HAVE_NETCDF
 integer :: ncid
#endif
 character(len=fnlen) :: path
 character(len=500) :: msg
 !integer,allocatable :: nband(:,:), wfd_istwfk(:)
 !logical,allocatable :: bks_mask(:,:,:), keep_ur(:,:,:)

 my_rank = xmpi_comm_rank(comm)
 call wrtout(std_out, 'Transport computation driver')

 sigmaph = sigmaph_read(dtset,dtfil,xmpi_comm_self,msg,ierr,keep_open=.true.)
 if (ierr/=0) MSG_ERROR(msg)

 ! Intialize transport
 transport_rta = transport_rta_new(dtset,sigmaph,cryst,ebands)
 sigmaph%ncid = nctk_noid
 call sigmaph_free(sigmaph)

 ! Compute transport
 call transport_rta_compute(transport_rta,cryst,dtset,comm)

 ! Write transport to stdout (for the test suite)
 call transport_rta_write(transport_rta)

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
!!***

!----------------------------------------------------------------------

!!****f* m_transport/transport_rta_new
!! NAME
!! transport_rta_new
!!
!! FUNCTION
!! Compute transport quantities in the relaxation time approximation
!!
!! INPUTS
!!
!! SOURCE

type(transport_rta_t) function transport_rta_new(dtset,sigmaph,cryst,ebands) result (new)

!Arguments -------------------------------------
 type(sigmaph_t),intent(in) :: sigmaph
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(dataset_type),intent(in) :: dtset

!Local variables ------------------------------
 integer,parameter :: occopt3=3, ndop=10
 real(dp), parameter :: doplist(ndop) = [-1d20,-1d19,-1d18,-1d17,-1d16,1d16,1d17,1d18,1d19,1d20]
 type(ebands_t) :: tmp_ebands
 character(len=500) :: msg
 integer :: itemp, ii, ierr, spin
 !integer,allocatable :: indq2ebands(:)
 real(dp) :: extrael

 ! Allocate temperature arrays
 new%ntemp = sigmaph%ntemp
 ABI_MALLOC(new%kTmesh,(new%ntemp))
 new%kTmesh = sigmaph%kTmesh

 ! Information about the Gaps
 new%nsppol = ebands%nsppol
 new%nspinor = ebands%nspinor
 ABI_MALLOC(new%eminmax_spin,(2,ebands%nsppol))
 new%eminmax_spin = get_minmax(ebands, "eig")

 ierr = get_gaps(ebands,new%gaps)
 if (ierr /= 0) then
   do spin=1, ebands%nsppol
     MSG_WARNING(trim(new%gaps%errmsg_spin(spin)))
   end do
   MSG_WARNING("get_gaps returned non-zero exit status. See above warning messages...")
 end if

 ! Read lifetimes to ebands object
 new%ebands = sigmaph_ebands(sigmaph,cryst,ebands,new%linewidth_serta,new%linewidth_mrta,new%velocity,xmpi_comm_self,ierr)
 ABI_CHECK(allocated(new%velocity),'Could not read velocities from SIGEPH.nc file')

 ! Compute Fermi for different dopings
#if 0
 new%ndop = ndop
 ABI_MALLOC(new%mu_carriers,(new%ntemp,new%ndop))
 ABI_MALLOC(new%carriers,(new%ndop))
#define ucvolcm3 (cryst%ucvol / (Bohr_Ang * 1d-8)**3)
 new%carriers = doplist * ucvolcm3
 if (new%ndop > 1) then
   call ebands_copy(ebands, tmp_ebands)
   do itemp=1,new%ntemp
     call ebands_set_scheme(tmp_ebands, occopt3, new%kTmesh(itemp), dtset%spinmagntarget, dtset%prtvol)
     do ii=1,new%ndop
       extrael = new%carriers(ii)
       call ebands_set_nelect(tmp_ebands, ebands%nelect + extrael, dtset%spinmagntarget, msg)
       call wrtout(ab_out,msg)
       new%mu_carriers(itemp,ii) = tmp_ebands%fermie
     end do
   end do
   call ebands_free(tmp_ebands)
 end if
#else
 ! Same doping case as sigmaph
 new%ndop = 1
 ABI_MALLOC(new%mu_carriers,(new%ntemp,new%ndop))
 ABI_MALLOC(new%carriers,(new%ndop))
 new%carriers = [dtset%eph_extrael]
 new%mu_carriers(:,1) = sigmaph%mu_e
#endif

end function transport_rta_new
!!***

!----------------------------------------------------------------------

!!****f* m_transport/transport_rta_compute
!! NAME
!! transport_rta_compute
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
 integer :: nsppol, nkpt, mband, ib, ik, iw, ispin, ii, jj, itemp
 integer :: ntens, nvecs, nvals, edos_intmeth
 real(dp) :: vr(3)
 real(dp) :: emin, emax, edos_broad, edos_step, max_occ, e, kT
 real(dp) :: linewidth, fact, n0, n0_dy
 real(dp) :: dummy_vals(1,1,1,1), dummy_vecs(1,1,1,1,1)
 real(dp),allocatable :: vv_tens(:,:,:,:,:,:) !vvdos_mesh(:), vvdos_tens(:,:,:,:,:,:),
 real(dp),allocatable :: dummy_dosvals(:,:,:,:), dummy_dosvecs(:,:,:,:,:)

 ! create alias for dimensions
 nsppol = self%ebands%nsppol
 nkpt   = self%ebands%nkpt
 mband  = self%ebands%mband
 nvals = 0; nvecs = 0

 ! Allocate vv tensors with and without the lifetimes
 ntens = 1+self%ntemp
 ABI_MALLOC(vv_tens, (3, 3, ntens, mband, nkpt, nsppol))
 do ispin=1,nsppol
   do ik=1,nkpt
     do ib=1,mband
       vr(:) = self%velocity(:,ib,ik,ispin)
       ! Store in vv_tens
       do ii=1,3
         do jj=1,3
           vv_tens(ii, jj, 1, ib, ik, ispin) = vr(ii) * vr(jj)
         end do
       end do
       ! Multiply by the lifetime
       do itemp=1,self%ntemp
         linewidth = abs(self%linewidth_serta(itemp, ib, ik, ispin))
         vv_tens(:, :, 1+itemp, ib, ik, ispin) = 0
         if (linewidth > tol12) vv_tens(:, :, 1+itemp, ib, ik, ispin) = vv_tens(:, :, 1, ib, ik, ispin) / linewidth
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
 do ispin=1,self%ebands%nsppol
   if (dtset%sigma_erange(1) >= zero) emin = self%gaps%vb_max(ispin) + tol2 * eV_Ha - dtset%sigma_erange(1)
   if (dtset%sigma_erange(2) >= zero) emax = self%gaps%cb_min(ispin) - tol2 * eV_Ha + dtset%sigma_erange(2)
 end do

 ! Compute dos and vvdos multiplied by lifetimes
 self%edos = ebands_get_dos_matrix_elements(self%ebands, cryst, &
                                            dummy_vals, nvals, dummy_vecs, nvecs, vv_tens, ntens, &
                                            edos_intmeth, edos_step, edos_broad, comm, &
                                            self%vvdos_mesh, &
                                            dummy_dosvals, dummy_dosvecs, self%vvdos, emin=emin, emax=emax)
 self%nw = self%edos%nw
 edos_step = self%edos%step

 ! Free memory
 ABI_SFREE(dummy_dosvals)
 ABI_SFREE(dummy_dosvecs)
 ABI_FREE(vv_tens)

 ! Transport coeficients computation
 ABI_MALLOC(self%n,   (self%nw,self%ntemp))

 ABI_MALLOC(self%l0,(self%nw,self%nsppol,3,3,self%ntemp+1))
 ABI_MALLOC(self%l1,(self%nw,self%nsppol,3,3,self%ntemp+1))
 ABI_MALLOC(self%l2,(self%nw,self%nsppol,3,3,self%ntemp+1))

 ABI_MALLOC(self%sigma,   (self%nw,self%nsppol,3,3,self%ntemp+1))
 ABI_MALLOC(self%mobility,(self%nw,self%nsppol,3,3,self%ntemp+1))
 ABI_MALLOC(self%seebeck, (self%nw,self%nsppol,3,3,self%ntemp+1))
 ABI_MALLOC(self%kappa,   (self%nw,self%nsppol,3,3,self%ntemp+1))
 ABI_MALLOC(self%pi,      (self%nw,self%nsppol,3,3,self%ntemp+1))

 ! Compute onsager coefficients
 call onsager(self%vvdos(:,1,1:,:,:,:),0,self%l0)
 call onsager(self%vvdos(:,1,1:,:,:,:),1,self%l1)
 call onsager(self%vvdos(:,1,1:,:,:,:),2,self%l2)

 ! Compute transport quantities
#define siemens (e_Cb**2/Ha_J/Time_Sec)
#define meter   (Bohr_Ang * 1d-10)
#define second  (Time_Sec)
#define volt    (Ha_J/e_Cb)
#define fact0   (siemens / (meter*second) / cryst%ucvol)
#define fact1   (volt * fact0)
#define fact2   (volt**2 * fact0)
 self%sigma = fact0 * self%l0
 self%pi(:,:,:,:,itemp) = (fact0 * self%l0(:,:,:,:,itemp)) / (fact1 * self%l1(:,:,:,:,itemp))
 do itemp=1,self%ntemp
   kT = self%kTmesh(itemp) / kb_HaK
   self%seebeck(:,:,:,:,itemp) = 1/kT * (fact1 * self%l1(:,:,:,:,itemp))/ &
                                        (fact0 * self%l0(:,:,:,:,itemp))
   self%kappa(:,:,:,:,itemp) = 1/kT * (-(fact1 * self%l1(:,:,:,:,itemp))**2 / &
                                        (fact0 * self%l0(:,:,:,:,itemp)) + &
                                        (fact2 * self%l2(:,:,:,:,itemp)))
 end do

 self%mobility = 0
 max_occ = two/(self%nspinor*self%nsppol)
 do ispin=1,self%nsppol
   do itemp=1,self%ntemp
     ! Compute carrier density
     kT = self%kTmesh(itemp)
     do iw=1,self%nw !doping
       self%n(iw,itemp) = carriers(self%vvdos_mesh,self%edos%dos(:,ispin)*max_occ,kT,self%vvdos_mesh(iw)) / &
                          cryst%ucvol / (Bohr_Ang * 1.0d-10)**3
     end do
     ! Compute carrier density at n0
     call polyn_interp(self%vvdos_mesh,self%n(:,itemp),self%ebands%fermie,n0,n0_dy)
     self%n(:,itemp) = self%n(:,itemp) - n0

     ! Compute mobility
     do ii=1,3
       do jj=1,3
         do iw=1,self%nw
           if (abs(self%n(iw,itemp)) > tol12) then
             self%mobility(iw,ispin,ii,jj,itemp) = self%sigma(iw,ispin,ii,jj,itemp) / &
                                                   ( e_Cb * self%n(iw,itemp) ) * 100**2
           end if
         end do
       end do
     end do
   end do
 end do

 contains
 real(dp) function carriers(wmesh,dos,kT,mu)

 !Arguments -------------------------------------------
 real(dp),intent(in) :: kT, mu
 real(dp),intent(in) :: wmesh(self%nw)
 real(dp),intent(in) :: dos(self%nw)

 !Local variables -------------------------------------
 integer :: ispin, iw, imu
 real(dp) :: kernel(self%nw)
 real(dp) :: integral(self%nw)

 do iw=1,self%nw
   kernel(iw) = dos(iw) * occ_fd(wmesh(iw),kT,mu)
 end do
 call simpson_int(self%nw,edos_step,kernel,integral)
 carriers = integral(self%nw)

 end function carriers

 subroutine onsager(vvdos,order,lorder)

 !Arguments -------------------------------------------
 integer,intent(in) :: order
 real(dp),intent(in) ::  vvdos(self%nw,self%nsppol,3,3,self%ntemp+1)
 real(dp),intent(out) :: lorder(self%nw,self%nsppol,3,3,self%ntemp+1)

 !Local variables -------------------------------------
 integer :: ispin, iw, imu
 real(dp) :: spin_degen, fact, mu, ee, kT
 real(dp) :: kernel(self%nw,self%nsppol,3,3)
 real(dp) :: integral(self%nw)

 ! Get spin degeneracy
 spin_degen = 2
 if (self%nsppol == 2) spin_degen = 1

 fact = spin_degen / ( 2*Ha_s )

 do itemp=1,self%ntemp
   kT = self%kTmesh(itemp)
   do imu=1,self%nw
     mu = self%vvdos_mesh(imu)
     do iw=1,self%nw
       ee = self%vvdos_mesh(iw)
       kernel(iw,:,:,:) = fact * vvdos(iw,:,:,:,1+itemp) * (mu - ee)**order * occ_dfd(ee,kT,mu)
     end do
     do ispin=1,self%nsppol
       do ii=1,3
         do jj=1,3
           call simpson_int(self%nw,edos_step,kernel(:,ispin,ii,jj),integral)
           lorder(imu,ispin,ii,jj,itemp) = integral(self%nw)
         end do
       end do
     end do
   end do
 end do

 end subroutine onsager

end subroutine transport_rta_compute
!!***

!----------------------------------------------------------------------

!!****f* m_transport/transport_rta_ncwrite
!! NAME
!! transport_rta_ncwrite
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

#ifdef HAVE_NETCDF
 ! Write to netcdf file
 ncerr = nctk_def_dims(ncid, [ nctkdim_t("ntemp", self%ntemp), &
                               nctkdim_t("ndop", self%ndop), &
                               nctkdim_t("nsppol", self%nsppol) ], defmode=.True.)
 NCF_CHECK(ncerr)
 NCF_CHECK(self%edos%ncwrite(ncid))
 NCF_CHECK(ebands_ncwrite(self%ebands,ncid))
 NCF_CHECK(cryst%ncwrite(ncid))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('vvdos_mesh', "dp", "edos_nw")], defmode=.True.))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('kTmesh', "dp", "ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('vvdos_vals', "dp", "edos_nw, nsppol_plus1, three, three")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('vvdos_tau', "dp", "edos_nw, nsppol_plus1, three, three, ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('mu_carriers', "dp", "ntemp, ndop")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('carriers', "dp", "ndop")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('N',  "dp", "edos_nw, ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('L0', "dp", "edos_nw, nsppol, three, three, ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('L1', "dp", "edos_nw, nsppol, three, three, ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('L2', "dp", "edos_nw, nsppol, three, three, ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('sigma',   "dp", "edos_nw, nsppol, three, three, ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('mobility',"dp", "edos_nw, nsppol, three, three, ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('kappa',   "dp", "edos_nw, nsppol, three, three, ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('seebeck', "dp", "edos_nw, nsppol, three, three, ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('pi',      "dp", "edos_nw, nsppol, three, three, ntemp")]))
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kTmesh"), self%kTmesh))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "mu_carriers"), self%mu_carriers))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "carriers"),    self%carriers))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vvdos_mesh"), self%vvdos_mesh))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vvdos_vals"), self%vvdos(:,1,:,:,:,1)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vvdos_tau"),  self%vvdos(:,1,:,:,:,2:)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "N"), self%n))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "L0"), self%l0(:,:,:,:,:self%ntemp)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "L1"), self%l1(:,:,:,:,:self%ntemp)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "L2"), self%l2(:,:,:,:,:self%ntemp)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "sigma"),   self%sigma(:,:,:,:,:self%ntemp)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "mobility"),self%mobility(:,:,:,:,:self%ntemp)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kappa"),   self%kappa(:,:,:,:,:self%ntemp)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "seebeck"), self%seebeck(:,:,:,:,:self%ntemp)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "pi"),      self%pi(:,:,:,:,:self%ntemp)))

#endif

end subroutine transport_rta_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_transport/transport_rta_write
!! NAME
!! transport_rta_write
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE

subroutine transport_rta_write(self)

!Arguments --------------------------------------
 type(transport_rta_t),intent(in) :: self

!Local variables --------------------------------
 integer :: ii, jj, itemp, ispin, idop
 real(dp) :: sigma, sigma_dy

 ! Write output to the terminal for the test suite
 write(std_out,*) "Electric conductivity [S/ms]"
 do itemp=1,self%ntemp
   write(std_out,"(a,f8.2,a)") "T= ", self%kTmesh(itemp) / kb_HaK, " [K]"
   do idop=1,self%ndop
   do ispin=1,self%nsppol
     do ii=1,3
       do jj=1,3
         call polyn_interp(self%vvdos_mesh,self%mobility(:,ispin,jj,ii,itemp),&
                           self%mu_carriers(itemp,idop),sigma,sigma_dy)
         write(std_out,"(f18.8)",advance="no") sigma
       end do
       write(std_out,*)
     end do
     write(std_out,*)
   end do !spin
   end do !dop
 end do !temp

end subroutine transport_rta_write
!!***

!----------------------------------------------------------------------

!!****f* m_transport/transport_rta_free
!! NAME
!! transport_rta_free
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE

subroutine transport_rta_free(self)

!Arguments --------------------------------------
 type(transport_rta_t),intent(inout) :: self

 call ebands_free(self%ebands)
 call self%gaps%free()
 call self%edos%free()

 ! free the allocated arrays and datastructure
 ABI_SFREE(self%n)

 ABI_SFREE(self%vvdos)
 ABI_SFREE(self%vvdos_mesh)
 ABI_SFREE(self%kTmesh)
 ABI_SFREE(self%eminmax_spin)
 ABI_SFREE(self%carriers)
 ABI_SFREE(self%mu_carriers)
 ABI_SFREE(self%velocity)
 ABI_SFREE(self%linewidth_mrta)
 ABI_SFREE(self%linewidth_serta)

 ABI_SFREE(self%l0)
 ABI_SFREE(self%l1)
 ABI_SFREE(self%l2)

 ABI_SFREE(self%sigma)
 ABI_SFREE(self%mobility)
 ABI_SFREE(self%seebeck)
 ABI_SFREE(self%kappa)
 ABI_SFREE(self%pi)

end subroutine transport_rta_free
!!***

end module m_transport
!!***
