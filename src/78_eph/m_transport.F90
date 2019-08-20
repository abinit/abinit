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

 use defs_datatypes,   only : ebands_t
 use m_crystal,        only : crystal_t
 use m_numeric_tools,  only : bisect, simpson_int, safe_div !polyn_interp,
 use m_fstrings,       only : strcat, sjoin, ltoa
 use m_kpts,           only : listkk
 use m_occ,            only : occ_fd, occ_dfd

 implicit none

 private
!!****

 public :: transport ! main entry point for transport calculations
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
   ! number of spinorial components

   integer :: ntemp
   ! number of temperatures

   integer :: ndop
   ! number of carrier concentrarions at which to evaluate chemical potential energy

   integer :: nw
   ! number of frequencies at which trasport quantities are computed

   real(dp) :: eph_extrael
   ! extra electrons per unit cell from sigeph (lifetimes)

   real(dp) :: eph_fermie
   ! fermi level from input file from sigeph (lifetimes)

   real(dp) :: transport_extrael
   ! extra electrons per unit cell

   real(dp) :: transport_fermie
   ! fermi level from input file

   real(dp),allocatable :: kTmesh(:)
   ! a list of temperatures at which to compute the transport

   real(dp),allocatable :: eph_mu_e(:)
   ! (ntemp, ndop)
   ! Chemical potential at this carrier concentrarion and temperature from sigeph (lifetime)

   real(dp),allocatable :: transport_mu_e(:)
   ! (ntemp, ndop)
   ! Chemical potential at this carrier concentrarion and temperature

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

   real(dp),allocatable :: ne(:)
   ! (ntemp) number of electrons at transport_mu_e(ntemp)

   real(dp),allocatable :: nh(:)
   ! (ntemp) number of holes at transport_mu_e(ntemp)

   real(dp),allocatable :: l0(:,:,:,:,:)
   real(dp),allocatable :: l1(:,:,:,:,:)
   real(dp),allocatable :: l2(:,:,:,:,:)
   ! (nw, nsppol, 3, 3, ntemp)
   ! onsager coeficients

   real(dp),allocatable :: sigma(:,:,:,:,:)
   real(dp),allocatable :: seebeck(:,:,:,:,:)
   real(dp),allocatable :: kappa(:,:,:,:,:)
   real(dp),allocatable :: pi(:,:,:,:,:)
   ! (nw, nsppol, 3, 3, ntemp)
   ! transport coefficients

   real(dp),allocatable :: mobility(:,:,:,:,:,:)
   ! Mobility

   real(dp),allocatable :: n(:,:,:)
   ! (nw,ntemp) carrier density (n/cm^3)

   real(dp),allocatable :: mobility_mu(:,:,:,:,:)
   ! (2, nsppol, 3, 3, ntemp)
   ! mobility for electrons and holes (first dimension) at transport_mu_e(ntemp)

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
!!
!! INPUTS
!! dtfil<datafiles_type>=variables related to files.
!! dtset<dataset_type>=All input variables for this dataset.
!! ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!! cryst<crystal_t>=Crystalline structure
!! comm=MPI communicator.
!!
!! SOURCE

subroutine transport(dtfil, dtset, ebands, cryst, comm)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
!arrays

!Local variables ------------------------------
 type(sigmaph_t) :: sigmaph
 type(transport_rta_t) :: transport_rta
 integer :: ierr, my_rank
 real(dp) :: extrael_fermie(2)
#ifdef HAVE_NETCDF
 integer :: ncid
#endif
 character(len=fnlen) :: path, sigeph_path
 character(len=500) :: msg

! *************************************************************************

 ABI_UNUSED((/comm, cryst%natom/))

 my_rank = xmpi_comm_rank(comm)
 call wrtout(std_out, ' Transport computation driver')

 sigeph_path = strcat(dtfil%filnam_ds(4), "_SIGEPH.nc")
 sigmaph = sigmaph_read(sigeph_path, dtset, xmpi_comm_self, msg, ierr, keep_open=.true., extrael_fermie=extrael_fermie)
 ABI_CHECK(ierr == 0, msg)
 ! if dtset%sigma_ngkpt /= sigeph_path

 ! Initialize transport
 transport_rta = transport_rta_new(dtset, sigmaph, cryst, ebands, extrael_fermie, comm)
 sigmaph%ncid = nctk_noid
 call sigmaph%free()

 ! Compute transport
 call transport_rta_compute(transport_rta, cryst, dtset, comm)

 ! Compute mobility
 call transport_rta_compute_mobility(transport_rta, cryst, dtset, comm)

 ! Write transport to stdout (for the test suite)
 call transport_rta_write(transport_rta, cryst)

 ! Master creates the netcdf file used to store the results of the calculation.
#ifdef HAVE_NETCDF
 if (my_rank == 0) then
   path = strcat(dtfil%filnam_ds(4), "_TRANSPORT.nc")
   NCF_CHECK(nctk_open_create(ncid, path, xmpi_comm_self))
   call transport_rta_ncwrite(transport_rta, cryst, dtset, ncid)
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
!!  dtset<dataset_type>=All input variables for this dataset.
!!  sigmaph<sigmaph_t>=Object with e-ph self-energy results.
!!  cryst<crystal_t>=Crystalline structure
!!  ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!!  extrael_fermie
!!  comm=MPI communicator.
!!
!!
!! SOURCE

type(transport_rta_t) function transport_rta_new(dtset, sigmaph, cryst, ebands, extrael_fermie, comm) result (new)

!Arguments -------------------------------------
 integer, intent(in) :: comm
 type(sigmaph_t),intent(in) :: sigmaph
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(dataset_type),intent(in) :: dtset
 real(dp),intent(in) :: extrael_fermie(2)

!Local variables ------------------------------
 type(ebands_t) :: tmp_ebands
 integer,parameter :: occopt3=3, timrev1=1, sppoldbl1=1
 integer :: ierr, itemp, spin
 integer :: nprocs, my_rank
 integer :: kptrlatt(3,3)
 integer,allocatable :: indkk(:,:)
 real(dp) :: nelect, dksqmax
 character(len=500) :: msg

!************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 ! Allocate temperature arrays
 new%ntemp = sigmaph%ntemp
 ABI_MALLOC(new%kTmesh, (new%ntemp))
 new%kTmesh = sigmaph%kTmesh

 ! Information about the Gaps
 new%nsppol = ebands%nsppol
 new%nspinor = ebands%nspinor
 ABI_MALLOC(new%eminmax_spin, (2,ebands%nsppol))
 new%eminmax_spin = get_minmax(ebands, "eig")

 ierr = get_gaps(ebands, new%gaps)
 if (ierr /= 0) then
   do spin=1, ebands%nsppol
     MSG_WARNING(trim(new%gaps%errmsg_spin(spin)))
     new%gaps%vb_max(spin) = ebands%fermie - 1 * eV_Ha
     new%gaps%cb_min(spin) = ebands%fermie + 1 * eV_Ha
   end do
   MSG_WARNING("get_gaps returned non-zero exit status. See above warning messages...")
 end if

 ! Read lifetimes to ebands object
 if (any(dtset%sigma_ngkpt /= 0)) then
   ! If integrals are computed with sigma_ngkpt k-mesh, we need to downsample ebands.
   call wrtout(std_out, sjoin(" Computing integrals with downsampled sigma_ngkpt:", ltoa(dtset%sigma_ngkpt)))
   kptrlatt = 0
   kptrlatt(1,1) = dtset%sigma_ngkpt(1)
   kptrlatt(2,2) = dtset%sigma_ngkpt(2)
   kptrlatt(3,3) = dtset%sigma_ngkpt(3)
   tmp_ebands = ebands_downsample(ebands, cryst, kptrlatt, 1, [zero,zero,zero])
   new%ebands = sigmaph%get_ebands(cryst, tmp_ebands, new%linewidth_serta, new%linewidth_mrta, &
                               new%velocity, xmpi_comm_self, ierr)
   call ebands_free(tmp_ebands)
 else
   new%ebands = sigmaph%get_ebands(cryst, ebands, new%linewidth_serta, new%linewidth_mrta, &
                                  new%velocity, xmpi_comm_self, ierr)
 end if

 ! Perform further downsampling (usefull for debugging purposes)
 ! TODO: introduce transport_ngkpt variable
 if (any(dtset%transport_ngkpt/=0)) then

   call wrtout(std_out, sjoin(" Downsampling the k-point mesh before computing transport:", ltoa(dtset%transport_ngkpt)))
   kptrlatt = 0
   kptrlatt(1,1) = dtset%transport_ngkpt(1)
   kptrlatt(2,2) = dtset%transport_ngkpt(2)
   kptrlatt(3,3) = dtset%transport_ngkpt(3)
   tmp_ebands = ebands_downsample(new%ebands, cryst, kptrlatt, 1, [zero,zero,zero])

   ! Map the points of downsampled bands to dense ebands
   ABI_MALLOC(indkk,(tmp_ebands%nkpt, 6))
   call listkk(dksqmax, cryst%gmet, indkk, new%ebands%kptns, tmp_ebands%kptns, &
               new%ebands%nkpt, tmp_ebands%nkpt, cryst%nsym, &
               sppoldbl1, cryst%symafm, cryst%symrec, timrev1, comm, exit_loop=.True., use_symrec=.True.)

   if (dksqmax > tol12) then
      write(msg, '(3a,es16.6,a)' ) &
       "Error downsampling ebands in transport driver",ch10,&
       "the k-point could not be generated from a symmetrical one. dksqmax: ",dksqmax, ch10
      MSG_ERROR(msg)
   end if

   ! linewidths serta
   if (allocated(new%linewidth_serta)) then
     call downsample_array(new%linewidth_serta,indkk,tmp_ebands%nkpt)
   end if

   ! linewidths mrta
   if (allocated(new%linewidth_mrta)) then
     call downsample_array(new%linewidth_mrta,indkk,tmp_ebands%nkpt)
   end if

   ! velocities
   if (allocated(new%linewidth_serta)) then
     call downsample_array(new%velocity,indkk,tmp_ebands%nkpt)
   end if

   ABI_SFREE(indkk)
   call ebands_copy(tmp_ebands, new%ebands)
   call ebands_free(tmp_ebands)
 end if

 ABI_CHECK(allocated(new%velocity), 'Could not read velocities from SIGEPH.nc file')

 ! Same doping case as sigmaph
 new%ndop = 1
 ABI_MALLOC(new%eph_mu_e, (new%ntemp))
 ABI_MALLOC(new%transport_mu_e, (new%ntemp))
 new%eph_extrael = extrael_fermie(1)
 new%eph_fermie = extrael_fermie(2)
 new%transport_fermie = dtset%eph_fermie
 new%transport_extrael = dtset%eph_extrael
 new%eph_mu_e = sigmaph%mu_e
 new%transport_mu_e = sigmaph%mu_e

 if (new%transport_fermie/=zero) then
   new%transport_mu_e = new%transport_fermie
 end if

 if (new%transport_fermie == zero .and. new%transport_extrael /= new%eph_extrael) then
   if (new%transport_extrael /= new%eph_extrael) then
     write(msg,'(a,e18.8,a,e18.8,a)') 'extrael from SIGEPH ',new%transport_extrael, &
                                      ' and input file ',new%eph_extrael, &
                                      ' differ. Will recompute the chemical potential'
   endif
   call wrtout(std_out, msg)
   call ebands_copy(ebands, tmp_ebands)

   ! We only need mu_e so MPI parallelize the T-loop.
   new%transport_mu_e = zero
   do itemp=1,new%ntemp
     if (mod(itemp, nprocs) /= my_rank) cycle ! MPI parallelism.
     ! Use Fermi-Dirac occopt
     call ebands_set_scheme(tmp_ebands, occopt3, new%kTmesh(itemp), dtset%spinmagntarget, dtset%prtvol)
     new%transport_mu_e(itemp) = tmp_ebands%fermie

     ! Check that the total number of electrons is correct
     ! This is to trigger problems as the routines that calculate the occupations in ebands_set_nelect
     ! are different from the occ_fd that will be used in the rest of the code.
     nelect = ebands_calc_nelect(tmp_ebands, new%kTmesh(itemp), new%transport_mu_e(itemp))

     if (abs(nelect - tmp_ebands%nelect) > tol6) then
       ! For T = 0 the number of occupied states goes in discrete steps (according to the k-point sampling)
       ! for finite doping it's hard to find nelect that exactly matches ebands%nelect.
       ! in this case we print a warning
       write(msg,'(3(a,f10.6))')&
         'Calculated number of electrons nelect = ',nelect,&
         ' does not correspond with ebands%nelect = ',tmp_ebands%nelect,' for T = ',new%kTmesh(itemp)
       if (new%kTmesh(itemp) == zero) then
         MSG_WARNING(msg)
       else
         MSG_ERROR(msg)
       end if
     end if
   end do ! it

   call ebands_free(tmp_ebands)
   call xmpi_sum(new%transport_mu_e, comm, ierr)
 endif

 contains
 subroutine downsample_array(array,indkk,nkpt)

   real(dp),allocatable,intent(inout) :: array(:,:,:,:)
   integer,allocatable,intent(in) :: indkk(:,:)

   real(dp),allocatable :: tmp_array(:,:,:,:)
   integer :: ikpt,nkpt
   integer :: tmp_shape(4)

   ABI_MOVE_ALLOC(array, tmp_array)
   tmp_shape = shape(array)
   tmp_shape(3) = nkpt
   ABI_MALLOC(array,(tmp_shape(1),tmp_shape(2),tmp_shape(3),tmp_shape(4)))
   do ikpt=1,nkpt
     array(:,:,ikpt,:) = tmp_array(:,:,indkk(ikpt,1),:)
   end do
   ABI_FREE(tmp_array)

 end subroutine downsample_array

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
!! cryst<crystal_t>=Crystalline structure
!! dtset<dataset_type>=All input variables for this dataset.
!! comm=MPI communicator.
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
 integer :: ntens, nvecs, nvals, edos_intmeth, ifermi, iel
 real(dp) :: vr(3)
 real(dp) :: emin, emax, edos_broad, edos_step, max_occ, kT
 real(dp) :: linewidth, fact0
 real(dp) :: dummy_vals(1,1,1,1), dummy_vecs(1,1,1,1,1)
 real(dp),allocatable :: vv_tens(:,:,:,:,:,:)
 real(dp),allocatable :: dummy_dosvals(:,:,:,:), dummy_dosvecs(:,:,:,:,:)
 !character(len=500) :: msg

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
       ! Store outer product in vv_tens
       do ii=1,3
         do jj=1,3
           vv_tens(ii, jj, 1, ib, ik, ispin) = vr(ii) * vr(jj)
         end do
       end do
       ! Multiply by the lifetime
       do itemp=1,self%ntemp
         linewidth = abs(self%linewidth_serta(itemp, ib, ik, ispin))
         call safe_div( vv_tens(:, :, 1, ib, ik, ispin), linewidth, zero, vv_tens(:, :, 1+itemp, ib, ik, ispin))
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

 ! Free memory
 ABI_SFREE(dummy_dosvals)
 ABI_SFREE(dummy_dosvecs)
 ABI_FREE(vv_tens)

 ! Transport coeficients computation
 ABI_MALLOC(self%l0,(self%nw,self%nsppol,3,3,self%ntemp+1))
 ABI_MALLOC(self%l1,(self%nw,self%nsppol,3,3,self%ntemp+1))
 ABI_MALLOC(self%l2,(self%nw,self%nsppol,3,3,self%ntemp+1))

 ABI_MALLOC(self%sigma,   (self%nw,self%nsppol,3,3,self%ntemp+1))
 ABI_CALLOC(self%seebeck, (self%nw,self%nsppol,3,3,self%ntemp+1))
 ABI_CALLOC(self%kappa,   (self%nw,self%nsppol,3,3,self%ntemp+1))
 ABI_MALLOC(self%pi,      (self%nw,self%nsppol,3,3,self%ntemp+1))

 ! Mobility
 ABI_MALLOC(self%n,   (self%nw,self%ntemp,2))
 ABI_MALLOC(self%mobility,(self%nw,self%nsppol,3,3,self%ntemp+1,2))

 ! Compute onsager coefficients
 call onsager(0, self%l0)
 call onsager(1, self%l1)
 call onsager(2, self%l2)

 ! Compute transport quantities
 fact0 = (Time_Sec * siemens_SI / Bohr_meter / cryst%ucvol)
 self%sigma = fact0 * self%l0
 call safe_div(volt_SI * self%l1, self%l0, zero, self%pi)
 do itemp=1,self%ntemp
   kT = self%kTmesh(itemp) / kb_HaK
   call safe_div(volt_SI * self%l1(:,:,:,:,itemp), &
                 kT * self%l0(:,:,:,:,itemp), zero, self%seebeck(:,:,:,:,itemp))

   ! HM: to write it as a single division I do:
   ! kappa = L1^2/L0 + L2 = (L1^2 + L2*L0)/L0
   ! Check why do we need minus sign here to get consistent results with Boltztrap!
   call safe_div( - volt_SI**2 * fact0 * (self%l1(:,:,:,:,itemp)**2 - self%l2(:,:,:,:,itemp)*self%l0(:,:,:,:,itemp)), &
                 kT * self%l0(:,:,:,:,itemp), zero, self%kappa(:,:,:,:,itemp))
 end do

 ! Compute the index of the fermi level
 ifermi = bisect(self%vvdos_mesh, self%ebands%fermie)

 ! Handle out of range condition.
 if (ifermi == 0 .or. ifermi == self%nw) then
   MSG_ERROR("Bisection could not find energy index of the Fermi level!")
 end if

 self%mobility = 0
 max_occ = two/(self%nspinor*self%nsppol)
 do ispin=1,self%nsppol
   do itemp=1,self%ntemp
     ! Compute carrier density
     kT = self%kTmesh(itemp)

     ! Compute carrier density of electrons
     do iw=1,self%nw !doping
       self%n(iw,itemp,1) = carriers(self%vvdos_mesh,self%edos%dos(:,ispin)*max_occ,ifermi,self%nw, &
                                     kT,self%vvdos_mesh(iw)) / &
                                     cryst%ucvol / Bohr_meter**3
     end do

     ! Compute carrier density of holes
     do iw=1,self%nw !doping
       self%n(iw,itemp,2) = carriers(self%vvdos_mesh,self%edos%dos(:,ispin)*max_occ,1,ifermi, &
                                     kT,self%vvdos_mesh(iw)) / &
                                     cryst%ucvol / Bohr_meter**3
     end do

     self%n(:,itemp,2) = self%n(self%nw,itemp,2) - self%n(:,itemp,2)

     ! Compute mobility
     do iel=1,2
     do ii=1,3
       do jj=1,3
         do iw=1,self%nw
           call safe_div( self%sigma(iw,ispin,ii,jj,itemp) * 100**2, &
                          e_Cb * self%n(iw,itemp,iel), &
                          zero, self%mobility(iw,ispin,ii,jj,itemp,iel) )
         end do
       end do
     end do
     end do
   end do !itemp
 end do !spin

contains
 real(dp) function carriers(wmesh,dos,istart,istop,kT,mu)

 !Arguments -------------------------------------------
 real(dp),intent(in) :: kT, mu
 real(dp),intent(in) :: wmesh(self%nw)
 real(dp),intent(in) :: dos(self%nw)

 !Local variables -------------------------------------
 integer :: iw, istart, istop
 real(dp) :: kernel(self%nw)
 real(dp) :: integral(self%nw)

 kernel = 0
 do iw=istart,istop
   kernel(iw) = dos(iw) * occ_fd(wmesh(iw),kT,mu)
 end do
 call simpson_int(self%nw, edos_step, kernel, integral)
 carriers = integral(self%nw)

 end function carriers

 subroutine onsager(order,lorder)

 !Arguments -------------------------------------------
 integer,intent(in) :: order
 real(dp),intent(out) :: lorder(self%nw,self%nsppol,3,3,self%ntemp+1)

 !Local variables -------------------------------------
 integer :: ispin, iw, imu
 real(dp) :: fact, mu, ee, kT
 real(dp) :: kernel(self%nw,self%nsppol,3,3)
 real(dp) :: integral(self%nw)

 ! Get spin degeneracy
 max_occ = two/(self%nspinor*self%nsppol)
 ! 2 comes from linewidth-lifetime relation
 fact = max_occ / two

 do itemp=1,self%ntemp
   kT = self%kTmesh(itemp)
   do imu=1,self%nw
     mu = self%vvdos_mesh(imu)
     do iw=1,self%nw
       ee = self%vvdos_mesh(iw)
       if (order > 0) then
         kernel(iw,:,:,:) = fact * self%vvdos(iw,1,1:,:,:,1+itemp) * &
                            (mu - ee)**order * occ_dfd(ee,kT,mu)
       else
         kernel(iw,:,:,:) = fact * self%vvdos(iw,1,1:,:,:,1+itemp) * &
                            occ_dfd(ee,kT,mu)
       end if
     end do
     do ispin=1,self%nsppol
       do ii=1,3
         do jj=1,3
           call simpson_int(self%nw, edos_step, kernel(:,ispin,ii,jj), integral)
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

!!****f* m_transport/transport_rta_compute_mobility
!! NAME
!! transport_rta_compute_mobility
!!
!! FUNCTION
!!
!! INPUTS
!! cryst<crystal_t>=Crystalline structure
!! dtset<dataset_type>=All input variables for this dataset.
!! comm=MPI communicator.
!!
!! SOURCE

subroutine transport_rta_compute_mobility(self, cryst, dtset, comm)

!Arguments ------------------------------------
 integer,intent(in) :: comm
 type(transport_rta_t),intent(inout) :: self
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst

!Local variables ------------------------------
 integer :: nsppol, nkpt, mband, ib, ik, ispin, ii, jj, itemp
 integer :: ielhol, nvalence
 integer :: bmin(2), bmax(2)
 real(dp) :: vr(3), vv_tens(3,3), vv_tenslw(3,3)
 real(dp) :: eig_nk, mu_e, linewidth, fact, fact0
 real(dp) :: max_occ, kT, wtk

 ABI_UNUSED((/dtset%natom, comm/))

 ABI_MALLOC(self%mobility_mu,(2,self%nsppol,3,3,self%ntemp+1))

 ! create alias for dimensions
 nsppol = self%ebands%nsppol
 nkpt   = self%ebands%nkpt
 mband  = self%ebands%mband

 ! Compute index of valence band
 max_occ = two/(self%nspinor*self%nsppol)
 ! TODO: should add nelect0 to ebands to keep track of intrinsic
 nvalence = nint(self%ebands%nelect - self%eph_extrael)/max_occ

 ABI_CALLOC(self%ne,(self%ntemp))
 ABI_CALLOC(self%nh,(self%ntemp))
 do ispin=1,nsppol
   do ik=1,nkpt
     wtk = self%ebands%wtk(ik)
     ! number of holes
     do ib=1,nvalence
       eig_nk = self%ebands%eig(ib, ik, ispin)
       do itemp=1,self%ntemp
         kT = self%kTmesh(itemp)
         mu_e = self%transport_mu_e(itemp)
         self%nh(itemp) = self%nh(itemp) + wtk*(1-occ_fd(eig_nk,kT,mu_e))*max_occ
       end do
     end do
     ! number of electrons
     do ib=nvalence+1,mband
       eig_nk = self%ebands%eig(ib, ik, ispin)
       do itemp=1,self%ntemp
         kT = self%kTmesh(itemp)
         mu_e = self%transport_mu_e(itemp)
         self%ne(itemp) = self%ne(itemp) + wtk*occ_fd(eig_nk,kT,mu_e)*max_occ
       end do
     end do
   end do
 end do

 ! Get units conversion factor and spin degeneracy
 fact0 = (Time_Sec * siemens_SI / Bohr_meter / cryst%ucvol)
 fact = max_occ * fact0 / e_Cb / two * 100**2

 ! Compute mobility
 self%mobility_mu = 0
 do ispin=1,nsppol
   do ik=1,nkpt
     wtk = self%ebands%wtk(ik)
     bmin = [nvalence+1,1]
     bmax = [mband, nvalence]
     do ib=1,mband
       ielhol = 2
       if (ib > nvalence) ielhol = 1
       eig_nk = self%ebands%eig(ib, ik, ispin)
       vr(:) = self%velocity(:,ib,ik,ispin)
       ! Store outer product in vv_tens
       do ii=1,3
         do jj=1,3
           vv_tens(ii, jj) = vr(ii) * vr(jj)
         end do
       end do
       ! Symmetrize tensor
       vv_tens = symmetrize_tensor(cryst,vv_tens)
       ! Multiply by the lifetime
       do itemp=1,self%ntemp
         mu_e = self%transport_mu_e(itemp)
         kT = self%kTmesh(itemp)
         linewidth = abs(self%linewidth_serta(itemp, ib, ik, ispin))
         call safe_div( wtk * vv_tens(:, :) * occ_dfd(eig_nk,kT,mu_e), linewidth, zero, vv_tenslw(:, :))
         self%mobility_mu(ielhol, ispin, :, :, itemp) = self%mobility_mu(ielhol, ispin, :, :, itemp) + vv_tenslw(:, :)
       end do
     end do
   end do !kpt
 end do !spin

 ! Scale by the carrier concentration
 do itemp=1,self%ntemp
   ! for electrons
   call safe_div( fact * self%mobility_mu(1,:,:,:,itemp), &
                  self%ne(itemp) / cryst%ucvol / Bohr_meter**3, &
                  zero, self%mobility_mu(1,:,:,:,itemp) )
   ! for holes
   call safe_div( fact * self%mobility_mu(2,:,:,:,itemp), &
                  self%nh(itemp) / cryst%ucvol / Bohr_meter**3, &
                  zero, self%mobility_mu(2,:,:,:,itemp) )
 end do

contains
 function symmetrize_tensor(cryst,t) result(tsum)
  integer :: isym
  type(crystal_t) :: cryst
  real(dp) :: tsym(3,3), tsum(3,3), t(3,3)

  !symmetrize
  tsum = 0
  do isym=1, cryst%nsym
    tsym = matmul( (cryst%symrel_cart(:,:,isym)), matmul(t, transpose(cryst%symrel_cart(:,:,isym))) )
    tsum = tsum + tsym
  end do
  tsum = tsum / cryst%nsym

 end function symmetrize_tensor

end subroutine transport_rta_compute_mobility
!!***

!----------------------------------------------------------------------

!!****f* m_transport/transport_rta_ncwrite
!! NAME
!! transport_rta_ncwrite
!!
!! FUNCTION
!!
!! INPUTS
!! cryst<crystal_t>=Crystalline structure
!! ncid=Netcdf file handle.
!!
!! SOURCE

subroutine transport_rta_ncwrite(self, cryst, dtset, ncid)

!Arguments --------------------------------------
 type(transport_rta_t),intent(in) :: self
 type(crystal_t),intent(in) :: cryst
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: ncid

!Local variables --------------------------------
 integer :: ncerr

#ifdef HAVE_NETCDF
 ! Write to netcdf file
 ncerr = nctk_def_dims(ncid, [ nctkdim_t("ntemp", self%ntemp), &
                               nctkdim_t("nsppol", self%nsppol) ], defmode=.True.)
 NCF_CHECK(ncerr)
 NCF_CHECK(self%edos%ncwrite(ncid))
 NCF_CHECK(ebands_ncwrite(self%ebands,ncid))
 NCF_CHECK(cryst%ncwrite(ncid))

 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('vvdos_mesh', "dp", "edos_nw")], defmode=.True.))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('transport_ngkpt', "int", "three")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('kTmesh', "dp", "ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('vvdos_vals', "dp", "edos_nw, nsppol_plus1, three, three")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('vvdos_tau', "dp", "edos_nw, nsppol_plus1, three, three, ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('eph_mu_e', "dp", "ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('transport_mu_e', "dp", "ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('L0', "dp", "edos_nw, nsppol, three, three, ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('L1', "dp", "edos_nw, nsppol, three, three, ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('L2', "dp", "edos_nw, nsppol, three, three, ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('sigma',   "dp", "edos_nw, nsppol, three, three, ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('kappa',   "dp", "edos_nw, nsppol, three, three, ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('seebeck', "dp", "edos_nw, nsppol, three, three, ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('pi',      "dp", "edos_nw, nsppol, three, three, ntemp")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('mobility',"dp", "edos_nw, nsppol, three, three, ntemp, two")]))

 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('N',  "dp", "edos_nw, ntemp, two")]))
 NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t('mobility_mu',"dp", "two, nsppol, three, three, ntemp")]))

 ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "eph_extrael", "eph_fermie", &
                                   "transport_extrael", "transport_fermie"])
 NCF_CHECK(ncerr)
 NCF_CHECK(nctk_set_datamode(ncid))

 ncerr = nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: "eph_extrael", "eph_fermie", &
                                     "transport_extrael", "transport_fermie"], &
                                    [self%eph_extrael, self%eph_fermie, self%transport_extrael, self%transport_fermie])
 NCF_CHECK(ncerr)

 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "transport_ngkpt"), dtset%transport_ngkpt))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kTmesh"), self%kTmesh))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "eph_mu_e"), self%eph_mu_e))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "transport_mu_e"), self%transport_mu_e))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vvdos_mesh"), self%vvdos_mesh))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vvdos_vals"), self%vvdos(:,1,:,:,:,1)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vvdos_tau"),  self%vvdos(:,1,:,:,:,2:)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "L0"), self%l0(:,:,:,:,:self%ntemp)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "L1"), self%l1(:,:,:,:,:self%ntemp)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "L2"), self%l2(:,:,:,:,:self%ntemp)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "sigma"),   self%sigma(:,:,:,:,:self%ntemp)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kappa"),   self%kappa(:,:,:,:,:self%ntemp)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "seebeck"), self%seebeck(:,:,:,:,:self%ntemp)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "pi"),      self%pi(:,:,:,:,:self%ntemp)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "N"), self%n))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "mobility"),self%mobility(:,:,:,:,:self%ntemp,:)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "mobility_mu"),self%mobility_mu(:,:,:,:,:self%ntemp)))
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
!! cryst<crystal_t>=Crystalline structure
!!
!! SOURCE

subroutine transport_rta_write(self, cryst)

!Arguments --------------------------------------
 type(transport_rta_t),intent(in) :: self
 type(crystal_t),intent(in) :: cryst

!Local variables --------------------------------
 integer :: itemp, ispin
 character(len=500) :: msg

 call wrtout(ab_out, 'Transport calculation results')
 write(msg,"(a16,a32,a32)") 'Temperature [K]', 'e/h density [cm^-3]', 'e/h mobility [cm^2/Vs]'
 call wrtout([std_out, ab_out], msg)

 do ispin=1,self%nsppol
   do itemp=1,self%ntemp
     write(msg,"(f16.2,2e16.2,2f16.2)") self%kTmesh(itemp) / kb_HaK, &
                            self%ne(itemp) / cryst%ucvol / (Bohr_meter * 100)**3, &
                            self%nh(itemp) / cryst%ucvol / (Bohr_meter * 100)**3, &
                            self%mobility_mu(1,ispin,1,1,itemp), self%mobility_mu(2,ispin,1,1,itemp)
     call wrtout([std_out, ab_out], msg)
   end do !temp
 end do !spin

end subroutine transport_rta_write
!!***

!----------------------------------------------------------------------

!!****f* m_transport/transport_rta_free
!! NAME
!! transport_rta_free
!!
!! FUNCTION
!!  Free dynamic memory.
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
 ABI_SFREE(self%eph_mu_e)
 ABI_SFREE(self%transport_mu_e)
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

 ABI_SFREE(self%mobility_mu)
 ABI_SFREE(self%nh)
 ABI_SFREE(self%ne)
 ABI_SFREE(self%n)

end subroutine transport_rta_free
!!***

end module m_transport
!!***
