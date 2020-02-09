!!****m* ABINIT/m_outscfcv
!! NAME
!!  m_outscfcv
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2020 ABINIT group (XG)
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

module m_outscfcv

 use defs_basis
 use defs_wvltypes
 use m_abicore
 use m_sort
 use m_efield
 use m_errors
 use m_xmpi
 use m_mpinfo
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_nctk
 use m_hdr
 use m_plowannier
 use m_splines
 use m_ebands
 use m_dtset
 use m_dtfil

 use defs_datatypes,     only : pseudopotential_type, ebands_t
 use defs_abitypes,      only : MPI_type
 use m_time,             only : timab
 use m_io_tools,         only : open_file
 use m_fstrings,         only : strcat, endswith
 use m_geometry,         only : bonds_lgth_angles
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype
 use m_oper,             only : oper_type,init_oper,destroy_oper
 use m_crystal,          only : crystal_init, crystal_t, prt_cif
 use m_results_gs,       only : results_gs_type, results_gs_ncwrite
 use m_ioarr,            only : ioarr, fftdatar_write
 use m_nucprop,          only : calc_efg,calc_fc
 use m_outwant,          only : outwant
 use m_pawang,           only : pawang_type
 use m_pawrad,           only : pawrad_type, simp_gen, bound_deriv
 use m_pawtab,           only : pawtab_type
 use m_paw_an,           only : paw_an_type
 use m_paw_ij,           only : paw_ij_type
 use m_paw_mkrho,        only : denfgr
 use m_pawfgrtab,        only : pawfgrtab_type
 use m_pawrhoij,         only : pawrhoij_type, pawrhoij_nullify, pawrhoij_copy, pawrhoij_free
 use m_pawcprj,          only : pawcprj_type
 use m_pawfgr,           only : pawfgr_type
 use m_paw_dmft,         only : paw_dmft_type,init_dmft,destroy_dmft,print_dmft
 use m_paw_optics,       only : optics_paw,optics_paw_core
 use m_paw_tools,        only : pawprt
 use m_numeric_tools,    only : simpson_int
 use m_epjdos,           only : dos_calcnwrite, partial_dos_fractions, partial_dos_fractions_paw, &
                                epjdos_t, epjdos_new, prtfatbands, fatbands_ncwrite
 use m_paral_atom,       only : get_my_atmtab, free_my_atmtab
 use m_io_kss,           only : outkss
 use m_multipoles,       only : multipoles_out, out1dm
 use m_mlwfovlp_qp,      only : mlwfovlp_qp
 use m_paw_mkaewf,       only : pawmkaewf
 use m_dens,             only : mag_penalty_e, calcdenmagsph, prtdenmagsph
 use m_mlwfovlp,         only : mlwfovlp
 use m_datafordmft,      only : datafordmft
 use m_mkrho,            only : read_atomden
 use m_positron,         only : poslifetime, posdoppler
 use m_optics_vloc,      only : optics_vloc
 use m_green,            only : green_type,compute_green,&
                                fourier_green,print_green,init_green,destroy_green,init_green_tau
 use m_self,             only : self_type,initialize_self,rw_self,destroy_self,destroy_self,selfreal2imag_self

 implicit none

 private
!!***

 public :: outscfcv
!!***

contains
!!***

!!****f* ABINIT/outscfcv
!! NAME
!! outscfcv
!!
!! FUNCTION
!! Output routine for the scfcv.F90 routine
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  cg(2,mcg)=planewave coefficients of wavefunctions (see also side effects)
!!  compch_fft=compensation charge, from FFT grid
!!  compch_sph=compensation charge, from sphere
!!  cprj(natom,mcprj*usecprj)=<p_lmn|Cnk> coefficients for each WF |Cnk>
!!          and each |p_lmn> non-local projector. See also side effects
!!  dimcprj(natom*usecprj)=array of dimensions of array cprj (not ordered)
!!  dmatpawu= fixed occupation matrix of correlated orbitals (DFT+U or DMFT only)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  elfr(nfft,nspden(+1))=electron localization function, real space.
!!   (+1) if spin-polarized in order to get total, spin up and spin down elf
!!  etotal=total energy
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  grhor(nfft,nspden,3)= gradient of electron density in electrons/bohr**4, real space
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  intgres(nspden,natom)=integrated residuals from constrained DFT. They are also Lagrange parameters, or gradients with respect to constraints.
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  lrhor(nfft,nspden)= Laplacian of electron density in electrons/bohr**5, real space
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*my_nspinor*mband*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mgfftc=maximum size of 1D FFTs for the PAW coarse grid
!!  mkmem =number of k points treated by this node.
!!  mpi_enreg=information about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw.
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nfft=(effective) number of FFT grid points (for this processor) (see NOTES at beginning of scfcv)
!!  ngfft(18)=contain all needed information about 3D FFT (see NOTES at beginning of scfcv)
!!  nhat(nfft,nspden*usepaw)= compensation charge density  (PAW)
!!  nkpt=number of k points.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of symmetries in space group
!!  ntypat=number of types of atoms in unit cell.
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  occ(mband*nkpt*nsppol)=occupation number for each band (usually 2) for each k.
!!  paw_an(my_natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr(natom) <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(my_natom) <type(pawfgrtab_type)> tables on PAW fine grid
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  paw_ij(my_natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!            note:structure factors are given on the coarse grid for PAW
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!     forces and its components, the stress tensor) of a ground-state computation
!!  rhor(nfft,nspden)=total electron density in electrons/bohr**3, real space.
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  taur(nfft,nspden)=total kinetic energy density in bohr**(-5), real space.
!!  ucvol=unit cell volume (bohr**3)
!!  usecprj=1 if cprj datastructure has been allocated
!!  vhartr(nfft)=Hartree potential
!!  vxc(nfft,nspden)=xc potential
!!  vtrial(nfft,nspden)=the trial potential = vxc + vpsp + vhartr, roughly speaking
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  (only writing, printing)
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  If prtwant==3 the following quantitities are updated using the unitary transformation
!!  defining the QP amplitudes in terms of the KS basis set:
!!   cg(2,mcg)=planewave coefficients of wavefunctions.
!!   cprj(natom,mcprj*usecpyj)=<p_lmn|Cnk> coefficients for each WF |Cnk> and each |p_lmn> non-local projector
!!
!! NOTES
!!   The string passed to fftdatar_write (first argument) gives the name used to store the data in the netcdf file
!!   The function  varname_from_fname defined in the module m_hdr.F90 gives the mapping between the Abinit
!!   file extension and the netcdf name e.g. foo_VHXC.nc --> vxc
!!   This function is used in cut3d so that we can immediately select the data to analyze without having
!!   to prompt the user. Remember to update varname_from_fname if you add a new file or if you change the
!!   name of the variable.
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      bonds_lgth_angles,bound_deriv,calc_efg,calc_fc,calcdenmagsph
!!      compute_coeff_plowannier,crystal_free,crystal_init,datafordmft,denfgr
!!      destroy_dmft,destroy_oper,destroy_plowannier,dos_calcnwrite,ebands_free
!!      ebands_init,ebands_interpolate_kpath,ebands_prtbltztrp,ebands_write
!!      fatbands_ncwrite,fftdatar_write,free_my_atmtab
!!      get_my_atmtab,init_dmft,init_oper,init_plowannier,ioarr,mag_penalty_e
!!      mlwfovlp,mlwfovlp_qp,multipoles_out,optics_paw,optics_paw_core
!!      optics_vloc,out1dm,outkss,outwant,partial_dos_fractions
!!      partial_dos_fractions_paw,pawmkaewf,pawprt,pawrhoij_copy
!!      pawrhoij_nullify,posdoppler,poslifetime,print_dmft,prt_cif,prtfatbands
!!      read_atomden,simpson_int,sort_dp,spline,splint,timab,wrtout,xmpi_sum
!!      xmpi_sum_master
!!
!! SOURCE

subroutine outscfcv(atindx1,cg,compch_fft,compch_sph,cprj,dimcprj,dmatpawu,dtfil,dtset,&
& ecut,eigen,electronpositron,elfr,etotal,gmet,gprimd,grhor,hdr,intgres,kg,&
& lrhor,mband,mcg,mcprj,mgfftc,mkmem,mpi_enreg,mpsang,mpw,my_natom,natom,&
& nattyp,nfft,ngfft,nhat,nkpt,npwarr,nspden,nsppol,nsym,ntypat,n3xccc,occ,&
& paw_dmft,pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,paw_an,paw_ij,&
& prtvol,psps,results_gs,rhor,rprimd,&
& taur,ucvol,usecprj,vhartr,vpsp,vtrial,vxc,wvl_den,xccc3d,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcg,mcprj,mgfftc,mkmem,mpsang,mpw,n3xccc,my_natom,natom,nfft
 integer,intent(in) :: nkpt,nspden,nsppol,nsym,ntypat,prtvol,usecprj
 real(dp),intent(in) :: compch_fft,compch_sph,ecut,ucvol
 real(dp),intent(inout) :: etotal
 type(electronpositron_type),pointer :: electronpositron
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(paw_dmft_type), intent(inout)  :: paw_dmft
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(results_gs_type),intent(in) :: results_gs
 type(wvl_denspot_type), intent(in) :: wvl_den
!arrays
 integer,intent(in) :: atindx1(natom),dimcprj(natom*usecprj)
 integer,intent(in) :: kg(3,mpw*mkmem),nattyp(ntypat),ngfft(18),npwarr(nkpt)
 real(dp),intent(in) :: dmatpawu(:,:,:,:),eigen(mband*nkpt*nsppol)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: intgres(:,:) ! (nspden,natom) if constrainedDFT otherwise (nspden,0)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol)
 real(dp),intent(in) :: rprimd(3,3),vhartr(nfft),xccc3d(n3xccc)
 real(dp),intent(in) :: vpsp(nfft)
 real(dp),intent(inout) :: cg(2,mcg)
 real(dp),intent(inout) :: nhat(nfft,nspden*psps%usepaw)
 real(dp),intent(inout),target :: rhor(nfft,nspden),vtrial(nfft,nspden)
 real(dp),intent(inout) :: vxc(nfft,nspden),xred(3,natom)
 real(dp),pointer :: elfr(:,:),grhor(:,:,:),lrhor(:,:),taur(:,:)
 type(pawcprj_type),intent(inout) :: cprj(natom,mcprj*usecprj)
 type(paw_an_type),intent(inout) :: paw_an(my_natom*psps%usepaw)
 type(pawfgrtab_type),intent(in) :: pawfgrtab(my_natom*psps%usepaw)
 type(paw_ij_type),intent(inout) :: paw_ij(my_natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawrhoij_type),target,intent(inout) :: pawrhoij(my_natom*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,cplex1=1,fform_den=52,rdwr2=2,rdwrpaw0=0
 integer :: bantot,fform,collect,timrev
 integer :: accessfil,coordn
 integer :: ii,ierr,ifft,ikpt,ispden,isppol,itypat
 integer :: me_fft,n1,n2,n3
 integer :: ifgd, iatom, iatom_tot,nradint
 integer :: me,my_natom_tmp
 integer :: occopt
 integer :: prtnabla
 integer :: pawprtden
 integer :: iband,nocc,spacecomm,comm_fft,tmp_unt,nfft_tot
 integer :: my_comm_atom
 integer :: opt_imagonly
#ifdef HAVE_NETCDF
 integer :: ncid
#endif
 real(dp) :: norm,occ_norm,unocc_norm
 real(dp) :: rate_dum,rate_dum2
 real(dp) :: yp1, ypn, dr
 character(len=500) :: message
 character(len=fnlen) :: fname
!arrays
 integer, allocatable :: isort(:)
 integer, pointer :: my_atmtab(:)
 real(dp) :: tsec(2),nt_ntone_norm(nspden),rhomag(2,nspden)
 real(dp),allocatable :: eigen2(:)
 real(dp),allocatable :: elfr_down(:,:),elfr_up(:,:),intgden(:,:)
 real(dp),allocatable :: rhor_paw(:,:),rhor_paw_core(:,:),rhor_paw_val(:,:),vpaw(:,:),vwork(:,:)
 real(dp),allocatable :: rhor_n_one(:,:),rhor_nt_one(:,:),ps_norms(:,:,:)
 real(dp), allocatable :: doccde(:)
 real(dp), allocatable :: vh1spl(:)
 real(dp), allocatable :: vh1_interp(:)
 real(dp), allocatable :: vh1_integ(:)
 real(dp), allocatable :: vh1_corrector(:)
 real(dp), allocatable :: radii(:)
 real(dp), ABI_CONTIGUOUS pointer :: rho_ptr(:,:)
 type(pawrhoij_type) :: pawrhoij_dum(0)
 type(pawrhoij_type),pointer :: pawrhoij_all(:)
 logical :: remove_inv
 logical :: paral_atom, paral_fft, my_atmtab_allocated
 real(dp) :: e_fermie
 type(oper_type) :: lda_occup
 type(crystal_t) :: crystal
 type(ebands_t) :: ebands
 type(epjdos_t) :: dos
 type(plowannier_type) :: wan
 type(self_type) :: selfr
 type(self_type) :: self
 type(green_type) :: greenr

! *************************************************************************

 DBG_ENTER("COLL")
 call timab(950,1,tsec) ! outscfcv

 if ((usecprj==0.or.mcprj==0).and.psps%usepaw==1.and. &
& (dtset%prtwant==2.or.dtset%prtwant==3.or.dtset%prtnabla>0.or.dtset%prtdos==3 &
& .or.dtset%kssform==3.or.dtset%pawfatbnd>0.or.dtset%pawprtwf>0)) then
   write (message,'(5a)')&
&   'cprj datastructure must be allocated',ch10,&
&   'with options prtwant=2,3, prtnabla>0, prtdos>3, kssform==3, pawfatbnd>0, pawprtwf>0',ch10,&
&   'Action: change pawusecp input keyword.'
   MSG_ERROR(message)
 end if

!Initialize two objects to facilitate the propagation of info.
!These objects should used more frequently, actually they should
!become basic objects used in abinit.

!Crystalline structure.
 remove_inv=.false.
 if(dtset%nspden==4 .and. dtset%usedmft==1) remove_inv=.true. ! MG: why this?

 timrev = 2; if (any(dtset%kptopt == [3, 4])) timrev= 1
 call crystal_init(dtset%amu_orig(:,1),crystal,dtset%spgroup,natom,dtset%npsp,ntypat, &
& dtset%nsym,rprimd,dtset%typat,xred,dtset%ziontypat,dtset%znucl,timrev,&
& dtset%nspden==2.and.dtset%nsppol==1,remove_inv,hdr%title,&
& dtset%symrel,dtset%tnons,dtset%symafm)

!Electron band energies.
 bantot= dtset%mband*dtset%nkpt*dtset%nsppol
 ABI_MALLOC(doccde,(bantot))
 doccde=zero

 call ebands_init(bantot,ebands,dtset%nelect,doccde,eigen,hdr%istwfk,hdr%kptns,hdr%nband,&
& hdr%nkpt,hdr%npwarr,hdr%nsppol,hdr%nspinor,hdr%tphysel,hdr%tsmear,hdr%occopt,hdr%occ,hdr%wtk,&
& hdr%charge, hdr%kptopt, hdr%kptrlatt_orig, hdr%nshiftk_orig, hdr%shiftk_orig, &
& hdr%kptrlatt, hdr%nshiftk, hdr%shiftk)

 ABI_FREE(doccde)

 ebands%fermie  = results_gs%energies%e_fermie
 e_fermie = results_gs%energies%e_fermie
 ebands%entropy = results_gs%energies%entropy
 !write(std_out,*)"ebands%efermi in outscfcv",ebands%fermie
 !write(std_out,*)"results_gs%energies%e_fermie in outscfcv",e_fermie
 !write(std_out,*)"results_gs%fermie in outscfcv",results_gs%fermie
 !write(std_out,*)"hdr%efermi in outscfcv",hdr%fermie

 ! Parameters for MPI-FFT
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3); nfft_tot = product(ngfft(1:3))
 comm_fft = mpi_enreg%comm_fft
 me_fft = xmpi_comm_rank(comm_fft)
 paral_fft = (mpi_enreg%paral_kgb==1)

 spacecomm = mpi_enreg%comm_cell
 me = xmpi_comm_rank(spacecomm)

 paral_atom=(my_natom/=natom)
 my_comm_atom = mpi_enreg%comm_atom
 nullify(my_atmtab)
 if (paral_atom) then
   call get_my_atmtab(mpi_enreg%comm_atom, my_atmtab, my_atmtab_allocated, paral_atom,natom,my_natom_ref=my_natom)
 else
   ABI_ALLOCATE(my_atmtab, (natom))
   my_atmtab = (/ (iatom, iatom=1, natom) /)
   my_atmtab_allocated = .true.
 end if

 ! YAML output
 if (me == master) then
  call results_gs%yaml_write(ab_out, dtset, crystal, comment="Summary of ground state results")
 end if

!wannier interface
 call timab(951,1,tsec)
 if (dtset%prtwant==2) then

   call mlwfovlp(crystal, ebands, hdr, atindx1,cg,cprj,dtset,dtfil,eigen,gprimd,kg,&
&   mband,mcg,mcprj,mgfftc,mkmem,mpi_enreg,mpw,natom,&
&   nattyp,nfft,ngfft,nkpt,npwarr,nsppol,ntypat,occ,&
&   pawang,pawrad,pawtab,prtvol,psps,rprimd,ucvol,xred)

 else if (dtset%prtwant==3) then

!  Convert cg and eigen to GW quasiparticle wave functions and eigenvalues in mlwfovlp_qp
   ABI_ALLOCATE(eigen2,(mband*nkpt*nsppol))
   eigen2=eigen

   call mlwfovlp_qp(cg,cprj,dtset,dtfil,eigen2,mband,mcg,mcprj,mkmem,mpw,natom,&
&   nkpt,npwarr,nspden,nsppol,ntypat,Hdr,pawtab,rprimd,MPI_enreg)

!  Call Wannier90
   call mlwfovlp(crystal, ebands, hdr, atindx1,cg,cprj,dtset,dtfil,eigen2,gprimd,kg,&
&   mband,mcg,mcprj,mgfftc,mkmem,mpi_enreg,mpw,natom,&
&   nattyp,nfft,ngfft,nkpt,npwarr,nsppol,ntypat,occ,&
&   pawang,pawrad,pawtab,prtvol,psps,rprimd,ucvol,xred)

!  this is the old implementation, risky due to unpredictable size effects
!  now eigen is not overwritten, one should use other ways to print the GW corrections
!  eigen=eigen2
   ABI_DEALLOCATE(eigen2)
 end if !prtwant
 call timab(951,2,tsec)

 occopt=dtset%occopt

 prtnabla=dtset%prtnabla
 pawprtden=dtset%prtden-1

 call timab(952,1,tsec)

 spacecomm=mpi_enreg%comm_cell; me=xmpi_comm_rank(spacecomm)
 comm_fft=mpi_enreg%comm_fft
 paral_atom=(my_natom/=natom)

!Warnings :
!- core charge is excluded from the charge density;
!- the potential is the INPUT vtrial.

 if (iwrite_fftdatar(mpi_enreg) .and. dtset%usewvl==0) then

   ! output the density.
   if (dtset%prtden/=0) then
     if (dtset%positron/=1) rho_ptr => rhor
     if (dtset%positron==1) rho_ptr => electronpositron%rhor_ep
     call fftdatar_write("density",dtfil%fnameabo_app_den,dtset%iomode,hdr,&
     crystal,ngfft,cplex1,nfft,nspden,rho_ptr,mpi_enreg,ebands=ebands)

     if (dtset%positron/=0) then
       if (dtset%positron/=1) rho_ptr => electronpositron%rhor_ep
       if (dtset%positron==1) rho_ptr => rhor
       fname = trim(dtfil%fnameabo_app_den)//'_POSITRON'
       if (dtset%iomode == IO_MODE_ETSF) fname = strcat(fname, ".nc")
       call fftdatar_write("positron_density",fname,dtset%iomode,hdr,&
       crystal,ngfft,cplex1,nfft,nspden,rho_ptr,mpi_enreg,ebands=ebands)
     end if
   end if

 else if (dtset%usewvl == 1 .and. dtset%prtden /= 0) then
   !if iomode == 2 then set all outputs to netcdf format
   !if iomode == 3 then set all outputs to ETSF format
   accessfil = 0
   if (dtset%iomode == IO_MODE_ETSF) accessfil = 3
   if (dtset%iomode == IO_MODE_MPI) accessfil = 4
   fform = fform_den
    ! Write wavelet DEN. Note however that this should be delegate to separated Bigdft routines.
    ! a lot of stuff written in outscf does not make sense if usewvl==0
   call ioarr(accessfil,rhor,dtset,etotal,fform,dtfil%fnameabo_app_den, &
   hdr,mpi_enreg,ngfft,cplex1,nfft,pawrhoij_dum,rdwr2,rdwrpaw0,wvl_den)
 end if ! if master

!! MS - Printing of PAWDEN parallellised and several possible options included
!We output the total electron density in the PAW case
!this requires removing nhat from rhor and making PAW on-site corrections
 if (pawprtden>0 .and. psps%usepaw==1) then
!  pawprtden 1 --> output PAW valence density
!  "     2 --> output PAW valence+core density
!  "     3 --> output core, valence and full atomic protodensity
!  "     4 --> options 1+3
!  "     5 --> options 2+3
!  "     6 --> output all individual PAW density contributions
   if (pawprtden/=3) then ! calc PAW valence density
     ABI_ALLOCATE(rhor_paw,(pawfgr%nfft,nspden))
     ABI_ALLOCATE(rhor_n_one,(pawfgr%nfft,nspden))
     ABI_ALLOCATE(rhor_nt_one,(pawfgr%nfft,nspden))
!    If the communicator used for denfgr is kpt_comm, it is not compatible with paral_atom
     if (mpi_enreg%paral_kgb==0.and.my_natom/=natom) then
       my_natom_tmp=natom
       ABI_DATATYPE_ALLOCATE(pawrhoij_all,(natom))
       call pawrhoij_nullify(pawrhoij_all)
       call pawrhoij_copy(pawrhoij,pawrhoij_all,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&       keep_cplex=.false.,keep_qphase=.false.,keep_itypat=.false.,keep_nspden=.false.)
     else
       my_natom_tmp=my_natom
       pawrhoij_all => pawrhoij
     end if
     if (pawprtden/=6) then
       call denfgr(atindx1,gmet,comm_fft,my_natom_tmp,natom,nattyp,ngfft,nhat,dtset%nspinor,nsppol,nspden,&
&       ntypat,pawfgr,pawrad,pawrhoij_all,pawtab,prtvol,rhor,rhor_paw,rhor_n_one,&
&       rhor_nt_one,rprimd,dtset%typat,ucvol,xred,&
&       comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
     else
       call denfgr(atindx1,gmet,comm_fft,my_natom_tmp,natom,nattyp,ngfft,nhat,dtset%nspinor,nsppol,nspden,&
&       ntypat,pawfgr,pawrad,pawrhoij_all,pawtab,prtvol,rhor,rhor_paw,rhor_n_one,&
&       rhor_nt_one,rprimd,dtset%typat,ucvol,xred,&
&       abs_n_tilde_nt_diff=nt_ntone_norm,znucl=dtset%znucl,&
&       comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
     end if
     if (mpi_enreg%paral_kgb==0.and.my_natom/=natom) then
       call pawrhoij_free(pawrhoij_all)
       ABI_DATATYPE_DEALLOCATE(pawrhoij_all)
     end if

     if (prtvol>9) then  ! Check normalisation
       norm = SUM(rhor_paw(:,1))*ucvol/PRODUCT(pawfgr%ngfft(1:3))
       call xmpi_sum(norm,comm_fft,ierr)
       write(message,'(a,F8.4)') '  PAWDEN - NORM OF DENSITY: ',norm
       call wrtout(std_out,message,'COLL')
     end if
   end if

   if (pawprtden>1.AND.pawprtden<6) then ! We will need the core density
     ABI_ALLOCATE(rhor_paw_core,(pawfgr%nfft,nspden))
     call read_atomden(mpi_enreg,natom,pawfgr%nfft,pawfgr%ngfft,nspden,ntypat,rhor_paw_core,&
&     dtset%typat,rprimd,xred,prtvol,file_prefix='core   ')

     if (prtvol>9) then  ! Check normalisation
       norm = SUM(rhor_paw_core(:,1))*ucvol/PRODUCT(pawfgr%ngfft(1:3))
       call xmpi_sum(norm,comm_fft,ierr)
       write(message,'(a,F8.4)') '  ATMDEN - NORM OF CORE DENSITY: ', norm
       call wrtout(std_out,message,'COLL')
     end if
   end if

   if (pawprtden>2.AND.pawprtden<6) then ! We will need the valence protodensity
     ABI_ALLOCATE(rhor_paw_val,(pawfgr%nfft,nspden))
     call read_atomden(mpi_enreg,natom,pawfgr%nfft,pawfgr%ngfft,nspden,ntypat,rhor_paw_val,&
&     dtset%typat,rprimd,xred,prtvol,file_prefix='valence')

     if (prtvol>9) then ! Check normalisation
       norm = SUM(rhor_paw_val(:,1))*ucvol/PRODUCT(pawfgr%ngfft(1:3))
       call xmpi_sum(norm,comm_fft,ierr)
       write(message,'(a,F8.4)') '  ATMDEN - NORM OF VALENCE PROTODENSITY: ', norm
       call wrtout(std_out,message,'COLL')
     end if
   end if

   if (iwrite_fftdatar(mpi_enreg)) then
     if (pawprtden/=3) then
       if (pawprtden==2.or.pawprtden==5) rhor_paw = rhor_paw + rhor_paw_core
!      PAWDEN
       call fftdatar_write("pawrhor",dtfil%fnameabo_app_pawden,dtset%iomode,hdr,&
       crystal,ngfft,cplex1,nfft,nspden,rhor_paw,mpi_enreg,ebands=ebands)
     end if

     if (pawprtden>2.AND.pawprtden<6) then
       ! ATMDEN_CORE
       call fftdatar_write("pawrhor_core",dtfil%fnameabo_app_atmden_core,dtset%iomode,hdr,&
       crystal,ngfft,cplex1,nfft,nspden,rhor_paw_core,mpi_enreg,ebands=ebands)

       ! valence protodensity. ATMDEN_VAL
       call fftdatar_write("pawrhor_val",dtfil%fnameabo_app_atmden_val,dtset%iomode,hdr,&
       crystal,ngfft,cplex1,nfft,nspden,rhor_paw_val,mpi_enreg,ebands=ebands)

       ! full protodensity. ATMDEN_FULL
       rhor_paw_val = rhor_paw_val + rhor_paw_core
       call fftdatar_write("pawrhor_full",dtfil%fnameabo_app_atmden_full,dtset%iomode,hdr,&
       crystal,ngfft,cplex1,nfft,nspden,rhor_paw_val,mpi_enreg,ebands=ebands)
     end if

     if (pawprtden==6) then ! Print all individual contributions to the density
       ! N_TILDE - N_HAT
       ! Use rhor_paw_val as temporary array
       if (.not.allocated(rhor_paw_val))  then
         ABI_ALLOCATE(rhor_paw_val,(pawfgr%nfft,nspden))
       end if
       rhor_paw_val = rhor - nhat

       call fftdatar_write("pawrhor_ntilde_minus_nhat",dtfil%fnameabo_app_n_tilde,dtset%iomode,hdr,&
       crystal,ngfft,cplex1,nfft,nspden,rhor_paw_val,mpi_enreg,ebands=ebands)

!      N_ONSITE
       call fftdatar_write("pawrhor_n_one",dtfil%fnameabo_app_n_one,dtset%iomode,hdr,&
       crystal,ngfft,cplex1,nfft,nspden,rhor_n_one,mpi_enreg,ebands=ebands)

!      N_TILDE_ONSITE
       call fftdatar_write("pawrhor_nt_one",dtfil%fnameabo_app_nt_one,dtset%iomode,hdr,&
       crystal,ngfft,cplex1,nfft,nspden,rhor_nt_one,mpi_enreg,ebands=ebands)

     end if ! All indivdual density cont.
   end if ! if master

   ABI_SFREE(rhor_paw)
   ABI_SFREE(rhor_paw_core)
   ABI_SFREE(rhor_paw_val)
   ABI_SFREE(rhor_n_one)
   ABI_SFREE(rhor_nt_one)

 end if ! if paw+pawprtden

 call timab(952,2,tsec)
 call timab(953,1,tsec)

 ! Output of the GSR file (except when we are inside mover)
#ifdef HAVE_NETCDF
 if (me == master .and. dtset%prtgsr == 1 .and. dtset%usewvl == 0) then
   !.and. (dtset%ionmov /= 0 .or. dtset%optcell /= 0)) then
   fname = strcat(dtfil%filnam_ds(4), "_GSR.nc")

   ! Write crystal and band structure energies.
   NCF_CHECK(nctk_open_create(ncid, fname, xmpi_comm_self))
   NCF_CHECK(hdr%ncwrite(ncid, fform_den, nc_define=.True.))
   NCF_CHECK(crystal%ncwrite(ncid))
   NCF_CHECK(ebands_ncwrite(ebands, ncid))
   ! Add energy, forces, stresses
   NCF_CHECK(results_gs_ncwrite(results_gs, ncid, dtset%ecut, dtset%pawecutdg))
   NCF_CHECK(nf90_close(ncid))
 end if
#endif

 ! Output of VCLMB file
 ! The PAW correction has to be computed here (all processors contribute)
 if (psps%usepaw > 0 .AND. dtset%prtvclmb>0) then
   nradint = 1000 ! radial integration grid density
   ABI_ALLOCATE(vpaw,(nfft,nspden))
   vpaw(:,:)=zero
   if (me == master .and. my_natom > 0) then
     if (paw_an(1)%cplex > 1) then
       MSG_WARNING('cplex = 2 : complex hartree potential in PAW spheres. This is not coded yet. Imag part ignored')
     end if
   end if

   do ispden=1,nspden
     ! for points inside spheres, replace with full AE hartree potential.
     ! In principle the correction could be more subtle (not spherical)
     do iatom=1,my_natom
       iatom_tot=iatom;if (paral_atom) iatom_tot=mpi_enreg%my_atmtab(iatom)
       itypat=dtset%typat(iatom_tot)

       ABI_ALLOCATE(vh1spl,(paw_an(iatom)%mesh_size))
       ABI_ALLOCATE(vh1_corrector,(paw_an(iatom)%mesh_size))
       ABI_ALLOCATE(vh1_interp,(pawfgrtab(iatom)%nfgd))
       ABI_ALLOCATE(radii,(pawfgrtab(iatom)%nfgd))
       ABI_ALLOCATE(isort,(pawfgrtab(iatom)%nfgd))
       ! vh1 vht1 contain the spherical first moments of the Hartree potentials, so re-divide by Y_00 = sqrt(four_pi)
       vh1_corrector(:) = (paw_an(iatom)%vh1(:,1,ispden)-paw_an(iatom)%vht1(:,1,ispden)) / sqrt(four_pi)

       ! get end point derivatives
       call bound_deriv(vh1_corrector, pawrad(itypat), pawrad(itypat)%mesh_size, yp1, ypn)
       ! spline the vh1 function
       ! NB for second argument of vh1: only first moment lm_size appears to be used
       ! NB2: vh1 can in principle be complex - not sure what to do with the imaginary part. Ignored for now.
       call spline(pawrad(itypat)%rad, vh1_corrector, paw_an(iatom)%mesh_size, yp1, ypn, vh1spl)

       do ifgd = 1, pawfgrtab(iatom)%nfgd
         ! get radii for this point
         isort(ifgd) = ifgd
         radii(ifgd) = sqrt(sum(pawfgrtab(iatom)%rfgd(:,ifgd)**2))
       end do

       if (pawfgrtab(iatom)%nfgd/=0) then
       ! spline interpolate the vh1 value for current radii
         call sort_dp(pawfgrtab(iatom)%nfgd, radii, isort, tol12)
         call splint(pawrad(itypat)%mesh_size, pawrad(itypat)%rad, &
&         vh1_corrector, vh1spl, pawfgrtab(iatom)%nfgd, radii,  vh1_interp, ierr)
       end if

       norm=SUM(vh1_interp)*ucvol/PRODUCT(ngfft(1:3))
       call xmpi_sum(norm,comm_fft,ierr)
       write(message,'(a,i6,a,E20.10)') ' sum of Hartree correction term on fft grid of atom : ', iatom, &
&       ' = ', norm
       call wrtout(std_out,message,'COLL')

       if (pawfgrtab(iatom)%nfgd/=0) then
         vpaw(pawfgrtab(iatom)%ifftsph(isort(1:pawfgrtab(iatom)%nfgd)),ispden) = &
&         vpaw(pawfgrtab(iatom)%ifftsph(isort(1:pawfgrtab(iatom)%nfgd)),ispden) + &
&         vh1_interp(1:pawfgrtab(iatom)%nfgd)
       end if

       ! get integral of correction term in whole sphere
       ABI_DEALLOCATE(radii)
       ABI_DEALLOCATE(vh1_interp)

       ABI_ALLOCATE(radii,(nradint))
       ABI_ALLOCATE(vh1_interp,(nradint))

       ABI_ALLOCATE(vh1_integ,(nradint))
       dr = pawrad(itypat)%rad(paw_an(iatom)%mesh_size) / dble(nradint)
       do ifgd = 1, nradint
         radii(ifgd) = dble(ifgd-1)*dr
       end do

       ! spline interpolate the vh1 value for current radii
       call splint(pawrad(itypat)%mesh_size, pawrad(itypat)%rad, &
&       vh1_corrector, vh1spl, nradint, radii,  vh1_interp, ierr)

       do ifgd = 1, nradint
         vh1_interp(ifgd) = vh1_interp(ifgd)*radii(ifgd)**2
       end do

       call simpson_int(nradint, dr, vh1_interp, vh1_integ)
       write(message,'(a,i6,a,E20.10)') ' integral of Hartree correction term in sphere of atom: ', iatom, &
&       ' = ', vh1_integ(nradint)*four*pi
       call wrtout(std_out,message,'COLL')

       ABI_DEALLOCATE(vh1spl)
       ABI_DEALLOCATE(vh1_corrector)
       ABI_DEALLOCATE(vh1_interp)
       ABI_DEALLOCATE(vh1_integ)
       ABI_DEALLOCATE(radii)
       ABI_DEALLOCATE(isort)
     end do ! iatom
   end do !ispden
   call xmpi_sum_master(vpaw,master,mpi_enreg%comm_atom,ierr)
   if (.not.iwrite_fftdatar(mpi_enreg)) then
     ABI_DEALLOCATE(vpaw)
   end if
 end if ! if paw - add all electron vhartree in spheres

 if (iwrite_fftdatar(mpi_enreg)) then

   ! output the electron localization function ELF
   if (dtset%prtelf/=0) then
     call fftdatar_write("elfr",dtfil%fnameabo_app_elf,dtset%iomode,hdr,&
     crystal,ngfft,cplex1,nfft,nspden,elfr,mpi_enreg,ebands=ebands)

     if (nspden==2)then
       ABI_ALLOCATE(elfr_up,(nfft,nspden))
       elfr_up(:,:) = zero
       do ifft=1,nfft
         elfr_up(ifft,1) = elfr(ifft,2)
       end do
!      ELF_UP
       call fftdatar_write("elfr_up",dtfil%fnameabo_app_elf_up,dtset%iomode,hdr,&
       crystal,ngfft,cplex1,nfft,nspden,elfr_up,mpi_enreg,ebands=ebands)

       ABI_ALLOCATE(elfr_down,(nfft,nspden))
       elfr_down(:,:) = zero
       do ifft=1,nfft
         elfr_down(ifft,1) = elfr(ifft,3)
       end do
!      ELF_DOWN'
       call fftdatar_write("elfr_down",dtfil%fnameabo_app_elf_down,dtset%iomode,hdr,&
       crystal,ngfft,cplex1,nfft,nspden,elfr_down,mpi_enreg,ebands=ebands)

       ABI_DEALLOCATE(elfr_up)
       ABI_DEALLOCATE(elfr_down)
     end if
   end if

   call timab(953,2,tsec)
   call timab(954,1,tsec)

!  We output the gradient of density
   if (dtset%prtgden/=0) then

     call fftdatar_write("grhor_1",dtfil%fnameabo_app_gden1,dtset%iomode,hdr,&
     crystal,ngfft,cplex1,nfft,nspden,grhor(:,:,1),mpi_enreg,ebands=ebands)

     call fftdatar_write("grhor_2",dtfil%fnameabo_app_gden2,dtset%iomode,hdr,&
     crystal,ngfft,cplex1,nfft,nspden,grhor(:,:,2),mpi_enreg,ebands=ebands)

     call fftdatar_write("grhor_3",dtfil%fnameabo_app_gden3,dtset%iomode,hdr,&
     crystal,ngfft,cplex1,nfft,nspden,grhor(:,:,3),mpi_enreg,ebands=ebands)
   end if

!  We output the total kinetic energy density KDEN
   if (dtset%prtkden/=0) then
     call fftdatar_write("kinedr",dtfil%fnameabo_app_kden,dtset%iomode,hdr,&
     crystal,ngfft,cplex1,nfft,nspden,taur,mpi_enreg,ebands=ebands)
   end if

!  We output the Laplacian of density
   if (dtset%prtlden/=0) then
     call fftdatar_write("laprhor",dtfil%fnameabo_app_lden,dtset%iomode,hdr,&
     crystal,ngfft,cplex1,nfft,nspden,lrhor,mpi_enreg,ebands=ebands)
   end if

   call timab(954,2,tsec)
   call timab(955,1,tsec)

   call timab(955,2,tsec)
   call timab(956,1,tsec)

!  POT
   if (dtset%prtpot>0) then
     call fftdatar_write("vtrial",dtfil%fnameabo_app_pot,dtset%iomode,hdr,&
     crystal,ngfft,cplex1,nfft,nspden,vtrial,mpi_enreg,ebands=ebands)
   end if

   call timab(956,2,tsec)
   call timab(957,1,tsec)

   if (dtset%prtgeo>0) then
     coordn=dtset%prtgeo
     call bonds_lgth_angles(coordn,dtfil%fnameabo_app_geo,natom,psps%ntypat,&
      rprimd,dtset%typat,xred,dtset%znucl)
   end if

   if (dtset%prtcif > 0) then
     call prt_cif(dtset%brvltt, dtfil%fnameabo_app_cif, natom, dtset%nsym, dtset%ntypat, rprimd, &
      dtset%spgaxor, dtset%spgroup, dtset%spgorig, dtset%symrel, dtset%tnons, dtset%typat, xred, dtset%znucl)
   end if

   call timab(957,2,tsec)
   call timab(958,1,tsec)

!  STM
   if (dtset%prtstm>0) then
     call fftdatar_write("stm",dtfil%fnameabo_app_stm,dtset%iomode,hdr,&
     crystal,ngfft,cplex1,nfft,nspden,rhor,mpi_enreg,ebands=ebands)
   end if

   if (dtset%prt1dm>0) then
     call out1dm(dtfil%fnameabo_app_1dm,mpi_enreg,natom,nfft,ngfft,nspden,psps%ntypat,&
&     rhor,rprimd,dtset%typat,ucvol,vtrial,xred,dtset%znucl)
   end if

!  VHA
   if (dtset%prtvha>0) then
     ABI_ALLOCATE(vwork,(nfft,nspden))
     do ispden=1,nspden
       vwork(:,ispden)=vhartr(:)
     end do

     call fftdatar_write("vhartree",dtfil%fnameabo_app_vha,dtset%iomode,hdr,&
     crystal,ngfft,cplex1,nfft,nspden,vwork,mpi_enreg,ebands=ebands)

     ABI_DEALLOCATE(vwork)
   end if

!  VPSP
   if (dtset%prtvpsp>0) then
     ABI_ALLOCATE(vwork,(nfft,nspden))
     do ispden=1,nspden
       vwork(:,ispden)=vpsp(:)
     end do

     call fftdatar_write("vpsp",dtfil%fnameabo_app_vpsp,dtset%iomode,hdr,&
     crystal,ngfft,cplex1,nfft,nspden,vwork,mpi_enreg,ebands=ebands)

     ABI_DEALLOCATE(vwork)
   end if

! VCouLoMB
   if (dtset%prtvclmb>0) then

     ABI_ALLOCATE(vwork,(nfft,nspden))
     do ispden=1,nspden
       vwork(:,ispden)=vpsp(:)+vhartr(:)
     end do
     if (psps%usepaw==1) then
       do ispden=1,nspden
         vwork(:,ispden)=vwork(:,ispden)+vpaw(:,ispden)
       end do
       ABI_DEALLOCATE(vpaw)
     end if

     call fftdatar_write("vhartree_vloc",dtfil%fnameabo_app_vclmb,dtset%iomode,hdr,&
&     crystal,ngfft,cplex1,nfft,nspden,vwork,mpi_enreg,ebands=ebands)

!TODO: find out why this combination of calls with fftdatar_write then out1dm fails on buda with 4 mpi-fft procs (npkpt 1).
!      For the moment comment it out. Only DS2 of mpiio test 27 fails
!     call out1dm(dtfil%fnameabo_app_vclmb_1dm,mpi_enreg,natom,nfft,ngfft,nspden,psps%ntypat,&
!&         rhor,rprimd,dtset%typat,ucvol,vwork,xred,dtset%znucl)

! TODO: add TEM phase with CE = (2 pi / lambda) (E+E0)/(E(E+2E0)) from p.49 of RE Dunin Borkowski 2004 encyclopedia of nanoscience volume 3 pp 41-99
!   where E is energy of electron, E0 rest mass, lambda the relativistic wavelength
!   values of CE at 200 300 and 1000 kV:  7.29e6  6.53e6   5.39e6 rad / V / m
!   vertical integral of vclmb * c / ngfft(3) / cross sectional area factor (= sin(gamma))
!      * 0.5291772083e-10*27.2113834 to get to SI
!      * CE factor above
!   should be done for each plane perpendicular to the axes...
     ABI_DEALLOCATE(vwork)
   end if ! prtvclmb


!  VHXC
   if (dtset%prtvhxc>0) then
     ABI_ALLOCATE(vwork,(nfft,nspden))
     do ispden=1,nspden
       vwork(:,ispden)=vhartr(:)+vxc(:,ispden)
     end do

     call fftdatar_write("vhxc",dtfil%fnameabo_app_vhxc,dtset%iomode,hdr,&
&     crystal,ngfft,cplex1,nfft,nspden,vwork,mpi_enreg,ebands=ebands)

     ABI_DEALLOCATE(vwork)
   end if

!  VXC
   if (dtset%prtvxc>0) then
     call fftdatar_write("exchange_correlation_potential",dtfil%fnameabo_app_vxc,dtset%iomode,hdr,&
     crystal,ngfft,cplex1,nfft,nspden,vxc,mpi_enreg,ebands=ebands)
   end if

   call timab(958,2,tsec)
 end if ! if iwrite_fftdatar

 call timab(959,1,tsec)

!Generate DOS using the tetrahedron method or using Gaussians
!FIXME: Should centralize all calculations of DOS here in outscfcv
 if (dtset%prtdos>=2.or.dtset%pawfatbnd>0) then
   dos = epjdos_new(dtset, psps, pawtab)

   if (dos%partial_dos_flag>=1 .or. dos%fatbands_flag==1)then
     ! Generate fractions for partial DOSs if needed partial_dos 1,2,3,4  give different decompositions
     collect = 1 !; if (psps%usepaw==1 .and. dos%partial_dos_flag /= 2) collect = 0
     if ((psps%usepaw==0.or.dtset%pawprtdos/=2) .and. dos%partial_dos_flag>=1) then
       call partial_dos_fractions(dos,crystal,dtset,eigen,occ,npwarr,kg,cg,mcg,collect,mpi_enreg)
     end if

     if (psps%usepaw==1 .and. dos%partial_dos_flag /= 2) then
!      TODO: update partial_dos_fractions_paw for extra atoms - no PAW contribution normally, but check bounds and so on.
       call partial_dos_fractions_paw(dos,cprj,dimcprj,dtset,mcprj,mkmem,mpi_enreg,pawrad,pawtab)
     end if

   else
     dos%fractions(:,:,:,1)=one
   end if

!  Here, print out fatbands for the k-points given in file appended _FATBANDS
   if (me == master .and. dtset%pawfatbnd>0 .and. dos%fatbands_flag==1) then
     call prtfatbands(dos,dtset,ebands,dtfil%fnameabo_app_fatbands,dtset%pawfatbnd,pawtab)
   end if

!  Here, computation and output of DOS and partial DOS  _DOS
   if (dos%fatbands_flag == 0 .and. dos%prtdos /= 4) then
     call dos_calcnwrite(dos,dtset,crystal,ebands,dtfil%fnameabo_app_dos,spacecomm)
   end if

#ifdef HAVE_NETCDF
   ! Write netcdf file with dos% results.
   if (me == master) then
     fname = trim(dtfil%filnam_ds(4))//'_FATBANDS.nc'
     NCF_CHECK(nctk_open_create(ncid, fname, xmpi_comm_self))
     call fatbands_ncwrite(dos, crystal, ebands, hdr, dtset, psps, pawtab, ncid)
     NCF_CHECK(nf90_close(ncid))
   end if
#endif

   call dos%free()
 end if ! prtdos > 1

 call timab(959,2,tsec)
 call timab(960,1,tsec)

!Output of integrated density inside atomic spheres
 if (dtset%prtdensph==1.and.dtset%usewvl==0)then
   ABI_MALLOC(intgden,(nspden,natom))
   call calcdenmagsph(gmet,mpi_enreg,natom,nfft,ngfft,nspden,&
&   ntypat,dtset%ratsm,dtset%ratsph,rhor,rprimd,dtset%typat,ucvol,xred,1,cplex1,intgden=intgden,rhomag=rhomag)
   if(all(dtset%constraint_kind(:)==0))then
     call prtdenmagsph(cplex1,intgden,natom,nspden,ntypat,ab_out,1,dtset%ratsm,dtset%ratsph,rhomag,dtset%typat)
   else
     call prtdenmagsph(cplex1,intgden,natom,nspden,ntypat,ab_out,1,dtset%ratsm,dtset%ratsph,rhomag,dtset%typat,dtset%ziontypat)
   endif
   if(any(dtset%constraint_kind(:)/=0))then
     call prtdenmagsph(cplex1,intgres,natom,nspden,ntypat,ab_out,11,dtset%ratsm,dtset%ratsph,rhomag,dtset%typat)
   endif
   ABI_SFREE(intgden)
 end if

 call timab(960,2,tsec)

 if (dtset%magconon /= 0) then
!  calculate final value of terms for magnetic constraint: "energy" term, lagrange multiplier term, and atomic contributions
   call mag_penalty_e(dtset%magconon,dtset%magcon_lambda,mpi_enreg,&
&   natom,nfft,ngfft,nspden,ntypat,dtset%ratsm,dtset%ratsph,rhor,rprimd,dtset%spinat,dtset%typat,xred)
 end if

 call timab(961,1,tsec)

!If PAW, provide additional outputs
 if (psps%usepaw==1) then
!  Output of compensation charge
   if (dtset%nstep>0.or.dtfil%ireadwf/=0) then
     write(message, '(4a)' )ch10,' PAW TEST:',ch10,&
&     ' ==== Compensation charge inside spheres ============'
     if (compch_sph>-1.d4.and.compch_fft>-1.d4) &
&     write(message, '(3a)' ) trim(message),ch10,' The following values must be close to each other ...'
     if (compch_sph>-1.d4) write(message, '(3a,f22.15)' ) trim(message),ch10,&
&     ' Compensation charge over spherical meshes = ',compch_sph
     if (compch_fft>-1.d4) then
       if (pawfgr%usefinegrid==1) then
         write(message, '(3a,f22.15)' ) trim(message),ch10,&
&         ' Compensation charge over fine fft grid    = ',compch_fft
       else
         write(message, '(3a,f22.15)' ) trim(message),ch10,&
&         ' Compensation charge over fft grid         = ',compch_fft
       end if
     end if
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if
!  Output of pseudopotential strength Dij and augmentation occupancies Rhoij
   call pawprt(dtset,my_natom,paw_ij,pawrhoij,pawtab,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&   electronpositron=electronpositron)
 end if

 call timab(961,2,tsec)
 call timab(962,1,tsec)


!PAW + output for optical conductivity   _OPT and _OPT2
 if (psps%usepaw==1.and.prtnabla>0) then
   if (prtnabla==1.or.prtnabla==2) then
     call optics_paw(atindx1,cg,cprj,dimcprj,dtfil,dtset,eigen,gprimd,hdr,kg,&
&     mband,mcg,mcprj,mkmem,mpi_enreg,mpsang,mpw,natom,nkpt,npwarr,nsppol,pawrad,pawtab)
   end if
   if (prtnabla==2.or.prtnabla==3) then
     call optics_paw_core(atindx1,cprj,dimcprj,dtfil,dtset,eigen,psps%filpsp,hdr,&
&     mband,mcprj,mkmem,mpi_enreg,mpsang,natom,nkpt,nsppol,pawrad,pawtab)
   end if
 end if
 if (prtnabla<0) then
   ! TODO: This routine is not tested but it's used in production.
   call optics_vloc(cg,dtfil,dtset,eigen,gprimd,hdr,kg,&
&   mband,mcg,mkmem,mpi_enreg,mpw,nkpt,npwarr,nsppol)
 end if

 call timab(962,2,tsec)
 call timab(963,1,tsec)

!Optionally provide output for AE wavefunctions (only for PAW)
 if (psps%usepaw==1 .and. dtset%pawprtwf==1) then
   ABI_ALLOCATE(ps_norms,(nsppol,nkpt,mband))

   call pawmkaewf(dtset,crystal,ebands,my_natom,mpw,mband,mcg,mcprj,nkpt,mkmem,nsppol,Dtset%nband,&
&   Dtset%istwfk,npwarr,Dtset%kptns,Dtset%ngfftdg,kg,dimcprj,pawfgrtab,&
&   Pawrad,Pawtab,Hdr,Dtfil,cg,Cprj,&
&   MPI_enreg,ierr,pseudo_norms=ps_norms,set_k=dtset%pawprt_k,set_band=dtset%pawprt_b,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

   if (dtset%pawprt_b==0) then
     fname = strcat(dtfil%filnam_ds(4), '_PAWSTAT')
     if (open_file(fname, message,newunit=tmp_unt,status='unknown',form='formatted') /= 0) then
       MSG_ERROR(message)
     end if
     write(tmp_unt,'(5a)') '# This file contains the statistics on the cancellation of',ch10,&
&     '# the onsite pseudo component of the all-electron wavefunction',ch10,&
&     '# with the plane wave part'
     ii = 0
     do isppol=1,nsppol
       write(tmp_unt,'(a,i0)') '# isppol = ',isppol
       do ikpt=1,nkpt
         write(tmp_unt,'(a,i0)') '# ikpt = ',ikpt
         write(tmp_unt,'(a)') '#    band      norm'
         occ_norm = zero; unocc_norm = zero; nocc = 0
         do iband=1,dtset%nband(ikpt + (isppol-1)*nkpt)
           ii = ii + 1
           write(tmp_unt,'(i8,ES16.6)') iband,ps_norms(isppol,ikpt,iband)
           if (abs(occ(ii)) <= tol16) then
             unocc_norm = unocc_norm + ps_norms(isppol,ikpt,iband)
           else
             occ_norm = occ_norm + ps_norms(isppol,ikpt,iband)
             nocc = nocc + 1
           end if
         end do
         if(mband/=nocc)then
           write(tmp_unt,'(2(a,ES16.6))') '# occ average: ',occ_norm/real(nocc),&
&           ' unocc average: ',unocc_norm/real(mband-nocc)
         else
           write(tmp_unt,'(2(a,ES16.6))') '# occ average: ',occ_norm/real(nocc)
         end if
       end do
     end do
     close(tmp_unt)
   end if
   ABI_DEALLOCATE(ps_norms)
 end if

 call timab(963,2,tsec)
 if(dtset%plowan_compute>0) then
   write(message,'(2a,i3)') ch10,&
&   ' ====================================================================================== '
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   write(message,'(2a,i3)') ch10,&
&   ' == Start computation of Projected Local Orbitals Wannier functions == ',dtset%nbandkss
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

!  ==  compute psichi

   call init_plowannier(dtset,wan)
   call compute_coeff_plowannier(crystal,cprj,dimcprj,dtset,eigen,e_fermie,&
&   mpi_enreg,occ,wan,pawtab,psps,usecprj,dtfil%unpaw,pawrad,dtfil)
   call destroy_plowannier(wan)
 end if

!Optionally provide output for the GW part of ABINIT
 if (dtset%nbandkss/=0) then
   ! Use DMFT to compute wannier function for cRPA calculation.
   if(dtset%usedmft==1) then
     write(message,'(2a,i3)') ch10,&
&     '  Warning: Psichi are renormalized in datafordmft because nbandkss is used',dtset%nbandkss
     call wrtout(std_out,message,'COLL')
     call init_dmft(dmatpawu,dtset,e_fermie,dtfil%fnameabo_app,&
&     dtfil%filnam_ds(3),dtset%nspinor,paw_dmft,pawtab,psps,dtset%typat)
     call print_dmft(paw_dmft,dtset%pawprtvol)

!    ==  compute psichi
     call init_oper(paw_dmft,lda_occup)

     call datafordmft(crystal,cprj,dimcprj,dtset,eigen,e_fermie &
&     ,lda_occup,dtset%mband,dtset%mband,dtset%mkmem,mpi_enreg,&
&     dtset%nkpt,dtset%nspinor,dtset%nsppol,occ,&
&     paw_dmft,paw_ij,pawang,pawtab,psps,usecprj,dtfil%unpaw,dtset%nbandkss)

     opt_imagonly=0
     if(paw_dmft%dmft_solv>=5) opt_imagonly=1


     ! Compute k-resolved spectral function in DMFT.
     if(dtset%dmft_kspectralfunc==1) then
      ! Initialize self on real axis
       call initialize_self(selfr,paw_dmft,wtype='real')

      ! Initialize self on  imag axis
       call initialize_self(self,paw_dmft)

      ! Initialize green on real axis
       call init_green(greenr,paw_dmft,opt_oper_ksloc=3,wtype='real')

      ! Read self energy in imag. Matsubara freq (for double counting
      ! and limit at high frequency)
       call rw_self(self,paw_dmft,prtopt=5,opt_rw=1,opt_stop=1)

      ! Read self energy on real axis obtained from Maxent
       call rw_self(selfr,paw_dmft,prtopt=5,opt_rw=1,opt_imagonly=opt_imagonly, &
     & opt_selflimit=self%oper(self%nw)%matlu,opt_hdc=self%hdc%matlu,pawang=pawang,cryst_struc=crystal)

      ! Check: from self on real axis, recompute self on Imaginary axis.
       call selfreal2imag_self(selfr,self,paw_dmft%filapp)

      !  paw_dmft%fermie=hdr%fermie ! for tests
       write(std_out,*) "    Fermi level is",paw_dmft%fermie

       ! selfr does not have any double couting in self%hdc
       ! hdc from self%hdc has been put in real part of self in rw_self.
       ! For the LDA BS: use opt_self=0 and fermie=fermie_lda

      ! Compute green  function on real axis
       call compute_green(crystal,greenr,paw_dmft,pawang,1,selfr,&
&       opt_self=1,opt_nonxsum=0)

      !write(6,*) "compute green done"
       if(me==master) then
         call print_green("from_realaxisself",greenr,5,paw_dmft,&
&         pawprtvol=3,opt_wt=1)
        !write(6,*) "print green done"
       endif

       call destroy_green(greenr)
       call destroy_self(selfr)
       call destroy_self(self)
     endif
     call destroy_dmft(paw_dmft)
     call destroy_oper(lda_occup)
   end if

   call timab(964,1,tsec) ! outscfcv(outkss)

   call outkss(crystal,dtfil,dtset,ecut,gmet,gprimd,hdr,&
&   dtset%kssform,mband,mcg,mcprj,mgfftc,mkmem,mpi_enreg,mpsang,mpw,natom,natom,&
&   nfft,nkpt,npwarr,nspden,nsppol,nsym,psps%ntypat,occ,pawtab,pawfgr,paw_ij,&
&   prtvol,psps,rprimd,vtrial,xred,cg,usecprj,cprj,eigen,ierr)
   call timab(964,2,tsec) ! outscfcv(outkss)
   if (ierr/=0) then
     MSG_WARNING("outkss returned a non zero status error, check log")
   end if
 end if

 if (electronpositron_calctype(electronpositron)/=0) then

!  Optionally provide output for  positron life time calculation
   call timab(965,1,tsec)
   call poslifetime(dtset,electronpositron,gprimd,my_natom,&
&   mpi_enreg,n3xccc,nfft,ngfft,nhat,1,pawang,&
&   pawrad,pawrhoij,pawtab,rate_dum,rate_dum2,&
&   rhor,ucvol,xccc3d)
   call timab(965,2,tsec)

!  Optionally provide output for momentum distribution of annihilation radiation
   if (dtset%posdoppler>0) then
     call posdoppler(cg,cprj,crystal,dimcprj,dtfil,dtset,electronpositron,psps%filpsp,&
&     kg,mcg,mcprj,mpi_enreg,my_natom,n3xccc,nfft,ngfft,nhat,npwarr,&
&     occ,pawang,pawrad,pawrhoij,pawtab,rhor,xccc3d)
   end if
 end if

!Optionally provide output for WanT
 if (dtset%prtwant==1) then
   call timab(966,1,tsec)
   ! WARNING: mpi_enreg not used --> MPI is not supported
   call outwant(dtset,eigen,cg,kg,npwarr,mband,mcg,nkpt,nsppol,mkmem,mpw,dtset%prtwant)
   call timab(966,2,tsec)
 end if

!Optionally provide output for electric field gradient calculation
 if (dtset%prtefg > 0) then
   call timab(967,1,tsec)
   call calc_efg(mpi_enreg,my_natom,natom,nfft,ngfft,nspden,dtset%nsym,ntypat,&
&   paw_an,pawang,pawrad,pawrhoij,pawtab,&
&   dtset%ptcharge,dtset%prtefg,dtset%quadmom,rhor,rprimd,dtset%symrel,&
&   dtset%tnons,dtset%typat,ucvol,psps%usepaw,xred,psps%zionpsp,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   call timab(967,2,tsec)
 end if

!Optionally provide output for Fermi-contact term at nuclear positions
 if (dtset%prtfc > 0) then
   call timab(967,1,tsec)
   call calc_fc(my_natom,natom,nspden,ntypat,pawrad,pawrhoij,pawtab,dtset%typat,psps%usepaw,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   call timab(967,2,tsec)
 end if

 ! Output electron bands.
 if (me == master .and. dtset%tfkinfunc==0) then
   if (size(dtset%kptbounds, dim=2) > 0) then
     call ebands_write(ebands, dtset%prtebands, dtfil%filnam_ds(4), kptbounds=dtset%kptbounds)
   else
     call ebands_write(ebands, dtset%prtebands, dtfil%filnam_ds(4))
   end if
 end if

!Optionally provide Xcrysden output for the Fermi surface (Only master writes)
 if (me == master .and. dtset%prtfsurf == 1) then
   if (ebands_write_bxsf(ebands,crystal,dtfil%fnameabo_app_bxsf) /= 0) then
     message = "Cannot produce BXSF file with Fermi surface, see log file for more info"
     MSG_WARNING(message)
     call wrtout(ab_out, message)
   end if
 end if ! prtfsurf==1

!output nesting factor for Fermi surface (requires ph_nqpath)
 if (me == master .and. dtset%prtnest>0 .and. dtset%ph_nqpath > 0) then
   call timab(968,1,tsec)
   ierr = ebands_write_nesting(ebands,crystal,dtfil%fnameabo_app_nesting,dtset%prtnest,&
&   dtset%tsmear,dtset%fermie_nest,dtset%ph_qpath(:,1:dtset%ph_nqpath),message)
   if (ierr /=0) then
     MSG_WARNING(message)
     call wrtout(ab_out,message,'COLL')
   end if
   call timab(968,2,tsec)
 end if ! prtnest=1

 call timab(969,1,tsec)

 if (dtset%prtdipole == 1) then
   call multipoles_out(rhor,mpi_enreg,natom,nfft,ngfft,dtset%nspden,dtset%ntypat,rprimd,&
&   dtset%typat,ucvol,ab_out,xred,dtset%ziontypat)
 end if

 ! BoltzTraP output files in GENEric format
 if (dtset%prtbltztrp == 1 .and. me==master) then
   call ebands_prtbltztrp(ebands, crystal, dtfil%filnam_ds(4))
 end if

 ! Band structure interpolation from eigenvalues computed on the k-mesh.
 if (nint(dtset%einterp(1)) /= 0) then
   call ebands_interpolate_kpath(ebands, dtset, crystal, [0, 0], dtfil%filnam_ds(4), spacecomm)
 end if

 if(associated(elfr))then
   ABI_DEALLOCATE(elfr)
 end if

 if(associated(grhor))then
   ABI_DEALLOCATE(grhor)
 end if

 if(associated(lrhor))then
   ABI_DEALLOCATE(lrhor)
 end if

 call crystal%free()
 call ebands_free(ebands)

 ! Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 call timab(969,2,tsec)
 call timab(950,2,tsec) ! outscfcv

 DBG_EXIT("COLL")

end subroutine outscfcv
!!***

end module m_outscfcv
!!***
