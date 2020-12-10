!!****m* ABINIT/m_io_kss
!! NAME
!!  m_io_kss
!!
!! FUNCTION
!!  This module contains procedured dealing with the IO of the KSS file.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2020 ABINIT group (MG, MT, VO, AR, LR, RWG, MM, XG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

MODULE m_io_kss

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_hdr
 use m_wfk
 use m_cgtools
 use m_hamiltonian
 use m_electronpositron
 use m_pawtab
 use m_paw_ij
 use m_pawcprj
 use m_pawfgr
 use m_dtfil
 use m_dtset

 use defs_datatypes,     only : pseudopotential_type
 use defs_abitypes,      only : MPI_type
 use m_time,             only : timab
 use m_io_tools,         only : open_file
 use m_fstrings,         only : sjoin, itoa, strcat
 use m_hide_lapack,      only : xheevx_cplex, xhegvx_cplex
 use m_geometry,         only : metric, remove_inversion
 use m_mpinfo,           only : destroy_mpi_enreg, proc_distrb_cycle
 use m_fftcore,          only : get_kg, sphere
 use m_fft,              only : fftpac
 use m_crystal ,         only : crystal_t
 use m_gsphere,          only : table_gbig2kg, merge_and_sort_kg
 use m_kg,               only : mkkin, mkkpg
 use m_ksdiago,          only : ksdiago, init_ddiago_ctl, ddiago_ctl_type
 use m_mkffnl,           only : mkffnl
 use m_getghc,           only : getghc
 use m_fourier_interpol, only : transgrid

 implicit none

 private

 public :: write_kss_header    ! Writes the header of the KSS file.
 !private :: write_vkb         ! Writes the KB form factors and derivates on file for a single k-point.
 public :: write_kss_wfgk      ! Write the Gamma-centered wavefunctions and energies on the KSS file for a single k-point.
 public :: k2gamma_centered    ! Convert a set of wavefunctions from the k-centered to the gamma-centered basis set.
 public :: make_gvec_kss       ! Build the list of G-vectors for the KSS file.
 public :: outkss              ! Generate KSS file

CONTAINS  !===========================================================
!!***

!!****f* m_io_kss/write_kss_header
!! NAME
!!  write_kss_header
!!
!! FUNCTION
!!  Write the header of the KSS file either using plain Fortran-IO or netcdf with ETSF-IO format.
!!  Returns the unit number to be used for further writing.
!!  It should be executed by master node only.
!!
!! INPUTS
!!  filekss(len=fnlen)=The name of the KSS file.
!!  kss_npw=Number of planewaves used for the wavefunctions in the KSS files.
!!  ishm=Max number of shells written on file
!!  shlim(ishm)=The cumulative number of G"s in each shell.
!!  nbandksseff=Number of bands to be written.
!!  mband=The maximum number of bands treated by abinit.
!!  nsym2=Number of symmetry operations to be written on the header.
!!  symrel2(3,3,nsym2)=The symmetry operations in real space to be written.
!!  tnons2(3,nsym2)=The fractional translations associateed to symrel2.
!!  gbig(3,kss_npw)=The set of G-vectors for the KSS wavefunctions (Gamma-centered)
!!  Hdr<hdr_type>=The abinit header.
!!  Dtset <dataset_type>=all input variables for this dataset
!!  Psps<pseudopotential_type>=Structure gathering info on the pseudopotentials.
!!  iomode=Input variables specifying the fileformat. (0-->Fortran,3-->netcdf with ETSF-IO format).
!!  occ(mband*nkpt*nsppol)=The occupation factors.
!!
!! OUTPUT
!!  kss_unt=The unit number of the opened file.
!!
!! SIDE EFFECTS
!!  The KSS Header is written on file.
!!
!! PARENTS
!!      m_io_kss
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine write_kss_header(filekss,kss_npw,ishm,nbandksseff,mband,nsym2,symrel2,tnons2,occ,gbig,shlim,&
                            crystal,Dtset,Hdr,Psps,iomode,kss_unt)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iomode,kss_npw,nbandksseff,ishm,nsym2,mband
 integer,intent(out) :: kss_unt
 character(len=fnlen),intent(in) :: filekss
 type(crystal_t),intent(in) :: crystal
 type(pseudopotential_type),intent(in) :: Psps
 type(Hdr_type),intent(in) :: Hdr
 type(Dataset_type),intent(in) :: Dtset
!arrays
 integer,intent(in) :: symrel2(3,3,nsym2)
 integer,intent(in) :: gbig(3,kss_npw),shlim(ishm)
 real(dp),intent(in) :: tnons2(3,nsym2)
 real(dp),intent(in) :: occ(mband*Dtset%nkpt*Dtset%nsppol)

!Local variables-------------------------------
!scalars
 integer :: nspinor,nsppol,nkpt,itypat,ierr
 integer :: nb,isppol,ik,fform,ii,jj,kk,ig
 integer :: il,il0,ilmn,in,ind1,ind2
 character(len=80) :: title
 character(len=500) :: msg
 type(hdr_type) :: my_Hdr
 type(dataset_type) :: Dtset_cpy
#ifdef HAVE_NETCDF
 integer :: ncerr
#endif
!arrays
 integer,allocatable :: vkbsign_int(:,:,:)
 real(dp),allocatable :: vkbsign(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

 nsppol = Dtset%nsppol
 nkpt   = Dtset%nkpt
 nspinor= Dtset%nspinor

 write(msg,'(3a)')ch10,' Opening file for KS structure output: ',TRIM(filekss)
 call wrtout(std_out,msg,'COLL')

 write(msg,'(a,i6)') ' number of Gamma centered plane waves ',kss_npw
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i6)') ' number of Gamma centered shells ',ishm
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i6)') ' number of bands ',nbandksseff
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i6)') ' maximum angular momentum components ',Psps%mpsang
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i2,a)')' number of symmetry operations ',nsym2,' (without inversion)'
 call wrtout(std_out,msg,'COLL')

!Copy the header so that we can change some basic dimensions using the KSS values:
!(bantot, npwarr, nband) and the occupation factors

!Note that nsym and symrel might have been changed this has to be fixed
!carefully in the next patch since in the new implementation symmorphy=0 should be dafault
 call hdr_copy(Hdr,my_Hdr)

 my_Hdr%npwarr =kss_npw
 my_Hdr%nband  =nbandksseff
 my_hdr%mband = maxval(my_hdr%nband)
 my_Hdr%bantot =nbandksseff*nkpt*nsppol

 my_Hdr%istwfk = 1  ! KSS file does not support istwfk/=1 even though the GS run
                    ! can take advantage of time-reversal symmetry.

!Copy the occ number in the new header with correct dimensions
!fill with zero the rest since mband can be < nbandksseff
 !write(std_out,*)associated(my_Hdr%occ)
 ABI_FREE(my_Hdr%occ)
 ABI_MALLOC(my_Hdr%occ,(my_Hdr%bantot))
 !mband = MAXVAL(Hdr%nband)

 my_Hdr%occ=zero; nb=MIN(mband,nbandksseff)
 do isppol=1,nsppol
   do ik=1,nkpt
     ind1=1+(ik-1)*nbandksseff+(isppol-1)*nkpt*nbandksseff
     ind2=1+(ik-1)*mband      +(isppol-1)*nkpt*mband
     my_Hdr%occ(ind1:ind1+nb-1) = occ(ind2:ind2+nb-1)
   end do
 end do

!Change dimension in the local Dtset_cpy as well.
 dtset_cpy = Dtset%copy()
 Dtset_cpy%mpw   = kss_npw
 Dtset_cpy%mband = nbandksseff

 fform=502

 SELECT CASE (iomode)

 CASE (IO_MODE_FORTRAN)

   if (open_file(filekss, msg, newunit=kss_unt, form="unformatted") /= 0) then
     ABI_ERROR(msg)
   end if

   call my_hdr%fort_write(kss_unt, fform, ierr)
   ABI_CHECK(ierr == 0, "hdr_Fort_write returned ierr != 0")

   title='Results from ABINIT code';          write(kss_unt) title(1:80)
   title='Ab-initio plane waves calculation'; write(kss_unt) title(1:80)

   write(kss_unt) nsym2,nbandksseff,kss_npw,ishm,Psps%mpsang ! To be modified to deal with more than one projector
   write(kss_unt) (((symrel2(ii,jj,kk),ii=1,3),jj=1,3),kk=1,nsym2)
   write(kss_unt) ((tnons2(ii,kk),ii=1,3),kk=1,nsym2)
   write(kss_unt) ((gbig(ii,ig),ii=1,3),ig=1,kss_npw)
   write(kss_unt) (shlim(in),in=1,ishm)

   ! Write vkbsign for NC pseudos with Fortran IO
   ! MG FIXME: only one projector in each angular channel is treated.
   ! Moreover the allocation is done in the wrong order for dimensions...
   ! but if I change this code, compatibility with external codes is broken.
   if (Psps%usepaw==0) then
     ABI_MALLOC(vkbsign,(Psps%ntypat,Psps%mpsang))
     vkbsign(:,:)=zero
     do itypat=1,Psps%ntypat
       il0=0
       do ilmn=1,Psps%lmnmax
         il=1+Psps%indlmn(1,ilmn,itypat)
         in=Psps%indlmn(3,ilmn,itypat)
         if (il/=il0 .and. in==1) then
           il0=il
           vkbsign(itypat,il)=DSIGN(one,Psps%ekb(ilmn,itypat))
         end if
       end do
     end do
     write(kss_unt) ((vkbsign(itypat,il),il=1,Psps%mpsang),itypat=1,Psps%ntypat)
     ABI_FREE(vkbsign)
   end if

#ifdef HAVE_NETCDF
 CASE (IO_MODE_ETSF)

   ! Create file.
   NCF_CHECK(nctk_open_create(kss_unt, nctk_ncify(filekss), xmpi_comm_self))

   ! Add additional info from abinit header.
   NCF_CHECK(my_hdr%ncwrite(kss_unt, fform, nc_define=.True.))

   ! Add info on crystalline structure
   ! FIXME: Check symmorphi trick and crystal%symrel!
   ! We currently use the dataset symmetries, as defined in the Hdr structure
   ! instead of the symmetries recomputed in outkss.
   NCF_CHECK(crystal%ncwrite(kss_unt))

   ! Defined G-vectors and wavefunctions.
   call wfk_ncdef_dims_vars(kss_unt, my_hdr, fform, iskss=.True.)
   !call abi_etsf_init(Dtset_cpy, filekss, 4, .false., my_Hdr%lmn_size, Psps, Dummy_wfs)

   ! If NC pseudos, write vkbsign.
   ! Here multi-projectors are supported, array is dimensioned according to etsf-io standard.
   if (psps%usepaw == 0) then

     ! Define dims and variables needed for KB matrix elements.
     ncerr = nctk_def_dims(kss_unt, [ &
       nctkdim_t("max_number_of_angular_momenta", psps%mpsang), &
       nctkdim_t("max_number_of_projectors", psps%mproj) &
     ])
     NCF_CHECK(ncerr)

     ncerr = nctk_def_arrays(kss_unt, [ &
       nctkarr_t("kb_formfactor_sign", "int", &
&"max_number_of_projectors, max_number_of_angular_momenta, number_of_atom_species"), &
       nctkarr_t("kb_formfactors", "dp", &
&"max_number_of_coefficients, number_of_kpoints, max_number_of_projectors,&
&max_number_of_angular_momenta, number_of_atom_species"), &
       nctkarr_t("kb_formfactor_derivative", "dp", &
&"max_number_of_coefficients, number_of_kpoints, max_number_of_projectors,&
&max_number_of_angular_momenta, number_of_atom_species") &
     ])
     NCF_CHECK(ncerr)

     ABI_MALLOC(vkbsign_int, (psps%mproj, Psps%mpsang, Psps%ntypat))
     vkbsign_int=0
     do itypat=1,Psps%ntypat
       do ilmn=1,Psps%lmnmax
         il=1+Psps%indlmn(1,ilmn,itypat)
         in=Psps%indlmn(3,ilmn,itypat)
         vkbsign_int(in,il,itypat)=NINT(DSIGN(one,Psps%ekb(ilmn,itypat)))
       end do
     end do

     NCF_CHECK(nctk_set_datamode(kss_unt))

     ! Write KB sign here
     NCF_CHECK(nf90_put_var(kss_unt, nctk_idname(kss_unt, "kb_formfactor_sign"), vkbsign_int))
     ABI_FREE(vkbsign_int)
   end if

   NCF_CHECK(nctk_set_datamode(kss_unt))
#endif

 CASE DEFAULT
   ABI_ERROR(sjoin("Unsupported value for iomode:", itoa(iomode)))
 END SELECT

 call Dtset_cpy%free()
 call my_Hdr%free()

 DBG_EXIT("COLL")

end subroutine write_kss_header
!!***

!----------------------------------------------------------------------

!!****f* m_io_kss/write_vkb
!! NAME
!!  write_vkb
!!
!! FUNCTION
!!  Writes the KB form factors and derivates on file for a single k-point.
!!  Supports plain Fortran IO and netcdf with ETSF-IO format
!!
!! INPUTS
!!  kss_unt=The unit number of the file
!!  ikpt=The index of the k-point
!!  kpoint(3)=The k-point in reduced coordinates.
!!  kss_npw=Number of planewaves used for the wavefunctions in the KSS files.
!!  npw_k=Number of planewaves at this k-point in the k-centered basis set used in abinit (ecut).
!!  ecut=cutoff energy used in abinit.
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr).
!!  Psps<Pseudopotential_type>=Datatype gathering data on the Pseudopotentials.
!!  iomode=Input variables specifying the fileformat. (0-->Fortran,3-->netcdf with ETSF-IO).
!!  gbig(3,kss_npw)=Set of G-vectors used in the KSS file.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      m_io_kss
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine write_vkb(kss_unt,ikpt,kpoint,kss_npw,gbig,rprimd,Psps,iomode)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikpt,iomode,kss_npw,kss_unt
 type(Pseudopotential_type),intent(in) :: Psps
!arrays
 integer,intent(in) :: gbig(3,kss_npw)
 real(dp),intent(in) :: kpoint(3),rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: itypat,il,ig,mpsang,ntypat
!array
 real(dp),allocatable :: vkb(:,:,:),vkbd(:,:,:)
 real(dp),allocatable :: dum_vkbsign(:,:)
#ifdef HAVE_NETCDF
 integer :: ncerr,varid
 real(dp),allocatable,target :: vkb_tgt(:,:,:,:), vkbd_tgt(:,:,:,:)
#endif

! *********************************************************************

 mpsang = Psps%mpsang; ntypat = Psps%ntypat

 ABI_MALLOC(vkb ,(kss_npw,ntypat,mpsang))
 ABI_MALLOC(vkbd,(kss_npw,ntypat,mpsang))
 ABI_MALLOC(dum_vkbsign,(ntypat,mpsang))

 call kss_calc_vkb(Psps,kpoint,kss_npw,gbig,rprimd,dum_vkbsign,vkb,vkbd)
 ABI_FREE(dum_vkbsign)

 SELECT CASE (iomode)

 CASE (IO_MODE_FORTRAN)
  do itypat=1,ntypat
    do il=1,mpsang
      write(kss_unt) (vkb (ig,itypat,il),ig=1,kss_npw)
      write(kss_unt) (vkbd(ig,itypat,il),ig=1,kss_npw)
    end do
  end do

#ifdef HAVE_NETCDF
 CASE (IO_MODE_ETSF)
   ABI_MALLOC(vkb_tgt ,(kss_npw,1,mpsang,ntypat))
   ABI_MALLOC(vkbd_tgt,(kss_npw,1,mpsang,ntypat))
   do itypat=1,ntypat
     do il=1,mpsang
       do ig=1,kss_npw
         vkb_tgt (ig,1,il,itypat)=vkb (ig,itypat,il)
         vkbd_tgt(ig,1,il,itypat)=vkbd(ig,itypat,il)
       end do
     end do
   end do

   ! FIXME: Multiple projectors
   !ABI_MALLOC(vkbd, (npw, psps%lnmax, cryst%ntypat))
   !call calc_vkb(cryst,psps,kpoint,npw,npw,gvec,vkbsign,vkb,vkbd)
   !ABI_FREE(vkbsign)
   !ABI_FREE(vkb)
   !ABI_FREE(vkbd)

   ! The shape of the variable on disk is (Fortran API):
   !  (max_number_of_coefficients, number_of_kpoints, max_number_of_projectors,
   !   max_number_of_angular_momenta, number_of_atom_species)
   varid = nctk_idname(kss_unt, "kb_formfactors")
   ncerr = nf90_put_var(kss_unt, varid, vkb_tgt, start=[1,ikpt,1,1,1], count=[kss_npw,1,1,mpsang,ntypat])
   NCF_CHECK(ncerr)

   varid = nctk_idname(kss_unt, "kb_formfactor_derivative")
   ncerr = nf90_put_var(kss_unt, varid, vkbd_tgt, start=[1,ikpt,1,1,1], count=[kss_npw,1,1,mpsang,ntypat])
   NCF_CHECK(ncerr)

   ABI_FREE(vkb_tgt)
   ABI_FREE(vkbd_tgt)
#endif

 CASE DEFAULT
   ABI_ERROR(sjoin("Unsupported value for iomode:", itoa(iomode)))
 END SELECT

 ABI_FREE(vkb)
 ABI_FREE(vkbd)

end subroutine write_vkb
!!***

!----------------------------------------------------------------------

!!****f* m_io_kss/write_kss_wfgk
!! NAME
!!  write_kss_wfgk
!!
!! FUNCTION
!!  Write the Gamma-centered wavefunctions and energies on the KSS file for a single k-point.
!!  (Only the master node should call this routine).
!!
!! INPUTS
!!  kss_unt=The unit number of the file
!!  ikpt=The index of the k-point
!!  isppol=The spin index.
!!  nspinor=number of spinorial components (on current proc)
!!  kss_npw=Number of planewaves used for the wavefunctions in the KSS files.
!!  npw_k=Number of plane-waves in the k-centered basis set.
!!  nbandksseff=Number of bands to be written.
!!  natom=Number of atoms.
!!  Psps<Pseudopotential_type>=Structure gathering pseudopotential data.
!!  kpoint(3)=The k-points in reduced coordinates.
!!  ene_k(nbandksseff)=Energies at this k-point
!!  occ_k(nbandksseff)=Occupation factors at this k-point.
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr).
!!  gbig(3,kss_npw)=The set of G-vectors for the KSS wavefunctions (Gamma-centered)
!!  wfg(2,kss_npw*nspinor,nbandksseff)=The wavefunction Fourier coefficients.
!!  iomode=Input variables specifying the fileformat. (0-->Fortran,3--> netcdf with ETSF-IO format).
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      m_io_kss
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine write_kss_wfgk(kss_unt,ikpt,isppol,kpoint,nspinor,kss_npw,&
&          nbandksseff,natom,Psps,ene_k,occ_k,rprimd,gbig,wfg,Cprjnk_k,iomode)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikpt,isppol,iomode,kss_npw,nspinor,kss_unt,nbandksseff
 integer,intent(in) :: natom
 type(pseudopotential_type),intent(in) :: Psps
!arrays
 integer,intent(in) :: gbig(3,kss_npw)
 real(dp),intent(in) :: kpoint(3),rprimd(3,3)
 real(dp),intent(in) :: ene_k(nbandksseff),occ_k(nbandksseff)
 real(dp),intent(in) ::  wfg(2,kss_npw*nspinor,nbandksseff)
 type(pawcprj_type),intent(in) :: Cprjnk_k(natom,nspinor*nbandksseff*Psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: ib,ibsp,ig,ispinor,iatom,ii !,ierr
#ifdef HAVE_NETCDF
 integer :: kg_varid,cg_varid,ncerr
 character(len=nctk_slen) :: kdep
#endif

! *********************************************************************

 ! Calculate and write KB form factors and derivative at this k-point.
 if (Psps%usepaw==0) then
   call write_vkb(kss_unt,ikpt,kpoint,kss_npw,gbig,rprimd,Psps,iomode)
 end if

 ! ============================================================
 ! ==== Write wavefunctions and PAW matrix elements on disk ====
 ! ============================================================
 SELECT CASE (iomode)

 CASE (IO_MODE_FORTRAN)
   write(kss_unt) (ene_k(ib),ib=1,nbandksseff)

   ibsp=0
   do ib=1,nbandksseff
     write(kss_unt) (wfg(:,ig,ib),ig=1,kss_npw*nspinor)
     if (Psps%usepaw==1) then ! Remember that cprj are unsorted.
       do ispinor=1,nspinor
         ibsp=ibsp+1
         do iatom=1,natom
           ii=Cprjnk_k(iatom,ibsp)%nlmn
           write(kss_unt) (Cprjnk_k(iatom,ibsp)%cp(:,1:ii))
         end do
       end do
     end if
   end do

#ifdef HAVE_NETCDF
 CASE (IO_MODE_ETSF)
   if (Psps%usepaw==1) then
     ABI_WARNING("PAW output with ETSF-IO netcdf: cprj won't be written")
   end if

   ! Write G-vectors (gbig because it's not k-dependent)
   NCF_CHECK(nf90_inq_varid(kss_unt, "reduced_coordinates_of_plane_waves", kg_varid))
   NCF_CHECK(nf90_get_att(kss_unt, kg_varid, "k_dependent", kdep))
   if (kdep == "no") then
     ncerr = nf90_put_var(kss_unt, kg_varid, gbig, start=[1,1], count=[3,kss_npw])
   else
     ncerr = nf90_put_var(kss_unt, kg_varid, gbig, start=[1,1,ikpt], count=[3,kss_npw,1])
   end if
   NCF_CHECK_MSG(ncerr, "putting gibg")

   ! Write wavefunctions
   ! The coefficients_of_wavefunctions on file have shape [cplex, mpw, nspinor, mband, nkpt, nsppol]
   NCF_CHECK(nf90_inq_varid(kss_unt, "coefficients_of_wavefunctions", cg_varid))
   ncerr = nf90_put_var(kss_unt, cg_varid, wfg, start=[1,1,1,1,ikpt,isppol], &
     count=[2,kss_npw,nspinor,nbandksseff,1,1])
   NCF_CHECK_MSG(ncerr, "putting cg_k")

   ! Write eigenvalues and occupations
   NCF_CHECK(nf90_put_var(kss_unt, nctk_idname(kss_unt, "eigenvalues"), ene_k, start=[1,ikpt,isppol]))
   NCF_CHECK(nf90_put_var(kss_unt, nctk_idname(kss_unt, "occupations"), occ_k, start=[1,ikpt,isppol]))
#endif

 CASE DEFAULT
   ABI_ERROR(sjoin("Unsupported iomode:", itoa(iomode)))
 END SELECT

end subroutine write_kss_wfgk
!!***

!----------------------------------------------------------------------

!!****f* m_io_kss/k2gamma_centered
!! NAME
!!  k2gamma_centered
!!
!! FUNCTION
!!  Helper function to translate a set of wavefunctions from the k-centered G-sphere
!!  to the Gamma-centered G-sphere used for GW calculations.
!!
!! INPUTS
!!  npw_k=Number of planewaves in the k-centered basis set.
!!  kss_npw=Number of planewaves in the Gamma-centered G-sphere.
!!  nspinor=Number of spinorial component.
!!  nbandksseff=Number of bands in input-output arrays.
!!  [icg]=Shift to be used when accessing the cg array. 0 if not specified (usually k_index).
!!  [eig_vec(2,npw_k*nspinor,nbandksseff)]=wavefunctions defined on the k-centered G-sphere.
!!  [cg(2,ikg+1:ikg+npw_k*nspinor*nbandksseff)]=wavefunctions defined on the k-centered G-sphere.
!!  ngfft(18)=Info on the FFT.
!!  MPI_enreg<MPI_type>=Structure gathering info on the parallelization.
!!  istwf_k
!!  ecut
!!  gbig(3,kss_npw)
!!  kg_k(3,npw_k)
!!  gmet(3,3)
!!  kpoint(3)
!!
!! OUTPUT
!!  wfg(2,kss_npw*nspinor,nbandksseff)=Wavefunctions in the Gamma-centered representation.
!!
!! NOTES
!!  1) icg is used only if cg is present.
!!  2) cg and eig_vec are mutually exclusive. One and only one can be passed to the routine.
!!
!! PARENTS
!!      m_io_kss
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine k2gamma_centered(kpoint,npw_k,istwf_k,ecut,kg_k,kss_npw,nspinor,nbandksseff,ngfft,gmet,&
&  MPI_enreg,gbig,ug,icg,cg,eig_vec)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbandksseff,nspinor,kss_npw,npw_k,istwf_k
 integer,optional,intent(in) :: icg
 real(dp),intent(in) :: ecut
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 integer,intent(in) :: gbig(3,kss_npw)
 integer,intent(in) :: kg_k(3,npw_k)
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gmet(3,3),kpoint(3)
 real(dp),intent(out) :: ug(2,kss_npw*nspinor,nbandksseff)
 real(dp),optional,intent(in) :: eig_vec(2,npw_k*nspinor,nbandksseff)
 real(dp),optional,intent(in) :: cg(:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: tobox=1,tosph=-1
 integer :: band,ispinor,spinor_shift2,spinor_shift1,ig,my_icg,ierr
 integer :: n1,n2,n3,n4,n5,n6,ndat,full_npw_k,ii
 character(len=500) :: msg
!arrays
 integer :: identity(3,3)=RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))
 integer :: no_shift(3)=(/0,0,0/)
 integer,allocatable :: trsl(:),full_kg_k(:,:)
 real(dp),allocatable :: cfft(:,:,:,:)
 real(dp),allocatable :: full_cg(:,:),tmp_cg(:,:)

! *********************************************************************

 if (PRESENT(cg).and.PRESENT(eig_vec)) then
   ABI_ERROR("Both cg and eig_vec are present!")
 end if

! Mapping between the gamma-centered basis set and the k-centered one.
! trsl(ig)=npw_k+1 if vector ig is not inside the k-centered G-sphere.
 ABI_MALLOC(trsl,(kss_npw))

 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 n4=ngfft(4); n5=ngfft(5); n6=ngfft(6)

 if (istwf_k==1) then ! Full k-centered G-sphere.
   call table_gbig2kg(npw_k,kg_k,kss_npw,gbig,trsl,ierr)
   if (ierr/=0.and.(kss_npw>=npw_k)) then
     ABI_ERROR(' The set of G vectors is inconsistent')
   end if

 else  ! Calculate full kg with istwf_k=1 then do the mapping.
   call get_kg(kpoint,1,ecut,gmet,full_npw_k,full_kg_k)

   call table_gbig2kg(full_npw_k,full_kg_k,kss_npw,gbig,trsl,ierr)
   if (ierr/=0.and.(kss_npw>=npw_k)) then
     ABI_ERROR(' The set of G vectors is inconsistent')
   end if
 end if
 !
 ! Branching, depending on optional arguments.
 if (PRESENT(cg)) then
   my_icg=0; if (PRESENT(icg)) my_icg=icg

   SELECT CASE (istwf_k)

   CASE (1)
     do band=1,nbandksseff
       do ispinor=1,nspinor
         spinor_shift1=(ispinor-1)*kss_npw
         spinor_shift2=(ispinor-1)*npw_k
         do ig=1,kss_npw ! Retrieve the correct components
           if (trsl(ig)<=npw_k) then
             ug(:,ig+spinor_shift1,band)=cg(:,trsl(ig)+spinor_shift2+(band-1)*npw_k*nspinor+my_icg)
           else
             ug(:,ig+spinor_shift1,band)=zero
           end if
         end do
       end do
     end do

   CASE (2:9)

     ABI_CHECK(nspinor==1,"nspinor/=1!")
     !
     ! Convert input wfs from reduced to full G-sphere.
     ndat=1
     ABI_MALLOC(cfft,(2,n4,n5,n6*ndat))
     ABI_MALLOC(full_cg,(2,full_npw_k*ndat))
     ABI_MALLOC(tmp_cg,(2,npw_k*ndat))

     !write(std_out,*)"npw_k, full_kg_k",npw_k,full_npw_k

     do band=1,nbandksseff
       ii = (band-1)*npw_k
       tmp_cg = cg(:,my_icg+ii+1:my_icg+ii+npw_k)
       !write(776,*)"band= ",band,tmp_cg !cg(1:,my_icg+1+ii:my_icg+ii+npw_k)

       call sphere(tmp_cg,ndat,npw_k,cfft,n1,n2,n3,n4,n5,n6,kg_k,istwf_k,tobox,MPI_enreg%me_g0,no_shift,identity,one)

       call sphere(full_cg,ndat,full_npw_k,cfft,n1,n2,n3,n4,n5,n6,full_kg_k,1,tosph,MPI_enreg%me_g0,no_shift,identity,one)
       !write(777,*)"band= ",band,full_cg(:,:)

       do ig=1,kss_npw ! Retrieve the correct components
         if (trsl(ig)<=full_npw_k) then
           ug(:,ig,band)=full_cg(:,trsl(ig))
         else
           ug(:,ig,band)=zero
         end if
       end do
     end do !band

     ABI_FREE(cfft)
     ABI_FREE(tmp_cg)
     ABI_FREE(full_cg)

   CASE DEFAULT
     ABI_BUG("Wrong istwf_k")
   END SELECT

 else if (PRESENT(eig_vec)) then

   SELECT CASE (istwf_k)

   CASE (1)
     do band=1,nbandksseff
       do ispinor=1,nspinor
         spinor_shift1=(ispinor-1)*kss_npw
         spinor_shift2=(ispinor-1)*npw_k
         do ig=1,kss_npw ! Retrieve the correct components
           if (trsl(ig)<=npw_k) then
             ug(:,ig+spinor_shift1,band)=eig_vec(:,trsl(ig)+spinor_shift2,band)
           else
             ug(:,ig+spinor_shift1,band)=zero
           end if
         end do
       end do
     end do

   CASE DEFAULT
     write(msg,'(a,i0)')" Unsupported value for istwf_k: ",istwf_k
     ABI_ERROR(msg)
   END SELECT

 else
   ABI_ERROR("neither cg not eig_vec are in input")
 end if

 ABI_FREE(trsl)
 if (allocated(full_kg_k))  then
   ABI_FREE(full_kg_k)
 end if

end subroutine k2gamma_centered
!!***

!----------------------------------------------------------------------

!!****f* m_io_kss/make_gvec_kss
!! NAME
!! make_gvec_kss
!!
!! FUNCTION
!!   Build the list of G-vectors using the KSS convention.
!!
!! INPUTS
!!  nkpt=Number of k-points.
!!  nsym=Number of symmetries.
!!  prtvol=Verbosity option.
!!  symmorphi=
!!    0 : Old (Obsolete) implementation => Suppress inversion from symmetries list
!!    1 : Use input symrel, tnons.
!!  ecut_eff=Effective cutoff
!!  symrel(3,3,nsym)= Symmetry operation in real space.
!!  tnons(3,nsym)=Fractional translations
!!  kptns(3,nkpt)=K-points in reduced coordinates.
!!
!! OUTPUT
!!  npwkss = Input: Initial guess for the number of G-vectors required. Use 0 to have the
!!           full list of G-vectors that form a closed shell.
!!           Output: Actual number of G-vectors that form a set of closed shells
!!  gvec_kss(:,:) = Input: null pointer. Output: gvec_kss(3,npwkss), list of G-vectors (closed shells)
!!  ierr=Status error
!!
!! PARENTS
!!      m_gwls_hamiltonian,m_screening_driver,m_sigma_driver
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine make_gvec_kss(nkpt,kptns,ecut_eff,symmorphi,nsym,symrel,tnons,gprimd,prtvol,npwkss,gvec_kss,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,nsym,prtvol,symmorphi
 integer,intent(out) :: ierr
 integer,intent(inout) :: npwkss
 real(dp),intent(in) :: ecut_eff
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 integer,pointer :: gvec_kss(:,:)
 real(dp),intent(in) :: tnons(3,nsym),kptns(3,nkpt)
 real(dp),intent(in) :: gprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,ishm,maxpw,nbase
 integer :: nrst1,nrst2,nsym2,pinv
 integer,pointer :: gbig(:,:)
 character(len=500) :: msg
!arrays
 integer,pointer :: symrel2(:,:,:),shlim(:)
 real(dp),pointer :: tnons2(:,:)
! *********************************************************************

 ierr = 0
 write(msg,'(2a)')ch10,' Sorting g-vecs for an output of states on an unique "big" PW basis.'
 call wrtout(std_out,msg,'COLL')

 !ecut_eff = ecut * Dtset%dilatmx**2  ! Use ecut_eff instead of ecut_eff since otherwise
 !
 !============================================================
 !=== Prepare set containing all G-vectors sorted by stars ===
 !============================================================
 !
 !=== Analyze symmetry operations ===
 if (symmorphi==0) then  ! Old (Obsolete) implementation: Suppress inversion from symmetries list:
   nullify(symrel2,tnons2)
   call remove_inversion(nsym,symrel,tnons,nsym2,symrel2,tnons2,pinv)
   if (ANY(ABS(tnons2(:,1:nsym2))>tol8)) then
     write(msg,'(3a)')&
&     ' Non-symmorphic operations still remain in the symmetries list ',ch10,&
&     ' Program does not stop but _KSS file will not be created...'
     ABI_WARNING(msg)
     ierr=ierr+1 ; RETURN
   end if
 else if (symmorphi==1) then
!  If in the input file symmorphi==1 all the symmetry operations are retained:
!  both identity and inversion (if any) as well as non-symmorphic operations.
   nsym2=nsym ; pinv=1
   ABI_MALLOC(symrel2,(3,3,nsym))
   ABI_MALLOC(tnons2,(3,nsym))
   symrel2(:,:,:)=symrel(:,:,1:nsym)
   tnons2(:,:)   =tnons(:,1:nsym)
 else
   write(msg,'(a,i4,3a)')&
&   ' symmorphi = ',symmorphi,' while it must be 0 or 1',ch10,&
&   ' Program does not stop but KSS file will not be created...'
   ABI_WARNING(msg)
   ierr=ierr+1 ; RETURN
 end if
 !
 !===================================================================
 !==== Merge the set of k-centered G-spheres into a big set gbig ====
 !===================================================================
 !* Vectors in gbig are ordered by shells
 !
 nullify(gbig,shlim)
 call merge_and_sort_kg(nkpt,kptns,ecut_eff,nsym2,pinv,symrel2,gprimd,gbig,prtvol,shlim_p=shlim)

 nbase = SIZE(shlim)   ! Number of independent G in the big sphere.
 maxpw = shlim(nbase)  ! Total number of G"s in the big sphere.
 !
 ! * Determine optimal number of bands and G"s to be written.
 !npwkss=Dtset%npwkss
 if ((npwkss==0).or.(npwkss>=maxpw)) then
   npwkss=maxpw
   write(msg,'(5a)')&
&   ' Since the number of g''s to be written on file',ch10,&
&   ' was 0 or too large, it has been set to the max. value.,',ch10,&
&   ' computed from the union of the sets of G vectors for the different k-points.'
   call wrtout(std_out,msg,'COLL')
 end if

 ishm=0
 do ii=1,nbase
   if (shlim(ii)<=npwkss) then
     ishm=ii
   else
     EXIT
   end if
 end do
 !ishm=bisect(shlim,npwkss)

 if (shlim(ishm)/=npwkss) then
   nrst1=shlim(ishm)
   nrst2=MIN0(shlim(MIN0(ishm+1,nbase)),maxpw)
   if (IABS(npwkss-nrst2)<IABS(npwkss-nrst1)) nrst1=nrst2
   npwkss=nrst1
   if (shlim(ishm)<npwkss) ishm=ishm+1
   write(msg,'(3a)')&
&   ' The number of G''s to be written on file is not a whole number of stars ',ch10,&
&   ' the program set it to the nearest star limit.'
   call wrtout(std_out,msg,'COLL')
 end if

 write(msg,'(a,i5)')' Number of G-vectors is: ',npwkss
 call wrtout(std_out,msg,'COLL')

 ABI_MALLOC(gvec_kss,(3,npwkss))
 gvec_kss = gbig(:,1:npwkss)

 ABI_FREE(gbig)
 ABI_FREE(symrel2)
 ABI_FREE(tnons2)
 ABI_FREE(shlim)

end subroutine make_gvec_kss
!!***

!----------------------------------------------------------------------

!!****f* m_io_kss/kss_calc_vkb
!! NAME
!!  kss_calc_vkb
!!
!! FUNCTION
!!  This routine calculates the Kleynman-Bylander form factors and its derivatives
!!  needed for the evaluation of the matrix elements of the dipole operator <phi1|r|phi2>.
!!
!! INPUTS
!!  npw_k=Number of plane waves for this k-point.
!!  Psps<pseudopotential_type>=Structured datatype gathering information on the pseudopotentials.
!!  kg_k(3,npw_k)=Reduced coordinates of the G-vectors.
!!  kpoint(3)=The k-point in reduced coordinates.
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!
!! OUTPUT
!!  vkb (npw_k,Psps%ntypat,Psps%mpsang)=KB form factors.
!!  vkbd(npw_k,Psps%ntypat,Psps%mpsang)=KB form factor derivatives.
!!  vkbsign(Psps%mpsang,Psps%ntypat)   =KS dyadic sign.
!!
!! NOTES
!!  This piece of code has been extracted from outkss.F90. The implementation is consistent
!!  with the KSS file formata (Fortran version) but it presents two design flaws.
!!
!!   1) Pseudo with more that one projector per l-channel are not supported.
!!   2) Ordering of dimensions in vkb and vkbd is not optimal. We are not programming C!!!
!!
!! TODO
!!  *) Spinorial case is not implemented.
!!
!! PARENTS
!!      m_io_kss
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine kss_calc_vkb(Psps,kpoint,npw_k,kg_k,rprimd,vkbsign,vkb,vkbd)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k
 type(Pseudopotential_type),intent(in) :: Psps
!arrays
 integer,intent(in) :: kg_k(3,npw_k)
 real(dp),intent(in) :: kpoint(3),rprimd(3,3)
 real(dp),intent(out) :: vkb (npw_k,Psps%ntypat,Psps%mpsang)
 real(dp),intent(out) :: vkbd(npw_k,Psps%ntypat,Psps%mpsang)
 real(dp),intent(out) :: vkbsign(Psps%mpsang,Psps%ntypat)

!Local variables ------------------------------
!scalars
 integer :: dimffnl,ider,idir,itypat,nkpg,il0,in
 integer :: il,ilmn,ig,is
 real(dp) :: ucvol,effmass_free,ecutsm,ecut
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: ffnl(:,:,:,:),kpg_dum(:,:),modkplusg(:)
 real(dp),allocatable :: ylm(:,:),ylm_gr(:,:,:),ylm_k(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(Psps%usepaw==0,"You should not be here!")
 ABI_CHECK(Psps%useylm==0,"useylm/=0 not considered!")

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 !
 ! === Save KB dyadic sign (integer-valued) ===
 vkbsign=zero
 do itypat=1,Psps%ntypat
   il0=0
   do ilmn=1,Psps%lmnmax
     il=1+Psps%indlmn(1,ilmn,itypat)
     in=Psps%indlmn(3,ilmn,itypat)
     if (il/=il0 .and. in==1) then
       il0=il
       vkbsign(il,itypat)=DSIGN(one,Psps%ekb(ilmn,itypat))
     end if
   end do
 end do

 ! === Allocate KB form factor and derivative wrt k+G ===
 ! * Here we do not use correct ordering for dimensions

 ider=1; dimffnl=2 ! To retrieve the first derivative.
 idir=0; nkpg=0
 !
 ! Quantities used only if useylm==1
 ABI_MALLOC(ylm,(npw_k,Psps%mpsang**2*Psps%useylm))
 ABI_MALLOC(ylm_gr,(npw_k,3+6*(ider/2),Psps%mpsang**2*Psps%useylm))
 ABI_MALLOC(ylm_k,(npw_k,Psps%mpsang**2*Psps%useylm))
 ABI_MALLOC(kpg_dum,(npw_k,nkpg))

 ABI_MALLOC(ffnl,(npw_k,dimffnl,Psps%lmnmax,Psps%ntypat))

 call mkffnl(Psps%dimekb,dimffnl,Psps%ekb,ffnl,Psps%ffspl,gmet,gprimd,ider,idir,Psps%indlmn,&
   kg_k,kpg_dum,kpoint,Psps%lmnmax,Psps%lnmax,Psps%mpsang,Psps%mqgrid_ff,nkpg,npw_k,&
   Psps%ntypat,Psps%pspso,Psps%qgrid_ff,rmet,Psps%usepaw,Psps%useylm,ylm_k,ylm_gr)

 ABI_FREE(kpg_dum)
 ABI_FREE(ylm)
 ABI_FREE(ylm_gr)
 ABI_FREE(ylm_k)

 ABI_MALLOC(modkplusg,(npw_k))

 effmass_free=one; ecutsm=zero; ecut=HUGE(one)
! call mkkin(ecut,ecutsm,effmass_free,gmet,kg_k,modkplusg,kpoint,npw_k)
 call mkkin(ecut,ecutsm,effmass_free,gmet,kg_k,modkplusg,kpoint,npw_k,0,0)
 modkplusg(:)=SQRT(half/pi**2*modkplusg(:))
 modkplusg(:)=MAX(modkplusg(:),tol10)

 !do ig=1,npw_k
 ! kpg(:)= kpoint(:)+kg_k(:,ig)
 ! modkplusg(ig) = normv(kpg,gmet,"G")
 !end do

 ! Calculate matrix elements.
 vkb=zero; vkbd=zero

 do is=1,Psps%ntypat
   il0=0
   do ilmn=1,Psps%lmnmax
     il=1+Psps%indlmn(1,ilmn,is)
     in=Psps%indlmn(3,ilmn,is)
     if ((il/=il0).and.(in==1)) then
       il0=il
       if (ABS(Psps%ekb(ilmn,is))>1.0d-10) then
         if (il==1) then
           vkb (1:npw_k,is,il) = ffnl(:,1,ilmn,is)
           vkbd(1:npw_k,is,il) = ffnl(:,2,ilmn,is)*modkplusg(:)/two_pi
         else if (il==2) then
           vkb(1:npw_k,is,il)  = ffnl(:,1,ilmn,is)*modkplusg(:)
           do ig=1,npw_k
             vkbd(ig,is,il) = ((ffnl(ig,2,ilmn,is)*modkplusg(ig)*modkplusg(ig))+&
              ffnl(ig,1,ilmn,is) )/two_pi
           end do
         else if (il==3) then
           vkb (1:npw_k,is,il) =  ffnl(:,1,ilmn,is)*modkplusg(:)**2
           vkbd(1:npw_k,is,il) = (ffnl(:,2,ilmn,is)*modkplusg(:)**3+&
            2*ffnl(:,1,ilmn,is)*modkplusg(:) )/two_pi
         else if (il==4) then
           vkb (1:npw_k,is,il) =  ffnl(:,1,ilmn,is)*modkplusg(:)**3
           vkbd(1:npw_k,is,il) = (ffnl(:,2,ilmn,is)*modkplusg(:)**4+&
            3*ffnl(:,1,ilmn,is)*modkplusg(:)**2 )/two_pi
         end if
         vkb (:,is,il) = SQRT(4*pi/ucvol*(2*il-1)*ABS(Psps%ekb(ilmn,is)))*vkb (:,is,il)
         vkbd(:,is,il) = SQRT(4*pi/ucvol*(2*il-1)*ABS(Psps%ekb(ilmn,is)))*vkbd(:,is,il)
       else
         vkb (:,is,il)=zero
         vkbd(:,is,il)=zero
       end if
     end if
   end do
 end do

 ABI_FREE(ffnl)
 ABI_FREE(modkplusg)

 DBG_EXIT("COLL")

end subroutine kss_calc_vkb
!!***

!!****f* m_io_kss/outkss
!! NAME
!! outkss
!!
!! FUNCTION
!!  This routine creates an output file containing the Kohn-Sham electronic Structure
!!  for a large number of eigenstates (energies and eigen-functions).
!!  The resulting file (_KSS) is needed for a GW post-treatment.
!!
!! The routine drives the following operations:
!!  - Re-ordering G-vectors according to stars (sets of Gs related by symmetry operations).
!!    A set of g for all k-points is created.
!!  - Creating and opening the output "_KSS'" file
!!  - Printing out output file header information...
!! ... and, for each k-point:
!!    According to 'kssform', either
!!      - Re-computing <G|H|G_prim> matrix elements for all (G, G_prim).
!!        Diagonalizing H in the plane-wave basis.
!!   or - Taking eigenvalues/vectors from congugate-gradient ones.
!!  - Writing out eigenvalues and eigenvectors.
!!
!! INPUTS
!!  cg(2,mcg)=planewave coefficients of wavefunctions.
!!  usecprj=1 if cprj datastructure has been allocated (ONLY PAW)
!!  Cprj(natom,mcprj*usecprj) <type(pawcprj_type)>=
!!    projected input wave functions <Proj_i|Cnk> with all NL projectors (only for PAW)
!!    NOTE that Cprj are unsorted, see ctoprj.F90
!!  Dtfil <type(datafiles_type)>=variables related to files
!!  Dtset <type(dataset_type)>=all input variables for this dataset
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  Hdr <type(hdr_type)>=the header of wf, den and pot files
!!  kssform=govern the Kohn-Sham Structure file format
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mgfft=maximum size of 1D FFTs
!!  mkmem =number of k points treated by this node.
!!  MPI_enreg=information about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw.
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nkpt=number of k points.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nspden=number of density components
!!  nsym=number of symmetries in space group
!!  ntypat=number of types of atoms in unit cell.
!!  occ(mband*nkpt*nsppol)=occupation number for each band (usually 2) for each k.
!!  Pawtab(Psps%ntypat*Psps%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  Pawfgr<pawfgr_type>=fine grid parameters and related data
!!  prtvol=control print volume and debugging output
!!  Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  vtrial(nfft,nspden)=the trial potential
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  Output is written on file.
!!  ierr=Status error.
!!
!! NOTES
!!
!! * This routine is maintained for legacy reasons. Abinit8 is not able to read KSS files
!!   anymore hence KSS files should be used only to interface Abinit with external codes
!!   that are still using the old KSS format.
!!
!! * The routine can be time consuming (in particular when computing
!!   <G|H|G_prim> elements for all (G, G_prim)) (kssform=1).
!!   So, it is recommended to call it once per run...
!!
!! * The IO code is not parallelized and this represents a serious bottleneck when np is large.
!!
!! * when kssform==1, the routine RE-computes all Hamiltonian terms.
!!   So it is equivalent to an additional electronic SC cycle.
!!   (This has no effect is convergence was reach...
!!   If not, eigenvalues/vectors may differs from the congugaste gradient ones)
!!
!! *  The KB form factors and derivatives are not calculated correctly if there are
!!    pseudos with more than one projector in an angular momentum channel.
!!
!! * In the ETSF output format (Dtset%iomode == 3), the complete symmetry set
!!   is output. So, if reading programs need only the symmorphic symmetries, they
!!   will need to remove themselves the non-symmorphic ones.
!!
!! * There exists two file formats:
!!    kssform==1 diagonalized file _KSS in real(dp) is generated.
!!    kssform==3 same as kssform=1 but the wavefunctions are not diagonalized
!!               (they are taken from conjugate-gradient ones)
!!    Old kssform=0 and kssform=2 are obsolete and no longer available
!!
!! TESTS
!! * ETSF_IO output is tested in tests/etsf_io/t02.
!!
!! PARENTS
!!      m_outscfcv
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine outkss(crystal,Dtfil,Dtset,ecut,gmet,gprimd,Hdr,&
& kssform,mband,mcg,mcprj,mgfft,mkmem,MPI_enreg,mpsang,mpw,my_natom,natom,&
& nfft,nkpt,npwarr,nspden,nsppol,nsym,ntypat,occ,Pawtab,Pawfgr,Paw_ij,&
& prtvol,Psps,rprimd,vtrial,xred,cg,usecprj,Cprj,eigen,ierr)

 use m_linalg_interfaces

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: kssform,mband,mcg,mcprj,mgfft,mkmem,mpsang,mpw,my_natom,natom,usecprj
 integer,intent(in) :: nfft,nkpt,nsppol,nspden,nsym,ntypat,prtvol
 integer,intent(out) :: ierr
 real(dp),intent(in) :: ecut
 type(MPI_type),intent(inout) :: MPI_enreg
 type(Datafiles_type),intent(in) :: Dtfil
 type(Dataset_type),intent(in) :: Dtset
 type(Hdr_type),intent(inout) :: Hdr
 type(Pseudopotential_type),intent(in) :: Psps
 type(pawfgr_type), intent(in) :: Pawfgr
 type(crystal_t),intent(in) :: crystal
!arrays
 integer,intent(in),target :: npwarr(nkpt)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),occ(mband*nkpt*nsppol)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: vtrial(nfft,nspden)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(in) :: cg(2,mcg),eigen(mband*nkpt*nsppol)
 type(pawcprj_type),intent(in) :: Cprj(natom,mcprj*usecprj)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
 type(paw_ij_type),intent(inout),target :: Paw_ij(my_natom*Psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_rwwf=0,bufnb=20
 integer :: untkss,onband_diago
 integer :: bdtot_index,i,ib,ibp,iomode
 integer :: ibsp,ibsp1,ibsp2,ibg,ig,ii,ikpt
 integer :: master,receiver,sender,spinor_shift1,shift
 integer :: ishm,ispinor,isppol,istwf_k,my_rank,j
 integer :: k_index,maxpw,mproj,n1,n2,n2dim,n3,n4,n5,n6,nband_k
 integer :: nbandkss_k,nbandksseff,nbase,nprocs,npw_k,onpw_k,npwkss
 integer :: nrst1,nrst2,nsym2,ntemp,pinv,sizepw,spaceComm,comm_self
 integer :: pad1,pad2
 integer :: bufrt,bufsz
 real(dp) :: cinf=1.0e24,csup=zero,einf=1.0e24,esup=zero
 real(dp) :: norm,cfact,ecut_eff
 logical :: do_diago,found,ltest,lhack
 logical,parameter :: skip_test_ortho=.FALSE.
 character(len=500) :: msg
 character(len=80) :: frmt1,frmt2
 character(len=10) :: stag(2)=(/'          ','          '/)
!arrays
 integer :: nbandkssk(nkpt)
 integer,pointer :: symrel2(:,:,:)
 integer,pointer :: gbig(:,:)
 integer,pointer :: shlim(:)
 integer,allocatable :: kg_k(:,:)
 integer,allocatable :: dimlmn(:)
 integer :: nattyp_dum(0)
 real(dp) :: ovlp(2),kpoint(3),tsec(2)
 real(dp),pointer :: tnons2(:,:)
 real(dp),allocatable :: ene(:)
 real(dp),pointer :: eig_ene(:),eig_vec(:,:,:)
 real(dp),allocatable :: occ_k(:)
 real(dp),allocatable,target :: wfg(:,:,:)
 real(dp),ABI_CONTIGUOUS pointer :: ug1(:,:),ug2(:,:)
 type(pawcprj_type),allocatable :: Cprjnk_k(:,:)
 type(pawcprj_type),pointer :: Cprj_diago_k(:,:)
 type(ddiago_ctl_type) :: Diago_ctl
 type(paw_ij_type),pointer :: Paw_ij_all(:)
! *********************************************************************

 ABI_UNUSED(mkmem)

 DBG_ENTER("COLL")

 call timab(933,1,tsec) ! outkss
 call timab(934,1,tsec) ! outkss(Gsort+hd)

 spaceComm=MPI_enreg%comm_cell
 my_rank=xmpi_comm_rank(spaceComm)
 nprocs=xmpi_comm_size(spaceComm)
 master=0

 iomode = Dtset%iomode
 nullify(eig_ene)
 nullify(eig_vec)
 nullify(Cprj_diago_k)

 ! JB: Valgrind complains about non initialized value. Set to -1 so if an array
 ! should be allocated with this "unintialized value" it crashes
 onband_diago = -1

!MG: since in seq case MPI_enreg%proc_distrb is not defined
!we hack a bit the data type in order to get rid of MPI preprocessing options.
!The previous status of %proc_distrb is restored before exiting.
!Note that in case of seq run MPI_enreg%proc_distrb is nullified at the very beginning of abinit.F90
!
!FIXME this is a design flaw that should be solved: proc_distrb should always
!be allocated and filled with my_rank in case of sequential run otherwise checks like
!if (nprocs>1.and.MPI_enreg%proc_distrb(ii)==me) leads to SIGFAULT under gfortran.
!as the second array is not allocated.
 lhack=.FALSE.
 if (nprocs==1) then
   ltest=allocated(MPI_enreg%proc_distrb)
   if (.not.ltest) then
     ABI_MALLOC(MPI_enreg%proc_distrb,(nkpt,mband,nsppol))
     MPI_enreg%proc_distrb=my_rank
     lhack=.TRUE.
   end if
   ltest=ALL(MPI_enreg%proc_distrb==my_rank)
   ABI_CHECK(ltest,'wrong values in %proc_distrb')
 end if
!
!============================
!==== Perform some tests ====
!============================
 ierr=0

 if (iomode==IO_MODE_ETSF) then
   write(msg,'(3a)')&
&   'when iomode==3 in outkss, support for netcdf ',ch10,&
&   'must be compiled. Use --enable-netcdf when configuring '
#ifndef HAVE_NETCDF
   ABI_WARNING(msg)
   ierr = ierr + 1
#endif
 end if

 if (kssform==3) then
   write(msg,'(a,70("="),4a)')ch10,ch10,&
&   ' Calculating and writing out Kohn-Sham electronic Structure file',ch10, &
&   ' Using conjugate gradient wavefunctions and energies (kssform=3)'
 else if (kssform==1) then
   write(msg,'(a,70("="),4a,i1,a)') ch10,ch10, &
&   ' Calculating and writing out Kohn-Sham electronic Structure file',ch10, &
&   ' Using diagonalized wavefunctions and energies (kssform=',kssform,')'
 else
   write(msg,'(a,i0,2a)')&
&   " Unsupported value for kssform: ",kssform,ch10,&
&   "  Program does not stop but _KSS file will not be created..."
   ierr=ierr+1
 end if
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
!
!* Check whether nband is constant in metals
 if ( (Dtset%occopt>=2.and.Dtset%occopt<=8) .and. (ANY(Dtset%nband(1:nkpt*nsppol)/=Dtset%nband(1))) ) then
   write(msg,'(3a,i4,a,i3,a,i4,3a)')&
&   ' The number of bands must be the same for all k-points ',ch10,&
&   ' but nband(1)=',Dtset%nband(1),' is different of nband(',&
&   ikpt+(isppol-1)*nkpt,')=',Dtset%nband(ikpt+(isppol-1)*nkpt),'.',ch10,&
&   '  Program does not stop but _KSS file will not be created...'
   ABI_WARNING(msg)
   ierr=ierr+1
 end if
!* istwfk must be 1 for each k-point
 if (ANY(Dtset%istwfk(1:nkpt)/=1).and.kssform/=3) then
   write(msg,'(7a)')&
&   ' istwfk/=1 not allowed when kssform/=3 :',ch10,&
&   ' States output not programmed for time-reversal symmetry.',ch10,&
&   ' Action : change istwfk in input file (put it to 1 for all kpt).',ch10,&
&   ' Program does not stop but _KSS file will not be created...'
   ABI_WARNING(msg)
   ierr=ierr+1
 end if
!* Check spin-orbit
 if (Psps%mpssoang/=mpsang) then
   write(msg,'(3a)')&
&   ' Variable mpspso should be 1 !',ch10,&
&   ' Program does not stop but _KSS file will not be created...'
   ABI_WARNING(msg)
   ierr=ierr+1
 end if
!* Check mproj
 mproj=MAXVAL(Psps%indlmn(3,:,:))
 if (mproj>1.and.Psps%usepaw==0) then ! TODO One has to derive the expression for [Vnl,r], in particular HGH and GTH psps
   write(msg,'(8a)')ch10,&
&   ' outkss : COMMENT - ',ch10,&
&   ' At least one NC pseudopotential has more that one projector per angular channel',ch10,&
&   ' Note that inclvkb==0 should be used in screening, since the evaluation of the commutator',ch10,&
&   ' for this particular case is not implemented yet'
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if
!* Check max angular momentum
 if (MAXVAL(Psps%indlmn(1,:,:))+1 >= 5) then
   write(msg,'(3a)')&
&   ' Pseudopotentials with f-projectors not implemented',ch10,&
&   ' Program does not stop but _KSS file will not be created...'
   ABI_WARNING(msg)
   ierr=ierr+1
 end if
!* Check useylm
 if (Psps%useylm/=0.and.Psps%usepaw==0) then
   write(msg,'(3a)')&
&   ' The present version of outkss does not work with useylm/=0 !',ch10,&
&   ' Program does not stop but _KSS file will not be created...'
   ABI_WARNING(msg)
   ierr=ierr+1
 end if
!* Check PAW and kssform value
 if (Psps%usepaw/=0) then
   if (nprocs>1.and.kssform==1) then
     write(msg,'(3a)')&
&     ' Parallel PAW with kssform=1, not yet allowed',ch10,&
&     ' Program does not stop but _KSS file will not be created...'
     ABI_WARNING(msg)
     ierr=ierr+1
   end if
   if (kssform==3.and.usecprj/=1) then
     write(msg,'(3a)')&
&     ' If PAW and kssform=3, usecprj must be 1',ch10,&
&     ' Program does not stop but _KSS file will not be created...'
     ABI_WARNING(msg)
     ierr=ierr+1
   end if
 end if
!* Check parallelization
 if (MPI_enreg%paralbd/=0) then
   write(msg,'(3a)')&
&   ' outkss cannot be used with parallelization on bands (paralbd/=0) !',ch10,&
&   ' Program does not stop but _KSS file will not be created...'
   ABI_WARNING(msg)
   ierr=ierr+1
 end if
 if (MPI_enreg%paral_spinor/=0) then
   write(msg,'(3a)')&
&   ' outkss cannot be used yet with parallelization on nspinors !',ch10,&
&   ' Program does not stop but _KSS file will not be created...'
   ABI_WARNING(msg)
   ierr=ierr+1

 endif
 if (ierr/=0) then
   write(msg,'(3a)')&
&   ' outkss: Not allowed options found !',ch10,&
&   ' Program does not stop but _KSS file will not be created...'
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(a)')' outkss: see the log file for more information.'
   call wrtout(ab_out,msg,'COLL')
   RETURN ! Houston we have a problem!
 end if
!
!Estimate required memory in case of diagonalization.
!TODO to be modified to take into account the case nsppol=2
 if (kssform/=3) then
   call memkss(mband,mgfft,mproj,Psps%mpssoang,mpw,natom,Dtset%ngfft,nkpt,dtset%nspinor,nsym,ntypat)
 end if
!
!=== Initialize some variables ===
 if (nsppol==2) stag(:)=(/'SPIN UP:  ','SPIN DOWN:'/)
 n1=Dtset%ngfft(1); n2=Dtset%ngfft(2); n3=Dtset%ngfft(3)
 n4=Dtset%ngfft(4); n5=Dtset%ngfft(5); n6=Dtset%ngfft(6)
 ecut_eff = ecut * Dtset%dilatmx**2  ! Use ecut_eff instead of ecut_eff since otherwise
!one cannot restart from a previous density file
 sizepw=2*mpw ; do_diago=(kssform/=3)
 ABI_MALLOC(dimlmn,(natom*Psps%usepaw))
 if (Psps%usepaw==1) then
   call pawcprj_getdim(dimlmn,natom,nattyp_dum,ntypat,Dtset%typat,pawtab,'R')
 end if
!
!============================================================
!=== Prepare set containing all G-vectors sorted by stars ===
!============================================================
 write(msg,'(2a)')ch10,' Sorting g-vecs for an output of states on an unique "big" PW basis.'
 call wrtout(std_out,msg,'COLL')
!
!=== Analyze symmetry operations ===
 if (Dtset%symmorphi==0) then  ! Old (Obsolete) implementation: Suppress inversion from symmetries list:
   nullify(symrel2,tnons2)
   call remove_inversion(nsym,Dtset%symrel,Dtset%tnons,nsym2,symrel2,tnons2,pinv)
   if (ANY(ABS(tnons2(:,1:nsym2))>tol8)) then
     write(msg,'(3a)')&
&     ' Non-symmorphic operations still remain in the symmetries list ',ch10,&
&     ' Program does not stop but _KSS file will not be created...'
     ABI_WARNING(msg)
     ierr=ierr+1 ; RETURN
   end if
 else if (Dtset%symmorphi==1) then
!  If in the input file symmorphi==1 all the symmetry operations are retained:
!  both identity and inversion (if any) as well as non-symmorphic operations.
   nsym2=nsym ; pinv=1
   ABI_MALLOC(symrel2,(3,3,nsym))
   ABI_MALLOC(tnons2,(3,nsym))
   symrel2(:,:,:)=Dtset%symrel(:,:,1:nsym)
   tnons2(:,:)   =Dtset%tnons(:,1:nsym)
 else
   write(msg,'(a,i4,3a)')&
&   ' symmorphi = ',Dtset%symmorphi,' while it must be 0 or 1',ch10,&
&   ' Program does not stop but KSS file will not be created...'
   ABI_WARNING(msg)
   ierr=ierr+1 ; RETURN
 end if
!
!===================================================================
!==== Merge the set of k-centered G-spheres into a big set gbig ====
!===================================================================
!* Vectors in gbig are ordered by shells
!
 nullify(gbig,shlim)
 call merge_and_sort_kg(nkpt,Dtset%kptns,ecut_eff,nsym2,pinv,symrel2,gprimd,gbig,prtvol,shlim_p=shlim)

 nbase = SIZE(shlim)   ! Number of independent G in the big sphere.
 maxpw = shlim(nbase)  ! Total number of G"s in the big sphere.
!
!* Determine optimal number of bands and G"s to be written.
 npwkss=Dtset%npwkss
 if ((npwkss==0).or.(npwkss>=maxpw)) then
   npwkss=maxpw
   write(msg,'(5a)')&
&   ' Since the number of g''s to be written on file',ch10,&
&   ' was 0 or too large, it has been set to the max. value.,',ch10,&
&   ' computed from the union of the sets of G vectors for the different k-points.'
   call wrtout(std_out,msg,'COLL')
 end if

 ishm=0
 do ii=1,nbase
   if (shlim(ii)<=npwkss) then
     ishm=ii
   else
     EXIT
   end if
 end do

 if (shlim(ishm)/=npwkss) then
   nrst1=shlim(ishm)
   nrst2=MIN0(shlim(MIN0(ishm+1,nbase)),maxpw)
   if (IABS(npwkss-nrst2)<IABS(npwkss-nrst1)) nrst1=nrst2
   npwkss=nrst1
   if (shlim(ishm)<npwkss) ishm=ishm+1
   write(msg,'(3a)')&
&   ' The number of G''s to be written on file is not a whole number of stars ',ch10,&
&   ' the program set it to the nearest star limit.'
   call wrtout(std_out,msg,'COLL')
 end if

 write(msg,'(a,i5)')' Number of g-vectors written on file is: ',npwkss
 call wrtout(std_out,msg,'COLL')
!
!=== Check on the number of stored bands ===
 if (do_diago) then

   if (Dtset%nbandkss==-1.or.Dtset%nbandkss>=maxpw) then
     nbandkssk(1:nkpt)=npwarr(1:nkpt)
     write(msg,'(6a)')ch10,&
&     ' Since the number of bands to be computed was (-1) or',ch10,&
&     ' too large, it has been set to the max. value. allowed for each k,',ch10,&
&     ' thus, the minimum of the number of plane waves for each k point.'
     call wrtout(std_out,msg,'COLL')
   else
     nbandkssk(1:nkpt)=Dtset%nbandkss
     found=.FALSE.
     do ikpt=1,nkpt
       if (Dtset%nbandkss>npwarr(ikpt)) then
         nbandkssk(ikpt)=npwarr(ikpt)
         found=.TRUE.
       end if
     end do
     if (found) then
       write(msg,'(7a)')&
&       ' The value choosen for the number of bands in file',ch10,&
&       ' (nbandkss) was greater than at least one number of plane waves ',ch10,&
&       ' for a given k-point (npw_k).',ch10,' It has been modified consequently.'
       ABI_WARNING(msg)
     end if
   end if
   found=.FALSE.
   do ikpt=1,nkpt
     if (nbandkssk(ikpt)>npwkss) then
       nbandkssk(ikpt)=npwkss
       found=.TRUE.
     end if
   end do
   if (found) then
     write(msg,'(5a)')&
&     ' The number of bands to be computed (for one k) was',ch10,&
&     ' greater than the number of g-vectors to be written.',ch10,&
&     ' It has been modified consequently.'
     ABI_WARNING(msg)
   end if
   nbandksseff=MINVAL(nbandkssk)

 else ! .not. do_diago
   do ikpt=1,nkpt
     do isppol=1,nsppol
       nbandkssk(ikpt)=Dtset%nband(ikpt+(isppol-1)*nkpt)
     end do
   end do
   nbandksseff=MINVAL(nbandkssk)
   if (Dtset%nbandkss>0 .and. Dtset%nbandkss<nbandksseff) then
     write(msg,'(a,i5,a,i5,2a)')&
&     ' Number of bands calculated=',nbandksseff,', greater than nbandkss=',Dtset%nbandkss,ch10,&
&     ' will write nbandkss bands on the KSS file'
     ABI_COMMENT(msg)
     nbandksseff=Dtset%nbandkss
   end if
 end if

 write(msg,'(a,i5)')' Number of bands written on file is: ',nbandksseff
 call wrtout(std_out,msg,'COLL')

 found= ANY(nbandkssk(1:nkpt)<npwarr(1:nkpt))

 if (do_diago) then
   if (found) then
     write(msg,'(6a)')ch10,&
&     ' Since the number of bands to be computed',ch10,&
&     ' is less than the number of G-vectors found,',ch10,&
&     ' the program will perform partial diagonalizations.'
   else
     write(msg,'(6a)')ch10,&
&     ' Since the number of bands to be computed',ch10,&
&     ' is equal to the nb of G-vectors found for each k-pt,',ch10,&
&     ' the program will perform complete diagonalizations.'
   end if
   call wrtout(std_out,msg,'COLL')
 end if
!
!==========================================================================
!=== Open KSS file for output, write header with dimensions and kb sign ===
!==========================================================================
!
!* Output required disk space.
 call dsksta(ishm,Psps%usepaw,nbandksseff,mpsang,natom,ntypat,npwkss,nkpt,dtset%nspinor,nsppol,nsym2,dimlmn)

 if (my_rank==master) then
   call write_kss_header(dtfil%fnameabo_kss,npwkss,ishm,nbandksseff,mband,nsym2,symrel2,tnons2,occ,gbig,shlim,&
&   crystal,Dtset,Hdr,Psps,iomode,untkss)
 end if

 ABI_FREE(shlim)

 if (     do_diago) msg = ' Diagonalized eigenvalues'
 if (.not.do_diago) msg = ' Conjugate gradient eigenvalues'
 call wrtout(ab_out,msg,'COLL')

 if (Dtset%enunit==1) then
   msg='   k    eigenvalues [eV]'
 else
   msg='   k    eigenvalues [Hartree]'
 end if
 call wrtout(ab_out,msg,'COLL')
!
!=== In case of PAW distributed atomic sites, need to retrieve the full paw_ij%dij ===
 if (do_diago.and.Psps%usepaw==1.and.MPI_enreg%nproc_atom>1) then
   ABI_MALLOC(Paw_ij_all,(dtset%natom))
   call paw_ij_gather(Paw_ij,Paw_ij_all,-1,MPI_enreg%comm_atom)
 else
   paw_ij_all => paw_ij
 end if


 call timab(934,2,tsec) ! outkss(Gsort+hd)
!

 k_index=0; bdtot_index=0; ibg=0

 do isppol=1,nsppol ! Loop over spins
!
   do ikpt=1,nkpt ! Loop over k-points.
     call timab(935,1,tsec) ! outkss(k-Loop)

     nband_k   =Dtset%nband(ikpt+(isppol-1)*nkpt)
     npw_k     =npwarr(ikpt)
     istwf_k   =Dtset%istwfk(ikpt)
     kpoint    =Dtset%kptns(:,ikpt)
     nbandkss_k=nbandkssk(ikpt)

     ! Get G-vectors, for this k-point.
     call get_kg(kpoint,istwf_k,ecut_eff,gmet,onpw_k,kg_k)
     ABI_CHECK(onpw_k==npw_k,"Mismatch in npw_k")
!
!    ============================================
!    ==== Parallelism over k-points and spin ====
!    ============================================
     if (MPI_enreg%proc_distrb(ikpt,1,isppol)==my_rank) then

       write(msg,'(2a,i3,3x,a)')ch10,' k-point ',ikpt,stag(isppol)
       call wrtout(std_out,msg,'PERS')

       if (do_diago) then
         ! Direct diagonalization of the KS Hamiltonian.
         ABI_SFREE_PTR(eig_ene)
         ABI_SFREE_PTR(eig_vec)
         comm_self = xmpi_comm_self

         call timab(936,1,tsec)

         call init_ddiago_ctl(Diago_ctl,"Vectors",isppol,dtset%nspinor,ecut_eff,Dtset%kptns(:,ikpt),Dtset%nloalg,gmet,&
           nband_k=nbandkssk(ikpt),effmass_free=Dtset%effmass_free,istwf_k=Dtset%istwfk(ikpt),prtvol=Dtset%prtvol)

         call ksdiago(Diago_ctl,nbandkssk(ikpt),Dtset%nfft,mgfft,Dtset%ngfft,natom,&
           Dtset%typat,nfft,dtset%nspinor,nspden,nsppol,Pawtab,Pawfgr,Paw_ij_all,&
           Psps,rprimd,vtrial,xred,onband_diago,eig_ene,eig_vec,Cprj_diago_k,comm_self,ierr)

         call timab(936,2,tsec)
       end if

     end if ! END of kpt+spin parallelism.
!
!    ===========================================================
!    ==== Transfer data between master and the working proc ====
!    ===========================================================
     call timab(937,1,tsec) !outkss(MPI_exch)
     ABI_MALLOC(Cprjnk_k,(0,0))
     if (nprocs==1) then

       if (Psps%usepaw==1) then ! Copy projectors for this k-point
         n2dim=min(nbandksseff*dtset%nspinor,onband_diago)
         if (kssform==3) n2dim=nband_k*dtset%nspinor
         ABI_FREE(Cprjnk_k)
         ABI_MALLOC(Cprjnk_k,(natom,n2dim))
         call pawcprj_alloc(Cprjnk_k,0,dimlmn)
         if (kssform==3) then
           call pawcprj_copy(Cprj(:,ibg+1:ibg+dtset%nspinor*nband_k),Cprjnk_k)
         else
           !ABI_WARNING("Here I have to use onband_diago") !FIXME
           call pawcprj_copy(Cprj_diago_k(:,1:n2dim),Cprjnk_k)
         end if
       end if

     else
       !parallel case

       receiver=master; sender=MPI_enreg%proc_distrb(ikpt,1,isppol)

       bufsz=nbandksseff/bufnb; bufrt=nbandksseff-bufnb*bufsz

       if (my_rank==receiver.or.my_rank==sender) then

         if (do_diago.and.(my_rank==receiver.and.my_rank/=sender)) then ! Alloc arrays if not done yet.
           ABI_MALLOC(eig_ene,(npw_k*dtset%nspinor))
           ABI_MALLOC(eig_vec,(2,npw_k*dtset%nspinor,nbandkssk(ikpt)))
         end if

         if (.not.do_diago) then

           ABI_MALLOC(eig_vec,(2,npw_k*dtset%nspinor,nbandkssk(ikpt)))

           if (my_rank==sender) then
             do ib=1,nbandksseff
               shift = k_index + (ib-1)*npw_k*dtset%nspinor
               do ig=1,npw_k*dtset%nspinor
                 eig_vec(:,ig,ib)=cg(:,ig+shift)
               end do
             end do
           end if
!
!          In case of PAW and kssform==3, retrieve matrix elements of the PAW projectors for this k-point
           if (Psps%usepaw==1) then
             n2dim=min(nbandksseff*dtset%nspinor,onband_diago)
             if (kssform==3) n2dim=nband_k*dtset%nspinor
             ABI_FREE(Cprjnk_k)
             ABI_MALLOC(Cprjnk_k,(natom,n2dim))
             call pawcprj_alloc(Cprjnk_k,0,dimlmn)
             if (my_rank==sender) then
               if (kssform==3) then
                 call pawcprj_copy(Cprj(:,ibg+1:ibg+dtset%nspinor*nband_k),Cprjnk_k)
               else
                !ABI_WARNING("Here I have to use onband_diago") !FIXME
                 call pawcprj_copy(Cprj_diago_k(:,1:n2dim),Cprjnk_k)
               end if
             end if
             if (sender/=receiver) then
               call pawcprj_mpi_exch(natom,n2dim,dimlmn,0,Cprjnk_k,Cprjnk_k,sender,receiver,spaceComm,ierr)
             end if
           end if ! usepaw

         else ! do_diago
           call xmpi_exch(eig_ene,nbandksseff,sender,eig_ene,receiver,spaceComm,ierr)
         end if

!        Exchange eigenvectors.
         if (bufsz>0) then
           do i=0,bufnb-1
             call xmpi_exch(eig_vec(:,:,i*bufsz+1:(i+1)*bufsz),2*npw_k*dtset%nspinor*bufsz,&
&             sender,eig_vec(:,:,i*bufsz+1:(i+1)*bufsz),receiver,spaceComm,ierr)
           end do
         end if
         if (bufrt>0) then
           call xmpi_exch(eig_vec(:,:,bufnb*bufsz+1:bufnb*bufsz+bufrt),2*npw_k*dtset%nspinor*bufrt,&
&           sender,eig_vec(:,:,bufnb*bufsz+1:bufnb*bufsz+bufrt),receiver,spaceComm,ierr)
         end if

       end if
     end if !nprocs > 1
     call timab(937,2,tsec) !outkss(MPI_exch)

     call timab(938,1,tsec) !outkss(write)

     if (my_rank==master) then ! Prepare data for writing on disk.
       ABI_MALLOC(ene,(nbandksseff))
       ABI_MALLOC(wfg,(2,npwkss*dtset%nspinor,nbandksseff))
       ene=zero; wfg=zero

       if (.not.do_diago) then
         ene(1:nbandksseff)=eigen(1+bdtot_index:nbandksseff+bdtot_index)

         if (nprocs>1) then
           call k2gamma_centered(kpoint,npw_k,istwf_k,ecut_eff,kg_k,npwkss,dtset%nspinor,nbandksseff,Dtset%ngfft,gmet,&
&           MPI_enreg,gbig,wfg,eig_vec=eig_vec)
         else
           call k2gamma_centered(kpoint,npw_k,istwf_k,ecut_eff,kg_k,npwkss,dtset%nspinor,nbandksseff,Dtset%ngfft,gmet,&
&           MPI_enreg,gbig,wfg,icg=k_index,cg=cg)
         end if

       else ! Direct diagonalization.
         ene(1:nbandksseff)=eig_ene(1:nbandksseff)

!        FIXME: receiver does not know Diago_ctl%npw_k
         call k2gamma_centered(kpoint,npw_k,istwf_k,ecut_eff,kg_k,npwkss,dtset%nspinor,nbandksseff,Dtset%ngfft,gmet,&
&         MPI_enreg,gbig,wfg,eig_vec=eig_vec)

!        * Check diagonalized eigenvalues with respect to conjugate gradient ones
         ntemp=MIN(nbandksseff,nband_k)
         if (ANY(ABS(ene(1:ntemp)-eigen(1+bdtot_index:ntemp+bdtot_index))>tol3)) then
           write(msg,'(3a)')&
&           ' The diagonalized eigenvalues differ by more than 10^-3 Hartree',ch10,&
&           ' with respect to the conjugated gradient values.'
           ABI_WARNING(msg)
         end if
       end if
!
!      * Write out energies
       if (Dtset%enunit==1) then
         cfact=Ha_eV ; frmt1='(i4,4x,9(1x,f7.2))' ; frmt2='(8x,9(1x,f7.2))'
         write(msg,'(a,i3,3x,a)')' Eigenvalues in eV for ikpt= ',ikpt,stag(isppol)
       else
         cfact=one   ; frmt1='(i4,4x,9(1x,f7.4))' ; frmt2='(8x,9(1x,f7.4))'
         write(msg,'(a,i3,3x,a)')' Eigenvalues in Hartree for ikpt= ',ikpt,stag(isppol)
       end if
       call wrtout(std_out,msg,'COLL')

       write(msg,frmt1)ikpt,(ene(ib)*cfact,ib=1,MIN(9,nbandksseff))
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out,msg,'COLL')

       if (nbandksseff>9) then
         do j=10,nbandksseff,9
           write(msg,frmt2) (ene(ib)*cfact,ib=j,MIN(j+8,nbandksseff))
           call wrtout(std_out,msg,'COLL')
           call wrtout(ab_out,msg,'COLL')
         end do
       end if

       if (skip_test_ortho) then ! Set this if to FALSE to skip test below
         einf=one; esup=one; cinf=zero; csup=zero
       else
         !
         ! Test on the normalization of wavefunctions.
         ibsp=0
         do ib=1,nbandksseff
           norm=zero
           do ispinor=1,dtset%nspinor
             ibsp=ibsp+1
             spinor_shift1=(ispinor-1)*npwkss
             ug1 => wfg(:,1+spinor_shift1:npwkss+spinor_shift1,ib)

             !ovlp(1) =ddot(npwkss,ug1(1,:),1,ug1(1,:),1) + ddot(npwkss,ug1(2,:),1,ug1(2,:),1)
             ovlp(1) = cg_dznrm2(npwkss,ug1)
             ovlp(1) = ovlp(1)**2
             ovlp(2) = zero
             if (Psps%usepaw==1) ovlp = ovlp &
&               + paw_overlap(Cprjnk_k(:,ibsp:ibsp),Cprjnk_k(:,ibsp:ibsp),Dtset%typat,Pawtab,&
&                             spinor_comm=MPI_enreg%comm_spinor)
             norm = norm + DABS(ovlp(1))
           end do
           if (norm<einf) einf=norm
           if (norm>esup) esup=norm
         end do
!
!        Test on the orthogonalization of wavefunctions.
         do ib=1,nbandksseff
           pad1=(ib-1)*dtset%nspinor
           do ibp=ib+1,nbandksseff
             pad2=(ibp-1)*dtset%nspinor
             ovlp(:)=zero
             do ispinor=1,dtset%nspinor
               ibsp1=pad1+ispinor
               ibsp2=pad2+ispinor
               spinor_shift1=(ispinor-1)*npwkss
               ug1 => wfg(:,1+spinor_shift1:npwkss+spinor_shift1,ib )
               ug2 => wfg(:,1+spinor_shift1:npwkss+spinor_shift1,ibp)

               !ovlp(1)=ddot(npwkss,ug1(1,:),1,ug2(1,:),1) + ddot(npwkss,ug1(2,:),1,ug2(2,:),1)
               !ovlp(2)=ddot(npwkss,ug1(1,:),1,ug2(2,:),1) - ddot(npwkss,ug1(2,:),1,ug2(1,:),1)
               ovlp = cg_zdotc(npwkss,ug1,ug2)

               if (Psps%usepaw==1) ovlp= ovlp &
&                 + paw_overlap(Cprjnk_k(:,ibsp1:ibsp1),Cprjnk_k(:,ibsp2:ibsp2),Dtset%typat,Pawtab,&
&                               spinor_comm=MPI_enreg%comm_spinor)
             end do
             norm = DSQRT(ovlp(1)**2+ovlp(2)**2)
             if (norm<cinf) cinf=norm
             if (norm>csup) csup=norm
           end do
         end do
       end if

       write(msg,'(a,i3,3x,a)')' Writing out eigenvalues/vectors for ikpt=',ikpt,stag(isppol)
       call wrtout(std_out,msg,'COLL')
!
!      * Write occupation numbers on std_out.
       ABI_MALLOC(occ_k,(MAX(nband_k,nbandksseff)))
       occ_k(1:nband_k)=occ(1+bdtot_index:nband_k+bdtot_index)
       if (nband_k < nbandksseff) occ_k(nband_k+1:nbandksseff)=zero

       write(msg,'(a,i3,3x,a)')' Occupation numbers for ikpt=',ikpt,stag(isppol)
       call wrtout(std_out,msg,'COLL')
       write(msg,'(i4,4x,9(1x,f7.4))')ikpt,(occ_k(ib),ib=1,MIN(9,nbandksseff))
       call wrtout(std_out,msg,'COLL')
       if (nbandksseff>9) then
         do j=10,nbandksseff,9
           write(msg,'(8x,9(1x,f7.4))') (occ_k(ib),ib=j,MIN(j+8,nbandksseff))
           call wrtout(std_out,msg,'COLL')
         end do
       end if
!
!      =================================================================
!      ==== Write wavefunctions, KB and PAW matrix elements on disk ====
!      =================================================================
       call write_kss_wfgk(untkss,ikpt,isppol,kpoint,dtset%nspinor,npwkss,&
             nbandksseff,natom,Psps,ene,occ_k,rprimd,gbig,wfg,Cprjnk_k,iomode)

       ABI_FREE(occ_k)
       ABI_FREE(ene)
       ABI_FREE(wfg)

     end if ! my_rank==master
     call timab(938,2,tsec) !outkss(write)

     if (my_rank==master.or.my_rank==MPI_enreg%proc_distrb(ikpt,1,isppol)) then
       ABI_SFREE_PTR(eig_ene)
       ABI_SFREE_PTR(eig_vec)
       if (Psps%usepaw==1) call pawcprj_free(Cprjnk_k)
     end if
     ABI_FREE(Cprjnk_k)
     ABI_SFREE(kg_k)

!    if (MPI_enreg%paral_compil_kpt==1) then !cannot be used in seq run!
     if (.not.(proc_distrb_cycle(MPI_enreg%proc_distrb,ikpt,1,nband_k,isppol,my_rank))) then
       k_index=k_index+npw_k*nband_k*dtset%nspinor
       ibg=ibg+dtset%nspinor*nband_k
     end if
     bdtot_index=bdtot_index+nband_k

     call xmpi_barrier(spaceComm) ! FIXME this barrier is detrimental in the case of direct diago!

     call timab(935,2,tsec) !outkss(k-loop)
   end do ! ! End loop over k-points.
 end do ! spin

 write(msg,'(3a,f9.6,2a,f9.6,4a,f9.6,2a,f9.6,a)')&
& ' Test on the normalization of the wavefunctions',ch10,&
& '  min sum_G |a(n,k,G)| = ',einf,ch10,&
& '  max sum_G |a(n,k,G)| = ',esup,ch10,&
& ' Test on the orthogonalization of the wavefunctions',ch10,&
& '  min sum_G a(n,k,G)a(n'',k,G) = ',cinf,ch10,&
& '  max sum_G a(n,k,G)a(n'',k,G) = ',csup,ch10
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 ABI_FREE(gbig)
 ABI_FREE(symrel2)
 ABI_FREE(tnons2)
 if (Psps%usepaw==1)  then
   ABI_FREE(dimlmn)
   if (do_diago.and.MPI_enreg%nproc_atom>1) then
     ABI_FREE(Paw_ij_all)
   end if
 end if
!
!* Close file
 if (my_rank==master) then
   if (iomode==IO_MODE_FORTRAN) close(unit=untkss)
#if defined HAVE_NETCDF
   if (iomode==IO_MODE_ETSF) then
     NCF_CHECK(nf90_close(untkss))
   end if
#endif
 end if

 if (associated(Cprj_diago_k)) then
   call pawcprj_free(Cprj_diago_k)
   ABI_FREE(Cprj_diago_k)
 end if

 if (lhack)  then
   ABI_FREE(MPI_enreg%proc_distrb)
 end if

 call wrtout(std_out, "outkss done", "COLL")
 call xmpi_barrier(spaceComm)

 DBG_EXIT("COLL")
 call timab(933,2,tsec) ! outkss

contains
!!***

!!****f* ABINIT/memkss
!! NAME
!! memkss
!!
!! FUNCTION
!! This routine evaluates the additional amount of memory required
!! by routine 'outkss'.
!!
!! INPUTS
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mproj=maximum dimension for number of projection operators for each
!!   angular momentum for nonlocal pseudopotential
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw.
!!  natom=number of atoms in cell.
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkpt=number of k points.
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsym=number of symmetries in space group
!!  ntypat=number of types of atoms in unit cell.
!!
!! NOTES
!! This routine is not available for paw calculations
!!
!! PARENTS
!!      m_io_kss
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine memkss(mband,mgfft,mproj,mpsang,mpw,natom,ngfft,nkpt,nspinor,nsym,ntypat)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mgfft,mproj,mpsang,mpw,natom,nkpt,nspinor
 integer,intent(in) :: nsym,ntypat
!arrays
 integer,intent(in) :: ngfft(18)

!Local variables-------------------------------
!scalars
 integer(i8b) :: isize,memsize
 character(len=500) :: msg

! *********************************************************************
!
 isize=580+fnlen+4*(81+nkpt+9*nsym)+8*15    !non allocatable var.
 if(xmpi_paral==1)then
   isize=isize+4*4                           !kpt_distrb
 end if
 memsize=isize
 isize=isize+4*nkpt+12*mpw+20*nkpt*mpw      !nbasek,gbasek,cnormk,gcurr
 memsize=max(memsize,isize)
 if(xmpi_paral==1)then
   isize=isize+12*mpw*nkpt                   !ibuf1,ibuf2,rbuf1
   memsize=max(memsize,isize)
   isize=isize-12*mpw*nkpt                   !ibuf1,ibuf2,rbuf1
 end if
 isize=isize+40*mpw                         !gbase,cnorm
 memsize=max(memsize,isize)
 isize=isize-4*nkpt-20*mpw*nkpt             !nbasek,gbasek,cnormk
 isize=isize+4*mpw                          !insort
 memsize=max(memsize,isize)
 isize=isize-16*mpw                         !cnorm
 isize=isize+28*mpw+24*nsym                 !gbig,nshell,gshell
 memsize=max(memsize,isize)
 isize=isize+4*mpw                          !shlim
 memsize=max(memsize,isize)
 isize=isize-44*mpw-24*nsym                 !gcurr,gbase,gshell,insort,nshell
 isize=isize-4*mpw                          !shlim
 isize=isize+8*mpw*nspinor&
& +16*mpw*nspinor*(mpw*nspinor+1)       !eigval,eigvec
 memsize=max(memsize,isize)
 isize=isize+8*mpw+8*ngfft(4)&
& *ngfft(5)*ngfft(6)      !ts,vlocal
 memsize=max(memsize,isize)
 isize=isize+8*mgfft+4+28*mpw               !gbound,indpw_k,kg_k
 memsize=max(memsize,isize)
 isize=isize+8*natom&
& +24*mpw*ntypat*mpsang*mproj      !phkxred,ffnl,kinpw
 memsize=max(memsize,isize)
 isize=isize+16*mpw*natom                   !ph3d
 memsize=max(memsize,isize)
 isize=isize+48*mpw*nspinor&
& +8*mpw*nspinor*(mpw*nspinor+1)            !pwave,subghg,gvnlg
 if (nspinor==2)&
& isize=isize+40*mpw*nspinor                !pwave_so,subghg_so
 memsize=max(memsize,isize)
 isize=isize+8*mpw*nspinor*(mpw*nspinor+1)  !ghg
 memsize=max(memsize,isize)
 isize=isize+8*ngfft(4)*ngfft(5)*ngfft(6)   !work
 memsize=max(memsize,isize)
 isize=isize-8*ngfft(4)*ngfft(5)*ngfft(6)   !work
 isize=isize-8*mgfft+4+28*mpw&              !gbound,indpw_k,kg_k
&-8*natom-24*mpw*ntypat*mpsang*mproj&  !phkxred,ffnl,kinpw
&-16*mpw*natom                        !ph3d
 isize=isize-48*mpw*nspinor&
& -8*mpw*nspinor*(mpw*nspinor+1)        !pwave,subghg,gvnlg
 if (nspinor==2)&
& isize=isize-40*mpw*nspinor                !pwave_so,subghg_so

 isize=isize+56*mpw*nspinor                 !cwork,rwork
 memsize=max(memsize,isize)
 isize=isize-56*mpw*nspinor                 !cwork,rwork
 isize=isize+112*mpw*nspinor                !cwork,rwork,iwork,ifail
 memsize=max(memsize,isize)
 isize=isize-112*mpw*nspinor                !cwork,rwork,iwork,ifail
 isize=isize-8*mpw*nspinor*(mpw*nspinor+1)  !ghg
 isize=isize+8*mband                        !occ_k
 memsize=max(memsize,isize)
 isize=isize-8*mband                        !occ_k
 isize=isize-8*mpw*nspinor&
& -16*mpw*nspinor*(mpw*nspinor+1)       !eigval,eigvec
 isize=isize-32*mpw-8*ngfft(4)&
& *ngfft(5)*ngfft(6)     !gbig,ts,vlocal
 if(xmpi_paral==1)then
   isize=isize-4*4                           !kpt_distrb
 end if
 isize=isize-580-fnlen-4*(81+nkpt+9*nsym)-8*15   !non allocatable var.
!
 write(msg,'(2a,f8.2,a)')ch10,&
& ' Additional amount of memory required by "outkss" routine=',memsize*b2Mb,' Mbytes.'
 call wrtout(std_out,msg,'COLL')
!
end subroutine memkss
!!***

!!****f* ABINIT/dsksta
!! NAME
!! dsksta
!!
!! FUNCTION
!! This routine evaluates the amount of disk space required by the _KSS file.
!!
!! INPUTS
!!  dimlmn(natom*usepaw)=Number of nlm partial waves for each atom.
!!  ishm=Number of G-shells to be saved in _KSS file.
!!  mpsang=Max angular momentum +1 for pseudos.
!!  natom=Number of atoms in the unit cell.
!!  nbandkss=Number of desired bands to be saved in _KSS file
!!  nkpt=Number of k points.
!!  npwkss=Number of desired G-vectors to be saved in _KSS file.
!!  nspinor=Number of spinorial components.
!!  nsppol=Number of independent spin polarizations.
!!  ntypat=Number of type of atoms.
!!  nsym2=Number of symmetries in space group, without INV
!!  usepaw=1 if PAW.
!!
!! OUTPUT
!!  Writes on standard output
!!
!! PARENTS
!!      m_io_kss
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine dsksta(ishm,usepaw,nbandkss,mpsang,natom,ntypat,npwkss,nkpt,nspinor,nsppol,nsym2,dimlmn)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: usepaw,ishm,nbandkss,mpsang,natom,ntypat,nkpt
 integer,intent(in) :: npwkss,nspinor,nsppol,nsym2
!arrays
 integer,intent(in) :: dimlmn(natom*usepaw)

!Local variables-------------------------------
!scalars
 integer :: bsize_tot,bsize_hdr,bsize_kb,bsize_wf,bsize_cprj
 character(len=500) :: msg

! *********************************************************************

!The Abinit header is not considered.
 bsize_hdr= 80*2 + & !title
&5*4 + & !nsym2,nbandksseff,npwkss,ishm,mpsang
&nsym2*9*4 + & !symrel2
&nsym2*3*8 + & !tnons
&npwkss*3*4 + & !gbig
&ishm*4     !shlim

!NOTE: vkb does not depend on nsppol, however the elements are written for each spin.
 bsize_kb=0
 if (usepaw==0) then
   bsize_kb= nsppol* &
&   (         mpsang*ntypat       *8 + & !vkbsign
&  2*(nkpt*mpsang*ntypat*npwkss*8)  & !vkbd,vkbd
&  )
 end if

 bsize_wf= nsppol* &
& ( nkpt*nbandkss                 *8 + & !energies
&nkpt*nbandkss*nspinor*npwkss*2*8   & !wfg
&)

!For PAW add space required by projectors.
 bsize_cprj=0
 if (usepaw==1) then
   bsize_cprj=SUM(dimlmn(:))*(nsppol*nkpt*nspinor*nbandkss*2*8)
 end if

 bsize_tot = bsize_hdr + bsize_kb + bsize_wf + bsize_cprj
 write(msg,'(2a,f8.2,4a,4(a,f8.2,2a))')ch10,&
& ' Total amount of disk space required by _KSS file = ',bsize_tot*b2Mb,' Mb.',ch10,&
& '  Subdivided into : ',ch10,&
& '  Header             = ',bsize_hdr *b2Mb,' Mb.',ch10,&
& '  KB elements        = ',bsize_kb  *b2Mb,' Mb.',ch10,&
& '  Wavefunctions (PW) = ',bsize_wf  *b2Mb,' Mb.',ch10,&
& '  PAW projectors     = ',bsize_cprj*b2Mb,' Mb.',ch10
 call wrtout(std_out,msg,'COLL')

end subroutine dsksta
!!***

end subroutine outkss
!!***

END MODULE m_io_kss
!!***
