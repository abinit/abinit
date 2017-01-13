!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_io_kss
!! NAME
!!  m_io_kss
!!
!! FUNCTION
!!  This module contains procedured dealing with the IO of the KSS file.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (MG, GMR, VO, LR, RWG, MM, XG, RShaltaf)
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
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_wffile
 use m_xmpi
 use m_errors
 use m_abi_etsf
 use m_nctk
#ifdef HAVE_ETSF_IO
 use etsf_io
#endif
 use m_hdr
 use m_pawcprj

 use m_io_tools,         only : open_file
 use m_dtset,            only : dtset_copy, dtset_free
 use m_blas,             only : xdotc
 use m_mpinfo,           only : destroy_mpi_enreg
 use m_fftcore,          only : get_kg, sphere
 use m_crystal ,         only : crystal_t
 use m_crystal_io,       only : crystal_ncwrite
 use m_gsphere,          only : table_gbig2kg, merge_and_sort_kg
 use m_wfd,              only : wfd_t, wfd_ihave_ug, wfd_mybands, wfd_update_bkstab, &
&                               WFD_STORED, WFD_ALLOCATED, wfd_set_mpicomm

 implicit none

 private

 public :: testkss             ! Test a KSS file reporting basic quantities and dimensions.
 public :: write_kss_header    ! Writes the header of the KSS file.
 public :: read_kss_header     ! Read the head of the KSS file.
 public :: write_vkb           ! Writes the KB form factors and derivates on file for a single k-point.
 public :: write_kss_wfgk      ! Write the Gamma-centered wavefunctions and energies on the KSS file for a single k-point.
 public :: k2gamma_centered    ! Convert a set of wavefunctions from the k-centered to the gamma-centered basis set.
 public :: make_gvec_kss       ! Build the list of G-vectors for the KSS file.
 public :: gshgg_mkncwrite     ! Debugging tool used to build <G|H|G'> and dump data in netcdf format.

CONTAINS  !===========================================================
!!***

!!****f* m_io_kss/testkss
!! NAME
!! testkss
!!
!! FUNCTION
!! Test the KSS type.
!!
!! INPUTS
!!  iomode=Define the access mode (plain FORTRAN or NETCDF with ETSF-IO).
!!  filkss=Name of the KSS file. ".nc" will be appended in case of ETSF-IO.
!!  (TODO: remove this kind of a hack, using a module to store units and filenames
!!  Comm=MPI Communicator.
!!
!! OUTPUT
!!  mpsang=1+maximum angular momentum for nonlocal pseudopotential.
!!  nbnds_kss=Number of bands contained in the KSS file
!!  ng_kss=Number of plane waves in KSS file.
!!  nsym_out=Number of symmetries reported in the KSS, warning might differ from the symmetries found by abinit.
!!  Hdr<Hdr_type>=The abinit header.
!!
!! SIDE EFFECTS
!!  gvec_p(3,ng_kss)=
!!   In input : pointer to integers (not associated)
!!   In output: allocated array containing the G vectors reported in the KSS file.
!!  energies_p(nbnds_kss,Hdr%nkpt,Hdr%nsppol)=
!!   In input : pointer to real (not associated)
!!   In output: allocated array with the KS energies.
!!
!! NOTES
!!  Starting version 5.6, KSS files in single precision are not supported anymore.
!!
!! PARENTS
!!
!! CHILDREN
!!      metric,mkffnl,mkkin
!!
!! SOURCE

subroutine testkss(filkss,iomode,nsym_out,nbnds_kss,ng_kss,mpsang,gvec_p,energies_p,Hdr,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'testkss'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iomode,comm
 integer,intent(out) :: mpsang,nbnds_kss,ng_kss,nsym_out
 character(len=fnlen),intent(in) :: filkss
 type(Hdr_type),intent(inout) :: Hdr !vz_i
!arrays
 integer,pointer :: gvec_p(:,:)
 real(dp),pointer :: energies_p(:,:,:)

!Local variables-------------------------------
#ifdef HAVE_ETSF_IO
 type(ETSF_io_low_error) :: error
 type(ETSF_electrons),target :: Electrons_folder
 logical :: ok
 real(dp),allocatable,target :: eigen(:)
#endif
!scalars
 integer :: iatom,iband,ierr,ikpt,il,ispinor,isppol,itypat,master
 integer :: my_rank,kss_unt,nprocs,my_prtvol
 real(dp) :: nelect
 character(len=500) :: msg
 character(len=fnlen) :: fname
!arrays
 real(dp),allocatable :: tmp_enek(:)

! *************************************************************************

 DBG_ENTER("COLL")

 if (iomode==IO_MODE_ETSF) then
   write(msg,'(3a)')&
&   ' when iomode==3, support for the ETSF I/O library ',ch10,&
&   ' must be compiled. Use --enable-etsf-io when configuring '
#if !defined HAVE_ETSF_IO
   MSG_ERROR(msg)
#endif
 end if

 ! === MPI info ===
 my_rank = xmpi_comm_rank(comm)
 nprocs  = xmpi_comm_size(comm)
 master=0

 if (my_rank==master) then
!  === Read the header of the GS wavefunction file ===
!  TODO: remove this kind of a hack, using a module to store units and filenames.
   if (iomode==IO_MODE_FORTRAN) fname=filkss
   if (iomode==IO_MODE_ETSF   ) fname=nctk_ncify(filkss)

   my_prtvol=0
   call read_kss_header(kss_unt,fname,iomode,my_prtvol,nsym_out,nbnds_kss,ng_kss,mpsang,nelect,gvec_p,Hdr)

!  In case of parallelism over bands or Adler-Wiser with time reversal find the band
!  index separating the occupied and partially occupied from the empty states (for each spin)
!  Each processor will store in memory the occupied states while the conduction
!  states will be shared between different processors

!  TODO this part can be completely avoided if we read Hdr%occ.
!  The only problem is that Hdr%occ has to be correctly calculated in outkss in case of NSCF run.

   call wrtout(std_out,' testkss: reading occupation numbers ...','COLL')

!  NOTE : In old version of the code, the number of bands defined in the header was different
!  from the value reported in in the first section of a KSS file generated using kssform 3 (if nbandkss<nband).
!  NOW the BUG has been fixed but we keep these tests.
   ABI_CHECK(ALL(Hdr%nband==Hdr%nband(1)),'nband must be constant')
   ABI_CHECK(ALL(Hdr%nband==nbnds_kss),'nband must be equal to nbnds_kss')

   ABI_ALLOCATE(energies_p,(nbnds_kss,Hdr%nkpt,Hdr%nsppol))
   energies_p(:,:,:)=zero
   ABI_ALLOCATE(tmp_enek,(1:nbnds_kss))

   if (iomode==IO_MODE_FORTRAN) then  !  Read eigenvalues from the KSS file in the FORTRAN format

     do isppol=1,Hdr%nsppol
       do ikpt=1,Hdr%nkpt

         if (Hdr%usepaw==0) then
           do itypat=1,Hdr%ntypat
            do il=1,mpsang
              read(kss_unt) !vkbdb(:,itypat,il)
              read(kss_unt) !vkbdd(:,itypat,il)
            end do
           end do
         end if

         read(kss_unt) tmp_enek(1:nbnds_kss)
         energies_p(1:nbnds_kss,ikpt,isppol)=tmp_enek(1:nbnds_kss)

         do iband=1,nbnds_kss
           read(kss_unt) !kss_ugd(kss_npw*nspinor)
           if (Hdr%usepaw==1) then
             do ispinor=1,Hdr%nspinor
               do iatom=1,Hdr%natom
                read(kss_unt) !(cprjnk_k(iatom,ibsp)%cp(:,1:cprjnk_k(iatom,ibsp)%nlmn))
               end do
             end do
           end if
         end do

       end do !ikpt
     end do !isppol

#ifdef HAVE_ETSF_IO
   else if (iomode==IO_MODE_ETSF) then
     ABI_ALLOCATE(eigen,(nbnds_kss))
     eigen(:)=zero
!    allocate(occ_vec(nbnds_kss))
     do isppol=1,Hdr%nsppol
       do ikpt=1,Hdr%nkpt
!        NOTE : occupation numbers have been read from Hdr, and are recalculated in fermi.
         Electrons_folder%eigenvalues%data1D         => eigen
         Electrons_folder%eigenvalues__kpoint_access =  ikpt
         Electrons_folder%eigenvalues__spin_access   =  isppol
         call etsf_io_electrons_get(kss_unt,Electrons_folder,ok,error)
         ETSF_CHECK_ERROR(ok,error)

         energies_p(1:nbnds_kss,ikpt,isppol)=eigen(1:nbnds_kss)
         !write(std_out,*)isppol,ikpt,eigen(:)*Ha_eV
       end do
     end do
     nullify(Electrons_folder%eigenvalues%data1D)
     ABI_DEALLOCATE(eigen)
!    nullify(Electrons_folder%occupations%data1D)
!    deallocate(occ_vec)
#endif
   end if

   ABI_DEALLOCATE(tmp_enek)

   if (iomode==IO_MODE_FORTRAN) close(kss_unt)
   if (iomode==IO_MODE_ETSF) then
#ifdef HAVE_ETSF_IO
     NCF_CHECK(nf90_close(kss_unt))
#endif
   end if
 end if ! (my_rank==master)
 !
 !==========================================
 !=== Cast data if KSS file is not local ===
 !==========================================
 if (nprocs>1) then
   call xmpi_bcast(nsym_out, master,comm,ierr)
   call xmpi_bcast(nbnds_kss,master,comm,ierr)
   call xmpi_bcast(ng_kss,   master,comm,ierr)
   call xmpi_bcast(mpsang,   master,comm,ierr)
   call xmpi_bcast(nelect,   master,comm,ierr)
   call hdr_bcast(Hdr,master,my_rank,comm)
   if (my_rank/=master) then ! this proc did not read.
     ABI_ALLOCATE(gvec_p,(3,ng_kss))
     ABI_ALLOCATE(energies_p,(nbnds_kss,Hdr%nkpt,Hdr%nsppol))
   end if
   call xmpi_bcast(gvec_p,    master,comm,ierr)
   call xmpi_bcast(energies_p,master,comm,ierr)
   call xmpi_barrier(comm)
 end if

 DBG_EXIT("COLL")

end subroutine testkss
!!***

!----------------------------------------------------------------------

!!****f* m_io_kss/write_kss_header
!! NAME
!!  write_kss_header
!!
!! FUNCTION
!!  Write the header of the KSS file either using plain Fortran-IO or ETSF-IO.
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
!!  iomode=Input variables specifying the fileformat. (0-->Fortran,3-->ETSF-IO).
!!  occ(mband*nkpt*nsppol)=The occupation factors.
!!
!! OUTPUT
!!  kss_unt=The unit number of the opened file.
!!
!! SIDE EFFECTS
!!  The KSS Header is written on file.
!!
!! PARENTS
!!      outkss
!!
!! CHILDREN
!!      metric,mkffnl,mkkin
!!
!! SOURCE

subroutine write_kss_header(filekss,kss_npw,ishm,nbandksseff,mband,nsym2,symrel2,tnons2,occ,gbig,shlim,&
&  crystal,Dtset,Hdr,Psps,iomode,kss_unt)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'write_kss_header'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

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
#ifdef HAVE_ETSF_IO
 logical :: ok
 type(etsf_gwdata) :: GW_data
 type(wvl_wf_type) :: Dummy_wfs
 type(etsf_io_low_error) :: error
#endif
!arrays
 integer,allocatable,target :: vkbsign_int(:,:)
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
 ABI_DEALLOCATE(my_Hdr%occ)
 ABI_ALLOCATE(my_Hdr%occ,(my_Hdr%bantot))
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
 call dtset_copy(Dtset_cpy, Dtset)
 Dtset_cpy%mpw   = kss_npw
 Dtset_cpy%mband = nbandksseff

 fform=502

 SELECT CASE (iomode)

 CASE (IO_MODE_FORTRAN)

   if (open_file(filekss, msg, newunit=kss_unt, form="unformatted") /= 0) then
     MSG_ERROR(msg)
   end if 

   call hdr_fort_write(my_hdr, kss_unt, fform, ierr)
   ABI_CHECK(ierr == 0, "hdr_Fort_write returned ierr != 0")

   title='Results from ABINIT code';          write(kss_unt) title(1:80)
   title='Ab-initio plane waves calculation'; write(kss_unt) title(1:80)

   write(kss_unt) nsym2,nbandksseff,kss_npw,ishm,Psps%mpsang ! To be modified to deal with more than one projector
   write(kss_unt) (((symrel2(ii,jj,kk),ii=1,3),jj=1,3),kk=1,nsym2)
   write(kss_unt) ((tnons2(ii,kk),ii=1,3),kk=1,nsym2)
   write(kss_unt) ((gbig(ii,ig),ii=1,3),ig=1,kss_npw)
   write(kss_unt) (shlim(in),in=1,ishm)

   ! Write vkbsign for NC pseudos.
   ! FIXME : only one projector in each angular is treated.
   ! Moreover the allocation is done in the wrong order for dimensions...
   if (Psps%usepaw==0) then
     ABI_ALLOCATE(vkbsign,(Psps%ntypat,Psps%mpsang))
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
     ABI_DEALLOCATE(vkbsign)
   end if

#ifdef HAVE_ETSF_IO
 CASE (IO_MODE_ETSF)
!  We currently use the dataset symmetries, as defined in the Hdr structure instead of the symmetries recomputed in outkss.
   call abi_etsf_init(Dtset_cpy, filekss, 4, .false., my_Hdr%lmn_size, Psps, Dummy_wfs)

   ! Open again for further additions
   NCF_CHECK(nctk_open_modify(kss_unt, nctk_ncify(filekss), xmpi_comm_self))
   NCF_CHECK(crystal_ncwrite(crystal, kss_unt))

   ! Add additional info from abinit header.
   NCF_CHECK(hdr_ncwrite(my_hdr, kss_unt, fform, nc_define=.True.))

   ! If NC pseudos, write vkbsign. Here ordering of dimensions is OK but multi-projectors not supported.
   if (Psps%usepaw==0) then
     ABI_ALLOCATE(vkbsign_int,(Psps%mpsang,Psps%ntypat))
     vkbsign_int=0
     do itypat=1,Psps%ntypat
       il0=0
       do ilmn=1,Psps%lmnmax
         il=1+Psps%indlmn(1,ilmn,itypat)
         in=Psps%indlmn(3,ilmn,itypat)
         if ((il/=il0).and.(in==1)) then
           il0=il
           vkbsign_int(il,itypat)=NINT(DSIGN(one,Psps%ekb(ilmn,itypat)))
         end if
       end do
     end do
     ! Write it now to be able to deallocate quickly.
     GW_data%kb_formfactor_sign%data2D => vkbsign_int
     call etsf_io_gwdata_put(kss_unt, GW_data, ok, error)
     ETSF_CHECK_ERROR(ok,error)
     nullify(GW_data%kb_formfactor_sign%data2D)
     ABI_DEALLOCATE(vkbsign_int)
   end if
#endif

 CASE DEFAULT
   write(msg,'(a,i0)')" Unsupported value for iomode: ",iomode
   MSG_ERROR(msg)
 END SELECT

 call dtset_free(Dtset_cpy)
 call hdr_free(my_Hdr)

 DBG_EXIT("COLL")

end subroutine write_kss_header
!!***

!----------------------------------------------------------------------

!!****f* m_io_kss/read_kss_header
!! NAME
!!  read_kss_header
!!
!! FUNCTION
!!  Read the header of the KSS file either using plain Fortran-IO or ETSF-IO.
!!
!! INPUTS
!!  filkss(len=fnlen)=The name of the KSS file.
!!  iomode=Input variables specifying the fileformat. (0-->Fortran,3-->ETSF-IO).
!!  prtvol=Flag governing verbosity output.
!!
!! OUTPUT
!!  kss_unt=The unit number of the opened file.
!!  nsym_out=Number of symmetry operations read from the KSS file (not from the abinit header)
!!  nbnds_kss=Number of bands stored.
!!  ng_kss=Number of planewaves used for the wavefunctions in the KSS files.
!!  mpsang=Max angular momentum + 1
!!  nelect=Number of electrons (including a possible charge in the unit cell)
!!  Hdr<hdr_type>=The abinit header.
!!
!! SIDE EFFECTS
!!  gvec_p(3,ng_kss)=Input:  Pointer to null()
!!                   Output: The set of G-vectors for the KSS wavefunctions (Gamma-centered)
!!
!! PARENTS
!!      m_io_kss
!!
!! CHILDREN
!!      metric,mkffnl,mkkin
!!
!! SOURCE

subroutine read_kss_header(kss_unt,filkss,iomode,prtvol,nsym_out,nbnds_kss,ng_kss,mpsang,nelect,gvec_p,Hdr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'read_kss_header'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iomode,prtvol
 integer,intent(out) :: nsym_out,nbnds_kss,ng_kss,mpsang,kss_unt
 real(dp),intent(out) :: nelect
 character(len=fnlen),intent(in) :: filkss
 type(Hdr_type),intent(inout) :: Hdr !vz_i
!arrays
 integer,pointer :: gvec_p(:,:)

!Local variables-------------------------------
!scalars
 integer :: nshells
 integer :: fform,ii,ig,nsym_kss
 character(len=80) :: title(2)
 character(len=500) :: msg
 logical :: ltest
#ifdef HAVE_ETSF_IO
 logical :: ok
 type(etsf_io_low_error) :: error
 type(ETSF_dims) :: Dims
 type(ETSF_basisdata),target :: Wave_folder
!arrays
 integer,allocatable,target :: kg_k(:,:)
#endif

! *************************************************************************

 DBG_ENTER("COLL")

#if !defined HAVE_ETSF_IO
 if (iomode==IO_MODE_ETSF) then
   write(msg,'(3a)')&
&   ' When iomode==3, support for the ETSF I/O library ',ch10,&
&   ' must be compiled. Use --enable-etsf-io when configuring '
   MSG_ERROR(msg)
 end if
#endif

 ! ===============================================
 ! ==== Read file according to the fileformat ====
 ! ===============================================
 SELECT CASE (iomode)

 CASE (IO_MODE_FORTRAN) ! * Formatted Fortran File
   write(msg,'(2a)')' reading Fortran Kohn-Sham structure file ',TRIM(filkss)
   call wrtout(std_out,msg,'COLL')

   if (open_file(filkss, msg, newunit=kss_unt, form='unformatted', status='old') /= 0) then
     MSG_ERROR(msg)
   end if

   call hdr_fort_read(hdr, kss_unt, fform)
   ABI_CHECK(fform /= 0, "hdr_fort_read returned fform == 0")

 CASE (IO_MODE_ETSF) ! * NETCDF-ETSF file format

#ifdef HAVE_ETSF_IO
   write(msg,'(2a)')' read_kss_header : reading NETCDF-ETSF Kohn-Sham structure file ',TRIM(filkss)
   call wrtout(std_out,msg,'COLL')

   NCF_CHECK(nctk_open_read(kss_unt, filkss, xmpi_comm_self))
   call hdr_ncread(Hdr, kss_unt, fform)

   ABI_CHECK(fform /= 0, "hdr_ncread returned fform == 0")
   ABI_CHECK(fform/=602,' Single precision KSS + ETSF-IO not implemented')
   ABI_CHECK(Hdr%usepaw==0,'PAW+ETSF-IO not yet coded')
#endif

 CASE DEFAULT
   write(msg,'(a,i4)')'Wrong value for iomode = ',iomode
   MSG_ERROR(msg)
 END SELECT

 if (fform==602) then ! * fform must be 502.
   write(msg,'(3a)')&
&   ' Starting v5.6, KSS files in single precision are not supported anymore,',ch10,&
&   ' Please, use an older version of abinit.'
   MSG_ERROR(msg)
 end if
 if (fform>=1.and.fform<=2) then
   write(msg,'(a,i4)')' (STA|QPLDA) format are obsolete and not supported anymore; fform= ',fform
   MSG_ERROR(msg)
 end if
 if (fform/=502) then
   write(msg,'(a,i4)')' Found unknown file format; fform= ',fform
   MSG_ERROR(msg)
 end if
 !
 ! === Output the header of the GS wavefunction file ===
 if (prtvol>0) call hdr_echo(hdr, fform, 4, unit=std_out)

 write(msg,'(1x,47a)')('-',ii=1,47)
 call wrtout(std_out,msg,'COLL')
 write(msg,'(3a,a6,a,i3)')&
&  ' KSS abinit double precision form',ch10,&
&  ' generated by ABINIT ',Hdr%codvsn,' header version ',Hdr%headform
 call wrtout(std_out,msg,'COLL')

! === Test spin-orbit characteristic ===
 if (Hdr%headform<53) then  ! Format previous to version 5.5, now pspo is obsolete and has been substituted by so_psp.
   ltest=ALL(Hdr%pspso(1:Hdr%ntypat)==1)
   ABI_CHECK(ltest,'pspso/=1 value not programmed')
 else  ! New format containing so_psp
   ltest=ALL(Hdr%so_psp(1:Hdr%npsp)==1)
   ABI_CHECK(ltest,'so_psp/=1 value not programmed')
 end if

 ! There might be extra charge thus we use occ_out. factors to calculate nelect.
 ! Besides note that nelect is real
 !nelect = hdr_get_nelect_byocc(Hdr)
 nelect = Hdr%nelect
 !
 ! === Abinit Header successfully read ===
 ! * Now read basic dimensions.
 ! * Note that, in the case of Fortran file, nsym_out is read from the second record
 nsym_out=Hdr%nsym

 SELECT CASE (iomode)

 CASE (IO_MODE_FORTRAN)
   read(kss_unt) title(1)
   read(kss_unt) title(2)
   write(msg,'(2a,1x,a79,a,1x,a79,a)')' title of file: ',ch10,title(1)(:79),ch10,title(2)(:79),ch10
   call wrtout(std_out,msg,'COLL')
   read(kss_unt) nsym_kss,nbnds_kss,ng_kss,nshells,mpsang
   read(kss_unt) !(((symrel2(jj,ii,isym),ii=1,3),jj=1,3),isym=1,nsym_kss)
   read(kss_unt) !((tnons(i,isym),i=1,3),isym=1,nsym_kss)

   ABI_ALLOCATE(gvec_p,(3,ng_kss))
   read(kss_unt)((gvec_p(ii,ig),ii=1,3),ig=1,ng_kss)
   nsym_out=nsym_kss

   read(kss_unt)                      !(shlim(i),i=1,nshells)
   if (Hdr%usepaw==0) read(kss_unt)   !((vkbsignd(il,is),il=1,mpsang),is=1,Hdr%ntypat)

 CASE (IO_MODE_ETSF) ! TODO spin-orbit not treated, number of projectors not treated
#ifdef HAVE_ETSF_IO
   call etsf_io_dims_get(kss_unt,Dims,ok,error)
   nsym_kss =Dims%number_of_symmetry_operations
   nbnds_kss=Dims%max_number_of_states
   ng_kss   =Dims%max_number_of_coefficients
   mpsang   =Dims%max_number_of_angular_momenta

   ABI_ALLOCATE(gvec_p,(3,ng_kss))
   ABI_ALLOCATE(kg_k,(3,ng_kss))
   Wave_folder%reduced_coordinates_of_plane_waves%data2D => kg_k(:,:)
   call etsf_io_basisdata_get(kss_unt,Wave_folder,ok,error)
   gvec_p(:,:)=kg_k(:,:)
   ABI_DEALLOCATE(kg_k)
   nshells=0 ! nshells is not defined in the ETSF spefications but it is not used
#endif

 CASE DEFAULT
   MSG_ERROR("Unsupported value for iomode")
 END SELECT

 ABI_CHECK(ALL(gvec_p(1:3,1)==0),'First G must be 0')

 if (prtvol>0) then ! Output important dimensions on the log file.
   write(msg,'(a,f8.2)')' number of electrons                    ',nelect
   call wrtout(std_out,msg,'COLL')
   write(msg,'(a,i8)')' number of symmetries without inversion ',nsym_out
   call wrtout(std_out,msg,'COLL')
   write(msg,'(a,i8)')' number of bands                        ',nbnds_kss
   call wrtout(std_out,msg,'COLL')
   write(msg,'(a,i8)')' number of plane waves                  ',ng_kss
   call wrtout(std_out,msg,'COLL')
   write(msg,'(a,i8)')' number of shells                       ',nshells
   call wrtout(std_out,msg,'COLL')
   write(msg,'(a,i8,a)')' maximum angular momentum +1          ',mpsang,ch10
   call wrtout(std_out,msg,'COLL')
   write(msg,'(1x,47a)')('-',ii=1,47)
   call wrtout(std_out,msg,'COLL')
 end if
 !
 ! === Check the value of some variables ===
 ! This is due to the awful treatment of symmetries done in outkss
 if (Hdr%nsym/=nsym_kss) then
   write(msg,'(2a,2(a,i3))')&
&    ' Code does not use the original set of symmetries.',ch10,&
&    ' Hdr%nsym= ',Hdr%nsym,' /= nsym_kss= ',nsym_kss
   MSG_COMMENT(msg)
 end if

 DBG_EXIT("COLL")

end subroutine read_kss_header
!!***

!----------------------------------------------------------------------

!!****f* m_io_kss/write_vkb
!! NAME
!!  write_vkb
!!
!! FUNCTION
!!  Writes the KB form factors and derivates on file for a single k-point.
!!  Supports plain Fortran IO and ETSF-IO.
!!
!! INPUTS
!!  kss_unt=The unit number of the file
!!  ikpt=The index of the k-point
!!  kpoint(3)=The k-point in reduced coordinates.
!!  kss_npw=Number of planewaves used for the wavefunctions in the KSS files.
!!  npw_k=Number of planewaves at this k-point in the k-centered basis set used in abinit (ecut).
!!  trsl(kss_npw)=Mapping between the G-sphere used for the KSS wavefunctions and the
!!   Abinit G-sphere (npw_k). As kss_npw>=npw_k, trsl=npw_k+1 is the "KSS" G-vector
!!   is not contained in the abinit one.
!!  ecut=cutoff energy used in abinit.
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr).
!!  Psps<Pseudopotential_type>=Datatype gathering data on the Pseudopotentials.
!!  iomode=Input variables specifying the fileformat. (0-->Fortran,3-->ETSF-IO).
!!  gbig(3,kss_npw)=Set of G-vectors used in the KSS file.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      m_io_kss
!!
!! CHILDREN
!!      metric,mkffnl,mkkin
!!
!! SOURCE

subroutine write_vkb(kss_unt,ikpt,kpoint,kss_npw,gbig,trsl,rprimd,Psps,iomode)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'write_vkb'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikpt,iomode,kss_npw,kss_unt
 type(Pseudopotential_type),intent(in) :: Psps
!arrays
 integer,intent(in) :: trsl(kss_npw),gbig(3,kss_npw)
 real(dp),intent(in) :: kpoint(3),rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: itypat,il,ig,mpsang,ntypat
 character(len=500) :: msg
!array
 real(dp),allocatable :: vkb(:,:,:),vkbd(:,:,:)
 real(dp),allocatable :: dum_vkbsign(:,:)
#ifdef HAVE_ETSF_IO
 logical :: ok
 real(dp),allocatable,target :: vkb_tgt(:,:,:), vkbd_tgt(:,:,:)
 type(etsf_io_low_error) :: error
 type(etsf_gwdata) :: GW_data
#endif

! *********************************************************************

 mpsang = Psps%mpsang
 ntypat = Psps%ntypat

 ABI_ALLOCATE(vkb ,(kss_npw,ntypat,mpsang))
 ABI_ALLOCATE(vkbd,(kss_npw,ntypat,mpsang))
 ABI_ALLOCATE(dum_vkbsign,(ntypat,mpsang))

 call kss_calc_vkb(Psps,kpoint,kss_npw,gbig,rprimd,dum_vkbsign,vkb,vkbd)
 ABI_DEALLOCATE(dum_vkbsign)

 SELECT CASE (iomode)

 CASE (IO_MODE_FORTRAN)
  do itypat=1,ntypat
    do il=1,mpsang
      !write(kss_unt) (vkb (trsl(ig),itypat,il),ig=1,kss_npw)
      !write(kss_unt) (vkbd(trsl(ig),itypat,il),ig=1,kss_npw)
      write(kss_unt) (vkb (ig,itypat,il),ig=1,kss_npw)
      write(kss_unt) (vkbd(ig,itypat,il),ig=1,kss_npw)
    end do
  end do

#ifdef HAVE_ETSF_IO
 CASE (IO_MODE_ETSF)
   ABI_ALLOCATE(vkb_tgt ,(kss_npw,mpsang,ntypat))
   ABI_ALLOCATE(vkbd_tgt,(kss_npw,mpsang,ntypat))
   do itypat=1,ntypat
     do il=1,mpsang
       do ig=1,kss_npw
         !vkb_tgt (ig,il,itypat)=vkb (trsl(ig),itypat,il)
         !vkbd_tgt(ig,il,itypat)=vkbd(trsl(ig),itypat,il)
         vkb_tgt (ig,il,itypat)=vkb (ig,itypat,il)
         vkbd_tgt(ig,il,itypat)=vkbd(ig,itypat,il)
       end do
     end do
   end do
   GW_data%kb_coeff__kpoint_access         = ikpt
   GW_data%kb_coeff_der__kpoint_access     = ikpt
   GW_data%kb_formfactors%data3D           => vkb_tgt
   GW_data%kb_formfactor_derivative%data3D => vkbd_tgt

   call etsf_io_gwdata_put(kss_unt, GW_data, ok, error)
   ETSF_CHECK_ERROR(ok,error)

   nullify(GW_data%kb_formfactors%data3D          ) ! Avoid dangling pointers
   nullify(GW_data%kb_formfactor_derivative%data3D)
   ABI_DEALLOCATE(vkb_tgt)
   ABI_DEALLOCATE(vkbd_tgt)
#endif

 CASE DEFAULT
   write(msg,'(a,i0)')" Unsupported value for iomode: ",iomode
   MSG_ERROR(msg)
 END SELECT

 ABI_DEALLOCATE(vkb)
 ABI_DEALLOCATE(vkbd)

 RETURN
 ABI_UNUSED(trsl(1)) ! just to keep trsl as an argument while in development

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
!!  kg_k(3,npw_k)=The G-vectors in the k-centered basis set.
!!  gbig(3,kss_npw)=The set of G-vectors for the KSS wavefunctions (Gamma-centered)
!!  wfg(2,kss_npw*nspinor,nbandksseff)=The wavefunction Fourier coefficients.
!!  iomode=Input variables specifying the fileformat. (0-->Fortran,3-->ETSF-IO).
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      outkss
!!
!! CHILDREN
!!      metric,mkffnl,mkkin
!!
!! SOURCE

subroutine write_kss_wfgk(kss_unt,ikpt,isppol,kpoint,nspinor,kss_npw,npw_k,kg_k,&
&          nbandksseff,natom,Psps,ene_k,occ_k,rprimd,gbig,wfg,Cprjnk_k,iomode)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'write_kss_wfgk'
 use interfaces_51_manage_mpi
 use interfaces_56_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikpt,isppol,iomode,kss_npw,nspinor,kss_unt,nbandksseff
 integer,intent(in) :: natom,npw_k
 type(pseudopotential_type),intent(in) :: Psps
!arrays
 integer,intent(in) :: gbig(3,kss_npw),kg_k(3,npw_k)
 real(dp),intent(in) :: kpoint(3),rprimd(3,3)
 real(dp),intent(in) :: ene_k(nbandksseff),occ_k(nbandksseff)
 real(dp),intent(in) ::  wfg(2,kss_npw*nspinor,nbandksseff)
 type(pawcprj_type),intent(in) :: Cprjnk_k(natom,nspinor*nbandksseff*Psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: ib,ibsp,ig,ispinor,iatom,ii,master,my_rank,ierr
 type(MPI_type) :: MPI_enreg_seq
 character(len=500) :: msg
#ifdef HAVE_ETSF_IO
 type(wffile_type) :: Wff
#endif
!arrays
 integer,allocatable :: trsl(:)

! *********************************************************************

!* Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)
 master=0; my_rank=master

 if (Psps%usepaw==0) then ! Calculate and write KB form factors and derivative at this k-point.
  ! The array trsl translates the index of gbig into the corresponding
  ! index in array kg_k. If gbig(ig) does not exist in kg_k, trsl(ig) is set to npw_k+1.
  ABI_ALLOCATE(trsl,(kss_npw))
  call table_gbig2kg(npw_k,kg_k,kss_npw,gbig,trsl,ierr)
  if (ierr/=0.and.(kss_npw>=npw_k)) then
    MSG_ERROR(' The set of g vectors is inconsistent ! Check source.')
  end if
  call write_vkb(kss_unt,ikpt,kpoint,kss_npw,gbig,trsl,rprimd,Psps,iomode)
  ABI_DEALLOCATE(trsl)
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

#ifdef HAVE_ETSF_IO
 CASE (IO_MODE_ETSF) ! When ETSF, use the rwwf routine (should always be done like that)
   if (Psps%usepaw==1) then
     MSG_WARNING("PAW output with ETSF-IO: cprj won't be written")
   end if
   Wff%master   =master
   Wff%me       =my_rank
   Wff%unwff    =kss_unt
   Wff%iomode=IO_MODE_ETSF
!  MG Tue Oct 23 occ_k is passed to writewf instead of occ to treat correctly metals
!  MG FIXME The RESHAPE below is ugly. Some compiler will likely create a buffer on the stack and then BOOM!
  ! call writewf(RESHAPE(wfg(:, 1:kss_npw, :),(/2,nbandksseff*kss_npw /)),ene_k,0,&
   call writewf(wfg(1,1,1),ene_k,0,&
&   0,ikpt,isppol,gbig(:,1:kss_npw),nbandksseff,nbandksseff*kss_npw,MPI_enreg_seq,&
&   nbandksseff,nbandksseff,kss_npw,nspinor,occ_k(1:nbandksseff),2,1,Wff)
#endif

 CASE DEFAULT
   write(msg,'(a,i0)')" Unsupported iomode: ",iomode
   MSG_ERROR(msg)
 END SELECT

 call destroy_mpi_enreg(MPI_enreg_seq)

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
!!      outkss
!!
!! CHILDREN
!!      metric,mkffnl,mkkin
!!
!! SOURCE

subroutine k2gamma_centered(kpoint,npw_k,istwf_k,ecut,kg_k,kss_npw,nspinor,nbandksseff,ngfft,gmet,&
&  MPI_enreg,gbig,ug,icg,cg,eig_vec)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'k2gamma_centered'
!End of the abilint section

 implicit none

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
   MSG_ERROR("Both cg and eig_vec are present!")
 end if

! Mapping between the gamma-centered basis set and the k-centered one.
! trsl(ig)=npw_k+1 if vector ig is not inside the k-centered G-sphere.
 ABI_ALLOCATE(trsl,(kss_npw))

 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 n4=ngfft(4); n5=ngfft(5); n6=ngfft(6)

 if (istwf_k==1) then ! Full k-centered G-sphere.
   call table_gbig2kg(npw_k,kg_k,kss_npw,gbig,trsl,ierr)
   if (ierr/=0.and.(kss_npw>=npw_k)) then
     MSG_ERROR(' The set of G vectors is inconsistent')
   end if

 else  ! Calculate full kg with istwf_k=1 then do the mapping.
   call get_kg(kpoint,1,ecut,gmet,full_npw_k,full_kg_k)

   call table_gbig2kg(full_npw_k,full_kg_k,kss_npw,gbig,trsl,ierr)
   if (ierr/=0.and.(kss_npw>=npw_k)) then
     MSG_ERROR(' The set of G vectors is inconsistent')
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
     ABI_ALLOCATE(cfft,(2,n4,n5,n6*ndat))
     ABI_ALLOCATE(full_cg,(2,full_npw_k*ndat))
     ABI_ALLOCATE(tmp_cg,(2,npw_k*ndat))

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

     ABI_DEALLOCATE(cfft)
     ABI_DEALLOCATE(tmp_cg)
     ABI_DEALLOCATE(full_cg)

   CASE DEFAULT
     MSG_BUG("Wrong istwf_k")
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
     MSG_ERROR(msg)
   END SELECT

 else
   MSG_ERROR("neither cg not eig_vec are in input")
 end if

 ABI_DEALLOCATE(trsl)
 if (allocated(full_kg_k))  then
   ABI_DEALLOCATE(full_kg_k)
 end if

end subroutine k2gamma_centered
!!***

!----------------------------------------------------------------------

!!****f* m_io_kss/skip_kss_record
!! NAME
!!  skip_kss_record
!!
!! FUNCTION
!!  Skip one or more wavefunction records of the KSS file.
!!
!! INPUTS
!!  kss_unt=Unit of the KSS file.
!!  nrec=Number of records to be skipped.
!!  usepaw=1 if PAW
!!  nspinor=Number of spinor components.
!!  natom=Nuber of atom
!!
!!  OUTPUT
!!   Error status reported by Fortran read statement.
!!
!! SIDE EFFECTS
!!  See description.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function skip_kss_record(kss_unt,nrec,usepaw,nspinor,natom) result(ios)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'skip_kss_record'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: kss_unt,usepaw,nspinor,natom,nrec
 integer :: ios

!Local variables-------------------------------
 integer :: ispinor,iatom,irec

! *************************************************************************

 ios=0
 do irec=1,nrec
   read(kss_unt,err=10) ! kss_ugd(1:kss_npw*nspinor)
   if (usepaw==1) then
     do ispinor=1,nspinor
       do iatom=1,natom
         read(kss_unt,err=10) !(Cprj_ibz(iatom,ibsp)%cp
       end do
     end do
   end if
 end do

 RETURN

 10 ios=-1

end function skip_kss_record
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
!!      gwls_hamiltonian,setup_screening,setup_sigma
!!
!! CHILDREN
!!      metric,mkffnl,mkkin
!!
!! SOURCE

subroutine make_gvec_kss(nkpt,kptns,ecut_eff,symmorphi,nsym,symrel,tnons,gprimd,prtvol,npwkss,gvec_kss,ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_gvec_kss'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
!End of the abilint section

 implicit none

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
     MSG_WARNING(msg)
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
   MSG_WARNING(msg)
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

!!****f* ABINIT/gshgg_mkncwrite
!! NAME
!! gshgg_mkncwrite
!!
!! FUNCTION
!!  This routine builds <G|H|G'> matrix elements for all (G, G').
!!  starting from the knowledge of the local potential on the real-space FFT mesh.
!!  It can also perform the direct diagonalization of the Kohn-Sham Hamiltonian
!!  for a given k-point and spin. This a debugging tool
!!
!! INPUTS
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  kpoint(3)
!!  natom=number of atoms in cell.
!!  nfftf=(effective) number of FFT grid points in the dense FFT mesh (for this processor)
!!         (nfftf=nfft for norm-conserving potential runs)
!!  nspinor=number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in unit cell.
!!  pawtab(psps%ntypat*psps%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  pawfgr<pawfgr_type>=fine grid parameters and related data
!!  Paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  vtrial(nfftf,nspden)=the trial potential
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  eigen(mband,nkpt,nsppol)=array for holding eigenvalues (hartree)
!!  ngfftc(18)=Info about 3D FFT for the coarse mesh, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  [Electronpositron] <electronpositron_type>=quantities for the electron-positron annihilation.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      metric,mkffnl,mkkin
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine gshgg_mkncwrite(istep, dtset, dtfil, psps, hdr, pawtab, pawfgr, paw_ij, mpi_enreg, &
  rprimd, xred, eigen, npwarr, kg, ylm, ngfftc, nfftc, ngfftf, nfftf, vtrial,&
  electronpositron) ! Optional arguments

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_hamiltonian
 use m_nctk
 use m_wfk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_electronpositron
 use m_pawtab
 use m_paw_ij
 use m_pawcprj
 use m_pawfgr

 use m_fstrings,    only : strcat
 use m_abilasi,     only : xheevx, xhegvx

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gshgg_mkncwrite'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_65_paw
 use interfaces_66_nonlocal
 use interfaces_66_wfs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep,nfftc,nfftf
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(pseudopotential_type),intent(in) :: psps
 type(pawfgr_type),intent(in) :: pawfgr
 type(hdr_type),intent(in) :: hdr
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfftc(18),ngfftf(18)
 integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem),npwarr(dtset%nkpt)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: eigen(dtset%mband,dtset%nkpt,dtset%nsppol)
 real(dp),intent(inout) :: vtrial(nfftf,dtset%nspden)
 real(dp),intent(in) :: xred(3,dtset%natom)
 real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(paw_ij_type),intent(in) :: Paw_ij(dtset%natom*psps%usepaw)
 type(electronpositron_type),optional,pointer :: electronpositron

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_getghc=4,master=0,fform=102
 integer :: comm,cpopt,dimffnl,ib,ider,idir,isppol,npw_k,npws,ierr
 integer :: ikg,istwf_k,ntypat,paral_kgb,ikpt,ilm,natom,nband_k,ncerr
 integer :: jj,n1,n2,n3,n4,n5,n6,negv,nkpg,nprocs,nspden,nspinor,mgfftc,ncid
 integer :: my_rank,ispden,ndat,type_calc,sij_opt,igsp2,cplex_ghg
 integer :: kg_varid,ghg_varid,pawgtg_varid,eighist_varid,eigdiago_varid
 real(dp) :: cfact,ucvol,lambda,size_mat
 character(len=50) :: jobz,range
 character(len=80) :: frmt1,frmt2
 character(len=10) :: stag(2)
 character(len=500) :: msg
 character(len=fnlen) :: path
 type(gs_hamiltonian_type) :: gs_hamk
!arrays
 integer,allocatable :: kg_k(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rhodum(1),kpoint(3),ylmgr_dum(1,1,1)
 real(dp),allocatable :: eig_ene(:),eig_vec(:,:,:),ph3d(:,:,:),pwave(:,:)
 real(dp),allocatable :: ffnl(:,:,:,:),kinpw(:),kpg_k(:,:)
 real(dp),allocatable :: vlocal(:,:,:,:),ylm_k(:,:),vlocal_tmp(:,:,:)
 real(dp),allocatable :: ghc(:,:),gvnlc(:,:),gsc(:,:),ghg_mat(:,:,:),gtg_mat(:,:,:),cgrvtrial(:,:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)

! *********************************************************************

 DBG_ENTER("COLL")
 ABI_UNUSED(ngfftf(1))

 call wrtout(std_out, "Constructing H_{GG')(k,spin) and dumping matrices to file")

 comm = mpi_enreg%comm_cell; nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 paral_kgb = mpi_enreg%paral_kgb
 ABI_CHECK(paral_kgb == 0 ,"paral_kgb == 1 not coded")
 ABI_CHECK(all(dtset%nband == dtset%mband) ,"nband must be constant")

 natom = dtset%natom; ntypat = dtset%ntypat
 nspden = dtset%nspden; nspinor = dtset%nspinor
 if (dtset%nsppol==1) stag=['          ','          ']
 if (dtset%nsppol==2) stag=['SPIN UP:  ','SPIN DOWN:']

 if (nprocs > 1 .and. .not. nctk_has_mpiio) then
   MSG_ERROR("Netcdf without MPI-IO support. Cannot produce HGG.nc file with nprocs > 1")
 end if

 path = strcat(dtfil%filnam_ds(4), "_HGG.nc")
#ifdef HAVE_ETSF_IO
 if (istep == 1) then
   ! Open netcdf file, define dims and variables.
   NCF_CHECK(nctk_open_create(ncid, path, comm))

   ncerr = nctk_def_dims(ncid, [&
     nctkdim_t("complex", 2), nctkdim_t("max_number_of_coefficients", dtset%mpw),&
     nctkdim_t("max_number_of_states", dtset%mband),&
     nctkdim_t("number_of_kpoints", dtset%nkpt), nctkdim_t("number_of_spins", dtset%nsppol), &
     nctkdim_t("nstep", dtset%nstep) ])
   NCF_CHECK(ncerr)

   call wfk_ncdef_dims_vars(ncid, hdr, fform)

   NCF_CHECK(nctk_def_iscalars(ncid, ["last_step"]))

   ncerr = nctk_def_arrays(ncid, nctkarr_t('eigenvalues_history', "dp", &
     "max_number_of_states, number_of_kpoints, number_of_spins, nstep"))
   NCF_CHECK(ncerr)

   ncerr = nctk_def_arrays(ncid, nctkarr_t('eigenvalues_diago', "dp", &
     "max_number_of_states, number_of_kpoints, number_of_spins, nstep"))
   NCF_CHECK(ncerr)

   ncerr = nctk_def_arrays(ncid, nctkarr_t('ghg', "dp", &
     "complex, max_number_of_coefficients, max_number_of_coefficients, number_of_kpoints, number_of_spins, nstep"))
   NCF_CHECK(ncerr)

   if (psps%usepaw == 1) then
     ncerr = nctk_def_arrays(ncid, nctkarr_t('paw_gtg', "dp", &
       "complex, max_number_of_coefficients, max_number_of_coefficients, number_of_kpoints, number_of_spins, nstep"))
     NCF_CHECK(ncerr)
   end if
 else
   NCF_CHECK(nctk_open_modify(ncid, path, comm))
 end if

 ! Use individual IO (default) for the G-vectors [3, mpw, nkpt]
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_inq_varid(ncid, "reduced_coordinates_of_plane_waves", kg_varid))
 NCF_CHECK(nf90_inq_varid(ncid, "eigenvalues_diago", eigdiago_varid))
 NCF_CHECK(nf90_inq_varid(ncid, "eigenvalues_history", eighist_varid))
 NCF_CHECK(nf90_inq_varid(ncid, "ghg", ghg_varid))
 if (psps%usepaw == 1) then
   NCF_CHECK(nf90_inq_varid(ncid, "paw_gtg", pawgtg_varid))
 end if

 NCF_CHECK(nctk_write_iscalars(ncid, ["last_step"], [istep]))
#endif

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 ! The coarse FFT mesh.
 mgfftc = maxval(ngfftc(:3))
 n1=ngfftc(1); n2=ngfftc(2); n3=ngfftc(3)
 n4=ngfftc(4); n5=ngfftc(5); n6=ngfftc(6)

!============================================
!==== Initialize most of the Hamiltonian ====
!============================================
!1) Allocate all arrays and initialize quantities that do not depend on k and spin.
!2) Perform the setup needed for the non-local factors:
!* Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
!* PAW: Initialize the overlap coefficients and allocate the Dij coefficients.

 if (PRESENT(Electronpositron)) then
   call init_hamiltonian(gs_hamk,psps,pawtab,nspinor,dtset%nsppol,nspden,natom,dtset%typat,xred,nfftc,&
&   mgfftc,ngfftc,rprimd,dtset%nloalg,paw_ij=paw_ij,usecprj=0,electronpositron=electronpositron,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab)
 else
   call init_hamiltonian(gs_hamk,psps,pawtab,nspinor,dtset%nsppol,nspden,natom,dtset%typat,xred,nfftc,&
&   mgfftc,ngfftc,rprimd,dtset%nloalg,paw_ij=paw_ij,usecprj=0,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab)
 end if

 ! Set up local potential vlocal with proper dimensioning, from vtrial.
 ! Select spin component of interest if nspden<=2 as nvloc==1, for nspden==4, nvloc==4
 ! option=2: vtrial(n1*n2*n3,ispden) --> vlocal(nd1,nd2,nd3) real case
 ABI_MALLOC(vlocal,(n4,n5,n6,gs_hamk%nvloc))
 ABI_MALLOC(kg_k,(3,dtset%mpw))

 do isppol=1,dtset%nsppol
   ikg = 0

   ! Set up local potential vlocal with proper dimensioning, from vtrial
   ! Also take into account the spin.
   if (nspden/=4)then
     if (psps%usepaw==0 .or. pawfgr%usefinegrid==0) then
       call fftpac(isppol,mpi_enreg,nspden,n1,n2,n3,n4,n5,n6,ngfftc,vtrial,vlocal,2)
     else
       ! Move from fine to coarse FFT mesh (PAW)
       ABI_MALLOC(cgrvtrial,(nfftc,nspden))
       call transgrid(1,mpi_enreg,nspden,-1,0,0,paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)
       call fftpac(isppol,mpi_enreg,nspden,n1,n2,n3,n4,n5,n6,ngfftc,cgrvtrial,vlocal,2)
       ABI_FREE(cgrvtrial)
     end if
   else
     ABI_MALLOC(vlocal_tmp,(n4,n5,n6))
     if (psps%usepaw==0 .or. pawfgr%usefinegrid==0) then
       do ispden=1,nspden
         call fftpac(ispden,mpi_enreg,nspden,n1,n2,n3,n4,n5,n6,ngfftc,vtrial,vlocal_tmp,2)
         vlocal(:,:,:,ispden)=vlocal_tmp(:,:,:)
       end do
     else
       ! Move from fine to coarse FFT mesh (PAW)
       ABI_MALLOC(cgrvtrial,(nfftc,nspden))
       call transgrid(1,mpi_enreg,nspden,-1,0,0,paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)
       do ispden=1,nspden
         call fftpac(ispden,mpi_enreg,nspden,n1,n2,n3,n4,n5,n6,ngfftc,cgrvtrial,vlocal_tmp,2)
         vlocal(:,:,:,ispden)=vlocal_tmp(:,:,:)
       end do
       ABI_FREE(cgrvtrial)
     end if
     ABI_FREE(vlocal_tmp)
   end if

   !Continue to initialize the Hamiltonian
   call load_spin_hamiltonian(gs_hamk,isppol,vlocal=vlocal,with_nonlocal=.true.)

   do ikpt=1,dtset%nkpt
     nband_k = dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
     istwf_k = dtset%istwfk(ikpt)
     npw_k = npwarr(ikpt)
     kpoint = dtset%kptns(:,ikpt)

     ! Skip the rest of the k-point loop
     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,my_rank)) cycle

     ABI_MALLOC(ylm_k, (npw_k,psps%mpsang**2*psps%useylm))

     kg_k(:,1:npw_k) = kg(:,1+ikg:npw_k+ikg)

     if (psps%useylm==1) then
       do ilm=1,psps%mpsang**2
         ylm_k(1:npw_k,ilm) = ylm(1+ikg:npw_k+ikg,ilm)
       end do
     end if

     ! Set up remaining of the Hamiltonian
     ! Compute (1/2) (2 Pi)**2 (k+G)**2:
     ABI_ALLOCATE(kinpw,(npw_k))
!     call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,kg_k,kinpw,kpoint,npw_k)
     call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,kg_k,kinpw,kpoint,npw_k,0,0)

     ! Compute (k+G) vectors (only if useylm=1)
     nkpg=3*dtset%nloalg(3)
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     if (nkpg>0) then
       call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
     end if

     ! Compute nonlocal form factors ffnl at all (k+G):
     ider=0; idir=0; dimffnl=1
     ABI_ALLOCATE(ffnl, (npw_k,dimffnl,psps%lmnmax,ntypat))

     call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
       gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
       psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
       npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,&
       psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)

     ABI_FREE(kpg_k)
     ABI_FREE(ylm_k)

!    Load k-dependent part in the Hamiltonian datastructure
     ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamk%matblk))
     call load_k_hamiltonian(gs_hamk,kpt_k=dtset%kptns(:,ikpt),npw_k=npw_k,istwf_k=istwf_k,&
&                            kinpw_k=kinpw,kg_k=kg_k,ffnl_k=ffnl,ph3d_k=ph3d,&
&                            compute_ph3d=.true.,compute_gbound=.true.)

     ! Prepare the call to getghc.
     ndat=1; lambda=zero; type_calc=0         ! For applying the whole Hamiltonian
     sij_opt=0; if (psps%usepaw==1) sij_opt=1 ! For PAW, <k+G|1+S|k+G"> is also needed.
     cpopt=-1                                 ! <p_lmn|in> (and derivatives) are computed here (and not saved)

     cplex_ghg=2; size_mat = cplex_ghg*(npw_k*nspinor)**2*dp*b2Mb
     ABI_STAT_MALLOC(ghg_mat,(cplex_ghg,npw_k*nspinor,npw_k*nspinor), ierr)
     write(msg,'(a,f0.3,a)')" Out-of-memory in ghg_mat. Memory required by the Hamiltonian matrix: ",size_mat," [Mb]."
     ABI_CHECK(ierr==0, msg)
     ghg_mat = huge(one)

     ABI_STAT_MALLOC(gtg_mat,(cplex_ghg,npw_k*nspinor,npw_k*nspinor*psps%usepaw), ierr)
     write(msg,'(a,f0.3,a)')" Out-of-memory in gtg_mat. Memory required by the PAW overlap operator: ",size_mat," [Mb]."
     ABI_CHECK(ierr==0, msg)

     ABI_DT_MALLOC(cwaveprj,(natom,nspinor*(1+cpopt)*gs_hamk%usepaw))
     if (cpopt==0) then
       call pawcprj_alloc(cwaveprj,0,gs_hamk%dimcprj)
     end if

     ABI_MALLOC(pwave,(2,npw_k*nspinor))
     pwave=zero ! Initialize plane-wave array:

     ABI_MALLOC(ghc  ,(2,npw_k*nspinor*ndat))
     ABI_MALLOC(gvnlc,(2,npw_k*nspinor*ndat))
     ABI_MALLOC(gsc  ,(2,npw_k*nspinor*ndat*(sij_opt+1)/2))

     if (dtset%prtvol > 0) call wrtout(std_out,' Calculating <G|H|G''> elements','PERS')

     ! Loop over the |beta,G''> component.
     ! Get <:|H|beta,G''> and <:|T_{PAW}|beta,G''>
     do igsp2=1,npw_k*nspinor 
       ghc = zero
       pwave = zero
       pwave(1,igsp2)=one     

       call getghc(cpopt,pwave,cwaveprj,ghc,gsc,gs_hamk,gvnlc,lambda,mpi_enreg,ndat,&
&                  dtset%prtvol,sij_opt,tim_getghc,type_calc)

       ! Fill the upper triangle.
       ghg_mat(:,1:igsp2,igsp2) = ghc(:,1:igsp2)
       if (psps%usepaw==1) gtg_mat(:,1:igsp2,igsp2) = gsc(:,1:igsp2)

       pwave(1,igsp2)=zero ! Reset the |G,beta> component that has been treated.
     end do

     ABI_FREE(kinpw)
     ABI_FREE(ffnl)
     ABI_FREE(ph3d)
     ABI_FREE(pwave)
     ABI_FREE(ghc)
     ABI_FREE(gvnlc)
     ABI_FREE(gsc)

     if (psps%usepaw==1.and.cpopt==0) call pawcprj_free(cwaveprj)
     ABI_DT_FREE(cwaveprj)

     if (dtset%prtvol > 0) then 
       !===========================================
       !=== Diagonalization of <G|H|G''> matrix ===
       !===========================================
       ABI_MALLOC(eig_ene,(nband_k))
       ABI_MALLOC(eig_vec,(cplex_ghg,npw_k*nspinor,nband_k))

       ! * Partial diagonalization
       write(msg,'(2a,3es16.8,3a,i5,a,i0)')ch10,&
       ' Begin partial diagonalization for kpt= ',kpoint,stag(isppol),ch10,&
       ' - Size of mat.=',npw_k*nspinor,' - # nband_k=',nband_k
       if (dtset%prtvol>0) call wrtout(std_out,msg,'PERS')

       jobz = "V"  ! Compute eigenvalues and eigenvectors.
       range = "I" ! the IL-th through IU-th eigenvalues will be found.

       if (psps%usepaw==0) then
         call xheevx(jobz,range,"Upper",cplex_ghg,npw_k*nspinor,ghg_mat,zero,zero,&
         1,nband_k,-tol8,negv,eig_ene,eig_vec,npw_k*nspinor)
       else
         call xhegvx(1,jobz,range,"Upper",cplex_ghg,npw_k*nspinor,ghg_mat,gtg_mat,zero,zero,&
         1,nband_k,-tol8,negv,eig_ene,eig_vec,npw_k*nspinor)
       end if

       ! Write eigenvalues.
       cfact=one !cfact=Ha_eV
       frmt1='(2x,8(1x,es9.2))' ; frmt2='(2x,8(1x,es9.2))'
       write(msg,'(a,3es16.8,3x,a)')' Eigenvalues in Ha for kpt= ',kpoint,stag(isppol)
       call wrtout(std_out,msg,'COLL')

       write(msg,frmt1)(eig_ene(ib)*cfact,ib=1,MIN(8,nband_k))
       call wrtout(std_out,msg,'COLL')
       if (nband_k>8) then
         do jj=9,nband_k,8
           write(msg,frmt2) (eig_ene(ib)*cfact,ib=jj,MIN(jj+7,nband_k))
           call wrtout(std_out,msg,'COLL')
         end do
       end if

#ifdef HAVE_ETSF_IO
       ncerr = nf90_put_var(ncid, eigdiago_varid, eig_ene, start=[1,ikpt,isppol,istep], count=[nband_k,1,1,1])
       NCF_CHECK(ncerr)
#endif

       ABI_FREE(eig_ene)
       ABI_FREE(eig_vec)
     end if

#ifdef HAVE_ETSF_IO
     !write(std_out,*)"Writing H_GG slice for ikpt",ikpt
     if (isppol == 1) then
       ncerr = nf90_put_var(ncid, kg_varid, kg_k, start=[1,1,ikpt], count=[3,npw_k,1])
       NCF_CHECK(ncerr)
     end if

     ncerr = nf90_put_var(ncid, eighist_varid, eigen(:,ikpt,isppol), &
       start=[1,ikpt,isppol,istep], count=[nband_k,1,1,1])
     NCF_CHECK(ncerr)

     npws = nspinor * npw_k
     ncerr = nf90_put_var(ncid, ghg_varid, ghg_mat, start=[1,1,1,ikpt,isppol,istep], &
       count=[2,npws,npws,1,1,1])
     NCF_CHECK(ncerr)

     if (psps%usepaw == 1) then
       ncerr = nf90_put_var(ncid, pawgtg_varid, gtg_mat, start=[1,1,1,ikpt,isppol,istep], &
         count=[2,npws,npws,1,1,1])
       NCF_CHECK(ncerr)
     end if

     !write(666,*)"{istep: ",istep,", ikpt: ", ikpt, ", isppol: ", isppol, ", npws: ",npws,"}"
     !do jj=1,9
     !  write(666, "(18(es11.3))")ghg_mat(:,jj,:9)
     !  write(666, *)
     !end do
     !write(666, *) ghg_mat
#endif

     ABI_FREE(ghg_mat)
     ABI_FREE(gtg_mat)

     ikg = ikg + npw_k
   end do ! ikpt
 end do ! isppol

 ABI_FREE(kg_k)
 ABI_FREE(vlocal)

 call destroy_hamiltonian(gs_hamk)

#ifdef HAVE_ETSF_IO
 NCF_CHECK(nf90_close(ncid))
#endif

 DBG_EXIT("COLL")

end subroutine gshgg_mkncwrite
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
!!      metric,mkffnl,mkkin
!!
!! SOURCE

subroutine kss_calc_vkb(Psps,kpoint,npw_k,kg_k,rprimd,vkbsign,vkb,vkbd)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kss_calc_vkb'
 use interfaces_41_geometry
 use interfaces_56_recipspace
 use interfaces_66_nonlocal
!End of the abilint section

 implicit none

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
 real(dp) :: ucvol,effmass,ecutsm,ecut
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

 effmass=one; ecutsm=zero; ecut=HUGE(one)
! call mkkin(ecut,ecutsm,effmass,gmet,kg_k,modkplusg,kpoint,npw_k)
 call mkkin(ecut,ecutsm,effmass,gmet,kg_k,modkplusg,kpoint,npw_k,0,0)
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

END MODULE m_io_kss
!!***
