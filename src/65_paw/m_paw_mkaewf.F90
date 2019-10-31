!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paw_mkaewf
!! NAME
!!  m_paw_mkaewf
!!
!! FUNCTION
!! Construct complete AE wave functions on the fine FFT grid adding onsite PAW corrections.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (MG)
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

module m_paw_mkaewf

 use defs_basis
 use defs_wvltypes
 use m_abicore
 use m_xmpi
 use m_hide_blas
 use m_splines
 use m_errors
 use m_nctk
 use m_hdr
 use m_dtset
 use m_dtfil
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_datatypes,   only : ebands_t
 use defs_abitypes,    only : MPI_type
 use m_io_tools,       only : flush_unit
 use m_numeric_tools,  only : wrap2_zero_one
 use m_fftcore,        only : sphereboundary
 use m_geometry,       only : xcart2xred
 use m_crystal,        only : crystal_t
 use m_ebands,         only : ebands_ncwrite
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type, pawtab_get_lsize
 use m_pawfgrtab,      only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_free, pawfgrtab_print
 use m_pawcprj,        only : pawcprj_type, pawcprj_alloc, pawcprj_get, pawcprj_free
 use m_paw_pwaves_lmn, only : paw_pwaves_lmn_t, paw_pwaves_lmn_init, paw_pwaves_lmn_free
 use m_paral_atom,     only : get_my_atmtab, free_my_atmtab
 use m_paw_nhat,       only : nhatgrid
 use m_mpinfo,         only : proc_distrb_cycle
 use m_fft,            only : fourwf

 implicit none

 private

 public :: pawmkaewf

CONTAINS  !========================================================================================
!!***

!!****f* m_paw_mkaewf/pawmkaewf
!! NAME
!! pawmkaewf
!!
!! FUNCTION
!! Construct complete AE wave functions on the fine FFT grid adding onsite PAW corrections.
!!
!! INPUTS
!! crystal<crystal_t>=Crystalline structure
!! ebands<ebands_t>=Electronic energies
!! dimcprj(natom)=array of dimensions of array cprj (not ordered)
!! mband=maximum number of bands
!! mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!! mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!! mkmem=number of k points treated by this node.
!! mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!! [comm_atom]= MPI communicator over atoms
!! mpw=maximum dimensioned size of npw.
!! my_natom=number of atoms treated by current processor
!! natom=number of atoms in cell
!! ntypat=number of types of atoms in the cell
!! nkpt=Total number of k-points
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! unks=unit number for G vectors.
!! nband(nkpt*nsppol)=Number of bands for each k-point and spin.
!! istwfk(nkpt)=Storage mode at each k-point.
!! Pawfgrtab(natom) <type(pawfgrtab_type)> : data about the fine grid around each atom
!! Pawrad(ntypat) <type(pawrad_type)> : radial mesh data for each type of atom
!! Pawtab(ntypat) <type(pawtab_type)> : PAW functions around each type of atom
!! Dtfil <type(datafiles_type)>=variables related to files
!! cg(2,mcg)=planewave coefficients of wavefunctions.
!! Cprj(natom,nspinor*mband*mkmem*nsppol)=<p_lmn|Cnk> coefficients for each WF |Cnk>
!!   and each |p_lmn> non-local projector
!! npwarr(nkpt)=Number of plane waves at each k-point
!! ngfftf(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  Note that ngfftf refers to the fine mesh.
!! kg(3,mpw*mkmem)=reduced planewave coordinates
!! Hdr<hdr_type>=the header of wf, den and pot files
!! kpt(3,nkpt)=reduced coordinates of k points.
!!
!! OUTPUT
!!  ierr=Status error
!!  Main output is written on file (ETSF_IO file format).
!!
!! NOTES
!! In PAW calculations, the pseudized wavefunction us represented
!! on a relatively small plane wave basis set and is not normalized
!! as it does not include the on-site PAW contributions which is described
!! in terms of real spherical harmonics and radial functions.
!! For post-processing and proper visualization, it is necessary
!! to use the full electronic wave function, which is what this subroutine constructs.
!! Specifically, it computes the pseudo part by doing an FFT from G- to r-space
!! using the dense mesh defined by pawecutdg. The on-site PAW terms are also
!! computed in real space inside each sphere and added to the pseudo part.
!! Notice that this formula is expressed on the fine grid, and requires
!! interpolating the PAW radial functions onto this grid, as well as calling
!! initylmr in order to get the angular functions on the grid points.
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      flush_unit,fourwf,free_my_atmtab,get_my_atmtab,nhatgrid
!!      paw_pwaves_lmn_free,paw_pwaves_lmn_init,pawcprj_alloc,pawcprj_free
!!      pawfgrtab_free,pawfgrtab_init,pawfgrtab_print,pawtab_get_lsize
!!      sphereboundary,wrap2_zero_one,wrtout,xcart2xred,xmpi_barrier,xmpi_max
!!      xmpi_sum
!!
!! SOURCE

subroutine pawmkaewf(Dtset,crystal,ebands,my_natom,mpw,mband,mcg,mcprj,nkpt,mkmem,nsppol,nband,&
& istwfk,npwarr,kpt,ngfftf,kg,dimcprj,Pawfgrtab,Pawrad,Pawtab,&
& Hdr,Dtfil,cg,Cprj,MPI_enreg,ierr,pseudo_norms,set_k,set_band , &
& mpi_atmtab,comm_atom) ! Optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_natom,mband,mcg,mcprj,mkmem,mpw,nsppol,nkpt
 integer,intent(in),optional :: comm_atom,set_k,set_band
 integer,intent(out) :: ierr
 type(Datafiles_type),intent(in) :: Dtfil
 type(MPI_type),intent(in) :: MPI_enreg
 type(hdr_type),intent(inout) :: Hdr
 type(dataset_type),intent(in) :: Dtset
 type(crystal_t),intent(in) :: crystal
 type(ebands_t),intent(in) :: ebands
!arrays
 integer,intent(in) :: nband(nkpt*nsppol),istwfk(nkpt),npwarr(nkpt),dimcprj(crystal%natom)
 integer,intent(in) :: ngfftf(18),kg(3,mpw*mkmem)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: kpt(3,nkpt)
 real(dp),optional,intent(out) :: pseudo_norms(nsppol,nkpt,mband)
 type(pawfgrtab_type),intent(in) :: Pawfgrtab(my_natom)
 type(pawrad_type),intent(in) :: Pawrad(crystal%ntypat)
 type(pawtab_type),intent(in) :: Pawtab(crystal%ntypat)
 type(pawcprj_type),intent(in) :: Cprj(crystal%natom,mcprj)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourwf0=0,tim_rwwf0=0,master=0
 integer :: bdtot_index,iband,icg,mgfftf,paral_kgb
 integer :: iatom,iatom_tot,ifgd,ifftsph,ifft,itypat,ispinor,ipw,ndat,ii,i1,i2,i3
 integer :: jl,jm,jlmn,natom
 integer :: max_nfgd,nfgd,ln_size,lmn_size,my_comm_atom,option
 integer :: iorder_cprj,comm_cell,me_kpt,ibsp,ibg,isppol,ikpt,nband_k,cplex
 integer :: n1,n2,n3,n4,n5,n6,ikg,npwout,istwf_k,npw_k
 integer :: nfftot,nprocs
 integer :: optcut,optgr0,optgr1,optgr2,optrad,start_band,start_kpt,stop_kpt,stop_band
 logical :: my_atmtab_allocated,paral_atom
 real(dp),parameter :: weight1=one
 real(dp) :: phj,tphj,re_p,im_p,norm,norm_rerr,max_rerr,imur,reur,arg
 character(len=500) :: msg
 character(len=nctk_slen) :: shape_str
!arrays
 integer,allocatable :: l_size_atm(:)
 integer, pointer :: my_atmtab(:)
 integer,allocatable :: gbound(:,:),kg_k(:,:)
 real(dp) :: red(3),shift(3),rfft(3),kpoint(3),cp_fact(2)
 real(dp),allocatable :: r0shift(:,:,:),phk_atm(:,:,:)
 real(dp),allocatable :: buf_tmp(:,:,:),fofgin(:,:),fofgin_down(:,:),fofgout(:,:)
 real(dp),allocatable :: denpot(:,:,:),fofr(:,:,:,:),fofr_down(:,:,:,:),phkr(:,:)
 real(dp),allocatable :: ur_ae(:,:), ur_pw(:,:),ur_ae_onsite(:,:),ur_ps_onsite(:,:)
 real(dp),allocatable :: ur_mask(:),dummy_1d(:),rsph_red(:,:),rsph_cart(:,:)
 type(pawcprj_type),allocatable :: Cprj_k(:,:)
 type(pawfgrtab_type) :: local_pawfgrtab(my_natom)
 type(paw_pwaves_lmn_t),allocatable :: Paw_onsite(:)
#ifdef HAVE_NETCDF
 integer :: fform,ncerr,ncid,ae_ncid,pw_ncid,aeons_ncid,psons_ncid
 character(len=fnlen) :: fname
#endif

! ************************************************************************

 DBG_ENTER("COLL")

!Init parallelism
 comm_cell = MPI_enreg%comm_cell; nprocs = xmpi_comm_size(comm_cell)
 me_kpt = MPI_enreg%me_kpt; paral_kgb=mpi_enreg%paral_kgb

!Compatibility tests
 ABI_CHECK(mkmem/=0, "mkmem==0 not supported anymore!")
 ABI_CHECK(MPI_enreg%paral_kgb == 0, "paral_kgb/=0 not coded")
 ABI_CHECK(SIZE(dimcprj)>0, "dimcprj should be allocated")
 ABI_CHECK(mpi_enreg%paral_spinor==0, "parallelisation over spinors not implemented")
 ABI_CHECK(nprocs==1, "k spin parallelism not yet active")
 ABI_CHECK(dtset%nspinor==1, "nspinor == 2 is buggy")

 natom = crystal%natom

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!If collection of pseudo norms is enabled, make sure the array is initialised
 if (present(pseudo_norms)) pseudo_norms = zero

!use a local copy of pawfgrtab to make sure we use the correction in the paw spheres
!the usual pawfgrtab uses r_shape which may not be the same as r_paw
 if (my_natom>0) then
   if (paral_atom) then
     call pawtab_get_lsize(pawtab,l_size_atm,my_natom,Dtset%typat,mpi_atmtab=my_atmtab)
     call pawfgrtab_init(local_pawfgrtab,Pawfgrtab(1)%cplex,l_size_atm,Dtset%nspden,Dtset%typat,&
&     mpi_atmtab=my_atmtab,comm_atom=my_comm_atom)
   else
     call pawtab_get_lsize(pawtab,l_size_atm,my_natom,Dtset%typat)
     call pawfgrtab_init(local_pawfgrtab,Pawfgrtab(1)%cplex,l_size_atm,Dtset%nspden,Dtset%typat)
   end if
   ABI_FREE(l_size_atm)
 end if
 optcut = 1 ! use rpaw to construct local_pawfgrtab
 optgr0 = 0; optgr1 = 0; optgr2 = 0 ! dont need gY terms locally
 optrad = 1 ! do store r-R

 if (paral_atom) then
   call nhatgrid(crystal%atindx1,crystal%gmet,my_natom,natom,crystal%nattyp,ngfftf,crystal%ntypat,&
&   optcut,optgr0,optgr1,optgr2,optrad,local_pawfgrtab,pawtab,crystal%rprimd,Dtset%typat,crystal%ucvol,Hdr%xred,&
&   comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
 else
   call nhatgrid(crystal%atindx1,crystal%gmet,my_natom,natom,crystal%nattyp,ngfftf,crystal%ntypat,&
&   optcut,optgr0,optgr1,optgr2,optrad,local_pawfgrtab,pawtab,crystal%rprimd,Dtset%typat,crystal%ucvol,Hdr%xred)
 end if
!now local_pawfgrtab is ready to use

 max_nfgd=MAXVAL(local_pawfgrtab(:)%nfgd) ! MAX no. of points in the fine grid for this PAW sphere
 ABI_MALLOC(r0shift,(3,max_nfgd,my_natom))
 ABI_MALLOC(phk_atm,(2,max_nfgd,my_natom))

 do iatom=1,my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)

   nfgd=local_pawfgrtab(iatom)%nfgd ! no. of points in the fine grid for this PAW sphere
   ABI_MALLOC(rsph_red,(3,nfgd))
   ABI_MALLOC(rsph_cart,(3,nfgd))
   do ifgd=1,nfgd
     rsph_cart(:,ifgd) = local_pawfgrtab(iatom)%rfgd(:,ifgd) + crystal%xcart(:,iatom_tot)
   end do
   call xcart2xred(nfgd,crystal%rprimd,rsph_cart,rsph_red) ! we work in reduced coordinates.
   do ifgd=1,nfgd
     call wrap2_zero_one(rsph_red(1,ifgd),red(1),shift(1)) ! num = red + shift
     call wrap2_zero_one(rsph_red(2,ifgd),red(2),shift(2))
     call wrap2_zero_one(rsph_red(3,ifgd),red(3),shift(3))
     r0shift(:,ifgd,iatom) = shift
     !if (ANY( ABS(shift) > tol12)) then
     !  MSG_WARNING("rmR_red is outside the first unit cell.")
     !  write(std_out,*)rsph_red(:,ifgd),shift
     !end if
   end do
   ABI_FREE(rsph_red)
   ABI_FREE(rsph_cart)
 end do

 if (.not.paral_atom .and. my_natom>0) then
   call pawfgrtab_print(local_pawfgrtab,natom=natom,unit=std_out,&
&   prtvol=Dtset%prtvol,mode_paral="COLL")
 end if

 ierr=0
#ifndef HAVE_NETCDF
 ierr = -1
 write(msg,'(3a)')&
& "netcdf support must be enabled in order to output AE PAW wavefunction. ",ch10,&
& "No output will be produced, use --enable-netcdf at configure-time. "
 MSG_WARNING(msg)
 return
!These statements are necessary to avoid the compiler complain about unused variables:
 ii=Dtset%usepaw;ii=Dtfil%unpaw;ii=Hdr%usepaw
#endif

!FIXME check ordering in cprj and Eventually in external file
!why is iorder_cprj not stored in the file for crosschecking purpose?
!Here Im assuming cprj are not ordered!
 iorder_cprj=0

!n4,n5,n6 are FFT dimensions, modified to avoid cache trashing
 n1=ngfftf(1); n2=ngfftf(2); n3=ngfftf(3)
 n4=ngfftf(4); n5=ngfftf(5); n6=ngfftf(6)
 nfftot=PRODUCT(ngfftf(1:3))
 mgfftf=MAXVAL(ngfftf(1:3))

 ABI_MALLOC(phkr,(2,nfftot))
 ABI_MALLOC(gbound,(2*mgfftf+8,2))

#ifdef HAVE_NETCDF
!=== Initialize ETSF_IO files ===
! FIXME: nspinor == 2 is buggy

 fname = trim(dtfil%filnam_ds(4))//'_PAWAVES.nc'
 write(msg,'(2a)')' Opening file for AE PAW wave functions: ',trim(fname)
 call wrtout([std_out, ab_out], msg, 'PERS')

 if (xmpi_comm_rank(comm_cell) == master) then
   NCF_CHECK(nctk_open_create(ncid, fname, xmpi_comm_self))

   fform = 602
   NCF_CHECK(hdr%ncwrite(ncid, fform, nc_define=.True.))

   ! Define wavefunctions in real space on the dense FFT mesh
   ! Fortran layout:
   !real_space_wavefunctions: double 8d array with shape:
   !  [real_or_complex_wavefunctions]
   !  [number_of_grid_points_vector1][number_of_grid_points_vector2][number_of_grid_points_vector3]
   !  [number_of_spinor_components]
   !  [max_number_of_states][number_of_kpoints][number_of_spins]

   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("real_or_complex_wavefunctions", 2),  &
     nctkdim_t("number_of_grid_points_vector1", n1), &
     nctkdim_t("number_of_grid_points_vector2", n2), &
     nctkdim_t("number_of_grid_points_vector3", n3)  &
   ], defmode=.True.)
   NCF_CHECK(ncerr)

   shape_str = "real_or_complex_wavefunctions, &
&   number_of_grid_points_vector1, number_of_grid_points_vector2, number_of_grid_points_vector3, &
&   number_of_spinor_components, &
&   max_number_of_states, number_of_kpoints, number_of_spins"

   ! Define wavefunctions in real space.
   ncerr = nctk_def_arrays(ncid, [&
     nctkarr_t('ur_ae', "dp", shape_str),&
     nctkarr_t('ur_pw', "dp", shape_str),&
     nctkarr_t('ur_ae_onsite', "dp", shape_str),&
     nctkarr_t('ur_ps_onsite', "dp", shape_str) &
   ], defmode=.True.)
   NCF_CHECK(ncerr)

   ! Complete the geometry information.
   NCF_CHECK(crystal%ncwrite(ncid))
   NCF_CHECK(ebands_ncwrite(ebands, ncid))

   NCF_CHECK(nf90_close(ncid))
 end if

 call xmpi_barrier(comm_cell)

 ! Reopen the file in parallel inside comm_cell
 ! Note that we use individual IO thus there's no need to handle idle processes
 ! if paral_kgb == 0 and nprocs > nkpt * nsppol
 NCF_CHECK(nctk_open_modify(ncid, fname, comm_cell))
 ae_ncid = nctk_idname(ncid, "ur_ae")
 pw_ncid = nctk_idname(ncid, "ur_pw")
 aeons_ncid = nctk_idname(ncid, "ur_ae_onsite")
 psons_ncid = nctk_idname(ncid, "ur_ps_onsite")

 NCF_CHECK(nctk_set_datamode(ncid))
#endif

!Init structure storing phi_{nlm} and tphi_(nlm} on the dense FFT points located in the PAW spheres.
 ABI_DT_MALLOC(Paw_onsite,(natom))
 if (paral_atom) then
   call paw_pwaves_lmn_init(Paw_onsite,my_natom,natom,crystal%ntypat,crystal%rprimd,crystal%xcart,&
   Pawtab,Pawrad,local_pawfgrtab, comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
 else
   call paw_pwaves_lmn_init(Paw_onsite,my_natom,natom,crystal%ntypat,crystal%rprimd,crystal%xcart,&
   Pawtab,Pawrad,local_pawfgrtab)
 end if

 bdtot_index=0; icg=0; ibg=0; norm_rerr=smallest_real

 ! === Loop over spin ===
 do isppol=1,nsppol
   ikg=0; start_kpt=1; stop_kpt=nkpt

   ! Check if k-point was specified (only serial)
   if (present(set_k) .and. nprocs==1) then
     if (set_k/=0) then
       start_kpt = set_k
       stop_kpt = set_k
       !MSG_ERROR("set_k")
     end if
   end if

   ! === Loop over k points ===
   do ikpt=start_kpt,stop_kpt
     kpoint  = kpt(:,ikpt)
     nband_k = nband(ikpt+(isppol-1)*nkpt)
     npw_k   = npwarr(ikpt)
     istwf_k = istwfk(ikpt)

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_kpt)) then
       bdtot_index=bdtot_index+nband_k
       !MSG_ERROR("cycle in seq!")
       cycle
     end if

     do i3=0,n3-1
       rfft(3)=DBLE(i3)/n3
       do i2=0,n2-1
         rfft(2)=DBLE(i2)/n2
         do i1=0,n1-1
           rfft(1)=DBLE(i1)/n1
           ifft = 1 +i1 +i2*n1 +i3*n1*n2
           phkr(1,ifft) = COS(two_pi*dot_product(kpoint,rfft))
           phkr(2,ifft) = SIN(two_pi*dot_product(kpoint,rfft))
         end do
       end do
     end do
     ! phkr(1,:)=one; phkr(2,:)=zero

!    Calculate the phase for the onsite PAW contributions.
     do iatom=1,my_natom
       nfgd=local_pawfgrtab(iatom)%nfgd ! no. of points in the fine grid for this PAW sphere
       do ifgd=1,nfgd
         arg = -two_pi* dot_product(r0shift(:,ifgd,iatom),kpoint)
         phk_atm(1,ifgd,iatom) = COS(arg)
         phk_atm(2,ifgd,iatom) = SIN(arg)
       end do
     end do

     ABI_DT_MALLOC(Cprj_k,(natom,dtset%nspinor*nband_k))
     call pawcprj_alloc(Cprj_k,0,dimcprj)

!    Extract cprj for this k-point.
     ibsp=0
     do iband=1,nband_k
       do ispinor=1,dtset%nspinor
         ibsp=ibsp+1
         do iatom=1,natom
           Cprj_k(iatom,ibsp)%cp(:,:)=Cprj(iatom,ibsp+ibg)%cp(:,:)
         end do
       end do
     end do

     ABI_MALLOC(kg_k,(3,npw_k))

     ! Extract G-vectors.
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     call sphereboundary(gbound,istwf_k,kg_k,mgfftf,npw_k)

     ! If a single band is requested, neuter the loop (only serial)
     start_band = 1; stop_band = nband_k
     if (present(set_band).AND.nprocs==1) then
       if (set_band/=0) then
         start_band = set_band
         stop_band = set_band
         !MSG_ERROR("set_band")
       end if
     end if

     ! Loop over bands.
     do iband=start_band,stop_band

       ! Fourier transform on the real fft box of the smooth part.
       ndat=Dtset%nspinor
       ABI_MALLOC(fofgin,(2,npw_k*ndat))
       ABI_MALLOC(fofr,(2,n4,n5,n6*ndat))

       do ipw=1,npw_k*dtset%nspinor
         fofgin(:,ipw)=cg(:,ipw+(iband-1)*npw_k*dtset%nspinor+icg)
       end do

       ! Complex can be set to 0 with this option(0) of fourwf
       option=0; cplex=0; npwout=1
       ABI_MALLOC(denpot,(cplex*n4,n5,n6))
       ABI_MALLOC(fofgout,(2,npwout*ndat))

       call fourwf(cplex,denpot,fofgin(:,1:npw_k),fofgout,fofr(:,:,:,1:n6),gbound,gbound,istwf_k,kg_k,kg_k,&
         mgfftf,MPI_enreg,1,ngfftf,npw_k,npwout,n4,n5,n6,option,tim_fourwf0,weight1,weight1,&
         use_gpu_cuda=Dtset%use_gpu_cuda)

!      Here I do not know if fourwf works in the case of spinors,
!      It seems that not all fftalg option support ndata! should check!
!      Do not forget to declare real(dp)::fofgin_down(:,:) to use the following statements
       if (Dtset%nspinor==2) then
         ABI_MALLOC(fofgin_down,(2,npw_k))
         ABI_MALLOC(fofr_down,(2,n4,n5,n6))
         fofgin_down(:,:)=fofgin(:,1+npw_k:2*npw_k)
!        Complex can be set to 0 with this option(0) of fourwf
!        cplex=1; option=1; npwout=1; ndat=1
!        NOTE: fofr_down can NOT be replaced by fofr(:,:,:,n6+1:2*n6), or else
!        the data in fofr(:,:,:,1:n6) will be the same with fofr(:,:,:,n6+1:2*n6)
         call fourwf(cplex,denpot,fofgin_down,fofgout,fofr_down,gbound,gbound,istwf_k,kg_k,kg_k,&
           mgfftf,MPI_enreg,1,ngfftf,npw_k,npwout,n4,n5,n6,option,tim_fourwf0,weight1,weight1)
         ABI_FREE(fofgin_down)
       end if

       ABI_MALLOC(ur_ae,(2,n1*n2*n3*ndat))
       ABI_MALLOC(ur_ae_onsite,(2,n1*n2*n3))
       ABI_MALLOC(ur_ps_onsite,(2,n1*n2*n3))
       ABI_MALLOC(ur_pw,(2,n1*n2*n3*ndat))
       ABI_MALLOC(ur_mask,(n1*n2*n3))

       ur_ae=zero;ur_ae_onsite=zero;ur_ps_onsite=zero;ur_pw=zero;ur_mask=zero

       ! * Add phase e^{ikr} since it is contained in cprj.
       do i3=1,n3
         do i2=1,n2
           do i1=1,n1
             ii = i1 + n1*(i2-1)+ n1*n2*(i3-1)
             ur_pw(:,ii)=fofr(:,i1,i2,i3) ! Save pw part separately without the phase.
             ur_ae(1,ii)= fofr(1,i1,i2,i3) * phkr(1,ii) - fofr(2,i1,i2,i3) * phkr(2,ii)
             ur_ae(2,ii)= fofr(1,i1,i2,i3) * phkr(2,ii) + fofr(2,i1,i2,i3) * phkr(1,ii)
             if(Dtset%nspinor==2) then
               ur_pw(:,ii+n1*n2*n3)=fofr_down(:,i1,i2,i3) ! Save pw part separately without the phase.
               ur_ae(1,ii+n1*n2*n3)= fofr_down(1,i1,i2,i3) * phkr(1,ii) - fofr_down(2,i1,i2,i3) * phkr(2,ii)
               ur_ae(2,ii+n1*n2*n3)= fofr_down(1,i1,i2,i3) * phkr(2,ii) + fofr_down(2,i1,i2,i3) * phkr(1,ii)
             end if
           end do
         end do
       end do
       ABI_FREE(fofr)

       if(Dtset%nspinor==2) then
         ABI_FREE(fofr_down)
       end if

       ! === Add onsite term on the augmented FFT mesh ===
       do iatom=1,my_natom
         itypat  =local_pawfgrtab(iatom)%itypat
         lmn_size=Pawtab(itypat)%lmn_size
         ln_size =Pawtab(itypat)%basis_size   ! no. of nl elements in PAW basis
         nfgd    =local_pawfgrtab(iatom)%nfgd ! no. of points in the fine grid for this PAW sphere

         ibsp=(iband-1)*dtset%nspinor
         do ispinor=1,dtset%nspinor
           ibsp=ibsp+1
           do jlmn=1,lmn_size
             jl=Pawtab(itypat)%indlmn(1,jlmn)
             jm=Pawtab(itypat)%indlmn(2,jlmn)
             cp_fact(1) = Cprj_k(iatom,ibsp)%cp(1,jlmn) *sqrt(crystal%ucvol) ! Magic factor
             cp_fact(2) = Cprj_k(iatom,ibsp)%cp(2,jlmn) *sqrt(crystal%ucvol)

             do ifgd=1,nfgd ! loop over fine grid points in current PAW sphere.
               ifftsph = local_pawfgrtab(iatom)%ifftsph(ifgd) ! index of the point on the grid
               phj  = Paw_onsite(iatom)% phi(ifgd,jlmn)
               tphj = Paw_onsite(iatom)%tphi(ifgd,jlmn)
               ! old code
               !re_p = cp_fact(1); im_p = cp_fact(2)
               ! apply the phase
               re_p = cp_fact(1) * phk_atm(1,ifgd,iatom) - cp_fact(2) * phk_atm(2,ifgd,iatom)
               im_p = cp_fact(1) * phk_atm(2,ifgd,iatom) + cp_fact(2) * phk_atm(1,ifgd,iatom)

               ur_ae(1,ifftsph+(ispinor-1)*nfftot) = ur_ae(1,ifftsph+(ispinor-1)*nfftot) + re_p * (phj-tphj)
               ur_ae(2,ifftsph+(ispinor-1)*nfftot) = ur_ae(2,ifftsph+(ispinor-1)*nfftot) + im_p * (phj-tphj)
               ur_ae_onsite(1,ifftsph) = ur_ae_onsite(1,ifftsph) + re_p * phj
               ur_ae_onsite(2,ifftsph) = ur_ae_onsite(2,ifftsph) + im_p * phj
               ur_ps_onsite(1,ifftsph) = ur_ps_onsite(1,ifftsph) + re_p * tphj
               ur_ps_onsite(2,ifftsph) = ur_ps_onsite(2,ifftsph) + im_p * tphj
               ur_mask(ifftsph) = one
             end do

           end do !jlmn
         end do !ispinor
       end do !iatom

       if (paral_atom) then
         ABI_MALLOC(buf_tmp,(2,n1*n2*n3,3))
         buf_tmp(:,:,1) = ur_ae
         buf_tmp(:,:,2) = ur_ae_onsite
         buf_tmp(:,:,3) = ur_ps_onsite
         call xmpi_sum(buf_tmp,my_comm_atom,ierr)
         ur_ae = buf_tmp(:,:,1)
         ur_ae_onsite= buf_tmp(:,:,2)
         ur_ps_onsite= buf_tmp(:,:,3)
         ABI_FREE(buf_tmp)
       end if

!      * Remove the phase e^{ikr}, we store u(r).
       do i3=1,n3
         do i2=1,n2
           do i1=1,n1
             ii = i1 + n1*(i2-1)+ n1*n2*(i3-1)
             reur=ur_ae(1,ii)
             imur=ur_ae(2,ii)
             ur_ae(1,ii)=  reur * phkr(1,ii) + imur * phkr(2,ii)
             ur_ae(2,ii)= -reur * phkr(2,ii) + imur * phkr(1,ii)
             if(Dtset%nspinor==2) then
               reur=ur_ae(1,ii+nfftot)    ! Important!
               imur=ur_ae(2,ii+nfftot)
               ur_ae(1,ii+nfftot)=  reur * phkr(1,ii) + imur * phkr(2,ii)
               ur_ae(2,ii+nfftot)= -reur * phkr(2,ii) + imur * phkr(1,ii)
             end if
             reur=ur_ae_onsite(1,ii)
             imur=ur_ae_onsite(2,ii)
             ur_ae_onsite(1,ii)=  reur * phkr(1,ii) + imur * phkr(2,ii)
             ur_ae_onsite(2,ii)= -reur * phkr(2,ii) + imur * phkr(1,ii)
             reur=ur_ps_onsite(1,ii)
             imur=ur_ps_onsite(2,ii)
             ur_ps_onsite(1,ii)=  reur * phkr(1,ii) + imur * phkr(2,ii)
             ur_ps_onsite(2,ii)= -reur * phkr(2,ii) + imur * phkr(1,ii)
           end do
         end do
       end do

       norm=zero
       do ii=1,npw_k*Dtset%nspinor
         norm=norm+fofgin(1,ii)**2+fofgin(2,ii)**2
       end do
       write(std_out,'(a,2i5,f22.16)',advance='no') 'ikpt,iband, norm (G,PSWF)=',ikpt,iband,norm
       norm=zero
       do ifft=1,nfftot*Dtset%nspinor
         norm = norm + ur_ae(1,ifft)**2+ur_ae(2,ifft)**2
       end do
       norm=norm/nfftot
       norm_rerr = MAX((ABS(norm-one))*100,norm_rerr)
       write(std_out,*)"norm (R,AEWF)= ",norm
       call flush_unit(std_out)

!      MS: Various testing and debugging options
       if (.TRUE..and.nprocs==1) then
         if (present(pseudo_norms)) then
!          Check the supposedly zero overlap |\tilde{Psi_n}-\tilde{Psi_n^1}|^2
           ABI_MALLOC(dummy_1d,(n1*n2*n3))
           dummy_1d = zero
           norm = zero
           do ifft = 1, nfftot
             dummy_1d(ifft) = ((ur_pw(1,ifft)-ur_ps_onsite(1,ifft))**2 &
              +  (ur_pw(2,ifft)-ur_ps_onsite(2,ifft))**2) * ur_mask(ifft)
             norm = norm + dummy_1d(ifft)
           end do
           norm = norm / nfftot
           pseudo_norms(isppol,ikpt,iband) = norm
           ABI_FREE(dummy_1d)
         end if

       else
         write(msg,'(5a)')&
          "The option to print PAW all-electron wavefunctions is on, but execution ",ch10,&
          "is in parallel on two or more processors. XcrysDen files with individual con-",ch10,&
          "tributions will not be written. In order to enable this you must run in serial."
         MSG_WARNING(msg)
       end if ! Check if serial run

#ifdef HAVE_NETCDF
       ncerr = nf90_put_var(ncid, ae_ncid, ur_ae, &
          start=[1,1,1,1,1,iband,ikpt,isppol], count=[2,n1,n2,n3,1,1,1,1])
       NCF_CHECK(ncerr)

       ncerr = nf90_put_var(ncid, pw_ncid, ur_pw, &
         start=[1,1,1,1,1,iband,ikpt,isppol], count=[2,n1,n2,n3,1,1,1,1])
       NCF_CHECK(ncerr)

       ncerr = nf90_put_var(ncid, aeons_ncid, ur_ae_onsite, &
         start=[1,1,1,1,1,iband,ikpt,isppol], count=[2,n1,n2,n3,1,1,1,1])
       NCF_CHECK(ncerr)

       ncerr = nf90_put_var(ncid, psons_ncid, ur_ps_onsite, &
         start=[1,1,1,1,1,iband,ikpt,isppol], count=[2,n1,n2,n3,1,1,1,1])
       NCF_CHECK(ncerr)
#endif

       ABI_FREE(ur_ae)
       ABI_FREE(ur_ae_onsite)
       ABI_FREE(ur_ps_onsite)
       ABI_FREE(ur_pw)
       ABI_FREE(ur_mask)
       ABI_FREE(fofgin)
       ABI_FREE(fofgout)
       ABI_FREE(denpot)
     end do !nband_k

     bdtot_index=bdtot_index+nband_k

     if (mkmem/=0) then
       ibg=ibg+dtset%nspinor*nband_k
       icg=icg+npw_k*dtset%nspinor*nband_k
       ikg=ikg+npw_k
     end if

     ABI_FREE(kg_k)

     call pawcprj_free(Cprj_k)
     ABI_DT_FREE(Cprj_k)

   end do !ikpt
 end do !nsppol

 ABI_FREE(phkr)
 ABI_FREE(gbound)

 ! Free augmentation waves.
 call paw_pwaves_lmn_free(Paw_onsite)
 ABI_DT_FREE(Paw_onsite)

 ! Maximum relative error over CPUs.
 call xmpi_max(norm_rerr,max_rerr,comm_cell,ierr)
 write(std_out,*)"max_rerr=",max_rerr

 if (max_rerr > ten) then
   write(msg,'(7a)')&
    "Inaccuracy on the normalization of the wave funtions exceeds 10%. ",ch10,&
    "Likely due to the use of a too coarse FFT mesh or unconverged wavefunctions. ",ch10,&
    "Numerical values inside the augmentation regions might be inaccurate. ",ch10,&
    "Action: increase pawecutdg in your input file. "
   MSG_COMMENT(msg)
 end if

 ABI_FREE(r0shift)
 ABI_FREE(phk_atm)
 call pawfgrtab_free(local_pawfgrtab)

 ! Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

end subroutine pawmkaewf
!!***

end module m_paw_mkaewf
!!***
