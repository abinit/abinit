!!****m* ABINIT/m_paw_optics
!! NAME
!!  m_paw_optics
!!
!! FUNCTION
!!  This module contains several routines related to conductivity:
!!    optical conductivity, X spectroscopy, linear susceptibility, ...
!!
!! COPYRIGHT
!! Copyright (C) 2018-2022 ABINIT group (SM,VR,FJ,MT,NB,PGhosh)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_optics

 use defs_basis
 use m_xmpi
 use m_errors
 use m_wffile
 use m_abicore
 use m_hdr
 use m_dtset
 use m_dtfil
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes,  only : MPI_type
 use m_time,         only : timab
 use m_io_tools,     only : open_file,get_unit,close_unit
 use m_pawpsp,       only : pawpsp_read_corewf
 use m_pawrad,       only : pawrad_type,pawrad_deducer0,simp_gen,nderiv_gen,poisson
 use m_pawtab,       only : pawtab_type
 use m_pawcprj,      only : pawcprj_type,pawcprj_alloc,pawcprj_get, &
&                           pawcprj_free,pawcprj_mpi_allgather
 use m_pawang,       only : pawang_type
 use m_paw_denpot,   only : pawdensities,pawkindensities,pawdenpot
 use m_paw_an,       only : paw_an_type,paw_an_init,paw_an_free,paw_an_copy
 use m_pawrhoij,     only : pawrhoij_type
 use m_paw_ij,       only : paw_ij_type
 use m_paw_onsite,   only : pawnabla_init,pawnabla_core_init
 use m_paw_sphharm,  only : setnabla_ylm
 use m_pawxc,        only : pawxc,pawxcm,pawxc_get_xclevel,pawxc_get_usekden
 use m_mpinfo,       only : destroy_mpi_enreg,nullify_mpi_enreg,initmpi_seq,proc_distrb_cycle
 use m_numeric_tools,only : kramerskronig
 use m_geometry,     only : metric
 use m_hide_lapack,  only : matrginv
 use m_paral_atom,       only : get_my_atmtab,free_my_atmtab

 implicit none

 private

!public procedures.
 public :: optics_paw
 public :: optics_paw_core
 public :: linear_optics_paw

!I/O parameters
!Set to true to use netcdf-MPIIO when available
 logical,parameter :: use_netcdf_mpiio=.true.
!Set to true to use unlimited dimensions in netCDF file
!  This is not mandatory because we know exactly the amount of data to write
!    and seems to impact performances negatively
 logical,parameter :: use_netcdf_unlimited=.false.

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_optics/optics_paw
!! NAME
!! optics_paw
!!
!! FUNCTION
!! Compute matrix elements need for optical conductivity (in the PAW context) and store them in a file
!!  Matrix elements = <Phi_i|Nabla|Phi_j>
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  cg(2,mcg)=planewave coefficients of wavefunctions.
!!  cprj(natom,mcprj)= <p_lmn|Cnk> coefficients for each WF |Cnk>
!!                                          and each |p_lmn> non-local projector
!!  dimcprj(natom)=array of dimensions of array cprj (not ordered)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mkmem =number of k points treated by this node.
!!  mpi_enreg=information about MPI parallelization
!!  mpsang =1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw.
!!  natom=number of atoms in cell.
!!  nkpt=number of k points.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  pawang <type(pawang_type)>= PAW ANGular mesh discretization and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= PAW rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  znucl(ntypat)=atomic number of atom type
!!
!! OUTPUT
!!  (only writing in a file)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_outscfcv
!!
!! CHILDREN
!!      destroy_mpi_enreg,hdr%free,hdr_fort_read,hdr%ncread
!!      kramerskronig,matrginv,metric
!!
!! SOURCE

 subroutine optics_paw(atindx1,cg,cprj,dimcprj,dtfil,dtset,eigen0,gprimd,hdr,kg,&
&               mband,mcg,mcprj,mkmem,mpi_enreg,mpsang,mpw,natom,nkpt,npwarr,nsppol,&
&               pawang,pawrad,pawrhoij,pawtab,znucl)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcg,mcprj,mkmem,mpsang,mpw,natom,nkpt,nsppol
 type(MPI_type),intent(in) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: atindx1(natom),dimcprj(natom),npwarr(nkpt)
 integer,intent(in),target :: kg(3,mpw*mkmem)
 real(dp),intent(in) :: eigen0(mband*nkpt*nsppol)
 real(dp),intent(in) :: gprimd(3,3),znucl(dtset%ntypat)
 real(dp),intent(inout) :: cg(2,mcg)
 type(pawcprj_type),target,intent(inout) :: cprj(natom,mcprj)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
 type(pawrhoij_type),intent(inout) :: pawrhoij(mpi_enreg%my_natom)
 type(pawtab_type),target,intent(inout) :: pawtab(dtset%ntypat)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: iomode,bdtot_index,cplex,etiq,fformopt,iatom,ib,ibmax,ibg,ibsp
 integer :: icg,ierr,ikg,ikpt,ilmn,ount,ncid,varid,idir
 integer :: iorder_cprj,ipw,ispinor,isppol,istwf_k,itypat,iwavef
 integer :: jb,jbsp,my_jb,jlmn,jwavef,lmn_size,mband_cprj,option_core
 integer :: my_nspinor,nband_k,nband_cprj_k,npw_k,sender,me,master_spfftband
 integer :: spaceComm_band,spaceComm_bandspinorfft,spaceComm_fft,spaceComm_kpt
 integer :: spaceComm_spinor,spaceComm_bandspinor,spaceComm_spinorfft,spaceComm_w
 logical :: already_has_nabla,cprj_paral_band,myband,mykpt,iomode_etsf_mpiio
 logical :: i_am_master,i_am_master_kpt,i_am_master_band,i_am_master_spfft,nc_unlimited
 real(dp) :: cgnm1,cgnm2,cpnm1,cpnm2,cpnm11,cpnm22,cpnm12,cpnm21,cpnm_11m22,cpnm_21p12,cpnm_21m12
 character(len=500) :: msg
 type(nctkdim_t) :: nctkdim
!arrays
 integer :: nc_count(6),nc_start(6),nc_stride(6),tmp_shape(3)
 integer, ABI_CONTIGUOUS pointer :: kg_k(:,:)
 real(dp) :: kpoint(3),tsec(2),nabla_ij(3)
 real(dp),allocatable :: kpg_k(:,:),psinablapsi(:,:,:,:)
 real(dp),allocatable :: psinablapsi_paw(:,:,:,:),psinablapsi_soc(:,:,:,:)
 real(dp),pointer :: soc_ij(:,:,:)
 type(coeff5_type),allocatable,target :: phisocphj(:)
 type(pawcprj_type),pointer :: cprj_k(:,:),cprj_k_loc(:,:)
 type(nctkarr_t) :: nctk_arrays(1)

! ************************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 ABI_CHECK(mkmem/=0,"mkmem==0 not supported anymore!")
 ABI_CHECK(mpi_enreg%paral_spinor==0.or.dtset%pawspnorb==0,"spinor parallelization not supported with SOC!")
!  MJV 6/12/2008: looks like mpi_enreg may not be completely initialized here
 if (xmpi_paral==1) then
   tmp_shape = shape(mpi_enreg%proc_distrb)
   if (nkpt > tmp_shape(1)) then
     ABI_BUG('problem with proc_distrb!')
   end if
 end if

!Init parallelism
 spaceComm_w=mpi_enreg%comm_cell
 if (mpi_enreg%paral_kgb==1) then
   spaceComm_kpt=mpi_enreg%comm_kpt
   spaceComm_fft=mpi_enreg%comm_fft
   spaceComm_band=mpi_enreg%comm_band
   spaceComm_spinor=mpi_enreg%comm_spinor
   spaceComm_bandspinor=mpi_enreg%comm_bandspinor
   spaceComm_spinorfft=mpi_enreg%comm_spinorfft
   spaceComm_bandspinorfft=mpi_enreg%comm_bandspinorfft
 else
   spaceComm_kpt=mpi_enreg%comm_kpt
   spaceComm_fft=xmpi_comm_self
   spaceComm_band=mpi_enreg%comm_band
   spaceComm_spinor=xmpi_comm_self
   spaceComm_bandspinor=spaceComm_band
   spaceComm_spinorfft=xmpi_comm_self
   spaceComm_bandspinorfft=xmpi_comm_self
 end if
 me=xmpi_comm_rank(spaceComm_w)
 i_am_master=(me==master)
 i_am_master_kpt=(xmpi_comm_rank(spaceComm_kpt)==master)
 i_am_master_band=(xmpi_comm_rank(spaceComm_band)==master)
 i_am_master_spfft=(xmpi_comm_rank(spaceComm_spinorfft)==master)
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

!----------------------------------------------------------------------------------
!1- Opening of OPT file and header writing
!----------------------------------------------------------------------------------

!(master proc only)
 if (i_am_master) then
   fformopt=610
!  ====> NETCDF format
   if (dtset%iomode==IO_MODE_ETSF) then
     iomode=IO_MODE_ETSF
#ifdef HAVE_NETCDF
!    Open/create nc file
     NCF_CHECK(nctk_open_create(ncid,nctk_ncify(dtfil%fnameabo_app_opt),xmpi_comm_self))
!    Write header data
     NCF_CHECK(hdr%ncwrite(ncid,fformopt,nc_define=.true.))
!    Define dims and array for dipole variables
     nctk_arrays(1)%name="dipole_valence_valence"
     nctk_arrays(1)%dtype="dp"
     nc_unlimited=(use_netcdf_unlimited.and.(.not.(nctk_has_mpiio.and.use_netcdf_mpiio)))
     if (nc_unlimited) then
       nctkdim%name="unlimited_bands"
       nctkdim%value=NF90_UNLIMITED
       NCF_CHECK(nctk_def_dims(ncid,nctkdim))
       nctk_arrays(1)%shape_str=&
&     "complex,number_of_cartesian_directions,max_number_of_states,number_of_kpoints,number_of_spins,unlimited_bands"
     else
       nctk_arrays(1)%shape_str=&
&     "complex,number_of_cartesian_directions,max_number_of_states,max_number_of_states,number_of_kpoints,number_of_spins"
     end if
     NCF_CHECK(nctk_def_arrays(ncid, nctk_arrays))
     NCF_CHECK(nctk_set_atomic_units(ncid, "dipole_valence_valence"))
!    Write eigenvalues
     NCF_CHECK(nctk_set_datamode(ncid))
     varid=nctk_idname(ncid,"eigenvalues")
     NCF_CHECK(nf90_put_var(ncid,varid,reshape(eigen0,[mband,nkpt,nsppol])))
     !Close file here because the rest has possibly to be written with collective I/O
     NCF_CHECK(nf90_close(ncid))
#else
     msg = "In order to use iomode=3 and prtnabla>0 together, NetCDF support must be enabled!"
     ABI_ERROR(msg)
#endif
!  ====> Standard FORTRAN binary format
   else
     iomode=IO_MODE_FORTRAN_MASTER
     if (open_file(dtfil%fnameabo_app_opt,msg,newunit=ount,form="unformatted",status="unknown")/= 0) then
       ABI_ERROR(msg)
     end if
     call hdr%fort_write(ount,fformopt,ierr,rewind=.true.)
     write(ount)(eigen0(ib),ib=1,mband*nkpt*nsppol)
   end if ! File format
 end if ! master node
 call xmpi_bcast(iomode,master,spaceComm_w,ierr)
 iomode_etsf_mpiio=(iomode==IO_MODE_ETSF.and.nctk_has_mpiio.and.use_netcdf_mpiio)
 nc_unlimited=(iomode==IO_MODE_ETSF.and.use_netcdf_unlimited.and.(.not.iomode_etsf_mpiio)) ! UNLIMITED not compatible with mpi-io

!----------------------------------------------------------------------------------
!2- Computation of on-site contribution: <phi_i|nabla|phi_j>-<tphi_i|nabla|tphi_j>
!----------------------------------------------------------------------------------

 already_has_nabla=all(pawtab(:)%has_nabla==2)
 call pawnabla_init(mpsang,dtset%ntypat,pawrad,pawtab)
 
!Compute spin-orbit contributions if necessary
 if (dtset%pawspnorb==1) then
   option_core=0
   call pawnabla_soc_init(phisocphj,option_core,dtset%ixc,mpi_enreg%my_natom,natom,&
&       dtset%nspden,dtset%ntypat,pawang,pawrad,pawrhoij,pawtab,dtset%pawxcdev,&
&       dtset%spnorbscl,dtset%typat,dtset%xc_denpos,znucl,&
&       comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
 end if

!----------------------------------------------------------------------------------
!3- Computation of <psi_n|-i.nabla|psi_m> for each k
!----------------------------------------------------------------------------------

!Prepare valence-valence dipoles writing
!In case of netCDF access to OPT file, prepare collective I/O
 if (iomode == IO_MODE_ETSF) then
   if (iomode_etsf_mpiio) then
     if (i_am_master_spfft) then
       NCF_CHECK(nctk_open_modify(ncid,nctk_ncify(dtfil%fnameabo_app_opt),spaceComm_band))
       varid=nctk_idname(ncid,"dipole_valence_valence")
       if (xmpi_comm_size(spaceComm_w)>1) then
         NCF_CHECK(nctk_set_collective(ncid,varid))
       end if
       NCF_CHECK(nctk_set_datamode(ncid))
     end if
   else if (i_am_master) then
     NCF_CHECK(nctk_open_modify(ncid,nctk_ncify(dtfil%fnameabo_app_opt),xmpi_comm_self))
     varid=nctk_idname(ncid,"dipole_valence_valence")
     if (nctk_has_mpiio.and.(.not.use_netcdf_mpiio)) then
       NCF_CHECK(nctk_set_collective(ncid,varid))
     end if
     NCF_CHECK(nctk_set_datamode(ncid))
   end if     
 end if
 if (iomode_etsf_mpiio) then
   !If MPI-IO, store only ib elements for each jb
   ABI_MALLOC(psinablapsi,(2,3,mband,1))
   ABI_MALLOC(psinablapsi_paw,(2,3,mband,1))
   if (dtset%pawspnorb==1) then
     ABI_MALLOC(psinablapsi_soc,(2,3,mband,1))
   end if
 else
   !If not, store all (ib,jb) pairs
   ABI_MALLOC(psinablapsi,(2,3,mband,mband))
   ABI_MALLOC(psinablapsi_paw,(2,3,mband,mband))
   if (dtset%pawspnorb==1) then
     ABI_MALLOC(psinablapsi_soc,(2,3,mband,mband))
   end if
 end if

!Determine if cprj datastructure is distributed over bands
 mband_cprj=mcprj/(my_nspinor*mkmem*nsppol)
 cprj_paral_band=(mband_cprj<mband)

!LOOP OVER SPINS
 ibg=0;icg=0
 bdtot_index=0
 do isppol=1,nsppol

!  LOOP OVER k POINTS
   ikg=0
   do ikpt=1,nkpt

     etiq=ikpt+(isppol-1)*nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
     master_spfftband=minval(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol))

!    Select k-points for current proc
     mykpt=.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,mpi_enreg%me_kpt))
     if (mykpt) then

!      Data depending on k-point
       npw_k=npwarr(ikpt)
       istwf_k=dtset%istwfk(ikpt)
       cplex=2;if (istwf_k>1) cplex=1
       kpoint(:)=dtset%kptns(:,ikpt)

!      Extract cprj for this k-point
       nband_cprj_k=nband_k;if (cprj_paral_band) nband_cprj_k=nband_k/mpi_enreg%nproc_band
       if (mkmem*nsppol/=1) then
         iorder_cprj=0
         ABI_MALLOC(cprj_k_loc,(natom,my_nspinor*nband_cprj_k))
         call pawcprj_alloc(cprj_k_loc,0,dimcprj)
         call pawcprj_get(atindx1,cprj_k_loc,cprj,natom,1,ibg,ikpt,iorder_cprj,isppol,&
&         mband_cprj,mkmem,natom,nband_cprj_k,nband_cprj_k,my_nspinor,nsppol,dtfil%unpaw,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
       else
         cprj_k_loc => cprj
       end if

!      if cprj are distributed over bands, gather them (because we need to mix bands)
       if (cprj_paral_band) then
         ABI_MALLOC(cprj_k,(natom,my_nspinor*nband_k))
         call pawcprj_alloc(cprj_k,0,dimcprj)
         call pawcprj_mpi_allgather(cprj_k_loc,cprj_k,natom, &
&             my_nspinor*nband_cprj_k,my_nspinor*mpi_enreg%bandpp,&
&         dimcprj,0,mpi_enreg%nproc_band,mpi_enreg%comm_band,ierr,rank_ordered=.false.)
       else
         cprj_k => cprj_k_loc
       end if

!      Compute k+G in cartesian coordinates
       ABI_MALLOC(kpg_k,(3,npw_k*dtset%nspinor))
       kg_k => kg(:,1+ikg:npw_k+ikg)
       do ipw=1,npw_k
         kpg_k(1,ipw)=(kpoint(1)+kg_k(1,ipw))*gprimd(1,1) &
&                    +(kpoint(2)+kg_k(2,ipw))*gprimd(1,2) &
&                    +(kpoint(3)+kg_k(3,ipw))*gprimd(1,3)
         kpg_k(2,ipw)=(kpoint(1)+kg_k(1,ipw))*gprimd(2,1) &
&                    +(kpoint(2)+kg_k(2,ipw))*gprimd(2,2) &
&                    +(kpoint(3)+kg_k(3,ipw))*gprimd(2,3)
         kpg_k(3,ipw)=(kpoint(1)+kg_k(1,ipw))*gprimd(3,1) &
&                    +(kpoint(2)+kg_k(2,ipw))*gprimd(3,2) &
&                    +(kpoint(3)+kg_k(3,ipw))*gprimd(3,3)
       end do
       kpg_k=two_pi*kpg_k
       if (dtset%nspinor==2) kpg_k(1:3,npw_k+1:2*npw_k)=kpg_k(1:3,1:npw_k)

!      Loops over bands
       do jb=1,nband_k
         jwavef=(jb-1)*npw_k*my_nspinor+icg
         !If MPI-IO, compute all (ib,jb) pairs ; if not, compute only ib<=jb
         ibmax=merge(nband_k,jb,iomode_etsf_mpiio)
         !If MPI-IO, store only ib elements for each jb ; if not, store all (ib,jb) pairs
         my_jb=merge(1,jb,iomode_etsf_mpiio)

         psinablapsi(:,:,:,my_jb)=zero
         psinablapsi_paw(:,:,:,my_jb)=zero
         if (dtset%pawspnorb==1) psinablapsi_soc(:,:,:,my_jb)=zero

!        2-A Computation of <psi_tild_n|-i.nabla|psi_tild_m>
!        ----------------------------------------------------------------------------------
!        Note: <psi_nk|-i.nabla|psi_mk> => Sum_g[ <G|Psi_nk>^* <G|Psi_mk> G ]

!        Select bands for current proc
         myband=.true.
         if (xmpi_paral==1.and.mpi_enreg%paral_kgb/=1.and.(.not.iomode_etsf_mpiio)) then
           myband=(abs(mpi_enreg%proc_distrb(ikpt,jb,isppol)-mpi_enreg%me_kpt)==0)
         end if
         if (myband) then

           do ib=1,ibmax
             iwavef=(ib-1)*npw_k*my_nspinor+icg

!            (C_nk^*)*C_mk*(k+g) is expressed in cartesian coordinates
             if (istwf_k>1) then
               !In this case (istwfk>1): Sum_g>=g0[ 2i.Im(<G|Psi_nk>^* <G|Psi_mk> G) ]
               !G=k+g=0 term is included but do not contribute
               do ipw=1,npw_k*my_nspinor
                 cgnm2=two*(cg(1,ipw+iwavef)*cg(2,ipw+jwavef)-cg(2,ipw+iwavef)*cg(1,ipw+jwavef))
                 psinablapsi(2,1:3,ib,my_jb)=psinablapsi(2,1:3,ib,my_jb)+cgnm2*kpg_k(1:3,ipw)
               end do
             else
               do ipw=1,npw_k*my_nspinor
                 cgnm1=cg(1,ipw+iwavef)*cg(1,ipw+jwavef)+cg(2,ipw+iwavef)*cg(2,ipw+jwavef)
                 cgnm2=cg(1,ipw+iwavef)*cg(2,ipw+jwavef)-cg(2,ipw+iwavef)*cg(1,ipw+jwavef)
                 psinablapsi(1,1:3,ib,my_jb)=psinablapsi(1,1:3,ib,my_jb)+cgnm1*kpg_k(1:3,ipw)
                 psinablapsi(2,1:3,ib,my_jb)=psinablapsi(2,1:3,ib,my_jb)+cgnm2*kpg_k(1:3,ipw)
               end do
             end if

!            Second half of the (n,m) matrix
             if ((ib/=jb).and.(.not.iomode_etsf_mpiio)) then
               psinablapsi(1,1:3,jb,ib)= psinablapsi(1,1:3,ib,jb)
               psinablapsi(2,1:3,jb,ib)=-psinablapsi(2,1:3,ib,jb)
             end if

           end do ! ib

!          Reduction in case of parallelism
           if (iomode_etsf_mpiio) then
             call timab(48,1,tsec)
             if (mpi_enreg%paral_kgb==1) then
               call xmpi_sum_master(psinablapsi,master,spaceComm_spinorfft,ierr)
             end if
             if (i_am_master_spfft) then
               call xmpi_sum(psinablapsi,spaceComm_band,ierr)
             end if
             call timab(48,2,tsec)
           end if

         end if ! myband

       
!        2-B Computation of <psi_n|p_i><p_j|psi_m>(<phi_i|-i.nabla|phi_j>-<tphi_i|-i.nabla|tphi_j>)
!        ----------------------------------------------------------------------------------
!        Non relativistic contribution
!        Note: <psi|-i.nabla|psi_mk>
!              => -i Sum_ij[ <p_i|Psi_nk>^* <p_j|Psi_mk> (<Phi_i|Nabla|Phi_j>-<tPhi_i|Nabla|tPhi_j>) ]

!        Select bands for current proc
         myband=.true.
         if (mpi_enreg%paral_kgb==1) then
           myband=(mod(jb-1,mpi_enreg%nproc_band)==mpi_enreg%me_band)
         else if (xmpi_paral==1) then
           myband=(abs(mpi_enreg%proc_distrb(ikpt,jb,isppol)-mpi_enreg%me_kpt)==0)
         end if
         if (myband) then

           do ib=1,ibmax          

             ibsp=(ib-1)*my_nspinor ; jbsp=(jb-1)*my_nspinor
             do ispinor=1,my_nspinor
               ibsp=ibsp+1;jbsp=jbsp+1
               if (cplex==1) then
                 do iatom=1,natom
                   itypat=dtset%typat(iatom)
                   lmn_size=pawtab(itypat)%lmn_size
                   do jlmn=1,lmn_size
                     do ilmn=1,lmn_size                  
                       nabla_ij(:)=pawtab(itypat)%nabla_ij(:,ilmn,jlmn)
                       if (ib>jb) nabla_ij(:)=-pawtab(itypat)%nabla_ij(:,jlmn,ilmn) ! We should take 1/2(nabla_ij-nabla_ji)
                       cpnm1=cprj_k(iatom,ibsp)%cp(1,ilmn)*cprj_k(iatom,jbsp)%cp(1,jlmn)
                       if (dtset%nspinor==2) cpnm1=cpnm1+cprj_k(iatom,ibsp)%cp(2,ilmn)*cprj_k(iatom,jbsp)%cp(2,jlmn)
                       psinablapsi_paw(2,:,ib,my_jb)=psinablapsi_paw(2,:,ib,my_jb)-cpnm1*nabla_ij(:)
                     end do !ilmn
                   end do !jlmn
                 end do !iatom
               else
                 do iatom=1,natom
                   itypat=dtset%typat(iatom)
                   lmn_size=pawtab(itypat)%lmn_size
                   do jlmn=1,lmn_size
                     do ilmn=1,lmn_size
                       nabla_ij(:)=pawtab(itypat)%nabla_ij(:,ilmn,jlmn)
                       if (ib>jb) nabla_ij(:)=-pawtab(itypat)%nabla_ij(:,jlmn,ilmn) ! We should take 1/2(nabla_ij-nabla_ji)
                       cpnm1=(cprj_k(iatom,ibsp)%cp(1,ilmn)*cprj_k(iatom,jbsp)%cp(1,jlmn) &
&                            +cprj_k(iatom,ibsp)%cp(2,ilmn)*cprj_k(iatom,jbsp)%cp(2,jlmn))
                       cpnm2=(cprj_k(iatom,ibsp)%cp(1,ilmn)*cprj_k(iatom,jbsp)%cp(2,jlmn) &
&                            -cprj_k(iatom,ibsp)%cp(2,ilmn)*cprj_k(iatom,jbsp)%cp(1,jlmn))
                       psinablapsi_paw(1,:,ib,my_jb)=psinablapsi_paw(1,:,ib,my_jb)+cpnm2*nabla_ij(:)
                       psinablapsi_paw(2,:,ib,my_jb)=psinablapsi_paw(2,:,ib,my_jb)-cpnm1*nabla_ij(:)
                     end do !ilmn
                   end do !jlmn
                 end do !iatom
               end if
             end do !ispinor

!        2-C Computation of  Spin-orbit coupling contribution:
!             Sum_ij,ss'[<psi_n,s|p_i><p_j|psi_m,s'>(<phi_i|1/4 Alpha^2 (Sigma^ss' X dV/dr)|phi_j>]
!        ----------------------------------------------------------------------------------
             if (dtset%pawspnorb==1) then
!              Add: Sum_ij,ss'[<Psi^s_n|p_i><p_j|Psi^s'_m> (Sigma^ss' X g_ij)]
!                  =Sum_ij[ (<Psi^up_n|p_i><p_j|Psi^up_m>-<Psi^dn_n|p_i><p_j|Psi^dn_m>) (Sigma^up-up     X g_ij)
!                          +(<Psi^dn_n|p_i><p_j|Psi^up_m>+<Psi^up_n|p_i><p_j|Psi^dn_m>) (Re{Sigma^dn-up} X g_ij)
!                          +(<Psi^dn_n|p_i><p_j|Psi^up_m>-<Psi^up_n|p_i><p_j|Psi^dn_m>) (Im{Sigma^dn-up} X g_ij) ]
!               where: g_ij = <Phi_i| 1/4 Alpha^2 dV/dr vec(r)/r) |Phi_j>
!              Note that:
!                phisocphj(:)%value(re:im,1,idir,ilmn,jlmn) is (Sigma^up-up X g_ij)
!                phisocphj(:)%value(re:im,2,idir,ilmn,jlmn) is (Sigma^dn-up X g_ij)
!              Not compatible with parallelization over spinors
               ibsp=1+(ib-1)*dtset%nspinor ; jbsp=1+(jb-1)*dtset%nspinor
               do iatom=1,natom
                 itypat=dtset%typat(iatom)
                 lmn_size=pawtab(itypat)%lmn_size
                 do jlmn=1,lmn_size
                   do ilmn=1,lmn_size
                     soc_ij => phisocphj(iatom)%value(:,:,:,ilmn,jlmn)
                     !Contribution from real part of <Psi^s_n|p_i><p_j|Psi^s'_m>
                     cpnm11=cprj_k(iatom,ibsp  )%cp(1,ilmn)*cprj_k(iatom,jbsp  )%cp(1,jlmn) &
&                          +cprj_k(iatom,ibsp  )%cp(2,ilmn)*cprj_k(iatom,jbsp  )%cp(2,jlmn)
                     cpnm22=cprj_k(iatom,ibsp+1)%cp(1,ilmn)*cprj_k(iatom,jbsp+1)%cp(1,jlmn) &
&                          +cprj_k(iatom,ibsp+1)%cp(2,ilmn)*cprj_k(iatom,jbsp+1)%cp(2,jlmn)
                     cpnm12=cprj_k(iatom,ibsp  )%cp(1,ilmn)*cprj_k(iatom,jbsp+1)%cp(1,jlmn) &
&                          +cprj_k(iatom,ibsp  )%cp(2,ilmn)*cprj_k(iatom,jbsp+1)%cp(2,jlmn)
                     cpnm21=cprj_k(iatom,ibsp+1)%cp(1,ilmn)*cprj_k(iatom,jbsp  )%cp(1,jlmn) &
&                          +cprj_k(iatom,ibsp+1)%cp(2,ilmn)*cprj_k(iatom,jbsp  )%cp(2,jlmn)
                     cpnm_11m22=cpnm11-cpnm22
                     cpnm_21p12=cpnm21+cpnm12
                     cpnm_21m12=cpnm21-cpnm12
                     do idir=1,3
                       psinablapsi_soc(1,idir,ib,my_jb)=psinablapsi_soc(1,idir,ib,my_jb) &
&                             +soc_ij(1,1,idir)*cpnm_11m22+soc_ij(1,2,idir)*cpnm_21p12
                       psinablapsi_soc(2,idir,ib,my_jb)=psinablapsi_soc(2,idir,ib,my_jb) &
&                             +soc_ij(2,2,idir)*cpnm_21m12
                     end do
                     !Contribution from imaginary part of <Psi^s_n|p_i><p_j|Psi^s'_m>
                     cpnm11=cprj_k(iatom,ibsp  )%cp(1,ilmn)*cprj_k(iatom,jbsp  )%cp(2,jlmn) &
&                          -cprj_k(iatom,ibsp  )%cp(2,ilmn)*cprj_k(iatom,jbsp  )%cp(1,jlmn)
                     cpnm22=cprj_k(iatom,ibsp+1)%cp(1,ilmn)*cprj_k(iatom,jbsp+1)%cp(2,jlmn) &
&                          -cprj_k(iatom,ibsp+1)%cp(2,ilmn)*cprj_k(iatom,jbsp+1)%cp(1,jlmn)
                     cpnm12=cprj_k(iatom,ibsp  )%cp(1,ilmn)*cprj_k(iatom,jbsp+1)%cp(2,jlmn) &
&                          -cprj_k(iatom,ibsp  )%cp(2,ilmn)*cprj_k(iatom,jbsp+1)%cp(1,jlmn)
                     cpnm21=cprj_k(iatom,ibsp+1)%cp(1,ilmn)*cprj_k(iatom,jbsp  )%cp(2,jlmn) &
&                          -cprj_k(iatom,ibsp+1)%cp(2,ilmn)*cprj_k(iatom,jbsp  )%cp(1,jlmn)
                     cpnm_11m22=cpnm11-cpnm22
                     cpnm_21p12=cpnm21+cpnm12
                     cpnm_21m12=cpnm21-cpnm12
                     do idir=1,3
                       psinablapsi_soc(1,idir,ib,my_jb)=psinablapsi_soc(1,idir,ib,my_jb) &
&                          -soc_ij(2,2,idir)*cpnm_21m12
                       psinablapsi_soc(2,idir,ib,my_jb)=psinablapsi_soc(2,idir,ib,my_jb) &
&                          +soc_ij(1,1,idir)*cpnm_11m22+soc_ij(1,2,idir)*cpnm_21p12
                     end do
                   end do ! ilmn
                 end do ! jlmn
               end do ! iatom
             end if ! pawspnorb

!            Second half of the (n,m) matrix
             if ((ib/=jb).and.(.not.iomode_etsf_mpiio)) then
               psinablapsi_paw(1,1:3,jb,ib)= psinablapsi_paw(1,1:3,ib,jb)
               psinablapsi_paw(2,1:3,jb,ib)=-psinablapsi_paw(2,1:3,ib,jb)
             end if
           end do ! ib

           if (iomode_etsf_mpiio.and.mpi_enreg%paral_spinor==1) then
             call timab(48,1,tsec)
             call xmpi_sum_master(psinablapsi_paw,master,spaceComm_spinor,ierr)
             !call xmpi_sum_master(psinablapsi_soc,master,spaceComm_spinor,ierr)
             call timab(48,2,tsec)
           end if

         end if ! myband

!        Write to OPT file in case of MPI-IO
         if (iomode_etsf_mpiio.and.i_am_master_spfft) then
           nc_start=[1,1,1,jb,ikpt,isppol];nc_stride=[1,1,1,1,1,1] 
           if (myband) then
             psinablapsi=psinablapsi+psinablapsi_paw
             if (dtset%pawspnorb==1) psinablapsi=psinablapsi+psinablapsi_soc
             nc_count=[2,3,mband,1,1,1]
           else
             nc_count=[0,0,0,0,0,0]
           end if
#ifdef HAVE_NETCDF
           NCF_CHECK(nf90_put_var(ncid,varid,psinablapsi,start=nc_start,stride=nc_stride,count=nc_count))
#endif
         end if

       end do ! jb

       if (mkmem/=0) then
         ibg = ibg +       my_nspinor*nband_cprj_k
         icg = icg + npw_k*my_nspinor*nband_k
         ikg = ikg + npw_k
       end if

       if (cprj_paral_band) then
         call pawcprj_free(cprj_k)
         ABI_FREE(cprj_k)
       end if
       if (mkmem*nsppol/=1) then
         call pawcprj_free(cprj_k_loc)
         ABI_FREE(cprj_k_loc)
       end if
       ABI_FREE(kpg_k)

!      Write to OPT file if not MPI-IO

!      >>> Reduction in case of parallelism
       if (.not.iomode_etsf_mpiio) then
         call timab(48,1,tsec)
         call xmpi_sum_master(psinablapsi,master,spaceComm_bandspinorfft,ierr)
         call xmpi_sum_master(psinablapsi_paw,master,spaceComm_bandspinor,ierr)
         call timab(48,2,tsec)
         psinablapsi=psinablapsi+psinablapsi_paw
         if (dtset%pawspnorb==1) then
           call xmpi_sum_master(psinablapsi_soc,master,spaceComm_band,ierr)
           psinablapsi=psinablapsi+psinablapsi_soc
         end if
       end if

!      >>> This my kpt and I am the master node: I write the data    
       if (.not.iomode_etsf_mpiio) then
         if (i_am_master) then
           if (iomode==IO_MODE_ETSF) then
             nc_stride=[1,1,1,1,1,1] 
             if (nc_unlimited) then
               nc_start=[1,1,1,ikpt,isppol,1] ; nc_count=[2,3,mband,1,1,mband]
             else
               nc_start=[1,1,1,1,ikpt,isppol] ; nc_count=[2,3,mband,mband,1,1]
             end if
#ifdef HAVE_NETCDF
             NCF_CHECK(nf90_put_var(ncid,varid,psinablapsi,start=nc_start,stride=nc_stride,count=nc_count))
#endif
           else
             write(ount)((psinablapsi(1:2,1,ib,jb),ib=1,nband_k),jb=1,nband_k)
             write(ount)((psinablapsi(1:2,2,ib,jb),ib=1,nband_k),jb=1,nband_k)
             write(ount)((psinablapsi(1:2,3,ib,jb),ib=1,nband_k),jb=1,nband_k)
           end if

!        >>> This my kpt and I am not the master node: I send the data    
         else if (i_am_master_band.and.i_am_master_spfft) then
           if (mpi_enreg%me_kpt/=master_spfftband) then
             ABI_BUG('Problem with band communicator!')
           end if
           call xmpi_exch(psinablapsi,etiq,mpi_enreg%me_kpt,psinablapsi,master,spaceComm_kpt,ierr)
         end if
       end if  

!    >>> This is not my kpt and I am the master node: I receive the data and I write    
     elseif ((.not.iomode_etsf_mpiio).and.i_am_master) then ! mykpt
       sender=master_spfftband
       call xmpi_exch(psinablapsi,etiq,sender,psinablapsi,master,spaceComm_kpt,ierr)
       if (iomode==IO_MODE_ETSF) then
         nc_stride=[1,1,1,1,1,1] 
         if (nc_unlimited) then
           nc_start=[1,1,1,ikpt,isppol,1] ; nc_count=[2,3,mband,1,1,mband]
         else
           nc_start=[1,1,1,1,ikpt,isppol] ; nc_count=[2,3,mband,mband,1,1]
         end if
#ifdef HAVE_NETCDF
         NCF_CHECK(nf90_put_var(ncid,varid,psinablapsi,start=nc_start,stride=nc_stride,count=nc_count))
#endif
       else
         write(ount)((psinablapsi(1:2,1,ib,jb),ib=1,nband_k),jb=1,nband_k)
         write(ount)((psinablapsi(1:2,2,ib,jb),ib=1,nband_k),jb=1,nband_k)
         write(ount)((psinablapsi(1:2,3,ib,jb),ib=1,nband_k),jb=1,nband_k)
       end if
     end if ! mykpt

     bdtot_index=bdtot_index+nband_k

!    End loop on spin,kpt
   end do ! ikpt
 end do !isppol

!Close file
 if (i_am_master.or.(iomode_etsf_mpiio.and.i_am_master_spfft)) then
   if (iomode==IO_MODE_ETSF) then
#ifdef HAVE_NETCDF
     NCF_CHECK(nf90_close(ncid))
#endif
   else
     ierr=close_unit(ount,msg)
     ABI_CHECK(ierr==0,"Error while closing OPT file")
   end if
 end if

!Datastructures deallocations
 ABI_FREE(psinablapsi)
 ABI_FREE(psinablapsi_paw)
 if (dtset%pawspnorb==1) then
   ABI_FREE(psinablapsi_soc)
   do iatom=1,natom
     ABI_FREE(phisocphj(iatom)%value)
   end do
   ABI_FREE(phisocphj)
 end if
 if (.not.already_has_nabla) then
   do itypat=1,dtset%ntypat
     if (allocated(pawtab(itypat)%nabla_ij)) then
       ABI_FREE(pawtab(itypat)%nabla_ij)
       pawtab(itypat)%has_nabla=0
     end if
   end do
 end if

 DBG_EXIT("COLL")

 end subroutine optics_paw
!!***

!----------------------------------------------------------------------

!!****f* m_paw_optics/optics_paw_core
!! NAME
!! optics_paw_core
!!
!! FUNCTION
!! Compute matrix elements need for X spectr. (in the PAW context) and store them in a file
!!  Matrix elements = <Phi_core|Nabla|Phi_j>
!!
!! COPYRIGHT
!! Copyright (C) 2005-2022 ABINIT group (SM,MT,NB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  cprj(natom,mcprj)= <p_lmn|Cnk> coefficients for each WF |Cnk> and each |p_lmn> non-local projector
!!  dimcprj(natom)=array of dimensions of array cprj (not ordered)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  filpsp(ntypat)=name(s) of the pseudopotential file(s)
!!  mband=maximum number of bands
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mkmem =number of k points treated by this node.
!!  mpi_enreg=information about MPI parallelization
!!  mpsang =1+maximum angular momentum for nonlocal pseudopotentials
!!  natom=number of atoms in cell.
!!  nkpt=number of k points.
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  pawang <type(pawang_type)>= PAW ANGular mesh discretization and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= PAW rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  znucl(ntypat)=atomic number of atom type
!!
!! OUTPUT
!!  (only writing in a file)
!!
!! PARENTS
!!      m_outscfcv
!!
!! CHILDREN
!!      destroy_mpi_enreg,hdr%free,hdr_io,hdr_read_from_fname,initmpi_seq
!!      kramerskronig,matrginv,metric
!!
!! SOURCE

 subroutine optics_paw_core(atindx1,cprj,dimcprj,dtfil,dtset,eigen0,filpsp,hdr,&
&               mband,mcprj,mkmem,mpi_enreg,mpsang,natom,nkpt,nsppol,&
&               pawang,pawrad,pawrhoij,pawtab,znucl)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcprj,mkmem,mpsang,natom,nkpt,nsppol
 type(MPI_type),intent(in) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: atindx1(natom),dimcprj(natom)
 character(len=fnlen),intent(in) :: filpsp(dtset%ntypat)
 real(dp),intent(in) :: eigen0(mband*nkpt*nsppol),znucl(dtset%ntypat)
 type(pawcprj_type),target,intent(inout) :: cprj(natom,mcprj)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
 type(pawrhoij_type),intent(inout) :: pawrhoij(mpi_enreg%my_natom)
 type(pawtab_type),target,intent(inout) :: pawtab(dtset%ntypat)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: bdtot_index,cplex,etiq,iatom,ic,ibg,idir
 integer :: ierr,ikpt,ilmn,iln,ount,is,my_jb
 integer :: iorder_cprj,ispinor,isppol,istwf_k,itypat
 integer :: jb,jbsp,jlmn,lmn_size,lmncmax,mband_cprj,ncid,varid
 integer :: me,my_nspinor,nband_cprj_k,option_core
 integer :: nband_k,nphicor,ncorespinor,sender,iomode,fformopt,master_spfftband
 integer :: spaceComm_band,spaceComm_bandspinorfft,spaceComm_fft,spaceComm_kpt
 integer :: spaceComm_spinor,spaceComm_bandspinor,spaceComm_spinorfft,spaceComm_w
 logical :: already_has_nabla,cprj_paral_band,ex,mykpt,myband
 logical :: iomode_etsf_mpiio,abinitcorewf,use_spinorbit,xmlcorewf
 logical :: i_am_master,i_am_master_band,i_am_master_spfft
 character(len=fnlen) :: filecore
 real(dp) :: cpnm1,cpnm2
 character(len=500) :: msg
!arrays
 integer :: nc_count(7),nc_start(7),nc_stride(7),tmp_shape(3)
 integer,allocatable :: indlmn_cor(:,:),lcor(:),ncor(:),kappacor(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: energy_cor(:),phi_cor(:,:)
 real(dp),allocatable :: psinablapsi(:,:,:,:,:),psinablapsi_soc(:,:,:,:,:)
 real(dp),pointer :: soc_ij(:,:,:)
 type(coeff5_type),allocatable,target :: phisocphj(:)
 type(pawcprj_type),pointer :: cprj_k(:,:),cprj_k_loc(:,:)
 type(nctkdim_t) :: ncdims(3)
 type(nctkarr_t) :: nctk_arrays(5)

! ************************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 msg="mkmem==0 not supported anymore!"
 ABI_CHECK(mkmem/=0,msg)
!Probably should check for spinor parallelism, because thats likely to not work correctly
 msg="Spinor parallelism not implemented for optics_paw_core!"
 ABI_CHECK(dtset%npspinor==1,msg)
!Is mpi_enreg initialized?
 if (xmpi_paral==1) then
   tmp_shape = shape(mpi_enreg%proc_distrb)
   if (nkpt > tmp_shape(1)) then
     ABI_BUG('problem with proc_distrb!')
   end if
 end if

!Init parallelism
 spaceComm_w=mpi_enreg%comm_cell
 if (mpi_enreg%paral_kgb==1) then
   spaceComm_kpt=mpi_enreg%comm_kpt
   spaceComm_fft=mpi_enreg%comm_fft
   spaceComm_band=mpi_enreg%comm_band
   spaceComm_spinor=mpi_enreg%comm_spinor
   spaceComm_bandspinor=mpi_enreg%comm_bandspinor
   spaceComm_spinorfft=mpi_enreg%comm_spinorfft
   spaceComm_bandspinorfft=mpi_enreg%comm_bandspinorfft
 else
   spaceComm_kpt=mpi_enreg%comm_kpt
   spaceComm_fft=xmpi_comm_self
   spaceComm_band=mpi_enreg%comm_band
   spaceComm_spinor=xmpi_comm_self
   spaceComm_bandspinor=spaceComm_band
   spaceComm_spinorfft=xmpi_comm_self
   spaceComm_bandspinorfft=xmpi_comm_self
 end if
 me=xmpi_comm_rank(spaceComm_w)
 i_am_master=(me==master)
 i_am_master_band=(xmpi_comm_rank(spaceComm_band)==master)
 i_am_master_spfft=(xmpi_comm_rank(spaceComm_spinorfft)==master)
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

!------------------------------------------------------------------------------------------------
!1- Reading of core wavefunctions
!------------------------------------------------------------------------------------------------

!Note: core WF is read for itypat=1
!At present, impose 2-spinor simulataneously for core anf valence WF
 ncorespinor=merge(2,1,dtset%pawspnorb==1.or.dtset%nspinor==2)
 filecore=trim(filpsp(1)) ; iln=len(trim(filecore))
 abinitcorewf=.false. ; if (iln>3) abinitcorewf=(filecore(iln-6:iln)=='.abinit')
 xmlcorewf=.false. ; if (iln>3) xmlcorewf=(filecore(iln-3:iln)=='.xml')
 if ((.not.xmlcorewf).and.(.not.abinitcorewf)) filecore=filecore(1:iln)//'.corewf'
 if (abinitcorewf) filecore=filecore(1:iln-6)//'corewf.abinit'
 if (xmlcorewf) filecore=filecore(1:iln-3)//'corewf.xml'
 inquire(file=filecore,exist=ex)

!Relativistic case
 if (ncorespinor==2) then
   if (ex) then
     !Use <filepsp>.corewf.xml or <filepsp>.corewf.abinit
     call pawpsp_read_corewf(energy_cor,indlmn_cor,lcor,lmncmax,ncor,nphicor,pawrad(1),phi_cor,&
&                            filename=filecore,kappacor=kappacor)
   else
     !Use default name
     call pawpsp_read_corewf(energy_cor,indlmn_cor,lcor,lmncmax,ncor,nphicor,pawrad(1),phi_cor,&
&                            kappacor=kappacor)
   end if

!Non-relativistic case
 else
   if (ex) then
     !Use <filepsp>.corewf.xml or <filepsp>.corewf.abinit
     call pawpsp_read_corewf(energy_cor,indlmn_cor,lcor,lmncmax,ncor,nphicor,pawrad(1),phi_cor,&
&                            filename=filecore)
   else
     !Use default name
     call pawpsp_read_corewf(energy_cor,indlmn_cor,lcor,lmncmax,ncor,nphicor,pawrad(1),phi_cor)
   end if
   ABI_MALLOC(kappacor,(nphicor))
   kappacor(:)=zero
 endif

!----------------------------------------------------------------------------------
!2- Computation of phipphj=<phi_i|nabla|phi_core>
!----------------------------------------------------------------------------------

 already_has_nabla=all(pawtab(:)%has_nabla==3)
 if (ncorespinor==2) then
   already_has_nabla=all(pawtab(:)%has_nabla==4)
!  Should check whether this would work with spinor parallelism
   call pawnabla_core_init(mpsang,dtset%ntypat,pawrad,pawtab,phi_cor,indlmn_cor,diracflag=1)
 else
   call pawnabla_core_init(mpsang,dtset%ntypat,pawrad,pawtab,phi_cor,indlmn_cor)
 endif

!Compute spin-orbit contributions if necessary
 use_spinorbit=(dtset%pawspnorb==1.and.dtset%userie/=111) ! For testing purpose
 if (use_spinorbit) then
   option_core=1
   call pawnabla_soc_init(phisocphj,option_core,dtset%ixc,mpi_enreg%my_natom,natom,&
&       dtset%nspden,dtset%ntypat,pawang,pawrad,pawrhoij,pawtab,dtset%pawxcdev,&
&       dtset%spnorbscl,dtset%typat,dtset%xc_denpos,znucl,&
&       phi_cor=phi_cor,indlmn_cor=indlmn_cor,&
&       comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
 end if

!----------------------------------------------------------------------------------
!3- Opening of OPT2 file and header writing
!----------------------------------------------------------------------------------

!(master proc only)
 if (i_am_master) then
!  ====> NETCDF format
   if (dtset%iomode==IO_MODE_ETSF) then
     iomode=IO_MODE_ETSF
     fformopt=611 
#ifdef HAVE_NETCDF
!    Open/create nc file
     NCF_CHECK(nctk_open_create(ncid,nctk_ncify(dtfil%fnameabo_app_opt2),xmpi_comm_self))
!    Write header data
     NCF_CHECK(hdr%ncwrite(ncid,fformopt,nc_define=.true.))
!    Define additional dimensions
     ncdims(1)%name="number_of_core_states"
     ncdims(1)%value=nphicor
     ncdims(2)%name="number_of_core_spinor_components"
     ncdims(2)%value=ncorespinor
     ncdims(3)%name="number_of_core_spins"
     ncdims(3)%value=3-ncorespinor
     NCF_CHECK(nctk_def_dims(ncid,ncdims))
!    Define additional variables
     nctk_arrays(1)%name="eigenvalues_core"
     nctk_arrays(1)%dtype="dp"
     nctk_arrays(1)%shape_str="number_of_core_states"
     nctk_arrays(2)%name="dipole_core_valence"
     nctk_arrays(2)%dtype="dp"
     nctk_arrays(2)%shape_str=&
&     "complex, number_of_cartesian_directions,number_of_core_states,"// &
&     "number_of_atoms,max_number_of_states,number_of_kpoints,number_of_spins"
     nctk_arrays(3)%name="n_quantum_number_core"
     nctk_arrays(3)%dtype="int"
     nctk_arrays(3)%shape_str="number_of_core_states"
     nctk_arrays(4)%name="l_quantum_number_core"
     nctk_arrays(4)%dtype="int"
     nctk_arrays(4)%shape_str="number_of_core_states"
     nctk_arrays(5)%name="kappa_core"
     nctk_arrays(5)%dtype="int"
     nctk_arrays(5)%shape_str="number_of_core_states"
     NCF_CHECK(nctk_def_arrays(ncid, nctk_arrays))
     NCF_CHECK(nctk_set_atomic_units(ncid, "eigenvalues_core"))
     NCF_CHECK(nctk_set_atomic_units(ncid, "dipole_core_valence"))
!    Write core states
     NCF_CHECK(nctk_set_datamode(ncid))
     varid=nctk_idname(ncid,"eigenvalues_core")
     NCF_CHECK(nf90_put_var(ncid,varid,energy_cor))
     varid=nctk_idname(ncid,"n_quantum_number_core")
     NCF_CHECK(nf90_put_var(ncid,varid,ncor))
     varid=nctk_idname(ncid,"l_quantum_number_core")
     NCF_CHECK(nf90_put_var(ncid,varid,lcor))
     varid=nctk_idname(ncid,"kappa_core")
     NCF_CHECK(nf90_put_var(ncid,varid,kappacor))
!    Write eigenvalues
     varid=nctk_idname(ncid,"eigenvalues")
     NCF_CHECK(nf90_put_var(ncid,varid,reshape(eigen0,[mband,nkpt,nsppol])))
     !Close file here because the rest has possibly to be written with collective I/O
     NCF_CHECK(nf90_close(ncid))
#else
     msg = "In order to use iomode=3 and prtnabla>0 together, NetCDF support must be enabled!"
     ABI_ERROR(msg)
#endif
!  ====> Standard FORTRAN binary file format
   else
     iomode=IO_MODE_FORTRAN_MASTER
     fformopt=612  ! MT 12sept21: change the OPT2 Fortran file format
     if (2*nphicor*natom*mband>2**30) fformopt=613 ! Format for large file records
     if (open_file(dtfil%fnameabo_app_opt2,msg,newunit=ount,form="unformatted",status="unknown")/= 0) then
       ABI_ERROR(msg)
     end if
     call hdr%fort_write(ount,fformopt,ierr,rewind=.true.)
     write(ount)(eigen0(jb),jb=1,mband*nkpt*nsppol)
     write(ount) nphicor
     do iln=1,nphicor
       write(ount) ncor(iln),lcor(iln),kappacor(iln),energy_cor(iln)
     end do
   end if ! File format
 end if ! master node
 call xmpi_bcast(iomode,master,spaceComm_w,ierr)
 call xmpi_bcast(fformopt,master,spaceComm_w,ierr)
 iomode_etsf_mpiio=(iomode==IO_MODE_ETSF.and.nctk_has_mpiio.and.use_netcdf_mpiio)

 ABI_FREE(ncor)
 ABI_FREE(lcor)
 ABI_FREE(kappacor)
 ABI_FREE(phi_cor)
 ABI_FREE(energy_cor)

!----------------------------------------------------------------------------------
!4- Computation of <psi_n|p_i>(<phi_i|-i.nabla|phi_core>)
!----------------------------------------------------------------------------------

!Prepare core-valence dipoles writing
!In case of netCDF access to OPT2 file, prepare collective I/O
 if (iomode == IO_MODE_ETSF) then
   if (iomode_etsf_mpiio) then
     if (i_am_master_spfft) then
       NCF_CHECK(nctk_open_modify(ncid,nctk_ncify(dtfil%fnameabo_app_opt2),spaceComm_band))
       varid=nctk_idname(ncid,"dipole_core_valence")
       if (xmpi_comm_size(spaceComm_w)>1) then
         NCF_CHECK(nctk_set_collective(ncid,varid))
       end if
       NCF_CHECK(nctk_set_datamode(ncid))
     end if
   else if (i_am_master) then
     NCF_CHECK(nctk_open_modify(ncid,nctk_ncify(dtfil%fnameabo_app_opt2),xmpi_comm_self))
     varid=nctk_idname(ncid,"dipole_core_valence")
     if (nctk_has_mpiio.and.(.not.use_netcdf_mpiio)) then
       NCF_CHECK(nctk_set_collective(ncid,varid))
     end if
     NCF_CHECK(nctk_set_datamode(ncid))
   end if     
 end if
 if (iomode_etsf_mpiio) then
   !If MPI-IO, store only elements for one band
   ABI_MALLOC(psinablapsi,(2,3,nphicor,natom,1))
   if (use_spinorbit) then
     ABI_MALLOC(psinablapsi_soc,(2,3,nphicor,natom,1))
   end if
 else
   !If not, store the elements for all bands
   ABI_MALLOC(psinablapsi,(2,3,nphicor,natom,mband))
   if (use_spinorbit) then
     ABI_MALLOC(psinablapsi_soc,(2,3,nphicor,natom,mband))
   end if
 end if

!Determine if cprj datastructure is distributed over bands
 mband_cprj=mcprj/(my_nspinor*mkmem*nsppol)
 cprj_paral_band=(mband_cprj<mband)

!LOOP OVER SPINS
 ibg=0
 bdtot_index=0
 do isppol=1,nsppol

!  LOOP OVER k POINTS
   do ikpt=1,nkpt

     etiq=ikpt+(isppol-1)*nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
     master_spfftband=minval(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol))

!    Select k-points for current proc
     mykpt=.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,mpi_enreg%me_kpt))
     if (mykpt) then

!      Data depending on k-point
       istwf_k=dtset%istwfk(ikpt)
       cplex=2;if (istwf_k>1) cplex=1

!      Extract cprj for this k-point
       nband_cprj_k=nband_k;if (cprj_paral_band) nband_cprj_k=nband_k/mpi_enreg%nproc_band
       if (mkmem*nsppol/=1) then
         iorder_cprj=0
         ABI_MALLOC(cprj_k_loc,(natom,my_nspinor*nband_cprj_k))
         call pawcprj_alloc(cprj_k_loc,0,dimcprj)
         call pawcprj_get(atindx1,cprj_k_loc,cprj,natom,1,ibg,ikpt,iorder_cprj,isppol,&
&         mband_cprj,mkmem,natom,nband_cprj_k,nband_cprj_k,my_nspinor,nsppol,dtfil%unpaw,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
       else
         cprj_k_loc => cprj
       end if

!      if cprj are distributed over bands, gather them (because we need to mix bands)
       if (cprj_paral_band) then
         ABI_MALLOC(cprj_k,(natom,my_nspinor*nband_k))
         call pawcprj_alloc(cprj_k,0,dimcprj)
         call pawcprj_mpi_allgather(cprj_k_loc,cprj_k,natom,my_nspinor*nband_cprj_k, &
&                                   my_nspinor*mpi_enreg%bandpp,&
&         dimcprj,0,mpi_enreg%nproc_band,mpi_enreg%comm_band,ierr,rank_ordered=.false.)
       else
         cprj_k => cprj_k_loc
       end if

!      Loops over bands
       do jb=1,nband_k
         !If MPI-IO, store only ib elements for each jb ; if not, store all (ib,jb) pairs
         my_jb=merge(1,jb,iomode_etsf_mpiio)

         psinablapsi(:,:,:,:,my_jb)=zero
         if (use_spinorbit) psinablapsi_soc(:,:,:,:,my_jb)=zero

!        Computation of <psi_n|p_i><phi_i|-i.nabla|phi_core>
!        ----------------------------------------------------------------------------------

!        Select bands for current proc
         myband=.true.
         if (mpi_enreg%paral_kgb==1) then
           myband=(mod(jb-1,mpi_enreg%nproc_band)==mpi_enreg%me_band)
         else if (xmpi_paral==1) then
           myband=(abs(mpi_enreg%proc_distrb(ikpt,jb,isppol)-mpi_enreg%me_kpt)==0)
         end if
         if (myband) then

           jbsp=(jb-1)*my_nspinor

!          1-spinor case
           if (dtset%nspinor==1.and.ncorespinor==1) then
             jbsp=jbsp+1
             if (cplex==1) then ! Real WF case
               do iatom=1,natom
                 itypat=dtset%typat(iatom)
                 lmn_size=pawtab(itypat)%lmn_size
                 do jlmn=1,lmn_size
                   do ilmn=1,lmncmax
                     ic=indlmn_cor(5,ilmn)
                     cpnm1=cprj_k(iatom,jbsp)%cp(1,jlmn)
                     psinablapsi(2,:,ic,iatom,my_jb)=psinablapsi(2,:,ic,iatom,my_jb) &
&                        -cpnm1*pawtab(itypat)%nabla_ij(:,jlmn,ilmn)
                   end do !ilmn
                 end do !jlmn
               end do !iatom
             else ! Complex WF case
               do iatom=1,natom
                 itypat=dtset%typat(iatom)
                 lmn_size=pawtab(itypat)%lmn_size
                 do jlmn=1,lmn_size
                   do ilmn=1,lmncmax
                     ic=indlmn_cor(5,ilmn)
                     cpnm1=cprj_k(iatom,jbsp)%cp(1,jlmn)
                     cpnm2=cprj_k(iatom,jbsp)%cp(2,jlmn)
                     psinablapsi(1,:,ic,iatom,my_jb)=psinablapsi(1,:,ic,iatom,my_jb) &
&                        -cpnm2*pawtab(itypat)%nabla_ij(:,jlmn,ilmn)
                     psinablapsi(2,:,ic,iatom,my_jb)=psinablapsi(2,:,ic,iatom,my_jb) &
&                        -cpnm1*pawtab(itypat)%nabla_ij(:,jlmn,ilmn)
                   end do !ilmn
                 end do !jlmn
               end do !iatom
             end if

!          2-spinor case
           else if (dtset%nspinor==2.and.ncorespinor==2) then
             do ispinor=1,my_nspinor
               jbsp=jbsp+1
               do iatom=1,natom
                 itypat=dtset%typat(iatom)
                 lmn_size=pawtab(itypat)%lmn_size
                 do jlmn=1,lmn_size
                   do ilmn=1,lmncmax
                     is=indlmn_cor(6,ilmn)
                     if (modulo(jbsp,2)==modulo(is,2)) then ! Nabla is a spin-diagonal operator
                       ic=indlmn_cor(5,ilmn)
                       if (ic>0) then
                         cpnm1=cprj_k(iatom,jbsp)%cp(1,jlmn)
                         cpnm2=cprj_k(iatom,jbsp)%cp(2,jlmn)
                         psinablapsi(1,:,ic,iatom,my_jb)=psinablapsi(1,:,ic,iatom,my_jb) &
&                            +cpnm1*pawtab(itypat)%nabla_im_ij(:,jlmn,ilmn) &
&                            -cpnm2*pawtab(itypat)%nabla_ij(:,jlmn,ilmn)
                         psinablapsi(2,:,ic,iatom,my_jb)=psinablapsi(2,:,ic,iatom,my_jb) &
&                            -cpnm1*pawtab(itypat)%nabla_ij(:,jlmn,ilmn) &
&                            -cpnm2*pawtab(itypat)%nabla_im_ij(:,jlmn,ilmn)
                       end if
                     end if
                   end do ! ilmn
                 end do !jlmn
               end do !iatom
             end do !ispinor
           else
             msg="Core and valence WF should have the same spinor representation!"
             ABI_BUG(msg)
             !N. Brouwer initial coding: should be justified!
             !if (dtset%nspinor==1.and.ncorespinor==2) then  ! Average core spinors
             !  psinablapsi(1,:,ic,iatom,my_jb)=psinablapsi(1,:,ic,iatom,my_jb) &
             ! &    +half_sqrt2*(cpnm1*pawtab(itypat)%nabla_im_ij(:,jlmn,ilmn) &
             ! &                +cpnm2*pawtab(itypat)%nabla_ij(:,jlmn,ilmn))
             !  psinablapsi(2,:,ic,iatom,my_jb)=psinablapsi(2,:,ic,iatom,my_jb) &
             ! &    +half_sqrt2*(cpnm2*pawtab(itypat)%nabla_im_ij(:,jlmn,ilmn) &
             ! &                -cpnm1*pawtab(itypat)%nabla_ij(:,jlmn,ilmn))
             !else if (dtset%nspinor==2.and.ncorespinor==1) then ! Average valence spinors
             !  psinablapsi(1,:,ic,iatom,my_jb)=psinablapsi(1,:,ic,iatom,my_jb) &
             ! &    +half_sqrt2*cpnm2*pawtab(itypat)%nabla_ij(:,jlmn,ilmn)
             !  psinablapsi(2,:,ic,iatom,my_jb)=psinablapsi(2,:,ic,iatom,my_jb) &
             ! &    -half_sqrt2*cpnm1*pawtab(itypat)%nabla_ij(:,jlmn,ilmn)
             !endif
           endif

!          Spin-orbit coupling contribution:
!             Sum_i,ss'[<psi^s_n|p_i><phi_i|1/4 Alpha^2 (Sigma^ss' X dV/dr)|phj_core^s'>]
           if (use_spinorbit) then
!            Add: Sum_i,ss'[<Psi^s_n|p_i> (Sigma^ss' X g_ij_core^s']
!             where: g_ij_core^s = <Phi_i| 1/4 Alpha^2 dV/dr vec(r)/r) |Phj_core^s>
!            Note that:
!             if phi_cor_jlmn_cor(:) is up:
!              phisocphj(:)%value(re:im,1,idir,ilmn,jlmn_cor) is (Sigma^up-up X g_ij_core^up)
!              phisocphj(:)%value(re:im,2,idir,ilmn,jlmn_cor) is (Sigma^dn-up X g_ij_core^up)
!             if phi_cor_jlmn_cor(:) is down:
!              phisocphj(:)%value(re:im,1,idir,ilmn,jlmn_cor) is (Sigma^up-dn X g_ij_core^dn)
!              phisocphj(:)%value(re:im,2,idir,ilmn,jlmn_cor) is (Sigma^dn-dn X g_ij_core^dn)
!            Not compatible with parallelization over spinors
             jbsp=1+(jb-1)*dtset%nspinor
             do iatom=1,natom
               itypat=dtset%typat(iatom)
               lmn_size=pawtab(itypat)%lmn_size
               do jlmn=1,lmn_size
                 do ilmn=1,lmncmax
                   ic=indlmn_cor(5,ilmn)
                   if (ic>0) then
                     soc_ij => phisocphj(iatom)%value(:,:,:,jlmn,ilmn)
                     !Contribution from real part of <Psi^s_n|p_i>
                     cpnm1=cprj_k(iatom,jbsp  )%cp(1,jlmn)
                     cpnm2=cprj_k(iatom,jbsp+1)%cp(1,jlmn)
                     do idir=1,3
                       psinablapsi_soc(1,idir,ic,iatom,my_jb)=psinablapsi_soc(1,idir,ic,iatom,my_jb) &
&                             +soc_ij(1,1,idir)*cpnm1+soc_ij(1,2,idir)*cpnm2
                       psinablapsi_soc(2,idir,ic,iatom,my_jb)=psinablapsi_soc(2,idir,ic,iatom,my_jb) &
&                             +soc_ij(2,1,idir)*cpnm1+soc_ij(2,2,idir)*cpnm2
                     end do
                     !Contribution from imaginary part of <Psi^s_n|p_i>
                     cpnm1=cprj_k(iatom,jbsp  )%cp(2,jlmn)
                     cpnm2=cprj_k(iatom,jbsp+1)%cp(2,jlmn)
                     do idir=1,3
                       psinablapsi_soc(1,idir,ic,iatom,my_jb)=psinablapsi_soc(1,idir,ic,iatom,my_jb) &
&                             +soc_ij(2,1,idir)*cpnm1+soc_ij(2,2,idir)*cpnm2
                       psinablapsi_soc(2,idir,ic,iatom,my_jb)=psinablapsi_soc(2,idir,ic,iatom,my_jb) &
&                             -soc_ij(1,1,idir)*cpnm1-soc_ij(1,2,idir)*cpnm2
                     end do
                   end if
                 end do ! ilmn
               end do ! jlmn
             end do ! iatom
           end if ! use_spinorbit

           if (iomode_etsf_mpiio.and.mpi_enreg%paral_spinor==1) then
             call timab(48,1,tsec)
             call xmpi_sum_master(psinablapsi,master,spaceComm_spinor,ierr)
             call timab(48,2,tsec)
           end if

         end if ! myband

!        Write to OPT2 file in case of MPI-IO
         if (iomode_etsf_mpiio.and.i_am_master_spfft) then
           nc_start=[1,1,1,1,jb,ikpt,isppol];nc_stride=[1,1,1,1,1,1,1] 
           if (myband) then
             if (use_spinorbit) psinablapsi=psinablapsi+psinablapsi_soc
             nc_count=[2,3,nphicor,natom,1,1,1]
           else
             nc_count=[0,0,0,0,0,0,0]
           end if
#ifdef HAVE_NETCDF
           NCF_CHECK(nf90_put_var(ncid,varid,psinablapsi,start=nc_start,stride=nc_stride,count=nc_count))
#endif
         end if

       end do ! jb
       
       if (mkmem/=0) then
         ibg = ibg +  my_nspinor*nband_cprj_k
       end if

       if (cprj_paral_band) then
         call pawcprj_free(cprj_k)
         ABI_FREE(cprj_k)
       end if
       if (mkmem*nsppol/=1) then
         call pawcprj_free(cprj_k_loc)
         ABI_FREE(cprj_k_loc)
       end if

!      Write to OPT2 file if not MPI-IO

!      >>> Reduction in case of parallelism
       if (.not.iomode_etsf_mpiio) then
         call timab(48,1,tsec)
         call xmpi_sum_master(psinablapsi,master,spaceComm_bandspinor,ierr)
         call timab(48,2,tsec)
         if (use_spinorbit) then
           call xmpi_sum_master(psinablapsi_soc,master,spaceComm_band,ierr)
           psinablapsi=psinablapsi+psinablapsi_soc
         end if
       end if

!      >>> This my kpt and I am the master node: I write the data    
       if (.not.iomode_etsf_mpiio) then
         if (i_am_master) then
           if (iomode==IO_MODE_ETSF) then
             nc_start=[1,1,1,1,1,ikpt,isppol];nc_stride=[1,1,1,1,1,1,1] 
             nc_count=[2,3,nphicor,natom,mband,1,1]
#ifdef HAVE_NETCDF
             NCF_CHECK(nf90_put_var(ncid,varid,psinablapsi,start=nc_start,stride=nc_stride,count=nc_count))
#endif
           else
             if (fformopt==612) then ! New OPT2 file format
               write(ount) (((psinablapsi(1:2,1,ic,iatom,jb),ic=1,nphicor),iatom=1,natom),jb=1,nband_k)
               write(ount) (((psinablapsi(1:2,2,ic,iatom,jb),ic=1,nphicor),iatom=1,natom),jb=1,nband_k)
               write(ount) (((psinablapsi(1:2,3,ic,iatom,jb),ic=1,nphicor),iatom=1,natom),jb=1,nband_k)
             else if (fformopt==613) then ! Large OPT2 file format
               do jb=1,nband_k
                 write(ount) ((psinablapsi(1:2,1,ic,iatom,jb),ic=1,nphicor),iatom=1,natom)
                 write(ount) ((psinablapsi(1:2,2,ic,iatom,jb),ic=1,nphicor),iatom=1,natom)
                 write(ount) ((psinablapsi(1:2,3,ic,iatom,jb),ic=1,nphicor),iatom=1,natom)
               end do
             else ! Old OPT2 file format
               !The old writing was not efficient (indexes order is bad)
               do iatom=1,natom
                 write(ount) ((psinablapsi(1:2,1,ic,iatom,jb),jb=1,nband_k),ic=1,nphicor)
                 write(ount) ((psinablapsi(1:2,2,ic,iatom,jb),jb=1,nband_k),ic=1,nphicor)
                 write(ount) ((psinablapsi(1:2,3,ic,iatom,jb),jb=1,nband_k),ic=1,nphicor)
               end do
             end if
           end if

!        >>> This my kpt and I am not the master node: I send the data    
         else if (i_am_master_band.and.i_am_master_spfft) then
           if (mpi_enreg%me_kpt/=master_spfftband) then
             ABI_BUG('Problem with band communicator!')
           end if
           call xmpi_exch(psinablapsi,etiq,mpi_enreg%me_kpt,psinablapsi,master,spaceComm_kpt,ierr)
         end if
       end if  

!    >>> This is not my kpt and I am the master node: I receive the data and I write    
     elseif ((.not.iomode_etsf_mpiio).and.i_am_master) then ! mykpt
       sender=master_spfftband
       call xmpi_exch(psinablapsi,etiq,sender,psinablapsi,master,spaceComm_kpt,ierr)
       if (iomode==IO_MODE_ETSF) then
         nc_start=[1,1,1,1,1,ikpt,isppol];nc_stride=[1,1,1,1,1,1,1] 
         nc_count=[2,3,nphicor,natom,mband,1,1]
#ifdef HAVE_NETCDF
         NCF_CHECK(nf90_put_var(ncid,varid,psinablapsi,start=nc_start,stride=nc_stride,count=nc_count))
#endif
       else
         if (fformopt==612) then ! New OPT2 file format
           write(ount) (((psinablapsi(1:2,1,ic,iatom,jb),ic=1,nphicor),iatom=1,natom),jb=1,nband_k)
           write(ount) (((psinablapsi(1:2,2,ic,iatom,jb),ic=1,nphicor),iatom=1,natom),jb=1,nband_k)
           write(ount) (((psinablapsi(1:2,3,ic,iatom,jb),ic=1,nphicor),iatom=1,natom),jb=1,nband_k)
         else if (fformopt==613) then ! Large OPT2 file format
           do jb=1,nband_k
             write(ount) ((psinablapsi(1:2,1,ic,iatom,jb),ic=1,nphicor),iatom=1,natom)
             write(ount) ((psinablapsi(1:2,2,ic,iatom,jb),ic=1,nphicor),iatom=1,natom)
             write(ount) ((psinablapsi(1:2,3,ic,iatom,jb),ic=1,nphicor),iatom=1,natom)
           end do
         else ! Old OPT2 file format
           !The old writing was not efficient (indexes order is bad)
           do iatom=1,natom
             write(ount) ((psinablapsi(1:2,1,ic,iatom,jb),jb=1,nband_k),ic=1,nphicor)
             write(ount) ((psinablapsi(1:2,2,ic,iatom,jb),jb=1,nband_k),ic=1,nphicor)
             write(ount) ((psinablapsi(1:2,3,ic,iatom,jb),jb=1,nband_k),ic=1,nphicor)
           end do
         end if
       end if
     end if ! mykpt

     bdtot_index=bdtot_index+nband_k

!    End loop on spin,kpt
   end do ! ikpt
 end do !isppol

!Close file
 if (i_am_master.or.(iomode_etsf_mpiio.and.i_am_master_spfft)) then
   if (iomode==IO_MODE_ETSF) then
#ifdef HAVE_NETCDF
     NCF_CHECK(nf90_close(ncid))
#endif
   else
     ierr=close_unit(ount,msg)
     ABI_CHECK(ierr==0,"Error while closing OPT2 file")
   end if
 end if

!Datastructures deallocations
 ABI_FREE(indlmn_cor)
 ABI_FREE(psinablapsi)
 if (use_spinorbit) then
   ABI_FREE(psinablapsi_soc)
   do iatom=1,natom
     ABI_FREE(phisocphj(iatom)%value)
   end do
   ABI_FREE(phisocphj)
 end if
 if (.not.already_has_nabla) then
   do itypat=1,dtset%ntypat
     if (allocated(pawtab(itypat)%nabla_ij)) then
       ABI_FREE(pawtab(itypat)%nabla_ij)
       pawtab(itypat)%has_nabla=0
     end if
     if (allocated(pawtab(itypat)%nabla_im_ij)) then
       ABI_FREE(pawtab(itypat)%nabla_im_ij)
     end if
   end do
 end if

 DBG_EXIT("COLL")

 end subroutine optics_paw_core
!!***

!----------------------------------------------------------------------

!!****f* m_paw_optics/linear_optics_paw
!! NAME
!! linear_optics_paw
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! linear susceptiblity using matrix elements <-i Nabla> obtained from a
!! PAW ground state calculation. It uses formula 17 from Gadoc et al,
!! Phys. Rev. B 73, 045112 (2006) [[cite:Gajdo2006]] together with a scissors correction. It uses
!! a Kramers-Kronig transform to compute the real part from the imaginary part, and
!! it will work on all types of unit cells. It outputs all tensor elements of
!! both the real and imaginary parts.
!!
!! INPUTS
!!  filnam: base of file names to read data from
!!  mpi_enreg: mpi set up variable, not used in this code
!!
!! OUTPUT
!!  _real and _imag output files
!!
!! NOTES
!!  This routine is not tested
!!
!! PARENTS
!!      conducti
!!
!! CHILDREN
!!      hdr%free,hdr_io,hdr_read_from_fname,initmpi_seq
!!      kramerskronig,matrginv,metric,wffopen
!!
!! SOURCE

 subroutine linear_optics_paw(filnam,filnam_out)

!Arguments -----------------------------------
!scalars
 character(len=fnlen),intent(in) :: filnam,filnam_out

!Local variables-------------------------------
 integer,parameter :: master=0
 integer :: iomode,bantot,bdtot_index,fform1,headform
 integer :: iband,ierr,ii,ikpt,iom,iout,isppol,isym,jband,jj,me,mband
 integer :: method,mom,nband_k,nkpt,nspinor,nsppol,nsym,occopt,only_check
 integer :: rdwr,spaceComm,inpunt,reunt,imunt,wfunt
 integer,allocatable :: nband(:),symrel(:,:,:)
 real(dp) :: del,dom,fij,gdelta,omin,omax,paijpbij(2),mbpt_sciss,wij,ucvol
 real(dp) :: diffwp, diffwm
 real(dp) :: e2rot(3,3),gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),rprimdinv(3,3),symd(3,3),symdinv(3,3)
 real(dp),allocatable :: e1(:,:,:),e2(:,:,:,:),epsilon_tot(:,:,:,:),eigen0(:),eig0_k(:)
 real(dp),allocatable :: kpts(:,:),occ(:),occ_k(:),oml1(:),wtk(:)
 complex(dpc),allocatable :: eps_work(:)
 character(len=fnlen) :: filnam1,filnam_gen
 character(len=500) :: msg
 type(hdr_type) :: hdr
 type(wffile_type) :: wff1
!arrays
 real(dp),allocatable :: psinablapsi(:,:,:,:)

! *********************************************************************************

 DBG_ENTER("COLL")

!write(std_out,'(a)')' Give the name of the output file ...'
!read(std_in, '(a)') filnam_out
!write(std_out,'(a)')' The name of the output file is :',filnam_out

!Read data file
 if (open_file(filnam,msg,newunit=inpunt,form='formatted') /= 0 ) then
   ABI_ERROR(msg)
 end if

 rewind(inpunt)
 read(inpunt,*)
 read(inpunt,'(a)')filnam_gen       ! generic name for the files
 filnam1=trim(filnam_gen)//'_OPT' ! nabla matrix elements file

!Open the Wavefunction and optic files
!These default values are typical of sequential use
 iomode=IO_MODE_FORTRAN ; spaceComm=xmpi_comm_self; me=0

! Read the header of the optic files
 call hdr_read_from_fname(hdr, filnam1, fform1, spaceComm)
 call hdr%free()
 if (fform1 /= 610) then
   ABI_ERROR("Abinit8 requires an OPT file with fform = 610")
 end if

!Open the conducti optic files
 wfunt = get_unit()
 call WffOpen(iomode,spaceComm,filnam1,ierr,wff1,master,me,wfunt)

!Read the header from Ground state file
 rdwr=1
 call hdr_io(fform1,hdr,rdwr,wff1)

!Extract info from the header
 headform=hdr%headform
 bantot=hdr%bantot
 nkpt=hdr%nkpt
 ABI_MALLOC(kpts,(3,nkpt))
 ABI_MALLOC(wtk,(nkpt))
 kpts(:,:)=hdr%kptns(:,:)
 wtk(:)=hdr%wtk(:)
 nspinor=hdr%nspinor
 nsppol=hdr%nsppol
 occopt=hdr%occopt
 rprimd(:,:)=hdr%rprimd(:,:)
 rprimdinv(:,:) = rprimd(:,:)
 call matrginv(rprimdinv,3,3) ! need the inverse of rprimd to symmetrize the tensors
 ABI_MALLOC(nband,(nkpt*nsppol))
 ABI_MALLOC(occ,(bantot))
 occ(1:bantot)=hdr%occ(1:bantot)
 nband(1:nkpt*nsppol)=hdr%nband(1:nkpt*nsppol)
 nsym=hdr%nsym
 ABI_MALLOC(symrel,(3,3,nsym))
 symrel(:,:,:)=hdr%symrel(:,:,:)

!Get mband, as the maximum value of nband(nkpt)
 mband=maxval(nband(:))

!get ucvol etc.
 iout = -1
 call metric(gmet,gprimd,iout,rmet,rprimd,ucvol)

 write(std_out,*)
 write(std_out,'(a,3f10.5,a)' )' rprimd(bohr)      =',rprimd(1:3,1)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,2)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,3)
 write(std_out,*)
 write(std_out,'(a,3f10.5,a)' )' rprimdinv         =',rprimdinv(1:3,1)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimdinv(1:3,2)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimdinv(1:3,3)
 write(std_out,'(a,2i8)')      ' nkpt,mband        =',nkpt,mband

!get eigen0
 ABI_MALLOC(eigen0,(mband*nkpt*nsppol))
 read(wfunt)(eigen0(iband),iband=1,mband*nkpt*nsppol)

 read(inpunt,*)mbpt_sciss
 read(inpunt,*)dom,omin,omax,mom
 close(inpunt)

 ABI_MALLOC(oml1,(mom))
 ABI_MALLOC(e1,(3,3,mom))
 ABI_MALLOC(e2,(2,3,3,mom))
 ABI_MALLOC(epsilon_tot,(2,3,3,mom))
 ABI_MALLOC(eps_work,(mom))
 del=(omax-omin)/(mom-1)
 do iom=1,mom
   oml1(iom)=omin+dble(iom-1)*del
 end do
 write(std_out,'(a,i8,4f10.5,a)')' npts,omin,omax,width,mbpt_sciss      =',mom,omin,omax,dom,mbpt_sciss,' Ha'

 ABI_MALLOC(psinablapsi,(2,3,mband,mband))

!loop over spin components
 do isppol=1,nsppol
   bdtot_index = 0
!  loop over k points
   do ikpt=1,nkpt
!
!    number of bands for this k point
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     ABI_MALLOC(eig0_k,(nband_k))
     ABI_MALLOC(occ_k,(nband_k))
!    eigenvalues for this k-point
     eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
!    occupation numbers for this k-point
     occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
!    values of -i*nabla matrix elements for this k point
     psinablapsi=zero
     read(wfunt)((psinablapsi(1:2,1,iband,jband),iband=1,nband_k),jband=1,nband_k)
     read(wfunt)((psinablapsi(1:2,2,iband,jband),iband=1,nband_k),jband=1,nband_k)
     read(wfunt)((psinablapsi(1:2,3,iband,jband),iband=1,nband_k),jband=1,nband_k)

!    occupation numbers for k-point
     occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
!    accumulate e2 for this k point, Eq. 17 from PRB 73, 045112 (2006) [[cite:Gajdo2006]]
     do iband = 1, nband_k
       do jband = 1, nband_k
         fij = occ_k(iband) - occ_k(jband) !occ number difference
         wij = eig0_k(iband) - eig0_k(jband) !energy difference
         if (abs(fij) > zero) then ! only consider states of differing occupation numbers
           do ii = 1, 3
             do jj = 1, 3
               paijpbij(1) = psinablapsi(1,ii,iband,jband)*psinablapsi(1,jj,iband,jband) + &
&               psinablapsi(2,ii,iband,jband)*psinablapsi(2,jj,iband,jband)
               paijpbij(2) = psinablapsi(2,ii,iband,jband)*psinablapsi(1,jj,iband,jband) - &
&               psinablapsi(1,ii,iband,jband)*psinablapsi(2,jj,iband,jband)
               do iom = 1, mom
!                original version
!                diffw = wij + mbpt_sciss - oml1(iom) ! apply scissors term here
!                gdelta = exp(-diffw*diffw/(4.0*dom*dom))/(2.0*dom*sqrt(pi)) ! delta fnc resolved as Gaussian
!                e2(1,ii,jj,iom) = e2(1,ii,jj,iom) - (4.0*pi*pi/ucvol)*wtk(ikpt)*fij*paijpbij(1)*gdelta/(oml1(iom)*oml1(iom))
!                e2(2,ii,jj,iom) = e2(2,ii,jj,iom) - (4.0*pi*pi/ucvol)*wtk(ikpt)*fij*paijpbij(2)*gdelta/(oml1(iom)*oml1(iom))
                 diffwm = wij - mbpt_sciss + oml1(iom) ! apply scissors term here
                 diffwp = wij + mbpt_sciss - oml1(iom) ! apply scissors term here
                 gdelta = exp(-diffwp*diffwp/(4.0*dom*dom))/(2.0*dom*sqrt(pi))
                 e2(1,ii,jj,iom) = e2(1,ii,jj,iom) - (4.0*pi*pi/ucvol)*wtk(ikpt)*fij*paijpbij(1)*gdelta/(wij*wij)
                 e2(2,ii,jj,iom) = e2(2,ii,jj,iom) - (4.0*pi*pi/ucvol)*wtk(ikpt)*fij*paijpbij(2)*gdelta/(wij*wij)
               end do ! end loop over spectral points
             end do ! end loop over jj = 1, 3
           end do ! end loop over ii = 1, 3
         end if ! end selection on fij /= 0
       end do ! end loop over jband
     end do ! end loop over iband

     ABI_FREE(eig0_k)
     ABI_FREE(occ_k)
     bdtot_index=bdtot_index+nband_k
   end do ! end loop over k points
 end do ! end loop over spin polarizations

!here apply nsym symrel transformations to reconstruct full tensor from IBZ part
 epsilon_tot(:,:,:,:) = zero
 do isym = 1, nsym
   symd(:,:)=matmul(rprimd(:,:),matmul(symrel(:,:,isym),rprimdinv(:,:)))
   symdinv(:,:)=symd(:,:)
   call matrginv(symdinv,3,3)
   do iom = 1, mom
     e2rot(:,:)=matmul(symdinv(:,:),matmul(e2(1,:,:,iom),symd(:,:)))
     epsilon_tot(2,:,:,iom) = epsilon_tot(2,:,:,iom)+e2rot(:,:)/nsym
   end do
 end do

!generate e1 from e2 via KK transforma
 method=0 ! use naive integration ( = 1 for simpson)
 only_check=0 ! compute real part of eps in kk routine
 do ii = 1, 3
   do jj = 1, 3
     eps_work(:) = cmplx(0.0,epsilon_tot(2,ii,jj,:), kind=dpc)
     call kramerskronig(mom,oml1,eps_work,method,only_check)
     epsilon_tot(1,ii,jj,:) = real(eps_work(:))
     if (ii /= jj) epsilon_tot(1,ii,jj,:) = epsilon_tot(1,ii,jj,:)- 1.0
   end do ! end loop over jj
 end do ! end loop over ii

 if (open_file(trim(filnam_out)//'_imag',msg,newunit=reunt,form='formatted') /= 0) then
   ABI_ERROR(msg)
 end if

 if (open_file(trim(filnam_out)//'_real',msg,unit=imunt,form='formatted') /= 0) then
   ABI_ERROR(msg)
 end if

 write(reunt,'(a12,6a13)')' # Energy/Ha ','eps_2_xx','eps_2_yy','eps_2_zz',&
& 'eps_2_yz','eps_2_xz','eps_2_xy'
 write(imunt,'(a12,6a13)')' # Energy/Ha ','eps_1_xx','eps_1_yy','eps_1_zz',&
& 'eps_1_yz','eps_1_xz','eps_1_xy'

 do iom = 1, mom
   write(reunt,'(ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4)') oml1(iom),' ',&
&   epsilon_tot(2,1,1,iom),' ',epsilon_tot(2,2,2,iom),' ',epsilon_tot(2,3,3,iom),' ',&
&   epsilon_tot(2,2,3,iom),' ',epsilon_tot(2,1,3,iom),' ',epsilon_tot(2,1,2,iom)
   write(imunt,'(ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4)') oml1(iom),' ',&
&   epsilon_tot(1,1,1,iom),' ',epsilon_tot(1,2,2,iom),' ',epsilon_tot(1,3,3,iom),' ',&
&   epsilon_tot(1,2,3,iom),' ',epsilon_tot(1,1,3,iom),' ',epsilon_tot(1,1,2,iom)
 end do

 close(reunt)
 close(imunt)

 ABI_FREE(nband)
 ABI_FREE(oml1)
 ABI_FREE(e2)
 ABI_FREE(e1)
 ABI_FREE(occ)
 ABI_FREE(psinablapsi)
 ABI_FREE(eigen0)
 ABI_FREE(wtk)
 ABI_FREE(kpts)

 call hdr%free()

 DBG_EXIT("COLL")

 end subroutine linear_optics_paw
!!***

!----------------------------------------------------------------------

!!****f* m_paw_optics/pawnabla_soc_init
!! NAME
!! pawnabla_soc_init
!!
!! FUNCTION
!! Compute the PAW SOC contribution(s) to the momentum PAW matrix elements,
!!  i.e. <Phi_i|1/4 Alpha^2 dV/dr (Sigma X vec(r)/r) |Phi_j>
!!   where:
!!    {Phi_i}= AE partial waves
!!    Alpha = inverse of fine structure constant
!!    Sigma^alpha,beta= Pauli matrices
!!    X = cross product
!!
!! There are 2 typical uses:
!!   - Valence-valence terms: Phi_i and Phi_j are PAW AE partial waves (unpolarized)
!!   - Core-valence terms: Phi_i is are AE partial waves and Phi_j are core AE wave-functions (spinors)
!!
!! In practice we compute:
!!  (Sigma^up-up X g_ij) and (Sigma^up-dn X g_ij)    (X = vector cross product)
!!  (Sigma^dn-up X g_ij) and (Sigma^dn-dn X g_ij)
!!   where:
!!    g_ij= 1/4 Alpha^2 Int_[Phi_i(r)/r Phi_j(r)/r dV(r)/dr r^2 dr] . Gvec_ij
!!        and Gvec_ij= Int[S_limi S_ljmj vec(r)/r dOmega] (Gaunt coefficients)
!!
!! COPYRIGHT
!! Copyright (C) 2021-2022 ABINIT group (NBrouwer,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ixc= choice of exchange-correlation scheme (see above, and below)
!!  my_natom=number of atoms treated by current processor
!!  natom=total number of atoms in cell
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms in unit cell.
!!  option_core=Type of calculation: 0=valence-valence, 1=core-valence
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  spnorbscl=scaling factor for spin-orbit coupling
!!  typat(natom) =Type of each atoms
!!  xc_denpos= lowest allowe density (usually for the computation of the XC functionals)
!!  znucl(ntypat)=gives the nuclear charge for all types of atoms
!!  [phi_cor(mesh_size,nphicor)]=--optional-- core wave-functions for the current type of atoms;
!!                               only needed when option_core=1
!!  [indlmn_cor(6,nlmn_core)]=--optional-- array giving l,m,n,lm,ln,s for i=lmn, for the core wave functions;
!!                            only needed when option_core=1
!!  [mpi_atmtab(:)]=--optional-- indexes of the atoms treated by current proc
!!  [comm_atom]=--optional-- MPI communicator over atoms
!!
!! OUTPUT
!    phisocphj(dtset%natom) <type(coeff4_type)>= stores soc coefficients:
!!!  If  option_core==0 or nspinor_cor==1:
!!     phisocphj(iat)%value(1,1,idir,ilmn,jlmn) is real part of (Sigma^up-up X g_ij)
!!     phisocphj(iat)%value(2,1,idir,ilmn,jlmn) is imaginary part of (Sigma^up-up X g_ij)
!!     phisocphj(iat)%value(1,2,idir,ilmn,jlmn) is real part of (Sigma^up-dn X g_ij)
!!     phisocphj(iat)%value(2,2,idir,ilmn,jlmn) is imaginary part of (Sigma^up-dn X g_ij)
!!   If option_core==1 and nspinor_cor==2 (core-valence with spinorial core WF):
!!     phisocphj(iat)%value(1,1,idir,ilmn,2*jlmn-1) is real part of (Sigma^up-up X g_ij^up)
!!     phisocphj(iat)%value(2,1,idir,ilmn,2*jlmn-1) is imaginary part of (Sigma^up-up X g_ij^up)
!!     phisocphj(iat)%value(1,2,idir,ilmn,2*jlmn-1) is real part of (Sigma^dn-up X g_ij^dn)
!!     phisocphj(iat)%value(2,2,idir,ilmn,2*jlmn-1) is imaginary part of (Sigma^dn-up X g_ij^dn)
!!     phisocphj(iat)%value(1,1,idir,ilmn,2*jlmn  ) is real part of (Sigma^up-dn X g_ij^up)
!!     phisocphj(iat)%value(2,1,idir,ilmn,2*jlmn  ) is imaginary part of (Sigma^up-dn X g_ij^up)
!!     phisocphj(iat)%value(1,2,idir,ilmn,2*jlmn  ) is real part of (Sigma^dn-dn X g_ij^dn)
!!     phisocphj(iat)%value(2,2,idir,ilmn,2*jlmn  ) is imaginary part of (Sigma^dn-dn X g_ij^dn)
!! (idir=cartesian direction)
!!
!! SIDE EFFECTS
!!
!! NOTES
!! If Phi_j is not polarized,
!!    (Sigma^dn-dn X g_ij)=-(Sigma^up-up X g_ij)
!!    (Sigma^dn-up X g_ij)= (Sigma^up-dn X g_ij)^*
!!    So, we store only 2 components.
!! If Phi_j is polarized, the spin component is included in the last dimension of
!!   phisocphj(iat)%value, i.e. lmn_size_cor=2*lmn_size
!!
!! PARENTS
!! optics_paw,optics_paw_core
!!
!! CHILDREN
!! nderiv_gen,simp_gen,pawrad_deducer0
!!
!! SOURCE

 subroutine pawnabla_soc_init(phisocphj,option_core,ixc,my_natom,natom,nspden,ntypat,pawang, &
&           pawrad,pawrhoij,pawtab,pawxcdev,spnorbscl,typat,xc_denpos,znucl, &
&           phi_cor,indlmn_cor,mpi_atmtab,comm_atom) ! Optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,my_natom,natom,nspden,ntypat,option_core,pawxcdev
 integer,optional,intent(in) :: comm_atom
 real(dp),intent(in) :: spnorbscl,xc_denpos
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in),target,optional :: indlmn_cor(:,:)
 integer,intent(in) :: typat(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in),target,optional :: phi_cor(:,:)
 real(dp),intent(in) :: znucl(ntypat)
 type(coeff5_type),allocatable,target,intent(inout) :: phisocphj(:)
 type(pawrad_type),target,intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)
 type(pawtab_type),target,intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 real(dp),parameter :: one_over_fourpi   = one/sqrt(four_pi)
 real(dp),parameter :: sqr_fourpi_over_3 = sqrt(four_pi/3)
 real(dp),parameter :: QuarterFineStruct2=(half/InvFineStruct)**2
 integer :: iatom,iatom_tot,itypat,ierr,ipts,ignt,sgnkappa
 integer :: idum,option,usenhat,usekden,usecore,xclevel,nkxc,my_comm_atom
 integer :: mesh_size,mesh_size_cor,lmn_size,lmn2_size,lmn_size_j,lmn_size_cor
 integer :: lm_size,ln_size,ln_size_j,ln_size_cor,nspinor_cor
 integer :: ilmn,ilm,iln,jl,jm,jm_re,jm_im,jlmn,jlm,jlm_re,jlm_im,jln,js,klm_re,klm_im
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: cgc,compch_sph_dum,eexc_dum,eexcdc_dum,jmj
 real(dp) :: fact_re,fact_im,gx_re,gx_im,gy_re,gy_im,gz_re,gz_im,if3
 character(len=500) :: msg
!arrays
 integer,pointer :: my_atmtab(:),indlmn(:,:),indlmn_j(:,:)
 logical,allocatable :: lmselect(:)
 real(dp) :: nhat_dum(1,1,1),trho_dum(1,1,1),kxc_dum(1,1,1),k3xc_dum(1,1,1)
 real(dp),allocatable :: rho1(:,:,:),tau1(:,:,:),rhosph(:),vhartree(:),vxc(:,:,:)
 real(dp),allocatable :: intf3(:,:),potsph(:),dVdr(:),func(:)
 real(dp),pointer :: phi_j(:,:),soc_ij(:,:,:,:,:)
 type(pawrad_type),pointer :: pawrd
 type(pawtab_type),pointer :: pawtb

! ************************************************************************

!Some checks in case of core-valence (option_core=1)
 if (option_core/=0.and.option_core/=1) then
   msg='Wrong option_core value!'
   ABI_BUG(msg)
 end if
 if (option_core==1) then
!  Check if we have the optional arguments
   if ((.not.present(phi_cor)).or.(.not.present(indlmn_cor))) then
     msg='For core-valence calculation, need phi_cor and indlmn_cor!'
     ABI_BUG(msg)
   end if
!  Check if we have relativistic core wave functions
   if (size(indlmn_cor,1)<8) then
     write(msg,'(a)') 'Wrong 1st dim. of indlmn_cor in pawnabla_soc_init (need spinors)!'
     ABI_BUG(msg)
   end if
 endif

!Some useful variables
 usekden=pawxc_get_usekden(ixc)
 usecore=1 ; nkxc=0 ; usenhat=0
 xclevel=pawxc_get_xclevel(ixc)
 if (option_core==1) then
   mesh_size_cor=size(phi_cor,1)
   lmn_size_cor=size(indlmn_cor,2) !Includes spinors
   ln_size_cor=size(phi_cor,2)
   nspinor_cor=maxval(indlmn_cor(6,1:lmn_size_cor))
 end if

!Prepare output arrays
 if (allocated(phisocphj)) then
   do iatom=1,natom
      if (allocated(phisocphj(iatom)%value)) then
        ABI_FREE(phisocphj(iatom)%value)
      end if
   end do
   ABI_FREE(phisocphj)
 end if
 ABI_MALLOC(phisocphj,(natom))
 do iatom=1,natom
   lmn_size=pawtab(typat(iatom))%lmn_size
   if (option_core==0) then
     ABI_MALLOC(phisocphj(iatom)%value,(2,2,3,lmn_size,lmn_size))
   else
     !lmn_size_cor is double because it contains the spin component
     ABI_MALLOC(phisocphj(iatom)%value,(2,2,3,lmn_size,lmn_size_cor))
   end if
   phisocphj(iatom)%value=zero
 end do

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,&
&                   my_natom_ref=my_natom)

!-------------------------------------------------------------------
!Loop over atoms

 do iatom=1,my_natom

!  Atom-dependent data
   itypat=typat(iatom)
   pawrd => pawrad(itypat)
   pawtb  => pawtab(itypat)
   lmn_size=pawtb%lmn_size
   lmn2_size=pawtb%lmn2_size
   lm_size=pawtb%lcut_size**2
   ln_size=pawtb%basis_size
   mesh_size=pawtb%mesh_size
   indlmn => pawtb%indlmn
   ABI_MALLOC(lmselect,(lm_size))
   lmselect(:)=.true.

!  Distinguish valence-valence and core-valence cases
   if (option_core==0) then
     ln_size_j=ln_size
     lmn_size_j=lmn_size
     indlmn_j => pawtb%indlmn
     phi_j => pawtb%phi
   else
     ln_size_j=ln_size_cor
     lmn_size_j=lmn_size_cor
     indlmn_j => indlmn_cor
     phi_j => phi_cor
     if (mesh_size_cor<mesh_size) then
       msg='mesh_size and mesh_size_cor not compatible!'
       ABI_BUG(msg)
     end if
   endif

!  Manage parallelism
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
   soc_ij => phisocphj(iatom_tot)%value(:,:,:,:,:)

!-------------------------------------------------------------------
!Compute all-electron density

   ABI_MALLOC(rho1,(mesh_size,lm_size,nspden))
   option=2 ; rho1=zero
   call pawdensities(compch_sph_dum,1,iatom_tot,lmselect,lmselect,lm_size,&
&       nhat_dum,nspden,-1,0,option,-1,0,pawang,0,pawrd,&
&       pawrhoij(iatom),pawtb,rho1,trho_dum)
   if (usekden==1) then
     ABI_MALLOC(tau1,(mesh_size,lm_size,nspden))
     tau1=zero
     call pawkindensities(1,lmselect,lm_size,nspden,-1,option,-1,&
&         pawang,pawrd,pawrhoij(iatom),pawtb,tau1,trho_dum)
   end if

!-------------------------------------------------------------------
!Compute spherical potential and compute its first derivative dV/dr

   ABI_MALLOC(potsph,(mesh_size))
   ABI_MALLOC(dVdr,(mesh_size))
   potsph=zero ; dVdr=zero

!  Compute XC potential
   option=1
   if (pawxcdev/=0) then
     ABI_MALLOC(vxc,(mesh_size,lm_size,nspden))
     vxc=zero
     call pawxcm(pawtb%coredens,eexc_dum,eexcdc_dum,idum,ixc,kxc_dum,lm_size,&
&         lmselect,nhat_dum,nkxc,.false.,mesh_size,nspden,option,pawang,pawrd,&
&         pawxcdev,rho1,usecore,usenhat,vxc,xclevel,xc_denpos)
     potsph(1:mesh_size)=half*(vxc(1:mesh_size,1,1)+vxc(1:mesh_size,1,nspden))
   else
     ABI_MALLOC(vxc,(mesh_size,pawang%angl_size,nspden))
     vxc=zero
     call pawxc(pawtb%coredens,eexc_dum,eexcdc_dum,ixc,kxc_dum,k3xc_dum,&
&         lm_size,lmselect,nhat_dum,nkxc,nkxc,.false.,mesh_size,nspden,option,pawang,&
&         pawrd,rho1,usecore,usenhat,vxc,xclevel,xc_denpos,&
&         coretau=pawtb%coretau,taur=tau1)
     potsph(1:mesh_size)=zero
     do ipts=1,pawang%angl_size
       potsph(1:mesh_size)=potsph(1:mesh_size) &
&        +half*(vxc(1:mesh_size,ipts,1)+vxc(1:mesh_size,ipts,nspden))*pawang%angwgth(ipts)
     end do
     potsph(1:mesh_size)=sqrt(four_pi)*potsph(1:mesh_size)
   end if

!  Compute Hartree potential
   ABI_MALLOC(vhartree,(mesh_size))
   ABI_MALLOC(rhosph,(mesh_size))
   vhartree=zero ; rhosph=zero
   rhosph(1:mesh_size)=rho1(1:mesh_size,1,1)
   if (usecore==1) rhosph(1:mesh_size)=rhosph(1:mesh_size)+sqrt(four_pi)*pawtb%coredens(1:mesh_size)
   rhosph(1:mesh_size)=rhosph(1:mesh_size)*four_pi*pawrd%rad(1:mesh_size)**2
   call poisson(rhosph,0,pawrd,vhartree)
   vhartree(2:mesh_size)=(vhartree(2:mesh_size)-sqrt(four_pi)*znucl(itypat))/pawrd%rad(2:mesh_size)
   call pawrad_deducer0(vhartree,mesh_size,pawrd)
   potsph(1:mesh_size)=potsph(1:mesh_size)+vhartree(1:mesh_size)

!  Apply angular scaling factor
   potsph(1:mesh_size)=one_over_fourpi*potsph(1:mesh_size)

!  Compute 1st derivative of potential
   call nderiv_gen(dVdr,potsph,pawrd)

!  Multiply by relativistic factor
   dVdr(1:mesh_size)=dVdr(1:mesh_size)*(one/(one-potsph(1:mesh_size)/InvFineStruct**2))

   ABI_FREE(vxc)
   ABI_FREE(vhartree)
   ABI_FREE(potsph)
   ABI_FREE(rhosph)
   ABI_FREE(lmselect)
   ABI_FREE(rho1)
   if (usekden==1) then
     ABI_FREE(tau1)
   end if

!-------------------------------------------------------------------
!Compute radial and angular contributions

   ABI_MALLOC(intf3,(ln_size,ln_size_j))
   intf3=zero

!  >>>> Calculate f_3= alpha^2/4 int[dr ui*(r) uj(r) dV/dr]
   ABI_MALLOC(func,(mesh_size))
   do jln=1,ln_size_j
     do iln=1,ln_size
       func(1:mesh_size)=dVdr(1:mesh_size)*pawtb%phi(1:mesh_size,iln)*phi_j(1:mesh_size,jln)
       call simp_gen(intf3(iln,jln),func,pawrd)
     end do
   end do
   intf3(:,:)=QuarterFineStruct2*spnorbscl*intf3(:,:)
   ABI_FREE(func)

!  Loop over initial states (valence or core, according to option_core)
   do jlmn=1,lmn_size_j
     jl=indlmn_j(1,jlmn)
     jm=indlmn_j(2,jlmn)
     jlm=indlmn_j(4,jlmn)
     jln=indlmn_j(5,jlmn)

!    In case of spinorial core wave function, we have to handle imaginary
!      spherical harmonics as linear combination of real spherical harmonics
!    See Brouwer et al, CPC 266, 108029 (2021), equation (21)
     if (option_core==1) then
       jm_re= abs(jm)
       jm_im=-abs(jm)
       jlm_re=jl*(jl+1)+jm_re+1
       jlm_im=jl*(jl+1)+jm_im+1
!      Calculate spinor dependent coefficients
       sgnkappa=indlmn_j(3,jlmn)   !sign of kappa
       jmj=half*indlmn_j(8,jlmn)   !2mj is stored in indlmn_cor
       js=indlmn_cor(6,jlmn)       !1 is up, 2 is down
       if (sgnkappa==1) then
         if(js==1) then
           cgc= sqrt((dble(jl)-jmj+half)/dble(2*jl+1))
         else
           cgc=-sqrt((dble(jl)+jmj+half)/dble(2*jl+1))
         endif
       else
         if(js==1) then
           cgc= sqrt((dble(jl)+jmj+half)/dble(2*jl+1))
         else
           cgc= sqrt((dble(jl)-jmj+half)/dble(2*jl+1))
         endif
       endif

!      Calculate factors to convert from complex to real sph. harm.
       if (jm<0) then
         fact_re=sqr_fourpi_over_3*half_sqrt2*cgc
         fact_im=-fact_re
       else if (jm>0) then
         fact_re=sqr_fourpi_over_3*half_sqrt2*cgc*(-1)**jm
         fact_im=fact_re
       else
         fact_re=sqr_fourpi_over_3
         fact_im=0
       end if
     else ! valence-valence case (real)
       js=1
       jlm_re=jlm ; jlm_im=jlm
       fact_re=sqr_fourpi_over_3
       fact_im=0
     end if

!    Loop over final states
     do ilmn=1,lmn_size
       ilm=indlmn(4,ilmn)
       iln=indlmn(5,ilmn)

!      >>>> Calculate g_ij=(g_x,g_y,g_z) = sqrt(4pi/3) int dOmega Ylm Ylm' S1-1,0,1
!              using real Gaunt coefficients
       gx_re=zero;gy_re=zero;gz_re=zero         
       gx_im=zero;gy_im=zero;gz_im=zero
       if3=zero

!      jl was set as a flag for invalid combinations
!        i.e. m=-(l+1) or m=(l+1)
!      In these cases, cgc=0 ; so gx=gy=gz=0 
       if (jl/=-1) then
         if3=intf3(iln,jln)
         klm_re=merge((jlm_re*(jlm_re-1))/2+ilm,(ilm*(ilm-1))/2+jlm_re,ilm<=jlm_re)

!        Real parts
!        M=-1
         ignt=pawang%gntselect(2,klm_re) !get index for L=1 M =-1 ilm jlm_re
         if (ignt/=0) gy_re=fact_re*pawang%realgnt(ignt)
!        M=0
         ignt=pawang%gntselect(3,klm_re) !get index for L=1 M = 0 ilm jlm_re
         if (ignt/=0) gz_re=fact_re*pawang%realgnt(ignt)
!        M=1
         ignt=pawang%gntselect(4,klm_re) !get index for L=1 M = 1 ilm jlm_re
         if (ignt/=0) gx_re=fact_re*pawang%realgnt(ignt)

!        Imaginary parts
         if (option_core==1) then
           klm_im=merge((jlm_im*(jlm_im-1))/2+ilm,(ilm*(ilm-1))/2+jlm_im,ilm<=jlm_im)
!          M=-1
           ignt=pawang%gntselect(2,klm_im) !get index for L=1 M =-1 ilm jlm_im
           if (ignt/=0) gy_im=fact_im*pawang%realgnt(ignt)
!          M=0
           ignt=pawang%gntselect(3,klm_im) !get index for L=1 M = 0 ilm jlm_im
           if (ignt/=0) gz_im=fact_im*pawang%realgnt(ignt)
!          M=1
           ignt=pawang%gntselect(4,klm_im) !get index for L=1 M = 1 ilm jlm_im
           if (ignt/=0) gx_im=fact_im*pawang%realgnt(ignt)
         end if
       end if

!      >>>> Calculate Sigma X g_ij

       if (option_core==0.or.js==1) then
         !(Sigma^up-up X gij)_x = -gy*f_3
         soc_ij(1,1,1,ilmn,jlmn)=-if3*gy_re           ! real part
         soc_ij(2,1,1,ilmn,jlmn)=-if3*gy_im           ! imag part
         !(Sigma^up-up X gij)_y = gx*f_3
         soc_ij(1,1,2,ilmn,jlmn)= if3*gx_re           ! real part
         soc_ij(2,1,2,ilmn,jlmn)= if3*gx_im           ! imag part
         !(Sigma^up-up X gij)_z = 0
         soc_ij(1,1,3,ilmn,jlmn)= zero                ! real part
         soc_ij(2,1,3,ilmn,jlmn)= zero                ! imag part

         !(Sigma^dn-up X gij)_x = i.gz*f_3
         soc_ij(1,2,1,ilmn,jlmn)=-if3*gz_im           ! real part
         soc_ij(2,2,1,ilmn,jlmn)= if3*gz_re           ! imag part
         !(Sigma^dn-up X gij)_y = -gz*f_3
         soc_ij(1,2,2,ilmn,jlmn)=-if3*gz_re           ! real part
         soc_ij(2,2,2,ilmn,jlmn)=-if3*gz_im           ! imag part
         !(Sigma^dn-up X gij)_z = (gy-i.gx)*f_3
         soc_ij(1,2,3,ilmn,jlmn)= if3*(gy_re+gx_im)   ! real part
         soc_ij(2,2,3,ilmn,jlmn)= if3*(gy_im-gx_re)   ! imag part

       else if (option_core==1.and.js==2) then
         !(Sigma^up-dn X gij^dn)_x = -i.gz^dn*f_3
         soc_ij(1,1,1,ilmn,jlmn)= if3*gz_im           ! real part
         soc_ij(2,1,1,ilmn,jlmn)=-if3*gz_re           ! imag part
         !(Sigma^up-dn X gij^dn)_y = -gz^dn*f_3
         soc_ij(1,1,2,ilmn,jlmn)=-if3*gz_re           ! real part
         soc_ij(2,1,2,ilmn,jlmn)=-if3*gz_im           ! imag part
         !(Sigma^up-dn X gij^dn)_z = (gy^dn+i.gx^dn)*f_3
         soc_ij(1,1,3,ilmn,jlmn)= if3*(gy_re-gx_im)   ! real part
         soc_ij(2,1,3,ilmn,jlmn)= if3*(gy_im+gx_re)   ! imag part

         !(Sigma^dn-dn X gij^dn)_x =  gy^dn*f_3
         soc_ij(1,2,1,ilmn,jlmn)= if3*gy_re           ! real part
         soc_ij(2,2,1,ilmn,jlmn)= if3*gy_im           ! imag part
         !(Sigma^dn-dn X gij^dn)_y = -gx^dn*f_3
         soc_ij(1,2,2,ilmn,jlmn)=-if3*gx_re           ! real part
         soc_ij(2,2,2,ilmn,jlmn)=-if3*gx_im           ! imag part
         !(Sigma^dn-dn X gij^dn)_z = 0
         soc_ij(1,2,3,ilmn,jlmn)= zero                ! real part
         soc_ij(2,2,3,ilmn,jlmn)= zero                ! imag part
       end if

     end do ! ilmn
   end do ! jlmn
     
   ABI_FREE(dVdr)
   ABI_FREE(intf3)

 end do  ! iatom
 
!Reduction in case of parallelism
 if (paral_atom) then
   call xmpi_sum(phisocphj,my_comm_atom,ierr)
 end if
 
!Destroy atom table(s) used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 end subroutine pawnabla_soc_init
!!***

!----------------------------------------------------------------------

END MODULE m_paw_optics
!!***
