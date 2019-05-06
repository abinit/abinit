!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paw_optics
!! NAME
!!  m_paw_optics
!!
!! FUNCTION
!!  This module contains several routines related to conductivity:
!!    optical conductivity, X spectroscopy, linear susceptibility, ...
!!
!! COPYRIGHT
!! Copyright (C) 2018-2019 ABINIT group (SM,VR,FJ,MT,PGhosh)
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
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_wffile
 use m_abicore
 use m_hdr

 use m_time,         only : timab
 use m_io_tools,     only : open_file,get_unit
 use m_pawpsp,       only : pawpsp_read_corewf
 use m_pawrad,       only : pawrad_type, pawrad_deducer0, simp_gen
 use m_pawtab,       only : pawtab_type
 use m_pawcprj,      only : pawcprj_type, pawcprj_alloc, pawcprj_get, &
&                           pawcprj_free,pawcprj_mpi_allgather
 use m_paw_onsite,   only : pawnabla_init,pawnabla_core_init
 use m_mpinfo,       only : destroy_mpi_enreg,nullify_mpi_enreg,initmpi_seq,proc_distrb_cycle
 use m_numeric_tools,only : kramerskronig
 use m_geometry,     only : metric
 use m_hide_lapack,  only : matrginv

 implicit none

 private

!public procedures.
 public :: optics_paw
 public :: optics_paw_core
 public :: linear_optics_paw

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
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data:
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  (only writing in a file)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      hdr_io,setnabla_ylm,pawcprj_alloc,pawcprj_free,pawcprj_get
!!      pawcprj_mpi_allgather,pawrad_deducer0,simp_gen,timab,wffclose,wffopen
!!      xmpi_exch,xmpi_sum,xmpi_sum_master
!!
!! SOURCE

 subroutine optics_paw(atindx1,cg,cprj,dimcprj,dtfil,dtset,eigen0,gprimd,hdr,kg,&
&               mband,mcg,mcprj,mkmem,mpi_enreg,mpsang,mpw,natom,nkpt,npwarr,nsppol,&
&               pawrad,pawtab)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcg,mcprj,mkmem,mpsang,mpw,natom,nkpt,nsppol
 type(MPI_type),intent(in) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
!arrays
 integer,intent(in) :: atindx1(natom),dimcprj(natom)
 integer,intent(in) :: kg(3,mpw*mkmem),npwarr(nkpt)
 real(dp),intent(in) :: eigen0(mband*nkpt*nsppol)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(inout) :: cg(2,mcg)
 type(pawcprj_type),target,intent(inout) :: cprj(natom,mcprj)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
 type(pawtab_type),target,intent(inout) :: pawtab(dtset%ntypat)

!Local variables-------------------------------
!scalars
 integer :: iomode,bdtot_index,cplex,etiq,fformopt,iatom,ib,ibg,ibsp
 integer :: icg,ierr,ikg,ikpt,ilmn,ount
 integer :: iorder_cprj,ipw,ispinor,isppol,istwf_k,itypat,iwavef
 integer :: jb,jbsp,jlmn,jwavef,lmn_size,mband_cprj
 integer :: me,me_kpt,my_nspinor,nband_k,nband_cprj_k,npw_k,sender
 integer :: spaceComm_band,spaceComm_bandfftspin,spaceComm_fft,spaceComm_k,spaceComm_spin,spaceComm_w
 logical :: already_has_nabla,cprj_paral_band,mykpt
 real(dp) :: cgnm1,cgnm2,cpnm1,cpnm2
 character(len=500) :: message
!arrays
 integer :: tmp_shape(3)
 integer,allocatable :: kg_k(:,:)
 real(dp) :: kpoint(3),tsec(2)
 real(dp),allocatable :: kpg_k(:,:),psinablapsi(:,:,:,:),tnm(:,:,:,:)
 type(pawcprj_type),pointer :: cprj_k(:,:),cprj_k_loc(:,:)
 type(wffile_type) :: wff1

! ************************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 ABI_CHECK(mkmem/=0,"mkmem==0 not supported anymore!")

!----------------------------------------------------------------------------------
!1- Computation of on-site contribution: <phi_i|nabla|phi_j>-<tphi_i|nabla|tphi_j>
!----------------------------------------------------------------------------------

 already_has_nabla=all(pawtab(:)%has_nabla==2)
 call pawnabla_init(mpsang,dtset%ntypat,pawrad,pawtab)

!----------------------------------------------------------------------------------
!2- Computation of <psi_n|-i.nabla|psi_m> for each k
!----------------------------------------------------------------------------------

!Init parallelism
 spaceComm_w=mpi_enreg%comm_cell
 me=xmpi_comm_rank(spaceComm_w)
 if (mpi_enreg%paral_kgb==1) then
   spaceComm_k=mpi_enreg%comm_kpt
   spaceComm_fft=mpi_enreg%comm_fft
   spaceComm_band=mpi_enreg%comm_band
   spaceComm_spin=mpi_enreg%comm_spinor
   spaceComm_bandfftspin=mpi_enreg%comm_bandspinorfft
   me_kpt=mpi_enreg%me_kpt
 else
   spaceComm_k=spaceComm_w
   spaceComm_band=0;spaceComm_fft=0;spaceComm_spin=xmpi_comm_self;spaceComm_bandfftspin=0
   me_kpt=me
 end if
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

!Determine if cprj datastructure is distributed over bands
 mband_cprj=mcprj/(my_nspinor*mkmem*nsppol)
 cprj_paral_band=(mband_cprj<mband)

!PAW file
 iorder_cprj=0

!Initialize main variables
 ABI_ALLOCATE(psinablapsi,(2,3,mband,mband))
 psinablapsi=zero

!open _OPT file for proc 0
 iomode=IO_MODE_FORTRAN_MASTER
 fformopt=610
 ount = get_unit()
 call WffOpen(iomode,spaceComm_k,dtfil%fnameabo_app_opt,ierr,wff1,0,me,ount)
 call hdr_io(fformopt,hdr,2,wff1)
 if (me==0) then
   write(ount)(eigen0(ib),ib=1,mband*nkpt*nsppol)
 end if

!LOOP OVER SPINS
 ibg=0;icg=0
 bdtot_index=0
 do isppol=1,nsppol

!  LOOP OVER k POINTS
   ikg=0
   do ikpt=1,nkpt

     nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
     etiq=ikpt+(isppol-1)*nkpt

!    Select k-points for current proc
     mykpt=.true.
     mykpt=.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_kpt))
     if (mykpt) then

!      Allocations depending on k-point
       kpoint(:)=dtset%kptns(:,ikpt)
       istwf_k=dtset%istwfk(ikpt)
       npw_k=npwarr(ikpt)
       cplex=2;if (istwf_k>1) cplex=1
       ABI_ALLOCATE(kg_k,(3,npw_k))
       ABI_ALLOCATE(kpg_k,(npw_k*dtset%nspinor,3))

!      Extract cprj for this k-point according to mkmem
       nband_cprj_k=nband_k;if (cprj_paral_band) nband_cprj_k=nband_k/mpi_enreg%nproc_band
       if (mkmem*nsppol/=1) then
         ABI_DATATYPE_ALLOCATE(cprj_k_loc,(natom,my_nspinor*nband_cprj_k))
         call pawcprj_alloc(cprj_k_loc,0,dimcprj)
         call pawcprj_get(atindx1,cprj_k_loc,cprj,natom,1,ibg,ikpt,iorder_cprj,isppol,&
&         mband_cprj,mkmem,natom,nband_cprj_k,nband_cprj_k,my_nspinor,nsppol,dtfil%unpaw,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
       else
         cprj_k_loc => cprj
       end if

!      if cprj are distributed over bands, gather them (because we need to mix bands)
       if (cprj_paral_band) then
         ABI_DATATYPE_ALLOCATE(cprj_k,(natom,my_nspinor*nband_k))
         call pawcprj_alloc(cprj_k,0,dimcprj)
         call pawcprj_mpi_allgather(cprj_k_loc,cprj_k,natom,my_nspinor*nband_cprj_k,dimcprj,0,&
&         mpi_enreg%nproc_band,mpi_enreg%comm_band,ierr,rank_ordered=.false.)
       else
         cprj_k => cprj_k_loc
       end if

!      Get G-vectors for this k-point.
       kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
       ikg=ikg+npw_k

!      Calculation of k+G in cartesian coordinates
       do ipw=1,npw_k
         kpg_k(ipw,1)=(kpoint(1)+kg_k(1,ipw))*gprimd(1,1)&
&         +(kpoint(2)+kg_k(2,ipw))*gprimd(1,2)&
&         +(kpoint(3)+kg_k(3,ipw))*gprimd(1,3)
         kpg_k(ipw,2)=(kpoint(1)+kg_k(1,ipw))*gprimd(2,1)&
&         +(kpoint(2)+kg_k(2,ipw))*gprimd(2,2)&
&         +(kpoint(3)+kg_k(3,ipw))*gprimd(2,3)
         kpg_k(ipw,3)=(kpoint(1)+kg_k(1,ipw))*gprimd(3,1)&
&         +(kpoint(2)+kg_k(2,ipw))*gprimd(3,2)&
&         +(kpoint(3)+kg_k(3,ipw))*gprimd(3,3)
       end do !ipw
       kpg_k=two_pi*kpg_k
       if (dtset%nspinor==2) kpg_k(npw_k+1:2*npw_k,1:3)=kpg_k(1:npw_k,1:3)
       ABI_DEALLOCATE(kg_k)

!      2-A Computation of <psi_tild_n|-i.nabla|psi_tild_m>
!      ----------------------------------------------------------------------------------
!      Computation of (C_nk^*)*C_mk*(k+g) in cartesian coordinates

       ABI_ALLOCATE(tnm,(2,3,nband_k,nband_k))
       tnm=zero

!      Loops on bands
       do jb=1,nband_k
         jwavef=(jb-1)*npw_k*my_nspinor+icg
         if (xmpi_paral==1.and.mpi_enreg%paral_kgb/=1) then
           tmp_shape = shape(mpi_enreg%proc_distrb)
           if (ikpt > tmp_shape(1)) then
             message = ' optics_paw : ikpt out of bounds '
             MSG_ERROR(message)
           end if
!          MJV 6/12/2008: looks like mpi_enreg may not be completely initialized here
           if (abs(mpi_enreg%proc_distrb(ikpt,jb,isppol)-me_kpt)/=0) cycle
         end if
         do ib=1,jb
           iwavef=(ib-1)*npw_k*my_nspinor+icg

!          Computation of (C_nk^*)*C_mk*(k+g) in cartesian coordinates
           if (cplex==1) then
             do ipw=1,npw_k*my_nspinor
               cgnm1=cg(1,ipw+iwavef)*cg(1,ipw+jwavef)
               tnm(1,1:3,ib,jb)=tnm(1,1:3,ib,jb)+cgnm1*kpg_k(ipw,1:3)
             end do
           else
             do ipw=1,npw_k*my_nspinor
               cgnm1=cg(1,ipw+iwavef)*cg(1,ipw+jwavef)+cg(2,ipw+iwavef)*cg(2,ipw+jwavef)
               cgnm2=cg(1,ipw+iwavef)*cg(2,ipw+jwavef)-cg(2,ipw+iwavef)*cg(1,ipw+jwavef)
               tnm(1,1:3,ib,jb)=tnm(1,1:3,ib,jb)+cgnm1*kpg_k(ipw,1:3)
               tnm(2,1:3,ib,jb)=tnm(2,1:3,ib,jb)+cgnm2*kpg_k(ipw,1:3)
             end do
           end if

!          Second half of the (n,m) matrix
           if (ib/=jb) then
             tnm(1,1:3,jb,ib)= tnm(1,1:3,ib,jb)
             tnm(2,1:3,jb,ib)=-tnm(2,1:3,ib,jb)
           end if

         end do ! ib
       end do ! jb

!      ATTEMPT: USING BLAS (not successful)
!      allocate(dg(2,npw_k*my_nspinor,nband_k),tnm_tmp(2,nband_k,nband_k))
!      do ii=1,3
!      do ib=1,nband_k
!      iwavef=icg+(ib-1)*npw_k*my_nspinor
!      do ipw=1,my_nspinor*npw_k
!      dg(:,ipw,nband_k)=cg(:,ipw+iwavef)*kpg_k(ipw,ii)
!      end do
!      end do
!      call zgemm('C','N',nband_k,nband_k,npw_k*my_nspinor,one,dg,npw_k*my_nspinor,&
!      &                cg(:,1+icg:npw_k*my_nspinor*nband_k+icg),npw_k*my_nspinor,zero,tnm_tmp,nband_k)
!      tnm(:,ii,:,:)=tnm_tmp(:,:,:)
!      deallocate(dg,tnm_tmp)
!      end do

!      Reduction in case of parallelism
       if (mpi_enreg%paral_kgb == 1) then
         call timab(48,1,tsec)
         call xmpi_sum_master(tnm,0,spaceComm_bandfftspin,ierr)
         call timab(48,2,tsec)
       end if

       psinablapsi(:,:,:,:)=tnm(:,:,:,:)

!      2-B Computation of <psi_n|p_i><p_j|psi_m>(<phi_i|-i.nabla|phi_j>-<tphi_i|-i.nabla|tphi_j>)
!      ----------------------------------------------------------------------------------

       tnm=zero

!      Loops on bands
       do jb=1,nband_k

         if (mpi_enreg%paral_kgb==1) then
           if (mod(jb-1,mpi_enreg%nproc_band)/=mpi_enreg%me_band) cycle
           if (abs(mpi_enreg%proc_distrb(ikpt,jb,isppol)-me_kpt)/=0) cycle
         end if
         do ib=1,jb
           ibsp=(ib-1)*my_nspinor;jbsp=(jb-1)*my_nspinor

           if (cplex==1) then
             do ispinor=1,my_nspinor
               ibsp=ibsp+1;jbsp=jbsp+1
               do iatom=1,natom
                 itypat=dtset%typat(iatom)
                 lmn_size=pawtab(itypat)%lmn_size
                 do jlmn=1,lmn_size
                   do ilmn=1,lmn_size
                     cpnm1=cprj_k(iatom,ibsp)%cp(1,ilmn)*cprj_k(iatom,jbsp)%cp(1,jlmn)
                     tnm(2,:,ib,jb)=tnm(2,:,ib,jb)+cpnm1*pawtab(itypat)%nabla_ij(:,ilmn,jlmn)
                   end do !ilmn
                 end do !jlmn
               end do !iatom
             end do !ispinor
           else
             do ispinor=1,my_nspinor
               ibsp=ibsp+1;jbsp=jbsp+1
               do iatom=1,natom
                 itypat=dtset%typat(iatom)
                 lmn_size=pawtab(itypat)%lmn_size
                 do jlmn=1,lmn_size
                   do ilmn=1,lmn_size
                     cpnm1=(cprj_k(iatom,ibsp)%cp(1,ilmn)*cprj_k(iatom,jbsp)%cp(1,jlmn) &
&                     +cprj_k(iatom,ibsp)%cp(2,ilmn)*cprj_k(iatom,jbsp)%cp(2,jlmn))
                     cpnm2=(cprj_k(iatom,ibsp)%cp(1,ilmn)*cprj_k(iatom,jbsp)%cp(2,jlmn) &
&                     -cprj_k(iatom,ibsp)%cp(2,ilmn)*cprj_k(iatom,jbsp)%cp(1,jlmn))
                     tnm(1,:,ib,jb)=tnm(1,:,ib,jb)+cpnm2*pawtab(itypat)%nabla_ij(:,ilmn,jlmn)
                     tnm(2,:,ib,jb)=tnm(2,:,ib,jb)-cpnm1*pawtab(itypat)%nabla_ij(:,ilmn,jlmn)
                   end do !ilmn
                 end do !jlmn
               end do !iatom
             end do !ispinor
           end if

!          Second half of the (n,m) matrix
           if (ib/=jb) then
             tnm(1,1:3,jb,ib)=-tnm(1,1:3,ib,jb)
             tnm(2,1:3,jb,ib)= tnm(2,1:3,ib,jb)
           end if

!          End loops on bands
         end do ! jb
       end do !ib

!      Reduction in case of parallelism
       if (mpi_enreg%paral_kgb==1) then
         call timab(48,1,tsec)
         call xmpi_sum(tnm,mpi_enreg%comm_bandspinor,ierr)
         call timab(48,2,tsec)
       else
!        In that case parallelisation on nspinor not implemented yet;put this line just in case
         call xmpi_sum(tnm,spacecomm_spin,ierr)
       end if

       psinablapsi(:,:,:,:)=psinablapsi(:,:,:,:)+tnm(:,:,:,:)
       ABI_DEALLOCATE(tnm)

       if (mkmem/=0) then
         ibg = ibg +       my_nspinor*nband_cprj_k
         icg = icg + npw_k*my_nspinor*nband_k
       end if

       if (cprj_paral_band) then
         call pawcprj_free(cprj_k)
         ABI_DATATYPE_DEALLOCATE(cprj_k)
       end if
       if (mkmem*nsppol/=1) then
         call pawcprj_free(cprj_k_loc)
         ABI_DATATYPE_DEALLOCATE(cprj_k_loc)
       end if
       ABI_DEALLOCATE(kpg_k)

       if (me==0) then
         write(ount)((psinablapsi(1:2,1,ib,jb),ib=1,nband_k),jb=1,nband_k)
         write(ount)((psinablapsi(1:2,2,ib,jb),ib=1,nband_k),jb=1,nband_k)
         write(ount)((psinablapsi(1:2,3,ib,jb),ib=1,nband_k),jb=1,nband_k)
       elseif (mpi_enreg%me_band==0.and.mpi_enreg%me_fft==0) then
         call xmpi_exch(psinablapsi,etiq,me_kpt,psinablapsi,0,spaceComm_k,ierr)
       end if

     elseif (me==0) then
       sender=minval(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol))
       call xmpi_exch(psinablapsi,etiq,sender,psinablapsi,0,spaceComm_k,ierr)
       write(ount)((psinablapsi(1:2,1,ib,jb),ib=1,nband_k),jb=1,nband_k)
       write(ount)((psinablapsi(1:2,2,ib,jb),ib=1,nband_k),jb=1,nband_k)
       write(ount)((psinablapsi(1:2,3,ib,jb),ib=1,nband_k),jb=1,nband_k)
     end if ! mykpt

     bdtot_index=bdtot_index+nband_k
!    End loop on spin,kpt
   end do ! ikpt
 end do !isppol

!Close file
 call WffClose(wff1,ierr)

!Datastructures deallocations
 ABI_DEALLOCATE(psinablapsi)
 if (.not.already_has_nabla) then
   do itypat=1,dtset%ntypat
     if (allocated(pawtab(itypat)%nabla_ij)) then
       ABI_DEALLOCATE(pawtab(itypat)%nabla_ij)
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
!! Copyright (C) 2005-2019 ABINIT group (SM,MT)
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
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data:
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  (only writing in a file)
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      hdr_io,setnabla_ylm,pawcprj_alloc,pawcprj_free,pawcprj_get
!!      pawcprj_mpi_allgather,pawpsp_read_corewf,pawrad_deducer0,simp_gen,timab
!!      wffclose,wffopen,xmpi_exch,xmpi_sum_master
!!
!! SOURCE

 subroutine optics_paw_core(atindx1,cprj,dimcprj,dtfil,dtset,eigen0,filpsp,hdr,&
&               mband,mcprj,mkmem,mpi_enreg,mpsang,natom,nkpt,nsppol,pawrad,pawtab)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcprj,mkmem,mpsang,natom,nkpt,nsppol
 type(MPI_type),intent(in) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
!arrays
 integer,intent(in) :: atindx1(natom),dimcprj(natom)
 character(len=fnlen),intent(in) :: filpsp(dtset%ntypat)
 real(dp),intent(in) :: eigen0(mband*nkpt*nsppol)
 type(pawcprj_type),target,intent(inout) :: cprj(natom,mcprj)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
 type(pawtab_type),target,intent(inout) :: pawtab(dtset%ntypat)

!Local variables-------------------------------
!scalars
 integer :: bdtot_index,cplex,etiq,iatom,ib,ibg
 integer :: ierr,ikpt,ilmn,iln,ount
 integer :: iorder_cprj,ispinor,isppol,istwf_k,itypat
 integer :: jb,jbsp,jlmn,lmn_size,lmncmax,mband_cprj
 integer :: me,me_kpt,my_nspinor,nband_cprj_k,nband_k,nphicor
 integer :: sender,spaceComm_bandspin,spaceComm_k,spaceComm_w
 integer :: iomode,fformopt
 logical :: already_has_nabla,cprj_paral_band,ex,mykpt
 character(len=fnlen) :: filecore
 real(dp) :: cpnm1,cpnm2
!arrays
 integer,allocatable :: indlmn_core(:,:),lcor(:),ncor(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: energy_cor(:),phi_cor(:,:),psinablapsi(:,:,:,:,:),tnm(:,:,:,:,:)
 type(pawcprj_type),pointer :: cprj_k(:,:),cprj_k_loc(:,:)
 type(wffile_type) :: wff1

! ************************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 ABI_CHECK(mkmem/=0,"mkmem==0 not supported anymore!")

!------------------------------------------------------------------------------------------------
!0-Reading core wavefunctions
!------------------------------------------------------------------------------------------------

!Note: core WF is read for itypat=1
 filecore=trim(filpsp(1))//'.corewf'
 inquire(file=filecore,exist=ex)
 if (ex) then
   !Use <filepsp>.corewf
   call pawpsp_read_corewf(energy_cor,indlmn_core,lcor,lmncmax,ncor,nphicor,pawrad(1),phi_cor,&
&   filename=filecore)
 else
   !Use default name
   call pawpsp_read_corewf(energy_cor,indlmn_core,lcor,lmncmax,ncor,nphicor,pawrad(1),phi_cor)
 end if

!----------------------------------------------------------------------------------
!1-Computation of phipphj=<phi_i|nabla|phi_core>
!----------------------------------------------------------------------------------

 already_has_nabla=all(pawtab(:)%has_nabla==3)
 call pawnabla_core_init(mpsang,dtset%ntypat,pawrad,pawtab,phi_cor,indlmn_core)

!----------------------------------------------------------------------------------
!2- Computation of <psi_n|p_i>(<phi_i|-i.nabla|phi_core>)
!----------------------------------------------------------------------------------

!Init parallelism
 spaceComm_w=mpi_enreg%comm_cell
 me=xmpi_comm_rank(spaceComm_w)
 if (mpi_enreg%paral_kgb==1) then
   spaceComm_k=mpi_enreg%comm_kpt
   spaceComm_bandspin=mpi_enreg%comm_bandspinor
   me_kpt=mpi_enreg%me_kpt
 else
   spaceComm_k=spaceComm_w
   spaceComm_bandspin=0
   me_kpt=me
 end if
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

!Determine if cprj datastructure is distributed over bands
 mband_cprj=mcprj/(my_nspinor*mkmem*nsppol)
 cprj_paral_band=(mband_cprj<mband)

!Prepare temporary PAW file if mkmem==0
 iorder_cprj=0
!Open _OPT2 file for proc 0
 iomode=IO_MODE_FORTRAN_MASTER
 fformopt=611
 ount = get_unit()
 call WffOpen(iomode,spaceComm_k,dtfil%fnameabo_app_opt2,ierr,wff1,0,me,ount)
 call hdr_io(fformopt,hdr,2,wff1)
 if (me==0) then
   write(ount)(eigen0(ib),ib=1,mband*nkpt*nsppol)
   write(ount) nphicor
   do iln=1,nphicor
     write(ount) ncor(iln),lcor(iln),energy_cor(iln)
   end do
 end if

 ABI_DEALLOCATE(ncor)
 ABI_DEALLOCATE(lcor)
 ABI_DEALLOCATE(phi_cor)
 ABI_DEALLOCATE(energy_cor)

 ABI_ALLOCATE(psinablapsi,(2,3,mband,nphicor,natom))

!LOOP OVER SPINS
 ibg=0
 bdtot_index=0
 do isppol=1,nsppol

!  LOOP OVER k POINTS
   do ikpt=1,nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
     etiq=ikpt+(isppol-1)*nkpt
     psinablapsi=zero

     mykpt=.true.
     mykpt=.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_kpt))
     if (mykpt) then

!      Constants depending on k-point
       istwf_k=dtset%istwfk(ikpt)
       cplex=2;if (istwf_k>1) cplex=1

!      Extract cprj for this k-point.
       nband_cprj_k=nband_k;if (cprj_paral_band) nband_cprj_k=nband_k/mpi_enreg%nproc_band
       if (mkmem*nsppol/=1) then
         ABI_DATATYPE_ALLOCATE(cprj_k_loc,(natom,my_nspinor*nband_cprj_k))
         call pawcprj_alloc(cprj_k_loc,0,dimcprj)
         call pawcprj_get(atindx1,cprj_k_loc,cprj,natom,1,ibg,ikpt,iorder_cprj,isppol,&
&         mband_cprj,mkmem,natom,nband_cprj_k,nband_cprj_k,my_nspinor,nsppol,dtfil%unpaw,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
       else
         cprj_k_loc => cprj
       end if

!      if cprj are distributed over bands, gather them (because we need to mix bands)
       if (cprj_paral_band) then
         ABI_DATATYPE_ALLOCATE(cprj_k,(natom,my_nspinor*nband_k))
         call pawcprj_alloc(cprj_k,0,dimcprj)
         call pawcprj_mpi_allgather(cprj_k_loc,cprj_k,natom,my_nspinor*nband_cprj_k,dimcprj,0,&
&         mpi_enreg%nproc_band,mpi_enreg%comm_band,ierr,rank_ordered=.false.)
       else
         cprj_k => cprj_k_loc
       end if

       ABI_ALLOCATE(tnm,(2,3,nband_k,nphicor,natom))
       tnm=zero

!      Loops on bands
       do jb=1,nband_k

         if (mpi_enreg%paral_kgb==1) then
           if (mod(jb-1,mpi_enreg%nproc_band)/=mpi_enreg%me_band) cycle
         elseif (xmpi_paral==1) then
           if (abs(mpi_enreg%proc_distrb(ikpt,jb,isppol)-me_kpt)/=0) cycle
         end if
         jbsp=(jb-1)*my_nspinor

         if (cplex==1) then
           do ispinor=1,my_nspinor
             jbsp=jbsp+1
             do iatom=1,natom
               itypat=dtset%typat(iatom)
               lmn_size=pawtab(itypat)%lmn_size
               do jlmn=1,lmn_size
                 do ilmn=1,lmncmax
                   ib=indlmn_core(5,ilmn)
                   cpnm1=cprj_k(iatom,jbsp)%cp(1,jlmn)
                   tnm(2,:,jb,ib,iatom)=tnm(2,:,jb,ib,iatom)+cpnm1*pawtab(itypat)%nabla_ij(:,jlmn,ilmn)
                 end do !ilmn
               end do !jlmn
             end do !iatom
           end do !ispinor
         else
           do ispinor=1,my_nspinor
             jbsp=jbsp+1
             do iatom=1,natom
               itypat=dtset%typat(iatom)
               lmn_size=pawtab(itypat)%lmn_size
               do jlmn=1,lmn_size
                 do ilmn=1,lmncmax
                   ib=indlmn_core(5,ilmn)
                   cpnm1=cprj_k(iatom,jbsp)%cp(1,jlmn)
                   cpnm2=cprj_k(iatom,jbsp)%cp(2,jlmn)
                   tnm(1,:,jb,ib,iatom)=tnm(1,:,jb,ib,iatom)+cpnm1*pawtab(itypat)%nabla_ij(:,jlmn,ilmn)
                   tnm(2,:,jb,ib,iatom)=tnm(2,:,jb,ib,iatom)+cpnm2*pawtab(itypat)%nabla_ij(:,jlmn,ilmn)
                 end do !ilmn
               end do !jlmn
             end do !iatom
           end do !ispinor
         end if

!        End loops on bands
       end do ! jb

!      Reduction in case of parallelism
       if (mpi_enreg%paral_kgb==1) then
         call timab(48,1,tsec)
         call xmpi_sum_master(tnm,0,spaceComm_bandspin,ierr)
         call timab(48,2,tsec)
       end if

       psinablapsi(:,:,:,:,:)=psinablapsi(:,:,:,:,:)+tnm(:,:,:,:,:)
       ABI_DEALLOCATE(tnm)

       if (mkmem/=0) ibg = ibg + my_nspinor*nband_cprj_k

       if (cprj_paral_band) then
         call pawcprj_free(cprj_k)
         ABI_DATATYPE_DEALLOCATE(cprj_k)
       end if
       if (mkmem*nsppol/=1) then
         call pawcprj_free(cprj_k_loc)
         ABI_DATATYPE_DEALLOCATE(cprj_k_loc)
       end if

       if (me==0) then
         do iatom=1,natom
           write(ount) ((psinablapsi(1:2,1,ib,jb,iatom),ib=1,nband_k),jb=1,nphicor)
           write(ount) ((psinablapsi(1:2,2,ib,jb,iatom),ib=1,nband_k),jb=1,nphicor)
           write(ount) ((psinablapsi(1:2,3,ib,jb,iatom),ib=1,nband_k),jb=1,nphicor)
         end do
!DEBUG
!         do iatom=1,natom
!           write(138,*) ((psinablapsi(1:2,1,ib,jb,iatom),ib=1,nband_k),jb=1,nphicor)
!           write(138,*) ((psinablapsi(1:2,2,ib,jb,iatom),ib=1,nband_k),jb=1,nphicor)
!           write(138,*) ((psinablapsi(1:2,3,ib,jb,iatom),ib=1,nband_k),jb=1,nphicor)
!         end do
!DEBUG

       elseif (mpi_enreg%me_band==0.and.mpi_enreg%me_fft==0) then
         call xmpi_exch(psinablapsi,etiq,me_kpt,psinablapsi,0,spaceComm_k,ierr)
       end if

     elseif (me==0) then
       sender=minval(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol))
       call xmpi_exch(psinablapsi,etiq,sender,psinablapsi,0,spaceComm_k,ierr)
       do iatom=1,natom
         write(ount) ((psinablapsi(1:2,1,ib,jb,iatom),ib=1,nband_k),jb=1,nphicor)
         write(ount) ((psinablapsi(1:2,2,ib,jb,iatom),ib=1,nband_k),jb=1,nphicor)
         write(ount) ((psinablapsi(1:2,3,ib,jb,iatom),ib=1,nband_k),jb=1,nphicor)
       end do

!DEBUG
!      do iatom=1,natom
!         write(138,*) ((psinablapsi(1:2,1,ib,jb,iatom),ib=1,nband_k),jb=1,nphicor)
!         write(138,*) ((psinablapsi(1:2,2,ib,jb,iatom),ib=1,nband_k),jb=1,nphicor)
!         write(138,*) ((psinablapsi(1:2,3,ib,jb,iatom),ib=1,nband_k),jb=1,nphicor)
!       end do
!DEBUG

     end if ! mykpt

     bdtot_index=bdtot_index+nband_k

!    End loop on spin,kpt
   end do ! ikpt
 end do !isppol

!Close file
 call WffClose(wff1,ierr)

!Datastructures deallocations
 ABI_DEALLOCATE(indlmn_core)
 ABI_DEALLOCATE(psinablapsi)
 if (.not.already_has_nabla) then
   do itypat=1,dtset%ntypat
     if (allocated(pawtab(itypat)%nabla_ij)) then
       ABI_DEALLOCATE(pawtab(itypat)%nabla_ij)
       pawtab(itypat)%has_nabla=0
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
!!      destroy_mpi_enreg,hdr_free,hdr_io,hdr_read_from_fname,initmpi_seq
!!      kramerskronig,matrginv,metric,wffopen
!!
!! SOURCE

 subroutine linear_optics_paw(filnam,filnam_out,mpi_enreg_seq)

!Arguments -----------------------------------
!scalars
 character(len=fnlen),intent(in) :: filnam,filnam_out
 type(MPI_type),intent(inout) :: mpi_enreg_seq

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
 complex,allocatable :: eps_work(:)
 character(len=fnlen) :: filnam1,filnam_gen
 character(len=500) :: msg
 type(hdr_type) :: hdr
 type(wffile_type) :: wff1
!arrays
 real(dp),allocatable :: psinablapsi(:,:,:,:)

! *********************************************************************************

 DBG_ENTER("COLL")

!* Fake MPI_type for the sequential part.
!This routine should not be parallelized as communicating gbig and other
!tables takes more time than recalculating them in sequential.
 call initmpi_seq(MPI_enreg_seq)

!write(std_out,'(a)')' Give the name of the output file ...'
!read(std_in, '(a)') filnam_out
!write(std_out,'(a)')' The name of the output file is :',filnam_out

!Read data file
 if (open_file(filnam,msg,newunit=inpunt,form='formatted') /= 0 ) then
   MSG_ERROR(msg)
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
 call hdr_free(hdr)
 if (fform1 /= 610) then
   MSG_ERROR("Abinit8 requires an OPT file with fform = 610")
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
 ABI_ALLOCATE(kpts,(3,nkpt))
 ABI_ALLOCATE(wtk,(nkpt))
 kpts(:,:)=hdr%kptns(:,:)
 wtk(:)=hdr%wtk(:)
 nspinor=hdr%nspinor
 nsppol=hdr%nsppol
 occopt=hdr%occopt
 rprimd(:,:)=hdr%rprimd(:,:)
 rprimdinv(:,:) = rprimd(:,:)
 call matrginv(rprimdinv,3,3) ! need the inverse of rprimd to symmetrize the tensors
 ABI_ALLOCATE(nband,(nkpt*nsppol))
 ABI_ALLOCATE(occ,(bantot))
 occ(1:bantot)=hdr%occ(1:bantot)
 nband(1:nkpt*nsppol)=hdr%nband(1:nkpt*nsppol)
 nsym=hdr%nsym
 ABI_ALLOCATE(symrel,(3,3,nsym))
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
 ABI_ALLOCATE(eigen0,(mband*nkpt*nsppol))
 read(wfunt)(eigen0(iband),iband=1,mband*nkpt*nsppol)

 read(inpunt,*)mbpt_sciss
 read(inpunt,*)dom,omin,omax,mom
 close(inpunt)

 ABI_ALLOCATE(oml1,(mom))
 ABI_ALLOCATE(e1,(3,3,mom))
 ABI_ALLOCATE(e2,(2,3,3,mom))
 ABI_ALLOCATE(epsilon_tot,(2,3,3,mom))
 ABI_ALLOCATE(eps_work,(mom))
 del=(omax-omin)/(mom-1)
 do iom=1,mom
   oml1(iom)=omin+dble(iom-1)*del
 end do
 write(std_out,'(a,i8,4f10.5,a)')' npts,omin,omax,width,mbpt_sciss      =',mom,omin,omax,dom,mbpt_sciss,' Ha'

 ABI_ALLOCATE(psinablapsi,(2,3,mband,mband))

!loop over spin components
 do isppol=1,nsppol
   bdtot_index = 0
!  loop over k points
   do ikpt=1,nkpt
!
!    number of bands for this k point
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     ABI_ALLOCATE(eig0_k,(nband_k))
     ABI_ALLOCATE(occ_k,(nband_k))
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

     ABI_DEALLOCATE(eig0_k)
     ABI_DEALLOCATE(occ_k)
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
     eps_work(:) = cmplx(0.0,epsilon_tot(2,ii,jj,:))
     call kramerskronig(mom,oml1,eps_work,method,only_check)
     epsilon_tot(1,ii,jj,:) = real(eps_work(:))
     if (ii /= jj) epsilon_tot(1,ii,jj,:) = epsilon_tot(1,ii,jj,:)- 1.0
   end do ! end loop over jj
 end do ! end loop over ii

 if (open_file(trim(filnam_out)//'_imag',msg,newunit=reunt,form='formatted') /= 0) then
   MSG_ERROR(msg)
 end if

 if (open_file(trim(filnam_out)//'_real',msg,unit=imunt,form='formatted') /= 0) then
   MSG_ERROR(msg)
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

 ABI_DEALLOCATE(nband)
 ABI_DEALLOCATE(oml1)
 ABI_DEALLOCATE(e2)
 ABI_DEALLOCATE(e1)
 ABI_DEALLOCATE(occ)
 ABI_DEALLOCATE(psinablapsi)
 ABI_DEALLOCATE(eigen0)
 ABI_DEALLOCATE(wtk)
 ABI_DEALLOCATE(kpts)

 call hdr_free(hdr)
 call destroy_mpi_enreg(MPI_enreg_seq)

 DBG_EXIT("COLL")

 end subroutine linear_optics_paw
!!***

!----------------------------------------------------------------------

END MODULE m_paw_optics
!!***
