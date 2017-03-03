!{\src2tex{textfont=tt}}
!!****f* ABINIT/optics_paw
!! NAME
!! optics_paw
!!
!! FUNCTION
!! Compute matrix elements need for optical conductivity (in the PAW context) and store them in a file
!!  Matrix elements = <Phi_i|Nabla|Phi_j>
!!
!! COPYRIGHT
!! Copyright (C) 2005-2017 ABINIT group (SM,VR,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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
!!      hdr_io,int_ang,nderiv_gen,pawcprj_alloc,pawcprj_free,pawcprj_get
!!      pawcprj_mpi_allgather,pawrad_deducer0,simp_gen,timab,wffclose,wffopen
!!      xmpi_exch,xmpi_sum,xmpi_sum_master
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine optics_paw(atindx1,cg,cprj,dimcprj,dtfil,dtset,eigen0,gprimd,hdr,kg,&
&               mband,mcg,mcprj,mkmem,mpi_enreg,mpsang,mpw,natom,nkpt,npwarr,nsppol,&
&               pawrad,pawtab)

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_wffile
 use m_profiling_abi
 use m_hdr

 use m_io_tools,  only : get_unit
 use m_pawrad,    only : pawrad_type, pawrad_deducer0, simp_gen, nderiv_gen
 use m_pawtab,    only : pawtab_type
 use m_pawcprj,   only : pawcprj_type, pawcprj_alloc, pawcprj_get, &
&                        pawcprj_free,pawcprj_mpi_allgather

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'optics_paw'
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_65_paw, except_this_one => optics_paw
!End of the abilint section

 implicit none

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
 type(pawtab_type),target,intent(in) :: pawtab(dtset%ntypat)

!Local variables-------------------------------
!scalars
 integer :: iomode,basis_size,bdtot_index,cplex,etiq,fformopt,iatom,ib,ibg,ibsp
 integer :: icg,ierr,ij_size,ikg,ikpt,il,ilm,ilmn,iln,ount
 integer :: iorder_cprj,ipw,ispinor,isppol,istwf_k,itypat,iwavef
 integer :: jb,jbsp,jl,jlm,jlmn,jln,jwavef,lmn_size,mband_cprj
 integer :: me,me_kpt,mesh_size,my_nspinor,nband_k,nband_cprj_k,npw_k,sender
 integer :: spaceComm_band,spaceComm_bandfftspin,spaceComm_fft,spaceComm_k,spaceComm_spin,spaceComm_w
 logical :: cprj_paral_band,mykpt
 real(dp) :: cgnm1,cgnm2,cpnm1,cpnm2,intg
 character(len=500) :: message
!arrays
 integer :: tmp_shape(3)
 integer,allocatable :: kg_k(:,:)
 integer,ABI_CONTIGUOUS pointer :: indlmn(:,:)
 real(dp) :: ang_phipphj(mpsang**2,mpsang**2,8),kpoint(3)
 real(dp) :: tsec(2)
 real(dp),allocatable :: dphi(:),dtphi(:)
 real(dp),allocatable ::ff(:),int1(:,:),int2(:,:)
 real(dp),allocatable :: kpg_k(:,:),phidphj(:,:),psinablapsi(:,:,:,:)
 real(dp),allocatable :: rad(:),tnm(:,:,:,:),tphidtphj(:,:)
 type(coeff3_type), allocatable :: phipphj(:)
 type(pawcprj_type),pointer :: cprj_k(:,:),cprj_k_loc(:,:)
 type(wffile_type) :: wff1

! ************************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 ABI_CHECK(mkmem/=0,"mkmem==0 not supported anymore!")

!----------------------------------------------------------------------------------
!1- Computation of phipphj=<phi_i|nabla|phi_j>-<tphi_i|nabla|tphi_j>
!----------------------------------------------------------------------------------

!1-A Integration of the angular part : all angular integrals have been
!computed outside Abinit and tabulated for each (l,m) value
!----------------------------------------------------------------------------------

 call int_ang(ang_phipphj,mpsang)

 ABI_DATATYPE_ALLOCATE(phipphj,(dtset%ntypat))

!loop on atoms type
 do itypat=1,dtset%ntypat

   mesh_size=pawtab(itypat)%mesh_size
   lmn_size=pawtab(itypat)%lmn_size
   basis_size=pawtab(itypat)%basis_size
   ij_size=lmn_size*lmn_size

   ABI_ALLOCATE(ff,(mesh_size))
   ABI_ALLOCATE(rad,(mesh_size))
   ABI_ALLOCATE(int2,(lmn_size,lmn_size))
   ABI_ALLOCATE(int1,(lmn_size,lmn_size))
   ABI_ALLOCATE(dphi,(mesh_size))
   ABI_ALLOCATE(dtphi,(mesh_size))
   ABI_ALLOCATE(phidphj,(mesh_size,ij_size))
   ABI_ALLOCATE(tphidtphj,(mesh_size,ij_size))
   ABI_ALLOCATE(phipphj(itypat)%value,(3,lmn_size,lmn_size))

   indlmn => pawtab(itypat)%indlmn
   rad(1:mesh_size)=pawrad(itypat)%rad(1:mesh_size)

!  1-B  Computation of int1=\int phi phj /r dr - \int tphi tphj /r dr
!  ----------------------------------------------------------------------------------
   do jln=1,basis_size
     do iln=1,basis_size
       ff(2:mesh_size)=(pawtab(itypat)%phi(2:mesh_size,iln)*pawtab(itypat)%phi(2:mesh_size,jln)&
&       -pawtab(itypat)%tphi(2:mesh_size,iln)*pawtab(itypat)%tphi(2:mesh_size,jln))/rad(2:mesh_size)
       call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intg,ff,pawrad(itypat))
       int1(iln,jln)=intg
     end do
   end do

!  1-C Computation of int2=\int phi/r d/dr(phj/r) r^2dr - \int tphi/r d/dr(tphj/r)r^2 dr
!  ----------------------------------------------------------------------------------
   do jln=1,basis_size

     ff(1:mesh_size)=pawtab(itypat)%phi(1:mesh_size,jln)
     call nderiv_gen(dphi,ff,pawrad(itypat))
     ff(1:mesh_size)=pawtab(itypat)%tphi(1:mesh_size,jln)
     call nderiv_gen(dtphi,ff,pawrad(itypat))

     do iln=1,basis_size
       ff(2:mesh_size)=pawtab(itypat)%phi(2:mesh_size,iln)*dphi(2:mesh_size) &
&       -pawtab(itypat)%phi (2:mesh_size,iln)*pawtab(itypat)%phi(2:mesh_size,jln)/ &
&       rad(2:mesh_size)-(pawtab(itypat)%tphi(2:mesh_size,iln)*dtphi(2:mesh_size) &
&       -pawtab(itypat)%tphi (2:mesh_size,iln)*pawtab(itypat)%tphi(2:mesh_size,jln)/rad(2:mesh_size))
       call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intg,ff,pawrad(itypat))
       int2(iln,jln)=intg
     end do
   end do

!  1-D Integration of the radial part
!  ----------------------------------------------------------------------------------
   do jlmn=1,lmn_size
     jlm=indlmn(4,jlmn)
     jl=indlmn(5,jlmn)
     do ilmn=1,lmn_size
       ilm=indlmn(4,ilmn)
       il=indlmn(5,ilmn)
       phipphj(itypat)%value(1,ilmn,jlmn)= int2(il,jl)*ang_phipphj(ilm,jlm,1)&
&       + int1(il,jl)*(ang_phipphj(ilm,jlm,2)+ang_phipphj(ilm,jlm,3))
       phipphj(itypat)%value(2,ilmn,jlmn)= int2(il,jl)*ang_phipphj(ilm,jlm,4)&
&       + int1(il,jl)*(ang_phipphj(ilm,jlm,5)+ang_phipphj(ilm,jlm,6))
       phipphj(itypat)%value(3,ilmn,jlmn)= int2(il,jl)*ang_phipphj(ilm,jlm,7)&
&       + int1(il,jl)*ang_phipphj(ilm,jlm,8)
     end do
   end do

   ABI_DEALLOCATE(ff)
   ABI_DEALLOCATE(rad)
   ABI_DEALLOCATE(int2)
   ABI_DEALLOCATE(int1)
   ABI_DEALLOCATE(dphi)
   ABI_DEALLOCATE(dtphi)
   ABI_DEALLOCATE(phidphj)
   ABI_DEALLOCATE(tphidtphj)

!  end loop on atoms type
 end do

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
                     tnm(2,:,ib,jb)=tnm(2,:,ib,jb)+cpnm1*phipphj(itypat)%value(:,ilmn,jlmn)
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
                     tnm(1,:,ib,jb)=tnm(1,:,ib,jb)+cpnm2*phipphj(itypat)%value(:,ilmn,jlmn)
                     tnm(2,:,ib,jb)=tnm(2,:,ib,jb)-cpnm1*phipphj(itypat)%value(:,ilmn,jlmn)
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
 do itypat=1,dtset%ntypat
   ABI_DEALLOCATE(phipphj(itypat)%value)
 end do
 ABI_DATATYPE_DEALLOCATE(phipphj)
 ABI_DEALLOCATE(psinablapsi)

 DBG_EXIT("COLL")

 end subroutine optics_paw
!!***
