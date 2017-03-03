!{\src2tex{textfont=tt}}
!!****f* ABINIT/optics_paw_core
!! NAME
!! optics_paw_core
!!
!! FUNCTION
!! Compute matrix elements need for X spectr. (in the PAW context) and store them in a file
!!  Matrix elements = <Phi_core|Nabla|Phi_j>
!!
!! COPYRIGHT
!! Copyright (C) 2005-2017 ABINIT group (SM,MT)
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
!!      hdr_io,int_ang,nderiv_gen,pawcprj_alloc,pawcprj_free,pawcprj_get
!!      pawcprj_mpi_allgather,pawpsp_read_corewf,pawrad_deducer0,simp_gen,timab
!!      wffclose,wffopen,xmpi_exch,xmpi_sum_master
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine optics_paw_core(atindx1,cprj,dimcprj,dtfil,dtset,eigen0,filpsp,hdr,&
&               mband,mcprj,mkmem,mpi_enreg,mpsang,natom,nkpt,nsppol,pawrad,pawtab)


 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_wffile
 use m_hdr

 use m_io_tools,  only : get_unit
 use m_pawpsp,    only : pawpsp_read_corewf
 use m_pawtab,    only : pawtab_type
 use m_pawcprj,   only : pawcprj_type,  pawcprj_alloc, pawcprj_get, &
&                        pawcprj_free, pawcprj_mpi_allgather
 use m_pawrad,    only : pawrad_type, pawrad_init, pawrad_free, pawrad_ifromr, &
&                        pawrad_deducer0, simp_gen, nderiv_gen, bound_deriv
 use m_splines,   only : spline,splint

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'optics_paw_core'
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_65_paw, except_this_one => optics_paw_core
!End of the abilint section

 implicit none

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
 type(pawtab_type),target,intent(in) :: pawtab(dtset%ntypat)

!Local variables-------------------------------
!scalars
 integer :: basis_size,bdtot_index,cplex,etiq,iatom,ib,ibg
 integer :: ierr,ij_size,ikpt,il,ilm,ilmn,iln,ount
 integer :: iorder_cprj,ispinor,isppol,istwf_k,itypat
 integer :: jb,jbsp,jl,jlm,jlmn,jln,lmn_size,lmncmax,mband_cprj
 integer :: me,me_kpt,mesh_size,my_nspinor,nband_cprj_k,nband_k,nphicor
 integer :: sender,spaceComm_bandspin,spaceComm_k,spaceComm_w
 integer :: iomode,fformopt
 logical :: cprj_paral_band,ex,mykpt
 character(len=fnlen) :: filecore
 real(dp) :: cpnm1,cpnm2,intg
!arrays
 integer,allocatable :: indlmn_core(:,:),lcor(:),ncor(:)
 integer,ABI_CONTIGUOUS pointer :: indlmn(:,:)
 real(dp) :: ang_phipphj(mpsang**2,mpsang**2,8)
 real(dp) :: tsec(2)
 real(dp),allocatable :: dphi(:),energy_cor(:),ff(:),int1(:,:),int2(:,:)
 real(dp),allocatable :: rad(:),phi_cor(:,:),psinablapsi(:,:,:,:,:)
 real(dp),allocatable :: tnm(:,:,:,:,:)
 type(coeff3_type), allocatable :: phipphj(:)
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

!1-A Integration of the angular part : all angular integrals have been
!computed outside Abinit and tabulated for each (l,m) value
!----------------------------------------------------------------------------------

 call int_ang(ang_phipphj,mpsang)

 ABI_DATATYPE_ALLOCATE(phipphj,(dtset%ntypat))

!We consider the impurity to be the first atom
!loop on atoms type
 do itypat=1,dtset%ntypat

   mesh_size=pawtab(itypat)%mesh_size
   lmn_size=pawtab(itypat)%lmn_size
   basis_size=pawtab(itypat)%basis_size
   ij_size=lmn_size*lmn_size
   indlmn => pawtab(itypat)%indlmn

   ABI_ALLOCATE(ff,(mesh_size))
   ABI_ALLOCATE(rad,(mesh_size))
   ABI_ALLOCATE(int2,(lmn_size,lmncmax))
   ABI_ALLOCATE(int1,(lmn_size,lmncmax))
   ABI_ALLOCATE(dphi,(mesh_size))
   ABI_ALLOCATE(phipphj(itypat)%value,(3,lmn_size,lmncmax))

   rad(1:mesh_size)=pawrad(itypat)%rad(1:mesh_size)

!  1-B Computation of int1=\int phi phi_core /r dr
!  ----------------------------------------------------------------------------------
   do jln=1,nphicor
     do iln=1,basis_size
       ff(2:mesh_size)=(pawtab(itypat)%phi(2:mesh_size,iln)*phi_cor(2:mesh_size,jln))/rad(2:mesh_size)
       call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intg,ff,pawrad(itypat))
       int1(iln,jln)=intg
     end do
   end do

!  1-C Computation of int2=\int phi/r d/dr(phi_core/r) - phi phi_core/r dr
!  ----------------------------------------------------------------------------------
   do jln=1,nphicor
     ff(1:mesh_size)=phi_cor(1:mesh_size,jln)
     call nderiv_gen(dphi,ff,pawrad(itypat))

     do iln=1,basis_size
       ff(2:mesh_size)=pawtab(itypat)%phi(2:mesh_size,iln)*dphi(2:mesh_size) &
&       -pawtab(itypat)%phi(2:mesh_size,iln)*phi_cor(2:mesh_size,jln)/ &
&       rad(2:mesh_size)
       call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intg,ff,pawrad(itypat))
       int2(iln,jln)=intg
     end do
   end do

!  1-D Integration of the radial part
!  ----------------------------------------------------------------------------------
   do jlmn=1,lmncmax
     jlm=indlmn_core(4,jlmn)
     jl=indlmn_core(5,jlmn)
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

!  end loop on atoms type
 end do

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
     write(ount) ncor(iln),lcor(iln),half*energy_cor(iln)
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
                   tnm(2,:,jb,ib,iatom)=tnm(2,:,jb,ib,iatom)+cpnm1*phipphj(itypat)%value(:,jlmn,ilmn)
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
                   tnm(1,:,jb,ib,iatom)=tnm(1,:,jb,ib,iatom)+cpnm1*phipphj(itypat)%value(:,jlmn,ilmn)
                   tnm(2,:,jb,ib,iatom)=tnm(2,:,jb,ib,iatom)+cpnm2*phipphj(itypat)%value(:,jlmn,ilmn)
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
 do itypat=1,dtset%ntypat
   ABI_DEALLOCATE(phipphj(itypat)%value)
 end do
 ABI_DATATYPE_DEALLOCATE(phipphj)
 ABI_DEALLOCATE(indlmn_core)
 ABI_DEALLOCATE(psinablapsi)

 DBG_EXIT("COLL")

 end subroutine optics_paw_core
!!***
