!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_nsteltwf
!! NAME
!! dfpt_nsteltwf
!!
!! FUNCTION
!! This routine computes the non-local and kinetic contribution to the
!! 2DTE matrix elements, in the non-stationary formulation
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (DRH,XG,AR,MB,MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of RF
!!    wavefunctions at k,q.
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)  (NOT NEEDED !)
!!  effmass=effective mass for electrons (1. in common case)
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cg
!!  icg1=shift to be applied on the location of data in the array cg1
!!  ikpt=number of the k-point
!!  isppol=1 for unpolarized, 2 for spin-polarized
!!  istwf_k=flag controlling the storage of WFs
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kg1_k(3,npw1_k)=reduced planewave coordinates at k+q, with RF k points
!!  kpoint(3)=k-point in reduced coordinates
!!  mband=maximum number of bands
!!  mkmem =number of k points treated by this node.
!!  mk1mem =number of k points treated by this node (RF data).
!!  mpert =maximum number of ipert
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  natom=number of atoms in cell.
!!  nband_k=number of bands at this k point for that spin polarization
!!  npw_k=number of plane waves at this k point
!!  npw1_k=number of plane waves at this k+q point
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rmet(3,3)=real space metric (bohr**2)
!!  wtk_k=weight assigned to the k point.
!!  ylm(npw_k,mpsang*mpsang)= real spherical harmonics for each G and k point
!!  ylmgr(npw_k,3,mpsang*mpsang*useylm)= gradients of real spherical for each G and k point
!!
!! OUTPUT
!!  d2nl_k(2,3,mpert)=non-local contributions to
!!   non-stationary 2DTE, for the present k point, and perturbation idir, ipert
!!
!! PARENTS
!!      dfpt_nselt
!!
!! CHILDREN
!!      dotprod_g,kpgstr,load_k_hamiltonian,mkffnl,mkkin,nonlop
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_nsteltwf(cg,cg1,d2nl_k,ecut,ecutsm,effmass,gs_hamk,icg,icg1,ikpt,isppol,&
&  istwf_k,kg_k,kg1_k,kpoint,mband,mkmem,mk1mem,mpert,mpi_enreg,mpw,mpw1,natom,nband_k,&
&  npw_k,npw1_k,nspinor,nsppol,ntypat,occ_k,psps,rmet,wtk_k,ylm,ylmgr)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_wffile
 use m_xmpi
 use m_profiling_abi
 use m_cgtools

 use m_pawcprj, only : pawcprj_type
 use m_hamiltonian, only : gs_hamiltonian_type,load_k_hamiltonian
 use m_nonlop

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_nsteltwf'
 use interfaces_56_recipspace
 use interfaces_66_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,icg1,ikpt,isppol,istwf_k,mband,mk1mem,mkmem,mpert,mpw,mpw1,natom
 integer,intent(in) :: nspinor,nsppol,ntypat
 integer,intent(inout) :: nband_k,npw1_k,npw_k
 real(dp),intent(in) :: ecut,ecutsm,effmass,wtk_k
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: kg1_k(3,npw1_k),kg_k(3,npw_k)
 real(dp),intent(in) :: kpoint(3)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: occ_k(nband_k),rmet(3,3)
 real(dp),intent(in) :: ylm(npw_k,psps%mpsang*psps%mpsang)
 real(dp),intent(in) :: ylmgr(npw_k,3,psps%mpsang*psps%mpsang)
 real(dp),intent(inout) :: d2nl_k(2,3,mpert)

!Local variables-------------------------------
!scalars
 integer :: choice,cpopt,dimffnl,dimffnl2,iband
 integer :: ider,idir0,idir1,ilmn,ipert1,ipw,ipws,ispinor,istr1,itypat
 integer :: nkpg,nnlout,paw_opt,signs,tim_nonlop
 real(dp) :: doti,dotr
 type(gs_hamiltonian_type) :: gs_hamk
!arrays
 real(dp) :: enlout(6),dum_svectout(1,1),dum(1),kpg_dum(0,0)
 real(dp),allocatable :: cwave0(:,:),cwavef(:,:),dkinpw(:),eig2_k(:)
 real(dp),allocatable :: ffnl(:,:,:,:),ffnl_ylm(:,:,:,:),ghc(:,:)
 real(dp),allocatable :: gvnl1(:,:),gvnlc(:,:),kinpw1(:),ph3d(:,:,:)
 type(pawcprj_type) :: cprj_dum(0,0)

! *********************************************************************


!DEBUG
!write(std_out,*)' dfpt_nstwf : enter '
!ENDDEBUG

!Init me
 ABI_ALLOCATE(ghc,(2,npw1_k*nspinor))
 ABI_ALLOCATE(gvnlc,(2,npw1_k*nspinor))
 ABI_ALLOCATE(gvnl1,(2,npw1_k*nspinor))
 ABI_ALLOCATE(eig2_k,(2*nsppol*mband**2))
 ABI_ALLOCATE(kinpw1,(npw1_k))
 ABI_ALLOCATE(dkinpw,(npw_k))
 nkpg=0

!Compute nonlocal form factors ffnl at (k+G), for all atoms
 dimffnl=2
 ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,ntypat))
 if (psps%useylm==0) then
   ider=1;idir0=0
   call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,gs_hamk%gmet,gs_hamk%gprimd,ider,idir0,&
&   psps%indlmn,kg_k,kpg_dum,kpoint,psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
&   npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm,ylmgr)
 else
   ider=1;idir0=-7;dimffnl2=7
   ABI_ALLOCATE(ffnl_ylm,(npw_k,dimffnl2,psps%lmnmax,ntypat))
   call mkffnl(psps%dimekb,dimffnl2,psps%ekb,ffnl_ylm,psps%ffspl,gs_hamk%gmet,gs_hamk%gprimd,&
&   ider,idir0,psps%indlmn,kg_k,kpg_dum,kpoint,psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,&
&   nkpg,npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm,ylmgr)
   do itypat=1,ntypat
     do ilmn=1,psps%lmnmax
       ffnl(:,1,ilmn,itypat)=ffnl_ylm(:,1,ilmn,itypat)
     end do
   end do
 end if

!Compute kinetic contributions (1/2) (2 Pi)**2 (k+G)**2:
! call mkkin(ecut,ecutsm,effmass,gs_hamk%gmet,kg1_k,kinpw1,kpoint,npw1_k)
 call mkkin(ecut,ecutsm,effmass,gs_hamk%gmet,kg1_k,kinpw1,kpoint,npw1_k,0,0)

!Load k/k+q-dependent part in the Hamiltonian datastructure
 ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamk%matblk))
 call load_k_hamiltonian(gs_hamk,kpt_k=kpoint,npw_k=npw_k,istwf_k=istwf_k,kg_k=kg_k,ffnl_k=ffnl,&
& ph3d_k=ph3d,compute_ph3d=.true.)

 ABI_ALLOCATE(cwave0,(2,npw_k*nspinor))
 ABI_ALLOCATE(cwavef,(2,npw1_k*nspinor))

!Loop over bands
 do iband=1,nband_k

   if(mpi_enreg%proc_distrb(ikpt, iband, isppol) /= mpi_enreg%me_kpt) then
!    Skip the eigenvalue and the gvnl records of this band
     cycle
   end if

!  Get ground-state and first-order wavefunctions
   cwave0(:,:)=cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)
   cwavef(:,:)=cg1(:,1+(iband-1)*npw1_k*nspinor+icg1:iband*npw1_k*nspinor+icg1)

!  Double loop over strain perturbations
   do ipert1=natom+3,natom+4
     do idir1=1,3
       if (ipert1==natom+3) istr1=idir1
       if (ipert1==natom+4) istr1=idir1+3

!      Compute the derivative of the kinetic operator vs strain in dkinpw
       call kpgstr(dkinpw,ecut,ecutsm,effmass,gs_hamk%gmet,gs_hamk%gprimd,istr1,&
&       kg1_k,kpoint,npw1_k)

!      Get |vnon-locj1|u0> :
!      first-order non-local, applied to zero-order wavefunction
!      (??) this routine gives MINUS the non-local contribution

!      When using Ylms, load the correct ffnl derivative
       if (psps%useylm==1) then
         do itypat=1,ntypat
           do ilmn=1,psps%lmnmax
             ffnl(:,2,ilmn,itypat)=ffnl_ylm(:,1+istr1,ilmn,itypat)
           end do
         end do
       end if

       signs=2 ; choice=3 ; nnlout=6 ; paw_opt=0 ; cpopt=-1 ; tim_nonlop=5
       call nonlop(choice,cpopt,cprj_dum,enlout,gs_hamk,istr1,dum,mpi_enreg,1,nnlout,paw_opt,&
&       signs,dum_svectout,tim_nonlop,cwave0,gvnl1)
!      <G|Vnl1|Cnk> is contained in gvnl1

!      Kinetic contribution
       do ispinor=1,nspinor
         do ipw=1,npw1_k
           ipws=ipw+npw1_k*(ispinor-1)
           if(kinpw1(ipw)<huge(0.0_dp)*1.d-11)then
             gvnl1(1,ipws)=gvnl1(1,ipws)+dkinpw(ipw)*cwave0(1,ipws)
             gvnl1(2,ipws)=gvnl1(2,ipws)+dkinpw(ipw)*cwave0(2,ipws)
           else
             gvnl1(1,ipws)=0.0_dp
             gvnl1(2,ipws)=0.0_dp
           end if
         end do
       end do

!      construct the matrix element (<uj2|vj1|u0>)complex conjug.
!      and add it to the 2nd-order matrix
!      imaginary term should be zero for strain-strain 2nd derivatives,
!      but keep it as a test for now
       call dotprod_g(dotr,doti,gs_hamk%istwf_k,npw1_k*nspinor,2,cwavef,gvnl1,&
&       mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

       d2nl_k(1,idir1,ipert1)= d2nl_k(1,idir1,ipert1)+wtk_k*occ_k(iband)*2.0_dp*dotr
       d2nl_k(2,idir1,ipert1)= d2nl_k(2,idir1,ipert1)-wtk_k*occ_k(iband)*2.0_dp*doti

     end do !idir1
   end do !ipert1

!  UNTIL NOW, DO NOT TAKE INTO ACCOUNT istwf_k

!  End loop over bands
 end do

 ABI_DEALLOCATE(cwave0)
 ABI_DEALLOCATE(cwavef)

!###################################################################

 ABI_DEALLOCATE(eig2_k)
 ABI_DEALLOCATE(ghc)
 ABI_DEALLOCATE(gvnlc)
 ABI_DEALLOCATE(gvnl1)
 ABI_DEALLOCATE(kinpw1)
 ABI_DEALLOCATE(dkinpw)
 ABI_DEALLOCATE(ffnl)
 ABI_DEALLOCATE(ph3d)
 if (psps%useylm==1)  then
   ABI_DEALLOCATE(ffnl_ylm)
 end if

end subroutine dfpt_nsteltwf
!!***
