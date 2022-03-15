!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dfptlw_pert
!! NAME
!!  m_dfptlw_pert
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2022 ABINIT group (MR)
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

module m_dfptlw_pert
    
 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use m_dtset
 use m_dtfil
 use m_errors
 use m_profiling_abi
 use m_hamiltonian
 use m_cgtools
 use m_pawcprj
 use m_pawfgr
 use m_wfk
 use m_xmpi
 use m_getgh1c
 use m_mklocl
 use m_initylmg,   only : initylmg
 use m_fstrings, only : itoa, sjoin
 use m_io_tools, only : file_exists
 use m_time, only : cwtime
 use m_kg, only : mkkpg
 use m_mpinfo, only : proc_distrb_cycle
 use m_dfptlw_wf

 implicit none

 public :: dfptlw_pert

 private

! *************************************************************************

contains 
!!***

!!****f* ABINIT/m_dfptlw_pert/dfptlw_pert
!! NAME
!!  dfptlw_pert
!!
!! FUNCTION
!! Compute first-order response function contributions to the spatial-dispersion
!! 3rd order energy derivatives of the longwave driver.
!! The main inputs are :
!!   - GS WFs and Hamiltonian (cg,gs_hamkq)
!!   - 1st-order WFs for two perturbations i1pert/i1dir,i2pert/i2dir (cg1,cg2)
!!   - 1st-order Local+SCF potentials for i1pert and i2pert (vtrial1_i1pert,vtrial1_i2pert)
!!   - 1st-order WFs DDK and 2nd-order WF D2_DKDK (d2_dkdk_f)
!!
!! COPYRIGHT
!! Copyright (C) 2018-2021 ABINIT group (MR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cg(2,mpw*nspinor*mband*mkmem_rbz*nsppol) = array for planewave
!!                                          coefficients of wavefunctions
!!  cg1 = first derivative of cg with respect the perturbation i1pert
!!  cg2 = first derivative of cg with respect the perturbation i2pert
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!          if 2, COMPLEX
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  gsqcut=large sphere cut-off
!!  i1dir,i2dir,i3dir=directions of the corresponding perturbations
!!  i1pert,i2pert,i3pert = type of perturbation that has to be computed
!!  kg(3,mpw*mkmem_rbz)=reduced planewave coordinates
!!  kxc(nfft,nkxc)=exchange and correlation kernel
!!  mband = maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem_rbz = maximum number of k points which can fit in core memory
!!  mk1mem = maximum number of k points for first-order WF
!!           which can fit in core memory
!!  mpert =maximum number of ipert
!!  mpi_enreg=MPI-parallelisation information
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw   = maximum number of planewaves in basis sphere (large number)
!!  natom = number of atoms in unit cell
!!  nattyp(ntypat)= # atoms of each type.
!!  n1dq= third dimension of vlocal1_i1pertdq
!!  n2dq= third dimension of vlocal1_i2pertdq
!!  nfft= number of FFT grid points (for this proc) 
!!  ngfft(1:18)=integer array with FFT box dimensions and other 
!!  nkpt = number of k points
!!  nkxc=second dimension of the kxc array. If /=0, the XC kernel must be computed.
!!  nspden = number of spin-density components
!!  nspinor = number of spinorial components of the wavefunctions
!!  nsppol = number of channels for spin-polarization (1 or 2)
!!  npwarr(nkpt) = array holding npw for each k point
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  rhog(2,nfft)=array for Fourier transform of GS electron density
!!  rho1g1(2,nfft)=G-space RF electron density in electrons/bohr**3 (i1pert)
!!  rhor(nfft,nspden)=array for GS electron density in electrons/bohr**3.
!!  rho2r1(cplex*nfft,nspden)=RF electron density in electrons/bohr**3 (i2pert)
!!  rmet(3,3)=real space metric tensor in bohr**2
!!  rprimd(3,3) = dimensional primitive translations (bohr)
!!  ucvol=volume of the unit cell
!!  vpsp1_i1pertdq(cplex*nfft,nspden,n1dq)= local potential of first-order
!!          gradient Hamiltonian for i1pert
!!  vpsp1_i2pertdq(cplex*nfft,nspden,n2dq)= local potential of first-order
!!          gradient Hamiltonian for i2pert
!!  vtrial1_i1pert(cplex*nfft,nspden)=firs-order local potential
!!  vtrial1_i2pert(cplex*nfft,nspden)=firs-order local potential
!!  ddk_f = wf files
!!  d2_dkdk_f = wf files
!!  xccc3d1(cplex*n3xccc)=3D change in core charge density (dummy) 
!!  xred(3,natom) = reduced atomic coordinates
!!
!! OUTPUT
!!  d3etot(2,3,mpert,3,mpert,3,mpert) = third derivatives of the energy tensor
!!
!! SIDE EFFECTS
!!  TO DO!
!!
!! PARENTS
!!      m_dfptlw_loop
!!
!! CHILDREN
!!      dotprod_vn
!!
!! SOURCE

subroutine dfptlw_pert(atindx,cg,cg1,cg2,cplex,dtfil,dtset,d3etot,gs_hamkq,gsqcut,i1dir,i2dir,i3dir,&
& i1pert,i2pert,i3pert,kg,kxc,mband,mgfft,mkmem_rbz,mk1mem,mpert,mpi_enreg,mpsang,mpw,natom,nattyp,&
& n1dq,n2dq,nfft,ngfft,nkpt,nkxc,&
& nspden,nspinor,nsppol,npwarr,occ,pawfgr,ph1d,psps,rhog,rho1g1,rhor,rho2r1,rmet,rprimd,&
& ucvol,vpsp1_i1pertdq,vpsp1_i2pertdq,vtrial1_i1pert,vtrial1_i2pert,ddk_f,d2_dkdk_f,xccc3d1,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,mband,mgfft
 integer,intent(in) :: mk1mem,mkmem_rbz,mpert,mpsang,mpw,natom,n1dq,n2dq,nfft,nkpt,nkxc,nspden
 integer,intent(in) :: nspinor,nsppol
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(pawfgr_type),intent(in) :: pawfgr
 type(wfk_t),intent(inout) :: ddk_f,d2_dkdk_f

!arrays
 integer,intent(in) :: atindx(natom),kg(3,mpw*mkmem_rbz),nattyp(psps%ntypat),ngfft(18),npwarr(nkpt)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem_rbz*nsppol)
 real(dp),intent(in) :: cg1(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: cg2(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: kxc(nfft,nkxc)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: rhog(2,nfft),rhor(nfft,dtset%nspden)
 real(dp),intent(in) :: rho1g1(2,nfft),rho2r1(cplex*nfft,dtset%nspden)
 real(dp),intent(in) :: rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: xccc3d1(cplex*nfft),xred(3,natom)
 real(dp),intent(in) :: vpsp1_i1pertdq(2*nfft,nspden,n1dq)
 real(dp),intent(in) :: vpsp1_i2pertdq(2*nfft,nspden,n2dq)
 real(dp),intent(in) :: vtrial1_i1pert(cplex*nfft,nspden)
 real(dp),intent(in) :: vtrial1_i2pert(cplex*nfft,nspden)
 real(dp),intent(inout) :: d3etot(2,3,mpert,3,mpert,3,mpert)

!Variables ------------------------------------
!scalars
 integer :: bandtot,icg,idq,ierr,ii,ikg,ikpt,ilm,isppol,istwf_k,me,n1,n2,n3,n4,n5,n6 
 integer :: nband_k,npw_k,nylmgr,option,spaceworld,tim_getgh1c
 integer :: usepaw,useylmgr
 real(dp) :: wtk_k
 character(len=1000) :: msg
 logical :: with_nonlocal_i1pert, with_nonlocal_i2pert
!arrays
 integer,allocatable :: kg_k(:,:)
 real(dp) :: buffer(10)
 real(dp) :: d3etot_k(2)
 real(dp) :: d3etot_t1(2),d3etot_t1_k(2)
 real(dp) :: d3etot_t2(2),d3etot_t2_k(2)
 real(dp) :: d3etot_t3(2),d3etot_t3_k(2)
 real(dp) :: d3etot_t4(2,n2dq),d3etot_t4_k(2,n2dq)
 real(dp) :: d3etot_t5(2,n1dq),d3etot_t5_k(2,n1dq)
 real(dp) :: d3etot_telec(2),d3etot_telec_k(2)
 real(dp) :: e3tot(2),kpt(3)
 real(dp),allocatable :: occ_k(:)
 real(dp),allocatable :: dum_vlocal(:,:,:,:),dum_vpsp(:)
 real(dp),allocatable :: vlocal1dq(:,:,:,:)
 real(dp),allocatable :: vlocal1(:,:,:,:)
 real(dp),allocatable :: vpsp1(:)
 real(dp),allocatable :: ylm(:,:),ylmgr(:,:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr_k(:,:,:)
 type(pawcprj_type),allocatable :: dum_cwaveprj(:,:)
 
! *************************************************************************

 DBG_ENTER("COLL")

write(msg,'(2a,3(a,i2,a,i1))') ch10,'LONGWAVE : ',&
 ' perts : ',i1pert,'.',i1dir,' / ',i2pert,'.',i2dir,' / ',i3pert,'.',i3dir
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

!Init parallelism
 spaceworld=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt 

!Additional definitions
 tim_getgh1c=0
 usepaw=dtset%usepaw
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 n4=ngfft(4) ; n5=ngfft(5) ; n6=ngfft(6)
 with_nonlocal_i1pert=.true. ; if (i1pert==natom+2) with_nonlocal_i1pert=.false.
 with_nonlocal_i2pert=.true. ; if (i2pert==natom+2) with_nonlocal_i2pert=.false.

!Additional allocations
 ABI_MALLOC(dum_vpsp,(nfft))
 ABI_MALLOC(dum_vlocal,(n4,n5,n6,gs_hamkq%nvloc))
 ABI_MALLOC(vlocal1,(cplex*n4,n5,n6,gs_hamkq%nvloc))
 ABI_MALLOC(vlocal1dq,(2*n4,n5,n6,gs_hamkq%nvloc))
 ABI_MALLOC(vpsp1,(cplex*nfft))
 ABI_MALLOC(dum_cwaveprj,(0,0))

!Set up the spherical harmonics (Ylm) and gradients at each k point 
 useylmgr=1; option=2 ; nylmgr=9
 ABI_MALLOC(ylm,(mpw*mkmem_rbz,psps%mpsang*psps%mpsang*psps%useylm))               
 ABI_MALLOC(ylmgr,(mpw*mkmem_rbz,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr))
 if (psps%useylm==1) then
   call initylmg(gs_hamkq%gprimd,kg,dtset%kptns,mkmem_rbz,mpi_enreg,&
 & psps%mpsang,mpw,dtset%nband,dtset%nkpt,npwarr,dtset%nsppol,option,&
 & rprimd,ylm,ylmgr)                                   
 end if

!Initialize d3etot parts
d3etot_t1=zero
d3etot_t2=zero
d3etot_t3=zero
d3etot_t4=zero
d3etot_t5=zero
d3etot_telec=zero

!Loop over spins
 bandtot = 0
 icg=0
 do isppol = 1, nsppol

!  Loop over k-points
   ikg = 0
   do ikpt = 1, nkpt

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,mband,isppol,mpi_enreg%me)) then
       cycle ! Skip the rest of the k-point loop
     end if

     nband_k = dtset%nband(ikpt+(isppol-1)*nkpt)
     npw_k = npwarr(ikpt)
     istwf_k = dtset%istwfk(ikpt)
     ABI_MALLOC(occ_k,(nband_k))
     occ_k(:) = occ(1+bandtot:nband_k+bandtot)
     wtk_k    = dtset%wtk(ikpt)
     kpt(:) = dtset%kptns(:,ikpt)

     ABI_MALLOC(kg_k,(3,npw_k))
     ABI_MALLOC(ylm_k,(npw_k,mpsang*mpsang*psps%useylm))
     ABI_MALLOC(ylmgr_k,(npw_k,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr))
     !Get plane-wave vectors and related data at k
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     if (psps%useylm==1) then
       do ilm=1,psps%mpsang*psps%mpsang
         ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       end do
       if (useylmgr==1) then
         do ilm=1,psps%mpsang*psps%mpsang
           do ii=1,nylmgr
             ylmgr_k(1:npw_k,ii,ilm)=ylmgr(1+ikg:npw_k+ikg,ii,ilm)
           end do
         end do
       end if
     end if

     !Compute the stationary terms of d3etot depending on response functions
     call dfpt_1wf(atindx,cg,cg1,cg2,cplex,ddk_f,d2_dkdk_f,d3etot_t1_k,d3etot_t2_k,d3etot_t3_k,& 
     & d3etot_t4_k,d3etot_t5_k,dtset,gs_hamkq,gsqcut,icg,&
     & i1dir,i2dir,i3dir,i1pert,i2pert,i3pert,ikpt,isppol,istwf_k,&
     & kg_k,kpt,kxc,mkmem_rbz,mpi_enreg,mpw,natom,nattyp,nband_k,&
     & n1dq,n2dq,nfft,ngfft,nkxc,npw_k,nspden,nsppol,nylmgr,occ_k,&
     & pawfgr,ph1d,psps,rhog,rhor,rmet,ucvol,useylmgr,&
     & vpsp1_i1pertdq,vpsp1_i2pertdq,vtrial1_i1pert,vtrial1_i2pert,&
     & wtk_k,xred,ylm_k,ylmgr_k)

!    Add the contribution from each k-point. 
     d3etot_t1=d3etot_t1 + d3etot_t1_k
 
!    Keep track of total number of bands
     bandtot = bandtot + nband_k

!    Shift arrays memory
     icg=icg+npw_k*dtset%nspinor*nband_k
     ikg=ikg+npw_k

     ABI_FREE(occ_k)
     ABI_FREE(kg_k)
     ABI_FREE(ylm_k)
     ABI_FREE(ylmgr_k)

   end do !ikpt

 end do !isppol

!=== MPI communications ==================
 if (xmpi_paral==1) then

 ! Real parts
   buffer(1)=d3etot_t1(1)

 ! Imaginary parts
   buffer(6)=d3etot_t1(2)

   call xmpi_sum(buffer,spaceworld,ierr)

 ! Real parts
   d3etot_t1(1)=buffer(1)

 ! Imaginary parts
   d3etot_t1(2)=buffer(6)

 end if

!Join all the contributions in e3tot. Apply here the two factor to 
!the stationary wf1 contributions (see PRB 105, 064101 (2022))
 d3etot_t1(:)=two*d3etot_t1(:)
 e3tot(:)=d3etot_t1(:)

!Before printing, set small contributions to zero
 if (dtset%kptopt==3) then

   !Real parts
   if (abs(d3etot_t1(1))<tol8) d3etot_t1(1)= zero
   if (abs(e3tot(1))    <tol8)     e3tot(1)= zero
   !Imaginary parts
   if (abs(d3etot_t1(2))<tol8) d3etot_t1(2)= zero
   if (abs(e3tot(2))    <tol8)     e3tot(2)= zero

 else if (dtset%kptopt==2) then

   !Real parts
   d3etot_t1(1)= zero
   e3tot(1)   = zero
   !Imaginary parts
   if (abs(d3etot_t1(2))<tol8) d3etot_t1(2)= zero
   if (abs(e3tot(2))    <tol8)     e3tot(2)= zero

 else

   write(msg,"(1a)") 'kptopt must be 2 or 3 for the longwave calculation'
   ABI_BUG(msg)

 end if

 if (dtset%prtvol>=10) then
   write(msg,'(2(a,2(a,f18.8)),a)') &
   ch10,'        d3etot_t1 = ',d3etot_t1(1),  ',',d3etot_t1(2),&
   ch10,'           d3etot = ',e3tot(1),      ',',e3tot(2), ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if

 d3etot(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)=e3tot(:)

 DBG_EXIT("COLL")

end subroutine dfptlw_pert
!!***

end module m_dfptlw_pert
!!***
