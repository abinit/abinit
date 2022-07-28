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
 use m_fstrings,   only : itoa, sjoin
 use m_io_tools,   only : file_exists
 use m_time,       only : cwtime
 use m_kg,         only : mkkpg
 use m_mpinfo,     only : proc_distrb_cycle
 use m_dfptlw_wf
 use m_dfpt_mkvxc, only : dfpt_mkvxcggadq
 use m_dfptlw_nv,  only : dfptlw_geom
 use m_spacepar,   only : hartredq
 use m_cgtools,    only : dotprod_vn
 use m_mkffnl,     only : mkffnl

 implicit none

 public :: dfptlw_pert
 public :: preca_ffnl

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
!!   - 1st-order Local+SCF potentials for i1pert and i2pert 
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
!!  dimffnl= third dimension of ffnl
!!  d3e_pert1(mpert)=array with the i1pert cases to calculate
!!  d3e_pert2(mpert)=array with the i2pert cases to calculate
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen1(2*mband*mband*nkpt*nsppol)=1st-order eigenvalues for i1pert,i1dir (hartree)
!!  eigen2(2*mband*mband*nkpt*nsppol)=1st-order eigenvalues for i2pert,i2dir (hartree)
!!  ffnl(dtset%mkmem,dtset%mpw,dimffnl,psps%lmnmax,psps%ntypat)= Nonlocal projectors and their derivatives
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
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
!!  nylmgr=second dimension of ylmgr_k
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  rhog(2,nfft)=array for Fourier transform of GS electron density
!!  rho1g1(2,nfft)=G-space RF electron density in electrons/bohr**3 (i1pert)
!!  rhor(nfft,nspden)=array for GS electron density in electrons/bohr**3.
!!  rho1r1(cplex*nfft,nspden)=RF electron density in electrons/bohr**3 (i1pert)
!!  rho2r1(cplex*nfft,nspden)=RF electron density in electrons/bohr**3 (i2pert)
!!  rmet(3,3)=real space metric tensor in bohr**2
!!  rprimd(3,3) = dimensional primitive translations (bohr)
!!  samepert= .true. if i1pert=i2pert and i1dir=i2dir
!!  ucvol=volume of the unit cell
!!  useylmgr= if 1 use the derivative of spherical harmonics
!!  vpsp1_i1pertdq(cplex*nfft,nspden,n1dq)= local potential of first-order
!!          gradient Hamiltonian for i1pert
!!  vpsp1_i1pertdq(cplex*nfft,nspden,n1dq)= local potential of second-order
!!          gradient Hamiltonian for i1pert
!!  vpsp1_i1pertdq_geom(cplex*nfft,nspden,3)= local potential of first-order
!!          gradient Hamiltonian for i1pert wrp to i3dir and i2dir
!!  vpsp1_i2pertdq(cplex*nfft,nspden,n2dq)= local potential of first-order
!!          gradient Hamiltonian for i2pert
!!  ddk_f = wf files
!!  d2_dkdk_f = wf files
!!  xccc3d1(cplex*n3xccc)=3D change in core charge density (dummy) 
!!  xred(3,natom) = reduced atomic coordinates
!!  ylm(mpw*mkmem,psps%mpsang*psps%mpsang*psps%useylm)=real spherical harmonics
!!  ylmgr(mpw*mkmem,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr)= k-gradients of real spherical harmonics
!!
!! OUTPUT
!!  d3etot(2,3,mpert,3,mpert,3,mpert) = third derivatives of the energy tensor
!!  d3etot_t4(2,n2dq)= t4 term which might need to be converted to type-II
!!  d3etot_t5(2,n1dq)= t5 term which might need to be converted to type-II
!!  d3etot_tgeom(2,2)= Geometric term which needs to be converted to type-II
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

subroutine dfptlw_pert(atindx,cg,cg1,cg2,cplex,d3e_pert1,d3e_pert2,d3etot,d3etot_t4,d3etot_t5,d3etot_tgeom,&
& dimffnl,dtfil,dtset,eigen1,eigen2,ffnl,gmet,gs_hamkq,gsqcut,i1dir,i2dir,i3dir,&
& i1pert,i2pert,i3pert,kg,kxc,mband,mgfft,mkmem_rbz,mk1mem,mpert,mpi_enreg,mpsang,mpw,natom,nattyp,&
& n1dq,n2dq,nfft,ngfft,nkpt,nkxc,&
& nspden,nspinor,nsppol,npwarr,nylmgr,occ,pawfgr,ph1d,psps,rhog,rho1g1,rhor,rho1r1,rho2r1,rmet,rprimd,samepert,&
& ucvol,useylmgr,vpsp1_i1pertdq,vpsp1_i1pertdqdq,vpsp1_i1pertdq_geom,vpsp1_i2pertdq,ddk_f,d2_dkdk_f,xccc3d1,xred,ylm,ylmgr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,dimffnl,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,mband,mgfft
 integer,intent(in) :: mk1mem,mkmem_rbz,mpert,mpsang,mpw,natom,n1dq,n2dq,nfft,nkpt,nkxc,nspden
 integer,intent(in) :: nspinor,nsppol,nylmgr,useylmgr
 real(dp),intent(in) :: gsqcut,ucvol
 logical,intent(in) :: samepert
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(pawfgr_type),intent(in) :: pawfgr
 type(wfk_t),intent(inout) :: ddk_f,d2_dkdk_f

!arrays
 integer,intent(in) :: atindx(natom),kg(3,mpw*mkmem_rbz),nattyp(psps%ntypat),ngfft(18),npwarr(nkpt)
 integer,intent(in) :: d3e_pert1(mpert),d3e_pert2(mpert)
 real(dp),intent(in) :: eigen1(2*mband*mband*nkpt*nsppol)
 real(dp),intent(in) :: eigen2(2*mband*mband*nkpt*nsppol)
 real(dp),intent(in) :: ffnl(mkmem_rbz,mpw,dimffnl,psps%lmnmax,psps%ntypat)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem_rbz*nsppol)
 real(dp),intent(in) :: cg1(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: cg2(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: gmet(3,3),kxc(nfft,nkxc)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: rhog(2,nfft),rhor(nfft,dtset%nspden)
 real(dp),intent(in) :: rho1g1(2,nfft),rho1r1(cplex*nfft,dtset%nspden)
 real(dp),intent(in) :: rho2r1(cplex*nfft,dtset%nspden)
 real(dp),intent(in) :: rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: xccc3d1(cplex*nfft),xred(3,natom)
 real(dp),intent(in) :: vpsp1_i1pertdq(2*nfft,nspden,n1dq)
 real(dp),intent(in) :: vpsp1_i1pertdqdq(2*nfft,nspden,n2dq)
 real(dp),intent(in) :: vpsp1_i1pertdq_geom(2*nfft,nspden,3)
 real(dp),intent(in) :: vpsp1_i2pertdq(2*nfft,nspden,n2dq)
 real(dp),intent(inout) :: d3etot(2,3,mpert,3,mpert,3,mpert)
 real(dp),intent(out) :: d3etot_t4(2,n2dq),d3etot_t5(2,n1dq)
 real(dp),intent(out) :: d3etot_tgeom(2,n2dq)
 real(dp),intent(in) :: ylm(mpw*mk1mem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(mpw*mk1mem,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr)

!Variables ------------------------------------
!scalars
 integer :: bandtot,bd2tot,icg,idq,ierr,ii,ikc,ikg,ikpt,ilm,isppol,istwf_k,me,n1,n2,n3,n4,n5,n6 
 integer :: nband_k,npw_k,spaceworld,tim_getgh1c
 integer :: usepaw
 real(dp) :: tmpim,tmpre,wtk_k
 real(dp) :: cpu, wall, gflops
 character(len=1000) :: msg
 logical :: with_nonlocal_i1pert, with_nonlocal_i2pert
!arrays
 integer,allocatable :: kg_k(:,:)
 real(dp) :: d3etot_t1(2),d3etot_t1_k(2)
 real(dp) :: d3etot_t2(2),d3etot_t2_k(2)
 real(dp) :: d3etot_t3(2),d3etot_t3_k(2)
 real(dp) :: d3etot_t4_k(2,n2dq)
 real(dp) :: d3etot_t5_k(2,n1dq)
 real(dp) :: d3etot_tgeom_k(2,n2dq)
 real(dp) :: d3etot_telec(2)
 real(dp) :: e3tot(2),kpt(3)
 real(dp),allocatable :: eig1_k(:),eig2_k(:),occ_k(:)
 real(dp),allocatable :: dum_vlocal(:,:,:,:),dum_vpsp(:)
 real(dp),allocatable :: vlocal1dq(:,:,:,:)
 real(dp),allocatable :: vlocal1(:,:,:,:)
 real(dp),allocatable :: vpsp1(:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr_k(:,:,:)
 real(dp),allocatable :: ffnl_k(:,:,:,:)
 type(pawcprj_type),allocatable :: dum_cwaveprj(:,:)

 !AZ_test_ini**************************************************************************
 character(40) :: i1dir_text, i2dir_text, i3dir_text, iband_text, ikpt_text, file_name
 !AZ_test_fin**************************************************************************
 
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

!Initialize d3etot parts
 d3etot_t1=zero
 d3etot_t2=zero
 d3etot_t3=zero
 d3etot_t4=zero
 d3etot_t5=zero
 d3etot_telec=zero
 d3etot_tgeom=zero

!Calculate the electrostatic contribution 
 call lw_elecstic(cplex,d3etot_telec,gmet,gs_hamkq%gprimd,gsqcut,&
& i1dir,i2dir,i3dir,i1pert,i2pert,i3pert,&
& kxc,mpi_enreg,nfft,ngfft,nkxc,nspden,rho1g1,rho1r1,rho2r1,ucvol)
    
 !AZ_test_ini******************************************************************
 write(i1dir_text,'(i8)') i1dir
 write(i2dir_text,'(i8)') i2dir
 write(i3dir_text,'(i8)') i3dir
 file_name = 'T_elec_i1dir_'//trim(adjustl(i1dir_text))//'_i2dir_'//&
 & trim(adjustl(i2dir_text))//'_i3dir_'//trim(adjustl(i3dir_text))//'.dat'
 open(unit=999,file=file_name,action='write',status='replace')
 write(999,'(2f12.6)') d3etot_telec(1)*four*pi, d3etot_telec(2)*four*pi
 close(999)
 !AZ_test_fin******************************************************************
 
!Loop over spins
 bandtot = 0
 bd2tot = 0
 icg=0
 do isppol = 1, nsppol

!  Loop over k-points
   ikg = 0
   ikc = 0
   do ikpt = 1, nkpt

     nband_k = dtset%nband(ikpt+(isppol-1)*nkpt)
     npw_k = npwarr(ikpt)
     istwf_k = dtset%istwfk(ikpt)

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,mband,isppol,mpi_enreg%me)) then
       bandtot = bandtot + nband_k
       bd2tot = bd2tot + 2*nband_k**2
       cycle ! Skip the rest of the k-point loop
     end if
     ikc= ikc + 1

     ABI_MALLOC(occ_k,(nband_k))
     occ_k(:) = occ(1+bandtot:nband_k+bandtot)
     wtk_k    = dtset%wtk(ikpt)
     kpt(:) = dtset%kptns(:,ikpt)

     ABI_MALLOC(eig1_k,(2*nband_k**2))
     ABI_MALLOC(eig2_k,(2*nband_k**2))
     ABI_MALLOC(kg_k,(3,npw_k))
     ABI_MALLOC(ylm_k,(npw_k,mpsang*mpsang*psps%useylm))
     ABI_MALLOC(ylmgr_k,(npw_k,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr))
     ABI_MALLOC(ffnl_k,(npw_k,dimffnl,psps%lmnmax,psps%ntypat))

     !Get plane-wave vectors and related data at k
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     if (dtset%ffnl_lw==1) then
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
     else if (dtset%ffnl_lw==0) then
       ffnl_k(1:npw_k,:,:,:)=ffnl(ikc,1:npw_k,:,:,:)
     end if

     !Get matrix elements for uniform perturbations
     eig1_k(:)=eigen1(1+bd2tot:2*nband_k**2+bd2tot)
     eig2_k(:)=eigen2(1+bd2tot:2*nband_k**2+bd2tot)

     
     !Compute the stationary terms of d3etot depending on response functions
     call dfpt_1wf(atindx,cg,cg1,cg2,cplex,ddk_f,d2_dkdk_f,d3etot_t1_k,d3etot_t2_k,d3etot_t3_k,& 
     & d3etot_t4_k,d3etot_t5_k,dimffnl,dtset,eig1_k,eig2_k,ffnl_k,gs_hamkq,gsqcut,icg,&
     & i1dir,i2dir,i3dir,i1pert,i2pert,i3pert,ikpt,isppol,istwf_k,&
     & kg_k,kpt,kxc,mkmem_rbz,mpi_enreg,mpw,natom,nattyp,nband_k,&
     & n1dq,n2dq,nfft,ngfft,nkxc,npw_k,nspden,nsppol,nylmgr,occ_k,&
     & pawfgr,ph1d,psps,rhog,rhor,rmet,rprimd,samepert,ucvol,useylmgr,&
     & vpsp1_i1pertdq,vpsp1_i2pertdq,&
     & wtk_k,xred,ylm_k,ylmgr_k)

!    Add the contribution from each k-point. 
     d3etot_t1=d3etot_t1 + d3etot_t1_k
     d3etot_t2=d3etot_t2 + d3etot_t2_k
     d3etot_t3=d3etot_t3 + d3etot_t3_k
     d3etot_t4=d3etot_t4 + d3etot_t4_k
     d3etot_t5=d3etot_t5 + d3etot_t5_k
 
     !Compute the nonvariational geometric term
     call cwtime(cpu, wall, gflops, "start")
     if (i1pert<=natom.and.(i2pert==natom+3.or.i2pert==natom+4)) then
       call dfptlw_geom(atindx,cg,d3etot_tgeom_k,dimffnl,dtset, &
       &  ffnl_k,gs_hamkq,gsqcut,icg, &
       &  i1dir,i2dir,i3dir,i1pert,i2pert,ikpt, &
       &  isppol,istwf_k,kg_k,kpt,mkmem_rbz,mpi_enreg,natom,mpw,nattyp,nband_k,n2dq,nfft, &
       &  ngfft,npw_k,nspden,nsppol,nylmgr,occ_k, &
       &  ph1d,psps,rmet,rprimd,ucvol,useylmgr,vpsp1_i1pertdqdq,vpsp1_i1pertdq_geom,wtk_k,ylm_k,ylmgr_k)

       !Add the contribution from each k-point
       d3etot_tgeom=d3etot_tgeom + d3etot_tgeom_k
     end if
     call cwtime(cpu, wall, gflops, "stop")

!    Keep track of total number of bands
     bandtot = bandtot + nband_k
     bd2tot = bd2tot + 2*nband_k**2

!    Shift arrays memory
     icg=icg+npw_k*dtset%nspinor*nband_k
     ikg=ikg+npw_k

     ABI_FREE(eig1_k)
     ABI_FREE(eig2_k)
     ABI_FREE(occ_k)
     ABI_FREE(kg_k)
     ABI_FREE(ylm_k)
     ABI_FREE(ylmgr_k)
     ABI_FREE(ffnl_k)

   end do !ikpt

 end do !isppol


!=== MPI communications ==================
 if (xmpi_paral==1) then
   call xmpi_sum(d3etot_t1,spaceworld,ierr)
   call xmpi_sum(d3etot_t2,spaceworld,ierr)
   call xmpi_sum(d3etot_t3,spaceworld,ierr)
   call xmpi_sum(d3etot_t4,spaceworld,ierr)
   call xmpi_sum(d3etot_t5,spaceworld,ierr)
   call xmpi_sum(d3etot_tgeom,spaceworld,ierr)
 end if

!Apply +i or -i in case of strain perturbation.
 if (i1pert==natom+3.or.i1pert==natom+4) then
   tmpre=d3etot_telec(1);tmpim=d3etot_telec(2) ; d3etot_telec(1)=tmpim;d3etot_telec(2)=-tmpre
   tmpre=d3etot_t1(1);tmpim=d3etot_t1(2) ; d3etot_t1(1)=tmpim;d3etot_t1(2)=-tmpre
   tmpre=d3etot_t2(1);tmpim=d3etot_t2(2) ; d3etot_t2(1)=tmpim;d3etot_t2(2)=-tmpre
   tmpre=d3etot_t3(1);tmpim=d3etot_t3(2) ; d3etot_t3(1)=tmpim;d3etot_t3(2)=-tmpre
   do idq=1,n2dq
     tmpre=d3etot_t4(1,idq);tmpim=d3etot_t4(2,idq) ; d3etot_t4(1,idq)=tmpim;d3etot_t4(2,idq)=-tmpre
   end do
   do idq=1,n1dq
     tmpre=d3etot_t5(1,idq);tmpim=d3etot_t5(2,idq) ; d3etot_t5(1,idq)=tmpim;d3etot_t5(2,idq)=-tmpre
   end do
 end if
 if (i2pert==natom+3.or.i2pert==natom+4) then
   tmpre=d3etot_telec(1);tmpim=d3etot_telec(2) ; d3etot_telec(1)=-tmpim;d3etot_telec(2)=tmpre
   tmpre=d3etot_t1(1);tmpim=d3etot_t1(2) ; d3etot_t1(1)=-tmpim;d3etot_t1(2)=tmpre
   tmpre=d3etot_t2(1);tmpim=d3etot_t2(2) ; d3etot_t2(1)=-tmpim;d3etot_t2(2)=tmpre
   tmpre=d3etot_t3(1);tmpim=d3etot_t3(2) ; d3etot_t3(1)=-tmpim;d3etot_t3(2)=tmpre
   do idq=1,n2dq
     tmpre=d3etot_t4(1,idq);tmpim=d3etot_t4(2,idq) ; d3etot_t4(1,idq)=-tmpim;d3etot_t4(2,idq)=tmpre
     if (i1pert<=natom) then
       tmpre=d3etot_tgeom(1,idq);tmpim=d3etot_tgeom(2,idq) ; d3etot_tgeom(1,idq)=-tmpim;d3etot_tgeom(2,idq)=tmpre
     end if
   end do
   do idq=1,n1dq
     tmpre=d3etot_t5(1,idq);tmpim=d3etot_t5(2,idq) ; d3etot_t5(1,idq)=-tmpim;d3etot_t5(2,idq)=tmpre
   end do
 end if

!Join all the contributions in e3tot except t4 and t5 which may need to be
!converted to type-II in case of strain perturbation. 
!Apply here the two factor to the stationary wf1 contributions 
!(see PRB 105, 064101 (2022))
 d3etot_t1(:)=two*d3etot_t1(:)
 d3etot_t2(:)=two*d3etot_t2(:)
 d3etot_t3(:)=two*d3etot_t3(:)
 d3etot_t4(:,:)=two*d3etot_t4(:,:)
 d3etot_t5(:,:)=two*d3etot_t5(:,:)
 e3tot(:)=d3etot_t1(:)+d3etot_t2(:)+d3etot_t3(:)+d3etot_telec(:)


!Before printing, set small contributions to zero
 !Real parts
 if (abs(d3etot_t1(1))<tol8) d3etot_t1(1)= zero
 if (abs(d3etot_t2(1))<tol8) d3etot_t2(1)= zero
 if (abs(d3etot_t3(1))<tol8) d3etot_t3(1)= zero
 do idq=1,n2dq
   if (abs(d3etot_t4(1,idq))<tol8) d3etot_t4(1,idq)= zero
   if (abs(d3etot_tgeom(1,idq))<tol8) d3etot_tgeom(1,idq)= zero
 end do 
 do idq=1,n1dq
   if (abs(d3etot_t5(1,idq))<tol8) d3etot_t5(1,idq)= zero
 end do 
 if (abs(d3etot_telec(1))<tol8) d3etot_telec(1)= zero
 if (abs(e3tot(1))    <tol8)     e3tot(1)= zero

 !Imaginary parts
 if (abs(d3etot_t1(2))<tol8) d3etot_t1(2)= zero
 if (abs(d3etot_t2(2))<tol8) d3etot_t2(2)= zero
 if (abs(d3etot_t3(2))<tol8) d3etot_t3(2)= zero
 do idq=1,n2dq
   if (abs(d3etot_t4(2,idq))<tol8) d3etot_t4(2,idq)= zero
   if (abs(d3etot_tgeom(2,idq))<tol8) d3etot_tgeom(2,idq)= zero
 end do
 do idq=1,n1dq
   if (abs(d3etot_t5(2,idq))<tol8) d3etot_t5(2,idq)= zero
 end do
 if (abs(d3etot_telec(2))<tol8) d3etot_telec(2)= zero
 if (abs(e3tot(2))    <tol8)     e3tot(2)= zero

 if (dtset%prtvol>=10) then
   write(msg,'(4(a,2(a,f18.8)),a)') &
   ch10,'          d3etot_telec = ',d3etot_telec(1),  ',',d3etot_telec(2),&
   ch10,'             d3etot_t1 = ',d3etot_t1(1),  ',',d3etot_t1(2),&
   ch10,'             d3etot_t2 = ',d3etot_t2(1),  ',',d3etot_t2(2),&
   ch10,'             d3etot_t3 = ',d3etot_t3(1),  ',',d3etot_t3(2)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   if (n2dq==1) then
     write(msg,'(2(a,f18.8))') &
     '             d3etot_t4 = ',d3etot_t4(1,1),  ',',d3etot_t4(2,1)
   else if (n2dq==2) then
     write(msg,'(2(2(a,f18.8)a))') &
     '   d3etot_t4(dw shear) = ',d3etot_t4(1,1),  ',',d3etot_t4(2,1),ch10,&
     '   d3etot_t4(up shear) = ',d3etot_t4(1,2),  ',',d3etot_t4(2,2)
   end if
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   if (n1dq==1) then
     write(msg,'(2(a,f18.8))') &
     '             d3etot_t5 = ',d3etot_t5(1,1),  ',',d3etot_t5(2,1)
   else if (n1dq==2) then
     write(msg,'(2(2(a,f18.8)a))') &
     '   d3etot_t5(dw shear) = ',d3etot_t5(1,1),  ',',d3etot_t5(2,1),ch10,&
     '   d3etot_t5(up shear) = ',d3etot_t5(1,2),  ',',d3etot_t5(2,2)
   end if
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   if (i1pert<=natom.and.(i2pert==natom+3.or.i2pert==natom+4)) then
     if (n2dq==1) then
       write(msg,'(2(a,f18.8))') &
       '          d3etot_tgeom = ',d3etot_tgeom(1,1),  ',',d3etot_tgeom(2,1)
     else if (n2dq==2) then
       write(msg,'(2(2(a,f18.8)a))') &
       'd3etot_tgeom(dw shear) = ',d3etot_tgeom(1,1),  ',',d3etot_tgeom(2,1),ch10,&
       'd3etot_tgeom(up shear) = ',d3etot_tgeom(1,2),  ',',d3etot_tgeom(2,2)
     end if
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
   end if
 end if

 d3etot(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)=e3tot(:)

 DBG_EXIT("COLL")

end subroutine dfptlw_pert
!!***

!!****f* ABINIT/m_dfptlw_pert/lw_elecstic
!! NAME
!!  lw_elecstic
!!
!! FUNCTION
!!  This routine calculates the electrostatic term of the spatial-dispersion
!!  third-order energy derivative for a couple of perturbations and a gradient
!!  direction.
!!
!! COPYRIGHT
!!  Copyright (C) 2022 ABINIT group (MR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!          if 2, COMPLEX
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  gsqcut=large sphere cut-off
!!  i1dir,i2dir,i3dir=directions of the corresponding perturbations
!!  i1pert,i2pert,i3pert = type of perturbation that has to be computed
!!  kxc(nfft,nkxc)=exchange and correlation kernel
!!  mpi_enreg=information about MPI parallelization
!!  nfft= number of FFT grid points (for this proc) 
!!  ngfft(1:18)=integer array with FFT box dimensions and other 
!!  nkxc=second dimension of the kxc array. If /=0, the XC kernel must be computed.
!!  nspden = number of spin-density components
!!  rho1g1(2,nfft)=G-space RF electron density in electrons/bohr**3 (i1pert)
!!  rho1r1(cplex*nfft,nspden)=RF electron density in electrons/bohr**3 (i1pert)
!!  rho2r1(cplex*nfft,nspden)=RF electron density in electrons/bohr**3 (i2pert)
!!  ucvol=volume of the unit cell
!!
!! OUTPUT
!!  d3etot_telec(2)= Electrostatic term of the third-order energy derivative
!!
!! SIDE EFFECTS
!!
!! NOTES
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


subroutine lw_elecstic(cplex,d3etot_telec,gmet,gprimd,gsqcut,&
& i1dir,i2dir,i3dir,i1pert,i2pert,i3pert,&
& kxc,mpi_enreg,nfft,ngfft,nkxc,nspden,rho1g1,rho1r1,rho2r1,ucvol)
    
 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: cplex,i1dir,i2dir,i3dir,i1pert,i2pert,i3pert
 integer,intent(in) :: nfft,nkxc,nspden
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(in) :: rho1g1(2,nfft),rho1r1(cplex*nfft,nspden)
 real(dp),intent(in) :: rho2r1(cplex*nfft,nspden),kxc(nfft,nkxc)
 real(dp),intent(out) :: d3etot_telec(2),gprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,nfftot,qcar
 real(dp) :: doti,dotr

!arrays
 real(dp),allocatable :: rhor1_cplx(:,:)
 real(dp),allocatable :: vxc1dq(:,:),vxc1dq_car(:,:,:),vqgradhart(:)
 
! *************************************************************************

 DBG_ENTER("COLL")
 
!If GGA xc first calculate the Cartesian q gradient of the xc kernel
 if (nkxc == 7) then
   ABI_MALLOC(vxc1dq,(2*nfft,nspden))
   ABI_MALLOC(vxc1dq_car,(2*nfft,nspden,3))
   do qcar=1,3
     call dfpt_mkvxcggadq(cplex,gprimd,kxc,mpi_enreg,nfft,ngfft,nkxc,nspden,qcar,rho1r1,vxc1dq)
     vxc1dq_car(:,:,qcar)=vxc1dq(:,:)
   end do
 end if

!Calculate the q gradient of the Hartree potential
 ABI_MALLOC(vqgradhart,(2*nfft))
 call hartredq(2,gmet,gsqcut,mpi_enreg,nfft,ngfft,i3dir,rho1g1,vqgradhart)

!If GGA convert the gradient of xc kernel to reduced coordinates and incorporate it to the Hartree part
 if (nkxc == 7) then
   vxc1dq=zero
   do qcar=1,3
     vxc1dq(:,:)=vxc1dq(:,:) + gprimd(qcar,i3dir) * vxc1dq_car(:,:,qcar)
   end do
   vqgradhart(:)=vqgradhart(:)+vxc1dq(:,1)
   ABI_FREE(vxc1dq_car)
 end if

!Calculate the electrostatic energy term with the i2pert density response
!I need a complex density for the dotprod_vn
 ABI_MALLOC(rhor1_cplx,(2*nfft,nspden))
 rhor1_cplx=zero
 do ii=1,nfft
   jj=ii*2
   rhor1_cplx(jj-1,:)=rho2r1(ii,:)
 end do

 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 call dotprod_vn(2,rhor1_cplx,dotr,doti,nfft,nfftot,nspden,2,vqgradhart,ucvol)
 
 d3etot_telec(1)=dotr
 d3etot_telec(2)=doti

!Deallocations
 if (nkxc == 7) ABI_FREE(vxc1dq)
 ABI_FREE(vqgradhart)
 ABI_FREE(rhor1_cplx)

 DBG_EXIT("COLL")

end subroutine lw_elecstic
!!***

!!****f* ABINIT/m_dfptlw_pert/preca_ffnl
!! NAME
!!  preca_ffnl
!!
!! FUNCTION
!!  Calculates the nonlocal form factors and derivatives for all the atoms
!!  and k points.
!!
!! COPYRIGHT
!!  Copyright (C) 2022 ABINIT group (MR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dimffnl= second dimension of ffnl
!!  gmet(3,3)= reciprocal-space metric tensor
!!  gprimd(3,3)= dimensional reciprocal space primitive translations (b^-1)
!!  ider= if 1 first order derivatives of ffnl are calculated
!!        if 2 first and second order derivatives of ffnl are calculated
!!  idir0= variable that controls the way in which the derivatives of ffnl are
!!         calculated and saved
!!  kg(3,mpw)=integer coordinates of G vectors in basis sphere
!!  kptns(3,nkpt)=k points in terms of reciprocal translations
!!  mband= masimum number of bands
!!  mkmem= maximum number of k points which can fit in core memory
!!  mpi_enreg=information about MPI parallelization
!!  mpw   = maximum number of planewaves in basis sphere (large number)
!!  nkpt = number of k point
!!  npwarr(nkpt)=array holding npw for each k point
!!  nylmgr=second dimension of ylmgr
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  rmet(3,3)= real-space metric tensor
!!  useylmgr= if 1 use the derivative of spherical harmonics
!!  ylm(mpw*mkmem,psps%mpsang*psps%mpsang*psps%useylm)=real spherical harmonics
!!  ylmgr(mpw*mkmem,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr)= k-gradients of real spherical harmonics
!!
!! OUTPUT
!!  ffnl(mkmem,npw_k,dimffnl,psps%lmnmax,psps%ntypat)= Nonlocal projectors and their derivatives
!!
!! SIDE EFFECTS
!!
!! NOTES
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


subroutine preca_ffnl(dimffnl,ffnl,gmet,gprimd,ider,idir0,kg,kptns,mband,mkmem,mpi_enreg,mpw,nkpt, &
& npwarr,nylmgr,psps,rmet,useylmgr,ylm,ylmgr)
    
 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

!Arguments ------------------------------------
!scalars
 integer , intent(in)  :: dimffnl,ider,idir0,mband,mkmem,mpw,nkpt,nylmgr,useylmgr
 type(pseudopotential_type),intent(in) :: psps
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: kg(3,mpw*mkmem)
 integer,intent(in) :: npwarr(nkpt)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),kptns(3,nkpt),rmet(3,3)
 real(dp),intent(in) :: ylm(mpw*mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(mpw*mkmem,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr)
 real(dp),intent(out) :: ffnl(mkmem,mpw,dimffnl,psps%lmnmax,psps%ntypat)                       

!Local variables-------------------------------
!scalars
 integer :: ii,ikc,ikg,ikpt,ilm,nkpg,npw_k
 character(len=500) :: msg
!arrays
 integer,allocatable :: kg_k(:,:)
 real(dp) :: kpt(3)
 real(dp),allocatable :: ffnl_k(:,:,:,:),kpg_k(:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr_k(:,:,:),ylmgr_k_part(:,:,:)
 
! *************************************************************************

 DBG_ENTER("COLL")

 !Loop over k-points
 ikg=0
 ikc=0
 do ikpt = 1, nkpt
 
   npw_k = npwarr(ikpt)
   if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,mband,1,mpi_enreg%me)) then
     cycle ! Skip the rest of the k-point loop
   end if
   ikc= ikc + 1

   ABI_MALLOC(kg_k,(3,npw_k))
   ABI_MALLOC(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
   ABI_MALLOC(ylmgr_k,(npw_k,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr))

 
   kpt(:)= kptns(:,ikpt)

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

   if (dimffnl==4) then
     ABI_MALLOC(ylmgr_k_part,(npw_k,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr))
     ylmgr_k_part(:,:,:)=ylmgr_k(:,1:3,:)
   else if (dimffnl==10) then
     ABI_MALLOC(ylmgr_k_part,(npw_k,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr))
     ylmgr_k_part(:,:,:)=ylmgr_k(:,:,:)
   else 
     msg='wrong size for ffnl via dimffnl!'
     ABI_BUG(msg)
   end if

   nkpg=0
   ABI_MALLOC(kpg_k,(npw_k,nkpg))
   ABI_MALLOC(ffnl_k,(npw_k,dimffnl,psps%lmnmax,psps%ntypat))
   call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl_k,psps%ffspl,gmet,gprimd,ider,idir0,&
 & psps%indlmn,kg_k,kpg_k,kpt,psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
 & npw_k,psps%ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_k_part)

   ffnl(ikc,1:npw_k,:,:,:)=ffnl_k(1:npw_k,:,:,:)

   ABI_FREE(kg_k)
   ABI_FREE(ylm_k)
   ABI_FREE(ylmgr_k)
   ABI_FREE(ylmgr_k_part)
   ABI_FREE(ffnl_k)
   ABI_FREE(kpg_k)

   !Shift arrays memory
   ikg=ikg+npw_k

 end do

 DBG_EXIT("COLL")

end subroutine preca_ffnl
!!***

end module m_dfptlw_pert
!!***
