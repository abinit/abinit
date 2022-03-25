!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dfptlw_nv
!! NAME
!!  m_dfptlw_nv
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2022 ABINIT group (MR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_dfptlw_nv
    
 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use m_abicore
 use m_xmpi
 use m_errors
 use m_mpinfo
 use m_dtset
 use m_hamiltonian
 use m_cgtools
 use m_wfk
 use m_xmpi
 use m_getgh1c
 use m_mklocl
 use m_pawcprj
 use m_pawfgr

 use m_dfpt_elt,    only : dfpt_ewalddq, dfpt_ewalddqdq
 use m_kg,          only : mkkpg
 use m_dynmat,      only : cart39

 implicit none

 private
!!***

 public :: dfptlw_nv
 public :: dfptlw_geom
!!***

! *************************************************************************

contains 
!!***

!!****f* ABINIT/m_dfptlw_nv/dfptlw_nv
!! NAME
!!  dfptlw_nv
!!
!! FUNCTION
!!  This routine calculates the nonvariational Ewald contributions to the 
!!  spatial-dispersion third-order energy derivatives.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2 
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1)
!!  mpert=maximum number of ipert
!!  my_natom=number of atoms treated by current processor
!!  rmet(3,3)=metric tensor in real space (length units squared)
!!  rprimd(3,3)=dimensional primitive translations (bohr)
!!  rfpert(3,mpert,3,mpert,3,mpert) = array defining the type of perturbations
!!       that have to be computed
!!       1   ->   element has to be computed explicitely
!!      -1   ->   use symmetry operations to obtain the corresponding element
!!  ucvol=unit cell volume in (whatever length scale units)**3
!!  xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!!  zion(ntypat)=charge on each type of atom (real number)
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!
!! OUTPUT
!!  d3etot_nv(2,3,mpert,3,mpert,3,mpert)= array with the nonvariational
!!              contributions of d3etot
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

subroutine dfptlw_nv(d3etot_nv,dtset,gmet,gprimd,mpert,my_natom,rfpert,rmet,rprimd,ucvol,xred,zion, &
&                 mpi_atmtab,comm_atom ) ! optional arguments (parallelism))
    
 implicit none

!Arguments ------------------------------------
!scalars
 integer , intent(in)  :: mpert,my_natom
 real(dp) :: ucvol
 type(dataset_type),intent(in) :: dtset
 integer,optional,intent(in) :: comm_atom

!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 integer,intent(in) :: rfpert(3,mpert,3,mpert,3,mpert)
 real(dp), intent(out) :: d3etot_nv(2,3,mpert,3,mpert,3,mpert)
 real(dp), intent(in) :: gmet(3,3),rmet(3,3),xred(3,dtset%natom),zion(*)
 real(dp), intent(in) :: gprimd(3,3),rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: alpha,beta,delta,gamma,i1dir,i2dir,i3dir,ii,i1pert,i2pert,i3pert,istr,natom,sumg0 
 real(dp) :: tmpim,tmpre
 character(len=500) :: msg

!arrays
 integer,save :: idx(18)=(/1,1,2,2,3,3,3,2,3,1,2,1,2,3,1,3,1,2/)
 integer :: flg1(3),flg2(3)
 real(dp),allocatable :: dyewdq(:,:,:,:,:,:),dyewdqdq(:,:,:,:,:,:)
 real(dp),allocatable :: dyewdqdq_tII(:,:,:,:,:,:)
 real(dp) :: qphon(3),vec1(3),vec2(3)
 
! *************************************************************************

 DBG_ENTER("COLL")

!Initialiations
 natom=dtset%natom
 d3etot_nv(:,:,:,:,:,:,:)=zero
 
 if (dtset%lw_flexo==1.or.dtset%lw_flexo==3) then

   !1st q-gradient of Ewald contribution to the IFCs
   ABI_MALLOC(dyewdq,(2,3,natom,3,natom,3))
   sumg0=0;qphon(:)=zero
   call dfpt_ewalddq(dyewdq,gmet,my_natom,natom,qphon,rmet,sumg0,dtset%typat,ucvol,xred,zion,&
 & mpi_atmtab=mpi_atmtab,comm_atom=comm_atom)

   i3pert=natom+8
   do i1pert=1,natom
     do i1dir=1,3
       do i2pert=1,natom
         do i2dir=1,3
           do i3dir=1,3
             tmpre=dyewdq(1,i1dir,i1pert,i2dir,i2pert,i3dir)
             tmpim=dyewdq(2,i1dir,i1pert,i2dir,i2pert,i3dir)
             if (abs(tmpre)>=tol8) d3etot_nv(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)= tmpre
             if (abs(tmpim)>=tol8) d3etot_nv(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)= tmpim
           end do
         end do
       end do
     end do
   end do 
   ABI_FREE(dyewdq)

 end if

 if (dtset%lw_flexo==1.or.dtset%lw_flexo==4) then

   !2nd q-gradient of Ewald contribution to the IFCs
   ABI_MALLOC(dyewdqdq,(2,3,natom,3,3,3))
   ABI_MALLOC(dyewdqdq_tII,(2,3,natom,3,3,3))
   sumg0=1;qphon(:)=zero
   call dfpt_ewalddqdq(dyewdqdq,gmet,my_natom,natom,qphon,rmet,sumg0,dtset%typat,ucvol,xred,zion,&
& mpi_atmtab=mpi_atmtab,comm_atom=comm_atom)

   !Convert the indexes labelling the strain perturbation into cartesian coordinates
   !Transform the metric perturbation direction 
   !(treat it as an atomic displacement)
   flg1(:)=1
   do i1pert=1,natom
     do i1dir=1,3
       do gamma=1,3
         do ii=1,2
           do delta=1,3
             do beta=1,3
               vec1(beta)=dyewdqdq(ii,i1dir,i1pert,beta,delta,gamma)
             end do
             call cart39(flg1,flg2,gprimd,i1pert,natom,rprimd,vec1,vec2)
             do beta=1,3
               dyewdqdq(ii,i1dir,i1pert,beta,delta,gamma)=vec2(beta)
             end do
           end do
         end do
       end do
     end do
   end do
               
   !Transform the second q-gradient direction 
   !(treat it as an electric field)
   do i1pert=1,natom
     do i1dir=1,3
       do gamma=1,3
         do ii=1,2
           do beta=1,3
             do delta=1,3
               vec1(delta)=dyewdqdq(ii,i1dir,i1pert,beta,delta,gamma)
             end do
             call cart39(flg1,flg2,gprimd,natom+2,natom,rprimd,vec1,vec2)
             do delta=1,3
               dyewdqdq(ii,i1dir,i1pert,beta,delta,gamma)=vec2(delta)
             end do
           end do

         end do
       end do
     end do
   end do

   !Convert to type-II quantity
   dyewdqdq_tII(:,:,:,:,:,:)=zero
   do i1pert=1,natom
     do alpha=1,3
       do gamma=1,3
         do beta=1,3
           do delta=1,3
             dyewdqdq_tII(:,alpha,i1pert,gamma,beta,delta)= &
           & dyewdqdq(:,alpha,i1pert,beta,delta,gamma) + &
           & dyewdqdq(:,alpha,i1pert,delta,gamma,beta) - &
           & dyewdqdq(:,alpha,i1pert,gamma,beta,delta)
           end do
         end do
       end do
     end do
   end do

   i3pert=natom+8
   do i1pert=1,natom
     do i1dir=1,3
       do i2pert=natom+3,natom+4
         do i2dir=1,3
           istr=(i2pert-natom-3)*3+i2dir
           beta=idx(2*istr-1); delta=idx(2*istr)
           do i3dir=1,3
             tmpre=dyewdqdq_tII(1,i1dir,i1pert,i3dir,beta,delta)
             tmpim=dyewdqdq_tII(2,i1dir,i1pert,i3dir,beta,delta)
             if (abs(tmpre)>=tol8) d3etot_nv(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)= half*tmpre
             if (abs(tmpim)>=tol8) d3etot_nv(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)= half*tmpim
           end do
         end do
       end do
     end do
   end do 

 end if

 !Print results
 if (dtset%prtvol>=10) then
   write(msg,'(3a)') ch10,'LONGWAVE NONVARIATIONAL EWALD D3ETOT: ',ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   do i1pert=1,mpert
     do i1dir=1,3
       do i2pert=1,mpert
         do i2dir=1,3
           do i3pert=1,mpert
             do i3dir=1,3
               if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then 
                 tmpre=d3etot_nv(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
                 tmpim=d3etot_nv(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
                 if (abs(tmpre)>zero.or.abs(tmpim)>zero) then 
                   write(msg,'(3(a,i2,a,i1),2f18.8)') &
           ' perts : ',i1pert,'.',i1dir,' / ',i2pert,'.',i2dir,' / ',i3pert,'.',i3dir,&
                 & tmpre, tmpim
                   call wrtout(std_out,msg,'COLL')
                   call wrtout(ab_out,msg,'COLL')
                 end if
               end if
             end do
           end do
         end do
       end do
     end do
   end do 
   write(msg,'(a)') ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if

 DBG_EXIT("COLL")

end subroutine dfptlw_nv
!!***

!!****f* ABINIT/dfptlw_geom
!! NAME
!!  dfptlw_geom
!!
!! FUNCTION
!!  This routine computes the nonvariational geometric contribution to the 
!!  third-order energy derivative of the flexoelectric force-response tensor.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions at k
!!  cplex: if 1, several magnitudes are REAL, if 2, COMPLEX
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  gsqcut=large sphere cut-off
!!  icg=shift to be applied on the location of data in the array cg
!!  i1dir,i2dir,i3dir=directions of the corresponding perturbations
!!  i1pert,i2pert = type of perturbation that has to be computed
!!  ikpt=number of the k-point
!!  isppol=1 for unpolarized, 2 for spin-polarized
!!  istwf_k=parameter that describes the storage of wfs
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kpt(3)=reduced coordinates of k point
!!  natom= number of atoms in the cell
!!  mkmem =number of k points treated by this node
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  natpert=number of atomic displacement perturbations
!!  nattyp(ntypat)= # atoms of each type.
!!  nband_k=number of bands at this k point for that spin polarization
!!  n2dq= second dimension of d3etot_tgeom_k
!!  nfft=(effective) number of FFT grid points (for this proc)
!!  ngfft(1:18)=integer array with FFT box dimensions and other
!!  npw_k=number of plane waves at this k point
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nylmgr=second dimension of ylmgr_k
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)=1-dimensional phases
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rmet(3,3)=real space metric (bohr**2)
!!  ucvol=unit cell volume in bohr**3.
!!  useylmgr= if 1 use the derivative of spherical harmonics
!!  wtk_k=weight assigned to the k point.
!!  ylm_k(npw_k,psps%mpsang*psps%mpsang*psps%useylm)=real spherical harmonics for the k point
!!  ylmgr_k(npw_k,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr)= k-gradients of real spherical
!!                                                                      harmonics for the k point
!!
!! OUTPUT
!!  d3etot_tgeom_k(2,n2dq)= nonvariational geometric contribution to d3etot for
!     this kpt.
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_dfpt_lw
!!
!! CHILDREN
!!      dfpt_vlocaldq,dfpt_vlocaldqdq,dotprod_g,getgh1dqc,getgh1dqc_setup
!!      init_rf_hamiltonian,mkkpg,rf_hamkq%free,rf_hamkq%load_spin
!!      rf_transgrid_and_pack
!!
!! SOURCE

subroutine dfptlw_geom(atindx,cg,d3etot_tgeom_k,dtset,gs_hamkq,gsqcut,icg, &
       &  i1dir,i2dir,i3dir,i1pert,i2pert,ikpt, &
       &  isppol,istwf_k,kg_k,kpt,mkmem,mpi_enreg,natom,mpw,nattyp,nband_k,n2dq,nfft, &
       &  ngfft,npw_k,nspden,nsppol,nylmgr,occ_k, &
       &  ph1d,psps,rmet,ucvol,useylmgr,wtk_k,ylm_k,ylmgr_k)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,ikpt,isppol,istwf_k
 integer,intent(in) :: i1dir,i1pert,i2dir,i2pert,i3dir
 integer,intent(in) :: natom,mkmem,mpw,nband_k,nfft
 integer,intent(in) :: npw_k,n2dq,nspden,nsppol,nylmgr
 integer,intent(in) :: useylmgr
 real(dp),intent(in) :: gsqcut,ucvol,wtk_k
 type(dataset_type),intent(in) :: dtset
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(MPI_type),intent(in) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps

!arrays
 integer,intent(in) :: atindx(dtset%natom)
 integer,intent(in) :: kg_k(3,npw_k),nattyp(dtset%ntypat),ngfft(18)
 real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*dtset%mband*mkmem*nsppol)
 real(dp),intent(in) :: kpt(3),occ_k(nband_k)
 real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
 real(dp),intent(in) :: rmet(3,3)
 real(dp),intent(in) :: ylm_k(npw_k,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr_k(npw_k,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr)
 real(dp),intent(out) :: d3etot_tgeom_k(2,n2dq)

!Local variables-------------------------------
!scalars
 integer :: beta,delta,gamma,iband,idq,ii,ipw,istr,ka,nkpg,nkpg1,nylmgrpart
 integer :: optlocal,optnl,q1dir,q2dir,tim_getgh1c,useylmgr1
 real(dp) :: doti,dotr
 type(pawfgr_type) :: pawfgr
 type(rf_hamiltonian_type) :: rf_hamkq

!arrays
 integer,save :: idx(18)=(/1,1,2,2,3,3,3,2,3,1,2,1,2,3,1,3,1,2/)
 real(dp) :: q1dirs(2),q2dirs(2)
 real(dp),allocatable :: cwave0i(:,:)
 real(dp),allocatable :: dkinpw(:)
 real(dp),allocatable :: ffnlk(:,:,:,:),ffnl1(:,:,:,:)
 real(dp),allocatable :: gh1dqc(:,:),gh1dqpkc(:,:),gvloc1dqc(:,:),gvnl1dqc(:,:)
 real(dp),allocatable :: kinpw1(:),kpg_k(:,:),kpg1_k(:,:),kpg_pk(:,:),ph3d(:,:,:),ph3d1(:,:,:)
 real(dp),allocatable :: dum_vlocal(:,:,:,:),vlocal1dq(:,:,:,:), dum_vpsp(:)
 real(dp),allocatable :: vpsp1dq(:),part_ylmgr_k(:,:,:)
 type(pawcprj_type),allocatable :: dum_cwaveprj(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

!Definitions
 tim_getgh1c=0
 useylmgr1=1;optlocal=1;optnl=1
 nylmgrpart=3
 nkpg=3
 d3etot_tgeom_k(:,:)=zero

!Allocations
 ABI_MALLOC(cwave0i,(2,npw_k*dtset%nspinor))
 ABI_MALLOC(dum_vpsp,(nfft))
 ABI_MALLOC(dum_vlocal,(ngfft(4),ngfft(5),ngfft(6),gs_hamkq%nvloc))
 ABI_MALLOC(dum_cwaveprj,(0,0))
 ABI_MALLOC(vpsp1dq,(2*nfft))
 ABI_MALLOC(vlocal1dq,(2*ngfft(4),ngfft(5),ngfft(6),gs_hamkq%nvloc))
 ABI_MALLOC(gh1dqc,(2,npw_k*dtset%nspinor))
 ABI_MALLOC(gvloc1dqc,(2,npw_k*dtset%nspinor))
 ABI_MALLOC(gvnl1dqc,(2,npw_k*dtset%nspinor))
 ABI_MALLOC(part_ylmgr_k,(npw_k,3, psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))
 part_ylmgr_k(:,:,:)=ylmgr_k(:,1:3,:)
 ABI_MALLOC(gh1dqpkc,(2,npw_k*dtset%nspinor))
 ABI_MALLOC(kpg_pk,(npw_k,nkpg))

!Generate k+G vectors
 call mkkpg(kg_k,kpg_pk,kpt,nkpg,npw_k)

!Since this is a type-I term, it has to be done for both up and down
!extradiagonal shear strains
 gamma=i3dir
 do idq=1, n2dq
   if (i2pert==natom+3) then                                   
     istr=i2dir
   else                                                        
     istr=idq*3+i2dir                                          
   endif                                                       
   beta=idx(2*istr-1); delta=idx(2*istr)

   !-----------------------------------------------------------------------------------------------
   !  q1-gradient of atomic displacement 1st order hamiltonian:
   !  < u_{i,k}^{(0)} | H^{\tau_{\kappa\alpha}_{\q1dir} \delta_{\beta\q2dir}| u_{i,k}^{(0)} >
   !-----------------------------------------------------------------------------------------------
   q1dirs=(/gamma,delta/)
   q2dirs=(/delta,gamma/)
   do ii=1,2
     q1dir=q1dirs(ii)
     q2dir=q2dirs(ii)

     if (beta==q2dir) then

       !Get q-gradient of first-order local part of the pseudopotential
       call dfpt_vlocaldq(atindx,2,gs_hamkq%gmet,gsqcut,i1dir,i1pert,mpi_enreg, &
       &  psps%mqgrid_vl,dtset%natom,&
       &  nattyp,nfft,ngfft,dtset%ntypat,ngfft(1),ngfft(2),ngfft(3), &
       &  ph1d,q1dir,psps%qgrid_vl,&
       &  dtset%qptn,ucvol,psps%vlspl,vpsp1dq)

       !Set up q-gradient of local potential vlocal1dq with proper dimensioning
       call rf_transgrid_and_pack(isppol,nspden,psps%usepaw,2,nfft,nfft,ngfft,&
       &  gs_hamkq%nvloc,pawfgr,mpi_enreg,dum_vpsp,vpsp1dq,dum_vlocal,vlocal1dq)

       !Initialize rf_hamiltonian (the k-dependent part is prepared in getgh1c_setup)
       call init_rf_hamiltonian(2,gs_hamkq,i1pert,rf_hamkq,&
       & comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab)
       call rf_hamkq%load_spin(isppol,vlocal1=vlocal1dq,with_nonlocal=.true.)

       !Set up the ground-state Hamiltonian, and some parts of the 1st-order Hamiltonian
       call getgh1dqc_setup(gs_hamkq,rf_hamkq,dtset,psps,kpt,kpt,i1dir,i1pert,q1dir, &
     & dtset%natom,rmet,gs_hamkq%gprimd,gs_hamkq%gmet,istwf_k,npw_k,npw_k,nylmgrpart,useylmgr1,kg_k, &
     & ylm_k,kg_k,ylm_k,part_ylmgr_k,nkpg,nkpg1,kpg_k,kpg1_k,dkinpw,kinpw1,ffnlk,ffnl1,ph3d,ph3d1)

       !LOOP OVER BANDS
       do iband=1,nband_k

         if(mpi_enreg%proc_distrb(ikpt,iband,isppol) /= mpi_enreg%me_kpt) cycle

         !Read ket ground-state wavefunctions
         cwave0i(:,:)=cg(:,1+(iband-1)*npw_k*dtset%nspinor+icg:iband*npw_k*dtset%nspinor+icg)

         !Compute < g |H^{\tau_{\kappa\alpha}}_{\q1dir} | u_{i,k}^{(0)} >
         call getgh1dqc(cwave0i,dum_cwaveprj,gh1dqc,gvloc1dqc,gvnl1dqc,gs_hamkq, &
         & i1dir,i1pert,mpi_enreg,optlocal,optnl,q1dir,rf_hamkq)

         !Calculate:
         !<u_{i,k}^{(0)} | H^{\tau_{\kappa\alpha}}_{\q1dir} | u_{i,k}^{(0)} >
         call dotprod_g(dotr,doti,istwf_k,npw_k*dtset%nspinor,2,cwave0i,gh1dqc, &
       & mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

         !Take into account the two pi factor from the term
         !(\hat{p}_{k\beta + \frac{q_{\beta}}{2}}) appearing before the double q-derivation
         !Take also into account here the -i factor and the complex conjugate
         d3etot_tgeom_k(1,idq)=d3etot_tgeom_k(1,idq)-occ_k(iband)*half*doti*two_pi
         d3etot_tgeom_k(2,idq)=d3etot_tgeom_k(2,idq)-occ_k(iband)*half*dotr*two_pi

       end do !iband

       !Clean the rf_hamiltonian
       call rf_hamkq%free()

       !Deallocations
       ABI_FREE(kpg_k)
       ABI_FREE(kpg1_k)
       ABI_FREE(dkinpw)
       ABI_FREE(kinpw1)
       ABI_FREE(ffnlk)
       ABI_FREE(ffnl1)
       ABI_FREE(ph3d)

     end if  

   end do !ii

   !-----------------------------------------------------------------------------------------------
   !  2nd q-gradient of atomic displacement 1st order hamiltonian * momentum operator :
   !  <u_{i,k}^{(0)} | H^{\tau_{\kappa\alpha}}_{\gamma\delta} (k+G)_{\beta} | u_{i,k}^{(0)} >
   !-----------------------------------------------------------------------------------------------

   !Get q-gradient of first-order local part of the pseudopotential
   call dfpt_vlocaldqdq(atindx,2,gs_hamkq%gmet,gsqcut,i1dir,i1pert,mpi_enreg, &
   &  psps%mqgrid_vl,dtset%natom,&
   &  nattyp,nfft,ngfft,dtset%ntypat,ngfft(1),ngfft(2),ngfft(3), &
   &  ph1d,gamma,delta,psps%qgrid_vl,&
   &  dtset%qptn,ucvol,psps%vlspl,vpsp1dq)

   !Set up q-gradient of local potential vlocal1dq with proper dimensioning
   call rf_transgrid_and_pack(isppol,nspden,psps%usepaw,2,nfft,dtset%nfft,dtset%ngfft,&
   &  gs_hamkq%nvloc,pawfgr,mpi_enreg,dum_vpsp,vpsp1dq,dum_vlocal,vlocal1dq)

   !Initialize rf_hamiltonian (the k-dependent part is prepared in getgh1c_setup)
   call init_rf_hamiltonian(2,gs_hamkq,i1pert,rf_hamkq,&
   & comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab)
   call rf_hamkq%load_spin(isppol,vlocal1=vlocal1dq,with_nonlocal=.true.)

   !Set up the ground-state Hamiltonian, and some parts of the 1st-order Hamiltonian
   call getgh1dqc_setup(gs_hamkq,rf_hamkq,dtset,psps,kpt,kpt,i1dir,i1pert,gamma, &
 & dtset%natom,rmet,gs_hamkq%gprimd,gs_hamkq%gmet,istwf_k,npw_k,npw_k,nylmgr,useylmgr1,kg_k, &
 & ylm_k,kg_k,ylm_k,ylmgr_k,nkpg,nkpg1,kpg_k,kpg1_k,dkinpw,kinpw1,ffnlk,ffnl1,ph3d,ph3d1, &
 & qdir2=delta)

   !LOOP OVER BANDS
   do iband=1,nband_k

     if(mpi_enreg%proc_distrb(ikpt,iband,isppol) /= mpi_enreg%me_kpt) cycle

     !Read ket ground-state wavefunctions
     cwave0i(:,:)=cg(:,1+(iband-1)*npw_k*dtset%nspinor+icg:iband*npw_k*dtset%nspinor+icg)

     !Compute < g |H^{\tau_{\kappa\alpha}}_{\gamma\delta} | u_{i,k}^{(0)} >
     call getgh1dqc(cwave0i,dum_cwaveprj,gh1dqc,gvloc1dqc,gvnl1dqc,gs_hamkq, &
     & i1dir,i1pert,mpi_enreg,optlocal,optnl,gamma,rf_hamkq,qdir2=delta)

     !LOOP OVER ONE OF THE STRAIN DIRECTIONS
     do ka=1,3
       do ipw=1,npw_k
         gh1dqpkc(:,ipw)=gh1dqc(:,ipw)*kpg_pk(ipw,ka)
       end do

       call dotprod_g(dotr,doti,istwf_k,npw_k*dtset%nspinor,2,cwave0i,gh1dqpkc, &
     & mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

         !Take into account the two pi factor from the term
         !(\hat{p}_{k\beta + \frac{q_{\beta}}{2}}) appearing before the double q-derivation
       d3etot_tgeom_k(1,idq)=d3etot_tgeom_k(1,idq)-occ_k(iband)*doti*two_pi
       d3etot_tgeom_k(2,idq)=d3etot_tgeom_k(2,idq)-occ_k(iband)*dotr*two_pi

     end do

   end do !iband

   !Clean the rf_hamiltonian
   call rf_hamkq%free()

   !Deallocations
   ABI_FREE(kpg_k)
   ABI_FREE(kpg1_k)
   ABI_FREE(dkinpw)
   ABI_FREE(kinpw1)
   ABI_FREE(ffnlk)
   ABI_FREE(ffnl1)
   ABI_FREE(ph3d)

 end do !idq

!scale by the k-point weight
 d3etot_tgeom_k(:,:)=d3etot_tgeom_k(:,:)*wtk_k

!Deallocations
 ABI_FREE(dum_cwaveprj)
 ABI_FREE(gh1dqc)
 ABI_FREE(gh1dqpkc)
 ABI_FREE(gvloc1dqc)
 ABI_FREE(gvnl1dqc)
 ABI_FREE(vpsp1dq)
 ABI_FREE(vlocal1dq)
 ABI_FREE(dum_vpsp)
 ABI_FREE(dum_vlocal)
 ABI_FREE(kpg_pk)
 ABI_FREE(cwave0i)
 ABI_FREE(part_ylmgr_k)

 DBG_EXIT("COLL")

 end subroutine dfptlw_geom
!!***
end module m_dfptlw_nv
!!***
