!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_wfkfermi
!! NAME
!! dfpt_wfkfermi
!!
!! FUNCTION
!! This routine computes the partial Fermi-level density at a given k-point,
!! and the fixed contribution to the 1st-order Fermi energy (nonlocal and kinetic)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (DRH, XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions
!!  cgq(2,mcgq)=array for planewave coefficients of wavefunctions.
!!  cplex=1 if rhoaug is real, 2 if rhoaug is complex
!!  cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)= wave functions at k
!!              projected with non-local projectors: cprj=<p_i|Cnk>
!!  cprjq(natom,nspinor*mband*mkqmem*nsppol*usecprj)= wave functions at k+q
!!              projected with non-local projectors: cprjq=<p_i|Cnk+q>
!!  dtfil <type(datafiles_type)>=variables related to files
!!  eig0_k(nband_k)=GS eigenvalues at k (hartree)
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  ibg=shift to be applied on the location of data in the array cprj
!!  ibgq=shift to be applied on the location of data in the array cprjq
!!  icg=shift to be applied on the location of data in the array cg
!!  icgq=shift to be applied on the location of data in the array cgq
!!  idir=direction of the current perturbation
!!  ikpt=number of the k-point
!!  ipert=type of the perturbation
!!  isppol=1 for unpolarized, 2 for spin-polarized
!!  kptopt=option for the generation of k points
!!  mband=maximum number of bands
!!  mcgq=second dimension of the cgq array
!!  mcprjq=second dimension of the cprjq array
!!  mkmem =number of k points treated by this node.
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  nband_k=number of bands at this k point for that spin polarization
!!  ncpgr=number of gradients stored in cprj array (cprj=<p_i|Cnk>)
!!  npw_k=number of plane waves at this k point
!!  npw1_k=number of plane waves at this k+q point
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  prtvol=control print volume and debugging output
!!  rf_hamkq <type(gs_hamiltonian_type)>=all data for the 1st-order Hamiltonian at k,q
!!  rhoaug(cplex*n4,n5,n6)= density in electrons/bohr**3,
!!   on the augmented fft grid. (cumulative, so input as well as output)
!!  rocceig(nband_k,nband_k)= (occ_kq(m)-occ_k(n))/(eig0_kq(m)-eig0_k(n)),
!!   if this ratio has been attributed to the band n (second argument), zero otherwise
!!  wtk_k=weight assigned to the k point.
!!
!! OUTPUT
!!  eig1_k(2*nband_k**2)=first-order eigenvalues (hartree)
!!  fe1fixed_k(nband_k)=contribution to 1st-order Fermi energy
!!      from changes of occupation from all bands at this k point.
!!  fe1norm_k(nband_k)=contribution to normalization for above
!!  rhoaug(cplex*n4,n5,n6)= Fermi-level density in electrons/bohr**3,
!!   on the augmented fft grid. (cumulative, so input as well as output).
!!  ==== if (gs_hamkq%usepaw==1) ====
!!    pawrhoijfermi(natom) <type(pawrhoij_type)>= paw rhoij occupancies
!!       at Fermi level (cumulative, so input as well as output)
!!
!! PARENTS
!!      dfpt_rhofermi
!!
!! CHILDREN
!!      dfpt_accrho,dotprod_g,getgh1c,pawcprj_alloc,pawcprj_axpby,pawcprj_copy
!!      pawcprj_free,pawcprj_get,status,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_wfkfermi(cg,cgq,cplex,cprj,cprjq,&
&          dtfil,eig0_k,eig1_k,fe1fixed_k,fe1norm_k,gs_hamkq,&
&          ibg,ibgq,icg,icgq,idir,ikpt,ipert,isppol,&
&          kptopt,mband,mcgq,mcprjq,mkmem,mpi_enreg,mpw,nband_k,ncpgr,&
&          npw_k,npw1_k,nspinor,nsppol,occ_k,pawrhoijfermi,prtvol,&
&          rf_hamkq,rhoaug,rocceig,wtk_k)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_cgtools

 use m_pawrhoij,    only : pawrhoij_type
 use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_get, pawcprj_copy, pawcprj_axpby, pawcprj_free
 use m_hamiltonian, only : gs_hamiltonian_type,rf_hamiltonian_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_wfkfermi'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_66_wfs
 use interfaces_72_response, except_this_one => dfpt_wfkfermi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,ibg,ibgq,icg,icgq,idir,ikpt
 integer,intent(in) :: ipert,isppol,kptopt,mband,mcgq,mcprjq,mkmem,mpw,ncpgr
 integer,intent(in) :: npw1_k,nspinor,nsppol,prtvol
 integer,intent(inout) :: nband_k,npw_k
 real(dp),intent(in) :: wtk_k
 type(MPI_type),intent(in) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(rf_hamiltonian_type),intent(inout) :: rf_hamkq
!arrays
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),cgq(2,mcgq)
 real(dp),intent(in) :: eig0_k(nband_k),occ_k(nband_k),rocceig(nband_k,nband_k)
 real(dp),intent(inout) :: rhoaug(cplex*gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,gs_hamkq%nvloc)
 real(dp),intent(inout) :: eig1_k(2*nband_k**2)
 real(dp),intent(out) :: fe1fixed_k(nband_k)
 real(dp),intent(out) :: fe1norm_k(nband_k)
 type(pawcprj_type),intent(in) :: cprj(gs_hamkq%natom,nspinor*mband*mkmem*nsppol*gs_hamkq%usecprj)
 type(pawcprj_type),intent(in) :: cprjq(gs_hamkq%natom,mcprjq)
 type(pawrhoij_type),intent(inout) :: pawrhoijfermi(gs_hamkq%natom*gs_hamkq%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=18
 integer :: berryopt,counter,iband,iexit,ii,indx,iorder_cprj
 integer :: ipw,me,nkpt_max,optlocal,optnl,opt_accrho,opt_corr
 integer :: opt_gvnl1,sij_opt,tim_fourwf,tim_getgh1c,usevnl
 real(dp) :: dotr,lambda,wtband
 character(len=500) :: msg
!arrays
 real(dp) :: dum_grad_berry(1,1),dum_gvnl1(1,1),dum_gs1(1,1),tsec(2)
 real(dp),allocatable :: cwave0(:,:),cwaveq(:,:),gh1(:,:)
 type(pawcprj_type),allocatable :: cwaveprj0(:,:),cwaveprjq(:,:),cwaveprj_tmp(:,:)

! *********************************************************************

 DBG_ENTER('COLL')

!Check arguments validity
 if (ipert>gs_hamkq%natom.and.ipert/=gs_hamkq%natom+3.and.ipert/=gs_hamkq%natom+4.and.ipert/=gs_hamkq%natom+5) then !SPr rfmagn deb
   msg='  wrong ipert argument !'
   MSG_BUG(msg)
 end if
 if (cplex/=1) then
   MSG_BUG('wrong cplex/=1 argument !')
 end if

!Debugging statements
 call status(0,dtfil%filstat,iexit,level,'enter dfpt_wfkfermi')
 if(prtvol==-level)then
   write(msg,'(80a,a,a)') ('=',ii=1,80),ch10,'dfpt_wfkfermi : enter'
   call wrtout(std_out,msg,'PERS')
 end if
 nkpt_max=50;if(xmpi_paral==1)nkpt_max=-1

 if(prtvol>2 .or. ikpt<=nkpt_max)then
   write(msg, '(a,a,i5,2x,a,3f9.5,2x,a)' ) ch10,&
&   ' Non-SCF iterations; k pt #',ikpt,'k=',gs_hamkq%kpt_k(:),' band residuals:'
   call wrtout(std_out,msg,'PERS')
 end if

!Retrieve parallelism data
 me=mpi_enreg%me_kpt
!Initializations and allocations

 ABI_ALLOCATE(gh1,(2,npw1_k*nspinor))
 ABI_ALLOCATE(cwave0,(2,npw_k*nspinor))
 ABI_ALLOCATE(cwaveq,(2,npw1_k*nspinor))
 iorder_cprj=0 ; eig1_k(:)=zero
 if (gs_hamkq%usepaw==1.and.gs_hamkq%usecprj==1) then
   ABI_DATATYPE_ALLOCATE(cwaveprj0,(gs_hamkq%natom,nspinor))
   ABI_DATATYPE_ALLOCATE(cwaveprjq,(gs_hamkq%natom,nspinor))
   call pawcprj_alloc(cwaveprj0,1,gs_hamkq%dimcprj)
   call pawcprj_alloc(cwaveprjq,0,gs_hamkq%dimcprj)
 else
   ABI_DATATYPE_ALLOCATE(cwaveprj0,(0,0))
   ABI_DATATYPE_ALLOCATE(cwaveprjq,(0,0))
 end if
!Arguments of getgh1c routine (want only (NL+kin) frozen H(1))
 berryopt=0;usevnl=0;sij_opt=-gs_hamkq%usepaw;tim_getgh1c=3
 optlocal=0;optnl=1;opt_gvnl1=0
 if(ipert==gs_hamkq%natom+5) optnl=0;    ! no 1st order NL in H(1), also no kin, but this will be taken into account later
!if(ipert==gs_hamkq%natom+5) optlocal=0; ! 1st order LOCAL potential present

!Arguments of the dfpt_accrho routine
 tim_fourwf=5 ; opt_accrho=1 ; opt_corr=0
!Null potentially unassigned output variables
 fe1fixed_k(:)=zero; fe1norm_k(:)=zero

!Read the npw and kg records of wf files
 call status(0,dtfil%filstat,iexit,level,'before WffRead')

 call timab(139,1,tsec)


!Loop over bands
 do iband=1,nband_k
   counter=100*iband+1

!  Skip bands not treated by current proc
   if(mpi_enreg%proc_distrb(ikpt, iband,isppol)/=me) cycle

   if(prtvol>=10)then
     call status(counter,dtfil%filstat,iexit,level,'loop iband    ')
   end if

!  Select occupied bands
   if(abs(occ_k(iband))>tol8.and.abs(rocceig(iband,iband))>tol8)then

     wtband=rocceig(iband,iband)/occ_k(iband)
!    Get ground-state wavefunctions at k
     do ipw=1,npw_k*nspinor
       cwave0(1,ipw)=cg(1,ipw+(iband-1)*npw_k*nspinor+icg)
       cwave0(2,ipw)=cg(2,ipw+(iband-1)*npw_k*nspinor+icg)
     end do

     if (gs_hamkq%usepaw==1.and.gs_hamkq%usecprj==1) then
!      Read PAW ground state projected WF (cprj)
       call pawcprj_get(gs_hamkq%atindx1,cwaveprj0,cprj,gs_hamkq%natom,iband,ibg,ikpt,iorder_cprj,&
&       isppol,mband,mkmem,gs_hamkq%natom,1,nband_k,nspinor,nsppol,dtfil%unpaw,&
&       mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb,&
&       icpgr=idir,ncpgr=ncpgr)
     end if

!    Read ground-state wavefunctions at k+q
     indx=npw1_k*nspinor*(iband-1)+icgq
     cwaveq(:,1:npw_k*nspinor)=wtband*cgq(:,1+indx:npw_k*nspinor+indx)
     if (gs_hamkq%usepaw==1.and.gs_hamkq%usecprj==1) then
!      Read PAW ground state projected WF (cprj)
       indx=nspinor*(iband-1)+ibgq
       call pawcprj_copy(cprjq(:,1+indx:nspinor+indx),cwaveprjq)
       call pawcprj_axpby(zero,wtband,cwaveprj_tmp,cwaveprjq)
     end if

     if(prtvol>=10)then
       call status(0,dtfil%filstat,iexit,level,'after wf read ')
     end if


!    Apply H^(1)-Esp.S^(1) to Psi^(0) (H(^1)=only (NL+kin) frozen part)
     lambda=eig0_k(iband)
     call getgh1c(berryopt,cwave0,cwaveprj0,gh1,dum_grad_berry,dum_gs1,gs_hamkq,dum_gvnl1,&
&     idir,ipert,lambda,mpi_enreg,optlocal,optnl,opt_gvnl1,rf_hamkq,sij_opt,&
&     tim_getgh1c,usevnl)
!    Compute Eig1=<Psi^(0)|H^(1)-Eps.S^(1)|Psi(0)>
     call dotprod_g(dotr,lambda,gs_hamkq%istwf_k,npw_k*nspinor,1,cwave0,gh1,mpi_enreg%me_g0, &
&     mpi_enreg%comm_spinorfft)
     indx=2*iband-1+(iband-1)*2*nband_k
     eig1_k(indx)=dotr
!    Compute the fixed contribution to the 1st-order Fermi energy
     fe1fixed_k(iband)=two*wtband*eig1_k(indx)
     fe1norm_k(iband) =two*wtband
     
!    Accumulate contribution to density and PAW occupation matrix
  
     call dfpt_accrho(counter,cplex,cwave0,cwaveq,cwaveq,cwaveprj0,cwaveprjq,dotr,&
&     dtfil%filstat,gs_hamkq,iband,0,0,isppol,kptopt,mpi_enreg,gs_hamkq%natom,nband_k,ncpgr,&
&     npw_k,npw1_k,nspinor,occ_k,opt_accrho,pawrhoijfermi,prtvol,rhoaug,tim_fourwf,&
&     opt_corr,wtk_k)

   end if ! End of non-zero occupation and rocceig

 end do ! End loop over bands

 call timab(139,2,tsec)
 call timab(130,1,tsec)
 call status(0,dtfil%filstat,iexit,level,'after loops   ')

 ABI_DEALLOCATE(cwave0)
 ABI_DEALLOCATE(cwaveq)
 ABI_DEALLOCATE(gh1)
 if (gs_hamkq%usepaw==1.and.gs_hamkq%usecprj==1) then
   call pawcprj_free(cwaveprj0)
   call pawcprj_free(cwaveprjq)
 end if
 ABI_DATATYPE_DEALLOCATE(cwaveprj0)
 ABI_DATATYPE_DEALLOCATE(cwaveprjq)

!Structured debugging : if prtvol=-level, stop here.
 if(prtvol==-level)then
   write(msg,'(a,a1,a,i2,a)')' fermie3 : exit prtvol=-',level,', debugging mode => stop '
   MSG_ERROR(msg)
 end if

 call status(0,dtfil%filstat,iexit,level,'exit dfpt_wfkfermi')

 call timab(130,2,tsec)

 DBG_EXIT('COLL')

end subroutine dfpt_wfkfermi
!!***
