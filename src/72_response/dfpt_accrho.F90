!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_accrho
!!
!! NAME
!! dfpt_accrho
!!
!! FUNCTION
!! Response function calculation only:
!!  Accumulate contribution to first-order density due do current (k,band)
!!  Also accumulate zero-order potential part of the 2nd-order total energy (if needed)
!!
!! COPYRIGHT
!! Copyright (C) 2009-2016 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  counter=counter for status file
!!  cplex=1 if 1st-order density is real, 2 if 1st-order density is complex
!!  cwave0(2,npw*nspinor)=GS wavefunction at k, in reciprocal space
!!  cwave1(2,npw1*nspinor)=1st-order wavefunction at k,q, in reciprocal space
!!  cwavef(2,npw1*nspinor)=1st-order wavefunction at k,q, in reciprocal space, without correction due to occupation change
!!  cwaveprj0(natom,nspinor*usecprj)= GS wave function at k projected with nl projectors
!!  cwaveprj1(natom,nspinor*usecprj)= 1st-order wave function at k,q projected with nl projectors
!!  filstat=name of the status file
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  iband=index of current band
!!  idir=direction of the current perturbation
!!  ipert=type of the perturbation
!!  isppol=1 index of current spin component
!!  kptopt=option for the generation of k points
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  nband_k=number of bands at this k point for that spin polarization
!!  ncpgr=number of gradients stored in cprj array (cprj=<p_i|Cnk>)
!!  npw_k=number of planewaves in basis sphere at k
!!  npw1_k=number of planewaves in basis sphere at k+q
!!  nspinor=number of spinorial components of the wavefunctions
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  option= 1: accumulate 1st-order density,
!!          2: accumulate 0-order potential part of the 2nd-order total energy
!!          3: accumulate both
!!  prtvol=control print volume and debugging output
!!  tim_fourwf= timing code for fourwf (5 from dfpt_vtowfk, 18 from dfpt_nstwf)
!!  wf_corrected=flag put to 1 if cwave1 is different from cwavef (if there is a contribution from occ. change)
!!  wtk_k=weight assigned to the k point.
!!
!! OUTPUT
!!  ====== if option=2 or option=3 =====
!!  eloc0_k=zero-order local contribution to 2nd-order total energy for current band and k
!!
!! SIDE EFFECTS
!!  ====== if option=1 or option=3 =====
!!    rhoaug1(cplex*n4,n5,n6,nvloc)= density in electrons/bohr**3,
!!    ==== if gs_hamkq%usepaw=1 =====
!!    pawrhoij1(natom) <type(pawrhoij_type)>= 1st-order paw rhoij occupancies and related data
!!                                            (cumulative, so input as well as output)
!!
!! NOTES
!!  In this part of the treatment of one band, one has to
!!  perform Fourier transforms, and to treat separately the
!!  two spinorial components of the wavefunction.
!!  Was part of dfpt_vtowfk before.
!!
!! PARENTS
!!      dfpt_nstpaw,dfpt_vtowfk,dfpt_wfkfermi
!!
!! CHILDREN
!!      fourwf,get_my_atmtab,getcprj,pawaccrhoij,pawcprj_alloc,pawcprj_free
!!      status
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dfpt_accrho(counter,cplex,cwave0,cwave1,cwavef,cwaveprj0,cwaveprj1,&
&                  eloc0_k,filstat,gs_hamkq,iband,idir,ipert,isppol,kptopt,&
&                  mpi_enreg,natom,nband_k,ncpgr,npw_k,npw1_k,nspinor,occ_k,&
&                  option,pawrhoij1,prtvol,rhoaug1,tim_fourwf,wf_corrected,&
&                  wtk_k,comm_atom,mpi_atmtab)


 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_hamiltonian,   only : gs_hamiltonian_type
 use m_pawrhoij,      only : pawrhoij_type
 use m_pawcprj,       only : pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_paral_atom,    only : get_my_atmtab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_accrho'
 use interfaces_32_util
 use interfaces_53_ffts
 use interfaces_65_paw
 use interfaces_66_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: counter,cplex,iband,idir,ipert,isppol,kptopt,natom,nband_k
 integer,intent(in) :: ncpgr,npw_k,npw1_k,nspinor,option,prtvol,tim_fourwf,wf_corrected
 integer,optional,intent(in) :: comm_atom
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: wtk_k
 real(dp),intent(out) :: eloc0_k
 character(len=fnlen),intent(in) :: filstat
 type(gs_hamiltonian_type),intent(inout),target :: gs_hamkq
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in),target :: cwave0(2,npw_k*nspinor),cwave1(2,npw1_k*nspinor),cwavef(2,npw1_k*nspinor)
 real(dp),intent(in) :: occ_k(nband_k)
 real(dp),intent(inout) :: rhoaug1(cplex*gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,gs_hamkq%nvloc)
 type(pawcprj_type),intent(in) :: cwaveprj0(natom,nspinor*gs_hamkq%usecprj)
 type(pawcprj_type),intent(in) :: cwaveprj1(natom,nspinor*gs_hamkq%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij1(:)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=14
 integer :: choice,cplex_cprj,i1,i2,i3,iexit,ispinor,my_comm_atom,my_natom,n1,n2,n3,option_rhoij
 logical :: my_atmtab_allocated,paral_atom
 logical :: usetimerev
 real(dp) :: im0,im1,re0,re1,valuer,diag,offdiag,weight
 real(dp) :: im0_up,im1_up,re0_up,re1_up,im0_down,im1_down,re0_down,re1_down
!arrays
 integer,pointer :: my_atmtab(:)
 real(dp) :: dummy(2,1)
 real(dp),allocatable :: rhoaug(:,:,:,:),wfraug(:,:,:,:),wfraug1(:,:,:,:)
 real(dp),allocatable :: wfraug1_up(:,:,:,:),wfraug1_down(:,:,:,:)
 real(dp),allocatable :: wfraug_up(:,:,:,:),wfraug_down(:,:,:,:)
 real(dp),pointer :: cwavef_sp(:,:),cwavef_up(:,:),cwavef_down(:,:)
 real(dp),pointer :: cwave0_up(:,:),cwave0_down(:,:),cwave1_up(:,:),cwave1_down(:,:)
 real(dp),pointer :: vlocal(:,:,:,:)=>null()
 type(pawcprj_type),allocatable :: cwaveprj_tmp(:,:)


! *********************************************************************

 DBG_ENTER("COLL")

 if (option/=1.and.option/=2.and.option/=3) return

!Initializations
 ABI_ALLOCATE(rhoaug,(gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,gs_hamkq%nvloc))
 ABI_ALLOCATE(wfraug1,(2,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))
 n1=gs_hamkq%ngfft(1);n2=gs_hamkq%ngfft(2);n3=gs_hamkq%ngfft(3)
 if (option==2.or.option==3) eloc0_k=zero
 if (option==2.or.option==3) vlocal => gs_hamkq%vlocal

!Loop on spinorial components
! TODO : double loop on spinors for full rhoaug1 matrix if nspden =4
 if (gs_hamkq%nvloc/=4) then  ! see later EB FR
   do ispinor=1,nspinor

     if (prtvol>=10) then
       call status(counter,filstat,iexit,level,'density update')
     end if

!  Part devoted to the accumulation of the 0-order potential part of the 2nd-order total energy
!  --------------------------------------------------------------------------------------------

!  Fourier transform of cwavef. Here, rhoaug is a dummy variable.
     if (wf_corrected==0.or.option==2.or.option==3) then
       if (ispinor==1) then
         cwavef_sp => cwavef(:,1:npw1_k)
       else
         cwavef_sp => cwavef(:,1+npw1_k:2*npw1_k)
       end if
       call fourwf(cplex,rhoaug,cwavef_sp,dummy,wfraug1,gs_hamkq%gbound_kp,gs_hamkq%gbound_kp,&
&       gs_hamkq%istwf_k,gs_hamkq%kg_kp,gs_hamkq%kg_kp,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&       gs_hamkq%npw_kp,1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,mpi_enreg%paral_kgb,tim_fourwf,&
&       weight,weight,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
       nullify(cwavef_sp)
     end if

!  Compute contribution of this band to zero-order potential part of the 2nd-order total energy
!  NB: this is spinor diagonal
     if (option==2.or.option==3) then
       valuer=zero
       do i3=1,n3
         do i2=1,n2
           do i1=1,n1
             valuer=valuer+vlocal(i1,i2,i3,1)*(wfraug1(1,i1,i2,i3)**2+wfraug1(2,i1,i2,i3)**2)
           end do
         end do
       end do

!    Local potential energy of this band
       eloc0_k=eloc0_k+two*valuer/dble(gs_hamkq%nfft)
     end if ! option

!  Part devoted to the accumulation of the 1st-order density
!  ---------------------------------------------------------
     if (option==1.or.option==3) then

!    Compute 1st-order WF in real space
!    One needs the Fourier transform of cwave1. However, only the one of
!    cwavef is available. If cwavef and cwave1 differs, this Fourier
!    transform must be computed. In both case the result is in wfraug1.
       if (wf_corrected==1) then
         if (ispinor==1) then
           cwavef_sp => cwave1(:,1:npw1_k)
         else
           cwavef_sp => cwave1(:,1+npw1_k:2*npw1_k)
         end if
         call fourwf(cplex,rhoaug,cwavef_sp,dummy,wfraug1,gs_hamkq%gbound_kp,gs_hamkq%gbound_kp,&
&         gs_hamkq%istwf_k,gs_hamkq%kg_kp,gs_hamkq%kg_kp,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&         gs_hamkq%npw_kp,1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,mpi_enreg%paral_kgb,tim_fourwf,&
&         weight,weight,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
         nullify(cwavef_sp)
       end if

!    Compute 0-order WF in real space
! TODO: add loop over ispinor_prime here
       ABI_ALLOCATE(wfraug,(2,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))
       if (ispinor==1) then
         cwavef_sp => cwave0(:,1:npw_k)
       else
         cwavef_sp => cwave0(:,1+npw_k:2*npw_k)
       end if
       call fourwf(1,rhoaug,cwavef_sp,dummy,wfraug,gs_hamkq%gbound_k,gs_hamkq%gbound_k,&
&       gs_hamkq%istwf_k,gs_hamkq%kg_k,gs_hamkq%kg_k,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&       gs_hamkq%npw_k,1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,mpi_enreg%paral_kgb,tim_fourwf,&
&       weight,weight,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
       nullify(cwavef_sp)

!    The factor 2 is not the spin factor (see Eq.44 of PRB55,10337 (1997))
       weight=two*occ_k(iband)*wtk_k/gs_hamkq%ucvol

!    Accumulate 1st-order density
       if (cplex==2) then
         do i3=1,n3
           do i2=1,n2
             do i1=1,n1
               re0=wfraug(1,i1,i2,i3)  ; im0=wfraug(2,i1,i2,i3)
               re1=wfraug1(1,i1,i2,i3) ; im1=wfraug1(2,i1,i2,i3)
! TODO: check which terms (ispinor ispinorp) enter a given element of rhoaug1
               rhoaug1(2*i1-1,i2,i3,1)=rhoaug1(2*i1-1,i2,i3,1)+weight*(re0*re1+im0*im1)
               rhoaug1(2*i1  ,i2,i3,1)=rhoaug1(2*i1  ,i2,i3,1)+weight*(re0*im1-im0*re1)
             end do
           end do
         end do
       else
         do i3=1,n3
           do i2=1,n2
             do i1=1,n1
               rhoaug1(i1,i2,i3,1)=rhoaug1(i1,i2,i3,1) &
&               +weight*(wfraug(1,i1,i2,i3)*wfraug1(1,i1,i2,i3) &
&               +wfraug(2,i1,i2,i3)*wfraug1(2,i1,i2,i3))
             end do
           end do
         end do
       end if
       ABI_DEALLOCATE(wfraug)
     end if ! option

   end do ! Loop on spinorial components if nspden=1 or 2
   ABI_DEALLOCATE(rhoaug)
   ABI_DEALLOCATE(wfraug1)

 else ! nvloc = 4
! The same lines of code are in 72_response/dfpt_mkrho.F90
! TODO merge these lines in a single routine??!!
   if (prtvol>=10) then
     call status(counter,filstat,iexit,level,'density update')
   end if
   ABI_ALLOCATE(wfraug1_up,(2,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))
   ABI_ALLOCATE(wfraug1_down,(2,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))

!  Part devoted to the accumulation of the 0-order potential part of the 2nd-order total energy
!  --------------------------------------------------------------------------------------------

!  Fourier transform of cwavef. Here, rhoaug is a dummy variable.
   if (wf_corrected==0.or.option==2.or.option==3) then
     cwavef_up => cwavef(:,1:npw1_k) ! wfs up spin-polarized
     call fourwf(cplex,rhoaug(:,:,:,1),cwavef_up,dummy,wfraug1_up,gs_hamkq%gbound_kp,gs_hamkq%gbound_kp,&
&     gs_hamkq%istwf_k,gs_hamkq%kg_kp,gs_hamkq%kg_kp,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&     gs_hamkq%npw_kp,1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,mpi_enreg%paral_kgb,tim_fourwf,&
&     weight,weight,use_gpu_cuda=gs_hamkq%use_gpu_cuda)

     cwavef_down => cwavef(:,1+npw1_k:2*npw1_k) ! wfs down spin-polarized
     call fourwf(cplex,rhoaug(:,:,:,1),cwavef_down,dummy,wfraug1_down,gs_hamkq%gbound_kp,gs_hamkq%gbound_kp,&
&     gs_hamkq%istwf_k,gs_hamkq%kg_kp,gs_hamkq%kg_kp,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&     gs_hamkq%npw_kp,1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,mpi_enreg%paral_kgb,tim_fourwf,&
&     weight,weight,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
     nullify(cwavef_up)
     nullify(cwavef_down)
   end if
   if (option==2.or.option==3) then
     valuer=zero
     diag=zero
     offdiag=zero
   ! EB FR 2nd term in Eq. 91 PRB52,1096 for non-collinear magnetism
     do i3=1,n3
       do i2=1,n2
         do i1=1,n1
           diag=vlocal(i1,i2,i3,1)*(wfraug1_up(1,i1,i2,i3)**2+wfraug1_up(2,i1,i2,i3)**2)&
&           +vlocal(i1,i2,i3,2)*(wfraug1_down(1,i1,i2,i3)**2+wfraug1_down(2,i1,i2,i3)**2)
           offdiag=(2*vlocal(i1,i2,i3,3)*((wfraug1_up(1,i1,i2,i3)*wfraug1_down(1,i1,i2,i3))+&
&           (wfraug1_up(2,i1,i2,i3)*wfraug1_down(2,i1,i2,i3))))+&
&           (2*vlocal(i1,i2,i3,4)*((-wfraug1_down(2,i1,i2,i3)*wfraug1_up(1,i1,i2,i3))+&
&           (wfraug1_down(1,i1,i2,i3)*wfraug1_up(2,i1,i2,i3))))
           valuer=valuer+diag+offdiag
         end do
       end do
     end do
     ABI_DEALLOCATE(wfraug1_up)
     ABI_DEALLOCATE(wfraug1_down)
     
!    Local potential energy of this band
     eloc0_k=eloc0_k+two*valuer/dble(gs_hamkq%nfft)
   end if ! option
!  Part devoted to the accumulation of the 1st-order density
!  ---------------------------------------------------------
   if (option==1.or.option==3) then
!        Build the four components of rho. We use only norm quantities and, so fourwf.
     ABI_ALLOCATE(wfraug_up,(2,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))
     ABI_ALLOCATE(wfraug_down,(2,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))    
     ABI_ALLOCATE(wfraug1_up,(2,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))
     ABI_ALLOCATE(wfraug1_down,(2,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))
! EB FR build spinorial wavefunctions
! zero order
     cwave0_up => cwave0(:,1:npw_k)
     cwave0_down => cwave0(:,1+npw_k:2*npw_k)
! first order
     cwave1_up => cwave1(:,1:npw_k)
     cwave1_down => cwave1(:,1+npw_k:2*npw_k)
   end if ! option
! The factor 2 is not the spin factor (see Eq.44 of PRB55,10337 (1997))
   weight=two*occ_k(iband)*wtk_k/gs_hamkq%ucvol
!density components
!GS wfk Fourrier Tranform
! EB FR in the fourwf calls rhoaug(:,:,:,2) is a dummy argument
   call fourwf(1,rhoaug(:,:,:,2),cwave0_up,dummy,wfraug_up,gs_hamkq%gbound_k,gs_hamkq%gbound_k,&
&   gs_hamkq%istwf_k,gs_hamkq%kg_k,gs_hamkq%kg_k,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&   gs_hamkq%npw_k,1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,mpi_enreg%paral_kgb,tim_fourwf,&
&   weight,weight,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
   call fourwf(1,rhoaug(:,:,:,2),cwave0_down,dummy,wfraug_down,gs_hamkq%gbound_k,gs_hamkq%gbound_k,&
&   gs_hamkq%istwf_k,gs_hamkq%kg_k,gs_hamkq%kg_k,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&   gs_hamkq%npw_k,1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,mpi_enreg%paral_kgb,tim_fourwf,&
&   weight,weight,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
!1st order wfk Fourrier Transform
   call fourwf(1,rhoaug1(:,:,:,2),cwave1_up,dummy,wfraug1_up,gs_hamkq%gbound_k,gs_hamkq%gbound_k,&
&   gs_hamkq%istwf_k,gs_hamkq%kg_k,gs_hamkq%kg_k,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&   gs_hamkq%npw_k,1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,mpi_enreg%paral_kgb,tim_fourwf,&
&   weight,weight,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
   call fourwf(1,rhoaug1(:,:,:,2),cwave1_down,dummy,wfraug1_down,gs_hamkq%gbound_k,gs_hamkq%gbound_k,&
&   gs_hamkq%istwf_k,gs_hamkq%kg_k,gs_hamkq%kg_k,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&   gs_hamkq%npw_k,1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,mpi_enreg%paral_kgb,tim_fourwf,&
&   weight,weight,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
! Accumulate 1st-order density (x component)
   if (cplex==2) then
     re0_up=zero;im0_up=zero;re1_up=zero;im1_up=zero;re0_down=zero;im0_down=zero
     re1_down=zero;im1_down=zero
     do i3=1,n3
       do i2=1,n2
         do i1=1,n1
           re0_up=wfraug_up(1,i1,i2,i3)  ;     im0_up=wfraug_up(2,i1,i2,i3)
           re1_up=wfraug1_up(1,i1,i2,i3) ;     im1_up=wfraug1_up(2,i1,i2,i3)
           re0_down=wfraug_down(1,i1,i2,i3)  ; im0_down=wfraug_down(2,i1,i2,i3)
           re1_down=wfraug1_down(1,i1,i2,i3) ; im1_down=wfraug1_down(2,i1,i2,i3)
           rhoaug1(2*i1-1,i2,i3,1)=rhoaug1(2*i1-1,i2,i3,1)+weight*((re0_up*re1_up+im0_up*im1_up)+&
&           (re0_down*re1_down+im0_down*im1_down)) ! trace
           rhoaug1(2*i1  ,i2,i3,1)=zero ! imag part of rho at k
           rhoaug1(2*i1-1,i2,i3,2)=rhoaug1(2*i1-1,i2,i3,2)+weight*((re0_up*re1_up+im0_up*im1_up)-&
&           (re0_down*re1_down+im0_down*im1_down)) ! m_z
           rhoaug1(2*i1  ,i2,i3,2)=zero ! imag part of rho at k
           rhoaug1(2*i1-1,i2,i3,3)=rhoaug1(2*i1-1,i2,i3,3)+weight*((re0_up*re1_down+im0_up*im1_down)+&
&           re0_down*re1_up+im0_down*im1_up) ! m_x
           rhoaug1(2*i1  ,i2,i3,3)=zero ! imag part of rho at k
           rhoaug1(2*i1-1,i2,i3,4)=rhoaug1(2*i1-1,i2,i3,4)+weight*((-re1_up*im0_down+im1_up*re0_down)&
&           -re0_up*im1_down+im0_up*re1_down) ! m_y
           rhoaug1(2*i1  ,i2,i3,4)=zero ! imag part of rho at k
         end do
       end do
     end do
   else !cplex
     re0_up=zero;im0_up=zero;re1_up=zero;im1_up=zero;re0_down=zero;im0_down=zero
     re1_down=zero;im1_down=zero
     do i3=1,n3
       do i2=1,n2
         do i1=1,n1
           re0_up=wfraug_up(1,i1,i2,i3)  ;     im0_up=wfraug_up(2,i1,i2,i3)
           re1_up=wfraug1_up(1,i1,i2,i3) ;     im1_up=wfraug1_up(2,i1,i2,i3)
           re0_down=wfraug_down(1,i1,i2,i3)  ; im0_down=wfraug_down(2,i1,i2,i3)
           re1_down=wfraug1_down(1,i1,i2,i3) ; im1_down=wfraug1_down(2,i1,i2,i3)
           rhoaug1(i1,i2,i3,1)=rhoaug1(i1,i2,i3,1)+weight*((re0_up*re1_up+im0_up*im1_up)+&
&           (re0_down*re1_down+im0_down*im1_down)) ! n
           rhoaug1(i1,i2,i3,2)=rhoaug1(i1,i2,i3,2)+weight*((re0_up*re1_up+im0_up*im1_up)-&
&           (re0_down*re1_down+im0_down*im1_down)) ! m_z
           rhoaug1(i1,i2,i3,3)=rhoaug1(i1,i2,i3,3)+weight*((re0_up*re1_down+im0_up*im1_down)+&
&           re0_up*re1_down+im0_down*im1_up)     ! m_x
           rhoaug1(i1,i2,i3,4)=rhoaug1(i1,i2,i3,4)+weight*((-re1_up*im0_down+im1_up*re0_down)+&
&           (-re0_up*im1_down+im0_up*re1_down)) ! m_y
         end do
       end do
     end do
   end if !cplex
   ABI_DEALLOCATE(wfraug_up)
   ABI_DEALLOCATE(wfraug_down)
   ABI_DEALLOCATE(wfraug1_up)
   ABI_DEALLOCATE(wfraug1_down)

 end if ! nvloc /= 4

!Part devoted to the accumulation of the 1st-order occupation matrix in PAW case
! TODO: parse for more nspden 4 dependencies on spinors
! EB FR CHECK: to be modified for non-collinear?????
!-------------------------------------------------------------------------------

 if ((option==1.or.option==3).and.gs_hamkq%usepaw==1) then

!  Set up parallelism over atoms
   my_natom=natom; if(gs_hamkq%usepaw==1) my_natom=size(pawrhoij1)
   paral_atom=(present(comm_atom).and.(my_natom/=natom))
   my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
   nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
   call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

   cplex_cprj=2;if (gs_hamkq%istwf_k>1) cplex_cprj=1
   option_rhoij=2;usetimerev=(kptopt>0.and.kptopt<3)

   if (gs_hamkq%usecprj==1) then
     call pawaccrhoij(gs_hamkq%atindx,cplex_cprj,cwaveprj0,cwaveprj1,ipert,isppol,&
&     my_natom,natom,nspinor,occ_k(iband),option_rhoij,pawrhoij1,usetimerev,wtk_k,&
&     comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
   else
     ABI_DATATYPE_ALLOCATE(cwaveprj_tmp,(natom,nspinor))
     call pawcprj_alloc(cwaveprj_tmp,ncpgr,gs_hamkq%dimcprj)
     choice=2
     call getcprj(choice,0,cwave0,cwaveprj_tmp,&
&     gs_hamkq%ffnl_k,idir,gs_hamkq%indlmn,gs_hamkq%istwf_k,&
&     gs_hamkq%kg_k,gs_hamkq%kpg_k,gs_hamkq%kpt_k,gs_hamkq%lmnmax,&
&     gs_hamkq%mgfft,mpi_enreg,gs_hamkq%natom,gs_hamkq%nattyp,gs_hamkq%ngfft,&
&     gs_hamkq%nloalg,gs_hamkq%npw_k,gs_hamkq%nspinor,gs_hamkq%ntypat,gs_hamkq%phkxred,&
&     gs_hamkq%ph1d,gs_hamkq%ph3d_k,gs_hamkq%ucvol,gs_hamkq%useylm)
     call pawaccrhoij(gs_hamkq%atindx,cplex_cprj,cwaveprj_tmp,cwaveprj1,ipert,isppol,&
&     my_natom,gs_hamkq%natom,nspinor,occ_k(iband),option_rhoij,pawrhoij1,usetimerev,wtk_k, &
&     comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
     call pawcprj_free(cwaveprj_tmp)
     ABI_DATATYPE_DEALLOCATE(cwaveprj_tmp)
   end if

 end if

 DBG_EXIT("COLL")

end subroutine dfpt_accrho
!!***
