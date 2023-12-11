!!****m* ABINIT/m_opernlc_ylm_ompgpu
!! NAME
!!  m_opernlc_ylm_ompgpu
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2022 ABINIT group (MT)
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

module m_opernlc_ylm_ompgpu

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi

 use defs_abitypes, only : MPI_type

 implicit none

 private
!!***

 public :: opernlc_ylm_ompgpu
!!***

contains
!!***

!!****f* ABINIT/opernlc_ylm_ompgpu
!! NAME
!! opernlc_ylm_ompgpu
!!
!! FUNCTION
!! * Operate with the non-local part of the hamiltonian,
!!   in order to reduce projected scalars
!! * Operate with the non-local projectors and the overlap matrix,
!!   in order to reduce projected scalars
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms (gives the absolute index of
!!                 an atom from its rank in a block of atoms)
!!  cplex=1 if <p_lmn|c> scalars are real (equivalent to istwfk>1)
!!        2 if <p_lmn|c> scalars are complex
!!  cplex_dgxdt(ndgxdt) = used only when cplex = 1
!!             cplex_dgxdt(i) = 1 if dgxdt(1,i,:,:)   is real, 2 if it is pure imaginary
!!  cplex_enl=1 if enl factors are real, 2 if they are complex
!!  cplex_fac=1 if gxfac scalars are real, 2 if gxfac scalars are complex
!!  dgxdt(cplex,ndgxdt,nlmn,nincat)=grads of projected scalars (only if optder>0)
!!  dimenl1,dimenl2=dimensions of enl (see enl)
!!  dimekbq=1 if enl factors do not contain a exp(-iqR) phase, 2 is they do
!!  enl(cplex_enl*dimenl1,dimenl2,nspinortot**2,dimekbq)=
!!  ->Norm conserving : ==== when paw_opt=0 ====
!!                      (Real) Kleinman-Bylander energies (hartree)
!!                      dimenl1=lmnmax  -  dimenl2=ntypat
!!                      dimekbq is 2 if Enl contains a exp(-iqR) phase, 1 otherwise
!!  ->PAW :             ==== when paw_opt=1, 2 or 4 ====
!!                      (Real or complex, hermitian) Dij coefs to connect projectors
!!                      dimenl1=cplex_enl*lmnmax*(lmnmax+1)/2  -  dimenl2=natom
!!                      These are complex numbers if cplex_enl=2
!!                        enl(:,:,1) contains Dij^up-up
!!                        enl(:,:,2) contains Dij^dn-dn
!!                        enl(:,:,3) contains Dij^up-dn (only if nspinor=2)
!!                        enl(:,:,4) contains Dij^dn-up (only if nspinor=2)
!!                      dimekbq is 2 if Dij contains a exp(-iqR) phase, 1 otherwise
!!  gx(cplex,nlmn,nincat*abs(enl_opt))= projected scalars
!!  iatm=absolute rank of first atom of the current block of atoms
!!  indlmn(6,nlmn)= array giving l,m,n,lm,ln,s for i=lmn
!!  itypat=type of atoms
!!  lambda=factor to be used when computing (Vln-lambda.S) - only for paw_opt=2
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell
!!  ndgxdt=second dimension of dgxdt
!!  ndgxdtfac=second dimension of dgxdtfac
!!  nincat=number of atoms in the subset here treated
!!  nlmn=number of (l,m,n) numbers for current type of atom
!!  nspinor= number of spinorial components of the wavefunctions (on current proc)
!!  nspinortot=total number of spinorial components of the wavefunctions
!!  optder=0=only gxfac is computed, 1=both gxfac and dgxdtfac are computed
!!         2=gxfac, dgxdtfac and d2gxdtfac are computed
!!  paw_opt= define the nonlocal operator concerned with:
!!           paw_opt=0 : Norm-conserving Vnl (use of Kleinman-Bylander ener.)
!!           paw_opt=1 : PAW nonlocal part of H (use of Dij coeffs)
!!           paw_opt=2 : PAW: (Vnl-lambda.Sij) (Sij=overlap matrix)
!!           paw_opt=3 : PAW overlap matrix (Sij)
!!           paw_opt=4 : both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
!!  sij(nlm*(nlmn+1)/2)=overlap matrix components (only if paw_opt=2, 3 or 4)
!!
!! OUTPUT
!!  if (paw_opt=0, 1, 2 or 4)
!!    gxfac(cplex_fac,nlmn,nincat,nspinor)= reduced projected scalars related to Vnl (NL operator)
!!  if (paw_opt=3 or 4)
!!    gxfac_sij(cplex,nlmn,nincat,nspinor)= reduced projected scalars related to Sij (overlap)
!!  if (optder==1.and.paw_opt=0, 1, 2 or 4)
!!    dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)= gradients of gxfac related to Vnl (NL operator)
!!  if (optder==1.and.paw_opt=3 or 4)
!!    dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat,nspinor)= gradients of gxfac related to Sij (overlap)
!!
!! NOTES
!! This routine operates for one type of atom, and within this given type of atom,
!! for a subset of at most nincat atoms.
!!
!! About the non-local factors symmetry:
!!   - The lower triangular part of the Dij matrix can be deduced from the upper one
!!     with the following relation: D^s2s1_ji = (D^s1s2_ij)^*
!!     where s1,s2 are spinor components
!!   - The Dij factors can contain a exp(-iqR) phase
!!     This phase does not have to be included in the symmetry rule
!!     For that reason, we first apply the real part (cos(qR).D^s1s2_ij)
!!     then, we apply the imaginary part (-sin(qR).D^s1s2_ij)
!!
!! PARENTS
!!      m_gpu_nonlop,m_nonlop_ylm
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

subroutine opernlc_ylm_ompgpu(atindx1,cplex,cplex_dgxdt,cplex_d2gxdt,cplex_enl,cplex_fac,&
&          dgxdt,dgxdtfac,dgxdtfac_sij,d2gxdt,d2gxdtfac,d2gxdtfac_sij,dimenl1,dimenl2,dimekbq,enl,&
&          gx,gxfac,gxfac_sij,iatm,indlmn,itypat,lambda,mpi_enreg,natom,ndgxdt,ndgxdtfac,&
&          nd2gxdt,nd2gxdtfac,nincat,nlmn,nspinor,nspinortot,optder,paw_opt,sij,ndat,ibeg,iend,nprojs,ntypat)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,cplex_enl,cplex_fac,dimenl1,dimenl2,dimekbq,iatm,itypat
 integer,intent(in) :: natom,ndgxdt,ndgxdtfac,nd2gxdt,nd2gxdtfac,nincat,nspinor,nspinortot,optder,paw_opt
 integer,intent(inout) :: nlmn
 integer,intent(in) :: ndat,ibeg,iend,nprojs,ntypat
 real(dp) :: lambda(ndat)
 type(MPI_type) , intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom),indlmn(6,nlmn,ntypat),cplex_dgxdt(ndgxdt),cplex_d2gxdt(nd2gxdt)
 real(dp),intent(in) :: dgxdt(cplex,ndgxdt,nlmn,nincat,nspinor)
 real(dp),intent(in) :: d2gxdt(cplex,nd2gxdt,nlmn,nincat,nspinor)
 real(dp),intent(in),target :: enl(dimenl1,dimenl2,nspinortot**2,dimekbq)
 real(dp),intent(inout) :: gx(cplex,nprojs,nspinor*ndat)
 real(dp),intent(in) :: sij(:,:)
 real(dp),intent(out),target :: dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)
 real(dp),intent(out) :: dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat,nspinor*(paw_opt/3))
 real(dp),intent(out),target :: d2gxdtfac(cplex_fac,nd2gxdtfac,nlmn,nincat,nspinor)
 real(dp),intent(out) :: d2gxdtfac_sij(cplex,nd2gxdtfac,nlmn,nincat,nspinor*(paw_opt/3))
 real(dp),intent(out),target :: gxfac(cplex_fac,nprojs,nspinor*ndat)
 real(dp),intent(out) :: gxfac_sij(cplex,nprojs,nspinor*ndat)

!Tested usecases :
! - Nvidia GPUs : FC_NVHPC + CUDA
! An eventual Intel implementation would use the OneAPI LLVM compiler.
! Homemade CUDA/HIP interfaces would allow the use of GCC.
! But it is likely that OpenMP performance won't be optimal outside GPU vendors compilers.
#ifdef HAVE_OPENMP_OFFLOAD

!Local variables-------------------------------
!Arrays
!scalars
 integer :: cplex_,ia,ierr,ijlmn,ijspin,ilm,ilmn,i0lmn,iln,index_enl,iphase,ispinor,ispinor_index,idat
 integer :: j0lmn,jilmn,jispin,jjlmn,jlm,jlmn,jspinor,jspinor_index,mu,shift
 real(dp) :: sijr
!arrays
 real(dp) :: enl_(2),gxfi(2),gxi(cplex),gxj(cplex)
 real(dp),allocatable :: d2gxdtfac_offdiag(:,:,:,:,:),dgxdtfac_offdiag(:,:,:,:,:)
 real(dp),allocatable :: gxfac_offdiag(:,:,:,:),gxfj(:,:)
 real(dp),pointer :: d2gxdtfac_(:,:,:,:,:),dgxdtfac_(:,:,:,:,:),gxfac_(:,:,:)
 real(dp),pointer :: enl_ptr(:,:,:)

! *************************************************************************

 if (dimekbq==2) then
   ABI_ERROR('dimekbq=2 not yet allowed with GPU !')
 end if

 DBG_ENTER("COLL")

!Parallelization over spinors treatment
 shift=0;if (mpi_enreg%paral_spinor==1) shift=mpi_enreg%me_spinor
 iphase=1
 gxfac_ => gxfac
 dgxdtfac_ => dgxdtfac



!Accumulate gxfac related to non-local operator (Norm-conserving)
!-------------------------------------------------------------------
  if (paw_opt==0) then
   !Enl is E(Kleinman-Bylander)
   ABI_CHECK(cplex_enl/=2,"BUG: invalid cplex_enl=2!")
   ABI_CHECK(cplex_fac==cplex,"BUG: invalid cplex_fac/=cplex!")
   !$OMP TARGET TEAMS DISTRIBUTE &
   !$OMP& PARALLEL DO COLLAPSE(3) &
   !$OMP& MAP(to:gxfac_,gx,enl_ptr) &
   !$OMP& IS_DEVICE_PTR(indlmn) &
   !$OMP& PRIVATE(ispinor,ispinor_index,ia,ilmn,iln,enl_)
   do ispinor=1,nspinor
     do ia=1,nincat
       do ilmn=1,nlmn
         ispinor_index=ispinor+shift
         iln=indlmn(5,ilmn,itypat)
         enl_(1)=enl_ptr(iln,itypat,ispinor_index)
         !TODO gxfac_(1:cplex,ilmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
         !& enl_(1)*gx(1:cplex,ilmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
       end do
     end do
   end do
   !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
  end if

!Accumulate gxfac related to nonlocal operator (PAW)
!-------------------------------------------------------------------
  if (paw_opt==1.or.paw_opt==2.or.paw_opt==4) then
   !Enl is psp strength Dij or (Dij-lambda.Sij)

!  === Diagonal term(s) (up-up, down-down)

    !$OMP TASKWAIT
!  1-Enl is real
   if (cplex_enl==1) then
     if (paw_opt==2) then
       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
       !$OMP& MAP(to:gxfac,enl,atindx1,gx,sij,lambda) &
       !$OMP& PRIVATE(idat,ispinor,ispinor_index,ia,index_enl,jlmn,j0lmn,jjlmn,ilmn,ijlmn)
       do idat=1,ndat
       do ispinor=1,nspinor
         do ia=1,nincat
           do jlmn=1,nlmn
             ispinor_index=ispinor+shift
             index_enl=atindx1(iatm+ia)
             j0lmn=jlmn*(jlmn-1)/2
             jjlmn=j0lmn+jlmn
             gxfac(1:cplex,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
&               gxfac(1:cplex,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) + &
&               (enl(jjlmn,index_enl,ispinor_index,iphase)-lambda(idat) * sij(jjlmn,itypat)) * &
&               gx(1:cplex,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
             do ilmn=1,jlmn-1
               ijlmn=j0lmn+ilmn
               gxfac(1:cplex,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
&                 gxfac(1:cplex,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) + &
&                 (enl(ijlmn,index_enl,ispinor_index,iphase)-lambda(idat) * sij(ijlmn,itypat)) * &
&                 gx(1:cplex,ilmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
             end do
             if(jlmn<nlmn) then
               do ilmn=jlmn+1,nlmn
                 i0lmn=(ilmn*(ilmn-1)/2)
                 ijlmn=i0lmn+jlmn
                 gxfac(1:cplex,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
&                   gxfac(1:cplex,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) + &
&                   (enl(ijlmn,index_enl,ispinor_index,iphase)-lambda(idat) * sij(ijlmn,itypat)) *&
&                   gx(1:cplex,ilmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
               end do
             end if
           end do
         end do
       end do
       end do

       !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
     else
       !$OMP TARGET TEAMS DISTRIBUTE &
       !$OMP& PARALLEL DO COLLAPSE(4) &
       !$OMP& MAP(to:enl,atindx1,gx,gxfac) &
       !$OMP& PRIVATE(idat,ispinor,ispinor_index,ia,index_enl,jlmn,j0lmn,jjlmn,ilmn,ijlmn) nowait
       do idat=1,ndat
       do ispinor=1,nspinor
         do ia=1,nincat
           do jlmn=1,nlmn
             ispinor_index=ispinor+shift
             index_enl=atindx1(iatm+ia)
             j0lmn=jlmn*(jlmn-1)/2
             jjlmn=j0lmn+jlmn
             gxfac(1:cplex,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
&               gxfac(1:cplex,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) + &
&               enl(jjlmn,index_enl,ispinor_index,iphase) * &
&               gx(1:cplex,ibeg+jlmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor)
             do ilmn=1,jlmn-1
               ijlmn=j0lmn+ilmn
               gxfac(1:cplex,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
&                 gxfac(1:cplex,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) + &
&                 enl(ijlmn,index_enl,ispinor_index,iphase)*gx(1:cplex,ibeg+ilmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor)
             end do
             if(jlmn<nlmn) then
               do ilmn=jlmn+1,nlmn
                 i0lmn=(ilmn*(ilmn-1)/2)
                 ijlmn=i0lmn+jlmn
                 gxfac(1:cplex,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
&                   gxfac(1:cplex,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) + &
&                   enl(ijlmn,index_enl,ispinor_index,iphase)*gx(1:cplex,ibeg+ilmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor)
               end do
             end if
           end do
         end do
       end do
       end do
       !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
     endif
    !$OMP TASKWAIT


!    2-Enl is complex  ===== D^ss'_ij=D^s's_ji^*
   else
     ABI_CHECK(cplex_fac==cplex_enl,"BUG: invalid cplex_fac/=cplex_enl!")

     if (nspinortot==1) then ! -------------> NO SPINORS
       if(paw_opt==2) then
         !$OMP TARGET TEAMS DISTRIBUTE &
         !$OMP& PARALLEL DO COLLAPSE(3) &
         !$OMP& MAP(to:gxfac_,gx,gxi,atindx1,gxj,enl_ptr) &
         !$OMP& IS_DEVICE_PTR(sij) &
         !$OMP& PRIVATE(idat,ia,index_enl,jlmn,j0lmn,jjlmn,enl_,gxj,ilmn,ijlmn,gxi)
         do idat=1,ndat
         do ia=1,nincat
           do jlmn=1,nlmn
             index_enl=atindx1(iatm+ia)
             j0lmn=jlmn*(jlmn-1)/2
             jjlmn=j0lmn+jlmn
             enl_(1)=enl_ptr(2*jjlmn-1,index_enl,1)-lambda(idat)*sij(jjlmn,itypat)
             gxj(1:cplex)=gx(1:cplex,jlmn+(ia-1)*nlmn+ibeg,1)
             gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&               gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(1)*gxj(1)
             if (cplex==2) then
               gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                 gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(1)*gxj(2)
             end if
             do ilmn=1,jlmn-1
               ijlmn=j0lmn+ilmn
               enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,1)
               enl_(1)=enl_(1)-lambda(idat)*sij(ijlmn,itypat)
               gxi(1:cplex)=gx(1:cplex,ilmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))
               gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                 gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(1)*gxi(1)
               gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                 gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))-enl_(2)*gxi(1)
               gxfac_(1,ilmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                 gxfac_(1,ilmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(1)*gxj(1)
               gxfac_(2,ilmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                 gxfac_(2,ilmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(2)*gxj(1)
               if (cplex==2) then
                 gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                   gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(2)*gxi(2)
                 gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                   gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(1)*gxi(2)
                 gxfac_(1,ilmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                   gxfac_(1,ilmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))-enl_(2)*gxj(2)
                 gxfac_(2,ilmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                   gxfac_(2,ilmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(1)*gxj(2)
               end if
             end do
             if(jlmn<nlmn) then
               do ilmn=jlmn+1,nlmn
                 i0lmn=ilmn*(ilmn-1)/2
                 ijlmn=i0lmn+jlmn
                 enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,1)
                 enl_(1)=enl_(1)-lambda(idat)*sij(ijlmn,itypat)
                 gxi(1:cplex)=gx(1:cplex,ilmn+(ia-1)*nlmn+ibeg,1)
                 gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                   gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(1)*gxi(1)
                 gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))= &
&                   gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(2)*gxi(1)
                 if (cplex==2) then
                   gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                     gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))-enl_(2)*gxi(2)
                   gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                     gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(1)*gxi(2)
                 end if
               end do
             end if
           end do
         end do
         end do
         !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
       else
         !$OMP TARGET TEAMS DISTRIBUTE &
         !$OMP& PARALLEL DO COLLAPSE(3) &
         !$OMP& MAP(to:gxfac_,gx,gxi,atindx1,gxj,enl_ptr) &
         !$OMP& IS_DEVICE_PTR(sij) &
         !$OMP& PRIVATE(idat,ia,index_enl,jlmn,j0lmn,jjlmn,enl_,gxj,ilmn,ijlmn,gxi)
         do idat=1,ndat
         do ia=1,nincat
           do jlmn=1,nlmn
           index_enl=atindx1(iatm+ia)
             j0lmn=jlmn*(jlmn-1)/2
             jjlmn=j0lmn+jlmn
             enl_(1)=enl_ptr(2*jjlmn-1,index_enl,1)
             gxj(1:cplex)=gx(1:cplex,jlmn+(ia-1)*nlmn+ibeg,1)
             gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&               gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(1)*gxj(1)
             if (cplex==2) then
               gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) + &
&                 enl_(1)*gxj(2)
             end if
             do ilmn=1,jlmn-1
               ijlmn=j0lmn+ilmn
               enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,1)
               gxi(1:cplex)=gx(1:cplex,ilmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))
               gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                 gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(1)*gxi(1)
               gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                 gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))-enl_(2)*gxi(1)
               gxfac_(1,ilmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                 gxfac_(1,ilmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(1)*gxj(1)
               gxfac_(2,ilmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                 gxfac_(2,ilmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(2)*gxj(1)
               if (cplex==2) then
                 gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                   gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(2)*gxi(2)
                 gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                   gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(1)*gxi(2)
                 gxfac_(1,ilmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                   gxfac_(1,ilmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))-enl_(2)*gxj(2)
                 gxfac_(2,ilmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                   gxfac_(2,ilmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(1)*gxj(2)
               end if
             end do
             if(jlmn<nlmn) then
               do ilmn=jlmn+1,nlmn
                 i0lmn=ilmn*(ilmn-1)/2
                 ijlmn=i0lmn+jlmn
                 enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,1)
   !TODO               gxi(1:cplex)=gx(1:cplex,ilmn+(ia-1)*nlmn+ibeg,1)
                 gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                   gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(1)*gxi(1)
                 gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                   gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(2)*gxi(1)
                 if (cplex==2) then
                   gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) - &
&                     enl_(2)*gxi(2)
                   gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) + &
&                     enl_(1)*gxi(2)
                 end if
               end do
             end if
           end do
         end do
         end do
         !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
       end if

     else ! -------------> SPINORIAL CASE

!  === Diagonal term(s) (up-up, down-down)

       !$OMP TARGET TEAMS DISTRIBUTE &
       !$OMP& PARALLEL DO COLLAPSE(4) &
       !$OMP& MAP(to:gxfac_,gx,gxi,atindx1,gxj,enl_ptr) &
       !$OMP& IS_DEVICE_PTR(sij) &
       !$OMP& PRIVATE(idat,ispinor,ispinor_index,ia,index_enl), &
       !$OMP& PRIVATE(jlmn,j0lmn,jjlmn,enl_,gxj,ilmn,ijlmn,gxi)
       do idat=1,ndat
       do ispinor=1,nspinor
         do ia=1,nincat
           do jlmn=1,nlmn
             ispinor_index=ispinor+shift
             index_enl=atindx1(iatm+ia)
             j0lmn=jlmn*(jlmn-1)/2
             jjlmn=j0lmn+jlmn
             enl_(1)=enl_ptr(2*jjlmn-1,index_enl,ispinor_index)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda(idat)*sij(jjlmn,itypat)
!TODO              gxj(1:cplex)=gx(1:cplex,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
             gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) = &
&                 gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)+enl_(1)*gxj(1)
             if (cplex==2) then
               gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) = &
&                 gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)+enl_(1)*gxj(2)
             end if
             do ilmn=1,jlmn-1
               ijlmn=j0lmn+ilmn
               enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,ispinor_index)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda(idat)*sij(ijlmn,itypat)
 !TODO               gxi(1:cplex)=gx(1:cplex,ilmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
               gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) = &
&                 gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)+enl_(1)*gxi(1)
               gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) = &
&                 gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)-enl_(2)*gxi(1)
               if (cplex==2) then
                 gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) = &
&                   gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)+enl_(2)*gxi(2)
                 gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) = &
&                   gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)+enl_(1)*gxi(2)
               end if
             end do
             if(jlmn<nlmn) then
               do ilmn=jlmn+1,nlmn
                 i0lmn=ilmn*(ilmn-1)/2
                 ijlmn=i0lmn+jlmn
                 enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,ispinor_index)
                 if (paw_opt==2) enl_(1)=enl_(1)-lambda(idat)*sij(ijlmn,itypat)
!TODO                  gxi(1:cplex)=gx(1:cplex,ilmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
                 gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) = &
&                   gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)+enl_(1)*gxi(1)
                 gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) = &
&                   gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)+enl_(2)*gxi(1)
                 if (cplex==2) then
                   gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) = &
&                     gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)-enl_(2)*gxi(2)
                   gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) = &
&                     gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)+enl_(1)*gxi(2)
                 end if
               end do
             end if
           end do
         end do
       end do
       end do
       !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
     end if !nspinortot
   end if !complex_enl

!  === Off-diagonal term(s) (up-down, down-up)

!  --- No parallelization over spinors ---
   if (nspinortot==2.and.nspinor==nspinortot) then
     ABI_CHECK(cplex_enl==2,"BUG: invalid cplex_enl/=2!")
     ABI_CHECK(cplex_fac==cplex,"BUG: invalid cplex_fac/=cplex)!")
     !$OMP TARGET TEAMS DISTRIBUTE &
     !$OMP& PARALLEL DO COLLAPSE(4) &
     !$OMP& MAP(to:gxfac_,gx,gxi,atindx1,gxj,enl_ptr) &
     !$OMP& PRIVATE(idat,ispinor,jspinor,ia,index_enl), &
     !$OMP& PRIVATE(jlmn,j0lmn,jjlmn,enl_,gxi,gxj,ilmn,ijlmn)
     do idat=1,ndat
     do ispinor=1,nspinortot
       do ia=1,nincat
         do jlmn=1,nlmn
       jspinor=3-ispinor
         index_enl=atindx1(iatm+ia)
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1:2)=enl_ptr(2*jjlmn-1:2*jjlmn,index_enl,2+ispinor )
!TODO            gxi(1:cplex)=gx(1:cplex,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
           gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,jspinor)=gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,jspinor)+enl_(1)*gxi(1)
           gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,jspinor)=gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,jspinor)-enl_(2)*gxi(1)
           if (cplex==2) then
             gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,jspinor)=gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,jspinor)+enl_(2)*gxi(2)
             gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,jspinor)=gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,jspinor)+enl_(1)*gxi(2)
           end if
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor+(idat-1)*nspinor)
!TODO              gxi(1:cplex)=gx(1:cplex,ilmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
             gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,jspinor)=gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,jspinor)+enl_(1)*gxi(1)
             gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,jspinor)=gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,jspinor)-enl_(2)*gxi(1)
             if (cplex==2) then
               gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,jspinor)=gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,jspinor)+enl_(2)*gxi(2)
               gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,jspinor)=gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,jspinor)+enl_(1)*gxi(2)
             end if
           end do
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=ilmn*(ilmn-1)/2
               ijlmn=i0lmn+jlmn
               enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor+(idat-1)*nspinor)
!TODO                gxi(1:cplex)=gx(1:cplex,ilmn+(ia-1)*nlmn+ibeg,jspinor)
               gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) = &
&                   gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)+enl_(1)*gxi(1)
               gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) = &
&                   gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)+enl_(2)*gxi(1)
               if (cplex==2) then
                 gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) = &
&                     gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)-enl_(2)*gxi(2)
                 gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) = &
&                     gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)+enl_(1)*gxi(2)
               end if
             end do
           end if
         end do
       end do
     end do
     end do
     !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

!    --- Parallelization over spinors ---
   else if (nspinortot==2.and.nspinor/=nspinortot) then
     ABI_CHECK(cplex_enl==2,"BUG: invalid cplex_enl/=2!")
     ABI_CHECK(cplex_fac==2,"BUG: invalid cplex_fac/=2!")
     ABI_MALLOC(gxfac_offdiag,(cplex_fac,nlmn,nincat,nspinortot))
!$OMP WORKSHARE
     gxfac_offdiag(:,:,:,:)=zero
!$OMP END WORKSHARE
     ispinor_index=mpi_enreg%me_spinor+1
     jspinor_index=3-ispinor_index
     if (ispinor_index==1) then
       ijspin=3;jispin=4
     else
       ijspin=4;jispin=3
     end if
     !$OMP TARGET TEAMS DISTRIBUTE &
     !$OMP& PARALLEL DO COLLAPSE(4) &
     !$OMP& MAP(to:gx,gxi,atindx1,enl_ptr) &
     !$OMP PRIVATE(idat,ia,index_enl,jlmn,j0lmn,ilmn,i0lmn,ijlmn,enl_,jilmn,gxi)
     do idat=1,ndat
     do ia=1,nincat
       do jlmn=1,nlmn
         do ilmn=1,nlmn
       index_enl=atindx1(iatm+ia)
         j0lmn=jlmn*(jlmn-1)/2
           i0lmn=ilmn*(ilmn-1)/2
           if (ilmn<=jlmn) then
             ijlmn=j0lmn+ilmn
             enl_(1)= enl_ptr(2*ijlmn-1,index_enl,ijspin)
             enl_(2)=-enl_ptr(2*ijlmn  ,index_enl,ijspin)
           else
             jilmn=i0lmn+jlmn
             enl_(1)= enl_ptr(2*jilmn-1,index_enl,jispin)
             enl_(2)= enl_ptr(2*jilmn  ,index_enl,jispin)
           end if
!TODO            gxi(1:cplex)=gx(1:cplex,ilmn,ia,1+nspinor*(idat-1))
           gxfac_offdiag(1,jlmn,ia,jspinor_index)= &
&            gxfac_offdiag(1,jlmn,ia,jspinor_index)+enl_(1)*gxi(1)
           gxfac_offdiag(2,jlmn,ia,jspinor_index)= &
&            gxfac_offdiag(2,jlmn,ia,jspinor_index)+enl_(2)*gxi(1)
           if (cplex==2) then
             gxfac_offdiag(1,jlmn,ia,jspinor_index)= &
&              gxfac_offdiag(1,jlmn,ia,jspinor_index)-enl_(2)*gxi(2)
             gxfac_offdiag(2,jlmn,ia,jspinor_index)= &
&              gxfac_offdiag(2,jlmn,ia,jspinor_index)+enl_(1)*gxi(2)
           end if
         end do !ilmn
       end do !jlmn
     end do !iat
     end do !idat
     !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
     call xmpi_sum(gxfac_offdiag,mpi_enreg%comm_spinor,ierr)
     !$OMP TARGET UPDATE FROM(gxfac_)
     !gxfac_(:,:,:,1)=gxfac_(:,:,:,1)+gxfac_offdiag(:,:,:,ispinor_index)
     ABI_FREE(gxfac_offdiag)
   end if

  end if !paw_opt


!Accumulate gxfac related to overlap (Sij) (PAW)
!------------------------------------------- ------------------------
 !$OMP TASKWAIT
 if (paw_opt==3.or.paw_opt==4) then ! Use Sij, overlap contribution
   !$OMP TARGET TEAMS DISTRIBUTE &
   !$OMP& PARALLEL DO COLLAPSE(4) &
   !$OMP& MAP(to:sij,gx,gxfac_sij) &
   !$OMP PRIVATE(idat, ispinor,ia,jlmn,j0lmn,jjlmn,jlm,ilmn,ilm,ijlmn) nowait
   do idat=1,ndat
   do ispinor=1,nspinor
     do ia=1,nincat
       do jlmn=1,nlmn
         j0lmn=jlmn*(jlmn-1)/2
         jjlmn=j0lmn+jlmn
         gxfac_sij(1:cplex,ibeg+jlmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor)= &
           gxfac_sij(1:cplex,ibeg+jlmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor) &
           + sij(jjlmn,itypat) * gx(1:cplex,ibeg+jlmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor)
         do ilmn=1,jlmn-1
           ijlmn=j0lmn+ilmn
           gxfac_sij(1:cplex,ibeg+jlmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor)= &
             gxfac_sij(1:cplex,ibeg+jlmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor) &
             + sij(ijlmn,itypat) * gx(1:cplex,ibeg+ilmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor)
         end do
         if(jlmn<nlmn) then
           do ilmn=jlmn+1,nlmn
             i0lmn=ilmn*(ilmn-1)/2
             ijlmn=i0lmn+jlmn
             gxfac_sij(1:cplex,ibeg+jlmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor)=&
               gxfac_sij(1:cplex,ibeg+jlmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor) &
               + sij(ijlmn,itypat) * gx(1:cplex,ibeg+ilmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor)
           end do
         end if
       end do
     end do
   end do
   end do
   !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
 end if
 !$OMP TASKWAIT

#else

 ABI_UNUSED((/cplex,cplex_enl,cplex_fac,dimenl1,dimenl2,dimekbq,iatm,itypat/))
 ABI_UNUSED((/natom,ndgxdt,ndgxdtfac,nd2gxdt,nd2gxdtfac,nincat,nspinor,nspinortot,optder,paw_opt/))
 ABI_UNUSED((/nlmn,ndat,ibeg,iend,nprojs,ntypat/))
 ABI_UNUSED((/atindx1,indlmn,cplex_dgxdt,cplex_d2gxdt/))
 ABI_UNUSED((/dgxdtfac,dgxdtfac_sij,d2gxdtfac,d2gxdtfac_sij,gxfac,gxfac_sij,dgxdt,d2gxdt,lambda,enl,gx,sij/))
 ABI_UNUSED_A(mpi_enreg)
 ABI_BUG("Unhandled configuration for OpenMP GPU immplementation")

#endif

end subroutine opernlc_ylm_ompgpu
!!***

end module m_opernlc_ylm_ompgpu
!!***
