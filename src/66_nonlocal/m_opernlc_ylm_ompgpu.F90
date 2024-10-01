!!****m* ABINIT/m_opernlc_ylm_ompgpu
!! NAME
!!  m_opernlc_ylm_ompgpu
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2024 ABINIT group (MT)
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
 use m_abi_linalg
 use, intrinsic :: iso_c_binding

 use defs_abitypes, only : MPI_type


 implicit none

 private
!!***

 public :: opernlc_ylm_ompgpu
!!***

! Work buffers to be used when iphase==2
 real(dp), allocatable, target :: d2gxdtfac_2ndphase(:,:,:,:)
 real(dp), allocatable, target :: dgxdtfac_2ndphase(:,:,:,:)
 real(dp), allocatable, target :: gxfac_2ndphase(:,:,:)

!----------------------------------------------------------------------

contains
!!***

!----------------------------------------------------------------------

!!****f* m_opernlc_ylm_ompgpu/alloc_work_arrays
!! NAME
!! alloc_work_arrays
!!
!! FUNCTION
!! Allocation of work arrays
!!
!! INPUTS
!!
!! SOURCE
 subroutine alloc_work_arrays(optder,cplex_fac,ndgxdtfac,nd2gxdtfac,nprojs,ndat)

  integer,intent(in) :: optder,cplex_fac,ndgxdtfac,nd2gxdtfac,nprojs,ndat

! *************************************************************************

   ABI_MALLOC(gxfac_2ndphase,(cplex_fac,nprojs,ndat))
   !$OMP TARGET ENTER DATA MAP(alloc:gxfac_2ndphase)
   call gpu_set_to_zero(gxfac_2ndphase, int(cplex_fac,c_size_t)*nprojs*ndat)
   if (optder>=1) then
     ABI_MALLOC(dgxdtfac_2ndphase,(cplex_fac,ndgxdtfac,nprojs,ndat))
     !$OMP TARGET ENTER DATA MAP(alloc:dgxdtfac_2ndphase)
     call gpu_set_to_zero(dgxdtfac_2ndphase, int(cplex_fac,c_size_t)*ndgxdtfac*nprojs*ndat)
   end if
   if (optder>=2) then
     ABI_MALLOC(d2gxdtfac_2ndphase,(cplex_fac,nd2gxdtfac,nprojs,ndat))
     !$OMP TARGET ENTER DATA MAP(alloc:dgxdtfac_2ndphase)
     call gpu_set_to_zero(d2gxdtfac_2ndphase, int(cplex_fac,c_size_t)*nd2gxdtfac*nprojs*ndat)
   end if

 end subroutine alloc_work_arrays
!!***

!----------------------------------------------------------------------

!!****f* m_opernlc_ylm_ompgpu/destroy_work_arrays
!! NAME
!! destroy_work_arrays
!!
!! FUNCTION
!! Destruction of work arrays
!!
!! INPUTS
!!
!! SOURCE
 subroutine destroy_work_arrays()

! *************************************************************************

  if(allocated(gxfac_2ndphase)) then
    !$OMP TARGET EXIT DATA MAP(delete:gxfac_2ndphase)
    ABI_FREE(gxfac_2ndphase)
  end if
  if(allocated(dgxdtfac_2ndphase)) then
    !$OMP TARGET EXIT DATA MAP(delete:dgxdtfac_2ndphase)
    ABI_FREE(dgxdtfac_2ndphase)
  end if
  if(allocated(d2gxdtfac_2ndphase)) then
    !$OMP TARGET EXIT DATA MAP(delete:d2gxdtfac_2ndphase)
    ABI_FREE(d2gxdtfac_2ndphase)
  end if

 end subroutine destroy_work_arrays
!!***

!----------------------------------------------------------------------

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
&          nd2gxdt,nd2gxdtfac,nincat,nlmn,nspinor,nspinortot,optder,paw_opt,sij,ndat,ibeg,iend,nprojs,ntypat,gpu_option)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,cplex_enl,cplex_fac,dimenl1,dimenl2,dimekbq,iatm,itypat
 integer,intent(in) :: natom,ndgxdt,ndgxdtfac,nd2gxdt,nd2gxdtfac,nincat,nspinor,nspinortot,optder,paw_opt,gpu_option
 integer,intent(inout) :: nlmn
 integer,intent(in) :: ndat,ibeg,iend,nprojs,ntypat
 real(dp) :: lambda(ndat)
 type(MPI_type) , intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom),indlmn(6,nlmn,ntypat),cplex_dgxdt(ndgxdt),cplex_d2gxdt(nd2gxdt)
 real(dp),intent(in) :: dgxdt(cplex,ndgxdt,nprojs,nspinor*ndat)
 real(dp),intent(in) :: d2gxdt(cplex,nd2gxdt,nlmn,nincat,nspinor*ndat)
 real(dp),intent(in),target :: enl(dimenl1,dimenl2,nspinortot**2,dimekbq)
 real(dp),intent(inout) :: gx(cplex,nprojs,nspinor*ndat)
 real(dp),intent(in) :: sij(:)
 real(dp),intent(out),target :: dgxdtfac(cplex_fac,ndgxdtfac,nprojs,nspinor*ndat)
 real(dp),intent(out) :: dgxdtfac_sij(cplex,ndgxdtfac,nprojs,nspinor*ndat*(paw_opt/3))
 real(dp),intent(out),target :: d2gxdtfac(cplex_fac,nd2gxdtfac,nprojs,nspinor*ndat)
 real(dp),intent(out) :: d2gxdtfac_sij(cplex,nd2gxdtfac,nprojs,nspinor*ndat*(paw_opt/3))
 real(dp),intent(out),target :: gxfac(cplex_fac,nprojs,nspinor*ndat)
 real(dp),intent(out) :: gxfac_sij(cplex,nprojs,nspinor*ndat)

!Tested usecases :
! - Nvidia GPUs : FC_NVHPC + CUDA
! - AMD GPUs    : FC_LLVM + HIP
! An eventual Intel implementation would use the OneAPI LLVM compiler.
! Homemade CUDA/HIP interfaces would allow the use of GCC.
! But it is likely that OpenMP performance won't be optimal outside GPU vendors compilers.
#ifdef HAVE_OPENMP_OFFLOAD

!Local variables-------------------------------
!Arrays
!scalars
 integer :: cplex_,ia,ierr,ijlmn,ijspin,ilm,ilmn,i0lmn,iln,index_enl,iphase,ispinor,ispinor_index,idat
 integer :: j0lmn,jilmn,jispin,jjlmn,jlm,jlmn,jspinor,jspinor_index,mu,shift,ii
 real(dp) :: sijr
!arrays
 real(dp) :: enl_(2),gxfi(2),gxi(cplex),gxj(cplex),tmp(2)
 real(dp),allocatable :: gxfj(:,:)
 real(dp), ABI_CONTIGUOUS pointer :: d2gxdtfac_(:,:,:,:),dgxdtfac_(:,:,:,:),gxfac_(:,:,:)
 real(dp), ABI_CONTIGUOUS pointer :: enl_ptr(:,:,:)

! *************************************************************************

 if (nspinor==2) then
   ABI_ERROR('nspinor=2 not yet allowed with GPU !')
 end if

 if (nspinortot==2) then
   ABI_ERROR('nspinortot=2 not yet allowed with GPU !')
 end if

 DBG_ENTER("COLL")

!Parallelization over spinors treatment
 shift=0;if (mpi_enreg%paral_spinor==1) shift=mpi_enreg%me_spinor
!When Enl factors contain a exp(-iqR) phase:
! - We loop over the real and imaginary parts
! - We need an additional memory space
 do iphase=1,dimekbq
  if (paw_opt==3) cycle
  if (iphase==1) then
   gxfac_ => gxfac ; dgxdtfac_ => dgxdtfac ; d2gxdtfac_ => d2gxdtfac
  else
   ABI_CHECK(cplex_fac==2,"BUG: invalid cplex_fac==1 when dimekbq=2!")
   call alloc_work_arrays(optder,cplex_fac,ndgxdtfac,nd2gxdtfac,nprojs,ndat)
   gxfac_ => gxfac_2ndphase
   if(optder>=1) dgxdtfac_ => dgxdtfac_2ndphase
   if(optder>=2) d2gxdtfac_ => d2gxdtfac_2ndphase
  end if
  enl_ptr => enl(:,:,:,iphase)





!Accumulate gxfac related to non-local operator (Norm-conserving)
!-------------------------------------------------------------------
  if (paw_opt==0) then
   !Enl is E(Kleinman-Bylander)
   ABI_CHECK(cplex_enl/=2,"BUG: invalid cplex_enl=2!")
   ABI_CHECK(cplex_fac==cplex,"BUG: invalid cplex_fac/=cplex!")
   !$OMP TARGET TEAMS DISTRIBUTE &
   !$OMP& MAP(to:gxfac_,gx,enl_ptr,indlmn) &
   !$OMP& PRIVATE(idat,ispinor)
   do idat=1,ndat
   do ispinor=1,nspinor
     !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(ia,ilmn,iln,ii)
     do ia=1,nincat
       do ilmn=1,nlmn
         do ii=1,cplex
           iln=indlmn(5,ilmn,itypat)
           gxfac_(ii,ilmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
           & enl_ptr(iln,itypat,ispinor+shift)*gx(ii,ilmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
         end do
       end do
     end do
   end do
   end do
  end if

!Accumulate gxfac related to nonlocal operator (PAW)
!-------------------------------------------------------------------
  if (paw_opt==1.or.paw_opt==2.or.paw_opt==4) then
   !Enl is psp strength Dij or (Dij-lambda.Sij)

!  === Diagonal term(s) (up-up, down-down)

!  1-Enl is real
   if (cplex_enl==1) then
     if (paw_opt==2) then
       !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) &
       !$OMP& MAP(to:gxfac,enl_ptr,atindx1,gx,sij,lambda) &
       !$OMP& PRIVATE(idat,ispinor) &
       !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
       do idat=1,ndat
       do ispinor=1,nspinor
         !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(ia,jlmn,ispinor_index,index_enl,j0lmn,i0lmn,ii)
         do ia=1,nincat
           do jlmn=1,nlmn
             ispinor_index=ispinor+shift
             index_enl=atindx1(iatm+ia)
             j0lmn=jlmn*(jlmn-1)/2
             do ii=1,cplex
               gxfac_(ii,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
&                 gxfac_(ii,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) + &
&                 (enl_ptr(j0lmn+jlmn,index_enl,ispinor_index)-lambda(idat) * sij(j0lmn+jlmn)) * &
&                 gx(ii,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
             end do
             do ilmn=1,jlmn-1
               do ii=1,cplex
                 gxfac_(ii,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
&                   gxfac_(ii,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) + &
&                   (enl_ptr(j0lmn+ilmn,index_enl,ispinor_index)-lambda(idat) * sij(j0lmn+ilmn)) * &
&                   gx(ii,ilmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
               end do
             end do
             if(jlmn<nlmn) then
               do ilmn=jlmn+1,nlmn
                 i0lmn=(ilmn*(ilmn-1)/2)
                 do ii=1,cplex
                   gxfac_(ii,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
&                     gxfac_(ii,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) + &
&                     (enl_ptr(i0lmn+jlmn,index_enl,ispinor_index)-lambda(idat) * sij(i0lmn+jlmn)) *&
&                     gx(ii,ilmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
                 end do
               end do
             end if
           end do
         end do
       end do
       end do

     else
       !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) &
       !$OMP& MAP(to:enl_ptr,atindx1,gx,gxfac_) &
       !$OMP& PRIVATE(idat,ispinor,ispinor_index,ia,index_enl,jlmn) &
       !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
       do idat=1,ndat
       do ispinor=1,nspinor
         !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(jlmn,j0lmn,ii,ia,ispinor_index,index_enl)
         do ia=1,nincat
           do jlmn=1,nlmn
             ispinor_index=ispinor+shift
             index_enl=atindx1(iatm+ia)
             j0lmn=jlmn*(jlmn-1)/2
             do ii=1,cplex
               gxfac_(ii,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
&                    gxfac_(ii,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) + &
&                    enl_ptr(jlmn*(jlmn-1)/2+jlmn,index_enl,ispinor+shift) * &
&                    gx(ii,ibeg+jlmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor)
             end do
           end do
         end do
       end do
       end do
       !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(3) &
       !$OMP& MAP(to:enl_ptr,atindx1,gx,gxfac_) &
       !$OMP& PRIVATE(idat,ispinor,ia,ispinor_index,index_enl)
       do idat=1,ndat
       do ispinor=1,nspinor
         do ia=1,nincat
           ispinor_index=ispinor+shift
           index_enl=atindx1(iatm+ia)
           !$OMP PARALLEL DO PRIVATE(tmp,j0lmn,jlmn,ilmn,i0lmn,ii)
           do jlmn=1,nlmn
             j0lmn=jlmn*(jlmn-1)/2
             tmp(1)=0
             tmp(2)=0
             do ilmn=1,jlmn-1
               do ii=1,cplex
                 tmp(ii)=tmp(ii) + enl_ptr(j0lmn+ilmn,index_enl,ispinor_index) &
&                      * gx(ii,ibeg+ilmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor)
               end do
             end do
             if(jlmn<nlmn) then
               do ilmn=jlmn+1,nlmn
                 i0lmn=(ilmn*(ilmn-1)/2)
                 do ii=1,cplex
                   tmp(ii)=tmp(ii) + enl_ptr(i0lmn+jlmn,index_enl,ispinor_index) &
&                       *gx(ii,ibeg+ilmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor)
                 end do
               end do
             end if
             do ii=1,cplex
               gxfac_(ii,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)= &
                   gxfac_(ii,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor) + tmp(ii)
             end do
           end do
         end do
       end do
       end do
     endif


!    2-Enl is complex  ===== D^ss'_ij=D^s's_ji^*
   else
     ABI_CHECK(cplex_fac==cplex_enl,"BUG: invalid cplex_fac/=cplex_enl!")

     if (nspinortot==1) then ! -------------> NO SPINORS
       if(paw_opt==2) then
         !$OMP TARGET TEAMS DISTRIBUTE &
         !$OMP& PARALLEL DO COLLAPSE(3) &
         !$OMP& MAP(to:gxfac_,gx,gxi,atindx1,gxj,enl_ptr) &
         !$OMP& PRIVATE(idat,ia,index_enl,jlmn,j0lmn,jjlmn,enl_,gxj,ilmn,i0lmn,ijlmn,gxi) &
         !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
         do idat=1,ndat
         do ia=1,nincat
           do jlmn=1,nlmn
             index_enl=atindx1(iatm+ia)
             j0lmn=jlmn*(jlmn-1)/2
             enl_(1)=enl_ptr(2*j0lmn+jlmn-1,index_enl,1)-lambda(idat)*sij(j0lmn+jlmn)
             gxj(1:cplex)=gx(1:cplex,jlmn+(ia-1)*nlmn+ibeg,1)
             gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&               gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(1)*gxj(1)
             if (cplex==2) then
               gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&                 gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(1)*gxj(2)
             end if
             do ilmn=1,jlmn-1
               enl_(1:2)=enl_ptr(2*j0lmn+ilmn-1:2*j0lmn+ilmn,index_enl,1)
               enl_(1)=enl_(1)-lambda(idat)*sij(j0lmn+ilmn)
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
                 enl_(1:2)=enl_ptr(2*i0lmn+jlmn-1:2*i0lmn+jlmn,index_enl,1)
                 enl_(1)=enl_(1)-lambda(idat)*sij(i0lmn+jlmn)
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
       else
         !$OMP TARGET TEAMS DISTRIBUTE &
         !$OMP& PARALLEL DO COLLAPSE(3) &
         !$OMP& MAP(to:gxfac_,gx,gxi,atindx1,gxj,enl_ptr) &
         !$OMP& PRIVATE(idat,ia,index_enl,jlmn,j0lmn,jjlmn,enl_,gxj,ilmn,i0lmn,ijlmn,gxi) &
         !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
         do idat=1,ndat
         do ia=1,nincat
           do jlmn=1,nlmn
           index_enl=atindx1(iatm+ia)
             j0lmn=jlmn*(jlmn-1)/2
             enl_(1)=enl_ptr(2*j0lmn+jlmn-1,index_enl,1)
             gxj(1:cplex)=gx(1:cplex,jlmn+(ia-1)*nlmn+ibeg,1)
             gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = &
&               gxfac_(1,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1))+enl_(1)*gxj(1)
             if (cplex==2) then
               gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) = gxfac_(2,jlmn+(ia-1)*nlmn+ibeg,1+nspinor*(idat-1)) + &
&                 enl_(1)*gxj(2)
             end if
             do ilmn=1,jlmn-1
               enl_(1:2)=enl_ptr(2*j0lmn+ilmn-1:2*j0lmn+ilmn,index_enl,1)
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
                 enl_(1:2)=enl_ptr(2*i0lmn+jlmn-1:2*i0lmn+jlmn,index_enl,1)
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
       end if

     else ! -------------> SPINORIAL CASE

       ABI_BUG("nspinor==2 not supported with OpenMP GPU")
     end if !nspinortot
   end if !complex_enl

!  === Off-diagonal term(s) (up-down, down-up)

!  --- No parallelization over spinors ---
   if (nspinortot==2.and.nspinor==nspinortot) then

     ABI_BUG("nspinor==2 not supported with OpenMP GPU")
!    --- Parallelization over spinors ---
   else if (nspinortot==2.and.nspinor/=nspinortot) then
     ABI_BUG("nspinor==2 not supported with OpenMP GPU")
   end if

  end if !paw_opt


!Accumulate dgxdtfac related to nonlocal operator (Norm-conserving)
!-------------------------------------------------------------------
  if (optder>=1.and.paw_opt==0) then
   !Enl is E(Kleinman-Bylander)
   ABI_CHECK(cplex_enl==1,"BUG: invalid cplex_enl/=1!")
   ABI_CHECK(cplex_fac==cplex,"BUG: invalid cplex_fac/=cplex!")
   ABI_WARNING("DEBUG: Code section not checked!")
   !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) &
   !$OMP& MAP(to:dgxdtfac_,enl_ptr,atindx1,dgxdt,indlmn) &
   !$OMP& PRIVATE(idat,ispinor,ispinor_index,ia,index_enl,jlmn,j0lmn,jjlmn,ilmn,ijlmn) &
   !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
   do idat=1,ndat
   do ispinor=1,nspinor
     ispinor_index = ispinor + shift
     !$OMP PARALLEL DO
     do ia=1,nincat
       do ilmn=1,nlmn
         do mu=1,ndgxdtfac
           dgxdtfac_(1:cplex,mu,ilmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
           &    enl_ptr(indlmn(5,ilmn,itypat),itypat,ispinor_index)&
           &    * dgxdt(1:cplex,mu,ilmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
         end do
       end do
     end do
   end do
   end do
  end if

!Accumulate dgxdtfac related to nonlocal operator (PAW)
!-------------------------------------------------------------------
  if (optder>=1.and.(paw_opt==1.or.paw_opt==2.or.paw_opt==4)) then
   !Enl is psp strength Dij or (Dij-lambda.Sij)

!  === Diagonal term(s) (up-up, down-down)

!  1-Enl is real
   if (cplex_enl==1) then
     if (paw_opt/=2) then
       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
       !$OMP& MAP(to:dgxdtfac_,enl_ptr,atindx1,dgxdt) &
       !$OMP& PRIVATE(idat,ispinor,ispinor_index,ia,index_enl,jlmn,j0lmn,ilmn,i0lmn) &
       !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
       do idat=1,ndat
       do ispinor=1,nspinor
         do ia=1,nincat
           do jlmn=1,nlmn
             ispinor_index=ispinor+shift
             index_enl=atindx1(iatm+ia)
             j0lmn=jlmn*(jlmn-1)/2
             do mu=1,ndgxdtfac
               dgxdtfac_(1:cplex,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
               &    dgxdtfac_(1:cplex,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)&
               &    + enl_ptr(j0lmn+jlmn,index_enl,ispinor_index) &
               &    * dgxdt(1:cplex,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
             end do
             do ilmn=1,jlmn-1
               do mu=1,ndgxdtfac
                 dgxdtfac_(1:cplex,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
                 &    dgxdtfac_(1:cplex,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)&
                 &    + enl_ptr(j0lmn+ilmn,index_enl,ispinor_index) &
                 &    * dgxdt(1:cplex,mu,ilmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
               end do
             end do
             if(jlmn<nlmn) then
               do ilmn=jlmn+1,nlmn
                 i0lmn=ilmn*(ilmn-1)/2
                 do mu=1,ndgxdtfac
                   dgxdtfac_(1:cplex,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
                   &    dgxdtfac_(1:cplex,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)&
                   &    + enl_ptr(i0lmn+jlmn,index_enl,ispinor_index) &
                   &    * dgxdt(1:cplex,mu,ilmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
                 end do
               end do
             end if
           end do
         end do
       end do
       end do
     else
       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
       !$OMP& MAP(to:dgxdtfac_,enl_ptr,atindx1,dgxdt,sij,lambda) &
       !$OMP& PRIVATE(idat,ispinor,ispinor_index,ia,index_enl,jlmn,j0lmn,ilmn,i0lmn) &
       !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
       do idat=1,ndat
       do ispinor=1,nspinor
         do ia=1,nincat
           do jlmn=1,nlmn
             ispinor_index=ispinor+shift
             index_enl=atindx1(iatm+ia)
             j0lmn=jlmn*(jlmn-1)/2
             do mu=1,ndgxdtfac
               dgxdtfac_(1:cplex,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
               &    dgxdtfac_(1:cplex,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)&
               &    + (enl_ptr(j0lmn+jlmn,index_enl,ispinor_index)-lambda(idat)*sij(j0lmn+jlmn)) &
               &    * dgxdt(1:cplex,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
             end do
             do ilmn=1,jlmn-1
               do mu=1,ndgxdtfac
                 dgxdtfac_(1:cplex,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
                 &    dgxdtfac_(1:cplex,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)&
                 &    + (enl_ptr(j0lmn+ilmn,index_enl,ispinor_index)-lambda(idat)*sij(j0lmn+ilmn)) &
                 &    * dgxdt(1:cplex,mu,ilmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
               end do
             end do
             if(jlmn<nlmn) then
               do ilmn=jlmn+1,nlmn
                 i0lmn=ilmn*(ilmn-1)/2
                 do mu=1,ndgxdtfac
                   dgxdtfac_(1:cplex,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
                   &    dgxdtfac_(1:cplex,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)&
                   &    + (enl_ptr(i0lmn+jlmn,index_enl,ispinor_index)-lambda(idat)*sij(i0lmn+jlmn)) &
                   &    * dgxdt(1:cplex,mu,ilmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
                 end do
               end do
             end if
           end do
         end do
       end do
       end do
     end if

!    2-Enl is complex  ===== D^ss'_ij=D^s's_ji^*
   else
     ABI_CHECK(cplex_fac==cplex_enl,"BUG: invalid cplex_fac/=cplex_enl!")

     if (nspinortot==1) then ! -------------> NO SPINORS

       ABI_WARNING("DEBUG: Code section not checked!")
       !if (paw_opt==2) ABI_BUG("paw_opt==2 not yet handled with OpenMP GPU")

       ABI_MALLOC(gxfj,(cplex,ndgxdtfac))
       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
       !$OMP& MAP(to:dgxdtfac_,enl_,atindx1,dgxdt,sij,lambda,gxfj,enl,gxfi) &
       !$OMP& PRIVATE(idat,ispinor,ispinor_index,ia,index_enl,jlmn,j0lmn,jjlmn,ilmn,i0lmn,ijlmn) &
       !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
       do idat=1,ndat
       do ia=1,nincat
         do jlmn=1,nlmn
           index_enl=atindx1(iatm+ia)
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1)=enl(2*jjlmn-1,index_enl,1,iphase)
           if (paw_opt==2) enl_(1)=enl_(1)-lambda(ndat)*sij(jjlmn)
           do mu=1,ndgxdtfac
             if(cplex_dgxdt(mu)==2)then
               cplex_ = 2 ; gxfj(1,mu) = zero ; gxfj(2,mu) = dgxdt(1,mu,jlmn+(ia-1)*nlmn+ibeg,1)
             else
               cplex_ = cplex ; gxfj(1:cplex,mu)=dgxdt(1:cplex,mu,jlmn+(ia-1)*nlmn+ibeg,1)
             end if
             dgxdtfac_(1,mu,jlmn+(ia-1)*nlmn+ibeg,idat)=dgxdtfac_(1,mu,jlmn+(ia-1)*nlmn+ibeg,idat)+enl_(1)*gxfj(1,mu)
             if (cplex_==2) dgxdtfac_(2,mu,jlmn+(ia-1)*nlmn+ibeg,idat)=dgxdtfac_(2,mu,jlmn+(ia-1)*nlmn+ibeg,idat)+enl_(1)*gxfj(2,mu)
           end do
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,1,iphase)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda(ndat)*sij(ijlmn)
             do mu=1,ndgxdtfac
               if(cplex_dgxdt(mu)==2)then
                 cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = dgxdt(1,mu,ilmn+(ia-1)*nlmn+ibeg,idat)
               else
                 cplex_ = cplex ; gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn+(ia-1)*nlmn+ibeg,idat)
               end if
               dgxdtfac_(1,mu,jlmn+(ia-1)*nlmn+ibeg,idat)=dgxdtfac_(1,mu,jlmn+(ia-1)*nlmn+ibeg,idat)+enl_(1)*gxfi(1)
               dgxdtfac_(2,mu,jlmn+(ia-1)*nlmn+ibeg,idat)=dgxdtfac_(2,mu,jlmn+(ia-1)*nlmn+ibeg,idat)-enl_(2)*gxfi(1)
               if (cplex_==2) then
                 dgxdtfac_(1,mu,jlmn+(ia-1)*nlmn+ibeg,idat)=dgxdtfac_(1,mu,jlmn+(ia-1)*nlmn+ibeg,idat)+enl_(2)*gxfi(2)
                 dgxdtfac_(2,mu,jlmn+(ia-1)*nlmn+ibeg,idat)=dgxdtfac_(2,mu,jlmn+(ia-1)*nlmn+ibeg,idat)+enl_(1)*gxfi(2)
               end if
             end do
           end do
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=ilmn*(ilmn-1)/2
               ijlmn=i0lmn+jlmn
               enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,1,iphase)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda(ndat)*sij(ijlmn)
               do mu=1,ndgxdtfac
                 if(cplex_dgxdt(mu)==2)then
                   cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = dgxdt(1,mu,ilmn+(ia-1)*nlmn+ibeg,idat)
                 else
                   cplex_ = cplex ; gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn+(ia-1)*nlmn+ibeg,idat)
                 end if
                 dgxdtfac_(1,mu,jlmn+(ia-1)*nlmn+ibeg,idat)=dgxdtfac_(1,mu,jlmn+(ia-1)*nlmn+ibeg,idat)+enl_(1)*gxfi(1)
                 dgxdtfac_(2,mu,jlmn+(ia-1)*nlmn+ibeg,idat)=dgxdtfac_(2,mu,jlmn+(ia-1)*nlmn+ibeg,idat)+enl_(2)*gxfi(1)
                 if (cplex_==2) then
                   dgxdtfac_(1,mu,jlmn+(ia-1)*nlmn+ibeg,idat)=dgxdtfac_(1,mu,jlmn+(ia-1)*nlmn+ibeg,idat)-enl_(2)*gxfi(2)
                   dgxdtfac_(2,mu,jlmn+(ia-1)*nlmn+ibeg,idat)=dgxdtfac_(2,mu,jlmn+(ia-1)*nlmn+ibeg,idat)+enl_(1)*gxfi(2)
                 end if
               end do
             end do
           end if
         end do
       end do
       end do
       ABI_FREE(gxfj)

     else ! -------------> SPINORIAL CASE

     ABI_BUG("nspinor==2 not supported with OpenMP GPU")
     end if !nspinortot
   end if !complex

!  === Off-diagonal term(s) (up-down, down-up)

!  --- No parallelization over spinors ---
   if (nspinortot==2.and.nspinor==nspinortot) then

     ABI_BUG("nspinor==2 not supported with OpenMP GPU")
!    --- Parallelization over spinors ---
   else if (nspinortot==2.and.nspinor/=nspinortot) then
     ABI_BUG("nspinor==2 not supported with OpenMP GPU")
   end if !nspinortot

  end if ! pawopt & optder

!End of loop when a exp(-iqR) phase is present
!------------------------------------------- ------------------------

!When iphase=1, gxfac and gxfac_ point to the same memory space
!When iphase=2, we add i.gxfac_ to gxfac
  if (iphase==2) then
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
    !$OMP& PRIVATE(idat,ia,ilmn) MAP(to:gxfac,gxfac_) &
    !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
    do idat=1,ndat*nspinor
      do ia=1,nincat
        do ilmn=1,nlmn
          gxfac(1,ilmn+(ia-1)*nlmn+ibeg,idat)=&
          &    gxfac(1,ilmn+(ia-1)*nlmn+ibeg,idat)-gxfac_(2,ilmn+(ia-1)*nlmn+ibeg,idat)
          gxfac(2,ilmn+(ia-1)*nlmn+ibeg,idat)=&
          &    gxfac(2,ilmn+(ia-1)*nlmn+ibeg,idat)+gxfac_(1,ilmn+(ia-1)*nlmn+ibeg,idat)
        end do
      end do
    end do
    if (optder>=1) then
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
      !$OMP& PRIVATE(idat,ia,ilmn,mu) MAP(to:dgxdtfac,dgxdtfac_) &
      !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
      do idat=1,ndat*nspinor
        do ia=1,nincat
          do ilmn=1,nlmn
            do mu=1,ndgxdtfac
              dgxdtfac(1,mu,ilmn+(ia-1)*nlmn+ibeg,idat)=&
              &    dgxdtfac(1,mu,ilmn+(ia-1)*nlmn+ibeg,idat)-dgxdtfac_(2,mu,ilmn+(ia-1)*nlmn+ibeg,idat)
              dgxdtfac(2,mu,ilmn+(ia-1)*nlmn+ibeg,idat)=&
              &    dgxdtfac(2,mu,ilmn+(ia-1)*nlmn+ibeg,idat)+dgxdtfac_(1,mu,ilmn+(ia-1)*nlmn+ibeg,idat)
            end do
          end do
        end do
      end do
    end if
    if (optder>=2) then
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
      !$OMP& PRIVATE(idat,ia,ilmn,mu) MAP(to:d2gxdtfac,d2gxdtfac_) &
      !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
      do idat=1,ndat*nspinor
        do ia=1,nincat
          do ilmn=1,nlmn
            do mu=1,nd2gxdtfac
              d2gxdtfac(1,mu,ilmn+(ia-1)*nlmn+ibeg,ispinor)=&
              &    d2gxdtfac(1,mu,ilmn+(ia-1)*nlmn+ibeg,ispinor)-d2gxdtfac_(2,mu,ilmn+(ia-1)*nlmn+ibeg,ispinor)
              d2gxdtfac(2,mu,ilmn+(ia-1)*nlmn+ibeg,ispinor)=&
              &    d2gxdtfac(2,mu,ilmn+(ia-1)*nlmn+ibeg,ispinor)+d2gxdtfac_(1,mu,ilmn+(ia-1)*nlmn+ibeg,ispinor)
            end do
          end do
        end do
      end do
    end if
    call destroy_work_arrays()
  end if

!End loop over real/imaginary part of the exp(-iqR) phase
 end do


!Accumulate gxfac related to overlap (Sij) (PAW)
!------------------------------------------- ------------------------
 if (paw_opt==3.or.paw_opt==4) then ! Use Sij, overlap contribution
   !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) &
   !$OMP& MAP(to:sij,gx,gxfac_sij) &
   !$OMP& PRIVATE(idat,ispinor) &
   !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
   do idat=1,ndat
   do ispinor=1,nspinor
     !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(ia,jlmn,j0lmn,jjlmn,jlm,ilmn,ilm,i0lmn,ijlmn,ii)
     do ia=1,nincat
       do jlmn=1,nlmn
         j0lmn=jlmn*(jlmn-1)/2
         jjlmn=j0lmn+jlmn
         do ii=1,cplex
           gxfac_sij(ii,ibeg+jlmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor)= &
             gxfac_sij(ii,ibeg+jlmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor) &
             + sij(jjlmn) * gx(ii,ibeg+jlmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor)
         end do
         do ilmn=1,jlmn-1
           ijlmn=j0lmn+ilmn
           do ii=1,cplex
             gxfac_sij(ii,ibeg+jlmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor)= &
               gxfac_sij(ii,ibeg+jlmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor) &
               + sij(ijlmn) * gx(ii,ibeg+ilmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor)
           end do
         end do
         if(jlmn<nlmn) then
           do ilmn=jlmn+1,nlmn
             i0lmn=ilmn*(ilmn-1)/2
             ijlmn=i0lmn+jlmn
             do ii=1,cplex
               gxfac_sij(ii,ibeg+jlmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor)=&
                 gxfac_sij(ii,ibeg+jlmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor) &
                 + sij(ijlmn) * gx(ii,ibeg+ilmn+(ia-1)*nlmn,ispinor+(idat-1)*nspinor)
             end do
           end do
         end if
       end do
     end do
   end do
   end do
 end if

!Accumulate dgxdtfac related to overlap (Sij) (PAW)
!-------------------------------------------------------------------
 if (optder>=1.and.(paw_opt==3.or.paw_opt==4)) then ! Use Sij, overlap contribution
   !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) &
   !$OMP& MAP(to:sij,dgxdt,dgxdtfac_sij) &
   !$OMP& PRIVATE(idat,ispinor) &
   !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
   do idat=1,ndat
   do ispinor=1,nspinor
     !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(ia,jlmn,j0lmn,jjlmn,jlm,ilmn,ilm,i0lmn,ijlmn,ii)
     do ia=1,nincat
       do jlmn=1,nlmn
         j0lmn=jlmn*(jlmn-1)/2
         jjlmn=j0lmn+jlmn
         do mu=1,ndgxdtfac
           do ii=1,cplex
             dgxdtfac_sij(ii,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
             &    dgxdtfac_sij(ii,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)&
             &    + sij(jjlmn) * dgxdt(ii,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
           end do
         end do
         do ilmn=1,jlmn-1
           ijlmn=j0lmn+ilmn
           do mu=1,ndgxdtfac
             do ii=1,cplex
               dgxdtfac_sij(ii,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
               &    dgxdtfac_sij(ii,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)&
               &    + sij(ijlmn) * dgxdt(ii,mu,ilmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
             end do
           end do
         end do
         if(jlmn<nlmn) then
           do ilmn=jlmn+1,nlmn
             i0lmn=ilmn*(ilmn-1)/2
             ijlmn=i0lmn+jlmn
             do mu=1,ndgxdtfac
               do ii=1,cplex
                 dgxdtfac_sij(ii,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)=&
                 &    dgxdtfac_sij(ii,mu,jlmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)&
                 &    + sij(ijlmn) * dgxdt(ii,mu,ilmn+(ia-1)*nlmn+ibeg,ispinor+(idat-1)*nspinor)
               end do
             end do
           end do
         end if
       end do
     end do
   end do
   end do
 end if

#else

 ABI_UNUSED((/cplex,cplex_enl,cplex_fac,dimenl1,dimenl2,dimekbq,iatm,itypat,gpu_option/))
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
