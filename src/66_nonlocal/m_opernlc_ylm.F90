!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_opernlc_ylm
!! NAME
!!  m_opernlc_ylm
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (MT)
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

module m_opernlc_ylm

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi

 use defs_abitypes, only : MPI_type

 implicit none

 private
!!***

 public :: opernlc_ylm
!!***

contains
!!***

!!****f* ABINIT/opernlc_ylm
!! NAME
!! opernlc_ylm
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
!!  mpi_enreg=informations about MPI parallelization
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
!!      m_gemm_nonlop,nonlop_ylm
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

subroutine opernlc_ylm(atindx1,cplex,cplex_dgxdt,cplex_d2gxdt,cplex_enl,cplex_fac,&
&          dgxdt,dgxdtfac,dgxdtfac_sij,d2gxdt,d2gxdtfac,d2gxdtfac_sij,dimenl1,dimenl2,dimekbq,enl,&
&          gx,gxfac,gxfac_sij,iatm,indlmn,itypat,lambda,mpi_enreg,natom,ndgxdt,ndgxdtfac,&
&          nd2gxdt,nd2gxdtfac,nincat,nlmn,nspinor,nspinortot,optder,paw_opt,sij)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,cplex_enl,cplex_fac,dimenl1,dimenl2,dimekbq,iatm,itypat
 integer,intent(in) :: natom,ndgxdt,ndgxdtfac,nd2gxdt,nd2gxdtfac,nincat,nspinor,nspinortot,optder,paw_opt
 integer,intent(inout) :: nlmn
 real(dp) :: lambda
 type(MPI_type) , intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom),indlmn(6,nlmn),cplex_dgxdt(ndgxdt),cplex_d2gxdt(nd2gxdt)
 real(dp),intent(in) :: dgxdt(cplex,ndgxdt,nlmn,nincat,nspinor)
 real(dp),intent(in) :: d2gxdt(cplex,nd2gxdt,nlmn,nincat,nspinor)
 real(dp),intent(in),target :: enl(dimenl1,dimenl2,nspinortot**2,dimekbq)
 real(dp),intent(inout) :: gx(cplex,nlmn,nincat,nspinor)
 real(dp),intent(in) :: sij(((paw_opt+1)/3)*nlmn*(nlmn+1)/2)
 real(dp),intent(out),target :: dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)
 real(dp),intent(out) :: dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat,nspinor*(paw_opt/3))
 real(dp),intent(out),target :: d2gxdtfac(cplex_fac,nd2gxdtfac,nlmn,nincat,nspinor)
 real(dp),intent(out) :: d2gxdtfac_sij(cplex,nd2gxdtfac,nlmn,nincat,nspinor*(paw_opt/3))
 real(dp),intent(out),target :: gxfac(cplex_fac,nlmn,nincat,nspinor)
 real(dp),intent(out) :: gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))

!Local variables-------------------------------
!Arrays
!scalars
 integer :: cplex_,ia,ierr,ijlmn,ijspin,ilm,ilmn,i0lmn,iln,index_enl,iphase,ispinor,ispinor_index
 integer :: j0lmn,jilmn,jispin,jjlmn,jlm,jlmn,jspinor,jspinor_index,mu,shift
 real(dp) :: sijr
!arrays
 real(dp) :: enl_(2),gxfi(2),gxi(cplex),gxj(cplex)
 real(dp),allocatable :: d2gxdtfac_offdiag(:,:,:,:,:),dgxdtfac_offdiag(:,:,:,:,:)
 real(dp),allocatable :: gxfac_offdiag(:,:,:,:),gxfj(:,:)
 real(dp),pointer :: d2gxdtfac_(:,:,:,:,:),dgxdtfac_(:,:,:,:,:),gxfac_(:,:,:,:)
 real(dp),pointer :: enl_ptr(:,:,:)

! *************************************************************************

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
   ABI_ALLOCATE(gxfac_,(cplex_fac,nlmn,nincat,nspinor))
   if (optder>=1) then
    ABI_ALLOCATE(dgxdtfac_,(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor))
   end if
   if (optder>=2) then
    ABI_ALLOCATE(d2gxdtfac_,(cplex_fac,nd2gxdtfac,nlmn,nincat,nspinor))
   end if
  end if
  gxfac_=zero
  if (optder>=1) dgxdtfac_=zero
  if (optder>=2) d2gxdtfac_=zero
  enl_ptr => enl(:,:,:,iphase)


!Accumulate gxfac related to non-local operator (Norm-conserving)
!-------------------------------------------------------------------
  if (paw_opt==0) then
   !Enl is E(Kleinman-Bylander)
   ABI_CHECK(cplex_enl/=2,"BUG: invalid cplex_enl=2!")
   ABI_CHECK(cplex_fac==cplex,"BUG: invalid cplex_fac/=cplex!")
!$OMP PARALLEL &
!$OMP PRIVATE(ispinor,ispinor_index,ia,ilmn,iln,enl_)
!$OMP DO COLLAPSE(3)
   do ispinor=1,nspinor
     do ia=1,nincat
       do ilmn=1,nlmn
         ispinor_index=ispinor+shift
         iln=indlmn(5,ilmn)
         enl_(1)=enl_ptr(iln,itypat,ispinor_index)
         gxfac_(1:cplex,ilmn,ia,ispinor)=enl_(1)*gx(1:cplex,ilmn,ia,ispinor)
       end do
     end do
   end do
!$OMP END DO
!$OMP END PARALLEL
  end if

!Accumulate gxfac related to nonlocal operator (PAW)
!-------------------------------------------------------------------
  if (paw_opt==1.or.paw_opt==2.or.paw_opt==4) then
   !Enl is psp strength Dij or (Dij-lambda.Sij)

!  === Diagonal term(s) (up-up, down-down)

!  1-Enl is real
   if (cplex_enl==1) then
!$OMP PARALLEL &
!$OMP PRIVATE(ispinor,ispinor_index,ia,index_enl), &
!$OMP PRIVATE(jlmn,j0lmn,jjlmn,enl_,gxj,ilmn,ijlmn,gxi)
!$OMP DO COLLAPSE(3)
     do ispinor=1,nspinor
       do ia=1,nincat
         do jlmn=1,nlmn
           ispinor_index=ispinor+shift
           index_enl=atindx1(iatm+ia)
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1)=enl_ptr(jjlmn,index_enl,ispinor_index)
           if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
           gxj(1:cplex)=gx(1:cplex,jlmn,ia,ispinor)
           gxfac_(1:cplex,jlmn,ia,ispinor)=gxfac_(1:cplex,jlmn,ia,ispinor)+enl_(1)*gxj(1:cplex)
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1)=enl_ptr(ijlmn,index_enl,ispinor_index)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
             gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
             gxfac_(1:cplex,jlmn,ia,ispinor)=gxfac_(1:cplex,jlmn,ia,ispinor)+enl_(1)*gxi(1:cplex)
#if !defined HAVE_OPENMP
             gxfac_(1:cplex,ilmn,ia,ispinor)=gxfac_(1:cplex,ilmn,ia,ispinor)+enl_(1)*gxj(1:cplex)
#endif
           end do
#if defined HAVE_OPENMP
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=(ilmn*(ilmn-1)/2)
               ijlmn=i0lmn+jlmn
               enl_(1)=enl_ptr(ijlmn,index_enl,ispinor_index)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
               gxfac_(1:cplex,jlmn,ia,ispinor)=gxfac_(1:cplex,jlmn,ia,ispinor)+enl_(1)*gxi(1:cplex)
             end do
           end if
#endif
         end do
       end do
     end do
!$OMP END DO
!$OMP END PARALLEL

!    2-Enl is complex  ===== D^ss'_ij=D^s's_ji^*
   else
     ABI_CHECK(cplex_fac==cplex_enl,"BUG: invalid cplex_fac/=cplex_enl!")

     if (nspinortot==1) then ! -------------> NO SPINORS

!$OMP PARALLEL &
!$OMP PRIVATE(ia,index_enl,jlmn,j0lmn,jjlmn,enl_,gxj,ilmn,ijlmn,gxi)
       do ia=1,nincat
         index_enl=atindx1(iatm+ia)
!$OMP DO
         do jlmn=1,nlmn
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1)=enl_ptr(2*jjlmn-1,index_enl,1)
           if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
           gxj(1:cplex)=gx(1:cplex,jlmn,ia,1)
           gxfac_(1,jlmn,ia,1)=gxfac_(1,jlmn,ia,1)+enl_(1)*gxj(1)
           if (cplex==2) gxfac_(2,jlmn,ia,1)=gxfac_(2,jlmn,ia,1)+enl_(1)*gxj(2)
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,1)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
             gxi(1:cplex)=gx(1:cplex,ilmn,ia,1)
             gxfac_(1,jlmn,ia,1)=gxfac_(1,jlmn,ia,1)+enl_(1)*gxi(1)
             gxfac_(2,jlmn,ia,1)=gxfac_(2,jlmn,ia,1)-enl_(2)*gxi(1)
#if !defined HAVE_OPENMP
             gxfac_(1,ilmn,ia,1)=gxfac_(1,ilmn,ia,1)+enl_(1)*gxj(1)
             gxfac_(2,ilmn,ia,1)=gxfac_(2,ilmn,ia,1)+enl_(2)*gxj(1)
#endif
             if (cplex==2) then
               gxfac_(1,jlmn,ia,1)=gxfac_(1,jlmn,ia,1)+enl_(2)*gxi(2)
               gxfac_(2,jlmn,ia,1)=gxfac_(2,jlmn,ia,1)+enl_(1)*gxi(2)
#if !defined HAVE_OPENMP
               gxfac_(1,ilmn,ia,1)=gxfac_(1,ilmn,ia,1)-enl_(2)*gxj(2)
               gxfac_(2,ilmn,ia,1)=gxfac_(2,ilmn,ia,1)+enl_(1)*gxj(2)
#endif
             end if
           end do
#if defined HAVE_OPENMP
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=ilmn*(ilmn-1)/2
               ijlmn=i0lmn+jlmn
               enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,1)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               gxi(1:cplex)=gx(1:cplex,ilmn,ia,1)
               gxfac_(1,jlmn,ia,1)=gxfac_(1,jlmn,ia,1)+enl_(1)*gxi(1)
               gxfac_(2,jlmn,ia,1)=gxfac_(2,jlmn,ia,1)+enl_(2)*gxi(1)
               if (cplex==2) then
                 gxfac_(1,jlmn,ia,1)=gxfac_(1,jlmn,ia,1)-enl_(2)*gxi(2)
                 gxfac_(2,jlmn,ia,1)=gxfac_(2,jlmn,ia,1)+enl_(1)*gxi(2)
               end if
             end do
           end if
#endif
         end do
!$OMP END DO
       end do
!$OMP END PARALLEL

     else ! -------------> SPINORIAL CASE

!  === Diagonal term(s) (up-up, down-down)

!$OMP PARALLEL &
!$OMP PRIVATE(ispinor,ispinor_index,ia,index_enl), &
!$OMP PRIVATE(jlmn,j0lmn,jjlmn,enl_,gxj,ilmn,ijlmn,gxi)
       do ispinor=1,nspinor
         ispinor_index=ispinor+shift
         do ia=1,nincat
           index_enl=atindx1(iatm+ia)
!$OMP DO
           do jlmn=1,nlmn
             j0lmn=jlmn*(jlmn-1)/2
             jjlmn=j0lmn+jlmn
             enl_(1)=enl_ptr(2*jjlmn-1,index_enl,ispinor_index)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
             gxj(1:cplex)=gx(1:cplex,jlmn,ia,ispinor)
             gxfac_(1,jlmn,ia,ispinor)=gxfac_(1,jlmn,ia,ispinor)+enl_(1)*gxj(1)
             if (cplex==2)  gxfac_(2,jlmn,ia,ispinor)=gxfac_(2,jlmn,ia,ispinor)+enl_(1)*gxj(2)
             do ilmn=1,jlmn-1
               ijlmn=j0lmn+ilmn
               enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,ispinor_index)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
               gxfac_(1,jlmn,ia,ispinor)=gxfac_(1,jlmn,ia,ispinor)+enl_(1)*gxi(1)
               gxfac_(2,jlmn,ia,ispinor)=gxfac_(2,jlmn,ia,ispinor)-enl_(2)*gxi(1)
#if !defined HAVE_OPENMP
               gxfac_(1,ilmn,ia,ispinor)=gxfac_(1,ilmn,ia,ispinor)+enl_(1)*gxj(1)
               gxfac_(2,ilmn,ia,ispinor)=gxfac_(2,ilmn,ia,ispinor)+enl_(2)*gxj(1)
#endif
               if (cplex==2) then
                 gxfac_(1,jlmn,ia,ispinor)=gxfac_(1,jlmn,ia,ispinor)+enl_(2)*gxi(2)
                 gxfac_(2,jlmn,ia,ispinor)=gxfac_(2,jlmn,ia,ispinor)+enl_(1)*gxi(2)
#if !defined HAVE_OPENMP
                 gxfac_(1,ilmn,ia,ispinor)=gxfac_(1,ilmn,ia,ispinor)-enl_(2)*gxj(2)
                 gxfac_(2,ilmn,ia,ispinor)=gxfac_(2,ilmn,ia,ispinor)+enl_(1)*gxj(2)
#endif
               end if
             end do
#if defined HAVE_OPENMP
             if(jlmn<nlmn) then
               do ilmn=jlmn+1,nlmn
                 i0lmn=ilmn*(ilmn-1)/2
                 ijlmn=i0lmn+jlmn
                 enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,ispinor_index)
                 if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
                 gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
                 gxfac_(1,jlmn,ia,ispinor)=gxfac_(1,jlmn,ia,ispinor)+enl_(1)*gxi(1)
                 gxfac_(2,jlmn,ia,ispinor)=gxfac_(2,jlmn,ia,ispinor)+enl_(2)*gxi(1)
                 if (cplex==2) then
                   gxfac_(1,jlmn,ia,ispinor)=gxfac_(1,jlmn,ia,ispinor)-enl_(2)*gxi(2)
                   gxfac_(2,jlmn,ia,ispinor)=gxfac_(2,jlmn,ia,ispinor)+enl_(1)*gxi(2)
                 end if
               end do
             end if
#endif
           end do
!$OMP END DO
         end do
       end do
!$OMP END PARALLEL
     end if !nspinortot
   end if !complex_enl

!  === Off-diagonal term(s) (up-down, down-up)

!  --- No parallelization over spinors ---
   if (nspinortot==2.and.nspinor==nspinortot) then
     ABI_CHECK(cplex_enl==2,"BUG: invalid cplex_enl/=2!")
     ABI_CHECK(cplex_fac==cplex,"BUG: invalid cplex_fac/=cplex)!")
!$OMP PARALLEL &
!$OMP PRIVATE(ispinor,jspinor,ia,index_enl), &
!$OMP PRIVATE(jlmn,j0lmn,jjlmn,enl_,gxi,gxj,ilmn,ijlmn)
     do ispinor=1,nspinortot
       jspinor=3-ispinor
       do ia=1,nincat
         index_enl=atindx1(iatm+ia)
!$OMP DO
         do jlmn=1,nlmn
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1:2)=enl_ptr(2*jjlmn-1:2*jjlmn,index_enl,2+ispinor )
           gxi(1:cplex)=gx(1:cplex,jlmn,ia,ispinor)
           gxfac_(1,jlmn,ia,jspinor)=gxfac_(1,jlmn,ia,jspinor)+enl_(1)*gxi(1)
           gxfac_(2,jlmn,ia,jspinor)=gxfac_(2,jlmn,ia,jspinor)-enl_(2)*gxi(1)
           if (cplex==2) then
             gxfac_(1,jlmn,ia,jspinor)=gxfac_(1,jlmn,ia,jspinor)+enl_(2)*gxi(2)
             gxfac_(2,jlmn,ia,jspinor)=gxfac_(2,jlmn,ia,jspinor)+enl_(1)*gxi(2)
           end if
#if !defined HAVE_OPENMP
           gxj(1:cplex)=gx(1:cplex,jlmn,ia,jspinor)
#endif
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor)
             gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
             gxfac_(1,jlmn,ia,jspinor)=gxfac_(1,jlmn,ia,jspinor)+enl_(1)*gxi(1)
             gxfac_(2,jlmn,ia,jspinor)=gxfac_(2,jlmn,ia,jspinor)-enl_(2)*gxi(1)
#if !defined HAVE_OPENMP
             gxfac_(1,ilmn,ia,ispinor)=gxfac_(1,ilmn,ia,ispinor)+enl_(1)*gxj(1)
             gxfac_(2,ilmn,ia,ispinor)=gxfac_(2,ilmn,ia,ispinor)+enl_(2)*gxj(1)
#endif
             if (cplex==2) then
               gxfac_(1,jlmn,ia,jspinor)=gxfac_(1,jlmn,ia,jspinor)+enl_(2)*gxi(2)
               gxfac_(2,jlmn,ia,jspinor)=gxfac_(2,jlmn,ia,jspinor)+enl_(1)*gxi(2)
#if !defined HAVE_OPENMP
               gxfac_(1,ilmn,ia,ispinor)=gxfac_(1,ilmn,ia,ispinor)-enl_(2)*gxj(2)
               gxfac_(2,ilmn,ia,ispinor)=gxfac_(2,ilmn,ia,ispinor)+enl_(1)*gxj(2)
#endif
             end if
           end do
#if defined HAVE_OPENMP
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=ilmn*(ilmn-1)/2
               ijlmn=i0lmn+jlmn
               enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor)
               gxi(1:cplex)=gx(1:cplex,ilmn,ia,jspinor)
               gxfac_(1,jlmn,ia,ispinor)=gxfac_(1,jlmn,ia,ispinor)+enl_(1)*gxi(1)
               gxfac_(2,jlmn,ia,ispinor)=gxfac_(2,jlmn,ia,ispinor)+enl_(2)*gxi(1)
               if (cplex==2) then
                 gxfac_(1,jlmn,ia,ispinor)=gxfac_(1,jlmn,ia,ispinor)-enl_(2)*gxi(2)
                 gxfac_(2,jlmn,ia,ispinor)=gxfac_(2,jlmn,ia,ispinor)+enl_(1)*gxi(2)
               end if
             end do
           end if
#endif
         end do
!$OMP END DO
       end do
     end do
!$OMP END PARALLEL

!    --- Parallelization over spinors ---
   else if (nspinortot==2.and.nspinor/=nspinortot) then
     ABI_CHECK(cplex_enl==2,"BUG: invalid cplex_enl/=2!")
     ABI_CHECK(cplex_fac==2,"BUG: invalid cplex_fac/=2!")
     ABI_ALLOCATE(gxfac_offdiag,(cplex_fac,nlmn,nincat,nspinortot))
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
!$OMP PARALLEL &
!$OMP PRIVATE(ia,index_enl,jlmn,j0lmn,ilmn,i0lmn,ijlmn,enl_,jilmn,gxi)
     do ia=1,nincat
       index_enl=atindx1(iatm+ia)
!$OMP DO
       do jlmn=1,nlmn
         j0lmn=jlmn*(jlmn-1)/2
         do ilmn=1,nlmn
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
           gxi(1:cplex)=gx(1:cplex,ilmn,ia,1)
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
!$OMP END DO
     end do !iat
!$OMP END PARALLEL
     call xmpi_sum(gxfac_offdiag,mpi_enreg%comm_spinor,ierr)
     gxfac_(:,:,:,1)=gxfac_(:,:,:,1)+gxfac_offdiag(:,:,:,ispinor_index)
     ABI_DEALLOCATE(gxfac_offdiag)
   end if

  end if !paw_opt

!Accumulate dgxdtfac related to nonlocal operator (Norm-conserving)
!-------------------------------------------------------------------
  if (optder>=1.and.paw_opt==0) then
   !Enl is E(Kleinman-Bylander)
   ABI_CHECK(cplex_enl==1,"BUG: invalid cplex_enl/=1!")
   ABI_CHECK(cplex_fac==cplex,"BUG: invalid cplex_fac/=cplex!")
!$OMP PARALLEL &
!$OMP PRIVATE(ispinor,ispinor_index,ia,ilmn,iln,enl_,mu)
   do ispinor=1,nspinor
     ispinor_index = ispinor + shift
     do ia=1,nincat
!$OMP DO
       do ilmn=1,nlmn
         iln=indlmn(5,ilmn)
         enl_(1)=enl_ptr(iln,itypat,ispinor_index)
         do mu=1,ndgxdtfac
           dgxdtfac_(1:cplex,mu,ilmn,ia,ispinor)=enl_(1)*dgxdt(1:cplex,mu,ilmn,ia,ispinor)
         end do
       end do
!$OMP END DO
     end do
   end do
!$OMP END PARALLEL
  end if

!Accumulate dgxdtfac related to nonlocal operator (PAW)
!-------------------------------------------------------------------
  if (optder>=1.and.(paw_opt==1.or.paw_opt==2.or.paw_opt==4)) then
   !Enl is psp strength Dij or (Dij-lambda.Sij)

!  === Diagonal term(s) (up-up, down-down)

!  1-Enl is real
   if (cplex_enl==1) then
!$OMP PARALLEL &
!$OMP PRIVATE(ispinor,ispinor_index,ia,index_enl,jlmn,j0lmn,jjlmn,enl_,mu,gxfj,ilmn,ijlmn,gxfi)
     ABI_ALLOCATE(gxfj,(cplex,ndgxdtfac))
     do ispinor=1,nspinor
       ispinor_index=ispinor+shift
       do ia=1,nincat
         index_enl=atindx1(iatm+ia)
!$OMP DO
         do jlmn=1,nlmn
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1)=enl_ptr(jjlmn,index_enl,ispinor_index)
           if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
           do mu=1,ndgxdtfac
             gxfj(1:cplex,mu)=dgxdt(1:cplex,mu,jlmn,ia,ispinor)
             dgxdtfac_(1:cplex,mu,jlmn,ia,ispinor)=dgxdtfac_(1:cplex,mu,jlmn,ia,ispinor)+enl_(1)*gxfj(1:cplex,mu)
           end do
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1)=enl_ptr(ijlmn,index_enl,ispinor_index)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
             do mu=1,ndgxdtfac
               gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
               dgxdtfac_(1:cplex,mu,jlmn,ia,ispinor)=dgxdtfac_(1:cplex,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1:cplex)
#if !defined HAVE_OPENMP
               dgxdtfac_(1:cplex,mu,ilmn,ia,ispinor)=dgxdtfac_(1:cplex,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(1:cplex,mu)
#endif
             end do
           end do
#if defined HAVE_OPENMP
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=ilmn*(ilmn-1)/2
               ijlmn=i0lmn+jlmn
               enl_(1)=enl_ptr(ijlmn,index_enl,ispinor_index)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               do mu=1,ndgxdtfac
                 gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
                 dgxdtfac_(1:cplex,mu,jlmn,ia,ispinor)=dgxdtfac_(1:cplex,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1:cplex)
               end do
             end do
           end if
#endif
         end do
!$OMP END DO
       end do
     end do
     ABI_DEALLOCATE(gxfj)
!$OMP END PARALLEL

!    2-Enl is complex  ===== D^ss'_ij=D^s's_ji^*
   else
     ABI_CHECK(cplex_fac==cplex_enl,"BUG: invalid cplex_fac/=cplex_enl!")

     if (nspinortot==1) then ! -------------> NO SPINORS

!$OMP PARALLEL &
!$OMP PRIVATE(ia,index_enl,jlmn,j0lmn,jjlmn,enl_,mu,gxfj,ilmn,ijlmn,gxfi)
       ABI_ALLOCATE(gxfj,(cplex,ndgxdtfac))
       do ia=1,nincat
         index_enl=atindx1(iatm+ia)
!$OMP DO
         do jlmn=1,nlmn
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1)=enl_ptr(2*jjlmn-1,index_enl,1)
           if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
           do mu=1,ndgxdtfac
             if(cplex_dgxdt(mu)==2)then
               cplex_ = 2 ; gxfj(1,mu) = zero ; gxfj(2,mu) = dgxdt(1,mu,jlmn,ia,1)
             else
               cplex_ = cplex ; gxfj(1:cplex,mu)=dgxdt(1:cplex,mu,jlmn,ia,1)
             end if
             dgxdtfac_(1,mu,jlmn,ia,1)=dgxdtfac_(1,mu,jlmn,ia,1)+enl_(1)*gxfj(1,mu)
             if (cplex_==2) dgxdtfac_(2,mu,jlmn,ia,1)=dgxdtfac_(2,mu,jlmn,ia,1)+enl_(1)*gxfj(2,mu)
           end do
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,1)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
             do mu=1,ndgxdtfac
               if(cplex_dgxdt(mu)==2)then
                 cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = dgxdt(1,mu,ilmn,ia,1)
               else
                 cplex_ = cplex ; gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,1)
               end if
               dgxdtfac_(1,mu,jlmn,ia,1)=dgxdtfac_(1,mu,jlmn,ia,1)+enl_(1)*gxfi(1)
               dgxdtfac_(2,mu,jlmn,ia,1)=dgxdtfac_(2,mu,jlmn,ia,1)-enl_(2)*gxfi(1)
#if !defined HAVE_OPENMP
               dgxdtfac_(1,mu,ilmn,ia,1)=dgxdtfac_(1,mu,ilmn,ia,1)+enl_(1)*gxfj(1,mu)
               dgxdtfac_(2,mu,ilmn,ia,1)=dgxdtfac_(2,mu,ilmn,ia,1)+enl_(2)*gxfj(1,mu)
#endif
               if (cplex_==2) then
                 dgxdtfac_(1,mu,jlmn,ia,1)=dgxdtfac_(1,mu,jlmn,ia,1)+enl_(2)*gxfi(2)
                 dgxdtfac_(2,mu,jlmn,ia,1)=dgxdtfac_(2,mu,jlmn,ia,1)+enl_(1)*gxfi(2)
#if !defined HAVE_OPENMP
                 dgxdtfac_(1,mu,ilmn,ia,1)=dgxdtfac_(1,mu,ilmn,ia,1)-enl_(2)*gxfj(2,mu)
                 dgxdtfac_(2,mu,ilmn,ia,1)=dgxdtfac_(2,mu,ilmn,ia,1)+enl_(1)*gxfj(2,mu)
#endif
               end if
             end do
           end do
#if defined HAVE_OPENMP
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=ilmn*(ilmn-1)/2
               ijlmn=i0lmn+jlmn
               enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,1)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               do mu=1,ndgxdtfac
                 if(cplex_dgxdt(mu)==2)then
                   cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = dgxdt(1,mu,ilmn,ia,1)
                 else
                   cplex_ = cplex ; gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,1)
                 end if
                 dgxdtfac_(1,mu,jlmn,ia,1)=dgxdtfac_(1,mu,jlmn,ia,1)+enl_(1)*gxfi(1)
                 dgxdtfac_(2,mu,jlmn,ia,1)=dgxdtfac_(2,mu,jlmn,ia,1)+enl_(2)*gxfi(1)
                 if (cplex_==2) then
                   dgxdtfac_(1,mu,jlmn,ia,1)=dgxdtfac_(1,mu,jlmn,ia,1)-enl_(2)*gxfi(2)
                   dgxdtfac_(2,mu,jlmn,ia,1)=dgxdtfac_(2,mu,jlmn,ia,1)+enl_(1)*gxfi(2)
                 end if
               end do
             end do
           end if
#endif
         end do
!$OMP END DO
       end do
       ABI_DEALLOCATE(gxfj)
!$OMP END PARALLEL

     else ! -------------> SPINORIAL CASE

!  === Diagonal term(s) (up-up, down-down)

!$OMP PARALLEL &
!$OMP PRIVATE(ispinor,ispinor_index,ia,index_enl), &
!$OMP PRIVATE(jlmn,j0lmn,jjlmn,enl_,mu,gxfj,ilmn,ijlmn,gxfi)
       ABI_ALLOCATE(gxfj,(cplex,ndgxdtfac))
       do ispinor=1,nspinor
         ispinor_index = ispinor + shift
         do ia=1,nincat
           index_enl=atindx1(iatm+ia)
!$OMP DO
           do jlmn=1,nlmn
             j0lmn=jlmn*(jlmn-1)/2
             jjlmn=j0lmn+jlmn
             enl_(1)=enl_ptr(2*jjlmn-1,index_enl,ispinor_index)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
             do mu=1,ndgxdtfac
               if(cplex_dgxdt(mu)==2)then
                 cplex_ = 2 ; gxfj(1,mu) = zero ; gxfj(2,mu) = dgxdt(1,mu,jlmn,ia,ispinor)
               else
                 cplex_ = cplex ; gxfj(1:cplex,mu)=dgxdt(1:cplex,mu,jlmn,ia,ispinor)
               end if
               dgxdtfac_(1,mu,jlmn,ia,ispinor)=dgxdtfac_(1,mu,jlmn,ia,ispinor)+enl_(1)*gxfj(1,mu)
               if (cplex_==2) dgxdtfac_(2,mu,jlmn,ia,ispinor)=dgxdtfac_(2,mu,jlmn,ia,ispinor)+enl_(1)*gxfj(2,mu)
             end do
             do ilmn=1,jlmn-1
               ijlmn=j0lmn+ilmn
               enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,ispinor_index)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               do mu=1,ndgxdtfac
                 if(cplex_dgxdt(mu)==2)then
                   cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = dgxdt(1,mu,ilmn,ia,ispinor)
                 else
                   cplex_ = cplex ; gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
                 end if
                 dgxdtfac_(1,mu,jlmn,ia,ispinor)=dgxdtfac_(1,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1)
                 dgxdtfac_(2,mu,jlmn,ia,ispinor)=dgxdtfac_(2,mu,jlmn,ia,ispinor)-enl_(2)*gxfi(1)
#if !defined HAVE_OPENMP
                 dgxdtfac_(1,mu,ilmn,ia,ispinor)=dgxdtfac_(1,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(1,mu)
                 dgxdtfac_(2,mu,ilmn,ia,ispinor)=dgxdtfac_(2,mu,ilmn,ia,ispinor)+enl_(2)*gxfj(1,mu)
#endif
                 if (cplex_==2) then
                   dgxdtfac_(1,mu,jlmn,ia,ispinor)=dgxdtfac_(1,mu,jlmn,ia,ispinor)+enl_(2)*gxfi(2)
                   dgxdtfac_(2,mu,jlmn,ia,ispinor)=dgxdtfac_(2,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(2)
#if !defined HAVE_OPENMP
                   dgxdtfac_(1,mu,ilmn,ia,ispinor)=dgxdtfac_(1,mu,ilmn,ia,ispinor)-enl_(2)*gxfj(2,mu)
                   dgxdtfac_(2,mu,ilmn,ia,ispinor)=dgxdtfac_(2,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(2,mu)
#endif
                 end if
               end do
             end do
#if defined HAVE_OPENMP
             if(jlmn<nlmn) then
               do ilmn=jlmn+1,nlmn
                 i0lmn=ilmn*(ilmn-1)/2
                 ijlmn=i0lmn+jlmn
                 enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,ispinor_index)
                 if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
                 do mu=1,ndgxdtfac
                   if(cplex_dgxdt(mu)==2)then
                     cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = dgxdt(1,mu,ilmn,ia,ispinor)
                   else
                     cplex_ = cplex ; gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
                   end if
                   dgxdtfac_(1,mu,jlmn,ia,ispinor)=dgxdtfac_(1,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1)
                   dgxdtfac_(2,mu,jlmn,ia,ispinor)=dgxdtfac_(2,mu,jlmn,ia,ispinor)+enl_(2)*gxfi(1)
                   if (cplex_==2) then
                     dgxdtfac_(1,mu,jlmn,ia,ispinor)=dgxdtfac_(1,mu,jlmn,ia,ispinor)-enl_(2)*gxfi(2)
                     dgxdtfac_(2,mu,jlmn,ia,ispinor)=dgxdtfac_(2,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(2)
                   end if
                 end do
               end do
             end if
#endif
           end do
!$OMP END DO
         end do
       end do
       ABI_DEALLOCATE(gxfj)
!$OMP END PARALLEL
     end if !nspinortot
   end if !complex

!  === Off-diagonal term(s) (up-down, down-up)

!  --- No parallelization over spinors ---
   if (nspinortot==2.and.nspinor==nspinortot) then
     ABI_CHECK(cplex_enl==2,"BUG: invalid cplex_enl/=2!")
     ABI_CHECK(cplex_fac==2,"BUG: invalid cplex_fac/=2!")
!$OMP PARALLEL &
!$OMP PRIVATE(ispinor,jspinor,ia,index_enl), &
!$OMP PRIVATE(jlmn,j0lmn,jjlmn,enl_,mu,gxfi,gxfj,ilmn,ijlmn)
     ABI_ALLOCATE(gxfj,(cplex,ndgxdtfac))
     do ispinor=1,nspinor
       jspinor=3-ispinor
       do ia=1,nincat
         index_enl=atindx1(iatm+ia)
!$OMP DO
         do jlmn=1,nlmn
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1:2)=enl_ptr(2*jjlmn-1:2*jjlmn,index_enl,2+ispinor)
           do mu=1,ndgxdtfac
             if(cplex_dgxdt(mu)==2)then
               cplex_ = 2 ;
               gxfi(1)    = zero ; gxfi(2)    = dgxdt(1,mu,jlmn,ia,ispinor)
               gxfj(1,mu) = zero ; gxfj(2,mu) = dgxdt(1,mu,jlmn,ia,jspinor)
             else
               cplex_ = cplex ;
               gxfi(1:cplex)   =dgxdt(1:cplex,mu,jlmn,ia,ispinor)
               gxfj(1:cplex,mu)=dgxdt(1:cplex,mu,jlmn,ia,jspinor)
             end if
             dgxdtfac_(1,mu,jlmn,ia,jspinor)=dgxdtfac_(1,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(1)
             dgxdtfac_(2,mu,jlmn,ia,jspinor)=dgxdtfac_(2,mu,jlmn,ia,jspinor)-enl_(2)*gxfi(1)
             if (cplex_==2) then
               dgxdtfac_(1,mu,jlmn,ia,jspinor)=dgxdtfac_(1,mu,jlmn,ia,jspinor)+enl_(2)*gxfi(2)
               dgxdtfac_(2,mu,jlmn,ia,jspinor)=dgxdtfac_(2,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(2)
             end if
           end do
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor)
             do mu=1,ndgxdtfac
               if(cplex_dgxdt(mu)==2)then
                 cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = dgxdt(1,mu,ilmn,ia,ispinor)
               else
                 cplex_ = cplex ; gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
               end if
               dgxdtfac_(1,mu,jlmn,ia,jspinor)=dgxdtfac_(1,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(1)
               dgxdtfac_(2,mu,jlmn,ia,jspinor)=dgxdtfac_(2,mu,jlmn,ia,jspinor)-enl_(2)*gxfi(1)
#if !defined HAVE_OPENMP
               dgxdtfac_(1,mu,ilmn,ia,ispinor)=dgxdtfac_(1,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(1,mu)
               dgxdtfac_(2,mu,ilmn,ia,ispinor)=dgxdtfac_(2,mu,ilmn,ia,ispinor)+enl_(2)*gxfj(1,mu)
#endif
               if (cplex_==2) then
                 dgxdtfac_(1,mu,jlmn,ia,jspinor)=dgxdtfac_(1,mu,jlmn,ia,jspinor)+enl_(2)*gxfi(2)
                 dgxdtfac_(2,mu,jlmn,ia,jspinor)=dgxdtfac_(2,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(2)
#if !defined HAVE_OPENMP
                 dgxdtfac_(1,mu,ilmn,ia,ispinor)=dgxdtfac_(1,mu,ilmn,ia,ispinor)-enl_(2)*gxfj(2,mu)
                 dgxdtfac_(2,mu,ilmn,ia,ispinor)=dgxdtfac_(2,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(2,mu)
#endif
               end if
             end do !mu
           end do !ilmn
#if defined HAVE_OPENMP
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=ilmn*(ilmn-1)/2
               ijlmn=i0lmn+jlmn
               enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor)
               do mu=1,ndgxdtfac
                 if(cplex_dgxdt(mu)==2)then
                   cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = dgxdt(1,mu,ilmn,ia,jspinor)
                 else
                   cplex_ = cplex ; gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,jspinor)
                 end if
                 dgxdtfac_(1,mu,jlmn,ia,ispinor)=dgxdtfac_(1,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1)
                 dgxdtfac_(2,mu,jlmn,ia,ispinor)=dgxdtfac_(2,mu,jlmn,ia,ispinor)+enl_(2)*gxfi(1)
                 if (cplex_==2) then
                   dgxdtfac_(1,mu,jlmn,ia,ispinor)=dgxdtfac_(1,mu,jlmn,ia,ispinor)-enl_(2)*gxfi(2)
                   dgxdtfac_(2,mu,jlmn,ia,ispinor)=dgxdtfac_(2,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(2)
                 end if
               end do !mu
             end do !ilmn
           end if
#endif
         end do !jmln
!$OMP END DO
       end do !ia
     end do !ispinor
     ABI_DEALLOCATE(gxfj)
!$OMP END PARALLEL

!    --- Parallelization over spinors ---
   else if (nspinortot==2.and.nspinor/=nspinortot) then
     ABI_CHECK(cplex_enl==2,"BUG: opernlc_ylm: invalid cplex_enl/=2!")
     ABI_CHECK(cplex_fac==2,"BUG: opernlc_ylm: invalid cplex_fac/=2!")
     ABI_ALLOCATE(dgxdtfac_offdiag,(cplex_fac,ndgxdtfac,nlmn,nincat,nspinortot))
!$OMP PARALLEL &
!$OMP PRIVATE(ia,index_enl), &
!$OMP PRIVATE(jlmn,j0lmn,ilmn,i0lmn,ijlmn,enl_,jilmn,mu,gxfi)
!$OMP WORKSHARE
     dgxdtfac_offdiag(:,:,:,:,:)=zero
!$OMP END WORKSHARE
!$OMP SINGLE
     ispinor_index=mpi_enreg%me_spinor+1
     jspinor_index=3-ispinor_index
     if (ispinor_index==1) then
       ijspin=3;jispin=4
     else
       ijspin=4;jispin=3
     end if
!$OMP END SINGLE
     do ia=1,nincat
       index_enl=atindx1(iatm+ia)
!$OMP DO
       do jlmn=1,nlmn
         j0lmn=jlmn*(jlmn-1)/2
         do ilmn=1,nlmn
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
           do mu=1,ndgxdtfac
             if(cplex_dgxdt(mu)==2)then
               cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = dgxdt(1,mu,ilmn,ia,1)
             else
               cplex_ = cplex ; gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,1)
             end if
             dgxdtfac_offdiag(1,mu,jlmn,ia,jspinor_index)= &
&                 dgxdtfac_offdiag(1,mu,jlmn,ia,jspinor_index)+enl_(1)*gxfi(1)
             dgxdtfac_offdiag(2,mu,jlmn,ia,jspinor_index)= &
&                 dgxdtfac_offdiag(2,mu,jlmn,ia,jspinor_index)+enl_(2)*gxfi(1)
             if (cplex_==2) then
               dgxdtfac_offdiag(1,mu,jlmn,ia,jspinor_index)= &
&                   dgxdtfac_offdiag(1,mu,jlmn,ia,jspinor_index)-enl_(2)*gxfi(2)
               dgxdtfac_offdiag(2,mu,jlmn,ia,jspinor_index)= &
&                   dgxdtfac_offdiag(2,mu,jlmn,ia,jspinor_index)+enl_(1)*gxfi(2)
             end if
           end do
         end do !ilmn
       end do !jlmn
!$OMP END DO
     end do !iat
!$OMP SINGLE
     call xmpi_sum(dgxdtfac_offdiag,mpi_enreg%comm_spinor,ierr)
!$OMP END SINGLE
!$OMP WORKSHARE
     dgxdtfac_(:,:,:,:,1)=dgxdtfac_(:,:,:,:,1)+dgxdtfac_offdiag(:,:,:,:,ispinor_index)
!$OMP END WORKSHARE
!$OMP END PARALLEL
     ABI_DEALLOCATE(dgxdtfac_offdiag)
   end if !nspinortot

  end if ! pawopt & optder

!Accumulate d2gxdtfac related to nonlocal operator (Norm-conserving)
!-------------------------------------------------------------------
  if (optder==2.and.paw_opt==0) then
   !Enl is E(Kleinman-Bylander)
!$OMP PARALLEL &
!$OMP PRIVATE(ispinor,ispinor_index,ia,ilmn,iln,enl_,mu)
   do ispinor=1,nspinor
     ispinor_index = ispinor + shift
     do ia=1,nincat
!$OMP DO
       do ilmn=1,nlmn
         iln=indlmn(5,ilmn)
         enl_(1)=enl_ptr(iln,itypat,ispinor_index)
         do mu=1,nd2gxdtfac
           d2gxdtfac_(1:cplex,mu,ilmn,ia,ispinor)=enl_(1)*d2gxdt(1:cplex,mu,ilmn,ia,ispinor)
         end do
       end do
!$OMP END DO
     end do
   end do
!$OMP END PARALLEL
 end if

 DBG_EXIT("COLL")

!Accumulate d2gxdtfac related to nonlocal operator (PAW)
!-------------------------------------------------------------------
  if (optder==2.and.(paw_opt==1.or.paw_opt==2.or.paw_opt==4)) then
   !Enl is psp strength Dij or (Dij-lambda.Sij)

!  === Diagonal term(s) (up-up, down-down)

!  1-Enl is real
   if (cplex_enl==1) then
!$OMP PARALLEL &
!$OMP PRIVATE(ispinor,ispinor_index,ia,index_enl,jlmn,j0lmn,jjlmn,enl_,mu,gxfj,ilmn,ijlmn,gxfi)
     ABI_ALLOCATE(gxfj,(cplex,nd2gxdtfac))
     do ispinor=1,nspinor
       ispinor_index=ispinor+shift
       do ia=1,nincat
         index_enl=atindx1(iatm+ia)
!$OMP DO
         do jlmn=1,nlmn
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1)=enl_ptr(jjlmn,index_enl,ispinor_index)
           if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
           do mu=1,nd2gxdtfac
             gxfj(1:cplex,mu)=d2gxdt(1:cplex,mu,jlmn,ia,ispinor)
             d2gxdtfac_(1:cplex,mu,jlmn,ia,ispinor)=d2gxdtfac_(1:cplex,mu,jlmn,ia,ispinor)+enl_(1)*gxfj(1:cplex,mu)
           end do
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1)=enl_ptr(ijlmn,index_enl,ispinor_index)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
             do mu=1,nd2gxdtfac
               gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,ispinor)
               d2gxdtfac_(1:cplex,mu,jlmn,ia,ispinor)=d2gxdtfac_(1:cplex,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1:cplex)
#if !defined HAVE_OPENMP
               d2gxdtfac_(1:cplex,mu,ilmn,ia,ispinor)=d2gxdtfac_(1:cplex,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(1:cplex,mu)
#endif
             end do
           end do
#if defined HAVE_OPENMP
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=ilmn*(ilmn-1)/2
               ijlmn=i0lmn+jlmn
               enl_(1)=enl_ptr(ijlmn,index_enl,ispinor_index)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               do mu=1,nd2gxdtfac
                 gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,ispinor)
                 d2gxdtfac_(1:cplex,mu,jlmn,ia,ispinor)=d2gxdtfac_(1:cplex,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1:cplex)
               end do
             end do
           end if
#endif
         end do
!$OMP END DO
       end do
     end do
     ABI_DEALLOCATE(gxfj)
!$OMP END PARALLEL

!    2-Enl is complex  ===== D^ss'_ij=D^s's_ji^*
   else
     ABI_CHECK(cplex_fac==cplex_enl,"BUG: invalid cplex_fac/=cplex_enl!")

     if (nspinortot==1) then ! -------------> NO SPINORS

!$OMP PARALLEL &
!$OMP PRIVATE(ia,index_enl,jlmn,j0lmn,jjlmn,enl_,mu,gxfj,ilmn,ijlmn,gxfi)
       ABI_ALLOCATE(gxfj,(cplex,nd2gxdtfac))
       do ia=1,nincat
         index_enl=atindx1(iatm+ia)
!$OMP DO
         do jlmn=1,nlmn
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1)=enl_ptr(2*jjlmn-1,index_enl,1)
           if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
           do mu=1,nd2gxdtfac
             if(cplex_d2gxdt(mu)==2)then
               cplex_ = 2 ; gxfj(1,mu) = zero ; gxfj(2,mu) = d2gxdt(1,mu,jlmn,ia,1)
             else
               cplex_ = cplex ; gxfj(1:cplex,mu)=d2gxdt(1:cplex,mu,jlmn,ia,1)
             end if
             d2gxdtfac_(1,mu,jlmn,ia,1)=d2gxdtfac_(1,mu,jlmn,ia,1)+enl_(1)*gxfj(1,mu)
             if (cplex_==2) d2gxdtfac_(2,mu,jlmn,ia,1)=d2gxdtfac_(2,mu,jlmn,ia,1)+enl_(1)*gxfj(2,mu)
           end do
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,1)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
             do mu=1,nd2gxdtfac
               if(cplex_d2gxdt(mu)==2)then
                 cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = d2gxdt(1,mu,ilmn,ia,1)
               else
                 cplex_ = cplex ; gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,1)
               end if
               d2gxdtfac_(1,mu,jlmn,ia,1)=d2gxdtfac_(1,mu,jlmn,ia,1)+enl_(1)*gxfi(1)
               d2gxdtfac_(2,mu,jlmn,ia,1)=d2gxdtfac_(2,mu,jlmn,ia,1)-enl_(2)*gxfi(1)
#if !defined HAVE_OPENMP
               d2gxdtfac_(1,mu,ilmn,ia,1)=d2gxdtfac_(1,mu,ilmn,ia,1)+enl_(1)*gxfj(1,mu)
               d2gxdtfac_(2,mu,ilmn,ia,1)=d2gxdtfac_(2,mu,ilmn,ia,1)+enl_(2)*gxfj(1,mu)
#endif
               if (cplex_==2) then
                 d2gxdtfac_(1,mu,jlmn,ia,1)=d2gxdtfac_(1,mu,jlmn,ia,1)+enl_(2)*gxfi(2)
                 d2gxdtfac_(2,mu,jlmn,ia,1)=d2gxdtfac_(2,mu,jlmn,ia,1)+enl_(1)*gxfi(2)
#if !defined HAVE_OPENMP
                 d2gxdtfac_(1,mu,ilmn,ia,1)=d2gxdtfac_(1,mu,ilmn,ia,1)-enl_(2)*gxfj(2,mu)
                 d2gxdtfac_(2,mu,ilmn,ia,1)=d2gxdtfac_(2,mu,ilmn,ia,1)+enl_(1)*gxfj(2,mu)
#endif
               end if
             end do
           end do
#if defined HAVE_OPENMP
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=ilmn*(ilmn-1)/2
               ijlmn=i0lmn+jlmn
               enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,1)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               do mu=1,nd2gxdtfac
                 if(cplex_d2gxdt(mu)==2)then
                   cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = d2gxdt(1,mu,ilmn,ia,1)
                 else
                   cplex_ = cplex ; gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,1)
                 end if
                 d2gxdtfac_(1,mu,jlmn,ia,1)=d2gxdtfac_(1,mu,jlmn,ia,1)+enl_(1)*gxfi(1)
                 d2gxdtfac_(2,mu,jlmn,ia,1)=d2gxdtfac_(2,mu,jlmn,ia,1)+enl_(2)*gxfi(1)
                 if (cplex_==2) then
                   d2gxdtfac_(1,mu,jlmn,ia,1)=d2gxdtfac_(1,mu,jlmn,ia,1)-enl_(2)*gxfi(2)
                   d2gxdtfac_(2,mu,jlmn,ia,1)=d2gxdtfac_(2,mu,jlmn,ia,1)+enl_(1)*gxfi(2)
                 end if
               end do
             end do
           end if
#endif
         end do
!$OMP END DO
       end do
       ABI_DEALLOCATE(gxfj)
!$OMP END PARALLEL

     else ! -------------> SPINORIAL CASE

!  === Diagonal term(s) (up-up, down-down)

!$OMP PARALLEL &
!$OMP PRIVATE(ispinor,ispinor_index,ia,index_enl), &
!$OMP PRIVATE(jlmn,j0lmn,jjlmn,enl_,mu,gxfj,ilmn,ijlmn,gxfi)
       ABI_ALLOCATE(gxfj,(cplex,nd2gxdtfac))
       do ispinor=1,nspinor
         ispinor_index = ispinor + shift
         do ia=1,nincat
           index_enl=atindx1(iatm+ia)
!$OMP DO
           do jlmn=1,nlmn
             j0lmn=jlmn*(jlmn-1)/2
             jjlmn=j0lmn+jlmn
             enl_(1)=enl_ptr(2*jjlmn-1,index_enl,ispinor_index)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
             do mu=1,nd2gxdtfac
               if(cplex_d2gxdt(mu)==2)then
                 cplex_ = 2 ; gxfj(1,mu) = zero ; gxfj(2,mu) = d2gxdt(1,mu,jlmn,ia,ispinor)
               else
                 cplex_ = cplex ; gxfj(1:cplex,mu)=d2gxdt(1:cplex,mu,jlmn,ia,ispinor)
               end if
               d2gxdtfac_(1,mu,jlmn,ia,ispinor)=d2gxdtfac_(1,mu,jlmn,ia,ispinor)+enl_(1)*gxfj(1,mu)
               if (cplex_==2) d2gxdtfac_(2,mu,jlmn,ia,ispinor)=d2gxdtfac_(2,mu,jlmn,ia,ispinor)+enl_(1)*gxfj(2,mu)
             end do
             do ilmn=1,jlmn-1
               ijlmn=j0lmn+ilmn
               enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,ispinor_index)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               do mu=1,nd2gxdtfac
                 if(cplex_d2gxdt(mu)==2)then
                   cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = d2gxdt(1,mu,ilmn,ia,ispinor)
                 else
                   cplex_ = cplex ; gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,ispinor)
                 end if
                 d2gxdtfac_(1,mu,jlmn,ia,ispinor)=d2gxdtfac_(1,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1)
                 d2gxdtfac_(2,mu,jlmn,ia,ispinor)=d2gxdtfac_(2,mu,jlmn,ia,ispinor)-enl_(2)*gxfi(1)
#if !defined HAVE_OPENMP
                 d2gxdtfac_(1,mu,ilmn,ia,ispinor)=d2gxdtfac_(1,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(1,mu)
                 d2gxdtfac_(2,mu,ilmn,ia,ispinor)=d2gxdtfac_(2,mu,ilmn,ia,ispinor)+enl_(2)*gxfj(1,mu)
#endif
                 if (cplex_==2) then
                   d2gxdtfac_(1,mu,jlmn,ia,ispinor)=d2gxdtfac_(1,mu,jlmn,ia,ispinor)+enl_(2)*gxfi(2)
                   d2gxdtfac_(2,mu,jlmn,ia,ispinor)=d2gxdtfac_(2,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(2)
#if !defined HAVE_OPENMP
                   d2gxdtfac_(1,mu,ilmn,ia,ispinor)=d2gxdtfac_(1,mu,ilmn,ia,ispinor)-enl_(2)*gxfj(2,mu)
                   d2gxdtfac_(2,mu,ilmn,ia,ispinor)=d2gxdtfac_(2,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(2,mu)
#endif
                 end if
               end do
             end do
#if defined HAVE_OPENMP
             if(jlmn<nlmn) then
               do ilmn=jlmn+1,nlmn
                 i0lmn=ilmn*(ilmn-1)/2
                 ijlmn=i0lmn+jlmn
                 enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,ispinor_index)
                 if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
                 do mu=1,nd2gxdtfac
                   if(cplex_d2gxdt(mu)==2)then
                     cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = d2gxdt(1,mu,ilmn,ia,ispinor)
                   else
                     cplex_ = cplex ; gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,ispinor)
                   end if
                   d2gxdtfac_(1,mu,jlmn,ia,ispinor)=d2gxdtfac_(1,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1)
                   d2gxdtfac_(2,mu,jlmn,ia,ispinor)=d2gxdtfac_(2,mu,jlmn,ia,ispinor)+enl_(2)*gxfi(1)
                   if (cplex_==2) then
                     d2gxdtfac_(1,mu,jlmn,ia,ispinor)=d2gxdtfac_(1,mu,jlmn,ia,ispinor)-enl_(2)*gxfi(2)
                     d2gxdtfac_(2,mu,jlmn,ia,ispinor)=d2gxdtfac_(2,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(2)
                   end if
                 end do
               end do
             end if
#endif
           end do
!$OMP END DO
         end do
       end do
       ABI_DEALLOCATE(gxfj)
!$OMP END PARALLEL
     end if !nspinortot
   end if !complex

!  === Off-diagonal term(s) (up-down, down-up)

!  --- No parallelization over spinors ---
   if (nspinortot==2.and.nspinor==nspinortot) then
     ABI_CHECK(cplex_enl==2,"BUG: invalid cplex_enl/=2!")
     ABI_CHECK(cplex_fac==2,"BUG: invalid cplex_fac/=2!")
!$OMP PARALLEL &
!$OMP PRIVATE(ispinor,jspinor,ia,index_enl), &
!$OMP PRIVATE(jlmn,j0lmn,jjlmn,enl_,mu,gxfi,gxfj,ilmn,ijlmn)
     ABI_ALLOCATE(gxfj,(cplex,nd2gxdtfac))
     do ispinor=1,nspinor
       jspinor=3-ispinor
       do ia=1,nincat
         index_enl=atindx1(iatm+ia)
!$OMP DO
         do jlmn=1,nlmn
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1:2)=enl_ptr(2*jjlmn-1:2*jjlmn,index_enl,2+ispinor)
           do mu=1,nd2gxdtfac
             if(cplex_d2gxdt(mu)==2)then
               cplex_ = 2 ;
               gxfi(1)    = zero ; gxfi(2)    = d2gxdt(1,mu,jlmn,ia,ispinor)
               gxfj(1,mu) = zero ; gxfj(2,mu) = d2gxdt(1,mu,jlmn,ia,jspinor)
             else
               cplex_ = cplex ;
               gxfi(1:cplex)   =d2gxdt(1:cplex,mu,jlmn,ia,ispinor)
               gxfj(1:cplex,mu)=d2gxdt(1:cplex,mu,jlmn,ia,jspinor)
             end if
             d2gxdtfac_(1,mu,jlmn,ia,jspinor)=d2gxdtfac_(1,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(1)
             d2gxdtfac_(2,mu,jlmn,ia,jspinor)=d2gxdtfac_(2,mu,jlmn,ia,jspinor)-enl_(2)*gxfi(1)
             if (cplex_==2) then
               d2gxdtfac_(1,mu,jlmn,ia,jspinor)=d2gxdtfac_(1,mu,jlmn,ia,jspinor)+enl_(2)*gxfi(2)
               d2gxdtfac_(2,mu,jlmn,ia,jspinor)=d2gxdtfac_(2,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(2)
             end if
           end do
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor)
             do mu=1,nd2gxdtfac
               if(cplex_d2gxdt(mu)==2)then
                 cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = d2gxdt(1,mu,ilmn,ia,ispinor)
               else
                 cplex_ = cplex ; gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,ispinor)
               end if
               d2gxdtfac_(1,mu,jlmn,ia,jspinor)=d2gxdtfac_(1,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(1)
               d2gxdtfac_(2,mu,jlmn,ia,jspinor)=d2gxdtfac_(2,mu,jlmn,ia,jspinor)-enl_(2)*gxfi(1)
#if !defined HAVE_OPENMP
               d2gxdtfac_(1,mu,ilmn,ia,ispinor)=d2gxdtfac_(1,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(1,mu)
               d2gxdtfac_(2,mu,ilmn,ia,ispinor)=d2gxdtfac_(2,mu,ilmn,ia,ispinor)+enl_(2)*gxfj(1,mu)
#endif
               if (cplex_==2) then
                 d2gxdtfac_(1,mu,jlmn,ia,jspinor)=d2gxdtfac_(1,mu,jlmn,ia,jspinor)+enl_(2)*gxfi(2)
                 d2gxdtfac_(2,mu,jlmn,ia,jspinor)=d2gxdtfac_(2,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(2)
#if !defined HAVE_OPENMP
                 d2gxdtfac_(1,mu,ilmn,ia,ispinor)=d2gxdtfac_(1,mu,ilmn,ia,ispinor)-enl_(2)*gxfj(2,mu)
                 d2gxdtfac_(2,mu,ilmn,ia,ispinor)=d2gxdtfac_(2,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(2,mu)
#endif
               end if
             end do !mu
           end do !ilmn
#if defined HAVE_OPENMP
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=ilmn*(ilmn-1)/2
               ijlmn=i0lmn+jlmn
               enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor)
               do mu=1,nd2gxdtfac
                 if(cplex_d2gxdt(mu)==2)then
                   cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = d2gxdt(1,mu,ilmn,ia,jspinor)
                 else
                   cplex_ = cplex ; gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,jspinor)
                 end if
                 d2gxdtfac_(1,mu,jlmn,ia,ispinor)=d2gxdtfac_(1,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1)
                 d2gxdtfac_(2,mu,jlmn,ia,ispinor)=d2gxdtfac_(2,mu,jlmn,ia,ispinor)+enl_(2)*gxfi(1)
                 if (cplex_==2) then
                   d2gxdtfac_(1,mu,jlmn,ia,ispinor)=d2gxdtfac_(1,mu,jlmn,ia,ispinor)-enl_(2)*gxfi(2)
                   d2gxdtfac_(2,mu,jlmn,ia,ispinor)=d2gxdtfac_(2,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(2)
                 end if
               end do !mu
             end do !ilmn
           end if
#endif
         end do !jmln
!$OMP END DO
       end do !ia
     end do !ispinor
     ABI_DEALLOCATE(gxfj)
!$OMP END PARALLEL

!    --- Parallelization over spinors ---
   else if (nspinortot==2.and.nspinor/=nspinortot) then
     ABI_CHECK(cplex_enl==2,"BUG: opernlc_ylm: invalid cplex_enl/=2!")
     ABI_CHECK(cplex_fac==2,"BUG: opernlc_ylm: invalid cplex_fac/=2!")
     ABI_ALLOCATE(d2gxdtfac_offdiag,(cplex_fac,nd2gxdtfac,nlmn,nincat,nspinortot))
!$OMP PARALLEL &
!$OMP PRIVATE(ia,index_enl), &
!$OMP PRIVATE(jlmn,j0lmn,ilmn,i0lmn,ijlmn,enl_,jilmn,mu,gxfi)
!$OMP WORKSHARE
     d2gxdtfac_offdiag(:,:,:,:,:)=zero
!$OMP END WORKSHARE
!$OMP SINGLE
     ispinor_index=mpi_enreg%me_spinor+1
     jspinor_index=3-ispinor_index
     if (ispinor_index==1) then
       ijspin=3;jispin=4
     else
       ijspin=4;jispin=3
     end if
!$OMP END SINGLE
     do ia=1,nincat
       index_enl=atindx1(iatm+ia)
!$OMP DO
       do jlmn=1,nlmn
         j0lmn=jlmn*(jlmn-1)/2
         do ilmn=1,nlmn
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
           do mu=1,nd2gxdtfac
             if(cplex_d2gxdt(mu)==2)then
               cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = d2gxdt(1,mu,ilmn,ia,1)
             else
               cplex_ = cplex ; gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,1)
             end if
             d2gxdtfac_offdiag(1,mu,jlmn,ia,jspinor_index)= &
&                  d2gxdtfac_offdiag(1,mu,jlmn,ia,jspinor_index)+enl_(1)*gxfi(1)
             d2gxdtfac_offdiag(2,mu,jlmn,ia,jspinor_index)= &
&                  d2gxdtfac_offdiag(2,mu,jlmn,ia,jspinor_index)+enl_(2)*gxfi(1)
             if (cplex_==2) then
               d2gxdtfac_offdiag(1,mu,jlmn,ia,jspinor_index)= &
&                    d2gxdtfac_offdiag(1,mu,jlmn,ia,jspinor_index)-enl_(2)*gxfi(2)
               d2gxdtfac_offdiag(2,mu,jlmn,ia,jspinor_index)= &
&                    d2gxdtfac_offdiag(2,mu,jlmn,ia,jspinor_index)+enl_(1)*gxfi(2)
             end if
           end do
         end do !ilmn
       end do !jlmn
!$OMP END DO
     end do !iat
!$OMP SINGLE
     call xmpi_sum(d2gxdtfac_offdiag,mpi_enreg%comm_spinor,ierr)
!$OMP END SINGLE
!$OMP WORKSHARE
     d2gxdtfac_(:,:,:,:,1)=d2gxdtfac_(:,:,:,:,1)+d2gxdtfac_offdiag(:,:,:,:,ispinor_index)
!$OMP END WORKSHARE
!$OMP END PARALLEL
     ABI_DEALLOCATE(d2gxdtfac_offdiag)
   end if !nspinortot

  end if ! pawopt & optder

!End of loop when a exp(-iqR) phase is present
!------------------------------------------- ------------------------

!When iphase=1, gxfac and gxfac_ point to the same memory space
!When iphase=2, we add i.gxfac_ to gxfac
  if (iphase==2) then
!$OMP PARALLEL PRIVATE(ispinor,ia,ilmn,mu)
!$OMP DO COLLAPSE(3)
    do ispinor=1,nspinor
      do ia=1,nincat
        do ilmn=1,nlmn
          gxfac(1,ilmn,ia,ispinor)=gxfac(1,ilmn,ia,ispinor)-gxfac_(2,ilmn,ia,ispinor)
          gxfac(2,ilmn,ia,ispinor)=gxfac(2,ilmn,ia,ispinor)+gxfac_(1,ilmn,ia,ispinor)
        end do
      end do
    end do
    ABI_DEALLOCATE(gxfac_)
    if (optder>=1) then
!$OMP DO COLLAPSE(4)
      do ispinor=1,nspinor
        do ia=1,nincat
          do ilmn=1,nlmn
            do mu=1,ndgxdtfac
              dgxdtfac(1,mu,ilmn,ia,ispinor)=dgxdtfac(1,mu,ilmn,ia,ispinor)-dgxdtfac_(2,mu,ilmn,ia,ispinor)
              dgxdtfac(2,mu,ilmn,ia,ispinor)=dgxdtfac(2,mu,ilmn,ia,ispinor)+dgxdtfac_(1,mu,ilmn,ia,ispinor)
            end do
          end do
        end do
      end do
      ABI_DEALLOCATE(dgxdtfac_)
    end if
    if (optder>=2) then
!$OMP DO COLLAPSE(4)
      do ispinor=1,nspinor
        do ia=1,nincat
          do ilmn=1,nlmn
            do mu=1,ndgxdtfac
              d2gxdtfac(1,mu,ilmn,ia,ispinor)=d2gxdtfac(1,mu,ilmn,ia,ispinor)-d2gxdtfac_(2,mu,ilmn,ia,ispinor)
              d2gxdtfac(2,mu,ilmn,ia,ispinor)=d2gxdtfac(2,mu,ilmn,ia,ispinor)+d2gxdtfac_(1,mu,ilmn,ia,ispinor)
            end do
          end do
        end do
      end do
      ABI_DEALLOCATE(d2gxdtfac_)
    end if
!$OMP END PARALLEL
  end if

!End loop over real/imaginary part of the exp(-iqR) phase
 end do


!Accumulate gxfac related to overlap (Sij) (PAW)
!------------------------------------------- ------------------------
 if (paw_opt==3.or.paw_opt==4) then ! Use Sij, overlap contribution
!$OMP PARALLEL &
!$OMP PRIVATE(ispinor,ia,jlmn,j0lmn,jjlmn,jlm,sijr,ilmn,ilm,ijlmn,gxi,gxj)
!$OMP WORKSHARE
   gxfac_sij(1:cplex,1:nlmn,1:nincat,1:nspinor)=zero
!$OMP END WORKSHARE
!$OMP DO COLLAPSE(3)
   do ispinor=1,nspinor
     do ia=1,nincat
       do jlmn=1,nlmn
         j0lmn=jlmn*(jlmn-1)/2
         jjlmn=j0lmn+jlmn
         jlm=indlmn(4,jlmn)
         sijr=sij(jjlmn)
         gxj(1:cplex)=gx(1:cplex,jlmn,ia,ispinor)
         gxfac_sij(1:cplex,jlmn,ia,ispinor)=gxfac_sij(1:cplex,jlmn,ia,ispinor)+sijr*gxj(1:cplex)
         do ilmn=1,jlmn-1
           ilm=indlmn(4,ilmn)
          !if (ilm==jlm) then
           ijlmn=j0lmn+ilmn
           sijr=sij(ijlmn)
           gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
           gxfac_sij(1:cplex,jlmn,ia,ispinor)=gxfac_sij(1:cplex,jlmn,ia,ispinor)+sijr*gxi(1:cplex)
#if !defined HAVE_OPENMP
           gxfac_sij(1:cplex,ilmn,ia,ispinor)=gxfac_sij(1:cplex,ilmn,ia,ispinor)+sijr*gxj(1:cplex)
#endif
          !end if
         end do
#if defined HAVE_OPENMP
         if(jlmn<nlmn) then
           do ilmn=jlmn+1,nlmn
             ilm=indlmn(4,ilmn)
             !if (ilm==jlm) then
             i0lmn=ilmn*(ilmn-1)/2
             ijlmn=i0lmn+jlmn
             sijr=sij(ijlmn)
             gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
             gxfac_sij(1:cplex,jlmn,ia,ispinor)=gxfac_sij(1:cplex,jlmn,ia,ispinor)+sijr*gxi(1:cplex)
             !end if
           end do
         end if
#endif
       end do
     end do
   end do
!$OMP END DO
!$OMP END PARALLEL
 end if

!Accumulate dgxdtfac related to overlap (Sij) (PAW)
!-------------------------------------------------------------------
 if (optder>=1.and.(paw_opt==3.or.paw_opt==4)) then ! Use Sij, overlap contribution
!$OMP PARALLEL &
!$OMP PRIVATE(ispinor,ia), &
!$OMP PRIVATE(jlmn,j0lmn,jjlmn,sijr,mu,gxfj,ilmn,ijlmn,gxfi)
   ABI_ALLOCATE(gxfj,(cplex,ndgxdtfac))
!$OMP WORKSHARE
   dgxdtfac_sij(1:cplex,1:ndgxdtfac,1:nlmn,1:nincat,1:nspinor)=zero
!$OMP END WORKSHARE
   do ispinor=1,nspinor
     do ia=1,nincat
!$OMP DO
       do jlmn=1,nlmn
         j0lmn=jlmn*(jlmn-1)/2
         jjlmn=j0lmn+jlmn
         sijr=sij(jjlmn)
         do mu=1,ndgxdtfac
           gxfj(1:cplex,mu)=dgxdt(1:cplex,mu,jlmn,ia,ispinor)
           dgxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)=dgxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)+sijr*gxfj(1:cplex,mu)
         end do
         do ilmn=1,jlmn-1
           ijlmn=j0lmn+ilmn
           sijr=sij(ijlmn)
           do mu=1,ndgxdtfac
             gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
             dgxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)=dgxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)+sijr*gxfi(1:cplex)
#if !defined HAVE_OPENMP
             dgxdtfac_sij(1:cplex,mu,ilmn,ia,ispinor)=dgxdtfac_sij(1:cplex,mu,ilmn,ia,ispinor)+sijr*gxfj(1:cplex,mu)
#endif
           end do
         end do
#if defined HAVE_OPENMP
         if(jlmn<nlmn) then
           do ilmn=jlmn+1,nlmn
             i0lmn=ilmn*(ilmn-1)/2
             ijlmn=i0lmn+jlmn
             sijr=sij(ijlmn)
             do mu=1,ndgxdtfac
               gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
               dgxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)=dgxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)+sijr*gxfi(1:cplex)
             end do
           end do
         end if
#endif
       end do
!$OMP END DO
     end do
   end do
   ABI_DEALLOCATE(gxfj)
!$OMP END PARALLEL
 end if

!Accumulate d2gxdtfac related to overlap (Sij) (PAW)
!-------------------------------------------------------------------
 if (optder==2.and.(paw_opt==3.or.paw_opt==4)) then ! Use Sij, overlap contribution
!$OMP PARALLEL &
!$OMP PRIVATE(ispinor,ia), &
!$OMP PRIVATE(jlmn,j0lmn,jjlmn,sijr,mu,gxfj,ilmn,ijlmn,gxfi)
   ABI_ALLOCATE(gxfj,(cplex,nd2gxdtfac))
!$OMP WORKSHARE
   d2gxdtfac_sij(1:cplex,1:nd2gxdtfac,1:nlmn,1:nincat,1:nspinor)=zero
!$OMP END WORKSHARE
   do ispinor=1,nspinor
     do ia=1,nincat
!$OMP DO
       do jlmn=1,nlmn
         j0lmn=jlmn*(jlmn-1)/2
         jjlmn=j0lmn+jlmn
         sijr=sij(jjlmn)
         do mu=1,nd2gxdtfac
           gxfj(1:cplex,mu)=d2gxdt(1:cplex,mu,jlmn,ia,ispinor)
           d2gxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)= &
&                d2gxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)+sijr*gxfj(1:cplex,mu)
         end do
         do ilmn=1,jlmn-1
           ijlmn=j0lmn+ilmn
           sijr=sij(ijlmn)
           do mu=1,nd2gxdtfac
             gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,ispinor)
             d2gxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)= &
&                  d2gxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)+sijr*gxfi(1:cplex)
#if !defined HAVE_OPENMP
             d2gxdtfac_sij(1:cplex,mu,ilmn,ia,ispinor)= &
&                  d2gxdtfac_sij(1:cplex,mu,ilmn,ia,ispinor)+sijr*gxfj(1:cplex,mu)
#endif
           end do
         end do
#if defined HAVE_OPENMP
         if(jlmn<nlmn) then
           do ilmn=jlmn+1,nlmn
             i0lmn=ilmn*(ilmn-1)/2
             ijlmn=i0lmn+jlmn
             sijr=sij(ijlmn)
             do mu=1,nd2gxdtfac
               gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,ispinor)
               d2gxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)= &
&                    d2gxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)+sijr*gxfi(1:cplex)
             end do
           end do
         end if
#endif
       end do
!$OMP END DO
     end do
   end do
   ABI_DEALLOCATE(gxfj)
!$OMP END PARALLEL
 end if

end subroutine opernlc_ylm
!!***

end module m_opernlc_ylm
!!***
