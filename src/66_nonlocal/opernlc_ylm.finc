!{\src2tex{textfont=tt}}
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
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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
!!  enl(cplex_enl*dimenl1,dimenl2,nspinortot**2)=
!!  ->Norm conserving : ==== when paw_opt=0 ====
!!                      (Real) Kleinman-Bylander energies (hartree)
!!                      dimenl1=lmnmax  -  dimenl2=ntypat
!!  ->PAW :             ==== when paw_opt=1, 2 or 4 ====
!!                      (Real or complex, hermitian) Dij coefs to connect projectors
!!                      dimenl1=cplex_enl*lmnmax*(lmnmax+1)/2  -  dimenl2=natom
!!                      These are complex numbers if cplex_enl=2
!!                        enl(:,:,1) contains Dij^up-up
!!                        enl(:,:,2) contains Dij^dn-dn
!!                        enl(:,:,3) contains Dij^up-dn (only if nspinor=2)
!!                        enl(:,:,4) contains Dij^dn-up (only if nspinor=2)
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
!!  paw_opt= define the nonlocal operator concerned with:
!!           paw_opt=0 : Norm-conserving Vnl (use of Kleinman-Bylander ener.)
!!           paw_opt=1 : PAW nonlocal part of H (use of Dij coeffs)
!!           paw_opt=2 : PAW: (Vnl-lambda.Sij) (Sij=overlap matrix)
!!           paw_opt=3 : PAW overlap matrix (Sij)
!!           paw_opt=4 : both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
!!  sij(nlm*(nlmn+1)/2)=overlap matrix components (only if paw_opt=2, 3 or 4)
!!  hermdij=optional logical argument governing whether Dij are to be applied explicitly Hermitian (true)
!!          or just symmetry Dij=Dji (false,default)
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
!! PARENTS
!!      m_gemm_nonlop,nonlop_ylm
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine opernlc_ylm(atindx1,cplex,cplex_dgxdt,cplex_d2gxdt,cplex_enl,cplex_fac,dgxdt,dgxdtfac,dgxdtfac_sij,&
&                      d2gxdt,d2gxdtfac,d2gxdtfac_sij,dimenl1,dimenl2,enl,gx,gxfac,gxfac_sij,iatm,indlmn,itypat,&
&                      lambda,mpi_enreg,natom,ndgxdt,ndgxdtfac,nd2gxdt,nd2gxdtfac,nincat,nlmn,&
&                      nspinor,nspinortot,optder,paw_opt,sij,hermdij)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'opernlc_ylm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,cplex_enl,cplex_fac,dimenl1,dimenl2,iatm,itypat
 integer,intent(in) :: natom,ndgxdt,ndgxdtfac,nd2gxdt,nd2gxdtfac,nincat,nspinor,nspinortot,optder,paw_opt
 integer,intent(inout) :: nlmn
 real(dp) :: lambda
 logical,optional,intent(in) :: hermdij
 type(MPI_type) , intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom),indlmn(6,nlmn),cplex_dgxdt(ndgxdt),cplex_d2gxdt(nd2gxdt)
 real(dp),intent(in) :: dgxdt(cplex,ndgxdt,nlmn,nincat,nspinor)
 real(dp),intent(in) :: d2gxdt(cplex,nd2gxdt,nlmn,nincat,nspinor)
 real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinortot**2)
 real(dp),intent(inout) :: gx(cplex,nlmn,nincat,nspinor)
 real(dp),intent(in) :: sij(((paw_opt+1)/3)*nlmn*(nlmn+1)/2)
 real(dp),intent(out) :: dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)
 real(dp),intent(out) :: dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat,nspinor*(paw_opt/3))
 real(dp),intent(out) :: d2gxdtfac(cplex_fac,nd2gxdtfac,nlmn,nincat,nspinor)
 real(dp),intent(out) :: d2gxdtfac_sij(cplex,nd2gxdtfac,nlmn,nincat,nspinor*(paw_opt/3))
 real(dp),intent(out) :: gxfac(cplex_fac,nlmn,nincat,nspinor)
 real(dp),intent(out) :: gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))

!Local variables-------------------------------
!Arrays
!scalars
 integer :: cplex_,ia,ierr,ijlmn,ijspin,ilm,ilmn,i0lmn,iln,index_enl,ispinor,ispinor_index
 integer :: j0lmn,jilmn,jispin,jjlmn,jlm,jlmn,jspinor,jspinor_index,mu,shift
 real(dp) :: sijr
 logical :: hermdij_
!arrays
 real(dp) :: enl_(2),gxfi(2),gxi(cplex),gxj(cplex)
 real(dp),allocatable :: d2gxdtfac_(:,:,:,:,:),dgxdtfac_(:,:,:,:,:),gxfac_(:,:,:,:),gxfj(:,:) 

! *************************************************************************

 DBG_ENTER("COLL")

!Parallelization over spinors treatment
 shift=0;if (mpi_enreg%paral_spinor==1) shift=mpi_enreg%me_spinor

! hermdij_ .false. invokes original coding, with Dij=Dji in complex case
! hermdij_ .true. invokes Dij = Dji*, needed for magnetic field cases
 hermdij_=.false.;if(present(hermdij)) hermdij_=hermdij

!Accumulate gxfac related to non-local operator (Norm-conserving)
!-------------------------------------------------------------------
 if (paw_opt==0) then                  ! Enl is E(Kleinman-Bylander)
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
         enl_(1)=enl(iln,itypat,ispinor_index)
         gxfac(1:cplex,ilmn,ia,ispinor)=enl_(1)*gx(1:cplex,ilmn,ia,ispinor)
       end do
     end do
   end do
!$OMP END DO
!$OMP END PARALLEL 
 end if

!Accumulate gxfac related to nonlocal operator (PAW)
!-------------------------------------------------------------------
 if (paw_opt==1.or.paw_opt==2.or.paw_opt==4) then        ! Enl is psp strength Dij
   gxfac(1:cplex_fac,1:nlmn,1:nincat,1:nspinor)=zero      ! or (Dij-lambda.Sij)
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
           enl_(1)=enl(jjlmn,index_enl,ispinor_index)
           if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
           gxj(1:cplex)=gx(1:cplex,jlmn,ia,ispinor)
           gxfac(1:cplex,jlmn,ia,ispinor)=gxfac(1:cplex,jlmn,ia,ispinor)+enl_(1)*gxj(1:cplex)
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1)=enl(ijlmn,index_enl,ispinor_index)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
             gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
             gxfac(1:cplex,jlmn,ia,ispinor)=gxfac(1:cplex,jlmn,ia,ispinor)+enl_(1)*gxi(1:cplex)
#if !defined HAVE_OPENMP
             gxfac(1:cplex,ilmn,ia,ispinor)=gxfac(1:cplex,ilmn,ia,ispinor)+enl_(1)*gxj(1:cplex)
#endif
           end do
#if defined HAVE_OPENMP
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=(ilmn*(ilmn-1)/2)
               ijlmn=i0lmn+jlmn
               enl_(1)=enl(ijlmn,index_enl,ispinor_index)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
               gxfac(1:cplex,jlmn,ia,ispinor)=gxfac(1:cplex,jlmn,ia,ispinor)+enl_(1)*gxi(1:cplex)
             end do
           end if
#endif
         end do
       end do
     end do
!$OMP END DO 
!$OMP END PARALLEL
!    2-Enl is complex
   else
     ABI_CHECK(cplex_fac==cplex_enl,"BUG: invalid cplex_fac/=cplex_enl!")
     if (nspinortot==1) then    !===== when nspinor=1, D_ij=D_ji except in hermdij case
!$OMP PARALLEL &
!$OMP PRIVATE(ia,index_enl,jlmn,j0lmn,jjlmn,enl_,gxj,ilmn,ijlmn,gxi)
       do ia=1,nincat
         index_enl=atindx1(iatm+ia)
!$OMP DO
         do jlmn=1,nlmn
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1:2)=enl(2*jjlmn-1:2*jjlmn,index_enl,1)
           if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
           gxj(1:cplex)=gx(1:cplex,jlmn,ia,1)
           if(hermdij_) then
             gxfac(1:cplex,jlmn,ia,1)=gxfac(1:cplex,jlmn,ia,1)+enl_(1)*gxj(1:cplex)
           else
             gxfac(1,jlmn,ia,1)=gxfac(1,jlmn,ia,1)+enl_(1)*gxj(1)
             gxfac(2,jlmn,ia,1)=gxfac(2,jlmn,ia,1)+enl_(2)*gxj(1)
             if (cplex==2) then
               gxfac(1,jlmn,ia,1)=gxfac(1,jlmn,ia,1)-enl_(2)*gxj(2)
               gxfac(2,jlmn,ia,1)=gxfac(2,jlmn,ia,1)+enl_(1)*gxj(2)
             end if
           end if
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,1)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
             gxi(1:cplex)=gx(1:cplex,ilmn,ia,1)
             gxfac(1,jlmn,ia,1)=gxfac(1,jlmn,ia,1)+enl_(1)*gxi(1)
             if(hermdij_) then
               gxfac(2,jlmn,ia,1)=gxfac(2,jlmn,ia,1)-enl_(2)*gxi(1)
             else
               gxfac(2,jlmn,ia,1)=gxfac(2,jlmn,ia,1)+enl_(2)*gxi(1)
             end if
#if !defined HAVE_OPENMP
             gxfac(1,ilmn,ia,1)=gxfac(1,ilmn,ia,1)+enl_(1)*gxj(1)
             gxfac(2,ilmn,ia,1)=gxfac(2,ilmn,ia,1)+enl_(2)*gxj(1)
#endif
             if (cplex==2) then
               if(hermdij_) then
                 gxfac(1,jlmn,ia,1)=gxfac(1,jlmn,ia,1)+enl_(2)*gxi(2)
               else
                 gxfac(1,jlmn,ia,1)=gxfac(1,jlmn,ia,1)-enl_(2)*gxi(2)
               end if
               gxfac(2,jlmn,ia,1)=gxfac(2,jlmn,ia,1)+enl_(1)*gxi(2)
#if !defined HAVE_OPENMP
               gxfac(1,ilmn,ia,1)=gxfac(1,ilmn,ia,1)-enl_(2)*gxj(2)
               gxfac(2,ilmn,ia,1)=gxfac(2,ilmn,ia,1)+enl_(1)*gxj(2)
#endif
             end if
           end do
#if defined HAVE_OPENMP
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=ilmn*(ilmn-1)/2
               ijlmn=i0lmn+jlmn
               enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,1)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               gxi(1:cplex)=gx(1:cplex,ilmn,ia,1)
               gxfac(1,jlmn,ia,1)=gxfac(1,jlmn,ia,1)+enl_(1)*gxi(1)
               gxfac(2,jlmn,ia,1)=gxfac(2,jlmn,ia,1)+enl_(2)*gxi(1)
               if (cplex==2) then
                 gxfac(1,jlmn,ia,1)=gxfac(1,jlmn,ia,1)-enl_(2)*gxi(2)
                 gxfac(2,jlmn,ia,1)=gxfac(2,jlmn,ia,1)+enl_(1)*gxi(2)
               end if
             end do
           end if
#endif
         end do
!$OMP END DO
       end do
!$OMP END PARALLEL
     else                    !===== when nspinor=2, D_ij=D_ji^*
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
             enl_(1)=enl(2*jjlmn-1,index_enl,ispinor_index)  ! enl_ii is real for up-up and dn-dn
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
             gxj(1:cplex)=gx(1:cplex,jlmn,ia,ispinor)
             gxfac(1:cplex,jlmn,ia,ispinor)=gxfac(1:cplex,jlmn,ia,ispinor)+enl_(1)*gxj(1:cplex)
             do ilmn=1,jlmn-1
               ijlmn=j0lmn+ilmn
               enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,ispinor_index)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
               gxfac(1,jlmn,ia,ispinor)=gxfac(1,jlmn,ia,ispinor)+enl_(1)*gxi(1)
               gxfac(2,jlmn,ia,ispinor)=gxfac(2,jlmn,ia,ispinor)-enl_(2)*gxi(1)
#if !defined HAVE_OPENMP
               gxfac(1,ilmn,ia,ispinor)=gxfac(1,ilmn,ia,ispinor)+enl_(1)*gxj(1)
               gxfac(2,ilmn,ia,ispinor)=gxfac(2,ilmn,ia,ispinor)+enl_(2)*gxj(1)
#endif
               if (cplex==2) then
                 gxfac(1,jlmn,ia,ispinor)=gxfac(1,jlmn,ia,ispinor)+enl_(2)*gxi(2)
                 gxfac(2,jlmn,ia,ispinor)=gxfac(2,jlmn,ia,ispinor)+enl_(1)*gxi(2)
#if !defined HAVE_OPENMP
                 gxfac(1,ilmn,ia,ispinor)=gxfac(1,ilmn,ia,ispinor)-enl_(2)*gxj(2)
                 gxfac(2,ilmn,ia,ispinor)=gxfac(2,ilmn,ia,ispinor)+enl_(1)*gxj(2)
#endif
               end if
             end do
#if defined HAVE_OPENMP
             if(jlmn<nlmn) then
               do ilmn=jlmn+1,nlmn
                 i0lmn=ilmn*(ilmn-1)/2
                 ijlmn=i0lmn+jlmn
                 enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,ispinor_index)
                 if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
                 gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
                 gxfac(1,jlmn,ia,ispinor)=gxfac(1,jlmn,ia,ispinor)+enl_(1)*gxi(1)
                 gxfac(2,jlmn,ia,ispinor)=gxfac(2,jlmn,ia,ispinor)+enl_(2)*gxi(1)
                 if (cplex==2) then
                   gxfac(1,jlmn,ia,ispinor)=gxfac(1,jlmn,ia,ispinor)-enl_(2)*gxi(2)
                   gxfac(2,jlmn,ia,ispinor)=gxfac(2,jlmn,ia,ispinor)+enl_(1)*gxi(2)
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
           enl_(1:2)=enl(2*jjlmn-1:2*jjlmn,index_enl,2+ispinor )
           gxi(1:cplex)=gx(1:cplex,jlmn,ia,ispinor)
           gxfac(1,jlmn,ia,jspinor)=gxfac(1,jlmn,ia,jspinor)+enl_(1)*gxi(1)
           gxfac(2,jlmn,ia,jspinor)=gxfac(2,jlmn,ia,jspinor)-enl_(2)*gxi(1)
           if (cplex==2) then
             gxfac(1,jlmn,ia,jspinor)=gxfac(1,jlmn,ia,jspinor)+enl_(2)*gxi(2)
             gxfac(2,jlmn,ia,jspinor)=gxfac(2,jlmn,ia,jspinor)+enl_(1)*gxi(2)
           end if
#if !defined HAVE_OPENMP
           gxj(1:cplex)=gx(1:cplex,jlmn,ia,jspinor)
#endif
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor)
             gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
             gxfac(1,jlmn,ia,jspinor)=gxfac(1,jlmn,ia,jspinor)+enl_(1)*gxi(1)
             gxfac(2,jlmn,ia,jspinor)=gxfac(2,jlmn,ia,jspinor)-enl_(2)*gxi(1)
#if !defined HAVE_OPENMP
             gxfac(1,ilmn,ia,ispinor)=gxfac(1,ilmn,ia,ispinor)+enl_(1)*gxj(1)
             gxfac(2,ilmn,ia,ispinor)=gxfac(2,ilmn,ia,ispinor)+enl_(2)*gxj(1)
#endif
             if (cplex==2) then
               gxfac(1,jlmn,ia,jspinor)=gxfac(1,jlmn,ia,jspinor)+enl_(2)*gxi(2)
               gxfac(2,jlmn,ia,jspinor)=gxfac(2,jlmn,ia,jspinor)+enl_(1)*gxi(2)
#if !defined HAVE_OPENMP
               gxfac(1,ilmn,ia,ispinor)=gxfac(1,ilmn,ia,ispinor)-enl_(2)*gxj(2)
               gxfac(2,ilmn,ia,ispinor)=gxfac(2,ilmn,ia,ispinor)+enl_(1)*gxj(2)
#endif
             end if
           end do
#if defined HAVE_OPENMP
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=ilmn*(ilmn-1)/2
               ijlmn=i0lmn+jlmn
               enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor)
               gxi(1:cplex)=gx(1:cplex,ilmn,ia,jspinor)
               gxfac(1,jlmn,ia,ispinor)=gxfac(1,jlmn,ia,ispinor)+enl_(1)*gxi(1)
               gxfac(2,jlmn,ia,ispinor)=gxfac(2,jlmn,ia,ispinor)+enl_(2)*gxi(1)
               if (cplex==2) then
                 gxfac(1,jlmn,ia,ispinor)=gxfac(1,jlmn,ia,ispinor)-enl_(2)*gxi(2)
                 gxfac(2,jlmn,ia,ispinor)=gxfac(2,jlmn,ia,ispinor)+enl_(1)*gxi(2)
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
     ABI_ALLOCATE(gxfac_,(cplex_fac,nlmn,nincat,nspinortot))
     gxfac_(:,:,:,:)=zero
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
             enl_(1)= enl(2*ijlmn-1,index_enl,ijspin)
             enl_(2)=-enl(2*ijlmn  ,index_enl,ijspin)
           else
             jilmn=i0lmn+jlmn
             enl_(1:2)=enl(2*jilmn-1:2*jilmn,index_enl,jispin)
           end if
           gxi(1:cplex)=gx(1:cplex,ilmn,ia,1)
           gxfac_(1,jlmn,ia,jspinor_index)=gxfac_(1,jlmn,ia,jspinor_index)+enl_(1)*gxi(1)
           gxfac_(2,jlmn,ia,jspinor_index)=gxfac_(2,jlmn,ia,jspinor_index)+enl_(2)*gxi(1)
           if (cplex==2) then
             gxfac_(1,jlmn,ia,jspinor_index)=gxfac_(1,jlmn,ia,jspinor_index)-enl_(2)*gxi(2)
             gxfac_(2,jlmn,ia,jspinor_index)=gxfac_(2,jlmn,ia,jspinor_index)+enl_(1)*gxi(2)
           end if
         end do !ilmn
       end do !jlmn
!$OMP END DO
     end do !iat
!$OMP END PARALLEL
     call xmpi_sum(gxfac_,mpi_enreg%comm_spinor,ierr)
     gxfac(:,:,:,1)=gxfac(:,:,:,1)+gxfac_(:,:,:,ispinor_index)
     ABI_DEALLOCATE(gxfac_)
   end if

 end if !paw_opt

!Accumulate gxfac related to overlap (Sij) (PAW)
!------------------------------------------- ------------------------
 if (paw_opt==3.or.paw_opt==4) then                    ! Use Sij, overlap contribution
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

!Accumulate dgxdtfac related to nonlocal operator (Norm-conserving)
!-------------------------------------------------------------------
 if (optder>=1.and.paw_opt==0) then    ! Enl is E(Kleinman-Bylander)
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
         enl_(1)=enl(iln,itypat,ispinor_index)
         do mu=1,ndgxdtfac
           dgxdtfac(1:cplex,mu,ilmn,ia,ispinor)=enl_(1)*dgxdt(1:cplex,mu,ilmn,ia,ispinor)
         end do
       end do
!$OMP END DO
     end do
   end do
!$OMP END PARALLEL
 end if

!Accumulate dgxdtfac related to nonlocal operator (PAW)
!-------------------------------------------------------------------
 if (optder>=1.and.(paw_opt==1.or.paw_opt==2.or.paw_opt==4)) then  ! Enl is psp strength Dij
   dgxdtfac(1:cplex,1:ndgxdtfac,1:nlmn,1:nincat,1:nspinor)=zero
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
           enl_(1)=enl(jjlmn,index_enl,ispinor_index)
           if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
           do mu=1,ndgxdtfac
             gxfj(1:cplex,mu)=dgxdt(1:cplex,mu,jlmn,ia,ispinor)
             dgxdtfac(1:cplex,mu,jlmn,ia,ispinor)=dgxdtfac(1:cplex,mu,jlmn,ia,ispinor)+enl_(1)*gxfj(1:cplex,mu)
           end do
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1)=enl(ijlmn,index_enl,ispinor_index)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
             do mu=1,ndgxdtfac
               gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
               dgxdtfac(1:cplex,mu,jlmn,ia,ispinor)=dgxdtfac(1:cplex,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1:cplex)
#if !defined HAVE_OPENMP
               dgxdtfac(1:cplex,mu,ilmn,ia,ispinor)=dgxdtfac(1:cplex,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(1:cplex,mu)
#endif
             end do
           end do
#if defined HAVE_OPENMP
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=ilmn*(ilmn-1)/2
               ijlmn=i0lmn+jlmn
               enl_(1)=enl(ijlmn,index_enl,ispinor_index)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               do mu=1,ndgxdtfac
                 gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
                 dgxdtfac(1:cplex,mu,jlmn,ia,ispinor)=dgxdtfac(1:cplex,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1:cplex)
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
!    2-Enl is complex
   else
     ABI_CHECK(cplex_fac==cplex_enl,"BUG: invalid cplex_fac/=cplex_enl!")
     if (nspinortot==1) then    !===== when nspinor=1, D_ij=D_ji
!$OMP PARALLEL &
!$OMP PRIVATE(ia,index_enl,jlmn,j0lmn,jjlmn,enl_,mu,gxfj,ilmn,ijlmn,gxfi)
       ABI_ALLOCATE(gxfj,(cplex,ndgxdtfac))
       do ia=1,nincat
         index_enl=atindx1(iatm+ia)
!$OMP DO
         do jlmn=1,nlmn
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1:2)=enl(2*jjlmn-1:2*jjlmn,index_enl,1)
           if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
           do mu=1,ndgxdtfac
             if(cplex_dgxdt(mu)==2)then
               cplex_ = 2 ; gxfj(1,mu) = zero ; gxfj(2,mu) = dgxdt(1,mu,jlmn,ia,1)
             else
               cplex_ = cplex ; gxfj(1:cplex,mu)=dgxdt(1:cplex,mu,jlmn,ia,1)
             end if
             dgxdtfac(1,mu,jlmn,ia,1)=dgxdtfac(1,mu,jlmn,ia,1)+enl_(1)*gxfj(1,mu)
             dgxdtfac(2,mu,jlmn,ia,1)=dgxdtfac(2,mu,jlmn,ia,1)+enl_(2)*gxfj(1,mu)
             if (cplex_==2) then
               dgxdtfac(1,mu,jlmn,ia,1)=dgxdtfac(1,mu,jlmn,ia,1)-enl_(2)*gxfj(2,mu)
               dgxdtfac(2,mu,jlmn,ia,1)=dgxdtfac(2,mu,jlmn,ia,1)+enl_(1)*gxfj(2,mu)
             end if
           end do
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,1)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
             do mu=1,ndgxdtfac
               if(cplex_dgxdt(mu)==2)then
                 cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = dgxdt(1,mu,ilmn,ia,1)
               else
                 cplex_ = cplex ; gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,1)
               end if
               dgxdtfac(1,mu,jlmn,ia,1)=dgxdtfac(1,mu,jlmn,ia,1)+enl_(1)*gxfi(1)
               dgxdtfac(2,mu,jlmn,ia,1)=dgxdtfac(2,mu,jlmn,ia,1)+enl_(2)*gxfi(1)
#if !defined HAVE_OPENMP
               dgxdtfac(1,mu,ilmn,ia,1)=dgxdtfac(1,mu,ilmn,ia,1)+enl_(1)*gxfj(1,mu)
               dgxdtfac(2,mu,ilmn,ia,1)=dgxdtfac(2,mu,ilmn,ia,1)+enl_(2)*gxfj(1,mu)
#endif
               if (cplex_==2) then
                 dgxdtfac(1,mu,jlmn,ia,1)=dgxdtfac(1,mu,jlmn,ia,1)-enl_(2)*gxfi(2)
                 dgxdtfac(2,mu,jlmn,ia,1)=dgxdtfac(2,mu,jlmn,ia,1)+enl_(1)*gxfi(2)
#if !defined HAVE_OPENMP
                 dgxdtfac(1,mu,ilmn,ia,1)=dgxdtfac(1,mu,ilmn,ia,1)-enl_(2)*gxfj(2,mu)
                 dgxdtfac(2,mu,ilmn,ia,1)=dgxdtfac(2,mu,ilmn,ia,1)+enl_(1)*gxfj(2,mu)
#endif
               end if
             end do
           end do
#if defined HAVE_OPENMP
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=ilmn*(ilmn-1)/2
               ijlmn=i0lmn+jlmn
               enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,1)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               do mu=1,ndgxdtfac
                 if(cplex_dgxdt(mu)==2)then
                   cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = dgxdt(1,mu,ilmn,ia,1)
                 else
                   cplex_ = cplex ; gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,1)
                 end if
                 dgxdtfac(1,mu,jlmn,ia,1)=dgxdtfac(1,mu,jlmn,ia,1)+enl_(1)*gxfi(1)
                 dgxdtfac(2,mu,jlmn,ia,1)=dgxdtfac(2,mu,jlmn,ia,1)+enl_(2)*gxfi(1)
                 if (cplex_==2) then
                   dgxdtfac(1,mu,jlmn,ia,1)=dgxdtfac(1,mu,jlmn,ia,1)-enl_(2)*gxfi(2)
                   dgxdtfac(2,mu,jlmn,ia,1)=dgxdtfac(2,mu,jlmn,ia,1)+enl_(1)*gxfi(2)
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
     else                    !===== when nspinor=2, D_ij=D_ji^*
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
             enl_(1)=enl(2*jjlmn-1,index_enl,ispinor_index)  ! enl_ii is real for up-up and dn-dn
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
             do mu=1,ndgxdtfac
               if(cplex_dgxdt(mu)==2)then
                 cplex_ = 2 ; gxfj(1,mu) = zero ; gxfj(2,mu) = dgxdt(1,mu,jlmn,ia,ispinor)
               else
                 cplex_ = cplex ; gxfj(1:cplex,mu)=dgxdt(1:cplex,mu,jlmn,ia,ispinor)
               end if
               dgxdtfac(1:cplex_,mu,jlmn,ia,ispinor)=dgxdtfac(1:cplex_,mu,jlmn,ia,ispinor)+enl_(1)*gxfj(1:cplex_,mu)
             end do
             do ilmn=1,jlmn-1
               ijlmn=j0lmn+ilmn
               enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,ispinor_index)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               do mu=1,ndgxdtfac
                 if(cplex_dgxdt(mu)==2)then
                   cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = dgxdt(1,mu,ilmn,ia,ispinor)
                 else
                   cplex_ = cplex ; gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
                 end if
                 dgxdtfac(1,mu,jlmn,ia,ispinor)=dgxdtfac(1,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1)
                 dgxdtfac(2,mu,jlmn,ia,ispinor)=dgxdtfac(2,mu,jlmn,ia,ispinor)-enl_(2)*gxfi(1)
#if !defined HAVE_OPENMP
                 dgxdtfac(1,mu,ilmn,ia,ispinor)=dgxdtfac(1,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(1,mu)
                 dgxdtfac(2,mu,ilmn,ia,ispinor)=dgxdtfac(2,mu,ilmn,ia,ispinor)+enl_(2)*gxfj(1,mu)
#endif
                 if (cplex_==2) then
                   dgxdtfac(1,mu,jlmn,ia,ispinor)=dgxdtfac(1,mu,jlmn,ia,ispinor)+enl_(2)*gxfi(2)
                   dgxdtfac(2,mu,jlmn,ia,ispinor)=dgxdtfac(2,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(2)
#if !defined HAVE_OPENMP
                   dgxdtfac(1,mu,ilmn,ia,ispinor)=dgxdtfac(1,mu,ilmn,ia,ispinor)-enl_(2)*gxfj(2,mu)
                   dgxdtfac(2,mu,ilmn,ia,ispinor)=dgxdtfac(2,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(2,mu)
#endif
                 end if
               end do
             end do
#if defined HAVE_OPENMP
             if(jlmn<nlmn) then
               do ilmn=jlmn+1,nlmn
                 i0lmn=ilmn*(ilmn-1)/2
                 ijlmn=i0lmn+jlmn
                 enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,ispinor_index)
                 if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
                 do mu=1,ndgxdtfac
                   if(cplex_dgxdt(mu)==2)then
                     cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = dgxdt(1,mu,ilmn,ia,ispinor)
                   else
                     cplex_ = cplex ; gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
                   end if
                   dgxdtfac(1,mu,jlmn,ia,ispinor)=dgxdtfac(1,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1)
                   dgxdtfac(2,mu,jlmn,ia,ispinor)=dgxdtfac(2,mu,jlmn,ia,ispinor)+enl_(2)*gxfi(1)
                   if (cplex_==2) then
                     dgxdtfac(1,mu,jlmn,ia,ispinor)=dgxdtfac(1,mu,jlmn,ia,ispinor)-enl_(2)*gxfi(2)
                     dgxdtfac(2,mu,jlmn,ia,ispinor)=dgxdtfac(2,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(2)
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
           enl_(1:2)=enl(2*jjlmn-1:2*jjlmn,index_enl,2+ispinor)
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
             dgxdtfac(1,mu,jlmn,ia,jspinor)=dgxdtfac(1,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(1)
             dgxdtfac(2,mu,jlmn,ia,jspinor)=dgxdtfac(2,mu,jlmn,ia,jspinor)-enl_(2)*gxfi(1)
             if (cplex_==2) then
               dgxdtfac(2,mu,jlmn,ia,jspinor)=dgxdtfac(2,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(2)
               dgxdtfac(1,mu,jlmn,ia,jspinor)=dgxdtfac(1,mu,jlmn,ia,jspinor)+enl_(2)*gxfi(2)
             end if
           end do
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor)
             do mu=1,ndgxdtfac
               if(cplex_dgxdt(mu)==2)then
                 cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = dgxdt(1,mu,ilmn,ia,ispinor)
               else
                 cplex_ = cplex ; gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
               end if
               dgxdtfac(1,mu,jlmn,ia,jspinor)=dgxdtfac(1,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(1)
               dgxdtfac(2,mu,jlmn,ia,jspinor)=dgxdtfac(2,mu,jlmn,ia,jspinor)-enl_(2)*gxfi(1)
#if !defined HAVE_OPENMP
               dgxdtfac(1,mu,ilmn,ia,ispinor)=dgxdtfac(1,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(1,mu)
               dgxdtfac(2,mu,ilmn,ia,ispinor)=dgxdtfac(2,mu,ilmn,ia,ispinor)+enl_(2)*gxfj(1,mu)
#endif
               if (cplex_==2) then
                 dgxdtfac(1,mu,jlmn,ia,jspinor)=dgxdtfac(1,mu,jlmn,ia,jspinor)+enl_(2)*gxfi(2)
                 dgxdtfac(2,mu,jlmn,ia,jspinor)=dgxdtfac(2,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(2)
#if !defined HAVE_OPENMP
                 dgxdtfac(1,mu,ilmn,ia,ispinor)=dgxdtfac(1,mu,ilmn,ia,ispinor)-enl_(2)*gxfj(2,mu)
                 dgxdtfac(2,mu,ilmn,ia,ispinor)=dgxdtfac(2,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(2,mu)
#endif
               end if
             end do !mu
           end do !ilmn
#if defined HAVE_OPENMP
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=ilmn*(ilmn-1)/2
               ijlmn=i0lmn+jlmn
               enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor)
               do mu=1,ndgxdtfac
                 if(cplex_dgxdt(mu)==2)then
                   cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = dgxdt(1,mu,ilmn,ia,jspinor)
                 else
                   cplex_ = cplex ; gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,jspinor)
                 end if
                 dgxdtfac(1,mu,jlmn,ia,ispinor)=dgxdtfac(1,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1)
                 dgxdtfac(2,mu,jlmn,ia,ispinor)=dgxdtfac(2,mu,jlmn,ia,ispinor)+enl_(2)*gxfi(1)
                 if (cplex_==2) then
                   dgxdtfac(1,mu,jlmn,ia,ispinor)=dgxdtfac(1,mu,jlmn,ia,ispinor)-enl_(2)*gxfi(2)
                   dgxdtfac(2,mu,jlmn,ia,ispinor)=dgxdtfac(2,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(2)
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
     ABI_ALLOCATE(dgxdtfac_,(cplex_fac,ndgxdtfac,nlmn,nincat,nspinortot))
!$OMP PARALLEL &
!$OMP PRIVATE(ia,index_enl), &
!$OMP PRIVATE(jlmn,j0lmn,ilmn,i0lmn,ijlmn,enl_,jilmn,mu,gxfi)
!$OMP WORKSHARE
     dgxdtfac_(:,:,:,:,:)=zero
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
             enl_(1)= enl(2*ijlmn-1,index_enl,ijspin)
             enl_(2)=-enl(2*ijlmn  ,index_enl,ijspin)
           else
             jilmn=i0lmn+jlmn
             enl_(1:2)=enl(2*jilmn-1:2*jilmn,index_enl,jispin)
           end if
           do mu=1,ndgxdtfac
             if(cplex_dgxdt(mu)==2)then
               cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = dgxdt(1,mu,ilmn,ia,1)
             else
               cplex_ = cplex ; gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,1)
             end if
             dgxdtfac_(1,mu,jlmn,ia,jspinor_index)=dgxdtfac_(1,mu,jlmn,ia,jspinor_index)+enl_(1)*gxfi(1)
             dgxdtfac_(2,mu,jlmn,ia,jspinor_index)=dgxdtfac_(2,mu,jlmn,ia,jspinor_index)+enl_(2)*gxfi(1)
             if (cplex_==2) then
               dgxdtfac_(1,mu,jlmn,ia,jspinor_index)=dgxdtfac_(1,mu,jlmn,ia,jspinor_index)-enl_(2)*gxfi(2)
               dgxdtfac_(2,mu,jlmn,ia,jspinor_index)=dgxdtfac_(2,mu,jlmn,ia,jspinor_index)+enl_(1)*gxfi(2)
             end if
           end do
         end do !ilmn
       end do !jlmn
!$OMP END DO
     end do !iat
!$OMP SINGLE
     call xmpi_sum(dgxdtfac_,mpi_enreg%comm_spinor,ierr)
!$OMP END SINGLE
!$OMP WORKSHARE
     dgxdtfac(:,:,:,:,1)=dgxdtfac(:,:,:,:,1)+dgxdtfac_(:,:,:,:,ispinor_index)
!$OMP END WORKSHARE
!$OMP END PARALLEL
     ABI_DEALLOCATE(dgxdtfac_)
   end if !nspinortot

 end if ! pawopt & optder

!Accumulate dgxdtfac related to overlap (Sij) (PAW)
!-------------------------------------------------------------------
 if (optder>=1.and.(paw_opt==3.or.paw_opt==4)) then  ! Use Sij, overlap contribution
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

!Accumulate d2gxdtfac related to nonlocal operator (Norm-conserving)
!-------------------------------------------------------------------
 if (optder==2.and.paw_opt==0) then    ! Enl is E(Kleinman-Bylander)
!$OMP PARALLEL &
!$OMP PRIVATE(ispinor,ispinor_index,ia,ilmn,iln,enl_,mu)
   do ispinor=1,nspinor
     ispinor_index = ispinor + shift
     do ia=1,nincat
!$OMP DO
       do ilmn=1,nlmn
         iln=indlmn(5,ilmn)
         enl_(1)=enl(iln,itypat,ispinor_index)
         do mu=1,nd2gxdtfac
           d2gxdtfac(1:cplex,mu,ilmn,ia,ispinor)=enl_(1)*d2gxdt(1:cplex,mu,ilmn,ia,ispinor)
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
 if (optder==2.and.(paw_opt==1.or.paw_opt==2.or.paw_opt==4)) then  ! Enl is psp strength Dij
   d2gxdtfac(1:cplex,1:nd2gxdtfac,1:nlmn,1:nincat,1:nspinor)=zero
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
           enl_(1)=enl(jjlmn,index_enl,ispinor_index)
           if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
           do mu=1,nd2gxdtfac
             gxfj(1:cplex,mu)=d2gxdt(1:cplex,mu,jlmn,ia,ispinor)
             d2gxdtfac(1:cplex,mu,jlmn,ia,ispinor)=d2gxdtfac(1:cplex,mu,jlmn,ia,ispinor)+enl_(1)*gxfj(1:cplex,mu)
           end do
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1)=enl(ijlmn,index_enl,ispinor_index)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
             do mu=1,nd2gxdtfac
               gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,ispinor)
               d2gxdtfac(1:cplex,mu,jlmn,ia,ispinor)=d2gxdtfac(1:cplex,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1:cplex)
#if !defined HAVE_OPENMP
               d2gxdtfac(1:cplex,mu,ilmn,ia,ispinor)=d2gxdtfac(1:cplex,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(1:cplex,mu)
#endif
             end do
           end do
#if defined HAVE_OPENMP
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=ilmn*(ilmn-1)/2
               ijlmn=i0lmn+jlmn
               enl_(1)=enl(ijlmn,index_enl,ispinor_index)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               do mu=1,nd2gxdtfac
                 gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,ispinor)
                 d2gxdtfac(1:cplex,mu,jlmn,ia,ispinor)=d2gxdtfac(1:cplex,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1:cplex)
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
!    2-Enl is complex
   else
     ABI_CHECK(cplex_fac==cplex_enl,"BUG: invalid cplex_fac/=cplex_enl!")
     if (nspinortot==1) then    !===== when nspinor=1, D_ij=D_ji
!$OMP PARALLEL &
!$OMP PRIVATE(ia,index_enl,jlmn,j0lmn,jjlmn,enl_,mu,gxfj,ilmn,ijlmn,gxfi)
       ABI_ALLOCATE(gxfj,(cplex,nd2gxdtfac))
       do ia=1,nincat
         index_enl=atindx1(iatm+ia)
!$OMP DO
         do jlmn=1,nlmn
           j0lmn=jlmn*(jlmn-1)/2
           jjlmn=j0lmn+jlmn
           enl_(1:2)=enl(2*jjlmn-1:2*jjlmn,index_enl,1)
           if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
           do mu=1,nd2gxdtfac
             if(cplex_d2gxdt(mu)==2)then
               cplex_ = 2 ; gxfj(1,mu) = zero ; gxfj(2,mu) = d2gxdt(1,mu,jlmn,ia,1)
             else
               cplex_ = cplex ; gxfj(1:cplex,mu)=d2gxdt(1:cplex,mu,jlmn,ia,1)
             end if
             d2gxdtfac(1,mu,jlmn,ia,1)=d2gxdtfac(1,mu,jlmn,ia,1)+enl_(1)*gxfj(1,mu)
             d2gxdtfac(2,mu,jlmn,ia,1)=d2gxdtfac(2,mu,jlmn,ia,1)+enl_(2)*gxfj(1,mu)
             if (cplex_==2) then
               d2gxdtfac(1,mu,jlmn,ia,1)=d2gxdtfac(1,mu,jlmn,ia,1)-enl_(2)*gxfj(2,mu)
               d2gxdtfac(2,mu,jlmn,ia,1)=d2gxdtfac(2,mu,jlmn,ia,1)+enl_(1)*gxfj(2,mu)
             end if
           end do
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,1)
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
             do mu=1,nd2gxdtfac
               if(cplex_d2gxdt(mu)==2)then
                 cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = d2gxdt(1,mu,ilmn,ia,1)
               else
                 cplex_ = cplex ; gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,1)
               end if
               d2gxdtfac(1,mu,jlmn,ia,1)=d2gxdtfac(1,mu,jlmn,ia,1)+enl_(1)*gxfi(1)
               d2gxdtfac(2,mu,jlmn,ia,1)=d2gxdtfac(2,mu,jlmn,ia,1)+enl_(2)*gxfi(1)
#if !defined HAVE_OPENMP
               d2gxdtfac(1,mu,ilmn,ia,1)=d2gxdtfac(1,mu,ilmn,ia,1)+enl_(1)*gxfj(1,mu)
               d2gxdtfac(2,mu,ilmn,ia,1)=d2gxdtfac(2,mu,ilmn,ia,1)+enl_(2)*gxfj(1,mu)
#endif
               if (cplex_==2) then
                 d2gxdtfac(1,mu,jlmn,ia,1)=d2gxdtfac(1,mu,jlmn,ia,1)-enl_(2)*gxfi(2)
                 d2gxdtfac(2,mu,jlmn,ia,1)=d2gxdtfac(2,mu,jlmn,ia,1)+enl_(1)*gxfi(2)
#if !defined HAVE_OPENMP
                 d2gxdtfac(1,mu,ilmn,ia,1)=d2gxdtfac(1,mu,ilmn,ia,1)-enl_(2)*gxfj(2,mu)
                 d2gxdtfac(2,mu,ilmn,ia,1)=d2gxdtfac(2,mu,ilmn,ia,1)+enl_(1)*gxfj(2,mu)
#endif
               end if
             end do
           end do
#if defined HAVE_OPENMP
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=ilmn*(ilmn-1)/2
               ijlmn=i0lmn+jlmn
               enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,1)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               do mu=1,nd2gxdtfac
                 if(cplex_d2gxdt(mu)==2)then
                   cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = d2gxdt(1,mu,ilmn,ia,1)
                 else
                   cplex_ = cplex ; gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,1)
                 end if
                 d2gxdtfac(1,mu,jlmn,ia,1)=d2gxdtfac(1,mu,jlmn,ia,1)+enl_(1)*gxfi(1)
                 d2gxdtfac(2,mu,jlmn,ia,1)=d2gxdtfac(2,mu,jlmn,ia,1)+enl_(2)*gxfi(1)
                 if (cplex_==2) then
                   d2gxdtfac(1,mu,jlmn,ia,1)=d2gxdtfac(1,mu,jlmn,ia,1)-enl_(2)*gxfi(2)
                   d2gxdtfac(2,mu,jlmn,ia,1)=d2gxdtfac(2,mu,jlmn,ia,1)+enl_(1)*gxfi(2)
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
     else                    !===== when nspinor=2, D_ij=D_ji^*
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
             enl_(1)=enl(2*jjlmn-1,index_enl,ispinor_index)  ! enl_ii is real for up-up and dn-dn
             if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
             do mu=1,nd2gxdtfac
               if(cplex_d2gxdt(mu)==2)then
                 cplex_ = 2 ; gxfj(1,mu) = zero ; gxfj(2,mu) = d2gxdt(1,mu,jlmn,ia,ispinor)
               else
                 cplex_ = cplex ; gxfj(1:cplex,mu)=d2gxdt(1:cplex,mu,jlmn,ia,ispinor)
               end if
               d2gxdtfac(1:cplex_,mu,jlmn,ia,ispinor)=d2gxdtfac(1:cplex_,mu,jlmn,ia,ispinor)+enl_(1)*gxfj(1:cplex_,mu)
             end do
             do ilmn=1,jlmn-1
               ijlmn=j0lmn+ilmn
               enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,ispinor_index)
               if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
               do mu=1,nd2gxdtfac
                 if(cplex_d2gxdt(mu)==2)then
                   cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = d2gxdt(1,mu,ilmn,ia,ispinor)
                 else
                   cplex_ = cplex ; gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,ispinor)
                 end if
                 d2gxdtfac(1,mu,jlmn,ia,ispinor)=d2gxdtfac(1,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1)
                 d2gxdtfac(2,mu,jlmn,ia,ispinor)=d2gxdtfac(2,mu,jlmn,ia,ispinor)-enl_(2)*gxfi(1)
#if !defined HAVE_OPENMP
                 d2gxdtfac(1,mu,ilmn,ia,ispinor)=d2gxdtfac(1,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(1,mu)
                 d2gxdtfac(2,mu,ilmn,ia,ispinor)=d2gxdtfac(2,mu,ilmn,ia,ispinor)+enl_(2)*gxfj(1,mu)
#endif
                 if (cplex_==2) then
                   d2gxdtfac(1,mu,jlmn,ia,ispinor)=d2gxdtfac(1,mu,jlmn,ia,ispinor)+enl_(2)*gxfi(2)
                   d2gxdtfac(2,mu,jlmn,ia,ispinor)=d2gxdtfac(2,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(2)
#if !defined HAVE_OPENMP
                   d2gxdtfac(1,mu,ilmn,ia,ispinor)=d2gxdtfac(1,mu,ilmn,ia,ispinor)-enl_(2)*gxfj(2,mu)
                   d2gxdtfac(2,mu,ilmn,ia,ispinor)=d2gxdtfac(2,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(2,mu)
#endif
                 end if
               end do
             end do
#if defined HAVE_OPENMP
             if(jlmn<nlmn) then
               do ilmn=jlmn+1,nlmn
                 i0lmn=ilmn*(ilmn-1)/2
                 ijlmn=i0lmn+jlmn
                 enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,ispinor_index)
                 if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
                 do mu=1,nd2gxdtfac
                   if(cplex_d2gxdt(mu)==2)then
                     cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = d2gxdt(1,mu,ilmn,ia,ispinor)
                   else
                     cplex_ = cplex ; gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,ispinor)
                   end if
                   d2gxdtfac(1,mu,jlmn,ia,ispinor)=d2gxdtfac(1,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1)
                   d2gxdtfac(2,mu,jlmn,ia,ispinor)=d2gxdtfac(2,mu,jlmn,ia,ispinor)+enl_(2)*gxfi(1)
                   if (cplex_==2) then
                     d2gxdtfac(1,mu,jlmn,ia,ispinor)=d2gxdtfac(1,mu,jlmn,ia,ispinor)-enl_(2)*gxfi(2)
                     d2gxdtfac(2,mu,jlmn,ia,ispinor)=d2gxdtfac(2,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(2)
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
           enl_(1:2)=enl(2*jjlmn-1:2*jjlmn,index_enl,2+ispinor)
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
             d2gxdtfac(1,mu,jlmn,ia,jspinor)=d2gxdtfac(1,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(1)
             d2gxdtfac(2,mu,jlmn,ia,jspinor)=d2gxdtfac(2,mu,jlmn,ia,jspinor)-enl_(2)*gxfi(1)
             if (cplex_==2) then
               d2gxdtfac(2,mu,jlmn,ia,jspinor)=d2gxdtfac(2,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(2)
               d2gxdtfac(1,mu,jlmn,ia,jspinor)=d2gxdtfac(1,mu,jlmn,ia,jspinor)+enl_(2)*gxfi(2)
             end if
           end do
           do ilmn=1,jlmn-1
             ijlmn=j0lmn+ilmn
             enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor)
             do mu=1,nd2gxdtfac
               if(cplex_d2gxdt(mu)==2)then
                 cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = d2gxdt(1,mu,ilmn,ia,ispinor)
               else
                 cplex_ = cplex ; gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,ispinor)
               end if
               d2gxdtfac(1,mu,jlmn,ia,jspinor)=d2gxdtfac(1,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(1)
               d2gxdtfac(2,mu,jlmn,ia,jspinor)=d2gxdtfac(2,mu,jlmn,ia,jspinor)-enl_(2)*gxfi(1)
#if !defined HAVE_OPENMP
               d2gxdtfac(1,mu,ilmn,ia,ispinor)=d2gxdtfac(1,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(1,mu)
               d2gxdtfac(2,mu,ilmn,ia,ispinor)=d2gxdtfac(2,mu,ilmn,ia,ispinor)+enl_(2)*gxfj(1,mu)
#endif
               if (cplex_==2) then
                 d2gxdtfac(1,mu,jlmn,ia,jspinor)=d2gxdtfac(1,mu,jlmn,ia,jspinor)+enl_(2)*gxfi(2)
                 d2gxdtfac(2,mu,jlmn,ia,jspinor)=d2gxdtfac(2,mu,jlmn,ia,jspinor)+enl_(1)*gxfi(2)
#if !defined HAVE_OPENMP
                 d2gxdtfac(1,mu,ilmn,ia,ispinor)=d2gxdtfac(1,mu,ilmn,ia,ispinor)-enl_(2)*gxfj(2,mu)
                 d2gxdtfac(2,mu,ilmn,ia,ispinor)=d2gxdtfac(2,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(2,mu)
#endif
               end if
             end do !mu
           end do !ilmn
#if defined HAVE_OPENMP
           if(jlmn<nlmn) then
             do ilmn=jlmn+1,nlmn
               i0lmn=ilmn*(ilmn-1)/2
               ijlmn=i0lmn+jlmn
               enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor)
               do mu=1,nd2gxdtfac
                 if(cplex_d2gxdt(mu)==2)then
                   cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = d2gxdt(1,mu,ilmn,ia,jspinor)
                 else
                   cplex_ = cplex ; gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,jspinor)
                 end if
                 d2gxdtfac(1,mu,jlmn,ia,ispinor)=d2gxdtfac(1,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1)
                 d2gxdtfac(2,mu,jlmn,ia,ispinor)=d2gxdtfac(2,mu,jlmn,ia,ispinor)+enl_(2)*gxfi(1)
                 if (cplex_==2) then
                   d2gxdtfac(1,mu,jlmn,ia,ispinor)=d2gxdtfac(1,mu,jlmn,ia,ispinor)-enl_(2)*gxfi(2)
                   d2gxdtfac(2,mu,jlmn,ia,ispinor)=d2gxdtfac(2,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(2)
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
     ABI_ALLOCATE(d2gxdtfac_,(cplex_fac,nd2gxdtfac,nlmn,nincat,nspinortot))
!$OMP PARALLEL &
!$OMP PRIVATE(ia,index_enl), &
!$OMP PRIVATE(jlmn,j0lmn,ilmn,i0lmn,ijlmn,enl_,jilmn,mu,gxfi)
!$OMP WORKSHARE
     d2gxdtfac_(:,:,:,:,:)=zero
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
             enl_(1)= enl(2*ijlmn-1,index_enl,ijspin)
             enl_(2)=-enl(2*ijlmn  ,index_enl,ijspin)
           else
             jilmn=i0lmn+jlmn
             enl_(1:2)=enl(2*jilmn-1:2*jilmn,index_enl,jispin)
           end if
           do mu=1,nd2gxdtfac
             if(cplex_d2gxdt(mu)==2)then
               cplex_ = 2 ; gxfi(1) = zero ; gxfi(2) = d2gxdt(1,mu,ilmn,ia,1)
             else
               cplex_ = cplex ; gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,1)
             end if
             d2gxdtfac_(1,mu,jlmn,ia,jspinor_index)=d2gxdtfac_(1,mu,jlmn,ia,jspinor_index)+enl_(1)*gxfi(1)
             d2gxdtfac_(2,mu,jlmn,ia,jspinor_index)=d2gxdtfac_(2,mu,jlmn,ia,jspinor_index)+enl_(2)*gxfi(1)
             if (cplex_==2) then
               d2gxdtfac_(1,mu,jlmn,ia,jspinor_index)=d2gxdtfac_(1,mu,jlmn,ia,jspinor_index)-enl_(2)*gxfi(2)
               d2gxdtfac_(2,mu,jlmn,ia,jspinor_index)=d2gxdtfac_(2,mu,jlmn,ia,jspinor_index)+enl_(1)*gxfi(2)
             end if
           end do
         end do !ilmn
       end do !jlmn
!$OMP END DO
     end do !iat
!$OMP SINGLE
     call xmpi_sum(d2gxdtfac_,mpi_enreg%comm_spinor,ierr)
!$OMP END SINGLE
!$OMP WORKSHARE
     d2gxdtfac(:,:,:,:,1)=d2gxdtfac(:,:,:,:,1)+d2gxdtfac_(:,:,:,:,ispinor_index)
!$OMP END WORKSHARE
!$OMP END PARALLEL
     ABI_DEALLOCATE(d2gxdtfac_)
   end if !nspinortot

 end if ! pawopt & optder

!Accumulate d2gxdtfac related to overlap (Sij) (PAW)
!-------------------------------------------------------------------
 if (optder==2.and.(paw_opt==3.or.paw_opt==4)) then  ! Use Sij, overlap contribution
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
           d2gxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)=d2gxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)+sijr*gxfj(1:cplex,mu)
         end do
         do ilmn=1,jlmn-1
           ijlmn=j0lmn+ilmn
           sijr=sij(ijlmn)
           do mu=1,nd2gxdtfac
             gxfi(1:cplex)=d2gxdt(1:cplex,mu,ilmn,ia,ispinor)
             d2gxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)=d2gxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)+sijr*gxfi(1:cplex)
#if !defined HAVE_OPENMP
             d2gxdtfac_sij(1:cplex,mu,ilmn,ia,ispinor)=d2gxdtfac_sij(1:cplex,mu,ilmn,ia,ispinor)+sijr*gxfj(1:cplex,mu)
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
               d2gxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)=d2gxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)+sijr*gxfi(1:cplex)
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
