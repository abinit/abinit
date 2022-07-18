!!****m* ABINIT/m_opernlc_ylm_allwf_stdpar
!! NAME
!!  m_opernlc_ylm_allwf_stdpar
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

module m_opernlc_ylm_allwf_stdpar

  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi

  use defs_abitypes, only : MPI_type

  implicit none

  private
  !!***

  public :: opernlc_ylm_allwf_stdpar
  !!***

contains
  !!***

  !!****f* ABINIT/opernlc_ylm_allwf_stdpar
  !! NAME
  !! opernlc_ylm_allwf_stdpar
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
  !!      m_gemm_nonlop,m_nonlop_ylm
  !!
  !! CHILDREN
  !!      xmpi_sum
  !!
  !! SOURCE

  subroutine opernlc_ylm_allwf_stdpar(atindx1,cplex,cplex_enl,cplex_fac, &
    &          dimenl1,dimenl2,dimekbq,enl, &
    &          gx,gxfac,gxfac_sij, &
    &          iatm,indlmn,itypat,ndat,lambda,mpi_enreg,natom, &
    &          nincat,nlmn,nspinor,nspinortot,paw_opt,sij)

    !Arguments ------------------------------------
    !scalars
    integer,intent(in)    :: cplex,cplex_enl,cplex_fac,dimenl1,dimenl2,dimekbq,iatm,itypat
    integer,intent(in)    :: natom,nincat,nspinor,nspinortot,paw_opt
    integer,intent(inout) :: nlmn
    type(MPI_type) , intent(in) :: mpi_enreg
    integer,intent(in)    :: ndat
    !arrays
    integer, intent(in)         :: atindx1(natom),indlmn(6,nlmn)
    real(dp),intent(in),target  :: enl(dimenl1,dimenl2,nspinortot**2,dimekbq)
    real(dp),intent(inout)      :: gx(cplex,nlmn,nincat,nspinor*ndat)
    real(dp),intent(in)         :: sij(((paw_opt+1)/3)*nlmn*(nlmn+1)/2)
    real(dp),intent(out),target :: gxfac(cplex_fac,nlmn,nincat,nspinor*ndat)
    real(dp),intent(out)        :: gxfac_sij(cplex,nlmn,nincat,nspinor*ndat*(paw_opt/3))
    real(dp),intent(in)         :: lambda(ndat)

    !Local variables-------------------------------
    !Arrays
    !scalars
    integer  :: cplex_, ia, idat, idat_ispinor, idat_jspinor, ierr
    integer  :: ijlmn,ijspin,ilm,ilmn,i0lmn,iln,index_enl,iphase,ispinor,ispinor_index
    integer  :: j0lmn,jilmn,jispin,jjlmn,jlm,jlmn,jspinor,jspinor_index,mu,shift
    real(dp) :: sijr
    !arrays
    real(dp)                    :: enl_(2), gxfi(2), gxi(cplex), gxj(cplex)
    real(dp),allocatable        :: gxfac_offdiag(:,:,:,:),gxfj(:,:)
    real(dp),pointer            :: gxfac_(:,:,:,:)
    real(dp),pointer            :: enl_ptr(:,:,:)

    ! *************************************************************************

    DBG_ENTER("COLL")

    !Parallelization over spinors treatment
    shift=0
    if (mpi_enreg%paral_spinor==1) then
      shift = mpi_enreg%me_spinor
    end if

    !When Enl factors contain a exp(-iqR) phase:
    ! - We loop over the real and imaginary parts
    ! - We need an additional memory space
    do iphase=1,dimekbq

      if (paw_opt==3) cycle

      if (iphase==1) then
        gxfac_ => gxfac
      else
        ABI_CHECK(cplex_fac==2,"BUG: invalid cplex_fac==1 when dimekbq=2!")
        ABI_MALLOC(gxfac_,(cplex_fac,nlmn,nincat,ndat*nspinor))
      end if

      gxfac_ = zero
      enl_ptr => enl(:,:,:,iphase)


      !Accumulate gxfac related to non-local operator (Norm-conserving)
      !-------------------------------------------------------------------
      if (paw_opt==0) then
        !Enl is E(Kleinman-Bylander)
        ABI_CHECK(cplex_enl/=2,"BUG: invalid cplex_enl=2!")
        ABI_CHECK(cplex_fac==cplex,"BUG: invalid cplex_fac/=cplex!")

        do concurrent(idat = 1:ndat, ispinor = 1:nspinor, ia = 1:nincat, ilmn = 1:nlmn) &
          & local(ispinor_index, idat_ispinor, iln, enl_) shared(shift)

          ispinor_index = ispinor+shift
          idat_ispinor = nspinor*(idat-1) + ispinor
          iln = indlmn(5,ilmn)
          enl_(1) = enl_ptr(iln,itypat,ispinor_index)
          gxfac_(1:cplex,ilmn,ia,idat_ispinor) = enl_(1)*gx(1:cplex,ilmn,ia,idat_ispinor)

        end do

      end if ! paw_opt == 0

      !Accumulate gxfac related to nonlocal operator (PAW)
      !-------------------------------------------------------------------
      if (paw_opt==1.or.paw_opt==2.or.paw_opt==4) then

        !Enl is psp strength Dij or (Dij-lambda.Sij)

        !  === Diagonal term(s) (up-up, down-down)

        if (cplex_enl==1) then
          !  1-Enl is real

          do concurrent(idat = 1:ndat, ispinor=1:nspinor, ia=1:nincat, jlmn=1:nlmn) &
            & local(ispinor_index, index_enl, i0lmn, j0lmn, ijlmn, jjlmn, enl_, idat_ispinor, gxi, gxj) &
            & shared(paw_opt, iatm, lambda, sij)

            ispinor_index = ispinor+shift
            index_enl = atindx1(iatm+ia)
            j0lmn = jlmn*(jlmn-1)/2
            jjlmn = j0lmn+jlmn
            enl_(1) = enl_ptr(jjlmn,index_enl,ispinor_index)
            idat_ispinor = nspinor*(idat-1) + ispinor

            if (paw_opt==2) enl_(1) = enl_(1)-lambda(idat)*sij(jjlmn)
            gxj   (1:cplex)                      = gx    (1:cplex,jlmn,ia,idat_ispinor)
            gxfac_(1:cplex,jlmn,ia,idat_ispinor) = gxfac_(1:cplex,jlmn,ia,idat_ispinor) + enl_(1)*gxj(1:cplex)

            do ilmn=1,nlmn
              if (ilmn < jlmn) then
                ijlmn = j0lmn+ilmn
                enl_(1) = enl_ptr(ijlmn,index_enl,ispinor_index)
                if (paw_opt==2) enl_(1)=enl_(1)-lambda(idat)*sij(ijlmn)
                gxi   (1:cplex)                      = gx    (1:cplex,ilmn,ia,idat_ispinor)
                gxfac_(1:cplex,jlmn,ia,idat_ispinor) = gxfac_(1:cplex,jlmn,ia,idat_ispinor) + enl_(1)*gxi(1:cplex)
              else if (ilmn > jlmn) then
                if(jlmn<nlmn) then
                  i0lmn = (ilmn*(ilmn-1)/2)
                  ijlmn = i0lmn+jlmn
                  enl_(1) = enl_ptr(ijlmn,index_enl,ispinor_index)
                  if (paw_opt==2) enl_(1)=enl_(1)-lambda(idat)*sij(ijlmn)
                  gxi   (1:cplex)                      = gx    (1:cplex,ilmn,ia,idat_ispinor)
                  gxfac_(1:cplex,jlmn,ia,idat_ispinor) = gxfac_(1:cplex,jlmn,ia,idat_ispinor) + enl_(1)*gxi(1:cplex)
                end if
              end if
            end do  ! ilmn

          end do ! concurrent

        else
          !  2-Enl is complex  ===== D^ss'_ij=D^s's_ji^*

          ABI_CHECK(cplex_fac==cplex_enl,"BUG: invalid cplex_fac/=cplex_enl!")

          if (nspinortot==1) then ! -------------> NO SPINORS

            do concurrent(idat = 1:ndat, ia=1:nincat, jlmn=1:nlmn) &
              & local(idat_ispinor, index_enl, i0lmn, j0lmn, ijlmn, jjlm, enl_, gxi, gxj) &
              & shared(iatm, nlmn)

              ! ispinor is 1
              ! nspinor should be 1 too (maybe an ABI_CHECK could ensure that)
              idat_ispinor = nspinor*(idat-1) + 1

              index_enl = atindx1(iatm+ia)

              j0lmn = jlmn*(jlmn-1)/2
              jjlmn = j0lmn+jlmn

              enl_(1) = enl_ptr(2*jjlmn-1,index_enl,1)
              if (paw_opt==2) then
                enl_(1) = enl_(1)-lambda(idat)*sij(jjlmn)
              end if

              gxj(1:cplex) = gx(1:cplex,jlmn,ia,idat_ispinor)

              gxfac_(1,jlmn,ia,idat_ispinor) = gxfac_(1,jlmn,ia,idat_ispinor) + enl_(1)*gxj(1)
              if (cplex==2) then
                gxfac_(2,jlmn,ia,idat_ispinor) = gxfac_(2,jlmn,ia,idat_ispinor) + enl_(1)*gxj(2)
              end if

              do ilmn=1,jlmn-1

                if (ilmn < jlmn) then

                  ijlmn = j0lmn + ilmn

                  enl_(1:2) = enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,1)
                  if (paw_opt==2) enl_(1)=enl_(1)-lambda(idat)*sij(ijlmn)

                  gxi(1:cplex) = gx(1:cplex,ilmn,ia,idat_ispinor)

                  gxfac_(1,jlmn,ia,idat_ispinor) = gxfac_(1,jlmn,ia,idat_ispinor) + enl_(1)*gxi(1)
                  gxfac_(2,jlmn,ia,idat_ispinor) = gxfac_(2,jlmn,ia,idat_ispinor) - enl_(2)*gxi(1)

                  if (cplex==2) then
                    gxfac_(1,jlmn,ia,idat_ispinor) = gxfac_(1,jlmn,ia,idat_ispinor) + enl_(2)*gxi(2)
                    gxfac_(2,jlmn,ia,idat_ispinor) = gxfac_(2,jlmn,ia,idat_ispinor) + enl_(1)*gxi(2)
                  end if

                else if (ilmn > jlmn) then

                  if(jlmn<nlmn) then
                    i0lmn=ilmn*(ilmn-1)/2
                    ijlmn=i0lmn+jlmn

                    enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,1)
                    if (paw_opt==2) enl_(1)=enl_(1)-lambda(idat)*sij(ijlmn)

                    gxi(1:cplex) = gx(1:cplex,ilmn,ia,idat_ispinor)

                    gxfac_(1,jlmn,ia,idat_ispinor) = gxfac_(1,jlmn,ia,idat_ispinor) + enl_(1)*gxi(1)
                    gxfac_(2,jlmn,ia,idat_ispinor) = gxfac_(2,jlmn,ia,idat_ispinor) + enl_(2)*gxi(1)

                    if (cplex==2) then
                      gxfac_(1,jlmn,ia,idat_ispinor) = gxfac_(1,jlmn,ia,idat_ispinor) - enl_(2)*gxi(2)
                      gxfac_(2,jlmn,ia,idat_ispinor) = gxfac_(2,jlmn,ia,idat_ispinor) + enl_(1)*gxi(2)
                    end if ! cplex == 2

                  end if ! jlmn

                end if ! ilmn

              end do ! ilmn

            end do ! concurrent

          else
            ! -------------> SPINORIAL CASE

            !  === Diagonal term(s) (up-up, down-down)

            do concurrent(idat = 1:ndat, ispinor = 1:nspinor, ia=1:nincat, jlmn = 1:nlmn) &
              & local(ispinor_index, idat_ispinor, index_enl, enl_, gxi, gxj, i0lmn, j0lmn, ijlmn, jjlmn) &
              & shared(shift, nspinor, iatm)

              ispinor_index = ispinor+shift
              idat_ispinor = nspinor*(idat-1) + ispinor

              index_enl = atindx1(iatm+ia)

              j0lmn = jlmn*(jlmn-1)/2
              jjlmn = j0lmn+jlmn

              enl_(1) = enl_ptr(2*jjlmn-1,index_enl,ispinor_index)
              if (paw_opt==2) then
                enl_(1) = enl_(1)-lambda(idat)*sij(jjlmn)
              end if

              gxj(1:cplex) = gx(1:cplex,jlmn,ia,idat_ispinor)

              gxfac_(1,jlmn,ia,idat_ispinor) = gxfac_(1,jlmn,ia,idat_ispinor) + enl_(1)*gxj(1)
              if (cplex==2) then
                gxfac_(2,jlmn,ia,idat_ispinor) = gxfac_(2,jlmn,ia,idat_ispinor) + enl_(1)*gxj(2)
              end if

              do ilmn = 1,nlmn

                if (ilmn < jlmn) then
                  ijlmn = j0lmn+ilmn

                  enl_(1:2) = enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,ispinor_index)
                  if (paw_opt==2) then
                    enl_(1) = enl_(1)-lambda(idat)*sij(ijlmn)
                  end if

                  gxi(1:cplex) = gx(1:cplex,ilmn,ia,idat_ispinor)

                  gxfac_(1,jlmn,ia,idat_ispinor) = gxfac_(1,jlmn,ia,idat_ispinor) + enl_(1)*gxi(1)
                  gxfac_(2,jlmn,ia,idat_ispinor) = gxfac_(2,jlmn,ia,idat_ispinor) - enl_(2)*gxi(1)

                  if (cplex==2) then
                    gxfac_(1,jlmn,ia,idat_ispinor) = gxfac_(1,jlmn,ia,idat_ispinor) + enl_(2)*gxi(2)
                    gxfac_(2,jlmn,ia,idat_ispinor) = gxfac_(2,jlmn,ia,idat_ispinor) + enl_(1)*gxi(2)
                  end if

                else if (ilmn > jlmn) then

                  if(jlmn<nlmn) then
                    i0lmn=ilmn*(ilmn-1)/2
                    ijlmn=i0lmn+jlmn
                    enl_(1:2) = enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,ispinor_index)

                    if (paw_opt==2) enl_(1)=enl_(1)-lambda(idat)*sij(ijlmn)
                    gxi(1:cplex) = gx(1:cplex,ilmn,ia,idat_ispinor)
                    gxfac_(1,jlmn,ia,idat_ispinor) = gxfac_(1,jlmn,ia,idat_ispinor) + enl_(1)*gxi(1)
                    gxfac_(2,jlmn,ia,idat_ispinor) = gxfac_(2,jlmn,ia,idat_ispinor) + enl_(2)*gxi(1)
                    if (cplex==2) then
                      gxfac_(1,jlmn,ia,idat_ispinor) = gxfac_(1,jlmn,ia,idat_ispinor) - enl_(2)*gxi(2)
                      gxfac_(2,jlmn,ia,idat_ispinor) = gxfac_(2,jlmn,ia,idat_ispinor) + enl_(1)*gxi(2)
                    end if ! cplex == 2
                  end if

                end if ! ilmn

              end do ! ilmn

            end do  ! concurrent

          end if ! nspinortot

        end if ! complex_enl

        !  === Off-diagonal term(s) (up-down, down-up)

        !  --- No parallelization over spinors ---
        if (nspinortot==2 .and. nspinor==nspinortot) then

          ABI_CHECK(cplex_enl==2,"BUG: invalid cplex_enl/=2!")
          ABI_CHECK(cplex_fac==cplex,"BUG: invalid cplex_fac/=cplex)!")

          do concurrent(idat = 1:ndat, ispinor=1:nspinortot, ia=1:nincat, jlmn=1:nlmn) &
            & local(jspinor, idat_ispinor, idat_jspinor, index_enl, gxi, gxj, i0lmn, j0lmn, ijlmn, jjlmn) &
            & shared(nspinor,iatm)

            jspinor = 3-ispinor

            idat_ispinor = nspinor*(idat-1) + ispinor

            ! watch out : ispinor replaced by jspinor
            idat_jspinor = nspinor*(idat-1) + jspinor

            index_enl=atindx1(iatm+ia)

            j0lmn=jlmn*(jlmn-1)/2
            jjlmn=j0lmn+jlmn
            enl_(1:2) = enl_ptr(2*jjlmn-1:2*jjlmn,index_enl,2+ispinor )
            gxi(1:cplex) = gx(1:cplex,jlmn,ia,idat_ispinor)
            gxfac_(1,jlmn,ia,idat_jspinor) = gxfac_(1,jlmn,ia,idat_jspinor) + enl_(1)*gxi(1)
            gxfac_(2,jlmn,ia,idat_jspinor) = gxfac_(2,jlmn,ia,idat_jspinor) - enl_(2)*gxi(1)
            if (cplex==2) then
              gxfac_(1,jlmn,ia,idat_jspinor) = gxfac_(1,jlmn,ia,idat_jspinor) + enl_(2)*gxi(2)
              gxfac_(2,jlmn,ia,idat_jspinor) = gxfac_(2,jlmn,ia,idat_jspinor) + enl_(1)*gxi(2)
            end if

            do ilmn=1,nlmn

              if (ilmn < jlmn) then

                ijlmn=j0lmn+ilmn
                enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor)
                gxi(1:cplex) = gx(1:cplex,ilmn,ia,idat_ispinor)
                gxfac_(1,jlmn,ia,idat_jspinor) = gxfac_(1,jlmn,ia,idat_jspinor) + enl_(1)*gxi(1)
                gxfac_(2,jlmn,ia,idat_jspinor) = gxfac_(2,jlmn,ia,idat_jspinor) - enl_(2)*gxi(1)
                if (cplex==2) then
                  gxfac_(1,jlmn,ia,idat_jspinor) = gxfac_(1,jlmn,ia,idat_jspinor) + enl_(2)*gxi(2)
                  gxfac_(2,jlmn,ia,idat_jspinor) = gxfac_(2,jlmn,ia,idat_jspinor) + enl_(1)*gxi(2)
                end if

              else if (ilmn > jlmn) then

                if(jlmn<nlmn) then
                  i0lmn=ilmn*(ilmn-1)/2
                  ijlmn=i0lmn+jlmn
                  enl_(1:2)=enl_ptr(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor)
                  gxi(1:cplex)=gx(1:cplex,ilmn,ia,idat_jspinor)
                  gxfac_(1,jlmn,ia,idat_ispinor) = gxfac_(1,jlmn,ia,idat_ispinor) + enl_(1)*gxi(1)
                  gxfac_(2,jlmn,ia,idat_ispinor) = gxfac_(2,jlmn,ia,idat_ispinor) + enl_(2)*gxi(1)
                  if (cplex==2) then
                    gxfac_(1,jlmn,ia,idat_ispinor) = gxfac_(1,jlmn,ia,idat_ispinor) - enl_(2)*gxi(2)
                    gxfac_(2,jlmn,ia,idat_ispinor) = gxfac_(2,jlmn,ia,idat_ispinor) + enl_(1)*gxi(2)
                  end if ! cplex==2

                end if ! if jlmn

              end if ! if ilmn

            end do ! ilmn

          end do ! concurrent

          !    --- Parallelization over spinors ---

          !
          ! this case is removed for GPU parallelization, npspinor>1 is forbidden when use_gpu_cuda is ON
          !

        ! else if (nspinortot==2 .and. nspinor/=nspinortot) then

        !   ABI_CHECK(cplex_enl==2,"BUG: invalid cplex_enl/=2!")
        !   ABI_CHECK(cplex_fac==2,"BUG: invalid cplex_fac/=2!")
        !   ABI_MALLOC(gxfac_offdiag,(cplex_fac,nlmn,nincat,nspinortot))

        !   gxfac_offdiag(:,:,:,:) = zero

        !   ispinor_index = mpi_enreg%me_spinor+1
        !   jspinor_index = 3-ispinor_index
        !   if (ispinor_index==1) then
        !     ijspin=3
        !     jispin=4
        !   else
        !     ijspin=4
        !     jispin=3
        !   end if

        !   do ia=1,nincat
        !     index_enl = atindx1(iatm+ia)

        !     do jlmn=1,nlmn
        !       j0lmn=jlmn*(jlmn-1)/2

        !       do ilmn=1,nlmn
        !         i0lmn=ilmn*(ilmn-1)/2

        !         if (ilmn<=jlmn) then
        !           ijlmn=j0lmn+ilmn
        !           enl_(1) =  enl_ptr(2*ijlmn-1,index_enl,ijspin)
        !           enl_(2) = -enl_ptr(2*ijlmn  ,index_enl,ijspin)
        !         else
        !           jilmn=i0lmn+jlmn
        !           enl_(1) = enl_ptr(2*jilmn-1,index_enl,jispin)
        !           enl_(2) = enl_ptr(2*jilmn  ,index_enl,jispin)
        !         end if

        !         gxi(1:cplex) = gx(1:cplex,ilmn,ia,1)

        !         gxfac_offdiag(1,jlmn,ia,jspinor_index) = &
        !           &            gxfac_offdiag(1,jlmn,ia,jspinor_index) + enl_(1)*gxi(1)

        !         gxfac_offdiag(2,jlmn,ia,jspinor_index) = &
        !           &            gxfac_offdiag(2,jlmn,ia,jspinor_index) + enl_(2)*gxi(1)

        !         if (cplex==2) then
        !           gxfac_offdiag(1,jlmn,ia,jspinor_index) = &
        !             &              gxfac_offdiag(1,jlmn,ia,jspinor_index) - enl_(2)*gxi(2)
        !           gxfac_offdiag(2,jlmn,ia,jspinor_index) = &
        !             &              gxfac_offdiag(2,jlmn,ia,jspinor_index) + enl_(1)*gxi(2)
        !         end if

        !       end do !ilmn

        !     end do !jlmn

        !   end do !iat

        !   !call xmpi_sum(gxfac_offdiag,mpi_enreg%comm_spinor,ierr)

        !   gxfac_(:,:,:,1) = gxfac_(:,:,:,1) + gxfac_offdiag(:,:,:,ispinor_index)
        !   ABI_FREE(gxfac_offdiag)

        end if

      end if !paw_opt

      !End of loop when a exp(-iqR) phase is present
      !------------------------------------------- ------------------------

      !When iphase=1, gxfac and gxfac_ point to the same memory space
      !When iphase=2, we add i.gxfac_ to gxfac
      if (iphase==2) then

        do concurrent(idat = 1:ndat, ispinor=1:nspinor, ia=1:nincat, ilmn=1:nlmn) &
          & local(idat_ispinor) &
          & shared(nspinor)

          idat_ispinor = nspinor*(idat-1) + ispinor
          gxfac(1,ilmn,ia,idat_ispinor) = gxfac(1,ilmn,ia,idat_ispinor) - gxfac_(2,ilmn,ia,idat_ispinor)
          gxfac(2,ilmn,ia,idat_ispinor) = gxfac(2,ilmn,ia,idat_ispinor) + gxfac_(1,ilmn,ia,idat_ispinor)

        end do ! concurrent

        ABI_FREE(gxfac_)

      end if ! iphase==2

      !End loop over real/imaginary part of the exp(-iqR) phase
    end do


    !Accumulate gxfac related to overlap (Sij) (PAW)
    !------------------------------------------- ------------------------
    if (paw_opt==3.or.paw_opt==4) then ! Use Sij, overlap contribution

      gxfac_sij(1:cplex,1:nlmn,1:nincat,1:nspinor) = zero

      do concurrent(idat = 1:ndat, ispinor = 1:nspinor, ia = 1:nincat, jlmn = 1:nlmn) &
        & local(idat_ispinor, i0lmn, j0lmn, ijlmn, jjlmn, sijr, gxi, gxj, ilm) &
        & shared(nspinor)

        idat_ispinor = nspinor*(idat-1) + ispinor

        j0lmn=jlmn*(jlmn-1)/2
        jjlmn=j0lmn+jlmn
        jlm=indlmn(4,jlmn)
        sijr=sij(jjlmn)
        gxj(1:cplex) = gx(1:cplex,jlmn,ia,idat_ispinor)
        gxfac_sij(1:cplex,jlmn,ia,idat_ispinor) = gxfac_sij(1:cplex,jlmn,ia,idat_ispinor) + sijr*gxj(1:cplex)

        do ilmn = 1,jlmn-1
          ilm=indlmn(4,ilmn)
          !if (ilm==jlm) then
          ijlmn=j0lmn+ilmn
          sijr=sij(ijlmn)
          gxi(1:cplex) = gx(1:cplex,ilmn,ia,idat_ispinor)
          gxfac_sij(1:cplex,jlmn,ia,idat_ispinor) = gxfac_sij(1:cplex,jlmn,ia,idat_ispinor) + sijr*gxi(1:cplex)
          !end if
        end do

        if(jlmn<nlmn) then

          do ilmn = jlmn+1,nlmn
            ilm=indlmn(4,ilmn)
            !if (ilm==jlm) then
            i0lmn=ilmn*(ilmn-1)/2
            ijlmn=i0lmn+jlmn
            sijr=sij(ijlmn)
            gxi(1:cplex) = gx(1:cplex,ilmn,ia,idat_ispinor)
            gxfac_sij(1:cplex,jlmn,ia,idat_ispinor) = gxfac_sij(1:cplex,jlmn,ia,idat_ispinor) + sijr*gxi(1:cplex)
            !end if
          end do ! ilmn

        end if ! jlmn

      end do ! concurrent

    end if ! paw_opt

  end subroutine opernlc_ylm_allwf_stdpar
  !!***

end module m_opernlc_ylm_allwf_stdpar
!!***
