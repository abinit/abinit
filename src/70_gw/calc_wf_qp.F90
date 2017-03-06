!{\src2tex{textfont=tt}}
!!****f* ABINIT/update_cprj
!! NAME
!! update_cprj
!!
!! FUNCTION
!!  Update the matrix elements of the PAW projectors in case of self-consistent GW.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2017 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dimlmn(natom)=number of (l,m,n) components for each atom (only for PAW)
!!  nkibz=number of k-points
!!  nsppol=number of spin
!!  nbnds=number of bands in the present GW calculation
!!  m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)= expansion of the QP amplitudes in terms of KS wavefunctions
!!  natom=number of atomd in unit cell
!!
!! OUTPUT
!!  Cprj_ibz(natom,nspinor*nkibz*nbnds*nsppol) <type(pawcprj_type)>=projected wave functions 
!!   <Proj_i|Cnk> with all NL projectors. On exit, it contains the projections onto the 
!!   QP amplitudes.
!!
!! TODO 
!! To be moved to cprj_utils, although here we use complex variables.
!!
!! PARENTS
!!      mlwfovlp_qp
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine update_cprj(natom,nkibz,nbnds,nsppol,nspinor,m_lda_to_qp,dimlmn,Cprj_ibz)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_pawcprj, only : pawcprj_type, pawcprj_alloc, pawcprj_free

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'update_cprj'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nbnds,nkibz,nsppol,nspinor
!arrays
 integer,intent(in) :: dimlmn(natom)
 complex(dpc),intent(in) :: m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)
 type(pawcprj_type),intent(inout) :: Cprj_ibz(natom,nspinor*nbnds*nkibz*nsppol)

!Local variables-------------------------------
!scalars
 integer :: iat,ib,ik,is,shift,indx_kibz,ilmn,nlmn,ispinor,ibsp,spad,ibdx
!arrays
 real(dp),allocatable :: re_p(:),im_p(:),vect(:,:),umat(:,:,:)
 type(pawcprj_type),allocatable :: Cprj_ks(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 ABI_DATATYPE_ALLOCATE(Cprj_ks,(natom,nspinor*nbnds))
 call pawcprj_alloc(Cprj_ks,0,dimlmn)

 ABI_ALLOCATE(re_p,(nbnds))
 ABI_ALLOCATE(im_p,(nbnds))
 ABI_ALLOCATE(vect,(2,nbnds))
 ABI_ALLOCATE(umat,(2,nbnds,nbnds))
 !
 ! $ \Psi^{QP}_{r,b} = \sum_n \Psi^{KS}_{r,n} M_{n,b} $ 
 !
 ! therefore the updated PAW projections are given by:
 !
 ! $ \<\tprj_j|\Psi^{QP}_a\> = sum_b M_{b,a} <\tprj_j|\Psi^{KS}_b\> $.
 !
 do is=1,nsppol
   do ik=1,nkibz

    shift=nspinor*nbnds*nkibz*(is-1)
    indx_kibz=nspinor*nbnds*(ik-1)+shift
    ibsp=0
    do ib=1,nbnds
      do ispinor=1,nspinor
        ibsp=ibsp+1
        do iat=1,natom
          Cprj_ks(iat,ibsp)%cp(:,:)=Cprj_ibz(iat,indx_kibz+ibsp)%cp(:,:)
        end do
      end do
    end do
    
    umat(1,:,:)=TRANSPOSE( REAL (m_lda_to_qp(:,:,ik,is)) )
    umat(2,:,:)=TRANSPOSE( AIMAG(m_lda_to_qp(:,:,ik,is)) )

    do iat=1,natom
      nlmn=dimlmn(iat)
      do ilmn=1,nlmn

        do ispinor=1,nspinor
           ! * Retrieve projections for this spinor component, at fixed atom and ilmn.
           spad=(ispinor-1)
           ibdx=0
           do ib=1,nbnds*nspinor,nspinor 
            ibdx=ibdx+1
            vect(1,ibdx)=Cprj_ks(iat,ib+spad)%cp(1,ilmn)
            vect(2,ibdx)=Cprj_ks(iat,ib+spad)%cp(2,ilmn)
           end do

           re_p(:)= &
&            MATMUL(umat(1,:,:),vect(1,:)) &
&           -MATMUL(umat(2,:,:),vect(2,:))

           im_p(:)= &
&            MATMUL(umat(1,:,:),vect(2,:)) &
&           +MATMUL(umat(2,:,:),vect(1,:))
           
           ! === Save values ===
           ibdx=0
           do ib=1,nbnds*nspinor,nspinor 
            ibdx=ibdx+1
            Cprj_ibz(iat,indx_kibz+spad+ib)%cp(1,ilmn)=re_p(ibdx)
            Cprj_ibz(iat,indx_kibz+spad+ib)%cp(2,ilmn)=im_p(ibdx)
           end do
        end do !ispinor

      end do !ilmn
    end do !iat

   end do !ik
 end do !is

 ABI_DEALLOCATE(re_p)
 ABI_DEALLOCATE(im_p)
 ABI_DEALLOCATE(vect)
 ABI_DEALLOCATE(umat)

 call pawcprj_free(Cprj_ks) 
 ABI_DATATYPE_DEALLOCATE(Cprj_ks)

 DBG_EXIT("COLL")

end subroutine update_cprj
!!***
