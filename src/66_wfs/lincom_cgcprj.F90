!{\src2tex{textfont=tt}}
!!****f* ABINIT/lincom_cgcprj
!!
!! NAME
!! lincom_cgcprj
!!
!! FUNCTION
!! For one k point and spinpol, compute a set (size nband_out) of linear combinations of nband_in wavefunctions,
!! that are known in the cg+cprj representation :
!! cgout_n(:,:) <--- Sum_m [ cg_m(:,:) . alpha_mn ]
!! cprjout_n(:,:) <--- Sum_n [ cprj_m(:,:) . alpha_mn ]
!! If nband_out is smaller or equal to nband_in, the result might be in-place (output in cg instead of cgout, and in cprj instead of cprjout).
!! Otherwise, it is contained in the optional cgout+cprjout pair.

!! In the present status, the cg and cgout relates to all the k points and spins, and rely on the icg index,
!! while it is assumed that cprj and cprjout refer to the specific k point and spin.
!! This is not coherent.
!! THIS MIGHT BE CHANGED IN THE FUTURE !

!! This implementation is NOT band-parallelized
!! Also, it is far of being optimal at the level of linear algebra, and involves extra copying
!! that are detrimental for performance...
!!
!! COPYRIGHT
!! Copyright (C) 2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  alpha_mn(2,nband_out,nband_in)=complex matrix of coefficients of the linear combinations to be computed
!!  dimcprj(natom)=number of lmn components in the <p_{lmn}^i|\psi> for the i-th atom
!!  icg=shift in cg array to locate current k-point and spinpol (for input, and possibly for in-place output)
!!  inplace= if 0, output in cgout and cprjout ; if 1, output in cg and cprj
!!  mcg=second dimension of cg array (mpw*nspinor*mband*mkmem*nsppol)
!!  mcprj=second dimension of cprj array 
!!  natom=number of atoms
!!  nband_in=number of bands, size of the input set of wavefunctions
!!  nband_out=number of bands, size of the output set of wavefunctions (should be equal to nband_in if inplace==1)
!!  npw=number of planewaves in basis at this k point
!!  nspinor=number of spinor components
!!  usepaw=1 if PAW is activated
!!  [icgout= shift in cgout array to locate current k-point and spinpol (for output)]
!!  [mcgout=second dimension of cgout array (mpw*nspinor*mband*mkmem*nsppol)]
!!  [mcprjout=second dimension of cprjout array] 
!!
!! OUTPUT
!!  [cgout(2,mcgout)= plane wave wavefunction coefficients for the set of output wavefunctions]
!!  [cprjout(natom,mcprjout) <type(pawcprj_type)>= projected output wave functions <Proj_i|Cnk> with NL projectors]
!!
!! SIDE EFFECTS
!!  (this quantities are input, and possibly updated output when inplace==1)
!!  cg(2,mcg)= plane wave wavefunction coefficients for the set of input wavefunctions (all k points and spinpol)
!!  cprj(natom,mcprj) <type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk> with NL projectors 
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

 subroutine lincom_cgcprj(alpha_mn,cg,cprj,dimcprj,&
& icg,inplace,mcg,mcprj,natom,nband_in,nband_out,npw,nspinor,usepaw, & 
& cgout,cprjout,icgout) ! optional args

 use defs_basis
 use m_errors 
 use m_profiling_abi
 use m_pawcprj, only : pawcprj_type, pawcprj_alloc, pawcprj_lincom, pawcprj_free

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lincom_cgcprj'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: icg,inplace,mcg,mcprj
 integer, intent(in) :: natom,nband_in,nband_out,npw,nspinor,usepaw
 integer, intent(in),optional :: icgout
!arrays
 integer, intent(in) :: dimcprj(natom)
 real(dp), intent(inout) :: cg(2,mcg)
 real(dp), intent(in) :: alpha_mn(2,nband_out,nband_in)
 real(dp), intent(out),optional :: cgout(:,:)
 type(pawcprj_type),intent(inout) :: cprj(natom,mcprj)
 type(pawcprj_type),intent(out),optional :: cprjout(:,:)

!Local variables-------------------------------
!scalars
 integer :: iband_in,iband_out,ii
!arrays
 real(dp),allocatable :: al(:,:),cgout_(:,:)
 type(pawcprj_type),allocatable :: cprjout_(:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' lincom_cgcprj : enter '
!write(std_out,*)' lincom_cgcprj : npw, nspinor=',npw,nspinor
!ENDDEBUG

 if(inplace==0)then
   if(.not.present(cgout))then
     MSG_ERROR(' inplace==0 while .not.present(cgout) is not permitted ')
   endif
   if(usepaw==1) then
     if(.not.present(cprjout))then 
       MSG_ERROR(' inplace==0 and usepaw==1 while .not.present(cprjout) is not permitted ')
     endif
   endif
 endif

!Take care of the plane wave part
 ABI_ALLOCATE(cgout_,(2,npw*nspinor*nband_out))

 call zgemm('N','N',npw*nspinor,nband_out,nband_in,dcmplx(1._dp), &
&  cg(:,icg+1:icg+npw*nspinor*nband_in),npw*nspinor, &
&  alpha_mn,nband_in,dcmplx(0._dp),cgout_,npw*nspinor)

 if(inplace==1)then
   cg(:,icg+1:icg+npw*nspinor*nband_out)=cgout_
 else
   cgout(:,icgout+1:icgout+npw*nspinor*nband_out)=cgout_
 endif
 ABI_DEALLOCATE(cgout_)

!Take care of the cprj part
 if(usepaw==1) then

   ABI_DATATYPE_ALLOCATE(cprjout_,(natom,nspinor*nband_out))
   call pawcprj_alloc(cprjout_,cprj(1,1)%ncpgr,dimcprj)
   ABI_ALLOCATE(al,(2,nband_in))
   do iband_out=1,nband_out
     ii=(iband_out-1)*nspinor
     do iband_in=1,nband_in
       al(1,iband_in)=alpha_mn(1,iband_out,iband_in)
       al(2,iband_in)=alpha_mn(2,iband_out,iband_in)
     enddo
     call pawcprj_lincom(al,cprj,cprjout_(:,ii+1:ii+nspinor),nband_in)
   enddo
   ABI_DEALLOCATE(al)

   if(inplace==1)then
     cprj=cprjout_
   else
     cprjout=cprjout_
   endif
   call pawcprj_free(cprjout_)
   ABI_DATATYPE_DEALLOCATE(cprjout_)

 endif

end subroutine lincom_cgcprj
!!***
