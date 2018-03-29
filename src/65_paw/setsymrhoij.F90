!{\src2tex{textfont=tt}}
!!****f* ABINIT/setsymrhoij
!! NAME
!! setsymrhoij
!!
!! FUNCTION
!! PAW only
!! Compute rotation matrices expressed in the basis of real spherical harmonics
!! This coefficients are used later to symmetrize rhoij quantities (augmentation occupancies)
!! and other similar quantities.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (NH, FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gprimd(3,3)==dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  lmax=value of lmax mentioned at the second line of the psp file
!!  nsym=number of symmetry elements in space group
!!  pawprtvol=control print volume and debugging output for PAW
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  sym(3,3,nsym)=symmetries of group in terms of operations on primitive translations
!!
!! OUTPUT
!!  zarot(2*lmax+1,2*lmax+1,lmax+1,nsym)=coefficients of the
!!      transformation of real spherical harmonics
!!      under the symmetry operations
!!
!! NOTES
!!  Typical use: sym(:,:,:) is symrec(:,:,:) (rotations in reciprocal space)
!!               because we need symrel^-1 (=transpose[symrec])
!!               to symmetrize quantities.
!!
!!  - This file comes from the file crystal_symmetry.f
!!    by N.A.W. Holzwarth and A. Tackett for the code pwpaw
!!  - Uses sign & phase convension of  M. E. Rose, Elementary Theory of Angular
!!    Momentum, John Wiley & Sons,. inc. 1957)
!!    zalpha = exp(-i*alpha)   zgamma = exp (-i*gamma)
!!  - Assumes each transformation  can be expressed in terms of 3 Euler
!!    angles with or without inversion
!!
!!  Reference for evaluation of rotation matrices in the basis of real SH:
!!  Blanco M.A., Florez M. and Bermejo M.
!!  Journal of Molecular Structure: THEOCHEM, Volume 419, Number 1, 8 December 1997 , pp. 19-27(9)
!!  http://www.unioviedo.es/qcg/art/Theochem419-19-ov-BF97-rotation-matrices.pdf
!!
!! PARENTS
!!      bethe_salpeter,dfpt_looppert,gstate,initberry,initorbmag,respfn
!!      screening,sigma,wfk_analyze
!!
!! CHILDREN
!!      mkeuler,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine setsymrhoij(gprimd,lmax,nsym,pawprtvol,rprimd,sym,zarot)

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_special_funcs, only : phim
 use m_angles,        only : mkeuler, dbeta

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setsymrhoij'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lmax,nsym,pawprtvol
!arrays
 integer,intent(in) :: sym(3,3,nsym)
 real(dp),intent(in) :: gprimd(3,3),rprimd(3,3)
 real(dp),intent(out) :: zarot(2*lmax+1,2*lmax+1,lmax+1,nsym)

!Local variables ------------------------------
!scalars
 integer :: i1,ii,il,irot,isn,j1,jj,k1,ll,mm,mp
 real(dp) :: cosalp,cosbeta,cosgam,sinalp,singam
 character(len=1000) :: message
!arrays
 real(dp) :: prod(3,3),rot(3,3)
!************************************************************************

 DBG_ENTER("COLL")

 if (abs(pawprtvol)>=3) then
   write(message,'(8a,i4)') ch10,&
&   ' PAW TEST:',ch10,&
&   ' ==== setsymrhoij: rotation matrices in the basis ============',ch10,&
&   ' ====              of real spherical harmonics    ============',ch10,&
&   '  > Number of symmetries (nsym)=',nsym
   call wrtout(std_out,message,'COLL')
 end if

 zarot=zero

 do irot=1,nsym

   if (abs(pawprtvol)>=3) then
     write(message,'(a,i2,a,9i2,a)') '   >For symmetry ',irot,' (',sym(:,:,irot),')'
     call wrtout(std_out,message,'COLL')
   end if

!  === l=0 case ===
   zarot(1,1,1,irot)=one

!  === l>0 case ===
   if (lmax>0) then
!    Calculate the rotations in the cartesian basis
     rot=zero;prod=zero
     do k1=1,3
       do j1=1,3
         do i1=1,3
           prod(i1,j1)=prod(i1,j1)+sym(i1,k1,irot)*rprimd(j1,k1)
         end do
       end do
     end do
     do j1=1,3
       do i1=1,3
         do k1=1,3
           rot(i1,j1)=rot(i1,j1)+gprimd(i1,k1)*prod(k1,j1)
         end do
         if(abs(rot(i1,j1))<tol10) rot(i1,j1)=zero
       end do
     end do
     call mkeuler(rot,cosbeta,cosalp,sinalp,cosgam,singam,isn)
     do ll=1,lmax
       il=(isn)**ll
       do mp=-ll,ll
         jj=mp+ll+1
         do mm=-ll,ll
           ii=mm+ll+1

!          Formula (47) from the paper of Blanco et al
           zarot(ii,jj,ll+1,irot)=il&
&           *(phim(cosalp,sinalp,mm)*phim(cosgam,singam,mp)*sign(1,mp)&
           *(dbeta(cosbeta,ll,abs(mp),abs(mm))&
&           +(-1._dp)**mm*dbeta(cosbeta,ll,abs(mm),-abs(mp)))*half&
&           -phim(cosalp,sinalp,-mm)*phim(cosgam,singam,-mp)*sign(1,mm)&
           *(dbeta(cosbeta,ll,abs(mp),abs(mm))&
&           -(-1._dp)**mm*dbeta(cosbeta,ll,abs(mm),-abs(mp)))*half)
         end do
       end do
     end do
   end if   ! lmax case

   if (abs(pawprtvol)>=3) then
     if(lmax>0) then
       write(message,'(2a,3(3(2x,f7.3),a))') &
&       '    Rotation matrice for l=1:',ch10,&
&       (zarot(1,jj,2,irot),jj=1,3),ch10,&
&       (zarot(2,jj,2,irot),jj=1,3),ch10,&
&       (zarot(3,jj,2,irot),jj=1,3)
       call wrtout(std_out,message,'COLL')
     end if
     if(lmax>1) then
       write(message,'(2a,5(5(2x,f7.3),a))') &
&       '    Rotation matrice for l=2:',ch10,&
&       (zarot(1,jj,3,irot),jj=1,5),ch10,&
&       (zarot(2,jj,3,irot),jj=1,5),ch10,&
&       (zarot(3,jj,3,irot),jj=1,5),ch10,&
&       (zarot(4,jj,3,irot),jj=1,5),ch10,&
&       (zarot(5,jj,3,irot),jj=1,5)
       call wrtout(std_out,message,'COLL')
     end if
     if(lmax>2) then
       write(message,'(2a,7(7(2x,f7.3),a))') &
&       '    Rotation matrice for l=3:',ch10,&
&       (zarot(1,jj,4,irot),jj=1,7),ch10,&
&       (zarot(2,jj,4,irot),jj=1,7),ch10,&
&       (zarot(3,jj,4,irot),jj=1,7),ch10,&
&       (zarot(4,jj,4,irot),jj=1,7),ch10,&
&       (zarot(5,jj,4,irot),jj=1,7),ch10,&
&       (zarot(6,jj,4,irot),jj=1,7),ch10,&
&       (zarot(7,jj,4,irot),jj=1,7)
       call wrtout(std_out,message,'COLL')
     end if
   end if

 end do  ! isym loop

 DBG_EXIT("COLL")

end subroutine setsymrhoij
!!***
