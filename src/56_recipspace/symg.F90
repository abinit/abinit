!{\src2tex{textfont=tt}}
!!****f* ABINIT/symg
!! NAME
!! symg
!!
!! FUNCTION
!! Treat symmetries applied to the G vectors, in view of the application
!! to symmetrization of the dielectric matrix.
!! Generate a list of time-reversed G vectors, as well as a list
!! of spatially-symmetric G vectors.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!! npwdiel=number of planewaves for the dielectric matrix
!! nsym=number of symmetry
!! symrel(3,3,nsym)=symmetry matrices in real space (integers)
!! tnons(3,nsym)=reduced nonsymmorphic translations
!! (symrel and tnons are in terms of real space primitive translations)
!!
!! OUTPUT
!! phdiel(2,npwdiel,nsym)=phase associated with point symmetries applied to G
!! sym_g(npwdiel,nsym)=index list of symmetric G vectors
!! (could save a bit of space by suppressing isym=1, since the
!! corresponding symmetry is the identity)
!! tmrev_g(npwdiel)=index list of inverted G vectors (time-reversed)
!!
!! PARENTS
!!      suscep_stat
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine symg(kg_diel,npwdiel,nsym,phdiel,sym_g,symrel,tmrev_g,tnons)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symg'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwdiel,nsym
!arrays
 integer,intent(in) :: kg_diel(3,npwdiel),symrel(3,3,nsym)
 integer,intent(out) :: sym_g(npwdiel,nsym),tmrev_g(npwdiel)
 real(dp),intent(in) :: tnons(3,nsym)
 real(dp),intent(out) :: phdiel(2,npwdiel,nsym)

!Local variables-------------------------------
!scalars
 integer :: g1,g2,g3,ipw,isym,j1,j2,j3,m1m,m1p,m2m,m2p,m3m,m3p,symmg,trevg
 real(dp) :: arg,tau1,tau2,tau3
 character(len=500) :: message
!arrays
 integer,allocatable :: grid(:,:,:)

! *************************************************************************

!Determines maximal bounds of the zone spanned by the planewaves
 m1m=0 ; m2m=0 ; m3m=0 ; m1p=0 ; m2p=0 ; m3p=0
 do ipw=1,npwdiel
   g1=kg_diel(1,ipw)
   g2=kg_diel(2,ipw)
   g3=kg_diel(3,ipw)
   if(g1<m1m)m1m=g1 ; if(g1>m1p)m1p=g1
   if(g2<m2m)m2m=g2 ; if(g2>m2p)m2p=g2
   if(g3<m3m)m3m=g3 ; if(g3>m3p)m3p=g3
 end do

!Set up grid, that associate to each point the index of the
!corresponding planewave, if there is one
 ABI_ALLOCATE(grid,(m1m:m1p,m2m:m2p,m3m:m3p))
 grid(:,:,:)=0
 do ipw=1,npwdiel
   g1=kg_diel(1,ipw)
   g2=kg_diel(2,ipw)
   g3=kg_diel(3,ipw)
   grid(g1,g2,g3)=ipw
 end do

!Set up tmrev_g and sym_g arrays
 do ipw=1,npwdiel
   g1=kg_diel(1,ipw)
   g2=kg_diel(2,ipw)
   g3=kg_diel(3,ipw)

!  Treat first time-reversal symmetry
   trevg=grid(-g1,-g2,-g3)
   if(trevg==0)then
     message = ' Do not find the time-reversed symmetric of a G-vector.'
     MSG_BUG(message)
   end if
   tmrev_g(ipw)=trevg

!  Treat now spatial symmetries
   do isym=1,nsym

!    Get rotated G vector Gj for each symmetry element
!    -- here we use the TRANSPOSE of symrel; assuming symrel expresses
!    the rotation in real space, the transpose is then appropriate
!    for G space symmetrization (according to Doug : see routine irrzg.f)
     j1=symrel(1,1,isym)*g1+&
&     symrel(2,1,isym)*g2+symrel(3,1,isym)*g3
     j2=symrel(1,2,isym)*g1+&
&     symrel(2,2,isym)*g2+symrel(3,2,isym)*g3
     j3=symrel(1,3,isym)*g1+&
&     symrel(2,3,isym)*g2+symrel(3,3,isym)*g3
     symmg=grid(j1,j2,j3)
     if(symmg==0)then
       message = ' Do not find the spatially symmetric of a G-vector.'
       MSG_BUG(message)
     end if
     sym_g(ipw,isym)=symmg

!    Get associated phase
     tau1=tnons(1,isym)
     tau2=tnons(2,isym)
     tau3=tnons(3,isym)
     if (abs(tau1)>tol12.or.abs(tau2)>tol12.or.abs(tau3)>tol12) then
!      compute exp(-2*Pi*I*G dot tau) using original G
       arg=two_pi*(dble(g1)*tau1+dble(g2)*tau2+dble(g3)*tau3)
       phdiel(1,ipw,isym)=cos(arg)
       phdiel(2,ipw,isym)=-sin(arg)
     else
       phdiel(1,ipw,isym)=1._dp
       phdiel(2,ipw,isym)=0._dp
     end if

   end do

 end do

 ABI_DEALLOCATE(grid)

end subroutine symg
!!***
