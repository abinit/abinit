!!****m* ABINIT/m_paw_lmn
!! NAME
!!  m_paw_lmn
!!
!! FUNCTION
!!  This module provides tools to calculate tables commonly used to iterate
!!  over the the (l,m,n) channels of the PAW partial waves.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_lmn

 use defs_basis
 use m_abicore
 use m_errors

 implicit none

 private

! Public procedures.
 public :: ilm2lm          ! Returns the value of l and m in from the index ilm
 public :: make_indlmn     ! Calculates the indlmn(6,lmn_size) table giving l,m,n,lm,ln,spin for i=lmn.
 public :: make_indklmn    ! Calculates the indklmn(8,lmn2_size) table giving
                           !   klm, kln, abs(il-jl), (il+jl), ilm and jlm, ilmn and jlmn for each symmetric klmn=(ilmn,jlmn)
 public :: make_kln2ln     ! Calculates the kln2ln(6,ln2_size) table giving
                           !   il, jl ,in, jn, iln, jln for each symmetric kln=(iln,jln)
 public :: make_klm2lm     ! Calculates the klm2lm(6,lm2_size) table giving
                           !   il, jl ,im, jm, ilm, jlm for each symmetric klm=(ilm,jlm)
 public :: klmn2ijlmn      ! Calculates ilmn and jlmn from klmn.
 public :: make_indln      ! Calculates indln(2,ln_size) giving l and n for i=ln
 public :: uppert_index    ! The sequential index of an element in the upper triangle of a matrix 

CONTAINS  !========================================================================================
!!***

!!****f* m_paw_lmn/ilm2lm
!! NAME
!!  ilm2lm
!!
!! FUNCTION
!!  Returns the value of l and m in from the index ilm
!!
!! INPUTS
!!  ilm=The contracted index for (l,m).
!!
!! OUTPUT
!!  ll=The angular momentun defined in [0,1,2,3,4].
!!  mm=The magnetic number in the interval [-l,....+l].
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ilm2lm(ilm,ll,mm)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ilm
 integer,intent(out) :: ll,mm

!Local variables-------------------------------
 integer :: ii
! *********************************************************************

 if (ilm<1) then 
   ABI_ERROR("Wrong ilm")
 end if

 ll = -1
 do ii=0,100
  if ( (ii+1)**2 >= ilm) then
   ll = ii
   EXIT
  end if
 end do

 mm = ilm - ll**2 -ll-1

 if (ll==-1) then
   ABI_ERROR("l>100 not programmed!")
 end if

end subroutine ilm2lm
!!***

!----------------------------------------------------------------------

!!****f* m_paw_lmn/make_indlmn
!! NAME
!!  make_indlmn
!!
!! FUNCTION
!!  Performs the setup of the indlmn table for PAW calculations (indices for (l,m,n) basis)
!!
!! INPUTS
!!  ln_size= Total number of nl components
!!  lmn_size= Second dimension in indlmn. Total number of (l,m,n) components for this atom.
!!  orbitals(ln_size)=Give the value of l for each element of the augmented basis set.
!!
!! OUTPUT
!!  indlmn(6,lmn_size)=array giving l,m,n,lm,ln,s for i=lmn
!!
!! PARENTS
!!      m_paw_atomorb
!!
!! CHILDREN
!!
!! SOURCE

subroutine make_indlmn(ln_size,lmn_size,orbitals,indlmn)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ln_size,lmn_size
 integer,intent(in) :: orbitals(ln_size)
!scalars
 integer,intent(out) :: indlmn(6,lmn_size)

!Local variables ------------------------------
!scalars
 integer :: ilmn,ib,il,iln,ilm
!arrays
 integer,allocatable :: nprj(:)

!************************************************************************

 ABI_MALLOC(nprj,(0:MAXVAL(orbitals)))

 ilmn=0; iln=0; nprj=0
 do ib=1,ln_size
   il=orbitals(ib)
   nprj(il)=nprj(il)+1
   iln=iln+1
   do ilm=1,2*il+1
     indlmn(1,ilmn+ilm)=il           ! l
     indlmn(2,ilmn+ilm)=ilm-(il+1)   ! m
     indlmn(3,ilmn+ilm)=nprj(il)     ! n
     indlmn(4,ilmn+ilm)=il*il+ilm    ! lm index
     indlmn(5,ilmn+ilm)=iln          ! ln index
     indlmn(6,ilmn+ilm)=1            ! spin (not yet used!)
   end do
   ilmn=ilmn+2*il+1
 end do

 ABI_FREE(nprj)

end subroutine make_indlmn
!!***

!----------------------------------------------------------------------

!!****f* m_paw_lmn/make_indklmn
!! NAME
!!  make_indklmn
!!
!! FUNCTION
!!  Performs the setup of the indklmn table for PAW calculations.
!!  Compute the indklmn indexes giving klm, kln, abs(il-jl) and (il+jl), ilm and jlm, ilmn and jlmn
!!  for each klmn=(ilmn,jlmn) with jlmn >= ilmn
!!
!! INPUTS
!!  lcutdens=Maximum l for densities/potentials moments computations
!!  lmn_size=Number of (l,m,n) elements for the PAW basis set.
!!  lmn2_size=Number of elements in the symmetric basis set: lmn2_size=lmn_size*(lmn_size+1)/2
!!  indlmn(6,lmn_size)=Array giving l,m,n,lm,ln,spin for i=lmn.
!!
!! OUTPUT
!!  indklmn(6,lmn2_size)=Array giving klm, kln, abs(il-jl), (il+jl), ilm and jlm, ilmn and jlmn
!!    for each klmn=(ilmn,jlmn). Note: ilmn=(il,im,in) and ilmn<=jlmn
!!  klm_diag(lmn2_size)=1 il==jl and im==jm, 0 otherwise.
!!
!! PARENTS
!!      m_paw_atomorb
!!
!! CHILDREN
!!
!! SOURCE

subroutine make_indklmn(lcutdens,lmn_size,lmn2_size,indlmn,indklmn,klm_diag)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmn2_size,lmn_size
 integer,intent(in) ::  lcutdens
!scalars
 integer,intent(in) :: indlmn(6,lmn_size)
 integer,intent(out) :: indklmn(8,lmn2_size)
 integer,intent(out) :: klm_diag(lmn2_size)

!Local variables ------------------------------
!scalars
 integer :: i0lm,i0ln,il,ilm,ilmn,iln
 integer :: j0lm,j0lmn,j0ln,jl,jlm,jlmn,jln,klmn

!************************************************************************

 klm_diag=0

 do jlmn=1,lmn_size
  jl= indlmn(1,jlmn); jlm=indlmn(4,jlmn); jln=indlmn(5,jlmn)
  j0lmn=jlmn*(jlmn-1)/2
  j0lm =jlm *(jlm -1)/2
  j0ln =jln *(jln -1)/2
  do ilmn=1,jlmn
   il=indlmn(1,ilmn); ilm=indlmn(4,ilmn); iln=indlmn(5,ilmn)
   klmn=j0lmn+ilmn
   if (ilm<=jlm) then
    indklmn(1,klmn)=j0lm+ilm     !klm
   else
    i0lm=ilm*(ilm-1)/2
    indklmn(1,klmn)=i0lm+jlm
   end if
   if (iln<=jln) then
    indklmn(2,klmn)=j0ln+iln     !kln
   else
    i0ln=iln*(iln-1)/2
    indklmn(2,klmn)=i0ln+jln
   end if
   !MG This is not safe, what happens if lcutdens < |il-jl|?
   indklmn(3,klmn)=MIN(ABS(il-jl),lcutdens)  ! abs(li-lj) NB this is a l-value, not an index >=1.
   indklmn(4,klmn)=MIN(il+jl,lcutdens)       ! abs(li+lj) NB this is a l-value, not an index >=1.

   indklmn(5,klmn)=ilm                       ! ilm
   indklmn(6,klmn)=jlm                       ! jlm

   indklmn(7,klmn)=ilmn                      ! ilmn
   indklmn(8,klmn)=jlmn                      ! jlmn


   if (ilm==jlm) klm_diag(klmn)=1
  end do
 end do

end subroutine make_indklmn
!!***

!----------------------------------------------------------------------

!!****f* m_paw_lmn/make_kln2ln
!! NAME
!!  make_kln2ln
!!
!! FUNCTION
!!  Performs the setup of the kln2ln table for PAW calculations.
!!
!! INPUTS
!!  lmn_size=Number of (l,m,n) elements for the PAW basis set.
!!  lmn2_size=Number of elements in the symmetric basis set: lmn2_size=lmn_size*(lmn_size+1)/2
!!  ln2_size=Number of symmetric (l,n) channels i.e. ln_size*(ln_size+1)/2
!!  indlmn(6,lmn_size)=Array giving l,m,n,lm,ln,spin for i=lmn.
!!  indklmn(8,lmn2_size)=Array giving klm, kln, abs(il-jl), (il+jl), ilm and jlm, ilmn and jlmn
!!   for each klmn=(ilmn,jlmn). Note: ilmn=(il,im,in) and ilmn<=jlmn
!!
!! OUTPUT
!!  kln2ln(6,ln2_size)=Table giving il, jl ,in, jn, iln, jln for each kln=(iln,jln)
!!  where iln=(il,in) and iln<=jln. NB: kln2ln is an application and not a bijection
!!
!! PARENTS
!!      m_paw_atomorb,m_paw_slater
!!
!! CHILDREN
!!
!! SOURCE

subroutine make_kln2ln(lmn_size,lmn2_size,ln2_size,indlmn,indklmn,kln2ln)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmn2_size,lmn_size,ln2_size
!arrays
 integer,intent(in) :: indlmn(6,lmn_size)
 integer,intent(in) :: indklmn(8,lmn2_size)
 integer,intent(out) :: kln2ln(6,ln2_size)

!Local variables ------------------------------
!scalars
 integer :: il,in,ilmn,iln
 integer :: jl,jn,jlmn,jln
 integer :: klmn,kln,kln_old

!************************************************************************

 kln2ln = 0

 kln_old = -1
 do klmn=1,lmn2_size
   kln = indklmn(2,klmn)
   if (kln /= kln_old) then
     kln_old = kln

     call klmn2ijlmn(klmn,lmn_size,ilmn,jlmn)

     il= indlmn(1,ilmn); iln=indlmn(5,ilmn)
     jl= indlmn(1,jlmn); jln=indlmn(5,jlmn)
     !$in = indklmn(9 ,klmn) !=indlmn(3,ilmn)
     !$jn = indklmn(10,klmn) !=indlmn(3,jlmn)

     in = indlmn(3,ilmn)
     jn = indlmn(3,jlmn)

     ! Shift il to have the index instead of the value of l.
     kln2ln(1,kln) = il +1
     kln2ln(2,kln) = jl +1
     kln2ln(3,kln) = in
     kln2ln(4,kln) = jn
     kln2ln(5,kln) = iln
     kln2ln(6,kln) = jln
   end if
 end do !klmn

end subroutine make_kln2ln
!!***

!----------------------------------------------------------------------

!!****f* m_paw_lmn/make_klm2lm
!! NAME
!!  make_klm2lm
!!
!! FUNCTION
!!  Performs the setup of the klm2lm table for PAW calculations.
!!
!! INPUTS
!!  lmn_size=Number of (l,m,n) elements for the PAW basis set.
!!  lmn2_size=Number of elements in the symmetric basis set: lmn2_size=lmn_size*(lmn_size+1)/2
!!  lm2_size)=Number of (l.m) elements in the symmetric basis set.
!!  indlmn(6,lmn_size)=Array giving l,m,n,lm,ln,spin for i=lmn.
!!  indklmn(8,lmn2_size)=Array giving klm, kln, abs(il-jl), (il+jl), ilm and jlm
!!   for each klmn=(ilmn,jlmn). Note: ilmn=(il,im,in) and ilmn<=jlmn
!!
!! OUTPUT
!!  klm2lm(6,lm2_size)=Table giving il, jl ,im, jm, ilm, jlm for each klm=(ilm,jlm)
!!  where ilm=(il,im) and ilm<=jlm. NB: klm2lm is an application and not a bijection.
!!
!! NOTES
!!  klm2lm can be calculated easily if we assume that all (l,m) channels
!!  are ordered by increasing l and m. This is the standard convention
!!  used in most of the PAW datasets. This routines, howevever, works
!!  works also in the unlikely case in with (l,m) are not ordered.
!!
!! PARENTS
!!      m_paw_slater
!!
!! CHILDREN
!!
!! SOURCE

subroutine make_klm2lm(lmn_size,lmn2_size,lm2_size,indlmn,indklmn,klm2lm)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmn2_size,lmn_size,lm2_size
!arrays
 integer,intent(in) :: indlmn(6,lmn_size)
 integer,intent(in) :: indklmn(8,lmn2_size)
 integer,intent(out) :: klm2lm(6,lm2_size)

!Local variables ------------------------------
!scalars
 integer :: il,ilm,ilmn,jl,jlm,jlmn,im,jm
 integer :: klmn,klm,klm_old !,iklm

!************************************************************************

 klm2lm = 0
 klm_old = -1
 do klmn=1,lmn2_size
   klm = indklmn(1,klmn)
   if (klm /= klm_old) then
     klm_old = klm

     ilm = indklmn(5,klmn)
     jlm = indklmn(6,klmn)

     call klmn2ijlmn(klmn,lmn_size,ilmn,jlmn)

     il=indlmn(1,ilmn); im=indlmn(2,ilmn)
     jl=indlmn(1,jlmn); jm=indlmn(2,jlmn)

     !shift to have the index instead of l or m.
     klm2lm(1,klm) = il +1       ! il
     klm2lm(2,klm) = jl +1       ! jl
     klm2lm(3,klm) = im + il +1  ! im
     klm2lm(4,klm) = jm + jl +1  ! jm
     klm2lm(5,klm) = ilm         ! ilm
     klm2lm(6,klm) = jlm         ! jlm
   end if
 end do !klmn

 if (ANY(klm2lm==0)) then
   ABI_BUG("check klm2lm")
 end if

!DEBUG
#if 0
 write(std_out,*)"Debugging make_klm2lm:"
 do iklm=1,lm2_size
  il  = klm2lm(1,iklm) -1      ! li
  jl  = klm2lm(2,iklm) -1      ! lj
  im  = klm2lm(3,iklm) -il -1  ! mi
  jm  = klm2lm(4,iklm) -jl -1  ! mj
  ilm = klm2lm(5,iklm)         ! ilm
  jlm = klm2lm(6,iklm)         ! jlm
  write(std_out,'(i3,2(a,2i3,a))')"iklm ",iklm," l m (",il,im,")"," l m (",jl,jm,")"
 end do
#endif

end subroutine make_klm2lm
!!***

!----------------------------------------------------------------------

!!****f* m_paw_lmn/klmn2ijlmn
!! NAME
!!  klmn2ijlmn
!!
!! FUNCTION
!!  Find ilmn and jlmn from klmn and lmn_size.
!!
!! INPUTS
!!  lmn_size=Number of (l,m,n) elements for the PAW basis set.
!!  klmn=The index corresponding to (ilmn,jlmn) in packed form.
!!
!! OUTPUT
!!  jlmn, ilmn=The two symmetrix indices corresponding to klmn. NB: jlmn >= ilmn
!!
!! PARENTS
!!      m_paw_lmn,m_paw_slater
!!
!! CHILDREN
!!
!! SOURCE

subroutine klmn2ijlmn(klmn,lmn_size,ilmn,jlmn)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: klmn,lmn_size
 integer,intent(out) :: ilmn,jlmn

!Local variables ------------------------------
!scalars
 integer :: ii,jj,k0
!************************************************************************

 ilmn=-1; jlmn=-1

 do jj=1,lmn_size
   k0=jj*(jj-1)/2
   do ii=1,jj
     if (klmn==ii+k0) then
       ilmn = ii
       jlmn = jj; RETURN
     end if
   end do
 end do

 ABI_BUG("Not able to found ilmn and jlmn")

end subroutine klmn2ijlmn
!!***

!----------------------------------------------------------------------

!!****f* m_paw_lmn/make_indln
!! NAME
!!  make_indln
!!
!! FUNCTION
!!  Performs the setup of the indln table for PAW calculations.
!!  Compute the indln indexes giving ilmn=(l,m,n)
!!
!! INPUTS
!!  lmn_size=Number of (l,m,n) elements for the PAW basis set.
!!  ln_size=Number of (l,n) elements
!!  indlmn(6,lmn_size)=Array giving l,m,n,lm,ln,spin for i=lmn.
!!
!! OUTPUT
!!  indln(2,ln_size)=Array giving l and n for i=ln
!!
!! PARENTS
!!      m_paw_slater
!!
!! CHILDREN
!!
!! SOURCE

subroutine make_indln(lmn_size,ln_size,indlmn,indln)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmn_size,ln_size
!arrays
 integer,intent(in) :: indlmn(6,lmn_size)
 integer,intent(out) :: indln(2,ln_size)

!Local variables ------------------------------
!scalars
 integer :: ilmn,ll,nn,ll_old,nn_old,ii
!************************************************************************

 ii=0; ll_old=-1; nn_old=-1
 do ilmn=1,lmn_size
  ll=indlmn(1,ilmn)
  nn=indlmn(3,ilmn)
  if (ll/=ll_old.or.nn/=nn_old) then
   ll_old=ll
   nn_old=nn
   ii=ii+1
   indln(1,ii)=ll
   indln(2,ii)=nn
  end if
 end do

 ABI_CHECK(ii==ln_size,"ii/=ln_size")

end subroutine make_indln
!!***

!----------------------------------------------------------------------

!!****f* m_paw_lmn/uppert_index
!! NAME
!!   uppert_index
!!
!! FUNCTION
!!  Helper function returning the sequential index of an element in the upper triangle of a matrix 
!!  given the row-index ii and the column-index jj. If ii>jj the index of the element a_{jj,ii} is returned.
!!
!! INPUTS
!!  ii=Row index
!!  jj=column index
!!
!! PARENTS
!!
!! SOURCE

function uppert_index(ii,jj)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ii,jj
 integer :: uppert_index
!scalars

!************************************************************************

 if (jj>=jj) then
   uppert_index = ii + jj*(jj-1)/2
 else
   uppert_index = jj + ii*(ii-1)/2
 end if

end function uppert_index
!!***

END MODULE m_paw_lmn
!!***
