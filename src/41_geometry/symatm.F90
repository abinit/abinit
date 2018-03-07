!{\src2tex{textfont=tt}}
!!****f* ABINIT/symatm
!! NAME
!! symatm
!!
!! FUNCTION
!! This routine has been improved using ideas of p. 649 of notes,
!! implementing suggestion of Andrew Horsfield: replace search for
!! equivalent atoms using direct primitive cell translations by
!! use of dot product relation which must produce an integer.
!! Relation: $[inv(S(i))*(x(a)-tnons(i)) - x(inv(S)(i,a))] = integer$
!! where $S(i) =$ symmetry matrix in real space, tnons=nonsymmorphic translation
!! (may be 0 0 0), and $x(inv(S)(i,a))$ is sought atom into which $x(a)$ gets
!! rotated by $inv(S)$.  Integer gives primitive translation coordinates to get
!! back to original unit cell.
!! Equivalent to $S*t(b)+tnons-x(a)=another$ $integer$ for $x(b)=x(inv(S))$.
!! For each symmetry operation, find the number of the position to
!! which each atom is sent in the unit cell by the INVERSE of the
!! symmetry operation inv(symrel); i.e. this is the atom which, when acted
!! upon by the given symmetry element isym, gets transformed into atom iatom.
!! This routine uses the fact that inv(symrel)=trans(symrec),
!! the inverse of the symmetry operation expressed in the basis of real
!! space primitive translations equals the transpose of the same symmetry
!! operation expressed in the basis of reciprocal space primitive transl.
!! $xred(nu,indsym(4,isym,ia))=symrec(mu,nu,isym)*(xred(mu,ia)-tnons(mu,isym))
!! - transl(mu)$ where $transl$ is also a set of integers and
!! where translation transl places coordinates within unit cell (note sign).
!! Note that symrec is the set of arrays which are actually input here.
!! These arrays have integer elements.
!! tnons is the nonsymmorphic translation or else is zero.
!! If nsym=1 (i.e. only the identity symmetry is present) then
!! indsym merely takes each atom into itself.
!! The array of integer translations "transl" gets included within array
!! "indsym" as seen below.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! natom=number of atoms in cell.
!! nsym=number of space group symmetries.
!! symrec(3,3,nsym)=symmetries expressed in terms of their action on
!!                  reciprocal space primitive translations (integer).
!! tnons(3,nsym)=nonsymmorphic translations for each symmetry (would
!!               be 0 0 0 each for a symmorphic space group)
!! typat(natom)=integer identifying type of atom.
!! xred(3,natom)=reduced coordinates of atoms in terms of real space
!!               primitive translations
!! tolsym=tolerance for the symmetries
!!
!! OUTPUT
!! indsym(4,nsym,natom)=indirect indexing array described above: for each
!!                      isym,iatom, fourth element is label of atom into
!!                      which iatom is sent by INVERSE of symmetry operation
!!                      isym; first three elements are the primitive
!!                      translations which must be subtracted after the
!!                      transformation to get back to the original unit cell.
!!
!! PARENTS
!!      get_npert_rbz,ingeo,initberry,initorbmag,m_ab7_symmetry,m_crystal,m_ddb
!!      m_polynomial_coeff,m_tdep_sym,setsym,thmeig
!!
!! CHILDREN
!!      symchk,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine symatm(indsym,natom,nsym,symrec,tnons,tolsym,typat,xred)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symatm'
 use interfaces_14_hidewrite
 use interfaces_41_geometry, except_this_one => symatm
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsym
 real(dp), intent(in) :: tolsym
!arrays
 integer,intent(in) :: symrec(3,3,nsym),typat(natom)
 integer,intent(out) :: indsym(4,nsym,natom)
 real(dp),intent(in) :: tnons(3,nsym),xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: eatom,errout,iatom,ii,isym,mu
 real(dp) :: difmax,err
 character(len=500) :: message
!arrays
 integer :: transl(3)
 real(dp) :: difmin(3),tratom(3)

! *************************************************************************

 err=zero
 errout=0

 do isym=1,nsym
   do iatom=1,natom

     do mu=1,3 ! Apply inverse transformation to original coordinates. Note transpose of symrec.
       tratom(mu) = dble(symrec(1,mu,isym))*(xred(1,iatom)-tnons(1,isym))&
&       +dble(symrec(2,mu,isym))*(xred(2,iatom)-tnons(2,isym))&
&       +dble(symrec(3,mu,isym))*(xred(3,iatom)-tnons(3,isym))
     end do
!    
!    Find symmetrically equivalent atom
     call symchk(difmin,eatom,natom,tratom,transl,typat(iatom),typat,xred)
!    
!    Put information into array indsym: translations and label
     indsym(1,isym,iatom)=transl(1)
     indsym(2,isym,iatom)=transl(2)
     indsym(3,isym,iatom)=transl(3)
     indsym(4,isym,iatom)=eatom
!    
!    Keep track of maximum difference between transformed coordinates and
!    nearest "target" coordinate
     difmax=max(abs(difmin(1)),abs(difmin(2)),abs(difmin(3)))
     err=max(err,difmax)

     if (difmax>tolsym) then ! Print warnings if differences exceed tolerance
       write(message, '(3a,i3,a,i6,a,i3,a,a,3es12.4,3a)' )&
&       'Trouble finding symmetrically equivalent atoms',ch10,&
&       'Applying inv of symm number',isym,' to atom number',iatom,'  of typat',typat(iatom),ch10,&
&       'gives tratom=',tratom(1:3),'.',ch10,&
&       'This is further away from every atom in crystal than the allowed tolerance.'
       MSG_WARNING(message)

       write(message, '(a,3i3,a,a,3i3,a,a,3i3)' ) &
&       '  The inverse symmetry matrix is',symrec(1,1:3,isym),ch10,&
&       '                                ',symrec(2,1:3,isym),ch10,&
&       '                                ',symrec(3,1:3,isym)
       call wrtout(std_out,message,'COLL')
       write(message, '(a,3f13.7)' )'  and the nonsymmorphic transl. tnons =',(tnons(mu,isym),mu=1,3)

       call wrtout(std_out,message,'COLL')
       write(message, '(a,1p,3e11.3,a,a,i5)' ) &
&       '  The nearest coordinate differs by',difmin(1:3),ch10,&
&       '  for indsym(nearest atom)=',indsym(4,isym,iatom)
       call wrtout(std_out,message,'COLL')
!      
!      Use errout to reduce volume of error diagnostic output
       if (errout==0) then
         write(message,'(6a)') ch10,&
&         '  This indicates that when symatm attempts to find atoms  symmetrically',ch10, &
&         '  related to a given atom, the nearest candidate is further away than some',ch10,&
&         '  tolerance.  Should check atomic coordinates and symmetry group input data.'
         call wrtout(std_out,message,'COLL')
         errout=1
       end if

     end if !difmax>tol
   end do !iatom
 end do !isym

 if (.FALSE.) then
   do iatom=1,natom
     write(message, '(a,i0,a)' )' symatm: atom number ',iatom,' is reached starting at atom'
     call wrtout(std_out,message,'COLL')
     do ii=1,(nsym-1)/24+1
       if(natom<100)then
         write(message, '(1x,24i3)' ) (indsym(4,isym,iatom),isym=1+(ii-1)*24,min(nsym,ii*24))
       else 
         write(message, '(1x,24i6)' ) (indsym(4,isym,iatom),isym=1+(ii-1)*24,min(nsym,ii*24))
       end if
       call wrtout(std_out,message,'COLL')
     end do
   end do
 end if

 if (err>tolsym) then
   write(message, '(1x,a,1p,e14.5,a,e12.4)' )'symatm: maximum (delta t)=',err,' is larger than tol=',tolsym
   call wrtout(std_out,message,'COLL')
 end if

!Stop execution if error is really big
 if (err>0.01d0) then
   write(message,'(5a)')&
&   'Largest error (above) is so large (0.01) that either input atomic coordinates (xred)',ch10,&
&   'are wrong or space group symmetry data is wrong.',ch10,&
&   'Action : correct your input file.'
   MSG_ERROR(message)
 end if

end subroutine symatm
!!***
