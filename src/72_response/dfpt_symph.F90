!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_symph
!! NAME
!! dfpt_symph
!!
!! FUNCTION
!! Determine the symmetry character of the different phonon modes.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2017 ABINIT group (GMR,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3)=length scales of primitive translations (bohr)
!! eigvec(2*3*natom*3*natom)=eigenvectors of the dynamical matrix
!! indsym(4,nsym,natom)=indirect indexing array : for each
!!   isym,iatom, fourth element is label of atom into which iatom is sent by
!!   INVERSE of symmetry operation isym; first three elements are the primitive
!!   translations which must be subtracted after the transformation to get back
!!   to the original unit cell.
!! iout=unit number to which output is written
!! natom=number of atoms in unit cell
!! nsym=number of space group symmetries
!! phfrq(3*natom)=phonon frequencies (Hartree units)
!! rprim(3,3)=dimensionless primitive translations in real space
!! symrel(3,3,nsym)=matrices of the group symmetries (real space)
!!
!! OUTPUT
!!
!! PARENTS
!!      anaddb,m_phonons
!!
!! CHILDREN
!!      matr3inv,mkrdim,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_symph(iout,acell,eigvec,indsym,natom,nsym,phfrq,rprim,symrel)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_symph'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,natom,nsym
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrel(3,3,nsym)
 real(dp),intent(in) :: acell(3),eigvec(2*3*natom*3*natom),phfrq(3*natom)
 real(dp),intent(in) :: rprim(3,3)

!Local variables -------------------------
!scalars
 integer :: iad1,iad2,iad3,iatom,idir,ii1,ii2,ii3,imode,isym,itol,jad,jatom,jj
 integer :: jmode,kk,ntol
 character(len=500) :: message
!arrays
 integer,allocatable :: degeneracy(:),integer_characters(:),symind(:,:)
 real(dp) :: gprimd(3,3),rprimd(3,3)
 real(dp),allocatable :: eigvtr(:),redvec(:),redvtr(:),symph(:,:)

!******************************************************************

!Compute dimensional primitive translations rprimd and its inverse gprimd
 call mkrdim(acell,rprim,rprimd)
 call matr3inv(rprimd,gprimd)

!Build the symmetry index (inverse of indsym(4,:,:))
 ABI_ALLOCATE(symind,(nsym,natom))
 do isym=1,nsym
   do iatom=1,natom
     symind(isym,indsym(4,isym,iatom))=iatom
   end do
 end do

 ABI_ALLOCATE(symph,(nsym,3*natom))
 ABI_ALLOCATE(redvec,(2*3*natom))
 ABI_ALLOCATE(redvtr,(2*3*natom))
 ABI_ALLOCATE(eigvtr,(2*3*natom))

!Loop over the vibration modes
 do imode=1,3*natom

!  Compute eigvec for this mode in reduced coordinates redvec
   do iatom=1,natom
     iad1=3*(iatom-1)+1
     ii1=2*3*natom*(imode-1)+2*(iad1-1)+1
     iad2=3*(iatom-1)+2
     ii2=2*3*natom*(imode-1)+2*(iad2-1)+1
     iad3=3*(iatom-1)+3
     ii3=2*3*natom*(imode-1)+2*(iad3-1)+1
     do idir=1,3
       jad=3*(iatom-1)+idir
       jj=2*(jad-1)+1
       redvec(jj)=gprimd(1,idir)*eigvec(ii1)+&
&       gprimd(2,idir)*eigvec(ii2)+&
&       gprimd(3,idir)*eigvec(ii3)
       redvec(jj+1)=gprimd(1,idir)*eigvec(ii1+1)+&
&       gprimd(2,idir)*eigvec(ii2+1)+&
&       gprimd(3,idir)*eigvec(ii3+1)
     end do !idir
   end do !iatom

!  Apply each transformation to redvec and store at the correct location in redvtr (iatom -> jatom)
   do isym=1,nsym
     do iatom=1,natom
       jatom=symind(isym,iatom)
       iad1=3*(iatom-1)+1
       ii1=2*(iad1-1)+1
       iad2=3*(iatom-1)+2
       ii2=2*(iad2-1)+1
       iad3=3*(iatom-1)+3
       ii3=2*(iad3-1)+1
       do idir=1,3
         jad=3*(jatom-1)+idir
         jj=2*(jad-1)+1
         redvtr(jj)=dble(symrel(idir,1,isym))*redvec(ii1)+&
&         dble(symrel(idir,2,isym))*redvec(ii2)+&
&         dble(symrel(idir,3,isym))*redvec(ii3)
         redvtr(jj+1)=dble(symrel(idir,1,isym))*redvec(ii1+1)+&
&         dble(symrel(idir,2,isym))*redvec(ii2+1)+&
&         dble(symrel(idir,3,isym))*redvec(ii3+1)

       end do !idir
     end do !iatom

!    Compute redvtr in cartesian coordinates eigvtr
     do iatom=1,natom
       iad1=3*(iatom-1)+1
       ii1=2*(iad1-1)+1
       iad2=3*(iatom-1)+2
       ii2=2*(iad2-1)+1
       iad3=3*(iatom-1)+3
       ii3=2*(iad3-1)+1
       do idir=1,3
         jad=3*(iatom-1)+idir
         jj=2*(jad-1)+1
         eigvtr(jj)=rprimd(idir,1)*redvtr(ii1)+&
&         rprimd(idir,2)*redvtr(ii2)+&
&         rprimd(idir,3)*redvtr(ii3)
         eigvtr(jj+1)=rprimd(idir,1)*redvtr(ii1+1)+&
&         rprimd(idir,2)*redvtr(ii2+1)+&
&         rprimd(idir,3)*redvtr(ii3+1)
       end do !idir
     end do !iatom

!    Compute scalar product...
     symph(isym,imode)=0.0_dp
     do jad=1,3*natom
       jj=2*(jad-1)+1
       kk=2*3*natom*(imode-1)+2*(jad-1)+1
       symph(isym,imode)=symph(isym,imode)+eigvtr(jj)*eigvec(kk)+eigvtr(jj+1)*eigvec(kk+1)
     end do

   end do !isym
 end do !imode

!Treat degeneracies (different tolerances will be tried)
!Compute the order of the degeneracy, and
!attribute it to the lowest of the degenerate modes
!Also attribute the characters to the lowest mode
!When all the characters are integers, consider that the
!mode is non-degenerate. The maximum difference in frequency
!that is tolerated is on the order of 4cm-1 (which is large...)
 ABI_ALLOCATE(degeneracy,(3*natom))
 ABI_ALLOCATE(integer_characters,(3*natom))
 degeneracy(:)=1
 integer_characters(:)=0
 do itol=1,20
   ntol=itol
   do imode=3*natom,2,-1
     if(integer_characters(imode)==0)then
       do jmode=imode-1,1,-1
         if(integer_characters(jmode)==0)then
           if(abs(phfrq(imode)-phfrq(jmode))<itol*tol6)then
             degeneracy(jmode)=degeneracy(jmode)+degeneracy(imode)
             degeneracy(imode)=0
             symph(:,jmode)=symph(:,jmode)+symph(:,imode)
             symph(:,imode)=0.0_dp
           end if
         end if !integer_characters(jmode)==0
       end do !jmode
     end if !integer_characters(imode)==0
   end do !imode
   do imode=1,3*natom
     if(maxval(abs( symph(:,imode)-nint(symph(:,imode)) ))<0.05_dp)then
       integer_characters(imode)=1
     end if
   end do
   if(sum(integer_characters(:))==3*natom)exit
 end do !itol

!DEBUG
!write(std_out,*)' dfpt_symph : degeneracy=',degeneracy(:)
!ENDDEBUG

 write(message,'(a,a,es8.2,a)')ch10,&
& ' Analysis of degeneracies and characters (maximum tolerance=',ntol*tol6,' a.u.)'
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 do imode=1,3*natom
   if(degeneracy(imode)/=0)then
     write(message,'(a,i4)') ' Symmetry characters of vibration mode #',imode
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,message,'COLL')
     if(degeneracy(imode)>=2)then
       if(degeneracy(imode)==2) write(message,'(a,i4)') &
&       '        degenerate with vibration mode #',imode+1
       if(degeneracy(imode)>=3) write(message,'(a,i4,a,i4)') &
&       '       degenerate with vibration modes #',imode+1,' to ',imode+degeneracy(imode)-1
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if
     do jj=1,(nsym-1)/16+1
       write(message,'(16f5.1)') (symph(isym,imode),isym=(jj-1)*16+1,min(nsym,jj*16))
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end do
   end if
 end do !imode

 ABI_DEALLOCATE(degeneracy)
 ABI_DEALLOCATE(integer_characters)
 ABI_DEALLOCATE(eigvtr)
 ABI_DEALLOCATE(redvtr)
 ABI_DEALLOCATE(redvec)
 ABI_DEALLOCATE(symph)
 ABI_DEALLOCATE(symind)

end subroutine dfpt_symph
!!***
