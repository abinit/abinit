!{\src2tex{textfont=tt}}
!!****f* ABINIT/alignph
!!
!! NAME
!! alignph
!!
!! FUNCTION
!! Construct linear combinations of the phonon eigendisplacements
!! of degenerate modes in order to align the mode effective charges
!! along the axes of the cartesian frame.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! amu(ntypat)=mass of the atoms (atomic mass unit)
!! displ(2,3*natom,3*natom)=
!! the displacements of atoms in cartesian coordinates.
!! The first index means either the real or the imaginary part,
!! The second index runs on the direction and the atoms displaced
!! The third index runs on the modes.
!! d2cart(2,3,mpert,3,mpert)=
!!  dynamical matrix, effective charges, dielectric tensor,....
!!  all in cartesian coordinates
!! mpert =maximum number of ipert
!! natom=number of atoms in unit cell
!! ntypat=number of types of atoms
!! phfrq(3*natom)=phonon frequencies (square root of the dynamical
!!  matrix eigenvalues, except if these are negative, and in this
!!  case, give minus the square root of the absolute value
!!  of the matrix eigenvalues). Hartree units.
!! typat(natom)=integer label of each type of atom (1,2,...)
!!
!! OUTPUT
!! displ(2,3*natom,3*natom)=
!! the displacements of atoms in cartesian coordinates.
!! The eigendisplacements of degenerate modes have been aligned along
!! the cartesian axes.
!!
!! PARENTS
!!      ddb_diel
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine alignph(amu,displ,d2cart,mpert,natom,ntypat,phfrq,typat)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'alignph'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,ntypat
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: amu(ntypat),d2cart(2,3,mpert,3,mpert),phfrq(3*natom)
 real(dp),intent(inout) :: displ(2,3*natom,3*natom)

!Local variables -------------------------
!scalars
 integer :: i1,idir1,idir2,ii,imode,imodex,imodey,imodez,ipert1
 real(dp) :: theta
!arrays
 integer,allocatable :: deg(:)
 real(dp) :: zvec(3,3),zvect(3,3)
 real(dp),allocatable :: modez(:,:,:),modezabs(:),oscstr(:,:,:),vec(:,:),vect(:,:)

! *********************************************************************

!DEBUG
! write(std_out,*)'alignph : enter'
!ENDDEBUG

!Get the oscillator strength and mode effective charge for each mode
 ABI_ALLOCATE(oscstr,(2,3,3*natom))
 ABI_ALLOCATE(modez,(2,3,3*natom))
 ABI_ALLOCATE(modezabs,(3*natom))
 ABI_ALLOCATE(vec,(3*natom,3))
 ABI_ALLOCATE(vect,(3*natom,3))
 ABI_ALLOCATE(deg,(3*natom))

 write(std_out,'(a,a)')ch10,' alignph : before modifying the eigenvectors, mode number and mode effective charges :'
 do imode=1,3*natom
   modezabs(imode)=zero
   do ii=1,2
     do idir2=1,3
       oscstr(ii,idir2,imode)=zero
       modez(ii,idir2,imode)=zero
       do idir1=1,3
         do ipert1=1,natom
           i1=idir1+(ipert1-1)*3
           oscstr(ii,idir2,imode)=oscstr(ii,idir2,imode)+&
&           displ(ii,i1,imode)*&
&           d2cart(1,idir1,ipert1,idir2,natom+2)
           modez(ii,idir2,imode)=modez(ii,idir2,imode)+&
&           displ(ii,i1,imode)*&
&           d2cart(1,idir1,ipert1,idir2,natom+2)*&
&           sqrt(amu(typat(ipert1))*amu_emass)
         end do
       end do
       if(abs(modez(ii,idir2,imode))>modezabs(imode))modezabs(imode)=abs(modez(ii,idir2,imode))
     end do
   end do
   write(std_out,'(i4,3f16.6)')imode,modez(1,:,imode)
 end do

!Find degenerate modes with non-zero mode effective charge
 imode = 0
 do while (imode < 3*natom)
   imode = imode + 1
   if (imode == 3*natom) then
     deg(imode) = 1
   else if (abs(phfrq(imode) - phfrq(imode+1)) > tol6 .or. modezabs(imode)<tol8 .or. modezabs(imode+1)<tol8) then
!    Differ by phonon frequency or zero mode effective charge
     deg(imode) = 1
   else
     deg(imode) = 2
     if (imode < 3*natom - 1) then
       if (abs(phfrq(imode) - phfrq(imode+2)) < tol6 .and. modezabs(imode+2)>tol8) then
         deg(imode) = 3
         imode = imode + 1
       end if
     end if
     imode = imode + 1
   end if
 end do


!In case of a degenerate mode, with non-zero mode effective charge, align the mode effective charge vector along
!the axes of the cartesian frame
 imode = 1
 do while (imode <= 3*natom)

   write(std_out,'(a,a,i4,a,i2)')ch10,' Mode number ',imode,' has degeneracy ',deg(imode)
   write(std_out,'(a,3es16.6)') ' Mode effective charge of this mode =',modez(1,:,imode)

   if (deg(imode) == 2) then

!    Optimize on the x direction
     write(std_out,'(a,3es16.6)') ' Mode effective charge of next mode =',modez(1,:,imode+1)
     if (abs(modez(1,1,imode)) > tol8) then
       theta = atan(-modez(1,1,imode+1)/modez(1,1,imode))
       vec(:,1) = displ(1,:,imode)
       vec(:,2) = displ(1,:,imode+1)
       displ(1,:,imode) = cos(theta)*vec(:,1) - sin(theta)*vec(:,2)
       displ(1,:,imode+1) = sin(theta)*vec(:,1) + cos(theta)*vec(:,2)
     end if

   else if (deg(imode) == 3) then

     write(std_out,'(a,3es16.6)') ' Mode effective charge of next mode =',modez(1,:,imode+1)
     write(std_out,'(a,3es16.6)') ' Mode effective charge of next-next mode =',modez(1,:,imode+2)

!    Before mixing them, select the mode-effective charge vectors as being predominently "x", "y" or "z" type.
     if(abs(modez(1,1,imode))>abs(modez(1,2,imode))-tol12 .and. &
&     abs(modez(1,1,imode))>abs(modez(1,3,imode))-tol12) then
       imodex=imode
       if(abs(modez(1,2,imode+1))>abs(modez(1,3,imode+1))-tol12)then
         imodey=imode+1 ; imodez=imode+2
       else
         imodez=imode+1 ; imodey=imode+2
       end if
     else if(abs(modez(1,2,imode))>abs(modez(1,1,imode))-tol12 .and. &
&       abs(modez(1,2,imode))>abs(modez(1,3,imode))-tol12) then
       imodey=imode
       if(abs(modez(1,1,imode+1))>abs(modez(1,3,imode+1))-tol12)then
         imodex=imode+1 ; imodez=imode+2
       else
         imodez=imode+1 ; imodex=imode+2
       end if
     else
       imodez=imode
       if(abs(modez(1,1,imode+1))>abs(modez(1,2,imode+1))-tol12)then
         imodex=imode+1 ; imodey=imode+2
       else
         imodey=imode+1 ; imodex=imode+2
       end if
     end if
     vec(:,1)=displ(1,:,imodex)
     vec(:,2)=displ(1,:,imodey)
     vec(:,3)=displ(1,:,imodez)
     zvec(:,1)=modez(1,:,imodex)
     zvec(:,2)=modez(1,:,imodey)
     zvec(:,3)=modez(1,:,imodez)


!    Optimize along x : does the first vector has a component along x ?
     if (abs(zvec(1,1)) > tol8) then
!      Optimize on the (1,2) pair of modes along x
       theta = atan(-zvec(1,2)/zvec(1,1))
       zvect(:,:)=zvec(:,:)
       zvec(:,1) = cos(theta)*zvect(:,1) - sin(theta)*zvect(:,2)
       zvec(:,2) = sin(theta)*zvect(:,1) + cos(theta)*zvect(:,2)
       vect(:,:)=vec(:,:)
       vec(:,1) = cos(theta)*vect(:,1) - sin(theta)*vect(:,2)
       vec(:,2) = sin(theta)*vect(:,1) + cos(theta)*vect(:,2)
!      Optimize on the (1,3) pair of modes along x
       theta = atan(-zvec(1,3)/zvec(1,1))
       zvect(:,:)=zvec(:,:)
       zvec(:,1) = cos(theta)*zvect(:,1) - sin(theta)*zvect(:,3)
       zvec(:,3) = sin(theta)*zvect(:,1) + cos(theta)*zvect(:,3)
       vect(:,:)=vec(:,:)
       vec(:,1) = cos(theta)*vect(:,1) - sin(theta)*vect(:,3)
       vec(:,3) = sin(theta)*vect(:,1) + cos(theta)*vect(:,3)
       if (abs(zvec(2,2)) > tol8) then
!        Optimize on the (2,3) pair of modes along y
         theta = atan(-zvec(2,3)/zvec(2,2))
         zvect(:,:)=zvec(:,:)
         zvec(:,2) = cos(theta)*zvect(:,2) - sin(theta)*zvect(:,3)
         zvec(:,3) = sin(theta)*zvect(:,2) + cos(theta)*zvect(:,3)
         vect(:,:)=vec(:,:)
         vec(:,2) = cos(theta)*vect(:,2) - sin(theta)*vect(:,3)
         vec(:,3) = sin(theta)*vect(:,2) + cos(theta)*vect(:,3)
       end if
!    Likely, the remaining is not needed ... because the vectors have been ordered in x, y, and z major component ...
!    Optimize along x : does the second vector has a component along x ?
     else if(abs(zvec(1,2)) > tol8) then
!      Optimize on the (2,3) pair of modes along x
       theta = atan(-zvec(1,3)/zvec(1,2))
       zvect(:,:)=zvec(:,:)
       zvec(:,2) = cos(theta)*zvect(:,2) - sin(theta)*zvect(:,3)
       zvec(:,3) = sin(theta)*zvect(:,2) + cos(theta)*zvect(:,3)
       vect(:,:)=vec(:,:)
       vec(:,2) = cos(theta)*vect(:,2) - sin(theta)*vect(:,3)
       vec(:,3) = sin(theta)*vect(:,2) + cos(theta)*vect(:,3)
!      Optimize on the (1,3) pair of modes along y
       if (abs(zvec(2,1)) > tol8) then
         theta = atan(-zvec(2,3)/zvec(2,1))
         zvect(:,:)=zvec(:,:)
         zvec(:,1) = cos(theta)*zvect(:,1) - sin(theta)*zvect(:,3)
         zvec(:,3) = sin(theta)*zvect(:,1) + cos(theta)*zvect(:,3)
         vect(:,:)=vec(:,:)
         vec(:,1) = cos(theta)*vect(:,1) - sin(theta)*vect(:,3)
         vec(:,3) = sin(theta)*vect(:,1) + cos(theta)*vect(:,3)
       end if
!    We are left with the pair of vectors (2,3)
     else if (abs(zvec(2,2)) > tol8) then
!      Optimize on the (2,3) pair of modes along y
       theta = atan(-zvec(2,3)/zvec(2,2))
       zvect(:,:)=zvec(:,:)
       zvec(:,2) = cos(theta)*zvect(:,2) - sin(theta)*zvect(:,3)
       zvec(:,3) = sin(theta)*zvect(:,2) + cos(theta)*zvect(:,3)
       vect(:,:)=vec(:,:)
       vec(:,2) = cos(theta)*vect(:,2) - sin(theta)*vect(:,3)
       vec(:,3) = sin(theta)*vect(:,2) + cos(theta)*vect(:,3)
     end if

     displ(1,:,imodex)=vec(:,1)
     displ(1,:,imodey)=vec(:,2)
     displ(1,:,imodez)=vec(:,3)

!    Previous coding, from Marek. Apparently, break the orthogonalization of vectors ...
!    do ii = 1,3
!      coeff(:) = 0._dp
!      if (ii == 1) then
!        jj = 2 ; kk = 3
!      else if (ii == 2) then
!        jj = 1 ; kk = 3
!      else
!        jj = 1 ; kk = 2
!      end if
!      coeff(ii) = 1._dp
!      c1 = modez(1,jj,imode+ii-1)
!      c2 = modez(1,jj,imode+jj-1)
!      c3 = modez(1,jj,imode+kk-1)
!      c4 = modez(1,kk,imode+ii-1)
!      c5 = modez(1,kk,imode+jj-1)
!      c6 = modez(1,kk,imode+kk-1)
!      dtm = c2*c6 - c3*c5
!      if (abs(dtm) > tol8) then
!        coeff(jj) = (c3*c4 - c1*c6)/dtm
!        coeff(kk) = (c1*c5 - c2*c4)/dtm
!      end if
!      mod_ = sqrt(1._dp + coeff(jj)*coeff(jj) + coeff(kk)*coeff(kk))
!      coeff(:) = coeff(:)/mod_
!      displ(1,:,imode+ii-1) = coeff(1)*vec(1,:) + coeff(2)*vec(2,:) + &
!&       coeff(3)*vec(3,:)
!    end do

   end if ! if deg mode

   imode = imode + deg(imode)

 end do

 write(std_out,'(a,a)')ch10,' alignph : after modifying the eigenvectors, mode number and mode effective charges :'
 do imode=1,3*natom
   do ii=1,2
     do idir2=1,3
       modez(ii,idir2,imode)=zero
       do idir1=1,3
         do ipert1=1,natom
           i1=idir1+(ipert1-1)*3
           modez(ii,idir2,imode)=modez(ii,idir2,imode)+&
&           displ(ii,i1,imode)*&
&           d2cart(1,idir1,ipert1,idir2,natom+2)*&
&           sqrt(amu(typat(ipert1))*amu_emass)
         end do
       end do
     end do
   end do
   write(std_out,'(i4,3f16.6)')imode,modez(1,:,imode)
 end do

 ABI_DEALLOCATE(deg)
 ABI_DEALLOCATE(oscstr)
 ABI_DEALLOCATE(modez)
 ABI_DEALLOCATE(modezabs)
 ABI_DEALLOCATE(vec)
 ABI_DEALLOCATE(vect)

end subroutine alignph
!!***
