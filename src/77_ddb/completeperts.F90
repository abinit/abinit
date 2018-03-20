!{\src2tex{textfont=tt}}
!!****f* ABINIT/completeperts
!!
!! NAME
!! completeperts
!!
!! FUNCTION
!!  Complete perturbations wrt atoms and reduced directions
!!  for a fixed qpoint. Normally there is a test in read_gkk which guarantees
!!  that enough irreducible perturbations are present to generate everything.
!!  h1_mat_el is first squared, making a (ipert,jpert) matrix which has the same
!!  symmetry properties as the dynamical matrix.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2018 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  Cryst<crystal_t>=Info on the unit cell and symmetries.
!!   nbranch=Number of phonon branches.
!!   nFSband=Number of bands in H1 matrix elements.
!!   nkpt=Number of k-points in matrix elements.
!!   nsppol=Number of independent spin polarizations.
!!   gkk_flag = flags for presence of gkk matrix elements
!!   h1_mat_el = irreducible matrix elements to be completed and squared
!!   qpt = qpoint
!!   symq = flags for symmetry elements conserving the present qpoint
!!   tnons = translation vectors associated with symops
!!
!! OUTPUT
!!   h1_mat_el_sq = irreducible matrix elements squared and completed
!!   gkk_flag = changed on output
!!
!! PARENTS
!!      read_gkk
!!
!! CHILDREN
!!      d2sym3
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine completeperts(Cryst,nbranch,nFSband,nkpt,nsppol,gkk_flag,h1_mat_el,h1_mat_el_sq,&
&   qpt,symq,qtimrev)

 use defs_basis
 use defs_elphon
 use m_profiling_abi
 use m_errors

 use m_dynmat,     only : d2sym3
 use m_crystal,    only : crystal_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'completeperts'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: qtimrev,nbranch,nFSband,nkpt,nsppol
 type(crystal_t),intent(in) :: Cryst
!arrays
 integer,intent(in) :: symq(4,2,Cryst%nsym)
 integer,intent(inout) :: gkk_flag(nbranch,nbranch,nkpt,nsppol)
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(in) :: h1_mat_el(2,nFSband**2,nbranch,nkpt,nsppol)
 real(dp),intent(out) :: h1_mat_el_sq(2,nFSband**2,nbranch**2,nkpt,nsppol)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,iatom1,iatom2,ibb,idir1,idir2,ipert1,ipert2,isppol,mpert,natom
 real(dp) :: im1,im2,re1,re2,res 
 character(len=500) :: msg
!arrays
 integer,allocatable :: tmpflg(:,:,:,:)
 real(dp),allocatable :: tmpval(:,:,:,:,:)

! *************************************************************************

!WARNING! Stupid patch in d2sym3 imposes these matrices to have size natom+2
 natom = Cryst%natom
 mpert = natom+2

 ABI_ALLOCATE(tmpflg,(3,mpert,3,mpert))
 ABI_ALLOCATE(tmpval,(2,3,mpert,3,mpert))

 h1_mat_el_sq = zero

 write (std_out,*) ' completeperts: shape(h1_mat_el_sq) = ', shape(h1_mat_el_sq)

 do isppol=1,nsppol
   write(std_out,*)'completeperts: isppol = ', isppol
!  
   do ikpt_phon=1,nkpt
     do ibb=1,nFSband**2
!      
       tmpval = zero
       tmpflg = 0
       ! for a fixed k (q) band and sppol construct the gamma matrix for (3 natom)^2 perturbation pairs
       do iatom1=1,natom
         do idir1=1,3
           ipert1 = (iatom1-1)*3+idir1
           if (gkk_flag(ipert1,ipert1,ikpt_phon,isppol) < 0) cycle
           re1 = h1_mat_el(1,ibb,ipert1,ikpt_phon,isppol)
           im1 = h1_mat_el(2,ibb,ipert1,ikpt_phon,isppol)

           do iatom2=1,natom
             do idir2=1,3
               ipert2 = (iatom2-1)*3+idir2
               if (gkk_flag(ipert2,ipert2,ikpt_phon,isppol) < 0) cycle
               tmpflg(idir1,iatom1,idir2,iatom2) = 1
               re2 = h1_mat_el(1,ibb,ipert2,ikpt_phon,isppol)
               im2 = h1_mat_el(2,ibb,ipert2,ikpt_phon,isppol)
!              
!              conjg(h1_mat_el_2) * h1_mat_el_1
               res =  re1*re2 + im1*im2
               tmpval(1,idir1,iatom1,idir2,iatom2) =  res
               res =  re1*im2 - im1*re2
               tmpval(2,idir1,iatom1,idir2,iatom2) = res

             end do !idir2 
           end do !iatom2
         end do !idir1
       end do !iatom1

       ! matrix is symmetrized like a dynamical matrix. No change of band or k
       !  in here. This should be checked (if we have to restrict further the symmetry operations)
       call d2sym3(tmpflg,tmpval,Cryst%indsym,mpert,natom,Cryst%nsym,qpt,symq,Cryst%symrec,Cryst%symrel,qtimrev)

       if (sum(tmpflg(:,1:natom,:,1:natom)) /= 3*natom*3*natom) then
         write(msg,'(3a,4i0)')&
&         'A perturbation is missing after completion with d2sym3',ch10,&
&         'tmpflg, ikpt_phon, isppol: ',tmpflg,ikpt_phon,isppol
         MSG_ERROR(msg)
       end if
!      
!      Save values for calculation of |gkk|^2
       do iatom1=1,natom
         do idir1=1,3
           ipert1 = (iatom1-1)*3+idir1
           do iatom2=1,natom
             do idir2=1,3
!              
!              mjv 29/10/2007 ipert2 now contains the composite index ip1*nperts+ip2
               ipert2 = (iatom2-1)*3 + idir2 + (ipert1-1)*3*natom
               h1_mat_el_sq(1,ibb,ipert2,ikpt_phon,isppol) = pi*tmpval(1,idir2,iatom2,idir1,iatom1)
               h1_mat_el_sq(2,ibb,ipert2,ikpt_phon,isppol) = pi*tmpval(2,idir2,iatom2,idir1,iatom1)
             end do
           end do
         end do
       end do
!      
     end do !end ibb band dos
!    
!    Set flags.
     do ipert1=1,3*natom
       do ipert2=1,3*natom
         if (gkk_flag(ipert2,ipert1,ikpt_phon,isppol) < 0) gkk_flag(ipert2,ipert1,ikpt_phon,isppol) = 1
       end do
     end do

   end do !end kpt_phon do
 end do !end sppol do

 ABI_DEALLOCATE(tmpflg)
 ABI_DEALLOCATE(tmpval)

end subroutine completeperts
!!***
