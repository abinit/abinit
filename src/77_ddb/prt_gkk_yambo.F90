!{\src2tex{textfont=tt}}
!!****f* ABINIT/prt_gkk_yambo
!!
!! NAME
!! prt_gkk_yambo
!!
!! FUNCTION
!! This routine outputs el-phon related quantities for the yambo code at 1
!!   q-point
!!
!! COPYRIGHT
!! Copyright (C) 2009-2016 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt
!.
!!
!! INPUTS
!!  displ_cart = phonon displacement vectors for this q-point in Cartesian coordinates.
!!  displ_red = phonon displacement vectors for this q-point, in reduced coordinates
!!  elph_ds = datastructure containing elphon matrix elements
!!  h1_mat_el = matrix elements of first order hamiltonian for present q-point,
!!     all perturbations
!!  iqptfull = index of present q-point in full array of q-points
!!  irredpert = index of irreducible perturbation (atom displaced)
!!  natom = number of atoms
!!  phfrq = phonon frequencies at present q-point
!!  qptn = q-point we will print for
!!
!! OUTPUT
!!  only writes to a file
!!
!! NOTES
!!
!! PARENTS
!!      read_gkk
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine prt_gkk_yambo(displ_cart,displ_red,kpt_phon,h1_mat_el,iqpt,&
&       natom,nFSband,nkpt_phon,phfrq,qptn)

 use defs_basis
 use m_profiling_abi

 use m_io_tools,   only : get_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prt_gkk_yambo'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,iqpt
 integer,intent(in) :: nFSband,nkpt_phon
 !arrays
 real(dp),intent(in) :: kpt_phon(3,nkpt_phon)
 real(dp),intent(in) :: h1_mat_el(2,nFSband*nFSband,3*natom,nkpt_phon,1)
 real(dp),intent(in) :: phfrq(3*natom)
 real(dp),intent(in) :: displ_cart(2,3*natom,3*natom)
 real(dp),intent(in) :: displ_red(2,3*natom,3*natom)
 real(dp),intent(in) :: qptn(3)

!Local variables-------------------------------
 !scalars
 integer, save :: firsttime=1
 integer :: outunit,ikpt,imode,iband,ibandp,iatom,idir,ibandindex
 integer :: jmode, outunit2, outunit3
 !arrays
 real(dp) :: gkk_mode_dep(2)
! *************************************************************************



!if first time round:
 if (firsttime==1) then
   firsttime=0
!  squash file
   outunit=get_unit()
   open (unit=outunit,file="yambo_elphon_data",status="REPLACE")
   outunit2=get_unit()
   open (unit=outunit2,file="yambo_elphon_gkk_bymode",status="replace")
   outunit3=get_unit()
   open (unit=outunit3,file="yambo_elphon_gkksqtw_bymode",status="replace")

!  write dimensions
   write (outunit,'(a,I6)') 'number of el atoms ', natom
   write (outunit2,'(a,I6)') 'number of el atoms ', natom
   write (outunit3,'(a,I6)') 'number of el atoms ', natom
   write (outunit,'(a,I6)') 'number of ph modes ', 3*natom
   write (outunit2,'(a,I6)') 'number of ph modes ', 3*natom
   write (outunit3,'(a,I6)') 'number of ph modes ', 3*natom
   write (outunit,'(a,I6)') 'number of el bands ', nFSband
   write (outunit2,'(a,I6)') 'number of el bands ', nFSband
   write (outunit3,'(a,I6)') 'number of el bands ', nFSband

!  write k-points
   write (outunit,'(a,I6)') 'number of k-points ', nkpt_phon
   write (outunit2,'(a,I6)') 'number of k-points ', nkpt_phon
   write (outunit3,'(a,I6)') 'number of k-points ', nkpt_phon
   do ikpt=1,nkpt_phon
     write (outunit,'(a,I6,3E20.10)') 'reduced coord kpoint no ', ikpt, kpt_phon(:,ikpt)
     write (outunit2,'(a,I6,3E20.10)') 'reduced coord kpoint no ', ikpt, kpt_phon(:,ikpt)
     write (outunit3,'(a,I6,3E20.10)') 'reduced coord kpoint no ', ikpt, kpt_phon(:,ikpt)
   end do

!  band energies are not accessible this deep in the code: simpler to get them
!  from elsewhere

   close (outunit)
   close (outunit2)
   close (outunit3)
 end if ! first time round

!open file
 outunit=get_unit()
 open (unit=outunit,file="yambo_elphon_data",status="unknown",position="append")

!qpoint
 write (outunit,'(a,I6,3E20.10)') 'reduced coord qpoint no ', iqpt, qptn(:)

!frequencies
 do imode=1,3*natom
   write (outunit,'(a,I6,3E20.10)') 'phonon freq no ', imode, phfrq(imode)
 end do

!displacement vector
 do imode=1,3*natom
   write (outunit,'(a,I6,3E20.10)') 'phonon displ vec no ', imode
   do iatom=1,natom
     write (outunit,'(3(2E20.10,2x))') displ_cart(:,(iatom-1)*3+1:iatom*3,imode)
   end do
 end do

!the beef: matrix elements of the first order hamiltonian for displacement of
!all atoms along all reduced directions
 write (outunit,'(a)') ' matrix elements of all perturbations for this q-point'
 do ikpt=1,nkpt_phon 
   write (outunit,'(a,I6)') ' kpoint ', ikpt
   imode=0
   do iatom=1,natom
     do idir=1,3
       imode=imode+1
       write (outunit,'(a,I6,I6)') ' atom, direction = ', iatom,idir
       ibandindex=0
       do iband=1,nFSband
         do ibandp=1,nFSband
           ibandindex=ibandindex+1
           write (outunit,'(a,I6,I6,2E20.10)') ' mat el for n,np ', iband,ibandp,&
&           h1_mat_el(:,ibandindex,imode,ikpt,1)
         end do !bandp
       end do !band
     end do !dir
   end do !atom
 end do

!blank line
 write (outunit,*)
 close (outunit)

 outunit2=get_unit()
 open (unit=outunit2,file="yambo_elphon_gkk_bymode",status="unknown",position="append")
 outunit3=get_unit()
 open (unit=outunit3,file="yambo_elphon_gkksqtw_bymode",status="unknown",position="append")

!qpoint
 write (outunit2,'(a,I6,3E20.10)') 'reduced coord qpoint no ', iqpt, qptn(:)
 write (outunit3,'(a,I6,3E20.10)') 'reduced coord qpoint no ', iqpt, qptn(:)

!print out mode-dependent matrix elements
 write (outunit2,'(a)') ' matrix elements of all phonon modes for this q-point'
 write (outunit3,'(a)') ' 1/w**1/2 times matrix elements of all phonon modes for this q-point'
 do ikpt=1,nkpt_phon
   write (outunit2,'(a,I6)') ' kpoint ', ikpt
   write (outunit3,'(a,I6)') ' kpoint ', ikpt
   ibandindex=0
   do iband=1,nFSband
     do ibandp=1,nFSband
       ibandindex=ibandindex+1
       write (outunit2,'(a,I6,I6)') ' el bands n,np ', iband,ibandp
       write (outunit3,'(a,I6,I6)') ' el bands n,np ', iband,ibandp
       do imode=1,3*natom
!        gkk_mode_dep = cg_zdotc(3*natom,displ_red(:,:,imode),h1_mat_el(:,ibandindex,:,ikpt,1))
         gkk_mode_dep = zero
         do jmode=1,3*natom
           gkk_mode_dep(1) = gkk_mode_dep(1) &
&           + displ_red(1,jmode,imode)*h1_mat_el(1,ibandindex,jmode,ikpt,1) &
&           + displ_red(2,jmode,imode)*h1_mat_el(2,ibandindex,jmode,ikpt,1)
           gkk_mode_dep(2) = gkk_mode_dep(2) &
&           + displ_red(1,jmode,imode)*h1_mat_el(2,ibandindex,jmode,ikpt,1) &
&           - displ_red(2,jmode,imode)*h1_mat_el(1,ibandindex,jmode,ikpt,1)
         end do
         write (outunit2,'(a,I6,2E20.10)') ' mat el for phonon mode num = ', imode, gkk_mode_dep
         write (outunit3,'(a,I6,2E20.10)') ' 1/w**1/2 * mat el for phonon mode num = ', &
&         imode, gkk_mode_dep/sqrt(two*abs(phfrq(imode))+tol10)
       end do !imode
     end do !bandp
   end do !band
 end do
!blank line
 write (outunit2,*)
 write (outunit3,*)

 close (outunit2)
 close (outunit3)

end subroutine prt_gkk_yambo
!!***
