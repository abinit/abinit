!{\src2tex{textfont=tt}}
!!****f* ABINIT/symgamma
!!
!! NAME
!! symgamma
!!
!! FUNCTION
!!  Symmetrize perturbations wrt atoms and reduced directions
!!  for a fixed qpoint.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2017 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   elph_ds = datastructure for elph data (dimensions and eventually data)
!!   kphon_full2full = mapping btw kpoints under symops
!!   h1_mat_el = irreducible matrix elements to be completed
!!   indsym = mapping of atoms under symops
!!   natom = number of atoms
!!   nsym = number of syms
!!   symq = flags for symmetry elements conserving the present qpoint
!!   symrec = symmetry operations in reciprocal space
!!
!! OUTPUT
!!   h1_mat_el = changed on output
!!
!! NOTES
!!  scheduled for destruction: not called at present (2/2010) so should eventually dissappear
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


subroutine symgamma(elph_ds,kphon_full2full,h1_mat_el,&
&   indsym,natom,nsym,symq,symrec)

 use defs_basis
 use defs_elphon
 use m_profiling_abi
 use m_io_tools

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symgamma'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsym
 type(elph_type),intent(in) :: elph_ds
!arrays
 integer,intent(in) :: kphon_full2full(2,nsym,elph_ds%k_phon%nkpt),indsym(4,nsym,natom)
 integer,intent(in) :: symq(4,2,nsym),symrec(3,3,nsym)
 real(dp),intent(inout) :: h1_mat_el(2,elph_ds%nFSband,elph_ds%nFSband,elph_ds%nbranch,elph_ds%k_phon%nkpt)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,ipreatom,isym,isymatom,itim
 integer :: nsymperts,symikpt_phon
 real(dp) :: timsign
!arrays
 real(dp) :: dsymrec(3,3)
 real(dp) :: sym_mat_el(2,elph_ds%nFSband,elph_ds%nFSband,elph_ds%nbranch,elph_ds%k_phon%nkpt)

! *************************************************************************

!2. symmetrize the whole set of perturbations
 sym_mat_el(:,:,:,:,:) = zero

 nsymperts = 0
!
!symrel(isym) sends ipreatom onto isymatom
!symrec(isym) sends ipredir onto isymdir
!symrec(isym) sends ikpt_phon onto symikpt_phon
!
 do isym=1,nsym

   do itim=0,1
     if (symq(4,itim+1,isym) == 0) cycle
     timsign = one - two*itim

     dsymrec(:,:) = timsign*dble(symrec(:,:,isym))

     nsymperts = nsymperts+1

!    loop over image perturbations
     do isymatom=1,natom
       ipreatom = indsym(4,isym,isymatom)
       write(std_out,*) ' symgamma : ', isym,itim,isymatom,ipreatom,dsymrec


       do symikpt_phon=1,elph_ds%k_phon%nkpt
         ikpt_phon = kphon_full2full(itim+1,isym,symikpt_phon)

!        Real part
         sym_mat_el (1,:,:,3*(isymatom-1)+1,symikpt_phon) = &
&         sym_mat_el (1,:,:,3*(isymatom-1)+1,symikpt_phon) &
&         + dsymrec(1,1)*h1_mat_el(1,:,:,3*(ipreatom-1)+1,ikpt_phon) &
&         + dsymrec(1,2)*h1_mat_el(1,:,:,3*(ipreatom-1)+2,ikpt_phon) &
&         + dsymrec(1,3)*h1_mat_el(1,:,:,3*(ipreatom-1)+3,ikpt_phon)
         sym_mat_el (1,:,:,3*(isymatom-1)+2,symikpt_phon) = &
&         sym_mat_el (1,:,:,3*(isymatom-1)+2,symikpt_phon) &
&         + dsymrec(2,1)*h1_mat_el(1,:,:,3*(ipreatom-1)+1,ikpt_phon) &
&         + dsymrec(2,2)*h1_mat_el(1,:,:,3*(ipreatom-1)+2,ikpt_phon) &
&         + dsymrec(2,3)*h1_mat_el(1,:,:,3*(ipreatom-1)+3,ikpt_phon)
         sym_mat_el (1,:,:,3*(isymatom-1)+3,symikpt_phon) = &
&         sym_mat_el (1,:,:,3*(isymatom-1)+3,symikpt_phon) &
&         + dsymrec(3,1)*h1_mat_el(1,:,:,3*(ipreatom-1)+1,ikpt_phon) &
&         + dsymrec(3,2)*h1_mat_el(1,:,:,3*(ipreatom-1)+2,ikpt_phon) &
&         + dsymrec(3,3)*h1_mat_el(1,:,:,3*(ipreatom-1)+3,ikpt_phon)
!        Imag part
         sym_mat_el (2,:,:,3*(isymatom-1)+1,symikpt_phon) = &
&         sym_mat_el (2,:,:,3*(isymatom-1)+1,symikpt_phon) &
&         +         (dsymrec(1,1)*h1_mat_el(2,:,:,3*(ipreatom-1)+1,ikpt_phon) &
&         +          dsymrec(1,2)*h1_mat_el(2,:,:,3*(ipreatom-1)+2,ikpt_phon) &
&         +          dsymrec(1,3)*h1_mat_el(2,:,:,3*(ipreatom-1)+3,ikpt_phon))
         sym_mat_el (2,:,:,3*(isymatom-1)+2,symikpt_phon) = &
&         sym_mat_el (2,:,:,3*(isymatom-1)+2,symikpt_phon) &
&         +         (dsymrec(2,1)*h1_mat_el(2,:,:,3*(ipreatom-1)+1,ikpt_phon) &
&         +          dsymrec(2,2)*h1_mat_el(2,:,:,3*(ipreatom-1)+2,ikpt_phon) &
&         +          dsymrec(2,3)*h1_mat_el(2,:,:,3*(ipreatom-1)+3,ikpt_phon))
         sym_mat_el (2,:,:,3*(isymatom-1)+3,symikpt_phon) = &
&         sym_mat_el (2,:,:,3*(isymatom-1)+3,symikpt_phon) &
&         +         (dsymrec(3,1)*h1_mat_el(2,:,:,3*(ipreatom-1)+1,ikpt_phon) &
&         +          dsymrec(3,2)*h1_mat_el(2,:,:,3*(ipreatom-1)+2,ikpt_phon) &
&         +          dsymrec(3,3)*h1_mat_el(2,:,:,3*(ipreatom-1)+3,ikpt_phon))
       end do
     end do
   end do
 end do
!end isym and itim do

!commented to use un-symmetrized version of h1_mat_el
 h1_mat_el(:,:,:,:,:) = sym_mat_el(:,:,:,:,:) / nsymperts

end subroutine symgamma
!!***
