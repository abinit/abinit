!{\src2tex{textfont=tt}}
!!****f* ABINIT/graph
!! NAME
!! graph
!!
!! FUNCTION
!! Writing  of the gnuplot script to show the computed part
!! of Bader surface with lines
!!
!! COPYRIGHT
!! Copyright (C) 2002-2018 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  untg = unit number of the file on which the info is written
!!  unts = unit number of the file from which the Bader surface is read
!!
!! OUTPUT
!!  (written in the untg file)
!!
!! PARENTS
!!      drvaim
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine graph(unts,untg)

 use m_profiling_abi

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'graph'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: untg,unts

!Local variables ------------------------------
!scalars
 integer :: ii,indx,jj,nphi,nth
 real(dp),parameter :: snull=1.d-6
 real(dp) :: phimax,phimin,ss,thmax,thmin
!arrays
 real(dp) :: xorig(3)
 real(dp),allocatable :: phi(:),rr(:,:),th(:)

! *********************************************************************

 rewind(unts)
 read(unts,*) indx, xorig(1:3)
 read(unts,*) nth, thmin, thmax
 read(unts,*) nphi, phimin, phimax
 ABI_ALLOCATE(th,(nth))
 ABI_ALLOCATE(phi,(nphi))
 ABI_ALLOCATE(rr,(nth,nphi))
 do ii=1,nth
   do jj=1,nphi
     read(unts,*) th(ii),phi(jj),rr(ii,jj),ss
   end do
 end do

!end of reading

 write(untg,*) 'reset'
 write(untg,*) 'set st d l'
 write(untg,*) 'set ticslevel 0'
 write(untg,*) 'set title ''Bader surface'' '
 write(untg,*) 'splot ''-'' using ($3*sin($1)*cos($2)):($3*sin($1)*sin($2)):($3*cos($1)) notitle'
 do ii=1,nth
   do jj=1,nphi
     write(untg,'(2F12.8,E16.8)') th(ii),phi(jj),rr(ii,jj)
   end do
   if ((ii==nth).and.(jj==nphi)) then
     cycle
   else
     write(untg,*)
   end if
 end do

end subroutine graph
!!***
