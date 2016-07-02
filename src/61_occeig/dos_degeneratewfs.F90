!{\src2tex{textfont=tt}}
!!****f* ABINIT/dos_degeneratewfs
!! NAME
!! dos_degeneratewfs
!!
!! FUNCTION
!! Average the contribution to a generalised DOS (that includes weighting factors from matrix elements)
!! from degenerate wavefunctions.
!! Warning : with a negative value of degeneracy_tol, this routine is desactivated. See KNOWN_PROBLEMS.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2016 ABINIT group (MVer,XG,SM,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dos_fractions_in(nkpt,mband,nsppol,ndos)= contribution to generalized DOS from each wavefunction
!!  mband =maximum number of electronic bands for each kpoint
!!  nband(nkpt*nsppol)  =number of electronic bands for each kpoint
!!  nkpt         =number of irreducible kpoints
!!  eigen(mband*nkpt*nsppol)=eigenvalues at irred kpoints
!!  ndos= number of components of dos_fractions_in and dos_fractions_out
!!
!! OUTPUT
!!  dos_fractions_out= contribution to generalized DOS from each wavefunction
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dos_degeneratewfs(dos_fractions_in,dos_fractions_out,eigen,mband,nband,ndos,nkpt,nsppol)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dos_degeneratewfs'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,ndos,nkpt,nsppol
 integer,intent(in) :: nband(nkpt*nsppol)
!arrays
 real(dp),intent(in) :: dos_fractions_in(nkpt,mband,nsppol,ndos)
 real(dp),intent(out) :: dos_fractions_out(nkpt,mband,nsppol,ndos)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol)

!Local variables-------------------------------
!scalars
 integer :: degeneracy_index,degeneracy_min,iband,idos,ikpt,index_eig,isppol,nband_ksp
!integer :: test_done
 real(dp) :: average,degeneracy_tol=-tol6  ! With a negative value, the routine is desactivated ...
 logical :: degeneracy_flag


! *********************************************************************

!DEBUG
!write(std_out,*) 'dos_degeneratewfs : enter'
!test_done=0
!ENDDEBUG
 

!Initialization of dos_fractions_out
 dos_fractions_out(:,:,:,:)=dos_fractions_in(:,:,:,:)

!Average on the degenerate contributions
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband_ksp=nband(ikpt+(isppol-1)*nkpt)
     index_eig=mband*(ikpt-1 + nkpt*(isppol-1))
     if(nband_ksp>1)then
       degeneracy_index=0 ; degeneracy_min=0
       do iband=2,nband_ksp
!        DEBUG
!        if(test_done==20)exit
!        ENDDEBUG
         degeneracy_flag=(eigen(iband+index_eig)-eigen(iband-1+index_eig) < degeneracy_tol)
         if(degeneracy_flag)degeneracy_index=degeneracy_index+1

!        A non-zero value of degeneracy_min will signal that it is time to perform the averaging
         if(degeneracy_index/=0 .and. (.not.degeneracy_flag))degeneracy_min=iband-degeneracy_index-1
         if(degeneracy_index/=0 .and. iband==nband_ksp)      degeneracy_min=iband-degeneracy_index

!        Perform average over all degenerate states at the end of a series of degenerate states
         if(degeneracy_min/=0)then
           do idos=1,ndos
             average=sum(dos_fractions_in(ikpt,degeneracy_min:degeneracy_min+degeneracy_index,isppol,idos))&
&             /real(degeneracy_index+1)
             dos_fractions_out(ikpt,degeneracy_min:degeneracy_min+degeneracy_index,isppol,idos)=average
           end do
!          DEBUG     
!          test_done=test_done+1
!          write(std_out,*)' dos_degeneratewfs '
!          write(std_out,*)' ikpt,isppol,iband,degeneracy_min,degeneracy_max=',&
!          &                      ikpt,isppol,iband,degeneracy_min,degeneracy_min+degeneracy_index
!          write(std_out,*)' eigen(1+index_eig:nband_ksp+index_eig)=',eigen(1+index_eig:nband_ksp+index_eig)
!          write(std_out,*)' dos_fractions_in(ikpt,:,isppol,idos)=',dos_fractions_in(ikpt,:,isppol,idos)
!          write(std_out,*)' dos_fractions_out(ikpt,:,isppol,idos)=',dos_fractions_out(ikpt,:,isppol,idos)
!          ENDDEBUG   

!          Reset degeneracy_index and degeneracy_min
           degeneracy_index=0 ; degeneracy_min=0
         end if
       end do ! iband
       
     end if
   end do ! ikpt
 end do ! isppol

!DEBUG
!write(std_out,*)' dos_degeneratewfs : exit '
!ENDDEBUG

end subroutine dos_degeneratewfs
!!***
