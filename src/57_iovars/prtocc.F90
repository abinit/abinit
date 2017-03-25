!{\src2tex{textfont=tt}}
!!****f* ABINIT/prtocc
!!
!! NAME
!! prtocc
!!
!! FUNCTION
!! Print the content of occ.
!! Due to the need to distinguish between different k-points and
!! different spin polarisations, prttagm.f cannot be used.
!! So, need a dedicated routine.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables
!!  iout=unit number for echoed output
!!  jdtset_(0:ndtset_alloc)=list of dataset indices.
!!  ndtset_alloc=govern second dimension of intarr and dprarr
!!  prtvol_glob= if 0, minimal output volume, if 1, no restriction.
!!  results_out(0:ndtset_alloc)=<type results_out_type>contains the results
!!   needed for outvars, including occ, an evolving variable
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      outvar_o_z
!!
!! CHILDREN
!!      appdig
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine prtocc(dtsets,iout,jdtset_,ndtset_alloc,prtvol_glob,results_out)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_results_out

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtocc'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,ndtset_alloc,prtvol_glob
!arrays
 integer,intent(in) :: jdtset_(0:ndtset_alloc)
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
 type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)

!Local variables-------------------------------
 character(len=*), parameter :: f_occ    ="(1x,a16,1x,(t22,6f10.6))"
 character(len=*), parameter :: f_occa   ="(1x,a16,a,1x,(t22,6f10.6))"
 character(len=*), parameter :: token='occ'
!scalars
 integer,parameter :: nkpt_max=50
 integer :: generic,iban,idtset,ikpsp,ikpt,isppol,jdtset,multi,multi_nband
 integer :: multi_nkpt,multi_nsppol,multi_occopt,nban,nkpt,nkpt_eff
 integer :: multi_tsmear
 integer :: print,tnkpt
 character(len=4) :: appen
 character(len=500) :: message

! *************************************************************************

 if(ndtset_alloc<1)then
   write(message, '(a,i0,a)' )' ndtset_alloc=',ndtset_alloc,', while it should be >= 1.'
   MSG_BUG(message)
 end if

 if(ndtset_alloc>9999)then
   write(message, '(a,i0,a)' )' ndtset_alloc=',ndtset_alloc,', while it must be lower than 100.'
   MSG_BUG(message)
 end if

!It is important to take iscf into account, since when it is -2, occupation numbers must be ignored

 multi_occopt=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%occopt/=dtsets(idtset)%occopt .and. &
&     dtsets(idtset)%iscf/=-2 )multi_occopt=1
   end do
 end if

 multi_tsmear=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%tsmear/=dtsets(idtset)%tsmear .and. &
&     dtsets(idtset)%iscf/=-2 )multi_tsmear=1
   end do
 end if

 multi_nkpt=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nkpt/=dtsets(idtset)%nkpt .and. dtsets(idtset)%iscf/=-2 )multi_nkpt=1
   end do
 end if
 if(multi_nkpt==0)nkpt=dtsets(1)%nkpt

 multi_nsppol=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nsppol/=dtsets(idtset)%nsppol .and. &
&     dtsets(idtset)%iscf/=-2 )multi_nsppol=1
   end do
 end if

 if(multi_nsppol==0 .and. multi_nkpt==0)then
   multi_nband=0
   if(ndtset_alloc>1)then
     do idtset=1,ndtset_alloc
       if(dtsets(idtset)%iscf/=-2)then
         do ikpsp=1,dtsets(1)%nkpt*dtsets(1)%nsppol
           if(dtsets(1)%nband(ikpsp)/=dtsets(idtset)%nband(ikpsp))multi_nband=1
         end do
       end if
     end do
   end if
 else
   multi_nband=1
 end if

!There is a possibility of a generic occupation-number set if
!multi_occopt==0 and multi_nband==0
 multi=1
 if(multi_occopt==0 .and. multi_nband==0) then
   nban=sum(dtsets(1)%nband(1:dtsets(1)%nsppol*dtsets(1)%nkpt))
   multi=0
   if(ndtset_alloc>1)then
     do idtset=1,ndtset_alloc
       if(dtsets(idtset)%iscf/=-2)then
!        nban counts all bands and kpoints and spins: see above
         do iban=1,nban
!          Use of tol8, because the format for multi=1 is f16.6, so will not
!          discriminate between relative values, or absolute values that
!          agree within more than 6 digits
           if( abs(results_out(1)%occ(iban,1)-results_out(idtset)%occ(iban,1)) > tol8) multi=1
         end do
       end if
     end do
   end if
 end if

!At this stage, if multi==1, the occ must be printed
!if multi==0, then it might be that we have the default values.
!Since the default is all zeros, it only happens when iscf=-2
!Also initialize the number of a idtset that can be used as generic
!(this might not be the case for idtset=1 !)

 generic=0
 print=0
 do idtset=1,ndtset_alloc
   if(dtsets(idtset)%iscf/=-2)then
     print=1
     generic=idtset
   end if
 end do

!Now, print in the generic occupation-number set case.
 if(print==1 .and. multi==0)then
!  Might restrict the number of k points to be printed
   tnkpt=0
   nkpt_eff=dtsets(1)%nkpt
   if(prtvol_glob==0 .and. nkpt_eff>nkpt_max)then
     nkpt_eff=nkpt_max
     tnkpt=1
   end if
!  The quantity of data to be output vary with occopt
   if(dtsets(1)%occopt>=2)then
     iban=1
     do isppol=1,dtsets(1)%nsppol
       do ikpt=1,nkpt_eff
         ikpsp=ikpt+dtsets(1)%nkpt*(isppol-1)
         nban=dtsets(generic)%nband(ikpsp)
         if(ikpsp==1)then
           write(iout, '(1x,a16,1x,(t22,6f10.6))' )&
&           token,results_out(generic)%occ(iban:iban+nban-1,1)
         else
           write(iout, '((t22,6f10.6))' )results_out(generic)%occ(iban:iban+nban-1,1)
         end if
         iban=iban+nban
       end do
       if(tnkpt==1) write(iout,'(23x,a)' ) &
&       'prtocc : prtvol=0, do not print more k-points.'
     end do
   else
!    The number of bands is identical for all k points and spin
     nban=dtsets(generic)%nband(1)
     write(iout, '(1x,a16,1x,(t22,6f10.6))' )&
&     token,results_out(generic)%occ(1:nban,1)
!    if occopt==1, the occ might differ with the spin
     if(dtsets(1)%nsppol/=1 .and. dtsets(1)%occopt==1)then
       write(iout,'((t22,6f10.6))')results_out(generic)%occ(nban*dtsets(1)%nkpt+1:&
&       nban*dtsets(1)%nkpt+nban,1)
     end if
   end if
 end if

!Now, print in the other cases
 if(print==1 .and. multi==1)then
   do idtset=1,ndtset_alloc
!    Might restrict the number of k points to be printed
     tnkpt=0
     nkpt_eff=dtsets(idtset)%nkpt
     if(prtvol_glob==0 .and. nkpt_eff>nkpt_max)then
       nkpt_eff=nkpt_max
       tnkpt=1
     end if
     if(dtsets(idtset)%iscf/=-2)then
       jdtset=jdtset_(idtset)
       call appdig(jdtset,'',appen)
!      The quantity of data to be output vary with occopt
       if(dtsets(idtset)%occopt>=2)then
         iban=1
         do isppol=1,dtsets(idtset)%nsppol
           do ikpt=1,nkpt_eff
             ikpsp=ikpt+dtsets(idtset)%nkpt*(isppol-1)
             nban=dtsets(idtset)%nband(ikpsp)
             if(ikpsp==1)then
               write(iout, '(1x,a16,a,1x,(t22,6f10.6))' )&
&               token,appen,results_out(idtset)%occ(iban:iban+nban-1,1)
             else
               write(iout, '((t22,6f10.6))' )results_out(idtset)%occ(iban:iban+nban-1,1)
             end if
             iban=iban+nban
           end do
           if(tnkpt==1) write(iout,'(23x,a)' ) &
&           'prtocc : prtvol=0, do not print more k-points.'
         end do
       else
!        The number of bands is identical for all k points and spin
         nban=dtsets(idtset)%nband(1)
         write(iout, '(1x,a16,a,1x,(t22,6f10.6))' )&
&         token,appen,results_out(idtset)%occ(1:nban,1)
!        if occopt==1, the occ might differ with the spin
         if(dtsets(idtset)%nsppol/=1 .and. dtsets(idtset)%occopt==1)then
           write(iout, '((t22,6f10.6))' ) &
&           results_out(idtset)%occ(nban*dtsets(idtset)%nkpt+1:nban*dtsets(idtset)%nkpt+nban,1)
         end if
       end if
     end if
!    Endloop on idtset
   end do
 end if

end subroutine prtocc
!!***
