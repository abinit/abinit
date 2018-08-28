!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_symkpt
!! NAME
!!  m_symkpt
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group ()
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! TODO
!!  Move it to m_kpts
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_symkpt

 implicit none

 private
!!***

 public :: symkpt
!!***

contains
!!***

!!****f* ABINIT/symkpt
!! NAME
!! symkpt
!!
!! FUNCTION
!! Determines the weights of the k-points for sampling the Brillouin Zone, starting from a first set
!! of weights wtk, and folding it to a new set, by taking into account the symmetries described
!! by symrec, and eventually the time-reversal symmetry.
!! Also compute the number of k points in the reduced set
!! This routine is also used for sampling the q vectors in the Brillouin zone for the computation
!! of thermodynamical properties (from the routine thm9).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (XG,LSI)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! chksymbreak= if 1, will check whether the k point grid is symmetric, and stop if not.
!! gmet(3,3)=reciprocal space metric (bohr**-2).
!! iout=if non-zero, output the new number of kpoints on unit iout
!! kbz(3,nkbz)= k vectors in the BZ.
!! nkbz = number of k-points whose weights are wtk
!! nsym=number of space group symmetries
!! symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!! timrev: if 1, the time reversal operation has to be taken into account
!!         if 0, no time reversal symmetry.
!! wtk(nkbz)=weight assigned to each k point.
!!
!! OUTPUT
!! ibz2bz(nkbz)=non-symmetrized indices of the k-points a.k.a. ibz2bz mapping
!!   The correspondence beween the iq_ibz point in IBZ and the iq_bz point in the full BZ is obtained via:
!!
!!       do ik_ibz=1,nkibz
!!         ik_bz = ibz2bz(ik_ibz)
!!       end do
!!
!! nkibz = number of k-points in the irreducible set
!! wtk_folded(nkbz)=weight assigned to each k point, taking into account the symmetries
!!
!! NOTES
!! The decomposition of the symmetry group in its primitives might speed up the execution.
!! The output variables are stored only in the range 1:nkbz
!!
!! PARENTS
!!      dfpt_looppert,elphon,ep_setupqpt,get_npert_rbz,getkgrid,harmonic_thermo
!!      m_bz_mesh,m_ifc,m_sigmaph
!!
!! CHILDREN
!!      sort_dp,wrtout
!!
!! SOURCE

subroutine symkpt(chksymbreak,gmet,ibz2bz,iout,kbz,nkbz,nkibz,nsym,&
& symrec,timrev,wtk,wtk_folded)

 use defs_basis
 use m_abicore
 use m_errors
 use m_sort

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symkpt'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: chksymbreak,iout,nkbz,nsym,timrev
 integer,intent(out) :: nkibz
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 integer,intent(inout) :: ibz2bz(nkbz) !vz_i
 real(dp),intent(in) :: gmet(3,3),kbz(3,nkbz),wtk(nkbz)
 real(dp),intent(out) :: wtk_folded(nkbz)

!Local variables -------------------------
!scalars
 integer :: identi,ii,ikpt,ikpt2,ind_ikpt,ind_ikpt2
 integer :: ikpt_current_length,isym,itim,jj,nkpout,quit,tident
 real(dp) :: difk,difk1,difk2,difk3,length2trial,reduce,reduce1,reduce2,reduce3
 character(len=500) :: message
!arrays
 integer,allocatable :: list(:)
 real(dp) :: gmetkpt(3),ksym(3)
 real(dp),allocatable :: length2(:)

! *********************************************************************

 if (timrev/=1 .and. timrev/=0) then
   write(message,'(a,i0)')' timrev should be 0 or 1, while it is equal to ',timrev
   MSG_BUG(message)
 end if

 if(nsym/=1)then
!  Find the identity symmetry operation
   do isym=1,nsym
     tident=1
     do jj=1,3
       if(symrec(jj,jj,isym)/=1)tident=0
       do ii=1,3
         if( ii/=jj .and. symrec(ii,jj,isym)/=0)tident=0
       end do
     end do
     if(tident==1)then
       identi=isym
       !write(message, '(a,i0)' )' symkpt : found identity, with number',identi
       !call wrtout(std_out,message,'COLL')
       exit
     end if
   end do
   if(tident==0)then
     MSG_BUG(' Did not find the identity operation')
   end if
 else
   identi=1
 end if

!Initialise the wtk_folded array using the wtk array :
 do ikpt=1,nkbz
   wtk_folded(ikpt)=wtk(ikpt)
 end do

!Here begins the serious business

!If there is some possibility for a change (otherwise, wtk_folded is
!correctly initialized to give no change) :
 if(nkbz/=1 .and. (nsym/=1 .or. timrev==1) )then

!  Store the length of vectors, but take into account umklapp
!  processes by selecting the smallest length of all symmetric vectors
   ABI_ALLOCATE(length2,(nkbz))


!MG FIXME:
! Here there's a possible problem with the order of symmetries because
! in listkk, time-reversal is the outermost loop. This can create inconsistencies in the symmetry tables.
   do ikpt=1,nkbz
     do isym=1,nsym
       do itim=1,(1-2*timrev),-2
!        Get the symmetric of the vector
         do ii=1,3
           ksym(ii)=itim*( kbz(1,ikpt)*symrec(ii,1,isym)&
&           +kbz(2,ikpt)*symrec(ii,2,isym)&
&           +kbz(3,ikpt)*symrec(ii,3,isym) )
           ksym(ii)=ksym(ii)-anint(ksym(ii)+tol8*half)
         end do
         gmetkpt(:)=gmet(:,1)*ksym(1)+gmet(:,2)*ksym(2)+gmet(:,3)*ksym(3)
         length2trial=ksym(1)*gmetkpt(1)+ksym(2)*gmetkpt(2)+ksym(3)*gmetkpt(3)
         if(isym==1 .and. itim==1)then
           length2(ikpt)=length2trial
         else
           if(length2(ikpt)>length2trial)length2(ikpt)=length2trial
         end if
       end do
     end do
   end do

!  Sort the lengths
   ABI_ALLOCATE(list,(nkbz))
   list(:)=(/ (ikpt,ikpt=1,nkbz) /)
   call sort_dp(nkbz,length2,list,tol14)
!  do ikpt=1,nkbz
!  write(std_out,*)ikpt,length2(ikpt),list(ikpt),kbz(1:3,list(ikpt))
!  end do

!  Examine whether the k point grid is symmetric or not
   if(chksymbreak==1)then
     ikpt_current_length=1
!    Loop on all k points
     do ikpt=1,nkbz
       ind_ikpt=list(ikpt)
!      Keep track of the current length, to avoid doing needless comparisons
       if(length2(ikpt)-length2(ikpt_current_length)>tol8)then
         ikpt_current_length=ikpt
       end if

       do isym=1,nsym
         do itim=1,(1-2*timrev),-2
           if(isym/=identi .or. itim/=1 )then

!            Get the symmetric of the vector
             do ii=1,3
               ksym(ii)=itim*( kbz(1,ind_ikpt)*symrec(ii,1,isym)&
&               +kbz(2,ind_ikpt)*symrec(ii,2,isym)&
&               +kbz(3,ind_ikpt)*symrec(ii,3,isym) )
             end do

!            Search over k-points with the same length, to find whether there is a connecting symmetry operation
             quit=0
             do ikpt2=ikpt_current_length,nkbz
!              The next line skip all ikpt2 vectors, as soon as one becomes larger than length2(ikpt)
!              Indeed, one is already supposed to have found a symmetric k point before this happens ...
               if(length2(ikpt2)-length2(ikpt)>tol8)exit
!              Ordered index
               ind_ikpt2=list(ikpt2)
               difk1= ksym(1)-kbz(1,ind_ikpt2)
               reduce1=difk1-anint(difk1)
               difk2= ksym(2)-kbz(2,ind_ikpt2)
               reduce2=difk2-anint(difk2)
               difk3= ksym(3)-kbz(3,ind_ikpt2)
               reduce3=difk3-anint(difk3)
               if(abs(reduce1)+abs(reduce2)+abs(reduce3)<tol8)then
!                The symmetric was found
                 quit=1
                 exit
               end if
             end do
             if(quit==0)then
               write(message,'(3a,i4,2a,9i3,2a,i6,1a,3es16.6,6a)' )&
&               'Chksymbreak=1 . It has been observed that the k point grid is not symmetric :',ch10,&
&               'for the symmetry number ',isym,ch10,&
&               'with symrec=',symrec(1:3,1:3,isym),ch10,&
&               'the symmetric of the k point number ',ind_ikpt2,' with components', kbz(1:3,ind_ikpt2),ch10,&
&               'does not belong to the k point grid.',ch10,&
&               'Read the description of the input variable chksymbreak,',ch10,&
&               'You might switch it to zero, or change your k point grid to one that is symmetric.'
               MSG_ERROR(message)
             end if

           end if ! End condition of non-identity symmetry
         end do ! itim
       end do ! isym

     end do ! ikpt
   end if ! chksymbreak==1


!  Eliminate the k points that are symmetric of another one
   do ikpt=1,nkbz-1

!    Ordered index
     ind_ikpt=list(ikpt)

!    Not worth to examine a k point that is a symmetric of another,
!    which is the case if its weight has been set to 0 by previous folding
     if(wtk_folded(ind_ikpt)<tol16)cycle

!    Loop on the remaining k-points
     do ikpt2=ikpt+1,nkbz

!      The next line eliminates pairs of vectors that differs by their length.
!      Moreover, since the list is ordered according to the length,
!      one can skip all other ikpt2 vectors, as soon as one becomes larger than length2(ikpt)
       if(length2(ikpt2)-length2(ikpt)>tol8)exit

!      Ordered index
       ind_ikpt2=list(ikpt2)

!      If the second vector is already empty, no interest to treat it
       if(wtk_folded(ind_ikpt2)<tol16)cycle

       quit=0
       do isym=1,nsym
         do itim=1,(1-2*timrev),-2
           if(isym/=identi .or. itim/=1 )then

!            Get the symmetric of the vector
             do ii=1,3
               ksym(ii)=itim*( kbz(1,ind_ikpt)*symrec(ii,1,isym)&
&               +kbz(2,ind_ikpt)*symrec(ii,2,isym)&
&               +kbz(3,ind_ikpt)*symrec(ii,3,isym) )
             end do

!            The do-loop was expanded to speed up the execution
             difk= ksym(1)-kbz(1,ind_ikpt2)
             reduce=difk-anint(difk)
             if(abs(reduce)>tol8)cycle
             difk= ksym(2)-kbz(2,ind_ikpt2)
             reduce=difk-anint(difk)
             if(abs(reduce)>tol8)cycle
             difk= ksym(3)-kbz(3,ind_ikpt2)
             reduce=difk-anint(difk)
             if(abs(reduce)>tol8)cycle

!            Here, have successfully found a symmetrical k-vector
!            Assign all the weight of the k-vector to its symmetrical
             wtk_folded(ind_ikpt)=wtk_folded(ind_ikpt)+wtk_folded(ind_ikpt2)
             wtk_folded(ind_ikpt2)=0._dp

!            Go to the next ikpt2 if the symmetric was found
             quit=1
             exit

           end if ! End condition of non-identity symmetry
         end do ! End loop on itim

         if(quit==1)exit
       end do !  End loop on isym
     end do ! End secondary loop over k-points
   end do ! End primary loop over k-points

   ABI_DEALLOCATE(length2)
   ABI_DEALLOCATE(list)
 end if ! End check on possibility of change

!Create the indexing array ibz2bz
 nkibz=0
 do ikpt=1,nkbz
   if(wtk_folded(ikpt)>tol8)then
     nkibz=nkibz+1
     ibz2bz(nkibz)=ikpt
   end if
 end do

 if(iout/=0)then
   if(nkbz/=nkibz)then
     write(message, '(a,a,a,i6,a)' )&
&     ' symkpt : the number of k-points, thanks to the symmetries,',ch10,&
&     ' is reduced to',nkibz,' .'
     call wrtout(iout,message,'COLL')
     if(iout/=std_out) then
       call wrtout(std_out,message,'COLL')
     end if

     nkpout=nkibz
     !if(nkibz>80)then
     !  write(message,'(a)' )' greater than 80, so only write 20 of them '
     !  call wrtout(std_out,message,'COLL')
     !  nkpout=20
     !end if
     !do ii=1,nkpout
     !  write(message, '(1x,i2,a2,3es16.8)' ) ii,') ',kbz(1:3,ibz2bz(ii))
     !  call wrtout(std_out,message,'COLL')
     !end do

!    DEBUG
!    write(message, '(a)' )'   Here are the new weights :'
!    call wrtout(std_out,message,'COLL')
!    do ikpt=1,nkbz,6
!    write(message, '(6f12.6)' ) wtk_folded(ikpt:min(nkbz,ikpt+5))
!    call wrtout(std_out,message,'COLL')
!    end do
!    ENDDEBUG
   else
     write(message, '(a)' )' symkpt : not enough symmetry to change the number of k points.'
     call wrtout(iout,message,'COLL')
     if(iout/=std_out) then
       call wrtout(std_out,message,'COLL')
     end if
   end if
 end if

!DEBUG
!write(std_out,*)' exit symkpt '
!write(std_out,*)' nkbz,nsym',nkbz,nsym
!write(std_out,*)' wtk',wtk
!write(std_out,*)' timrev',timrev
!if(timrev==0)stop
!if(option==1)stop
!ENDDEBUG

end subroutine symkpt
!!***

end module m_symkpt
!!***
