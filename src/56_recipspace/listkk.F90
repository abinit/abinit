!{\src2tex{textfont=tt}}
!!****f* ABINIT/listkk
!! NAME
!! listkk
!!
!! FUNCTION
!! Given a list of nkpt1 initial k points kptns1 and a list of nkpt2
!! final k points kptns2, associates each final k pt with a "closest"
!! initial k point (or symmetric thereof, also taking possible umklapp)
!! as determined by a metric gmet, that commutes with the symmetry operations.
!! The algorithm does not scale as nkpt1 times nkpt2, thanks
!! to the ordering of the kptns1 and kptns2 vectors according to their
!! lengths, and comparison first between vectors of similar lengths.
!! Returns indirect indexing list indkk.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!  kptns1(3,nkpt1)=list of initial k points (reduced coordinates)
!!  kptns2(3,nkpt2)=list of final k points
!!  nkpt1=number of initial k points
!!  nkpt2=number of final k points
!!  nsym=number of symmetry elements in space group
!!  sppoldbl=if 1, no spin-polarisation doubling
!!           if 2, spin-polarisation doubling using symafm
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symmat(3,3,nsym)=symmetry operations (symrel or symrec, depending on
!!                   value of use_symrec
!!  timrev=1 if the use of time-reversal is allowed; 0 otherwise
!!  use_symrec :: if present and true, symmat assumed to be symrec, otherwise assumed to be symrel (default)
!!
!! OUTPUT
!!  dksqmax=maximal value of the norm**2 of the difference between
!!    a kpt2 vector and the closest k-point found from the kptns1 set, using symmetries.
!!  indkk(nkpt2*sppoldbl,6)=describe k point number of kpt1 that allows to
!!    generate wavefunctions closest to given kpt2
!!    if sppoldbl=2, use symafm to generate spin down wfs from spin up wfs
!!
!!    indkk(:,1)=k point number of kptns1
!!    indkk(:,2)=symmetry operation to be applied to kpt1, to give kpt1a
!!      (if 0, means no symmetry operation, equivalent to identity )
!!    indkk(:,3:5)=shift in reciprocal space to be given to kpt1a,
!!      to give kpt1b, that is the closest to kpt2.
!!    indkk(:,6)=1 if time-reversal was used to generate kpt1a from kpt1, 0 otherwise
!!
!! NOTES
!!  The tolerances tol12 and tol8 aims at giving a machine-independent ordering.
!!  (this trick is used in bonds.f, listkk.f, prtrhomxmn.f and rsiaf9.f)
!!  The tolerance tol12 is used for each component of the k vectors,
!!  and for the length of the vectors
!!  while the tolerance tol8 is used for the comparison of the squared lengths
!!  of the separate vectors.
!!
!! PARENTS
!!      initberry,inwffil,m_dvdb,m_ebands,m_fock,m_fstab,m_ifc,m_kpts,m_phgamma
!!      m_sigmaph,mlwfovlp_qp,printbxsf
!!
!! CHILDREN
!!      sort_dp,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine listkk(dksqmax,gmet,indkk,kptns1,kptns2,nkpt1,nkpt2,nsym,&
& sppoldbl,symafm,symmat,timrev,use_symrec)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_sort

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'listkk'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt1,nkpt2,nsym,sppoldbl,timrev
 real(dp),intent(out) :: dksqmax
 logical,optional,intent(in) :: use_symrec
!arrays
 integer,intent(in) :: symafm(nsym),symmat(3,3,nsym)
 integer,intent(out) :: indkk(nkpt2*sppoldbl,6)
 real(dp),intent(in) :: gmet(3,3),kptns1(3,nkpt1),kptns2(3,nkpt2)

!Local variables-------------------------------
!scalars
 integer :: l3,ig1,ig2,ig3,ii,ikpg1,ikpt1,ikpt2,ikpt2_done
 integer :: ilarger,ismaller,itrial
 integer :: isppol,isym,itimrev,jkpt1,jsym,jtime,limit
 integer :: nsym_used,timrev_used,usesym
 real(dp) :: dksq,dksqmn,lk2,llarger,ldiff,lsmaller,ltrial,min_l
 character(len=500) :: message
!arrays
 integer :: dkint(3),jdkint(3),k1int(3),k2int(3)
 integer, allocatable :: isort(:)
 real(dp) :: tsec(2)
 real(dp) :: dk(3),kpg1(3),kpt1a(3),k1(3),k2(3)
!real(dp) :: kasq,ka(3)
 real(dp),allocatable :: lkpg1(:),lkpg1_sorted(:)

! *************************************************************************

!write(std_out,*)' listkk : nkpt1,nkpt2,nsym=',nkpt1,nkpt2,nsym
 call timab(1021,1,tsec)

 if(sppoldbl<1 .or. sppoldbl>2)then
   write(message, '(a,i4,3a)' )&
&   'The value of sppoldbl is',sppoldbl,',',ch10,&
&   'but it should be either 1 or 2.'
   MSG_BUG(message)
 end if

!When usesym=0, the old way of converting the wavefunctions (without
!using the symmetries), is recovered.
 usesym=1

 nsym_used=nsym
 timrev_used=timrev
 if(usesym==0)nsym_used=1
 if(usesym==0)timrev_used=0

!Precompute the length of the kpt1 vectors, also taking into account
!possible umpklapp vectors
 limit=1 ; l3 = (2*limit+1)**3
 ABI_ALLOCATE(lkpg1,(l3*nkpt1))
 ABI_ALLOCATE(lkpg1_sorted,(l3*nkpt1))
 ABI_ALLOCATE(isort,(l3*nkpt1))
!write(std_out,*)' List of kpt1 vectors '
!write(std_out,*)' Length of the kpt1 vectors :'

 do ikpt1=1,nkpt1
   k1(:)=kptns1(:,ikpt1)
!  write(std_out,*)ikpt1,k1(:)
   k1int(:)=nint(k1(:)+tol12)
   k1(:)=k1(:)-k1int(:)
   do ig1=-limit,limit
     kpg1(1)=k1(1)+ig1
     do ig2=-limit,limit
       kpg1(2)=k1(2)+ig2
       do ig3=-limit,limit
         kpg1(3)=k1(3)+ig3

         ikpg1=ig1+limit+1 + (2*limit+1)*(ig2+limit) + (2*limit+1)**2*(ig3+limit) + l3*(ikpt1-1)
!        Compute the norm of the vector (also taking into account possible umklapp)
         lkpg1(ikpg1)=sqrt(gmet(1,1)*kpg1(1)**2+gmet(2,2)*kpg1(2)**2+&
&         gmet(3,3)*kpg1(3)**2+two*(gmet(2,1)*kpg1(2)*kpg1(1)+&
&         gmet(3,2)*kpg1(3)*kpg1(2)+gmet(3,1)*kpg1(3)*kpg1(1)))
         lkpg1_sorted(ikpg1)=lkpg1(ikpg1)
         isort(ikpg1)=ikpg1
!        write(std_out,*)' ikpt1,ig1,ig2,ig3,lkpg1=',ikpt1,ig1,ig2,ig3,lkpg1(ikpg1)
       end do
     end do
   end do
 end do

 call sort_dp( l3*nkpt1,lkpg1_sorted,isort,tol12)

!DEBUG
!write(std_out,*)' listkk : output list of kpt1 for checking purposes '
!write(std_out,*)' ii,ikpt1,isort(ii)-l3*(ikpt1-1),lkpg1_sorted(ii),lkpg1(isort(ii)) '
!do ii=1,l3*nkpt1
!ikpt1=(isort(ii)-1)/l3+1
!write(std_out,*)ii,ikpt1,isort(ii)-l3*(ikpt1-1),lkpg1_sorted(ii),lkpg1(isort(ii))
!enddo
!stop
!ENDDEBUG

 dksqmax=zero
 do isppol=1,sppoldbl
   do ikpt2=1,nkpt2

     ikpt2_done=0
!    Precompute the length of the kpt2 vector, with the Umklapp vector such that it is the closest to the Gamma point
     k2(:)=kptns2(:,ikpt2)
     k2int(:)=nint(k2(:)+tol12)
     k2(:)=k2(:)-k2int(:)
     lk2=sqrt(gmet(1,1)*k2(1)**2+gmet(2,2)*k2(2)**2+&
&     gmet(3,3)*k2(3)**2+two*(gmet(2,1)*k2(2)*k2(1)+&
&     gmet(3,2)*k2(3)*k2(2)+gmet(3,1)*k2(3)*k2(1)))

!    DEBUG
!    write(std_out, '(a,i4,7es16.6)' )' listkk : ikpt2,kptns2(:,ikpt2),k2(:),lk2=',ikpt2,kptns2(:,ikpt2),k2(:),lk2
!    if(ikpt2/=17)cycle
!    ENDDEBUG

!    Find the kpt1 vector whose length is the most similar to the length of lk2
!    up to a tolerance. Use a bissection algorithm.
     ismaller=0       ; lsmaller=zero
     ilarger=l3*nkpt1+1 ; llarger=huge(one)
!    This loop should never reach l3*nkpt1, since this is a bissection algorithm
     do ii=1,l3*nkpt1
       if((ilarger-ismaller)<2 .or. (llarger-lsmaller)<2*tol12)exit
       itrial=(ilarger+ismaller)/2 ; ltrial=lkpg1_sorted(itrial)
       if((ltrial-lk2)>tol12)then
         ilarger=itrial ; llarger=ltrial
       else if((ltrial-lk2)<-tol12)then
         ismaller=itrial ; lsmaller=ltrial
       else
         ismaller=itrial ; lsmaller=ltrial
         ilarger=itrial ; llarger=ltrial
       end if
     end do
     itrial=ismaller
     if(abs(llarger-lk2)<abs(lsmaller-lk2)-tol12)itrial=ilarger
     if(itrial==0)itrial=ilarger
     ismaller=itrial ; ilarger=itrial
!    write(std_out,*)' listkk : starting search at itrial=',itrial

     dksqmn=huge(one)

!    The ii index is dummy. This avoids an infinite loop.
     do ii=1,l3*nkpt1
!      do ikpt1=1,nkpt1

!      If the difference in length between the trial vector and the target vector is bigger
!      than the already achieved distance, the search is finished ...
       ldiff=abs(lkpg1_sorted(itrial)-lk2)


!      DEBUG
!      write(std_out,*)' listkk : ii,itrial,lkpg1_sorted(itrial),lk2,ldiff,dksqmn=',ii,itrial,lkpg1_sorted(itrial),lk2,ldiff,dksqmn
!      ENDDEBUG
       if(ldiff**2>dksqmn+tol8)exit

!      If this k-point has already been examined in a previous batch, skip it
!      First, compute the minimum of the difference of length of the sets of associated vectors thanks to Umklapp vectors
!      with the target vector
       ikpt1=(isort(itrial)-1)/l3+1
       min_l=minval(abs(lkpg1((ikpt1-1)*l3+1:(ikpt1-1)*l3+l3)-lk2))
!      Then compare with the current ldiff

!      DEBUG
!      write(std_out,*)' listkk : ikpt1,min_l,ldiff=',ikpt1,min_l,ldiff
!      ENDDEBUG

       if(min_l > ldiff-tol12)then

!        Now, will examine the trial vector, and the symmetric ones
!MG FIXME:
! Here there's a possible problem with the order of symmetries because
! in symkpt, time-reversal is the innermost loop. This can create inconsistencies in the symmetry tables.
         do itimrev=0,timrev_used
           do isym=1,nsym_used

!            Select magnetic characteristic of symmetries
             if(isppol==1 .and. symafm(isym)==-1)cycle
             if(isppol==2 .and. symafm(isym)==1)cycle

!            Compute symmetric point to kpt1
             if(usesym==1)then
!              original code only used transpose(symrel)
!              kpt1a(:)=symrel(1,:,isym)*kptns1(1,ikpt1)+&
!              &             symrel(2,:,isym)*kptns1(2,ikpt1)+&
!              &             symrel(3,:,isym)*kptns1(3,ikpt1)
               if (present(use_symrec)) then
                 if (use_symrec) then
                   kpt1a(:) = MATMUL(symmat(:,:,isym),kptns1(:,ikpt1))
                 else
                   kpt1a(:) = MATMUL(TRANSPOSE(symmat(:,:,isym)),kptns1(:,ikpt1))
                 end if
               else
                 kpt1a(:) = MATMUL(TRANSPOSE(symmat(:,:,isym)),kptns1(:,ikpt1))
               end if
               kpt1a(:)=(1-2*itimrev)*kpt1a(:)
             else
               kpt1a(:)=kptns1(:,ikpt1)
             end if

!            Compute difference with respect to kpt2, modulo a lattice vector
             dk(:)=kptns2(:,ikpt2)-kpt1a(:)
             if(usesym==1)then
!              The tolerance insure similar behaviour on different platforms
!              XG120418 : Actually, *assumes* that the closest point will have reduced
!              coordinates differing by less than 1/2 . There might be elongated
!              cells where this is not correct ...
               dkint(:)=nint(dk(:)+tol12)
               dk(:)=dk(:)-dkint(:)
             else
               dkint(:)=0
             end if

!            Compute norm of the difference vector, and update kpt1 if better.
             dksq=gmet(1,1)*dk(1)**2+gmet(2,2)*dk(2)**2+&
&             gmet(3,3)*dk(3)**2+two*(gmet(2,1)*dk(2)*dk(1)+&
&             gmet(3,2)*dk(3)*dk(2)+gmet(3,1)*dk(3)*dk(1))

             if (dksq<dksqmn+tol8) then

!              If exactly the right point (without using symmetries neither umklapp vector), will exit the search
!              Note that in this condition, each coordinate is tested separately, without squaring. So, it is a much stronger
!              condition than dksqmn<tol12
               if(sum(abs(kptns2(:,ikpt2)-kptns1(:,ikpt1)))<3*tol12)then
                 ikpt2_done=1
               end if

!              Update in three cases : either if succeeded to have exactly the vector, or the distance is better,
!              or the distance is only slightly worsened so select the lowest itimrev, isym or ikpt1, in order to respect previous ordering
               if(  ikpt2_done==1 .or. &
&               dksq+tol12<dksqmn .or. &
&               ( abs(dksq-dksqmn)<tol12 .and. &
&               ((itimrev<jtime) .or. &
&               (itimrev==jtime .and. isym<jsym) .or. &
&               (itimrev==jtime .and. isym==jsym .and. ikpt1<jkpt1))))then

                 dksqmn=dksq
                 jkpt1=ikpt1
                 jsym=isym
                 jtime=itimrev
                 jdkint(:)=dkint(:)

!                DEBUG
!                write(std_out,*)' ikpt1,ikpt2=',ikpt1,ikpt2
!                write(std_out,*)' timrev_used=',timrev_used
!                write(std_out,*)' Succeeded to lower dskmn,ikpt2_done=',dksqmn,ikpt2_done
!                write(std_out,*)' ikpt1,isym,dkint(:),itimrev=',ikpt1,isym,dkint(:),itimrev
!                ka(:)=kpt1a(:)+dkint(:)
!                write(std_out,*)'        k1=',kpt1a(:)
!                write(std_out,*)'     dkint=',dkint(:)
!                write(std_out,*)' Actual k1=',ka(:)
!                write(std_out,*)'        k2=',kptns2(:,ikpt2)
!                kasq=gmet(1,1)*ka(1)**2+gmet(2,2)*ka(2)**2+&
!                &                  gmet(3,3)*ka(3)**2+two*(gmet(2,1)*ka(2)*ka(1)+&
!                &                  gmet(3,2)*ka(3)*ka(2)+gmet(3,1)*ka(3)*ka(1))
!                write(std_out,*)' Actual k1sq=',kasq
!                ENDDEBUG
               end if

             end if
             if(ikpt2_done==1)exit

           end do ! isym
           if(ikpt2_done==1)exit

         end do ! itimrev
         if(ikpt2_done==1)exit

       end if

!      Update the interval that has been explored
       if(itrial<ismaller)ismaller=itrial
       if(itrial>ilarger)ilarger=itrial

!      Select the next index to be tried (preferably the smaller indices, but this is a bit arbitrary).

!      DEBUG
!      write(std_out,*)' before choosing the next index :'
!      write(std_out,*)' ismaller,itrial,ilarger=',ismaller,itrial,ilarger
!      write(std_out,*)' lkpg1_sorted(ismaller-1),lk2,lkpg1_sorted(ilarger+1)=',lkpg1_sorted(ismaller-1),lk2,lkpg1_sorted(ilarger+1)
!      ENDDEBUG
       if(ismaller>1 .and. ilarger<l3*nkpt1)then
         if(abs(lkpg1_sorted(ismaller-1)-lk2)<abs(lkpg1_sorted(ilarger+1)-lk2)+tol12)then
           itrial=ismaller-1
         else
           itrial=ilarger+1
         end if
       end if
       if(ismaller==1 .and. ilarger<l3*nkpt1)itrial=ilarger+1
       if(ismaller>1 .and. ilarger==l3*nkpt1)itrial=ismaller-1
!      if(ismaller==1 .and. ilarger==l3*nkpt1), we are done with the loop !

     end do ! ikpt1

     indkk(ikpt2+(isppol-1)*nkpt2,1)=jkpt1
     indkk(ikpt2+(isppol-1)*nkpt2,2)=jsym
     indkk(ikpt2+(isppol-1)*nkpt2,3:5)=jdkint(:)
     indkk(ikpt2+(isppol-1)*nkpt2,6)=jtime
     dksqmax=max(dksqmax,dksqmn)

     if(dksqmn<-tol12)then
       write(message, '(a,es16.6)' )'  The minimum square of dk has negative norm: dksqmn=',dksqmn
       MSG_BUG(message)
     end if

!    DEBUG
!    write(std_out,'(a,i6,i2,2x,i6,5i3,es24.14)' )' listkk: ikpt2,isppol,indkk(ikpt2+(isppol-1)*nkpt2,:)=',ikpt2,isppol,indkk(ikpt2+(isppol-1)*nkpt2,:),dksqmn
!    if(nkpt1==17)stop
!    ENDDEBUG

   end do ! ikpt2
 end do ! isppol

 ABI_DEALLOCATE(isort)
 ABI_DEALLOCATE(lkpg1)
 ABI_DEALLOCATE(lkpg1_sorted)

 call timab(1021,2,tsec)

end subroutine listkk
!!***
