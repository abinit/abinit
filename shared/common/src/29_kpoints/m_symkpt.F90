!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_symkpt
!! NAME
!!  m_symkpt
!!
!! COPYRIGHT
!!  Copyright (C) 1999-2019 ABINIT group (XG,LSI,HM)
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

 use defs_basis
 use m_abicore
 use m_errors
 use m_sort
 use m_krank
 use m_numeric_tools
 use m_time

 implicit none

 private
!!***

 public :: symkpt
 public :: symkpt_new
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
!! comm=MPI communicator.
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
!! bz2ibz_smap(nkbz, 6)= Mapping BZ --> IBZ.
!!
!! NOTES
!! The decomposition of the symmetry group in its primitives might speed up the execution.
!! The output variables are stored only in the range 1:nkbz
!!
!! TODO
!!  Bad scaling wrt nkbz. Should try to MPI parallelize or implement more efficient algorithm
!!
!! PARENTS
!!      dfpt_looppert,elphon,ep_setupqpt,get_npert_rbz,getkgrid,harmonic_thermo
!!      m_bz_mesh,m_ifc,m_sigmaph
!!
!! CHILDREN
!!      sort_dp,wrtout
!!
!! SOURCE

subroutine symkpt(chksymbreak,gmet,ibz2bz,iout,kbz,nkbz,nkibz,nsym,symrec,timrev,wtk,wtk_folded, bz2ibz_smap, comm)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: chksymbreak,iout,nkbz,nsym,timrev,comm
 integer,intent(out) :: nkibz
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 integer,intent(inout) :: ibz2bz(nkbz) !vz_i
 real(dp),intent(in) :: gmet(3,3),kbz(3,nkbz),wtk(nkbz)
 real(dp),intent(out) :: wtk_folded(nkbz)
 integer,intent(out) :: bz2ibz_smap(6, nkbz)

!Local variables -------------------------
!scalars
 integer :: identi,ii,ikpt,ikpt2,ind_ikpt,ind_ikpt2,ierr
 integer :: ikpt_current_length,isym,itim,jj,nkpout,quit,tident
 real(dp) :: difk,difk1,difk2,difk3,length2trial,reduce,reduce1,reduce2,reduce3
 !real(dp) :: cpu,wall,gflops
 character(len=500) :: message
!arrays
 integer,allocatable :: list(:),bz2ibz_idx(:)
 real(dp) :: gmetkpt(3),ksym(3)
 real(dp),allocatable :: length2(:)

! *********************************************************************

 ABI_UNUSED((/comm/))

 if (timrev/=1 .and. timrev/=0) then
   write(message,'(a,i0)')' timrev should be 0 or 1, while it is equal to ',timrev
   MSG_BUG(message)
 end if

 ! Find the identity symmetry operation
 identi = 1
 tident = -1
 if (nsym/=1) then
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
   ABI_CHECK(tident == 1, 'Did not find the identity operation')
 end if

 ! Initialise the wtk_folded array using the wtk array
 do ikpt=1,nkbz
   wtk_folded(ikpt)=wtk(ikpt)
 end do

 ! Initialize bz2ibz_smap
 bz2ibz_smap = 0
 do ikpt=1,nkbz
   bz2ibz_smap(1, ikpt) = ikpt
   bz2ibz_smap(2, ikpt) = 1
 end do

 ! Here begins the serious business

 ! If there is some possibility for a change (otherwise, wtk_folded is correctly initialized to give no change)
 if(nkbz/=1 .and. (nsym/=1 .or. timrev==1) )then
   !call cwtime(cpu, wall, gflops, "start")

   ! Store the length of vectors, but take into account umklapp
   ! processes by selecting the smallest length of all symmetric vectors
   ABI_ALLOCATE(length2,(nkbz))

   do ikpt=1,nkbz
     do isym=1,nsym
       do itim=1,(1-2*timrev),-2
         ! Get the symmetric of the vector
         do ii=1,3
           ksym(ii)=itim*( kbz(1,ikpt)*symrec(ii,1,isym)&
            +kbz(2,ikpt)*symrec(ii,2,isym)&
            +kbz(3,ikpt)*symrec(ii,3,isym) )
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

   !call cwtime_report("symkpt: length", cpu, wall, gflops)

   ! Sort the lengths
   ABI_ALLOCATE(list,(nkbz))
   list(:)=(/ (ikpt,ikpt=1,nkbz) /)
   call sort_dp(nkbz,length2,list,tol14)
   ! do ikpt=1,nkbz; write(std_out,*)ikpt,length2(ikpt),list(ikpt),kbz(1:3,list(ikpt)); end do

   !call cwtime_report("symkpt: sort", cpu, wall, gflops)

   ! Examine whether the k point grid is symmetric or not
   ! This check scales badly with nkbz hence it's disabled for dense meshes.
   if (chksymbreak == 1 .and. nkbz < 40**3) then
     ikpt_current_length=1
     ! Loop on all k points
     do ikpt=1,nkbz
       ind_ikpt=list(ikpt)
       ! Keep track of the current length, to avoid doing needless comparisons
       if(length2(ikpt)-length2(ikpt_current_length)>tol8) ikpt_current_length=ikpt

       do isym=1,nsym
         do itim=1,(1-2*timrev),-2
           if(isym/=identi .or. itim/=1 )then
             ! Get the symmetric of the vector
             do ii=1,3
               ksym(ii)=itim*( kbz(1,ind_ikpt)*symrec(ii,1,isym)&
                 +kbz(2,ind_ikpt)*symrec(ii,2,isym)&
                 +kbz(3,ind_ikpt)*symrec(ii,3,isym) )
             end do

             ! Search over k-points with the same length, to find whether there is a connecting symmetry operation
             quit=0
             do ikpt2=ikpt_current_length,nkbz
               ! The next line skip all ikpt2 vectors, as soon as one becomes larger than length2(ikpt)
               ! Indeed, one is already supposed to have found a symmetric k point before this happens ...
               if(length2(ikpt2)-length2(ikpt)>tol8)exit
               ! Ordered index
               ind_ikpt2=list(ikpt2)
               difk1= ksym(1)-kbz(1,ind_ikpt2)
               reduce1=difk1-anint(difk1)
               difk2= ksym(2)-kbz(2,ind_ikpt2)
               reduce2=difk2-anint(difk2)
               difk3= ksym(3)-kbz(3,ind_ikpt2)
               reduce3=difk3-anint(difk3)
               if (abs(reduce1)+abs(reduce2)+abs(reduce3) < tol8) then
                 ! The symmetric was found
                 quit=1; exit
               end if
             end do

             if (quit==0) then
               write(message,'(3a,i4,2a,9i3,2a,i6,1a,3es16.6,6a)' )&
               'Chksymbreak=1. It has been observed that the k point grid is not symmetric:',ch10,&
               'for the symmetry number: ',isym,ch10,&
               'with symrec= ',symrec(1:3,1:3,isym),ch10,&
               'the symmetric of the k point number: ',ind_ikpt2,' with components: ', kbz(1:3,ind_ikpt2),ch10,&
               'does not belong to the k point grid.',ch10,&
               'Read the description of the input variable chksymbreak,',ch10,&
               'You might switch it to zero, or change your k point grid to one that is symmetric.'
               MSG_ERROR(message)
             end if

           end if ! End condition of non-identity symmetry
         end do ! itim
       end do ! isym

     end do ! ikpt
   end if ! chksymbreak==1

   ! Eliminate the k points that are symmetric of another one
   do ikpt=1,nkbz-1
     ! Ordered index
     ind_ikpt=list(ikpt)

     ! Not worth to examine a k point that is a symmetric of another,
     ! which is the case if its weight has been set to 0 by previous folding
     if (wtk_folded(ind_ikpt) < tol16) cycle

     ! Loop on the remaining k-points
     do ikpt2=ikpt+1,nkbz

       ! The next line eliminates pairs of vectors that differs by their length.
       ! Moreover, since the list is ordered according to the length,
       ! one can skip all other ikpt2 vectors, as soon as one becomes larger than length2(ikpt)
       if (length2(ikpt2) - length2(ikpt) > tol8) exit

       ! Ordered index
       ind_ikpt2=list(ikpt2)

       ! If the second vector is already empty, no interest to treat it
       if (wtk_folded(ind_ikpt2) < tol16) cycle

       quit = 0
       ! MG Dec 16 2018, Invert isym, itim loop to be consistent with listkk and GW routines
       ! Should always use this convention when applying symmetry operations in k-space.
       ! TODO: Postponed to v9 because it won't be possible to read old WFK files.
       do isym=1,nsym
         do itim=1,(1-2*timrev),-2
           if (isym/=identi .or. itim/=1) then
             ! Get the symmetric of the vector
             do ii=1,3
               ksym(ii)=itim*( kbz(1,ind_ikpt)*symrec(ii,1,isym) &
               +kbz(2,ind_ikpt)*symrec(ii,2,isym)&
               +kbz(3,ind_ikpt)*symrec(ii,3,isym) )
             end do

             ! The do-loop was expanded to speed up the execution
             difk= ksym(1)-kbz(1,ind_ikpt2)
             reduce=difk-anint(difk)
             if (abs(reduce)>tol8) cycle
             difk= ksym(2)-kbz(2,ind_ikpt2)
             reduce=difk-anint(difk)
             if (abs(reduce)>tol8) cycle
             difk= ksym(3)-kbz(3,ind_ikpt2)
             reduce=difk-anint(difk)
             if (abs(reduce)>tol8) cycle

             ! Here, have successfully found a symmetrical k-vector
             ! Assign all the weight of the k-vector to its symmetrical
             wtk_folded(ind_ikpt) = wtk_folded(ind_ikpt) + wtk_folded(ind_ikpt2)
             wtk_folded(ind_ikpt2) = zero

             ! Fill entries following listkk convention.
             ! Note however that here we always use symrec whereas listkk uses symrel^T by default
             ! so pay attention when using these tables to symmetrize wavefunctions.
             bz2ibz_smap(1, ind_ikpt2) = ind_ikpt
             bz2ibz_smap(2, ind_ikpt2) = isym
             ! Compute difference with respect to kpt2, modulo a lattice vector
             ! Sk1 + G0 = k2
             bz2ibz_smap(3:5, ind_ikpt2) = nint(-ksym(:) + kbz(:, ind_ikpt2) + tol12)
             ii = 0; if (itim == -1) ii = 1
             bz2ibz_smap(6, ind_ikpt2) = ii

             ! Go to the next ikpt2 if the symmetric was found
             quit = 1; exit
           end if ! End condition of non-identity symmetry
         end do ! isym
         if (quit == 1) exit
       end do ! itim

     end do ! ikpt2
   end do ! ikpt

   ABI_DEALLOCATE(length2)
   ABI_DEALLOCATE(list)
   !call cwtime_report("symkpt: loop", cpu, wall, gflops)
 end if ! End check on possibility of change

 ! Create the indexing array ibz2bz
 ABI_MALLOC(bz2ibz_idx, (nkbz))
 bz2ibz_idx = 0
 nkibz = 0
 do ikpt=1,nkbz
   if (wtk_folded(ikpt) > tol8) then
     nkibz = nkibz+1
     ibz2bz(nkibz) = ikpt
     bz2ibz_idx(ikpt) = nkibz
   end if
 end do

 ! bz2ibz_smap stores the index in the BZ. Here we replace the BZ index with the IBZ index.
 ierr = 0
 do ikpt=1,nkbz
   ind_ikpt = bz2ibz_idx(bz2ibz_smap(1, ikpt))
   if (ind_ikpt /= 0) then
     bz2ibz_smap(1, ikpt) = ind_ikpt
   else
     ierr = ierr + 1
   end if
 end do
 ABI_CHECK(ierr == 0, "Error while remapping bz2ibz_smap array")

 ABI_FREE(bz2ibz_idx)

 if(iout/=0)then
   if(nkbz/=nkibz)then
     write(message, '(a,a,a,i6,a)' )&
     ' symkpt : the number of k-points, thanks to the symmetries,',ch10,' is reduced to',nkibz,' .'
     call wrtout(iout,message,'COLL')
     if(iout/=std_out) call wrtout(std_out,message,'COLL')

     nkpout=nkibz
     !if (nkibz>80) then
     !  call wrtout(std_out,' greater than 80, so only write 20 of them ')
     !  nkpout=20
     !end if
     !do ii=1,nkpout
     !  write(message, '(1x,i2,a2,3es16.8)' ) ii,') ',kbz(1:3,ibz2bz(ii))
     !  call wrtout(std_out,message,'COLL')
     !end do

     !DEBUG
     !call wrtout(std_out,'   Here are the new weights :','COLL')
     !do ikpt=1,nkbz,6
     !  write(message, '(6f12.6)' ) wtk_folded(ikpt:min(nkbz,ikpt+5))
     !  call wrtout(std_out,message,'COLL')
     !end do
     !ENDDEBUG
   else
     write(message, '(a)' )' symkpt : not enough symmetry to change the number of k points.'
     call wrtout(iout,message,'COLL')
     if (iout/=std_out) call wrtout(std_out,message,'COLL')
   end if
 end if

end subroutine symkpt
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/symkpt_new
!! NAME
!! symkpt_new
!!
!! FUNCTION
!! Same routine as above but with an algorithm with better scalling than before.
!! From a few tests it produces the same IBZ as before but avoids computing the lengths and sorting.
!! Instead it uses the krank datatype to map k-points onto each other.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine symkpt_new(chksymbreak,gmet,ibz2bz,iout,kbz,nkbz,nkibz,nsym,symrec,timrev,bz2ibz_smap, comm)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: chksymbreak,iout,nkbz,nsym,timrev,comm
 integer,intent(out) :: nkibz
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 integer,intent(out) :: ibz2bz(nkbz)
 real(dp),intent(in) :: gmet(3,3),kbz(3,nkbz)
 integer,intent(out) :: bz2ibz_smap(6,nkbz)

!Local variables -------------------------
!scalars
 type(krank_t) :: krank
 integer :: identi,ii,ikpt,ikibz,ikpt_found
 integer :: isym,itim,jj,nkpout,tident
 !real(dp) :: cpu, gflops, wall
 character(len=500) :: message
!arrays
 real(dp) :: ksym(3),kpt1(3)

! *********************************************************************

 ABI_UNUSED(comm)
 ABI_UNUSED(gmet)

 if (timrev/=1 .and. timrev/=0) then
   write(message,'(a,i0)')' timrev should be 0 or 1, while it is equal to ',timrev
   MSG_BUG(message)
 end if

 ! Find the identity symmetry operation
 identi = 1
 tident = -1
 if (nsym/=1) then
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
   ABI_CHECK(tident == 1, 'Did not find the identity operation')
 end if

 ! Initialize
 ibz2bz = 0
 bz2ibz_smap = 0
 do ikpt=1,nkbz
   bz2ibz_smap(1, ikpt) = ikpt
   bz2ibz_smap(2, ikpt) = 1
   bz2ibz_smap(4, ikpt) = 1 ! We will use this as wtk_folded
 end do

 ! Start krank
 krank = krank_new(nkbz, kbz)

 ! Here begins the serious business
 !call cwtime(cpu, wall, gflops, "start")

 ! If there is some possibility for a change
 if(nkbz/=1 .and. (nsym/=1 .or. timrev==1) )then

   ! Examine whether the k point grid is symmetric or not
   ! This check scales badly with nkbz hence it's disabled for dense meshes.
   if (chksymbreak == 1 .and. nkbz < 40**3) then
     do ikpt=1,nkbz
       kpt1 = kbz(:,ikpt)

       do isym=1,nsym
         do itim=0,timrev
           ! Skip identity symmetry
           if (isym==identi .and. itim==0) cycle

           ! Get the symmetric of the vector
           do ii=1,3
             ksym(ii)=(1-2*itim)*( kpt1(1)*symrec(ii,1,isym)+&
                                   kpt1(2)*symrec(ii,2,isym)+&
                                   kpt1(3)*symrec(ii,3,isym) )
           end do

           !find this point
           ikpt_found = krank%get_index(ksym)
           !if (sum(abs(mod(ksym-kbz(:,ikpt_found),one)))>tol8) then
           !  MSG_ERROR('Wrong k-point mapping found by krank')
           !end if
           !if k-point not found
           if (ikpt_found < 0) then
             write(message,'(3a,i4,2a,9i3,2a,i6,1a,3es16.6,6a)' )&
             'Chksymbreak=1. It has been observed that the k point grid is not symmetric:',ch10,&
             'for the symmetry number: ',isym,ch10,&
             'with symrec= ',symrec(1:3,1:3,isym),ch10,&
             'the symmetric of the k point number: ',ikpt,' with components: ',kpt1(:),ch10,&
             'does not belong to the k point grid.',ch10,&
             'Read the description of the input variable chksymbreak,',ch10,&
             'You might switch it to zero, or change your k point grid to one that is symmetric.'
             MSG_ERROR(message)
           end if
         end do ! itim
       end do ! isym
     end do ! ikpt
   end if

   ! Here I generate the IBZ
   ikpt_loop: do ikpt=1,nkbz

     if (bz2ibz_smap(4, ikpt)==0) cycle
     kpt1 = kbz(:,ikpt)

     ! MG Dec 16 2018, Invert isym, itim loop to be consistent with listkk and GW routines
     ! Should always use this convention when applying symmetry operations in k-space.
     ! TODO: Postponed to v9 because it won't be possible to read old WFK files.
     do isym=1,nsym
       do itim=0,timrev
         ! Skip identity symmetry
         if (isym==identi .and. itim==0) cycle

         ! Get the symmetric of the vector
         do ii=1,3
           ksym(ii)=(1-2*itim)*( kpt1(1)*symrec(ii,1,isym)+&
                                 kpt1(2)*symrec(ii,2,isym)+&
                                 kpt1(3)*symrec(ii,3,isym) )
         end do

         !find this point
         ikpt_found = krank%get_index(ksym)
         !if k-point not found just cycle
         if (ikpt_found < 0) cycle
         !if (sum(abs(mod(ksym-kbz(:,ikpt_found),one)))>tol8) then
         !  MSG_ERROR('Wrong k-point mapping found by krank')
         !end if
         if (ikpt_found >= ikpt) cycle
         bz2ibz_smap(:3, ikpt)  = [ikpt_found,isym,itim]
         bz2ibz_smap(4,ikpt) = bz2ibz_smap(4,ikpt) + bz2ibz_smap(4,ikpt_found)
         bz2ibz_smap(4,ikpt_found) = 0
       end do ! itim
     end do ! isym
   end do ikpt_loop! ikpt

 end if ! End check on possibility of change
 !call cwtime_report(" ibz", cpu, wall, gflops)

 nkibz = 0
 do ikpt=1,nkbz
   ikibz = bz2ibz_smap(1,ikpt)
   if (ikibz /= ikpt) cycle
   nkibz = nkibz + 1
   ibz2bz(nkibz) = ikpt
 end do

 !do ikpt=1,nkbz
 !  write(*,*) ikpt, ibz2bz(ikpt), bz2ibz_smap(1,ikpt)
 !end do

 ! Initialize again
 bz2ibz_smap = 0

 ! Now I loop again over the points in the IBZ to find the mapping to the BZ
 do ikpt=1,nkibz
   ikibz = ibz2bz(ikpt)
   kpt1 = kbz(:,ikibz)

   ! HM: Here I invert the itim and isym loop to generate the same mapping as listkk
   do itim=0,timrev
     do isym=1,nsym
       ! Get the symmetric of the vector
       do ii=1,3
         ksym(ii)=(1-2*itim)*( kpt1(1)*symrec(ii,1,isym)+&
                               kpt1(2)*symrec(ii,2,isym)+&
                               kpt1(3)*symrec(ii,3,isym) )
       end do

       !find this point
       ikpt_found = krank%get_index(ksym)
       if (ikpt_found < 0) cycle
       if (bz2ibz_smap(1, ikpt_found) /= 0) cycle
       bz2ibz_smap(:3, ikpt_found) = [ikpt,isym,itim]
       bz2ibz_smap(4:, ikpt_found) = nint(kbz(:,ikpt_found)-ksym)
     end do
   end do
 end do
 !call cwtime_report(" map", cpu, wall, gflops)

 call krank%free()

 !Here I make a check if the mapping was sucessfull
 if (any(bz2ibz_smap(1, :) == 0)) then
   MSG_ERROR('Could not find mapping BZ to IBZ')
 end if

 !do ikpt=1,nkbz
 !  write(*,*) ikpt, ibz2bz(ikpt), bz2ibz_smap(1,ikpt)
 !end do

 if(iout/=0)then
   if(nkbz/=nkibz)then
     write(message, '(a,a,a,i6,a)' )&
     ' symkpt : the number of k-points, thanks to the symmetries,',ch10,' is reduced to',nkibz,' .'
     call wrtout(iout,message,'COLL')
     if(iout/=std_out) call wrtout(std_out,message,'COLL')

     nkpout=nkibz
   else
     write(message, '(a)' )' symkpt : not enough symmetry to change the number of k points.'
     call wrtout(iout,message,'COLL')
     if (iout/=std_out) call wrtout(std_out,message,'COLL')
   end if
 end if

end subroutine symkpt_new
!!***

end module m_symkpt
!!***
