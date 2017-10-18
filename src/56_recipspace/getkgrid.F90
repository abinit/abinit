!{\src2tex{textfont=tt}}
!!****f* ABINIT/getkgrid
!! NAME
!! getkgrid
!!
!! FUNCTION
!! Compute the grid of k points in the irreducible Brillouin zone.
!! Note that nkpt (and nkpthf) can be computed by calling this routine with nkpt=0, provided that kptopt/=0.
!! If downsampling is present, also compute a downsampled k grid.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR, MM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! chksymbreak= if 1, will check whether the k point grid is symmetric (for kptopt=1,2 and 4), and stop if not.
!! iout=unit number for echoed output . 0 if no output is wished.
!! iscf= ( <= 0 =>non-SCF), >0 => SCF)  MG: FIXME I don't understand why we have to pass the value iscf.
!! kptopt=option for the generation of k points
!!   (defines whether spatial symmetries and/or time-reversal can be used)
!! msym=default maximal number of symmetries
!! nsym=number of symmetries
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym)=symmetry operations in real space in terms of primitive translations
!! vacuum(3)=for each direction, 0 if no vacuum, 1 if vacuum
!! [downsampling(3) = input variable that governs the downsampling]
!!
!! OUTPUT
!! kptrlen=length of the smallest real space supercell vector associated with the lattice of k points.
!! nkpt_computed=number of k-points in the IBZ computed in the present routine
!! If nkpt/=0  the following are also output :
!!   kpt(3,nkpt)=reduced coordinates of k points.
!!   wtk(nkpt)=weight assigned to each k point.
!! [fullbz(3,nkpt_fullbz)]=k-points generated in the full Brillouin zone.
!!   In output: allocated array with the list of k-points in the BZ.
!! [kpthf(3,nkpthf)]=k-points generated in the full Brillouin zone, possibly downsampled (for Fock).
!!
!! NOTES
!!  msym not needed since nsym is the last index.
!!
!! SIDE EFFECTS
!! Input/Output
!! nkpt=number of k points (might be zero, see output description)
!! kptrlatt(3,3)=k-point lattice specification
!! nshiftk=actual number of k-point shifts in shiftk
!! shiftk(3,210)=shift vectors for k point generation
!! [nkpthf = number of k points in the full BZ, for the Fock operator]
!!
!! PARENTS
!!      ep_setupqpt,getshell,inkpts,inqpt,m_ab7_kpoints,m_bz_mesh,m_kpts
!!      nonlinear,testkgrid,thmeig
!!
!! CHILDREN
!!      mati3inv,matr3inv,metric,smallprim,smpbz,symkpt
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine getkgrid(chksymbreak,iout,iscf,kpt,kptopt,kptrlatt,kptrlen,&
& msym,nkpt,nkpt_computed,nshiftk,nsym,rprimd,shiftk,symafm,symrel,vacuum,wtk,&
& fullbz,nkpthf,kpthf,downsampling) ! optional

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getkgrid'
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_56_recipspace, except_this_one => getkgrid
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: chksymbreak,iout,iscf,kptopt,msym,nkpt,nsym
 integer,intent(inout),optional :: nkpthf
 integer,intent(inout) :: nshiftk
 integer,intent(inout) :: nkpt_computed !vz_i
 real(dp),intent(out) :: kptrlen
!arrays
 integer,intent(in) :: symafm(msym),symrel(3,3,msym),vacuum(3)
 integer,optional,intent(in) :: downsampling(3)
 integer,intent(inout) :: kptrlatt(3,3)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: shiftk(3,210)
 real(dp),intent(inout) :: kpt(3,nkpt) !vz_i
 real(dp),intent(inout) :: wtk(nkpt)
 real(dp),optional,allocatable,intent(out) :: fullbz(:,:)
 real(dp),optional,intent(out) :: kpthf(:,:)

!Local variables-------------------------------
!scalars
 integer, parameter :: max_number_of_prime=47,mshiftk=210
 integer :: brav,decreased,found,ii,ikpt,iprime,ishiftk,isym,jshiftk,kshiftk,mkpt,mult
 integer :: nkpt_fullbz,nkptlatt,nshiftk2,nsym_used,option
 integer :: test_prime,timrev
 real(dp) :: length2,ucvol,ucvol_super
 character(len=500) :: message
!arrays
 integer, parameter :: prime_factor(max_number_of_prime)=(/2,3,5,7,9, 11,13,17,19,23,&
&  29,31,37,41,43, 47,53,59,61,67,&
&  71,73,79,83,89, 97,101,103,107,109,&
&  113,127,131,137,139, 149,151,157,163,167,&
&  173,179,181,191,193, 197,199/)
 integer :: kptrlatt2(3,3)
 integer,allocatable :: belong_chain(:),generator(:),indkpt(:),number_in_chain(:)
 integer,allocatable :: repetition_factor(:),symrec(:,:,:)
! real(dp) :: cart(3,3)
 real(dp) :: dijk(3),delta_dmult(3),dmult(3),fact_vacuum(3),gmet(3,3)
 real(dp) :: gmet_super(3,3),gprimd(3,3),gprimd_super(3,3),klatt2(3,3)
 real(dp) :: klatt3(3,3),kptrlattr(3,3),ktransf(3,3),ktransf_invt(3,3)
 real(dp) :: metmin(3,3),minim(3,3),rmet(3,3),rmet_super(3,3),rprimd_super(3,3)
 real(dp),allocatable :: deltak(:,:),kpt_fullbz(:,:),shiftk2(:,:),shiftk3(:,:),spkpt(:,:),wtk_folded(:),wtk_fullbz(:)

! *************************************************************************

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 if (kptopt==1.or.kptopt==4) then
   ! Cannot use antiferromagnetic symmetry operations to decrease the number of k points
   nsym_used=0
   do isym=1,nsym
     if(symafm(isym)==1)nsym_used=nsym_used+1
   end do
   ABI_ALLOCATE(symrec,(3,3,nsym_used))
   nsym_used=0
   do isym=1,nsym ! Get the symmetry matrices in terms of reciprocal basis
     if(symafm(isym)==1)then
       nsym_used=nsym_used+1
       call mati3inv(symrel(:,:,isym),symrec(:,:,nsym_used))
     end if
   end do
 else if (kptopt==2) then
   !Use only the time-reversal
   nsym_used=1
   ABI_ALLOCATE(symrec,(3,3,1))
   symrec(1:3,1:3,1)=0
   do ii=1,3
     symrec(ii,ii,1)=1
   end do
 end if

 kptrlatt2(:,:)=kptrlatt(:,:)
 nshiftk2=nshiftk
 ABI_ALLOCATE(shiftk2,(3,mshiftk))
 ABI_ALLOCATE(shiftk3,(3,mshiftk))
 shiftk2(:,:)=shiftk(:,:)

!Find a primitive k point lattice, if possible, by decreasing the number of shifts.
 if(nshiftk2/=1)then

   do ! Loop to be repeated if there has been a successful reduction of nshiftk2

     ABI_ALLOCATE(deltak,(3,nshiftk2))
     ABI_ALLOCATE(repetition_factor,(nshiftk2))
     ABI_ALLOCATE(generator,(nshiftk2))
     ABI_ALLOCATE(belong_chain,(nshiftk2))
     ABI_ALLOCATE(number_in_chain,(nshiftk2))

     decreased=0
     deltak(1,1:nshiftk2)=shiftk2(1,1:nshiftk2)-shiftk2(1,1)
     deltak(2,1:nshiftk2)=shiftk2(2,1:nshiftk2)-shiftk2(2,1)
     deltak(3,1:nshiftk2)=shiftk2(3,1:nshiftk2)-shiftk2(3,1)
     deltak(:,:)=deltak(:,:)-floor(deltak(:,:)+tol8)

!    Identify for each shift, the smallest repetition prime factor that yields a reciprocal lattice vector.
     repetition_factor(:)=0
     repetition_factor(1)=1
     do ishiftk=2,nshiftk2
       do iprime=1,max_number_of_prime
         test_prime=prime_factor(iprime)
         dmult(:)=test_prime*deltak(:,ishiftk)
         if(sum(abs( dmult(:)-nint(dmult(:)) ))<tol8)then
           repetition_factor(ishiftk)=test_prime
           exit
         end if
       end do
     end do

!    Initialize the selection of tentative generators
     generator(:)=1
     do ishiftk=1,nshiftk2
       if(repetition_factor(ishiftk)==0 .or. repetition_factor(ishiftk)==1)generator(ishiftk)=0
     end do

!    Try different shifts as generators, by order of increasing repetition factor,
!    provided they are equal or bigger than 2
     do iprime=1,max_number_of_prime
       do ishiftk=2,nshiftk2
         ! Note that ishiftk=1 is never a generator. It is the reference starting point.
         if(generator(ishiftk)==1 .and. repetition_factor(ishiftk)==prime_factor(iprime))then
!          Test the generator : is it indeed closed ?
           if(prime_factor(iprime)/=2)then
             do mult=2,prime_factor(iprime)-1
               dmult(:)=mult*deltak(:,ishiftk)
               found=0
               do jshiftk=1,nshiftk2
                 delta_dmult(:)=deltak(:,jshiftk)-dmult(:)
                 if(sum(abs(delta_dmult(:)-nint(delta_dmult(:)) ))<tol8)then
                   found=1
                   exit
                 end if
               end do
               if(found==0)exit
             end do
             if(found==0)generator(ishiftk)=0
           end if
           if(generator(ishiftk)==0)cycle
         else
           cycle
         end if
!        Now, test whether all k points can be found in all possible chains
         belong_chain(:)=0
         do jshiftk=1,nshiftk2
!          Initialize a chain starting from a k point not yet in a chain
           if(belong_chain(jshiftk)==0)then
             number_in_chain(:)=0   ! Not a member of the chain (yet)
             number_in_chain(jshiftk)=1   ! The first point in chain
             do mult=1,prime_factor(iprime)-1
               dmult(:)=mult*deltak(:,ishiftk)
               found=0
               do kshiftk=jshiftk+1,nshiftk2
                 delta_dmult(:)=deltak(:,kshiftk)-deltak(:,jshiftk)-dmult(:)
                 if(sum(abs(delta_dmult(:)-nint(delta_dmult(:)) ))<tol8)then
                   found=1
                   number_in_chain(kshiftk)=mult+1
                   exit
                 end if
               end do
               if(found==0)then
                 generator(ishiftk)=0
                 exit
               end if
             end do
             if(generator(ishiftk)==1)then
!              Store the chain
               do kshiftk=1,nshiftk2
                 if(number_in_chain(kshiftk)/=0)belong_chain(kshiftk)=number_in_chain(kshiftk)
               end do
             else
               exit
             end if
           end if
         end do

         if(generator(ishiftk)==0)cycle

!        For the generator based on ishiftk, all the k points have been found to belong to one chain.
!        All the initializing k points in the different chains have belong_chain(:)=1 .
!        They must be kept, and the others thrown away.
         ktransf(:,:)=0.0_dp
         ktransf(1,1)=1.0_dp
         ktransf(2,2)=1.0_dp
         ktransf(3,3)=1.0_dp
!        Replace one of the unit vectors by the shift vector deltak(:,ishiftk).
!        However, must pay attention not to make linear combinations.
!        Also, choose positive sign for first-non-zero value.
         if(abs(deltak(1,ishiftk)-nint(deltak(1,ishiftk)))>tol8)then
           if(deltak(1,ishiftk)>0)ktransf(:,1)= deltak(:,ishiftk)
           if(deltak(1,ishiftk)<0)ktransf(:,1)=-deltak(:,ishiftk)
         else if(abs(deltak(2,ishiftk)-nint(deltak(2,ishiftk)))>tol8)then
           if(deltak(2,ishiftk)>0)ktransf(:,2)= deltak(:,ishiftk)
           if(deltak(2,ishiftk)<0)ktransf(:,2)=-deltak(:,ishiftk)
         else if(abs(deltak(3,ishiftk)-nint(deltak(3,ishiftk)))>tol8)then
           if(deltak(3,ishiftk)>0)ktransf(:,3)= deltak(:,ishiftk)
           if(deltak(3,ishiftk)<0)ktransf(:,3)=-deltak(:,ishiftk)
         end if
!        Copy the integers to real(dp)
         kptrlattr(:,:)=kptrlatt2(:,:)
!        Go to reciprocal space
         call matr3inv(kptrlattr,klatt2)
!        Make the transformation
         do ii=1,3
           klatt3(:,ii)=ktransf(1,ii)*klatt2(:,1)+ktransf(2,ii)*klatt2(:,2)+ktransf(3,ii)*klatt2(:,3)
         end do
!        Back to real space
         call matr3inv(klatt3,kptrlattr)
!        real(dp) to integer
         kptrlatt2(:,:)=nint(kptrlattr(:,:))
!        Prepare the transformation of the shifts
         call matr3inv(ktransf,ktransf_invt)
         decreased=1
         kshiftk=0
         do jshiftk=1,nshiftk2
           if(belong_chain(jshiftk)==1)then
             kshiftk=kshiftk+1
!            Place the shift with index jshiftk in place of the one in kshiftk,
!            also transform the shift from the old to the new coordinate system
             shiftk3(:,kshiftk)=ktransf_invt(1,:)*shiftk2(1,jshiftk)+&
&             ktransf_invt(2,:)*shiftk2(2,jshiftk)+&
&             ktransf_invt(3,:)*shiftk2(3,jshiftk)
           end if
         end do
         nshiftk2=nshiftk2/prime_factor(iprime)
         shiftk2(:,1:nshiftk2)=shiftk3(:,1:nshiftk2)-floor(shiftk3(:,1:nshiftk2)+tol8)
         if(kshiftk/=nshiftk2)then
           MSG_BUG('The search for a primitive k point lattice contains a bug.')
         end if

!        If this trial shift was successful, must exit the loop on trial ishiftk,
!        and reinitialize the global loop
         if(decreased==1)exit
       end do ! ishiftk
       if(decreased==1)exit
     end do ! iprime

     ABI_DEALLOCATE(belong_chain)
     ABI_DEALLOCATE(deltak)
     ABI_DEALLOCATE(number_in_chain)
     ABI_DEALLOCATE(repetition_factor)
     ABI_DEALLOCATE(generator)

     if(decreased==0 .or. nshiftk2==1)exit

   end do ! Infinite loop

 end if !  End nshiftk being 1 or larger

!Impose shiftk coordinates to be in [0,1[
 do ishiftk=1,nshiftk2
   do ii=1,3
     if(shiftk2(ii,ishiftk)>one-tol8) shiftk2(ii,ishiftk)=shiftk2(ii,ishiftk)-1.0_dp
     if(shiftk2(ii,ishiftk)<-tol8)    shiftk2(ii,ishiftk)=shiftk2(ii,ishiftk)+1.0_dp
   end do
 end do

!Compute the number of k points in the G-space unit cell
 nkptlatt=kptrlatt2(1,1)*kptrlatt2(2,2)*kptrlatt2(3,3) &
& +kptrlatt2(1,2)*kptrlatt2(2,3)*kptrlatt2(3,1) &
& +kptrlatt2(1,3)*kptrlatt2(2,1)*kptrlatt2(3,2) &
& -kptrlatt2(1,2)*kptrlatt2(2,1)*kptrlatt2(3,3) &
& -kptrlatt2(1,3)*kptrlatt2(2,2)*kptrlatt2(3,1) &
& -kptrlatt2(1,1)*kptrlatt2(2,3)*kptrlatt2(3,2)

!Check whether the number of k points is positive,
!otherwise, change the handedness of kptrlatt2
 if(nkptlatt<=0)then
!  write(std_out,*)' getkgrid : nkptlatt is negative !'
   kptrlatt2(:,3)=-kptrlatt2(:,3)
   nkptlatt=-nkptlatt
   do ishiftk=1,nshiftk2
     shiftk2(3,ishiftk)=-shiftk2(3,ishiftk)
   end do
 end if

!Determine the smallest supercell R-vector whose contribution
!is not taken correctly into account in the k point integration.
!Increase enormously the size of the cell when vacuum is present.
 fact_vacuum(:)=1
 if(vacuum(1)==1)fact_vacuum(1)=1000.0_dp
 if(vacuum(2)==1)fact_vacuum(2)=1000.0_dp
 if(vacuum(3)==1)fact_vacuum(3)=1000.0_dp
 do ii=1,3
   rprimd_super(:,ii)=fact_vacuum(1)*rprimd(:,1)*kptrlatt2(1,ii)+&
&   fact_vacuum(2)*rprimd(:,2)*kptrlatt2(2,ii)+&
&   fact_vacuum(3)*rprimd(:,3)*kptrlatt2(3,ii)
 end do

 call metric(gmet_super,gprimd_super,-1,rmet_super,rprimd_super,ucvol_super)
 call smallprim(metmin,minim,rprimd_super)
 length2=min(metmin(1,1),metmin(2,2),metmin(3,3))
 kptrlen=sqrt(length2)

 !write(message,'(a,es16.6)' )' getkgrid : length of smallest supercell vector (bohr)=',kptrlen
 !call wrtout(std_out,message,'COLL')
! If the number of shifts has been decreased, determine the set of kptrlatt2 vectors
! with minimal length (without using fact_vacuum)
! It is worth to determine the minimal set of vectors so that the kptrlatt that is output
! does not seem screwy, although correct but surprising.
 if(nshiftk/=nshiftk2)then
   do ii=1,3
     rprimd_super(:,ii)=rprimd(:,1)*kptrlatt2(1,ii)+rprimd(:,2)*kptrlatt2(2,ii)+rprimd(:,3)*kptrlatt2(3,ii)
   end do
   call metric(gmet_super,gprimd_super,-1,rmet_super,rprimd_super,ucvol_super)
!  Shift vectors in cartesian coordinates (reciprocal space)
   do ishiftk=1,nshiftk2
     shiftk3(:,ishiftk)=gprimd_super(:,1)*shiftk2(1,ishiftk)+&
&     gprimd_super(:,2)*shiftk2(2,ishiftk)+&
&     gprimd_super(:,3)*shiftk2(3,ishiftk)
   end do
   call smallprim(metmin,minim,rprimd_super)
   call metric(gmet_super,gprimd_super,-1,rmet_super,minim,ucvol_super)
!  This is the new kptrlatt2
   do ii=1,3
     dijk(:)=gprimd(1,:)*minim(1,ii)+&
&     gprimd(2,:)*minim(2,ii)+&
&     gprimd(3,:)*minim(3,ii)
     kptrlatt2(:,ii)=nint(dijk(:))
   end do
!  Shifts in the new set of kptrlatt vectors
   do ishiftk=1,nshiftk2
     shiftk2(:,ishiftk)=minim(1,:)*shiftk3(1,ishiftk)+&
&     minim(2,:)*shiftk3(2,ishiftk)+&
&     minim(3,:)*shiftk3(3,ishiftk)
   end do
 end if

!brav=1 is able to treat all bravais lattices.
 brav=1
 mkpt=nkptlatt*nshiftk2

 ABI_ALLOCATE(spkpt,(3,mkpt))
 option=0
 if(iout/=0)option=1

 if (PRESENT(downsampling))then
   call smpbz(brav,iout,kptrlatt2,mkpt,nkpthf,nshiftk2,option,shiftk2,spkpt,downsampling=downsampling)
   if (PRESENT(kpthf)) then ! Returns list of k-points in the Full BZ, possibly downsampled for Fock
     kpthf = spkpt(:,1:nkpthf)
   end if
 endif

 call smpbz(brav,iout,kptrlatt2,mkpt,nkpt_fullbz,nshiftk2,option,shiftk2,spkpt)

 if (PRESENT(fullbz)) then ! Returns list of k-points in the Full BZ.
   ABI_ALLOCATE(fullbz,(3,nkpt_fullbz))
   fullbz = spkpt(:,1:nkpt_fullbz)
 end if

 if(kptopt==1 .or. kptopt==2 .or. kptopt==4)then

   ABI_ALLOCATE(indkpt,(nkpt_fullbz))
   ABI_ALLOCATE(kpt_fullbz,(3,nkpt_fullbz))
   ABI_ALLOCATE(wtk_fullbz,(nkpt_fullbz))
   ABI_ALLOCATE(wtk_folded,(nkpt_fullbz))

   kpt_fullbz(:,:)=spkpt(:,1:nkpt_fullbz)
   wtk_fullbz(1:nkpt_fullbz)=1.0_dp/dble(nkpt_fullbz)

   timrev=1;if (kptopt==4) timrev=0

   call symkpt(chksymbreak,gmet,indkpt,iout,kpt_fullbz,nkpt_fullbz,&
&   nkpt_computed,nsym_used,symrec,timrev,wtk_fullbz,wtk_folded)

   ABI_DEALLOCATE(symrec)
   ABI_DEALLOCATE(wtk_fullbz)

 else if(kptopt==3)then
   nkpt_computed=nkpt_fullbz
 end if

!The number of k points has been computed from kptopt, kptrlatt, nshiftk, shiftk,
!and the eventual symmetries, it is presently called nkpt_computed.

!Check that the argument nkpt is coherent with nkpt_computed, if nkpt/=0.
 if(nkpt/=nkpt_computed .and. nkpt/=0)then
   write(message, '(a,i6,5a,i6,7a)') &
&   'The argument nkpt=',nkpt,', does not match',ch10,&
&   'the number of k points generated by kptopt, kptrlatt, shiftk,',ch10,&
&   'and the eventual symmetries, that is, nkpt=',nkpt_computed,'.',ch10,&
&   'However, note that it might be due to the user,',ch10,&
&   'if nkpt is explicitely defined in the input file.',ch10,&
&   'In this case, please check your input file.'
   MSG_BUG(message)
 end if

 if(kptopt==1 .or. kptopt==2 .or. kptopt==4)then

   if(nkpt/=0)then
     do ikpt=1,nkpt
       kpt(:,ikpt)=kpt_fullbz(:,indkpt(ikpt))
       if(iscf>=0 .or. iscf==-3 .or. iscf==-1.or.iscf==-2)wtk(ikpt)=wtk_folded(indkpt(ikpt))
     end do
   end if

   ABI_DEALLOCATE(indkpt)
   ABI_DEALLOCATE(kpt_fullbz)
   ABI_DEALLOCATE(spkpt)
   ABI_DEALLOCATE(wtk_folded)

 else if(kptopt==3)then

   if(nkpt/=0)then
     kpt(:,1:nkpt)=spkpt(:,1:nkpt)
     if(iscf>1 .or. iscf==-3 .or. iscf==-1.or.iscf==-2)wtk(1:nkpt)=1.0_dp/dble(nkpt)
   end if
   ABI_DEALLOCATE(spkpt)

 end if

 kptrlatt(:,:)=kptrlatt2(:,:)
 nshiftk=nshiftk2
 shiftk(:,1:nshiftk)=shiftk2(:,1:nshiftk)
 ABI_DEALLOCATE(shiftk2)
 ABI_DEALLOCATE(shiftk3)

end subroutine getkgrid
!!***
