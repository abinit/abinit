!{\src2tex{textfont=tt}}
!!****f* ABINIT/irrzg
!!
!! NAME
!! irrzg
!!
!! FUNCTION
!! Find the irreducible zone in reciprocal space under the
!! symmetry group with real space rotations in symrel(3,3,nsym).
!! The (integer) rotation matrices symrel(3,3,nsym) express the new
!! real space positions (e.g. rotated atom positions) in REDUCED
!! coordinates, i.e. in coordinates expressed as fractions of real space
!! primitive translations (atomic coordinates xred).  tnons(3,nsym) express
!! the associated nonsymmorphic translations, again in reduced coordinates.
!! Special data structure created in irrzon.
!! First half holds mapping from irr zone to full zone;
!! part of second half holds repetition number info.
!! work1 is a work array to keep track of grid points found so far.
!! In case nspden=2 and nsppol=1, one has to take care of antiferromagnetic
!! operations. The subgroup of non-magnetic operations is used
!! to generate irrzon(:,:,2) and phnons(:,:,2), while the
!! full group is used to generate irrzon(:,:,1) and phnons(:,:,1)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of symmetry elements in group
!!  n1,n2,n3=box dimensions of real space grid (or fft box)
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrel(3,3,nsym)=symmetry matrices in real space (integers)
!!  tnons(3,nsym)=reduced nonsymmorphic translations
!! (symrel and tnons are in terms of real space primitive translations)
!!
!! OUTPUT
!!  irrzon(n1*n2*n3,2+(nspden/4),(nspden/nsppol)-3*(nspden/4))=integer array which contains the locations of related
!!   grid points and the repetition number for each symmetry class.
!!  phnons(2,n1*n2*n3,(nspden/nsppol)-3*(nspden/4))=phases associated with nonsymmorphic translations
!!
!! PARENTS
!!      m_ab7_kpoints,setsym,wfd_mkrho
!!
!! CHILDREN
!!      sort_int,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine irrzg(irrzon,nspden,nsppol,nsym,n1,n2,n3,phnons,symafm,symrel,tnons)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_sort

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'irrzg'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1,n2,n3,nspden,nsppol,nsym
!arrays
 integer,intent(in) :: symafm(nsym),symrel(3,3,nsym)
 integer,intent(out) :: irrzon(n1*n2*n3,2,(nspden/nsppol)-3*(nspden/4))
 real(dp),intent(in) :: tnons(3,nsym)
 real(dp),intent(out) :: phnons(2,n1*n2*n3,(nspden/nsppol)-3*(nspden/4))

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,id1,id2,id3,ifft,imagn,ind1,ind2,ipt,irep,isym,izone
 integer :: izonemax,j1,j2,j3,jj,k1,k2,k3,l1,l2,l3,nfftot,npt,nsym_used
 integer :: nzone,setzer,sppoldbl
 real(dp) :: arg,ph1i,ph1r,ph2i,ph2r,tau1,tau2,tau3
 logical,parameter :: afm_noncoll=.true. ! TRUE if antiferro symmetries are used in non-collinear magnetism
 character(len=500) :: message
!arrays
 integer,allocatable :: class(:),iperm(:),symafm_used(:),symrel_used(:,:,:)
 integer,allocatable :: work1(:)
 real(dp),allocatable :: tnons_used(:,:),work2(:,:)

! *************************************************************************

 ABI_ALLOCATE(class,(nsym))
 ABI_ALLOCATE(iperm,(nsym))
 ABI_ALLOCATE(work1,(n1*n2*n3))
 ABI_ALLOCATE(work2,(2,n1*n2*n3))

 nfftot=n1*n2*n3

 id1=n1/2+2
 id2=n2/2+2
 id3=n3/2+2

 sppoldbl=nspden/nsppol;if (nspden==4) sppoldbl=1

 do imagn=1,sppoldbl

!  Treat in a similar way the case of the full group and the non-magnetic subgroup
   nsym_used=0
   do isym=1,nsym
     if( (imagn==1 .and. sppoldbl==2) .or. symafm(isym)==1 .or. &
&     ((nspden==4).and.afm_noncoll) )then
       nsym_used=nsym_used+1
     end if
   end do

   if(imagn==2 .and. nsym_used/=nsym/2)then
     write(message, '(a,a,a,a,a,i4,a,i0)' )&
&     '  The number of ferromagnetic symmetry operations must be',ch10,&
&     '  half the total number of operations, while it is observed that',ch10,&
&     '  nsym=',nsym,' and nsym_magn=',nsym_used
     MSG_BUG(message)
   end if

   ABI_ALLOCATE(symafm_used,(nsym_used))
   ABI_ALLOCATE(symrel_used,(3,3,nsym_used))
   ABI_ALLOCATE(tnons_used,(3,nsym_used))

   nsym_used=0
   do isym=1,nsym
     if( (imagn==1 .and. sppoldbl==2) .or. symafm(isym)==1 .or.  &
&     ((nspden==4).and.afm_noncoll) ) then
       nsym_used=nsym_used+1
       symrel_used(:,:,nsym_used)=symrel(:,:,isym)
       tnons_used(:,nsym_used)=tnons(:,isym)
       symafm_used(nsym_used)=symafm(isym)
     end if
   end do
   if ((nspden/=4).or.(.not.afm_noncoll)) symafm_used=1


!  Zero out work array--later on, a zero entry will mean that
!  a given grid point has not yet been assigned to an ibz point
   work1(1:nfftot)=0
   irrzon(:,2,imagn)=0

!  Initialize at G=0 (always in irreducible zone)
   nzone=1
   irrzon(1,1,imagn)=1
   irrzon(1,2,imagn)=nsym_used
!  Set phase exp(2*Pi*I*G dot tnons) for G=0
   phnons(1,1,imagn)=one
   phnons(2,1,imagn)=zero
   npt=1
!  setting work1(1)=1 indicates that first grid point (G=0) is
!  in the iz (irreducible zone)
   work1(1)=1

   ind1=0

!  Loop over reciprocal space grid points:
   do i3=1,n3
     do i2=1,n2
       do i1=1,n1

         ind1=ind1+1

!        Check to see whether present grid point is equivalent to
!        any previously identified ibz point--if not, a new ibz point
!        has been found

         if (work1(ind1)==0) then

!          A new class has been found.

!          Get location of G vector (grid point) centered at 0 0 0
           l3=i3-(i3/id3)*n3-1
           l2=i2-(i2/id2)*n2-1
           l1=i1-(i1/id1)*n1-1

           do isym=1,nsym_used

!            Get rotated G vector Gj for each symmetry element
!            -- here we use the TRANSPOSE of symrel_used; assuming symrel_used expresses
!            the rotation in real space, the transpose is then appropriate
!            for G space symmetrization (p. 1172d,e of notes, 2 June 1995).
             j1=symrel_used(1,1,isym)*l1+&
&             symrel_used(2,1,isym)*l2+symrel_used(3,1,isym)*l3
             j2=symrel_used(1,2,isym)*l1+&
&             symrel_used(2,2,isym)*l2+symrel_used(3,2,isym)*l3
             j3=symrel_used(1,3,isym)*l1+&
&             symrel_used(2,3,isym)*l2+symrel_used(3,3,isym)*l3

!            Map into [0,n-1] and then add 1 for array index in [1,n]
             k1=1+mod(n1+mod(j1,n1),n1)
             k2=1+mod(n2+mod(j2,n2),n2)
             k3=1+mod(n3+mod(j3,n3),n3)
!            k1=1+map(j1,n1)
!            k2=1+map(j2,n2)
!            k3=1+map(j3,n3)

!            Get linear index of rotated point Gj
             ind2=k1+n1*((k2-1)+n2*(k3-1))

!            Store info for new class:
             class(isym)=ind2
             iperm(isym)=isym

!            Setting work array element to 1 indicates grid point has been
!            identified with iz point
             work1(ind2)=1

!            End of loop on isym
           end do

!          Sort integers into ascending order in each class
!          (this lumps together G vectors with the same linear index, i.e.
!          groups together symmetries which land on the same Gj)
           call sort_int(nsym_used,class,iperm)

!          Check repetition factor (how many identical copies of Gj occur
!          from all symmetries applied to G)
           irep=0
           do isym=1,nsym_used
             if (class(isym)==class(1)) then
               irep=irep+1
             end if
           end do
           ipt=nsym_used/irep

!          Repetition factor must be divisor of nsym_used:
           if (nsym_used/=(ipt*irep)) then
             write(message, '(a,i5,a,i6,a,a,a,a,a,a)' )&
&             '  irep=',irep,' not a divisor of nsym_used=',nsym_used,ch10,&
&             ' This usually indicates that',&
&             ' the input symmetries do not form a group.',ch10,&
&             ' Action : check the input symmetries carefully do they',&
&             ' form a group ? If they do, there is a code bug.'
             MSG_ERROR(message)
           end if

!          Compute phases for any nonsymmorphic symmetries
!          exp(-2*Pi*I*G dot tau(j)) for each symmetry j with
!          (possibly zero) nonsymmorphic translation tau(j)
           do jj=1,nsym_used
!            First get nonsymmorphic translation and see if nonzero
!            (iperm grabs the symmetries in the new order after sorting)
             isym=iperm(jj)
             tau1=tnons_used(1,isym)
             tau2=tnons_used(2,isym)
             tau3=tnons_used(3,isym)
             if (abs(tau1)>tol12.or.abs(tau2)>tol12&
&             .or.abs(tau3)>tol12) then
!              compute exp(-2*Pi*I*G dot tau) using original G
               arg=two_pi*(dble(l1)*tau1+dble(l2)*tau2+dble(l3)*tau3)
               work2(1,jj)=cos(arg)
               work2(2,jj)=-sin(arg)
             else
               work2(1,jj)=one
               work2(2,jj)=zero
             end if
           end do

!          All phases arising from symmetries which map to the same
!          G vector must actually be the same because
!          rho(Strans*G)=exp(2*Pi*I*(G) dot tau_S) rho(G)
!          must be satisfied; if exp(2*Pi*I*(G) dot tau_S) can be different
!          for two different symmetries S which both take G to the same St*G,
!          then the related Fourier components rho(St*G) must VANISH.
!          Thus: set "phase" to ZERO here in that case.
!          The G mappings occur in sets of irep members; if irep=1 then
!          all the G are unique.
!          MT 090212:
!          In the case of antiferromagn. symetries, one can have
!          rho(Strans*G)= -exp(2*Pi*I*(G) dot tau_S) rho(G)
!          (look at the minus !)
!          A special treatment is then operated on phons.
!          The later must be consistent with the use of phnons array
!          in symrhg.F90 routine.
!          XG 001108 :
!          Note that there is a tolerance on the
!          accuracy of tnons, especially when they are found from
!          the symmetry finder (with xred that might be a bit inaccurate)
           if (irep > 1) then
             do jj=1,nsym_used,irep
               setzer=0
               ph1r=work2(1,jj);ph1i=work2(2,jj)
               do j1=jj,jj+irep-1
                 ph2r=work2(1,j1);ph2i=work2(2,j1)
                 if (((ph2r+ph1r)**2+(ph2i+ph1i)**2) <= tol14) then
                   if (setzer/=1) setzer=-1
                 else if (((ph2r-ph1r)**2+(ph2i-ph1i)**2) > tol14) then
                   setzer=1
                 end if
               end do
!              Setzer= 0: phnons are all equal
!              Setzer=-1: phnons are equal in absolute value
!              Setzer= 1: some phnons are different
               if (setzer/=0) then
                 if (setzer==-1) then
                   if (afm_noncoll.and.nspden==4) then
                     arg=symafm_used(iperm(jj))
                     if (all(symafm_used(iperm(jj:jj+irep-1))==arg)) then
                       setzer=1
                     else
                       do j1=jj,jj+irep-1
                         work2(:,j1)=work2(:,j1)*dble(symafm_used(iperm(j1)))
                       end do
                     end if
                   else
                     setzer=1
                   end if
                 end if
                 if (setzer==1) work2(:,jj:jj+irep-1)=zero
               end if
             end do
!            Compress data if irep>1:
             jj=0
             do isym=1,nsym_used,irep
               jj=jj+1
               class(jj)=class(isym)
               work2(1,jj)=work2(1,isym)
               work2(2,jj)=work2(2,isym)
             end do
           end if

!          Put new unique points into irrzon array:
           irrzon(1+npt:ipt+npt,1,imagn)=class(1:ipt)

!          Put repetition number into irrzon array:
           irrzon(1+nzone,2,imagn)=irep

!          DEBUG
!          write(std_out,'(a,6i7)' )' irrzg : izone,i1,i2,i3,imagn,irrzon(859,2,1)=',&
!          &      1+nzone,i1,i2,i3,imagn,irrzon(859,2,1)
!          ENDDEBUG

!          Put phases (or 0) in phnons array:
           phnons(:,1+npt:ipt+npt,imagn)=work2(:,1:ipt)

!          Update number of points in irrzon array:
!          (irep must divide evenly into nsym_used !)
           npt=npt+ipt

!          Update number of classes:
           nzone=nzone+1

         end if
!        
!        End of loop on reciprocal space points, with indices i1, i2, i3
       end do
     end do
   end do

   if (allocated(symafm_used))  then
     ABI_DEALLOCATE(symafm_used)
   end if
   if (allocated(symrel_used))  then
     ABI_DEALLOCATE(symrel_used)
   end if
   if (allocated(tnons_used))  then
     ABI_DEALLOCATE(tnons_used)
   end if

 end do ! imagn

!Make sure number of real space points accounted for equals actual number of grid points
 if (npt/=n1*n2*n3) then
   write(message, '(a,a,a,a,i10,a,i10,a,a,a,a,a,a,a,a,a)' ) ch10,&
&   ' irrzg : ERROR -',ch10,&
&   '  npt=',npt,' and n1*n2*n3=',n1*n2*n3,' are not equal',ch10,&
&   '  This says that the total of all points in the irreducible',&
&   '  sector in real space',ch10,&
&   '  and all symmetrically equivalent',&
&   '  points, npt, does not equal the actual number',ch10,&
&   '  of real space grid points.'
   call wrtout(std_out,message,'COLL')
   write(message,'(3a)') &
&   ' This may mean that the input symmetries do not form a group',ch10,&
&   ' Action : check input symmetries carefully for errors.'
   MSG_ERROR(message)
 end if

!Perform some checks
 do imagn=1,sppoldbl

   do ifft=1,nfftot
     if (irrzon(ifft,1,imagn)<1.or.irrzon(ifft,1,imagn)>nfftot) then
       write(message,'(a,4i0,a,a)')&
&       '  ifft,irrzon(ifft,1,imagn),nfftot,imagn=',ifft,irrzon(ifft,1,imagn),nfftot,imagn,ch10,&
&       '  =>irrzon goes outside acceptable bounds.'
       MSG_BUG(message)
     end if
   end do

   izonemax=0
   do izone=1,nfftot
!    Global bounds
     if (irrzon(izone,2,imagn)<0.or.irrzon(izone,2,imagn)>(nsym/imagn)) then
       write(message, '(a,5i7,a,a)' )&
&       ' izone,nzone,irrzon(izone,2,imagn),nsym,imagn =',izone,nzone,irrzon(izone,2,imagn),nsym,imagn,ch10,&
&       '  =>irrzon goes outside acceptable bounds.'
       MSG_BUG(message)
     end if
!    Second index only goes up to nzone
     if(izonemax==0)then
       if (irrzon(izone,2,imagn)==0)izonemax=izone-1
     end if
     if(izonemax/=0)then
       if (irrzon(izone,2,imagn)/=0) then
         message = ' beyond izonemax, irrzon(izone,2,imagn) should be zero'
         MSG_BUG(message)
       end if
     end if
   end do

 end do ! imagn

 ABI_DEALLOCATE(class)
 ABI_DEALLOCATE(iperm)
 ABI_DEALLOCATE(work1)
 ABI_DEALLOCATE(work2)

end subroutine irrzg
!!***
