!!****m* ABINIT/m_symfind
!! NAME
!!  m_symfind
!!
!! FUNCTION
!!  Symmetry finder high-level API.
!!
!! COPYRIGHT
!!  Copyright (C) 2000-2020 ABINIT group (XG, RC)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_symfind

 use defs_basis
 use m_errors
 use m_abicore
 use m_symlist

 use m_symtk,     only : chkgrp, chkprimit, matr3inv, symrelrot, symdet, symcharac, holocell, smallprim, print_symmetries
 use m_geometry,  only : acrossb, xred2xcart
 use m_spgdata,   only : getptgroupma, symptgroup, spgdata

 implicit none

 private
!!***

 public :: symfind     ! From the symmetries of the Bravais lattice,
                       ! select those that leave invariant the system, and generate tnons
 public :: symanal     ! Find the space group from the list of symmetries and lattice parameters
 public :: symbrav     ! Determine the Bravais information from the list of symmetry operations, and the lattice vectors.
 public :: symlatt     ! Find the Bravais lattice and its symmetry operations (ptsymrel).
                       ! From the unit cell vectors (rprimd) and the corresponding metric tensor.

contains
!!***

!!****f* m_symfind/symfind
!! NAME
!! symfind
!!
!! FUNCTION
!! Symmetry finder.
!! From the symmetries of the Bravais lattice (ptsymrel),
!! select those that leave invariant the system, and generate
!! the corresponding tnons vectors.
!! The algorithm is explained in T.G. Worlton and J.L. Warren, Comp. Phys. Comm. 3, 88 (1972) [[cite:Worton1972]]
!!
!! INPUTS
!! berryopt    =  4/14, 6/16, 7/17: electric or displacement field
!! chrgat(natom) (optional)=target charge for each atom. Not always used, it depends on the value of constraint_kind
!! efield=cartesian coordinates of the electric field
!! gprimd(3,3)=dimensional primitive translations for reciprocal space
!! msym=default maximal number of symmetries
!! natom=number of atoms in cell.
!! noncoll=1 if non-collinear magnetism is activated
!          (3 components of spinat are taken into account)
!!         else 0
!! nptsym=number of point symmetries of the Bravais lattice
!! nucdipmom(3,natom) (optional) array of nuclear dipole moments
!! nzchempot=if non-zero, means that a z-spatially varying chemical potential is added
!! ptsymrel(3,3,1:msym)= nptsym point-symmetry operations
!!   of the Bravais lattice in real space in terms
!!   of primitive translations.
!! spinat(3,natom)=initial spin of each atom, in unit of hbar/2.
!! tolsym=tolerance for the symmetries
!! typat(natom)=integer identifying type of atom.
!! use_inversion=1 if inversion and improper rotations can be included in set of symmetries
!! xred(3,natom)=reduced coordinates of atoms in terms of real space
!!   primitive translations
!!
!! OUTPUT
!! ierr (optional)=if non-zero, the symmetry operations do not form a group
!! nsym=actual number of symmetries
!! symafm(1:msym)=(anti)ferromagnetic part of nsym symmetry operations
!! symrel(3,3,1:msym)= nsym symmetry operations in real space in terms
!!  of primitive translations
!! tnons(3,1:msym)=nonsymmorphic translations for each symmetry (would
!!  be 0 0 0 each for a symmorphic space group)
!!
!! PARENTS
!!      m_ab7_symmetry,m_effective_potential_file,m_ingeo,m_inkpts
!!      m_mover_effpot,m_tdep_sym,m_thmeig,m_use_ga
!!
!! CHILDREN
!!      holocell,matr3inv,smallprim,symrelrot,wrtout
!!
!! SOURCE

 subroutine symfind(berryopt,efield,gprimd,jellslab,msym,natom,noncoll,nptsym,nsym,&
&  nzchempot,prtvol, ptsymrel,spinat,symafm,symrel,tnons,tolsym,typat,use_inversion,xred,&
&  chrgat,ierr,nucdipmom)  ! Optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: berryopt,jellslab,msym,natom,noncoll,nptsym,nzchempot,use_inversion
 integer,intent(in) :: prtvol
 integer,optional,intent(out) :: ierr
 integer,intent(out) :: nsym
 real(dp),intent(in) :: tolsym
!arrays
 integer,intent(in) :: ptsymrel(3,3,msym),typat(natom)
 integer,intent(inout) :: symafm(msym),symrel(3,3,msym) !vz_i
 real(dp),intent(in) :: efield(3),gprimd(3,3),spinat(3,natom),xred(3,natom)
 real(dp),optional,intent(in) :: chrgat(natom)
 real(dp),optional, intent(in) :: nucdipmom(3,natom)
 real(dp),intent(inout) :: tnons(3,msym) !vz_i

!Local variables-------------------------------
!scalars
 integer :: found3,foundcl,iatom,iatom0,iatom1,iatom2,iatom3,iclass,iclass0,ierr_,ii
 integer :: isym,jj,kk,natom0,nclass,ntrial,printed,trialafm,trialok
 real(dp) :: det,ndnorm,nucdipmomcl2,nucdipmomcl20
 real(dp) :: spinat2,spinatcl2,spinatcl20
! TRUE if antiferro symmetries are used with non-collinear magnetism.
 integer :: afm_noncoll=1 !For noncoll==1.  If 1, all symops are permitted ; if 0 symafm must be 1.
!For noncoll=1. If noncoll_orthorhombic1, require the symmetry operations to be a subset of the orthorhombic symmetries, except if all spinat=0..
 integer :: noncoll_orthorhombic=0
 logical :: test_sameabsspin,test_samechrg
 logical :: test_samenucdipmom
 character(len=500) :: message
!arrays
 integer,allocatable :: class(:,:),natomcl(:),typecl(:)
 real(dp) :: diff(3),efieldrot(3),hand2(3),hand3(3),ndtest(3),rprimd(3,3),spinat0(3),xred0(3)
 !real(dp) :: symnucdipmom2(3)
 real(dp) :: symnucdipmom2cart(3,3),symnucdipmom2red(3,3)
 real(dp) :: symspinat1(3),symspinat2(3),symxred2(3),trialnons(3)
 real(dp),allocatable :: chrgat_(:)
 real(dp),allocatable :: chrgatcl(:)
 real(dp),allocatable :: local_nucdipmom(:,:,:),nucdipmomcl(:,:),nucdipmomred(:,:,:)
 real(dp),allocatable :: spinatcl(:,:),spinatred(:,:)

!**************************************************************************

!DEBUG
 if (prtvol>1) message="remove me later"
! write(std_out,*)' symfind : enter'
! call flush(6)
! write(std_out,*)' symfind : nzchempot= ',nzchempot
! write(std_out,*)'   ptsymrel matrices are :'
! do isym=1,nptsym
! write(std_out,'(i4,4x,9i4)' )isym,ptsymrel(:,:,isym)
! end do
! write(std_out,*)' symfind : natom=',natom
! do iatom=1,natom
! write(std_out,*)'  atom number',iatom
! write(std_out,*)'   typat   =',typat(iatom)
! write(std_out,*)'   spinat  =',spinat(:,iatom)
! write(std_out,*)'   xred    =',xred(:,iatom)
! end do
! write(std_out,*)' '
! call flush(6)
!ENDDEBUG

!Find the number of classes of atoms (type, chrg and spinat must be identical,
!spinat might differ by a sign, if aligned with the z direction, or,
! type and nucdipmom must be identical)
!natomcl(iclass) will contain the number of atoms in the class
!typecl(iclass) will contain the type of the atoms in the class
!chrgcl(iclass) will contain the charge of the atoms in the class
!spinatcl(1:3,iclass) will contain the spinat of the atoms in the class
!class(1:natomclass(iclass),iclass) will contain the index of the
!atoms belonging to the class
 ABI_ALLOCATE(class,(natom+3,natom))
 ABI_ALLOCATE(natomcl,(natom))
 ABI_ALLOCATE(typecl,(natom))
 ABI_ALLOCATE(chrgat_,(natom))
 ABI_ALLOCATE(chrgatcl,(natom))
 ABI_ALLOCATE(spinatcl,(3,natom))
 ABI_ALLOCATE(local_nucdipmom,(3,3,natom))
 ABI_ALLOCATE(nucdipmomcl,(3,natom))

 chrgat_(:)=zero
 if(present(chrgat))then
   chrgat_(:)=chrgat(:)
 endif

 local_nucdipmom(:,:,:) = zero
 if(present(nucdipmom)) then
    local_nucdipmom(1:3,1,:) = nucdipmom(1:3,:)
 end if
 ! for each nuclear dipole we need a local right handed coord system, so we can
 ! test later for whether a symmetry operation preserves the circulation induced
 ! by the dipole
 do iatom=1, natom
    ndnorm=sqrt(DOT_PRODUCT(local_nucdipmom(1:3,1,iatom),local_nucdipmom(1:3,1,iatom)))

    ! if nuclear dipole has effectively zero size, move on to the next atom
    if (ndnorm < tol8) cycle

    ! for testing purposes, we care only about direction so renormalize to unity
    local_nucdipmom(1:3,1,iatom) = local_nucdipmom(1:3,1,iatom)/ndnorm

    ! make a random vector, each component is (0,1]
    call random_number(ndtest)

    ! vector 2 is constructed to be orthogonal to original nuclear dipole moment vector
    call acrossb(local_nucdipmom(1:3,1,iatom),ndtest(1:3),local_nucdipmom(1:3,2,iatom))

    ! vector 3 is orthogonal to 1 and 2, and 1,2,3 form a right-handed set
    call acrossb(local_nucdipmom(1:3,1,iatom),local_nucdipmom(1:3,2,iatom),local_nucdipmom(1:3,3,iatom))
 end do

 ! need rprimd later to transform back to cart coords
 call matr3inv(gprimd,rprimd)

!Initialise with the first atom
 nclass=1
 natomcl(1)=1
 typecl(1)=typat(1)
 chrgatcl(1)=chrgat_(1)
 spinatcl(:,1)=spinat(:,1)
 nucdipmomcl(:,1)=local_nucdipmom(:,1,1)
 class(1,1)=1
 if(natom>1)then
   do iatom=2,natom
!    DEBUG
!    write(std_out,*)' '
!    write(std_out,*)' symfind : examine iatom=',iatom
!    ENDDEBUG
     foundcl=0
     do iclass=1,nclass
!      Compare the typat, chrg and spinat of atom iatom with existing ones.
!      At this stage, admit either identical spinat, or spin-flip spinat.
       if( typat(iatom)==typecl(iclass)) then
         test_samechrg= (abs(chrgat_(iatom)-chrgatcl(iclass))<tolsym)
         if(noncoll==0)then
           test_sameabsspin=(abs(abs(spinat(3,iatom))-abs(spinatcl(3,iclass)))<tolsym)
         else if(noncoll==1)then
           spinat2  =spinat(1,iatom)**2+spinat(2,iatom)**2+spinat(3,iatom)**2
           spinatcl2=spinatcl(1,iclass)**2+spinatcl(2,iclass)**2+spinatcl(3,iclass)**2
           test_sameabsspin=abs(spinat2-spinatcl2)<tolsym
         endif
         test_samenucdipmom= &
&             abs(local_nucdipmom(1,1,iatom)-nucdipmomcl(1,iclass))<tolsym .and. &
&             abs(local_nucdipmom(2,1,iatom)-nucdipmomcl(2,iclass))<tolsym .and. &
&             abs(local_nucdipmom(3,1,iatom)-nucdipmomcl(3,iclass))<tolsym
         ! note in the following test, m_chkinp/chkinp has already prevented nucdipmom to be
         ! nonzero when spinat is nonzero
         if( test_samechrg .and. test_sameabsspin .and. test_samenucdipmom ) then
!          DEBUG
!          write(std_out,*)' symfind : find it belongs to class iclass=',iclass
!          write(std_out,*)' symfind : spinat(:,iatom)=',spinat(:,iatom)
!          write(std_out,*)' symfind : spinatcl(:,iclass)=',spinatcl(:,iclass)
!          write(std_out,*)' symfind : test_sameabsspin=',test_sameabsspin
!          write(std_out,*)' '
!          ENDDEBUG
           natomcl(iclass)=natomcl(iclass)+1
           class(natomcl(iclass),iclass)=iatom
           foundcl=1
           exit
         end if
       end if
     end do
!    If no class with these characteristics exist, create one
     if(foundcl==0)then
       nclass=nclass+1
       natomcl(nclass)=1
       typecl(nclass)=typat(iatom)
       chrgatcl(nclass)=chrgat_(iatom)
       spinatcl(:,nclass)=spinat(:,iatom)
       nucdipmomcl(:,nclass)=local_nucdipmom(:,1,iatom)
       class(1,nclass)=iatom
     end if
   end do
 end if

!DEBUG
!write(std_out,*)' '
!write(std_out,*)' symfind : found ',nclass,' nclass of atoms'
!do iclass=1,nclass
!write(std_out,*)'  class number',iclass
!write(std_out,*)'   natomcl =',natomcl(iclass)
!write(std_out,*)'   typecl  =',typecl(iclass)
!write(std_out,*)'   spinatcl=',spinatcl(:,iclass)
!write(std_out,*)'   class   =',(class(iatom,iclass),iatom=1,natomcl(iclass))
!end do
!write(std_out,*)' '
!ENDDEBUG

!Select the class with the least number of atoms, and non-zero spinat if any
!It is important to select a magnetic class of atom, if any, otherwise
!the determination of the initial (inclusive) set of symmetries takes only
!non-magnetic symmetries, and not both magnetic and non-magnetic ones, see later.
!On the contrary, the chrgat_ data does not play any role, it is invariant upon atomic-centered  symmetries
 iclass0=1
 natom0=natomcl(1)
 spinatcl20=spinatcl(1,1)**2+spinatcl(2,1)**2+spinatcl(3,1)**2
 nucdipmomcl20=nucdipmomcl(1,1)**2+nucdipmomcl(2,1)**2+nucdipmomcl(3,1)**2
 if(nclass>1)then
   do iclass=2,nclass
     spinatcl2=spinatcl(1,iclass)**2+spinatcl(2,iclass)**2+spinatcl(3,iclass)**2
     nucdipmomcl2=nucdipmomcl(1,iclass)**2+nucdipmomcl(2,iclass)**2+nucdipmomcl(3,iclass)**2
     if( (natomcl(iclass)<natom0 &
&          .and. .not. (spinatcl20>tolsym .and. spinatcl2<tolsym) &
&          .and. .not. (nucdipmomcl20>tolsym .and. nucdipmomcl2<tolsym) )  &
&     .or. (spinatcl20<tolsym .and. spinatcl2>tolsym) &
&     .or. (nucdipmomcl20<tolsym .and. nucdipmomcl2>tolsym)                        )then
       iclass0=iclass
       natom0=natomcl(iclass)
       spinatcl20=spinatcl2
       nucdipmomcl20=nucdipmomcl2
     end if
   end do
 end if

 printed=0

!If non-collinear spinat have to be used, transfer them in reduced coordinates
 if (noncoll==1) then
   ABI_ALLOCATE(spinatred,(3,natom))
   do iatom=1,natom
     do ii=1,3
       spinatred(1:3,iatom)=MATMUL(TRANSPOSE(gprimd),spinat(1:3,iatom))
     end do
   end do
 end if

!DEBUG
!write(std_out,*)' '
!write(std_out,*)' symfind : has selected iclass0=',iclass0
!write(std_out,*)'  # iatom  xred                          spinat       (spinatred if noncoll=1) '
!do iatom0=1,natomcl(iclass0)
!iatom=class(iatom0,iclass0)
!if(noncoll==0)then
!  write(std_out,'(2i4,6f10.4)' )iatom0,iatom,xred(:,iatom),spinat(:,iatom)
!else if(noncoll==1)then
!  write(std_out,'(2i4,9f10.4)' )iatom0,iatom,xred(:,iatom),spinat(:,iatom),spinatred(:,iatom)
!endif
!end do
!write(std_out,*)' '
!ENDDEBUG


 !represent nuclear dipole moments in reduced coords
 ABI_ALLOCATE(nucdipmomred,(3,3,natom))
 do iatom=1,natom
    do ii=1,3
       nucdipmomred(1:3,ii,iatom)=MATMUL(TRANSPOSE(gprimd),local_nucdipmom(1:3,ii,iatom))
    end do
 end do

!Big loop over each symmetry operation of the Bravais lattice
 nsym=0
 do isym=1,nptsym

!  ji: Check whether symmetry operation leaves efield invariant
   if (berryopt==4 .or. berryopt==6 .or. berryopt==7 .or. &
&   berryopt==14 .or. berryopt==16 .or. berryopt==17) then
     efieldrot(:) = ptsymrel(:,1,isym)*efield(1) +  &
&     ptsymrel(:,2,isym)*efield(2) +  &
&     ptsymrel(:,3,isym)*efield(3)
     diff(:)=efield(:)-efieldrot(:)
     if( (diff(1)**2+diff(2)**2+diff(3)**2) > tolsym**2 ) cycle
   end if

   if (use_inversion==0) then
     det=ptsymrel(1,1,isym)*ptsymrel(2,2,isym)*ptsymrel(3,3,isym)+&
&     ptsymrel(2,1,isym)*ptsymrel(3,2,isym)*ptsymrel(1,3,isym)+&
&     ptsymrel(1,2,isym)*ptsymrel(2,3,isym)*ptsymrel(3,1,isym) - &
&     (ptsymrel(3,1,isym)*ptsymrel(2,2,isym)*ptsymrel(1,3,isym)+&
&     ptsymrel(2,1,isym)*ptsymrel(1,2,isym)*ptsymrel(3,3,isym)+&
&     ptsymrel(3,2,isym)*ptsymrel(2,3,isym)*ptsymrel(1,1,isym))
     if(det==-1) cycle
   end if

!  jellium slab and spatially varying chemical potential cases:
!  (actually, an inversion symmetry/mirror plane perpendicular to z symmetry operation might still be allowed... TO BE DONE !)
   if (jellslab/=0 .or. nzchempot/=0) then
!    check whether symmetry operation produce a rotation only in the xy plane
     if( ptsymrel(1,3,isym)/=0 .or. ptsymrel(2,3,isym)/=0 .or. &
&     ptsymrel(3,1,isym)/=0 .or. ptsymrel(3,2,isym)/=0 ) cycle
!    check whether symmetry operation does not change the z
     if( ptsymrel(3,3,isym)/=1 ) cycle
   end if

!  If noncoll_orthorhombic=1, require orthorhombic operations of symmetries, except if spinat=0.
   if (noncoll==1 .and. noncoll_orthorhombic==1)then
     if(sum(abs(spinat(:,:)))>tol14)then
       if( ptsymrel(1,3,isym)/=0 .or. ptsymrel(2,3,isym)/=0 .or. &
&          ptsymrel(1,2,isym)/=0 .or. ptsymrel(3,2,isym)/=0 .or. &
&          ptsymrel(2,1,isym)/=0 .or. ptsymrel(3,2,isym)/=0 ) cycle
     endif
   endif

!  Select a tentative set of associated translations
!  First compute the symmetric of the first atom in the smallest class,
!  using the point symmetry, and also the symmetric of spinat(red).
   iatom0=class(1,iclass0)
   xred0(:)=ptsymrel(:,1,isym)*xred(1,iatom0)+ &
&   ptsymrel(:,2,isym)*xred(2,iatom0)+ &
&   ptsymrel(:,3,isym)*xred(3,iatom0)
   if (noncoll==0) then
     spinat0(:)=spinat(:,iatom0)
   else
     spinat0(:)=ptsymrel(:,1,isym)*spinatred(1,iatom0)+ &
&           ptsymrel(:,2,isym)*spinatred(2,iatom0)+ &
&           ptsymrel(:,3,isym)*spinatred(3,iatom0)
   endif

!  From the set of possible images, deduce tentative translations,
!  and magnetic factor then test whether it send each atom on a symmetric one
   ntrial=0
   do ii=1,natom0
!DEBUG
!    write(std_out,'(a,2i4)')' symfind : loop isym,ii=',isym,ii
!ENDDEBUG
     iatom1=class(ii,iclass0)

!    The tentative translation is found
     trialnons(:)=xred(:,iatom1)-xred0(:)
!    Compare the spinat vectors
     if (noncoll==0) then
       symspinat1(:)=spinat(:,iatom1)
     else
       symspinat1(:)=spinatred(:,iatom1)
     end if

!DEBUG
!    write(std_out,'(a,6f10.4)')' symspinat1,spinat0=',symspinat1(:),spinat0(:)
!ENDDEBUG

     trialafm=1
     if(sum(abs(symspinat1(:)-spinat0(:)))>tolsym)then
       trialafm=-1
       if(noncoll==1 .and. afm_noncoll==0)cycle
       if(sum(abs(symspinat1(:)+spinat0(:)))>tolsym)cycle
     endif

     if(sum(abs(local_nucdipmom(:,1,iatom1)-local_nucdipmom(:,1,iatom0)))>tolsym)then
       write(message,'(3a,3i5)')&
&       'Problem with matching the nuclear dipole moment within a class.',ch10,&
&       'isym,iatom0,iatom1=',isym,iatom0,iatom1
       MSG_ERROR_CLASS(message, "TolSymError")
     end if
!    jellium slab case: check whether symmetry operation has no translational
!    component along z
     if( jellslab/=0 .and. abs(trialnons(3)) > tolsym ) cycle
     trialok=1

!    DEBUG
!    write(std_out, '(a,i3,a,i3,a,i3,a,3f12.4,i3)') ' Try isym=',isym,' sending iatom0 ',iatom0,' to iatom1 ',iatom1,' with trialnons(:),trialafm =',trialnons(:),trialafm
!    ENDDEBUG

!    Loop over all classes, then all atoms in the class,
!    to find whether they have a symmetric
     do iclass=1,nclass
       do jj=1,natomcl(iclass)

         iatom2=class(jj,iclass)
!        Generate the tentative symmetric position of iatom2
         symxred2(:)=ptsymrel(:,1,isym)*xred(1,iatom2)+ &
&         ptsymrel(:,2,isym)*xred(2,iatom2)+ &
&         ptsymrel(:,3,isym)*xred(3,iatom2)+ trialnons(:)
!        Generate the tentative symmetric spinat of iatom2
         if (noncoll==0) then
           symspinat2(:)=trialafm*spinat(:,iatom2)
         else
           symspinat2(:)=trialafm*(ptsymrel(:,1,isym)*spinatred(1,iatom2)+ &
&           ptsymrel(:,2,isym)*spinatred(2,iatom2)+ &
&           ptsymrel(:,3,isym)*spinatred(3,iatom2))
         end if
         !        Generate the tentative symmetric nucdipmom of iatom2
         do kk = 1, 3
            symnucdipmom2red(:,kk)=ptsymrel(:,1,isym)*nucdipmomred(1,kk,iatom2)+ &
                 &           ptsymrel(:,2,isym)*nucdipmomred(2,kk,iatom2)+ &
                 &           ptsymrel(:,3,isym)*nucdipmomred(3,kk,iatom2)
            ! transform back to cart coords for final comparison to nucdipmom
            symnucdipmom2cart(:,kk)=MATMUL(rprimd,symnucdipmom2red(:,kk))
         end do

!        DEBUG
!        write(std_out,'(a,i4,a,3f8.4,a,3f8.4,a,3f8.4)')&
!&          ' Test iatom2=',iatom2,' at xred=',xred(:,iatom2),'. Is sent to',symxred2(:),' with symspinat2=',symspinat2(:)
!        ENDDEBUG

!        Check whether there exists an atom of the same class at the
!        same location, with the correct spinat and nuclear dipole moment circulation
         do kk=1,natomcl(iclass)

           found3=1
           iatom3=class(kk,iclass)
!          Check the location
           diff(:)=xred(:,iatom3)-symxred2(:)
           diff(:)=diff(:)-nint(diff(:))
           if( (diff(1)**2+diff(2)**2+diff(3)**2) > tolsym**2 )found3=0
!          Check the spinat
           if (noncoll==0) then
             diff(:)=spinat(:,iatom3)-symspinat2(:)
           else
             diff(:)=spinatred(:,iatom3)-symspinat2(:)
           end if
           if( (diff(1)**2+diff(2)**2+diff(3)**2) > tolsym**2 )found3=0
           !          Check the nucdipmom
           ! hand3 gives original circulation sense of nuclear dipole
           call acrossb(local_nucdipmom(1:3,2,iatom3),local_nucdipmom(1:3,3,iatom3),hand3)

           ! hand2 gives circulation sense of tentative, symmetry equivalent nuclear dipole
           call acrossb(symnucdipmom2cart(1:3,2),symnucdipmom2cart(1:3,3),hand2)

           diff(:)=hand3(:)-hand2(:)
           if( any(abs(diff)>tolsym) )found3=0

           if(found3==1)exit
         end do ! End loop over iatom3

         if(found3==0)then
           trialok=0
           exit
         end if
       end do ! End loop over iatom2

       if(trialok==0)exit
     end do ! End loop over all classes

!DEBUG
!    write(std_out,*)' For trial isym=',isym,', trialok = ',trialok
!    write(std_out,*)' '
!ENDDEBUG

     if(trialok==1)then
       nsym=nsym+1
       if(nsym>msym)then
         write(message,'(a,i0,2a,i0,4a)')&
         'The number of symmetries (including non-symmorphic translations) is:', nsym, ch10,&
         'is larger than maxnsym: ',msym,ch10,&
         'Action: increase maxnsym in the input, or take a cell that is primitive, ',ch10,&
         'or at least smaller than the present one.'
        MSG_ERROR(message)
       end if
       ntrial=ntrial+1
       symrel(:,:,nsym)=ptsymrel(:,:,isym)
       symafm(nsym)=trialafm
       tnons(:,nsym)=trialnons(:)-nint(trialnons(:)-tolsym)
     end if

   end do ! End the loop on tentative translations
 end do ! End big loop over each symmetry operation of the Bravais lattice

 ABI_DEALLOCATE(class)
 ABI_DEALLOCATE(natomcl)
 ABI_DEALLOCATE(chrgat_)
 ABI_DEALLOCATE(chrgatcl)
 ABI_DEALLOCATE(spinatcl)
 ABI_DEALLOCATE(typecl)
 ABI_DEALLOCATE(local_nucdipmom)
 ABI_DEALLOCATE(nucdipmomcl)
 ABI_DEALLOCATE(nucdipmomred)
 if (noncoll==1)   then
   ABI_DEALLOCATE(spinatred)
 end if

 call chkgrp(nsym,symafm,symrel,ierr_)
 if (ierr_/=0) then
   call print_symmetries(nsym,symrel,tnons,symafm)
 end if

 if(.not.present(ierr))then
   ABI_CHECK(ierr_==0,"Error in group closure")
 else
   ierr=ierr_
 endif

!DEBUG
! write(message,'(a,I0,a)')' symfind : exit, nsym=',nsym,ch10
! write(message,'(2a)') trim(message),'   symrel matrices, symafm and tnons are :'
! call wrtout(std_out,message,'COLL')
! do isym=1,nsym
!   write(message,'(i4,4x,3i4,2x,3i4,2x,3i4,4x,i4,4x,3f8.4)' ) isym,symrel(:,:,isym),&
!&   symafm(isym),tnons(:,isym)
!   call wrtout(std_out,message,'COLL')
! end do
!stop
!ENDDEBUG

end subroutine symfind
!!***

!!****f* m_symfind/symanal
!! NAME
!! symanal
!!
!! FUNCTION
!! Find the space group, Bravais lattice, including Shubnikov characteristics
!! from the list of symmetries (including magnetic characteristics), and lattice parameters
!! Warning: the recognition of the space group might not yet work for the
!! Shubnikov group of type IV
!!
!! INPUTS
!! chkprim= if 1 then stop if the cell is not primitive
!! msym=default maximal number of symmetries
!! nsym=actual number of symmetries
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!! symafm(1:msym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,1:msym)=symmetry operations in real space in terms
!!  of primitive translations
!! tnons(3,1:msym)=nonsymmorphic translations for symmetry operations
!! tolsym=tolerance for the symmetry operations
!! verbose= if true, will list the symmetry operation labels
!!
!! OUTPUT
!! bravais(11)=characteristics of Bravais lattice (see symlatt.F90)
!! genafm(3)=magnetic translation generator (in case of Shubnikov group type IV)
!! ptgroupma = magnetic point group number
!! spgroup=symmetry space group
!!
!! PARENTS
!!      m_ab7_symmetry,m_ingeo,m_mover_effpot,m_phonons,m_tdep_sym,m_use_ga
!!
!! CHILDREN
!!      holocell,matr3inv,smallprim,symrelrot,wrtout
!!
!! SOURCE

subroutine symanal(bravais,chkprim,genafm,msym,nsym,ptgroupma,rprimd,spgroup,symafm,symrel,tnons,tolsym,verbose)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: chkprim,msym,nsym
 integer,intent(out) :: ptgroupma,spgroup
 real(dp),intent(in) :: tolsym
 logical,optional,intent(in) :: verbose
!arrays
 integer,intent(out) :: bravais(11)
 integer,intent(in) :: symafm(msym),symrel(3,3,msym)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: tnons(3,msym)
 real(dp),intent(out) :: genafm(3)

!Local variables-------------------------------
!scalars
 integer, parameter :: maxsym=192
! In this routine, maxsym is used either to determine the ptsymrel from the rprimd (routine symlatt),
! so, it might be up to 192 = 4*48 for FCC, and also to define the maximum number of symmetry operation labels,
! but only in case the cell is primitive, which gives the same upper bound. Thus in this routine, msym might
! be equal to nsym.
 integer :: iholohedry_nomagn,isym,isym_nomagn,multi
 integer :: nptsym,nsym_nomagn,shubnikov
 logical :: verbose_
 character(len=5) :: ptgroup,ptgroupha
 character(len=500) :: message
!arrays
 integer :: identity(3,3)
 integer,allocatable :: ptsymrel(:,:,:),symrel_nomagn(:,:,:)
 real(dp),allocatable :: tnons_nomagn(:,:)
 character(len=128) :: labels(maxsym)

! *************************************************************************

!DEBUG
!write(std_out,*)' symanal : enter'
!ENDDEBUG

 verbose_=.false.
 if(present(verbose))then
   verbose_=verbose
 endif

!This routine finds the Bravais characteristics, without actually
!looking at the symmetry operations.
 ABI_ALLOCATE(ptsymrel,(3,3,maxsym))
 call symlatt(bravais,maxsym,nptsym,ptsymrel,rprimd,tolsym)
 ABI_DEALLOCATE(ptsymrel)

!Check whether the cell is primitive or not.
 call chkprimit(chkprim,multi,nsym,symafm,symrel)

 spgroup=0 ; ptgroupma=0 ; genafm(:)=zero

 if(multi>1)then !  Modify bravais if the cell is not primitive ; no determination of the space group
   bravais(1)=-bravais(1)
 else

!  The cell is primitive, so that the space group can be
!  determined. Need to distinguish Fedorov and Shubnikov groups.
!  Do not distinguish Shubnikov types I and II.
!  Also identify genafm, in case of Shubnikov type IV
   identity(:,:)=reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
   shubnikov=1
   do isym=1,nsym
     if(symafm(isym)==-1)then
       shubnikov=3
       if(sum(abs(symrel(:,:,isym)-identity(:,:)))==0)then
         shubnikov=4
         genafm(:)=tnons(:,isym)
!        DEBUG
!        write(std_out,*)' isym=',isym
!        write(std_out,*)' symrel(:,:,isym)',symrel(:,:,isym)
!        write(std_out,*)' tnons(:,isym)',tnons(:,isym)
!        write(std_out,*)' symafm(isym)',symafm(isym)
!        ENDDEBUG
         exit
       end if
     end if
   end do

   if(shubnikov/=1)then
     if(shubnikov==3)write(message, '(a)' )' Shubnikov space group type III'
     if(shubnikov==4)write(message, '(a)' )' Shubnikov space group type IV'
     call wrtout(std_out,message,'COLL')
   end if

   if(shubnikov==1 .or. shubnikov==3)then
!    Find the correct Bravais characteristics and point group
!    Should also be used for Shubnikov groups of type IV ...
     call symbrav(bravais,msym,nsym,ptgroup,rprimd,symrel,tolsym)

!    Find the space group
     call symspgr(bravais,labels,nsym,spgroup,symrel,tnons,tolsym)

     if(verbose_)then
       do isym=1,nsym
         write(message,'(a,i3,2a)')' symanal : the symmetry operation no. ',isym,' is ',trim(labels(isym))
         call wrtout(std_out,message,'COLL')
       enddo
     endif

   end if

   if(shubnikov/=1)then

!    Determine nonmagnetic symmetry operations
     nsym_nomagn=nsym/2
     ABI_ALLOCATE(symrel_nomagn,(3,3,nsym_nomagn))
     ABI_ALLOCATE(tnons_nomagn,(3,nsym_nomagn))
     isym_nomagn=0
     do isym=1,nsym
       if(symafm(isym)==1)then
         isym_nomagn=isym_nomagn+1
         symrel_nomagn(:,:,isym_nomagn)=symrel(:,:,isym)
         tnons_nomagn(:,isym_nomagn)=tnons(:,isym)
       end if
     end do

     if(shubnikov==3)then

!      DEBUG
!      write(std_out,*)' symanal : will enter symbrav with halved symmetry set'
!      write(std_out,*)' Describe the different symmetry operations (index,symrel,tnons,symafm)'
!      do isym=1,nsym_nomagn
!      write(std_out,'(i3,2x,9i3,3es12.2,i3)')isym,symrel_nomagn(:,:,isym),tnons_nomagn(:,isym)
!      end do
!      ENDDEBUG

!      Find the point group of the halved symmetry set
       call symptgroup(iholohedry_nomagn,nsym_nomagn,ptgroupha,symrel_nomagn)

!      Deduce the magnetic point group (ptgroupma) from ptgroup and ptgroupha
       call getptgroupma(ptgroup,ptgroupha,ptgroupma)

     else if(shubnikov==4)then

!      Find the Fedorov space group of the halved symmetry set
       call symspgr(bravais,labels,nsym_nomagn,spgroup,symrel_nomagn,tnons_nomagn,tolsym)

!      The magnetic translation generator genafm has already been determined
!      write(std_out,*)' genafm =',genafm, ' spgroup=',spgroup

       if(verbose_)then

         write(message, '(a)' )' Select only the non-magnetic symmetry operations '
         call wrtout(std_out,message,'COLL')

         do isym=1,nsym
           if(symafm(isym)==1)then
             isym_nomagn=isym_nomagn+1
             write(message,'(a,i3,2a)')' symspgr : the symmetry operation no. ',isym,' is ',trim(labels(isym_nomagn))
             call wrtout(std_out,message,'COLL')
           endif
         enddo
       endif

     end if

     ABI_DEALLOCATE(symrel_nomagn)
     ABI_DEALLOCATE(tnons_nomagn)
   end if ! Shubnikov groups

 end if

!DEBUG
!write(std_out,'(a)') ' symanal : exit '
!ENDDEBUG

end subroutine symanal
!!***

!!****f* m_symfind/symbrav
!! NAME
!! symbrav
!!
!! FUNCTION
!! From the list of symmetry operations, and the lattice vectors,
!! determine the Bravais information (including the holohedry, the centering,
!! the coordinate of the primitive vectors in the conventional vectors),
!! as well as the point group.
!!
!! INPUTS
!! msym=dimension of symrel
!! nsym=actual number of symmetries
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!! symrel(3,3,msym)=symmetry operations in real space in terms
!!                  of primitive translations
!! tolsym=tolerance for the symmetries
!!
!! OUTPUT
!! bravais(11): bravais(1)=iholohedry
!!              bravais(2)=center
!!              bravais(3:11)=coordinates of rprimd in the axes
!!              of the conventional bravais lattice (*2 if center/=0)
!! ptgroup=symmetry point group
!! [axis(3)]=Invariant axis in the conventional vector coordinates
!!   Set to (/0,0,0/) if the lattice belongs to the same holohedry as the lattice+atoms (+electric field + ...).
!!
!! PARENTS
!!      m_esymm,m_symfind
!!
!! CHILDREN
!!      holocell,matr3inv,smallprim,symrelrot,wrtout
!!
!! SOURCE

subroutine symbrav(bravais,msym,nsym,ptgroup,rprimd,symrel,tolsym,axis)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: msym,nsym
 real(dp),intent(in) :: tolsym
 character(len=5),intent(out) :: ptgroup
!arrays
 integer,intent(in) :: symrel(3,3,msym)
 integer,optional,intent(out) :: axis(3)
 integer,intent(out) :: bravais(11)
 real(dp),intent(in) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: iaxis,ii,bravais1now,ideform,iholohedry,invariant,isym
 integer :: jaxis,next_stage,nptsym,problem,maxsym
 real(dp) :: norm,scprod
 character(len=500) :: message
!arrays
 integer :: identity(3,3),axis_trial(3),hexa_axes(3,7),ortho_axes(3,13)
 integer,allocatable :: ptsymrel(:,:,:),symrelconv(:,:,:)
 real(dp) :: axes(3,3),axis_cart(3),axis_red(3)
 real(dp) :: rprimdconv(3,3),rprimdtry(3,3),rprimdnow(3,3)
 real(dp) :: rprimdconv_invt(3,3)

!**************************************************************************

!DEBUG
!write(std_out,*)' symbrav : enter '
!call flush(std_out)
!ENDDEBUG

 identity(:,:)=0
 identity(1,1)=1 ; identity(2,2)=1 ; identity(3,3)=1

 ortho_axes(:,:)=0
 ortho_axes(1,1)=1
 ortho_axes(2,2)=1
 ortho_axes(3,3)=1
 ortho_axes(:,4)=(/0,1,1/)
 ortho_axes(:,5)=(/1,0,1/)
 ortho_axes(:,6)=(/1,1,0/)
 ortho_axes(:,7)=(/0,1,-1/)
 ortho_axes(:,8)=(/-1,0,1/)
 ortho_axes(:,9)=(/1,-1,0/)
 ortho_axes(:,10)=(/1,1,1/)
 ortho_axes(:,11)=(/-1,1,1/)
 ortho_axes(:,12)=(/1,-1,1/)
 ortho_axes(:,13)=(/1,1,-1/)

 hexa_axes(:,:)=0
 hexa_axes(1,1)=1
 hexa_axes(2,2)=1
 hexa_axes(3,3)=1
 hexa_axes(:,4)=(/1,-1,0/)
 hexa_axes(:,5)=(/2,1,0/)
 hexa_axes(:,6)=(/1,1,0/)
 hexa_axes(:,7)=(/1,2,0/)

!Determine the point group from the list of symmetry operations.
!Also determine the holohedry, up to one undeterminacy : hR versus hP
 call symptgroup(iholohedry,nsym,ptgroup,symrel)

!DEBUG
!write(std_out,*)' symbrav, after symptgroup: nsym=',nsym
!call flush(std_out)
!write(std_out,*)' symbrav: symrel='
!do isym=1,nsym
!  write(std_out,'(9i4)')symrel(:,:,isym)
!enddo
!write(std_out,*)' symbrav: iholohedry=',iholohedry
!call flush(std_out)
!ENDDEBUG

!Loop over trial deformations
!This is needed in case the Bravais lattice determination from the lattice vectors
!has a higher holohedry than the real one, in which the symmetry
!operations for the atoms (or electric field, etc) are taken into account
 iaxis=0
 invariant=0
 next_stage=0
 rprimdnow(:,:)=rprimd(:,:)
 rprimdtry(:,:)=rprimd(:,:)
 ABI_ALLOCATE(symrelconv,(3,3,nsym))

!At most will have to try 65 deformations (13 axes, five stages)
 do ideform=1,65

!DEBUG
!write(std_out,*)' symbrav: inside loop with ideform=',ideform
!call flush(std_out)
!write(std_out,'(a,9f10.4)')' rprimdtry=',rprimdtry(:,:)
!ENDDEBUG

   maxsym=max(192,msym)
   ABI_ALLOCATE(ptsymrel,(3,3,maxsym))
   call symlatt(bravais,maxsym,nptsym,ptsymrel,rprimdtry,tolsym)
   ABI_DEALLOCATE(ptsymrel)

!  Examine the agreement with bravais(1)
!  Warning : might change Bravais lattice hR to hP, if hexagonal axes
   problem=0
   select case (bravais(1))
   case(7)
     if(iholohedry<6)problem=1
     if(iholohedry==6)problem=2
   case(6)
     if(iholohedry<4)problem=1
     if(iholohedry==7 .or. iholohedry==4)problem=2
!      Here, change hR into hP
     if(iholohedry==5)iholohedry=6
   case(5)
     if(iholohedry<4)problem=1
     if(iholohedry==7 .or. iholohedry==6 .or. iholohedry==4)problem=2
   case(4)
     if(iholohedry<4)problem=1
     if(iholohedry>4)problem=2
   case(3)
     if(iholohedry<3)problem=1
     if(iholohedry>3)problem=2
   case(2)
     if(iholohedry<2)problem=1
     if(iholohedry>2)problem=2
   case(1)
     if(iholohedry>1)problem=2
   end select

!  This is the usual situation, in which the lattice belong to the same holohedry
!  as the lattice+atoms (+electric field + ...)
   if(problem==0)exit

   if(problem==2)then
     if(iaxis==0)then
       write(message, '(3a,i3,3a,i3,7a)' )&
&       'The Bravais lattice determined only from the primitive',ch10,&
&       'vectors (rprim or angdeg), bravais(1)=',bravais(1),', is not compatible',ch10,&
&       'with the real one, iholohedry=',iholohedry,', obtained by taking into',ch10,&
&       'account the symmetry operations. This might be due to an insufficient',ch10,&
&       'number of digits in the specification of rprim (at least 10),',ch10,&
&       'or to an erroneous rprim or angdeg. If this is not the case, then ...'
       MSG_BUG(message)
     end if
     if(iaxis==1)then
       write(message, '(3a,3i3,2a,i3,2a,i3)' )&
&       'Could not succeed to determine the bravais lattice',ch10,&
&       'problem,iaxis,invariant=',problem,iaxis,invariant,ch10,&
&       'bravais(1)=',bravais(1),ch10,&
&       'iholohedry=',iholohedry
       MSG_BUG(message)
     end if
   end if

   if(problem==1)then  ! One is left with the problem=1 case, basically iholohedry is lower than bravais(1)
     if(iaxis==0)then
       write(message, '(a,a,a,i3,a,a,a,i3,a,a,a)' )&
&       'The Bravais lattice determined only from the primitive',ch10,&
&       'vectors, bravais(1)=',bravais(1),', is more symmetric',ch10,&
&       'than the real one, iholohedry=',iholohedry,', obtained by taking into',ch10,&
&       'account the atomic positions. Start deforming the primitive vector set.'
       MSG_COMMENT(message)
       next_stage=1
     else if(iaxis/=0)then
       if(bravais(1)<bravais1now)then
         write(message, '(3a,i3,3a,i3,2a)' )&
&         'The Bravais lattice determined from modified primitive',ch10,&
&         'vectors, bravais(1)=',bravais(1),', has a lower symmetry than before,',ch10,&
&         'but is still more symmetric than the real one, iholohedry=',iholohedry,ch10,&
&         'obtained by taking into account the atomic positions.'
         MSG_COMMENT(message)
         next_stage=1
       else if(iaxis==1)then
         write(message, '(3a,3i3,2a,i3,2a,i3)' )&
&         'Could not succeed to determine the bravais lattice',ch10,&
&         'problem,iaxis,invariant=',problem,iaxis,invariant,ch10,&
&         'bravais(1)=',bravais(1),ch10,&
&         'iholohedry=',iholohedry
         MSG_BUG(message)
       end if
     end if
   end if ! problem==1

   if(next_stage==1)then
     bravais1now=bravais(1)
     rprimdnow(:,:)=rprimdtry(:,:)
!    Generate the symmetry operations in the conventional vector coordinates
     rprimdconv(:,1)=bravais(3:5)
     rprimdconv(:,2)=bravais(6:8)
     rprimdconv(:,3)=bravais(9:11)
     axes(:,:)=zero
     axes(1,1)=one ; axes(2,2)=one ; axes(3,3)=one
     symrelconv(:,:,1:nsym)=symrel(:,:,1:nsym)
     call symrelrot(nsym,rprimdconv,axes,symrelconv,tolsym)
     if(bravais(1)/=6)then
       iaxis=14
     else
       iaxis=8
     end if
     next_stage=0
   end if

   iaxis=iaxis-1
   do jaxis=iaxis,1,-1
     if(bravais(1)/=6)then
       axis_trial(:)=ortho_axes(:,jaxis)
     else
       axis_trial(:)=hexa_axes(:,jaxis)
     end if
!    DEBUG
!    write(std_out,*)' symbrav : try jaxis=',jaxis
!    write(std_out,*)' axis_trial=',axis_trial
!    ENDDEBUG
     invariant=1
!    Examine whether all symmetry operations leave the axis invariant (might be reversed, though)
     do isym=1,nsym
       if(sum(abs(matmul(symrelconv(:,:,isym),axis_trial)+(-axis_trial(:))))/=0 .and. &
&       sum(abs(matmul(symrelconv(:,:,isym),axis_trial)+axis_trial(:)))/=0 )invariant=0
     end do
     if(invariant==1)then
       iaxis=jaxis
!      write(message, '(2a,i3)' )ch10,' symbrav : found invariant axis, jaxis=',iaxis
!      call wrtout(std_out,message,'COLL')
       exit
     end if
   end do

   if(invariant==0)then
!    Not a single axis was invariant with respect to all operations ?!
!    do isym=1,nsym; write(std_out, '(a,10i4)' )' isym,symrelconv=',isym,symrelconv(:,:,isym); enddo
     write(message, '(3a,3i3,2a,i3,2a,i3)' )&
&     'Could not succeed to determine the bravais lattice (not a single invariant)',ch10,&
&     'problem,iaxis,invariant=',problem,iaxis,invariant,ch10,&
&     'bravais(1)=',bravais(1),ch10,&
&     'iholohedry=',iholohedry
     MSG_BUG(message)
   end if

   call matr3inv(rprimdconv,rprimdconv_invt)
   axis_red(:)=axis_trial(1)*rprimdconv_invt(1,:)+ &
&   axis_trial(2)*rprimdconv_invt(2,:)+ &
&   axis_trial(3)*rprimdconv_invt(3,:)
   axis_cart(:)=axis_red(1)*rprimdnow(:,1)+ &
&   axis_red(2)*rprimdnow(:,2)+ &
&   axis_red(3)*rprimdnow(:,3)
   norm=sum(axis_cart(:)**2)
!  Expand by a uniform, quite arbitrary, dilatation, along the invariant axis
!  Note : make these dilatation different, according to ideform
!  XG 20151221  : Still, the interplay between the size of the deformation and the tolsym is not easy to address.
!  Indeed the deformation must be sufficiently large to be perceived by symlatt as a real breaking of the
!  symmetry of the lattice. In order to deal with all the small values od tolsym, it has been set at a minimum of tol3,
!  but it must also be larger than tolsym. Moreover, for some axis choice, the deformation is not aligned with the axis, decreasing
!  the effective deformation length. An additional factor of three is thus included, actually increased to six just to be sure...
   do ii=1,3
     scprod=axis_cart(1)*rprimdnow(1,ii)+axis_cart(2)*rprimdnow(2,ii)+axis_cart(3)*rprimdnow(3,ii)
     rprimdtry(:,ii)=rprimdnow(:,ii)+ideform*(max(tol3,six*tolsym)-tol6)*scprod/norm*axis_cart(:)
   end do

 end do ! ideform

 if(bravais(1)/=iholohedry)then
   write(message, '(3a,3i3,2a,i3,2a,i3)' )&
&   'Despite efforts, Could not succeed to determine the bravais lattice :',ch10,&
&   'bravais(1)=',bravais(1),ch10,&
&   'iholohedry=',iholohedry
   MSG_BUG(message)
 end if

 ABI_DEALLOCATE(symrelconv)

 if (PRESENT(axis)) then  ! Return symmetry axis.
   axis=(/0,0,0/)
   if (iaxis/=0) then
     if(bravais(1)/=6)then
       axis=ortho_axes(:,iaxis)
     else
       axis=hexa_axes(:,iaxis)
     end if
   end if
 end if

!DEBUG
!write(std_out,'(a)')' symbrav : exit '
!ENDDEBUG

end subroutine symbrav
!!***

!!****f* m_symfind/symspgr
!! NAME
!! symspgr
!!
!! FUNCTION
!! Find the type of each symmetry operation (calling symcharac):
!!   proper symmetries 1,2,2_1,3,3_1,3_2,4,4_1,4_2,4_3,6,6_1,...6_5
!!   improper symmetries -1,m,a,b,c,d,n,g,-3,-4,-6 ,
!! Then, build an array with the number of such operations.
!! Then, call symlist to identify the space group.
!!
!! INPUTS
!! bravais(11): bravais(1)=iholohedry
!!              bravais(2)=center
!!              bravais(3:11)=coordinates of rprimd in the axes
!!              of the conventional bravais lattice (*2 if center/=0)
!! nsym=actual number of symmetries
!! symrel(3,3,nsym)= nsym symmetry operations in real space in terms
!!   of primitive translations
!! tnons(3,nsym)=nonsymmorphic translations for each symmetry (would
!!   be 0 0 0 each for a symmorphic space group)
!!
!! OUTPUT
!! labels(maxsym=192)= labels of the symmetry operations
!! spgroup=symmetry space group number
!!
!! NOTES
!! It is assumed that the symmetry operations will be entered in the
!! symrel tnons arrays, for the PRIMITIVE cell. The matrix of transformation
!! from the primitive cell to the conventional cell is described
!! in the array "bravais" (see symlatt.F90).
!! The present routine first make the transformation from the
!! primitive coordinates to the conventional ones, then eventually
!! generate additional symmetries, taking into account the
!! centering translations.
!! Then, the order and determinant of each symmetry operation
!! is determined.
!!
!! For proper symmetries (rotations), the
!! associated translation is also determined.
!! However, left or right handed screw rotations are
!! not (presently) distinguished, and will be attributed equally
!! to left or right.
!!
!! For the detailed description of the labelling of the axes,
!! see symaxes.f and symplanes.f
!!
!! PARENTS
!!      m_symfind
!!
!! CHILDREN
!!      holocell,matr3inv,smallprim,symrelrot,wrtout
!!
!! SOURCE

subroutine symspgr(bravais,labels,nsym,spgroup,symrel,tnons,tolsym)

 use m_numeric_tools, only : OPERATOR(.x.)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(out) :: spgroup
 real(dp),intent(in) :: tolsym
!arrays
 integer,intent(in) :: bravais(11),symrel(3,3,nsym)
 real(dp),intent(in) :: tnons(3,nsym)
 character(len=128),intent(out) :: labels(192) ! 192 = maxsym

!Local variables-------------------------------
!scalars
! logical,parameter :: verbose=.FALSE.
 integer :: additional_info,brvltt,center,direction=0,found,iholohedry,ii
 integer :: ishift,isym,jj,nshift,nsymconv,spgaxor,spgorig,sporder
 character(len=1) :: brvsb
 character(len=15) :: intsb,ptintsb,ptschsb,schsb
 character(len=35) :: intsbl
 character(len=500) :: message
!arrays
 integer :: ivec1(3), ivec2(3)
 integer :: n_axes(31),n_axest(31),prime(5),test_direction(3),symrel_uni(3,3)
 integer :: uniaxis(3),uniaxis_try(3)
 integer,allocatable :: determinant(:),symrelconv(:,:,:),t_axes(:)
 real(dp) :: axes(3,3),rprimdconv(3,3),trialt(3),vect(3,3)
 real(dp),allocatable :: shift(:,:),tnonsconv(:,:)

!**************************************************************************

 DBG_ENTER("COLL")

!Initialize brvltt, from bravais(2) and bravais(1)
 center=bravais(2)
 iholohedry=bravais(1)
 brvltt=1
 if(center==-1)brvltt=2  ! Inner centering
 if(center==-3)brvltt=3  ! Face centering
 if(center==1)brvltt=5  ! A-Face centering
 if(center==2)brvltt=6  ! B-Face centering
 if(center==3)brvltt=4  ! C-Face centering
 if(iholohedry==5)brvltt=7  ! Rhombohedral

!Produce the symmetry operations, in the axis of the conventional cell
 nsymconv=nsym
 if(center/=0)nsymconv=2*nsymconv
 if(center==-3)nsymconv=4*nsym
 ABI_ALLOCATE(symrelconv,(3,3,nsymconv))
 ABI_ALLOCATE(tnonsconv,(3,nsymconv))

!Produce symrel and tnons in conventional axes,
!name them symrelconv and tnonsconv
 rprimdconv(:,1)=bravais(3:5)
 rprimdconv(:,2)=bravais(6:8)
 rprimdconv(:,3)=bravais(9:11)

 if(center/=0)rprimdconv(:,:)=rprimdconv(:,:)*half

 axes(:,:)=zero
 axes(1,1)=one ; axes(2,2)=one ; axes(3,3)=one
 symrelconv(:,:,1:nsym)=symrel(:,:,1:nsym)
!Note that the number of symmetry operations is still nsym
 call symrelrot(nsym,rprimdconv,axes,symrelconv,tolsym)
 call xred2xcart(nsym,rprimdconv,tnonsconv,tnons)
!Gives the associated translation, with components in the
!interval ]-0.5,0.5] .
 tnonsconv(:,1:nsym)=tnonsconv(:,1:nsym)-nint(tnonsconv(:,1:nsym)-tol6)

!If the Bravais lattice is centered, duplicate or quadruplicate
!the number of symmetry operations, using the Bravais
!lattice shifts
 nshift=1
 if(center/=0)nshift=2
 if(center==-3)nshift=4
 ABI_ALLOCATE(shift,(3,nshift))
 shift(:,1)=zero
 if(center/=0 .and. center/=-3)then
   shift(:,2)=half
   if(center==1)shift(1,2)=zero
   if(center==2)shift(2,2)=zero
   if(center==3)shift(3,2)=zero
 else if(center==-3)then
   shift(:,2)=half ; shift(1,2)=zero
   shift(:,3)=half ; shift(2,3)=zero
   shift(:,4)=half ; shift(3,4)=zero
 end if ! center/=0 or -3
 if(nshift/=1)then
   do ishift=2,nshift
     symrelconv(:,:,(ishift-1)*nsym+1:ishift*nsym)=symrelconv(:,:,1:nsym)
     do isym=1,nsym
       tnonsconv(:,(ishift-1)*nsym+isym)=tnonsconv(:,isym)+shift(:,ishift)
     end do
   end do ! ishift
 end if ! nshift/=1

!At this stage, all the symmetry operations are available,
!expressed in the conventional axis, and also include
!the Bravais lattive translations, and associated operations...

 n_axes(:)=0

 ABI_ALLOCATE(determinant,(nsymconv))

!Get the determinant
 call symdet(determinant,nsymconv,symrelconv)

!Get the order of each the symmetry operation, as well as the maximal order
!Also, examine whether each symmetry operation is the inversion, or a root
!of the inversion (like -3)
!Decide which kind of point symmetry operation it is
!Finally assign tnonsconv order and decide the space symmetry operation

 ABI_ALLOCATE(t_axes,(nsymconv))

 do isym=1,nsymconv

!  Note : nsymconv might be bigger than 192, but only for non-primitive cells, in which case labels will not be echoed anywhere.
!  192 is the fixed dimension of labels, so this avoids possible memory problems.
   call symcharac(center, determinant(isym), iholohedry, isym, labels(mod(isym-1,192)+1), &
   symrelconv(:,:,isym), tnonsconv(:,isym), t_axes(isym))
   if (t_axes(isym) == -1) then
     write(message, '(a,a,i3,a,3(a,3i4,a),a,3es22.12,a,a,3es22.12)' )ch10,&
&     ' symspgr: problem with isym=',isym,ch10,&
&     '  symrelconv(:,1,isym)=',symrelconv(:,1,isym),ch10,&
&     '  symrelconv(:,2,isym)=',symrelconv(:,2,isym),ch10,&
&     '  symrelconv(:,3,isym)=',symrelconv(:,3,isym),ch10,&
&     '  tnonsconv(:,isym)=',tnonsconv(:,isym),ch10,&
&     '  trialt(:)=',trialt(:)
     call wrtout(std_out,message,'COLL')
     write(message, '(a,i4,2a)' )&
&     'The space symmetry operation number',isym,ch10,'is not a (translated) root of unity'
     MSG_BUG(message)
   else if (t_axes(isym) == -2) then
     write(message, '(a,i0,a)' )'The symmetry operation number ',isym,' is not a root of unity'
     MSG_BUG(message)
   end if

   n_axes(t_axes(isym))=n_axes(t_axes(isym))+1

 end do ! isym=1,nsymconv

 if (sum(n_axes)-nsymconv/=0) then
   write(message, '(7a)' )&
&   'Not all the symmetries have been recognized. ',ch10,&
&   'This might be due either to an error in the input file',ch10,&
&   'or to a BUG in ABINIT',ch10,&
&   'Please contact the ABINIT group.'
   MSG_WARNING(message)
 end if

!DEBUG
!write(std_out,*)' symspgr : brvltt,nsymconv=',brvltt,nsymconv
!write(std_out,*)' n_axes(1:10)=',n_axes(1:10)
!write(std_out,*)' n_axes(11:20)=',n_axes(11:20)
!write(std_out,*)' n_axes(21:31)=',n_axes(21:31)
!ENDDEBUG

!Treat cases in which the space group cannot be identified on the
!basis of n_axes one need additional informations
 if(brvltt==1)then
!  If the bravais lattice is primitive
   if(nsymconv==4)then
     n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,0,0,0,0,2,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0)then    ! Spgroup 27 (Pcc2) or 32 (Pba2)
       write(std_out,*)' symspgr: 27 or 32'
       additional_info=2
!      Select binary axis
       do isym=1,nsymconv
         if(t_axes(isym)==8)then
!          Find direction of binary axis
           if(symrelconv(1,1,isym)==1)direction=1
           if(symrelconv(2,2,isym)==1)direction=2
           if(symrelconv(3,3,isym)==1)direction=3
         end if
       end do
!      Examine the projection of the translation vector of the a, b or c mirror planes
!      onto the binary axis
       do isym=1,nsymconv
         if(t_axes(isym)==16)then
           if(abs(tnonsconv(direction,isym))>tol8)additional_info=1
         end if
       end do
     end if
   else if(nsymconv==8)then
     n_axest=(/0,0,0,0,1,0,0,1,1,0,  0,0,0,0,1,2,0,0,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0)then    ! Spgroup 55 (Pbam) or 57 (Pbcm)
       write(std_out,*)' symspgr: 55 or 57'
       additional_info=1
!      Select mirror plane m
       do isym=1,nsymconv
         if(t_axes(isym)==15)then
!          Find direction of mirror plane
           if(symrelconv(1,1,isym)==-1)direction=1
           if(symrelconv(2,2,isym)==-1)direction=2
           if(symrelconv(3,3,isym)==-1)direction=3
         end if
       end do
!      Examine the projection of the translation vector of the a, b, or c mirror planes
!      onto the binary axis
       do isym=1,nsymconv
         if(t_axes(isym)==16)then
           if(abs(tnonsconv(direction,isym))>tol8)additional_info=2
         end if
       end do
     end if
     n_axest=(/0,0,0,0,1,0,0,1,1,0,  0,0,0,0,0,2,0,1,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0)then    ! Spgroup 56 (Pccn) or 60 (Pbcn)
       write(std_out,*)' symspgr: 56 or 60'
       additional_info=1
!      Select mirror plane n
       do isym=1,nsymconv
         if(t_axes(isym)==18)then
!          Find direction of mirror plane
           if(symrelconv(1,1,isym)==-1)direction=1
           if(symrelconv(2,2,isym)==-1)direction=2
           if(symrelconv(3,3,isym)==-1)direction=3
         end if
       end do
!      Examine the projection of the translation vector of the a, b, or c mirror planes
!      onto the binary axis
       do isym=1,nsymconv
         if(t_axes(isym)==16)then
           if(abs(tnonsconv(direction,isym))<tol8)additional_info=2
         end if
       end do
     end if
   end if
 else if(brvltt==2)then
!  In the few next lines, use additional_info as a flag
   additional_info=0
!  If the bravais lattice is inner-centered
   if(nsymconv==8)then
!    Test spgroup 23 (I222) or 24 (I2_{1}2_{1}2_{1})
     n_axest=(/0,0,0,0,0,0,1,1,3,0,  0,0,0,0,0,0,0,0,0,3,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) additional_info=1
   else if(nsymconv==24)then
!    Test spgroup 197 (I23) or 199 (I2_{1}3)
     n_axest=(/0,0,0,0,0,0,1,1,3,16, 0,0,0,0,0,0,0,0,0,3,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) additional_info=1
   end if
   if(additional_info==1)then
     write(std_out,*)' symspgr: (23 or 24) or (197 or 199)'
!    Select the three binary axes (they might be 2 or 2_1 !)
     test_direction(:)=0
     do isym=1,nsymconv
       if(t_axes(isym)==20)then
!        Find direction of axis
         do direction=1,3
           if(symrelconv(direction,direction,isym)==1)then
             test_direction(direction)=1
             if(abs(tnonsconv(direction,isym))<tol8)then
               vect(:,direction)=tnonsconv(:,isym)
             else
               vect(:,direction)=tnonsconv(:,isym)+half
             end if
             vect(:,direction)=vect(:,direction)-nint(vect(:,direction)-tol8)
             vect(direction,direction)=zero
           end if
         end do ! direction=1,3
       end if ! if binary axis
     end do ! isym
     if(test_direction(1)/=1 .or. test_direction(2)/=1 .and. test_direction(3)/=1)then
       write(message, '(5a,3i4)' )&
&       'For space groups 23, 24, 197 or 197, the three binary axes',ch10,&
&       'are not equally partitioned along the x, y and z directions',ch10,&
&       'test_direction(1:3)=',test_direction(:)
       MSG_BUG(message)
     end if
     additional_info=1
     if(abs(vect(1,2)-vect(1,3))>tol8 .or. &
&     abs(vect(2,1)-vect(2,3))>tol8 .or. &
&     abs(vect(3,1)-vect(3,2))>tol8) additional_info=2
   end if ! additional informations are needed
 end if ! brvltt==1

 if (brvltt==0 .or. brvltt==1) then ! Primitive
   call symlist_prim(additional_info,nsymconv,n_axes,spgroup)
 else if(brvltt==2)then
   call symlist_bcc(additional_info,nsymconv,n_axes,spgroup)
 else if(brvltt==3)then
   call symlist_fcc(nsymconv,n_axes,spgroup)
 else
   call symlist_others(brvltt,nsymconv,n_axes,spgroup)
 end if

 if(spgroup==0) then
   write(message, '(a,a,a,a,a)' )&
&   'Could not find the space group.',ch10,&
&   'This often happens when the user selects a restricted set of symmetries ',ch10,&
&   'in the input file, instead of letting the code automatically find symmetries.'
   MSG_WARNING(message)
 end if

 spgorig=1 ; spgaxor=1
 call spgdata(brvsb,intsb,intsbl,ptintsb,ptschsb,schsb,spgaxor,spgroup,sporder,spgorig)

 if(spgroup/=0)then
   write(message, '(a,i4,2x,a,a,a,a,a)' ) ' symspgr: spgroup=',spgroup,trim(brvsb),trim(intsb),'   (=',trim(schsb),')'
   call wrtout(std_out,message,'COLL')
 end if

 if(bravais(1)==7)then
   write(message, '(a)' ) ' symspgr: optical characteristics = isotropic '
   call wrtout(std_out,message,'COLL')
 else if(bravais(1)==4 .or. bravais(1)==5 .or. bravais(1)==6)then
   write(message, '(a)' ) ' symspgr: optical characteristics = uniaxial '
   call wrtout(std_out,message,'COLL')
!  Identify the first symmetry operation that is order 3, 4 or 6
   found=0
   do isym=1,nsym
!    Proper rotations
     if( minval( abs( t_axes(isym)-(/10,12,14,22,23,24,25,26,27,28,29,30,31/) ))==0) then
       found=1 ; exit
!   Improper symmetry operations
     else if( minval( abs( t_axes(isym)-(/1,2,3/) ))==0) then
       found=-1 ; exit
     end if
   end do
   if(found==-1 .or. found==1)then
     symrel_uni=symrel(:,:,isym)
     if(found==-1)symrel_uni=-symrel_uni
!    Now, symrel_uni is a rotation of order 3, 4, 6, for which the axis must be identified
!    It is actually the only eigenvector with eigenvalue 1. It can be found by cross products
!    Subtract the unit matrix.
     do ii=1,3
       symrel_uni(ii,ii)=symrel_uni(ii,ii)-1
     end do
     found=0
     do ii=1,3
       jj=ii+1 ; if(jj==4)jj=1
!      Cross product
       ivec1 = symrel_uni(ii,:); ivec2 = symrel_uni(jj,:)
       uniaxis = ivec1 .x. ivec2
       if(sum(uniaxis**2)/=0)then
         found=1 ; exit
       end if
     end do
     if(found==1)then
!      Try to reduce the length, by an integer factor (try only primes 2, 3, 5, 7, 11)
       prime=(/2,3,5,7,11/)
       ii=1
       do while (ii<6)
         uniaxis_try=uniaxis/prime(ii)
         if(sum(abs(uniaxis_try*prime(ii)-uniaxis))==0)then
           uniaxis=uniaxis_try
         else
           ii=ii+1
         end if
       end do
       write(message, '(a,3i4)' ) ' Optical axis (in reduced coordinates, real space ) :',uniaxis
     end if
   end if
   if(found==0)then
     write(message, '(a)' ) ' However, the axis has not been found. Sorry for this.'
   end if
   call wrtout(std_out,message,'COLL')
 end if

 ABI_DEALLOCATE(determinant)
 ABI_DEALLOCATE(shift)
 ABI_DEALLOCATE(symrelconv)
 ABI_DEALLOCATE(tnonsconv)
 ABI_DEALLOCATE(t_axes)

 DBG_EXIT("COLL")

end subroutine symspgr
!!***

!!****f* m_symfind/symlatt
!! NAME
!! symlatt
!!
!! FUNCTION
!! From the unit cell vectors (rprimd) and the corresponding metric tensor,
!! find the Bravais lattice and its symmetry operations (ptsymrel).
!! 1) Find the shortest possible primitive vectors for the lattice
!! 2) Determines the holohedral group of the lattice, and the
!!    axes to be used for the conventional cell
!!    (this is a delicate part, in which the centering of the
!!    reduced cell must be taken into account)
!!    The idea is to determine the basis vectors of the conventional
!!    cell from the reduced cell basis vectors.
!! 3) Generate the symmetry operations of the holohedral group
!!
!! INPUTS
!! msym=default maximal number of symmetries. WARNING : cannot be simply set to nsym, because
!!   the number of symmetries found here will likely be bigger than sym !
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!! tolsym=tolerance for the symmetries
!!
!! OUTPUT
!!  bravais(11): bravais(1)=iholohedry
!!               bravais(2)=center
!!               bravais(3:11)=coordinates of rprim in the axes
!!               of the conventional bravais lattice (*2 if center/=0)
!! nptsym=number of point symmetries of the Bravais lattice
!! ptsymrel(3,3,1:msym)= nptsym point-symmetry operations
!! of the Bravais lattice in real space in terms of primitive translations.
!!
!! NOTES
!! WARNING: bravais(1) might be given a negative value in another
!! routine, if the cell is non-primitive.
!! The holohedral groups are numbered as follows
!! (see international tables for crystallography (1983), p. 13)
!! iholohedry=1   triclinic      1bar
!! iholohedry=2   monoclinic     2/m
!! iholohedry=3   orthorhombic   mmm
!! iholohedry=4   tetragonal     4/mmm
!! iholohedry=5   trigonal       3bar m
!! iholohedry=6   hexagonal      6/mmm
!! iholohedry=7   cubic          m3bar m
!! Centering
!! center=0        no centering
!! center=-1       body-centered
!! center=-3       face-centered
!! center=1        A-face centered
!! center=2        B-face centered
!! center=3        C-face centered
!!
!! PARENTS
!!      m_ab7_symmetry,m_effective_potential_file,m_ingeo,m_inkpts
!!      m_mover_effpot,m_symfind,m_tdep_sym,m_thmeig,m_use_ga
!!
!! CHILDREN
!!      holocell,matr3inv,smallprim,symrelrot,wrtout
!!
!! SOURCE

subroutine symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: msym
 integer,intent(out) :: nptsym
 real(dp),intent(in) :: tolsym
!arrays
 integer,intent(out) :: bravais(11),ptsymrel(3,3,msym)
 real(dp),intent(in) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer,parameter :: mgen=4
 integer :: center,fact,found,foundc,ia,ib,icase,igen,iholohedry,ii,index,isym
 integer :: itrial,jj,jsym,ngen=0,orthogonal,sign12,sign13,sign23,sumsign
 real(dp) :: determinant,norm2a,norm2b,norm2c,norm2trial,reduceda,reducedb,sca
 real(dp) :: scalarprod,scb,trace,val
 character(len=500) :: message
!arrays
 integer,parameter :: list_holo(7)=(/7,6,4,3,5,2,1/)
 integer :: ang90(3),equal(3),gen(3,3,mgen),gen2xy(3,3),gen2y(3,3),gen2z(3,3)
 integer :: gen3(3,3),gen6(3,3),icoord(3,3),identity(3,3),nvecta(3),nvectb(3)
 integer :: order(mgen)
 real(dp) :: axes(3,3),axesinvt(3,3),cell_base(3,3),coord(3,3),metmin(3,3)
 real(dp) :: minim(3,3),scprods(3,3),vecta(3),vectb(3),vectc(3),vin1(3),vin2(3),vext(3)

!**************************************************************************

!DEBUG
!write(std_out,'(a)') ' symlatt : enter '
!ENDDEBUG

 identity(:,:)=0 ; identity(1,1)=1 ; identity(2,2)=1 ; identity(3,3)=1
 nvecta(1)=2 ; nvectb(1)=3
 nvecta(2)=1 ; nvectb(2)=3
 nvecta(3)=1 ; nvectb(3)=2

!--------------------------------------------------------------------------
!Reduce the input vectors to a set of minimal vectors
 call smallprim(metmin,minim,rprimd)

!DEBUG
!write(std_out,*)' symlatt : minim(:,1)=',minim(:,1)
!write(std_out,*)' symlatt : minim(:,2)=',minim(:,2)
!write(std_out,*)' symlatt : minim(:,3)=',minim(:,3)
!ENDDEBUG

!--------------------------------------------------------------------------
!Examine the angles and vector lengths
 ang90(:)=0
 if(metmin(1,2)**2<tolsym**2*metmin(1,1)*metmin(2,2))ang90(3)=1
 if(metmin(1,3)**2<tolsym**2*metmin(1,1)*metmin(3,3))ang90(2)=1
 if(metmin(2,3)**2<tolsym**2*metmin(2,2)*metmin(3,3))ang90(1)=1
 equal(:)=0
 if(abs(metmin(1,1)-metmin(2,2))<tolsym*half*(metmin(1,1)+metmin(2,2)))equal(3)=1
 if(abs(metmin(1,1)-metmin(3,3))<tolsym*half*(metmin(1,1)+metmin(3,3)))equal(2)=1
 if(abs(metmin(2,2)-metmin(3,3))<tolsym*half*(metmin(2,2)+metmin(3,3)))equal(1)=1

!DEBUG
!write(std_out,*)' ang90=',ang90(:)
!write(std_out,*)' equal=',equal(:)
!ENDDEBUG

!-----------------------------------------------------------------------
!Identification of the centering

 foundc=0
!Default values
 fact=1 ; center=0
 cell_base(:,:)=minim(:,:)

!Examine each holohedral group
!This search is ordered : should not be happy with tetragonal,
!while there is FCC ...
 do index=1,6

!  If the holohedry is already found, exit
   if(foundc==1)exit

!  Initialize the target holohedry
   iholohedry=list_holo(index)

!  DEBUG
!  write(std_out,*)' symlatt : trial holohedry',iholohedry
!  ENDDEBUG

   orthogonal=0
   if(iholohedry==7 .or. iholohedry==4 .or. iholohedry==3)orthogonal=1

!  Now, will examine different working hypothesis.
!  The set of these hypothesis is thought to cover all possible cases ...

!  Working hypothesis : the basis is orthogonal
   if(ang90(1)+ang90(2)+ang90(3)==3 .and. orthogonal==1)then
     fact=1 ; center=0
     cell_base(:,:)=minim(:,:)
!    Checks that the basis vectors are OK for the target holohedry
     call holocell(cell_base,0,foundc,iholohedry,tolsym)
   end if

!  Select one trial direction
   do itrial=1,3

!    If the holohedry is already found, exit
     if(foundc==1)exit

     ia=nvecta(itrial) ; ib=nvectb(itrial)

!    This is in case of hexagonal holohedry
     if(foundc==0 .and. iholohedry==6 .and. ang90(ia)==1 .and. ang90(ib)==1 .and. equal(itrial)==1 )then
       reduceda=metmin(ib,ia)/metmin(ia,ia)
       fact=1 ; center=0
       if(abs(reduceda+0.5d0)<tolsym)then
         cell_base(:,1)=minim(:,ia)
         cell_base(:,2)=minim(:,ib)
         cell_base(:,3)=minim(:,itrial)
!        Checks that the basis vectors are OK for the target holohedry
         call holocell(cell_base,0,foundc,iholohedry,tolsym)
       else if(abs(reduceda-0.5d0)<tolsym)then
         cell_base(:,1)=minim(:,ia)
         cell_base(:,2)=minim(:,ib)-minim(:,ia)
         cell_base(:,3)=minim(:,itrial)
!        Checks that the basis vectors are OK for the target holohedry
         call holocell(cell_base,0,foundc,iholohedry,tolsym)
       end if
     end if

!    Working hypothesis : the conventional cell is orthogonal,
!    and the two other vectors are axes of the conventional cell
     if(foundc==0 .and. orthogonal==1 .and. ang90(itrial)==1)then

!      Compute the reduced coordinate of trial vector in the basis
!      of the two other vectors
       reduceda=metmin(itrial,ia)/metmin(ia,ia)
       reducedb=metmin(itrial,ib)/metmin(ib,ib)
       cell_base(:,ia)=minim(:,ia)
       cell_base(:,ib)=minim(:,ib)
       if( (abs(abs(reduceda)-0.5d0)<tolsym .and. abs(reducedb)<tolsym ) .or. &
&       ( abs(reduceda)<tolsym .and. abs(abs(reducedb)-0.5d0)<tolsym)       )then
         if(abs(abs(reduceda)-0.5d0)<tolsym)center=ib
         if(abs(abs(reducedb)-0.5d0)<tolsym)center=ia
         fact=2
         cell_base(:,itrial)= &
&         (minim(:,itrial)-reduceda*minim(:,ia)-reducedb*minim(:,ib) )*2.0d0
         call holocell(cell_base,0,foundc,iholohedry,tolsym)
       else if( abs(abs(reduceda)-0.5d0)<tolsym .and.&
&         abs(abs(reducedb)-0.5d0)<tolsym       ) then
         fact=2 ; center=-1
         cell_base(:,itrial)= &
&         (minim(:,itrial)-reduceda*minim(:,ia)-reducedb*minim(:,ib) )*2.0d0
         call holocell(cell_base,0,foundc,iholohedry,tolsym)
       end if
     end if

!    Working hypothesis : the conventional cell is orthogonal, and
!    the trial vector is one of the future axes, and the face perpendicular to it is centered
     if(foundc==0 .and. iholohedry==3 .and. &
&     ang90(ia)==1 .and. ang90(ib)==1 .and. equal(itrial)==1 )then
       fact=2 ; center=itrial
       cell_base(:,ia)=minim(:,ia)+minim(:,ib)
       cell_base(:,ib)=minim(:,ia)-minim(:,ib)
       cell_base(:,itrial)=minim(:,itrial)
!      Checks that the basis vectors are OK for the target holohedry
       call holocell(cell_base,0,foundc,iholohedry,tolsym)
     end if

!    DEBUG
!    write(std_out,*)' after test_b, foundc=',foundc
!    ENDDEBUG

!    Working hypothesis : the conventional cell is orthogonal, and
!    the trial vector is one of the future axes
     if(foundc==0 .and. orthogonal==1)then
!      Compute the projection of the two other vectors on the trial vector
       reduceda=metmin(itrial,ia)/metmin(itrial,itrial)
       reducedb=metmin(itrial,ib)/metmin(itrial,itrial)
!      If both projections are half-integer, one might have found an axis
       if( abs(abs(reduceda)-0.5d0)<tolsym .and.&
&       abs(abs(reducedb)-0.5d0)<tolsym       ) then
         vecta(:)=minim(:,ia)-reduceda*minim(:,itrial)
         vectb(:)=minim(:,ib)-reducedb*minim(:,itrial)
         norm2a=vecta(1)**2+vecta(2)**2+vecta(3)**2
         norm2b=vectb(1)**2+vectb(2)**2+vectb(3)**2
         scalarprod=vecta(1)*vectb(1)+vecta(2)*vectb(2)+vecta(3)*vectb(3)
!        Note the order of selection : body-centered is prefered
!        over face centered, which is correct for the tetragonal case
         if(abs(norm2a-norm2b)<tolsym*half*(norm2a+norm2b))then
!          The lattice is body centered
           fact=2 ; center=-1
           cell_base(:,ia)=vecta(:)+vectb(:)
           cell_base(:,ib)=vecta(:)-vectb(:)
           cell_base(:,itrial)=minim(:,itrial)
           call holocell(cell_base,0,foundc,iholohedry,tolsym)
         else if(abs(scalarprod)<tolsym*half*(norm2a+norm2b))then
!          The lattice is face centered
           fact=2 ; center=-3
           cell_base(:,ia)=2.0d0*vecta(:)
           cell_base(:,ib)=2.0d0*vectb(:)
           cell_base(:,itrial)=minim(:,itrial)
           call holocell(cell_base,0,foundc,iholohedry,tolsym)
         end if
       end if
     end if

!    DEBUG
!    write(std_out,*)' after test_c, foundc=',foundc
!    ENDDEBUG

!    Working hypothesis : the conventional cell is orthogonal,
!    and body centered with no basis vector being an axis,
!    in which case the basis vectors must be equal (even for orthorhombic)
     if(foundc==0 .and. orthogonal==1 .and. &
&     equal(1)==1 .and. equal(2)==1 .and. equal(3)==1 )then
!      Compute the combination of the two other vectors
       vecta(:)=minim(:,ia)+minim(:,ib)
       vectb(:)=minim(:,ia)-minim(:,ib)
       norm2a=vecta(1)**2+vecta(2)**2+vecta(3)**2
       norm2b=vectb(1)**2+vectb(2)**2+vectb(3)**2
!      Project the trial vector on the first of the two vectors
       reduceda=( minim(1,itrial)*vecta(1)+       &
&       minim(2,itrial)*vecta(2)+       &
&       minim(3,itrial)*vecta(3) )/norm2a
       reducedb=( minim(1,itrial)*vectb(1)+       &
&       minim(2,itrial)*vectb(2)+       &
&       minim(3,itrial)*vectb(3) )/norm2b
       if( abs(abs(reduceda)-0.5d0)<tolsym )then
!        The first vector is an axis
         fact=2 ; center=-1
         cell_base(:,ia)=vecta(:)
         vecta(:)=minim(:,itrial)-reduceda*vecta(:)
         vectb(:)=0.5d0*vectb(:)
         cell_base(:,ib)=vecta(:)+vectb(:)
         cell_base(:,itrial)=vecta(:)-vectb(:)
         call holocell(cell_base,0,foundc,iholohedry,tolsym)
       else if( abs(abs(reducedb)-0.5d0)<tolsym )then
!        The second vector is an axis
         fact=2 ; center=-1
         cell_base(:,ib)=vectb(:)
         vectb(:)=minim(:,itrial)-reducedb*vectb(:)
         vecta(:)=0.5d0*vecta(:)
         cell_base(:,ia)=vectb(:)+vecta(:)
         cell_base(:,itrial)=vectb(:)-vecta(:)
         call holocell(cell_base,0,foundc,iholohedry,tolsym)
       end if
     end if

!    Working hypothesis : the conventional cell is orthogonal,
!    and face centered, in the case where two minimal vectors are equal
     if(foundc==0 .and. orthogonal==1 .and. equal(itrial)==1 ) then
!      Compute the combination of these two vectors
       vecta(:)=minim(:,ia)+minim(:,ib)
       vectb(:)=minim(:,ia)-minim(:,ib)
       norm2a=vecta(1)**2+vecta(2)**2+vecta(3)**2
       norm2b=vectb(1)**2+vectb(2)**2+vectb(3)**2
!      Project the trial vector on the two vectors
       reduceda=( minim(1,itrial)*vecta(1)+       &
&       minim(2,itrial)*vecta(2)+       &
&       minim(3,itrial)*vecta(3) )/norm2a
       reducedb=( minim(1,itrial)*vectb(1)+       &
&       minim(2,itrial)*vectb(2)+       &
&       minim(3,itrial)*vectb(3) )/norm2b
       if( (abs(abs(reduceda)-0.5d0)<tolsym .and. abs(reducedb)<tolsym ) .or. &
&       ( abs(reduceda)<tolsym .and. abs(abs(reducedb)-0.5d0)<tolsym)       )then
         fact=2 ; center=-3
         cell_base(:,itrial)= &
&         (minim(:,itrial)-reduceda*vecta(:)-reducedb*vectb(:) )*2.0d0
         cell_base(:,ia)=vecta(:)
         cell_base(:,ib)=vectb(:)
         call holocell(cell_base,0,foundc,iholohedry,tolsym)
       end if
     end if

!    Working hypothesis : the conventional cell is orthogonal,
!    face centered, but no two vectors are on the same "square"
     if(foundc==0 .and. orthogonal==1)then
!      Compute the combination of these two vectors
       vecta(:)=minim(:,ia)+minim(:,ib)
       vectb(:)=minim(:,ia)-minim(:,ib)
       norm2a=vecta(1)**2+vecta(2)**2+vecta(3)**2
       norm2b=vectb(1)**2+vectb(2)**2+vectb(3)**2
!      The trial vector length must be equal to one of these lengths
       if(abs(metmin(itrial,itrial)-norm2a)<tolsym*norm2a)then
         fact=2 ; center=-3
         cell_base(:,ia)=vecta(:)+minim(:,itrial)
         cell_base(:,ib)=vecta(:)-minim(:,itrial)
!        Project vectb perpendicular to cell_base(:,ia) and cell_base(:,ib)
         norm2a=cell_base(1,ia)**2+cell_base(2,ia)**2+cell_base(3,ia)**2
         norm2b=cell_base(1,ib)**2+cell_base(2,ib)**2+cell_base(3,ib)**2
         reduceda=( cell_base(1,ia)*vectb(1)+       &
&         cell_base(2,ia)*vectb(2)+       &
&         cell_base(3,ia)*vectb(3) )/norm2a
         reducedb=( cell_base(1,ib)*vectb(1)+       &
&         cell_base(2,ib)*vectb(2)+       &
&         cell_base(3,ib)*vectb(3) )/norm2b
         if( abs(abs(reduceda)-0.5d0)<tolsym .and.         &
&         abs(abs(reducedb)-0.5d0)<tolsym      )then
           cell_base(:,itrial)=vectb(:)-reduceda*cell_base(:,ia)-reducedb*cell_base(:,ib)
           call holocell(cell_base,0,foundc,iholohedry,tolsym)
         end if
       else if(abs(metmin(itrial,itrial)-norm2b)<tolsym*norm2b)then
         fact=2 ; center=-3
         cell_base(:,ia)=vectb(:)+minim(:,itrial)
         cell_base(:,ib)=vectb(:)-minim(:,itrial)
!        Project vecta perpendicular to cell_base(:,ia) and cell_base(:,ib)
         norm2a=cell_base(1,ia)**2+cell_base(2,ia)**2+cell_base(3,ia)**2
         norm2b=cell_base(1,ib)**2+cell_base(2,ib)**2+cell_base(3,ib)**2
         reduceda=( cell_base(1,ia)*vecta(1)+       &
&         cell_base(2,ia)*vecta(2)+       &
&         cell_base(3,ia)*vecta(3) )/norm2a
         reducedb=( cell_base(1,ib)*vecta(1)+       &
&         cell_base(2,ib)*vecta(2)+       &
&         cell_base(3,ib)*vecta(3) )/norm2b
         if( abs(abs(reduceda)-0.5d0)<tolsym .and.         &
&         abs(abs(reducedb)-0.5d0)<tolsym      )then
           cell_base(:,itrial)=vecta(:)-reduceda*cell_base(:,ia)-reducedb*cell_base(:,ib)
           call holocell(cell_base,0,foundc,iholohedry,tolsym)
         end if
       end if
     end if

!    Working hypothesis : the cell is rhombohedral, and
!    the three minimal vectors have same length and same absolute scalar product
     if(foundc==0 .and. iholohedry==5 .and. &
&     equal(1)==1 .and. equal(2)==1 .and. equal(3)==1 )then
       if( abs(abs(metmin(1,2))-abs(metmin(1,3)))<tolsym*metmin(1,1) .and.     &
&       abs(abs(metmin(1,2))-abs(metmin(2,3)))<tolsym*metmin(2,2)      )then
         fact=1 ; center=0
         cell_base(:,:)=minim(:,:)
!        One might have to change the sign of one of the vectors
         sign12=1 ; sign13=1 ; sign23=1
         if(metmin(1,2)<0.0d0)sign12=-1
         if(metmin(1,3)<0.0d0)sign13=-1
         if(metmin(2,3)<0.0d0)sign23=-1
         sumsign=sign12+sign13+sign23
         if(sumsign==-1)then
           if(sign12==1)cell_base(:,3)=-cell_base(:,3)
           if(sign13==1)cell_base(:,2)=-cell_base(:,2)
           if(sign23==1)cell_base(:,1)=-cell_base(:,1)
         else if(sumsign==1)then
           if(sign12==-1)cell_base(:,3)=-cell_base(:,3)
           if(sign13==-1)cell_base(:,2)=-cell_base(:,2)
           if(sign23==-1)cell_base(:,1)=-cell_base(:,1)
         end if
         call holocell(cell_base,0,foundc,iholohedry,tolsym)
       end if
     end if

!    DEBUG
!    write(std_out,*)' after test_3a, foundc=',foundc
!    write(std_out,*)' after test_3a, itrial=',itrial
!    write(std_out,*)' after test_3a, equal(:)=',equal(:)
!    ENDDEBUG

!    Working hypothesis : the cell is rhombohedral, one vector
!    is parallel to the trigonal axis
     if(foundc==0 .and. iholohedry==5 .and. equal(itrial)==1 )then
       vecta(:)=minim(:,ia) ; vectb(:)=minim(:,ib)
       norm2trial=minim(1,itrial)**2+minim(2,itrial)**2+minim(3,itrial)**2
       reduceda=( minim(1,itrial)*vecta(1)+       &
&       minim(2,itrial)*vecta(2)+       &
&       minim(3,itrial)*vecta(3) )/norm2trial
       reducedb=( minim(1,itrial)*vectb(1)+       &
&       minim(2,itrial)*vectb(2)+       &
&       minim(3,itrial)*vectb(3) )/norm2trial
!      DEBUG
!      write(std_out,*)' reduceda,reducedb=',reduceda,reducedb
!      ENDDEBUG
       if(abs(abs(reduceda)-1.0d0/3.0d0)<tolsym .and.      &
&       abs(abs(reducedb)-1.0d0/3.0d0)<tolsym      ) then
!        Possibly change of sign to make positive the scalar product with
!        the vector parallel to the trigonal axis
         if(reduceda<zero)vecta(:)=-vecta(:)
         if(reducedb<zero)vectb(:)=-vectb(:)
!        Projection on the orthogonal plane
         vecta(:)=vecta(:)-abs(reduceda)*cell_base(:,itrial)
         vectb(:)=vectb(:)-abs(reducedb)*cell_base(:,itrial)
!        These two vectors should have an angle of 120 degrees
         norm2a=vecta(1)**2+vecta(2)**2+vecta(3)**2
         scalarprod=vecta(1)*vectb(1)+vecta(2)*vectb(2)+vecta(3)*vectb(3)
!        DEBUG
!        write(std_out,*)' norm2a,scalarprod=',norm2a,scalarprod
!        ENDDEBUG
         if(abs(two*scalarprod+norm2a)<tolsym*norm2a)then
           fact=1 ; center=0
           if(scalarprod>0.0d0)vectb(:)=-vectb(:)
!          Now vecta and vectb have an angle of 120 degrees
           cell_base(:,1)=cell_base(:,itrial)/3.0d0+vecta(:)
           cell_base(:,2)=cell_base(:,itrial)/3.0d0+vectb(:)
           cell_base(:,3)=cell_base(:,itrial)/3.0d0-vecta(:)-vectb(:)
!          DEBUG
!          write(std_out,*)' cell_base(:,1)=',cell_base(:,1)
!          write(std_out,*)' cell_base(:,2)=',cell_base(:,2)
!          write(std_out,*)' cell_base(:,3)=',cell_base(:,3)
!          ENDDEBUG
           call holocell(cell_base,0,foundc,iholohedry,tolsym)
         end if
       end if
     end if

!    Working hypothesis : the cell is rhombohedral, one vector
!    is in the plane perpendicular to the trigonal axis
     if(foundc==0 .and. iholohedry==5 .and. equal(itrial)==1 ) then
       vecta(:)=minim(:,ia)+minim(:,ib)
       vectb(:)=minim(:,ia)-minim(:,ib)
       norm2trial=cell_base(1,itrial)**2+cell_base(2,itrial)**2+cell_base(3,itrial)**2
       norm2a=vecta(1)**2+vecta(2)**2+vecta(3)**2
       norm2b=vecta(1)**2+vecta(2)**2+vecta(3)**2
       reduceda=( cell_base(1,itrial)*vecta(1)+       &
&       cell_base(2,itrial)*vecta(2)+       &
&       cell_base(3,itrial)*vecta(3) )/norm2trial
       reducedb=( cell_base(1,itrial)*vectb(1)+       &
&       cell_base(2,itrial)*vectb(2)+       &
&       cell_base(3,itrial)*vectb(3) )/norm2trial
       if(abs(norm2trial-norm2a)<tolsym*norm2a .and. &
&       abs(abs(2*reduceda)-norm2trial)<tolsym*norm2trial    )then
         fact=1 ; center=0
         cell_base(:,1)=minim(:,ia)
         cell_base(:,2)=-minim(:,ib)
         cell_base(:,3)=-minim(:,ib)+2*reduceda*minim(:,itrial)
         call holocell(cell_base,0,foundc,iholohedry,tolsym)
       else if (abs(norm2trial-norm2b)<tolsym*norm2b .and. &
&         abs(abs(2*reducedb)-norm2trial)<tolsym*norm2trial    )then
         fact=1 ; center=0
         cell_base(:,1)=minim(:,ia)
         cell_base(:,2)=minim(:,ib)
         cell_base(:,3)=minim(:,ib)+2*reducedb*minim(:,itrial)
         call holocell(cell_base,0,foundc,iholohedry,tolsym)
       end if
     end if

!    Working hypothesis : the cell is rhombohedral, two vectors
!    are in the plane perpendicular to the trigonal axis
     if(foundc==0 .and. iholohedry==5 .and. equal(itrial)==1 ) then
       vecta(:)=minim(:,ia) ; vectb(:)=minim(:,ib)
       norm2a=vecta(1)**2+vecta(2)**2+vecta(3)**2
       norm2b=vectb(1)**2+vectb(2)**2+vectb(3)**2
       scalarprod=vecta(1)*vectb(1)+vecta(2)*vectb(2)+vecta(3)*vectb(3)
       if(abs(abs(2*scalarprod)-norm2a)<tolsym*norm2a)then
!        This is in order to have 120 angle between vecta and vectb
         if(scalarprod>0.0d0)vectb(:)=-vectb(:)
         reduceda=( cell_base(1,itrial)*vecta(1)+        &
&         cell_base(2,itrial)*vecta(2)+        &
&         cell_base(3,itrial)*vecta(3) )/norm2a
         reducedb=( cell_base(1,itrial)*vectb(1)+        &
&         cell_base(2,itrial)*vectb(2)+        &
&         cell_base(3,itrial)*vectb(3) )/norm2b
         fact=1 ; center=0
         cell_base(:,1)=minim(:,itrial)
         if(abs(reduceda-0.5d0)<tolsym .and. abs(reducedb)<tolsym )then
           cell_base(:,2)=minim(:,itrial)-vecta(:)
           cell_base(:,3)=minim(:,itrial)-vecta(:)-vectb(:)
           call holocell(cell_base,0,foundc,iholohedry,tolsym)
         else if(abs(reduceda-0.5d0)<tolsym.and. abs(reducedb+0.5d0)<tolsym )then
           cell_base(:,2)=minim(:,itrial)-vecta(:)
           cell_base(:,3)=minim(:,itrial)+vectb(:)
           call holocell(cell_base,0,foundc,iholohedry,tolsym)
         else if(abs(reduceda)<tolsym .and. abs(reducedb+0.5d0)<tolsym )then
           cell_base(:,2)=minim(:,itrial)+vectb(:)
           cell_base(:,3)=minim(:,itrial)+vecta(:)+vectb(:)
           call holocell(cell_base,0,foundc,iholohedry,tolsym)
         else if(abs(reduceda+0.5d0)<tolsym .and. abs(reducedb)<tolsym)then
           cell_base(:,2)=minim(:,itrial)+vecta(:)
           cell_base(:,3)=minim(:,itrial)+vecta(:)+vectb(:)
           call holocell(cell_base,0,foundc,iholohedry,tolsym)
         else if(abs(reduceda+0.5d0)<tolsym .and. abs(reducedb-0.5d0)<tolsym)then
           cell_base(:,2)=minim(:,itrial)+vecta(:)
           cell_base(:,3)=minim(:,itrial)-vectb(:)
           call holocell(cell_base,0,foundc,iholohedry,tolsym)
         else if(abs(reduceda)<tolsym .and. abs(reducedb-0.5d0)<tolsym )then
           cell_base(:,2)=minim(:,itrial)-vectb(:)
           cell_base(:,3)=minim(:,itrial)-vecta(:)-vectb(:)
           call holocell(cell_base,0,foundc,iholohedry,tolsym)
         end if
       end if
     end if

!    Working hypothesis : monoclinic holohedry, primitive. Then, two angles are 90 degrees
     if(foundc==0 .and. iholohedry==2 .and. &
&     ang90(ia)==1 .and. ang90(ib)==1 ) then
       fact=1 ; center=0
       cell_base(:,1)=minim(:,ia)
       cell_base(:,2)=minim(:,itrial)
       cell_base(:,3)=minim(:,ib)
!      Checks that the basis vectors are OK for the target holohedry
       call holocell(cell_base,0,foundc,iholohedry,tolsym)
     end if

!    Monoclinic holohedry, one-face-centered cell
!    Working hypothesis, two vectors have equal length.
     do icase=1,5
       if(foundc==0 .and. iholohedry==2 .and. equal(itrial)==1 ) then
         vecta(:)=cell_base(:,ia)+cell_base(:,ib)
         vectb(:)=cell_base(:,ia)-cell_base(:,ib)
         norm2a=vecta(1)**2+vecta(2)**2+vecta(3)**2
         norm2b=vectb(1)**2+vectb(2)**2+vectb(3)**2
!        The minim(:,trial) vector belongs to the
!        plane parallel to the cell_base(:,ia),cell_base(:,ib) plane
!        In that plane, must try minim(:,itrial),
!        as well as the 4 different combinations of
!        minim(:,itrial) with the vectors in the plane
         if(icase==1)vectc(:)=minim(:,itrial)
         if(icase==2)vectc(:)=minim(:,itrial)+cell_base(:,ia)
         if(icase==3)vectc(:)=minim(:,itrial)+cell_base(:,ib)
         if(icase==4)vectc(:)=minim(:,itrial)-cell_base(:,ia)
         if(icase==5)vectc(:)=minim(:,itrial)-cell_base(:,ib)
         norm2c=vectc(1)**2+vectc(2)**2+vectc(3)**2
         sca=vectc(1)*vecta(1)+&
&         vectc(2)*vecta(2)+&
&         vectc(3)*vecta(3)
         scb=vectc(1)*vectb(1)+&
&         vectc(2)*vectb(2)+&
&         vectc(3)*vectb(3)
!        DEBUG
!        write(std_out,*)' symlatt : test iholohedry=2, sca,scb=',sca,scb
!        ENDDEBUG
         if(abs(sca)<tolsym*sqrt(norm2c*norm2a) .or. abs(scb)<tolsym*sqrt(norm2c*norm2b))then
           fact=2 ; center=3
!          The itrial direction is centered
           cell_base(:,3)=vectc(:)
           if(abs(sca)<tolsym*sqrt(norm2c*norm2a))then
             cell_base(:,2)=vecta(:)
             cell_base(:,1)=vectb(:)
             call holocell(cell_base,0,foundc,iholohedry,tolsym)
           else if(abs(scb)<tolsym*sqrt(norm2c*norm2b))then
             cell_base(:,2)=vectb(:)
             cell_base(:,1)=vecta(:)
             call holocell(cell_base,0,foundc,iholohedry,tolsym)
           end if
         end if
       end if
     end do ! icase=1,5

!    Monoclinic holohedry, one-face-centered cell, but non equivalent.
!    This case, one pair of vectors is orthogonal
     if(foundc==0 .and. iholohedry==2 .and. ang90(itrial)==1) then
       vecta(:)=minim(:,ia)
       vectb(:)=minim(:,ib)
       norm2a=vecta(1)**2+vecta(2)**2+vecta(3)**2
       norm2b=vectb(1)**2+vectb(2)**2+vectb(3)**2
!      Project the trial vector on the two vectors
       reduceda=( minim(1,itrial)*vecta(1)+       &
&       minim(2,itrial)*vecta(2)+       &
&       minim(3,itrial)*vecta(3) )/norm2a
       reducedb=( minim(1,itrial)*vectb(1)+       &
&       minim(2,itrial)*vectb(2)+       &
&       minim(3,itrial)*vectb(3) )/norm2b
       if(abs(abs(reduceda)-0.5d0)<tolsym .or. abs(abs(reducedb)-0.5d0)<tolsym) then
         fact=2 ; center=3
         if(abs(abs(reduceda)-0.5d0)<tolsym)then
           cell_base(:,2)=vecta(:)
           cell_base(:,3)=vectb(:)
           cell_base(:,1)=2*(minim(:,itrial)-reduceda*vecta(:))
           call holocell(cell_base,0,foundc,iholohedry,tolsym)
         else if(abs(abs(reducedb)-0.5d0)<tolsym)then
           cell_base(:,2)=vectb(:)
           cell_base(:,3)=vecta(:)
           cell_base(:,1)=2*(minim(:,itrial)-reducedb*vectb(:))
           call holocell(cell_base,0,foundc,iholohedry,tolsym)
         end if
       end if
     end if

!    Monoclinic holohedry, one-face-centered cell, but non equivalent.
!    This case, no pair of vectors is orthogonal, no pair of vector of equal lentgh
     if(foundc==0 .and. iholohedry==2)then
!      Try to find a vector that belongs to the mediator plane, or the binary vector.
!      There must be one such vector, if centered monoclinic and no pair of vectors of equal length,
!      either among the three vectors, or among one of their differences or sums.
!      And there must be, among the two other vectors, one vector whose projection
!      on this vector is half the length of this vector.
       vecta(:)=minim(:,ia)
       vectb(:)=minim(:,ib)
!      Try the different possibilities for the vector on which the projection will be half ...
       do ii=1,5
         if(ii==1)vectc(:)=minim(:,itrial)
         if(ii==2)vectc(:)=minim(:,itrial)+vecta(:)
         if(ii==3)vectc(:)=minim(:,itrial)-vecta(:)
         if(ii==4)vectc(:)=minim(:,itrial)+vectb(:)
         if(ii==5)vectc(:)=minim(:,itrial)-vectb(:)
         norm2trial=vectc(1)**2+vectc(2)**2+vectc(3)**2
!        Project the two vectors on the trial vector
         reduceda=( vectc(1)*vecta(1)+vectc(2)*vecta(2)+vectc(3)*vecta(3) )/norm2trial
         reducedb=( vectc(1)*vectb(1)+vectc(2)*vectb(2)+vectc(3)*vectb(3) )/norm2trial
         found=0
         if(abs(abs(reduceda)-0.5d0)<tolsym)then
           vin1(:)=vectc(:)
           vin2(:)=2.0d0*(vecta(:)-reduceda*vectc(:))
           vext(:)=vectb(:)
           found=1
         else if(abs(abs(reducedb)-0.5d0)<tolsym)then
           vin1(:)=vectc(:)
           vin2(:)=2.0d0*(vectb(:)-reduceda*vectc(:))
           vext(:)=vecta(:)
           found=1
         end if
         if(found==1)exit
       end do
!      Now, vin1 and vin2 are perpendicular to each other, and in the plane that contains the binary vector.
!      One of them must be the binary vector if any.
!      On the other hand, vext is out-of-plane. Might belong to the mediator plane or not.
!      If C monoclinc, then the projection of this vext on the binary vector will be either 0 or +1/2 or -1/2.
!      The binary axis must be stored in cell_base(:,2) for conventional C-cell
       if(found==1)then
         found=0

!        Test vin1 being the binary axis
         norm2trial=vin1(1)**2+vin1(2)**2+vin1(3)**2
         reduceda=(vext(1)*vin1(1)+vext(2)*vin1(2)+vext(3)*vin1(3))/norm2trial
         if(abs(reduceda)<tolsym)then  ! vin1 is the binary axis and vext is in the mediator plane
           found=1
           cell_base(:,1)=vin2(:)
           cell_base(:,2)=vin1(:)
           cell_base(:,3)=vext(:)
         else if(abs(abs(reduceda)-0.5d0)<tolsym)then  ! vin1 is the binary axis and vext has +1/2 or -1/2 as projection
           found=1
           cell_base(:,1)=vin2(:)
           cell_base(:,2)=vin1(:)
           cell_base(:,3)=vext(:)-reduceda*vin1(:)+vin2(:)*half
         else
!          Test vin2 being the binary axis
           norm2trial=vin2(1)**2+vin2(2)**2+vin2(3)**2
           reduceda=(vext(1)*vin2(1)+vext(2)*vin2(2)+vext(3)*vin2(3))/norm2trial
           if(abs(reduceda)<tolsym)then  ! vin2 is the binary axis and vext is in the mediator plane
             found=1
             cell_base(:,1)=vin1(:)
             cell_base(:,2)=vin2(:)
             cell_base(:,3)=vext(:)
           else if(abs(abs(reduceda)-0.5d0)<tolsym)then  ! vin2 is the binary axis and vext has +1/2 or -1/2 as projection
             found=1
             cell_base(:,1)=vin1(:)
             cell_base(:,2)=vin2(:)
             cell_base(:,3)=vext(:)-reduceda*vin2(:)+vin1(:)*half
           end if
         end if

         if(found==1)then
           fact=2 ; center=3
           call holocell(cell_base,0,foundc,iholohedry,tolsym)
         end if
       end if
     end if

   end do ! Do-loop on three different directions
 end do !  Do-loop on different target holohedries

 if(foundc==0)then
   iholohedry=1 ; fact=1 ; center=0
   cell_base(:,:)=minim(:,:)
 end if

!DEBUG
!write(std_out,*)' symlatt : done with centering tests, foundc=',foundc
!write(std_out,*)'  center=',center
!write(std_out,*)'  iholohedry=',iholohedry
!ENDDEBUG

!--------------------------------------------------------------------------
!Final check on the Bravais lattice, using the basis vectors

!Recompute the metric tensor
 if(foundc==1)then
   do ii=1,3
     metmin(:,ii)=cell_base(1,:)*cell_base(1,ii)+&
&     cell_base(2,:)*cell_base(2,ii)+&
&     cell_base(3,:)*cell_base(3,ii)
   end do
 end if

!Examine the angles and vector lengths
 ang90(:)=0
 if(metmin(1,2)**2<tolsym**2*metmin(1,1)*metmin(2,2))ang90(3)=1
 if(metmin(1,3)**2<tolsym**2*metmin(1,1)*metmin(3,3))ang90(2)=1
 if(metmin(2,3)**2<tolsym**2*metmin(2,2)*metmin(3,3))ang90(1)=1
 equal(:)=0
 if(abs(metmin(1,1)-metmin(2,2))<tolsym*half*(metmin(1,1)+metmin(2,2)))equal(3)=1
 if(abs(metmin(1,1)-metmin(3,3))<tolsym*half*(metmin(1,1)+metmin(3,3)))equal(2)=1
 if(abs(metmin(2,2)-metmin(3,3))<tolsym*half*(metmin(2,2)+metmin(3,3)))equal(1)=1

!DEBUG
!write(std_out,*)' symlatt : recompute the  metric tensor '
!write(std_out,*)'  ang90=',ang90
!write(std_out,*)'  equal=',equal
!ENDDEBUG

!The axes will be aligned with the previously determined
!basis vectors, EXCEPT for the tetragonal cell, see later
 axes(:,:)=cell_base(:,:)

 found=0
!Check orthogonal conventional cells
 if(ang90(1)+ang90(2)+ang90(3)==3)then

!  Cubic system
   if(equal(1)+equal(2)+equal(3)==3)then
!    However, one-face centered is not admitted
     if(center==0 .or. center==-1 .or. center==-3)then
       iholohedry=7 ; found=1
       if(center==0)then
         write(message,'(a,a)')ch10,' symlatt: the Bravais lattice is cP (primitive cubic)'
       else if(center==-1)then
         write(message,'(a,a)')ch10,' symlatt: the Bravais lattice is cI (body-centered cubic)'
       else if(center==-3)then
         write(message,'(a,a)')ch10,' symlatt: the Bravais lattice is cF (face-centered cubic)'
       end if
     end if
   end if

!  Tetragonal system
   if(found==0 .and. &
&   (equal(1)==1 .or. equal(2)==1 .or. equal(3)==1) )then
!    However, one-face centered or face-centered is not admitted
     if(center==0 .or. center==-1)then
       iholohedry=4 ; found=1
       if(equal(1)==1)then
         axes(:,3)=cell_base(:,1) ; axes(:,1)=cell_base(:,2) ; axes(:,2)=cell_base(:,3)
       else if(equal(2)==1)then
         axes(:,3)=cell_base(:,2) ; axes(:,2)=cell_base(:,1) ; axes(:,1)=cell_base(:,3)
       else if(equal(3)==1)then
         axes(:,:)=cell_base(:,:)
       end if
       if(center==0)then
         write(message,'(a,a)')ch10,' symlatt: the Bravais lattice is tP (primitive tetragonal)'
       else if(center==-1)then
         write(message,'(a,a)')ch10,' symlatt: the Bravais lattice is tI (body-centered tetragonal)'
       end if
     end if
   end if

!  Orthorhombic system
   if(found==0)then
     iholohedry=3 ; found=1
     axes(:,:)=cell_base(:,:)
     if(center==0)then
       write(message,'(a,a)')ch10,' symlatt: the Bravais lattice is oP (primitive orthorhombic)'
     else if(center==-1)then
       write(message,'(a,a)')ch10,' symlatt: the Bravais lattice is oI (body-centered orthorhombic)'
     else if(center==1 .or. center==2 .or. center==3)then
       write(message,'(a,a)')ch10,' symlatt: the Bravais lattice is oC (one-face-centered orthorhombic)'
     else if(center==-3)then
       write(message,'(a,a)')ch10,' symlatt: the Bravais lattice is oF (face-centered orthorhombic)'
     end if
   end if

 else

!  Hexagonal system
   if(found==0 .and. ang90(1)==1 .and. ang90(2)==1 .and. equal(3)==1 .and. (2*metmin(2,1)+metmin(1,1))<tolsym*metmin(1,1))then
     iholohedry=6 ; found=1
     write(message,'(a,a)')ch10,' symlatt: the Bravais lattice is hP (primitive hexagonal)'
   end if

!  Rhombohedral system
   if(found==0 .and. equal(1)+equal(2)+equal(3)==3 .and.       &
&   abs(metmin(2,1)-metmin(3,2))<tolsym*metmin(2,2)             .and.       &
&   abs(metmin(2,1)-metmin(3,1))<tolsym*metmin(1,1) )then
     iholohedry=5 ; found=1
     write(message,'(a,a)')ch10,' symlatt: the Bravais lattice is hR (rhombohedral)'
   end if

!  Monoclinic system
   if(found==0 .and. ang90(1)+ang90(2)+ang90(3)==2 )then
     iholohedry=2 ; found=1
     if(center==0)then
       write(message,'(a,a)')ch10,' symlatt: the Bravais lattice is mP (primitive monoclinic)'
     else if(center==3)then
       write(message,'(a,a)')ch10,' symlatt: the Bravais lattice is mC (one-face-centered monoclinic)'
     end if
   end if

!  Triclinic system
   if(found==0)then
     iholohedry=1 ; found=1
     write(message,'(a,a)')ch10,' symlatt: the Bravais lattice is aP (primitive triclinic)'
   end if

 end if

 call wrtout(std_out,message,'COLL')

!--------------------------------------------------------------------------
!Make sure that axes form a right-handed coordinate system
!(Note : this should be done in the body of the routine,
!by making changes that leave the sign of the mixed product of the three
!vectors invariant)
 determinant=axes(1,1)*axes(2,2)*axes(3,3) &
& +axes(1,2)*axes(2,3)*axes(3,1) &
& +axes(1,3)*axes(3,2)*axes(2,1) &
& -axes(1,1)*axes(3,2)*axes(2,3) &
& -axes(1,3)*axes(2,2)*axes(3,1) &
& -axes(1,2)*axes(2,1)*axes(3,3)
 if(determinant<0.0d0)then
   axes(:,:)=-axes(:,:)
 end if

!--------------------------------------------------------------------------
!Prefer symmetry axes on the same side as the primitive axes,
!when the changes are allowed
 do
   do ia=1,3
     scprods(ia,:)=axes(1,ia)*rprimd(1,:)+&
&     axes(2,ia)*rprimd(2,:)+&
&     axes(3,ia)*rprimd(3,:)
     norm2trial=sum(axes(:,ia)**2)
     scprods(ia,:)=scprods(ia,:)/sqrt(norm2trial)
   end do
   do ia=1,3
     norm2trial=sum(rprimd(:,ia)**2)
     scprods(:,ia)=scprods(:,ia)/sqrt(norm2trial)
   end do

!  One should now try all the generators of the
!  proper rotations of each Bravais lattice, coupled with change of
!  signs of each vector. This is not done systematically in what follows ...
!  Here, the third axis is left unchanged
   if(iholohedry/=5)then
     if(scprods(1,1)<-tolsym .and. scprods(2,2)<-tolsym)then
       axes(:,1)=-axes(:,1) ; axes(:,2)=-axes(:,2)
       cycle
     end if
   end if
!  The first (or second) axis is left unchanged
   if(iholohedry/=5 .and. iholohedry/=6)then
     if(scprods(2,2)<-tolsym .and. scprods(3,3)<-tolsym)then
       axes(:,2)=-axes(:,2) ; axes(:,3)=-axes(:,3)
       cycle
     end if
     if(scprods(1,1)<-tolsym .and. scprods(3,3)<-tolsym)then
       axes(:,1)=-axes(:,1) ; axes(:,3)=-axes(:,3)
       cycle
     end if
   end if
!  Permutation of the three axis
   if(iholohedry==5 .or. iholohedry==7)then
     trace=scprods(1,1)+scprods(2,2)+scprods(3,3)
     if(trace+tolsym< scprods(1,2)+scprods(2,3)+scprods(3,1))then
       vecta(:)=axes(:,1) ; axes(:,1)=axes(:,3)
       axes(:,3)=axes(:,2); axes(:,2)=vecta(:)
       cycle
     end if
     if(trace+tolsym < scprods(1,3)+scprods(2,1)+scprods(3,2))then
       vecta(:)=axes(:,1) ; axes(:,1)=axes(:,2)
       axes(:,2)=axes(:,3); axes(:,3)=vecta(:)
       cycle
     end if
!    This case is observed when the three new vectors
!    are pointing opposite to the three original vectors
!    One takes their opposite, then switch to of them, then process
!    them again in the loop
     if(sum(scprods(:,:))<tolsym)then
       axes(:,1)=-axes(:,1)
       vecta(:)=-axes(:,2)
       axes(:,2)=-axes(:,3)
       axes(:,3)=vecta(:)
       cycle
     end if
   end if
   exit
!  Other cases might be coded ...
 end do

!--------------------------------------------------------------------------

!DEBUG
!write(std_out,'(a,3es14.6,a,3es14.6,a,3es14.6)')' rprimd=',&
!&  rprimd(:,1),ch10,rprimd(:,2),ch10,rprimd(:,3)
!write(std_out,'(a,3es14.6,a,3es14.6,a,3es14.6)')' axes  =',&
!&  axes(:,1),ch10,axes(:,2),ch10,axes(:,3)
!ENDDEBUG

!Compute the coordinates of rprimd in the system defined by axes(:,:)
 call matr3inv(axes,axesinvt)
 do ii=1,3
   coord(:,ii)=rprimd(1,ii)*axesinvt(1,:)+ &
&   rprimd(2,ii)*axesinvt(2,:)+ &
&   rprimd(3,ii)*axesinvt(3,:)
 end do

!Check that the coordinates are integers, or half-integer in
!the case there is a centering, and generate integer coordinates
 do ii=1,3
   do jj=1,3
     val=coord(ii,jj)*fact
     if(abs(val-nint(val))>fact*two*tolsym)then
       write(message,'(4a,a,3es18.10,a,a,3es18.10,a,a,3es18.10,a,a,i4)')&
&       'One of the coordinates of rprimd in axes is non-integer,',ch10,&
&       'or non-half-integer (if centering), within 2*tolsym.',ch10,&
&       'coord=',coord(:,1),ch10,&
&       '      ',coord(:,2),ch10,&
&       '      ',coord(:,3),ch10,&
&       'fact=',fact
       MSG_BUG(message)
     end if
     icoord(ii,jj)=nint(val)
   end do
 end do

!Store the bravais lattice characteristics
 bravais(1)=iholohedry
 bravais(2)=center
 bravais(3:5)=icoord(1:3,1)
 bravais(6:8)=icoord(1:3,2)
 bravais(9:11)=icoord(1:3,3)

!--------------------------------------------------------------------------
!Initialize the set of symmetries
!Bravais lattices are always invariant under identity and inversion

!Identity and inversion
 ptsymrel(:,:,1)=identity(:,:) ; ptsymrel(:,:,2)=-identity(:,:)
 nptsym=2

!Keep this for IFCv70 compiler
 if(nptsym/=2)then
   write(message,'(a,a,a,a)')ch10,&
&   ' symlatt : BUG -',ch10,&
&   '  Crazy error, compiler bug '
   call wrtout(std_out,message,'COLL')
 end if

!--------------------------------------------------------------------------
!Initialize some generators
!gen6 is defined in a coordinated system with gamma=120 degrees
 gen6(:,:)=0  ; gen6(3,3)=1  ; gen6(1,1)=1  ; gen6(1,2)=-1 ; gen6(2,1)=1
 gen3(:,:)=0  ; gen3(1,2)=1  ; gen3(2,3)=1  ; gen3(3,1)=1
 gen2xy(:,:)=0 ; gen2xy(2,1)=1 ; gen2xy(1,2)=1; gen2xy(3,3)=1
 gen2y(:,:)=0 ; gen2y(1,1)=-1; gen2y(2,2)=1 ; gen2y(3,3)=-1
 gen2z(:,:)=0 ; gen2z(1,1)=-1; gen2z(2,2)=-1; gen2z(3,3)=1

!--------------------------------------------------------------------------

!Define the generators for each holohedry (inversion is already included)
 if(iholohedry==6)then
   ngen=2
   gen(:,:,1)=gen2xy(:,:) ; order(1)=2
   gen(:,:,2)=gen6(:,:)   ; order(2)=6
 else if(iholohedry==5)then
   ngen=2
   gen(:,:,1)=gen2xy(:,:) ; order(1)=2
   gen(:,:,2)=gen3(:,:)   ; order(2)=3
 else
   gen(:,:,1)=gen2y(:,:)  ; order(1)=2
   gen(:,:,2)=gen2z(:,:)  ; order(2)=2
   gen(:,:,3)=gen2xy(:,:) ; order(3)=2
   gen(:,:,4)=gen3(:,:)   ; order(4)=3
   if(iholohedry<=4)ngen=iholohedry-1
   if(iholohedry==7)ngen=4
 end if

!Build the point symmetry operations from generators, in the reduced system
!of coordinates defined by axes(:,:)
 if(ngen/=0)then
   do igen=1,ngen
     do isym=1+nptsym,order(igen)*nptsym
       jsym=isym-nptsym
       do ii=1,3
         ptsymrel(:,ii,isym)=gen(:,1,igen)*ptsymrel(1,ii,jsym)+ &
&         gen(:,2,igen)*ptsymrel(2,ii,jsym)+ &
&         gen(:,3,igen)*ptsymrel(3,ii,jsym)
       end do
     end do
     nptsym=order(igen)*nptsym

   end do
 end if

!--------------------------------------------------------------------------

!Transform symmetry matrices in the system defined by rprimd
 call symrelrot(nptsym,axes,rprimd,ptsymrel,tolsym)

!DEBUG
!write(std_out,'(a)') ' symlatt : exit '
!ENDDEBUG

end subroutine symlatt
!!***

end module m_symfind
!!***
