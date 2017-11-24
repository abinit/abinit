!{\src2tex{textfont=tt}}
!!****f* ABINIT/symfind
!! NAME
!! symfind
!!
!! FUNCTION
!! Symmetry finder.
!! From the symmetries of the Bravais lattice (ptsymrel),
!! select those that leave invariant the system, and generate
!! the corresponding tnons vectors.
!! The algorithm is explained in T.G. Worlton and J.L. Warren,
!! Comp. Phys. Comm. 3, 88 (1972)
!!
!! COPYRIGHT
!! Copyright (C) 2000-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! berryopt    =  4/14, 6/16, 7/17: electric or displacement field 
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
!! use_inversion=1 if inversion can be included in set of symmetries
!! xred(3,natom)=reduced coordinates of atoms in terms of real space
!!   primitive translations
!!
!! OUTPUT
!! nsym=actual number of symmetries
!! symafm(1:msym)=(anti)ferromagnetic part of nsym symmetry operations
!! symrel(3,3,1:msym)= nsym symmetry operations in real space in terms
!!  of primitive translations
!! tnons(3,1:msym)=nonsymmorphic translations for each symmetry (would
!!  be 0 0 0 each for a symmorphic space group)
!!
!! PARENTS
!!      ingeo,inqpt,m_ab7_symmetry,m_effective_potential_file,m_tdep_sym
!!      m_use_ga,thmeig
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine symfind(berryopt,efield,gprimd,jellslab,msym,natom,noncoll,nptsym,nsym,&
&  nzchempot,ptsymrel,spinat,symafm,symrel,tnons,tolsym,typat,use_inversion,xred,&
&  nucdipmom)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symfind'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: berryopt,jellslab,msym,natom,noncoll,nptsym,nzchempot,use_inversion
 integer,intent(out) :: nsym
 real(dp),intent(in) :: tolsym
!arrays
 integer,intent(in) :: ptsymrel(3,3,msym),typat(natom)
 integer,intent(inout) :: symafm(msym),symrel(3,3,msym) !vz_i
 real(dp),intent(in) :: efield(3),gprimd(3,3),spinat(3,natom),xred(3,natom)
 real(dp),optional :: nucdipmom(3,natom)
 real(dp),intent(inout) :: tnons(3,msym) !vz_i

!Local variables-------------------------------
! TRUE if antiferro symmetries are used with non-collinear magnetism
!scalars
 integer :: found3,foundcl,iatom,iatom0,iatom1,iatom2,iatom3,iclass,iclass0,ii
 integer :: isym,jj,kk,natom0,nclass,ntrial,printed,trialafm,trialok
 integer :: itrialafm
 real(dp) :: spinatcl2,spinatcl20,det
 logical,parameter :: afm_noncoll=.true.
 logical :: test_sameabscollin,test_sameabsnoncoll,test_samespin, test_afmspin_noncoll
 character(len=500) :: message
!arrays
 integer,allocatable :: class(:,:),natomcl(:),typecl(:)
 real(dp) :: diff(3),efieldrot(3),sxred0(3),symnucdipmom2(3)
 real(dp) :: symspinat2(3),symxred2(3),trialnons(3)
 real(dp),allocatable :: spinatcl(:,:),spinatred(:,:)

!**************************************************************************

!DEBUG
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
! call flush(6)
!ENDDEBUG

!Find the number of classes of atoms (type and spinat must be identical,
!spinat might differ by a sign, if aligned with the z direction)
!natomcl(iclass) will contain the number of atoms in the class
!typecl(iclass) will contain the type of the atoms in the class
!spinatcl(1:3,iclass) will contain the spinat of the atoms in the class
!class(1:natomclass(iclass),iclass) will contain the index of the
!atoms belonging to the class
 ABI_ALLOCATE(class,(natom+3,natom))
 ABI_ALLOCATE(natomcl,(natom))
 ABI_ALLOCATE(typecl,(natom))
 ABI_ALLOCATE(spinatcl,(3,natom))

!Initialise with the first atom
 nclass=1
 natomcl(1)=1
 typecl(1)=typat(1)
 spinatcl(:,1)=spinat(:,1)
 class(1,1)=1
 if(natom>1)then
   do iatom=2,natom
!    DEBUG
!    write(std_out,*)' symfind : examine iatom=',iatom
!    ENDDEBUG
     foundcl=0
     do iclass=1,nclass
!      Compare the typat and spinat of atom iatom with existing ones.
!      Admit either identical spinat, or z-aligned spinat with same
!      absolute magnitude
       if( typat(iatom)==typecl(iclass)) then
! spins are vector identical
         test_samespin=  &
&         abs(spinat(1,iatom)-spinatcl(1,iclass))<tolsym .and. &
&         abs(spinat(2,iatom)-spinatcl(2,iclass))<tolsym .and. &
&         abs(spinat(3,iatom)-spinatcl(3,iclass))<tolsym
! spins are vector identical except z component which can have a sign flip
         test_sameabscollin= &
&         noncoll==0 .and.&
&         abs(spinat(1,iatom))<tolsym .and. abs(spinatcl(1,iclass))<tolsym .and.&
&         abs(spinat(2,iatom))<tolsym .and. abs(spinatcl(2,iclass))<tolsym .and.&
&         abs(abs(spinat(3,iatom))-abs(spinatcl(3,iclass)))<tolsym
! spins are vector identical with an AFM sign flip to +
         test_afmspin_noncoll= &
&         noncoll==1 .and. afm_noncoll .and. &
&         abs(spinat(1,iatom)+spinatcl(1,iclass))<tolsym .and. &
&         abs(spinat(2,iatom)+spinatcl(2,iclass))<tolsym .and. &
&         abs(spinat(3,iatom)+spinatcl(3,iclass))<tolsym
! spin vectors have the same norm... 
! TODO: This has to be improved as it is not sufficient to assert they have the same class.
         test_sameabsnoncoll= &
&         noncoll==1 .and. afm_noncoll .and. &
&         ( spinat(1,iatom)**2+spinat(2,iatom)**2+spinat(3,iatom)**2 - &
&          (spinatcl(1,iclass)**2+spinatcl(2,iclass)**2+spinatcl(3,iclass)**2) < tolsym)

         if( test_samespin .or. test_sameabscollin .or. test_afmspin_noncoll .or. test_sameabsnoncoll) then
!          DEBUG
!          write(std_out,*)' symfind : find it belongs to class iclass=',iclass
!          write(std_out,*)' symfind : spinat(:,iatom)=',spinat(:,iatom)
!          write(std_out,*)' symfind : spinatcl(:,iclass)=',spinatcl(:,iclass)
!          write(std_out,*)' symfind : test_samespin,test_sameabscollin,test_afmspin_noncoll,test_sameabsnoncoll=',&
!          &      test_samespin,test_sameabscollin,test_afmspin_noncoll,test_sameabsnoncoll
!          ENDDEBUG
           natomcl(iclass)=natomcl(iclass)+1
           class(natomcl(iclass),iclass)=iatom
           foundcl=1
           exit
         end if
       end if
     end do !loop over iclass
!    If no class with these characteristics exist, create one
     if(foundcl==0)then
       nclass=nclass+1
       natomcl(nclass)=1
       typecl(nclass)=typat(iatom)
       spinatcl(:,nclass)=spinat(:,iatom)
       class(1,nclass)=iatom
     end if
   end do ! loop over atoms
 end if

!DEBUG
!write(std_out,*)' symfind : found ',nclass,' nclass of atoms'
!do iclass=1,nclass
!write(std_out,*)'  class number',iclass
!write(std_out,*)'   natomcl =',natomcl(iclass)
!write(std_out,*)'   typecl  =',typecl(iclass)
!write(std_out,*)'   spinatcl=',spinatcl(:,iclass)
!write(std_out,*)'   class   =',(class(iatom,iclass),iatom=1,natomcl(iclass))
!end do
!ENDDEBUG

!Select the class with the least number of atoms, and non-zero spinat if any
!It is important to select a magnetic class of atom, if any, otherwise
!the determination of the initial (inclusive) set of symmetries takes only
!non-magnetic symmetries, and not both magnetic and non-magnetic ones, see later.
 iclass0=1
 natom0=natomcl(1)
 spinatcl20=spinatcl(1,1)**2+spinatcl(2,1)**2+spinatcl(3,1)**2
 if(nclass>1)then
   do iclass=2,nclass
     spinatcl2=spinatcl(1,iclass)**2+spinatcl(2,iclass)**2+spinatcl(3,iclass)**2
     if( (natomcl(iclass)<natom0 .and. (spinatcl20<tolsym .or. spinatcl2>tolsym))  &
&     .or. (spinatcl20<tolsym .and. spinatcl2>tolsym)                         )then
       iclass0=iclass
       natom0=natomcl(iclass)
       spinatcl20=spinatcl2
     end if
   end do
 end if

 printed=0

!DEBUG
!write(std_out,*)' symfind : has selected iclass0=',iclass0
!write(std_out,*)' #    iatom     xred             spinat '
!do iatom0=1,natomcl(iclass0)
!iatom=class(iatom0,iclass0)
!write(std_out,'(2i4,6f10.4)' )iatom0,iatom,xred(:,iatom),spinat(:,iatom)
!end do
!ENDDEBUG

!If non-collinear spinat have to be used, transfer them in reduced coordinates
 if (noncoll==1) then
   ABI_ALLOCATE(spinatred,(3,natom))
   do iatom=1,natom
     do ii=1,3
       spinatred(ii,iatom)=gprimd(1,ii)*spinat(1,iatom) &
&       +gprimd(2,ii)*spinat(2,iatom) &
&       +gprimd(3,ii)*spinat(3,iatom)
     end do
   end do
 end if

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

!  Select a tentative set of associated translations
!  First compute the symmetric of the first atom in the smallest class,
!  using the point symmetry
   iatom0=class(1,iclass0)
   sxred0(:)=ptsymrel(:,1,isym)*xred(1,iatom0)+ &
&   ptsymrel(:,2,isym)*xred(2,iatom0)+ &
&   ptsymrel(:,3,isym)*xred(3,iatom0)

!  From the set of possible images, deduce tentative translations,
!  and magnetic factor then test whether it sends each atom on a symmetric one
   ntrial=0
   do ii=1,natom0
     iatom1=class(ii,iclass0)

!    The tentative translation is found
     trialnons(:)=xred(:,iatom1)-sxred0(:)
!     trialafm=1
!     if(sum(abs(spinat(:,iatom1)-spinat(:,iatom0))) > tolsym) then
!       trialafm=-1
!     end if

!     if(sum(abs(spinat(:,iatom1)*trialafm-spinat(:,iatom0))) > tolsym)then
! Just impose the norms are equal. The turning spinat is dealt with below
     if(sum(spinat(:,iatom1)**2)-sum(spinat(:,iatom0)**2) > tolsym)then
       write(message,'(3a,3i5)')&
&       'Problem with matching the spin part within a class.',ch10,&
&       'isym,iatom0,iatom1=',isym,iatom0,iatom1
       MSG_ERROR_CLASS(message, "TolSymError")
     end if

!    jellium slab case: check whether symmetry operation has no translational
!    component along z
     if( jellslab/=0 .and. abs(trialnons(3)) > tolsym ) cycle
     trialok=1

!    DEBUG
!    write(std_out,*)' isym,trialnons(:)=',isym,trialnons(:)
!    ENDDEBUG

!    Loop over all classes, then all atoms in the class,
!    to find whether they have a symmetric
     do iclass=1,nclass
       do jj=1,natomcl(iclass)

         iatom2=class(jj,iclass)
!        Generate the tentative symmetric position of iatom2
         symxred2(:)=ptsymrel(:,1,isym)*xred(1,iatom2)+ &
&                    ptsymrel(:,2,isym)*xred(2,iatom2)+ &
&                    ptsymrel(:,3,isym)*xred(3,iatom2)+ trialnons(:)

!        added an explicit loop over FM and AFM - for non-collinear spins you
!        can not tell ahead of time which to choose
!        start with FM to default to that if there is no spin on this atom
         do itrialafm = 1,0,-1
           trialafm = itrialafm*2-1
!          Generate the tentative symmetric spinat of iatom2
           if (noncoll==0) then
             symspinat2(:)=trialafm*spinat(:,iatom2)
           else
             symspinat2(:)=trialafm*(ptsymrel(:,1,isym)*spinatred(1,iatom2)+ &
&                                    ptsymrel(:,2,isym)*spinatred(2,iatom2)+ &
&                                    ptsymrel(:,3,isym)*spinatred(3,iatom2))
           end if
           if(present(nucdipmom)) then
!          Generate the tentative symmetric nuclear dipole moment of iatom2
! TODO: add trialafm? Does it make sense for the nuclear dipole moment, or do
! these phases exist, which are AFM wrt nuclear moments?
             symnucdipmom2(:)=(ptsymrel(:,1,isym)*nucdipmom(1,iatom2)+ &
&                              ptsymrel(:,2,isym)*nucdipmom(2,iatom2)+ &
&                              ptsymrel(:,3,isym)*nucdipmom(3,iatom2))
           end if

!          DEBUG
!          write(std_out,'(a,3i6,3f12.4,a,3f12.4)') ' iclass jj trialafm Send atom at xred=',&
!&            iclass, jj, trialafm, xred(:,iatom2),' to ',symxred2(:)
!          ENDDEBUG

!          Check whether there exists an atom of the same class at the
!          same location, with the correct spinat and nuclear dipole moment
           do kk=1,natomcl(iclass)

             found3=1
             iatom3=class(kk,iclass)
!            Check the location
             diff(:)=xred(:,iatom3)-symxred2(:)
             diff(:)=diff(:)-nint(diff(:))
             if( (diff(1)**2+diff(2)**2+diff(3)**2) > tolsym**2 )found3=0
!            Check the spinat
             if (noncoll==0) then
               diff(:)=spinat(:,iatom3)-symspinat2(:)
             else
               diff(:)=spinatred(:,iatom3)-symspinat2(:)
             end if
             if( (diff(1)**2+diff(2)**2+diff(3)**2) > tolsym**2 )found3=0
             if(present(nucdipmom)) then
!              Check the nuclear dipole moment
               diff(:) = nucdipmom(:,iatom3) - symnucdipmom2(:)
               if(any(diff>tolsym))found3=0
             end if
             
             if(found3==1)exit

!            End loop over iatom3
           end do
           if(found3==1)exit
         end do ! itrialafm

         if(found3==0)then
           trialok=0
           exit
         end if

!        End loop over iatom2
       end do

       if(trialok==0)exit

!      End loop over all classes
     end do

     if(trialok==1)then
       nsym=nsym+1
       if(nsym>msym)then
         write(message,'(3a,i0,4a)')&
&         'The number of symmetries (including non-symmorphic translations)',ch10,&
&         'is larger than maxnsym=',msym,ch10,&
&         'Action: increase maxnsym in the input, or take a cell that is primitive, ',ch10,&
&         'or at least smaller than the present one.'
         MSG_ERROR(message)
       end if
       ntrial=ntrial+1
       symrel(:,:,nsym)=ptsymrel(:,:,isym)
! TODO: fix potential confusion if atoms in a class are both FM and AFM coupled.
! This might lead to arbitrary values for trialafm = +-1
! to be tested
       symafm(nsym)=trialafm
       tnons(:,nsym)=trialnons(:)-nint(trialnons(:)-tolsym)
     end if

!    End the loop on tentative translations
   end do

!  End big loop over each symmetry operation of the Bravais lattice
 end do

 ABI_DEALLOCATE(class)
 ABI_DEALLOCATE(natomcl)
 ABI_DEALLOCATE(spinatcl)
 ABI_DEALLOCATE(typecl)
 if (noncoll==1)   then
   ABI_DEALLOCATE(spinatred)
 end if

!DEBUG
! write(message,'(a,I0,a)')' symfind : exit, nsym=',nsym,ch10
! write(message,'(2a)') trim(message),'   symrel matrices, symafm and tnons are :'
! call wrtout(std_out,message,'COLL')
! do isym=1,nsym
!   write(message,'(i4,4x,3i4,2x,3i4,2x,3i4,4x,i4,4x,3f8.4)' ) isym,symrel(:,:,isym),&
!&   symafm(isym),tnons(:,isym)
!   call wrtout(std_out,message,'COLL')
! end do
!
!stop
!ENDDEBUG

end subroutine symfind
!!***
