!{\src2tex{textfont=tt}}
!!****f* ABINIT/symlatt
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
!! COPYRIGHT
!! Copyright (C) 2000-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! msym=default maximal number of symmetries
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
!! of the Bravais lattice in real space in terms
!! of primitive translations.
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
!!      ingeo,inqpt,m_ab7_symmetry,m_effective_potential_file,m_use_ga,symanal
!!      symbrav,thmeig
!!
!! CHILDREN
!!      holocell,matr3inv,smallprim,symrelrot,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symlatt'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry, except_this_one => symlatt
!End of the abilint section

 implicit none

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
         write(message,'(a,a)')ch10,' symlatt : the Bravais lattice is cP (primitive cubic)'
       else if(center==-1)then
         write(message,'(a,a)')ch10,' symlatt : the Bravais lattice is cI (body-centered cubic)'
       else if(center==-3)then
         write(message,'(a,a)')ch10,' symlatt : the Bravais lattice is cF (face-centered cubic)'
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
         write(message,'(a,a)')ch10,' symlatt : the Bravais lattice is tP (primitive tetragonal)'
       else if(center==-1)then
         write(message,'(a,a)')ch10,' symlatt : the Bravais lattice is tI (body-centered tetragonal)'
       end if
     end if
   end if

!  Orthorhombic system
   if(found==0)then
     iholohedry=3 ; found=1
     axes(:,:)=cell_base(:,:)
     if(center==0)then
       write(message,'(a,a)')ch10,' symlatt : the Bravais lattice is oP (primitive orthorhombic)'
     else if(center==-1)then
       write(message,'(a,a)')ch10,' symlatt : the Bravais lattice is oI (body-centered orthorhombic)'
     else if(center==1 .or. center==2 .or. center==3)then
       write(message,'(a,a)')ch10,' symlatt : the Bravais lattice is oC (one-face-centered orthorhombic)'
     else if(center==-3)then
       write(message,'(a,a)')ch10,' symlatt : the Bravais lattice is oF (face-centered orthorhombic)'
     end if
   end if

 else

!  Hexagonal system
   if(found==0 .and. ang90(1)==1 .and. ang90(2)==1 .and. equal(3)==1 .and. (2*metmin(2,1)+metmin(1,1))<tolsym*metmin(1,1))then
     iholohedry=6 ; found=1
     write(message,'(a,a)')ch10,' symlatt : the Bravais lattice is hP (primitive hexagonal)'
   end if

!  Rhombohedral system
   if(found==0 .and. equal(1)+equal(2)+equal(3)==3 .and.       &
&   abs(metmin(2,1)-metmin(3,2))<tolsym*metmin(2,2)             .and.       &
&   abs(metmin(2,1)-metmin(3,1))<tolsym*metmin(1,1) )then
     iholohedry=5 ; found=1
     write(message,'(a,a)')ch10,' symlatt : the Bravais lattice is hR (rhombohedral)'
   end if

!  Monoclinic system
   if(found==0 .and. ang90(1)+ang90(2)+ang90(3)==2 )then
     iholohedry=2 ; found=1
     if(center==0)then
       write(message,'(a,a)')ch10,' symlatt : the Bravais lattice is mP (primitive monoclinic)'
     else if(center==3)then
       write(message,'(a,a)')ch10,' symlatt : the Bravais lattice is mC (one-face-centered monoclinic)'
     end if
   end if

!  Triclinic system
   if(found==0)then
     iholohedry=1 ; found=1
     write(message,'(a,a)')ch10,' symlatt : the Bravais lattice is aP (primitive triclinic)'
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
     if(abs(val-nint(val))>fact*tolsym)then
       write(message,'(4a,a,3es18.10,a,a,3es18.10,a,a,3es18.10,a,a,i4)')&
&       'One of the coordinates of rprimd in axes is non-integer,',ch10,&
&       'or non-half-integer (if centering).',ch10,&
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

end subroutine symlatt
!!***
