!{\src2tex{textfont=tt}}
!!****f* ABINIT/ingeobld
!!
!! NAME
!! ingeobld
!!
!! FUNCTION
!! The geometry builder.
!! Start from the types and coordinates of the primitive atoms
!! and produce the completed set of atoms, by using the definition
!! of objects, then application of rotation, translation and repetition.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! iout=unit number of output file
!! jdtset=number of the dataset looked for
!! lenstr=actual length of the string
!! natrd=number of atoms that have been read in the calling routine
!! natom=number of atoms
!! nobj=the number of objects
!! string*(*)=character string containing all the input data, used
!!  only if choice=1 or 3. Initialized previously in instrng.
!! typat_read(natrd)=type integer for each atom in the primitive set
!! xcart_read(3,natrd)=cartesian coordinates of atoms (bohr), in the primitive set
!!
!! OUTPUT
!! typat(natom)=type integer for each atom in cell
!! xcart(3,natom)=cartesian coordinates of atoms (bohr)
!!
!! PARENTS
!!      ingeo
!!
!! CHILDREN
!!      intagm,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine ingeobld (iout,jdtset,lenstr,natrd,natom,nobj,string,typat,typat_read,xcart,xcart_read)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ingeobld'
 use interfaces_14_hidewrite
 use interfaces_42_parser
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,jdtset,lenstr,natom,natrd,nobj
 character(len=*),intent(in) :: string
!arrays
 integer,intent(in) :: typat_read(natrd)
 integer,intent(out) :: typat(natom)
 real(dp),intent(in) :: xcart_read(3,natrd)
 real(dp),intent(out) :: xcart(3,natom)

!Local variables-------------------------------
 character(len=*), parameter :: format01110 ="(1x,a6,1x,(t9,8i8) )"
 character(len=*), parameter :: format01160 ="(1x,a6,1x,1p,(t9,3g18.10)) "
!scalars
 integer :: belonga,belongb,iatom,iatrd,ii,irep,irep1,irep2,irep3,ivac,marr
 integer :: natom_toberead,nread,objan,objbn,rotate,shift,tread,vacnum
 real(dp) :: angle,cosine,norm2per,norma,normb,normper,project,sine
 character(len=500) :: message
!arrays
 integer :: objarf(3),objbrf(3)
 integer,allocatable :: objaat(:),objbat(:),vaclst(:)
 real(dp) :: axis2(3),axis3(3),axisa(3),axisb(3),objaax(6),objaro(4),objatr(12)
 real(dp) :: objbax(6),objbro(4),objbtr(12),parall(3),perpen(3),rotated(3)
 real(dp) :: vectora(3),vectorb(3)
 real(dp),allocatable :: typat_full(:),xcart_full(:,:)
!no_abirules
!Dummy arguments for subroutine 'intagm' to parse input file
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:)

! *************************************************************************

 marr=max(12,3*natom)
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

!1) Set up the number of vacancies.

!This is the default
 vacnum=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vacnum',tread,'INT')
 if(tread==1) vacnum=intarr(1)

 if (vacnum>0)then
   ABI_ALLOCATE(vaclst,(vacnum))
!  Read list of atoms to be suppressed to create vacancies
   call intagm(dprarr,intarr,jdtset,marr,vacnum,string(1:lenstr),'vaclst',tread,'INT')
   if(tread==1) vaclst(:)=intarr(1:vacnum)
   if(tread/=1)then
     write(message, '(a,a,a,a,a)' )&
&     'The array vaclst MUST be initialized in the input file',ch10,&
&     'when vacnum is non-zero.',ch10,&
&     'Action: initialize vaclst in your input file.'
     MSG_ERROR(message)
   end if
 end if

 natom_toberead=natom+vacnum

!2) Set up list and number of atoms in objects, and the --------------
!operations to be performed on objects.

 write(message,'(80a,a)')('=',ii=1,80),ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(iout,message,'COLL')

 write(message, '(a,a)' )&
& '--ingeobld: echo values of variables connected to objects --------',ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(iout,message,'COLL')

 if(vacnum>0)then
   write(iout,format01110) 'vacnum',vacnum
   write(std_out,format01110) 'vacnum',vacnum
   write(iout,'(1x,a6,1x,(t9,20i3))') 'vaclst',vaclst(:)
   write(std_out,'(1x,a6,1x,(t9,20i3))') 'vaclst',vaclst(:)
   write(iout, '(a)' ) ' '
   write(std_out,'(a)' ) ' '
 end if

 write(iout,format01110) 'nobj',nobj
 write(std_out,format01110) 'nobj',nobj

 if(nobj/=1 .and. nobj/=2)then
   write(message, '(a,a,a,i8,a,a,a)' )&
&   'The number of object (nobj) must be either 1 or 2,',ch10,&
&   'while the input file has  nobj=',nobj,'.',ch10,&
&   'Action: correct nobj in your input file.'
   MSG_ERROR(message)
 end if

 if(nobj==1 .or. nobj==2)then

!  Read the number of atoms of the object a
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'objan',tread,'INT')
   if(tread==1) objan=intarr(1)

   if(tread/=1)then
     write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&     'The number of atoms in object a (objan) must be initialized',ch10,&
&     'in the input file, when nobj=',nobj,'.',ch10,&
&     'This is not the case.',ch10,&
&     'Action: correct objan in your input file.'
     MSG_ERROR(message)
   end if

   write(iout, '(a)' ) ' '
   write(std_out,'(a)' ) ' '
   write(iout,format01110) 'objan',objan
   write(std_out,format01110) 'objan',objan

   if(objan<=1 .or. objan>natom)then
     write(message, '(a,a,a,a,a,i8,a,a,a)' )&
&     'The number of atoms in object a (objan) must be larger than 0',ch10,&
&     'and smaller than natom.',ch10,&
&     'It is equal to ',objan,', an unacceptable value.',ch10,&
&     'Action: correct objan in your input file.'
     MSG_ERROR(message)
   end if

!  Read list of atoms in object a
   call intagm(dprarr,intarr,jdtset,marr,objan,string(1:lenstr),'objaat',tread,'INT')
   ABI_ALLOCATE(objaat,(objan))
   if(tread==1) objaat(1:objan)=intarr(1:objan)

   if(tread/=1)then
     write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&     'The list of atoms in object a (objaat) must be initialized',ch10,&
&     'in the input file, when nobj=',nobj,'.',ch10,&
&     'This is not the case.',ch10,&
&     'Action: initialize objaat in your input file.'
     MSG_ERROR(message)
   end if

   write(iout,'(1x,a6,1x,(t9,20i3))') 'objaat',objaat(:)
   write(std_out,'(1x,a6,1x,(t9,20i3))') 'objaat',objaat(:)

   do iatom=1,objan
     if(objaat(iatom)<1 .or. objaat(iatom)>natom)then
       write(message, '(a,i8,a,a,i8,4a)' )&
&       'The input value of objaat for atom number ',iatom,ch10,&
&       'is equal to ',objaat(iatom),', an unacceptable value :',ch10,&
&       'it should be between 1 and natom. ',&
&       'Action: correct the array objaat in your input file.'
       MSG_ERROR(message)
     end if
   end do

   if(objan>1)then
     do iatom=1,objan-1
       if( objaat(iatom)>=objaat(iatom+1) )then
         write(message, '(a,i8,a,a,a,a,a,a)' )&
&         'The input value of objaat for atom number ',iatom,ch10,&
&         'is larger or equal to the one of the next atom,',ch10,&
&         'while this list should be ordered, and an atom cannot be repeated.',ch10,&
&         'Action: correct the array objaat in your input file.'
         MSG_ERROR(message)
       end if
     end do
   end if

!  Read repetition factors
   objarf(1:3)=1
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'objarf',tread,'INT')
   if(tread==1) objarf(1:3)=intarr(1:3)
   write(iout,'(1x,a6,1x,(t9,20i3))') 'objarf',objarf(:)
   write(std_out,'(1x,a6,1x,(t9,20i3))') 'objarf',objarf(:)

   if(tread==1)then
     do irep=1,3
       if(objarf(irep)<1)then
         write(message, '(a,a,a,3i8,a,a,a)' )&
&         'The input values of objarf(1:3) must be positive,',ch10,&
&         'while it is ',objarf(1:3),'.',ch10,&
&         'Action: correct objarf in your input file.'
         MSG_ERROR(message)
       end if
     end do
   end if

!  Modify the number of atoms to be read
   natom_toberead=natom_toberead-objan*(objarf(1)*objarf(2)*objarf(3)-1)

!  Read rotations angles and translations
   objaro(1:4)=0.0_dp
   objatr(1:12)=0.0_dp
   if (objarf(1)*objarf(2)*objarf(3) ==1) then
     nread=1
   else if (objarf(2)*objarf(3) ==1) then
     nread=2
   else if (objarf(3) ==1) then
     nread=3
   else
     nread=4
   end if
   call intagm(dprarr,intarr,jdtset,marr,nread,string(1:lenstr),'objaro',tread,'DPR')
   if(tread==1) objaro(1:nread)=dprarr(1:nread)

   call intagm(dprarr,intarr,jdtset,marr,3*nread,string(1:lenstr),'objatr',tread,'LEN')

   if(tread==1) objatr(1:3*nread)=dprarr(1:3*nread)
   write(iout,format01160) 'objaro',objaro(1:4)
   write(std_out,format01160) 'objaro',objaro(1:4)
   write(iout,format01160) 'objatr',objatr(1:12)
   write(std_out,format01160) 'objatr',objatr(1:12)
!  If needed, read axes, but default to the x-axis to avoid errors later
   objaax(1:6)=0.0_dp ; objaax(4)=1.0_dp

   if(abs(objaro(1))+abs(objaro(2))+abs(objaro(3))+abs(objaro(4)) > 1.0d-10) then
     call intagm(dprarr,intarr,jdtset,marr,6,string(1:lenstr),'objaax',tread,'LEN')
     if(tread==1) objaax(1:6)=dprarr(1:6)
     if(tread/=1)then
       write(message, '(a,a,a,a,a,a,a)' )&
&       'The axis of object a (objaax) must be initialized',ch10,&
&       'in the input file, when rotations (objaro) are present.',ch10,&
&       'This is not the case.',ch10,&
&       'Action: initialize objaax in your input file.'
       MSG_ERROR(message)
     end if
     write(iout,format01160) 'objaax',objaax(1:6)
     write(std_out,format01160) 'objaax',objaax(1:6)
   end if

   axisa(1:3)=objaax(4:6)-objaax(1:3)
   norma=axisa(1)**2+axisa(2)**2+axisa(3)**2

   if(norma<1.0d-10)then
     write(message, '(5a)' )&
&     'The two points defined by the input array objaax are too',ch10,&
&     'close to each other, and will not be used to define an axis.',ch10,&
&     'Action: correct objaax in your input file.'
     MSG_ERROR(message)
   end if
   axisa(1:3)=axisa(1:3)/sqrt(norma)

!  End condition of existence of a first object
 end if

 if(nobj==2)then

!  Read the number of atoms of the object b
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'objbn',tread,'INT')
   if(tread==1) objbn=intarr(1)

   if(tread/=1)then
     write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&     'The number of atoms in object b (objbn) must be initialized',ch10,&
&     'in the input file, when nobj=',nobj,'.',ch10,&
&     'This is not the case.',ch10,&
&     'Action: initialize objbn in your input file.'
     MSG_ERROR(message)
   end if

   write(iout, '(a)' ) ' '
   write(std_out,'(a)' ) ' '
   write(iout,format01110) 'objbn',objbn
   write(std_out,format01110) 'objbn',objbn

   if(objbn<=1 .or. objbn>natom)then
     write(message, '(a,a,a,a,a,i8,a,a,a)' )&
&     'The number of atoms in object b (objbn) must be larger than 0',ch10,&
&     'and smaller than natom.',ch10,&
&     'It is equal to ',objbn,', an unacceptable value.',ch10,&
&     'Action: correct objbn in your input file.'
     MSG_ERROR(message)
   end if

!  Read list of atoms in object b
   call intagm(dprarr,intarr,jdtset,marr,objbn,string(1:lenstr),'objbat',tread,'INT')
   ABI_ALLOCATE(objbat,(objbn))

   if(tread==1) objbat(1:objbn)=intarr(1:objbn)
   if(tread/=1)then
     write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&     'The list of atoms in object b (objbat) must be initialized',ch10,&
&     'in the input file, when nobj=',nobj,'.',ch10,&
&     'This is not the case.',ch10,&
&     'Action: initialize objbat in your input file.'
     MSG_ERROR(message)
   end if

   write(iout,'(1x,a6,1x,(t9,20i3))') 'objbat',objbat(:)
   write(std_out,'(1x,a6,1x,(t9,20i3))') 'objbat',objbat(:)

   do iatom=1,objbn
     if(objbat(iatom)<1 .or. objbat(iatom)>natom)then
       write(message, '(a,i8,a,a,i8,a,a,a,a,a)' )&
&       'The input value of objbat for atom number ',iatom,ch10,&
&       'is equal to ',objbat(iatom),', an unacceptable value :',ch10,&
&       'it should be between 1 and natom. ',ch10,&
&       'Action: correct objbat in your input file.'
       MSG_ERROR(message)
     end if
   end do

   if(objbn>1)then
     do iatom=1,objbn-1
       if( objbat(iatom)>=objbat(iatom+1) )then
         write(message, '(a,i8,a,a,a,a,a,a)' )&
&         'The input value of objbat for atom number ',iatom,ch10,&
&         'is larger or equal to the one of the next atom,',ch10,&
&         'while this list should be ordered, and an atom cannot be repeated.',ch10,&
&         'Action: correct the array objbat in the input file.'
         MSG_ERROR(message)
       end if
     end do
   end if

!  Read repetition factors
   objbrf(1:3)=1
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'objbrf',tread,'INT')
   if(tread==1) objbrf(1:3)=intarr(1:3)
   write(iout,'(1x,a6,1x,(t9,20i3))') 'objbrf',objbrf(:)
   write(std_out,'(1x,a6,1x,(t9,20i3))') 'objbrf',objbrf(:)

   if(tread==1)then
     do irep=1,3
       if(objbrf(irep)<1)then
         write(message, '(a,a,a,3i8,a,a,a)' )&
&         'The input values of objbrf(1:3) must be positive,',ch10,&
&         'while it is ',objbrf(1:3),'.',ch10,&
&         'Action: correct objbrf in your input file.'
         MSG_ERROR(message)
       end if
     end do
   end if

!  Modify the number of atoms to be read
   natom_toberead=natom_toberead-objbn*(objbrf(1)*objbrf(2)*objbrf(3)-1)
!  Read rotations angles and translations
   objbro(1:4)=0.0_dp
   objbtr(1:12)=0.0_dp
   if (objbrf(1)*objbrf(2)*objbrf(3) ==1) then
     nread=1
   else if (objbrf(2)*objbrf(3) ==1) then
     nread=2
   else if (objbrf(3) ==1) then
     nread=3
   else
     nread=4
   end if
   call intagm(dprarr,intarr,jdtset,marr,nread,string(1:lenstr),'objbro',tread,'DPR')
   if(tread==1) objbro(1:nread)=dprarr(1:nread)

   call intagm(dprarr,intarr,jdtset,marr,3*nread,string(1:lenstr),'objbtr',tread,'LEN')
   if(tread==1) objbtr(1:3*nread)=dprarr(1:3*nread)

   write(iout,format01160) 'objbro',objbro(1:4)
   write(std_out,format01160) 'objbro',objbro(1:4)
   write(iout,format01160) 'objbtr',objbtr(1:12)
   write(std_out,format01160) 'objbtr',objbtr(1:12)

!  If needed, read axes, but default to the x-axis to avoid errors later
   objbax(1:6)=0.0_dp ; objbax(4)=1.0_dp
   if(abs(objbro(1))+abs(objbro(2))+abs(objbro(3))+abs(objbro(4)) > 1.0d-10) then
     call intagm(dprarr,intarr,jdtset,marr,6,string(1:lenstr),'objbax',tread,'LEN')
     if(tread==1) objbax(1:6)=dprarr(1:6)
     if(tread/=1)then
       write(message, '(a,a,a,a,a,a,a)' )&
&       'The axis of object b (objbax) must be initialized',ch10,&
&       'in the input file, when rotations (objbro) are present.',ch10,&
&       'This is not the case.',ch10,&
&       'Action: initialize objbax in your input file.'
       MSG_ERROR(message)
     end if
     write(iout,format01160) 'objbax',objbax(1:6)
     write(std_out,format01160) 'objbax',objbax(1:6)
   end if
   axisb(1:3)=objbax(4:6)-objbax(1:3)
   normb=axisb(1)**2+axisb(2)**2+axisb(3)**2
   if(normb<1.0d-10)then
     write(message, '(5a)' )&
&     'The two points defined by the input array objbax are too',ch10,&
&     'close to each other, and will not be used to define an axis.',ch10,&
&     'Action: correct objbax in your input file.'
     MSG_ERROR(message)
   end if
   axisb(1:3)=axisb(1:3)/sqrt(normb)

!  Check whether both lists are disjoints. Use a very primitive algorithm.
   do iatom=1,objan
     do ii=1,objbn
       if(objaat(iatom)==objbat(ii))then
         write(message, '(6a,i8,a,i8,3a)' )&
&         'The objects a and b cannot have a common atom, but it is',ch10,&
&         'found that the values of objaat and objbat ',&
&         ' are identical, for their',ch10,&
&         'atoms number ',iatom,' and ',ii,'.',ch10,&
&         'Action: change objaat and/or objbat so that they have no common atom anymore.'
         MSG_ERROR(message)
       end if
     end do
   end do

!  End condition of existence of a second object
 end if

!Check whether the number of atoms to be read obtained by relying
!on natom, vacnum and the object definitions, or from natrd coincide
 if(natrd/=natom_toberead)then
   write(message,'(11a,i0,a,i0,2a,i0,a)' )&
&   ' ingeobld : ERROR -- ',ch10,&
&   '  The number of atoms to be read (natrd) must be equal',ch10,&
&   '  to the total number of atoms (natom), plus',ch10,&
&   '  the number of vacancies (vacnum), minus',ch10,&
&   '  the number of atoms added by the repetition of objects.',ch10,&
&   '  This is not the case : natrd= ',natrd,', natom= ',natom,ch10,&
&   ', vacnum= ',vacnum,';'
   call wrtout(std_out,message,"COLL")

   if(nobj==1 .or. nobj==2) then
     write(message,'(a,i3,a,3i3,a,i5,a)' )&
&     '   object a : objan=',objan,', objarf(1:3)=',objarf(1:3),&
&     ' => adds ',objan*(objarf(1)*objarf(2)*objarf(3)-1),' atoms.'
     call wrtout(std_out,message,"COLL")
   end if

   if(nobj==2) then
     write(message,'(a,i3,a,3i3,a,i5,a)' )&
&     '   object b : objbn=',objbn,', objbrf(1:3)=',objbrf(1:3),&
&     ' => adds ',objbn*(objbrf(1)*objbrf(2)*objbrf(3)-1),' atoms.'
     call wrtout(std_out,message,"COLL")
   end if

   write(message,'(3a)' )&
&   '  Action : check the correspondence between natom+vacnum on one side,',ch10,&
&   '           and natrd, objan, objbn, objarf and objbrf on the other side.'
   MSG_ERROR(message)
 end if

!6) Produce full set of atoms

!Print the initial atom coordinates if the geometry builder is used
 write(iout, '(/,a)' )  ' Cartesian coordinates of the primitive atoms '
 write(std_out,'(/,a)' )' Cartesian coordinates of the primitive atoms '
 write(iout,format01160) '      ',xcart_read(:,:)
 write(std_out,format01160) '      ',xcart_read(:,:)

 ABI_ALLOCATE(typat_full,(natom+vacnum))
 ABI_ALLOCATE(xcart_full,(3,natom+vacnum))

!Use the work array xcart_full to produce full set of atoms,
!including those coming from repeated objects.
 iatom=1
 do iatrd=1,natrd

   belonga=0 ; belongb=0
   if(nobj==1 .or. nobj==2)then
!    Determine whether the atom belongs to object a
     do ii=1,objan
       if(iatrd==objaat(ii))belonga=ii
     end do
   end if
   if(nobj==2)then
!    Determine whether the atom belong to object b
     do ii=1,objbn
       if(iatrd==objbat(ii))belongb=ii
     end do
   end if

   write(std_out,'(a,i5,a,i2,i2,a)' ) &
&   ' ingeobld : treating iatrd=',iatrd,', belong(a,b)=',belonga,belongb,'.'

!  In case it does not belong to an object
   if(belonga==0 .and. belongb==0)then
     xcart_full(1:3,iatom)=xcart_read(1:3,iatrd)
     typat_full(iatom)=typat_read(iatrd)
     iatom=iatom+1
   else

!    Repeat, rotate and translate this atom
     if(belonga/=0)then

!      Treat object a
!      Compute the relative coordinate of atom with respect to first point of axis
       vectora(1:3)=xcart_read(1:3,iatrd)-objaax(1:3)
!      Project on axis
       project=vectora(1)*axisa(1)+vectora(2)*axisa(2)+vectora(3)*axisa(3)
!      Get the parallel part
       parall(1:3)=project*axisa(1:3)
!      Get the perpendicular part, to be rotated
       perpen(1:3)=vectora(1:3)-parall(1:3)
!      Compute the norm of the perpendicular part
       norm2per=perpen(1)**2+perpen(2)**2+perpen(3)**2
!      Initialisation to avoid warnings even if used behind if rotate == 1.
       normper = 0
!      It the norm is too small, there is not need to rotate
       rotate=0
       if(norm2per>=1.0d-18)then
         rotate=1
         normper=sqrt(norm2per)
         axis2(1:3)=perpen(1:3)/normper
!        Get the vector perpendicular to axisa and axisa2
         axis3(1)=axisa(2)*axis2(3)-axisa(3)*axis2(2)
         axis3(2)=axisa(3)*axis2(1)-axisa(1)*axis2(3)
         axis3(3)=axisa(1)*axis2(2)-axisa(2)*axis2(1)
       end if

!      Here the repetition loop
       do irep3=1,objarf(3)
         do irep2=1,objarf(2)
           do irep1=1,objarf(1)
!            Here the rotation
             if(rotate==1)then
!              Compute the angle of rotation
               angle=objaro(1)+(irep1-1)*objaro(2)+                     &
&               (irep2-1)*objaro(3)+(irep3-1)*objaro(4)
               cosine=cos(angle/180.0*pi)
               sine=sin(angle/180.0*pi)
               rotated(1:3)=objaax(1:3)+parall(1:3)+&
&               normper*(cosine*axis2(1:3)+sine*axis3(1:3))
             else
               rotated(1:3)=vectora(1:3)
             end if
!            Here the translation
             xcart_full(1:3,iatom)=rotated(1:3)+objatr(1:3)+&
&             (irep1-1)*objatr(4:6)+(irep2-1)*objatr(7:9)+(irep3-1)*objatr(10:12)
             typat_full(iatom)=typat_read(iatrd)
             iatom=iatom+1
           end do
         end do
!        End the repetition loop
       end do

     else
!      If the atom belong to object b
!      Compute the relative coordinate of atom with respect to first point of axis
       vectorb(1:3)=xcart_read(1:3,iatrd)-objbax(1:3)
!      Project on axis
       project=vectorb(1)*axisb(1)+vectorb(2)*axisb(2)+vectorb(3)*axisb(3)
!      Get the parallel part
       parall(1:3)=project*axisb(1:3)
!      Get the perpendicular part, to be rotated
       perpen(1:3)=vectorb(1:3)-parall(1:3)
!      Compute the norm of the perpendicular part
       norm2per=perpen(1)**2+perpen(2)**2+perpen(3)**2
!      Initialisation to avoid warnings even if used behind if rotate == 1.
       normper = 0
!      It the norm is too small, there is not need to rotate
       rotate=0
       if(norm2per>=1.0d-18)then
         rotate=1
         normper=sqrt(norm2per)
         axis2(1:3)=perpen(1:3)/normper
!        Get the vector perpendicular to axisb and axis2
         axis3(1)=axisb(2)*axis2(3)-axisb(3)*axis2(2)
         axis3(2)=axisb(3)*axis2(1)-axisb(1)*axis2(3)
         axis3(3)=axisb(1)*axis2(2)-axisb(2)*axis2(1)
       end if
!      Here the repetition loop
       do irep3=1,objbrf(3)
         do irep2=1,objbrf(2)
           do irep1=1,objbrf(1)
!            Here the rotation
             if(rotate==1)then
!              Compute the angle of rotation
               angle=objbro(1)+(irep1-1)*objbro(2)+                      &
&               (irep2-1)*objbro(3)+ (irep3-1)*objbro(4)
               cosine=cos(angle/180.0*pi)
               sine=sin(angle/180.0*pi)
               rotated(1:3)=objbax(1:3)+parall(1:3)+&
&               normper*(cosine*axis2(1:3)+sine*axis3(1:3))
             else
               rotated(1:3)=vectorb(1:3)
             end if
!            Here the translation
             xcart_full(1:3,iatom)=rotated(1:3)+objbtr(1:3)+&
&             (irep1-1)*objbtr(4:6)+(irep2-1)*objbtr(7:9)+(irep3-1)*objbtr(10:12)
             typat_full(iatom)=typat_read(iatrd)
             iatom=iatom+1
           end do
         end do
!        End the repetition loop
       end do

!      End the condition of belonging to object b
     end if

!    End the condition of belonging to an object
   end if

!  End the loop on atoms
 end do

!Create the vacancies here
 if(vacnum/=0)then
!  First label the vacant atoms as belonging to typat 0
   do ivac=1,vacnum
     typat_full(vaclst(ivac))=0
   end do
!  Then compact the arrays
   shift=0
   do iatom=1,natom
     if(typat_full(iatom+shift)==0) shift=shift+1
     if(shift/=0)then
       xcart_full(1:3,iatom)=xcart_full(1:3,iatom+shift)
       typat_full(iatom)=typat_full(iatom+shift)
     end if
   end do
 end if

!Transfer the content of xcart_full and typat_full to the proper
!location
 xcart(:,1:natom)=xcart_full(:,1:natom)
 typat(1:natom)=typat_full(1:natom)

 ABI_DEALLOCATE(typat_full)
 ABI_DEALLOCATE(xcart_full)
 if(allocated(objaat)) then
   ABI_DEALLOCATE(objaat)
 end if
 if(allocated(objbat)) then
   ABI_DEALLOCATE(objbat)
 end if

 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)
 if (vacnum>0)  then
   ABI_DEALLOCATE(vaclst)
 end if

end subroutine ingeobld
!!***
