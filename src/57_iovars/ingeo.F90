!{\src2tex{textfont=tt}}
!!****f* ABINIT/ingeo
!!
!! NAME
!! ingeo
!!
!! FUNCTION
!! Initialize geometry variables for the ABINIT code.
!! 1) set up unit cell : acell, rprim and rprimd ; deduce Bravais lattice
!! 2) (removed)
!! 3) Set up the number of atoms (natrd) in the primitive set, to be read.
!! 4) Read the type of each atom in the primitive set
!! 5) Read coordinates for each atom in the primitive set
!! 6) Eventually read the symmetries
!! 7) Checks whether the geometry builder must be used,
!!    and call it if needed. Call eventually the symmetry builder and analyser
!!    Make the adequate transfers if the geometry
!!    builder is not needed.
!! 8) Initialize the fixing of atoms,
!!    the initial velocities, and the initial atomic spin
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (XG, RC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! berryopt == 4/14: electric field is on; berryopt = 6/7/16/17: electric displacement field is on
!! iimage= index of the current image
!! iout=unit number of output file
!! jdtset=number of the dataset looked for
!! lenstr=actual length of the string
!! msym=default maximal number of symmetries
!! natom=number of atoms
!! nimage=number of images
!! npsp=number of pseudopotentials (needed for the dimension of znucl)
!! nspden=number of spin-density components
!! nsppol=number of independent spin polarizations
!! ntypat=number of type of atoms
!! nzchempot=defines the use of a spatially-varying chemical potential along z
!! pawspnorb=1 when spin-orbit is activated within PAW
!! ratsph(1:ntypat)=radius of the atomic sphere
!! string*(*)=character string containing all the input data, used
!!  only if choice=1 or 3. Initialized previously in instrng.
!!
!! OUTPUT
!! acell(3)=length of primitive vectors
!! amu(ntypat)=mass of each atomic type
!! bravais(11)=characteristics of Bravais lattice (see symlatt.F90)
!! genafm(3)=magnetic translation generator (in case of Shubnikov group type IV)
!! iatfix(3,natom)=indices for atoms fixed along some (or all) directions
!! jellslab=not zero if jellslab keyword is activated
!! slabzbeg, slabzend= the z coordinates of beginning / end of the jellium slab
!! mixalch(npspalch,ntypalch)=alchemical mixing factors
!! nsym=actual number of symmetries
!! nucdipmom(3,natom)=nuclear magnetic dipole moment of each atom in atomic units
!! ptgroupma = magnetic point group number
!! rprim(3,3)=dimensionless real space primitive translations
!! spgroup=symmetry space group
!! spinat(3,natom)=initial spin of each atom, in unit of hbar/2.
!! symafm(1:msym)=(anti)ferromagnetic part of symmetry operations
!! symmorphi=if 0, only allows symmorphic symmetry operations
!! symrel(3,3,1:msym)=symmetry operations in real space in terms
!!  of primitive translations
!! tnons(3,1:msym)=nonsymmorphic translations for symmetry operations
!! tolsym=tolerance for the symmetry operations
!! typat(natom)=type integer for each atom in cell
!! vel(3,natom)=initial velocity of atoms in bohr/atomic time units
!! vel_cell(3,3)=initial velocity of cell parameters in bohr/atomic time units
!! xred(3,natom)=reduced dimensionless atomic coordinates
!! znucl(1:npsp)=nuclear number of atom as specified in psp file
!!
!! SIDE EFFECTS
!!
!! NOTES
!! the parameters ntypat and natom have already been read in indims,
!! and were used to dimension the arrays needed here.
!!
!! TODO
!! The dtset datastructure should NOT be an argument of this routine ... !
!!
!! MG: I completely agree. Abinit developers must learn that Fortran does not allow for aliasing!
!!
!! PARENTS
!!      invars1
!!
!! CHILDREN
!!      atomdata_from_znucl,chkorthsy,fillcell,gensymshub,gensymshub4
!!      gensymspgr,ingeo_img,ingeobld,intagm,mati3inv,metric,mkradim,mkrdim
!!      randomcellpos,symanal,symatm,symfind,symlatt,symmetrize_rprimd
!!      symmetrize_xred,symrelrot,wrtout,xcart2xred,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine ingeo (acell,amu,dtset,bravais,&
& genafm,iatfix,icoulomb,iimage,iout,jdtset,jellslab,lenstr,mixalch,&
& msym,natom,nimage,npsp,npspalch,nspden,nsppol,nsym,ntypalch,ntypat,&
& nucdipmom,nzchempot,pawspnorb,&
& ptgroupma,ratsph,rprim,slabzbeg,slabzend,spgroup,spinat,string,symafm,&
& symmorphi,symrel,tnons,tolsym,typat,vel,vel_cell,xred,znucl)

 use defs_basis
 use defs_abitypes
 use m_ingeo_img
 use m_profiling_abi
 use m_errors
 use m_atomdata

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ingeo'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_42_parser
 use interfaces_57_iovars, except_this_one => ingeo
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iimage,iout,jdtset,lenstr,msym
 integer,intent(in) :: nimage,npsp,npspalch,nspden,nsppol
 integer,intent(in) :: ntypalch,ntypat,nzchempot,pawspnorb
 integer,intent(inout) :: natom,symmorphi
 integer,intent(out) :: icoulomb,jellslab,ptgroupma,spgroup !vz_i
 integer,intent(inout) :: nsym !vz_i
 real(dp),intent(out) :: slabzbeg,slabzend,tolsym
 character(len=*),intent(in) :: string
!arrays
 integer,intent(out) :: bravais(11),iatfix(3,natom) !vz_i
 integer,intent(inout) :: symafm(msym) !vz_i
 integer,intent(inout) :: symrel(3,3,msym) !vz_i
 integer,intent(out) :: typat(natom)
 real(dp),intent(inout) :: nucdipmom(3,natom)
 real(dp),intent(in) :: ratsph(ntypat)
 real(dp),intent(inout) :: spinat(3,natom)
 real(dp),intent(out) :: acell(3),amu(ntypat),genafm(3),mixalch(npspalch,ntypalch)
 real(dp),intent(inout) :: rprim(3,3),tnons(3,msym) !vz_i
 real(dp),intent(out) :: vel(3,natom),vel_cell(3,3),xred(3,natom)
 real(dp),intent(in) :: znucl(npsp)
 type(dataset_type),intent(in) :: dtset

!Local variables-------------------------------
 character(len=*), parameter :: format01110 ="(1x,a6,1x,(t9,8i8) )"
 character(len=*), parameter :: format01160 ="(1x,a6,1x,1p,(t9,3g18.10)) "
!scalars
 integer :: bckbrvltt,brvltt,chkprim,iatom,idir,iexit,ii,ipsp,irreducible,isym,itypat
 integer :: jsym,marr,mu,natfix,natrd,nobj,noncoll
 integer :: nptsym,nsym_now,ntyppure,random_atpos,shubnikov,spgaxor,spgorig
 integer :: spgroupma,tacell,tangdeg,tgenafm,tnatrd,tread,trprim,tscalecart,tspgroupma
 integer :: txangst,txcart,txred,txrandom,use_inversion
 real(dp) :: amu_default,a2,aa,cc,cosang,ucvol,sumalch
 character(len=500) :: message
 type(atomdata_t) :: atom
!arrays
 integer,allocatable :: ptsymrel(:,:,:),typat_read(:),symrec(:,:,:),indsym(:,:,:)
 integer,allocatable :: intarr(:)
 real(dp) :: angdeg(3), field_xred(3),gmet(3,3),gprimd(3,3),rmet(3,3),rcm(3)
 real(dp) :: rprimd(3,3),rprimd_new(3,3),scalecart(3)
!real(dp) :: tsec(2)
 real(dp),allocatable :: mass_psp(:)
 real(dp),allocatable :: tnons_cart(:,:),xangst_read(:,:)
 real(dp),allocatable :: xcart(:,:),xcart_read(:,:),xred_read(:,:),dprarr(:)

! *************************************************************************

 marr=max(12,3*natom,9*msym)
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

!1) set up unit cell : acell, rprim and rprimd ---------------------

 acell(1:3)=one
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'acell',tacell,'LEN')
 if(tacell==1) acell(1:3)=dprarr(1:3)
 call ingeo_img(acell,iimage,jdtset,lenstr,nimage,3,string,"acell",tacell,'LEN')

 scalecart(1:3)=one
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'scalecart',tscalecart,'LEN')
 if(tscalecart==1) scalecart(1:3)=dprarr(1:3)
 call ingeo_img(scalecart,iimage,jdtset,lenstr,nimage,3,string,"scalecart",tscalecart,'LEN')

!Check that input length scales acell(3) are > 0
 do mu=1,3
   if(acell(mu)<=zero) then
     write(message, '(a,i0,a, 1p,e14.6,a,a,a,a)' )&
&     'Length scale ',mu,' is input as acell=',acell(mu),ch10,&
&     'However, length scales must be > 0 ==> stop',ch10,&
&     'Action: correct acell in input file.'
     MSG_ERROR(message)
   end if
 end do

!Initialize rprim, or read the angles
 tread=0
 call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),'rprim',trprim,'DPR')
 if(trprim==1)rprim(:,:)=reshape( dprarr(1:9) , (/3,3/) )
 call ingeo_img(rprim,iimage,jdtset,lenstr,nimage,3,3,string,"rprim",trprim,'DPR')

!If none of the rprim were read ...
 if(trprim==0)then
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'angdeg',tangdeg,'DPR')
   angdeg(:)=dprarr(1:3)
   call ingeo_img(angdeg,iimage,jdtset,lenstr,nimage,3,string,"angdeg",tangdeg,'DPR')

   if(tangdeg==1)then
     call wrtout(std_out,' ingeo : use angdeg to generate rprim.',"COLL")

!    Check that input angles are positive
     do mu=1,3
       if(angdeg(mu)<=0.0_dp) then
         write(message, '(a,i0,a,1p,e14.6,a,a,a,a)' )&
&         'Angle number ',mu,' is input as angdeg=',angdeg(mu),ch10,&
&         'However, angles must be > 0 ==> stop',ch10,&
&         'Action: correct angdeg in input file.'
         MSG_ERROR(message)
       end if
     end do

!    Check that the sum of angles is smaller than 360 degrees
     if(angdeg(1)+angdeg(2)+angdeg(3)>=360.0_dp) then
       write(message, '(a,a,a,es14.4,a,a,a)' )&
&       'The sum of input angles (angdeg(1:3)) must be lower than 360 degrees',ch10,&
&       'while it is ',angdeg(1)+angdeg(2)+angdeg(3),'.',ch10,&
&       'Action: correct angdeg in input file.'
       MSG_ERROR(message)
     end if

     if( abs(angdeg(1)-angdeg(2))<tol12 .and. &
&     abs(angdeg(2)-angdeg(3))<tol12 .and. &
&     abs(angdeg(1)-90._dp)+abs(angdeg(2)-90._dp)+abs(angdeg(3)-90._dp)>tol12 )then
!      Treat the case of equal angles (except all right angles):
!      generates trigonal symmetry wrt third axis
       cosang=cos(pi*angdeg(1)/180.0_dp)
       a2=2.0_dp/3.0_dp*(1.0_dp-cosang)
       aa=sqrt(a2)
       cc=sqrt(1.0_dp-a2)
       rprim(1,1)=aa        ; rprim(2,1)=0.0_dp                 ; rprim(3,1)=cc
       rprim(1,2)=-0.5_dp*aa ; rprim(2,2)= sqrt(3.0_dp)*0.5_dp*aa ; rprim(3,2)=cc
       rprim(1,3)=-0.5_dp*aa ; rprim(2,3)=-sqrt(3.0_dp)*0.5_dp*aa ; rprim(3,3)=cc
!      DEBUG
!      write(std_out,*)' ingeo : angdeg=',angdeg(1:3)
!      write(std_out,*)' ingeo : aa,cc=',aa,cc
!      ENDDEBUG
     else
!      Treat all the other cases
       rprim(:,:)=0.0_dp
       rprim(1,1)=1.0_dp
       rprim(1,2)=cos(pi*angdeg(3)/180.0_dp)
       rprim(2,2)=sin(pi*angdeg(3)/180.0_dp)
       rprim(1,3)=cos(pi*angdeg(2)/180.0_dp)
       rprim(2,3)=(cos(pi*angdeg(1)/180.0_dp)-rprim(1,2)*rprim(1,3))/rprim(2,2)
       rprim(3,3)=sqrt(1.0_dp-rprim(1,3)**2-rprim(2,3)**2)
     end if

   end if
!  No problem if neither rprim nor angdeg are defined : use default rprim
 end if

!Rescale rprim using scalecart (and set scalecart to one)
 rprim(:,1)=scalecart(:)*rprim(:,1)
 rprim(:,2)=scalecart(:)*rprim(:,2)
 rprim(:,3)=scalecart(:)*rprim(:,3)
 scalecart(:)=one

!Compute different matrices in real and reciprocal space, also checks whether ucvol is positive.
 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 tolsym=tol8
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'tolsym',tread,'DPR')
 if(tread==1) tolsym=dprarr(1)

!Find a tentative Bravais lattice and its point symmetries (might not use them)
!Note that the Bravais lattice might not be the correct one yet (because the
!actual atomic locations might lower the symattry obtained from the lattice parameters only)
 ABI_ALLOCATE(ptsymrel,(3,3,msym))
 call symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)

!3) Possibly, initialize a jellium slab
 jellslab=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'jellslab',tread,'INT')
 if(tread==1) jellslab=intarr(1)

 slabzbeg=zero
 slabzend=zero
 if(jellslab/=0)then
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'slabzbeg',tread,'DPR')
   if(tread==1) slabzbeg=dprarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'slabzend',tread,'DPR')
   if(tread==1) slabzend=dprarr(1)
 end if

!4) Set up the number of atoms in the primitive set, to be read.

!This is the default
 natrd=natom

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natrd',tnatrd,'INT')
 if(tnatrd==1) natrd=intarr(1)

 if(natrd<1 .or. natrd>natom)then
   write(message, '(3a,i0,a,i0,2a,a)' )&
&   'The number of atoms to be read (natrd) must be positive and not bigger than natom.',ch10,&
&   'This is not the case : natrd=',natrd,', natom=',natom,ch10,&
&   'Action: correct natrd or natom in your input file.'
   MSG_ERROR(message)
 end if


!5) Read the type and initial spin of each atom in the primitive set--------

!Check for the use of the old name of this variable
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'type',tread,'INT')
 if(tread==1) then
   write(message,'(a,a,a)')&
&   'The use of the "type" input variable is forbidden since version 4.1 .',ch10,&
&   'Action: replace "type" by "typat".'
   MSG_ERROR(message)
 end if

 ABI_ALLOCATE(typat_read,(natrd))
 typat_read(1)=1

 call intagm(dprarr,intarr,jdtset,marr,natrd,string(1:lenstr),'typat',tread,'INT')

!If not read, try the XYZ data
 if(tread==0)then
   call intagm(dprarr,intarr,jdtset,marr,natrd,string(1:lenstr),'_typat',tread,'INT')
 end if
 if(tread==1) typat_read(1:natrd)=intarr(1:natrd)

 do iatom=1,natrd
   if(typat_read(iatom)<1 .or. typat_read(iatom)>ntypat )then
     write(message,'(a,i0,a,i0,a,a,a,i0,a,a,a)')&
&     'The input type of atom number ',iatom,' is equal to ',typat_read(iatom),',',ch10,&
&     'while it should be between 1 and ntypat= ',ntypat,'.',ch10,&
&     'Action: change either the variable typat or the variable ntypat.'
     MSG_ERROR(message)
   end if
 end do

!6) Read coordinates for each atom in the primitive set--------

 ABI_ALLOCATE(xangst_read,(3,natrd))
 ABI_ALLOCATE(xcart_read,(3,natrd))
 ABI_ALLOCATE(xred_read,(3,natrd))

 random_atpos=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'random_atpos',txrandom,'INT')
 if(txrandom==1) random_atpos=intarr(1)
 if (random_atpos < 0 .or. random_atpos > 5) then
   write(message,'(a,a,a)')&
&   'Random positions is a variable defined between 0 and 5. Error in the input file. ',ch10,&
&   'Action: define one of these in your input file.'
   MSG_ERROR(message)
 end if
!if(nimage/=1 .and. iimage/=1)then
!FIXME : should this be called outside the above end if?
 call randomcellpos(natom,npsp,ntypat,random_atpos,ratsph,rprim,rprimd,typat_read,xred_read(:,1:natrd),znucl,acell)
!This should not be printed if randomcellpos did nothing - it contains garbage. Spurious output anyway
!end if

 call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),'xred',txred,'DPR')
 if(txred==1 .and. txrandom == 0) xred_read(:,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )
 call ingeo_img(xred_read,iimage,jdtset,lenstr,nimage,3,natrd,string,"xred",txred,'DPR')

 call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),'xangst',txangst,'DPR')
 if(txangst==1 .and. txrandom==0) xangst_read(:,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )
 call ingeo_img(xangst_read,iimage,jdtset,lenstr,nimage,3,natrd,string,"xangst",txangst,'DPR')

 call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),'xcart',txcart,'LEN')
 if(txcart==1 .and. txrandom==0)xcart_read(:,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )
 call ingeo_img(xcart_read,iimage,jdtset,lenstr,nimage,3,natrd,string,"xcart",txcart,'LEN')

!Might initialize xred from XYZ file
 if(txred+txcart+txangst+txrandom==0)then
   call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),'_xred',txred,'DPR')
   if(txred==1 .and. txrandom==0) xred_read(:,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )

   call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),'_xangst',txangst,'DPR')
   if(txangst==1 .and. txrandom==0) xangst_read(:,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )
 end if

 if (txred+txcart+txangst+txrandom==0) then
   write(message, '(a,a,a)' )&
&   'Neither xred nor xangst nor xcart are present in input file. ',ch10,&
&   'Action: define one of these in your input file.'
   MSG_ERROR(message)
 end if

 if (txred==1)   write(message, '(a)' ) '  xred   is defined in input file'
 if (txangst==1) write(message, '(a)' ) '  xangst is defined in input file'
 if (txcart ==1) write(message, '(a)' ) '  xcart  is defined in input file'
 if (txrandom ==1) write(message, '(a)' ) '  xred  as random positions in the unit cell'
 if (txrandom ==1) write(message, '(a)' ) '  xcart  are defined from a random distribution '
 call wrtout(std_out,message,'COLL')

 if (txred+txcart+txangst+txrandom>1)then
   write(message, '(3a)' )&
&   'Too many input channels for atomic positions are defined.',ch10,&
&   'Action: choose to define only one of these.'
   MSG_ERROR(message)
 end if

 if(txred==1 .or. txrandom /=0 )then
   call wrtout(std_out,' ingeo : takes atomic coordinates from input array xred ','COLL')
   call xred2xcart(natrd,rprimd,xcart_read,xred_read)
 else
   if(txangst==1)then
     call wrtout(std_out,' ingeo : takes atomic coordinates from input array xangst','COLL')
     xcart_read(:,:)=xangst_read(:,:)/Bohr_Ang
   else
     call wrtout(std_out,' ingeo : takes atomic coordinates from input array xcart','COLL')
   end if
   txred=1
 end if
!At this stage, the cartesian coordinates are known, for the atoms whose coordinates where read.

!Here, allocate the variable that will contain the completed
!sets of xcart, after the use of the geometry builder or the symmetry builder
 ABI_ALLOCATE(xcart,(3,natom))

!7) Eventually read the symmetries

!Take care of the symmetries
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nsym',tread,'INT')

 if(tread==1) nsym=intarr(1)

!Check that nsym is not negative
 if (nsym<0) then
   write(message, '(a,i0,a,a,a,a)' )&
&   'Input nsym must be positive or 0, but was ',nsym,ch10,&
&   'This is not allowed.',ch10,&
&   'Action: correct nsym in your input file.'
   MSG_ERROR(message)
 end if
!Check that nsym is not bigger than msym
 if (nsym>msym) then
   write(message, '(a,i0,a,i0,a,a,a,a,a)' )&
&   'Input nsym=',nsym,' exceeds msym=',msym,'.',ch10,&
&   'This is not allowed.',ch10,&
&   'Action: correct nsym in your input file.'
   MSG_ERROR(message)
 end if

!Read symmorphi
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'symmorphi',tread,'INT')
 if(tread==1) symmorphi=intarr(1)

!Now, read the symmetry operations
 if(nsym>0)then
   call intagm(dprarr,intarr,jdtset,marr,9*nsym,string(1:lenstr),'symrel',tread,'INT')
   if(nsym>1 .and. tread==0)then
     write(message,'(3a)')&
&     'When nsym>1, symrel must be defined in the input file.',ch10,&
&     'Action : either change nsym, or define symrel in your input file.'
     MSG_ERROR(message)
   end if
   if(tread==1) symrel(:,:,1:nsym)=reshape( intarr(1:9*nsym) , (/3,3,nsym/) )

!  Take care of tnons
   tnons(:,1:nsym)=zero
   call intagm(dprarr,intarr,jdtset,marr,3*nsym,string(1:lenstr),'tnons',tread,'DPR')
   if(tread==1) tnons(:,1:nsym)=reshape( dprarr(1:3*nsym) , (/3,nsym/) )

   if(symmorphi==0)then
     do isym=1,nsym
       if(sum(tnons(:,isym)**2)>tol6)then
         write(message, '(5a,i0,a,3f8.4,3a)' )&
&         'When symmorph/=1, the vectors of translation (tnons)',ch10,&
&         'a symmetry operation must vanish.',ch10,&
&         'However, for the symmetry operation number ',isym,', tnons =',tnons(:,isym),'.',ch10,&
&         'Action : either change your list of allowed symmetry operations, or use the symmetry finder (nsym=0).'
         MSG_ERROR(message)
       end if
     end do
   end if

!  Take care of symafm
   call intagm(dprarr,intarr,jdtset,marr,nsym,string(1:lenstr),'symafm',tread,'INT')
   if(tread==1) symafm(1:nsym)=intarr(1:nsym)

 end if


!8) Checks whether the geometry builder must be used, and call it if needed.
!Call the symmetry builder and analyzer if needed.

!At this stage, nsym might still contain the default 0, msym contains the default dtset%maxnsym.
!The cartesian coordinates of the atoms of the primitive set are contained in xcart_read.

 nobj=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nobj',tread,'INT')
 if(tread==1) nobj=intarr(1)

!If there are objects, chkprim will not be used immediately
!But, if there are no objects, but a space group, it will be used directly.
 chkprim=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'chkprim',tread,'INT')
 if(tread==1) chkprim=intarr(1)

 if(nobj/=0)then

!  Spinat is read for each atom, from 1 to natom
   call intagm(dprarr,intarr,jdtset,marr,3*natom,string(1:lenstr),'spinat',tread,'DPR')
   if(tread==1) then
     spinat(1:3,1:natom) = reshape( dprarr(1:3*natom) , (/3,natom/) )
   else if (nspden==4.or.(nspden==2.and.nsppol==1)) then
     write(message, '(5a)' )&
&     'When nspden=4 or (nspden==2 and nsppol==1), the input variable spinat must be',ch10,&
&     'defined in the input file, which is apparently not the case.',ch10,&
&     'Action: define spinat or use nspden=1 in your input file.'
     MSG_ERROR(message)
   end if

!  nucdipmom is read for each atom, from 1 to natom
   call intagm(dprarr,intarr,jdtset,marr,3*natom,string(1:lenstr),'nucdipmom',tread,'DPR')
   if(tread==1) then
     nucdipmom(1:3,1:natom) = reshape( dprarr(1:3*natom) , (/3,natom/) )
   end if

!  Will use the geometry builder
   if(tnatrd/=1 .and. nobj/=0)then
     write(message, '(a,a,a,i0,a,a,a,a,a)' )&
&     'The number of atoms to be read (natrd) must be initialized',ch10,&
&     'in the input file, when nobj= ',nobj,'.',ch10,&
&     'This is not the case.',ch10,&
&     'Action: initialize natrd in your input file.'
     MSG_ERROR(message)
   end if

   if(jellslab/=0)then
     write(message, '(a,i0,3a)' )&
&     'A jellium slab cannot be used when nobj= ',nobj,'.',ch10,&
&     'Action: change one of the input variables jellslab or nobj in your input file.'
     MSG_ERROR(message)
   end if

   call ingeobld (iout,jdtset,lenstr,natrd,natom,nobj,string,typat,typat_read,xcart,xcart_read)

!  Finalize the computation of coordinates : produce xred.
   call xcart2xred(natom,rprimd,xcart,xred)

 else ! nobj==0

!  Spinat is read for each irreducible atom, from 1 to natrd
   call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),'spinat',tread,'DPR')
   if(tread==1)spinat(1:3,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )

!  nucdipmom is read for each irreducible atom, from 1 to natrd
   call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),'nucdipmom',tread,'DPR')
   if(tread==1)nucdipmom(1:3,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )

!  Get xred
   call xcart2xred(natrd,rprimd,xcart_read,xred)

   spgroup=0
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spgroup',tread,'INT')
   if(tread==1) spgroup=intarr(1)

   if(spgroup/=0 .or. nsym/=0)then

     if(jellslab/=0 .and. nsym/=1 .and. spgroup/=1)then
       write(message, '(5a)' )&
&       'For the time being, a jellium slab can only be used',ch10,&
&       'either with the symmetry finder (nsym=0) or with the space group 1 (nsym=1)',ch10,&
&       'Action: change one of the input variables jellslab or nsym or spgroup in your input file.'
       MSG_ERROR(message)
     end if

     if(nzchempot/=0 .and. nsym/=1 .and. spgroup/=1)then
       write(message, '(5a)' )&
&       'For the time being, a spatially-varying chemical potential can only be used',ch10,&
&       'either with the symmetry finder (nsym=0) or with the space group 1 (nsym=1)',ch10,&
&       'Action: change one of the input variables nzchempot or nsym or spgroup in your input file.'
       MSG_ERROR(message)
     end if

     typat(1:natrd)=typat_read(1:natrd)

     if(spgroup/=0 .and. nsym/=0)then
       write(message, '(a,i0,a,a,i0,a,a,a,a,a,a,a,a)' )&
&       'The spatial group number spgroup= ',spgroup,ch10,&
&       'is specified, as well as the number of symmetries nsym= ',nsym,ch10,&
&       'This is not allowed, as you can define the symmetries',ch10,&
&       'either using spgroup OR using nsym, but not both.',ch10,&
&       'Action: modify your input file',ch10,&
&       '(either set spgroup to 0, or nsym to 0)'
       MSG_ERROR(message)
     end if

     brvltt=0

     if(spgroup/=0)then

!      Will generate the spatial group using spgroup
!      Assign default values
       spgaxor=1
       spgorig=1
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'brvltt',tread,'INT')
       if(tread==1) brvltt=intarr(1)
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spgaxor',tread,'INT')
       if(tread==1) spgaxor=intarr(1)
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spgorig',tread,'INT')
       if(tread==1) spgorig=intarr(1)

!      Treat the case of magnetic groups
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spgroupma',tspgroupma,'INT')
       if(tspgroupma==1) spgroupma=intarr(1)
       call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'genafm',tgenafm,'DPR')
       if(tgenafm==1) genafm(1:3)=dprarr(1:3)
       if(tspgroupma/=0 .and. tgenafm/=0)then
         write(message, '(a,i0,a,a,3es9.2,a,a,a,a,a,a,a,a)' )&
&         'The spatial group number spgroupma= ',spgroupma,ch10,&
&         'is specified, as well as the antiferromagnetic generator genafm=',genafm(1:3),ch10,&
&         'This is not allowed, as you can define the magnetic space group',ch10,&
&         'either using spgroupma OR using genafm, but not both.',ch10,&
&         'Action: modify your input file',ch10,&
&         '(either define spgroupma or genafm)'
         MSG_ERROR(message)
       end if

!      TODO : all the symmetry generation operations should be in one big routine

!      If spgroupma is defined, check whether it is consistent
!      with spgroup, determine the Shubnikov type,
!      and, for type IV, find the corresponding genafm
       shubnikov=1
       if(tspgroupma==1)then
         call gensymshub(genafm,spgroup,spgroupma,shubnikov)
       else if(tgenafm==1)then
         shubnikov=4
       end if

!      Generate the spatial group of symmetries in a conventional cell
!      In case of Shubnikov space group type IV, only generate the
!      Fedorov (non-magnetic) group. For Shubnikov type III space group,
!      the magnetic part is generated here.
       bckbrvltt=brvltt
       if(brvltt==-1)brvltt=0
       call gensymspgr(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,spgroup,spgroupma,symafm,symrel,tnons)

!      For shubnikov type IV groups,
!      double the space group, using the antiferromagnetic translation generator
       if(shubnikov==4)then
         call gensymshub4(genafm,msym,nsym,symafm,symrel,tnons)
       end if

!      DEBUG
!      write(std_out,*)' after gensymshub4, nsym =',nsym
!      write(std_out,*)' Describe the different symmetry operations (index,symrel,tnons,symafm)'
!      do ii=1,nsym
!      write(std_out,'(i3,2x,9i3,3es12.2,i3)')ii,symrel(:,:,ii),tnons(:,ii),symafm(ii)
!      end do
!      ENDDEBUG

!      If brvltt was -1 at input, one should now change the conventional cell
!      to a primitive one, if brvltt/=1
       if(bckbrvltt==-1 .and. brvltt/=1)then
!        Will work with rprim only
         rprim(:,:)=rprimd(:,:)
         rprimd_new(:,:)=rprimd(:,:)
         acell(:)=1.0_dp
         select case(brvltt)
         case(5)
           rprimd_new(:,2)=(rprim(:,2)+rprim(:,3))*0.5_dp
           rprimd_new(:,3)=(rprim(:,3)-rprim(:,2))*0.5_dp
         case(6)
           rprimd_new(:,1)=(rprim(:,1)+rprim(:,3))*0.5_dp
           rprimd_new(:,3)=(rprim(:,3)-rprim(:,1))*0.5_dp
         case(4)
           rprimd_new(:,1)=(rprim(:,1)+rprim(:,2))*0.5_dp
           rprimd_new(:,2)=(rprim(:,2)-rprim(:,1))*0.5_dp
         case(3)
           rprimd_new(:,1)=(rprim(:,2)+rprim(:,3))*0.5_dp
           rprimd_new(:,2)=(rprim(:,1)+rprim(:,3))*0.5_dp
           rprimd_new(:,3)=(rprim(:,1)+rprim(:,2))*0.5_dp
         case(2)
           rprimd_new(:,1)=(-rprim(:,1)+rprim(:,2)+rprim(:,3))*0.5_dp
           rprimd_new(:,2)=( rprim(:,1)-rprim(:,2)+rprim(:,3))*0.5_dp
           rprimd_new(:,3)=( rprim(:,1)+rprim(:,2)-rprim(:,3))*0.5_dp
         case(7)
           rprimd_new(:,1)=( rprim(:,1)*2.0_dp+rprim(:,2)+rprim(:,3))/3.0_dp
           rprimd_new(:,2)=(-rprim(:,1)      +rprim(:,2)+rprim(:,3))/3.0_dp
           rprimd_new(:,3)=(-rprim(:,1)-rprim(:,2)*2.0_dp+rprim(:,3))/3.0_dp
         end select
         call symrelrot(nsym,rprimd,rprimd_new,symrel,tolsym)
!        Produce xred in the new system of coordinates
         call xred2xcart(natrd,rprimd,xcart,xred)
         call xcart2xred(natrd,rprimd_new,xcart,xred)
!        Produce tnons in the new system of coordinates
         ABI_ALLOCATE(tnons_cart,(3,nsym))
         call xred2xcart(nsym,rprimd,tnons_cart,tnons)
         call xcart2xred(nsym,rprimd_new,tnons_cart,tnons)
         ABI_DEALLOCATE(tnons_cart)

!        DEBUG
!        write(std_out,*)' after change of coordinates, nsym =',nsym
!        write(std_out,*)' Describe the different symmetry operations (index,symrel,tnons,symafm)'
!        do ii=1,nsym
!        write(std_out,'(i3,2x,9i3,3es12.2,i3)')ii,symrel(:,:,ii),tnons(:,ii),symafm(ii)
!        end do
!        ENDDEBUG

!        Prune the symmetry operations : suppress those with
!        exactly the same point and magnetic part
         nsym_now=1
         do isym=2,nsym
           irreducible=1
           do jsym=1,nsym_now
             if(sum(abs(symrel(:,:,isym)-symrel(:,:,jsym)))==0 .and. symafm(isym)==symafm(jsym)) then
               irreducible=0
               exit
             end if
           end do
           if(irreducible==1)then
             nsym_now=nsym_now+1
             symrel(:,:,nsym_now)=symrel(:,:,isym)
             tnons(:,nsym_now)=tnons(:,isym)
             symafm(nsym_now)=symafm(isym)
           end if
         end do
         nsym=nsym_now
!        Translate tnons in the ]-0.5,0.5] interval
         tnons(:,1:nsym)=tnons(:,1:nsym)-nint(tnons(:,1:nsym)-1.0d-8)

!        DEBUG
!        write(std_out,*)' after reduction, nsym =',nsym
!        write(std_out,*)' Describe the different symmetry operations (index,symrel,tnons,symafm)'
!        do ii=1,nsym
!        write(std_out,'(i3,2x,9i3,3es12.2,i3)')ii,symrel(:,:,ii),tnons(:,ii),symafm(ii)
!        end do
!        ENDDEBUG

!        Now that symrel, tnons and xred are expressed in the primitive
!        axis system, update the geometric quantities
         rprimd(:,:)=rprimd_new(:,:)
         rprim(:,:)=rprimd_new(:,:)
         call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
         call symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)
       end if

     end if

     if(natom/=natrd)then
!      Generate the full set of atoms from its knowledge in the irreducible part.
       call fillcell(natom,natrd,nsym,nucdipmom,spinat,symafm,symrel,tnons,tolsym,typat,xred)
     end if

!    Check whether the symmetry operations are consistent with the lattice vectors
     iexit=0

     call chkorthsy(gprimd,iexit,nsym,rmet,rprimd,symrel)

   else ! spgroup==0 and nsym==0

!    Here, spgroup==0 as well as nsym==0, so must generate
!    the spatial group of symmetry. However, all the atom
!    positions must be known, so the number
!    of atoms to be read must equal the total number of atoms.
     if(natrd/=natom)then

       write(message, '(a,i0,a,a,i0,a,a,a,a,a,a,a,a,a)' )&
&       'The number of atoms to be read (natrd)= ',natrd,ch10,&
&       'differs from the total number of atoms (natom)= ',natom,ch10,&
&       'while spgroup=0 and nsym=0.',&
&       'This is not allowed, since the information needed to',ch10,&
&       'generate the missing atomic coordinates is not available.',ch10,&
&       'Action: modify your input file',ch10,&
&       '(either natrd, or natom, or spgroup, or nsym)'
       MSG_ERROR(message)
     else

       typat(:)=typat_read(:)
!      Find the symmetry operations : nsym, symafm, symrel and tnons.
!      Use nptsym and ptsymrel, as determined by symlatt
       noncoll=0;if (nspden==4) noncoll=1
       use_inversion=1;if (nspden==4.or.pawspnorb>0) use_inversion=0
       !use_inversion=1

!      Get field in reduced coordinates (reduced e/d field)

       field_xred(:)=zero
       if (dtset%berryopt ==4) then
         do ii=1,3
           field_xred(ii)=dot_product(dtset%efield(:),gprimd(:,ii))
         end do
       else if (dtset%berryopt == 6 ) then
         do ii=1,3
           field_xred(ii)=dot_product(dtset%dfield(:),gprimd(:,ii))
           field_xred(ii)=field_xred(ii)+ dot_product(dtset%efield(:),gprimd(:,ii)) ! note: symmetry broken by D and E
         end do
       else if (dtset%berryopt == 14) then
         do ii=1,3
           field_xred(ii)=dot_product(dtset%red_efieldbar(:),gmet(:,ii))
         end do
       else if (dtset%berryopt == 16) then
         do ii=1,3
           field_xred(ii)=dtset%red_dfield(ii)+dtset%red_efield(ii)  ! symmetry broken by reduced d and e
         end do
       else if (dtset%berryopt == 17) then
         do ii=1,3
           field_xred(ii)=dot_product(dtset%red_efieldbar(:),gmet(:,ii))
           if(dtset%jfielddir(ii)==2) field_xred(ii)=dtset%red_dfield(ii)
         end do
       end if


       call symfind(dtset%berryopt,field_xred,gprimd,jellslab,msym,natom,noncoll,nptsym,nsym,&
&       nzchempot,ptsymrel,spinat,symafm,symrel,tnons,tolsym,typat,use_inversion,xred,nucdipmom)

!      If the tolerance on symmetries is bigger than 1.e-8, symmetrize the atomic positions
       if(tolsym>1.00001e-8)then
         ABI_ALLOCATE(indsym,(4,natom,nsym))
         ABI_ALLOCATE(symrec,(3,3,nsym))
         do isym=1,nsym
           call mati3inv(symrel(:,:,isym),symrec(:,:,isym))
         end do
         call symatm(indsym,natom,nsym,symrec,tnons,tolsym,typat,xred)
         call symmetrize_xred(indsym,natom,nsym,symrel,tnons,xred)
         ABI_DEALLOCATE(indsym)
         ABI_DEALLOCATE(symrec)

         write(message,'(a,es14.6,10a)')&
&         'The tolerance on symmetries =',tolsym,ch10,&
&         'is bigger than the usual tolerance, i.e. 1.0e-8 .',ch10,&
&         'In order to avoid spurious effect, the atomic coordinates have been',ch10,&
&         'symmetrized before storing them in the dataset internal variable.',ch10,&
&         'So, do not be surprised by the fact that your input variables (xcart, xred, ...)',ch10,&
&         'do not correspond to the ones echoed by ABINIT, the latter being used to do the calculations.'
         MSG_WARNING(message)
       end if

     end if

   end if

!  Finalize the computation of coordinates : produce xcart
   call xred2xcart(natom,rprimd,xcart,xred)

 end if ! check of existence of an object

 ABI_DEALLOCATE(ptsymrel)
 ABI_DEALLOCATE(xangst_read)
 ABI_DEALLOCATE(xcart_read)
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(xred_read)
 ABI_DEALLOCATE(typat_read)

!Correct the default nsym value, if a symmetry group has not been generated.
 if(nsym==0)nsym=1

!--------------------------------------------------------------------------------------------------------

 icoulomb=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'icoulomb',tread,'INT')
 if(tread==1)icoulomb=intarr(1)

!calculate the center of the atomic system such as to put the
!atoms in the middle of the simulation box for the free BC case.
 if (icoulomb == 1) then
   rcm(:)=zero
   do iatom=1,natom
     rcm(:)=rcm(:)+xred(:,iatom)
   end do
   rcm(:)=rcm(:)/real(natom,dp)-half
   do iatom=1,natom
     xred(:,iatom)=xred(:,iatom)-rcm(:)
   end do
!  Also modify the tnons
   do isym=1,nsym
     tnons(:,isym)=matmul(symrel(:,:,isym),rcm(:))-rcm(:)+tnons(:,isym)
   end do

   message = ' Because icoulomb is 1, the average center of coordinates of the system has been translated to (0.5,0.5,0.5) '
   MSG_WARNING(message)
 end if

!========================================================================================================
!
!At this stage, the cell parameters and atomic coordinates are known, as well as the symmetry operations
!There has been a preliminary analysis of the holohedry (not definitive, though ...)
!
!========================================================================================================

!Here, determine correctly the Bravais lattice and other space group or shubnikov group characteristics
 call symanal(bravais,chkprim,genafm,msym,nsym,ptgroupma,rprimd,spgroup,symafm,symrel,tnons,tolsym)

!If the tolerance on symmetries is bigger than 1.e-8, symmetrize the rprimd. Keep xred fixed.
 if(tolsym>1.00001e-8)then
!  Check whether the symmetry operations are consistent with the lattice vectors
   iexit=1

   call chkorthsy(gprimd,iexit,nsym,rmet,rprimd,symrel)
   if(iexit==-1)then
     call symmetrize_rprimd(bravais,nsym,rprimd,symrel,tolsym)
     call mkradim(acell,rprim,rprimd)
     write(message,'(a,es14.6,10a)')&
&     'The tolerance on symmetries =',tolsym,ch10,&
&     'is bigger than the usual tolerance, i.e. 1.0e-8 .',ch10,&
&     'In order to avoid spurious effect, the primitive vectors have been',ch10,&
&     'symmetrized before storing them in the dataset internal variable.',ch10,&
&     'So, do not be surprised by the fact that your input variables (acell, rprim, xcart, xred, ...)',ch10,&
&     'do not correspond to the ones echoed by ABINIT, the latter being used to do the calculations.'
     MSG_WARNING(message)
   end if

 end if

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 angdeg(1)=180.0_dp/pi * acos(rmet(2,3)/sqrt(rmet(2,2)*rmet(3,3)))
 angdeg(2)=180.0_dp/pi * acos(rmet(1,3)/sqrt(rmet(1,1)*rmet(3,3)))
 angdeg(3)=180.0_dp/pi * acos(rmet(1,2)/sqrt(rmet(1,1)*rmet(2,2)))
 write(std_out,'(a,3f14.8)') ' ingeo : angdeg(1:3)=',angdeg(1:3)

!--------------------------------------------------------------------------------------

!Finally prune the set of symmetry in case non-symmorphic operations must be excluded
 if(symmorphi==0)then
   jsym=0
   do isym=1,nsym
     if(sum(tnons(:,isym)**2)<tol6)then
       jsym=jsym+1
!      This symmetry operation is non-symmorphic, and can be kept
       if(isym/=jsym)then
         symrel(:,:,jsym)=symrel(:,:,isym)
         tnons(:,jsym)=tnons(:,isym)
         symafm(jsym)=symafm(isym)
       end if
     end if
   end do
   nsym=jsym
 end if


!DEBUG
!call symmultsg(nsym,symafm,symrel,tnons)
!ENDDEBUG

!9) initialize the list of fixed atoms, and initial velocities -----------------
!Note : these inputs do not influence the previous generation of
!symmetry operations. This might be changed in the future

!idir=0 is for iatfix , idir=1 is for iatfixx,
!idir=2 is for iatfixy, idir=3 is for iatfixz
 iatfix(:,:)=0

 do idir=0,3

   if(idir==0)then
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natfix',tread,'INT')
   else if(idir==1)then
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natfixx',tread,'INT')
   else if(idir==2)then
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natfixy',tread,'INT')
   else if(idir==3)then
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natfixz',tread,'INT')
   end if

!  Use natfix also for natfixx,natfixy,natfixz
   natfix=0
   if(tread==1) natfix=intarr(1)


!  Checks the validity of natfix
   if (natfix<0 .or. natfix>natom) then
     write(message, '(a,a,a,i0,a,i4,a,a,a)' )&
&     'The input variables natfix, natfixx, natfixy and natfixz must be',ch10,&
&     'between 0 and natom (= ',natom,'), while one of them is ',natfix,'.',ch10,&
&     'Action: correct that occurence in your input file.'
     MSG_ERROR(message)
   end if

!  Read iatfix
   if(idir==0)then
     call intagm(dprarr,intarr,jdtset,marr,natfix,string(1:lenstr),'iatfix',tread,'INT')
   else if(idir==1)then
     call intagm(dprarr,intarr,jdtset,marr,natfix,string(1:lenstr),'iatfixx',tread,'INT')
   else if(idir==2)then
     call intagm(dprarr,intarr,jdtset,marr,natfix,string(1:lenstr),'iatfixy',tread,'INT')
   else if(idir==3)then
     call intagm(dprarr,intarr,jdtset,marr,natfix,string(1:lenstr),'iatfixz',tread,'INT')
   end if

!  If some iatfix was read, natfix must vanish
   if (natfix==0 .and. tread==1)then
     write(message, '(a,i1,5a)' )&
&     'For direction ',idir,' the corresponding natfix is zero,',ch10,&
&     'while iatfix specifies some atoms to be fixed.',ch10,&
&     'Action: either specify a non-zero natfix(x,y,z) or suppress iatfix(x,y,z).'
     MSG_ERROR(message)
   end if

!  If natfix is non-zero, iatfix must be defined
   if (natfix>0 .and. tread==0)then
     write(message, '(a,i1,3a,i0,3a)' )&
&     'For direction ',idir,' no iatfix has been specified,',ch10,&
&     'while natfix specifies that some atoms to be fixed, natfix= ',natfix,'.',ch10,&
&     'Action: either set natfix(x,y,z) to zero or define iatfix(x,y,z).'
     MSG_ERROR(message)
   end if

   if(tread==1)then
     do ii=1,natfix
!      Checks the validity of the input iatfix
       if (intarr(ii)<1 .or. intarr(ii)>natom) then
         write(message, '(a,a,a,i0,a,a,a)' )&
&         'The input variables iatfix, iatfixx, iatfixy and iatfixz must be',ch10,&
&         'between 1 and natom, while one of them is ',intarr(ii),'.',ch10,&
&         'Action: correct that occurence in your input file.'
         MSG_ERROR(message)
       end if
!      Finally set the value of the internal iatfix array
       do iatom=1,natom
         if(intarr(ii)==iatom)then
           if(idir==0)iatfix(1:3,iatom)=1
           if(idir/=0)iatfix(idir,iatom)=1
         end if
       end do
     end do
   end if

 end do

 vel(:,:)=zero
 call intagm(dprarr,intarr,jdtset,marr,3*natom,string(1:lenstr),'vel',tread,'DPR')
 if(tread==1)vel(:,:)=reshape( dprarr(1:3*natom) , (/3,natom/) )
 call ingeo_img(vel,iimage,jdtset,lenstr,nimage,3,natom,string,"vel",tread,'DPR')

 vel_cell(:,:)=zero
 call intagm(dprarr,intarr,jdtset,marr,3*3,string(1:lenstr),'vel_cell',tread,'DPR')
 if(tread==1)vel_cell(:,:)=reshape( dprarr(1:9) , (/3,3/) )

!mixalch
 if(ntypalch>0)then
   call intagm(dprarr,intarr,jdtset,marr,npspalch*ntypalch,string(1:lenstr),'mixalch',tread,'DPR')
   if(tread==1) mixalch(1:npspalch,1:ntypalch)=&
&   reshape(dprarr(1:npspalch*ntypalch),(/npspalch,ntypalch/))
   do itypat=1,ntypalch
     sumalch=sum(mixalch(1:npspalch,itypat))
     if(abs(sumalch-one)>tol10)then
       write(message, '(a,i0,2a,f8.2,4a)' )&
&       'For the alchemical atom number ',itypat,ch10,&
&       'the sum of the pseudopotential coefficients is',sumalch,ch10,&
&       'while it should be one.',ch10,&
&       'Action : check the content of the input variable mixalch.'
       MSG_ERROR(message)
     end if
   end do
   call ingeo_img(mixalch,iimage,jdtset,lenstr,nimage,npspalch,ntypalch,string,"mixalch",tread,'DPR')
 end if

!amu (needs mixalch to be initialized ...)
!Find the default mass
 ABI_ALLOCATE(mass_psp,(npsp))
 do ipsp=1,npsp
   call atomdata_from_znucl(atom,znucl(ipsp))
   amu_default = atom%amu
   mass_psp(ipsp)=amu_default
 end do
!When the pseudo-atom is pure, simple copy
 ntyppure=ntypat-ntypalch
 if(ntyppure>0)then
   amu(1:ntyppure)=mass_psp(1:ntyppure)
 end if
!When the pseudo-atom is alchemical, must make mixing
 if(ntypalch>0)then
   do itypat=ntyppure+1,ntypat
     amu(itypat)=zero
     do ipsp=ntyppure+1,npsp
       amu(itypat)=amu(itypat)+mixalch(ipsp-ntyppure,itypat-ntyppure)*mass_psp(ipsp)
     end do
   end do
 end if
 ABI_DEALLOCATE(mass_psp)

 call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),'amu',tread,'DPR')
 if(tread==1)amu(:)=dprarr(1:ntypat)
 call ingeo_img(amu,iimage,jdtset,lenstr,nimage,ntypat,string,"amu",tread,'DPR')


 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)

end subroutine ingeo
!!***
