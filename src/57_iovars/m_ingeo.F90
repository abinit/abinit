!!****m* ABINIT/m_ingeo
!! NAME
!!  m_ingeo
!!
!! FUNCTION
!! Initialize geometry variables for the ABINIT code.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (XG, RC)
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

module m_ingeo

 use defs_basis
 use m_intagm_img
 use m_abicore
 use m_errors
 use m_atomdata
 use m_sort
 use m_dtset

 use m_symtk,      only : mati3inv, chkorthsy, symrelrot, mati3det, symmetrize_rprimd, symmetrize_xred, symatm
 use m_spgbuilder, only : gensymspgr, gensymshub, gensymshub4
 use m_symfind,    only : symfind, symanal, symlatt
 use m_geometry,   only : mkradim, mkrdim, xcart2xred, xred2xcart, randomcellpos, metric
 use m_parser,     only : intagm

 implicit none

 private
!!***

 public :: ingeo        ! Initialize geometry variables for the ABINIT code.
 public :: invacuum     ! Determine whether there is vacuum along some of the primitive directions
!!***

contains
!!***

!!****f* m_ingeo/ingeo
!!
!! NAME
!! ingeo
!!
!! FUNCTION
!! Initialize geometry variables for the ABINIT code.
!! 1) set up unit cell: acell, rprim and rprimd ; deduce Bravais lattice
!! 2) (removed)
!! 3) Set up the number of atoms (natrd) in the primitive set, to be read.
!! 4) Read the type of each atom in the primitive set
!! 5) Read coordinates for each atom in the primitive set
!! 6) Eventually read the symmetries
!! 7) Checks whether the geometry builder must be used,
!!    and call it if needed. Call eventually the symmetry builder and analyser
!!    Make the adequate transfers if the geometry
!!    builder is not needed.
!! 8) Initialize the fixing of atoms, the initial velocities, and the initial atomic spin
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
!! supercell_latt(3,3)=supercell lattice
!!
!! OUTPUT
!! acell(3)=length of primitive vectors
!! amu(ntypat)=mass of each atomic type
!! bravais(11)=characteristics of Bravais lattice (see symlatt.F90)
!! chrgat(natom)=target charge for each atom. Not always used, it depends on the value of constraint_kind
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
!!      gensymspgr,intagm_img,ingeobld,intagm,mati3inv,metric,mkradim,mkrdim
!!      randomcellpos,symanal,symatm,symfind,symlatt,symmetrize_rprimd
!!      symmetrize_xred,symrelrot,wrtout,xcart2xred,xred2xcart
!!
!! SOURCE

subroutine ingeo (acell,amu,bravais,chrgat,dtset,&
& genafm,iatfix,icoulomb,iimage,iout,jdtset,jellslab,lenstr,mixalch,&
& msym,natom,nimage,npsp,npspalch,nspden,nsppol,nsym,ntypalch,ntypat,&
& nucdipmom,nzchempot,pawspnorb,&
& ptgroupma,ratsph,rprim,slabzbeg,slabzend,spgroup,spinat,string,supercell_lattice,symafm,&
& symmorphi,symrel,tnons,tolsym,typat,vel,vel_cell,xred,znucl)

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
 integer,intent(in) :: supercell_lattice(3,3)
 integer,intent(out) :: bravais(11),iatfix(3,natom) !vz_i
 integer,intent(inout) :: symafm(msym) !vz_i
 integer,intent(inout) :: symrel(3,3,msym) !vz_i
 integer,intent(out) :: typat(natom)
 real(dp),intent(inout) :: chrgat(natom)
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
 integer :: bckbrvltt,brvltt,chkprim,i1,i2,i3,iatom,iatom_supercell,idir,iexit,ii
 integer :: ipsp,irreducible,isym,itypat,jsym,marr,mu,multiplicity,natom_uc,natfix,natrd
 integer :: nobj,noncoll,nptsym,nsym_now,ntyppure,random_atpos,shubnikov,spgaxor,spgorig
 integer :: spgroupma,tacell,tangdeg,tgenafm,tnatrd,tread,trprim,tscalecart,tspgroupma
 integer :: txangst,txcart,txred,txrandom,use_inversion
 real(dp) :: amu_default,a2,aa,cc,cosang,ucvol,sumalch
 character(len=500) :: message
 type(atomdata_t) :: atom
!arrays
 integer,allocatable :: ptsymrel(:,:,:),typat_read(:),symrec(:,:,:),indsym(:,:,:)
 integer,allocatable :: intarr(:)
 real(dp) :: angdeg(3), field_xred(3),gmet(3,3),gprimd(3,3),rmet(3,3),rcm(3)
 real(dp) :: rprimd(3,3),rprimd_read(3,3),rprimd_new(3,3),scalecart(3)
 real(dp),allocatable :: mass_psp(:)
 real(dp),allocatable :: tnons_cart(:,:),xangst_read(:,:)
 real(dp),allocatable :: xcart(:,:),xcart_read(:,:),xred_read(:,:),dprarr(:)

! *************************************************************************

 marr=max(12,3*natom,9*msym)
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

!1) set up unit cell: acell, rprim and rprimd ---------------------
 acell(1:3)=one
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'acell',tacell,'LEN')
 if(tacell==1) acell(1:3)=dprarr(1:3)
 call intagm_img(acell,iimage,jdtset,lenstr,nimage,3,string,"acell",tacell,'LEN')

 scalecart(1:3)=one
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'scalecart',tscalecart,'LEN')
 if(tscalecart==1) scalecart(1:3)=dprarr(1:3)
 call intagm_img(scalecart,iimage,jdtset,lenstr,nimage,3,string,"scalecart",tscalecart,'LEN')

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
 call intagm_img(rprim,iimage,jdtset,lenstr,nimage,3,3,string,"rprim",trprim,'DPR')

!If none of the rprim were read ...
 if(trprim==0)then
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'angdeg',tangdeg,'DPR')
   angdeg(:)=dprarr(1:3)
   call intagm_img(angdeg,iimage,jdtset,lenstr,nimage,3,string,"angdeg",tangdeg,'DPR')

   if(tangdeg==1)then
     call wrtout(std_out,' ingeo: use angdeg to generate rprim.',"COLL")

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
!      write(std_out,*)' ingeo: angdeg=',angdeg(1:3), aa,cc=',aa,cc
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
!  No problem if neither rprim nor angdeg are defined: use default rprim
 end if

!Rescale rprim using scalecart (and set scalecart to one)
 rprim(:,1)=scalecart(:)*rprim(:,1)
 rprim(:,2)=scalecart(:)*rprim(:,2)
 rprim(:,3)=scalecart(:)*rprim(:,3)
 scalecart(:)=one

!Compute the multiplicity of the supercell
 call mati3det(supercell_lattice,multiplicity)
!Get the number of atom in the unit cell
!Read natom from string
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natom',tread,'INT')
!Might initialize natom from XYZ file
 if(tread==0)then
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'_natom',tread,'INT')
 end if
 if(tread==1)natom_uc=intarr(1)

!Store the rprimd of the unit cell
 call mkrdim(acell,rprim,rprimd_read)
!Multiply the rprim to get the rprim of the supercell
 if(multiplicity > 1)then
   rprim(:,1) = rprim(:,1) * supercell_lattice(1,1)
   rprim(:,2) = rprim(:,2) * supercell_lattice(2,2)
   rprim(:,3) = rprim(:,3) * supercell_lattice(3,3)
 end if

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
 if(multiplicity > 1) natrd = natom_uc

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natrd',tnatrd,'INT')
 if(tnatrd==1) natrd=intarr(1)

 if(natrd<1 .or. natrd>natom)then
   if(natrd>1 .and. multiplicity > 1) then
     if(natrd < natom)then
       write(message, '(3a)' )&
&       'The number of atoms to be read (natrd) can not be used with supercell_latt.',ch10,&
&       'Action: Remove natrd or supercell_latt in your input file.'
       MSG_ERROR(message)
     else
       write(message,'(3a,I0,a,I0,a,I0,2a)')&
       'The input variable supercell_latt is present',ch10,&
&      'thus a supercell of ',supercell_lattice(1,1),' ',supercell_lattice(2,2),&
&      ' ',supercell_lattice(3,3),' is generated',ch10
       MSG_WARNING(message)
     end if
   else
     write(message, '(3a,i0,a,i0,2a,a)' )&
&   'The number of atoms to be read (natrd) must be positive and not bigger than natom.',ch10,&
&   'This is not the case: natrd=',natrd,', natom=',natom,ch10,&
&   'Action: correct natrd or natom in your input file.'
     MSG_ERROR(message)
   end if
 end if

!5) Read the type and initial spin of each atom in the primitive set--------
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
   write(message,'(3a)')&
&   'Random positions is a variable defined between 0 and 5. Error in the input file. ',ch10,&
&   'Action: define one of these in your input file.'
   MSG_ERROR(message)
 end if
!if(nimage/=1 .and. iimage/=1)then
!FIXME: should this be called outside the above end if?
 call randomcellpos(natom,npsp,ntypat,random_atpos,ratsph,rprim,rprimd_read,typat_read,&
&                   xred_read(:,1:natrd),znucl,acell)
!This should not be printed if randomcellpos did nothing - it contains garbage. Spurious output anyway
!end if

 call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),'xred',txred,'DPR')
 if(txred==1 .and. txrandom == 0) xred_read(:,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )
 call intagm_img(xred_read,iimage,jdtset,lenstr,nimage,3,natrd,string,"xred",txred,'DPR')

 call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),'xangst',txangst,'DPR')
 if(txangst==1 .and. txrandom==0) xangst_read(:,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )
 call intagm_img(xangst_read,iimage,jdtset,lenstr,nimage,3,natrd,string,"xangst",txangst,'DPR')

 call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),'xcart',txcart,'LEN')
 if(txcart==1 .and. txrandom==0)xcart_read(:,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )
 call intagm_img(xcart_read,iimage,jdtset,lenstr,nimage,3,natrd,string,"xcart",txcart,'LEN')

!Might initialize xred from XYZ file
 if(txred+txcart+txangst+txrandom==0)then
   call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),'_xred',txred,'DPR')
   if(txred==1 .and. txrandom==0) xred_read(:,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )

   call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),'_xangst',txangst,'DPR')
   if(txangst==1 .and. txrandom==0) xangst_read(:,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )
 end if

 if (txred+txcart+txangst+txrandom==0) then
   write(message, '(3a)' )&
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
   call wrtout(std_out,' ingeo: takes atomic coordinates from input array xred ','COLL')
   call xred2xcart(natrd,rprimd_read,xcart_read,xred_read)
 else
   if(txangst==1)then
     call wrtout(std_out,' ingeo: takes atomic coordinates from input array xangst','COLL')
     xcart_read(:,:)=xangst_read(:,:)/Bohr_Ang
   else
     call wrtout(std_out,' ingeo: takes atomic coordinates from input array xcart','COLL')
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
   write(message, '(2(a,i0),5a)')&
&   'Input nsym = ',nsym,' exceeds msym = ',msym,'.',ch10,&
&   'This is not allowed.',ch10,&
&   'Action: correct nsym in your input file.'
   MSG_ERROR(message)
 end if
 if (multiplicity>1) then
   nsym = 1
   MSG_WARNING('Input nsym is now set to one due to the supercell_latt input')
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
&     'Action: either change nsym, or define symrel in your input file.'
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
&         'Action: either change your list of allowed symmetry operations, or use the symmetry finder (nsym=0).'
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
 if(nobj /= 0 .and. multiplicity > 1)then
   write(message, '(3a)' )&
&    'nobj can not be used with supercell_latt.',ch10,&
&    'Action: Remove nobj or supercell_latt in your input file.'
   MSG_ERROR(message)
 end if

!If there are objects, chkprim will not be used immediately
!But, if there are no objects, but a space group, it will be used directly.
 chkprim=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'chkprim',tread,'INT')
 if(tread==1) chkprim=intarr(1)

 if(nobj/=0)then

!  chrgat is read for each atom, from 1 to natom
   call intagm(dprarr,intarr,jdtset,marr,natom,string(1:lenstr),'chrgat',tread,'DPR')
   if(tread==1) then
     chrgat(1:natom) = dprarr(1:natom)
   end if

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

!  Finalize the computation of coordinates: produce xred.
   call xcart2xred(natom,rprimd,xcart,xred)

 else ! nobj==0

!  chrgat is read for each irreducible atom, from 1 to natrd
   call intagm(dprarr,intarr,jdtset,marr,natrd,string(1:lenstr),'chrgat',tread,'DPR')
   if(tread==1)chrgat(1:natrd) = dprarr(1:natrd)

!  Spinat is read for each irreducible atom, from 1 to natrd
   call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),'spinat',tread,'DPR')
   if(tread==1)spinat(1:3,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )

!  nucdipmom is read for each irreducible atom, from 1 to natrd
   call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),'nucdipmom',tread,'DPR')
   if(tread==1)nucdipmom(1:3,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )

!  Compute xred/typat and spinat for the supercell
   if(multiplicity > 1)then
     iatom_supercell = 0
     do i1 = 1, supercell_lattice(1,1)
       do i2 = 1, supercell_lattice(2,2)
         do i3 = 1, supercell_lattice(3,3)
           do iatom = 1, natom_uc
             iatom_supercell = iatom_supercell + 1
             xcart(:,iatom_supercell) = xcart_read(:,iatom) &
&            + matmul(rprimd_read,(/i1-1,i2-1,i3-1/))
             chrgat(iatom_supercell) = chrgat(iatom)
             spinat(1:3,iatom_supercell) = spinat(1:3,iatom)
             typat(iatom_supercell) = typat_read(iatom)
           end do
         end do
       end do
     end do
     call xcart2xred(natom,rprimd,xcart,xred)
   else
!    No supercell
     call xcart2xred(natrd,rprimd,xcart_read,xred)
   end if


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

!      TODO: all the symmetry generation operations should be in one big routine

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

!        Prune the symmetry operations: suppress those with
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

     if(natom/=natrd.and.multiplicity == 1)then
!      Generate the full set of atoms from its knowledge in the irreducible part.
       call fillcell(chrgat,natom,natrd,nsym,nucdipmom,spinat,symafm,symrel,tnons,tolsym,typat,xred)
     end if

!    Check whether the symmetry operations are consistent with the lattice vectors
     iexit=0

     call chkorthsy(gprimd,iexit,nsym,rmet,rprimd,symrel)

   else ! spgroup==0 and nsym==0

!    Here, spgroup==0 as well as nsym==0, so must generate
!    the spatial group of symmetry. However, all the atom
!    positions must be known, so the number
!    of atoms to be read must equal the total number of atoms.
     if(natrd/=natom .and. multiplicity== 1)then
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

       if (multiplicity==1) typat(:)=typat_read(:)
!      Find the symmetry operations: nsym, symafm, symrel and tnons.
!      Use nptsym and ptsymrel, as determined by symlatt
       noncoll=0;if (nspden==4) noncoll=1
       use_inversion=1
       if (dtset%usepaw == 1 .and. (nspden==4.or.pawspnorb>0)) then
         MSG_COMMENT("Removing inversion and improper rotations from initial space group because of PAW + SOC")
         use_inversion=0
       end if

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
&       nzchempot,dtset%prtvol,ptsymrel,spinat,symafm,symrel,tnons,tolsym,typat,use_inversion,xred,&
&       chrgat=chrgat,nucdipmom=nucdipmom)

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

!  Finalize the computation of coordinates: produce xcart
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
 !write(std_out,'(a,3f14.8)') ' ingeo: angdeg(1:3)=',angdeg(1:3)

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

 !call symmultsg(nsym,symafm,symrel,tnons)

!9) initialize the list of fixed atoms, and initial velocities -----------------
!Note: these inputs do not influence the previous generation of
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
 call intagm_img(vel,iimage,jdtset,lenstr,nimage,3,natom,string,"vel",tread,'DPR')

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
&       'Action: check the content of the input variable mixalch.'
       MSG_ERROR(message)
     end if
   end do
   call intagm_img(mixalch,iimage,jdtset,lenstr,nimage,npspalch,ntypalch,string,"mixalch",tread,'DPR')
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
 call intagm_img(amu,iimage,jdtset,lenstr,nimage,ntypat,string,"amu",tread,'DPR')


 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)

end subroutine ingeo
!!***

!!****f* m_ingeo/ingeobld
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

subroutine ingeobld (iout,jdtset,lenstr,natrd,natom,nobj,string,typat,typat_read,xcart,xcart_read)

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
 integer,allocatable :: objaat(:),objbat(:),typat_full(:),vaclst(:)
 real(dp) :: axis2(3),axis3(3),axisa(3),axisb(3),objaax(6),objaro(4),objatr(12)
 real(dp) :: objbax(6),objbro(4),objbtr(12),parall(3),perpen(3),rotated(3)
 real(dp) :: vectora(3),vectorb(3)
 real(dp),allocatable :: xcart_full(:,:)
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

 write(message, '(a,a)' )'--ingeobld: echo values of variables connected to objects --------',ch10
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
 end if !  End condition of existence of a first object

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
 end if !  End condition of existence of a second object

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

   !write(std_out,'(a,i5,a,i2,i2,a)' )' ingeobld : treating iatrd=',iatrd,', belong(a,b)=',belonga,belongb,'.'

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
               angle=objaro(1)+(irep1-1)*objaro(2) + &
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
       end do ! End the repetition loop

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
               angle=objbro(1)+(irep1-1)*objbro(2) + &
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
       end do ! End the repetition loop
     end if ! Condition of belonging to object b
   end if ! Condition of belonging to an object
 end do ! Loop on atoms

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

!Transfer the content of xcart_full and typat_full to the proper location
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

!!****f* m_ingeo/fillcell
!! NAME
!! fillcell
!!
!! FUNCTION
!! Computes the atomic position of all the atoms in the unit cell starting
!! with the symmetry operations and the atoms from the asymetric unit cell.
!!
!! INPUTS
!!  chrgat(natom)=target charge for each atom. Not always used, it depends on the value of constraint_kind
!!  natrd = number of atoms in the assymetric unit cell
!!  natom = total number of atoms (to be checked)
!!  nsym = number of symmetry operations
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrel(3,3,nsym)=symmetry operations in real space in terms
!!   of primitive translations
!!  tnons(3,nsym)=nonsymmorphic translations for symmetry operations
!!  tolsym=tolerance on symmetries
!!  typat(1:natrd)=type integer for each atom in cell
!!  xred(3,1:natrd)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  At input, for the assymetric unit cell
!!  nucdipmom(3,1:natrd)=nuclear magnetic dipole moments of the atoms
!!  spinat(3,1:natrd)=spin-magnetization of the atoms
!!  typat(1:natrd)=type integer for each atom in cell
!!  xred(3,1:natrd)=reduced dimensionless atomic coordinates
!!
!!  At output, for the complete unit cell
!!  nucdipmom(3,1:natom)=nuclear magnetic dipole moments of the atoms
!!  spinat(3,1:natom)=spin-magnetization of the atoms
!!  typat(1:natom)=type integer for each atom in cell
!!  xred(3,1:natom)=reduced dimensionless atomic coordinates
!!
!! PARENTS
!!      ingeo
!!
!! CHILDREN
!!
!! SOURCE

subroutine fillcell(chrgat,natom,natrd,nsym,nucdipmom,spinat,symafm,symrel,tnons,tolsym,typat,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,natrd,nsym
!arrays
 integer,intent(in) :: symafm(nsym),symrel(3,3,nsym)
 integer,intent(inout) :: typat(natom)
 real(dp),intent(in) :: tolsym
 real(dp),intent(in) :: tnons(3,nsym)
 real(dp),intent(inout) :: chrgat(natom),nucdipmom(3,natom),spinat(3,natom),xred(3,natom)

!Local variables ------------------------------
!scalars
 integer :: curat,flagch,flageq,ii,iij,jj,kk
 character(len=500) :: message
!arrays
 integer :: bcktypat(nsym*natrd)
 real(dp) :: bckat(3),bcknucdipmom(3,nsym*natrd)
 real(dp) :: bckchrgat(nsym*natrd),bckspinat(3,nsym*natrd),bckxred(3,nsym*natrd)

! *************************************************************************

!DEBUG
!write(std_out,*)' fillcell : enter with nsym, natrd= ',nsym,natrd
!write(std_out,*)' Describe the different symmetry operations (index,symrel,tnons,symafm)'
!do ii=1,nsym
!write(std_out,'(i3,2x,9i3,3es12.2,i3)')ii,symrel(:,:,ii),tnons(:,ii),symafm(ii)
!end do
!write(std_out,*)' Describe the input atoms (index,typat,xred,spinat)'
!do jj=1,natrd
!write(std_out,'(i3,2x,i3,6es12.2)')jj,typat(jj),xred(:,jj),spinat(:,jj)
!end do
!ENDDEBUG

 curat=0

!Cycle over all the symmetry operations
 do ii=1,nsym

!  Cycle over all the atoms in the assymetric unit cell
   do jj=1,natrd

!    Symmetry operation application
     bckat(:)=matmul(symrel(:,:,ii),xred(:,jj))+tnons(:,ii)

!    Normalization of the coordinates in [0,1)
     do iij=1,3
       do while (bckat(iij)<-tolsym)
         bckat(iij)=bckat(iij)+1.0d0
       end do
       do while (bckat(iij)>=1.0d0-tolsym)
         bckat(iij)=bckat(iij)-1.0d0
       end do
     end do

!    Check for duplicate atoms
     flagch=0
     do kk=1,curat
       flageq=0
       if ( abs(bckxred(1,kk)-bckat(1))<tolsym  .and. &
&       abs(bckxred(2,kk)-bckat(2))<tolsym  .and. &
&       abs(bckxred(3,kk)-bckat(3))<tolsym       ) exit
       flagch=flagch+1
     end do

     if (flagch==curat) then
!      Add the obtained atom to the bckxred list
       curat=curat+1
       bckxred(:,curat)=bckat
       bcktypat(curat)=typat(jj)
       bckchrgat(curat)=chrgat(jj)
       bcknucdipmom(:,curat)=nucdipmom(:,jj)
       bckspinat(:,curat)=spinat(:,jj)*symafm(ii)
     end if

   end do
 end do

!DEBUG
!write(std_out,*)' fillcell : Proposed coordinates ='
!do ii=1,curat
!write(std_out,'(i4,3es16.6)' )ii,bckxred(:,ii)
!end do
!ENDDEBUG

 if (curat>natom) then
   write(message, '(a,i3,a,a,i7,a,a,a,a)' )&
&   'The number of atoms obtained from symmetries, ',curat,ch10,&
&   'is greater than the input number of atoms, natom=',natom,ch10,&
&   'This is not allowed.',ch10,&
&   'Action: modify natom or the symmetry data in the input file.'
   MSG_ERROR(message)
 end if

 if (curat<natom) then
   write(message, '(a,i3,a,a,i7,a,a,a,a)' )&
&   'The number of atoms obtained from symmetries, ',curat,ch10,&
&   'is lower than the input number of atoms, natom=',natom,ch10,&
&   'This is not allowed.',ch10,&
&   'Action: modify natom or the symmetry data in the input file.'
   MSG_ERROR(message)
 end if

!Assignment of symmetry to xred
 xred(:,1:natom)=bckxred(:,1:natom)
 typat(1:natom)=bcktypat(1:natom)
 chrgat(1:natom)=bckchrgat(1:natom)
 nucdipmom(1:3,1:natom)=bcknucdipmom(1:3,1:natom)
 spinat(1:3,1:natom)=bckspinat(1:3,1:natom)

!DEBUG
!write(std_out,*)' fillcell : exit with natom=',natom
!write(std_out,*)' Describe the output atoms (index,typat,xred,spinat)'
!do jj=1,natom
!write(std_out,'(i3,2x,i3,6es12.2)')jj,typat(jj),xred(:,jj),spinat(:,jj)
!end do
!ENDDEBUG

end subroutine fillcell
!!***

!!****f* m_ingeo/invacuum
!!
!! NAME
!! invacuum
!!
!! FUNCTION
!! Determine whether there is vacuum along some of the primitive directions in real space.
!!
!! INPUTS
!! jdtset=number of the dataset looked for
!! lenstr=actual length of the string
!! natom=number of atoms
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! string*(*)=character string containing all the input data.
!!  Initialized previously in instrng.
!! xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!! vacuum(3)= for each direction, 0 if no vacuum, 1 if vacuum
!!
!! PARENTS
!!      invars1,invars2
!!
!! CHILDREN
!!      intagm,metric,sort_dp
!!
!! SOURCE

subroutine invacuum(jdtset,lenstr,natom,rprimd,string,vacuum,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: jdtset,lenstr,natom
 character(len=*),intent(in) :: string
!arrays
 integer,intent(out) :: vacuum(3)
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: ia,ii,marr,tread
 real(dp) :: max_diff_xred,ucvol,vacwidth,vacxred
!arrays
 integer,allocatable :: list(:)
 integer,allocatable :: intarr(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: xred_sorted(:)
 real(dp),allocatable :: dprarr(:)

! *************************************************************************

!Compute the maximum size of arrays intarr and dprarr
 marr=3
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

!Get metric quantities
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Read vacwidth, or set the default
 vacwidth=10.0_dp
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vacwidth',tread,'LEN')
 if(tread==1) vacwidth=dprarr(1)

!Read vacuum, or compute it using the atomic coordinates and vacwidth.
 vacuum(1:3)=0
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'vacuum',tread,'INT')

 if(tread==1)then
   vacuum(1:3)=intarr(1:3)
 else
!  For each direction, determine whether a vacuum space exists
   ABI_ALLOCATE(list,(natom))
   ABI_ALLOCATE(xred_sorted,(natom))
   do ii=1,3
!    This is the minimum xred difference needed to have vacwidth
     vacxred=vacwidth*sqrt(sum(gprimd(:,ii)**2))
!    Project the reduced coordinate in the [0.0_dp,1.0_dp[ interval
     xred_sorted(:)=mod(xred(ii,:),1.0_dp)
!    list is dummy
     list(:)=0
!    Sort xred_sorted
     call sort_dp(natom,xred_sorted,list,tol14)
     if(natom==1)then
       max_diff_xred=1.0_dp
     else
!      Compute the difference between each pair of atom in the sorted order
       max_diff_xred=0.0_dp
       do ia=1,natom-1
         max_diff_xred=max(max_diff_xred,xred_sorted(ia+1)-xred_sorted(ia))
       end do
!      Do not forget the image of the first atom in the next cell
       max_diff_xred=max(max_diff_xred,1.0_dp+xred_sorted(1)-xred_sorted(ia))
     end if
     if(vacxred<max_diff_xred+tol10)vacuum(ii)=1
   end do
   ABI_DEALLOCATE(list)
   ABI_DEALLOCATE(xred_sorted)
 end if

!DEBUG
!write(std_out,*)' invacuum : vacuum=',vacuum(1:3)
!ENDDEBUG

 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)

end subroutine invacuum
!!***

end module m_ingeo
!!***
