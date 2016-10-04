!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_crystal
!! NAME
!! m_crystal
!!
!! FUNCTION
!! Module containing the definition of the crystal_t data type and methods used to handle it. 
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2016 ABINIT group (MG, YP)
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

MODULE m_crystal

 use defs_basis
 use m_errors 
 use m_profiling_abi
 use m_atomdata

 use m_numeric_tools,  only : set2unit

 implicit none

 private 
!!***

!----------------------------------------------------------------------

!!****t* m_crystal/crystal_t
!! NAME
!! crystal_t
!! 
!! FUNCTION
!! Structure defining the unit cell (geometry, atomic positions and symmetry operations in real and reciprocal space)
!! 
!! SOURCE

 type,public :: crystal_t

!scalars
  !integer :: point_group                    ! Point group
  !integer :: bravais,crystsys               ! Bravais lattice, Crystal system
  !integer :: nptsym                         ! No of point symmetries of the Bravais lattice
  !integer :: bravais(11)                    ! bravais(1)=iholohedry, bravais(2)=center
                                             ! bravais(3:11)=coordinates of rprim in the axes of the conventional
                                             ! bravais lattice (*2 if center/=0)

  !integer :: vacuum(3)
  !integer,pointer ptsymrel(:,:,:)
  !ptsymrel(3,3,nptsym)
  ! nptsym point-symmetry operations of the Bravais lattice in real space in terms of primitive translations.

  integer :: natom     
  ! Number of atoms

  integer :: nsym
  ! Number of symmetry operations

  integer :: ntypat   
  ! Number of type of atoms

  !$integer :: ntypalch,ntyppure

  integer :: npsp  
  ! No. of pseudopotentials

  integer :: space_group 
  ! Space group

  integer :: timrev  
  ! TODO BE CAREFUL here, as the convention used in abinit is different.
  ! 1 => do not use time-reversal symmetry.
  ! 2 => take advantage of time-reversal symmetry. 

  real(dp) :: ucvol                          
  ! Real space unit cell volume.

  logical :: use_antiferro
  ! .TRUE. if AFM symmetries are present and used.

!arrays
  real(dp) :: angdeg(3)
  ! Angles among rprim (degree).

  real(dp) :: gmet(3,3)
  ! Reciprocal space metric ($\textrm{bohr}^{-2}$).

  real(dp) :: gprimd(3,3)
  ! Dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)

  real(dp) :: rmet(3,3)
  ! Metric in real space.

  real(dp) :: rprimd(3,3)
  ! Direct lattice vectors, Bohr units.

  integer,allocatable :: indsym(:,:,:) 
  ! indsym(4,nsym,natom)
  ! indirect indexing array for atoms, see symatm.F90.

  integer,allocatable :: symafm(:)  
  ! symafm(nsym)
  ! (Anti)Ferromagnetic symmetries.

  integer,allocatable :: symrec(:,:,:)  
  ! symrec(3,3,nsym)
  ! Symmetry operation in reciprocal space (reduced coordinates)

  integer,allocatable :: symrel(:,:,:)
  ! symrel(3,3,nsym)
  ! Symmetry operations in direct space (reduced coordinates).

  integer,allocatable :: atindx(:)
  integer,allocatable :: atindx1(:)
  ! atindx(natom), atindx1(natom)
  ! Index tables for atoms useful to treat atoms type after type.

  integer,allocatable :: typat(:) 
  integer,allocatable :: nattyp(:)
  ! typat(natom), nattyp(ntypat)
  ! Type of each natom and number of atoms of each type.

  real(dp),allocatable :: tnons(:,:) 
  ! tnons(3,nsym)
  ! Fractional translations (reduced coordinates)

  real(dp),allocatable :: xcart(:,:) 
  ! xcart(3,natom)
  ! Cartesian coordinates.

  real(dp),allocatable :: xred(:,:)
  ! xred(3,natom)
  ! Reduced coordinates.

  real(dp),allocatable :: spinrot(:,:) 
  ! spinrot(4,nsym)
  ! spinor rotation matrices.

  ! Useful quantities that might be added in the future
  real(dp),allocatable :: amu(:)                
  !  amu(ntypat)
  !  mass of the atoms (atomic mass unit)

  real(dp),allocatable :: zion(:) 
  ! zion(ntypat)
  ! Charge of the pseudo-ion 
  ! (No of valence electrons needed to screen exactly the pseudopotential).

  !real(dp),allocatable :: znucltypat(:)
   ! znucltypat(ntypat)
   ! The atomic number of each type of atom (might be alchemy wrt psps)

  real(dp),allocatable :: znucl(:) 
  ! znucl(npsp)
  ! Nuclear charge for each type of pseudopotential

  character(len=132),allocatable :: title(:) 
   ! title(ntypat)
   ! The content of first line read from the psp file

 end type crystal_t

 public :: crystal_init            ! Main Creation method.
 public :: crystal_free            ! Free memory.
 public :: crystal_print           ! Print dimensions and basic info stored in the object
 public :: print_symmetries        ! Helper function to print symmetries in a nice format.
 public :: idx_spatial_inversion   ! Return the index of the spatial inversion, 0 if not present.
 public :: isymmorphic             ! .TRUE. if space group is symmorphic.
 public :: isalchemical            ! .TRUE. if we are using alchemical pseudopotentials.
 public :: adata_type              ! Return atomic data from the itypat index
 public :: symbol_type             ! Return the atomic symbol from the itypat index
!!***

CONTAINS  !====================================================================================================
!!***

!!****f* m_crystal/crystal_init
!! NAME
!!  crystal_init 
!!
!! FUNCTION
!!  Initialize a crystal_t data type.
!!  Ideally the routine should work in two different modes:
!!  Either the symmetries are directly supplied or the space group
!!  is determined starting from the definition of the unit cell.
!!  Only the first method is implemented, the second one should be
!!  a wrapper for the symmetry finder library. To implement the 
!!  second case I have to add additional entries in the object
!!  and I have also to pass an object describing the (optional) geometry builder.
!!
!! INPUTS
!!  natom=number of atom
!!  ntypat=number of type of atoms
!!  nsym=number of symmetry operations
!!  rprimd(3,3)=dimensional lattive vector (real space)
!!  typat(natom)=type of each atom
!!  xred(3,natom)=reduced coordinates of each atom
!!  symrel(3,3,nsym) [optional]=symmetry operations in real space
!!  space_group=Space group (0 if not available)
!!  tnons(3,nsym) [optional]=fractional Translations
!!  symafm(nsym) [optional]=  ferromagnetic symmetries
!!  remove_inv [optional]= if .TRUE. the inversion is removed from the set of symmetries
!!  timrev ==2 => take advantage of time-reversal symmetry
!!         ==1 ==> do not use time-reversal symmetry 
!!
!! OUTPUT
!!  Cryst<crystal_t>= the object completely initialized.
!!
!! TODO
!!  Add additional entries in the class:
!!  1) Info on space and point group (generators?).
!!  2) alchemy
!!  3) masses and nuclear (pseudo&AE) charge
!!  4) forces stresses, velocities.
!!  5) constraints for the relaxation
!!  6) Likely I will need also info on the electric field and berryopt
!!
!! PARENTS
!!      dfpt_looppert,eig2tot,gwls_hamiltonian,m_crystal_io,m_ddb,mover
!!      outddbnc,outscfcv,vtorho
!!
!! CHILDREN
!!      atomdata_from_znucl,wrtout
!!
!! SOURCE

subroutine crystal_init(Cryst,space_group,natom,npsp,ntypat,nsym,rprimd,typat,xred,&
& zion,znucl,timrev,use_antiferro,remove_inv,title,&
& symrel,tnons,symafm) ! Optional


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'crystal_init'
 use interfaces_32_util
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat,nsym,timrev,space_group,npsp
 type(crystal_t),intent(inout) :: Cryst
 logical,intent(in) :: remove_inv,use_antiferro
!arrays
 integer,intent(in) :: typat(natom)
 integer,optional,intent(in) :: symrel(3,3,nsym),symafm(nsym)
 real(dp),intent(in) :: xred(3,natom),rprimd(3,3),zion(ntypat),znucl(npsp)
 real(dp),optional,intent(in) :: tnons(3,nsym)
 character(len=*),intent(in) :: title(ntypat)

!Local variables-------------------------------
!scalars
 integer :: iat,indx,itypat,pinv,isym,nsym_noI
 real(dp) :: tolsym8,ucvol
 character(len=500) :: msg      
!arrays
 integer :: symrec(3,3),inversion(3,3)
 real(dp) :: gprimd(3,3),gmet(3,3),rmet(3,3)
 integer,pointer :: symrel_noI(:,:,:)
 integer,allocatable :: indsym(:,:,:)
 real(dp),pointer :: tnons_noI(:,:)
! *************************************************************************

 !@crystal_t
 Cryst%natom  = natom 
 Cryst%ntypat = ntypat
 Cryst%npsp   = npsp

 Cryst%space_group = space_group

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 
 Cryst%ucvol  = ucvol
 Cryst%rprimd = rprimd 
 Cryst%rmet   = rmet
 Cryst%gmet   = gmet
 Cryst%gprimd = gprimd

 Cryst%angdeg(1)=ACOS(Cryst%rmet(2,3)/SQRT(Cryst%rmet(2,2)*Cryst%rmet(3,3)))/two_pi*360.0d0
 Cryst%angdeg(2)=ACOS(Cryst%rmet(1,3)/SQRT(Cryst%rmet(1,1)*Cryst%rmet(3,3)))/two_pi*360.0d0
 Cryst%angdeg(3)=ACOS(Cryst%rmet(1,2)/SQRT(Cryst%rmet(1,1)*Cryst%rmet(2,2)))/two_pi*360.0d0

 ABI_MALLOC(Cryst%typat,(natom))
 ABI_MALLOC(Cryst%xred,(3,natom))
 ABI_MALLOC(Cryst%xcart,(3,natom))
 ABI_MALLOC(Cryst%zion,(ntypat))
 ABI_MALLOC(Cryst%znucl,(npsp))
 ABI_MALLOC(Cryst%amu, (ntypat))

 Cryst%typat = typat 
 Cryst%xred  = xred 
 Cryst%zion  = zion
 Cryst%znucl = znucl

 call xred2xcart(natom,rprimd,Cryst%xcart,Cryst%xred)

 ABI_MALLOC(Cryst%title,(ntypat))
 Cryst%title=title
 !
 ! === Generate index table of atoms, in order for them to be used type after type ===
 ABI_MALLOC(Cryst%atindx,(natom))
 ABI_MALLOC(Cryst%atindx1,(natom))
 ABI_MALLOC(Cryst%nattyp,(ntypat))

 indx=1
 do itypat=1,ntypat
   Cryst%nattyp(itypat)=0
   do iat=1,natom
     if (Cryst%typat(iat)==itypat) then
       Cryst%atindx (iat )=indx 
       Cryst%atindx1(indx)=iat
       indx=indx+1
       Cryst%nattyp(itypat)=Cryst%nattyp(itypat)+1
     end if
   end do
 end do

 ! Be careful when we have alchemy.
 !if isalchemical(Cryst)
 !ABI_MALLOC(Cryst%amu, (ntypat))
 !call atomdata_from_znucl(atom,znucl)
 !atom%amu

 Cryst%timrev = timrev
 call set2unit(inversion) ; inversion=-inversion

 if (PRESENT(symrel).and.PRESENT(tnons).and.PRESENT(symafm)) then 
   if (.not.remove_inv) then
     ! * Just a copy 
     Cryst%nsym= nsym
     ABI_MALLOC(Cryst%symrel,(3,3,nsym))
     ABI_MALLOC(Cryst%symrec,(3,3,nsym))
     ABI_MALLOC(Cryst%tnons,(3,nsym))
     ABI_MALLOC(Cryst%symafm,(nsym))
     Cryst%symrel=symrel 
     Cryst%tnons=tnons
     Cryst%symafm=symafm
     Cryst%use_antiferro = use_antiferro
     do isym=1,nsym
       call mati3inv(symrel(:,:,isym),symrec)
       Cryst%symrec(:,:,isym)=symrec
     end do
   else 
     ! * Remove inversion, just to be compatible with old GW implementation
     ! TODO should be removed!
     call remove_inversion(nsym,symrel,tnons,nsym_noI,symrel_noI,tnons_noI,pinv)
     Cryst%nsym=nsym_noI
     ABI_MALLOC(Cryst%symrel,(3,3,nsym_noI))
     ABI_MALLOC(Cryst%symrec,(3,3,nsym_noI))
     ABI_MALLOC(Cryst%tnons,(3,nsym_noI))
     ABI_MALLOC(Cryst%symafm,(nsym_noI))
     Cryst%symrel=symrel_noI
     Cryst%tnons=tnons_noI
     if (ANY(symafm==-1)) then 
       msg = 'Solve the problem with inversion before adding ferromagnetic symmetries '
       MSG_BUG(msg)
     end if
     Cryst%symafm=1
     Cryst%use_antiferro=use_antiferro 
     do isym=1,nsym_noI
       call mati3inv(symrel_noI(:,:,isym),symrec)
       Cryst%symrec(:,:,isym)=symrec
     end do
     ABI_FREE(symrel_noI)
     ABI_FREE(tnons_noI)
   end if

 else
   ! * Find symmetries symrec,symrel,tnons,symafm
   ! TODO This should be a wrapper around the abinit library whose usage is not so straightforward
   MSG_BUG('NotImplememented: symrel, symrec and tnons should be specied')
 end if

 ! === Obtain a list of rotated atoms ===
 ! $ R^{-1} (xred(:,iat)-\tau) = xred(:,iat_sym) + R_0 $ 
 ! * indsym(4,  isym,iat) gives iat_sym in the original unit cell.
 ! * indsym(1:3,isym,iat) gives the lattice vector $R_0$.
 ! 
 ABI_MALLOC(indsym,(4,Cryst%nsym,natom))
 tolsym8=tol8
 call symatm(indsym,natom,Cryst%nsym,Cryst%symrec,Cryst%tnons,tolsym8,Cryst%typat,Cryst%xred)

 ABI_MALLOC(Cryst%indsym,(4,Cryst%nsym,natom))
 Cryst%indsym=indsym  
 ABI_FREE(indsym)

 ! === Rotation in spinor space ===
 ABI_MALLOC(Cryst%spinrot,(4,Cryst%nsym))
 do isym=1,Cryst%nsym
   call getspinrot(Cryst%rprimd,Cryst%spinrot(:,isym),Cryst%symrel(:,:,isym))
 end do
end subroutine crystal_init
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/crystal_free
!! NAME
!!  crystal_free 
!!
!! FUNCTION
!!  Destroy the dynamic arrays in a crystal_t data type.
!!
!! PARENTS
!!      anaddb,bethe_salpeter,dfpt_looppert,eig2tot,eph,gstate,gwls_hamiltonian
!!      m_ddk,m_dvdb,m_ioarr,m_iowf,m_wfd,m_wfk,mlwfovlp_qp,mover,mrgscr
!!      outddbnc,outscfcv,screening,sigma,vtorho,wfk_analyze
!!
!! CHILDREN
!!      atomdata_from_znucl,wrtout
!!
!! SOURCE

subroutine crystal_free(Cryst)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'crystal_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(crystal_t),intent(inout) :: Cryst

! *********************************************************************

 DBG_ENTER("COLL")
 
 !@crystal_t

!integer
 if (allocated(Cryst%indsym))  then
   ABI_FREE(Cryst%indsym)
 end if
 if (allocated(Cryst%symafm))  then
   ABI_FREE(Cryst%symafm)
 end if
 if (allocated(Cryst%symrec))  then
   ABI_FREE(Cryst%symrec)
 end if
 if (allocated(Cryst%symrel))  then
   ABI_FREE(Cryst%symrel)
 end if
 if (allocated(Cryst%atindx))  then
   ABI_FREE(Cryst%atindx)
 end if
 if (allocated(Cryst%atindx1))  then
   ABI_FREE(Cryst%atindx1)
 end if
 if (allocated(Cryst%typat  ))  then
   ABI_FREE(Cryst%typat)
 end if
 if (allocated(Cryst%nattyp))  then
   ABI_FREE(Cryst%nattyp)
 end if

!real
 if (allocated(Cryst%tnons))  then
   ABI_FREE(Cryst%tnons)
 end if
 if (allocated(Cryst%xcart))  then
   ABI_FREE(Cryst%xcart)
 end if
 if (allocated(Cryst%xred))  then
   ABI_FREE(Cryst%xred)
 end if
 if (allocated(Cryst%zion))  then
   ABI_FREE(Cryst%zion)
 end if
 if (allocated(Cryst%znucl))  then
   ABI_FREE(Cryst%znucl)
 end if
 if (allocated(Cryst%amu))  then
   ABI_FREE(Cryst%amu)
 end if
 if (allocated(Cryst%spinrot)) then
    ABI_FREE(Cryst%spinrot)
 end if

!character
 if (allocated(Cryst%title))  then
   ABI_FREE(Cryst%title)
 end if

 DBG_EXIT("COLL")

end subroutine crystal_free
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/crystal_print
!! NAME
!!  crystal_print 
!!
!! FUNCTION
!!  Print the content of crystal_t data type
!!
!! INPUTS
!!  Cryst<crystal_t>=The structure.
!!  [unit]=Unit number for output. Defaults to std_out
!!  [prtvol]=Verbosity level
!!  [mode_paral]=Either "COLL" or "PERS"
!!  [header]=String to be printed as header for additional info.
!!
!! OUTPUT
!!  Only printing 
!!
!! PARENTS
!!      eph,gwls_hamiltonian,setup_bse,setup_screening,setup_sigma,wfk_analyze
!!
!! CHILDREN
!!      atomdata_from_znucl,wrtout
!!
!! SOURCE

subroutine crystal_print(Cryst,header,unit,mode_paral,prtvol) 


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'crystal_print'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit,prtvol
 character(len=*),optional,intent(in) :: mode_paral 
 character(len=*),optional,intent(in) :: header
 type(crystal_t),intent(in) :: Cryst

!Local variables-------------------------------
 integer :: my_unt,my_prtvol,nu,iatom
 character(len=4) :: my_mode
 character(len=500) :: msg      
! ********************************************************************* 

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol 
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=' ==== Info on the Cryst% object ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(a)')' Real(R)+Recip(G) space primitive vectors, cartesian coordinates (Bohr,Bohr^-1):'
 call wrtout(my_unt,msg,my_mode)
 do nu=1,3
   write(msg,'(1x,a,i1,a,3f11.7,2x,a,i1,a,3f11.7)')&
    'R(',nu,')=',Cryst%rprimd(:,nu)+tol10,&
    'G(',nu,')=',Cryst%gprimd(:,nu)+tol10 !tol10 is used to be consistent with metric.F90
   call wrtout(my_unt,msg,my_mode)
 end do

 write(msg,'(a,1p,e15.7,a)')' Unit cell volume ucvol=',Cryst%ucvol+tol10,' bohr^3'
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(a,3es16.8,a)')' Angles (23,13,12)=',Cryst%angdeg(1:3),' degrees'
 call wrtout(my_unt,msg,my_mode)

 if (Cryst%timrev==1) then 
   msg = ' Time-reversal symmetry is not present '
 else if (Cryst%timrev==2) then 
   msg = ' Time-reversal symmetry is present '
 else 
   MSG_BUG('Wrong value for timrev') 
 end if
 call wrtout(my_unt,msg,my_mode)

 call print_symmetries(Cryst%nsym,Cryst%symrel,Cryst%tnons,Cryst%symafm,unit=my_unt,mode_paral=my_mode)

 if (Cryst%use_antiferro) call wrtout(my_unt,' System has magnetic symmetries ',my_mode)

 call wrtout(my_unt,"Reduced atomic positions [iatom, xred, symbol]:",my_mode)
 do iatom=1,cryst%natom
   write(msg,"(i5,a,2x,3f11.7,2x,a)")iatom,")",cryst%xred(:,iatom),symbol_type(cryst,cryst%typat(iatom))
   call wrtout(my_unt,msg,my_mode)
 end do

end subroutine crystal_print
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/print_symmetries
!! NAME
!! print_symmetries
!!
!! FUNCTION
!!  Helper function to print the set of symmetries.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gensymspgr,hdr_vs_dtset,m_crystal
!!
!! CHILDREN
!!      atomdata_from_znucl,wrtout
!!
!! SOURCE

subroutine print_symmetries(nsym,symrel,tnons,symafm,unit,mode_paral)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_symmetries'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,optional,intent(in) :: unit
 character(len=4),optional,intent(in) :: mode_paral
!arrays
 integer,intent(in) :: symrel(3,3,nsym),symafm(nsym)
 real(dp),intent(in) :: tnons(3,nsym)

!Local variables-------------------------------
 integer :: my_unt,isym,isymin,isymend,ii,jj
 character(len=500) :: msg      
 character(len=4) :: my_mode
! *********************************************************************

 my_unt =std_out; if (PRESENT(unit      )) my_unt =unit
 my_mode='COLL' ; if (PRESENT(mode_paral)) my_mode=mode_paral

 !write(msg,'(2a)')ch10,' Rotations                           Translations     Symafm '
 !do isym=1,nsym
 ! write(msg,'(1x,3(3i3,1x),4x,3(f11.7,1x),6x,i2)')symrel(:,:,isym),tnons(:,isym),symafm(isym)
 ! call wrtout(my_unt,msg,my_mode)
 !end do 

 write(msg,'(2a)')ch10,' Symmetry operations in real space (Rotation tnons AFM)'
 call wrtout(my_unt,msg,my_mode)
 do isymin=1,nsym,4
   isymend=isymin+3
   if (isymend>nsym) isymend=nsym
   do ii=1,3
     write(msg,'(4(3i3,f8.3,i3,3x))')((symrel(ii,jj,isym),jj=1,3),tnons(ii,isym),symafm(isym),isym=isymin,isymend)
     call wrtout(my_unt,msg,my_mode)
   end do
   write(msg,'(a)')ch10
   call wrtout(my_unt,msg,my_mode)
 end do

end subroutine print_symmetries 
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/idx_spatial_inversion
!! NAME
!!  idx_spatial_inversion
!!
!! FUNCTION
!!  Return the index of the spatial inversion, 0 if not present
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure function idx_spatial_inversion(Cryst) result(inv_idx)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'idx_spatial_inversion'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: inv_idx
 type(crystal_t),intent(in) :: Cryst

!Local variables-------------------------------
!scalars
 integer :: isym
!arrays
 integer :: inversion(3,3)

! *************************************************************************

 inversion=RESHAPE((/-1,0,0,0,-1,0,0,0,-1/),(/3,3/))

 inv_idx=0
 do isym=1,Cryst%nsym
   if ( ALL(Cryst%symrel(:,:,isym)==inversion) ) then 
    inv_idx=isym; RETURN
   end if
 end do

end function idx_spatial_inversion
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/isymmorphic
!! NAME
!!  isymmorphic
!!
!! FUNCTION
!!  Returns .TRUE. if space group is symmorphic, i.e. all fractional translations are zero.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure function isymmorphic(Cryst) result(ans)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'isymmorphic'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical :: ans
 type(crystal_t),intent(in) :: Cryst

! *************************************************************************

 ans = ALL (ABS(Cryst%tnons)<tol6) 

end function isymmorphic
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/isalchemical
!! NAME
!!  isalchemical
!!
!! FUNCTION
!!  Returns .TRUE. if we are using alchemical pseudopotentials
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure function isalchemical(Cryst) result(ans)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'isalchemical'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical :: ans
 type(crystal_t),intent(in) :: Cryst

! *************************************************************************

 ans = (Cryst%npsp /= Cryst%ntypat)

end function isalchemical
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/adata_type
!! NAME
!!  adata_type
!!
!! FUNCTION
!!  Return atomic data from the itypat index
!!
!! PARENTS
!!
!! SOURCE

function adata_type(crystal,itypat) result(atom)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'adata_type'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itypat
 type(crystal_t),intent(in) :: crystal
 type(atomdata_t) :: atom

! *************************************************************************

 call atomdata_from_znucl(atom,crystal%znucl(itypat)) 

end function adata_type
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/symbol_type
!! NAME
!!  symbol_type
!!
!! FUNCTION
!!  Return the atomic symbol from the itypat index
!!
!! PARENTS
!!
!! SOURCE

function symbol_type(crystal,itypat) result(symbol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symbol_type'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itypat
 character(len=2) :: symbol
 type(crystal_t),intent(in) :: crystal

!Local variables-------------------------------
!scalars
 type(atomdata_t) :: atom

! *************************************************************************

 atom = adata_type(crystal, itypat)
 symbol = atom%symbol

end function symbol_type
!!***

!----------------------------------------------------------------------

END MODULE m_crystal
!!***
