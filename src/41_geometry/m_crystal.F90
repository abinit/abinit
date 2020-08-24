!!****m* ABINIT/m_crystal
!! NAME
!! m_crystal
!!
!! FUNCTION
!! Module containing the definition of the crystal_t data type and methods used to handle it.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (MG, YP, MJV)
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
 use m_abicore
 use m_atomdata
 use m_xmpi
 use m_nctk
 use iso_c_binding
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_io_tools,       only : file_exists
 use m_numeric_tools,  only : set2unit
 use m_fstrings,       only : int2char10, sjoin, yesno, itoa
 use m_symtk,          only : mati3inv, sg_multable, symatm, print_symmetries
 use m_spgdata,        only : spgdata
 use m_geometry,       only : metric, xred2xcart, xcart2xred, remove_inversion, getspinrot, symredcart
 use m_io_tools,       only : open_file

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
  !integer,pointer ptsymrel(:,:,:)
  !ptsymrel(3,3,nptsym)
  ! nptsym point-symmetry operations of the Bravais lattice in real space in terms of primitive translations.

  integer :: natom
  ! Number of atoms

  integer :: nsym
  ! Number of symmetry operations

  integer :: ntypat
  ! Number of type of atoms

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
  ! (Anti)Ferromagnetic symmetries. +1/-1

  integer,allocatable :: symrec(:,:,:)
  ! symrec(3,3,nsym)
  ! Symmetry operation in reciprocal space (reduced coordinates)

  integer,allocatable :: symrel(:,:,:)
  ! symrel(3,3,nsym)
  ! Symmetry operations in direct space (reduced coordinates).

  real(dp),allocatable :: symrel_cart(:,:,:)
  ! symrel_cart(3,3,nsym)
  ! Symmetry operations in cartesian coordinates (same order as symrel)

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

 contains

   procedure :: ncwrite => crystal_ncwrite
   ! Write the object in netcdf format

   procedure :: ncwrite_path => crystal_ncwrite_path
   ! Dump the object to netcdf file.

   !procedure :: ncread => crystal_ncread
   ! TODO: Should add routine to read crystal from structure without hdr

   procedure :: isymmorphic
   ! True if space group is symmorphic.

   procedure :: idx_spatial_inversion
   ! Return the index of the spatial inversion, 0 if not present.

   procedure :: isalchemical
   ! True if we are using alchemical pseudopotentials.

   procedure :: free => crystal_free
   ! Free memory.

   procedure :: new_without_symmetries => crystal_without_symmetries
   ! Return new object without symmetries (actually nsym = 1 and identity operation)

   procedure :: get_point_group => crystal_point_group
   ! Return the symmetries of the point group of the crystal.

   procedure :: symbol_type
   ! Return the atomic symbol from the itypat index.

   procedure :: symbol_iatom
   ! Return the atomic symbol from the iatom index.

   procedure :: adata_type
   ! Return atomic data from the itypat index.

   !procedure :: compare => crystal_compare
   ! Compare two structures, write warning messages if they differ

   procedure :: print => crystal_print
   ! Print dimensions and basic info stored in the object

   procedure :: symmetrize_cart_vec3 => crystal_symmetrize_cart_vec3
   ! Symmetrize a 3d cartesian vector

   procedure :: symmetrize_cart_tens33 => crystal_symmetrize_cart_tens33
   ! Symmetrize a cartesian 3x3 tensor

 end type crystal_t

 public :: crystal_init            ! Main Creation method.

 public :: symbols_crystal         ! Return an array with the atomic symbol:["Sr","Ru","O1","O2","O3"]
 public :: prt_cif                 ! Print CIF file.
 public :: prtposcar               ! output VASP style POSCAR and FORCES files.
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
!!      m_crystal,m_ddb,m_dfpt_looppert,m_effective_potential
!!      m_effective_potential_file,m_eig2d,m_gwls_hamiltonian,m_hdr,m_outscfcv
!!      m_precpred_1geo,m_respfn_driver,m_tdep_abitypes,m_unittests,m_vtorho
!!      optic
!!
!! CHILDREN
!!      atomdata_from_znucl
!!
!! SOURCE

subroutine crystal_init(amu,Cryst,space_group,natom,npsp,ntypat,nsym,rprimd,typat,xred,&
   zion,znucl,timrev,use_antiferro,remove_inv,title,&
   symrel,tnons,symafm) ! Optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat,nsym,timrev,space_group,npsp
 type(crystal_t),intent(inout) :: Cryst
 logical,intent(in) :: remove_inv,use_antiferro
!arrays
 integer,intent(in) :: typat(natom)
 integer,optional,intent(in) :: symrel(3,3,nsym),symafm(nsym)
 real(dp),intent(in) :: amu(ntypat),xred(3,natom),rprimd(3,3),zion(ntypat),znucl(npsp)
 real(dp),optional,intent(in) :: tnons(3,nsym)
 character(len=*),intent(in) :: title(ntypat)

!Local variables-------------------------------
!scalars
 integer :: iat,indx,itypat,pinv,isym,nsym_noI
 real(dp) :: tolsym8,ucvol
 !character(len=500) :: msg
!arrays
 integer :: symrec(3,3)
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

 Cryst%amu   = amu
 Cryst%typat = typat
 Cryst%xred  = xred
 Cryst%zion  = zion
 Cryst%znucl = znucl

 call xred2xcart(natom,rprimd,Cryst%xcart,Cryst%xred)

 ABI_MALLOC(Cryst%title,(ntypat))
 Cryst%title = title

 ! Generate index table of atoms, in order for them to be used type after type.
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

 Cryst%timrev = timrev

 if (PRESENT(symrel).and.PRESENT(tnons).and.PRESENT(symafm)) then
   if (.not.remove_inv) then
     ! Just a copy
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
     ! Remove inversion, just to be compatible with old GW implementation
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
       MSG_BUG('Solve the problem with inversion before adding ferromagnetic symmetries')
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
   ! Find symmetries symrec,symrel,tnons,symafm
   ! TODO This should be a wrapper around the abinit library whose usage is not so straightforward
   MSG_BUG('NotImplememented: symrel, symrec and tnons should be specied')
 end if

 ! Get symmetries in cartesian coordinates
 ABI_MALLOC(cryst%symrel_cart, (3, 3, cryst%nsym))
 do isym =1,cryst%nsym
   call symredcart(cryst%rprimd, cryst%gprimd, cryst%symrel_cart(:,:,isym), cryst%symrel(:,:,isym))
   ! purify operations in cartesian coordinates.
   where (abs(cryst%symrel_cart(:,:,isym)) < tol14)
     cryst%symrel_cart(:,:,isym) = zero
   end where
 end do

 ! === Obtain a list of rotated atoms ===
 ! $ R^{-1} (xred(:,iat)-\tau) = xred(:,iat_sym) + R_0 $
 ! * indsym(4,  isym,iat) gives iat_sym in the original unit cell.
 ! * indsym(1:3,isym,iat) gives the lattice vector $R_0$.
 !
 ABI_MALLOC(indsym,(4,Cryst%nsym,natom)); indsym = 0
 tolsym8=tol8
 call symatm(indsym,natom,Cryst%nsym,Cryst%symrec,Cryst%tnons,tolsym8,Cryst%typat,Cryst%xred)

 ABI_MALLOC(Cryst%indsym,(4,Cryst%nsym,natom))
 Cryst%indsym=indsym
 ABI_FREE(indsym)

 ! Rotations in spinor space
 ABI_MALLOC(Cryst%spinrot,(4,Cryst%nsym))
 do isym=1,Cryst%nsym
   call getspinrot(Cryst%rprimd,Cryst%spinrot(:,isym),Cryst%symrel(:,:,isym))
 end do

end subroutine crystal_init
!!***

!!****f* m_crystal/crystal_without_symmetries
!! NAME
!!  crystal_without_symmetries
!!
!! FUNCTION
!!  Return new crystal_t object without symmetries (actually nsym = 1 and identity operation)
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

type(crystal_t) function crystal_without_symmetries(self) result(new)

!Arguments ------------------------------------
 class(crystal_t), intent(in) :: self

!Local variables-------------------------------
 integer,parameter :: timrev1 = 1, new_symafm(1) = 1
 real(dp),parameter :: new_tnons(3,1) = zero
! *************************************************************************

 call crystal_init(self%amu, new, 1, self%natom, self%npsp, self%ntypat, 1, self%rprimd, self%typat, &
  self%xred, self%zion, self%znucl, timrev1, .False., .False., self%title, &
  symrel=identity_3d, tnons=new_tnons, symafm=new_symafm)

end function crystal_without_symmetries
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
!!
!! CHILDREN
!!      atomdata_from_znucl
!!
!! SOURCE

subroutine crystal_free(Cryst)

!Arguments ------------------------------------
 class(crystal_t),intent(inout) :: Cryst

! *********************************************************************

!integer
 ABI_SFREE(Cryst%indsym)
 ABI_SFREE(Cryst%symafm)
 ABI_SFREE(Cryst%symrec)
 ABI_SFREE(Cryst%symrel)
 ABI_SFREE(Cryst%symrel_cart)
 ABI_SFREE(Cryst%atindx)
 ABI_SFREE(Cryst%atindx1)
 ABI_SFREE(Cryst%typat)
 ABI_SFREE(Cryst%nattyp)

!real
 ABI_SFREE(Cryst%tnons)
 ABI_SFREE(Cryst%xcart)
 ABI_SFREE(Cryst%xred)
 ABI_SFREE(Cryst%zion)
 ABI_SFREE(Cryst%znucl)
 ABI_SFREE(Cryst%amu)
 ABI_SFREE(Cryst%spinrot)

!character
 ABI_SFREE(Cryst%title)

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
!!  [prtvol]=Verbosity level. If prtvol== -1, only lattice parameters are printed. Defaults to 0
!!  [mode_paral]=Either "COLL" or "PERS"
!!  [header]=String to be printed as header for additional info.
!!
!! OUTPUT
!!  Only printing
!!
!! PARENTS
!!
!! CHILDREN
!!      atomdata_from_znucl
!!
!! SOURCE

subroutine crystal_print(Cryst, header, unit, mode_paral, prtvol)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit,prtvol
 character(len=*),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header
 class(crystal_t),intent(in) :: Cryst

!Local variables-------------------------------
 integer :: my_unt,my_prtvol,nu,iatom, isym, ii, nsym
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
    'R(',nu,')=',Cryst%rprimd(:,nu)+tol10, &
    'G(',nu,')=',Cryst%gprimd(:,nu)+tol10  ! tol10 is used to be consistent with metric.F90
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
 if (my_prtvol == -1) return

 if (my_prtvol > 0) then
   call print_symmetries(Cryst%nsym,Cryst%symrel,Cryst%tnons,Cryst%symafm,unit=my_unt,mode_paral=my_mode)
   if (Cryst%use_antiferro) call wrtout(my_unt,' System has magnetic symmetries ',my_mode)

   ! Print indsym using the same format as in symatm
   nsym = cryst%nsym
   do iatom=1,cryst%natom
     write(msg, '(a,i0,a)' )' symatm: atom number ',iatom,' is reached starting at atom'
     call wrtout(std_out, msg)
     do ii=1,(nsym-1)/24+1
       if (cryst%natom<100) then
         write(msg, '(1x,24i3)' ) (cryst%indsym(4,isym,iatom),isym=1+(ii-1)*24,min(nsym,ii*24))
       else
         write(msg, '(1x,24i6)' ) (cryst%indsym(4,isym,iatom),isym=1+(ii-1)*24,min(nsym,ii*24))
       end if
       call wrtout(std_out, msg)
     end do
   end do

 end if

 call wrtout(my_unt, " Reduced atomic positions [iatom, xred, symbol]:", my_mode)
 do iatom=1,cryst%natom
   write(msg,"(i5,a,2x,3f11.7,2x,a)")iatom,")",cryst%xred(:,iatom),symbol_type(cryst,cryst%typat(iatom))
   call wrtout(my_unt,msg,my_mode)
 end do

end subroutine crystal_print
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/symbols_crystal
!!
!! NAME
!! symbols_crystal
!!
!! FUNCTION
!! Return a array with the symbol of each atoms with indexation e.g.
!! ["Sr","Ru","O1","O2","O3"]
!!
!! INPUTS
!! natom = number of atoms
!! ntypat = number of typat
!! npsp =  number of pseudopotentials
!! znucl = Nuclear charge for each type of pseudopotential
!!
!! OUTPUT
!! symbols = array with the symbol of each atoms
!!
!! PARENTS
!!      m_effective_potential_file,m_fit_polynomial_coeff,m_opt_effpot
!!      m_polynomial_coeff
!!
!! CHILDREN
!!      atomdata_from_znucl
!!
!! SOURCE

subroutine symbols_crystal(natom,ntypat,npsp,symbols,typat,znucl)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat,npsp
!arrays
 real(dp),intent(in):: znucl(npsp)
 integer,intent(in) :: typat(natom)
 character(len=5),intent(out) :: symbols(natom)
 character(len=3) :: powerchar

!Local variables-------------------------------
!scalar
 integer :: ia,ii,itypat,jj
! *************************************************************************

 ! Fill the symbols array
 do ia=1,natom
   symbols(ia) = adjustl(znucl2symbol(znucl(typat(ia))))
 end do
 itypat = 0
 do itypat =1,ntypat
   ii = 0
   do ia=1,natom
     if(typat(ia)==itypat) then
       ii = ii + 1
     end if
   end do
   if(ii>1)then
     jj=1
     do ia=1,natom
       if(typat(ia)==itypat) then
         write(powerchar,'(I0)') jj
         symbols(ia) = trim(symbols(ia))//trim(powerchar)
         jj=jj+1
       end if
     end do
   end if
 end do

end subroutine symbols_crystal
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

!Arguments ------------------------------------
!scalars
 integer :: inv_idx
 class(crystal_t),intent(in) :: Cryst

!Local variables-------------------------------
!scalars
 integer :: isym

! *************************************************************************

 inv_idx=0
 do isym=1,cryst%nsym
   if (all(cryst%symrel(:,:,isym) == inversion_3d)) then
    inv_idx=isym; return
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

!Arguments ------------------------------------
!scalars
 logical :: ans
 class(crystal_t),intent(in) :: Cryst

! *************************************************************************

 ans = ALL(ABS(Cryst%tnons) < tol6)

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

pure logical function isalchemical(Cryst) result(ans)

!Arguments ------------------------------------
 class(crystal_t),intent(in) :: Cryst

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

type(atomdata_t) function adata_type(crystal, itypat) result(atom)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itypat
 class(crystal_t),intent(in) :: crystal

! *************************************************************************

 call atomdata_from_znucl(atom, crystal%znucl(itypat))

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

function symbol_type(crystal, itypat) result(symbol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itypat
 character(len=2) :: symbol
 class(crystal_t),intent(in) :: crystal

!Local variables-------------------------------
!scalars
 type(atomdata_t) :: atom

! *************************************************************************

 atom = crystal%adata_type(itypat)
 symbol = atom%symbol

end function symbol_type
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/symbol_iatom
!! NAME
!!  symbol_iatom
!!
!! FUNCTION
!!  Return the atomic symbol from the iatom index
!!
!! PARENTS
!!
!! SOURCE

function symbol_iatom(crystal, iatom) result(symbol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iatom
 character(len=2) :: symbol
 class(crystal_t),intent(in) :: crystal

! *************************************************************************

 symbol = crystal%symbol_type(crystal%typat(iatom))

end function symbol_iatom
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/crystal_point_group
!! NAME
!!  crystal_point_group
!!
!! FUNCTION
!!  Return the symmetries of the point group of the crystal.
!!
!! INPUTS
!!  [include_timrev]=If True, time-reversal symmetry is included in the point group unless
!!    the system has spatial inversion. Default: False
!!
!! OUTPUT
!!  ptg_nsym=Number of symmetries in the point group
!!  ptg_symrel(3,3,ptg_nsym)=Rotations in real space
!!  ptg_symrec(3,3,ptg_nsym)=Rotations in reciprocal space
!!  has_inversion=True if spatial inversion is present in the point group.
!!
!! PARENTS
!!
!! CHILDREN
!!      atomdata_from_znucl
!!
!! SOURCE

subroutine crystal_point_group(cryst, ptg_nsym, ptg_symrel, ptg_symrec, has_inversion, include_timrev)

!Arguments ------------------------------------
!scalars
 class(crystal_t),intent(in) :: cryst
 integer,intent(out) :: ptg_nsym
 logical,optional,intent(in) :: include_timrev
 logical,intent(out) :: has_inversion
!arrays
 integer,allocatable,intent(out) :: ptg_symrel(:,:,:),ptg_symrec(:,:,:)

!Local variables-------------------------------
!scalars
 logical :: debug
 integer :: isym,search,tmp_nsym,ierr
 logical :: found,my_include_timrev
!arrays
 integer :: work_symrel(3,3,cryst%nsym)
 integer,allocatable :: symafm(:)
 real(dp),allocatable :: tnons(:,:)

! *************************************************************************

 my_include_timrev = .False.; if (present(include_timrev)) my_include_timrev = include_timrev

 tmp_nsym = 1; work_symrel(:,:,1) = cryst%symrel(:,:,1)
 do isym=2,cryst%nsym
   if (cryst%symafm(isym) == -1) cycle
   do search=1,tmp_nsym
     found = all(work_symrel(:,:,search) == cryst%symrel(:,:,isym))
     if (found) exit
   end do
   if (.not. found) then
     tmp_nsym = tmp_nsym + 1
     work_symrel(:,:,tmp_nsym) = cryst%symrel(:,:,isym)
   end if
 end do

 has_inversion = .False.
 do isym=1,tmp_nsym
   if (all(work_symrel(:,:,isym) == inversion_3d) ) then
     has_inversion = .True.; exit
   end if
 end do

 ! Now we know the symmetries of the point group.
 ptg_nsym = tmp_nsym; if (.not. has_inversion .and. my_include_timrev) ptg_nsym = 2 * tmp_nsym
 ABI_MALLOC(ptg_symrel, (3,3,ptg_nsym))
 ABI_MALLOC(ptg_symrec, (3,3,ptg_nsym))

 ptg_symrel(:,:,1:tmp_nsym) = work_symrel(:,:,1:tmp_nsym)
 do isym=1,tmp_nsym
   call mati3inv(ptg_symrel(:,:,isym), ptg_symrec(:,:,isym))
 end do

 if (.not. has_inversion .and. my_include_timrev) then
   ptg_symrel(:,:,tmp_nsym+1:) = -work_symrel(:,:,1:tmp_nsym)
   do isym=tmp_nsym+1,ptg_nsym
     call mati3inv(ptg_symrel(:,:,isym), ptg_symrec(:,:,isym))
   end do
 end if

 debug = .False.
 if (debug) then
   ABI_CALLOC(tnons, (3, ptg_nsym))
   ABI_MALLOC(symafm, (ptg_nsym))
   symafm = 1
   call sg_multable(ptg_nsym, symafm, ptg_symrel, tnons, tol12, ierr)
   ABI_CHECK(ierr == 0, "point group is not a group! See messages above")
   ABI_FREE(tnons)
   ABI_FREE(symafm)
 end if

end subroutine crystal_point_group
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/crystal_ncwrite
!! NAME
!! crystal_ncwrite
!!
!! FUNCTION
!! Output system geometry to a file, using the NETCDF file format and ETSF I/O.
!! Data are taken from the crystal_t object.
!!
!! INPUTS
!!  cryst<crystal_t>=Object defining the unit cell and its symmetries.
!!  ncid=NC file handle.
!!
!! OUTPUT
!!  Only writing
!!
!! NOTES
!!  Alchemy not treated, since crystal should be initialized at the beginning of the run.
!!
!! PARENTS
!!      anaddb,eig2tot,exc_spectra,dfpt_looppert,m_haydock,m_phonons,m_shirley
!!      outscfcv,sigma
!!
!! CHILDREN
!!      atomdata_from_znucl
!!
!! SOURCE

integer function crystal_ncwrite(cryst, ncid) result(ncerr)

!Arguments ------------------------------------
!scalars
 class(crystal_t),intent(in) :: cryst
 integer,intent(in) :: ncid

#ifdef HAVE_NETCDF
!Local variables-------------------------------
!scalars
 integer :: itypat
 character(len=500) :: msg
 character(len=etsfio_charlen) :: symmorphic
 type(atomdata_t) :: atom
!arrays
 character(len=2) :: symbols(cryst%ntypat)
 character(len=80) :: psp_desc(cryst%ntypat),symbols_long(cryst%ntypat)

! *************************************************************************

 ! TODO alchemy not treated correctly by ETSF_IO specs.
 if (cryst%isalchemical()) then
   write(msg,"(3a)")&
    "Alchemical crystals are not fully supported by the netcdf format",ch10,&
    "Important parameters (e.g. znucl, symbols) are not written with the correct value"
   MSG_WARNING(msg)
 end if

 symmorphic = yesno(cryst%isymmorphic())

 ! Define dimensions.
 ! npsp added in v9.
 ncerr = nctk_def_dims(ncid, [ &
   nctkdim_t("complex", 2), nctkdim_t("symbol_length", 2),&
   nctkdim_t("character_string_length", 80), nctkdim_t("number_of_cartesian_directions", 3),&
   nctkdim_t("number_of_reduced_dimensions", 3), nctkdim_t("number_of_vectors", 3),&
   nctkdim_t("number_of_atoms", cryst%natom), nctkdim_t("number_of_atom_species", cryst%ntypat),&
   nctkdim_t("number_of_atom_pseudopotentials", cryst%npsp),&
   nctkdim_t("number_of_symmetry_operations", cryst%nsym)], defmode=.True.)
 NCF_CHECK(ncerr)

 ! Define variables
 NCF_CHECK(nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "space_group"]))

 ncerr = nctk_def_arrays(ncid, [ &
  ! Atomic structure and symmetry operations
  nctkarr_t("primitive_vectors", "dp", "number_of_cartesian_directions, number_of_vectors"), &
  nctkarr_t("reduced_symmetry_matrices", "int", &
    "number_of_reduced_dimensions, number_of_reduced_dimensions, number_of_symmetry_operations"), &
  nctkarr_t("reduced_symmetry_translations", "dp", "number_of_reduced_dimensions, number_of_symmetry_operations"), &
  nctkarr_t("atom_species", "int", "number_of_atoms"), &
  nctkarr_t("reduced_atom_positions", "dp", "number_of_reduced_dimensions, number_of_atoms"), &
  nctkarr_t("atomic_numbers", "dp", "number_of_atom_species"), &
  nctkarr_t("atom_species_names", "char", "character_string_length, number_of_atom_species"), &
  nctkarr_t("chemical_symbols", "char", "symbol_length, number_of_atom_species"), &
  nctkarr_t('atomic_mass_units', "dp", "number_of_atom_species"), &
  ! Atomic information.
  nctkarr_t("valence_charges", "dp", "number_of_atom_species"), &  ! NB: This variable is not written if alchemical
  nctkarr_t("pseudopotential_types", "char", "character_string_length, number_of_atom_species") &
 ])
 NCF_CHECK(ncerr)

 ! Some variables require the "symmorphic" attribute.
 NCF_CHECK(nf90_put_att(ncid, vid("reduced_symmetry_matrices"), "symmorphic", symmorphic))
 NCF_CHECK(nf90_put_att(ncid, vid("reduced_symmetry_translations"), "symmorphic", symmorphic))

 ! At this point we have an ETSF-compliant file. Add additional data for internal use in abinit.
 ncerr = nctk_def_arrays(ncid, [ &
   nctkarr_t('symafm', "int", "number_of_symmetry_operations"), &
   nctkarr_t('symrel_cart', "dp", "three, three, number_of_symmetry_operations"), &
   nctkarr_t('indsym', "int", "four, number_of_symmetry_operations, number_of_atoms") &
 ])
 NCF_CHECK(ncerr)

 ! Set-up atomic symbols.
 do itypat=1,cryst%ntypat
   call atomdata_from_znucl(atom, cryst%znucl(itypat))
   symbols(itypat) = atom%symbol
   write(symbols_long(itypat),'(a2,a78)') symbols(itypat),REPEAT(CHAR(0),78)
   write(psp_desc(itypat),'(2a)') &
     cryst%title(itypat)(1:MIN(80,LEN_TRIM(cryst%title(itypat)))),REPEAT(CHAR(0),MAX(0,80-LEN_TRIM(cryst%title(itypat))))
 end do

 ! Write data.
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, vid("space_group"), cryst%space_group))
 NCF_CHECK(nf90_put_var(ncid, vid("primitive_vectors"), cryst%rprimd))
 NCF_CHECK(nf90_put_var(ncid, vid("reduced_symmetry_matrices"), cryst%symrel))
 NCF_CHECK(nf90_put_var(ncid, vid("reduced_symmetry_translations"), cryst%tnons))
 NCF_CHECK(nf90_put_var(ncid, vid("atom_species"), cryst%typat))
 NCF_CHECK(nf90_put_var(ncid, vid("reduced_atom_positions"), cryst%xred))
 NCF_CHECK(nf90_put_var(ncid, vid("atomic_numbers"), cryst%znucl(1:cryst%ntypat)))
 NCF_CHECK(nf90_put_var(ncid, vid("atom_species_names"), symbols_long))
 NCF_CHECK(nf90_put_var(ncid, vid("chemical_symbols"), symbols))
 NCF_CHECK(nf90_put_var(ncid, vid('atomic_mass_units'), cryst%amu))
 NCF_CHECK(nf90_put_var(ncid, vid("pseudopotential_types"), psp_desc))
 if (cryst%npsp == cryst%ntypat) then
   NCF_CHECK(nf90_put_var(ncid, vid("valence_charges"), cryst%zion))
 end if

 NCF_CHECK(nf90_put_var(ncid, vid("symafm"), cryst%symafm))
 NCF_CHECK(nf90_put_var(ncid, vid("symrel_cart"), cryst%symrel_cart))
 NCF_CHECK(nf90_put_var(ncid, vid("indsym"), cryst%indsym))

#else
 MSG_ERROR("netcdf library not available")
#endif

contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end function crystal_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/crystal_ncwrite_path
!! NAME
!! crystal_ncwrite_path
!!
!! FUNCTION
!! Output system geometry to a file, using the NETCDF file format and ETSF I/O.
!!
!! INPUTS
!!  crystal<crystal_t>=Object defining the unit cell and its symmetries.
!!  path=filename
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function crystal_ncwrite_path(crystal, path) result(ncerr)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: path
 class(crystal_t),intent(in) :: crystal

#ifdef HAVE_NETCDF
!Local variables-------------------------------
!scalars
 integer :: ncid

! *************************************************************************

 ncerr = nf90_noerr
 if (file_exists(path)) then
   NCF_CHECK(nctk_open_modify(ncid, path, xmpi_comm_self))
 else
   ncerr = nctk_open_create(ncid, path, xmpi_comm_self)
   NCF_CHECK_MSG(ncerr, sjoin("creating:", path))
 end if

 NCF_CHECK(crystal_ncwrite(crystal, ncid))
 NCF_CHECK(nf90_close(ncid))
#endif

end function crystal_ncwrite_path
!!***

!!****f* m_crystal/prt_cif
!! NAME
!! prt_cif
!!
!! FUNCTION
!!   print out CIF format file
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      m_outscfcv
!!
!! CHILDREN
!!      atomdata_from_znucl
!!
!! SOURCE

subroutine prt_cif(brvltt, ciffname, natom, nsym, ntypat, rprimd, &
                   spgaxor, spgroup, spgorig, symrel, tnon, typat, xred, znucl)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom, ntypat, nsym
 integer, intent(in) :: brvltt, spgaxor, spgroup, spgorig
!arrays
 integer, intent(in) :: typat(natom)
 integer, intent(in) :: symrel(3,3,nsym)
 character(len=*), intent(in) :: ciffname
 real(dp), intent(in) :: tnon(3,nsym)
 real(dp), intent(in) :: rprimd(3,3)
 real(dp), intent(in) :: xred(3,natom)
 real(dp), intent(in) :: znucl(ntypat)

!Local variables -------------------------------
!scalars
 integer :: unitcif, iatom, isym, sporder, itypat, nat_this_type
 real(dp) :: ucvol
 type(atomdata_t) :: atom
!arrays
 character(len=80) :: tmpstring
 character(len=1) :: brvsb
 character(len=15) :: intsb,ptintsb,ptschsb,schsb
 character(len=35) :: intsbl
 character(len=10) :: str_nat_type
 character(len=100) :: chemformula
 character(len=500) :: msg
 real(dp) :: angle(3), gprimd(3,3), rmet(3,3), gmet(3,3)

!*************************************************************************

 ! open file in append mode xlf and other compilers refuse append mode
 if (open_file(ciffname,msg,newunit=unitcif) /=0) then
   MSG_WARNING(msg)
   return
 end if

 ! print title for dataset
 write (unitcif,'(a)') 'data_set'

 ! print cell parameters a,b,c, angles, volume
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 angle(1)=acos(rmet(2,3)/sqrt(rmet(2,2)*rmet(3,3)))/two_pi*360.0_dp
 angle(2)=acos(rmet(1,3)/sqrt(rmet(1,1)*rmet(3,3)))/two_pi*360.0_dp
 angle(3)=acos(rmet(1,2)/sqrt(rmet(1,1)*rmet(2,2)))/two_pi*360.0_dp

 write (unitcif,'(a,E20.10)') '_cell_length_a                     ', sqrt(rmet(1,1))*Bohr_Ang
 write (unitcif,'(a,E20.10)') '_cell_length_b                     ', sqrt(rmet(2,2))*Bohr_Ang
 write (unitcif,'(a,E20.10)') '_cell_length_c                     ', sqrt(rmet(3,3))*Bohr_Ang
 write (unitcif,'(a,E20.10)') '_cell_angle_alpha                  ', angle(1)
 write (unitcif,'(a,E20.10)') '_cell_angle_beta                   ', angle(2)
 write (unitcif,'(a,E20.10)') '_cell_angle_gamma                  ', angle(3)
 write (unitcif,'(a,E20.10)') '_cell_volume                       ', ucvol*(Bohr_Ang)**3

 ! print reduced positions
 write (unitcif,'(a)') 'loop_'
 write (unitcif,'(a,E20.10)') '  _atom_site_label                   '
 write (unitcif,'(a,E20.10)') '  _atom_site_fract_x                 '
 write (unitcif,'(a,E20.10)') '  _atom_site_fract_y                 '
 write (unitcif,'(a,E20.10)') '  _atom_site_fract_z                 '
 do iatom = 1, natom
   call atomdata_from_znucl(atom,znucl(typat(iatom)))
   write (unitcif,'(2a,3E20.10)') '  ', atom%symbol, xred(:,iatom)
 end do

!other specs in CIF dictionary which may be useful:
!GEOM_BOND GEOM_ANGLE GEOM_TORSION

 ! print chemical composition in simplest form
 chemformula = "'"
 do itypat = 1, ntypat
   nat_this_type = 0
   do iatom = 1, natom
     if (typat(iatom) == itypat) nat_this_type = nat_this_type+1
   end do
   call atomdata_from_znucl(atom,znucl(itypat))
   call int2char10(nat_this_type, str_nat_type)
   chemformula = trim(chemformula) // atom%symbol // trim(str_nat_type) // "  "
 end do
 chemformula = trim(chemformula) // "'"
 write (unitcif,'(2a)') '_chemical_formula_analytical              ', chemformula

 !FIXME: check that brvltt is correctly used here - is it equal to bravais(1) in the invars routines?
 if (brvltt==1) then
   write (unitcif,'(a)') '_symmetry_cell_setting             triclinic'
 else if(brvltt==2)then
   write (unitcif,'(a)') '_symmetry_cell_setting             monoclinic'
 else if(brvltt==3)then
   write (unitcif,'(a)') '_symmetry_cell_setting             orthorhombic'
 else if(brvltt==4)then
   write (unitcif,'(a)') '_symmetry_cell_setting             tetragonal'
 else if(brvltt==5)then
   write (unitcif,'(a)') '_symmetry_cell_setting             rhombohedral'
 else if(brvltt==6)then
   write (unitcif,'(a)') '_symmetry_cell_setting             hexagonal'
 else if(brvltt==7)then
   write (unitcif,'(a)') '_symmetry_cell_setting             cubic'
 end if

 call spgdata(brvsb,intsb,intsbl,ptintsb,ptschsb,schsb,spgaxor,spgroup,sporder,spgorig)

 ! print symmetry operations
 write (unitcif,'(a,I6)') "_symmetry_Int_Tables_number          ", spgroup
 write (unitcif,'(5a)') "_symmetry_space_group_name_H-M        '", brvsb, " ", trim(intsb), "'"
 write (unitcif,'(a)') ''
 write (unitcif,'(a)') 'loop_'
 write (unitcif,'(a)') '  _symmetry_equiv_pos_as_xyz           '
 do isym = 1, nsym
   call  symrel2string(symrel(:,:,isym), tnon(:,isym), tmpstring)
   write (unitcif,'(2a)') '  ', trim(tmpstring)
 end do

 close(unitcif)

end subroutine prt_cif
!!***

!!****f* m_crystal/symrel2string
!! NAME
!! symrel2string
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      m_crystal
!!
!! CHILDREN
!!      atomdata_from_znucl
!!
!! SOURCE

subroutine symrel2string(symrel1, tnon, string)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: symrel1(3,3)
 real(dp), intent(in) :: tnon(3)
 character(len=80), intent(out) :: string

!Local variables -------------------------
!scalars
 integer :: i1,i2
 character(len=1) :: xyz(3)

! *********************************************************************

 xyz(1) = 'x'
 xyz(2) = 'y'
 xyz(3) = 'z'

 string = ''
 do i1=1,3
   if (abs(tnon(i1)) > tol10) then
     ! find fraction 1/n for tnon, otherwise do not know what to print
     if (abs(one-two*tnon(i1)) < tol10) string = trim(string)//'1/2'
     if (abs(one+two*tnon(i1)) < tol10) string = trim(string)//'-1/2'

     if (abs(one-three*tnon(i1)) < tol10) string = trim(string)//'1/3'
     if (abs(one+three*tnon(i1)) < tol10) string = trim(string)//'-1/3'
     if (abs(two-three*tnon(i1)) < tol10) string = trim(string)//'2/3'
     if (abs(two+three*tnon(i1)) < tol10) string = trim(string)//'-2/3'

     if (abs(one-six*tnon(i1)) < tol10) string = trim(string)//'1/6'
     if (abs(one+six*tnon(i1)) < tol10) string = trim(string)//'-1/6'
     if (abs(five-six*tnon(i1)) < tol10) string = trim(string)//'5/6'
     if (abs(five+six*tnon(i1)) < tol10) string = trim(string)//'-5/6'
   end if
   do i2=1,3
     ! FIXME: check if this is correct ordering for symrel(i1,i2) looks ok
     if (symrel1(i1,i2) == 1)  string = trim(string)//'+'//xyz(i2)
     if (symrel1(i1,i2) == -1) string = trim(string)//'-'//xyz(i2)
   end do
   if (i1 /= 3) string = trim(string)//','
 end do

end subroutine symrel2string
!!***

!!****f* m_crystal/prtposcar
!! NAME
!!  prtposcar
!!
!! FUNCTION
!!  output VASP style POSCAR and FORCES files for use with frozen phonon codes, like
!!  PHON from Dario Alfe' or frophon
!!  IMPORTANT: the order of atoms is fixed such that typat is re-grouped.
!!  First typat=1 then typat=2, etc...
!!  Only master should call this routine in MPI-mode.
!!
!! INPUTS
!!  fcart = forces on atoms in cartesian coordinates
!!  natom = number of atoms
!!  ntypat = number of types of atoms
!!  rprimd = lattice vectors for the primitive cell
!!  typat = type for each of the natom atoms
!!  ucvol = unit cell volume
!!  xred = reduced positions of the atoms
!!  znucl = nuclear charge of each atomic type
!!
!! OUTPUTS
!!   Only files written
!!
!! PARENTS
!!      m_afterscfloop
!!
!! CHILDREN
!!      atomdata_from_znucl
!!
!! SOURCE

subroutine prtposcar(fcart, fnameradix, natom, ntypat, rprimd, typat, ucvol, xred, znucl)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom, ntypat
 real(dp), intent(in) :: ucvol
!arrays
 integer, intent(in) :: typat(natom)
 real(dp), intent(in) :: fcart(3,natom)
 real(dp), intent(in) :: rprimd(3,3)
 real(dp), intent(in) :: xred(3,natom)
 real(dp), intent(in) :: znucl(ntypat)
 character(len=fnlen), intent(in) :: fnameradix

!Local variables-------------------------------
!scalars
 integer :: iatom, itypat, iout
 type(atomdata_t) :: atom
! arrays
 integer :: natoms_this_type(ntypat)
 character(len=2) :: symbol
 character(len=7) :: natoms_this_type_str
 character(len=100) :: chem_formula, natoms_all_types
 character(len=500) :: msg

!************************************************************************

 ! Output POSCAR file for positions, atom types etc
 if (open_file(trim(fnameradix)//"_POSCAR", msg, newunit=iout) /= 0) then
   MSG_ERROR(msg)
 end if

 natoms_this_type = 0
 do itypat=1,ntypat
   do iatom=1,natom
     if (typat(iatom) == itypat) natoms_this_type(itypat) = natoms_this_type(itypat) + 1
   end do
 end do

 chem_formula = ""
 do itypat=1, ntypat
   call atomdata_from_znucl(atom, znucl(itypat))
   symbol = atom%symbol
   if (natoms_this_type(itypat) < 10) then
     write (natoms_this_type_str, '(I1)') natoms_this_type(itypat)
   else if (natoms_this_type(itypat) < 100) then
     write (natoms_this_type_str, '(I2)') natoms_this_type(itypat)
   else if (natoms_this_type(itypat) < 1000) then
     write (natoms_this_type_str, '(I3)') natoms_this_type(itypat)
   end if
   chem_formula = trim(chem_formula) // symbol // trim(natoms_this_type_str)
 end do

 write (iout,'(2a)') "ABINIT generated POSCAR file. Title string - should be chemical formula... ",trim(chem_formula)

 write (iout,'(E24.14)') -ucvol*Bohr_Ang*Bohr_Ang*Bohr_Ang
 write (iout,'(3E24.14,1x)') Bohr_Ang*rprimd(:,1) ! (angstr? bohr?)
 write (iout,'(3E24.14,1x)') Bohr_Ang*rprimd(:,2)
 write (iout,'(3E24.14,1x)') Bohr_Ang*rprimd(:,3)

 natoms_all_types = "   "
 do itypat=1, ntypat
   write (natoms_this_type_str, '(I7)') natoms_this_type(itypat)
   natoms_all_types = trim(natoms_all_types) // "   " // trim(natoms_this_type_str)
 end do

 write (iout,'(a)') trim(natoms_all_types)
 write (iout,'(a)') "Direct"

 do itypat=1, ntypat
   do iatom=1,natom
     if (typat(iatom) /= itypat) cycle
     write (iout,'(3(E24.14,1x))') xred(:,iatom)
   end do
 end do
 close (iout)

 ! output FORCES file for forces in same order as positions above
 if (open_file(trim(fnameradix)//"_FORCES", msg, newunit=iout) /= 0 ) then
   MSG_ERROR(msg)
 end if

 !ndisplacements
 !iatom_displaced displacement_red_coord(3)
 !forces_cart_ev_Angstr(3)
 !...
 !<repeat for other displaced atoms>
 write (iout,'(I7)') 1
 write (iout,'(a)') '1 0 0 0        ! TO BE FILLED IN '
 do itypat=1, ntypat
   do iatom=1,natom
     if (typat(iatom) /= itypat) cycle
     write (iout,'(3(E24.14,1x))') Ha_eV/Bohr_Ang*fcart(:,iatom)
   end do
 end do

 close (iout)

end subroutine prtposcar
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/crystal_symmetrize_cart_vec3
!! NAME
!!  crystal_symmetrize_cart_vec3
!!
!! FUNCTION
!!  Symmetrize a cartesian vector.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function crystal_symmetrize_cart_vec3(cryst, v) result(vsum)

!Arguments ------------------------------------
 class(crystal_t),intent(in) :: cryst
 real(dp),intent(in) :: v(3)
 real(dp) :: vsum(3)

!Local variables-------------------------------
 integer :: isym
 real(dp) :: vsym(3)

 !symmetrize
 vsum = zero
 do isym=1, cryst%nsym
   vsym = matmul( (cryst%symrel_cart(:,:,isym)), v)
   vsum = vsum + vsym
 end do
 vsum = vsum / cryst%nsym

end function crystal_symmetrize_cart_vec3
!!***

!----------------------------------------------------------------------

!!****f* m_crystal/crystal_symmetrize_cart_tens33
!! NAME
!!  crystal_symmetrize_cart_tens33
!!
!! FUNCTION
!!  Symmetrize a cartesian 3x3 tensor
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function crystal_symmetrize_cart_tens33(cryst, t) result(tsum)

!Arguments ------------------------------------
 class(crystal_t),intent(in) :: cryst
 real(dp),intent(in) :: t(3,3)
 real(dp) :: tsum(3,3)

!Local variables-------------------------------
 integer :: isym
 real(dp) :: tsym(3,3)

 !symmetrize
 tsum = zero
 do isym=1, cryst%nsym
   tsym = matmul( (cryst%symrel_cart(:,:,isym)), matmul(t, transpose(cryst%symrel_cart(:,:,isym))) )
   tsum = tsum + tsym
 end do
 tsum = tsum / cryst%nsym

end function crystal_symmetrize_cart_tens33
!!***

END MODULE m_crystal
!!***
