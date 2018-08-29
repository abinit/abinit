!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_crystal
!! NAME
!! m_crystal
!!
!! FUNCTION
!! Module containing the definition of the crystal_t data type and methods used to handle it.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (MG, YP, MJV)
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

 use m_numeric_tools,  only : set2unit
 use m_symtk,          only : mati3inv, sg_multable, symatm, print_symmetries
 use m_spgdata,        only : spgdata
 use m_geometry,       only : metric, xred2xcart, remove_inversion, getspinrot
 use m_io_tools,       only : open_file
 use m_fstrings,       only : int2char10

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
  ! (Anti)Ferromagnetic symmetries. +1/-1

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
 public :: idx_spatial_inversion   ! Return the index of the spatial inversion, 0 if not present.
 public :: isymmorphic             ! True if space group is symmorphic.
 public :: isalchemical            ! True if we are using alchemical pseudopotentials.
 public :: adata_type              ! Return atomic data from the itypat index.
 public :: symbol_type             ! Return the atomic symbol from the itypat index.
 public :: symbols_crystal         ! Return an array with the atomic symbol:["Sr","Ru","O1","O2","O3"]
 public :: crystal_point_group     ! Return the symmetries of the point group of the crystal.
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
!!      dfpt_looppert,eig2tot,gwls_hamiltonian,m_crystal_io,m_ddb
!!      m_effective_potential,m_effective_potential_file,m_tdep_abitypes,mover
!!      optic,outscfcv,respfn,vtorho
!!
!! CHILDREN
!!      mati3inv,sg_multable
!!
!! SOURCE

subroutine crystal_init(amu,Cryst,space_group,natom,npsp,ntypat,nsym,rprimd,typat,xred,&
& zion,znucl,timrev,use_antiferro,remove_inv,title,&
& symrel,tnons,symafm) ! Optional


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'crystal_init'
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
 real(dp),intent(in) :: amu(ntypat),xred(3,natom),rprimd(3,3),zion(ntypat),znucl(npsp)
 real(dp),optional,intent(in) :: tnons(3,nsym)
 character(len=*),intent(in) :: title(ntypat)

!Local variables-------------------------------
!scalars
 integer :: iat,indx,itypat,pinv,isym,nsym_noI
 real(dp) :: tolsym8,ucvol
 character(len=500) :: msg
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

 Cryst%timrev = timrev

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
 ABI_MALLOC(indsym,(4,Cryst%nsym,natom)); indsym(:,:,:)=zero
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
!!      anaddb,bethe_salpeter,cut3d,dfpt_looppert,eig2tot,eph,fold2Bloch,gstate
!!      gwls_hamiltonian,m_ddk,m_dvdb,m_effective_potential
!!      m_effective_potential_file,m_gruneisen,m_ioarr,m_iowf,m_wfd,m_wfk
!!      mlwfovlp_qp,mover,mrgscr,optic,outscfcv,respfn,screening,sigma,vtorho
!!      wfk_analyze
!!
!! CHILDREN
!!      mati3inv,sg_multable
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
!!  [prtvol]=Verbosity level. If prtvol== -1, only lattice parameters are printed. Defaults to 0
!!  [mode_paral]=Either "COLL" or "PERS"
!!  [header]=String to be printed as header for additional info.
!!
!! OUTPUT
!!  Only printing
!!
!! PARENTS
!!      eph,gwls_hamiltonian,m_dvdb,m_gruneisen,setup_bse,setup_screening
!!      setup_sigma,wfk_analyze
!!
!! CHILDREN
!!      mati3inv,sg_multable
!!
!! SOURCE

subroutine crystal_print(Cryst,header,unit,mode_paral,prtvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'crystal_print'
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
 if (my_prtvol == -1) return

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

!!****f* m_crystal/symbols_crystal
!!
!! NAME
!! symbols_crystal
!!
!! FUNCTION
!! Return a array with the symbol of each atoms
!! with indexation
!! ["Sr","Ru","O1","O2","O3"] for example
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
!!      m_effective_potential_file,m_fit_polynomial_coeff,m_polynomial_coeff
!!
!! CHILDREN
!!      mati3inv,sg_multable
!!
!! SOURCE

subroutine symbols_crystal(natom,ntypat,npsp,symbols,typat,znucl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symbols_crystal'
!End of the abilint section

 implicit none

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
!arrays
! *************************************************************************

!  Fill the symbols array
   do ia=1,natom
     symbols(ia) = adjustl(znucl2symbol(znucl(typat(ia))))
   end do
   itypat = zero
   do itypat =1,ntypat
     ii = zero
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
!!      m_skw
!!
!! CHILDREN
!!      mati3inv,sg_multable
!!
!! SOURCE

subroutine crystal_point_group(cryst, ptg_nsym, ptg_symrel, ptg_symrec, has_inversion, include_timrev)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'crystal_point_group'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst
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
!!      outscfcv
!!
!! CHILDREN
!!
!! SOURCE

subroutine prt_cif(brvltt, ciffname, natom, nsym, ntypat, rprimd, &
&   spgaxor, spgroup, spgorig, symrel, tnon, typat, xred, znucl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prt_cif'
!End of the abilint section

 implicit none

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
 integer :: unitcif
 integer :: iatom, isym
 integer :: sporder
 integer :: itypat, nat_this_type
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

 real(dp) :: angle(3)
 real(dp) :: gprimd(3,3)
 real(dp) :: rmet(3,3), gmet(3,3)

!*************************************************************************

!open file in append mode xlf and other compilers refuse append mode
 if (open_file(ciffname,msg,newunit=unitcif) /=0) then
   MSG_WARNING(msg)
   return
 end if

!print title for dataset
 write (unitcif,'(a)') 'data_set'

!print cell parameters a,b,c, angles, volume
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

!print reduced positions
 write (unitcif,'(a)') 'loop_'
 write (unitcif,'(a,E20.10)') '  _atom_site_label                   '
 write (unitcif,'(a,E20.10)') '  _atom_site_fract_x                 '
 write (unitcif,'(a,E20.10)') '  _atom_site_fract_y                 '
 write (unitcif,'(a,E20.10)') '  _atom_site_fract_z                 '
 do iatom = 1, natom
   call atomdata_from_znucl(atom,znucl(typat(iatom)))
   write (unitcif,'(2a,3E20.10)') '  ', atom%symbol, xred(:,iatom)
 end do

!
!other specs in CIF dictionary which may be useful:
!GEOM_BOND GEOM_ANGLE GEOM_TORSION
!

!print chemical composition in simplest form
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
 if     (brvltt==1)then
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

!print symmetry operations
 write (unitcif,'(a,I6)') "_symmetry_Int_Tables_number          ", spgroup
 write (unitcif,'(5a)') "_symmetry_space_group_name_H-M        '", brvsb, " ", trim(intsb), "'"
 write (unitcif,'(a)') ''
 write (unitcif,'(a)') 'loop_'
 write (unitcif,'(a)') '  _symmetry_equiv_pos_as_xyz           '
 do isym = 1, nsym
   call  symrel2string(symrel(:,:,isym), tnon(:,isym), tmpstring)
   write (unitcif,'(2a)') '  ', trim(tmpstring)
 end do

 close (unitcif)

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
!!      prt_cif
!!
!! CHILDREN
!!
!! SOURCE

subroutine symrel2string(symrel1, tnon, string)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symrel2string'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
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
!    find fraction 1/n for tnon, otherwise do not know what to print
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
!    FIXME: check if this is correct ordering for symrel(i1,i2) looks ok
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
!!  IMPORTANT: the order of atoms is fixed such that typat is re-grouped. First typat=1
!!  then typat=2, etc...
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
!!      afterscfloop
!!
!! CHILDREN
!!      atomdata_from_znucl
!!
!! SOURCE

subroutine prtposcar(fcart, fnameradix, natom, ntypat, rprimd, typat, ucvol, xred, znucl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtposcar'
!End of the abilint section

 implicit none

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
 character(len=fnlen) :: fname

!************************************************************************

!output POSCAR file for positions, atom types etc
 fname = trim(fnameradix)//"_POSCAR"
 if (open_file(fname,msg,newunit=iout) /= 0) then
   MSG_ERROR(msg)
 end if

 natoms_this_type = 0
 do itypat=1, ntypat
   do iatom=1,natom
     if (typat(iatom) == itypat) natoms_this_type(itypat) = natoms_this_type(itypat)+1
   end do
 end do

 chem_formula = ""
 do itypat=1, ntypat
   call atomdata_from_znucl(atom,znucl(itypat))
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
 write (iout,'(2a)') "ABINIT generated POSCAR file. Title string - should be chemical formula... ",&
& trim(chem_formula)

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


!output FORCES file for forces in same order as positions above
 fname = trim(fnameradix)//"_FORCES"
 if (open_file(fname,msg,newunit=iout) /= 0 ) then
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

END MODULE m_crystal
!!***
