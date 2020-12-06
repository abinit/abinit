
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_latt

 use defs_basis
 use m_abicore
 use m_errors

 use m_symtk,            only : matr3inv
 use m_geometry,         only : metric
 use m_tdep_readwrite,   only : Input_Variables_type

 implicit none

  type Lattice_Variables_type

    double precision :: acell_unitcell(3)
    double precision :: angle_alpha
    integer :: brav
    integer :: bravais(11)
    integer :: line
    double precision :: metmin        (3,3)
    double precision :: minim         (3,3)
    double precision :: multiplicity  (3,3)
    double precision :: multiplicitym1(3,3)
    double precision :: gmet          (3,3)
    double precision :: rmet          (3,3)
    double precision :: gprimd        (3,3)
    double precision :: gprim         (3,3)
    double precision :: gprimt        (3,3)
    double precision :: rprim         (3,3)
    double precision :: rprimt        (3,3)
    double precision :: rprimm1       (3,3)
    double precision :: rprimd        (3,3)
    double precision :: rprimdm1      (3,3)
    double precision :: rprimdt       (3,3)
    double precision :: rprimdtm1     (3,3)
    double precision :: rprimd_MD     (3,3)
    double precision :: Sij           (6,6)
    double precision :: ucvol
    double precision :: BulkModulus_T
    double precision :: BulkModulus_S
    double precision :: HeatCapa_V
    double precision :: HeatCapa_P
    double precision :: Shear
    double precision :: Density

  end type Lattice_Variables_type

  public :: tdep_make_inbox
  public :: tdep_make_latt

contains

!=====================================================================================================
subroutine tdep_make_inbox(tab,natom,tol,&
&                temp) !Optional

  implicit none
  integer :: natom,ii,jj,iatom
  double precision :: tol
  double precision :: tab(3,natom) !Input (and ouput if temp is not present)
  double precision,optional :: temp(3,natom) !Input/Output (if present)

  do iatom=1,natom
    do ii=1,3
      if (tab(ii,iatom).lt.0.d0) then
        jj=int(tab(ii,iatom)-0.5+tol)
      else if (tab(ii,iatom).ge.0.d0) then
        jj=int(tab(ii,iatom)+0.5+tol)
      end if
      if (present(temp)) then
        temp(ii,iatom)=temp(ii,iatom)-real(jj)
      else
        tab(ii,iatom)=tab(ii,iatom)-real(jj)
      end if
    end do
  end do

end subroutine tdep_make_inbox

!=====================================================================================================
 subroutine tdep_make_latt(InVar,Lattice)

  implicit none

  integer :: brav,ii,jj,kk,INFO,line
  integer, allocatable :: IPIV(:)
  double precision :: acell_unitcell(3),multiplicity(3,3),multiplicitym1(3,3),temp2(3,3)
  double precision :: rprimd(3,3),rprimdt(3,3),rprimd_MD(3,3),rprim_tmp(3,3),rprimdm1(3,3)
  double precision :: rprim(3,3),temp(3,3),rprimm1(3,3),rprimt(3,3),rprimdtm1(3,3)
  double precision :: xi,hh
  double precision, allocatable :: WORK(:)
  type(Input_Variables_type) :: InVar
  type(Lattice_Variables_type),intent(out) :: Lattice

! For bravais(1):
! The holohedral groups are numbered as follows
! (see international tables for crystallography (1983), p. 13)
! iholohedry=1   triclinic      1bar
! iholohedry=2   monoclinic     2/m
! iholohedry=3   orthorhombic   mmm
! iholohedry=4   tetragonal     4/mmm
! iholohedry=5   trigonal       3bar m
! iholohedry=6   hexagonal      6/mmm
! iholohedry=7   cubic          m3bar m

! For bravais(2):
! Centering
! center=0        no centering
! center=-1       body-centered
! center=-3       face-centered
! center=1        A-face centered
! center=2        B-face centered
! center=3        C-face centered

! Correspondency between brav and bravais: brav=1-S.C., 2-F.C., 3-B.C., 4-Hex.)

! Initialize some local variables
  acell_unitcell(:)=zero; multiplicity(:,:)=zero; multiplicitym1(:,:)=zero; temp2(:,:)=zero
  rprimd(:,:)=zero; rprimdm1(:,:)=zero; rprimdt(:,:)=zero; rprimd_MD(:,:)=zero; rprim_tmp(:,:)=zero
  rprim(:,:)=zero; temp(:,:)=zero; rprimm1(:,:)=zero; rprimt(:,:)=zero; rprimdtm1(:,:)=zero

! Echo some (re)computed quantities
  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '########################## Computed quantities ##############################'
  write(InVar%stdout,*) '#############################################################################'

! Define inverse of the multiplicity tab
  ABI_MALLOC(IPIV,(3)); IPIV(:)=0
  ABI_MALLOC(WORK,(3)); WORK(:)=0.d0
  multiplicitym1(:,:)=InVar%multiplicity(:,:)
  call DGETRF(3,3,multiplicitym1,3,IPIV,INFO)
  call DGETRI(3,multiplicitym1,3,IPIV,WORK,3,INFO)
  ABI_FREE(IPIV)
  ABI_FREE(WORK)

  do ii=1,3
    do jj=1,3
      do kk=1,3
        temp(ii,jj)=temp(ii,jj)+multiplicitym1(ii,kk)*InVar%rprimd_MD(kk,jj)
      end do
    end do
!FB    write(InVar%stdout,'(a,x,3(f16.10,x))') 'temp=',(temp(ii,jj),jj=1,3)
  end do

!=============================================================================================
! Here, rprim defines a Primitive Lattice and NOT a Conventional lattice
! The lattice parameters have to multiply rprim on:
! 0/ line or column         --> line 0
!                               Cubic, Fcc, Bcc, Ortho, Tetra, Rhombo, Hexa
! 1/ line only              --> line=1 (see rprim using acell in ABINIT) :
!                               Mono, Tri
! 2/ column only            --> line=2 (see rprim using scalecart in ABINIT) :
!                               Bct, Face-centered-Ortho, Base-centered-Ortho, C-centered-Ortho
! 3/ neither line or column --> line=3 (rprim has to be directly dimensioned in ABINIT) :
!                               C-centered-Mono
!=============================================================================================
  line=0
! For monoclinic: bravais(1)=2
  if (InVar%bravais(1).eq.2.and.InVar%bravais(2).eq.0) then !monoclinic
    brav=1
    line=1
!FB See the m_tdep_qpt.F90 routine for the other modifications
!FB    rprim(1,1)= 1.0d0 ; rprim(1,2)= 0.0d0                             ; rprim(1,3)= 0.0d0
!FB    rprim(2,1)= 0.0d0 ; rprim(2,2)= 1.0d0                             ; rprim(2,3)= 0.0d0
!FB    rprim(3,1)= 0.0d0 ; rprim(3,2)= dcos(InVar%angle_alpha*pi/180.d0) ; rprim(3,3)= dsin(InVar%angle_alpha*pi/180.d0)
    rprim(1,1)= 1.0d0                             ; rprim(1,2)= 0.0d0 ; rprim(1,3)= 0.0d0
    rprim(2,1)= 0.0d0                             ; rprim(2,2)= 1.0d0 ; rprim(2,3)= 0.0d0
    rprim(3,1)= dcos(InVar%angle_alpha*pi/180.d0) ; rprim(3,2)= 0.0d0 ; rprim(3,3)= dsin(InVar%angle_alpha*pi/180.d0)
! For orthorhombic: bravais(1)=3
  else if (InVar%bravais(1).eq.3.and.InVar%bravais(2).eq.0) then !orthorhombic
    brav=1
    line=0
    rprim(1,1)=1.0d0 ; rprim(1,2)=0.0d0 ; rprim(1,3)=0.0d0
    rprim(2,1)=0.0d0 ; rprim(2,2)=1.0d0 ; rprim(2,3)=0.0d0
    rprim(3,1)=0.0d0 ; rprim(3,2)=0.0d0 ; rprim(3,3)=1.0d0
  else if (InVar%bravais(1).eq.3.and.InVar%bravais(2).eq.3) then !orthorombic C-face centered
    brav=1
    line=2
    rprim(1,1)= 0.5d0 ; rprim(1,2)=-0.5d0 ; rprim(1,3)=0.0d0
    rprim(2,1)= 0.5d0 ; rprim(2,2)= 0.5d0 ; rprim(2,3)=0.0d0
    rprim(3,1)= 0.0d0 ; rprim(3,2)= 0.0d0 ; rprim(3,3)=1.0d0
! For tetragonal: bravais(1)=4
  else if (InVar%bravais(1).eq.4.and.InVar%bravais(2).eq.0) then !tetragonal
    brav=1
    line=2
    rprim(1,1)=1.0d0 ; rprim(1,2)=0.0d0 ; rprim(1,3)=0.0d0
    rprim(2,1)=0.0d0 ; rprim(2,2)=1.0d0 ; rprim(2,3)=0.0d0
    rprim(3,1)=0.0d0 ; rprim(3,2)=0.0d0 ; rprim(3,3)=1.0d0
  else if (InVar%bravais(1).eq.4.and.InVar%bravais(2).eq.-1) then !body centered tetragonal
    brav=3
    line=2
    rprim(1,1)=-0.5d0 ; rprim(1,2)= 0.5d0 ; rprim(1,3)= 0.5d0
    rprim(2,1)= 0.5d0 ; rprim(2,2)=-0.5d0 ; rprim(2,3)= 0.5d0
    rprim(3,1)= 0.5d0 ; rprim(3,2)= 0.5d0 ; rprim(3,3)=-0.5d0
! For trigonal: bravais(1)=5
  else if (InVar%bravais(1).eq.5.and.InVar%bravais(2).eq.0) then !rhombo
    brav=1
    line=0
    xi=dsin(InVar%angle_alpha*pi/180.d0/2d0)
    hh=dsqrt(1d0-4d0/3d0*xi**2)
    rprim(1,1)=xi    ; rprim(1,2)=-xi/dsqrt(3d0)    ; rprim(1,3)= hh
    rprim(2,1)= 0.d0 ; rprim(2,2)=2d0*xi/dsqrt(3d0) ; rprim(2,3)= hh
    rprim(3,1)=-xi   ; rprim(3,2)=-xi/dsqrt(3d0)    ; rprim(3,3)= hh

! For hexagonal: bravais(1)=6
  else if (InVar%bravais(1).eq.6.and.InVar%bravais(2).eq.0) then !hexagonal
    brav=4
    line=0
    rprim(1,1)= 1.0d0 ; rprim(1,2)= 0.0d0 ; rprim(1,3)= 0.0d0
    rprim(2,1)=-0.5d0 ; rprim(2,2)= dsqrt(3.d0)/2.d0 ; rprim(2,3)= 0.0d0
    rprim(3,1)= 0.0d0 ; rprim(3,2)= 0.0d0 ; rprim(3,3)= 1.0d0
!FB    rprim(1,1)= 1.0d0 ; rprim(1,2)=-0.5d0 ; rprim(1,3)= 0.0d0
!FB    rprim(2,1)= 0.0d0 ; rprim(2,2)= dsqrt(3.d0)/2.d0 ; rprim(2,3)= 0.0d0
!FB    rprim(3,1)= 0.0d0 ; rprim(3,2)= 0.0d0 ; rprim(3,3)= 1.0d0
!FB    rprim(1,1)= dsqrt(3.d0)/2.d0 ; rprim(1,2)= 0.5d0 ; rprim(1,3)= 0.0d0
!FB    rprim(2,1)=-dsqrt(3.d0)/2.d0 ; rprim(2,2)= 0.5d0 ; rprim(2,3)= 0.0d0
!FB    rprim(3,1)= 0.0d0            ; rprim(3,2)= 0.0d0 ; rprim(3,3)= 1.0d0
! For cubic: bravais(1)=7
  else if (InVar%bravais(1).eq.7.and.InVar%bravais(2).eq.0) then !simple cubic
    brav=1
    line=0
    rprim(1,1)=1.0d0 ; rprim(1,2)=0.0d0 ; rprim(1,3)=0.0d0
    rprim(2,1)=0.0d0 ; rprim(2,2)=1.0d0 ; rprim(2,3)=0.0d0
    rprim(3,1)=0.0d0 ; rprim(3,2)=0.0d0 ; rprim(3,3)=1.0d0
  else if (InVar%bravais(1).eq.7.and.InVar%bravais(2).eq.-3) then !face centered cubic
    brav=2
    line=0
    rprim(1,1)=0.0d0 ; rprim(1,2)=0.5d0 ; rprim(1,3)=0.5d0
    rprim(2,1)=0.5d0 ; rprim(2,2)=0.0d0 ; rprim(2,3)=0.5d0
    rprim(3,1)=0.5d0 ; rprim(3,2)=0.5d0 ; rprim(3,3)=0.0d0
  else if (InVar%bravais(1).eq.7.and.InVar%bravais(2).eq.-1) then !body centered cubic
    brav=3
    line=0
    rprim(1,1)=-0.5d0 ; rprim(1,2)= 0.5d0 ; rprim(1,3)= 0.5d0
    rprim(2,1)= 0.5d0 ; rprim(2,2)=-0.5d0 ; rprim(2,3)= 0.5d0
    rprim(3,1)= 0.5d0 ; rprim(3,2)= 0.5d0 ; rprim(3,3)=-0.5d0
  else
    MSG_ERROR('THIS BRAVAIS IS NOT DEFINED')
  end if
! Compute gprim and (transpose of gprim) gprimt
  call matr3inv(rprim,Lattice%gprimt)
  do ii=1,3
    do jj=1,3
      Lattice%gprim(ii,jj)=Lattice%gprimt(jj,ii)
    end do
  end do
! Define transpose and inverse of rprim
  do ii=1,3
    do jj=1,3
      rprimt(ii,jj)=rprim(jj,ii)
    end do
  end do
  ABI_MALLOC(IPIV,(3)); IPIV(:)=0
  ABI_MALLOC(WORK,(3)); WORK(:)=0.d0
  rprimm1(:,:)=rprim(:,:)
  call DGETRF(3,3,rprimm1,3,IPIV,INFO)
  call DGETRI(3,rprimm1,3,IPIV,WORK,3,INFO)
  ABI_FREE(IPIV)
  ABI_FREE(WORK)

! Compute the acell_unitcell according to the definition of
! the unitcell and multiplicity
  do ii=1,3
    do jj=1,3
      do kk=1,3
        temp2(ii,jj)=temp2(ii,jj)+rprimm1(ii,kk)*temp(kk,jj)
      end do
    end do
    acell_unitcell(ii)=temp2(ii,ii)
  end do

! Check the precision and the order of the lattice parameters
  if ((InVar%bravais(1).eq.2).and.(InVar%bravais(2).eq.0)) then !monoclinic
!FB    if ((acell_unitcell(1).gt.acell_unitcell(3)).or.&
!FB&       (acell_unitcell(2).gt.acell_unitcell(3)).or.&
!FB&       (InVar%angle_alpha.ge.90)) then
!FB      write(InVar%stdout,*) ' STOP: You must set a,b <= c and alpha<90 in the conventional lattice'
!FB      stop -1
!FB    end if
    if ((acell_unitcell(1).gt.acell_unitcell(3)).or.&
&       (acell_unitcell(2).gt.acell_unitcell(3))) then
      MSG_ERROR('You must set a,b <= c in the conventional lattice')
    end if
  else if ((InVar%bravais(1).eq.3).and.(InVar%bravais(2).eq.3)) then !C face centered orthorombique
    if (acell_unitcell(1).ge.acell_unitcell(2)) then
      MSG_ERROR('You must set a < b in the conventional lattice')
    end if
  else if (InVar%bravais(1).eq.6) then !hexagonal
    if(abs(acell_unitcell(1)-acell_unitcell(2)).gt.tol8) then
      MSG_ERROR(' STOP: THE PRECISION ON THE LATTICE PARAMETERS IS NOT SUFFICIENT')
    end if
    acell_unitcell(2)=acell_unitcell(1)
  else if (InVar%bravais(1).eq.7) then !cubic
    if((abs(acell_unitcell(1)-acell_unitcell(2)).gt.tol8).or.&
&      (abs(acell_unitcell(2)-acell_unitcell(3)).gt.tol8)) then
      MSG_ERROR('THE PRECISION ON THE LATTICE PARAMETERS IS NOT SUFFICIENT')
    end if
    acell_unitcell(2)=acell_unitcell(1)
    acell_unitcell(3)=acell_unitcell(1)
  end if
  write(InVar%stdout,'(a,1x,3(f16.10,1x))') ' acell_unitcell=',acell_unitcell(:)
! TODO: Check also the off-diagonal elements

! Redefine rprimd_MD (in order to have a precision higher than 1.d-8)
  rprimd_MD(:,:)=0.d0
  rprim_tmp(:,:)=rprim(:,:)
  do ii=1,3
    do jj=1,3
      do kk=1,3
        rprimd_MD(ii,jj)=rprimd_MD(ii,jj)+acell_unitcell(ii)*Invar%multiplicity(ii,kk)*rprim_tmp(kk,jj)
      end do
    end do
    write(InVar%stdout,'(a,1x,3(f16.10,1x))') ' rprimd_MD=',(rprimd_MD(ii,jj),jj=1,3)
  end do

! Define the rprimd with respect to the acell_unitcell, according to the value of "line"
  if (line==2) then
    rprimd(1,1)=rprim(1,1)*acell_unitcell(1) ; rprimd(1,2)=rprim(1,2)*acell_unitcell(2) ; rprimd(1,3)=rprim(1,3)*acell_unitcell(3)
    rprimd(2,1)=rprim(2,1)*acell_unitcell(1) ; rprimd(2,2)=rprim(2,2)*acell_unitcell(2) ; rprimd(2,3)=rprim(2,3)*acell_unitcell(3)
    rprimd(3,1)=rprim(3,1)*acell_unitcell(1) ; rprimd(3,2)=rprim(3,2)*acell_unitcell(2) ; rprimd(3,3)=rprim(3,3)*acell_unitcell(3)
  else if (line==0.or.line==1) then
    rprimd(1,1)=rprim(1,1)*acell_unitcell(1) ; rprimd(1,2)=rprim(1,2)*acell_unitcell(1) ; rprimd(1,3)=rprim(1,3)*acell_unitcell(1)
    rprimd(2,1)=rprim(2,1)*acell_unitcell(2) ; rprimd(2,2)=rprim(2,2)*acell_unitcell(2) ; rprimd(2,3)=rprim(2,3)*acell_unitcell(2)
    rprimd(3,1)=rprim(3,1)*acell_unitcell(3) ; rprimd(3,2)=rprim(3,2)*acell_unitcell(3) ; rprimd(3,3)=rprim(3,3)*acell_unitcell(3)
  else
    write(InVar%stdout,'(a)') ' STOP : CALCULATION OF RPRIMD NOT IMPLEMENTED'
  end if

! Starting from rprimd, compute gmet, rmet, gprimd  
!FB  call metric(Lattice%gmet,Lattice%gprimd,InVar%stdlog,Lattice%rmet,rprimd,Lattice%ucvol)

! Define transpose and inverse of rprimd
  do ii=1,3
    do jj=1,3
      rprimdt(ii,jj)=rprimd(jj,ii)
    end do
  end do  
  call metric(Lattice%gmet,Lattice%gprimd,InVar%stdlog,Lattice%rmet,rprimdt,Lattice%ucvol)
  ABI_MALLOC(IPIV,(3)); IPIV(:)=0
  ABI_MALLOC(WORK,(3)); WORK(:)=0.d0
  rprimdtm1(:,:)=rprimdt(:,:)
  call DGETRF(3,3,rprimdtm1,3,IPIV,INFO)
  call DGETRI(3,rprimdtm1,3,IPIV,WORK,3,INFO)
  ABI_FREE(IPIV)
  ABI_FREE(WORK)

  ABI_MALLOC(IPIV,(3)); IPIV(:)=0
  ABI_MALLOC(WORK,(3)); WORK(:)=0.d0
  rprimdm1(:,:)=rprimd(:,:)
  call DGETRF(3,3,rprimdm1,3,IPIV,INFO)
  call DGETRI(3,rprimdm1,3,IPIV,WORK,3,INFO)
  ABI_FREE(IPIV)
  ABI_FREE(WORK)

! Store all these values in the 'Lattice' datatype
  Lattice%acell_unitcell(:)  =acell_unitcell(:)
  Lattice%angle_alpha        =InVar%angle_alpha
  Lattice%brav               =brav
  Lattice%bravais(:)         =InVar%bravais(:)
  Lattice%line               =line
  Lattice%multiplicity  (:,:)=Invar%multiplicity(:,:)
  Lattice%multiplicitym1(:,:)=multiplicitym1(:,:)
  Lattice%rprim         (:,:)=rprim         (:,:)
  Lattice%rprimt        (:,:)=rprimt        (:,:)
  Lattice%rprimm1       (:,:)=rprimm1       (:,:)
  Lattice%rprimd        (:,:)=rprimd        (:,:)
  Lattice%rprimdm1      (:,:)=rprimdm1      (:,:)
  Lattice%rprimdt       (:,:)=rprimdt       (:,:)
  Lattice%rprimdtm1     (:,:)=rprimdtm1     (:,:)
  Lattice%rprimd_MD     (:,:)=rprimd_MD     (:,:)

 end subroutine tdep_make_latt
!=====================================================================================================

end module m_tdep_latt
