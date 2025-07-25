
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_latt

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi

 use m_matrix,           only : matr3inv, mat33det
 use m_geometry,         only : metric
 use m_tdep_readwrite,   only : Input_type, MPI_enreg_type

 implicit none

  type Lattice_type

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
    double precision :: rprimd_md     (3,3)
    double precision :: Sij           (6,6)
    double precision :: ucvol
    double precision :: BulkModulus_T
    double precision :: BulkModulus_S
    double precision :: HeatCapa_V
    double precision :: HeatCapa_P
    double precision :: Shear
    double precision :: Density

  end type Lattice_type

  public :: tdep_make_inbox
  public :: tdep_make_latt
  public :: tdep_shift_xred

contains

!=====================================================================================================
subroutine tdep_make_inbox(tab,natom,tol,&
&                          temp) !Optional

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
 subroutine tdep_make_latt(Invar,Lattice)

  integer :: brav,ii,jj,line
  double precision :: acell_unitcell(3),multiplicity(3,3),multiplicitym1(3,3)
  double precision :: rprimd(3,3),rprimdt(3,3),rprimd_md(3,3),rprimdm1(3,3)
  double precision :: rprim(3,3),rprimm1(3,3),rprimt(3,3),rprimdtm1(3,3)
  double precision :: rprimd_unitcell(3,3)
  double precision :: rprimd_tmp(3,3), rprimdm1_tmp(3,3)
  double precision :: RAmat(3,3),Amat2(3,3),Rmat(3,3)
  double precision :: rotation(3,3), rotation_cart(3,3)
  double precision :: xi,hh
  character(len=500) :: msg
  type(Input_type),intent(inout) :: Invar
  type(Lattice_type),intent(out) :: Lattice

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
  acell_unitcell(:)=zero; multiplicity(:,:)=zero; multiplicitym1(:,:)=zero
  rprimd(:,:)=zero; rprimdm1(:,:)=zero; rprimdt(:,:)=zero; rprimd_md(:,:)=zero
  rprim(:,:)=zero; rprimm1(:,:)=zero; rprimt(:,:)=zero; rprimdtm1(:,:)=zero
  rprimd_unitcell(:,:)=zero; RAmat(:,:)=zero; Amat2(:,:)=zero; Rmat(:,:)=zero


! ---------------------------------------------------------------------------- !
! Here, rprim defines a Primitive Lattice and NOT a Conventional lattice
! The lattice parameters have to multiply rprim on:
! 0/ line or column         --> line 0
!                               Cubic, Fcc, Bcc, Ortho, Tetra, Rhombo, Hexa
! 1/ line only              --> line=1 (see rprim using acell in ABINIT) :
!                               Mono, Tri
! 2/ column only            --> line=2 (see rprim using scalecart in ABINIT) :
!                               Bct, Face-centered-Ortho, Body-centered-Ortho, C-centered-Ortho
! 3/ neither line or column --> line=3 (rprim has to be directly dimensioned in ABINIT) :
!                               C-centered-Mono
! ---------------------------------------------------------------------------- !
  line=0
! For monoclinic: bravais(1)=2
  if (Invar%bravais(1).eq.2.and.Invar%bravais(2).eq.0) then !monoclinic
    brav=1
    line=1
    rprim(1,1)= 1.0d0                             ; rprim(1,2)= 0.0d0 ; rprim(1,3)= 0.0d0
    rprim(2,1)= 0.0d0                             ; rprim(2,2)= 1.0d0 ; rprim(2,3)= 0.0d0
    rprim(3,1)= dcos(Invar%angle_alpha*pi/180.d0) ; rprim(3,2)= 0.0d0 ; rprim(3,3)= dsin(Invar%angle_alpha*pi/180.d0)
! For orthorhombic: bravais(1)=3
  else if (Invar%bravais(1).eq.3.and.Invar%bravais(2).eq.0) then !orthorhombic
    brav=1
    line=0
    rprim(1,1)=1.0d0 ; rprim(1,2)=0.0d0 ; rprim(1,3)=0.0d0
    rprim(2,1)=0.0d0 ; rprim(2,2)=1.0d0 ; rprim(2,3)=0.0d0
    rprim(3,1)=0.0d0 ; rprim(3,2)=0.0d0 ; rprim(3,3)=1.0d0
  else if (Invar%bravais(1).eq.3.and.Invar%bravais(2).eq.-3) then !face centered orthorhombic
    brav=1
    line=2
    rprim(1,1)=0.0d0 ; rprim(1,2)=0.5d0 ; rprim(1,3)=0.5d0
    rprim(2,1)=0.5d0 ; rprim(2,2)=0.0d0 ; rprim(2,3)=0.5d0
    rprim(3,1)=0.5d0 ; rprim(3,2)=0.5d0 ; rprim(3,3)=0.0d0
  else if (Invar%bravais(1).eq.3.and.Invar%bravais(2).eq.-1) then !body centered orthorhombic
    brav=1
    line=2
    rprim(1,1)=-0.5d0 ; rprim(1,2)= 0.5d0 ; rprim(1,3)= 0.5d0
    rprim(2,1)= 0.5d0 ; rprim(2,2)=-0.5d0 ; rprim(2,3)= 0.5d0
    rprim(3,1)= 0.5d0 ; rprim(3,2)= 0.5d0 ; rprim(3,3)=-0.5d0
  else if (Invar%bravais(1).eq.3.and.Invar%bravais(2).eq.3) then !orthorombic C-face centered
    brav=1
    line=2
    rprim(1,1)= 0.5d0 ; rprim(1,2)=-0.5d0 ; rprim(1,3)=0.0d0
    rprim(2,1)= 0.5d0 ; rprim(2,2)= 0.5d0 ; rprim(2,3)=0.0d0
    rprim(3,1)= 0.0d0 ; rprim(3,2)= 0.0d0 ; rprim(3,3)=1.0d0
! For tetragonal: bravais(1)=4
  else if (Invar%bravais(1).eq.4.and.Invar%bravais(2).eq.0) then !tetragonal
    brav=1
    line=2
    rprim(1,1)=1.0d0 ; rprim(1,2)=0.0d0 ; rprim(1,3)=0.0d0
    rprim(2,1)=0.0d0 ; rprim(2,2)=1.0d0 ; rprim(2,3)=0.0d0
    rprim(3,1)=0.0d0 ; rprim(3,2)=0.0d0 ; rprim(3,3)=1.0d0
  else if (Invar%bravais(1).eq.4.and.Invar%bravais(2).eq.-1) then !body centered tetragonal
    brav=3
    line=2
    rprim(1,1)=-0.5d0 ; rprim(1,2)= 0.5d0 ; rprim(1,3)= 0.5d0
    rprim(2,1)= 0.5d0 ; rprim(2,2)=-0.5d0 ; rprim(2,3)= 0.5d0
    rprim(3,1)= 0.5d0 ; rprim(3,2)= 0.5d0 ; rprim(3,3)=-0.5d0
! For trigonal: bravais(1)=5
  else if (Invar%bravais(1).eq.5.and.Invar%bravais(2).eq.0) then !rhombo
    brav=1
    line=0
    xi=dsin(Invar%angle_alpha*pi/180.d0/2d0)
    hh=dsqrt(1d0-4d0/3d0*xi**2)
    rprim(1,1)=xi    ; rprim(1,2)=-xi/dsqrt(3d0)    ; rprim(1,3)= hh
    rprim(2,1)= 0.d0 ; rprim(2,2)=2d0*xi/dsqrt(3d0) ; rprim(2,3)= hh
    rprim(3,1)=-xi   ; rprim(3,2)=-xi/dsqrt(3d0)    ; rprim(3,3)= hh
!  double precision :: AA,BB,CC,DD
!    AA=dcos(Invar%angle_alpha*pi/180.d0/2d0)
!    BB=dsin(Invar%angle_alpha*pi/180.d0/2d0)
!    CC=dcos(Invar%angle_alpha*pi/180.d0)
!    DD=dsqrt(1-CC**2/AA**2)
!    rprim(1,1)= AA   ; rprim(1,2)= AA  ; rprim(1,3)= CC/AA
!    rprim(2,1)=-BB   ; rprim(2,2)= BB  ; rprim(2,3)= 0.d0
!    rprim(3,1)= 0.d0 ; rprim(3,2)= 0.d0; rprim(3,3)= DD
!    rprim(1,1)= AA    ; rprim(1,2)= -BB ; rprim(1,3)= 0.d0
!    rprim(2,1)= AA    ; rprim(2,2)=  BB ; rprim(2,3)= 0.d0
!    rprim(3,1)= CC/AA ; rprim(3,2)= 0.d0; rprim(3,3)= DD
! For hexagonal: bravais(1)=6
  else if (Invar%bravais(1).eq.6.and.Invar%bravais(2).eq.0) then !hexagonal
    brav=4
    line=0
! The following definition of the hcp is fixed in m_dynmat (chkrp9)
    rprim(1,1)= 1.0d0 ; rprim(1,2)= 0.0d0            ; rprim(1,3)= 0.0d0
    rprim(2,1)=-0.5d0 ; rprim(2,2)= dsqrt(3.d0)/2.d0 ; rprim(2,3)= 0.0d0
    rprim(3,1)= 0.0d0 ; rprim(3,2)= 0.0d0            ; rprim(3,3)= 1.0d0
! For cubic: bravais(1)=7
  else if (Invar%bravais(1).eq.7.and.Invar%bravais(2).eq.0) then !simple cubic
    brav=1
    line=0
    rprim(1,1)=1.0d0 ; rprim(1,2)=0.0d0 ; rprim(1,3)=0.0d0
    rprim(2,1)=0.0d0 ; rprim(2,2)=1.0d0 ; rprim(2,3)=0.0d0
    rprim(3,1)=0.0d0 ; rprim(3,2)=0.0d0 ; rprim(3,3)=1.0d0
  else if (Invar%bravais(1).eq.7.and.Invar%bravais(2).eq.-3) then !face centered cubic
    brav=2
    line=0
    rprim(1,1)=0.0d0 ; rprim(1,2)=0.5d0 ; rprim(1,3)=0.5d0
    rprim(2,1)=0.5d0 ; rprim(2,2)=0.0d0 ; rprim(2,3)=0.5d0
    rprim(3,1)=0.5d0 ; rprim(3,2)=0.5d0 ; rprim(3,3)=0.0d0
  else if (Invar%bravais(1).eq.7.and.Invar%bravais(2).eq.-1) then !body centered cubic
    brav=3
    line=0
    rprim(1,1)=-0.5d0 ; rprim(1,2)= 0.5d0 ; rprim(1,3)= 0.5d0
    rprim(2,1)= 0.5d0 ; rprim(2,2)=-0.5d0 ; rprim(2,3)= 0.5d0
    rprim(3,1)= 0.5d0 ; rprim(3,2)= 0.5d0 ; rprim(3,3)=-0.5d0
  else
    ABI_ERROR('THIS BRAVAIS IS NOT DEFINED')
  end if

! ---------------------------------------------------------------------------- !

! Define inverse of the multiplicity
  multiplicity = Invar%multiplicity
  call matr3inv(multiplicity,multiplicitym1)
  multiplicitym1 = TRANSPOSE(multiplicitym1)

! Compute gprim and (transpose of gprim) gprimt
  call matr3inv(rprim, Lattice%gprimt)
  Lattice%gprim = TRANSPOSE(Lattice%gprimt)

! Define transpose and inverse of rprim
  rprimt = TRANSPOSE(rprim)
  call matr3inv(rprimt, rprimm1)

! Compute acell_unitcell and a rotation matrix, given that
! the supercell rprimd is related to the unitcell rprim according to
!     rprimd  = R * A * multiplicity * rprim
! where R is unitary and A is a diagonal matrix containing acell.
  RAmat = MATMUL(MATMUL(Invar%rprimd_md, rprimm1), multiplicitym1)

! This matrix should be diagonal
  Amat2 = MATMUL(TRANSPOSE(RAmat), RAmat)
  do ii=1,3
    acell_unitcell(ii) = sqrt(Amat2(ii,ii))
  end do
  do ii=1,3
    do jj=1,3
      Rmat(ii,jj) = RAmat(ii,jj) / acell_unitcell(jj)
    end do
  end do
  rotation = TRANSPOSE(Rmat)

! Compute rotation matrix in cartesian coordinates
  call matr3inv(Invar%rprimd_md, rprimdm1_tmp)
  rprimdm1_tmp = TRANSPOSE(rprimdm1_tmp)
  rotation_cart = MATMUL(rotation, Invar%rprimd_md)
  rotation_cart = MATMUL(rprimdm1_tmp, rotation_cart)
  rotation_cart = TRANSPOSE(rotation_cart)

  ! Perform some checks
  if ((mat33det(rotation) - 1) .gt. tol8) then
     rprimd_tmp = MATMUL(multiplicitym1, Invar%rprimd_md)
     write(msg, '(6a,3(3f16.10,1x,a),2a,3(3f16.10,1x,a))')&
      'The input primitive vectors cannot be aligned',ch10,&
      'with the expected primitive vectors through rotation.',ch10,&
      'Input unitcell primitive vectors:',ch10,&
      (rprimd_tmp(1,jj),jj=1,3),ch10,&
      (rprimd_tmp(2,jj),jj=1,3),ch10,&
      (rprimd_tmp(3,jj),jj=1,3),ch10,&
      'Expected primitive vectors:',ch10,&
      (rprim(1,jj),jj=1,3),ch10,&
      (rprim(2,jj),jj=1,3),ch10,&
      (rprim(3,jj),jj=1,3),ch10
    ABI_ERROR(msg)
  end if

! Apply rotation to rprim_md
  Invar%rprimd_md = MATMUL(rotation, Invar%rprimd_md)

! Apply rotation to fcart
  do jj=1,Invar%my_nstep
    do ii=1,Invar%natom
      Invar%fcart(:,ii,jj) = MATMUL(rotation_cart, Invar%fcart(:,ii,jj))
    end do
  end do

! ---------------------------------------------------------------------------- !

! Recompute dimensioned primitive vectors
  rprimd_md = MATMUL(multiplicity, rprim)
  if (line==0.or.line==1) then
    do ii=1,3
      do jj=1,3
        rprimd_md(ii,jj) = acell_unitcell(ii) * rprimd_md(ii,jj)
      end do
    end do
  else if (line==2) then
    do ii=1,3
      do jj=1,3
        rprimd_md(ii,jj) = acell_unitcell(jj) * rprimd_md(ii,jj)
      end do
    end do
  end if

! Echo some (re)computed quantities
  write(Invar%stdout,*) ' '
  write(Invar%stdout,*) '#############################################################################'
  write(Invar%stdout,*) '########################## Computed quantities ##############################'
  write(Invar%stdout,*) '#############################################################################'

! Check the off-diagonal elements
  do ii=1,3
    write(Invar%stdlog,'(a,1x,3(f16.10,1x))') 'The rprimd_md (computed)=',(rprimd_md(ii,jj),jj=1,3)
  end do
  if((abs(rprimd_md(1,2)-Invar%rprimd_md(1,2)).gt.tol5).or.(abs(rprimd_md(1,3)-Invar%rprimd_md(1,3)).gt.tol5).or.&
&    (abs(rprimd_md(2,1)-Invar%rprimd_md(2,1)).gt.tol5).or.(abs(rprimd_md(3,1)-Invar%rprimd_md(3,1)).gt.tol5).or.&
&    (abs(rprimd_md(2,3)-Invar%rprimd_md(2,3)).gt.tol5).or.(abs(rprimd_md(3,2)-Invar%rprimd_md(3,2)).gt.tol5).or.&
&    (abs(rprimd_md(1,1)-Invar%rprimd_md(1,1)).gt.tol5).or.(abs(rprimd_md(2,2)-Invar%rprimd_md(2,2)).gt.tol5).or.&
&    (abs(rprimd_md(3,3)-Invar%rprimd_md(3,3)).gt.tol5)) then
    do ii=1,3
      write(Invar%stdlog,'(a,1x,3(f16.10,1x))') 'The rprimd (from the input file or NetCDF file) is=',&
&                                             (Invar%rprimd_md(ii,jj),jj=1,3)
    end do
    do ii=1,3
      write(Invar%stdlog,'(a,1x,3(f16.10,1x))') 'However, using multiplicity (from the input file)=',&
&                                             (multiplicity(ii,jj),jj=1,3)
    end do
    do ii=1,3
      write(Invar%stdlog,'(a,1x,3(f16.10,1x))') 'rprim (from the aTDEP code)=',(rprim(ii,jj),jj=1,3)
    end do
    write(Invar%stdlog,'(a,1x,3(f16.10,1x))') 'and acell (from the calculation)=',(acell_unitcell(ii),ii=1,3)
    ABI_ERROR(' RPRIMD IS NOT RELATED TO RPRIM AND MULTIPLICITY. MODIFY YOUR RPRIMD.')
  end if

! Check the diagonal elements
  if ((Invar%bravais(1).eq.2).and.(Invar%bravais(2).eq.0)) then !monoclinic
    if ((acell_unitcell(1).gt.acell_unitcell(3)).or.&
&       (acell_unitcell(2).gt.acell_unitcell(3))) then
      ABI_ERROR('You must set a,b <= c in the conventional lattice')
    end if
  else if ((Invar%bravais(1).eq.3).and.(Invar%bravais(2).eq.3)) then !C face centered orthorombique
    if (acell_unitcell(1).ge.acell_unitcell(2)) then
      ABI_ERROR('You must set a < b in the conventional lattice')
    end if
  else if (Invar%bravais(1).eq.6) then !hexagonal
    if(abs(acell_unitcell(1)-acell_unitcell(2)).gt.tol8) then
      ABI_ERROR(' STOP: THE PRECISION ON THE LATTICE PARAMETERS IS NOT SUFFICIENT')
    end if
    acell_unitcell(2)=acell_unitcell(1)
  else if (Invar%bravais(1).eq.7) then !cubic
    if((abs(acell_unitcell(1)-acell_unitcell(2)).gt.tol8).or.&
&      (abs(acell_unitcell(2)-acell_unitcell(3)).gt.tol8)) then
      ABI_ERROR('THE PRECISION ON THE LATTICE PARAMETERS IS NOT SUFFICIENT')
    end if
    acell_unitcell(2)=acell_unitcell(1)
    acell_unitcell(3)=acell_unitcell(1)
  end if
  write(Invar%stdout,'(a,1x,3(f16.10,1x))') ' acell_unitcell=',acell_unitcell(:)

! Recompute rprimd_md with the symmetrized acell,
! in order to have a precision higher than 1.d-8.
  rprimd_md = MATMUL(multiplicity, rprim)
  if (line==0.or.line==1) then
    do ii=1,3
      do jj=1,3
        rprimd_md(ii,jj) = acell_unitcell(ii) * rprimd_md(ii,jj)
      end do
    end do
  else if (line==2) then
    do ii=1,3
      do jj=1,3
        rprimd_md(ii,jj) = acell_unitcell(jj) * rprimd_md(ii,jj)
      end do
    end do
  end if
  do ii=1,3
    write(Invar%stdout,'(a,1x,3(f16.10,1x))') ' rprimd_md=',(rprimd_md(ii,jj),jj=1,3)
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
  end if

! Define transpose and inverse of rprimd
  rprimdt = TRANSPOSE(rprimd)

  ! Compute gmet, rmet, gprimd
  call metric(Lattice%gmet,Lattice%gprimd,Invar%stdlog,Lattice%rmet,rprimdt,Lattice%ucvol)

  call matr3inv(rprimd, rprimdtm1)
  call matr3inv(rprimdt, rprimdm1)

! Store all these values in the 'Lattice' datatype
  Lattice%acell_unitcell(:)  =acell_unitcell(:)
  Lattice%angle_alpha        =Invar%angle_alpha
  Lattice%brav               =brav
  Lattice%bravais(:)         =Invar%bravais(:)
  Lattice%line               =line
  Lattice%multiplicity  (:,:)=multiplicity(:,:)
  Lattice%multiplicitym1(:,:)=multiplicitym1(:,:)
  Lattice%rprim         (:,:)=rprim         (:,:)
  Lattice%rprimt        (:,:)=rprimt        (:,:)
  Lattice%rprimm1       (:,:)=rprimm1       (:,:)
  Lattice%rprimd        (:,:)=rprimd        (:,:)
  Lattice%rprimdm1      (:,:)=rprimdm1      (:,:)
  Lattice%rprimdt       (:,:)=rprimdt       (:,:)
  Lattice%rprimdtm1     (:,:)=rprimdtm1     (:,:)
  Lattice%rprimd_md     (:,:)=rprimd_md     (:,:)

 end subroutine tdep_make_latt

!=====================================================================================================

! Shift xred to keep atoms in the same unit cell at each step.
subroutine tdep_shift_xred(Invar,MPIdata)

  type(Input_type), intent(inout) :: Invar
  type(MPI_enreg_type), intent(in) :: MPIdata
  integer :: natom,ii,iatom,istep,ierr
  integer :: shift,shift_max,shift_best
  double precision :: xi, dist, best_dist
  double precision, allocatable :: x0(:,:)

  natom = Invar%natom
  ABI_MALLOC(x0,(3,natom))

  ! Communicate xred at the first step
  x0(:,:) = zero
  if (MPIdata%my_step(1)) then
    x0(:,:) = Invar%xred(:,:,1)
  end if
  call xmpi_sum(x0,MPIdata%comm_step,ierr)

  ! Shift xred from all steps in the same unitcell as the first step
  shift_max = 1
  do istep=1, Invar%my_nstep
    do iatom=1,natom
      do ii=1,3
        best_dist = abs(Invar%xred(ii,iatom,istep) - x0(ii,iatom))
        shift_best = 0
        do shift=-shift_max,shift_max
          xi = Invar%xred(ii,iatom,istep) + shift
          dist = abs(xi - x0(ii,iatom))
          if (dist < best_dist) then
            best_dist = dist
            shift_best = shift
          end if
        end do
        Invar%xred(ii,iatom,istep) = Invar%xred(ii,iatom,istep) + shift_best
      end do
    end do
  end do

  ABI_FREE(x0)

end subroutine tdep_shift_xred

!=====================================================================================================

end module m_tdep_latt
