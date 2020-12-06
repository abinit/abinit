
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_qpt

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_tdep_readwrite,   only : Input_Variables_type, MPI_enreg_type
 use m_tdep_latt,        only : Lattice_Variables_type

 implicit none

  type QptBound_type

    integer :: ihol,center
    character (len=5) :: letter
    double precision :: x,y,z

  end type QptBound_type

  type Qpoints_type

    integer :: nqpt,qpt_tot,qptbound_tot
    integer, allocatable :: lgth_segments(:)
    double precision, allocatable :: qpt_red(:,:),qpt_cart(:,:)
    double precision, allocatable :: special_red(:,:),special_cart(:,:)
    character (len=5), allocatable :: special_qpt(:)

  end type Qpoints_type

  public :: tdep_make_qptpath
  public :: tdep_make_specialqpt

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_make_specialqpt(InVar,Lattice,MPIdata,Qpt,QptBound)

  implicit none
  integer :: qpt_tot,qptbound_tot
  double precision :: zeta,eta,nu,angle_alpha
  type(Input_Variables_type),intent(in) :: InVar
  type(QptBound_type), allocatable,intent(out) :: QptBound(:)
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Qpoints_type),intent(out) :: Qpt
  type(MPI_enreg_type), intent(in) :: MPIdata

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

! Define the special Q points IN GENERAL
! Here we use the definitions of special Q points in reduced coordinates
! as defined in the article: Setyawan and Curtarolo CMS 49, 299 (2010)
  if ((InVar%bravais(1).eq.2).and.(InVar%bravais(2).eq.0)) then
!FB    qptbound_tot=16
    qptbound_tot=10
    angle_alpha=Lattice%angle_alpha
    eta=(1.-Lattice%acell_unitcell(2)*dcos(angle_alpha*pi/180.d0)/Lattice%acell_unitcell(3))/(2.*dsin(angle_alpha*pi/180.d0)**2)
    nu = 1./2.-eta*Lattice%acell_unitcell(3)*dcos(angle_alpha*pi/180.d0)/Lattice%acell_unitcell(2)
    ABI_MALLOC(QptBound,(qptbound_tot))
    QptBound(:)=(/ QptBound_type (2, 0,'G ', 0.000, 0.000, 0.000),&
&                  QptBound_type (2, 0,'A ', 0.500, 0.500, 0.000),&
&                  QptBound_type (2, 0,'C ', 0.000, 0.500, 0.500),&
&                  QptBound_type (2, 0,'D ', 0.500, 0.000, 0.500),&
&                  QptBound_type (2, 0,'D1', 0.500, 0.000,-0.500),&
&                  QptBound_type (2, 0,'E ', 0.500, 0.500, 0.500),&
!FB&                  QptBound_type (2, 0,'H ', 0.000, eta  , 1-nu ),&
!FB&                  QptBound_type (2, 0,'H1', 0.000, 1-eta, nu   ),&
!FB&                  QptBound_type (2, 0,'H2', 0.000, eta  , -nu  ),&
!FB&                  QptBound_type (2, 0,'M ', 0.500, eta  , 1-nu ),&
!FB&                  QptBound_type (2, 0,'M1', 0.500, 1-eta, nu   ),&
!FB&                  QptBound_type (2, 0,'M2', 0.500, eta  , -nu  ),&
&                  QptBound_type (2, 0,'X ', 0.000, 0.500, 0.000),&
&                  QptBound_type (2, 0,'Y ', 0.000, 0.000, 0.500),&
&                  QptBound_type (2, 0,'Y1', 0.000, 0.000,-0.500),&
&                  QptBound_type (2, 0,'Z ', 0.500, 0.000, 0.000) /)
  else if ((InVar%bravais(1).eq.3).and.(InVar%bravais(2).eq.0)) then
    qptbound_tot=8
    ABI_MALLOC(QptBound,(qptbound_tot))
    QptBound(:)=(/ QptBound_type (3, 0,'G ', 0.000, 0.000, 0.000),&
&                  QptBound_type (3, 0,'R ', 0.500, 0.500, 0.500),&
&                  QptBound_type (3, 0,'S ', 0.500, 0.500, 0.000),&
&                  QptBound_type (3, 0,'T ', 0.000, 0.500, 0.500),&
&                  QptBound_type (3, 0,'U ', 0.500, 0.000, 0.500),&
&                  QptBound_type (3, 0,'X ', 0.500, 0.000, 0.000),&
&                  QptBound_type (3, 0,'Y ', 0.000, 0.500, 0.000),&
&                  QptBound_type (3, 0,'Z ', 0.000, 0.000, 0.500) /)
  else if ((InVar%bravais(1).eq.3).and.(InVar%bravais(2).eq.3)) then
    zeta=(1.d0+Lattice%acell_unitcell(1)**2/Lattice%acell_unitcell(2)**2)/4.d0
    qptbound_tot=13
    ABI_MALLOC(QptBound,(qptbound_tot))
    QptBound(:)=(/ QptBound_type (3, 3,'G ', 0.000, 0.000, 0.000),&
&                  QptBound_type (3, 3,'Gp', 1.000, 0.000, 0.000),&
&                  QptBound_type (3, 3,'A ', zeta , zeta , 0.500),&
&                  QptBound_type (3, 3,'A1',-zeta ,1-zeta, 0.500),&
&                  QptBound_type (3, 3,'R ', 0.000, 0.500, 0.500),&
&                  QptBound_type (3, 3,'S ', 0.000, 0.500, 0.000),&
&                  QptBound_type (3, 3,'T ',-0.500, 0.500, 0.500),&
&                  QptBound_type (3, 3,'X ', zeta , zeta , 0.000),&
&                  QptBound_type (3, 3,'X1',-zeta ,1-zeta, 0.000),&
&                  QptBound_type (3, 3,'Y ',-0.500, 0.500, 0.000),&
&                  QptBound_type (3, 3,'Yp', 0.500, 0.500, 0.000),&
&                  QptBound_type (3, 3,'Z ', 0.000, 0.000, 0.500),&
&                  QptBound_type (3, 3,'Zp', 1.000, 0.000, 0.500) /)
  else if ((InVar%bravais(1).eq.4).and.(InVar%bravais(2).eq.0)) then
    qptbound_tot=6
    ABI_MALLOC(QptBound,(qptbound_tot))
    QptBound(:)=(/ QptBound_type (4, 0,'G ', 0.000, 0.000, 0.000),&
&                  QptBound_type (4, 0,'A ', 0.500, 0.500, 0.500),&
&                  QptBound_type (4, 0,'M ', 0.500, 0.500, 0.000),&
&                  QptBound_type (4, 0,'R ', 0.000, 0.500, 0.500),&
&                  QptBound_type (4, 0,'X ', 0.000, 0.500, 0.000),&
&                  QptBound_type (4, 0,'Z ', 0.000, 0.000, 0.500) /)
  else if ((InVar%bravais(1).eq.4).and.(InVar%bravais(2).eq.-1)) then
    if (Lattice%acell_unitcell(3).lt.Lattice%acell_unitcell(1)) then
      qptbound_tot=7
      eta=(1.d0+Lattice%acell_unitcell(3)**2/Lattice%acell_unitcell(1)**2)/4.d0
      ABI_MALLOC(QptBound,(qptbound_tot))
      QptBound(:)=(/ QptBound_type (4,-1,'G ', 0.000, 0.000, 0.000),&
&                    QptBound_type (4,-1,'M ',-0.500, 0.500, 0.500),&
&                    QptBound_type (4,-1,'N ', 0.000, 0.500, 0.500),&
&                    QptBound_type (4,-1,'P ', 0.250, 0.250, 0.250),&
&                    QptBound_type (4,-1,'X ', 0.000, 0.000, 0.500),&
&                    QptBound_type (4,-1,'Z ', eta  , eta  ,-eta  ),&
&                    QptBound_type (4,-1,'Z1',-eta  ,1.-eta, eta  ) /)
    else
      qptbound_tot=9
      eta=(1.d0+Lattice%acell_unitcell(1)**2/Lattice%acell_unitcell(3)**2)/4.d0
      zeta=Lattice%acell_unitcell(1)**2/Lattice%acell_unitcell(3)**2/2.d0
      ABI_MALLOC(QptBound,(qptbound_tot))
      QptBound(:)=(/ QptBound_type (4,-1,'G ', 0.000, 0.000, 0.000),&
&                    QptBound_type (4,-1,'N ', 0.000, 0.500, 0.000),&
&                    QptBound_type (4,-1,'P ', 0.250, 0.250, 0.250),&
&                    QptBound_type (4,-1,'S ',-eta  , eta  , eta  ),&
&                    QptBound_type (4,-1,'S1', eta  ,1.-eta,-eta  ),&
&                    QptBound_type (4,-1,'X ', 0.000, 0.000, 0.500),&
&                    QptBound_type (4,-1,'Y ',-zeta , zeta , 0.500),&
&                    QptBound_type (4,-1,'Y1', 0.500, 0.500,-zeta ),&
&                    QptBound_type (4,-1,'Z ', 0.500, 0.500,-0.500) /)
    end if
  else if ((InVar%bravais(1).eq.5).and.(InVar%bravais(2).eq.0)) then
    qptbound_tot=9
    angle_alpha=Lattice%angle_alpha
    eta=(1d0+4*dcos(angle_alpha*pi/180d0))/(2d0+4*dcos(angle_alpha*pi/180d0))
    nu=3d0/4d0-eta/2d0
    ABI_MALLOC(QptBound,(qptbound_tot))
    QptBound(:)=(/ QptBound_type (5, 0,'G ', 0.000, 0.000, 0.000),&
&                  QptBound_type (5, 0,'F ', 0.500, 0.500, 0.000),&
&                  QptBound_type (5, 0,'F1 ', 0.500, 0.000, -0.500),&
&                  QptBound_type (5, 0,'L ', 0.500, 0.000, 0.000),&
&                  QptBound_type (5, 0,'Z ', 0.500, 0.500, 0.500),&
&                  QptBound_type (5, 0,'Q ', 1-nu,nu,0),&
&                  QptBound_type (5, 0,'X ', nu,0,-nu),&
&                  QptBound_type (5, 0,'B1 ', 0.500,1-eta,eta-1),&
&                  QptBound_type (5, 0,'B ', eta,0.500, 1-eta) /)
  else if ((InVar%bravais(1).eq.6).and.(InVar%bravais(2).eq.0)) then
    qptbound_tot=6
    ABI_MALLOC(QptBound,(qptbound_tot))
    QptBound(:)=(/ QptBound_type (6, 0,'G ', 0.000, 0.000, 0.000),&
&                  QptBound_type (6, 0,'A ', 0.000, 0.000, 0.500),&
&                  QptBound_type (6, 0,'H ', 0.333, 0.333, 0.500),&
&                  QptBound_type (6, 0,'K ', 0.333, 0.333, 0.000),&
&                  QptBound_type (6, 0,'L ', 0.500, 0.000, 0.500),&
&                  QptBound_type (6, 0,'M ', 0.500, 0.000, 0.000) /)
  else if ((InVar%bravais(1).eq.7).and.(InVar%bravais(2).eq.0)) then
    qptbound_tot=8
    ABI_MALLOC(QptBound,(qptbound_tot))
    QptBound(:)=(/ QptBound_type (7, 0,'G ', 0.000, 0.000, 0.000),&
&                  QptBound_type (7, 0,'M ', 0.500, 0.500, 0.000),&
&                  QptBound_type (7, 0,'R ', 0.500, 0.500, 0.500),&
&                  QptBound_type (7, 0,'X ', 0.000, 0.500, 0.000),&
! For testing purpose only!!!!!!!!
&                  QptBound_type (7, 0,'A ', 0.000, 1.000, 0.000),&
&                  QptBound_type (7, 0,'B ', 0.500, 1.000, 0.000),&
&                  QptBound_type (7, 0,'C ', 1.000, 1.000, 0.000),&
&                  QptBound_type (7, 0,'D ', 0.750, 0.750, 0.000) /)
  else if ((InVar%bravais(1).eq.7).and.(InVar%bravais(2).eq.-1)) then
    qptbound_tot=4
    ABI_MALLOC(QptBound,(qptbound_tot))
    QptBound(:)=(/ QptBound_type (7,-1,'G ', 0.000, 0.000, 0.000),&
&                  QptBound_type (7,-1,'H ', 0.500,-0.500, 0.500),&
&                  QptBound_type (7,-1,'P ', 0.250, 0.250, 0.250),&
&                  QptBound_type (7,-1,'N ', 0.000, 0.000, 0.500) /)
  else if ((InVar%bravais(1).eq.7).and.(InVar%bravais(2).eq.-3)) then
    qptbound_tot=8
    ABI_MALLOC(QptBound,(qptbound_tot))
    QptBound(:)=(/ QptBound_type (7,-3,'G ', 0.000, 0.000, 0.000),&
&                  QptBound_type (7,-3,'K ', 0.375, 0.375, 0.750),&
&                  QptBound_type (7,-3,'L ', 0.500, 0.500, 0.500),&
&                  QptBound_type (7,-3,'U ', 0.625, 0.250, 0.625),&
&                  QptBound_type (7,-3,'W ', 0.500, 0.250, 0.750),&
&                  QptBound_type (7,-3,'X ', 0.500, 0.000, 0.500),&
&                  QptBound_type (7,-3,'M ', 0.500, 0.500, 0.000),&
&                  QptBound_type (7,-3,'Xp', 0.500, 0.500, 1.000) /)
  end if
  Qpt%qptbound_tot=qptbound_tot

! Define the special Q points USED IN THE CALCULATIONS
! Two cases of generation: default (0) or by hand (>=1)
  if (InVar%BZpath.eq.0) then
    write(InVar%stdout,*) 'Generate the BZ path using the Q points defined by default'
    if (MPIdata%iam_master) then
      write(40,*)         'Generate the BZ path using the Q points defined by default'
    end if  
    if ((InVar%bravais(1).eq.2).and.(InVar%bravais(2).eq.0)) then
!     MONO: G-Y-H-C-E-M1-A-X-H1
!FB      qpt_tot=9
      qpt_tot=5
      ABI_MALLOC(Qpt%special_qpt,(qpt_tot))
      Qpt%special_qpt(1)="X "
      Qpt%special_qpt(2)="G "
      Qpt%special_qpt(3)="Y "
      Qpt%special_qpt(4)="G "
      Qpt%special_qpt(5)="Z "
!FB      Qpt%special_qpt(1)="G "
!FB      Qpt%special_qpt(2)="Y "
!FB      Qpt%special_qpt(3)="H "
!FB      Qpt%special_qpt(4)="C "
!FB      Qpt%special_qpt(5)="E "
!FB      Qpt%special_qpt(6)="M1"
!FB      Qpt%special_qpt(7)="A "
!FB      Qpt%special_qpt(8)="X "
!FB      Qpt%special_qpt(9)="H1"
    else if ((InVar%bravais(1).eq.3).and.(InVar%bravais(2).eq.0)) then
!     ORTH: G-X-S-Y-G-Z
      qpt_tot=6
      ABI_MALLOC(Qpt%special_qpt,(qpt_tot))
      Qpt%special_qpt(1)="G "
      Qpt%special_qpt(2)="X "
      Qpt%special_qpt(3)="S "
      Qpt%special_qpt(4)="Y "
      Qpt%special_qpt(5)="G "
      Qpt%special_qpt(6)="Z "
    else if ((InVar%bravais(1).eq.3).and.(InVar%bravais(2).eq.3)) then
!     ORTH-C: G-Yp-Gp-Z
      qpt_tot=4
      ABI_MALLOC(Qpt%special_qpt,(qpt_tot))
      Qpt%special_qpt(1) ="G "
      Qpt%special_qpt(2) ="Yp"
      Qpt%special_qpt(3) ="Gp"
      Qpt%special_qpt(4) ="Zp"
    else if ((InVar%bravais(1).eq.4).and.(InVar%bravais(2).eq.0)) then
!     TET: G-X-M-G-Z-R-A-Z
      qpt_tot=8
      ABI_MALLOC(Qpt%special_qpt,(qpt_tot))
      Qpt%special_qpt(1) ="G "
      Qpt%special_qpt(2) ="X "
      Qpt%special_qpt(3) ="M "
      Qpt%special_qpt(4) ="G "
      Qpt%special_qpt(5) ="Z "
      Qpt%special_qpt(6) ="R "
      Qpt%special_qpt(7) ="A "
      Qpt%special_qpt(8) ="Z "
    else if ((InVar%bravais(1).eq.4).and.(InVar%bravais(2).eq.-1)) then
      if (Lattice%acell_unitcell(3).lt.Lattice%acell_unitcell(1)) then
!       BCT1: G-X-M-G-Z-P-N-Z1-M
        qpt_tot=9
        ABI_MALLOC(Qpt%special_qpt,(qpt_tot))
        Qpt%special_qpt(1) ="G "
        Qpt%special_qpt(2) ="X "
        Qpt%special_qpt(3) ="M "
        Qpt%special_qpt(4) ="G "
        Qpt%special_qpt(5) ="Z "
        Qpt%special_qpt(6) ="P "
        Qpt%special_qpt(7) ="N "
        Qpt%special_qpt(8) ="Z1"
        Qpt%special_qpt(9) ="M "
      else
!       BCT2: G-X-Y-S-G-Z-S1-N-P-Y1-Z
        qpt_tot=11
        ABI_MALLOC(Qpt%special_qpt,(qpt_tot))
        Qpt%special_qpt(1) ="G "
        Qpt%special_qpt(2) ="X "
        Qpt%special_qpt(3) ="Y "
        Qpt%special_qpt(4) ="S "
        Qpt%special_qpt(5) ="G "
        Qpt%special_qpt(6) ="Z "
        Qpt%special_qpt(7) ="S1"
        Qpt%special_qpt(8) ="N "
        Qpt%special_qpt(9) ="P "
        Qpt%special_qpt(10)="Y1"
        Qpt%special_qpt(11)="Z "
      end if
    else if ((InVar%bravais(1).eq.5).and.(InVar%bravais(2).eq.0)) then
!     RHOMBO:F1-Q-G-Z-B-B1-L-G-F
      qpt_tot=9
      ABI_MALLOC(Qpt%special_qpt,(qpt_tot))
      Qpt%special_qpt(1)="F1"
      Qpt%special_qpt(2)="X"
      Qpt%special_qpt(3)="G "
      Qpt%special_qpt(4)="Z "
      Qpt%special_qpt(5)="B "
      Qpt%special_qpt(6)="B1 "
      Qpt%special_qpt(7)="L "
      Qpt%special_qpt(8)="G "
      Qpt%special_qpt(9)="F "
    else if ((InVar%bravais(1).eq.6).and.(InVar%bravais(2).eq.0)) then
!     HEX: G-M-K-G-A-L-H-A
      qpt_tot=8
      ABI_MALLOC(Qpt%special_qpt,(qpt_tot))
      Qpt%special_qpt(1)="G "
      Qpt%special_qpt(2)="M "
      Qpt%special_qpt(3)="K "
      Qpt%special_qpt(4)="G "
      Qpt%special_qpt(5)="A "
      Qpt%special_qpt(6)="L "
      Qpt%special_qpt(7)="H "
      Qpt%special_qpt(8)="A "
    else if ((InVar%bravais(1).eq.7).and.(InVar%bravais(2).eq.0)) then
!     SC: G-X-M-G-R
      qpt_tot=5
      ABI_MALLOC(Qpt%special_qpt,(qpt_tot))
      Qpt%special_qpt(1)="G "
      Qpt%special_qpt(2)="X "
      Qpt%special_qpt(3)="M "
      Qpt%special_qpt(4)="G "
      Qpt%special_qpt(5)="R "
    else if ((InVar%bravais(1).eq.7).and.(InVar%bravais(2).eq.-1)) then
!     BCC: G-P-H-G-N
      qpt_tot=5
      ABI_MALLOC(Qpt%special_qpt,(qpt_tot))
      Qpt%special_qpt(1)="G "
      Qpt%special_qpt(2)="P "
      Qpt%special_qpt(3)="H "
      Qpt%special_qpt(4)="G "
      Qpt%special_qpt(5)="N "
    else if ((InVar%bravais(1).eq.7).and.(InVar%bravais(2).eq.-3)) then
!     FCC: G-X-W-Xp-K-G-L
      qpt_tot=7
      ABI_MALLOC(Qpt%special_qpt,(qpt_tot))
      Qpt%special_qpt(1)="G "
      Qpt%special_qpt(2)="X "
      Qpt%special_qpt(3)="W "
      Qpt%special_qpt(4)="Xp"
      Qpt%special_qpt(5)="K "
      Qpt%special_qpt(6)="G "
      Qpt%special_qpt(7)="L "
    end if
  else if (InVar%BZpath.ge.1) then
    write(InVar%stdout,*) 'Generate the BZ path using the Q points given in the input file'
    if (MPIdata%iam_master) then
      write(40,*)         'Generate the BZ path using the Q points given in the input file'
    end if  
    qpt_tot=InVar%BZpath
    ABI_MALLOC(Qpt%special_qpt,(qpt_tot))
    Qpt%special_qpt(:)=InVar%special_qpt(:)
  end if
  Qpt%qpt_tot=qpt_tot

 end subroutine tdep_make_specialqpt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_make_qptpath(InVar,Lattice,MPIdata,Qpt)

  implicit none
  integer :: ii,jj,kk,nqpt,iqpt,qpt_tot,tmp_int
  logical :: IsThisAllowed
  type(Input_Variables_type),intent(in) :: InVar
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Qpoints_type),intent(out) :: Qpt
  type(MPI_enreg_type), intent(in) :: MPIdata
  type(QptBound_type), allocatable :: QptBound(:)

  nqpt=0
  if (MPIdata%iam_master) open(unit=40,file=trim(InVar%output_prefix)//'qpt.dat')
  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '########################## Q points generation  #############################'
  write(InVar%stdout,*) '#############################################################################'
! Define the special Q points
  call tdep_make_specialqpt(InVar,Lattice,MPIdata,Qpt,QptBound)

! Define the path in the BZ
! Two cases of generation: default (0) or by hand (>=1)
  qpt_tot     =Qpt%qpt_tot
  ABI_MALLOC(Qpt%special_red ,(qpt_tot,3)); Qpt%special_red (:,:)=zero
  ABI_MALLOC(Qpt%special_cart,(qpt_tot,3)); Qpt%special_cart(:,:)=zero

  if (InVar%BZpath.ge.0) then
!   If "by hand", verify that the Letter defining the Qpt bound is allowed for this
!   crystallographic group
    do jj=1,qpt_tot
      IsThisAllowed=.false.
      do ii=1,Qpt%qptbound_tot
        if ((QptBound(ii)%ihol.eq.InVar%bravais(1)).and.(QptBound(ii)%center.eq.InVar%bravais(2))) then
          if (QptBound(ii)%letter.eq.Qpt%special_qpt(jj)) then
            IsThisAllowed=.true.
            Qpt%special_red(jj,1)=QptBound(ii)%x
            Qpt%special_red(jj,2)=QptBound(ii)%y
            Qpt%special_red(jj,3)=QptBound(ii)%z
          end if
        end if
      end do
      if (.not.IsThisAllowed) then
        MSG_ERROR('One of the Qpt bound (letter) is not allowed.')
      end if
    end do
!   Compute the cartesian coordinates of the special Q points in the reciprocical lattice
    do ii=1,qpt_tot
      do jj=1,3
        do kk=1,3
          if (Lattice%line==2) then
            Qpt%special_cart(ii,jj)=Qpt%special_cart(ii,jj)+Lattice%gprimt(kk,jj)*Qpt%special_red(ii,kk)/Lattice%acell_unitcell(jj)
          else if (Lattice%line==0.or.Lattice%line==1) then
            Qpt%special_cart(ii,jj)=Qpt%special_cart(ii,jj)+Lattice%gprimt(kk,jj)*Qpt%special_red(ii,kk)/Lattice%acell_unitcell(kk)
          end if
        end do
      end do
    end do

    if (qpt_tot.gt.1) then
      ABI_MALLOC(Qpt%lgth_segments,(qpt_tot-1)); Qpt%lgth_segments(:)=0
      do ii=1,qpt_tot-1
        Qpt%lgth_segments(ii)=int(dsqrt((Qpt%special_cart(ii,1)-Qpt%special_cart(ii+1,1))**2+&
&                                       (Qpt%special_cart(ii,2)-Qpt%special_cart(ii+1,2))**2+&
&                                       (Qpt%special_cart(ii,3)-Qpt%special_cart(ii+1,3))**2)*100*2*pi)
      end do

      tmp_int=Qpt%lgth_segments(1)
      do ii=1,qpt_tot-1
        if (InVar%BZlength.eq.0) then
          Qpt%lgth_segments(ii)=int(real(Qpt%lgth_segments(ii))/real(tmp_int)*100)
        else if (InVar%BZlength.gt.0) then
          Qpt%lgth_segments(ii)=InVar%lgth_segments(ii)
        else if (InVar%BZlength.lt.0) then
          Qpt%lgth_segments(ii)=0
        end if
      end do

!     Allocate and define the qpt points along the segments
      do ii=1,qpt_tot-1
        nqpt=nqpt+Qpt%lgth_segments(ii)
      end do
      nqpt=nqpt+1
      ABI_MALLOC(Qpt%qpt_red ,(3,nqpt)); Qpt%qpt_red (:,:)=zero
      ABI_MALLOC(Qpt%qpt_cart,(3,nqpt)); Qpt%qpt_cart(:,:)=zero
      iqpt=0
      do ii=1,qpt_tot-1
        if (Qpt%lgth_segments(ii).eq.0) cycle
        do jj=1,Qpt%lgth_segments(ii)
          iqpt=iqpt+1
          Qpt%qpt_red (:,iqpt)=((jj-1)*Qpt%special_red (ii+1,:)+(Qpt%lgth_segments(ii)-jj+1)*Qpt%special_red (ii,:))&
&           /Qpt%lgth_segments(ii)
          Qpt%qpt_cart(:,iqpt)=((jj-1)*Qpt%special_cart(ii+1,:)+(Qpt%lgth_segments(ii)-jj+1)*Qpt%special_cart(ii,:))&
&           /Qpt%lgth_segments(ii)
        end do
      end do
      Qpt%qpt_red (:,nqpt)=Qpt%special_red (qpt_tot,:)
      Qpt%qpt_cart(:,nqpt)=Qpt%special_cart(qpt_tot,:)
    else if (qpt_tot.eq.1) then
      nqpt=1
      ABI_MALLOC(Qpt%qpt_red ,(3,nqpt)); Qpt%qpt_red (:,:)=zero
      ABI_MALLOC(Qpt%qpt_cart,(3,nqpt)); Qpt%qpt_cart(:,:)=zero
      Qpt%qpt_red (:,1)=Qpt%special_red (1,:)
      Qpt%qpt_cart(:,1)=Qpt%special_cart(1,:)
    end if !qpt_tot.gt.1
    if (MPIdata%iam_master) then
      write(40,*) '  In reduced coordinates:'
      do ii=1,qpt_tot
        write(40,'(a,1x,3(f10.5,1x))') Qpt%special_qpt(ii),Qpt%special_red(ii,1),Qpt%special_red(ii,2),Qpt%special_red(ii,3)
      end do
      write(40,*) ' '
      write(40,*) '  In cartesian coordinates:'
      do ii=1,qpt_tot
        write(40,'(a,1x,3(f10.5,1x))') Qpt%special_qpt(ii),Qpt%special_cart(ii,1),Qpt%special_cart(ii,2),Qpt%special_cart(ii,3)
      end do
      write(40,*) ' '
      write(40,*) '  Using gprimt='
      write(40,'(3(f10.5,1x))') Lattice%gprimt(1,1),Lattice%gprimt(1,2),Lattice%gprimt(1,3)
      write(40,'(3(f10.5,1x))') Lattice%gprimt(2,1),Lattice%gprimt(2,2),Lattice%gprimt(2,3)
      write(40,'(3(f10.5,1x))') Lattice%gprimt(3,1),Lattice%gprimt(3,2),Lattice%gprimt(3,3)
      write(40,*) ' '
      write(40,*) '  The number of points along each direction in the BZ='
      if (qpt_tot.gt.1) then
        do ii=1,qpt_tot-1
          write(40,'(a2,a,a2,1x,i4)') Qpt%special_qpt(ii),'-',Qpt%special_qpt(ii+1),Qpt%lgth_segments(ii)
        end do
      end if	
    end if

  else
    write(InVar%stdout,*) 'The Q points path is defined in the input file'
    if (MPIdata%iam_master) write(40,*)           'The Q points path is defined in the input file'
    nqpt=abs(InVar%BZpath)
    ABI_MALLOC(Qpt%qpt_red ,(3,nqpt)); Qpt%qpt_red (:,:)=zero
    ABI_MALLOC(Qpt%qpt_cart,(3,nqpt)); Qpt%qpt_cart(:,:)=zero
    Qpt%qpt_red(:,:)=InVar%qpt(:,:)
    MSG_ERROR('The indices in the loop below are not consistent')
    do ii=1,nqpt
      do jj=1,3
        do kk=1,3
          if (Lattice%line==2) then
            Qpt%qpt_cart(ii,jj)=Qpt%qpt_cart(ii,jj)+Lattice%gprimt(kk,jj)*Qpt%qpt_red(ii,kk)/Lattice%acell_unitcell(jj)
          else if (Lattice%line==0.or.Lattice%line==1) then
            Qpt%qpt_cart(ii,jj)=Qpt%qpt_cart(ii,jj)+Lattice%gprimt(kk,jj)*Qpt%qpt_red(ii,kk)/Lattice%acell_unitcell(kk)
          end if
        end do
      end do
    end do
  end if  !BZpath>=0
  write(InVar%stdout,*) 'See the qpt.dat file'

! Write the q-points along the path defined in the BZ, in the qpt.dat file
  if (MPIdata%iam_master) then
    write(40,*) ' '
    write(40,*) '  Q-points path (in reduced coordinates) and (in cartesian coordinates)='
    do  iqpt=1,nqpt
      write(40,'(i4,1x,6(f10.5,1x))') iqpt,Qpt%qpt_red(1,iqpt),Qpt%qpt_red(2,iqpt),Qpt%qpt_red(3,iqpt),&
&       Qpt%qpt_cart(1,iqpt),Qpt%qpt_cart(2,iqpt),Qpt%qpt_cart(3,iqpt)
    end do
    close(40)
  end if
  Qpt%nqpt=nqpt
 end subroutine tdep_make_qptpath
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module m_tdep_qpt
