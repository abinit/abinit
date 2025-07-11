
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_abitypes

  use defs_basis
  use m_abicore
  use m_nctk
  use m_xmpi
  use m_errors
  use netcdf
  use m_geometry,         only : xred2xcart
  use m_dynmat,           only : asrif9, d2cart_to_red
  use m_tdep_readwrite,   only : Input_type, MPI_enreg_type
  use m_tdep_latt,        only : Lattice_type
  use m_tdep_phi2,        only : Eigen_type, tdep_build_phi2_33, Phi2_type
  use m_tdep_qpt,         only : Qpoints_type
  use m_tdep_sym,         only : Symetries_type
  use m_tdep_shell,       only : Shell_type
  use m_ifc,              only : ifc_type
  use m_crystal,          only : crystal_t
  use m_ddb,              only : ddb_type
  use m_ddb_hdr,          only : ddb_hdr_type
  use m_dynmat,           only : gtdyn9, d2cart_to_red
  use m_kpts,             only : smpbz, kpts_ibz_from_kptrlatt
  use m_copy,             only : alloc_copy
  use m_symkpt,           only : symkpt
  use m_matrix,           only : mat33det

  implicit none

  type Qbz_type

    integer :: nqbz,nqibz
    integer, allocatable :: ibz2bz(:)
    double precision, allocatable :: wtq(:),wtqibz(:)
    double precision, allocatable :: qbz(:,:), qbz_cart(:,:)
    double precision, allocatable :: qibz(:,:),qibz_cart(:,:)

  end type Qbz_type

  public :: tdep_init_crystal
  public :: tdep_init_ifc
  public :: tdep_read_ifc
  public :: tdep_write_ifc
  public :: tdep_ifc2phi2
  public :: tdep_init_ddb
  public :: tdep_write_ddb
  public :: tdep_destroy_qbz

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_init_crystal(Crystal,Invar,Lattice,Sym)

  implicit none
  type(crystal_t),intent(out) :: Crystal
  type(Input_type),intent(in) :: Invar
  type(Lattice_type),intent(in) :: Lattice
  type(Symetries_type),intent(in) :: Sym

  integer :: timrev,npsp
  logical :: remove_inv,use_antiferro
  double precision, allocatable :: zion(:),znucl(:)
  character(len=132), allocatable :: title(:)

! Initialize the Crystal datatype
! ===============================
  timrev=1
  use_antiferro=.false.
  remove_inv=.false.
  npsp=Invar%ntypat
  ABI_MALLOC(title,(Invar%ntypat))
  ABI_MALLOC(zion,(Invar%ntypat)) ; zion(:)=one
  ABI_MALLOC(znucl,(Invar%ntypat)) ; znucl(:)=one
  if (allocated(Invar%znucl)) znucl = Invar%znucl
  call crystal%init(Invar%amu,Sym%spgroup,Invar%natom_unitcell,npsp,&
&   Invar%ntypat,Sym%nsym,Lattice%rprimdt,Invar%typat_unitcell,Sym%xred_zero,&
&  zion,znucl,timrev,use_antiferro,remove_inv,title,&
!BUG&  Sym%ptsymrel(:,:,1:Sym%nsym),Sym%tnons(:,1:Sym%nsym),Sym%symafm(1:Sym%nsym)) ! Optional
&  Sym%symrel(:,:,1:Sym%nsym),Sym%tnons(:,1:Sym%nsym),Sym%symafm(1:Sym%nsym)) ! Optional
  ABI_FREE(title)
  ABI_FREE(znucl)
  ABI_FREE(zion)

 end subroutine tdep_init_crystal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine tdep_init_ifc(Crystal,DDB,Ifc,Invar,Lattice,MPIdata,Phi2,Rlatt_cart,Shell2at,Sym)

  implicit none
  type(crystal_t),intent(in) :: Crystal
  type(ifc_type),intent(out) :: Ifc
  type(ddb_type),intent(in) :: DDB
  type(Input_type),intent(in) :: Invar
  type(Lattice_type),intent(in) :: Lattice
  type(Symetries_type),intent(in) :: Sym
  type(Shell_type),intent(in) :: Shell2at
  type(MPI_enreg_type), intent(in) :: MPIdata
  type(Phi2_type),intent(inout) :: Phi2
  double precision,intent(in) :: Rlatt_cart(3,Invar%natom_unitcell,Invar%natom)

  integer :: iatcell,dipdip,ii,asr,symdynmat,rfmeth
  integer :: iatom,nsphere,prtsrlr,enunit,nqshft
  integer :: ngqpt_in(3)
  double precision :: rifcsph
  double precision :: dielt(3,3)
  double precision, allocatable :: zeff(:,:,:)
  double precision, allocatable :: q1shft(:,:)
  real(dp),allocatable :: qdrp_cart(:,:,:,:)

! Define matrices for LO-TO
! =========================
  ABI_MALLOC(zeff, (3,3,Invar%natom_unitcell)); zeff (:,:,:)=zero
  dielt(:,:)=    zero
  dielt(1,1)=    1.d0
  dielt(2,2)=    1.d0
  dielt(3,3)=    1.d0
  if (Invar%loto) then
    dipdip=1
    do iatcell=1,Invar%natom_unitcell
      do ii=1,3
        zeff(ii,ii,iatcell)=Invar%born_charge(Invar%typat_unitcell(iatcell))
        dielt(ii,ii)       =Invar%dielec_constant
      end do
    end do
  else
    dipdip=0
  end if

! Initialize the Ifc datatype
! ===========================
  asr=         1
  symdynmat=   1
  rfmeth=      1
  nsphere=     0
!LOTO Seems to not change strongly the results.
!LOTO  rifcsph=     Invar%rcut
  rifcsph=     0
  prtsrlr=     0
  enunit=      0
  ngqpt_in(:)= Invar%ngqpt1(:)
  nqshft=      1
  ABI_MALLOC(qdrp_cart,(3,3,3,Invar%natom_unitcell))
  ABI_MALLOC(q1shft,(3,nqshft))

!LOTO Keep the Gamma point in the q-point mesh
  if (Invar%loto) then
    q1shft(:,:)=0.0d0
  else
    q1shft(:,:)=0.5d0
  end if

  call ifc%init(Crystal,DDB,Lattice%brav,asr,symdynmat,dipdip,&
!LOTO Keep the correct definition of the Lattice
!LOTO  call ifc_init(Ifc,Crystal,DDB,1,asr,symdynmat,dipdip,&
  rfmeth,ngqpt_in,nqshft,q1shft,dielt,zeff,qdrp_cart,nsphere,rifcsph,&
  prtsrlr,enunit,XMPI_WORLD, prtout=.false.)

  ABI_FREE(q1shft)
  ABI_FREE(qdrp_cart)
  ABI_FREE(zeff)

! Read an IFC from ifc_in.dat input file, write it in the ifc_check.dat file and copy to Phi2
! =================================================================================
  if (Invar%readifc.eq.1) then
    write(Invar%stdout,*) ' '
    write(Invar%stdout,*) '#############################################################################'
    write(Invar%stdout,*) '################ Read IFCs from input file ##################################'
    write(Invar%stdout,*) '#############################################################################'
!   Read IFC from ifc_in.dat (readifc=1)
    call tdep_read_ifc(Ifc,Invar,Invar%natom_unitcell)
!   Copy Ifc%atmfrc to Phi2
    call tdep_ifc2phi2(Ifc%dipdip,Ifc,Invar,Lattice,Invar%natom_unitcell,1,Phi2,Rlatt_cart,Shell2at,Sym)
!   Copy Phi2 to Ifc%atmfrc
    call tdep_ifc2phi2(Ifc%dipdip,Ifc,Invar,Lattice,Invar%natom_unitcell,0,Phi2,Rlatt_cart,Shell2at,Sym)
!   Write IFC in ifc_check.dat (for check)
    if (MPIdata%iam_master) call tdep_write_ifc(Crystal,Ifc,Invar,Invar%natom_unitcell,1)

!   Write the Phi2.dat file
    if (Invar%debug.and.MPIdata%iam_master) then
      write(Invar%stdout,'(a)') ' See the Phi2.dat file corresponding to the ifc_in.dat file'
      open(unit=55,file=trim(Invar%output_prefix)//'_Phi2.dat')
      do iatom=1,3*Invar%natom
        write(55,'(10000(f10.6,1x))') Phi2%SR(iatom,:)
      end do
      close(55)
    end if
  end if

!TMP!LOTO
!TMP! If the IFC is read (as above), we assume that there is no residual contribution of the
!TMP! "LR part" within or that the decomposition between "LR" and "SR" parts is done (as it is
!TMP! performed in the ABINIT output of IFC).
!TMP! ============================================================================================
!TMP  if (Invar%loto) then
!TMP    if (MPIdata%iam_master) call tdep_write_ifc(Crystal,Ifc,Invar,Invar%natom_unitcell,1)
!TMP    call tdep_ifc2phi2(Ifc%dipdip,Ifc,Invar,Lattice,Invar%natom_unitcell,1,Phi2,Rlatt_cart,Shell2at,Sym)
!TMP  end if

 end subroutine tdep_init_ifc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine tdep_init_ddb(Crystal,DDB,Invar,Lattice,MPIdata,Qbz)

  implicit none
  type(crystal_t),intent(in) :: Crystal
  type(ddb_type),intent(out) :: DDB
  type(Qbz_type),intent(out) :: Qbz
  type(MPI_enreg_type), intent(in) :: MPIdata
  type(Input_type),intent(in) :: Invar
  type(Lattice_type),intent(in) :: Lattice

  integer :: option,nqshft,mqpt,nqbz,iq_ibz
  integer :: mpert,msize
  integer :: qptrlatt(3,3),ngqpt_in(3)
  integer, allocatable :: bz2ibz_smap(:,:)
  double precision, allocatable :: spqpt(:,:)
  double precision, allocatable :: q1shft(:,:)
  double precision, allocatable :: wtq_folded(:)

! Compute the number of points in the Brillouin zone: nqbz and qbz
! ================================================================
  option=1
  nqshft=1
  ngqpt_in(:)=Invar%ngqpt1(:)
  qptrlatt(:,:)=0
  qptrlatt(1,1)=ngqpt_in(1)
  qptrlatt(2,2)=ngqpt_in(2)
  qptrlatt(3,3)=ngqpt_in(3)
  mqpt=ngqpt_in(1)*ngqpt_in(2)*ngqpt_in(3)*nqshft
  if(Lattice%brav==2)mqpt=mqpt/2
  if(Lattice%brav==3)mqpt=mqpt/4
  ABI_MALLOC(spqpt,(3,mqpt)); spqpt(:,:)=zero
  ABI_MALLOC(q1shft,(3,nqshft))
!LOTO Keep the Gamma point in the q-point mesh
  if (Invar%loto) then
    q1shft(:,:)=0.0d0
  else
    q1shft(:,:)=0.5d0
  end if
!LOTO Keep the correct definition of the Lattice
!LOTO  call smpbz(1,Invar%stdlog,qptrlatt,mqpt,nqbz,nqshft,option,q1shft,spqpt)
  call smpbz(Lattice%brav,Invar%stdlog,qptrlatt,mqpt,nqbz,nqshft,option,q1shft,spqpt)
  ABI_FREE(q1shft)

! Initialize the DDB datatype (only needed values)
! ================================================
!! flg(nsize,nqbz)= flag of existence for each element of the DDB:
!!                   integer,intent(in) :: blkflg(3,mpert,3,mpert,nqbz)
!! nrm(1,nqbz)=norm of qpt providing normalization
!!                   real(dp),intent(in) :: blknrm(3,nqbz)
!! qpt(1<ii<9,nqbz)=q vector of a phonon mode (ii=1,2,3)
!!                   real(dp),intent(in) :: blkqpt(9,nqbz)
!! typ(nqbz)=1 or 2 depending on non-stationary or stationary block 3 for third order derivatives
!!                   integer,intent(in) :: blktyp(nqbz)
!! val(2,3*mpert*3*mpert,nqbz)= all the dynamical matrices
!!                   real(dp),intent(in) :: blkval(2,3*mpert*3*mpert,nqbz)
!! nqbz=number of blocks in the DDB
  mpert=Invar%natom_unitcell+6
  msize=3*mpert*3*mpert
  DDB%mpert = mpert
  DDB%msize = msize
  ABI_MALLOC(DDB%flg,(msize,nqbz))  ; DDB%flg(:,:)=1
  ABI_MALLOC(DDB%nrm,(3,nqbz))      ; DDB%nrm(:,:)  =zero
  ABI_MALLOC(DDB%qpt,(9,nqbz))      ; DDB%qpt(:,:)  =zero
  ABI_MALLOC(DDB%typ,(nqbz))        ; DDB%typ(:)    =1
  ABI_MALLOC(DDB%val,(2,msize,nqbz)); DDB%val(:,:,:)=zero
  DDB%nblok=nqbz
  DDB%nrm(1,:)=1.d0
  DDB%qpt(1:3,:)=spqpt(:,:)
  ABI_FREE(spqpt)
  call alloc_copy(Invar%amu,DDB%amu)
  DDB%rprim=Lattice%rprimt
  DDB%gprim=Lattice%gprim
  DDB%acell=Lattice%acell_unitcell

! Initialize useful quantities stored in the Qbz datatype
! =======================================================
  Qbz%nqbz=nqbz
  ABI_MALLOC(Qbz%qbz     ,(3,nqbz)); Qbz%qbz     (:,:)=zero
  ABI_MALLOC(Qbz%qbz_cart,(3,nqbz)); Qbz%qbz_cart(:,:)=zero
  Qbz%qbz(:,:)=DDB%qpt(1:3,:)
  do iq_ibz=1,Qbz%nqbz
    Qbz%qbz_cart(:,iq_ibz)=matmul(Crystal%gprimd,Qbz%qbz(:,iq_ibz))
  end do

! Reduce the number of such points by symmetrization.
  ABI_MALLOC(Qbz%ibz2bz,(nqbz)); Qbz%ibz2bz=zero
  ABI_MALLOC(Qbz%wtq,(nqbz)); Qbz%wtq=zero
  ABI_MALLOC(wtq_folded,(nqbz)); wtq_folded=zero
  ABI_MALLOC(bz2ibz_smap,(6,nqbz)); wtq_folded=zero
  Qbz%wtq(:)=one/nqbz         ! Weights sum up to one
! Set nsym=1 in order to compute this quantity in the full BZ.
  call symkpt(0,Crystal%gmet,Qbz%ibz2bz,Invar%stdlog,Qbz%qbz,Qbz%nqbz,Qbz%nqibz,&
&             1,Crystal%symrec,1,Qbz%wtq,wtq_folded,bz2ibz_smap,xmpi_comm_self)
  ABI_MALLOC(Qbz%wtqibz   ,(Qbz%nqibz))
  ABI_MALLOC(Qbz%qibz     ,(3,Qbz%nqibz))
  ABI_MALLOC(Qbz%qibz_cart,(3,Qbz%nqibz))
  do iq_ibz=1,Qbz%nqibz
    Qbz%wtqibz(iq_ibz)=wtq_folded(Qbz%ibz2bz(iq_ibz))
    Qbz%qibz(:,iq_ibz)=Qbz%qbz(:,Qbz%ibz2bz(iq_ibz))
    Qbz%qibz_cart(:,iq_ibz)=matmul(Crystal%gprimd,Qbz%qibz(:,iq_ibz))
  end do
  if (MPIdata%iam_master) then
    open(unit=40,file=trim(Invar%output_prefix)//'_qbz.dat')
    open(unit=41,file=trim(Invar%output_prefix)//'_iqbz.dat')
    do iq_ibz=1,Qbz%nqbz
      write(40,'(i4,7(1x,f10.6))') iq_ibz,Qbz%qbz(1:3,iq_ibz),Qbz%wtq(iq_ibz),Qbz%qbz_cart(1:3,iq_ibz)
    end do
    do iq_ibz=1,Qbz%nqibz
      write(41,'(i4,7(1x,f10.6))') iq_ibz,Qbz%qibz(1:3,iq_ibz),Qbz%wtqibz(iq_ibz),Qbz%qibz_cart(1:3,iq_ibz)
    end do
    close(40)
    close(41)
  end if
  ABI_FREE(wtq_folded)
  ABI_FREE(bz2ibz_smap)

 end subroutine tdep_init_ddb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine tdep_write_ddb(DDB,Crystal,Invar)

  implicit none
  type(crystal_t),intent(inout) :: Crystal
  type(ddb_type),intent(inout) :: DDB
  type(Input_type),intent(in) :: Invar

  character (len=fnlen):: filename
  type(ddb_hdr_type) :: ddb_hdr

  ! Initialize ddb_hdr manually
  !write(Invar%stdout,'(a)') ' Writing DDB'
  call ddb_hdr%init_from_crystal(Crystal)
  ddb_hdr%dscrpt = "DDB from TDEP"

  ddb_hdr%mpert = DDB%mpert
  ddb_hdr%nblok = DDB%nblok
  ABI_MALLOC(ddb_hdr%typ,(ddb_hdr%nblok))
  ddb_hdr%typ = DDB%typ

  filename = trim(Invar%output_prefix)//'_DDB'
  call DDB%write(ddb_hdr, filename)
  call ddb_hdr%free()

 end subroutine tdep_write_ddb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tdep_read_ifc(Ifc,Invar,natom_unitcell)

  implicit none

  integer,intent(in) :: natom_unitcell
  type(Input_type),intent(in) :: Invar
  type(ifc_type),intent(inout) :: Ifc

  integer :: iatcell,ifcout,ii,irpt,jatcell,jj,mu
  double precision :: atmfrctot1,atmfrctot2,atmfrctot3,dist
  double precision :: atmfrclr1,atmfrclr2,atmfrclr3
  double precision :: atmfrcsr1,atmfrcsr2,atmfrcsr3
  character (len=25):: string,string1,string2,string3

! Read IFC from ifc_in.dat or ifc_out.dat (according to readifc)
! =======================================================
  ifcout=min(Ifc%nrpt,200)
  Ifc%atmfrc(:,:,:,:,:)=zero
  if (Invar%loto) then
    Ifc%short_atmfrc(:,:,:,:,:)=zero
    Ifc%ewald_atmfrc(:,:,:,:,:)=zero
  end if
  if (Invar%readifc==1) then
    write(Invar%stdout,'(a)') ' Read the IFC from ifc_in.dat'
    open(unit=40,file=trim(Invar%input_prefix)//'_ifc_in.dat')
  else if (Invar%readifc==2) then
    write(Invar%stdout,'(a)') ' Read the IFC from ifc_out.dat'
    open(unit=40,file=trim(Invar%output_prefix)//'_ifc_out.dat')
  end if
  read(40,*) string
  read(40,*) string
  read(40,*) string
  read(40,*) string
  read(40,*) string
  read(40,*) string
  read(40,*) string
  read(40,*) string
  read(40,*) string
  if (Invar%loto) then
    read(40,*) string
    read(40,*) string
    read(40,*) string
  end if
  do iatcell=1,natom_unitcell
    read(40,*) string
    read(40,*) string
    read(40,*) string
    read(40,*) string
    do ii=1,ifcout
      read(40,*) jj,string1,string2,string3,jatcell,string,irpt
      read(40,*) string
      read(40,*) string1,string2,dist
      do mu=1,3
        if (Invar%loto) then
          read(40,*) atmfrctot1,atmfrctot2,atmfrctot3,atmfrclr1,atmfrclr2,atmfrclr3,atmfrcsr1,atmfrcsr2,atmfrcsr3
        else
          read(40,*) atmfrcsr1,atmfrcsr2,atmfrcsr3
        end if
        if (dist.lt.Invar%rcut*0.99) then
!FB       See ifc_getiaf.F90 (for wghatm)
!WARNING : atmfrc != atmfrctot
          Ifc%atmfrc      (1,iatcell,mu,jatcell,irpt)=atmfrcsr1 /Ifc%wghatm(iatcell,jatcell,irpt)
          Ifc%atmfrc      (2,iatcell,mu,jatcell,irpt)=atmfrcsr2 /Ifc%wghatm(iatcell,jatcell,irpt)
          Ifc%atmfrc      (3,iatcell,mu,jatcell,irpt)=atmfrcsr3 /Ifc%wghatm(iatcell,jatcell,irpt)
          if (Invar%loto) then
            Ifc%ewald_atmfrc(1,iatcell,mu,jatcell,irpt)=atmfrclr1 /Ifc%wghatm(iatcell,jatcell,irpt)
            Ifc%ewald_atmfrc(2,iatcell,mu,jatcell,irpt)=atmfrclr2 /Ifc%wghatm(iatcell,jatcell,irpt)
            Ifc%ewald_atmfrc(3,iatcell,mu,jatcell,irpt)=atmfrclr3 /Ifc%wghatm(iatcell,jatcell,irpt)
            Ifc%short_atmfrc(1,iatcell,mu,jatcell,irpt)=atmfrcsr1 /Ifc%wghatm(iatcell,jatcell,irpt)
            Ifc%short_atmfrc(2,iatcell,mu,jatcell,irpt)=atmfrcsr2 /Ifc%wghatm(iatcell,jatcell,irpt)
            Ifc%short_atmfrc(3,iatcell,mu,jatcell,irpt)=atmfrcsr3 /Ifc%wghatm(iatcell,jatcell,irpt)
          end if
        else
          Ifc%atmfrc(1:3,iatcell,mu,jatcell,irpt)=zero
          if (Invar%loto) then
            Ifc%ewald_atmfrc(1,iatcell,mu,jatcell,irpt)=zero
            Ifc%short_atmfrc(1,iatcell,mu,jatcell,irpt)=zero
          end if
        end if
      end do
      read(40,*) string
      read(40,*) string
      read(40,*) string
      read(40,*) string
      read(40,*) string
      read(40,*) string
      read(40,*) string
      read(40,*) string
      if (Invar%loto) then
        read(40,*) string
        read(40,*) string
        read(40,*) string
        read(40,*) string
        read(40,*) string
        read(40,*) string
      end if
    end do
  end do
  close(40)
  write(Invar%stdout,'(a)') ' ------- achieved'
! When the IFC come from other calculations (readifc=1), then impose acoustic sum rule
  call asrif9(Ifc%asr,Ifc%atmfrc,natom_unitcell,Ifc%nrpt,Ifc%rpt,Ifc%wghatm)

end subroutine tdep_read_ifc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tdep_write_ifc(Crystal,Ifc,Invar,natom_unitcell,unitfile)

  implicit none

  integer,intent(in) :: natom_unitcell,unitfile
  type(Input_type),intent(in) :: Invar
  type(ifc_type),intent(inout) :: Ifc
  type(crystal_t),intent(in) :: Crystal

  integer :: ifcana,ifcout,ncid,prt_ifc
  integer :: atifc(natom_unitcell)
  character(len=500) :: message

! Write IFC in ifc_out.dat or ifc_check.dat (according to readifc)
! ================================================================
  ifcana  = 1
  atifc(:)= 1
  ifcout  = min(Ifc%nrpt,200)
  prt_ifc = 1
  if (unitfile.eq.0) then
    write(Invar%stdout,'(a)') ' Write the IFC of TDEP in ifc_out.dat (and ifc_out.nc)'
    open(unit=77,file=trim(Invar%output_prefix)//'_ifc_out.dat')
  else if (unitfile.eq.1) then
    write(Invar%stdout,'(a)') ' Write in ifc_check.dat (and ifc_check.nc) the IFC read previously'
    open(unit=77,file=trim(Invar%output_prefix)//'_ifc_check.dat')
  else if (unitfile.eq.2) then
    write(Invar%stdout,'(a)') ' Write in ifc_ddb.dat (and ifc_ddb.nc) the IFC read from DDB file'
    open(unit=77,file=trim(Invar%output_prefix)//'_ifc_ddb.dat')
  else
    write(message, '(a,i3,a)' )&
&   ' The value of unitfile ',unitfile,' is not allowed.'
    ABI_ERROR(message)
  end if
  if (unitfile.eq.0) then
    NCF_CHECK_MSG(nctk_open_create(ncid, trim(Invar%output_prefix)//"_ifc_out.nc", xmpi_comm_self), "Creating ifc_out.nc")
  else if (unitfile.eq.1) then
    NCF_CHECK_MSG(nctk_open_create(ncid, trim(Invar%output_prefix)//"_ifc_check.nc", xmpi_comm_self), "Creating ifc_check.nc")
  else if (unitfile.eq.2) then
    NCF_CHECK_MSG(nctk_open_create(ncid, trim(Invar%output_prefix)//"_ifc_ddb.nc", xmpi_comm_self), "Creating ifc_ddb.nc")
  end if
  NCF_CHECK(nctk_def_basedims(ncid))
  NCF_CHECK(nctk_defnwrite_ivars(ncid, ["anaddb_version"], [1]))
  NCF_CHECK(crystal%ncwrite(ncid))
  call ifc%write(ifcana,atifc,ifcout,prt_ifc,ncid,Invar%output_prefix,77)
  close(77)
  write(Invar%stdout,'(a)') ' ------- achieved'

end subroutine tdep_write_ifc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tdep_ifc2phi2(dipdip,Ifc,Invar,Lattice,natom_unitcell,option,Phi2,Rlatt4abi,Shell2at,Sym)

  implicit none

  integer,intent(in) :: dipdip,natom_unitcell,option
  type(Input_type),intent(in) :: Invar
  type(Lattice_type),intent(in) :: Lattice
  type(Symetries_type),intent(in) :: Sym
  type(Shell_type),intent(in) :: Shell2at
  type(ifc_type),intent(inout) :: Ifc
  type(Phi2_type), intent(inout) :: Phi2
  double precision, intent(in) :: Rlatt4abi(3,natom_unitcell,Invar%natom)

  integer :: eatom,fatom,iatcell,ii,irpt,jatcell,jatom,jj,isym,kk
  integer :: ishell,iatref,jatref,iatshell,trans
  double precision :: tol,dist
  double precision :: tmp(3)
  double precision :: vect_ddb(3),vect_tdep(3),xcart(3,natom_unitcell)
  double precision,allocatable :: Phi2_33(:,:),Phi2_ref(:,:,:)

! Link the atmfrc to the Phi2
! ==============================
  if (dipdip==0.or.Invar%readifc.ne.0) then
    if (option==0) then
      Ifc%atmfrc      (:,:,:,:,:)=zero
    end if
    if (option==1) then
      Phi2%SR (:,:)=zero
    end if
  end if
  call xred2xcart(natom_unitcell,Lattice%rprimt(:,:),xcart(:,:),Sym%xred_zero(:,:))
  tol=1.d-8

  do iatcell=1,natom_unitcell
    do jatcell=1,natom_unitcell
      do jatom=jatcell,Invar%natom,natom_unitcell
        vect_tdep(:)=Rlatt4abi(:,iatcell,jatom)+xcart(:,jatcell)-xcart(:,iatcell)
        tmp(:)=(vect_tdep(:)*Lattice%acell_unitcell(:))**2
        dist=dsqrt(sum(tmp(:)))
        if (dist.gt.(Invar%rcut*0.99).and.option==1) then
          Phi2%SR ((iatcell-1)*3+1:(iatcell-1)*3+3,(jatom-1)*3+1:(jatom-1)*3+3)=zero
          if (Invar%loto) then
            Phi2%LR ((iatcell-1)*3+1:(iatcell-1)*3+3,(jatom-1)*3+1:(jatom-1)*3+3)=zero
            Phi2%Tot((iatcell-1)*3+1:(iatcell-1)*3+3,(jatom-1)*3+1:(jatom-1)*3+3)=zero
          end if
          cycle
        end if
        do irpt=1,Ifc%nrpt
          if (Ifc%wghatm(iatcell,jatcell,irpt).lt.tol8) cycle
          vect_ddb (:)=Ifc%rpt(:,irpt)+Ifc%rcan(:,jatcell)-Ifc%rcan(:,iatcell)
          tmp(:)=(vect_ddb(:)*Lattice%acell_unitcell(:))**2
          dist=dsqrt(sum(tmp(:)))
          if (dist.gt.(Invar%rcut*0.99).and.option==0) then
            Ifc%atmfrc      (:,iatcell,:,jatcell,irpt)=zero
            if (Invar%loto) then
              Ifc%short_atmfrc(:,iatcell,:,jatcell,irpt)=zero
              Ifc%ewald_atmfrc(:,iatcell,:,jatcell,irpt)=zero
            end if
            cycle
          end if
          if (abs(vect_ddb(1)-vect_tdep(1)).lt.tol.and.&
&             abs(vect_ddb(2)-vect_tdep(2)).lt.tol.and.&
&             abs(vect_ddb(3)-vect_tdep(3)).lt.tol) then
            do ii=1,3
              do jj=1,3
                if (option==0) then
                  if (Invar%loto) then
                    Ifc%atmfrc      (ii,iatcell,jj,jatcell,irpt)=Ifc%atmfrc      (ii,iatcell,jj,jatcell,irpt)+&
&                          Phi2%SR  (ii+(iatcell-1)*3,(jatom-1)*3+jj)/Ifc%wghatm(iatcell,jatcell,irpt)
!LOTO                    Ifc%ewald_atmfrc(ii,iatcell,jj,jatcell,irpt)=Ifc%ewald_atmfrc(ii,iatcell,jj,jatcell,irpt)+&
!LOTO&                          Phi2%LR  (ii+(iatcell-1)*3,(jatom-1)*3+jj)/Ifc%wghatm(iatcell,jatcell,irpt)
                    Ifc%ewald_atmfrc(ii,iatcell,jj,jatcell,irpt)=Ifc%ewald_atmfrc(ii,iatcell,jj,jatcell,irpt)+&
&                          Phi2%LR  (ii+(iatcell-1)*3,(jatom-1)*3+jj)
                    Ifc%short_atmfrc(ii,iatcell,jj,jatcell,irpt)=Ifc%short_atmfrc(ii,iatcell,jj,jatcell,irpt)+&
&                          Phi2%SR  (ii+(iatcell-1)*3,(jatom-1)*3+jj)/Ifc%wghatm(iatcell,jatcell,irpt)
                  else
                    Ifc%atmfrc      (ii,iatcell,jj,jatcell,irpt)=Ifc%atmfrc      (ii,iatcell,jj,jatcell,irpt)+&
&                          Phi2%SR  (ii+(iatcell-1)*3,(jatom-1)*3+jj)/Ifc%wghatm(iatcell,jatcell,irpt)
                  end if
                else if (option==1) then
                  if (Invar%loto) then
!LOTO                    Phi2%LR (ii+(iatcell-1)*3,(jatom-1)*3+jj)=Ifc%ewald_atmfrc(ii,iatcell,jj,jatcell,irpt)*&
!LOTO&                                                               Ifc%wghatm(iatcell,jatcell,irpt)
                    Phi2%LR (ii+(iatcell-1)*3,(jatom-1)*3+jj)=Ifc%ewald_atmfrc(ii,iatcell,jj,jatcell,irpt)
                    Phi2%SR (ii+(iatcell-1)*3,(jatom-1)*3+jj)=Ifc%short_atmfrc(ii,iatcell,jj,jatcell,irpt)*&
&                                                               Ifc%wghatm(iatcell,jatcell,irpt)
                    Phi2%Tot(ii+(iatcell-1)*3,(jatom-1)*3+jj)=Phi2%SR (ii+(iatcell-1)*3,(jatom-1)*3+jj)+&
&                                                             Phi2%LR (ii+(iatcell-1)*3,(jatom-1)*3+jj)
                  else
                    Phi2%SR (ii+(iatcell-1)*3,(jatom-1)*3+jj)=Ifc%atmfrc(ii,iatcell,jj,jatcell,irpt)*&
&                                                               Ifc%wghatm(iatcell,jatcell,irpt)
                  end if
                end if
              enddo
            enddo
          endif
        end do
      end do
    end do
  end do

  if (option==0) call asrif9(Ifc%asr,Ifc%atmfrc,natom_unitcell,Ifc%nrpt,Ifc%rpt,Ifc%wghatm)

!TODO : symetrization of the LR part
  if (option==1) then
!   Build the Phi2
    ABI_MALLOC(Phi2_ref,(3,3,Shell2at%nshell)); Phi2_ref(:,:,:)=zero
    ABI_MALLOC(Phi2_33,(3,3)) ; Phi2_33(:,:)=zero
    do ishell=1,Shell2at%nshell
!     Build the 3x3 IFC per shell
      iatref=Shell2at%iatref(ishell)
      jatref=Shell2at%jatref(ishell)
      Phi2_ref(:,:,ishell)=Phi2%SR((iatref-1)*3+1:(iatref-1)*3+3,(jatref-1)*3+1:(jatref-1)*3+3)
    end do
    Phi2%SR(:,:)=zero
    do ishell=1,Shell2at%nshell
      do eatom=1,Invar%natom
!       Build the 3x3 IFC of an atom in this shell
        if (Shell2at%neighbours(eatom,ishell)%n_interactions.eq.0) cycle
        do iatshell=1,Shell2at%neighbours(eatom,ishell)%n_interactions
          fatom=Shell2at%neighbours(eatom,ishell)%atomj_in_shell(iatshell)
          isym =Shell2at%neighbours(eatom,ishell)%sym_in_shell(iatshell)
          trans=Shell2at%neighbours(eatom,ishell)%transpose_in_shell(iatshell)
          if (fatom.lt.eatom) cycle
          call tdep_build_phi2_33(isym,Phi2_ref(:,:,ishell),Phi2_33,Sym,trans)
!         Symetrization of the Phi2 matrix
          Phi2%SR((eatom-1)*3+1:(eatom-1)*3+3,3*(fatom-1)+1:3*(fatom-1)+3)=Phi2_33(:,:)
          do ii=1,3
            do jj=1,3
              Phi2%SR((fatom  -1)*3+ii,3*(eatom-1)+jj)=Phi2_33(jj,ii)
            end do
          end do
        end do !iatshell
      end do !eatom
    end do !ishell
!   Acoustic sum rule
    do eatom=1,Invar%natom
      Phi2%SR((eatom-1)*3+1:(eatom-1)*3+3,(eatom-1)*3+1:(eatom-1)*3+3)=zero
      do jj=1,3
        do kk=1,3
          do fatom=1,Invar%natom
            if (fatom==eatom) cycle
            Phi2%SR((eatom-1)*3+jj,(eatom-1)*3+kk)=Phi2%SR((eatom-1)*3+jj,3*(eatom-1)+kk)&
&                                                 -Phi2%SR((eatom-1)*3+jj,3*(fatom-1)+kk)
          enddo !fatom
        enddo !kk
      enddo !jj
    end do !eatom
    ABI_FREE(Phi2_33)
    ABI_FREE(Phi2_ref)
  end if

end subroutine tdep_ifc2phi2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine tdep_destroy_qbz(Qbz)

  implicit none
  type(Qbz_type),intent(inout) :: Qbz

  ABI_FREE(Qbz%qbz)
  ABI_FREE(Qbz%qbz_cart)
  ABI_FREE(Qbz%ibz2bz)
  ABI_FREE(Qbz%wtq)
  ABI_FREE(Qbz%wtqibz)
  ABI_FREE(Qbz%qibz)
  ABI_FREE(Qbz%qibz_cart)

 end subroutine tdep_destroy_qbz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module m_tdep_abitypes
