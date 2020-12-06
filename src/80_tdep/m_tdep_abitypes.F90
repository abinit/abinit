
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
#ifdef HAVE_NETCDF
  use netcdf
#endif
  use m_ddb_hdr,          only : ddb_hdr_type
  use m_geometry,         only : xred2xcart
  use m_dynmat,           only : asrif9, d2cart_to_red
  use m_tdep_readwrite,   only : Input_Variables_type, MPI_enreg_type
  use m_tdep_latt,        only : Lattice_Variables_type
  use m_tdep_phi2,        only : Eigen_Variables_type, tdep_build_phi2_33
  use m_tdep_qpt,         only : Qpoints_type
  use m_tdep_sym,         only : Symetries_Variables_type
  use m_tdep_shell,       only : Shell_Variables_type
  use m_ifc,              only : ifc_type, ifc_init, ifc_init_fromFile
  use m_crystal,          only : crystal_t, crystal_init
  use m_ddb,              only : ddb_type
  use m_kpts,             only : smpbz
  use m_copy,             only : alloc_copy
  use m_symkpt,           only : symkpt 

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

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_init_crystal(Crystal,InVar,Lattice,Sym)

  implicit none
  type(crystal_t),intent(out) :: Crystal
  type(Input_Variables_type),intent(in) :: InVar
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Symetries_Variables_type),intent(in) :: Sym

  integer :: timrev,npsp
  logical :: remove_inv,use_antiferro
  double precision, allocatable :: zion(:),znucl(:)
  character(len=132), allocatable :: title(:)

! Initialize the Crystal datatype
! ===============================
  timrev=1
  use_antiferro=.false.
  remove_inv=.false.
  npsp=InVar%ntypat
  ABI_MALLOC(zion,(InVar%ntypat)) ; zion(:)=zero
  ABI_MALLOC(znucl,(npsp)) ; znucl(:)=zero
  ABI_MALLOC(title,(InVar%ntypat))
  call crystal_init(InVar%amu,Crystal,Sym%spgroup,InVar%natom_unitcell,npsp,&
&   InVar%ntypat,Sym%nsym,Lattice%rprimdt,InVar%typat_unitcell,Sym%xred_zero,&
&  zion,znucl,timrev,use_antiferro,remove_inv,title,&
&  Sym%ptsymrel(:,:,1:Sym%nsym),Sym%tnons(:,1:Sym%nsym),Sym%symafm(1:Sym%nsym)) ! Optional
  ABI_FREE(title)
  ABI_FREE(znucl)
  ABI_FREE(zion)

 end subroutine tdep_init_crystal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine tdep_init_ifc(Crystal,DDB,Ifc,InVar,Lattice,MPIdata,Phi2,Rlatt_cart,Shell2at,Sym)

  implicit none
  type(crystal_t),intent(in) :: Crystal
  type(ifc_type),intent(out) :: Ifc
  type(ddb_type),intent(in) :: DDB
  type(Input_Variables_type),intent(in) :: InVar
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Shell_Variables_type),intent(in) :: Shell2at
  type(MPI_enreg_type), intent(in) :: MPIdata
  double precision,intent(out) :: Phi2(3*InVar%natom,3*InVar%natom)
  double precision,intent(in) :: Rlatt_cart(3,InVar%natom_unitcell,InVar%natom)

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
  ABI_MALLOC(zeff, (3,3,InVar%natom_unitcell)); zeff (:,:,:)=zero
  dielt(:,:)=    zero
  dielt(1,1)=    1.d0
  dielt(2,2)=    1.d0
  dielt(3,3)=    1.d0
  if (InVar%loto) then
    dipdip=1
    do iatcell=1,InVar%natom_unitcell
      do ii=1,3
        zeff(ii,ii,iatcell)=InVar%born_charge(InVar%typat_unitcell(iatcell))
        dielt(ii,ii)       =InVar%dielec_constant
      end do
    end do
  else
    dipdip=0
  end if

! Initialize the Ifc datatype
! ===========================
  asr=         1
!  asr=         2
  symdynmat=   1
  rfmeth=      1
  nsphere=     0
  rifcsph=     0
!FB  rifcsph=     InVar%Rcut
!  prtsrlr=     1
  prtsrlr=     0
  enunit=      0

  ngqpt_in(:)= InVar%ngqpt1(:)
  nqshft=      1
  ABI_MALLOC(q1shft,(3,nqshft)); q1shft(:,:)=0.5d0
!FB  ABI_MALLOC(q1shft,(3,nqshft)); q1shft(:,:)=0.0d0
  ABI_ALLOCATE(qdrp_cart,(3,3,3,InVar%natom_unitcell))

  call ifc_init(Ifc,Crystal,DDB,Lattice%brav,asr,symdynmat,dipdip,&
!FB  call ifc_init(Ifc,Crystal,DDB,1,asr,symdynmat,dipdip,&
  rfmeth,ngqpt_in,nqshft,q1shft,dielt,zeff,qdrp_cart,nsphere,rifcsph,&
  prtsrlr,enunit,XMPI_WORLD)

  ABI_DEALLOCATE(qdrp_cart)

  ABI_FREE(q1shft)
  ABI_FREE(zeff)

! Read an IFC from ifc_in.dat input file, write it in the ifc_check.dat file and copy to Phi2
! =================================================================================
  if (InVar%ReadIFC.eq.1) then
    write(InVar%stdout,*) ' '
    write(InVar%stdout,*) '#############################################################################'
    write(InVar%stdout,*) '################ Read IFCs from input file ##################################'
    write(InVar%stdout,*) '#############################################################################'
!   Read IFC from ifc_in.dat (ReadIFC=1)
    call tdep_read_ifc(Ifc,InVar,InVar%natom_unitcell)
!   Copy Ifc%atmfrc to Phi2
    call tdep_ifc2phi2(Ifc%dipdip,Ifc,InVar,Lattice,InVar%natom_unitcell,1,Phi2,Rlatt_cart,Shell2at,Sym)
!   Copy Phi2 to Ifc%atmfrc
    call tdep_ifc2phi2(Ifc%dipdip,Ifc,InVar,Lattice,InVar%natom_unitcell,0,Phi2,Rlatt_cart,Shell2at,Sym)
!   Write IFC in ifc_out.dat (for check)
    if (MPIdata%iam_master) call tdep_write_ifc(Crystal,Ifc,InVar,MPIdata,InVar%natom_unitcell,1)

!   Write the Phi2.dat file
    if (InVar%debug.and.MPIdata%iam_master) then
      write(InVar%stdout,'(a)') ' See the Phi2.dat file corresponding to the ifc_in.dat file'
      open(unit=55,file=trim(InVar%output_prefix)//'Phi2.dat')
      do iatom=1,3*InVar%natom
        write(55,'(10000(f10.6,1x))') Phi2(iatom,:)
      end do
      close(55)
    end if
  end if

 end subroutine tdep_init_ifc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine tdep_init_ddb(Crystal,DDB,InVar,Lattice,MPIdata,Qbz)

  implicit none
  type(crystal_t),intent(in) :: Crystal
  type(ddb_type),intent(out) :: DDB
  type(Qbz_type),intent(out) :: Qbz
  type(MPI_enreg_type), intent(in) :: MPIdata
  type(Input_Variables_type),intent(in) :: InVar
  type(Lattice_Variables_type),intent(in) :: Lattice

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
  ngqpt_in(:)=InVar%ngqpt1(:)
  qptrlatt(:,:)=0
  qptrlatt(1,1)=ngqpt_in(1)
  qptrlatt(2,2)=ngqpt_in(2)
  qptrlatt(3,3)=ngqpt_in(3)
  mqpt=ngqpt_in(1)*ngqpt_in(2)*ngqpt_in(3)*nqshft
  if(Lattice%brav==2)mqpt=mqpt/2
  if(Lattice%brav==3)mqpt=mqpt/4
  ABI_MALLOC(spqpt,(3,mqpt)); spqpt(:,:)=zero
  ABI_MALLOC(q1shft,(3,nqshft)); q1shft(:,:)=0.5d0
!FB  ABI_MALLOC(q1shft,(3,nqshft)); q1shft(:,:)=0.0d0
  call smpbz(Lattice%brav,InVar%stdlog,qptrlatt,mqpt,nqbz,nqshft,option,q1shft,spqpt)
!FB  call smpbz(1,InVar%stdlog,qptrlatt,mqpt,nqbz,nqshft,option,q1shft,spqpt)
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
  mpert=InVar%natom_unitcell+6
  msize=3*mpert*3*mpert
  ABI_MALLOC(DDB%flg,(msize,nqbz))  ; DDB%flg(:,:)=1
  ABI_MALLOC(DDB%nrm,(3,nqbz))      ; DDB%nrm(:,:)  =zero
  ABI_MALLOC(DDB%qpt,(9,nqbz))      ; DDB%qpt(:,:)  =zero
  ABI_MALLOC(DDB%typ,(nqbz))        ; DDB%typ(:)    =1
  ABI_MALLOC(DDB%val,(2,msize,nqbz)); DDB%val(:,:,:)=zero
  DDB%nblok=nqbz
  DDB%nrm(1,:)=1.d0
  DDB%qpt(1:3,:)=spqpt(:,:)
  ABI_FREE(spqpt)
  call alloc_copy(InVar%amu,DDB%amu)
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
  ABI_ALLOCATE(Qbz%ibz2bz,(nqbz)); Qbz%ibz2bz=zero
  ABI_ALLOCATE(Qbz%wtq,(nqbz)); Qbz%wtq=zero
  ABI_ALLOCATE(wtq_folded,(nqbz)); wtq_folded=zero
  ABI_ALLOCATE(bz2ibz_smap,(6,nqbz)); wtq_folded=zero
  Qbz%wtq(:)=one/nqbz         ! Weights sum up to one
! Set nsym=1 in order to compute this quantity in the full BZ. 
  call symkpt(0,Crystal%gmet,Qbz%ibz2bz,InVar%stdlog,Qbz%qbz,Qbz%nqbz,Qbz%nqibz,&
&             1,Crystal%symrec,1,Qbz%wtq,wtq_folded,bz2ibz_smap,xmpi_comm_self)
  ABI_ALLOCATE(Qbz%wtqibz   ,(Qbz%nqibz))
  ABI_ALLOCATE(Qbz%qibz     ,(3,Qbz%nqibz))
  ABI_ALLOCATE(Qbz%qibz_cart,(3,Qbz%nqibz))
  do iq_ibz=1,Qbz%nqibz
    Qbz%wtqibz(iq_ibz)=wtq_folded(Qbz%ibz2bz(iq_ibz))
    Qbz%qibz(:,iq_ibz)=Qbz%qbz(:,Qbz%ibz2bz(iq_ibz))
    Qbz%qibz_cart(:,iq_ibz)=matmul(Crystal%gprimd,Qbz%qibz(:,iq_ibz))
  end do
  if (MPIdata%iam_master) then
    open(unit=40,file='qbz.dat')
    open(unit=41,file='iqbz.dat')
    do iq_ibz=1,Qbz%nqbz
      write(40,'(i4,7(1x,f10.6))') iq_ibz,Qbz%qbz(1:3,iq_ibz),Qbz%wtq(iq_ibz),Qbz%qbz_cart(1:3,iq_ibz)
    end do
    do iq_ibz=1,Qbz%nqibz
      write(41,'(i4,7(1x,f10.6))') iq_ibz,Qbz%qibz(1:3,iq_ibz),Qbz%wtqibz(iq_ibz),Qbz%qibz_cart(1:3,iq_ibz)
    end do
    close(40)
    close(41)
  end if  
  ABI_DEALLOCATE(wtq_folded)
  ABI_DEALLOCATE(bz2ibz_smap)

 end subroutine tdep_init_ddb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine tdep_write_ddb(Crystal,DDB,Eigen2nd,InVar,Lattice,MPIdata,Qbz,Sym)

  implicit none
  type(crystal_t),intent(in) :: Crystal
  type(ddb_type),intent(inout) :: DDB
  type(Eigen_Variables_type),intent(in) :: Eigen2nd
  type(Input_Variables_type),intent(in) :: InVar
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Qbz_type),intent(in) :: Qbz
  type(Symetries_Variables_type),intent(in) :: Sym
  type(MPI_enreg_type), intent(in) :: MPIdata

  integer :: choice,mband,mpert,msize,unddb,iq_ibz,iatom
  integer :: jatom,ii,jj,nblok,counter
  type(ddb_hdr_type) :: ddb_hdr
  character(len=fnlen) :: dscrpt,filename,filename_bis
  character(len=500) :: message

! Test !!!!!!!!!!!!!
  integer :: natom_uc,comm,unitfile
  integer :: nqshft
  integer :: ngqpt(3)
  real(dp) :: dielt(3,3)
  real(dp),allocatable :: zeff(:,:,:)
  real(dp),allocatable :: qshft(:,:)
  real(dp),allocatable :: d2cart(:,:,:,:,:),d2red(:,:,:,:,:)
  real(dp),allocatable :: qdrp_cart(:,:,:,:)
  type(ifc_type) :: Ifc_tmp
  type(crystal_t) :: Crystal_tmp

!FB  nblok=Qbz%nqibz
  nblok=Qbz%nqbz
  choice=2
  mband=1
!FB  mpert=InVar%natom_unitcell+6
  mpert=InVar%natom_unitcell
  msize=3*mpert*3*mpert

!FB  write(6,*) 'BEFORE INIT VAL'
  ABI_MALLOC(d2cart,(2,3,InVar%natom_unitcell,3,InVar%natom_unitcell));d2cart=zero
  ABI_MALLOC(d2red ,(2,3,InVar%natom_unitcell,3,InVar%natom_unitcell));d2red =zero
! Fill the DDB%val and DDB%flg datatype
  ddb%flg(:,:)=zero
!FB  write(6,*) ' In cartesian coordinates'
  do iq_ibz=1,Qbz%nqibz
!FB  do iq_ibz=1,Qbz%nqbz
!FB    write(6,*) 'qpt=',iq_ibz
    d2cart=zero ; d2red=zero
    do jatom=1,InVar%natom_unitcell
      do jj=1,3
        do iatom=1,InVar%natom_unitcell
          do ii=1,3
            counter=counter+1
            d2cart(1,ii,iatom,jj,iatom) =Eigen2nd%dynmat(1,ii,iatom,jj,jatom,iq_ibz)
            d2cart(2,ii,iatom,jj,iatom) =Eigen2nd%dynmat(2,ii,iatom,jj,jatom,iq_ibz)
!FB            write(6,'(a,4i4,2d22.14)')'In cart :',ii,iatom,jj,jatom,d2cart(1,ii,iatom,jj,iatom),d2cart(2,ii,iatom,jj,iatom)
          end do
        end do
      end do
    end do
!FB    write(6,*) ' '
    call d2cart_to_red(d2cart,d2red,Crystal%gprimd,Crystal%rprimd,mpert, &
&                      InVar%natom_unitcell,InVar%ntypat,Crystal%typat,Crystal%ucvol,Crystal%zion)
    counter=0
    do jatom=1,InVar%natom_unitcell
      do jj=1,3
        do iatom=1,InVar%natom_unitcell
          do ii=1,3
            counter=counter+1
            DDB%val(1,counter,iq_ibz) =d2red(1,ii,iatom,jj,jatom)
            DDB%val(2,counter,iq_ibz) =d2red(2,ii,iatom,jj,jatom)
            DDB%flg(counter,iq_ibz)=1
!FB            write(6,'(a,4i4,2d22.14)')'In red :',ii,iatom,jj,jatom,DDB%val(1,counter,iq_ibz),DDB%val(2,counter,iq_ibz)
          end do
        end do
      end do
    end do
  end do
  ABI_FREE(d2cart)
  ABI_FREE(d2red)

!FB  write(6,*) 'BEFORE ERROR'
! Print the header of the DDB
  unddb=16
  filename = trim(InVar%output_prefix)//'_DDB'
  filename_bis = 'input.DDB'
  ddb_hdr%dscrpt = "DDB from TDEP"
  ddb_hdr%nblok = nblok
  ddb_hdr%mblktyp = choice

  ddb_hdr%matom=20
  if (InVar%natom_unitcell.gt.ddb_hdr%matom) then
    write(message, '(a,i3,a,i3,a)' )&
&   ' natom ',InVar%natom_unitcell,' is greater than matom ',ddb_hdr%matom,' !'
    MSG_ERROR(message)
  end if
  ddb_hdr%mband=mband
  ddb_hdr%mkpt=1
  ddb_hdr%msym=48
  if (Sym%nsym.gt.ddb_hdr%msym) then
    write(message, '(a,i3,a,i3,a)' )&
&   ' nsym ',Sym%nsym,' is greater than msym ',ddb_hdr%msym,' !'
    MSG_ERROR(message)
  end if
  ddb_hdr%mtypat=5
  if (InVar%ntypat.gt.ddb_hdr%mtypat) then
    write(message, '(a,i3,a,i3,a)' )&
&   ' nsym ',InVar%ntypat,' is greater than mtypat ',ddb_hdr%mtypat,' !'
    MSG_ERROR(message)
  end if
  ddb_hdr%ddb_version=100401

!FB  write(6,*) 'BEFORE CONSTANT'
  ddb_hdr%dilatmx=0
  ddb_hdr%ecut=0
  ddb_hdr%ecutsm=0
  ddb_hdr%intxc=0
  ddb_hdr%iscf=0
  ddb_hdr%ixc=0
  ddb_hdr%kptnrm=0
  ddb_hdr%ngfft=1
  ddb_hdr%nkpt=1
  ddb_hdr%nspden=0
  ddb_hdr%nspinor=0
  ddb_hdr%nsppol=1
  ddb_hdr%occopt=0
  ddb_hdr%pawecutdg=0
  ddb_hdr%dfpt_sciss=0
  ddb_hdr%tolwfr=0
  ddb_hdr%tphysel=0
  ddb_hdr%tsmear=0
  ddb_hdr%usepaw=0

  ddb_hdr%natom=InVar%natom_unitcell
  ddb_hdr%nsym=Sym%nsym
  ddb_hdr%ntypat=InVar%ntypat

!FB  write(6,*) 'BEFORE MALLOC'
!FB TO REMOVE  call ddb_hdr_malloc(ddb_hdr)
  call ddb_hdr%malloc()

!FB  write(6,*) 'BEFORE TABS'
  ddb_hdr%acell=Lattice%acell_unitcell
  ddb_hdr%rprim=Lattice%rprimt
  ddb_hdr%amu=InVar%amu
  ddb_hdr%nband=1
  ddb_hdr%symafm=0
  ddb_hdr%symrel=Sym%ptsymrel
  ddb_hdr%typat=InVar%typat_unitcell
  ddb_hdr%kpt=0
  ddb_hdr%occ=0
  ddb_hdr%spinat=0
  ddb_hdr%tnons=Sym%tnons
  ddb_hdr%wtk=0
  ddb_hdr%xred=Sym%xred_zero
  ddb_hdr%zion=0
  ddb_hdr%znucl=0

  if (MPIdata%iam_master) then
!FB  write(6,*) 'BEFORE OPEN DDB'
!FB TO REMOVE    call ddb_hdr_open_write(ddb_hdr, filename, unddb)
    call ddb_hdr%open_write(filename, unddb)

!FB  write(6,*) 'BEFORE WRITE BLOK'
! Print each blok of the DDB
    do iq_ibz=1,Qbz%nqibz
!FB  do iq_ibz=1,Qbz%nqbz
      call ddb%write_block(iq_ibz,choice,mband,mpert,msize,Qbz%nqibz,unddb)
!FB    call ddb%write_block(iq_ibz,choice,mband,mpert,msize,Qbz%nqbz,unddb)
    end do
    close(unddb)
  end if  

!!!!!!!!! Test !!!!!!!!  
  nqshft=1
!FB  ngqpt=InVar%ngqpt1/2.d0
  ngqpt=InVar%ngqpt1
  allocate(qshft(3,nqshft))
  qshft(:,:)=0.5d0
  ABI_ALLOCATE(qdrp_cart,(3,3,3,InVar%natom_unitcell))
!FB  qshft(:,:)=0.0d0
  if (MPIdata%iam_master) then
    comm=xmpi_world
!FB  call ifc_init_fromFile(dielt,trim(filename_bis),Ifc_tmp,natom_uc,ngqpt,nqshft,qshft,Crystal_tmp,zeff,comm)
    call ifc_init_fromFile(dielt,trim(filename),Ifc_tmp,natom_uc,ngqpt,nqshft,qshft,Crystal_tmp,zeff,qdrp_cart,comm)
    unitfile=2
    call tdep_write_ifc(Crystal_tmp,Ifc_tmp,InVar,MPIdata,InVar%natom_unitcell,unitfile)
  end if  
  ABI_DEALLOCATE(qdrp_cart)

 end subroutine tdep_write_ddb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tdep_read_ifc(Ifc,InVar,natom_unitcell)

  implicit none

  integer,intent(in) :: natom_unitcell
  type(Input_Variables_type),intent(in) :: InVar
  type(ifc_type),intent(inout) :: Ifc

  integer :: iatcell,ifcout,ii,irpt,jatcell,jj,mu
  double precision :: atmfrctot1,atmfrctot2,atmfrctot3,dist
  double precision :: atmfrclr1,atmfrclr2,atmfrclr3
  double precision :: atmfrcsr1,atmfrcsr2,atmfrcsr3
  character (len=25):: string,string1,string2,string3

! Read IFC from ifc_in.dat or ifc_out.dat (according to ReadIFC)
! =======================================================
  ifcout=          min(Ifc%nrpt,200)
  Ifc%atmfrc(:,:,:,:,:)=zero
  if (InVar%ReadIFC==1) then
    write(InVar%stdout,'(a)') ' Read the IFC from ifc_in.dat'
    open(unit=40,file='ifc_in.dat')
  else if (InVar%ReadIFC==2) then
    write(InVar%stdout,'(a)') ' Read the IFC from ifc_out.dat'
    open(unit=40,file=trim(InVar%output_prefix)//'ifc_out.dat')
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
  if (InVar%loto) then
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
        if (InVar%loto) then
          read(40,*) atmfrctot1,atmfrctot2,atmfrctot3,atmfrclr1,atmfrclr2,atmfrclr3,atmfrcsr1,atmfrcsr2,atmfrcsr3
        else
          read(40,*) atmfrcsr1,atmfrcsr2,atmfrcsr3
        end if
        if (dist.lt.InVar%Rcut*0.99) then
!FB       See ifc_getiaf.F90 (for wghatm)
          Ifc%atmfrc(1,iatcell,mu,jatcell,irpt)=atmfrcsr1/Ifc%wghatm(iatcell,jatcell,irpt)
          Ifc%atmfrc(2,iatcell,mu,jatcell,irpt)=atmfrcsr2/Ifc%wghatm(iatcell,jatcell,irpt)
          Ifc%atmfrc(3,iatcell,mu,jatcell,irpt)=atmfrcsr3/Ifc%wghatm(iatcell,jatcell,irpt)
        else
          Ifc%atmfrc(1:3,iatcell,mu,jatcell,irpt)=zero
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
      if (InVar%loto) then
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
  write(InVar%stdout,'(a)') ' ------- achieved'
! When the IFC come from other calculations (ReadIFC=1), then impose acoustic sum rule
  call asrif9(Ifc%asr,Ifc%atmfrc,natom_unitcell,Ifc%nrpt,Ifc%rpt,Ifc%wghatm)

end subroutine tdep_read_ifc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tdep_write_ifc(Crystal,Ifc,InVar,MPIdata,natom_unitcell,unitfile)

  implicit none

  integer,intent(in) :: natom_unitcell,unitfile
  type(Input_Variables_type),intent(in) :: InVar
  type(ifc_type),intent(inout) :: Ifc
  type(crystal_t),intent(in) :: Crystal
  type(MPI_enreg_type), intent(in) :: MPIdata

  integer :: ifcana,ifcout,ncid,prt_ifc
  integer :: atifc(natom_unitcell)
  character(len=500) :: message

! Write IFC in ifc_out.dat or ifc_check.dat (according to ReadIFC)
! =======================================================
  ifcana=            1
  atifc(:)=          1
  ifcout=          min(Ifc%nrpt,200)
  prt_ifc=           1
  if (unitfile.eq.0) then
    write(InVar%stdout,'(a)') ' Write the IFC of TDEP in ifc_out.dat (and ifc_out.nc)'
    open(unit=77,file=trim(InVar%output_prefix)//'ifc_out.dat')
  else if (unitfile.eq.1) then
    write(InVar%stdout,'(a)') ' Write in ifc_check.dat (and ifc_check.nc) the IFC read previously'
    open(unit=77,file=trim(InVar%output_prefix)//'ifc_check.dat')
  else if (unitfile.eq.2) then
    write(InVar%stdout,'(a)') ' Write in ifc_ddb.dat (and ifc_ddb.nc) the IFC read from DDB file'
    open(unit=77,file=trim(InVar%output_prefix)//'ifc_ddb.dat')
  else
    write(message, '(a,i3,a)' )&
&   ' The value of unitfile ',unitfile,' is not allowed.'
    MSG_ERROR(message)
  end if
#ifdef HAVE_NETCDF
  if (unitfile.eq.0) then
    NCF_CHECK_MSG(nctk_open_create(ncid, trim(InVar%output_prefix)//"ifc_out.nc", xmpi_comm_self), "Creating ifc_out.nc")
  else if (unitfile.eq.1) then
    NCF_CHECK_MSG(nctk_open_create(ncid, trim(InVar%output_prefix)//"ifc_check.nc", xmpi_comm_self), "Creating ifc_check.nc")
  else if (unitfile.eq.2) then
    NCF_CHECK_MSG(nctk_open_create(ncid, trim(InVar%output_prefix)//"ifc_ddb.nc", xmpi_comm_self), "Creating ifc_ddb.nc")
  end if
  NCF_CHECK(nctk_def_basedims(ncid))
  NCF_CHECK(nctk_defnwrite_ivars(ncid, ["anaddb_version"], [1]))
  NCF_CHECK(crystal%ncwrite(ncid))
  call ifc%write(ifcana,atifc,ifcout,prt_ifc,ncid,77)
#endif
  close(77)
  write(InVar%stdout,'(a)') ' ------- achieved'

end subroutine tdep_write_ifc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tdep_ifc2phi2(dipdip,Ifc,InVar,Lattice,natom_unitcell,option,Phi2,Rlatt4abi,Shell2at,Sym)

  implicit none

  integer,intent(in) :: dipdip,natom_unitcell,option
  type(Input_Variables_type),intent(in) :: InVar
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Shell_Variables_type),intent(in) :: Shell2at
  type(ifc_type),intent(inout) :: Ifc
  double precision, intent(inout) :: Phi2(3*InVar%natom,3*InVar%natom)
  double precision, intent(in) :: Rlatt4abi(3,natom_unitcell,InVar%natom)

  integer :: eatom,fatom,iatcell,ii,irpt,jatcell,jatom,jj,isym,kk
  integer :: ishell,iatref,jatref,iatshell,trans
  double precision :: tol,dist
  double precision :: tmp(3)
  double precision :: vect_ddb(3),vect_tdep(3),xcart(3,natom_unitcell)
  double precision,allocatable :: Phi2_33(:,:),Phi2_ref(:,:,:)

! Link the atmfrc to the Phi2
! ==============================
  if (dipdip==0.or.InVar%ReadIFC.ne.0) then
    if (option==0) Ifc%atmfrc(:,:,:,:,:)=zero
    if (option==1) Phi2(:,:)=zero
  end if  
!FB  if (dipdip==0.and.option==0) Ifc%atmfrc(:,:,:,:,:)=zero
!FB  if (dipdip==0.and.option==1) Phi2(:,:)=zero
  call xred2xcart(natom_unitcell,Lattice%rprimt(:,:),xcart(:,:),Sym%xred_zero(:,:))
  tol=1.d-8

  do iatcell=1,natom_unitcell
    do jatcell=1,natom_unitcell
      do jatom=jatcell,InVar%natom,natom_unitcell
        vect_tdep(:)=Rlatt4abi(:,iatcell,jatom)+xcart(:,jatcell)-xcart(:,iatcell)
        tmp(:)=(vect_tdep(:)*Lattice%acell_unitcell(:))**2
        dist=dsqrt(sum(tmp(:)))
        if (dist.gt.(InVar%Rcut*0.99).and.option==1) then
          Phi2((iatcell-1)*3+1:(iatcell-1)*3+3,(jatom-1)*3+1:(jatom-1)*3+3)=zero
          cycle
        end if
        do irpt=1,Ifc%nrpt
          vect_ddb (:)=Ifc%rpt(:,irpt)+Ifc%rcan(:,jatcell)-Ifc%rcan(:,iatcell)
          tmp(:)=(vect_ddb(:)*Lattice%acell_unitcell(:))**2
          dist=dsqrt(sum(tmp(:)))
          if (dist.gt.(InVar%Rcut*0.99).and.option==0) then
            Ifc%atmfrc(:,iatcell,:,jatcell,irpt)=zero
            cycle
          end if
          if (abs(vect_ddb(1)-vect_tdep(1)).lt.tol.and.&
&             abs(vect_ddb(2)-vect_tdep(2)).lt.tol.and.&
&             abs(vect_ddb(3)-vect_tdep(3)).lt.tol) then
            do ii=1,3
              do jj=1,3
                if (option==0) then
                  Ifc%atmfrc(ii,iatcell,jj,jatcell,irpt)=Ifc%atmfrc(ii,iatcell,jj,jatcell,irpt)+&
&                         Phi2(ii+(iatcell-1)*3,(jatom-1)*3+jj)/Ifc%wghatm(iatcell,jatcell,irpt)
                else if (option==1) then
                  Phi2(ii+(iatcell-1)*3,(jatom-1)*3+jj)=Ifc%atmfrc(ii,iatcell,jj,jatcell,irpt)*Ifc%wghatm(iatcell,jatcell,irpt)
                  if (abs(Phi2(ii+(iatcell-1)*3,(jatom-1)*3+jj)-0.000040).lt.1.d-6) then
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

  if (option==1) then
!   Build the Phi2
    ABI_MALLOC(Phi2_ref,(3,3,Shell2at%nshell)); Phi2_ref(:,:,:)=zero
    ABI_MALLOC(Phi2_33,(3,3)) ; Phi2_33(:,:)=zero
    do ishell=1,Shell2at%nshell
!     Build the 3x3 IFC per shell
      iatref=Shell2at%iatref(ishell)
      jatref=Shell2at%jatref(ishell)
      Phi2_ref(:,:,ishell)=Phi2((iatref-1)*3+1:(iatref-1)*3+3,(jatref-1)*3+1:(jatref-1)*3+3)
    end do
    Phi2(:,:)=zero
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
          Phi2((eatom-1)*3+1:(eatom-1)*3+3,3*(fatom-1)+1:3*(fatom-1)+3)=Phi2_33(:,:)
          do ii=1,3
            do jj=1,3
              Phi2((fatom  -1)*3+ii,3*(eatom-1)+jj)=Phi2_33(jj,ii)
            end do
          end do
        end do !iatshell
      end do !eatom
    end do !ishell
!   Acoustic sum rule
    do eatom=1,InVar%natom
      Phi2((eatom-1)*3+1:(eatom-1)*3+3,(eatom-1)*3+1:(eatom-1)*3+3)=zero
      do jj=1,3
        do kk=1,3
          do fatom=1,InVar%natom
            if (fatom==eatom) cycle
            Phi2((eatom-1)*3+jj,(eatom-1)*3+kk)=Phi2((eatom-1)*3+jj,3*(eatom-1)+kk)&
&                                                 -Phi2((eatom-1)*3+jj,3*(fatom-1)+kk)
          enddo !fatom
        enddo !kk
      enddo !jj
    end do !eatom
    ABI_FREE(Phi2_33)
    ABI_FREE(Phi2_ref)
  end if

end subroutine tdep_ifc2phi2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module m_tdep_abitypes
