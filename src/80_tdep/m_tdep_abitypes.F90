
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

  use m_geometry,         only : xred2xcart
  use m_dynmat,           only : asrif9
  use m_tdep_readwrite,   only : Input_Variables_type
  use m_tdep_latt,        only : Lattice_Variables_type
  use m_tdep_phij,        only : tdep_build_phij33
  use m_tdep_sym,         only : Symetries_Variables_type
  use m_tdep_shell,       only : Shell_Variables_type
  use m_ifc,              only : ifc_type, ifc_init
  use m_crystal,          only : crystal_t, crystal_init
  use m_ddb,              only : ddb_type
  use m_kpts,             only : smpbz

  implicit none

  public :: tdep_init_crystal
  public :: tdep_init_ifc
  public :: tdep_read_ifc
  public :: tdep_write_ifc
  public :: tdep_ifc2phij
  public :: tdep_init_ddb

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

 subroutine tdep_init_ifc(Crystal,DDB,Ifc,InVar,Lattice,Phij_NN,Rlatt_cart,Shell2at,Sym)

  implicit none
  type(crystal_t),intent(in) :: Crystal
  type(ifc_type),intent(out) :: Ifc
  type(ddb_type),intent(in) :: DDB
  type(Input_Variables_type),intent(in) :: InVar
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Shell_Variables_type),intent(in) :: Shell2at
  double precision,intent(out) :: Phij_NN(3*InVar%natom,3*InVar%natom)
  double precision,intent(in) :: Rlatt_cart(3,InVar%natom_unitcell,InVar%natom)

  integer :: iatcell,dipdip,ii,asr,symdynmat,rfmeth
  integer :: iatom,nsphere,prtsrlr,enunit,nqshft
  integer :: ngqpt_in(3)
  double precision :: rifcsph
  double precision :: dielt(3,3)
  double precision, allocatable :: zeff(:,:,:)
  double precision, allocatable :: qdrp_cart(:,:,:,:)
  double precision, allocatable :: q1shft(:,:)

! Define matrices for LO-TO
! =========================
  ABI_MALLOC(zeff, (3,3,InVar%natom_unitcell)); zeff (:,:,:)=zero
  ABI_MALLOC(qdrp_cart, (3,3,3,InVar%natom_unitcell))
  dielt(:,:)=    zero
  dielt(1,1)=    1.d0
  dielt(2,2)=    1.d0
  dielt(3,3)=    1.d0
  qdrp_cart = zero
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
!  prtsrlr=     1
  prtsrlr=     0
  enunit=      0

  ngqpt_in(:)= InVar%ngqpt1(:)
  nqshft=      1
  ABI_MALLOC(q1shft,(3,nqshft)); q1shft(:,:)=0.0d0

  call ifc_init(Ifc,Crystal,DDB,Lattice%brav,asr,symdynmat,dipdip,&
  rfmeth,ngqpt_in,nqshft,q1shft,dielt,zeff,qdrp_cart,nsphere,rifcsph,&
  prtsrlr,enunit,XMPI_WORLD)

  ABI_FREE(q1shft)
  ABI_FREE(zeff)
  ABI_FREE(qdrp_cart)

! Read an IFC from ifc.in input file, write it in the ifc.out file and copy to Phij
! =================================================================================
  if (InVar%ReadIFC.eq.1) then
    write(InVar%stdout,*) ' '
    write(InVar%stdout,*) '#############################################################################'
    write(InVar%stdout,*) '################ Read IFCs from input file ##################################'
    write(InVar%stdout,*) '#############################################################################'
!   Read IFC from ifc.in (ReadIFC=1)
    call tdep_read_ifc(Ifc,InVar,InVar%natom_unitcell)
!   Copy Ifc%atmfrc to Phij_NN
    call tdep_ifc2phij(Ifc%dipdip,Ifc,InVar,Lattice,InVar%natom_unitcell,1,Phij_NN,Rlatt_cart,Shell2at,Sym)
!   Copy Phij_NN to Ifc%atmfrc
    call tdep_ifc2phij(Ifc%dipdip,Ifc,InVar,Lattice,InVar%natom_unitcell,0,Phij_NN,Rlatt_cart,Shell2at,Sym)
!   Write IFC in ifc.out (for check)
    call tdep_write_ifc(Crystal,Ifc,InVar,InVar%natom_unitcell,1)

!   Write the Phij_NN.dat file
    if (InVar%debug) then
      write(InVar%stdout,'(a)') ' See the Phij_NN.dat file corresponding to the ifc.in file'
      open(unit=55,file=trim(InVar%output_prefix)//'Phij_NN.dat')
      do iatom=1,3*InVar%natom
        write(55,'(10000(f10.6,1x))') Phij_NN(iatom,:)
      end do
      close(55)
    end if
  end if

 end subroutine tdep_init_ifc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine tdep_init_ddb(DDB,InVar,Lattice)

  use m_copy,             only : alloc_copy
  implicit none
  type(ddb_type),intent(out) :: DDB
  type(Input_Variables_type),intent(in) :: InVar
  type(Lattice_Variables_type),intent(in) :: Lattice

  integer :: option,nqshft,mqpt,nqbz
  integer :: mpert,nblok,msize
  integer :: qptrlatt(3,3),ngqpt_in(3)
  double precision, allocatable :: spqpt(:,:)
  double precision, allocatable :: q1shft(:,:)

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
  ABI_MALLOC(q1shft,(3,nqshft)); q1shft(:,:)=0.0d0
  call smpbz(Lattice%brav,InVar%stdlog,qptrlatt,mqpt,nqbz,nqshft,option,q1shft,spqpt)

! Initialize the DDB datatype (only needed values)
! ================================================
!! flg(nsize,nblok)= flag of existence for each element of the DDB:
!!                   integer,intent(in) :: blkflg(3,mpert,3,mpert,nblok)
!! nrm(1,nblok)=norm of qpt providing normalization
!!                   real(dp),intent(in) :: blknrm(3,nblok)
!! qpt(1<ii<9,nblok)=q vector of a phonon mode (ii=1,2,3)
!!                   real(dp),intent(in) :: blkqpt(9,nblok)
!! typ(nblok)=1 or 2 depending on non-stationary or stationary block 3 for third order derivatives
!!                   integer,intent(in) :: blktyp(nblok)
!! val(2,3*mpert*3*mpert,nblok)= all the dynamical matrices
!!                   real(dp),intent(in) :: blkval(2,3*mpert*3*mpert,nblok)
!! nblok=number of blocks in the DDB
  mpert=InVar%natom_unitcell+MPERT_MAX
  msize=3*mpert*3*mpert
  nblok=nqbz
  ABI_MALLOC(DDB%flg,(msize,nblok))  ; DDB%flg(:,:)  =1
  ABI_MALLOC(DDB%nrm,(3,nblok))      ; DDB%nrm(:,:)  =zero
  ABI_MALLOC(DDB%qpt,(9,nblok))      ; DDB%qpt(:,:)  =zero
  ABI_MALLOC(DDB%typ,(nblok))        ; DDB%typ(:)    =1
  ABI_MALLOC(DDB%val,(2,msize,nblok)); DDB%val(:,:,:)=zero
  DDB%nblok=nblok
  DDB%nrm(1,:)=1.d0
  DDB%qpt(1:3,:)=spqpt(:,:)
  call alloc_copy(InVar%amu,DDB%amu)
  DDB%rprim=Lattice%rprimt
  DDB%gprim=Lattice%gprim
  DDB%acell=Lattice%acell_unitcell

  ABI_FREE(spqpt)
  ABI_FREE(q1shft)

 end subroutine tdep_init_ddb
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

! Read IFC from ifc.in or ifc.tdep (according to ReadIFC)
! =======================================================
  ifcout=          min(Ifc%nrpt,200)
  Ifc%atmfrc(:,:,:,:,:)=zero
  if (InVar%ReadIFC==1) then
    write(InVar%stdout,'(a)') ' Read the IFC from ifc.in'
    open(unit=40,file='ifc.in')
  else if (InVar%ReadIFC==2) then
    write(InVar%stdout,'(a)') ' Read the IFC from ifc.tdep'
    open(unit=40,file='ifc.tdep')
  end if
  Ifc%atmfrc(:,:,:,:,:)=zero
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

subroutine tdep_write_ifc(Crystal,Ifc,InVar,natom_unitcell,unitfile)

  implicit none

  integer,intent(in) :: natom_unitcell,unitfile
  type(Input_Variables_type),intent(in) :: InVar
  type(ifc_type),intent(inout) :: Ifc
  type(crystal_t),intent(in) :: Crystal

  integer :: ifcana,ifcout,ncid,prt_ifc
  integer :: atifc(natom_unitcell)

! Write IFC in ifc.tdep or ifc.out (according to ReadIFC)
! =======================================================
  ifcana=            1
  atifc(:)=          1
  ifcout=          min(Ifc%nrpt,200)
  prt_ifc=           1
#ifdef HAVE_NETCDF
  if (unitfile.eq.0) then
    write(InVar%stdout,'(a)') ' Write the IFC of TDEP in ifc.tdep and anaddb.nc'
  else
    write(InVar%stdout,'(a)') ' Write in ifc.out and anaddb.nc the IFC read previously'
  end if
  NCF_CHECK_MSG(nctk_open_create(ncid, trim(InVar%output_prefix)//"anaddb.nc", xmpi_comm_self), "Creating anaddb.nc")
  NCF_CHECK(nctk_def_basedims(ncid))
  NCF_CHECK(nctk_defnwrite_ivars(ncid, ["anaddb_version"], [1]))
  NCF_CHECK(crystal%ncwrite(ncid))
!JB  call ifc%print(Ifc%dielt,Ifc%zeff,ifcana,atifc,ifcout,prt_ifc,ncid)
  call ifc%write(ifcana,atifc,ifcout,prt_ifc,ncid)
  write(InVar%stdout,'(a)') ' ------- achieved'
#else
  if (unitfile.eq.0) then
    open(unit=77,file=trim(InVar%output_prefix)//'ifc.tdep')
  else
    open(unit=77,file=trim(InVar%output_prefix)//'ifc.out')
  end if
  if (unitfile.eq.0) then
    write(InVar%stdout,'(a)') ' Write the IFC of TDEP in ifc.tdep'
  else
    write(InVar%stdout,'(a)') ' Write in ifc.out the IFC read previously'
  end if
  call ifc%print("TDEP",77,prt_ifc)
#endif
  close(77)

end subroutine tdep_write_ifc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tdep_ifc2phij(dipdip,Ifc,InVar,Lattice,natom_unitcell,option,Phij_NN,Rlatt4abi,Shell2at,Sym)

  implicit none

  integer,intent(in) :: dipdip,natom_unitcell,option
  type(Input_Variables_type),intent(in) :: InVar
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Shell_Variables_type),intent(in) :: Shell2at
  type(ifc_type),intent(inout) :: Ifc
  double precision, intent(inout) :: Phij_NN(3*InVar%natom,3*InVar%natom)
  double precision, intent(in) :: Rlatt4abi(3,natom_unitcell,InVar%natom)

  integer :: eatom,fatom,iatcell,ii,irpt,jatcell,jatom,jj,isym,kk
  integer :: ishell,iatref,jatref,iatshell,trans
  double precision :: tol,dist
  double precision :: tmp(3)
  double precision :: vect_ddb(3),vect_tdep(3),xcart(3,natom_unitcell)
  double precision,allocatable :: Phij_33(:,:),Phij_ref(:,:,:)

! Link the atmfrc to the Phij_NN
! ==============================
  if (dipdip==0.and.option==0) Ifc%atmfrc(:,:,:,:,:)=zero
  if (dipdip==0.and.option==1) Phij_NN(:,:)=zero
  call xred2xcart(natom_unitcell,Lattice%rprimt(:,:),xcart(:,:),Sym%xred_zero(:,:))
  tol=1.d-8

  do iatcell=1,natom_unitcell
    do jatcell=1,natom_unitcell
      do jatom=jatcell,InVar%natom,natom_unitcell
        vect_tdep(:)=Rlatt4abi(:,iatcell,jatom)+xcart(:,jatcell)-xcart(:,iatcell)
        tmp(:)=(vect_tdep(:)*Lattice%acell_unitcell(:))**2
        dist=dsqrt(sum(tmp(:)))
        if (dist.gt.(InVar%Rcut*0.99).and.option==1) then
          Phij_NN((iatcell-1)*3+1:(iatcell-1)*3+3,(jatom-1)*3+1:(jatom-1)*3+3)=zero
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
&                         Phij_NN(ii+(iatcell-1)*3,(jatom-1)*3+jj)/Ifc%wghatm(iatcell,jatcell,irpt)
                else if (option==1) then
                  Phij_NN(ii+(iatcell-1)*3,(jatom-1)*3+jj)=Ifc%atmfrc(ii,iatcell,jj,jatcell,irpt)*Ifc%wghatm(iatcell,jatcell,irpt)
                  if (abs(Phij_NN(ii+(iatcell-1)*3,(jatom-1)*3+jj)-0.000040).lt.1.d-6) then
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
!   Build the Phij_NN
    ABI_MALLOC(Phij_ref,(3,3,Shell2at%nshell)); Phij_ref(:,:,:)=zero
    ABI_MALLOC(Phij_33,(3,3)) ; Phij_33(:,:)=zero
    do ishell=1,Shell2at%nshell
!     Build the 3x3 IFC per shell
      iatref=Shell2at%iatref(ishell)
      jatref=Shell2at%jatref(ishell)
      Phij_ref(:,:,ishell)=Phij_NN((iatref-1)*3+1:(iatref-1)*3+3,(jatref-1)*3+1:(jatref-1)*3+3)
    end do
    Phij_NN(:,:)=zero
    do ishell=1,Shell2at%nshell
      do eatom=1,Invar%natom
!       Build the 3x3 IFC of an atom in this shell
        if (Shell2at%neighbours(eatom,ishell)%n_interactions.eq.0) cycle
        do iatshell=1,Shell2at%neighbours(eatom,ishell)%n_interactions
          fatom=Shell2at%neighbours(eatom,ishell)%atomj_in_shell(iatshell)
          isym =Shell2at%neighbours(eatom,ishell)%sym_in_shell(iatshell)
          trans=Shell2at%neighbours(eatom,ishell)%transpose_in_shell(iatshell)
          if (fatom.lt.eatom) cycle
          call tdep_build_phij33(isym,Phij_ref(:,:,ishell),Phij_33,Sym,trans)
!         Symetrization of the Phij_NN matrix
          Phij_NN((eatom-1)*3+1:(eatom-1)*3+3,3*(fatom-1)+1:3*(fatom-1)+3)=Phij_33(:,:)
          do ii=1,3
            do jj=1,3
              Phij_NN((fatom  -1)*3+ii,3*(eatom-1)+jj)=Phij_33(jj,ii)
            end do
          end do
        end do !iatshell
      end do !eatom
    end do !ishell
!   Acoustic sum rule
    do eatom=1,InVar%natom
      Phij_NN((eatom-1)*3+1:(eatom-1)*3+3,(eatom-1)*3+1:(eatom-1)*3+3)=zero
      do jj=1,3
        do kk=1,3
          do fatom=1,InVar%natom
            if (fatom==eatom) cycle
            Phij_NN((eatom-1)*3+jj,(eatom-1)*3+kk)=Phij_NN((eatom-1)*3+jj,3*(eatom-1)+kk)&
&                                                 -Phij_NN((eatom-1)*3+jj,3*(fatom-1)+kk)
          enddo !fatom
        enddo !kk
      enddo !jj
    end do !eatom
    ABI_FREE(Phij_33)
    ABI_FREE(Phij_ref)
  end if

end subroutine tdep_ifc2phij
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module m_tdep_abitypes
