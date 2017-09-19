#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdepdos

  use defs_basis
  use m_nctk
  use m_xmpi
  use m_errors
  use m_phonons
  use m_ifc
  use m_crystal
  use m_ddb,         only : ddb_type
  use m_qpt,         only : Qpoints_type
  use m_crystal_io,  only : crystal_ncwrite
  use m_dynmat,      only : make_bigbox, canat9, wght9, q0dy3_calc, asrif9, ftifc_q2r, ftifc_r2q
  use m_copy,        only : alloc_copy
  use m_readwrite,   only : Input_Variables_type
  use m_latt,        only : Lattice_Variables_type
  use m_sym,         only : make_sym, SearchMatR_1at, SearchMatR_2at, Symetries_Variables_type
#ifdef HAVE_NETCDF
  use netcdf
#endif

  implicit none

  public :: make_phdos
  public :: thermo
  public :: elastic

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine make_phdos(Phij_NN,Ifc,InVar,Lattice,natom,natom_unitcell,PHdos,Qpt,Rlatt_cart,Sym)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_phdos'
 use interfaces_41_geometry
 use interfaces_56_recipspace
!End of the abilint section

  implicit none 

  integer :: prtdos,option,nqpt,ii,jj,timrev,npsp,nqbz,iqpt
  integer :: mrpt,mqpt,msym,iatom,natom,natom_unitcell,irpt_new
  integer :: irpt,iatcell,jatcell,ifcana,ifcout,prt_ifc,ncid,mu,nu
  integer :: asr,symdynmat,dipdip,rfmeth,nqshft,nsphere,prtsrlr,enunit
  integer :: mpert,nblok,msize,iomega
  integer :: dos_ngqpt(3),qptrlatt(3,3),ngqpt_in(3)
  integer :: atifc(natom_unitcell)
  logical :: remove_inv,use_antiferro
  character (len=25):: string,string1,string2,string3,phdos_fname
  double precision :: dossmear,tol,rifcsph,dist,integ,domega
  double precision :: atmfrctot1,atmfrctot2,atmfrctot3
  double precision :: atmfrclr1,atmfrclr2,atmfrclr3
  double precision :: atmfrcsr1,atmfrcsr2,atmfrcsr3
  double precision :: Phij_NN(3*natom_unitcell,3*natom)
  double precision :: Rlatt_cart(3,natom_unitcell,natom)
  double precision :: Rlatt_tmp(3,natom_unitcell,natom)
  double precision :: dos_qshift(3),dielt(3,3)
  double precision :: vect_ddb(3),vect_tdep(3),xcart(3,natom_unitcell)
  double precision :: tmp(3)
  double precision, allocatable :: zion(:),znucl(:),zeff(:,:,:)
  double precision, allocatable :: spqpt(:,:)
  double precision, allocatable :: displ(:,:),omega(:,:)
  double precision, allocatable :: q1shft(:,:)
!FB  double precision, allocatable :: atmfrc_loto(:,:,:,:,:,:)
!FB  double precision, allocatable :: dynmat(:,:,:,:,:,:)
  character(len=132), allocatable :: title(:)
  type(Input_Variables_type),intent(in) :: InVar
  type(phonon_dos_type),intent(out) :: PHdos
  type(ifc_type),intent(out) :: Ifc
  type(ifc_type) :: Ifc_tmp
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Symetries_Variables_type),intent(in) :: Sym
  type(crystal_t) :: Crystal
  type(Qpoints_type) :: Qpt
  type(ddb_type) :: ddb

  write(InVar%stdout,*)' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '################### vibrational Density OF States (vDOS) ####################'
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,'(a)') ' See the vdos.dat and TDEP_PHDOS* files'

! Initialize the Crystal datatype
! ===============================
  timrev=1
  use_antiferro=.false.
  remove_inv=.false.
  npsp=InVar%ntypat
  ABI_MALLOC(zion,(InVar%ntypat)) ; zion(:)=zero
  ABI_MALLOC(znucl,(npsp)) ; znucl(:)=zero
  ABI_MALLOC(title,(InVar%ntypat))
  call crystal_init(InVar%amu,Crystal,Sym%spgroup,natom_unitcell,npsp,InVar%ntypat,Sym%nsym,Lattice%rprimdt,InVar%typat_unitcell,Sym%xred_zero,&
&  zion,znucl,timrev,use_antiferro,remove_inv,title,&
&  Sym%ptsymrel(:,:,1:Sym%nsym),Sym%tnons(:,1:Sym%nsym),Sym%symafm(1:Sym%nsym)) ! Optional

! Define matrices for LO-TO  
! =========================
  ABI_MALLOC(zeff, (3,3,natom_unitcell)); zeff (:,:,:)=zero
  dielt(:,:)=    zero
  dielt(1,1)=    1.d0
  dielt(2,2)=    1.d0
  dielt(3,3)=    1.d0
  if (InVar%loto) then
    dipdip=1
    do iatcell=1,natom_unitcell
      do ii=1,3
        zeff(ii,ii,iatcell)=InVar%born_charge(InVar%typat_unitcell(iatcell))
        dielt(ii,ii)       =InVar%dielec_constant
      end do
    end do  
  else  
    dipdip=0
  end if

! Compute the number of points in the Brilloiun zone: nqbz and qbz
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
  mpert=natom_unitcell+6
  msize=3*mpert*3*mpert
  nblok=nqbz
  ABI_MALLOC(ddb%flg,(msize,nblok))  ; ddb%flg(:,:)  =1
  ABI_MALLOC(ddb%nrm,(3,nblok))      ; ddb%nrm(:,:)  =zero
  ABI_MALLOC(ddb%qpt,(9,nblok))      ; ddb%qpt(:,:)  =zero
  ABI_MALLOC(ddb%typ,(nblok))        ; ddb%typ(:)    =1
  ABI_MALLOC(ddb%val,(2,msize,nblok)); ddb%val(:,:,:)=zero
  ddb%nblok=nblok
  ddb%nrm(1,:)=1.d0
  ddb%qpt(1:3,:)=spqpt(:,:)
  call alloc_copy(InVar%amu,ddb%amu)
  ddb%rprim=Lattice%rprimt
  ddb%gprim=Lattice%gprim
  ddb%acell=Lattice%acell_unitcell

! Initialize the Ifc datatype
! ===========================
  asr=         1
!FB  asr=         2
  symdynmat=   1
  rfmeth=      1
  nsphere=     0
  rifcsph=     0
!  prtsrlr=     1
  prtsrlr=     0
  enunit=      0
  call ifc_init(Ifc,Crystal,ddb,Lattice%brav,asr,symdynmat,dipdip,&
  rfmeth,ngqpt_in,nqshft,q1shft,dielt,zeff,nsphere,rifcsph,&
  prtsrlr,enunit,xmpi_world)

  do iatcell=1,natom_unitcell
    do iatom=1,natom
!FB      Rlatt_tmp(:,iatcell,iatom)=Rlatt_cart(:,iatcell,iatom)/Lattice%acell_unitcell(:)
      Rlatt_tmp(:,iatcell,iatom)=Rlatt_cart(:,iatcell,iatom)
    end do   
  end do  

! Link the atmfrc to the Phij_NN  
! ==============================
!FB
!FB  ABI_MALLOC(atmfrc_loto,(2,3,natom_unitcell,3,natom_unitcell,Ifc%nrpt)); atmfrc_loto(:,:,:,:,:,:)=zero 
!FB  atmfrc_loto(:,:,:,:,:,:)=Ifc%atmfrc(:,:,:,:,:,:)
!FB  Ifc%atmfrc(:,:,:,:,:,:)=zero
!FB! From atmfrc to dynmat  
!FB  ABI_MALLOC(dynmat,(2,3,natom_unitcell,3,natom_unitcell,nqbz)); dynmat(:,:,:,:,:,:)=zero
!FB  call ftifc_r2q(Ifc%atmfrc,dynmat,Lattice%gprim,natom_unitcell,nqbz,Ifc%nrpt,Ifc%rpt,spqpt,Ifc%wghatm)
!FB  do iqpt=1,nqbz
!FB    write(InVar%stdout,*) ddb%qpt(1:3,iqpt)
!FB!   Calculation of the k coordinates in Normalized Reciprocal coordinates
!FB    tmp(1)=(spqpt(1,iqpt)*Lattice%gprim(1,1)+spqpt(2,iqpt)*Lattice%gprim(1,2)+spqpt(3,iqpt)*Lattice%gprim(1,3))**2
!FB    tmp(2)=(spqpt(1,iqpt)*Lattice%gprim(2,1)+spqpt(2,iqpt)*Lattice%gprim(2,2)+spqpt(3,iqpt)*Lattice%gprim(2,3))**2
!FB    tmp(3)=(spqpt(1,iqpt)*Lattice%gprim(3,1)+spqpt(2,iqpt)*Lattice%gprim(3,2)+spqpt(3,iqpt)*Lattice%gprim(3,3))**2
!FB    dist=sum(tmp(:))
!FB    write(InVar%stdout,*) 'norm**2=',dist,dynmat(:,:,:,:,:,iqpt)
!FB    dynmat(:,:,:,:,:,iqpt)=dynmat(:,:,:,:,:,iqpt)*exp(-dist/1.**2)
!FB  end do
!FB! From dynmat to atmfrc
!FB  Ifc%atmfrc(:,:,:,:,:,:)=zero
!FB  call ftifc_q2r(Ifc%atmfrc,dynmat,Lattice%gprim,natom_unitcell,nqbz,Ifc%nrpt,Ifc%rpt,spqpt)
  ABI_FREE(spqpt)
!FB
  if (dipdip==0) Ifc%atmfrc(:,:,:,:,:,:)=zero
  call xred2xcart(natom_unitcell,Lattice%rprimt(:,:),xcart(:,:),Sym%xred_zero(:,:))
  tol=1.d-8
  do iatcell=1,natom_unitcell
    do jatcell=1,natom_unitcell
      do iatom=jatcell,natom,natom_unitcell
        vect_tdep(:)=Rlatt_tmp(:,iatcell,iatom)+xcart(:,jatcell)-xcart(:,iatcell)
        do irpt=1,Ifc%nrpt
          vect_ddb (:)=Ifc%rpt(:,irpt)+Ifc%rcan(:,jatcell)-Ifc%rcan(:,iatcell)
          tmp(:)=(vect_ddb(:)*Lattice%acell_unitcell(:))**2
          dist=dsqrt(sum(tmp(:)))
          if (dist.gt.InVar%Rcut.and.dipdip==1) then  
            Ifc%atmfrc(1,:,iatcell,:,jatcell,irpt)=zero
            cycle
          end if
          if (abs(vect_ddb(1)-vect_tdep(1)).lt.tol.and.&
&             abs(vect_ddb(2)-vect_tdep(2)).lt.tol.and.&
&             abs(vect_ddb(3)-vect_tdep(3)).lt.tol) then
            do ii=1,3
              do jj=1,3
                Ifc%atmfrc(1,ii,iatcell,jj,jatcell,irpt)=Ifc%atmfrc(1,ii,iatcell,jj,jatcell,irpt)+Phij_NN(ii+(iatcell-1)*3,(iatom-1)*3+jj)
              enddo
            enddo
          endif
        end do
      end do
    end do
  end do
!FB  if (InVar%loto) then
    call asrif9(Ifc%asr,Ifc%atmfrc,natom_unitcell,Ifc%nrpt,Ifc%rpt,Ifc%wghatm)
!FB  end if  

! Write IFC in ifc.tdep
! =====================
  ifcana=            1
  atifc(:)=          1
  ifcout=          min(Ifc%nrpt,200)
  prt_ifc=           1
  open(unit=7,file='ifc.tdep')
#ifdef HAVE_NETCDF
  write(InVar%stdout,'(a)') 'Write the IFC of TDEP in ifc.tdep and anaddb.nc' 
  NCF_CHECK_MSG(nctk_open_create(ncid, "anaddb.nc", xmpi_comm_self), "Creating anaddb.nc")
  NCF_CHECK(nctk_def_basedims(ncid))
  NCF_CHECK(nctk_defnwrite_ivars(ncid, ["anaddb_version"], [1]))
  NCF_CHECK(crystal_ncwrite(crystal,ncid))
  call ifc_write(Ifc,ifcana,atifc,ifcout,prt_ifc,ncid)
  write(InVar%stdout,'(a)') '------- achieved'
#else  
  write(InVar%stdout,'(a)') 'Write the IFC of TDEP in ifc.tdep' 
  call ifc_print(Ifc,unit=7)
#endif
  close(7)

! Read IFC from ifc.in (ifc.tdep)
! ===============================
  if (InVar%ReadIFC==1) then
    write(InVar%stdout,'(a)') 'Read the IFC from ifc.in' 
    open(unit=40,file='ifc.in')
  else if (InVar%ReadIFC==2) then  
    write(InVar%stdout,'(a)') 'Read the IFC from ifc.tdep' 
    open(unit=40,file='ifc.tdep')
  end if
  if ((InVar%ReadIFC.eq.1).or.(InVar%ReadIFC.eq.2)) then
    Ifc%atmfrc(:,:,:,:,:,:)=zero
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
        read(40,*) string
        do mu=1,3
          if (InVar%loto) then
            read(40,*) atmfrctot1,atmfrctot2,atmfrctot3,atmfrclr1,atmfrclr2,atmfrclr3,atmfrcsr1,atmfrcsr2,atmfrcsr3
          else
            read(40,*) atmfrcsr1,atmfrcsr2,atmfrcsr3
          end if  
!FB       See ifc_getiaf.F90 (for wghatm)          
          Ifc%atmfrc(1,1,iatcell,mu,jatcell,irpt)=atmfrcsr1/Ifc%wghatm(iatcell,jatcell,irpt)
          Ifc%atmfrc(1,2,iatcell,mu,jatcell,irpt)=atmfrcsr2/Ifc%wghatm(iatcell,jatcell,irpt)
          Ifc%atmfrc(1,3,iatcell,mu,jatcell,irpt)=atmfrcsr3/Ifc%wghatm(iatcell,jatcell,irpt)
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
!FB    if (InVar%loto) then
      call asrif9(Ifc%asr,Ifc%atmfrc,natom_unitcell,Ifc%nrpt,Ifc%rpt,Ifc%wghatm)
!FB    end if  

! Write IFC in ifc.out (for check)
! ================================
    open(unit=7,file='ifc.out')
#ifdef HAVE_NETCDF
    write(InVar%stdout,'(a)') 'Write the IFC of ANADDB in ifc.out and anaddb.nc' 
    NCF_CHECK_MSG(nctk_open_create(ncid, "anaddb.nc", xmpi_comm_self), "Creating anaddb.nc")
    NCF_CHECK(nctk_def_basedims(ncid))
    NCF_CHECK(nctk_defnwrite_ivars(ncid, ["anaddb_version"], [1]))
    NCF_CHECK(crystal_ncwrite(crystal,ncid))
    call ifc_write(Ifc,ifcana,atifc,ifcout,prt_ifc,ncid)
    write(InVar%stdout,'(a)') '------- achieved'
#else  
    write(InVar%stdout,'(a)') 'Write the IFC of ANADDB in ifc.out' 
    call ifc_print(Ifc,unit=7)
#endif
    close(7)
  end if !ReadIFC

! Compute the DOS
! ===============
  prtdos=1 !Gaussian
!  prtdos=2 !Tetra
  dossmear=4.5d-6
!  dossmear=4.5d-5
  dos_qshift(:)=     zero
  dos_ngqpt(:)=InVar%ngqpt2(:)
  write(InVar%stdout,'(a)') 'Compute the vDOS'
  call mkphdos(PHdos,Crystal,Ifc,prtdos,InVar%dosdeltae,dossmear,dos_ngqpt,dos_qshift,xmpi_world)
  write(InVar%stdout,'(a)') '------- achieved'

! Compute the frequencies
! =======================
  ABI_MALLOC(displ,(2*3*natom_unitcell*3*natom_unitcell,Qpt%nqpt)); displ(:,:)=zero
  ABI_MALLOC(omega,(3*natom_unitcell,Qpt%nqpt)); omega(:,:)=zero
  open(unit=53,file='omega-abinit.dat')
  write(6,*) 'My_qpt='
  do iqpt=1,Qpt%nqpt 
    call ifc_fourq(Ifc,Crystal,Qpt%qpt_red(:,iqpt),omega(:,iqpt),displ(:,iqpt))
    if (iqpt.le.Qpt%nqpt) write(53,'(i5,x,100(f15.6,x))') iqpt,(omega(ii,iqpt)*Ha_cmm1,ii=1,3*natom_unitcell)
  end do  
  close(53)

! Print the DOS
! =============
  phdos_fname = "TDEP_PHDOS"
  call phdos_print(PHdos,phdos_fname)
  domega=(InVar%dosdeltae*Ha_meV)
  integ=0.d0
  do iomega=1,PHdos%nomega
    integ=integ + domega*PHdos%phdos(iomega)
  end do
  PHdos%phdos(:)=PHdos%phdos(:)/integ
  open(unit=56,file='vdos.dat')
  do iomega=1,PHdos%nomega
    write(56,'(2(f18.6,x))') PHdos%omega(iomega)*Ha_eV*8065.544,PHdos%phdos(iomega)
  end do
  close(56)

end subroutine  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine thermo(DeltaFree_AH2,InVar,PHdos,U0)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'thermo'
!End of the abilint section

  implicit none 

  integer :: iomega,itemp
  double precision :: k_B,wovert,integ,expm2x,ln2shx,cothx,xx
  double precision :: Ftot,Fvib,domega
  double precision, intent(in) :: U0,DeltaFree_AH2
  type(Input_Variables_type),intent(in) :: InVar
  type(phonon_dos_type),intent(in) :: PHdos

  write(InVar%stdout,*)' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '################# Thermodynamic quantities: Free energy,...##################'
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,'(a)') ' See the thermo.dat file'

! The heat capacity
! =================
  open(unit=20,file='thermo.dat')
  k_B=8.617343d-5 !in eV/K
  domega=(InVar%dosdeltae*Ha_meV)
  wovert=(domega*1.d-3)/(2*InVar%temperature*k_B)
  integ=0.d0
  do iomega=1,PHdos%nomega
    xx=PHdos%omega(iomega)*Ha_eV*1.d+03/domega
    if (xx.lt.0) cycle
    integ=integ + (wovert*xx/sinh(wovert*xx))**2*PHdos%phdos(iomega)*domega
  end do
  integ=integ*3
  write(20,'(a)')'============= Direct results (without any inter/extrapolation) =================='
  write(20,'(x,a,f10.5)')'For present temperature (in Kelvin): T= ',InVar%temperature
  write(20,'(x,a,f10.5)')'  The specific heat (in k_b/atom): C_v=',integ

! The entropy
! ===========
  integ=0.d0
  do iomega=1,PHdos%nomega
    xx=PHdos%omega(iomega)*Ha_eV*1.d+03/domega
    if (xx.lt.0) cycle
    expm2x=exp(-2.d0*wovert*xx)
    ln2shx=wovert*xx+log(1.d0-expm2x)
    cothx=(1.d0+expm2x)/(1.d0-expm2x)
    integ=integ + ((wovert*xx)*cothx-ln2shx)*PHdos%phdos(iomega)*domega
  end do
  integ=integ*3
  write(20,'(x,a,f10.5)')'  The vibrational entropy (in k_b/atom): S_vib =',integ

! The free energy (for present temperature)
! =========================================
  integ=0.d0
  do iomega=1,PHdos%nomega
    xx=PHdos%omega(iomega)*Ha_eV*1.d+03/domega
    if (xx.lt.0) cycle
    integ=integ + log(2*sinh(wovert*xx))*PHdos%phdos(iomega)*domega
  end do
  Fvib=integ*3*k_B*InVar%temperature
  write(20,'(x,a,f12.5)')'  The cold contribution (in eV/atom): U_0 =',U0*Ha_eV
  write(20,'(x,a,f10.5)')'  The vibrational contribution (in eV/atom): F_vib = -T.S_vib =',Fvib
  write(20,'(x,a,f10.5)')'  The anharmonic contribution (in eV/atom): DeltaF_AH =',DeltaFree_AH2*Ha_eV
  Ftot=U0*Ha_eV+Fvib+DeltaFree_AH2*Ha_eV
  write(20,'(x,a)')'  So the free energy (in eV/atom) is equal to:'
  write(20,'(x,a,f12.5)')'     Harmonic only -->  F_tot^HA = U_0 + F_vib =',U0*Ha_eV+Fvib
  write(20,'(x,a,f12.5)')'     With anharmonic contribution -->  F_tot^AH = U_0 + F_vib + DeltaF_AH =',Ftot
  write(20,'(a)')' '

! The free energy (extrapolation)
! ===============================
  write(20,'(a)')'============= Quasi-Harmonic Approximation (QHA) =================='
  write(20,'(x,a)')'  Note that the following results come from an EXTRAPOLATION:'
  write(20,'(x,a,i5,a)')'    1/ F_vib^QHA(T) is computed for each T using vDOS(T=',int(InVar%temperature),')'
  write(20,'(x,a)')'    2/ F_tot^QHA(T) = F_vib^QHA(T) + U_0'
  write(20,'(x,a)')'    3/ We assume that DeltaF_AH^QHA(T)=a(V)*T**2'
  write(20,'(x,a)')'    4/ F_tot^QHA+AH(T) = U_0 + F_vib^QHA(T) + DeltaF_AH^QHA(T)'
  write(20,'(a)')'   T      F_vib^QHA(T)   F_tot^QHA(T)    DeltaF_AH^QHA(T)   F_tot^QHA+AH(T)'   
  do itemp=1,100
    wovert=(domega*1.d-3)/(2*dfloat(itemp)*100*k_B)
    integ=0.d0
    do iomega=1,PHdos%nomega
      xx=PHdos%omega(iomega)*Ha_eV*1.d+03/domega
      if (xx.lt.0) cycle
      integ=integ + log(2*sinh(wovert*xx))*PHdos%phdos(iomega)*domega
    end do
    integ=integ*3*k_B
    Ftot=U0*Ha_eV+integ*itemp*100+DeltaFree_AH2*Ha_eV*(itemp*100)**2/(InVar%temperature*100)**2
    write(20,'(x,i5,4(x,f15.5))') itemp*100,integ*itemp*100,U0*Ha_eV+integ*itemp*100,&
&                          DeltaFree_AH2*Ha_eV*(itemp*100)**2/(InVar%temperature)**2,Ftot
  end do  

  close(20)
end subroutine  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine elastic(Phij_NN,distance,InVar,Lattice)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elastic'
!End of the abilint section

  implicit none 

  integer :: iatom,ii,jj,kk,ll,iatcell,itypat
  double precision :: BH,BR,BV,GR,GV,GH,Eaverage,Nuaverage,Laverage,Vp,Vs,Vphi
  double precision :: rho,E1,E2,E3,Nu12,Nu13,Nu23,Nu21,Nu31,Nu32,G23,G13,G12
  double precision :: mass_amu,bohr
  double precision, allocatable :: Sij(:,:),Cij(:,:),aijkl(:,:,:,:),cijkl(:,:,:,:)
  type(Input_Variables_type), intent(in) :: InVar
  type(Lattice_Variables_type), intent(in) :: Lattice
  double precision, intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision, intent(in) :: Phij_NN(3*InVar%natom,3*InVar%natom)

  write(InVar%stdout,*)' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '######################### Elastic constants #################################'
  write(InVar%stdout,*) '################ Bulk and Shear modulus--Sound velocities ###################'
  write(InVar%stdout,*) '#############################################################################'

  bohr=0.5291772108e-10
! Define atomic mass average
  mass_amu=zero
  do iatom=1,InVar%natom_unitcell
    itypat=InVar%typat_unitcell(iatom) 
    mass_amu=mass_amu+InVar%amu(itypat)
  end do  
  mass_amu=mass_amu/dfloat(InVar%natom_unitcell)

  rho=(mass_amu/1e3)*InVar%natom_unitcell/Lattice%ucvol/bohr**3/6.022e23

!==========================================================================================
!===================== Elastic constants ==================================================
!==========================================================================================
! New calculation of elastic constants using the formula (12.28 and 12.29 of
! Wallace, Statistical physics of crystals and liquids, Worl Scientific)  
  ABI_MALLOC(aijkl,(3,3,3,3)); aijkl(:,:,:,:)=0.d0
  ABI_MALLOC(cijkl,(3,3,3,3)); cijkl(:,:,:,:)=0.d0
  do ii=1,3
    do jj=1,3
      do kk=1,3
        do ll=1,3
          do iatom=1,InVar%natom
            do iatcell=1,InVar%natom_unitcell
              aijkl(ii,jj,kk,ll)=aijkl(ii,jj,kk,ll)-Phij_NN(ii+(iatcell-1)*3,3*(iatom-1)+jj)*distance(iatcell,iatom,kk+1)*distance(iatcell,iatom,ll+1)/2.d0/Lattice%ucvol
            end do  
          end do
        enddo
      enddo
    enddo
  enddo

  do ii=1,3
    do jj=1,3
      do kk=1,3
        do ll=1,3
          cijkl(ii,jj,kk,ll)=aijkl(ii,kk,jj,ll)+aijkl(jj,kk,ii,ll)-aijkl(ii,jj,kk,ll)
        enddo
      enddo
    enddo
  enddo
  ABI_FREE(aijkl)
  
  cijkl(:,:,:,:)=cijkl(:,:,:,:)*29421.033d0

  ABI_MALLOC(Cij,(6,6)) ; Cij(:,:)=0.d0
  Cij(1,1)=cijkl(1,1,1,1) ; Cij(1,2)=cijkl(1,1,2,2) ; Cij(1,3)=cijkl(1,1,3,3) ; Cij(1,4)=cijkl(1,1,2,3) ; Cij(1,5)=cijkl(1,1,1,3) ; Cij(1,6)=cijkl(1,1,1,2)
  Cij(2,1)=cijkl(2,2,1,1) ; Cij(2,2)=cijkl(2,2,2,2) ; Cij(2,3)=cijkl(2,2,3,3) ; Cij(2,4)=cijkl(2,2,2,3) ; Cij(2,5)=cijkl(2,2,1,3) ; Cij(2,6)=cijkl(2,2,1,2)
  Cij(3,1)=cijkl(3,3,1,1) ; Cij(3,2)=cijkl(3,3,2,2) ; Cij(3,3)=cijkl(3,3,3,3) ; Cij(3,4)=cijkl(3,3,2,3) ; Cij(3,5)=cijkl(3,3,1,3) ; Cij(3,6)=cijkl(3,3,1,2)
  Cij(4,1)=cijkl(2,3,1,1) ; Cij(4,2)=cijkl(2,3,2,2) ; Cij(4,3)=cijkl(2,3,3,3) ; Cij(4,4)=cijkl(2,3,2,3) ; Cij(4,5)=cijkl(2,3,1,3) ; Cij(4,6)=cijkl(2,3,1,2)
  Cij(5,1)=cijkl(1,3,1,1) ; Cij(5,2)=cijkl(1,3,2,2) ; Cij(5,3)=cijkl(1,3,3,3) ; Cij(5,4)=cijkl(1,3,2,3) ; Cij(5,5)=cijkl(1,3,1,3) ; Cij(5,6)=cijkl(1,3,1,2)
  Cij(6,1)=cijkl(1,2,1,1) ; Cij(6,2)=cijkl(1,2,2,2) ; Cij(6,3)=cijkl(1,2,3,3) ; Cij(6,4)=cijkl(1,2,2,3) ; Cij(6,5)=cijkl(1,2,1,3) ; Cij(6,6)=cijkl(1,2,1,2)
! Remove the rounding errors before writing (for non regression testing purposes)
  do ii=1,6
    do jj=1,6
      if (abs(Cij(ii,jj)).lt.tol8) Cij(ii,jj)=zero
    end do
  end do  
  write(InVar%stdout,'(a)') ' '
  write(InVar%stdout,'(a)') '========== Using the formulation proposed by Wallace (using the IFC) ========='
  write(InVar%stdout,'(a)') 'Cijkl='
  write(InVar%stdout,'(a,6(f8.3,x))') '| C11 C12 C13 C14 C15 C16 |   ',Cij(1,1),Cij(1,2),Cij(1,3),Cij(1,4),Cij(1,5),Cij(1,6)
  write(InVar%stdout,'(a,6(f8.3,x))') '| C21 C22 C23 C24 C25 C26 |   ',Cij(2,1),Cij(2,2),Cij(2,3),Cij(2,4),Cij(2,5),Cij(2,6)
  write(InVar%stdout,'(a,6(f8.3,x))') '| C31 C32 C33 C34 C35 C36 |   ',Cij(3,1),Cij(3,2),Cij(3,3),Cij(3,4),Cij(3,5),Cij(3,6)
  write(InVar%stdout,'(a,6(f8.3,x))') '| C41 C42 C43 C44 C45 C46 | = ',Cij(4,1),Cij(4,2),Cij(4,3),Cij(4,4),Cij(4,5),Cij(4,6)
  write(InVar%stdout,'(a,6(f8.3,x))') '| C51 C52 C53 C54 C55 C56 |   ',Cij(5,1),Cij(5,2),Cij(5,3),Cij(5,4),Cij(5,5),Cij(5,6)
  write(InVar%stdout,'(a,6(f8.3,x))') '| C61 C62 C63 C64 C65 C66 |   ',Cij(6,1),Cij(6,2),Cij(6,3),Cij(6,4),Cij(6,5),Cij(6,6)

! Mean value of the off-diagonal elements
  Cij(1,2)=(Cij(1,2)+Cij(2,1))/2.d0 ; Cij(2,1)=Cij(1,2)
  Cij(1,3)=(Cij(1,3)+Cij(3,1))/2.d0 ; Cij(3,1)=Cij(1,3)
  Cij(2,3)=(Cij(2,3)+Cij(3,2))/2.d0 ; Cij(3,2)=Cij(2,3)
  
! Young's modulus
  E1=(Cij(1,1)*Cij(2,2)*Cij(3,3)+2.d0*Cij(2,3)*Cij(1,2)*Cij(1,3)-Cij(1,1)*Cij(2,3)**2-Cij(2,2)*Cij(1,3)**2-Cij(3,3)*Cij(1,2)**2)/(Cij(2,2)*Cij(3,3)-Cij(2,3)**2)
  E2=(Cij(1,1)*Cij(2,2)*Cij(3,3)+2.d0*Cij(2,3)*Cij(1,2)*Cij(1,3)-Cij(1,1)*Cij(2,3)**2-Cij(2,2)*Cij(1,3)**2-Cij(3,3)*Cij(1,2)**2)/(Cij(1,1)*Cij(3,3)-Cij(1,3)**2)
  E3=(Cij(1,1)*Cij(2,2)*Cij(3,3)+2.d0*Cij(2,3)*Cij(1,2)*Cij(1,3)-Cij(1,1)*Cij(2,3)**2-Cij(2,2)*Cij(1,3)**2-Cij(3,3)*Cij(1,2)**2)/(Cij(1,1)*Cij(2,2)-Cij(1,2)**2)
  write(InVar%stdout,'(a,3(f8.3,x))') 'Young modulus E1, E2 and E3=',E1,E2,E3

! Poisson Ratio
  Nu21=(Cij(1,2)*Cij(3,3)-Cij(1,3)*Cij(2,3))/(Cij(1,1)*Cij(3,3)-Cij(1,3)**2)
  Nu31=(Cij(1,3)*Cij(2,2)-Cij(1,2)*Cij(2,3))/(Cij(1,1)*Cij(2,2)-Cij(1,2)**2)
  Nu23=(Cij(1,1)*Cij(2,3)-Cij(1,2)*Cij(1,3))/(Cij(1,1)*Cij(3,3)-Cij(1,3)**2)
  Nu12=(Cij(1,2)*Cij(3,3)-Cij(1,3)*Cij(2,3))/(Cij(2,2)*Cij(3,3)-Cij(2,3)**2)
  Nu13=(Cij(2,2)*Cij(1,3)-Cij(1,2)*Cij(2,3))/(Cij(2,2)*Cij(3,3)-Cij(2,3)**2)
  Nu32=(Cij(1,1)*Cij(2,3)-Cij(1,2)*Cij(1,3))/(Cij(1,1)*Cij(2,2)-Cij(1,2)**2)
  write(InVar%stdout,'(a,6(f8.3,x))') 'Poisson ratio Nu21, Nu31, Nu23, Nu12, Nu13 and Nu32=',Nu21,Nu31,Nu23,Nu12,Nu13,Nu32
  
! Shear modulus  
  G23=Cij(4,4) ; G13=Cij(5,5) ; G12=Cij(6,6)
  write(InVar%stdout,'(a,3(f8.3,x))') 'Shear modulus G23, G13 and G12=',G23,G13,G12
  
! Compliance matrix  
  ABI_MALLOC(Sij,(6,6)) ; Sij(:,:)=0.d0
  Sij(1,1)= 1.d0/E1 ; Sij(1,2)=-Nu21/E2 ; Sij(1,3)=-Nu31/E3 ; Sij(1,4)=0.d0     ; Sij(1,5)=0.d0     ; Sij(1,6)=0.d0 
  Sij(2,1)=-Nu12/E1 ; Sij(2,2)= 1.d0/E2 ; Sij(2,3)=-Nu32/E3 ; Sij(2,4)=0.d0     ; Sij(2,5)=0.d0     ; Sij(2,6)=0.d0 
  Sij(3,1)=-Nu13/E1 ; Sij(3,2)=-Nu23/E2 ; Sij(3,3)= 1.d0/E3 ; Sij(3,4)=0.d0     ; Sij(3,5)=0.d0     ; Sij(3,6)=0.d0 
  Sij(4,1)= 0.d0    ; Sij(4,2)= 0.d0    ; Sij(4,3)= 0.d0    ; Sij(4,4)=1.d0/G23 ; Sij(4,5)=0.d0     ; Sij(4,6)=0.d0 
  Sij(5,1)= 0.d0    ; Sij(5,2)= 0.d0    ; Sij(5,3)= 0.d0    ; Sij(5,4)=0.d0     ; Sij(5,5)=1.d0/G13 ; Sij(5,6)=0.d0 
  Sij(6,1)= 0.d0    ; Sij(6,2)= 0.d0    ; Sij(6,3)= 0.d0    ; Sij(6,4)=0.d0     ; Sij(6,5)=0.d0     ; Sij(6,6)=1.d0/G12 
! Remove the rounding errors before writing (for non regression testing purposes)
  do ii=1,6
    do jj=1,6
      if (abs(Sij(ii,jj)).lt.tol8) Sij(ii,jj)=zero
    end do
  end do  
  write(InVar%stdout,'(a)') ' '
  write(InVar%stdout,'(a)') 'Sijkl='
  write(InVar%stdout,'(a,6(f8.3,x))') '| S11 S12 S13 S14 S15 S16 |   ',Sij(1,1),Sij(1,2),Sij(1,3),Sij(1,4),Sij(1,5),Sij(1,6)
  write(InVar%stdout,'(a,6(f8.3,x))') '| S21 S22 S23 S24 S25 S26 |   ',Sij(2,2),Sij(2,2),Sij(2,3),Sij(2,4),Sij(2,5),Sij(2,6)
  write(InVar%stdout,'(a,6(f8.3,x))') '| S31 S32 S33 S34 S35 S36 |   ',Sij(3,1),Sij(3,2),Sij(3,3),Sij(3,4),Sij(3,5),Sij(3,6)
  write(InVar%stdout,'(a,6(f8.3,x))') '| S41 S42 S43 S44 S45 S46 | = ',Sij(4,1),Sij(4,2),Sij(4,3),Sij(4,4),Sij(4,5),Sij(4,6)
  write(InVar%stdout,'(a,6(f8.3,x))') '| S51 S52 S53 S54 S55 S56 |   ',Sij(5,1),Sij(5,2),Sij(5,3),Sij(5,4),Sij(5,5),Sij(5,6)
  write(InVar%stdout,'(a,6(f8.3,x))') '| S61 S62 S63 S64 S65 S66 |   ',Sij(6,1),Sij(6,2),Sij(6,3),Sij(6,4),Sij(6,5),Sij(6,6)

!==========================================================================================
!===================== Bulk and Shear modulus--Sound velocities ===========================
!==========================================================================================
! Voigt notation  
  write(InVar%stdout,'(a,f9.3)')'For density rho=',rho
  write(InVar%stdout,*)' '
  write(InVar%stdout,*) '========================= Voigt average (constant strain) ===================' 
  BV=((Cij(1,1)+Cij(2,2)+Cij(3,3))+2.d0*(Cij(1,2)+Cij(1,3)+Cij(2,3)))/9.d0
  GV=((Cij(1,1)+Cij(2,2)+Cij(3,3))-     (Cij(1,2)+Cij(1,3)+Cij(2,3))+3.d0*(Cij(4,4)+Cij(5,5)+Cij(6,6)))/15.d0
  write(InVar%stdout,'(2(a,f9.3))')'ISOTHERMAL modulus: Bulk Kt=',BV,' and Shear G=',GV
  Eaverage=9.d0*BV*GV/(3*BV+GV)
  Nuaverage=0.5*(1.d0-(3.d0*GV)/(3.d0*BV+GV) )
  Laverage=(3.d0*BV-2.d0*GV)/3.d0
  write(InVar%stdout,'(3(a,f9.3))')'Average of Young modulus E=',Eaverage,' Lame modulus Lambda=',Laverage,' and Poisson ratio Nu=',Nuaverage
  Vp=dsqrt(1.d9*(BV+4.d0*GV/3.d0)/rho)
  Vs=dsqrt(1.d9*GV/rho)
  Vphi=dsqrt(1.d9*BV/rho)
  write(InVar%stdout,'(3(a,f9.3,x))')'Velocities: compressional Vp=',Vp,' shear Vs=',Vs,' and bulk Vphi=',Vphi

! Reuss notation  
  write(InVar%stdout,*)' '
  write(InVar%stdout,*) '========================= Reuss average (constant stress) ===================' 
  BR=1.d0/(Sij(1,1)+Sij(2,2)+Sij(3,3)+2.d0*(Sij(1,2)+Sij(1,3)+Sij(2,3)))
  GR=15.d0/(4.d0*(Sij(1,1)+Sij(2,2)+Sij(3,3))-4.d0*(Sij(1,2)+Sij(1,3)+Sij(2,3))+3.d0*(Sij(4,4)+Sij(5,5)+Sij(6,6)))
  write(InVar%stdout,'(2(a,f9.3))')'ISOTHERMAL modulus: Bulk Kt=',BR,' and Shear G=',GR
  Eaverage=9.d0*BR*GR/(3*BR+GR)
  Nuaverage=0.5*(1.d0-(3.d0*GR)/(3.d0*BR+GR) )
  Laverage=(3.d0*BR-2.d0*GR)/3.d0
  write(InVar%stdout,'(3(a,f9.3))')'Average of Young modulus E=',Eaverage,' Lame modulus Lambda=',Laverage,' and Poisson ratio Nu=',Nuaverage
  Vp=dsqrt(1.d9*(BR+4.d0*GR/3.d0)/rho)
  Vs=dsqrt(1.d9*GR/rho)
  Vphi=dsqrt(1.d9*BR/rho)
  write(InVar%stdout,'(3(a,f9.3,x))')'Velocities: compressional Vp=',Vp,' shear Vs=',Vs,' and bulk Vphi=',Vphi

! Voigt-Reuss-Hill notation
  write(InVar%stdout,*)' '
  write(InVar%stdout,*) '============================== Hill average =================================' 
  BH=(BR+BV)/2.d0
  GH=(GR+GV)/2.d0
  write(InVar%stdout,'(2(a,f9.3))')'ISOTHERMAL modulus: Bulk Kt=',BH,' and Shear G=',GH
  Eaverage=9.d0*BH*GH/(3*BH+GH)
  Nuaverage=0.5*(1.d0-(3.d0*GH)/(3.d0*BH+GH) )
  Laverage=(3.d0*BH-2.d0*GH)/3.d0
  write(InVar%stdout,'(3(a,f9.3))')'Average of Young modulus E=',Eaverage,' Lame modulus Lambda=',Laverage,' and Poisson ratio Nu=',Nuaverage
  Vp=dsqrt(1.d9*(BH+4.d0*GH/3.d0)/rho)
  Vs=dsqrt(1.d9*GH/rho)
  Vphi=dsqrt(1.d9*BH/rho)
  write(InVar%stdout,'(3(a,f9.3,x))')'Velocities: compressional Vp=',Vp,' shear Vs=',Vs,' and bulk Vphi=',Vphi

  ABI_FREE(cijkl)
  ABI_FREE(Cij)
  ABI_FREE(Sij)
end subroutine  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module m_tdepdos
