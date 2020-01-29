
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_readwrite
  
  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  use m_abihist
  use m_abimover, only : abimover

 implicit none

  type Input_Variables_type

    integer :: Impose_Symetry=0
    integer :: natom
    integer :: natom_unitcell
    integer :: nstep_max
    integer :: nstep_min
    integer :: nstep
    integer :: ntypat
    integer :: Use_ideal_positions
    integer :: stdout
    integer :: stdlog
    integer :: BZpath
    integer :: Order
    integer :: Slice
    integer :: Enunit
    integer :: ReadIFC
    integer :: firstqptseg
    integer :: ngqpt1(3)
    integer :: ngqpt2(3)
    integer :: bravais(11)
    integer, allocatable ::typat_unitcell(:)
    integer, allocatable ::typat(:)
    logical :: debug
    logical :: loto
    logical :: netcdf
    double precision :: angle_alpha
    double precision :: dielec_constant
    double precision :: dosdeltae
    double precision :: Rcut
    double precision :: Rcut3
    double precision :: temperature
    double precision :: tolread
    double precision :: tolinbox
    double precision :: tolmatch
    double precision :: tolmotif
    double precision :: rprimd_MD(3,3)
    double precision :: multiplicity(3,3)
    double precision, allocatable :: amu(:)
    double precision, allocatable :: born_charge(:)
    double precision, allocatable :: qpt(:,:)
    double precision, allocatable :: xred_ideal(:,:)
    double precision, allocatable :: xred_unitcell(:,:)
    double precision, allocatable :: xred(:,:,:)
    double precision, allocatable :: fcart(:,:,:)
    double precision, allocatable :: etot(:)
!FB    double precision, allocatable :: sigma(:,:)
    character (len=2), allocatable :: special_qpt(:)
    character (len=200) :: output_prefix
    
  end type Input_Variables_type

!FB  type, public :: Hist_type
!FB
!FB  ! scalars
!FB    ! Index of the last element on all records
!FB    integer :: ihist = 0
!FB    ! Maximun size of the historical records
!FB    integer :: mxhist = 0
!FB    ! Booleans to know if some arrays are changing
!FB    logical :: isVused  ! If velocities are changing
!FB    logical :: isARused ! If Acell and Rprimd are changing
!FB
!FB  ! arrays
!FB    ! Vector of (x,y,z)X(mxhist)
!FB    real(dp), allocatable :: histA(:,:)
!FB    ! Vector of (mxhist) values of energy
!FB    real(dp), allocatable :: histE(:)
!FB    ! Vector of (mxhist) values of ionic kinetic energy
!FB    real(dp), allocatable :: histEk(:)
!FB    ! Vector of (mxhist) values of Entropy
!FB    real(dp), allocatable :: histEnt(:)
!FB    ! Vector of (mxhist) values of time (relevant
!FB    ! for MD calculations)
!FB    real(dp), allocatable :: histT(:)
!FB    ! Vector of (x,y,z)X(x,y,z)X(mxhist)
!FB    real(dp), allocatable :: histR(:,:,:)
!FB    ! Vector of (stress [6])X(mxhist)
!FB    real(dp), allocatable :: histS(:,:)
!FB    ! Vector of (x,y,z)X(natom)X(mxhist) values of velocity
!FB    real(dp), allocatable :: histV(:,:,:)
!FB    ! Vector of (x,y,z)X(natom)X(xcart,xred,fcart,fred)X(mxhist)
!FB    real(dp), allocatable :: histXF(:,:,:,:)
!FB
!FB  end type Hist_type

 public :: tdep_print_Aknowledgments
 public :: tdep_ReadEcho

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_print_Aknowledgments(InVar)

  type(Input_Variables_type) :: InVar
  integer :: stdout
  stdout = InVar%stdout

  write(stdout,*) ' '
  write(stdout,'(a)') ' #############################################################################'
  write(stdout,'(a)') ' ######################### CALCULATION COMPLETED #############################'
  write(stdout,'(a)') ' #############################################################################'
  write(stdout,'(a)') ' Suggested references for the acknowledgment of ABINIT usage.'
  write(stdout,'(a)') ' '
  write(stdout,'(a)') ' The users of ABINIT have little formal obligations with respect to the ABINIT group'
  write(stdout,'(a)') ' (those specified in the GNU General Public License, http://www.gnu.org/copyleft/gpl.txt).'
  write(stdout,'(a)') ' However, it is common practice in the scientific literature,'
  write(stdout,'(a)') ' to acknowledge the efforts of people that have made the research possible.'
  write(stdout,'(a)') ' In this spirit, please find below suggested citations of work written by ABINIT developers,'
  write(stdout,'(a)') ' corresponding to implementations inside of ABINIT that you have used in the present run.'
  write(stdout,'(a)') ' Note also that it will be of great value to readers of publications presenting these results,'
  write(stdout,'(a)') ' to read papers enabling them to understand the theoretical formalism and details'
  write(stdout,'(a)') ' of the ABINIT implementation.'
  write(stdout,'(a)') ' For information on why they are suggested, see also https://docs.abinit.org/theory/acknowledgments.'
  write(stdout,'(a)') ' '
  write(stdout,'(a)') '.[1] Thermal evolution of vibrational properties of $\\alpha$-U' 
  write(stdout,'(a)') ' J. Bouchet and F. Bottin, Phys. Rev. B 92, 174108 (2015).' ! [[cite:Bouchet2015]]
  write(stdout,'(a)') ' Strong suggestion to cite this paper in your publications.'
  write(stdout,'(a)') ' This paper is also available at http://www.arxiv.org/abs/xxxx'
  write(stdout,'(a)') ' '
  write(stdout,'(a)') ' [2] Lattice dynamics of anharmonic solids from first principles'
  write(stdout,'(a)') ' O. Hellman and I.A. Abrikosov and S.I. Simak, Phys. Rev. B 84, 180301(R) (2011).' ! [[cite:Hellman2011]]
  write(stdout,'(a)') ' Strong suggestion to cite this paper in your publications.'
  write(stdout,'(a)') ' '
  write(stdout,'(a)') ' [3] Temperature dependent effective potential method for accurate free energy calculations of solids'
  write(stdout,'(a)') ' O. Hellman and P. Steneteg and I.A. Abrikosov and S.I. Simak, Phys. Rev. B 87, 104111 (2013).' ! [[cite:Hellman2013]]
  write(stdout,'(a)') ' Strong suggestion to cite this paper in your publications.'

 end subroutine tdep_print_Aknowledgments 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_ReadEcho(InVar)

#if defined HAVE_NETCDF
 use netcdf
#endif

  integer :: ii,jj,tmp,istep,iatom,this_istep
  character (len=30):: string,NormalMode,DebugMode,Impose_Symetry,Use_ideal_positions
  character (len=30):: Born_charge,Dielec_constant,tolmotifinboxmatch,TheEnd,BZpath
  character (len=30):: Order,Slice,Enunit,ReadIFC,firstqptseg,Ngqpt1,Ngqpt2,DosDeltae
  double precision :: version_value,tmp1,tmp2,tmp3,nstep_int,nstep_float
  character (len=8) :: date
  character (len=10) :: time
  character (len=5) :: zone
  character(len=3),parameter :: month_names(12)=(/'Jan','Feb','Mar','Apr','May','Jun',&
&                                                 'Jul','Aug','Sep','Oct','Nov','Dec'/)
  character(len=500) :: filename
  character(len=500) :: inputfilename
  integer :: values(8)  
  type(Input_Variables_type),intent(out) :: InVar
  type(abihist) :: Hist
  ! Temp variable to get dimensions of HIST file
  integer :: ncid, ncerr
  integer :: nimage, mdtime, natom_id,nimage_id,time_id,xyz_id,six_id
  integer :: ntypat_id
  logical :: has_nimage
  real(dp) :: dtion
  real(dp), allocatable :: znucl(:)

! Define output files  
  InVar%stdout=7
  InVar%stdlog=6
  !open(unit=InVar%stdlog,file='data.log')

! Define Keywords
  NormalMode='NormalMode'
  DebugMode='DebugMode'
  Impose_Symetry='Impose_Symetry'
  Use_Ideal_Positions='Use_Ideal_Positions'
  Born_Charge='Born_Charge'
  Dielec_Constant='Dielec_Constant'
  DosDeltae='DosDeltae'
  BZpath='BZpath'
  Firstqptseg='Firstqptseg'
  Order='Order'
  Slice='Slice'
  Enunit='Enunit'
  ReadIFC='ReadIFC'
  Ngqpt1='Ngqpt1'
  Ngqpt2='Ngqpt2'
  TolMotifInboxMatch='TolMotifInboxMatch'
  TheEnd='TheEnd'
! Define default values
  InVar%angle_alpha=90.d0
  InVar%BZpath=0
  InVar%Order=2
  InVar%Slice=1
  InVar%Enunit=0
  InVar%ReadIFC=0
  InVar%firstqptseg=100
  InVar%tolread=1.d-8
  InVar%tolmotif=5.d-2
  InVar%tolinbox=5.d-2
  InVar%tolmatch=5.d-2
  InVar%dosdeltae=4.5d-6
  InVar%debug=.false.
  InVar%loto=.false.
  InVar%netcdf=.false.
  InVar%Use_ideal_positions=0
  version_value=2.d0
! In order to have an accuracy better than 1meV  
  InVar%ngqpt1(:)=8
  InVar%ngqpt2(:)=32

! Check if a NetCDF file is available
  filename='HIST.nc'
  inputfilename='input.in'
  write(InVar%stdlog,'(a)',err=10) ' Give name for input file '
  read(*, '(a)',err=10) inputfilename
  write(InVar%stdlog, '(a)',err=10) '.'//trim(inputfilename)
10 continue
  write(InVar%stdlog,'(a)',err=11) ' Give name for HIST file '
  read(*, '(a)',err=11) filename
  write(InVar%stdlog, '(a)',err=11) '.'//trim(filename)
11 continue
  write(InVar%stdlog,'(a)', err=12)' Give root name for generic output files:'
  read (*, '(a)', err=12) InVar%output_prefix
  write (InVar%stdlog, '(a)', err=12 ) InVar%output_prefix 
12 continue
  if ( inputfilename == "" ) inputfilename='input.in'
  if ( filename == "" ) filename='HIST.nc'

  open(unit=InVar%stdout,file=trim(InVar%output_prefix)//'.out')


#if defined HAVE_NETCDF
 !Open netCDF file
  ncerr=nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=ncid)
  if(ncerr /= NF90_NOERR) then
    write(InVar%stdout,'(3a)') '-'//'Could no open ',trim(filename),', starting from scratch'
    InVar%netcdf=.false.
  else
    write(InVar%stdout,'(3a)') '-'//'Succesfully open ',trim(filename),' for reading'
    write(InVar%stdout,'(a)') ' Extracting information from NetCDF file...'
    InVar%netcdf=.true.
  end if

  if ( InVar%netcdf) then
    call get_dims_hist(ncid,InVar%natom,InVar%ntypat,nimage,mdtime,&
&       natom_id,ntypat_id,nimage_id,time_id,xyz_id,six_id,has_nimage)
    ABI_MALLOC(InVar%amu,(InVar%ntypat)); InVar%amu(:)=zero
    ABI_MALLOC(InVar%typat,(InVar%natom)); InVar%typat(:)=zero
    ABI_MALLOC(znucl,(InVar%ntypat)) ; znucl(:)=zero
    call read_csts_hist(ncid,dtion,InVar%typat,znucl,InVar%amu)
    ABI_FREE(znucl)

    ! Need to close NetCDF file because it is going to be reopened by read_md_hist
    ncerr = nf90_close(ncid)
    ! .false. -> Velocities note used
    ! .true. -> acell and rprimd may change (2017_04 only NVT/isoK used but maybe
    ! .false. -> read all times
    ! NPT one day ?)
    call read_md_hist(filename,Hist,.false.,.true.,.false.)
  end if
#endif

! Write version, copyright, date...
  write(InVar%stdout,*) ' '
  open(unit=40,file=inputfilename)
  read(40,*) string
  if (string.eq.NormalMode) then
    write(InVar%stdout,'(a,f6.1,a)') '.Version ', version_value,' of PHONONS'
  else if (string.eq.DebugMode) then
    InVar%debug=.true.
    write(InVar%stdout,'(a,f6.1,a)') '.Version ', version_value,' of PHONONS (Debug)'
  else
    MSG_ERROR('Please use recent format for the input file')
  end if  
  write(InVar%stdout,'(a)') '.Copyright (C) 1998-2020 ABINIT group (FB,JB).'
  write(InVar%stdout,'(a)') ' ABINIT comes with ABSOLUTELY NO WARRANTY.'
  write(InVar%stdout,'(a)') ' It is free software, and you are welcome to redistribute it'
  write(InVar%stdout,'(a)') ' under certain conditions (GNU General Public License,'
  write(InVar%stdout,'(a)') ' see ~abinit/COPYING or http://www.gnu.org/copyleft/gpl.txt).'
  write(InVar%stdout,*) ' '
  write(InVar%stdout,'(a)') ' ABINIT is a project of the Universite Catholique de Louvain,'
  write(InVar%stdout,'(a)') ' Corning Inc. and other collaborators, see'
  write(InVar%stdout,'(a)') ' ~abinit/doc/developers/contributors.txt .'
  write(InVar%stdout,'(a)') ' Please read https://docs.abinit.org/theory/acknowledgments for suggested'
  write(InVar%stdout,'(a)') ' acknowledgments of the ABINIT effort.'
  write(InVar%stdout,'(a)') ' For more information, see http://www.abinit.org .'

  call date_and_time(date,time,zone,values)
  write(InVar%stdout,'(/,a,i2,1x,a,1x,i4,a)') '.Starting date : ',values(3),month_names(values(2)),values(1),'.'

! Read (and echo) of input variables from the input.in input file
  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '######################### ECHO OF INPUT FILE ################################'
  write(InVar%stdout,*) '#############################################################################'
! Define unit cell  
  read(40,*) string
  write(InVar%stdout,'(a)') ' ======================= Define the unitcell =================================' 
  read(40,*) string,InVar%bravais(1),InVar%bravais(2)
  write(InVar%stdout,'(1x,a20,1x,i4,1x,i4)') string,InVar%bravais(1),InVar%bravais(2)
  if ((InVar%bravais(1).eq.2).or.(InVar%bravais(1).eq.5)) then
    read(40,*) string,InVar%angle_alpha
    write(InVar%stdout,'(1x,a20,1x,f15.10)') string,InVar%angle_alpha
  else
    !read(40,*)
    InVar%angle_alpha=90.d0
  end if
  read(40,*) string,InVar%natom_unitcell
  write(InVar%stdout,'(1x,a20,1x,i4)') string,InVar%natom_unitcell
  ABI_MALLOC(InVar%xred_unitcell,(3,InVar%natom_unitcell)); InVar%xred_unitcell(:,:)=zero
  read(40,*) string,InVar%xred_unitcell(:,:)
  write(InVar%stdout,'(1x,a20)') string
  do ii=1,InVar%natom_unitcell
    write(InVar%stdout,'(22x,3(f15.10,1x))') (InVar%xred_unitcell(jj,ii), jj=1,3)
  end do  
  ABI_MALLOC(InVar%typat_unitcell,(InVar%natom_unitcell)); InVar%typat_unitcell(:)=0 
  read(40,*) string,InVar%typat_unitcell(:)
  write(InVar%stdout,'(1x,a20,20(1x,i4))') string,(InVar%typat_unitcell(jj),jj=1,InVar%natom_unitcell)
  if (InVar%netcdf) then
    string='ntypat'
  else
    read(40,*) string,InVar%ntypat
  end if  
  write(InVar%stdout,'(1x,a20,1x,i4)') string,InVar%ntypat
  if (InVar%netcdf) then
    string='amu'
  else  
    ABI_MALLOC(InVar%amu,(InVar%ntypat)); InVar%amu(:)=zero
    read(40,*) string,InVar%amu(:)
  end if  
  write(InVar%stdout,'(1x,a20,20(1x,f15.10))') string,(InVar%amu(jj),jj=1,InVar%ntypat)
! Define supercell (as a function of the unitcell defined above)
  read(40,*) string
  write(InVar%stdout,'(a)') ' ======================= Define the supercell ================================' 
  if (InVar%netcdf) then
    InVar%rprimd_MD(:,:)=Hist%rprimd(:,:,Hist%ihist)
    string='rprimd'
  else
    read(40,*) string,InVar%rprimd_MD(1,:),InVar%rprimd_MD(2,:),InVar%rprimd_MD(3,:)
  end if  
  write(InVar%stdout,'(1x,a20)') string
  do ii=1,3
    write(InVar%stdout,'(22x,3(f15.10,1x))') (InVar%rprimd_MD(ii,jj),jj=1,3)
  end do  
  read(40,*) string,InVar%multiplicity(1,:),InVar%multiplicity(2,:),InVar%multiplicity(3,:)
  write(InVar%stdout,'(1x,a20)') string
  do ii=1,3
    write(InVar%stdout,'(22x,3(f15.10,1x))') (InVar%multiplicity(ii,jj),jj=1,3)
  end do  
  if (InVar%netcdf) then
    string='natom'
  else
    read(40,*) string,InVar%natom
  end if
  write(InVar%stdout,'(1x,a20,1x,i4)') string,InVar%natom
  if (InVar%netcdf) then
    string='typat'
  else
    ABI_MALLOC(InVar%typat,(InVar%natom)); InVar%typat(:)=0 
    read(40,*) string,InVar%typat(:)
  end if  
  write(InVar%stdout,'(1x,a20)') string
  do ii=1,InVar%natom,10
    if (ii+9.lt.InVar%natom) then
      write(InVar%stdout,'(22x,10(i4,1x))') (InVar%typat(ii+jj-1),jj=1,10)
    else
      write(InVar%stdout,'(22x,10(i4,1x))') (InVar%typat(jj),jj=ii,InVar%natom)
    end if  
  end do  
  read(40,*) string,InVar%temperature
  write(InVar%stdout,'(1x,a20,1x,f15.10)') string,InVar%temperature
! Define phonons computational details
  read(40,*) string
  write(InVar%stdout,'(a)') ' ======================= Define computational details ========================' 
  read(40,*) string,InVar%nstep_max
  write(InVar%stdout,'(1x,a20,1x,i5)') string,InVar%nstep_max
  read(40,*) string,InVar%nstep_min
  write(InVar%stdout,'(1x,a20,1x,i5)') string,InVar%nstep_min
  read(40,*) string,InVar%Rcut
  write(InVar%stdout,'(1x,a20,1x,f15.10)') string,InVar%Rcut
! Optional input variables  
  read(40,*) string
  write(InVar%stdout,'(a)') ' ======================= Optional input variables ============================' 
  do ii=1,100
    read(40,*) string
    backspace(40)
    if (string.eq.DosDeltae) then
      read(40,*) string,InVar%dosdeltae
      write(InVar%stdout,'(1x,a20,1x,f15.10)') string,InVar%dosdeltae
    else if (string.eq.Impose_Symetry) then
      read(40,*) string,InVar%Impose_Symetry
      write(InVar%stdout,'(1x,a20,1x,i4)') string,InVar%Impose_Symetry
    else if (string.eq.Use_ideal_positions) then  
      read(40,*) string,InVar%Use_ideal_positions
      write(InVar%stdout,'(1x,a20,1x,i4)') string,InVar%Use_ideal_positions
    else if (string.eq.Born_charge) then  
      ABI_MALLOC(InVar%born_charge,(InVar%ntypat)); InVar%born_charge(:)=0.d0
      InVar%loto=.true.
      read(40,*) string,InVar%born_charge(:)
      write(InVar%stdout,'(1x,a20,20(1x,f15.10))') string,(InVar%born_charge(jj),jj=1,InVar%ntypat)
    else if (string.eq.Dielec_constant) then  
      read(40,*) string,InVar%dielec_constant
      write(InVar%stdout,'(1x,a20,1x,f15.10)') string,InVar%dielec_constant
    else if (string.eq.BZpath) then  
      read(40,*) string,InVar%BZpath
      write(InVar%stdout,'(1x,a20,1x,i4)') string,InVar%BZpath
      if (InVar%BZpath.lt.0) then
        ABI_MALLOC(InVar%qpt,(3,abs(InVar%BZpath))); InVar%qpt(:,:)=zero
        write(InVar%stdout,'(a)') ' Q points as given in the input file:'
        do jj=1,abs(InVar%BZpath)
          read(40,*) InVar%qpt(:,jj)
          write(InVar%stdout,'(22x,3(f15.10,1x))') InVar%qpt(:,jj)
        end do  
      else if (InVar%BZpath.gt.0) then
        ABI_MALLOC(InVar%special_qpt,(InVar%BZpath))
        backspace(40)
        read(40,*) string,tmp,(InVar%special_qpt(jj),jj=1,InVar%BZpath)
        write(InVar%stdout,'(a,1x,10(a2,"-"))') ' Special q-points: ',InVar%special_qpt(:)
      end if
    else if (string.eq.Order) then  
      read(40,*) string,InVar%Order,InVar%Rcut3
      write(InVar%stdout,'(1x,a20,1x,i4,1x,f15.10)') string,InVar%Order,InVar%Rcut3
      if (InVar%Rcut3.gt.InVar%Rcut) then
        MSG_ERROR('The cutoff radius of the third order cannot be greater than the second order one.')
      end if  
    else if (string.eq.Slice) then  
      read(40,*) string,InVar%Slice
      write(InVar%stdout,'(1x,a20,1x,i4)') string,InVar%Slice
      nstep_float=float(InVar%nstep_max-InVar%nstep_min+1)/float(InVar%Slice)
      nstep_int  =float(int(nstep_float))
      write(InVar%stdout,*) nstep_int,nstep_float
      if (abs(nstep_float-nstep_int).gt.tol8) then
        MSG_ERROR('Change nstep_min. (nstep_max-nstep_min+1)/Slice has to be an integer.')
      end if  
    else if (string.eq.Enunit) then  
      read(40,*) string,InVar%Enunit
      if (InVar%Enunit.eq.0) write(InVar%stdout,'(1x,a20,1x,i4,1x,a)') string,InVar%Enunit,'(energy in meV)'
      if (InVar%Enunit.eq.1) write(InVar%stdout,'(1x,a20,1x,i4,1x,a)') string,InVar%Enunit,'(energy in cm-1)'
      if (InVar%Enunit.eq.2) write(InVar%stdout,'(1x,a20,1x,i4,1x,a)') string,InVar%Enunit,'(energy in Ha)'
    else if (string.eq.ReadIFC) then  
      read(40,*) string,InVar%ReadIFC
      if (InVar%ReadIFC.eq.1) then
        backspace(40)
        read(40,*) string,InVar%ReadIFC,InVar%tolread
        write(InVar%stdout,'(1x,a20,1x,i4,1x,f15.10)') string,InVar%ReadIFC,InVar%tolread
      else  
        write(InVar%stdout,'(1x,a20,1x,i4)') string,InVar%ReadIFC
      end if  
    else if (string.eq.Firstqptseg) then  
      read(40,*) string,InVar%firstqptseg
      write(InVar%stdout,'(1x,a20,1x,i4)') string,InVar%firstqptseg
    else if (string.eq.Ngqpt1) then  
      read(40,*) string,InVar%ngqpt1(:)
      write(InVar%stdout,'(1x,a20,1x,3(i4,1x))') string,InVar%ngqpt1(:)
    else if (string.eq.Ngqpt2) then  
      read(40,*) string,InVar%ngqpt2(:)
      write(InVar%stdout,'(1x,a20,1x,3(i4,1x))') string,InVar%ngqpt2(:)
    else if (string.eq.tolmotifinboxmatch) then  
      read(40,*) string,InVar%tolmotif,InVar%tolinbox,InVar%tolmatch
      write(InVar%stdout,'(1x,a20,f10.5)') 'tolmotif            ',InVar%tolmotif
      write(InVar%stdout,'(1x,a20,f10.5)') 'tolinbox            ',InVar%tolinbox
      write(InVar%stdout,'(1x,a20,f10.5)') 'tolmatch            ',InVar%tolmatch
    else if (string.eq.TheEnd) then
      exit
    else 
      write(InVar%stdout,'(a,1x,a)') 'This keyword is not allowed',string
      MSG_ERROR('A keyword is not allowed. See the log file.')
    end if  
  end do
! Output very important informations 
  write(InVar%stdout,'(a)') ' '
  if (InVar%Impose_Symetry.eq.0) then
    write(InVar%stdout,'(a)') ' STOP IF THE DIJ (IFC) MATRIX IS NOT HERMITIAN (SYMETRIC)'
  else if (InVar%Impose_Symetry.eq.1) then
    write(InVar%stdout,'(a)') ' SYMETRIZE THE DIJ MATRIX (HERMITIAN)'
  else if (InVar%Impose_Symetry.eq.2) then
    write(InVar%stdout,'(a)') ' SYMETRIZE THE IFC MATRIX (SYMETRIC)'
  else if (InVar%Impose_Symetry.eq.3) then
    write(InVar%stdout,'(a)') ' SYMETRIZE THE DIJ MATRIX (HERMITIAN) AND THE IFC MATRIX (SYMETRIC)'
  else
    write(InVar%stdout,'(a)') ' STOP: THIS VALUE IS NOT ALLOWED FOR Impose_Symetry'
  end if
  if (InVar%Use_ideal_positions.eq.0) then
    write(InVar%stdout,'(a)') ' USE AVERAGE POSITIONS TO COMPUTE SPECTRUM'
  else if (InVar%Use_ideal_positions.eq.1) then
    write(InVar%stdout,'(a)') ' USE IDEAL POSITIONS TO COMPUTE SPECTRUM'
  else
    write(InVar%stdout,'(a)') ' STOP: THIS VALUE IS NOT ALLOWED FOR Use_Ideal_Positions'
  end if
  if (InVar%loto) write(InVar%stdout,'(a)') ' USE NON-ANALYTICAL CORRECTIONS (LO-TO)'
  InVar%nstep=(InVar%nstep_max-InVar%nstep_min+1)/InVar%Slice
  write(InVar%stdout,'(a)') ' '
  write(InVar%stdout,'(a)') ' WARNING: ALL the quantities are now computed :'
  write(InVar%stdout,'(a,1x,i4)') '                                      from nstep_min=',InVar%nstep_min
  write(InVar%stdout,'(a,1x,i4)') '                                        to nstep_max=',InVar%nstep_max
  if (InVar%Slice.ne.1) then
    write(InVar%stdout,'(a,1x,i4)') '                                    by using a slice=',InVar%Slice
  end if  
  write(InVar%stdout,'(a,1x,i4)') '          So, the real number of time steps is nstep=',InVar%nstep
! End of read and echo  
  close(40)

! Read xred.dat, fcart.dat and etot.dat ASCII files or extract them from the HIST.nc netcdf file.
  write(InVar%stdout,'(a)') ' '
  ABI_MALLOC(InVar%xred,(3,InVar%natom,InVar%nstep))  ; InVar%xred(:,:,:)=0.d0
  ABI_MALLOC(InVar%fcart,(3,InVar%natom,InVar%nstep)) ; InVar%fcart(:,:,:)=0.d0
  ABI_MALLOC(InVar%etot,(InVar%nstep))                ; InVar%etot(:)=0.d0
!FB  ABI_MALLOC(InVar%sigma,(6,InVar%nstep))             ; InVar%sigma(:,:)=0.d0
  if (InVar%netcdf) then
    this_istep=0
    do istep=1,InVar%nstep_max
      if ((istep.lt.InVar%nstep_min).or.(mod(istep-InVar%nstep_min,InVar%Slice).ne.0)) then
        cycle
      else
        this_istep=this_istep+1
        InVar%xred(:,:,this_istep) =Hist%xred (:,:,istep)
        InVar%fcart(:,:,this_istep)=Hist%fcart(:,:,istep)
        InVar%etot(this_istep)     =Hist%etot     (istep)
!FB        InVar%sigma(:,this_istep)  =Hist%strten (:,istep)
      end if
    end do !istep  
    write(InVar%stdout,'(a)') ' The Xred, Fcart, Etot and Stress data are extracted from the NetCDF file: HIST.nc '
  else
    open(unit=60,file='fcart.dat')
    open(unit=50,file='xred.dat')
    open(unit=40,file='etot.dat')
    this_istep=0
    do istep=1,InVar%nstep_max
      if ((istep.lt.InVar%nstep_min).or.(mod(istep-InVar%nstep_min,InVar%Slice).ne.0)) then
        read(40,*) tmp1
      else 
        this_istep=this_istep+1
        read(40,*) InVar%etot(this_istep)
      end if  
      do iatom=1,InVar%natom
        if ((istep.lt.InVar%nstep_min).or.(mod(istep-InVar%nstep_min,InVar%Slice).ne.0)) then
          read(50,*) tmp1,tmp2,tmp3
          read(60,*) tmp1,tmp2,tmp3
        else 
          read(50,*) InVar%xred (1,iatom,this_istep),InVar%xred (2,iatom,this_istep),InVar%xred (3,iatom,this_istep)
          read(60,*) InVar%fcart(1,iatom,this_istep),InVar%fcart(2,iatom,this_istep),InVar%fcart(3,iatom,this_istep)
        end if
      end do
    end do !istep 
    close(40)
    close(50)
    close(60)
    write(InVar%stdout,'(a)') ' The Xred, Fcart and Etot data are extracted from the ASCII files: xred.dat, fcart.dat \& etot.dat'
  end if  
 end subroutine tdep_ReadEcho
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module m_tdep_readwrite
