
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

    integer :: natom
    integer :: natom_unitcell
    integer :: nstep_max
    integer :: nstep_min
    integer :: nstep_tot
    integer :: my_nstep
    integer :: ntypat
    integer :: Use_ideal_positions
    integer :: stdout
    integer :: stdlog
    integer :: BZpath
    integer :: Order
    integer :: Slice
    integer :: Enunit
    integer :: ReadIFC
    integer :: together
    integer :: Nproc(2)
    integer :: BZlength
    integer :: ngqpt1(3)
    integer :: ngqpt2(3)
    integer :: bravais(11)
    integer, allocatable :: typat_unitcell(:)
    integer, allocatable :: typat(:)
    integer, allocatable :: lgth_segments(:)
    logical :: debug
    logical :: loto
    logical :: netcdf
    double precision :: angle_alpha
    double precision :: dielec_constant
    double precision :: dosdeltae
    double precision :: Rcut
    double precision :: Rcut3
    double precision :: Rcut4
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
    character (len=2), allocatable :: special_qpt(:)
    character (len=200) :: output_prefix
    
  end type Input_Variables_type

  type MPI_enreg_type

    integer :: comm_shell
    integer :: comm_step
    integer :: comm_shellstep
    integer :: nproc
    integer :: nproc_shell
    integer :: nproc_step
    integer :: master
    integer :: me_shell
    integer :: me_step
    integer, allocatable :: my_nshell(:)
    integer :: my_nstep
    logical :: iam_master
    integer, allocatable :: nstep_all(:)
    integer, allocatable :: shft_step(:)
    logical, allocatable :: my_shell(:)
    logical, allocatable :: my_step(:)

  end type MPI_enreg_type

 public :: tdep_print_Aknowledgments
 public :: tdep_read_input
 public :: tdep_distrib_data
 public :: tdep_init_MPIdata
 public :: tdep_init_MPIshell
 public :: tdep_clean_MPI

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
  write(stdout,'(a)') ' [1] Thermal evolution of vibrational properties of alpha-U' 
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
 subroutine tdep_read_input(Hist,InVar)

#if defined HAVE_NETCDF
 use netcdf
#endif

  type(Input_Variables_type),intent(out) :: InVar
  type(abihist), intent(out) :: Hist

  integer :: values(8)  
  integer :: ncid, ncerr,ierr
  integer :: nimage, mdtime, natom_id,nimage_id,time_id,xyz_id,six_id
  integer :: ntypat_id
  integer :: ii,jj,tmp,istep,iatom,this_istep
  double precision :: version_value,tmp1,tmp2,tmp3,dtion
  character (len=30):: string,NormalMode,DebugMode,Use_ideal_positions
  character (len=30):: Born_charge,Dielec_constant,tolmotifinboxmatch,TheEnd,BZpath
  character (len=30):: Order,Slice,Enunit,ReadIFC,together,Nproc,BZlength,Ngqpt1,Ngqpt2,DosDeltae
  character (len=8) :: date
  character (len=10) :: time
  character (len=5) :: zone
  character(len=3),parameter :: month_names(12)=(/'Jan','Feb','Mar','Apr','May','Jun',&
&                                                 'Jul','Aug','Sep','Oct','Nov','Dec'/)
  character(len=500) :: filename,inputfilename,message
  logical :: has_nimage
  real(dp), allocatable :: znucl(:)

! Define output files  
  InVar%stdout=8
  InVar%stdlog=6
  !open(unit=InVar%stdlog,file='data.log')

! Define Keywords
  NormalMode='NormalMode'
  DebugMode='DebugMode'
  Use_Ideal_Positions='Use_Ideal_Positions'
  Born_Charge='Born_Charge'
  Dielec_Constant='Dielec_Constant'
  DosDeltae='DosDeltae'
  BZpath='BZpath'
  BZlength='BZlength'
  Order='Order'
  Slice='Slice'
  Enunit='Enunit'
  ReadIFC='ReadIFC'
  together='together'
  Nproc='Nproc'
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
  InVar%together=1
  InVar%Nproc(:)=1
  InVar%BZlength=0
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
  if (xmpi_comm_rank(xmpi_world).eq.0) open(unit=InVar%stdout,file=trim(InVar%output_prefix)//'.out')

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
  write(InVar%stdout,'(a)') '.Copyright (C) 1998-2019 ABINIT group (FB,JB).'
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
! To avoid some troubles and inconsistency between the .in and .nc files, when we use NetCDF or not.    
  if (InVar%netcdf) then
    read(40,*) string
    backspace(40)
    if (string.eq.'ntypat') then
      write(InVar%stdlog,'(1x,a)') 'When the NetCDF file .nc is used, the ntypat keywork is not allowed.' 
      MSG_ERROR('ACTION : Please modify your input file')
    end if  
    string='ntypat'
  else
    read(40,*) string
    backspace(40)
    if (string.ne.'ntypat') then
      write(InVar%stdlog,'(1x,a)') 'The NetCDF file .nc is not used.' 
      write(InVar%stdlog,'(1x,a,1x,a)') 'In your input file, the code search the ntypat keywork but found :',string
      MSG_ERROR('ACTION : Please modify your input file')
    end if  
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
        write(InVar%stdout,'(a,1x,a2,10("-",a2))') ' Special q-points: ',InVar%special_qpt(:)
      end if
    else if (string.eq.Order) then  
      read(40,*) string,InVar%Order
      backspace(40)
      if (InVar%Order.eq.3) then
        read(40,*) string,InVar%Order,InVar%Rcut3
        write(InVar%stdout,'(1x,a20,1x,i4,1x,f15.10)') string,InVar%Order,InVar%Rcut3
        if (InVar%Rcut3.gt.InVar%Rcut) then
          MSG_ERROR('The cutoff radius of the third order cannot be greater than the second order one.')
        end if  
      else if (InVar%Order.eq.4) then
        read(40,*) string,InVar%Order,InVar%Rcut3,InVar%Rcut4
        write(InVar%stdout,'(1x,a20,1x,i4,2(1x,f15.10))') string,InVar%Order,InVar%Rcut3,InVar%Rcut4
        if (InVar%Rcut4.gt.InVar%Rcut) then
          MSG_ERROR('The cutoff radius of the fourth order cannot be greater than the second order one.')
        end if  
      else
        MSG_ERROR('Only the 3rd and 4th orders are allowed. Change your input file.')
      end if
    else if (string.eq.Slice) then  
      read(40,*) string,InVar%Slice
      write(InVar%stdout,'(1x,a20,1x,i4)') string,InVar%Slice
    else if (string.eq.Enunit) then  
      read(40,*) string,InVar%Enunit
      if (InVar%Enunit.eq.0) write(InVar%stdout,'(1x,a20,1x,i4,1x,a)') string,InVar%Enunit,'(Phonon frequencies in meV)'
      if (InVar%Enunit.eq.1) write(InVar%stdout,'(1x,a20,1x,i4,1x,a)') string,InVar%Enunit,'(Phonon frequencies in cm-1)'
      if (InVar%Enunit.eq.2) write(InVar%stdout,'(1x,a20,1x,i4,1x,a)') string,InVar%Enunit,'(Phonon frequencies in Ha)'
    else if (string.eq.Nproc) then  
      read(40,*) string,InVar%Nproc(1),InVar%Nproc(2)
      write(InVar%stdout,'(1x,a20,1x,i4,1x,i4)') string,InVar%Nproc(1),InVar%Nproc(2)
    else if (string.eq.ReadIFC) then  
      read(40,*) string,InVar%ReadIFC
      if (InVar%ReadIFC.eq.1) then
        backspace(40)
        read(40,*) string,InVar%ReadIFC,InVar%tolread
        write(InVar%stdout,'(1x,a20,1x,i4,1x,f15.10)') string,InVar%ReadIFC,InVar%tolread
      else  
        write(InVar%stdout,'(1x,a20,1x,i4)') string,InVar%ReadIFC
      end if  
    else if (string.eq.together) then  
      read(40,*) string,InVar%together
      write(InVar%stdout,'(1x,a20,1x,i4)') string,InVar%together
    else if (string.eq.BZlength) then  
      read(40,*) string,InVar%BZlength
      write(InVar%stdout,'(1x,a20,1x,i4)') string,InVar%BZlength
      ABI_MALLOC(InVar%lgth_segments,(InVar%BZlength))
      backspace(40)
      read(40,*) string,tmp,(InVar%lgth_segments(jj),jj=1,InVar%BZlength)
      write(InVar%stdout,'(a,1x,i3,10("-",i3))') ' Length of BZ : ',InVar%lgth_segments(:)
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
  if (InVar%Use_ideal_positions.eq.0) then
    write(InVar%stdout,'(a)') ' USE AVERAGE POSITIONS TO COMPUTE SPECTRUM'
  else if (InVar%Use_ideal_positions.eq.1) then
    write(InVar%stdout,'(a)') ' USE IDEAL POSITIONS TO COMPUTE SPECTRUM'
  else
    write(InVar%stdout,'(a)') ' STOP: THIS VALUE IS NOT ALLOWED FOR Use_Ideal_Positions'
  end if
  if (InVar%loto) write(InVar%stdout,'(a)') ' USE NON-ANALYTICAL CORRECTIONS (LO-TO)'

! Allowed values
  if ((InVar%together.ne.1).and.(InVar%together.ne.0)) then
    MSG_ERROR('STOP: The value of input variable TOGETHER is not allowed') 
  end if  
! Incompatible variables :
  if ((InVar%ReadIFC.eq.1).and.(InVar%together.eq.1).and.(InVar%Order.gt.2)) then
    MSG_ERROR('STOP: ReadIFC=1, together=1 and Order=3 or 4 are incompatible')
  end if  

! Compute Nstep as a function of the slice
  InVar%nstep_tot=int(float(InVar%nstep_max-InVar%nstep_min)/float(InVar%Slice)+1)
  write(6,*) 'nstep_tot=',InVar%nstep_tot
  close(40)


 end subroutine tdep_read_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_distrib_data(Hist,InVar,MPIdata)

  implicit none 

  type(Input_Variables_type), intent(inout) :: InVar
  type(MPI_enreg_type), intent(in) :: MPIdata
  type(abihist), intent(in) :: Hist

  integer :: this_istep,istep,iatom,jstep
  double precision :: tmp1,tmp2,tmp3
  
  InVar%my_nstep=MPIdata%my_nstep
 
! Read xred.dat, fcart.dat and etot.dat ASCII files or extract them from the HIST.nc netcdf file.
  write(InVar%stdout,'(a)') ' '
  ABI_MALLOC(InVar%xred,(3,InVar%natom,InVar%my_nstep))  ; InVar%xred(:,:,:)=0.d0
  ABI_MALLOC(InVar%fcart,(3,InVar%natom,InVar%my_nstep)) ; InVar%fcart(:,:,:)=0.d0
  ABI_MALLOC(InVar%etot,(InVar%my_nstep))                ; InVar%etot(:)=0.d0
  this_istep=0
  jstep=0
  if (InVar%netcdf) then
    do istep=InVar%nstep_min,InVar%nstep_max
      if (mod(istep-InVar%nstep_min,InVar%Slice).ne.0) then
        cycle
      else
        jstep=jstep+1
        if (.not.MPIdata%my_step(jstep)) cycle
        this_istep=this_istep+1
        InVar%xred(:,:,this_istep) =Hist%xred (:,:,istep)
        InVar%fcart(:,:,this_istep)=Hist%fcart(:,:,istep)
        InVar%etot(this_istep)     =Hist%etot     (istep)
      end if
    end do !istep  
    write(InVar%stdout,'(a)') ' The Xred, Fcart, Etot and Stress data are extracted from the NetCDF file: HIST.nc '
  else
    open(unit=60,file='fcart.dat')
    open(unit=50,file='xred.dat')
    open(unit=40,file='etot.dat')
    do istep=1,InVar%nstep_min-1
      read(40,*) tmp1
      do iatom=1,InVar%natom
        read(50,*) tmp1,tmp2,tmp3
        read(60,*) tmp1,tmp2,tmp3
      end do
    end do 
    do istep=InVar%nstep_min,InVar%nstep_max
      if (mod(istep-InVar%nstep_min,InVar%Slice).ne.0) then
        read(40,*) tmp1
        do iatom=1,InVar%natom
          read(50,*) tmp1,tmp2,tmp3
          read(60,*) tmp1,tmp2,tmp3
        end do
      else
        jstep=jstep+1
        if (.not.MPIdata%my_step(jstep)) then
          read(40,*) tmp1
          do iatom=1,InVar%natom
            read(50,*) tmp1,tmp2,tmp3
            read(60,*) tmp1,tmp2,tmp3
          end do
        else
          this_istep=this_istep+1
          read(40,*) InVar%etot(this_istep)
          do iatom=1,InVar%natom
            read(50,*) InVar%xred (1,iatom,this_istep),InVar%xred (2,iatom,this_istep),InVar%xred (3,iatom,this_istep)
            read(60,*) InVar%fcart(1,iatom,this_istep),InVar%fcart(2,iatom,this_istep),InVar%fcart(3,iatom,this_istep)
          end do
        end if !my_step  
      end if !slice
    end do !istep
    close(40)
    close(50)
    close(60)
    write(InVar%stdout,'(a)') ' The Xred, Fcart and Etot data are extracted from the ASCII files: xred.dat, fcart.dat \& etot.dat'
  end if !netcdf

 end subroutine tdep_distrib_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine tdep_init_MPIshell(InVar,MPIdata,Shelllllll)
 subroutine tdep_init_MPIshell(InVar,MPIdata)

  implicit none 

  type(Input_Variables_type), intent(in) :: InVar
  type(MPI_enreg_type), intent(in) :: MPIdata

 end subroutine tdep_init_MPIshell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_init_MPIdata(InVar,MPIdata)

  implicit none 

  type(Input_Variables_type), intent(in) :: InVar
  type(MPI_enreg_type), intent(out) :: MPIdata
  integer :: ii,remain,ierr,iproc,istep
  integer, allocatable :: nstep_acc(:)
  integer, allocatable :: tab_step(:)
  character(len=500) :: message
  
#if defined HAVE_MPI
  integer :: dimcart,commcart_2d,me_cart_2d
  logical :: reorder
  integer,allocatable :: coords(:),sizecart(:)
  logical,allocatable :: periode(:), keepdim(:)
#endif

! Check the number of processors
  MPIdata%nproc_shell=InVar%Nproc(1)
  MPIdata%nproc_step =InVar%Nproc(2)
  MPIdata%nproc = xmpi_comm_size(xmpi_world)
  if (MPIdata%nproc_step*MPIdata%nproc_shell.ne.MPIdata%nproc) then
    MSG_WARNING('The parallelization is performed over steps')
    MPIdata%nproc_step = xmpi_comm_size(xmpi_world)
  end if  

  MPIdata%master         = 0
  MPIdata%iam_master     =.false.
! Initialize the MPIdata datastructure for sequential calculation
  if (MPIdata%nproc.eq.1) then
    MPIdata%comm_shell     = xmpi_comm_null
    MPIdata%comm_step      = xmpi_comm_null
    MPIdata%comm_shellstep = xmpi_comm_null
    MPIdata%me_shell       = 0
    MPIdata%me_step        = 0
    MPIdata%my_nstep       = InVar%nstep_tot
    MPIdata%iam_master     = (MPIdata%me_step == MPIdata%master)
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! Parallel calculation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! Definition of the processor grid !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined HAVE_MPI
!FB  if (MPIdata%nproc.eq.1) return

! Create the global cartesian 2D-communicator
  dimcart=2
  ABI_ALLOCATE(sizecart,(dimcart))
  ABI_ALLOCATE(periode,(dimcart))
  sizecart(1)=MPIdata%nproc_shell ! MPIdata%nproc_shell
  sizecart(2)=MPIdata%nproc_step  ! MPIdata%nproc_step
  periode(:)=.false.;reorder=.false.
  call MPI_CART_CREATE(xmpi_world,dimcart,sizecart,periode,reorder,commcart_2d,ierr)
  ABI_DEALLOCATE(periode)
  ABI_DEALLOCATE(sizecart)

! Find the index and coordinates of the current processor
  call MPI_COMM_RANK(commcart_2d,me_cart_2d,ierr)
  ABI_ALLOCATE(coords,(dimcart))
  call MPI_CART_COORDS(commcart_2d,me_cart_2d,dimcart,coords,ierr)
  MPIdata%me_shell=coords(1)
  MPIdata%me_step =coords(2)
  ABI_DEALLOCATE(coords)
  if ((MPIdata%me_shell == MPIdata%master).and.(MPIdata%me_step == MPIdata%master)) then
    MPIdata%iam_master = .true.
  end if  

  ABI_ALLOCATE(keepdim,(dimcart))
! Create the communicator for shell distribution
  keepdim(1)=.true.
  keepdim(2)=.false.
  call MPI_CART_SUB(commcart_2d,keepdim,MPIdata%comm_shell,ierr)
! Create the communicator for step distribution
  keepdim(1)=.false.
  keepdim(2)=.true.
  call MPI_CART_SUB(commcart_2d,keepdim,MPIdata%comm_step,ierr)
! Create the communicator for shellstep distribution
  keepdim(1)=.true.
  keepdim(2)=.true.
  call MPI_CART_SUB(commcart_2d,keepdim,MPIdata%comm_shellstep,ierr)
  ABI_DEALLOCATE(keepdim)
  call xmpi_comm_free(commcart_2d)

! Write some data
  write(message,'(1x,a20,2(1x,i4))') 'Number of processors :',MPIdata%nproc_shell,MPIdata%nproc_step
  call wrtout(InVar%stdout,message,'COLL')
!FB  write(message,'(a,2i5)') 'me_shell and me_step : ',MPIdata%me_shell,MPIdata%me_step
!FB  call wrtout(InVar%stdout,message,'COLL')
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! Distribution over STEP processors !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MPIdata%my_nstep =int(InVar%nstep_tot/MPIdata%nproc_step)
  remain=InVar%nstep_tot-MPIdata%nproc_step*MPIdata%my_nstep
  do ii=1,remain
    if ((ii-1).eq.MPIdata%me_step) MPIdata%my_nstep=MPIdata%my_nstep+1
  end do
  ABI_MALLOC(MPIdata%nstep_all,(MPIdata%nproc_step)); MPIdata%nstep_all(:)=zero
  call xmpi_allgather(MPIdata%my_nstep,MPIdata%nstep_all,MPIdata%comm_step,ierr)
  write(InVar%stdout,'(a)') ' '
  write(InVar%stdout,'(a)') ' WARNING: ALL the quantities are now computed :'
  write(InVar%stdout,'(a,1x,i4)') '                                      from nstep_min=',InVar%nstep_min
  write(InVar%stdout,'(a,1x,i4)') '                                        to nstep_max=',InVar%nstep_max
  if (InVar%Slice.ne.1) then
    write(InVar%stdout,'(a,1x,i4)') '                                    by using a slice=',InVar%Slice
  end if  
  write(InVar%stdout,'(a,1x,i4)') '          So, the real number of time steps is nstep=',InVar%nstep_tot
  if (MPIdata%nproc_step.gt.1) then
    write(Invar%stdout,'(a,1000(1x,i5))') 'Distribution of number of steps wrt the number of processors=',MPIdata%nstep_all(:)
  end if

  ABI_MALLOC(nstep_acc,(MPIdata%nproc_step+1)); nstep_acc(:)=zero
  nstep_acc(1)=0
  do ii=2,MPIdata%nproc_step+1
    nstep_acc(ii)=nstep_acc(ii-1)+MPIdata%nstep_all(ii-1)
  end do
  if (nstep_acc(MPIdata%nproc_step+1).ne.InVar%nstep_tot) then
    write(6,*) 'STOP : pb in nstep_acc'
    stop
  end if

  ABI_MALLOC(tab_step,(InVar%nstep_tot)); tab_step(:)=zero
  ABI_MALLOC(MPIdata%my_step ,(InVar%nstep_tot)); MPIdata%my_step (:)=.false.
  do iproc=1,MPIdata%nproc_step
    do istep=1,InVar%nstep_tot
      if ((istep.gt.nstep_acc(iproc)).and.(istep.le.nstep_acc(iproc+1))) then
        tab_step(istep)=iproc-1
      end if
    end do
  end do
  do istep=1,InVar%nstep_tot
    MPIdata%my_step(istep) = (tab_step(istep) == MPIdata%me_step)
  end do

  ABI_MALLOC(MPIdata%shft_step,(MPIdata%nproc_step)); MPIdata%shft_step(:)=zero
  MPIdata%shft_step(1)=0
  do ii=2,MPIdata%nproc_step
    MPIdata%shft_step(ii)=MPIdata%shft_step(ii-1)+MPIdata%nstep_all(ii-1)
  end do
  ABI_FREE(nstep_acc)
  ABI_FREE(tab_step)

 end subroutine tdep_init_MPIdata

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_clean_mpi(MPIdata)

  implicit none 

  type(MPI_enreg_type), intent(inout) :: MPIdata

  ABI_FREE(MPIdata%shft_step)
  ABI_FREE(MPIdata%nstep_all)
  ABI_FREE(MPIdata%my_step)

 end subroutine tdep_clean_mpi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module m_tdep_readwrite
