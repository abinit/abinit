
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
    integer :: use_ideal_positions
    integer :: stdout
    integer :: stdlog
    integer :: bzpath
    integer :: order
    integer :: slice
    integer :: enunit
    integer :: readifc
    integer :: together
    integer :: nproc(2)
    integer :: bzlength
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
    double precision :: rcut
    double precision :: rcut3
    double precision :: rcut4
    double precision :: temperature
    double precision :: tolread
    double precision :: tolinbox
    double precision :: tolmatch
    double precision :: tolmotif
    double precision :: rprimd_md(3,3)
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
    character (len=200) :: input_prefix
    
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
 subroutine tdep_print_Aknowledgments(Invar)

  type(Input_Variables_type) :: Invar
  integer :: stdout
  stdout = Invar%stdout

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
  write(stdout,'(a)') ' [1] a-TDEP: Temperature Dependent Effective Potential for Abinit '
  write(stdout,'(a)') ' -- Lattice dynamic properties including anharmonicity' 
  write(stdout,'(a)') ' F. Bottin, J. Bieder and J. Bouchet, Comput. Phys. Comm. 254, 107301 (2020).' ! [[cite:Bottin2020]]
  write(stdout,'(a)') ' Strong suggestion to cite this paper in your publications.'
  write(stdout,'(a)') ' '
  write(stdout,'(a)') ' [2] Thermal evolution of vibrational properties of alpha-U' 
  write(stdout,'(a)') ' J. Bouchet and F. Bottin, Phys. Rev. B 92, 174108 (2015).' ! [[cite:Bouchet2015]]
  write(stdout,'(a)') ' Strong suggestion to cite this paper in your publications.'
  write(stdout,'(a)') ' '
  write(stdout,'(a)') ' [3] Lattice dynamics of anharmonic solids from first principles'
  write(stdout,'(a)') ' O. Hellman, I.A. Abrikosov and S.I. Simak, Phys. Rev. B 84, 180301(R) (2011).' ! [[cite:Hellman2011]]
  write(stdout,'(a)') ' '
  write(stdout,'(a)') ' [4] Temperature dependent effective potential method for accurate free energy calculations of solids'
  write(stdout,'(a)') ' O. Hellman, P. Steneteg, I.A. Abrikosov and S.I. Simak, Phys. Rev. B 87, 104111 (2013).' ! [[cite:Hellman2013]]

 end subroutine tdep_print_Aknowledgments 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_read_input(Hist,Invar)

#if defined HAVE_NETCDF
 use netcdf
#endif

  type(Input_Variables_type),intent(out) :: Invar
  type(abihist), intent(out) :: Hist

  integer :: values(8)  
  integer :: ncid, ncerr,ierr
  integer :: nimage, mdtime, natom_id,nimage_id,time_id,xyz_id,six_id
  integer :: ntypat_id
  integer :: ii,jj,tmp,istep,iatom,this_istep
  double precision :: version_value,tmp1,tmp2,tmp3,dtion
  character (len=30):: string,NormalMode,DebugMode,use_ideal_positions
  character (len=30):: born_charge,dielec_constant,tolmotifinboxmatch,TheEnd,bzpath
  character (len=30):: order,slice,enunit,readifc,together,nproc,bzlength,ngqpt1,ngqpt2,dosdeltae
  character (len=8) :: date
  character (len=10) :: time
  character (len=5) :: zone
  character(len=3),parameter :: month_names(12)=(/'Jan','Feb','Mar','Apr','May','Jun',&
&                                                 'Jul','Aug','Sep','Oct','Nov','Dec'/)
  character(len=500) :: filename,inputfilename,message
  logical :: has_nimage
  real(dp), allocatable :: znucl(:)

! Define output files  
  Invar%stdout=8
  Invar%stdlog=6
  !open(unit=Invar%stdlog,file='data.log')

! Define Keywords
  NormalMode='NormalMode'
  DebugMode='DebugMode'
  use_ideal_positions='use_ideal_positions'
  born_charge='born_charge'
  dielec_constant='dielec_constant'
  dosdeltae='dosdeltae'
  bzpath='bzpath'
  bzlength='bzlength'
  order='order'
  slice='slice'
  enunit='enunit'
  readifc='readifc'
  together='together'
  nproc='nproc'
  ngqpt1='ngqpt1'
  ngqpt2='ngqpt2'
  tolmotifinboxmatch='tolmotifinboxmatch'
  TheEnd='TheEnd'
! Define default values
  Invar%angle_alpha=90.d0
  Invar%bzpath=0
  Invar%order=2
  Invar%slice=1
  Invar%enunit=0
  Invar%together=1
  Invar%nproc(:)=1
  Invar%bzlength=0
  Invar%tolread=1.d-8
  Invar%tolmotif=5.d-2
  Invar%tolinbox=5.d-2
  Invar%tolmatch=5.d-2
  Invar%dosdeltae=4.5d-6
  Invar%debug=.false.
  Invar%loto=.false.
  Invar%netcdf=.false.
  Invar%use_ideal_positions=0
  version_value=3.d0
! In order to have an accuracy better than 1meV  
  Invar%ngqpt1(:)=8
  Invar%ngqpt2(:)=32

! Check if a NetCDF file is available
  write(Invar%stdlog,'(a)',err=10) ' Give name for input file '
  read(*, '(a)',err=10) inputfilename
  if ( inputfilename == "" ) inputfilename='input.in'
  write(Invar%stdlog, '(a)',err=10) '.'//trim(inputfilename)
10 continue
  write(Invar%stdlog,'(a)',err=11) ' Give root name for generic input files (NetCDF or ASCII)'
  read(*, '(a)',err=11) Invar%input_prefix
  if ( Invar%input_prefix == "" ) then
    filename='HIST.nc'
  else
    filename=trim(Invar%input_prefix)//'HIST.nc'  
  end if  
  write(Invar%stdlog, '(a)',err=11) '.'//trim(Invar%input_prefix)
11 continue
  write(Invar%stdlog,'(a)', err=12)' Give root name for generic output files:'
  read (*, '(a)', err=12) Invar%output_prefix
  if ( Invar%output_prefix == "" ) then
    if (xmpi_comm_rank(xmpi_world).eq.0) open(unit=Invar%stdout,file='atdep.out')
  else
    if (xmpi_comm_rank(xmpi_world).eq.0) open(unit=Invar%stdout,file=trim(Invar%output_prefix)//'.out')
  end if  
  write (Invar%stdlog, '(a)', err=12 ) '.'//trim(Invar%output_prefix)
12 continue
  if ( inputfilename == "" ) inputfilename='input.in'
  if ( filename == "" ) filename='HIST.nc'

  open(unit=InVar%stdout,file=trim(InVar%output_prefix)//'.abo')

#if defined HAVE_NETCDF
 !Open netCDF file
  ncerr=nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=ncid)
  if(ncerr /= NF90_NOERR) then
    write(Invar%stdout,'(3a)') '-'//'Could not open ',trim(filename),', starting from scratch'
    Invar%netcdf=.false.
  else
    write(Invar%stdout,'(3a)') '-'//'Succesfully open ',trim(filename),' for reading'
    write(Invar%stdout,'(a)') ' Extracting information from NetCDF file...'
    Invar%netcdf=.true.
  end if

  if ( Invar%netcdf) then
    call get_dims_hist(ncid,Invar%natom,Invar%ntypat,nimage,mdtime,&
&       natom_id,ntypat_id,nimage_id,time_id,xyz_id,six_id,has_nimage)
    ABI_MALLOC(Invar%amu,(Invar%ntypat)); Invar%amu(:)=zero
    ABI_MALLOC(Invar%typat,(Invar%natom)); Invar%typat(:)=zero
    ABI_MALLOC(znucl,(Invar%ntypat)) ; znucl(:)=zero
    call read_csts_hist(ncid,dtion,Invar%typat,znucl,Invar%amu)
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
  write(Invar%stdout,*) ' '
  open(unit=40,file=inputfilename)
  read(40,*) string
  if (string.eq.NormalMode) then
    write(Invar%stdout,'(a,f6.1,a)') '.Version ', version_value,' of PHONONS'
  else if (string.eq.DebugMode) then
    Invar%debug=.true.
    write(Invar%stdout,'(a,f6.1,a)') '.Version ', version_value,' of PHONONS (Debug)'
  else
    MSG_ERROR('Please use recent format for the input file')
  end if  
  write(Invar%stdout,'(a)') '.Copyright (C) 1998-2020 ABINIT group (FB,JB).'
  write(Invar%stdout,'(a)') ' ABINIT comes with ABSOLUTELY NO WARRANTY.'
  write(Invar%stdout,'(a)') ' It is free software, and you are welcome to redistribute it'
  write(Invar%stdout,'(a)') ' under certain conditions (GNU General Public License,'
  write(Invar%stdout,'(a)') ' see ~abinit/COPYING or http://www.gnu.org/copyleft/gpl.txt).'
  write(Invar%stdout,*) ' '
  write(Invar%stdout,'(a)') ' ABINIT is a project of the Universite Catholique de Louvain,'
  write(Invar%stdout,'(a)') ' Corning Inc. and other collaborators, see'
  write(Invar%stdout,'(a)') ' ~abinit/doc/developers/contributors.txt .'
  write(Invar%stdout,'(a)') ' Please read https://docs.abinit.org/theory/acknowledgments for suggested'
  write(Invar%stdout,'(a)') ' acknowledgments of the ABINIT effort.'
  write(Invar%stdout,'(a)') ' For more information, see http://www.abinit.org .'

  call date_and_time(date,time,zone,values)
  write(Invar%stdout,'(/,a,i2,1x,a,1x,i4,a)') '.Starting date : ',values(3),month_names(values(2)),values(1),'.'

! Read (and echo) of input variables from the input.in input file
  write(Invar%stdout,*) ' '
  write(Invar%stdout,*) '#############################################################################'
  write(Invar%stdout,*) '######################### ECHO OF INPUT FILE ################################'
  write(Invar%stdout,*) '#############################################################################'
! Define unit cell  
  read(40,*) string
  write(Invar%stdout,'(a)') ' ======================= Define the unitcell =================================' 
  read(40,*) string,Invar%bravais(1),Invar%bravais(2)
  write(Invar%stdout,'(1x,a20,1x,i4,1x,i4)') string,Invar%bravais(1),Invar%bravais(2)
  if ((Invar%bravais(1).eq.2).or.(Invar%bravais(1).eq.5)) then
    read(40,*) string,Invar%angle_alpha
    write(Invar%stdout,'(1x,a20,1x,f15.10)') string,Invar%angle_alpha
  else
    !read(40,*)
    Invar%angle_alpha=90.d0
  end if
  read(40,*) string,Invar%natom_unitcell
  write(Invar%stdout,'(1x,a20,1x,i4)') string,Invar%natom_unitcell
  ABI_MALLOC(Invar%xred_unitcell,(3,Invar%natom_unitcell)); Invar%xred_unitcell(:,:)=zero
  read(40,*) string,Invar%xred_unitcell(:,:)
  write(Invar%stdout,'(1x,a20)') string
  do ii=1,Invar%natom_unitcell
    write(Invar%stdout,'(22x,3(f15.10,1x))') (Invar%xred_unitcell(jj,ii), jj=1,3)
  end do  
  ABI_MALLOC(Invar%typat_unitcell,(Invar%natom_unitcell)); Invar%typat_unitcell(:)=0 
  read(40,*) string,Invar%typat_unitcell(:)
  write(Invar%stdout,'(1x,a20,20(1x,i4))') string,(Invar%typat_unitcell(jj),jj=1,Invar%natom_unitcell)
! To avoid some troubles and inconsistency between the .in and .nc files, when we use NetCDF or not.    
  if (Invar%netcdf) then
    read(40,*) string
    backspace(40)
    if (string.eq.'ntypat') then
      write(Invar%stdlog,'(1x,a)') 'When the NetCDF file .nc is used, the ntypat keywork is not allowed.' 
      MSG_ERROR('ACTION : Please modify your input file')
    end if  
    string='ntypat'
  else
    read(40,*) string
    backspace(40)
    if (string.ne.'ntypat') then
      write(Invar%stdlog,'(1x,a)') 'The NetCDF file .nc is not used.' 
      write(Invar%stdlog,'(1x,a,1x,a)') 'In your input file, the code search the ntypat keywork but found :',string
      MSG_ERROR('ACTION : Please modify your input file')
    end if  
    read(40,*) string,Invar%ntypat
  end if  
  write(Invar%stdout,'(1x,a20,1x,i4)') string,Invar%ntypat
  if (Invar%netcdf) then
    string='amu'
  else  
    ABI_MALLOC(Invar%amu,(Invar%ntypat)); Invar%amu(:)=zero
    read(40,*) string,Invar%amu(:)
  end if  
  write(Invar%stdout,'(1x,a20,20(1x,f15.10))') string,(Invar%amu(jj),jj=1,Invar%ntypat)
! Define supercell (as a function of the unitcell defined above)
  read(40,*) string
  write(Invar%stdout,'(a)') ' ======================= Define the supercell ================================' 
  if (Invar%netcdf) then
    Invar%rprimd_md(:,:)=Hist%rprimd(:,:,Hist%ihist)
    string='rprimd'
  else
    read(40,*) string,Invar%rprimd_md(1,:),Invar%rprimd_md(2,:),Invar%rprimd_md(3,:)
  end if  
  write(Invar%stdout,'(1x,a20)') string
  do ii=1,3
    write(Invar%stdout,'(22x,3(f15.10,1x))') (Invar%rprimd_md(ii,jj),jj=1,3)
  end do  
  read(40,*) string,Invar%multiplicity(1,:),Invar%multiplicity(2,:),Invar%multiplicity(3,:)
  write(Invar%stdout,'(1x,a20)') string
  do ii=1,3
    write(Invar%stdout,'(22x,3(f15.10,1x))') (Invar%multiplicity(ii,jj),jj=1,3)
  end do  
  if (Invar%netcdf) then
    string='natom'
  else
    read(40,*) string,Invar%natom
  end if
  write(Invar%stdout,'(1x,a20,1x,i4)') string,Invar%natom
  if (Invar%netcdf) then
    string='typat'
  else
    ABI_MALLOC(Invar%typat,(Invar%natom)); Invar%typat(:)=0 
    read(40,*) string,Invar%typat(:)
  end if  
  write(Invar%stdout,'(1x,a20)') string
  do ii=1,Invar%natom,10
    if (ii+9.lt.Invar%natom) then
      write(Invar%stdout,'(22x,10(i4,1x))') (Invar%typat(ii+jj-1),jj=1,10)
    else
      write(Invar%stdout,'(22x,10(i4,1x))') (Invar%typat(jj),jj=ii,Invar%natom)
    end if  
  end do  
  read(40,*) string,Invar%temperature
  write(Invar%stdout,'(1x,a20,1x,f15.10)') string,Invar%temperature
! Define phonons computational details
  read(40,*) string
  write(Invar%stdout,'(a)') ' ======================= Define computational details ========================' 
  read(40,*) string,Invar%nstep_max
  write(Invar%stdout,'(1x,a20,1x,i5)') string,Invar%nstep_max
  read(40,*) string,Invar%nstep_min
  write(Invar%stdout,'(1x,a20,1x,i5)') string,Invar%nstep_min
  read(40,*) string,Invar%rcut
  write(Invar%stdout,'(1x,a20,1x,f15.10)') string,Invar%rcut
! Optional input variables  
  read(40,*) string
  write(Invar%stdout,'(a)') ' ======================= Optional input variables ============================' 
  do ii=1,100
    read(40,*) string
    backspace(40)
    if (string.eq.dosdeltae) then
      read(40,*) string,Invar%dosdeltae
      write(Invar%stdout,'(1x,a20,1x,f15.10)') string,Invar%dosdeltae
    else if (string.eq.use_ideal_positions) then  
      read(40,*) string,Invar%use_ideal_positions
      write(Invar%stdout,'(1x,a20,1x,i4)') string,Invar%use_ideal_positions
    else if (string.eq.born_charge) then  
      ABI_MALLOC(Invar%born_charge,(Invar%ntypat)); Invar%born_charge(:)=0.d0
      Invar%loto=.true.
      read(40,*) string,Invar%born_charge(:)
      write(Invar%stdout,'(1x,a20,20(1x,f15.10))') string,(Invar%born_charge(jj),jj=1,Invar%ntypat)
    else if (string.eq.dielec_constant) then  
      read(40,*) string,Invar%dielec_constant
      write(Invar%stdout,'(1x,a20,1x,f15.10)') string,Invar%dielec_constant
    else if (string.eq.bzpath) then  
      read(40,*) string,Invar%bzpath
      write(Invar%stdout,'(1x,a20,1x,i4)') string,Invar%bzpath
      if (Invar%bzpath.lt.0) then
        ABI_MALLOC(Invar%qpt,(3,abs(Invar%bzpath))); Invar%qpt(:,:)=zero
        write(Invar%stdout,'(a)') ' Q points as given in the input file:'
        do jj=1,abs(Invar%bzpath)
          read(40,*) Invar%qpt(:,jj)
          write(Invar%stdout,'(22x,3(f15.10,1x))') Invar%qpt(:,jj)
        end do  
      else if (Invar%bzpath.gt.0) then
        ABI_MALLOC(Invar%special_qpt,(Invar%bzpath))
        backspace(40)
        read(40,*) string,tmp,(Invar%special_qpt(jj),jj=1,Invar%bzpath)
        write(Invar%stdout,'(a,1x,a2,10("-",a2))') ' Special q-points: ',Invar%special_qpt(:)
      end if
    else if (string.eq.order) then  
      read(40,*) string,Invar%order
      backspace(40)
      if (Invar%order.eq.3) then
        read(40,*) string,Invar%order,Invar%rcut3
        write(Invar%stdout,'(1x,a20,1x,i4,1x,f15.10)') string,Invar%order,Invar%rcut3
        if (Invar%rcut3.gt.Invar%rcut) then
          MSG_ERROR('The cutoff radius of the third order cannot be greater than the second order one.')
        end if  
      else if (Invar%order.eq.4) then
        read(40,*) string,Invar%order,Invar%rcut3,Invar%rcut4
        write(Invar%stdout,'(1x,a20,1x,i4,2(1x,f15.10))') string,Invar%order,Invar%rcut3,Invar%rcut4
        if (Invar%rcut4.gt.Invar%rcut) then
          MSG_ERROR('The cutoff radius of the fourth order cannot be greater than the second order one.')
        end if  
      else
        MSG_ERROR('Only the 3rd and 4th orders are allowed. Change your input file.')
      end if
    else if (string.eq.slice) then  
      read(40,*) string,Invar%slice
      write(Invar%stdout,'(1x,a20,1x,i4)') string,Invar%slice
    else if (string.eq.enunit) then  
      read(40,*) string,Invar%enunit
      if (Invar%enunit.eq.0) write(Invar%stdout,'(1x,a20,1x,i4,1x,a)') string,Invar%enunit,'(Phonon frequencies in meV)'
      if (Invar%enunit.eq.1) write(Invar%stdout,'(1x,a20,1x,i4,1x,a)') string,Invar%enunit,'(Phonon frequencies in cm-1)'
      if (Invar%enunit.eq.2) write(Invar%stdout,'(1x,a20,1x,i4,1x,a)') string,Invar%enunit,'(Phonon frequencies in Ha)'
    else if (string.eq.nproc) then  
      read(40,*) string,Invar%nproc(1),Invar%nproc(2)
      write(Invar%stdout,'(1x,a20,1x,i4,1x,i4)') string,Invar%nproc(1),Invar%nproc(2)
    else if (string.eq.readifc) then  
      read(40,*) string,Invar%readifc
      if (Invar%readifc.eq.1) then
        backspace(40)
        read(40,*) string,Invar%readifc,Invar%tolread
        write(Invar%stdout,'(1x,a20,1x,i4,1x,f15.10)') string,Invar%readifc,Invar%tolread
      else  
        write(Invar%stdout,'(1x,a20,1x,i4)') string,Invar%readifc
      end if  
    else if (string.eq.together) then  
      read(40,*) string,Invar%together
      write(Invar%stdout,'(1x,a20,1x,i4)') string,Invar%together
    else if (string.eq.bzlength) then  
      read(40,*) string,Invar%bzlength
      write(Invar%stdout,'(1x,a20,1x,i4)') string,Invar%bzlength
      ABI_MALLOC(Invar%lgth_segments,(Invar%bzlength))
      backspace(40)
      read(40,*) string,tmp,(Invar%lgth_segments(jj),jj=1,Invar%bzlength)
      write(Invar%stdout,'(a,1x,i3,10("-",i3))') ' Length of BZ : ',Invar%lgth_segments(:)
    else if (string.eq.ngqpt1) then  
      read(40,*) string,Invar%ngqpt1(:)
      write(Invar%stdout,'(1x,a20,1x,3(i4,1x))') string,Invar%ngqpt1(:)
    else if (string.eq.ngqpt2) then  
      read(40,*) string,Invar%ngqpt2(:)
      write(Invar%stdout,'(1x,a20,1x,3(i4,1x))') string,Invar%ngqpt2(:)
    else if (string.eq.tolmotifinboxmatch) then  
      read(40,*) string,Invar%tolmotif,Invar%tolinbox,Invar%tolmatch
      write(Invar%stdout,'(1x,a20,f10.5)') 'tolmotif            ',Invar%tolmotif
      write(Invar%stdout,'(1x,a20,f10.5)') 'tolinbox            ',Invar%tolinbox
      write(Invar%stdout,'(1x,a20,f10.5)') 'tolmatch            ',Invar%tolmatch
    else if (string.eq.TheEnd) then
      exit
    else 
      write(Invar%stdout,'(a,1x,a)') 'This keyword is not allowed',string
      MSG_ERROR('A keyword is not allowed. See the log file.')
    end if  
  end do
! Output very important informations 
  if (Invar%use_ideal_positions.eq.0) then
    write(Invar%stdout,'(a)') ' USE AVERAGE POSITIONS TO COMPUTE SPECTRUM'
  else if (Invar%use_ideal_positions.eq.1) then
    write(Invar%stdout,'(a)') ' USE IDEAL POSITIONS TO COMPUTE SPECTRUM'
  else
    write(Invar%stdout,'(a)') ' STOP: THIS VALUE IS NOT ALLOWED FOR use_ideal_positions'
  end if
  if (Invar%loto) write(Invar%stdout,'(a)') ' USE NON-ANALYTICAL CORRECTIONS (LO-TO)'

! Allowed values
  if ((Invar%together.ne.1).and.(Invar%together.ne.0)) then
    MSG_ERROR('STOP: The value of input variable TOGETHER is not allowed') 
  end if  
! Incompatible variables :
  if ((Invar%readifc.eq.1).and.(Invar%together.eq.1).and.(Invar%order.gt.2)) then
    MSG_ERROR('STOP: readifc=1, together=1 and order=3 or 4 are incompatible')
  end if  

! Compute Nstep as a function of the slice
  Invar%nstep_tot=int(float(Invar%nstep_max-Invar%nstep_min)/float(Invar%slice)+1)
  close(40)

  write(6,*) 'nstep_tot=',Invar%nstep_tot


 end subroutine tdep_read_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_distrib_data(Hist,Invar,MPIdata)

  implicit none 

  type(Input_Variables_type), intent(inout) :: Invar
  type(MPI_enreg_type), intent(in) :: MPIdata
  type(abihist), intent(in) :: Hist

  integer :: this_istep,istep,iatom,jstep
  double precision :: tmp1,tmp2,tmp3
  
  Invar%my_nstep=MPIdata%my_nstep
 
! Read xred.dat, fcart.dat and etot.dat ASCII files or extract them from the HIST.nc netcdf file.
  write(Invar%stdout,'(a)') ' '
  ABI_MALLOC(Invar%xred,(3,Invar%natom,Invar%my_nstep))  ; Invar%xred(:,:,:)=0.d0
  ABI_MALLOC(Invar%fcart,(3,Invar%natom,Invar%my_nstep)) ; Invar%fcart(:,:,:)=0.d0
  ABI_MALLOC(Invar%etot,(Invar%my_nstep))                ; Invar%etot(:)=0.d0
  this_istep=0
  jstep=0
  if (Invar%netcdf) then
    do istep=Invar%nstep_min,Invar%nstep_max
      if (mod(istep-Invar%nstep_min,Invar%slice).ne.0) then
        cycle
      else
        jstep=jstep+1
        if (.not.MPIdata%my_step(jstep)) cycle
        this_istep=this_istep+1
        Invar%xred(:,:,this_istep) =Hist%xred (:,:,istep)
        Invar%fcart(:,:,this_istep)=Hist%fcart(:,:,istep)
        Invar%etot(this_istep)     =Hist%etot     (istep)
      end if
    end do !istep  
    write(Invar%stdout,'(a)') ' The Xred, Fcart, Etot and Stress data are extracted from the NetCDF file: HIST.nc '
  else
    open(unit=60,file=trim(Invar%input_prefix)//'fcart.dat')
    open(unit=50,file=trim(Invar%input_prefix)//'xred.dat')
    open(unit=40,file=trim(Invar%input_prefix)//'etot.dat')
    do istep=1,Invar%nstep_min-1
      read(40,*) tmp1
      do iatom=1,Invar%natom
        read(50,*) tmp1,tmp2,tmp3
        read(60,*) tmp1,tmp2,tmp3
      end do
    end do 
    do istep=Invar%nstep_min,Invar%nstep_max
      if (mod(istep-Invar%nstep_min,Invar%slice).ne.0) then
        read(40,*) tmp1
        do iatom=1,Invar%natom
          read(50,*) tmp1,tmp2,tmp3
          read(60,*) tmp1,tmp2,tmp3
        end do
      else
        jstep=jstep+1
        if (.not.MPIdata%my_step(jstep)) then
          read(40,*) tmp1
          do iatom=1,Invar%natom
            read(50,*) tmp1,tmp2,tmp3
            read(60,*) tmp1,tmp2,tmp3
          end do
        else
          this_istep=this_istep+1
          read(40,*) Invar%etot(this_istep)
          do iatom=1,Invar%natom
            read(50,*) Invar%xred (1,iatom,this_istep),Invar%xred (2,iatom,this_istep),Invar%xred (3,iatom,this_istep)
            read(60,*) Invar%fcart(1,iatom,this_istep),Invar%fcart(2,iatom,this_istep),Invar%fcart(3,iatom,this_istep)
          end do
        end if !my_step  
      end if !slice
    end do !istep
    close(40)
    close(50)
    close(60)
    write(Invar%stdout,'(a)') ' The Xred, Fcart and Etot data are extracted from the ASCII files: xred.dat, fcart.dat \& etot.dat'
  end if !netcdf

 end subroutine tdep_distrib_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine tdep_init_MPIshell(Invar,MPIdata,Shelllllll)
 subroutine tdep_init_MPIshell(Invar,MPIdata)

  implicit none 

  type(Input_Variables_type), intent(in) :: Invar
  type(MPI_enreg_type), intent(in) :: MPIdata

 end subroutine tdep_init_MPIshell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_init_MPIdata(Invar,MPIdata)

  implicit none 

  type(Input_Variables_type), intent(in) :: Invar
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
  MPIdata%nproc_shell=Invar%nproc(1)
  MPIdata%nproc_step =Invar%nproc(2)
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
    MPIdata%my_nstep       = Invar%nstep_tot
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
  call wrtout(Invar%stdout,message,'COLL')
!FB  write(message,'(a,2i5)') 'me_shell and me_step : ',MPIdata%me_shell,MPIdata%me_step
!FB  call wrtout(Invar%stdout,message,'COLL')
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! Distribution over STEP processors !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MPIdata%my_nstep =int(Invar%nstep_tot/MPIdata%nproc_step)
  remain=Invar%nstep_tot-MPIdata%nproc_step*MPIdata%my_nstep
  do ii=1,remain
    if ((ii-1).eq.MPIdata%me_step) MPIdata%my_nstep=MPIdata%my_nstep+1
  end do
  ABI_MALLOC(MPIdata%nstep_all,(MPIdata%nproc_step)); MPIdata%nstep_all(:)=zero
  call xmpi_allgather(MPIdata%my_nstep,MPIdata%nstep_all,MPIdata%comm_step,ierr)
  write(Invar%stdout,'(a)') ' '
  write(Invar%stdout,'(a)') ' WARNING: ALL the quantities are now computed :'
  write(Invar%stdout,'(a,1x,i4)') '                                      from nstep_min=',Invar%nstep_min
  write(Invar%stdout,'(a,1x,i4)') '                                        to nstep_max=',Invar%nstep_max
  if (Invar%slice.ne.1) then
    write(Invar%stdout,'(a,1x,i4)') '                                    by using a slice=',Invar%slice
  end if  
  write(Invar%stdout,'(a,1x,i4)') '          So, the real number of time steps is nstep=',Invar%nstep_tot
  if (MPIdata%nproc_step.gt.1) then
    write(Invar%stdout,'(a,1000(1x,i5))') 'Distribution of number of steps wrt the number of processors=',MPIdata%nstep_all(:)
  end if

  ABI_MALLOC(nstep_acc,(MPIdata%nproc_step+1)); nstep_acc(:)=zero
  nstep_acc(1)=0
  do ii=2,MPIdata%nproc_step+1
    nstep_acc(ii)=nstep_acc(ii-1)+MPIdata%nstep_all(ii-1)
  end do
  if (nstep_acc(MPIdata%nproc_step+1).ne.Invar%nstep_tot) then
    write(6,*) 'STOP : pb in nstep_acc'
    stop
  end if

  ABI_MALLOC(tab_step,(Invar%nstep_tot)); tab_step(:)=zero
  ABI_MALLOC(MPIdata%my_step ,(Invar%nstep_tot)); MPIdata%my_step (:)=.false.
  do iproc=1,MPIdata%nproc_step
    do istep=1,Invar%nstep_tot
      if ((istep.gt.nstep_acc(iproc)).and.(istep.le.nstep_acc(iproc+1))) then
        tab_step(istep)=iproc-1
      end if
    end do
  end do
  do istep=1,Invar%nstep_tot
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
