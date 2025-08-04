
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
  use m_parser
  use m_fstrings,  only : inupper,ljust,next_token
  use m_abimover, only : abimover

 implicit none

  character(len=6),public,parameter :: version_string = '   4.0'

  type Input_type

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
    integer :: alloy
    integer :: ityp_alloy1
    integer :: ityp_alloy2
    integer :: nproc(2)
    integer :: bzlength
    integer :: ngqpt1(3)
    integer :: ngqpt2(3)
    integer :: bravais(11)
    integer :: use_weights
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
    double precision, allocatable :: znucl(:)
    double precision, allocatable :: qpt(:,:)
    double precision, allocatable :: xred_ideal(:,:)
    double precision, allocatable :: xred_unitcell(:,:)
    double precision, allocatable :: xred(:,:,:)
    double precision, allocatable :: fcart(:,:,:)
    double precision, allocatable :: etot(:)
    double precision, allocatable :: weights(:)
    character (len=2), allocatable :: special_qpt(:)
    character (len=fnlen) :: output_prefix
    character (len=fnlen) :: input_prefix
    character (len=fnlen) :: output_file

  end type Input_type

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
    integer, allocatable :: nstep_acc(:)
    integer, allocatable :: nstep_all(:)
    integer, allocatable :: shft_step(:)
    logical, allocatable :: my_shell(:)
    logical, allocatable :: my_step(:)

  end type MPI_enreg_type

 public :: tdep_print_Aknowledgments
 public :: tdep_read_input
 public :: tdep_distrib_data
 public :: tdep_init_MPIdata
!FB public :: tdep_init_MPIshell
 public :: tdep_destroy_mpidata
 public :: tdep_destroy_invar

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_print_Aknowledgments(Invar)

  type(Input_type) :: Invar
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

!----------------------------------------------------------------------

 subroutine tdep_read_input(input_path,Hist,Invar)

#if defined HAVE_NETCDF
 use netcdf
#endif

! Arguments-------------------------------
  character(len=*), intent(in):: input_path
  type(Input_type),intent(out) :: Invar
  type(abihist), intent(out) :: Hist

! Local variables-------------------------
! scalars
  integer :: values(8)
  integer :: ncid, ncerr, me,ierr,master
  integer :: nimage, mdtime, natom_id,nimage_id,time_id,xyz_id,six_id
  integer :: ntypat_id,iatcell
  integer :: ii,jj,shift,iatom,itypat,sum_alloy1,sum_alloy2
  integer:: lenstr, marr, jdtset, tread
  logical :: has_nimage
  double precision :: dtion,amu_average,born_average
  character (len=8) :: date
  character (len=10) :: time
  character (len=5) :: zone
  character(len=500) :: msg
  character(len=500) :: ncfilename,inputfilename
  character(len=strlen):: string, raw_string
! arrays
  character(len=3),parameter :: month_names(12)=(/'Jan','Feb','Mar','Apr','May','Jun',&
&                                                 'Jul','Aug','Sep','Oct','Nov','Dec'/)
  integer, allocatable:: intarr(:)
  integer, allocatable :: typat_unitcell_tmp(:)
  real(dp), allocatable :: xred_unitcell_tmp(:,:),amu_tmp(:),born_charge_tmp(:),znucl_tmp(:)
  real(dp), allocatable:: dprarr(:)

! *********************************************************************

! Define output files
  Invar%stdout=ab_out
  Invar%stdlog=std_out

! Define default values
  Invar%angle_alpha=90.d0
  Invar%bzpath=0
  Invar%order=2
  Invar%slice=1
  Invar%enunit=0
  Invar%together=1
  Invar%alloy=0
  Invar%ityp_alloy1=0
  Invar%ityp_alloy2=0
  Invar%nproc(:)=1
  Invar%bzlength=0
  Invar%tolread=1.d-8
  Invar%tolmotif=5.d-2
  Invar%tolinbox=5.d-2
  Invar%tolmatch=5.d-2
  Invar%dosdeltae=0.2_dp/Ha_cmm1
  Invar%debug=.false.
  Invar%loto=.false.
  Invar%netcdf=.false.
  Invar%use_ideal_positions=0
  Invar%use_weights=0
! In order to have an accuracy better than 1meV
  Invar%ngqpt1(:)=8
  Invar%ngqpt2(:)=32

  me = xmpi_comm_rank(xmpi_world)
  if (me==0) then

    if (len_trim(input_path) == 0) then

      write(std_out, "(2a)")" DeprecationWarning: ",ch10
      write(std_out, "(a)") "     The files file has been deprecated in Abinit10 and will be removed in Abinit11."
      write(std_out, "(2a)")"     Use the syntax `atdep t01.abi` where t01.abi is an atdep input,",ch10
      write(std_out, "(2a)")"     and use input variables output_file, indata_prefix, outdata_prefix.",ch10

      write(Invar%stdlog,'(a)',err=10) ' Give name for input file '
      read(*, '(a)',err=10) inputfilename
      if ( inputfilename == "" ) inputfilename='input.in'
      write(Invar%stdlog, '(a)',err=10) '.'//trim(inputfilename)
10     continue
!     Check if a NetCDF file is available
      write(Invar%stdlog,'(a)',err=11) ' Give root name for generic input files (NetCDF or ASCII)'
      read(*, '(a)',err=11) Invar%input_prefix
      if ( Invar%input_prefix == "" ) then
        ncfilename='HIST.nc'
      else
        ncfilename=trim(Invar%input_prefix)//'_HIST.nc'
      end if
      write(Invar%stdlog, '(a)',err=11) '.'//trim(Invar%input_prefix)
11     continue
      write(Invar%stdlog,'(a)', err=12)' Give root name for generic output files:'
      read (*, '(a)', err=12) Invar%output_prefix
      if ( Invar%output_prefix == "" ) Invar%output_prefix = 'atdep'
      write (Invar%stdlog, '(a)', err=12 ) '.'//trim(Invar%output_prefix)
12     continue
      Invar%output_file = trim(Invar%output_prefix)//'.abo'
    else
      inputfilename = input_path

      ! Read input
      string = repeat(" ", strlen)
      raw_string = repeat(" ", strlen)
      call instrng(input_path, lenstr, 1, strlen, string, raw_string)
      ! To make case-insensitive, map characters to upper case.
      call inupper(string(1:lenstr))

      marr = 3
      ABI_MALLOC(intarr, (marr))
      ABI_MALLOC(dprarr, (marr))
      jdtset = 0

      call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), &
                  "indata_prefix", tread, 'KEY', key_value=Invar%input_prefix)
      if (tread == 0) then
        Invar%input_prefix = ''
        ncfilename='HIST.nc'
      else
        ncfilename=trim(Invar%input_prefix)//'_HIST.nc'
      end if
      write(Invar%stdlog, "(2a)")"- Root name for input files: ", trim(Invar%input_prefix)

      call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), &
                  "outdata_prefix", tread, 'KEY', key_value=Invar%output_prefix)
      if (tread == 0) then
        Invar%output_prefix = 'atdep'
      end if
      write(Invar%stdlog, "(2a)")"- Root name for output files: ", trim(Invar%output_prefix)

      call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), &
                  "output_file", tread, 'KEY', key_value=Invar%output_file)
      if (tread == 0) then
        Invar%output_file = trim(Invar%output_prefix) // '.abo'
      end if
      write(Invar%stdlog, "(2a)")"- Main output file: ", trim(Invar%output_file)

      ABI_FREE(intarr)
      ABI_FREE(dprarr)

    end if

    open(unit=Invar%stdout,file=trim(Invar%output_file))

  end if !me

  master = 0
  call xmpi_bcast(inputfilename,master,xmpi_world,ierr)
  call xmpi_bcast(ncfilename,master,xmpi_world,ierr)
  call xmpi_bcast(Invar%output_prefix,master,xmpi_world,ierr)
  call xmpi_bcast(Invar%input_prefix,master,xmpi_world,ierr)

#if defined HAVE_NETCDF
 !Open netCDF file
  ncerr=nf90_open(path=trim(ncfilename),mode=NF90_NOWRITE,ncid=ncid)
  if(ncerr /= NF90_NOERR) then
    write(Invar%stdlog,'(3a)') '-'//'Could not open ',trim(ncfilename),', starting from scratch'
    Invar%netcdf=.false.
  else
    write(Invar%stdlog,'(3a)') '-'//'Succesfully open ',trim(ncfilename),' for reading'
    write(Invar%stdlog,'(a)') ' Extracting information from NetCDF file...'
    Invar%netcdf=.true.
  end if

  if ( Invar%netcdf) then
    call get_dims_hist(ncid,Invar%natom,Invar%ntypat,nimage,mdtime,&
&       natom_id,ntypat_id,nimage_id,time_id,xyz_id,six_id,has_nimage)
    ABI_MALLOC(Invar%amu,(Invar%ntypat)); Invar%amu(:)=zero
    ABI_MALLOC(Invar%typat,(Invar%natom)); Invar%typat(:)=zero
    ABI_MALLOC(Invar%znucl,(Invar%ntypat)) ; Invar%znucl(:)=zero
    call read_csts_hist(ncid,dtion,Invar%typat,Invar%znucl,Invar%amu)

    ! Need to close NetCDF file because it is going to be reopened by read_md_hist
    ncerr = nf90_close(ncid)
    ! .false. -> Velocities note used
    ! .true. -> acell and rprimd may change (2017_04 only NVT/isoK used but maybe
    ! .false. -> read all times
    ! NPT one day ?)
    call read_md_hist(ncfilename,Hist,.false.,.true.,.false.)
  end if
#endif

! =========================================================================== !
! Read input file

  string = repeat(" ", strlen)
  raw_string = repeat(" ", strlen)
  call instrng(input_path, lenstr, 1, strlen, string, raw_string)
  ! To make case-insensitive, map characters to upper case.
  call inupper(string(1:lenstr))

! marr is the aximum array size. It is thus a hard-coded maximum value
! for 3 * (number of atoms)
  marr = 9000
  ABI_MALLOC(intarr, (marr))
  ABI_MALLOC(dprarr, (marr))
  jdtset = 0

! Mandatory input variables
! -------------------------

! Bravais lattice
  call intagm(dprarr, intarr, jdtset, marr, 2, string(1:lenstr), 'brav', tread, 'INT')
  if (tread == 0) then
    write(msg,*)&
     'Variable "brav" is mandatory, but was not found in input file.'
    ABI_ERROR(msg)
  end if
  Invar%bravais(1:2) = intarr(1:2)

! Angle alpha
  if ((Invar%bravais(1).eq.2).or.(Invar%bravais(1).eq.5)) then
    call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'angle', tread, 'DPR')
    if (tread == 0) then
      write(msg,*)&
       'Variable "angle" is mandatory for this bravais lattice,',ch10,&
       'but was not found in input file.'
      ABI_ERROR(msg)
    end if
    Invar%angle_alpha = dprarr(1)
  else
    Invar%angle_alpha=90.d0
  end if

! natom_unitcell
  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'natom_unitcell', tread, 'INT')
  if (tread == 0) then
    write(msg,*)&
     'Variable "natom_unitcell" is mandatory, but was not found in input file.'
    ABI_ERROR(msg)
  end if
  Invar%natom_unitcell = intarr(1)
  if (3*Invar%natom_unitcell.gt.marr) then
    write(msg,*)&
     'Maximum number of atoms exceeded. Modify source code to circumvent problem.'
    ABI_ERROR(msg)
  end if

! xred_unitcell
  ABI_MALLOC(Invar%xred_unitcell,(3,Invar%natom_unitcell)); Invar%xred_unitcell(:,:)=zero
  call intagm(dprarr, intarr, jdtset, marr, 3*Invar%natom_unitcell, string(1:lenstr), 'xred_unitcell', tread, 'DPR')
  if (tread == 0) then
    write(msg,*)&
     'Variable "xred_unitcell" is mandatory, but was not found in input file.'
    ABI_ERROR(msg)
  end if
  Invar%xred_unitcell(:,:) = reshape(dprarr(1:3*Invar%natom_unitcell),(/3,Invar%natom_unitcell/))
  do ii=1,3
    do iatcell=1,Invar%natom_unitcell
      if ((Invar%xred_unitcell(ii,iatcell).le.(-0.5)).or.(Invar%xred_unitcell(ii,iatcell).gt.(0.5))) then
        do while (Invar%xred_unitcell(ii,iatcell).le.(-0.5))
          Invar%xred_unitcell(ii,iatcell)=Invar%xred_unitcell(ii,iatcell)+1.d0
        end do
        do while (Invar%xred_unitcell(ii,iatcell).gt.(0.5))
          Invar%xred_unitcell(ii,iatcell)=Invar%xred_unitcell(ii,iatcell)-1.d0
        end do
      end if
    end do
  end do

! typat_unitcell
  ABI_MALLOC(Invar%typat_unitcell,(Invar%natom_unitcell)); Invar%typat_unitcell(:)=0
  call intagm(dprarr, intarr, jdtset, marr, Invar%natom_unitcell, string(1:lenstr), 'typat_unitcell', tread, 'INT')
  if (tread == 0) then
    write(msg,*)&
     'Variable "typat_unitcell" is mandatory, but was not found in input file.'
    ABI_ERROR(msg)
  end if
  Invar%typat_unitcell(:) = intarr(1:Invar%natom_unitcell)

  if (Invar%netcdf) then
    Invar%rprimd_md(:,:)=TRANSPOSE(Hist%rprimd(:,:,Hist%ihist))
  else
! ntypat
    call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'ntypat', tread, 'INT')
    if (tread == 1) then
      Invar%ntypat = intarr(1)
    else
      write(msg,*)&
       'The NetCDF file .nc is not used.',ch10,&
       'The variable "ntypat" is thus mandatory.',ch10,&
       'ACTION : Please modify your input file'
      ABI_ERROR(msg)
    end if

! amu
    ABI_MALLOC(Invar%amu,(Invar%ntypat)); Invar%amu(:)=zero
    call intagm(dprarr, intarr, jdtset, marr, Invar%ntypat, string(1:lenstr), 'amu', tread, 'DPR')
    if (tread == 1) then
      Invar%amu = dprarr(1:Invar%ntypat)
    else
      write(msg,*)&
       'The NetCDF file .nc is not used.',ch10,&
       'The variable "amu" is thus mandatory.',ch10,&
       'ACTION : Please modify your input file'
      ABI_ERROR(msg)
    end if

! rprimd_md
    call intagm(dprarr, intarr, jdtset, marr, 9, string(1:lenstr), 'rprimd', tread, 'LEN')
    if (tread == 1) then
      Invar%rprimd_md(:,:) = TRANSPOSE(reshape(dprarr(1:9),(/3, 3/)))
    else
      write(msg,*)&
       'The NetCDF file .nc is not used.',ch10,&
       'The variable "rprimd" is thus mandatory.',ch10,&
       'ACTION : Please modify your input file'
      ABI_ERROR(msg)
    end if

! natom
    call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'natom', tread, 'INT')
    if (tread == 1) then
      Invar%natom = intarr(1)
    else
      write(msg,*)&
       'The NetCDF file .nc is not used.',ch10,&
       'The variable "natom" is thus mandatory.',ch10,&
       'ACTION : Please modify your input file'
      ABI_ERROR(msg)
    end if

! typat
    ABI_MALLOC(Invar%typat,(Invar%natom)); Invar%typat(:)=0
    call intagm(dprarr, intarr, jdtset, marr, Invar%natom, string(1:lenstr), 'typat', tread, 'INT')
    if (tread == 1) then
      Invar%typat = intarr(1:Invar%natom)
    else
      write(msg,*)&
       'The NetCDF file .nc is not used.',ch10,&
       'The variable "typat" is thus mandatory.',ch10,&
       'ACTION : Please modify your input file'
      ABI_ERROR(msg)
    end if

  end if

! multiplicity
  call intagm(dprarr, intarr, jdtset, marr, 9, string(1:lenstr), 'multiplicity', tread, 'DPR')
  if (tread == 1) then
    Invar%multiplicity(:,:) = TRANSPOSE(reshape(dprarr(1:9),(/3, 3/)))
  else
    write(msg,*)&
     'Variable "multiplicity" is mandatory, but was not found in input file.'
    ABI_ERROR(msg)
  end if

! temperature
  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'temperature', tread, 'DPR')
  if (tread == 1) then
    Invar%temperature = dprarr(1)
  else
    write(msg,*)&
     'Variable "temperature" is mandatory, but was not found in input file.'
    ABI_ERROR(msg)
  end if

! nstep_max
  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'nstep_max', tread, 'INT')
  if (tread == 1) then
    Invar%nstep_max = intarr(1)
  else
    write(msg,*)&
     'Variable "nstep_max" is mandatory, but was not found in input file.'
    ABI_ERROR(msg)
  end if

! nstep_min
  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'nstep_min', tread, 'INT')
  if (tread == 1) then
    Invar%nstep_min = intarr(1)
  else
    write(msg,*)&
     'Variable "nstep_min" is mandatory, but was not found in input file.'
    ABI_ERROR(msg)
  end if

! rcut
  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'rcut', tread, 'LEN')
  if (tread == 1) then
    Invar%rcut = dprarr(1)
  else
    write(msg,*)&
     'Variable "rcut" is mandatory, but was not found in input file.'
    ABI_ERROR(msg)
  end if

! debug_mode (optional)
  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'debug_mode', tread, 'INT')
  if (tread == 1) then
    if (intarr(1).ne.0) Invar%debug = .true.
  end if

! =========================================================================== !
! Output header and mandatory input variables

! Write version, copyright, date...
  write(Invar%stdout,*) ' '

  if (Invar%debug) then
    write(Invar%stdout,'(a,a,a)') '.Version ', version_string,' of ATDEP (Debug)'
  else
    write(Invar%stdout,'(a,a,a)') '.Version ', version_string,' of ATDEP'
  end if

  write(Invar%stdout,'(a)') '.Copyright (C) 1998-2025 ABINIT group (FB,JB,GA).'
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

  write(Invar%stdout,*) ' '
  write(Invar%stdout,*) '#############################################################################'
  write(Invar%stdout,*) '######################### ECHO OF INPUT FILE ################################'
  write(Invar%stdout,*) '#############################################################################'

  write(Invar%stdout,'(a)') ' ======================= Define the unitcell ================================='
  write(Invar%stdout,'(1x,a20,1x,i4,1x,i4)') ljust('brav',20),Invar%bravais(1),Invar%bravais(2)
  if ((Invar%bravais(1).eq.2).or.(Invar%bravais(1).eq.5)) then
    write(Invar%stdout,'(1x,a20,1x,f15.10)') 'angle',Invar%angle_alpha
  end if
  write(Invar%stdout,'(1x,a20,1x,i4)') ljust('natom_unitcell',20),Invar%natom_unitcell
  write(Invar%stdout,'(1x,a20)') ljust('xred_unitcell',20)
  do ii=1,Invar%natom_unitcell
    write(Invar%stdout,'(22x,3(f15.10,1x))') (Invar%xred_unitcell(jj,ii), jj=1,3)
  end do
  write(Invar%stdout,'(1x,a20,20(1x,i4))') ljust('typat_unitcell',20),(Invar%typat_unitcell(jj),jj=1,Invar%natom_unitcell)
  write(Invar%stdout,'(1x,a20,1x,i4)') ljust('ntypat',20),Invar%ntypat
  write(Invar%stdout,'(1x,a20,20(1x,f15.10))') ljust('amu',20),(Invar%amu(jj),jj=1,Invar%ntypat)

! znucl (optional)
  if (.not.Invar%netcdf) then
    call intagm(dprarr, intarr, jdtset, marr, Invar%ntypat, string(1:lenstr), 'znucl', tread, 'DPR')
    if (tread == 1) then
      ABI_MALLOC(Invar%znucl,(Invar%ntypat)); Invar%znucl(:)=zero
      Invar%znucl = dprarr(1:Invar%ntypat)
      write(Invar%stdout,'(1x,a20,20(1x,f15.10))') ljust('znucl',20),(Invar%znucl(jj),jj=1,Invar%ntypat)
    end if
  end if


  write(Invar%stdout,'(a)') ' ======================= Define the supercell ================================'
  write(Invar%stdout,'(1x,a20)') ljust('rprimd',20)
  do ii=1,3
    write(Invar%stdout,'(22x,3(f15.10,1x))') (Invar%rprimd_md(ii,jj),jj=1,3)
  end do
  write(Invar%stdout,'(1x,a20)') ljust('multiplicity',20)
  do ii=1,3
    write(Invar%stdout,'(22x,3(f15.10,1x))') (Invar%multiplicity(ii,jj),jj=1,3)
  end do
  write(Invar%stdout,'(1x,a20,1x,i4)') ljust('natom',20),Invar%natom
  write(Invar%stdout,'(1x,a20)') ljust('typat',20)
  do ii=1,Invar%natom,10
    if (ii+9.lt.Invar%natom) then
      write(Invar%stdout,'(22x,10(i4,1x))') (Invar%typat(ii+jj-1),jj=1,10)
    else
      write(Invar%stdout,'(22x,10(i4,1x))') (Invar%typat(jj),jj=ii,Invar%natom)
    end if
  end do

  write(Invar%stdout,'(a)') ' ======================= Define computational details ========================'
  write(Invar%stdout,'(1x,a20,1x,i5)') ljust('nstep_max',20),Invar%nstep_max
  write(Invar%stdout,'(1x,a20,1x,i5)') ljust('nstep_min',20),Invar%nstep_min
  write(Invar%stdout,'(1x,a20,1x,f15.10)') ljust('rcut',20),Invar%rcut
  write(Invar%stdout,'(1x,a20,1x,f15.10)') ljust('temperature',20),Invar%temperature


! =========================================================================== !
! Optional input variables

  write(Invar%stdout,'(a)') ' ======================= Optional input variables ============================'

! dosdeltae
  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'dosdeltae', tread, 'ENE')
  if (tread == 1) then
    Invar%dosdeltae = dprarr(1)
    write(Invar%stdout,'(1x,a20,1x,f15.10)') ljust('dosdeltae',20),Invar%dosdeltae
  end if

! use_ideal_positions
  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'use_ideal_positions', tread, 'INT')
  if (tread == 1) then
    Invar%use_ideal_positions = intarr(1)
    write(Invar%stdout,'(1x,a20,1x,i4)') ljust('use_ideal_positions',20),Invar%use_ideal_positions
  end if

! born_charge
  call intagm(dprarr, intarr, jdtset, marr, Invar%ntypat, string(1:lenstr), 'born_charge', tread, 'DPR')
  if (tread == 1) then
    ABI_MALLOC(Invar%born_charge,(Invar%ntypat)); Invar%born_charge(:)=zero
    Invar%loto=.true.
    Invar%born_charge = dprarr(1:Invar%ntypat)
    write(Invar%stdout,'(1x,a20,20(1x,f15.10))') ljust('born_charge',20),(Invar%born_charge(jj),jj=1,Invar%ntypat)
  end if

! dielec_constant
  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'dielec_constant', tread, 'DPR')
  if (tread == 1) then
    Invar%dielec_constant = dprarr(1)
    write(Invar%stdout,'(1x,a20,1x,f15.10)') ljust('dielec_constant',20),Invar%dielec_constant
  end if

! bzpath
  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'bzpath', tread, 'INT')
  if (tread == 1) then
    Invar%bzpath = intarr(1)
    write(Invar%stdout,'(1x,a20,1x,i4)') ljust('bzpath',20),Invar%bzpath

    if (Invar%bzpath.lt.0) then
      ABI_MALLOC(Invar%qpt,(3,abs(Invar%bzpath))); Invar%qpt(:,:)=zero

      call intagm(dprarr, intarr, jdtset, marr, 1+3*abs(Invar%bzpath), string(1:lenstr), 'bzpath', tread, 'DPR')
      Invar%qpt(:,:) = reshape(dprarr(2:1+3*abs(Invar%bzpath)), (/3,abs(Invar%bzpath)/))

      write(Invar%stdout,'(a)') ' Q points as given in the input file:'
      do jj=1,abs(Invar%bzpath)
        write(Invar%stdout,'(22x,3(f15.10,1x))') Invar%qpt(:,jj)
      end do

    else if (Invar%bzpath.gt.0) then
      ABI_MALLOC(Invar%special_qpt,(Invar%bzpath))
      call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'special_qpt', tread, 'KEY', key_value=msg)
      if (tread == 1) then
        jj = 1
        do ii=1,Invar%bzpath
          Invar%special_qpt(ii) = "  "
          ierr = next_token(msg,jj,Invar%special_qpt(ii))
          !jj = 2*ii-1
          !Invar%special_qpt(ii) = msg(jj:jj)
        end do
        write(Invar%stdout,'(a,1x,a2,10("-",a2))') ' Special q-points: ',Invar%special_qpt(:)
      else
        write(msg,*)&
         'Variable bzpath > 0, but variable "special_qpt" was not found in input file.'
        ABI_ERROR(msg)
      end if
    end if
  end if

! order
  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'order', tread, 'INT')
  if (tread == 1) then
    Invar%order = intarr(1)
    if (Invar%order.eq.3) then
      call intagm(dprarr, intarr, jdtset, marr, 2, string(1:lenstr), 'order', tread, 'DPR')
      Invar%rcut3 = dprarr(2)
      write(Invar%stdout,'(1x,a20,1x,i4,1x,f15.10)') ljust('order',20),Invar%order,Invar%rcut3
      if (Invar%rcut3.gt.Invar%rcut) then
        ABI_ERROR('The cutoff radius of the third order cannot be greater than the second order one.')
      end if
    else if (Invar%order.eq.4) then
      call intagm(dprarr, intarr, jdtset, marr, 3, string(1:lenstr), 'order', tread, 'DPR')
      Invar%rcut3 = dprarr(2)
      Invar%rcut4 = dprarr(3)
      write(Invar%stdout,'(1x,a20,1x,i4,2(1x,f15.10))') ljust('order',20),Invar%order,Invar%rcut3,Invar%rcut4
      if (Invar%rcut4.gt.Invar%rcut) then
        ABI_ERROR('The cutoff radius of the fourth order cannot be greater than the second order one.')
      end if
    else
      ABI_ERROR('Only the 3rd and 4th orders are allowed. Change your input file.')
    end if
  end if

! slice
  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'slice', tread, 'INT')
  if (tread == 1) then
    Invar%slice = intarr(1)
    write(Invar%stdout,'(1x,a20,1x,i4)') ljust('slice',20),Invar%slice
  end if

! enunit
  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'enunit', tread, 'INT')
  if (tread == 1) then
    Invar%enunit = intarr(1)
    if (Invar%enunit.eq.0) write(Invar%stdout,'(1x,a20,1x,i4,1x,a)') ljust('enunit',20),Invar%enunit,'(Phonon frequencies in meV)'
    if (Invar%enunit.eq.1) write(Invar%stdout,'(1x,a20,1x,i4,1x,a)') ljust('enunit',20),Invar%enunit,'(Phonon frequencies in cm-1)'
    if (Invar%enunit.eq.2) write(Invar%stdout,'(1x,a20,1x,i4,1x,a)') ljust('enunit',20),Invar%enunit,'(Phonon frequencies in mHa)'
    if (Invar%enunit.eq.3) write(Invar%stdout,'(1x,a20,1x,i4,1x,a)') ljust('enunit',20),Invar%enunit,'(Phonon frequencies in THz)'
  end if

! nproc
  call intagm(dprarr, intarr, jdtset, marr, 2, string(1:lenstr), 'nproc', tread, 'INT')
  if (tread == 1) then
    Invar%nproc(1:2) = intarr(1:2)
    write(Invar%stdout,'(1x,a20,1x,i4,1x,i4)') ljust('nproc',20),Invar%nproc(1),Invar%nproc(2)
  end if

! readifc
  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'readifc', tread, 'INT')
  if (tread == 1) then
    Invar%readifc = intarr(1)
    if (Invar%readifc.eq.1) then
      call intagm(dprarr, intarr, jdtset, marr, 2, string(1:lenstr), 'readifc', tread, 'DPR')
      Invar%tolread = dprarr(1)
      write(Invar%stdout,'(1x,a20,1x,i4,1x,f15.10)') ljust('readifc',20),Invar%readifc,Invar%tolread
    else
      write(Invar%stdout,'(1x,a20,1x,i4)') ljust('readifc',20),Invar%readifc
    end if
  end if

! alloy
  call intagm(dprarr, intarr, jdtset, marr, 3, string(1:lenstr), 'alloy', tread, 'INT')
  if (tread == 1) then
    Invar%alloy = intarr(1)
    Invar%ityp_alloy1 = intarr(2)
    Invar%ityp_alloy2 = intarr(3)
    write(Invar%stdout,'(1x,a20,1x,3(i4,1x))') ljust('alloy',20),Invar%alloy,Invar%ityp_alloy1,Invar%ityp_alloy2
  end if

! together
  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'together', tread, 'INT')
  if (tread == 1) then
    Invar%together = intarr(1)
    write(Invar%stdout,'(1x,a20,1x,i4)') ljust('together',20),Invar%together
  end if

! bzlength
  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'bzlength', tread, 'INT')
  if (tread == 1) then
    Invar%bzlength = intarr(1)
    write(Invar%stdout,'(1x,a20,1x,i4)') ljust('bzlength',20),Invar%bzlength
    ABI_MALLOC(Invar%lgth_segments,(Invar%bzlength))
    call intagm(dprarr, intarr, jdtset, marr, 1+Invar%bzlength, string(1:lenstr), 'bzlength', tread, 'DPR')
    Invar%lgth_segments = dprarr(2:1+Invar%bzlength)
    write(Invar%stdout,'(a,1x,i3,10("-",i3))') ' Length of BZ : ',Invar%lgth_segments(:)
  end if

! ngqpt1
  call intagm(dprarr, intarr, jdtset, marr, 3, string(1:lenstr), 'ngqpt1', tread, 'INT')
  if (tread == 1) then
    Invar%ngqpt1 = intarr(1:3)
    write(Invar%stdout,'(1x,a20,1x,3(i4,1x))') ljust('ngqpt1',20),Invar%ngqpt1(:)
  end if

! ngqpt2
  call intagm(dprarr, intarr, jdtset, marr, 3, string(1:lenstr), 'ngqpt2', tread, 'INT')
  if (tread == 1) then
    Invar%ngqpt2 = intarr(1:3)
    write(Invar%stdout,'(1x,a20,1x,3(i4,1x))') ljust('ngqpt2',20),Invar%ngqpt2(:)
  end if

! tolmotifinboxmatch
  !call intagm(dprarr, intarr, jdtset, marr, 3, string(1:lenstr), 'tolmotifinboxmatch', tread, 'DPR')
  !if (tread == 1) then
  !  Invar%tolmotif = dprarr(1)
  !  Invar%tolinbox = dprarr(2)
  !  Invar%tolmatch = dprarr(3)
  !  write(Invar%stdout,'(1x,a20,f10.5)') ljust('tolmotif',20),Invar%tolmotif
  !  write(Invar%stdout,'(1x,a20,f10.5)') ljust('tolinbox',20),Invar%tolinbox
  !  write(Invar%stdout,'(1x,a20,f10.5)') ljust('tolmatch',20),Invar%tolmatch
  !end if

  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'tolmotif', tread, 'DPR')
  if (tread == 1) then
    Invar%tolmotif = dprarr(1)
    write(Invar%stdout,'(1x,a20,f10.5)') ljust('tolmotif',20),Invar%tolmotif
  end if

  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'tolinbox', tread, 'DPR')
  if (tread == 1) then
    Invar%tolinbox = dprarr(1)
    write(Invar%stdout,'(1x,a20,f10.5)') ljust('tolinbox',20),Invar%tolinbox
  end if

  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'tolmatch', tread, 'DPR')
  if (tread == 1) then
    Invar%tolmatch = dprarr(1)
    write(Invar%stdout,'(1x,a20,f10.5)') ljust('tolmatch',20),Invar%tolmatch
  end if

! use_weights
  call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'use_weights', tread, 'INT')
  if (tread == 1) then
    Invar%use_weights = intarr(1)
    write(Invar%stdout,'(1x,a20,1x,i4)') ljust('use_weights',20),Invar%use_weights
  end if

  ABI_FREE(intarr)
  ABI_FREE(dprarr)

! =========================================================================== !
! Output other information and perform some checks

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
    ABI_ERROR('STOP: The value of input variable TOGETHER is not allowed')
  end if
  if ((Invar%alloy.ne.1).and.(Invar%alloy.ne.0)) then
    ABI_ERROR('STOP: The value of input variable ALLOY is not allowed')
  end if
  if (Invar%alloy.ge.1) then
    if ((Invar%ityp_alloy1.lt.1).or.(Invar%ityp_alloy2.lt.1).or.&
&       (Invar%ityp_alloy1.gt.Invar%natom_unitcell).or.(Invar%ityp_alloy2.gt.Invar%natom_unitcell)) then
      ABI_ERROR('STOP: The value of input variables IALLOY are not allowed')
    end if
  end if

! Incompatible variables :
  if ((Invar%readifc.eq.1).and.(Invar%together.eq.1).and.(Invar%order.gt.2)) then
    ABI_ERROR('STOP: readifc=1, together=1 and order=3 or 4 are incompatible')
  end if

! =========================================================================== !
! Treat virtual crystal approximation (VCA) when alloy=1.
! Redefine all the data depending on (n)typat(_unitcell) and natom_unitcell

  if (Invar%alloy.eq.1) then
    sum_alloy1=0
    sum_alloy2=0
    do iatom=1,Invar%natom
      if (Invar%typat(iatom).eq.Invar%ityp_alloy1) then
        sum_alloy1=sum_alloy1+1
      end if
      if (Invar%typat(iatom).eq.Invar%ityp_alloy2) then
        sum_alloy2=sum_alloy2+1
      end if
    end do
    amu_average   =(Invar%amu        (Invar%ityp_alloy1)*sum_alloy1+&
&                   Invar%amu        (Invar%ityp_alloy2)*sum_alloy2)/(sum_alloy1+sum_alloy2)
    if (Invar%loto) then
      born_average=(Invar%born_charge(Invar%ityp_alloy1)*sum_alloy1+&
&                   Invar%born_charge(Invar%ityp_alloy2)*sum_alloy2)/(sum_alloy1+sum_alloy2)
    end if
    shift=0
    do iatom=1,Invar%natom_unitcell
      if (Invar%typat_unitcell(iatom).lt.max(Invar%ityp_alloy1,Invar%ityp_alloy2)) then
        Invar%typat_unitcell (iatom-shift)=Invar%typat_unitcell (iatom)
        Invar%xred_unitcell(:,iatom-shift)=Invar%xred_unitcell(:,iatom)
      else if (Invar%typat_unitcell(iatom).eq.max(Invar%ityp_alloy1,Invar%ityp_alloy2)) then
        shift=shift+1
      else if (Invar%typat_unitcell(iatom).gt.max(Invar%ityp_alloy1,Invar%ityp_alloy2)) then
        Invar%typat_unitcell (iatom-shift)=Invar%typat_unitcell (iatom) - 1
        Invar%xred_unitcell(:,iatom-shift)=Invar%xred_unitcell(:,iatom)
      end if
    end do
    Invar%natom_unitcell=Invar%natom_unitcell-shift
    do iatom=1,Invar%natom
      if (Invar%typat(iatom).ge.max(Invar%ityp_alloy1,Invar%ityp_alloy2)) then
        Invar%typat(iatom)=Invar%typat(iatom) - 1
      end if
    end do
    do itypat=1,Invar%ntypat
      if (itypat.eq.min(Invar%ityp_alloy1,Invar%ityp_alloy2)) then
        Invar%amu          (itypat)=amu_average
        if (Invar%loto) then
          Invar%born_charge(itypat)=born_average
        end if
      else if (itypat.gt.max(Invar%ityp_alloy1,Invar%ityp_alloy2)) then
        Invar%amu          (itypat-1)=Invar%amu        (itypat)
        if (Invar%loto) then
          Invar%born_charge(itypat-1)=Invar%born_charge(itypat)
        end if
      end if
    end do
    Invar%ntypat=Invar%ntypat-1
!    write(6,*) 'ntypat=',Invar%ntypat
    ABI_MALLOC(typat_unitcell_tmp,(  Invar%natom_unitcell))
    ABI_MALLOC(xred_unitcell_tmp ,(3,Invar%natom_unitcell))
    ABI_MALLOC(amu_tmp           ,(  Invar%ntypat))
    if (Invar%loto) then
      ABI_MALLOC(born_charge_tmp ,(  Invar%ntypat))
    end if
    typat_unitcell_tmp (:)=Invar%typat_unitcell(1:Invar%natom_unitcell)
    xred_unitcell_tmp(:,:)=Invar%xred_unitcell(:,1:Invar%natom_unitcell)
    amu_tmp            (:)=Invar%amu(1:Invar%ntypat)
    if (Invar%loto) then
      born_charge_tmp  (:)=Invar%born_charge(1:Invar%ntypat)
    end if
    ABI_REMALLOC(Invar%typat_unitcell,(  Invar%natom_unitcell))
    ABI_REMALLOC(Invar%xred_unitcell ,(3,Invar%natom_unitcell))
    ABI_REMALLOC(Invar%amu           ,(  Invar%ntypat))
    if (Invar%loto) then
      ABI_REMALLOC(Invar%born_charge ,(  Invar%ntypat))
    end if
    Invar%typat_unitcell (:)=typat_unitcell_tmp (:)
    Invar%xred_unitcell(:,:)=xred_unitcell_tmp(:,:)
    Invar%amu            (:)=amu_tmp            (:)
    if (Invar%loto) then
      Invar%born_charge  (:)=born_charge_tmp (:)
    end if
    if (allocated(Invar%znucl)) then
      ABI_MALLOC(znucl_tmp,(Invar%ntypat))
      znucl_tmp(:) = Invar%znucl(1:Invar%ntypat)
      ABI_REMALLOC(Invar%znucl,(Invar%ntypat))
      Invar%znucl(:) = znucl_tmp(:)
      ABI_FREE(znucl_tmp)
    end if
    ABI_FREE(typat_unitcell_tmp)
    ABI_FREE(xred_unitcell_tmp)
    ABI_FREE(amu_tmp)
    if (Invar%loto) then
      ABI_FREE(born_charge_tmp)
    end if

    write(Invar%stdout,'(a)') ' ==================== Virtual Crystal Approximation =========================='
    write(Invar%stdout,'(a)') ' ================ Several input variables are modified ======================='
    write(Invar%stdout,'(a)') ' --> Beginning of the modifications'
    write(Invar%stdout,'(1x,a20,1x,i4)') ljust('ntypat',20),Invar%ntypat
    write(Invar%stdout,'(1x,a20,1x,i4)') ljust('natom_unitcell',20),Invar%natom_unitcell
    write(Invar%stdout,'(1x,a20,20(1x,f15.10))') ljust('amu',20),(Invar%amu(jj),jj=1,Invar%ntypat)
    write(Invar%stdout,'(1x,a20,20(1x,i4))') ljust('typat_unitcell',20),(Invar%typat_unitcell(jj),jj=1,Invar%natom_unitcell)
    write(Invar%stdout,'(1x,a20)') ljust('xred_unitcell',20)
    do ii=1,Invar%natom_unitcell
      write(Invar%stdout,'(22x,3(f15.10,1x))') (Invar%xred_unitcell(jj,ii), jj=1,3)
    end do
    if (Invar%loto) then
      write(Invar%stdout,'(1x,a20,20(1x,f15.10))') ljust('born_charge',20),(Invar%born_charge(jj),jj=1,Invar%ntypat)
    end if
    write(Invar%stdout,'(1x,a20)') ljust('typat',20)
    do ii=1,Invar%natom,10
      if (ii+9.lt.Invar%natom) then
        write(Invar%stdout,'(22x,10(i4,1x))') (Invar%typat(ii+jj-1),jj=1,10)
      else
        write(Invar%stdout,'(22x,10(i4,1x))') (Invar%typat(jj),jj=ii,Invar%natom)
      end if
    end do
    write(Invar%stdout,'(a)') ' --> End of the modifications'
    write(Invar%stdout,'(a)') ' '
  end if

! Compute Nstep as a function of the slice
  Invar%nstep_tot=int(float(Invar%nstep_max-Invar%nstep_min)/float(Invar%slice)+1)
  write(Invar%stdlog,*) 'nstep_tot=',Invar%nstep_tot


 end subroutine tdep_read_input

! ---------------------------------------------------------------------------

 subroutine tdep_distrib_data(Hist,Invar,MPIdata)

  type(Input_type), intent(inout) :: Invar
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
  ABI_MALLOC(Invar%weights,(Invar%my_nstep))             ; Invar%weights(:)=0.d0
  this_istep=0
  jstep=0
  if (Invar%use_weights.eq.1) then
    open(unit=30,file=trim(Invar%input_prefix)//'_weights.dat')
  else if (Invar%use_weights.eq.0) then
    Invar%weights=1.0d0/real(Invar%nstep_tot)
  endif
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
    write(Invar%stdout,'(a)') ' The positions, forces and energies are extracted from the NetCDF file: HIST.nc'
  else
    open(unit=60,file=trim(Invar%input_prefix)//'_fcart.dat')
    open(unit=50,file=trim(Invar%input_prefix)//'_xred.dat')
    open(unit=40,file=trim(Invar%input_prefix)//'_etot.dat')
    do istep=1,Invar%nstep_min-1
      if (Invar%use_weights.eq.1) then
         read(30,*) tmp1
      endif
      read(40,*) tmp1
      do iatom=1,Invar%natom
        read(50,*) tmp1,tmp2,tmp3
        read(60,*) tmp1,tmp2,tmp3
      end do
    end do
    do istep=Invar%nstep_min,Invar%nstep_max
      if (mod(istep-Invar%nstep_min,Invar%slice).ne.0) then
        if (Invar%use_weights.eq.1) then
           read(30,*) tmp1
        endif
        read(40,*) tmp1
        do iatom=1,Invar%natom
          read(50,*) tmp1,tmp2,tmp3
          read(60,*) tmp1,tmp2,tmp3
        end do
      else
        jstep=jstep+1
        if (.not.MPIdata%my_step(jstep)) then
          if (Invar%use_weights.eq.1) then
             read(30,*) tmp1
          endif
          read(40,*) tmp1
          do iatom=1,Invar%natom
            read(50,*) tmp1,tmp2,tmp3
            read(60,*) tmp1,tmp2,tmp3
          end do
        else
          this_istep=this_istep+1
          if (Invar%use_weights.eq.1) then
             read(30,*) Invar%weights(this_istep)
          endif
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
    write(Invar%stdout,'(2a)') ' The positions, forces and energies are extracted from the ASCII files:',&
&                              ' xred.dat, fcart.dat & etot.dat'
  end if !netcdf
  if (Invar%use_weights.eq.1) then
    close(30)
  end if

 end subroutine tdep_distrib_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!FB subroutine tdep_init_MPIshell(Invar,MPIdata)
!FB
!FB
!FB  type(Input_type), intent(in) :: Invar
!FB  type(MPI_enreg_type), intent(in) :: MPIdata
!FB
!FB end subroutine tdep_init_MPIshell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_init_MPIdata(Invar,MPIdata)

  type(Input_type), intent(in) :: Invar
  type(MPI_enreg_type), intent(out) :: MPIdata
  integer :: ii,remain,ierr,iproc,istep
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
    ABI_WARNING('The parallelization is performed over steps')
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
  ABI_MALLOC(sizecart,(dimcart))
  ABI_MALLOC(periode,(dimcart))
  sizecart(1)=MPIdata%nproc_shell ! MPIdata%nproc_shell
  sizecart(2)=MPIdata%nproc_step  ! MPIdata%nproc_step
  periode(:)=.false.;reorder=.false.
  call MPI_CART_CREATE(xmpi_world,dimcart,sizecart,periode,reorder,commcart_2d,ierr)
  ABI_FREE(periode)
  ABI_FREE(sizecart)

! Find the index and coordinates of the current processor
  call MPI_COMM_RANK(commcart_2d,me_cart_2d,ierr)
  ABI_MALLOC(coords,(dimcart))
  call MPI_CART_COORDS(commcart_2d,me_cart_2d,dimcart,coords,ierr)
  MPIdata%me_shell=coords(1)
  MPIdata%me_step =coords(2)
  ABI_FREE(coords)
  if ((MPIdata%me_shell == MPIdata%master).and.(MPIdata%me_step == MPIdata%master)) then
    MPIdata%iam_master = .true.
  end if

  ABI_MALLOC(keepdim,(dimcart))
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
  ABI_FREE(keepdim)
  call xmpi_comm_free(commcart_2d)

! Write some data
  write(message,'(a21,2(1x,i4))') '-Number of processors',MPIdata%nproc_shell,MPIdata%nproc_step
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
  write(Invar%stdout,'(a,1x,i4)') ' All quantities are computed from nstep_min=',Invar%nstep_min
  write(Invar%stdout,'(a,1x,i4)') '                               to nstep_max=',Invar%nstep_max
  if (Invar%slice.ne.1) then
    write(Invar%stdout,'(a,1x,i4)') '                                    by using a slice=',Invar%slice
  end if
  write(Invar%stdout,'(a,1x,i4)') ' So, the real number of time steps is nstep=',Invar%nstep_tot
  if (MPIdata%nproc_step.gt.1) then
    write(Invar%stdout,'(a,1000(1x,i5))') '-Distribution of number of steps wrt the number of processors=',MPIdata%nstep_all(:)
  end if

  ABI_MALLOC(MPIdata%nstep_acc,(MPIdata%nproc_step+1)); MPIdata%nstep_acc(:)=zero
  MPIdata%nstep_acc(1)=0
  do ii=2,MPIdata%nproc_step+1
    MPIdata%nstep_acc(ii)=MPIdata%nstep_acc(ii-1)+MPIdata%nstep_all(ii-1)
  end do
  if (MPIdata%nstep_acc(MPIdata%nproc_step+1).ne.Invar%nstep_tot) then
    write(Invar%stdlog,*) 'STOP : pb in nstep_acc'
    stop
  end if

  ABI_MALLOC(tab_step,(Invar%nstep_tot)); tab_step(:)=zero
  ABI_MALLOC(MPIdata%my_step ,(Invar%nstep_tot)); MPIdata%my_step (:)=.false.
  do iproc=1,MPIdata%nproc_step
    do istep=1,Invar%nstep_tot
      if ((istep.gt.MPIdata%nstep_acc(iproc)).and.(istep.le.MPIdata%nstep_acc(iproc+1))) then
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
  ABI_FREE(tab_step)

 end subroutine tdep_init_MPIdata

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_destroy_mpidata(MPIdata)

  type(MPI_enreg_type), intent(inout) :: MPIdata

  ABI_FREE(MPIdata%shft_step)
  ABI_FREE(MPIdata%nstep_acc)
  ABI_FREE(MPIdata%nstep_all)
  ABI_FREE(MPIdata%my_step)

 end subroutine tdep_destroy_mpidata

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_destroy_invar(Invar)

  type(Input_type), intent(inout) :: Invar

  ABI_FREE(Invar%amu)
  ABI_FREE(Invar%typat)
  ABI_FREE(Invar%xred_unitcell)
  ABI_FREE(Invar%typat_unitcell)
  if (Invar%bzpath.lt.0) then
    ABI_FREE(Invar%qpt)
  else if (Invar%bzpath.gt.0) then
    ABI_FREE(Invar%special_qpt)
    end if
  if (Invar%bzlength.gt.0) then
    ABI_FREE(Invar%lgth_segments)
  end if
  ABI_FREE(Invar%xred)
  ABI_FREE(Invar%fcart)
  ABI_FREE(Invar%etot)
  ABI_FREE(Invar%weights)
  ABI_FREE(Invar%xred_ideal)
  if (Invar%loto) then
    ABI_FREE(Invar%born_charge)
  end if
  if (allocated(Invar%znucl)) then
    ABI_FREE(Invar%znucl)
  end if

 end subroutine tdep_destroy_invar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module m_tdep_readwrite
