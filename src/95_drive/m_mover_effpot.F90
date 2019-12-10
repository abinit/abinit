!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_mover_effpot
!! NAME
!!  m_mover_effpot
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (AM)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_mover_effpot

 use defs_basis
 use m_errors
 use m_abicore
 use m_dtset
 use m_dtfil
 use m_abimover
 use m_build_info
 use m_scf_history
 use defs_wvltypes
 use m_xmpi
 use m_phonons
 use m_strain
 use m_effective_potential_file
 use m_supercell
 use m_psps
 use m_args_gs
 use m_ifc

 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes, only : MPI_type
 use m_geometry, only : xcart2xred, xred2xcart
 use m_multibinit_dataset, only : multibinit_dtset_type
 use m_effective_potential,only : effective_potential_type
 use m_fit_polynomial_coeff, only : polynomial_coeff_writeXML, &
   fit_polynomial_coeff_fit, genereList, fit_polynomial_coeff_getPositive,fit_polynomial_coeff_getCoeffBound
 use m_polynomial_coeff,only : polynomial_coeff_getNorder
! use m_pawang,       only : pawang_type, pawang_free
! use m_pawrad,       only : pawrad_type, pawrad_free
! use m_pawtab,       only : pawtab_type, pawtab_nullify, pawtab_free
! use m_pawxmlps, only : paw_setup, ipsp2xml, rdpawpsxml, &
!&                       paw_setup_copy, paw_setup_free, getecutfromxml
 use m_abihist, only : abihist
 use m_ewald
 use m_mpinfo,           only : init_mpi_enreg,destroy_mpi_enreg
 use m_copy            , only : alloc_copy
 use m_electronpositron, only : electronpositron_type
 use m_scfcv,            only : scfcv_t, scfcv_run,scfcv_destroy
 use m_results_gs,       only : results_gs_type,init_results_gs,destroy_results_gs
 use m_mover,            only : mover
 use m_io_tools,         only : get_unit, open_file

 implicit none

 private
!!***

 public :: mover_effpot
!!***

contains
!!***

!!****f* ABINIT/mover_effpot
!! NAME
!! mover_effpot
!!
!! FUNCTION
!! this routine is driver for using mover with effective potential
!!
!! INPUTS
!!  inp = input of multibinit
!!  effective_potential =  effective potential of the reference structure
!!  option = flag for the option:
!!                         -1  = Bound the anharmonic part
!!                          12 = NVT simulation
!!                          13 = NPT simulation
!! OUTPUT
!!
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!!
!! PARENTS
!!      multibinit
!!
!! CHILDREN
!!      alloc_copy,destroy_mpi_enreg,destroy_results_gs,dtset_free
!!      effective_potential_setcoeffs,effective_potential_setsupercell
!!      fit_polynomial_coeff_fit,fit_polynomial_coeff_getpositive,generelist
!!      init_results_gs,mover,polynomial_coeff_free,polynomial_coeff_getnorder
!!      polynomial_coeff_init,polynomial_coeff_setcoefficient
!!      polynomial_coeff_writexml,scfcv_destroy,wrtout,xcart2xred,xmpi_barrier
!!      xred2xcart
!!
!! SOURCE

subroutine mover_effpot(inp,filnam,effective_potential,option,comm,hist)


!Arguments --------------------------------
!scalar
 integer, intent(in) :: option,comm
!array
 type(multibinit_dtset_type),intent(in) :: inp
 type(effective_potential_type),intent(inout)  :: effective_potential
 character(len=fnlen),intent(in) :: filnam(15)
 type(abihist),optional,intent(inout):: hist
!Local variables-------------------------------
!scalar
 integer :: filetype,icoeff_bound,ii, unit_out
!integer :: iexit,initialized
 integer :: jj,kk,nproc,icoeff,ncoeff,nmodels,ncoeff_bound,ncoeff_bound_tot,ncoeff_max
 integer :: model_bound,model_ncoeffbound,my_rank
!integer :: mtypalch,,npsp,paw_size,type
!integer,save :: paw_size_old=-1
 real(dp):: cutoff,freq_q,freq_b,qmass,bmass
! real(dp):: time_q,time_b
 logical :: iam_master
 integer, parameter:: master=0
 logical :: verbose,writeHIST
!real(dp) :: cpui
!character(len=6) :: codvsn

!TEST_AM
! integer :: ia,mu,rand_seed = 5
! real(dp):: mass_ia,rescale_vel,sum_mass,v2gauss
!TEST_AM
! Set array dimensions
 character(len=500) :: message
 type(MPI_type),target :: mpi_enreg
 type(dataset_type),target :: dtset
 type(scfcv_t) :: scfcv_args
 type(datafiles_type),target :: dtfil
 integer,target :: zero_integer
 type(ab_xfh_type) :: ab_xfh
 type(results_gs_type),target :: results_gs
 type(pseudopotential_type),target :: psps
!arrays
!no_abirules
 integer :: sc_size(3),sc_size_TS(3)
 integer,pointer :: indsym(:,:,:)
 integer,allocatable :: listcoeff(:),listcoeff_bound(:,:),list_tmp(:),list_bound(:,:)
 integer,allocatable :: isPositive(:)
 integer,allocatable :: symrel(:,:,:)
!integer,allocatable :: npwtot(:)
 real(dp) :: acell(3)
!real(dp) :: ecut_tmp(3,2,10)
 real(dp),allocatable :: coeff_values(:,:)
 real(dp),pointer :: rhog(:,:),rhor(:,:)
 real(dp),allocatable :: tnons(:,:)
 real(dp),allocatable :: xred(:,:),xred_old(:,:),xcart(:,:)
 real(dp),allocatable :: fred(:,:),fcart(:,:)
 real(dp),allocatable :: vel(:,:)
 real(dp) :: vel_cell(3,3),rprimd(3,3)
 type(polynomial_coeff_type),dimension(:),allocatable :: coeffs_all,coeffs_tmp,coeffs_bound
 character(len=fnlen) :: filename
 character(len=50) :: name_file
 character(len=200):: term_name
!character(len=fnlen) :: filename_psp(3)
 type(electronpositron_type),pointer :: electronpositron
! type(pspheader_type),allocatable :: pspheads(:)
! type(pawrad_type),allocatable :: pawrad(:)
! type(pawtab_type),allocatable :: pawtab(:)
! type(args_gs_type) :: args_gs
!type(wvl_data) :: wvl
!type(pawang_type) :: pawang
!type(scf_history_type) :: scf_history

!******************************************************************

!MPI variables
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

 write(message, '(a,(80a),a)') ch10,&
& ('=',ii=1,80),ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

!*******************************************************************
! 1 Generate and check supercell for the dynamics
!*******************************************************************

!a new supercell is compute
!Initialisaton of variable
 if(option == -1.or.option == -2) then
!  Bound process option
   sc_size(:) = inp%bound_cell
 else if(option == -3) then
!  Heff option
   sc_size(:) = (/1,1,1/)
 else
!  Normal dynamics
   sc_size(:) = inp%ncell
 end if

 if(option/=0)then

!  acell is always set to one, only rprimd is used for the effective potential
   acell = one
   rprimd = effective_potential%crystal%rprimd

   ABI_ALLOCATE(xred,(3,effective_potential%crystal%natom))
   ABI_ALLOCATE(xcart,(3,effective_potential%crystal%natom))

!  convert new xcart
   call xcart2xred(effective_potential%crystal%natom,effective_potential%crystal%rprimd,&
&   effective_potential%crystal%xcart,xred)
   call xred2xcart(effective_potential%crystal%natom, rprimd, xcart, xred)
!  Generate supercell for the simulation
   call effective_potential_setSupercell(effective_potential,comm,ncell=sc_size)

   ABI_DEALLOCATE(xred)
   ABI_DEALLOCATE(xcart)

!***************************************************************
!1 Convert some parameters into the structures used by mover.F90
!***************************************************************
!NOTE:ARGUMENTS OF MOVER SHOULD BE CLEAN
!     We may just need to provide AB_MOVER wich is the main object
!     for mover and set  scfcv_args as an optional and depending on
!     the kind of calculation (abinit or multibinit), we provide
!     to mover scfcv_args or effective_potential...
!***************************************************************
!  Free dtset
   call dtset%free()

!  Set mpi_eng
   mpi_enreg%comm_cell  = comm
   mpi_enreg%me = my_rank

!  Set the abinit dataset for mover with fake values
!  Scalar
   dtset%dmft_entropy = 0
   dtset%nctime = inp%nctime ! NetCdf TIME between output of molecular dynamics informations
   dtset%delayperm = 0  ! DELAY between trials to PERMUTE atoms
   dtset%dilatmx = one  ! DILATation : MaXimal value
   dtset%chkdilatmx = 0 ! No check on dilatmx is needed in multibilint
   dtset%diismemory = 8 ! Direct Inversion in the Iterative Subspace MEMORY
   dtset%friction = 0.0001d0 ! internal FRICTION coefficient
   dtset%goprecon = 0   ! Geometry Optimization PREconditioner equations
   dtset%istatr = 0     ! Integer for STATus file SHiFT
   dtset%jellslab = 0   ! include a JELLium SLAB in the cell
   dtset%mqgrid = 0     ! Maximum number of Q-space GRID points for pseudopotentials
   dtset%mqgriddg = 0   ! Maximum number of Q-wavevectors for the 1-dimensional GRID
                        ! for the Double Grid in PAW
   dtset%mdwall = 10000d0 ! Molecular Dynamics WALL location
   dtset%ntypalch = 0   ! Number of TYPe of atoms that are "ALCHemical"
   dtset%natom = effective_potential%supercell%natom
   dtset%ntypat = effective_potential%crystal%ntypat
   dtset%npspalch = effective_potential%crystal%ntypat
   dtset%nconeq = 0     ! Number of CONstraint EQuations
   dtset%noseinert = 1.d-5 ! NOSE INERTia factor
   dtset%nnos = inp%nnos   ! Number of nose masses Characteristic
   do ii=1,3 ! Only copy the diagonal part
             ! might be adapt with ngqpt(9) instead of ngqpt(3) because is needed in wght9.f
     dtset%ph_ngqpt(ii) = inp%ngqpt((ii))
   end do
   dtset%ph_nqshift = inp%nqshft
   dtset%nsym = 1       ! Number of SYMmetry operations
   dtset%prtxml = 0     ! print the xml
   dtset%signperm = 1   ! SIGN of PERMutation potential
   dtset%strprecon = 1  ! STRess PRECONditioner
   dtset%supercell_latt(:,:) = 0
   do ii=1,3
     dtset%supercell_latt(ii,ii) = sc_size(ii)
   end do
   dtset%tolmxf = 2.0d-5
   dtset%tsmear = 0.009500446 !
   dtset%vis = 100      ! VIScosity
   dtset%usewvl = 0     !
   dtset%useylm = 0     !

!    if(option == -3) then
!      write(message,'(a)')' Read the DDB file to fill the dtset array'
!      call wrtout(std_out,message,"COLL")
! !    Copy real informtions from the ddb
!      call effective_potential_file_getType(filnam(3),type)
!      if(type /= 1) then
!        write(message, '(5a)' )&
! &         ' You need to provide DDB file in the input to compute ahnarmonic',ch10,&
! &         ' part of effective Hamiltionian',ch10,&
! &         'Action: add DDB file in the inputs'
!        MSG_BUG(message)
!      end if
!      call ddb_to_dtset(comm, dtset,filnam(3),psps)
!      ABI_ALLOCATE(dtset%kptns,(3,dtset%nkpt))
!      dtset%kptns(:,:) = dtset%kpt(:,:)
!      ABI_ALLOCATE(dtset%istwfk,(dtset%nkpt))
!      dtset%istwfk(:) = 1

!    else
!    Need to init some values
   ABI_ALLOCATE(symrel,(3,3,dtset%nsym))
   symrel = 1
   call alloc_copy(symrel,dtset%symrel)
   ABI_ALLOCATE(tnons,(3,dtset%nsym))
   tnons = zero
   call alloc_copy(tnons,dtset%tnons)
   call alloc_copy(effective_potential%supercell%typat,dtset%typat)
   call alloc_copy(effective_potential%crystal%znucl,dtset%znucl)
   ABI_DEALLOCATE(symrel)
   ABI_DEALLOCATE(tnons)
!   end if

   !array
   ABI_ALLOCATE(dtset%iatfix,(3,dtset%natom)) ! Indices of AToms that are FIXed
   dtset%iatfix = 0
   dtset%goprecprm(:) = zero !Geometry Optimization PREconditioner PaRaMeters equations
   ABI_ALLOCATE(dtset%prtatlist,(dtset%natom)) !PRinT by ATom LIST of ATom
   dtset%prtatlist(:) = 0
   ABI_ALLOCATE(dtset%mixalch_orig,(dtset%npspalch,dtset%ntypalch,1))
   dtset%mixalch_orig(:,:,:)=zero
   ABI_ALLOCATE(dtset%ph_qshift,(3,dtset%ph_nqshift))
   dtset%ph_qshift = inp%q1shft
   if(option  > 0)then
     verbose = .TRUE.
     writeHIST = .TRUE.
     dtset%dtion = inp%dtion  ! Delta Time for IONs
     dtset%ionmov = inp%dynamics  ! Number for the dynamic
     dtset%ntime = inp%ntime  ! Number of TIME steps
     dtset%optcell = inp%optcell    ! OPTimize the CELL shape and dimensions Characteristic
     dtset%restartxf = inp%restartxf  ! RESTART from (X,F) history
     dtset%mdtemp(1) = inp%temperature   !Molecular Dynamics Temperatures
     dtset%mdtemp(2) = inp%temperature   !Molecular Dynamics Temperatures
     dtset%strtarget(1:6) = -1 * inp%strtarget(1:6) / 29421.033d0 ! STRess TARGET
   else if(option == -1.or.option == -2) then
!    Set default for the fit
     verbose = .false.
     writeHIST = .false.
     dtset%restartxf = 0  ! RESTART from (X,F) history
     dtset%dtion = inp%dtion  ! Delta Time for IONs
     dtset%ionmov = 13  ! Number for the dynamic
     dtset%ntime = inp%bound_step  ! Number of TIME steps
     dtset%optcell = 2    ! OPTimize the CELL shape and dimensions Characteristic
     dtset%mdtemp(1) = inp%bound_temp   !Molecular Dynamics Temperatures
     dtset%mdtemp(2) = inp%bound_temp !Molecular Dynamics Temperatures
     dtset%strtarget(1:6) = zero
   end if

!  Set the barostat and thermonstat if ionmov == 13
   if(dtset%ionmov == 13)then

!    Select frequency of the barostat as a function of temperature
!    For small temperature, we need huge barostat and inversely
     ! if(dtset%mdtemp(1) <= 10) then
     !   freq_q = 0.002
     !   freq_b = 0.0002
     ! else  if(dtset%mdtemp(1) <= 50) then
     !   freq_q = 0.02
     !   freq_b = 0.002
     ! else  if(dtset%mdtemp(1) > 50.and.dtset%mdtemp(1) < 300) then
     !   freq_q = 0.1
     !   freq_b = 0.01
     ! else
     !   freq_q = 0.2
     !   freq_b = 0.02
     ! end if

!TEST_AM_old way
!     freq_q = 0.1
!     freq_b = 0.01
!     qmass = dtset%natom* kb_THzK * dtset%mdtemp(1) / (freq_q**2)
!     bmass = dtset%natom* kb_THzK * dtset%mdtemp(1) / (freq_b**2)
!TEST_AM


!TEST_AM
     freq_q = 800  / Ha_cmm1
     freq_b = 800  / Ha_cmm1
     qmass = 10 * dtset%natom * kb_HaK * dtset%mdtemp(1) / (freq_q**2)
     bmass = 10000*qmass
!TEST_AM


     if(dtset%nnos==0) then
       dtset%nnos = 1
       ABI_ALLOCATE(dtset%qmass,(dtset%nnos))
       dtset%qmass(:)  = qmass
       write(message,'(3a,F30.10,a)')&
&       ' WARNING: nnos is set to zero in the input',ch10,&
&       '          value by default for qmass: ',dtset%qmass(:),ch10
       if(verbose)call wrtout(std_out,message,"COLL")
     else
       ABI_ALLOCATE(dtset%qmass,(dtset%nnos)) ! Q thermostat mass
       dtset%qmass(:) = inp%qmass(:)
     end if
     if (abs(inp%bmass) < tol10) then
       dtset%bmass = bmass
       write(message,'(3a,F30.10,a)')&
&       ' WARNING: bmass is set to zero in the input',ch10,&
&       '          value by default for bmass: ',dtset%bmass,ch10
       if(verbose)call wrtout(std_out,message,"COLL")
     else
       dtset%bmass = inp%bmass  ! Barostat mass
     end if
   end if

!  Do a check
   if(dtset%ionmov == 27)then
     call effective_potential_file_getType(filnam(3),filetype)
     if(filetype /= 1)then
       write(message, '(5a)' )&
&       ' The file ',trim(filnam(3)),' is not a DDB',ch10,&
&       ' It is not compatible with ionmov 27'
       MSG_ERROR(message)
     end if

   end if

!  set psps
   psps%useylm = dtset%useylm
!    if(option == -3)then
!      mtypalch = 0
!      npsp = dtset%ntypat
!      call psps_free(psps)
!      filename_psp(1) = "/home/alex/calcul/psp/Sr.LDA_PW-JTH.xml"
!      filename_psp(2) = "/home/alex/calcul/psp/Ti.LDA_PW-JTH.xml"
!      filename_psp(3) = "/home/alex/calcul/psp/O.LDA_PW-JTH.xml"
!      ABI_DATATYPE_ALLOCATE(pspheads,(npsp))
!      call inpspheads(filename_psp,npsp,pspheads,ecut_tmp)
!      call psps_init_global(mtypalch, npsp, psps, pspheads)
!      call psps_init_from_dtset(dtset, 1, psps, pspheads)
!    end if

! !  The correct dimension of pawrad/tab is ntypat. In case of alchemical psps
! !  pawrad/tab(ipsp) is invoked with ipsp<=npsp. So, in order to avoid any problem,
! !  declare pawrad/tab at paw_size=max(ntypat,npsp).
!    paw_size=0;if (psps%usepaw==1) paw_size=max(dtset%ntypat,dtset%npsp)
!    if (paw_size/=paw_size_old) then
!      if (paw_size_old/=-1) then
!        call pawrad_free(pawrad)
!        call pawtab_free(pawtab)
!        ABI_DATATYPE_DEALLOCATE(pawrad)
!        ABI_DATATYPE_DEALLOCATE(pawtab)
!      end if
!      ABI_DATATYPE_ALLOCATE(pawrad,(paw_size))
!      ABI_DATATYPE_ALLOCATE(pawtab,(paw_size))
!      call pawtab_nullify(pawtab)
!      paw_size_old=paw_size
!    end if

!  set args_gs
!    if (option == -3 then
!      call args_gs_init(args_gs, &
! &       effective_potential%crystal%amu(:),dtset%mixalch_orig(:,:,1),&
! &       dtset%dmatpawu(:,:,:,:,1),dtset%upawu(:,1),dtset%jpawu(:,1),&
! &       dtset%rprimd_orig(:,:,1))
!      ABI_ALLOCATE(npwtot,(dtset%nkpt))
!    end if
!  initialisation of results_gs
   call init_results_gs(dtset%natom,1,1,results_gs)

!  Set the pointers of scfcv_args
   zero_integer = 0
   scfcv_args%dtset     => dtset
   ABI_ALLOCATE(indsym,(4,dtset%nsym,dtset%natom))
   indsym = 0
   scfcv_args%indsym => indsym
   scfcv_args%mpi_enreg => mpi_enreg
   scfcv_args%ndtpawuj  => zero_integer
   scfcv_args%results_gs => results_gs
   scfcv_args%psps => psps
!  Set other arguments of the mover.F90 routines

!  Set the dffil structure
   dtfil%filnam_ds(1:2)=filnam(1:2)
   dtfil%filnam_ds(3)=""
   dtfil%filnam_ds(4)=filnam(2)
   dtfil%filstat='_STATUS'
   nullify (electronpositron)
   ABI_ALLOCATE(rhog,(2,1))
   ABI_ALLOCATE(rhor,(2,1))

!  Initialize xf history (should be put in inwffil)
!  Not yet implemented for ionmov 2 3 10 11 22 (memory problem...)
!  ab_xfh%mxfh=(ab_xfh%nxfh-dtset%restartxf+1)+dtset%ntime+5
   ab_xfh%nxfh = 0
   ab_xfh%mxfh = 1
   ABI_ALLOCATE(ab_xfh%xfhist,(3,dtset%natom+4,2,ab_xfh%mxfh))
   if (any((/2,3,10,11,22/)==dtset%ionmov)) then
     write(message, '(3a)' )&
&     ' This dynamics can not be used with effective potential',ch10,&
&     'Action: correct dynamics input'
     MSG_BUG(message)
   end if

!***************************************************************
!2  initialization of the structure for the dynamics
!***************************************************************

   if (allocated(dtset%rprimd_orig)) then
     ABI_DEALLOCATE(dtset%rprimd_orig)
   end if
   ABI_ALLOCATE(dtset%rprimd_orig,(3,3,1))
   dtset%rprimd_orig(:,:,1) = effective_potential%supercell%rprimd


   ABI_ALLOCATE(xred,(3,dtset%natom))
   ABI_ALLOCATE(xred_old,(3,dtset%natom))
   ABI_ALLOCATE(vel,(3,dtset%natom))
   ABI_ALLOCATE(fred,(3,dtset%natom))
   ABI_ALLOCATE(fcart,(3,dtset%natom))

   call xcart2xred(dtset%natom,effective_potential%supercell%rprimd,&
&   effective_potential%supercell%xcart,xred)

   xred_old = xred
   vel_cell(:,:) = zero
   vel(:,:)      = zero

!*********************************************************
!4   Call main routine for the bound process,
!    monte carlo / molecular dynamics / project
!*********************************************************
   if(option > 0)then
     !*************************************************************
     !  call mover in case of NPT or NVT simulation
     !*************************************************************
     write(message, '((80a),3a)' ) ('-',ii=1,80), ch10,&
&     '-Monte Carlo / Molecular Dynamics ',ch10

     ! Marcus: if wanted analyze anharmonic terms of effective potential &&
     ! and print anharmonic contribution to file anharmonic_energy_terms.out
     ! Open File and write header
     ncoeff = effective_potential%anharmonics_terms%ncoeff
     name_file='anharmonic_energy_terms.out'
     unit_out = get_unit()
     if(inp%analyze_anh_pot == 1)then
       open(unit=unit_out,file=name_file,status='replace',form='formatted')
       write(unit_out,*) '#---------------------------------------------#'
       write(unit_out,*) '#    Anharmonic Terms Energy Contribution     #'
       write(unit_out,*) '#---------------------------------------------#'
       write(unit_out,*) ''
       write(unit_out,'(A,I5)') 'Number of Terms: ', ncoeff
       write(unit_out,*) ''
       write(unit_out,'(A)') 'Terms     Names'
       do icoeff=1,ncoeff
         term_name = effective_potential%anharmonics_terms%coefficients(icoeff)%name
         write(unit_out,'(I5,A,A)') icoeff,'     ',trim(term_name)
       enddo
       write(unit_out,*) ''
       write(unit_out,'(A)',advance='no')  'Cycle/Terms'
       do icoeff=1,ncoeff
         if(icoeff<ncoeff)then
         write(unit_out,'(I5)',advance='no') icoeff
         else
         write(unit_out,'(I5)',advance='yes') icoeff
         endif
       enddo
     end if

     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     call mover(scfcv_args,ab_xfh,acell,effective_potential%crystal%amu,dtfil,electronpositron,&
&     rhog,rhor,dtset%rprimd_orig,vel,vel_cell,xred,xred_old,&
&     effective_potential=effective_potential,filename_ddb=filnam(3),&
&     verbose=verbose,writeHIST=writeHIST)
     close(unit_out)
   else if(option== -1.or.option==-2)then
     !*************************************************************
     !   Try to bound the model
     !*************************************************************
     write(message, '((80a),4a)' ) ('-',ii=1,80), ch10,&
&     ' Try to bound the model',ch10,' Check if the model is bounded or not'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')

!    Try the model
     call mover(scfcv_args,ab_xfh,acell,effective_potential%crystal%amu,dtfil,electronpositron,&
&     rhog,rhor,dtset%rprimd_orig,vel,vel_cell,xred,xred_old,&
&     effective_potential=effective_potential,verbose=verbose,writeHIST=writeHIST)

     write(message, '(a)' ) ' => The model'
     if(effective_potential%anharmonics_terms%bounded)then
       write(message, '(2a)' ) trim(message),' is bound'
       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,message,'COLL')
     else
       write(message, '(2a)' ) trim(message),' is not bound'
       call wrtout(std_out,message,'COLL')


       if(option==-2)then

!        Fill the list for the fixcoeff input of the fit_polynomial_coeff_fit routine
!        Store the number of coefficients before adding other coeff for the bounding
         ncoeff = effective_potential%anharmonics_terms%ncoeff
         ABI_ALLOCATE(listcoeff,(ncoeff))
         do ii=1,ncoeff
           listcoeff(ii)=ii
         end do

         write(message, '(2a)')ch10,' Generate the list of addionnal terms with the fit process...'
         call wrtout(std_out,message,'COLL')

!        Get the additional coeff
         call fit_polynomial_coeff_fit(effective_potential,(/0/),listcoeff,hist,1,&
&         inp%bound_rangePower,0,inp%bound_maxCoeff,ncoeff,1,comm,cutoff_in=inp%bound_cutoff,&
&         max_power_strain=2,verbose=.true.,positive=.true.,spcoupling=inp%bound_SPCoupling==1,&
&         anharmstr=inp%bound_anhaStrain==1,only_even_power=.true.)

!        Store the max number of coefficients after the fit process
         ncoeff_max = effective_potential%anharmonics_terms%ncoeff
!        Store all the coefficients in coeffs_all
         ABI_DATATYPE_ALLOCATE(coeffs_all,(ncoeff_max))
         ABI_DATATYPE_ALLOCATE(coeffs_tmp,(ncoeff_max))
         do ii=1,ncoeff_max
           call polynomial_coeff_init(&
&           effective_potential%anharmonics_terms%coefficients(ii)%coefficient,&
&           effective_potential%anharmonics_terms%coefficients(ii)%nterm,&
&           coeffs_all(ii),&
&           effective_potential%anharmonics_terms%coefficients(ii)%terms,&
&           effective_potential%anharmonics_terms%coefficients(ii)%name,&
&           check=.false.)
         end do

         model_ncoeffbound = 0
         do ii=1,ncoeff_max-ncoeff

           write(message, '(5a,I0,a)')ch10,'--',ch10,' Try to bound the model ',&
&           'with ', ii,' additional term'
           if(ii>1)write(message,'(2a)') trim(message),'s'
           call wrtout(std_out,message,'COLL')
           call wrtout(ab_out,message,'COLL')

!          Copy the new model in coeffs_tmp(jj)
!          Free the coeffs_tmp array before
           do jj=1,ncoeff_max
             call polynomial_coeff_free(coeffs_tmp(jj))
           end do

           do jj=1,ncoeff+ii
             call polynomial_coeff_init(&
&             coeffs_all(jj)%coefficient,&
&             coeffs_all(jj)%nterm,&
&             coeffs_tmp(jj),&
&             coeffs_all(jj)%terms,&
&             coeffs_all(jj)%name,&
&             check=.false.)
           end do

           model_ncoeffbound = ii

!          Reset the simulation and set the coefficients of the model
           call effective_potential_setCoeffs(coeffs_tmp(1:ncoeff+ii),effective_potential,ncoeff+ii)
           call fit_polynomial_coeff_fit(effective_potential,(/0/),(/0/),hist,0,(/0,0/),0,0,&
&           -1,1,comm,verbose=.true.,positive=.false.)
           call effective_potential_setSupercell(effective_potential,comm,ncell=sc_size)
           dtset%rprimd_orig(:,:,1) = effective_potential%supercell%rprimd
           acell(:) = one
           call xcart2xred(dtset%natom,effective_potential%supercell%rprimd,&
&           effective_potential%supercell%xcart,xred)
           xred_old = xred
           vel_cell(:,:) = zero
           vel(:,:)      = zero
           fred(:,:)     = zero
           fcart(:,:)    = zero

!          Run mover to check if the model is bound
           call mover(scfcv_args,ab_xfh,acell,effective_potential%crystal%amu,dtfil,electronpositron,&
&           rhog,rhor,dtset%rprimd_orig,vel,vel_cell,xred,xred_old,&
&           effective_potential=effective_potential,verbose=verbose,writeHIST=writeHIST)
           if(.not.effective_potential%anharmonics_terms%bounded)then
             write(message, '(2a)' ) ' => The model is not bounded'
           else
             write(message, '(2a)' ) ' => The model is bounded'
           end if
           call wrtout(std_out,message,'COLL')

!          Exit if the model is bounded
           if(effective_potential%anharmonics_terms%bounded)  exit

         end do

       else
!      Get the list of possible coefficients to bound the model
         cutoff = zero
         do ii=1,3
           cutoff = cutoff + effective_potential%crystal%rprimd(ii,ii)
         end do
         cutoff = cutoff / 3.0

!        call fit_polynomial_coeff_getCoeffBound(effective_potential,coeffs_bound,&
!&                                              hist,ncoeff_bound,comm,verbose=.true.)

!TEST_AM!
         sc_size_TS = (/2,2,2/)
         call polynomial_coeff_getNorder(coeffs_bound,effective_potential%crystal,cutoff,&
&         ncoeff_bound,ncoeff_bound_tot,inp%bound_rangePower,inp%bound_rangePower(2),2,sc_size_TS,&
&         comm,anharmstr=inp%bound_anhaStrain==1,&
&         spcoupling=inp%bound_SPCoupling==1,verbose=.false.,distributed=.false.,&
&         only_even_power=.true.,only_odd_power=.false.)

         if(iam_master)then
           filename=trim(filnam(2))//"_boundcoeff.xml"
           call polynomial_coeff_writeXML(coeffs_bound,ncoeff_bound,filename=filename,newfile=.true.)
         end if
!      wait
         call xmpi_barrier(comm)
!      Store all the initial coefficients
         ncoeff = effective_potential%anharmonics_terms%ncoeff
         ABI_DATATYPE_ALLOCATE(coeffs_all,(ncoeff+ncoeff_bound))
         do ii=1,ncoeff
           call polynomial_coeff_init(effective_potential%anharmonics_terms%coefficients(ii)%coefficient,&
&           effective_potential%anharmonics_terms%coefficients(ii)%nterm,&
&           coeffs_all(ii),&
&           effective_potential%anharmonics_terms%coefficients(ii)%terms,&
&           effective_potential%anharmonics_terms%coefficients(ii)%name,&
&           check=.false.)
         end do

         do ii=1,ncoeff_bound
           call polynomial_coeff_init(coeffs_bound(ii)%coefficient,&
&           coeffs_bound(ii)%nterm,&
&           coeffs_all(ncoeff+ii),&
&           coeffs_bound(ii)%terms,&
&           coeffs_bound(ii)%name,&
&           check=.false.)
         end do

!      Copy the fixed coefficients from the model (without bound coeff)
         ncoeff = effective_potential%anharmonics_terms%ncoeff
         ABI_DATATYPE_ALLOCATE(coeffs_tmp,(ncoeff+ncoeff_bound))
         do ii=1,ncoeff
           call polynomial_coeff_init(effective_potential%anharmonics_terms%coefficients(ii)%coefficient,&
&           effective_potential%anharmonics_terms%coefficients(ii)%nterm,&
&           coeffs_tmp(ii),&
&           effective_potential%anharmonics_terms%coefficients(ii)%terms,&
&           effective_potential%anharmonics_terms%coefficients(ii)%name,&
&           check=.false.)
         end do

         ncoeff_max = ncoeff+ncoeff_bound
         ABI_ALLOCATE(listcoeff,(ncoeff_max))
         listcoeff = 0
         do jj=1,ncoeff
           listcoeff(jj) = jj
         end do

         model_bound = 0
         model_ncoeffbound = 0

         do ii=2,inp%bound_maxCoeff
!        Compute the number of possible combination
           nmodels = 1
           ABI_ALLOCATE(list_bound,(nmodels,ii))
           ABI_ALLOCATE(list_tmp,(ii))
           list_bound = 0; list_tmp = 0; kk = 0;  jj = 1

!        Generate the list of possible combinaison 1st count
           call genereList(kk,jj,ii,ncoeff_bound,list_tmp,list_bound,nmodels,.false.)
           nmodels = kk

           write(message, '(5a,I0,a,I0,a)')ch10,'--',ch10,' Try to bound the model ',&
&           'with ', ii,' additional positive terms (',nmodels,') possibilities'
           call wrtout(std_out,message,'COLL')

!        allocate and generate combinaisons
           ABI_DEALLOCATE(list_bound)
           ABI_DEALLOCATE(list_tmp)
           ABI_ALLOCATE(coeff_values,(nmodels,ncoeff+ii))
           ABI_ALLOCATE(listcoeff_bound,(nmodels,ncoeff+ii))
           ABI_ALLOCATE(list_bound,(nmodels,ii))
           ABI_ALLOCATE(list_tmp,(ii))
           ABI_ALLOCATE(isPositive,(nmodels))
           list_bound = 0;  listcoeff_bound = 0;  list_tmp = 0; isPositive = 0; kk = 0; jj = 1
           call genereList(kk,jj,ii,ncoeff_bound,list_tmp,list_bound,nmodels,.true.)
!        Generate the models
           do jj=1,nmodels
             listcoeff_bound(jj,1:ncoeff) = listcoeff(1:ncoeff)
             listcoeff_bound(jj,ncoeff+1:ncoeff+ii) = list_bound(jj,:) + ncoeff
           end do

!        Reset the simulation
           call effective_potential_setCoeffs(coeffs_all,effective_potential,ncoeff+ncoeff_bound)
           call fit_polynomial_coeff_getPositive(effective_potential,hist,coeff_values,&
&           isPositive,listcoeff_bound,ncoeff+ii,ncoeff,nmodels,comm,verbose=.false.)
           if(all(isPositive == 0)) then
             write(message, '(5a,I0,a)')ch10,'--',ch10,' No possible model ',&
&             'with ', ii,' additional terms found'
             call wrtout(std_out,message,'COLL')
           else

             do jj=1,nmodels
               if(isPositive(jj) == 1 .and. all(abs(coeff_values(jj,:)) < 1.0E5)) then
                 write(message, '(2a,I0,a)') ch10,' The model number ',jj,' ['
                 do kk=1,ncoeff+ii
                   if(kk<ncoeff+ii)then
                     write(message, '(a,I0,a)') trim(message),listcoeff_bound(jj,kk),','
                   else
                     write(message, '(a,I0)') trim(message),listcoeff_bound(jj,kk)
                   end if
                 end do
                 write(message, '(2a)') trim(message),'] is positive'
                 call wrtout(std_out,message,'COLL')
                 write(message, '(2a,I0,a)') ' Check if the model ',&
&                 'number ', jj,' is bounded...'
                 call wrtout(std_out,message,'COLL')

!             Set the coefficients of the model
                 do kk=1,ncoeff+ii
                   if(kk<=ncoeff)then
!                 just set the values of the coefficient
                     call polynomial_coeff_setCoefficient(coeff_values(jj,kk),coeffs_tmp(kk))
                     write(message, '(a,I0,a,ES19.10,2a)') ' Set the value of the coefficient ',kk,&
&                     ' =>',coeff_values(jj,kk),'     ',trim(coeffs_tmp(kk)%name)
                     call wrtout(std_out,message,'COLL')

                   else
!                 Set the good coefficient
                     icoeff_bound = listcoeff_bound(jj,kk)-ncoeff ! need to remove ncoeff value
                     write(message, '(a,I0,a,I0,a,ES19.10,2a)')&
&                     ' Set the value of the coefficient ',kk,' (',icoeff_bound,&
&                     ') =>',coeff_values(jj,kk),&
&                     '     ',trim(coeffs_bound(icoeff_bound)%name)
                     call wrtout(std_out,message,'COLL')
                     call polynomial_coeff_free(coeffs_tmp(kk))
                     call polynomial_coeff_init(coeff_values(jj,kk),&
&                     coeffs_bound(icoeff_bound)%nterm,&
&                     coeffs_tmp(kk),&
&                     coeffs_bound(icoeff_bound)%terms,&
&                     coeffs_bound(icoeff_bound)%name,&
&                     check=.false.)

                   end if
                 end do

!              Reset the simulation and set the coefficients of the model
                 call effective_potential_setCoeffs(coeffs_tmp(1:ncoeff+ii),effective_potential,&
&                 ncoeff+ii)
                 call fit_polynomial_coeff_fit(effective_potential,(/0/),(/0/),hist,0,(/0,0/),1,0,&
&                 -1,1,comm,verbose=.false.,positive=.false.)
                 call effective_potential_setSupercell(effective_potential,comm,ncell=sc_size)
                 dtset%rprimd_orig(:,:,1) = effective_potential%supercell%rprimd
                 acell(:) = one
                 call xcart2xred(dtset%natom,effective_potential%supercell%rprimd,&
&                 effective_potential%supercell%xcart,xred)
                 xred_old = xred
                 vel_cell(:,:) = zero
                 vel(:,:)      = zero
                 fred(:,:)     = zero
                 fcart(:,:)    = zero

!              Run mover
                 call mover(scfcv_args,ab_xfh,acell,effective_potential%crystal%amu,dtfil,electronpositron,&
&                 rhog,rhor,dtset%rprimd_orig,vel,vel_cell,xred,xred_old,&
&                 effective_potential=effective_potential,verbose=verbose,writeHIST=writeHIST)

                 if(.not.effective_potential%anharmonics_terms%bounded)then
                   write(message, '(2a)' ) ' => The model is not bounded'
                 else
                   write(message, '(2a)' ) ' => The model is bounded'
                 end if
                 call wrtout(std_out,message,'COLL')
!             Exit if the model is bounded
                 if(effective_potential%anharmonics_terms%bounded) then
                   model_bound = jj
                   model_ncoeffbound = ii
                   exit
                 end if
               end if
             end do
           end if

           ABI_DEALLOCATE(list_tmp)
           ABI_DEALLOCATE(list_bound)
           ABI_DEALLOCATE(isPositive)

!        Exit if the model is bounded
           if(effective_potential%anharmonics_terms%bounded) then
!         Final transfert
             write(message, '(3a)' ) ch10,' => The model is now bounded'
             call wrtout(ab_out,message,'COLL')
             call wrtout(std_out,message,'COLL')
             do kk=ncoeff+1,ncoeff+model_ncoeffbound
               icoeff_bound = listcoeff_bound(model_bound,kk)-ncoeff ! need to remove ncoeff value
               call polynomial_coeff_free(coeffs_tmp(kk))
               call polynomial_coeff_init(coeff_values(model_bound,kk),&
&               coeffs_bound(icoeff_bound)%nterm,&
&               coeffs_tmp(kk),&
&               coeffs_bound(icoeff_bound)%terms,&
&               coeffs_bound(icoeff_bound)%name,&
&               check=.false.)
             end do
             ABI_DEALLOCATE(coeff_values)
             ABI_DEALLOCATE(listcoeff_bound)
             exit
           end if
           ABI_DEALLOCATE(coeff_values)

         end do

         do ii=1,ncoeff_bound
           call polynomial_coeff_free(coeffs_bound(ii))
         end do
         if(allocated(coeffs_bound)) ABI_DEALLOCATE(coeffs_bound)

       end if

       if(.not.effective_potential%anharmonics_terms%bounded)then
         write(message, '(3a)' ) ch10,' => The model cannot be bounded'
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')
         model_ncoeffbound = 0
         model_bound = 0
       end if

!      Fit the final model
       call effective_potential_setCoeffs(coeffs_tmp(1:ncoeff+model_ncoeffbound),effective_potential,&
&       ncoeff+model_ncoeffbound)

       call fit_polynomial_coeff_fit(effective_potential,(/0/),(/0/),hist,0,(/0,0/),0,0,&
&       -1,1,comm,verbose=.false.)

       write(message, '(3a)') ch10,' Fitted coefficients at the end of the fit bound process: '
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')

       do ii = 1,ncoeff+model_ncoeffbound
         write(message, '(a,I0,a,ES19.10,2a)') " ",ii," =>",&
&         effective_potential%anharmonics_terms%coefficients(ii)%coefficient,&
&         " ",trim(effective_potential%anharmonics_terms%coefficients(ii)%name)
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')
       end do

!      Deallocation
       ABI_DEALLOCATE(listcoeff)
       do ii=1,ncoeff_max
         call polynomial_coeff_free(coeffs_tmp(ii))
       end do
       if(allocated(coeffs_tmp)) ABI_DEALLOCATE(coeffs_tmp)

       do ii=1,ncoeff_max
         call polynomial_coeff_free(coeffs_all(ii))
       end do
       if(allocated(coeffs_all)) ABI_DEALLOCATE(coeffs_all)

     end if

   else  if (option == -3) then

!*************************************************************
!   Call the routine for calculation of the energy for specific
!   partern of displacement or strain for the effective
!   Hamiltonian
!*************************************************************
!     write(message, '((80a),4a)' ) ('-',ii=1,80), ch10,&
!&     ' Effective Hamiltonian calculation'
!     call wrtout(ab_out,message,'COLL')
!     call wrtout(std_out,message,'COLL')

!     acell = one
!     call gstate(args_gs,acell,codvsn,cpui,dtfil,dtset,iexit,initialized,&
!&                mpi_enreg,npwtot,dtset%occ_orig,pawang,pawrad,pawtab,&
!&                psps,results_gs,dtset%rprimd_orig,scf_history,vel,vel_cell,wvl,xred)

   end if

!***************************************************************
! 5   Deallocation of array
!***************************************************************

   ABI_DEALLOCATE(fred)
   ABI_DEALLOCATE(fcart)
   ABI_DEALLOCATE(indsym)
   ABI_DEALLOCATE(rhog)
   ABI_DEALLOCATE(rhor)
   ABI_DEALLOCATE(vel)
   ABI_DEALLOCATE(xred)
   ABI_DEALLOCATE(xred_old)
   ABI_DEALLOCATE(ab_xfh%xfhist)

   ! if(option == -3)then
   !   call args_gs_free(args_gs)
   !   call psps_free(psps)
   !   do ii = 1,npsp
   !     call paw_setup_free(paw_setup(ii))
   !   end do
   !   ABI_DEALLOCATE(paw_setup)
   !   ABI_DEALLOCATE(ipsp2xml)
   !   ABI_DEALLOCATE(pspheads)
   !   call pawrad_free(pawrad)
   !   call pawtab_free(pawtab)
   !   ABI_DATATYPE_DEALLOCATE(pawrad)
   !   ABI_DATATYPE_DEALLOCATE(pawtab)
   !   ABI_DEALLOCATE(npwtot)
   ! end if
   call dtset%free()
   call destroy_results_gs(results_gs)
   call scfcv_destroy(scfcv_args)
   call destroy_mpi_enreg(mpi_enreg)

 end if

 write(message, '(a,(80a),a,a)' ) ch10,&
& ('=',ii=1,80),ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

end subroutine mover_effpot
!!***

end module m_mover_effpot
!!***
