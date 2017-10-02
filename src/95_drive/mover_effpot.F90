!{\src2tex{textfont=tt}}
!!****f* ABINIT/mover_effpot
!! NAME
!! mover_effpot
!!
!! FUNCTION
!! this routine is driver for using mover with effective potential
!! 
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (AM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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
!!      fit_polynomial_coeff_fit,fit_polynomial_coeff_getnorder
!!      fit_polynomial_coeff_getpositive,generelist,init_results_gs,mover
!!      polynomial_coeff_free,polynomial_coeff_init
!!      polynomial_coeff_setcoefficient,polynomial_coeff_writexml,scfcv_destroy
!!      wrtout,xcart2xred,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mover_effpot(inp,filnam,effective_potential,option,comm,hist)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use defs_datatypes
 use m_errors
 use m_abimover
 use m_build_info
 use m_xmpi
 use m_abimover
 use m_phonons
 use m_strain
 use m_effective_potential_file
 use m_supercell
 use m_multibinit_dataset, only : multibinit_dataset_type
 use m_effective_potential,only : effective_potential_type
 use m_fit_polynomial_coeff, only : polynomial_coeff_writeXML
 use m_fit_polynomial_coeff, only : fit_polynomial_coeff_fit,genereList
 use m_fit_polynomial_coeff, only : fit_polynomial_coeff_getPositive
 use m_polynomial_coeff,only : polynomial_coeff_getNorder
 use m_dtset,  only : dtset_free
 use m_abihist, only : abihist
 use m_ifc
 use m_ewald
 use m_mpinfo,           only : init_mpi_enreg,destroy_mpi_enreg
 use m_copy            , only : alloc_copy 
 use m_electronpositron, only : electronpositron_type
 use m_scfcv,            only : scfcv_t, scfcv_run,scfcv_destroy
 use m_results_gs,       only : results_gs_type,init_results_gs,destroy_results_gs

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mover_effpot'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_95_drive, except_this_one => mover_effpot
!End of the abilint section

implicit none

!Arguments --------------------------------
!scalar
 integer, intent(in) :: option,comm
!array
 type(multibinit_dataset_type),intent(in) :: inp
 type(effective_potential_type),intent(inout)  :: effective_potential
 character(len=fnlen),intent(in) :: filnam(15)
 type(abihist),optional,intent(inout):: hist
!Local variables-------------------------------
!scalar
 integer :: icoeff_bound,ii,jj,kk,nproc,ncoeff,nmodels,ncoeff_bound,ncoeff_max
 integer :: model_bound,model_ncoeffbound,my_rank,option_fit = zero
 real(dp):: cutoff,freq_q,freq_b,qmass,bmass
 logical :: iam_master
 integer, parameter:: master = 0
 logical :: verbose,writeHIST
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
 type(electronpositron_type), pointer :: electronpositron
 integer,target :: zero_integer
 type(ab_xfh_type) :: ab_xfh
 type(results_gs_type),target :: results_gs
 type(pseudopotential_type),target :: psps
!arrays
!no_abirules
 integer :: sc_size(3)
 integer,pointer :: indsym(:,:,:)
 integer,allocatable :: listcoeff(:),listcoeff_bound(:,:),list_tmp(:),list_bound(:,:)
 integer,allocatable :: isPositive(:)
 integer,allocatable :: symrel(:,:,:)
 real(dp) :: acell(3)
 real(dp),allocatable :: amass(:) ,coeff_values(:,:)
 real(dp),pointer :: rhog(:,:),rhor(:,:)
 real(dp),allocatable :: tnons(:,:)
 real(dp),allocatable :: xred(:,:),xred_old(:,:),xcart(:,:)
 real(dp),allocatable :: fred(:,:),fcart(:,:)
 real(dp),allocatable :: vel(:,:)
 real(dp) :: vel_cell(3,3),rprimd(3,3)
 type(polynomial_coeff_type),dimension(:),allocatable :: coeffs_all,coeffs_tmp,coeffs_bound
 character(len=fnlen) :: filename
!TEST_AM
 !real(dp) :: tsec(2),tcpu,tcpui,twall,twalli
 !real(dp),allocatable :: energy(:)
 !integer :: option
 !integer :: funit = 1,ii,kk
 !real(dp):: energy_harmonic
 !real(dp):: ener1(inp%ntime,2),ener2(inp%ntime,2)
 !character (len=500000) :: line,readline
 !character(len=500) :: message
 !character(len=fnlen) :: filename,filename2
! real(dp),allocatable :: disp(:,:,:),disp_tmp(:,:)
!TEST_AM

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

 if(option == 1) then
   option_fit = 1
   sc_size(:) = inp%fit_boundCell
 else
   option_fit = 0
   sc_size(:) = inp%n_cell
 end if

 acell = one
 rprimd = effective_potential%crystal%rprimd

 ABI_ALLOCATE(xred,(3,effective_potential%crystal%natom))
 ABI_ALLOCATE(xcart,(3,effective_potential%crystal%natom))

!convert new xcart
 call xcart2xred(effective_potential%crystal%natom,effective_potential%crystal%rprimd,&
& effective_potential%crystal%xcart,xred)
 call xred2xcart(effective_potential%crystal%natom, rprimd, xcart, xred)
!Generate supercell for the simulation
 call effective_potential_setSupercell(effective_potential,comm,n_cell=sc_size)

 ABI_DEALLOCATE(xred)
 ABI_DEALLOCATE(xcart)

 if(option_fit==1.or.&
& inp%dynamics==12.or.&
& inp%dynamics==13.or.&
& inp%dynamics==24.or.&
& inp%dynamics==25)then

!***************************************************************
!1 Convert some parameters into the structures used by mover.F90
!***************************************************************

!Set mpi_eng
   mpi_enreg%comm_cell  = comm
   mpi_enreg%me = my_rank

!Set the fake abinit dataset 
!Scalar
   dtset%nctime = 0     ! NetCdf TIME between output of molecular dynamics informations 
   dtset%delayperm = 0  ! DELAY between trials to PERMUTE atoms
   dtset%dilatmx = 1.0  ! DILATation : MaXimal value
   dtset%diismemory = 8 ! Direct Inversion in the Iterative Subspace MEMORY
   dtset%friction = 0.0001 ! internal FRICTION coefficient
   dtset%goprecon = 0   ! Geometry Optimization PREconditioner equations
   dtset%jellslab = 0   ! include a JELLium SLAB in the cell
   dtset%mdwall = 10000 ! Molecular Dynamics WALL location
   dtset%natom = effective_potential%supercell%natom
   dtset%ntypat = effective_potential%crystal%ntypat
   dtset%nconeq = 0     ! Number of CONstraint EQuations
   dtset%noseinert = 1.d-5 ! NOSE INERTia factor
   dtset%nnos = inp%nnos   ! Number of nose masses Characteristic
   dtset%nsym = 1       ! Number of SYMmetry operations
   dtset%prtxml = 0     ! print the xml
   dtset%signperm = 1   ! SIGN of PERMutation potential      
   dtset%strprecon = 1  ! STRess PRECONditioner
   dtset%tolmxf = 2.0d-5
   dtset%tsmear = 0.009500446 !
   dtset%vis = 100      ! VIScosity
   dtset%usewvl = 0     !
   dtset%useylm = 0     !

!array
   ABI_ALLOCATE(dtset%iatfix,(3,dtset%natom)) ! Indices of AToms that are FIXed
   dtset%iatfix = zero
   dtset%goprecprm(:) = zero !Geometry Optimization PREconditioner PaRaMeters equations
   ABI_ALLOCATE(dtset%prtatlist,(dtset%natom)) !PRinT by ATom LIST of ATom
   dtset%prtatlist(:) = zero

   if(option_fit==0)then
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
   else
!    Set default for the fit
     verbose = .FALSE.
     writeHIST = .FALSE.
     dtset%restartxf = 0  ! RESTART from (X,F) history
     dtset%dtion = 100  ! Delta Time for IONs
     dtset%ionmov = 13  ! Number for the dynamic
     dtset%ntime = inp%fit_boundStep  ! Number of TIME steps 
     dtset%optcell = 2    ! OPTimize the CELL shape and dimensions Characteristic 
     dtset%mdtemp(1) = inp%fit_boundTemp   !Molecular Dynamics Temperatures 
     dtset%mdtemp(2) = inp%fit_boundTemp !Molecular Dynamics Temperatures 
     dtset%strtarget(1:6) = zero
   end if


!  Set the barostat and thermonstat if ionmov == 13
   if(dtset%ionmov == 13)then

!    Select frequency of the barostat as a function of temperature
!    For small temperature, we need huge barostat and inversely
     if(dtset%mdtemp(1) <= 10) then
       freq_q = 0.002
       freq_b = 0.0002
     else  if(dtset%mdtemp(1) <= 50) then
       freq_q = 0.02
       freq_b = 0.002 
     else  if(dtset%mdtemp(1) > 50.and.dtset%mdtemp(1) < 300) then
       freq_q = 0.1
       freq_b = 0.01 
     else
       freq_q = 0.2
       freq_b = 0.02
     end if

!TEST_AM
     freq_q = 0.1
     freq_b = 0.01
!TEST_AM
     
     qmass=(abs(1+product(inp%strtarget(1:3)/3))*dtset%natom* kb_THzK * dtset%mdtemp(1)) / (freq_q**2)
     bmass=(abs(1+product(inp%strtarget(1:3)/3))*dtset%natom* kb_THzK * dtset%mdtemp(1)) / (freq_b**2)

     if(dtset%nnos==0) then
       dtset%nnos = 1
       ABI_ALLOCATE(dtset%qmass,(dtset%nnos))
       dtset%qmass(:)  = qmass
       write(message,'(3a,F20.1,a)')&
&       ' WARNING: nnos is set to zero in the input',ch10,&
&       '          value by default for qmass: ',dtset%qmass(:),ch10
       if(verbose)call wrtout(std_out,message,"COLL")
     else
       ABI_ALLOCATE(dtset%qmass,(dtset%nnos)) ! Q thermostat mass
       dtset%qmass(:) = inp%qmass(:)
     end if
     if (inp%bmass == zero) then
       dtset%bmass = bmass
       write(message,'(3a,F20.4,a)')&
&       ' WARNING: bmass is set to zero in the input',ch10,&
&       '          value by default for bmass: ',dtset%bmass,ch10
       if(verbose)call wrtout(std_out,message,"COLL")
     else
       dtset%bmass = inp%bmass  ! Barostat mass
     end if
   end if
   
   ABI_ALLOCATE(symrel,(3,3,dtset%nsym))
   symrel = one
   call alloc_copy(symrel,dtset%symrel)
   ABI_ALLOCATE(tnons,(3,dtset%nsym))
   tnons = zero
   call alloc_copy(tnons,dtset%tnons)
   call alloc_copy(effective_potential%supercell%typat,dtset%typat)
   call alloc_copy(effective_potential%crystal%znucl,dtset%znucl)   

!  set psps 
   psps%useylm = dtset%useylm
   
!  initialisation of results_gs
   call init_results_gs(dtset%natom,1,results_gs)

!  Set the pointers of scfcv_args
   zero_integer = 0
   scfcv_args%dtset     => dtset
   ABI_ALLOCATE(indsym,(4,dtset%nsym,dtset%natom))
   indsym = zero
   scfcv_args%indsym => indsym
   scfcv_args%mpi_enreg => mpi_enreg
   scfcv_args%ndtpawuj  => zero_integer
   scfcv_args%results_gs => results_gs
   scfcv_args%psps => psps
!  Set other arguments of the mover.F90 routines

   ABI_ALLOCATE(amass,(dtset%natom))
!  Assign masses to each atom (for MD)
   do jj = 1,dtset%natom
     amass(jj)=amu_emass*&
&     effective_potential%crystal%amu(effective_potential%supercell%typat(jj))
   end do
!  Set the dffil structure
   dtfil%filnam_ds(1:2)=filnam(1:2)
   dtfil%filnam_ds(3)=""
   dtfil%filnam_ds(4)=filnam(2)
   
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

   ABI_ALLOCATE(dtset%rprimd_orig,(3,3,1))
   dtset%rprimd_orig(:,:,1) = effective_potential%supercell%rprimd
   
   acell(1) = dtset%rprimd_orig(1,1,1)
   acell(2) = dtset%rprimd_orig(2,2,1)
   acell(3) = dtset%rprimd_orig(3,3,1)

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
!    monte carlo / molecular dynamics
!*********************************************************
   if(option_fit==1)then
     write(message, '((80a),4a)' ) ('-',ii=1,80), ch10,&
&     ' Try to bound the model',ch10,' Check if the model is bounded or not'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')

!    Try the model             
     call mover(scfcv_args,ab_xfh,acell,amass,dtfil,electronpositron,&
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


!      Get the list of possible coefficients to bound the model
       cutoff = zero
       do ii=1,3
         cutoff = cutoff + effective_potential%crystal%rprimd(ii,ii)
       end do
       cutoff = cutoff / 3.0
       
      call polynomial_coeff_getNorder(coeffs_bound,effective_potential%crystal,cutoff,&
&                                     ncoeff_bound,inp%fit_boundPower,1,comm,&
&                                     anharmstr=inp%fit_anhaStrain==1,&
&                                     spcoupling=inp%fit_SPCoupling==1)

       if(iam_master)then
         filename=trim(filnam(2))//"_boundcoeff.xml"
         call polynomial_coeff_writeXML(coeffs_bound,ncoeff_bound,filename=filename,newfile=.true.)
       end if
       
!     Store all the initial coefficients
       ncoeff = effective_potential%anharmonics_terms%ncoeff
       ABI_DATATYPE_ALLOCATE(coeffs_all,(ncoeff+ncoeff_bound))
       do ii=1,ncoeff
         call polynomial_coeff_init(effective_potential%anharmonics_terms%coefficients(ii)%coefficient,&
&         effective_potential%anharmonics_terms%coefficients(ii)%nterm,&
&         coeffs_all(ii),&
&         effective_potential%anharmonics_terms%coefficients(ii)%terms,&
&         effective_potential%anharmonics_terms%coefficients(ii)%name,&
&         check=.false.) 
       end do
       
       do ii=1,ncoeff_bound
         call polynomial_coeff_init(coeffs_bound(ii)%coefficient,&
&         coeffs_bound(ii)%nterm,&
&         coeffs_all(ncoeff+ii),&
&         coeffs_bound(ii)%terms,&
&         coeffs_bound(ii)%name,&
&         check=.false.)
       end do

!     Copy the fixed coefficients from the model (without bound coeff)
       ncoeff = effective_potential%anharmonics_terms%ncoeff
       ABI_DATATYPE_ALLOCATE(coeffs_tmp,(ncoeff+ncoeff_bound))
       do ii=1,ncoeff
         call polynomial_coeff_init(effective_potential%anharmonics_terms%coefficients(ii)%coefficient,&
&         effective_potential%anharmonics_terms%coefficients(ii)%nterm,&
&         coeffs_tmp(ii),&
&         effective_potential%anharmonics_terms%coefficients(ii)%terms,&
&         effective_potential%anharmonics_terms%coefficients(ii)%name,&
&         check=.false.) 
       end do
       
       ncoeff_max = ncoeff+ncoeff_bound              
       ABI_ALLOCATE(listcoeff,(ncoeff_max))
       listcoeff = zero
       do jj=1,ncoeff
         listcoeff(jj) = jj
       end do

       model_bound = zero
       model_ncoeffbound = zero
       
       do ii=3,inp%fit_boundTerm
!       Compute the number of possible combination         
         nmodels = factorial(ncoeff_bound) / (factorial(ii)*factorial(ncoeff_bound-ii))

         write(message, '(5a,I0,a,I0,a)')ch10,'--',ch10,' Try to bound the model ',&
&         'with ', ii,' additional positive terms (',nmodels,') possibilities'
         call wrtout(std_out,message,'COLL')

         ABI_ALLOCATE(coeff_values,(nmodels,ncoeff+ii))         
         ABI_ALLOCATE(listcoeff_bound,(nmodels,ncoeff+ii))
         ABI_ALLOCATE(list_bound,(nmodels,ii))
         ABI_ALLOCATE(list_tmp,(ii))
         ABI_ALLOCATE(isPositive,(nmodels))
         list_bound = zero
         listcoeff_bound = zero
         list_tmp = zero
         isPositive = zero
         kk = one
         jj = one
         
!       Generate the list of possible combinaison
         call genereList(kk,jj,ii,ncoeff_bound,list_tmp,list_bound,nmodels)
         do jj=1,nmodels
           listcoeff_bound(jj,1:ncoeff) = listcoeff(1:ncoeff)
           listcoeff_bound(jj,ncoeff+1:ncoeff+ii) = list_bound(jj,:) + ncoeff
         end do
         
!       Reset the simulation
         call effective_potential_setCoeffs(coeffs_all,effective_potential,ncoeff+ncoeff_bound)

         call fit_polynomial_coeff_getPositive(effective_potential,hist,coeff_values,&
&         isPositive,listcoeff_bound,ncoeff+ii,&
&         ncoeff,nmodels,comm,verbose=.false.)
         if(all(isPositive == zero)) then
           write(message, '(5a,I0,a)')ch10,'--',ch10,' No possible model ',&
&           'with ', ii,' additional terms found'
           call wrtout(std_out,message,'COLL')
         else

           do jj=1,nmodels
             if(isPositive(jj) == one.and.all(coeff_values(jj,:) < 10)) then
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
&               'number ', jj,' is bounded...'
               call wrtout(std_out,message,'COLL')

!             Set the coefficients of the model
               do kk=1,ncoeff+ii
                 if(kk<=ncoeff)then
!                 just set the values of the coefficient
                   call polynomial_coeff_setCoefficient(coeff_values(jj,kk),coeffs_tmp(kk))
                   write(message, '(a,I0,a,ES19.10,2a)') ' Set the value of the coefficient ',kk,&
&                   ' =>',coeff_values(jj,kk),'     ',trim(coeffs_tmp(kk)%name)
                   call wrtout(std_out,message,'COLL')

                 else
!                 Set the good coefficient
                   icoeff_bound = listcoeff_bound(jj,kk)-ncoeff ! need to remove ncoeff value
                   write(message, '(a,I0,a,I0,a,ES19.10,2a)')&
&                   ' Set the value of the coefficient ',kk,' (',icoeff_bound,&
&                   ') =>',coeff_values(jj,kk),&
&                   '     ',trim(coeffs_bound(icoeff_bound)%name)
                   call wrtout(std_out,message,'COLL')
                   call polynomial_coeff_free(coeffs_tmp(kk))
                   call polynomial_coeff_init(coeff_values(jj,kk),&
&                   coeffs_bound(icoeff_bound)%nterm,&
&                   coeffs_tmp(kk),&
&                   coeffs_bound(icoeff_bound)%terms,&
&                   coeffs_bound(icoeff_bound)%name,&
&                   check=.false.)
                   
                 end if
               end do

!             Reset the simulation and set the coefficients of the model 
               call effective_potential_setCoeffs(coeffs_tmp(1:ncoeff+ii),effective_potential,&
&               ncoeff+ii)
               call fit_polynomial_coeff_fit(effective_potential,&
&               (/0/),(/0/),hist,(/0,0/),1,0,&
&               -1,comm,verbose=.false.,positive=.false.) 
               call effective_potential_setSupercell(effective_potential,comm,n_cell=sc_size)
               dtset%rprimd_orig(:,:,1) = effective_potential%supercell%rprimd
               acell(1) = dtset%rprimd_orig(1,1,1)
               acell(2) = dtset%rprimd_orig(2,2,1)
               acell(3) = dtset%rprimd_orig(3,3,1)
               call xcart2xred(dtset%natom,effective_potential%supercell%rprimd,&
&               effective_potential%supercell%xcart,xred)
               xred_old = xred
               vel_cell(:,:) = zero
               vel(:,:)      = zero
               fred(:,:)     = zero
               fcart(:,:)    = zero
               
!             Run mover
               call mover(scfcv_args,ab_xfh,acell,amass,dtfil,electronpositron,&
&               rhog,rhor,dtset%rprimd_orig,vel,vel_cell,xred,xred_old,&
&               effective_potential=effective_potential,verbose=verbose,writeHIST=.false.)

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
         
!       Exit if the model is bounded
         if(effective_potential%anharmonics_terms%bounded) then
!         Final transfert
           write(message, '(3a)' ) ch10,' => The model is now bounded'
           call wrtout(ab_out,message,'COLL')
           call wrtout(std_out,message,'COLL')
           do kk=ncoeff+1,ncoeff+model_ncoeffbound
             icoeff_bound = listcoeff_bound(model_bound,kk)-ncoeff ! need to remove ncoeff value
             call polynomial_coeff_free(coeffs_tmp(kk))
             call polynomial_coeff_init(coeff_values(model_bound,kk),&
&             coeffs_bound(icoeff_bound)%nterm,&
&             coeffs_tmp(kk),&
&             coeffs_bound(icoeff_bound)%terms,&
&             coeffs_bound(icoeff_bound)%name,&
&             check=.false.)
           end do
           ABI_DEALLOCATE(coeff_values)
           ABI_DEALLOCATE(listcoeff_bound)
           exit
         end if
         ABI_DEALLOCATE(coeff_values)
         ABI_DEALLOCATE(listcoeff_bound)
       end do

       if(.not.effective_potential%anharmonics_terms%bounded)then
         write(message, '(3a)' ) ch10,' => The model cannot be bounded'
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')
         model_ncoeffbound = zero
         model_bound = zero
       end if
       
!     Fit the final model
       call effective_potential_setCoeffs(coeffs_tmp(1:ncoeff+model_ncoeffbound),effective_potential,&
&       ncoeff+model_ncoeffbound)

       call fit_polynomial_coeff_fit(effective_potential,&
&       (/0/),(/0/),hist,(/0,0/),1,0,&
&       -1,comm,verbose=.false.,positive=.false.,anharmstr=.false.)
       
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

!     Deallocation
       ABI_DEALLOCATE(listcoeff)
       do ii=1,ncoeff+ncoeff_bound
         call polynomial_coeff_free(coeffs_tmp(ii))
       end do
       if(allocated(coeffs_tmp)) ABI_DEALLOCATE(coeffs_tmp)
       
       do ii=1,ncoeff_bound
         call polynomial_coeff_free(coeffs_bound(ii))
       end do
       if(allocated(coeffs_bound)) ABI_DEALLOCATE(coeffs_bound)
       
       do ii=1,ncoeff+ncoeff_bound
         call polynomial_coeff_free(coeffs_all(ii))
       end do
       if(allocated(coeffs_all)) ABI_DEALLOCATE(coeffs_all)

     end if
     
   else
!    just call mover in case of NPT or NVT simulation
     write(message, '((80a),3a)' ) ('-',ii=1,80), ch10,&
&     '-Monte Carlo / Molecular Dynamics ',ch10
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     call mover(scfcv_args,ab_xfh,acell,amass,dtfil,electronpositron,&
&     rhog,rhor,dtset%rprimd_orig,vel,vel_cell,xred,xred_old,&
&     effective_potential=effective_potential,verbose=verbose,writeHIST=writeHIST)
   end if

!***************************************************************
! 5   Deallocation of array   
!***************************************************************

   ABI_DEALLOCATE(amass)
   ABI_DEALLOCATE(fred)
   ABI_DEALLOCATE(fcart)
   ABI_DEALLOCATE(indsym)
   ABI_DEALLOCATE(rhog)
   ABI_DEALLOCATE(rhor)
   ABI_DEALLOCATE(symrel)
   ABI_DEALLOCATE(tnons)
   ABI_DEALLOCATE(vel)
   ABI_DEALLOCATE(xred)
   ABI_DEALLOCATE(xred_old)
   ABI_DEALLOCATE(ab_xfh%xfhist)

   call dtset_free(dtset)
   call destroy_results_gs(results_gs)
   call scfcv_destroy(scfcv_args)
   call destroy_mpi_enreg(mpi_enreg)
 end if

 write(message, '(a,(80a),a,a)' ) ch10,&
& ('=',ii=1,80),ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')


!  if(inp%dynamics == 3) then
! !*************************************************************
! !   Call the routine for calculation of the energy for specific 
! !   partern of displacement or strain for the effective 
! !   Hamiltonian
! !*************************************************************
!    write(message, '(3a,(80a),a)' ) &
! &    ch10,'Read displacement for effective hamiltonian ',ch10,&
! &    ('-',ii=1,80),ch10
!    call wrtout(ab_out,message,'COLL')
!    call wrtout(std_out,message,'COLL')

!    ABI_ALLOCATE(disp,(inp%ntime,3,supercell%natom))
!    ABI_ALLOCATE(disp_tmp,(inp%ntime,3*supercell%natom))
!    ABI_ALLOCATE(energy,(inp%ntime))
!    ABI_ALLOCATE(fcart,(3,supercell%natom))
!    ABI_ALLOCATE(fred,(3,supercell%natom))
! !  Read displacement
!    call effective_potential_file_readDisplacement(filnam(4),disp,inp%ntime,&
! &                              supercell%natom)

!    do ii=1,inp%ntime
!      write(message, '(a,(80a),2a,I3)' ) ch10,&
! &    ('-',jj=1,80),ch10,' Displacements ',ii
!      call wrtout(ab_out,message,'COLL')
!      call wrtout(std_out,message,'COLL')

!      write(111,*) disp(ii,:,:)
!      call ifc_contribution(effective_potential,disp(ii,:,:),energy(ii),fcart)
!      write(112,'(I3,es23.14)') ii , energy(ii)
! !SHOULD ADD STRTEN
! !     call effective_potential_getHarmonicContributions(effective_potential,energy_harmonic,fcart,fred,&
! !&                                       supercell%natom,&
! !&                                       supercell%rprimd,&
! !&                                       supercell%xcart,1,disp(1,:,:))
 
!    end do

!    ABI_DEALLOCATE(disp)
!    ABI_DEALLOCATE(disp_tmp)
!    ABI_DEALLOCATE(energy)
!    ABI_DEALLOCATE(fcart)
!    ABI_DEALLOCATE(fred)
 
!  else if(.false.) then

!    ABI_ALLOCATE(disp,(inp%ntime,3,effective_potential%supercell%natom))
!    ABI_ALLOCATE(disp_tmp,(inp%ntime,3*effective_potential%supercell%natom))

! !filename="/home/alex/Desktop/dev/test/BaTiO3/heff/config/cfg_6.dat"
! !    filename="/home/alex/Desktop/dev/test/CaTiO3/script/disp.dat"
!    filename="/home/alex/Desktop/dev/test/CaTiO3/spld_output/result_112/fort.111"
! !   filename="/home/alex/Desktop/dev/test/CaTiO3/spld_output/result_112/xcart_9.dat"
! !   filename="/home/alex/Desktop/dev/test/CaTiO3/spld_output/result_112/fort.111"
! !   filename="/home/alex/Desktop/dev/test/CaTiO3/spld_output/result_222/fort.111"
! !    filename="/home/alex/Desktop/dev/test/CaTiO3/spld_output/CaTiO3_NC.disp"


!    option = 2

!    do ii=1,inp%ntime
!      write(std_out,*),"Step :", ii
!      disp = zero;
! !  First version:
!  ! end if
!       if (open_file(filename,message,unit=funit,form="formatted",&
!                  status="old",action="read") /= 0) then
!         MSG_ERROR(message)
!       end if


!       if(option ==1)then
!         write(std_out,*),"read displacement "
!         do jj=1,(effective_potential%supercell%natom)
! !          do kk=1,3
!             read(funit,'(a)',err=10,end=10) readline
!             line=adjustl(readline)
!             read(unit=line,fmt=*) (disp(1,kk,jj),kk=1,3)
!             write(667,*),disp(1,:,jj)
! !         end do
!        end do
!      else if (option ==2) then
! ! !second version: (fort.111)
!        write(std_out,*),"read displacement 2"
!        read(funit,'(a)',err=10,end=10) readline
!        line=adjustl(readline)  
!        read(unit=line,fmt=*) (disp_tmp(1,jj),jj=1,3*effective_potential%supercell%natom)
!        disp(1,:,:) = reshape(disp_tmp(1,:),(/3,effective_potential%supercell%natom/))
 
!         do jj=1,(effective_potential%supercell%natom)
!           do kk=1,3
! !              disp(1,kk,jj) = disp(1,kk,jj) *  effective_potential%supercell%rprimd(kk,kk)
!           end do
!         end do

!      else if (option==3)then
! ! third verison
!        write(std_out,*),"read displacement 3"
!        disp = zero
!        do jj=1,(effective_potential%supercell%natom)
!          do kk=1,3
!            if(effective_potential%supercell%typat(jj) == 1) then
!              disp(1,1,1) = 2*1*0.0005 * effective_potential%supercell%rprimd(1,1)
!            end if
!          end do
!        end do
!      else if(option==4)then
!        ABI_ALLOCATE(xred,(3,effective_potential%supercell%natom))
!        ABI_ALLOCATE(xcart,(3,effective_potential%supercell%natom))
 
!        !four version
!        write(std_out,*),"read position" 
!        do jj=1,(effective_potential%supercell%natom)
!          read(funit,'(a)',err=10,end=10) readline
!          line=adjustl(readline)
!          read(unit=line,fmt=*)  (xcart(kk,jj),kk=1,3) 
!        end do
!        close(funit)
!        do jj = 1, effective_potential%supercell%natom
!          do kk=1,3
!            disp(1,kk,jj) =  (xcart(kk,jj) - effective_potential%supercell%xcart(kk,jj)) / &
! &          effective_potential%supercell%rprimd(kk,kk) 
!          end do
!        end do
!        ABI_DEALLOCATE(xred)
!        ABI_DEALLOCATE(xcart)

!      end if
 
! !!!!!!!!!
! 10   continue

!      ABI_ALLOCATE(energy,(inp%ntime))
!      ABI_ALLOCATE(fcart,(3,effective_potential%supercell%natom))
!      ABI_ALLOCATE(fred,(3,effective_potential%supercell%natom))
 
!      write(111,*) disp(1,:,:)
!      write(std_out,*),"compute energy",ii
!      call ifc_contribution(effective_potential,disp(1,:,:),energy_harmonic,fcart)
!      write(112,*) ii , energy_harmonic
!      write(std_out,*),"harmonic energy :",ii,energy_harmonic
! !SHOULD ADD STRTEN
! !     call effective_potential_getHarmonicContributions(effective_potential,energy_harmonic,fcart,fred,&
! !&                                                      effective_potential%supercell%natom,&
! !&                                                      effective_potential%supercell%rprimd,&
! !&                                                      effective_potential%supercell%xcart,&
! !&                                                      1,disp(1,:,:))

!      write(std_out,*),"forces cart:",fcart(1,1)
!      write(std_out,*),"forces red :",fred(1,1)

!      ABI_DEALLOCATE(energy)
!      ABI_DEALLOCATE(fcart)
!      ABI_DEALLOCATE(fred)
 

! !xred_old is xcart!!
! !    xred_old = effective_potential%supercell%xcart + disp
! !    call xcart2xred(dtset%natom,rprimd,xred_old,xred)
! !call effective_potential_getEnergy(effective_potential,energy_harmonic,dtset%natom,rprimd,xred,comm)
!      write(std_out,*),"done"
! !end version
! !  call effective_potential_getEnergy(effective_potential,energy_harmonic,dtset%natom,rprimd,xred,&
! !&                                   comm,displacement=disp)

!    end do
 
!    close(112)


!    if(.false.)then
! !  check the calculation
!      write(std_out,*) "final check:"
!      filename="fort.112"
!      filename2="/home/alex/Desktop/dev/test/CaTiO3/spld_output/result_112/fort.112"
 
!      if (open_file(filename,message,unit=1,form="formatted",&
!        status="old",action="read") /= 0) then
!        MSG_ERROR(message)
!      end if
 
!      if (open_file(filename2,message,unit=2,form="formatted",&
!        status="old",action="read") /= 0) then
!        MSG_ERROR(message)
!      end if
 
!      do ii=1,inp%ntime
!        read(1,'(a)',err=20,end=20) readline
!        line=adjustl(readline)
!        read(unit=line,fmt=*)  (ener1(ii,kk),kk=1,2) 
 
!        read(2,'(a)',err=20,end=20) readline
!        line=adjustl(readline)
!        read(unit=line,fmt=*)  (ener2(ii,kk),kk=1,2) 
!        if(abs(100*(ener1(ii,2)-ener2(ii,2)) / ener1(ii,2)) > tol5)then
!          write(std_out,*) "step", ii,abs(100*(ener1(ii,2)-ener2(ii,2)) / ener1(ii,2))
!        end if
!      end do
! 20   continue 
!    end if
!  else 

end subroutine mover_effpot
!!***

