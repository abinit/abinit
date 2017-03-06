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
!!      alloc_copy,copy_supercell,destroy_mpi_enreg,destroy_results_gs
!!      destroy_supercell,dtset_free,effective_potential_initmpi_supercell
!!      effective_potential_printsupercell,init_results_gs,init_supercell,mover
!!      scfcv_destroy,strain_apply,strain_get,strain_init,strain_print,wrtout
!!      xcart2xred,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mover_effpot(inp,filnam,effective_potential,comm)

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
 use m_dtset,  only : dtset_free
 use m_strain
 use m_effective_potential
 use m_effective_potential_file
 use m_multibinit_dataset
 use m_phonon_supercell
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
 use interfaces_28_numeric_noabirule
 use interfaces_41_geometry
 use interfaces_95_drive, except_this_one => mover_effpot
!End of the abilint section

implicit none

!Arguments --------------------------------
!scalar
 integer, intent(in) :: comm
!array
 type(multibinit_dataset_type),intent(in) :: inp
 type(effective_potential_type),intent(inout)  :: effective_potential
 character(len=fnlen),intent(in) :: filnam(15)
!Local variables-------------------------------
!scalar
 integer :: ii,jj,nproc,my_rank
 real(dp):: qmass,bmass
 logical :: iam_master
 integer, parameter:: master = 0
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
 type(strain_type) :: strain
!arrays
!no_abirules
 integer,pointer :: indsym(:,:,:)
 integer,allocatable :: symrel(:,:,:)
 real(dp) :: acell(3)
 real(dp),allocatable :: amass(:) 
 real(dp),pointer :: rhog(:,:),rhor(:,:)
 real(dp),allocatable :: tnons(:,:)
 real(dp),allocatable :: xred(:,:),xred_old(:,:),xcart(:,:)
 real(dp),allocatable :: fred(:,:),fcart(:,:)
 real(dp),allocatable :: vel(:,:)
 real(dp) :: vel_cell(3,3),rprimd(3,3)
 real(dp) :: mat_strain(3,3)
 type(supercell_type) :: supercell
!TEST_AM
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

!*******************************************************************
! 1 Generate supercell and print information
!*******************************************************************

 write(message, '(a,(80a),a)') ch10,&
& ('=',ii=1,80),ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

!if special structure is specified in the input,
!a new supercell is compute
 acell = one
 rprimd = effective_potential%crystal%rprimd

 if(all(inp%acell > one))then   
   acell = inp%acell
 end if
 if(any(inp%rprim > zero))then
   do jj=1,3
     rprimd(:,jj)=inp%rprim(:,jj)*acell(jj)
   end do
 end if

! check new rprimd
 if(all(rprimd(1,:)==zero).or.&
& all(rprimd(2,:)==zero).or.all(rprimd(3,:)==zero)) then
   write(message, '(3a)' )&
&   ' There is a problem with rprim',ch10,&
&   'Action: correct rprim'
   MSG_BUG(message)
 end if

!if rprim is different from initial structure, we print warning
 if(any(rprimd-effective_potential%crystal%rprimd>tol10))then
   write(message,'(a)')&
&   ' WARNING: the structure for the dynamics is different than initiale structure'
   call wrtout(std_out,message,"COLL")
   call wrtout(ab_out,message,"COLL")
 end if

 ABI_ALLOCATE(xred,(3,effective_potential%crystal%natom))
 ABI_ALLOCATE(xcart,(3,effective_potential%crystal%natom))

!convert new xcart
 call xcart2xred(effective_potential%crystal%natom,effective_potential%crystal%rprimd,&
& effective_potential%crystal%xcart,xred)
 call xred2xcart(effective_potential%crystal%natom, rprimd, xcart, xred)

 call init_supercell(effective_potential%crystal%natom, 0,&
& real(inp%n_cell,dp),&
& rprimd,&
& effective_potential%crystal%typat,&
& xcart,&
& supercell)

!Store the information of the supercell of the reference structure into effective potential
 call copy_supercell(supercell,effective_potential%supercell)
!Set new MPI for the new supercell 
 call effective_potential_initmpi_supercell(effective_potential,comm)
!Deallocation of useless array
 call destroy_supercell(supercell)   

 ABI_DEALLOCATE(xred)
 ABI_DEALLOCATE(xcart)

 call effective_potential_printSupercell(effective_potential)

 if(inp%dynamics==12.or.inp%dynamics==13) then
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
   dtset%dtion = inp%dtion  ! Delta Time for IONs
   dtset%diismemory = 8 ! Direct Inversion in the Iterative Subspace MEMORY
   dtset%friction = 0.0001 ! internal FRICTION coefficient
   dtset%goprecon = 0   ! Geometry Optimization PREconditioner equations
   dtset%ionmov = inp%dynamics  ! Number for the dynamic
   dtset%jellslab = 0   ! include a JELLium SLAB in the cell
   dtset%mdwall = 10000 ! Molecular Dynamics WALL location
   dtset%natom = effective_potential%supercell%natom_supercell
   dtset%ntypat = effective_potential%crystal%ntypat
   dtset%nconeq = 0     ! Number of CONstraint EQuations
   dtset%noseinert = 1.d-5 ! NOSE INERTia factor
   dtset%nnos = inp%nnos       ! Number of nose masses Characteristic
   dtset%ntime = inp%ntime  ! Number of TIME steps 
   dtset%nsym = 1       ! Number of SYMmetry operations
   dtset%prtxml = 0     ! print the xml
   dtset%optcell = inp%optcell    ! OPTimize the CELL shape and dimensions Characteristic
   dtset%restartxf = 0  ! RESTART from (X,F) history
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
   dtset%mdtemp(1) = inp%temperature   !Molecular Dynamics Temperatures 
   dtset%mdtemp(2) = inp%temperature   !Molecular Dynamics Temperatures 
   ABI_ALLOCATE(dtset%prtatlist,(dtset%natom)) !PRinT by ATom LIST of ATom
   dtset%prtatlist(:) = zero


!  Set the barostat and thermonstat if ionmov == 13
   if(dtset%ionmov == 13)then

     qmass = (abs(1+product(inp%strtarget(1:3)/3))*dtset%natom* kb_THzK * dtset%mdtemp(1)) / (0.1**2)
     bmass = (abs(1+product(inp%strtarget(1:3)/3))*dtset%natom* kb_THzK * dtset%mdtemp(1)) / (0.01**2)

     if(dtset%nnos==0) then
       dtset%nnos = 1
       ABI_ALLOCATE(dtset%qmass,(dtset%nnos))
       dtset%qmass(:)  = qmass
       write(message,'(3a,F12.1,a)')&
&      ' WARNING: nnos is set to zero in the input',ch10,&
&      '          value by default for qmass: ',dtset%qmass(:),ch10
       call wrtout(std_out,message,"COLL")
     else
       ABI_ALLOCATE(dtset%qmass,(dtset%nnos)) ! Q thermostat mass
        dtset%qmass(:) = inp%qmass(:)
     end if
     if (inp%bmass == zero) then
       dtset%bmass = bmass
       write(message,'(3a,F12.1,a)')&
&      ' WARNING: bmass is set to zero in the input',ch10,&
&       '          value by default for bmass: ',dtset%bmass,ch10
       call wrtout(std_out,message,"COLL")
     else
       dtset%bmass = inp%bmass  ! Barostat mass
     end if
   end if
   
   dtset%strtarget(1:6) = -1 * inp%strtarget(1:6) / 29421.033d0 ! STRess TARGET
   ABI_ALLOCATE(symrel,(3,3,dtset%nsym))
   symrel = one
   call alloc_copy(symrel,dtset%symrel)
   ABI_ALLOCATE(tnons,(3,dtset%nsym))
   tnons = zero
   dtset%tnons = tnons
   call alloc_copy(tnons,dtset%tnons)
   call alloc_copy(effective_potential%supercell%typat_supercell,dtset%typat)
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
&     effective_potential%crystal%amu(effective_potential%supercell%typat_supercell(jj))
   end do
!  Set the dffil structure
   dtfil%filnam_ds(1:2)=filnam(1:2)
   dtfil%filnam_ds(3)=""
   dtfil%filnam_ds(4)=filnam(2)
   
   nullify (electronpositron)
   ABI_ALLOCATE(rhog,(2,1))
   ABI_ALLOCATE(rhor,(2,1))
   
!  Initialize xf history (should be put in inwffil)
   ab_xfh%nxfh=0
   ab_xfh%mxfh=(ab_xfh%nxfh-dtset%restartxf+1)+dtset%ntime+5 
   ABI_ALLOCATE(ab_xfh%xfhist,(3,dtset%natom+4,2,ab_xfh%mxfh))

!***************************************************************
!2  initialization of the structure for the dynamics
!***************************************************************

   ABI_ALLOCATE(dtset%rprimd_orig,(3,3,1))
   dtset%rprimd_orig(:,:,1) = effective_potential%supercell%rprimd_supercell
   
   acell(1) = dtset%rprimd_orig(1,1,1)
   acell(2) = dtset%rprimd_orig(2,2,1)
   acell(3) = dtset%rprimd_orig(3,3,1)

   ABI_ALLOCATE(xred,(3,dtset%natom))
   ABI_ALLOCATE(xred_old,(3,dtset%natom))
   ABI_ALLOCATE(vel,(3,dtset%natom))
   ABI_ALLOCATE(fred,(3,dtset%natom))
   ABI_ALLOCATE(fcart,(3,dtset%natom))

! Fill the strain from input file
   call strain_init(strain)
   if (any(inp%strain /= zero)) then
     write(message,'(2a)') ch10, ' Strain is imposed during the simulation'
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')
!    convert strain into matrix
     mat_strain(1,1) = inp%strain(1); mat_strain(2,2) = inp%strain(2); mat_strain(3,3) = inp%strain(3)
     mat_strain(3,2) = half * inp%strain(4) ; mat_strain(2,3) = half * inp%strain(4) 
     mat_strain(3,1) = half * inp%strain(5) ; mat_strain(1,3) = half * inp%strain(5)
     mat_strain(1,2) = half * inp%strain(6) ; mat_strain(2,1) = half * inp%strain(6)
     call strain_get(strain,mat_delta = mat_strain)
     effective_potential%strain = strain
     effective_potential%has_strain = .FALSE.
     call strain_print(effective_potential%strain)
     call strain_apply(effective_potential%supercell%rprimd_supercell,dtset%rprimd_orig(:,:,1),&
&     effective_potential%strain)
   end if
   
   
   call xcart2xred(dtset%natom,effective_potential%supercell%rprimd_supercell,&
&   effective_potential%supercell%xcart_supercell,xred)

   xred_old = xred
   vel_cell(:,:) = zero
   vel(:,:)      = zero

!TEST_AM
!  Random initilisation of the velocitie and scale to the temperature 
!  with Maxwell-Boltzman distribution
!     do ia=1,dtset%natom
!       do mu=1,3
!         vel(mu,ia)=sqrt(kb_HaK*dtset%mdtemp(1)/amass(ia))*cos(two_pi*uniformrandom(rand_seed))
!         vel(mu,ia)=vel(mu,ia)*sqrt(-2._dp*log(uniformrandom(rand_seed)))
!       end do
!     end do

! ! !  Get rid of center-of-mass velocity
!     sum_mass=sum(amass(:))
!     do mu=1,3
!       mass_ia=sum(amass(:)*vel(mu,:))
!       vel(mu,:)=vel(mu,:)-mass_ia/sum_mass
!     end do

! ! !  Compute v2gauss
!     v2gauss = zero
!     do ia=1,dtset%natom
!       do mu=1,3
!         v2gauss=v2gauss+vel(mu,ia)*vel(mu,ia)*amass(ia)
!       end do
!     end do
!  !  Now rescale the velocities to give the exact temperature
!     rescale_vel=sqrt(3._dp*dtset%natom*kb_HaK*dtset%mdtemp(1)/v2gauss)
!     vel(:,:)=vel(:,:)*rescale_vel

   vel_cell(:,:) = zero
   vel(:,:)      = zero

!TEST_AM

!*********************************************************
!3   Call main routine for monte carlo / molecular dynamics
!*********************************************************
   write(message, '((80a),3a)' ) ('-',ii=1,80), ch10,&
&   '-Monte Carlo / Molecular Dynamics ',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   call mover(scfcv_args,ab_xfh,acell,amass,dtfil,electronpositron,&
&   rhog,rhor,dtset%rprimd_orig,vel,vel_cell,xred,xred_old,effective_potential)
   
!***************************************************************
! 4   Deallocation of array   
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

!    ABI_ALLOCATE(disp,(inp%ntime,3,supercell%natom_supercell))
!    ABI_ALLOCATE(disp_tmp,(inp%ntime,3*supercell%natom_supercell))
!    ABI_ALLOCATE(energy,(inp%ntime))
!    ABI_ALLOCATE(fcart,(3,supercell%natom_supercell))
!    ABI_ALLOCATE(fred,(3,supercell%natom_supercell))
! !  Read displacement
!    call effective_potential_file_readDisplacement(filnam(4),disp,inp%ntime,&
! &                              supercell%natom_supercell)

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
! !&                                       supercell%natom_supercell,&
! !&                                       supercell%rprimd_supercell,&
! !&                                       supercell%xcart_supercell,1,disp(1,:,:))
 
!    end do

!    ABI_DEALLOCATE(disp)
!    ABI_DEALLOCATE(disp_tmp)
!    ABI_DEALLOCATE(energy)
!    ABI_DEALLOCATE(fcart)
!    ABI_DEALLOCATE(fred)
 
!  else if(.false.) then

!    ABI_ALLOCATE(disp,(inp%ntime,3,effective_potential%supercell%natom_supercell))
!    ABI_ALLOCATE(disp_tmp,(inp%ntime,3*effective_potential%supercell%natom_supercell))

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
!         do jj=1,(effective_potential%supercell%natom_supercell)
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
!        read(unit=line,fmt=*) (disp_tmp(1,jj),jj=1,3*effective_potential%supercell%natom_supercell)
!        disp(1,:,:) = reshape(disp_tmp(1,:),(/3,effective_potential%supercell%natom_supercell/))
 
!         do jj=1,(effective_potential%supercell%natom_supercell)
!           do kk=1,3
! !              disp(1,kk,jj) = disp(1,kk,jj) *  effective_potential%supercell%rprimd_supercell(kk,kk)
!           end do
!         end do

!      else if (option==3)then
! ! third verison
!        write(std_out,*),"read displacement 3"
!        disp = zero
!        do jj=1,(effective_potential%supercell%natom_supercell)
!          do kk=1,3
!            if(effective_potential%supercell%typat_supercell(jj) == 1) then
!              disp(1,1,1) = 2*1*0.0005 * effective_potential%supercell%rprimd_supercell(1,1)
!            end if
!          end do
!        end do
!      else if(option==4)then
!        ABI_ALLOCATE(xred,(3,effective_potential%supercell%natom_supercell))
!        ABI_ALLOCATE(xcart,(3,effective_potential%supercell%natom_supercell))
 
!        !four version
!        write(std_out,*),"read position" 
!        do jj=1,(effective_potential%supercell%natom_supercell)
!          read(funit,'(a)',err=10,end=10) readline
!          line=adjustl(readline)
!          read(unit=line,fmt=*)  (xcart(kk,jj),kk=1,3) 
!        end do
!        close(funit)
!        do jj = 1, effective_potential%supercell%natom_supercell
!          do kk=1,3
!            disp(1,kk,jj) =  (xcart(kk,jj) - effective_potential%supercell%xcart_supercell(kk,jj)) / &
! &          effective_potential%supercell%rprimd_supercell(kk,kk) 
!          end do
!        end do
!        ABI_DEALLOCATE(xred)
!        ABI_DEALLOCATE(xcart)

!      end if
 
! !!!!!!!!!
! 10   continue

!      ABI_ALLOCATE(energy,(inp%ntime))
!      ABI_ALLOCATE(fcart,(3,effective_potential%supercell%natom_supercell))
!      ABI_ALLOCATE(fred,(3,effective_potential%supercell%natom_supercell))
 
!      write(111,*) disp(1,:,:)
!      write(std_out,*),"compute energy",ii
!      call ifc_contribution(effective_potential,disp(1,:,:),energy_harmonic,fcart)
!      write(112,*) ii , energy_harmonic
!      write(std_out,*),"harmonic energy :",ii,energy_harmonic
! !SHOULD ADD STRTEN
! !     call effective_potential_getHarmonicContributions(effective_potential,energy_harmonic,fcart,fred,&
! !&                                                      effective_potential%supercell%natom_supercell,&
! !&                                                      effective_potential%supercell%rprimd_supercell,&
! !&                                                      effective_potential%supercell%xcart_supercell,&
! !&                                                      1,disp(1,:,:))

!      write(std_out,*),"forces cart:",fcart(1,1)
!      write(std_out,*),"forces red :",fred(1,1)

!      ABI_DEALLOCATE(energy)
!      ABI_DEALLOCATE(fcart)
!      ABI_DEALLOCATE(fred)
 

! !xred_old is xcart!!
! !    xred_old = effective_potential%supercell%xcart_supercell + disp
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
