!{\src2tex{textfont=tt}}
!!****f* ABINIT/mover_effpot
!! NAME
!! mover_effpot
!!
!! FUNCTION
!! this routine is driver for using mover with effective potential
!! 
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (AM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  inp = input of epigene
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
!!      epigine
!!
!! CHILDREN
!!      mover
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
 use m_errors
 use defs_datatypes
 use m_abimover
 use m_build_info
 use m_xmpi
 use m_abimover
 use m_phonons
 use m_effective_potential
 use m_effective_potential_file
 use m_epigene_dataset
 use m_phonon_supercell

!TEST
 use m_ifc
 use m_ewald

 use m_electronpositron, only : electronpositron_type
 use m_scfcv,            only : scfcv_t, scfcv_run
 use m_results_gs,       only : results_gs_type,init_results_gs

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mover_effpot'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_95_drive, except_this_one => mover_effpot
!End of the abilint section

implicit none

!Arguments --------------------------------
!scalar
 integer, intent(inout) :: comm
!array
 type(epigene_dataset_type),intent(in) :: inp
 type(effective_potential_type),intent(inout)  :: effective_potential
 character(len=fnlen),intent(in) :: filnam(15)
!Local variables-------------------------------
!scalar
 integer :: jj
 integer :: nproc,my_rank
 logical :: iam_master
 integer, parameter:: master = 0
! Set array dimensions
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
 real(dp) :: acell(3)
 real(dp),allocatable :: energy(:)
 real(dp),allocatable :: amass(:) 
 real(dp),pointer :: rhog(:,:),rhor(:,:)
 real(dp),allocatable :: xred(:,:),xred_old(:,:),xcart(:,:)
 real(dp),allocatable :: fred(:,:),fcart(:,:)
 real(dp),allocatable :: vel(:,:)
 real(dp) ::gprimd(3,3),vel_cell(3,3),rprimd(3,3)
!TEST_AM
 integer :: option
 integer :: funit = 1,ii,kk
 real(dp):: energy_harmonic
 character (len=500000) :: line,readline
 character(len=500) :: message
 character(len=fnlen) :: filename
 real(dp),allocatable :: disp(:,:,:),disp_tmp(:,:)
!TEST_AM

!******************************************************************

!MPI variables
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

!*********************************
! 1 print  supercell informatoins
!*********************************
 write(message, '(a,(80a),a)'),ch10,&
&    ('=',ii=1,80),ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 call effective_potential_printSupercell(effective_potential)

 if(inp%monte_carlo == 2) then
!*************************************************************
!   Call the routine for calculation of the energy for specific 
!   partern of displacement or strain for the effective 
!   Hamiltonian
!*************************************************************
   write(message, '(3a,(80a),a)' ) &
&    ch10,'Read displacement for effective hamiltonian ',ch10,&
&    ('-',ii=1,80),ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   ABI_ALLOCATE(disp,(inp%ntime,3,effective_potential%supercell%natom_supercell))
   ABI_ALLOCATE(disp_tmp,(inp%ntime,3*effective_potential%supercell%natom_supercell))
   ABI_ALLOCATE(energy,(inp%ntime))
!  Read displacement
   call effective_potential_file_readDisplacement(filnam(4),disp,inp%ntime,&
&                              effective_potential%supercell%natom_supercell)

   do ii=1,inp%ntime
     write(message, '(a,(80a),2a,I3)' ) ch10,&
&    ('-',jj=1,80),ch10,' Displacements ',ii
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')

     write(111,*) disp(ii,:,:)
     energy(ii) = harmonic_energy(effective_potential,disp(ii,:,:))
     write(112,'(I3,es23.14)') ii , energy(ii)
     call effective_potential_getEnergy(effective_potential,energy_harmonic,&
&                                        effective_potential%supercell%natom_supercell,&
&                                        effective_potential%supercell%rprimd_supercell,&
&                                        xred,1,disp(ii,:,:))
     
   end do
   
   ABI_DEALLOCATE(disp)
   ABI_DEALLOCATE(disp_tmp)
   ABI_DEALLOCATE(energy)

 else if(.false.) then

   ABI_ALLOCATE(disp,(inp%ntime,3,effective_potential%supercell%natom_supercell))
   ABI_ALLOCATE(disp_tmp,(inp%ntime,3*effective_potential%supercell%natom_supercell))

!filename="/home/alex/Desktop/dev/test/BaTiO3/heff/config/cfg_6.dat"
!    filename="/home/alex/Desktop/dev/test/CaTiO3/script/disp.dat"
   filename="/home/alex/Desktop/dev/test/CaTiO3/spld_output/result_112/fort.111"
!   filename="/home/alex/Desktop/dev/test/CaTiO3/spld_output/result_112/xcart_9.dat"
!   filename="/home/alex/Desktop/dev/test/CaTiO3/spld_output/fort.111"
!   filename="/home/alex/Desktop/dev/test/CaTiO3/spld_output/result_222/fort.111"
!    filename="/home/alex/Desktop/dev/test/CaTiO3/spld_output/CaTiO3_NC.disp"


   option = 2

   do ii=1,inp%ntime
     print*,"Step :", ii
     disp = zero;
!  First version:
 ! end if
      if (open_file(filename,message,unit=funit,form="formatted",&
                 status="old",action="read") /= 0) then
        MSG_ERROR(message)
      end if


      if(option ==1)then
        print*,"read displacement "
        do jj=1,(effective_potential%supercell%natom_supercell)
!          do kk=1,3
            read(funit,'(a)',err=10,end=10) readline
            line=adjustl(readline)
            read(unit=line,fmt=*) (disp(1,kk,jj),kk=1,3)
            write(667,*),disp(1,:,jj)
!         end do
       end do
       close(funit)
     else if (option ==2) then
! !second version: (fort.111)
       print*,"read displacement 2"
       read(funit,'(a)',err=10,end=10) readline
       line=adjustl(readline)  
       read(unit=line,fmt=*) (disp_tmp(1,jj),jj=1,3*effective_potential%supercell%natom_supercell)
       disp(1,:,:) = reshape(disp_tmp(1,:),(/3,effective_potential%supercell%natom_supercell/))
       
       ! do jj=1,(effective_potential%supercell%natom_supercell)
       !   do kk=1,3
       !       disp(1,kk,jj) = disp(1,kk,jj) *  effective_potential%supercell%rprimd_supercell(kk,kk)
       !   end do
       ! end do
     else if (option==3)then
! third verison
       print*,"read displacement 3"
       disp = zero
       do jj=1,(effective_potential%supercell%natom_supercell)
         do kk=1,3
           if(effective_potential%supercell%typat_supercell(jj) == 1) then
             disp(1,1,1) = 2*1*0.0005 * effective_potential%supercell%rprimd_supercell(1,1)
           end if
         end do
       end do
     else if(option==4)then
       !four version
       print*,"read position" 
       do jj=1,(effective_potential%supercell%natom_supercell)
         read(funit,'(a)',err=10,end=10) readline
         line=adjustl(readline)
         read(unit=line,fmt=*)  (xcart(kk,jj),kk=1,3) 
         print*,"tata",xcart(:,jj)
         write(456,*) effective_potential%supercell%xcart_supercell(:,jj)
       end do
       close(funit)
       do jj = 1, effective_potential%supercell%natom_supercell
         do kk=1,3
           disp(1,kk,jj) = (xcart(kk,jj) - effective_potential%supercell%xcart_supercell(kk,jj))! / &
!&          effective_potential%supercell%rprimd_supercell(kk,kk) 
         end do
       end do
     end if
      
!!!!!!!!!
10   continue

     write(111,*) disp(1,:,:)
     print*,"compute energy",ii
     energy_harmonic = harmonic_energy(effective_potential,disp(1,:,:))
     write(112,*) ii , energy_harmonic
     print*,"harmonic energy :",ii,energy_harmonic
     call effective_potential_getEnergy(effective_potential,energy_harmonic,&
&                                       effective_potential%supercell%natom_supercell,&
&                                       effective_potential%supercell%rprimd_supercell,&
&                                       effective_potential%supercell%xcart_supercell,1,disp(1,:,:))

     call effective_potential_getForces(effective_potential,fcart,fred,&
&                                       effective_potential%supercell%natom_supercell,&
&                                       effective_potential%supercell%rprimd_supercell,&
&                                       effective_potential%supercell%xcart_supercell,1,disp(1,:,:))
     print*,"forces cart:",fcart(1,1)
     print*,"forces red :",fred(1,1)

!xred_old is xcart!!
!    xred_old = effective_potential%supercell%xcart_supercell + disp
!    call xcart2xred(dtset%natom,rprimd,xred_old,xred)
!call effective_potential_getEnergy(effective_potential,energy_harmonic,dtset%natom,rprimd,xred,comm)
     print*,"done"
!end version
!  call effective_potential_getEnergy(effective_potential,energy_harmonic,dtset%natom,rprimd,xred,&
!&                                   comm,displacement=disp)

   end do
   
   close(funit)

 else 

!*************************************************************
!1 Convert some parameters into the structures used by mover.F90
!************************************************************
!Set mpi_eng
   mpi_enreg%comm_cell = comm

!Set the fake abinit dataset 
!Scalar
   dtset%bmass = 10     ! Barostat mass
   dtset%nctime = 0     ! NetCdf TIME between output of molecular dynamics informations 
   dtset%delayperm = 0  ! DELAY between trials to PERMUTE atoms
   dtset%dilatmx = 1    ! DILATation : MaXimal value
   dtset%dtion = 100     ! Delta Time for IONs
   dtset%diismemory = 8 ! Direct Inversion in the Iterative Subspace MEMORY
   dtset%friction = 0.0001 ! internal FRICTION coefficient
   dtset%goprecon = 0   ! Geometry Optimization PREconditioner equations
   if(inp%monte_carlo==1)then
     dtset%ionmov = 12    ! Number for the montecarlo
   end if
   dtset%jellslab = 0   ! include a JELLium SLAB in the cell
   dtset%mdwall = 10000 ! Molecular Dynamics WALL location
   dtset%natom = effective_potential%supercell%natom_supercell
   dtset%nconeq = 0     ! Number of CONstraint EQuations
   dtset%noseinert = 1.d-5 ! NOSE INERTia factor
   dtset%nnos = 5       ! Number of nose masses Characteristic
   dtset%ntime = inp%ntime  ! Number of TIME steps 
   dtset%nsym = 1       ! Number of SYMmetry operations
   dtset%optcell = 0    ! OPTimize the CELL shape and dimensions Characteristic
   dtset%restartxf = 0  ! RESTART from (X,F) history
   dtset%signperm = 1   ! SIGN of PERMutation potential      
   dtset%strprecon = 1  ! STRess PRECONditioner
   dtset%tolmxf = 2.0d-5
   dtset%tsmear = 0.009500446 !
   dtset%vis = 100      ! VIScosity
 
   dtset%usewvl = 0     !
 
!array
   dtset%goprecprm(:) = zero !Geometry Optimization PREconditioner PaRaMeters equations
   dtset%mdtemp(1) = 500   !Molecular Dynamics Temperatures 
   dtset%mdtemp(2) = 500  !Molecular Dynamics Temperatures 
   dtset%ntypat = effective_potential%ntypat
   ABI_ALLOCATE(dtset%prtatlist,(dtset%natom)) !PRinT by ATom LIST of ATom
   dtset%prtatlist(:) = zero
   ABI_ALLOCATE(dtset%iatfix,(3,dtset%natom)) ! Indices of AToms that are FIXed
   dtset%iatfix = zero
   if(dtset%nnos>0) then
     ABI_ALLOCATE(dtset%qmass,(dtset%nnos)) ! Q thermostat mass
     dtset%qmass = dtset%nnos * 10 
   end if
   dtset%strtarget = zero ! STRess TARGET
   !      dtset%symrel ! SYMmetry in REaL space
   dtset%typat  = effective_potential%supercell%typat_supercell
   dtset%useylm = 0
   dtset%znucl  = effective_potential%znucl
   
!  set psps 
   psps%useylm = 0
   
!  initialisation of results_gs
   call init_results_gs(dtset%natom,1,results_gs)

!  Set the pointers of scfcv_args
   zero_integer = 0
   scfcv_args%dtset     => dtset
   scfcv_args%mpi_enreg => mpi_enreg
   scfcv_args%ndtpawuj  => zero_integer
   scfcv_args%results_gs => results_gs
   scfcv_args%psps => psps
!  Set other arguments of the mover.F90 routines
   ABI_ALLOCATE(amass,(dtset%natom))
!  Assign masses to each atom (for MD)

   do jj = 1,dtset%natom
     amass(jj)=amu_emass*&
&      effective_potential%amu(effective_potential%supercell%typat_supercell(jj))
   end do
   acell  = one
   rprimd = effective_potential%supercell%rprimd_supercell
 
   ABI_ALLOCATE(dtset%rprimd_orig,(3,3,1))
   dtset%rprimd_orig(:,:,1) = rprimd
 
   ABI_ALLOCATE(xred,(3,dtset%natom))
   ABI_ALLOCATE(xcart,(3,dtset%natom))
   ABI_ALLOCATE(xred_old,(3,dtset%natom))
   ABI_ALLOCATE(vel,(3,dtset%natom))
   ABI_ALLOCATE(fred,(3,dtset%natom))
   ABI_ALLOCATE(fcart,(3,dtset%natom))
   
   call matr3inv(effective_potential%supercell%rprimd_supercell,gprimd)
   call xcart2xred(dtset%natom,rprimd,effective_potential%supercell%xcart_supercell,xred)
   
   xred_old = xred
   vel = zero
   vel_cell = zero
   
!  Set the dffil structure
   dtfil%filnam_ds(1:2)=filnam(1:2)
   dtfil%filnam_ds(3)=""
   dtfil%filnam_ds(4)="test"
   
   nullify (electronpositron)
   ABI_ALLOCATE(rhog,(2,1))
   ABI_ALLOCATE(rhor,(2,1))
   
!  Initialize xf history (should be put in inwffil)
   ab_xfh%nxfh=0
   ab_xfh%mxfh=(ab_xfh%nxfh-dtset%restartxf+1)+dtset%ntime+5 
   ABI_ALLOCATE(ab_xfh%xfhist,(3,dtset%natom+4,2,ab_xfh%mxfh))
   ab_xfh%xfhist(:,:,:,:) = zero
   
!*********************************************************
!2   Call main routine for monte carlo / molecular dynamics
!*********************************************************
   write(message, '(2a,(80a),a)' ) &
&    'Monte Carlo / Molecular Dynamics ',ch10,&
&    ('-',ii=1,80),ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   call mover(scfcv_args,ab_xfh,acell,amass,dtfil,electronpositron,&
&    rhog,rhor,rprimd,vel,vel_cell,xred,xred_old,effective_potential)
 

!***********************************
! 3   Deallocation of array   
!***********************************
   ABI_DEALLOCATE(amass)
   ABI_DEALLOCATE(dtset%prtatlist)
   ABI_DEALLOCATE(dtset%iatfix)
   ABI_DEALLOCATE(dtset%rprimd_orig)
   if(dtset%nnos>0) then
     ABI_DEALLOCATE(dtset%qmass)
   end if
   ABI_DEALLOCATE(rhog)
   ABI_DEALLOCATE(rhor)
   ABI_DEALLOCATE(vel)
   ABI_DEALLOCATE(xred)
   ABI_DEALLOCATE(xcart)
   ABI_DEALLOCATE(xred_old)
   ABI_DEALLOCATE(ab_xfh%xfhist)

 end if

 write(message, '(a,(80a),a,a)' ) ch10,&
&     ('=',ii=1,80),ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

end subroutine mover_effpot
!!***
