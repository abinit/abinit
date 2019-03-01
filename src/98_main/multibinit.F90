!{\src2tex{textfont=tt}}
!!****p* ABINIT/multibinit
!! NAME
!! multibinit
!!
!! FUNCTION
!! Main routine MULTIBINIT.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (AM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,abihist_bcast,abihist_free
!!      abimem_init,abinit_doctor,compute_anharmonics
!!      effective_potential_file_getdimsystem,effective_potential_file_gettype
!!      effective_potential_file_maphisttoref,effective_potential_file_read
!!      effective_potential_file_readmdfile,effective_potential_free
!!      effective_potential_setconfinement,effective_potential_writenetcdf
!!      effective_potential_writexml,fit_polynomial_coeff_fit
!!      fit_polynomial_printsystemfiles,flush_unit,herald,init10,instrng
!!      inupper,invars10,isfile,mover_effpot,multibinit_dtset_free
!!      outvars_multibinit,timein,wrtout,xmpi_bcast,xmpi_init,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program multibinit

 use defs_basis
 use defs_abitypes
 use m_build_info
 use m_xmpi
 use m_xomp
 use m_abicore
 use m_errors
 use m_argparse
 use m_effective_potential
 use m_fit_polynomial_coeff
 use m_opt_effpot
 use m_multibinit_dataset
 use m_effective_potential_file
 use m_spin_model
 use m_abihist

 use m_specialmsg, only : specialmsg_getcount, herald
 use m_io_tools,   only : flush_unit, open_file
 use m_fstrings,   only : replace, inupper
 use m_time,       only : asctime, timein
 use m_parser,     only : instrng
 use m_dtset,      only : chkvars
 use m_dtfil,      only : isfile
 use m_mover_effpot, only : mover_effpot
#if defined DEV_MS_SCALEUP 
 use scup_global, only : global_init_model
#endif 
 !use m_generate_training_set, only : generate_training_set
 use m_compute_anharmonics, only : compute_anharmonics
 use m_init10,              only : init10
 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
! Set array dimensions
 integer,parameter :: master=0 ! FIXME: these should not be reserved unit numbers!
 integer :: comm,filetype,ii,ierr,lenstr
 integer :: natom,nph1l,nrpt,ntypat,nproc,my_rank
 integer :: option
 logical :: iam_master
 logical :: need_analyze_anh_pot
 logical :: err_init_elec
 real(dp) :: tcpu,tcpui,twall,twalli
 real(dp) :: tsec(2)
 character(len=24) :: codename,start_datetime
 character(len=strlen) :: string
 character(len=fnlen) :: filnam(18),tmpfilename,name
 character(len=500) :: message
 type(multibinit_dtset_type) :: inp
 type(effective_potential_type) :: reference_effective_potential
 type(abihist) :: hist, hist_tes ! TODO MARCUS: test set implementation hist_tes: hist_TestSet
 type(args_t) :: args

!TODO hexu: add types for spin here.
 type(spin_model_t) :: spin_model
!TEST_AM
! integer :: natom_sp
! real(dp),allocatable :: dynmat(:,:,:,:,:)
!TEST_AM
!******************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

!Initialize MPI
 call xmpi_init()

!MPI variables
 comm = xmpi_world
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)


! Parse command line arguments.
 args = args_parser(); if (args%exit /= 0) goto 100

!Initialize memory profiling if it is activated !if a full abimem.mocc report is desired,
!set the argument of abimem_init to "2" instead of "0"
!note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(2)
#endif

!Initialisation of the timing
 call timein(tcpui,twalli)

 if (iam_master) then
   codename='MULTIBINIT'//repeat(' ',14)
   call herald(codename,abinit_version,std_out)
 end if

 start_datetime = asctime()

!Print the number of cpus in log file
 write(message,'(a,i5,a)') '-  nproc =',nproc,ch10
 call wrtout(std_out,message,'COLL')

!Initialise the code : write heading, and read names of files.
 call init10(filnam,comm)

!******************************************************************

 call timein(tcpu,twall)

 write(message, '(a,f11.3,a,f11.3,a)' )'-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')

! Open output files and ab_out (might change its name if needed)
! MJV 1/2010 : now output file is open, but filnam(2) continues unmodified
! so the other output files are overwritten instead of accumulating.
 if (iam_master) then
   tmpfilename = filnam(2)
   call isfile(tmpfilename,'new')
   if (open_file(tmpfilename,message,unit=ab_out,form="formatted",status="new",&
&   action="write") /= 0) then
     MSG_ERROR(message)
   end if
!  Call open_file(unit=ab_out,file=tmpfilename,form='formatted',status='new')
   rewind (unit=ab_out)
   call herald(codename,abinit_version,ab_out)
!  Print the number of cpus in output
   write(message,'(a,i5,a)') '-  nproc =',nproc
   call wrtout(ab_out,message,'COLL')
 else
   ab_out = dev_null
 end if

 write(message, '(a,(80a),a)' ) ch10,&
& ('=',ii=1,80),ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')


!To automate a maximum calculation, multibinit reads the number of atoms
!in the file (ddb or xml). If DDB file is present in input, the ifc calculation
!will be initilaze array to the maximum of atoms (natifc=natom,atifc=1,natom...) in invars10
 write(message, '(6a)' )' Read the information in the reference structure in ',ch10,&
& '-',trim(filnam(3)),ch10,' to initialize the multibinit input'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 call effective_potential_file_getDimSystem(filnam(3),natom,ntypat,nph1l,nrpt)

!Read the input file, and store the information in a long string of characters
!strlen from defs_basis module
 option=1
 if (iam_master) then
   call instrng (filnam(1),lenstr,option,strlen,string)
   !To make case-insensitive, map characters to upper case:
   call inupper(string(1:lenstr))

   !Check whether the string only contains valid keywords
   call chkvars(string)

 end if

 call xmpi_bcast(string,master, comm, ierr)
 call xmpi_bcast(lenstr,master, comm, ierr)

!Read the input file
 call invars10(inp,lenstr,natom,string)

 if (iam_master) then
!  Echo the inputs to console and main output file
   call outvars_multibinit(inp,std_out)
   call outvars_multibinit(inp,ab_out)
 end if

! Read and treat the reference structure
!****************************************************************************************
  if (inp%spin_dynamics>0) then
   if (iam_master) then

     write(message,'(a,(80a),3a)') ch10,('=',ii=1,80),ch10,ch10,&
&     'reading spin terms.'
     call spin_model_t_initialize(spin_model, filnam, inp )
   end if
 else
   !  Read the model (from DDB or XML)
   call effective_potential_file_read(filnam(3),reference_effective_potential,inp,comm)

 !Read the coefficient from fit
   if(filnam(4)/=''.and.filnam(4)/='no')then
     call effective_potential_file_getType(filnam(4),filetype)
   ! TODO hexu: filetype==(33?)
     if(filetype==3.or.filetype==23) then
       call effective_potential_file_read(filnam(4),reference_effective_potential,inp,comm)
     else
       write(message,'(a,(80a),3a)') ch10,('=',ii=1,80),ch10,ch10,&
&         ' There is no specific file for the coefficients from polynomial fitting'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if
   else
     if(inp%ncoeff/=0) then
       write(message, '(5a)' )&
&         'ncoeff is specified in the input but,',ch10,&
&         'there is no file for the coefficients ',ch10,&
&         'Action: add coefficients.xml file'
       MSG_ERROR(message)

     else
       write(message,'(a,(80a),3a)') ch10,('=',ii=1,80),ch10,ch10,&
&         ' There is no file for the coefficients from polynomial fitting'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if
   end if
 end if

 !Initialized the electronic model (If scale-up is available)
#if defined DEV_MS_SCALEUP 
  err_init_elec = global_init_model(filnam(4),inp%ncell,needlattice=.FALSE.,needelectrons=.TRUE.,didi=.FALSE.,harm_der=.FALSE.,tcharge=0,ismagnetic=.FALSE.,istddft=.FALSE.) 
#endif 

!****************************************************************************************

! Compute the third order derivative with finite differences
!****************************************************************************************
   if (inp%strcpling > 0) then
     call compute_anharmonics(reference_effective_potential,filnam,inp,comm)
   end if
!****************************************************************************************

! If needed, fit the anharmonic part and compute the confinement potential
!****************************************************************************************
   if (inp%fit_coeff/=0.or.inp%confinement==2.or.inp%bound_model/=0 .or. inp%opt_effpot/=0) then

     if(iam_master) then
!    Read the MD file
       write(message,'(a,(80a),7a)')ch10,('=',ii=1,80),ch10,ch10,&
&       '-Reading the training-set file :',ch10,&
&       '-',trim(filnam(5)),ch10

       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,message,'COLL')
       if(filnam(5)/=''.and.filnam(5)/='no')then
         call effective_potential_file_readMDfile(filnam(5),hist,option=inp%ts_option)
         if (hist%mxhist == 0)then
           write(message, '(5a)' )&
&           'The trainig-set ',trim(filnam(5)),' file is not correct ',ch10,&
&           'Action: add training-set file'
           MSG_ERROR(message)
         end if
       else
         if (inp%fit_coeff/=0) then
           write(message, '(3a)' )&
&           'There is no training-set file to fit the lattice model ',ch10,&
&           'Action: add trainings-set-file'
           MSG_ERROR(message)
         else if (inp%bound_model/=0) then
             write(message, '(3a)' )&
&             'There is no  training-set file to bound the model ',ch10,&
&             'Action: add training-set file '
             MSG_ERROR(message)
         else if(inp%confinement==2) then
             write(message, '(3a)' )&
&             'There is no training-set file to compute the confinement',ch10,&
&             'Action: add training-set file '
             MSG_ERROR(message) 
         else if(inp%opt_effpot==2) then
             write(message, '(3a)' )&
&             'There is no training-set file to optimize the latice model',ch10,&
&             'Action: add training-set file '
             MSG_ERROR(message)
         end if
       end if
     end if
!  MPI BROADCAST the history of the MD
     call abihist_bcast(hist,master,comm)
!  Map the hist in order to be consistent with the supercell into reference_effective_potential
     call effective_potential_file_mapHistToRef(reference_effective_potential,hist,comm)
   end if

!TEST_AM
! call effective_potential_checkDEV(reference_effective_potential,hist,size(hist%xred,2),hist%mxhist)
! stop
!TEST_AM

!Generate the confinement polynome (not working yet)
   if(inp%confinement/=0)then
     option=inp%confinement
     select case(option)
     case(1)
       call effective_potential_setConfinement(inp%conf_cutoff_disp,inp%conf_cutoff_strain,&
&       reference_effective_potential,inp%conf_power_fact_disp,&
&       inp%conf_power_fact_strain,inp%conf_power_disp,&
&       inp%conf_power_disp,inp%conf_power_strain,&
&       need_confinement=.TRUE.)

       write(message,'(a,(80a),3a)') ch10,('=',ii=1,80),ch10,ch10,&
&       ' The confinement potential is active.'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')

     case(2)
       write(message,'(a,(80a),3a)') ch10,('=',ii=1,80),ch10,ch10,&
&       ' The confinement potential is computed from the MD file and actived.'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')

     end select
   end if


!Fit the coeff
   if (inp%fit_coeff/=0)then
     option=inp%fit_coeff
     if(hist%mxhist >0)then
       if (option==-1)then
!      option == -1
!      Print the file in the specific format for the script of carlos
!      Born_Charges
!      Dielectric_Tensor
!      harmonic.xml
!      Reference_structure
!      Strain_Tensor
!      symmetry_operations (only cubic)
         if (iam_master) then
           call fit_polynomial_printSystemFiles(reference_effective_potential,hist)
         end if
       else if (option==1.or.option==2)then
!      option = 1
         call fit_polynomial_coeff_fit(reference_effective_potential,&
&         inp%fit_bancoeff,inp%fit_fixcoeff,hist,inp%fit_generateCoeff,&
&         inp%fit_rangePower,inp%fit_nbancoeff,inp%fit_ncoeff,&
&         inp%fit_nfixcoeff,option,comm,cutoff_in=inp%fit_cutoff,&
&         initialize_data=inp%fit_initializeData==1,&
&         fit_tolMSDF=inp%fit_tolMSDF,fit_tolMSDS=inp%fit_tolMSDS,fit_tolMSDE=inp%fit_tolMSDE,&
&         fit_tolMSDFS=inp%fit_tolMSDFS,&
&         verbose=.true.,positive=.false.,&
&         anharmstr=inp%fit_anhaStrain==1,&
&         spcoupling=inp%fit_SPCoupling==1,prt_names=inp%prt_names,prt_anh=inp%analyze_anh_pot)
       end if
     else
       write(message, '(3a)' )&
&       'There is no step in the MD file ',ch10,&
&       'Action: add correct MD file'
       MSG_ERROR(message)
     end if
   end if


!try to bound the model with mover_effpot
!we need to use the molecular dynamics
   if(inp%bound_model>0.and.inp%bound_model<=2)then
     call mover_effpot(inp,filnam,reference_effective_potential,-1*inp%bound_model,comm,hist=hist)
   !Marcus: New option for bound_model: use optimize routine for generting specific high order terms
   elseif(inp%bound_model == 3)then
    write(message,'(a,(80a),4a)')ch10,('=',ii=1,80),ch10,ch10,&
&    'Bound Process 3: Generate equivalent high order terms',ch10            
     call opt_effpotbound(reference_effective_potential,inp%bound_rangePower,hist,comm) 
   end if





!****************************************************************************************
! OPTIMIZE SECTION, Optimize selected coefficients of effective potential while
! keeping the others constant
!****************************************************************************************
 
 if(inp%opt_effpot == 1)then
    write(message,'(a,(80a),4a)')ch10,('=',ii=1,80),ch10,ch10,&
&    'Optimizing Effective Potential',ch10            

    call wrtout(std_out,message,'COLL')
    call wrtout(ab_out,message,'COLL')

     need_analyze_anh_pot = .FALSE.
     if(inp%analyze_anh_pot == 1) need_analyze_anh_pot = .TRUE. 

    call opt_effpot(reference_effective_potential,inp%opt_ncoeff,inp%opt_coeff,hist,comm,& 
&        print_anh=need_analyze_anh_pot)
 end if 

 

!****************************************************************************************
! TEST SECTION test effective potential with regard to test-set 
!****************************************************************************************
   if(inp%test_effpot == 1)then 
     if(iam_master) then
!    Read the test-set .nc file
       write(message,'(a,(80a),9a)')ch10,('=',ii=1,80),ch10,ch10,&
&       'TEST - SET Option',ch10,&               
&       '-Reading the test-set file :',ch10,&
&       '-',trim(filnam(6)),ch10

       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,message,'COLL')
       if(filnam(6)/=''.and.filnam(6)/='no')then
         call effective_potential_file_readMDfile(filnam(6),hist_tes,option=inp%ts_option)
         if (hist_tes%mxhist == 0)then
           write(message, '(5a)' )&
&           'The test-set ',trim(filnam(6)),' file is empty ',ch10,&
&           'Action: add non-empty test-set'
           MSG_ERROR(message)
         end if
       else
           write(message, '(3a)' )&
&           'There is no test-set file ',ch10,&
&           'Action: add test-set file'
           MSG_ERROR(message)
       end if
     end if
!  MPI BROADCAST the history of the MD
     call abihist_bcast(hist_tes,master,comm)
!  Map the hist in order to be consistent with the supercell into reference_effective_potential
     call effective_potential_file_mapHistToRef(reference_effective_potential,hist_tes,comm)
     !  Initialize if to print anharmonic contribution to energy or not   
     need_analyze_anh_pot = .FALSE.
     if(inp%analyze_anh_pot == 1) need_analyze_anh_pot = .TRUE.  
!  Call to test routine 
     call fit_polynomial_coeff_testEffPot(reference_effective_potential,hist_tes,master,comm,print_anharmonic=need_analyze_anh_pot)

   end if ! End if(inp%test_effpot == 1)then 

!****************************************************************************************

!****************************************************************************************
!Print the effective potential system + coefficients (only master CPU)
! TODO hexu: add print spin model.
   if(iam_master) then
     if (inp%prt_model >= 1) then
       write(message, '(a,(80a),a)' ) ch10,&
&       ('=',ii=1,80)
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
       name = replace(trim(filnam(2)),".out","")
       call effective_potential_writeXML(reference_effective_potential,inp%prt_model,filename=name,&
&       prt_dipdip=inp%dipdip_prt==1)
     else if (inp%prt_model == -2)then
!    NetCDF case, in progress
       name = trim(filnam(2))//"_sys.nc"
       call effective_potential_writeNETCDF(reference_effective_potential,1,filename=name)
     end if
   end if
!****************************************************************************************

!TEST_AM SECTION
! Print the Phonon dos/spectrum
! if(inp%prt_phfrq > 0) then
!     call effective_potential_printPDOS(reference_effective_potential,filnam(2),&
!&           inp%ncell,inp%nph1l,inp%prt_phfrq,inp%qph1l)
!   end if

!Intialisation of the effective potential type
!   call effective_potential_file_read(filnam(4),reference_effective_potential,inp,comm)
!   name = "test.xml"
!   call effective_potential_writeXML(reference_effective_potential,1,filename=name)
! just for TEST
!   if(inp%prt_phfrq > 0) then
!      natom_sp = reference_effective_potential%supercell%natom_supercell
!      ABI_ALLOCATE(dynmat,(2,3,natom_sp,3,natom_sp))
!      call effective_potential_effpot2dynmat(dynmat,inp%delta_df,reference_effective_potential,&
! &                                           reference_effective_potential%supercell%natom_supercell,&
! &                                           int(reference_effective_potential%supercell%qphon),3)

!      ABI_DEALLOCATE(dynmat)
!    end if
! end if
!TEST_AM SECTION


! Compute the monte carlo, molecular dynamics of compute specific energy
!****************************************************************************************
   if(inp%dynamics>=1) then
     call mover_effpot(inp,filnam,reference_effective_potential,inp%dynamics,comm)
   end if

!****************************************************************************************


! Run spin dynamics
!****************************************************************************************
   if(inp%spin_dynamics>0) then
  ! TODO hexu: no mpi yet.
     if(iam_master) then
       call spin_model_t_run(spin_model)
     end if
   end if
!****************************************************************************************


!Free the effective_potential and dataset
!****************************************************************************************
   call effective_potential_free(reference_effective_potential)
   call multibinit_dtset_free(inp)
   call abihist_free(hist)
   call abihist_free(hist_tes)
   call spin_model_t_finalize(spin_model)
!****************************************************************************************

   write(message,'(a,a,a,(80a))') ch10,('=',ii=1,80),ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   call timein(tcpu,twall)
   tsec(1)=tcpu-tcpui
   tsec(2)=twall-twalli

   write(message, '(a,i4,a,f13.1,a,f13.1)' )' Proc.',my_rank,' individual time (sec): cpu=',&
&   tsec(1),'  wall=',tsec(2)
   call wrtout(std_out,message,"COLL")

   if (iam_master) then
     write(ab_out, '(a,a,a,i4,a,f13.1,a,f13.1)' )'-',ch10,&
&     '- Proc.',my_rank,' individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)
   end if

   call xmpi_sum(tsec,comm,ierr)

   write(message, '(a,(80a),a,a,a,f11.3,a,f11.3,a,a,a,a)' ) ch10,&
&   ('=',ii=1,80),ch10,ch10,&
&   '+Total cpu time',tsec(1),&
&   '  and wall time',tsec(2),' sec',ch10,ch10,&
&   ' multibinit : the run completed succesfully.'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   if (iam_master) then
   ! Write YAML document with the final summary.
   ! we use this doc to test whether the calculation is completed.
     write(std_out,"(a)")"--- !FinalSummary"
     write(std_out,"(a)")"program: multibinit"
     write(std_out,"(2a)")"version: ",trim(abinit_version)
     write(std_out,"(2a)")"start_datetime: ",start_datetime
     write(std_out,"(2a)")"end_datetime: ",asctime()
     write(std_out,"(a,f13.1)")"overall_cpu_time: ",tsec(1)
     write(std_out,"(a,f13.1)")"overall_wall_time: ",tsec(2)
     write(std_out,"(a,i0)")"mpi_procs: ",xmpi_comm_size(xmpi_world)
     write(std_out,"(a,i0)")"omp_threads: ",xomp_get_num_threads(open_parallel=.True.)
   !write(std_out,"(a,i0)")"num_warnings: ",nwarning
   !write(std_out,"(a,i0)")"num_comments: ",ncomment
     write(std_out,"(a)")"..."
     call flush_unit(std_out)
   end if

!Write information on file about the memory before ending mpi module, if memory profiling is enabled
   call abinit_doctor("__multibinit")

   call flush_unit(ab_out)
   call flush_unit(std_out)

   if (iam_master) close(ab_out)

   100 call xmpi_end()

   end program multibinit
!!***
