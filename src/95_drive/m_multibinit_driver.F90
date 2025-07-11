
!!****m* ABINIT/m_multibinit_main
!! NAME
!! m_multibinit_main
!!
!! FUNCTION
!! Main routine MULTIBINIT.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2025 ABINIT group (AM, hexu)
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
!! NOTES
!! Should be
!! 1 moved to somewhere else
!! 2 be replaced with the new implementation multibinit_main2.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

! FIXME: This module is temporarily here. It should be removed once we have the new lattice mover.
! Move the main calculation part into a subroutine
! TODO: And this subroutine should be replaced by a one liner then: call multibinit_manager%run_all().
! TODO: The module is here because it uses lv95_drive mover_effpot (very hacky).
! TODO: The mover_effpot is in 95_drive due to that it use some >78 level module
! TODO: which should be changed when we implement the new lattice mover.
module m_multibinit_driver
  use defs_basis
  use defs_abitypes
  use m_xmpi
  use m_xomp
  use m_abicore
  use m_errors

  use m_effective_potential
  use m_fit_polynomial_coeff
 use m_opt_effpot
  use m_multibinit_dataset
  use m_effective_potential_file
 use m_scup_dataset
  !use m_spin_model, only: spin_model_t
  use m_abihist

  use m_build_info,         only: abinit_version
  use m_multibinit_manager, only: mb_manager_t
  use m_multibinit_main2,   only: multibinit_main2

  use m_mover_effpot, only : mover_effpot
#if defined DEV_MS_SCALEUP
 use scup_global
#endif
  !use m_generate_training_set, only : generate_training_set
  use m_compute_anharmonics, only : compute_anharmonics
  use m_init10,              only : init10, postfix_fnames
  use m_parser,     only : instrng
  use m_fstrings,   only : replace, inupper
  use m_dtset,      only : chkvars
  implicit none
  !!***
contains
  !!****f* m_multbinit_main/multibinit_main
  !!
  !! NAME
  !! multibinit_main
  !!
  !! FUNCTION
  !! The main function of multibinit
  !!
  !! INPUTS
  !! filnam: The filenames from the files file. 17 files in total.
  !!
  !! OUTPUT
  !!
  !! SOURCE
  subroutine multibinit_main(input_path, filnam, dry_run)
    character(len=fnlen), intent(inout) :: input_path
    character(len=fnlen), intent(inout) :: filnam(18)
    integer, intent(in) :: dry_run
    type(multibinit_dtset_type), target :: inp
    type(effective_potential_type) :: reference_effective_potential, read_effective_potential
    type(abihist) :: hist, hist_tes

    !type(spin_model_t) :: spin_model
    character(len=strlen) :: string, raw_string
    character(len=500) :: message
    character(len=fnlen) :: name
    character(len=fnlen) :: sys_fname

    integer :: filetype,ii,lenstr,iiter,niter
    integer :: natom,nph1l,nrpt,ntypat
    integer :: option
    logical :: need_analyze_anh_pot,need_prt_files

    ! Whether the "new" MULTIBNIT framework should be used.
    logical :: need_new_multibinit
! MS
! temporary variables for testing SCALE-UP with Multibinit
  !Variable to pass to effpot_evaluate routine of multibinit
  !To declare evaluation of electronice model
  logical  :: elec_eval
#if defined DEV_MS_SCALEUP
  !Variables needed to call SCALE-UP
  logical :: err_init_elec
  logical*1 :: needlattice = .FALSE.
  logical*1 :: needelectrons = .TRUE.
  logical*1 :: didi = .FALSE.
  logical*1 :: harm_der = .FALSE.
  logical*1 :: initorbocc = .FALSE.
  logical*1 :: ismagnetic = .FALSE.
  logical*1 :: istddft = .FALSE.
  logical*4 :: printgeom = .FALSE.
  logical*4 :: printeigv = .FALSE.
  logical*4 :: printeltic = .FALSE.
  logical*4 :: printorbocc = .FALSE.
  integer :: ksamp(3)
  real*8 :: tcharge
#endif
!TEST_AM
! integer :: natom_sp
! real(dp),allocatable :: dynmat(:,:,:,:,:)
!TEST_AM
!******************************************************************

    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master


!MPI variables
    master = 0
    comm = xmpi_world
    nproc = xmpi_comm_size(comm)
    my_rank = xmpi_comm_rank(comm)
    iam_master = (my_rank == master)

    !Read the input file, only to the the name of the file which contains the ddb file or xml file
    ! for getting the number of atoms.
    option=1
    if (iam_master) then
       call instrng (filnam(1),lenstr,option,strlen,string,raw_string)
       !To make case-insensitive, map characters to upper case:
       call inupper(string(1:lenstr))
       !Check whether the string only contains valid keywords
       call chkvars(string)
    end if
    call xmpi_bcast(string,master, comm, ierr)
    call xmpi_bcast(raw_string,master, comm, ierr)
    call xmpi_bcast(lenstr,master, comm, ierr)
    !To automate a maximum calculation, multibinit reads the number of atoms
    !in the file (ddb or xml). If DDB file is present in input, the ifc calculation
    !will be initilaze array to the maximum of atoms (natifc=natom,atifc=1,natom...) in invars10


    !Read the input file assuming natom=1 so that the invars10 can work.

    !call invars10(inp,lenstr,natom,string)
    if(trim(filnam(3)) /='') then
      sys_fname=filnam(3)
    else
      call invars_multibinit_filenames(string, lenstr,  sys_fname=sys_fname)
    end if



    ! read the reference structure to get natom
    if (iam_master) then
      write(message, '(6a)' )' Read the information in the reference structure in ',ch10,&
            & '-',trim(sys_fname),ch10,' to initialize the multibinit input'
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')
    end if

    call effective_potential_file_getDimSystem(sys_fname,comm,natom,ntypat,nph1l,nrpt)
    !call effective_potential_file_getDimSystem(filnam(3),natom,ntypat,nph1l,nrpt)


    ! read the input again to use the right natom
    call invars10(inp,lenstr,natom,string)
    call postfix_fnames(input_path, filnam, inp)

    need_new_multibinit= inp%spin_dynamics > 0 .or. inp%lwf_dynamics > 0 .or. inp%dynamics >= 100

    if (iam_master) then
        if(need_new_multibinit) then
            ABI_ERROR("The new MULTINIT mode should be enabled with --F03 option. ")
        end if
       !  Echo the inputs to console and main output file
       call outvars_multibinit(inp,std_out)
       call outvars_multibinit(inp,ab_out)
    end if

    if(dry_run/=0) then
       call wrtout([std_out, ab_out], "Multibinit in dry_run mode. Exiting after input parser")
       call xmpi_end()
       !goto 100
    endif

    !  Read the model (from DDB or XML)
    call effective_potential_file_read(filnam(3),reference_effective_potential,inp,comm)
      ! TODO: coeff_file_rw is currently set to 0 as this functionality needs more testing. 
       if (inp%coeff_file_rw==2) then
          read_effective_potential=reference_effective_potential
          !call effective_potential_file_read(filnam(3),read_effective_potential,inp,comm)
          if(filnam(6)/=''.and.filnam(6)/='no')then
              call effective_potential_file_getType(filnam(6),filetype)
              if(filetype==3.or.filetype==23) then
                  call effective_potential_file_read(filnam(6),read_effective_potential,inp,comm)
              end if
          end if
       end if

    if(filnam(4)/=''.and.filnam(4)/='no') then
       call effective_potential_file_getType(filnam(4),filetype)
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
          ABI_ERROR(message)

       else
          write(message,'(a,(80a),3a)') ch10,('=',ii=1,80),ch10,ch10,&
               &         ' There is no file for the coefficients from polynomial fitting'
          call wrtout(ab_out,message,'COLL')
          call wrtout(std_out,message,'COLL')
       end if
    end if

!****************************************************************************************
!SCALE UP Initialize the electronic model (If scale-up is available)
!****************************************************************************************
elec_eval = .FALSE.

#if defined DEV_MS_SCALEUP
 if(inp%scup_dtset%scup_elec_model)then
   write(message,'(a,(80a),4a)') ch10,('=',ii=1,80),ch10,ch10,&
        ' Initializing Electronic Model with SCALE-UP',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   !Set Variables
   elec_eval = .TRUE.
   ksamp = inp%scup_dtset%scup_ksamp
   tcharge = inp%scup_dtset%scup_tcharge
   if(inp%scup_dtset%scup_ismagnetic)ismagnetic=.TRUE.
   if(inp%scup_dtset%scup_istddft)istddft=.TRUE.
   if(inp%scup_dtset%scup_initorbocc)initorbocc=.TRUE.

   ! Call to Scale-Up
   err_init_elec = global_init_model(filnam(3),inp%ncell,needlattice,needelectrons,didi,&
&                               harm_der,tcharge,ksamp,ismagnetic,istddft,initorbocc)

   !Set Print variables
   if(inp%scup_dtset%scup_printgeom)printgeom=.TRUE.
   if(inp%scup_dtset%scup_printeigv)printeigv=.TRUE.
   if(inp%scup_dtset%scup_printeltic)printeltic=.TRUE.
   if(inp%scup_dtset%scup_printorbocc)printorbocc=.TRUE.

   !Set Print Parameters within scaleup
   call global_set_print_parameters(printgeom,printeigv,printeltic,printorbocc,&
&                                  inp%scup_dtset%scup_printbands)

   !Set SCF controling variables (values contain defaults, if not specified in the input)
   call global_set_scf_parameters(inp%scup_dtset%scup_scfmixing,inp%scup_dtset%scup_scfthresh,&
&                                 inp%scup_dtset%scup_smearing,inp%scup_dtset%scup_maxscfstep,&
&                                 inp%scup_dtset%scup_startpulay,inp%scup_dtset%scup_freezden)


   !Create kpath if printbands=true and pass it to SCALE UP
   if(inp%scup_dtset%scup_printbands)then

           call scup_kpath_new(inp%scup_dtset%scup_speck,&
&                                      reference_effective_potential%supercell%rprimd,&
&                                inp%scup_dtset%scup_ndivsm,inp%scup_dtset%scup_kpath)
           call scup_kpath_print(inp%scup_dtset%scup_kpath)

           call global_set_print_bands(inp%scup_dtset%scup_printbands,&
&               inp%scup_dtset%scup_nspeck,inp%scup_dtset%scup_kpath%ndivs,&
&               inp%scup_dtset%scup_speck)
   endif
 endif
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
           ABI_ERROR(message)
         end if
       else
         if (inp%fit_coeff/=0) then
           write(message, '(3a)' )&
&           'There is no training-set file to fit the lattice model ',ch10,&
&           'Action: add trainings-set-file'
           ABI_ERROR(message)
         else if (inp%bound_model/=0) then
             write(message, '(3a)' )&
&             'There is no  training-set file to bound the model ',ch10,&
&             'Action: add training-set file '
             ABI_ERROR(message)
         else if(inp%confinement==2) then
             write(message, '(3a)' )&
&             'There is no training-set file to compute the confinement',ch10,&
&             'Action: add training-set file '
             ABI_ERROR(message)
         else if(inp%opt_effpot==2) then
             write(message, '(3a)' )&
&             'There is no training-set file to optimize the latice model',ch10,&
&             'Action: add training-set file '
             ABI_ERROR(message)
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
             if(inp%fit_iatom/=0)then
               call fit_polynomial_coeff_fit(reference_effective_potential,&
                    &         inp%fit_bancoeff,inp%fit_fixcoeff,hist,inp%fit_generateCoeff,&
                    &         inp%fit_rangePower,inp%fit_nbancoeff,inp%fit_ncoeff,&
                    &         inp%fit_nfixcoeff,inp%fit_nimposecoeff,inp%fit_imposecoeff,&
                    &         option,comm,cutoff_in=inp%fit_cutoff,&
                    &         max_power_strain=inp%fit_SPC_maxS,initialize_data=inp%fit_initializeData==1,&
                    &         fit_tolMSDF=inp%fit_tolMSDF,fit_tolMSDS=inp%fit_tolMSDS,fit_tolMSDE=inp%fit_tolMSDE,&
                    &         fit_tolMSDFS=inp%fit_tolMSDFS,fit_tolGF=inp%fit_tolGF,&
                    &         verbose=.true.,positive=.false.,&
                    &         anharmstr=inp%fit_anhaStrain==1,&
                    &         spcoupling=inp%fit_SPCoupling==1,prt_anh=inp%analyze_anh_pot,&
                    &         fit_iatom=inp%fit_iatom,prt_files=.TRUE.,fit_on=inp%fit_on,sel_on=inp%sel_on,&
                    &         fit_factors=inp%fit_factors,prt_GF_csv=inp%prt_GF_csv,dispterms=inp%fit_dispterms==1,&
                    &         coeff_file_rw=inp%coeff_file_rw,read_effective_potential=read_effective_potential, &
                    &         max_nbody=inp%fit_max_nbody,  &
                    &         drop_rate=inp%fit_drop_rate, ncoeff_per_cycle=inp%fit_ncoeff_per_cycle, fit_weight_T=inp%fit_weight_T)
             else
                if (inp%fit_ncoeff_per_iatom/=0)then
                   if (mod(inp%fit_ncoeff,inp%fit_ncoeff_per_iatom) /= 0)then
                      write(message,'(2a,I3,2a,I3,3a)') ch10,&
                           & 'fit_ncoeff_per_iatom = ', inp%fit_ncoeff_per_iatom,ch10,&
                           & 'is not a divider of fit_ncoeff = ', inp%fit_ncoeff,ch10,&
                           & 'Action: Change fit_ncoeff and/or fit_ncoeff_per_iatom',ch10
                      ABI_ERROR(message)
                   endif
                   niter = inp%fit_ncoeff/inp%fit_ncoeff_per_iatom
                   if (mod(niter,reference_effective_potential%crystal%nirredat) /= 0)then
                      write(message,'(2a,I3,2a,I3,2a,I3,3a)') ch10,&
                           & 'fit_ncoeff_per_iatom = ', inp%fit_ncoeff_per_iatom,ch10,&
                           & 'times the number of irreducible atoms = ',reference_effective_potential%crystal%nirredat,ch10,&
                           & 'is not a divider of fit_ncoeff = ', inp%fit_ncoeff,ch10,&
                           & 'Action: Change fit_ncoeff and/or fit_ncoeff_per_iatom',ch10
                      ABI_ERROR(message)
                   endif
                   niter = niter/reference_effective_potential%crystal%nirredat
                else if (inp%fit_ncoeff_per_iatom == 0)then
                   if (mod(inp%fit_ncoeff,reference_effective_potential%crystal%nirredat) /= 0)then
                      write(message,'(2a,I3,2a,I3,3a)') ch10,&
                           & 'The number of irreducible atoms = ',reference_effective_potential%crystal%nirredat,ch10,&
                           & 'is not a divider of fit_ncoeff = ', inp%fit_ncoeff,ch10,&
                           & 'Action: Change fit_ncoeff',ch10
                      ABI_ERROR(message)
                   endif
                   inp%fit_ncoeff_per_iatom = inp%fit_ncoeff/reference_effective_potential%crystal%nirredat
                   niter = 1
                endif
                write(message,'(a,(80a),7a,I3,3a,I3,3a,I3,3a,I3,2a)') ch10,('=',ii=1,80),ch10,ch10,&
                     & '  Starting Fit Iterations  ',ch10,&
                     & '  -----------------------  ',ch10,&
                     & '  Select in total fit_ncoeff = ', inp%fit_ncoeff,' coefficients',ch10,&
                     & '  In ', niter,' iterations',ch10,&
                     & '  Over ', reference_effective_potential%crystal%nirredat, ' irreducible atoms',ch10,&
                     & '  Selecting ', inp%fit_ncoeff_per_iatom, ' coefficients per atom in each iteration',ch10
                call wrtout(std_out,message,'COLL')
                call wrtout(ab_out,message,'COLL')
                need_prt_files=.FALSE.
                do iiter=1,niter
                  write(message,'(a,(80a),3a,I3,a,I3,2a)') ch10,('-',ii=1,80),ch10,ch10,&
                          &    ' Start Iteration (',iiter,'/',niter,')',ch10
                  call wrtout(std_out,message,'COLL')
                  call wrtout(ab_out,message,'COLL')
                  do ii=1,reference_effective_potential%crystal%nirredat
                    if(ii == reference_effective_potential%crystal%nirredat .and. iiter==niter)need_prt_files=.TRUE.
                    if(ii > 1 .or. iiter > 1)inp%fit_nfixcoeff = -1
                       call fit_polynomial_coeff_fit(reference_effective_potential,&
                          &         inp%fit_bancoeff,inp%fit_fixcoeff,hist,inp%fit_generateCoeff,&
                          &         inp%fit_rangePower,inp%fit_nbancoeff,inp%fit_ncoeff_per_iatom,&
                          &         inp%fit_nfixcoeff,inp%fit_nimposecoeff,inp%fit_imposecoeff,&
                          &         option,comm,cutoff_in=inp%fit_cutoff,&
                          &         max_power_strain=inp%fit_SPC_maxS,initialize_data=inp%fit_initializeData==1,&
                          &         fit_tolMSDF=inp%fit_tolMSDF,fit_tolMSDS=inp%fit_tolMSDS,fit_tolMSDE=inp%fit_tolMSDE,&
                          &         fit_tolMSDFS=inp%fit_tolMSDFS,fit_tolGF=inp%fit_tolGF,&
                          &         verbose=.true.,positive=.false.,&
                          &         anharmstr=inp%fit_anhaStrain==1,&
                          &         spcoupling=inp%fit_SPCoupling==1,prt_anh=inp%analyze_anh_pot,&
                          &         fit_iatom=reference_effective_potential%crystal%irredatindx(ii),&
                          &         prt_files=need_prt_files,fit_on=inp%fit_on,sel_on=inp%sel_on,&
                          &         fit_factors=inp%fit_factors,prt_GF_csv=inp%prt_GF_csv,dispterms=inp%fit_dispterms==1, &
   &                                coeff_file_rw=inp%coeff_file_rw,read_effective_potential=read_effective_potential, &
   &                      max_nbody=inp%fit_max_nbody,  &
   &                      drop_rate=inp%fit_drop_rate, ncoeff_per_cycle=inp%fit_ncoeff_per_cycle, fit_weight_T=inp%fit_weight_T)
                  enddo
                enddo
             endif
          end if
       else
          write(message, '(3a)' )&
               &       'There is no step in the MD file ',ch10,&
               &       'Action: add correct MD file'
          ABI_ERROR(message)
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
    call wrtout(std_out,message,'COLL')
    call wrtout(ab_out,message,'COLL')

    call opt_effpotbound(reference_effective_potential,inp%bound_rangePower,hist, inp%bound_EFS,&
&                       inp%bound_factors,inp%bound_penalty,comm, fit_weight_T=inp%fit_weight_T, bound_option=inp%bound_option)

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

    call opt_effpot(reference_effective_potential,inp%opt_ncoeff,inp%opt_coeff,hist,inp%opt_on,&
&                   inp%opt_factors,comm,fit_weight_T=inp%fit_weight_T, print_anh=need_analyze_anh_pot)
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
           ABI_ERROR(message)
         end if
       else
           write(message, '(3a)' )&
&           'There is no test-set file ',ch10,&
&           'Action: add test-set file'
           ABI_ERROR(message)
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
     call fit_polynomial_coeff_testEffPot(reference_effective_potential,hist_tes,master,comm,&
&                                   print_anharmonic=need_analyze_anh_pot,scup_dtset=inp%scup_dtset,&
&                                         prt_ph=inp%test_prt_ph)



   end if ! End if(inp%test_effpot == 1)then

    !TEST_AM
    !Effective Hamiltonian, compute the energy for given patern
    ! call mover_effpot(inp,filnam,reference_effective_potential,-2,comm,hist=hist)
    !TEST_AM

    !****************************************************************************************

    !****************************************************************************************
    !Print the effective potential system + coefficients (only master CPU)
    if(iam_master) then
       if (inp%prt_model >= 1) then
          write(message, '(a,(80a),a)' ) ch10,&
               &       ('=',ii=1,80)
          call wrtout(ab_out,message,'COLL')
          call wrtout(std_out,message,'COLL')
          !name = replace(trim(filnam(2)),".out","")
          ! Assume new .abo convention
          name = replace(trim(filnam(2)),".abo","")
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
    !      ABI_MALLOC(dynmat,(2,3,natom_sp,3,natom_sp))
    !      call effective_potential_effpot2dynmat(dynmat,inp%delta_df,reference_effective_potential,&
    ! &                                           reference_effective_potential%supercell%natom_supercell,&
    ! &                                           int(reference_effective_potential%supercell%qphon),3)

    !      ABI_FREE(dynmat)
    !    end if
    ! end if
    !TEST_AM SECTION


    ! Run lattice dynamics (relaxation or molecular dynamics, most of abinits ionmovs are allowed)
    !****************************************************************************************
    if(inp%dynamics>=1) then
       call mover_effpot(inp,filnam,reference_effective_potential,inp%dynamics,comm)
    end if

    !****************************************************************************************



    !Free the effective_potential and dataset
    !****************************************************************************************
    call effective_potential_free(reference_effective_potential)
    call multibinit_dtset_free(inp)
    call abihist_free(hist)
    call abihist_free(hist_tes)
!****************************************************************************************

  end subroutine multibinit_main
  !!***

end module m_multibinit_driver
