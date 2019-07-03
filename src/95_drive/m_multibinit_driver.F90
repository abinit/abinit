
!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_multibinit_main
!! NAME
!! m_multibinit_main
!!
!! FUNCTION
!! Main routine MULTIBINIT.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2019 ABINIT group (AM, hexu)
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
!! PARENTS
!!
!! CHILDREN
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
  use m_build_info
  use m_xmpi
  use m_xomp
  use m_abicore
  use m_errors

  use m_effective_potential
  use m_fit_polynomial_coeff
  use m_multibinit_dataset
  use m_effective_potential_file
  !use m_spin_model, only: spin_model_t
  use m_abihist

  use m_multibinit_manager, only: mb_manager_t
  use m_multibinit_main2, only: multibinit_main2

  use m_mover_effpot, only : mover_effpot
  !use m_generate_training_set, only : generate_training_set
  use m_compute_anharmonics, only : compute_anharmonics
  use m_init10,              only : init10
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
  !! PARENTS
  !!
  !!
  !! CHILDREN
  !!
  !!
  !! SOURCE
  subroutine multibinit_main(filnam, dry_run)
    character(len=fnlen), intent(inout) :: filnam(17)
    integer, intent(in) :: dry_run
    type(multibinit_dtset_type), target :: inp
    type(effective_potential_type) :: reference_effective_potential
    type(abihist) :: hist

    ! data for spin
    !type(spin_model_t) :: spin_model
    type(mb_manager_t) :: manager
    character(len=strlen) :: string
    character(len=500) :: message
    character(len=fnlen) :: name

    integer :: filetype,ii,lenstr
    integer :: natom,nph1l,nrpt,ntypat
    integer :: option

    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master


    master = 0
    comm = xmpi_world
    nproc = xmpi_comm_size(comm)
    my_rank = xmpi_comm_rank(comm)
    iam_master = (my_rank == master)


    !To automate a maximum calculation, multibinit reads the number of atoms
    !in the file (ddb or xml). If DDB file is present in input, the ifc calculation
    !will be initilaze array to the maximum of atoms (natifc=natom,atifc=1,natom...) in invars10

    !TODO: hexu comment: why no if(iam_master) ?
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

    if(dry_run/=0) then
       call wrtout([std_out, ab_out], "Multibinit in dry_run mode. Exiting after input parser")
       call xmpi_end()
       !goto 100
    endif

    ! Read and treat the reference structure
    !****************************************************************************************
    if (inp%spin_dynamics>0) then
       if (iam_master) then
          write(message,'(a,(80a),3a)') ch10,('=',ii=1,80),ch10,ch10,&
               &     'reading spin terms.'
       end if
       !call spin_model%initialize( filnam, inp )
       call  manager%initialize(filnam, params=inp)
    else
       !  Read the model (from DDB or XML)
       call effective_potential_file_read(filnam(3),reference_effective_potential,inp,comm)

       !Read the coefficient from fit
       !FIXME: hexu: on test farm, it is not no but $path/no
       if(filnam(4)/=''.and.filnam(4)/='no')then
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
             MSG_ERROR(message)

          else
             write(message,'(a,(80a),3a)') ch10,('=',ii=1,80),ch10,ch10,&
                  &         ' There is no file for the coefficients from polynomial fitting'
             call wrtout(ab_out,message,'COLL')
             call wrtout(std_out,message,'COLL')
          end if
       end if
    end if

    !****************************************************************************************

    ! Compute the third order derivative with finite differences
    !****************************************************************************************
    if (inp%strcpling > 0) then
       call compute_anharmonics(reference_effective_potential,filnam,inp,comm)
    end if
    !****************************************************************************************

    ! If needed, fit the anharmonic part and compute the confinement potential
    !****************************************************************************************
    if (inp%fit_coeff/=0.or.inp%confinement==2.or.inp%bound_model/=0) then

       if(iam_master) then
          !    Read the MD file
          write(message,'(a,(80a),7a)')ch10,('=',ii=1,80),ch10,ch10,&
               &       '-Reading the file the HIST file :',ch10,&
               &       '-',trim(filnam(5)),ch10

          call wrtout(std_out,message,'COLL')
          call wrtout(ab_out,message,'COLL')
          if(filnam(5)/=''.and.filnam(5)/='no')then
             call effective_potential_file_readMDfile(filnam(5),hist,option=inp%ts_option)
             if (hist%mxhist == 0)then
                write(message, '(5a)' )&
                     &           'The MD ',trim(filnam(5)),' file is not correct ',ch10,&
                     &           'Action: add MD file'
                MSG_ERROR(message)
             end if
          else
             if (inp%fit_coeff/=0) then
                write(message, '(3a)' )&
                     &           'There is no MD file to fit the coefficients ',ch10,&
                     &           'Action: add MD file'
                MSG_ERROR(message)
             else
                if (inp%bound_model/=0) then
                   write(message, '(3a)' )&
                        &             'There is no MD file to bound the model ',ch10,&
                        &             'Action: add MD file'
                   MSG_ERROR(message)
                else if(inp%confinement==2) then
                   write(message, '(3a)' )&
                        &             'There is no MD file to compute the confinement',ch10,&
                        &             'Action: add MD file'
                   MSG_ERROR(message)
                end if
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
                  &         spcoupling=inp%fit_SPCoupling==1,prt_names=inp%prt_names)
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
    end if

    !****************************************************************************************

    !****************************************************************************************
    !TEST_AM
    !Effective Hamiltonian, compute the energy for given patern
    ! call mover_effpot(inp,filnam,reference_effective_potential,-2,comm,hist=hist)
    !TEST_AM

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
    if(inp%spin_dynamics/=0) then
       !call spin_model%run()
       call manager%run()
    end if
    !****************************************************************************************


    !Free the effective_potential and dataset
    !****************************************************************************************
    if(inp%spin_dynamics/=0) then
       call manager%finalize()
    end if

    call effective_potential_free(reference_effective_potential)
    call multibinit_dtset_free(inp)
    call abihist_free(hist)

  end subroutine multibinit_main
  !!***

end module m_multibinit_driver

