!!****m* ABINIT/m_compute_anharmonics
!! NAME
!!  m_compute_anharmonics
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group ()
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

module m_compute_anharmonics

 implicit none

 private
!!***

 public :: compute_anharmonics
!!***

contains
!!***

!!****f* ABINIT/compute_anharmonics
!!
!! NAME
!! compute_anharmonics
!!
!! FUNCTION
!! Compute strain phonon coupling by finite differences
!! Return the effective_potential with the third order
!!
!! INPUTS
!! filenames(17) = path with all name files
!! inp <type(multibinit_dtset_type)> = datatype with all the input variables
!! comm=MPI communicator
!!
!! OUTPUT
!! eff_pot<type(effective_potential_type)> = effective_potential datatype to be initialized
!!
!! PARENTS
!!      m_multibinit_driver
!!
!! CHILDREN
!!      effective_potential_file_read,effective_potential_free
!!      effective_potential_setelastic3rd,effective_potential_setelastic4th
!!      effective_potential_setelasticdispcoupling
!!      effective_potential_setstrainphononcoupling
!!      effective_potential_writeabiinput,harmonics_terms_applysumrule
!!      phonon_strain,strain_free,strain_get,strain_init,wrtout,xmpi_bcast
!!
!! SOURCE

subroutine compute_anharmonics(eff_pot,filenames,inp,comm)

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_io_tools, only : open_file

 use m_ifc
 use m_anharmonics_terms
 use m_effective_potential
 use m_effective_potential_file
 use m_multibinit_dataset, only : multibinit_dtset_type
 use m_strain
 use m_fstrings, only : itoa,int2char4,ftoa
  implicit none

 !Arguments ------------------------------------
 !scalars
  integer, intent(in) :: comm
  character(len=fnlen),intent(in) :: filenames(17)
  type(effective_potential_type),target, intent(inout) :: eff_pot
  type(multibinit_dtset_type),intent(in) :: inp
 !arrays

 !Local variables-------------------------------
 !scalar
  integer :: ia,ii,ierr,irpt,jj,kk,my_rank,natom
  integer :: nfile,nrpt,nproc
  real(dp) :: delta,delta1,delta2
  character(len=500) :: message
  character(len=fnlen):: name
  logical :: files_availables = .True.,has_any_strain = .False.
  logical :: has_all_strain = .True.
  logical :: iam_master=.FALSE.
  integer,parameter :: master=0
 !arrays
  integer  :: have_strain(6)
  real(dp) :: deformation(6,2),elastics3rd(6,6,6)
  real(dp) :: elastics4th(6,6,6,6),rprimd_def(3,3)
  type(strain_type) :: strain
  type(ifc_type) :: phonon_strain(6)
  logical, allocatable :: file_usable(:)
  real(dp),allocatable :: elastic_displacement(:,:,:,:)
  type(effective_potential_type),dimension(:),allocatable :: eff_pots
  type(strain_type),dimension(:),allocatable :: effpot_strain
  type(effective_potential_type),pointer :: ref_eff_pot

 ! *************************************************************************

  write(message,'(a,(80a),a)') ch10,('=',ii=1,80),ch10
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

  write(message, '(a,a,a)' )' Compute the third order derivative by finite differences',ch10
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  write(message, '(a,a,a)' )' The following files will be used :'
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')

 !==========================================
 !0)Initialisation of variables:
! Set MPI local varibaless
  nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
  iam_master = (my_rank == master)

 !==========================================
 !1) Get the list of files
  nfile = 0
  jj=6
  do while (jj < 18)
    if (filenames(jj)/="") then
      if(jj==6) nfile = 0
      write(message, '(a,a)' )'  - ',trim(filenames(jj))
      call wrtout(std_out,message,'COLL')
      call wrtout(ab_out,message,'COLL')
      jj = jj + 1
      nfile = nfile + 1
    else
      exit
    end if
  end do

  if(nfile==0) then
    write(message,'(a)') '  - No file found -'
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')
  end if

  write(message,'(a,(80a),a)') ch10,('-',ii=1,80),ch10
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

 !============================================
 !2) Read the effectives potential from files"
 !   - store the reference effective potential
 !   - Also get the strain
 !   - perform some checks
  ABI_DATATYPE_ALLOCATE(eff_pots,(nfile))
  ABI_DATATYPE_ALLOCATE(effpot_strain,(nfile))
  ABI_ALLOCATE(file_usable,(nfile))

  ref_eff_pot => eff_pot
  file_usable(:) = .True.

  ii = 1 ! Start at the index 1
  jj = 6 ! Start at the index 6
  do while (jj < 18)
    if (filenames(jj)/="".and.filenames(jj)/="no") then
      !Read and Intialisation of the effective potential type
      call effective_potential_file_read(filenames(jj),eff_pots(ii),inp,comm)
      !Eventualy print the xml file
!      if(inp%prt_model==-1.or.inp%prt_model>=3) then
!        call int2char4(ii,message)
!        name = 'structure_'//trim(itoa(ii-1))//'.xml'
!        call isfile(name,'new')
!        call effective_potential_writeXML(eff_pots(ii),1,filename=name)
!      end if

      !Fill the eff_pots with the conresponding strain
      call strain_get(effpot_strain(ii),rprim=eff_pot%crystal%rprimd,&
&                     rprim_def=eff_pots(ii)%crystal%rprimd)

      jj = jj + 1; ii = ii + 1

      write(message,'(a,(80a))') ch10,('-',ia=1,80)
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')
    else
      exit
    end if
  end do

  !Do some checks
  if(iam_master)then
    do ii=1,size(eff_pots)
      if (eff_pots(ii)%harmonics_terms%ifcs%nrpt/=ref_eff_pot%harmonics_terms%ifcs%nrpt) then
        write(message,'(a,I0,a,a,a,a,a,I0,a,a,a,a)' )&
&      'the number of cell in reference  (',ref_eff_pot%harmonics_terms%ifcs%nrpt,&
&       ') is not equal to the  ',ch10,'the number of cell  in ',trim(filenames(ii+5)),&
&      ' (',eff_pots(ii)%harmonics_terms%ifcs%nrpt,')',ch10,'this files cannot be used',ch10
        ABI_WARNING(message)
        file_usable(ii) = .False.
      end if
      if (eff_pots(ii)%crystal%natom/=ref_eff_pot%crystal%natom) then
        write(message, '(a,I0,a,a,a,a,a,I0,a,a,a,a)' )&
&      'the number of atoms in reference  (',ref_eff_pot%crystal%natom,') is not equal to the  ',ch10,&
&      'the number of atoms  in ',trim(filenames(ii+5)),' (',eff_pots(ii)%crystal%natom,')',ch10,&
&      'this files cannot be used',ch10
        ABI_WARNING(message)
        file_usable(ii) = .False.
      end if
      if (eff_pots(ii)%crystal%ntypat/=ref_eff_pot%crystal%ntypat) then
        write(message, '(a,I0,a,a,a,a,a,I0,a,a,a,a)' )&
&      'the number of type of atoms in reference  (',ref_eff_pot%crystal%ntypat,&
&       ') is not equal to the  ',&
&       ch10,'the number of type of atoms  in ',trim(filenames(ii+5)),&
&       ' (',eff_pots(ii)%crystal%ntypat,')',&
&       ch10,'this files can not be used',ch10
        ABI_WARNING(message)
        file_usable(ii) = .False.
      end if
    end do
  end if

! MPI BROADCAST
  do ii=1,size(eff_pots)
    call xmpi_bcast (file_usable(ii), master, comm, ierr)
  end do

  if (count((effpot_strain%name=="reference"))>1) then
    write(message, '(2a)' )&
&    ' There is several file corresponding to the reference ',ch10
    ABI_BUG(message)
  end if

  have_strain = 0

  write(message,'(a)') ' Strains available after reading the files:'
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')
  do ii=1,size(eff_pots)
    if(effpot_strain(ii)%name /= "".and.file_usable(ii)) then
      write(message,'(a,a,a,I2,a,(ES10.2),a)')&
&       ' A ',trim(effpot_strain(ii)%name),' strain in the direction ',&
&       effpot_strain(ii)%direction,' with delta of ',effpot_strain(ii)%delta
      has_any_strain = .True.
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')
    end if
  end do


  if(nfile>1.and.has_any_strain) then
    write(message,'(a,a,a)') ch10, ' ---analize in more details these files---',ch10
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')
  else
    write(message,'(a)') '  - No strain found -'
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')
    write(message,'(a,(80a),a)') ch10,('-',ia=1,80),ch10
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')
  end if

 !First check the strain
  do ii =1,6
    jj = 0
    jj = count(effpot_strain%direction==ii)
    if(jj>2) then
      write(message, '(a,I1,a)' )&
 &    ' There is several file corresponding to strain uniaxial in direction ',ii,ch10
      ABI_ERROR(message)
    else
      name = 'uniaxial'
      if(ii>=4) name = 'shear'
      if (jj==1) then
        write(message, '(a,a,a,I1,a,a)' )&
&       ' WARNING: There is only one strain ',trim(name),' in direction ',ii,ch10,&
&       '          the finate diferences will not be centering'
        call wrtout(std_out,message,"COLL")
        has_all_strain = .False.
        have_strain(ii)=jj
      else
        if(jj==2)then
          write(message, '(a,a,a,I1,a)' )&
&          ' There is two files corresponding to strain ',trim(name),' in direction ',ii,ch10
          call wrtout(ab_out,message,'COLL')
          call wrtout(std_out,message,'COLL')
          have_strain(ii)=jj
        else
          write(message, '(a,a,a,I1,a,a)' )&
&        ' WARNING: There is no strain ',trim(name),' in direction ',ii,ch10
          call wrtout(std_out,message,"COLL")
          has_all_strain = .False.
          if (inp%strcpling == 2) then
            do kk = 1,2
              delta = inp%delta_df
              if (kk==1) delta = -1 * delta
              call strain_init(strain,name=name,direction=ii,delta=delta)
              rprimd_def = matmul(eff_pot%crystal%rprimd,strain%strain)
              if(kk==1) then
                write(message, '(a,a,a,a,a,I1,a,a,a,a)' )&
&                 ' if you want to get the correct structure, please run dfpt calculation with',ch10,&
&                 ' strain ',trim(name),' in the direction',ii,' with delta=',trim(ftoa(delta)),ch10,&
&                 ' The corresponding primitive vectors are:'
              else
                write(message, '(a,a,a,I1,a,a,a,a)' )&
&                 ' And a strain ',trim(name),' in the direction',ii,' with delta = ',&
&                 trim(ftoa(delta)),ch10,' The corresponding primitive vectors are:'
              end if
              call wrtout(ab_out,message,'COLL')
              call wrtout(std_out,message,'COLL')
              write(message,'(3(F20.10),a,3(F20.10),a,3(F20.10))')&
&               rprimd_def(:,1),ch10, rprimd_def(:,2), ch10,rprimd_def(:,3)
              call wrtout(ab_out,message,'COLL')
              call wrtout(std_out,message,'COLL')
              if(iam_master)then
                call effective_potential_writeAbiInput(eff_pot,strain=strain)
              end if
              call strain_free(strain)
            end do
          end if
        end if
      end if
    end if
  end do

! check if strain exist
  if(all(have_strain==0).and.inp%strcpling /= 2) then
      write(message, '(6a)' )&
&    ' WARNING: There is no file corresponding to strain',&
&    ' to compute 3rd order derivatives.',ch10,&
&    '          In this case the 3rd order derivatives are not set',ch10,&
&    '          (add files or set strcpling to 0)'
      call wrtout(std_out,message,"COLL")
  end if

! check if the existing strains have opposite deformation
  deformation = zero
  do ii=1,6
    if(have_strain(ii)/=0) then
      ia = 1
      do jj=1,size(eff_pots)
        if (effpot_strain(jj)%direction==ii)then
          deformation(ii,ia) = effpot_strain(jj)%delta
          ia = ia + 1
        end if
      end do
      if (have_strain(ii)==2)  then
        delta1 = deformation(ii,1)
        delta2 = deformation(ii,2)
        if (delta1+delta2 > tol15) then
          write(message, '(a,I1,a,a)' )&
&             ' The deformations for strain ',ii,&
&             ' are not the opposite',ch10
          ABI_ERROR(message)
        end if
      end if
    end if
  end do

  write(message,'(a,(80a))') ch10,('-',ia=1,80)
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

  write(message,'(a,a)') ch10, ' After analyzing, the strains available are:'
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')
  if(has_any_strain) then
    do ii=1,6
      if(have_strain(ii)/=0) then
        do jj=1,size(eff_pots)
          if (effpot_strain(jj)%direction==ii)then
            write(message,'(a,a,a,I2,a,(ES10.2),a)')&
&             ' A ',trim(effpot_strain(jj)%name),' strain in the direction ',&
&             effpot_strain(jj)%direction,' with delta of ',effpot_strain(jj)%delta
            call wrtout(ab_out,message,'COLL')
            call wrtout(std_out,message,'COLL')
          end if
        end do
      else
        files_availables = .False.
      end if
    end do
  else
    files_availables = .False.
    write(message,'(a)') '  - No strain available -'
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')
  end if
  write(message,'(a,(80a))') ch10,('-',ia=1,80)
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

  if(has_all_strain) then
    write(message,'(3a)') ch10, ' The computation of the third order derivative ',&
&    'is possible'
  else
    if (inp%strcpling /= 2) then
      if(ref_eff_pot%has_anharmonicsTerms)then
        write(message,'(10a)') ch10, ' The computation of the third order derivative ',&
&        'is not possible',ch10,' somes files are missing please use strcpling 2 to generate',&
&        ' inputs files',ch10,' usable by abinit. The third order derivatives  present in  ',&
&        trim(filenames(3)),' will be used'
      else
        write(message,'(9a)') ch10, ' The computation of the third order derivative ',&
&        'is not possible',ch10,' somes files are missing please use strcpling 2 to generate',&
&        ' inputs files',ch10,' usable by abinit. The third order derivative will not be set in',&
&        ' the XML file'
      end if
    else
      if(ref_eff_pot%has_anharmonicsTerms)then
        write(message,'(10a)') ch10, ' The computation of the third order derivative ',&
&      'is not possible',ch10,' somes files are missing, the input files usable by abinit have been',&
&      ' generate.',ch10,' The third order derivatives present in ',trim(filenames(3)),' will be used'
      else
        write(message,'(8a)') ch10, ' The computation of the third order derivative ',&
&      'is not possible',ch10,' somes files are missing, the input files usable by abinit have been',&
&      ' generate.',ch10,' The third order derivatives will be not set in the XML file'
      end if
    end if
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')
  end if

 !================================================
 !3) Compute finate differences
  if(has_all_strain) then

!   Allocation of array and set some values
    nrpt  = ref_eff_pot%harmonics_terms%ifcs%nrpt
    natom = ref_eff_pot%crystal%natom
    ABI_ALLOCATE(elastic_displacement,(6,6,3,natom))

    elastics3rd = zero
    elastics4th = zero

    do ii=1,6
      if(have_strain(ii)/=0) then
 !      We want the find the index of the perturbation ii in eff_pots(ii)
 !      And store in delta1 and delta2
        delta1 = zero
        delta2 = zero
        do jj=1,size(eff_pots)
          if (effpot_strain(jj)%direction==ii.and.(effpot_strain(jj)%direction/=0))then
            if (abs(delta1)<tol16) then
              delta1 = jj
            else
              delta2 = jj
            end if
          end if
        end do
        if (abs(delta1)>tol16.and.abs(delta1)>tol16)then
 !        check if delta1 < delta2, in this case, inverse delta1 and delta2
          if (effpot_strain(int(delta1))%delta < effpot_strain(int(delta2))%delta) then
            delta = delta1
            delta1 = delta2
            delta2 = delta
          end if
!         Compute strain phonon-coupling
          phonon_strain(ii)%nrpt =  nrpt
          ABI_ALLOCATE(phonon_strain(ii)%atmfrc,(3,natom,3,natom,nrpt))
          ABI_ALLOCATE(phonon_strain(ii)%cell,(3,nrpt))
          phonon_strain(ii)%atmfrc = zero
          phonon_strain(ii)%cell =  eff_pots(int(delta1))%harmonics_terms%ifcs%cell

          do irpt=1,phonon_strain(ii)%nrpt
            phonon_strain(ii)%atmfrc(:,:,:,:,irpt) =&
&           (eff_pots(int(delta1))%harmonics_terms%ifcs%atmfrc(:,:,:,:,irpt)&
&          - eff_pots(int(delta2))%harmonics_terms%ifcs%atmfrc(:,:,:,:,irpt)) / &
&            (2 * abs(effpot_strain(int(delta1))%delta))
          end do

          if(inp%asr >= 0) then
!           Impose sum rule
            call harmonics_terms_applySumRule(inp%asr,phonon_strain(ii),&
&                                                 eff_pot%crystal%natom)
          end if

!         Compute elastic constants
          elastics3rd(ii,:,:) = (eff_pots(int(delta1))%harmonics_terms%elastic_constants(:,:)&
&          - eff_pots(int(delta2))%harmonics_terms%elastic_constants(:,:)) / &
&            (2 * abs(effpot_strain(int(delta1))%delta))

!         Compute elastic-displacement coupling
          elastic_displacement(ii,:,:,:)=(eff_pots(int(delta1))%harmonics_terms%strain_coupling(:,:,:)&
&          - eff_pots(int(delta2))%harmonics_terms%strain_coupling(:,:,:)) / &
&            (2 * abs(effpot_strain(int(delta1))%delta))

!         Compute elastic constants
          elastics4th(ii,ii,:,:) = (eff_pots(int(delta1))%harmonics_terms%elastic_constants(:,:)&
&            - 2*ref_eff_pot%harmonics_terms%elastic_constants(:,:)&
&          + eff_pots(int(delta2))%harmonics_terms%elastic_constants(:,:)) / &
&            (abs(effpot_strain(int(delta1))%delta)**2)
        end if

      end if
    end do

!   Set all the values in the effective potential type
    call effective_potential_setStrainPhononCoupling(eff_pot,natom,phonon_strain)
    call effective_potential_setElastic3rd(eff_pot,elastics3rd)
    call effective_potential_setElastic4th(eff_pot,elastics4th)
    call effective_potential_setElasticDispCoupling(eff_pot,natom,elastic_displacement)


!   Free the phonon-strain coupling array
    do ii = 1,6
      call phonon_strain(ii)%free()
    end do
    ABI_DEALLOCATE(elastic_displacement)

    write(message,'(4a)') ch10, ' The computation of the 3rd order elastics constants, ',ch10,&
&    ' the phonon-strain coupling and the elastic-displacement coupling is done'
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')

  end if


 !===============================================
 !4) Free the array of effective potential

  do jj=1,nfile
 !  Free the effective potential type
    call effective_potential_free(eff_pots(jj))
  end do

  ABI_DATATYPE_DEALLOCATE(effpot_strain)
  ABI_DATATYPE_DEALLOCATE(eff_pots)
  ABI_DEALLOCATE(file_usable)


  write(message,'(a,a,a,(80a))') ch10,('=',ii=1,80),ch10
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

end subroutine compute_anharmonics
!!***

end module m_compute_anharmonics
!!***
