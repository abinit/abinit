!{\src2tex{textfont=tt}}
!!****f* ABINIT/strain_phonon_coupling
!!
!! NAME
!! strain_phonon_coupling
!!
!! FUNCTION
!! Compute strain phonon coupling by finite differences
!! Return the effective_potential with the third order
!!
!! INPUTS
!! filenames = input with all name files
!! inp  = input of epigene
!! comm=MPI communicator
!!
!! OUTPUT
!! effective_potantial = effective_potential structure with 3rd orders
!!
!! PARENTS
!!    epigene
!!
!! CHILDREN
!!    check_effpot,effective_potential_check_strain,file_to_effective_potential,get_strain
!!    free_effective_potential,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine strain_phonon_coupling(effective_potential,filenames,inp,comm)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_xmpi

 use m_ddb
 use m_effective_potential
 use m_effective_potential_file
 use m_epigene_dataset, only : epigene_dataset_type
 use m_strain
 use m_fstrings, only : itoa,int2char4,ftoa

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'strain_phonon_coupling'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(inout) :: comm
 character(len=fnlen),intent(in) :: filenames(15)
 type(effective_potential_type), intent(inout) :: effective_potential
 type(epigene_dataset_type),intent(in) :: inp
!arrays

!Local variables-------------------------------
!scalar
 integer :: ia,ii,irpt,jj,kk,nfile
 real(dp) :: delta,delta1,delta2
 character(len=500) :: message
 character(len=fnlen):: name
 logical :: files_availables = .True.
 logical :: has_any_strain = .False.
 logical :: has_all_strain = .True.
!arrays
 type(effective_potential_type),dimension(:),allocatable :: eff_pots
 type(strain_type) :: strain
 integer  :: have_strain(6)
 real(dp) :: deformation(6,2),rprimd_def(3,3)
 logical, allocatable  :: file_usable(:)
! *************************************************************************

 write(message,'(a,a,(80a))') ('=',ii=1,80),ch10
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


!==========================================
!1) Get the list of files
 nfile = 1
 jj = 4
 do while (jj < 16) 
   if (filenames(jj)/="") then 
     if(jj==4) nfile = 1
     write(message, '(a,a)' )'  - ',trim(filenames(jj))
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL') 
     jj = jj + 1
     nfile = nfile + 1
   else 
     exit 
   end if
 end do
 
 if(nfile==1) then
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
 ii = 2
 jj = 4

 ABI_DATATYPE_ALLOCATE(eff_pots,(nfile))
 ABI_ALLOCATE(file_usable,(nfile))

 eff_pots(1) = effective_potential
 file_usable(:) = .True.

 do while (jj < 16) 
   if (filenames(jj)/="") then
     call effective_potential_file_read(filenames(jj),eff_pots(ii),inp,comm)
     !Intialisation of the effective potential type

     !Eventualy print the xml file
     if(inp%prt_effpot==-1.or.inp%prt_effpot>=3) then
       call int2char4(ii,message)
       name = 'structure_'//trim(itoa(ii-1))//'.xml'
       call isfile(name,'new')
       call effective_potential_writeXML(eff_pots(ii),1,filename=name)
     end if

     !Fill the eff_pots with the conresponding strain
     call strain_get(eff_pots(ii)%strain,rprim=effective_potential%rprimd,rprim_def=eff_pots(ii)%rprimd)
   
     jj = jj + 1; ii = ii + 1
   
     write(message,'(a,(80a))') ch10,('-',ia=1,80)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
   else
     exit
   end if
 end do

 !Do some checks
 do ii=1,size(eff_pots)
   if (eff_pots(ii)%ifcs%nrpt/=eff_pots(1)%ifcs%nrpt) then
     write(message,'(a,I5,a,a,a,a,a,I5,a,a,a,a)' )&
&   'the number of cell in reference  (',eff_pots(1)%ifcs%nrpt,') is not equal to the  ',ch10,&
&   'the number of cell  in',trim(filenames(ii+2)),' (',eff_pots(ii)%ifcs%nrpt,')',ch10,&
&   'this files cannot be used',ch10
     MSG_WARNING(message)
     file_usable(ii) = .False.
   end if
   if (eff_pots(ii)%natom/=eff_pots(1)%natom) then
     write(message, '(a,I5,a,a,a,a,a,I5,a,a,a,a)' )&
&   'the number of atoms in reference  (',eff_pots(1)%natom,') is not equal to the  ',ch10,&
&   'the number of atoms  in',trim(filenames(ii+2)),' (',eff_pots(ii)%natom,')',ch10,&
&   'this files cannot be used',ch10
     MSG_WARNING(message)
     file_usable(ii) = .False.
   end if
   if (eff_pots(ii)%ntypat/=eff_pots(1)%ntypat) then
     write(message, '(a,I5,a,a,a,a,a,I5,a,a,a,a)' )&
&   'the number of type of atoms in reference  (',eff_pots(1)%ntypat,') is not equal to the  ',ch10,&
&   'the number of type of atoms  in',trim(filenames(ii+2)),' (',eff_pots(ii)%ntypat,')',ch10,&
&   'this files can not be used',ch10
     MSG_WARNING(message)
     file_usable(ii) = .False.
   end if
 end do
 
 if (count((eff_pots%strain%name=="reference"))>1) then
   write(message, '(a,a)' )&
&   ' There is several file corresponding to the reference ',ch10
   MSG_BUG(message)
 end if  

 have_strain = zero

 write(message,'(a,a)') ch10, ' Strains available after reading the files after:'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 do ii=1,size(eff_pots)
   if(eff_pots(ii)%strain%name /= "".and.file_usable(ii)) then 
     write(message,'(a,a,a,I2,a,(ES10.2),a)')&
&      ' A ',trim(eff_pots(ii)%strain%name),' strain in the direction ',&
&      eff_pots(ii)%strain%direction,' with delta of ',eff_pots(ii)%strain%delta
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
   write(message,'(a,(80a))') ch10,('-',ia=1,80)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

!First check the strain  
 do ii =1,6
   jj = zero
   jj = count(eff_pots%strain%direction==ii)
   if(jj>2) then
     write(message, '(a,I1,a)' )&
 &   ' There is several file corresponding to strain uniaxial in direction ',ii,ch10
     MSG_ERROR(message)
   else
     name = 'uniaxial'
     if(ii>=4) name = 'shear'
     if (jj==1) then
       write(message, '(a,a,a,I1,a,a)' )&
&   ' There is only one strain ',trim(name),' in direction ',ii,ch10,&
&   'the finate diferences will not be centering'
       MSG_WARNING(message)
       has_all_strain = .False.
       have_strain(ii)=jj
     else
       if(jj==2)then 
         write(message, '(a,a,a,I1,a)' )&
&         ' There is two files corresponding to strain ',trim(name),' in direction ',ii,ch10
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')
         have_strain(ii)=jj
       else
         write(message, '(a,a,a,I1,a,a)' )&
&       ' There is no strain ',trim(name),' in direction ',ii,ch10
         MSG_WARNING(message)
         has_all_strain = .False.
         if (inp%prt_3rd == 2) then
           do kk = 1,2 
             delta = inp%delta_df
             if (kk==1) delta = -1 * delta 
             call strain_init(strain,name=name,direction=ii,delta=delta)
             rprimd_def = matmul(effective_potential%rprimd,strain%strain)
             if(kk==1) then
               write(message, '(a,a,a,a,a,I1,a,a,a,a)' ) &
&                ' if you want to get the correct structure, please run dfpt calculation with',ch10,&
&                ' strain ',trim(name),' in the direction',ii,' with delta=',trim(ftoa(delta)),ch10,&
&                ' The corresponding primitive vectors are:'
             else
               write(message, '(a,a,a,I1,a,a,a,a)' ) &
&                ' And a strain ',trim(name),' in the direction',ii,' with delta = ',&
&                trim(ftoa(delta)),ch10,' The corresponding primitive vectors are:'               
             end if
             call wrtout(ab_out,message,'COLL')
             call wrtout(std_out,message,'COLL')
             write(message,'(3(F20.10),a,3(F20.10),a,3(F20.10))')&
&              rprimd_def(:,1),ch10, rprimd_def(:,2), ch10,rprimd_def(:,3)
             call wrtout(ab_out,message,'COLL')
             call wrtout(std_out,message,'COLL')
             call effective_potential_writeAbiInput(effective_potential,strain=strain)
             call strain_free(strain)
           end do
         end if
       end if
     end if
   end if
 end do

!check if strain exist
 if(all(have_strain==zero).and.inp%prt_3rd /= 2) then
     write(message, '(a,a,a,a,a,a,a)' )&
&   ' There is no file corresponding to strain',&
&   ' to compute 3rd order derivatives.',ch10,&
&   'In this case the 3rd order derivatives are not set',ch10,&
&   '(add files or set prt_3rd to 0)',ch10
     MSG_WARNING(message)   
 end if

!check if the existing strains have opposite deformation
 deformation = zero
 do ii=1,6
   if(have_strain(ii)/=0) then
     ia = 1
     do jj=1,size(eff_pots)
       if (eff_pots(jj)%strain%direction==ii)then
         deformation(ii,ia) = eff_pots(jj)%strain%delta
         ia = ia + 1
       end if
     end do
     if (have_strain(ii)==2)  then 
        delta1 = deformation(ii,1) 
        delta2 = deformation(ii,2)
        if (delta1/=-delta2) then
            write(message, '(a,I1,a,a)' )&
&            ' The deformations for strain ',ii,&
&            ' are not the opposite',ch10
            MSG_ERROR(message)   
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
         if (eff_pots(jj)%strain%direction==ii)then
           write(message,'(a,a,a,I2,a,(ES10.2),a)')&
&            ' A ',trim(eff_pots(jj)%strain%name),' strain in the direction ',&
&            eff_pots(jj)%strain%direction,' with delta of ',eff_pots(jj)%strain%delta
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
   write(message,'(a,a,a)') ch10, ' The computation of the third order derivative ',&
&   'is possible'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 else
   if (inp%prt_3rd /= 2) then
     write(message,'(a,a,a,a,a,a,a,a,a)') ch10, ' The computation of the third order derivative ',&
&     'is not possible',ch10,' somes files are missing please use prt_3rd 2 to generate',&
&     ' inputs files',ch10,' usable by abinit The third order derivative will not be set in',&
&     ' the XML file'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')   
   else
     write(message,'(a,a,a,a,a,a,a,a,a,a)') ch10, ' The computation of the third order derivative ',&
&   'is not possible',ch10,' somes files are missing, the input files usable by abinit have been',&
&   ' generate.',ch10,' The third order derivative will be not set in the XML file'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')   
   end if
 end if

!================================================
!3) Compute finate differences
 if(has_all_strain) then
   do ii=1,6
     if(have_strain(ii)/=0) then
!      We want the find the index of the perturbation ii in eff_pots(ii)
!      And store in delta1 and delta2 
       delta1 = zero
       delta2 = zero
       do jj=1,size(eff_pots)
         if (eff_pots(jj)%strain%direction==ii.and.(eff_pots(jj)%strain%direction/=zero))then
           if (delta1==zero) then
             delta1 = jj 
           else
             delta2 = jj
           end if
         end if
       end do
       if (delta1/=0.and.delta1/=0)then
!        check if delta1 < delta2, in this case, inverse delta1 and delta2
         if (eff_pots(int(delta1))%strain%delta < eff_pots(int(delta2))%strain%delta) then
           delta = delta1
           delta1 = delta2
           delta2 = delta
         end if
         effective_potential%phonon_strain_coupling(ii)%nrpt = eff_pots(int(delta1))%ifcs%nrpt
         effective_potential%phonon_strain_coupling(ii)%cell = eff_pots(int(delta1))%ifcs%cell
         do irpt=1,effective_potential%phonon_strain_coupling(ii)%nrpt
           effective_potential%phonon_strain_coupling(ii)%atmfrc(:,:,:,:,:,irpt) = &
&           (eff_pots(int(delta1))%ifcs%short_atmfrc(:,:,:,:,:,irpt)&
&          - eff_pots(int(delta2))%ifcs%short_atmfrc(:,:,:,:,:,irpt) ) / &
&           (2 * abs(eff_pots(int(delta1))%strain%delta))
         end do
       end if
       cycle
     end if
   end do
   effective_potential%has_3rd = .True.
 end if

!===============================================
!4) Free the array of effective potential

 do jj=1,size(eff_pots)
!  Free the effective potential type
   call effective_potential_free(eff_pots(jj))
 end do

 ABI_DATATYPE_DEALLOCATE(eff_pots)
 ABI_DEALLOCATE(file_usable)
 write(message,'(a,a,a,(80a))') ch10,('=',ii=1,80),ch10  
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

end subroutine strain_phonon_coupling
!!***
