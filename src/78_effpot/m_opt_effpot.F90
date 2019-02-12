!!****m* ABINIT/m_opt_effpot
!!
!! NAME
!! m_opt_effpot
!!
!! FUNCTION
!!      
!! 
!!
!! COPYRIGHT
!! Copyright (C) 2010-2018 ABINIT group (AM)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE



#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_opt_effpot

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_abicore
 use m_xmpi
 use m_effective_potential 
 use m_effective_potential_file, only : effective_potential_file_mapHistToRef
 use m_fit_data
 use m_fit_polynomial_coeff
 use m_polynomial_coeff
 use m_crystal,only : symbols_crystal

 implicit none

 public :: opt_effpot 
 public :: opt_effpotbound

!!****
CONTAINS 
      
!!****f* m_opt_effpot/opt_effpot 
!!
!! NAME
!! opt_effpot
!!
!! FUNCTION
!! Optimize Effective Potential by fitting the value of certain
!! coefficients while keeping the values of the others     
!!
!! INPUTS
!! eff_pot<type(effective_potential)> = effective potential
!!
!! opt_coeff(opt_ncoeff) = list of terms whose coefficients are to be
!! optimized
!!
!! hist<type(abihist)> = Training set Data(or snapshot of DFT)
!! comm = MPI communicator
!!
!! OUTPUT
!! eff_pot<type(effective_potential)> = effective potential datatype with new fitted coefficients
!!
!! PARENTS
!! multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine opt_effpot(eff_pot,opt_ncoeff,opt_coeff,hist,comm,print_anh) 

 implicit none  

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,opt_ncoeff
 type(effective_potential_type),target,intent(inout) :: eff_pot
 type(abihist),intent(inout) :: hist
!arrays 
 integer,intent(in) :: opt_coeff(opt_ncoeff)
!Logicals
 logical,optional,intent(in) :: print_anh 
!Strings 
!Local variables ------------------------------
!scalars
 integer :: ii, info,natom_sc,ntime,unit_anh1,unit_anh2
 real(dp) :: factor,mse,msef,mses
 real(dp),parameter :: HaBohr_meVAng = 27.21138386 / 0.529177249
!arrays 
 integer :: sc_size(3)
 integer :: coeff_inds(opt_ncoeff)
 type(fit_data_type) :: fit_data
 type(polynomial_coeff_type) :: my_coeffs(opt_ncoeff)
 real(dp) :: coeff_values(opt_ncoeff), coeff_init_values(opt_ncoeff)
 real(dp), allocatable :: energy_coeffs(:,:),fcart_coeffs(:,:,:,:)
 real(dp), allocatable :: strten_coeffs(:,:,:)
!Logicals
 logical :: need_print_anh,file_opened 
!Strings 
 character(len=1000) :: message
 character(len=1000) :: frmt
 character(len=fnlen) :: fn_bf='before_opt_diff', fn_af='after_opt_diff'
! *************************************************************************

 !Setting/Initializing Variables
  ntime = hist%mxhist
  natom_sc = size(hist%xred,2)
  factor   = 1._dp/natom_sc
  if(present(print_anh)) need_print_anh = print_anh
 
 !if the number of atoms in reference supercell into effpot is not correct,
 !wrt to the number of atom in the hist, we set map the hist and set the good supercell
  if (natom_sc /= eff_pot%supercell%natom) then
    call effective_potential_file_mapHistToRef(eff_pot,hist,comm,verbose=.TRUE.)
  end if

 !we get the size of the supercell in the hist file
  do ii=1,3
    sc_size(ii) = int(anint(sqrt(eff_pot%supercell%rprimd(ii,1)**2+&
 &                               eff_pot%supercell%rprimd(ii,2)**2+&
 &                               eff_pot%supercell%rprimd(ii,3)**2) / &
 &                          sqrt(eff_pot%crystal%rprimd(ii,1)**2+&
 &                               eff_pot%crystal%rprimd(ii,2)**2+&
 &                               eff_pot%crystal%rprimd(ii,3)**2)))
  end do


 !Before the fit, compute constants with fit_data_compute.
 !Conpute the strain of each configuration.
 !Compute the displacmeent of each configuration.
 !Compute the variation of the displacement due to strain of each configuration.
 !Compute fixed forces and stresse and get the standard deviation.
 !Compute Sheppard and al Factors  \Omega^{2} see J.Chem Phys 136, 074103 (2012) [[cite:Sheppard2012]].
  call fit_data_compute(fit_data,eff_pot,hist,comm,verbose=.FALSE.)
 

  if(need_print_anh) call effective_potential_writeAnhHead(eff_pot%anharmonics_terms%ncoeff,&
&                            fn_bf,eff_pot%anharmonics_terms)                  

 !Before deleting coefficients calculate MSD of initial model  
  call fit_polynomial_coeff_computeMSD(eff_pot,hist,mse,msef,mses,&
 &                                     natom_sc,ntime,fit_data%training_set%sqomega,&
 &                                     compute_anharmonic=.TRUE.,print_file=.TRUE.,filename=fn_bf)


 !  Print the standard devition of initial model 
      write(message,'(6a,ES24.16,6a,ES24.16,2a,ES24.16,2a,ES24.16,a)' )ch10,&
 &                    ' Mean Standard Deviation values of the effective-potential',ch10,&
 &                    ' with respect to the training-set before optimization (meV/atm):',&
 &               ch10,'   Energy          : ',&
 &               mse*Ha_EV*1000*factor ,ch10,&
 &                    ' Goal function values of the effective.potential',ch10,& 
 &                    ' with respect to the test-set (eV^2/A^2):',ch10,&
 &                    '   Forces+Stresses : ',&
 &               (msef+mses)*(HaBohr_meVAng)**2,ch10,&
 &                    '   Forces          : ',&
 &               msef*(HaBohr_meVAng)**2,ch10,&
 &                    '   Stresses        : ',&
 &               mses*(HaBohr_meVAng)**2,ch10
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')


 ! Write terms to my_coeffs(ii) and zero them in eff_pot      
 do ii=1,opt_ncoeff
   !Store indices for later 
   coeff_inds(ii) = ii
   !Initialize coefficients for optimizing
   call polynomial_coeff_init(coeff_values(ii),eff_pot%anharmonics_terms%coefficients(opt_coeff(ii))%nterm,&
 &                            my_coeffs(ii), eff_pot%anharmonics_terms%coefficients(opt_coeff(ii))%terms, & 
 &                            check=.TRUE.)
   !Store initial values of coefficients  
   coeff_init_values(ii) = eff_pot%anharmonics_terms%coefficients(opt_coeff(ii))%coefficient
   !Put them temporarely to zero 
   eff_pot%anharmonics_terms%coefficients(opt_coeff(ii))%coefficient = zero    
 end do  

 !Before the fit, compute constants with fit_data_compute.
 !And coefficients to be optimized put to zero 
 !Conpute the strain of each configuration.
 !Compute the displacmeent of each configuration.
 !Compute the variation of the displacement due to strain of each configuration.
 !Compute fixed forces and stresse and get the standard deviation.
 !Compute Sheppard and al Factors  \Omega^{2} see J.Chem Phys 136, 074103 (2012) [[cite:Sheppard2012]].
  call fit_data_compute(fit_data,eff_pot,hist,comm,verbose=.TRUE.)

 !After deleting coefficients calculate MSD  
  call fit_polynomial_coeff_computeMSD(eff_pot,hist,mse,msef,mses,&
 &                                     natom_sc,ntime,fit_data%training_set%sqomega,&
 &                                      compute_anharmonic=.TRUE.)


 !  Print the standard deviation after deleting
      write(message,'(6a,ES24.16,6a,ES24.16,2a,ES24.16,2a,ES24.16,a)' )ch10,&
 &                    ' Mean Standard Deviation values of the effective-potential',ch10,&
 &                    ' with respect to the training-set after deleting selected terms (meV/atm):',&
 &               ch10,'   Energy          : ',&
 &               mse*Ha_EV*1000*factor ,ch10,&
 &                    ' Goal function values of the effective.potential',ch10,& 
 &                    ' with respect to the test-set (eV^2/A^2):',ch10,&
 &                    '   Forces+Stresses : ',&
 &               (msef+mses)*(HaBohr_meVAng)**2,ch10,&
 &                    '   Forces          : ',&
 &               msef*(HaBohr_meVAng)**2,ch10,&
 &                    '   Stresses        : ',&
 &               mses*(HaBohr_meVAng)**2,ch10
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')


      
   ! Allocate necessary arrays for the fit-data 
   ABI_ALLOCATE(energy_coeffs,(opt_ncoeff,ntime))
   ABI_ALLOCATE(fcart_coeffs,(3,natom_sc,opt_ncoeff,ntime))
   ABI_ALLOCATE(strten_coeffs,(6,ntime,opt_ncoeff))
   ! Calculate forces and stresses per coefficient, which are to be optimized 
   call fit_polynomial_coeff_getFS(my_coeffs,fit_data%training_set%du_delta,&
&                                 fit_data%training_set%displacement,&
&                                 energy_coeffs,fcart_coeffs,natom_sc,eff_pot%crystal%natom,&
&                                 opt_ncoeff,ntime,sc_size,fit_data%training_set%strain,&
&                                 strten_coeffs,fit_data%training_set%ucvol,coeff_inds,opt_ncoeff)


!  call the fit process routine
!  This routine solves the linear system proposed 
!  by C.Escorihuela-Sayalero see PRB95,094115(2017) [[cite:Escorihuela-Sayalero2017]]
   call fit_polynomial_coeff_solve(coeff_values(1:opt_ncoeff),fcart_coeffs,fit_data%fcart_diff,&
&                                  energy_coeffs,fit_data%energy_diff,info,&
&                                  coeff_inds,natom_sc,opt_ncoeff,opt_ncoeff,ntime,&
&                                  strten_coeffs,fit_data%strten_diff,&
&                                  fit_data%training_set%sqomega)

  if (info /= 0 .and. all(coeff_values < tol16))then
    write(frmt,*) opt_ncoeff  
    write(message, '(2a,'//ADJUSTR(frmt)//'I4,8a)' ) ch10,&
&        '     The attempt to optimize the terms: ', opt_coeff ,ch10,&
&        '     , returned a singular solution', ch10,&
&        '     The terms could not be optimized ',ch10,&
&        '     and the effective potential has not been altered.', ch10,&
&        '     Action: Change training set or coefficients to be optimized.'
    MSG_WARNING(message)
    do ii=1,opt_ncoeff 
      eff_pot%anharmonics_terms%coefficients(opt_coeff(ii))%coefficient = coeff_init_values(ii)       
      call polynomial_coeff_free(my_coeffs(ii))
    end do 
  else
  ! Transfer new fitted values to coefficients and write them into effective potential
  ! Deallcoate temporary coefficients my_coeffs
    do ii=1,opt_ncoeff 
       eff_pot%anharmonics_terms%coefficients(opt_coeff(ii))%coefficient = coeff_values(ii)
       call polynomial_coeff_free(my_coeffs(ii))
    end do
    !Recalculate MSD of Final Model 
    
    !Conpute the strain of each configuration.
    !Compute the displacmeent of each configuration.
    !Compute the variation of the displacement due to strain of each configuration.
    !Compute fixed forces and stresse and get the standard deviation.
    !Compute Sheppard and al Factors  \Omega^{2} see J.Chem Phys 136, 074103 (2012) [[cite:Sheppard2012]].
     call fit_data_compute(fit_data,eff_pot,hist,comm,verbose=.TRUE.)
   
  
     if(need_print_anh) call effective_potential_writeAnhHead(eff_pot%anharmonics_terms%ncoeff,&
&                            fn_af,eff_pot%anharmonics_terms)                  

    !After optimization of coefficients opt_coeff recalculate MSD
     call fit_polynomial_coeff_computeMSD(eff_pot,hist,mse,msef,mses,&
     &                                     natom_sc,ntime,fit_data%training_set%sqomega,&
     &                                     compute_anharmonic=.TRUE.,print_file=.TRUE.,filename=fn_af)
     
     
     !  Print the standard deviation after optimization
          write(message,'(6a,ES24.16,6a,ES24.16,2a,ES24.16,2a,ES24.16,a)' )ch10,&
     &                    ' Mean Standard Deviation values of the effective-potential',ch10,&
     &                    ' with respect to the training-set after optimizing selected terms (meV/atm):',&
     &               ch10,'   Energy          : ',&
     &               mse*Ha_EV*1000*factor ,ch10,&
     &                    ' Goal function values of the effective.potential',ch10,& 
     &                    ' with respect to the test-set (eV^2/A^2):',ch10,&
     &                    '   Forces+Stresses : ',&
     &               (msef+mses)*(HaBohr_meVAng)**2,ch10,&
     &                    '   Forces          : ',&
     &               msef*(HaBohr_meVAng)**2,ch10,&
     &                    '   Stresses        : ',&
     &               mses*(HaBohr_meVAng)**2,ch10
          call wrtout(ab_out,message,'COLL')
          call wrtout(std_out,message,'COLL')
  end if 

 !Deallocation of fitting variables
 ABI_DEALLOCATE(energy_coeffs) 
 ABI_DEALLOCATE(fcart_coeffs)
 ABI_DEALLOCATE(strten_coeffs)
  
 if(need_print_anh)then 
  INQUIRE(FILE='before_opt_diff_anharmonic_terms_energy.dat',OPENED=file_opened,number=unit_anh1)
  if(file_opened) close(unit_anh1)
  INQUIRE(FILE='after_opt_diff_anharmonic_terms_energy.dat',OPENED=file_opened,number=unit_anh2)
  if(file_opened) close(unit_anh2)
 end if 
 ! Deallocate and delete the fit-date 
 call fit_data_free(fit_data)
end subroutine opt_effpot


!!****f* m_opt_effpot/opt_effpotbound
!!
!! NAME
!! opt_effpotbound
!!
!! FUNCTION
!! Compute and add high order terms to existing odd or negative even anharmonic terms
!! Fix the coefficient of the added new high order terms to a value such that it 
!! doesn't influence the precision of the existing anharmonic potential with respect 
!! to a relevent training set (ATTENTIION: A user must know what a relevant training set 
!! is for the system he want's to study. Typically something oscillating around its ground-state.) 
!! Finally optimize the coefficients of the orignal anharmonic terms under the presence of the 
!! added high order terms. 
!!
!! INPUTS
!! eff_pot: existing effective potential 
!! order: order for which bounding terms are generated
!!
!! 
!! OUTPUT
!! eff_pot new effective potential 
!! 
!!
!! PARENTS
!! multibinit
!!
!! CHILDREN
!! opt_effpot 
!!
!! SOURCE

subroutine opt_effpotbound(eff_pot,order,hist,comm,print_anh)

 implicit none 
         
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 type(effective_potential_type),target,intent(inout) :: eff_pot
 type(abihist),intent(inout) :: hist
!arrays 
 integer,intent(in) :: order(2)
!Logicals
 logical,optional,intent(in) :: print_anh 
!Strings 
!Local variables ------------------------------
!scalars
 integer :: i,ii, info,natom_sc,ntime,iterm,nterm,iorder,idisp,ndisp
 integer :: nterm_of_term, nterm_inloop, iterm_of_term
 real(dp) :: factor,mse_ini,msef_ini,mses_ini,mse,msef,mses,coeff_ini=0.0001
 real(dp),parameter :: HaBohr_meVAng = 27.21138386 / 0.529177249
!arrays 
 integer :: sc_size(3)
 integer,allocatable :: terms(:) 
 type(fit_data_type) :: fit_data
 !real(dp), allocatable :: energy_coeffs(:,:),fcart_coeffs(:,:,:,:)
 !real(dp), allocatable :: strten_coeffs(:,:,:)
 type(polynomial_coeff_type),allocatable :: my_coeffs(:)
!Logicals
 logical :: need_print_anh=.FALSE.,file_opened ! MARCUS FOR THE MOMENT PRINT NO FILES
!Strings
 character(len=5),allocatable :: symbols(:)
 character(len=200):: name
 character(len=1000) :: message
 character(len=1000) :: frmt
 character(len=fnlen) :: fn_bf='before_opt_diff', fn_af='after_opt_diff'
!*************************************************************************
 
  !Setting/Initializing Variables
  ntime = hist%mxhist
  natom_sc = size(hist%xred,2)
  factor   = 1._dp/natom_sc
  nterm =eff_pot%anharmonics_terms%ncoeff
  if(present(print_anh)) need_print_anh = print_anh
  ABI_ALLOCATE(symbols,(eff_pot%crystal%natom))
  ABI_ALLOCATE(terms,(nterm))
  call symbols_crystal(eff_pot%crystal%natom,eff_pot%crystal%ntypat,eff_pot%crystal%npsp,&
 &                     symbols,eff_pot%crystal%typat,eff_pot%crystal%znucl)
 
 
 !if the number of atoms in reference supercell into effpot is not correct,
 !wrt to the number of atom in the hist, we set map the hist and set the good supercell
  if (natom_sc /= eff_pot%supercell%natom) then
    call effective_potential_file_mapHistToRef(eff_pot,hist,comm,verbose=.TRUE.)
  end if

 !Check if input of order is correct
 ! TODO write error message here  
  if (any(mod(order,2) /= 0)) return 

 !we get the size of the supercell in the hist file
  do ii=1,3
    sc_size(ii) = int(anint(sqrt(eff_pot%supercell%rprimd(ii,1)**2+&
 &                               eff_pot%supercell%rprimd(ii,2)**2+&
 &                               eff_pot%supercell%rprimd(ii,3)**2) / &
 &                          sqrt(eff_pot%crystal%rprimd(ii,1)**2+&
 &                               eff_pot%crystal%rprimd(ii,2)**2+&
 &                               eff_pot%crystal%rprimd(ii,3)**2)))
  end do


 !Before the fit, compute constants with fit_data_compute.
 !Conpute the strain of each configuration.
 !Compute the displacmeent of each configuration.
 !Compute the variation of the displacement due to strain of each configuration.
 !Compute fixed forces and stresse and get the standard deviation.
 !Compute Sheppard and al Factors  \Omega^{2} see J.Chem Phys 136, 074103 (2012) [[cite:Sheppard2012]].
 !call fit_data_compute(fit_data,eff_pot,hist,comm,verbose=.FALSE.)
  call fit_data_compute(fit_data,eff_pot,hist,comm,verbose=.FALSE.)

  if(need_print_anh) call effective_potential_writeAnhHead(eff_pot%anharmonics_terms%ncoeff,&
&                            fn_bf,eff_pot%anharmonics_terms)                  

 !Before adding bound coefficients calculate MSD of initial model  
 !MS FOR THE MOMENT PRINT NO FILE 
  call fit_polynomial_coeff_computeMSD(eff_pot,hist,mse_ini,msef_ini,mses_ini,&
 &                                     natom_sc,ntime,fit_data%training_set%sqomega,&
 &                                     compute_anharmonic=.TRUE.,print_file=.FALSE.)


 !  Print the standard devition of initial model 
      write(message,'(6a,ES24.16,6a,ES24.16,2a,ES24.16,2a,ES24.16,a)' )ch10,&
 &                    ' Mean Standard Deviation values of the effective-potential',ch10,&
 &                    ' with respect to the training-set before attempted bounding (meV/atm):',&
 &               ch10,'   Energy          : ',&
 &               mse_ini*Ha_EV*1000*factor ,ch10,&
 &                    ' Goal function values of the effective.potential',ch10,& 
 &                    ' with respect to the test-set (eV^2/A^2):',ch10,&
 &                    '   Forces+Stresses : ',&
 &               (msef_ini+mses_ini)*(HaBohr_meVAng)**2,ch10,&
 &                    '   Forces          : ',&
 &               msef_ini*(HaBohr_meVAng)**2,ch10,&
 &                    '   Stresses        : ',&
 &               mses_ini*(HaBohr_meVAng)**2,ch10
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')

 ! Loop over orders givin in input in steps of two !TODO CHECK IF STEPS ARE WRITTEN LIKE THIS IN FORTRAN      
 !For the moment order loop commented
 !do iorder=order(1),order(2),2 
   !Loop over all original terms
   do iterm =1,nterm
     terms(iterm) = iterm
     ! Get current number of terms
     ! Add one for new bounding term 
     nterm_inloop = eff_pot%anharmonics_terms%ncoeff + 1 
     ABI_DATATYPE_ALLOCATE(my_coeffs,(nterm_inloop)) 
     !Copy original terms to my_coeffs 
     my_coeffs(1:nterm_inloop-1) =  eff_pot%anharmonics_terms%coefficients(:)
     !Message 
      write(message,'(2a,I3,a,I3,3a)' )ch10,&
&     ' Check term (',iterm,'/',nterm,'): ', trim(my_coeffs(iterm)%name),ch10 
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')
     ! Copy current term to last in array change the value of its coefficient to a start value
     ! The start is estimed to not be larger then half the initial term's value
     ! This is because hihger order terms should have smaller coefficients 
     my_coeffs(nterm_inloop) = eff_pot%anharmonics_terms%coefficients(iterm)  
     coeff_ini = my_coeffs(nterm_inloop)%coefficient / 2
     my_coeffs(nterm_inloop)%coefficient = coeff_ini 
     !Get number of displacements and equivalent terms for this term
     !Chose term one to get ndisp. ndisp is equal for all terms of the term 
     ndisp = eff_pot%anharmonics_terms%coefficients(iterm)%terms(1)%ndisp  
     nterm_of_term = my_coeffs(iterm)%nterm
     ! If the term is even and its coefficient positive we skip it. Also here we take terms(1) as example for all terms of term
     if(my_coeffs(iterm)%coefficient > 0 .and. .not. any(mod(my_coeffs(iterm)%terms(1)%power_disp(:),2) /= 0))then 
             ABI_DATATYPE_DEALLOCATE(my_coeffs)
              ! Message to Output 
              write(message,'(3a)' )ch10,&
&             ' ==>No need for high order bounding term',ch10
              call wrtout(ab_out,message,'COLL')
              call wrtout(std_out,message,'COLL')
             cycle 
     end if 
     ! Change all weights of equivalent terms of this term to 1 
     !(Even terms have allways a weight of 1...)
     ! Change all orders of the term to an even order
     ! TODO The sum of the orders may not be greater than the order specified in the input
     ! TODO generate more flexible routine that does stuff like ^4*^2 and ^2*^4 etc... 
     ! TODO compare generate term with terms already in potential, if already contained skip!
     do iterm_of_term=1,nterm_of_term 
       my_coeffs(nterm_inloop)%terms(iterm_of_term)%weight = 1
       do idisp=1,ndisp
         if(mod(my_coeffs(iterm)%terms(iterm_of_term)%power_disp(idisp),2) == 1)then 
           my_coeffs(nterm_inloop)%terms(iterm_of_term)%power_disp(idisp) = my_coeffs(iterm)%terms(iterm_of_term)%power_disp(idisp) + 1 
          else
           my_coeffs(nterm_inloop)%terms(iterm_of_term)%power_disp(idisp) = my_coeffs(iterm)%terms(iterm_of_term)%power_disp(idisp) + 2 
         endif
       enddo !idisp
     enddo !iterm_of_term 

     ! Get new name of term and set new terms to potential 
     call polynomial_coeff_getName(name,my_coeffs(nterm_inloop),symbols,recompute=.TRUE.)
     call polynomial_coeff_SetName(name,my_coeffs(nterm_inloop))
     call effective_potential_setCoeffs(my_coeffs,eff_pot,nterm_inloop)
        
     ! Message to Output 
      write(message,'(5a)' )ch10,&
&     ' ==> high order term: ', trim(my_coeffs(nterm_inloop)%name),' created',ch10
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')

     ABI_DATATYPE_DEALLOCATE(my_coeffs)  
     

     call fit_polynomial_coeff_computeMSD(eff_pot,hist,mse,msef,mses,&
 &                                 natom_sc,ntime,fit_data%training_set%sqomega,&
 &                                 compute_anharmonic=.TRUE.,print_file=.FALSE.)
     i = 0 
     do while(msef/msef_ini>=2)
        i = i + 1 
        write(*,*) i
        eff_pot%anharmonics_terms%coefficients(nterm_inloop)%coefficient = coeff_ini / 2**i
        call fit_polynomial_coeff_computeMSD(eff_pot,hist,mse,msef,mses,&
 &                                     natom_sc,ntime,fit_data%training_set%sqomega,&
 &                                     compute_anharmonic=.TRUE.,print_file=.FALSE.)
     enddo ! while mse/mse_ini>10  
   enddo !iterm
 !enddo ! order 

 
 write(message,'(3a)' )ch10,&
&     'Finished creating high-order terms',ch10
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL')

  
  call opt_effpot(eff_pot,nterm,terms,hist,comm,print_anh=.FALSE.)
  !DEALLOCATION 
  ABI_DEALLOCATE(terms)
  ABI_DEALLOCATE(symbols)
  call fit_data_free(fit_data)
end subroutine opt_effpotbound

end module m_opt_effpot

      
