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
 use m_polynomial_term
 use m_crystal,only : symbols_crystal

 implicit none

 public :: opt_effpot 
 public :: opt_effpotbound
 public :: opt_getHOforterm
 public :: opt_getCombisforterm
 public :: opt_getHoTerms
 public :: opt_getHOstrain
 public :: opt_getHOcrossdisp
 public :: opt_filterdisp
 public :: opt_getSingleDispTerms
 public :: opt_getHOSingleDispTerms
 private :: opt_boundcoeff
 private :: check_to_skip
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
 type(effective_potential_type),intent(inout) :: eff_pot
 type(abihist),intent(inout) :: hist
!arrays 
 integer,intent(in) :: opt_coeff(opt_ncoeff)
!Logicals
 logical,optional,intent(in) :: print_anh 
!Strings 
!Local variables ------------------------------
!scalars
 integer :: ii, info,natom_sc,ntime,unit_anh1,unit_anh2
 integer :: master,nproc,my_rank 
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
 logical :: need_print_anh,file_opened,iam_master
!Strings 
 character(len=1000) :: message
 character(len=1000) :: frmt
 character(len=fnlen) :: fn_bf='before_opt_diff', fn_af='after_opt_diff'
! *************************************************************************
 !MPI
 master = 0
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)


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
 &                                     natom_sc,ntime,fit_data%training_set%sqomega,comm,&
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
 &                                     natom_sc,ntime,fit_data%training_set%sqomega,comm,&
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
     &                                     natom_sc,ntime,fit_data%training_set%sqomega,comm,&
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
!!***

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

subroutine opt_effpotbound(eff_pot,order_ran,hist,comm,print_anh)

 implicit none 
         
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 type(effective_potential_type),target,intent(inout) :: eff_pot
 type(abihist),intent(inout) :: hist
!arrays 
 integer,intent(in) :: order_ran(2)
!Logicals
 logical,optional,intent(in) :: print_anh 
!Strings 
!Local variables ------------------------------
!scalars
 integer :: i,ii,natom_sc,ntime,iterm,nterm,j
 integer :: jterm, ncombi,ncombi1,ncombi2
 integer :: icombi
 integer :: nterm_start,nterm2
 integer :: nproc,my_rank,master
 !1406
 real(dp) :: factor,mse_ini,msef_ini,mses_ini,mse,msef,mses,coeff_ini=0.1
 real(dp) :: coeff_tmp
 real(dp),parameter :: HaBohr_meVAng = 27.21138386 / 0.529177249
!arrays 
 integer :: sc_size(3)
 integer,allocatable :: terms(:)
 logical,allocatable :: exists(:) 
 type(fit_data_type) :: fit_data
 real(dp) :: msefs_arr(2),coeff_opt(2)
 !real(dp), allocatable :: energy_coeffs(:,:),fcart_coeffs(:,:,:,:)
 !real(dp), allocatable :: strten_coeffs(:,:,:)
 !1406 strain_temrs_tmp
 type(polynomial_coeff_type),target,allocatable :: my_coeffs(:),my_coeffs_tmp(:)
 type(polynomial_coeff_type),allocatable :: singledisp_terms(:),HOsingledisp_terms(:)
 type(polynomial_coeff_type),allocatable :: HOcrossdisp_terms(:)
!Logicals
 logical :: need_print_anh=.FALSE. ! MARCUS FOR THE MOMENT PRINT NO FILES
 logical :: to_skip,iam_master
!Strings
 character(len=5),allocatable :: symbols(:)
 character(len=200):: name
 character(len=1000) :: message
 character(len=fnlen) :: fn_bf='before_opt_diff'!, fn_af='after_opt_diff'
!*************************************************************************
   !MPI variables
   master = 0
   nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
   iam_master = (my_rank == master)
  
  ! Say hello to the world!
  write(message, '(3a)' )'-Start Bound optimization of Anharmonic Potential ',ch10
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

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
  if (any(mod(order_ran,2) /= 0)) return 

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
 &                                     natom_sc,ntime,fit_data%training_set%sqomega,comm,&
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

 !MS DEV 
 call opt_getSingleDispTerms(singledisp_terms,eff_pot%crystal,sc_size,comm) 

 !For the moment order loop commented
 !do iorder=order(1),order(2),2, Order will be done per term 
 !Loop over all original terms + 1 
 ! + 1 to bound pure strain
  do iterm =1,nterm  !+1 
     if(iterm <=nterm)then
       ncombi1=0 
       ncombi2=0
       !Store for optimization
       terms(iterm) = iterm
       !Message: The world wants to know where we stand Batman
       write(message, '(a,(80a),a)' ) ch10,&
&      ('_',ii=1,80),ch10
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
        write(message,'(2a,I3,a,I3,3a)' )ch10,&
&       ' Check term (',iterm,'/',nterm,'): ', trim(eff_pot%anharmonics_terms%coefficients(iterm)%name),ch10 
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')                                     
       to_skip = .FALSE.
       to_skip = check_to_skip(eff_pot%anharmonics_terms%coefficients(iterm))
       !Skip term if it doesn't need bounding 
       if(.not. to_skip)then
          !Get List of high order single Terms for terms 
          call opt_getHOSingleDispTerms(eff_pot%anharmonics_terms%coefficients(iterm),&
&                                      HOsingledisp_terms,symbols,singledisp_terms,order_ran,ncombi1)
          !Get List of high order cross Terms for term if ndisp > 1
          if(eff_pot%anharmonics_terms%coefficients(iterm)%terms(1)%ndisp>1 .or. &
&            eff_pot%anharmonics_terms%coefficients(iterm)%terms(1)%ndisp /= 0 .and. &
&            eff_pot%anharmonics_terms%coefficients(iterm)%terms(1)%nstrain /= 0)then 
             call opt_getHOcrossdisp(HOcrossdisp_terms,ncombi2,eff_pot%anharmonics_terms%coefficients(iterm),order_ran)
          endif  
          !Copy everything together
          ABI_DATATYPE_ALLOCATE(my_coeffs,(size(eff_pot%anharmonics_terms%coefficients)+size(HOsingledisp_terms)))
          !my_coeffs = eff_pot%anharmonics_terms%coefficients + HOsingledisp_terms
          !Test 
          do i=1,size(eff_pot%anharmonics_terms%coefficients)+size(HOsingledisp_terms)
             if(i<=size(eff_pot%anharmonics_terms%coefficients))then 
                call polynomial_coeff_init(coeff_ini,eff_pot%anharmonics_terms%coefficients(i)%nterm,my_coeffs(i),&
& eff_pot%anharmonics_terms%coefficients(i)%terms,eff_pot%anharmonics_terms%coefficients(i)%name,check=.TRUE.)
            else 
                j=i-size(eff_pot%anharmonics_terms%coefficients)
                call polynomial_coeff_init(coeff_ini,HOsingledisp_terms(j)%nterm,my_coeffs(i),HOsingledisp_terms(j)%terms,&
                                           HOsingledisp_terms(j)%name,check=.TRUE.)
             endif 
          enddo 

          if(ncombi2 > 0)then
            ABI_DATATYPE_ALLOCATE(my_coeffs_tmp,(size(my_coeffs)))
            my_coeffs_tmp = my_coeffs
            if(allocated(my_coeffs)) call polynomial_coeff_list_free(my_coeffs)
            ABI_DATATYPE_ALLOCATE(my_coeffs,(size(my_coeffs_tmp)+size(HOcrossdisp_terms)))
            !my_coeffs = my_coeffs_tmp + HOcrossdisp_terms
            do i=1,size(my_coeffs_tmp)+size(HOcrossdisp_terms)
               if(i<=size(my_coeffs_tmp))then 
                  call polynomial_coeff_init(coeff_ini,my_coeffs_tmp(i)%nterm,my_coeffs(i),my_coeffs_tmp(i)%terms,&
                                             my_coeffs_tmp(i)%name,check=.TRUE.)
              else 
                  j=i-size(my_coeffs_tmp)
                  call polynomial_coeff_init(coeff_ini,HOcrossdisp_terms(j)%nterm,my_coeffs(i),HOcrossdisp_terms(j)%terms,&
                                             HOcrossdisp_terms(j)%name,check=.TRUE.)
               endif 
            enddo 
            if(allocated(my_coeffs_tmp)) call polynomial_coeff_list_free(my_coeffs_tmp)
          endif
       else 
         ncombi2=0
         ncombi1=0
       endif 
       ncombi = ncombi1 + ncombi2 
       nterm_start = eff_pot%anharmonics_terms%ncoeff
     else ! if iterm = nterm + 1 => Take care about strain 
       call opt_getHOstrain(my_coeffs,ncombi,nterm_start,eff_pot,order_ran,comm)     
     endif ! 

     do icombi=1,ncombi
            ! Copy all the terms in eff pot 
            ! Get new name of term and set new terms to potential 
            !write(*,*) 'ndisp of term', my_coeffs(nterm_start+icombi)%nterm
            !write(*,*) 'and wath is nterm_start', nterm_start,'and icomb btw', icombi
            call polynomial_coeff_getName(name,my_coeffs(nterm_start+icombi),symbols,recompute=.TRUE.)
            call polynomial_coeff_SetName(name,my_coeffs(nterm_start+icombi))
          
            ! Set dimensions of temporary my_coeffs array 
            nterm2 = eff_pot%anharmonics_terms%ncoeff + 1
            ABI_DATATYPE_ALLOCATE(my_coeffs_tmp,(nterm2))
            ! Copy terms of previous cycle
            my_coeffs_tmp(1:nterm2-1) = eff_pot%anharmonics_terms%coefficients 
            !Put new term to my_coeffs_tmp
            !my_coeffs_tmp(nterm2) = my_coeffs(nterm_start+icombi)
            call polynomial_coeff_init(my_coeffs(nterm_start+icombi)%coefficient,my_coeffs(nterm_start+icombi)%nterm,&
&                                      my_coeffs_tmp(nterm2),my_coeffs(nterm_start+icombi)%terms,my_coeffs(nterm_start+icombi)%name)
            ! If order is greater than specified cycle
            if(sum(my_coeffs_tmp(nterm2)%terms(1)%power_disp) & 
&             +sum(my_coeffs_tmp(nterm2)%terms(1)%power_strain) > maxval(order_ran))then 
               call polynomial_coeff_list_free(my_coeffs_tmp)
               cycle
            endif 
            ! Message to Output 
             write(message,'(5a)' )ch10,&
&            ' ==> high order term: ', trim(my_coeffs_tmp(nterm2)%name),' created',ch10
             call wrtout(ab_out,message,'COLL')
             call wrtout(std_out,message,'COLL')
            ! Check if generated term is not already contained in effpot
            ! If yes cycle 
            ABI_ALLOCATE(exists,(nterm2))
            exists=.FALSE.
            do jterm=1,nterm2-1
!               write(std_out,*) 'what is jterm here', jterm 
!               write(std_out,*) "c1%nterm", my_coeffs_tmp(jterm), "c2%nterm", my_coeffs_tmp(nterm2)%nterm
               exists(jterm) = coeffs_compare(my_coeffs_tmp(jterm),my_coeffs_tmp(nterm2))
            enddo !jterm
            if(any(exists))then 
               write(message,'(3a)' )ch10,&
&              '   ==> Term exists already. We cycle',ch10
               call wrtout(ab_out,message,'COLL')
               call wrtout(std_out,message,'COLL')
               call polynomial_coeff_list_free(my_coeffs_tmp)
               ABI_DEALLOCATE(exists)
               cycle 
            endif 
          
            ! Set new term into effective potential 
            call effective_potential_setCoeffs(my_coeffs_tmp,eff_pot,nterm2)
            
            ! Tell the world what we do, They want to know.     
            write(message,'(3a)' )ch10,&
&           '   ==> Optimizing coefficient',ch10
            call wrtout(ab_out,message,'COLL')
            call wrtout(std_out,message,'COLL')
             
            ! Deallocation in loop 
            call polynomial_coeff_list_free(my_coeffs_tmp)
            !Optimizing coefficient old style 

            if(iterm>nterm)then
               ! MS 2006 Decomment for old style optimization  
               call fit_polynomial_coeff_computeMSD(eff_pot,hist,mse,msef,mses,&
 &                                           natom_sc,ntime,fit_data%training_set%sqomega,comm,&
 &                                           compute_anharmonic=.TRUE.,print_file=.FALSE.)
               i = 0
               write(message,'(a,I2,a,F12.7)') "cycle ", i ," (msef+mses)/(msef_ini+mses_ini): ", (msef+mses)/(msef_ini+mses_ini)
               call wrtout(std_out,message,'COLL')
               write(message,'(a,I2,a,ES24.16)') "cycle ", i ," (msef+mses): ", (msef+mses)
               call wrtout(std_out,message,'COLL')
               do  while((msef+mses)/(msef_ini+mses_ini) >= 1.001)
                  i = i + 1 
                  eff_pot%anharmonics_terms%coefficients(nterm2)%coefficient = coeff_ini / 2**i
                  call fit_polynomial_coeff_computeMSD(eff_pot,hist,mse,msef,mses,&
 &                                              natom_sc,ntime,fit_data%training_set%sqomega,comm,&
 &                                              compute_anharmonic=.TRUE.,print_file=.FALSE.)
 
                  write(message,'(a,I2,a,F12.7)') "cycle ", i ," (msef+mses)/(msef_ini+mses_ini): ", (msef+mses)/(msef_ini+mses_ini)
                  call wrtout(std_out,message,'COLL')
                  write(message,'(a,I2,a,ES24.16)') "cycle ", i ," (msef+mses): ", (msef+mses)
                  call wrtout(std_out,message,'COLL')
               enddo ! while mse/mse_ini>1.0001 
               write(message,'(a,F12.7)') "coeff after opt:",   eff_pot%anharmonics_terms%coefficients(nterm2)%coefficient 
               call wrtout(std_out,message,'COLL')
               msef_ini = msef
               mses_ini = mses
!           !Optimize coefficient with opt routine 
!           optterm(1)= nterm2 
!           nterm_opt = 1
!           call opt_effpot(eff_pot,nterm_opt,optterm,hist,comm,print_anh=.FALSE.)
!           write(*,*) "coeff after opt:",   eff_pot%anharmonics_terms%coefficients(nterm2)%coefficient 
!            !Store new "inital precision for next coefficient
!            msef_ini = msef
!            mses_ini = mses

            else 
            !Optimizing coefficient precisely ?
              coeff_opt = 0 
              msefs_arr = 0 
                i = 1 
                do while(i<=2)
                  eff_pot%anharmonics_terms%coefficients(nterm2)%coefficient = &
&                 eff_pot%anharmonics_terms%coefficients(nterm2)%coefficient/ 2**(i-1)
                  call fit_polynomial_coeff_computeMSD(eff_pot,hist,mse,msef,mses,&
&                                              natom_sc,ntime,fit_data%training_set%sqomega,comm,&
&                                              compute_anharmonic=.TRUE.,print_file=.FALSE.)
 
                  write(message,'(a,I2,a,ES24.16)') "cycle ",i," (msef+mses)/(msef_ini+mses_ini): ",(msef+mses)/(msef_ini+mses_ini)
                  call wrtout(std_out,message,'COLL')
                  write(message,'(a,I2,a,ES24.16)') "cycle ", i ," (msef+mses): ", (msef+mses)
                  call wrtout(std_out,message,'COLL')
                  coeff_opt(i) =  eff_pot%anharmonics_terms%coefficients(nterm2)%coefficient
                  msefs_arr(i) =  (msef+mses)/(msef_ini+mses_ini)
                  if(i==2 .and. abs(msefs_arr(1)-msefs_arr(2)) < tol8)then 
                     eff_pot%anharmonics_terms%coefficients(nterm2)%coefficient =& 
                     eff_pot%anharmonics_terms%coefficients(nterm2)%coefficient*10d5 
                     write(message,'(5a)') ch10,"Differences between test-cycles to small increase",ch10, & 
&                                          "test coefficient value by factor 1000",ch10
                     call wrtout(std_out,message,'COLL')
                     i = 1
                  else 
                     i=i+1
                  end if 
                enddo ! while mse/mse_ini>10
                eff_pot%anharmonics_terms%coefficients(nterm2)%coefficient = opt_boundcoeff(msefs_arr,coeff_opt)
                write(message,'(a,ES24.16)') "coeff after opt1:",   eff_pot%anharmonics_terms%coefficients(nterm2)%coefficient 
                call wrtout(std_out,message,'COLL')
                coeff_tmp = ANINT(eff_pot%anharmonics_terms%coefficients(nterm2)%coefficient*10d10)
                eff_pot%anharmonics_terms%coefficients(nterm2)%coefficient = coeff_tmp/10d10
                write(message,'(a,ES24.16)') "coeff after opt2:",   eff_pot%anharmonics_terms%coefficients(nterm2)%coefficient 
                call wrtout(std_out,message,'COLL')
                call fit_polynomial_coeff_computeMSD(eff_pot,hist,mse,msef,mses,&
 &                                               natom_sc,ntime,fit_data%training_set%sqomega,comm,&
 &                                               compute_anharmonic=.TRUE.,print_file=.FALSE.)
                write(message,'(a,ES24.16)') "(msef+mses)/(msef_ini+mses_ini) after_opt: ", (msef+mses)/(msef_ini+mses_ini)
                call wrtout(std_out,message,'COLL')
                msef_ini = msef
                mses_ini = mses
           endif
           !DEALLOCATION 
           ABI_DEALLOCATE(exists)
     enddo ! icombi

       if(allocated(my_coeffs)) call polynomial_coeff_list_free(my_coeffs)
       if(allocated(HOcrossdisp_terms)) call polynomial_coeff_list_free(HOcrossdisp_terms)  
       if(allocated(HOsingledisp_terms)) call polynomial_coeff_list_free(HOsingledisp_terms)  
   enddo !iterm
 !enddo ! order 

 if(allocated(singledisp_terms)) call polynomial_coeff_list_free(singledisp_terms)  
 write(message, '(a,(80a),a)' ) ch10,&
&('_',ii=1,80),ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 
 write(message,'(3a)' )ch10,&
&     ' Finished creating high-order terms',ch10  !,&
!&     ' Optimize initial anharmonic terms !NOT IS COMMENTED NOW!',ch10
      call wrtout(ab_out,message,'COLL')
      call wrtout(std_out,message,'COLL') 
 !  call opt_effpot(eff_pot,nterm,terms,hist,comm,print_anh=.FALSE.)

!  Print the standard devition of final model 
      write(message,'(6a,ES24.16,6a,ES24.16,2a,ES24.16,2a,ES24.16,a)' )ch10,&
 &                    ' Mean Standard Deviation values of the effective-potential',ch10,&
 &                    ' with respect to the training-set after attempted bounding (meV/atm):',&
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


  !DEALLOCATION 
  ABI_DEALLOCATE(symbols)
  ABI_DEALLOCATE(terms)
  
  !ABI_DATATYPE_DEALLOCATE(my_coeffs)  
  call fit_data_free(fit_data)
end subroutine opt_effpotbound
!!***


!!****f* m_opt_effpot/opt_getHOforterm
!!
!! NAME
!! opt_effpotbound
!!
!! FUNCTION
!! Compute possible high orders for a given anharmonic term. 
!! In the range of order_start,order_stop 
!!
!!
!! INPUTS
!! term<polynomial_coeff_type>:anharmonic term 
!! order_range(2):start and stop order desired by user
!! 
!! 
!! OUTPUT
!! order_start:possible start order 
!! order_stop: possible stop order 
!!
!! PARENTS
!! opt_effpotbound
!!
!! CHILDREN
!!
!! SOURCE

subroutine opt_getHOforterm(term,order_range,order_start,order_stop)

 implicit none 
         
!Arguments ------------------------------------
!scalars
 type(polynomial_coeff_type),intent(in) :: term
!arrays 
 integer,intent(in) :: order_range(2)
 integer,intent(out) :: order_start, order_stop 
!Logicals
!Strings 
!Local variables ------------------------------
!scalars
 integer :: idisp,ndisp,nterm_of_term,power_tot 
!arrays 
 integer,allocatable :: powers(:) 
!Logicals
!Strings
!*************************************************************************

     !Get/Initialize variables
     ndisp = term%terms(1)%ndisp
     nterm_of_term = term%nterm

     ABI_ALLOCATE(powers,(ndisp)) 
     powers = term%terms(1)%power_disp

     power_tot = 0 


     !Get rid off odd displacements
     do idisp=1,ndisp
       if(mod(term%terms(1)%power_disp(idisp),2) == 1)then 
          powers(idisp) = powers(idisp) + 1
       end if 
     enddo !idisp
     ! Count order
     do idisp=1,ndisp
        power_tot = power_tot + powers(idisp)
     enddo 

     ! Get start and stop order for this term 
     ! If term doesn't fit in order range give back order_start = ordre_stop = 0 
     if(power_tot >= order_range(1) .and. power_tot <=order_range(2))then
        order_start = power_tot 
        order_stop  = order_range(2) 
     elseif(power_tot < order_range(1))then 
        order_start = order_range(1)
        order_stop  = order_range(2)   
     elseif(power_tot > order_range(2))then          
        order_start = 0 
        order_stop  = 0 
     endif 

     ABI_DEALLOCATE(powers) 

end subroutine opt_getHOforterm
!!***


!!****f* m_opt_effpot/opt_getCombisforterm
!!
!! NAME
!! opt_getCombisforterm
!!
!! FUNCTION
!! For a given order range: order_start, order_stop 
!! calculate number of total possible combinations ncombi
!! and calculat combinations per order ncombi_order(i)
!!
!!
!! INPUTS
!! order_start: start order for bounding terms  
!! order_end: end order for bounding terms 
!! ndisp: number of displacements a given term contains 
!! 
!! 
!! OUTPUT
!! ncombi: total number of combinations 
!! ncombi_order: array with number of combinations per order
!! 
!!
!! PARENTS
!! opt_effpotbound
!!
!! CHILDREN
!!
!! SOURCE

subroutine opt_getCombisforterm(order_start,order_end,ndisp,ncombi,ncombi_order)

 implicit none 
         
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndisp
 integer,intent(in) :: order_start, order_end
!arrays 
 integer,intent(out) :: ncombi
 integer,intent(out) :: ncombi_order(:)
!Logicals
!Strings 
!Local variables ------------------------------
!scalars
 integer :: i
 integer :: order,iorder1,iorder2
!arrays 
!integer 
!Logicals
!Strings
 character(len=1000) :: message
!*************************************************************************

!Test 
if(mod(order_start,2) /= 0 .or. mod(order_end,2) /= 0)then 
   ! Message to Output 
   write(message,'(4a)' )ch10,&
&  'Either start or stop order are not even numbers',ch10,&
&  'Action: change bound_range in input',ch10 
   MSG_ERROR(message)
endif 

!Initialize Variables 
i = 0
ncombi = 0
ncombi_order = 0  

!Calculate Combinations 
do order=order_start,order_end,2
   i = i+1
   if(ndisp == 1)then
      ncombi = ncombi + 1 
      ncombi_order(i) = 1
   else 
      do iorder1 = 2,order-2*(ndisp-1),2 
         if(ndisp*iorder1 == order)then
            ncombi = ncombi + 1 
            ncombi_order(i) = ncombi_order(i) + 1
            cycle
         endif
         do iorder2=iorder1+2,order-2*(ndisp-1),2    
           if( iorder1 + (ndisp-1)*iorder2 == order)then 
             ncombi = ncombi + ndisp
             ncombi_order(i) = ncombi_order(i) + ndisp
           elseif(iorder1 * (ndisp-1) + iorder2 == order)then 
             ncombi = ncombi + ndisp 
             ncombi_order(i) = ncombi_order(i) + ndisp
           endif
         enddo !iorder2 
      enddo !iorder1 !
   endif 
   !write(*,*) 'ncombi(',i,') for order',order,' is:', ncombi_order(i), 'are we happy?'
enddo !order

!write(*,*) ncombi_order(:)
!write(*,*) 'ncombi for term is:', ncombi, 'are we happy?'

end subroutine opt_getCombisforterm
!!***

!!****f* m_opt_effpot/opt_getHoTerms
!!
!! NAME
!! opt_geHoTerms
!!
!! FUNCTION
!! For a term give all possible high order terms 
!! Attention order_start, order_stop, ncombi and 
!! ncombi_order have to be calculated before! 
!!
!!
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

subroutine opt_getHoTerms(terms,order_start,order_stop,ndisp,ncombi_order)

 implicit none 
         
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndisp
 integer,intent(in) :: order_start, order_stop
!arrays 
 integer,intent(in) :: ncombi_order(:)
!Logicals
 type(polynomial_coeff_type),intent(inout) :: terms(:)
!Strings 
!Local variables ------------------------------
!scalars
 integer :: i,icombi,icombi2,icombi_start,icombi_stop,idisp,nterm_of_term
 integer :: order,iterm_of_term,jdisp,power_tot
 integer :: jdisp1,jdisp2,sec
 real(sp) :: to_divide,divider1,divider2,divided
!arrays 
!integer 
!Logicals
 logical :: equal_term_done
!Strings
 character(len=1000) :: message
!*************************************************************************
     !Get Variables 
     nterm_of_term = terms(1)%nterm

     ! Create all possible combinaions for the specified orders
     ! If anybody has ever, ever to read and understand this I'm terribly sorry 
     ! ---this will be quite hell---
     power_tot = 0
     icombi_start = 1
     i = 0
     !write(*,*) "what is ncombi?", ncombi
     !write(*,*) "what is ncmobi_order", ncombi_order 
     !write(*,*) "order_start", order_start
     !write(*,*) "order_stop", order_stop
     order = order_start
     ! TODO work here icombi start and order conting does not work yet. 
     do order=order_start,order_stop,2
        !write(*,*) 'Was I now here?!(order-loop)'
        i = i + 1
        icombi_stop = icombi_start + ncombi_order(i) - 1
        equal_term_done = .FALSE.
        !write(*,*) 'order', order
        !write(*,*) 'icombi_start', icombi_start
        !write(*,*) 'icombi_stop', icombi_stop
        icombi=icombi_start 
        do while (icombi<=icombi_stop)
           power_tot = 0 
           do idisp=1,ndisp
              !write(*,*) "what is icombi here?", icombi 
              !write(*,*) "what is icombi_stop here?", icombi_stop
              power_tot = power_tot + terms(icombi)%terms(1)%power_disp(idisp)
           enddo
           ! Probably have to increase order at same time 
           ! If the term already has the right order we cycle 
           ! we increase icombi_start and icombi to go to the next term in the array
           if(power_tot == order)then
              icombi_start = icombi_start + 1
              icombi = icombi + 1
              !write(*,*) 'icombi-start in if', icombi_start 
              cycle
           ! If the term is not already in the right order, we manipulate it until 
           ! it is 
           else        
              jdisp1 = 1 
              sec = 1
              !write(*,*) 'what is ndisp actually', ndisp
              ! Treat single permutations from the bottom of the order 
              ! so start from ^2^2^2 and get ^6^2^2,^2^6^2, and ^2^2^6 f.E.
              ! do loop over displacements 
              do while(jdisp1<=ndisp .and. sec < 100)                    
                 !write(*,*) "what is jdisp1 here?", jdisp1
                 !write(*,*) "what is icombi here?", icombi
                 !write(*,*) "what is icombi_sotp?", icombi_stop
                 ! Increase the order of the displacements 
                 do iterm_of_term=1,nterm_of_term
                    terms(icombi)%terms(iterm_of_term)%power_disp(jdisp1) =&
&                            terms(icombi)%terms(iterm_of_term)%power_disp(jdisp1) + 2 
                 enddo
                 ! Check total order of term after increse
                 power_tot=0
                 do jdisp2=1,ndisp 
                    power_tot = power_tot + terms(icombi)%terms(1)%power_disp(jdisp2)
                 enddo
                 ! If the term is at the right order do next displacement 
                 ! increase icombi and icombi_start to go to next term in array 
                 if(power_tot == order)then 
                    icombi = icombi + 1
                    icombi_start = icombi_start +1
                    !write(*,*) 'what is icombi_start here', icombi_start
                    jdisp1 = jdisp1 + 1
                    !write(*,*) 'and what is jdisp1?', jdisp1

                 endif  
               enddo!jdisp1
               !Treat permutations in terms ndisp > 2 and order >= 10 
               !Start f.E. from ^4^4^4 and get ^4^2^4, ^2^4^4,^4^4^2
               if(icombi_stop - icombi > 1 .and. ndisp >2 )then    
                  do icombi2=icombi,icombi+ndisp-1
                      do jdisp=1,ndisp
                         do iterm_of_term=1,nterm_of_term
                            terms(icombi2)%terms(iterm_of_term)%power_disp(jdisp) = &
&                                    terms(icombi2)%terms(iterm_of_term)%power_disp(jdisp) +2
                            !write(*,*) "What's the power now?", terms(icombi2)%terms(iterm_of_term)%power_disp(jdisp)
                         enddo                     
                      enddo !jdisp
                  enddo 
                  jdisp1 = 1                  
                  do while(jdisp1<=ndisp .and. sec < 100) 
                      sec = sec + 1  
                      !write(*,*) "how often did I go here, hu ?"                  
                     !write(*,*) "I did at least one displacement" 
                     do iterm_of_term=1,nterm_of_term
                        terms(icombi)%terms(iterm_of_term)%power_disp(jdisp1) = &
&                                terms(icombi)%terms(iterm_of_term)%power_disp(jdisp1) - 2 
                     enddo
                     power_tot=0
                     do jdisp2=1,ndisp 
                        power_tot = power_tot + terms(icombi)%terms(1)%power_disp(jdisp2)
                     enddo
                     if(power_tot == order)then 
                        icombi = icombi + 1
                        icombi_start = icombi_start +1
                        !write(*,*) 'what is icombi_start here', icombi_start
                        jdisp1 = jdisp1 + 1
                        !write(*,*) 'and what is jdisp1?', jdisp1
                     endif  
                   enddo!jdisp1  
                   ! Message to Output 
                   if(sec>100)then 
                      write(message,'(4a)' )ch10,&
&                     "You're stuck in a while loop.",ch10,&
&                     'Action: Contact Abinit Group',ch10 
                      MSG_ERROR(message)   
                   endif                      
               endif! (icombi_stop - icombi)  
               !write(*,*) 'I was here!'
               to_divide = real(order)
               divider1 = real(ndisp)
               divided = real(to_divide/divider1)
               divider2 = real(2)
               !write(*,*) 'divided', divided, 'divider2', divider2
               !Treat terms with even power f.E. ^2^2^2^2, ^4^4 etc...
               if(mod(divided,divider2) == 0 .and. .not. equal_term_done .and. ndisp > 1)then 
                  !write(*,*) "Sometimes I should be here sometimes I shouldn't" 
                  do jdisp=1,ndisp
                     do iterm_of_term=1,nterm_of_term
                        terms(icombi)%terms(iterm_of_term)%power_disp(jdisp) = order/ndisp 
                     enddo                     
                  enddo !jdisp
                  if(order < order_stop)then
                     do icombi2=icombi+1,icombi+ndisp
                         do jdisp=1,ndisp
                            do iterm_of_term=1,nterm_of_term
                               terms(icombi2)%terms(iterm_of_term)%power_disp(jdisp) = order/ndisp 
                            enddo                     
                         enddo !jdisp
                     enddo
                  endif 
                  equal_term_done = .TRUE.
                  icombi = icombi + 1
                  icombi_start = icombi_start +1
               endif ! equal term if 
           endif ! power_tot == order  
        enddo !icombination
      enddo !order

end subroutine opt_getHoTerms
!!***

!!****f* m_opt_effpot/opt_filterdisp
!!
!! NAME
!! opt_opt_filterdisp
!!
!! FUNCTION
!! If a anharmonic term represents a strain-phonon coupling 
!! delete the strain and only keep the displacement part.
!!
!!
!! INPUTS
!! term<polynomial_coeff_type>: anharmonic term to check 
!! nterm_of_term: number of symmetry equivalent terms for term 
!! 
!! OUTPUT
!! term<polynomial_coeff_type>: only the displacement part of original term 
!!
!! PARENTS
!! opt_effpotbound
!!
!! CHILDREN
!! m_polynomial_coeff.F90/polynomial_coeff_init 
!!
!! SOURCE

subroutine opt_filterdisp(term,nterm_of_term)

 implicit none 
         
!Arguments ------------------------------------
!scalars
 type(polynomial_coeff_type),intent(inout) :: term
 integer,intent(in) :: nterm_of_term
 !arrays 
!Logicals
!Strings 
!Local variables ------------------------------
!scalars
 integer :: iterm_of_term
!reals 
 real(dp) :: coeff
!arrays 
 type(polynomial_term_type) :: terms(nterm_of_term)
!Logicals
!Strings
!*************************************************************************

!Initialize/Get Variables 
!nterm_of_term = term%nterm

! Set strain to zero in terms 
do iterm_of_term = 1, nterm_of_term 
  !Set strain in all terms to zero 
  !terms(iterm_of_term)%nstrain = 0
  !Free initial strain array 
  !ABI_DEALLOCATE(term%terms(iterm_of_term)%strain)
  !ABI_DEALLOCATE(term%terms(iterm_of_term)%power_strain)
  !Reallocate them with size zero 
enddo ! iterm_of term 

do iterm_of_term=1,nterm_of_term
   !terms(iterm_of_term) = term%terms(iterm_of_term)
   call polynomial_term_init(term%terms(iterm_of_term)%atindx,term%terms(iterm_of_term)%cell,term%terms(iterm_of_term)%direction,&
&                            term%terms(iterm_of_term)%ndisp,term%terms(iterm_of_term)%nstrain,terms(iterm_of_term),&
&                            term%terms(iterm_of_term)%power_disp,term%terms(iterm_of_term)%power_strain,&
&                            term%terms(iterm_of_term)%strain,term%terms(iterm_of_term)%weight,check=.TRUE.)
   terms(iterm_of_term)%nstrain = 0
   terms(iterm_of_term)%power_strain = 0
   terms(iterm_of_term)%strain = 0
enddo

call polynomial_coeff_free(term)
!Reinitial term 
!check=.TRUE. checks for duplicate terms
call polynomial_coeff_init(coeff,nterm_of_term,term,terms,check=.TRUE.)

do iterm_of_term=1,nterm_of_term
  call polynomial_term_free(terms(iterm_of_term))
enddo 
!if(nterm_of_term /= term%nterm)then 
!  write(*,*) "nterm_of_term changed after deleting strain"
!endif


end subroutine opt_filterdisp
!!***

!!****f* m_opt_effpot/opt_getHOstrain
!!
!! NAME
!! opt_getHOstrain
!!
!! FUNCTION
!! Get HO anharmnonic strain terms for bounding and add them to list 
!! of existing terms in effective potential
!!
!! INPUTS
!! eff_pot<effective_potential_type>: datatype with all the informations 
!! 				      about the effective potential 
!! power_strain(2): start and stop order for strain terms
!! comm: mpi communicator (at the moment only sequential tested) 
!! 
!! OUTPUT
!! terms<polynomial_coeff_type>: list with original terms in effective
!!                               potential + HO even strain terms
!!
!! PARENTS
!! opt_effpotbound
!!
!! CHILDREN
!! m_polynomial_coeff.F90/polynomial_coeff_init 
!!
!! SOURCE

subroutine opt_getHOstrain(terms,ncombi,nterm_start,eff_pot,power_strain,comm)

 implicit none 
         
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 type(polynomial_coeff_type),allocatable,intent(inout) :: terms(:)
 type(effective_potential_type), intent(in) :: eff_pot 
 integer,intent(in) :: power_strain(2) 
 integer,intent(out) :: ncombi,nterm_start
 !arrays 
!Logicals
!Strings 
!Local variables ------------------------------
!scalars
 integer ::  nterm_tot_tmp,icombi 
 integer :: ii 
!reals 
 type(crystal_t) :: crystal
!arrays
 type(polynomial_coeff_type),allocatable :: strain_terms_tmp(:)
!Logicals
!Strings
 character(len=1000) :: message
!*************************************************************************
    !Get variables 
    crystal = eff_pot%crystal

    write(message, '(a,(80a),a)' ) ch10,&
&   ('_',ii=1,80),ch10
    call wrtout(ab_out,message,'COLL')
    call wrtout(std_out,message,'COLL')
     write(message,'(3a)' )ch10,&
&    ' Chreate high order strain terms ',ch10 
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')

    !1406 get count of high order even anharmonic strain terms and the strain terms itself  
    call polynomial_coeff_getEvenAnhaStrain(strain_terms_tmp,crystal,ncombi,power_strain,comm)
    ! Allocate my_coeffs with ncombi free space to work with 

    nterm_start = eff_pot%anharmonics_terms%ncoeff
    nterm_tot_tmp = eff_pot%anharmonics_terms%ncoeff + ncombi
    ABI_DATATYPE_ALLOCATE(terms,(nterm_tot_tmp)) 
    do icombi=1,ncombi
       terms(nterm_start+icombi) = strain_terms_tmp(icombi)
       terms(nterm_start+icombi)%coefficient = 10000000      ! eff_pot%harmonics_terms%elastic_constants(1,1)
    enddo

end subroutine opt_getHOstrain
!!***

!!****f* m_opt_effpot/opt_getHOcrossdisp
!!
!! NAME
!! opt_getHOcrossdisp
!!
!! FUNCTION
!! Get even high order displacement terms for a given input term and 
!! add them to an excisting list of terms. If the term is strain phonon type 
!! the strain part gets deletet and the high order terms for the displacement 
!! part are computed. 
!! Example: for a tree linear term with three displacements x*y*z the even high order 
!!          possibilites are computed and stored.
!!          For range 6 to 8: x^2*y^2*z^2,x^4*y^2*z^2,x^2*y^4*z^2,x^2*y^2*z^4. 
!!
!! INPUTS
!! eff_pot<effective_potential_type>: datatype with all the informations 
!! 				      about the effective potential 
!! power_disp(2): start and stop order for disp terms
!! comm: mpi communicator (at the moment only sequential tested) 
!! 
!! OUTPUT
!! terms<polynomial_coeff_type>: list with original terms in effective
!!                               potential + HO even disp terms
!!
!! PARENTS
!! opt_effpotbound
!!
!! CHILDREN
!! m_polynomial_coeff.F90/polynomial_coeff_init 
!!
!! SOURCE

subroutine opt_getHOcrossdisp(terms_out,ncombi,term_in,power_disp)

 implicit none 
         
!Arguments ------------------------------------
!scalars
 type(polynomial_coeff_type),allocatable,intent(inout) :: terms_out(:)
 type(polynomial_coeff_type),intent(inout) :: term_in
 integer,intent(in) :: power_disp(2) 
 integer,intent(out) :: ncombi
 !arrays 
!Logicals
!Strings 
!Local variables ------------------------------
!scalars
 integer ::  ndisp,nterm_of_term
 integer ::  order_start,order_stop,norder
 integer ::  icombi,idisp,iterm_of_term
 integer :: ncombi_tot 
!reals
 real(dp) :: coeff_ini 
!arrays
 type(polynomial_coeff_type) :: term
 integer,allocatable :: ncombi_order(:),dummy(:)
!Logicals
 logical :: had_strain 
!Strings
 character(len=1000) :: message
!*************************************************************************
       !Get/Set Variables         
       norder = abs(((power_disp(2)-power_disp(1))/2)) + 1 
       ABI_ALLOCATE(ncombi_order,(norder))
       ncombi_order = 0 

       ncombi = 0 
       !Get this term (iterm) and infromations about it 
       !Get number of displacements and equivalent terms for this term
       !Chose term one to get ndisp. ndisp is equal for all terms of the term
       !Get minimum oder for this term
       !Get total number of terms in effpot for message
       ndisp = term_in%terms(1)%ndisp
       nterm_of_term = term_in%nterm
       ABI_ALLOCATE(dummy,(5))      
       call polynomial_coeff_init(coeff_ini,nterm_of_term,term,term_in%terms,term_in%name,check=.true.)
       ABI_DEALLOCATE(dummy)      
       ! Check if term has strain component. 
       ! If yes filter strain and fit high order atomic displacement terms
       had_strain = .FALSE.
       ABI_ALLOCATE(dummy,(5))      
       if(term%terms(1)%nstrain /= 0)then
                ! Message to Output 
                write(message,'(5a)' )ch10,&
&               '- Term has strain compenent',ch10,&
&               ' -> Filter Displacement',ch10
                call wrtout(ab_out,message,'COLL')
                call wrtout(std_out,message,'COLL')
                call opt_filterdisp(term,nterm_of_term)
                !Get new value of symmetry equivalent term nterm_of_term 
                nterm_of_term = term%nterm
                !Remember if this term had strain
                had_strain = .TRUE.
               !cycle 
       endif 
       ABI_DEALLOCATE(dummy)      
       ! Ok we want it. Let's go. 
      
       ! get start and stop order for this term 
       call opt_getHOforterm(term,power_disp,order_start,order_stop)
       if(order_start == 0)then 
                ! Message to Output 
                write(message,'(5a,I2,a,I2,3a)' )ch10,&
&               " ==> High order cross product terms for term ", trim(term_in%name),ch10,&
&               " ==> do not fit into specified order range from ", power_disp(1),' to ',power_disp(2),ch10,&        
&               " ==> Can not construct high order cross product bounding term",ch10
                call wrtout(ab_out,message,'COLL')
                call wrtout(std_out,message,'COLL')
               return  
       end if
       
       ! get total amount of combinations and combinations per order for the term
       call opt_getCombisforterm(order_start,order_stop,ndisp,ncombi,ncombi_order)
             
       ! Allocate terms with ncombi free space to work with 
       if(had_strain)then 
          ABI_DATATYPE_ALLOCATE(terms_out,(3*ncombi))
          ncombi_tot = 3*ncombi  
       else 
          ABI_DATATYPE_ALLOCATE(terms_out,(ncombi)) 
          ncombi_tot = ncombi 
       endif     
       ! Copy current term to the ncombination elemenst a the end in array terms 
       ! change the value of their coefficient to a start value
       ! The start is estimed to not be larger then half the initial term's value
       ! This is because higher order terms should have smaller coefficients 
       do icombi=1,ncombi_tot
           if(icombi <= ncombi)then
              coeff_ini = 1 !abs(terms_out(icombi)%coefficient / 2)
              nterm_of_term = term%nterm
              call polynomial_coeff_init(coeff_ini,nterm_of_term,terms_out(icombi),term%terms,check=.true.)
           else
              coeff_ini = 10d3
              nterm_of_term = term_in%nterm
              call polynomial_coeff_init(coeff_ini,nterm_of_term,terms_out(icombi),term_in%terms,check=.true.)
           endif
           ! Set the power of all terms we want to add to two. We find the correct power later
           ! Change the weight of the term to 1 (even terms have allways weight=1)
           do iterm_of_term=1,nterm_of_term 
             terms_out(icombi)%terms(iterm_of_term)%weight = 1
             if(icombi > ncombi .and. icombi <= 2*ncombi)then 
                terms_out(icombi)%terms(iterm_of_term)%power_strain = 2
                coeff_ini = 10d3 !abs(terms_out(icombi)%coefficient / 2)
                terms_out(icombi)%coefficient = coeff_ini
             elseif(icombi>2*ncombi)then
                terms_out(icombi)%terms(iterm_of_term)%power_strain = 4
                coeff_ini = 10d5 !abs(terms_out(icombi)%coefficient / 2)
                terms_out(icombi)%coefficient = coeff_ini
             endif 
             do idisp=1,ndisp
                terms_out(icombi)%terms(iterm_of_term)%power_disp(idisp) = 2
             enddo !idisp
           enddo !iterm_of_term 
       enddo !icombi
       
       !If term had strain we had to reinitialize it in the process 
       !Refree memory
       !if(had_strain) 
       call polynomial_coeff_free(term)  
      
       ! Get high order combinations 
       if(had_strain)then           
          call opt_getHoTerms(terms_out(:ncombi),order_start,order_stop,ndisp,ncombi_order) 
          call opt_getHoTerms(terms_out(ncombi+1:2*ncombi),order_start,order_stop,ndisp,ncombi_order) 
          call opt_getHoTerms(terms_out(2*ncombi+1:),order_start,order_stop,ndisp,ncombi_order) 
          ncombi = ncombi_tot
       else 
          call opt_getHoTerms(terms_out,order_start,order_stop,ndisp,ncombi_order) 
       endif
       !DEALLOCATION 
       ABI_DEALLOCATE(ncombi_order)

end subroutine opt_getHOcrossdisp
!!***
 
!!****f* m_opt_effpot/opt_getSingleDispTerms
!!
!! NAME
!! opt_getSingleDispTerms
!!
!! FUNCTION
!! Get polynomial terms (<polynomial_coeff_type>) with single displacements 
!! at second order inside a given range defined by the supercell size 
!!
!! INPUTS
!! sc_size(3): supercell size
!! crystal<type(crystal_t)>: all information about the crystal
!! comm: mpi communicator (at the moment only sequential tested) 
!! 
!! OUTPUT
!! terms<polynomial_coeff_type>: list single displacement polynomial_coeffs
!!
!! PARENTS
!! opt_effpotbound
!!
!! CHILDREN
!! m_polynomial_coeff.F90/polynomial_coeff_init 
!!
!! SOURCE

subroutine opt_getSingleDispTerms(terms,crystal,sc_size,comm)

 implicit none 
         
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 type(polynomial_coeff_type),allocatable,intent(inout) :: terms(:)
 type(crystal_t),intent(inout) :: crystal
 !arrays
 integer :: sc_size(3) 
!Logicals
!Strings 
!Local variables ------------------------------
!scalars
 integer :: natom,nsym,nrpt,ncoeff_sym,nstr_sym
 integer :: ncoeff,ncoeff_out,power_strph,option_GN,option
 integer :: nterms_out
 integer :: lim1,lim2,lim3
 integer :: ii,ia,ib,r1,r2,r3
 integer :: irpt,irpt_ref
 integer :: master,nproc,my_rank
 real(dp) :: norm, cutoff
!arrays
 integer :: ncell(3),power_disp(2)
 integer, allocatable :: cell(:,:)
 integer,allocatable :: list_symcoeff(:,:,:),list_symstr(:,:,:)
 !type(polynomial_coeff_type),allocatable :: terms(:)
 real(dp),allocatable :: xcart(:,:),xred(:,:),rpt(:,:)
 real(dp) :: rprimd(3,3),range_ifc(3)
 real(dp),allocatable :: dist(:,:,:,:)
 character(len=5),allocatable :: symbols(:)
!Logicals
 logical :: iam_master,need_verbose 
!Strings
 character(len=1000) :: message
!*************************************************************************

option = 2 

if(option == 1)then 
   !MPI variables
   master = 0
   nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
   iam_master = (my_rank == master)
   
   !Set/Get other variables
   natom  = crystal%natom
   nsym   = crystal%nsym
   rprimd = crystal%rprimd
   need_verbose = .TRUE.
  
   ABI_ALLOCATE(xcart,(3,natom))
   ABI_ALLOCATE(xred,(3,natom))
   xcart(:,:) = crystal%xcart(:,:)
   xred(:,:)  = crystal%xred(:,:)
  
   !Compute the max range of the ifc with respect to the trainning set
   range_ifc(:) = zero
   do ii=1,3
     norm = sqrt(rprimd(ii,1)**2+ rprimd(ii,2)**2+rprimd(ii,3)**2)
     range_ifc(ii) = range_ifc(ii) + norm * sc_size(ii) / 2.0
   end do
  
   !compute new ncell
   ncell = sc_size
   lim1=((ncell(1)/2)) + 1
   lim2=((ncell(2)/2)) + 1
   lim3=((ncell(3)/2)) + 1
   if(mod(ncell(1),2)/=0) lim1=lim1+1
   if(mod(ncell(2),2)/=0) lim2=lim2+1
   if(mod(ncell(3),2)/=0) lim3=lim3+1
   nrpt=(2*lim1+1)*(2*lim2+1)*(2*lim3+1)
  
   ncell(1) = 2*lim1+1
   ncell(2) = 2*lim2+1
   ncell(3) = 2*lim3+1
  
   !Build the rpt point
   ABI_ALLOCATE(rpt,(3,nrpt))
   ABI_ALLOCATE(cell,(3,nrpt))
  
   !WARNING:
   !Put the reference cell into the first element
   !the code will first deal with the atoms of the first cell
   irpt = 1
   irpt_ref = 1
   rpt(:,1) = zero
   cell(:,irpt)=0
   !Fill other rpt:
   do r1=lim1,-lim1,-1
     do r2=lim2,-lim2,-1
       do r3=lim3,-lim3,-1
         if(r1==0.and.r2==0.and.r3==0) then
           cycle
         end if
         irpt=irpt+1
         rpt(1,irpt)=r1*rprimd(1,1)+r2*rprimd(1,2)+r3*rprimd(1,3)
         rpt(2,irpt)=r1*rprimd(2,1)+r2*rprimd(2,2)+r3*rprimd(2,3)
         rpt(3,irpt)=r1*rprimd(3,1)+r2*rprimd(3,2)+r3*rprimd(3,3)
         cell(1,irpt)=r1;cell(2,irpt)=r2;cell(3,irpt)=r3
       end do
     end do
   end do
  
   ABI_ALLOCATE(symbols,(natom))
   call symbols_crystal(crystal%natom,crystal%ntypat,crystal%npsp,&
&                     symbols,crystal%typat,crystal%znucl)

   !Compute the distances between atoms
   !Now dist(3,ia,ib,irpt) contains the distance from atom ia to atom ib in unit cell irpt.
   ABI_ALLOCATE(dist,(3,natom,natom,nrpt))
   dist = zero
   do ia=1,natom
     do ib=1,natom
       do irpt=1,nrpt
         dist(1,ia,ib,irpt) = xcart(1,ib)-xcart(1,ia)+rpt(1,irpt)
         dist(2,ia,ib,irpt) = xcart(2,ib)-xcart(2,ia)+rpt(2,irpt)
         dist(3,ia,ib,irpt) = xcart(3,ib)-xcart(3,ia)+rpt(3,irpt)
       end do
     end do
   end do
  
   cutoff = rprimd(1,1)
  
   if(need_verbose)then
     write(message,'(1a)')' Generation of the list of all the possible coefficients'
     call wrtout(std_out,message,'COLL')
   end if
   call polynomial_coeff_getList(cell,crystal,dist,list_symcoeff,list_symstr,&
&                              natom,nstr_sym,ncoeff_sym,nrpt,range_ifc,cutoff,sc_size=sc_size)

   ABI_DEALLOCATE(dist)
   ABI_DEALLOCATE(rpt)
 
   !write(*,*) "polynomial_getList worked"

   !Get Order1 term s
   !Check difference between ncoeff_out, ncoeff 
   call polynomial_coeff_getOrder1(cell,terms,list_symcoeff,natom,nterms_out,ncoeff_sym,nrpt,nsym,symbols)
  
   !do ii=1,nterms_out
   !   write(*,*) "Term(",ii,"/",nterms_out,"): ", terms(ii)%name 
   !enddo

elseif(option == 2)then 
   !Get/set variables 
   power_disp = (/2,2/)
   power_strph = zero 
   option_GN = 0 
   sc_size = (/2,2,2/)
   cutoff = 0  

   do ii=1,3
      cutoff = cutoff + sqrt(crystal%rprimd(ii,1)**2 + &
&                            crystal%rprimd(ii,2)**2 + &
&                            crystal%rprimd(ii,3)**2)
   enddo 


call polynomial_coeff_getNorder(terms,crystal,cutoff,ncoeff,ncoeff_out,power_disp,& 
&                               power_strph,option_GN,sc_size,comm,anharmstr=.false.,spcoupling=.false.,&
&                               only_odd_power=.false.,only_even_power=.true.,verbose=.false.,&
&                               compute_symmetric=.false.) 

 !TEST MS 
 !  write(*,*) "behind call getNorder"
 !  write(*,*) "ncoeff_out: ", ncoeff_out
 !  do ii=1,ncoeff_out 
 !     write(*,*) "Term(",ii,"/",ncoeff_out,"): ", terms(ii)%name 
 !  enddo
 !TEST MS

endif !option

end subroutine opt_getSingleDispTerms
!!***

!!****f* m_opt_effpot/opt_getHOSingleDispTerms
!!
!! NAME
!! opt_getHOSingleDispTerms
!!
!! FUNCTION
!! For a givne anharmonic term that might consits of product of 
!! terms find all even high order single displacement terms
!! Example: Input term: x*y*z 
!!          Input HO-range: 6 8
!!          Output terms: x^6,x^8,y^6,y^8,z^6,z^8
!!          Caution only atomic displacements are taken into account
!!
!! INPUTS
!! term_in<polynomial_coeff_type>: input anharmonic term 
!! crystal<type(crystal_t)>: all information about the crystal
!! single_disp_terms<polynomial_coeff_out>: list of single disp terms at 
!!                                          second order to select terms from 
!! power_disp(2): Start and stop power for HO terms
!! comm: mpi communicator (at the moment only sequential tested) 
!! 
!! OUTPUT
!! terms_out<polynomial_coeff_out>: output high order even terms 
!! ncoeff: number of coefficients
!!
!! PARENTS
!! opt_effpotbound
!!
!! CHILDREN
!! m_polynomial_coeff.F90/polynomial_coeff_init 
!!
!! SOURCE

subroutine opt_getHOSingleDispTerms(term_in,terms_out,symbols,single_disp_terms,power_disp,ncoeff)

 implicit none 
         
!Arguments ------------------------------------
!scalars
 integer,intent(out) :: ncoeff
 type(polynomial_coeff_type),intent(in) :: term_in
 type(polynomial_coeff_type),intent(in) :: single_disp_terms(:)
 type(polynomial_coeff_type),allocatable,intent(out) :: terms_out(:) 
 !type(crystal_t),intent(inout) :: crystal
 !arrays
 integer, intent(in) :: power_disp(2)
 character(len=5),intent(in) :: symbols(:)
!Logicals
!Strings 
!Local variables ------------------------------
!scalars
integer :: ndisp,norder, nterm_of_term
integer :: icoeff,iorder,idisp, iterm1,iterm2,iterm3
real(dp) :: coeff_ini = 1
!Strings
 character(len=200):: name
!arrays
 type(polynomial_coeff_type),allocatable :: terms_out_tmp(:) 
!Logicals
 logical,allocatable :: found(:) 
!*************************************************************************

!Get/Set Variables 
!Number of output terms
ndisp = term_in%terms(1)%ndisp  
norder = abs(power_disp(2)-power_disp(1))/2 + 1
ncoeff = norder * ndisp
!Allocate output terms 
ABI_DATATYPE_ALLOCATE(terms_out_tmp,(ncoeff))

!find equivalent second order terms in list of single disp terms 
!for each displacement in input term
!Transfer to output term and increase order 

icoeff = 0
do idisp=1,ndisp 
   do iterm1=1,size(single_disp_terms)
      do iterm2=1,single_disp_terms(iterm1)%nterm 
         if(all(term_in%terms(1)%atindx(:,idisp) == single_disp_terms(iterm1)%terms(iterm2)%atindx(:,1)))then 
            if(all(term_in%terms(1)%cell(:,:,idisp) == single_disp_terms(iterm1)%terms(iterm2)%cell(:,:,1)))then 
               if(term_in%terms(1)%direction(idisp) == single_disp_terms(iterm1)%terms(iterm2)%direction(1))then 
                  do iorder=1,norder
                     icoeff = icoeff + 1
                     terms_out_tmp(icoeff) = single_disp_terms(iterm1)
                     nterm_of_term = single_disp_terms(iterm1)%nterm 
                     !Change order of term 
                     call polynomial_coeff_init(coeff_ini,nterm_of_term,terms_out_tmp(icoeff),&
&                                               single_disp_terms(iterm1)%terms(:),check=.true.)
                     do iterm3=1,nterm_of_term
                        terms_out_tmp(icoeff)%terms(iterm3)%power_disp = power_disp(1) + (iorder-1)*2  
                     enddo !iterm3
                  enddo !iorder
               endif
            endif
         endif
      enddo !iterm2 
   enddo !iterm1
enddo!idisp

!Change Name 
do icoeff=1,ncoeff 
!   write(std_out,*) "DEBUG icoeff: ", icoeff  
!   write(std_out,*) "Term(",icoeff,"/",ncoeff,"): ", terms_out_tmp(icoeff)%name, "name before set name" 
!   write(std_out,*) "Term(",icoeff,"/",ncoeff,") nterm: ", terms_out_tmp(icoeff)%nterm
  call polynomial_coeff_getName(name,terms_out_tmp(icoeff),symbols,recompute=.TRUE.)
  call polynomial_coeff_SetName(name,terms_out_tmp(icoeff))
!   write(*,*) "Term(",icoeff,"/",ncoeff,"): ", terms_out_tmp(icoeff)%name, "after set name" 
enddo 


!Check for doubles and delete them
!First count irreducible terms 
ABI_ALLOCATE(found,(ncoeff)) 
found = .FALSE.  
iterm3 = 0
do iterm1=1,ncoeff 
  do iterm2=iterm1+1,ncoeff 
     if(terms_out_tmp(iterm1) == terms_out_tmp(iterm2) .and. .not. found(iterm2))then 
        found(iterm2) = .TRUE.
        iterm3 = iterm3 + 1
     endif
  enddo 
enddo 

iterm3 = ncoeff - iterm3
ABI_DATATYPE_ALLOCATE(terms_out,(iterm3))

!Second copy them 
iterm3 = 0
do iterm1=1,ncoeff 
   if(.not. found(iterm1))then 
      iterm3 = iterm3 + 1
!      terms_out(iterm3) = terms_out_tmp(iterm1)
      call polynomial_coeff_init(coeff_ini,terms_out_tmp(iterm1)%nterm,terms_out(iterm3),&
&                                terms_out_tmp(iterm1)%terms,terms_out_tmp(iterm1)%name,check=.TRUE.)
   endif 
enddo

!ABI_DEALLOCATE(terms_out_tmp)
call polynomial_coeff_list_free(terms_out_tmp)
ABI_DEALLOCATE(found)
ncoeff = iterm3

 !TEST MS 
 !  write(*,*) "behind call getNorder"
 !  write(*,*) "ncoeff_out: ", ncoeff_out
 !  do ii=1,ncoeff_out 
 !     write(*,*) "Term(",ii,"/",ncoeff_out,"): ", terms(ii)%name 
 !  enddo
 !TEST MS

!!Check after reduction
!do icoeff=1,ncoeff 
!   write(std_out,*) "DEBUG icoeff: ", icoeff  
!   write(*,*) "Term(",icoeff,"/",ncoeff,"): ", terms_out(icoeff)%name, "name before set name" 
!  call polynomial_coeff_getName(name,terms_out(icoeff),symbols,recompute=.TRUE.)
!  call polynomial_coeff_SetName(name,terms_out(icoeff))
!   write(*,*) "Term(",icoeff,"/",ncoeff,"): ", terms_out(icoeff)%name, "after set name" 
!enddo 

end subroutine opt_getHOSingleDispTerms
!!***



!!****f* m_opt_effpot/opt_boundcoeff
!! NAME
!! opt_boundcoeff
!!
!! FUNCTION
!!
!! optimize a bound coefficient if optimized value is negative 
!! put an positive value that respects a precision penalty  
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

function opt_boundcoeff(yvalues,cvalues) result (coeff)
!Arguments ------------------------------------
 implicit none

!Arguments ------------------------------------
  real(dp),intent(in) :: yvalues(2),cvalues(2)
  real(dp) :: coeff
!local
!variable
 real(dp) :: a,b,coeff_tmp,x1,x2,penalty
!array
! *************************************************************************
 
 a = ( (yvalues(1) - 1) - (yvalues(2)-1)*(cvalues(1)/cvalues(2))) / (cvalues(1)**2 - cvalues(1)*cvalues(2)) 
 
 b = ( (yvalues(2) - 1)/cvalues(2) ) - ( (yvalues(1) -1)*cvalues(2) - (yvalues(2) - 1)*cvalues(1) )&
&    / (cvalues(1)**2 - cvalues(1)*cvalues(2)) 
 
 !write(*,*) "a", a
 !write(*,*) "b", b
 penalty = 0.001 
 coeff_tmp = -b/(2*a)
 !write(*,*) "coeff_tmp", coeff_tmp 
 if(coeff_tmp > 0)then
   coeff = coeff_tmp 
 elseif(coeff_tmp <= 0)then 
   x1 = (-b + sqrt(b**2 + 4*a*penalty)) / (2*a) ! 1.001 penalty value 
   x2 = (-b - sqrt(b**2 + 4*a*penalty)) / (2*a)
   !write(*,*) "x1", x1 
   !write(*,*) "x2", x2
   if(x1>0)then 
     coeff = x1 
   else 
     coeff = x2 
   endif
 endif

end function opt_boundcoeff
!!***

!!****f* m_opt_effpot/check_to_skip
!! NAME
!! check_to_skip
!!
!! FUNCTION
!!
!! Check if term contains only bodies with even power
!! and has a positive coefficient. If yes term doesn't need 
!! a bounding high order equivalent and we can skip it. 
!! Function retursn logical to_skip.
!! 
!! INPUTS
!!
!! term<polynomial_coeff_type>:anharmonic term 
!!
!! OUTPUT
!!
!! logical: to_skip
!!
!! SOURCE

function check_to_skip(term) result (to_skip)
!Arguments ------------------------------------
 implicit none
 type(polynomial_coeff_type),intent(in) :: term
 logical :: to_skip
! ------------------------------------
!local
!variable
 character(len=1000) :: message
!array
! *************************************************************************

to_skip = .FALSE. 
! Let's check if we really want all this mess 
! If the term is even and its coefficient positive we skip it. Also here we take terms(1) as example for all equivalent terms of term
if(term%coefficient > 0 .and. .not. any(mod(term%terms(1)%power_disp(:),2) /= 0))then
   if(.not. any(mod(term%terms(1)%power_strain(:),2) /= 0))then 
         ! Message to Output 
         write(message,'(3a)' )ch10,&
&         ' ==> No need for high order bounding term',ch10
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')
        to_skip = .TRUE. 
        return  
   end if
end if

end function check_to_skip
!!***

end module m_opt_effpot
!!***

