!!****m* ABINIT/m_fit_polynomial_coeff
!!
!! NAME
!! m_fit_polynomial_coeff
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2010-2019 ABINIT group (AM)
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

module m_fit_polynomial_coeff

 use defs_basis
 use m_errors
 use m_abicore
 use m_polynomial_coeff
 use m_atomdata
 use m_xmpi
 use m_supercell

 use m_special_funcs,only : factorial
 use m_geometry,       only : xred2xcart
 use m_crystal,only : symbols_crystal
 use m_strain,only : strain_type,strain_get
 use m_effective_potential,only : effective_potential_type, effective_potential_evaluate
 use m_effective_potential,only : effective_potential_freeCoeffs,effective_potential_setCoeffs
 use m_effective_potential,only : effective_potential_getDisp, effective_potential_writeAnhHead
 use m_effective_potential_file, only : effective_potential_file_mapHistToRef
 use m_io_tools,   only : open_file,get_unit
 use m_abihist, only : abihist,abihist_free,abihist_init,abihist_copy,write_md_hist,var2hist
 use m_random_zbq
 use m_fit_data
 use m_geometry, only: metric
 use m_scup_dataset 
#if defined DEV_MS_SCALEUP 
 use scup_global, only : global_set_parent_iter,global_set_print_parameters 
#endif 

 implicit none

 public :: fit_polynomial_coeff_computeGF
 public :: fit_polynomial_coeff_computeMSD
 public :: fit_polynomial_coeff_fit
 public :: fit_polynomial_coeff_getFS
 public :: fit_polynomial_coeff_getPositive
 public :: fit_polynomial_coeff_getCoeffBound
 public :: fit_polynomial_coeff_solve
 public :: fit_polynomial_coeff_testEffPot
 public :: fit_polynomial_printSystemFiles
 public :: genereList
!!***

CONTAINS  
!===========================================================================================


!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_fit
!!
!! NAME
!! fit_polynomial_coeff_fit
!!
!! FUNCTION
!! Fit the list of coefficients included in eff_pot,
!! if the coefficients are not set in eff_pot, this routine will genenerate
!! a list of coefficients by taking into acount the symmetries of the system
!! and the cutoff
!!
!! INPUTS
!! eff_pot<type(effective_potential)> = effective potential
!! bancoeff(nbancoeff) = list of bannned coeffcients, these coefficients will NOT be
!!                       used during the fit process
!! fixcoeff(nfixcoeff) = list of fixed coefficient, these coefficients will be
!!                       imposed during the fit process
!! hist<type(abihist)> = The history of the MD (or snapshot of DFT)
!! generateterm = term to activate the generation of the term set
!! power_disps(2) = array with the minimal and maximal power_disp to be computed
!! nbancoeff = number of banned coeffcients
!! ncycle_in = number of maximum cycle (maximum coefficient to be fitted)
!! nfixcoeff = Number of coefficients imposed during the fit process
!! option = option of the fit process : 1 - selection of the coefficient one by one
!!                                      2 - selection of the coefficients with Monte Carlo(testversion)
!! comm = MPI communicator
!! cutoff_in = optional,cut off to apply to the range of interation if
!!           the coefficient are genereted in this routine
!! max_power_strain = maximum order of the strain of the strain phonon coupling
!! fit_initializeData = optional, logical !If true, we store all the informations for the fit,
!!                      it will reduce the computation time but increase a lot the memory...
!! fit_tolMSDF = optional, tolerance in eV^2/A^2 on the Forces for the fit process
!! fit_tolMSDS = optional, tolerance in eV^2/A^2 on the Stresses for the fit process
!! fit_tolMSDE = optional, tolerance in meV^2/A^2 on the Energy for the fit process
!! fit_tolMSDFS= optional, tolerance in eV^2/A^2 on the Forces+stresses for the fit process
!! positive = optional, TRUE will return only positive coefficients
!!                      FALSE, default
!! verbose  = optional, flag for the verbose mode
!! anhstr = logical, optional : TRUE, the anharmonic strain are computed
!!                              FALSE, (default) the anharmonic strain are not computed
!! only_odd_power = logical, optional : if TRUE generate only odd power
!! only_even_power= logical, optional : if TRUE generate only even power
!!
!! OUTPUT
!! eff_pot<type(effective_potential)> = effective potential datatype with new fitted coefficients
!!
!! PARENTS
!!      m_fit_polynomial_coeff,mover_effpot,multibinit
!!
!! CHILDREN
!!      destroy_supercell,generelist,init_supercell,xred2xcart
!!
!! SOURCE

subroutine fit_polynomial_coeff_fit(eff_pot,bancoeff,fixcoeff,hist,generateterm,power_disps,&
&                                   nbancoeff,ncycle_in,nfixcoeff,option,comm,cutoff_in,&
&                                   max_power_strain,initialize_data,&
&                                   fit_tolMSDF,fit_tolMSDS,fit_tolMSDE,fit_tolMSDFS,&
&                                   positive,verbose,anharmstr,spcoupling,&
&                                   only_odd_power,only_even_power,prt_names,prt_anh,& 
&                                   fit_iatom)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncycle_in,nfixcoeff,comm
 integer,intent(in) :: generateterm,nbancoeff,option
!arrays
 integer,intent(in) :: fixcoeff(nfixcoeff), bancoeff(nbancoeff)
 integer,intent(in) :: power_disps(2)
 type(effective_potential_type),target,intent(inout) :: eff_pot
 type(abihist),intent(inout) :: hist
 integer,optional,intent(in) :: max_power_strain,prt_names,prt_anh,fit_iatom
 real(dp),optional,intent(in) :: cutoff_in,fit_tolMSDF,fit_tolMSDS,fit_tolMSDE,fit_tolMSDFS
 logical,optional,intent(in) :: verbose,positive,anharmstr,spcoupling
 logical,optional,intent(in) :: only_odd_power,only_even_power
 logical,optional,intent(in) :: initialize_data
!Local variables-------------------------------
!scalar
 integer :: ii,icoeff,my_icoeff,icycle,icycle_tmp,ierr,info,index_min,iproc,isweep,jcoeff
 integer :: master,max_power_strain_in,my_rank,my_ncoeff,ncoeff_model,ncoeff_tot,natom_sc,ncell,ncycle
 integer :: ncycle_tot,ncycle_max,need_prt_names,nproc,ntime,nsweep,size_mpi,ncoeff_fix
 integer :: rank_to_send,unit_names,unit_anh,fit_iatom_in
 real(dp) :: cutoff,factor,time,tolMSDF,tolMSDS,tolMSDE,tolMSDFS
 real(dp),parameter :: HaBohr_meVAng = 27.21138386 / 0.529177249
 logical :: iam_master,need_verbose,need_positive,converge,file_opened
 logical :: need_anharmstr,need_spcoupling,ditributed_coefficients,need_prt_anh
 logical :: need_only_odd_power,need_only_even_power,need_initialize_data
!arrays
 real(dp) :: mingf(4)
 integer :: sc_size(3)
 integer,allocatable  :: buffsize(:),buffdisp(:),buffin(:)
 integer,allocatable  :: list_coeffs(:),list_coeffs_tmp(:),list_coeffs_tmp2(:)
 integer,allocatable  :: my_coeffindexes(:),singular_coeffs(:)
 integer,allocatable  :: my_coefflist(:) ,stat_coeff(:)
 real(dp),allocatable :: buffGF(:,:),coeff_values(:),energy_coeffs(:,:)
 real(dp),allocatable :: energy_coeffs_tmp(:,:)
 real(dp),allocatable :: fcart_coeffs(:,:,:,:),gf_values(:,:),gf_mpi(:,:)
 real(dp),allocatable :: fcart_coeffs_tmp(:,:,:,:),strten_coeffs_tmp(:,:,:)
 real(dp),allocatable :: strten_coeffs(:,:,:)
 type(polynomial_coeff_type),allocatable :: my_coeffs(:)
 type(polynomial_coeff_type),target,allocatable :: coeffs_tmp(:)
 type(polynomial_coeff_type),pointer :: coeffs_in(:)
 type(fit_data_type) :: fit_data
 character(len=1000) :: message
 character(len=fnlen) :: filename
 character(len=5) :: powerstr,rangestr
 character(len=200) :: namefile
 character(len=3)  :: i_char
 character(len=7)  :: j_char
 character(len=5),allocatable :: symbols(:)
! *************************************************************************

!MPI variables
 master = 0
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

!Initialisation of optional arguments
 need_verbose = .TRUE.
 if(present(verbose)) need_verbose = verbose
 need_initialize_data = .TRUE.
 if(present(initialize_data)) need_initialize_data = initialize_data
 need_positive = .FALSE.
 if(present(positive)) need_positive = positive
 need_anharmstr = .FALSE.
 if(present(anharmstr)) need_anharmstr = anharmstr
 need_spcoupling = .TRUE.
 if(present(spcoupling)) need_spcoupling = spcoupling
 need_only_odd_power = .FALSE.
 if(present(only_odd_power)) need_only_odd_power = only_odd_power
 need_prt_names = 0
 if(present(prt_names)) need_prt_names = prt_names
 need_prt_anh = .FALSE. 
 if(present(prt_anh))then
   if(prt_anh == 1) need_prt_anh = .TRUE.
 end if
 need_only_even_power = .FALSE.
 if(present(only_even_power)) need_only_even_power = only_even_power
 if(need_only_odd_power.and.need_only_even_power)then
      write(message, '(3a)' )&
&       'need_only_odd_power and need_only_even_power are both true',ch10,&
&       'Action: contact abinit group'
   MSG_ERROR(message)
 end if
 max_power_strain_in = 1
 if(present(max_power_strain))then
   max_power_strain_in = max_power_strain
 end if
 if(max_power_strain_in <= 0)then
      write(message, '(3a)' )&
&       'max_power_strain can not be inferior or equal to zero',ch10,&
&       'Action: contact abinit group'
   MSG_ERROR(message)
 end if
 !Check which atom to fit, if not present do all atoms 
 if(present(fit_iatom))then 
    fit_iatom_in = fit_iatom 
 else 
    fit_iatom_in = -1 
 endif


 ABI_ALLOCATE(symbols,(eff_pot%crystal%natom))
 call symbols_crystal(eff_pot%crystal%natom,eff_pot%crystal%ntypat,eff_pot%crystal%npsp,&
&                     symbols,eff_pot%crystal%typat,eff_pot%crystal%znucl)

!Set the tolerance for the fit
 tolMSDF=zero;tolMSDS=zero;tolMSDE=zero;tolMSDFS=zero
 if(present(fit_tolMSDF)) tolMSDF  = fit_tolMSDF
 if(present(fit_tolMSDS)) tolMSDS  = fit_tolMSDS
 if(present(fit_tolMSDE)) tolMSDE  = fit_tolMSDE
 if(present(fit_tolMSDFS))tolMSDFS = fit_tolMSDFS

 if(need_verbose) then
   write(message,'(a,(80a))') ch10,('=',ii=1,80)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   if(fit_iatom_in > 0)then
     write(message,'(2a,I3,2a)') ch10,' Starting Fit process around atom', fit_iatom_in,": ", symbols(fit_iatom)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
   else 
     write(message,'(2a)') ch10,' Starting Fit process with all possible cross-terms'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
   endif
   write(message,'(a,(80a))') ch10,('-',ii=1,80)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

 ditributed_coefficients = .true.
 if(option==2) ditributed_coefficients = .false.

!if the number of atoms in reference supercell into effpot is not correct,
!wrt to the number of atom in the hist, we set map the hist and set the good supercell
 if (size(hist%xred,2) /= eff_pot%supercell%natom) then
   call effective_potential_file_mapHistToRef(eff_pot,hist,comm,verbose=need_verbose)
 end if

!Set the cut off
 cutoff = zero
 if(present(cutoff_in))then
   cutoff = cutoff_in
 end if
!If the cutoff is set to zero, we define a default value
 if(abs(cutoff)<tol16)then
   do ii=1,3
     cutoff = cutoff + sqrt(eff_pot%supercell%rprimd(ii,1)**2+&
&                           eff_pot%supercell%rprimd(ii,2)**2+&
&                           eff_pot%supercell%rprimd(ii,3)**2)
   end do
   cutoff = cutoff / 3.0_dp
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


!Get the list of coefficients to fit:
!get from the eff_pot type (from the input)
!or
!regenerate the list
 my_ncoeff = 0
 ncoeff_tot = 0
 ncoeff_model = eff_pot%anharmonics_terms%ncoeff

!Reset ncoeff_tot
 if(ncoeff_model > 0)then
   if(need_verbose)then
     write(message, '(4a)' )ch10,' The coefficients present in the effective',&
&    ' potential will be used for the fit'
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')
   end if
 end if

 if(generateterm == 1)then
! we need to regerate them
   if(need_verbose)then
     write(message, '(4a)' )ch10,' The coefficients for the fit will be generated'
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')

     write(message,'(a,F6.3,a)') " Cutoff of ",cutoff," Bohr is imposed"
     call wrtout(std_out,message,'COLL')
   end if

   call polynomial_coeff_getNorder(coeffs_tmp,eff_pot%crystal,cutoff,my_ncoeff,ncoeff_tot,power_disps,&
&                                  max_power_strain_in,0,sc_size,comm,anharmstr=need_anharmstr,&
&                                  spcoupling=need_spcoupling,distributed=.true.,&
&                                  only_odd_power=need_only_odd_power,&
&                                  only_even_power=need_only_even_power,& 
&                                  fit_iatom=fit_iatom_in)
 end if
!Copy the initial coefficients from the model on the CPU 0
 ncoeff_tot = ncoeff_tot + ncoeff_model
 if(iam_master .and. ncoeff_model > 0) my_ncoeff = my_ncoeff + ncoeff_model

!Get number of fixed coeff
 ncoeff_fix = 0  
 if(nfixcoeff /=0) then 
   if(nfixcoeff == -1)then 
      ncoeff_fix = ncoeff_model 
   else 
      ncoeff_fix = nfixcoeff 
   endif 
 endif 

!Get the list with the number of coeff on each CPU
!In order to be abble to compute the my_coeffindexes array which is for example:
! if CPU0 has 200  Coeff and CPU1 has 203 Coeff then
! for CPU0:my_coeffindexes=>1-200 and for CPU1:my_coeffindexes=>201-403
!Also fill the my_coeffs array with the generated coefficients and/or the coefficient from the input xml
 ABI_ALLOCATE(buffin,(nproc))
 buffin = 0
 buffin(my_rank+1) = my_ncoeff
 call xmpi_sum(buffin,comm,ierr)
 ABI_ALLOCATE(my_coeffindexes,(my_ncoeff))
 ABI_ALLOCATE(my_coefflist,(my_ncoeff))
 ABI_DATATYPE_ALLOCATE(my_coeffs,(my_ncoeff))
 do icoeff=1,my_ncoeff

   jcoeff = icoeff
   my_coefflist(icoeff) = icoeff

   if(my_rank==0) then
     my_coeffindexes(icoeff) = icoeff
   else
     my_coeffindexes(icoeff) = sum(buffin(1:my_rank)) + icoeff
   end if

!  Only copy the input coefficients on the CPU0
   if(my_rank==0) then
     if(icoeff <= ncoeff_model)then
       coeffs_in => eff_pot%anharmonics_terms%coefficients
     else
       coeffs_in => coeffs_tmp
       jcoeff = jcoeff-ncoeff_model
     end if
   else
     coeffs_in => coeffs_tmp
   end if
   call polynomial_coeff_init(one,coeffs_in(jcoeff)%nterm,&
&                             my_coeffs(icoeff),coeffs_in(jcoeff)%terms,&
&                             coeffs_in(jcoeff)%name,&
&                             check=.true.)
   call polynomial_coeff_free(coeffs_in(jcoeff))
 end do

!Deallocation
 if(allocated(coeffs_tmp)) then
   ABI_DATATYPE_DEALLOCATE(coeffs_tmp)
 end if
 NULLIFY(coeffs_in)
 ABI_DEALLOCATE(buffin)

 !wait everybody
 call xmpi_barrier(comm)
  
 if(need_prt_names == 1 .and. nproc == 1)then
   unit_names = get_unit()
   write (powerstr,'(I0,A1,I0)') power_disps(1),'-',power_disps(2)
   write (rangestr,'(F4.2)') cutoff 
   namefile='name-of-terms_range-'//trim(rangestr)//'_power-'//trim(powerstr)//'.out'
   namefile=trim(namefile)
   write(message,'(a)') " Printing of list of terms is asked"
   call wrtout(std_out,message,'COLL')
   write(message,'(a,a)') " Write list of generated terms to file: ",namefile
   call wrtout(std_out,message,'COLL')
   open(unit_names,file=namefile,status='replace')
   do icoeff=1,ncoeff_tot
       write(unit_names,*) icoeff, trim(my_coeffs(icoeff)%name ) ! Marcus Write name of coefficient to file
   enddo
   close(unit_names)
 else if(need_prt_names == 1 .and. nproc /= 1)then
   write(message, '(15a)' )ch10,&
&        ' --- !WARNING',ch10,&
&        '     The printing of the list of generated Terms has been requested.',ch10,&
&        '     This option is currently limited to serial execution of multibinit ',ch10,&
&        '     The terms are not printed.',ch10,&
&        '     Action: Rerun in serial.',ch10,&
&        ' ---',ch10
     call wrtout(std_out,message,"COLL")
 endif


!Reset the output (we free the memory)
 call effective_potential_freeCoeffs(eff_pot)


!Check if ncycle_in is not zero or superior to ncoeff_tot
 if(need_verbose.and.(ncycle_in > ncoeff_tot).or.(ncycle_in<0.and.nfixcoeff /= -1)) then
   write(message, '(6a,I0,3a)' )ch10,&
&        ' --- !WARNING',ch10,&
&        '     The number of cycle requested in the input is not correct.',ch10,&
&        '     This number will be set to the maximum of coefficients: ',ncoeff_tot,ch10,&
&        ' ---',ch10
     call wrtout(std_out,message,"COLL")
   end if

!Use fixcoeff
!ncycle_tot store the curent number of coefficient in the model
!Do not reset this variable...
 ncycle_tot = 0
 if (nfixcoeff == -1)then
   write(message, '(3a)')' nfixcoeff is set to -1, the coefficients present in the model',&
&                        ' are imposed.',ch10
   ncycle_tot = ncycle_tot + ncoeff_model
 else
   if (nfixcoeff > 0)then
     if(maxval(fixcoeff(:)) > ncoeff_tot) then
       write(message, '(4a,I0,6a)' )ch10,&
&        ' --- !WARNING',ch10,&
&        '     The value ',maxval(fixcoeff(:)),' is not in the list.',ch10,&
&        '     Start from scratch...',ch10,&
&        ' ---',ch10
     else
       ncycle_tot = ncycle_tot + nfixcoeff
       write(message, '(2a)')' Some coefficients are imposed from the input.',ch10
     end if
   else
     write(message, '(4a)')' There is no coefficient imposed from the input.',ch10,&
&                        ' Start from scratch',ch10
   end if
 end if

 if(need_verbose) call wrtout(std_out,message,'COLL')

!Compute the number of cycle:
 ncycle     = ncycle_in
!Compute the maximum number of cycle
 ncycle_max = ncycle_in + ncycle_tot

!Check if the number of request cycle + the initial number of coeff is superior to
 !the maximum number of coefficient allowed
 if(ncycle_max > ncoeff_tot) then
   ncycle = ncoeff_tot - ncycle_tot
   ncycle_max = ncoeff_tot
   write(message, '(4a,I0,2a,I0,2a,I0,3a)' )ch10,&
&      ' --- !WARNING',ch10,&
&      '     The number of cycle + the number of imposed coefficients: ',ncycle_max,ch10,&
&      '     is superior to the maximum number of coefficients in the initial list: ',ncoeff_tot,ch10,&
&      '     The number of cycle is set to ',ncycle,ch10,&
&      ' ---',ch10
   if(need_verbose) call wrtout(std_out,message,'COLL')
 else if (option==2)then
!  Always set to the maximum
   ncycle_max = ncoeff_tot
 end if

!Initialisation of constants
 ntime    = hist%mxhist
 natom_sc = eff_pot%supercell%natom
 ncell    = eff_pot%supercell%ncells
 factor   = 1._dp/natom_sc

!Initialisation of arrays:
 ABI_ALLOCATE(energy_coeffs_tmp,(ncycle_max,ntime))
 ABI_ALLOCATE(list_coeffs,(ncycle_max))
 ABI_ALLOCATE(fcart_coeffs_tmp,(3,natom_sc,ncycle_max,ntime))
 ABI_ALLOCATE(strten_coeffs_tmp,(6,ntime,ncycle_max))
 list_coeffs  = 0

!if ncycle_tot > 0 fill list_coeffs with the fixed coefficients
 if(ncycle_tot > 0)then
   do ii = 1,ncycle_tot
     if(nfixcoeff == -1)then
       if(ii <= ncoeff_model)then
         list_coeffs(ii) = ii
       end if
     else
       list_coeffs(ii) = fixcoeff(ii)
     end if
   end do
 end if

!Get the decomposition for each coefficients of the forces and stresses for
!each atoms and each step  equations 11 & 12 of  PRB95,094115(2017) [[cite:Escorihuela-Sayalero2017]]
 if(need_verbose)then
   write(message, '(a)' ) ' Initialisation of the fit process...'
   call wrtout(std_out,message,'COLL')
 end if
!Before the fit, compute constants with fit_data_compute.
!Conpute the strain of each configuration.
!Compute the displacmeent of each configuration.
!Compute the variation of the displacement due to strain of each configuration.
!Compute fixed forces and stresse and get the standard deviation.
!Compute Sheppard and al Factors  \Omega^{2} see J.Chem Phys 136, 074103 (2012) [[cite:Sheppard2012]].
 call fit_data_compute(fit_data,eff_pot,hist,comm,verbose=need_verbose)

!Get the decomposition for each coefficients of the forces,stresses and energy for
!each atoms and each step  (see equations 11 & 12 of  
!PRB95,094115(2017)) [[cite:Escorihuela-Sayalero2017]]+ allocation
!If the user does not turn off this initialization, we store all the informations for the fit,
!it will reduce the computation time but increase a lot the memory...
 if(need_initialize_data)then
   ABI_ALLOCATE(energy_coeffs,(my_ncoeff,ntime))
   ABI_ALLOCATE(fcart_coeffs,(3,natom_sc,my_ncoeff,ntime))
   ABI_ALLOCATE(strten_coeffs,(6,ntime,my_ncoeff))
   call fit_polynomial_coeff_getFS(my_coeffs,fit_data%training_set%du_delta,&
&                                 fit_data%training_set%displacement,&
&                                 energy_coeffs,fcart_coeffs,natom_sc,eff_pot%crystal%natom,&
&                                 my_ncoeff,ntime,sc_size,fit_data%training_set%strain,&
&                                 strten_coeffs,fit_data%training_set%ucvol,my_coefflist,my_ncoeff)
 else
!  Allocate just 1 dimension ! Save MEMORY !
   ABI_ALLOCATE(energy_coeffs,(1,ntime))
   ABI_ALLOCATE(fcart_coeffs,(3,natom_sc,1,ntime))
   ABI_ALLOCATE(strten_coeffs,(6,ntime,1))   
 end if

!Allocation of arrays
 ABI_DATATYPE_ALLOCATE(coeffs_tmp,(ncycle_max))
 ABI_ALLOCATE(singular_coeffs,(max(1,my_ncoeff)))
 ABI_ALLOCATE(coeff_values,(ncycle_max))
 ABI_ALLOCATE(gf_values,(4,max(1,my_ncoeff)))
 ABI_ALLOCATE(list_coeffs_tmp,(ncycle_max))
 ABI_ALLOCATE(list_coeffs_tmp2,(ncycle_max))
 ABI_ALLOCATE(stat_coeff,(ncoeff_tot))
 coeff_values = zero
 singular_coeffs = 0
 stat_coeff = 0
!Set mpi buffer
!Set the bufsize for mpi allgather
 ABI_ALLOCATE(buffsize,(nproc))
 ABI_ALLOCATE(buffdisp,(nproc))
 ABI_ALLOCATE(buffGF,(5,1))
 ABI_ALLOCATE(gf_mpi,(5,nproc))
 buffsize(:) = 0
 buffdisp(1) = 0
 do ii= 1,nproc
   buffsize(ii) =  5
 end do
 do ii = 2,nproc
   buffdisp(ii) = buffdisp(ii-1) + buffsize(ii-1)
 end do
 size_mpi = 5*nproc
!If some coeff are imposed by the input, we need to fill the arrays
!with this coeffs and broadcast to the others CPUs :
 if(ncycle_tot>=1)then
   do icycle = 1,ncycle_tot
     list_coeffs_tmp(icycle) = icycle
     rank_to_send = 0
     do icoeff=1,my_ncoeff
       if((my_coeffindexes(icoeff)==list_coeffs(icycle)))then

         if(need_initialize_data)then
           my_icoeff = icoeff
         else
           my_icoeff = 1
!          Need to initialized the data for the fit for this coefficient 
           call fit_polynomial_coeff_getFS(my_coeffs,fit_data%training_set%du_delta,&
&                                          fit_data%training_set%displacement,&
&                                          energy_coeffs,fcart_coeffs,natom_sc,eff_pot%crystal%natom,&
&                                          my_ncoeff,ntime,sc_size,fit_data%training_set%strain,&
&                                          strten_coeffs,fit_data%training_set%ucvol,&
&                                          my_coefflist(icoeff),1)
         end if
         
         energy_coeffs_tmp(icycle,:)    = energy_coeffs(my_icoeff,:)
         fcart_coeffs_tmp(:,:,icycle,:) = fcart_coeffs(:,:,my_icoeff,:)
         strten_coeffs_tmp(:,:,icycle)  = strten_coeffs(:,:,my_icoeff)
         rank_to_send = my_rank
         call polynomial_coeff_free(coeffs_tmp(icycle))
         call polynomial_coeff_init(coeff_values(icycle),my_coeffs(icoeff)%nterm,&
&                                   coeffs_tmp(icycle),my_coeffs(icoeff)%terms,&
&                                   my_coeffs(icoeff)%name,&
&                                   check=.false.)
         exit
       end if
     end do
!    Need to send the rank with the chosen coefficient
     call xmpi_sum(rank_to_send, comm, ierr)
!    Boadcast the coefficient
     call xmpi_bcast(energy_coeffs_tmp(icycle,:), rank_to_send, comm, ierr)
     call xmpi_bcast(fcart_coeffs_tmp(:,:,icycle,:) , rank_to_send, comm, ierr)
     call xmpi_bcast(strten_coeffs_tmp(:,:,icycle), rank_to_send, comm, ierr)
     call polynomial_coeff_broadcast(coeffs_tmp(icycle), rank_to_send, comm)
   end do
 end if

!Waiting for all
 if(nproc > 1)  then
   if(need_verbose)then
     write(message, '(a)') ' Initialisation done... waiting for all the CPU'
     call wrtout(std_out,message,'COLL')
   end if
   call xmpi_barrier(comm)
 end if

!Compute GF, coeff_values,strten_coeffs and fcart_coeffs are set to zero
!it means that only the harmonic part wiil be computed
 coeff_values = zero

 call fit_polynomial_coeff_computeGF(coeff_values,energy_coeffs,fit_data%energy_diff,fcart_coeffs,&
&                                    fit_data%fcart_diff,gf_values(:,1),int((/1/)),natom_sc,&
&                                    0,my_ncoeff,ntime,strten_coeffs,fit_data%strten_diff,&
&                                    fit_data%training_set%sqomega)

!Print the standard deviation before the fit
 write(message,'(3a,ES24.16,4a,ES24.16,2a,ES24.16,2a,ES24.16,a)' ) &
&                   ' Mean Standard Deviation values at the begining of the fit process (meV/atm):',&
&               ch10,'   Energy          : ',&
&               gf_values(4,1)*Ha_EV*1000*factor  ,ch10,&
&                    ' Goal function values at the begining of the fit process (eV^2/A^2):',ch10,&
&                    '   Forces+Stresses : ',&
&               gf_values(1,1)*(HaBohr_meVAng)**2,ch10,&
&                    '   Forces          : ',&
&               gf_values(2,1)*(HaBohr_meVAng)**2,ch10,&
&                    '   Stresses        : ',&
&               gf_values(3,1)*(HaBohr_meVAng)**2,ch10
 if(need_verbose)then
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

 select case(option)

 case(1)
   !Option 1, we select the coefficients one by one
   if(need_verbose.and.ncycle > 0)then
     write(message,'(a,3x,a,10x,a,14x,a,14x,a,14x,a)') " N","Selecting","MSDE","MSDFS","MSDF","MSDS"
     call wrtout(ab_out,message,'COLL')
     write(message,'(4x,a,6x,a,8x,a,8x,a,8x,a)') "Coefficient","(meV/atm)","(eV^2/A^2)","(eV^2/A^2)",&
&                                            "(eV^2/A^2)"
     call wrtout(ab_out,message,'COLL')
   end if

!  Start fit process
   do icycle_tmp = 1,ncycle
     icycle = ncycle_tot + 1
     list_coeffs_tmp(icycle) = icycle
     if(need_verbose)then
       write(message, '(4a,I0,a)')ch10,'--',ch10,' Try to find the best model with ',&
&                                 icycle,' coefficient'
       if(icycle > 1)  write(message, '(2a)') trim(message),'s'
       if(nproc > 1)  then
         if(my_ncoeff>=1) then
           write(message, '(2a,I0,a)')trim(message), ' (only the ',my_ncoeff,&
&                                                ' first are printed for this CPU)'
         else
           write(message, '(2a)')trim(message), ' (no coefficient treated by this CPU)'
         end if
       end if
       call wrtout(std_out,message,'COLL')
       if(icycle>1 .or. any(list_coeffs(:) > zero))then
         write(message, '(3a)') ' The coefficient numbers from the previous cycle are:',ch10,' ['
         do ii=1,icycle-1
           if(ii<icycle-1)then
             write(message, '(a,I0,a)') trim(message),list_coeffs(ii),','
           else
             write(message, '(a,I0)') trim(message),list_coeffs(ii)
           end if
         end do
         write(message, '(3a)') trim(message),']',ch10
         call wrtout(std_out,message,'COLL')
       end if

       write(message,'(2x,a,12x,a,14x,a,13x,a,14x,a)') " Testing","MSDE","MSDFS","MSDF","MSDS"
       call wrtout(std_out,message,'COLL')
       write(message,'(a,7x,a,8x,a,8x,a,8x,a)') " Coefficient","(meV/atm)","(eV^2/A^2)","(eV^2/A^2)",&
&                                            "(eV^2/A^2)"
       call wrtout(std_out,message,'COLL')
     end if!End if verbose

!    Reset gf_values
     gf_values(:,:) = zero
     do icoeff=1,my_ncoeff
!    cycle if this coefficient is not allowed
       if(any(list_coeffs==my_coeffindexes(icoeff)).or.singular_coeffs(icoeff) == 1) cycle
       if(nbancoeff >= 1)then
         if(any(bancoeff==my_coeffindexes(icoeff))) cycle
       end if
       list_coeffs(icycle) = my_coeffindexes(icoeff)

       if(need_initialize_data)then
         my_icoeff = icoeff
       else
!        Need to initialized the data for the fit for this coefficient
         my_icoeff = 1        
         call fit_polynomial_coeff_getFS(my_coeffs,fit_data%training_set%du_delta,&
&                                        fit_data%training_set%displacement,&
&                                        energy_coeffs,fcart_coeffs,natom_sc,eff_pot%crystal%natom,&
&                                        my_ncoeff,ntime,sc_size,fit_data%training_set%strain,&
&                                        strten_coeffs,fit_data%training_set%ucvol,&
&                                        my_coefflist(icoeff),1)
       end if
       
!      Fill the temporary arrays
       energy_coeffs_tmp(icycle,:)    = energy_coeffs(my_icoeff,:)
       fcart_coeffs_tmp(:,:,icycle,:) = fcart_coeffs(:,:,my_icoeff,:)
       strten_coeffs_tmp(:,:,icycle)  = strten_coeffs(:,:,my_icoeff)

!      call the fit process routine
!      This routine solves the linear system proposed 
!      by C.Escorihuela-Sayalero see PRB95,094115(2017) [[cite:Escorihuela-Sayalero2017]]
       call fit_polynomial_coeff_solve(coeff_values(1:icycle),fcart_coeffs_tmp,fit_data%fcart_diff,&
&                                      energy_coeffs_tmp,fit_data%energy_diff,info,&
&                                      list_coeffs_tmp(1:icycle),natom_sc,icycle,ncycle_max,ntime,&
&                                      strten_coeffs_tmp,fit_data%strten_diff,&
&                                      fit_data%training_set%sqomega)
       if(info==0)then
         if (need_positive.and.any(coeff_values(ncoeff_fix+1:icycle) < zero)) then
           write(message, '(a)') ' Negative value detected...'
           gf_values(:,icoeff) = zero
           coeff_values = zero
         else
           call fit_polynomial_coeff_computeGF(coeff_values(1:icycle),energy_coeffs_tmp,&
&                                            fit_data%energy_diff,fcart_coeffs_tmp,fit_data%fcart_diff,&
&                                            gf_values(:,icoeff),list_coeffs_tmp(1:icycle),natom_sc,&
&                                            icycle,ncycle_max,ntime,strten_coeffs_tmp,&
&                                            fit_data%strten_diff,fit_data%training_set%sqomega)


           write (j_char, '(i7)') my_coeffindexes(icoeff)
           write(message, '(4x,a,3x,4ES18.10)') adjustl(j_char),&
&                                   gf_values(4,icoeff)* 1000*Ha_ev *factor,&
&                                   gf_values(1,icoeff)*HaBohr_meVAng**2,&
&                                   gf_values(2,icoeff)*HaBohr_meVAng**2,&
&                                   gf_values(3,icoeff)*HaBohr_meVAng**2
         end if
       else!In this case the matrix is singular
         gf_values(:,icoeff) = zero
         singular_coeffs(icoeff) = 1
         write(message, '(a)') ' The matrix is singular...'
       end if
       if(need_verbose) call wrtout(std_out,message,'COLL')
     end do

!    find the best coeff on each CPU
     mingf(:)  = 9D99
     index_min = 0
     do icoeff=1,my_ncoeff
       if(gf_values(1,icoeff) < zero) cycle
       if(abs(gf_values(1,icoeff)) <tol16) cycle
       if(gf_values(1,icoeff) < mingf(1) ) then
         mingf(:) = gf_values(:,icoeff)
         index_min = my_coeffindexes(icoeff)
       end if
     end do

!    MPI GATHER THE BEST COEFF ON EACH CPU
     if(nproc > 1)then
       buffGF(1,1) = index_min
       buffGF(2:5,1) =  mingf(:)

       call xmpi_allgatherv(buffGF,5,gf_mpi,buffsize,buffdisp, comm, ierr)
!      find the best coeff
       mingf(:)    = 9D99
       index_min= 0
       do icoeff=1,nproc
         if(gf_mpi(2,icoeff) < zero) cycle
         if(abs(gf_mpi(2,icoeff)) < tol16) cycle
         if(gf_mpi(2,icoeff) < mingf(1) ) then
           mingf(:) = gf_mpi(2:5,icoeff)
           index_min = int(gf_mpi(1,icoeff))
         end if
       end do
     end if

!    Check if there is still coefficient
     if(index_min==0) then
       exit
     else
       list_coeffs(icycle) = index_min
     end if

!    Check if this coeff is treat by this cpu and fill the
!    temporary array before broadcast
     rank_to_send = 0
     do icoeff=1,my_ncoeff


       if((my_coeffindexes(icoeff)==list_coeffs(icycle)))then

         if(need_initialize_data)then
           my_icoeff = icoeff
         else
!          Need to initialized the data for the fit for this coefficient
           my_icoeff = 1           
           call fit_polynomial_coeff_getFS(my_coeffs,fit_data%training_set%du_delta,&
&                                          fit_data%training_set%displacement,&
&                                          energy_coeffs,fcart_coeffs,natom_sc,eff_pot%crystal%natom,&
&                                          my_ncoeff,ntime,sc_size,fit_data%training_set%strain,&
&                                          strten_coeffs,fit_data%training_set%ucvol,&
&                                          my_coefflist(icoeff),1)
         end if

         energy_coeffs_tmp(icycle,:)    = energy_coeffs(my_icoeff,:)
         fcart_coeffs_tmp(:,:,icycle,:) = fcart_coeffs(:,:,my_icoeff,:)
         strten_coeffs_tmp(:,:,icycle)  = strten_coeffs(:,:,my_icoeff)
         call polynomial_coeff_free(coeffs_tmp(icycle))
         call polynomial_coeff_init(coeff_values(icycle),my_coeffs(icoeff)%nterm,&
&                                   coeffs_tmp(icycle),my_coeffs(icoeff)%terms,&
&                                   my_coeffs(icoeff)%name,&
&                                   check=.false.)

         rank_to_send = my_rank
         exit
       end if
     end do
!    Need to send the rank with the chosen coefficient
     call xmpi_sum(rank_to_send, comm, ierr)

!    Boadcast the coefficient
     call xmpi_bcast(energy_coeffs_tmp(icycle,:), rank_to_send, comm, ierr)
     call xmpi_bcast(fcart_coeffs_tmp(:,:,icycle,:) , rank_to_send, comm, ierr)
     call xmpi_bcast(strten_coeffs_tmp(:,:,icycle), rank_to_send, comm, ierr)
     call polynomial_coeff_broadcast(coeffs_tmp(icycle), rank_to_send, comm)

     if(need_verbose) then
       write(message, '(a,I0,2a)' )' Selecting the coefficient number ',list_coeffs(icycle),&
&                                   ' ===> ',trim(coeffs_tmp(icycle)%name)
       call wrtout(std_out,message,'COLL')

       write(message, '(2a,I0,a,ES24.16)' )' Standard deviation of the energy for',&
&                                        ' the iteration ',icycle_tmp,' (meV/atm): ',&
&                         mingf(4)* Ha_eV *1000 *factor
       call wrtout(std_out,message,'COLL')

       write (i_char, '(i3)') icycle
       write (j_char, '(i7)') list_coeffs(icycle)
       write(message, '(a,a,3x,a,3x,4ES18.10)') " ",adjustl(i_char),adjustl(j_char),&
&                                    mingf(4)* 1000*Ha_eV *factor,&
&                                    mingf(1)*HaBohr_meVAng**2,&
&                                    mingf(2)*HaBohr_meVAng**2,&
&                                    mingf(3)*HaBohr_meVAng**2
       call wrtout(ab_out,message,'COLL')
     end if

     ncycle_tot = ncycle_tot + 1

!    Check the stopping criterion
     converge = .false.
     if(tolMSDE  > zero)then
       if(abs(tolMSDE) > abs(mingf(4)*1000*Ha_eV *factor))then
         write(message,'(2a,ES18.10,a,ES18.10,a)') ch10," Fit process complete =>",&
&                                                mingf(4)*1000*Ha_eV * factor ," < ",tolMSDE,&
&                                              ' for MSDE'
         converge = .true.
       end if
     end if
     if(tolMSDF  > zero) then
       if(abs(tolMSDF) > abs(mingf(2)*HaBohr_meVAng**2))then
         write(message,'(2a,ES18.10,a,ES18.10,a)') ch10," Fit process complete =>",&
&                                                  mingf(2)*HaBohr_meVAng**2 ," < ",tolMSDF,&
&                                              ' for MSDF'
         converge = .true.
       end if
     end if
     if(tolMSDS  > zero) then
       if(abs(tolMSDS) > abs(mingf(3)*HaBohr_meVAng**2))then
         write(message,'(2a,ES18.10,a,ES18.10,a)') ch10," Fit process complete =>",&
&                                                  mingf(3)*HaBohr_meVAng**2 ," < ",tolMSDS,&
&                                              ' for MSDS'
         converge = .true.
       end if
     end if
     if(tolMSDFS > zero)then
       if(abs(tolMSDFS) > abs(mingf(1)*HaBohr_meVAng**2))then
         write(message,'(2a,ES18.10,a,ES18.10,a)') ch10," Fit process complete =>",&
&                                                  mingf(1)*HaBohr_meVAng**2 ," < ",tolMSDFS,&
&                                              ' for MSDFS'
         converge = .true.
       end if
     end if
     if(converge)then
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
       exit
     else
       if(any((/abs(tolMSDE),abs(tolMSDF),abs(tolMSDS),abs(tolMSDFS)/) > tol20) .and.&
&         icycle_tmp == ncycle)then
         write(message,'(2a,I0,a)') ch10," WARNING: ",ncycle,&
&                                   " cycles was not enougth to converge the fit process"
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')
       end if
     end if

   end do

 case(2)

!  Monte Carlo selection
   nsweep = 10000
!  If no coefficient imposed in the inputs we reset the goal function
   if (ncycle_tot == 0) then
     gf_values(:,:) = zero
     mingf(:) = 9D99
   else
     mingf = gf_values(:,1)
   end if
   call cpu_time(time)
   call ZBQLINI(int(time*1000000/(my_rank+1)))
   if(need_verbose)then
     write(message,'(a,I0,a)') " Start Monte Carlo simulations on ", nproc," CPU"
     if(nproc>1) write(message,'(2a)') trim(message)," (only print result of the master)"
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,2x,a,9x,a,14x,a,13x,a,14x,a)') ch10," Iteration ","MSDE","MSDFS","MSDF","MSdS"
     call wrtout(std_out,message,'COLL')
     write(message,'(a,5x,a,8x,a,8x,a,8x,a)') "              ","(meV/atm)","(eV^2/A^2)","(eV^2/A^2)",&
&                                            "(eV^2/A^2)"
     call wrtout(std_out,message,'COLL')

   end if

   do ii = 1,1!nyccle
     do isweep =1,nsweep
       write (j_char, '(i7)') isweep
!TEST_AM
       icycle_tmp = int(ZBQLU01(zero)*(ncycle+1-1))+1
       icycle_tmp = ncycle
       do icycle=1,icycle_tmp
         icoeff = int(ZBQLU01(zero)*(my_ncoeff))+1
!         icycle = int(ZBQLU01(zero)*(ncycle))+1
         list_coeffs_tmp2(icycle) = icoeff
         list_coeffs_tmp(icycle)= icycle
!        Fill the temporary arrays
         energy_coeffs_tmp(icycle,:)    = energy_coeffs(icoeff,:)
         fcart_coeffs_tmp(:,:,icycle,:) = fcart_coeffs(:,:,icoeff,:)
         strten_coeffs_tmp(:,:,icycle)  = strten_coeffs(:,:,icoeff)
       end do
!TEST_AM

!      call the fit process routine
!      This routine solves the linear system proposed by 
!      C.Escorihuela-Sayalero see PRB95,094115(2017) [[cite:Escorihuela-Sayalero2017]]
       call fit_polynomial_coeff_solve(coeff_values(1:icycle_tmp),fcart_coeffs_tmp,fit_data%fcart_diff,&
&                                      energy_coeffs_tmp,fit_data%energy_diff,info,&
&                                      list_coeffs_tmp(1:icycle_tmp),natom_sc,icycle_tmp,ncycle_max,&
&                                      ntime,strten_coeffs_tmp,fit_data%strten_diff,&
&                                      fit_data%training_set%sqomega)
       if(info==0)then
         call fit_polynomial_coeff_computeGF(coeff_values(1:icycle_tmp),energy_coeffs_tmp,&
&                                            fit_data%energy_diff,fcart_coeffs_tmp,fit_data%fcart_diff,&
&                                            gf_values(:,1),list_coeffs_tmp(1:icycle_tmp),natom_sc,&
&                                            icycle_tmp,ncycle_max,ntime,strten_coeffs_tmp,&
&                                            fit_data%strten_diff,fit_data%training_set%sqomega)

       else!In this case the matrix is singular
         gf_values(:,icoeff) = zero
         singular_coeffs(icoeff) = 1
       end if

       if(gf_values(1,1) > zero.and.abs(gf_values(1,1))>tol16.and.&
&         gf_values(1,1) < mingf(1) ) then
         mingf = gf_values(:,1)
         list_coeffs(1:icycle_tmp) = list_coeffs_tmp2(1:icycle_tmp)
         ncycle_tot = icycle_tmp

         write(message, '(4x,a,3x,4ES18.10)') adjustl(j_char),&
&                                   gf_values(4,1)* 1000*Ha_ev *factor,&
&                                   gf_values(1,1)*HaBohr_meVAng**2,&
&                                   gf_values(2,1)*HaBohr_meVAng**2,&
&                                   gf_values(3,1)*HaBohr_meVAng**2
         if(need_verbose) call wrtout(std_out,message,'COLL')
       else
         list_coeffs_tmp2(1:icycle_tmp) = list_coeffs(1:icycle_tmp)
       end if
     end do

     if(nproc > 1) then
!TEST_AM
       do iproc=1,ncycle_tot
         stat_coeff(list_coeffs(iproc)) =  stat_coeff(list_coeffs(iproc)) + 1
       end do
!TEST_AM

!    Find the best model on all the CPUs
       buffGF(1,1) = zero
       buffGF(2:5,1) =  mingf(:)
       call xmpi_allgatherv(buffGF,5,gf_mpi,buffsize,buffdisp, comm, ierr)
!      find the best coeff
       mingf(:) = 9D99
       index_min= 0
       do iproc=1,nproc
         if(gf_mpi(2,iproc) < zero) cycle
         if(abs(gf_mpi(2,iproc)) <tol16) cycle
         if(gf_mpi(2,iproc) < mingf(1) ) then
           mingf(:) = gf_mpi(2:5,iproc)
           index_min = int(gf_mpi(1,iproc))
           rank_to_send = iproc-1
         end if
       end do
       write(message, '(2a,I0)') ch10,' Best model found on the CPU: ', rank_to_send
       call wrtout(std_out,message,'COLL')
     end if
   end do

!TEST_AM
   call xmpi_sum(stat_coeff, comm, ierr)
   do ii=1,ncoeff_tot
     write(100,*) ii,stat_coeff(ii)
   end do
   close(100)
!TEST_AM

!  Transfert final model
   if(nproc>1)then
     call xmpi_bcast(ncycle_tot,rank_to_send,comm,ierr)
     call xmpi_bcast(list_coeffs(1:ncycle_tot),rank_to_send,comm,ierr)
   end if
   do ii=1,ncycle_tot
     icoeff = list_coeffs(ii)
     list_coeffs_tmp(ii) = ii
!    Fill the temporary arrays
     energy_coeffs_tmp(ii,:)    = energy_coeffs(icoeff,:)
     fcart_coeffs_tmp(:,:,ii,:) = fcart_coeffs(:,:,icoeff,:)
     strten_coeffs_tmp(:,:,ii)  = strten_coeffs(:,:,icoeff)
     call polynomial_coeff_free(coeffs_tmp(ii))
     call polynomial_coeff_init(one,my_coeffs(icoeff)%nterm,&
&                               coeffs_tmp(ii),my_coeffs(icoeff)%terms,&
&                               my_coeffs(icoeff)%name,&
&                               check=.false.)
   end do
 end select

!This routine solves the linear system proposed by 
! C.Escorihuela-Sayalero see PRB95,094115(2017) [[cite:Escorihuela-Sayalero2017]]
 if(ncycle_tot > 0)then

   call fit_polynomial_coeff_solve(coeff_values(1:ncycle_tot),fcart_coeffs_tmp,fit_data%fcart_diff,&
&                                  energy_coeffs_tmp,fit_data%energy_diff,info,&
&                                  list_coeffs_tmp(1:ncycle_tot),natom_sc,&
&                                  ncycle_tot,ncycle_max,ntime,strten_coeffs_tmp,&
&                                  fit_data%strten_diff,fit_data%training_set%sqomega)

   if(need_verbose) then
     write(message, '(3a)') ch10,' Fitted coefficients at the end of the fit process: '
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if
   do ii = 1,ncycle_tot
     if(list_coeffs(ii) ==0) cycle
!    Set the value of the coefficient
     coeffs_tmp(ii)%coefficient = coeff_values(ii)
     if(need_verbose) then
       write(message, '(a,I0,a,ES19.10,2a)') " ",list_coeffs(ii)," =>",coeff_values(ii),&
         &                                " ",trim(coeffs_tmp(ii)%name)
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if
   end do

   call fit_polynomial_coeff_computeGF(coeff_values(1:ncycle_tot),energy_coeffs_tmp,&
&                                      fit_data%energy_diff,fcart_coeffs_tmp,fit_data%fcart_diff,&
&                                      gf_values(:,1),list_coeffs_tmp(1:ncycle_tot),natom_sc,&
&                                      ncycle_tot,ncycle_max,ntime,strten_coeffs_tmp,&
&                                      fit_data%strten_diff,fit_data%training_set%sqomega)

   if(need_verbose) then
!  Print the standard deviation after the fit
     write(message,'(4a,ES24.16,4a,ES24.16,2a,ES24.16,2a,ES24.16,a)' )ch10,&
&                    ' Mean Standard Deviation values at the end of the fit process (meV/atm): ',ch10,& 
&                    '   Energy          : ',&
&               gf_values(4,1)*Ha_EV*1000*factor ,ch10,&
&                    ' Goal function values at the end of the fit process (eV^2/A^2):',ch10,&
&                    '   Forces+Stresses : ',&
&               gf_values(1,1)*(HaBohr_meVAng)**2,ch10,&
&                    '   Forces          : ',&
&               gf_values(2,1)*(HaBohr_meVAng)**2,ch10,&
&                    '   Stresses        : ',&
&               gf_values(3,1)*(HaBohr_meVAng)**2,ch10
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if
    

!  Set the final set of coefficients into the eff_pot type
   call effective_potential_setCoeffs(coeffs_tmp(1:ncycle_tot),eff_pot,ncycle_tot)

   ! If Wanted open the anharmonic_terms_file and write header
   filename = "TRS_fit_diff"
   ncoeff_model = eff_pot%anharmonics_terms%ncoeff
   if(need_prt_anh .and. ncoeff_model > 0 )then 
     call effective_potential_writeAnhHead(ncoeff_model,filename,&
&                                     eff_pot%anharmonics_terms) 
   else if (need_prt_anh)then
     write(message, '(6a,I3,3a)' )ch10,&
&          ' --- !WARNING',ch10,&
&          '     Printing of anharmonic terms has been asked,but',ch10,&
&          '     there are',ncoeff_model,'anharmonic terms in the potential',ch10,&
&          ' ---',ch10
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if 

! Calculate MSD values for final model 
   call fit_polynomial_coeff_computeMSD(eff_pot,hist,gf_values(4,1),gf_values(2,1),gf_values(1,1),&
&                                       natom_sc,ntime,fit_data%training_set%sqomega,&
&                                       compute_anharmonic=.TRUE.,print_file=.TRUE.,filename=filename)


   INQUIRE(FILE='TRS_fit_diff_anharmonic_terms_energy.dat',OPENED=file_opened,number=unit_anh)
   if(file_opened) close(unit_anh)

!    if(need_verbose) then
! !  Print the standard deviation after the fit
!      write(message,'(4a,ES24.16,4a,ES24.16,2a,ES24.16,2a,ES24.16,a)' )ch10,&
! &                    ' Mean Standard Deviation values at the end of the fit process (meV/f.u.):',&
! &               ch10,'   Energy          : ',&
! &               gf_values(4,1)*Ha_EV*1000/ ncell ,ch10,&
! &                    ' Goal function values at the end of the fit process (eV^2/A^2):',ch10,&
! &                    '   Forces+Stresses : ',&
! &               (gf_values(1,1)+gf_values(2,1))*(HaBohr_meVAng)**2,ch10,&
! &                    '   Forces          : ',&
! &               gf_values(2,1)*(HaBohr_meVAng)**2,ch10,&
! &                    '   Stresses        : ',&
! &               gf_values(3,1)*(HaBohr_meVAng)**2,ch10
!      call wrtout(ab_out,message,'COLL')
!      call wrtout(std_out,message,'COLL')
!    end if

 else
   if(need_verbose) then
     write(message, '(9a)' )ch10,&
&          ' --- !WARNING',ch10,&
&          '     The fit process does not provide possible terms.',ch10,&
&          '     Please make sure that the terms set is correct',ch10,&
&          ' ---',ch10
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if
 end if
!Deallocation of arrays
 call fit_data_free(fit_data)

!Deallocate the temporary coefficient
 do ii=1,ncycle_max
   call polynomial_coeff_free(coeffs_tmp(ii))
 end do
 ABI_DATATYPE_DEALLOCATE(coeffs_tmp)
 do ii=1,my_ncoeff
   call polynomial_coeff_free(my_coeffs(ii))
 end do

 ABI_DATATYPE_DEALLOCATE(my_coeffs)
 ABI_DEALLOCATE(buffsize)
 ABI_DEALLOCATE(buffdisp)
 ABI_DEALLOCATE(buffGF)
 ABI_DEALLOCATE(coeff_values)
 ABI_DEALLOCATE(energy_coeffs)
 ABI_DEALLOCATE(energy_coeffs_tmp)
 ABI_DEALLOCATE(fcart_coeffs)
 ABI_DEALLOCATE(fcart_coeffs_tmp)
 ABI_DEALLOCATE(gf_mpi)
 ABI_DEALLOCATE(gf_values)
 ABI_DEALLOCATE(list_coeffs)
 ABI_DEALLOCATE(list_coeffs_tmp)
 ABI_DEALLOCATE(list_coeffs_tmp2)
 ABI_DEALLOCATE(my_coeffindexes)
 ABI_DEALLOCATE(my_coefflist)
 ABI_DEALLOCATE(singular_coeffs)
 ABI_DEALLOCATE(strten_coeffs)
 ABI_DEALLOCATE(strten_coeffs_tmp)
 ABI_DEALLOCATE(stat_coeff)

end subroutine fit_polynomial_coeff_fit
!!***

!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_getPositive
!!
!! NAME
!! fit_polynomial_coeff_getPositive
!!
!! FUNCTION
!! This routine fit a list of possible model.
!! Return in the isPositive array:
!!   0 if the model ii does not contain possive coefficients
!!   1 if the model ii contain possive coefficients
!!
!! INPUTS
!! eff_pot<type(effective_potential)> = effective potential
!! hist<type(abihist)> = The history of the MD (or snapshot of DFT
!! coeff_values(nmodel,ncoeff) = values of the coefficients for each model
!! isPositive(nmodel) = see description below
!! list_coeff(nmodel,ncoeff) = list of the models
!! ncoeff = number of coeff per model
!! nfixcoeff = will not test the nfixcoeff first coeffcients
!! nmodel = number of model
!! comm = MPI communicator
!! verbose  = optional, flag for the verbose mode
!!
!! OUTPUT
!! eff_pot = effective potential datatype with new fitted coefficients
!!
!! PARENTS
!!      mover_effpot
!!
!! CHILDREN
!!      destroy_supercell,generelist,init_supercell,xred2xcart
!!
!! SOURCE

subroutine fit_polynomial_coeff_getPositive(eff_pot,hist,coeff_values,isPositive,list_coeff,ncoeff,&
&                                           nfixcoeff,nmodel,comm,verbose)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncoeff,nfixcoeff,nmodel,comm
!arrays
 integer,intent(in)  :: list_coeff(nmodel,ncoeff)
 integer,intent(out) :: isPositive(nmodel)
 real(dp),intent(out) :: coeff_values(nmodel,ncoeff)
 type(effective_potential_type),intent(inout) :: eff_pot
 type(abihist),intent(inout) :: hist
 logical,optional,intent(in) :: verbose
!Local variables-------------------------------
!scalar
 integer :: ierr,ii,info,imodel,my_nmodel,nmodel_alone
 integer :: master,my_rank,ncoeff_tot,natom_sc,ncell
 integer :: nproc,ntime
 logical :: iam_master,need_verbose
!arrays
 integer :: sc_size(3)
 integer,allocatable  :: list_coeffs(:),my_modelindexes(:),my_modellist(:)
 real(dp),allocatable :: energy_coeffs(:,:),fcart_coeffs(:,:,:,:), strten_coeffs(:,:,:)
 type(polynomial_coeff_type),allocatable :: coeffs_in(:)
 type(fit_data_type) :: fit_data
 character(len=500) :: message
! *************************************************************************

!MPI variables
 master = 0
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

!Initialisation of optional arguments
 need_verbose = .TRUE.
 if(present(verbose)) need_verbose = verbose

!Get the list of coefficients from the eff_pot
 if(eff_pot%anharmonics_terms%ncoeff > 0)then
!  Copy the initial coefficients array
   ncoeff_tot = eff_pot%anharmonics_terms%ncoeff
   ABI_DATATYPE_ALLOCATE(coeffs_in,(ncoeff_tot))
   do ii=1,ncoeff_tot
     call polynomial_coeff_init(eff_pot%anharmonics_terms%coefficients(ii)%coefficient,&
&                               eff_pot%anharmonics_terms%coefficients(ii)%nterm,&
&                               coeffs_in(ii),&
&                               eff_pot%anharmonics_terms%coefficients(ii)%terms,&
&                               eff_pot%anharmonics_terms%coefficients(ii)%name,&
&                               check=.false.)
   end do
 end if

!Reset the output (we free the memory)
 call effective_potential_freeCoeffs(eff_pot)

!if the number of atoms in reference supercell into effpot is not corret,
!wrt to the number of atom in the hist, we set map the hist and set the good
!supercell
 if (size(hist%xred,2) /= eff_pot%supercell%natom) then
   call effective_potential_file_mapHistToRef(eff_pot,hist,comm,verbose=need_verbose)
 end if

!Initialisation of constants
 natom_sc   = eff_pot%supercell%natom
 ncell      = eff_pot%supercell%ncells
 ntime      = hist%mxhist
 do ii = 1, 3
   sc_size(ii) = eff_pot%supercell%rlatt(ii,ii)
 end do

!Initialisation of arrays:
 ABI_ALLOCATE(list_coeffs,(ncoeff_tot))
 list_coeffs  = 0
 do ii = 1,ncoeff_tot
   list_coeffs(ii) = ii
 end do

!Get the decomposition for each coefficients of the forces and stresses for
!each atoms and each step  equations 11 & 12 of  PRB95,094115(2017) [[cite:Escorihuela-Sayalero2017]]
 if(need_verbose)then
   write(message, '(a)' ) ' Initialisation of the fit process...'
   call wrtout(std_out,message,'COLL')
 end if
!Before the fit, compute constants with fit_data_compute.
!Conpute the strain of each configuration.
!Compute the displacmeent of each configuration.
!Compute the variation of the displacement due to strain of each configuration.
!Compute fixed forces and stresse and get the standard deviation.
!Compute Sheppard and al Factors  \Omega^{2} see J.Chem Phys 136, 074103 (2012) [[cite:Sheppard2012]].
 call fit_data_compute(fit_data,eff_pot,hist,comm,verbose=need_verbose)

!Get the decomposition for each coefficients of the forces,stresses and energy for
!each atoms and each step  (see equations 11 & 12 of  
! PRB95,094115(2017)) [[cite:Escorihuela-Sayalero2017]] + allocation
 ABI_ALLOCATE(energy_coeffs,(ncoeff_tot,ntime))
 ABI_ALLOCATE(fcart_coeffs,(3,natom_sc,ncoeff_tot,ntime))
 ABI_ALLOCATE(strten_coeffs,(6,ntime,ncoeff_tot))

 call fit_polynomial_coeff_getFS(coeffs_in,fit_data%training_set%du_delta,&
&                                fit_data%training_set%displacement,&
&                                energy_coeffs,fcart_coeffs,natom_sc,eff_pot%crystal%natom,&
&                                ncoeff_tot,ntime,sc_size,&
&                                fit_data%training_set%strain,strten_coeffs,&
&                                fit_data%training_set%ucvol,list_coeffs,ncoeff_tot)


!set MPI, really basic stuff...
 nmodel_alone = mod(nmodel,nproc)
 my_nmodel = int(aint(real(nmodel,sp)/(nproc)))

 if(my_rank >= (nproc-nmodel_alone)) then
   my_nmodel = my_nmodel  + 1
 end if

 ABI_ALLOCATE(my_modelindexes,(my_nmodel))
 ABI_ALLOCATE(my_modellist,(my_nmodel))

!2:compute the number of model and the list of the corresponding for each CPU.
 do imodel=1,my_nmodel
   if(my_rank >= (nproc-nmodel_alone))then
     my_modelindexes(imodel)=(int(aint(real(nmodel,sp)/nproc)))*(my_rank)+&
&                              (my_rank - (nproc-nmodel_alone)) + imodel
     my_modellist(imodel) = imodel
   else
     my_modelindexes(imodel)=(my_nmodel)*(my_rank)  + imodel
     my_modellist(imodel) = imodel
  end if
 end do


!Start fit process
 isPositive   = 0
 coeff_values = zero
 do ii=1,my_nmodel
   imodel = my_modelindexes(ii)
   call fit_polynomial_coeff_solve(coeff_values(imodel,1:ncoeff),fcart_coeffs,fit_data%fcart_diff,&
&                                  energy_coeffs,fit_data%energy_diff,info,&
&                                  list_coeff(imodel,1:ncoeff),natom_sc,ncoeff,&
&                                  ncoeff_tot,ntime,strten_coeffs,fit_data%strten_diff,&
&                                  fit_data%training_set%sqomega)

   if(info==0)then

     if (any(coeff_values(imodel,nfixcoeff+1:ncoeff) < zero))then
!       coeff_values(imodel,:) = zero
       isPositive(imodel) = 0
     else
       isPositive(imodel) = 1
     end if
   end if
 end do

 call xmpi_sum(isPositive, comm, ierr)
 call xmpi_sum(coeff_values, comm, ierr)

!Deallocation of arrays
 do ii=1,ncoeff_tot
   call polynomial_coeff_free(coeffs_in(ii))
 end do
 call fit_data_free(fit_data)
 ABI_DATATYPE_DEALLOCATE(coeffs_in)
 ABI_DEALLOCATE(energy_coeffs)
 ABI_DEALLOCATE(fcart_coeffs)
 ABI_DEALLOCATE(list_coeffs)
 ABI_DEALLOCATE(my_modelindexes)
 ABI_DEALLOCATE(my_modellist)
 ABI_DEALLOCATE(strten_coeffs)

end subroutine fit_polynomial_coeff_getPositive
!!***

!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_getCoeffBound
!!
!! NAME
!! fit_polynomial_coeff_getCoeffBound
!!
!! FUNCTION
!! This routine fit a list of possible model.
!! Return in the isPositive array:
!!
!! INPUTS
!! NEED TO UPDATE
!! eff_pot<type(effective_potential)> = effective potential
!! hist<type(abihist)> = The history of the MD (or snapshot of DFT
!! comm = MPI communicator
!! verbose  = optional, flag for the verbose mode
!!
!! OUTPUT
!!
!!
!! PARENTS
!!
!! CHILDREN
!!      destroy_supercell,generelist,init_supercell,xred2xcart
!!
!! SOURCE

subroutine fit_polynomial_coeff_getCoeffBound(eff_pot,coeffs_out,hist,ncoeff_bound,comm,verbose)

 implicit none

!Arguments ------------------------------------
 !scalars
 integer,intent(in) :: comm
 integer,intent(out) :: ncoeff_bound
 logical,optional,intent(in) :: verbose
!arrays
 type(abihist),intent(inout) :: hist
 type(effective_potential_type),target,intent(inout) :: eff_pot
 type(polynomial_coeff_type),allocatable,intent(out) :: coeffs_out(:)
!Local variables-------------------------------
!scalar
 integer :: counter,icoeff,icoeff_bound,idisp,istrain,ii
 integer :: istart,iterm,ndisp,nstrain,nterm,ncoeff_model,ncoeff_in,ncoeff_max
 real(dp):: weight
 logical :: need_verbose
!arrays
 integer,allocatable :: atindx(:,:),cells(:,:,:),direction(:)
 integer,allocatable :: power_disps(:),power_strain(:),strain(:)
 type(polynomial_term_type),dimension(:),allocatable :: terms
 integer,allocatable :: odd_coeff(:),need_bound(:)
 type(polynomial_coeff_type),pointer :: coeffs_in(:)
 type(polynomial_coeff_type),allocatable :: coeffs_test(:)
 character(len=5),allocatable :: symbols(:)
 character(len=200):: name
 character(len=500) :: msg
! *************************************************************************


!set the inputs varaibles
 ncoeff_model =  eff_pot%anharmonics_terms%ncoeff
 coeffs_in => eff_pot%anharmonics_terms%coefficients

!Do check
 if(ncoeff_model == 0)then
   write(msg,'(a)')'ncoeff_model must be different to 0'
   MSG_BUG(msg)
 end if

!Map the hist in order to be consistent with the supercell into reference_effective_potential
 call effective_potential_file_mapHistToRef(eff_pot,hist,comm)

!Initialisation of optional arguments
 need_verbose = .TRUE.
 if(present(verbose)) need_verbose = verbose

 write(msg, '(a)' ) ' Detection of the unbound coefficients'
 if(need_verbose)call wrtout(std_out,msg,'COLL')

!Allocation
 ncoeff_max = 2 * ncoeff_model
 ABI_ALLOCATE(odd_coeff,(ncoeff_max))
 ABI_ALLOCATE(need_bound,(ncoeff_max))

 ABI_ALLOCATE(symbols,(eff_pot%crystal%natom))
 call symbols_crystal(eff_pot%crystal%natom,eff_pot%crystal%ntypat,eff_pot%crystal%npsp,&
&                     symbols,eff_pot%crystal%typat,eff_pot%crystal%znucl)


 ABI_DATATYPE_ALLOCATE(coeffs_test,(ncoeff_max))

 do icoeff=1,ncoeff_model
   call polynomial_coeff_init(coeffs_in(icoeff)%coefficient,coeffs_in(icoeff)%nterm,&
&                             coeffs_test(icoeff),coeffs_in(icoeff)%terms,&
&                             coeffs_in(icoeff)%name,check=.false.)
 end do

!array to know which coeff has to be bound
 need_bound(:) = 1
 counter = 0
 ncoeff_in = ncoeff_model

 do while(.not.all(need_bound == 0).and.counter<1)
!  Get the coefficients with odd coefficient
   odd_coeff = 0
   if(counter>0) then
     need_bound(1:ncoeff_in) = 0
     icoeff_bound = ncoeff_in
   else
     icoeff_bound = 1
   end if

   do icoeff=icoeff_bound,ncoeff_model
     if(any(mod(coeffs_in(icoeff)%terms(1)%power_disp(:),2)/=0))then
       odd_coeff(icoeff) = 1
     end if
     if(any(mod(coeffs_in(icoeff)%terms(1)%power_strain(:),2)/=0))then
       odd_coeff(icoeff) = 1
     end if
     if(odd_coeff(icoeff) == 0 .and. coeffs_in(icoeff)%coefficient > zero) then
       need_bound(icoeff) = 0
     else
        need_bound(icoeff) = 1
     end if
   end do
   if(need_verbose)then
     write(msg, '(a)' ) ' The following coefficients need to be bound:'
     call wrtout(std_out,msg,'COLL')
     do icoeff=1,ncoeff_model
       if(need_bound(icoeff) == 1)then
         write(msg, '(2a)' ) ' =>',trim(coeffs_in(icoeff)%name)
         call wrtout(std_out,msg,'COLL')
       end if
     end do
   end if


   icoeff_bound = ncoeff_in + 1
   if(counter==0)then
     istart = 1
   else
     istart = ncoeff_in
   end if

   ncoeff_bound = count(need_bound(istart:ncoeff_model)==1)

   do icoeff=istart,ncoeff_model
     if(need_bound(icoeff)==1)then

       nterm = coeffs_in(icoeff)%nterm
       ndisp = coeffs_in(icoeff)%terms(1)%ndisp
       nstrain = coeffs_in(icoeff)%terms(1)%nstrain

       ABI_ALLOCATE(terms,(nterm))
       ABI_ALLOCATE(atindx,(2,ndisp))
       ABI_ALLOCATE(cells,(3,2,ndisp))
       ABI_ALLOCATE(direction,(ndisp))
       ABI_ALLOCATE(power_disps,(ndisp))
       ABI_ALLOCATE(power_strain,(nstrain))
       ABI_ALLOCATE(strain,(nstrain))

       do iterm=1,coeffs_in(icoeff)%nterm
         atindx(:,:) = coeffs_in(icoeff)%terms(iterm)%atindx(:,:)
         cells(:,:,:) = coeffs_in(icoeff)%terms(iterm)%cell(:,:,:)
         direction(:) = coeffs_in(icoeff)%terms(iterm)%direction(:)
         power_strain(:) = coeffs_in(icoeff)%terms(iterm)%power_strain(:)
         power_disps(:) = coeffs_in(icoeff)%terms(iterm)%power_disp(:)
         strain(:) =  coeffs_in(icoeff)%terms(iterm)%strain(:)
         weight =  1
         do idisp=1,ndisp
           if(mod(power_disps(idisp),2) /= 0) then
             power_disps(idisp) = power_disps(idisp) + 1
           else
             power_disps(idisp) = power_disps(idisp) + 2
           end if
         end do
         do istrain=1,nstrain
           if(mod(power_strain(istrain),2) /= 0)then
             power_strain(istrain) = power_strain(istrain) + 1
           else
             if(power_strain(istrain) < 4 ) power_strain(istrain) = power_strain(istrain) + 2
           end if
         end do

         call polynomial_term_init(atindx,cells,direction,ndisp,nstrain,terms(iterm),&
&                                  power_disps,power_strain,strain,weight,check=.true.)
       end do

       name = ""
       call polynomial_coeff_init(one,nterm,coeffs_test(icoeff_bound),terms,name,check=.true.)
       call polynomial_coeff_getName(name,coeffs_test(icoeff_bound),symbols,recompute=.TRUE.)
       call polynomial_coeff_SetName(name,coeffs_test(icoeff_bound))

!      Deallocate the terms
       do iterm=1,nterm
         call polynomial_term_free(terms(iterm))
       end do
       ABI_DEALLOCATE(terms)
       ABI_DEALLOCATE(atindx)
       ABI_DEALLOCATE(cells)
       ABI_DEALLOCATE(direction)
       ABI_DEALLOCATE(power_disps)
       ABI_DEALLOCATE(power_strain)
       ABI_DEALLOCATE(strain)

       icoeff_bound = icoeff_bound  + 1

     end if
   end do


   if(counter==0)ncoeff_model = ncoeff_model + ncoeff_bound
!   call effective_potential_setCoeffs(coeffs_test,eff_pot,ncoeff_model)
   call fit_polynomial_coeff_fit(eff_pot,(/0/),(/0/),hist,0,(/0,0/),1,0,&
&             -1,1,comm,verbose=.true.,positive=.false.)

   coeffs_in => eff_pot%anharmonics_terms%coefficients

   counter = counter + 1
 end do

 ABI_ALLOCATE(coeffs_out,(ncoeff_bound))
 do ii=1,ncoeff_bound
   icoeff_bound = ncoeff_in + ii
   call polynomial_coeff_init(one,coeffs_test(icoeff_bound)%nterm,coeffs_out(ii),&
&                    coeffs_test(icoeff_bound)%terms,coeffs_test(icoeff_bound)%name,check=.true.)
 end do
!Deallocation
 do ii=ncoeff_model,ncoeff_max
   call polynomial_coeff_free(coeffs_test(ii))
 end do

!Deallocation
 do icoeff=1,ncoeff_max
   call polynomial_coeff_free(coeffs_test(icoeff))
 end do
 ABI_DATATYPE_DEALLOCATE(coeffs_test)
 ABI_DEALLOCATE(odd_coeff)
 ABI_DEALLOCATE(need_bound)
 ABI_DEALLOCATE(symbols)


end subroutine fit_polynomial_coeff_getCoeffBound
!!***


!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_solve
!!
!! NAME
!! fit_polynomial_coeff_solve
!!
!! FUNCTION
!! Build and the solve the system to get the values of the coefficients
!! This routine solves the linear system proposed by 
!! C.Escorihuela-Sayalero see PRB95,094115(2017) [[cite:Escorihuela-Sayalero2017]]
!!
!! INPUTS
!! fcart_coeffs(3,natom_sc,ncoeff_max,ntime) = List of the values of the contribution to the
!!                                             cartesian forces for all coefficients
!!                                             for each direction and each time
!! fcart_diff(3,natom,ntime) = Difference of cartesian forces between DFT calculation and
!!                             fixed part of the model (more often harmonic part)
!! energy_coeffs(ncoeff,ntime)   = value of the energy for each  coefficient (Ha)
!! energy_diff(ntime) = Difference of energ ybetween DFT calculation and fixed part
!!                             of the model (more often harmonic part)
!! list_coeffs(ncoeff_fit) = List with the index of the coefficients used for this model
!! natom = Number of atoms
!! ncoeff_fit = Number of coeff for the fit (dimension of the system)
!! ncoeff_max = Maximum number of coeff in the list
!! ntime = Number of time (number of snapshot, number of md step...)
!! strten_coeffs(6,ntime,ncoeff_max) = List of the values of the contribution to the stress tensor
!!                                      of  the coefficients for each direction,time
!! strten_diff(6,natom) = Difference of stress tensor between DFT calculation and
!!                        fixed part of the model (more often harmonic part)
!! sqomega(ntime) =  Sheppard and al Factors \Omega^{2} see J.Chem Phys 136, 074103 (2012) [[cite:Sheppard2012]]
!!
!! OUTPUT
!! coefficients(ncoeff_fit) = Values of the coefficients
!! info_out = 0:  successful exit
!!          < 0:  if INFO = -i, the i-th argument had an illegal value
!!          > 0:  if INFO = i, U(i,i) computed in DOUBLE PRECISION is
!!                exactly zero.  The factorization has been completed,
!!                but the factor U is exactly singular, so the solution
!!                could not be computed.  = 0:  successful exit
!!          information from the subroutine dsgesv in LAPACK
!!
!! PARENTS
!!      m_fit_polynomial_coeff
!!
!! CHILDREN
!!      destroy_supercell,generelist,init_supercell,xred2xcart
!!
!! SOURCE

subroutine fit_polynomial_coeff_solve(coefficients,fcart_coeffs,fcart_diff,energy_coeffs,energy_diff,&
&                                     info_out,list_coeffs,natom,ncoeff_fit,ncoeff_max,ntime,&
&                                     strten_coeffs,strten_diff,sqomega)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: natom,ncoeff_fit,ncoeff_max,ntime
 integer,intent(out) :: info_out
!arrays
 real(dp),intent(in) :: energy_coeffs(ncoeff_max,ntime)
 real(dp),intent(in) :: energy_diff(ntime)
 integer,intent(in)  :: list_coeffs(ncoeff_fit)
 real(dp),intent(in) :: fcart_coeffs(3,natom,ncoeff_max,ntime)
 real(dp),intent(in) :: fcart_diff(3,natom,ntime)
 real(dp),intent(in) :: strten_coeffs(6,ntime,ncoeff_max)
 real(dp),intent(in) :: strten_diff(6,ntime),sqomega(ntime)
 real(dp),intent(out):: coefficients(ncoeff_fit)
!Local variables-------------------------------
!scalar
 integer :: ia,itime,icoeff,jcoeff,icoeff_tmp,jcoeff_tmp,mu,LDA,LDB,LDX,LDAF,N,NRHS
 real(dp):: efact,ffact,sfact,ftmpA,stmpA,ftmpB,stmpB,etmpA,etmpB,fmu,fnu,smu,snu,emu,enu
 integer :: INFO,ITER
 real(dp):: RCOND
 real(dp):: fcart_coeffs_tmp(3,natom,ntime)
 real(dp),allocatable:: AF(:,:),BERR(:),FERR(:),WORK(:),C(:),R(:)
 integer,allocatable :: IPIV(:),IWORK(:),SWORK(:)
!arrays
 real(dp),allocatable :: A(:,:),B(:,:)
 character(len=1) :: FACT,EQUED,TRANS
! character(len=500) :: message
! *************************************************************************

!0-Set variables for the
 N    = ncoeff_fit; NRHS = 1; LDA  = ncoeff_fit; LDB  = ncoeff_fit; LDX  = ncoeff_fit
 LDAF = ncoeff_fit;  RCOND = zero; INFO  = 0; TRANS='N'; EQUED='N'; FACT='N'

!Set the factors
 ffact = one/(3*natom*ntime)
 sfact = one/(6*ntime)
 efact = one/(ntime)

!0-Allocation
 ABI_ALLOCATE(A,(LDA,N))
 ABI_ALLOCATE(B,(LDB,NRHS))
 ABI_ALLOCATE(AF,(LDAF,N))
 ABI_ALLOCATE(IPIV,(N))
 ABI_ALLOCATE(R,(N))
 ABI_ALLOCATE(C,(N))
 ABI_ALLOCATE(FERR,(NRHS))
 ABI_ALLOCATE(BERR,(NRHS))
 ABI_ALLOCATE(WORK,(4*N))
 ABI_ALLOCATE(IWORK,(N))
 ABI_ALLOCATE(SWORK,(N*(N+NRHS)))
 A=zero; B=zero;
 AF = zero; IPIV = 1;
 R = one; C = one;
 FERR = zero; BERR = zero
 IWORK = 0; WORK = 0

!1-Get forces and stresses from the model and fill A
!  Fill alsor B with the forces and stresses from
!  the DFT snapshot and the model
!  See equation 17 of PRB95 094115 (2017) [[cite:Escorihuela-Sayalero2017]]
 do icoeff=1,ncoeff_fit
   icoeff_tmp = list_coeffs(icoeff)
   fcart_coeffs_tmp(:,:,:) = fcart_coeffs(:,:,icoeff_tmp,:)
   ftmpA= zero; ftmpB = zero
   stmpA= zero; stmpB = zero
   etmpA= zero; etmpB = zero
!  loop over the configuration
   do itime=1,ntime
!    Fill energy
     emu = energy_coeffs(icoeff_tmp,itime)
     do jcoeff=1,ncoeff_fit
       jcoeff_tmp = list_coeffs(jcoeff)
       enu = energy_coeffs(jcoeff_tmp,itime)
!       etmpA =  emu*enu
!       A(icoeff,jcoeff) = A(icoeff,jcoeff) + efact*etmpA
     end do
     etmpB = etmpB + energy_diff(itime)*emu / (sqomega(itime)**3)
     etmpB = zero ! REMOVE THIS LINE TO TAKE INTO ACOUNT THE ENERGY     

!    Fill forces
     do ia=1,natom
       do mu=1,3
         fmu = fcart_coeffs_tmp(mu,ia,itime)
         do jcoeff=1,ncoeff_fit
           jcoeff_tmp = list_coeffs(jcoeff)
           fnu = fcart_coeffs(mu,ia,jcoeff_tmp,itime)
           ftmpA =  fmu*fnu
           A(icoeff,jcoeff) = A(icoeff,jcoeff) + ffact*ftmpA
         end do
         ftmpB = ftmpB + fcart_diff(mu,ia,itime)*fmu
       end do !End loop dir
     end do !End loop natom
!    Fill stresses
     do mu=1,6
       smu = strten_coeffs(mu,itime,icoeff_tmp)
       do jcoeff=1,ncoeff_fit
         jcoeff_tmp = list_coeffs(jcoeff)
         snu = strten_coeffs(mu,itime,jcoeff_tmp)
         stmpA =  sqomega(itime)*smu*snu
         A(icoeff,jcoeff) = A(icoeff,jcoeff) + sfact*stmpA
       end do
       stmpB = stmpB + sqomega(itime)*strten_diff(mu,itime)*smu
     end do !End loop stress dir
   end do ! End loop time
   B(icoeff,1) = B(icoeff,1) + ffact*ftmpB + sfact*stmpB + efact*etmpB
 end do ! End loop icoeff

!2-Solve Ax=B
!OLD VERSION..
! call dgesvx(FACT,TRANS,N,NRHS,A,LDA,AF,LDAF,IPIV,EQUED,R,C,B,LDB,coefficients,LDX,&
!             RCOND,FERR,BERR,WORK,IWORK,INFO)
!U is nonsingular
! if (INFO==N+1) then
!   coefficients = zero
! end if
 call DSGESV(N,NRHS,A,LDA,IPIV,B,LDB,coefficients,LDX,WORK,SWORK,ITER,INFO)

!other routine
! call dgesv(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
! coefficients = B(:,NRHS)
 !U is nonsingular
 if (INFO==N+2) then
   coefficients = zero
 end if

 if(any(abs(coefficients)>1.0E10))then
   INFO = 1
   coefficients = zero
 end if

 info_out = INFO

 ABI_DEALLOCATE(AF)
 ABI_DEALLOCATE(IPIV)
 ABI_DEALLOCATE(R)
 ABI_DEALLOCATE(C)
 ABI_DEALLOCATE(FERR)
 ABI_DEALLOCATE(BERR)
 ABI_DEALLOCATE(WORK)
 ABI_DEALLOCATE(IWORK)
 ABI_DEALLOCATE(SWORK)
 ABI_DEALLOCATE(A)
 ABI_DEALLOCATE(B)

end subroutine fit_polynomial_coeff_solve
!!***

!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_computeGF
!!
!! NAME
!! fit_polynomial_coeff_computeGF
!!
!! FUNCTION
!! Compute the values of the goal function (Mean squared error) for
!!   gf_value(1) = forces (Ha/Bohr)**2
!!   gf_value(2) = stresses (Ha/Bohr)**2
!!   gf_value(3) = stresses+forces (Ha/Bohr)**2
!!   gf_value(4) = energy (Ha)
!!
!! INPUTS
!! coefficients(ncoeff)          = type(polynomial_coeff_type)
!! energy_coeffs(ncoeff,ntime)   = value of the energy for each  coefficient (Ha)
!! energy_diff(ntime) = Difference of energ ybetween DFT calculation and fixed part
!!                             of the model (more often harmonic part)
!!                             fixed part of the model (more often harmonic part)
!! fcart_coeffs(ncoeff,3,natom,ntime) = value of the forces for each coefficient
!!                                      (-1 factor is taking into acount) (Ha/Bohr)
!! fcart_diff(3,natom,ntime) = Difference of cartesian forces between DFT calculation and
!!                             fixed part of the model (more often harmonic part)
!! list_coeffs(ncoeff_fit) = List with the indexes of the coefficients used for this model
!! natom = Number of atoms
!! ncoeff_fit = Number of coefficients fitted
!! ncoeff_max = Maximum number of coeff in the list
!! ntime = Number of time in the history
!! strten_coeffs(ncoeff,3,natom,ntime)= value of the stresses for each coefficient
!!                                      (1/ucvol factor is taking into acount) (Ha/Bohr^3)
!! strten_diff(6,natom) = Difference of stress tensor between DFT calculation and
!!                        fixed part of the model (more often harmonic part)
!! sqomega =  Sheppard and al Factors \Omega^{2} see J.Chem Phys 136, 074103 (2012) [[cite:Sheppard2012]]
!!
!! OUTPUT
!! gf_value(4) = Goal function
!!
!! PARENTS
!!      m_fit_polynomial_coeff
!!
!! CHILDREN
!!      destroy_supercell,generelist,init_supercell,xred2xcart
!!
!! SOURCE

subroutine fit_polynomial_coeff_computeGF(coefficients,energy_coeffs,energy_diff,&
&                                         fcart_coeffs,fcart_diff,gf_value,list_coeffs,&
&                                         natom,ncoeff_fit,ncoeff_max,ntime,strten_coeffs,&
&                                         strten_diff,sqomega)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: natom,ncoeff_fit,ncoeff_max,ntime
!arrays
 integer,intent(in)  :: list_coeffs(ncoeff_fit)
 real(dp),intent(in) :: energy_coeffs(ncoeff_max,ntime)
 real(dp),intent(in) :: energy_diff(ntime)
 real(dp),intent(in) :: fcart_coeffs(3,natom,ncoeff_max,ntime)
 real(dp),intent(in) :: fcart_diff(3,natom,ntime)
 real(dp),intent(in) :: strten_coeffs(6,ntime,ncoeff_max)
 real(dp),intent(in) :: strten_diff(6,ntime),sqomega(ntime)
 real(dp),intent(in) :: coefficients(ncoeff_fit)
 real(dp),intent(out) :: gf_value(4)
!Local variables-------------------------------
!scalar
 integer :: ia,icoeff,icoeff_tmp,itime,mu
 real(dp):: etmp,emu,fmu,ftmp,smu,stmp
 real(dp) :: ffact,sfact,efact
!arrays
! *************************************************************************

!1-Compute the value of the goal function
! see equation 9 of PRB 95 094115(2017) [[cite:Escorihuela-Sayalero2017]]
 gf_value = zero
 etmp     = zero
 ftmp     = zero
 stmp     = zero

!Compute factors
 ffact = one/(3*natom*ntime)
 sfact = one/(6*ntime)
 efact = one/(ntime)

! loop over the configuration
 do itime=1,ntime
! Fill energy
   emu = zero
   do icoeff=1,ncoeff_fit
     icoeff_tmp = list_coeffs(icoeff)
     emu = emu + coefficients(icoeff)*energy_coeffs(icoeff_tmp,itime)
   end do
!   uncomment the next line to be consistent with the definition of the goal function   
!   etmp = etmp + (energy_diff(itime)-emu)**2
   etmp = etmp + abs(energy_diff(itime)-emu)
!  Fill forces
   do ia=1,natom
     do mu=1,3
       fmu  = zero
       do icoeff=1,ncoeff_fit
         icoeff_tmp = list_coeffs(icoeff)
         fmu =  fmu + coefficients(icoeff)*fcart_coeffs(mu,ia,icoeff_tmp,itime)
       end do
       ftmp = ftmp + (fcart_diff(mu,ia,itime)-fmu)**2
     end do !End loop dir
   end do !End loop natom
   do mu=1,6
     smu = zero
     do icoeff=1,ncoeff_fit
       icoeff_tmp = list_coeffs(icoeff)
       smu = smu + coefficients(icoeff)*strten_coeffs(mu,itime,icoeff_tmp)
     end do
     stmp = stmp + sqomega(itime)*(strten_diff(mu,itime)-smu)**2
   end do !End loop stress dir
 end do ! End loop time

 gf_value(1)   =  ffact*ftmp + sfact*stmp !+ efact*etmp !Stresses + Forces
 gf_value(2)   =  ffact*ftmp ! only Forces
 gf_value(3)   =  sfact*stmp ! only Stresses
 gf_value(4)   =  efact*etmp !abs(Energy)

end subroutine fit_polynomial_coeff_computeGF
!!***


!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_getFS
!!
!! NAME
!! fit_polynomial_coeff_getFS
!!
!! FUNCTION
!! Compute all the matrix elements of eq.11 and 12 in PRB95,094115 (2017) [[cite:Escorihuela-Sayalero2017]]
!!
!! INPUTS
!! coefficients(ncoeff)          = type(polynomial_coeff_type)
!! du_delta(6,3,natom_sc,ntime)  = Variation to displacements wrt to the strain (Bohr)
!! displacement(3,natom_sc,ntime)= Atomic displacement wrt to the reference (Bohr)
!! natom_sc = Number of atoms in the supercell
!! natom_uc = Number of atoms in the unit cell
!! ncoeff = Number of coefficients
!! ntime = Number of time in the history
!! sc_size(3) = Size of the supercell
!! strain(6,ntime) = Strain
!! ucvol(ntime) = Volume of the supercell for each time (Bohr^3)
!! cells(ncell) = Indexes of the cell treat by this CPU
!! ncell = Number of cell treat by this CPU
!! index_cells(ncell,3) = Indexes of the cells (1 1 1, 0 0 0 for instance) treat by this CPU
!! comm  = MPI communicator
!!
!! OUTPUT
!! fcart_out(ncoeff,3,natom,ntime) = value of the forces for each coefficient
!!                                   (-1 factor is taking into acount) (Ha/Bohr)
!! strten_out(ncoeff,3,natom,ntime)= value of the stresses for each coefficient
!!                                   (-1/ucvol factor is taking into acount) (Ha/Bohr^3)
!! energy_out(ncoeff,ntime)        = value of the energy for each  coefficient (Ha)
!!
!! PARENTS
!!      m_fit_polynomial_coeff
!!
!! CHILDREN
!!      destroy_supercell,generelist,init_supercell,xred2xcart
!!
!! SOURCE

subroutine fit_polynomial_coeff_getFS(coefficients,du_delta,displacement,energy_out,fcart_out,&
&                                     natom_sc,natom_uc,ncoeff_max,ntime,sc_size,strain,strten_out,&
&                                     ucvol,coeffs,ncoeff)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom_sc,natom_uc,ncoeff_max,ntime
 integer,intent(in) :: ncoeff
!arrays
 integer,intent(in) :: sc_size(3)
 integer,intent(in) :: coeffs(ncoeff_max)
 real(dp),intent(in) :: du_delta(6,3,natom_sc,ntime)
 real(dp),intent(in) :: displacement(3,natom_sc,ntime)
 real(dp),intent(in) :: strain(6,ntime),ucvol(ntime)
 real(dp),intent(out):: energy_out(ncoeff,ntime)
 real(dp),intent(out) :: fcart_out(3,natom_sc,ncoeff,ntime)
 real(dp),intent(out) :: strten_out(6,ntime,ncoeff)
 type(polynomial_coeff_type), intent(in) :: coefficients(ncoeff_max)
!Local variables-------------------------------
!scalar
 integer :: i1,i2,i3,ia1,ia2,ib1,ib2,ii,icell,icoeff,icoeff_tmp
 integer :: idir1,idir2,idisp1,idisp2,idisp1_strain,idisp2_strain
 integer :: iterm,itime,ndisp,ndisp_tot,nstrain,power_disp,power_strain
 real(dp):: disp1,disp2,tmp1,tmp2,tmp3,weight
!arrays
 integer :: cell_atoma1(3),cell_atoma2(3)
 integer :: cell_atomb1(3),cell_atomb2(3)

! *************************************************************************


!1-Get forces and stresses from the model
!  Initialisation of variables
 fcart_out(:,:,:,:) = zero
 strten_out(:,:,:)  = zero
 energy_out(:,:)    = zero

 icell = 0; ib1=0; ia1=0
 do i1=1,sc_size(1)
   do i2=1,sc_size(2)
     do i3=1,sc_size(3)
       ii = icell*natom_uc
       icell = icell + 1
!      Loop over configurations
       do itime=1,ntime
!       Loop over coefficients
         do icoeff_tmp=1,ncoeff
           icoeff = coeffs(icoeff_tmp)
!          Loop over terms of this coefficient
           do iterm=1,coefficients(icoeff)%nterm
             ndisp = coefficients(icoeff)%terms(iterm)%ndisp
             nstrain = coefficients(icoeff)%terms(iterm)%nstrain
             ndisp_tot = ndisp + nstrain
!            Set the weight of this term
             weight =coefficients(icoeff)%terms(iterm)%weight
             tmp1 = one
!            Loop over displacement and strain
             do idisp1=1,ndisp_tot

!              Set to one the acculation of forces and strain
               tmp2 = one
               tmp3 = one
!              Strain case idir => -6, -5, -4, -3, -2 or -1
               if (idisp1 > ndisp)then
                 idisp1_strain = idisp1 - ndisp
                 power_strain = coefficients(icoeff)%terms(iterm)%power_strain(idisp1_strain)
!                Get the direction of the displacement or strain
                 idir1 = coefficients(icoeff)%terms(iterm)%strain(idisp1_strain)
                 if(abs(strain(idir1,itime)) > tol10)then
!                  Accumulate energy fo each displacement (\sum ((A_x-O_x)^Y(A_y-O_c)^Z))
                   tmp1 = tmp1 * (strain(idir1,itime))**power_strain
                   if(power_strain > 1) then
!                    Accumulate stress for each strain (\sum (Y(eta_2)^Y-1(eta_2)^Z+...))
                     tmp3 = tmp3 *  power_strain*(strain(idir1,itime))**(power_strain-1)
                   end if
                 else
                   tmp1 = zero
                   if(power_strain > 1) then
                     tmp3 = zero
                   end if
                 end if
               else
!                Set the power_disp of the displacement:
                 power_disp = coefficients(icoeff)%terms(iterm)%power_disp(idisp1)
!                Get the direction of the displacement or strain
                 idir1 = coefficients(icoeff)%terms(iterm)%direction(idisp1)
!                Displacement case idir = 1, 2  or 3
!                indexes of the cell of the atom a
                 cell_atoma1 = coefficients(icoeff)%terms(iterm)%cell(:,1,idisp1)
                 if(cell_atoma1(1)/=0.or.cell_atoma1(2)/=0.or.cell_atoma1(3)/=0) then
!                  if the cell is not 0 0 0 we apply PBC:
                   cell_atoma1(1) =  i1 + cell_atoma1(1)
                   cell_atoma1(2) =  i2 + cell_atoma1(2)
                   cell_atoma1(3) =  i3 + cell_atoma1(3)
                   call getPBCIndexes_supercell(cell_atoma1(1:3),sc_size(1:3))
!                  index of the first atom (position in the supercell if the cell is not 0 0 0)
                   ia1 = (cell_atoma1(1)-1)*sc_size(2)*sc_size(3)*natom_uc+&
&                        (cell_atoma1(2)-1)*sc_size(3)*natom_uc+&
&                        (cell_atoma1(3)-1)*natom_uc+&
&                        coefficients(icoeff)%terms(iterm)%atindx(1,idisp1)
                 else
!                  index of the first atom (position in the supercell if the cell is 0 0 0)
                   ia1 = ii + coefficients(icoeff)%terms(iterm)%atindx(1,idisp1)
                 end if

!                indexes of the cell of the atom b  (with PBC) same as ia1
                 cell_atomb1 = coefficients(icoeff)%terms(iterm)%cell(:,2,idisp1)
                 if(cell_atomb1(1)/=0.or.cell_atomb1(2)/=0.or.cell_atomb1(3)/=0) then
                   cell_atomb1(1) =  i1 + cell_atomb1(1)
                   cell_atomb1(2) =  i2 + cell_atomb1(2)
                   cell_atomb1(3) =  i3 + cell_atomb1(3)
                   call getPBCIndexes_supercell(cell_atomb1(1:3),sc_size(1:3))

!                  index of the second atom in the (position in the supercell  if the cell is not 0 0 0)
                   ib1 = (cell_atomb1(1)-1)*sc_size(2)*sc_size(3)*natom_uc+&
&                        (cell_atomb1(2)-1)*sc_size(3)*natom_uc+&
&                        (cell_atomb1(3)-1)*natom_uc+&
&                        coefficients(icoeff)%terms(iterm)%atindx(2,idisp1)
                 else
!                  index of the first atom (position in the supercell if the cell is 0 0 0)
                   ib1 = ii + coefficients(icoeff)%terms(iterm)%atindx(2,idisp1)
                 end if

!                Get the displacement for the both atoms
                 disp1 = displacement(idir1,ia1,itime)
                 disp2 = displacement(idir1,ib1,itime)

                 if(abs(disp1) > tol10 .or. abs(disp2)> tol10)then
!                  Accumulate energy fo each displacement (\sum ((A_x-O_x)^Y(A_y-O_c)^Z))
                   tmp1 = tmp1 * (disp1-disp2)**power_disp
                   if(power_disp > 1) then
!                    Accumulate forces for each displacement (\sum (Y(A_x-O_x)^Y-1(A_y-O_c)^Z+...))
                     tmp2 = tmp2 * power_disp*(disp1-disp2)**(power_disp-1)
                   end if
                 else
                   tmp1 = zero
                   if(power_disp > 1) then
                     tmp2 = zero
                   end if
                 end if
               end if

               do idisp2=1,ndisp_tot
                 if(idisp2 /= idisp1) then

!                  Strain case
                   if (idisp2 > ndisp)then
                     idisp2_strain = idisp2 - ndisp
                     idir2 = coefficients(icoeff)%terms(iterm)%strain(idisp2_strain)
!                    Set the power_strain of the strain:
                     power_strain = coefficients(icoeff)%terms(iterm)%power_strain(idisp2_strain)
!                    Accumulate energy forces
                     tmp2 = tmp2 * (strain(idir2,itime))**power_strain
!                    Accumulate stress for each strain (\sum (Y(eta_2)^Y-1(eta_2)^Z+...))
                     tmp3 = tmp3 * (strain(idir2,itime))**power_strain
!                  Atomic displacement case
                   else
!                    Set the power_disp of the displacement:
                     power_disp = coefficients(icoeff)%terms(iterm)%power_disp(idisp2)
!                    Set the direction of the displacement:
                     idir2 = coefficients(icoeff)%terms(iterm)%direction(idisp2)

                     cell_atoma2=coefficients(icoeff)%terms(iterm)%cell(:,1,idisp2)
                     if(cell_atoma2(1)/=0.or.cell_atoma2(2)/=0.or.cell_atoma2(3)/=0) then
                       cell_atoma2(1) =  i1 + cell_atoma2(1)
                       cell_atoma2(2) =  i2 + cell_atoma2(2)
                       cell_atoma2(3) =  i3 + cell_atoma2(3)
                       call getPBCIndexes_supercell(cell_atoma2(1:3),sc_size(1:3))
!                      index of the first atom (position in the supercell and direction)
!                      if the cell of the atom a is not 0 0 0 (may happen)
                       ia2 = (cell_atoma2(1)-1)*sc_size(2)*sc_size(3)*natom_uc+&
&                            (cell_atoma2(2)-1)*sc_size(3)*natom_uc+&
&                            (cell_atoma2(3)-1)*natom_uc+&
&                        coefficients(icoeff)%terms(iterm)%atindx(1,idisp2)
                     else
!                      index of the first atom (position in the supercell and direction)
                       ia2 = ii + coefficients(icoeff)%terms(iterm)%atindx(1,idisp2)
                     end if

                     cell_atomb2 = coefficients(icoeff)%terms(iterm)%cell(:,2,idisp2)

                     if(cell_atomb2(1)/=0.or.cell_atomb2(2)/=0.or.cell_atomb2(3)/=0) then
!                      indexes of the cell2 (with PBC)
                       cell_atomb2(1) =  i1 + cell_atomb2(1)
                       cell_atomb2(2) =  i2 + cell_atomb2(2)
                       cell_atomb2(3) =  i3 + cell_atomb2(3)
                       call getPBCIndexes_supercell(cell_atomb2(1:3),sc_size(1:3))

!                      index of the second atom in the (position in the supercell)
                       ib2 = (cell_atomb2(1)-1)*sc_size(2)*sc_size(3)*natom_uc+&
&                            (cell_atomb2(2)-1)*sc_size(3)*natom_uc+&
&                            (cell_atomb2(3)-1)*natom_uc+&
&                            coefficients(icoeff)%terms(iterm)%atindx(2,idisp2)
                     else
                       ib2 = ii + coefficients(icoeff)%terms(iterm)%atindx(2,idisp2)
                     end if

                     disp1 = displacement(idir2,ia2,itime)
                     disp2 = displacement(idir2,ib2,itime)

                     tmp2 = tmp2 * (disp1-disp2)**power_disp
                     tmp3 = tmp3 * (disp1-disp2)**power_disp

                   end if
                 end if
               end do

               if(idisp1 > ndisp)then
!                Accumule stress tensor
                 strten_out(idir1,itime,icoeff_tmp) = strten_out(idir1,itime,icoeff_tmp) + &
&                                                      weight * tmp3 / ucvol(itime)
               else
!                Accumule  forces
                 fcart_out(idir1,ia1,icoeff_tmp,itime)=fcart_out(idir1,ia1,icoeff_tmp,itime)+weight*tmp2
                 fcart_out(idir1,ib1,icoeff_tmp,itime)=fcart_out(idir1,ib1,icoeff_tmp,itime)-weight*tmp2
               end if
             end do

!            accumule energy
             energy_out(icoeff_tmp,itime) = energy_out(icoeff_tmp,itime) +  weight * tmp1

           end do!End do iterm
         end do!End do coeff
       end do!End time
     end do!End do i3
   end do!End do i2
 end do!End do i1

!ADD variation of the atomic displacement due to the strain
 do icoeff=1,ncoeff
   do itime=1,ntime
     do ia1=1,natom_sc
       do idir1=1,3
         do idir2=1,6
           strten_out(idir2,itime,icoeff) = strten_out(idir2,itime,icoeff) + &
&                     du_delta(idir2,idir1,ia1,itime)*fcart_out(idir1,ia1,icoeff,itime)/ucvol(itime)
         end do
       end do
     end do
   end do
 end do

! multiply by -1
 fcart_out(:,:,:,:) = -1 * fcart_out(:,:,:,:)

end subroutine fit_polynomial_coeff_getFS
!!***


!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_computeMSD
!!
!! NAME
!! fit_polynomial_coeff_computeMSD
!!
!! FUNCTION
!! Compute the Mean square error of the energy, forces and stresses
!!
!! INPUTS
!! eff_pot<type(effective_potential)> = effective potential
!! hist<type(abihist)> = The history of the MD
!! natom = number of atom
!! ntime = number of time in the hist
!! sqomega =  Sheppard and al Factors \Omega^{2} see J.Chem Phys 136, 074103 (2012) [[cite:Sheppard2012]]
!! compute_anharmonic = TRUE if the anharmonic part of the effective potential
!!                           has to be taking into acount
!! print_file = if True, a ASCII file with the difference in energy will be print
!!
!! OUTPUT
!! mse  =  Mean square error of the energy   (Hatree)
!! msef =  Mean square error of the forces   (Hatree/Bohr)**2
!! mses =  Mean square error of the stresses (Hatree/Bohr)**2
!!
!! PARENTS
!!      m_fit_polynomial_coeff
!!
!! CHILDREN
!!      destroy_supercell,generelist,init_supercell,xred2xcart
!!
!! SOURCE

subroutine fit_polynomial_coeff_computeMSD(eff_pot,hist,mse,msef,mses,natom,ntime,sqomega,&
&                                          compute_anharmonic,print_file,filename,scup_dtset)

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,ntime
 real(dp),intent(out):: mse,msef,mses
 logical,optional,intent(in) :: compute_anharmonic,print_file
!arrays
 real(dp) :: sqomega(ntime)
 type(effective_potential_type),intent(in) :: eff_pot
 type(abihist),intent(in) :: hist
!Strings/Characters
 character(len=fnlen),optional,intent(in) :: filename
 type(scup_dtset_type),optional,intent(inout) :: scup_dtset
!Local variables-------------------------------
!scalar
integer :: ii,ia,mu,unit_energy,unit_stress,unit_anh,ifirst,itime
! integer :: ifirst
 real(dp):: energy,energy_harm
 logical :: need_anharmonic = .TRUE.,need_print=.FALSE., anh_opened,need_elec_eval
 !arrays
 real(dp):: fcart(3,natom),fred(3,natom),strten(6),rprimd(3,3),xred(3,natom)
!Strings/Characters 
 character(len=fnlen) :: file_energy, file_stress, file_anh, name_file
 character(len=500) :: msg
 type(abihist) :: hist_out
 character(len=200) :: filename_hist

! *************************************************************************

 !Do some checks
 if(ntime /= hist%mxhist)then
   write(msg,'(a)')'ntime is not correct'
   MSG_BUG(msg)
 end if

 if(natom /= size(hist%xred,2)) then
   write(msg,'(a)')'natom is not correct'
   MSG_BUG(msg)
 end if

 if(present(compute_anharmonic))then
   need_anharmonic = compute_anharmonic
 end if

 name_file=''
 if(present(filename))name_file = filename
 
 need_print=.FALSE. 
 if(present(print_file))need_print=print_file
 
 need_elec_eval = .FALSE. 
 if(present(scup_dtset))need_elec_eval=scup_dtset%scup_elec_model


 if(need_print .and. present(filename))then
   !MS hist out uncommented for PHONOPY test
   !call abihist_init(hist_out,natom,ntime,.false.,.false.)
   file_energy=trim(name_file)//'_energy.dat'
   unit_energy = get_unit()
   if (open_file(file_energy,msg,unit=unit_energy,form="formatted",&
&     status="unknown",action="write") /= 0) then
     MSG_ERROR(msg)
   end if
   unit_stress = get_unit()
   file_stress=trim(name_file)//'_stress.dat'
   if (open_file(file_stress,msg,unit=unit_stress,form="formatted",&
&     status="unknown",action="write") /= 0) then
     MSG_ERROR(msg)
   end if
 else if(need_print .and. .not. present(filename))then 
   write(msg,'(3a)')' You asked for printing of the MSD-values',ch10,& 
&        ' without specifying a filename'
   MSG_ERROR(msg) 
 end if
 
 file_anh=trim(name_file)//'_anharmonic_terms_energy.dat'
 anh_opened=.FALSE.
 INQUIRE(FILE=file_anh,OPENED=anh_opened,number=unit_anh)

 mse  = zero
 msef = zero
 mses = zero
 do ii=1,ntime ! Loop over configurations
   xred(:,:)   = hist%xred(:,:,ii)
   rprimd(:,:) = hist%rprimd(:,:,ii)
   if(anh_opened .eqv. .TRUE.)then
     write(unit_anh,'(I7)',advance='no') ii !If wanted Write cycle to anharmonic_energy_contribution file
   end if
#if defined DEV_MS_SCALEUP 
   !Pass print options to scale-up
   itime = ii 
   if(need_elec_eval)then
        call global_set_parent_iter(itime)
        ! Set all print options to false. 
        call global_set_print_parameters(geom=.FALSE.,eigvals=.FALSE.,eltic=.FALSE.,&
&                orbocc=.FALSE.,bands=.FALSE.)
        if(ii == 1 .or. modulo(ii,scup_dtset%scup_printniter) == 0)then 
           call global_set_print_parameters(scup_dtset%scup_printgeom,scup_dtset%scup_printeigv,scup_dtset%scup_printeltic,& 
&                   scup_dtset%scup_printorbocc,scup_dtset%scup_printbands)
        end if 
   end if 
#endif
   call effective_potential_evaluate(eff_pot,energy_harm,fcart,fred,strten,natom,rprimd,&
&                                    xred=xred,compute_anharmonic=.False.,verbose=.false.,&
&                                    elec_eval=need_elec_eval)

   call effective_potential_evaluate(eff_pot,energy,fcart,fred,strten,natom,rprimd,&
&                                    xred=xred,compute_anharmonic=need_anharmonic,verbose=.false.,&
&                                    filename=file_anh,elec_eval=need_elec_eval)

   if(need_print)then
     WRITE(unit_energy ,'(I10,7(F23.14))') ii,hist%etot(ii),energy_harm,energy,&
&                                       abs(hist%etot(ii) - energy_harm),abs(hist%etot(ii) - energy)
     WRITE(unit_stress,'(I10,12(F23.14))') ii,hist%strten(:,ii),strten(:)
   end if
    
    !MS Uncommented for abihist test 
    !ifirst=merge(0,1,(ii>1))
    !filename_hist = trim("test.nc")
    !hist_out%fcart(:,:,hist_out%ihist) = fcart(:,:)
    !hist_out%strten(:,hist_out%ihist)  = strten(:)
    !hist_out%etot(hist_out%ihist)      = energy
    !hist_out%entropy(hist_out%ihist)   = hist%entropy(ii)
    !hist_out%time(hist_out%ihist)      = real(ii,kind=dp)
    !call vel2hist(ab_mover%amass,hist,vel,vel_cell)
    !call var2hist(hist%acell(:,ii),hist_out,natom,hist%rprimd(:,:,ii),hist%xred(:,:,ii),.false.)
    !call write_md_hist(hist_out,filename_hist,ifirst,ii,natom,1,eff_pot%crystal%ntypat,&
 !&                    eff_pot%supercell%typat,eff_pot%crystal%amu,eff_pot%crystal%znucl,&
 !&                    real(100,dp),(/real(100,dp),real(100,dp)/))

   mse  = mse  + abs(hist%etot(ii) - energy)
   do ia=1,natom ! Loop over atoms
     do mu=1,3   ! Loop over cartesian directions
       msef = msef + (hist%fcart(mu,ia,ii)  - fcart(mu,ia))**2
     end do
   end do
   do mu=1,6 ! Loop over stresses
     mses = mses + sqomega(ii)*(hist%strten(mu,ii) - strten(mu))**2
   end do
   end do ! End loop itime 

 mse  = mse  /  ntime
 msef = msef / (3*natom*ntime)
 mses = mses / (6*ntime)

 if(need_print)then
   close(unit_energy)
   close(unit_stress)
 end if

 !MS uncommented for PHONOPY TEST
 !call abihist_free(hist_out)

end subroutine fit_polynomial_coeff_computeMSD
!!***

!MARCUS_EXPERIMENTAL_SECTION 
!!****f* m_fit_polynomiaL_coeff/testEffPot
!! NAME
!!  testEffPot
!!
!! FUNCTION
!!  Calculate the energy, forces for displacements provided 
!!  in an test-set (input:hist) within a given effective potential 
!!  (input: eff_pot)
!!  If the test set is from DFT and contains DFT energies and forces 
!!  calculate the Goal Function values and the MSD of the Energy with 
!!  respect to the DFT energies 
!!
!! INPUTS
!! eff_pot = effective_potential datatype
!! hist = abihist datatype
!!
!! OUTPUT
!!
!! SOURCE

subroutine fit_polynomial_coeff_testEffPot(eff_pot,hist,master,comm,print_anharmonic,scup_dtset)

       
  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in) :: master,comm
!logicals
  logical,optional,intent(in) :: print_anharmonic
!array
  type(effective_potential_type),intent(inout) :: eff_pot
  type(abihist),intent(in) :: hist
  type(scup_dtset_type),optional,intent(inout) :: scup_dtset
!Local variables-------------------------------
!reals 
  real(dp) :: factor,mse,msef,mses
  type(fit_data_type) :: test_data
  real(dp),allocatable :: sqomega(:),ucvol(:)
  real(dp),parameter :: HaBohr_meVAng = 27.21138386 / 0.529177249
!scalar
  integer :: itime, test,unit_anh
  integer :: natom,ntime,ncoeff,my_rank
!logicals 
  logical :: iam_master, need_print_anharmonic,file_opened
!strings/characters
 character(len=fnlen) :: filename 
 character(len=1000) :: message
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
! *************************************************************************
  
  !MPI variables
  my_rank=xmpi_comm_rank(comm)
  iam_master = (my_rank == master)

  !Initialisation of optional arguments
  need_print_anharmonic = .FALSE. 
  if(present(print_anharmonic)) need_print_anharmonic = print_anharmonic


  !Setting/Allocating other Variables 
  natom = size(hist%xred,2)   
  factor   = 1._dp/natom
  ntime = hist%mxhist 
  ABI_ALLOCATE(sqomega,(ntime))
  ABI_ALLOCATE(ucvol,(ntime))
  sqomega = zero 
  filename = 'TES_fit_diff'
  ncoeff = eff_pot%anharmonics_terms%ncoeff
   
  do itime=1,ntime 
!  Compute \Omega^{2} and ucvol for each time
   call metric(gmet,gprimd,-1,rmet,hist%rprimd(:,:,itime),ucvol(itime))
!  Formula: sqomega(itime) = (((ucvol(itime)**(-2.))* ((natom)**(0.5)))**(-1.0/3.0))**2
!   Compact form:
   sqomega(itime) = ((ucvol(itime)**(4.0/3.0)) / ((natom)**(1/3.0)))
  end do 

       
  if(need_print_anharmonic) call effective_potential_writeAnhHead(ncoeff,&
&                            filename,eff_pot%anharmonics_terms)                  

  call fit_polynomial_coeff_computeMSD(eff_pot,hist,mse,msef,mses,natom,ntime,&
&                                      sqomega,&
&                 compute_anharmonic=.TRUE.,print_file=.TRUE.,filename=filename,scup_dtset=scup_dtset)


!  Print the standard deviation after the fit
     write(message,'(6a,ES24.16,6a,ES24.16,2a,ES24.16,2a,ES24.16,a)' )ch10,&
&                    ' Mean Standard Deviation values of the effective-potential',ch10,&
&                    ' with respect to the test-set (meV/atm):',&
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


  !Deallocating
  ABI_DEALLOCATE(sqomega)
  ABI_DEALLOCATE(ucvol)

  INQUIRE(FILE='TES_fit_diff_anharmonic_terms_energy.dat',OPENED=file_opened,number=unit_anh)
  if(file_opened) close(unit_anh)


end subroutine fit_polynomial_coeff_testEffPot
!!***

!!      m_fit_polynomial_coeff,multibinit
!!      generelist,polynomial_coeff_free,polynomial_coeff_getname
!!      polynomial_coeff_init,polynomial_term_free,polynomial_term_init,wrtout

!!****f* m_fit_polynomial_coeff/fit_polynomial_printSystemFiles
!!
!! NAME
!! fit_polynomial_printSystemFiles
!!
!! FUNCTION
!! Print the files for the fitting script
!!
!! INPUTS
!! eff_pot<type(effective_potential)> = effective potential
!! hist<type(abihist)> = datatype with the  history of the MD
!!
!! OUTPUT
!!
!! PARENTS
!!      multibinit
!!
!! CHILDREN
!!      destroy_supercell,generelist,init_supercell,xred2xcart
!!
!! SOURCE

subroutine fit_polynomial_printSystemFiles(eff_pot,hist)

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(effective_potential_type), intent(in) :: eff_pot
 type(abihist),intent(in) :: hist
!Local variables-------------------------------
!scalar
 integer :: ia,ib,ib1,ii,jj,irpt,kk,ll,mu,nu,nstep,nshift
 integer :: natom_uc
 integer :: unit_born=22,unit_epsiloninf=23,unit_md=24
 integer :: unit_harmonic=25,unit_ref=26,unit_strain=27,unit_sym=28
!arrays
 integer,allocatable :: typat_order(:),typat_order_uc(:)
 integer, dimension(3)  :: A,ncell
 real(dp), allocatable :: xcart(:,:),fcart(:,:)
 character(len=500) :: msg
 type(supercell_type) :: supercell
! *************************************************************************

!Create new supercell corresponding to the MD
 ncell = (/2,2,2/)
 call init_supercell(eff_pot%crystal%natom, (/ncell(1),0,0,  0,ncell(2),0,  0,0,ncell(3)/),&
&                    eff_pot%crystal%rprimd,eff_pot%crystal%typat,&
&                    eff_pot%crystal%xcart,eff_pot%crystal%znucl, supercell)

!allocation of array
 ABI_ALLOCATE(xcart,(3,supercell%natom))
 ABI_ALLOCATE(fcart,(3,supercell%natom))
 ABI_ALLOCATE(typat_order,(supercell%natom))
 ABI_ALLOCATE(typat_order_uc,(eff_pot%crystal%natom))

 A = (/ 2, 3, 1/)

 nshift = product(ncell)
 natom_uc = eff_pot%crystal%natom
!Fill the typat_order array:
!In the fit script the atom must be in the order 11111 222222 33333 ..
!and the order of the atom can not be change in the fit script,
!we transform into the format of the script
 ib = 1
 ib1= 1
 do ii=1,eff_pot%crystal%ntypat
   jj = A(ii)
   do kk=1,natom_uc
     if(supercell%typat(kk)==jj)then
       typat_order_uc(ib1) = kk
       ib1 = ib1 + 1
       do ll=1,nshift
         ia = (ll-1)*natom_uc + kk
         typat_order(ib) = ia
         ib = ib + 1
       end do
     end if
   end do
 end do

! BORN CHARGES FILE
 if (open_file('system/Born_Charges',msg,unit=unit_born,form="formatted",&
&    status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 do ii=1,eff_pot%crystal%ntypat
   jj = A(ii)
   do ia=1,eff_pot%crystal%natom
     if(eff_pot%crystal%typat(ia)==jj)then
       write(unit_born,'(i2,a,1F10.5)') ia,"    ",eff_pot%crystal%amu(eff_pot%crystal%typat(ia))
       do mu=1,3
         WRITE(unit_born,'(a,3(F23.14))') "     ",eff_pot%harmonics_terms%zeff(:,mu,ia)
       end do
     end if
   end do
 end do

!DIELECTRIC TENSOR FILE
 if (open_file('system/Dielectric_Tensor',msg,unit=unit_epsiloninf,form="formatted",&
&    status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 do mu=1,3
   WRITE(unit_epsiloninf,'(3(F23.14))') eff_pot%harmonics_terms%epsilon_inf(:,mu)
 end do


!REFERENCE STRUCTURE FILE
 if (open_file('system/Reference_structure',msg,unit=unit_ref,form="formatted",&
&    status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 write(unit_ref,'("Energy (Hartree)")')
 write(unit_ref,'("================")')
 write(unit_ref,'(F23.14)') (hist%etot(1)/nshift)
 write(unit_ref,'("")')
 write(unit_ref,'("Cell vectors")')
 write(unit_ref,'("============")')
 do jj=1,3
   write(unit_ref,'(3(F22.14))') (supercell%rprimd(:,jj))
 end do

 write(unit_ref,'("")')
 write(unit_ref,'("Atomic positions (Bohr radius)")')
 write(unit_ref,'("==============================")')

 do ia=1,supercell%natom
   write(unit_ref,'(3(F23.14))') supercell%xcart(:,typat_order(ia))
 end do

!Harmonic XML file
 if (open_file('system/harmonic.xml',msg,unit=unit_harmonic,form="formatted",&
&     status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if

!Write header
 write(unit_harmonic,'("<?xml version=""1.0"" ?>")')
 write(unit_harmonic,'("<name>")')

 do irpt=1,eff_pot%harmonics_terms%ifcs%nrpt
   if(any(abs(eff_pot%harmonics_terms%ifcs%short_atmfrc(:,:,:,:,irpt))>tol9)) then
     write(unit_harmonic,'("  <local_force_constant units=""hartree/bohrradius**2"">")')
     write(unit_harmonic,'("    <data>")')
     do ia=1,eff_pot%crystal%natom
       do mu=1,3
         do ib=1,eff_pot%crystal%natom
           do  nu=1,3
             write(unit_harmonic,'(F22.14)', advance="no")&
&                 (eff_pot%harmonics_terms%ifcs%short_atmfrc(mu,typat_order_uc(ia),&
&                                                              nu,typat_order_uc(ib),irpt))
           end do
         end do
         write(unit_harmonic,'(a)')''
       end do
     end do
     write(unit_harmonic,'("    </data>")')
     write(unit_harmonic,'("    <cell>")')
     write(unit_harmonic,'(3(I4))') (eff_pot%harmonics_terms%ifcs%cell(:,irpt))
     write(unit_harmonic,'("    </cell>")')
     write(unit_harmonic,'("  </local_force_constant>")')
   end if
 end do
 write(unit_harmonic,'("</name>")')

!STRAIN FILE
 if (open_file('system/Strain_Tensor',msg,unit=unit_strain,form="formatted",&
&     status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write(unit_strain,'(6(F23.14))') (eff_pot%harmonics_terms%elastic_constants)

! SYM FILE
 if (open_file('system/symmetry_operations',msg,unit=unit_sym,form="formatted",&
&     status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write(unit_sym,'("(x,y,z)  (y,-x,z) (z,x,y) (y,z,x) (x,z,y) (y,x,z) (z,y,x) (x,-y,-z) (z,-x,-y)",&
&                " (y,-z,-x) (x,-z,-y) (y,-x,-z) (z,-y,-x) (-x,y,-z) (-z,x,-y) (-y,z,-x) (-x,z,-y)",&
&                " (-y,x,-z) (-z,y,-x) (-x,-y,z) (-z,-x,y) (-y,-z,x) (-x,-z,y) (-y,-x,z) (-z,-y,x)",&
&                " (-x,-y,-z) (-z,-x,-y) (-y,-z,-x) (-x,-z,-y) (-y,-x,-z) (-z,-y,-x) (-x,y,z)",&
&                " (-z,x,y) (-y,z,x) (-x,z,y) (-y,x,z) (-z,y,x) (x,-y,z) (z,-x,y) (y,-z,x) (x,-z,y)",&
&                " (z,-y,x) (x,y,-z) (z,x,-y) (y,z,-x) (x,z,-y) (y,x,-z) (z,y,-x)")')


!MD file
 nstep = hist%mxhist
 if (open_file('system/Molecular_dynamic',msg,unit=unit_md,form="formatted",&
&     status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 do ii=1,nstep
   write(unit_md,'(I5)') ii-1
   write(unit_md,'(F22.14)') hist%etot(ii)/nshift
   do jj=1,3
     write(unit_md,'(3(F22.14))') (hist%rprimd(:,jj,ii))
   end do
!  Set xcart and fcart for this step
   call xred2xcart(supercell%natom,hist%rprimd(:,:,ii),&
&                  xcart,hist%xred(:,:,ii))

   fcart(:,:) = hist%fcart(:,:,ii)

   do ia=1,supercell%natom
     write(unit_md,'(3(E22.14),3(E22.14))') xcart(:,typat_order(ia)),fcart(:,typat_order(ia))
   end do
   write(unit_md,'(6(E22.14))') hist%strten(:,ii)
 end do

!Close files
 close(unit_ref)
 close(unit_born)
 close(unit_harmonic)
 close(unit_epsiloninf)
 close(unit_md)
 close(unit_strain)
 close(unit_sym)

!Deallocation array
 ABI_DEALLOCATE(typat_order)
 ABI_DEALLOCATE(typat_order_uc)
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(fcart)
 call destroy_supercell(supercell)

end subroutine fit_polynomial_printSystemFiles
!!***

recursive subroutine genereList(i,m,m_max,n_max,list,list_out,size,compute)

 implicit none

!Arguments ---------------------------------------------
!scalar
 integer, intent(in) :: m_max,n_max,m,size
 integer, intent(inout) :: i
 logical,intent(in) :: compute
!arrays
 integer, intent(out) :: list(m_max),list_out(size,m_max)
!Local variables ---------------------------------------
!scalar
 integer n
!arrays

! *************************************************************************
 if (m > m_max) then
   i = i + 1
   if(compute)list_out(i,:) = list(:)
 else
   do n = 1, n_max
     if (m == 1)then
       list(m) = n
       call genereList (i, m + 1,m_max,n_max,list,list_out,size,compute)
     else if (n > list(m - 1)) then
       list(m) = n
       call genereList (i, m + 1,m_max,n_max,list,list_out,size,compute)
     end if
   end do
 end if

end subroutine genereList
!!***

end module m_fit_polynomial_coeff
!!***
