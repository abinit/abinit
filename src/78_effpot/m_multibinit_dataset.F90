!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_multibinit_dataset
!! NAME
!!  m_multibinit_dataset
!!
!! FUNCTION
!!  module with the type for the input of multibinit (should be clean)
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2019 ABINIT group (AM)
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

module m_multibinit_dataset

 use defs_basis
 use m_abicore
 use m_errors

 use m_parser, only : intagm
 use m_ddb,    only : DDB_QTOL

 implicit none

 private

 public :: multibinit_dtset_type
 public :: multibinit_dtset_init
 public :: multibinit_dtset_free
 public :: outvars_multibinit
 public :: invars10
!!***

!----------------------------------------------------------------------

!!****t* m_multibinit_dataset/multibinit_dtset_type
!! NAME
!! multibinit_dtset_type
!!
!! FUNCTION
!! The multibinit_dtset_type structured datatype
!! gather all the input variables for the multibinit code.
!!
!! SOURCE

 type multibinit_dtset_type
    
! Integer
  integer :: asr
  integer :: analyze_anh_pot
  integer :: brav
  integer :: chneut
  integer :: confinement
  integer :: conf_power_disp
  integer :: conf_power_strain
  integer :: dipdip
  integer :: eivec
  integer :: elphflag
  integer :: enunit
  integer :: bound_model
  integer :: bound_maxCoeff
  integer :: bound_SPCoupling
  integer :: bound_AnhaStrain
  integer :: bound_step
  integer :: fit_anhaStrain
  integer :: fit_SPCoupling
  integer :: fit_generateCoeff
  integer :: fit_initializeData
  integer :: fit_coeff
  integer :: fit_option
  integer :: fit_ncoeff
  integer :: fit_nbancoeff
  integer :: fit_nfixcoeff
  integer :: ts_option
  integer :: ifcana
  integer :: ifcflag
  integer :: ifcout
  ! TODO hexu: why integer dtion?
  integer :: dtion
  integer :: dynamics
  integer :: natifc
  integer :: natom
  integer :: ncoeff
  integer :: nctime
  integer :: ntime
  integer :: nnos
  integer :: nph1l
  integer :: nph2l
  integer :: nqshft
  integer :: nsphere
  integer :: optcell
  integer :: prt_model
  integer :: prt_names
  integer :: dipdip_prt
  integer :: prt_phfrq
  integer :: prt_ifc
  integer :: strcpling  ! Print the 3rd order in xml file
  integer :: prtsrlr  ! print the short-range/long-range decomposition of phonon freq.
  integer :: rfmeth
  integer :: restartxf
  integer :: symdynmat
  integer :: dipdip_range(3)
  integer :: fit_grid(3)
  integer :: fit_rangePower(2)
  integer :: bound_rangePower(2)
  integer :: bound_cell(3)
  integer :: ncell(3)
  integer :: ngqpt(9)             ! ngqpt(9) instead of ngqpt(3) is needed in wght9.f
  integer :: ng2qpt(3)
  integer :: kptrlatt(3,3)
  integer :: kptrlatt_fine(3,3)
  integer :: qrefine(3)


  ! TODO hexu: add parameters for spin.
  integer :: spin_calc_traj_obs
  integer :: spin_calc_thermo_obs
  integer :: spin_calc_correlation_obs
  integer :: spin_dipdip
  integer :: spin_dynamics
  integer :: spin_init_state
  integer :: spin_nctime
  integer :: spin_ntime_pre
  integer :: spin_ntime
  integer :: spin_nmatom !TODO hexu: is it needed?
  integer :: spin_n1l
  integer :: spin_n2l
  integer :: spin_sia_add
  integer :: spin_temperature_nstep    ! var temperature number of steps
  integer :: spin_var_temperature
  integer :: spin_write_traj

  ! parameters for spin-lattice coupling
  integer :: slc_coupling

! Real(dp)
  real(dp) :: bmass
  real(dp) :: conf_power_fact_disp
  real(dp) :: conf_power_fact_strain
  real(dp) :: delta_df
  real(dp) :: energy_reference
  real(dp) :: bound_cutoff
  real(dp) :: bound_Temp
  real(dp) :: fit_cutoff
  real(dp) :: fit_tolMSDF
  real(dp) :: fit_tolMSDS
  real(dp) :: fit_tolMSDE
  real(dp) :: fit_tolMSDFS
  real(dp) :: temperature
  real(dp) :: rifcsph
  real(dp) :: conf
  real(dp) :: acell(3)
  real(dp) :: strten_reference(6)
  real(dp) :: strtarget(6)
  real(dp) :: conf_cutoff_strain(6)
  real(dp) :: rprim(3,3)

  ! lattice (new) related 
  real(dp) :: latt_friction ! langevin dynamics friction
  real(dp) :: latt_taut     ! Berendsen taut
  real(dp) :: latt_taup     ! 
  real(dp) :: latt_compressibility
  integer :: latt_mask(3)

  ! TODO hexu:add parameters for spin
  real(dp) :: spin_dt
  real(dp) :: spin_damping
  real(dp) :: spin_sia_k1amp
  real(dp) :: spin_temperature
  ! TODO hexu: add spin convergence tol. (or remove it)
  real(dp) :: spin_temperature_start   ! var temperature start
  real(dp) :: spin_temperature_end     ! var temperature end
  real(dp) :: spin_tolavg !average
  real(dp) :: spin_tolvar !covariance

  real(dp) :: spin_mag_field(3)  ! external magnetic field
  real(dp) :: spin_qpoint(3)
  real(dp) :: spin_sia_k1dir(3)
  real(dp) :: spin_init_qpoint(3) ! qpoint to specify initial spin configuration
  real(dp) :: spin_init_rotate_axis(3) ! rotation axis to specify initial spin configuration  
  real(dp) :: spin_init_orientation(3) ! spin orientation in primitive cell which is then rotated

! Integer arrays
  integer, allocatable :: atifc(:)
  ! atifc(natom)
  integer, allocatable :: fit_fixcoeff(:)
  ! fit_fixcoeffs(fit_nfixcoeff)

  integer, allocatable :: fit_bancoeff(:)
  ! fit_bancoeffs(fit_nbancoeff)

  !integer, allocatable :: spin_sublattice(:) ! TODO hexu: difficult to use, better in xml?

  real(dp), allocatable :: qmass(:)
  ! qmass(nnos)


! Real arrays
  real(dp), allocatable :: coefficients(:)
  ! coefficients(ncoeff)

  real(dp), allocatable :: conf_cutoff_disp(:)
  ! conf_cuttoff(natom)
  
  real(dp),allocatable  :: q1shft(:,:)
  !q1shft(3,nqshft)  SHIFT for Q point

  real(dp), allocatable :: qnrml1(:)
  ! qnrml1(nph1l)

  real(dp), allocatable :: qnrml2(:)
  ! qnrml1(nph1l)

  real(dp), allocatable :: qph1l(:,:)
  ! qph1l(3,nph1l)

  real(dp), allocatable :: qph2l(:,:)
  ! qph2l(3,nph2l)

  ! spin part
  !real(dp), allocatable :: gilbert_damping ! if not provided in xml or override is needed. 
  !real(dp), allocatable :: gyro_ratio(:) ! if not provided in xml

  !real(dp), allocatable :: qspin1l(:,:)
  !real(dp), allocatable :: qspin2l(:,:)

 end type multibinit_dtset_type
!!***

contains
!!***

!!****f* m_multibinit_dataset/multibinit_dtset_init
!!
!! NAME
!! multibinit_dtset_init
!!
!! FUNCTION
!! Init the dtset datatype
!!
!! INPUTS
!! natom=number of atoms, needed for atifc
!!
!! OUTPUT
!! multibinit_dtset <type(multibinit_dtset_type)> = datatype with all the input variables
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine multibinit_dtset_init(multibinit_dtset,natom)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom
 type(multibinit_dtset_type),intent(inout) :: multibinit_dtset
!Local variables -------------------------
!scalars
!arrays

!*********************************************************************

!copy natom to multibinit_dtset
 multibinit_dtset%natom=natom

!=====================================================================
!Scalars
!=====================================================================
 multibinit_dtset%asr=2
 multibinit_dtset%analyze_anh_pot=0 
 multibinit_dtset%brav=1
 multibinit_dtset%bmass=0
 multibinit_dtset%chneut=0
 multibinit_dtset%confinement=0
 multibinit_dtset%conf_power_disp=0
 multibinit_dtset%conf_power_strain=0
 multibinit_dtset%conf_power_fact_disp=100
 multibinit_dtset%conf_power_fact_strain=100
 multibinit_dtset%delta_df= 1d-02
 multibinit_dtset%dipdip=1
 multibinit_dtset%dipdip_prt=0
 multibinit_dtset%dtion=100
 multibinit_dtset%dynamics=0
 multibinit_dtset%eivec=0
 multibinit_dtset%energy_reference= zero
 multibinit_dtset%enunit=0
 multibinit_dtset%fit_anhaStrain=0
 multibinit_dtset%bound_model=0
 multibinit_dtset%bound_anhaStrain=0
 multibinit_dtset%bound_cutoff=0 
 multibinit_dtset%bound_maxCoeff=4
 multibinit_dtset%bound_temp=325
 multibinit_dtset%bound_step=1000
 multibinit_dtset%bound_SPCoupling=1
 multibinit_dtset%fit_coeff=0
 multibinit_dtset%fit_cutoff=0
 multibinit_dtset%fit_nbancoeff=0
 multibinit_dtset%fit_ncoeff=0
 multibinit_dtset%ts_option=0
 multibinit_dtset%fit_nfixcoeff=0
 multibinit_dtset%fit_option=0
 multibinit_dtset%fit_SPCoupling=1
 multibinit_dtset%fit_generateCoeff=1
 multibinit_dtset%fit_initializeData=1
 multibinit_dtset%fit_tolMSDE=zero
 multibinit_dtset%fit_tolMSDS=zero
 multibinit_dtset%fit_tolMSDF=zero
 multibinit_dtset%fit_tolMSDFS=zero
 multibinit_dtset%ifcana=0
 multibinit_dtset%ifcflag=1
 multibinit_dtset%ifcout=-1
 multibinit_dtset%prtsrlr=0
 ! Langevin friction
 multibinit_dtset%latt_friction=1d-4
 ! Berendsen taut
 multibinit_dtset%latt_taut=1000.0
 multibinit_dtset%latt_taup=1000.0
 multibinit_dtset%latt_compressibility=0.0

 multibinit_dtset%ntime=200
 multibinit_dtset%nctime=1
 multibinit_dtset%natifc=natom
 multibinit_dtset%ncoeff=0
 multibinit_dtset%nph1l=1
 multibinit_dtset%nph2l=0
 multibinit_dtset%nqshft=1
 multibinit_dtset%nnos=0
 multibinit_dtset%nsphere=0
 multibinit_dtset%optcell=0
 multibinit_dtset%prt_model=0
 multibinit_dtset%prt_names=0
 multibinit_dtset%prt_phfrq=0
 multibinit_dtset%prt_ifc = 0
 multibinit_dtset%strcpling = -1
 multibinit_dtset%qrefine=1
 multibinit_dtset%restartxf=0
 multibinit_dtset%rfmeth=1
 multibinit_dtset%rifcsph=zero
 multibinit_dtset%symdynmat=1
 multibinit_dtset%temperature=325

 multibinit_dtset%spin_calc_traj_obs=0
 multibinit_dtset%spin_calc_thermo_obs=1
 multibinit_dtset%spin_calc_correlation_obs=0
 multibinit_dtset%spin_dipdip=0
 multibinit_dtset%spin_dynamics=0
 multibinit_dtset%spin_init_state=1
 multibinit_dtset%spin_ntime_pre=0
 multibinit_dtset%spin_ntime=10000
 multibinit_dtset%spin_nctime=100
 multibinit_dtset%spin_nmatom=0
 multibinit_dtset%spin_n1l=1
 multibinit_dtset%spin_n2l=0

 multibinit_dtset%spin_dt=100

 multibinit_dtset%spin_damping=-1.0
 multibinit_dtset%spin_sia_add=0
 multibinit_dtset%spin_sia_k1amp=zero
multibinit_dtset%spin_temperature=325
multibinit_dtset%spin_temperature_start=0.0
multibinit_dtset%spin_temperature_end= 0.0
multibinit_dtset%spin_temperature_nstep= 0
multibinit_dtset%spin_tolavg=1d-2 ! TODO hexu: to be decided. should it be a function of temperature?
multibinit_dtset%spin_tolvar=1d-3 ! TODO hexu: as above. 

multibinit_dtset%spin_var_temperature=0 
multibinit_dtset%spin_write_traj=1 

multibinit_dtset%slc_coupling=0

!=======================================================================
!Arrays
!=======================================================================
 multibinit_dtset%acell(:) = one
 multibinit_dtset%conf_cutoff_strain(1:6) = zero
 multibinit_dtset%dipdip_range(:)= (/0,0,0/)
 multibinit_dtset%fit_grid(:)= 1
 multibinit_dtset%fit_rangePower(:)= (/3,4/)
 multibinit_dtset%bound_rangePower(:)= (/6,6/)
 multibinit_dtset%bound_cell(:)= (/6,6,6/)
 multibinit_dtset%ncell(:)= 0
 multibinit_dtset%ngqpt(:) = 0
 multibinit_dtset%ng2qpt(:)= 0
 multibinit_dtset%strtarget(1:6) = zero
 multibinit_dtset%qmass(:)= zero
 multibinit_dtset%rprim(:,:)= zero
 multibinit_dtset%strten_reference(:)= zero

 multibinit_dtset%spin_mag_field(:)=zero
 multibinit_dtset%spin_qpoint(:)=zero
 multibinit_dtset%spin_init_qpoint(:)=zero
 multibinit_dtset%spin_init_rotate_axis(:)=(/1.0, 0.0, 0.0/)
 multibinit_dtset%spin_init_orientation(:)=(/0.0, 0.0, 1.0/)
 
 multibinit_dtset%spin_sia_k1dir(:)=(/0.0,0.0,1.0/)


 ABI_ALLOCATE(multibinit_dtset%atifc,(natom))
 multibinit_dtset%atifc(:)=0
 ABI_ALLOCATE(multibinit_dtset%conf_cutoff_disp,(multibinit_dtset%natom))
 multibinit_dtset%conf_cutoff_disp(:)=zero
 ABI_ALLOCATE(multibinit_dtset%q1shft,(3,multibinit_dtset%nqshft))
 multibinit_dtset%q1shft(:,:) = zero

 multibinit_dtset%latt_mask(:) = 0

end subroutine multibinit_dtset_init
!!***

!!****f* m_multibinit_dataset/multibinit_dtset_free
!!
!! NAME
!!  multibinit_dtset_free
!!
!! FUNCTION
!!  deallocate remaining arrays in the multibinit_dtset datastructure
!!
!! INPUTS
!!  multibinit_dtset <type(multibinit_dtset_type)> = multibinit_dataset structure
!!
!! OUTPUTS
!!  multibinit_dtset <type(multibinit_dtset_type)> = multibinit_dataset structure
!!
!! PARENTS
!!      multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine multibinit_dtset_free(multibinit_dtset)

!Arguments ------------------------------------
!scalars
 type(multibinit_dtset_type), intent(inout) :: multibinit_dtset

! *************************************************************************

 if (allocated(multibinit_dtset%atifc))  then
   ABI_DEALLOCATE(multibinit_dtset%atifc)
 end if
 if (allocated(multibinit_dtset%conf_cutoff_disp))  then
   ABI_DEALLOCATE(multibinit_dtset%conf_cutoff_disp)
 end if
 if (allocated(multibinit_dtset%fit_fixcoeff))  then
   ABI_DEALLOCATE(multibinit_dtset%fit_fixcoeff)
 end if
  if (allocated(multibinit_dtset%fit_bancoeff))  then
   ABI_DEALLOCATE(multibinit_dtset%fit_bancoeff)
 end if
 if (allocated(multibinit_dtset%qmass))  then
   ABI_DEALLOCATE(multibinit_dtset%qmass)
 end if
 if (allocated(multibinit_dtset%coefficients))  then
   ABI_DEALLOCATE(multibinit_dtset%coefficients)
 end if
 if (allocated(multibinit_dtset%qnrml1))  then
   ABI_DEALLOCATE(multibinit_dtset%qnrml1)
 end if
 if (allocated(multibinit_dtset%qnrml2))  then
   ABI_DEALLOCATE(multibinit_dtset%qnrml2)
 end if
 if (allocated(multibinit_dtset%qph1l))  then
   ABI_DEALLOCATE(multibinit_dtset%qph1l)
 end if
 if (allocated(multibinit_dtset%qph2l))  then
   ABI_DEALLOCATE(multibinit_dtset%qph2l)
 end if
 if(allocated(multibinit_dtset%q1shft))then
   ABI_DEALLOCATE(multibinit_dtset%q1shft)
 end if

 !if (allocated(multibinit_dtset%gilbert_damping))  then
 !  ABI_DEALLOCATE(multibinit_dtset%gilbert_damping)
 !end if

 !if (allocated(multibinit_dtset%gyro_ratio))  then
 !  ABI_DEALLOCATE(multibinit_dtset%gyro_ratio)
 !end if

 !if (allocated(multibinit_dtset%qph1l_spin))  then
 !  ABI_DEALLOCATE(multibinit_dtset%qph1l_spin)
 !end if
 !if (allocated(multibinit_dtset%qph2l_spin))  then
 !  ABI_DEALLOCATE(multibinit_dtset%qph2l_spin)
 !end if


end subroutine multibinit_dtset_free
!!***

!----------------------------------------------------------------------

!!****f* m_multibinit_dataset/invars10
!!
!! NAME
!! invars10
!!
!! FUNCTION
!! Open input file for the multibinit code, then reads or echoes the input information.
!!
!! INPUTS
!! lenstr=actual length of string
!! natom=number of atoms, needed for atifc
!! string*(*)=string of characters containing all input variables and data
!!
!! OUTPUT
!! multibinit_dtset <type(multibinit_dtset_type)> = datatype with all the input variables
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! PARENTS
!!      multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine invars10(multibinit_dtset,lenstr,natom,string)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: lenstr,natom
 character(len=*),intent(in) :: string
 type(multibinit_dtset_type),intent(inout) :: multibinit_dtset

!Local variables -------------------------
!Dummy arguments for subroutine 'intagm' to parse input file
!Set routine version number here:
!scalars
 integer :: iatifc,ii,iph1,iph2,jdtset,jj,marr,tread
 character(len=500) :: message
!arrays
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:),work(:)

!*********************************************************************
 marr=30
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

 jdtset=1

!copy natom to multibinit_dtset
 multibinit_dtset%natom=natom

!=====================================================================
!start reading in dimensions and non-dependent variables
!=====================================================================

!A
 multibinit_dtset%asr=2
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'asr',tread,'INT')
 if(tread==1) multibinit_dtset%asr=intarr(1)
 if(multibinit_dtset%asr<-2.or.multibinit_dtset%asr>5)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'asr is',multibinit_dtset%asr,', but the only allowed values',ch10,&
&   'are 0, 1, 2, 3, 4, 5, -1 or -2 .',ch10,&
&   'Action: correct asr in your input file.'
!  Note : negative values are allowed when the acoustic sum rule
!  is to be applied after the analysis of IFCs
!  3,4 are for rotational invariance (under development)
!  5 is for hermitian imposition of the ASR
   MSG_ERROR(message)
 end if

 multibinit_dtset%analyze_anh_pot=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'analyze_anh_pot',tread,'INT')
 if(tread==1) multibinit_dtset%analyze_anh_pot=intarr(1)
 if(multibinit_dtset%analyze_anh_pot < 0 .or. multibinit_dtset%analyze_anh_pot > 1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'analyze_anh_pot is',multibinit_dtset%analyze_anh_pot,', but the only allowed values',ch10,&
&   'are 0 and 1 .',ch10,&
&   'Action: correct analyze_anh_pot in your input file.'
   MSG_ERROR(message)
 end if

!B
 multibinit_dtset%brav=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'brav',tread,'INT')
 if(tread==1) multibinit_dtset%brav=intarr(1)
 if(multibinit_dtset%brav/=1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'brav is',multibinit_dtset%brav,', but the only allowed values',ch10,&
&   'are 1 for multibinit (not implemented) .',ch10,&
&   'Action: correct brav in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%bmass=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bmass',tread,'DPR')
 if(tread==1) multibinit_dtset%bmass=dprarr(1)
 if(multibinit_dtset%bmass<0)then
   write(message, '(a,f10.2,a,a,a,a,a)' )&
&   'bmass is',multibinit_dtset%bmass,', but the only allowed values',ch10,&
&   'is superior to 0.',ch10,&
&   'Action: correct bmass in your input file.'
   MSG_ERROR(message)
 end if


!C
 multibinit_dtset%chneut=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'chneut',tread,'INT')
 if(tread==1) multibinit_dtset%chneut=intarr(1)
 if(multibinit_dtset%chneut<0.or.multibinit_dtset%chneut>2)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'chneut is',multibinit_dtset%chneut,', but the only allowed values',ch10,&
&   'are 0, 1 or 2 .',ch10,&
&   'Action: correct chneut in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%confinement=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'confinement',tread,'INT')
 if(tread==1) multibinit_dtset%confinement=intarr(1)
 if(multibinit_dtset%confinement<0.or.multibinit_dtset%confinement>2)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'confinement is',multibinit_dtset%confinement,', but the only allowed values',ch10,&
&   'are 0, 1 or 2 .',ch10,&
&   'Action: correct confinement in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%conf_power_disp=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'conf_power_disp',tread,'INT')
 if(tread==1) multibinit_dtset%conf_power_disp=intarr(1)
 if(multibinit_dtset%conf_power_disp<0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'conf_power_disp is',multibinit_dtset%conf_power_disp,', but the only allowed values',ch10,&
&   'positive .',ch10,&
&   'Action: correct conf_power_disp in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%conf_power_strain=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'conf_power_strain',tread,'INT')
 if(tread==1) multibinit_dtset%conf_power_strain=intarr(1)
 if(multibinit_dtset%conf_power_strain<0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'conf_power_strain is',multibinit_dtset%conf_power_strain,', but the only allowed values',ch10,&
&   'are positive .',ch10,&
&   'Action: correct conf_power_strain in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%conf_power_fact_disp=100
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'conf_power_fact_disp',tread,'DPR')
 if(tread==1) multibinit_dtset%conf_power_fact_disp=dprarr(1)

 multibinit_dtset%conf_power_fact_strain=100
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'conf_power_fact_strain',tread,'DPR')
 if(tread==1) multibinit_dtset%conf_power_fact_strain=dprarr(1)

!D
 multibinit_dtset%dipdip=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dipdip',tread,'INT')
 if(tread==1) multibinit_dtset%dipdip=intarr(1)
 if(multibinit_dtset%dipdip>1.or.multibinit_dtset%dipdip<0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'dipdip is',multibinit_dtset%dipdip,', but the only allowed values',ch10,&
&   'is 1.',ch10,&
&   'Action: correct dipdip in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%dipdip_prt=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dipdip_prt',tread,'INT')
 if(tread==1) multibinit_dtset%dipdip_prt=intarr(1)
 if(multibinit_dtset%dipdip_prt<0.or.multibinit_dtset%dipdip_prt>1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'dipdip_prt is',multibinit_dtset%prtsrlr,', but the only allowed values',ch10,&
    'are 0 or 1.',ch10,&
&   'Action: correct dipdip_prt in your input file.'
   MSG_ERROR(message)
 end if


 multibinit_dtset%dtion=100
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dtion',tread,'INT')
 if(tread==1) multibinit_dtset%dtion=intarr(1)
 if(multibinit_dtset%dtion<1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'dtion is',multibinit_dtset%dtion,', but the only allowed values',ch10,&
&   'is superior to 1.',ch10,&
&   'Action: correct dtion in your input file.'
   MSG_ERROR(message)
 end if


 multibinit_dtset%delta_df= 1d-02
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'delta_df',tread,'DPR')
 if(tread==1) multibinit_dtset%delta_df=dprarr(1)
 if(multibinit_dtset%delta_df<0)then
   write(message, '(a,es10.2,a,a,a,a,a)' )&
&   'delta_df is',multibinit_dtset%delta_df,', but the only allowed values',ch10,&
&   'are superior to 0  .',ch10,&
&   'Action: correct delta_df in your input file.'
   MSG_ERROR(message)
 end if

!E
 multibinit_dtset%energy_reference= zero
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'energy_reference',tread,'DPR')
 if(tread==1) multibinit_dtset%energy_reference=dprarr(1)

 multibinit_dtset%enunit=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'enunit',tread,'INT')
 if(tread==1) multibinit_dtset%enunit=intarr(1)
 if(multibinit_dtset%enunit<0.or.multibinit_dtset%enunit>2)then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'enunit is',multibinit_dtset%enunit,', but the only allowed values',ch10,&
&   'are 0, 1 or 2.',ch10,&
&   'Action: correct enunit in your input file.'
   MSG_ERROR(message)
 end if

!F
 multibinit_dtset%fit_option=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_option',tread,'INT')
 if(tread==1) multibinit_dtset%fit_option=intarr(1)
 if(multibinit_dtset%fit_option<0.or.multibinit_dtset%fit_option>2)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'fit_option is',multibinit_dtset%fit_option,', but the only allowed values',ch10,&
&   'are 0, 1 or 2 for multibinit.',ch10,&
&   'Action: correct fit_option in your input file.'
   MSG_ERROR(message)
 end if


 multibinit_dtset%fit_ncoeff=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_ncoeff',tread,'INT')
 if(tread==1) multibinit_dtset%fit_ncoeff=intarr(1)
 if(multibinit_dtset%fit_ncoeff<0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'fit_ncoeff is',multibinit_dtset%fit_ncoeff,', but the only allowed values',ch10,&
&   'are positives for multibinit.',ch10,&
&   'Action: correct fit_ncoeff in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%fit_nbancoeff=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_nbancoeff',tread,'INT')
 if(tread==1) multibinit_dtset%fit_nbancoeff=intarr(1)
 if(multibinit_dtset%fit_nbancoeff<0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'fit_nbancoeff is',multibinit_dtset%fit_nbancoeff,', but the only allowed values',ch10,&
&   'are 0 or positive values for multibinit.',ch10,&
&   'Action: correct fit_nbancoeff in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%fit_nfixcoeff=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_nfixcoeff',tread,'INT')
 if(tread==1) multibinit_dtset%fit_nfixcoeff=intarr(1)
 if(multibinit_dtset%fit_nfixcoeff<-2)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'fit_nfixcoeff is',multibinit_dtset%fit_nfixcoeff,', but the only allowed values',ch10,&
&   'are -1 or positives for multibinit.',ch10,&
&   'Action: correct fit_nfixcoeff in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%ts_option=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ts_option',tread,'INT')
 if(tread==1) multibinit_dtset%ts_option=intarr(1)
 if(multibinit_dtset%ts_option<0.or.multibinit_dtset%ts_option>1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'ts_option is',multibinit_dtset%ts_option,', but the only allowed values',ch10,&
&   'are positives for multibinit.',ch10,&
&   'Action: correct ts_option in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%ifcana=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ifcana',tread,'INT')
 if(tread==1) multibinit_dtset%ifcana=intarr(1)
 if(multibinit_dtset%ifcana<0.or.multibinit_dtset%ifcana>1)then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'ifcana is',multibinit_dtset%ifcana,', but the only allowed values',ch10,&
&   'are 0 or 1.',ch10,&
&   'Action: correct ifcana in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%ifcflag=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ifcflag',tread,'INT')
 if(tread==1) multibinit_dtset%ifcflag=intarr(1)
 if(multibinit_dtset%ifcflag<0.or.multibinit_dtset%ifcflag>1)then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'ifcflag is',multibinit_dtset%ifcflag,', but the only allowed values',ch10,&
&   'are 0 or 1.',ch10,&
&   'Action: correct ifcflag in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%prtsrlr=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtsrlr',tread,'INT')
 if(tread==1) multibinit_dtset%prtsrlr=intarr(1)
 if(multibinit_dtset%prtsrlr<0.or.multibinit_dtset%prtsrlr>1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'prtsrlr is',multibinit_dtset%prtsrlr,', but the only allowed values',ch10,&
&   'are 0 or 1.',ch10,&
&   'Action: correct prtsrlr in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%ifcout=2000000 ! or -1 -> max number of ifc
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ifcout',tread,'INT')
 if(tread==1) multibinit_dtset%ifcout=intarr(1)
 if(multibinit_dtset%ifcout<-1)then
   write(message, '(a,i0,a,a,a)' )&
&   'ifcout is',multibinit_dtset%ifcout,', which is lower than -1 (default = all ifc) .',ch10,&
&   'Action: correct ifcout in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%nctime=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nctime',tread,'INT')
 if(tread==1) multibinit_dtset%nctime=intarr(1)
 if(multibinit_dtset%nctime<=0)then
   write(message, '(a,i0,a,a,a)' )&
&   'nctime is',multibinit_dtset%ntime,', which is not positive .',ch10,&
&   'Action: correct nctime in your input file.'
   MSG_ERROR(message)
 end if


 multibinit_dtset%ntime=200
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ntime',tread,'INT')
 if(tread==1) multibinit_dtset%ntime=intarr(1)
 if(multibinit_dtset%ntime<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'ntime is',multibinit_dtset%ntime,', which is lower than 0 .',ch10,&
&   'Action: correct ntime in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%dynamics=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dynamics',tread,'INT')
 if(tread==1) multibinit_dtset%dynamics=intarr(1)
 ! >100: the builtin multibinit lattice movers
 if(multibinit_dtset%dynamics/=0.and.multibinit_dtset%dynamics/=6.and.&
      &   multibinit_dtset%dynamics/=12.and.multibinit_dtset%dynamics/=13.and.&
      &   multibinit_dtset%dynamics/=27.and.&
      &   multibinit_dtset%dynamics/=24.and.multibinit_dtset%dynamics/=25 .and. &
      &   multibinit_dtset%dynamics/=101.and.multibinit_dtset%dynamics/=102 .and. & 
      &   multibinit_dtset%dynamics/=103.and.multibinit_dtset%dynamics/=120    &
    ) then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'dynamics is',multibinit_dtset%dynamics,', but the only allowed values',ch10,&
&   'are 6,12,24,25 or  13, 101, 102, 103 120 (see ionmov in abinit documentation).',ch10,&
&   'Action: correct dynamics in your input file.'
   MSG_ERROR(message)
 end if

 if(multibinit_dtset%dynamics==120) then
    write(message, '(a,i8,a)' )&
         &   'dynamics is',multibinit_dtset%dynamics,'The atoms will not move. For test only!'
    MSG_WARNING(message)
 end if


!L
 multibinit_dtset%latt_compressibility=0.0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'latt_compressibility',tread,'DPR')
 if(tread==1) multibinit_dtset%latt_compressibility=dprarr(1)

 multibinit_dtset%latt_friction=1e-4
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'latt_friction',tread,'DPR')
 if(tread==1) multibinit_dtset%latt_friction=dprarr(1)

 multibinit_dtset%latt_taut=1000
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'latt_taut',tread,'DPR')
 if(tread==1) multibinit_dtset%latt_taut=dprarr(1)

 multibinit_dtset%latt_taup=1000
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'latt_taup',tread,'DPR')
 if(tread==1) multibinit_dtset%latt_taup=dprarr(1)


!N
 multibinit_dtset%natifc=natom
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natifc',tread,'INT')
 if(tread==1) multibinit_dtset%natifc=intarr(1)
 if(multibinit_dtset%natifc<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'natifc is',multibinit_dtset%natifc,', which is lower than 0 .',ch10,&
&   'Action: correct natifc in your input file.'
   MSG_ERROR(message)
 end if

 if(multibinit_dtset%natifc>natom)then
   write(message, '(a,i0,a,a,a,i0,a,a,a)' )&
&   'The number of atom ifc in the input files',multibinit_dtset%natifc,',',ch10,&
&   'is larger than the number of atoms',natom,'.',ch10,&
&   'Action: change natifc in the input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%ncoeff=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ncoeff',tread,'INT')
 if(tread==1) multibinit_dtset%ncoeff=intarr(1)
 if(multibinit_dtset%ncoeff<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'ncoeff is',multibinit_dtset%ncoeff,', which is lower than 0 .',ch10,&
&   'Action: correct natifc in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%ng2qpt(:)=0
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'ng2qpt',tread,'INT')
 if(tread==1) multibinit_dtset%ng2qpt(:)=intarr(1:3)
 do ii=1,3
   if(multibinit_dtset%ng2qpt(ii)<0)then
     write(message, '(a,i0,a,i0,a,a,a,i0,a)' )&
&     'ng2qpt(',ii,') is',multibinit_dtset%ng2qpt(ii),', which is lower than 0 .',ch10,&
&     'Action: correct ng2qpt(',ii,') in your input file.'
     MSG_ERROR(message)
   end if
 end do

 multibinit_dtset%ncell(:)= 1
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'ncell',tread,'INT')
 if(tread==1) multibinit_dtset%ncell(1:3)=intarr(1:3)
 do ii=1,3
   if(multibinit_dtset%ncell(ii)<0.or.multibinit_dtset%ncell(ii)>100)then
     write(message, '(a,i0,a,i0,3a,i0,a)' )&
&     'ncell(',ii,') is ',multibinit_dtset%ncell(ii),', which is lower than 0 of superior than 50.',&
&     ch10,'Action: correct ncell(',ii,') in your input file.'
     MSG_ERROR(message)
   end if
 end do

 multibinit_dtset%ngqpt(:)= 1
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'ngqpt',tread,'INT')
 if(tread==1) multibinit_dtset%ngqpt(1:3)=intarr(1:3)
 do ii=1,3
   if(multibinit_dtset%ngqpt(ii)<0)then
     write(message, '(a,i0,a,i0,a,a,a,i0,a)' )&
&     'ngqpt(',ii,') is',multibinit_dtset%ngqpt(ii),', which is lower than 0 .',ch10,&
&     'Action: correct ngqpt(',ii,') in your input file.'
     MSG_ERROR(message)
   end if
 end do

 multibinit_dtset%nph1l=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nph1l',tread,'INT')
 if(tread==1) multibinit_dtset%nph1l=intarr(1)
 if(multibinit_dtset%nph1l<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'nph1l is',multibinit_dtset%nph1l,', which is lower than 0 .',ch10,&
&   'Action: correct nph1l in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%nph2l=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nph2l',tread,'INT')
 if(tread==1) multibinit_dtset%nph2l=intarr(1)
 if(multibinit_dtset%nph2l<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'nph2l is',multibinit_dtset%nph2l,', which is lower than 0 .',ch10,&
&   'Action: correct nph2l in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%nqshft=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nqshft',tread,'INT')
 if(tread==1) multibinit_dtset%nqshft=intarr(1)
 if(multibinit_dtset%nqshft<0 .or. multibinit_dtset%nqshft==3 .or.&
& multibinit_dtset%nqshft>=5 )then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'nqshft is',multibinit_dtset%nqshft,', but the only allowed values',ch10,&
&   'are 1, 2 or 4 .',ch10,&
&   'Action: correct nqshft in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%nnos=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nnos',tread,'INT')
 if(tread==1) multibinit_dtset%nnos=intarr(1)
 if(multibinit_dtset%nnos<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'nnos is',multibinit_dtset%nnos,', which is lower than 0',ch10,&
&   'Action: correct nnos in your input file.'
   MSG_ERROR(message)
 end if


 multibinit_dtset%nsphere=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nsphere',tread,'INT')
 if(tread==1) multibinit_dtset%nsphere=intarr(1)
 if(multibinit_dtset%nsphere<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'nsphere is',multibinit_dtset%nsphere,', which is lower than 0',ch10,&
&   'Action: correct nsphere in your input file.'
   MSG_ERROR(message)
 end if

!O
 multibinit_dtset%optcell=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'optcell',tread,'INT')
 if(tread==1) multibinit_dtset%optcell=intarr(1)
 if(multibinit_dtset%optcell<0.or.multibinit_dtset%optcell>2)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'optcell is',multibinit_dtset%prtsrlr,', but the only allowed values',ch10,&
&   'are 0, 1 or 2.',ch10,&
&   'Action: correct optcell in your input file.'
   MSG_ERROR(message)
 end if


!P
 multibinit_dtset%prt_model=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prt_model',tread,'INT')
 if(tread==1) multibinit_dtset%prt_model=intarr(1)
 if(multibinit_dtset%prt_model<0.or.multibinit_dtset%prt_model>4)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'prt_model is',multibinit_dtset%prtsrlr,', but the only allowed values',ch10,&
&   'are 0, 1 or 2.',ch10,&
&   'Action: correct prt_model in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%prt_names=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prt_names',tread,'INT')
 if(tread==1) multibinit_dtset%prt_names=intarr(1)
 if(multibinit_dtset%prt_names<0.or.multibinit_dtset%prt_names>1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'prt_names is',multibinit_dtset%prt_names,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct prt_names in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%prt_phfrq=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prt_phfrq',tread,'INT')
 if(tread==1) multibinit_dtset%prt_phfrq=intarr(1)
 if(multibinit_dtset%prt_phfrq<0.or.multibinit_dtset%prt_phfrq>2)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'prt_phfrq is',multibinit_dtset%prtsrlr,', but the only allowed values',ch10,&
&   'are 0, 1 or 2.',ch10,&
&   'Action: correct prt_phfrq in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%fit_initializeData=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_initializeData',tread,'INT')
 if(tread==1) multibinit_dtset%fit_initializeData=intarr(1)
 if(multibinit_dtset%fit_initializeData<0.or.multibinit_dtset%fit_initializeData>1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'fit_initializeData is',multibinit_dtset%prtsrlr,', but the only allowed values',ch10,&
&   'are 0, 1 or 2.',ch10,&
&   'Action: correct fit_initializeData in your input file.'
   MSG_ERROR(message)
 end if

 
 multibinit_dtset%fit_generateCoeff=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_generateCoeff',tread,'INT')
 if(tread==1) multibinit_dtset%fit_generateCoeff=intarr(1)
 if(multibinit_dtset%fit_generateCoeff<0.or.multibinit_dtset%fit_generateCoeff>1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'fit_generateCoeff is',multibinit_dtset%prtsrlr,', but the only allowed values',ch10,&
&   'are 0, 1 or 2.',ch10,&
&   'Action: correct fit_generateCoeff in your input file.'
   MSG_ERROR(message)
 end if

!Default is no output of the real space IFC to file
 multibinit_dtset%prt_ifc = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prt_ifc',tread,'INT')
 if(tread==1) multibinit_dtset%prt_ifc = intarr(1)
 if(multibinit_dtset%prt_ifc < 0 .or. multibinit_dtset%prt_ifc > 1) then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'prtf_ifc is',multibinit_dtset%prt_ifc,'. The only allowed values',ch10,&
&   'are 0 (no output) or 1 (AI2PS format)',ch10,  &
&   'Action: correct prt_ifc in your input file.'
   MSG_ERROR(message)
 end if

!Default is no output of the 3rd derivative
 multibinit_dtset%strcpling = -1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'strcpling',tread,'INT')
 if(tread==1) multibinit_dtset%strcpling = intarr(1)
 if(multibinit_dtset%strcpling < -1 .or. multibinit_dtset%strcpling > 2) then
   write(message, '(a,i0,a,a,a,a,a,a,a)' )&
&   'prtf_3rd is ',multibinit_dtset%strcpling,'. The only allowed values',ch10,&
&   'are 0 (no computation), 1 (only computation)',ch10,&
&   'or 2 (computation and print in xml file)',ch10,  &
&   'Action: correct strcpling in your input file.'
   MSG_ERROR(message)
 end if

!Q
 multibinit_dtset%qrefine=1 ! default is no refinement
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'qrefine',tread,'INT')
 if(tread==1) multibinit_dtset%qrefine = intarr(1:3)
 do ii=1,3
   if(multibinit_dtset%qrefine(ii) < 1) then
     write(message, '(a,3i0,a,a,a,a,a)' )&
&     'qrefine is',multibinit_dtset%qrefine,' The only allowed values',ch10,&
&     'are integers >= 1 giving the refinement of the ngqpt grid',ch10,&
&     'Action: correct qrefine in your input file.'
     MSG_ERROR(message)
   end if
 end do

!R
 multibinit_dtset%restartxf=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'restartxf',tread,'INT')
 if(tread==1) multibinit_dtset%restartxf=intarr(1)
 if(multibinit_dtset%restartxf < -3 .or. multibinit_dtset%restartxf > 0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'restartxf is',multibinit_dtset%restartxf,', but the only allowed values',ch10,&
&   'is -2 or 0.',ch10,&
&   'Action: correct restartxf in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%rfmeth=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfmeth',tread,'INT')
 if(tread==1) multibinit_dtset%rfmeth=intarr(1)
 if(multibinit_dtset%rfmeth<1.or.multibinit_dtset%rfmeth>2)then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'rfmeth is',multibinit_dtset%rfmeth,', but the only allowed values',ch10,&
&   'are 1 or 2 . ',ch10,&
&   'Action: correct rfmeth in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%rifcsph=zero
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rifcsph',tread,'DPR')
 if(tread==1) multibinit_dtset%rifcsph=dprarr(1)
 if(multibinit_dtset%rifcsph<-tol12)then
   write(message, '(a,f10.3,a,a,a)' )&
&   'rifcsph is',multibinit_dtset%rifcsph,', which is lower than zero.',ch10,&
&   'Action: correct rifcsph in your input file.'
   MSG_ERROR(message)
 end if

!S
 multibinit_dtset%spin_damping=-1.0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_damping',tread,'DPR')
 if(tread==1) multibinit_dtset%spin_damping=dprarr(1)

 multibinit_dtset%spin_calc_correlation_obs=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_calc_correlation_obs',tread,'INT')
 if(tread==1) multibinit_dtset%spin_calc_correlation_obs=intarr(1)
 if(multibinit_dtset%spin_calc_correlation_obs>1.or.multibinit_dtset%spin_calc_correlation_obs<0)then
    write(message, '(a,i8,a,a,a,a,a)' )&
         &   'spin_calc_correlation_obs is',multibinit_dtset%spin_calc_correlation_obs,', but the only allowed values',ch10,&
         &   'is 0 or 1.',ch10,&
         &   'Action: correct spin_calc_correlation_obs in your input file.'
    MSG_ERROR(message)
 end if

 multibinit_dtset%spin_calc_thermo_obs=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_calc_thermo_obs',tread,'INT')
 if(tread==1) multibinit_dtset%spin_calc_thermo_obs=intarr(1)
 if(multibinit_dtset%spin_calc_thermo_obs>1.or.multibinit_dtset%spin_calc_thermo_obs<0)then
    write(message, '(a,i8,a,a,a,a,a)' )&
         &   'spin_calc_thermo_obs is',multibinit_dtset%spin_calc_thermo_obs,', but the only allowed values',ch10,&
         &   'is 0 or 1.',ch10,&
         &   'Action: correct spin_calc_thermo_obs in your input file.'
    MSG_ERROR(message)
 end if


 multibinit_dtset%spin_calc_traj_obs=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_calc_traj_obs',tread,'INT')
 if(tread==1) multibinit_dtset%spin_calc_traj_obs=intarr(1)
 if(multibinit_dtset%spin_calc_traj_obs>1.or.multibinit_dtset%spin_calc_traj_obs<0)then
    write(message, '(a,i8,a,a,a,a,a)' )&
         &   'spin_calc_traj_obs is',multibinit_dtset%spin_calc_traj_obs,', but the only allowed values',ch10,&
         &   'is 0 or 1.',ch10,&
         &   'Action: correct spin_calc_traj_obs in your input file.'
    MSG_ERROR(message)
 end if


 multibinit_dtset%spin_dipdip=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_dipdip',tread,'INT')
 if(tread==1) multibinit_dtset%spin_dipdip=intarr(1)
 if(multibinit_dtset%spin_dipdip>1.or.multibinit_dtset%spin_dipdip<0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'spin_dipdip is',multibinit_dtset%spin_dipdip,', but the only allowed values',ch10,&
&   'is 0 or 1.',ch10,&
&   'Action: correct spin_dipdip in your input file.'
   MSG_ERROR(message)
 end if


 multibinit_dtset%spin_dt= 1d-16
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_dt',tread,'TIM')
 if(tread==1) multibinit_dtset%spin_dt=dprarr(1)
 if(multibinit_dtset%spin_dt<0)then
    write(message, '(a,es10.2,a,a,a,a,a)' )&
         &   'spin_dt is',multibinit_dtset%spin_dt,', but the only allowed values',ch10,&
         &   'are superior to 0  .',ch10,&
         &   'Action: correct spin_dt in your input file.'
    MSG_ERROR(message)
 end if

 
 multibinit_dtset%spin_dynamics=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_dynamics',tread,'INT')
 if(tread==1) multibinit_dtset%spin_dynamics=intarr(1)
 if( .not. (multibinit_dtset%spin_dynamics <= 3 .or. multibinit_dtset%spin_dynamics==20) ) then
    write(message, '(a,i8,a,a,a,a,a)' )&
         &   'spin_dynamics is',multibinit_dtset%spin_dynamics,', but the only allowed values',ch10,&
         &   'are 0, 1, 2, 3 and 20 and negative values.',ch10,&
         &   'Action: correct spin_dynamics in your input file.'
    MSG_ERROR(message)
 end if

 if(multibinit_dtset%spin_dynamics == 20) then
    write(message, '(a,i8,a)' )&
         &   'spin_dynamics is',multibinit_dtset%spin_dynamics,', spins will not move. For test only!!'
    MSG_WARNING(message)
 end if


 
 multibinit_dtset%spin_init_state=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_init_state',tread,'INT')
 if(tread==1) multibinit_dtset%spin_init_state=intarr(1)
 if(multibinit_dtset%spin_init_state<1 .or. &
      &   multibinit_dtset%spin_init_state>5) then
    write(message, '(a,i8,a,a,a,a,a)' )&
         &   'spin_init_state is',multibinit_dtset%spin_init_state,', but the only allowed values',ch10,&
         &   'are 1, 2, 3, 4, and 5.',ch10,&
         &   'Action: correct spin_init_state in your input file.'
    MSG_ERROR(message)
 end if
 


 multibinit_dtset%spin_mag_field= zero
 if(3>marr)then
    marr=3
    ABI_DEALLOCATE(intarr)
    ABI_DEALLOCATE(dprarr)
    ABI_ALLOCATE(intarr,(marr))
    ABI_ALLOCATE(dprarr,(marr))
 end if
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'spin_mag_field',tread,'BFI')
 if(tread==1) multibinit_dtset%spin_mag_field(1:3)= dprarr(1:3)

 multibinit_dtset%spin_nctime=100
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_nctime',tread,'INT')
 if(tread==1) multibinit_dtset%spin_nctime=intarr(1)
 if(multibinit_dtset%spin_nctime<0)then
    write(message, '(a,i0,a,a,a)' )&
         &   'spin_nctime is',multibinit_dtset%spin_nctime,', which is lower than 0 .',ch10,&
         &   'Action: correct spin_nctime in your input file.'
    MSG_ERROR(message)
 end if
 
 multibinit_dtset%spin_ntime_pre=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_ntime_pre',tread,'INT')
 if(tread==1) multibinit_dtset%spin_ntime_pre=intarr(1)
 if(multibinit_dtset%spin_ntime_pre<0)then
    write(message, '(a,i0,a,a,a)' )&
         &   'spin_ntime_pre is',multibinit_dtset%spin_ntime_pre,', which is lower than 0 .',ch10,&
         &   'Action: correct spin_ntime_pre in your input file.'
    MSG_ERROR(message)
 end if



 multibinit_dtset%spin_ntime=10000
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_ntime',tread,'INT')
 if(tread==1) multibinit_dtset%spin_ntime=intarr(1)
 if(multibinit_dtset%spin_ntime<0)then
    write(message, '(a,i0,a,a,a)' )&
         &   'spin_ntime is',multibinit_dtset%spin_ntime,', which is lower than 0 .',ch10,&
         &   'Action: correct spin_ntime in your input file.'
    MSG_ERROR(message)
 end if


 multibinit_dtset%spin_n1l=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_n1l',tread,'INT')
 if(tread==1) multibinit_dtset%spin_n1l=intarr(1)
 if(multibinit_dtset%spin_n1l<0)then
    write(message, '(a,i0,a,a,a)' )&
         &   'spin_n1l is',multibinit_dtset%spin_n1l,', which is lower than 0 .',ch10,&
         &   'Action: correct spin_n1l in your input file.'
    MSG_ERROR(message)
 end if

 multibinit_dtset%spin_n2l=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_n2l',tread,'INT')
 if(tread==1) multibinit_dtset%spin_n2l=intarr(1)
 if(multibinit_dtset%spin_n2l<0)then
    write(message, '(a,i0,a,a,a)' )&
         &   'spin_n2l is',multibinit_dtset%spin_n2l,', which is lower than 0 .',ch10,&
         &   'Action: correct spin_n2l in your input file.'
    MSG_ERROR(message)
 end if
 
 multibinit_dtset%spin_init_orientation= [0.0, 0.0, 1.0]
 if(3>marr)then
    marr=3
    ABI_DEALLOCATE(intarr)
    ABI_DEALLOCATE(dprarr)
    ABI_ALLOCATE(intarr,(marr))
    ABI_ALLOCATE(dprarr,(marr))
 end if
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'spin_init_orientation',tread,'DPR')
 if(tread==1) multibinit_dtset%spin_init_orientation(1:3)= dprarr(1:3)

 multibinit_dtset%spin_qpoint= zero
 if(3>marr)then
    marr=3
    ABI_DEALLOCATE(intarr)
    ABI_DEALLOCATE(dprarr)
    ABI_ALLOCATE(intarr,(marr))
    ABI_ALLOCATE(dprarr,(marr))
 end if
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'spin_qpoint',tread,'DPR')
 if(tread==1) multibinit_dtset%spin_qpoint(1:3)= dprarr(1:3)

 multibinit_dtset%spin_init_qpoint= zero
 if(3>marr)then
    marr=3
    ABI_DEALLOCATE(intarr)
    ABI_DEALLOCATE(dprarr)
    ABI_ALLOCATE(intarr,(marr))
    ABI_ALLOCATE(dprarr,(marr))
 end if
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'spin_init_qpoint',tread,'DPR')
 if(tread==1) multibinit_dtset%spin_init_qpoint(1:3)= dprarr(1:3)

 multibinit_dtset%spin_init_rotate_axis= [1.0, 0.0, 0.0]
 if(3>marr)then
    marr=3
    ABI_DEALLOCATE(intarr)
    ABI_DEALLOCATE(dprarr)
    ABI_ALLOCATE(intarr,(marr))
    ABI_ALLOCATE(dprarr,(marr))
 end if
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'spin_init_rotate_axis',tread,'DPR')
 if(tread==1) multibinit_dtset%spin_init_rotate_axis(1:3)= dprarr(1:3)

 multibinit_dtset%spin_sia_add=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_sia_add',tread,'INT')
 if(tread==1) multibinit_dtset%spin_sia_add=intarr(1)
 if(multibinit_dtset%spin_sia_add <0 .or. multibinit_dtset%spin_sia_add>2 )then
    write(message, '(a,i0,a,a,a)' )&
         &   'spin_sia_add is',multibinit_dtset%spin_sia_add,', which is not 0, 1, or 2.',ch10,&
         &   'Action: correct spin_sia_add in your input file.'
    MSG_ERROR(message)
 end if

 multibinit_dtset%spin_sia_k1amp=0.0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_sia_k1amp',tread,'ENE')
 if(tread==1) multibinit_dtset%spin_sia_k1amp=dprarr(1)

 multibinit_dtset%spin_sia_k1dir(:)= [0.0,0.0,1.0]
 if(3>marr)then
    marr=3
    ABI_DEALLOCATE(intarr)
    ABI_DEALLOCATE(dprarr)
    ABI_ALLOCATE(intarr,(marr))
    ABI_ALLOCATE(dprarr,(marr))
 end if
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'spin_sia_k1dir',tread,'DPR')
 if(tread==1) then
    dprarr(1:3)=dprarr(1:3)/sqrt(sum(dprarr(1:3)**2))
    multibinit_dtset%spin_sia_k1dir(1:3)= dprarr(1:3)
 endif

 multibinit_dtset%spin_temperature=325
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_temperature',tread,'DPR')
 if(tread==1) multibinit_dtset%spin_temperature=dprarr(1)
 if(multibinit_dtset%spin_temperature<0)then
   write(message, '(a,f10.1,a,a,a,a,a)' )&
&   'spin_temperature is ',multibinit_dtset%spin_temperature,'. The only allowed values',ch10,&
&   'are non-negative values.',ch10,&
&   'Action: correct spin_temperature in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%spin_temperature_start=0.0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_temperature_start',tread,'DPR')
 if(tread==1) multibinit_dtset%spin_temperature_start=dprarr(1)
 if(multibinit_dtset%spin_temperature_start<0.0)then
    write(message, '(a,f10.1,a,a,a,a,a)' )&
         &   'spin_temperature_start is ',multibinit_dtset%spin_temperature_start,'. The only allowed values',ch10,&
         &   'are positives values.',ch10,&
         &   'Action: correct spin_semperature_start in your input file.'
    MSG_ERROR(message)
 end if

 multibinit_dtset%spin_temperature_end=0.0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_temperature_end',tread,'DPR')
 if(tread==1) multibinit_dtset%spin_temperature_end=dprarr(1)
 if(multibinit_dtset%spin_temperature_end<0)then
    write(message, '(a,f10.1,a,a,a,a,a)' )&
         &   'spin_temperature_end is ',multibinit_dtset%spin_temperature_end,'. The only allowed values',ch10,&
         &   'are positives values.',ch10,&
         &   'Action: correct spin_semperature_end in your input file.'
    MSG_ERROR(message)
 end if

 multibinit_dtset%spin_temperature_nstep=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_temperature_nstep',tread,'INT')
 if(tread==1) multibinit_dtset%spin_temperature_nstep=intarr(1)
 if(multibinit_dtset%spin_temperature_nstep<=0)then
    write(message, '(a,i0,a,a,a,a)' )&
         &   'spin_temperature_nstep is',multibinit_dtset%spin_temperature_nstep,', while it should be larger than 0',ch10,&
         &   'Action: correct spin_temperature_nstep in your input file.'
    MSG_ERROR(message)
 end if



 multibinit_dtset%spin_tolavg=1d-02
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_tolavg',tread,'DPR')
 if(tread==1) multibinit_dtset%spin_tolavg=dprarr(1)
 if(multibinit_dtset%spin_tolavg<=0)then
    write(message, '(a,f10.1,a,a,a,a,a)' )&
         &   'spin_tolavg is ',multibinit_dtset%spin_tolavg,'. The only allowed values',ch10,&
         &   'are positives values.',ch10,&
         &   'Action: correct spin_tolavg in your input file.'
    MSG_ERROR(message)
 end if

 multibinit_dtset%spin_tolvar=1d-02
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_tolvar',tread,'DPR')
 if(tread==1) multibinit_dtset%spin_tolvar=dprarr(1)
 if(multibinit_dtset%spin_tolvar<=0)then
    write(message, '(a,f10.1,a,a,a,a,a)' )&
         &   'spin_tolvar is ',multibinit_dtset%spin_tolvar,'. The only allowed values',ch10,&
         &   'are positives values.',ch10,&
         &   'Action: correct spin_tolvar in your input file.'
    MSG_ERROR(message)
 end if
 
 multibinit_dtset%spin_var_temperature=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_var_temperature',tread,'INT')
 if(tread==1) multibinit_dtset%spin_var_temperature=intarr(1)
 if(multibinit_dtset%spin_var_temperature/=0.and.multibinit_dtset%spin_var_temperature/=1)then
    write(message, '(a,i0,a,a,a,a,a)' )&
         &   'spin_var_temperature is',multibinit_dtset%spin_var_temperature,'. The only allowed values',ch10,&
         &   'are 0, or 1.',ch10,&
         &   'Action: correct spin_var_temperature in your input file.'
    MSG_ERROR(message)
 end if

 multibinit_dtset%spin_write_traj=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spin_write_traj',tread,'INT')
 if(tread==1) multibinit_dtset%spin_write_traj=intarr(1)
 if(multibinit_dtset%spin_write_traj/=0.and.multibinit_dtset%spin_write_traj/=1)then
    write(message, '(a,i0,a,a,a,a,a)' )&
         &   'spin_write_traj is',multibinit_dtset%spin_write_traj,'. The only allowed values',ch10,&
         &   'are 0, or 1.',ch10,&
         &   'Action: correct spin_write_traj in your input file.'
    MSG_ERROR(message)
 end if

 multibinit_dtset%slc_coupling=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'slc_coupling',tread,'INT')
 if(tread==1) multibinit_dtset%slc_coupling=intarr(1)
 if(multibinit_dtset%slc_coupling.ne. 1111 .and. & 
      &   multibinit_dtset%slc_coupling.ne. 1110 .and. &
      &   multibinit_dtset%slc_coupling.ne. 1101 .and. &
      &   multibinit_dtset%slc_coupling.ne. 1011 .and. &
      &   multibinit_dtset%slc_coupling.ne.  111 .and. &
      &   multibinit_dtset%slc_coupling.ne. 1100 .and. &
      &   multibinit_dtset%slc_coupling.ne. 1010 .and. &
      &   multibinit_dtset%slc_coupling.ne. 1001 .and. &
      &   multibinit_dtset%slc_coupling.ne.  110 .and. &
      &   multibinit_dtset%slc_coupling.ne.  101 .and. &
      &   multibinit_dtset%slc_coupling.ne.   11 .and. &
      &   multibinit_dtset%slc_coupling.ne. 1000 .and. &
      &   multibinit_dtset%slc_coupling.ne.  100 .and. &
      &   multibinit_dtset%slc_coupling.ne.   10 .and. &
      &   multibinit_dtset%slc_coupling.ne.    1 .and. &
      &   multibinit_dtset%slc_coupling.ne.    0) then
    write(message, '(a,i8,a,a,a,a,a)' )&
         &   'slc_coupling is',multibinit_dtset%slc_coupling,', but the only allowed values',ch10,&
         &   'are 1111, 1110, 1101, 1011, 111, 1100, 1010, 1001, 110, 101, 11, 1000, 100, 10, 1, and 0.',ch10,&
         &   'Action: correct slc_coupling in your input file.'
    MSG_ERROR(message)
 end if


 multibinit_dtset%symdynmat=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'symdynmat',tread,'INT')
 if(tread==1) multibinit_dtset%symdynmat=intarr(1)
 if(multibinit_dtset%symdynmat/=0.and.multibinit_dtset%symdynmat/=1)then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'symdynmat is',multibinit_dtset%symdynmat,'. The only allowed values',ch10,&
&   'are 0, or 1.',ch10,&
&   'Action: correct symdynmat in your input file.'
   MSG_ERROR(message)
 end if



!T


 multibinit_dtset%temperature=325
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'temperature',tread,'DPR')
 if(tread==1) multibinit_dtset%temperature=dprarr(1)
 if(multibinit_dtset%temperature<=0)then
   write(message, '(a,f10.1,a,a,a,a,a)' )&
&   'Temperature is ',multibinit_dtset%temperature,'. The only allowed values',ch10,&
&   'are positives values.',ch10,&
&   'Action: correct Temperature in your input file.'
   MSG_ERROR(message)
 end if



!U

!V

!W

!X

!Y

!Z

!=====================================================================
!end non-dependent variables
!=====================================================================

!=======================================================================
!Read in dependent variables (dependent on dimensions above)
!=======================================================================

!A
 multibinit_dtset%acell= one
 if(3>marr)then
   marr=3
   ABI_DEALLOCATE(intarr)
   ABI_DEALLOCATE(dprarr)
   ABI_ALLOCATE(intarr,(marr))
   ABI_ALLOCATE(dprarr,(marr))
 end if
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'acell',tread,'DPR')
 if(tread==1) multibinit_dtset%acell(1:3)= dprarr(1:3)
 if(any(multibinit_dtset%acell<=tol10))then
    write(message, '(3a)' )&
&       'There is negative or zero value for cell ',ch10,&
&       'Action: change acell in your input file.'
      MSG_ERROR(message)
 end if

 if(6>marr)then
   marr=6
   ABI_DEALLOCATE(intarr)
   ABI_DEALLOCATE(dprarr)
   ABI_ALLOCATE(intarr,(marr))
   ABI_ALLOCATE(dprarr,(marr))
 end if
 multibinit_dtset%strtarget(1:6) = zero
 call intagm(dprarr,intarr,jdtset,marr,6,string(1:lenstr),'strtarget',tread,'DPR')
 if(tread==1) multibinit_dtset%strtarget(1:6)=dprarr(1:6)


 ABI_ALLOCATE(multibinit_dtset%atifc,(natom))
 multibinit_dtset%atifc(:)=0
 if(multibinit_dtset%natifc>=1)then
   if(multibinit_dtset%natifc>marr)then
     marr=multibinit_dtset%natifc
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   call intagm(dprarr,intarr,jdtset,marr,multibinit_dtset%natifc,string(1:lenstr),'atifc',tread,'INT')
   if(tread==1) then
     multibinit_dtset%atifc(1:multibinit_dtset%natifc)= intarr(1:multibinit_dtset%natifc)
   else ! set to the maximum
     do iatifc=1,multibinit_dtset%natifc
       multibinit_dtset%atifc(iatifc) =  iatifc
     end do
   end if
   ABI_MALLOC(work,(natom))
   work(:)=0

   do iatifc=1,multibinit_dtset%natifc
     if(multibinit_dtset%atifc(iatifc)<=0.or.multibinit_dtset%atifc(iatifc)>natom)then
       write(message, '(a,i0,a,a,a,a,a,i0,a,a,a)' )&
&       'For iatifc=',iatifc,', the number of the atom ifc to be ',ch10,&
&       'analysed is not valid : either negative, ',ch10,&
&       'zero, or larger than natom =',natom,'.',ch10,&
&       'Action: change atifc in your input file.'
       MSG_ERROR(message)
     end if
     work(multibinit_dtset%atifc(iatifc))=1
   end do
   multibinit_dtset%atifc(1:natom)=int(work(:))
   ABI_FREE(work)
 end if

!B

!C
 ABI_ALLOCATE(multibinit_dtset%coefficients,(multibinit_dtset%ncoeff))
 if (multibinit_dtset%ncoeff/=0)then
   if(multibinit_dtset%ncoeff>marr)then
     marr=multibinit_dtset%ncoeff
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   multibinit_dtset%coefficients(:)=zero
   call intagm(dprarr,intarr,jdtset,marr,multibinit_dtset%ncoeff,&
&              string(1:lenstr),'coefficients',tread,'DPR')
   if(tread==1)then
     do ii=1,multibinit_dtset%ncoeff
       multibinit_dtset%coefficients(ii)=dprarr(ii)
     end do
   end if
 end if

 ABI_ALLOCATE(multibinit_dtset%conf_cutoff_disp,(multibinit_dtset%natom))
 if (multibinit_dtset%natom/=0)then
   if(multibinit_dtset%natom>marr)then
     marr=multibinit_dtset%natom
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   multibinit_dtset%conf_cutoff_disp(:)=zero
   call intagm(dprarr,intarr,jdtset,marr,multibinit_dtset%natom,&
&              string(1:lenstr),'conf_cutoff_disp',tread,'DPR')
   if(tread==1)then
     do ii=1,multibinit_dtset%natom
       multibinit_dtset%conf_cutoff_disp(ii)=dprarr(ii)
     end do
   end if
   if(any(multibinit_dtset%conf_cutoff_disp<zero))then
     write(message, '(3a)' )&
&       'There is negative value for conf_cutoff_disp ',ch10,&
&       'Action: change acell in your input file.'
     MSG_ERROR(message)
   end if
 end if

 if(6>marr)then
   marr=6
   ABI_DEALLOCATE(intarr)
   ABI_DEALLOCATE(dprarr)
   ABI_ALLOCATE(intarr,(marr))
   ABI_ALLOCATE(dprarr,(marr))
 end if
 multibinit_dtset%conf_cutoff_strain(1:6) = zero
 call intagm(dprarr,intarr,jdtset,marr,6,string(1:lenstr),'conf_cutoff_strain',tread,'DPR')
 if(tread==1) multibinit_dtset%conf_cutoff_strain(1:6)=dprarr(1:6)
 if(any(multibinit_dtset%conf_cutoff_disp<zero))then
   write(message, '(3a)' )&
&     'There is negative value for conf_cutoff_strain ',ch10,&
&     'Action: change acell in your input file.'
   MSG_ERROR(message)
 end if

!D
 multibinit_dtset%dipdip_range(:)= (/0,0,0/)
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'dipdip_range',tread,'INT')
 if(tread==1) multibinit_dtset%dipdip_range(1:3)=intarr(1:3)
 do ii=1,3
   if(multibinit_dtset%dipdip_range(ii)<0.or.multibinit_dtset%dipdip_range(ii)>50)then
     write(message, '(a,i0,a,i0,4a,i0,a)' )&
&     'dipdip_range(',ii,') is ',multibinit_dtset%dipdip_range(ii),', which is lower',&
&     ' than 0 of superior than 50.',&
&     ch10,'Action: correct dipdip_range(',ii,') in your input file.'
     MSG_ERROR(message)
   end if
 end do
!E
 multibinit_dtset%eivec=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'eivec',tread,'INT')
 if(tread==1) multibinit_dtset%eivec=intarr(1)
 if(multibinit_dtset%eivec<0.or.multibinit_dtset%eivec>4)then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'eivec is',multibinit_dtset%eivec,', but the only allowed values',ch10,&
&   'are 0, 1, 2, 3 or 4.',ch10,&
&   'Action: correct eivec in your input file.'
   MSG_ERROR(message)
 end if

!F
  multibinit_dtset%fit_anhaStrain=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_anhaStrain',tread,'INT')
 if(tread==1) multibinit_dtset%fit_anhaStrain=intarr(1)
 if(multibinit_dtset%fit_anhaStrain<0.and.multibinit_dtset%fit_anhaStrain>1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'fit_anhaStrain is',multibinit_dtset%fit_anhaStrain,', but the only allowed values',ch10,&
&   'are 0 or 1 for multibinit.',ch10,&
&   'Action: correct fit_anhaStrain in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%bound_model=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bound_model',tread,'INT')
 if(tread==1) multibinit_dtset%bound_model=intarr(1)
 if(multibinit_dtset%bound_model<0.and.multibinit_dtset%bound_model>1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'bound_model is',multibinit_dtset%bound_model,', but the only allowed values',ch10,&
&   'are 0 or 1 for multibinit.',ch10,&
&   'Action: correct bound_model in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%bound_anhaStrain=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_anhaStrain',tread,'INT')
 if(tread==1) multibinit_dtset%bound_anhaStrain=intarr(1)
 if(multibinit_dtset%bound_anhaStrain<0.and.multibinit_dtset%bound_anhaStrain>1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'fit_anhaStrain is',multibinit_dtset%bound_anhaStrain,', but the only allowed values',ch10,&
&   'are 0 or 1 for multibinit.',ch10,&
&   'Action: correct fit_anhaStrain in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%bound_SPCoupling=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bound_SPCoupling',tread,'INT')
 if(tread==1) multibinit_dtset%bound_SPCoupling=intarr(1)
 if(multibinit_dtset%bound_SPCoupling<0.and.multibinit_dtset%bound_SPCoupling>1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
     &   'fit_SPCoupling is',multibinit_dtset%bound_SPCoupling,&
     &   ', but the only allowed values',ch10,&
&   'are 0 or 1 for multibinit.',ch10,&
&   'Action: correct fit_SPCoupling in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%fit_SPCoupling=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_SPCoupling',tread,'INT')
 if(tread==1) multibinit_dtset%fit_SPCoupling=intarr(1)
 if(multibinit_dtset%fit_SPCoupling<0.and.multibinit_dtset%fit_SPCoupling>1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
     &   'fit_SPCoupling is',multibinit_dtset%fit_SPCoupling,&
     &   ', but the only allowed values',ch10,&
&   'are 0 or 1 for multibinit.',ch10,&
&   'Action: correct fit_SPCoupling in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%bound_cutoff=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bound_cutoff',tread,'DPR')
 if(tread==1) multibinit_dtset%bound_cutoff=dprarr(1)
 if(multibinit_dtset%bound_cutoff<0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'bound_cutoff is',multibinit_dtset%bound_cutoff,', but the only allowed values',ch10,&
&   'are positives for multibinit.',ch10,&
&   'Action: correct bound_cutoff in your input file.'
   MSG_ERROR(message)
 end if


  multibinit_dtset%bound_maxCoeff=4
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bound_maxCoeff',tread,'INT')
 if(tread==1) multibinit_dtset%bound_maxCoeff=intarr(1)
 if(multibinit_dtset%bound_maxCoeff<0.and.multibinit_dtset%bound_maxCoeff>1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'bound_maxCoeff is',multibinit_dtset%bound_maxCoeff,', but the only allowed values',ch10,&
&   'are 0 or 1 for multibinit.',ch10,&
&   'Action: correct bound_maxCoeff in your input file.'
   MSG_ERROR(message)
 end if

  multibinit_dtset%bound_temp=325
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bound_temp',tread,'DPR')
 if(tread==1) multibinit_dtset%bound_temp=dprarr(1)
 if(multibinit_dtset%bound_temp<=0)then
   write(message, '(a,f10.1,a,a,a,a,a)' )&
&   'Bound_Temp is ',multibinit_dtset%bound_temp,'. The only allowed values',ch10,&
&   'are positives values.',ch10,&
&   'Action: correct Bound_Temp in your input file.'
   MSG_ERROR(message)
 end if
 multibinit_dtset%bound_step=1000
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bound_step',tread,'INT')
 if(tread==1) multibinit_dtset%bound_step=intarr(1)
 if(multibinit_dtset%bound_step<0.and.multibinit_dtset%bound_step>1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'bound_step is',multibinit_dtset%bound_step,', but the only allowed values',ch10,&
&   'are 0 or 1 for multibinit.',ch10,&
&   'Action: correct bound_step in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%fit_coeff=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_coeff',tread,'INT')
 if(tread==1) multibinit_dtset%fit_coeff=intarr(1)
 if(multibinit_dtset%fit_coeff<0.and.multibinit_dtset%fit_coeff>1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'fit_coeff is',multibinit_dtset%fit_coeff,', but the only allowed values',ch10,&
&   'are 0 or 1 for multibinit.',ch10,&
&   'Action: correct fit_coeff in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%fit_cutoff=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_cutoff',tread,'DPR')
 if(tread==1) multibinit_dtset%fit_cutoff=dprarr(1)
 if(multibinit_dtset%fit_cutoff<0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'fit_cutoff is',multibinit_dtset%fit_cutoff,', but the only allowed values',ch10,&
&   'are positives for multibinit.',ch10,&
&   'Action: correct fit_cutoff in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%fit_grid(:)= 1
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'fit_grid',tread,'INT')
 if(tread==1) multibinit_dtset%fit_grid(1:3)=intarr(1:3)
 do ii=1,3
   if(multibinit_dtset%fit_grid(ii)<0.or.multibinit_dtset%fit_grid(ii)>20)then
     write(message, '(a,i0,a,i0,a,a,a,i0,a)' )&
&     'fit_grid(',ii,') is ',multibinit_dtset%fit_grid(ii),', which is lower',&
&     ' than 0 of superior than 20.',&
&     ch10,'Action: correct fit_grid(',ii,') in your input file.'
     MSG_ERROR(message)
   end if
 end do

 multibinit_dtset%fit_rangePower(:)= (/3,4/)
 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'fit_rangePower',tread,'INT')
 if(tread==1) multibinit_dtset%fit_rangePower(1:2)=intarr(1:2)
 do ii=1,2
   if(multibinit_dtset%fit_rangePower(ii)<0.or.multibinit_dtset%fit_rangePower(ii)>20)then
     write(message, '(a,i0,a,i0,a,a,a,i0,a)' )&
&     'fit_rangePower(',ii,') is ',multibinit_dtset%fit_rangePower(ii),', which is lower',&
&     ' than 0 of superior than 20.',&
&     ch10,'Action: correct fit_rangePower(',ii,') in your input file.'
     MSG_ERROR(message)
   end if
 end do

 multibinit_dtset%bound_rangePower(:)= (/6,6/)
 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'bound_rangePower',tread,'INT')
 if(tread==1) multibinit_dtset%bound_rangePower(1:2)=intarr(1:2)
 do ii=1,2
   if(multibinit_dtset%bound_rangePower(ii)<=0.or.multibinit_dtset%bound_rangePower(ii)>20)then
     write(message, '(a,i0,a,i0,a,a,a,i0,a)' )&
&     'bound_rangePower(',ii,') is ',multibinit_dtset%bound_rangePower(ii),', which is lower',&
&     ' than 0 of superior than 20.',&
&     ch10,'Action: correct bound_rangePower(',ii,') in your input file.'
     MSG_ERROR(message)
   end if
 end do

  multibinit_dtset%bound_cell(:)= (/6,6,6/)
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'bound_cell',tread,'INT')
 if(tread==1) multibinit_dtset%bound_cell(1:3)=intarr(1:3)
 do ii=1,3
   if(multibinit_dtset%bound_cell(ii)<=0.or.multibinit_dtset%bound_cell(ii)>20)then
     write(message, '(a,i0,a,i0,4a,i0,a)' )&
&     'bound_cell(',ii,') is ',multibinit_dtset%bound_cell(ii),', which is lower',&
&     ' than 0 of superior than 20.',&
&     ch10,'Action: correct bound_cell(',ii,') in your input file.'
     MSG_ERROR(message)
   end if
 end do

 multibinit_dtset%fit_tolMSDE=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_tolMSDE',tread,'DPR')
 if(tread==1) multibinit_dtset%fit_tolMSDE=dprarr(1)
 if(multibinit_dtset%fit_tolMSDE<0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'fit_tolMSDE is',multibinit_dtset%fit_tolMSDE,', but the only allowed values',ch10,&
&   'are positives for multibinit.',ch10,&
&   'Action: correct fit_tolMSDE in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%fit_tolMSDF=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_tolMSDF',tread,'DPR')
 if(tread==1) multibinit_dtset%fit_tolMSDF=dprarr(1)
 if(multibinit_dtset%fit_tolMSDF<0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'fit_tolMSDF is',multibinit_dtset%fit_tolMSDF,', but the only allowed values',ch10,&
&   'are positives for multibinit.',ch10,&
&   'Action: correct fit_tolMSDF in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%fit_tolMSDS=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_tolMSDS',tread,'DPR')
 if(tread==1) multibinit_dtset%fit_tolMSDS=dprarr(1)
 if(multibinit_dtset%fit_tolMSDS<0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'fit_tolMSDS is',multibinit_dtset%fit_tolMSDS,', but the only allowed values',ch10,&
&   'are positives for multibinit.',ch10,&
&   'Action: correct fit_tolMSDS in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%fit_tolMSDFS=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_tolMSDFS',tread,'DPR')
 if(tread==1) multibinit_dtset%fit_tolMSDFS=dprarr(1)
 if(multibinit_dtset%fit_tolMSDFS<0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'fit_tolMSDFS is',multibinit_dtset%fit_tolMSDFS,', but the only allowed values',ch10,&
&   'are positives for multibinit.',ch10,&
&   'Action: correct fit_tolMSDFS in your input file.'
   MSG_ERROR(message)
 end if

!G

!H

!I


!J

!K

!L
 multibinit_dtset%latt_mask(:)=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'latt_mask',tread,'INT')
 if(tread==1) then
    do ii=1, 3
       multibinit_dtset%latt_mask(ii)=intarr(ii)
       if(multibinit_dtset%latt_mask(ii) <0 .or. multibinit_dtset%latt_mask(ii) >1)then
          write(message, '(a)' )&
               &   ' latt_mask element should be 0 or 1.'
          MSG_ERROR(message)
       end if
    end do
 end if


!M
!N

 ABI_ALLOCATE(multibinit_dtset%fit_bancoeff,(multibinit_dtset%fit_nbancoeff))
 if (multibinit_dtset%fit_nbancoeff >0)then
   if(multibinit_dtset%fit_nbancoeff>marr)then
     marr=multibinit_dtset%fit_nbancoeff
     ABI_DEALLOCATE(intarr)
     ABI_ALLOCATE(intarr,(marr))
   end if
   multibinit_dtset%fit_bancoeff(:)=0
   call intagm(dprarr,intarr,jdtset,marr,multibinit_dtset%fit_nbancoeff,&
&              string(1:lenstr),'fit_bancoeff',tread,'INT')
   if(tread==1)then
     do ii=1,multibinit_dtset%fit_nbancoeff
       multibinit_dtset%fit_bancoeff(ii)=intarr(ii)
     end do
   end if
 end if

 ABI_ALLOCATE(multibinit_dtset%fit_fixcoeff,(multibinit_dtset%fit_nfixcoeff))
 if (multibinit_dtset%fit_nfixcoeff >0)then
   if(multibinit_dtset%fit_nfixcoeff>marr)then
     marr=multibinit_dtset%fit_nfixcoeff
     ABI_DEALLOCATE(intarr)
     ABI_ALLOCATE(intarr,(marr))
   end if
   multibinit_dtset%fit_fixcoeff(:)=0
   call intagm(dprarr,intarr,jdtset,marr,multibinit_dtset%fit_nfixcoeff,&
&              string(1:lenstr),'fit_fixcoeff',tread,'INT')
   if(tread==1)then
     do ii=1,multibinit_dtset%fit_nfixcoeff
       multibinit_dtset%fit_fixcoeff(ii)=intarr(ii)
     end do
   end if
 end if

!O

!P

!Q
 ABI_ALLOCATE(multibinit_dtset%qmass,(multibinit_dtset%nnos))
 multibinit_dtset%qmass(:)= zero
 if(multibinit_dtset%nnos>=1)then
   if(multibinit_dtset%nnos>marr)then
     marr=multibinit_dtset%nnos
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   call intagm(dprarr,intarr,jdtset,marr,multibinit_dtset%nnos,string(1:lenstr),'qmass',tread,'DPR')
   if(tread==1) multibinit_dtset%qmass(:)=dprarr(1:multibinit_dtset%nnos)
 end if

 if (multibinit_dtset%nqshft/=0)then
   if(3*multibinit_dtset%nqshft>marr)then
     marr=3*multibinit_dtset%nqshft
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   ABI_ALLOCATE(multibinit_dtset%q1shft,(3,multibinit_dtset%nqshft))
   multibinit_dtset%q1shft(:,:)=zero
   call intagm(dprarr,intarr,jdtset,marr,3*multibinit_dtset%nqshft, string(1:lenstr),'q1shft',tread,'DPR')
   if(tread==1) multibinit_dtset%q1shft(1:3,1:multibinit_dtset%nqshft)=&
&   reshape(dprarr(1:3*multibinit_dtset%nqshft),(/3,multibinit_dtset%nqshft/))
 end if

 ABI_ALLOCATE(multibinit_dtset%qph1l,(3,multibinit_dtset%nph1l))
 ABI_ALLOCATE(multibinit_dtset%qnrml1,(multibinit_dtset%nph1l))
 if (multibinit_dtset%nph1l/=0)then
   if(4*multibinit_dtset%nph1l>marr)then
     marr=4*multibinit_dtset%nph1l
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   multibinit_dtset%qph1l(:,:)=zero
   multibinit_dtset%qnrml1(:)=zero
   call intagm(dprarr,intarr,jdtset,marr,4*multibinit_dtset%nph1l,string(1:lenstr),'qph1l',tread,'DPR')
   if(tread==1)then
     do iph1=1,multibinit_dtset%nph1l
       do ii=1,3
         multibinit_dtset%qph1l(ii,iph1)=dprarr(ii+(iph1-1)*4)
       end do
       multibinit_dtset%qnrml1(iph1)=dprarr(4+(iph1-1)*4)
       if(abs(multibinit_dtset%qnrml1(iph1))<DDB_QTOL)then
         write(message, '(a,a,a,a,a)' )&
&         'The first list of wavevectors ','should not have non-analytical data.',ch10,&
&         'Action: correct the first list',' of wavevectors in the input file.'
         MSG_ERROR(message)
       end if
     end do
   end if
 end if

 ABI_ALLOCATE(multibinit_dtset%qph2l,(3,multibinit_dtset%nph2l))
 ABI_ALLOCATE(multibinit_dtset%qnrml2,(multibinit_dtset%nph2l))
 if (multibinit_dtset%nph2l/=0)then
   if(4*multibinit_dtset%nph2l>marr)then
     marr=4*multibinit_dtset%nph2l
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   multibinit_dtset%qph2l(:,:)=zero
   multibinit_dtset%qnrml2(:)=zero
   call intagm(dprarr,intarr,jdtset,marr,4*multibinit_dtset%nph2l,string(1:lenstr),'qph2l',tread,'DPR')
   if(tread==1)then
     do iph2=1,multibinit_dtset%nph2l
       do ii=1,3
         multibinit_dtset%qph2l(ii,iph2)=dprarr(ii+(iph2-1)*4)
       end do
       multibinit_dtset%qnrml2(iph2)=dprarr(4+(iph2-1)*4)
       if(abs(multibinit_dtset%qnrml2(iph2))>DDB_QTOL)then
         write(message, '(a,a,a,a,a)' )&
&         'The second list of wavevectors',' should have only non-analytical data.',ch10,&
&         'Action: correct the second list','of wavevectors in the input file.'
         MSG_ERROR(message)
       end if
     end do
   end if
 end if

!R
 if(9>marr)then
   marr=9
   ABI_DEALLOCATE(intarr)
   ABI_DEALLOCATE(dprarr)
   ABI_ALLOCATE(intarr,(marr))
   ABI_ALLOCATE(dprarr,(marr))
 end if
 multibinit_dtset%rprim(:,:)= zero
 call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),'rprim',tread,'DPR')
 if(tread==1) then
   multibinit_dtset%rprim(1:3,1:3)= reshape(dprarr(1:9),(/3,3/))
! check new rprimd
   if(all(abs(multibinit_dtset%rprim(1,:))<tol16).or.&
&    all(abs(multibinit_dtset%rprim(2,:))<tol16).or.all(abs(multibinit_dtset%rprim(3,:))<tol16)) then
     write(message, '(3a)' )&
&  ' There is a problem with rprim',ch10,&
&   'Action: correct rprim'
     MSG_BUG(message)
   end if
 end if
!S

 if(6>marr)then
   marr=6
   ABI_DEALLOCATE(intarr)
   ABI_DEALLOCATE(dprarr)
   ABI_ALLOCATE(intarr,(marr))
   ABI_ALLOCATE(dprarr,(marr))
 end if
 multibinit_dtset%strten_reference(:)= zero
 call intagm(dprarr,intarr,jdtset,marr,6,string(1:lenstr),'strten_reference',tread,'DPR')
 if(tread==1) multibinit_dtset%strten_reference(1:6)= dprarr(1:6)

!T

!U

!V

!W

!X

!Y

!Z

!=======================================================================
!Finished reading in variables - deallocate
!=======================================================================

 ABI_DEALLOCATE(dprarr)
 ABI_DEALLOCATE(intarr)

!=======================================================================
!Check consistency of input variables:
!=======================================================================

 if(multibinit_dtset%prtsrlr/=0 .and. multibinit_dtset%ifcflag/=1) then
   write(message, '(3a)' )&
&   'ifcflag must be 1 for the SR/LR decomposition of the phonon frequencies',ch10,&
&   'Action: correct ifcflag in your input file.'
   MSG_ERROR(message)
 end if

!FIXME: add check that if freeze_displ /= 0 then you need to be doing ifc and phonon interpolation

 if (multibinit_dtset%ifcflag > 0 .and. sum(abs(multibinit_dtset%ngqpt)) == 0) then
   write(message, '(3a)' )&
&   'if you want interatomic force constant output, multibinit needs ngqpt input variable ',ch10,&
&   'Action: set ngqpt in your input file.'
   MSG_ERROR(message)
 end if

!check that q-grid refinement is a divisor of ngqpt in each direction
 if(any(multibinit_dtset%qrefine(:) > 1) .and. &
&    any(abs(dmod(dble(multibinit_dtset%ngqpt(1:3))/dble(multibinit_dtset%qrefine(:)),one)) > tol10)) then
   write(message, '(a,3i0,a,a,a,3i8,a,a)' )&
&   'qrefine is',multibinit_dtset%qrefine,' The only allowed values',ch10,&
&   'are integers which are divisors of the ngqpt grid', multibinit_dtset%ngqpt,ch10,&
&   'Action: correct qrefine in your input file.'
   MSG_ERROR(message)
 end if

! check new rprimd
 if(all(multibinit_dtset%acell(:) > one).and.all(abs(multibinit_dtset%rprim(:,:))<tol16))then
   write(message, '(3a)' )&
&         ' acell is defined but there is no rprim',ch10,&
&         'Action: add rprim input'
   MSG_BUG(message)
 end if


!check the fit_bancoeff and fit_fixcoeff
 do ii=1,multibinit_dtset%fit_nbancoeff
   do jj=ii+1,multibinit_dtset%fit_nbancoeff
     if (multibinit_dtset%fit_bancoeff(ii) == multibinit_dtset%fit_bancoeff(jj))then
       write(message, '(a,I0,a,I0,2a)' )&
&           ' There is two similar numbers for fit_bancoeff: ',multibinit_dtset%fit_bancoeff(ii),&
&           ' and ', multibinit_dtset%fit_bancoeff(jj),ch10,&
&            'Action: change fit_bancoeff'
       MSG_BUG(message)
     end if
   end do
 end do

 do ii=1,multibinit_dtset%fit_nfixcoeff
   do jj=ii+1,multibinit_dtset%fit_nfixcoeff
     if (multibinit_dtset%fit_fixcoeff(ii) == multibinit_dtset%fit_fixcoeff(jj))then
       write(message, '(a,I0,a,I0,2a)' )&
&           ' There is two similar numbers for fit_fixcoeff: ',multibinit_dtset%fit_fixcoeff(ii),&
&           ' and ', multibinit_dtset%fit_fixcoeff(jj),ch10,&
&            'Action: change fit_fixcoeff'
       MSG_BUG(message)
     end if
   end do
 end do

 do ii=1,3
   if(multibinit_dtset%dipdip_range(ii) < multibinit_dtset%ncell(ii)) then
     write(message,'(4a,3I3,3a,3I3,6a)') ch10,&
&                 ' --- !WARNING',ch10,&
&                 '     The range of dipdip_range (',multibinit_dtset%dipdip_range(:),')',ch10,&
&                 '     But the range of the cell for the simulation is',&
&                       multibinit_dtset%ncell(:),')',ch10,&
&                 '     dipdip_range is set to ncell.',ch10,&
&                 ' ---',ch10
     multibinit_dtset%dipdip_range(:) =  multibinit_dtset%ncell(:)
     call wrtout(std_out,message,'COLL')
     exit
   end if
   if(multibinit_dtset%dipdip_range(ii) < multibinit_dtset%bound_cell(ii)) then
     write(message,'(4a,3I3,3a,3I3,6a)') ch10,&
&                 ' --- !WARNING',ch10,&
&                 '     The range of dipdip_range (',multibinit_dtset%dipdip_range(:),')',ch10,&
&                 '     But the range of the cell for the simulation is',&
&                       multibinit_dtset%ncell(:),')',ch10,&
&                 '     dipdip_range is set to bound_cell.',ch10,&
&                 ' ---',ch10
     multibinit_dtset%dipdip_range(:) =  multibinit_dtset%ncell(:)
     call wrtout(std_out,message,'COLL')
     exit
   end if
 end do

!Check if only one tolerance is specify
 if(abs(multibinit_dtset%fit_tolMSDF) >zero .and. abs(multibinit_dtset%fit_tolMSDS) >zero) then
   write(message, '(3a)' ) &
&           ' There is two tolerance flags for the fit: fit_tolMSDF and fit_tolMSDS',ch10,&
&            'Action: Put only one tolerance flag'
   MSG_BUG(message)
 end if
 if(abs(multibinit_dtset%fit_tolMSDF) >zero .and. abs(multibinit_dtset%fit_tolMSDE) >zero)then
   write(message, '(3a)' ) &
&           ' There is two tolerance flags for the fit: fit_tolMSDF and fit_tolMSDE',ch10,&
&            'Action: Put only one tolerance flag'
   MSG_BUG(message)
 end if
 if(abs(multibinit_dtset%fit_tolMSDF) >zero .and. abs(multibinit_dtset%fit_tolMSDFS) >zero)then
   write(message, '(3a)' ) &
&           ' There is two tolerance flags for the fit: fit_tolMSDF and fit_tolMSDFS',ch10,&
&            'Action: Put only one tolerance flag'
   MSG_BUG(message)
 end if
 if(abs(multibinit_dtset%fit_tolMSDS) >zero .and. abs(multibinit_dtset%fit_tolMSDE) >zero)then
   write(message, '(3a)' ) &
&           ' There is two tolerance flags for the fit: fit_tolMSDS and fit_tolMSDE',ch10,&
&            'Action: Put only one tolerance flag'
   MSG_BUG(message)
 end if
 if(abs(multibinit_dtset%fit_tolMSDS) >zero .and. abs(multibinit_dtset%fit_tolMSDFS) >zero)then
   write(message, '(3a)' ) &
&           ' There is two tolerance flags for the fit: fit_tolMSDS and fit_tolMSDFS',ch10,&
&            'Action: Put only one tolerance flag'
   MSG_BUG(message)
 end if
 if(abs(multibinit_dtset%fit_tolMSDE) >zero .and. abs(multibinit_dtset%fit_tolMSDFS) >zero)then
   write(message, '(3a)' ) &
&           ' There is two tolerance flags for the fit: fit_tolMSDE and fit_tolMSDFS',ch10,&
&            'Action: Put only one tolerance flag'
   MSG_BUG(message)
 end if

end subroutine invars10
!!***



!----------------------------------------------------------------------

!!****f* m_multibinit_dataset/outvars_multibinit
!!
!! NAME
!! outvars_multibinit
!!
!! FUNCTION
!! Open input file for the multibinit code, then
!! echoes the input information.
!!
!! INPUTS
!! multibinit_dtset <type(multibinit_dtset_type)> datatype with all the input variables
!! nunit=unit number for input or output
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! PARENTS
!!      multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine outvars_multibinit (multibinit_dtset,nunit)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nunit
 type(multibinit_dtset_type),intent(in) :: multibinit_dtset

!Local variables -------------------------
!Set routine version number here:
!scalars
 integer :: ii,iph1,iph2,iqshft

!*********************************************************************

!Write the heading
 write(nunit,'(a,80a,a)') ch10,('=',ii=1,80),ch10
 write(nunit, '(a,a)' )&
& ' -outvars_multibinit: echo values of input variables ----------------------',ch10

!The flags
 if(multibinit_dtset%ifcflag/=0)then
   write(nunit,'(a)')' Flags : '
   if(multibinit_dtset%ifcflag/=0)write(nunit,'(3x,a9,3i10)')'  ifcflag',multibinit_dtset%ifcflag
   if(multibinit_dtset%prt_model/=0)write(nunit,'(3x,a9,3i10)')'prt_model',multibinit_dtset%prt_model
   if(multibinit_dtset%prt_phfrq/=0)write(nunit,'(3x,a9,3i10)')'prt_phfrq',multibinit_dtset%prt_phfrq
   if(multibinit_dtset%strcpling/=0)write(nunit,'(3x,a9,3i10)')'  strcpling',multibinit_dtset%strcpling
   if(multibinit_dtset%strcpling==2)write(nunit,'(3x,a9,3es8.2)')'delta_df',multibinit_dtset%delta_df
 end if

 if(multibinit_dtset%dynamics/=0)then
   write(nunit,'(a)')' Molecular Dynamics :'
   write(nunit,'(3x,a9,3I10.1)')' dynamics',multibinit_dtset%dynamics
   write(nunit,'(3x,a9,3F10.1)')'     temp',multibinit_dtset%temperature
   write(nunit,'(3x,a9,3I10.1)')'    ntime',multibinit_dtset%ntime
   if (multibinit_dtset%nctime /=1)then
     write(nunit,'(3x,a9,3I10.1)')'   nctime',multibinit_dtset%nctime
   end if
   write(nunit,'(3x,a9,3i10)')  '    ncell',multibinit_dtset%ncell
   write(nunit,'(3x,a9,3i10)')  '    dtion',multibinit_dtset%dtion
   if (multibinit_dtset%restartxf/=0) then
     write(nunit,'(3x,a9,3i10)')  'restartxf',multibinit_dtset%restartxf
   end if
   if(multibinit_dtset%dynamics==13)then
     write(nunit,'(3x,a9,3i10)')'  optcell',multibinit_dtset%optcell
     write(nunit,'(3x,a9,3F12.1)')'    bmass',multibinit_dtset%bmass
     write(nunit,'(3x,a9,3I10)')'     nnos',multibinit_dtset%nnos
     write(nunit,'(3x,a12)',advance='no')'    qmass  '
     write(nunit,'(3x,15F12.10)') (multibinit_dtset%qmass(ii),ii=1,multibinit_dtset%nnos)
   end if

   if(multibinit_dtset%dynamics==101)then
   end if

   if(multibinit_dtset%dynamics==102)then
      write(nunit,'(a15,ES15.5)')'latt_friction',multibinit_dtset%latt_friction
   end if

   if(multibinit_dtset%dynamics==103)then
      write(nunit,'(a15,ES15.5)')'     latt_taut',multibinit_dtset%latt_taut
   end if

   if(multibinit_dtset%dynamics==104)then
      write(nunit,'(a15,ES15.5)')'     latt_taut',multibinit_dtset%latt_taut
      write(nunit,'(a15,ES15.5)')'     latt_taup',multibinit_dtset%latt_taup
      write(nunit,'(a15,ES15.5)')'compressibility',multibinit_dtset%latt_compressibility
   end if
   
   if(multibinit_dtset%dynamics==105)then
      write(nunit,'(a15,ES15.5)')'     latt_taut',multibinit_dtset%latt_taut
      write(nunit,'(a15,ES15.5)')'     latt_taup',multibinit_dtset%latt_taup
      write(nunit,'(a15,ES15.5)')'compressibility',multibinit_dtset%latt_compressibility
      write(nunit,'(a15,ES15.5)')'     latt_mask',(multibinit_dtset%latt_mask(ii), ii=1, 3)
   end if

 end if

 if(multibinit_dtset%spin_dynamics/=0) then
    write(nunit,'(a)')' Spin Dynamics :'

    write(nunit,'(3x,a25,I12.1)')'spin_calc_correlation_obs',multibinit_dtset%spin_calc_correlation_obs
    write(nunit,'(3x,a25,I12.1)')'spin_calc_thermo_obs',multibinit_dtset%spin_calc_thermo_obs
    write(nunit,'(3x,a25,I12.1)')'spin_calc_traj_obs',multibinit_dtset%spin_calc_traj_obs
    write(nunit,'(12x,a16,I12.1)')'spin_dynamics',multibinit_dtset%spin_dynamics
    write(nunit,'(10x, a18, 5x, F10.5)')'spin_temperature',multibinit_dtset%spin_temperature
    write(nunit,'(10x, a18, 5x, F10.5)')'spin_damping',multibinit_dtset%spin_damping
    write(nunit,'(9x,a19,I10.1)')'spin_ntime_pre',multibinit_dtset%spin_ntime_pre
    write(nunit,'(13x,a15,I10.1)')'spin_ntime',multibinit_dtset%spin_ntime
    write(nunit,'(13x,a15,3I10)')  'ncell',multibinit_dtset%ncell !TODO hexu: duplicate but dynamics can be 0.
    write(nunit,'(13x,a15,ES15.5, a8)')  'spin_dt',multibinit_dtset%spin_dt*Time_Sec , ' second' !TODO: use a.u.
    !write(nunit,'(3x,a14,3es10.5)')  '   spin_tolavg',multibinit_dtset%spin_tolavg
    !write(nunit,'(3x,a14,3es10.5)')  '   spin_tolvar',multibinit_dtset%spin_tolvar
    write(nunit,'(13x,a15)')   'spin_mag_field'
    write(nunit,'(31x,3es12.5, a8)')   (multibinit_dtset%spin_mag_field(ii)/Bfield_Tesla,ii=1,3), '   Tesla'
    write(nunit, '(13x, a15, I12.1)') 'spin_sia_add', multibinit_dtset%spin_sia_add
    write(nunit, '(13x, a15, ES15.5)') 'spin_sia_k1amp', multibinit_dtset%spin_sia_k1amp
    write(nunit, '(13x, a15, 3ES15.5)') 'spin_sia_k1dir', (multibinit_dtset%spin_sia_k1dir(ii), ii=1,3)
    write(nunit,'(13x,a15)')   'spin_qpoint'
    write(nunit,'(28x,3es12.5)')   (multibinit_dtset%spin_qpoint(ii),ii=1,3)
    write(nunit, '(13x, a15, I12.1)') 'spin_init_state', multibinit_dtset%spin_init_state
    if(multibinit_dtset%spin_init_state==4) then
      write(nunit,'(13x,a25)')   'spin_init_orientation'
      write(nunit,'(28x,3es12.5)')   (multibinit_dtset%spin_init_orientation(ii),ii=1,3)
      write(nunit,'(13x,a18)')   'spin_init_qpoint'
      write(nunit,'(28x,3es12.5)')   (multibinit_dtset%spin_init_qpoint(ii),ii=1,3)
      write(nunit,'(13x,a25)')   'spin_init_rotate_axis'
      write(nunit,'(28x,3es12.5)')   (multibinit_dtset%spin_init_rotate_axis(ii),ii=1,3)
    endif
    write(nunit, '(6x, a22, I12.1)') 'spin_var_temperature', multibinit_dtset%spin_var_temperature
    write(nunit, '(6x, a22, 5x, F10.5)') 'spin_temperature_start', multibinit_dtset%spin_temperature_start
    write(nunit, '(6x, a22, 5x, F10.5)') 'spin_temperature_end', multibinit_dtset%spin_temperature_end
    write(nunit, '(6x, a22, 5x, I12.1)') 'spin_temperature_nstep', multibinit_dtset%spin_temperature_nstep
    write(nunit, '(13x, a15, I12.1)') 'spin_write_traj', multibinit_dtset%spin_write_traj
 end if

 if(multibinit_dtset%slc_coupling/=0) then
   write(nunit,'(a)')' Spin-Lattice coupling :'
 endif

 if(multibinit_dtset%confinement==1)then
   write(nunit,'(a)')' Confinement information :'
   write(nunit,'(1x,a22,I5.1)')'       conf_power_disp',multibinit_dtset%conf_power_disp
   write(nunit,'(1x,a22,I5.1)')'     conf_power_strain',multibinit_dtset%conf_power_strain
   write(nunit,'(1x,a22,3es16.8)')'  conf_power_fact_disp',multibinit_dtset%conf_power_fact_disp
   write(nunit,'(1x,a22,3es16.8)')'conf_power_fact_strain',multibinit_dtset%conf_power_fact_strain
   write(nunit,'(1x,a22)')'     conf_cutoff_disp'
   write(nunit,'(19x,3es16.8)') (multibinit_dtset%conf_cutoff_disp(ii),ii=1,multibinit_dtset%natom)
   write(nunit,'(1x,a22)')'    conf_cutoff_strain'
   write(nunit,'(19x,3es16.8)') (multibinit_dtset%conf_cutoff_strain(ii),ii=1,6)
 end if

 if(multibinit_dtset%fit_coeff/=0)then
   write(nunit,'(a)')' Fit the coefficients :'
   write(nunit,'(1x,a17,I3.1)')'       fit_coeff',multibinit_dtset%fit_coeff
   write(nunit,'(1x,a17,I3.1)')'fit_generateCoeff',multibinit_dtset%fit_generateCoeff
   if(multibinit_dtset%fit_initializeData==0)then
     write(nunit,'(1x,a17,I3.1)')'fit_initializeData',multibinit_dtset%fit_initializeData
   end if
   if(multibinit_dtset%fit_tolMSDE > 0)then
     write(nunit,'(1x,a17,es16.8)')'    fit_tolMSDE',multibinit_dtset%fit_tolMSDE
   end if
   if(multibinit_dtset%fit_tolMSDF > 0)then
     write(nunit,'(1x,a17,es16.8)')'    fit_tolMSDF',multibinit_dtset%fit_tolMSDF
   end if
   if(multibinit_dtset%fit_tolMSDS > 0)then
     write(nunit,'(1x,a17,es16.8)')'    fit_tolMSDS',multibinit_dtset%fit_tolMSDS
   end if
   if(multibinit_dtset%fit_tolMSDFS > 0)then
     write(nunit,'(1x,a17,es16.8)')'   fit_tolMSDFS',multibinit_dtset%fit_tolMSDFS
   end if
   write(nunit,'(1x,a17,es16.8)')'      fit_cutoff',multibinit_dtset%fit_cutoff
   write(nunit,'(1x,a17,I3.1)')'      fit_option',multibinit_dtset%fit_option
   write(nunit,'(1x,a17,2x,I0)')'      fit_ncoeff',multibinit_dtset%fit_ncoeff   
   write(nunit,'(1x,a17,3i3)') '        fit_grid',multibinit_dtset%fit_grid
   write(nunit,'(1x,a17,I3.1)')'   ts_option',multibinit_dtset%ts_option
   write(nunit,'(1x,a17,2i3)') '  fit_rangePower',multibinit_dtset%fit_rangePower
   write(nunit,'(1x,a17,I3)')  '  fit_anhaStrain',multibinit_dtset%fit_anhaStrain
   write(nunit,'(1x,a17,I3)')  '  fit_SPCoupling',multibinit_dtset%fit_SPCoupling
   if(multibinit_dtset%fit_nbancoeff /= 0) then
     write(nunit,'(1x,a17,I3)')  '   fit_nbancoeff',multibinit_dtset%fit_nbancoeff
     write(nunit,'(1x,a17)',advance='no')'   fit_bancoeff'
     write(nunit,'(4x,9i7)') (multibinit_dtset%fit_bancoeff(ii),ii=1,multibinit_dtset%fit_nbancoeff)
   end if
   if(multibinit_dtset%fit_nfixcoeff /= 0) then    
     write(nunit,'(1x,a17,I3)')  '   fit_nfixcoeff',multibinit_dtset%fit_nfixcoeff
     write(nunit,'(1x,a17)',advance='no')'   fit_fixcoeff'
     write(nunit,'(4x,9i7)') (multibinit_dtset%fit_fixcoeff(ii),ii=1,multibinit_dtset%fit_nfixcoeff)
   end if
 end if

 if(multibinit_dtset%bound_model /=0)then
   write(nunit,'(a)')' Bound the coefficients :'
   write(nunit,'(1x,a16,I3)')    'bound_anhaStrain',multibinit_dtset%bound_anhaStrain   
   write(nunit,'(1x,a16,I3)')    'bound_SPCoupling',multibinit_dtset%bound_SPCoupling
   write(nunit,'(1x,a16,es16.8)')'    bound_cutoff',multibinit_dtset%bound_cutoff
   write(nunit,'(1x,a16,1x,3I3)')   '      bound_cell',multibinit_dtset%bound_cell
   write(nunit,'(1x,a16,1x,I3)')    '  bound_maxCoeff',multibinit_dtset%bound_maxCoeff
   write(nunit,'(1x,a16,es16.8)') '      bound_temp',multibinit_dtset%bound_temp
   write(nunit,'(1x,a16,I7)')   '      bound_step',multibinit_dtset%bound_step
   write(nunit,'(1x,a16,2I3.1)')'bound_rangePower',multibinit_dtset%bound_rangePower
 end if

!Write the general information
 if( multibinit_dtset%rfmeth/=1 .or. &
& multibinit_dtset%enunit/=0 .or. &
& multibinit_dtset%eivec/=0 .or. &
& multibinit_dtset%asr/=0 .or. &
& multibinit_dtset%chneut/=0)then
   write(nunit,'(a)')' Miscellaneous information :'
   if(multibinit_dtset%rfmeth/=1)write(nunit,'(3x,a9,3i10)')'   rfmeth',multibinit_dtset%rfmeth
   if(multibinit_dtset%enunit/=0)write(nunit,'(3x,a9,3i10)')'   enunit',multibinit_dtset%enunit
   if(multibinit_dtset%eivec/=0) write(nunit,'(3x,a9,3i10)')'    eivec',multibinit_dtset%eivec
   if(multibinit_dtset%asr/=0)   write(nunit,'(3x,a9,3i10)')'      asr',multibinit_dtset%asr
   if(multibinit_dtset%chneut/=0)write(nunit,'(3x,a9,3i10)')'   chneut',multibinit_dtset%chneut
 end if


!For interatomic force constant information
 if(multibinit_dtset%ifcflag/=0)then
   write(nunit,'(a)')' Interatomic Force Constants Inputs :'
   write(nunit,'(3x,a9,3i10)')'   dipdip',multibinit_dtset%dipdip
   if(multibinit_dtset%dipdip /= 0)then
     write(nunit,'(a12,3i10)') 'dipdip_range',multibinit_dtset%dipdip_range
   end if
   if(multibinit_dtset%dipdip_prt/=0)then
     write(nunit,'(a12,3i10)') 'dipdip_prt',multibinit_dtset%dipdip_prt
   end if
   if(multibinit_dtset%nsphere/=0)write(nunit,'(3x,a9,3i10)')'  nsphere',multibinit_dtset%nsphere
   if(abs(multibinit_dtset%rifcsph)>tol10)write(nunit,'(3x,a9,E16.6)')'  nsphere',multibinit_dtset%rifcsph
   write(nunit,'(3x,a9,3i10)')'   ifcana',multibinit_dtset%ifcana
   write(nunit,'(3x,a9,3i10)')'   ifcout',multibinit_dtset%ifcout
   if(multibinit_dtset%natifc>=1)then
     write(nunit,'(3x,a9,3i10)')'   natifc',multibinit_dtset%natifc
     write(nunit,'(3x,a12)',advance='no')'    atifc   '
     write(nunit,'(3x,15i4)') (multibinit_dtset%atifc(ii)*ii,ii=1,multibinit_dtset%natifc)

   end if
   write(nunit,'(a)')' Description of grid 1 :'
   write(nunit,'(3x,a9,3i10)')'     brav',multibinit_dtset%brav
   write(nunit,'(3x,a9,3i10)')'    ngqpt',multibinit_dtset%ngqpt(1:3)
   write(nunit,'(3x,a9,3i10)')'   nqshft',multibinit_dtset%nqshft
   if (multibinit_dtset%nqshft/=0)then
     write(nunit,'(3x,a9)')'   q1shft'
     do iqshft=1,multibinit_dtset%nqshft
       write(nunit,'(19x,4es16.8)') (multibinit_dtset%q1shft(ii,iqshft),ii=1,3)
     end do
   end if
   if (any(multibinit_dtset%qrefine(:) > 1)) then
     write(nunit,'(3x,a9,3i10)')'  qrefine', multibinit_dtset%qrefine
   end if
 end if


!List of vector 1  (reduced coordinates)
 if(multibinit_dtset%nph1l/=0)then
   write(nunit,'(a)')' First list of wavevector (reduced coord.) :'
   write(nunit,'(3x,a9,3i10)')'    nph1l',multibinit_dtset%nph1l
   write(nunit,'(3x,a9)')'    qph1l'
   do iph1=1,multibinit_dtset%nph1l
     write(nunit,'(19x,3es16.8,2x,es11.3)') &
&     (multibinit_dtset%qph1l(ii,iph1),ii=1,3),multibinit_dtset%qnrml1(iph1)
   end do
 end if

!List of vector 2  (cartesian coordinates)
 if(multibinit_dtset%nph2l/=0)then
   write(nunit,'(a)')' Second list of wavevector (cart. coord.) :'
   write(nunit,'(3x,a9,3i10)')'    nph2l',multibinit_dtset%nph2l
   write(nunit,'(3x,a9)')'    qph2l'
   do iph2=1,multibinit_dtset%nph2l
     write(nunit,'(19x,3es16.8,2x,es11.3)') &
&     (multibinit_dtset%qph2l(ii,iph2),ii=1,3),multibinit_dtset%qnrml2(iph2)
   end do
 end if

 write(nunit,'(a,80a,a)') ch10,('=',ii=1,80),ch10

end subroutine outvars_multibinit
!!***

end module m_multibinit_dataset
!!***
