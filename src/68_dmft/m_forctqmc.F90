!!****m* ABINIT/m_forctqmc
!! NAME
!!  m_forctqmc
!!
!! FUNCTION
!! Prepare CTQMC and call CTQMC
!!
!! COPYRIGHT
!! Copyright (C) 2006-2025 ABINIT group (BAmadon, VPlanes)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_forctqmc

 use defs_basis
 use m_abicore
 use m_Ctqmc
 use m_CtqmcInterface
 use m_Ctqmcoffdiag
 use m_CtqmcoffdiagInterface
 use m_data4entropyDMFT
 use m_errors
 use m_GreenHyb

 use m_crystal, only : crystal_t
 use m_datafordmft, only : compute_levels,hybridization_asymptotic_coefficient
 use m_energy, only : compute_migdal_energy,compute_trace_log_loc
 use m_fstrings, only : int2char4
 use m_green, only : compute_moments_loc,copy_green,destroy_green,green_type, &
    & init_green,int_fct,occup_green_tau,print_green
 use m_hide_lapack, only : matrginv,xginv
 use m_hu, only : copy_hu,destroy_hu,destroy_vee,hu_type,init_vee, &
    & rotatevee_hu,vee_type,vee_ndim2tndim_hu_r
 use m_io_tools, only : flush_unit,open_file
 use m_matlu, only : add_matlu,checkdiag_matlu,checkreal_matlu,chi_matlu,copy_matlu,destroy_matlu, &
     & diag_matlu,diff_matlu,fac_matlu,gather_matlu,init_matlu,magmomforb_matlu,magmomfspin_matlu, &
     & magmomfzeeman_matlu,matlu_type,print_matlu,printplot_matlu,prod_matlu,rotate_matlu,shift_matlu, &
     & slm2ylm_matlu,sym_matlu,symmetrize_matlu,xmpi_matlu,ylm2jmj_matlu,zero_matlu,magnfield_matlu
 use m_numeric_tools, only : coeffs_gausslegint
 use m_oper, only : destroy_oper,gather_oper,identity_oper,init_oper,inverse_oper,oper_type
 use m_paw_correlations, only : calc_vee
 use m_paw_dmft, only : paw_dmft_type
 use m_paw_numeric, only : jbessel => paw_jbessel
 use m_pawang, only : pawang_type
 use m_self, only : destroy_self,initialize_self,self_type
 use m_special_funcs, only : sbf8
 use m_splines, only : spline2_complex

 use netcdf !If calling TRIQS via python invocation, write a .nc file

 implicit none

 private

 public :: qmc_prep_ctqmc
 public :: testcode_ctqmc
 public :: testcode_ctqmc_b
 public :: ctqmcoutput_to_green
 public :: ctqmcoutput_printgreen
 public :: ctqmc_calltriqs
 public :: ctqmc_calltriqs_c
!!***

contains
!!****f* m_forctqmc/qmc_prep_ctqmc
!! NAME
!! qmc_prep_ctqmc
!!
!! FUNCTION
!! Prepare and call the qmc subroutines
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  self <type(self_type)>= self-energy
!!  hu <type(hu_type)>= U interaction
!!  paw_dmft <type(paw_dmft_type)>= DMFT data structure
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawprtvol = drive the amount of writed data.
!!  weiss <type(green_type)>= weiss function
!!
!! OUTPUT
!!  green <type(green_type)>= green function
!!
!! NOTES
!!
!! SOURCE

subroutine qmc_prep_ctqmc(cryst_struc,green,self,hu,paw_dmft,pawang,pawprtvol,weiss)

!Arguments ------------------------------------
!scalars
! type(pawang_type), intent(in) :: pawang
 type(crystal_t), intent(in) :: cryst_struc
 type(green_type), intent(inout) :: green  ! MGNAG: This fix the problem with v7[27:29] on nag@petrus
 type(hu_type), intent(inout) :: hu(cryst_struc%ntypat)
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(pawang_type), intent(in) :: pawang
 integer, intent(in) :: pawprtvol
 type(green_type), target, intent(inout) :: weiss
 type(self_type), intent(inout) :: self
!Local variables ------------------------------
 integer :: iatom,icomp,ierr,if1,if2,iflavor1,iflavor2,ifreq,im1,im2,ima,imb,ispa,ispb,ispinor
 integer :: ispinor1,ispinor2,isppol,itau,itypat,lpawu,myproc,natom,ndim,nflavor,nomega,nproc
 integer :: nspinor,nsppol,nsppol_imp,ntypat,nwlo,opt_diag,opt_fk,opt_nondiag
 integer :: opt_rot,rot_type_vee,testcode,testrot,tndim,unt,unt2,useylm
 integer, parameter :: optdb = 0
 logical :: nondiaglevels
 logical(kind=1) :: leg_measure = .true.
 real(dp) :: doccsum,EE,f4of2_sla,f6of2_sla,noise,omega
 type(green_type) :: weiss_for_rot
 type(oper_type) :: energy_level,level_diag
 type(CtqmcInterface) :: hybrid
 type(CtqmcoffdiagInterface) :: hybridoffdiag
 real(dp) :: umod(2,2)
 complex(dp) :: integral(2,2)
 real(dp), allocatable :: docc(:,:),gtmp(:,:),gtmp_nd(:,:,:),levels_ctqmc(:),vee(:,:,:,:)
 complex(dp), allocatable :: muorb,muspin,muzeem
 complex(dp), allocatable :: fw1(:,:),fw1_nd(:,:,:),gw_tmp(:,:),gw_tmp_nd(:,:,:)
 complex(dp), allocatable :: gw1_nd(:,:,:),hybri_limit(:,:),levels_ctqmc_nd(:,:),shift(:)
 type(coeff2c_type), allocatable :: magmom_orb(:),magmom_spin(:),magmom_tot(:)
 type(hu_type), allocatable :: hu_for_s(:)
 type(matlu_type), allocatable :: dmat_diag(:),eigvectmatlu(:),hybri_coeff(:),matlu1(:),matlu2(:),matlu3(:)
 type(matlu_type), allocatable :: matlu4(:),matlumag(:),matlumag_orb(:),matlumag_spin(:),matlumag_tot(:)
 type(matlu_type), allocatable :: udens_atoms(:),udens_atoms_for_s(:)
 type(matlu_type), allocatable :: levels_temp(:),magnfield(:)
 type(vee_type), allocatable :: vee_for_s(:),vee_rotated(:)
 character(len=13) :: tag
 character(len=500) :: message
! ************************************************************************

 !mbandc=paw_dmft%mbandc
 !nkpt=paw_dmft%nkpt
 natom   = paw_dmft%natom
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol
 ntypat  = paw_dmft%ntypat
 nwlo    = paw_dmft%dmft_nwlo
 !greendft%whichgreen="DFT"

 call init_green(weiss_for_rot,paw_dmft,opt_oper_ksloc=2)
 ! call init_green(gw_loc,paw_dmft)
 call copy_green(weiss,weiss_for_rot,opt_tw=2)
 !=======================================================================
 !== Use one QMC solver   ===============================================
 !=======================================================================
 write(message,'(3a)') ch10,'  ===  CT-QMC solver === ',ch10
 call wrtout(std_out,message,'COLL')

 ! Initialise for compiler
 ! omega_current=czero

 ! Initialise nproc
 nproc=paw_dmft%nproc
 myproc = paw_dmft%myproc

 ! ======================================
 ! Allocations: diagonalization and eigenvectors
 ! ======================================
 ABI_MALLOC(udens_atoms,(natom))
 ABI_MALLOC(eigvectmatlu,(natom))
 ABI_MALLOC(magmom_orb,(natom))
 ABI_MALLOC(matlumag_orb,(natom))
 ABI_MALLOC(magmom_spin,(natom))
 ABI_MALLOC(matlumag_spin,(natom))
 ABI_MALLOC(magmom_tot,(natom))
 ABI_MALLOC(matlumag_tot,(natom))
 if (paw_dmft%ientropy == 1) then
   ABI_MALLOC(udens_atoms_for_s,(natom))
 end if
 ABI_MALLOC(dmat_diag,(natom))
 call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),dmat_diag(:))
 call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),eigvectmatlu(:))
 call init_matlu(natom,2,1,paw_dmft%lpawu(:),udens_atoms(:))
 if (paw_dmft%ientropy == 1) then
   call init_matlu(natom,2,1,paw_dmft%lpawu(:),udens_atoms_for_s(:))
 end if
 do iatom=1,natom
   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   ABI_MALLOC(magmom_orb(iatom)%value,(2*(2*lpawu+1),2*(2*lpawu+1)))
   magmom_orb(iatom)%value=czero
   ABI_MALLOC(magmom_spin(iatom)%value,(2*(2*lpawu+1),2*(2*lpawu+1)))
   magmom_spin(iatom)%value=czero
   ABI_MALLOC(magmom_tot(iatom)%value,(2*(2*lpawu+1),2*(2*lpawu+1)))
   magmom_tot(iatom)%value=czero
 end do ! iatom

 ! =================================================================
 ! Impose diago of density matrix
 ! =================================================================

 ! =================================================================
 ! Impose diago of levels and Ylm basis if opt_nondiag=1
 ! =================================================================
 ! opt_diag=1 ! 1: diago the levels (The best choice).
 ! opt_diag=2 ! 2: diago density matrix (can be used for historical reasons)

 !  Need in the general case of two input variable for opt_diag and
 !  opt_nondiag!
 !  opt_nondiag should be 0 by default
 opt_diag    = 1
 opt_nondiag = 0
 if (paw_dmft%dmft_solv >= 6) opt_nondiag = 1 ! Use ctqmc in abinit with offdiag terms in F
 !else
 !  opt_nondiag = 0 ! use fast ctqmc in ABINIT without off diagonal terms in F
 !end if

 useylm = 0
 if (nspinor == 2) useylm = 1 ! to avoid complex G(tau)

 !write(6,*) "nspinor,useylm",nspinor,useylm
 if (useylm == 0) then
   write(std_out,*) " Slm (real spherical harmonics) basis is used (before a possible rotation)"
   rot_type_vee = 1 ! for rotatevee_hu
 else if (useylm == 1) then
   write(std_out,*) " Ylm (complex spherical harmonics) basis is used (before rotation)"
   rot_type_vee = 4 ! for rotatevee_hu
 end if ! useylm

 ! if(useylm==1.and.opt_diag/=1) ABI_ERROR("useylm==1 and opt_diag/=0 is not possible")
 if (hu(1)%jpawu_zero .and. nsppol == 2) nsppol_imp = 2 ! J=0 and nsppol=2
 if (.not. hu(1)%jpawu_zero .or. nsppol /= 2) nsppol_imp = 1  ! J/=0 ou nsppol=1
 ! =================================================================
 ! Compute DFT Green's function to compare to weiss_for_rot (check)
 ! =================================================================
 ! call init_green(greendft,paw_dmft,opt_oper_ksloc=3)
 ! call greendftcompute_green(cryst_struc,greendft,pawang,paw_dmft)
 !! call copy_green(greendft,weiss_for_rot,2)

 ! =================================================================
 ! Compute atomic levels
 ! =================================================================
 call init_oper(paw_dmft,energy_level,opt_ksloc=2)

 ! ----------------------------------
 ! Compute atomic levels in Slm basis
 ! ----------------------------------
 call compute_levels(energy_level,self%hdc,paw_dmft,nondiag=nondiaglevels)

 ! ------------------------------------------------
 ! If levels are not diagonal, then diagonalize it (according to
 ! dmftctqmc_basis)
 ! ------------------------------------------------
 if (paw_dmft%dmftctqmc_basis == 1) then
   if (nondiaglevels .or. useylm == 1) then
     opt_diag = 1
     write(message,'(3a)') ch10,"   == Hamiltonian in local basis is not diagonal: diagonalize it",ch10
   else
     opt_diag = 0
     write(message,'(5a)') ch10,"   == Hamiltonian in local basis is diagonal in the Slm basis ",ch10, &
        & "      CTQMC will use this basis",ch10
   end if ! nondiaglevels or useylm
 else if (paw_dmft%dmftctqmc_basis == 2) then
   if (nondiaglevels .or. useylm == 1) then
     write(message,'(7a)') ch10,"   == Hamiltonian in local basis is not diagonal",ch10, &
       & "   == According to dmftctqmc_basis: diagonalize density matrix",ch10, &
       & "   == Warning : Check that the Hamiltonian is diagonal !",ch10
     opt_diag = 2
   else
     write(message,'(5a)') ch10,"   == Hamiltonian in local basis is diagonal in the Slm basis ",ch10, &
        & "      CTQMC will use this basis",ch10
     opt_diag = 0
   end if ! nondiaglevels or useylm
 else if (paw_dmft%dmftctqmc_basis == 0) then
   if (nondiaglevels) then
     write(message,'(4a)') ch10,"   == Hamiltonian in local basis is not diagonal",ch10, &
       & "   == According to dmftctqmc_basis: keep this non diagonal basis for the calculation"
   else
     write(message,'(5a)') ch10,"   == Hamiltonian in local basis is diagonal in the Slm basis ",ch10, &
       & "      CTQMC will use this basis",ch10
   end if ! nondiaglevels
   opt_diag = 0
 end if ! dmftctqmc_basis
 call wrtout(std_out,message,'COLL')
 if (opt_diag == 1) then
   write(std_out,*) "  ==  The atomic levels are diagonalized"
 else if (opt_diag == 2) then
   write(std_out,*) "  ==  The correlated occupation matrix is diagonalized"
 end if ! opt_diag

 ! =================================================================
 ! Now, check if diagonalisation is necessary
 ! =================================================================


 ! =================================================================
 ! First rotate to Ylm basis the atomic levels
 ! =================================================================

 if (useylm == 1) then

   ! Rotate from Slm to Ylm the atomic levels
   ! ----------------------------------------
   call slm2ylm_matlu(energy_level%matlu(:),natom,paw_dmft,1,pawprtvol)

   ! Print atomic energy levels in Ylm basis
   ! --------------------------------
   if (pawprtvol >= 3) then
     write(message,'(2a)') ch10," == Print Energy levels in Ylm basis"
     call wrtout(std_out,message,'COLL')
     call print_matlu(energy_level%matlu(:),natom,1)
   end if ! pawprtvol>=3

   !==================================================================
   ! Add Zeeman contributions to local energy levels when nspinor = 2
   !==================================================================
   if(paw_dmft%dmft_magnfield .eq. 2 .and. nspinor .eq. 2) then

     ABI_MALLOC(magnfield,(natom))
     ABI_MALLOC(levels_temp,(natom))

     write(message,'(a,2x,a)') ch10, " == Add Zeeman contributions to local energy levels in Ylm"
     call wrtout(std_out,message,'COLL')

     call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,magnfield)
     call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,levels_temp)
     call copy_matlu(energy_level%matlu,levels_temp,natom)

     !Spin-Orbit case, not much tested so far (need to remove AFM sym)
     call magnfield_matlu(magnfield,natom,paw_dmft%dmft_magnfield_b,2)
     !call print_matlu(magnfield,natom,1)
     call add_matlu(levels_temp,magnfield,energy_level%matlu,natom,-1)
     call print_matlu(energy_level%matlu,natom,1)
     call destroy_matlu(magnfield,natom)
     call destroy_matlu(levels_temp,natom)

     ABI_FREE(magnfield)
     ABI_FREE(levels_temp)
   endif !dmft_magnfield
 end if ! useylm

 ABI_MALLOC(vee_rotated,(natom))
 call init_vee(paw_dmft,vee_rotated(:))

 ! ===========================================================================================
 ! Start for diagonalization of levels/density matrix according to opt_diag
 ! ===========================================================================================

 !opt_rot=2 ! do it one time before CTQMC
 opt_rot = 1 ! do all the rotations successively on all different quantities.
 if (opt_diag == 1 .or. opt_diag == 0) then

   if (opt_diag == 1) then
     ! =================================================================
     ! Diagonalize atomic levels
     ! =================================================================
     call init_oper(paw_dmft,level_diag,opt_ksloc=2)

     ! Diagonalize atomic levels (opt_real is necessary, because
     ! rotation must be real in order for the occupations and Green's
     ! function to be real)
     ! ---------------------------------------------------------------
     call diag_matlu(energy_level%matlu(:),level_diag%matlu(:),natom,pawprtvol,eigvectmatlu(:),&
                   & nsppol_imp=nsppol_imp,opt_real=1,test=paw_dmft%dmft_solv)  ! temporary: test should be extended to all cases.

 !     call rotate_matlu(energy_level%matlu,eigvectmatlu,natom,3,1)
 !       write(message,'(a,2x,a,f13.5)') ch10,&
 !&       " == Print first Diagonalized Energy levels for Fermi Level=",paw_dmft%fermie
 !       call wrtout(std_out,message,'COLL')
 !       call print_matlu(energy_level%matlu,natom,1,compl=1,opt_exp=1)

     if (opt_rot == 1) then
       call copy_matlu(level_diag%matlu(:),energy_level%matlu(:),natom)
     end if

     call destroy_oper(level_diag)

     ! Print diagonalized levels
     ! --------------------------
     write(tag,'(f13.5)') paw_dmft%fermie
     if (pawprtvol >= 3) then
       write(message,'(a,2x,2a)') ch10,&
         & " == Print Diagonalized Energy levels for Fermi Level=",adjustl(tag)
       call wrtout(std_out,message,'COLL')
       call print_matlu(energy_level%matlu(:),natom,1,compl=1,opt_exp=1)
     else
       write(message,'(a,2x,2a)') ch10,&
         & " == Energy levels Diagonalized for Fermi Level=",adjustl(tag)
       call wrtout(std_out,message,'COLL')
     end if ! pawprtvol>=3

   else if (opt_diag == 0) then
     do iatom = 1, natom
       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       itypat = paw_dmft%typat(iatom)
       !  write(6,*) size(udens_atoms(iatom)%value)
       !  write(6,*) size(hu(itypat)%udens)
       !  write(6,*) udens_atoms(iatom)%value
       !  write(6,*) hu(itypat)%udens
       udens_atoms(iatom)%mat(:,:,1)=hu(itypat)%udens(:,:)
       vee_rotated(iatom)%mat(:,:,:,:) = hu(itypat)%veeslm2(:,:,:,:)
     end do ! iatom
   end if ! opt_diag=0 or 1
 ! call rotatevee_hu(cryst_struc,hu,nspinor,nsppol,pawprtvol,eigvectmatlu,udens_atoms)

 else if (opt_diag == 2) then
   ! =================================================================
   ! Diagonalizes density matrix and keep eigenvectors in eigvectmatlu
   ! =================================================================

   ! Print density matrix before diagonalization
   ! -------------------------------------------
   if (pawprtvol >= 3) then
     write(message,'(a,2x,a)') ch10," == Density Matrix before diagonalization ="
     call wrtout(std_out,message,'COLL')
     !MGNAG: This call is wrong if green has intent(out), now we use intent(inout)
     call print_matlu(green%occup%matlu(:),natom,1)
   end if ! pawprtvol>=3

   !!  checkstop: we can have two different diagonalisation basis for the up and dn
   !!  but one use the same basis, unless the error is really to large(>0.1)

   ! Diagonalize density matrix
   ! ---------------------------
   call diag_matlu(green%occup%matlu(:),dmat_diag(:),natom,4,eigvectmatlu(:), &
                 & nsppol_imp=nsppol_imp,checkstop=.false.)

   ! Print diagonalized density matrix
   ! ----------------------------------
   if (pawprtvol >= 3) then
     write(message,'(a,2x,a)') ch10,&
       & " == Diagonalized Density Matrix in the basis used for QMC ="
     call wrtout(std_out,message,'COLL')
     call print_matlu(dmat_diag(:),natom,1)

     !write(message,'(2a,i3,13x,a)') ch10,'    ==  Rotation of interaction matrix =='
     !call wrtout(std_out,message,'COLL')
   end if ! pawprtvol>=3

   !if (.not.hu(1)%jpawu_zero) &
   !ABI_WARNING("In qmc_prep_ctqmc J/=0 and rotation matrix not rotated")
   !  Rotate interaction.
   !   call rotatevee_hu(cryst_struc,hu,nspinor,nsppol,pawprtvol,eigvectmatlu,udens_atoms)
   !   call rotatevee_hu(cryst_struc,hu,nspinor,nsppol,pawprtvol,eigvectmatlu,udens_atoms,rot_type_vee)
 end if
 ! ===========================================================================================
 ! END Of diagonalization
 ! ===========================================================================================

 if (paw_dmft%ientropy == 1) then
   ABI_MALLOC(hu_for_s,(ntypat))
       ! Usefull to compute interaction energy for U=1 J=J/U when U=0.
   call copy_hu(ntypat,hu(:),hu_for_s(:))
   f4of2_sla = - one
   f6of2_sla = - one
   do itypat=1,ntypat
     ndim = 2*hu(itypat)%lpawu + 1
     ABI_MALLOC(vee,(ndim,ndim,ndim,ndim))
     call calc_vee(f4of2_sla,f6of2_sla,paw_dmft%j_for_s/paw_dmft%u_for_s, &
                 & hu_for_s(itypat)%lpawu,pawang,one,vee(:,:,:,:))
     hu_for_s(itypat)%vee(:,:,:,:) = cmplx(vee(:,:,:,:),zero,kind=dp)
     ABI_FREE(vee)
   end do
   ABI_MALLOC(vee_for_s,(natom))
   call init_vee(paw_dmft,vee_for_s(:))
   call rotatevee_hu(hu_for_s(:),paw_dmft,pawprtvol,eigvectmatlu(:), &
                   & rot_type_vee,udens_atoms_for_s(:),vee_for_s(:))
   call destroy_hu(hu_for_s(:),ntypat)
!      udens_atoms_for_s will be used later.
   ABI_FREE(hu_for_s)
   call destroy_vee(paw_dmft,vee_for_s(:))
   ABI_FREE(vee_for_s)
 end if ! ientropy=1

 call flush_unit(std_out)

 ! ===========================================================================================
 ! Broadcast matrix of rotation from processor 0 to the other
 ! In case of degenerate levels, severals rotations are possible. Here we
 ! choose the rotation of proc 0. It is arbitrary.
 ! ===========================================================================================
 call xmpi_matlu(eigvectmatlu(:),natom,paw_dmft%spacecomm,master=0,option=2)

 if (opt_diag /= 0) then
   call rotatevee_hu(hu(:),paw_dmft,pawprtvol,eigvectmatlu(:), &
                   & rot_type_vee,udens_atoms(:),vee_rotated(:))
 end if

 !unitnb=300000+paw_dmft%myproc
 !call int2char4(paw_dmft%myproc,tag_proc)
 !tmpfil = 'eigvectmatluaftermpi'//tag_proc
 !open (unit=unitnb,file=trim(tmpfil),status='unknown',form='formatted')
 !do iflavor1=1,14
 !  do iflavor2=1,14
 !    write(unitnb,*) iflavor1,iflavor2,eigvectmatlu(1,1)%value(iflavor1,iflavor2)
 !  enddo
 !enddo

 ! ===========================================================================================
 ! Now rotate various quantities in the new basis
 ! ===========================================================================================

 !=======================================================
 ! Allocate, Compute, and Rotate atomic levels for CTQMC
 !=======================================================

 ! If levels not rotated, rotate them
 ! -----------------------------------
 if (opt_diag == 2 .and. opt_rot == 1) then
   call rotate_matlu(energy_level%matlu(:),eigvectmatlu(:),natom,1)
 end if

 ! Print atomic levels
 ! -------------------
 if (pawprtvol >= 3 .and. opt_diag == 2 .and. opt_rot == 1) then
   write(message,'(a,2x,a)') ch10," == Print Energy levels in CTQMC basis"
   call wrtout(std_out,message,'COLL')
   call print_matlu(energy_level%matlu(:),natom,1)
 else if (opt_diag == 2 .and. opt_rot == 1) then
   write(message,'(a,2x,a)') ch10," == CT-QMC Energy levels rotated"
   call wrtout(std_out,message,'COLL')
 end if ! pawprtvol>=3

 !====================================================================
 ! If levels were diagonalized before, then rotate density matrix for
 ! information.
 !====================================================================
 if (opt_diag == 1) then

   ABI_MALLOC(matlu1,(natom))
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),matlu1(:))
   call copy_matlu(green%occup%matlu(:),matlu1(:),natom)
   if (pawprtvol >= 3) then
     write(message,'(a,2x,a)') ch10," == Occupations before rotations"
     call wrtout(std_out,message,'COLL')
     call print_matlu(green%occup%matlu(:),natom,1)
   end if ! pawprtvol>=3

   ! 1) rotate density matrix to Ylm basis
   ! --------------------------------------
   if (useylm == 1) then
     call slm2ylm_matlu(matlu1(:),natom,paw_dmft,1,pawprtvol)
     if (pawprtvol >= 3) then
       write(message,'(2a)') ch10," == Print occupations in Ylm basis"
       call wrtout(std_out,message,'COLL')
       call print_matlu(matlu1(:),natom,1)
     end if
   end if ! useylm

   ! 2) rotate density matrix to rotated basis
   ! -------------------------------------------
   if (opt_rot == 1 .or. opt_rot == 2) then
     call rotate_matlu(matlu1(:),eigvectmatlu(:),natom,1)
   end if
   write(message,'(a,2x,a)') ch10," == Rotated occupations (for information)"
   call wrtout(std_out,message,'COLL')
   call print_matlu(matlu1(:),natom,1,compl=1)
   call checkreal_matlu(matlu1(:),natom,tol10)
   call destroy_matlu(matlu1(:),natom)
   ABI_FREE(matlu1)

 end if ! opt_diag=1

 call flush_unit(std_out)

 ! =================================================================
 ! Rotate weiss function according to eigenvectors.
 ! =================================================================
 !!!stop
 ! Rotate Weiss function first in Ylm basis
 ! -----------------------------------------------------------------
 if(useylm==1) then
   write(message,'(a,2x,a)') ch10, " == Rotation of weiss and greendft in the Ylm Basis="
   call wrtout(std_out,message,'COLL')
   do ifreq=1,nwlo
     call slm2ylm_matlu(weiss_for_rot%oper(ifreq)%matlu(:),natom,paw_dmft,1,0)
     call slm2ylm_matlu(weiss%oper(ifreq)%matlu(:),natom,paw_dmft,1,0)
     ! call slm2ylm_matlu(greendft%oper(ifreq)%matlu,natom,1,0)
   end do
 end if

 if (pawprtvol >= 3) then
   !   write(message,'(a,2x,a,f13.5)') ch10,& ! debug
   !   " == Print weiss for small freq 1 before rot" ! debug
   !   call wrtout(std_out,message,'COLL') ! debug
   !   call print_matlu(weiss_for_rot%oper(1)%matlu,natom,1) !  debug

   ! Print Weiss function
   ! --------------------
   write(message,'(a,2x,a)') ch10," == Print weiss for 1st freq before rot" ! debug
   call wrtout(std_out,message,'COLL') ! debug
   call print_matlu(weiss_for_rot%oper(1)%matlu(:),natom,1,compl=1) !  debug
   write(message,'(a,2x,a)') ch10," == Print weiss for last freq before rot" ! debug
   call wrtout(std_out,message,'COLL') ! debug
   call print_matlu(weiss_for_rot%oper(nwlo)%matlu(:),natom,1,compl=1) !  debug
!    write(message,'(a,2x,a,f13.5)') ch10,& ! debug
!&   " == Print DFT G for 1st freq before rot" ! debug
!    call wrtout(std_out,message,'COLL') ! debug
!    call print_matlu(greendft%oper(1)%matlu,natom,1,compl=1,opt_exp=2) !  debug
!    write(message,'(a,2x,a,f13.5)') ch10,& ! debug
!&   " == Print DFT G for last freq before rot" ! debug
!    call wrtout(std_out,message,'COLL') ! debug
!    call print_matlu(greendft%oper(paw_dmft%dmft_nwlo)%matlu,natom,1,compl=1,opt_exp=2) !  debug
 end if ! pawprtvol>=3

 if (opt_diag /= 0) then
   ! Rotate Weiss function from the Slm (or Ylm) to the basis of diagonalisation
   ! -------------------------------------------------------------------
   write(message,'(a,2x,a)') ch10, " == Rotation of weiss ="
   call wrtout(std_out,message,'COLL')

   do ifreq=1,nwlo
     if (opt_rot == 1) then
       call rotate_matlu(weiss_for_rot%oper(ifreq)%matlu(:),eigvectmatlu(:),natom,1)
       call rotate_matlu(weiss%oper(ifreq)%matlu(:),eigvectmatlu(:),natom,1)
     end if ! opt_rot=1
!    call checkdiag_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,tol6)
   end do ! ifreq

   if (myproc == mod(nproc+1,nproc)) then
     if (open_file(trim(paw_dmft%filapp)//"_atom_G0w_.dat",message,newunit=unt) /= 0) ABI_ERROR(message)
     ndim = 2*paw_dmft%lpawu(natom) + 1
     do ifreq=1,nwlo
       write(unt,'(29f21.14)') paw_dmft%omega_lo(ifreq),&
         & (((weiss_for_rot%oper(ifreq)%matlu(natom)%mat(im1+(ispinor-1)*ndim,im1+(ispinor-1)*ndim,isppol),&
         & im1=1,3),ispinor=1,nspinor),isppol=1,nsppol)
     end do ! ifreq
     close(unt)
   end if ! myproc=master

   call flush_unit(std_out)
   if (pawprtvol >= 3) then
     write(message,'(a,2x,a)') ch10," == Print weiss for small freq 1 after rot" ! debug
     call wrtout(std_out,message,'COLL') ! debug
     call print_matlu(weiss_for_rot%oper(1)%matlu(:),natom,1,compl=1) !  debug
     write(message,'(a,2x,a)') ch10," == Print weiss for last freq after rot"   ! debug
     call wrtout(std_out,message,'COLL')   ! debug
     call print_matlu(weiss_for_rot%oper(nwlo)%matlu(:),natom,1,compl=1) ! debug
   end if ! pawprtvol>=3

   !   ! Rotate DFT Green's function first in Ylm basis then in the rotated basis and compare to weiss_for_rot
   !   ! -----------------------------------------------------------------------------------------------------
   !   write(message,'(a,2x,a)') ch10, " == Rotation of greendft ="
   !   call wrtout(std_out,message,'COLL')
   !   do ifreq=1,paw_dmft%dmft_nwlo
   !     if(opt_rot==1) call rotate_matlu(greendft%oper(ifreq)%matlu,eigvectmatlu,natom,3,1)
   !     call diff_matlu("Weiss_for_rot","greendft",weiss_for_rot%oper(ifreq)%matlu,greendft%oper(ifreq)%matlu,natom,1,tol14)
   !!    call checkdiag_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,tol6)
   !   end do
   !   if(pawprtvol>=3) then
   !     write(message,'(a,2x,a,f13.5)') ch10,& ! debug
   !&    " == Print greendft for small freq 1 after rot" ! debug
   !     call wrtout(std_out,message,'COLL') ! debug
   !     call print_matlu(greendft%oper(1)%matlu,natom,1,compl=1,opt_exp=2) !  debug
   !     write(message,'(a,2x,a,f13.5)') ch10,&   ! debug
   !&    " == Print greendft for last freq after rot"   ! debug
   !     call wrtout(std_out,message,'COLL')   ! debug
   !     call print_matlu(greendft%oper(paw_dmft%dmft_nwlo)%matlu,natom,1,compl=1,opt_exp=2) ! debug
   !   end if
   !   call flush_unit(std_out)
 end if

 ! =================================================================
 ! Compute analytic limit of hybridization and rotate it
 ! =================================================================
 ABI_MALLOC(hybri_coeff,(natom))
 call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),hybri_coeff(:))
 !write(6,*)"hybri1",hybri_coeff(1)%mat(1,1,1,1,1),paw_dmft%natom,cryst_struc%natom

 ! Compute analytical C_ij such that F_ij -> C_ij/iw_n
 ! ---------------------------------------
 call hybridization_asymptotic_coefficient(cryst_struc,paw_dmft,hybri_coeff(:))
 write(message,'(a,2x,a)') ch10," == Coeff analytical C_ij such that F -> C_ij/iw_n for large frequency"
 call wrtout(std_out,message,'COLL')

 ! Print analytical C_ij (not rotated)
 ! ---------------------------------------
 call print_matlu(hybri_coeff(:),natom,1)

 ! Rotate analytical C_ij in Ylm basis
 ! ---------------------------------------
 if (useylm == 1) then
   call slm2ylm_matlu(hybri_coeff(:),natom,paw_dmft,1,pawprtvol)
 end if

 if (opt_diag /= 0)  then
   ! Rotate analytical C_ij in rotated basis
   ! ---------------------------------------
   if (opt_rot == 1 .or. opt_rot == 2) then
     call rotate_matlu(hybri_coeff(:),eigvectmatlu(:),natom,1)
   end if

   ! Print analytical C_ij (rotated)
   ! ---------------------------------------
   write(message,'(a,2x,a)') ch10," == Coeff analytical C_ij such that F -> C_ij/iw_n after rotation"
   call wrtout(std_out,message,'COLL')
   call print_matlu(hybri_coeff(:),natom,1,compl=1,opt_exp=1)
 end if

 ! =================================================================
 ! Check if rotation is properly done.
 ! =================================================================
 if(3 == 4) then
   write(message,'(a,2x,a)') ch10, " == Print  dmat before rot"
   call wrtout(std_out,message,'COLL')
   call print_matlu(green%occup%matlu(:),natom,1)
   if (useylm == 1) then
     call slm2ylm_matlu(green%occup%matlu(:),natom,paw_dmft,1,pawprtvol)
   end if
   if (opt_rot == 1) then
     call rotate_matlu(green%occup%matlu(:),eigvectmatlu(:),natom,1)
   end if
   write(message,'(a,2x,a)') ch10," == Print  dmat after rot"
   call wrtout(std_out,message,'COLL')
   call print_matlu(green%occup%matlu(:),natom,1)

   write(message,'(2a)') ch10,' QMC STOP: DEBUG'
   call wrtout(std_out,message,'COLL')
   ABI_ERROR(message)
 end if
 ! =================================================================
 ! Check
 ! =================================================================

 ! write(message,'(a,2x,a,f13.5)') ch10,&
 !&   " == Print weiss for small tau"
 ! call wrtout(std_out,message,'COLL')
 ! call print_matlu(weiss%oper(1)%matlu,natom,1)
 ! write(message,'(a,2x,a,f13.5)') ch10,&
 !&   " == Print weiss for large tau"
 ! call wrtout(std_out,message,'COLL')
 ! call print_matlu(weiss%oper(paw_dmft%dmft_nwlo)%matlu,natom,1)
 ! call flush_unit(std_out)
 ! write(message,'(2a)') ch10,' Check weiss_for_rot(last freq)'
 ! call wrtout(std_out,message,'COLL')
 ! call checkdiag_matlu(weiss_for_rot%oper(paw_dmft%dmft_nwlo)%matlu,natom,tol6,opt=nspinor)
 ! call flush_unit(std_out)
 ! write(message,'(2a)') ch10,' Check weiss_for_rot(ifreq=1)'
 ! call wrtout(std_out,message,'COLL')
 ! call checkdiag_matlu(weiss_for_rot%oper(1)%matlu,natom,tol6,opt=nspinor)
 ! call flush_unit(std_out)

 !master = 0

 ! =================================================================
 ! Print out
 ! =================================================================

! Print Weiss
! -------------
 if (paw_dmft%dmft_prgn == 1) then
   call print_green('Weiss_diag',weiss_for_rot,1,paw_dmft,opt_wt=1,opt_decim=1)
 end if

 write(message,'(a,2x,a)') ch10," == Preparing data for CTQMC"
 call wrtout(std_out,message,'COLL')

 ! Print Rotate Weiss for 1st and last frequencies
 ! ------------------------------------------------
 if (pawprtvol >= 3) then
   write(message,'(a,2x,a)') ch10," == Print rotated weiss function for small freq in the rotated basis"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(weiss_for_rot%oper(1)%matlu(:),natom,1,compl=1)  ! debug
   write(message,'(a,2x,a)') ch10," == Print rotated weiss function for largest freq in the rotated basis"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(weiss_for_rot%oper(nwlo)%matlu(:),natom,1,compl=1)  ! debug
 end if ! pawprtvol>=3

 ! =================================================================
 !  VARIABLES FOR CTQMC TESTS
 testcode = 0
 testrot = 0
 ! opt_fk = 0 ! for developpers to check Fourier transform and computes G0(tau)
 opt_fk = 1 ! usual case: for real calculations
 ! =================================================================

 ! _________________________________________________________________
 !
 !  SECOND PART : BUILT HYBRIDIZATION FROM G0
 ! _________________________________________________________________
 !
 ! =================================================================
 ! Compute inverse of weiss and compute hybridization
 ! =================================================================

 ! Compute inverse of weiss for each Frequency
 ! -----------------------------------------------------------------

  do ifreq=1,nwlo
    ABI_MALLOC(matlu1,(natom))
    ABI_MALLOC(matlu2,(natom))
    call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),matlu1(:))
    call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),matlu2(:))

    call copy_matlu(weiss_for_rot%oper(ifreq)%matlu(:),matlu1(:),natom)

    ! Print G_0(iw_n)
    ! ----------------
    if (optdb == 1) then
      call printplot_matlu(weiss_for_rot%oper(ifreq)%matlu(:),natom,paw_dmft%omega_lo(ifreq),"go",60000,imre=1)
    end if

    ! Compute G_0^-1
    ! -------------------------------------------
    ! if opt_fk=1 or testcode/=0  Do the inversion
    ! if opt_fk=0                 Do not inverse.
    ! If testcode=2 and opt_fk=0  Do the inversion
    ! If testcode=1 and opt_fk=0  Do the inversion but no effect, because it will nevertheless be erased
    ! If opt_fk=1                 Do the inversion
    ! -------------------------------------------
    if (optdb == 1) then
      call printplot_matlu(matlu1(:),natom,paw_dmft%omega_lo(ifreq),"weiss",12000,imre=1)
    end if
    if (opt_fk == 1 .or. testcode /= 0) then
      call inverse_oper(weiss_for_rot%oper(ifreq),2)
    end if

    ! Print G_0^-1(iw_n)
    ! ----------------
    if (optdb == 1) then
      call printplot_matlu(weiss_for_rot%oper(ifreq)%matlu(:),natom,paw_dmft%omega_lo(ifreq),"goinv",70000,imre=1)
    end if

    if (pawprtvol >= 4 .or. ifreq == nwlo) then
      if (opt_fk == 1 .or. testcode /= 0) then
        ! Check inversion : do the product
        ! ----------------------------------------------
        call prod_matlu(weiss_for_rot%oper(ifreq)%matlu(:),matlu1(:),matlu2(:),natom)
        write(message,'(a,2x,a,i7)') ch10," == Print product of  weiss times invers for freq",ifreq
        call wrtout(std_out,message,'COLL')  ! debug
        call print_matlu(matlu2(:),natom,1)  ! debug
      end if
    end if

    call destroy_matlu(matlu1(:),natom)
    call destroy_matlu(matlu2(:),natom)
    ABI_FREE(matlu1)
    ABI_FREE(matlu2)

  end do ! ifreq

 ! Copy weiss_for_rot into weiss
 ! -------------------------------
 !call copy_matlu(weiss_for_rot%oper(ifreq)%matlu,weiss%oper(ifreq)%matlu,natom)


 ! Print G_0^-1 for 1st and last frequencies.
 ! -----------------------------------------
 if (pawprtvol >= 3) then
   write(message,'(a,2x,a)') ch10," == Print G_0^-1 for small freq in the rotated basis"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(weiss_for_rot%oper(1)%matlu(:),natom,1)  ! debug
   write(message,'(a,2x,a,e18.10,a)') ch10,&   ! debug
     & " == Print G_0^-1 for last freq in the rotated basis (last freq=",paw_dmft%omega_lo(nwlo),")"  ! debug
   call wrtout(std_out,message,'COLL')   ! debug
   call print_matlu(weiss_for_rot%oper(paw_dmft%dmft_nwlo)%matlu(:),natom,1,compl=1) ! debug
 end if ! pawprtvol>=3

 ! Substract frequency from diagonal part
 ! ======================================

 ABI_MALLOC(shift,(natom))
 do ifreq=1,nwlo

   shift(:) = cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)

   !  write(5555,'(400e17.4)') paw_dmft%omega_lo(ifreq),((((((weiss_for_rot%oper(ifreq)%matlu(1)%mat&
   !  & (im,im1,isppol,ispinor,ispinor1)-cmplx(0.d0,paw_dmft%omega_lo(ifreq),kind=dp)),im=1,2*3+1),&
   !&      im1=1,2*3+1),isppol=1,nsppol),ispinor=1,nspinor),ispinor1=1,nspinor)

   if (opt_fk == 1) then
     ! Compute G_0^-1-iw_n
     ! --------------------
     call shift_matlu(weiss_for_rot%oper(ifreq)%matlu(:),natom,shift(:))

     ! Compute -G_0^-1+iw_n
     ! --------------------
     call fac_matlu(weiss_for_rot%oper(ifreq)%matlu(:),natom,-cone)
   end if

   ! Print -G_0^-1+iw_n
   ! --------------------
   if (optdb == 1) then
     call printplot_matlu(weiss_for_rot%oper(ifreq)%matlu(:),natom,paw_dmft%omega_lo(ifreq), &
                        & "G0inv_minus_omega",20000,imre=1)
   end if
 end do ! ifreq

 ! Print -G_0^+1-iw_n=(F-levels) for last freq in the rotated basis"
 ! ------------------------------------------------------------------
 ABI_FREE(shift)
 if (pawprtvol >= 3) then
   write(message,'(a,2x,a)') ch10,&  ! debug
     & " == Print G_0^-1-iw_n=-(F-levels) for last freq in the rotated basis"  ! debug
   call wrtout(std_out,message,'COLL')   ! debug
   call print_matlu(weiss_for_rot%oper(nwlo)%matlu(:),natom,1,compl=1) ! debug
 end if ! pawprtvol>=3

 ! Check numerical limit of F(i_wn)*iw_n (can be used also to compute F )
 ! ======================================

 if (opt_nondiag == 1) then

   ABI_MALLOC(matlu1,(natom))
   ABI_MALLOC(matlu2,(natom))
   ABI_MALLOC(matlu3,(natom))
   ABI_MALLOC(matlu4,(natom))
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),matlu1(:))
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),matlu2(:))
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),matlu3(:))
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),matlu4(:))

   write(message,'(a,2x,a)') ch10," == energy_levels"
   call wrtout(std_out,message,'COLL')
   call print_matlu(energy_level%matlu(:),natom,1,opt_exp=2,compl=1)

   do ifreq=nwlo,1,-1 ! necessary to have matlu4 computed for the max frequency and available for all frequency.
     ! do ifreq=paw_dmft%dmftqmc_l,1,-1 ! necessary to have matlu4 computed for the max frequency and available for all frequency.

     ! Compute F (substract levels) for max frequency
     ! -----------------------------------------------
     call add_matlu(weiss_for_rot%oper(ifreq)%matlu(:),energy_level%matlu(:),matlu1(:),natom,-1)

     ! Print F(iw_n)=-(G_0^-1-iw_n+levels)  for last frequency.
     ! --------------------------------------------------------
     if (ifreq == nwlo .or. ifreq == paw_dmft%dmftqmc_l) then
       write(message,'(a,2x,a,i4,a,f13.5,a)') ch10, &
         & " == Print F(iw_n)=-(G_0^-1-iw_n+levels) for freq nb",ifreq," (=",paw_dmft%omega_lo(ifreq),")"
       call wrtout(std_out,message,'COLL')
       call print_matlu(matlu1(:),natom,1,opt_exp=1,compl=1)
     end if
     if (optdb == 1) then
       call printplot_matlu(matlu1(:),natom,paw_dmft%omega_lo(ifreq),"Hybridization",10000,imre=1)
     end if

     ! Put F in weiss_for_rot -> CTQMC
     ! -------------------------------
     if (opt_rot == 2) then
       call rotate_matlu(weiss_for_rot%oper(ifreq)%matlu(:),eigvectmatlu(:),natom,1)
     end if
     ! The following line will produce directly the weiss function for the CTQMC code
     if (opt_fk == 1) then
       call copy_matlu(matlu1(:),weiss_for_rot%oper(ifreq)%matlu(:),natom)
     end if

     ! Multiply F by frequency
     ! ------------------------
     call copy_matlu(matlu1(:),matlu2(:),natom)
     call fac_matlu(matlu1(:),natom,cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp))
     if (ifreq == nwlo .or. ifreq == paw_dmft%dmftqmc_l) then
       write(message,'(a,2x,a,i4,a,f13.5,a)') ch10, &
         & " == Print numerical C_ij = F(iw_n)*iw_n for freq nb",ifreq," (=",paw_dmft%omega_lo(ifreq),")"
       call wrtout(std_out,message,'COLL')
       call print_matlu(matlu1(:),natom,1,opt_exp=1,compl=1)
     end if
     if (optdb == 1) then
       call printplot_matlu(matlu1(:),natom,paw_dmft%omega_lo(ifreq),"cij",72800,imre=1)
     end if
     ! call rotate_matlu(matlu1,eigvectmatlu,natom,3,1)

     if (ifreq == nwlo .or. ifreq == paw_dmft%dmftqmc_l) then
       write(message,'(a,2x,a,i4,a,f13.5,a)') ch10, &
          & " == Print numerical after back rotation C_ij = F(iw_n)*iw_n for freq nb",ifreq," (=",paw_dmft%omega_lo(ifreq),")"
       call wrtout(std_out,message,'COLL')
       call print_matlu(matlu1(:),natom,1,opt_exp=1,compl=1)
     end if
     if (optdb == 1) then
       call printplot_matlu(matlu1(:),natom,paw_dmft%omega_lo(ifreq),"cij_rotated",72900,imre=1)
     end if

     ! Built C_ij/iw_n
     ! ------------------------
     call copy_matlu(hybri_coeff(:),matlu1(:),natom)
     call fac_matlu(matlu1(:),natom,cone/cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp))
     if (optdb == 1) then
       call printplot_matlu(matlu1(:),natom,paw_dmft%omega_lo(ifreq),"cij_over_omega",72000)
     end if
     ! if(ifreq==paw_dmft%dmft_nwlo) then
     !   write(message,'(a,2x,a,f13.5)') ch10,  " == Print numerical C_ij/iw_n for frequency",paw_dmft%omega_lo(ifreq)
     !   call wrtout(std_out,message,'COLL')
     !   call print_matlu(matlu1,natom,1,opt_exp=1,compl=1)
     ! endif

     ! For test: put C_ij/i_wn into weiss_for_rot
     ! --------------------------------------------
     ! call copy_matlu(matlu1,weiss_for_rot%oper(ifreq)%matlu,natom,opt_non_diag=1)

     ! Compute Hybri - C_ij/iw_n
     ! ------------------------
     call add_matlu(matlu2(:),matlu1(:),matlu3(:),natom,-1)

     ! Print Hybri - C_ij/iw_n
     ! ------------------------
     if (optdb == 1) then
       call printplot_matlu(matlu3(:),natom,paw_dmft%omega_lo(ifreq),"hybri_minus_asymp",74000,imre=1)
     end if

     ! Multiply (F-C_ij/i_wn) by (iw_n)**2 to find D_ij such that (F-C_ij/i_wn) -> D_ij/(iw_n)^2 only for last frequency.
     ! ------------------------------------------------------------------------------------------------------------------
     call copy_matlu(matlu3(:),matlu2(:),natom)
     call fac_matlu(matlu2(:),natom,cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)**2)
     if (optdb == 1) then
       call printplot_matlu(matlu2(:),natom,paw_dmft%omega_lo(ifreq),"fminuscijtimesw2",75000,imre=1)
     end if
     if (ifreq == nwlo .or. ifreq == paw_dmft%dmftqmc_l) then
       call copy_matlu(matlu2(:),matlu4(:),natom)
       write(message,'(a,2x,a,i4,a,f13.5,a)') ch10, &
          & " == Print numerical (F(iw_n)-C_ij/iw_n)%iw_n^2 for freq nb",ifreq," (=",paw_dmft%omega_lo(ifreq),")"
       call wrtout(std_out,message,'COLL')
       call print_matlu(matlu4(:),natom,1)
     end if

     ! Built C_ij/iw_n+D_ij/(iw_n)^2
     ! ------------------------
     call copy_matlu(matlu4(:),matlu3(:),natom,opt_re=1)
     call fac_matlu(matlu3(:),natom,cone/cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)**2)
     call add_matlu(matlu1(:),matlu3(:),matlu2(:),natom,1)
     if (optdb == 1) then
       call printplot_matlu(matlu2(:),natom,paw_dmft%omega_lo(ifreq),"cij_w_plus_dij_w2",72700,imre=1)
     end if
     ! For test: put C_ij/i_wn +D_ij/(iw_n)^2 into weiss_for_rot
     ! --------------------------------------------
     ! call copy_matlu(matlu2,weiss_for_rot%oper(ifreq)%matlu,natom,opt_non_diag=1)


   end do ! ifreq

   call destroy_matlu(matlu1(:),natom)
   call destroy_matlu(matlu2(:),natom)
   call destroy_matlu(matlu3(:),natom)
   call destroy_matlu(matlu4(:),natom)
   ABI_FREE(matlu1)
   ABI_FREE(matlu2)
   ABI_FREE(matlu3)
   ABI_FREE(matlu4)
 end if ! if opt_nondiag=1

 ! =======================
 !
 ! Rotation of Magnetic moment for CT-QMC
 !
 ! =======================
 if(nspinor .eq. 2 .and. paw_dmft%dmftctqmc_localprop .gt. 1) then
   write(message,'(a,2x,2a)') ch10, " == Making rotation for magnetic moments", ch10
   call wrtout(std_out,message,'COLL')

   !create a rotation matrix for diagonal Hamiltonian
   if(opt_diag == 0) then
     write(message,'(a,2x,2a)') ch10, " --> Hamiltonian is already diagonal in Slm", ch10
     call wrtout(std_out,message,'COLL')
     do iatom = 1,paw_dmft%natom
       if(paw_dmft%lpawu(iatom) /= -1) then
         do iflavor1=1,tndim
           do iflavor2=1,tndim
             if(iflavor1==iflavor2) then
               eigvectmatlu(iatom)%mat(iflavor1,iflavor2,1)=cone
             else
               eigvectmatlu(iatom)%mat(iflavor1,iflavor2,1)=czero
             end if
           end do
         end do
       end if
     end do
   end if !end opt_diag=0

   ! == orbital angular momentum
   call init_matlu(natom=natom,nspinor=paw_dmft%nspinor,nsppol=paw_dmft%nsppol,lpawu_natom=paw_dmft%lpawu,matlu=matlumag_orb)
   call zero_matlu(matlumag_orb,natom=natom)
   call chi_matlu(matlumag_orb,natom=natom,option=1,optprt=0)
   call rotate_matlu(matlumag_orb,eigvectmatlu,natom=natom,inverse=1)
   !call print_matlu(matlumag_orb,iatom,prtopt=1)
   call gather_matlu(matlumag_orb,magmom_orb,natom=natom,option=1,prtopt=0)
   call destroy_matlu(matlumag_orb,natom=natom)

   ! == spin angular momentum
   call init_matlu(natom=natom,nspinor=paw_dmft%nspinor,nsppol=paw_dmft%nsppol,lpawu_natom=paw_dmft%lpawu,matlu=matlumag_spin)
   call zero_matlu(matlumag_spin,natom=natom)
   call chi_matlu(matlumag_spin,natom=natom,option=2,optprt=0)
   call rotate_matlu(matlumag_spin,eigvectmatlu,natom=natom,inverse=1)
   !call print_matlu(matlumag_spin,natom,prtopt=1)
   call gather_matlu(matlumag_spin,magmom_spin,natom=natom,option=1,prtopt=0)
   call destroy_matlu(matlumag_spin,natom=natom)

   ! == total angular momentum
   call init_matlu(natom=natom,nspinor=paw_dmft%nspinor,nsppol=paw_dmft%nsppol,lpawu_natom=paw_dmft%lpawu,matlu=matlumag_tot)
   call zero_matlu(matlumag_tot,natom=natom)
   call chi_matlu(matlumag_tot,natom=natom,option=3,optprt=0)
   call rotate_matlu(matlumag_tot,eigvectmatlu,natom=natom,inverse=1)
   !call print_matlu(matlumag_tot,natom=1,prtopt=1)
   call gather_matlu(matlumag_tot,magmom_tot,natom=natom,option=1,prtopt=0)
   call destroy_matlu(matlumag_tot,natom=natom)

   write(message,'(a,2x,2a)') ch10, " ==> Rotation done", ch10
   call wrtout(std_out,message,'COLL')

 end if ! dmftctqmc_localprop
 !======================

 ! =========================================================================================
 ! Start big loop over atoms to compute hybridization and do the CTQMC
 ! =========================================================================================

 do iatom=1,natom

   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle

   green%ecorr_qmc(iatom) = zero
   itypat = paw_dmft%typat(iatom)
   tndim = 2*lpawu + 1
   ! if(lpawu/=-1) then

   nflavor = 2 * tndim

   if (testcode >= 1) then
     nflavor = 2
     if (testcode == 2) then
       ispa = 1
       ispb = 2
       if (nspinor == 1) ispb = 1
       ima = 1
       imb = 1
       if (tndim > 4) then
         ima = 5 ! row
         imb = 4 ! column
       end if
     end if ! testcode=2
   end if ! testcode >=1

   ABI_MALLOC(fw1,(nwlo,nflavor))
   ABI_MALLOC(fw1_nd,(nwlo,nflavor,nflavor))
   ABI_MALLOC(levels_ctqmc,(nflavor))
   ABI_MALLOC(levels_ctqmc_nd,(nflavor,nflavor))
   levels_ctqmc_nd(:,:) = czero
   ABI_MALLOC(hybri_limit,(nflavor,nflavor))
   hybri_limit(:,:) = czero
   fw1_nd(:,:,:) = czero
   fw1(:,:) = czero

   ! =================================================================
   ! Put hybridization in new arrays for CTQMC
   ! =================================================================
   if (testcode == 0) then
     iflavor1 = 0
     iflavor2 = 0

     do isppol=1,nsppol
       do ispinor1=1,nspinor
         do ispinor2=1,nspinor
           do im1=1,tndim
             do im2=1,tndim

                 ! first diagonal terms whatever opt_nondiag
               iflavor1 = im1 + tndim*(ispinor1-1) + tndim*(isppol-1)
               iflavor2 = im2 + tndim*(ispinor2-1) + tndim*(isppol-1)

               if (iflavor1 == iflavor2 ) then

                   ! Put weiss_for_rot in fw1
                 do ifreq=1,nwlo
                   if (opt_fk == 1 .or. opt_fk == 0) fw1(ifreq,iflavor1) = &
                     & weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(im1+(ispinor1-1)*tndim,im1+(ispinor1-1)*tndim,isppol)
                 end do  ! ifreq
                 fw1_nd(:,iflavor1,iflavor1) = fw1(:,iflavor1)

                 levels_ctqmc(iflavor1) = &
                    & dble(energy_level%matlu(iatom)%mat(im1+(ispinor1-1)*tndim,im1+(ispinor1-1)*tndim,isppol))
                 hybri_limit(iflavor1,iflavor1) = hybri_coeff(iatom)%mat(im1+(ispinor1-1)*tndim,im1+(ispinor1-1)*tndim,isppol)


                   ! case nsppol=nspinor=1
                 if (nsppol == 1 .and. nspinor == 1) then
                   fw1(:,iflavor1+tndim) = fw1(:,iflavor1)
                   fw1_nd(:,iflavor1+tndim,iflavor1+tndim) = fw1(:,iflavor1)
                   levels_ctqmc(iflavor1+tndim) = levels_ctqmc(iflavor1)
                   hybri_limit(iflavor1+tndim,iflavor1+tndim) = hybri_limit(iflavor1,iflavor1)
                 end if

               ! off diagonal terms
               else

                 ! Put weiss_for_rot in fw1_nd
                 do ifreq=1,nwlo
                   if (opt_fk == 1 .or. opt_fk == 0) fw1_nd(ifreq,iflavor1,iflavor2) = &
                     & weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(im1+(ispinor1-1)*tndim,im2+(ispinor2-1)*tndim,isppol)
                 end do ! ifreq
                 hybri_limit(iflavor1,iflavor2) = hybri_coeff(iatom)%mat(im1+(ispinor1-1)*tndim,im2+(ispinor2-1)*tndim,isppol)

                 ! case nsppol=nspinor=1
                 if (nsppol == 1 .and. nspinor == 1) then
                   fw1_nd(:,iflavor1+tndim,iflavor2+tndim) = fw1_nd(:,iflavor1,iflavor2)
                   hybri_limit(iflavor1+tndim,iflavor2+tndim) = hybri_limit(iflavor1,iflavor2)
                 end if

               end if ! iflavor1=iflavor2

! <  / HACK >
             end do ! im2
           end do ! im1
         end do  ! ispinor2
       end do  ! ispinor1
     end do  ! isppol
! < HACK >
     ! JB. On 1000 cpus this can not work since all CPU try to open/write the files
     ! Action : Don't print it or check only one cpu does it.

     if (pawprtvol >= 10000000) then
       write(message,'(a,2x,a)') ch10,  " == Hybri for all flavors for CTQMC "
       call wrtout(std_out,message,'COLL')
       do iflavor1=1,nflavor
         write(message,'(4x,14(2e14.5,2x))') (hybri_limit(iflavor1,iflavor2),iflavor2=1,nflavor)
         call wrtout(std_out,message,'COLL')
       end do ! iflavor1

       if (open_file('Hybri_cijoveromega',message,newunit=unt,status='unknown',form='formatted') /= 0) &
         & ABI_ERROR(message)
       if (open_file('Hybri',message,newunit=unt2,status='unknown',form='formatted') /= 0) ABI_ERROR(message)
       do ifreq=1,nwlo
         !  weiss_for_rot is G_0^-1-iw_n=-(F-levels)
         if (optdb == 1) then
           call printplot_matlu(weiss_for_rot%oper(ifreq)%matlu(:),natom,paw_dmft%omega_lo(ifreq),"weissbefore112",30000)
         end if
       end do
       do iflavor1=1,nflavor
         do iflavor2=1,nflavor
           do ifreq=1,nwlo
             omega = pi * paw_dmft%temp * (two*float(ifreq)-1)
             ! fw1_nd is -G_0^+1-iw_n=(F-levels)
             write(unt,'(300e16.5)') paw_dmft%omega_lo(ifreq), &
               & fw1_nd(ifreq,iflavor1,iflavor2)-hybri_limit(iflavor1,iflavor2)/cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)
             write(unt2,'(300e16.5)') paw_dmft%omega_lo(ifreq),fw1_nd(ifreq,iflavor1,iflavor2)
           end do ! ifreq
           write(unt,*)
           write(unt2,*)
         end do ! iflavor2
       end do ! iflavor1
       close(unt)
       close(unt2)
     end if ! pawprtvol>=10000000
   end if ! testcode
! </ HACK >

     ! ====================================================================================
     !  TEST
     !  For testing purpose, built ultra simple hybridization (constant in
     !  imaginary time or very simple) or extract some part of the calculated hybridization
     ! ====================================================================================
   if (testcode >= 1) then
       !dmft_nwlo = paw_dmft%dmft_nwlo
     paw_dmft%dmft_nwlo = paw_dmft%dmftqmc_l
     ABI_MALLOC(gw1_nd,(paw_dmft%dmft_nwlo,nflavor,nflavor))
     gw1_nd(:,:,:) = czero

       !  Call testcode_ctqmc: built simple hybridization
       !--------------------------------------------------
     if (testcode == 1) then
       call testcode_ctqmc(paw_dmft%dmftqmc_l,fw1_nd(:,:,:),fw1(:,:),gtmp_nd(:,:,:),gw_tmp_nd(:,:,:),&
         & levels_ctqmc(:),hybri_limit(:,:),nflavor,1,paw_dmft%temp,testrot,testcode,umod(:,:))
       !  Select 2x2 hybridization matrix from the current larger matrix
       !  ima and imb are defined above.
       !----------------------------------------------------------------
     else if (testcode == 2) then
         !close(unt)
         !close(unt2)
       call testcode_ctqmc_b(energy_level,hybri_coeff,weiss_for_rot,paw_dmft%dmftqmc_l,fw1_nd(:,:,:),&
         & levels_ctqmc(:),levels_ctqmc_nd(:,:),hybri_limit(:,:),paw_dmft%temp,umod(:,:),opt_diag,opt_fk)
     end if

       ! Calculation of Inverse Green's function from hybridization
       !-------------------------------------------------------------
     do if1=1,2
       do if2=1,2
         do ifreq=1,paw_dmft%dmftqmc_l
           omega = pi * paw_dmft%temp * (two*dble(ifreq)-1)
           if (if1 == if2) then
             gw1_nd(ifreq,if1,if2) = (cmplx(zero,omega,kind=dp)-fw1_nd(ifreq,if1,if2))
           else
             gw1_nd(ifreq,if1,if2) = -fw1_nd(ifreq,if1,if2)
           end if
         end do ! ifreq
       end do ! if2
     end do ! if1
       ! Calculation of Green's function (call inverse)
       !-------------------------------------------------------------
     do ifreq=1,paw_dmft%dmftqmc_l
       call xginv(gw1_nd(ifreq,:,:),2)
     end do
     write(std_out,*) " testctqmc high frequency limit of hybridization",fw1_nd(paw_dmft%dmftqmc_l,:,:)

       ! Integrate Green's function
       !-------------------------------------------------------------
     do if1=1,2
       do if2=1,2
         call int_fct(gw1_nd(:,if1,if2),(if1==if2),2,paw_dmft,integral(if1,if2))  ! test_1
       end do
     end do
       ! Write Occupations
     write(std_out,*) "Occupation of model in matrix form"
     do if1=1,2
       write(std_out,'(2(2f13.5,3x))') ((integral(if1,if2)+conjg(integral(if2,if1)))/two,if2=1,2)
     end do
     write(std_out,*) "Limit of hybridization "
     do if1=1,2
       write(std_out,'(2(2f13.5,3x))') (hybri_limit(if1,if2),if2=1,2)
     end do

       ! If opt_fk=0, give Green's function to CTQMC code instead of
       ! hybridization
       !-------------------------------------------------------------
     if (opt_fk == 0) fw1_nd(:,:,:) = gw1_nd(:,:,:)

     ABI_FREE(gw1_nd)
     paw_dmft%dmft_nwlo = nwlo

     ! and testcode>1
   end if ! testcode>=1


   call flush_unit(std_out)
   ! =================================================================

   ! ___________________________________________________________________________________
   !
   !  THIRD PART : CALL CTQMC
   ! ___________________________________________________________________________________

   ! ==================================================================
   !    Main calls to CTQMC code in ABINIT (INITIALIZATION and OPTIONS)
   ! ==================================================================
   if (paw_dmft%dmft_solv == 5 .or. paw_dmft%dmft_solv == 8) then
     write(message,'(a,2x,a)') ch10," == Initializing CTQMC"
     call wrtout(std_out,message,'COLL')

     !    Initialisation
     ! =================================================================
     if (paw_dmft%dmft_solv == 5) then
       nomega = paw_dmft%dmftqmc_l
       call CtqmcInterface_init(hybrid,paw_dmft%dmftqmc_seed,paw_dmft%dmftqmc_n, &
         & paw_dmft%dmftqmc_therm,paw_dmft%dmftctqmc_meas,nflavor,paw_dmft%dmftqmc_l,&
         & one/paw_dmft%temp,zero,std_out,paw_dmft%spacecomm,paw_dmft%nspinor)
       !    options
       ! =================================================================
       call CtqmcInterface_setOpts(hybrid, &
          & opt_Fk       = opt_fk, &
          & opt_order    = paw_dmft%dmftctqmc_order, &
          & opt_histo    = paw_dmft%dmftctqmc_localprop, &
          & opt_movie    = paw_dmft%dmftctqmc_mov, &
          & opt_analysis = paw_dmft%dmftctqmc_correl, &
          & opt_check    = paw_dmft%dmftctqmc_check, &
          & opt_noise    = paw_dmft%dmftctqmc_grnns, &
          & opt_spectra  = paw_dmft%dmftctqmc_mrka, &
          & opt_gmove    = paw_dmft%dmftctqmc_gmove)
     end if

     if (paw_dmft%dmft_solv == 8) then
       nomega = paw_dmft%dmftqmc_l
       call CtqmcoffdiagInterface_init(hybridoffdiag,paw_dmft%dmftqmc_seed,&
         & paw_dmft%dmftqmc_n,paw_dmft%dmftqmc_therm,paw_dmft%dmftctqmc_meas,&
         & nflavor,paw_dmft%dmftqmc_l,one/paw_dmft%temp,zero,std_out,&
         & paw_dmft%spacecomm,opt_nondiag,paw_dmft%nspinor)
       !    options
       ! =================================================================
       call CtqmcoffdiagInterface_setOpts(hybridoffdiag,opt_Fk=opt_fk, &
           & opt_order    = paw_dmft%dmftctqmc_order, &
           & opt_histo    = paw_dmft%dmftctqmc_localprop, &
           & opt_movie    = paw_dmft%dmftctqmc_mov, &
           & opt_analysis = paw_dmft%dmftctqmc_correl, &
           & opt_check    = paw_dmft%dmftctqmc_check, &
           & opt_noise    = paw_dmft%dmftctqmc_grnns, &
           & opt_spectra  = paw_dmft%dmftctqmc_mrka, &
           & opt_gmove    = paw_dmft%dmftctqmc_gmove)
     end if

     write(message,'(a,2x,2a)') ch10, " == Initialization CTQMC done", ch10
     call wrtout(std_out,message,'COLL')

   end if ! dmft_solv=5 or dmft_solv=8

   if (paw_dmft%dmft_solv == 9) then
     ABI_MALLOC(gw_tmp_nd,(paw_dmft%dmft_nwli,nflavor,nflavor))
     ! because size allocation problem with TRIQS paw_dmft%dmft_nwlo must be >= paw_dmft%dmft_nwli
     open(unit=505,file=trim(paw_dmft%filapp)//"_Legendre_coefficients.dat",status='unknown',form='formatted')
   else
     if (paw_dmft%dmft_solv == 5) then
       ABI_MALLOC(gw_tmp,(paw_dmft%dmft_nwlo,nflavor+1))
     end if
     ABI_MALLOC(gw_tmp_nd,(paw_dmft%dmft_nwlo,nflavor,nflavor+1))
       !use  gw_tmp to put freq
     do ifreq=1,paw_dmft%dmft_nwlo
       if (paw_dmft%dmft_solv == 5) gw_tmp(ifreq,nflavor+1) = cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)
       gw_tmp_nd(ifreq,nflavor,nflavor+1) = cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)
     end do
   end if ! dmft_solv=9

   ABI_MALLOC(gtmp,(paw_dmft%dmftqmc_l,nflavor))
     ! THIS IS A BACKUP PLAN. USING paw_dmft%hybrid makes a segfault on TIKAL
     ! PSC with MPI only (and max2_open64). paw_dmf%hybrid is corrupted
     ! somewhere but I could not find the place in all DMFT routines
   ABI_MALLOC(gtmp_nd,(paw_dmft%dmftqmc_l,nflavor,nflavor))
   call flush_unit(std_out)

     ! =================================================================
     !    BEGIN CALL TO CTQMC SOLVERS
     ! =================================================================

   if (testcode == 0) then

       ! =================================================================
       !    CTQMC run Abinit
       ! =================================================================
     if (paw_dmft%dmft_solv == 5) then

       ABI_MALLOC(docc,(nflavor,nflavor))
       docc(:,:) = zero
       call CtqmcInterface_run(hybrid,fw1(1:paw_dmft%dmftqmc_l,:),Gtau=gtmp(:,:),Gw=gw_tmp(:,:),D=docc(:,:),E=green%ecorr_qmc(iatom),&
         ! & matU=hu(itypat)%udens,opt_levels=levels_ctqmc)
         & matU=dble(udens_atoms(iatom)%mat(:,:,1)),opt_levels=levels_ctqmc(:),Magmom_orb=REAL(magmom_orb(iatom)%value),&
         & Magmom_spin=REAL(magmom_spin(iatom)%value),Magmom_tot=REAL(magmom_tot(iatom)%value),Iatom=iatom,fname=paw_dmft%filapp)
       if (paw_dmft%dmft_entropy > 0) then
         call data4entropyDMFT_setDocc(paw_dmft%forentropyDMFT,iatom,docc(:,:))
       end if
       ABI_FREE(docc)
       !DO iflavor = 1, nflavor
       !  hybrid%Hybrid%Greens(iflavor)%oper(1:this%samples) = gtmp(1:this%samples,iflavor)
       !  CALL GreenHyb_forFourier(this%Greens(iflavor), Gomega=Gw(:,iflavor), omega=Gw(:,this%flavors+1))
       !END DO

       ! =================================================================
       !    CTQMC run Abinit off diagonal terms in hybridization
       ! =================================================================
     else if (paw_dmft%dmft_solv == 8) then
       ! =================================================================

       ABI_MALLOC(docc,(nflavor,nflavor))
       docc(:,:) = zero

       call CtqmcoffdiagInterface_run(hybridoffdiag,fw1_nd(1:paw_dmft%dmftqmc_l,:,:),Gtau=gtmp_nd(:,:,:),&
          & Gw=gw_tmp_nd(:,:,:),D=doccsum,E=green%ecorr_qmc(iatom),Noise=noise,matU=dble(udens_atoms(iatom)%mat(:,:,1)),&
          & Docc=docc(:,:),opt_levels=levels_ctqmc(:),hybri_limit=hybri_limit(:,:),Magmom_orb=REAL(magmom_orb(iatom)%value),&
          & Magmom_spin=REAL(magmom_spin(iatom)%value),Magmom_tot=REAL(magmom_tot(iatom)%value),Iatom=iatom,fname=paw_dmft%filapp)

       ! For entropy (alternative formulation)
       if (paw_dmft%ientropy == 1) then
         EE = zero
         do if1=1,nflavor
           do if2=if1+1,nflavor
             EE = EE + docc(if1,if2)*dble(udens_atoms_for_s(iatom)%mat(if1,if2,1))
             ! write(std_out,*) udens_atoms_for_s(iatom)%value(if1,if2),docc(if1,if2)
           end do
         end do
         ! Here in udens U=1, J=J/U, so we need to multiply bu U/Ha_eV
         write(message,'(a,3(f14.10,3x))') "For entropy calculation E_corr_qmc, u_for_s, j_for,s", &
             & paw_dmft%u_for_s*EE/Ha_eV,paw_dmft%u_for_s,paw_dmft%j_for_s
         call wrtout(std_out,message,'COLL')
         EE = zero
         do if1=1,nflavor
           do if2=if1+1,nflavor
             EE = EE + dble(docc(if1,if2)*udens_atoms(iatom)%mat(if1,if2,1))
             ! write(std_out,*) udens_atoms(iatom)%value(if1,if2),docc(if1,if2)
           end do
         end do
         ! Here in udens U=U, J=J, so we obtain directly the results
         write(message,'(a,3(f14.10,3x))') "Reference   calculation E_corr_qmc, upawu  , jpawu  ", &
             & EE,hu(itypat)%upawu*Ha_eV,hu(itypat)%jpawu*Ha_eV
         call wrtout(std_out,message,'COLL')
       end if
       ABI_FREE(docc)
       ! TODO: Handle de luj0 case for entropy

       ! =================================================================
       !    CTQMC run TRIQS
       ! =================================================================
     else if (paw_dmft%dmft_solv == 9) then
       ! =================================================================

       call ctqmc_calltriqs(paw_dmft,cryst_struc,hu(:),levels_ctqmc,gtmp_nd,gw_tmp_nd,fw1_nd,leg_measure,iatom)

     end if

   ! =================================================================
   !    CTQMC run for tests
   ! =================================================================
   else if (testcode >= 1) then
     call CtqmcInterface_run(hybrid,fw1(1:nomega,:),Gtau=gtmp(:,:),Gw=gw_tmp(:,:),E=green%ecorr_qmc(iatom),&
        & matU=umod(:,:),opt_levels=levels_ctqmc(:),Iatom=iatom,fname=paw_dmft%filapp)

     ! for non diagonal code
     !       call CtqmcInterface_run(hybrid,fw1_nd(1:nomega,:,:),Gtau=gtmp_nd,&
     !&       Gw=gw_tmp_nd,D=Doccsum,E=green%ecorr_qmc(iatom),&
     !&       Noise=Noise,matU=umod,opt_levels=levels_ctqmc,hybri_limit=hybri_limit)

     !  If test of the code is activated, and testrot =1 rotate back green's function   and stop the code.
     ! --------------------------------------------------------------------------------------------------
     if (testcode == 1) then

       call testcode_ctqmc(paw_dmft%dmftqmc_l,fw1_nd(:,:,:),fw1(:,:),gtmp_nd(:,:,:),gw_tmp_nd(:,:,:), &
            & levels_ctqmc(:),hybri_limit(:,:),nflavor,2,paw_dmft%temp,testrot,testcode,umod(:,:))

       write(message,'(2a)') ch10,' testcode end of test calculation'
       ABI_ERROR(message)
     end if

     if (testcode == 2) then
       write(message,'(2a)') ch10,' testcode 2 end of test calculation'
       ABI_ERROR(message)
     end if

   end if
   ! =================================================================
   !    END CALL TO CTQMC SOLVERS
   ! =================================================================


   ! Print green function is files directly from CTQMC
   ! --------------------------------------------------
   call ctqmcoutput_printgreen(paw_dmft,gtmp_nd,gw_tmp_nd,gtmp,gw_tmp,iatom)


   ! If the CTQMC code in ABINIT was used, then destroy it and deallocate arrays
   ! ----------------------------------------------------------------------------
   ! if(paw_dmft%dmft_solv<6.and.paw_dmft%dmft_solv>7) then
   ! Nothing just hybrid var problem
   ! else
   write(message,'(a,2x,a)') ch10," == Destroy CTQMC"
   call wrtout(std_out,message,'COLL')
   if (paw_dmft%dmft_solv == 5) then
     call CtqmcInterface_finalize(hybrid)
   end if
   if (paw_dmft%dmft_solv == 8) then
     call CtqmcoffdiagInterface_finalize(hybridoffdiag)
   end if
   write(message,'(a,2x,a)') ch10," == Destroy CTQMC done"
   call wrtout(std_out,message,'COLL')
   ABI_FREE(hybri_limit)
   ABI_FREE(levels_ctqmc_nd)
   ABI_FREE(levels_ctqmc)
   ABI_FREE(fw1)
   ABI_FREE(fw1_nd)

   ! ____________________________________________________________
   !
   !  FOURTH PART : USE OUTPUT OF CTQMC AND THEN DO BACK ROTATION
   ! ____________________________________________________________
   !

   ! Put green's function values from CTQMC into green structure
   !------------------------------------------------------------
   call ctqmcoutput_to_green(green,paw_dmft,gtmp_nd,gw_tmp_nd,gtmp,gw_tmp,iatom,leg_measure,opt_nondiag)

   ! Deallocate arrays for CTQMC
   !----------------------------
   if (paw_dmft%dmft_solv < 6) then
     ABI_FREE(gw_tmp)
   end if
   ABI_FREE(gw_tmp_nd)
   ABI_FREE(gtmp)
   ABI_FREE(gtmp_nd)

   ! Do Fourier transform if it was not done (ie if TRIQS is used without legendre measurement)
   !-------------------------------------------------------------------------------------------
   ! if(opt_nondiag==1) then  ! (As leg_measure is activated by defautl, this fourier is never done).
   !   if(paw_dmft%dmft_solv>=6.and..not.leg_measure.and.paw_dmft%dmft_solv<=7) then
   !     write(message,'(2a,i3,13x,a)') ch10,'   ===  Direct Fourier Transform t->w of Weiss Field'
   !     call wrtout(std_out,message,'COLL')
   !     call fourier_green(cryst_struc,green,paw_dmft,&
   !     & pawang,opt_ksloc=2,opt_tw=1)
   !     end if
   !   endif

  ! end if

 end do ! iatom
 ! ==================================================================
 !  End big loop over atoms to compute hybridization and do the CTQMC
 ! ==================================================================

 if (paw_dmft%dmft_prgn == 1) then
   call print_green('QMC_diag_notsym',green,1,paw_dmft,opt_wt=2)
   call print_green('QMC_diag_notsym',green,1,paw_dmft,opt_wt=1)
 end if
 ! write(message,'(i3,4x,2e21.14)') 6,weiss_for_rot%oper(1)%matlu(1)%mat(1,1,1,1,1)
 ! call wrtout(std_out,message,'COLL')  ! debug
 ! =================================================================
 ! Inverse Weiss, then
 ! Copy Weiss_for_rot into weiss and rotate back weiss to the original basis
 ! =================================================================

 ! ABI_MALLOC(shift,(natom))
 ! do ifreq=1,paw_dmft%dmft_nwlo
 !  ! First weiss_for_rot contains -G_0^-1+iw_n
 !  ! -------------------------------------------
 !  ! Compute G_0^-1-iw_n
 !  ! --------------------
 !       write(6,*) "1"
 !  if(opt_fk==1) call fac_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,-cone)
 !
 !
 !       write(6,*) "2"
 !  ! Compute G_0^-1
 !  ! --------------------
 !  shift(:)=cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)
 !  if(opt_fk==1) call shift_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,shift,signe=1)
 !
 !       write(6,*) "3"
 !  ! Compute G_0
 !  ! --------------------
 !   call inverse_oper(weiss_for_rot%oper(ifreq),option=1,prtopt=1)
 !   ! No need to copy if weiss_for_rot is a pointer to weiss ...
 !!   if(useylm==1) call slm2ylm_matlu(weiss%oper(ifreq)%matlu,natom,2,0)
 !!   if(opt_diag/=0) call rotate_matlu(weiss%oper(ifreq)%matlu,eigvectmatlu,natom,3,0)
 !
 !  ! Compute G_0 in the original basis
 !  ! --------------------
 !   call rotate_matlu(weiss_for_rot%oper(ifreq)%matlu,eigvectmatlu,natom,3,0)
 ! end do
 ! ABI_FREE(shift)

 ! =================================================================
 ! Here compute Self energy from Dyson and print it
 ! Warning : Weiss_for_rot is inversed inside dyson
 ! =================================================================
 ! call initialize_self(self,paw_dmft)
 ! call dyson(green,paw_dmft,self,weiss_for_rot,opt_weissself=2)
 ! call rw_self(self,mpi_enreg,paw_dmft,prtopt=2,opt_rw=2,opt_char="diag")
 ! call destroy_self(self)
  !write(message,'(i3,4x,2e21.14)') 7,weiss%oper(1)%matlu(1)%mat(1,1,1,1,1)
  !call wrtout(std_out,message,'COLL')  ! debug

! =================================================================
! Rotate back green function to original basis (non-diagonal)
!  (and Weiss for further use: might be useful if an back Fourier
!     transformation is done).
! =================================================================
 if (pawprtvol >= 3) then
   write(message,'(a,2x,a)') ch10, &  ! debug
      & " == Print Green's function for tau=0+ in the CTQMC basis"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(green%oper_tau(1)%matlu(:),natom,1)  ! debug
   write(message,'(a,2x,a)') ch10,&  ! debug
      & " == Print Green's function for smallest freq in the CTQMC basis"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(green%oper(1)%matlu(:),natom,1)  ! debug
 end if ! pawprtvol>=3

 !  === Compute rotated Occupations in green%occup_tau
 call occup_green_tau(green)

 if (pawprtvol >= 3) then
 ! === Compute non rotated Occupations in green%occup_tau
   write(message,'(a,2x,a)') ch10," == Occupations from G(tau=0-) in the CTQMC basis"
   call wrtout(std_out,message,'COLL')
   call print_matlu(green%occup_tau%matlu(:),natom,1)
 end if ! pawprtvol>=3

 ! =================================================================
 !
 !  === Compute magnetic moments from CT-QMC occupations for
 !  the x,y and z axes when SOC is activated
 !
 ! =================================================================
 if (paw_dmft%nspinor .eq. 2) then
  ABI_MALLOC(matlumag,(natom))
  write(message,'(a,2x,a)') ch10,"== Magnetic moments from CT-QMC occupation matrix "
  call wrtout(std_out,message,'COLL')

  do iatom=1,cryst_struc%natom
    lpawu=paw_dmft%lpawu(iatom)
    if(lpawu .ne. -1) then
      write(message,'(a,3x,a,i4)') ch10,"-------> For Correlated Atom",iatom
      call wrtout(std_out,message,'COLL')

      ! == orbital angular momentum
      do icomp=1,3 !x,y,z components
        muorb=czero
        call init_matlu(natom=natom,nspinor=paw_dmft%nspinor,nsppol=paw_dmft%nsppol,lpawu_natom=paw_dmft%lpawu,matlu=matlumag)
        call copy_matlu(green%occup_tau%matlu,matlumag,natom)
        call rotate_matlu(matlumag,eigvectmatlu,natom=natom,inverse=0)
        call magmomforb_matlu(matlumag,muorb,natom=natom,option=icomp,optprt=0)
        write(message,'(a,2x,a,i4,a,f8.4)') ch10," Orbital angular momentum for axis ", icomp, " is ", REAL(muorb)
        call wrtout(std_out,message,'COLL')
        call destroy_matlu(matlumag,(natom))
      end do

      ! == spin angular momentum
      do icomp=1,3 !x,y,z components
        muspin=czero
        call init_matlu(natom=natom,nspinor=paw_dmft%nspinor,nsppol=paw_dmft%nsppol,lpawu_natom=paw_dmft%lpawu,matlu=matlumag)
        call copy_matlu(green%occup_tau%matlu,matlumag,natom=natom)
        call rotate_matlu(matlumag,eigvectmatlu,natom=natom,inverse=0)
        call magmomfspin_matlu(matlumag,muspin,natom=natom,option=icomp,optprt=0)
        write(message,'(a,2x,a,i4,a,f8.4)') ch10," Spin angular momentum for axis ", icomp, " is ", REAL(muspin)
        call wrtout(std_out,message,'COLL')
        call destroy_matlu(matlumag,(natom))
      end do

      ! == total angular momentum (L_u + 2*S_u)
      do icomp=1,3 !x,y,z components
        muzeem=czero
        call init_matlu(natom=natom,nspinor=paw_dmft%nspinor,nsppol=paw_dmft%nsppol,lpawu_natom=paw_dmft%lpawu,matlu=matlumag)
        call copy_matlu(green%occup_tau%matlu,matlumag,natom=natom)
        call rotate_matlu(matlumag,eigvectmatlu,natom=natom,inverse=0)
        call magmomfzeeman_matlu(matlumag,muzeem,natom=natom,option=icomp,optprt=0)
        write(message,'(a,2x,a,i4,a,f8.4)') ch10," Zeeman angular momentum for axis ", icomp, " is ", REAL(muzeem)
        call wrtout(std_out,message,'COLL')
        call destroy_matlu(matlumag,(natom))
      end do
    endif !lpawu
  end do !iatom
  ABI_FREE(matlumag)
 end if !nspinor
 ! =================================================================

 if (opt_diag /= 0) then
   write(message,'(a,2x,a)') ch10," == Rotate Green's function back to original basis "
   call wrtout(std_out,message,'COLL')
 end if
 ! write(message,'(i3,4x,2e21.14)') 8,weiss%oper(1)%matlu(1)%mat(1,1,1,1,1)
 ! call wrtout(std_out,message,'COLL')  ! debug

 ! Rotate oper_tau into Ylm basis and then Slm basis
 ! -------------------------------------------------------------
 ! do itau=1,paw_dmft%dmftqmc_l
 !   if (opt_diag /= 0) call rotate_matlu(green%oper_tau(itau)%matlu(:),eigvectmatlu(:),natom,3,0)
 !   if (useylm == 1) call slm2ylm_matlu(green%oper_tau(itau)%matlu(:),natom,2,0)
 ! end do
 do itau=1,paw_dmft%dmftqmc_l
   if (opt_diag /= 0) then
     call rotate_matlu(green%oper_tau(itau)%matlu(:),eigvectmatlu(:),natom,0)
   end if
   if (useylm == 1) then
     call slm2ylm_matlu(green%oper_tau(itau)%matlu(:),natom,paw_dmft,2,0)
   end if
 end do ! itau

 ! Rotate occup_tau into Ylm basis and then Slm basis

 ! Rotate occup_tau into Ylm basis and then Slm basis
 !-------------------------------------------------------------
 if (opt_diag /= 0) then
   call rotate_matlu(green%occup_tau%matlu(:),eigvectmatlu(:),natom,0)
 end if
 if (useylm == 1) then
   write(message,'(a,2x,a)') ch10," == Occupations from G(tau=0-) in the Ylm basis"
   call wrtout(std_out,message,'COLL')
   call print_matlu(green%occup_tau%matlu(:),natom,1)
   call slm2ylm_matlu(green%occup_tau%matlu(:),natom,paw_dmft,2,0)
 end if
 ! write(message,'(a,2x,a)') ch10," == Occupations from G(tau=0-) in the Slm basis"
 ! call wrtout(std_out,message,'COLL')
 ! call print_matlu(green%occup_tau%matlu(:),natom,1)

 ! Put Weiss off diagonal terms to zero because Green function will not have any offdiag terms
 ! ------------------------------------------------------------------------------
 !   (if opt_nondiag=0 ie dmft_solv=5)
 if (opt_nondiag == 0) then
   do ifreq=1,nwlo
     call zero_matlu(weiss%oper(ifreq)%matlu(:),natom,onlynondiag=1)
   end do ! ifreq
 end if ! opt_nondiag=0
 !    ( if opt_nondiag=0, then:
 !       As Green's function is diagonal, one suppress off diag  terms in Weiss, if any.
 !      (If off diag are non zero in the density matrix and thus in the Green's function,
 !       there is a warning in checkreal_matlu above).)

 ! Rotate Green's and Weiss functions into Ylm basis and then Slm basis
 !-------------------------------------------------------------

 do ifreq=1,nwlo
   if (opt_diag /= 0) then
     call rotate_matlu(green%oper(ifreq)%matlu(:),eigvectmatlu(:),natom,0)
     call rotate_matlu(weiss%oper(ifreq)%matlu(:),eigvectmatlu(:),natom,0)
   end if
   if (useylm == 1) then
     call slm2ylm_matlu(green%oper(ifreq)%matlu(:),natom,paw_dmft,2,0)
     call slm2ylm_matlu(weiss%oper(ifreq)%matlu(:),natom,paw_dmft,2,0)
   end if
 end do ! ifreq

 !write(message,'(i3,4x,2e21.14)') 10,weiss%oper(1)%matlu(1)%mat(1,1,1,1,1)
 !call wrtout(std_out,message,'COLL')  ! debug

 if (pawprtvol >= 3) then
!< HACK >
   write(message,'(a,2x,a)') ch10,&  ! debug
       & " == Print diagonalized weiss_for_rot function after rotation for small freq in the ctqmc basis"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(weiss_for_rot%oper(1)%matlu(:),natom,1)  ! debug

!</ HACK >
   write(message,'(a,2x,a)') ch10,&  ! debug
     & " == Print Weiss function for smallest freq in the Slm basis"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(weiss%oper(1)%matlu(:),natom,1)  ! debug

   do ifreq=1,nwlo
     call sym_matlu(weiss%oper(ifreq)%matlu(:),paw_dmft)
   end do
   write(message,'(a,2x,a)') ch10,&  ! debug
     & " == Print symmetrized Weiss function for smallest freq in the Slm basis"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(weiss%oper(1)%matlu(:),natom,1)  ! debug
   write(message,'(a,2x,a)') ch10, &                  ! debug
     & " == Print Green's function for tau=0+ in the Slm basis" ! debug
   call wrtout(std_out,message,'COLL')                  ! debug
   call print_matlu(green%oper_tau(1)%matlu(:),natom,1)  ! debug
   write(message,'(a,2x,a)') ch10,&                  ! debug
     & " == Print Green's function for smallest freq in the Slm basis" ! debug
   call wrtout(std_out,message,'COLL')                  ! debug
   call print_matlu(green%oper(1)%matlu(:),natom,1)  ! debug
 end if ! pawprtvol>=3

 ABI_MALLOC(matlu1,(natom))
 call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),matlu1(:))
 call copy_matlu(green%occup_tau%matlu(:),matlu1(:),natom)
 call sym_matlu(matlu1(:),paw_dmft)

 write(message,'(a,2x,a)') ch10," == Occupations from G(tau=0-) in the Slm basis"
 call wrtout(std_out,message,'COLL')
 call print_matlu(green%occup_tau%matlu(:),natom,1)

 write(message,'(a,2x,a)') ch10," == Symmetrized occupations"
 call wrtout(std_out,message,'COLL')
 call print_matlu(matlu1(:),natom,1)

 call diff_matlu("CTQMC occupations","Symmetrized CTQMC occupations",green%occup_tau%matlu(:),matlu1(:),natom,0,tol4,ierr=ierr)
 call destroy_matlu(matlu1(:),natom)
 ABI_FREE(matlu1)

 ! =================================================================
 ! Symmetrize green function G(tau) and G(ifreq) to recover symmetry
 ! artificially broken by QMC
 ! =================================================================
 write(message,'(a,2x,a)') ch10," == Symmetrize Green's function after CTQMC "
 call wrtout(std_out,message,'COLL')

 do itau=1,1 !paw_dmft%dmftqmc_l
   call sym_matlu(green%oper_tau(itau)%matlu(:),paw_dmft)
 end do ! itau
 do ifreq=1,paw_dmft%dmft_nwlo
   call sym_matlu(green%oper(ifreq)%matlu(:),paw_dmft)
 end do ! ifreq
 if (pawprtvol >= 3) then
   write(message,'(a,2x,a)') ch10, &  ! debug
      & " == Print Green's function for tau=0+ after symmetrization"  !  debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(green%oper_tau(1)%matlu(:),natom,1)  ! debug
   write(message,'(a,2x,a)') ch10, &  ! debug
      & " == Print Green's function for smallest freq after symmetrization"  !  debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(green%oper(1)%matlu(:),natom,1)  ! debug
 end if ! pawprtvol>=3
 if (paw_dmft%dmft_prgn == 1) then
   call print_green('QMC_sym',green,1,paw_dmft,opt_wt=2)
   call print_green('QMC_sym',green,1,paw_dmft,opt_wt=1)
 end if

 ! === Compute Occupations  (Symmetrized from oper_tau)
 call occup_green_tau(green)

 ! === Print occupations
 ! call printocc_green(green,6,paw_dmft,3)

 call destroy_oper(energy_level)
 call destroy_matlu(dmat_diag(:),natom)
 call destroy_matlu(eigvectmatlu(:),natom)
 call destroy_matlu(udens_atoms(:),natom)
 ABI_FREE(dmat_diag)
 do iatom=1,natom
   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   ABI_FREE(magmom_orb(iatom)%value)
   ABI_FREE(magmom_spin(iatom)%value)
   ABI_FREE(magmom_tot(iatom)%value)
 end do
 ABI_FREE(udens_atoms)
 ABI_FREE(eigvectmatlu)
 ABI_FREE(magmom_orb)
 ABI_FREE(magmom_spin)
 ABI_FREE(magmom_tot)
 ABI_FREE(matlumag_orb)
 ABI_FREE(matlumag_spin)
 ABI_FREE(matlumag_tot)
 call destroy_green(weiss_for_rot)
 ! call destroy_green(gw_loc)
 ! call destroy_green(greendft)

 ! destroy limit of hybridization
 call destroy_matlu(hybri_coeff(:),paw_dmft%natom)
 ABI_FREE(hybri_coeff)

 call destroy_vee(paw_dmft,vee_rotated(:))
 ABI_FREE(vee_rotated)

end subroutine qmc_prep_ctqmc
!!***

!!****f* m_forctqmc/testcode_ctqmc_b
!! NAME
!! testcode_ctqmc_b
!!
!! FUNCTION
!! Setup ultra simple hybridization to test CTQMC in simple situations.
!!
!! INPUTS
!! temp = temperature
!! dmftqmc_l = number of times slices
!! levels_ctqmc_nd=level matrix
!!
!! OUTPUT
!! fw1_nd=hybridization matrix
!! umod = value of U
!! hybri_limit= limit of F
!! weiss_for_rot= weiss function
!! hybri_coeff
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

subroutine testcode_ctqmc_b(energy_level,hybri_coeff,weiss_for_rot,dmftqmc_l,fw1_nd,levels_ctqmc,&
&   levels_ctqmc_nd,hybri_limit,temp,umod,opt_diag,opt_fk)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: dmftqmc_l,opt_diag,opt_fk
 real(dp), intent(in) :: temp
 real(dp), intent(out) :: umod(2,2)
 real(dp), intent(inout) :: levels_ctqmc(:)
 complex(dp), intent(out) :: fw1_nd(:,:,:)
 complex(dp),  intent(inout) :: levels_ctqmc_nd(:,:)
 complex(dp),  intent(inout) :: hybri_limit(:,:)
 type(oper_type)  :: energy_level
 type(matlu_type), allocatable  :: hybri_coeff(:)
 type(green_type)  :: weiss_for_rot

!Local variables ------------------------------
 integer :: ifreq,iatom,ima,imb,ispa,ispb,ndim
 real(dp) :: omega
 real(dp) :: facnd, facd
 character(len=30) :: tmpfil
! ************************************************************************
 facnd=0.8d0
 facd=1.0d0
 ndim=2*energy_level%matlu(iatom)%lpawu+1
 !write(6,*) "fac",facnd,facd
 levels_ctqmc_nd(2,2)   = energy_level%matlu(iatom)%mat(imb+(ispb-1)*ndim,imb+(ispb-1)*ndim,1)
 levels_ctqmc_nd(1,1)   = energy_level%matlu(iatom)%mat(ima+(ispa-1)*ndim,ima+(ispa-1)*ndim,1)
 levels_ctqmc(2)   = real(energy_level%matlu(iatom)%mat(imb+(ispb-1)*ndim,imb+(ispb-1)*ndim,1),kind=dp)
 levels_ctqmc(1)   = real(energy_level%matlu(iatom)%mat(ima+(ispa-1)*ndim,ima+(ispa-1)*ndim,1),kind=dp)
 if(opt_diag/=1) then
   levels_ctqmc_nd(1,2)   = energy_level%matlu(iatom)%mat(ima+(ispa-1)*ndim,imb+(ispb-1)*ndim,1)
   levels_ctqmc_nd(2,1)   = energy_level%matlu(iatom)%mat(imb+(ispb-1)*ndim,ima+(ispa-1)*ndim,1)
 end if
 hybri_limit(1,1)  = facd*hybri_coeff(iatom)%mat(ima+(ispa-1)*ndim,ima+(ispa-1)*ndim,1)
 hybri_limit(2,2)  = facd*hybri_coeff(iatom)%mat(imb+(ispb-1)*ndim,imb+(ispb-1)*ndim,1)
 hybri_limit(1,2)  = facnd*hybri_coeff(iatom)%mat(ima+(ispa-1)*ndim,imb+(ispb-1)*ndim,1)
 hybri_limit(2,1)  = facnd*hybri_coeff(iatom)%mat(imb+(ispb-1)*ndim,ima+(ispa-1)*ndim,1)
 !write(6,*) "hybri_limit",hybri_limit
 !write(6,*) "levels_ctqmc",levels_ctqmc
 umod=zero

 tmpfil = 'fw1_nd_re'
 !if (open_file(newunit=unt,message,file=trim(tmpfil),status='unknown',form='formatted')/=0) then
 !  ABI_ERROR(message)
 !end if
 tmpfil = 'fw1_nd_im'
 !if (open_file(newunit=unt2,message,file=trim(tmpfil),status='unknown',form='formatted')/=0) then
 !  ABI_ERROR(message)
 !end if
 write(std_out,*) "testcode==2",ispa,ispb,ima,imb
 write(std_out,*) "opt_fk==",opt_fk
 do ifreq=1,dmftqmc_l
   if (opt_fk==1) then
     fw1_nd(ifreq,1,1) = facd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(ima+(ispa-1)*ndim,ima+(ispa-1)*ndim,1)
     fw1_nd(ifreq,2,2) = facd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(imb+(ispb-1)*ndim,imb+(ispb-1)*ndim,1)
     !fw1_nd(ifreq,1,2) =  weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(ima,imb,1,ispa,ispb)
     !fw1_nd(ifreq,2,1) =  weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(imb,ima,1,ispb,ispa)
     fw1_nd(ifreq,1,2) = facnd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(ima+(ispa-1)*ndim,imb+(ispb-1)*ndim,1)
     fw1_nd(ifreq,2,1) = facnd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(imb+(ispb-1)*ndim,ima+(ispa-1)*ndim,1)
     omega=pi*temp*(two*float(ifreq)-1)
   else if (opt_fk==0) then
     fw1_nd(ifreq,1,1) =  facd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(ima+(ispa-1)*ndim,ima+(ispa-1)*ndim,1)
     fw1_nd(ifreq,2,2) =  facd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(imb+(ispb-1)*ndim,imb+(ispb-1)*ndim,1)
     fw1_nd(ifreq,1,2) =  facnd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(ima+(ispa-1)*ndim,imb+(ispb-1)*ndim,1)
     fw1_nd(ifreq,2,1) =  facnd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(imb+(ispb-1)*ndim,ima+(ispa-1)*ndim,1)
     call xginv(fw1_nd(ifreq,:,:),2)
   end if
 end do
end subroutine testcode_ctqmc_b
!!***


!!****f* m_forctqmc/testcode_ctqmc
!! NAME
!! testcode_ctqmc
!!
!! FUNCTION
!! Setup ultra simple hybridization to test CTQMC in simple situations.
!!
!! INPUTS
!! gtmp_nd
!! gw_tmp_nd
!! temp = temperature
!! dmftqmc_l = number of times slices
!! nflavor = number of flavor
!! testrot = 0/1 if rotation of hybridization is tested or not
!! testcode = 1 if tests are activated.
!! opt = 1/2 if pre or postprocessing of CTQMC data.
!!
!! OUTPUT
!! fw1_nd = non diagonal hybridization
!! fw1 = hybridization
!! umod = value of U
!!
!!
!! SIDE EFFECTS
!!  gtmp_nd
!!  gw_tmp_nd
!!
!! NOTES
!!
!! SOURCE

subroutine testcode_ctqmc(dmftqmc_l,fw1_nd,fw1,gtmp_nd,gw_tmp_nd,levels_ctqmc,hybri_limit,&
&   nflavor,opt,temp,testrot,testcode,umod)


!Arguments ------------------------------------
!scalars
 integer, intent(in) :: dmftqmc_l,nflavor,testrot,testcode,opt
 real(dp), intent(in) :: temp
 real(dp), intent(out) :: umod(2,2)
 complex(dp), intent(inout) :: gw_tmp_nd(:,:,:)
 real(dp),  intent(inout) :: gtmp_nd(:,:,:)
 complex(dp), intent(out) :: fw1(:,:)
 complex(dp), intent(out) :: fw1_nd(:,:,:)
 real(dp),  intent(inout) :: levels_ctqmc(:)
 complex(dp),  intent(inout) :: hybri_limit(:,:)

!Local variables ------------------------------
 character(len=500) :: message
 integer :: ifreq, itau,realrot,simplehyb
 real(dp) :: omega
 real(dp) :: tbi1,tbi2,e2,tbi3,tbi4,e3,e4,tbi21,tbi12,e3b,e4b,tbi21b,tbi12b
 complex(dp) :: e1
! arrays
 complex(dp) :: RR(2,2)
 complex(dp) :: RR1(2,2)
 complex(dp) :: RRi(2,2)
 complex(dp) :: RRt(2,2)
! ************************************************************************
 if (testcode==0) return
 if (nflavor/=2) then
   write(message,'(2a)') ch10,' testcode nflavor.ne.2'
   ABI_ERROR(message)
 end if

 simplehyb=2
 simplehyb=1
 simplehyb=3
 !=========================
 ! Built rotation matrix
 !=========================
 realrot=0
 realrot=2
 if (realrot==1) then
   ! Real rotation
   !=========================
   RR(1,1)  =  SQRT(3.d0)/2.d0
   RR(1,2)  = -1.d0/2.d0
   RR(2,1)  =  1.d0/2.d0
   RR(2,2)  =  SQRT(3.d0)/2.d0
 else if (realrot==2) then
   ! Real rotation
   !=========================
   RR(1,1)  =  SQRT(1.d0/2.d0)
   RR(1,2)  = -SQRT(1.d0/2.d0)
   RR(2,1)  =  SQRT(1.d0/2.d0)
   RR(2,2)  =  SQRT(1.d0/2.d0)
 else
   ! Complex rotation
   !=========================
   RR(1,1)  =  CMPLX(one,two)
   RR(1,2)  =  CMPLX(one,one)
   RR(2,1)  =  CMPLX(one,-one)
   RR(2,2)  =  CMPLX(-one,two)
   RR=RR/sqrt(seven)
 end if
 ! Check rotation is unitary
 !==========================
 RRi(1,1) =  conjg(RR(1,1))
 RRi(1,2) =  conjg(RR(2,1))
 RRi(2,1) =  conjg(RR(1,2))
 RRi(2,2) =  conjg(RR(2,2))
 RR1(:,:)  = MATMUL ( RR(:,:) , RRi(:,:)          )
 !write(6,*) "RR1",RR1
 if(abs(RR1(1,1)-one).gt.tol7.or.abs(RR1(1,2)).gt.tol7.or.abs(RR1(2,2)-one).gt.tol7.or.abs(RR1(2,1)).gt.tol7) then
   write(message,'(2a)') ch10,' testcode error in rotation matrix'
   ABI_ERROR(message)
 end if


 !=================================
 ! Built hybridization  for CTQMC
 !=================================
 if (opt==1) then

 !  Parameters: tight-binding + U
 !  firt test of the code try umod=0, and (tbi1,tbi2,e1,e2)=(2,1,0.5,0.0) testrot=1
 !  second test of the code try umod=four, and (tbi1,tbi2,e1,e2)=(2,1,0.0,0.0) testrot=1
 !=======================================================================================
   fw1_nd(:,:,:)= czero
   tbi1=2.0_dp
   tbi2=1.0_dp
   tbi3=1.0_dp
   tbi4=1.0_dp
   tbi12=2.5_dp
   tbi12b=2.5_dp
   tbi21=2.5_dp
   tbi21b=2.5_dp
   e1=cmplx(0.0,0.0,8)
   e2=zero
   e3=0.2
   e4=0.3
   e3b=0.3
   e4b=-0.2
   umod(:,:)=0.d0

   if(testrot==1.and.(abs(tbi1-tbi2)<tol6)) then
     write(message,'(3a)') ch10,' testrot=1 with tbi1=tbi2 is equivalent' &
     ,'to testrot=0: change testrot'
     ABI_WARNING(message)
   end if
   ! Built fw1_nd
   !==============
   do ifreq=1,dmftqmc_l

     omega=pi*temp*(two*float(ifreq)-1)

     if(simplehyb==1) then
       fw1_nd(ifreq,1,1) =  -umod(1,1)/two+tbi1**2/(dcmplx(0.d0,omega)-e1)
       fw1_nd(ifreq,2,2) =  -umod(1,1)/two+tbi2**2/(dcmplx(0.d0,omega)-e2)
       fw1(ifreq,1)      =  -umod(1,1)/two+tbi1**2/(dcmplx(0.d0,omega)-e1)
       fw1(ifreq,2)      =  -umod(1,1)/two+tbi2**2/(dcmplx(0.d0,omega)-e2)
       hybri_limit(1,1)=tbi1**2
       hybri_limit(2,2)=tbi2**2
       hybri_limit(1,2)=0.d0
       hybri_limit(2,1)=0.d0
     else if(simplehyb==2) then
       fw1_nd(ifreq,1,1) =  -umod(1,1)/two+tbi1**2/(dcmplx(0.d0,omega)-e1)+tbi3**2/(dcmplx(0.d0,omega)-e3)
       fw1_nd(ifreq,2,2) =  -umod(1,1)/two+tbi2**2/(dcmplx(0.d0,omega)-e2)+tbi4**2/(dcmplx(0.d0,omega)-e4)
       fw1(ifreq,1)      =  -umod(1,1)/two+tbi1**2/(dcmplx(0.d0,omega)-e1)
       fw1(ifreq,2)      =  -umod(1,1)/two+tbi2**2/(dcmplx(0.d0,omega)-e2)
     else if(simplehyb==3) then
       fw1_nd(ifreq,1,1) =  -umod(1,1)/two+tbi1**2/(dcmplx(0.d0,omega)-e1)
       fw1_nd(ifreq,2,2) =  -umod(1,1)/two+tbi2**2/(dcmplx(0.d0,omega)-e2)
       fw1_nd(ifreq,1,2) =  tbi12**2/(dcmplx(0.d0,omega)-e3)+tbi12b**2/(dcmplx(0.d0,omega)-e3b)
       fw1_nd(ifreq,2,1) =  tbi21**2/(dcmplx(0.d0,omega)-e4)+tbi21b**2/(dcmplx(0.d0,omega)-e4b)
       fw1(ifreq,1)      =  -umod(1,1)/two+tbi1**2/(dcmplx(0.d0,omega)-e1)
       fw1(ifreq,2)      =  -umod(1,1)/two+tbi2**2/(dcmplx(0.d0,omega)-e2)
       hybri_limit(1,1)=tbi1**2
       hybri_limit(2,2)=tbi2**2
       hybri_limit(1,2)=tbi12**2+tbi12b**2
       hybri_limit(2,1)=tbi21**2+tbi21b**2
     end if
     write(132,*) omega,real(fw1_nd(ifreq,1,1)),aimag(fw1_nd(ifreq,1,1))
     write(133,*) omega,real(fw1_nd(ifreq,1,2)),aimag(fw1_nd(ifreq,1,2))
     write(134,*) omega,real(fw1_nd(ifreq,2,1)),aimag(fw1_nd(ifreq,2,1))
     write(135,*) omega,real(fw1_nd(ifreq,2,2)),aimag(fw1_nd(ifreq,2,2))
     write(1234,*) omega, real(fw1(ifreq,1)),aimag(fw1(ifreq,1))
   end do
   ! Built level and limit of hybridization
   !=======================================
   levels_ctqmc(1:nflavor)=-umod(1,1)/two

   write(std_out,*) "fw1_nd"
   write(std_out,*) fw1_nd(1,1,1), fw1_nd(1,1,2)
   write(std_out,*) fw1_nd(1,2,1), fw1_nd(1,2,2)
   write(std_out,*) "fw1"
   write(std_out,*) fw1(1,1), fw1(1,2)
   write(std_out,*) fw1(2,1), fw1(2,2)

 ! Rotate hybridization if testrot=1
 !==================================
   if(testrot==1) then

     do ifreq=1,dmftqmc_l
       RRt(:,:)  = MATMUL ( RR(:,:)  , fw1_nd(ifreq,:,:) )
   !write(6,*) "RRt"
   !write(6,*) RRt(1,1), RRt(1,2)
   !write(6,*) RRt(2,1), RRt(2,2)
       RR1(:,:)  = MATMUL ( RRt(:,:) , RRi(:,:)          )
   !write(6,*) "RR1"
   !write(6,*) RR1(1,1), RR1(1,2)
   !write(6,*) RR1(2,1), RR1(2,2)
       fw1_nd(ifreq,:,:)=RR1(:,:)
       omega=pi*temp*(two*float(ifreq)+1)
       write(3322,*) omega,real(fw1_nd(ifreq,1,1)),aimag(fw1_nd(ifreq,1,1))
       write(232,*) omega,real(fw1_nd(ifreq,1,1)),aimag(fw1_nd(ifreq,1,1))
       write(233,*) omega,real(fw1_nd(ifreq,1,2)),aimag(fw1_nd(ifreq,1,2))
       write(234,*) omega,real(fw1_nd(ifreq,2,1)),aimag(fw1_nd(ifreq,2,1))
       write(235,*) omega,real(fw1_nd(ifreq,2,2)),aimag(fw1_nd(ifreq,2,2))
     end do

     ! Rotate limit of hybridization
     !=======================================
     RRt(:,:)  = MATMUL ( RR(:,:)  , hybri_limit(:,:)  )
     RR1(:,:)  = MATMUL ( RRt(:,:) , RRi(:,:)          )
     hybri_limit(:,:)=RR1(:,:)

   end if
   ! rajouter test real(fw1_nd(1,:,:)) doit etre diagonale

 !======================================
 ! Rotate Green's function from CTQMC
 !======================================
 else if(opt==2) then

   write(std_out,*) "gw_tmp_nd"
   write(std_out,*) gw_tmp_nd(1,1,1), gw_tmp_nd(1,1,2)
   write(std_out,*) gw_tmp_nd(1,2,1), gw_tmp_nd(1,2,2)
   ! Rotate Green's function back
   !==============================
   if(testrot==1) then
     do ifreq=1,dmftqmc_l
       RRt(1:nflavor,1:nflavor) = MATMUL ( RRi(1:nflavor,1:nflavor),gw_tmp_nd(ifreq,1:nflavor,1:nflavor) )
       RR1(1:nflavor,1:nflavor) = MATMUL ( RRt(1:nflavor,1:nflavor),RR(1:nflavor,1:nflavor) )
       gw_tmp_nd(ifreq,1:nflavor,1:nflavor)=RR1(1:nflavor,1:nflavor)
     end do

     write(std_out,*) "gw_tmp_nd after rotation"
     write(std_out,*) gw_tmp_nd(1,1,1), gw_tmp_nd(1,1,2)
     write(std_out,*) gw_tmp_nd(1,2,1), gw_tmp_nd(1,2,2)

     do itau=1,dmftqmc_l
       RRt(1:nflavor,1:nflavor) = MATMUL ( RRi(1:nflavor,1:nflavor),gtmp_nd(itau,1:nflavor,1:nflavor) )
       RR1(1:nflavor,1:nflavor)  = MATMUL ( RRt(1:nflavor,1:nflavor),RR(1:nflavor,1:nflavor) )
       gtmp_nd(itau,1:nflavor,1:nflavor)=real(RR1(1:nflavor,1:nflavor))
     end do

   ! Rotate Green's function for comparison with testrot=1
   !======================================================
   else if (testrot==0) then ! produce rotated green's function to compare to testrot=1 case

     do itau=1,dmftqmc_l
       RRt(1:nflavor,1:nflavor) = MATMUL ( RR(1:nflavor,1:nflavor),gtmp_nd(itau,1:nflavor,1:nflavor) )
       RR1(1:nflavor,1:nflavor)  = MATMUL ( RRt(1:nflavor,1:nflavor),RRi(1:nflavor,1:nflavor) )
       write(444,*) real(itau-1)/(temp*real(dmftqmc_l)),real(RR1(1,1)),real(RR1(2,2)),real(RR1(1,2)),real(RR1(2,1))
     end do

   end if

   ! Print out rotated Green's function
   !=====================================
   do itau=1,dmftqmc_l
     write(555,'(e14.5,4(2e14.5,3x))') real(itau-1)/(temp*real(dmftqmc_l)),gtmp_nd(itau,1,1),&
&     gtmp_nd(itau,2,2),gtmp_nd(itau,1,2),gtmp_nd(itau,2,1)
   end do

   write(message,'(2a)') ch10,' testcode end of test calculation'
   ABI_ERROR(message)

 end if
 close(444)
 close(555)

end subroutine testcode_ctqmc
!!***

!!****f* m_forctqmc/ctqmcoutput_to_green
!! NAME
!! ctqmcoutput_to_green
!!
!! FUNCTION
!!  Put values of green function from ctqmc into green datatype
!!  Symetrize over spin if calculation is non magnetic
!!
!! INPUTS
!!  paw_dmft <type(paw_dmft_type)>= DMFT data structure
!!  gtmp_nd(dmftqmc_l,nflavor,nflavor) = Green's fct in imag time (with off diag terms)
!!  gw_tmp_nd(nb_of_frequency,nflavor,nflavor) = Green's fct in imag freq (with off diag terms)
!!  gtmp(dmftqmc_l,nflavor) = Green's fct in imag time (diag)
!!  gw_tmp(nb_of_frequency,nflavor+1) =Green's fct in imag freq (diag)
!!  iatom = atoms on which the calculation has been done
!!  leg_measure = logical, to Legendre Measurement or not (if done Green function is frequency is computed)
!!  opt_nondiag = integer, it activated, then
!!
!! OUTPUT
!!  green <type(green_type)>= green's function
!!
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

subroutine ctqmcoutput_to_green(green,paw_dmft,gtmp_nd,gw_tmp_nd,gtmp,gw_tmp,iatom,leg_measure,opt_nondiag)

!Arguments ------------------------------------
!scalars
 type(paw_dmft_type), intent(in)  :: paw_dmft
 type(green_type), intent(inout) :: green
 real(dp), allocatable, intent(in) :: gtmp_nd(:,:,:)
 complex(dp), allocatable, intent(in) :: gw_tmp(:,:)
 complex(dp), allocatable, intent(in) :: gw_tmp_nd(:,:,:)
 real(dp), allocatable, intent(in) :: gtmp(:,:)
 integer, intent(in) :: iatom,opt_nondiag
 logical(kind=1), intent(in) :: leg_measure
 character(len=500) :: message

!Local variables ------------------------------
 integer :: ifreq, itau,im1,im2,isppol,ispinor1,ispinor2,iflavor1
 integer :: iflavor2,tndim,ispinor,iflavor,im,nflavor
! ************************************************************************
 tndim=2*paw_dmft%lpawu(iatom)+1
 nflavor=2*(tndim)

 do itau=1,paw_dmft%dmftqmc_l
   green%oper_tau(itau)%matlu(iatom)%mat(:,:,:)=czero
 end do
 green%occup_tau%matlu(iatom)%mat(nflavor:,:,:)=czero

 do ifreq=1,paw_dmft%dmft_nwlo
   green%oper(ifreq)%matlu(iatom)%mat(:,:,:)=czero
 end do
 green%occup%matlu(iatom)%mat(:,:,:)=czero

!   built time and frequency green's function from output of CTQMC
! =================================================================
 if(opt_nondiag==1) then
   do isppol=1,paw_dmft%nsppol
     do ispinor1=1,paw_dmft%nspinor
       do im1=1,tndim
         iflavor1=im1+tndim*(ispinor1-1)+tndim*(isppol-1)
         do ispinor2=1,paw_dmft%nspinor
           do im2=1,tndim
             iflavor2=im2+tndim*(ispinor2-1)+tndim*(isppol-1)
             do itau=1,paw_dmft%dmftqmc_l
               green%oper_tau(itau)%matlu(iatom)%mat(im1+(ispinor1-1)*tndim,im2+(ispinor2-1)*tndim,isppol)=&
&               gtmp_nd(itau,iflavor1,iflavor2)
               ! symetrize over spin if nsppol=nspinor=1
               if(paw_dmft%nsppol==1.and.paw_dmft%nspinor==1) then
                 green%oper_tau(itau)%matlu(iatom)%mat(im1+(ispinor1-1)*tndim,im2+(ispinor2-1)*tndim,isppol)=&
&                 (gtmp_nd(itau,iflavor1,iflavor2)+gtmp_nd(itau,iflavor1+tndim,iflavor2+tndim))/two
               end if
             end do  !itau
             if(paw_dmft%dmft_solv<6.or.leg_measure) then
               do ifreq=1,paw_dmft%dmft_nwlo
                 green%oper(ifreq)%matlu(iatom)%mat(im1+(ispinor1-1)*tndim,im2+(ispinor2-1)*tndim,isppol)=&
&                 gw_tmp_nd(ifreq,iflavor1,iflavor2)
               ! symetrize over spin if nsppol=nspinor=1
                 if(paw_dmft%nsppol==1.and.paw_dmft%nspinor==1) then
                   green%oper(ifreq)%matlu(iatom)%mat(im1+(ispinor1-1)*tndim,im2+(ispinor2-1)*tndim,isppol)=&
&                   (gw_tmp_nd(ifreq,iflavor1,iflavor2)+&
&                   gw_tmp_nd(ifreq,iflavor1+tndim,iflavor2+tndim))/two
                 end if
               end do ! ifreq
             end if
           end do  ! im2
         end do  ! ispinor2
       end do  ! im1
     end do  ! ispinor
   end do ! isppol
 else
   iflavor=0
   do isppol=1,paw_dmft%nsppol
     do ispinor=1,paw_dmft%nspinor
       do im=1,tndim
         iflavor=iflavor+1
         do itau=1,paw_dmft%dmftqmc_l
           green%oper_tau(itau)%matlu(iatom)%mat(im+(ispinor-1)*tndim,im+(ispinor-1)*tndim,isppol)=gtmp(itau,iflavor)
           ! symetrize over spin if nsppol=paw_dmft%nspinor=1
           if(paw_dmft%nsppol==1.and.paw_dmft%nspinor==1) then
             green%oper_tau(itau)%matlu(iatom)%mat(im+(ispinor-1)*tndim,im+(ispinor-1)*tndim,isppol)=&
&             (gtmp(itau,iflavor)+gtmp(itau,iflavor+tndim))/two
           end if
         end do
!       ifreq2=0
         do ifreq=1,paw_dmft%dmft_nwlo
!         if(paw_dmft%select_log(ifreq)==1) then
!           ifreq2=ifreq2+1
           green%oper(ifreq)%matlu(iatom)%mat(im+(ispinor-1)*tndim,im+(ispinor-1)*tndim,isppol)=gw_tmp(ifreq,iflavor)
           ! symetrize over spin if nsppol=paw_dmft%nspinor=1
           if(paw_dmft%nsppol==1.and.paw_dmft%nspinor==1) then
             green%oper(ifreq)%matlu(iatom)%mat(im+(ispinor-1)*tndim,im+(ispinor-1)*tndim,isppol)=&
&             (gw_tmp(ifreq,iflavor)+gw_tmp(ifreq,iflavor+tndim))/two
           end if
         end do
       end do
     end do
   end do
 end if
 if(paw_dmft%nsppol==1.and.paw_dmft%nspinor==1) then
   write(message,'(a,2x,a,f13.5)') ch10,&
&   " == nsppol==1 and nspden==1: Green functions from CTQMC have been symetrized over spin"
   call wrtout(std_out,message,'COLL')
 end if

end subroutine ctqmcoutput_to_green
!!***

!!****f* m_forctqmc/ctqmcoutput_printgreen
!! NAME
!! ctqmcoutput_printgreen
!!
!! FUNCTION
!!  Print values of green function in files.
!!  Symetrize imaginary time Green's function in a peculiar case
!!  (dmft_solv=8 and natom=1). Should be moved later.
!!
!! INPUTS
!!  paw_dmft <type(paw_dmft_type)>= DMFT data structure
!!  gtmp_nd(dmftqmc_l,nflavor,nflavor) = Green's fct in imag time (with off diag terms)
!!  gw_tmp_nd(nb_of_frequency,nflavor,nflavor) = Green's fct in imag freq (with off diag terms)
!!  gtmp(dmftqmc_l,nflavor) = Green's fct in imag time (diag)
!!  gw_tmp(nb_of_frequency,nflavor+1) =Green's fct in imag freq (diag)
!!  iatom = atoms on which the calculation has been done
!!
!! OUTPUT
!!
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

subroutine ctqmcoutput_printgreen(paw_dmft,gtmp_nd,gw_tmp_nd,gtmp,gw_tmp,iatom)

!Arguments ------------------------------------
!scalars
 type(paw_dmft_type), intent(in)  :: paw_dmft
 real(dp), allocatable, intent(inout) :: gtmp_nd(:,:,:)
 complex(dp), allocatable, intent(in) :: gw_tmp(:,:)
 complex(dp), allocatable, intent(in) :: gw_tmp_nd(:,:,:)
 real(dp), allocatable, intent(in) :: gtmp(:,:)
 integer, intent(in) :: iatom

!Local variables ------------------------------
 character(len=500) :: message
 integer :: ifreq, itau,iflavor1
 integer :: tndim,iflavor,nflavor
 character(len=2) :: gtau_iter,iatomnb
 integer :: unt
! ************************************************************************
 tndim=2*paw_dmft%lpawu(iatom)+1
 nflavor=2*(tndim)
 !----------------------------------------
 ! <DEBUG>
 !----------------------------------------
 ! Construct UNIT
 if(paw_dmft%idmftloop < 10) then
   write(gtau_iter,'("0",i1)') paw_dmft%idmftloop
 elseif(paw_dmft%idmftloop >= 10 .and. paw_dmft%idmftloop < 100) then
   write(gtau_iter,'(i2)') paw_dmft%idmftloop
 else
   gtau_iter="xx"
 end if
 if(iatom < 10) then
   write(iatomnb,'("0",i1)') iatom
 elseif(iatom >= 10 .and. iatom < 100) then
   write(iatomnb,'(i2)') iatom
 else
   iatomnb='xx'
 end if

 if(paw_dmft%myproc .eq. mod(paw_dmft%nproc+1,paw_dmft%nproc)) then
! < HACK >
  if(paw_dmft%dmft_solv==6.or.paw_dmft%dmft_solv==7) then
    if (open_file(trim(paw_dmft%filapp)//"_atom_"//iatomnb//"_Gw_"//gtau_iter//".dat", message, newunit=unt) /=0) then
      ABI_ERROR(message)
    end if
    do ifreq=1,paw_dmft%dmft_nwli
      write(unt,'(29f21.14)') paw_dmft%omega_lo(ifreq),((gw_tmp_nd(ifreq,iflavor,iflavor)), iflavor=1, nflavor)
    end do
    close(unt)
  else
    if(paw_dmft%dmft_solv==5) then
      if (open_file(trim(paw_dmft%filapp)//"_atom_"//iatomnb//"_Gtau_"//gtau_iter//".dat", message, newunit=unt) /= 0) then
        ABI_ERROR(message)
      end if
      do itau=1,paw_dmft%dmftqmc_l
        write(unt,'(29f21.14)') float(itau-1)/float(paw_dmft%dmftqmc_l)/paw_dmft%temp,&
        (gtmp(itau,iflavor), iflavor=1, nflavor)
      end do
      write(unt,'(29f21.14)') 1/paw_dmft%temp, (-1_dp-gtmp(1,iflavor), iflavor=1, nflavor)
      close(unt)
    endif
    if(paw_dmft%dmft_solv==8) then
      if (open_file(trim(paw_dmft%filapp)//"_atom_"//iatomnb//"_Gtau_offdiag_unsym_"//gtau_iter//".dat",&
&      message, newunit=unt) /= 0) then
        ABI_ERROR(message)
      end if
      do itau=1,paw_dmft%dmftqmc_l
        write(unt,'(196f21.14)') float(itau-1)/float(paw_dmft%dmftqmc_l)/paw_dmft%temp,&
        ((gtmp_nd(itau,iflavor,iflavor1), iflavor=1, nflavor),iflavor1=1, nflavor)
      end do
      close(unt)
!      if(paw_dmft%natom==1) then ! If natom>1, it should be moved outside the loop over atoms
!        ABI_MALLOC(matlu1,(paw_dmft%natom))
!        call init_matlu(paw_dmft%natom,paw_dmft%nspinor,paw_dmft%nsppol,paw_dmft%lpawu,matlu1)
!        do itau=1,paw_dmft%dmftqmc_l
!          do isppol=1,paw_dmft%nsppol
!            do ispinor1=1,paw_dmft%nspinor
!              do im1=1,tndim
!                iflavor1=im1+tndim*(ispinor1-1)+tndim*(isppol-1)
!                do ispinor2=1,paw_dmft%nspinor
!                  do im2=1,tndim
!                    iflavor2=im2+tndim*(ispinor2-1)+tndim*(isppol-1)
!                    matlu1(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)=&
!&                     gtmp_nd(itau,iflavor1,iflavor2)
!                  end do  ! im2
!                end do  ! ispinor2
!              end do  ! im1
!            end do  ! ispinor
!          end do ! isppol
!          call rotate_matlu(matlu1,eigvectmatlu,paw_dmft%natom,3,0)
!          call slm2ylm_matlu(matlu1,paw_dmft%natom,2,0)
!          call sym_matlu(cryst_struc,matlu1,pawang,paw_dmft)
!          call slm2ylm_matlu(matlu1,paw_dmft%natom,1,0)
!          call rotate_matlu(matlu1,eigvectmatlu,paw_dmft%natom,3,1)
!          do isppol=1,paw_dmft%nsppol
!            do ispinor1=1,paw_dmft%nspinor
!              do im1=1,tndim
!                iflavor1=im1+tndim*(ispinor1-1)+tndim*(isppol-1)
!                do ispinor2=1,paw_dmft%nspinor
!                  do im2=1,tndim
!                    iflavor2=im2+tndim*(ispinor2-1)+tndim*(isppol-1)
!                    gtmp_nd(itau,iflavor1,iflavor2)=&
!                     matlu1(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)
!                  end do  ! im2
!                end do  ! ispinor2
!              end do  ! im1
!            end do  ! ispinor
!          end do ! isppol
!        end do  !itau
!        call destroy_matlu(matlu1,paw_dmft%natom)
!        ABI_FREE(matlu1)
!      endif ! if natom=1
      if (open_file(trim(paw_dmft%filapp)//"_atom_"//iatomnb//"_Gtau_offdiag_"//gtau_iter//".dat",&
&      message, newunit=unt) /= 0) then
        ABI_ERROR(message)
      end if
      do itau=1,paw_dmft%dmftqmc_l
        write(unt,'(196f21.14)') float(itau-1)/float(paw_dmft%dmftqmc_l)/paw_dmft%temp,&
        ((gtmp_nd(itau,iflavor,iflavor1), iflavor=1, nflavor),iflavor1=1, nflavor)
      end do
      close(unt)
    endif
    !open(unit=4243, file=trim(paw_dmft%filapp)//"_atom_"//iatomnb//"_F_"//gtau_iter//".dat")
    !call BathOperator_printF(paw_dmft%hybrid(iatom)%hybrid%bath,4243) !Already comment here
    !close(4243)
    if(paw_dmft%dmft_solv==5) then
      if (open_file(trim(paw_dmft%filapp)//"_atom_"//iatomnb//"_Gw_"//gtau_iter//".dat", message, newunit=unt) /= 0) then
        ABI_ERROR(message)
      end if
      do ifreq=1,paw_dmft%dmft_nwlo
        write(unt,'(29f21.14)') paw_dmft%omega_lo(ifreq), &
&        (gw_tmp(ifreq,iflavor), iflavor=1, nflavor)
      end do
    endif
    close(unt)
  end if
! </ HACK >
  end if


end subroutine ctqmcoutput_printgreen
!!***

!!****f* m_forctqmc/ctqmc_calltriqs
!! NAME
!! ctqmc_calltriqs
!!
!! FUNCTION
!!  Call TRIQS solver and perform calculation of Green's function using
!!  Legendre coefficients.
!!
!! INPUTS
!!  paw_dmft <type(paw_dmft_type)>= DMFT data structure
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  hu <type(hu_type)>= U interaction
!!  levels_ctqmc(nflavor) = atomic levels
!!  gw_tmp_nd(nb_of_frequency,nflavor,nflavor) = Green's fct in imag freq (with off diag terms)
!!  gtmp_nd(dmftqmc_l,nflavor,nflavor) = Green's fct in imag time (with off diag terms)
!!  fw1_nd(dmft_nwlo,nflavor,nflavor) = Hybridization fct in imag time (with off diag terms)
!!  leg_measure = logical, true is legendre measurement is activated
!!  iatom= index of atom
!!
!! OUTPUT
!!
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

subroutine ctqmc_calltriqs(paw_dmft,cryst_struc,hu,levels_ctqmc,gtmp_nd,gw_tmp_nd,fw1_nd,leg_measure,iatom)

#if defined HAVE_TRIQS_v4_0 || defined HAVE_TRIQS_v2_0 || defined HAVE_TRIQS_v1_4
 use TRIQS_CTQMC !Triqs module
#endif
#if defined HAVE_PYTHON_INVOCATION
 use m_invoke_python
#endif
 use, intrinsic :: iso_c_binding

!Arguments ------------------------------------
!scalars
 type(paw_dmft_type), intent(in)  :: paw_dmft
 type(crystal_t),intent(in) :: cryst_struc
 type(hu_type), intent(in) :: hu(cryst_struc%ntypat)
 real(dp), allocatable, target, intent(inout) :: gtmp_nd(:,:,:)
 complex(dp), allocatable, target, intent(inout) :: gw_tmp_nd(:,:,:)
 complex(dp), allocatable, target, intent(in) :: fw1_nd(:,:,:)
 real(dp), allocatable, target, intent(inout) ::  levels_ctqmc(:)
 logical(kind=1), intent(in) :: leg_measure
 integer, intent(in) :: iatom

!Local variables ------------------------------
 complex(dp), allocatable, target ::fw1_nd_tmp(:,:,:)
 complex(dp), allocatable, target :: g_iw(:,:,:)
 real(dp), allocatable, target :: u_mat_ij(:,:)
 real(dp), allocatable, target :: u_mat_ijkl(:,:,:,:)
 real(dp), allocatable, target :: u_mat_ijkl_tmp(:,:,:,:)
 real(dp), allocatable, target :: gl_nd(:,:,:)
 type(c_ptr) :: levels_ptr, fw1_nd_ptr, u_mat_ij_ptr, u_mat_ijkl_ptr, g_iw_ptr, gtau_ptr, gl_ptr
 real(dp), allocatable :: jbes(:)
 character(len=500) :: message
 integer :: ifreq, iflavor1
 integer :: iflavor2,iflavor,nflavor,iflavor3,itypat
 integer :: nfreq,ntau,nleg,ileg
 integer :: verbosity_solver ! min 0 -> max 3
 logical(kind=1) :: rot_inv = .false.
#if defined HAVE_TRIQS_v2_0 || defined HAVE_TRIQS_v1_4 || defined HAVE_PYTHON_INVOCATION
 logical(kind=1) :: hist = .false.
 logical(kind=1) :: wrt_files = .true.
 logical(kind=1) :: tot_not = .true.
#endif
 real(dp) :: beta,besp,bespp,xx
 complex(dp) :: u_nl

#if defined HAVE_PYTHON_INVOCATION
!----------
!Variables for writing out the NETCDF file when calling PYTHON_INVOCATION
!----------
 integer(kind=4) :: ncid
 integer(kind=4) :: dim_one_id, dim_nflavor_id, dim_nwlo_id, dim_nwli_id
 integer(kind=4) :: dim_qmc_l_id, dim_nleg_id
 integer(kind=4), dimension(2) :: dim_u_mat_ij_id
 integer(kind=4), dimension(3) :: dim_fw1_id, dim_g_iw_id, dim_gl_id, dim_gtau_id
 integer(kind=4), dimension(4) :: dim_u_mat_ijkl_id
 integer(kind=4) :: var_rot_inv_id, var_leg_measure_id, var_hist_id, var_wrt_files_id
 integer(kind=4) :: var_tot_not_id, var_n_orbitals_id, var_n_freq_id, var_n_tau_id, var_n_l_id, var_n_cycles_id
 integer(kind=4) :: var_cycle_length_id, var_ntherm_id, var_verbo_id, var_seed_id, var_beta_id
 integer(kind=4) :: var_levels_id, var_u_mat_ij_id, var_u_mat_ijkl_id, var_real_fw1_nd_id, var_imag_fw1_nd_id
 integer(kind=4) :: var_real_g_iw_id, var_imag_g_iw_id, var_gtau_id, var_gl_id, var_spacecomm_id

 integer :: itau

 integer(kind=4) :: varid
 logical :: file_exists
 complex :: i
 character(len=100) :: filename

 real(dp), allocatable, target :: new_re_g_iw(:,:,:), new_im_g_iw(:,:,:)
 real(dp), allocatable, target :: new_g_tau(:,:,:), new_gl(:,:,:)
!----------
#endif
! ************************************************************************

 ! fw1_nd: Hybridation
 ! levels_ctqmc: niveaux
 ! hu(itypat)%udens(:,:) : U_ij
 ! hu(itypat)%u(:,:,:,:) : uijkl
 ! temperature : paw_dmft%temp
 ! paw_dmft%dmftqmc_l: nombre de points en temps -1
 ! paw_dmft%dmftqmc_n: nombre de cycles
 ! ?? Quelles sorties: Les fonctions de Green
 ! frequence/temps/Legendre.
 ! Double occupations ?? <n_i n_j>
 ! test n_tau > 2*nfreq => ntau = 2*nfreq + 1
 !   for non diagonal code:
 !   call CtqmcInterface_run(hybrid,fw1_nd(1:paw_dmft%dmftqmc_l,:,:),Gtau=gtmp_nd,&
 !&  Gw=gw_tmp_nd,D=Doccsum,E=green%ecorr_qmc(iatom),&
 !&  Noise=Noise,matU=hu(itypat)%udens,opt_levels=levels_ctqmc,hybri_limit=hybri_limit)
 !Check choice of user to fix model bool var for the solver
 if (paw_dmft%dmft_solv==6) then
   rot_inv = .false.
 else !obviously paw_dmft%dmft_solv==7 with rot invariant terms
   rot_inv = .true.
 end if

 nfreq = paw_dmft%dmft_nwli
 !paw_dmft%dmft_nwlo = paw_dmft%dmft_nwli !transparent for user
 ntau  = paw_dmft%dmftqmc_l !(2*paw_dmft%dmftqmc_l)+1 !nfreq=paw_dmft%dmft_nwli
 nleg  = paw_dmft%dmft_triqs_nleg
 nflavor=2*(2*paw_dmft%lpawu(iatom)+1)
 itypat=cryst_struc%typat(iatom)


 verbosity_solver = paw_dmft%prtvol
 beta = 1.0/(paw_dmft%temp*Ha_eV)

 !Allocation in/output array phase:
 ABI_MALLOC(fw1_nd_tmp,(1:nflavor,1:nflavor,1:nfreq)) !column major
 ABI_MALLOC(g_iw,(1:nflavor,1:nflavor,1:nfreq)) !column major
 ABI_MALLOC(u_mat_ij,(1:nflavor,1:nflavor)) !column major
 ABI_MALLOC(u_mat_ijkl,(1:nflavor,1:nflavor,1:nflavor,1:nflavor)) !column major
 ABI_MALLOC(u_mat_ijkl_tmp,(1:nflavor,1:nflavor,1:nflavor,1:nflavor)) !column major

 if ( leg_measure ) then !only if functionality is enabled
   ABI_MALLOC(gl_nd,(1:nleg,1:nflavor,1:nflavor)) !column major !nl = 30 by default
 end if

 !Conversion datas Ha -> eV (some duplications for test...)
 !fw1_nd_tmp = fw1_nd(1:paw_dmft%dmftqmc_l,:,:) * Ha_eV !fw1_nd = fw1_nd * Ha_eV !Ok?

 do iflavor=1,nflavor
   do iflavor1=1,nflavor
     do ifreq=1,nfreq
       fw1_nd_tmp(iflavor,iflavor1,ifreq) = fw1_nd(ifreq,iflavor,iflavor1) * Ha_eV
!        WRITE(500,*) "[IN Fortran] F[ w= ",ifreq," l= ",iflavor," l_= ",iflavor1,"] = ",fw1_nd(ifreq,iflavor,iflavor1)
     end do
   end do
 end do

    !Report test
!    WRITE(502,*) hu(itypat)%udens
!    do ifreq=1,paw_dmft%dmftqmc_l
!      write(501,*) ((fw1_nd(ifreq,iflavor,iflavor1),iflavor=1,nflavor),iflavor1=1,nflavor)
!    enddo
    !write(866,*)paw_dmft%dmft_nwlo,paw_dmft%dmftqmc_l
    !write(866,*) u_mat_ij
!   do iflavor=1,nflavor+1
!     do iflavor1=1,nflavor+1
!       WRITE(502,*) "[OUT Fortran] U(i,j)[ l= ",iflavor," l_= ",iflavor1,"] = ",hu(itypat)%udens(iflavor,iflavor1)
!     enddo
!   enddo

!          if(paw_dmft%myproc==0) then
!          do iflavor=1,nflavor
!            do iflavor1=1,nflavor
!               do iflavor2=1,nflavor
!                  do iflavor3=1,nflavor
!                    write(490,*), hu(itypat)%vee(iflavor,iflavor1,iflavor2,iflavor3)
!                  enddo
!                 enddo
!                enddo
!              enddo
!          endif

!          if(paw_dmft%myproc==0) then
!          do iflavor=1,nflavor
!            do iflavor1=1,nflavor
!            write(491,*), hu(itypat)%udens(iflavor,iflavor1) !(1,1,1,1)
!            enddo
!          enddo
!          endif

!          do iflavor=1,nflavor
!            do iflavor1=1,nflavor
!          do iflavor2=1,nflavor
!            do iflavor3=1,nflavor
                 ! WRITE(552,*), hu(itypat)%vee!(iflavor,iflavor1,iflavor2,iflavor3)
!            enddo
!          enddo
!            enddo
!          enddo

 call vee_ndim2tndim_hu_r(paw_dmft%lpawu(iatom),dble(hu(itypat)%vee),u_mat_ijkl_tmp,1)
 do iflavor=1,nflavor
   do iflavor1=1,nflavor
     do iflavor2=1,nflavor
       do iflavor3=1,nflavor
         u_mat_ijkl(iflavor,iflavor1,iflavor2,iflavor3)   =  Ha_eV * u_mat_ijkl_tmp(iflavor,iflavor1,iflavor2,iflavor3)
       end do
     end do
   end do
 end do

 !u_mat_ijkl   =  Ha_eV * reshape( u_mat_ijkl , [nflavor,nflavor,nflavor,nflavor] )  !column -> row major + conversion
 u_mat_ij     = transpose( dble(hu(itypat)%udens) ) * Ha_eV !column -> row major + conversion
 levels_ctqmc = levels_ctqmc * Ha_eV

 !Location array in memory for C++ pointer args to pass
 !----------------------------------------------------
 g_iw_ptr       = C_LOC( gw_tmp_nd ) !C_LOC( g_iw )
 gtau_ptr       = C_LOC( gtmp_nd ) !C_LOC( gtau )
 gl_ptr         = C_LOC( gl_nd )
 fw1_nd_ptr     = C_LOC( fw1_nd_tmp )
 u_mat_ij_ptr   = C_LOC( u_mat_ij )
 u_mat_ijkl_ptr = C_LOC( u_mat_ijkl )
 levels_ptr     = C_LOC( levels_ctqmc )

 !Calling interfaced TRIQS solver subroutine from src/67_triqs_ext package
 if (paw_dmft%dmft_solv==9) then
#ifndef HAVE_PYTHON_INVOCATION
  write(message,'(23a)') ch10,' Python invocation flag requiered! You need to install ABINIT with ',&
   'enable_python_invocation = yes" in your "configure.ac" file.'
  call wrtout(std_out,message,'COLL')
  ABI_ERROR(message)
#else
  ! Creating the NETCDF file
  ! write(std_out, "(a)") trim(paw_dmft%filapp)
  write(filename, '(a, a)') trim(paw_dmft%filnamei), "_abinit_output_for_py.nc"
  write(std_out, '(3a)') ch10, "    Creating NETCDF file: ", trim(filename)
  NCF_CHECK(nf90_create(filename, NF90_CLOBBER, ncid))

  ! Defining the dimensions of the variables to write in the NETCDF file
  NCF_CHECK(nf90_def_dim(ncid, "one", 1, dim_one_id))
  NCF_CHECK(nf90_def_dim(ncid, "nflavor", nflavor, dim_nflavor_id))
  NCF_CHECK(nf90_def_dim(ncid, "nwlo", paw_dmft%dmft_nwlo, dim_nwlo_id))
  NCF_CHECK(nf90_def_dim(ncid, "nwli", paw_dmft%dmft_nwli, dim_nwli_id))
  NCF_CHECK(nf90_def_dim(ncid, "qmc_l", paw_dmft%dmftqmc_l, dim_qmc_l_id))
  NCF_CHECK(nf90_def_dim(ncid, "nleg", nleg, dim_nleg_id))

  dim_u_mat_ij_id = (/ dim_nflavor_id, dim_nflavor_id /)
  dim_u_mat_ijkl_id = (/ dim_nflavor_id, dim_nflavor_id, dim_nflavor_id, dim_nflavor_id /)
  dim_fw1_id = (/ dim_nflavor_id, dim_nflavor_id, dim_nwli_id /)
  dim_g_iw_id = (/ dim_nwli_id, dim_nflavor_id, dim_nflavor_id /)
  dim_gtau_id = (/ dim_qmc_l_id, dim_nflavor_id, dim_nflavor_id /)
  dim_gl_id = (/ dim_nleg_id, dim_nflavor_id, dim_nflavor_id /)

  ! Defining the variables
  NCF_CHECK(nf90_def_var(ncid, "rot_inv",         NF90_INT, dim_one_id,           var_rot_inv_id))
  NCF_CHECK(nf90_def_var(ncid, "leg_measure",     NF90_INT, dim_one_id,           var_leg_measure_id))
  NCF_CHECK(nf90_def_var(ncid, "hist",            NF90_INT, dim_one_id,           var_hist_id))
  NCF_CHECK(nf90_def_var(ncid, "wrt_files",       NF90_INT, dim_one_id,           var_wrt_files_id))
  NCF_CHECK(nf90_def_var(ncid, "tot_not",         NF90_INT, dim_one_id,           var_tot_not_id))
  NCF_CHECK(nf90_def_var(ncid, "n_orbitals",      NF90_INT, dim_one_id,           var_n_orbitals_id))
  NCF_CHECK(nf90_def_var(ncid, "n_freq",          NF90_INT, dim_one_id,           var_n_freq_id))
  NCF_CHECK(nf90_def_var(ncid, "n_tau",           NF90_INT, dim_one_id,           var_n_tau_id))
  NCF_CHECK(nf90_def_var(ncid, "n_l",             NF90_INT, dim_one_id,           var_n_l_id))
  NCF_CHECK(nf90_def_var(ncid, "n_cycles",        NF90_INT, dim_one_id,           var_n_cycles_id))
  NCF_CHECK(nf90_def_var(ncid, "cycle_length",    NF90_INT, dim_one_id,           var_cycle_length_id))
  NCF_CHECK(nf90_def_var(ncid, "ntherm",          NF90_INT, dim_one_id,           var_ntherm_id))
  NCF_CHECK(nf90_def_var(ncid, "verbo",           NF90_INT, dim_one_id,           var_verbo_id))
  NCF_CHECK(nf90_def_var(ncid, "seed",            NF90_INT, dim_one_id,           var_seed_id))
  NCF_CHECK(nf90_def_var(ncid, "beta",            NF90_FLOAT, dim_one_id,         var_beta_id))
  NCF_CHECK(nf90_def_var(ncid, "levels",          NF90_DOUBLE, dim_nflavor_id,    var_levels_id))
  NCF_CHECK(nf90_def_var(ncid, "u_mat_ij",        NF90_DOUBLE, dim_u_mat_ij_id,   var_u_mat_ij_id))
  NCF_CHECK(nf90_def_var(ncid, "u_mat_ijkl",      NF90_DOUBLE, dim_u_mat_ijkl_id, var_u_mat_ijkl_id))
  NCF_CHECK(nf90_def_var(ncid, "real_fw1_nd",     NF90_DOUBLE, dim_fw1_id,        var_real_fw1_nd_id))
  NCF_CHECK(nf90_def_var(ncid, "imag_fw1_nd",     NF90_DOUBLE, dim_fw1_id,        var_imag_fw1_nd_id))
  NCF_CHECK(nf90_def_var(ncid, "real_g_iw",       NF90_DOUBLE, dim_g_iw_id,       var_real_g_iw_id))
  NCF_CHECK(nf90_def_var(ncid, "imag_g_iw",       NF90_DOUBLE, dim_g_iw_id,       var_imag_g_iw_id))
  NCF_CHECK(nf90_def_var(ncid, "gtau",            NF90_DOUBLE, dim_gtau_id,       var_gtau_id))
  NCF_CHECK(nf90_def_var(ncid, "gl",              NF90_DOUBLE, dim_gl_id,         var_gl_id))
  NCF_CHECK(nf90_def_var(ncid, "spacecomm",       NF90_INT, dim_one_id,           var_spacecomm_id))
  NCF_CHECK(nf90_enddef(ncid))

  ! Filling the variables with actual data
  if (rot_inv) then
   NCF_CHECK(nf90_put_var(ncid, var_rot_inv_id,       1))
  else
   NCF_CHECK(nf90_put_var(ncid, var_rot_inv_id,       0))
  end if
  if (leg_measure) then
   NCF_CHECK(nf90_put_var(ncid, var_leg_measure_id,   1))
  else
   NCF_CHECK(nf90_put_var(ncid, var_leg_measure_id,   0))
  end if
  if (hist) then
   NCF_CHECK(nf90_put_var(ncid, var_hist_id,          1))
  else
   NCF_CHECK(nf90_put_var(ncid, var_hist_id,          0))
  end if
  if (wrt_files) then
   NCF_CHECK(nf90_put_var(ncid, var_wrt_files_id,     1))
  else
   NCF_CHECK(nf90_put_var(ncid, var_wrt_files_id,     0))
  end if
  if (tot_not) then
   NCF_CHECK(nf90_put_var(ncid, var_tot_not_id,       1))
  else
   NCF_CHECK(nf90_put_var(ncid, var_tot_not_id,       0))
  end if
  NCF_CHECK(nf90_put_var(ncid, var_n_orbitals_id,         nflavor))
  NCF_CHECK(nf90_put_var(ncid, var_n_freq_id,             nfreq))
  NCF_CHECK(nf90_put_var(ncid, var_n_tau_id,              ntau))
  NCF_CHECK(nf90_put_var(ncid, var_n_l_id,                nleg))
  NCF_CHECK(nf90_put_var(ncid, var_n_cycles_id,           int(paw_dmft%dmftqmc_n/paw_dmft%nproc)))
  NCF_CHECK(nf90_put_var(ncid, var_cycle_length_id,       paw_dmft%dmftctqmc_meas*2*2*nflavor))
  NCF_CHECK(nf90_put_var(ncid, var_ntherm_id,             paw_dmft%dmftqmc_therm))
  NCF_CHECK(nf90_put_var(ncid, var_verbo_id,              verbosity_solver))
  NCF_CHECK(nf90_put_var(ncid, var_seed_id,               paw_dmft%dmftqmc_seed))
  NCF_CHECK(nf90_put_var(ncid, var_beta_id,               beta))
  NCF_CHECK(nf90_put_var(ncid, var_levels_id,             levels_ctqmc))
  NCF_CHECK(nf90_put_var(ncid, var_u_mat_ij_id,           u_mat_ij))
  NCF_CHECK(nf90_put_var(ncid, var_u_mat_ijkl_id,         u_mat_ijkl))
  NCF_CHECK(nf90_put_var(ncid, var_real_fw1_nd_id,        real(fw1_nd_tmp)))
  NCF_CHECK(nf90_put_var(ncid, var_imag_fw1_nd_id,        aimag(fw1_nd_tmp)))
  NCF_CHECK(nf90_put_var(ncid, var_real_g_iw_id,          real(gw_tmp_nd)))
  NCF_CHECK(nf90_put_var(ncid, var_imag_g_iw_id,          aimag(gw_tmp_nd)))
  NCF_CHECK(nf90_put_var(ncid, var_gtau_id,               gtmp_nd))
  NCF_CHECK(nf90_put_var(ncid, var_gl_id,                 gl_nd))
  NCF_CHECK(nf90_put_var(ncid, var_spacecomm_id,          paw_dmft%spacecomm))
  NCF_CHECK(nf90_close(ncid))

  write(std_out, '(4a)') ch10, "    NETCDF file ", trim(filename), " written; Launching python invocation"

  ! Invoking python to execute the script
  call invoke_python_run_script (0, paw_dmft%myproc, trim(paw_dmft%filnamei), paw_dmft%spacecomm)
  ! call Invoke_python_triqs (paw_dmft%myproc, trim(paw_dmft%filnamei)//c_null_char)
  call xmpi_barrier(paw_dmft%spacecomm)
  call flush_unit(std_out)

  ! Allocating the fortran variables for the results
  ABI_MALLOC(new_re_g_iw,(nflavor,nflavor, paw_dmft%dmft_nwli))
  ABI_MALLOC(new_im_g_iw,(nflavor,nflavor, paw_dmft%dmft_nwli))
  ABI_MALLOC(new_g_tau,(nflavor,nflavor, paw_dmft%dmftqmc_l))
  ABI_MALLOC(new_gl,(nflavor,nflavor, nleg))
  i = (0, 1)

  ! Check if file exists
  write(filename, '(a, a)') trim(paw_dmft%filnamei), "_py_output_for_abinit.nc"

  INQUIRE(FILE=filename, EXIST=file_exists)
  if(.not. file_exists) then
   write(message,'(4a)') ch10,' Cannot find file ', trim(filename), '! Make sure the python script writes it with the right name and at the right place!'
   call wrtout(std_out,message,'COLL')
   ABI_ERROR(message)
  endif

  write(std_out, '(3a)') ch10, "    Reading NETCDF file ", trim(filename)

  ! Opening the NETCDF file
  NCF_CHECK(nf90_open(filename, nf90_nowrite, ncid))

  ! Read from the file
  ! Re{G_iw}
  write(std_out, '(2a)') ch10, "    -- Re[G(iw_n)]"
  NCF_CHECK(nf90_inq_varid(ncid, "re_g_iw", varid))
  NCF_CHECK(nf90_get_var(ncid, varid, new_re_g_iw))
  ! Im{G_iw}
  write(std_out, '(2a)') ch10, "    -- Im[G(iw_n)]"
  NCF_CHECK(nf90_inq_varid(ncid, "im_g_iw", varid))
  NCF_CHECK(nf90_get_var(ncid, varid, new_im_g_iw))
  ! G_tau
  write(std_out, '(2a)') ch10, "    -- G(tau)"
  NCF_CHECK(nf90_inq_varid(ncid, "g_tau", varid))
  NCF_CHECK(nf90_get_var(ncid, varid, new_g_tau))
  ! G_l
  write(std_out, '(2a)') ch10, "    -- G_l"
  NCF_CHECK(nf90_inq_varid(ncid, "gl", varid))
  NCF_CHECK(nf90_get_var(ncid, varid, new_gl))

  ! Assigning data
  do iflavor1=1, nflavor
   do iflavor2=1, nflavor
    do ifreq=1, paw_dmft%dmft_nwli
     gw_tmp_nd(ifreq, iflavor1, iflavor2) = new_re_g_iw(iflavor1, iflavor2, ifreq) &
&               + i*new_im_g_iw(iflavor1, iflavor2, ifreq)
    end do
    do itau=1, paw_dmft%dmftqmc_l
     gtmp_nd(itau, iflavor1, iflavor2) = new_g_tau(iflavor1, iflavor2, itau)
    end do
    do ileg=1, nleg
     gl_nd(ileg, iflavor1, iflavor2) = new_gl(iflavor1, iflavor2, ileg)
    end do
   end do
  end do

  ! Deallocating
  ABI_FREE(new_re_g_iw)
  ABI_FREE(new_im_g_iw)
  ABI_FREE(new_g_tau)
  ABI_FREE(new_gl)
#endif
 elseif(paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7) then
  !Calling interfaced TRIQS solver subroutine from src/01_triqs_ext package
  !----------------------------------------------------
#if defined HAVE_TRIQS_v2_0 || defined HAVE_TRIQS_v1_4
 call Ctqmc_triqs_run (     rot_inv, leg_measure, hist, wrt_files, tot_not,   &
&  nflavor, nfreq, ntau , nleg, int(paw_dmft%dmftqmc_n/paw_dmft%nproc),       &
&  paw_dmft%dmftctqmc_meas*2*2*nflavor, paw_dmft%dmftqmc_therm,               &
&  verbosity_solver, paw_dmft%dmftqmc_seed,beta,                              &
&  levels_ptr,  u_mat_ij_ptr, u_mat_ijkl_ptr, fw1_nd_ptr,                     &
!&  g_iw_ptr, gtau_ptr, gl_ptr, paw_dmft%spacecomm                             )
&  g_iw_ptr, gtau_ptr, gl_ptr, paw_dmft%myproc                             )
#endif
 endif

  !WRITE(*,*) "Hello Debug"
  !call xmpi_barrier(paw_dmft%spacecomm) !Resynch all processus after calling Impurity solver from TRIQS

  !Report output datas from TRIQS to Abinit
  !Interacting G(iw)
 ! OG Commented these loops because they are useless
 !do ifreq=1,nfreq
 !  do iflavor1=1,nflavor
 !    do iflavor=1,nflavor
 !   !   gw_tmp_nd(ifreq,iflavor,iflavor1) = g_iw(iflavor,iflavor1,ifreq) !* Ha_eV !because 1/ G0(eV)
 !   !  WRITE(503,*) "[OUT Fortran] G(iw)[ w= ",ifreq," l= ",iflavor," l_= ",iflavor1,"] = ",gw_tmp_nd(ifreq,iflavor,iflavor1)!g_iw(iflavor,iflavor1,ifreq)
 !    end do
 !  end do
 !end do

! Convert in Ha
 gw_tmp_nd = gw_tmp_nd*Ha_eV

!     do iflavor1=1,nflavor
!       do iflavor=1,nflavor
!
!        WRITE(510,*) "[OUT Fortran] U[ l= ",iflavor," l_= ",iflavor1,"] = ",u_mat_ij(iflavor,iflavor1)
!       enddo
!     enddo

! if(paw_dmft%myproc==0) write(6,*) "essai",paw_dmft%myproc, gw_tmp_nd(2,1,1)
! if(paw_dmft%myproc==1) write(6,*) "essai",paw_dmft%myproc,gw_tmp_nd(2,1,1)
! if(paw_dmft%myproc==0) write(621,*) "essai",paw_dmft%myproc, gw_tmp_nd(2,1,1)
! if(paw_dmft%myproc==1) write(622,*) "essai",paw_dmft%myproc,gw_tmp_nd(2,1,1)
! call flush_unit(621)
! call flush_unit(622)
! write(message,*) ch10, "essai",paw_dmft%myproc, paw_dmft%myproc,paw_dmft%dmftqmc_seed!gw_tmp_nd(2,1,1)
! call wrtout(555,message,'PERS',.true.)
! if(paw_dmft%myproc==0) write(499,*) "essai",paw_dmft%myproc, paw_dmft%dmftqmc_seed
! if(paw_dmft%myproc==1) write(498,*) "essai",paw_dmft%myproc,paw_dmft%dmftqmc_seed

  !Its associated G(tau): Problem of compatibility => paw_dmft%dmftqmc_l < (2*paw_dmft%dmftqmc_l)+1 => We report only  paw_dmft%dmftqmc_l =  first values of G(tau)...
!   do iflavor=1,nflavor
!     do iflavor1=1,nflavor
!       do itau=1,ntau
!         if ( modulo(itau,2) == 1 ) then !Problem of binding: paw_dmft%dmftqmc_l =! ntau => We take one value by 2 and Write in file all the G(tau) out function from TRIQS
          !gtmp_nd(itau,iflavor,iflavor1) = gtau(iflavor,iflavor1,itau)
!         endif
!         if(paw_dmft%myproc==0) then
!           WRITE(504,*) "[OUT Fortran] G[ tau= ",itau," l= ",iflavor," l_= ",iflavor1,"] = ",gtmp_nd(itau,iflavor,iflavor1) !gtmp_nd(itau,iflavor,iflavor1) !passage ok avec ntau/iflavor1/iflavor (iflavor,iflavor1,ntau)
!         endif
!       enddo
!     enddo
!   enddo

  ! Write Legendre Polynoms G(L) for extrapolation of Interacting G(iw) by FT, only if leg_measure == TRUE
  ! -------------------------------------------------------------------------------------------
 if (leg_measure) then
   do ileg=1,nleg
     WRITE(505,*) ileg,((gl_nd(ileg,iflavor,iflavor1),iflavor=1,nflavor),iflavor1=1,nflavor)
   end do
   close(505)
 end if
! f(paw_dmft%myproc==0) then
!  do itau=1,paw_dmft%dmftqmc_l
!    write(490,*) ((gtmp_nd(itau,iflavor,iflavor1),iflavor=1,nflavor),iflavor1=1,nflavor)
!  enddo
! ndif
 ABI_FREE( fw1_nd_tmp )
 ABI_FREE( g_iw )
 ABI_FREE( u_mat_ijkl )
 ABI_FREE( u_mat_ijkl_tmp )
 ABI_FREE( u_mat_ij )


  !  Compute Green's function in imaginary freq using Legendre coefficients
  ! -----------------------------------------------------------------------
 if (leg_measure) then
   call xmpi_barrier(paw_dmft%spacecomm)
   call flush_unit(std_out)
   write(message,'(2a)') ch10,"    ==  Compute G(iw_n) from Legendre coefficients"
   call wrtout(std_out,message,'COLL')
   ABI_MALLOC( jbes, (nleg))
   gw_tmp_nd=czero

  !   write(77,*) " TEST OF BESSEL S ROUTINES 0 0"

  !   xx=0_dp
  !   ileg=0
  !   call sbf8(ileg+1,xx,jbes)
  !   write(77,*) "T0 A",jbes(ileg+1)
  !   call jbessel(jbes(ileg+1),besp,bespp,ileg,1,xx)
  !   write(77,*) "T0 B",jbes(ileg+1)
  !   write(77,*) "T0 C",bessel_jn(ileg,xx)

  !   write(77,*) " TEST OF BESSEL S ROUTINES 1.5 0"

  !   xx=1.5_dp
  !   ileg=0
  !   call sbf8(ileg+1,xx,jbes)
  !   write(77,*) "T1 A",jbes(ileg+1)
  !   call jbessel(jbes(ileg+1),besp,bespp,ileg,1,xx)
  !   write(77,*) "T1 B",jbes(ileg+1)
  !   write(77,*) "T1 C",bessel_jn(ileg,xx)

  !   write(77,*) " TEST OF BESSEL S ROUTINES 1.5 1"

  !   xx=1.5_dp
  !   ileg=1
  !   call sbf8(ileg+1,xx,jbes)
  !   write(77,*) "T2 A",jbes(ileg+1)
  !   call jbessel(jbes(ileg+1),besp,bespp,ileg,1,xx)
  !   write(77,*) "T2 B",jbes(ileg+1)
  !   write(77,*) "T2 C",bessel_jn(ileg,xx)


   do ifreq=1,paw_dmft%dmft_nwli
     xx=real(2*ifreq-1,kind=dp)*pi/two
     if(xx<=100_dp) call sbf8(nleg,xx,jbes)
     do ileg=1,nleg
    ! write(77,*) "A",ifreq,jbes(ileg),xx

       if(xx>=99) call jbessel(jbes(ileg),besp,bespp,ileg-1,1,xx)
    ! write(77,*) "B",ifreq,jbes(ileg),xx

     !write(77,*) "C",ifreq,jbes(ileg),xx

       u_nl=sqrt(float(2*ileg-1))*(-1)**(ifreq-1)*cmplx(0_dp,one)**(ileg)*jbes(ileg)
      write(77,*) "----------",ileg,jbes(ileg), u_nl,gl_nd(ileg,1,1)

       do iflavor=1,nflavor
         do iflavor1=1,nflavor
           gw_tmp_nd(ifreq,iflavor,iflavor1)= gw_tmp_nd(ifreq,iflavor,iflavor1) + &
&           u_nl*gl_nd(ileg,iflavor,iflavor1)
         end do
       end do

  !    write(77,*) "------------------", gw_tmp_nd(ifreq,1,1)

     end do
  !  write(77,*) "------------------ sum ", gw_tmp_nd(ifreq,1,1)
   end do
   ABI_FREE( jbes )
   call xmpi_barrier(paw_dmft%spacecomm)
   call flush_unit(std_out)
 end if
 gw_tmp_nd = gw_tmp_nd*Ha_eV


 if ( leg_measure ) then !only if functionality is enabled
   ABI_FREE(gl_nd)
 end if


end subroutine ctqmc_calltriqs
!!***

!!****f* m_forctqmc/ctqmc_calltriqs_c
!! NAME
!! ctqmc_calltriqs_c
!!
!! FUNCTION
!! This routines calls TRIQS/CTHYB using the C++ API in order
!! to solve the impurity model.
!!
!! INPUTS
!!  paw_dmft <type(paw_dmft_type)>= DMFT data structure
!!  green <type(green_type)>= green's function
!!  self <type(self_type)>= self-energy
!!  hu <type(hu_type)>= U interaction
!!  weiss <type(green_type)>= inverse of weiss function
!!  self_new <type(self_type)>= impurity self-energy
!!  pawprtvol = flag for printing
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

subroutine ctqmc_calltriqs_c(paw_dmft,green,self,hu,weiss,self_new,pawprtvol)

#if defined HAVE_TRIQS_v4_0 || defined HAVE_TRIQS_v3_2
 use TRIQS_CTQMC
#endif
 use ISO_C_BINDING

!Arguments ------------------------------------
 integer, intent(in) :: pawprtvol
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(green_type), target, intent(inout) :: green,weiss
 type(self_type), intent(inout) :: self,self_new
 type(hu_type), intent(inout) :: hu(paw_dmft%ntypat)
!Local variables ------------------------------
 integer :: basis,i,iatom,iblock,iflavor,iflavor1,iflavor2,ifreq,ilam,ilam_prev,ileg,im,im1,integral,isppol,isub
 integer :: itau,itypat,iw,l,len_t,lpawu,myproc,natom,ncon,ndim,nflavor,nflavor_max,ngauss,nleg,nmoments
 integer :: nspinor,nsppol,nsub,ntau,ntot,nwlo,p,pad_elam,pad_lambda,read_data,rot_type_vee,tndim,unt,verbo,wdlr_size
 integer, target :: ndlr
 logical :: debug,density_matrix,entropy,leg_measure,lexist,nondiag,off_diag,rot_inv
 real(dp) :: besp,bespp,beta,dx,elam,emig_tot,err,err_,fact,fact2,tau,tol,xtau,xx
 complex(dp) :: mself_1,mself_2,occ_tmp,u_nl
 complex(dp), target :: eu
 type(oper_type) :: hdc_ctqmc
 type(oper_type), target :: energy_level
 type(self_type) :: hybmwdhyb
 type(c_ptr) :: block_ptr,eu_ptr,flavor_ptr,fname_data_ptr,fname_dataw_ptr,fname_histo_ptr,ftau_ptr,gl_ptr,gtau_ptr
 type(c_ptr) :: inner_ptr,levels_ptr,mself_1_ptr,mself_2_ptr,ndlr_ptr,occ_ptr,siz_ptr,udens_ptr,vee_ptr,wdlr_ptr
 integer, allocatable :: flavor_list(:,:,:),nblocks(:)
 integer, target, allocatable :: block_list(:,:),flavor_tmp(:,:),inner_list(:,:),siz_block(:,:)
 real(dp), allocatable :: adlr(:,:),bdlr(:),elam_list(:),emig(:),gl_dlr_re(:),gl_dlr_im(:),jbes(:),lam_list(:)
 real(dp), allocatable :: leg_array(:,:),moment_fit(:),t_lp(:,:),tpoints(:),tweights(:),wdlr(:),wdlr_beta(:,:)
 real(dp), target, allocatable :: wdlr_tmp(:)
 complex(dp), allocatable :: adlr_iw(:,:),gl_dlr(:,:,:,:),gl_tmp(:,:,:,:),gtau_dlr(:,:,:),gtau_leg(:,:,:),shift(:)
 complex(dp), target, allocatable :: gl(:,:,:),gtau(:,:,:),levels_ctqmc(:,:),moments_self_1(:),moments_self_2(:),occ(:)
 type(matlu_type), allocatable :: eigvectmatlu(:),matlu_tmp(:)
 type(matlu_type), target, allocatable :: dmat_ctqmc(:),ftau(:),udens_rot(:)
 type(matlu_type), pointer :: matlu_pt(:) => null()
 type(vee_type), target, allocatable :: vee_rot(:)
 character(len=1) :: tag_block4
 character(len=2) :: tag_block,tag_block3,tag_lam,tag_lam_prev
 character(len=4) :: tag_at
 character(len=14) :: tag_elam,tag_lambda
 character(len=500) :: stringfile,stringfile_prev,tag_block2,tag_lam2,tag_lam_prev2
 character(len=10000) :: message
 character(len=fnlen), target :: fname_data,fname_dataw,fname_histo
! ************************************************************************

 basis          = paw_dmft%dmftctqmc_basis
 beta           = one / paw_dmft%temp
 debug          = paw_dmft%dmft_triqs_debug
 density_matrix = paw_dmft%dmft_triqs_measure_density_matrix
 entropy        = (paw_dmft%dmft_triqs_entropy == 1)
 integral       = paw_dmft%dmft_triqs_compute_integral
 leg_measure    = paw_dmft%dmft_triqs_leg_measure
 myproc         = paw_dmft%myproc
 natom          = paw_dmft%natom
 nflavor_max    = 2 * (2*paw_dmft%maxlpawu+1)
 ngauss         = paw_dmft%dmft_triqs_gaussorder
 nleg           = paw_dmft%dmft_triqs_nleg
 nspinor        = paw_dmft%nspinor
 nsppol         = paw_dmft%nsppol
 nsub           = paw_dmft%dmft_triqs_nsubdivisions
 ntau           = paw_dmft%dmftqmc_l
 nwlo           = paw_dmft%dmft_nwlo
 off_diag       = paw_dmft%dmft_triqs_off_diag
 rot_inv        = (paw_dmft%dmft_solv == 7)
 tol            = paw_dmft%dmft_triqs_tol_block

 if (rot_inv) then
   write(message,'(a,3x,a)') ch10,"== Rotationally Invariant Terms Included"
 else
   write(message,'(a,3x,a)') ch10,"== Density-Density Terms Included"
 end if
 call wrtout(std_out,message,"COLL")

 ABI_MALLOC(block_list,(nflavor_max,natom))
 ABI_MALLOC(dmat_ctqmc,(natom))
 ABI_MALLOC(eigvectmatlu,(natom))
 ABI_MALLOC(flavor_list,(nflavor_max,nflavor_max,natom))
 ABI_MALLOC(ftau,(natom))
 ABI_MALLOC(inner_list,(nflavor_max,natom))
 ABI_MALLOC(matlu_tmp,(natom))
 ABI_MALLOC(nblocks,(natom))
 ABI_MALLOC(shift,(natom))
 ABI_MALLOC(siz_block,(nflavor_max,natom))
 ABI_MALLOC(udens_rot,(natom))
 ABI_MALLOC(vee_rot,(natom))

 call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),dmat_ctqmc(:))
 call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),eigvectmatlu(:))
 call init_matlu(natom,2,ntau,paw_dmft%lpawu(:),ftau(:))
 call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),matlu_tmp(:))
 call init_matlu(natom,2,1,paw_dmft%lpawu(:),udens_rot(:))

 call init_oper(paw_dmft,energy_level,opt_ksloc=2)

 call init_vee(paw_dmft,vee_rot(:))

 if (entropy .and. integral == 2) then
   call init_oper(paw_dmft,hdc_ctqmc,opt_ksloc=2)
   call copy_matlu(self%hdc%matlu(:),hdc_ctqmc%matlu(:),natom)
 end if

 call compute_levels(energy_level,self%hdc,paw_dmft)

 write(message,'(a,3x,a)') ch10,"== Print Occupation matrix in cubic basis"
 call wrtout(std_out,message,"COLL")
 call print_matlu(green%occup%matlu(:),natom,1)

 call compute_moments_loc(green,self,energy_level,weiss,0)

 ! Build hybridization and remove spurious 0th order moment
 do ifreq=1,nwlo
   shift(:) = cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)
   call shift_matlu(weiss%oper(ifreq)%matlu(:),natom,shift(:))
   call fac_matlu(weiss%oper(ifreq)%matlu(:),natom,-cone)
   call add_matlu(weiss%oper(ifreq)%matlu(:),energy_level%matlu(:),matlu_tmp(:),natom,-1)
   call add_matlu(matlu_tmp(:),weiss%moments(1)%matlu(:),weiss%oper(ifreq)%matlu(:),natom,-1)
 end do ! ifreq

 call zero_matlu(weiss%moments(1)%matlu(:),natom)

 write(message,'(a,3x,a)') ch10,"== Print Delta(iw) for first frequency in cubic basis"
 call wrtout(std_out,message,"COLL")
 call print_matlu(weiss%oper(1)%matlu(:),natom,1)

 if (basis == 0) then
   write(message,'(a,3x,a)') ch10,"== Switching to CTQMC basis: staying in cubic basis"
 else if (basis == 1) then
   write(message,'(a,3x,2a)') ch10,"== Switching to CTQMC basis: using basis that", &
                                 & " diagonalizes the electronic levels"
 else if (basis == 2) then
   write(message,'(a,3x,2a)') ch10,"== Switching to CTQMC basis: using basis that", &
                                 & " diagonalizes the occupation matrix"
 else if (basis == 3) then
   write(message,'(a,3x,a)') ch10,"== Switching to CTQMC basis: using Ylm basis"
 else if (basis == 4) then
   write(message,'(a,3x,a)') ch10,"== Switching to CTQMC basis: using JmJ basis"
 end if
 call wrtout(std_out,message,"COLL")

 if (basis == 1) then
   call checkdiag_matlu(energy_level%matlu(:),natom,tol,nondiag)
   if (.not. nondiag) then
     basis = 0
     write(message,'(a,3x,a)') ch10,"== Electronic levels are already diagonal: staying in the cubic basis"
   else
     write(message,'(a,3x,a)') ch10,"== Switching to Ylm basis first"
   end if ! nondiag
   call wrtout(std_out,message,"COLL")
 end if ! basis=1

 if (basis == 2) then
   call checkdiag_matlu(green%occup%matlu(:),natom,tol,nondiag)
   if (.not. nondiag) then
     basis = 0
     write(message,'(a,3x,a)') ch10,"== Occupation matrix is already diagonal: staying in the cubic basis"
   else
     write(message,'(a,3x,a)') ch10,"== Switching to Ylm basis first"
   end if ! not nondiag
   call wrtout(std_out,message,"COLL")
 end if ! basis=2

 call copy_matlu(green%occup%matlu(:),dmat_ctqmc(:),natom)

 if (basis > 0) then ! First switch to Ylm basis in every case
   call slm2ylm_matlu(energy_level%matlu(:),natom,paw_dmft,1,0)
   call slm2ylm_matlu(dmat_ctqmc(:),natom,paw_dmft,1,0)
   do i=2,weiss%nmoments-1
     call slm2ylm_matlu(weiss%moments(i)%matlu(:),natom,paw_dmft,1,0)
   end do ! i
   do ifreq=1,nwlo
     if (weiss%distrib%procf(ifreq) /= myproc) cycle
     call slm2ylm_matlu(weiss%oper(ifreq)%matlu(:),natom,paw_dmft,1,0)
   end do ! ifreq
   if (entropy .and. integral == 2) then
     call slm2ylm_matlu(hdc_ctqmc%matlu(:),natom,paw_dmft,1,0)
   end if
 end if ! basis>0

 if (basis == 1 .or. basis == 2) then
   ! Find block structure in Ylm basis and diagonalize for each block (extremely useful in the
   ! case of degenerate levels ; this ensures minimal mixing of Ylm and thus maximal number of subspaces)
   if (basis == 1) then
     matlu_pt => energy_level%matlu(:)
   else
     matlu_pt => dmat_ctqmc(:)
   end if ! basis

   if (pawprtvol >= 3) then
     if (basis == 1) then
       write(message,'(a,3x,a)') ch10,"== Print Energy levels in Ylm basis"
     else
       write(message,'(a,3x,a)') ch10,"== Print Occupation matrix in Ylm basis"
     end if
     call wrtout(std_out,message,"COLL")

     call print_matlu(matlu_pt(:),natom,1)
   end if ! pawprtvol>=3

   call find_block_structure(paw_dmft,block_list(:,:),inner_list(:,:), &
       & flavor_list(:,:,:),siz_block(:,:),nblocks(:),matlu_pt(:),natom,nflavor_max)
   call diag_block(matlu_pt)
   matlu_pt => null()

   ! Make sure every process has the same rotation matrix in case of degenerate levels
   call xmpi_matlu(eigvectmatlu(:),natom,paw_dmft%spacecomm,master=0,option=2)

   if (basis == 1) then
     call rotate_matlu(dmat_ctqmc(:),eigvectmatlu(:),natom,1)
   else
     call rotate_matlu(energy_level%matlu(:),eigvectmatlu(:),natom,1)
   end if ! basis
   do i=2,weiss%nmoments-1
     call rotate_matlu(weiss%moments(i)%matlu(:),eigvectmatlu(:),natom,1)
   end do ! i
   do ifreq=1,nwlo
     if (weiss%distrib%procf(ifreq) /= myproc) cycle
     call rotate_matlu(weiss%oper(ifreq)%matlu(:),eigvectmatlu(:),natom,1)
   end do ! ifreq
   if (entropy .and. integral == 2) then
     call rotate_matlu(hdc_ctqmc%matlu(:),eigvectmatlu(:),natom,1)
   end if
 end if ! basis=1 or 2

 if (basis == 4) then
   call ylm2jmj_matlu(energy_level%matlu(:),natom,1,paw_dmft)
   call ylm2jmj_matlu(dmat_ctqmc(:),natom,1,paw_dmft)
   do i=2,weiss%nmoments-1
     call ylm2jmj_matlu(weiss%moments(i)%matlu(:),natom,1,paw_dmft)
   end do
   do ifreq=1,nwlo
     if (weiss%distrib%procf(ifreq) /= myproc) cycle
     call ylm2jmj_matlu(weiss%oper(ifreq)%matlu(:),natom,1,paw_dmft)
   end do ! ifreq
   if (entropy .and. integral == 2) then
     call ylm2jmj_matlu(hdc_ctqmc%matlu(:),natom,1,paw_dmft)
   end if
 end if ! basis=4

 if (basis == 0) then
   do iatom=1,natom
     lpawu = paw_dmft%lpawu(iatom)
     if (lpawu == -1) cycle
     itypat = paw_dmft%typat(iatom)
     udens_rot(iatom)%mat(:,:,1) = hu(itypat)%udens(:,:)
     vee_rot(iatom)%mat(:,:,:,:) = hu(itypat)%veeslm2(:,:,:,:)
   end do ! iatom
 else
   call gather_oper(weiss%oper(:),weiss%distrib,paw_dmft,opt_ksloc=2)
   rot_type_vee = 4
   if (basis == 3) rot_type_vee = 2
   if (basis == 4) rot_type_vee = 3
   call rotatevee_hu(hu(:),paw_dmft,pawprtvol,eigvectmatlu(:),rot_type_vee,udens_rot(:),vee_rot(:))
 end if ! basis

 write(message,'(a,3x,a)') ch10,"== Print Energy levels in CTQMC basis"
 call wrtout(std_out,message,"COLL")
 call print_matlu(energy_level%matlu(:),natom,1)

 write(message,'(a,3x,a)') ch10,"== Print Occupation matrix in CTQMC basis"
 call wrtout(std_out,message,"COLL")
 call print_matlu(dmat_ctqmc(:),natom,1)

 write(message,'(a,3x,a)') ch10,"== Print Delta(iw) for first frequency in CTQMC basis"
 call wrtout(std_out,message,"COLL")
 call print_matlu(weiss%oper(1)%matlu(:),natom,1)

 ! Possibly set the imaginary part and off-diagonal elements to 0 now that we
 ! are in the CTQMC basis. This is extremely important to do it explicitly instead of
 ! simply sending the real part or the diagonal elements to TRIQS since this modifies
 ! the electronic levels and hybridization that are used in Dyson's equation later.

#ifndef HAVE_TRIQS_COMPLEX
 write(message,'(a,3x,2a)') ch10,"== The imaginary part of Delta(tau) and the ", &
                         & "electronic levels is now set to 0"
 call wrtout(std_out,message,"COLL")
 ! Symmetrizing Delta(iw) is equivalent to neglecting the imaginary part of Delta(tau)
 err = zero
 do ifreq=1,nwlo
   call symmetrize_matlu(weiss%oper(ifreq)%matlu(:),natom,err=err_)
   if (err_ > err) err = err_
 end do ! ifreq
 do i=2,weiss%nmoments-1
   call symmetrize_matlu(weiss%moments(i)%matlu(:),natom,err=err_)
   if (err_ > err) err = err_
 end do ! i
 call zero_matlu(energy_level%matlu(:),natom,onlyimag=1,err=err_)
 if (err_ > err) err = err_
 if (err > tol) then
   write(message,'(2a)') "WARNING: This is not a good approximation ; the imaginary ", &
                      & "part is non negligible !"
   ABI_WARNING(message)
 end if ! err>tol
 if (entropy .and. integral == 2) then
   call zero_matlu(hdc_ctqmc%matlu(:),natom,onlyimag=1)
 end if
#endif

 err = zero
 if ((.not. rot_inv) .or. (.not. off_diag)) then
   call zero_matlu(energy_level%matlu(:),natom,onlynondiag=1,err=err_)
   if (err_ > err) err = err_
   if (entropy .and. integral == 2) then
     call zero_matlu(hdc_ctqmc%matlu(:),natom,onlynondiag=1)
   end if
 end if

 if (.not. off_diag) then
   write(message,'(a,3x,2a)') ch10,"== The off-diagonal elements of the hybridization ", &
                            & "and the electronic levels are now set to 0"
   call wrtout(std_out,message,"COLL")
   do ifreq=1,nwlo
     call zero_matlu(weiss%oper(ifreq)%matlu(:),natom,onlynondiag=1,err=err_)
     if (err_ > err) err = err_
   end do ! ifreq
   do i=2,weiss%nmoments-1
     call zero_matlu(weiss%moments(i)%matlu(:),natom,onlynondiag=1,err=err_)
     if (err_ > err) err = err_
   end do ! i
   if (err > tol) then
     write(message,'(2a)') "WARNING: This is not a good approximation ; the off-diagonal ", &
                         & "elements are non negligible !"
     ABI_WARNING(message)
   end if ! err>tol
 end if ! not off_diag

 ! Prepare DLR frequencies
 if (.not. leg_measure) then
   wdlr_size = 1000 ! make sure this is big enough
   ABI_MALLOC(wdlr_tmp,(wdlr_size))
   ndlr_ptr = C_LOC(ndlr)
   wdlr_ptr = C_LOC(wdlr_tmp)
#if defined HAVE_TRIQS_v4_0 || defined HAVE_TRIQS_v3_2
   call build_dlr(wdlr_size,ndlr_ptr,wdlr_ptr,paw_dmft%dmft_triqs_lambda,paw_dmft%dmft_triqs_epsilon)
#endif
   if (ndlr > wdlr_size) then
     write(message,'(a,i4,2a)') "You have more than ",wdlr_size," DLR frequencies.", &
                    & " Something is wrong here."
     ABI_ERROR(message)
   end if

   ABI_MALLOC(wdlr,(ndlr))
   ABI_MALLOC(wdlr_beta,(ndlr,4))
   wdlr(:) = wdlr_tmp(1:ndlr)
   wdlr_beta(:,1) = wdlr(:) / beta
   do i=2,4
     wdlr_beta(:,i) = wdlr_beta(:,i-1) * wdlr_beta(:,1)
   end do
   ABI_FREE(wdlr_tmp)
   ABI_MALLOC(adlr,(ndlr,ntau))
   ABI_MALLOC(adlr_iw,(ndlr,nwlo))

   do ifreq=1,nwlo
     do iw=1,ndlr
       adlr_iw(iw,ifreq) = k_iw(paw_dmft%omega_lo(ifreq),wdlr_beta(iw,1))
     end do ! iw
   end do ! ifreq

   do itau=1,ntau
     do iw=1,ndlr
       adlr(iw,itau) = k_it(dble(itau-1)/dble(ntau-1),wdlr(iw))
     end do ! iw
   end do ! itau

   write(tag_at,'(i4)') ndlr
   write(message,'(a,3x,3a)') ch10,"== There are ",trim(adjustl(tag_at))," DLR frequencies"
   call wrtout(std_out,message,"COLL")
   write(message,'(3x,1000(e10.3,2x))') wdlr_beta(:,1)
   call wrtout(std_out,message,"COLL")
   call identity_oper(green%moments(1),2)
 end if ! not leg_measure

 ! ntot is total number of lambda pts, + 1 is because we add the case lambda = 1 (which has no reason to be included in the
 ! Gauss-Legendre grid), since we need it for the rest of the SCF calculation
 ntot = merge(ngauss*nsub+1,1,integral>0.and.entropy)

 ABI_MALLOC(elam_list,(ntot))
 ABI_MALLOC(lam_list,(ntot)) ! scaling factors of U matrix for thermodynamic integration
 lam_list(ntot) = one
 green%integral = zero
 green%ekin_imp = zero

 ! Prepare Gauss-Legendre quadrature for thermodynamic integration over U
 if (integral > 0 .and. entropy) then

   ABI_MALLOC(tweights,(ngauss))
   ABI_MALLOC(tpoints,(ngauss))

   ! Calculation of Gauss-Legendre grid (on [-1,1]) of size ngauss
   call coeffs_gausslegint(-one,one,tpoints(:),tweights(:),ngauss)

   dx = one / dble(nsub)

   ! Split [0,1] into nsub intervals [x_i,x_{i+1}],i=1,nsub with x_i=(i-1)*dx
   do isub=1,nsub
     ! For each interval, the Gauss-Legendre grid of size ngauss is mapped from the
     ! t-world where t in [-1,1] to the x-world where x in [x_{isub},x_{isub+1}]
     ! We want an array sorted in descending order to optimize restart.
     ! No need to flip the tweights, as they are symmetric.
     lam_list(ntot-isub*ngauss:ntot-(isub-1)*ngauss-1) = (dble(isub)*two+tpoints(ngauss:1:-1)-one) * dx * half
   end do ! ilam

 end if ! integral and entropy

 ! Build most optimal block structure in CTQMC basis
 write(message,'(a,3x,2a)') ch10,"== Searching for the most optimal block structure of", &
                           & " the electronic levels and the hybridization"
 call wrtout(std_out,message,"COLL")

 call find_block_structure(paw_dmft,block_list(:,:),inner_list(:,:),flavor_list(:,:,:), &
               & siz_block(:,:),nblocks(:),energy_level%matlu(:),natom,nflavor_max,hyb=weiss)

 nmoments = weiss%nmoments - 2

  ! Inverse Fourier transform of the hybridization
 call fourier_inv(paw_dmft,nmoments,ntau,ftau(:),weiss%oper(:),weiss%moments(2:nmoments+1))

 if (entropy) then
   ABI_MALLOC(emig,(natom))
   ! Cubic splines to compute the derivative of Delta(iw)
   call initialize_self(hybmwdhyb,paw_dmft,opt_moments=1)
   call cubic_spline()
   call copy_matlu(energy_level%matlu(:),hybmwdhyb%moments(1)%matlu(:),natom)
   do i=2,weiss%nmoments-1
     call copy_matlu(weiss%moments(i)%matlu(:),hybmwdhyb%moments(i)%matlu(:),natom)
     call fac_matlu(hybmwdhyb%moments(i)%matlu(:),natom,cmplx(dble(i),zero,kind=dp))
   end do ! i
 end if ! entropy

 ! Solve impurity model for each atom
 do iatom=1,natom

   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   ndim    = 2*lpawu + 1
   tndim   = nspinor * ndim
   nflavor = 2 * ndim

   write(tag_at,'(i4)') iatom
   write(tag_block,'(i2)') nblocks(iatom)
   write(message,'(a,3x,6a)') ch10,"== Solving impurity model for atom ",trim(adjustl(tag_at)), &
                            & ", where there are ",trim(adjustl(tag_block))," blocks",ch10
   call wrtout(std_out,message,'COLL')

   do iblock=1,nblocks(iatom)
     write(tag_block,'(i2)') iblock - 1
     tag_block2 = ""
     do iflavor=1,siz_block(iblock,iatom)
       write(tag_block3,'(i2)') flavor_list(iflavor,iblock,iatom)
       tag_block2 = trim(tag_block2) // " " // trim(adjustl(tag_block3))
     end do ! iflavor
     tag_block4 = ""
     if (siz_block(iblock,iatom) > 1) tag_block4 = "s"
     write(message,'(2x,4a,1x,a)') "--> Block ",trim(adjustl(tag_block))," contains flavor",trim(adjustl(tag_block4)),trim(adjustl(tag_block2))
     call wrtout(std_out,message,'COLL')
   end do ! iblock

   write(message,'(a,3x,2a)') ch10,"== Schematic of the block structure",ch10
   call wrtout(std_out,message,'COLL')

   iflavor = 1
   do iblock=1,nblocks(iatom)
     do iflavor1=1,siz_block(iblock,iatom)
       tag_block2 = ""
       do iflavor2=1,iflavor-1
         tag_block2 = trim(tag_block2) // "  ."
       end do ! iflavor2
       do iflavor2=1,iflavor1-1
         tag_block2 = trim(tag_block2) // "  x"
       end do ! iflavor2
       write(tag_block,'(i2)') flavor_list(iflavor1,iblock,iatom)
       i = merge(1,2,flavor_list(iflavor1,iblock,iatom)>=10)
       tag_block2 = trim(tag_block2) // repeat(" ",i) // trim(adjustl(tag_block))
       do iflavor2=iflavor1+1,siz_block(iblock,iatom)
         tag_block2 = trim(tag_block2) // "  x"
       end do ! iflavor2
       do iflavor2=iflavor+siz_block(iblock,iatom),nflavor
         tag_block2 = trim(tag_block2) // "  ."
       end do ! iflavor2
       write(message,'(4x,a)') tag_block2
       call wrtout(std_out,message,'COLL')
     end do ! iflavor1
     iflavor = iflavor + siz_block(iblock,iatom)
   end do ! iblock

   call int2char4(iatom,tag_at)
   ABI_CHECK((tag_at(1:1)/='#'),'Bug: string length too short!')

   if (myproc == 0 .and. off_diag) then

     if (open_file(trim(paw_dmft%filapp)//"_Hybridization_offdiag_iatom"//tag_at//".dat",message,newunit=unt) /= 0) ABI_ERROR(message)
     write(unt,'(6a)') "# Off-diagonal components of Delta(tau) in the CTQMC basis",ch10, &
                     & "# Columns are ordered this way:",ch10, &
                     & "# Imaginary Time     ((Re(Delta_{ij}) Im(Delta_{ij}),i=1,2*(2*l+1)),j=1,2*(2*l+1)) where the", &
                     & " leftmost index varies first"

     do itau=1,ntau
       write(unt,'(2x,393(e18.10e3,2x))') beta*dble(itau-1)/dble(ntau-1), &
          & ((dble(ftau(iatom)%mat(im,im1,itau)),aimag(ftau(iatom)%mat(im,im1,itau)),im=1,nflavor),im1=1,nflavor)
     end do ! itau
     close(unt)

   end if ! myproc=0

   if (myproc == 0) then

     if (open_file(trim(paw_dmft%filapp)//"_Hybridization_diag_iatom"//tag_at//".dat",message,newunit=unt) /= 0) ABI_ERROR(message)
     write(unt,'(5a)') "# Diagonal components of Delta(tau) in the CTQMC basis",ch10, &
                     & "# Columns are ordered this way:",ch10, &
                     & "# Imaginary Time     (Delta_{ii},i=1,2*(2*l+1))"

     do itau=1,ntau
       write(unt,'(2x,393(e25.17e3,2x))') beta*dble(itau-1)/dble(ntau-1),(dble(ftau(iatom)%mat(im,im,itau)),im=1,nflavor)
     end do ! itau
     close(unt)

   end if ! myproc=0

   ABI_MALLOC(levels_ctqmc,(nflavor,nflavor))
   ABI_MALLOC(gtau,(ntau,nflavor,nflavor))
   gtau(:,:,:) = czero
   if (density_matrix) then
     ABI_MALLOC(occ,(nflavor))
   end if
   if ((.not. leg_measure) .and. density_matrix) then
     ABI_MALLOC(moments_self_1,(nflavor))
     ABI_MALLOC(moments_self_2,(nflavor))
   end if
   if (leg_measure) then
     ABI_MALLOC(gl,(nleg,nflavor,nflavor))
     gl(:,:,:) = czero
   end if
   ABI_MALLOC(flavor_tmp,(nflavor,nflavor))

   levels_ctqmc(:,:) = czero
   do isppol=1,nsppol
     levels_ctqmc(1+(isppol-1)*ndim:tndim+(isppol-1)*ndim,1+(isppol-1)*ndim:tndim+(isppol-1)*ndim) = &
        & energy_level%matlu(iatom)%mat(:,:,isppol)
     if (nsppol == 1 .and. nspinor == 1) levels_ctqmc(1+ndim:2*ndim,1+ndim:2*ndim) = levels_ctqmc(1:ndim,1:ndim)
   end do ! isppol

   ! Need to slice flavor_list to make it size nflavor*nflavor instead of size nflavor_max*nflavor_max
   flavor_tmp(:,:) = flavor_list(1:nflavor,1:nflavor,iatom)

   block_ptr       = C_LOC(block_list(:,iatom))
   eu_ptr          = C_LOC(eu)
   flavor_ptr      = C_LOC(flavor_tmp(:,:))
   fname_data_ptr  = C_LOC(fname_data)
   fname_dataw_ptr = C_LOC(fname_dataw)
   fname_histo_ptr = C_LOC(fname_histo)
   ftau_ptr        = C_LOC(ftau(iatom)%mat(:,:,:))
   gl_ptr          = C_LOC(gl(:,:,:))
   gtau_ptr        = C_LOC(gtau(:,:,:))
   inner_ptr       = C_LOC(inner_list(:,iatom))
   levels_ptr      = C_LOC(levels_ctqmc(:,:))
   mself_1_ptr     = C_LOC(moments_self_1(:))
   mself_2_ptr     = C_LOC(moments_self_2(:))
   occ_ptr         = C_LOC(occ(:))
   siz_ptr         = C_LOC(siz_block(:,iatom))
   udens_ptr       = C_LOC(udens_rot(iatom)%mat(:,:,1))
   vee_ptr         = C_LOC(vee_rot(iatom)%mat(:,:,:,:))

   verbo = 1

   do ilam=1,ntot

     if (ilam /= ntot) then
       write(message,'(a,3x,a,f6.4,a)') ch10,"== Thermodynamic integration over interaction for lambda= ",lam_list(ilam),ch10
       call wrtout(std_out,message,'COLL')
       if (integral == 2) then
         do isppol=1,nsppol
           levels_ctqmc(1+(isppol-1)*ndim:tndim+(isppol-1)*ndim,1+(isppol-1)*ndim:tndim+(isppol-1)*ndim) = &
             & energy_level%matlu(iatom)%mat(:,:,isppol) + (one-lam_list(ilam))*hdc_ctqmc%matlu(iatom)%mat(:,:,isppol)
           if (nsppol == 1 .and. nspinor == 1) levels_ctqmc(1+ndim:2*ndim,1+ndim:2*ndim) = levels_ctqmc(1:ndim,1:ndim)
         end do ! isppol
       end if ! integral=2
     end if ! ilam/=ntot

     if (ilam == ntot .and. entropy .and. integral == 2) then
       do isppol=1,nsppol
         levels_ctqmc(1+(isppol-1)*ndim:tndim+(isppol-1)*ndim,1+(isppol-1)*ndim:tndim+(isppol-1)*ndim) = &
            & energy_level%matlu(iatom)%mat(:,:,isppol)
         if (nsppol == 1 .and. nspinor == 1) levels_ctqmc(1+ndim:2*ndim,1+ndim:2*ndim) = levels_ctqmc(1:ndim,1:ndim)
       end do ! isppol
     end if

     if (ilam == 2) verbo = 0

     if (ilam < 10) then
       write(tag_lam,'("0",i1)') ilam
     else
       write(tag_lam,'(i2)') ilam
     end if

     ilam_prev = merge(ntot,merge(1,ilam-1,ilam==ntot),ilam==1) ! Index of the previous integration point

     if (ilam_prev < 10) then
       write(tag_lam_prev,'("0",i1)') ilam_prev
     else
       write(tag_lam_prev,'(i2)') ilam_prev
     end if

     tag_lam2 = ""
     if (ilam /= ntot) tag_lam2 = "_ilam" // tag_lam

     tag_lam_prev2 = ""
     if (ilam_prev /= ntot) tag_lam_prev2 = "_ilam" // tag_lam_prev

     read_data = paw_dmft%dmft_triqs_read_ctqmcdata
     stringfile = "_iatom" // tag_at // trim(adjustl(tag_lam2)) // ".h5"
     stringfile_prev = "_iatom" // tag_at // trim(adjustl(tag_lam_prev2)) // ".h5"

     if (paw_dmft%idmftloop == 1) then
       read_data = 0
       if (paw_dmft%dmft_triqs_read_ctqmcdata == 1 .and. paw_dmft%ireadctqmcdata == 1) read_data = 1
       fname_data = trim(adjustl(paw_dmft%filctqmcdatain)) // stringfile
       inquire(file=trim(fname_data),exist=lexist)
       if ((.not. lexist) .and. ntot > 1) then ! try to restart from config of previous integration point instead
         if (read_data == 1) read_data = 2 ! to indicate not to use the value qmc_therm for warmup
         if (ilam == 1) then
           fname_data = trim(adjustl(paw_dmft%filctqmcdatain)) // stringfile_prev
         else
           fname_data = trim(adjustl(paw_dmft%filapp)) // "_CTQMC_DATA" // stringfile_prev
         end if
       end if
     else
       fname_data = trim(adjustl(paw_dmft%filapp)) // "_CTQMC_DATA" // stringfile
     end if
     len_t = len(trim(adjustl(fname_data))) + 1
     fname_data(len_t:len_t) = c_null_char

     fname_dataw = trim(adjustl(paw_dmft%filapp)) // "_CTQMC_DATA" // stringfile
     len_t = len(trim(adjustl(fname_dataw))) + 1
     fname_dataw(len_t:len_t) = c_null_char

     fname_histo = trim(adjustl(paw_dmft%filapp)) // "_CTQMC_HISTOGRAM_iatom" // tag_at // trim(adjustl(tag_lam2)) // ".dat"
     len_t = len(trim(adjustl(fname_histo))) + 1
     fname_histo(len_t:len_t) = c_null_char

     call flush_unit(std_out)

#if defined HAVE_TRIQS_v4_0 || defined HAVE_TRIQS_v3_2
     call Ctqmc_triqs_run(rot_inv,leg_measure,paw_dmft%dmft_triqs_move_shift,paw_dmft%dmft_triqs_move_double, &
                        & density_matrix,paw_dmft%dmft_triqs_time_invariance,paw_dmft%dmft_triqs_use_norm_as_weight, &
                        & debug,merge(integral,0,ilam/=ntot),paw_dmft%dmft_triqs_loc_n_min,paw_dmft%dmft_triqs_loc_n_max, &
                        & paw_dmft%dmft_triqs_seed_a,paw_dmft%dmft_triqs_seed_b,nflavor,ntau,nleg, &
                        & paw_dmft%dmft_triqs_n_cycles,paw_dmft%dmftctqmc_meas,paw_dmft%dmftqmc_therm, &
                        & paw_dmft%dmft_triqs_therm_restart,paw_dmft%dmft_triqs_det_init_size, &
                        & paw_dmft%dmft_triqs_det_n_operations_before_check,myproc,nblocks(iatom),read_data,verbo, &
                        & beta,paw_dmft%dmft_triqs_imag_threshold,paw_dmft%dmft_triqs_det_precision_warning, &
                        & paw_dmft%dmft_triqs_det_precision_error,paw_dmft%dmft_triqs_det_singular_threshold,lam_list(ilam), &
                        & paw_dmft%dmft_triqs_pauli_prob,block_ptr,flavor_ptr,inner_ptr,siz_ptr,ftau_ptr,gtau_ptr,gl_ptr, &
                        & udens_ptr,vee_ptr,levels_ptr,mself_1_ptr,mself_2_ptr,occ_ptr,eu_ptr,fname_data_ptr,fname_dataw_ptr, fname_histo_ptr)
#endif

     call flush_unit(std_out)

     if (ilam == ntot .or. debug) then

       do isppol=1,nsppol
         if (nsppol == 1 .and. nspinor == 1) then
           green%oper_tau(1)%matlu(iatom)%mat(:,:,isppol) = (gtau(1,1:ndim,1:ndim)+gtau(1,ndim+1:2*ndim,ndim+1:2*ndim)) * half
         else
           green%oper_tau(1)%matlu(iatom)%mat(:,:,isppol) = gtau(1,1+(isppol-1)*ndim:tndim+(isppol-1)*ndim,1+(isppol-1)*ndim:tndim+(isppol-1)*ndim)
         end if
       end do ! isppol

       if (ilam < ntot) then
         call occup_green_tau(green)

         write(message,'(a,3x,a)') ch10,"== Print Occupation matrix in CTQMC basis"
         call wrtout(std_out,message,"COLL")
         call print_matlu(green%occup_tau%matlu(:),natom,1)
       end if

       if ((.not. leg_measure) .and. density_matrix) then

         ! Constrain the occupations and high-frequency moments with the more accurate values sampled from the CTQMC

         do isppol=1,nsppol
           do im=1,tndim
             iflavor = im + (isppol-1)*ndim

             if (nsppol == 1 .and. nspinor == 1) then
               mself_1 = (moments_self_1(iflavor)+moments_self_1(iflavor+ndim)) * half
               mself_2 = (moments_self_2(iflavor)+moments_self_2(iflavor+ndim)) * half
             else
               mself_1 = moments_self_1(iflavor)
               mself_2 = moments_self_2(iflavor)
             end if

             green%moments(2)%matlu(iatom)%mat(im,im,isppol) = energy_level%matlu(iatom)%mat(im,im,isppol) + mself_1
             green%moments(3)%matlu(iatom)%mat(im,im,isppol) = weiss%moments(2)%matlu(iatom)%mat(im,im,isppol) + mself_2

           end do ! im

           ! Use matmul in prevision of the day where the off-diagonal density matrix will be available
           green%moments(3)%matlu(iatom)%mat(:,:,isppol) = green%moments(3)%matlu(iatom)%mat(:,:,isppol) + &
             & matmul(green%moments(2)%matlu(iatom)%mat(:,:,isppol),green%moments(2)%matlu(iatom)%mat(:,:,isppol))

         end do ! isppol

       end if ! not leg and density_matrix

       if (density_matrix) then

         green%ecorr_qmc(iatom) = dble(eu)

         do isppol=1,nsppol
           do im=1,tndim
             iflavor = im + (isppol-1)*ndim

             if (nsppol == 1 .and. nspinor == 1) then
               occ_tmp = (occ(iflavor)+occ(iflavor+ndim)) * half
             else
               occ_tmp = occ(iflavor)
             end if

             green%oper_tau(1)%matlu(iatom)%mat(im,im,isppol) = occ_tmp - cone

           end do ! im
         end do ! isppol

       end if ! density_matrix

       if (leg_measure) then

         ABI_MALLOC(gl_tmp,(nleg,tndim,tndim,nsppol))
         ABI_MALLOC(jbes,(nleg))

         gl_tmp(:,:,:,:) = czero

         do ifreq=1,nwlo
           green%oper(ifreq)%matlu(iatom)%mat(:,:,:) = czero
           xx = dble(2*ifreq-1) * pi / two
           if (xx <= dble(100)) then
             call sbf8(nleg,xx,jbes(:))
           end if
           do isppol=1,nsppol
             do im1=1,tndim
               iflavor1 = im1 + (isppol-1)*ndim
               do im=1,tndim
                 iflavor = im + (isppol-1)*ndim
                 do ileg=1,nleg
                   if (xx >= dble(99)) then
                     call jbessel(jbes(ileg),besp,bespp,ileg-1,1,xx)
                   end if
                   u_nl = sqrt(dble(2*ileg-1))*(-1)**(ifreq-1)*(j_dpc**(ileg))*jbes(ileg)
                   if (nsppol == 1 .and. nspinor == 1) then
                     gl_tmp(ileg,im,im1,isppol) = &
                       & (gl(ileg,iflavor,iflavor1)+gl(ileg,iflavor+ndim,iflavor1+ndim))*half
                   else
                     gl_tmp(ileg,im,im1,isppol) = gl(ileg,iflavor,iflavor1)
                   end if
                   green%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol) = &
                     & green%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol) + &
                     & u_nl*gl_tmp(ileg,im,im1,isppol)
                 end do ! ileg
               end do ! im
             end do ! im1
           end do ! isppol
         end do ! ifreq

         ABI_FREE(jbes)

         ABI_MALLOC(gtau_leg,(ntau,nflavor,nflavor))
         ABI_MALLOC(leg_array,(nleg,ntau))

         do itau=1,ntau
           tau  = dble(itau-1) * beta / (ntau-1)
           xtau = two*tau/beta - one
           leg_array(1,itau) = one
           leg_array(2,itau) = xtau
           do ileg=3,nleg
             leg_array(ileg,itau) = (dble(2*(ileg-2)+1)*xtau*leg_array(ileg-1,itau)- &
             & dble(ileg-2)*leg_array(ileg-2,itau)) / dble(ileg-1)
           end do ! ileg
         end do ! itau

         gtau_leg(:,:,:) = czero
         do iflavor1=1,nflavor
           do iflavor=1,nflavor
             do itau=1,ntau
               do ileg=1,nleg
                 gtau_leg(itau,iflavor,iflavor1) = gtau_leg(itau,iflavor,iflavor1) + &
                  & gl(ileg,iflavor,iflavor1)*leg_array(ileg,itau)*sqrt(dble(2*ileg-1))/beta
               end do ! ileg
             end do ! itau
           end do ! iflavor
         end do ! iflavor1

         if (myproc == 0 .and. off_diag) then

           if (open_file(trim(paw_dmft%filapp)//"_Gtau_offdiag_Leg_iatom"//tag_at//trim(adjustl(tag_lam2))//".dat",message,newunit=unt) /= 0) ABI_ERROR(message)
           write(unt,'(6a)') "# Off-diagonal components of Legendre-sampled G(tau) in the CTQMC basis",ch10, &
                           & "# Columns are ordered this way:",ch10, &
                           & "# Imaginary Time     ((Re(G_{ij}) Im(G_{ij}),i=1,2*(2*l+1)),j=1,2*(2*l+1)) where the", &
                           & " leftmost index varies first"
           do itau=1,ntau
             write(unt,'(2x,393(e18.10e3,2x))') beta*dble(itau-1)/dble(ntau-1), &
               & ((dble(gtau_leg(itau,im,im1)),aimag(gtau_leg(itau,im,im1)),im=1,nflavor),im1=1,nflavor)
           end do ! itau
           close(unt)

         end if ! myproc=0

         if (myproc == 0) then

           if (open_file(trim(paw_dmft%filapp)//"_Gtau_diag_Leg_iatom"//tag_at//trim(adjustl(tag_lam2))//".dat",message,newunit=unt) /= 0) ABI_ERROR(message)
           write(unt,'(5a)') "# Diagonal components of Legendre-sampled G(tau) in the CTQMC basis",ch10, &
                           & "# Columns are ordered this way:",ch10, &
                           & "# Imaginary Time     (G_{ii},i=1,2*(2*l+1))"
           do itau=1,ntau
             write(unt,'(2x,393(e25.17e3,2x))') beta*dble(itau-1)/dble(ntau-1),(dble(gtau_leg(itau,im,im)),im=1,nflavor)
           end do ! itau
           close(unt)

         end if ! myproc=0

         ABI_FREE(gtau_leg)
         ABI_FREE(leg_array)

         ABI_MALLOC(t_lp,(nleg,green%nmoments))

         fact = one ! this is equal to (p-1)!

         ! Compute analytical moments of Fourier transform of Legendre polynomial
         ! (equation (E2) of PRB, 84(7), 2011, Boehnke et al)
         do p=1,green%nmoments
           if (p > 1) fact = fact * dble(p-1)
           fact2 = fact
           do l=0,nleg-1
             if (l > 0) fact2 = fact2 * dble(l+p-1)
             if (p > l+1) then
               t_lp(l+1,p) = zero
             else
               if (l-p+1 > 0) fact2 = fact2 / dble(l-p+1)  ! fact2 is now equal to (l+p-1)...(l-p+2)
               if (mod(p+l,2) == 0) then
                 t_lp(l+1,p) = zero
               else
                 t_lp(l+1,p) = (-1)**p * two * sqrt(dble(2*l+1)) * fact2 / fact
               end if
             end if ! p>l+1
           end do ! l
           do isppol=1,nsppol
             do im1=1,tndim
               do im=1,tndim
                 ! Do not use DOT_PRODUCT
                 green%moments(p)%matlu(iatom)%mat(im,im1,isppol) = sum(t_lp(:,p)*gl_tmp(:,im,im1,isppol)) / beta**p
               end do ! im
             end do ! im1
           end do ! isppol
         end do ! p

         ABI_FREE(gl_tmp)
         ABI_FREE(t_lp)

       else

         ABI_MALLOC(gl_dlr,(ndlr,tndim,tndim,nsppol))
         ABI_MALLOC(gl_dlr_re,(ndlr))
         ABI_MALLOC(gl_dlr_im,(ndlr))
         ABI_MALLOC(bdlr,(ntau))
         ABI_MALLOC(gtau_dlr,(ntau,nflavor,nflavor))
         ABI_MALLOC(moment_fit,(green%nmoments))

         gtau_dlr(:,:,:) = czero

         call fit_dlr()

         if (myproc == 0 .and. off_diag) then

           if (open_file(trim(paw_dmft%filapp)//"_Gtau_offdiag_DLR_iatom"//tag_at//trim(adjustl(tag_lam2))//".dat",message,newunit=unt) /= 0) ABI_ERROR(message)
           write(unt,'(6a)') "# Off-diagonal components of DLR fit of G(tau) in the CTQMC basis",ch10, &
                           & "# Columns are ordered this way:",ch10, &
                           & "# Imaginary Time     ((Re(G_{ij}) Im(G_{ij}),i=1,2*(2*l+1)),j=1,2*(2*l+1)) where the", &
                           & " leftmost index varies first"

           do itau=1,ntau
             write(unt,'(2x,393(e18.10e3,2x))') beta*dble(itau-1)/dble(ntau-1), &
                & ((dble(gtau_dlr(itau,im,im1)),aimag(gtau_dlr(itau,im,im1)),im=1,nflavor),im1=1,nflavor)
           end do ! itau
           close(unt)

         end if ! myproc=0 and off_diag

         if (myproc == 0) then

           if (open_file(trim(paw_dmft%filapp)//"_Gtau_diag_DLR_iatom"//tag_at//trim(adjustl(tag_lam2))//".dat",message,newunit=unt) /= 0) ABI_ERROR(message)
           write(unt,'(5a)') "# Diagonal components of DLR fit of G(tau) in the CTQMC basis",ch10, &
                           & "# Columns are ordered this way:",ch10, &
                           & "# Imaginary Time     (G_{ii},i=1,2*(2*l+1))"
           do itau=1,ntau
             write(unt,'(2x,393(e25.17e3,2x))') beta*dble(itau-1)/dble(ntau-1),(dble(gtau_dlr(itau,im,im)),im=1,nflavor)
           end do ! itau
           close(unt)

         end if ! myproc

         ABI_FREE(gl_dlr)
         ABI_FREE(gl_dlr_re)
         ABI_FREE(gl_dlr_im)
         ABI_FREE(bdlr)
         ABI_FREE(gtau_dlr)
         ABI_FREE(moment_fit)

       end if ! leg_measure

       if (myproc == 0 .and. off_diag) then

         if (open_file(trim(paw_dmft%filapp)//"_Gtau_offdiag_iatom"//tag_at//trim(adjustl(tag_lam2))//".dat",message,newunit=unt) /= 0) ABI_ERROR(message)
         write(unt,'(6a)') "# Off-diagonal components of binned G(tau) in the CTQMC basis",ch10, &
                         & "# Columns are ordered this way:",ch10, &
                         & "# Imaginary Time     ((Re(G_{ij}) Im(G_{ij}),i=1,2*(2*l+1)),j=1,2*(2*l+1)) where the", &
                         & " leftmost index varies first"

         do itau=1,ntau
           write(unt,'(2x,393(e18.10e3,2x))') beta*dble(itau-1)/dble(ntau-1), &
                & ((dble(gtau(itau,im,im1)),aimag(gtau(itau,im,im1)),im=1,nflavor),im1=1,nflavor)
         end do ! itau
         close(unt)

       end if ! myproc

       if (myproc == 0) then

         if (open_file(trim(paw_dmft%filapp)//"_Gtau_diag_iatom"//tag_at//trim(adjustl(tag_lam2))//".dat",message,newunit=unt) /= 0) ABI_ERROR(message)
         write(unt,'(5a)') "# Diagonal components of binned G(tau) in the CTQMC basis",ch10, &
                         & "# Columns are ordered this way:",ch10, &
                         & "# Imaginary Time     (G_{ii},i=1,2*(2*l+1))"
         do itau=1,ntau
           write(unt,'(2x,393(e25.17e3,2x))') beta*dble(itau-1)/dble(ntau-1),(dble(gtau(itau,im,im)),im=1,nflavor)
         end do ! itau
         close(unt)

       end if ! myproc

       if (entropy .and. (.not. debug)) then

         call compute_migdal_energy(emig(:),emig_tot,green,paw_dmft,hybmwdhyb,iatom=iatom)
         green%ekin_imp = green%ekin_imp + two*emig_tot

       end if ! entropy

     end if ! ilam=ntot or debug

     if (integral > 0 .and. ilam < ntot .and. entropy) then

       elam = dble(eu)

       if (integral == 2) then
         ! CAREFUL: if one day the density matrix sampling is implemented for
         ! off-diagonal components, you would need to add them here
         do isppol=1,nsppol
           do im=1,tndim
             iflavor = im + (isppol-1)*ndim

             if (nsppol == 1 .and. nspinor == 1) then
               occ_tmp = occ(iflavor) + occ(iflavor+ndim)
             else
               occ_tmp = occ(iflavor)
             end if
             elam = elam - dble(hdc_ctqmc%matlu(iatom)%mat(im,im,isppol)*occ_tmp)
           end do ! im
         end do ! isppol
       end if ! integral=2

       i = mod(ilam-1,ngauss) + 1
       green%integral = green%integral + tweights(i)*elam*dx*half
       elam_list(ilam) = elam

     end if ! integral and ilam<ntot and entropy

     if (integral > 0 .and. ilam == ntot-1 .and. entropy) then
       write(message,'(a,3(3x,2a),a,12x,a,6x,2a,8x,a)') ch10,repeat("=",39),ch10,"== Summary of thermodynamic integration", &
            & ch10,repeat("=",39),ch10,ch10,"Lambda","<dH/dlambda>",ch10,repeat("-",29)
       call wrtout(std_out,message,'COLL')
       do i=ntot-1,1,-1
         write(tag_lambda,'(f14.4)') lam_list(i)
         write(tag_elam,'(f14.4)') elam_list(i)
         tag_lambda = adjustl(tag_lambda)
         tag_elam = adjustl(tag_elam)
         pad_lambda = (14-len_trim(tag_lambda)) / 2
         pad_elam = (14-len_trim(tag_elam)) / 2
         write(message,'(8x,2(3a,1x),a,8x,a)') repeat(" ",pad_lambda),trim(tag_lambda),repeat(" ",14-pad_lambda-len_trim(tag_lambda)), &
                                             & repeat(" ",pad_elam),trim(tag_elam),repeat(" ",14-pad_elam-len_trim(tag_elam)),ch10,repeat("-",29)
         call wrtout(std_out,message,'COLL')
       end do ! i
       write(message,'(a,3x,a,f10.4,a)') ch10,"--> Integral is: ",green%integral,ch10
       call wrtout(std_out,message,'COLL')
     end if ! integral and ilam=ntot-1 and entropy

   end do ! ilam

   ABI_FREE(flavor_tmp)
   ABI_FREE(gtau)
   ABI_FREE(levels_ctqmc)
   ABI_SFREE(gl)
   ABI_SFREE(moments_self_1)
   ABI_SFREE(moments_self_2)
   ABI_SFREE(occ)

 end do ! iatom

 ABI_FREE(elam_list)
 ABI_FREE(lam_list)

 ABI_SFREE(tweights)
 ABI_SFREE(tpoints)

 call occup_green_tau(green)

 write(message,'(a,3x,a)') ch10,"== Print Occupation matrix in CTQMC basis"
 call wrtout(std_out,message,"COLL")
 call print_matlu(green%occup_tau%matlu(:),natom,1)

 if (basis > 0) then
   write(message,'(a,3x,a)') ch10,"== Rotating back to cubic basis"
   call wrtout(std_out,message,"COLL")
 end if

 ! Build back Weiss field
 do ifreq=1,nwlo
   shift(:) = cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)
   call add_matlu(weiss%oper(ifreq)%matlu(:),energy_level%matlu(:),matlu_tmp(:),natom,1)
   call copy_matlu(matlu_tmp(:),weiss%oper(ifreq)%matlu(:),natom)
   call shift_matlu(weiss%oper(ifreq)%matlu(:),natom,shift(:))
   call fac_matlu(weiss%oper(ifreq)%matlu(:),natom,-cone)
 end do ! ifreq

 if (basis > 0) then
   if (basis == 1 .or. basis == 2) then
     call rotate_matlu(energy_level%matlu(:),eigvectmatlu(:),natom,0)
     if (entropy .and. integral == 2) then
       call rotate_matlu(hdc_ctqmc%matlu(:),eigvectmatlu(:),natom,0)
     end if
   else if (basis == 4) then
     call ylm2jmj_matlu(energy_level%matlu(:),natom,2,paw_dmft)
     if (entropy .and. integral == 2) then
       call ylm2jmj_matlu(hdc_ctqmc%matlu(:),natom,2,paw_dmft)
     end if
   end if ! basis /= 3
   call slm2ylm_matlu(energy_level%matlu(:),natom,paw_dmft,2,0)
   if (entropy .and. integral == 2) then
     call slm2ylm_matlu(hdc_ctqmc%matlu(:),natom,paw_dmft,2,0)
   end if
   do i=2,weiss%nmoments-1
     if (basis == 1 .or. basis == 2) then
       call rotate_matlu(weiss%moments(i)%matlu(:),eigvectmatlu(:),natom,0)
     else if (basis == 4) then
       call ylm2jmj_matlu(weiss%moments(i)%matlu(:),natom,2,paw_dmft)
     end if
     call slm2ylm_matlu(weiss%moments(i)%matlu(:),natom,paw_dmft,2,0)
   end do ! i
   do i=1,green%nmoments
     if (basis == 1 .or. basis == 2) then
       call rotate_matlu(green%moments(i)%matlu(:),eigvectmatlu(:),natom,0)
     else if (basis == 4) then
       call ylm2jmj_matlu(green%moments(i)%matlu(:),natom,2,paw_dmft)
     end if
     call slm2ylm_matlu(green%moments(i)%matlu(:),natom,paw_dmft,2,0)
   end do ! i
   do ifreq=1,nwlo
     if (green%distrib%procf(ifreq) /= myproc) cycle
     if (basis == 1 .or. basis == 2) then
       call rotate_matlu(weiss%oper(ifreq)%matlu(:),eigvectmatlu(:),natom,0)
       call rotate_matlu(green%oper(ifreq)%matlu(:),eigvectmatlu(:),natom,0)
     else if (basis == 4) then
       call ylm2jmj_matlu(weiss%oper(ifreq)%matlu(:),natom,2,paw_dmft)
       call ylm2jmj_matlu(green%oper(ifreq)%matlu(:),natom,2,paw_dmft)
     end if
     call slm2ylm_matlu(weiss%oper(ifreq)%matlu(:),natom,paw_dmft,2,0)
     call slm2ylm_matlu(green%oper(ifreq)%matlu(:),natom,paw_dmft,2,0)
   end do ! ifreq
   if (basis == 1 .or. basis == 2) then
     call rotate_matlu(green%occup_tau%matlu(:),eigvectmatlu(:),natom,0)
   else if (basis == 4) then
     call ylm2jmj_matlu(green%occup_tau%matlu(:),natom,2,paw_dmft)
   end if
   call slm2ylm_matlu(green%occup_tau%matlu(:),natom,paw_dmft,2,0)
 end if! basis > 0

 ! Since we possibly neglected some off-diagonal elements and imaginary part,
 ! the levels and hybridization might not be symmetrized anymore
 call sym_matlu(energy_level%matlu(:),paw_dmft)
 if (entropy .and. integral == 2) then
   call sym_matlu(hdc_ctqmc%matlu(:),paw_dmft)
 end if
 do i=2,weiss%nmoments-1
   call sym_matlu(weiss%moments(i)%matlu(:),paw_dmft)
 end do ! i
 do i=1,green%nmoments
   call sym_matlu(green%moments(i)%matlu(:),paw_dmft)
 end do ! i
 do ifreq=1,nwlo
   if (green%distrib%procf(ifreq) /= myproc) cycle
   call sym_matlu(weiss%oper(ifreq)%matlu(:),paw_dmft)
   call sym_matlu(green%oper(ifreq)%matlu(:),paw_dmft)
 end do ! ifreq
 call copy_matlu(green%occup_tau%matlu(:),dmat_ctqmc(:),natom)
 call sym_matlu(green%occup_tau%matlu(:),paw_dmft)

 call diff_matlu("CTQMC occupations","Symmetrized CTQMC occupations",dmat_ctqmc(:),green%occup_tau%matlu(:),natom,0,tol4)

 call gather_oper(weiss%oper(:),weiss%distrib,paw_dmft,opt_ksloc=2)
 call gather_oper(green%oper(:),green%distrib,paw_dmft,opt_ksloc=2)

 call compute_moments_loc(green,self_new,energy_level,weiss,1,opt_log=merge(max(1,integral),0,entropy),opt_hdc=hdc_ctqmc)

 if (entropy .and. integral > 0) then
   call compute_trace_log_loc(weiss,paw_dmft,green%fband_weiss,opt_inv=max(1,integral),opt_hdc=hdc_ctqmc)
 end if

 call destroy_matlu(dmat_ctqmc(:),natom)
 call destroy_matlu(eigvectmatlu(:),natom)
 call destroy_matlu(ftau(:),natom)
 call destroy_matlu(matlu_tmp(:),natom)
 call destroy_matlu(udens_rot(:),natom)

 if (entropy .and. integral == 2) then
   call destroy_oper(hdc_ctqmc)
 end if

 call destroy_oper(energy_level)

 call destroy_vee(paw_dmft,vee_rot(:))

 ABI_FREE(block_list)
 ABI_FREE(dmat_ctqmc)
 ABI_FREE(eigvectmatlu)
 ABI_FREE(flavor_list)
 ABI_FREE(ftau)
 ABI_FREE(inner_list)
 ABI_FREE(matlu_tmp)
 ABI_FREE(nblocks)
 ABI_FREE(shift)
 ABI_FREE(siz_block)
 ABI_FREE(udens_rot)
 ABI_FREE(vee_rot)

 ABI_SFREE(adlr)
 ABI_SFREE(adlr_iw)
 ABI_SFREE(emig)
 ABI_SFREE(wdlr)
 ABI_SFREE(wdlr_beta)

 if (entropy) then
   call destroy_self(hybmwdhyb)
 end if

contains

subroutine cubic_spline()

!Arguments ------------------------------------
!Local variables ------------------------------
 complex(dp), allocatable :: y(:),yp(:)
! ************************************************************************

 ABI_MALLOC(y,(nwlo))
 ABI_MALLOC(yp,(nwlo))

 do iatom=1,natom
   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   tndim = nspinor * (2*lpawu+1)
   do isppol=1,nsppol
     do im1=1,tndim
       do im=1,tndim
         do ifreq=1,nwlo
           y(ifreq) = weiss%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol)
         end do ! ifreq
         ! Boundary condition not-a-knot seems to get the best results
         call spline2_complex(paw_dmft%omega_lo(:),y(:),nwlo,yp(:),czero,czero,3,3)
         do ifreq=1,nwlo
           hybmwdhyb%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol) = &
            & energy_level%matlu(iatom)%mat(im,im1,isppol) + y(ifreq) - &
            & paw_dmft%omega_lo(ifreq)*yp(ifreq)
         end do ! ifreq
       end do ! im
     end do ! im1
   end do ! isppol
 end do ! iatom

 ABI_FREE(y)
 ABI_FREE(yp)

end subroutine cubic_spline

subroutine fit_dlr()

!Arguments ------------------------------------
!Local variables ------------------------------
! ************************************************************************

 do isppol=1,nsppol
   do im1=1,tndim
     iflavor1 = im1 + (isppol-1)*ndim
     do im=1,tndim
       iflavor = im + (isppol-1)*ndim

       ncon = merge(merge(4,3,iflavor==iflavor1),1,density_matrix)
       nmoments = merge(3,1,density_matrix)
       if (nsppol == 1 .and. nspinor == 1) gtau(:,iflavor,iflavor1) = &
             & (gtau(:,iflavor,iflavor1)+gtau(:,iflavor+ndim,iflavor1+ndim)) * half
       occ_tmp = green%oper_tau(1)%matlu(iatom)%mat(im,im1,isppol)

       bdlr(:) = dble(gtau(:,iflavor,iflavor1))

       do i=1,nmoments
         moment_fit(i) = dble(green%moments(i)%matlu(iatom)%mat(im,im1,isppol))
       end do ! i

       call slsqp_wrapper(ncon,ndlr,lsq_g,con_moments,jac_lsq_g,jac_con_moments,gl_dlr_re(:))

       if (density_matrix) ncon = 3

       bdlr(:) = aimag(gtau(:,iflavor,iflavor1))
       do i=1,nmoments
         moment_fit(i) = aimag(green%moments(i)%matlu(iatom)%mat(im,im1,isppol))
       end do ! i

       call slsqp_wrapper(ncon,ndlr,lsq_g,con_moments,jac_lsq_g,jac_con_moments,gl_dlr_im(:))

       gl_dlr(:,im,im1,isppol) = cmplx(gl_dlr_re(:),gl_dlr_im(:),kind=dp)

       ! Do not use DOT_PRODUCT
       green%moments(1)%matlu(iatom)%mat(im,im1,isppol) = sum(gl_dlr(:,im,im1,isppol))
       do i=2,green%nmoments
         green%moments(i)%matlu(iatom)%mat(im,im1,isppol) = sum(gl_dlr(:,im,im1,isppol)*wdlr_beta(:,i-1))
       end do ! i

       do itau=1,ntau
         ! Do not use DOT_PRODUCT
         gtau_dlr(itau,iflavor,iflavor1) = sum(gl_dlr(:,im,im1,isppol)*adlr(:,itau))
       end do ! itau
       if (nsppol == 1 .and. nspinor == 1) gtau_dlr(:,iflavor+ndim,iflavor1+ndim) = gtau_dlr(:,iflavor,iflavor1)

     end do ! im
   end do ! im1
 end do ! isppol

 do ifreq=1,nwlo
   do isppol=1,nsppol
     do im1=1,tndim
       do im=1,tndim
         ! Do not use DOT_PRODUCT
         green%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol) = sum(gl_dlr(:,im,im1,isppol)*adlr_iw(:,ifreq))
       end do ! im
     end do ! im1
   end do ! isppol
 end do ! ifreq

end subroutine fit_dlr

subroutine diag_block(matlu)

!Arguments ------------------------------------
 type(matlu_type), intent(inout) :: matlu(natom)
!Local variables ------------------------------
 integer :: i,iatom,iblock,iflavor,iflavor1,im,im1
 integer :: info,is,j,lpawu,lwork,ndim,sizb,tndim
 real(dp), allocatable :: eig(:),rwork(:)
 complex(dp), allocatable :: mat_tmp(:,:),work(:)
! ************************************************************************

 is = 1
 ABI_MALLOC(eig,(nflavor_max))
 ABI_MALLOC(mat_tmp,(nflavor_max,nflavor_max))
 ABI_MALLOC(rwork,(3*nflavor_max-2))
 ABI_MALLOC(work,(2*nflavor_max-1))

 do iatom=1,natom
   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   ndim  = 2*lpawu + 1
   tndim = nspinor * ndim
   do iblock=1,nblocks(iatom)
     sizb = siz_block(iblock,iatom)
     lwork = 2*sizb - 1
     do j=1,sizb
       iflavor1 = flavor_list(j,iblock,iatom) + 1
       im1 = mod(iflavor1-1,tndim) + 1
       if (nspinor == 1) then
         is = (iflavor1-1)/ndim + 1
         if (is > nsppol) exit
       end if
       do i=1,sizb
         iflavor = flavor_list(i,iblock,iatom) + 1
         im = mod(iflavor-1,tndim) + 1
         mat_tmp(i,j) = matlu(iatom)%mat(im,im1,is)
         matlu(iatom)%mat(im,im1,is) = czero
       end do ! i
     end do ! j

     if (is > nsppol) cycle

     call zheev('v','u',sizb,mat_tmp(:,1:sizb),nflavor_max,eig(1:sizb), &
              & work(1:lwork),lwork,rwork(1:3*sizb-2),info)

     do j=1,sizb
       iflavor1 = flavor_list(j,iblock,iatom) + 1
       im1 = mod(iflavor1-1,tndim) + 1
       matlu(iatom)%mat(im1,im1,is) = cmplx(eig(j),zero,kind=dp)
       do i=1,sizb
         iflavor = flavor_list(i,iblock,iatom) + 1
         im = mod(iflavor-1,tndim) + 1
         eigvectmatlu(iatom)%mat(im,im1,is) = mat_tmp(i,j)
       end do ! i
     end do ! j
   end do ! iblock

   if (pawprtvol >= 3) then
     write(tag_at,'(i4)') iatom
     do isppol=1,nsppol
       write(message,'(4a,i1)') ch10,"       EIGENVECTORS for atom ",trim(adjustl(tag_at))
       if (nspinor == 1) then
         write(tag_block4,'(i1)') isppol
         message = trim(message) // " and isppol " // tag_block4
       end if
       call wrtout(std_out,message,'COLL')
       do im=1,tndim
         write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))') &
            & (eigvectmatlu(iatom)%mat(im,im1,isppol),im1=1,tndim)
         call wrtout(std_out,message,'COLL')
       end do ! im1
     end do ! isppol
   end if ! pawprtvol>=3

 end do ! iatom

 ABI_FREE(eig)
 ABI_FREE(mat_tmp)
 ABI_FREE(rwork)
 ABI_FREE(work)

end subroutine diag_block

subroutine lsq_g(gl,err)

!Arguments ------------------------------------
 real(dp), intent(in) :: gl(:)
 real(dp), intent(out) :: err
!Local variables ------------------------------
! ************************************************************************

 err = zero
 do itau=1,ntau
   err = err + (dot_product(gl(:),adlr(:,itau))-bdlr(itau))**2
 end do

end subroutine lsq_g

subroutine jac_lsq_g(gl,jac)

!Arguments ------------------------------------
 real(dp), intent(in) :: gl(:)
 real(dp), intent(inout) :: jac(:)
!Local variables ------------------------------
! ************************************************************************

 jac(:) = zero
 do itau=1,ntau
   jac(1:ndlr) = jac(1:ndlr) + (dot_product(gl(:),adlr(:,itau))-bdlr(itau))*adlr(:,itau)
 end do
 jac = jac * two

end subroutine jac_lsq_g

subroutine con_moments(gl,con)

!Arguments ------------------------------------
 real(dp), intent(in) :: gl(:)
 real(dp), intent(inout) :: con(:)
!Local variables ------------------------------
! ************************************************************************

 con(1) = sum(gl(:)) - moment_fit(1)
 if (ncon > 1) then
   con(2) = dot_product(gl(:),wdlr_beta(:,1)) - moment_fit(2)
   con(3) = dot_product(gl(:),wdlr_beta(:,2)) - moment_fit(3)
   if (ncon == 4) con(4) = dot_product(gl(:),adlr(:,1)) - dble(occ_tmp)
 end if ! ncon>1

end subroutine con_moments

subroutine jac_con_moments(gl,jac_con)

!Arguments ------------------------------------
 real(dp), intent(in) :: gl(:)
 real(dp), intent(inout) :: jac_con(:,:)
!Local variables ------------------------------
! ************************************************************************

 ABI_UNUSED(gl(:))

 jac_con(1,1:ndlr) = one
 if (ncon > 1) then
   jac_con(2,1:ndlr) = wdlr_beta(:,1)
   jac_con(3,1:ndlr) = wdlr_beta(:,2)
   if (ncon == 4) jac_con(4,1:ndlr) = adlr(:,1)
 end if ! ncon>1
 jac_con(:,ndlr+1) = zero  ! not sure if this is necessary, but this is what they do in the SCIPY interface with SLSQP

end subroutine jac_con_moments

end subroutine ctqmc_calltriqs_c
!!***

!!****f* m_forctqmc/k_it
!! NAME
!! k_it
!!
!! FUNCTION
!! Computes the imaginary time kernel K(tau,omega).
!!
!! INPUTS
!! tau = imaginary time divided by beta
!! omega = real frequency multiplied by beta
!!
!! OUTPUT
!!
!! SOURCE

function k_it(tau,omega)

!Arguments ------------------------------------
 real(dp), intent(in) :: tau, omega
 real(dp) :: k_it
! *********************************************************************

 k_it = merge(-exp(-tau*omega)/(one+exp(-omega)),-exp((one-tau)*omega)/(one+exp(omega)),omega>=0)

end function k_it
!!***

!!****f* m_forctqmc/k_iw
!! NAME
!! k_iw
!!
!! FUNCTION
!! Computes the imaginary frequency kernel K(iom,omega)
!!
!! INPUTS
!! iom = imaginary part of the Matsubara frequency
!! omega = real frequency
!!
!! OUTPUT
!!
!! SOURCE

function k_iw(iom,omega)

!Arguments ------------------------------------
 real(dp), intent(in) :: iom,omega
 complex(dp) :: k_iw
! *********************************************************************

 k_iw = cone / (cmplx(zero,iom,kind=dp)-omega)

end function k_iw
!!***

!!****f* m_forctqmc/slsqp_wrapper
!! NAME
!! slsqp_wrapper
!!
!! FUNCTION
!! Optimizes a function with several variables
!! under several constraints, using the SLSQP algorithm.
!!
!! INPUTS
!! m = number of constraints
!! n = number of variables
!! fun = function to optimize
!! con = constraints
!! jac = jacobian of the function
!! jac_con = jacobian of the constraints
!!
!! OUTPUT
!! x(n) = minimizer
!!
!! SOURCE

subroutine slsqp_wrapper(m,n,fun,con,jac,jac_con,x)

 use m_slsqp, only : slsqp

!Arguments ------------------------------------
 integer, intent(in) :: m,n
 real(dp), intent(inout) :: x(n)

 interface

   subroutine fun(x,f)
     use defs_basis
     real(dp), intent(in) :: x(:)
     real(dp), intent(out) :: f
   end subroutine fun

   subroutine con(x,c)
     use defs_basis
     real(dp), intent(in) :: x(:)
     real(dp), intent(inout) :: c(:)
   end subroutine con

   subroutine jac(x,g)
     use defs_basis
     real(dp), intent(in) :: x(:)
     real(dp), intent(inout) :: g(:)
   end subroutine jac

   subroutine jac_con(x,a)
     use defs_basis
     real(dp), intent(in) :: x(:)
     real(dp), intent(inout) :: a(:,:)
   end subroutine jac_con

 end interface
!Local variables ------------------------------
 integer :: i,iexact,incons,ireset,iter,itermx,l_jw,l_w,la
 integer :: line,maxiter,meq,mineq,mode,n1,n2,n3
 real(dp) :: f,acc,alpha,f0,gs,h1,h2,h3,h4,t,t0,tol
 real(dp) :: c(max(m,1)),g(n+1),xl(n),xu(n)
 real(dp) :: a(max(m,1),n+1)
 real(dp), allocatable :: w(:)
 integer, allocatable :: jw(:)
! ************************************************************************

 maxiter = 10000
 meq = m   ! all constraints are equality constraints here
 la  = max(m,1)
 x(:)  = zero ; xl(:) = zero ; xu(:) = zero
 xl(:) = xl(:) / zero  ! set lower and upper bounds to NaN (very important, this is how slsqp recognizes that no bounds should be applied)
 xu(:) = xu(:) / zero

 call fun(x(:),f)
 call con(x(:),c(:))
 call jac(x(:),g(:))
 call jac_con(x(:),a(:,:))
 acc = tol8  ! best not to overconverge the result, in order to avoid overfitting
 iter = maxiter
 mode = 0
 n1 = n + 1
 mineq = m - meq + 2*n1
 l_w = (3*n1+m)*(n1+1) + (n1-meq+1)*(mineq+2) + 2*mineq &   ! as recommended
      & +(n1+mineq)*(n1-meq) + 2*meq + n1 + n1*n/2 + 2*m + 3*n + 3*n1 + 1
 ABI_MALLOC(w,(l_w))
 w = zero
 l_jw = mineq   ! as recommended
 ABI_MALLOC(jw,(l_jw))
 jw = 0 ; alpha = zero ; f0 = zero
 gs = zero ; h1 = zero ; h2 = zero
 h3 = zero ; h4 = zero ; t = zero
 t0 = zero ; tol = zero ; iexact = 0
 incons = 0 ; ireset = 0 ; itermx = 0
 line = 0 ; n1 = 0 ; n2 = 0 ; n3 = 0

 do i=1,maxiter

   call slsqp(m,meq,la,n,x(:),xl(:),xu(:),f,c(:),g(:),a(:,:),acc,iter,mode,w(:),l_w, &
            & jw(:),l_jw,alpha,f0,gs,h1,h2,h3,h4,t,t0,tol,iexact,incons,ireset,itermx, &
            & line,n1,n2,n3)

   if (abs(mode) /= 1) exit
   if (mode == -1) then
     call jac(x(:),g(:))
     call jac_con(x(:),a(:,:))
   end if
   if (mode == 1) then
     call fun(x(:),f)
     call con(x(:),c(:))
   end if
 end do ! i

 if (mode /= 0) ABI_ERROR("Error in the optimization procedure during the DLR fit")

 ABI_FREE(w)
 ABI_FREE(jw)

end subroutine slsqp_wrapper
!!***

!!****f* m_forctqmc/fourier_inv
!! NAME
!! fourier_inv
!!
!! FUNCTION
!!  Computes the inverse Fourier transform of a frequency-dependent operator,
!!  using analytical formulas for the asymptotic behavior. It is assumed the
!!  asymptotic behavior is moments(1)/(iw_n) + moments(2)/(iwn)**2 + ...,
!!  so if you want to use this routine for an operator with a 0th order moment,
!!  you need to subtract it before calling this routine.
!!
!! INPUTS
!!  paw_dmft <type(paw_dmft_type)>= DMFT data structure
!!  nmoments = number of moments
!!  ntau = number of (equidistant) imaginary time points on [0,beta]
!!  oper_freq = operator for each Matsubara frequency
!!  moments = high-frequency moments of the operator
!!
!! OUTPUT
!!  matlu_tau(2*(2*lpawu+1),2*(2*lpawu+1),ntau) = operator for each tau point
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

subroutine fourier_inv(paw_dmft,nmoments,ntau,matlu_tau,oper_freq,moments)

!Arguments ------------------------------------
 integer, intent(in) :: nmoments,ntau
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(oper_type), intent(in) :: oper_freq(paw_dmft%dmft_nwlo),moments(nmoments)
 type(matlu_type), intent(inout) :: matlu_tau(paw_dmft%natom)
!Local variables ------------------------------
 integer :: i,iatom,ibuf,ibuf_tau,ierr,ifreq,im1,isppol
 integer :: itau,itaub,itauf,lpawu,myproc,natom,ndim,nproc
 integer :: nspinor,nsppol,ntau_proc,nwlo,ratio,residu,siz_buf,tndim
 real(dp) :: beta,omegatau,tau
 complex(dp) :: fac
 integer, allocatable :: displs(:),recvcounts(:)
 complex(dp), allocatable :: buffer(:),buffer_tot(:),omega_fac(:)
! ************************************************************************

 beta    = one / paw_dmft%temp
 myproc  = paw_dmft%myproc
 natom   = paw_dmft%natom
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol
 nproc   = paw_dmft%nproc
 nwlo    = paw_dmft%dmft_nwlo

 call zero_matlu(matlu_tau(:),natom)

 ABI_MALLOC(displs,(nproc))
 ABI_MALLOC(recvcounts,(nproc))

 ratio  = ntau / nproc
 residu = ntau - ratio*nproc

 itau = 1
 do i=0,nproc-1
   ntau_proc = merge(ratio+1,ratio,i<residu)
   recvcounts(i+1) = ntau_proc
   if (myproc == i) itaub = itau
   itau = itau + ntau_proc
 end do ! i
 itauf = itaub + recvcounts(myproc+1) - 1

 siz_buf = 0
 do iatom=1,natom
   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   siz_buf = siz_buf + (2*lpawu+1)**2
 end do ! iatom

 siz_buf = siz_buf * (nspinor**2) * nsppol

 recvcounts(:) = recvcounts(:) * siz_buf

 displs(1) = 0
 do i=2,nproc
   displs(i) = displs(i-1) + recvcounts(i-1)
 end do ! i

 ABI_MALLOC(buffer,(recvcounts(myproc+1)))
 ABI_MALLOC(buffer_tot,(recvcounts(nproc)+displs(nproc)))
 ABI_MALLOC(omega_fac,(nmoments))

 buffer(:) = czero

 ibuf_tau = 0
 do itau=itaub,itauf

   tau = dble(itau-1) * beta / dble(ntau-1)
   omega_fac(:) = czero

   do ifreq=nwlo,1,-1 ! NEVER change this summation order and DON'T replace by the intrinsic SUM
     omegatau = mod(paw_dmft%omega_lo(ifreq)*tau,two_pi)
     fac = two * paw_dmft%temp * exp(-j_dpc*omegatau)
     do i=1,nmoments
       omega_fac(i) = omega_fac(i) - fac/(j_dpc*paw_dmft%omega_lo(ifreq))**i
     end do
     ibuf = 0
     do iatom=1,natom
       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       tndim = nspinor * (2*lpawu+1)
       do isppol=1,nsppol
         do im1=1,tndim
           buffer(ibuf_tau+ibuf+1:ibuf_tau+ibuf+tndim) = buffer(ibuf_tau+ibuf+1:ibuf_tau+ibuf+tndim) + &
                   fac*oper_freq(ifreq)%matlu(iatom)%mat(:,im1,isppol)
           ibuf = ibuf + tndim
         end do ! im
       end do ! isppol
     end do ! iatom
   end do ! ifreq

   omega_fac(1) = omega_fac(1) - half
   omega_fac(2) = omega_fac(2) + tau/two - beta/four
   omega_fac(3) = omega_fac(3) - (tau**2)/four + tau*beta/four

   do i=1,nmoments
     ibuf = 0
     do iatom=1,natom
       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       tndim = nspinor * (2*lpawu+1)
       do isppol=1,nsppol
         do im1=1,tndim
           buffer(ibuf_tau+ibuf+1:ibuf_tau+ibuf+tndim) = buffer(ibuf_tau+ibuf+1:ibuf_tau+ibuf+tndim) + &
                    & moments(i)%matlu(iatom)%mat(:,im1,isppol)*omega_fac(i)
           ibuf = ibuf + tndim
         end do ! im1
       end do ! isppol
     end do ! iatom
   end do ! i

   ibuf_tau = ibuf_tau + siz_buf

 end do ! itau

 ABI_FREE(omega_fac)

 call xmpi_allgatherv(buffer(:),recvcounts(myproc+1),buffer_tot(:),recvcounts(:),displs(:),paw_dmft%spacecomm,ierr)

 ABI_FREE(displs)
 ABI_FREE(recvcounts)

 ibuf = 0
 do itau=1,ntau
   do iatom=1,natom
     lpawu = paw_dmft%lpawu(iatom)
     if (lpawu == -1) cycle
     ndim  = 2*lpawu + 1
     tndim = ndim * nspinor
     do isppol=1,nsppol
       do im1=1,tndim
         matlu_tau(iatom)%mat(1+(isppol-1)*ndim:tndim+(isppol-1)*ndim,im1+(isppol-1)*ndim,itau) = buffer_tot(ibuf+1:ibuf+tndim)
         ibuf = ibuf + tndim
       end do ! im1
     end do ! isppol
     if (nsppol == 1 .and. nspinor == 1) matlu_tau(iatom)%mat(ndim+1:2*ndim,ndim+1:2*ndim,itau) = &
         & matlu_tau(iatom)%mat(1:ndim,1:ndim,itau)
     ndim = 2 * ndim
     matlu_tau(iatom)%mat(1:ndim,1:ndim,itau) = (matlu_tau(iatom)%mat(1:ndim,1:ndim,itau)+ &
         & conjg(transpose(matlu_tau(iatom)%mat(1:ndim,1:ndim,itau)))) * half
   end do ! iatom
 end do ! itau

 ABI_FREE(buffer)
 ABI_FREE(buffer_tot)

end subroutine fourier_inv
!!***

!!****f* m_forctqmc/find_block_structure
!! NAME
!! find_block_structure
!!
!! FUNCTION
!!  Find the most optimal block structure of a matlu
!!
!! INPUTS
!!  paw_dmft <type(paw_dmft_type)>= DMFT data structure
!!  matlu <type(oper_type)>= matrix for which the block structure is to be found
!!  natom = number of atoms
!!  nflavor_max = max number of orbitals
!!  hyb <type(green_type)>= hybridization ; if present, the block structure will
!!                          match both matlu and hyb
!!
!! OUTPUT
!!  block_list(nflavor,natom) = block index for each flavor and atom
!!  inner_list(nflavor,natom) = inner block index for each flavor and atom
!!  flavor_list(nflavor,nflavor,natom) = flavor for each block and inner indexes
!!  siz_block(nflavor,natom) = block size for each block and atom
!!  nblocks = number of blocks for each atom
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

subroutine find_block_structure(paw_dmft,block_list,inner_list,flavor_list, &
                              & siz_block,nblocks,matlu,natom,nflavor_max,hyb)

!Arguments ------------------------------------
 integer, intent(in) :: natom,nflavor_max
 integer, intent(inout) :: block_list(nflavor_max,natom),inner_list(nflavor_max,natom)
 integer, intent(inout) :: flavor_list(nflavor_max,nflavor_max,natom)
 integer, intent(inout) :: siz_block(nflavor_max,natom),nblocks(natom)
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(green_type), optional, intent(inout) :: hyb
 type(matlu_type), intent(inout) :: matlu(natom)
!Local variables ------------------------------
 integer :: i,iatom,iblock,iblock1,iblock2,iflavor,ifreq,lpawu
 integer :: nflavor,nspinor,nsppol,nwlo
 integer, allocatable :: found_block(:),label_block(:)
! ************************************************************************

 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol
 nwlo    = paw_dmft%dmft_nwlo

 do iflavor=1,nflavor_max
   block_list(iflavor,:) = iflavor - 1
 end do ! iflavor

 ABI_MALLOC(found_block,(nflavor_max))
 ABI_MALLOC(label_block,(nflavor_max))

 call find_block_structure_mat(matlu(:))

 if (present(hyb)) then
   do i=2,hyb%nmoments-1
     call find_block_structure_mat(hyb%moments(i)%matlu(:))
   end do ! i

   do ifreq=1,nwlo
     call find_block_structure_mat(hyb%oper(ifreq)%matlu(:))
   end do ! ifreq
 end if ! present(hyb)

 siz_block(:,:) = 0

 ! Rename the blocks from 0 to nblocks-1 and build lists
 do iatom=1,natom
   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   nflavor = 2 * (2*lpawu+1)

   found_block(:) = 0

   iblock1 = 0 ! number of blocks
   do iflavor=1,nflavor
     iblock = block_list(iflavor,iatom)
     if (found_block(iblock+1) == 0) then
       label_block(iblock+1) = iblock1
       iblock1 = iblock1 + 1 ! next block index
     end if
     iblock2 = label_block(iblock+1) ! new block index
     block_list(iflavor,iatom) = iblock2
     inner_list(iflavor,iatom) = found_block(iblock+1)
     found_block(iblock+1) = found_block(iblock+1) + 1
     siz_block(iblock2+1,iatom) = siz_block(iblock2+1,iatom) + 1
     flavor_list(found_block(iblock+1),iblock2+1,iatom) = iflavor - 1
   end do ! iflavor
   nblocks(iatom) = iblock1
 end do ! iatom

 ABI_FREE(found_block)
 ABI_FREE(label_block)

 ! Set to 0 the off-diagonal elements that are not kept in a block
 call apply_block_structure_mat(matlu(:))

 if (present(hyb)) then

   do i=2,hyb%nmoments-1
     call apply_block_structure_mat(hyb%moments(i)%matlu(:))
   end do ! i

   do ifreq=1,nwlo
     call apply_block_structure_mat(hyb%oper(ifreq)%matlu(:))
   end do ! ifreq

 end if ! present(hyb)

contains

subroutine find_block_structure_mat(mat)

!Arguments ------------------------------------
 type(matlu_type), intent(inout) :: mat(natom)
!Local variables ------------------------------
 integer :: i,iflavor1,im,im1,isppol,ndim,tndim
! ************************************************************************

  do iatom=1,natom
    lpawu = paw_dmft%lpawu(iatom)
    if (lpawu == -1) cycle
    ndim = 2*lpawu + 1
    tndim = ndim * nspinor
    nflavor = 2 * ndim
    do isppol=1,nsppol
      do im1=1,tndim
        iflavor1 = im1 + (isppol-1)*ndim
        iblock1  = block_list(iflavor1,iatom)
        do im=1,tndim
          iflavor = im + (isppol-1)*ndim
          iblock  = block_list(iflavor,iatom)
          if (iblock == iblock1) cycle
          if (abs(mat(iatom)%mat(im,im1,isppol)) > paw_dmft%dmft_triqs_tol_block) then ! Merge the two blocks
            do i=1,nflavor
              if (block_list(i,iatom) == iblock) block_list(i,iatom) = iblock1
            end do ! i
          end if
        end do ! im
      end do ! im1
    end do ! isppol
    if (nsppol == 1 .and. nspinor == 1) block_list(ndim+1:nflavor,iatom) = block_list(1:ndim,iatom) + ndim
  end do ! iatom

 end subroutine find_block_structure_mat

 subroutine apply_block_structure_mat(mat)

!Arguments ------------------------------------
 type(matlu_type), intent(inout) :: mat(natom)
!Local variables ------------------------------
 integer :: iflavor1,im,im1,isppol,ndim,tndim
! ************************************************************************

 do iatom=1,natom
   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   ndim  = 2*lpawu + 1
   tndim = nspinor * ndim
   do isppol=1,nsppol
     do im1=1,tndim
       iflavor1 = im1 + (isppol-1)*ndim
       iblock1  = block_list(iflavor1,iatom)
       do im=1,tndim
         iflavor = im + (isppol-1)*ndim
         iblock  = block_list(iflavor,iatom)
         if (iblock == iblock1) cycle
         mat(iatom)%mat(im,im1,isppol) = czero
       end do ! im
     end do ! im1
   end do ! isppol
 end do ! iatom

 end subroutine apply_block_structure_mat

 end subroutine find_block_structure
!!***

END MODULE m_forctqmc
!!***
