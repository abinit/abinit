!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_forctqmc
!! NAME
!!  m_forctqmc
!!
!! FUNCTION
!! Prepare CTQMC and call CTQMC
!!
!! COPYRIGHT
!! Copyright (C) 2006-2019 ABINIT group (BAmadon, VPlanes)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
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

MODULE m_forctqmc

 use defs_basis
 use m_errors
 use m_xmpi
 use m_abicore
 use m_Ctqmc
 use m_CtqmcInterface
 use m_Ctqmcoffdiag
 use m_CtqmcoffdiagInterface
 use m_GreenHyb
 use m_data4entropyDMFT

 use m_pawang, only : pawang_type
 use m_crystal, only : crystal_t
 use m_green, only : green_type,occup_green_tau,print_green,printocc_green,spline_fct,copy_green,init_green,destroy_green,&
& int_fct,greenldacompute_green,fourier_green
 use m_paw_dmft, only : paw_dmft_type
 use m_hide_lapack,         only : xginv
 use m_oper, only : oper_type,destroy_oper,init_oper,inverse_oper
 use m_self, only : self_type
 use m_matlu, only : matlu_type,sym_matlu, print_matlu, &
& diag_matlu,init_matlu,destroy_matlu,rotate_matlu,checkdiag_matlu,checkreal_matlu, &
& copy_matlu, diff_matlu, slm2ylm_matlu, shift_matlu, prod_matlu,fac_matlu,&
& add_matlu,printplot_matlu,identity_matlu,zero_matlu
 use m_hu, only : hu_type,rotatevee_hu,vee_ndim2tndim_hu_r
 use m_io_tools, only : flush_unit, open_file
 use m_datafordmft, only : hybridization_asymptotic_coefficient,compute_levels
 use m_special_funcs, only : sbf8
 use m_paw_numeric, only : jbessel=>paw_jbessel

 implicit none

 private

 public :: qmc_prep_ctqmc
 public :: testcode_ctqmc
 public :: testcode_ctqmc_b
 public :: ctqmcoutput_to_green
 public :: ctqmcoutput_printgreen
 public :: ctqmc_calltriqs
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
!! PARENTS
!!      impurity_solve
!!
!! CHILDREN
!!      add_matlu,checkreal_matlu,compute_levels,copy_green,copy_matlu
!!      ctqmc_triqs_run,ctqmcinterface_finalize,ctqmcinterface_init
!!      ctqmcinterface_run,ctqmcinterface_setopts,data4entropydmft_setdocc
!!      destroy_green,destroy_matlu,destroy_oper,diag_matlu,diff_matlu
!!      fac_matlu,flush_unit,fourier_green,hybridization_asymptotic_coefficient
!!      identity_matlu,init_green,init_matlu,init_oper,int_fct,inverse_oper
!!      jbessel,occup_green_tau,print_green,print_matlu,printocc_green
!!      printplot_matlu,prod_matlu,rotate_matlu,rotatevee_hu,sbf8,shift_matlu
!!      slm2ylm_matlu,sym_matlu,testcode_ctqmc,vee_ndim2tndim_hu_r,wrtout,xginv
!!      xmpi_barrier,xmpi_bcast
!!
!! SOURCE

subroutine qmc_prep_ctqmc(cryst_struc,green,self,hu,paw_dmft,pawang,pawprtvol,weiss)

!Arguments ------------------------------------
!scalars
! type(pawang_type), intent(in) :: pawang
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type), intent(inout) :: green  ! MGNAG: This fix the problem with v7[27:29] on nag@petrus
 type(hu_type), intent(in) :: hu(cryst_struc%ntypat)
 type(paw_dmft_type), intent(inout)  :: paw_dmft
 type(pawang_type), intent(in) :: pawang
 integer, intent(in) :: pawprtvol
 type(green_type), intent(inout) :: weiss
 type(self_type), intent(in) :: self

!Local variables ------------------------------
 character(len=500) :: message
 integer :: iatom,ierr,if1,if2,iflavor1,iflavor2,ifreq,im1,ispinor,ispinor1,isppol,itau,itypat,im2,ispinor2
 integer :: lpawu,master,mbandc,natom,nflavor,nkpt,nspinor,nsppol,nsppol_imp,tndim,ispa,ispb,ima,imb
 integer :: nproc,opt_diag,opt_nondiag,testcode,testrot,dmft_nwlo,opt_fk,useylm,nomega,opt_rot
 integer :: ier,rot_type_vee
 complex(dpc) :: omega_current,integral(2,2)
 real(dp) :: doccsum,noise,omega
 logical :: nondiaglevels
! arrays
 real(dp), allocatable :: docc(:,:)
 real(dp), allocatable, target :: gtmp(:,:), levels_ctqmc(:) !modif
 complex(dpc), allocatable :: levels_ctqmc_nd(:,:)
 complex(dpc), allocatable :: hybri_limit(:,:)
 real(dp), allocatable, target :: gtmp_nd(:,:,:)
 real(dp) :: umod(2,2)
 complex(dpc), allocatable :: fw1(:,:),gw_tmp(:,:)
 complex(dpc), allocatable, target :: gw_tmp_nd(:,:,:) !modif
 complex(dpc), allocatable, target :: fw1_nd(:,:,:) !modif
 complex(dpc), allocatable :: gw1_nd(:,:,:)
 complex(dpc), allocatable :: shift(:)
 integer,parameter :: optdb=0
 type(coeff2_type), allocatable :: udens_atoms(:)
! Type    -----------------------------------------
 type(coeff2c_type), allocatable :: eigvectmatlu(:,:)
 type(green_type)  :: weiss_for_rot
 type(matlu_type), allocatable :: dmat_diag(:)
 type(matlu_type), allocatable :: matlu1(:)
 type(matlu_type), allocatable :: matlu2(:)
 type(matlu_type), allocatable :: matlu3(:)
 type(matlu_type), allocatable :: matlu4(:)
 type(matlu_type), allocatable :: identity(:)
 type(matlu_type), allocatable :: level_diag(:)
 type(oper_type)  :: energy_level
 !type(self_type) :: self
! type(green_type) :: gw_loc
 type(CtqmcInterface) :: hybrid   !!! WARNING THIS IS A BACKUP PLAN
 type(CtqmcoffdiagInterface) :: hybridoffdiag   !!! WARNING THIS IS A BACKUP PLAN
 type(green_type) :: greenlda
 type(matlu_type), allocatable  :: hybri_coeff(:)
 integer :: unt,unt2
! Var added to the code for TRIQS_CTQMC test and default value -----------------------------------------------------------
 logical(kind=1) :: leg_measure = .true.

! ************************************************************************
 mbandc=paw_dmft%mbandc
 nkpt=paw_dmft%nkpt
 nsppol=paw_dmft%nsppol
 natom=paw_dmft%natom
 nspinor=paw_dmft%nspinor
 greenlda%whichgreen="LDA"

 call init_green(weiss_for_rot,paw_dmft,opt_oper_ksloc=2)
! weiss_for_rot=>weiss
! call init_green(gw_loc,paw_dmft)
 call copy_green(weiss,weiss_for_rot,opt_tw=2)
!=======================================================================
!== Use one QMC solver   ===============================================
!=======================================================================
 write(message,'(2a)') ch10,'  ===  CT-QMC solver === '
 call wrtout(std_out,message,'COLL')

! Initialise for compiler
 omega_current=czero

! Initialise nproc
 nproc=paw_dmft%nproc

! ======================================
! Allocations: diagonalization and eigenvectors
! ======================================
 ABI_DATATYPE_ALLOCATE(udens_atoms,(natom))
 ABI_DATATYPE_ALLOCATE(eigvectmatlu,(natom,nsppol))
 ABI_DATATYPE_ALLOCATE(dmat_diag,(natom))
 ABI_DATATYPE_ALLOCATE(identity,(natom))
 call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,dmat_diag)
 call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,identity)
 call identity_matlu(identity,natom)
 do iatom=1,cryst_struc%natom
   lpawu=paw_dmft%lpawu(iatom)
   if(lpawu/=-1) then
     tndim=nspinor*(2*lpawu+1)
     do isppol=1,nsppol
       ABI_ALLOCATE(eigvectmatlu(iatom,isppol)%value,(tndim,tndim))
     end do
     ABI_ALLOCATE(udens_atoms(iatom)%value,(2*(2*lpawu+1),2*(2*lpawu+1)))
     dmat_diag(iatom)%mat=czero
   end if
 end do

! ___________________________________________________________________________________
!
!  FIRST PART: DIAGONALISATION AND ROTATIONS.
! ___________________________________________________________________________________
!

! =================================================================
! Choose to diagonalize and how to do it
! =================================================================

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
!  The default value of opt_diag should be 2 for historical reasons (or
!  we decide to change the automatic tests)
!  opt_nondiag should be 0 by default
 opt_diag    = 1
 if(paw_dmft%dmft_solv>=6)  then
   opt_nondiag = 1 ! Use cthyb in triqs or ctqmc in abinit with offdiag terms in F
 else
   opt_nondiag = 0 ! use fast ctqmc in ABINIT without off diagonal terms in F
 end if

 useylm=0
 if(nspinor==2) then
   useylm=1      ! to avoid complex G(tau)
 end if

 !write(6,*) "nspinor,useylm",nspinor,useylm
 if(useylm==0) then
   write(std_out,*) " Slm basis is used (before rotation)"
   rot_type_vee=1 ! for rotatevee_hu
 else if(useylm==1) then
   write(std_out,*) " Ylm basis is used (before rotation)"
   rot_type_vee=4 ! for rotatevee_hu
 end if


! if(useylm==1.and.opt_diag/=1) MSG_ERROR("useylm==1 and opt_diag/=0 is not possible")
 if(hu(1)%jpawu_zero.and.nsppol==2) nsppol_imp=2 ! J=0 and nsppol=2
 if(.not.hu(1)%jpawu_zero.or.nsppol/=2) nsppol_imp=1  ! J/=0 ou nsppol=1
! =================================================================
! Compute LDA Green's function to compare to weiss_for_rot (check)
! =================================================================
! call init_green(greenlda,paw_dmft,opt_oper_ksloc=3)
! call greenldacompute_green(cryst_struc,greenlda,pawang,paw_dmft)
!! call copy_green(greenlda,weiss_for_rot,2)

! =================================================================
! Compute atomic levels
! =================================================================
 call init_oper(paw_dmft,energy_level,opt_ksloc=3)

 ! Compute atomic levels in Slm basis
 ! ----------------------------------
 call compute_levels(cryst_struc,energy_level,self%hdc,pawang,paw_dmft,nondiag=nondiaglevels)

 ! If levels are not diagonal, then diagonalize it (according to
 ! dmftctqmc_basis)
 ! ------------------------------------------------
 if(paw_dmft%dmftctqmc_basis==1) then
   if(nondiaglevels.or.useylm==1) then
     opt_diag=1
     write(message,'(3a)') ch10, "   == Hamiltonian in local basis is non diagonal: diagonalise it",ch10
   else
     opt_diag=0
     write(message,'(5a)') ch10, "   == Hamiltonian in local basis is diagonal in the Slm basis ",ch10 &
&     ,"      CTQMC will use this basis",ch10
   end if
 else if (paw_dmft%dmftctqmc_basis==2) then
   if(nondiaglevels.or.useylm==1) then
     write(message,'(7a)') ch10, "   == Hamiltonian in local basis is non diagonal",ch10, &
&     "   == According to dmftctqmc_basis: diagonalise density matrix",ch10, &
&     "   == Warning : Check that the Hamiltonian is diagonal !",ch10
     opt_diag=2
   else
     write(message,'(5a)') ch10, "   == Hamiltonian in local basis is diagonal in the Slm basis ",ch10 &
&     ,"      CTQMC will use this basis",ch10
     opt_diag=0
   end if
 else if (paw_dmft%dmftctqmc_basis==0) then
   if(nondiaglevels) then
     write(message,'(4a)') ch10, "   == Hamiltonian in local basis is non diagonal",ch10, &
&     "   == According to dmftctqmc_basis: keep this non diagonal basis for the calculation"
   else
     write(message,'(5a)') ch10, "   == Hamiltonian in local basis is diagonal in the Slm basis ",ch10 &
&     ,"      CTQMC will use this basis",ch10
   end if
   opt_diag=0
 end if
 call wrtout(std_out,message,'COLL')
 if(opt_diag==1) then
   write(std_out,*) "  ==  The atomic levels are diagonalized"
 else if(opt_diag==2) then
   write(std_out,*) "  ==  The correlated occupation matrix is diagonalized"
 end if

! =================================================================
! Now, check if diagonalisation is necessary
! =================================================================


! =================================================================
! First rotate to Ylm basis the atomic levels
! =================================================================

 if(useylm==1) then

   ! Rotate from Slm to Ylm the atomic levels
   ! ----------------------------------------
   call slm2ylm_matlu(energy_level%matlu,natom,1,pawprtvol)

   ! Print atomic energy levels in Ylm basis
   ! --------------------------------
   if(pawprtvol>=3) then
     write(message,'(a,a)') ch10, " == Print Energy levels in Ylm basis"
     call wrtout(std_out,message,'COLL')
     call print_matlu(energy_level%matlu,natom,1)
   end if

 end if ! useylm

! ===========================================================================================
! Start for diagonalization of levels/density matrix according to opt_diag
! ===========================================================================================
 !opt_rot=2 ! do it one time before CTQMC
 opt_rot=1 ! do all the rotations successively on all different quantities.
 if(opt_diag==1.or.opt_diag==0) then


   if(opt_diag==1) then
! =================================================================
! Diagonalize atomic levels
! =================================================================
     ABI_DATATYPE_ALLOCATE(level_diag,(natom))
     call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,level_diag)

     ! Diagonalise atomic levels (opt_real is necessary, because
     ! rotation must be real in order for the occupations and Green's
     ! function to be real)
     ! ---------------------------------------------------------------
     call diag_matlu(energy_level%matlu,level_diag,natom,&
&     prtopt=pawprtvol,eigvectmatlu=eigvectmatlu,nsppol_imp=nsppol_imp,optreal=1,&
&     test=paw_dmft%dmft_solv)  ! temporary: test should be extended to all cases.

!     call rotate_matlu(energy_level%matlu,eigvectmatlu,natom,3,1)
!       write(message,'(a,2x,a,f13.5)') ch10,&
!&       " == Print first Diagonalized Energy levels for Fermi Level=",paw_dmft%fermie
!       call wrtout(std_out,message,'COLL')
!       call print_matlu(energy_level%matlu,natom,1,compl=1,opt_exp=1)

     if(opt_rot==1) call copy_matlu(level_diag,energy_level%matlu,natom)


     call destroy_matlu(level_diag,natom)
     ABI_DATATYPE_DEALLOCATE(level_diag)

     ! Print diagonalized levels
     ! --------------------------
     if(pawprtvol>=3) then
       write(message,'(a,2x,a,f13.5)') ch10,&
&       " == Print Diagonalized Energy levels for Fermi Level=",paw_dmft%fermie
       call wrtout(std_out,message,'COLL')
       call print_matlu(energy_level%matlu,natom,1,compl=1,opt_exp=1)
     else
       write(message,'(a,2x,a,f13.5)') ch10,&
&       " == Energy levels Diagonalized for Fermi Level=",paw_dmft%fermie
       call wrtout(std_out,message,'COLL')
     end if

     call rotatevee_hu(cryst_struc,hu,nspinor,nsppol,pawprtvol,eigvectmatlu,udens_atoms,rot_type_vee)

   else if (opt_diag==0) then
     do iatom=1,cryst_struc%natom
       lpawu=paw_dmft%lpawu(iatom)
       itypat=cryst_struc%typat(iatom)
       if(lpawu/=-1) then
       !  write(6,*) size(udens_atoms(iatom)%value)
       !  write(6,*) size(hu(itypat)%udens)
       !  write(6,*) udens_atoms(iatom)%value
       !  write(6,*) hu(itypat)%udens
         udens_atoms(iatom)%value=hu(itypat)%udens
       end if
     end do
   end if
  ! call rotatevee_hu(cryst_struc,hu,nspinor,nsppol,pawprtvol,eigvectmatlu,udens_atoms)

 else if(opt_diag==2) then
! =================================================================
! Diagonalizes density matrix and keep eigenvectors in eigvectmatlu
! =================================================================

   ! Print density matrix before diagonalization
   ! -------------------------------------------
   if(pawprtvol>=3) then
     write(message,'(a,2x,a)') ch10,        " == Density Matrix before diagonalisation ="
     call wrtout(std_out,message,'COLL')
     !MGNAG: This call is wrong if green has intent(out), now we use intent(inout)
     call print_matlu(green%occup%matlu,natom,1)
   end if

!!  checkstop: we can have two different diagonalisation basis for the up and dn
!!  but one use the same basis, unless the error is really to large(>0.1)

   ! Diagonalize density matrix
   ! ---------------------------
   call diag_matlu(green%occup%matlu,dmat_diag,natom,&
&   prtopt=4,eigvectmatlu=eigvectmatlu,nsppol_imp=nsppol_imp,checkstop=.false.)

   ! Print diagonalized density matrix
   ! ----------------------------------
   if(pawprtvol>=3) then
     write(message,'(a,2x,a)') ch10,&
&     " == Diagonalized Density Matrix in the basis used for QMC ="
     call wrtout(std_out,message,'COLL')
     call print_matlu(dmat_diag,natom,1)

     !write(message,'(2a,i3,13x,a)') ch10,'    ==  Rotation of interaction matrix =='
     !call wrtout(std_out,message,'COLL')
   end if

   !if (.not.hu(1)%jpawu_zero) &
   !MSG_WARNING("In qmc_prep_ctqmc J/=0 and rotation matrix not rotated")
!  Rotate interaction.
!   call rotatevee_hu(cryst_struc,hu,nspinor,nsppol,pawprtvol,eigvectmatlu,udens_atoms)
   call rotatevee_hu(cryst_struc,hu,nspinor,nsppol,pawprtvol,eigvectmatlu,udens_atoms,rot_type_vee)

 end if
! ===========================================================================================
! END Of diagonalization
! ===========================================================================================

 call flush_unit(std_out)

! ===========================================================================================
! Broadcast matrix of rotation from processor 0 to the other
! In case of degenerate levels, severals rotations are possible. Here we
! choose the rotation of proc 0. It is arbitrary.
! ===========================================================================================
 do iatom=1,cryst_struc%natom
   lpawu=paw_dmft%lpawu(iatom)
   if(lpawu/=-1) then
     tndim=nspinor*(2*lpawu+1)
     do isppol=1,nsppol
       call xmpi_bcast(eigvectmatlu(iatom,isppol)%value,0,paw_dmft%spacecomm,ier)
     end do
   end if
 end do


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
 if(opt_diag==2.and.opt_rot==1) call rotate_matlu(energy_level%matlu,eigvectmatlu,natom,3,1)

   ! Print atomic levels
   ! -------------------
 if(pawprtvol>=3) then
   write(message,'(a,2x,a,f13.5)') ch10," == Print Energy levels after rotation"
   call wrtout(std_out,message,'COLL')
   call print_matlu(energy_level%matlu,natom,1)
 else
   write(message,'(a,2x,a,f13.5)') ch10," == CT-QMC Energy levels rotated"
   call wrtout(std_out,message,'COLL')
 end if

!====================================================================
! If levels were diagonalized before, then rotate density matrix for
! information.
!====================================================================
 if(opt_diag==1) then
   ABI_DATATYPE_ALLOCATE(matlu1,(natom))
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,matlu1)
   call copy_matlu(green%occup%matlu,matlu1,natom)
   if(pawprtvol>=3) then
     write(message,'(a,2x,a)') ch10,&
&     " == Occupations before rotations"
     call wrtout(std_out,message,'COLL')
     call print_matlu(green%occup%matlu,natom,1)
   end if

   ! 1) rotate density matrix to Ylm basis
   ! --------------------------------------
   if(useylm==1) then
     call slm2ylm_matlu(matlu1,natom,1,pawprtvol)
     if(pawprtvol>=3) then
       write(message,'(a,a)') ch10, " == Print occupations in Ylm basis"
       call wrtout(std_out,message,'COLL')
       call print_matlu(matlu1,natom,1)
     end if
   end if

   ! 2) rotate density matrix to rotated basis
   ! -------------------------------------------
   if(opt_rot==1.or.opt_rot==2) call rotate_matlu(matlu1,eigvectmatlu,natom,3,1)
   write(message,'(a,2x,a,f13.5)') ch10," == Rotated occupations (for information)"
   call wrtout(std_out,message,'COLL')
   call print_matlu(matlu1,natom,1,compl=1)
   call checkreal_matlu(matlu1,natom,tol10)
   call destroy_matlu(matlu1,natom)
   ABI_DATATYPE_DEALLOCATE(matlu1)

 end if

 call flush_unit(std_out)


! =================================================================
! Rotate weiss function according to eigenvectors.
! =================================================================
!!!stop
  ! Rotate Weiss function first in Ylm basis
  ! -------------------------------------------------------------------
 if(useylm==1) then
   write(message,'(a,2x,a)') ch10, " == Rotation of weiss and greenlda in the Ylm Basis="
   call wrtout(std_out,message,'COLL')
   do ifreq=1,paw_dmft%dmft_nwlo
     call slm2ylm_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,1,0)
     call slm2ylm_matlu(weiss%oper(ifreq)%matlu,natom,1,0)
     ! call slm2ylm_matlu(greenlda%oper(ifreq)%matlu,natom,1,0)
   end do
 end if

 if(pawprtvol>=3) then
   !   write(message,'(a,2x,a,f13.5)') ch10,& ! debug
   !   " == Print weiss for small freq 1 before rot" ! debug
   !   call wrtout(std_out,message,'COLL') ! debug
   !   call print_matlu(weiss_for_rot%oper(1)%matlu,natom,1) !  debug

    ! Print Weiss function
    ! --------------------
   write(message,'(a,2x,a,f13.5)') ch10,& ! debug
&  " == Print weiss for 1st freq before rot" ! debug
   call wrtout(std_out,message,'COLL') ! debug
   call print_matlu(weiss_for_rot%oper(1)%matlu,natom,1,compl=1) !  debug
   write(message,'(a,2x,a,f13.5)') ch10,& ! debug
&  " == Print weiss for last freq before rot" ! debug
   call wrtout(std_out,message,'COLL') ! debug
   call print_matlu(weiss_for_rot%oper(paw_dmft%dmft_nwlo)%matlu,natom,1,compl=1) !  debug
!    write(message,'(a,2x,a,f13.5)') ch10,& ! debug
!&   " == Print LDA G for 1st freq before rot" ! debug
!    call wrtout(std_out,message,'COLL') ! debug
!    call print_matlu(greenlda%oper(1)%matlu,natom,1,compl=1,opt_exp=2) !  debug
!    write(message,'(a,2x,a,f13.5)') ch10,& ! debug
!&   " == Print LDA G for last freq before rot" ! debug
!    call wrtout(std_out,message,'COLL') ! debug
!    call print_matlu(greenlda%oper(paw_dmft%dmft_nwlo)%matlu,natom,1,compl=1,opt_exp=2) !  debug
 end if

 if(opt_diag/=0) then
   ! Rotate Weiss function from the Slm (or Ylm) to the basis of diagonalisation
   ! -------------------------------------------------------------------
   write(message,'(a,2x,a)') ch10, " == Rotation of weiss ="
   call wrtout(std_out,message,'COLL')
   do ifreq=1,paw_dmft%dmft_nwlo
     if(opt_rot==1) call rotate_matlu(weiss_for_rot%oper(ifreq)%matlu,eigvectmatlu,natom,3,1)
     if(opt_rot==1) call rotate_matlu(weiss%oper(ifreq)%matlu,eigvectmatlu,natom,3,1)
!    call checkdiag_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,tol6)
   end do

   if(paw_dmft%myproc .eq. mod(nproc+1,nproc)) then
     if (open_file(trim(paw_dmft%filapp)//"_atom__G0w_.dat", message, newunit=unt) /= 0) then
       MSG_ERROR(message)
     end if
     do ifreq=1,paw_dmft%dmft_nwlo
       write(unt,'(29f21.14)') paw_dmft%omega_lo(ifreq),&
&       (((weiss_for_rot%oper(ifreq)%matlu(1)%mat(im1,im1,isppol,ispinor,ispinor),&
&       im1=1,3),ispinor=1,nspinor),isppol=1,nsppol)
     end do
     close(unt)
   end if

   call flush_unit(std_out)
   if(pawprtvol>=3) then
     write(message,'(a,2x,a,f13.5)') ch10,& ! debug
&    " == Print weiss for small freq 1 after rot" ! debug
     call wrtout(std_out,message,'COLL') ! debug
     call print_matlu(weiss_for_rot%oper(1)%matlu,natom,1,compl=1) !  debug
     write(message,'(a,2x,a,f13.5)') ch10,&   ! debug
&    " == Print weiss for last freq after rot"   ! debug
     call wrtout(std_out,message,'COLL')   ! debug
     call print_matlu(weiss_for_rot%oper(paw_dmft%dmft_nwlo)%matlu,natom,1,compl=1) ! debug
   end if

!   ! Rotate LDA Green's function first in Ylm basis then in the rotated basis and compare to weiss_for_rot
!   ! -----------------------------------------------------------------------------------------------------
!   write(message,'(a,2x,a)') ch10, " == Rotation of greenlda ="
!   call wrtout(std_out,message,'COLL')
!   do ifreq=1,paw_dmft%dmft_nwlo
!     if(opt_rot==1) call rotate_matlu(greenlda%oper(ifreq)%matlu,eigvectmatlu,natom,3,1)
!     call diff_matlu("Weiss_for_rot","greenlda",weiss_for_rot%oper(ifreq)%matlu,greenlda%oper(ifreq)%matlu,natom,1,tol14)
!!    call checkdiag_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,tol6)
!   end do
!   if(pawprtvol>=3) then
!     write(message,'(a,2x,a,f13.5)') ch10,& ! debug
!&    " == Print greenlda for small freq 1 after rot" ! debug
!     call wrtout(std_out,message,'COLL') ! debug
!     call print_matlu(greenlda%oper(1)%matlu,natom,1,compl=1,opt_exp=2) !  debug
!     write(message,'(a,2x,a,f13.5)') ch10,&   ! debug
!&    " == Print greenlda for last freq after rot"   ! debug
!     call wrtout(std_out,message,'COLL')   ! debug
!     call print_matlu(greenlda%oper(paw_dmft%dmft_nwlo)%matlu,natom,1,compl=1,opt_exp=2) ! debug
!   end if
!   call flush_unit(std_out)
 end if

! =================================================================
! Compute analytic limit of hybridization and rotate it
! =================================================================
 ABI_DATATYPE_ALLOCATE(hybri_coeff,(paw_dmft%natom))
 call init_matlu(paw_dmft%natom,paw_dmft%nspinor,paw_dmft%nsppol,paw_dmft%lpawu,hybri_coeff)
 !write(6,*)"hybri1",hybri_coeff(1)%mat(1,1,1,1,1),paw_dmft%natom,cryst_struc%natom

 ! Compute analytical C_ij such that F_ij -> C_ij/iw_n
 ! ---------------------------------------
 call hybridization_asymptotic_coefficient(cryst_struc,paw_dmft,pawang,hybri_coeff)
 write(message,'(a,2x,a)') ch10," == Coeff analytical C_ij such that F -> C_ij/iw_n for large frequency"
 call wrtout(std_out,message,'COLL')

 ! Print analytical C_ij (not rotated)
 ! ---------------------------------------
 call print_matlu(hybri_coeff,natom,1)

 ! Rotate analytical C_ij in Ylm basis
 ! ---------------------------------------
 if(useylm==1) call slm2ylm_matlu(hybri_coeff,natom,1,pawprtvol)
 if(opt_diag/=0)  then

 ! Rotate analytical C_ij in rotated basis
 ! ---------------------------------------
   if(opt_rot==1.or.opt_rot==2) call rotate_matlu(hybri_coeff,eigvectmatlu,natom,3,1)

 ! Print analytical C_ij (rotated)
 ! ---------------------------------------
   write(message,'(a,2x,a)') ch10," == Coeff analytical C_ij such that F -> C_ij/iw_n after rotation"
   call wrtout(std_out,message,'COLL')
   call print_matlu(hybri_coeff,natom,1,compl=1,opt_exp=1)
 end if

! =================================================================
! Check if rotation is properly done.
! =================================================================
 if(3==4) then
   write(message,'(a,2x,a)') ch10,&
&   " == Print  dmat before rot"
   call wrtout(std_out,message,'COLL')
   call print_matlu(green%occup%matlu,natom,1)
   if(useylm==1) call slm2ylm_matlu(green%occup%matlu,natom,1,pawprtvol)
   if(opt_rot==1) call rotate_matlu(green%occup%matlu,eigvectmatlu,natom,3,1)
   write(message,'(a,2x,a)') ch10,&
&   " == Print  dmat after rot"
   call wrtout(std_out,message,'COLL')
   call print_matlu(green%occup%matlu,natom,1)

   write(message,'(2a)') ch10,' QMC STOP: DEBUG'
   call wrtout(std_out,message,'COLL')
   MSG_ERROR(message)
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

 master=0

! =================================================================
! Print out
! =================================================================

! Print Weiss
! -------------
 if(paw_dmft%dmft_prgn==1) then
   call print_green('Weiss_diag',weiss_for_rot,1,paw_dmft,pawprtvol=1,opt_wt=1,opt_decim=1)
 end if

 write(message,'(a,2x,a,f13.5)') ch10,&
& " == Preparing data for CTQMC"
 call wrtout(std_out,message,'COLL')

! Print Rotate Weiss for 1st and last frequencies
! ------------------------------------------------
 if (pawprtvol>=3) then
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print rotated weiss function for small freq in the rotated basis"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(weiss_for_rot%oper(1)%matlu,natom,1,compl=1)  ! debug
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print rotated weiss function for largest freq in the rotated basis"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(weiss_for_rot%oper(paw_dmft%dmft_nwlo)%matlu,natom,1,compl=1)  ! debug
 end if

! =================================================================
!  VARIABLES FOR CTQMC TESTS
 testcode = 0
 testrot  = 0
 opt_fk=0 ! for developpers to check Fourier transform and computes G0(tau)
 opt_fk=1 ! usual case: for real calculations
! =================================================================

! ___________________________________________________________________________________
!
!  SECOND PART : BUILT HYBRIDIZATION FROM G0
! ___________________________________________________________________________________
!
! ===========================================================================================
! Compute inverse of weiss  and compute hybridization
! ===========================================================================================

! Compute inverse of weiss  for each Frequency
! ----------------------------------------------
 do ifreq=1,paw_dmft%dmft_nwlo
   ABI_DATATYPE_ALLOCATE(matlu1,(natom))
   ABI_DATATYPE_ALLOCATE(matlu2,(natom))
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,matlu1)
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,matlu2)

   call copy_matlu(weiss_for_rot%oper(ifreq)%matlu,matlu1,natom)

   ! Print G_0(iw_n)
   ! ----------------
   if(optdb==1) call printplot_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,paw_dmft%omega_lo(ifreq),"go",60000,imre=1)

   ! Compute G_0^-1
   ! -------------------------------------------
   ! if opt_fk=1 or testcode/=0  Do the inversion
   ! if opt_fk=0                 Do not inverse.
   ! If testcode=2 and opt_fk=0  Do the inversion
   ! If testcode=1 and opt_fk=0  Do the inversion but no effect, because it will nevertheless be erased
   ! If opt_fk=1                 Do the inversion
   ! -------------------------------------------
   if(optdb==1) call printplot_matlu(matlu1,natom,paw_dmft%omega_lo(ifreq),"weiss",12000,imre=1)
   if(opt_fk==1.or.testcode/=0) call inverse_oper(weiss_for_rot%oper(ifreq),option=1,prtopt=1)

   ! Print G_0^-1(iw_n)
   ! ----------------
   if(optdb==1) call printplot_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,paw_dmft%omega_lo(ifreq),"goinv",70000,imre=1)

   if(pawprtvol>=4.or.ifreq==paw_dmft%dmft_nwlo) then
     if(opt_fk==1.or.testcode/=0) then
      ! Check inversion : do the product
      ! ----------------------------------------------
       call prod_matlu(weiss_for_rot%oper(ifreq)%matlu,matlu1,matlu2,natom)
       write(message,'(a,2x,a,i7)') ch10,&  ! debug
&      " == Print product of  weiss times invers for freq",ifreq
       call wrtout(std_out,message,'COLL')  ! debug
       call print_matlu(matlu2,natom,1)  ! debug
     end if
   end if

   call destroy_matlu(matlu1,natom)
   call destroy_matlu(matlu2,natom)
   ABI_DATATYPE_DEALLOCATE(matlu1)
   ABI_DATATYPE_DEALLOCATE(matlu2)
 end do

 ! Copy weiss_for_rot into weiss
 ! -------------------------------
 !call copy_matlu(weiss_for_rot%oper(ifreq)%matlu,weiss%oper(ifreq)%matlu,natom)


 ! Print G_0^-1 for 1st and last frequencies.
 ! -----------------------------------------
 if(pawprtvol>=3) then
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print G_0^-1 for small freq in the rotated basis"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(weiss_for_rot%oper(1)%matlu,natom,1)  ! debug
   write(message,'(a,2x,a,e18.10,a)') ch10,&   ! debug
&  " == Print G_0^-1 for last freq in the rotated basis (last freq=", paw_dmft%omega_lo(paw_dmft%dmft_nwlo),")"  ! debug
   call wrtout(std_out,message,'COLL')   ! debug
   call print_matlu(weiss_for_rot%oper(paw_dmft%dmft_nwlo)%matlu,natom,1,compl=1) ! debug
 end if

! Substract frequency from diagonal part
! ======================================

 ABI_ALLOCATE(shift,(natom))
 do ifreq=1,paw_dmft%dmft_nwlo
   shift(:)=cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)

!  write(5555,'(400e17.4)') paw_dmft%omega_lo(ifreq),((((((weiss_for_rot%oper(ifreq)%matlu(1)%mat&
!  & (im,im1,isppol,ispinor,ispinor1)-cmplx(0.d0,paw_dmft%omega_lo(ifreq),kind=dp)),im=1,2*3+1),&
!&      im1=1,2*3+1),isppol=1,nsppol),ispinor=1,nspinor),ispinor1=1,nspinor)

  ! Compute G_0^-1-iw_n
  ! --------------------
   if(opt_fk==1) call shift_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,shift)

  ! Compute -G_0^-1+iw_n
  ! --------------------
   if(opt_fk==1) call fac_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,-cone)

  ! Print -G_0^-1+iw_n
  ! --------------------
   if(optdb==1) then
     call printplot_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,paw_dmft%omega_lo(ifreq),"G0inv_minus_omega",20000,imre=1)
   end if
 end do

 ! Print -G_0^+1-iw_n=(F-levels) for last freq in the rotated basis"
 ! ------------------------------------------------------------------
 ABI_DEALLOCATE(shift)
 if(pawprtvol>=3) then
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print G_0^-1-iw_n=-(F-levels) for last freq in the rotated basis"  ! debug
   call wrtout(std_out,message,'COLL')   ! debug
   call print_matlu(weiss_for_rot%oper(paw_dmft%dmft_nwlo)%matlu,natom,1,compl=1) ! debug
 end if

! Check numerical limit of F(i_wn)*iw_n (can be used also to compute F )
! ======================================

 if(opt_nondiag==1) then
   ABI_DATATYPE_ALLOCATE(matlu1,(natom))
   ABI_DATATYPE_ALLOCATE(matlu2,(natom))
   ABI_DATATYPE_ALLOCATE(matlu3,(natom))
   ABI_DATATYPE_ALLOCATE(matlu4,(natom))
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,matlu1)
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,matlu2)
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,matlu3)
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,matlu4)

   write(message,'(a,2x,a,f13.5)') ch10,  " == energy_levels"
   call wrtout(std_out,message,'COLL')
   call print_matlu(energy_level%matlu,natom,1,opt_exp=2,compl=1)

   do ifreq=paw_dmft%dmft_nwlo,1,-1 ! necessary to have matlu4 computed for the max frequency and available for all frequency.
   !do ifreq=paw_dmft%dmftqmc_l,1,-1 ! necessary to have matlu4 computed for the max frequency and available for all frequency.
      ! Compute F (substract levels) for max frequency
      ! -----------------------------------------------
     call add_matlu(weiss_for_rot%oper(ifreq)%matlu,energy_level%matlu,matlu1,natom,-1)

      ! Print F(iw_n)=-(G_0^-1-iw_n+levels)  for last frequency.
      ! --------------------------------------------------------
     if(ifreq==paw_dmft%dmft_nwlo.or.ifreq==paw_dmft%dmftqmc_l) then
       write(message,'(a,2x,a,i4,a,f13.5,a)') ch10, &
&       " == Print F(iw_n)=-(G_0^-1-iw_n+levels) for freq nb",ifreq," (=",paw_dmft%omega_lo(ifreq),")"
       call wrtout(std_out,message,'COLL')
       call print_matlu(matlu1,natom,1,opt_exp=1,compl=1)
     end if
     if(optdb==1) call printplot_matlu(matlu1,natom,paw_dmft%omega_lo(ifreq),"Hybridization",10000,imre=1)

      ! Put F in weiss_for_rot -> CTQMC
      ! -------------------------------
     if(opt_rot==2) call rotate_matlu(weiss_for_rot%oper(ifreq)%matlu,eigvectmatlu,natom,3,1)
!   The following line will produce directly the weiss function for the CTQMC code
     if(opt_fk==1) call copy_matlu(matlu1,weiss_for_rot%oper(ifreq)%matlu,natom)


      ! Multiply F by frequency
      ! ------------------------
     call copy_matlu(matlu1,matlu2,natom)
     call fac_matlu(matlu1,natom,cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp))
     if(ifreq==paw_dmft%dmft_nwlo.or.ifreq==paw_dmft%dmftqmc_l) then
       write(message,'(a,2x,a,i4,a,f13.5,a)') ch10, &
&       " == Print numerical C_ij = F(iw_n)*iw_n for freq nb",ifreq," (=",paw_dmft%omega_lo(ifreq),")"
       call wrtout(std_out,message,'COLL')
       call print_matlu(matlu1,natom,1,opt_exp=1,compl=1)
     end if
     if(optdb==1) call printplot_matlu(matlu1,natom,paw_dmft%omega_lo(ifreq),"cij",72800,imre=1)
     !call rotate_matlu(matlu1,eigvectmatlu,natom,3,1)

     if(ifreq==paw_dmft%dmft_nwlo.or.ifreq==paw_dmft%dmftqmc_l) then
       write(message,'(a,2x,a,i4,a,f13.5,a)') ch10, &
&       " == Print numerical after back rotation C_ij = F(iw_n)*iw_n for freq nb",ifreq," (=",paw_dmft%omega_lo(ifreq),")"
       call wrtout(std_out,message,'COLL')
       call print_matlu(matlu1,natom,1,opt_exp=1,compl=1)
     end if
     if(optdb==1) call printplot_matlu(matlu1,natom,paw_dmft%omega_lo(ifreq),"cij_rotated",72900,imre=1)

      ! Built C_ij/iw_n
      ! ------------------------
     call copy_matlu(hybri_coeff,matlu1,natom)
     call fac_matlu(matlu1,natom,1.d0/cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp))
     if(optdb==1) call printplot_matlu(matlu1,natom,paw_dmft%omega_lo(ifreq),"cij_over_omega",72000)
    ! if(ifreq==paw_dmft%dmft_nwlo) then
    !   write(message,'(a,2x,a,f13.5)') ch10,  " == Print numerical C_ij/iw_n for frequency",paw_dmft%omega_lo(ifreq)
    !   call wrtout(std_out,message,'COLL')
    !   call print_matlu(matlu1,natom,1,opt_exp=1,compl=1)
    ! endif

      ! For test: put C_ij/i_wn into weiss_for_rot
      ! --------------------------------------------
     !call copy_matlu(matlu1,weiss_for_rot%oper(ifreq)%matlu,natom,opt_non_diag=1)

      ! Compute Hybri - C_ij/iw_n
      ! ------------------------
     call add_matlu(matlu2,matlu1,matlu3,natom,-1)

      ! Print Hybri - C_ij/iw_n
      ! ------------------------
     if(optdb==1) call printplot_matlu(matlu3,natom,paw_dmft%omega_lo(ifreq),"hybri_minus_asymp",74000,imre=1)

      ! Multiply (F-C_ij/i_wn) by (iw_n)**2 to find D_ij such that (F-C_ij/i_wn) -> D_ij/(iw_n)^2 only for last frequency.
      ! ------------------------------------------------------------------------------------------------------------------
     call copy_matlu(matlu3,matlu2,natom)
     call fac_matlu(matlu2,natom,cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)**2)
     if(optdb==1) call printplot_matlu(matlu2,natom,paw_dmft%omega_lo(ifreq),"fminuscijtimesw2",75000,imre=1)
     if(ifreq==paw_dmft%dmft_nwlo.or.ifreq==paw_dmft%dmftqmc_l) then
       call copy_matlu(matlu2,matlu4,natom)
       write(message,'(a,2x,a,i4,a,f13.5,a)') ch10, &
&       " == Print numerical (F(iw_n)-C_ij/iw_n)%iw_n^2 for freq nb",ifreq," (=",paw_dmft%omega_lo(ifreq),")"
       call wrtout(std_out,message,'COLL')
       call print_matlu(matlu4,natom,1)
     end if

      ! Built C_ij/iw_n+D_ij/(iw_n)^2
      ! ------------------------
     call copy_matlu(matlu4,matlu3,natom,opt_re=1)
     call fac_matlu(matlu3,natom,1.d0/cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)**2)
     call add_matlu(matlu1,matlu3,matlu2,natom,1)
     if(optdb==1) call printplot_matlu(matlu2,natom,paw_dmft%omega_lo(ifreq),"cij_w_plus_dij_w2",72700,imre=1)
      ! For test: put C_ij/i_wn +D_ij/(iw_n)^2 into weiss_for_rot
      ! --------------------------------------------
     !call copy_matlu(matlu2,weiss_for_rot%oper(ifreq)%matlu,natom,opt_non_diag=1)


   end do

   call destroy_matlu(matlu1,natom)
   call destroy_matlu(matlu2,natom)
   call destroy_matlu(matlu3,natom)
   call destroy_matlu(matlu4,natom)
   ABI_DATATYPE_DEALLOCATE(matlu1)
   ABI_DATATYPE_DEALLOCATE(matlu2)
   ABI_DATATYPE_DEALLOCATE(matlu3)
   ABI_DATATYPE_DEALLOCATE(matlu4)
 end if ! if opt_nondiag=1

! =========================================================================================
! Start big loop over atoms to compute hybridization and do the CTQMC
! =========================================================================================
 do iatom=1,cryst_struc%natom
   green%ecorr_qmc(iatom)=zero
   itypat=cryst_struc%typat(iatom)
   lpawu=paw_dmft%lpawu(iatom)
   tndim=2*lpawu+1
   if(lpawu/=-1) then

     nflavor=2*(tndim)
     if(testcode>=1) then
       nflavor=2
       if(testcode==2) then
         ispa=1
         ispb=2
         if(nspinor==1) ispb=1
         ima=1
         imb=1
         if(tndim>4) then
           ima=5 ! row
           imb=4 ! column
         end if
       end if
     end if

     ABI_ALLOCATE(fw1,(paw_dmft%dmft_nwlo,nflavor))
     ABI_ALLOCATE(fw1_nd,(paw_dmft%dmft_nwlo,nflavor,nflavor))
     ABI_ALLOCATE(levels_ctqmc,(nflavor))
     ABI_ALLOCATE(levels_ctqmc_nd,(nflavor,nflavor))
     levels_ctqmc_nd=czero
     ABI_ALLOCATE(hybri_limit,(nflavor,nflavor))
     hybri_limit=czero
     fw1_nd=czero
     fw1=czero

     ! =================================================================
     ! Put hybridization in new arrays for CTQMC
     ! =================================================================
     if (testcode==0) then
       iflavor1=0
       iflavor2=0

       do isppol=1,nsppol
         do ispinor1=1,nspinor
           do ispinor2=1,nspinor
             do im1=1,tndim
               do im2=1,tndim

                 ! first diagonal terms whatever opt_nondiag
                 iflavor1=im1+tndim*(ispinor1-1)+tndim*(isppol-1)
                 iflavor2=im2+tndim*(ispinor2-1)+tndim*(isppol-1)

                 if ( iflavor1==iflavor2 ) then

                   ! Put weiss_for_rot in fw1
                   do ifreq=1,paw_dmft%dmft_nwlo
                     if(opt_fk==1) then
                       fw1(ifreq,iflavor1)= weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(im1,im1,isppol,ispinor1,ispinor1)
                     else if (opt_fk==0) then
                       fw1(ifreq,iflavor1)= weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(im1,im1,isppol,ispinor1,ispinor1)
                     end if
                   end do
                   fw1_nd(:,iflavor1,iflavor1)=fw1(:,iflavor1)

                   levels_ctqmc(iflavor1)=real(energy_level%matlu(iatom)%mat(im1,im1,isppol,ispinor1,ispinor1),kind=dp)
                   hybri_limit(iflavor1,iflavor1)=hybri_coeff(iatom)%mat(im1,im1,isppol,ispinor1,ispinor1)


                   ! case nsppol=nspinor=1
                   if(nsppol==1.and.nspinor==1) then
                     fw1(:,iflavor1+tndim)=fw1(:,iflavor1)
                     fw1_nd(:,iflavor1+tndim,iflavor1+tndim)=fw1(:,iflavor1)
                     levels_ctqmc(iflavor1+tndim)=levels_ctqmc(iflavor1)
                     hybri_limit(iflavor1+tndim,iflavor1+tndim)=hybri_limit(iflavor1,iflavor1)
                   end if

                 ! off diagonal terms
                 else

                   ! Put weiss_for_rot in fw1_nd
                   do ifreq=1,paw_dmft%dmft_nwlo
                     if(opt_fk==1) then
                       fw1_nd(ifreq,iflavor1,iflavor2)= &
&                       weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)
                     else if (opt_fk==0) then
                       fw1_nd(ifreq,iflavor1,iflavor2)= &
&                       weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)
                     end if
                   end do
                   hybri_limit(iflavor1,iflavor2)=hybri_coeff(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)

                   ! case nsppol=nspinor=1
                   if(nsppol==1.and.nspinor==1) then
                     fw1_nd(:,iflavor1+tndim,iflavor2+tndim) = fw1_nd(:,iflavor1,iflavor2)
                     hybri_limit(iflavor1+tndim,iflavor2+tndim)=hybri_limit(iflavor1,iflavor2)
                   end if

                 end if

! <  / HACK >
               end do !im2
             end do !im1
           end do  !ispinor2
         end do  !ispinor1
       end do  !isppol
! < HACK >
       ! JB. On 1000 cpus this can not work since all CPU try to open/write the files
       ! Action : Don't print it or check only one cpu does it.

       if(pawprtvol>=10000000) then
         write(message,'(a,2x,a)') ch10,  " == Hybri for all flavors for CTQMC "
         call wrtout(std_out,message,'COLL')
         do iflavor1=1,nflavor
           write(message,'(4x,14(2e14.5,2x))') (hybri_limit(iflavor1,iflavor2),iflavor2=1,nflavor)
           call wrtout(std_out,message,'COLL')
         end do

         if (open_file('Hybri_cijoveromega',message, newunit=unt, status='unknown', form='formatted') /= 0) then
           MSG_ERROR(message)
         end if
         if (open_file('Hybri',message,newunit=unt2,status='unknown',form='formatted') /= 0) then
           MSG_ERROR(message)
         end if
         do ifreq=1,paw_dmft%dmft_nwlo
           !  weiss_for_rot is G_0^-1-iw_n=-(F-levels)
           if(optdb==1) call printplot_matlu(weiss_for_rot%oper(ifreq)%matlu,natom,paw_dmft%omega_lo(ifreq),"weissbefore112",30000)
         end do
         do iflavor1=1,nflavor
           do iflavor2=1,nflavor
             do ifreq=1,paw_dmft%dmft_nwlo
               omega=pi*paw_dmft%temp*(two*float(ifreq)-1)
               ! fw1_nd is -G_0^+1-iw_n=(F-levels)
               write(unt,'(300e16.5)') paw_dmft%omega_lo(ifreq)&
&               ,fw1_nd(ifreq,iflavor1,iflavor2)-hybri_limit(iflavor1,iflavor2)/cmplx(0.d0,paw_dmft%omega_lo(ifreq),kind=dp)
               write(unt2,'(300e16.5)') paw_dmft%omega_lo(ifreq),fw1_nd(ifreq,iflavor1,iflavor2)
             end do
             write(unt,*)
             write(unt2,*)
           end do
         end do
         close(unt)
         close(unt2)
       end if
     end if ! testcode
   ! </ HACK >

     ! ====================================================================================
     !  TEST
     !  For testing purpose, built ultra simple hybridization (constant in
     !  imaginary time or very simple) or extract some part of the calculated hybridization
     ! ====================================================================================
     if(testcode>=1) then
       dmft_nwlo=paw_dmft%dmft_nwlo
       paw_dmft%dmft_nwlo=paw_dmft%dmftqmc_l
       ABI_ALLOCATE(gw1_nd,(paw_dmft%dmft_nwlo,nflavor,nflavor))
       gw1_nd=czero

       !  Call testcode_ctqmc: built simple hybridization
       !--------------------------------------------------
       if (testcode==1) then
         call testcode_ctqmc(paw_dmft%dmftqmc_l,fw1_nd,fw1,gtmp_nd,gw_tmp_nd,&
&         levels_ctqmc,hybri_limit,nflavor,1,paw_dmft%temp,testrot,testcode,umod)
       !  Select 2x2 hybridization matrix from the current larger matrix
       !  ima and imb are defined above.
       !----------------------------------------------------------------
       else if (testcode==2) then
         !close(unt)
         !close(unt2)
         call testcode_ctqmc_b(energy_level,hybri_coeff,weiss_for_rot,paw_dmft%dmftqmc_l,fw1_nd,&
&         levels_ctqmc,levels_ctqmc_nd,hybri_limit,paw_dmft%temp,umod,opt_diag,opt_fk)
       end if

       ! Calculation of Inverse Green's function from hybridization
       !-------------------------------------------------------------
       do if1=1,2
         do if2=1,2
           do ifreq=1,paw_dmft%dmftqmc_l
             omega=pi*paw_dmft%temp*(two*float(ifreq)-1)
             if(if1==if2) then
               gw1_nd(ifreq,if1,if2) =  (cmplx(0.d0,omega,kind=dp)-fw1_nd(ifreq,if1,if2))
             else
               gw1_nd(ifreq,if1,if2) =  (-fw1_nd(ifreq,if1,if2))
             end if
           end do
         end do
       end do
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
       if(opt_fk==0) then
         fw1_nd=gw1_nd
       end if

       ABI_DEALLOCATE(gw1_nd)
       paw_dmft%dmft_nwlo=dmft_nwlo

     ! and testcode>1
     end if


     call flush_unit(std_out)
! =================================================================

! ___________________________________________________________________________________
!
!  THIRD PART : CALL CTQMC
! ___________________________________________________________________________________

! =================================================================
!    Main calls to CTQMC code in ABINIT (INITIALIZATION and OPTIONS)
! =================================================================
     if(paw_dmft%dmft_solv==5.or.paw_dmft%dmft_solv==8) then

       write(message,'(a,2x,a)') ch10,&
&       " == Initializing CTQMC"
       call wrtout(std_out,message,'COLL')

!    Initialisation
! =================================================================
     if(paw_dmft%dmft_solv==5) then
       nomega=paw_dmft%dmftqmc_l
       call CtqmcInterface_init(hybrid,paw_dmft%dmftqmc_seed,paw_dmft%dmftqmc_n, &
&       paw_dmft%dmftqmc_therm, paw_dmft%dmftctqmc_meas,nflavor,paw_dmft%dmftqmc_l,one/paw_dmft%temp,zero,&
&       std_out,paw_dmft%spacecomm)
!    options
! =================================================================
       call CtqmcInterface_setOpts(hybrid,&
       opt_Fk      =opt_fk,&
&       opt_order   =paw_dmft%dmftctqmc_order ,&
&       opt_movie   =paw_dmft%dmftctqmc_mov   ,&
&       opt_analysis=paw_dmft%dmftctqmc_correl,&
&       opt_check   =paw_dmft%dmftctqmc_check ,&
&       opt_noise   =paw_dmft%dmftctqmc_grnns ,&
&       opt_spectra =paw_dmft%dmftctqmc_mrka  ,&
&       opt_gmove   =paw_dmft%dmftctqmc_gmove )
     endif

     if(paw_dmft%dmft_solv==8) then
       nomega=paw_dmft%dmftqmc_l
       call CtqmcoffdiagInterface_init(hybridoffdiag,paw_dmft%dmftqmc_seed,paw_dmft%dmftqmc_n, &
&        paw_dmft%dmftqmc_therm, paw_dmft%dmftctqmc_meas,nflavor,&
&        paw_dmft%dmftqmc_l,one/paw_dmft%temp,zero,&
&        std_out,paw_dmft%spacecomm,opt_nondiag)
!    options
! =================================================================
       call CtqmcoffdiagInterface_setOpts(hybridoffdiag,&
       opt_Fk      =opt_fk,&
&       opt_order   =paw_dmft%dmftctqmc_order ,&
&       opt_movie   =paw_dmft%dmftctqmc_mov   ,&
&       opt_analysis=paw_dmft%dmftctqmc_correl,&
&       opt_check   =paw_dmft%dmftctqmc_check ,&
&       opt_noise   =paw_dmft%dmftctqmc_grnns ,&
&       opt_spectra =paw_dmft%dmftctqmc_mrka  ,&
&       opt_gmove   =paw_dmft%dmftctqmc_gmove )
     endif

       write(message,'(a,2x,2a)') ch10, " == Initialization CTQMC done", ch10
       call wrtout(std_out,message,'COLL')
     end if

     if(paw_dmft%dmft_solv==6.or.paw_dmft%dmft_solv==7) then
       !because size allocation problem with TRIQS paw_dmft%dmft_nwlo must be >= paw_dmft%dmft_nwli
       ABI_ALLOCATE(gw_tmp_nd,(paw_dmft%dmft_nwli,nflavor,nflavor))
       open(unit=505,file=trim(paw_dmft%filapp)//"_Legendre_coefficients.dat", status='unknown',form='formatted')
     else
       if(paw_dmft%dmft_solv==5) then
         ABI_ALLOCATE(gw_tmp,(paw_dmft%dmft_nwlo,nflavor+1))
       end if
       ABI_ALLOCATE(gw_tmp_nd,(paw_dmft%dmft_nwlo,nflavor,nflavor+1))
       !use  gw_tmp to put freq
       do ifreq=1,paw_dmft%dmft_nwlo
         if(paw_dmft%dmft_solv==5) gw_tmp(ifreq,nflavor+1)=cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)
         gw_tmp_nd(ifreq,nflavor,nflavor+1)=cmplx(zero,paw_dmft%omega_lo(ifreq),kind=dp)
       end do
     end if

     ABI_ALLOCATE(gtmp,(paw_dmft%dmftqmc_l,nflavor))
     ! THIS IS A BACKUP PLAN. USING paw_dmft%hybrid makes a segfault on TIKAL
     ! PSC with MPI only (and max2_open64). paw_dmf%hybrid is corrupted
     ! somewhere but I could not find the place in all DMFT routines
     ABI_ALLOCATE(gtmp_nd,(paw_dmft%dmftqmc_l,nflavor,nflavor))
     call flush_unit(std_out)

     ! =================================================================
     !    BEGIN CALL TO CTQMC SOLVERS
     ! =================================================================

     if(testcode==0) then

       ! =================================================================
       !    CTQMC run Abinit
       ! =================================================================
       if(paw_dmft%dmft_solv==5) then

         ABI_ALLOCATE(docc,(1:nflavor,1:nflavor))
         docc(:,:) = zero
         call CtqmcInterface_run(hybrid,fw1(1:paw_dmft%dmftqmc_l,:),Gtau=gtmp,&
&         Gw=gw_tmp,D=docc(:,:),E=green%ecorr_qmc(iatom),&
!&       matU=hu(itypat)%udens,opt_levels=levels_ctqmc)
&         matU=udens_atoms(iatom)%value,opt_levels=levels_ctqmc)
         call data4entropyDMFT_setDocc(paw_dmft%forentropyDMFT,iatom,docc)
         ABI_DEALLOCATE(docc)
         !DO iflavor = 1, nflavor
         !  hybrid%Hybrid%Greens(iflavor)%oper(1:this%samples) = gtmp(1:this%samples,iflavor)
         !  CALL GreenHyb_forFourier(this%Greens(iflavor), Gomega=Gw(:,iflavor), omega=Gw(:,this%flavors+1))
         !END DO

       ! =================================================================
       !    CTQMC run Abinit off diagonal terms in hybridization
       ! =================================================================
       else if (paw_dmft%dmft_solv==8) then
       ! =================================================================

         ABI_ALLOCATE(docc,(1:nflavor,1:nflavor))
         docc(:,:) = zero
         call CtqmcoffdiagInterface_run(hybridoffdiag,fw1_nd(1:paw_dmft%dmftqmc_l,:,:),Gtau=gtmp_nd,&
&        Gw=gw_tmp_nd,D=doccsum,E=green%ecorr_qmc(iatom),&
&        Noise=noise,matU=udens_atoms(iatom)%value,Docc=docc,opt_levels=levels_ctqmc,hybri_limit=hybri_limit)
         ABI_DEALLOCATE(docc)
       ! TODO: Handle de luj0 case for entropy

       ! =================================================================
       !    CTQMC run TRIQS
       ! =================================================================
       else if (paw_dmft%dmft_solv>=6.and.paw_dmft%dmft_solv<=7) then
       ! =================================================================

         if ( paw_dmft%dmftqmc_l >= (2*paw_dmft%dmft_nwli)+1 ) then

           call ctqmc_calltriqs(paw_dmft,cryst_struc,hu,levels_ctqmc,gtmp_nd,gw_tmp_nd,fw1_nd,leg_measure,iatom)

         else
           write(message,'(2a)') ch10," Can't launch TRIQS CTHYB solver because dmftqmc_l must be >= 2*dmft_nwli + 1"
           MSG_ERROR(message)
         end if

       end if

     ! =================================================================
     !    CTQMC run for tests
     ! =================================================================
     else if (testcode>=1) then
       call CtqmcInterface_run(hybrid,fw1(1:nomega,:),Gtau=gtmp,&
&       Gw=gw_tmp,E=green%ecorr_qmc(iatom),&
&       matU=umod,opt_levels=levels_ctqmc)

      ! for non diagonal code
      !       call CtqmcInterface_run(hybrid,fw1_nd(1:nomega,:,:),Gtau=gtmp_nd,&
      !&       Gw=gw_tmp_nd,D=Doccsum,E=green%ecorr_qmc(iatom),&
      !&       Noise=Noise,matU=umod,opt_levels=levels_ctqmc,hybri_limit=hybri_limit)

      !  If test of the code is activated, and testrot =1 rotate back green's function   and stop the code.
      ! --------------------------------------------------------------------------------------------------
       if(testcode==1) then

         call testcode_ctqmc(paw_dmft%dmftqmc_l,fw1_nd,fw1,gtmp_nd,gw_tmp_nd,&
&         levels_ctqmc,hybri_limit,nflavor,2,paw_dmft%temp,testrot,testcode,umod)

         write(message,'(2a)') ch10,' testcode end of test calculation'
         MSG_ERROR(message)
       end if
       if(testcode==2) then
         write(message,'(2a)') ch10,' testcode 2 end of test calculation'
         MSG_ERROR(message)
       end if

     end if
     ! =================================================================
     !    END CALL TO CTQMC SOLVERS
     ! =================================================================


     ! Print green function is files directly from CTQMC
     ! --------------------------------------------------
     call ctqmcoutput_printgreen(cryst_struc,eigvectmatlu,pawang,paw_dmft,gtmp_nd,gw_tmp_nd,gtmp,gw_tmp,iatom)


     ! If the CTQMC code in ABINIT was used, then destroy it and deallocate arrays
     ! ----------------------------------------------------------------------------
     if(paw_dmft%dmft_solv<6.and.paw_dmft%dmft_solv>7) then
     !Nothing just hybrid var problem
     else
       write(message,'(a,2x,a)') ch10,&
&       " == Destroy CTQMC"
       call wrtout(std_out,message,'COLL')
       if(paw_dmft%dmft_solv==5) call CtqmcInterface_finalize(hybrid)
       if(paw_dmft%dmft_solv==8) call CtqmcoffdiagInterface_finalize(hybridoffdiag)
       write(message,'(a,2x,a)') ch10,&
&       " == Destroy CTQMC done"
       call wrtout(std_out,message,'COLL')
     end if
     ABI_DEALLOCATE(hybri_limit)
     ABI_DEALLOCATE(levels_ctqmc_nd)
     ABI_DEALLOCATE(levels_ctqmc)
     ABI_DEALLOCATE(fw1)
     ABI_DEALLOCATE(fw1_nd)

! ___________________________________________________________________________________
!
!  FOURTH PART : USE OUTPUT OF CTQMC AND THEN DO BACK ROTATION
! ___________________________________________________________________________________
!

     ! Put green's function values from CTQMC into green structure
     !-------------------------------------------------------------
     call ctqmcoutput_to_green(green,paw_dmft,gtmp_nd,gw_tmp_nd,gtmp,gw_tmp,iatom,leg_measure,opt_nondiag)

     ! Deallocate arrays for CTQMC
     !-----------------------------
     if(paw_dmft%dmft_solv<6) then
       ABI_DEALLOCATE(gw_tmp)
     endif
     ABI_DEALLOCATE(gw_tmp_nd)
     ABI_DEALLOCATE(gtmp)
     ABI_DEALLOCATE(gtmp_nd)



     ! Do Fourier transform if it was not done (ie if TRIQS is used without legendre measurement)
     !----------------------------------------------------------------------------------------------
     if(opt_nondiag==1) then  ! (As leg_measure is activated by defautl, this fourier is never done).
       if(paw_dmft%dmft_solv>=6.and..not.leg_measure.and.paw_dmft%dmft_solv<=7) then
         write(message,'(2a,i3,13x,a)') ch10,'   ===  Direct Fourier Transform t->w of Weiss Field'
         call wrtout(std_out,message,'COLL')
         call fourier_green(cryst_struc,green,paw_dmft,&
&         pawang,opt_ksloc=2,opt_tw=1)
       end if
     endif


   end if

 end do ! iatom
! =========================================================================================
!  End big loop over atoms to compute hybridization and do the CTQMC
! =========================================================================================


 if(paw_dmft%dmft_prgn==1) then
   call print_green('QMC_diag_notsym',green,1,paw_dmft,pawprtvol=1,opt_wt=2)
   call print_green('QMC_diag_notsym',green,1,paw_dmft,pawprtvol=1,opt_wt=1)
 end if
 !write(message,'(i3,4x,2e21.14)') 6,weiss_for_rot%oper(1)%matlu(1)%mat(1,1,1,1,1)
 !call wrtout(std_out,message,'COLL')  ! debug
! =================================================================
! Inverse Weiss, then
! Copy Weiss_for_rot into weiss and rotate back weiss to the original basis
! =================================================================

! ABI_ALLOCATE(shift,(natom))
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
! ABI_DEALLOCATE(shift)

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
 if(pawprtvol>=3) then
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print green function for small tau after CTQMC"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(green%oper_tau(1)%matlu,natom,1)  ! debug
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print green function for small freq after CTQMC"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(green%oper(1)%matlu,natom,1)  ! debug
 end if


!  === Compute rotated Occupations in green%occup_tau
 call occup_green_tau(green)

 if(pawprtvol>=3) then
!  === Compute non rotated Occupations in green%occup_tau
   write(message,'(a,2x,a,f13.5)') ch10," == Occupation from G(tau) in the ctqmc basis"
   call wrtout(std_out,message,'COLL')
   call print_matlu(green%occup_tau%matlu,natom,1)
 end if

 write(message,'(a,2x,a,f13.5)') ch10,&
& " == Rotate Green function to original basis "
 call wrtout(std_out,message,'COLL')
 !write(message,'(i3,4x,2e21.14)') 8,weiss%oper(1)%matlu(1)%mat(1,1,1,1,1)
 !call wrtout(std_out,message,'COLL')  ! debug

 ! Rotate oper_tau into Ylm basis and then Slm basis
 !-------------------------------------------------------------
 do itau=1,paw_dmft%dmftqmc_l
   if(opt_diag/=0) call rotate_matlu(green%oper_tau(itau)%matlu,eigvectmatlu,natom,3,0)
   if(useylm==1) call slm2ylm_matlu(green%oper_tau(itau)%matlu,natom,2,0)
 end do

 ! Rotate occup_tau into Ylm basis and then Slm basis

 ! Rotate occup_tau into Ylm basis and then Slm basis
 !-------------------------------------------------------------
 if(opt_diag/=0) call rotate_matlu(green%occup_tau%matlu,eigvectmatlu,natom,3,0)
 if(useylm==1) then
   write(message,'(a,2x,a,f13.5)') ch10," == Occupation from G(tau) in the Ylm basis"
   call wrtout(std_out,message,'COLL')
   call print_matlu(green%occup_tau%matlu,natom,1)
 end if

 if(useylm==1) call slm2ylm_matlu(green%occup_tau%matlu,natom,2,0)
 write(message,'(a,2x,a,f13.5)') ch10," == Occupation from G(tau) in the Slm basis"
 call wrtout(std_out,message,'COLL')
 call print_matlu(green%occup_tau%matlu,natom,1)

 ! Put Weiss off diagonal terms to zero because Green function will not have any offdiag terms
 !------------------------------------------------------------------------------
 !   (if opt_nondiag=0 ie dmft_solv=5)
 do ifreq=1,paw_dmft%dmft_nwlo
   if(opt_nondiag==0) call zero_matlu(weiss%oper(ifreq)%matlu,natom,onlynondiag=1)
 end do
 !    ( if opt_nondiag=0, then:
 !       As Green's function is diagonal, one suppress off diag  terms in Weiss, if any.
 !      (If off diag are non zero in the density matrix and thus in the Green's function,
 !       there is a warning in checkreal_matlu above).)

 ! Rotate Green's and Weiss functions into Ylm basis and then Slm basis
 !-------------------------------------------------------------
 do ifreq=1,paw_dmft%dmft_nwlo
   if(opt_diag/=0) call rotate_matlu(green%oper(ifreq)%matlu,eigvectmatlu,natom,3,0)
   if(useylm==1) call slm2ylm_matlu(green%oper(ifreq)%matlu,natom,2,0)
   if(opt_diag/=0) call rotate_matlu(weiss%oper(ifreq)%matlu,eigvectmatlu,natom,3,0)
   if(useylm==1) call slm2ylm_matlu(weiss%oper(ifreq)%matlu,natom,2,0)
 end do
 !write(message,'(i3,4x,2e21.14)') 10,weiss%oper(1)%matlu(1)%mat(1,1,1,1,1)
 !call wrtout(std_out,message,'COLL')  ! debug

 if(pawprtvol>=3) then
   write(message,'(a,2x,a,f13.5)') ch10,&                  ! debug
&  " == Print green function for small time after rotation (in the original basis)" ! debug
   call wrtout(std_out,message,'COLL')                  ! debug
   call print_matlu(green%oper_tau(1)%matlu,natom,1)  ! debug
   write(message,'(a,2x,a,f13.5)') ch10,&                  ! debug
&  " == Print green function for small freq after rotation (in the original basis)" ! debug
   call wrtout(std_out,message,'COLL')                  ! debug
   call print_matlu(green%oper(1)%matlu,natom,1)  ! debug
   !< HACK >
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print diagonalized weiss_for_rot function after rotation for small freq in the ctqmc basis"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(weiss_for_rot%oper(1)%matlu,natom,1)  ! debug
   !</ HACK >
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print weiss function for small freq in the original basis"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(weiss%oper(1)%matlu,natom,1)  ! debug

   do ifreq=1,paw_dmft%dmft_nwlo
     call sym_matlu(cryst_struc,weiss%oper(ifreq)%matlu,pawang,paw_dmft)
   end do
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print symetrized weiss function for small freq in the original basis"  ! debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(weiss%oper(1)%matlu,natom,1)  ! debug
 end if


 ABI_DATATYPE_ALLOCATE(matlu1,(natom))
 call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,matlu1)
 call copy_matlu(green%occup_tau%matlu,matlu1,natom)
 call sym_matlu(cryst_struc,matlu1,pawang,paw_dmft)

 write(message,'(a,2x,a,f13.5)') ch10," == Occupation from G(tau) in the original basis"
 call wrtout(std_out,message,'COLL')
 call print_matlu(green%occup_tau%matlu,natom,1)

 write(message,'(a,2x,a,f13.5)') ch10," == Symetrized occupations"
 call wrtout(std_out,message,'COLL')
 call print_matlu(matlu1,natom,1)

 call diff_matlu("CTQMC Occup","CTQMC Occup symetrized",green%occup_tau%matlu,matlu1,natom,0,tol4,ierr)
 call destroy_matlu(matlu1,natom)
 ABI_DATATYPE_DEALLOCATE(matlu1)

! =================================================================
! Symetrise green function G(tau) and G(ifreq) to recover symetry
! artificially broken by QMC
! =================================================================
 write(message,'(a,2x,a,f13.5)') ch10,&
& " == Symetrise green function after QMC "
 call wrtout(std_out,message,'COLL')
 do itau=1,paw_dmft%dmftqmc_l
   call sym_matlu(cryst_struc,green%oper_tau(itau)%matlu,pawang,paw_dmft)
 end do
 do ifreq=1,paw_dmft%dmft_nwlo
   call sym_matlu(cryst_struc,green%oper(ifreq)%matlu,pawang,paw_dmft)
 end do
 if(pawprtvol>=3) then
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print green function for small time after symetrisation"  !  debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(green%oper_tau(1)%matlu,natom,1)  ! debug
   write(message,'(a,2x,a,f13.5)') ch10,&  ! debug
&  " == Print green function for small freq after symetrisation"  !  debug
   call wrtout(std_out,message,'COLL')  ! debug
   call print_matlu(green%oper(1)%matlu,natom,1)  ! debug
 end if
 if(paw_dmft%dmft_prgn==1) then
   call print_green('QMC_sym',green,1,paw_dmft,pawprtvol=1,opt_wt=2)
   call print_green('QMC_sym',green,1,paw_dmft,pawprtvol=1,opt_wt=1)
 end if

!  === Compute Occupations  (Symetrized from oper_tau)
 call occup_green_tau(green)


!  === Print occupations
 call printocc_green(green,6,paw_dmft,3)

 call destroy_oper(energy_level)
 call destroy_matlu(dmat_diag,natom)
 ABI_DATATYPE_DEALLOCATE(dmat_diag)
 call destroy_matlu(identity,natom)
 ABI_DATATYPE_DEALLOCATE(identity)
 do iatom=1,cryst_struc%natom
   lpawu=paw_dmft%lpawu(iatom)
   if(lpawu/=-1) then
     do isppol=1,nsppol
       ABI_DEALLOCATE(eigvectmatlu(iatom,isppol)%value)
       !ABI_DEALLOCATE(udens_atoms(iatom))
     end do
     ABI_DEALLOCATE(udens_atoms(iatom)%value)
   end if
 end do
 ABI_DATATYPE_DEALLOCATE(udens_atoms)
 ABI_DATATYPE_DEALLOCATE(eigvectmatlu)
 call destroy_green(weiss_for_rot)
! call destroy_green(gw_loc)
! call destroy_green(greenlda)

!  destroy limit of hybridization
 call destroy_matlu(hybri_coeff,paw_dmft%natom)
 ABI_DATATYPE_DEALLOCATE(hybri_coeff)

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
!! PARENTS
!!      qmc_prep_ctqmc
!!
!! CHILDREN
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
 complex(dpc), intent(out) :: fw1_nd(:,:,:)
 complex(dpc),  intent(inout) :: levels_ctqmc_nd(:,:)
 complex(dpc),  intent(inout) :: hybri_limit(:,:)
 type(oper_type)  :: energy_level
 type(matlu_type), allocatable  :: hybri_coeff(:)
 type(green_type)  :: weiss_for_rot

!Local variables ------------------------------
 integer :: ifreq,iatom,ima,imb,ispa,ispb
 real(dp) :: omega
 real(dp) :: facnd, facd
 character(len=30) :: tmpfil
! ************************************************************************
 facnd=0.8d0
 facd=1.0d0
 !write(6,*) "fac",facnd,facd
 levels_ctqmc_nd(2,2)   = energy_level%matlu(iatom)%mat(imb,imb,1,ispb,ispb)
 levels_ctqmc_nd(1,1)   = energy_level%matlu(iatom)%mat(ima,ima,1,ispa,ispa)
 levels_ctqmc(2)   = real(energy_level%matlu(iatom)%mat(imb,imb,1,ispb,ispb),kind=dp)
 levels_ctqmc(1)   = real(energy_level%matlu(iatom)%mat(ima,ima,1,ispa,ispa),kind=dp)
 if(opt_diag/=1) then
   levels_ctqmc_nd(1,2)   = energy_level%matlu(iatom)%mat(ima,imb,1,ispa,ispb)
   levels_ctqmc_nd(2,1)   = energy_level%matlu(iatom)%mat(imb,ima,1,ispb,ispa)
 end if
 hybri_limit(1,1)  = facd*hybri_coeff(iatom)%mat(ima,ima,1,ispa,ispa)
 hybri_limit(2,2)  = facd*hybri_coeff(iatom)%mat(imb,imb,1,ispb,ispb)
 hybri_limit(1,2)  = facnd*hybri_coeff(iatom)%mat(ima,imb,1,ispa,ispb)
 hybri_limit(2,1)  = facnd*hybri_coeff(iatom)%mat(imb,ima,1,ispb,ispa)
 !write(6,*) "hybri_limit",hybri_limit
 !write(6,*) "levels_ctqmc",levels_ctqmc
 umod=zero

 tmpfil = 'fw1_nd_re'
 !if (open_file(newunit=unt,message,file=trim(tmpfil),status='unknown',form='formatted')/=0) then
 !  MSG_ERROR(message)
 !end if
 tmpfil = 'fw1_nd_im'
 !if (open_file(newunit=unt2,message,file=trim(tmpfil),status='unknown',form='formatted')/=0) then
 !  MSG_ERROR(message)
 !end if
 write(std_out,*) "testcode==2",ispa,ispb,ima,imb
 write(std_out,*) "opt_fk==",opt_fk
 do ifreq=1,dmftqmc_l
   if (opt_fk==1) then
     fw1_nd(ifreq,1,1) = facd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(ima,ima,1,ispa,ispa)
     fw1_nd(ifreq,2,2) = facd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(imb,imb,1,ispb,ispb)
     !fw1_nd(ifreq,1,2) =  weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(ima,imb,1,ispa,ispb)
     !fw1_nd(ifreq,2,1) =  weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(imb,ima,1,ispb,ispa)
     fw1_nd(ifreq,1,2) = facnd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(ima,imb,1,ispa,ispb)
     fw1_nd(ifreq,2,1) = facnd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(imb,ima,1,ispb,ispa)
     omega=pi*temp*(two*float(ifreq)-1)
   else if (opt_fk==0) then
     fw1_nd(ifreq,1,1) =  facd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(ima,ima,1,ispa,ispa)
     fw1_nd(ifreq,2,2) =  facd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(imb,imb,1,ispb,ispb)
     fw1_nd(ifreq,1,2) =  facnd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(ima,imb,1,ispa,ispb)
     fw1_nd(ifreq,2,1) =  facnd*weiss_for_rot%oper(ifreq)%matlu(iatom)%mat(imb,ima,1,ispb,ispa)
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
!! PARENTS
!!      qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE

subroutine testcode_ctqmc(dmftqmc_l,fw1_nd,fw1,gtmp_nd,gw_tmp_nd,levels_ctqmc,hybri_limit,&
&   nflavor,opt,temp,testrot,testcode,umod)


!Arguments ------------------------------------
!scalars
 integer, intent(in) :: dmftqmc_l,nflavor,testrot,testcode,opt
 real(dp), intent(in) :: temp
 real(dp), intent(out) :: umod(2,2)
 complex(dpc), intent(inout) :: gw_tmp_nd(:,:,:)
 real(dp),  intent(inout) :: gtmp_nd(:,:,:)
 complex(dpc), intent(out) :: fw1(:,:)
 complex(dpc), intent(out) :: fw1_nd(:,:,:)
 real(dp),  intent(inout) :: levels_ctqmc(:)
 complex(dpc),  intent(inout) :: hybri_limit(:,:)

!Local variables ------------------------------
 character(len=500) :: message
 integer :: ifreq, itau,realrot,simplehyb
 real(dp) :: omega
 real(dp) :: tbi1,tbi2,e2,tbi3,tbi4,e3,e4,tbi21,tbi12,e3b,e4b,tbi21b,tbi12b
 complex(dpc) :: e1
! arrays
 complex(dpc) :: RR(2,2)
 complex(dpc) :: RR1(2,2)
 complex(dpc) :: RRi(2,2)
 complex(dpc) :: RRt(2,2)
! ************************************************************************
 if (testcode==0) return
 if (nflavor/=2) then
   write(message,'(2a)') ch10,' testcode nflavor.ne.2'
   MSG_ERROR(message)
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
   MSG_ERROR(message)
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
     MSG_WARNING(message)
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
   MSG_ERROR(message)

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
!! PARENTS
!!      qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE

subroutine ctqmcoutput_to_green(green,paw_dmft,gtmp_nd,gw_tmp_nd,gtmp,gw_tmp,iatom,leg_measure,opt_nondiag)

!Arguments ------------------------------------
!scalars
 type(paw_dmft_type), intent(in)  :: paw_dmft
 type(green_type), intent(inout) :: green
 real(dp), allocatable, intent(in) :: gtmp_nd(:,:,:)
 complex(dpc), allocatable, intent(in) :: gw_tmp(:,:)
 complex(dpc), allocatable, intent(in) :: gw_tmp_nd(:,:,:)
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
   green%oper_tau(itau)%matlu(iatom)%mat(:,:,:,:,:)=czero
 end do
 green%occup_tau%matlu(iatom)%mat(nflavor:,:,:,:,:)=czero

 do ifreq=1,paw_dmft%dmft_nwlo
   green%oper(ifreq)%matlu(iatom)%mat(:,:,:,:,:)=czero
 end do
 green%occup%matlu(iatom)%mat(:,:,:,:,:)=czero

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
               green%oper_tau(itau)%matlu(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)=&
&               gtmp_nd(itau,iflavor1,iflavor2)
               ! symetrize over spin if nsppol=nspinor=1
               if(paw_dmft%nsppol==1.and.paw_dmft%nspinor==1) then
                 green%oper_tau(itau)%matlu(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)=&
&                 (gtmp_nd(itau,iflavor1,iflavor2)+gtmp_nd(itau,iflavor1+tndim,iflavor2+tndim))/two
               end if
             end do  !itau
             if(paw_dmft%dmft_solv<6.or.leg_measure) then
               do ifreq=1,paw_dmft%dmft_nwlo
                 green%oper(ifreq)%matlu(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)=&
&                 gw_tmp_nd(ifreq,iflavor1,iflavor2)
               ! symetrize over spin if nsppol=nspinor=1
                 if(paw_dmft%nsppol==1.and.paw_dmft%nspinor==1) then
                   green%oper(ifreq)%matlu(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)=&
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
           green%oper_tau(itau)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)=gtmp(itau,iflavor)
           ! symetrize over spin if nsppol=paw_dmft%nspinor=1
           if(paw_dmft%nsppol==1.and.paw_dmft%nspinor==1) then
             green%oper_tau(itau)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)=&
&             (gtmp(itau,iflavor)+gtmp(itau,iflavor+tndim))/two
           end if
         end do
!       ifreq2=0
         do ifreq=1,paw_dmft%dmft_nwlo
!         if(paw_dmft%select_log(ifreq)==1) then
!           ifreq2=ifreq2+1
           green%oper(ifreq)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)=gw_tmp(ifreq,iflavor)
           ! symetrize over spin if nsppol=paw_dmft%nspinor=1
           if(paw_dmft%nsppol==1.and.paw_dmft%nspinor==1) then
             green%oper(ifreq)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)=&
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
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  iatom = atoms on which the calculation has been done
!!
!! OUTPUT
!!
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE

subroutine ctqmcoutput_printgreen(cryst_struc,eigvectmatlu,pawang,paw_dmft,gtmp_nd,gw_tmp_nd,gtmp,gw_tmp,iatom)

!Arguments ------------------------------------
!scalars
 type(paw_dmft_type), intent(in)  :: paw_dmft
 real(dp), allocatable, intent(inout) :: gtmp_nd(:,:,:)
 complex(dpc), allocatable, intent(in) :: gw_tmp(:,:)
 complex(dpc), allocatable, intent(in) :: gw_tmp_nd(:,:,:)
 real(dp), allocatable, intent(in) :: gtmp(:,:)
 type(crystal_t),intent(in) :: cryst_struc
 integer, intent(in) :: iatom
 type(coeff2c_type), intent(inout) :: eigvectmatlu(:,:)
 type(pawang_type), intent(in) :: pawang

!Local variables ------------------------------
 character(len=500) :: message
 type(matlu_type), allocatable :: matlu1(:)
 integer :: ifreq, itau,im1,im2,isppol,ispinor1,ispinor2,iflavor1
 integer :: iflavor2,tndim,iflavor,nflavor
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
      MSG_ERROR(message)
    end if
    do ifreq=1,paw_dmft%dmft_nwli
      write(unt,'(29f21.14)') paw_dmft%omega_lo(ifreq),((gw_tmp_nd(ifreq,iflavor,iflavor)), iflavor=1, nflavor)
    end do
    close(unt)
  else
    if(paw_dmft%dmft_solv==5) then
      if (open_file(trim(paw_dmft%filapp)//"_atom_"//iatomnb//"_Gtau_"//gtau_iter//".dat", message, newunit=unt) /= 0) then
        MSG_ERROR(message)
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
        MSG_ERROR(message)
      end if
      do itau=1,paw_dmft%dmftqmc_l
        write(unt,'(196f21.14)') float(itau-1)/float(paw_dmft%dmftqmc_l)/paw_dmft%temp,&
        ((gtmp_nd(itau,iflavor,iflavor1), iflavor=1, nflavor),iflavor1=1, nflavor)
      end do
      close(unt)
      if(paw_dmft%natom==1) then ! If natom>1, it should be moved outside the loop over atoms
        ABI_DATATYPE_ALLOCATE(matlu1,(paw_dmft%natom))
        call init_matlu(paw_dmft%natom,paw_dmft%nspinor,paw_dmft%nsppol,paw_dmft%lpawu,matlu1)
        do itau=1,paw_dmft%dmftqmc_l
          do isppol=1,paw_dmft%nsppol
            do ispinor1=1,paw_dmft%nspinor
              do im1=1,tndim
                iflavor1=im1+tndim*(ispinor1-1)+tndim*(isppol-1)
                do ispinor2=1,paw_dmft%nspinor
                  do im2=1,tndim
                    iflavor2=im2+tndim*(ispinor2-1)+tndim*(isppol-1)
                    matlu1(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)=&
&                     gtmp_nd(itau,iflavor1,iflavor2)
                  end do  ! im2
                end do  ! ispinor2
              end do  ! im1
            end do  ! ispinor
          end do ! isppol
          call rotate_matlu(matlu1,eigvectmatlu,paw_dmft%natom,3,0)
          call slm2ylm_matlu(matlu1,paw_dmft%natom,2,0)
          call sym_matlu(cryst_struc,matlu1,pawang,paw_dmft)
          call slm2ylm_matlu(matlu1,paw_dmft%natom,1,0)
          call rotate_matlu(matlu1,eigvectmatlu,paw_dmft%natom,3,1)
          do isppol=1,paw_dmft%nsppol
            do ispinor1=1,paw_dmft%nspinor
              do im1=1,tndim
                iflavor1=im1+tndim*(ispinor1-1)+tndim*(isppol-1)
                do ispinor2=1,paw_dmft%nspinor
                  do im2=1,tndim
                    iflavor2=im2+tndim*(ispinor2-1)+tndim*(isppol-1)
                    gtmp_nd(itau,iflavor1,iflavor2)=&
                     matlu1(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)
                  end do  ! im2
                end do  ! ispinor2
              end do  ! im1
            end do  ! ispinor
          end do ! isppol
        end do  !itau
        call destroy_matlu(matlu1,paw_dmft%natom)
        ABI_DATATYPE_DEALLOCATE(matlu1)
      endif ! if natom=1
      if (open_file(trim(paw_dmft%filapp)//"_atom_"//iatomnb//"_Gtau_offdiag_"//gtau_iter//".dat",&
&      message, newunit=unt) /= 0) then
        MSG_ERROR(message)
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
        MSG_ERROR(message)
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
!! PARENTS
!!      qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE

subroutine ctqmc_calltriqs(paw_dmft,cryst_struc,hu,levels_ctqmc,gtmp_nd,gw_tmp_nd,fw1_nd,leg_measure,iatom)

#if defined HAVE_TRIQS_v2_0 || defined HAVE_TRIQS_v1_4
 use TRIQS_CTQMC !Triqs module
#endif
 use ISO_C_BINDING

!Arguments ------------------------------------
!scalars
 type(paw_dmft_type), intent(in)  :: paw_dmft
 type(crystal_t),intent(in) :: cryst_struc
 type(hu_type), intent(in) :: hu(cryst_struc%ntypat)
 real(dp), allocatable, target, intent(inout) :: gtmp_nd(:,:,:)
 complex(dpc), allocatable, target, intent(inout) :: gw_tmp_nd(:,:,:)
 complex(dpc), allocatable, target, intent(in) :: fw1_nd(:,:,:)
 real(dp), allocatable, target, intent(inout) ::  levels_ctqmc(:)
 logical(kind=1), intent(in) :: leg_measure
 integer, intent(in) :: iatom

!Local variables ------------------------------
 complex(dpc), allocatable, target ::fw1_nd_tmp(:,:,:)
 complex(dpc), allocatable, target :: g_iw(:,:,:)
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
#if defined HAVE_TRIQS_v2_0 || defined HAVE_TRIQS_v1_4
 logical(kind=1) :: hist = .false.
 logical(kind=1) :: wrt_files = .true.
 logical(kind=1) :: tot_not = .true.
#endif
 real(dp) :: beta,besp,bespp,xx
 complex(dpc) :: u_nl
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
 nleg  = paw_dmft%dmftctqmc_triqs_nleg
 nflavor=2*(2*paw_dmft%lpawu(iatom)+1)
 itypat=cryst_struc%typat(iatom)


 verbosity_solver = paw_dmft%prtvol
 beta = 1.0/(paw_dmft%temp*Ha_eV)

 !Allocation in/output array phase:
 ABI_ALLOCATE(fw1_nd_tmp,(1:nflavor,1:nflavor,1:nfreq)) !column major
 ABI_ALLOCATE(g_iw,(1:nflavor,1:nflavor,1:nfreq)) !column major
 ABI_ALLOCATE(u_mat_ij,(1:nflavor,1:nflavor)) !column major
 ABI_ALLOCATE(u_mat_ijkl,(1:nflavor,1:nflavor,1:nflavor,1:nflavor)) !column major
 ABI_ALLOCATE(u_mat_ijkl_tmp,(1:nflavor,1:nflavor,1:nflavor,1:nflavor)) !column major

 if ( leg_measure ) then !only if functionality is enabled
   ABI_ALLOCATE(gl_nd,(1:nleg,1:nflavor,1:nflavor)) !column major !nl = 30 by default
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

 call vee_ndim2tndim_hu_r(paw_dmft%lpawu(iatom),hu(itypat)%vee,u_mat_ijkl_tmp,1)
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
 u_mat_ij     = transpose( hu(itypat)%udens ) * Ha_eV !column -> row major + conversion
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

  !Calling interfaced TRIQS solver subroutine from src/01_triqs_ext package
  !----------------------------------------------------
#if defined HAVE_TRIQS_v2_0 || defined HAVE_TRIQS_v1_4
 call Ctqmc_triqs_run (     rot_inv, leg_measure, hist, wrt_files, tot_not,                            &
&  nflavor, nfreq, ntau , nleg, int(paw_dmft%dmftqmc_n/paw_dmft%nproc),       &
&  paw_dmft%dmftctqmc_meas*2*2*nflavor, paw_dmft%dmftqmc_therm,               &
&  verbosity_solver, paw_dmft%dmftqmc_seed,beta,                              &
&  levels_ptr,  u_mat_ij_ptr, u_mat_ijkl_ptr, fw1_nd_ptr,                     &
&  g_iw_ptr, gtau_ptr, gl_ptr, paw_dmft%spacecomm                             )
#endif

  !WRITE(*,*) "Hello Debug"
  !call xmpi_barrier(paw_dmft%spacecomm) !Resynch all processus after calling Impurity solver from TRIQS

  !Report output datas from TRIQS to Abinit
  !Interacting G(iw)
 do ifreq=1,nfreq
   do iflavor1=1,nflavor
     do iflavor=1,nflavor
    !   gw_tmp_nd(ifreq,iflavor,iflavor1) = g_iw(iflavor,iflavor1,ifreq) !* Ha_eV !because 1/ G0(eV)
    !  WRITE(503,*) "[OUT Fortran] G(iw)[ w= ",ifreq," l= ",iflavor," l_= ",iflavor1,"] = ",gw_tmp_nd(ifreq,iflavor,iflavor1)!g_iw(iflavor,iflavor1,ifreq)
     end do
   end do
 end do

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
 ABI_DEALLOCATE( fw1_nd_tmp )
 ABI_DEALLOCATE( g_iw )
 ABI_DEALLOCATE( u_mat_ijkl )
 ABI_DEALLOCATE( u_mat_ijkl_tmp )
 ABI_DEALLOCATE( u_mat_ij )


  !  Compute Green's function in imaginary freq using Legendre coefficients
  ! -----------------------------------------------------------------------
 if (leg_measure) then
   call xmpi_barrier(paw_dmft%spacecomm)
   call flush_unit(std_out)
   write(message,'(2a)') ch10,"    ==  Compute G(iw_n) from Legendre coefficients"
   call wrtout(std_out,message,'COLL')
   ABI_ALLOCATE( jbes, (nleg))
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
   ABI_DEALLOCATE( jbes )
   call xmpi_barrier(paw_dmft%spacecomm)
   call flush_unit(std_out)
 end if
 gw_tmp_nd = gw_tmp_nd*Ha_eV


 if ( leg_measure ) then !only if functionality is enabled
   ABI_DEALLOCATE(gl_nd)
 end if


end subroutine ctqmc_calltriqs
!!***

END MODULE m_forctqmc
!!***
