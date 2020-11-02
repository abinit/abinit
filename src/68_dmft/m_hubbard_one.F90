!!****m* ABINIT/m_hubbard_one
!! NAME
!!  m_hubbard_one
!!
!! FUNCTION
!!
!! Solve Anderson model with the density/density Hubbard one approximation
!!
!! COPYRIGHT
!! Copyright (C) 2006-2020 ABINIT group (BAmadon)
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

MODULE m_hubbard_one


 use defs_basis

 implicit none

 private

 public :: hubbard_one
!!***

contains

!!****f* m_hubbard_one/hubbard_one
!! NAME
!! hubbard_one
!!
!! FUNCTION
!! Solve the hubbard one approximation
!!
!! COPYRIGHT
!! Copyright (C) 1999-2020 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cryst_struc
!!  istep    =  step of iteration for DFT.
!!  dft_occup
!!  mpi_enreg=informations about MPI parallelization
!!  paw_dmft =  data for self-consistent DFT+DMFT calculations.
!!  pawang <type(pawang)>=paw angular mesh and related data
!!
!! OUTPUT
!!  paw_dmft =  data for self-consistent DFT+DMFT calculations.
!!
!! NOTES
!!
!! PARENTS
!!      m_dmft
!!
!! CHILDREN
!!      combin,destroy_green,init_green,wrtout
!!
!! SOURCE

subroutine hubbard_one(cryst_struc,green,hu,paw_dmft,pawang,pawprtvol,hdc,weiss)

 use defs_basis
 use m_errors
 use m_abicore

 use m_pawang, only : pawang_type
 use m_crystal, only : crystal_t
 use m_green, only : green_type,init_green,destroy_green
 use m_paw_dmft, only : paw_dmft_type
 use m_oper, only : oper_type,init_oper,destroy_oper,loc_oper,print_oper
 use m_matlu, only : matlu_type,sym_matlu, print_matlu, gather_matlu,&
& diag_matlu,init_matlu,destroy_matlu,rotate_matlu,copy_matlu,slm2ylm_matlu
 use m_hu, only : hu_type,rotatevee_hu
 use m_datafordmft, only : compute_levels
 implicit none

!Arguments ------------------------------------
!scalars
! type(pawang_type), intent(in) :: pawang
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type), intent(inout) :: green
 type(paw_dmft_type), intent(in)  :: paw_dmft
 type(hu_type), intent(inout) :: hu(cryst_struc%ntypat)
 type(pawang_type), intent(in) :: pawang
 type(oper_type), intent(inout) :: hdc
 integer, intent(in) :: pawprtvol
 type(green_type), intent(inout) :: weiss

!Local variables ------------------------------
 type  :: level2_type
  integer, pointer :: repart(:,:) => null()
  integer, pointer :: ocp(:,:) => null()
  integer, pointer :: transition(:,:) => null()
  integer, pointer :: transition_m(:,:) => null()
 end type level2_type
 type  :: level1_type
  real(dp), pointer :: config(:) => null()
 end type level1_type
! scalars
 character(len=500) :: message
 integer :: iatom,ifreq,im,im1,isppol,ispinor,ispinor1
 integer :: lpawu,mbandc,natom,nkpt,nspinor,nsppol,nsppol_imp,testblock,tndim,useylm
! complex(dpc) :: g,g0,w
! arrays
 complex(dp), allocatable :: Id(:,:,:,:)
 type(coeff2c_type), allocatable :: eigvectmatlu(:,:)
 type(coeff2_type), allocatable :: udens_atoms(:)
 type(oper_type)  :: energy_level
 type(green_type) :: green_hubbard
 type(matlu_type), allocatable :: level_diag(:)
 complex(dpc) :: omega_current
! ************************************************************************
 mbandc=paw_dmft%mbandc
 nkpt=paw_dmft%nkpt
 nsppol=paw_dmft%nsppol
 natom=paw_dmft%natom
 nspinor=paw_dmft%nspinor


!Initialise for compiler
 omega_current=czero

!======================================
!Allocations: levels and eigenvectors
!======================================
 ABI_DATATYPE_ALLOCATE(level_diag,(natom))
 ABI_DATATYPE_ALLOCATE(eigvectmatlu,(natom,nsppol))
 ABI_DATATYPE_ALLOCATE(udens_atoms,(natom))
 call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,level_diag)
 do iatom=1,cryst_struc%natom
   lpawu=paw_dmft%lpawu(iatom)
   if(lpawu/=-1) then
     tndim=nspinor*(2*lpawu+1)
     do isppol=1,nsppol
       ABI_ALLOCATE(eigvectmatlu(iatom,isppol)%value,(tndim,tndim))
     end do
     ABI_ALLOCATE(udens_atoms(iatom)%value,(2*(2*lpawu+1),2*(2*lpawu+1)))
     level_diag(iatom)%mat=czero
   end if
 end do

 call init_oper(paw_dmft,energy_level,opt_ksloc=3)
 call compute_levels(cryst_struc,energy_level,hdc,pawang,paw_dmft)
!!========================
!!Get KS eigenvalues
!!========================
! call init_oper(paw_dmft,energy_level,opt_ksloc=3)
! do iband=1,mbandc
!   do ikpt=1,nkpt
!     do isppol=1,nsppol
!!      Take \epsilon_{nks}
!!      ========================
!       energy_level%ks(isppol,ikpt,iband,iband)=paw_dmft%eigen_dft(isppol,ikpt,iband)
!     end do
!   end do
! end do
!
!
!!======================================================================
!!Compute atomic levels from projection of \epsilon_{nks} and symetrize
!!======================================================================
! call loc_oper(energy_level,paw_dmft,1)
! write(message,'(a,2x,a,f13.5)') ch10," == Print Energy levels before sym and only DFT"
! call wrtout(std_out,message,'COLL')
! call print_matlu(energy_level%matlu,natom,1)
! do iatom = 1 , natom
!   lpawu=paw_dmft%lpawu(iatom)
!   if(lpawu/=-1) then
!     do isppol=1,nsppol
!       do ispinor=1,nspinor
!         do im1=1,2*lpawu+1
!           energy_level%matlu(iatom)%mat(im1,im1,isppol,ispinor,ispinor)=&
!&           energy_level%matlu(iatom)%mat(im1,im1,isppol,ispinor,ispinor)&
!&           -hdc%matlu(iatom)%mat(im1,im1,isppol,ispinor,ispinor)-paw_dmft%fermie
!         end do
!       end do
!     end do
!!    write(std_out,*) "DC,fermie",hdc%matlu(iatom)%mat(1,1,1,1,1),paw_dmft%fermie
!   end if
! end do ! natom
! call sym_matlu(cryst_struc,energy_level%matlu,pawang)
!
! write(message,'(a,2x,a,f13.5)') ch10," == Print Energy levels for Fermi Level=",paw_dmft%fermie
! call wrtout(std_out,message,'COLL')
!!call print_oper(energy_level,1,paw_dmft,1)
! call print_matlu(energy_level%matlu,natom,1)

!========================
!Compute Weiss function
!========================
 ABI_ALLOCATE(Id,(20,20,nspinor,nspinor))
 do iatom = 1 , natom
   lpawu=paw_dmft%lpawu(iatom)
   if(lpawu/=-1) then
     Id=czero
     do im=1,2*lpawu+1
       do ispinor=1,nspinor
         Id(im,im,ispinor,ispinor)=cone
       end do
     end do ! ib
     do ifreq=1,weiss%nw
       if(weiss%w_type=="imag") then
         omega_current=cmplx(zero,weiss%omega(ifreq),kind=dp)
       else if(green%w_type=="real") then
         omega_current=cmplx(weiss%omega(ifreq),zero,kind=dp)
       end if
       do im=1,2*lpawu+1
         do im1=1,2*lpawu+1
           do isppol=1,nsppol
             do ispinor=1,nspinor
               do ispinor1=1,nspinor
                 weiss%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)=&
&                 ( omega_current*Id(im,im1,ispinor,ispinor1) - &
&                 energy_level%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
               end do ! ispinor1
             end do ! ispinor
           end do ! isppol
         end do ! im1
       end do ! im
     end do ! ifreq
   end if ! lpawu
 end do ! natom
 ABI_DEALLOCATE(Id)

!=================================================================
!Diagonalizes atomic levels and keep eigenvectors in eigvectmatlu
!=================================================================
!if jpawu=0, rotatevee_hu will have no effect so it is not necessary to
!have a single rotation matrix for up and dn spins.

 if(hu(1)%jpawu_zero.and.nsppol==2) nsppol_imp=2
 if(.not.hu(1)%jpawu_zero.or.nsppol/=2) nsppol_imp=1
!  Diagonalize energy levels
 useylm=paw_dmft%dmft_blockdiag
 if(useylm==1) call slm2ylm_matlu(energy_level%matlu,natom,1,pawprtvol)
 testblock=1
 if(useylm==1) testblock=8
 call diag_matlu(energy_level%matlu,level_diag,natom,&
& prtopt=pawprtvol,eigvectmatlu=eigvectmatlu,nsppol_imp=nsppol_imp,optreal=1,test=testblock)

!  Use rotation matrix to rotate interaction
 if(useylm==1) then 
   call rotatevee_hu(cryst_struc,hu,nspinor,nsppol,pawprtvol,eigvectmatlu,udens_atoms,4)
 else
   call rotatevee_hu(cryst_struc,hu,nspinor,nsppol,pawprtvol,eigvectmatlu,udens_atoms,1)
 endif 
!write(std_out,*)"udens after rotatevee", udens_atoms(1)%value
 write(message,'(a,2x,a,f13.5)') ch10,&
& " == Print Diagonalized Energy levels for Fermi Level=",paw_dmft%fermie
 call wrtout(std_out,message,'COLL')
 call print_matlu(level_diag,natom,1)

!  Print out
 if(nspinor==2) then
   write(message,'(a,2x,a,f13.5)') ch10,&
&   " == Print weiss for small freq"
   call wrtout(std_out,message,'COLL')
   call print_matlu(weiss%oper(1)%matlu,natom,1)
   write(message,'(a,2x,a,f13.5)') ch10,&
&   " == Print weiss for large freq"
   call wrtout(std_out,message,'COLL')
   call print_matlu(weiss%oper(weiss%nw)%matlu,natom,1)
 end if

!========================
!Compute Green function
!========================
 call init_green(green_hubbard,paw_dmft,opt_oper_ksloc=2,wtype=green%w_type) ! initialize only matlu
!write(std_out,*)"udens", udens_atoms(1)%value
! write(std_out,*)"levels",  level_diag(1)%mat
 call green_atomic_hubbard(cryst_struc,green_hubbard,hu,level_diag,paw_dmft,udens_atoms)
!call rotate_matlu(energy_level%matlu,natom,pawprtvol=3)
!write(81,*) "I1",paw_dmft%omega_lo(1), real(green%oper(1)%matlu(1)%mat(1,1,1,1,1)),imag(green%oper(1)%matlu(1)%mat(1,1,1,1,1))
!========================================================================
!Rotate back Green function in the original basis before diagonalization
!========================================================================
!call print_matlu(level_diag,natom,1)
!test scall rotate_matlu(level_diag,eigvectmatlu,natom,3)
!todo_ab: add check here for back rotation
!call print_matlu(level_diag,natom,1)
 write(message,'(2a,f13.5)') ch10," == Green function before rotation"
 call wrtout(std_out,message,'COLL')
 call print_matlu(green_hubbard%oper(1)%matlu,natom,1)
 do ifreq=1,green_hubbard%nw
   call rotate_matlu(green_hubbard%oper(ifreq)%matlu,eigvectmatlu,natom,3,0)
   if(useylm==1) call slm2ylm_matlu(green_hubbard%oper(ifreq)%matlu,natom,2,0)
   call copy_matlu(green_hubbard%oper(ifreq)%matlu,green%oper(ifreq)%matlu,natom)
 end do
 write(message,'(2a,f13.5)') ch10," == Green function after rotation"
 call wrtout(std_out,message,'COLL')
 call print_matlu(green%oper(1)%matlu,natom,1)
 if(nspinor==2) then
   write(message,'(a,2x,a,f13.5)') ch10,&
&   " == Print green for small freq"
   call wrtout(std_out,message,'COLL')
   call print_matlu(green%oper(1)%matlu,natom,1)
   write(message,'(a,2x,a,f13.5)') ch10,&
&   " == Print green for large freq"
   call wrtout(std_out,message,'COLL')
   call print_matlu(green%oper(green%nw)%matlu,natom,1)
 end if
!do ifreq=1,paw_dmft%dmft_nwlo
!g=green%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)
!g0=cone/weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)
!w=cmplx(0.d0,paw_dmft%omega_lo(ifreq),kind=dp)
!write(160,*) paw_dmft%omega_lo(ifreq),real(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),imag(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
!write(161,*) paw_dmft%omega_lo(ifreq),real(green%oper(ifreq)%matlu(1)%mat(1,1,1,1,1) ),imag(green%oper(ifreq)%matlu(1)%mat(1,1,1,1,1) )
!write(164,*) paw_dmft%omega_lo(ifreq),real(one/green%oper(ifreq)%matlu(1)%mat(1,1,1,1,1) ),imag(one/green%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
!write(166,*) paw_dmft%omega_lo(ifreq),real(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)-one/green%oper(ifreq)%matlu(1)%mat(1,1,1,1,1) ),imag(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)-one/green%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
!write(167,*) paw_dmft%omega_lo(ifreq),real((g-g0)/(g0*g)),imag((g-g0)/(g0*g))
!write(168,*) paw_dmft%omega_lo(ifreq),real(1/g-w),imag(1/g-w)
!write(169,*) paw_dmft%omega_lo(ifreq),real(1/g0-w),imag(1/g0-w)
!write(170,*) paw_dmft%omega_lo(ifreq),w
!write(171,*) paw_dmft%omega_lo(ifreq),real(1/g),imag(1/g)
!write(172,*) paw_dmft%omega_lo(ifreq),real(w),imag(w)

!! voir si en faisant GG0/(G-G0) cela reduit l'erreur
!enddo
!call abi_abort('COLL')


!write(message,'(2a,f13.5)') ch10," == Print Energy levels after diagonalisation"
!call wrtout(std_out,message,'COLL')
!call print_matlu(energy_level%matlu,natom,1)

!======================================
!Deallocations and destroys
!======================================
 call destroy_green(green_hubbard)
 call destroy_oper(energy_level)
 call destroy_matlu(level_diag,natom)
 ABI_DATATYPE_DEALLOCATE(level_diag)
 do iatom=1,cryst_struc%natom
   lpawu=paw_dmft%lpawu(iatom)
   if(lpawu/=-1) then
     ABI_DEALLOCATE(udens_atoms(iatom)%value)
     do isppol=1,nsppol
       ABI_DEALLOCATE(eigvectmatlu(iatom,isppol)%value)
     end do
   end if
 end do
 ABI_DATATYPE_DEALLOCATE(eigvectmatlu)
 ABI_DATATYPE_DEALLOCATE(udens_atoms)
!!***

contains

!!****f* m_hubbard_one/green_atomic_hubbard
!! NAME
!! green_atomic_hubbard
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!! Copyright (C) 1999-2020 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cryst_struc
!!  istep    =  step of iteration for DFT.
!!  dft_occup
!!  mpi_enreg=informations about MPI parallelization
!!  paw_dmft =  data for self-consistent DFT+DMFT calculations.
!!  pawang <type(pawang)>=paw angular mesh and related data
!!
!! OUTPUT
!!  paw_dmft =  data for self-consistent DFT+DMFT calculations.
!!
!! NOTES
!!
!! PARENTS
!!      m_hubbard_one
!!
!! CHILDREN
!!      combin,destroy_green,init_green,wrtout
!!
!! SOURCE

subroutine green_atomic_hubbard(cryst_struc,green_hubbard,hu,level_diag,paw_dmft,udens_atoms)


 use defs_basis
 use m_errors
 use m_abicore
 use m_crystal, only : crystal_t
 use m_special_funcs,  only : factorial, permutations
 use m_green, only : green_type,init_green,destroy_green
 use m_hu, only : hu_type
 use m_paw_dmft, only : paw_dmft_type
 implicit none

!Arguments ------------------------------------
!scalars
! type(pawang_type), intent(in) :: pawang
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type), intent(inout) :: green_hubbard
 type(paw_dmft_type), intent(in)  :: paw_dmft
 type(matlu_type),intent(in) :: level_diag(cryst_struc%natom)
 type(hu_type), intent(inout) :: hu(cryst_struc%ntypat)
 type(coeff2_type),intent(in) :: udens_atoms(cryst_struc%natom)

!Local variables ------------------------------
! scalars
 integer :: cnk,iacc,iatom,iconfig,ielec,ifreq,ilevel,im,im1,isppol,ispinor,itrans,jconfig,jelec
 integer :: lpawu,m_temp,nconfig,nelec,nlevels,nspinor,nsppol,occupied_level,sum_test
 integer, allocatable :: occup(:,:),nconfig_nelec(:)
 character(len=500) :: message
! arrays
 type(green_type) :: green_hubbard_realw
 type(level2_type), allocatable :: occ_level(:)
 type(level1_type), allocatable :: e_nelec(:)
 complex(dpc), allocatable :: green_temp(:,:)
 complex(dpc), allocatable :: green_temp_realw(:,:)
 complex(dpc) :: Z_part
 real(dp), allocatable :: maxener(:),minener(:)
 real(dp), allocatable :: elevels(:)
 real(dp) :: emax,emin,eshift,prtopt, Ej_np1, Ei_n,beta,maxarg_exp,tmp
!************************************************************************
   maxarg_exp=300

!  hu is not used anymore.
   if(hu(1)%lpawu==0) then
   end if
!  ======================================
!  General loop over atoms
!  ======================================
   nsppol=paw_dmft%nsppol
   nspinor=paw_dmft%nspinor
   prtopt=1
   beta=one/paw_dmft%temp
   call init_green(green_hubbard_realw,paw_dmft,opt_oper_ksloc=2,wtype=green_hubbard%w_type) ! initialize only matlu

   do iatom=1,cryst_struc%natom
     lpawu=paw_dmft%lpawu(iatom)
     if(lpawu/=-1) then
       nlevels=nsppol*nspinor*(2*lpawu+1)
       if(nsppol==1.and.nspinor==1) nlevels=2*nspinor*(2*lpawu+1)

!      ===================================
!      Allocations
!      ===================================
       ABI_DATATYPE_ALLOCATE(occ_level,(0:nlevels))
       ABI_ALLOCATE(maxener,(0:nlevels))
       ABI_ALLOCATE(minener,(0:nlevels))
       ABI_ALLOCATE(elevels,(nlevels))
       ABI_DATATYPE_ALLOCATE(e_nelec,(0:nlevels))
       do nelec=0,nlevels ! number of electrons
         cnk=nint(permutations(nlevels,nelec)/factorial(nelec))
         ABI_ALLOCATE(occ_level(nelec)%repart      ,(cnk,nelec))
         ABI_ALLOCATE(occ_level(nelec)%ocp         ,(cnk,nlevels))
         ABI_ALLOCATE(occ_level(nelec)%transition  ,(cnk,nlevels-nelec))
         ABI_ALLOCATE(occ_level(nelec)%transition_m,(cnk,nlevels))
         ABI_ALLOCATE(e_nelec  (nelec)%config      ,(cnk))
         e_nelec(nelec)%config(:)=zero
!        write(std_out,*) "permutations",nint(permutations(nlevels,nelec)/factorial(nelec))
!        write(std_out,*) "size",size(occ_level),size(occ_level(nelec)%repart,1)
!        write(std_out,*) "size",size(occ_level),size(occ_level(nelec)%repart,2)
!        for a given nb of electrons nelec, gives for a given repartition
!        of electron, the position of the ielec electron inside atomic
!        levels
!        levels
       end do
       ABI_ALLOCATE(occup,(0:nlevels,nlevels))
       ABI_ALLOCATE(nconfig_nelec,(0:nlevels))

!      ===================================
!      Initialization
!      ===================================
       nconfig_nelec=0
       nconfig=1
       occup=0
       nconfig_nelec(0)=1
       occup(0,:)=0
       iacc=0
       elevels=zero
       do isppol=1,nsppol
         do ispinor=1,nspinor
           do im1=1,(2*lpawu+1)
             iacc=iacc+1
             elevels(iacc)=level_diag(iatom)%mat(im1,im1,isppol,ispinor,ispinor)
             if(abs(aimag(level_diag(iatom)%mat(im1,im1,isppol,ispinor,ispinor)))>tol8) then
               message = " Hubbard I: levels are imaginary"
               write(std_out,*) level_diag(iatom)%mat(im1,im1,isppol,ispinor,ispinor)
               MSG_BUG(message)
             end if
             if(nsppol==1.and.nspinor==1) then
               elevels(nspinor*(2*lpawu+1)+iacc)=level_diag(iatom)%mat(im1,im1,isppol,ispinor,ispinor)
             end if
!            ! change it by: real(level_diag) with warning
           end do
         end do
       end do

!      ===================================
!      Compute possible occupations
!      ===================================
!      Value for nelec=0:
       nconfig_nelec(0)=1
       occ_level(0)%ocp(1,:)=0
!      Loop on possible occupation of levels with nelec
       do nelec=1,nlevels ! number of electrons
!        write(message,'(2a,i3,a)') ch10," For number of electrons",  &
!        &       nelec," positions of electrons are:"
!        call wrtout(std_out,message,'COLL')
!        write(std_out,*) "nelec",nelec
!        write(std_out,*) "nlevels",nlevels
         call combin(1,nconfig,nconfig_nelec,nelec,nlevels,occ_level,occup)
         if(nconfig_nelec(nelec)/=nint(permutations(nlevels,nelec)/factorial(nelec))) then
           message = " BUG in hubbard_one/combin"
           MSG_BUG(message)
         end if
         occ_level(nelec)%ocp=zero
         do iconfig=1,nconfig_nelec(nelec)
           do ielec=1,nelec
!            occ_level%repart: gives the place of electron ielec for the configuration iconfig (among the config for the total number of electron nelec
             occupied_level=occ_level(nelec)%repart(iconfig,ielec)
!            occ_level%ocp: gives if level occupied_level is occupied or not
             occ_level(nelec)%ocp(iconfig,occupied_level)=1
           end do
         end do
       end do

!      ===================================
!      Print possible occupations
!      ===================================
       if(prtopt>3) then
         do nelec=0,nlevels ! number of electrons f
           write(message,'(2a,i3,2a,i5,a)') ch10," For",nelec," electrons, ", &
&           "there are ",nconfig_nelec(nelec)," repartitions which are:"
           call wrtout(std_out,message,'COLL')
           do iconfig=1,nconfig_nelec(nelec)
             write(message,'(40i4)') (occ_level(nelec)%ocp(iconfig,ilevel),ilevel=1,nlevels),&
&             (occ_level(nelec)%repart(iconfig,ielec),ielec=1,nelec)
             call wrtout(std_out,message,'COLL')
           end do
         end do
       end if

!      ============================================
!      Compute energy for each of the occupations
!      ============================================
       do nelec=0,nlevels !
         e_nelec(nelec)%config=zero
         do iconfig=1,nconfig_nelec(nelec)
!          First compute energy level contribution
           do ielec=1,nelec
             e_nelec(nelec)%config(iconfig)= e_nelec(nelec)%config(iconfig) &
&             + elevels(occ_level(nelec)%repart(iconfig,ielec))
           end do
!          write(std_out,*) "Nelec",nelec,"iconfig",iconfig,"eleve",e_nelec(nelec)%config(iconfig)

!          Second: Compute interaction part
!          do ielec=1,nelec-1 ! compute interaction among the nelec electrons in the configuration iconfig
!          e_nelec(nelec)%config(iconfig)= e_nelec(nelec)%config(iconfig)   &
!          &             + hu(cryst_struc%typat(iatom))%udens(occ_level(nelec)%repart(iconfig,ielec), &
!          &               occ_level(nelec)%repart(iconfig,ielec+1))
!          enddo
           do ielec=1,nelec ! compute interaction among the nelec electrons in the configuration iconfig
             do jelec=1,nelec
               e_nelec(nelec)%config(iconfig)= e_nelec(nelec)%config(iconfig)   &
!              &               + hu(cryst_struc%typat(iatom))%udens(occ_level(nelec)%repart(iconfig,ielec), &
&               + udens_atoms(iatom)%value(occ_level(nelec)%repart(iconfig,ielec), &
&               occ_level(nelec)%repart(iconfig,jelec))/2.d0 ! udens(i,i)=0
!              write(std_out,*) ielec,occ_level(nelec)%repart(iconfig,ielec)
!              write(std_out,*) jelec,occ_level(nelec)%repart(iconfig,jelec)
!              write(std_out,*)hu(cryst_struc%typat(iatom))%udens(occ_level(nelec)%repart(iconfig,ielec), &
!              &                occ_level(nelec)%repart(iconfig,jelec))/2.d0
             end do ! jelec
           end do ! ielec
!          write(std_out,*) "Nelec",nelec,"iconfig",iconfig,"ecorr",e_nelec(nelec)%config(iconfig)

         end do ! iconfig
         maxener(nelec)=maxval(-e_nelec(nelec)%config(:))
         minener(nelec)=minval(-e_nelec(nelec)%config(:))
       end do
!      write(std_out,*) "maxener", maxener(:)
       emax=maxval(maxener(:))
       emin=minval(minener(:))
       eshift=zero
       eshift=emax/two
       eshift=emax-maxarg_exp/beta
!      eshift=emax
!      write(std_out,*)"emax",emax
!      write(std_out,*)"emin",emin
!      write(std_out,*)"eshift",eshift
       write(message,'(a,3x,3a,3x,a)') ch10," Hubbard I: Energies as a", &
&       " function of number of electrons",ch10,&
&       "     Nelec     Min. Ene.       Max. Ener."
       call wrtout(std_out,message,'COLL')
       do nelec=0,nlevels
         write(message,'(3x,a,i4,2f17.7)') "HI", nelec,&
&         minval(e_nelec(nelec)%config(:)),maxval(e_nelec(nelec)%config(:))
         call wrtout(std_out,message,'COLL')
       end do

!      ===================================
!      Print possibles occupations
!      ===================================
       if(prtopt>3) then
         do nelec=0,nlevels ! number of electrons
           write(message,'(2a,i3,2a,i5,3a)') ch10," For",nelec," electrons, ", &
&           "there are ",nconfig_nelec(nelec)," repartitions which are :", &
&           ch10,"Energy and Occupations"
           call wrtout(std_out,message,'COLL')
           do iconfig=1,nconfig_nelec(nelec)
             write(message,'(f12.6,20i4)') e_nelec(nelec)%config(iconfig),&
&             (occ_level(nelec)%repart(iconfig,ielec),ielec=1,nelec)
             call wrtout(std_out,message,'COLL')
           end do
         end do
       end if

!      sum_test=zero
!      do ielec=1,nelec+1
!      sum_test = sum_test + (occ_level(nelec)%repart(iconfig,ielec)  &
!      &              -occ_level(nelec)%repart(iconfig,ielec))
!      enddo
!      ===================================
!      Built transitions between configurations
!      ===================================
       do nelec=0,nlevels-1
         do iconfig=1,nconfig_nelec(nelec)
           itrans=0 ! transition from iconfig
           do jconfig=1, nconfig_nelec(nelec+1)
             sum_test=0
             do ilevel=1,nlevels
!              test if their is one electron added to the starting configuration
               sum_test=sum_test + &
&               (occ_level(nelec+1)%ocp(jconfig,ilevel)- occ_level(nelec)%ocp(iconfig,ilevel))**2
!              save the level for the electron added
               if(occ_level(nelec+1)%ocp(jconfig,ilevel)==1.and.occ_level(nelec)%ocp(iconfig,ilevel)==0) then
                 m_temp=ilevel
               end if
             end do ! ilevel
             if(sum_test==1) then
               itrans=itrans+1
               if(itrans>nlevels-nelec) then
                 write(message,'(a,4i4)') "BUG: itrans is to big in hubbard_one",itrans,iconfig,jconfig,ilevel
                 call wrtout(std_out,message,'COLL')
               end if
               occ_level(nelec)%transition(iconfig,itrans)=jconfig  ! jconfig=config(n+1) obtained after transition
               occ_level(nelec)%transition_m(iconfig,itrans)=m_temp  !  level to fill to do the transition
             end if
           end do ! jconfig
           if(prtopt>3) then
             write(std_out,'(a,2i5,a,18i5)') "occ_level", nelec,&
&             iconfig,"  :",(occ_level(nelec)%transition(iconfig,itrans),itrans=1,nlevels-nelec)
             write(std_out,'(a,2i5,a,18i5)') "electron added", nelec,iconfig,&
&             "  :",(occ_level(nelec)%transition_m(iconfig,itrans),itrans=1,nlevels-nelec)
           end if
         end do ! iconfig
       end do ! nelec

!      ===================================
!      Built Partition Function
!      ===================================
       Z_part=czero
!      do nelec=1,nlevels-1
       do nelec=0,nlevels
         do iconfig=1,nconfig_nelec(nelec)
           Ei_n    = e_nelec  (nelec  )%config(iconfig) + eshift
           Z_part=Z_part+dexp(-Ei_n*beta)
!          write(std_out,*) "fonction de partition",nelec,iconfig, Z_part,Ei_n*beta,Ei_n,eshift
         end do
       end do
!      write(std_out,*) "Z_part",Z_part

!      ===================================
!      Built Green Function
!      ===================================
       ABI_ALLOCATE(green_temp,(green_hubbard%nw,nlevels))
       ABI_ALLOCATE(green_temp_realw,(green_hubbard%nw,nlevels))
!      For each freq.

       green_temp=czero
       green_temp_realw=czero
       tmp=zero
       do nelec=0,nlevels-1
!        write(std_out,*) "For nelec    =",nelec
         do iconfig=1,nconfig_nelec(nelec)
!          write(std_out,*) "The config nb:",iconfig
           do itrans=1,nlevels-nelec
             jconfig = occ_level(nelec  )%transition(iconfig,itrans)
             m_temp  = occ_level(nelec  )%transition_m(iconfig,itrans)
             Ej_np1  = e_nelec  (nelec+1)%config(jconfig) + eshift
             Ei_n    = e_nelec  (nelec  )%config(iconfig) + eshift
!            write(std_out,'(a,i4,a)') "Transition nb:",itrans,"involve"
!            write(std_out,'(a,i4,a)') "                        jconfig=",jconfig
!            write(std_out,'(a,i4,a)') "                        m_temp=",m_temp
             do ifreq=1,green_hubbard%nw
               if(green_hubbard%w_type=="imag") then
                 omega_current=cmplx(zero,green_hubbard%omega(ifreq),kind=dp)
               else if(green_hubbard%w_type=="real") then
                 omega_current=cmplx(green_hubbard%omega(ifreq),zero,kind=dp)
               end if
               green_temp(ifreq,m_temp)=green_temp(ifreq,m_temp)+  &
&               (dexp(-Ej_np1*beta)+ dexp(-Ei_n*beta))/ &
&               ( omega_current +Ei_n-Ej_np1)
               if(ifreq==1.and.m_temp==1) tmp=tmp+(dexp(-Ej_np1*beta)+ dexp(-Ei_n*beta))

               green_temp_realw(ifreq,m_temp)=green_temp_realw(ifreq,m_temp)+  &
&               (dexp(-Ej_np1*beta)+ dexp(-Ei_n*beta))/ &
&               ( omega_current +Ei_n-Ej_np1)
             end do
!            green_temp_realw(m_temp)=green_temp_realw(m_temp)+  &
!            &           (dexp(-Ej_np1*beta)+ dexp(-Ei_n*beta)) -> will give one at the end
!            write(std_out,*) "green",-Ej_np1*beta,-Ei_n*beta,dexp(-Ej_np1*beta),dexp(-Ei_n*beta)
           end do
         end do
       end do
!      write(std_out,*) "tmp",tmp
       ilevel=0
       do ispinor=1,nspinor
         do isppol=1,nsppol
           do im=1,(2*lpawu+1)
             ilevel=ilevel+1
!            write(std_out,'(16e15.6)') paw_dmft%omega_lo(ifreq),(real(green_temp_realw(ilevel)/Z_part),ilevel=1,nlevels)
             do ifreq=1,green_hubbard%nw
               green_hubbard%oper(ifreq)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)=green_temp(ifreq,ilevel)/Z_part
               green_hubbard_realw%oper(ifreq)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)=green_temp_realw(ifreq,ilevel)/Z_part
             end do
           end do
         end do
       end do

!      End calculation for this frequency
       ABI_DEALLOCATE(green_temp)
       ABI_DEALLOCATE(green_temp_realw)

!      ===================================
!      Deallocations
!      ===================================
       do nelec=0,nlevels
         ABI_DEALLOCATE(occ_level(nelec)%repart)
         ABI_DEALLOCATE(occ_level(nelec)%ocp)
         ABI_DEALLOCATE(occ_level(nelec)%transition)
         ABI_DEALLOCATE(occ_level(nelec)%transition_m)
         ABI_DEALLOCATE(e_nelec(nelec)%config)
       end do
       ABI_DATATYPE_DEALLOCATE(occ_level)
       ABI_DEALLOCATE(occup)
       ABI_DEALLOCATE(nconfig_nelec)
       ABI_DATATYPE_DEALLOCATE(e_nelec)
       ABI_DEALLOCATE(elevels)
       ABI_DEALLOCATE(maxener)
       ABI_DEALLOCATE(minener)
     end if
   end do
   call destroy_green(green_hubbard_realw)

 end subroutine green_atomic_hubbard
!!***

!!****f* m_hubbard_one/combin
!! NAME
!! combin
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!! Copyright (C) 1999-2020 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 recursive subroutine combin(ielec,nconfig,nconfig_nelec,nelec,nlevels,occ_level,occup)

 use defs_basis
 implicit none

!Arguments ------------------------------------
!scalars
! type(pawang_type), intent(in) :: pawang
 integer, intent(in) :: ielec,nelec,nlevels
 integer, intent(inout) :: nconfig,nconfig_nelec(0:nlevels)
 integer, intent(inout) :: occup(0:nlevels,nlevels)
! type  :: level2_type
!  integer, pointer :: repart(:,:)
! end type
 type(level2_type), intent(inout) :: occ_level(0:nlevels)
! integer, intent(in) :: prtopt

!Local variables ------------------------------
! scalars
 integer :: max_ielec,pos,min_ielec,jelec,prtopt
 character(len=500) :: message
! arrays
!************************************************************************
   prtopt=1
   max_ielec=nlevels-nelec+ielec
!  write(std_out,*) "call to combin ielec,nelec,nlevels",ielec,nelec,nlevels
   select case (ielec)
   case (1)
     min_ielec=1
   case default
     min_ielec=occup(nelec,ielec-1)+1
   end select
!  write(std_out,*) "For ielec", ielec, "min_ielec,max_ielec",min_ielec,max_ielec
   do pos = min_ielec, max_ielec
     if(ielec==nelec) then
       occup(nelec,ielec)=pos
       nconfig=nconfig+1
       nconfig_nelec(nelec)=nconfig_nelec(nelec)+1
!      write(std_out,*) "size",size(occ_level),size(occ_level(nelec)%repart,1)
!      write(std_out,*) "size",size(occ_level),size(occ_level(nelec)%repart,2)
       do jelec=1,nelec
!        write(std_out,*) "nconfig",nconfig_nelec(nelec),nelec
!        write(std_out,*) "occup",occup(nelec,jelec)
         occ_level(nelec)%repart(nconfig_nelec(nelec),jelec)=occup(nelec,jelec)
       end do
!      write(std_out,*) "For ielec", ielec, "case nelec"
       if(prtopt>=3) then
         write(message,'(a,i3,a,30i5)') "For ielec",ielec," Occupf are", (occup(nelec,jelec),jelec=1,nelec)
         call wrtout(std_out,message,'COLL')
       end if
     else
       occup(nelec,ielec)=pos
!      write(std_out,*) "For ielec", ielec, "case 1 and default"
       call combin(ielec+1,nconfig,nconfig_nelec,nelec,nlevels,occ_level,occup)
       if(prtopt>=3) then
         write(message,'(a,i3,a,30i5)') "For ielec",ielec," Occup are", (occup(nelec,jelec),jelec=1,nelec)
         call wrtout(std_out,message,'COLL')
       end if
     end if
   end do

 end subroutine combin
!!***

end subroutine hubbard_one
!!***

END MODULE m_hubbard_one
!!***
