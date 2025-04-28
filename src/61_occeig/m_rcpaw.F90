!!****m* ABINIT/m_rcpaw
!! NAME
!!  m_rcpaw
!!
!! FUNCTION
!!  
!!
!! COPYRIGHT
!!  Copyright (C) 2019-2019 ABINIT group (NBrouwer,MT, JBoust)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_rcpaw
 use defs_basis
 use defs_abitypes
 use m_dtset
 use m_pawtab
 use m_pawrad
 use m_xmpi
 use m_errors
 use m_paw_atomorb
 use m_paw_atom
 use m_paral_atom
 use m_paw_relax
 use m_pawpsp,           only : pawpsp_init_core
 use m_extfpmd,          only : extfpmd_type
 use defs_datatypes,     only : pseudopotential_type
 use m_pawang,           only : pawang_type
 use m_pawrhoij,         only : pawrhoij_type
 use m_paw_an,           only : paw_an_type 
 use m_pawfgrtab,        only : pawfgrtab_type
 use m_paw_finegrid,     only : pawrfgd_fft

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_rcpaw/valdens_type
!! NAME
!! valdens_type
!!
!! FUNCTION
!!
!! SOURCE
 type,public :: valdens_type
   logical :: has_dens
   real(dp) :: compch_sph
   real(dp), allocatable :: rho1(:,:,:)
   real(dp), allocatable :: trho1(:,:,:)
   real(dp), allocatable :: nhat1(:,:,:)
 end type valdens_type

!----------------------------------------------------------------------

!!****t* m_rcpaw/rcpaw_type
!! NAME
!! rcpaw_type
!!
!! FUNCTION
!!
!! SOURCE
 type,public :: rcpaw_type
   real(dp) :: nelect_core
   real(dp) :: nelect_core_orig
   real(dp) :: ehnzc
   real(dp) :: ekinc
   real(dp) :: edcc
   real(dp) :: eeigc
   real(dp) :: entropy
   real(dp) :: tolnc
   integer :: ntypat
   integer :: istep
   integer :: nfrpaw
   integer :: nfrocc
   logical :: all_atoms_relaxed
   type(atomorb_type),allocatable :: atm(:)
   type(atompaw_type),allocatable :: atp(:)
   type(valdens_type),allocatable :: val(:)
   logical, allocatable :: eijkl_is_sym(:)
   integer :: nfrtnc
   logical :: frocc
 end type rcpaw_type
!!***

!----------------------------------------------------------------------

 public :: rcpaw_destroy       ! Destroy RCPAW
 public :: rcpaw_init          ! Initialize RCPAW
 public :: rcpaw_core_eig      ! Compute core eigenenergies
 public :: rcpaw_core_energies ! Compute total energy contributions from the core
!!***


CONTAINS !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_rcpaw/rcpaw_destroy
!! NAME
!! rcpaw_destroy
!!
!! FUNCTION
!!  Destroy RCPAW object
!!
!! INPUTS
!!
!!
!! OUTPUT
!!
!!
!! SOURCE

subroutine rcpaw_destroy(rcpaw)
!Arguments ------------------------------------
!scalars
 integer :: ii
 type(rcpaw_type), pointer, intent(inout) :: rcpaw

!******************************************************************************************

 if(allocated(rcpaw%atm)) then
   do ii=1,size(rcpaw%atm)
     call destroy_atomorb(rcpaw%atm(ii))
   enddo
   ABI_FREE(rcpaw%atm)
 endif
 if(allocated(rcpaw%atp)) then
   do ii=1,size(rcpaw%atp)
     call atompaw_destroy(rcpaw%atp(ii))
   enddo
   ABI_FREE(rcpaw%atp)
 endif
 if(allocated(rcpaw%val)) then
   do ii=1,size(rcpaw%val)
     call destroy_valdens(rcpaw%val(ii))
   enddo
   ABI_FREE(rcpaw%val)
 endif
 if(allocated(rcpaw%eijkl_is_sym) ABI_FREE(rcpaw%eijkl_is_sym)
 rcpaw=>null()

end subroutine rcpaw_destroy
!!***


!----------------------------------------------------------------------

!!****f* m_rcpaw/destroy_valdens
!! NAME
!! destroy_valdens
!!
!! FUNCTION
!!  Destroy valdens object
!!
!! INPUTS
!!
!!
!! OUTPUT
!!
!!
!! SOURCE

subroutine destroy_valdens(val)
!Arguments ------------------------------------
!scalars
 type(valdens_type), intent(inout) :: val

!******************************************************************************************

 val%compch_sph=zero
 val%has_dens=.false.
 if (allocated(val%rho1)) ABI_FREE(val%rho1)
 if (allocated(val%trho1)) ABI_FREE(val%trho1)
 if (allocated(val%nhat1)) ABI_FREE(val%nhat1)

end subroutine destroy_valdens
!!***



!----------------------------------------------------------------------

!!****f* m_rcpaw/rcpaw_init
!! NAME
!! rcpaw_init
!!
!! FUNCTION
!!  Initialize the RCPAW functionality
!!
!! INPUTS
!!
!!
!! OUTPUT
!!
!!
!! SOURCE

subroutine rcpaw_init(rcpaw,dtset,filpsp,pawrad,pawtab,psps,ntypat,paw_an,my_natom,comm_atom,mpi_atmtab)
!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ntypat,my_natom
 integer,optional,intent(in) :: comm_atom
 type(pseudopotential_type),intent(inout) :: psps
 type(rcpaw_type), pointer, intent(inout) :: rcpaw
 type(dataset_type), intent(in) :: dtset
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 character(len=fnlen), intent(in) :: filpsp(ntypat)
 type(pawrad_type), intent(in) :: pawrad(ntypat)
 type(pawtab_type), intent(inout) :: pawtab(ntypat)
 type(paw_an_type), intent(in) :: paw_an(my_natom)

!Local variables-------------------------------
!scalars
 integer :: itypat,iatom,lm_size,mesh_size,cplex,my_comm_atom,iat,iln
 logical :: my_atmtab_allocated,paral_atom
!arrays
 integer :: mult(ntypat)
 integer,pointer :: my_atmtab(:)

!******************************************************************************************

 !Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=dtset%natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,dtset%natom,my_natom_ref=my_natom)

 ABI_MALLOC(rcpaw%eijkl_is_sym,(dtset%ntypat))
 rcpaw%eijkl_is_sym(:)=.true.

 ! Initilizations
 rcpaw%istep=0
 mult=0
 do iatom=1,dtset%natom
   itypat=dtset%typat(iatom)
   mult(itypat)=mult(itypat)+1
 enddo
 rcpaw%ntypat=ntypat
 rcpaw%edcc=zero
 rcpaw%eeigc=zero
 rcpaw%ehnzc=zero
 rcpaw%ekinc=zero
 rcpaw%tolnc=dtset%rcpaw_tolnc

 ! Init arrays 
 if(.not.allocated(rcpaw%val)) ABI_MALLOC(rcpaw%val,(my_natom))
 if(.not.allocated(rcpaw%atm)) ABI_MALLOC(rcpaw%atm,(ntypat))
 if(.not.allocated(rcpaw%atp)) ABI_MALLOC(rcpaw%atp,(ntypat))
 rcpaw%all_atoms_relaxed=.true.
 do itypat=1,ntypat
   call pawpsp_init_core(rcpaw%atm(itypat),psp_filename=filpsp(itypat))
   ABI_MALLOC(rcpaw%atm(itypat)%vhtnzc_orig,(size(pawtab(itypat)%vhtnzc)))
   rcpaw%atm(itypat)%mode=dtset%rcpaw_frtypat(itypat)
   rcpaw%atm(itypat)%vhtnzc_orig=pawtab(itypat)%vhtnzc
   rcpaw%atm(itypat)%mult=mult(itypat)
   rcpaw%atm(itypat)%nspden=dtset%nspden
   if(rcpaw%atm(itypat)%mode(1,1)==orb_relaxed_core) rcpaw%all_atoms_relaxed=.false.
 enddo

 ! Init atp
 do itypat=1,ntypat
   if(rcpaw%atm(itypat)%zcore_orig>zero) then
     rcpaw%atp(itypat)%ixc=dtset%ixc
     rcpaw%atp(itypat)%xclevel=dtset%xclevel
     rcpaw%atp(itypat)%electrons=dtset%nelect
     call atompaw_init(itypat,pawtab(itypat),pawrad(itypat),rcpaw%atp(itypat),&
&    int(rcpaw%atm(itypat)%znucl),rcpaw%atm(itypat),dtset%rcpaw_verbosity,dtset%rcpaw_scenergy(itypat))
   endif
 enddo
 
 ! Init val
 do iat=1,my_natom
  iatom=iat;if (paral_atom) iatom=my_atmtab(iat)
  cplex=paw_an(iat)%cplex
  lm_size=paw_an(iat)%lm_size
  itypat=dtset%typat(iatom)
  mesh_size=pawtab(itypat)%mesh_size
  ABI_MALLOC(rcpaw%val(iat)%nhat1,(mesh_size*cplex,lm_size,dtset%nspden))
  ABI_MALLOC(rcpaw%val(iat)%rho1,(mesh_size*cplex,lm_size,dtset%nspden))
  ABI_MALLOC(rcpaw%val(iat)%trho1,(mesh_size*cplex,lm_size,dtset%nspden))
  rcpaw%val(iat)%nhat1=zero
  rcpaw%val(iat)%rho1=zero
  rcpaw%val(iat)%trho1=zero
  rcpaw%val(iat)%compch_sph=zero
  rcpaw%val(iat)%has_dens=.false.
 enddo

 ! Init non arrays
 rcpaw%nelect_core=zero
 do itypat=1,ntypat
   rcpaw%nelect_core=rcpaw%nelect_core+rcpaw%atm(itypat)%zcore*rcpaw%atm(itypat)%mult
 enddo
 rcpaw%nelect_core_orig=rcpaw%nelect_core
 if(dtset%rcpaw_frocc==1) then 
   rcpaw%frocc=.true.
 else
   rcpaw%frocc=.false.
 endif
 rcpaw%nfrpaw=dtset%rcpaw_nfrpaw
 if(rcpaw%frocc) then
   rcpaw%nfrocc=rcpaw%nfrpaw
 else
   rcpaw%nfrocc=dtset%nstep
 endif
 rcpaw%nfrtnc=dtset%rcpaw_nfrtnc

 call rcpaw_core_energies(rcpaw,pawtab,ntypat)

end subroutine rcpaw_init
!!***


!----------------------------------------------------------------------

!!****f* m_rcpaw/rcpaw_core_eig
!! NAME
!! rcpaw_core_eig
!!
!! FUNCTION
!! Computes the core eigenenergies
!!
!! INPUTS
!!
!!
!! OUTPUT
!!
!!
!! SOURCE

subroutine rcpaw_core_eig(pawtab,pawrad,ntypat,rcpaw,dtset,&
& nfft,pawfgrtab,vtrial,cplex,ucvol,paw_an,&
&                      gmet,rprimd,xred,ngfft,my_natom,&
&                      distribfft,comm_fft,mpi_atmtab,comm_atom)
!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ntypat,cplex
 integer,intent(in) :: nfft,my_natom
 integer,optional,intent(in) :: comm_atom
 integer,optional,intent(in) :: comm_fft
 real(dp), intent(in) :: ucvol
 type(distribfft_type),optional,target,intent(in)  :: distribfft
 type(rcpaw_type), intent(inout) :: rcpaw
 type(dataset_type), intent(in) :: dtset
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 integer,intent(in) :: ngfft(18)
 type(pawtab_type), target,intent(inout) :: pawtab(ntypat)
 type(pawrad_type), intent(in) :: pawrad(ntypat)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
 real(dp),intent(in),target :: vtrial(cplex*nfft,dtset%nspden)
 type(paw_an_type),intent(inout) :: paw_an(my_natom)
 real(dp),intent(in) :: gmet(3,3)
 real(dp), intent(in) :: rprimd(3,3)
 real(dp), intent(in) :: xred(3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer :: me_fft,iatom
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 integer :: ii,itypat
 integer :: mesh_size
 integer :: il,nfgd,ifft
 integer :: n1,n2,n3,i3,ispden
 integer :: my_comm_atom,iat,ierr
 logical :: my_atmtab_allocated,paral_atom,grid_found
 real(dp) :: eigshift,r2_tmp
!arrays
 integer,pointer :: my_atmtab(:) 
 real(dp), allocatable :: nt1hat0(:)
 integer,allocatable :: ifftsph(:)
 real(dp),allocatable :: rfgd(:,:)
 real(dp),allocatable :: vh_sph(:)

!******************************************************************************************
 
 if(.not.rcpaw%all_atoms_relaxed) then
   ! Initializations
   if(cplex.ne.1) ABI_ERROR('cplex not 1')
   paral_atom=(present(comm_atom).and.(my_natom/=dtset%natom))
   nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
   my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
   call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,dtset%natom,my_natom_ref=my_natom)
   me_fft=0
   if (present(comm_fft)) then
     me_fft=xmpi_comm_rank(comm_fft)
   end if
   n1=ngfft(1);n2=ngfft(2);n3=ngfft(3)
   if (present(distribfft)) then
     grid_found=.false.
     if (n2 == distribfft%n2_coarse) then
       if (n3== size(distribfft%tab_fftdp3_distrib)) then
         fftn3_distrib => distribfft%tab_fftdp3_distrib
         ffti3_local => distribfft%tab_fftdp3_local
         grid_found=.true.
       end if
     end if
     if (n2 == distribfft%n2_fine) then
       if (n3 == size(distribfft%tab_fftdp3dg_distrib)) then
         fftn3_distrib => distribfft%tab_fftdp3dg_distrib
         ffti3_local => distribfft%tab_fftdp3dg_local
         grid_found = .true.
       end if
     end if
     if (.not.(grid_found)) then
       ABI_BUG('Unable to find an allocated distrib for this fft grid!')
     end if
   else
     ABI_MALLOC(fftn3_distrib,(n3))
     ABI_MALLOC(ffti3_local,(n3))
     fftn3_distrib=0;ffti3_local=(/(i3,i3=1,n3)/)
   end if
 endif

   ! Loop on typat 
 do itypat=1,dtset%ntypat
   if(.not.rcpaw%atm(itypat)%nc_conv) then
     eigshift=zero
     do iat=1,my_natom
       iatom=iat;if (paral_atom) iatom=my_atmtab(iat)
       if(dtset%typat(iatom)==itypat) then
         mesh_size=pawtab(itypat)%mesh_size
         ABI_MALLOC(nt1hat0,(pawtab(itypat)%mesh_size))
         nt1hat0=zero
         do ispden=1,dtset%nspden
           nt1hat0(1:pawtab(itypat)%mesh_size)=nt1hat0(1:pawtab(itypat)%mesh_size)+&
&           rcpaw%val(iat)%trho1(1:pawtab(itypat)%mesh_size,1,ispden)*sqrt(four*pi)*&
&           pawrad(itypat)%rad(1:pawtab(itypat)%mesh_size)**2+&
&           rcpaw%val(iat)%nhat1(1:pawtab(itypat)%mesh_size,1,ispden)*sqrt(four*pi)*&
&           pawrad(itypat)%rad(1:pawtab(itypat)%mesh_size)**2
         end do
         ABI_MALLOC(vh_sph,(mesh_size))
         call poisson(nt1hat0,0,pawrad(itypat),vh_sph)
         do il=2,mesh_size
           vh_sph(il)=vh_sph(il)/pawrad(itypat)%rad(il)
         enddo
         call pawrad_deducer0(vh_sph,mesh_size,pawrad(itypat))
         ABI_FREE(nt1hat0)
         call pawrfgd_fft(ifftsph,gmet,n1,n2,n3,nfgd,0.5_dp,rfgd,rprimd,ucvol,xred(:,iatom),&
&                        fftn3_distrib,ffti3_local,me_fft) 
         ifft=ifftsph(1)
         r2_tmp=norm2(rfgd(:,1))
         do ii=2,nfgd
           if(norm2(rfgd(:,ii))<r2_tmp) then
             ifft=ifftsph(ii)
             r2_tmp=norm2(rfgd(:,ii))
           endif
         enddo 
         ABI_FREE(ifftsph)
         ABI_FREE(rfgd)
         eigshift=eigshift+vtrial(ifft,1)-vh_sph(1)-pawtab(itypat)%vhtnzc(1)-paw_an(iat)%vxct1(1,1,1)/sqrt(four_pi)
         ABI_FREE(vh_sph)
       endif
     enddo
     !!! mpi reduction
     if(paral_atom) then
       call xmpi_sum(eigshift,my_comm_atom,ierr)
       call xmpi_bcast(eigshift,0,my_comm_atom,ierr)
     endif
     rcpaw%atm(itypat)%eig=rcpaw%atm(itypat)%eig+eigshift/rcpaw%atm(itypat)%mult
     if(rcpaw%atm(itypat)%nresid_c<rcpaw%tolnc)rcpaw%atm(itypat)%nc_conv=.true.
   endif
 enddo

 if(.not.rcpaw%all_atoms_relaxed) then
   if (.not.present(distribfft)) then
     ABI_FREE(fftn3_distrib)
     ABI_FREE(ffti3_local)
   end if
 endif

 rcpaw%all_atoms_relaxed=.true.
 do itypat=1,dtset%ntypat 
   if(rcpaw%atm(itypat)%zcore_conv.and.rcpaw%atm(itypat)%nc_conv) then
     rcpaw%atm(itypat)%mode(:,:)=ORB_FROZEN
   else
     rcpaw%all_atoms_relaxed=.false.
   endif
 enddo


end subroutine rcpaw_core_eig
!!***



!----------------------------------------------------------------------

!!****f* m_rcpaw/rcpaw_core_energies
!! NAME
!! rcpaw_core_energies
!!
!! FUNCTION
!!
!! INPUTS
!!
!!
!! OUTPUT
!!
!!
!! SOURCE

subroutine rcpaw_core_energies(rcpaw,pawtab,ntypat)
!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ntypat
 type(rcpaw_type), pointer, intent(inout) :: rcpaw
!arrays
 type(pawtab_type), intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: itypat,iln

!******************************************************************************************

 rcpaw%ehnzc=zero
 rcpaw%edcc=zero
 rcpaw%eeigc=zero
 rcpaw%ekinc=zero
 do itypat=1,ntypat
   if(rcpaw%atm(itypat)%zcore_orig>zero) then
     rcpaw%edcc=rcpaw%edcc+rcpaw%atm(itypat)%edcc*rcpaw%atm(itypat)%mult
     rcpaw%ekinc=rcpaw%ekinc+rcpaw%atm(itypat)%ekinc*rcpaw%atm(itypat)%mult
     rcpaw%eeigc=rcpaw%eeigc+rcpaw%atm(itypat)%eeigc*rcpaw%atm(itypat)%mult
     rcpaw%ehnzc=rcpaw%ehnzc+rcpaw%atm(itypat)%ehnzc*rcpaw%atm(itypat)%mult
   endif
 enddo
 write(std_out,*) rcpaw%istep,'ekinc=',rcpaw%ekinc, 'eeig=',rcpaw%eeigc,'edcc=',rcpaw%edcc,'ehnzc',rcpaw%ehnzc
 write(std_out,*) rcpaw%atm(1)%mult,rcpaw%atm(1)%edcc,ntypat

end subroutine rcpaw_core_energies
!!***


end module m_rcpaw
!!***

