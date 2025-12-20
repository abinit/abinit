!!****m* ABINIT/m_rcpaw
!! NAME
!!  m_rcpaw
!!
!! FUNCTION
!! This module contains types and subroutines linked to the PAW core relaxation
!!  approach  
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
 use m_paw_atom_solve
 use m_pawpsp,           only : pawpsp_init_core
 use m_extfpmd,          only : extfpmd_type
 use defs_datatypes,     only : pseudopotential_type
 use m_pawang,           only : pawang_type
 use m_pawrhoij,         only : pawrhoij_type
 use m_paw_an,           only : paw_an_type 
 use m_pawfgrtab,        only : pawfgrtab_type
 use m_paw_finegrid,     only : pawrfgd_fft

#ifdef HAVE_MPI2
 use mpi
#endif


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
!!***

!----------------------------------------------------------------------

!!****t* m_rcpaw/rcpaw_type
!! NAME
!! rcpaw_type
!!
!! FUNCTION
!!
!! SOURCE
 type,public :: rcpaw_type
   integer :: ntypat
   integer :: istep
   integer :: updatepaw(2)
   integer :: updateocc
   integer :: updatetnc
   logical :: frocc
   logical :: all_atoms_relaxed
   real(dp) :: nelect_core
   real(dp) :: nelect_core_orig
   real(dp) :: ehnzc
   real(dp) :: ekinc
   real(dp) :: edcc
   real(dp) :: eeigc
   real(dp) :: entropy
   real(dp) :: tolnc
   logical, allocatable :: eijkl_is_sym(:)
   type(atomorb_type),allocatable :: atm(:)
   type(atompaw_type),allocatable :: atp(:)
   type(valdens_type),allocatable :: val(:)
 end type rcpaw_type
!!***

!----------------------------------------------------------------------

 public :: rcpaw_destroy       ! Destroy RCPAW
 public :: rcpaw_init          ! Initialize RCPAW
 public :: rcpaw_reinit        ! Re-Initialize RCPAW
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
 type(rcpaw_type), pointer,intent(inout) :: rcpaw

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
 ABI_SFREE(rcpaw%eijkl_is_sym)

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
 ABI_SFREE(val%rho1)
 ABI_SFREE(val%trho1)
 ABI_SFREE(val%nhat1)

end subroutine destroy_valdens
!!***


!----------------------------------------------------------------------

!!****f* m_rcpaw/rcpaw_reinit
!! NAME
!! rcpaw_reinit
!!
!! FUNCTION
!! Reinitialize rcpaw object
!!
!! INPUTS
!!
!!
!! OUTPUT
!!
!!
!! SOURCE

subroutine rcpaw_reinit(rcpaw)
!Arguments ------------------------------------
!scalars
 integer :: itypat
 type(rcpaw_type), pointer,intent(inout) :: rcpaw

!******************************************************************************************

 rcpaw%all_atoms_relaxed=.true.
 do itypat=1,size(rcpaw%atm)
   rcpaw%atm(itypat)%nresid_c=one
   rcpaw%atm(itypat)%nc_conv=.false.
   if(rcpaw%atm(itypat)%mode(1,1)==orb_relaxed_core) then
     rcpaw%all_atoms_relaxed=.false.
   else
     rcpaw%atm(itypat)%nc_conv=.true.
   endif
 enddo 

end subroutine rcpaw_reinit
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

subroutine rcpaw_init(rcpaw,dtset,filpsp,pawrad,pawtab,ntypat,cplex,my_natom,comm_atom,mpi_atmtab)
!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ntypat,my_natom,cplex
 integer,optional,intent(in) :: comm_atom
 type(rcpaw_type), pointer, intent(inout) :: rcpaw
 type(dataset_type), intent(in) :: dtset
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 character(len=fnlen), intent(in) :: filpsp(ntypat)
 type(pawrad_type), intent(in) :: pawrad(ntypat)
 type(pawtab_type), intent(inout) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: itypat,iatom,lm_size,mesh_size,my_comm_atom,iat
 logical :: my_atmtab_allocated,paral_atom
!arrays
 integer :: mult(ntypat)
 integer,pointer :: my_atmtab(:)

!******************************************************************************************

 write(std_out, * ) 'RCPAW initialization'

 !Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=dtset%natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,dtset%natom,my_natom_ref=my_natom)

 ! Set up eijkl_is_sym
 ABI_MALLOC(rcpaw%eijkl_is_sym,(dtset%ntypat))
 rcpaw%eijkl_is_sym(:)=.true.

 ! Set up multiplicity of atoms
 rcpaw%istep=0
 mult=0
 do iatom=1,dtset%natom
   itypat=dtset%typat(iatom)
   mult(itypat)=mult(itypat)+1
 enddo

 ! Allocate arrays 
 if(.not.allocated(rcpaw%val)) then
   ABI_MALLOC(rcpaw%val,(my_natom))
 endif
 if(.not.allocated(rcpaw%atm)) then
   ABI_MALLOC(rcpaw%atm,(ntypat))
 endif
 if(.not.allocated(rcpaw%atp)) then
   ABI_MALLOC(rcpaw%atp,(ntypat))
 endif

 ! Init atm
 rcpaw%all_atoms_relaxed=.true.
 do itypat=1,ntypat
   call pawpsp_init_core(rcpaw%atm(itypat),psp_filename=filpsp(itypat))
   ABI_MALLOC(rcpaw%atm(itypat)%vhtnzc_orig,(size(pawtab(itypat)%vhtnzc)))
   rcpaw%atm(itypat)%mode=dtset%rcpaw_rctypat(itypat)
   rcpaw%atm(itypat)%vhtnzc_orig=pawtab(itypat)%vhtnzc
   rcpaw%atm(itypat)%mult=mult(itypat)
   rcpaw%atm(itypat)%nspden=dtset%nspden
   rcpaw%atm(itypat)%eigshift=zero
   if(rcpaw%atm(itypat)%mode(1,1)==orb_relaxed_core) then
     rcpaw%all_atoms_relaxed=.false.
   else
     rcpaw%atm(itypat)%nc_conv=.true.
   endif
 enddo

 ! Init atp
 do itypat=1,ntypat
     rcpaw%atp(itypat)%ixc=dtset%ixc
     rcpaw%atp(itypat)%xclevel=dtset%xclevel
     rcpaw%atp(itypat)%electrons=dtset%nelect
     call atompaw_init(pawtab(itypat),pawrad(itypat),rcpaw%atp(itypat),&
&    int(rcpaw%atm(itypat)%znucl),rcpaw%atm(itypat),dtset%rcpaw_sc(itypat),&
&    dtset%rcpaw_elin,dtset%rcpaw_vhtnzc,dtset%rcpaw_tpaw)
 enddo
 
 ! Init val
 do iat=1,my_natom
  iatom=iat;if (paral_atom) iatom=my_atmtab(iat)
  itypat=dtset%typat(iatom)
  lm_size=Pawtab(itypat)%lcut_size**2
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
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 ! Init non arrays
 rcpaw%edcc=zero
 rcpaw%eeigc=zero
 rcpaw%ehnzc=zero
 rcpaw%ekinc=zero
 rcpaw%entropy=zero
 rcpaw%tolnc=dtset%rcpaw_tolnc
 rcpaw%ntypat=ntypat
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
 rcpaw%updatepaw=dtset%rcpaw_updatepaw
 if(rcpaw%frocc.and.rcpaw%updatepaw(2)>0) then
   rcpaw%updateocc=rcpaw%updatepaw(2)
 else
   rcpaw%updateocc=dtset%nstep
 endif
 rcpaw%updatetnc=dtset%rcpaw_updatetnc

 ! Init core energies
 call rcpaw_core_energies(rcpaw,ntypat)
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
& nfft,vtrial,cplex,ucvol,&
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
 real(dp),intent(in) :: gmet(3,3)
 real(dp), intent(in) :: rprimd(3,3)
 real(dp), intent(in) :: xred(3,dtset%natom)
 real(dp),intent(in),target :: vtrial(cplex*nfft)
 type(pawtab_type), target,intent(inout) :: pawtab(ntypat)
 type(pawrad_type), intent(in) :: pawrad(ntypat)

!Local variables-------------------------------
!scalars
 integer :: me_fft,iatom
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 integer :: ii,itypat,ifft_old
 integer :: mesh_size,ind1,ind2
 integer :: il,nfgd,ifft,iln
 integer :: n1,n2,n3,i3,ispden
 integer :: my_comm_atom,iat,ierr
 logical :: my_atmtab_allocated,paral_atom,grid_found
 real(dp) :: eigshift,r1,r2,est_err,vh1,vh2
!arrays
 integer,pointer :: my_atmtab(:) 
 integer,allocatable :: ifftsph(:)
 real(dp), allocatable :: nt1hat0(:)
 real(dp),allocatable :: rfgd(:,:)
 real(dp),allocatable :: vh_sph(:)
 integer :: local_arr(2),global_arr(2)

!******************************************************************************************

 ! FFT grid 
 if(.not.rcpaw%all_atoms_relaxed) then
   if(cplex.ne.1) then
     ABI_ERROR('cplex not 1')
   endif
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
     est_err=zero
     do iat=1,my_natom
       iatom=iat;if (paral_atom) iatom=my_atmtab(iat)
       if(dtset%typat(iatom)==itypat) then ! Average on atoms of same type
         mesh_size=pawtab(itypat)%mesh_size
         ABI_MALLOC(nt1hat0,(pawtab(itypat)%mesh_size)) ! spherical part of nt1+nhat
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
         r1=0.6_dp
         r2=r1
         ifft=1
         ifft_old=1
         do ii=1,nfgd
           if(norm2(rfgd(:,ii))<r1) then
             ifft_old=ifft
             ifft=ifftsph(ii)
             r2=r1
             r1=norm2(rfgd(:,ii))
           endif
         enddo 
         ABI_FREE(ifftsph)
         ABI_FREE(rfgd)
         global_arr(1)=r1
         global_arr(2)=0
         if(present(comm_fft)) then
           local_arr(1)=r1
           local_arr(2)=me_fft
           ! Check which processor is closest to minimum
#if defined HAVE_MPI
           call MPI_ALLREDUCE(local_arr,global_arr,1,MPI_2INT,MPI_MINLOC,comm_fft,ierr)
#endif
         endif
         if(me_fft==global_arr(2)) then
           ind1=pawrad_ifromr(pawrad(itypat),r1)
           ind2=pawrad_ifromr(pawrad(itypat),r2)
           vh1=vh_sph(ind1)+pawtab(itypat)%vhtnzc(ind1)+(r1-pawrad(itypat)%rad(ind1))*&
&              (pawtab(itypat)%vhtnzc(ind1+1)+vh_sph(ind1+1)-vh_sph(ind1)-pawtab(itypat)%vhtnzc(ind1))/&
&              (pawrad(itypat)%rad(ind1+1)-pawrad(itypat)%rad(ind1)) 
           vh2=vh_sph(ind2)+pawtab(itypat)%vhtnzc(ind2)+(r2-pawrad(itypat)%rad(ind2))*&
&              (pawtab(itypat)%vhtnzc(ind2+1)+vh_sph(ind2+1)-vh_sph(ind2)-pawtab(itypat)%vhtnzc(ind2))/&
&              (pawrad(itypat)%rad(ind2+1)-pawrad(itypat)%rad(ind2))
           eigshift=eigshift+vtrial(ifft)-vh1
           if(r1>zero) then
              est_err=est_err+abs((vtrial(ifft)-vh1)-&
&                                  (vtrial(ifft_old)-vh2))
            endif
         endif
         if(present(comm_fft)) then
           call xmpi_bcast(eigshift,global_arr(2),comm_fft,ierr)
           call xmpi_bcast(est_err,global_arr(2),comm_fft,ierr)
         endif
         ABI_FREE(vh_sph)
       endif
     enddo
     ! mpi reduction
     if(paral_atom) then
       call xmpi_sum(eigshift,my_comm_atom,ierr)
       call xmpi_bcast(eigshift,0,my_comm_atom,ierr)
       call xmpi_sum(est_err,my_comm_atom,ierr)
       call xmpi_bcast(est_err,0,my_comm_atom,ierr)
     endif
     rcpaw%atm(itypat)%eigshift=eigshift/rcpaw%atm(itypat)%mult
     write(std_out,*) 'ESTIMATED ERROR ON CORE EIGS OF TYPAT',itypat,' = ',est_err*Ha_eV, ' eV'
     if(allocated(rcpaw%atm(itypat)%eig)) rcpaw%atm(itypat)%eig=rcpaw%atm(itypat)%eig+eigshift/rcpaw%atm(itypat)%mult ! Average on atoms of same type
     if(rcpaw%atm(itypat)%nresid_c<rcpaw%tolnc.and.rcpaw%istep>rcpaw%updatepaw(2).and.rcpaw%updatepaw(2)/=0)then
        rcpaw%atm(itypat)%nc_conv=.true.
     endif
   endif
 enddo

 !Destroy atom table used for parallelism
 if(.not.rcpaw%all_atoms_relaxed) then
   call free_my_atmtab(my_atmtab,my_atmtab_allocated)
   if (.not.present(distribfft)) then
     ABI_FREE(fftn3_distrib)
     ABI_FREE(ffti3_local)
   end if
 endif

 ! Update convergence status of cores
 rcpaw%all_atoms_relaxed=.true.
 do itypat=1,dtset%ntypat 
   if(rcpaw%atm(itypat)%zcore_conv.and.rcpaw%atm(itypat)%nc_conv) then
     rcpaw%atm(itypat)%mode(:,:)=ORB_FROZEN
   else
     rcpaw%all_atoms_relaxed=.false.
   endif
 enddo

 ! Print core eigenenergies and occupations
 do itypat=1,dtset%ntypat
   write(std_out,*) 'RCPAW core eigenergies (Ha) and occupations for typat ',itypat
   do iln=1,rcpaw%atm(itypat)%ln_size
     write(std_out,*) rcpaw%atm(itypat)%eig(iln,1),rcpaw%atm(itypat)%occ(iln,1)
   enddo
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

subroutine rcpaw_core_energies(rcpaw,ntypat)
!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ntypat
 type(rcpaw_type), pointer, intent(inout) :: rcpaw
!arrays

!Local variables-------------------------------
!scalars
 integer :: itypat

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
 write(std_out,*)'RCPAW energies at step ',rcpaw%istep,':', ' ekinc = ',rcpaw%ekinc,' eeig = ',&
& rcpaw%eeigc,' edcc = ',rcpaw%edcc,' ehnzc = ',rcpaw%ehnzc

end subroutine rcpaw_core_energies
!!***


end module m_rcpaw
!!***

