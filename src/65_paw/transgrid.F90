!{\src2tex{textfont=tt}}
!!****f* ABINIT/transgrid
!! NAME
!! transgrid
!!
!! FUNCTION
!! Convert a given density (or potential) from the coarse to the fine rectangular grid and vice versa
!! Used in PAW calculations
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cplex=1 if rhor[f] is real, 2 if rhor[f] is complex
!!  mpi_enreg=informations about MPI parallelization
!!  nspden=number of spin-density components
!!  optgrid=+1 to go from the coarse grid towards the fine grid
!!          -1 to go from the fine grid towards the coarse grid
!!  optin= 0: input density/potential is taken from rhor(:,nspden)
!!         1: input density/potential is taken from rhog(:)     (ispden=1)
!!                                              and rhor(:,2:4) (ispden=2,3,4)
!!  optout= 0: output density/potential is given in r space in rhor(:,nspden)
!!          1: output density/potential is given in r space in rhor(:,nspden)
!!                                           and in g space in rhog(:)
!!  pawfgr <type(paw_fgr_type)>=fine rectangular grid parameters
!!    %nfftc=number of points in the coarse FFT box
!!    %nfft =number of points in the fine FFT box
!!    %ngfftc(18)=all needed information about 3D FFT, for the coarse grid
!!    %ngfft(18) =all needed information about 3D FFT, for the fine grid
!!    %coatofin(nfftc)=Index of the points of the coarse grid on the fine grid
!!    %fintocoa(nfft) =Index of the points of the fine grid on the coarse grid
!!    %usefinegrid= 1 if a fine FFT grid is used (0 otherwise)
!!  if optgrid=+1 and optin=1:
!!    rhog(2,nfftc)=Fourier transform of input density/potential on the coarse grid
!!  if optgrid=-1 and optin=1:
!!    rhogf(2,nfftf)=Fourier transform of input density/potential on the fine grid
!!  if optgrid=+1
!!    rhor(cplex*nfftc,nspden)=input density/potential in r space on the coarse grid
!!  if optgrid=-1:
!!    rhorf(cplex*nfftf,nspden)=input density/potential in r space on the fine grid
!!
!! OUTPUT
!!  if optgrid=-1 and optout=1:
!!    rhog(2,nfftc)=Fourier transform of output density/potential on the coarse grid
!!  if optgrid=+1 and optout=1:
!!    rhogf(2,nfftf)=Fourier transform of output density/potential on the fine grid
!!  if optgrid=-1
!!    rhor(cplex*nfftc,nspden)=output density/potential in r space on the coarse grid
!!  if optgrid=+1:
!!    rhorf(cplex*nfftf,nspden)=output density/potential in r space on the fine grid
!!
!! PARENTS
!!      dfpt_looppert,energy,fourier_interpol,getgh1c,gstate,ks_ddiago,m_io_kss
!!      pawmkrho,respfn,vtorho,vtorhorec
!!
!! CHILDREN
!!      fourdp,indirect_parallel_fourier,zerosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine transgrid(cplex,mpi_enreg,nspden,optgrid,optin,optout,paral_kgb,pawfgr,rhog,rhogf,rhor,rhorf)

 use m_profiling_abi
 use defs_basis
 use defs_abitypes
 use m_errors

 use m_pawfgr, only : pawfgr_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'transgrid'
 use interfaces_52_fft_mpi_noabirule
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,nspden,optgrid,optin,optout,paral_kgb
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawfgr_type),intent(in) :: pawfgr
!arrays
 real(dp),intent(inout) :: rhog(2,pawfgr%nfftc),rhogf(2,pawfgr%nfft)
 real(dp),intent(inout) :: rhor(cplex*pawfgr%nfftc,nspden),rhorf(cplex*pawfgr%nfft,nspden)

!Local variables ---------------------------------------
!scalars
 integer :: i1,ispden,nfftc,nfftctot,nfftf,nfftftot
 character(len=500) :: msg
!arrays
 integer :: ngfftc(18),ngfftf(18)
 real(dp),allocatable :: vectg(:,:),work(:,:),workfft(:)

! *************************************************************************

 DBG_ENTER("COLL")

!Tests
 if(pawfgr%nfft<pawfgr%nfftc) then
   write(msg,'(a,2(i0,1x))')' nfft (fine grid) must be >= nfft (coarse grid) while: ',pawfgr%nfft, pawfgr%nfftc
   MSG_ERROR(msg)
 end if

!Store FFT dimensions
 nfftc=pawfgr%nfftc;ngfftc(:)=pawfgr%ngfftc(:);nfftctot=ngfftc(1)*ngfftc(2)*ngfftc(3)
 nfftf=pawfgr%nfft ;ngfftf(:)=pawfgr%ngfft (:);nfftftot=ngfftf(1)*ngfftf(2)*ngfftf(3)

!If no fine FFT grid is used, this is only a simple transfer
 if (pawfgr%usefinegrid==0) then
   if (optgrid==1) then
     rhorf=rhor
     if (optout==1.and.optin==1) rhogf=rhog
     if (optout==1.and.optin/=1) then
       ABI_ALLOCATE(workfft,(cplex*nfftc))
       workfft(:)=rhor(:,1)
       call fourdp(cplex,rhogf,workfft,-1,mpi_enreg,nfftc,ngfftc,paral_kgb,0)
       ABI_DEALLOCATE(workfft)
     end if
   end if
   if (optgrid==-1) then
     rhor=rhorf
     if (optout==1.and.optin==1) rhog=rhogf
     if (optout==1.and.optin/=1) then
       ABI_ALLOCATE(workfft,(cplex*nfftc))
       workfft(:)=rhorf(:,1)
       call fourdp(cplex,rhog,workfft,-1,mpi_enreg,nfftc,ngfftc,paral_kgb,0)
       ABI_DEALLOCATE(workfft)
     end if
   end if
   return
 end if

!====== FROM THE COARSE GRID TOWARDS THE FINE GRID =============
!===============================================================
!Calculate the FT of rhor to have it in the g space on the coarse grid
!Transfer the FT of rhor on the coarse grid towards the fine grid
!Then calculate the FT back to get rhorf on the fine grid
 if (optgrid==1) then

   ABI_ALLOCATE(work,(2,nfftc))

!  First spin component
!  --------------------------------------------------------------
   if (optout==0) then
!    if optout=0, rhog on the fine grid is temporary (in vectg)
     ABI_ALLOCATE(vectg,(2,nfftf))
     vectg(:,:)=zero
     if (optin==1) then
       call zerosym(rhog,2,ngfftc(1),ngfftc(2),ngfftc(3),&
&       comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
       if(mpi_enreg%nproc_fft > 1 .and. mpi_enreg%paral_kgb==1) then
         call indirect_parallel_Fourier&
&         (pawfgr%coatofin,vectg,mpi_enreg,ngfftf,ngfftc,nfftf,nfftc,paral_kgb,rhog,nfftctot)
       else
         do i1=1,nfftc
           vectg(:,pawfgr%coatofin(i1))=rhog(:,i1)
         end do
       end if
     else
       ABI_ALLOCATE(workfft,(cplex*nfftc))
       workfft(:)=rhor(:,1)
       call fourdp(cplex,work,workfft,-1,mpi_enreg,nfftc,ngfftc,paral_kgb,0)
       ABI_DEALLOCATE(workfft)
       call zerosym(work,2,ngfftc(1),ngfftc(2),ngfftc(3),&
&       comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
       if(mpi_enreg%nproc_fft > 1 .and. mpi_enreg%paral_kgb==1) then
         call indirect_parallel_Fourier&
&         (pawfgr%coatofin,vectg,mpi_enreg,ngfftf,ngfftc,nfftf,nfftc,paral_kgb,work,nfftctot)
       else
         do i1=1,nfftc
           vectg(:,pawfgr%coatofin(i1))=work(:,i1)
         end do
       end if
     end if
!    call zerosym(vectg,2,ngfftf(1),ngfftf(2),ngfftf(3),&
!    &        comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
     ABI_ALLOCATE(workfft,(cplex*nfftf))
     call fourdp(cplex,vectg,workfft,1,mpi_enreg,nfftf,ngfftf,paral_kgb,0)
     rhorf(:,1)=workfft(:)
     ABI_DEALLOCATE(workfft)
     ABI_DEALLOCATE(vectg)
   else
!    if optout=1, rhog on the fine grid is saved
     call zerosym(rhog,2,ngfftc(1),ngfftc(2),ngfftc(3),&
&     comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
     rhogf(:,:)=zero
     if (optin==1) then
       if(mpi_enreg%nproc_fft > 1 .and. mpi_enreg%paral_kgb==1) then
         call indirect_parallel_Fourier&
&         (pawfgr%coatofin,rhogf,mpi_enreg,ngfftf,ngfftc,nfftf,nfftc,paral_kgb,rhog,nfftctot)
       else
         do i1=1,nfftc
           rhogf(:,pawfgr%coatofin(i1))=rhog(:,i1)
         end do
       end if
     else
       ABI_ALLOCATE(workfft,(cplex*nfftc))
       workfft(:)=rhor(:,1)
       call fourdp(cplex,work,workfft,-1,mpi_enreg,nfftc,ngfftc,paral_kgb,0)
       ABI_DEALLOCATE(workfft)
       call zerosym(work,2,ngfftc(1),ngfftc(2),ngfftc(3),&
&       comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
       if(mpi_enreg%nproc_fft > 1 .and. mpi_enreg%paral_kgb==1) then
         call indirect_parallel_Fourier&
&         (pawfgr%coatofin,rhogf,mpi_enreg,ngfftf,ngfftc,nfftf,nfftc,paral_kgb,work,nfftctot)
       else
         do i1=1,nfftc
           rhogf(:,pawfgr%coatofin(i1))=work(:,i1)
         end do
       end if
     end if
!    call zerosym(rhogf,2,ngfftf(1),ngfftf(2),ngfftf(3),&
!    &        comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
     ABI_ALLOCATE(workfft,(cplex*nfftf))
     call fourdp(cplex,rhogf,workfft,1,mpi_enreg,nfftf,ngfftf,paral_kgb,0)
     rhorf(:,1)=workfft(:)
     ABI_DEALLOCATE(workfft)
   end if

!  Additional spin components
!  ----------------------------------------------------
   if (nspden>=2) then
     ABI_ALLOCATE(vectg,(2,nfftf))
     do ispden=2,nspden
       vectg(:,:)=zero
       ABI_ALLOCATE(workfft,(cplex*nfftc))
       workfft(:)=rhor(:,ispden)
       call fourdp(cplex,work,workfft,-1,mpi_enreg,nfftc,ngfftc,paral_kgb,0)
       ABI_DEALLOCATE(workfft)
       call zerosym(work,2,ngfftc(1),ngfftc(2),ngfftc(3),&
&       comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
       if(mpi_enreg%nproc_fft > 1 .and. mpi_enreg%paral_kgb==1) then
         call indirect_parallel_Fourier&
&         (pawfgr%coatofin,vectg,mpi_enreg,ngfftf,ngfftc,nfftf,nfftc,paral_kgb,work,nfftctot)
       else
         do i1=1,nfftc
           vectg(:,pawfgr%coatofin(i1))=work(:,i1)
         end do
       end if
!      call zerosym(vectg,2,ngfftf(1),ngfftf(2),ngfftf(3),&
!      &          comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
       ABI_ALLOCATE(workfft,(cplex*nfftf))
       call fourdp(cplex,vectg,workfft,1,mpi_enreg,nfftf,ngfftf,paral_kgb,0)
       rhorf(:,ispden)=workfft(:)
       ABI_DEALLOCATE(workfft)
     end do
     ABI_DEALLOCATE(vectg)
   end if

   ABI_DEALLOCATE(work)


!  ====== FROM THE FINE GRID TOWARDS THE COARSE GRID =============
!  ==============================================================
!  Calculate the FT of rhorf to have it in the g space on the fine grid
!  Transfer the FT of rhorf on the fine grid towards the coarse grid
!  Then calculate the FT back to get rhor on the coarse grid
 else if (optgrid==-1) then

   ABI_ALLOCATE(work,(2,nfftf))

!  First spin component
!  --------------------------------------------------------------
   if (optout==0) then
!    if optout=0, rhog on the fine grid is temporary (in vectg)
     ABI_ALLOCATE(vectg,(2,nfftc))
     vectg(:,:)=zero
     if (optin==1) then
       do i1=1,nfftf
         if (pawfgr%fintocoa(i1)/=0) vectg(:,pawfgr%fintocoa(i1))=rhogf(:,i1)
       end do
     else
       ABI_ALLOCATE(workfft,(cplex*nfftf))
       workfft(:)=rhorf(:,1)
       call fourdp(cplex,work,workfft,-1,mpi_enreg,nfftf,ngfftf,paral_kgb,0)
       ABI_DEALLOCATE(workfft)
       if(mpi_enreg%nproc_fft > 1 .and. mpi_enreg%paral_kgb==1) then
         call indirect_parallel_Fourier&
&         (pawfgr%fintocoa,vectg,mpi_enreg,ngfftc,ngfftf,nfftc,nfftf,paral_kgb,work,nfftftot)
       else
         do i1=1,nfftf
           if (pawfgr%fintocoa(i1)/=0) vectg(:,pawfgr%fintocoa(i1))=work(:,i1)
         end do
       end if
     end if
     call zerosym(vectg,2,ngfftc(1),ngfftc(2),ngfftc(3),&
&     comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
     ABI_ALLOCATE(workfft,(cplex*nfftc))
     call fourdp(cplex,vectg,workfft,1,mpi_enreg,nfftc,ngfftc,paral_kgb,0)
     rhor(:,1)=workfft(:)
     ABI_DEALLOCATE(workfft)
     ABI_DEALLOCATE(vectg)
   else
!    if optout=1, rhog on the fine grid is saved
     rhog(:,:)=zero
     if (optin==1) then
       do i1=1,nfftf
         if (pawfgr%fintocoa(i1)/=0) rhog(:,pawfgr%fintocoa(i1))=rhogf(:,i1)
       end do
     else
       ABI_ALLOCATE(workfft,(cplex*nfftf))
       workfft(:)=rhorf(:,1)
       call fourdp(cplex,work,workfft,-1,mpi_enreg,nfftf,ngfftf,paral_kgb,0)
       ABI_DEALLOCATE(workfft)
       if(mpi_enreg%nproc_fft > 1 .and. mpi_enreg%paral_kgb==1) then
         call indirect_parallel_Fourier&
&         (pawfgr%fintocoa,rhog,mpi_enreg,ngfftc,ngfftf,nfftc,nfftf,paral_kgb,work,nfftftot)
       else
         do i1=1,nfftf
           if (pawfgr%fintocoa(i1)/=0) rhog(:,pawfgr%fintocoa(i1))=work(:,i1)
         end do
       end if
     end if
     call zerosym(rhog,2,ngfftc(1),ngfftc(2),ngfftc(3),&
&     comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
     ABI_ALLOCATE(workfft,(cplex*nfftc))
     call fourdp(cplex,rhog,workfft,1,mpi_enreg,nfftc,ngfftc,paral_kgb,0)
     rhor(:,1)=workfft(:)
     ABI_DEALLOCATE(workfft)
   end if

!  Additional spin components
!  ----------------------------------------------------
   if (nspden>=2) then
     ABI_ALLOCATE(vectg,(2,nfftc))
     do ispden=2,nspden
       vectg(:,:)=zero
       ABI_ALLOCATE(workfft,(cplex*nfftf))
       workfft(:)=rhorf(:,ispden)
       call fourdp(cplex,work,workfft,-1,mpi_enreg,nfftf,ngfftf,paral_kgb,0)
       ABI_DEALLOCATE(workfft)
       if(mpi_enreg%nproc_fft > 1 .and. mpi_enreg%paral_kgb==1) then
         call indirect_parallel_Fourier&
&         (pawfgr%fintocoa,vectg,mpi_enreg,ngfftc,ngfftf,nfftc,nfftf,paral_kgb,work,nfftftot)
       else
         do i1=1,nfftf
           if (pawfgr%fintocoa(i1)/=0) vectg(:,pawfgr%fintocoa(i1))=work(:,i1)
         end do
       end if
       call zerosym(vectg,2,ngfftc(1),ngfftc(2),ngfftc(3),&
&       comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
       ABI_ALLOCATE(workfft,(cplex*nfftc))
       call fourdp(cplex,vectg,workfft,1,mpi_enreg,nfftc,ngfftc,paral_kgb,0)
       rhor(:,ispden)=workfft(:)
       ABI_DEALLOCATE(workfft)
     end do
     ABI_DEALLOCATE(vectg)
   end if

   ABI_DEALLOCATE(work)

 end if

 DBG_EXIT("COLL")

end subroutine transgrid
!!***
