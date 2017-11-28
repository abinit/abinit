!{\src2tex{textfont=tt}}
!!****f* ABINIT/printmagvtk
!! NAME
!!  printmagvtk
!!
!! FUNCTION
!!  Auxiliary routine for printing out magnetization density in VTK format.
!!  Output file name is DEN.vtk
!!
!! COPYRIGHT
!!  Copyright (C) 2017 ABINIT group (SPr)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  mpi_enreg = information about adopted parallelization strategy
!!  nspden    = number of components of density matrix (possible values re 1,2, or 4)
!!              nspden:   1 -> rho
!!              nspden:   2 -> rho_up,rho_dwn
!!              nspden:   4 -> rho,mx,my,mz
!!  nfft      = number of fft points per FFT processor
!!  ngfft     = full information about FFT mesh
!!  rhor      = density array stored in the memory of current FFT processor
!!  rprimd    = array of lattice vectors
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  At the moment this routine is mainly used for development and debugging 
!!  of gs and dfpt calculations with non-collinear spins. If needed, can be used
!!  to print final density in vtk format.
!!  IMPORTANT: implementation is thoroughly checked only for npspinor = 1,
!!             for other case might need to change the part gathering
!!             the FFT mesh info
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_gather
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine printmagvtk(mpi_enreg,cplex,nspden,nfft,ngfft,rhor,rprimd,fname)
    
 use defs_basis
 use m_errors
 use m_profiling_abi
 use defs_abitypes
 use m_xmpi
 use m_io_tools,        only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'printmagvtk'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(in)   :: mpi_enreg
 integer,intent(in)          :: nfft,nspden,cplex
!arrays
 integer,intent(in)          :: ngfft(18)
 real(dp),intent(in)         :: rhor(cplex*nfft,nspden),rprimd(3,3)
 character(len=*),intent(in) :: fname

!Local variables-------------------------------
!scalars
 integer :: denvtk,denxyz,nfields
 integer :: nx,ny,nz,nfft_tot
 integer :: ii,jj,kk,ind,jfft,ispden
 integer :: mpi_comm,mpi_head,mpi_rank,ierr
 real    :: rx,ry,rz
 integer :: nproc_fft,ir
 character(len=500) :: msg
 character(len=10)  :: outformat
!arrays
 real(dp),allocatable :: rhorfull(:,:)

 
! *************************************************************************

 DBG_ENTER("COLL")
 
! if (option/=1 .and. option/=2 ) then
!   write(msg,'(3a,i0)')&
!&   'The argument option should be 1 or 2,',ch10,&
!&   'however, option=',option
!   MSG_BUG(msg)
! end if
!
! if (sizein<1) then
!   write(msg,'(3a,i0)')&
!&   'The argument sizein should be a positive number,',ch10,&
!&   'however, sizein=',sizein
!   MSG_ERROR(msg)
! end if

 DBG_EXIT("COLL")

!write(std_out,*) ' Writing out .vtk file: ',fname

  !if 1 or two component density then write out either 1 or 2 scalar density fields
  !if 4, then write one scalar field (density) and one vector field (magnetization density)
 if(nspden/=4)then
   nfields=nspden
 else
   nfields=2
 end if

 nfields=nfields*cplex

  ! FFT mesh specifications: full grid
 nx=ngfft(1)        ! number of points along 1st lattice vector
 ny=ngfft(2)        ! number of points along 2nd lattice vector
 nz=ngfft(3)        ! number of points along 3rd lattice vector
 nfft_tot=nx*ny*nz  ! total number of fft mesh points (can be different from nfft in case of distributed memory of nproc_fft processors)
 

  ! Gather information about memory distribution
 mpi_head=0
 mpi_comm = mpi_enreg%comm_fft
 mpi_rank = xmpi_comm_rank(mpi_comm)
 nproc_fft=ngfft(10)

  ! Create array to host full FFT mesh
 if(mpi_rank==mpi_head)then
   ABI_ALLOCATE(rhorfull,(cplex*nfft_tot,nspden))
 end if

  ! Fill in the full mesh
 if(nproc_fft==1)then
   rhorfull=rhor
 else
   do ir=1,nspden
     call xmpi_gather(rhor(:,ir),cplex*nfft,rhorfull(:,ir),cplex*nfft,mpi_head,mpi_comm,ierr)
   end do
 end if

 if(mpi_rank==mpi_head)then

    ! Open the output vtk file
   if (open_file(fname,msg,newunit=denvtk,status='replace',form='formatted') /=0) then
     MSG_WARNING(msg)
     RETURN
   end if

   if (open_file("DEN.xyz",msg,newunit=denxyz,status='replace',form='formatted') /=0) then
     MSG_WARNING(msg)
     RETURN
   end if

    ! Write the header of the output vtk file
   write(denvtk,"(a)") '# vtk DataFile Version 2.0'
   write(denvtk,"(a)") 'Electron density components'
   write(denvtk,"(a)") 'ASCII'
   write(denvtk,"(a)") 'DATASET STRUCTURED_GRID'
   write(denvtk,"(a,3i6)") 'DIMENSIONS ', nx,ny,nz
   write(denvtk,"(a,i18,a)") 'POINTS ',nfft_tot,' double'

   if (nspden==1) then
     outformat="(4e16.8)"
   else if (nspden==2) then
     outformat="(5e16.8)"
   else
     outformat="(7e16.8)"
   endif

    ! Write out information about grid points
   do kk=0,nz-1
     do jj=0,ny-1
       do ii=0,nx-1

         rx=(dble(ii)/nx)*rprimd(1,1)+(dble(jj)/ny)*rprimd(1,2)+(dble(kk)/nz)*rprimd(1,3)
         ry=(dble(ii)/nx)*rprimd(2,1)+(dble(jj)/ny)*rprimd(2,2)+(dble(kk)/nz)*rprimd(2,3)
         rz=(dble(ii)/nx)*rprimd(3,1)+(dble(jj)/ny)*rprimd(3,2)+(dble(kk)/nz)*rprimd(3,3)
         write(denvtk,'(3f16.8)') rx,ry,rz  !coordinates of the grid point
         ind=1+ii+nx*(jj+ny*kk)
         write(denxyz,outformat) rx,ry,rz,(rhorfull(ind,ispden),ispden=1,nspden)
       end do
     end do
   end do

   close(denxyz)

    ! Write out information about field defined on the FFT mesh
   write(denvtk,"(a,i18)") 'POINT_DATA ',nfft_tot
   write(denvtk,"(a,i6)")  'FIELD Densities ',nfields


    ! Write out different fields depending on the number of density matrix components
   if(nspden==1)then
     
      !single component, so just write out the density
     if(cplex==1) then
       write(denvtk,"(a,i18,a)") 'rho 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=1+ii+nx*(jj+ny*kk)
             write(denvtk,'(f16.8)') rhorfull(ind,1)
           end do
         end do
       end do
     else
       write(denvtk,"(a,i18,a)") 'Re_rho 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))-1
             write(denvtk,'(f16.8)') rhorfull(ind,1)
           end do
         end do
       end do
       write(denvtk,"(a,i18,a)") 'Im_rho 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))
             write(denvtk,'(f16.8)') rhorfull(ind,1)
           end do
         end do
       end do
     end if

   else if(nspden==2)then

      !two component, write the density for spin_up and spin_down channels
     if(cplex==1) then
       
       write(denvtk,"(a,i18,a)") 'rho 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=1+ii+nx*(jj+ny*kk)
             write(denvtk,'(f16.8)') rhorfull(ind,1)
           end do
         end do
       end do
       write(denvtk,"(a,i18,a)") 'mag 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=1+ii+nx*(jj+ny*kk)
             write(denvtk,'(f16.8)') 2*rhorfull(ind,2)-rhorfull(ind,1)
           end do
         end do
       end do
       
     else
       
       write(denvtk,"(a,i18,a)") 'Re_rho 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))-1
             write(denvtk,'(f16.8)') rhorfull(ind,1)
           end do
         end do
       end do
       write(denvtk,"(a,i18,a)") 'Im_rho 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))
             write(denvtk,'(f16.8)') rhorfull(ind,1)
           end do
         end do
       end do
       write(denvtk,"(a,i18,a)") 'Re_mag 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))-1
             write(denvtk,'(f16.8)') 2*rhorfull(ind,2)-rhorfull(ind,1)
           end do
         end do
       end do
       write(denvtk,"(a,i18,a)") 'Im_mag 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))
             write(denvtk,'(f16.8)') 2*rhorfull(ind,2)-rhorfull(ind,1)
           end do
         end do
       end do
       
     end if   

   else  !here is the last option: nspden==4

     if(cplex==1) then
       
        !four component, write the density (scalar field) and magnetization density (vector field)
       write(denvtk,"(a,i18,a)") 'rho 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=1+ii+nx*(jj+ny*kk)
             write(denvtk,'(f16.8)') rhorfull(ind,1)
           end do
         end do
       end do
       write(denvtk,"(a,i18,a)") 'mag 3 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=1+ii+nx*(jj+ny*kk)
             write(denvtk,'(3f16.8)') rhorfull(ind,2),rhorfull(ind,3),rhorfull(ind,4)
           end do
         end do
       end do
       
     else
       
       write(denvtk,"(a,i18,a)") 'Re_rho 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))-1
             write(denvtk,'(f16.8)') rhorfull(ind,1)
           end do
         end do
       end do
       write(denvtk,"(a,i18,a)") 'Im_rho 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))
             write(denvtk,'(f16.8)') rhorfull(ind,1)
           end do
         end do
       end do
       write(denvtk,"(a,i18,a)") 'Re_mag 3 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))-1
             write(denvtk,'(3f16.8)') rhorfull(ind,2),rhorfull(ind,3),rhorfull(ind,4)
           end do
         end do
       end do
       write(denvtk,"(a,i18,a)") 'Im_mag 3 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))
             write(denvtk,'(3f16.8)') rhorfull(ind,2),rhorfull(ind,3),rhorfull(ind,4)
           end do
         end do
       end do
       
     end if

   end if ! nspden options condition

   close (denvtk)

    !clean up the gathered FFT mesh
   ABI_DEALLOCATE(rhorfull)

 end if


end subroutine printmagvtk
!!***
