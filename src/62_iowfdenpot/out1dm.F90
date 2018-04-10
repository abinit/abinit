!{\src2tex{textfont=tt}}
!!****f* ABINIT/out1dm
!! NAME
!! out1dm
!!
!! FUNCTION
!! Output the 1 dimensional mean of potential and density
!! on the three reduced axis.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  fnameabo_app_1dm=name of the file in which the data is written, appended with _1DM
!!  natom=number of atoms in unit cell
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms in unit cell.
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=real space dimensional primitive translations (bohr)
!!  typat(natom)=type integer for each atom in cell
!!  ucvol=unit cell volume in bohr**3.
!!  vtrial(nfft,nspden)=INPUT Vtrial(r).
!!  xred(3,natom)=reduced coordinates of atoms
!!  znucl(ntypat)=real(dp), nuclear number of atom type
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  data written in file fnameabo_app_1dm
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      atomdata_from_znucl,ptabs_fourdp,wrtout,xmpi_sum,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine out1dm(fnameabo_app_1dm,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,&
&  rhor,rprimd,typat,ucvol,vtrial,xred,znucl)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_atomdata

 use m_io_tools, only : open_file
 use m_geometry,  only : xred2xcart
 use m_mpinfo,   only : ptabs_fourdp
 use m_xmpi,     only : xmpi_sum, xmpi_comm_size

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'out1dm'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,nspden,ntypat
 real(dp),intent(in) :: ucvol
 character(len=fnlen),intent(in) :: fnameabo_app_1dm
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18),typat(natom)
 real(dp),intent(in) :: rhor(nfft,nspden),rprimd(3,3),vtrial(nfft,nspden)
 real(dp),intent(in) :: znucl(ntypat),xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: ia,ib,idim,ierr,ifft,islice,ispden,na,nb,ndig,nslice,nu,temp_unit
 integer :: comm_fft, nproc_fft, me_fft
 real(dp) :: global_den,global_pot
 character(len=2) :: symbol
 character(len=500) :: message
 type(atomdata_t) :: atom
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp),allocatable :: lin_den(:),mean_pot(:),reduced_coord(:),xcart(:,:)
 character(len=8),allocatable :: iden(:)

! *************************************************************************

!Initialize the file
 write(message, '(a,a)' ) ' io1dm : about to open file ',fnameabo_app_1dm
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 comm_fft = mpi_enreg%comm_fft; nproc_fft = xmpi_comm_size(comm_fft); me_fft = mpi_enreg%me_fft

 if (me_fft == 0) then
   if (open_file(fnameabo_app_1dm,message,newunit=temp_unit,status='unknown',form='formatted') /= 0) then
     MSG_ERROR(message)
   end if
   rewind(temp_unit)
 end if

 write(message, '(a,a)' ) ch10,'# ABINIT package : 1DM file '
 if (me_fft == 0) write(temp_unit,'(a)') message

 write(message, '(a,a)' )ch10,'# Primitive vectors of the periodic cell (bohr)'
 if (me_fft == 0) write(temp_unit,'(a)') message
 do nu=1,3
   write(message, '(1x,a,i1,a,3f10.5)' ) '#  R(',nu,')=',rprimd(:,nu)
   if (me_fft == 0) write(temp_unit,'(a)') message
 end do

 write(message, '(a,a)' ) ch10,&
& '# Atom list        Reduced coordinates          Cartesian coordinates (bohr)'
 if (me_fft == 0) write(temp_unit,'(a)') message

!Set up a list of character identifiers for all atoms : iden(ia)
 ABI_ALLOCATE(iden,(natom))
 do ia=1,natom
   call atomdata_from_znucl(atom,znucl(typat(ia)))
   symbol = atom%symbol

   ndig=int(log10(dble(ia)+0.5_dp))+1
   if(ndig==1) write(iden(ia), '(a,a,i1,a)' )symbol,'(',ia,')   '
   if(ndig==2) write(iden(ia), '(a,a,i1,a)' )symbol,'(',ia,')  '
   if(ndig==3) write(iden(ia), '(a,a,i1,a)' )symbol,'(',ia,') '
   if(ndig==4) write(iden(ia), '(a,a,i1,a)' )symbol,'(',ia,')'
   if(ndig>4)then
     close(temp_unit)
     write(message, '(a,i0)' )&
&     ' out1dm : cannot handle more than 9999 atoms, while natom=',natom
     MSG_ERROR(message)
   end if
 end do

!Compute cartesian coordinates, and print reduced and cartesian coordinates
 ABI_ALLOCATE(xcart,(3,natom))
 call xred2xcart(natom,rprimd,xcart,xred)
 do ia=1,natom
   write(message, '(a,a,3f10.5,a,3f10.5)' ) &
&   '#   ',iden(ia),xred(1:3,ia),'    ',xcart(1:3,ia)
   if (me_fft == 0) write(temp_unit,'(a)') message
 end do
 ABI_DEALLOCATE(iden)
 ABI_DEALLOCATE(xcart)

!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,ngfft(2),ngfft(3),fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

 do idim=1,3

   nslice=ngfft(idim)

!  Dummy initialisation of na, nb and ifft.
   na=0 ; nb=0; ifft=0
   select case(idim)
   case(1)
     na=ngfft(2) ; nb=ngfft(3)
   case(2)
     na=ngfft(1) ; nb=ngfft(3)
   case(3)
     na=ngfft(1) ; nb=ngfft(2)
   end select

   ABI_ALLOCATE( reduced_coord,(nslice))
   ABI_ALLOCATE(mean_pot,(nslice))
   ABI_ALLOCATE(lin_den,(nslice))

   do ispden=1,nspden

     if(ispden==1)then
       write(message, '(a,a,a)' ) ch10,'#===========',&
&       '====================================================================='
       if (me_fft == 0) write(temp_unit,'(a)') message
     end if

     select case(idim)
     case(1)
       write(message, '(a)' )'# Projection along the first dimension '
     case(2)
       write(message, '(a)' )'# Projection along the second dimension '
     case(3)
       write(message, '(a)' )'# Projection along the third dimension '
     end select
     if (me_fft == 0) write(temp_unit,'(a)') message

     if(nspden==2)then
       select case(ispden)
       case(1)
         write(message, '(a)' )'# Spin up '
       case(2)
         write(message, '(a)' )'# Spin down '
       end select
       if (me_fft == 0) write(temp_unit,'(a)') message
     else if (nspden == 4) then
       select case(ispden)
       case(1)
         write(message, '(a)' )'# Spinor up up '
       case(2)
         write(message, '(a)' )'# Spinor down down'
       case(3)
         write(message, '(a)' )'# Spinor up down'
       case(4)
         write(message, '(a)' )'# Spinor down up'
       end select
       if (me_fft == 0) write(temp_unit,'(a)') message
     end if

     write(message, '(2a)' ) ch10,&
&     '#     Red. coord. Mean KS potential  Linear density  '
     if (me_fft == 0) write(temp_unit,'(a)') message

     write(message, '(a)' ) &
&     '#                  (Hartree unit)   (electron/red. unit)'
     if (me_fft == 0) write(temp_unit,'(a)') message

     global_pot=zero
     global_den=zero
     mean_pot(:)=zero
     lin_den(:)=zero
     do islice=1,nslice
       reduced_coord(islice)=(islice-1)/dble(nslice)
     end do
     if (idim == 1) then
       do islice=1,nslice
         do ib=1,nb
! skip z values not on this processor
           if (fftn3_distrib(ib)/=mpi_enreg%me_fft) cycle
           do ia=1,na
             ifft = islice + nslice*( ia    -1 + na    *(ffti3_local(ib)    -1) )
             mean_pot(islice)=mean_pot(islice)+vtrial(ifft,ispden)
             lin_den(islice)=lin_den(islice)+rhor(ifft,ispden)
           end do
         end do
       end do
     else if (idim == 2) then
       do islice=1,nslice
         do ib=1,nb
! skip z values not on this processor
           if (fftn3_distrib(ib)/=mpi_enreg%me_fft) cycle
           do ia=1,na
             ifft = ia     + na    *( islice-1 + nslice*(ffti3_local(ib)    -1) )
             mean_pot(islice)=mean_pot(islice)+vtrial(ifft,ispden)
             lin_den(islice)=lin_den(islice)+rhor(ifft,ispden)
           end do
         end do
       end do
     else if (idim == 3) then
! this full z slice is mine in parallel case
       do islice=1,nslice
         if (fftn3_distrib(islice)/=mpi_enreg%me_fft) cycle
         do ib=1,nb
           do ia=1,na
             ifft = ia     + na    *( ib    -1 + nb    *(ffti3_local(islice)-1) )
             mean_pot(islice)=mean_pot(islice)+vtrial(ifft,ispden)
             lin_den(islice)=lin_den(islice)+rhor(ifft,ispden)
           end do
         end do
       end do
     end if
     mean_pot(:)=mean_pot(:)/dble(na*nb)
     lin_den(:)=lin_den(:)/dble(na*nb)*ucvol
     global_pot=sum(mean_pot)/dble(nslice)
     global_den=sum(lin_den)/dble(nslice)

!    FFT parallelization
     if(mpi_enreg%nproc_fft>1)then
       call xmpi_sum(mean_pot  ,mpi_enreg%comm_fft,ierr)
       call xmpi_sum(global_pot,mpi_enreg%comm_fft,ierr)
       call xmpi_sum(lin_den   ,mpi_enreg%comm_fft,ierr)
       call xmpi_sum(global_den,mpi_enreg%comm_fft,ierr)
     end if

     do islice=1,ngfft(idim)
       write(message, '(f10.4,es20.6,es16.6)' )&
&       reduced_coord(islice),mean_pot(islice),lin_den(islice)
       if (me_fft == 0) write(temp_unit,'(a)') message
     end do

     write(message, '(a,a,es15.6,es16.6,a)' ) ch10,&
&     '# Cell mean       :',global_pot,global_den, ch10
     if (me_fft == 0) write(temp_unit,'(a)') message


!    End of the loop on spins
   end do

   ABI_DEALLOCATE(reduced_coord)
   ABI_DEALLOCATE(mean_pot)
   ABI_DEALLOCATE(lin_den)

!  End of the loops on the three dimensions
 end do

 if (me_fft == 0) then
   write(message, '(a,a,a)' ) ch10,'#===========',&
&   '====================================================================='
   call wrtout(temp_unit,message,'COLL')
   close(temp_unit)
 end if

end subroutine out1dm
!!***
