!!****m* ABINIT/m_multipoles
!! NAME
!!  m_multipoles
!!
!! FUNCTION
!!  Compute spatial multipole moments of input array on FFT grid
!!
!! COPYRIGHT
!!  Copyright (C) 2003-2020 ABINIT group (MJV, MT, XG)
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

module m_multipoles

 use defs_basis
 use m_errors
 use m_abicore
 use m_distribfft
 use m_xmpi
 use m_atomdata

 use defs_abitypes,    only : mpi_type
 use m_io_tools, only : open_file
 use m_geometry, only : xred2xcart
 use m_mpinfo,   only : ptabs_fourdp

 implicit none

 private
!!***

 public :: multipoles_out    ! Output multipole moments of input array on FFT grid, calculated with multipoles_fftr
 public :: out1dm            ! Output the 1 dimensional mean of potential and density on the three reduced axis.
!!***

contains
!!***

!!****f* m_multipoles/multipoles_fftr
!! NAME
!! multipoles_fftr
!!
!! FUNCTION
!!  Compute spatial multipole moments of input array on FFT grid
!!  Namely, the electrical dipole, quadrupole, etc... of the electron density
!!
!! INPUTS
!!  arraysp(nfft,nspden)=the array whose average has to be computed
!!  [distribfft<type(distribfft_type)>]= -optional- infos related to FFT parallelism
!!  [mpi_comm_grid]= -optional- MPI communicator over the grid
!!  origin(3)=vector defining the origin of the dipole (point of observation, in reduced coordinates)
!!  nfft=number of FFT points stored by one proc
!!  nfftot=total number of FFT points
!!  ngfft =number of subdivisions along each lattice vector
!!  nspden=number of spin-density components
!!  rprimd = dimensionful lattice vectors
!!
!! OUTPUT
!!  dipole(nspden)=mean value of the dipole of input array, for each nspden component
!!
!! PARENTS
!!      m_multipoles
!!
!! CHILDREN
!!      atomdata_from_znucl,ptabs_fourdp,wrtout,xmpi_sum,xred2xcart
!!
!! SOURCE

subroutine multipoles_fftr(arraysp,dipole,nfft,ngfft,nspden,rprimd,origin,&
                           distribfft,mpi_comm_grid)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nspden
 integer,intent(in),optional :: mpi_comm_grid
 type(distribfft_type),intent(in),optional,target :: distribfft
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: arraysp(nfft,nspden),origin(3),rprimd(3,3)
 real(dp),intent(out) :: dipole(3,nspden)
!Local variables-------------------------------
!scalars
 integer,parameter :: ishift=5
 integer :: ierr,ifft1,ifft2,ifft3,ifft0,ifft,ispden,i1,i2,i3,n1,n2,n3
 integer :: me_fft,my_mpi_comm,nfftot
 logical :: fftgrid_found
 real(dp) :: invn1,invn2,invn3,invnfftot,wrapfft1,wrapfft2,wrapfft3
 character(len=500) :: msg
 type(distribfft_type),pointer :: my_distribfft
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: dipole_tmp(3,nspden)

! *************************************************************************


!Several initializations
 my_mpi_comm=xmpi_comm_self;if (present(mpi_comm_grid)) my_mpi_comm=mpi_comm_grid
 n1=ngfft(1);n2=ngfft(2);n3=ngfft(3)
 nfftot=n1*n2*n3;invnfftot=one/(dble(nfftot))
 invn1=one/real(n1,kind=dp);invn2=one/real(n2,kind=dp);invn3=one/real(n3,kind=dp)
 me_fft=xmpi_comm_rank(my_mpi_comm)

 dipole(:,:)=zero

!Get the distrib associated with the FFT grid
 if (present(distribfft)) then
   my_distribfft => distribfft
 else
   ABI_MALLOC(my_distribfft,)
   call init_distribfft_seq(my_distribfft,'f',n2,n3,'fourdp')
 end if
 fftgrid_found=.false.
 if (n2 == my_distribfft%n2_coarse ) then
   if (n3 == size(my_distribfft%tab_fftdp3_distrib)) then
     fftn3_distrib => my_distribfft%tab_fftdp3_distrib
     ffti3_local => my_distribfft%tab_fftdp3_local
     fftgrid_found=.true.
   end if
 end if
 if (n2 == my_distribfft%n2_fine ) then
   if (n3 == size(my_distribfft%tab_fftdp3dg_distrib)) then
     fftn3_distrib => my_distribfft%tab_fftdp3dg_distrib
     ffti3_local => my_distribfft%tab_fftdp3dg_local
     fftgrid_found = .true.
   end if
 end if
 if (.not.(fftgrid_found)) then
   msg='Unable to find an allocated distrib for the FFT grid!'
   ABI_BUG(msg)
 end if

!Loop over FFT grid points
!$OMP PARALLEL PRIVATE(ifft,i1,i2,ifft1,ifft2,ispden,wrapfft1,wrapfft2)
 do ifft3=1,n3
   i3=mod(ifft3-1+ishift*n3,n3)
   if(fftn3_distrib(1+i3)==me_fft) then
     wrapfft3=mod(real(ifft3-1,kind=dp)*invn3-origin(3)+1.5_dp,one)-half
     ifft0=1+n1*n2*(ffti3_local(1+i3)-1)
!$OMP SINGLE
     dipole_tmp=zero
!$OMP END SINGLE
!$OMP DO COLLAPSE(2) REDUCTION(+:dipole_tmp)
     do ifft2=1,n2
       do ifft1=1,n1
         i2=mod(ifft2-1+ishift*n2,n2)
         i1=mod(ifft1-1+ishift*n1,n1)
         wrapfft2=mod(real(ifft2-1,kind=dp)*invn2-origin(2)+1.5_dp,one)-half
         wrapfft1=mod(real(ifft1-1,kind=dp)*invn1-origin(1)+1.5_dp,one)-half
         ifft=ifft0+n1*i2+i1

!        Computation of integral(s)
         do ispden=1,nspden
           dipole_tmp(1,ispden)=dipole_tmp(1,ispden)+wrapfft1*arraysp(ifft,ispden)
           dipole_tmp(2,ispden)=dipole_tmp(2,ispden)+wrapfft2*arraysp(ifft,ispden)
           dipole_tmp(3,ispden)=dipole_tmp(3,ispden)+wrapfft3*arraysp(ifft,ispden)
         end do

       end do
     end do
!$OMP END DO
!$OMP SINGLE
     dipole=dipole+dipole_tmp
!$OMP END SINGLE
   end if
 end do
!$OMP END PARALLEL

!MPI parallelization
 if (xmpi_comm_size(my_mpi_comm)>1) then
   call xmpi_sum(dipole,my_mpi_comm,ierr)
 end if

!From reduced to cartesian coordinates
 do ispden=1,nspden
   dipole(:,ispden)=matmul(rprimd,dipole(:,ispden))*invnfftot
 end do

 if (.not.present(distribfft)) then
   call destroy_distribfft(my_distribfft)
   ABI_FREE(my_distribfft)
 end if

end subroutine multipoles_fftr
!!***

!!****f* m_multipoles/multipoles_out
!! NAME
!! multipoles_out
!!
!! FUNCTION
!!  Output multipole moments of input array on FFT grid, calculated with multipoles_fftr
!!  Namely, the electrical dipole, quadrupole, etc... of the electron density
!!
!! INPUTS
!!  rhor(nfft,nspden)=electronic density
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms
!!  nfft=number of FFT points stored by one proc
!!  ngfft =number of subdivisions along each lattice vector
!!  nspden=number of spin-density components
!!  ntypat=number of atom types
!!  rprimd(3,3)=dimensionful lattice vectors
!!  typat(ntypat)=types of atoms
!!  ucvol=unit cell volume
!!  unit_out=file unit to print out
!!  ziontypat(ntypat)=ionic charge of each type of atom
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      m_outscfcv
!!
!! CHILDREN
!!      atomdata_from_znucl,ptabs_fourdp,wrtout,xmpi_sum,xred2xcart
!!
!! SOURCE

subroutine multipoles_out(rhor,mpi_enreg,natom,nfft,ngfft,nspden,&
                          ntypat,rprimd,typat,ucvol,unit_out,xred,ziontypat)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,nspden,ntypat,unit_out
 real(dp), intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer, intent(in) :: ngfft(18),typat(natom)
 real(dp),intent(in) :: rhor(nfft,nspden),rprimd(3,3),xred(3,natom),ziontypat(ntypat)
!Local variables ------------------------------
!scalars
 integer :: iatom,nspden_updn
 real(dp) :: ziontotal
 character(len=500) :: message
!arrays
 real(dp) :: center_of_charge(3),dipole_el(3,2),dipole_ions_cart(3)
 real(dp) :: dipole_ions_red(3),dipole_tot(3),tmp(3)

! *************************************************************************

!Separate spins only for nspden=2
 nspden_updn=merge(1,2,nspden/=2)

!Title

!Get nuclear part of dipole
 dipole_ions_red(:) = zero ; ziontotal = zero
 do iatom = 1, natom
   dipole_ions_red(:) = dipole_ions_red(:) + xred(:,iatom)*ziontypat(typat(iatom))
   ziontotal = ziontotal + ziontypat(typat(iatom))
 end do
 dipole_ions_cart(:) = matmul(rprimd,dipole_ions_red(:))

!Find coordinates of center of charge on FFT grid
 center_of_charge(1:3) = dipole_ions_red(1:3)/ziontotal

!Get electronic part of dipole with respect to center of charge of ions (in cart. coord.)
 dipole_el = zero
 call multipoles_fftr(rhor(:,1:nspden_updn),dipole_el(:,1:nspden_updn),nfft,ngfft,nspden_updn,&
& rprimd,center_of_charge,distribfft=mpi_enreg%distribfft,mpi_comm_grid=mpi_enreg%comm_fft)
 dipole_el(1:3,1:nspden_updn)=dipole_el(1:3,1:nspden_updn)*ucvol

!Take into account storage of rhor (up+dn,up)
 if (nspden==2) then
   tmp(1:3)=dipole_el(1:3,1)
   dipole_el(1:3,1)=dipole_el(1:3,2)
   dipole_el(1:3,2)=tmp(1:3)-dipole_el(1:3,2)
 end if

!Compute total dipole
!NOTE: wrt center of charge, dipole_ions is 0
 dipole_tot(1) = - sum(dipole_el(1,1:nspden_updn))
 dipole_tot(2) = - sum(dipole_el(2,1:nspden_updn))
 dipole_tot(3) = - sum(dipole_el(3,1:nspden_updn))

!Output
 write (message, '(2a)') ch10,' ----- Electric nuclear dipole wrt the center of ionic charge ----- '
 call wrtout(unit_out, message, 'COLL')
 write (message, '(a,3(1x,ES12.5))') &
& ' Center of charge for ionic distribution (red. coord.): ',center_of_charge(1:3)
 call wrtout(unit_out, message, 'COLL')
 write (message, '(3a,3(1x,E16.6),3a,3(1x,E16.6),a)') ' -----',ch10,&
& ' Ionic dipole (cart. coord.)     = ',dipole_ions_cart, ' (a.u.)',ch10, &
& '                                 = ',dipole_ions_cart/dipole_moment_debye,' (D)'
 call wrtout(unit_out, message, 'COLL')
 if (nspden/=2) then
   !This is compatible with nspden=4
   write (message, '(3a,3(1x,E16.6),3a,3(1x,E16.6),a)') ' -----',ch10,&
&   ' Electronic dipole (cart. coord.)= ',dipole_el(:,1),' (a.u.)',ch10,&
&   '                                 = ',dipole_el(:,1)/dipole_moment_debye,' (D)'
 else
   write (message, '(3a,3(1x,E16.6),a,3(2a,3(1x,E16.6),a))') ' -----',ch10,&
&   ' Electronic dipole (cart. coord.)= ',dipole_el(:,1),' up (a.u.)',ch10,&
&   '                                 = ',dipole_el(:,2),' dn (a.u.)',ch10,&
&   '                                 = ',dipole_el(:,1)/dipole_moment_debye,' up (D)',ch10,&
&   '                                 = ',dipole_el(:,2)/dipole_moment_debye,' dn (D)'
 end if
 call wrtout(unit_out, message, 'COLL')
 write (message, '(3a,3(1x,E16.6),3a,3(1x,E16.6),a)') ' -----',ch10,&
& ' Total dipole (cart. coord.)     = ',dipole_tot,' (a.u.)',ch10,&
& '                                 = ',dipole_tot/dipole_moment_debye,' (D)'
 call wrtout(unit_out, message, 'COLL')
 call wrtout(unit_out, ' ', 'COLL')

end subroutine multipoles_out
!!***

!!****f* m_multipoles/out1dm
!! NAME
!! out1dm
!!
!! FUNCTION
!! Output the 1 dimensional mean of potential and density on the three reduced axis.
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
!!      m_outscfcv
!!
!! CHILDREN
!!      atomdata_from_znucl,ptabs_fourdp,wrtout,xmpi_sum,xred2xcart
!!
!! SOURCE

subroutine out1dm(fnameabo_app_1dm,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,&
                  rhor,rprimd,typat,ucvol,vtrial,xred,znucl)

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
 write(message, '(a,a)' ) ' io1dm : about to open file ',trim(fnameabo_app_1dm)
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 comm_fft = mpi_enreg%comm_fft; nproc_fft = xmpi_comm_size(comm_fft); me_fft = mpi_enreg%me_fft

 if (me_fft == 0) then
   if (open_file(fnameabo_app_1dm,message,newunit=temp_unit,status='unknown',form='formatted') /= 0) then
     ABI_ERROR(message)
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
 ABI_MALLOC(iden,(natom))
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
     ABI_ERROR(message)
   end if
 end do

!Compute cartesian coordinates, and print reduced and cartesian coordinates
 ABI_MALLOC(xcart,(3,natom))
 call xred2xcart(natom,rprimd,xcart,xred)
 do ia=1,natom
   write(message, '(a,a,3f10.5,a,3f10.5)' ) &
&   '#   ',iden(ia),xred(1:3,ia),'    ',xcart(1:3,ia)
   if (me_fft == 0) write(temp_unit,'(a)') message
 end do
 ABI_FREE(iden)
 ABI_FREE(xcart)

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

   ABI_MALLOC( reduced_coord,(nslice))
   ABI_MALLOC(mean_pot,(nslice))
   ABI_MALLOC(lin_den,(nslice))

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

   ABI_FREE(reduced_coord)
   ABI_FREE(mean_pot)
   ABI_FREE(lin_den)

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

end module m_multipoles
!!***
