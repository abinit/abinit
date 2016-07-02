!{\src2tex{textfont=tt}}
!!****f* ABINIT/linear_optics_paw
!! NAME
!! linear_optics_paw
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! linear susceptiblity using matrix elements <-i Nabla> obtained from a
!! PAW ground state calculation. It uses formula 17 from Gadoc et al,
!! Phys. Rev. B 73, 045112 (2006) together with a scissors correction. It uses
!! a Kramers-Kronig transform to compute the real part from the imaginary part, and
!! it will work on all types of unit cells. It outputs all tensor elements of
!! both the real and imaginary parts.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (VR, PGhosh)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  filnam: base of file names to read data from
!!  mpi_enreg: mpi set up variable, not used in this code
!!
!! OUTPUT
!!  _real and _imag output files
!!
!! NOTES
!!  This routine is not tested
!!
!! PARENTS
!!      conducti
!!
!! CHILDREN
!!      destroy_mpi_enreg,hdr_free,hdr_io,hdr_read_from_fname,initmpi_seq
!!      kramerskronig,matrginv,metric,wffopen
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine linear_optics_paw(filnam,filnam_out,mpi_enreg_seq)


 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_wffile
 use m_profiling_abi
 use m_hdr

 use m_io_tools,   only : open_file, get_unit
 use m_mpinfo,     only : destroy_mpi_enreg, nullify_mpi_enreg

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'linear_optics_paw'
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 character(len=fnlen),intent(in) :: filnam,filnam_out
 type(MPI_type),intent(inout) :: mpi_enreg_seq

!Local variables-------------------------------
 integer,parameter :: master=0
 integer :: iomode,bantot,bdtot_index,fform1,headform
 integer :: iband,ierr,ii,ikpt,iom,iout,isppol,isym,jband,jj,me,mband
 integer :: method,mom,nband_k,nkpt,nspinor,nsppol,nsym,occopt,only_check
 integer :: rdwr,spaceComm,inpunt,reunt,imunt,wfunt
 integer,allocatable :: nband(:),symrel(:,:,:)
 real(dp) :: del,dom,fij,gdelta,omin,omax,paijpbij(2),mbpt_sciss,wij,ucvol
 real(dp) :: diffwp, diffwm
 real(dp) :: e2rot(3,3),gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),rprimdinv(3,3),symd(3,3),symdinv(3,3)
 real(dp),allocatable :: e1(:,:,:),e2(:,:,:,:),epsilon_tot(:,:,:,:),eigen0(:),eig0_k(:)
 real(dp),allocatable :: kpts(:,:),occ(:),occ_k(:),oml1(:),wtk(:)
 complex,allocatable :: eps_work(:)
 character(len=fnlen) :: filnam1,filnam_gen
 character(len=500) :: msg
 type(hdr_type) :: hdr
 type(wffile_type) :: wff1
!arrays
 real(dp),allocatable :: psinablapsi(:,:,:,:)

! *********************************************************************************

 DBG_ENTER("COLL")

!* Fake MPI_type for the sequential part.
!This routine should not be parallelized as communicating gbig and other
!tables takes more time than recalculating them in sequential.
 call initmpi_seq(MPI_enreg_seq)
 
!write(std_out,'(a)')' Give the name of the output file ...'
!read(std_in, '(a)') filnam_out
!write(std_out,'(a)')' The name of the output file is :',filnam_out

!Read data file
 if (open_file(filnam,msg,newunit=inpunt,form='formatted') /= 0 ) then
   MSG_ERROR(msg)
 end if

 rewind(inpunt)
 read(inpunt,*)
 read(inpunt,'(a)')filnam_gen       ! generic name for the files
 filnam1=trim(filnam_gen)//'_OPT' ! nabla matrix elements file

!Open the Wavefunction and optic files
!These default values are typical of sequential use
 iomode=IO_MODE_FORTRAN ; spaceComm=xmpi_comm_self; me=0

! Read the header of the optic files
 call hdr_read_from_fname(hdr, filnam1, fform1, spaceComm)
 call hdr_free(hdr)
 if (fform1 /= 610) then
   MSG_ERROR("Abinit8 requires an OPT file with fform = 610")
 end if

!Open the conducti optic files
 wfunt = get_unit()
 call WffOpen(iomode,spaceComm,filnam1,ierr,wff1,master,me,wfunt)

!Read the header from Ground state file
 rdwr=1
 call hdr_io(fform1,hdr,rdwr,wff1)

!Extract info from the header
 headform=hdr%headform
 bantot=hdr%bantot
 nkpt=hdr%nkpt
 ABI_ALLOCATE(kpts,(3,nkpt))
 ABI_ALLOCATE(wtk,(nkpt))
 kpts(:,:)=hdr%kptns(:,:)
 wtk(:)=hdr%wtk(:)
 nspinor=hdr%nspinor
 nsppol=hdr%nsppol
 occopt=hdr%occopt
 rprimd(:,:)=hdr%rprimd(:,:)
 rprimdinv(:,:) = rprimd(:,:)
 call matrginv(rprimdinv,3,3) ! need the inverse of rprimd to symmetrize the tensors
 ABI_ALLOCATE(nband,(nkpt*nsppol))
 ABI_ALLOCATE(occ,(bantot))
 occ(1:bantot)=hdr%occ(1:bantot)
 nband(1:nkpt*nsppol)=hdr%nband(1:nkpt*nsppol)
 nsym=hdr%nsym
 ABI_ALLOCATE(symrel,(3,3,nsym))
 symrel(:,:,:)=hdr%symrel(:,:,:)

!Get mband, as the maximum value of nband(nkpt)
 mband=maxval(nband(:))

!get ucvol etc.
 iout = -1
 call metric(gmet,gprimd,iout,rmet,rprimd,ucvol)

 write(std_out,*)
 write(std_out,'(a,3f10.5,a)' )' rprimd(bohr)      =',rprimd(1:3,1)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,2)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,3)
 write(std_out,*)
 write(std_out,'(a,3f10.5,a)' )' rprimdinv         =',rprimdinv(1:3,1)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimdinv(1:3,2)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimdinv(1:3,3)
 write(std_out,'(a,2i8)')      ' nkpt,mband        =',nkpt,mband

!get eigen0
 ABI_ALLOCATE(eigen0,(mband*nkpt*nsppol))
 read(wfunt)(eigen0(iband),iband=1,mband*nkpt*nsppol)

 read(inpunt,*)mbpt_sciss
 read(inpunt,*)dom,omin,omax,mom
 close(inpunt)

 ABI_ALLOCATE(oml1,(mom))
 ABI_ALLOCATE(e1,(3,3,mom))
 ABI_ALLOCATE(e2,(2,3,3,mom))
 ABI_ALLOCATE(epsilon_tot,(2,3,3,mom))
 ABI_ALLOCATE(eps_work,(mom))
 del=(omax-omin)/(mom-1)
 do iom=1,mom
   oml1(iom)=omin+dble(iom-1)*del
 end do
 write(std_out,'(a,i8,4f10.5,a)')' npts,omin,omax,width,mbpt_sciss      =',mom,omin,omax,dom,mbpt_sciss,' Ha'

 ABI_ALLOCATE(psinablapsi,(2,3,mband,mband))

!loop over spin components
 do isppol=1,nsppol
   bdtot_index = 0
!  loop over k points
   do ikpt=1,nkpt
!    
!    number of bands for this k point
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     ABI_ALLOCATE(eig0_k,(nband_k))
     ABI_ALLOCATE(occ_k,(nband_k))
!    eigenvalues for this k-point
     eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
!    occupation numbers for this k-point
     occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
!    values of -i*nabla matrix elements for this k point
     psinablapsi=zero
     read(wfunt)((psinablapsi(1:2,1,iband,jband),iband=1,nband_k),jband=1,nband_k)
     read(wfunt)((psinablapsi(1:2,2,iband,jband),iband=1,nband_k),jband=1,nband_k)
     read(wfunt)((psinablapsi(1:2,3,iband,jband),iband=1,nband_k),jband=1,nband_k) 

!    occupation numbers for k-point
     occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
!    accumulate e2 for this k point, Eq. 17 from PRB 73, 045112 (2006)
     do iband = 1, nband_k
       do jband = 1, nband_k
         fij = occ_k(iband) - occ_k(jband) !occ number difference
         wij = eig0_k(iband) - eig0_k(jband) !energy difference
         if (abs(fij) > zero) then ! only consider states of differing occupation numbers
           do ii = 1, 3
             do jj = 1, 3
               paijpbij(1) = psinablapsi(1,ii,iband,jband)*psinablapsi(1,jj,iband,jband) + &
&               psinablapsi(2,ii,iband,jband)*psinablapsi(2,jj,iband,jband)
               paijpbij(2) = psinablapsi(2,ii,iband,jband)*psinablapsi(1,jj,iband,jband) - &
&               psinablapsi(1,ii,iband,jband)*psinablapsi(2,jj,iband,jband)
               do iom = 1, mom
!                original version
!                diffw = wij + mbpt_sciss - oml1(iom) ! apply scissors term here
!                gdelta = exp(-diffw*diffw/(4.0*dom*dom))/(2.0*dom*sqrt(pi)) ! delta fnc resolved as Gaussian
!                e2(1,ii,jj,iom) = e2(1,ii,jj,iom) - (4.0*pi*pi/ucvol)*wtk(ikpt)*fij*paijpbij(1)*gdelta/(oml1(iom)*oml1(iom))
!                e2(2,ii,jj,iom) = e2(2,ii,jj,iom) - (4.0*pi*pi/ucvol)*wtk(ikpt)*fij*paijpbij(2)*gdelta/(oml1(iom)*oml1(iom))
                 diffwm = wij - mbpt_sciss + oml1(iom) ! apply scissors term here
                 diffwp = wij + mbpt_sciss - oml1(iom) ! apply scissors term here
                 gdelta = exp(-diffwp*diffwp/(4.0*dom*dom))/(2.0*dom*sqrt(pi))
                 e2(1,ii,jj,iom) = e2(1,ii,jj,iom) - (4.0*pi*pi/ucvol)*wtk(ikpt)*fij*paijpbij(1)*gdelta/(wij*wij)
                 e2(2,ii,jj,iom) = e2(2,ii,jj,iom) - (4.0*pi*pi/ucvol)*wtk(ikpt)*fij*paijpbij(2)*gdelta/(wij*wij)
               end do ! end loop over spectral points
             end do ! end loop over jj = 1, 3
           end do ! end loop over ii = 1, 3
         end if ! end selection on fij /= 0
       end do ! end loop over jband
     end do ! end loop over iband

     ABI_DEALLOCATE(eig0_k)
     ABI_DEALLOCATE(occ_k)
     bdtot_index=bdtot_index+nband_k
   end do ! end loop over k points
 end do ! end loop over spin polarizations

!here apply nsym symrel transformations to reconstruct full tensor from IBZ part
 epsilon_tot(:,:,:,:) = zero
 do isym = 1, nsym
   symd(:,:)=matmul(rprimd(:,:),matmul(symrel(:,:,isym),rprimdinv(:,:)))
   symdinv(:,:)=symd(:,:)
   call matrginv(symdinv,3,3)
   do iom = 1, mom
     e2rot(:,:)=matmul(symdinv(:,:),matmul(e2(1,:,:,iom),symd(:,:)))
     epsilon_tot(2,:,:,iom) = epsilon_tot(2,:,:,iom)+e2rot(:,:)/nsym
   end do
 end do

!generate e1 from e2 via KK transforma
 method=0 ! use naive integration ( = 1 for simpson)
 only_check=0 ! compute real part of eps in kk routine
 do ii = 1, 3
   do jj = 1, 3
     eps_work(:) = cmplx(0.0,epsilon_tot(2,ii,jj,:))
     call kramerskronig(mom,oml1,eps_work,method,only_check)
     epsilon_tot(1,ii,jj,:) = real(eps_work(:))
     if (ii /= jj) epsilon_tot(1,ii,jj,:) = epsilon_tot(1,ii,jj,:)- 1.0
   end do ! end loop over jj
 end do ! end loop over ii

 if (open_file(trim(filnam_out)//'_imag',msg,newunit=reunt,form='formatted') /= 0) then
   MSG_ERROR(msg)
 end if 

 if (open_file(trim(filnam_out)//'_real',msg,unit=imunt,form='formatted') /= 0) then
   MSG_ERROR(msg)
 end if

 write(reunt,'(a12,6a13)')' # Energy/Ha ','eps_2_xx','eps_2_yy','eps_2_zz',&
& 'eps_2_yz','eps_2_xz','eps_2_xy'
 write(imunt,'(a12,6a13)')' # Energy/Ha ','eps_1_xx','eps_1_yy','eps_1_zz',&
& 'eps_1_yz','eps_1_xz','eps_1_xy'

 do iom = 1, mom
   write(reunt,'(ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4)') oml1(iom),' ',&
&   epsilon_tot(2,1,1,iom),' ',epsilon_tot(2,2,2,iom),' ',epsilon_tot(2,3,3,iom),' ',&
&   epsilon_tot(2,2,3,iom),' ',epsilon_tot(2,1,3,iom),' ',epsilon_tot(2,1,2,iom)
   write(imunt,'(ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4,a,ES12.4)') oml1(iom),' ',&
&   epsilon_tot(1,1,1,iom),' ',epsilon_tot(1,2,2,iom),' ',epsilon_tot(1,3,3,iom),' ',&
&   epsilon_tot(1,2,3,iom),' ',epsilon_tot(1,1,3,iom),' ',epsilon_tot(1,1,2,iom)
 end do

 close(reunt)
 close(imunt)

 ABI_DEALLOCATE(nband)
 ABI_DEALLOCATE(oml1)
 ABI_DEALLOCATE(e2)
 ABI_DEALLOCATE(e1)
 ABI_DEALLOCATE(occ)
 ABI_DEALLOCATE(psinablapsi)
 ABI_DEALLOCATE(eigen0)
 ABI_DEALLOCATE(wtk)
 ABI_DEALLOCATE(kpts)

 call hdr_free(hdr)
 call destroy_mpi_enreg(MPI_enreg_seq)

 DBG_EXIT("COLL")

 end subroutine linear_optics_paw
!!***
