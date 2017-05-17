!{\src2tex{textfont=tt}}
!!****f* ABINIT/emispec_paw
!! NAME
!! emispec_paw
!!
!! FUNCTION
!! This program computes the elements of the emission spectra 
!! from the Kubo-Greenwood formula for PAW formalism
!! largely inspired from the conducti_core_paw routine
!!
!! COPYRIGHT
!! Copyright (C) 2002-2017 ABINIT group (SM, SVinko)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! NOTES
!!  bantot
!!  dom=frequency range
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree).
!!  ecut=kinetic energy planewave cutoff (hartree).
!!  fermie= fermi energy (Hartree)
!!  mom=number of frequency for conductivity computation
!!  mband=maximum number of bands.
!!  natom = number of atoms in the unit cell.
!!  nband(nkpt*nsppol)=number of bands at each RF k point for each spin.
!!  nkpt=number of k points in the IBZ for this perturbation
!!  ngfft(3)=integer fft box dimensions.
!!  nspinor=number of spinorial components of the wavefunctions.
!!  nsppol=1 for unpolarized, 2 for spin-polarized.
!!  ntypat = number of atom types.
!!  occ(mband*nkpt*nsppol)=occupation number for each band and k.
!!  occopt==option for occupancies
!!  psinablapsi2(2,3,mband,nphicor,natom)Matrix elements = <Phi_core|Nabla|Phi_i>
!!  rmet(3,3)=real space metric ($\textrm{bohr}^{2}$).sigx(mom,nphicor))
!!  rprimd(3,3)=real space primitive translations.
!!  of primitive translations.
!!  ucvol=unit cell volume in ($\textrm{bohr}^{3}$).
!!  wind=frequency windows for computations of sigma
!!  wtk(nkpt)=weight assigned to each k point.
!!
!! PARENTS
!!      conducti
!!
!! CHILDREN
!!      hdr_free,hdr_io,hdr_read_from_fname,metric,wffclose,wffopen
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine emispec_paw(filnam,filnam_out,mpi_enreg)

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_wffile
 use m_profiling_abi
 use m_errors
 use m_hdr

 use m_io_tools,     only : open_file, get_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'emispec_paw'
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 character(len=fnlen) :: filnam,filnam_out
 type(MPI_type),intent(in) :: mpi_enreg

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: iomode,atnbr,bantot,bdtot_index
 integer :: fform2,headform,iatom,iband,icor,ierr,ikpt
 integer :: iom,isppol,l1,mband,me,mom
 integer :: natom,nband_k,nkpt,nphicor,nspinor,nsppol,ntypat
 integer :: occopt,rdwr,spaceComm,iunt,ems_unt,opt2_unt
 real(dp) :: del,ecut,fermie
 real(dp) :: omin,omax,dom,oml
 real(dp) :: Tatm,tsmear,ucvol
 character(len=fnlen) :: filnam2,filnam_gen
 character(len=500) :: msg
 type(hdr_type) :: hdr
 type(wffile_type) :: wff2
!arrays
 integer,allocatable :: nband(:),ncor(:),lcor(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: dhdk2_g(:,:,:)
 real(dp),allocatable :: eig0_k(:),eigen0(:)
 real(dp),allocatable :: energy_cor(:)
 real(dp),allocatable :: occ(:),occ_k(:),oml1(:)
 real(dp),allocatable :: psinablapsi2(:,:,:,:,:)
 real(dp),allocatable :: sigx(:,:,:),sigx_av(:,:),wtk(:)

! *********************************************************************************
 ABI_UNUSED(mpi_enreg%paral_kgb)

!Read output file name
!write(std_out,'(a)')' emispec_paw : Please, give the name of the output file ...'
!read(std_in, '(a)') filnam_out
!write(std_out,'(a)')' The name of the output file is :',trim(filnam_out)
!Read data file
 if (open_file(filnam,msg,newunit=iunt,form='formatted',action="read") /=0) then
   MSG_ERROR(msg)
 end if
 rewind(iunt)
 read(iunt,*)
 read(iunt,'(a)')filnam_gen       ! generic name for the files
 filnam2=trim(filnam_gen)//'_OPT2'
!Read size of the frequency range
 read(iunt,*) dom,omin,omax,mom,atnbr
 close(iunt)

 write(std_out,'(a,i8,3f10.5,a)')' npts,omin,omax,width      =',mom,omin,omax,dom,' Ha'
 write(std_out,'(a)')'--------------------------------------------'
 write(std_out,'(a,i4)') 'selected atom for X ray emission',atnbr
 write(std_out,'(a)')'--------------------------------------------'

!Open the Wavefunction and optic files
!These default values are typical of sequential use
 iomode=IO_MODE_FORTRAN; spaceComm=xmpi_comm_self; me=0

! Read the header of the OPT2 file.
 call hdr_read_from_fname(hdr, filnam2, fform2, spaceComm)
 call hdr_free(hdr)

 if (fform2 /= 611) then
   MSG_ERROR("Abinit8 requires an OPT2 file with fform = 611")
 end if

!Open the optic files
 opt2_unt = get_unit()
 call WffOpen(iomode,spaceComm,filnam2,ierr,wff2,master,me,opt2_unt)

!Read the header 
 rdwr=1
 call hdr_io(fform2,hdr,rdwr,wff2)

!Extract info from the header
 headform=hdr%headform
 bantot=hdr%bantot
 ecut=hdr%ecut_eff
 natom=hdr%natom
 nkpt=hdr%nkpt
 nspinor=hdr%nspinor
 nsppol=hdr%nsppol
 ntypat=hdr%ntypat
 occopt=hdr%occopt
 rprimd(:,:)=hdr%rprimd(:,:)
 ABI_ALLOCATE(nband,(nkpt*nsppol))
 ABI_ALLOCATE(occ,(bantot))
 ABI_ALLOCATE(wtk,(nkpt))
 fermie=hdr%fermie
 tsmear=hdr%tsmear
 occ(1:bantot)=hdr%occ(1:bantot)
 wtk(1:nkpt)=hdr%wtk(1:nkpt)
 nband(1:nkpt*nsppol)=hdr%nband(1:nkpt*nsppol)

!Get mband, as the maximum value of nband(nkpt)
 mband=maxval(nband(:))

 write(std_out,*)
 write(std_out,'(a,3f10.5,a)' )' rprimd(bohr)      =',rprimd(1,1:3)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(2,1:3)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(3,1:3)
 write(std_out,'(a,i8)')       ' natom             =',natom
 write(std_out,'(a,3i8)')      ' nkpt,mband,nsppol        =',nkpt,mband,nsppol
 write(std_out, '(a, f10.5,a)' ) ' ecut              =',ecut,' Ha'
 write(std_out,'(a,f10.5,a,f10.5,a)' )' fermie            =',fermie,' Ha',fermie*Ha_eV,' eV'
 Tatm=tsmear*Ha_K
 write(std_out,'(a,f12.5,a,f12.5,a)') ' Temp              =',tsmear,' Ha ',Tatm,' Kelvin'

 ABI_ALLOCATE(eigen0,(mband*nkpt*nsppol))
 read(opt2_unt)(eigen0(iband),iband=1,mband*nkpt*nsppol)

 write(std_out,'(a)')'--------------------------------------------'
 read(opt2_unt) nphicor
 write(std_out,'(a,i4)') 'Number of core orbitals nc=',nphicor
 ABI_ALLOCATE(ncor,(nphicor))
 ABI_ALLOCATE(lcor,(nphicor))
 ABI_ALLOCATE(energy_cor,(nphicor))
 do icor=1,nphicor
   read(opt2_unt) ncor(icor),lcor(icor),energy_cor(icor)
   write(std_out,'(a,2i4,f10.5)') ' n, l, Energy (Ha): ',ncor(icor),lcor(icor),energy_cor(icor)
 end do
 write(std_out,'(a)')'--------------------------------------------'
 
!---------------------------------------------------------------------------------
!gmet inversion to get ucvol
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!---------------------------------------------------------------------------------
!frequency range fixed by the occupation of the valence states
 bdtot_index=0
 omax=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     do iband=1,nband_k
       if(eigen0(1+bdtot_index)>=omax) omax=eigen0(1+bdtot_index)
       bdtot_index=bdtot_index+1
     end do
   end do
 end do
 omin=eigen0(1)
 del=(omax-omin)/(mom-1)
 ABI_ALLOCATE(oml1,(mom))
 do iom=1,mom
   oml1(iom)=omin+dble(iom-1)*del
 end do
 write(std_out,*) 'Valence state orbital energies: omin,omax',omin,omax
 ABI_ALLOCATE(sigx,(natom,mom,nphicor))
 ABI_ALLOCATE(sigx_av,(mom,nphicor))
!---------------------------------------------------------------------------------
!emission X  -------
!
 ABI_ALLOCATE(psinablapsi2,(2,3,mband,nphicor,natom))
 sigx=zero
 sigx_av=zero
 bdtot_index = 0

!LOOP OVER SPINS
 do isppol=1,nsppol

!  LOOP OVER K-POINTS
   do ikpt=1,nkpt
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     ABI_ALLOCATE(eig0_k,(nband_k))
     ABI_ALLOCATE(occ_k,(nband_k))
     ABI_ALLOCATE(dhdk2_g,(natom,nband_k,nphicor))

     dhdk2_g   = zero
     psinablapsi2=zero

!    eigenvalue for k-point
     eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)

!    occupation numbers for k-point
     occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)

!    dipole matrix elements
     do iatom=1,natom
       read(opt2_unt) ((psinablapsi2(1:2,1,iband,icor,iatom),iband=1,nband_k),icor=1,nphicor)
       read(opt2_unt) ((psinablapsi2(1:2,2,iband,icor,iatom),iband=1,nband_k),icor=1,nphicor)
       read(opt2_unt) ((psinablapsi2(1:2,3,iband,icor,iatom),iband=1,nband_k),icor=1,nphicor)
     end do

!    loops over atom, bands, core states 
     do iatom=1,natom
       do iband=1,nband_k
         do icor=1,nphicor

           do l1=1,3
             dhdk2_g(iatom,iband,icor)=dhdk2_g(iatom,iband,icor)+( &
&             psinablapsi2(1,l1,iband,icor,iatom)*psinablapsi2(1,l1,iband,icor,iatom) &
&             +psinablapsi2(2,l1,iband,icor,iatom)*psinablapsi2(2,l1,iband,icor,iatom))
           end do
           
!          emission for each omega
           do iom=1,mom
             oml=-energy_cor(icor)-oml1(iom) 
             sigx(iatom,iom,icor)=sigx(iatom,iom,icor)+ wtk(ikpt)*dhdk2_g(iatom,iband,icor)&
&             *occ_k(iband)/oml*dexp(-((-energy_cor(icor)+eig0_k(iband)-oml)/dom)**2)        
           end do
         end do !icor
       end do  !iband
     end do !iatom
     bdtot_index=bdtot_index+nband_k
     ABI_DEALLOCATE(eig0_k)
     ABI_DEALLOCATE(occ_k)
     ABI_DEALLOCATE(dhdk2_g)
!    end loop over k
   end do
!  end loop over spins
 end do
 ABI_DEALLOCATE(psinablapsi2)

 do iatom=1,natom
   do icor=1,nphicor
     do iom=1,mom
       if(sigx(iatom,iom,icor)<=tol16) sigx(iatom,iom,icor)=zero
     end do
   end do 
 end do ! iatom

 sigx=sigx*two_pi*third*dble(natom)/(dom*ucvol)*half/dsqrt(pi)

 do icor=1,nphicor
   do iom=1,mom
     do iatom=1,natom
       sigx_av(iom,icor) =sigx_av(iom,icor)+sigx(iatom,iom,icor)/dble(natom)
     end do
   end do
 end do 

 if (open_file(trim(filnam_out)//'_emisX',msg,newunit=ems_unt,form='formatted', action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 do iom=1,mom
   write(ems_unt,'(9(1x,e15.8))') &
&   ((-energy_cor(icor)+oml1(iom)),sigx_av(iom,icor),sigx(atnbr,iom,icor),icor=1,nphicor)
 end do
 close(ems_unt)

 call WffClose(wff2,ierr)

 ABI_DEALLOCATE(sigx)
 ABI_DEALLOCATE(sigx_av)
 ABI_DEALLOCATE(ncor)
 ABI_DEALLOCATE(lcor)
 ABI_DEALLOCATE(energy_cor)
 ABI_DEALLOCATE(eigen0)
 ABI_DEALLOCATE(nband)
 ABI_DEALLOCATE(oml1)
 ABI_DEALLOCATE(occ)
 ABI_DEALLOCATE(wtk)

 call hdr_free(hdr)

end subroutine emispec_paw
!!***
