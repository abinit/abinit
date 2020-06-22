!!****m* ABINIT/m_conducti
!! NAME
!!  m_conducti
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! conductivity tensor and the conductivity along the three principal axes
!! from the Kubo-Greenwood formula.
!!
!! COPYRIGHT
!!  Copyright (C) 2002-2020 ABINIT group (VRecoules, PGhosh, SMazevet, SM, SVinko)
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

module m_conducti

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_wffile
 use m_wfk
 use m_hdr
 use m_nctk

 use defs_abitypes,  only : MPI_type
 use m_io_tools,     only : open_file, get_unit
 use m_fstrings,     only : sjoin
 use m_symtk,        only : matr3inv
 use m_hide_lapack,  only : jacobi
 use m_occ,          only : getnel
 use m_geometry,     only : metric
 use m_splines,      only : intrpl

 implicit none

 private
!!***

 public :: conducti_paw
 public :: conducti_paw_core
 public :: conducti_nc
 public :: emispec_paw
!!***

contains
!!***

!!****f* m_conducti/conducti_paw
!! NAME
!! conducti_paw
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! conductivity tensor and the conductivity along the three principal axes
!! from the Kubo-Greenwood formula for PAW formalism
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! NOTES
!!  bantot
!!  doccde(mband*nkpt_rbz*nsppol)=derivative of occ_rbz wrt the energy.
!!  dom=frequency range
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree).
!!  eigen11(2,nkpt,mband,mband,nsppol)=first-order eigenvalues (hartree)
!!  in direction x
!!  eigen12(2,nkpt,mband,mband,nsppol)=first-order eigenvalues (hartree)
!!  in direction y
!!  eigen13(2,nkpt,mband,mband,nsppol)=first-order eigenvalues (hartree)
!!  in direction z
!!  ecut=kinetic energy planewave cutoff (hartree).
!!  fermie= fermi energy (Hartree)
!!  gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{2}$).
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1).
!!  kin11= Onsager kinetic coeficient=optical conductivity
!!  kin12= Onsager kinetic coeficient
!!  kin21= Onsager kinetic coeficient
!!  kin22= Onsager kinetic coeficient
!!  Kth=thermal conductivity
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
!!  rmet(3,3)=real space metric ($\textrm{bohr}^{2}$).sigx(mom,nphicor))
!!  rprimd(3,3)=real space primitive translations.
!!  of primitive translations.
!!  Sth=thermopower
!!  tsmear=smearing width (or temperature) in Hartree
!!  ucvol=unit cell volume in ($\textrm{bohr}^{3}$).
!!  wind=frequency windows for computations of sigma
!!  wtk(nkpt)=weight assigned to each k point.
!!  znucl(natom)=atomic number of atoms
!!  np_sum=noziere-pines sumrule
!!
!! PARENTS
!!      conducti
!!
!! CHILDREN
!!      hdr_free,hdr_io,hdr_read_from_fname,metric,msig,wffclose,wffopen
!!
!! SOURCE

 subroutine conducti_paw(filnam,filnam_out,mpi_enreg)

!Arguments -----------------------------------
!scalars
 character(len=fnlen) :: filnam,filnam_out
 type(MPI_type),intent(in) :: mpi_enreg

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: iomode,bantot,bdtot_index
 integer :: fform1,headform,iband,ierr,ikpt
 integer :: iom,isppol,jband,l1,l2,mband,me,mom
 integer :: natom,nband_k,nkpt,nspinor,nsppol,ntypat
 integer :: occopt,rdwr,spaceComm,iunt,opt_unt
 integer :: lij_unt,sig_unt,kth_unt,ocond_unt
 real(dp) :: del,deltae,diff_occ,ecut,fermie,maxocc
 real(dp) :: np_sum,np_sum_k1,np_sum_k2,omin,omax,dom,oml,sig,socc,socc_k
 real(dp) :: Tatm,tphysel,tsmear,ucvol
 character(len=fnlen) :: filnam1,filnam_gen
 character(len=500) :: msg
 type(hdr_type) :: hdr
 type(wffile_type) :: wff1
!arrays
 integer,allocatable :: nband(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: cond_nd(:,:,:),dhdk2_r(:,:,:,:),dhdk2_g(:,:,:)
 real(dp),allocatable :: doccde(:),doccde_k(:),eig0_k(:),eigen0(:)
 real(dp),allocatable :: occ(:),occ_k(:),wtk(:),oml1(:)
 real(dp),allocatable :: kin11(:,:),kin12(:),kin21(:),kin22(:)
 real(dp),allocatable :: kin11_k(:),kin12_k(:),kin21_k(:),kin22_k(:),Kth(:),Stp(:)
 real(dp),allocatable :: psinablapsi(:,:,:,:),sig_abs(:)

! *********************************************************************************
 ABI_UNUSED(mpi_enreg%paral_kgb)

!write(std_out,'(a)')' The name of the output file is :',trim(filnam_out)
!Read data file
 if (open_file(filnam, msg, newunit=iunt, form='formatted', action="read", status="old") /=0 ) then
   MSG_ERROR(msg)
 end if
 rewind(iunt)
 read(iunt,*)
 read(iunt,'(a)')filnam_gen       ! generic name for the files
 filnam1=trim(filnam_gen)//'_OPT'
!Read size of the frequency range
 read(iunt,*) dom,omin,omax,mom
 close(iunt)
 write(std_out,'(a,i8,3f10.5,a)')' npts,omin,omax,width      =',mom,omin,omax,dom,' Ha'

!These default values are typical of sequential use
 iomode=IO_MODE_FORTRAN ; spaceComm=xmpi_comm_self; me=0

! Read the header of the optic files
 call hdr_read_from_fname(hdr, filnam1, fform1, spaceComm)
 call hdr%free()
 if (fform1 /= 610) then
   MSG_ERROR("Abinit8 requires an OPT file with fform = 610")
 end if

!Open the conducti and/or optic files
 opt_unt = get_unit()
 call WffOpen(iomode,spaceComm,filnam1,ierr,wff1,master,me,opt_unt)

!Read the header from Ground state file
 rdwr=1
 call hdr_io(fform1,hdr,rdwr,wff1)
 ABI_CHECK(fform1/=0,"Error opening wff1")

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
 write(std_out,'(a,3f10.5,a)' )' rprimd(bohr)      =',rprimd(1:3,1)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,2)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,3)
 write(std_out,'(a,i8)')       ' natom             =',natom
 write(std_out,'(a,3i8)')      ' nkpt,mband,nsppol        =',nkpt,mband,nsppol
 write(std_out, '(a, f10.5,a)' ) ' ecut              =',ecut,' Ha'
 write(std_out,'(a,f10.5,a,f10.5,a)' )' fermie            =',fermie,' Ha',fermie*Ha_eV,' eV'
 Tatm=tsmear*Ha_K
 write(std_out,'(a,f12.5,a,f12.5,a)') ' Temp              =',tsmear,' Ha ',Tatm,' Kelvin'

 ABI_ALLOCATE(eigen0,(mband*nkpt*nsppol))
 read(opt_unt)(eigen0(iband),iband=1,mband*nkpt*nsppol)
!
!
!---------------------------------------------------------------------------------
!gmet inversion to get ucvol
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!---------------------------------------------------------------------------------
!derivative of occupation wrt the energy.
 ABI_ALLOCATE(doccde,(mband*nkpt*nsppol))
 if (occopt==1) then
   write(std_out,'(a,i4)')  ' occopt            =',occopt
   doccde=zero
 else
   tphysel=zero
   maxocc=two/(nsppol*nspinor)
 end if
!---------------------------------------------------------------------------------
!size of the frequency range
 del=(omax-omin)/(mom-1)
 ABI_ALLOCATE(oml1,(mom))
 do iom=1,mom
   oml1(iom)=omin+dble(iom-1)*del
 end do
 ABI_ALLOCATE(cond_nd,(mom,3,3))
 ABI_ALLOCATE(kin11,(mom,nsppol))
 ABI_ALLOCATE(kin12,(mom))
 ABI_ALLOCATE(kin21,(mom))
 ABI_ALLOCATE(kin22,(mom))
 ABI_ALLOCATE(sig_abs,(mom))
 ABI_ALLOCATE(kin11_k,(mom))
 ABI_ALLOCATE(kin12_k,(mom))
 ABI_ALLOCATE(kin21_k,(mom))
 ABI_ALLOCATE(kin22_k,(mom))
 ABI_ALLOCATE(Kth,(mom))
 ABI_ALLOCATE(Stp,(mom))

!---------------------------------------------------------------------------------
!Conductivity -------
!
 ABI_ALLOCATE(psinablapsi,(2,3,mband,mband))
 kin11   = zero
 kin12   = zero
 kin21   = zero
 kin22   = zero
 np_sum  = zero
 socc    = zero
 sig_abs = zero

 bdtot_index = 0

!LOOP OVER SPINS/K
 deltae  = zero
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     ABI_ALLOCATE(eig0_k,(nband_k))
     ABI_ALLOCATE(occ_k,(nband_k))
     ABI_ALLOCATE(doccde_k,(nband_k))
     ABI_ALLOCATE(dhdk2_r,(nband_k,nband_k,3,3))
     ABI_ALLOCATE(dhdk2_g,(natom,nband_k,nband_k))

     cond_nd   = zero
     kin11_k   = zero
     kin12_k   = zero
     kin21_k   = zero
     kin22_k   = zero
     np_sum_k1 = zero
     np_sum_k2 = zero
     socc_k    = zero
     dhdk2_r   = zero
     dhdk2_g   = zero

!    eigenvalue for k-point
     eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
!    first derivative eigenvalues for k-point
     psinablapsi=zero
     read(opt_unt)((psinablapsi(1:2,1,iband,jband),iband=1,nband_k),jband=1,nband_k)
     read(opt_unt)((psinablapsi(1:2,2,iband,jband),iband=1,nband_k),jband=1,nband_k)
     read(opt_unt)((psinablapsi(1:2,3,iband,jband),iband=1,nband_k),jband=1,nband_k)
!    DEBUG
!    write(963,*)isppol,ikpt,((psinablapsi(1:2,1,iband,jband),iband=1,nband_k),jband=1,nband_k)
!    write(963,*)isppol,ikpt,((psinablapsi(1:2,2,iband,jband),iband=1,nband_k),jband=1,nband_k)
!    write(963,*)isppol,ikpt,((psinablapsi(1:2,3,iband,jband),iband=1,nband_k),jband=1,nband_k)
!    ENDDEBUG

!    occupation numbers for k-point
     occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
!    derivative of occupation number for k-point
     doccde_k(:)=doccde(1+bdtot_index:nband_k+bdtot_index)

!    LOOP OVER BANDS
     do iband=1,nband_k
       do jband=1,nband_k
         do l1=1,3
           do l2=1,3
             dhdk2_r(iband,jband,l1,l2)=dhdk2_r(iband,jband,l1,l2)+(&
&             psinablapsi(1,l1,iband,jband)*psinablapsi(1,l2,iband,jband)&
&             +psinablapsi(2,l1,iband,jband)*psinablapsi(2,l2,iband,jband))
           end do
         end do

         do l1=1,3
           dhdk2_g(1,iband,jband)=dhdk2_g(1,iband,jband)+( &
&           psinablapsi(1,l1,iband,jband)*psinablapsi(1,l1,iband,jband) &
&           +psinablapsi(2,l1,iband,jband)*psinablapsi(2,l1,iband,jband))
         end do

         diff_occ = occ_k(iband)-occ_k(jband)
         if (dabs(diff_occ)>=tol8) then

!          Conductivity for each omega
!          omin = zero
           do iom=1,mom
             oml=oml1(iom)
             if (jband>iband) then
               sig= dhdk2_g(1,iband,jband)&
&               *(diff_occ)/oml*(dexp(-((eig0_k(jband)-eig0_k(iband)-oml)/dom)**2)&
&               -dexp(-((eig0_k(iband)-eig0_k(jband)-oml)/dom)**2))
               kin11_k(iom)=kin11_k(iom)+sig
               kin12_k(iom)=kin12_k(iom)-sig*(eig0_k(jband)-fermie)
               kin21_k(iom)=kin21_k(iom)-sig*(eig0_k(iband)-fermie)
               kin22_k(iom)=kin22_k(iom) + &
&               sig*(eig0_k(iband)-fermie)*(eig0_k(jband)-fermie)
             end if
             do l1=1,3
               do l2=1,3
                 cond_nd(iom,l1,l2)=cond_nd(iom,l1,l2) +dhdk2_r(iband,jband,l1,l2)&
&                 *(diff_occ)/oml*dexp(-((eig0_k(jband)-eig0_k(iband)-oml)/dom)**2)
               end do
             end do
           end do

!          Sumrule start
           if (dabs(eig0_k(iband)-eig0_k(jband))>=tol10) then
             np_sum_k1=np_sum_k1 -dhdk2_g(1,iband,jband)&
&             *(diff_occ)/(eig0_k(iband)-eig0_k(jband))
           else
             np_sum_k2=np_sum_k2 - doccde_k(iband)*dhdk2_g(1,iband,jband)
           end if

!          end loop over band
         end if
       end do
       socc_k=socc_k+occ_k(iband)
     end do
     do iom=1,mom
       kin11(iom,isppol)=kin11(iom,isppol)+wtk(ikpt)*kin11_k(iom)
       kin12(iom)=kin12(iom)+wtk(ikpt)*kin12_k(iom)
       kin21(iom)=kin21(iom)+wtk(ikpt)*kin21_k(iom)
       kin22(iom)=kin22(iom)+wtk(ikpt)*kin22_k(iom)
     end do
     np_sum=np_sum + wtk(ikpt)*(np_sum_k1+np_sum_k2)
     socc=socc+wtk(ikpt)*socc_k

!    Validity limit
     deltae=deltae+(eig0_k(nband_k)-fermie)

     bdtot_index=bdtot_index+nband_k
     ABI_DEALLOCATE(eig0_k)
     ABI_DEALLOCATE(occ_k)
     ABI_DEALLOCATE(doccde_k)
     ABI_DEALLOCATE(dhdk2_r)
     ABI_DEALLOCATE(dhdk2_g)
!    End loop over k
   end do
!  End loop over Spin
 end do

 write(std_out,'(a,3f10.5)')' sumrule           =',np_sum/socc/three/dble(nsppol),socc
 write(std_out,'(a,f10.5,a,f10.5,a)')&
& ' Emax-Efermi       =',deltae/dble(nkpt*nsppol),' Ha',deltae/dble(nkpt*nsppol)*Ha_eV,' eV'

 if (open_file(trim(filnam_out)//'_Lij',msg, newunit=lij_unt, form='formatted', action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write(lij_unt,'(a)')' # omega(ua) L12 L21 L22 L22'

 if (open_file(trim(filnam_out)//'_sig', msg, newunit=sig_unt, form='formatted', action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 if (nsppol==1) then
   write(sig_unt,'(a)')' # omega(ua) hbar*omega(eV)    cond(ua)             cond(ohm.cm)-1'
 else
   write(sig_unt,'(2a)')' # omega(ua) hbar*omega(eV)      cond(ua)            cond(ohm.cm)-1',&
&   '      cond(ohm.cm)-1 UP      cond(ohm.cm)-1 DN'
 end if

 if (open_file(trim(filnam_out)//'_Kth', msg, newunit=kth_unt, form='formatted', action="write") /=0) then
   MSG_ERROR(msg)
 end if
 write(kth_unt,'(a)')&
& ' #omega(ua) hbar*omega(eV)  thermal cond(ua)   Kth(W/m/K)   thermopower(ua)   Stp(microohm/K)'

 if (open_file(trim(filnam_out)//'.out', msg, newunit=ocond_unt, form='formatted', action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write(ocond_unt,'(a)' )' #Conducti output file:'
 write(ocond_unt,'(a)' )' #Contains all results produced by conducti utility'
 write(ocond_unt,'(a)' )' '
 write(ocond_unt,'(a)')' # omega(ua)       cond(ua)             thermal cond(ua)       thermopower(ua)'

!call isfile(filnam_out,'new')

!Compute thermal conductivity and thermopower
 do iom=1,mom
   oml=oml1(iom)
   do isppol=1,nsppol
     kin11(iom,isppol)=kin11(iom,isppol)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
     if (dabs(kin11(iom,isppol))<10.0d-20) kin11(iom,isppol)=zero
     sig_abs(iom)=sig_abs(iom)+kin11(iom,isppol)
   end do
   kin21(iom)=kin21(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
   kin12(iom)=kin12(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
   kin22(iom)=kin22(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
   Kth(iom)=kin22(iom)
   Stp(iom)=zero
   if(sig_abs(iom)/=zero)  then
     Kth(iom)=Kth(iom)-(kin12(iom)*kin21(iom)/sig_abs(iom))
     Stp(iom)=kin12(iom)/(sig_abs(iom)*Tatm)
   end if
   if (dabs(Kth(iom))<10.0d-20) Kth(iom)=zero
   if (dabs(Stp(iom))<10.0d-20) Stp(iom)=zero
   if (abs(kin12(iom))<10.0d-80) kin12=zero
   if (abs(kin21(iom))<10.0d-80) kin21=zero
   if (abs(kin22(iom))<10.0d-80) kin22=zero

   write(lij_unt,'(f12.5,4es22.12)')oml,kin12(iom),kin21(iom),kin22(iom),kin22(iom)/Tatm*3.4057d9

   if (nsppol==1) then
     write(sig_unt,'(2f12.5,2es22.12)') oml,oml*Ha_eV,sig_abs(iom),sig_abs(iom)*Ohmcm
   else
     write(sig_unt,'(2f12.5,4es22.12)') oml,oml*Ha_eV,sig_abs(iom),sig_abs(iom)*Ohmcm,&
&     kin11(iom,1)*Ohmcm,kin11(iom,2)*Ohmcm
   end if
   write(kth_unt,'(2f12.5,4es22.12)') oml,oml*Ha_eV,Kth(iom),Kth(iom)*3.4057d9/Tatm,&
&   Stp(iom),Stp(iom)*3.6753d-2
   write(ocond_unt,'(1f12.5,3es22.12)') oml,sig_abs(iom),Kth(iom),Stp(iom)
 end do

!Calculate the imaginary part of the conductivity (principal value)
!+derived optical properties.
 call msig (sig_abs,mom,oml1,filnam_out)

 close(lij_unt)
 close(sig_unt)
 close(kth_unt)
 close(ocond_unt)

 write(std_out,'(2a)')ch10,'OUTPUT'
 write(std_out,'(a)')trim(filnam_out)//'_Lij : Onsager kinetic coefficients'
 write(std_out,'(a)')trim(filnam_out)//'_sig : Optical conductivity'
 write(std_out,'(a)')trim(filnam_out)//'_Kth : Thermal conductivity and thermopower'
 write(std_out,'(a)')trim(filnam_out)//'_eps : Dielectric fonction'
 write(std_out,'(a)')trim(filnam_out)//'_abs : n, k, reflectivity, absorption'

 call WffClose(wff1,ierr)

 ABI_DEALLOCATE(psinablapsi)
 ABI_DEALLOCATE(kin11)
 ABI_DEALLOCATE(kin22)
 ABI_DEALLOCATE(kin12)
 ABI_DEALLOCATE(kin21)
 ABI_DEALLOCATE(kin11_k)
 ABI_DEALLOCATE(kin22_k)
 ABI_DEALLOCATE(kin12_k)
 ABI_DEALLOCATE(kin21_k)
 ABI_DEALLOCATE(Stp)
 ABI_DEALLOCATE(Kth)
 ABI_DEALLOCATE(cond_nd)
 ABI_DEALLOCATE(sig_abs)
 ABI_DEALLOCATE(eigen0)
 ABI_DEALLOCATE(nband)
 ABI_DEALLOCATE(oml1)
 ABI_DEALLOCATE(occ)
 ABI_DEALLOCATE(doccde)
 ABI_DEALLOCATE(wtk)

 call hdr%free()

end subroutine conducti_paw
!!***

!!****f* m_conducti/conducti_paw_core
!! NAME
!! conducti_paw_core
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! conductivity tensor and the conductivity along the three principal axes
!! from the Kubo-Greenwood formula for PAW formalism
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

 subroutine conducti_paw_core(filnam,filnam_out,mpi_enreg)

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
 integer :: occopt,rdwr,spaceComm,iunt,opt2_unt,sigx_unt
 real(dp) :: del,diff_occ,ecut,fermie
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
!write(std_out,'(a)')' conducti_paw_core : Please, give the name of the output file ...'
!read(5, '(a)') filnam_out
!write(std_out,'(a)')' The name of the output file is :',trim(filnam_out)
!Read data file
 if (open_file(filnam,msg, newunit=iunt, form='formatted', action="read", status="old") /= 0) then
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
 write(std_out,'(a,i4)') 'selected atom for spectro X',atnbr
 write(std_out,'(a)')'--------------------------------------------'

!These default values are typical of sequential use
 iomode=IO_MODE_FORTRAN; spaceComm=xmpi_comm_self; me=0

! Read the header of the OPT2 file.
 call hdr_read_from_fname(hdr, filnam2, fform2, spaceComm)
 call hdr%free()

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
!
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
!size of the frequency range
 del=(omax-omin)/(mom-1)
 ABI_ALLOCATE(oml1,(mom))
 do iom=1,mom
   oml1(iom)=omin+dble(iom-1)*del
 end do

 ABI_ALLOCATE(sigx,(natom,mom,nphicor))
 ABI_ALLOCATE(sigx_av,(mom,nphicor))
!---------------------------------------------------------------------------------
!SpectroX  -------
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

     do iatom=1,natom
!      first derivative eigenvalues for k-point
       read(opt2_unt) ((psinablapsi2(1:2,1,iband,icor,iatom),iband=1,nband_k),icor=1,nphicor)
       read(opt2_unt) ((psinablapsi2(1:2,2,iband,icor,iatom),iband=1,nband_k),icor=1,nphicor)
       read(opt2_unt) ((psinablapsi2(1:2,3,iband,icor,iatom),iband=1,nband_k),icor=1,nphicor)
     end do

!    LOOP OVER ATOMS/BANDS
     do iatom=1,natom
       do iband=1,nband_k
         do icor=1,nphicor

           do l1=1,3
             dhdk2_g(iatom,iband,icor)=dhdk2_g(iatom,iband,icor)+( &
&             psinablapsi2(1,l1,iband,icor,iatom)*psinablapsi2(1,l1,iband,icor,iatom) &
&             +psinablapsi2(2,l1,iband,icor,iatom)*psinablapsi2(2,l1,iband,icor,iatom))
           end do

           diff_occ = (two/dble(nsppol))-occ_k(iband)
!          Spectro for each omega
           omin = -1.0
           do iom=1,mom
             oml=-energy_cor(icor)+oml1(iom)+omin
             sigx(iatom,iom,icor)=sigx(iatom,iom,icor)+ wtk(ikpt)*dhdk2_g(iatom,iband,icor)&
&             *(diff_occ)/oml*dexp(-((-energy_cor(icor)+eig0_k(iband)-oml)/dom)**2)
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

 if (open_file(trim(filnam_out)//'_sigX', msg, newunit=sigx_unt, form='formatted', action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write(sigx_unt,*) '# conducti: Xray core level conductivity, all in atomic units by default '
 write(sigx_unt,*) '# One block of 3 columns per core wavefunction'
 write(sigx_unt,*) '# energy, sigx_av, sigx, etc... '
 do iom=1,mom
   write(sigx_unt,'( 3(3(1x,e14.8),2x) )') &
&   ((-energy_cor(icor)+oml1(iom)+omin),sigx_av(iom,icor),sigx(atnbr,iom,icor),icor=1,nphicor)
 end do
 close(sigx_unt)

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

 call hdr%free()

end subroutine conducti_paw_core
!!***

!!****f* m_conducti/conducti_nc
!! NAME
!! conducti_nc
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! conductivity tensor and the conductivity along the three principal axes
!! from the Kubo-Greenwood formula.
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! NOTES
!!  bantot
!!  doccde(mband*nkpt_rbz*nsppol)=derivative of occ_rbz wrt the energy.
!!  dom=frequency range
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree).
!!  eigen11(2*mband*mband*nkpt_rbz*nsppol)=first-order eigenvalues (hartree)
!!  in reciprocal direction 100
!!  eigen12(2*mband*mband*nkpt_rbz*nsppol)=first-order eigenvalues (hartree)
!!  in reciprocal direction 010
!!  eigen13(2*mband*mband*nkpt_rbz*nsppol)=first-order eigenvalues (hartree)
!!  in reciprocal direction 001
!!  ecut=kinetic energy planewave cutoff (hartree).
!!  entropy= entropy associated with the smearing (adimensional)
!!  fermie= fermi energy (Hartree)
!!  gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{2}$).
!!  gmet_inv(3,3)=inverse of reciprocal space metric ($\textrm{bohr}^{2}$).
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1).
!!  kin11= Onsager kinetic coeficient=optical conductivity
!!  kin12= Onsager kinetic coeficient
!!  kin21= Onsager kinetic coeficient
!!  kin22= Onsager kinetic coeficient
!!  Kth=thermal conductivity
!!  mom=number of frequency for conductivity computation
!!  mband=maximum number of bands.
!!  natom = number of atoms in the unit cell.
!!  nband(nkpt*nsppol)=number of bands at each RF k point for each spin.
!!  nelect=number of electrons per unit cell
!!  nkpt=number of k points in the IBZ for this perturbation
!!  ngfft(3)=integer fft box dimensions.
!!  nspinor=number of spinorial components of the wavefunctions.
!!  nsppol=1 for unpolarized, 2 for spin-polarized.
!!  ntypat = number of atom types.
!!  occ(mband*nkpt*nsppol)=occupation number for each band and k.
!!  occopt==option for occupancies
!!  rmet(3,3)=real space metric ($\textrm{bohr}^{2}$).
!!  rprimd(3,3)=real space primitive translations.
!!  of primitive translations.
!!  Sth=thermopower
!!  tsmear=smearing width (or temperature) in Hartree
!!  ucvol=unit cell volume in ($\textrm{bohr}^{3}$).
!!  wind=frequency windows for computations of sigma
!!  wtk(nkpt)=weight assigned to each k point.
!!  znucl(natom)=atomic number of atoms
!!  np_sum=noziere-pines sumrule
!!  cond_kg(mom)=kubo-greenwood conductivity
!!
!! PARENTS
!!      conducti
!!
!! CHILDREN
!!      getnel,hdr_free,jacobi,matr3inv,metric,msig,nctk_fort_or_ncfile
!!      wfk_open_read,wfk_read_eigk
!!
!! SOURCE


subroutine conducti_nc(filnam,filnam_out,mpi_enreg)

!Arguments -----------------------------------
!scalars
 character(len=fnlen) :: filnam,filnam_out
 type(MPI_type),intent(in) :: mpi_enreg

!Local variables-------------------------------
!scalars
 integer,parameter :: formeig0=0,formeig1=1
 integer :: bantot,bd2tot_index,bdtot0_index,bdtot_index
 integer :: headform,iband,ii,jj,ikpt,iunt
 integer :: index_1,iom,isppol,jband,l1,l2,mband,mom,natom,nband1
 integer :: nrot,iomode
 integer :: nband_k,nkpt,nlign,nrest,nspinor,nsppol,ntypat
 integer :: occopt,comm
 integer :: tens_unt,lij_unt,sig_unt,kth_unt,ocond_unt
 real(dp) :: deltae,dosdeltae,diff_occ,dom,ecut,entropy,fermie,maxocc
 real(dp) :: nelect,np_sum,np_sum_k1,np_sum_k2,omin,oml,socc,socc_k,sig
 real(dp) :: tphysel,tsmear,ucvol,wind,Tatm
 character(len=fnlen) :: filnam0,filnam1,filnam2,filnam3
 character(len=500) :: msg
 type(hdr_type) :: hdr
 type(wfk_t) :: gswfk,ddk1,ddk2,ddk3
!arrays
 integer,allocatable :: nband(:)
 real(dp) :: gmet(3,3),gmet_inv(3,3),gprimd(3,3),gprimd_inv(3,3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: cond_kg(:,:,:),cond_kg_cart(:,:,:),cond_nd(:,:,:),dhdk2_r(:,:,:,:),dhdk2_g(:,:)
 real(dp),allocatable :: doccde(:),doccde_k(:),cond_kg_xx(:),cond_kg_yy(:),cond_kg_zz(:),trace(:)
 real(dp),allocatable :: eig0_k(:),eig0tmp(:),eig1_k(:,:),eigen0(:),eigen11(:)
 real(dp),allocatable :: eigen12(:),eigtmp(:)
 real(dp),allocatable :: eigen13(:),occ(:),occ_k(:),wtk(:),cond_tot(:),oml1(:)
 real(dp),allocatable :: kin11(:),kin12(:),kin21(:),kin22(:)
 real(dp),allocatable :: kin11_k(:),kin12_k(:),kin21_k(:),kin22_k(:),Kth(:),Stp(:)
 real(dp) :: cond_kg_w(3,3),z(3,3)
 real(dp) :: eig_cond(3)

! *********************************************************************************

 ABI_UNUSED(mpi_enreg%paral_kgb)

!Read data file
 if (open_file(filnam,msg,newunit=iunt,form='formatted', status="old") /= 0) then
   MSG_ERROR(msg)
 end if

 rewind(iunt)
 read(iunt,*)
 read(iunt,'(a)')filnam1       ! first ddk file
 read(iunt,'(a)')filnam2       ! second ddk file
 read(iunt,'(a)')filnam3       ! third ddk file
 read(iunt,'(a)')filnam0       ! ground-state data

!Open the GS Wavefunction file and the 3 DDK files.

! TODO: one should perform basic consistency tests for the GS WFK and the DDK files, e.g.
! k-points and their order, spins, number of bands could differ in the four files.
! Note indeed that we are assuming the same numer of bands in all the files.
 comm = xmpi_comm_self
 call nctk_fort_or_ncfile(filnam0, iomode, msg)
 if (len_trim(msg) /= 0) MSG_ERROR(msg)
 call wfk_open_read(gswfk,filnam0, formeig0, iomode, get_unit(), comm)

 call nctk_fort_or_ncfile(filnam1, iomode, msg)
 if (len_trim(msg) /= 0) MSG_ERROR(msg)
 call wfk_open_read(ddk1,filnam1, formeig1, iomode, get_unit(), comm, hdr_out=hdr)

 call nctk_fort_or_ncfile(filnam2, iomode, msg)
 if (len_trim(msg) /= 0) MSG_ERROR(msg)
 call wfk_open_read(ddk2,filnam2, formeig1, iomode, get_unit(), comm)

 call nctk_fort_or_ncfile(filnam3, iomode, msg)
 if (len_trim(msg) /= 0) MSG_ERROR(msg)
 call wfk_open_read(ddk3,filnam3, formeig1, iomode, get_unit(), comm)

 if (ddk1%compare(ddk2) /= 0) then
   MSG_ERROR("ddk1 and ddk2 are not consistent. see above messages")
 end if
 if (ddk1%compare(ddk3) /= 0) then
   MSG_ERROR("ddk1 and ddk3 are not consistent. see above messages")
 end if

!Extract params from the header of the first ddk file (might have been the GS file ?)

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
 fermie=hdr%fermie
 occ(1:bantot)=hdr%occ(1:bantot)
 nband(1:nkpt*nsppol)=hdr%nband(1:nkpt*nsppol)

!Get mband, as the maximum value of nband(nkpt)
 mband=maxval(nband(:))

 write(std_out,*)
 write(std_out,'(a,3f10.5,a)' )' rprimd(bohr)      =',rprimd(1:3,1)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,2)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,3)
 write(std_out,'(a,i8)')       ' natom             =',natom
 write(std_out,'(a,2i8)')      ' nkpt,mband        =',nkpt,mband
 write(std_out,'(a, f10.5,a)' ) ' ecut              =',ecut,' Ha'
 write(std_out,'(a,f10.5,a,f10.5,a)' )' fermie            =',fermie,' Ha',fermie*Ha_eV,' eV'

!Prepare the reading of ddk Wff files
 ABI_ALLOCATE(eigtmp,(2*mband*mband))
 ABI_ALLOCATE(eig0tmp,(mband))

!Read the eigenvalues of ground-state and ddk files
 ABI_ALLOCATE(eigen0,(mband*nkpt*nsppol))
 ABI_ALLOCATE(eigen11,(2*mband*mband*nkpt*nsppol))
 ABI_ALLOCATE(eigen12,(2*mband*mband*nkpt*nsppol))
 ABI_ALLOCATE(eigen13,(2*mband*mband*nkpt*nsppol))
 bdtot0_index=0 ; bdtot_index=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband1=nband(ikpt+(isppol-1)*nkpt)
     call gswfk%read_eigk(ikpt,isppol,xmpio_single,eig0tmp)
     eigen0(1+bdtot0_index:nband1+bdtot0_index)=eig0tmp(1:nband1)

     call ddk1%read_eigk(ikpt,isppol,xmpio_single,eigtmp)
     eigen11(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)

     call ddk2%read_eigk(ikpt,isppol,xmpio_single,eigtmp)
     eigen12(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)

     call ddk3%read_eigk(ikpt,isppol,xmpio_single,eigtmp)
     eigen13(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)

     bdtot0_index=bdtot0_index+nband1
     bdtot_index=bdtot_index+2*nband1**2
   end do
 end do

 ! Close files
 call gswfk%close()
 call ddk1%close()
 call ddk2%close()
 call ddk3%close()

 ABI_DEALLOCATE(eigtmp)
 ABI_DEALLOCATE(eig0tmp)

!---------------------------------------------------------------------------------
!gmet inversion
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 call matr3inv(gmet,gmet_inv)
 call matr3inv(gprimd,gprimd_inv)

!---------------------------------------------------------------------------------
!derivative of occupation wrt the energy.
 ABI_ALLOCATE(doccde,(mband*nkpt*nsppol))
 ABI_ALLOCATE(wtk,(nkpt))

 read(iunt,*)tsmear
 Tatm=tsmear*Ha_K
 write(std_out,'(a,f12.5,a,f12.5,a)') ' Temp              =',tsmear,' Ha ',Tatm,' Kelvin'
!
 nlign=nkpt/6
 nrest=nkpt-6*nlign
 index_1=0
 do ii=1,nlign
   read(iunt,*)wtk(1+index_1:6+index_1)
   index_1=index_1+6
 end do
 if (nrest/=0) then
   read(iunt,*)wtk(6*nlign+1:nkpt)
 end if
!
 if (occopt==1) then
   write(std_out,'(a,i4)')  ' occopt            =',occopt
   doccde=0.0d0
 else
   tphysel=zero
   maxocc=two/(nsppol*nspinor)
   dosdeltae=zero
   call getnel(doccde,dosdeltae,eigen0,entropy,fermie,maxocc,mband,nband,&
&   nelect,nkpt,nsppol,occ,occopt,1,tphysel,tsmear,11,wtk)
!  DEBUG
!  write(std_out,'(a,f10.5)')' getnel : nelect   =',nelect
!  ENDDEBUG
 end if
!---------------------------------------------------------------------------------
!size of the frequency range
 read(iunt,*)dom,wind
 close(iunt)
 mom=int(wind/dom)
 ABI_ALLOCATE(oml1,(mom))
 do iom=1,mom
   oml1(iom)=tol10*1000.0d0+dble(iom)*dom
 end do

 ABI_ALLOCATE(cond_nd,(mom,3,3))
 ABI_ALLOCATE(cond_kg,(mom,3,3))
 ABI_ALLOCATE(cond_kg_cart,(mom,3,3))
 ABI_ALLOCATE(cond_kg_xx,(mom))
 ABI_ALLOCATE(cond_kg_yy,(mom))
 ABI_ALLOCATE(trace,(mom))
 ABI_ALLOCATE(cond_kg_zz,(mom))
 ABI_ALLOCATE(cond_tot,(mom))
 ABI_ALLOCATE(kin11,(mom))
 ABI_ALLOCATE(kin12,(mom))
 ABI_ALLOCATE(kin21,(mom))
 ABI_ALLOCATE(kin22,(mom))
 ABI_ALLOCATE(kin11_k,(mom))
 ABI_ALLOCATE(kin12_k,(mom))
 ABI_ALLOCATE(kin21_k,(mom))
 ABI_ALLOCATE(kin22_k,(mom))
 ABI_ALLOCATE(Kth,(mom))
 ABI_ALLOCATE(Stp,(mom))
 write(std_out,'(a,i8,2f10.5,a)')' mom,wind,dom      =',mom,wind,dom,' Ha'

!---------------------------------------------------------------------------------

 kin11   = 0.0d0
 kin12   = 0.0d0
 kin21   = 0.0d0
 kin22   = 0.0d0
 np_sum  = 0.0d0
 socc    = 0.0d0
 cond_kg = 0.0d0


!LOOP OVER SPINS
 do isppol=1,nsppol
!
   bdtot_index = 0
   bd2tot_index = 0
!
   deltae  = 0.0d0
!
!  BIG FAT k POINT LOOP
!
   do ikpt=1,nkpt
!
     nband_k=nband(ikpt+(isppol-1)*nkpt)
!
     ABI_ALLOCATE(eig0_k,(nband_k))
     ABI_ALLOCATE(eig1_k,(2*nband_k**2,3))
     ABI_ALLOCATE(occ_k,(nband_k))
     ABI_ALLOCATE(doccde_k,(nband_k))
     ABI_ALLOCATE(dhdk2_r,(nband_k,nband_k,3,3))
     ABI_ALLOCATE(dhdk2_g,(nband_k,nband_k))

     cond_nd   = 0.0d0
     kin11_k   = 0.0d0
     kin12_k   = 0.0d0
     kin21_k   = 0.0d0
     kin22_k   = 0.0d0
     np_sum_k1 = 0.0d0
     np_sum_k2 = 0.0d0
     socc_k    = 0.0d0
     dhdk2_r   = 0.0d0
     dhdk2_g   = 0.0d0
!
!    eigenvalue for k-point
     eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
!    first derivative eigenvalues for k-point
     eig1_k(:,1)=eigen11(1+bd2tot_index:2*nband_k**2+bd2tot_index)
     eig1_k(:,2)=eigen12(1+bd2tot_index:2*nband_k**2+bd2tot_index)
     eig1_k(:,3)=eigen13(1+bd2tot_index:2*nband_k**2+bd2tot_index)
!    occupation numbers for k-point
     occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
!    derivative of occupation number for k-point
     doccde_k(:)=doccde(1+bdtot_index:nband_k+bdtot_index)
!
!    DEBUG
!    write(16,*)
!    write(16,*)' conducti : ikpt=',ikpt
!    do iband=1,nband_k
!    write(16, '(i4,4es22.12)' )iband,wtk(ikpt),occ_k(iband),&
!    &                            doccde_k(iband),eig0_k(iband)
!    end do
!    write(16,*)
!    ENDDEBUG
!
!    LOOP OVER BAND
     do iband=1,nband_k
       do jband=1,nband_k
!
! TODO : replace with BLAS calls
         do l1=1,3
           do l2=1,3
             do ii=1,3
               do jj=1,3
                 dhdk2_r(iband,jband,l1,l2)=dhdk2_r(iband,jband,l1,l2)+(rprimd(l1,ii)&
&                 *eig1_k(2*iband-1+(jband-1)*2*nband_k,ii)*&
&                 rprimd(l2,jj)*eig1_k(2*iband-1+(jband-1)*2*nband_k,jj)&
&                 +rprimd(l1,ii)*eig1_k(2*iband  +(jband-1)*2*nband_k,ii)*&
&                 rprimd(l2,jj)*eig1_k(2*iband+(jband-1)*2*nband_k,jj))
               end do
             end do
           end do
         end do

         do l1=1,3
           do l2=1,3
             dhdk2_r(iband,jband,l1,l2)=dhdk2_r(iband,jband,l1,l2)/two_pi/two_pi
           end do
         end do
!
! TODO: replace with BLAS calls
         do l1=1,3
           do l2=1,3
             dhdk2_g(iband,jband)=dhdk2_g(iband,jband)+gmet_inv(l1,l2)*( &
&             eig1_k(2*iband-1+(jband-1)*2*nband_k,l1)*&
&             eig1_k(2*iband-1+(jband-1)*2*nband_k,l2) &
&            +eig1_k(2*iband  +(jband-1)*2*nband_k,l1)*&
&             eig1_k(2*iband  +(jband-1)*2*nband_k,l2))
           end do
         end do
         dhdk2_g(iband,jband)=dhdk2_g(iband,jband)/two_pi/two_pi
!
         diff_occ = occ_k(iband)-occ_k(jband)
!        if (dabs(diff_occ)>=tol8) then
!
!        Conductivity for each omega
         omin = 0.0d0
         do iom=1,mom
           oml=oml1(iom)
           if (jband>iband) then
             sig= dhdk2_g(iband,jband)&
&             *(diff_occ)/oml*(dexp(-((eig0_k(jband)-eig0_k(iband)-oml)/dom)**2)&
&             -dexp(-((eig0_k(iband)-eig0_k(jband)-oml)/dom)**2))
             kin11_k(iom)=kin11_k(iom)+sig
             kin12_k(iom)=kin12_k(iom)-sig*(eig0_k(jband)-fermie)
             kin21_k(iom)=kin21_k(iom)-sig*(eig0_k(iband)-fermie)
             kin22_k(iom)=kin22_k(iom) + &
&             sig*(eig0_k(iband)-fermie)*(eig0_k(jband)-fermie)
           end if
           do l1=1,3
             do l2=1,3
               cond_nd(iom,l1,l2)=cond_nd(iom,l1,l2) +dhdk2_r(iband,jband,l1,l2)&
&               *(diff_occ)/oml*dexp(-((eig0_k(jband)-eig0_k(iband)-oml)/dom)**2)
             end do
           end do

         end do
!
!        Sumrule start
         if (dabs(eig0_k(iband)-eig0_k(jband))>=tol10) then
           np_sum_k1=np_sum_k1 -dhdk2_g(iband,jband)&
&           *(diff_occ)/(eig0_k(iband)-eig0_k(jband))
         else
           np_sum_k2=np_sum_k2 - doccde_k(iband)*dhdk2_g(iband,jband)
         end if
!

!        end loop over band
!        end if
       end do
       socc_k=socc_k+occ_k(iband)
     end do
!
     do iom=1,mom
       kin11(iom)=kin11(iom)+wtk(ikpt)*kin11_k(iom)
       kin12(iom)=kin12(iom)+wtk(ikpt)*kin12_k(iom)
       kin21(iom)=kin21(iom)+wtk(ikpt)*kin21_k(iom)
       kin22(iom)=kin22(iom)+wtk(ikpt)*kin22_k(iom)
       do l1=1,3
         do l2=1,3
           cond_kg(iom,l1,l2)=cond_kg(iom,l1,l2)+wtk(ikpt)*cond_nd(iom,l1,l2)
         end do
       end do
     end do

     np_sum=np_sum + wtk(ikpt)*(np_sum_k1+np_sum_k2)
     socc=socc+wtk(ikpt)*socc_k
!
!    validity limit
     deltae=deltae+(eig0_k(nband_k)-fermie)

     bd2tot_index=bd2tot_index+2*nband_k**2
     bdtot_index=bdtot_index+nband_k
     ABI_DEALLOCATE(eig0_k)
     ABI_DEALLOCATE(eig1_k)
     ABI_DEALLOCATE(occ_k)
     ABI_DEALLOCATE(doccde_k)
     ABI_DEALLOCATE(dhdk2_r)
     ABI_DEALLOCATE(dhdk2_g)
!    End loop over k
   end do

   write(std_out,'(a,3f10.5)')' sumrule           =',np_sum/socc/three,socc
   write(std_out,'(a,f10.5,a,f10.5,a)')&
&   ' Emax-Efermi       =',deltae/dble(nkpt),' Ha',deltae/dble(nkpt)*Ha_eV,' eV'

!  End loop over spins
 end do

 cond_kg=cond_kg*two_pi*third/(dom*ucvol)*half/dsqrt(pi)


!Check that new output file does NOT exist
!Keep this line : prevent silly (compiler ?) bug on HP 8000
 write(std_out,*)' conducti : call isfile '
!
 if (open_file(trim(filnam_out)//'_tens',msg,newunit=tens_unt,form='formatted',action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 if (open_file(trim(filnam_out)//'_Lij',msg,newunit=lij_unt,form='formatted',action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write(lij_unt,'(a)')' # omega(ua) L12 L21 L22 L22'

 if (open_file(trim(filnam_out)//'_sig',msg,newunit=sig_unt,form='formatted',action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write(sig_unt,'(a)')' # omega(ua) hbar*omega(eV)    cond(ua)             cond(ohm.cm)-1'

 if (open_file(trim(filnam_out)//'_Kth',msg,newunit=kth_unt,form='formatted',action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write(kth_unt,'(a)')&
& ' #omega(ua) hbar*omega(eV)  thermal cond(ua)   Kth(W/m/K)   thermopower(ua)   Stp(microohm/K)'

 if (open_file(trim(filnam_out)//'.out',msg,newunit=ocond_unt,form='formatted',action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write(ocond_unt,'(a)' )' Conducti output file:'
 write(ocond_unt,'(a)' )' Contains all results produced by conducti utility'
 write(ocond_unt,'(a)' )' '
 write(ocond_unt,'(a)')' # omega(ua)       cond(ua)             thermal cond(ua)       thermopower(ua)'
!
!call isfile(filnam_out,'new')

!Keep this line : prevent silly (compiler ?) bug on HP 8000
 write(std_out,*)' conducti : after call isfile '
!
!Compute thermal conductivity and thermopower
 do iom=1,mom
   oml=oml1(iom)
   kin11(iom)=kin11(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
   kin21(iom)=kin21(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
   kin12(iom)=kin12(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
   kin22(iom)=kin22(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
   if (dabs(kin11(iom))<10.0d-20) kin11(iom)=zero
   Kth(iom)=kin22(iom)
   Stp(iom)=zero
   if(kin11(iom)/=zero)  then
     Kth(iom)=Kth(iom)-(kin12(iom)*kin21(iom)/kin11(iom))
     Stp(iom)=kin12(iom)/(kin11(iom)*Tatm)
   end if
   if (dabs(Kth(iom))<10.0d-20) Kth(iom)=0.0d0
   if (dabs(Stp(iom))<10.0d-20) Stp(iom)=0.0d0
   if (dabs(kin12(iom))<10.0d-20) kin12(iom)=zero
   if (dabs(kin21(iom))<10.0d-20) kin21(iom)=zero
   if (dabs(kin22(iom))<10.0d-20) kin22(iom)=zero

   write(lij_unt,'(f12.5,4es22.12)')oml,kin12(iom),kin21(iom),kin22(iom),kin22(iom)/Tatm*3.4057d9
   write(sig_unt,'(2f12.5,2es22.12)') oml,oml*Ha_eV,kin11(iom),kin11(iom)*Ohmcm
   write(kth_unt,'(2f12.5,4es22.12)') oml,oml*Ha_eV,Kth(iom),Kth(iom)*3.4057d9/Tatm,Stp(iom),Stp(iom)*3.6753d-2
   write(ocond_unt,'(1f12.5,3es22.12)') oml,kin11(iom),Kth(iom),Stp(iom)
 end do
!

 write(tens_unt,'(a)' )' Conductivity file '
 write(tens_unt,'(a)' )' ----------------- '
 write(tens_unt,'(a)' )' Contain first the full conductivity tensor, for the desired set of energies,'
 write(tens_unt,'(a)' )' then, the three principal values, for the desired set of energies'
 write(tens_unt,'(a)' )' (note that eigenvalues are not directly associated with xx,yy,zz)'
 write(tens_unt,'(a)' )' '

 write(ocond_unt,'(a)' )' '
 write(ocond_unt,'(a)' )' full conductivity tensor, for the desired set of energies'
 write(ocond_unt,'(a)' )' then, the three principal values, for the desired set of energies:'

 do iom=1,mom
   oml=oml1(iom)*Ha_eV
   write(tens_unt, '(a,es16.6,a)' ) ' energy (in eV) =',oml,', conductivity tensor (in Ohm.cm-1) follows :'
   write(ocond_unt, '(a,es16.6,a)' ) ' energy (in eV) =',oml,', conductivity tensor (in Ohm.cm-1) follows :'
   do l1=1,3
     write(tens_unt,"(3f25.15)") (cond_kg(iom,l1,l2)*Ohmcm,l2=1,3)
     write(ocond_unt,"(3f25.15)") (cond_kg(iom,l1,l2)*Ohmcm,l2=1,3)
   end do
 end do

!Diagonalizing the conductivity matrix for sigma_xx,sigma_yy,sigma_zz
 cond_kg_xx=0d0
 cond_kg_yy=0d0
 cond_kg_zz=0d0
!trace=0d0    used for checking with the original version of the code
 do iom=1,mom
   oml=oml1(iom)*Ha_eV
   cond_kg_w=0d0
   do l1=1,3
     do l2=1,3
       cond_kg_w(l1,l2)=cond_kg(iom,l1,l2)
     end do
   end do
   call jacobi(cond_kg_w,3,3,eig_cond,z,nrot)

!  When the value is too small, set it to zero before printing
   if(abs(eig_cond(1))<tol10)eig_cond(1)=zero
   if(abs(eig_cond(2))<tol10)eig_cond(2)=zero
   if(abs(eig_cond(3))<tol10)eig_cond(3)=zero

   cond_kg_xx(iom)=eig_cond(1)
   cond_kg_yy(iom)=eig_cond(2)
   cond_kg_zz(iom)=eig_cond(3)
!  trace(iom)=cond_kg_xx(iom)+cond_kg_yy(iom)+cond_kg_zz(iom)
 end do

!DEBUG Keep this line : prevent silly (compiler ?) bug on HP 8000
!write(std_out,*)' conducti : after open '
!ENDDEBUG

 write(tens_unt,'(a,a)')ch10,' Now, print principal values of the conductivity tensor.'
 write(tens_unt,'(a)')' '
 write(tens_unt,'(a)')' #omega(ua)   cond_1(ua)     cond_2(ua) cond_3(ua)  cond_tot(ua)'

 write(ocond_unt,'(a)')' '
 write(ocond_unt,'(a,a)')ch10,' Now, print principal values of the conductivity tensor.'
 write(ocond_unt,'(a)')' '
 write(ocond_unt,'(a)')' #omega(ua)   cond_1(ua)     cond_2(ua) cond_3(ua)  cond_tot(ua)'


 do iom=1,mom
   cond_tot(iom)=cond_kg_xx(iom)+cond_kg_yy(iom)+cond_kg_zz(iom)
   write(tens_unt,'(f12.5,4es22.12)')oml1(iom),cond_kg_xx(iom),cond_kg_yy(iom),cond_kg_zz(iom),cond_tot(iom)
   write(ocond_unt,'(f12.5,4es22.12)')oml1(iom),cond_kg_xx(iom),cond_kg_yy(iom),cond_kg_zz(iom),cond_tot(iom)
 end do

 write(tens_unt,*)
 write(tens_unt,'(a)')' #hbar*omega(eV)    cond_1(ohm.cm)-1    cond_2(ohm.cm)-1    cond_3(ohm.cm)-1    cond_t(ohm.cm)-1'
 write(ocond_unt,*)
 write(ocond_unt,'(a)')' #hbar*omega(eV)    cond_1(ohm.cm)-1    cond_2(ohm.cm)-1    cond_3(ohm.cm)-1    cond_t(ohm.cm)-1'

 do iom=1,mom
   oml=oml1(iom)*Ha_eV
   cond_tot(iom)=cond_tot(iom)*Ohmcm
   cond_kg_xx(iom)=cond_kg_xx(iom)*Ohmcm
   cond_kg_yy(iom)=cond_kg_yy(iom)*Ohmcm
   cond_kg_zz(iom)=cond_kg_zz(iom)*Ohmcm
   write(tens_unt,'(f12.5,4es22.12)')oml,cond_kg_xx(iom),cond_kg_yy(iom),cond_kg_zz(iom),cond_tot(iom)
   write(ocond_unt,'(f12.5,4es22.12)')oml,cond_kg_xx(iom),cond_kg_yy(iom),cond_kg_zz(iom),cond_tot(iom)
 end do
!Calculate the imaginary part of the conductivity (principal value)
!+derived optical properties.

 call msig (kin11,mom,oml1,filnam_out)

 close(tens_unt)
 close(lij_unt)
 close(sig_unt)
 close(kth_unt)
 close(ocond_unt)

 ABI_DEALLOCATE(nband)
 ABI_DEALLOCATE(oml1)
 ABI_DEALLOCATE(occ)
 ABI_DEALLOCATE(eigen11)
 ABI_DEALLOCATE(eigen12)
 ABI_DEALLOCATE(eigen13)
 ABI_DEALLOCATE(eigen0)
 ABI_DEALLOCATE(doccde)
 ABI_DEALLOCATE(wtk)
 ABI_DEALLOCATE(cond_nd)
 ABI_DEALLOCATE(cond_kg)
 ABI_DEALLOCATE(cond_kg_cart)
 ABI_DEALLOCATE(cond_kg_xx)
 ABI_DEALLOCATE(cond_kg_yy)
 ABI_DEALLOCATE(trace)
 ABI_DEALLOCATE(cond_kg_zz)
 ABI_DEALLOCATE(cond_tot)
 ABI_DEALLOCATE(kin11)
 ABI_DEALLOCATE(kin22)
 ABI_DEALLOCATE(kin12)
 ABI_DEALLOCATE(kin21)
 ABI_DEALLOCATE(kin11_k)
 ABI_DEALLOCATE(kin22_k)
 ABI_DEALLOCATE(kin12_k)
 ABI_DEALLOCATE(kin21_k)
 ABI_DEALLOCATE(Stp)
 ABI_DEALLOCATE(Kth)

 call hdr%free()

 end subroutine conducti_nc
!!***

!!****f* m_conducti/msig
!! NAME
!! msig
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! conductivity tensor and the conductivity along the three principal axes
!! from the Kubo-Greenwood formula for PAW formalism
!!
!! INPUTS
!!  fcti(npti)=  conductivity, as calculated in conducti
!!  npti= number of points to calculate conductivity
!!  xi(npti)= energies where the conductivity is calculated
!!
!! OUTPUT
!!   no output, only files
!!
!! NOTES
!!     this program calculates the imaginary part of the conductivity (principal value)
!!     +derived optical properties.
!!     the calculation is performed on the same grid as the initial input
!!     to calculate the principal value, a trapezoidale integration +taylor expansion to
!!     third order is used (W.J. Thomson computer in physics vol 12 p94 1998)
!!    two input files are needed inppv.dat (parameters) and sigma.dat (energy,sigma_1)
!!     two output files ppsigma.dat (energy,sigma_1,sigma_2,epsilon_1,epsilon_2)
!!                      abs.dat     (energy,nomega,komega,romega,absomega)
!!     march 2002 s.mazevet
!!
!! PARENTS
!!      conducti_nc,conducti_paw
!!
!! CHILDREN
!!      intrpl
!!
!! SOURCE

subroutine msig(fcti,npti,xi,filnam_out_sig)

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: npti
!arrays
 real(dp),intent(in) :: fcti(npti),xi(npti)
 character(len=fnlen),intent(in) :: filnam_out_sig

!Local variables-------------------------------
!scalars
 integer :: npt = 10000
 integer :: ii,ip,npt1,npt2,eps_unt,abs_unt
 real(dp),parameter :: del=0.001_dp,ohmtosec=9.d11
 real(dp) :: dx,dx1,dx2,eps1,eps2,idel,komega,pole,refl,sigma2,xsum
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: abso(:),fct(:),fct1(:),fct2(:),fct3(:),fct4(:),fct5(:)
 real(dp),allocatable :: fctii(:),fp(:),fpp(:),fppp(:),nomega(:),ppsig(:)
 real(dp),allocatable :: x1(:),x2(:)

! *********************************************************************************
!BEGIN EXECUTABLE SECTION

 if (npti > 12000) then
   msg = "Sorry - the interpolator INTRPL is hard coded for maximum 12000 points." // &
&        ch10 // " Reduce the conducti input npti, or implement a better interpolator!"
   MSG_ERROR(msg)
 end if

 write(std_out,'(2a)')ch10,'Calculate the principal value and related optical properties'
 write(std_out,'(a)')'following W.J. Thomson computer in physics vol 12 p94 1998 for '
 write(std_out,'(a)')'the principal value. S. Mazevet'
 write(std_out,'(a)')'OPTIONS'
 write(std_out,'(a)')'use default number of integration pts: npt=10000'
 write(std_out,'(a)')'Use default value for delta interval: del=1e-3'

 if (open_file(trim(filnam_out_sig)//'_eps',msg, newunit=eps_unt,status='replace',action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write(eps_unt,'(a)')'#energy (eV),sigma_1(Ohm-1cm-1),sigma_2(Ohm-1cm-1),epsilon_1,epsilon_2'

 if (open_file(trim(filnam_out_sig)//'_abs', msg, newunit=abs_unt, status='replace',action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write(abs_unt,'(a)')'#energy(eV),nomega,komega,refl.,abso.(cm-1)'

 ABI_ALLOCATE(fct,(npt))
 ABI_ALLOCATE(fct2,(npt))
 ABI_ALLOCATE(fct3,(npt))
 ABI_ALLOCATE(fct4,(npt))
 ABI_ALLOCATE(fct5,(npt))
 ABI_ALLOCATE(fp,(npt))
 ABI_ALLOCATE(fpp,(npt))
 ABI_ALLOCATE(fppp,(npt))
 ABI_ALLOCATE(x1,(npt))
 ABI_ALLOCATE(x2,(npt))
 ABI_ALLOCATE(fct1,(npt))
 ABI_ALLOCATE(ppsig,(npt))
 ABI_ALLOCATE(fctii,(npt))
 ABI_ALLOCATE(abso,(npt))
 ABI_ALLOCATE(nomega,(npt))

!loop on the initial energy grid
 do ip=1,npti

!  adjust the interval before and after the pole to reflect range/npt interval
   xsum=zero
   dx=(xi(npti)-xi(1))/dble(npt-1)
   pole=xi(ip)
   npt1=int((pole-del)/dx)
   dx1=zero
   if(npt1/=1) dx1=(pole-del)/(npt1-1)
   npt2=int((xi(npti)-pole-del)/dx)
   dx2=(xi(npti)-pole-del)/(npt2-1)

!  for the moment skip the pp calculation when the pole if too close to the end of the range
   if(npt1<=1.or.npt2<=1) then

     xsum=zero
     ppsig(ip)=zero

   else

!    define the fct for which the pp calculation is needed using xi^2-pole^2 factorization
     fctii(1:npti) = zero
     fctii(1:npti)=fcti(1:npti)*pole/(xi(1:npti)+pole)

!    define the grid on each side of the pole x1 before x2 after
     do ii=1,npt1
       x1(ii)=dx1*dble(ii-1)
     end do
     do ii=1,npt2
       x2(ii)=pole+del+dx2*dble(ii-1)
     end do

!    interpolate the initial fct fii on the new grids x1 and x2 (cubic spline)
!    write(std_out,*) npti,npt1

!    MJV 6/12/2008:
!    for each use of fctii should ensure that npt1 npt2 etc... are less than
!    npt=len(fctii)
! TODO: move to spline/splint routines with no memory limitation
     call intrpl(npti,xi,fctii,npt1,x1,fct4,fct1,fct5,1)
     call intrpl(npti,xi,fctii,npt2,x2,fct3,fct2,fct5,1)

!    calculate the two integrals from 0-->pole-lamda and pole+lamda--> end range
!    trapezoidal integration
     do ii=1,npt1
       fct1(ii)=fct4(ii)/(x1(ii)-pole)
     end do
     do ii=1,npt2
       fct2(ii)=fct3(ii)/(x2(ii)-pole)
     end do

     do ii=2,npt1-1
       xsum=xsum+fct1(ii)*dx1
     end do
     do ii=2,npt2-1
       xsum=xsum+fct2(ii)*dx2
     end do
     xsum=xsum+half*(fct1(1)+fct1(npt1))*dx1+half*(fct2(1)+fct2(npt2))*dx2

!    calculate the first and third derivative at the pole and add the taylor expansion
! TODO: move to spline/splint routines with no memory limitation
     call intrpl(npti,xi,fctii,npti,xi,fct3,fct4,fct5,1)
     call intrpl(npti,xi,fct4,1,(/pole/),fp,fpp,fppp,1)

     idel=two*fp(1)*(del)+fppp(1)*(del**3)/nine
     xsum=xsum+idel

   end if

!  calculate the derivated optical quantities and output the value
   sigma2=(-two/pi)*xsum
   eps1=one-(four_pi*sigma2/(pole))
   eps2=four*fcti(ip)*pi/(pole)

!  A special treatment of the case where eps2 is very small compared to eps1 is needed
   if(eps2**2 > eps1**2 * tol12)then
     nomega(ip)=sqrt(half*(eps1 + sqrt(eps1**2 + eps2**2)))
     komega=sqrt(half*(-eps1 + sqrt(eps1**2 + eps2**2)))
     abso(ip)=four_pi*fcti(ip)*ohmtosec*Ohmcm/nomega(ip)/(Sp_Lt_SI*100._dp)
   else if(eps1>zero)then
     nomega(ip)=sqrt(half*(eps1 + sqrt(eps1**2 + eps2**2)))
     komega=half*abs(eps2/sqrt(eps1))
     abso(ip)=four_pi*fcti(ip)*ohmtosec*Ohmcm/nomega(ip)/(Sp_Lt_SI*100._dp)
   else if(eps1<zero)then
     nomega(ip)=half*abs(eps2/sqrt(-eps1))
     komega=sqrt(half*(-eps1 + sqrt(eps1**2 + eps2**2)))
     abso(ip)=two*sqrt(-eps1)*pole*ohmtosec*Ohmcm/(Sp_Lt_SI*100._dp)
   end if

   refl=((one-nomega(ip))**2 + komega**2)/ &
&   ((one+nomega(ip))**2 + komega**2)

   write(eps_unt,'(5e18.10)') Ha_eV*pole,fcti(ip)*Ohmcm,sigma2*Ohmcm,eps1,eps2
   write(abs_unt,'(5e18.10)') Ha_eV*pole,nomega(ip),komega,refl,abso(ip)

 end do

 close(eps_unt)
 close(abs_unt)

 ABI_DEALLOCATE(fct)
 ABI_DEALLOCATE(x1)
 ABI_DEALLOCATE(x2)
 ABI_DEALLOCATE(fct2)
 ABI_DEALLOCATE(fct3)
 ABI_DEALLOCATE(fct4)
 ABI_DEALLOCATE(fct5)
 ABI_DEALLOCATE(fp)
 ABI_DEALLOCATE(fpp)
 ABI_DEALLOCATE(fppp)
 ABI_DEALLOCATE(fct1)
 ABI_DEALLOCATE(ppsig)
 ABI_DEALLOCATE(fctii)
 ABI_DEALLOCATE(abso)
 ABI_DEALLOCATE(nomega)

end subroutine msig
!!***

!!****f* ABINIT/emispec_paw
!! NAME
!! emispec_paw
!!
!! FUNCTION
!! This program computes the elements of the emission spectra
!! from the Kubo-Greenwood formula for PAW formalism
!! largely inspired from the conducti_core_paw routine
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

 subroutine emispec_paw(filnam,filnam_out,mpi_enreg)

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
 call hdr%free()

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
 write(ems_unt,*) '# conducti: Xray emission spectrum, all in atomic units by default '
 write(ems_unt,*) '# One block of 3 columns per core wavefunction'
 write(ems_unt,*) '# energy, sigx_av, sigx, etc... '
 do iom=1,mom
   write(ems_unt,'( 3(3(1x,e15.8),2x) )') &
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

 call hdr%free()

end subroutine emispec_paw
!!***

end module m_conducti
!!***
