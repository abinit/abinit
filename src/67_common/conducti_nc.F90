!{\src2tex{textfont=tt}}
!!****f* ABINIT/conducti_nc
!! NAME
!! conducti_nc
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! conductivity tensor and the conductivity along the three principal axes
!! from the Kubo-Greenwood formula.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2018 ABINIT group (VRecoules, PGhosh)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!!      wfk_close,wfk_open_read,wfk_read_eigk
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine conducti_nc(filnam,filnam_out,mpi_enreg)

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_profiling_abi
 use m_errors
 use m_wfk
 use m_nctk
 use m_hdr

 use m_io_tools,     only : open_file, get_unit
 use m_geometry,     only : metric
 use m_fstrings,     only : sjoin
 use m_abilasi,      only : jacobi
 use m_occ,          only : getnel

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'conducti_nc'
 use interfaces_32_util
 use interfaces_61_occeig
 use interfaces_67_common, except_this_one => conducti_nc
!End of the abilint section

 implicit none

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

 if (wfk_compare(ddk1, ddk2) /= 0) then
   MSG_ERROR("ddk1 and ddk2 are not consistent. see above messages")
 end if
 if (wfk_compare(ddk1, ddk3) /= 0) then
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
     call wfk_read_eigk(gswfk,ikpt,isppol,xmpio_single,eig0tmp)
     eigen0(1+bdtot0_index:nband1+bdtot0_index)=eig0tmp(1:nband1)

     call wfk_read_eigk(ddk1,ikpt,isppol,xmpio_single,eigtmp)
     eigen11(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)

     call wfk_read_eigk(ddk2,ikpt,isppol,xmpio_single,eigtmp)
     eigen12(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)

     call wfk_read_eigk(ddk3,ikpt,isppol,xmpio_single,eigtmp)
     eigen13(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)

     bdtot0_index=bdtot0_index+nband1
     bdtot_index=bdtot_index+2*nband1**2
   end do
 end do

 ! Close files
 call wfk_close(gswfk)
 call wfk_close(ddk1)
 call wfk_close(ddk2)
 call wfk_close(ddk3)

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
 mom=dint(wind/dom)
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
         do l1=1,3
           do l2=1,3
             dhdk2_g(iband,jband)=dhdk2_g(iband,jband)+gmet_inv(l1,l2)*( &
&             eig1_k(2*iband-1+(jband-1)*2*nband_k,l1)*&
&             eig1_k(2*iband-1+(jband-1)*2*nband_k,l2) &
&             +eig1_k(2*iband  +(jband-1)*2*nband_k,l1)*&
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

 call hdr_free(hdr)

 end subroutine conducti_nc
!!***
