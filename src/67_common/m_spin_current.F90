!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spin_current
!! NAME
!!  m_spin_current
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (Mver)
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

module m_spin_current

 use defs_datatypes
 use m_errors
 use m_abicore
 use m_splines
 use m_hdr
 use m_dtset
 use m_dtfil


 use defs_abitypes, only : MPI_type
 use m_io_tools,   only : open_file
 use m_pptools,    only : printxsf
 use m_geometry,   only : xred2xcart
 use m_fftcore,    only : sphereboundary
 use m_special_funcs,   only : gamma_function
 use m_fft,            only : fourwf

 implicit none

 private
!!***

 public :: spin_current
!!***

contains
!!***

!!****f* m_spin_current/spin_current
!! NAME
!! spin_current
!!
!! FUNCTION
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  atindx1(natom)=inverse of atindx
!!  cg(2,mcg)=wavefunctions (may be read from disk instead of input)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  gmet = reciprocal space metric
!!  gprimd = dimensionful reciprocal space vectors
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  kg(3,mpw*mkmem)=reduced (integer) coordinates of G vecs in basis sphere
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mpi_enreg=information about MPI parallelization
!!  nattyp(dtset%ntypat)=number of atoms of each type
!!  nfftf = fft grid dimensions for fine grid
!!  ph1d = phase factors in 1 radial dimension
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum
!!  rhog(2,nfftf)=Fourier transform of total electron density (including compensation density in PAW)
!!  rhor(nfftf,nspden)=total electron density (including compensation density in PAW)
!!  rmet = real space metric tensor
!!  symrec(3,3,nsym)=symmetries in reciprocal space, reduced coordinates
!!  ucvol = unit cell volume
!!  wffnow=unit number for current wf disk file
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!   only output to file
!!
!! PARENTS
!!      afterscfloop
!!
!! CHILDREN
!!      fourwf,printxsf,sphereboundary,vso_realspace_local,xred2xcart
!!
!! SOURCE

subroutine spin_current(cg,dtfil,dtset,gprimd,hdr,kg,mcg,mpi_enreg,psps)

!Arguments ------------------------------------
!scalars
!integer,intent(in) :: nfftf
!real(dp),intent(in) :: ucvol
 integer,intent(in) :: mcg
 type(MPI_type),intent(in) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(pseudopotential_type),intent(in) :: psps
!type(wffile_type),intent(in) :: wffnow
!arrays
!integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
!integer,intent(in) :: nattyp(dtset%ntypat)
!integer,intent(in) :: symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: cg(2,mcg)
!real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol),gmet(3,3)
 real(dp),intent(in) :: gprimd(3,3)
!real(dp),intent(in) :: rhog(2,nfftf),rhor(nfftf,dtset%nspden)
!real(dp),intent(in) :: rmet(3,3)
!real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
!real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
!real(dp),intent(inout) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)

!Local variables-------------------------------
!scalars
 integer :: cplex,fft_option,i1
 integer :: i2,i3,iband,icartdir,icg,ig
 integer :: ikg,ikpt,iocc,irealsp,ispindir,ispinor,ispinorp
 integer :: npw
 integer :: icplex
 integer :: realrecip
 integer :: iatom,spcur_unit
 real(dp) :: prefact_nk
 real(dp) :: rescale_current
 character(len=500) :: message
 character(len=fnlen) :: filnam
!arrays
 integer,allocatable :: gbound(:,:),kg_k(:,:)
 real(dp),allocatable :: dpsidr(:,:,:,:,:,:)
 real(dp),allocatable :: density(:,:,:,:)
 real(dp),allocatable :: dummy_denpot(:,:,:)
 real(dp),allocatable :: gpsi(:,:,:,:),kgcart(:,:)
 real(dp),allocatable :: position_op(:,:,:,:)
 real(dp),allocatable :: psi(:,:,:),psi_r(:,:,:,:,:)
 real(dp),allocatable :: spincurrent(:,:,:,:,:)
 real(dp),allocatable :: vso_realspace(:,:,:,:,:),datagrid(:)
 real(dp) :: dummy_fofgout(0,0)
 real(dp),allocatable :: xcart(:,:)
 character :: spin_symbol(3)
 character :: spinor_sym(2)
 character(len=2) :: realimag(2)
!no_abirules
!real(dp),allocatable :: density_matrix(:,:,:,:,:)
!real(dp),allocatable :: vso_realspace_nl(:,:,:,:,:)

! *************************************************************************

!write(std_out,*) ' Entering subroutine spin_current '
!write(std_out,*) ' dtset%ngfft = ', dtset%ngfft
!write(std_out,*) ' hdr%istwfk = ', hdr%istwfk

!===================== init and checks ============================
!check if nspinor is 2
 if (dtset%nspinor /= 2) then
   write(message, '(a,i0)' )' nspinor must be 2, but it is ',dtset%nspinor
   MSG_ERROR(message)
 end if

 if (dtset%nsppol /= 1) then
   write(message, '(a,i0)' )' spin_current:  nsppol must be 1 but it is ',dtset%nsppol
   MSG_ERROR(message)
 end if

 if (dtset%mkmem /= dtset%nkpt) then
   write(message, '(a,i6,a,i6,a,a)' )&
&   ' mkmem =  ',dtset%mkmem,' must be equal to nkpt ',dtset%nkpt,ch10,&
&   ' keep all kpt in memory'
   MSG_ERROR(message)
 end if

 if (dtset%usepaw /= 0) then
   write(message, '(a,i0,a,a,a)' )&
&   'usepaw =  ',dtset%usepaw,' must be equal to 0 ',ch10,&
&   'Not functional for PAW case yet.'
   MSG_ERROR(message)
 end if

 cplex=2
 fft_option = 0 ! just do direct fft
 spin_symbol = (/'x','y','z'/)
 spinor_sym = (/'u','d'/)
 realimag = (/'Re','Im'/)

 write(std_out,*) ' psps%mpsang,psps%mpssoang ', psps%mpsang,psps%mpssoang

!======================= main code ================================
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!first get normal contribution to current, as psi tau dpsidr + dpsidr tau psi
!where tau are 1/2 the pauli matrices
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

!init plane wave coeff counter
 icg = 0
!init plane wave counter
 ikg = 0
!init occupation/band counter
 iocc = 1

!rspace point, cartesian direction, spin pol=x,y,z
 ABI_ALLOCATE(spincurrent,(dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),3,3))
 spincurrent = zero

 ABI_ALLOCATE(dummy_denpot,(cplex*dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6)))

 ABI_ALLOCATE(gbound,(2*dtset%mgfft+8,2))

!allocate (density_matrix(2,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),dtset%nspinor,&
!&                           dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),dtset%nspinor))
!density_matrix= zero
 ABI_ALLOCATE(density,(2,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),dtset%nspinor,dtset%nspinor))
 density= zero

 ABI_ALLOCATE(dpsidr,(2,dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),dtset%nspinor,3))
 ABI_ALLOCATE(psi_r,(2,dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),dtset%nspinor))

!loop over kpoints
 do ikpt=1,dtset%nkpt

!  number of plane waves for this kpt
   npw = hdr%npwarr(ikpt)

!  allocate arrays dep on number of pw
   ABI_ALLOCATE(kg_k,(3,npw))
   ABI_ALLOCATE(gpsi,(2,npw,dtset%nspinor,3))
   ABI_ALLOCATE(psi,(2,npw,dtset%nspinor))
   ABI_ALLOCATE(kgcart,(3,npw))

!  get cartesian coordinates of k+G vectors around this kpoint
   do ig=1,npw
     kgcart(:,ig) = matmul(gprimd(:,:),dtset%kpt(:,ikpt)+kg(:,ikg+ig))
     kg_k (:,ig) = kg(:,ikg+ig)
   end do

!  get gbound
   call sphereboundary(gbound,dtset%istwfk(ikpt),kg_k,dtset%mgfft,npw)

!  loop over bands
   do iband=1,dtset%nband(ikpt)

!    prefactor for sum over bands and kpoints
     prefact_nk = hdr%occ(iocc) * dtset%wtk(ikpt)

!    initialize this wf
     gpsi=zero
     psi=zero
     psi(:,1:npw,1) = cg(:,icg+1:icg+npw)

!    multiply psi by - i 2 pi G
     do ig=1,npw
       gpsi(1,ig,:,1) =  two_pi * kgcart(1,ig)*psi(2,ig,:)
       gpsi(2,ig,:,1) = -two_pi * kgcart(1,ig)*psi(1,ig,:)
       gpsi(1,ig,:,2) =  two_pi * kgcart(2,ig)*psi(2,ig,:)
       gpsi(2,ig,:,2) = -two_pi * kgcart(2,ig)*psi(1,ig,:)
       gpsi(1,ig,:,3) =  two_pi * kgcart(3,ig)*psi(2,ig,:)
       gpsi(2,ig,:,3) = -two_pi * kgcart(3,ig)*psi(1,ig,:)
     end do

!    loop over spinorial components
     do ispinor=1,dtset%nspinor
!      FT Gpsi_x to real space
       call fourwf(cplex,dummy_denpot,gpsi(:,:,ispinor,1),dummy_fofgout,&
&       dpsidr(:,:,:,:,ispinor,1),gbound,gbound,&
&       hdr%istwfk(ikpt),kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,npw,&
&       npw,dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),&
&       fft_option,0,one,one,use_gpu_cuda=dtset%use_gpu_cuda)

!      FT Gpsi_y to real space
       call fourwf(cplex,dummy_denpot,gpsi(:,:,ispinor,2),dummy_fofgout,&
&       dpsidr(:,:,:,:,ispinor,2),gbound,gbound,&
&       hdr%istwfk(ikpt),kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,npw,&
&       npw,dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),&
&       fft_option,0,one,one,use_gpu_cuda=dtset%use_gpu_cuda)

!      FT Gpsi_z to real space
       call fourwf(cplex,dummy_denpot,gpsi(:,:,ispinor,3),dummy_fofgout,&
&       dpsidr(:,:,:,:,ispinor,3),gbound,gbound,&
&       hdr%istwfk(ikpt),kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,npw,&
&       npw,dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),&
&       fft_option,0,one,one,use_gpu_cuda=dtset%use_gpu_cuda)

!      FT psi to real space
       call fourwf(cplex,dummy_denpot,psi(:,:,ispinor),dummy_fofgout,&
&       psi_r(:,:,:,:,ispinor),gbound,gbound,&
&       hdr%istwfk(ikpt),kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,npw,&
&       npw,dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),&
&       fft_option,0,one,one,use_gpu_cuda=dtset%use_gpu_cuda)

     end do ! ispinor

!    dpsidr now contains the full derivative of psi wrt space (gradient) in cartesian coordinates

!    get 3 pauli matrix contributions to the current: x,y,z, cart dir, spin dir
     do icartdir=1,3

!      x pauli spin matrix
!      sigma_x =  | 0   1 |
!      | 1   0 |
       spincurrent(:,:,:,icartdir,1) =  spincurrent(:,:,:,icartdir,1) + prefact_nk * &
!      Re(psi_r(up)^* dpsidr(down))
&       real(psi_r(1,:,:,:,1)*dpsidr(1,:,:,:,2,icartdir)  &
&       + psi_r(2,:,:,:,1)*dpsidr(2,:,:,:,2,icartdir)  &
!      Re(psi_r(down)^* dpsidr(up))
&       + psi_r(1,:,:,:,2)*dpsidr(1,:,:,:,1,icartdir)  &
&       + psi_r(2,:,:,:,2)*dpsidr(2,:,:,:,1,icartdir))

!      y pauli spin matrix
!      sigma_y =  | 0  -i |
!      | i   0 |
       spincurrent(:,:,:,icartdir,2) =  spincurrent(:,:,:,icartdir,2) + prefact_nk * &
!      Re(-i psi_r(up)^* dpsidr(down))
&       real(psi_r(1,:,:,:,1)*dpsidr(2,:,:,:,2,icartdir)  &
&       - psi_r(2,:,:,:,1)*dpsidr(1,:,:,:,2,icartdir)  &
!      Re(i psi_r(down)^* dpsidr(up))
&       - psi_r(1,:,:,:,2)*dpsidr(2,:,:,:,1,icartdir)  &
&       + psi_r(2,:,:,:,2)*dpsidr(1,:,:,:,1,icartdir))

!      z pauli spin matrix
!      sigma_z =  | 1   0 |
!      | 0  -1 |
       spincurrent(:,:,:,icartdir,3) =  spincurrent(:,:,:,icartdir,3) + prefact_nk * &
!      Re(psi_r(up)^* dpsidr(up))
&       real(psi_r(1,:,:,:,1)*dpsidr(1,:,:,:,1,icartdir)  &
&       - psi_r(2,:,:,:,1)*dpsidr(2,:,:,:,1,icartdir)  &
!      Re(-psi_r(down)^* dpsidr(down))
&       - psi_r(1,:,:,:,2)*dpsidr(1,:,:,:,2,icartdir)  &
&       + psi_r(2,:,:,:,2)*dpsidr(2,:,:,:,2,icartdir))
     end do ! end icartdir

!
!    accumulate non local density matrix in real space
!    NOTE: if we are only using the local part of the current, this becomes the
!    density spinor matrix! (much lighter to calculate) rho(r, sigma, sigmaprime)
!
     do ispinor=1,dtset%nspinor
       do i3=1,dtset%ngfft(3)
         do i2=1,dtset%ngfft(2)
           do i1=1,dtset%ngfft(1)
             irealsp = i1 + (i2-1)*dtset%ngfft(1) + (i3-1)*dtset%ngfft(2)*dtset%ngfft(1)

             do ispinorp=1,dtset%nspinor
               density(1,irealsp,ispinor,ispinorp) = &
&               density(1,irealsp,ispinor,ispinorp) + &
&               prefact_nk * (psi_r(1,i1,i2,i3,ispinor)*psi_r(1,i1,i2,i3,ispinorp)&
&               +  psi_r(2,i1,i2,i3,ispinor)*psi_r(2,i1,i2,i3,ispinorp))
               density(2,irealsp,ispinor,ispinorp) = &
&               density(2,irealsp,ispinor,ispinorp) + &
&               prefact_nk * (psi_r(1,i1,i2,i3,ispinor)*psi_r(2,i1,i2,i3,ispinorp)&
&               -  psi_r(2,i1,i2,i3,ispinor)*psi_r(1,i1,i2,i3,ispinorp))

!              do i3p=1,dtset%ngfft(3)
!              do i2p=1,dtset%ngfft(2)
!              do i1p=1,dtset%ngfft(1)
!              irealsp_p = i1p + (i2p-1)*dtset%ngfft(1) + (i3p-1)*dtset%ngfft(2)*dtset%ngfft(1)
!
!              NOTE : sign changes in second terms below because rho = psi*(r) psi(rprime)
!
!              density_matrix(1,irealsp,ispinor,irealsp_p,ispinorp) = &
!              &           density_matrix(1,irealsp,ispinor,irealsp_p,ispinorp) + &
!              &           prefact_nk * (psi_r(1,i1,i2,i3,ispinor)*psi_r(1,i1p,i2p,i3p,ispinorp)&
!              &           +  psi_r(2,i1,i2,i3,ispinor)*psi_r(2,i1p,i2p,i3p,ispinorp))
!              density_matrix(2,irealsp,ispinor,irealsp_p,ispinorp) = &
!              &           density_matrix(2,irealsp,ispinor,irealsp_p,ispinorp) + &
!              &           prefact_nk * (psi_r(1,i1,i2,i3,ispinor)*psi_r(2,i1p,i2p,i3p,ispinorp)&
!              &           -  psi_r(2,i1,i2,i3,ispinor)*psi_r(1,i1p,i2p,i3p,ispinorp))
!              end do
!              end do
!              end do ! end irealspprime

             end do !end ispinorp do

           end do
         end do
       end do ! end irealsp
     end do !end ispinor do

!    update pw counter
     icg=icg+npw
     iocc=iocc+1
   end do ! iband

   ikg=ikg+npw

!  deallocate arrays dep on npw for this kpoint
   ABI_DEALLOCATE(kg_k)
   ABI_DEALLOCATE(gpsi)
   ABI_DEALLOCATE(psi)
   ABI_DEALLOCATE(kgcart)

 end do ! ikpt

 ABI_DEALLOCATE(dpsidr)
 ABI_DEALLOCATE(psi_r)
 ABI_DEALLOCATE(dummy_denpot)
 ABI_DEALLOCATE(gbound)

!prefactor for contribution to spin current
!prefactor is 1/2 * 1/2 * 2 Re(.):
!1/2 from the formula for the current
!1/2 from the use of the normalized Pauli matrices
!2 from the complex conjugate part
!total = 1/2
 spincurrent = half * spincurrent

!make array of positions for all points on grid
 ABI_ALLOCATE(position_op,(3,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3)))
 do i3=1,dtset%ngfft(3)
   do i2=1,dtset%ngfft(2)
     do i1=1,dtset%ngfft(1)
       position_op(:,i1,i2,i3) = matmul(hdr%rprimd,(/i1-one,i2-one,i3-one/))&
&       /(/dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3)/)
     end do
   end do
 end do

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!add electric field term to current. Non local term in case of pseudopotential SO
!present theory is that it is equal to A(r,r') = (W_SO(r,r') + W_SO(r',r))
!For the strictly local part of the current, this becomes 2 W_SO(r,r)
!
!W_SO is the prefactor in the spinorbit part of the potential, such that it
!can be written V_SO = W_SO . p (momentum operator)
!decomposed from V_SO = v_SO(r,r') L.S = v_SO(r,r') (rxp).S = v_SO(r,r') (Sxr).p
!and ensuring symmetrization for the r operator wrt the two arguments of v_SO(r,r')
!Hence:
!W_SO(r,r) = v_SO(r,r) (Sxr)
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

!allocate (vso_realspace_nl(2,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),dtset%nspinor,&
!& dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),dtset%nspinor))

!call vso_realspace_nonlop(atindx,atindx1,dtfil,dtset,gmet,gprimd,hdr,kg,&
!& mpi_enreg,nattyp,ph1d,position_op,psps,rmet,ucvol,vso_realspace_nl,ylm,ylmgr)
!anticommutator of VSO with position operator
!--- not needed in local spin current case ---

 ABI_ALLOCATE(vso_realspace,(2,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),dtset%nspinor,dtset%nspinor,3))

 call vso_realspace_local(dtset,hdr,position_op,psps,vso_realspace)


!multiply by density (or density matrix for nonlocal case)
!and add to spin current



 ABI_DEALLOCATE(density)

 realrecip = 0 ! real space for xsf output
 ABI_ALLOCATE(xcart,(3,dtset%natom))
 call xred2xcart(dtset%natom,hdr%rprimd,xcart,hdr%xred)

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!output 3 components of current for each real space point
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
 do ispindir=1, 3
!  choose rescale_current such that the maximum current component printed out
!  is 1 percent of lattice distance
!  By default XCrysDen multiplies by 200 to get something comparable to a distance in real space.
   rescale_current = maxval(abs(spincurrent(:, :, :, :, ispindir)))
   if (abs(rescale_current) < tol8) then
     rescale_current = one
   else
     rescale_current = 0.001_dp * sqrt(max(sum(hdr%rprimd(:,1)**2), &
&     sum(hdr%rprimd(:,2)**2), sum(hdr%rprimd(:,3)**2)))
   end if

   filnam=trim(dtfil%fnameabo_spcur)//spin_symbol(ispindir)//".xsf"
   if (open_file(filnam,message,newunit=spcur_unit,status='unknown') /= 0) then
     MSG_ERROR(message)
   end if

!  print header
   write (spcur_unit,'(a)')          '#'
   write (spcur_unit,'(a)')          '#  Xcrysden format file'
   write (spcur_unit,'(a)')          '#  spin current density, for all real space points'
   write (spcur_unit,'(a,3(I5,1x))') '#  fft grid is ', dtset%ngfft(1), dtset%ngfft(2), dtset%ngfft(3)
   write (spcur_unit,'(a,a,a)')  '# ', spin_symbol(ispindir), '-spin current, full vector '

   write (spcur_unit,'(a)')  'ATOMS'
   do iatom = 1, dtset%natom
     write (spcur_unit,'(I4, 2x, 3(E16.6, 1x))') int(dtset%znucl(dtset%typat(iatom))), xcart(:,iatom)
   end do

   do i3=1,dtset%ngfft(3)
     do i2=1,dtset%ngfft(2)
       do i1=1,dtset%ngfft(1)
         write (spcur_unit,'(a, 3(E10.3),2x, 3(E20.10))') 'X ', &
&         position_op(:, i1, i2, i3), spincurrent(i1, i2, i3, :, ispindir)*rescale_current
       end do
     end do
   end do
   close (spcur_unit)

 end do ! end ispindir

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!output 3 spin components of V_SO matrices, for each real space point
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
 ABI_MALLOC(datagrid, (dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3)))

 do ispindir=1,3
   do icplex=1,2
     do ispinor=1,dtset%nspinor
       do ispinorp=1,dtset%nspinor

!        for the moment only print out if non zero
         if (abs(sum(vso_realspace(icplex, :, ispinor, ispinorp, ispindir))) < tol8) cycle

         filnam=trim(dtfil%fnameabo_vso)//"_spin_"//spin_symbol(ispindir)//"_"//&
&         spinor_sym(ispinor)//spinor_sym(ispinorp)//"_"//realimag(icplex)//".xsf"

         if (open_file(filnam,message,newunit=spcur_unit,status='unknown') /= 0) then
           MSG_ERROR(message)
         end if

         ! print header
         write (spcur_unit,'(a)')        '#'
         write (spcur_unit,'(a)')        '#  Xcrysden format file'
         write (spcur_unit,'(a)')        '#  spin-orbit potential (space diagonal), for all real space points'
         write (spcur_unit,'(a)')        '#    Real part first, then imaginary part'
         write (spcur_unit,'(a,3(I5,1x))') '#  fft grid is ', dtset%ngfft(1), dtset%ngfft(2),   dtset%ngfft(3)
         write (spcur_unit,'(a,a,a)')    '# ', spin_symbol(ispindir), '-spin contribution '
         write (spcur_unit,'(a,a,a)')    '# ', spinor_sym(ispinor)//spinor_sym(ispinorp), &
             '-spin element of the spinor 2x2 matrix '
         write (spcur_unit,'(a,a)')      '#  cart x     *  cart y    *  cart z    ***',&
&         ' up-up component    *  up-down           * down-up          * down-down          '

         ! Build contiguous array
         datagrid = vso_realspace(icplex,:,ispinor, ispinorp, ispindir)

         call printxsf(dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),&
&         datagrid,hdr%rprimd,(/zero,zero,zero/), dtset%natom, dtset%ntypat, &
&         dtset%typat, xcart, dtset%znucl, spcur_unit,realrecip)

!
!        NOTE: have chosen actual dims of grid (n123) instead of fft box, for which n45
!        may be different
!
!        do i3_dum=1,dtset%ngfft(3)+1
!        i3 = mod(i3_dum-1,dtset%ngfft(3)) + 1
!        do i2_dum=1,dtset%ngfft(2)+1
!        i2 = mod(i2_dum-1,dtset%ngfft(2)) + 1
!        do i1_dum=1,dtset%ngfft(1)+1
!        i1 = mod(i1_dum-1,dtset%ngfft(1)) + 1
!
!        irealsp = i1 + (i2-1)*dtset%ngfft(1) + (i3-1)*dtset%ngfft(2)*dtset%ngfft(1)
!        write (spcur_unit,'(E20.10,1x)')&
!        &      vso_realspace(icplex,irealsp,  ispinor, ispinorp, ispindir)
!        end do
!        end do
!        end do

         close (spcur_unit)

       end do ! ispinorp
     end do ! ispinor
   end do ! icplex
 end do ! end ispindir

 ABI_FREE(datagrid)
 ABI_DEALLOCATE(vso_realspace)
!deallocate (vso_realspace_nl)
 ABI_DEALLOCATE(position_op)
 ABI_DEALLOCATE(spincurrent)
 ABI_DEALLOCATE(xcart)

 write(std_out,*) ' Exiting subroutine spin_current '

end subroutine spin_current
!!***

!!****f* m_spin_current/vso_realspace_local
!! NAME
!!   vso_realspace_local
!!
!! FUNCTION
!!  Calculate real space (local - (r,r)) values of the SO part of the
!!  pseudopotential. Reconstructed explicitly in the HGH/GTH case.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      spin_current
!!
!! CHILDREN
!!      gamma_function,spline,splint,xred2xcart
!!
!! SOURCE

subroutine vso_realspace_local(dtset,hdr,position_op,psps,vso_realspace)

!Arguments -------------------------------
 type(hdr_type),intent(inout) :: hdr
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 real(dp),intent(in) :: position_op(3,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3))
 real(dp),intent(out) :: vso_realspace(2,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),&
& dtset%nspinor,dtset%nspinor,3)

!Local variables -------------------------
!scalars
 integer :: i,j,l, lmax,ipsp,iatom, ir1,ir2,ir3
 integer :: rcexponent,irealsp
 integer :: nradgrid,iradgrid
 real(dp) :: gammai, gammaj, relative_position(3), radial_cutoff, norm_rel_pos
 real(dp) :: expfact,lfact, vso_interpol, x,y,z
!arrays
 real(dp) :: xcart(3,dtset%natom),splint_x(1),splint_y(1)
 real(dp), allocatable :: radial_grid(:)
 real(dp), allocatable :: prefact_ijl(:,:,:,:),tmpvso(:),tmpvso_pp(:)
 real(dp), allocatable :: vso_radial(:,:),vso_radial_pp(:,:),tmp_spline(:)
 real(dp), allocatable :: offdiag_l_fact(:,:,:),kpar_matrix(:,:)

! *********************************************************************

!recalculate xcart (option = 1)
 call xred2xcart(dtset%natom,hdr%rprimd,xcart,hdr%xred)

 lmax = psps%mpsang-1

!content of gth pseudo type:
!These are {rloc, C(1...4)} coefficients for psppar(0, :, :) indices,
!Followed by the h coefficients for psppar(1:2, 1:, :) indices.
!size (0:2, 0:4, npsp)
!potential radius r_l is in psppar(l+1,0,ipsp)
!real(dp), pointer :: psppar(:, :, :)
!The covalence radii for each pseudo (?) size (npsp)
!real(dp), pointer :: radii_cov(:)
!Cut-off radii for core part and long-range part.
!radii_cf(:, 1) is for the long-range cut-off and
!radii_cf(:, 2) is for the core cut-off.
!size (npsp, 2)
!real(dp), pointer :: radii_cf(:, :)
!Spin orbit coefficients in HGH/GTH formats: k11p
!etc... see psp3ini.F90
!dimension = num l channels, 3 coeffs, num psp =
!(1:lmax+1,1:3,npsp)
!real(dp), pointer :: psp_k_par(:, :, :)

!v_SO^l (r,r') = sum_i sum_j sum_m Y_{lm} (\hat{r}) p_i^l (r) k_{ij}^l p_j^l(r') Y^{*}_lm (\hat{r'})
!
!v_SO^l (r,r)  = sum_ij  p_i^l (r) k_{ij}^l p_j^l(r) sum_m Y_{lm} (\hat{r}) Y^{*}_lm (\hat{r})
!= (2l+1)/4\pi sum_ij  p_i^l (r) k_{ij}^l p_j^l(r) (eq B.17 Patrick Rinke thesis)
!p are gaussian projectors (from HGH paper prb 58 3641) [[cite:Hartwigsen1998]]
!sum_l v_SO^l (r,r) is a purely radial quantity (function of |r|), so spline it

!maximum distance needed in unit cell
 radial_cutoff = four * maxval(psps%gth_params%psppar(:, 0, :))

!setup radial grid; Should we use a logarithmic grid? The spline functions can
!take it...
 nradgrid = 201 ! this is heuristic
 ABI_ALLOCATE(radial_grid,(nradgrid))
 do iradgrid=1,nradgrid
   radial_grid(iradgrid) = (iradgrid-1)*radial_cutoff/(nradgrid-1)
 end do

!calculate prefactors independent of r
 ABI_ALLOCATE(prefact_ijl,(3,3,0:lmax,psps%npsp))
 ABI_ALLOCATE(offdiag_l_fact,(3,3,0:lmax))
 ABI_ALLOCATE(kpar_matrix,(3,3))

!these factors complete the full 3x3 matrix of k (or h) parameters for the
!HGH pseudos
 offdiag_l_fact = zero
!l=0
 offdiag_l_fact(1,2,0) = -half*sqrt(three/five)
 offdiag_l_fact(1,3,0) = half*sqrt(five/21._dp)
 offdiag_l_fact(2,3,0) = -half*sqrt(100._dp/63._dp)
!l=1
 offdiag_l_fact(1,2,1) = -half*sqrt(five/seven)
 offdiag_l_fact(1,3,1) = sixth*sqrt(35._dp/11._dp)
 offdiag_l_fact(2,3,1) = -sixth*14._dp/sqrt(11._dp)
!l=2
 if (lmax >= 2) then
   offdiag_l_fact(1,2,2) = -half*sqrt(seven/nine)
   offdiag_l_fact(1,3,2) = half*sqrt(63._dp/143._dp)
   offdiag_l_fact(2,3,2) = -half*18._dp /sqrt(143._dp)
 end if
!l=3
 if (lmax >= 3) then
   offdiag_l_fact(1,2,3) = zero
   offdiag_l_fact(1,3,3) = zero
   offdiag_l_fact(2,3,3) = zero
 end if
!get prefactors for evaluation of V_SO: terms that do not depend on r
 prefact_ijl = zero
 do l=0,lmax
!  first the diagonal i=j term
   do i=1,3
     call gamma_function(l+(4._dp*i-1._dp)*0.5_dp, gammai)
     gammai = sqrt(gammai)
     rcexponent = 2*l+2*i+2*i-1
     do ipsp=1,psps%npsp
       prefact_ijl(i,i,l,ipsp) = psps%gth_params%psp_k_par(l+1,i,ipsp) &
&       / ( (psps%gth_params%psppar(l+1,0,ipsp))**(rcexponent) &
&       * gammai * gammai)
     end do
   end do
!  now the off diagonal elements
   do ipsp=1,psps%npsp
     kpar_matrix(1,2) = offdiag_l_fact (1,2,l)* psps%gth_params%psp_k_par(l+1,2,ipsp)
     kpar_matrix(2,1) = kpar_matrix(1,2)
     kpar_matrix(1,3) = offdiag_l_fact (1,3,l)* psps%gth_params%psp_k_par(l+1,3,ipsp)
     kpar_matrix(3,1) = kpar_matrix(1,3)
     kpar_matrix(2,3) = offdiag_l_fact (2,3,l)* psps%gth_params%psp_k_par(l+1,3,ipsp)
     kpar_matrix(3,2) = kpar_matrix(2,3)
   end do

!  for the f case only the 1,1 matrix element is non 0 - it is done above and
!  all these terms are actually 0
   if (l > 2) cycle

   do i=1,3
     call gamma_function(l+(4._dp*i-1._dp)*0.5_dp, gammai)
     gammai = sqrt(gammai)
     do j=1,3
       if (j==i) cycle
       rcexponent = 2*l+2*i+2*j-1
       call gamma_function(l+(4._dp*j-1._dp)*0.5_dp,gammaj)
       gammaj = sqrt(gammaj)
       do ipsp=1,psps%npsp
         prefact_ijl(i,j,l,ipsp) = kpar_matrix(i,j) &
&         / ( (psps%gth_params%psppar(l+1,0,ipsp))**rcexponent &
&         * gammai * gammaj )
       end do
     end do
   end do
 end do

 ABI_DEALLOCATE(kpar_matrix)
 ABI_DEALLOCATE(offdiag_l_fact)

 prefact_ijl = prefact_ijl * two

!calculate v_SO on radial grid
! MGNAG Runtime Error: *** Arithmetic exception: Floating invalid operation - aborting
 ABI_ALLOCATE(vso_radial,(nradgrid,psps%npsp))
 vso_radial = zero
 do l=0,lmax
   lfact=(2._dp*l+1._dp)/four/pi
   do iradgrid=1,nradgrid
     norm_rel_pos = radial_grid(iradgrid)
     do ipsp=1,psps%npsp
       expfact = exp(-norm_rel_pos**2 / (psps%gth_params%psppar(l+1,0,ipsp))**2)

       do i=1,3
         do j=1,3
           rcexponent = 2*l +2*i+2*j-4
           if(prefact_ijl(i,j,l,ipsp)/=0) then !vz_d 0**0
             vso_radial(iradgrid,ipsp) = vso_radial(iradgrid,ipsp) + &
&             prefact_ijl(i,j,l,ipsp)*(norm_rel_pos**rcexponent) * expfact
           end if  !vz_d
         end do ! j
       end do ! i
     end do ! ipsp
   end do ! iradgrid
 end do ! lmax

!spline v_SO(radial coord): get second derivative coefficients
 ABI_ALLOCATE(vso_radial_pp,(nradgrid,psps%npsp))

 ABI_ALLOCATE(tmp_spline,(nradgrid))
 ABI_ALLOCATE(tmpvso,(nradgrid))
 ABI_ALLOCATE(tmpvso_pp,(nradgrid))
 do ipsp=1,psps%npsp
   tmpvso = vso_radial(:,ipsp)
   call spline( radial_grid, tmpvso, nradgrid, zero, radial_grid(nradgrid), tmpvso_pp )
   vso_radial_pp(:,ipsp) = tmpvso_pp
 end do
 ABI_DEALLOCATE(tmp_spline)
 ABI_DEALLOCATE(tmpvso)
 ABI_DEALLOCATE(tmpvso_pp)

!to optimize this I should precalculate the distances which are actually needed by
!symmetry, or only sum over irreducible points in space and use weights

!for each physical atom present in unit cell
 vso_realspace = zero
 do iatom=1,dtset%natom
!  atom type will be dtset%typat(iatom)

!  for each point on grid
   do ir3=1,dtset%ngfft(3)
     do ir2=1,dtset%ngfft(2)
       do ir1=1,dtset%ngfft(1)
         irealsp = ir1 + (ir2-1)*dtset%ngfft(1) + (ir3-1)*dtset%ngfft(2)*dtset%ngfft(1)

!        relative position from atom to point
         relative_position = position_op(:,ir1,ir2,ir3) - xcart(:,iatom)
         x=relative_position(1)
         y=relative_position(2)
         z=relative_position(3)

!        calculate norm^2
         norm_rel_pos = relative_position(1)**2+relative_position(2)**2+relative_position(3)**2

!        if norm^2 is too large, skip this point
         if (norm_rel_pos > radial_cutoff*radial_cutoff) cycle

!        calculate norm
         splint_x(1) = sqrt(norm_rel_pos)

!        spline interpolate vso only depends on position (through pos - atomic position)
         call splint (nradgrid,radial_grid,vso_radial(:,dtset%typat(iatom)),&
&         vso_radial_pp(:,dtset%typat(iatom)),1,splint_x,splint_y)
         vso_interpol=splint_y(1)

!        multiply by vectorial spin factor (S x r)
!        NOTE: this r is taken relative to atom center. It could be that the r operator should
!        applied in an absolute way wrt the origin...
!
!        Is this correct: accumulated sum over atoms ?
         vso_realspace(1,irealsp,:,:,1) = vso_realspace(1,irealsp,:,:,1) + &
&         vso_interpol * reshape((/y,   zero,zero,-y/),(/2,2/))
         vso_realspace(2,irealsp,:,:,1) = vso_realspace(2,irealsp,:,:,1) + &
&         vso_interpol * reshape((/zero,z,  -z,    zero/),(/2,2/))

         vso_realspace(1,irealsp,:,:,2) = vso_realspace(1,irealsp,:,:,2) + &
&         vso_interpol * reshape((/-x,  z,   z,   x/),(/2,2/))
         vso_realspace(2,irealsp,:,:,2) = vso_realspace(2,irealsp,:,:,2) + &
&         vso_interpol * reshape((/zero,zero,zero,zero/),(/2,2/))

         vso_realspace(1,irealsp,:,:,3) = vso_realspace(1,irealsp,:,:,3) + &
&         vso_interpol * reshape((/zero,-y, -y,   zero/),(/2,2/))
         vso_realspace(2,irealsp,:,:,3) = vso_realspace(2,irealsp,:,:,3) + &
&         vso_interpol * reshape((/zero,-x, -x,    zero/),(/2,2/))

       end do  ! ir3
     end do  ! ir2
   end do  ! ir1
 end do ! iatom

 ABI_DEALLOCATE(prefact_ijl)
 ABI_DEALLOCATE(vso_radial)
 ABI_DEALLOCATE(vso_radial_pp)
 ABI_DEALLOCATE(radial_grid)

end subroutine vso_realspace_local
!!***

end module m_spin_current
!!***
