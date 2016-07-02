!{\src2tex{textfont=tt}}
!!****f* ABINIT/exc_den
!! NAME
!!  exc_den
!!
!! FUNCTION
!!  This routines calculates the electron-hole excited state density.
!!
!! COPYRIGHT
!! Copyright (C) 1992-2009 EXC group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida)
!! Copyright (C) 2009-2016 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  filbseig=Name of the file containing the excitonic eigenvectors and eigenvalues.
!!
!! OUTPUT
!!
!! PARENTS
!!      bethe_salpeter
!!
!! CHILDREN
!!      wfd_get_ur,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine exc_den(BSp,BS_files,ngfft,nfftot,Kmesh,ktabr,Wfd)

 use defs_basis
 use m_profiling_abi
 use m_bs_defs
 use m_errors

 use m_io_tools,        only : open_file
 use m_bz_mesh,         only : kmesh_t
 use m_wfd,             only : wfd_t, wfd_get_ur

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_den'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftot
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(kmesh_t),intent(in) :: Kmesh
 type(wfd_t),intent(inout) :: Wfd
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(in) :: ktabr(nfftot,BSp%nkbz)

!Local variables ------------------------------
!scalars
 integer :: nfft1,nfft2,nfft3,spin,ierr
 integer :: it,itp,ik_ibz,ikp_bz,ik_bz,band,iv,ivp,ic,icp,ir,i1,i2,i3,iik
 integer :: state,nstates,min_idx,eig_unt,den_unt,sden_unt,hsize,homo
 real(dp) :: min_ene,cost
 character(len=500) :: msg
 character(len=fnlen) :: filbseig
!arrays
 real(dp) :: n0(nfftot),rho_eh(nfftot),nexc(nfftot)
 real(dp),allocatable :: exc_ene(:)
 complex(dpc) :: rhor_h(nfftot),rhor_e(nfftot)
 complex(gwpc),allocatable :: wfr(:,:,:), wfrk(:,:)
 complex(dpc),allocatable :: exc_ene_cplx(:),exc_vec(:)

!************************************************************************

 ABI_CHECK(Wfd%nsppol==1,"nsppol==2 not coded")

 if (BS_files%in_eig /= BSE_NOFILE) then
   filbseig = BS_files%in_eig
 else 
   filbseig = BS_files%out_eig
 end if

 call wrtout(std_out,' Calculating electron-hole, excited state density',"COLL")
 call wrtout(std_out," Reading eigenstates from: "//TRIM(filbseig),"COLL")
 MSG_ERROR("Not tested")

 ! Prepare FFT tables to have u(r) on the ngfft_osc mesh.
 !if ( ANY(ngfftf(1:3) /= Wfd%ngfft(1:3)) ) then
 !  call wfd_change_ngfft(Wfd,Cryst,Psps,ngfftf) 
 !end if

 nfft1 = ngfft(1)
 nfft2 = ngfft(2)
 nfft3 = ngfft(3)
 ABI_CHECK(nfftot==PRODUCT(ngfft(1:3)),"Mismatch in FFT size")

!allocate and load wavefunctions in real space
 ABI_STAT_MALLOC(wfr,(nfftot*Wfd%nspinor,BSp%nbnds,Wfd%nkibz), ierr)
 ABI_CHECK(ierr==0, "out of memory: exc_den, wfr")

 spin=1 ! SPIN support is missing.

 do ik_ibz=1,Wfd%nkibz
   do band=1,BSp%nbnds
     call wfd_get_ur(Wfd,band,ik_ibz,spin,wfr(:,band,ik_ibz))
   end do
 end do

 if (open_file(filbseig,msg,newunit=eig_unt,form="unformatted", status="old", action="read") /= 0) then
   MSG_ERROR(msg)
 end if

 read(eig_unt) hsize,nstates

 if (nstates/=Bsp%nreh(1)) then 
   write(std_out,*) 'not resonant calculation'
   close(eig_unt)
   RETURN 
 end if

 ABI_MALLOC(exc_ene_cplx,(nstates))
 ABI_MALLOC(exc_ene,(nstates))

 read(eig_unt) exc_ene_cplx(1:nstates)

 ! Small imaginary part is always neglected.
 exc_ene(:) = exc_ene_cplx(:)
 ABI_FREE(exc_ene_cplx)
 !
 ! Find the lowest non-negative eigenvalue.
 min_ene = greatest_real 
 do state=1,nstates
   if (exc_ene(state) < zero) cycle
   if (exc_ene(state) < min_ene) then
     min_ene = exc_ene(state)
     min_idx = state
   end if
 end do
 ABI_FREE(exc_ene)

 write(std_out,*)' Lowest eigenvalue ', min_idx, min_ene*Ha_eV
 !
 ! skip other eigenvectors
 do state=1,min_idx-1
   read(eig_unt)
 end do
 !
 ! read "lowest" eigenvector
 ABI_MALLOC(exc_vec,(hsize))
 read(eig_unt) exc_vec

 close(eig_unt)
 
 ABI_STAT_MALLOC(wfrk,(nfftot,BSp%nbnds), ierr)
 ABI_CHECK(ierr==0, 'out of memory: exc_den, wfrk')

!calculate ground state density
 n0(:) = zero
 spin = 1
 homo = Bsp%homo_spin(spin)
 do ik_bz = 1, BSp%nkbz
   ik_ibz = Kmesh%tab(ik_bz)
   iik = (3-Kmesh%tabi(ik_bz))/2
    
   if (iik==1) then
     do ir=1,nfftot
       wfrk(ir,1:homo) = wfr(ktabr(ir,ik_bz),1:homo,ik_ibz)
     end do
   else
     do ir=1,nfftot
       wfrk(ir,1:homo) = conjg(wfr(ktabr(ir,ik_bz),1:homo,ik_ibz))
     end do
   end if

   do iv=1,homo
     n0(:) = n0(:) + conjg(wfrk(:,band)) * wfrk(:,band)
   end do
 end do
 !
 ! calculate electron and hole density
 rhor_h=czero
 rhor_e=czero
 ABI_CHECK(BSp%nsppol==1,"nsppol=2 not coded")

 do it=1,Bsp%nreh(1)
   ik_bz = Bsp%Trans(it,1)%k
   iv    = Bsp%Trans(it,1)%v
   ic    = Bsp%Trans(it,1)%c
   ik_ibz = Kmesh%tab(ik_bz)
   iik = (3-Kmesh%tabi(ik_bz))/2
    
   if (iik==1) then
     do ir = 1, nfftot
       wfrk(ir,:) = wfr(ktabr(ir,ik_bz),:,ik_ibz)
     end do
   else
     do ir = 1, nfftot
       wfrk(ir,:) = conjg(wfr(ktabr(ir,ik_bz),:,ik_ibz))
     end do
   end if
   !
   do itp=1,Bsp%nreh(1)
     ikp_bz = Bsp%Trans(itp,1)%k
     if (ikp_bz/=ik_bz) CYCLE
     icp = Bsp%Trans(itp,1)%c
     ivp = Bsp%Trans(itp,1)%v
     if(icp==ic) then
       rhor_h = rhor_h - conjg(exc_vec(it)) * exc_vec(itp) * wfrk(:,iv) * conjg(wfrk(:,ivp))
     end if
     if(ivp==iv) then
       rhor_e = rhor_e + conjg(exc_vec(it)) * exc_vec(itp) * conjg(wfrk(:,ic)) * wfrk(:,icp)
     end if
   end do
 end do

 ABI_FREE(exc_vec)
 !
 ! calculate excited state density minus ground state density
 ! n* - n0 = rhor_e + rhor_h
 rho_eh(:) = rhor_e + rhor_h
 !
 !calculate excited state density
 ! n* = n0 + rhor_e + rhor_h
 nexc(:) = n0(:) + rho_eh(:)

 if (open_file('out.den',msg,newunit=den_unt) /= 0) then
   MSG_ERROR(msg) 
 end if

 ! here conversion to cartesian through a1, a2, a3
 cost = zero
 do i1=0,nfft1-1
   do i2=0,nfft2-1
     do i3=0,nfft3-1
       ir=1 + i1 + i2*nfft1 + i3*nfft1*nfft2
       write(den_unt,'(3i3,2x,5e11.4)')i1,i2,i3,n0(ir),nexc(ir),rho_eh(ir),real(rhor_e(ir)),real(rhor_h(ir))
       cost = cost + n0(ir)
     end do
   end do
 end do

 write(std_out,*) 'density normalization constant ', cost
 close(den_unt)

 if (open_file('out.sden',msg,newunit=sden_unt) /= 0) then
   MSG_ERROR(msg) 
 end if

 ! we are looking for the plane between (100) (111)
 ! so it's the place where (111) [ (100) x v ] = 0 (mixed product)
 ! finally v2 - v3 = 0
 do i2=0, nfft2-1
   do i3=0, nfft3-1
     if(i2==i3) then
       write(std_out,*) i2, i3
       write(sden_unt,*) (rho_eh(1+i1+i2*nfft1+i3*nfft1*nfft2),i1=0,nfft1-1)
     end if
   end do
 end do

 close(sden_unt)

end subroutine exc_den      
!!***
