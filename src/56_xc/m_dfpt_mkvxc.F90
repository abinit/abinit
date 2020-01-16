!!****m* ABINIT/m_dfpt_mkvxc
!! NAME
!!  m_dfpt_mkvxc
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2019 ABINIT group (XG, DRH, FR, EB, SPr)
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

module m_dfpt_mkvxc

 use defs_basis
 use m_errors
 use m_abicore
 use m_xc_noncoll

 use defs_abitypes,     only : MPI_type
 use m_time,     only : timab
 use m_symtk,    only : matr3inv
 use m_xctk,     only : xcden, xcpot

 implicit none

 private
!!***

 public :: dfpt_mkvxc
 public :: dfpt_mkvxc_noncoll
!!***

contains
!!***

!!****f* ABINIT/dfpt_mkvxc
!! NAME
!! dfpt_mkvxc
!!
!! FUNCTION
!! Compute the first-order change of exchange-correlation potential
!! due to atomic displacement: assemble the first-order
!! density change with the frozen-core density change, then use
!! the exchange-correlation kernel.
!!
!! INPUTS
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!         if 2, COMPLEX
!!  ixc= choice of exchange-correlation scheme
!!  kxc(nfft,nkxc)=exchange and correlation kernel (see below)
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/variables/vargs.htm#ngfft
!!  nhat1(cplex*nfft,2nspden*nhat1dim)= -PAW only- 1st-order compensation density
!!  nhat1dim= -PAW only- 1 if nhat1 array is used ; 0 otherwise
!!  nhat1gr(cplex*nfft,nspden,3*nhat1grdim)= -PAW only- gradients of 1st-order compensation density
!!  nhat1grdim= -PAW only- 1 if nhat1gr array is used ; 0 otherwise
!!  nkxc=second dimension of the kxc array
!!  non_magnetic_xc= if true, handle density/potential as non-magnetic (even if it is)
!!  nspden=number of spin-density components
!!  n3xccc=dimension of xccc3d1 ; 0 if no XC core correction is used, otherwise, nfft
!!  option=if 0, work only with the XC core-correction,
!!         if 1, treat both density change and XC core correction
!!         if 2, treat only density change
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rhor1(cplex*nfft,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  usexcnhat= -PAW only- 1 if nhat density has to be taken into account in Vxc
!!  xccc3d1(cplex*n3xccc)=3D change in core charge density, see n3xccc
!!
!! OUTPUT
!!  vxc1(cplex*nfft,nspden)=change in exchange-correlation potential (including
!!   core-correction, if applicable)
!!
!! NOTES
!!  Content of Kxc array:
!!   ===== if LDA
!!    if nspden==1: kxc(:,1)= d2Exc/drho2
!!                 (kxc(:,2)= d2Exc/drho_up drho_dn)
!!    if nspden>=2: kxc(:,1)= d2Exc/drho_up drho_up
!!                  kxc(:,2)= d2Exc/drho_up drho_dn
!!                  kxc(:,3)= d2Exc/drho_dn drho_dn
!!   ===== if GGA
!!    if nspden==1:
!!       kxc(:,1)= d2Exc/drho2
!!       kxc(:,2)= 1/|grad(rho)| dExc/d|grad(rho)|
!!       kxc(:,3)= 1/|grad(rho)| d2Exc/d|grad(rho)| drho
!!       kxc(:,4)= 1/|grad(rho)| * d/d|grad(rho)| ( 1/|grad(rho)| dExc/d|grad(rho)| )
!!       kxc(:,5)= gradx(rho)
!!       kxc(:,6)= grady(rho)
!!       kxc(:,7)= gradz(rho)
!!    if nspden>=2:
!!       kxc(:,1)= d2Exc/drho_up drho_up
!!       kxc(:,2)= d2Exc/drho_up drho_dn
!!       kxc(:,3)= d2Exc/drho_dn drho_dn
!!       kxc(:,4)= 1/|grad(rho_up)| dEx/d|grad(rho_up)|
!!       kxc(:,5)= 1/|grad(rho_dn)| dEx/d|grad(rho_dn)|
!!       kxc(:,6)= 1/|grad(rho_up)| d2Ex/d|grad(rho_up)| drho_up
!!       kxc(:,7)= 1/|grad(rho_dn)| d2Ex/d|grad(rho_dn)| drho_dn
!!       kxc(:,8)= 1/|grad(rho_up)| * d/d|grad(rho_up)| ( 1/|grad(rho_up)| dEx/d|grad(rho_up)| )
!!       kxc(:,9)= 1/|grad(rho_dn)| * d/d|grad(rho_dn)| ( 1/|grad(rho_dn)| dEx/d|grad(rho_dn)| )
!!       kxc(:,10)=1/|grad(rho)| dEc/d|grad(rho)|
!!       kxc(:,11)=1/|grad(rho)| d2Ec/d|grad(rho)| drho_up
!!       kxc(:,12)=1/|grad(rho)| d2Ec/d|grad(rho)| drho_dn
!!       kxc(:,13)=1/|grad(rho)| * d/d|grad(rho)| ( 1/|grad(rho)| dEc/d|grad(rho)| )
!!       kxc(:,14)=gradx(rho_up)
!!       kxc(:,15)=gradx(rho_dn)
!!       kxc(:,16)=grady(rho_up)
!!       kxc(:,17)=grady(rho_dn)
!!       kxc(:,18)=gradz(rho_up)
!!       kxc(:,19)=gradz(rho_dn)
!!
!! PARENTS
!!      dfpt_dyxc1,dfpt_mkvxc_noncoll,dfpt_nstdy,dfpt_nstpaw,dfpt_rhotov
!!      dfptnl_loop,m_kxc,nres2vres
!!
!! CHILDREN
!!      dfpt_mkvxcgga,matr3inv,timab
!!
!! SOURCE

subroutine dfpt_mkvxc(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,nhat1,nhat1dim,nhat1gr,nhat1grdim,&
&          nkxc,non_magnetic_xc,nspden,n3xccc,option,qphon,rhor1,rprimd,usexcnhat,vxc1,xccc3d1)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,ixc,n3xccc,nfft,nhat1dim,nhat1grdim
 integer,intent(in) :: nkxc,nspden,option,usexcnhat
 logical,intent(in) :: non_magnetic_xc
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in),target :: nhat1(cplex*nfft,nspden*nhat1dim)
 real(dp),intent(in),target :: nhat1gr(cplex*nfft,nspden,3*nhat1grdim)
 real(dp),intent(in) :: kxc(nfft,nkxc),qphon(3)
 real(dp),intent(in),target :: rhor1(cplex*nfft,nspden)
 real(dp),intent(in) :: rprimd(3,3),xccc3d1(cplex*n3xccc)
 real(dp),intent(out) :: vxc1(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: ii,ir,ispden,nhat1dim_,nhat1rgdim_
 real(dp) :: rho1_dn,rho1_up,rho1im_dn,rho1im_up,rho1re_dn,rho1re_up
 real(dp) :: spin_scale
!arrays
 real(dp) :: gprimd(3,3),tsec(2)
 real(dp), ABI_CONTIGUOUS pointer :: nhat1_(:,:),nhat1gr_(:,:,:),rhor1_(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(181,1,tsec)

 if(nspden/=1 .and. nspden/=2) then
   MSG_BUG('For nspden==4 please use dfpt_mkvxc_noncoll!')
 end if

!Special case: no XC applied
 if (ixc==0.or.nkxc==0) then
   MSG_WARNING('Note that no xc is applied (ixc=0)')
   vxc1=zero
   return
 end if

!Treat first LDA
 if(nkxc==1.or.nkxc==3)then

!  PAW: eventually substract compensation density
   if (option/=0) then
     if ((usexcnhat==0.and.nhat1dim==1).or.(non_magnetic_xc)) then
       ABI_ALLOCATE(rhor1_,(cplex*nfft,nspden))
       if (usexcnhat==0.and.nhat1dim==1) then
         rhor1_(:,:)=rhor1(:,:)-nhat1(:,:)
       else
         rhor1_(:,:)=rhor1(:,:)
       end if
       if (non_magnetic_xc) then
         if(nspden==2) rhor1_(:,2)=rhor1_(:,1)*half
         if(nspden==4) rhor1_(:,2:4)=zero
       end if
     else
       rhor1_ => rhor1
     end if
   end if

!  Case without non-linear core correction
   if(n3xccc==0 .or. option==2)then

     if(option==0)then  ! No straight XC to compute

       vxc1(:,:)=zero

     else               ! XC, without non-linear XC correction

!      Non-spin-polarized
       if(nspden==1)then
         if(cplex==1)then
           do ir=1,nfft
             vxc1(ir,1)=kxc(ir,1)*rhor1_(ir,1)
           end do
         else
           do ir=1,nfft
             vxc1(2*ir-1,1)=kxc(ir,1)*rhor1_(2*ir-1,1)
             vxc1(2*ir  ,1)=kxc(ir,1)*rhor1_(2*ir  ,1)
           end do
         end if ! cplex==1

!        Spin-polarized
       else
         if(cplex==1)then
           do ir=1,nfft
             rho1_dn=rhor1_(ir,1)-rhor1_(ir,2)
             vxc1(ir,1)=kxc(ir,1)*rhor1_(ir,2)+kxc(ir,2)*rho1_dn
             vxc1(ir,2)=kxc(ir,2)*rhor1_(ir,2)+kxc(ir,3)*rho1_dn
           end do
         else
           do ir=1,nfft
             rho1re_dn=rhor1_(2*ir-1,1)-rhor1_(2*ir-1,2)
             rho1im_dn=rhor1_(2*ir  ,1)-rhor1_(2*ir  ,2)
             vxc1(2*ir-1,1)=kxc(ir,1)*rhor1_(2*ir-1,2)+kxc(ir,2)*rho1re_dn
             vxc1(2*ir  ,1)=kxc(ir,1)*rhor1_(2*ir  ,2)+kxc(ir,2)*rho1im_dn
             vxc1(2*ir-1,2)=kxc(ir,2)*rhor1_(2*ir-1,2)+kxc(ir,3)*rho1re_dn
             vxc1(2*ir  ,2)=kxc(ir,2)*rhor1_(2*ir  ,2)+kxc(ir,3)*rho1im_dn
           end do
         end if ! cplex==1
       end if ! nspden==1

     end if ! option==0

!    Treat case with non-linear core correction
   else

     if(option==0)then

       if(nspden==1)then
         if(cplex==1)then
           do ir=1,nfft
             vxc1(ir,1)=kxc(ir,1)*xccc3d1(ir)
           end do
         else
           do ir=1,nfft
             vxc1(2*ir-1,1)=kxc(ir,1)*xccc3d1(2*ir-1)
             vxc1(2*ir  ,1)=kxc(ir,1)*xccc3d1(2*ir  )
           end do
         end if ! cplex==1
       else
         if(cplex==1)then
           do ir=1,nfft
             vxc1(ir,1)=(kxc(ir,1)+kxc(ir,2))*xccc3d1(ir)*half
             vxc1(ir,2)=(kxc(ir,2)+kxc(ir,3))*xccc3d1(ir)*half
           end do
         else
           do ir=1,nfft
             vxc1(2*ir-1,1)=(kxc(ir,1)+kxc(ir,2))*xccc3d1(2*ir-1)*half
             vxc1(2*ir  ,1)=(kxc(ir,1)+kxc(ir,2))*xccc3d1(2*ir  )*half
             vxc1(2*ir-1,2)=(kxc(ir,2)+kxc(ir,3))*xccc3d1(2*ir-1)*half
             vxc1(2*ir  ,2)=(kxc(ir,2)+kxc(ir,3))*xccc3d1(2*ir  )*half
           end do
         end if ! cplex==1
       end if ! nspden==1

     else ! option/=0

       if(nspden==1)then
         if(cplex==1)then
           do ir=1,nfft
             vxc1(ir,1)=kxc(ir,1)*(rhor1_(ir,1)+xccc3d1(ir))
           end do
         else
           do ir=1,nfft
             vxc1(2*ir-1,1)=kxc(ir,1)*(rhor1_(2*ir-1,1)+xccc3d1(2*ir-1))
             vxc1(2*ir  ,1)=kxc(ir,1)*(rhor1_(2*ir  ,1)+xccc3d1(2*ir  ))
           end do
         end if ! cplex==1
       else
         if(cplex==1)then
           do ir=1,nfft
             rho1_dn=rhor1_(ir,1)-rhor1_(ir,2) + xccc3d1(ir)*half
             rho1_up=rhor1_(ir,2)             + xccc3d1(ir)*half
             vxc1(ir,1)=kxc(ir,1)*rho1_up+kxc(ir,2)*rho1_dn
             vxc1(ir,2)=kxc(ir,2)*rho1_up+kxc(ir,3)*rho1_dn
           end do
         else
           do ir=1,nfft
             rho1re_dn=rhor1_(2*ir-1,1)-rhor1_(2*ir-1,2) + xccc3d1(2*ir-1)*half
             rho1im_dn=rhor1_(2*ir  ,1)-rhor1_(2*ir  ,2) + xccc3d1(2*ir  )*half
             rho1re_up=rhor1_(2*ir-1,2)                 + xccc3d1(2*ir-1)*half
             rho1im_up=rhor1_(2*ir  ,2)                 + xccc3d1(2*ir  )*half
             vxc1(2*ir-1,1)=kxc(ir,1)*rho1re_up+kxc(ir,2)*rho1re_dn
             vxc1(2*ir  ,1)=kxc(ir,1)*rho1im_up+kxc(ir,2)*rho1im_dn
             vxc1(2*ir-1,2)=kxc(ir,2)*rho1re_up+kxc(ir,3)*rho1re_dn
             vxc1(2*ir  ,2)=kxc(ir,2)*rho1im_up+kxc(ir,3)*rho1im_dn
           end do
         end if ! cplex==1
       end if ! nspden==1

     end if ! option==0

   end if ! n3xccc==0

   if (option/=0.and.((usexcnhat==0.and.nhat1dim==1).or.(non_magnetic_xc))) then
     ABI_DEALLOCATE(rhor1_)
   end if

!  Treat GGA
 else if (nkxc==7.or.nkxc==19) then

!  Transfer the data to spin-polarized storage

!  Treat the density change
   ABI_ALLOCATE(rhor1_,(cplex*nfft,nspden))
   if (option==1 .or. option==2) then
     if (nspden==1) then
       do ir=1,cplex*nfft
         rhor1_(ir,1)=rhor1(ir,1)
       end do
     else
       if(non_magnetic_xc) then
         do ir=1,cplex*nfft
           rho1_dn=rhor1(ir,1)*half
           rhor1_(ir,1)=rho1_dn
           rhor1_(ir,2)=rho1_dn
         end do
       else
         do ir=1,cplex*nfft
           rho1_dn=rhor1(ir,1)-rhor1(ir,2)
           rhor1_(ir,1)=rhor1(ir,2)
           rhor1_(ir,2)=rho1_dn
         end do
       end if
     end if
   else
     do ispden=1,nspden
       do ir=1,cplex*nfft
         rhor1_(ir,ispden)=zero
       end do
     end do
   end if

   if( (option==0 .or. option==1) .and. n3xccc/=0)then
     spin_scale=one;if (nspden==2) spin_scale=half
     do ispden=1,nspden
       do ir=1,cplex*nfft
         rhor1_(ir,ispden)=rhor1_(ir,ispden)+xccc3d1(ir)*spin_scale
       end do
     end do
   end if

!  PAW: treat also compensation density (and gradients)
   nhat1dim_=nhat1dim ; nhat1rgdim_=nhat1grdim
   if (option/=0.and.nhat1dim==1.and.nspden==2) then
     ABI_ALLOCATE(nhat1_,(cplex*nfft,nspden))
     if (non_magnetic_xc) then
       do ir=1,cplex*nfft
         rho1_dn=nhat1(ir,1)*half
         nhat1_(ir,1:2)=rho1_dn
       end do
     else
       do ir=1,cplex*nfft
         rho1_dn=nhat1(ir,1)-nhat1(ir,2)
         nhat1_(ir,1)=nhat1(ir,2)
         nhat1_(ir,2)=rho1_dn
       end do
     end if
   else if (option==0) then
     ABI_ALLOCATE(nhat1_,(0,0))
     nhat1dim_=0
   else
     nhat1_ => nhat1
   end if
   if (option/=0.and.nhat1grdim==1.and.nspden==2) then
     ABI_ALLOCATE(nhat1gr_,(cplex*nfft,nspden,3))
     if (non_magnetic_xc) then
       do ii=1,3
         do ir=1,cplex*nfft
           rho1_dn=nhat1(ir,1)*half
           nhat1gr_(ir,1:2,ii)=rho1_dn
         end do
       end do
     else
       do ii=1,3
         do ir=1,cplex*nfft
           rho1_dn=nhat1gr(ir,1,ii)-nhat1gr(ir,2,ii)
           nhat1gr_(ir,1,ii)=nhat1gr(ir,2,ii)
           nhat1gr_(ir,2,ii)=rho1_dn
         end do
       end do
     end if
   else if (option==0) then
     ABI_ALLOCATE(nhat1gr_,(0,0,0))
     nhat1rgdim_=0
   else
     nhat1gr_ => nhat1gr
   end if

   call matr3inv(rprimd,gprimd)

   call dfpt_mkvxcgga(cplex,gprimd,kxc,mpi_enreg,nfft,ngfft,nhat1_,nhat1dim_,&
&   nhat1gr_,nhat1rgdim_,nkxc,nspden,qphon,rhor1_,usexcnhat,vxc1)

   ABI_DEALLOCATE(rhor1_)
   if ((option==0).or.(nhat1dim==1.and.nspden==2)) then
     ABI_DEALLOCATE(nhat1_)
   end if
   if ((option==0).or.(nhat1grdim==1.and.nspden==2)) then
     ABI_DEALLOCATE(nhat1gr_)
   end if

 else
   MSG_BUG('Invalid nkxc!')

 end if ! LDA or GGA

 call timab(181,2,tsec)

 DBG_EXIT("COLL")

end subroutine dfpt_mkvxc
!!***

!!****f* ABINIT/dfpt_mkvxcgga
!! NAME
!! dfpt_mkvxcgga
!!
!! FUNCTION
!! Compute the first-order change of exchange-correlation potential
!! in case of GGA functionals
!! Use the exchange-correlation kernel.
!!
!! INPUTS
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!    if 2, COMPLEX
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  gprimd(3,3)=dimensional primitive translations in reciprocal space (bohr^-1)
!!  gsqcut=cutoff value on G**2 for sphere inside fft box.
!!  kxc(nfft,nkxc)=exchange and correlation kernel (see below)
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT
!!  nhat1(cplex*nfft,2*nhat1dim)= -PAW only- 1st-order compensation density
!!  nhat1dim= -PAW only- 1 if nhat1 array is used ; 0 otherwise
!!  nhat1gr(cplex*nfft,2,3*nhat1grdim)= -PAW only- gradients of 1st-order compensation density
!!  nhat1grdim= -PAW only- 1 if nhat1gr array is used ; 0 otherwise
!!  nkxc=second dimension of the kxc array
!!  nspden=number of spin-density components
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rhor1tmp(cplex*nfft,2)=array for first-order electron spin-density
!!   in electrons/bohr**3 (second index corresponds to spin-up and spin-down)
!!  usexcnhat= -PAW only- 1 if nhat density has to be taken into account in Vxc
!!
!! OUTPUT
!!  vxc1(cplex*nfft,nspden)=change in exchange-correlation potential
!!
!! NOTES
!!  For the time being, a rather crude coding, to be optimized ...
!!  Content of Kxc array:
!!   ===== if GGA
!!    if nspden==1:
!!       kxc(:,1)= d2Exc/drho2
!!       kxc(:,2)= 1/|grad(rho)| dExc/d|grad(rho)|
!!       kxc(:,3)= 1/|grad(rho)| d2Exc/d|grad(rho)| drho
!!       kxc(:,4)= 1/|grad(rho)| * d/d|grad(rho)| ( 1/|grad(rho)| dExc/d|grad(rho)| )
!!       kxc(:,5)= gradx(rho)
!!       kxc(:,6)= grady(rho)
!!       kxc(:,7)= gradz(rho)
!!    if nspden>=2:
!!       kxc(:,1)= d2Exc/drho_up drho_up
!!       kxc(:,2)= d2Exc/drho_up drho_dn
!!       kxc(:,3)= d2Exc/drho_dn drho_dn
!!       kxc(:,4)= 1/|grad(rho_up)| dEx/d|grad(rho_up)|
!!       kxc(:,5)= 1/|grad(rho_dn)| dEx/d|grad(rho_dn)|
!!       kxc(:,6)= 1/|grad(rho_up)| d2Ex/d|grad(rho_up)| drho_up
!!       kxc(:,7)= 1/|grad(rho_dn)| d2Ex/d|grad(rho_dn)| drho_dn
!!       kxc(:,8)= 1/|grad(rho_up)| * d/d|grad(rho_up)| ( 1/|grad(rho_up)| dEx/d|grad(rho_up)| )
!!       kxc(:,9)= 1/|grad(rho_dn)| * d/d|grad(rho_dn)| ( 1/|grad(rho_dn)| dEx/d|grad(rho_dn)| )
!!       kxc(:,10)=1/|grad(rho)| dEc/d|grad(rho)|
!!       kxc(:,11)=1/|grad(rho)| d2Ec/d|grad(rho)| drho_up
!!       kxc(:,12)=1/|grad(rho)| d2Ec/d|grad(rho)| drho_dn
!!       kxc(:,13)=1/|grad(rho)| * d/d|grad(rho)| ( 1/|grad(rho)| dEc/d|grad(rho)| )
!!       kxc(:,14)=gradx(rho_up)
!!       kxc(:,15)=gradx(rho_dn)
!!       kxc(:,16)=grady(rho_up)
!!       kxc(:,17)=grady(rho_dn)
!!       kxc(:,18)=gradz(rho_up)
!!       kxc(:,19)=gradz(rho_dn)
!!
!! PARENTS
!!      dfpt_mkvxc
!!
!! CHILDREN
!!      xcden,xcpot
!!
!! SOURCE

subroutine dfpt_mkvxcgga(cplex,gprimd,kxc,mpi_enreg,nfft,ngfft,&
&                    nhat1,nhat1dim,nhat1gr,nhat1grdim,nkxc,&
&                    nspden,qphon,rhor1,usexcnhat,vxc1)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nhat1dim,nhat1grdim,nkxc,nspden,usexcnhat
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: kxc(nfft,nkxc)
 real(dp),intent(in) :: nhat1(cplex*nfft,nspden*nhat1dim)
 real(dp),intent(in) :: nhat1gr(cplex*nfft,nspden,3*nhat1grdim)
 real(dp),intent(in) :: qphon(3)
 real(dp),intent(in),target :: rhor1(cplex*nfft,nspden)
 real(dp),intent(out) :: vxc1(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: ii,ir,ishift,ngrad,nspgrad,use_laplacian
 logical :: test_nhat
 real(dp) :: coeff_grho,coeff_grho_corr,coeff_grho_dn,coeff_grho_up
 real(dp) :: coeffim_grho,coeffim_grho_corr,coeffim_grho_dn,coeffim_grho_up
 real(dp) :: gradrho_gradrho1,gradrho_gradrho1_dn,gradrho_gradrho1_up
 real(dp) :: gradrho_gradrho1im,gradrho_gradrho1im_dn,gradrho_gradrho1im_up
 character(len=500) :: msg
!arrays
 real(dp) :: r0(3),r0_dn(3),r0_up(3),r1(3),r1_dn(3),r1_up(3)
 real(dp) :: r1im(3),r1im_dn(3),r1im_up(3)
 real(dp),allocatable :: dnexcdn(:,:),rho1now(:,:,:)
 real(dp),ABI_CONTIGUOUS pointer :: rhor1_ptr(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 if (nkxc/=12*min(nspden,2)-5) then
   msg='Wrong nkxc value for GGA!'
   MSG_BUG(msg)
 end if

!metaGGA contributions are not taken into account here
 use_laplacian=0

!PAW: substract 1st-order compensation density from 1st-order density
 test_nhat=((nhat1dim==1).and.(usexcnhat==0.or.nhat1grdim==1))
 if (test_nhat) then
   ABI_ALLOCATE(rhor1_ptr,(cplex*nfft,nspden))
   rhor1_ptr(:,:)=rhor1(:,:)-nhat1(:,:)
 else
   rhor1_ptr => rhor1
 end if

!call filterpot(paral_kgb,cplex,gmet,gsqcut,nfft,ngfft,2,qphon,rhor1_ptr)

!Compute the gradients of the first-order density
!rho1now(:,:,1) contains the first-order density, and
!rho1now(:,:,2:4) contains the gradients of the first-order density
 ishift=0 ; ngrad=2
 ABI_ALLOCATE(rho1now,(cplex*nfft,nspden,ngrad*ngrad))
 call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden,qphon,rhor1_ptr,rho1now)

!PAW: add "exact" gradients of compensation density
 if (test_nhat.and.usexcnhat==1) then
   rho1now(:,1:nspden,1)=rho1now(:,1:nspden,1)+nhat1(:,1:nspden)
 end if
 if (nhat1grdim==1) then
   do ii=1,3
     rho1now(:,1:nspden,ii+1)=rho1now(:,1:nspden,ii+1)+nhat1gr(:,1:nspden,ii)
   end do
 end if
 if (test_nhat) then
   ABI_DEALLOCATE(rhor1_ptr)
 end if

!Apply the XC kernel
 nspgrad=2; if (nspden==2) nspgrad=5
 ABI_ALLOCATE(dnexcdn,(cplex*nfft,nspgrad))

 if (cplex==1) then  ! Treat real case first
   if (nspden==1) then
     do ir=1,nfft
       r0(:)=kxc(ir,5:7) ; r1(:)=rho1now(ir,1,2:4)
       gradrho_gradrho1=dot_product(r0,r1)
       dnexcdn(ir,1)=kxc(ir,1)*rho1now(ir,1,1) + kxc(ir,3)*gradrho_gradrho1
       coeff_grho=kxc(ir,3)*rho1now(ir,1,1) + kxc(ir,4)*gradrho_gradrho1
  !    Reuse the storage in rho1now
       rho1now(ir,1,2:4)=r1(:)*kxc(ir,2)+r0(:)*coeff_grho
     end do
   else
     do ir=1,nfft
       do ii=1,3  ! grad of spin-up ans spin_dn GS rho
         r0_up(ii)=kxc(ir,13+2*ii);r0_dn(ii)=kxc(ir,12+2*ii)-kxc(ir,13+2*ii)
       end do
       r0(:)=r0_up(:)+r0_dn(:)      ! grad of GS rho
       r1_up(:)=rho1now(ir,1,2:4)   ! grad of spin-up rho1
       r1_dn(:)=rho1now(ir,2,2:4)   ! grad of spin-down rho1
       r1(:)=r1_up(:)+r1_dn(:)      ! grad of GS rho1
       gradrho_gradrho1_up=dot_product(r0_up,r1_up)
       gradrho_gradrho1_dn=dot_product(r0_dn,r1_dn)
       gradrho_gradrho1   =dot_product(r0,r1)
       dnexcdn(ir,1)=kxc(ir, 1)*rho1now(ir,1,1)     &
&       +kxc(ir, 2)*rho1now(ir,2,1)     &
&       +kxc(ir, 6)*gradrho_gradrho1_up &
&       +kxc(ir,11)*gradrho_gradrho1
       dnexcdn(ir,2)=kxc(ir, 3)*rho1now(ir,2,1)     &
&       +kxc(ir, 2)*rho1now(ir,1,1)     &
&       +kxc(ir, 7)*gradrho_gradrho1_dn &
&       +kxc(ir,12)*gradrho_gradrho1
       coeff_grho_corr=kxc(ir,11)*rho1now(ir,1,1) &
&       +kxc(ir,12)*rho1now(ir,2,1) &
&       +kxc(ir,13)*gradrho_gradrho1
       coeff_grho_up=kxc(ir,6)*rho1now(ir,1,1)+kxc(ir,8)*gradrho_gradrho1_up
       coeff_grho_dn=kxc(ir,7)*rho1now(ir,2,1)+kxc(ir,9)*gradrho_gradrho1_dn
  !    Reuse the storage in rho1now
       rho1now(ir,1,2:4)=(kxc(ir,4)+kxc(ir,10))*r1_up(:) &
&       +kxc(ir,10)            *r1_dn(:) &
&       +coeff_grho_up         *r0_up(:) &
&       +coeff_grho_corr       *r0(:)
       rho1now(ir,2,2:4)=(kxc(ir,5)+kxc(ir,10))*r1_dn(:) &
&       +kxc(ir,10)            *r1_up(:) &
&       +coeff_grho_dn         *r0_dn(:) &
&       +coeff_grho_corr       *r0(:)
     end do
   end if ! nspden

 else ! if cplex==2
   if (nspden==1) then
     do ir=1,nfft
       r0(:)=kxc(ir,5:7)
       r1(:)  =rho1now(2*ir-1,1,2:4)
       r1im(:)=rho1now(2*ir  ,1,2:4)
       gradrho_gradrho1  =dot_product(r0,r1)
       gradrho_gradrho1im=dot_product(r0,r1im)
       dnexcdn(2*ir-1,1)=kxc(ir,1)*rho1now(2*ir-1,1,1) + kxc(ir,3)*gradrho_gradrho1
       dnexcdn(2*ir  ,1)=kxc(ir,1)*rho1now(2*ir  ,1,1) + kxc(ir,3)*gradrho_gradrho1im
       coeff_grho  =kxc(ir,3)*rho1now(2*ir-1,1,1) + kxc(ir,4)*gradrho_gradrho1
       coeffim_grho=kxc(ir,3)*rho1now(2*ir  ,1,1) + kxc(ir,4)*gradrho_gradrho1im
  !    Reuse the storage in rho1now
       rho1now(2*ir-1,1,2:4)=r1(:)  *kxc(ir,2)+r0(:)*coeff_grho
       rho1now(2*ir  ,1,2:4)=r1im(:)*kxc(ir,2)+r0(:)*coeffim_grho
     end do
   else
     do ir=1,nfft
       do ii=1,3  ! grad of spin-up ans spin_dn GS rho
         r0_up(ii)=kxc(ir,13+2*ii);r0_dn(ii)=kxc(ir,12+2*ii)-kxc(ir,13+2*ii)
       end do
       r0(:)=r0_up(:)+r0_dn(:)          ! grad of GS rho
       r1_up(:)=rho1now(2*ir-1,1,2:4)   ! grad of spin-up rho1
       r1im_up(:)=rho1now(2*ir,1,2:4)   ! grad of spin-up rho1 , im part
       r1_dn(:)=rho1now(2*ir-1,2,2:4)   ! grad of spin-down rho1
       r1im_dn(:)=rho1now(2*ir,2,2:4)   ! grad of spin-down rho1 , im part
       r1(:)=r1_up(:)+r1_dn(:)      ! grad of GS rho1
       r1im(:)=r1im_up(:)+r1im_dn(:)      ! grad of GS rho1, im part
       gradrho_gradrho1_up  =dot_product(r0_up,r1_up)
       gradrho_gradrho1_dn  =dot_product(r0_dn,r1_dn)
       gradrho_gradrho1     =dot_product(r0,r1)
       gradrho_gradrho1im_up=dot_product(r0_up,r1im_up)
       gradrho_gradrho1im_dn=dot_product(r0_dn,r1im_dn)
       gradrho_gradrho1im   =dot_product(r0,r1im)
       dnexcdn(2*ir-1,1)=kxc(ir, 1)*rho1now(2*ir-1,1,1) &
&       +kxc(ir, 2)*rho1now(2*ir-1,2,1) &
&       +kxc(ir, 6)*gradrho_gradrho1_up &
&       +kxc(ir,11)*gradrho_gradrho1
       dnexcdn(2*ir-1,2)=kxc(ir, 3)*rho1now(2*ir-1,2,1) &
&       +kxc(ir, 2)*rho1now(2*ir-1,1,1) &
&       +kxc(ir, 7)*gradrho_gradrho1_dn &
&       +kxc(ir,12)*gradrho_gradrho1
       dnexcdn(2*ir  ,1)=kxc(ir, 1)*rho1now(2*ir  ,1,1) &
&       +kxc(ir, 2)*rho1now(2*ir  ,2,1) &
&       +kxc(ir, 6)*gradrho_gradrho1im_up &
&       +kxc(ir,11)*gradrho_gradrho1im
       dnexcdn(2*ir  ,2)=kxc(ir, 3)*rho1now(2*ir  ,2,1) &
&       +kxc(ir, 2)*rho1now(2*ir  ,1,1) &
&       +kxc(ir, 7)*gradrho_gradrho1im_dn &
&       +kxc(ir,12)*gradrho_gradrho1im
       coeff_grho_corr  =kxc(ir,11)*rho1now(2*ir-1,1,1) &
&       +kxc(ir,12)*rho1now(2*ir-1,2,1) &
&       +kxc(ir,13)*gradrho_gradrho1
       coeffim_grho_corr=kxc(ir,11)*rho1now(2*ir  ,1,1) &
&       +kxc(ir,12)*rho1now(2*ir  ,2,1) &
&       +kxc(ir,13)*gradrho_gradrho1im
       coeff_grho_up  =kxc(ir,6)*rho1now(2*ir-1,1,1)+kxc(ir,8)*gradrho_gradrho1_up
       coeff_grho_dn  =kxc(ir,7)*rho1now(2*ir-1,2,1)+kxc(ir,9)*gradrho_gradrho1_dn
       coeffim_grho_up=kxc(ir,6)*rho1now(2*ir  ,1,1)+kxc(ir,8)*gradrho_gradrho1im_up
       coeffim_grho_dn=kxc(ir,7)*rho1now(2*ir  ,2,1)+kxc(ir,9)*gradrho_gradrho1im_dn
!      Reuse the storage in rho1now
       rho1now(2*ir-1,1,2:4)=(kxc(ir,4)+kxc(ir,10))*r1_up(:) &
&       +kxc(ir,10)            *r1_dn(:) &
&       +coeff_grho_up         *r0_up(:) &
&       +coeff_grho_corr*r0(:)
       rho1now(2*ir-1,2,2:4)=(kxc(ir,5)+kxc(ir,10))*r1_dn(:) &
&       +kxc(ir,10)            *r1_up(:) &
&       +coeff_grho_dn         *r0_dn(:) &
&       +coeff_grho_corr*r0(:)
       rho1now(2*ir  ,1,2:4)=(kxc(ir,4)+kxc(ir,10))*r1im_up(:) &
&       +kxc(ir,10)            *r1im_dn(:) &
&       +coeffim_grho_up       *r0_up(:)   &
&       +coeffim_grho_corr     *r0(:)
       rho1now(2*ir  ,2,2:4)=(kxc(ir,5)+kxc(ir,10))*r1im_dn(:) &
&       +kxc(ir,10)            *r1im_up(:) &
&       +coeffim_grho_dn       *r0_dn(:)   &
&       +coeffim_grho_corr     *r0(:)
     end do
   end if ! nspden

 end if

 vxc1(:,:)=zero
 call xcpot(cplex,gprimd,ishift,use_laplacian,mpi_enreg,nfft,ngfft,ngrad,nspden,&
& nspgrad,qphon,depsxc=dnexcdn,rhonow=rho1now,vxc=vxc1)

!call filterpot(paral_kgb,cplex,gmet,gsqcut,nfft,ngfft,nspden,qphon,vxc1)

 ABI_DEALLOCATE(dnexcdn)
 ABI_DEALLOCATE(rho1now)

 DBG_EXIT("COLL")

end subroutine dfpt_mkvxcgga
!!***

!!****f* ABINIT/dfpt_mkvxc_noncoll
!! NAME
!! dfpt_mkvxc_noncoll
!!
!! FUNCTION
!! Compute the first-order change of exchange-correlation potential
!! due to atomic displacement for non-collinear spins: assemble the first-order
!! density change with the frozen-core density change, then use
!! the exchange-correlation kernel.
!!
!! INPUTS
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!         if 2, COMPLEX
!!  ixc= choice of exchange-correlation scheme
!!  kxc(nfft,nkxc)=exchange and correlation kernel (see rhotoxc.F90)
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/variables/vargs.htm#ngfft
!!  nhat(nfft,nspden*nhatdim)= -PAW only- GS compensation density
!!  nhatdim= -PAW only- 1 if nhat array is used ; 0 otherwise
!!  nhat1(cplex*nfft,nspden*nhat1dim)= -PAW only- 1st-order compensation density
!!  nhat1dim= -PAW only- 1 if nhat1 array is used ; 0 otherwise
!!  nhat1gr(cplex*nfft,nspden,3*nhat1grdim)= -PAW only- gradients of 1st-order compensation density
!!  nhat1grdim= -PAW only- 1 if nhat1gr array is used ; 0 otherwise
!!  nkxc=second dimension of the kxc array
!!  non_magnetic_xc= if true, handle density/potential as non-magnetic (even if it is)
!!  nspden=number of spin-density components
!!  n3xccc=dimension of xccc3d1 ; 0 if no XC core correction is used, otherwise, nfft
!!  optnc=option for non-collinear magnetism (nspden=4):
!!       1: the whole 2x2 Vres matrix is computed
!!       2: only Vres^{11} and Vres^{22} are computed
!!  option=if 0, work only with the XC core-correction,
!!         if 1, treat both density change and XC core correction
!!         if 2, treat only density change
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rhor(nfft,nspden)=GS electron density in real space
!!  rhor1(cplex*nfft,nspden)=1st-order electron density in real space
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  usexcnhat= -PAW only- 1 if nhat density has to be taken into account in Vxc
!!  vxc(nfft,nspden)=GS XC potential
!!
!!
!! OUTPUT
!!  vxc1(cplex*nfft,nspden)=change in exchange-correlation potential (including
!!   core-correction, if applicable)
!!
!! PARENTS
!!      dfpt_dyxc1,dfpt_nstdy,dfpt_nstpaw,dfpt_rhotov,nres2vres
!!
!! CHILDREN
!!      dfpt_mkvxc,rotate_back_mag,rotate_back_mag_dfpt,rotate_mag,timab
!!
!! SOURCE

subroutine dfpt_mkvxc_noncoll(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,nhat,nhatdim,nhat1,nhat1dim,&
&          nhat1gr,nhat1grdim,nkxc,non_magnetic_xc,nspden,n3xccc,optnc,option,qphon,&
&          rhor,rhor1,rprimd,usexcnhat,vxc,vxc1,xccc3d1,ixcrot)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,ixc,n3xccc,nfft,nhatdim,nhat1dim,nhat1grdim,optnc
 integer,intent(in) :: nkxc,nspden,option,usexcnhat
 logical,intent(in) :: non_magnetic_xc
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: nhat1gr(cplex*nfft,nspden,3*nhat1grdim)
 real(dp),intent(in) :: kxc(nfft,nkxc)
 real(dp),intent(in) :: vxc(nfft,nspden)
 real(dp),intent(in) :: nhat(nfft,nspden*nhatdim),nhat1(cplex*nfft,nspden*nhat1dim)
 real(dp),intent(in),target :: rhor(nfft,nspden),rhor1(cplex*nfft,nspden)
 real(dp),intent(in) :: qphon(3),rprimd(3,3),xccc3d1(cplex*n3xccc)
 real(dp),intent(out) :: vxc1(cplex*nfft,nspden)
 integer,optional,intent(in) :: ixcrot
!Local variables-------------------------------
!scalars
!arrays
 real(dp) :: nhat1_zero(0,0),nhat1gr_zero(0,0,0),tsec(2)
 real(dp),allocatable :: m_norm(:),rhor1_diag(:,:),vxc1_diag(:,:)
 real(dp), ABI_CONTIGUOUS pointer :: mag(:,:),rhor_(:,:),rhor1_(:,:)
! *************************************************************************

!  Non-collinear magnetism
!  Has to locally "rotate" rho(r)^(1) (according to magnetization),
!  Compute Vxc(r)^(1) in the spin frame aligned with \vec{m} and rotate it back

 DBG_ENTER("COLL")
 ABI_UNUSED(nhat1gr)

 call timab(181,1,tsec)

 if(nspden/=4) then
   MSG_BUG('only for nspden=4!')
 end if

 if(nkxc/=2*min(nspden,2)-1) then
   MSG_BUG('nspden=4 works only with LSDA.')
 end if

!Special case: no XC applied
 if (ixc==0.or.nkxc==0) then
   MSG_WARNING('Note that no xc is applied (ixc=0)')
   vxc1(:,:)=zero
   return
 end if



!Treat first LDA
 if(nkxc==1.or.nkxc==3)then

   vxc1(:,:)=zero

!  PAW: possibly substract compensation density
   if ((usexcnhat==0.and.nhatdim==1).or.(non_magnetic_xc)) then
     ABI_ALLOCATE(rhor_,(nfft,nspden))
     if (usexcnhat==0.and.nhatdim==1) then
       rhor_(:,:) =rhor(:,:)-nhat(:,:)
     else
       rhor_(:,:) =rhor(:,:)
     end if
     if (non_magnetic_xc) then
       if(nspden==2) rhor_(:,2)=rhor_(:,1)*half
       if(nspden==4) rhor_(:,2:4)=zero
     end if
   else
     rhor_ => rhor
   end if
   if ((usexcnhat==0.and.nhat1dim==1).or.(non_magnetic_xc)) then
     ABI_ALLOCATE(rhor1_,(cplex*nfft,nspden))
     if (usexcnhat==0.and.nhatdim==1) then
       rhor1_(:,:)=rhor1(:,:)-nhat1(:,:)
     else
       rhor1_(:,:)=rhor1(:,:)
     end if
     if (non_magnetic_xc) then
       if(nspden==2) rhor1_(:,2)=rhor1_(:,1)*half
       if(nspden==4) rhor1_(:,2:4)=zero
     end if
   else
     rhor1_ => rhor1
   end if

!  Magnetization
   mag => rhor_(:,2:4)
   ABI_ALLOCATE(rhor1_diag,(cplex*nfft,2))
   ABI_ALLOCATE(vxc1_diag,(cplex*nfft,2))
   ABI_ALLOCATE(m_norm,(nfft))

!  -- Rotate rho(r)^(1)
!  SPr: for option=0 the rhor is not used, only core density xccc3d1
!       rotate_mag is only to compute the m_norm
   call rotate_mag(rhor1_,rhor1_diag,mag,nfft,cplex,mag_norm_out=m_norm,&
&   rho_out_format=2)

!  -- Compute Vxc(r)^(1)=Kxc(r).rho(r)^(1)_rotated
!  Note for PAW: nhat has already been substracted; don't use it in dfpt_mkvxc
!                 (put all nhat options to zero).
!  The collinear routine dfpt_mkvxc wants a general density built as (tr[rho],rho_upup)
   call dfpt_mkvxc(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,nhat1_zero,0,nhat1gr_zero,0,&
&   nkxc,non_magnetic_xc,2,n3xccc,option,qphon,rhor1_diag,rprimd,0,vxc1_diag,xccc3d1)

   !call test_rotations(0,1)

!  -- Rotate back Vxc(r)^(1)
   if (optnc==1) then
     if(present(ixcrot)) then
       call rotate_back_mag_dfpt(option,vxc1_diag,vxc1,vxc,kxc,rhor1_,mag,nfft,cplex,&
&       mag_norm_in=m_norm,rot_method=ixcrot)
     else
       call rotate_back_mag_dfpt(option,vxc1_diag,vxc1,vxc,kxc,rhor1_,mag,nfft,cplex,&
&       mag_norm_in=m_norm)
     end if
   else
     call rotate_back_mag(vxc1_diag,vxc1,mag,nfft,mag_norm_in=m_norm)
     vxc1(:,3:4)=zero
   end if

   ABI_DEALLOCATE(rhor1_diag)
   ABI_DEALLOCATE(vxc1_diag)
   ABI_DEALLOCATE(m_norm)
   if ((usexcnhat==0.and.nhatdim==1).or.(non_magnetic_xc)) then
     ABI_DEALLOCATE(rhor_)
   end if
   if ((usexcnhat==0.and.nhat1dim==1).or.(non_magnetic_xc)) then
     ABI_DEALLOCATE(rhor1_)
   end if

 end if ! nkxc=1 or nkxc=3

 call timab(181,2,tsec)

 DBG_EXIT("COLL")

end subroutine dfpt_mkvxc_noncoll
!!***

end module m_dfpt_mkvxc
!!***
