!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spacepar
!! NAME
!! m_spacepar
!!
!! FUNCTION
!!  Relatively Low-level procedures operating on arrays defined on the FFT box (G- or R- space)
!!  Unlike the procedures in m_cgtools, the routines declared in this module can use mpi_type.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (XG, BA, MT, DRH, DCA, GMR, MJV, JWZ)
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

module m_spacepar

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_sort
 use m_hightemp

 use m_time,            only : timab
 use defs_abitypes,     only : MPI_type
 use m_symtk,           only : mati3inv, chkgrp, symdet, symatm, matr3inv
 use m_geometry,        only : metric, symredcart
 use m_mpinfo,          only : ptabs_fourdp
 use m_fft,             only : zerosym, fourdp

 implicit none

 private
!!***

public :: hartre            ! Given rho(G), compute Hartree potential (=FFT of rho(G)/pi/(G+q)**2)
public :: make_vectornd     ! compute vector potential due to nuclear magnetic dipoles, in real space
public :: meanvalue_g       ! Compute <wf|op|wf> where op is real and diagonal in G-space.
public :: laplacian         ! Compute the laplacian of a function defined in real space
public :: redgr             ! Compute reduced gradients of a real function on the usual unshifted FFT grid.
public :: hartrestr         ! FFT of (rho(G)/pi)*[d(1/G**2)/d(strain) - delta(diagonal strain)*(1/G**2)]
public :: symrhg            ! Symmetrize rhor(r)
public :: irrzg             ! Find the irreducible zone in reciprocal space (used by symrhg)
public :: setsym            ! Set up irreducible zone in  G space by direct calculation.

! MG FIXME This routine is deprecated. Now the symmetrization of the **potentials** is done in the m_dvdb
public :: rotate_rho
!!***

contains
!!***

!!****f* m_spacepar/make_vectornd
!! NAME
!! make_vectornd
!!
!! FUNCTION
!! For nuclear dipole moments m, compute vector potential A(r) = (m x (r-R))/|r-R|^3
!! in r space. This is done by computing A(G) followed by FFT.
!!
!! NOTES
!! This code is copied and modified from m_spacepar/hartre where a very similar loop
!! over G is done followed by FFT to real space
!!
!! INPUTS
!!
!! OUTPUT
!!  vectornd(3,nfft)=Vector potential in real space, along Cartesian directions
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine make_vectornd(cplex,gsqcut,izero,mpi_enreg,natom,nfft,ngfft,nucdipmom,&
     & rprimd,vectornd,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,izero,natom,nfft
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: nucdipmom(3,natom),rprimd(3,3),xred(3,natom)
 real(dp),intent(out) :: vectornd(nfft,3)

!Local variables-------------------------------
 !scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i2_local,i23,i3,iatom,id1,id2,id3,ig,ig1,ig2,ig3,ig1max,ig2max,ig3max
 integer :: ig1min,ig2min,ig3min
 integer :: ii,ii1,ing,me_fft,n1,n2,n3,nd_atom,nd_atom_tot,nproc_fft
 real(dp),parameter :: tolfix=1.000000001e0_dp
 real(dp) :: cutoff,gqgm12,gqg2p3,gqgm23,gqgm13,gs2,gs3,gs
 real(dp) :: phase,prefac,precosph,presinph,prefacgs,ucvol
 !arrays
 integer :: id(3)
 integer,allocatable :: nd_list(:)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: gcart(3),gmet(3,3),gprimd(3,3),gred(3),mcg_cart(3),rmet(3,3)
 real(dp) :: AGre_red(3),AGre_cart(3),AGim_red(3),AGim_cart(3)
 real(dp),allocatable :: gq(:,:),nd_m(:,:),ndvecr(:),work1(:,:),work2(:,:),work3(:,:)

! *************************************************************************

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)


 ! make list of atoms with nonzero nuclear dipole moments
 ! in typical applications only 0 or 1 atoms have nonzero dipoles. This
 ! code shouldn't even be called if all dipoles are zero.
 nd_atom_tot = 0
 do iatom = 1, natom
    if (any(abs(nucdipmom(:,iatom))>tol8)) then
       nd_atom_tot = nd_atom_tot + 1
    end if
 end do

 ! note that nucdipmom is input as vectors in atomic units referenced
 ! to cartesian coordinates
 ABI_ALLOCATE(nd_list,(nd_atom_tot))
 ABI_ALLOCATE(nd_m,(3,nd_atom_tot))
 nd_atom_tot = 0
 do iatom = 1, natom
    if (any(abs(nucdipmom(:,iatom))>tol8)) then
       nd_atom_tot = nd_atom_tot + 1
       nd_list(nd_atom_tot) = iatom
       nd_m(:,nd_atom_tot) = nucdipmom(:,iatom)
    end if
 end do

 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 nproc_fft = mpi_enreg%nproc_fft; me_fft = mpi_enreg%me_fft

 prefac = -four_pi/(ucvol*two_pi)

 ! Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

 ! Initialize a few quantities
 cutoff=gsqcut*tolfix

 ! In order to speed the routine, precompute the components of g+q
 ! Also check if the booked space was large enough...
 ABI_ALLOCATE(gq,(3,max(n1,n2,n3)))
 do ii=1,3
   id(ii)=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
     ig=ing-(ing/id(ii))*ngfft(ii)-1
     gq(ii,ing)=ig
   end do
 end do
 ig1max=-1;ig2max=-1;ig3max=-1
 ig1min=n1;ig2min=n2;ig3min=n3

 ABI_ALLOCATE(work1,(2,nfft))
 ABI_ALLOCATE(work2,(2,nfft))
 ABI_ALLOCATE(work3,(2,nfft))
 work1=zero; work2=zero; work3=zero
 id1=n1/2+2;id2=n2/2+2;id3=n3/2+2

 ! Triple loop on each dimension
 do i3=1,n3
   ig3=i3-(i3/id3)*n3-1
   ! Precompute some products that do not depend on i2 and i1
   gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
   gqgm23=gq(3,i3)*gmet(2,3)*2
   gqgm13=gq(3,i3)*gmet(1,3)*2

   do i2=1,n2
     ig2=i2-(i2/id2)*n2-1
     if (fftn2_distrib(i2) == me_fft) then
       gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
       gqgm12=gq(2,i2)*gmet(1,2)*2
       gqg2p3=gqgm13+gqgm12

       i2_local = ffti2_local(i2)
       i23=n1*(i2_local-1 +(n2/nproc_fft)*(i3-1))
       ! Do the test that eliminates the Gamma point outside of the inner loop
       ii1=1
       if(i23==0 .and. ig2==0 .and. ig3==0)then
         ii1=2
         work1(re,1+i23)=zero
         work1(im,1+i23)=zero
       end if

       ! Final inner loop on the first dimension (note the lower limit)
       do i1=ii1,n1
          gs=gs2+ gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
          ig1 = i1 - (i1/id1)*n1 -1
          ii=i1+i23

          gred(1) = one*ig1; gred(2) = one*ig2; gred(3) = one*ig3
          ! obtain \vec{G} in cartesian coordinates
          gcart(1:3) = MATMUL(gprimd,gred)
          gs = DOT_PRODUCT(gred,MATMUL(gmet,gred))

         if(gs .LE. cutoff)then

            prefacgs = prefac/gs
            do iatom = 1, nd_atom_tot
               nd_atom = nd_list(iatom)
               phase = two_pi*DOT_PRODUCT(xred(:,nd_atom),gred(:))
               presinph=prefacgs*sin(phase)
               precosph=prefacgs*cos(phase)

               ! cross product m x G
               mcg_cart(1) =  nd_m(2,iatom)*gcart(3) - nd_m(3,iatom)*gcart(2)
               mcg_cart(2) = -nd_m(1,iatom)*gcart(3) + nd_m(3,iatom)*gcart(1)
               mcg_cart(3) =  nd_m(1,iatom)*gcart(2) - nd_m(2,iatom)*gcart(1)

               ! Re(A(G)), in cartesian coordinates
               AGre_cart = presinph*mcg_cart

               ! refer back to recip space
               AGre_red = MATMUL(TRANSPOSE(gprimd),AGre_cart)

               ! Im(A(G)), in cartesian coordinates
               AGim_cart = precosph*mcg_cart

               ! refer back to recip space
               AGim_red = MATMUL(TRANSPOSE(gprimd),AGim_cart)

               work1(re,ii) = AGre_red(1)
               work2(re,ii) = AGre_red(2)
               work3(re,ii) = AGre_red(3)
               work1(im,ii) = AGim_red(1)
               work2(im,ii) = AGim_red(2)
               work3(im,ii) = AGim_red(3)
            end do
         else
           ! gs>cutoff
           work1(re,ii)=zero
           work1(im,ii)=zero
           work2(re,ii)=zero
           work2(im,ii)=zero
           work3(re,ii)=zero
           work3(im,ii)=zero
         end if

       end do ! End loop on i1
     end if
   end do ! End loop on i2
 end do ! End loop on i3

 ABI_DEALLOCATE(gq)
 ABI_DEALLOCATE(nd_list)
 ABI_DEALLOCATE(nd_m)

 if ( izero .EQ. 1 ) then
   ! Set contribution of unbalanced components to zero

    call zerosym(work1,2,n1,n2,n3,comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
    call zerosym(work2,2,n1,n2,n3,comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
    call zerosym(work3,2,n1,n2,n3,comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)

 end if

 ! Fourier Transform
 ABI_ALLOCATE(ndvecr,(cplex*nfft))
 ndvecr=zero
 call fourdp(cplex,work1,ndvecr,1,mpi_enreg,nfft,1,ngfft,0)
 vectornd(:,1)=ndvecr(:)
 ABI_DEALLOCATE(work1)

 ndvecr=zero
 call fourdp(cplex,work2,ndvecr,1,mpi_enreg,nfft,1,ngfft,0)
 vectornd(:,2) = ndvecr(:)
 ABI_DEALLOCATE(work2)

 ndvecr=zero
 call fourdp(cplex,work3,ndvecr,1,mpi_enreg,nfft,1,ngfft,0)
 vectornd(:,3) = ndvecr(:)
 ABI_DEALLOCATE(work3)
 ABI_DEALLOCATE(ndvecr)

end subroutine make_vectornd
!!***

!!****f* m_spacepar/hartre
!! NAME
!! hartre
!!
!! FUNCTION
!! Given rho(G), compute Hartree potential (=FFT of rho(G)/pi/(G+q)**2)
!! When cplex=1, assume q=(0 0 0), and vhartr will be REAL
!! When cplex=2, q must be taken into account, and vhartr will be COMPLEX
!!
!! NOTES
!! *Modified code to avoid if statements inside loops to skip G=0.
!!  Replaced if statement on G^2>gsqcut to skip G s outside where
!!  rho(G) should be 0.  Effect is negligible but gsqcut should be
!!  used to be strictly consistent with usage elsewhere in code.
!! *The speed-up is provided by doing a few precomputations outside
!!  the inner loop. One variable size array is needed for this (gq).
!!
!! INPUTS
!!  cplex= if 1, vhartr is REAL, if 2, vhartr is COMPLEX
!!  gsqcut=cutoff value on G**2 for sphere inside fft box.
!!         (gsqcut=(boxcut**2)*ecut/(2.d0*(Pi**2))
!!  izero=if 1, unbalanced components of Vhartree(g) are set to zero
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  [qpt(3)=reduced coordinates for a wavevector to be combined with the G vectors (needed if cplex==2).]
!!  rhog(2,nfft)=electron density in G space
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  divgq0= [optional argument] value of the integration of the Coulomb singularity 4pi\int_BZ 1/q^2 dq
!!
!! OUTPUT
!!  vhartr(cplex*nfft)=Hartree potential in real space, either REAL or COMPLEX
!!
!! PARENTS
!!      dfpt_rhotov,dfptnl_loop,energy,fock_getghc,m_kxc,nonlinear,nres2vres
!!      odamix,prcref,prcref_PMA,respfn,rhotov,setup_positron,setvtr,tddft
!!
!! CHILDREN
!!      fourdp,metric,ptabs_fourdp,timab,zerosym
!!
!! SOURCE

subroutine hartre(cplex,gsqcut,izero,mpi_enreg,nfft,ngfft,rhog,rprimd,vhartr,&
&  divgq0,qpt) ! Optional argument

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,izero,nfft
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rprimd(3,3),rhog(2,nfft)
 real(dp),intent(in),optional :: divgq0
 real(dp),intent(in),optional :: qpt(3)
 real(dp),intent(out) :: vhartr(cplex*nfft)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i23,i2_local,i3,id1,id2,id3
 integer :: ig,ig1min,ig1,ig1max,ig2,ig2min,ig2max,ig3,ig3min,ig3max
 integer :: ii,ii1,ing,n1,n2,n3,qeq0,qeq05,me_fft,nproc_fft
 real(dp),parameter :: tolfix=1.000000001e0_dp
 real(dp) :: cutoff,den,gqg2p3,gqgm12,gqgm13,gqgm23,gs,gs2,gs3,ucvol
 character(len=500) :: message
!arrays
 integer :: id(3)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: gmet(3,3),gprimd(3,3),qpt_(3),rmet(3,3),tsec(2)
 real(dp),allocatable :: gq(:,:),work1(:,:)

! *************************************************************************

 ! Keep track of total time spent in hartre
 call timab(10,1,tsec)

 ! Check that cplex has an allowed value
 if(cplex/=1 .and. cplex/=2)then
   write(message, '(a,i0,a,a)' )&
   'From the calling routine, cplex=',cplex,ch10,&
   'but the only value allowed are 1 and 2.'
   MSG_BUG(message)
 end if

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 nproc_fft = mpi_enreg%nproc_fft; me_fft = mpi_enreg%me_fft

 ! Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

 ! Initialize a few quantities
 cutoff=gsqcut*tolfix
 if(present(qpt))then
   qpt_=qpt
 else
   qpt_=zero
 end if
 qeq0=0
 if(qpt_(1)**2+qpt_(2)**2+qpt_(3)**2<1.d-15) qeq0=1
 qeq05=0
 if (qeq0==0) then
   if (abs(abs(qpt_(1))-half)<tol12.or.abs(abs(qpt_(2))-half)<tol12.or.abs(abs(qpt_(3))-half)<tol12) qeq05=1
 end if

 ! If cplex=1 then qpt_ should be 0 0 0
 if (cplex==1.and. qeq0/=1) then
   write(message,'(a,3e12.4,a,a)')&
   'cplex=1 but qpt=',qpt_,ch10,&
   'qpt should be 0 0 0.'
   MSG_BUG(message)
 end if

 ! If FFT parallelism then qpt should not be 1/2
 if (nproc_fft>1.and.qeq05==1) then
   write(message, '(a,3e12.4,a,a)' )&
   'FFT parallelism selected but qpt',qpt_,ch10,&
   'qpt(i) should not be 1/2...'
   MSG_ERROR(message)
 end if

 ! In order to speed the routine, precompute the components of g+q
 ! Also check if the booked space was large enough...
 ABI_ALLOCATE(gq,(3,max(n1,n2,n3)))
 do ii=1,3
   id(ii)=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
     ig=ing-(ing/id(ii))*ngfft(ii)-1
     gq(ii,ing)=ig+qpt_(ii)
   end do
 end do
 ig1max=-1;ig2max=-1;ig3max=-1
 ig1min=n1;ig2min=n2;ig3min=n3

 ABI_ALLOCATE(work1,(2,nfft))
 id1=n1/2+2;id2=n2/2+2;id3=n3/2+2

 ! Triple loop on each dimension
 do i3=1,n3
   ig3=i3-(i3/id3)*n3-1
   ! Precompute some products that do not depend on i2 and i1
   gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
   gqgm23=gq(3,i3)*gmet(2,3)*2
   gqgm13=gq(3,i3)*gmet(1,3)*2

   do i2=1,n2
     ig2=i2-(i2/id2)*n2-1
     if (fftn2_distrib(i2) == me_fft) then
       gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
       gqgm12=gq(2,i2)*gmet(1,2)*2
       gqg2p3=gqgm13+gqgm12

       i2_local = ffti2_local(i2)
       i23=n1*(i2_local-1 +(n2/nproc_fft)*(i3-1))
       ! Do the test that eliminates the Gamma point outside of the inner loop
       ii1=1
       if(i23==0 .and. qeq0==1  .and. ig2==0 .and. ig3==0)then
         ii1=2
         work1(re,1+i23)=zero
         work1(im,1+i23)=zero
         ! If the value of the integration of the Coulomb singularity 4pi\int_BZ 1/q^2 dq is given, use it
         if (PRESENT(divgq0)) then
           work1(re,1+i23)=rhog(re,1+i23)*divgq0*piinv
           work1(im,1+i23)=rhog(im,1+i23)*divgq0*piinv
         end if
       end if

       ! Final inner loop on the first dimension (note the lower limit)
       do i1=ii1,n1
         gs=gs2+ gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
         ii=i1+i23

         if(gs<=cutoff)then
           ! Identify min/max indexes (to cancel unbalanced contributions later)
           ! Count (q+g)-vectors with similar norm
           if ((qeq05==1).and.(izero==1)) then
             ig1=i1-(i1/id1)*n1-1
             ig1max=max(ig1max,ig1); ig1min=min(ig1min,ig1)
             ig2max=max(ig2max,ig2); ig2min=min(ig2min,ig2)
             ig3max=max(ig3max,ig3); ig3min=min(ig3min,ig3)
           end if

           den=piinv/gs
           work1(re,ii)=rhog(re,ii)*den
           work1(im,ii)=rhog(im,ii)*den
         else
           ! gs>cutoff
           work1(re,ii)=zero
           work1(im,ii)=zero
         end if

       end do ! End loop on i1
     end if
   end do ! End loop on i2
 end do ! End loop on i3

 ABI_DEALLOCATE(gq)

 if (izero==1) then
   ! Set contribution of unbalanced components to zero

   if (qeq0==1) then !q=0
     call zerosym(work1,2,n1,n2,n3,comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)

   else if (qeq05==1) then
     !q=1/2; this doesn't work in parallel
     ig1=-1;if (mod(n1,2)==0) ig1=1+n1/2
     ig2=-1;if (mod(n2,2)==0) ig2=1+n2/2
     ig3=-1;if (mod(n3,2)==0) ig3=1+n3/2
     if (abs(abs(qpt_(1))-half)<tol12) then
       if (abs(ig1min)<abs(ig1max)) ig1=abs(ig1max)
       if (abs(ig1min)>abs(ig1max)) ig1=n1-abs(ig1min)
     end if
     if (abs(abs(qpt_(2))-half)<tol12) then
       if (abs(ig2min)<abs(ig2max)) ig2=abs(ig2max)
       if (abs(ig2min)>abs(ig2max)) ig2=n2-abs(ig2min)
     end if
     if (abs(abs(qpt_(3))-half)<tol12) then
       if (abs(ig3min)<abs(ig3max)) ig3=abs(ig3max)
       if (abs(ig3min)>abs(ig3max)) ig3=n3-abs(ig3min)
     end if
     call zerosym(work1,2,n1,n2,n3,ig1=ig1,ig2=ig2,ig3=ig3,&
       comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
   end if
 end if

 ! Fourier Transform Vhartree. Vh in reciprocal space was stored in work1
 call fourdp(cplex,work1,vhartr,1,mpi_enreg,nfft,1,ngfft,0)

 ABI_DEALLOCATE(work1)

 call timab(10,2,tsec)

end subroutine hartre
!!***

!!****f* m_spacepar/meanvalue_g
!! NAME
!! meanvalue_g
!!
!! FUNCTION
!!  Compute the mean value of one wavefunction, in reciprocal space,
!!  for an operator that is real, diagonal in G-space: <wf|op|wf>
!!  For the time being, only spin-independent operators are treated.
!!
!! INPUTS
!!  diag(npw)=diagonal operator (real, spin-independent!)
!!  filter= if 1, need to filter on the value of diag, that must be less than huge(0.0d0)*1.d-11
!!          otherwise, should be 0
!!  istwf_k=storage mode of the vectors
!!  npw=number of planewaves of the vector
!!  nspinor=number of spinor components
!!  vect(2,npw*nspinor)=vector
!!  vect1(2,npw*nspinor*use_ndo)=vector1 (=vector in most of the cases)
!!  use_ndo = says if vect=/vect1
!!
!! OUTPUT
!!  ar=mean value
!!
!! PARENTS
!!      dfpt_vtowfk,energy,forstrnps,vtowfk
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

subroutine meanvalue_g(ar,diag,filter,istwf_k,mpi_enreg,npw,nspinor,vect,vect1,use_ndo,ar_im)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: filter,istwf_k,npw,nspinor,use_ndo
 real(dp),intent(out) :: ar
 real(dp),intent(out),optional :: ar_im
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 real(dp),intent(in) :: diag(npw),vect(2,npw*nspinor)
 real(dp),intent(in) :: vect1(2,npw*nspinor)

!Local variables-------------------------------
!scalars
 integer :: i1,ierr,ipw,jpw,me_g0
 character(len=500) :: message
!arrays

! *************************************************************************
 me_g0 = mpi_enreg%me_g0

 DBG_CHECK(ANY(filter==(/0,1/)),"Wrong filter")
 DBG_CHECK(ANY(nspinor==(/1,2/)),"Wrong nspinor")
 DBG_CHECK(ANY(istwf_k==(/(ipw,ipw=1,9)/)),"Wrong istwf_k")

 if(nspinor==2 .and. istwf_k/=1)then
   write(message,'(a,a,a,i6,a,i6)')&
&   'When istwf_k/=1, nspinor must be 1,',ch10,&
&   'however, nspinor=',nspinor,', and istwf_k=',istwf_k
   MSG_BUG(message)
 end if

 if(use_ndo==1 .and. (istwf_k==2 .and.me_g0==1)) then
   MSG_BUG('use_ndo==1, not tested, use istwfk=1')
 end if

 ar=zero
 if(present(ar_im)) ar_im=zero

!Normal storage mode
 if(istwf_k==1)then

!  No filter
   if(filter==0)then
!$OMP PARALLEL DO REDUCTION(+:ar)
     do ipw=1,npw
       ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
     end do
     if(nspinor==2)then
!$OMP PARALLEL DO REDUCTION(+:ar) PRIVATE(jpw)
       do ipw=1+npw,2*npw
         jpw=ipw-npw
         ar=ar+diag(jpw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
       end do
     end if
     if(use_ndo==1)then
!$OMP PARALLEL DO REDUCTION(+:ar_im)
       do ipw=1,npw
         ar_im=ar_im+diag(ipw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
       end do
       if(nspinor == 2) then
!$OMP PARALLEL DO REDUCTION(+:ar_im) PRIVATE(jpw)
         do ipw=1+npw,2*npw
           jpw=ipw-npw
           ar_im=ar_im+diag(jpw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
         end do
       end if
     end if

!    !$OMP PARALLEL DO REDUCTION(+:ar,ar_im)
!    do ipw=1,npw
!    ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
!    if(use_ndo==1.and.nspinor==2) ar_im=ar_im+diag(ipw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
!    end do
!    if(nspinor==2)then
!    !$OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ar,ar_im)
!    do ipw=1+npw,2*npw
!    ar=ar+diag(ipw-npw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
!    if(use_ndo==1.and.nspinor==2) ar_im=ar_im+diag(ipw-npw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
!    end do
!    end if
   else ! will filter

!$OMP PARALLEL DO REDUCTION(+:ar)
     do ipw=1,npw
       if(diag(ipw)<huge(0.0d0)*1.d-11)then
         ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
       end if
     end do
     if(nspinor==2)then
!$OMP PARALLEL DO REDUCTION(+:ar) PRIVATE(jpw)
       do ipw=1+npw,2*npw
         jpw=ipw-npw
         if(diag(jpw)<huge(0.0d0)*1.d-11)then
           ar=ar+diag(jpw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
         end if
       end do
     end if
     if(use_ndo==1)then
       if(.not.present(ar_im)) then
         MSG_BUG("use_ndo true and ar_im not present")
       end if
!$OMP PARALLEL DO REDUCTION(+:ar_im)
       do ipw=1,npw
         if(diag(ipw)<huge(0.0d0)*1.d-11)then
           ar_im=ar_im+diag(ipw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
         end if
       end do
       if(nspinor == 2) then
!$OMP PARALLEL DO REDUCTION(+:ar_im) PRIVATE(jpw)
         do ipw=1+npw,2*npw
           jpw=ipw-npw
           if(diag(jpw)<huge(0.0d0)*1.d-11)then
             ar_im=ar_im+diag(jpw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
           end if
         end do
       end if
     end if


!    !$OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ar,ar_im)
!    do ipw=1,npw
!    if(diag(ipw)<huge(0.0d0)*1.d-11)then
!    ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
!    if(use_ndo==1.and.nspinor==2) ar_im=ar_im+diag(ipw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
!    end if
!    end do
!    if(nspinor==2)then
!    !$OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ar,ar_im)
!    do ipw=1+npw,2*npw
!    if(diag(ipw-npw)<huge(0.0d0)*1.d-11)then
!    ar=ar+diag(ipw-npw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
!    if(use_ndo==1.and.nspinor==2) ar_im=ar_im+diag(ipw-npw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
!    end if
!    end do
!    end if ! nspinor==2

   end if ! filter==0

 else if(istwf_k>=2)then

   if(filter==0)then
     i1=1
     if(istwf_k==2 .and. me_g0==1)then ! MPIWF need to know which proc has G=0
       ar=half*diag(1)*vect(1,1)*vect1(1,1) ; i1=2
     end if

!$OMP PARALLEL DO REDUCTION(+:ar)
     do ipw=i1,npw
       ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
     end do

   else ! filter/=0
     i1=1
     if(istwf_k==2 .and. me_g0==1)then
       if(diag(1)<huge(0.0d0)*1.d-11)then
         ar=half*diag(1)*vect(1,1)*vect1(1,1) ; i1=2
       end if
     end if

!$OMP PARALLEL DO REDUCTION(+:ar)
     do ipw=i1,npw
       if(diag(ipw)<huge(0.0d0)*1.d-11)then
         ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
       end if
     end do
   end if ! filter==0

   ar=two*ar

 end if ! istwf_k

!MPIWF need to make reduction on ar and ai .
 if(mpi_enreg%paral_kgb==1)then
   call xmpi_sum(ar,mpi_enreg%comm_bandspinorfft ,ierr)
   if(present(ar_im))then
     call xmpi_sum(ar_im,mpi_enreg%comm_bandspinorfft,ierr)
   end if
 end if

end subroutine meanvalue_g
!!***

!!****f* m_spacepar/laplacian
!! NAME
!! laplacian
!!
!! FUNCTION
!! compute the laplacian of a function defined in real space
!! the code is written in the way of /3xc/xcden.F90
!!
!! INPUTS
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  mpi_enreg=information about MPI parallelization
!!  nfft=number of points of the fft grid
!!  nfunc=number of functions on the grid for which the laplacian is to be calculated
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  (optional) rdfuncr(nfft,nfunc)=real(dp) discretized functions in real space
!!  rdfuncg_in TO BE DESCRIBED SB 090901
!!  laplacerdfuncg_in TO BE DESCRIBED SB 090901
!!  (optional) g2cart_in(nfft) = G**2 on the grid
!!
!! OUTPUT
!! (optional) laplacerdfuncr = laplacian in real space of the functions in rdfuncr
!! (optional) rdfuncg = real(dp) discretized functions in fourier space
!! (optional) laplacerdfuncg = real(dp) discretized laplacian of the functions in fourier space
!! (optional) g2cart_out(nfft) = G**2 on the grid
!!  rdfuncg_out TO BE DESCRIBED SB 090901
!!  laplacerdfuncg_out TO BE DESCRIBED SB 090901
!!
!! PARENTS
!!      frskerker1,frskerker2,moddiel_csrb,prcrskerker1,prcrskerker2
!!
!! CHILDREN
!!      fourdp,ptabs_fourdp
!!
!! SOURCE

subroutine laplacian(gprimd,mpi_enreg,nfft,nfunc,ngfft,rdfuncr,&
&  laplacerdfuncr,rdfuncg_out,laplacerdfuncg_out,g2cart_out,rdfuncg_in,g2cart_in)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nfunc
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(inout),optional :: laplacerdfuncr(nfft,nfunc)
 real(dp),intent(inout),optional,target :: rdfuncr(nfft,nfunc)
 real(dp),intent(in),optional,target :: g2cart_in(nfft) !vz_i
 real(dp),intent(out),optional,target :: g2cart_out(nfft)  !vz_i
 real(dp),intent(out),optional,target :: laplacerdfuncg_out(2,nfft,nfunc)
 real(dp),intent(in),optional,target :: rdfuncg_in(2,nfft,nfunc) !vz_i
 real(dp),intent(out),optional,target :: rdfuncg_out(2,nfft,nfunc)

!Local variables-------------------------------
!scalars
 integer :: count,i1,i2,i3,id1,id2,id3,ifft,ifunc,ig1,ig2,ig3,ii1,n1,n2
 integer :: n3
 real(dp) :: b11,b12,b13,b21,b22,b23,b31,b32,b33
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp),ABI_CONTIGUOUS pointer :: g2cart(:),laplacerdfuncg(:,:,:),rdfuncg(:,:,:)

! *************************************************************************

!Keep local copy of fft dimensions
 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)

 if(present(laplacerdfuncg_out)) then
   laplacerdfuncg => laplacerdfuncg_out
 else
   ABI_ALLOCATE(laplacerdfuncg,(2,nfft,nfunc))
 end if

 ! Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!change the real density rdfuncr on real space on the real density
!rdfuncg in reciprocal space
 if(.not.present(rdfuncg_in)) then
   if(present(rdfuncg_out)) then
     rdfuncg => rdfuncg_out
   else
     ABI_ALLOCATE(rdfuncg,(2,nfft,nfunc))
   end if
   if(present(rdfuncr)) then
     do ifunc=1,nfunc
       call fourdp(1,rdfuncg(:,:,ifunc),rdfuncr(:,ifunc),-1,mpi_enreg,nfft,1,ngfft,0)
     end do
   end if
 else
   rdfuncg => rdfuncg_in
 end if

!apply the laplacian on laplacerdfuncr
!code from /3xc/xcden.F90
!see also src/5common/hatre.F90 and src/5common/moddiel.F90
!Keep local copy of fft dimensions
!Initialize computation of G^2 in cartesian coordinates
 if(.not.present(g2cart_in)) then
   if(present(g2cart_out)) then
     g2cart => g2cart_out
   else
     ABI_ALLOCATE(g2cart,(nfft))
   end if
   id1=int(n1/2)+2
   id2=int(n2/2)+2
   id3=int(n3/2)+2
   count=0
   do i3=1,n3
     ifft=(i3-1)*n1*(n2/mpi_enreg%nproc_fft)
     ig3=i3-int(i3/id3)*n3-1
     do i2=1,n2
       if (fftn2_distrib(i2)==mpi_enreg%me_fft) then
         ig2=i2-int(i2/id2)*n2-1

         ii1=1
         do i1=ii1,n1
           ig1=i1-int(i1/id1)*n1-1
           ifft=ifft+1

           b11=gprimd(1,1)*real(ig1,dp)
           b21=gprimd(2,1)*real(ig1,dp)
           b31=gprimd(3,1)*real(ig1,dp)
           b12=gprimd(1,2)*real(ig2,dp)
           b22=gprimd(2,2)*real(ig2,dp)
           b32=gprimd(3,2)*real(ig2,dp)
           b13=gprimd(1,3)*real(ig3,dp)
           b23=gprimd(2,3)*real(ig3,dp)
           b33=gprimd(3,3)*real(ig3,dp)

           g2cart(ifft)=( &
&           (b11+b12+b13)**2&
&           +(b21+b22+b23)**2&
&           +(b31+b32+b33)**2&
&           )
           do ifunc=1,nfunc
!            compute the laplacian in Fourier space that is * (i x 2pi x G)**2
             laplacerdfuncg(1,ifft,ifunc) = -rdfuncg(1,ifft,ifunc)*g2cart(ifft)*two_pi*two_pi
             laplacerdfuncg(2,ifft,ifunc) = -rdfuncg(2,ifft,ifunc)*g2cart(ifft)*two_pi*two_pi
           end do
         end do
       end if
     end do
   end do
   if(.not.present(g2cart_out))  then
     ABI_DEALLOCATE(g2cart)
   end if
 else
   g2cart => g2cart_in
   do ifunc=1,nfunc
     do ifft=1,nfft
!      compute the laplacian in Fourier space that is * (i x 2pi x G)**2
       laplacerdfuncg(1,ifft,ifunc) = -rdfuncg(1,ifft,ifunc)*g2cart(ifft)*two_pi*two_pi
       laplacerdfuncg(2,ifft,ifunc) = -rdfuncg(2,ifft,ifunc)*g2cart(ifft)*two_pi*two_pi
     end do
   end do
 end if

!get the result back into real space
 if(present(laplacerdfuncr)) then
   do ifunc=1,nfunc
     call fourdp(1,laplacerdfuncg(:,:,ifunc),laplacerdfuncr(:,ifunc),1,mpi_enreg,nfft,1,ngfft,0)
   end do
 end if

!deallocate pointers
 if((.not.present(rdfuncg_in)).and.(.not.present(rdfuncg_in)))  then
   ABI_DEALLOCATE(rdfuncg)
 end if
 if(.not.present(laplacerdfuncg_out))  then
   ABI_DEALLOCATE(laplacerdfuncg)
 end if

end subroutine laplacian
!!***

!!****f* m_spacepar/redgr
!! NAME
!! redgr
!!
!! FUNCTION
!! Compute reduced gradients of a real function on the usual unshifted
!! fft grid. The gradient directions are the along the primitive
!! reciprocal lattice vectors.
!! The input function is intended to be a single spin component of
!! the valence charge density, the valence + core charge densities
!! or the first-order core charge density for use in frozen wf
!! elastic tensor calculations within the GGA.
!!
!! NOTES
!! Closely linked to xcden, but limited to Q=0, real charge densities,
!! and unshifted grids.
!!
!! INPUTS
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  frin(nfft)=real space input function
!!
!! OUTPUT
!!  frredgr(nfft,3)= reduced gradient of input function (same units as frin)
!!
!! PARENTS
!!      dfpt_eltfrxc
!!
!! CHILDREN
!!      fourdp
!!
!! SOURCE

subroutine redgr(frin,frredgr,mpi_enreg,nfft,ngfft)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: frin(nfft)
 real(dp),intent(out) :: frredgr(nfft,3)

!Local variables-------------------------------
!scalars
 integer :: cplex_tmp,i1,i2,i3,id,idir,ifft,ig,ii,ing,n1,n2,n3
!arrays
 real(dp),allocatable :: gg(:,:),wkcmpx(:,:),work(:),workgr(:,:)

! *************************************************************************

!Only real arrays are treated
 cplex_tmp=1

!Keep local copy of fft dimensions
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)

!In order to speed the routine, precompute the components of g, including 2pi factor
 ABI_ALLOCATE(gg,(max(n1,n2,n3),3))
 do ii=1,3
   id=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
     ig=ing-(ing/id)*ngfft(ii)-1
     gg(ing,ii)=two_pi*ig
   end do
!  Note that the G <-> -G symmetry must be maintained
   if(mod(ngfft(ii),2)==0)gg(ngfft(ii)/2+1,ii)=zero
 end do

 ABI_ALLOCATE(wkcmpx,(2,nfft))
 ABI_ALLOCATE(work,(nfft))
 ABI_ALLOCATE(workgr,(2,nfft))

!Obtain rho(G) in wkcmpx from input rho(r)
 work(:)=frin(:)

 call fourdp(cplex_tmp,wkcmpx,work,-1,mpi_enreg,nfft,1,ngfft,0)

!Gradient calculation for three reduced components in turn.
!Code duplicated to remove logic from loops.
 do idir=1,3
   if(idir==1) then
!$OMP PARALLEL DO PRIVATE(ifft)
     do i3=1,n3
       ifft=(i3-1)*n1*n2
       do i2=1,n2
         do i1=1,n1
           ifft=ifft+1
!          Multiply by i 2pi G(idir)
           workgr(2,ifft)= gg(i1,idir)*wkcmpx(1,ifft)
           workgr(1,ifft)=-gg(i1,idir)*wkcmpx(2,ifft)
         end do
       end do
     end do
   else if(idir==2) then
!$OMP PARALLEL DO PRIVATE(ifft)
     do i3=1,n3
       ifft=(i3-1)*n1*n2
       do i2=1,n2
         do i1=1,n1
           ifft=ifft+1
!          Multiply by i 2pi G(idir)
           workgr(2,ifft)= gg(i2,idir)*wkcmpx(1,ifft)
           workgr(1,ifft)=-gg(i2,idir)*wkcmpx(2,ifft)
         end do
       end do
     end do
   else
!$OMP PARALLEL DO PRIVATE(ifft)
     do i3=1,n3
       ifft=(i3-1)*n1*n2
       do i2=1,n2
         do i1=1,n1
           ifft=ifft+1
!          Multiply by i 2pi G(idir)
           workgr(2,ifft)= gg(i3,idir)*wkcmpx(1,ifft)
           workgr(1,ifft)=-gg(i3,idir)*wkcmpx(2,ifft)
         end do
       end do
     end do
   end if !idir

   call fourdp(cplex_tmp,workgr,work,1,mpi_enreg,nfft,1,ngfft,0)

!$OMP PARALLEL DO
   do ifft=1,nfft
     frredgr(ifft,idir)=work(ifft)
   end do

 end do !idir

 ABI_DEALLOCATE(gg)
 ABI_DEALLOCATE(wkcmpx)
 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(workgr)

end subroutine redgr
!!***

!!****f* m_spacepar/hartrestr
!! NAME
!! hartrestr
!!
!! FUNCTION
!! To be called for strain perturbation only
!! Compute the inhomogenous terms generated by the strain derivative of
!! Hartree potential due to the ground state charge rho(G)
!!
!!  FFT of (rho(G)/pi)*[d(1/G**2)/d(strain) - delta(diagonal strain)*(1/G**2)]
!!
!! NOTES
!! *based largely on hartre.f
!! *Modified code to avoid if statements inside loops to skip G=0.
!!  Replaced if statement on G^2>gsqcut to skip G s outside where
!!  rho(G) should be 0.  Effect is negligible but gsqcut should be
!!  used to be strictly consistent with usage elsewhere in code.
!! *The speed-up is provided by doing a few precomputations outside
!!  the inner loop. One variable size array is needed for this (gq).
!!
!! INPUTS
!!  gsqcut=cutoff value on G**2 for sphere inside fft box.
!!  idir=direction of the current perturbation
!!  ipert=type of the perturbation
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell.
!!  nfft=number of fft grid points (gsqcut=(boxcut**2)*ecut/(2._dp*(Pi**2))
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  rhog(2,nfft)=array for Fourier transform of GS electron density
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!
!! OUTPUT
!!  vhartr1(nfft)=Inhomogeneous term in strain-perturbation-induced Hartree
!!   potential in real space,
!!
!! PARENTS
!!      dfpt_nselt,dfpt_nstpaw,dfpt_rhotov
!!
!! CHILDREN
!!      fourdp,metric,ptabs_fourdp
!!
!! SOURCE

subroutine hartrestr(gsqcut,idir,ipert,mpi_enreg,natom,nfft,ngfft,rhog,rprimd,vhartr1)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idir,ipert,natom,nfft
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhog(2,nfft),rprimd(3,3)
 real(dp),intent(out) :: vhartr1(nfft)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i23,i3,id2,id3,ig,ig2,ig3,ii,ii1,ing,istr,ka,kb,n1,n2,n3
 real(dp),parameter :: tolfix=1.000000001_dp
 real(dp) :: cutoff,ddends,den,dgsds,gqg2p3,gqgm12,gqgm13,gqgm23,gs,gs2,gs3
 real(dp) :: term,ucvol
 character(len=500) :: message
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 integer :: id(3)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: dgmetds(3,3),gmet(3,3),gprimd(3,3),gqr(3),rmet(3,3)
 real(dp),allocatable :: gq(:,:),work1(:,:)

! *************************************************************************

 if( .not. (ipert==natom+3 .or. ipert==natom+4))then
   write(message, '(a,i0,a,a)' )&
&   'From the calling routine, ipert=',ipert,ch10,&
&   'so this routine for the strain perturbation should not be called.'
   MSG_BUG(message)
 end if

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)

!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!Initialize a few quantities
 cutoff=gsqcut*tolfix

 istr=idir + 3*(ipert-natom-3)

 if(istr<1 .or. istr>6)then
   write(message, '(a,i10,a,a,a)' )&
&   'Input dir gives istr=',istr,' not allowed.',ch10,&
&   'Possible values are 1,2,3,4,5,6 only.'
   MSG_BUG(message)
 end if

 ka=idx(2*istr-1);kb=idx(2*istr)
 do ii = 1,3
   dgmetds(:,ii)=-(gprimd(ka,:)*gprimd(kb,ii)+gprimd(kb,:)*gprimd(ka,ii))
 end do
!For historical reasons:
 dgmetds(:,:)=0.5_dp*dgmetds(:,:)

!In order to speed the routine, precompute the components of g+q
!Also check if the booked space was large enough...
 ABI_ALLOCATE(gq,(3,max(n1,n2,n3)))
 do ii=1,3
   id(ii)=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
     ig=ing-(ing/id(ii))*ngfft(ii)-1
     gq(ii,ing)=ig
   end do
 end do

 ABI_ALLOCATE(work1,(2,nfft))
 id2=n2/2+2
 id3=n3/2+2
!Triple loop on each dimension
 do i3=1,n3
   ig3=i3-(i3/id3)*n3-1
!  Precompute some products that do not depend on i2 and i1
   gqr(3)=gq(3,i3)
   gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
   gqgm23=gq(3,i3)*gmet(2,3)*2
   gqgm13=gq(3,i3)*gmet(1,3)*2

   do i2=1,n2
     if (fftn2_distrib(i2)==mpi_enreg%me_fft) then
       gqr(2)=gq(2,i2)
       gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
       gqgm12=gq(2,i2)*gmet(1,2)*2
       gqg2p3=gqgm13+gqgm12
       ig2=i2-(i2/id2)*n2-1
!      i23=n1*((i2-1)+n2*(i3-1))
       i23=n1*((ffti2_local(i2)-1)+(n2/mpi_enreg%nproc_fft)*(i3-1))
!      Do the test that eliminates the Gamma point outside
!      of the inner loop
       ii1=1
       if(i23==0  .and. ig2==0 .and. ig3==0)then
         ii1=2
         work1(re,1+i23)=0.0_dp
         work1(im,1+i23)=0.0_dp
       end if

!      Final inner loop on the first dimension
!      (note the lower limit)
       do i1=ii1,n1
         gs=gs2+ gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
         ii=i1+i23
         if(gs<=cutoff)then
           den=piinv/gs
           gqr(1)=gq(1,i1)
           dgsds=&
&           (gqr(1)*(dgmetds(1,1)*gqr(1)+dgmetds(1,2)*gqr(2)+dgmetds(1,3)*gqr(3))+  &
&           gqr(2)*(dgmetds(2,1)*gqr(1)+dgmetds(2,2)*gqr(2)+dgmetds(2,3)*gqr(3))+  &
&           gqr(3)*(dgmetds(3,1)*gqr(1)+dgmetds(3,2)*gqr(2)+dgmetds(3,3)*gqr(3)) )
           ddends=-piinv*dgsds/gs**2
           if(istr<=3)then
             term=2.0_dp*ddends-den
           else
             term=2.0_dp*ddends
           end if
           work1(re,ii)=rhog(re,ii)*term
           work1(im,ii)=rhog(im,ii)*term
         else
           work1(re,ii)=0.0_dp
           work1(im,ii)=0.0_dp
         end if

       end do ! End loop on i1
     end if
   end do ! End loop on i2
 end do !  End loop on i3

 ABI_DEALLOCATE(gq)

!Fourier Transform Vhartree.
!Vh in reciprocal space was stored in work1
 call fourdp(1,work1,vhartr1,1,mpi_enreg,nfft,1,ngfft,0)

 ABI_DEALLOCATE(work1)

end subroutine hartrestr
!!***

!!****f* m_spacepar/symrhg
!! NAME
!! symrhg
!!
!! FUNCTION
!! From rho(r), generate rho(G), symmetrize it, and
!! come back to the real space for a symmetrized rho(r).
!!
!! INPUTS
!! cplex=1 if rhor is real, 2 if rhor is complex
!! gprimd(3,3)=dimensional reciprocal space primitive translations
!! irrzon(nfft,2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!! mpi_enreg=information about MPI parallelization
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! nspden=number of spin-density components
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! nsym=number of symmetry elements.
!! phnons(2,nfft,(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!! rprimd(3,3)=dimensional real space primitive translations
!! symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym)=symmetry matrices in real space (integers)
!!
!! OUTPUT
!! rhog(2,nfft)=symmetrized rho(G) (total) electron density in G space
!!
!! SIDE EFFECTS
!! Input/Output
!! rhor(cplex*nfft,nspden)=array for electron density in electrons/bohr**3.
!! Input, but also output, if symmetrization is applied.
!! Also output if nspden > 1 (change spin components)
!!
!! NOTES
!! When using spin-polarization (nspden==2),
!! put total density in first half of rhor array and spin up in second half
!! If (nspden=2 and nsppol=2) the density is transformed as  (up,down) => (up+down,up)
!! If (nspden=2 and nsppol=1) anti-ferromagnetic symmetry operations
!!  must be used, such as to transform (2*up) => (up+down,up)
!! In spin-polarized, and if there is no symmetry to be
!! applied on the system, only the total density is generated in G space
!!
!! PARENTS
!!      dfpt_mkrho,dfpt_nstpaw,dfpt_rhofermi,dfpt_vtorho,m_dvdb,mkrho
!!      suscep_stat,vtorho,vtorhorec,vtorhotf,wfd_mkrho
!!
!! CHILDREN
!!      fourdp,matr3inv,ptabs_fourdp,symredcart,timab,xmpi_sum
!!
!! SOURCE

subroutine symrhg(cplex,gprimd,irrzon,mpi_enreg,nfft,nfftot,ngfft,nspden,nsppol,nsym,&
&                 phnons,rhog,rhor,rprimd,symafm,symrel)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nfftot,nspden,nsppol,nsym
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: irrzon(nfftot**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4)),ngfft(18)
 integer,intent(in) :: symafm(nsym),symrel(3,3,nsym)
 real(dp),intent(in) :: gprimd(3,3),phnons(2,nfftot**(1-1/nsym),(nspden/nsppol)-3*(nspden/4)),rprimd(3,3)
 real(dp),intent(inout) :: rhor(cplex*nfft,nspden)
 real(dp),intent(out) :: rhog(2,nfft)

!Local variables-------------------------------
!scalars
 integer :: ier,imagn,ind,ind2,indsy,ispden,isym,iup,izone,izone_max,j,j1,j2,j3,jsym
 integer :: k1,k2,k3,l1,l2,l3,me_fft
 integer :: n1,n2,n3,nd2,nproc_fft,nspden_eff,nsym_used,numpt,nup
 integer :: r2,rep,spaceComm
 logical,parameter :: afm_noncoll=.true.  ! TRUE if antiferro symmetries are used in non-collinear magnetism
 real(dp) :: magxsu1,magxsu2,magysu1,magysu2,magzsu1,magzsu2,mxi,mxr,myi,myr,mzi,mzr,phi,phr,rhosu1,rhosu2
 !character(len=500) :: message
!arrays
 integer,allocatable :: isymg(:)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: magngx(:,:),magngy(:,:),magngz(:,:)
 real(dp),allocatable :: rhosu1_arr(:),rhosu2_arr(:),work(:)
 real(dp),allocatable :: symafm_used(:),symrec_cart(:,:,:),symrel_cart(:,:,:)

!*************************************************************************
!
!Note the timing channel 17 excludes the different Fourier transforms

 ABI_ALLOCATE(work,(cplex*nfft))

!Special treatment for spin-polarized case
 if(nspden==2 .and. nsppol==2) then
!  When nspden=2 and nsppol=2, put total density in first half
!  of rhor array and spin up in second half  (up,down) => (up+down,up)
   call timab(17,1,tsec)
   work(:)=rhor(:,1)               ! up => work
   rhor(:,1)=rhor(:,1)+rhor(:,2)   ! up+down
   rhor(:,2)=work(:)               ! work => up
   call timab(17,2,tsec)
 end if

!Special treatment for antiferromagnetism case
 if(nspden==2 .and. nsppol==1) then
   call timab(17,1,tsec)
!  When nspden=2 and nsppol=1, (2*up) => (2*up,up)
!  Indeed, what was delivered to the present routine is a "total" density,
!  obtained from occupation numbers varying between 0 and 2,
!  but for spin up only potential.
   rhor(:,2)=half*rhor(:,1)
   call timab(17,2,tsec)
 end if

!Special treatment for non-collinear magnetism case
 if(nspden==4) then
   call timab(17,1,tsec)
!FR the half factors missing are recovered in dfpt_mkvxc_noncoll and dfpt_accrho
   rhor(:,1)=rhor(:,1)+rhor(:,4)     !nup+ndown
   rhor(:,2)=rhor(:,2)-rhor(:,1)     !mx (n+mx-n)
   rhor(:,3)=rhor(:,3)-rhor(:,1)     !my (n+my-n)
   rhor(:,4)=rhor(:,1)-two*rhor(:,4) !mz=n-2ndown
   call timab(17,2,tsec)
 end if


 if(nsym==1)then

   if(nspden==2 .and. nsppol==1) then ! There must be at least one anti-ferromagnetic operation
     MSG_BUG('In the antiferromagnetic case, nsym cannot be 1')
   end if

!  Blanchet Add free electron gas contribution
   if(associated(hightemp)) then
     if(hightemp%ioptden==1) then
       rhor(:,:)=rhor(:,:)+hightemp%nfreeel/hightemp%ucvol/nspden
     end if
   end if

!  If not using symmetry, still want total density in G space rho(G).
!  Fourier transform (incl normalization) to get rho(G)
   work(:)=rhor(:,1)
   call fourdp(cplex,rhog,work,-1,mpi_enreg,nfft,1,ngfft,0)
 else

!  Treat either full density, spin-up density or magnetization
!  Note the decrease of ispden to the value 1, in order to finish
!  with rhog of the total density (and not the spin-up density or magnetization)
   nspden_eff=nspden;if (nspden==4) nspden_eff=1
   do ispden=nspden_eff,1,-1

!    Prepare the density to be symmetrized, in the reciprocal space
     if(nspden==1 .or. nsppol==2 .or. (nspden==4.and.(.not.afm_noncoll)))then
       imagn=1
       nsym_used=0
       do isym=1,nsym
         if(symafm(isym)==1)nsym_used=nsym_used+1
!        DEBUG
!        write(std_out,*)' symrhg : isym,symafm(isym)',isym,symafm(isym)
!        ENDDEBUG
       end do
     else if(nspden==2 .and. nsppol==1)then   ! antiferromagnetic case
       imagn=ispden
       nsym_used=nsym/ispden
     else if (nspden==4) then
       imagn=1
       nsym_used=nsym/ispden
     end if

!    write(std_out,*)' symrhg : nsym_used=',nsym_used

!    rhor -fft-> rhog    (rhog is used as work space)
!    Note : it should be possible to reuse rhog in the antiferromagnetic case this would avoid one FFT
     work(:)=rhor(:,ispden)
     call fourdp(cplex,rhog,work,-1,mpi_enreg,nfft,1,ngfft,0)
     if (nspden==4) then
       ABI_ALLOCATE(magngx,(2,nfft))
       ABI_ALLOCATE(magngy,(2,nfft))
       ABI_ALLOCATE(magngz,(2,nfft))
       work(:)=rhor(:,2)
       call fourdp(cplex,magngx,work,-1,mpi_enreg,nfft,1,ngfft,0)
       work(:)=rhor(:,3)
       call fourdp(cplex,magngy,work,-1,mpi_enreg,nfft,1,ngfft,0)
       work(:)=rhor(:,4)
       call fourdp(cplex,magngz,work,-1,mpi_enreg,nfft,1,ngfft,0)
     end if

!    Begins the timing here only , to exclude FFTs
     call timab(17,1,tsec)

     n1=ngfft(1);n2=ngfft(2);n3=ngfft(3);nproc_fft=ngfft(10);me_fft=ngfft(11);nd2=n2/nproc_fft

!    Get the distrib associated with this fft_grid
     call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!    The following is only valid for total, up or dn density
!    -------------------------------------------------------

!    Get maxvalue of izone
     izone_max=count(irrzon(:,2,imagn)>0)
     ABI_ALLOCATE(rhosu1_arr,(izone_max))
     ABI_ALLOCATE(rhosu2_arr,(izone_max))

     numpt=0
     do izone=1,nfftot

!      Get repetition number
       rep=irrzon(izone,2,imagn)
       if(rep==0)exit

!      Compute number of unique points in this symm class:
       nup=nsym_used/rep

!      Accumulate charge over equivalent points
       rhosu1=zero
       rhosu2=zero
       do iup=1,nup
         ind=irrzon(iup+numpt,1,imagn)
         j=ind-1;j1=modulo(j,n1);j2=modulo(j/n1,n2);j3=j/(n1*n2);
         if(fftn2_distrib(j2+1)==me_fft)  then ! this ind is to be treated by me_fft
           r2=ffti2_local(j2+1) - 1
           ind=n1*(nd2*j3+r2)+j1+1 !this is ind in the current proc
           rhosu1=rhosu1+rhog(1,ind)*phnons(1,iup+numpt,imagn)&
&           -rhog(2,ind)*phnons(2,iup+numpt,imagn)
           rhosu2=rhosu2+rhog(2,ind)*phnons(1,iup+numpt,imagn)&
&           +rhog(1,ind)*phnons(2,iup+numpt,imagn)
         end if

       end do
       rhosu1=rhosu1/dble(nup)
       rhosu2=rhosu2/dble(nup)
       rhosu1_arr(izone)=rhosu1
       rhosu2_arr(izone)=rhosu2
!      Keep index of how many points have been considered:
       numpt=numpt+nup

     end do  ! End loop over izone

!    Reduction in case of FFT parallelization
     if(mpi_enreg%nproc_fft>1)then
       spaceComm=mpi_enreg%comm_fft
       call xmpi_sum(rhosu1_arr,spaceComm,ier)
       call xmpi_sum(rhosu2_arr,spaceComm,ier)
     end if

!    Now symmetrize the density
     numpt=0
     do izone=1,nfftot

!      Get repetition number
       rep=irrzon(izone,2,imagn)
       if(rep==0)exit

!      Compute number of unique points in this symm class:
       nup=nsym_used/rep

!      Define symmetrized rho(G) at equivalent points:
       do iup=1,nup
         ind=irrzon(iup+numpt,1,imagn)
!        decompose ind-1=n1(n2 j3+ j2)+j1
         j=ind-1;j1=modulo(j,n1);j2=modulo(j/n1,n2);j3=j/(n1*n2);
         if(fftn2_distrib(j2+1)==me_fft)  then ! this ind is to be treated by me_fft
           r2=ffti2_local(j2+1) - 1
!          ind in the proc ind-1=n1(nd2 j3+ r2)+j1
           ind=n1*(nd2*j3+r2)+j1+1 !this is ind in the current proc
           rhog(1,ind)=rhosu1_arr(izone)*phnons(1,iup+numpt,imagn)&
&           +rhosu2_arr(izone)*phnons(2,iup+numpt,imagn)
           rhog(2,ind)=rhosu2_arr(izone)*phnons(1,iup+numpt,imagn)&
&           -rhosu1_arr(izone)*phnons(2,iup+numpt,imagn)
         end if
       end do

!      Keep index of how many points have been considered:
       numpt=numpt+nup

     end do ! End loop over izone

     ABI_DEALLOCATE(rhosu1_arr)
     ABI_DEALLOCATE(rhosu2_arr)

!    The following is only valid for magnetization
!    ---------------------------------------------
     if (nspden==4) then

!      Transfer symmetries in cartesian coordinates
!      Compute symmetries in reciprocal space in cartesian coordinates
       ABI_ALLOCATE(symrec_cart,(3,3,nsym_used))
       ABI_ALLOCATE(symrel_cart,(3,3,nsym_used))
       ABI_ALLOCATE(symafm_used,(nsym_used))
       jsym=0
       do isym=1,nsym
         if (symafm(isym)/=1.and.(.not.afm_noncoll)) cycle
         jsym=jsym+1
         symafm_used(jsym)=dble(symafm(isym))
         call symredcart(rprimd,gprimd,symrel_cart(:,:,jsym),symrel(:,:,isym))
         call matr3inv(symrel_cart(:,:,jsym),symrec_cart(:,:,jsym))
       end do

       numpt=count(irrzon(:,1,imagn)>0)
       ABI_ALLOCATE(isymg,(numpt))
       isymg=0
       ABI_ALLOCATE(rhosu1_arr,(3*izone_max))
       ABI_ALLOCATE(rhosu2_arr,(3*izone_max))

!      Accumulate magnetization over equivalent points
!      Use all symmetries (not only those linking different g points)
!      Use Inverse[Transpose[symrel]]=symrec
       numpt=0
       do izone=1,izone_max
         magxsu1=zero;magxsu2=zero
         magysu1=zero;magysu2=zero
         magzsu1=zero;magzsu2=zero
         ind=irrzon(1+numpt,1,1)
         rep=irrzon(izone,2,1)
         nup=nsym_used/rep
         j=ind-1;l1=modulo(j,n1);l2=modulo(j/n1,n2);l3=j/(n1*n2)
         jsym=0
         do isym=1,nsym
           if (symafm(isym)/=1.and.(.not.afm_noncoll)) cycle
           jsym=jsym+1
           j1=symrel(1,1,isym)*l1+symrel(2,1,isym)*l2+symrel(3,1,isym)*l3
           j2=symrel(1,2,isym)*l1+symrel(2,2,isym)*l2+symrel(3,2,isym)*l3
           j3=symrel(1,3,isym)*l1+symrel(2,3,isym)*l2+symrel(3,3,isym)*l3
           k1=map_symrhg(j1,n1);k2=map_symrhg(j2,n2);k3=map_symrhg(j3,n3)
           indsy=1+k1+n1*(k2+n2*k3)
           ind2=-1;iup=numpt
           do while (ind2/=indsy.and.iup<numpt+nup)
             iup=iup+1;ind2=irrzon(iup,1,1)
           end do
           if (ind2/=indsy) then
             MSG_ERROR("ind2/=indsy in symrhg !")
           end if
           if (isymg(iup)==0) isymg(iup)=jsym
           if(fftn2_distrib(modulo((indsy-1)/n1,n2) + 1) == me_fft ) then  ! this is indsy is to be treated by me_fft
             indsy=n1*(nd2*k3+ ffti2_local(k2+1) -1)+k1+1        ! this is indsy in the current proc
             phr=phnons(1,iup,imagn);if (rep==1) phr=phr*symafm_used(jsym) !if rep==2, symafm is already included in phnons
             phi=phnons(2,iup,imagn);if (rep==1) phi=phi*symafm_used(jsym) !(see irrzg.F90)
             mxr=symrel_cart(1,1,jsym)*magngx(1,indsy)+symrel_cart(1,2,jsym)*magngy(1,indsy)+symrel_cart(1,3,jsym)*magngz(1,indsy)
             mxi=symrel_cart(1,1,jsym)*magngx(2,indsy)+symrel_cart(1,2,jsym)*magngy(2,indsy)+symrel_cart(1,3,jsym)*magngz(2,indsy)
             myr=symrel_cart(2,1,jsym)*magngx(1,indsy)+symrel_cart(2,2,jsym)*magngy(1,indsy)+symrel_cart(2,3,jsym)*magngz(1,indsy)
             myi=symrel_cart(2,1,jsym)*magngx(2,indsy)+symrel_cart(2,2,jsym)*magngy(2,indsy)+symrel_cart(2,3,jsym)*magngz(2,indsy)
             mzr=symrel_cart(3,1,jsym)*magngx(1,indsy)+symrel_cart(3,2,jsym)*magngy(1,indsy)+symrel_cart(3,3,jsym)*magngz(1,indsy)
             mzi=symrel_cart(3,1,jsym)*magngx(2,indsy)+symrel_cart(3,2,jsym)*magngy(2,indsy)+symrel_cart(3,3,jsym)*magngz(2,indsy)
             magxsu1=magxsu1+mxr*phr-mxi*phi;magxsu2=magxsu2+mxi*phr+mxr*phi
             magysu1=magysu1+myr*phr-myi*phi;magysu2=magysu2+myi*phr+myr*phi
             magzsu1=magzsu1+mzr*phr-mzi*phi;magzsu2=magzsu2+mzi*phr+mzr*phi
           end if
         end do
         rhosu1_arr(3*izone-2)=magxsu1/dble(nsym_used)
         rhosu1_arr(3*izone-1)=magysu1/dble(nsym_used)
         rhosu1_arr(3*izone  )=magzsu1/dble(nsym_used)
         rhosu2_arr(3*izone-2)=magxsu2/dble(nsym_used)
         rhosu2_arr(3*izone-1)=magysu2/dble(nsym_used)
         rhosu2_arr(3*izone  )=magzsu2/dble(nsym_used)
         numpt=numpt+nup
       end do

!      Reduction in case of FFT parallelization
       if(mpi_enreg%nproc_fft>1)then
         spaceComm=mpi_enreg%comm_fft
         call xmpi_sum(rhosu1_arr,spaceComm,ier)
         call xmpi_sum(rhosu2_arr,spaceComm,ier)
       end if

!      Now symmetrize the magnetization at equivalent points
!      Use Transpose[symrel]
       numpt=0
       do izone=1,izone_max
         rep=irrzon(izone,2,imagn)
         nup=nsym_used/rep
         do iup=1,nup
           ind=irrzon(iup+numpt,1,imagn)
           j=ind-1;j1=modulo(j,n1);j2=modulo(j/n1,n2);j3=j/(n1*n2)
           if(fftn2_distrib(j2+1)==me_fft)  then ! this ind is to be treated by me_fft
             r2=ffti2_local(j2+1) - 1
             ind=n1*(nd2*j3+r2)+j1+1  ! this is ind in the current proc
             jsym=isymg(iup+numpt)
             if (jsym==0) then
               MSG_ERROR("jsym=0 in symrhg !")
             end if
             magxsu1=rhosu1_arr(3*izone-2);magxsu2=rhosu2_arr(3*izone-2)
             magysu1=rhosu1_arr(3*izone-1);magysu2=rhosu2_arr(3*izone-1)
             magzsu1=rhosu1_arr(3*izone  );magzsu2=rhosu2_arr(3*izone  )
             phr=phnons(1,iup,imagn);if (rep==1) phr=phr*symafm_used(jsym) !if rep==2, symafm is already included in phnons
             phi=phnons(2,iup,imagn);if (rep==1) phi=phi*symafm_used(jsym) !(see irrzg.F90)
             mxr=symrec_cart(1,1,jsym)*magxsu1+symrec_cart(2,1,jsym)*magysu1+symrec_cart(3,1,jsym)*magzsu1
             mxi=symrec_cart(1,1,jsym)*magxsu2+symrec_cart(2,1,jsym)*magysu2+symrec_cart(3,1,jsym)*magzsu2
             myr=symrec_cart(1,2,jsym)*magxsu1+symrec_cart(2,2,jsym)*magysu1+symrec_cart(3,2,jsym)*magzsu1
             myi=symrec_cart(1,2,jsym)*magxsu2+symrec_cart(2,2,jsym)*magysu2+symrec_cart(3,2,jsym)*magzsu2
             mzr=symrec_cart(1,3,jsym)*magxsu1+symrec_cart(2,3,jsym)*magysu1+symrec_cart(3,3,jsym)*magzsu1
             mzi=symrec_cart(1,3,jsym)*magxsu2+symrec_cart(2,3,jsym)*magysu2+symrec_cart(3,3,jsym)*magzsu2
             magngx(1,ind)=mxr*phr+mxi*phi
             magngx(2,ind)=mxi*phr-mxr*phi
             magngy(1,ind)=myr*phr+myi*phi
             magngy(2,ind)=myi*phr-myr*phi
             magngz(1,ind)=mzr*phr+mzi*phi
             magngz(2,ind)=mzi*phr-mzr*phi
           end if
         end do
         numpt=numpt+nup
       end do
       ABI_DEALLOCATE(isymg)
       ABI_DEALLOCATE(rhosu1_arr)
       ABI_DEALLOCATE(rhosu2_arr)
       ABI_DEALLOCATE(symrec_cart)
       ABI_DEALLOCATE(symrel_cart)
       ABI_DEALLOCATE(symafm_used)

     end if ! nspden==4

!    Blanchet Add free electron gas contribution
     if(associated(hightemp)) then
       if(hightemp%ioptden==1) then
         rhog(1,1)=rhog(1,1)+hightemp%nfreeel/hightemp%ucvol/nspden
       end if
     end if

     call timab(17,2,tsec)

!    Pull out full or spin up density, now symmetrized
     call fourdp(cplex,rhog,work,1,mpi_enreg,nfft,1,ngfft,0)
     rhor(:,ispden)=work(:)
     if (nspden==4) then
       call fourdp(cplex,magngx,work,1,mpi_enreg,nfft,1,ngfft,0)
       rhor(:,2)=work(:)
       call fourdp(cplex,magngy,work,1,mpi_enreg,nfft,1,ngfft,0)
       rhor(:,3)=work(:)
       call fourdp(cplex,magngz,work,1,mpi_enreg,nfft,1,ngfft,0)
       rhor(:,4)=work(:)
       ABI_DEALLOCATE(magngx)
       ABI_DEALLOCATE(magngy)
       ABI_DEALLOCATE(magngz)
     end if

   end do ! ispden

 end if !  End on the condition nsym==1

 ABI_DEALLOCATE(work)

 contains

   function map_symrhg(j1,n1)

   integer :: map_symrhg
   integer,intent(in) :: j1,n1
   map_symrhg=mod(n1+mod(j1,n1),n1)
 end function map_symrhg

end subroutine symrhg
!!***

!!****f* m_spacepar/irrzg
!! NAME
!! irrzg
!!
!! FUNCTION
!! Find the irreducible zone in reciprocal space under the
!! symmetry group with real space rotations in symrel(3,3,nsym).
!! The (integer) rotation matrices symrel(3,3,nsym) express the new
!! real space positions (e.g. rotated atom positions) in REDUCED
!! coordinates, i.e. in coordinates expressed as fractions of real space
!! primitive translations (atomic coordinates xred).  tnons(3,nsym) express
!! the associated nonsymmorphic translations, again in reduced coordinates.
!! Special data structure created in irrzon.
!! First half holds mapping from irr zone to full zone;
!! part of second half holds repetition number info.
!! work1 is a work array to keep track of grid points found so far.
!! In case nspden=2 and nsppol=1, one has to take care of antiferromagnetic
!! operations. The subgroup of non-magnetic operations is used
!! to generate irrzon(:,:,2) and phnons(:,:,2), while the
!! full group is used to generate irrzon(:,:,1) and phnons(:,:,1)
!!
!! NOTES
!! for reference in the near future (2018), some notes: this routine should be duplicated for
!! magnetizations in spinorial formalism. The only main difference will
!! be that the density is not simply transported to the image point under symrel
!! but is a vector which has to be transformed as well
!! $S \vec{m} (\vec{x}) = \sigma \vec{m} (S \vec{x} + \tau)$
!! $S \vec{m} (\vec{G}) = \sigma exp(+ 2 \pi i \vec{G} \vec{tau}) \vec{m}(S^{-1 T} \vec{G})$
!! S is a symop, sigma the AFM sign flip if any, tau the partial non symmorphic translation
!! x a position, m the magnetization 3 vector
!!
!! The phase factor is the same as below for the density, but the collection of elements which
!! are equal is more complex: for a 3 or 6 axis the m vector could transform one component
!! to a linear combination of several others (I think). Things are not necessarily aligned
!! with the z axis.
!!
!! INPUTS
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of symmetry elements in group
!!  n1,n2,n3=box dimensions of real space grid (or fft box)
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrel(3,3,nsym)=symmetry matrices in real space (integers)
!!  tnons(3,nsym)=reduced nonsymmorphic translations
!! (symrel and tnons are in terms of real space primitive translations)
!!
!! OUTPUT
!!  irrzon(n1*n2*n3,2+(nspden/4),(nspden/nsppol)-3*(nspden/4))=integer array which contains the locations of related
!!   grid points and the repetition number for each symmetry class.
!!  phnons(2,n1*n2*n3,(nspden/nsppol)-3*(nspden/4))=phases associated with nonsymmorphic translations
!!
!! PARENTS
!!      m_ab7_kpoints,setsym,wfd_mkrho
!!
!! CHILDREN
!!      sort_int,wrtout
!!
!! SOURCE

subroutine irrzg(irrzon,nspden,nsppol,nsym,n1,n2,n3,phnons,symafm,symrel,tnons)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1,n2,n3,nspden,nsppol,nsym
!arrays
 integer,intent(in) :: symafm(nsym),symrel(3,3,nsym)
 integer,intent(out) :: irrzon(n1*n2*n3,2,(nspden/nsppol)-3*(nspden/4))
 real(dp),intent(in) :: tnons(3,nsym)
 real(dp),intent(out) :: phnons(2,n1*n2*n3,(nspden/nsppol)-3*(nspden/4))

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,id1,id2,id3,ifft,imagn,ind1,ind2,ipt,irep,isym,izone
 integer :: izonemax,j1,j2,j3,jj,k1,k2,k3,l1,l2,l3,nfftot,npt,nsym_used
 integer :: nzone,setzer,sppoldbl
 real(dp) :: arg,ph1i,ph1r,ph2i,ph2r,tau1,tau2,tau3
 logical,parameter :: afm_noncoll=.true. ! TRUE if antiferro symmetries are used in non-collinear magnetism
 character(len=500) :: message
!arrays
 integer,allocatable :: class(:),iperm(:),symafm_used(:),symrel_used(:,:,:)
 integer,allocatable :: work1(:)
 real(dp),allocatable :: tnons_used(:,:),work2(:,:)

! *************************************************************************

 ABI_ALLOCATE(class,(nsym))
 ABI_ALLOCATE(iperm,(nsym))
 ABI_ALLOCATE(work1,(n1*n2*n3))
 ABI_ALLOCATE(work2,(2,n1*n2*n3))

 nfftot=n1*n2*n3

 id1=n1/2+2
 id2=n2/2+2
 id3=n3/2+2

 sppoldbl=nspden/nsppol;if (nspden==4) sppoldbl=1

 do imagn=1,sppoldbl

!  Treat in a similar way the case of the full group and the non-magnetic subgroup
   nsym_used=0
   do isym=1,nsym
     if( (imagn==1 .and. sppoldbl==2) .or. symafm(isym)==1 .or. &
&     ((nspden==4).and.afm_noncoll) )then
       nsym_used=nsym_used+1
     end if
   end do

   if(imagn==2 .and. nsym_used/=nsym/2)then
     write(message, '(a,a,a,a,a,i4,a,i0)' )&
&     '  The number of ferromagnetic symmetry operations must be',ch10,&
&     '  half the total number of operations, while it is observed that',ch10,&
&     '  nsym=',nsym,' and nsym_magn=',nsym_used
     MSG_BUG(message)
   end if

   ABI_ALLOCATE(symafm_used,(nsym_used))
   ABI_ALLOCATE(symrel_used,(3,3,nsym_used))
   ABI_ALLOCATE(tnons_used,(3,nsym_used))

   nsym_used=0
   do isym=1,nsym
     if( (imagn==1 .and. sppoldbl==2) .or. symafm(isym)==1 .or.  &
&     ((nspden==4).and.afm_noncoll) ) then
       nsym_used=nsym_used+1
       symrel_used(:,:,nsym_used)=symrel(:,:,isym)
       tnons_used(:,nsym_used)=tnons(:,isym)
       symafm_used(nsym_used)=symafm(isym)
     end if
   end do
   if ((nspden/=4).or.(.not.afm_noncoll)) symafm_used=1


!  Zero out work array--later on, a zero entry will mean that
!  a given grid point has not yet been assigned to an ibz point
   work1(1:nfftot)=0
   irrzon(:,2,imagn)=0

!  Initialize at G=0 (always in irreducible zone)
   nzone=1
   irrzon(1,1,imagn)=1
   irrzon(1,2,imagn)=nsym_used
!  Set phase exp(2*Pi*I*G dot tnons) for G=0
   phnons(1,1,imagn)=one
   phnons(2,1,imagn)=zero
   npt=1
!  setting work1(1)=1 indicates that first grid point (G=0) is
!  in the iz (irreducible zone)
   work1(1)=1

   ind1=0

!  Loop over reciprocal space grid points:
   do i3=1,n3
     do i2=1,n2
       do i1=1,n1

         ind1=ind1+1

!        Check to see whether present grid point is equivalent to
!        any previously identified ibz point--if not, a new ibz point
!        has been found

         if (work1(ind1)==0) then

!          A new class has been found.

!          Get location of G vector (grid point) centered at 0 0 0
           l3=i3-(i3/id3)*n3-1
           l2=i2-(i2/id2)*n2-1
           l1=i1-(i1/id1)*n1-1

           do isym=1,nsym_used

!            Get rotated G vector Gj for each symmetry element
!            -- here we use the TRANSPOSE of symrel_used; assuming symrel_used expresses
!            the rotation in real space, the transpose is then appropriate
!            for G space symmetrization (p. 1172d,e of notes, 2 June 1995).
             j1=symrel_used(1,1,isym)*l1+&
&             symrel_used(2,1,isym)*l2+symrel_used(3,1,isym)*l3
             j2=symrel_used(1,2,isym)*l1+&
&             symrel_used(2,2,isym)*l2+symrel_used(3,2,isym)*l3
             j3=symrel_used(1,3,isym)*l1+&
&             symrel_used(2,3,isym)*l2+symrel_used(3,3,isym)*l3

!            Map into [0,n-1] and then add 1 for array index in [1,n]
             k1=1+mod(n1+mod(j1,n1),n1)
             k2=1+mod(n2+mod(j2,n2),n2)
             k3=1+mod(n3+mod(j3,n3),n3)
!            k1=1+map(j1,n1)
!            k2=1+map(j2,n2)
!            k3=1+map(j3,n3)

!            Get linear index of rotated point Gj
             ind2=k1+n1*((k2-1)+n2*(k3-1))

!            Store info for new class:
             class(isym)=ind2
             iperm(isym)=isym

!            Setting work array element to 1 indicates grid point has been
!            identified with iz point
             work1(ind2)=1

!            End of loop on isym
           end do

!          Sort integers into ascending order in each class
!          (this lumps together G vectors with the same linear index, i.e.
!          groups together symmetries which land on the same Gj)
           call sort_int(nsym_used,class,iperm)

!          Check repetition factor (how many identical copies of Gj occur
!          from all symmetries applied to G)
           irep=0
           do isym=1,nsym_used
             if (class(isym)==class(1)) then
               irep=irep+1
             end if
           end do
           ipt=nsym_used/irep

!          Repetition factor must be divisor of nsym_used:
           if (nsym_used/=(ipt*irep)) then
             write(message, '(a,i5,a,i6,a,a,a,a,a,a)' )&
&             '  irep=',irep,' not a divisor of nsym_used=',nsym_used,ch10,&
&             ' This usually indicates that',&
&             ' the input symmetries do not form a group.',ch10,&
&             ' Action : check the input symmetries carefully do they',&
&             ' form a group ? If they do, there is a code bug.'
             MSG_ERROR(message)
           end if

!          Compute phases for any nonsymmorphic symmetries
!          exp(-2*Pi*I*G dot tau(j)) for each symmetry j with
!          (possibly zero) nonsymmorphic translation tau(j)
           do jj=1,nsym_used
!            First get nonsymmorphic translation and see if nonzero
!            (iperm grabs the symmetries in the new order after sorting)
             isym=iperm(jj)
             tau1=tnons_used(1,isym)
             tau2=tnons_used(2,isym)
             tau3=tnons_used(3,isym)
             if (abs(tau1)>tol12.or.abs(tau2)>tol12&
&             .or.abs(tau3)>tol12) then
!              compute exp(-2*Pi*I*G dot tau) using original G
               arg=two_pi*(dble(l1)*tau1+dble(l2)*tau2+dble(l3)*tau3)
               work2(1,jj)=cos(arg)
               work2(2,jj)=-sin(arg)
             else
               work2(1,jj)=one
               work2(2,jj)=zero
             end if
           end do

!          All phases arising from symmetries which map to the same
!          G vector must actually be the same because
!          rho(Strans*G)=exp(2*Pi*I*(G) dot tau_S) rho(G)
!          must be satisfied; if exp(2*Pi*I*(G) dot tau_S) can be different
!          for two different symmetries S which both take G to the same St*G,
!          then the related Fourier components rho(St*G) must VANISH.
!          Thus: set "phase" to ZERO here in that case.
!          The G mappings occur in sets of irep members; if irep=1 then
!          all the G are unique.
!          MT 090212:
!          In the case of antiferromagn. symetries, one can have
!          rho(Strans*G)= -exp(2*Pi*I*(G) dot tau_S) rho(G)
!          (look at the minus !)
!          A special treatment is then operated on phons.
!          The later must be consistent with the use of phnons array
!          in symrhg.F90 routine.
!          XG 001108 :
!          Note that there is a tolerance on the
!          accuracy of tnons, especially when they are found from
!          the symmetry finder (with xred that might be a bit inaccurate)
           if (irep > 1) then
             do jj=1,nsym_used,irep
               setzer=0
               ph1r=work2(1,jj);ph1i=work2(2,jj)
               do j1=jj,jj+irep-1
                 ph2r=work2(1,j1);ph2i=work2(2,j1)
                 if (((ph2r+ph1r)**2+(ph2i+ph1i)**2) <= tol14) then
                   if (setzer/=1) setzer=-1
                 else if (((ph2r-ph1r)**2+(ph2i-ph1i)**2) > tol14) then
                   setzer=1
                 end if
               end do
!              Setzer= 0: phnons are all equal
!              Setzer=-1: phnons are equal in absolute value
!              Setzer= 1: some phnons are different
               if (setzer/=0) then
                 if (setzer==-1) then
                   if (afm_noncoll.and.nspden==4) then
                     arg=symafm_used(iperm(jj))
                     if (all(symafm_used(iperm(jj:jj+irep-1))==arg)) then
                       setzer=1
                     else
                       do j1=jj,jj+irep-1
                         work2(:,j1)=work2(:,j1)*dble(symafm_used(iperm(j1)))
                       end do
                     end if
                   else
                     setzer=1
                   end if
                 end if
                 if (setzer==1) work2(:,jj:jj+irep-1)=zero
               end if
             end do
!            Compress data if irep>1:
             jj=0
             do isym=1,nsym_used,irep
               jj=jj+1
               class(jj)=class(isym)
               work2(1,jj)=work2(1,isym)
               work2(2,jj)=work2(2,isym)
             end do
           end if

!          Put new unique points into irrzon array:
           irrzon(1+npt:ipt+npt,1,imagn)=class(1:ipt)

!          Put repetition number into irrzon array:
           irrzon(1+nzone,2,imagn)=irep

!          DEBUG
!          write(std_out,'(a,6i7)' )' irrzg : izone,i1,i2,i3,imagn,irrzon(859,2,1)=',&
!          &      1+nzone,i1,i2,i3,imagn,irrzon(859,2,1)
!          ENDDEBUG

!          Put phases (or 0) in phnons array:
           phnons(:,1+npt:ipt+npt,imagn)=work2(:,1:ipt)

!          Update number of points in irrzon array:
!          (irep must divide evenly into nsym_used !)
           npt=npt+ipt

!          Update number of classes:
           nzone=nzone+1

         end if
!
!        End of loop on reciprocal space points, with indices i1, i2, i3
       end do
     end do
   end do

   if (allocated(symafm_used))  then
     ABI_DEALLOCATE(symafm_used)
   end if
   if (allocated(symrel_used))  then
     ABI_DEALLOCATE(symrel_used)
   end if
   if (allocated(tnons_used))  then
     ABI_DEALLOCATE(tnons_used)
   end if

 end do ! imagn

!Make sure number of real space points accounted for equals actual number of grid points
 if (npt/=n1*n2*n3) then
   write(message, '(a,a,a,a,i10,a,i10,a,a,a,a,a,a,a,a,a)' ) ch10,&
&   ' irrzg : ERROR -',ch10,&
&   '  npt=',npt,' and n1*n2*n3=',n1*n2*n3,' are not equal',ch10,&
&   '  This says that the total of all points in the irreducible',&
&   '  sector in real space',ch10,&
&   '  and all symmetrically equivalent',&
&   '  points, npt, does not equal the actual number',ch10,&
&   '  of real space grid points.'
   call wrtout(std_out,message,'COLL')
   write(message,'(3a)') &
&   ' This may mean that the input symmetries do not form a group',ch10,&
&   ' Action : check input symmetries carefully for errors.'
   MSG_ERROR(message)
 end if

!Perform some checks
 do imagn=1,sppoldbl

   do ifft=1,nfftot
     if (irrzon(ifft,1,imagn)<1.or.irrzon(ifft,1,imagn)>nfftot) then
       write(message,'(a,4i0,a,a)')&
&       '  ifft,irrzon(ifft,1,imagn),nfftot,imagn=',ifft,irrzon(ifft,1,imagn),nfftot,imagn,ch10,&
&       '  =>irrzon goes outside acceptable bounds.'
       MSG_BUG(message)
     end if
   end do

   izonemax=0
   do izone=1,nfftot
!    Global bounds
     if (irrzon(izone,2,imagn)<0.or.irrzon(izone,2,imagn)>(nsym/imagn)) then
       write(message, '(a,5i7,a,a)' )&
&       ' izone,nzone,irrzon(izone,2,imagn),nsym,imagn =',izone,nzone,irrzon(izone,2,imagn),nsym,imagn,ch10,&
&       '  =>irrzon goes outside acceptable bounds.'
       MSG_BUG(message)
     end if
!    Second index only goes up to nzone
     if(izonemax==0)then
       if (irrzon(izone,2,imagn)==0)izonemax=izone-1
     end if
     if(izonemax/=0)then
       if (irrzon(izone,2,imagn)/=0) then
         message = ' beyond izonemax, irrzon(izone,2,imagn) should be zero'
         MSG_BUG(message)
       end if
     end if
   end do

 end do ! imagn

 ABI_DEALLOCATE(class)
 ABI_DEALLOCATE(iperm)
 ABI_DEALLOCATE(work1)
 ABI_DEALLOCATE(work2)

end subroutine irrzg
!!***

!!****f* m_spacepar/rotate_rho
!! NAME
!! rotate_rho
!!
!! FUNCTION
!! rotate density in real and reciprocal space
!!
!! INPUTS
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this proc) for the "fine" grid (see NOTES in respfn.F90)
!!  ngfft=array of dimensions for different FFT grids
!!  nspden=number of spin-density components
!!  rhor1(cplex*nfft,nspden)=array for Fourier transform of RF electron density
!!  === if psps%usepaw==1 TODO: extend to PAW
!!    pawrhoij1(natom) <type(pawrhoij_type)>= 1st-order paw rhoij occupancies and related data
!!  symrel1=single symmetry operation in real space to apply to rho
!!  tnon = eventual translation associated to symrel1
!!
!! OUTPUT
!!  rhog1_eq(2,nfft)= symmetric density in reciprocal space for equivalent perturbation
!!  rhor1_eq(cplex*nfft,nspden) = symmetric density in real space for equivalent perturbation
!!
!! PARENTS
!!      dfpt_looppert
!!
!! CHILDREN
!!      fourdp
!!
!! TODO
!!  This routine is deprecated. Now the symmetrization of the **potentials** is done in the m_dvdb
!!  Moreover one should take linear combinations of symmetrical components
!!
!! SOURCE

subroutine rotate_rho(cplex, itirev, mpi_enreg, nfft, ngfft, nspden, &
&   rhor1, rhog1_eq, rhor1_eq, symrel1, tnon)

!args
 integer,intent(in) :: cplex, nfft, nspden, itirev
 integer,intent(in) :: ngfft(18)

 integer, intent(in) :: symrel1(3,3)
 real(dp),intent(in) :: tnon(3)
 real(dp),intent(inout) :: rhor1(cplex*nfft,nspden)

 real(dp),intent(out) :: rhog1_eq(2,nfft)
 real(dp),intent(out) :: rhor1_eq(cplex*nfft,nspden)

 type(MPI_type),intent(in) :: mpi_enreg

! local vars
 integer :: id1, id2, id3
 integer :: n1, n2, n3, nd2
 integer :: l1, l2, l3
 integer :: i1, i2, i3, ind1, ind2
 integer :: j1, j2, j3
 integer :: k1, k2, k3
 integer :: nproc_fft, ispden, me_fft
 real(dp) :: arg
 logical :: t_tnon_nonzero

 real(dp) :: phnon1(2)
 real(dp), allocatable :: workg(:,:), workg_eq(:,:)
 character(len=500) :: message

! *************************************************************************

 n1=ngfft(1);n2=ngfft(2);n3=ngfft(3);nproc_fft=ngfft(10);me_fft=ngfft(11);nd2=n2/nproc_fft

 id1=n1/2+2
 id2=n2/2+2
 id3=n3/2+2

 rhog1_eq = zero
 rhor1_eq = zero

 if (itirev == 2) then
   write (message,'(3a,9I4,1a)') 'using time reversal. ',ch10,'Symrel1 = ', symrel1, ch10
 else
   write (message,'(3a,9I4,1a)') 'no time reversal. ',ch10,'Symrel1 = ', symrel1, ch10
 end if
 !call wrtout(std_out, message, 'COLL')

 t_tnon_nonzero = (any(abs(tnon) > tol12))

! eventually, for FFT parallelization
! call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

 ABI_ALLOCATE(workg,(2,nfft))
 ABI_ALLOCATE(workg_eq,(2,nfft))
 do ispden = 1, nspden

! fft input rhor1 to reciprocal space: uses work* as a buffer
   call fourdp(cplex,workg,rhor1(:,ispden),-1,mpi_enreg,nfft,1,ngfft,0)

! below taken from irrzg and setsym
!  Loop over reciprocal space grid points:
!  loop over local points in workg, and get back transform from rhog1,
!  which is presumed complete on each proc
   ind1=0
   do i3=1,n3
     do i2=1,n2
!       if(fftn2_distrib(i2)/=me_fft)  cycle ! this ind is not to be treated by me_fft
       do i1=1,n1

         ind1=ind1+1
!       r2=ffti2_local(i2+1) - 1
!       ind=n1*(nd2*i3+r2)+i1+1 !this is ind in the current proc

!      Get location of G vector (grid point) centered at 0 0 0
         l1=i1-(i1/id1)*n1-1
         l2=i2-(i2/id2)*n2-1
         l3=i3-(i3/id3)*n3-1

!      Get rotated G vector Gj for each symmetry element
!      -- here we use the TRANSPOSE of symrel1; assuming symrel1 expresses
!      the rotation in real space, the transpose is then appropriate
!      for G space symmetrization (p. 1172d,e of notes, 2 June 1995).
         j1=symrel1(1,1)*l1+symrel1(2,1)*l2+symrel1(3,1)*l3
         j2=symrel1(1,2)*l1+symrel1(2,2)*l2+symrel1(3,2)*l3
         j3=symrel1(1,3)*l1+symrel1(2,3)*l2+symrel1(3,3)*l3

!      Map into [0,n-1] and then add 1 for array index in [1,n]
         k1=1+mod(n1+mod(j1,n1),n1)
         k2=1+mod(n2+mod(j2,n2),n2)
         k3=1+mod(n3+mod(j3,n3),n3)

!      Get linear index of rotated point Gj
         ind2=k1+n1*((k2-1)+n2*(k3-1))
!       r2=ffti2_local(j2+1) - 1
!       ind=n1*(nd2*j3+r2)+j1+1 !this is ind may be in another proc!!

         phnon1(1) = one
         phnon1(2) = zero
         if (t_tnon_nonzero) then
!        compute exp(-2*Pi*I*G dot tau) using original G
! NB: this phase is same as that in irrzg and phnons1, and corresponds to complex conjugate of phase from G to Gj;
! we use it immediately below, to go _to_ workg(ind1)
! TODO : replace this with complex powers of exp(2pi tnon(1)) etc...
           arg=two_pi*(dble(l1)*tnon(1)+dble(l2)*tnon(2)+dble(l3)*tnon(3))
           phnon1(1) = cos(arg)
           phnon1(2) =-sin(arg)
         end if

!      rho(Strans*G)=exp(2*Pi*I*(G) dot tau_S) rho(G)
         workg_eq (1, ind1) = phnon1(1) * workg(1, ind2) &
&         - phnon1(2) * workg(2, ind2)
         workg_eq (2, ind1) = phnon1(1) * workg(2, ind2) &
&         + phnon1(2) * workg(1, ind2)

       end do
     end do
   end do

! accumulate rhog1_eq
   if (ispden == 1) rhog1_eq = workg_eq

! FFT back to real space to get rhor1_eq
!    Pull out full or spin up density, now symmetrized
   call fourdp(cplex,workg_eq,rhor1_eq(:,ispden),1,mpi_enreg,nfft,1,ngfft,0)

 end do !nspden

 ABI_DEALLOCATE(workg)
 ABI_DEALLOCATE(workg_eq)

end subroutine rotate_rho
!!***

!!****f* m_spacepar/setsym
!! NAME
!! setsym
!!
!! FUNCTION
!! Set up irreducible zone in  G space by direct calculation.
!! Do not call this routine if nsym=1 (only identity symmetry).
!! Only indsym and symrec get returned if iscf=0.
!! symrec needed to symmetrize coordinate gradients in sygrad.
!! (symrec is redundant and could be removed later in favor of symrel)
!!
!! INPUTS
!! iscf=(<= 0  =>non-SCF), >0 => SCF
!! natom=number of atoms in unit cell
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! nspden=number of spin-density components
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! nsym=number of symmetries in space group (at least 1)
!! symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym)=symmetry operations in terms of real space primitive translations
!! tnons(3,nsym)=nonsymmorphic translations of space group in terms
!! of real space primitive translations (may be 0)
!! typat(natom)=atom type (integer) for each atom
!! xred(3,natom)=atomic coordinates in terms of real space primitive translations
!!
!! OUTPUT
!! indsym(4,nsym,natom)=indirect indexing of atom labels--see subroutine
!!   symatm for definition (if nsym>1)
!! irrzon(nfft,2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!! phnons(2,nfft,(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!! symrec(3,3,nsym)=symmetry operations in terms of reciprocal
!!   space primitive translations (if nsym>1)
!!
!! NOTES
!! nsppol and nspden are needed in case of (anti)ferromagnetic symmetry operations
!!
!! PARENTS
!!      dfpt_looppert,gstate,m_dvdb,nonlinear,respfn,scfcv
!!
!! CHILDREN
!!      chkgrp,irrzg,mati3inv,symatm,symdet,timab
!!
!! SOURCE

subroutine setsym(indsym,irrzon,iscf,natom,nfft,ngfft,nspden,nsppol,nsym,phnons,&
& symafm,symrec,symrel,tnons,typat,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iscf,natom,nfft,nspden,nsppol,nsym
!arrays
 integer,intent(in) :: ngfft(18),symafm(nsym),symrel(3,3,nsym),typat(natom)
 integer,intent(out) :: indsym(4,nsym,natom)
 integer,intent(inout) :: irrzon(nfft,2,(nspden/nsppol)-3*(nspden/4)) !vz_i
 integer,intent(out) :: symrec(3,3,nsym)
 real(dp),intent(in) :: tnons(3,nsym),xred(3,natom)
 real(dp),intent(out) :: phnons(2,nfft,(nspden/nsppol)-3*(nspden/4))

!Local variables-------------------------------
!scalars
 integer :: isym,ierr
 real(dp) :: tolsym8
!arrays
 integer,allocatable :: determinant(:)
 real(dp) :: tsec(2)

! *************************************************************************

 call timab(6,1,tsec)

!Check that symmetries have unity determinant
 ABI_ALLOCATE(determinant,(nsym))
 call symdet(determinant,nsym,symrel)
 ABI_DEALLOCATE(determinant)


!Get the symmetry matrices in terms of reciprocal basis
 do isym=1,nsym
   call mati3inv(symrel(:,:,isym),symrec(:,:,isym))
 end do

!Check for group closure
 call chkgrp(nsym,symafm,symrel,ierr)
 ABI_CHECK(ierr==0,"Error in group closure")

 call chkgrp(nsym,symafm,symrec,ierr)
 ABI_CHECK(ierr==0,"Error in group closure")

!Obtain a list of rotated atom labels:
 tolsym8=tol8
 call symatm(indsym,natom,nsym,symrec,tnons,tolsym8,typat,xred)

!If non-SCF calculation, or nsym==1, do not need IBZ data
 if ( (iscf>0 .or. iscf==-3) .and. nsym>1 ) then
!  Locate irreducible zone in reciprocal space for symmetrization:
   call irrzg(irrzon,nspden,nsppol,nsym,ngfft(1),ngfft(2),ngfft(3),phnons,symafm,symrel,tnons)
 end if

 call timab(6,2,tsec)

end subroutine setsym
!!***

end module m_spacepar
!!***
