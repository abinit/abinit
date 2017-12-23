!{\src2tex{textfont=tt}}
!!****f* ABINIT/mklocl_recipspace
!! NAME
!! mklocl_recipspace
!!
!! FUNCTION
!! Optionally compute :
!!  option=1 : local ionic potential throughout unit cell
!!  option=2 : contribution of local ionic potential to E gradient wrt xred
!!  option=3 : contribution of local ionic potential to stress tensor
!!  option=4 : contribution of local ionic potential to
!!                second derivative of E wrt xred
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  if(option==3) eei=local pseudopotential part of total energy (hartree)
!!  gmet(3,3)=reciprocal space metric ($\textrm{Bohr}^{-2}$).
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  gsqcut=cutoff on $|G|^2$: see setup1 for definition (doubled sphere).
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=information about MPI parallelization
!!  mqgrid=number of grid pts in q array for f(q) spline.
!!  natom=number of atoms in unit cell.
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  ntypat=number of types of atoms.
!!  option= (see above)
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information.
!!  qgrid(mqgrid)=q grid for spline from 0 to qmax.
!!  qprtrb(3)= integer wavevector of possible perturbing potential
!!   in basis of reciprocal lattice translations
!!  rhog(2,nfft)=electron density rho(G) (electrons/$\textrm{Bohr}^3$)
!!    (needed if option==2 or if option==4)
!!  ucvol=unit cell volume ($\textrm{Bohr}^3$).
!!  vlspl(mqgrid,2,ntypat)=q^2 v(q) spline for each type of atom.
!!  vprtrb(2)=complex amplitude of possible perturbing potential; if nonzero,
!!   perturbing potential is added of the form
!!   $V(G)=(vprtrb(1)+I*vprtrb(2))/2$ at the values G=qprtrb and
!!   $(vprtrb(1)-I*vprtrb(2))/2$ at $G=-qprtrb$ (integers)
!!
!! OUTPUT
!!  (if option==1) vpsp(nfft)=local crystal pseudopotential in real space.
!!  (if option==2) grtn(3,natom)=grads of Etot wrt tn.
!!  (if option==3) lpsstr(6)=components of local psp part of stress tensor
!!   (Cartesian coordinates, symmetric tensor) in hartree/$\textrm{bohr}^3$
!!   Store 6 unique components in order 11, 22, 33, 32, 31, 21
!!  (if option==4) dyfrlo(3,3,natom)=d2 Eei/d tn(i)/d tn(j).  (Hartrees)
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Note that the present routine is tightly connected to the dfpt_vlocal.f routine,
!! that compute the derivative of the local ionic potential
!! with respect to one atomic displacement. The argument list
!! and the internal loops to be considered were sufficiently different
!! as to make the two routine different.
!!
!! PARENTS
!!      dfpt_dyfro,mklocl,stress
!!
!! CHILDREN
!!      fourdp,ptabs_fourdp,timab,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mklocl_recipspace(dyfrlo,eei,gmet,gprimd,grtn,gsqcut,lpsstr,mgfft,&
&  mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,ntypat,option,paral_kgb,ph1d,qgrid,qprtrb,&
&  rhog,ucvol,vlspl,vprtrb,vpsp)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_mpinfo,   only : ptabs_fourdp

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mklocl_recipspace'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,mqgrid,natom,nfft,ntypat,option,paral_kgb
 real(dp),intent(in) :: eei,gsqcut,ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: nattyp(ntypat),ngfft(18),qprtrb(3)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: qgrid(mqgrid),rhog(2,nfft),vlspl(mqgrid,2,ntypat)
 real(dp),intent(in) :: vprtrb(2)
 real(dp),intent(out) :: dyfrlo(3,3,natom),grtn(3,natom),lpsstr(6) !vz_i
 real(dp),intent(inout) :: vpsp(nfft) !vz_i

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i3,ia,ia1,ia2,id1,id2,id3,ierr,ig1,ig2,ig3,ii,itypat
 integer :: jj,me_fft,me_g0,n1,n2,n3,nproc_fft,shift1
 integer :: shift2,shift3
 real(dp),parameter :: tolfix=1.0000001_dp
 real(dp) :: aa,bb,cc,cutoff,dbl_ig1,dbl_ig2,dbl_ig3,dd,diff,dq,dq2div6,dqdiv6
 real(dp) :: dqm1,ee,ff,gmag,gsquar,ph12i,ph12r,ph1i,ph1r,ph2i,ph2r
 real(dp) :: ph3i,ph3r,phimag_igia,phre_igia,sfi,sfr
 real(dp) :: svion,svioni,svionr,term,vion1,vion2,xnorm
 character(len=500) :: message
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: gcart(3),tsec(2)
 real(dp),allocatable :: work1(:,:)

! *************************************************************************

!Define G^2 based on G space metric gmet.
! gsq(i1,i2,i3)=dble(i1*i1)*gmet(1,1)+dble(i2*i2)*gmet(2,2)+&
!& dble(i3*i3)*gmet(3,3)+dble(2*i1*i2)*gmet(1,2)+&
!& dble(2*i2*i3)*gmet(2,3)+dble(2*i3*i1)*gmet(3,1)

!Real and imaginary parts of phase--statment functions:
! phr(x1,y1,x2,y2,x3,y3)=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
! phi(x1,y1,x2,y2,x3,y3)=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
! ph1(nri,i1,ia)=ph1d(nri,i1+1+n1+(ia-1)*(2*n1+1))
! ph2(nri,i2,ia)=ph1d(nri,i2+1+n2+(ia-1)*(2*n2+1)+&
!& natom*(2*n1+1))
! ph3(nri,i3,ia)=ph1d(nri,i3+1+n3+(ia-1)*(2*n3+1)+&
!& natom*(2*n1+1+2*n2+1))
! phre(i1,i2,i3,ia)=phr(ph1(re,i1,ia),ph1(im,i1,ia),ph2(re,i2,ia),&
!& ph2(im,i2,ia),ph3(re,i3,ia),ph3(im,i3,ia))
! phimag(i1,i2,i3,ia)=phi(ph1(re,i1,ia),ph1(im,i1,ia),ph2(re,i2,ia),&
!& ph2(im,i2,ia),ph3(re,i3,ia),ph3(im,i3,ia))

!-----

!Keep track of total time spent in mklocl
 if(option==2)then
   call timab(72,1,tsec)
 end if
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 me_fft=ngfft(11)
 nproc_fft=ngfft(10)

!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!Zero out array to permit accumulation over atom types below:
 if(option==1)then
   ABI_ALLOCATE(work1,(2,nfft))
   work1(:,:)=zero
 end if
!
 dq=(qgrid(mqgrid)-qgrid(1))/dble(mqgrid-1)
 dqm1=1.0_dp/dq
 dqdiv6=dq/6.0_dp
 dq2div6=dq**2/6.0_dp
 cutoff=gsqcut*tolfix
 id1=n1/2+2
 id2=n2/2+2
 id3=n3/2+2
 grtn(:,:)=zero
 lpsstr(:)=zero
 dyfrlo(:,:,:)=zero
 me_g0=0
 ia1=1

 do itypat=1,ntypat
!  ia1,ia2 sets range of loop over atoms:
   ia2=ia1+nattyp(itypat)-1

   ii=0
   do i3=1,n3
     ig3=i3-(i3/id3)*n3-1
     do i2=1,n2
       ig2=i2-(i2/id2)*n2-1
       if(fftn2_distrib(i2) == me_fft ) then 
         do i1=1,n1
           ig1=i1-(i1/id1)*n1-1

           ii=ii+1
!          ***     GET RID OF THIS THESE IF STATEMENTS (if they slow code)
!          Skip G=0:
!          if (ii==1) cycle
           if (ig1==0 .and. ig2==0 .and. ig3==0) me_g0=1
           if (ig1==0 .and. ig2==0 .and. ig3==0) cycle

           gsquar=gsq_mk(ig1,ig2,ig3)
!          Skip G**2 outside cutoff:
           if (gsquar<=cutoff) then
             gmag=sqrt(gsquar)

!            Compute vion(G) for given type of atom
             jj=1+int(gmag*dqm1)
             diff=gmag-qgrid(jj)

!            Evaluate spline fit from q^2 V(q) to get V(q):
!            (p. 86 Numerical Recipes, Press et al;
!            NOTE error in book for sign
!            of "aa" term in derivative; also see splfit routine).

             bb = diff*dqm1
             aa = 1.0_dp-bb
             cc = aa*(aa**2-1.0_dp)*dq2div6
             dd = bb*(bb**2-1.0_dp)*dq2div6

             vion1 = (aa*vlspl(jj,1,itypat)+bb*vlspl(jj+1,1,itypat) +&
&             cc*vlspl(jj,2,itypat)+dd*vlspl(jj+1,2,itypat) ) / gsquar

             if(option==1)then

!              Assemble structure factor over all atoms of given type:
               sfr=zero
               sfi=zero
               do ia=ia1,ia2
                 sfr=sfr+phre_mk(ig1,ig2,ig3,ia)
                 sfi=sfi-phimag_mk(ig1,ig2,ig3,ia)
               end do
!              Multiply structure factor times vion:
               work1(re,ii)=work1(re,ii)+sfr*vion1
               work1(im,ii)=work1(im,ii)+sfi*vion1

             else if(option==2 .or. option==4)then

!              Compute Re and Im part of (2Pi)*Vion(G)*rho(G):
               svionr=(two_pi*vion1)*rhog(re,ii)
               svioni=(two_pi*vion1)*rhog(im,ii)

!              Loop over all atoms of this type:
               do ia=ia1,ia2
                 shift1=1+n1+(ia-1)*(2*n1+1)
                 shift2=1+n2+(ia-1)*(2*n2+1)+natom*(2*n1+1)
                 shift3=1+n3+(ia-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1)
                 ph1r=ph1d(1,ig1+shift1)
                 ph1i=ph1d(2,ig1+shift1)
                 ph2r=ph1d(1,ig2+shift2)
                 ph2i=ph1d(2,ig2+shift2)
                 ph3r=ph1d(1,ig3+shift3)
                 ph3i=ph1d(2,ig3+shift3)
                 ph12r=ph1r*ph2r-ph1i*ph2i
                 ph12i=ph1r*ph2i+ph1i*ph2r
                 phre_igia=ph12r*ph3r-ph12i*ph3i
                 phimag_igia=ph12r*ph3i+ph12i*ph3r

                 if(option==2)then

!                  Compute "Vion" part of gradient
!                  svion=svioni*phre(ig1,ig2,ig3,ia)+svionr*phimag(ig1,ig2,ig3,ia)
                   svion=svioni*phre_igia+svionr*phimag_igia

!                  Open loop over 3-index for speed:
                   grtn(1,ia)=grtn(1,ia)-dble(ig1)*svion
                   grtn(2,ia)=grtn(2,ia)-dble(ig2)*svion
                   grtn(3,ia)=grtn(3,ia)-dble(ig3)*svion

                 else

!                  Compute "Vion" part of the second derivative
!                  svion=two_pi*
!                  (svionr*phre(ig1,ig2,ig3,ia)-svioni*phimag(ig1,ig2,ig3,ia))
                   svion=two_pi*(svionr*phre_igia-svioni*phimag_igia)

!                  Open loop over 3-index for speed
                   dbl_ig1=dble(ig1) ; dbl_ig2=dble(ig2) ; dbl_ig3=dble(ig3)
                   dyfrlo(1,1,ia)=dyfrlo(1,1,ia)-dbl_ig1*dbl_ig1*svion
                   dyfrlo(1,2,ia)=dyfrlo(1,2,ia)-dbl_ig1*dbl_ig2*svion
                   dyfrlo(1,3,ia)=dyfrlo(1,3,ia)-dbl_ig1*dbl_ig3*svion
                   dyfrlo(2,2,ia)=dyfrlo(2,2,ia)-dbl_ig2*dbl_ig2*svion
                   dyfrlo(2,3,ia)=dyfrlo(2,3,ia)-dbl_ig2*dbl_ig3*svion
                   dyfrlo(3,3,ia)=dyfrlo(3,3,ia)-dbl_ig3*dbl_ig3*svion

                 end if

               end do

             else if(option==3)then

!              Also get (dV(q)/dq)/q:
!              (note correction of Numerical Recipes sign error
!              before (3._dp*aa**2-1._dp)
!              ee*dqm1 + ff*dqdiv6 is the best estimate of dV(q)/dq from splines
               ee= vlspl(jj+1,1,itypat)-vlspl(jj,1,itypat)
               ff=  (3._dp*bb**2-1._dp)*vlspl(jj+1,2,itypat) &
&               - (3._dp*aa**2-1._dp)*vlspl(jj,2,itypat)
               vion2 = ( ( ee*dqm1 + ff*dqdiv6 )/gmag&
&               - 2.0_dp*vion1                 ) / gsquar

               gcart(1)=gprimd(1,1)*dble(ig1)+gprimd(1,2)*dble(ig2)+&
&               gprimd(1,3)*dble(ig3)
               gcart(2)=gprimd(2,1)*dble(ig1)+gprimd(2,2)*dble(ig2)+&
&               gprimd(2,3)*dble(ig3)
               gcart(3)=gprimd(3,1)*dble(ig1)+gprimd(3,2)*dble(ig2)+&
&               gprimd(3,3)*dble(ig3)
!              Assemble structure over all atoms of given type
               sfr=zero
               sfi=zero
               do ia=ia1,ia2
                 sfr=sfr+phre_mk(ig1,ig2,ig3,ia)
                 sfi=sfi-phimag_mk(ig1,ig2,ig3,ia)
               end do

!              Compute Re( rho^*(G)* sf ) * [(dV(G)/dG)/|G|]
               term=(rhog(re,ii)*sfr+rhog(im,ii)*sfi)*vion2

!              Compute contribution to stress tensor
               lpsstr(1)=lpsstr(1)-term*gcart(1)*gcart(1)
               lpsstr(2)=lpsstr(2)-term*gcart(2)*gcart(2)
               lpsstr(3)=lpsstr(3)-term*gcart(3)*gcart(3)
               lpsstr(4)=lpsstr(4)-term*gcart(3)*gcart(2)
               lpsstr(5)=lpsstr(5)-term*gcart(3)*gcart(1)
               lpsstr(6)=lpsstr(6)-term*gcart(2)*gcart(1)

             else
               write(message, '(a,i0,a)' )' mklocl: Option=',option,' not allowed.'
               MSG_BUG(message)
             end if ! End option choice

!            End skip G**2 outside cutoff:
           end if

!          End loop on n1, n2, n3. There is a "cycle" inside the loop
         end do
       end if ! this plane is for me_fft
     end do
   end do

!  Symmetrize the dynamical matrix with respect to indices
   do ia=ia1,ia2
     dyfrlo(2,1,ia)=dyfrlo(1,2,ia)
     dyfrlo(3,1,ia)=dyfrlo(1,3,ia)
     dyfrlo(3,2,ia)=dyfrlo(2,3,ia)
   end do

   ia1=ia2+1

!  End loop on type of atoms
 end do

 if(option==1)then
!  Dont't change work1 on g=0 if Poisson solver is used since work1
!  hold not the potential but the density generated by the pseudo.
   if(me_g0 == 1) then
!    Set Vloc(G=0)=0:
     work1(re,1)=zero
     work1(im,1)=zero
   end if

!  DEBUG
!  write(std_out,*) ' mklocl_recipspace : will add potential with strength vprtrb(:)=',vprtrb(:)
!  ENDDEBUG

!  Allow for the addition of a perturbing potential
   if ((vprtrb(1)**2+vprtrb(2)**2) > 1.d-30) then
!    Find the linear indices which correspond with the input
!    wavevector qprtrb
!    The double modulus handles both i>=n and i<0, mapping into [0,n-1];
!    then add 1 to get range [1,n] for each
     i3=1+mod(n3+mod(qprtrb(3),n3),n3)
     i2=1+mod(n2+mod(qprtrb(2),n2),n2)
     i1=1+mod(n1+mod(qprtrb(1),n1),n1)
!    Compute the linear index in the 3 dimensional array
     ii=i1+n1*((ffti2_local(i2)-1)+(n2/nproc_fft)*(i3-1))
!    Add in the perturbation at G=qprtrb
     work1(re,ii)=work1(re,ii)+0.5_dp*vprtrb(1)
     work1(im,ii)=work1(im,ii)+0.5_dp*vprtrb(2)
!    Same thing for G=-qprtrb
     i3=1+mod(n3+mod(-qprtrb(3),n3),n3)
     i2=1+mod(n2+mod(-qprtrb(2),n2),n2)
     i1=1+mod(n1+mod(-qprtrb(1),n1),n1)
!    ii=i1+n1*((i2-1)+n2*(i3-1))
     work1(re,ii)=work1(re,ii)+0.5_dp*vprtrb(1)
     work1(im,ii)=work1(im,ii)-0.5_dp*vprtrb(2)
     write(message, '(a,1p,2e12.4,a,0p,3i4,a)' )&
&     ' mklocl: perturbation of vprtrb=', vprtrb,&
&     ' and q=',qprtrb,' has been added'
     call wrtout(std_out,message,'COLL')
   end if

!  Transform back to real space
   call fourdp(1,work1,vpsp,1,mpi_enreg,nfft,ngfft,paral_kgb,0)

!  Divide by unit cell volume
   xnorm=1.0_dp/ucvol
   vpsp(:)=vpsp(:)*xnorm

   ABI_DEALLOCATE(work1)

 end if

 if(option==2)then
!  Init mpi_comm
   if(mpi_enreg%nproc_fft>1)then
     call timab(48,1,tsec)
     call xmpi_sum(grtn,mpi_enreg%comm_fft ,ierr)
     call timab(48,2,tsec)
   end if
   call timab(72,2,tsec)
 end if

 if(option==3)then
!  Init mpi_comm
   if(mpi_enreg%nproc_fft>1)then
     call timab(48,1,tsec)
     call xmpi_sum(lpsstr,mpi_enreg%comm_fft ,ierr)
     call timab(48,2,tsec)
   end if

!  Normalize and add term -eei/ucvol on diagonal
!  (see page 802 of notes)
   lpsstr(1)=(lpsstr(1)-eei)/ucvol
   lpsstr(2)=(lpsstr(2)-eei)/ucvol
   lpsstr(3)=(lpsstr(3)-eei)/ucvol
   lpsstr(4)=lpsstr(4)/ucvol
   lpsstr(5)=lpsstr(5)/ucvol
   lpsstr(6)=lpsstr(6)/ucvol

 end if

 if(option==4)then
!  Init mpi_comm
   if(mpi_enreg%nproc_fft>1)then
     call timab(48,1,tsec)
     call xmpi_sum(dyfrlo,mpi_enreg%comm_fft ,ierr)
     call timab(48,2,tsec)
   end if
 end if

 contains

!Real and imaginary parts of phase--statment functions:
   function phr_mk(x1,y1,x2,y2,x3,y3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phr_mk'
!End of the abilint section

   real(dp) :: phr_mk,x1,x2,x3,y1,y2,y3
   phr_mk=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
 end function phr_mk

   function phi_mk(x1,y1,x2,y2,x3,y3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phi_mk'
!End of the abilint section

   real(dp):: phi_mk,x1,x2,x3,y1,y2,y3
   phi_mk=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
 end function phi_mk

   function ph1_mk(nri,ig1,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ph1_mk'
!End of the abilint section

   real(dp):: ph1_mk
   integer :: nri,ig1,ia
   ph1_mk=ph1d(nri,ig1+1+n1+(ia-1)*(2*n1+1))
 end function ph1_mk

   function ph2_mk(nri,ig2,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ph2_mk'
!End of the abilint section

   real(dp):: ph2_mk
   integer :: nri,ig2,ia
   ph2_mk=ph1d(nri,ig2+1+n2+(ia-1)*(2*n2+1)+natom*(2*n1+1))
 end function ph2_mk

   function ph3_mk(nri,ig3,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ph3_mk'
!End of the abilint section

   real(dp):: ph3_mk
   integer :: nri,ig3,ia
   ph3_mk=ph1d(nri,ig3+1+n3+(ia-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1))
 end function ph3_mk

   function phre_mk(ig1,ig2,ig3,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phre_mk'
!End of the abilint section

   real(dp):: phre_mk
   integer :: ig1,ig2,ig3,ia
   phre_mk=phr_mk(ph1_mk(re,ig1,ia),ph1_mk(im,ig1,ia),&
&   ph2_mk(re,ig2,ia),ph2_mk(im,ig2,ia),ph3_mk(re,ig3,ia),ph3_mk(im,ig3,ia))
 end function phre_mk

   function phimag_mk(ig1,ig2,ig3,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phimag_mk'
!End of the abilint section

   real(dp) :: phimag_mk
   integer :: ig1,ig2,ig3,ia
   phimag_mk=phi_mk(ph1_mk(re,ig1,ia),ph1_mk(im,ig1,ia),&
&   ph2_mk(re,ig2,ia),ph2_mk(im,ig2,ia),ph3_mk(re,ig3,ia),ph3_mk(im,ig3,ia))
 end function phimag_mk

   function gsq_mk(i1,i2,i3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gsq_mk'
!End of the abilint section

   real(dp) :: gsq_mk
   integer :: i1,i2,i3
   gsq_mk=dble(i1*i1)*gmet(1,1)+dble(i2*i2)*gmet(2,2)+&
&   dble(i3*i3)*gmet(3,3)+dble(2*i1*i2)*gmet(1,2)+&
&   dble(2*i2*i3)*gmet(2,3)+dble(2*i3*i1)*gmet(3,1)
 end function gsq_mk

end subroutine mklocl_recipspace
!!***
