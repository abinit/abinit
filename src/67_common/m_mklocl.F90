!!****m* ABINIT/m_mklocl
!! NAME
!!  m_mklocl
!!
!! FUNCTION
!!   Routines related to the local part of the pseudopotentials.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, MM, DRH)
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

module m_mklocl

 use defs_basis
 use defs_wvltypes
 use m_abicore
 use m_errors
 use m_xmpi
 use m_dtset

 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes, only : MPI_type
 use m_time,     only : timab
 use m_geometry, only : xred2xcart
 use m_mpinfo,   only : ptabs_fourdp
 use m_pawtab,   only : pawtab_type
 use m_mklocl_realspace, only : mklocl_realspace, mklocl_wavelets
 use m_fft,      only : fourdp
 use m_gtermcutoff,only : termcutoff

 use m_splines,  only : splfit

#if defined HAVE_BIGDFT
 use BigDFT_API, only : ELECTRONIC_DENSITY
 use m_abi2big, only : wvl_rho_abi2big
#endif

 implicit none

 private
!!***

 public :: mklocl
 public :: mklocl_recipspace
 public :: dfpt_vlocal           ! Local part of 1st-order potential due to atomic displacement.
 public :: vlocalstr             ! Compute strain derivatives of local ionic potential
 public :: dfpt_vlocaldq         ! Compute the first q-gradient of the 1st-order potential due to atomic displacement. 
 public :: dfpt_vlocaldqdq       ! Compute the second q-gradient of the 1st-order potential due to atomic displacement. 
 public :: dfpt_vmetdqdq       ! Compute the second q-gradient of the 1st-order potential due to a metric perturbation.
!!***

contains
!!***

!!****f* ABINIT/mklocl
!! NAME
!! mklocl
!!
!! FUNCTION
!! This method is a wrapper for mklocl_recipspace and mklocl_realspace.
!! It does some consistency checks before calling one of the two methods.
!!
!! Optionally compute :
!!  option=1 : local ionic potential throughout unit cell
!!  option=2 : contribution of local ionic potential to E gradient wrt xred
!!  option=3 : contribution of local ionic potential to
!!                stress tensor (only with reciprocal space computations)
!!  option=4 : contribution of local ionic potential to
!!                second derivative of E wrt xred  (only with reciprocal space computations)
!!
!! INPUTS
!!  if(option==3) eei=local pseudopotential part of total energy (hartree)
!!  gmet(3,3)=reciprocal space metric ($\textrm{Bohr}^{-2}$).
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  gsqcut=cutoff on $|G|^2$: see setup1 for definition (doubled sphere).
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in unit cell.
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms.
!!  option= (see above)
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  qprtrb(3)= integer wavevector of possible perturbing potential
!!   in basis of reciprocal lattice translations
!!  rhog(2,nfft)=electron density rho(G) (electrons/$\textrm{Bohr}^3$)
!!    (needed if option==2 or if option==4)
!!  rhor(nfft,nspden)=electron density in electrons/bohr**3.
!!    (needed if option==2 or if option==4)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  ucvol=unit cell volume ($\textrm{Bohr}^3$).
!!  vprtrb(2)=complex amplitude of possible perturbing potential; if nonzero,
!!   perturbing potential is added of the form
!!   $V(G)=(vprtrb(1)+I*vprtrb(2))/2$ at the values G=qprtrb and
!!   $(vprtrb(1)-I*vprtrb(2))/2$ at $G=-qprtrb$ (integers)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  (if option==1) vpsp(nfft)=local crystal pseudopotential in real space.
!!  (if option==2) grtn(3,natom)=grads of Etot wrt tn.
!!  (if option==3) lpsstr(6)=components of local psp part of stress tensor
!!   (Cartesian coordinates, symmetric tensor) in hartree/$\textrm{bohr}^3$
!!   Store 6 unique components in order 11, 22, 33, 32, 31, 21
!!  (if option==4) dyfrlo(3,3,natom)=d2 Eei/d tn(i)/d tn(j).  (Hartrees)
!!
!! NOTES
!! Note that the present routine is tightly connected to the dfpt_vlocal.f routine,
!! that compute the derivative of the local ionic potential
!! with respect to one atomic displacement. The argument list
!! and the internal loops to be considered were sufficiently different
!! as to make the two routine different.
!!
!! PARENTS
!!      forces,prcref,prcref_PMA,respfn,setvtr
!!
!! CHILDREN
!!      mklocl_realspace,mklocl_recipspace,mklocl_wavelets,wvl_rho_abi2big
!!      xred2xcart
!!
!! SOURCE

subroutine mklocl(dtset, dyfrlo,eei,gmet,gprimd,grtn,gsqcut,lpsstr,mgfft,&
&  mpi_enreg,natom,nattyp,nfft,ngfft,nspden,ntypat,option,pawtab,ph1d,psps,qprtrb,&
&  rhog,rhor,rprimd,ucvol,vprtrb,vpsp,wvl,wvl_den,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,natom,nfft,nspden,ntypat,option
 real(dp),intent(in) :: eei,gsqcut,ucvol
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type), intent(in) :: wvl
 type(wvl_denspot_type), intent(inout) :: wvl_den
 type(pawtab_type),intent(in)  :: pawtab(ntypat*psps%usepaw)
!arrays
 integer,intent(in) :: nattyp(ntypat),ngfft(18),qprtrb(3)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: rhog(2,nfft),rprimd(3,3)
 real(dp),intent(in) :: vprtrb(2),xred(3,natom)
 real(dp),intent(in),target :: rhor(nfft,nspden)
 real(dp),intent(out) :: dyfrlo(3,3,natom),grtn(3,natom),lpsstr(6)
 real(dp),intent(inout) :: vpsp(nfft)

!Local variables-------------------------------
!scalars
 character(len=500) :: message
!arrays
 real(dp),allocatable :: xcart(:,:)
#if defined HAVE_BIGDFT
 real(dp),pointer :: rhor_ptr(:,:)
#endif

! *************************************************************************

 if (option < 1 .or. option > 4) then
   write(message,'(a,i0,a,a)')&
&   'From the calling routine, option=',option,ch10,&
&   'The only allowed values are between 1 and 4.'
   MSG_ERROR(message)
 end if
 if (option > 2 .and. .not.psps%vlspl_recipSpace) then
   write(message,'(a,i0,a,a,a,a)')&
&   'From the calling routine, option=',option,ch10,&
&   'but the local part of the pseudo-potential is in real space.',ch10,&
&   'Action: set icoulomb = 0 to turn-off real space computations.'
   MSG_ERROR(message)
 end if
 if (option > 2 .and. dtset%usewvl == 1) then
   write(message,'(a,i0,a,a)')&
&   'From the calling routine, option=',option,ch10,&
&   'but this is not implemented yet from wavelets.'
   MSG_ERROR(message)
 end if

 if (dtset%usewvl == 0) then
!  Plane wave case
   if (psps%vlspl_recipSpace) then
     call mklocl_recipspace(dyfrlo,eei,dtset%icutcoul,dtset%vcutgeo,&
&      gmet,gprimd,grtn,gsqcut,lpsstr,mgfft, &
&     mpi_enreg,psps%mqgrid_vl,natom,nattyp,nfft,ngfft, &
&     ntypat,option,ph1d,psps%qgrid_vl,qprtrb,rhog,ucvol, &
&     psps%vlspl,vprtrb,vpsp)
   else
     call mklocl_realspace(grtn,dtset%icoulomb,mpi_enreg,natom,nattyp,nfft, &
&     ngfft,dtset%nscforder,nspden,ntypat,option,pawtab,psps,rhog,rhor, &
&     rprimd,dtset%typat,ucvol,dtset%usewvl,vpsp,xred)
   end if
 else
!  Store xcart for each atom
   ABI_ALLOCATE(xcart,(3, dtset%natom))
   call xred2xcart(dtset%natom, rprimd, xcart, xred)
!  Eventually retrieve density
#if defined HAVE_BIGDFT
   if (option>1.and.wvl_den%denspot%rhov_is/=ELECTRONIC_DENSITY) then
     rhor_ptr => rhor ! Just to bypass intent(inout)
     call wvl_rho_abi2big(1,rhor_ptr,wvl_den)
   end if
#endif
!  Wavelets case
   call mklocl_wavelets(dtset%efield, grtn, mpi_enreg, dtset%natom, &
&   nfft, nspden, option, rprimd, vpsp, &
&   wvl_den, wvl, xcart)
   ABI_DEALLOCATE(xcart)
 end if

end subroutine mklocl
!!***

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

subroutine mklocl_recipspace(dyfrlo,eei,icutcoul,vcutgeo,gmet,gprimd,grtn,gsqcut,lpsstr,mgfft,&
&  mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,ntypat,option,ph1d,qgrid,qprtrb,&
&  rhog,ucvol,vlspl,vprtrb,vpsp)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,mqgrid,natom,nfft,ntypat,option,icutcoul
 real(dp),intent(in) :: eei,gsqcut,ucvol,vcutgeo(3)
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
 real(dp),allocatable :: gcutoff(:)
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

 !Initialize Gcut-off array from m_gtermcutoff
 call termcutoff(icutcoul,vcutgeo,gmet,gprimd,nfft,ngfft,gsqcut,ucvol,gcutoff)
 !Print out the method used for cut-off
 !if (icutcoul==0) mode='SPHERE'
 !if (icutcoul==1) mode='CYLINDER'
 !if (icutcoul==2) mode='SURFACE'
 !if (icutcoul==3) mode='CRYSTAL'
 !if (icutcoul==4) mode='ERF'
 !if (icutcoul==5) mode='ERFC'

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
&             cc*vlspl(jj,2,itypat)+dd*vlspl(jj+1,2,itypat) ) / gsquar * gcutoff(ii)

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
!  write(std_out,*) ' mklocl_recipspace : will add potential with strength vprtrb(:)=',vprtrb(:)

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
   call fourdp(1,work1,vpsp,1,mpi_enreg,nfft,1,ngfft,0)

!  Divide by unit cell volume
   xnorm=1.0_dp/ucvol
   vpsp(:)=vpsp(:)*xnorm

   ABI_DEALLOCATE(work1)

 end if

 ABI_DEALLOCATE(gcutoff) 

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

   real(dp) :: phr_mk,x1,x2,x3,y1,y2,y3
   phr_mk=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
 end function phr_mk

   function phi_mk(x1,y1,x2,y2,x3,y3)

   real(dp):: phi_mk,x1,x2,x3,y1,y2,y3
   phi_mk=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
 end function phi_mk

   function ph1_mk(nri,ig1,ia)

   real(dp):: ph1_mk
   integer :: nri,ig1,ia
   ph1_mk=ph1d(nri,ig1+1+n1+(ia-1)*(2*n1+1))
 end function ph1_mk

   function ph2_mk(nri,ig2,ia)

   real(dp):: ph2_mk
   integer :: nri,ig2,ia
   ph2_mk=ph1d(nri,ig2+1+n2+(ia-1)*(2*n2+1)+natom*(2*n1+1))
 end function ph2_mk

   function ph3_mk(nri,ig3,ia)

   real(dp):: ph3_mk
   integer :: nri,ig3,ia
   ph3_mk=ph1d(nri,ig3+1+n3+(ia-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1))
 end function ph3_mk

   function phre_mk(ig1,ig2,ig3,ia)

   real(dp):: phre_mk
   integer :: ig1,ig2,ig3,ia
   phre_mk=phr_mk(ph1_mk(re,ig1,ia),ph1_mk(im,ig1,ia),&
&   ph2_mk(re,ig2,ia),ph2_mk(im,ig2,ia),ph3_mk(re,ig3,ia),ph3_mk(im,ig3,ia))
 end function phre_mk

   function phimag_mk(ig1,ig2,ig3,ia)

   real(dp) :: phimag_mk
   integer :: ig1,ig2,ig3,ia
   phimag_mk=phi_mk(ph1_mk(re,ig1,ia),ph1_mk(im,ig1,ia),&
&   ph2_mk(re,ig2,ia),ph2_mk(im,ig2,ia),ph3_mk(re,ig3,ia),ph3_mk(im,ig3,ia))
 end function phimag_mk

   function gsq_mk(i1,i2,i3)

   real(dp) :: gsq_mk
   integer :: i1,i2,i3
   gsq_mk=dble(i1*i1)*gmet(1,1)+dble(i2*i2)*gmet(2,2)+&
&   dble(i3*i3)*gmet(3,3)+dble(2*i1*i2)*gmet(1,2)+&
&   dble(2*i2*i3)*gmet(2,3)+dble(2*i3*i1)*gmet(3,1)
 end function gsq_mk

end subroutine mklocl_recipspace
!!***

!!****f* ABINIT/dfpt_vlocal
!! NAME
!! dfpt_vlocal
!!
!! FUNCTION
!! Compute local part of 1st-order potential from the appropriate
!! atomic pseudopotential with structure and derivative factor.
!! In case of derivative with respect to k or
!! electric (magnetic Zeeman) field perturbation, the 1st-order local potential vanishes.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cplex: if 1, real space 1-order functions on FFT grid
!!    are REAL, if 2, COMPLEX
!!  gmet(3,3)=reciprocal space metric (Bohr**-2)
!!  gsqcut=cutoff G**2 for included G s in fft box.
!!  idir=direction of atomic displacement (=1,2 or 3 : displacement of
!!    atom ipert along the 1st, 2nd or 3rd axis).
!!  ipert=number of the atom being displaced in the frozen-phonon
!!  mpi_enreg=information about MPI parallelization
!!  mqgrid=dimension of q grid for pseudopotentials
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  ntypat=number of types of atoms in cell.
!!  n1,n2,n3=fft grid.
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information.
!!  qgrid(mqgrid)=grid of q points from 0 to qmax.
!!  qphon(3)=wavevector of the phonon
!!  ucvol=unit cell volume (Bohr**3).
!!  vlspl(mqgrid,2,ntypat)=spline fit of q^2 V(q) for each type of atom.
!!  xred(3,natom)=reduced atomic coordinates
!!
!! OUTPUT
!!  vpsp1(cplex*nfft)=first-order local crystal pseudopotential in real space
!!    (including the minus sign, forgotten in the paper non-linear..
!!
!! PARENTS
!!      dfpt_looppert,dfpt_nstdy,dfpt_nstpaw,dfptnl_loop
!!
!! CHILDREN
!!      fourdp,ptabs_fourdp
!!
!! SOURCE

subroutine dfpt_vlocal(atindx,cplex,gmet,gsqcut,idir,ipert,&
& mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,&
& ntypat,n1,n2,n3,ph1d,qgrid,qphon,ucvol,vlspl,vpsp1,xred)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,mqgrid,n1,n2,n3,natom,nfft,ntypat
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx(natom),nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: gmet(3,3),ph1d(2,(2*n1+1+2*n2+1+2*n3+1)*natom)
 real(dp),intent(in) :: qgrid(mqgrid),qphon(3),vlspl(mqgrid,2,ntypat)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(out) :: vpsp1(cplex*nfft)

!Local variables -------------------------
!scalars
 integer :: i1,i2,i3,ia1,iatom,id1,id2,id3,ig1,ig2,ig3,ii,ii1,im=2
 integer :: itypat,jj,re=1
 real(dp),parameter :: tolfix=1.000000001_dp
 real(dp) :: aa,bb,cc,cutoff,dd,diff,dq,dq2div6,dqdiv6,dqm1,gmag,gq1
 real(dp) :: gq2,gq3,gsquar,phqim,phqre
 real(dp) :: qxred2pi,sfi,sfr,vion1,xnorm
 logical :: qeq0
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: gq(3)
 real(dp),allocatable :: work1(:,:)

! *********************************************************************

 iatom=ipert

 if(iatom==natom+1 .or. iatom==natom+2 .or. iatom==natom+10  .or. iatom==natom+11 .or. iatom==natom+5)then

!  (In case of d/dk or an electric field, or magnetic (Zeeman) field->[natom+5] SPr deb )
   vpsp1(1:cplex*nfft)=zero

 else

!  (In case of a phonon perturbation)
   ABI_ALLOCATE(work1,(2,nfft))
   work1(1:2,1:nfft)=0.0_dp

   dq=(qgrid(mqgrid)-qgrid(1))/dble(mqgrid-1)
   dqm1=1.0_dp/dq
   dqdiv6=dq/6.0_dp
   dq2div6=dq**2/6.0_dp
   cutoff=gsqcut*tolfix
   id1=n1/2+2
   id2=n2/2+2
   id3=n3/2+2

   ! Get the distrib associated with this fft_grid
   call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!  This is to allow q=0
   qeq0=.false.
   if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15)qeq0=.true.

!  Determination of the atom type
   ia1=0
   itypat=0
   do ii=1,ntypat
     ia1=ia1+nattyp(ii)
     if(atindx(iatom)<=ia1.and.itypat==0)itypat=ii
   end do

!  Determination of phase qxred*
   qxred2pi=2.0_dp*pi*(qphon(1)*xred(1,iatom)+ &
&   qphon(2)*xred(2,iatom)+ &
&   qphon(3)*xred(3,iatom) )
   phqre=cos(qxred2pi)
   phqim=sin(qxred2pi)
   ii=0

   do i3=1,n3
     ig3=i3-(i3/id3)*n3-1
     gq3=dble(ig3)+qphon(3)
     gq(3)=gq3
     do i2=1,n2
       if (fftn2_distrib(i2)==mpi_enreg%me_fft) then
         ig2=i2-(i2/id2)*n2-1
         gq2=dble(ig2)+qphon(2)
         gq(2)=gq2

!        Note the lower limit of the next loop
         ii1=1
         if(i3==1 .and. i2==1 .and. qeq0 .and. ig2==0 .and. ig3==0)then
           ii1=2
           ii=ii+1
         end if
         do i1=ii1,n1
           ig1=i1-(i1/id1)*n1-1
           gq1=dble(ig1)+qphon(1)
           gq(1)=gq1
           ii=ii+1
           gsquar=gsq_vl3(gq1,gq2,gq3)
!          Skip G**2 outside cutoff:
           if (gsquar<=cutoff) then
             gmag=sqrt(gsquar)

!            Compute vion(G) for given type of atom
             jj=1+int(gmag*dqm1)
             diff=gmag-qgrid(jj)

!            Evaluate spline fit from q^2 V(q) to get V(q):
!            (p. 86 Numerical Recipes, Press et al; NOTE error in book for sign
!            of "aa" term in derivative; also see splfit routine.
!            This bug fixed here 27 Jan 1992.)

             bb = diff*dqm1
             aa = 1.0_dp-bb
             cc = aa*(aa**2-1.0_dp)*dq2div6
             dd = bb*(bb**2-1.0_dp)*dq2div6
             vion1 = (aa*vlspl(jj,1,itypat)+bb*vlspl(jj+1,1,itypat) + &
&             cc*vlspl(jj,2,itypat)+dd*vlspl(jj+1,2,itypat) ) &
&             / gsquar

!            Phase   G*xred  (complex conjugate) * -i *2pi*(g+q)*vion
             sfr=-phimag_vl3(ig1,ig2,ig3,iatom)*2.0_dp*pi*gq(idir)*vion1
             sfi=-phre_vl3(ig1,ig2,ig3,iatom)*2.0_dp*pi*gq(idir)*vion1

!            Phase   q*xred  (complex conjugate)
             work1(re,ii)=sfr*phqre+sfi*phqim
             work1(im,ii)=-sfr*phqim+sfi*phqre
           end if

         end do
       end if
     end do
   end do

!  Transform back to real space
   call fourdp(cplex,work1,vpsp1,1,mpi_enreg,nfft,1,ngfft,0)

   xnorm=1.0_dp/ucvol
   vpsp1(1:cplex*nfft)=vpsp1(1:cplex*nfft)*xnorm

   ABI_DEALLOCATE(work1)

!  End the condition of non-electric-field
 end if

 contains

!Real and imaginary parts of phase.
   function phr_vl3(x1,y1,x2,y2,x3,y3)

   real(dp) :: phr_vl3
   real(dp),intent(in) :: x1,x2,x3,y1,y2,y3
   phr_vl3=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
 end function phr_vl3

   function phi_vl3(x1,y1,x2,y2,x3,y3)

   real(dp) :: phi_vl3
   real(dp),intent(in) :: x1,x2,x3,y1,y2,y3
   phi_vl3=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
 end function phi_vl3

!  Warning : this function differ from similar ones for ground-state calculations : note the atindx !!
   function ph1_vl3(nri,ig1,ia)

   real(dp) :: ph1_vl3
   integer,intent(in) :: nri,ig1,ia
   ph1_vl3=ph1d(nri,ig1+1+n1+(atindx(ia)-1)*(2*n1+1))
 end function ph1_vl3

!  Warning : this function differ from similar ones for ground-state calculations : note the atindx !!
   function ph2_vl3(nri,ig2,ia)

   real(dp) :: ph2_vl3
   integer,intent(in) :: nri,ig2,ia
   ph2_vl3=ph1d(nri,ig2+1+n2+(atindx(ia)-1)*(2*n2+1)+natom*(2*n1+1))
 end function ph2_vl3

!  Warning : this function differ from similar ones for ground-state calculations : note the atindx !!
   function ph3_vl3(nri,ig3,ia)

   real(dp) :: ph3_vl3
   integer,intent(in) :: nri,ig3,ia
   ph3_vl3=ph1d(nri,ig3+1+n3+(atindx(ia)-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1))
 end function ph3_vl3

   function phre_vl3(ig1,ig2,ig3,ia)

   real(dp) :: phre_vl3
   integer,intent(in) :: ig1,ig2,ig3,ia
   phre_vl3=phr_vl3(ph1_vl3(re,ig1,ia),ph1_vl3(im,ig1,ia),&
&   ph2_vl3(re,ig2,ia),ph2_vl3(im,ig2,ia),ph3_vl3(re,ig3,ia),ph3_vl3(im,ig3,ia))
 end function phre_vl3

   function phimag_vl3(ig1,ig2,ig3,ia)

   real(dp) :: phimag_vl3
   integer,intent(in) :: ig1,ig2,ig3,ia
   phimag_vl3=phi_vl3(ph1_vl3(re,ig1,ia),ph1_vl3(im,ig1,ia),&
&   ph2_vl3(re,ig2,ia),ph2_vl3(im,ig2,ia),ph3_vl3(re,ig3,ia),ph3_vl3(im,ig3,ia))
 end function phimag_vl3

   function gsq_vl3(g1,g2,g3)

   real(dp) :: gsq_vl3
   real(dp),intent(in) :: g1,g2,g3 ! Note that they are real, unlike in other similar function definitions
!Define G^2 based on G space metric gmet.
   gsq_vl3=g1*g1*gmet(1,1)+g2*g2*gmet(2,2)+&
&   g3*g3*gmet(3,3)+2.0_dp*g1*g2*gmet(1,2)+&
&   2.0_dp*g2*g3*gmet(2,3)+2.0_dp*g3*g1*gmet(3,1)
 end function gsq_vl3

end subroutine dfpt_vlocal
!!***

!!****f* ABINIT/vlocalstr
!! NAME
!! vlocalstr
!!
!! FUNCTION
!! Compute strain derivatives of local ionic potential
!!                second derivative of E wrt xred
!!
!! INPUTS
!!  gmet(3,3)=reciprocal space metric ($\textrm{Bohr}^{-2}$).
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  gsqcut=cutoff on $|G|^2$: see setup1 for definition (doubled sphere).
!!  istr=1,...6 specifies cartesian strain component 11,22,33,32,31,21
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=information about MPI parallelization
!!  mqgrid=number of grid pts in q array for f(q) spline.
!!  natom=number of atoms in unit cell.
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/variables/vargs.htm#ngfft
!!  ntypat=number of types of atoms.
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information.
!!  qgrid(mqgrid)=q grid for spline from 0 to qmax.
!!  ucvol=unit cell volume ($\textrm{Bohr}^3$).
!!  vlspl(mqgrid,2,ntypat)=q^2 v(q) spline for each type of atom.
!!  [g0term]= optional, if present an alternative treatment of the G=0 term,
!!            adopoted for the flexoelectric tensor calculation, is performed.
!!
!! OUTPUT
!!  vpsp1(nfft)=first-order local crystal pseudopotential in real space.
!!
!! NOTES
!! * Note that the present routine is tightly connected to the dfpt_vlocal.f routine,
!! that compute the derivative of the local ionic potential
!! with respect to one atomic displacement. The argument list
!! and the internal loops to be considered were sufficiently different
!! as to make the two routines different.
!! * The routine was adapted from mklocl.F90
!!
!! PARENTS
!!      dfpt_looppert,dfpt_nselt,dfpt_nstpaw
!!
!! CHILDREN
!!      fourdp,ptabs_fourdp
!!
!! SOURCE

subroutine vlocalstr(gmet,gprimd,gsqcut,istr,mgfft,mpi_enreg,&
&  mqgrid,natom,nattyp,nfft,ngfft,ntypat,ph1d,qgrid,&
&  ucvol,vlspl,vpsp1,g0term)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istr,mgfft,mqgrid,natom,nfft,ntypat
 integer,optional,intent(in) :: g0term
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: qgrid(mqgrid),vlspl(mqgrid,2,ntypat)
 real(dp),intent(out) :: vpsp1(nfft)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: g0term_
 integer :: i1,i2,i3,ia,ia1,ia2,id1,id2,id3,ig1,ig2,ig3,ii,itypat,jj
 integer :: ka,kb,n1,n2,n3
 real(dp),parameter :: tolfix=1.0000001_dp
 real(dp) :: aa,bb,cc,cutoff,dd,dgsquards,diff
 real(dp) :: dq,dq2div6,dqdiv6,dqm1,ee,ff,gmag,gsquar
 real(dp) :: sfi,sfr,term,vion1,vion2,vlocg0
 real(dp) :: xnorm
 character(len=500) :: message
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: dgmetds(3,3)
 real(dp),allocatable :: work1(:,:)

! *************************************************************************

!Define G^2 based on G space metric gmet.
! gsq_vl(i1,i2,i3)=dble(i1*i1)*gmet(1,1)+dble(i2*i2)*gmet(2,2)+&
!& dble(i3*i3)*gmet(3,3)+dble(2*i1*i2)*gmet(1,2)+&
!& dble(2*i2*i3)*gmet(2,3)+dble(2*i3*i1)*gmet(3,1)

!Define dG^2/ds based on G space metric derivative dgmetds.
! dgsqds_vl(i1,i2,i3)=dble(i1*i1)*dgmetds(1,1)+dble(i2*i2)*dgmetds(2,2)+&
!& dble(i3*i3)*dgmetds(3,3)+&
!& dble(i1*i2)*(dgmetds(1,2)+dgmetds(2,1))+&
!& dble(i1*i3)*(dgmetds(1,3)+dgmetds(3,1))+&
!& dble(i2*i3)*(dgmetds(2,3)+dgmetds(3,2))

!Real and imaginary parts of phase--statment functions:
! phr_vl(x1,y1,x2,y2,x3,y3)=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
! phi_vl(x1,y1,x2,y2,x3,y3)=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
! ph1_vl(nri,i1,ia)=ph1d(nri,i1+1+n1+(ia-1)*(2*n1+1))
! ph2_vl(nri,i2,ia)=ph1d(nri,i2+1+n2+(ia-1)*(2*n2+1)+&
!& natom*(2*n1+1))
! ph3_vl(nri,i3,ia)=ph1d(nri,i3+1+n3+(ia-1)*(2*n3+1)+&
!& natom*(2*n1+1+2*n2+1))
! phre_vl(i1,i2,i3,ia)=phr_vl(ph1_vl(re,i1,ia),ph1_vl(im,i1,ia),ph2_vl(re,i2,ia),&
!& ph2_vl(im,i2,ia),ph3_vl(re,i3,ia),ph3_vl(im,i3,ia))
! phimag_vl(i1,i2,i3,ia)=phi_vl(ph1_vl(re,i1,ia),ph1_vl(im,i1,ia),ph2_vl(re,i2,ia),&
!& ph2_vl(im,i2,ia),ph3_vl(re,i3,ia),ph3_vl(im,i3,ia))

!-----
!Compute derivative of metric tensor wrt strain component istr
 if(istr<1 .or. istr>6)then
   write(message, '(a,i10,a,a,a)' )&
&   ' Input istr=',istr,' not allowed.',ch10,&
&   ' Possible values are 1,2,3,4,5,6 only.'
   MSG_BUG(message)
 end if

 ka=idx(2*istr-1);kb=idx(2*istr)
 do ii = 1,3
   dgmetds(:,ii)=-(gprimd(ka,:)*gprimd(kb,ii)+gprimd(kb,:)*gprimd(ka,ii))
 end do
!For historical reasons:
 dgmetds(:,:)=0.5_dp*dgmetds(:,:)

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)

!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!Zero out array to permit accumulation over atom types below:
 ABI_ALLOCATE(work1,(2,nfft))
 work1(:,:)=0.0_dp
!
 dq=(qgrid(mqgrid)-qgrid(1))/dble(mqgrid-1)
 dqm1=1.0_dp/dq
 dqdiv6=dq/6.0_dp
 dq2div6=dq**2/6.0_dp
 cutoff=gsqcut*tolfix
 id1=n1/2+2
 id2=n2/2+2
 id3=n3/2+2

 ia1=1
 do itypat=1,ntypat
!  ia1,ia2 sets range of loop over atoms:
   ia2=ia1+nattyp(itypat)-1

   ii=0
   do i3=1,n3
     ig3=i3-(i3/id3)*n3-1
     do i2=1,n2
       if (fftn2_distrib(i2)==mpi_enreg%me_fft) then
         ig2=i2-(i2/id2)*n2-1
         do i1=1,n1
           ig1=i1-(i1/id1)*n1-1
           ii=ii+1
!          ***     GET RID OF THIS THESE IF STATEMENTS (if they slow code)
!          Skip G=0:
!          if (ii==1) cycle
           if (ig1==0 .and. ig2==0 .and. ig3==0) cycle
           gsquar=gsq_vl(ig1,ig2,ig3)

!          Skip G**2 outside cutoff:
           if (gsquar<=cutoff) then
             gmag=sqrt(gsquar)
             dgsquards=dgsqds_vl(ig1,ig2,ig3)
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
&             cc*vlspl(jj,2,itypat)+dd*vlspl(jj+1,2,itypat) ) &
&             / gsquar

!            Also get (dV(q)/dq)/q:
!            (note correction of Numerical Recipes sign error
!            before (3._dp*aa**2-1._dp)
             ee= vlspl(jj+1,1,itypat)-vlspl(jj,1,itypat)
             ff=  (3._dp*bb**2-1._dp)*vlspl(jj+1,2,itypat) &
&             - (3._dp*aa**2-1._dp)*vlspl(jj,2,itypat)
             vion2 = ( ( ee*dqm1 + ff*dqdiv6 )/gmag&
&             - 2.0_dp*vion1                 ) / gsquar


!            Assemble structure factor over all atoms of given type:
             sfr=0.0_dp
             sfi=0.0_dp
             do ia=ia1,ia2
               sfr=sfr+phre_vl(ig1,ig2,ig3,ia)
               sfi=sfi-phimag_vl(ig1,ig2,ig3,ia)
             end do

             term=dgsquards*vion2
!            Add potential for diagonal strain components
             if(istr <=3) then
               term=term-vion1
             end if

!            Multiply structure factor times vion derivatives:
             work1(re,ii)=work1(re,ii)+sfr*term
             work1(im,ii)=work1(im,ii)+sfi*term

!            End skip G**2 outside cutoff:
           end if
!          End loop on n1, n2, n3. There is a "cycle" inside the loop
         end do
       end if
     end do
   end do

   ia1=ia2+1

!  End loop on type of atoms
 end do


!Set Vloc(G=0)=0:
 work1(re,1)=0.0_dp
 work1(im,1)=0.0_dp

!Alternative treatment of Vloc(G=0) for the flexoelectric tensor calculation
 g0term_=0; if (present(g0term)) g0term_=g0term
 if (g0term_==1) then
   vlocg0=zero
   if (istr<=3) then
     ia1=1
     do itypat=1,ntypat
    !  ia1,ia2 sets range of loop over atoms:

       ia2=ia1+nattyp(itypat)-1
       do ia=ia1,ia2
         vlocg0=vlocg0+vlspl(1,2,itypat)
       end do
     end do
     work1(re,1)=-half*vlocg0
   end if
 end if

!Transform back to real space
 call fourdp(1,work1,vpsp1,1,mpi_enreg,nfft,1,ngfft,0)

!Divide by unit cell volume
 xnorm=1.0_dp/ucvol
 vpsp1(:)=vpsp1(:)*xnorm

 ABI_DEALLOCATE(work1)

 contains

!Real and imaginary parts of phase.
   function phr_vl(x1,y1,x2,y2,x3,y3)

   real(dp) :: phr_vl,x1,x2,x3,y1,y2,y3
   phr_vl=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
 end function phr_vl

   function phi_vl(x1,y1,x2,y2,x3,y3)

   real(dp):: phi_vl,x1,x2,x3,y1,y2,y3
   phi_vl=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
 end function phi_vl

   function ph1_vl(nri,ig1,ia)

   real(dp):: ph1_vl
   integer :: nri,ig1,ia
   ph1_vl=ph1d(nri,ig1+1+n1+(ia-1)*(2*n1+1))
 end function ph1_vl

   function ph2_vl(nri,ig2,ia)

   real(dp):: ph2_vl
   integer :: nri,ig2,ia
   ph2_vl=ph1d(nri,ig2+1+n2+(ia-1)*(2*n2+1)+natom*(2*n1+1))
 end function ph2_vl

   function ph3_vl(nri,ig3,ia)

   real(dp):: ph3_vl
   integer :: nri,ig3,ia
   ph3_vl=ph1d(nri,ig3+1+n3+(ia-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1))
 end function ph3_vl

   function phre_vl(ig1,ig2,ig3,ia)

   real(dp):: phre_vl
   integer :: ig1,ig2,ig3,ia
   phre_vl=phr_vl(ph1_vl(re,ig1,ia),ph1_vl(im,ig1,ia),&
&   ph2_vl(re,ig2,ia),ph2_vl(im,ig2,ia),ph3_vl(re,ig3,ia),ph3_vl(im,ig3,ia))
 end function phre_vl

   function phimag_vl(ig1,ig2,ig3,ia)

   real(dp) :: phimag_vl
   integer :: ig1,ig2,ig3,ia
   phimag_vl=phi_vl(ph1_vl(re,ig1,ia),ph1_vl(im,ig1,ia),&
&   ph2_vl(re,ig2,ia),ph2_vl(im,ig2,ia),ph3_vl(re,ig3,ia),ph3_vl(im,ig3,ia))
 end function phimag_vl

   function gsq_vl(i1,i2,i3)

   real(dp) :: gsq_vl
   integer :: i1,i2,i3
!Define G^2 based on G space metric gmet.
   gsq_vl=dble(i1*i1)*gmet(1,1)+dble(i2*i2)*gmet(2,2)+&
&   dble(i3*i3)*gmet(3,3)+dble(2*i1*i2)*gmet(1,2)+&
&   dble(2*i2*i3)*gmet(2,3)+dble(2*i3*i1)*gmet(3,1)
 end function gsq_vl

   function dgsqds_vl(i1,i2,i3)

   real(dp) :: dgsqds_vl
   integer :: i1,i2,i3
!Define dG^2/ds based on G space metric derivative dgmetds.
   dgsqds_vl=dble(i1*i1)*dgmetds(1,1)+dble(i2*i2)*dgmetds(2,2)+&
&   dble(i3*i3)*dgmetds(3,3)+&
&   dble(i1*i2)*(dgmetds(1,2)+dgmetds(2,1))+&
&   dble(i1*i3)*(dgmetds(1,3)+dgmetds(3,1))+&
&   dble(i2*i3)*(dgmetds(2,3)+dgmetds(3,2))
 end function dgsqds_vl

end subroutine vlocalstr
!!***

!!****f* ABINIT/dfpt_vlocaldq
!! NAME
!! dfpt_vlocaldq
!!
!! FUNCTION
!! Compute q-gradient (at q=0) of the local part of 1st-order 
!! atomic displacement potential from the appropriate
!! atomic pseudopotential with structure and derivative factor.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (MR,MS)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cplex: if 1, real space 1-order functions on FFT grid
!!    are REAL, if 2, COMPLEX
!!  gmet(3,3)=reciprocal space metric (Bohr**-2)
!!  gsqcut=cutoff G**2 for included G s in fft box.
!!  idir=direction of atomic displacement (=1,2 or 3 : displacement of
!!    atom ipert along the 1st, 2nd or 3rd axis).
!!  ipert=number of the atom being displaced in the frozen-phonon
!!  mpi_enreg=information about MPI parallelization
!!  mqgrid=dimension of q grid for pseudopotentials
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ntypat=number of types of atoms in cell.
!!  n1,n2,n3=fft grid.
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information.
!!  qdir=direction of the q-gradient 
!!  qgrid(mqgrid)=grid of q points from 0 to qmax.
!!  qphon(3)=wavevector of the phonon
!!  ucvol=unit cell volume (Bohr**3).
!!  vlspl(mqgrid,2,ntypat)=spline fit of q^2 V(q) for each type of atom.
!!
!! OUTPUT
!!  vpsp1dq(cplex*nfft)=q-gradient (at q=0) of the first-order local 
!!  crystal pseudopotential in real space
!!    (including the minus sign, forgotten in the paper non-linear..
!!
!! NOTES
!! * IMPORTANT: the formalism followed in this routine
!!   assumes a phase factor for the perturbation that 
!!   is different to the one used elsewhere in the code (See M.Stengel paper):
!!         
!!             here: e^{i q (R_l + \tau_{\kappa})}
!!   rest of ABINIT: e^{i q R_l}  
!!
!!  **A -i factor has been factorized out in all the contributions of the first
!!    q-gradient of the atomic displacement Hamiltonian. This is lately included 
!!    in the matrix element calculation.
!!
!! PARENTS
!!      dfpt_qdrpwf
!!
!! CHILDREN
!!      fourdp,ptabs_fourdp,splfit
!!
!! SOURCE

subroutine dfpt_vlocaldq(atindx,cplex,gmet,gsqcut,idir,ipert,&
& mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,&
& ntypat,n1,n2,n3,ph1d,qdir,qgrid,qphon,ucvol,vlspl,vpsp1dq)

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,mqgrid,n1,n2,n3,natom,nfft,ntypat
 integer,intent(in) :: qdir
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx(natom),nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: gmet(3,3),ph1d(2,(2*n1+1+2*n2+1+2*n3+1)*natom)
 real(dp),intent(in) :: qgrid(mqgrid),qphon(3),vlspl(mqgrid,2,ntypat)
 real(dp),intent(out) :: vpsp1dq(cplex*nfft)

!Local variables -------------------------
!scalars
 integer :: i1,i2,i3,ia1,iatom,id1,id2,id3,ig1,ig2,ig3,ii,ii1,im=2
 integer :: itypat,re=1
 real(dp),parameter :: tolfix=1.000000001_dp
 real(dp) :: cutoff,gfact,gmag,gq1
 real(dp) :: gq2,gq3,gsquar
 real(dp) :: sfi,sfr,xnorm
 logical :: qeq0
 character(len=500) :: msg
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: gq(3),gvec(3),vion1(1),vion1dq(1)
 real(dp),allocatable :: work1(:,:)

! *********************************************************************

 iatom=ipert

 if(iatom==natom+1 .or. iatom==natom+2 .or. iatom==natom+10  .or. iatom==natom+11 .or. iatom==natom+5)then

!  (In case of d/dk or an electric field, or magnetic (Zeeman) field->[natom+5] SPr deb )
   vpsp1dq(1:cplex*nfft)=zero

 else

!  (In case of a phonon perturbation)
   ABI_ALLOCATE(work1,(2,nfft))
   work1(1:2,1:nfft)=0.0_dp

   cutoff=gsqcut*tolfix
   id1=n1/2+2
   id2=n2/2+2
   id3=n3/2+2

   ! Get the distrib associated with this fft_grid
   call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!  This is to allow q=0
   qeq0=.false.
   if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15) then 
     qeq0=.true.
   else
     msg='This routine cannot be used for q/=0'
     MSG_BUG(msg)
   end if

!  Determination of the atom type
   ia1=0
   itypat=0
   do ii=1,ntypat
     ia1=ia1+nattyp(ii)
     if(atindx(iatom)<=ia1.and.itypat==0)itypat=ii
   end do

   ii=0

   do i3=1,n3
     ig3=i3-(i3/id3)*n3-1
     gq3=dble(ig3)+qphon(3)
     gq(3)=gq3
     do i2=1,n2
       if (fftn2_distrib(i2)==mpi_enreg%me_fft) then
         ig2=i2-(i2/id2)*n2-1
         gq2=dble(ig2)+qphon(2)
         gq(2)=gq2

!        Note the lower limit of the next loop
         ii1=1
         if(i3==1 .and. i2==1 .and. qeq0 .and. ig2==0 .and. ig3==0)then
           ii1=2
           ii=ii+1
         end if
         do i1=ii1,n1
           ig1=i1-(i1/id1)*n1-1
           gq1=dble(ig1)+qphon(1)
           gq(1)=gq1
           ii=ii+1
           gsquar=gsq_vl3(gq1,gq2,gq3)
!          Skip G**2 outside cutoff:
           if (gsquar<=cutoff) then
             gmag=sqrt(gsquar)

!            Evaluate spline fit to get V(q) and V(q)':
             call splfit(qgrid,vion1dq,vlspl(:,:,itypat),1,(/gmag/),vion1,mqgrid,1)

             vion1=vion1/gsquar
             vion1dq=(vion1dq-2.0_dp*gmag*vion1)/gsquar

             gvec=(/ig1,ig2,ig3/)
             gfact=dot_product(gmet(qdir,:),gvec(:))/gmag

!            Phase   G*xred  (complex conjugate) *2*pi*(g_idir)*vion1dq*gfact
             sfr=phre_vl3(ig1,ig2,ig3,iatom)*two_pi*gq(idir)*vion1dq(1)*gfact
             sfi=-phimag_vl3(ig1,ig2,ig3,iatom)*two_pi*gq(idir)*vion1dq(1)*gfact
!            Phase   G*xred  (complex conjugate) *2*pi*(\delta_{idir,qdir})*vion1
             if (idir==qdir) then
               sfr=sfr + phre_vl3(ig1,ig2,ig3,iatom)*two_pi*vion1(1)
               sfi=sfi - phimag_vl3(ig1,ig2,ig3,iatom)*two_pi*vion1(1)
             end if         

             work1(re,ii)=sfr
             work1(im,ii)=sfi
           end if

         end do
       end if
     end do
   end do

!  Transform back to real space
   call fourdp(cplex,work1,vpsp1dq,1,mpi_enreg,nfft,1,ngfft,0)

   xnorm=1.0_dp/ucvol
   vpsp1dq(1:cplex*nfft)=vpsp1dq(1:cplex*nfft)*xnorm

   ABI_DEALLOCATE(work1)

!  End the condition of non-electric-field
 end if

 contains

!Real and imaginary parts of phase.
 function phr_vl3(x1,y1,x2,y2,x3,y3)
   real(dp) :: phr_vl3
   real(dp),intent(in) :: x1,x2,x3,y1,y2,y3
   phr_vl3=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
 end function phr_vl3

 function phi_vl3(x1,y1,x2,y2,x3,y3)
   real(dp) :: phi_vl3
   real(dp),intent(in) :: x1,x2,x3,y1,y2,y3
   phi_vl3=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
 end function phi_vl3

!  Warning : this function differ from similar ones for ground-state calculations : note the atindx !!
 function ph1_vl3(nri,ig1,ia)
   real(dp) :: ph1_vl3
   integer,intent(in) :: nri,ig1,ia
   ph1_vl3=ph1d(nri,ig1+1+n1+(atindx(ia)-1)*(2*n1+1))
 end function ph1_vl3

!  Warning : this function differ from similar ones for ground-state calculations : note the atindx !!
 function ph2_vl3(nri,ig2,ia)
   real(dp) :: ph2_vl3
   integer,intent(in) :: nri,ig2,ia
   ph2_vl3=ph1d(nri,ig2+1+n2+(atindx(ia)-1)*(2*n2+1)+natom*(2*n1+1))
 end function ph2_vl3

!  Warning : this function differ from similar ones for ground-state calculations : note the atindx !!
 function ph3_vl3(nri,ig3,ia)
   real(dp) :: ph3_vl3
   integer,intent(in) :: nri,ig3,ia
   ph3_vl3=ph1d(nri,ig3+1+n3+(atindx(ia)-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1))
 end function ph3_vl3

 function phre_vl3(ig1,ig2,ig3,ia)
   real(dp) :: phre_vl3
   integer,intent(in) :: ig1,ig2,ig3,ia
   phre_vl3=phr_vl3(ph1_vl3(re,ig1,ia),ph1_vl3(im,ig1,ia),&
&   ph2_vl3(re,ig2,ia),ph2_vl3(im,ig2,ia),ph3_vl3(re,ig3,ia),ph3_vl3(im,ig3,ia))
 end function phre_vl3

 function phimag_vl3(ig1,ig2,ig3,ia)
   real(dp) :: phimag_vl3
   integer,intent(in) :: ig1,ig2,ig3,ia
   phimag_vl3=phi_vl3(ph1_vl3(re,ig1,ia),ph1_vl3(im,ig1,ia),&
&   ph2_vl3(re,ig2,ia),ph2_vl3(im,ig2,ia),ph3_vl3(re,ig3,ia),ph3_vl3(im,ig3,ia))
 end function phimag_vl3

 function gsq_vl3(g1,g2,g3)
   real(dp) :: gsq_vl3
   real(dp),intent(in) :: g1,g2,g3 ! Note that they are real, unlike in other similar function definitions
!Define G^2 based on G space metric gmet.
   gsq_vl3=g1*g1*gmet(1,1)+g2*g2*gmet(2,2)+&
&   g3*g3*gmet(3,3)+2.0_dp*g1*g2*gmet(1,2)+&
&   2.0_dp*g2*g3*gmet(2,3)+2.0_dp*g3*g1*gmet(3,1)
 end function gsq_vl3

end subroutine dfpt_vlocaldq
!!***

!!****f* ABINIT/dfpt_vlocaldqdq
!! NAME
!! dfpt_vlocaldqdq
!!
!! FUNCTION
!! Compute 2nd q-gradient (at q=0) of the local part of 1st-order 
!! atomic displacement potential from the appropriate
!! atomic pseudopotential with structure and derivative factor.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (MR,MS)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cplex: if 1, real space 1-order functions on FFT grid
!!    are REAL, if 2, COMPLEX
!!  gmet(3,3)=reciprocal space metric (Bohr**-2)
!!  gsqcut=cutoff G**2 for included G s in fft box.
!!  idir=direction of atomic displacement (=1,2 or 3 : displacement of
!!    atom ipert along the 1st, 2nd or 3rd axis).
!!  ipert=number of the atom being displaced in the frozen-phonon
!!  mpi_enreg=information about MPI parallelization
!!  mqgrid=dimension of q grid for pseudopotentials
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ntypat=number of types of atoms in cell.
!!  n1,n2,n3=fft grid.
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information.
!!  qdir1=direction of the first q-gradient 
!!  qdir2=direction of the second q-gradient 
!!  qgrid(mqgrid)=grid of q points from 0 to qmax.
!!  qphon(3)=wavevector of the phonon
!!  ucvol=unit cell volume (Bohr**3).
!!  vlspl(mqgrid,2,ntypat)=spline fit of q^2 V(q) for each type of atom.
!!
!! OUTPUT
!!  vpsp1dqdq(cplex*nfft)=2nd q-gradient (at q=0) of the first-order local 
!!  crystal pseudopotential in real space
!!
!! NOTES
!! * IMPORTANT: the formalism followed in this routine
!!   assumes a phase factor for the perturbation that 
!!   is different to the one used elsewhere in the code (See M.Stengel paper):
!!         
!!             here: e^{i q (R_l + \tau_{\kappa})}
!!   rest of ABINIT: e^{i q R_l}  
!!
!!  **A -i factor has been factorized out in all the contributions of the second
!!    q-gradient of the atomic displacement Hamiltonian. This is lately included 
!!    in the whole frozen contribution to the q-gradient of the 
!!    2nd order energy wrt an atomic displacement and a strain:
!!    \Delta E^{\tau_{\kappa\alpha}^* (\beta)}_{m\kvec,\gamma\delta}
!!     
!!
!! PARENTS
!!      dfpt_qdrpwf
!!
!! CHILDREN
!!      fourdp,ptabs_fourdp,splfit
!!
!! SOURCE

subroutine dfpt_vlocaldqdq(atindx,cplex,gmet,gsqcut,idir,ipert,&
& mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,&
& ntypat,n1,n2,n3,ph1d,qdir1,qdir2,qgrid,qphon,ucvol,vlspl,vpsp1dqdq)

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,mqgrid,n1,n2,n3,natom,nfft,ntypat
 integer,intent(in) :: qdir1,qdir2
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx(natom),nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: gmet(3,3),ph1d(2,(2*n1+1+2*n2+1+2*n3+1)*natom)
 real(dp),intent(in) :: qgrid(mqgrid),qphon(3),vlspl(mqgrid,2,ntypat)
 real(dp),intent(out) :: vpsp1dqdq(cplex*nfft)

!Local variables -------------------------
!scalars
 integer :: alpha, delta, gamma
 integer :: i1,i2,i3,ia1,iatom,id1,id2,id3,ig1,ig2,ig3,ii,ii1,im=2
 integer :: itypat,re=1
 real(dp),parameter :: tolfix=1.000000001_dp
 real(dp) :: cutoff,delad,delag,gfact,gfact1,gfact2,gmag,gq1
 real(dp) :: gq2,gq3,gsquar
 real(dp) :: sfi,sfr,term1,term2,xnorm
 logical :: qeq0
 character(len=500) :: msg
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: gq(3),gvec(3),vion1(1),vion1dq(1),vion1dqdq(1)
 real(dp),allocatable :: work1(:,:)

! *********************************************************************

 iatom=ipert

 if(iatom==natom+1 .or. iatom==natom+2 .or. iatom==natom+10  .or. iatom==natom+11 .or. iatom==natom+5)then

!  (In case of d/dk or an electric field, or magnetic (Zeeman) field->[natom+5] SPr deb )
   vpsp1dqdq(1:cplex*nfft)=zero

 else

   alpha=idir; delta=qdir2; gamma=qdir1

   !Kronecker deltas
   delad=0.0_dp; delag=0.0_dp
   if (alpha==delta) delad=1.0_dp
   if (alpha==gamma) delag=1.0_dp

!  (In case of a phonon perturbation)
   ABI_ALLOCATE(work1,(2,nfft))
   work1(1:2,1:nfft)=0.0_dp

   cutoff=gsqcut*tolfix
   id1=n1/2+2
   id2=n2/2+2
   id3=n3/2+2

   ! Get the distrib associated with this fft_grid
   call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!  This is to allow q=0
   qeq0=.false.
   if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15) then 
     qeq0=.true.
   else
     msg='This routine cannot be used for q/=0'
     MSG_BUG(msg)
   end if

!  Determination of the atom type
   ia1=0
   itypat=0
   do ii=1,ntypat
     ia1=ia1+nattyp(ii)
     if(atindx(iatom)<=ia1.and.itypat==0)itypat=ii
   end do

   ii=0

   do i3=1,n3
     ig3=i3-(i3/id3)*n3-1
     gq3=dble(ig3)+qphon(3)
     gq(3)=gq3
     do i2=1,n2
       if (fftn2_distrib(i2)==mpi_enreg%me_fft) then
         ig2=i2-(i2/id2)*n2-1
         gq2=dble(ig2)+qphon(2)
         gq(2)=gq2

!        Note the lower limit of the next loop
         ii1=1
         if(i3==1 .and. i2==1 .and. qeq0 .and. ig2==0 .and. ig3==0)then
           ii1=2
           ii=ii+1
         end if
         do i1=ii1,n1
           ig1=i1-(i1/id1)*n1-1
           gq1=dble(ig1)+qphon(1)
           gq(1)=gq1
           ii=ii+1
           gsquar=gsq_vl3(gq1,gq2,gq3)
!          Skip G**2 outside cutoff:
           if (gsquar<=cutoff) then
             gmag=sqrt(gsquar)

!            Evaluate spline fit to get first V(q) and V(q)' and later V(q)'':
             call splfit(qgrid,vion1dq,vlspl(:,:,itypat),1,(/gmag/),vion1,mqgrid,1)
             vion1=vion1/gsquar
             vion1dq=(vion1dq-2.0_dp*gmag*vion1)/gsquar

             call splfit(qgrid,vion1dqdq,vlspl(:,:,itypat),2,(/gmag/),vion1,mqgrid,1)
             vion1dqdq=(vion1dqdq-4.0_dp*gmag*vion1dq-2.0_dp*vion1)/gsquar

             gvec=(/ig1,ig2,ig3/)
             gfact1=dot_product(gmet(gamma,:),gvec(:))
             gfact2=dot_product(gmet(delta,:),gvec(:))
             gfact=gvec(alpha)*gfact1*gfact2/gsquar

             term1=delag*gfact2+delad*gfact1+gvec(alpha)*gmet(gamma,delta)
             term1=term1-gfact
             term1=term1*vion1dq(1)/gmag

             term2=vion1dqdq(1)*gfact

!            structure factors
             sfr=phre_vl3(ig1,ig2,ig3,iatom)
             sfi=-phimag_vl3(ig1,ig2,ig3,iatom)

!            Multiply structure factor times vion derivatives:
             work1(re,ii)=sfr*(term1+term2)*two_pi
             work1(im,ii)=sfi*(term1+term2)*two_pi

           end if

         end do
       end if
     end do
   end do

!  Transform back to real space
   call fourdp(cplex,work1,vpsp1dqdq,1,mpi_enreg,nfft,1,ngfft,0)

   xnorm=1.0_dp/ucvol
   vpsp1dqdq(1:cplex*nfft)=vpsp1dqdq(1:cplex*nfft)*xnorm

   ABI_DEALLOCATE(work1)

!  End the condition of non-electric-field
 end if

 contains

!Real and imaginary parts of phase.
 function phr_vl3(x1,y1,x2,y2,x3,y3)
   real(dp) :: phr_vl3
   real(dp),intent(in) :: x1,x2,x3,y1,y2,y3
   phr_vl3=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
 end function phr_vl3

 function phi_vl3(x1,y1,x2,y2,x3,y3)
   real(dp) :: phi_vl3
   real(dp),intent(in) :: x1,x2,x3,y1,y2,y3
   phi_vl3=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
 end function phi_vl3

!  Warning : this function differ from similar ones for ground-state calculations : note the atindx !!
 function ph1_vl3(nri,ig1,ia)
   real(dp) :: ph1_vl3
   integer,intent(in) :: nri,ig1,ia
   ph1_vl3=ph1d(nri,ig1+1+n1+(atindx(ia)-1)*(2*n1+1))
 end function ph1_vl3

!  Warning : this function differ from similar ones for ground-state calculations : note the atindx !!
 function ph2_vl3(nri,ig2,ia)
   real(dp) :: ph2_vl3
   integer,intent(in) :: nri,ig2,ia
   ph2_vl3=ph1d(nri,ig2+1+n2+(atindx(ia)-1)*(2*n2+1)+natom*(2*n1+1))
 end function ph2_vl3

!  Warning : this function differ from similar ones for ground-state calculations : note the atindx !!
 function ph3_vl3(nri,ig3,ia)
   real(dp) :: ph3_vl3
   integer,intent(in) :: nri,ig3,ia
   ph3_vl3=ph1d(nri,ig3+1+n3+(atindx(ia)-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1))
 end function ph3_vl3

 function phre_vl3(ig1,ig2,ig3,ia)
   real(dp) :: phre_vl3
   integer,intent(in) :: ig1,ig2,ig3,ia
   phre_vl3=phr_vl3(ph1_vl3(re,ig1,ia),ph1_vl3(im,ig1,ia),&
&   ph2_vl3(re,ig2,ia),ph2_vl3(im,ig2,ia),ph3_vl3(re,ig3,ia),ph3_vl3(im,ig3,ia))
 end function phre_vl3

 function phimag_vl3(ig1,ig2,ig3,ia)
   real(dp) :: phimag_vl3
   integer,intent(in) :: ig1,ig2,ig3,ia
   phimag_vl3=phi_vl3(ph1_vl3(re,ig1,ia),ph1_vl3(im,ig1,ia),&
&   ph2_vl3(re,ig2,ia),ph2_vl3(im,ig2,ia),ph3_vl3(re,ig3,ia),ph3_vl3(im,ig3,ia))
 end function phimag_vl3

 function gsq_vl3(g1,g2,g3)
   real(dp) :: gsq_vl3
   real(dp),intent(in) :: g1,g2,g3 ! Note that they are real, unlike in other similar function definitions
!Define G^2 based on G space metric gmet.
   gsq_vl3=g1*g1*gmet(1,1)+g2*g2*gmet(2,2)+&
&   g3*g3*gmet(3,3)+2.0_dp*g1*g2*gmet(1,2)+&
&   2.0_dp*g2*g3*gmet(2,3)+2.0_dp*g3*g1*gmet(3,1)
 end function gsq_vl3

end subroutine dfpt_vlocaldqdq
!!***


!!****f* ABINIT/dfpt_vmetdqdq	
!! NAME
!! dfpt_vmetdqdq
!!
!! FUNCTION
!! Compute second q-gradient (at q=0) of the local part of 1st-order 
!! metric potential from the appropriate atomic pseudopotential 
!! with structure and derivative factor. Additionaly, compute the 
!! second q-gradient (at q=0) of the Hartree potential of the metric 
!! perturbation.
!! Cartesian coordinates are employed to define the direction of the 
!! metric perturbation and the two q-gradients.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (MR,MS)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex: if 1, real space 1-order functions on FFT grid
!!    are REAL, if 2, COMPLEX
!!  gmet(3,3)=reciprocal space metric (Bohr**-2)
!!  gsqcut=cutoff G**2 for included G s in fft box.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  idir= strain perturbation direction
!!  ipert=number of the atom being displaced in the frozen-phonon
!!  mpi_enreg=information about MPI parallelization
!!  mqgrid=dimension of q grid for pseudopotentials
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ntypat=number of types of atoms in cell.
!!  n1,n2,n3=fft grid.
!!  opthartdqdq= if 1 activates the calculation 2nd q-gradient of the Hartree potential 
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information.
!!  qdir=direction of the q-gradient 
!!  qgrid(mqgrid)=grid of q points from 0 to qmax.
!!  qphon(3)=wavevector of the phonon
!!  rhog(2,nfft)=array for Fourier transform of GS electron density
!!  ucvol=unit cell volume (Bohr**3).
!!  vlspl(mqgrid,2,ntypat)=spline fit of q^2 V(q) for each type of atom.
!!
!! OUTPUT
!!  vhart1dqdq(cplex*nfft)=2nd q-gradient (at q=0) of the GS density Hartree potential from the metric perturbation
!!
!!  vpsp1dqdq(cplex*nfft)=2nd q-gradient (at q=0) of the first-order metric local 
!!  crystal pseudopotential in real space
!!
!! NOTES
!! ** IMPORTANT: the formalism followed in this routine
!!    assumes a phase factor for the perturbation that 
!!    is different to the one used elsewhere in the code (See M.Stengel paper):
!!         
!!             here: e^{i q (R_l + \tau_{\kappa})}
!!    rest of ABINIT: e^{i q R_l}  
!!
!!  **Since the 2nd derivative w.r.t q-vector is calculated along cartesian
!!    directions, the 1/twopi**2 factor (that in the rest of the code is applied
!!    in the reduced to cartesian derivative conversion process) is here 
!!    explicictly included in the formulas.
!!  
!!  **Notice that idir=1-9, in contrast to the strain perturbation (idir=1-6),
!!    because this term is not symmetric w.r.t permutations of the two strain
!!    indices.
!!
!!  **A -i factor has been factorized out in all the contributions of the second
!!    q-gradient of the metric Hamiltonian. This is lately included in the contribution
!!    of the corresponing term (T4) to the flexoelectric tensor in dfpt_flexoout.F90
!!
!! PARENTS
!!
!!      dfpt_flexowf
!!
!! CHILDREN
!!
!!      fourdp,ptabs_fourdp,splfit
!!
!! SOURCE

subroutine dfpt_vmetdqdq(cplex,gmet,gprimd,gsqcut,idir,ipert,&
& mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,&
& ntypat,n1,n2,n3,opthartdqdq,ph1d,qdir,qgrid,qphon,rhog,&
& ucvol,vlspl,vhart1dqdq,vpsp1dqdq)

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,mqgrid,n1,n2,n3,natom,nfft,ntypat
 integer,intent(in) :: opthartdqdq,qdir
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,(2*n1+1+2*n2+1+2*n3+1)*natom)
 real(dp),intent(in) :: qgrid(mqgrid),qphon(3), rhog(2,nfft)
 real(dp),intent(in) :: vlspl(mqgrid,2,ntypat)
 real(dp),intent(out) :: vhart1dqdq(cplex*nfft),vpsp1dqdq(cplex*nfft)

!Local variables -------------------------
!scalars
 integer :: beta, delta, gamma
 integer :: ia,i1,i2,i3,ia1,ia2,id1,id2,id3,ig1,ig2,ig3,ii,ii1,im=2
 integer :: itypat,re=1
 real(dp),parameter :: tolfix=1.000000001_dp
 real(dp) :: cutoff,delbd,delbg,deldg,gfact,gmag,gq1
 real(dp) :: gq2,gq3,gsquar,pisqrinv
 real(dp) :: sfi,sfr,term1,term2,uogsquar,work1re,xnorm
 logical :: qeq0
 character(len=500) :: msg
!arrays
 integer,save :: idx(18)=(/1,1,2,2,3,3,3,2,3,1,2,1,2,3,1,3,1,2/)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: gq(3),gqc(3),vion1(1),vion1dq(1),vion1dqdq(1)
 real(dp),allocatable :: work1(:,:)

! *********************************************************************


 if(ipert/=natom+3 .and. ipert/=natom+4)then

   vpsp1dqdq(1:cplex*nfft)=zero

 else

   beta=idx(2*idir-1); delta=idx(2*idir); gamma=qdir

   !Kronecker deltas
   delbd=0.0_dp; delbg=0.0_dp; deldg=0.0_dp
   if (beta==delta) delbd=1.0_dp
   if (beta==gamma) delbg=1.0_dp
   if (delta==gamma) deldg=1.0_dp

   ABI_ALLOCATE(work1,(2,nfft))
   work1(1:2,1:nfft)=0.0_dp

   cutoff=gsqcut*tolfix
   id1=n1/2+2
   id2=n2/2+2
   id3=n3/2+2

   !Get the distrib associated with this fft_grid
   call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!  This is to allow q=0
   qeq0=.false.
   if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15) then 
     qeq0=.true.
   else
     msg='This routine cannot be used for q/=0'
     MSG_BUG(msg)
   end if

   ia1=1
   do itypat=1,ntypat
  !  ia1,ia2 sets range of loop over atoms:
     ia2=ia1+nattyp(itypat)-1

     ii=0

     do i3=1,n3
       ig3=i3-(i3/id3)*n3-1
       gq3=dble(ig3)+qphon(3)
       gq(3)=gq3
       do i2=1,n2
         if (fftn2_distrib(i2)==mpi_enreg%me_fft) then
           ig2=i2-(i2/id2)*n2-1
           gq2=dble(ig2)+qphon(2)
           gq(2)=gq2

!          Note the lower limit of the next loop
           ii1=1
           if(i3==1 .and. i2==1 .and. qeq0 .and. ig2==0 .and. ig3==0)then
            ii1=2
             ii=ii+1
           end if
           do i1=ii1,n1
             ig1=i1-(i1/id1)*n1-1
             gq1=dble(ig1)+qphon(1)
             gq(1)=gq1
             ii=ii+1
             gsquar=gsq_vl(ig1,ig2,ig3)
!            Skip G**2 outside cutoff:
             if (gsquar<=cutoff) then
               gmag=sqrt(gsquar)

!              Obtain G in cartesian coordinates
               gqc(1)=gprimd(1,1)*gq(1)+gprimd(1,2)*gq(2)+gprimd(1,3)*gq(3)
               gqc(2)=gprimd(2,1)*gq(1)+gprimd(2,2)*gq(2)+gprimd(2,3)*gq(3)
               gqc(3)=gprimd(3,1)*gq(1)+gprimd(3,2)*gq(2)+gprimd(3,3)*gq(3)

!              Evaluate spline fit to get first V(q) and V(q)' and later V(q)'':
               call splfit(qgrid,vion1dq,vlspl(:,:,itypat),1,(/gmag/),vion1,mqgrid,1)
               vion1=vion1/gsquar
               vion1dq=(vion1dq-2.0_dp*gmag*vion1)/gsquar

               call splfit(qgrid,vion1dqdq,vlspl(:,:,itypat),2,(/gmag/),vion1,mqgrid,1)
               vion1dqdq=(vion1dqdq-4.0_dp*gmag*vion1dq-2.0_dp*vion1)/gsquar

!              Assemble structure factor over all atoms of given type:
               sfr=0.0_dp
               sfi=0.0_dp
               do ia=ia1,ia2
                 sfr=sfr+phre_vl(ig1,ig2,ig3,ia)
                 sfi=sfi-phimag_vl(ig1,ig2,ig3,ia)
               end do

               gfact=gqc(beta)*gqc(delta)*gqc(gamma)/gsquar
  
               term1=delbd*gqc(gamma)+delbg*gqc(delta)+deldg*gqc(beta)
               term1=term1-gfact
               term1=term1*vion1dq(1)/gmag

               term2=vion1dqdq(1)*gfact

!              Multiply structure factor times vion derivatives:
               work1(re,ii)=work1(re,ii)+sfr*(term1+term2)
               work1(im,ii)=work1(im,ii)+sfi*(term1+term2)

!              End skip G**2 outside cutoff:
             end if

           end do
         end if
       end do
     end do

     ia1=ia2+1

!  End loop on type of atoms
   end do

!  Set Vloc(G=0)=0:
   work1(re,1)=0.0_dp
   work1(im,1)=0.0_dp

!  Transform back to real space
   call fourdp(cplex,work1,vpsp1dqdq,1,mpi_enreg,nfft,1,ngfft,0)

   xnorm=1.0_dp/ucvol/two_pi
   vpsp1dqdq(1:cplex*nfft)=vpsp1dqdq(1:cplex*nfft)*xnorm

   work1=0.0_dp

!  Calculate the GS density Hartree contribution
   if (opthartdqdq==1) then

     pisqrinv=1.0_dp/pi**2
     
     ii=0
     do i3=1,n3
       ig3=i3-(i3/id3)*n3-1
       gq3=dble(ig3)+qphon(3)
       gq(3)=gq3
       do i2=1,n2
         if (fftn2_distrib(i2)==mpi_enreg%me_fft) then
           ig2=i2-(i2/id2)*n2-1
           gq2=dble(ig2)+qphon(2)
           gq(2)=gq2

!          Note the lower limit of the next loop
           ii1=1
           if(i3==1 .and. i2==1 .and. qeq0 .and. ig2==0 .and. ig3==0)then
             ii1=2
             ii=ii+1
           end if
           do i1=ii1,n1
             ig1=i1-(i1/id1)*n1-1
             gq1=dble(ig1)+qphon(1)
             gq(1)=gq1
             ii=ii+1
             gsquar=gsq_vl(ig1,ig2,ig3)
!            Skip G**2 outside cutoff:
             if (gsquar<=cutoff) then

!              Precalculate quotient of G powers
               uogsquar= 1.0_dp/gsquar

!              Obtain G in cartesian coordinates
               gqc(1)=gprimd(1,1)*gq(1)+gprimd(1,2)*gq(2)+gprimd(1,3)*gq(3)
               gqc(2)=gprimd(2,1)*gq(1)+gprimd(2,2)*gq(2)+gprimd(2,3)*gq(3)
               gqc(3)=gprimd(3,1)*gq(1)+gprimd(3,2)*gq(2)+gprimd(3,3)*gq(3)

               term1=4.0_dp*gqc(beta)*gqc(gamma)*gqc(delta)*uogsquar*uogsquar
               term2=delbd*gqc(gamma)+delbg*gqc(delta)+deldg*gqc(beta) 
               term2=-term2*uogsquar

               work1re=pisqrinv*uogsquar*(term1+term2)
               work1(re,ii)=rhog(re,ii)*work1re
               work1(im,ii)=rhog(im,ii)*work1re

!              End skip G**2 outside cutoff:
             end if

           end do
         end if
       end do
     end do

!    Set V(G=0)=0:
     work1(re,1)=0.0_dp
     work1(im,1)=0.0_dp

!    Transform back to real space
     call fourdp(cplex,work1,vhart1dqdq,1,mpi_enreg,nfft,1,ngfft,0)

!  End the calculation of the Hartree contribution 
   end if

   ABI_DEALLOCATE(work1)

!End the condition of non-electric-field
 end if

 contains

!Real and imaginary parts of phase.
   function phr_vl(x1,y1,x2,y2,x3,y3)
   real(dp) :: phr_vl,x1,x2,x3,y1,y2,y3
   phr_vl=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
 end function phr_vl

   function phi_vl(x1,y1,x2,y2,x3,y3)
   real(dp):: phi_vl,x1,x2,x3,y1,y2,y3
   phi_vl=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
 end function phi_vl

   function ph1_vl(nri,ig1,ia)
   real(dp):: ph1_vl
   integer :: nri,ig1,ia
   ph1_vl=ph1d(nri,ig1+1+n1+(ia-1)*(2*n1+1))
 end function ph1_vl

   function ph2_vl(nri,ig2,ia)
   real(dp):: ph2_vl
   integer :: nri,ig2,ia
   ph2_vl=ph1d(nri,ig2+1+n2+(ia-1)*(2*n2+1)+natom*(2*n1+1))
 end function ph2_vl

   function ph3_vl(nri,ig3,ia)
   real(dp):: ph3_vl
   integer :: nri,ig3,ia
   ph3_vl=ph1d(nri,ig3+1+n3+(ia-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1))
 end function ph3_vl

   function phre_vl(ig1,ig2,ig3,ia)
   real(dp):: phre_vl
   integer :: ig1,ig2,ig3,ia
   phre_vl=phr_vl(ph1_vl(re,ig1,ia),ph1_vl(im,ig1,ia),&
&   ph2_vl(re,ig2,ia),ph2_vl(im,ig2,ia),ph3_vl(re,ig3,ia),ph3_vl(im,ig3,ia))
 end function phre_vl

   function phimag_vl(ig1,ig2,ig3,ia)
   real(dp) :: phimag_vl
   integer :: ig1,ig2,ig3,ia
   phimag_vl=phi_vl(ph1_vl(re,ig1,ia),ph1_vl(im,ig1,ia),&
&   ph2_vl(re,ig2,ia),ph2_vl(im,ig2,ia),ph3_vl(re,ig3,ia),ph3_vl(im,ig3,ia))
 end function phimag_vl

   function gsq_vl(i1,i2,i3)
   real(dp) :: gsq_vl
   integer :: i1,i2,i3
!Define G^2 based on G space metric gmet.
   gsq_vl=dble(i1*i1)*gmet(1,1)+dble(i2*i2)*gmet(2,2)+&
&   dble(i3*i3)*gmet(3,3)+dble(2*i1*i2)*gmet(1,2)+&
&   dble(2*i2*i3)*gmet(2,3)+dble(2*i3*i1)*gmet(3,1)
 end function gsq_vl

end subroutine dfpt_vmetdqdq

end module m_mklocl
!!***
