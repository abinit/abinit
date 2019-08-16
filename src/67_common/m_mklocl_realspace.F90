!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_mklocl_realspace
!! NAME
!!  m_mklocl_realspace
!!
!! FUNCTION
!!   Routines related to the local part of the pseudopotentials.
!!   Computation is done in real space (useful for isolated systems).
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (TRangel, MT, DC)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!  This module could be merged with m_mklocl
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

module m_mklocl_realspace

 use defs_basis
 use defs_datatypes
 use defs_wvltypes
 use m_xmpi
 use m_abicore
 use m_errors

 use defs_abitypes, only : MPI_type
 use m_time,        only : timab
 use m_geometry,    only : xred2xcart
 use m_fft_mesh,    only : mkgrid_fft
 use m_mpinfo,      only : ptabs_fourdp
 use m_pawtab,      only : pawtab_type
 use m_paw_numeric, only : paw_splint, paw_splint_der
 use m_psolver,     only : psolver_hartree, psolver_kernel
 use m_abi2big,     only : wvl_rhov_abi2big
 use m_wvl_wfs,     only : derf_ab
 use m_fft,         only : fourdp

 implicit none

 private
!!***

 public :: mklocl_realspace
 public :: mklocl_wavelets
!!***

contains
!!***

!!****f* ABINIT/mklocl_realspace
!! NAME
!!  mklocl_realspace
!!
!! FUNCTION
!! This method is equivalent to mklocl_recipspace except that
!! it uses real space pseudo-potentials. It is useful for isolated
!! systems. Then the option 3 and 4 are not available for this implementation.
!!
!! Optionally compute:
!!  option=1 : local ionic potential throughout unit cell
!!  option=2 : contribution of local ionic potential to E gradient wrt xred
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in unit cell.
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms.
!!  option= (see above)
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfft)=electron density rho(G) (electrons/$\textrm{Bohr}^3$)
!!  rhor(nfft,nspden)=electron density in electrons/bohr**3.
!!    (needed if option==2 or if option==4)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  ucvol=unit cell volume ($\textrm{Bohr}^3$).
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  (if option==1) vpsp(nfft)=local crystal pseudopotential in real space.
!!  (if option==2) grtn(3,natom)=grads of Etot wrt tn. These gradients are in
!!                 reduced coordinates. Multiply them by rprimd to get
!!                 gradients in cartesian coordinates.
!!
!! PARENTS
!!      mklocl
!!
!! CHILDREN
!!
!! SOURCE

subroutine mklocl_realspace(grtn,icoulomb,mpi_enreg,natom,nattyp,nfft,ngfft,nscforder, &
&                           nspden,ntypat,option,pawtab,psps,rhog,rhor,rprimd,typat,&
&                           ucvol,usewvl,vpsp,xred)

#if defined HAVE_BIGDFT
 use BigDFT_API,    only : coulomb_operator,deallocate_coulomb_operator
 use defs_PSolver
#else
 use defs_wvltypes, only : coulomb_operator
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,nspden,ntypat,option
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
 type(pawtab_type),intent(in)  :: pawtab(ntypat*psps%usepaw)
!arrays
 integer,intent(in)  :: icoulomb,nscforder,usewvl
 integer,intent(in)  :: nattyp(ntypat),ngfft(18),typat(natom)
 real(dp),intent(in) :: rhog(2,nfft)
 real(dp),intent(in) :: rhor(nfft,nspden),rprimd(3,3)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(out) :: grtn(3,natom),vpsp(nfft)

!Local variables-------------------------------
 character(len=1) :: geocode
  !testing variables
 !scalars
 integer,parameter :: nStep=2
 integer :: comm_fft,countParSeconde,i1,i2,i3
 integer :: ia,ia1,ia2,igeo,ii,ind,itypat,ix,iy,iz,jj
 integer :: kk,me_fft,n1,n2,n3,n3d,n_interpol
 integer :: nproc_fft,tpsStart,tpsStop
 real(dp),parameter :: min_rho_value=1.0d-12
 real(dp) :: aa,bb,cc,dd,delta,deltaV,dr,dr2div6,invdr,r,vol_interpol,x,y,z,hgx,hgy,hgz,entmp
 logical,parameter :: customRho=.false.,finiteDiff=.false.,testing=.false.
 logical :: doIt
 character(len=500) :: message
!arrays
 integer :: ngfft_interpol(18)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: coord(3),coordXYZ(3),refValue(3),tsec(2)
 real(dp),allocatable :: coordCart_interpol(:,:),coordRed_interpol(:,:)
 real(dp),allocatable :: gridcart(:,:)
 real(dp),allocatable :: grtn_cart_interpol(:,:),grtn_diff(:,:)
 real(dp),allocatable :: rhog_interpol(:,:),rhog_testing(:,:),rhor_interpol(:)
 real(dp),allocatable :: rhor_testing(:),rhor_work(:),xcart(:,:),vhartr(:),gxyz(:,:)
 type(coulomb_operator):: kernel

! *************************************************************************

!Keep track of total time spent here
 if (option==2) then
   call timab(72,1,tsec)
 end if

!Several constants (FFT sizes and parallelism)
 n1 = ngfft(1) ; n2 = ngfft(2) ; n3 = ngfft(3)
 nproc_fft = ngfft(10) ;  me_fft = ngfft(11)
 n3d = ngfft(13)          !for parallel runs
 if (nproc_fft==1) n3d=n3  !for serial runs
 comm_fft=mpi_enreg%comm_fft
 if(me_fft /= mpi_enreg%me_fft .or. nproc_fft /= mpi_enreg%nproc_fft) then
   MSG_BUG("mpi_enreg%x_fft not equal to the corresponding values in ngfft")
 end if

!Conditions for periodicity in the three directions
 geocode='P'
 if (icoulomb==1) geocode='F'
 if (icoulomb==2) geocode='S'

!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!Store xcart for each atom
 ABI_ALLOCATE(xcart,(3, natom))
 call xred2xcart(natom, rprimd, xcart, xred)
!Store cartesian coordinates for each grid points
 ABI_ALLOCATE(gridcart,(3, nfft))
 call mkgrid_fft(ffti3_local,fftn3_distrib,gridcart,nfft,ngfft,rprimd)

!Check whether all the PSP considered are of type GTH-HGH or PAW
 doIt=.true.
 do ii=1,psps%npsp
   doIt=doIt .and.&
   (psps%pspcod(ii)==2.or.psps%pspcod(ii)==3.or.psps%pspcod(ii)==10.or.psps%pspcod(ii)==7)
 end do

!HGH-GTH/PAW treatment presumably starts here
 if (doIt) then

!  Definition of the grid spacings as in the kernel routine
   hgx = rprimd(1,1)/(n1)
   hgy = rprimd(2,2)/(n2)
   hgz = rprimd(3,3)/(n3)


!----------------------------------------------------------------------
! ----- Option 1: compute local ionic potential                   -----
!----------------------------------------------------------------------
   if (option==1) then

     call psolver_kernel( (/ hgx, hgy, hgz /), 2, icoulomb, me_fft, kernel, comm_fft, &
&     (/n1,n2,n3/), nproc_fft, nscforder)

     call createIonicPotential_new(fftn3_distrib,ffti3_local,&
&     geocode,me_fft, nproc_fft, natom, &
&     ntypat, typat, psps%gth_params%psppar, &
&     int(psps%ziontypat), xcart,gridcart, hgx,hgy,hgz, &
&     n1,n2,n3d,n3, kernel, vpsp, comm_fft,pawtab,psps%usepaw)

!----------------------------------------------------------------------
! ----- Option 2: compute forces induced by local ionic potential -----
!----------------------------------------------------------------------
   else if (option == 2) then

!    Compute Hartree potential from rhor
     ABI_ALLOCATE(vhartr,(nfft))
     call psolver_hartree(entmp, (/ hgx, hgy, hgz /), icoulomb, me_fft, comm_fft, nfft, &
&     (/n1,n2,n3/), nproc_fft, nscforder, nspden, rhor, vhartr, usewvl)

!    Allocate temporary array for forces
     ABI_ALLOCATE(gxyz,(3, natom))

!    Calculate local part of the forces grtn (inspired from BigDFT routine)
     call local_forces_new(fftn3_distrib,ffti3_local,geocode,me_fft, ntypat, natom, &
&     typat, xcart, gridcart, psps%gth_params%psppar, int(psps%ziontypat), &
&     hgx,hgy,hgz, n1,n2,n3,n3d, rhor,vhartr, gxyz, pawtab,psps%usepaw)

!    Forces should be in reduced coordinates.
     do ia = 1, natom, 1
       do igeo = 1, 3, 1
         grtn(igeo, ia) = - rprimd(1, igeo) * gxyz(1, ia) &
&         - rprimd(2, igeo) * gxyz(2, ia) &
&         - rprimd(3, igeo) * gxyz(3, ia)
       end do
     end do

!    Deallocate local variables
     ABI_DEALLOCATE(vhartr)
     ABI_DEALLOCATE(gxyz)
   end if

!----------------------------------------------------------------------
! ----- Section for the non-HGH/GTH/PAW pseudopotentials (testing) ----
!----------------------------------------------------------------------
 else

   if (testing) then
     call system_clock(count_rate = countParSeconde)
     call system_clock(tpsStart, count_rate = countParSeconde)
   end if

!  dr is the r step in the sampling psps%vlspl
   dr = psps%qgrid_vl(2)
   invdr = 1._dp / dr
   dr2div6 = dr * dr / 6._dp

   if (option == 1) then
!    Set 0 in vpsp before summing
     vpsp(:) = 0._dp
   else if (option == 2) then
!    Allocate array to store cartesian gradient computed with
!    an interpolation of rhor
     ABI_ALLOCATE(grtn_cart_interpol,(3, natom))
     grtn_cart_interpol(:, :) = 0._dp

     n_interpol = nStep ** 3
     ABI_ALLOCATE(coordRed_interpol,(3, nStep ** 3))
     ABI_ALLOCATE(coordCart_interpol,(3, nStep ** 3))

     if (testing .and. customRho) then
!      Use a custom rho instead of the self-consistent one.
       ABI_ALLOCATE(rhor_testing,(nfft))
       ABI_ALLOCATE(rhog_testing,(2, nfft))
     end if

     ABI_ALLOCATE(rhor_interpol,(nfft * n_interpol))
     ABI_ALLOCATE(rhor_work,(nfft * n_interpol))
     ABI_ALLOCATE(rhog_interpol,(2, nfft * n_interpol))

     if (testing .and. customRho) then
!      Testing only, changing rho with a centered gaussian
       do ii = 1, nfft, 1
!        using the position of the first atom as center.
         r = (gridcart(1, ii) - xcart(1, 1)) ** 2 + &
&         (gridcart(2, ii) - xcart(2, 1)) ** 2 + &
&         (gridcart(3, ii) - xcart(3, 1)) ** 2
         rhor_testing(ii) = exp(-r/4._dp)
       end do
!      Testing only, compute rhog_testing from rhor_testing
       call fourdp(1,rhog_testing,rhor_testing,-1,mpi_enreg,nfft,1,ngfft,0)
     end if

!    Compute the interpolation of rho, using a fourier transform
     rhog_interpol(:, :) = 0._dp
     ii = 0
     do i3 = 1, n3, 1
       if (i3 <= n3 / 2) then
         iz = i3
       else
         iz = n3 * nStep - n3 + i3
       end if
       do i2 = 1, n2, 1
         if (i2 <= n2 / 2) then
           iy = i2
         else
           iy = n2 * nStep - n2 + i2
         end if
         do i1 = 1, n1, 1
           ii = ii + 1
           if (i1 <= n1 / 2) then
             ix = i1
           else
             ix = n1 * nStep - n1 + i1
           end if
           jj = (iz - 1) * n2 * n1 * nStep ** 2 + (iy - 1) * n3 * nStep + ix
           if (testing .and. customRho) then
             rhog_interpol(:, jj) = rhog_testing(:, ii)
           else
             rhog_interpol(:, jj) = rhog(:, ii)
           end if
         end do
       end do
     end do

!    Compute the interpolation of rho from the Fourier transformation
     ngfft_interpol(:) = ngfft(:)
     ngfft_interpol(1:3) = (/ n1 * nStep, n2 * nStep, n3 * nStep /)
     ngfft_interpol(4:6) = (/ n1 * nStep + 1, n2 * nStep + 1, n3 * nStep /)
     call fourdp(1,rhog_interpol,rhor_work,1,mpi_enreg,nfft*n_interpol,1,ngfft_interpol,0)

!    Reorder rhor_interpol to be able to read it linearly
     jj = 0
     do i3 = 1, n3, 1
       do i2 = 1, n2, 1
         do i1 = 1, n1, 1
           do iz = 1, nStep, 1
             do iy = 1, nStep, 1
               do ix = 1, nStep, 1
                 jj = jj + 1
                 kk = ((i3 - 1) * nStep + iz - 1) ! z coordinate in the interpolated grid
                 kk = kk * n1 * n2 * nStep ** 2
                 kk = kk + ((i2 - 1) * nStep + iy - 1) * n1 * nStep ! adding y coordinate
                 kk = kk + (i1 - 1) * nStep + ix ! adding x coordinate
                 rhor_interpol(jj) = rhor_work(kk)
               end do
             end do
           end do
         end do
       end do
     end do
     ABI_DEALLOCATE(rhor_work)

!    Compute grid access in the interpolated volume
     ii = 0
     do iz = 1, nStep, 1
       z = real(iz - 1, dp) / real(nStep, dp)
       do iy = 1, nStep, 1
         y = real(iy - 1, dp) / real(nStep, dp)
         do ix = 1, nStep, 1
           x = real(ix - 1, dp) / real(nStep, dp)
           ii = ii + 1
           coordRed_interpol(:, ii) = (/ x, y, z /)
!          Assuming orthogonal box (should be change later)
           coordCart_interpol(:, ii) = (/ x * rprimd(1, 1) / real(n1, dp), &
&           y * rprimd(2, 2) / real(n2, dp), &
&           z * rprimd(3, 3) / real(n3, dp) /)
         end do
       end do
     end do

     vol_interpol = 1._dp / real(nStep, dp) ** 3
!    Compute the coordinates (integer) of each atom and deduce
!    the max extens of the integral summation.
!    !!  do ia = 1, natom, 1
!    !!   coordAtom(1, ia) = int(xred(1, ia) * n1) + 1
!    !!   coordAtom(2, ia) = int(xred(2, ia) * n2) + 1
!    !!   coordAtom(3, ia) = int(xred(3, ia) * n3) + 1
!    !!  end do
   end if

   if (testing .and. option == 2) then
     call system_clock(tpsStop, count_rate = countParSeconde)
     write(std_out,*) "Tps : ", real(tpsStop - tpsStart) / real(countParSeconde)
   end if

   ia1=1
   do itypat = 1, ntypat, 1
!    ia1,ia2 sets range of loop over atoms:
     ia2 = ia1 + nattyp(itypat) - 1

     do ii = 1, nfft, 1
       do ia = ia1, ia2, 1
         if (option == 1) then
!          Compute the potential
!          r is the distance between grid point and atom
           r = sqrt((gridcart(1, ii) - xcart(1, ia)) ** 2 + &
&           (gridcart(2, ii) - xcart(2, ia)) ** 2 + &
&           (gridcart(3, ii) - xcart(3, ia)) ** 2)

!          Coefficients needed to compute the spline.
           jj = int(r * invdr) + 1
           if (jj > psps%mqgrid_vl - 2) then
             write(message, '(3a,i0,a,i0,a,a)' )&
&             '  pseudo-potential local part sampling is not wide enough', ch10, &
&             '  want to access position ', jj, ' whereas mqgrid_vl = ', psps%mqgrid_vl, ch10, &
&             '  Action : no idea, contact developpers...'
             MSG_ERROR(message)
           end if
           delta = r - psps%qgrid_vl(jj)
           bb = delta * invdr
           aa = 1._dp - bb
           cc = aa * (aa ** 2 - 1._dp) * dr2div6
           dd = bb * (bb ** 2 - 1._dp) * dr2div6

!          compute V(r) from the spline, jj and jj + 1 is braketting r in
!          the sampling
           deltaV = aa * psps%vlspl(jj, 1, itypat) + bb * psps%vlspl(jj + 1, 1, itypat) + &
&           cc * psps%vlspl(jj, 2, itypat) + dd * psps%vlspl(jj + 1, 2, itypat)
!          Add on grid point ii the contribution of atom ia
           vpsp(ii) = vpsp(ii) + deltaV
         else if (option == 2) then
!          Compute the forces, as gradient of energy (V(r).rho(r))

!          Testing only - reference points
           if (.false.) then
!            r is the distance between grid point and atom
             r = sqrt((gridcart(1, ii) - xcart(1, ia)) ** 2 + &
&             (gridcart(2, ii) - xcart(2, ia)) ** 2 + &
&             (gridcart(3, ii) - xcart(3, ia)) ** 2)

!            Coefficients needed to compute the spline.
             jj = int(r * invdr) + 1
             delta = r - psps%qgrid_vl(jj)
             bb = delta * invdr
             aa = 1._dp - bb
             cc = aa * (aa ** 2 - 1._dp) * dr2div6
             dd = bb * (bb ** 2 - 1._dp) * dr2div6

!            When mesh position is on a node, forces are null.
             if (r /= 0._dp) then
!              This value deltaV is the first derivative of V(r) taken at r.
               deltaV = aa * psps%dvlspl(jj, 1, itypat) + bb * psps%dvlspl(jj + 1, 1, itypat) + &
&               cc * psps%dvlspl(jj, 2, itypat) + dd * psps%dvlspl(jj + 1, 2, itypat)
!              We multiply by rho(r) to have an energy.
               deltaV = deltaV * rhor(ii, 1) / r
               refValue(:) = - deltaV * (gridcart(:, ii) - xcart(:, ia))
               grtn_cart_interpol(:, ia) = grtn_cart_interpol(:, ia) + refValue(:)
             end if
           end if

!          Compute the interpolation for the point ii
           ind = (ii - 1) * n_interpol
           do kk = 1, n_interpol, 1
             ind = ind + 1

             if (rhor_interpol(ind) > min_rho_value) then
!              Assume orthogonal box...
               coordXYZ(1) = gridcart(1, ii) - xcart(1, ia) + coordCart_interpol(1, kk)
               coordXYZ(2) = gridcart(2, ii) - xcart(2, ia) + coordCart_interpol(2, kk)
               coordXYZ(3) = gridcart(3, ii) - xcart(3, ia) + coordCart_interpol(3, kk)
               r = coordXYZ(1) ** 2 + coordXYZ(2) ** 2 + coordXYZ(3) ** 2

               if (r /= 0._dp) then
                 r = sqrt(r)
!                Coefficients needed to compute the spline.
                 jj = int(r * invdr) + 1
                 delta = r - psps%qgrid_vl(jj)
                 bb = delta * invdr
                 aa = 1._dp - bb
                 cc = aa * (aa ** 2 - 1._dp) * dr2div6
                 dd = bb * (bb ** 2 - 1._dp) * dr2div6
                 deltaV = aa * psps%dvlspl(jj, 1, itypat) + &
&                 bb * psps%dvlspl(jj + 1, 1, itypat) + &
&                 cc * psps%dvlspl(jj, 2, itypat) + &
&                 dd * psps%dvlspl(jj + 1, 2, itypat)
                 deltaV = deltaV * rhor_interpol(ind) / r
                 grtn_cart_interpol(1, ia) = grtn_cart_interpol(1, ia) - deltaV * coordXYZ(1)
                 grtn_cart_interpol(2, ia) = grtn_cart_interpol(2, ia) - deltaV * coordXYZ(2)
                 grtn_cart_interpol(3, ia) = grtn_cart_interpol(3, ia) - deltaV * coordXYZ(3)
!                do igeo = 1, 3, 1
!                grtn_cart_interpol(igeo, ia) = grtn_cart_interpol(igeo, ia) - deltaV * coordXYZ(igeo)
!                end do
               end if
             end if
           end do

!          =============
!          Testing only
!          =============
!          use of finite differences
           if (finiteDiff) then
             do igeo = 1, 3, 1
               coord(:) = 0._dp
               coord(igeo) = dr / 2.0_dp
               r = sqrt((gridcart(1, ii) - xcart(1, ia) + coord(1)) ** 2 + &
&               (gridcart(2, ii) - xcart(2, ia) + coord(2)) ** 2 + &
&               (gridcart(3, ii) - xcart(3, ia) + coord(3)) ** 2)

!              Coefficients needed to compute the spline.
               jj = int(r * invdr) + 1
               delta = r - psps%qgrid_vl(jj)
               bb = delta * invdr
               aa = 1._dp - bb
               cc = aa * (aa ** 2 - 1._dp) * dr2div6
               dd = bb * (bb ** 2 - 1._dp) * dr2div6

               deltaV = aa * psps%vlspl(jj, 1, itypat) + bb * psps%vlspl(jj + 1, 1, itypat) + &
&               cc * psps%vlspl(jj, 2, itypat) + dd * psps%vlspl(jj + 1, 2, itypat)


               coord(:) = 0._dp
               coord(igeo) = -dr / 2.0_dp
               r = sqrt((gridcart(1, ii) - xcart(1, ia) + coord(1)) ** 2 + &
&               (gridcart(2, ii) - xcart(2, ia) + coord(2)) ** 2 + &
&               (gridcart(3, ii) - xcart(3, ia) + coord(3)) ** 2)

!              Coefficients needed to compute the spline.
               jj = int(r * invdr) + 1
               delta = r - psps%qgrid_vl(jj)
               bb = delta * invdr
               aa = 1._dp - bb
               cc = aa * (aa ** 2 - 1._dp) * dr2div6
               dd = bb * (bb ** 2 - 1._dp) * dr2div6

               deltaV = deltaV - (aa * psps%vlspl(jj, 1, itypat) + &
&               bb * psps%vlspl(jj + 1, 1, itypat) + &
&               cc * psps%vlspl(jj, 2, itypat) + &
&               dd * psps%vlspl(jj + 1, 2, itypat))
               grtn_diff(igeo, ia) = grtn_diff(igeo, ia) - deltaV * rhor(ii, 1) / dr
             end do
           end if
!          =============
!          Testing only
!          =============

         end if
       end do
!      End loop over atoms of type itypat
     end do
!    End loop over real space grid points

     ia1 = ia2 + 1
   end do
!  End loop over type of atoms

   if(option==2)then
!    multiply the forces by the volume of a single box mesh.
     grtn_cart_interpol(:, :) = grtn_cart_interpol(:, :) * &
&     ucvol / real(n1 * n2 * n3, dp) * vol_interpol
!    Transform cartesian forces to reduce coordinates
     do ia = 1, natom, 1
       do igeo = 1, 3, 1
         grtn(igeo, ia) = rprimd(1, igeo) * grtn_cart_interpol(1, ia) + &
&         rprimd(2, igeo) * grtn_cart_interpol(2, ia) + &
&         rprimd(3, igeo) * grtn_cart_interpol(3, ia)
       end do
     end do
     ABI_DEALLOCATE(rhor_interpol)
     ABI_DEALLOCATE(rhog_interpol)
     ABI_DEALLOCATE(coordRed_interpol)
     ABI_DEALLOCATE(coordCart_interpol)
     if (testing .and. customRho) then
       ABI_DEALLOCATE(rhor_testing)
       ABI_DEALLOCATE(rhog_testing)
     end if

     if (testing) then
       call system_clock(tpsStop, count_rate = countParSeconde)
       write(std_out,*) "Tps : ", real(tpsStop - tpsStart) / real(countParSeconde)
       write(std_out,*) grtn_cart_interpol
       MSG_ERROR("Testing section!")
     end if

   end if

!-------------------------
 end if ! GTH/HGH/PAW psps

!Release temporary memory
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(gridcart)

!Close timing counters
 if (option==2)then
   call timab(72,2,tsec)
 end if

end subroutine mklocl_realspace
!!***

!----------------------------------------------------------------------

!!****f* mklocl_realspace/createIonicPotential_new
!! NAME
!!  createIonicPotential_new
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      mklocl_realspace
!!
!! CHILDREN
!!
!! SOURCE

subroutine createIonicPotential_new(fftn3_distrib,ffti3_local,geocode,iproc,&
&  nproc,nat,ntypes,iatype,psppar,nelpsp,rxyz,gridcart,&
&  hxh,hyh,hzh,n1i,n2i,n3d,n3i,kernel,pot_ion,spaceworld,pawtab,usepaw)

 use defs_wvltypes, only : coulomb_operator

!Arguments -------------------------------
!scalars
 integer, intent(in) :: iproc,nproc,ntypes,nat,n1i,n2i,n3i,n3d,spaceworld,usepaw
 real(dp), intent(in) :: hxh,hyh,hzh
 character(len=1), intent(in) :: geocode
 type(coulomb_operator), intent(in) :: kernel
!arrays
 integer, dimension(nat), intent(in) :: iatype
 integer, dimension(ntypes), intent(in) :: nelpsp
 integer, dimension(*), intent(in) ::fftn3_distrib,ffti3_local
 real(dp), dimension(3,n1i*n2i*n3d), intent(in) :: gridcart
 real(dp), dimension(0:4,0:6,ntypes), intent(in) :: psppar
 real(dp), dimension(3,nat), intent(in) :: rxyz
 real(dp), dimension(*), intent(inout) :: pot_ion
 type(pawtab_type),intent(in)  :: pawtab(ntypes*usepaw)

!Local variables -------------------------
#if defined HAVE_BIGDFT
!scalars
 integer :: iat,i1,i2,i3,j1,j2,j3,isx,isy,isz,iex,iey,iez,ierr,ityp
 integer :: ind,nloc,iloc,i3loc,msz
 logical :: gox,goy,goz,perx,pery,perz
 real(dp) :: arg,charge,cutoff,ehart,eexcu,rholeaked,rholeaked_tot,rloc
 real(dp) :: rx,ry,rz,rzero,r2,tt,tt_tot,vexcu,vhgh,x,xp,y,z
!arrays
 real(dp) :: rr(1),vpaw(1)
 real(dp),pointer :: rad(:),vloc(:),d2vloc(:)
#endif

! *********************************************************************

#if defined HAVE_BIGDFT

 if(nproc<0)then
   MSG_ERROR('nproc should not be negative')
 end if

!Ionic charge (must be calculated for the PS active processes)
 rholeaked=0._dp
!Ionic energy (can be calculated for all the processors)

!here we should insert the calculation of the ewald energy for the periodic BC case
!!!  eion=0._dp
!!!  do iat=1,nat
!!!     ityp=iatype(iat)
!!!     rx=rxyz(1,iat)
!!!     ry=rxyz(2,iat)
!!!     rz=rxyz(3,iat)
!!!     !    ion-ion interaction
!!!     do jat=1,iat-1
!!!        dist=sqrt( (rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2 )
!!!        jtyp=iatype(jat)
!!!        eion=eion+real(nelpsp(jtyp)*nelpsp(ityp),kind=dp)/dist
!!!     enddo
!!!  end do
!!!  if (iproc.eq.0) write(std_out,'(1x,a,1pe22.14)') 'ion-ion interaction energy',eion

!Creates charge density arising from the ionic PSP cores
!the n3pi dimension indicates the number of planes trated by each processor in the FFT parallelisation
!for a plane wave treatment this value depends on whether the direct space is divided in planes or not
!I don't know this variable, which in the future must be inserted at the place of n3pi (LG)
!if n3pi=0 this means that the processors doesn't calculate anything
!if (n3pi >0 ) then

!conditions for periodicity in the three directions
 perx=(geocode /= 'F')
 pery=(geocode == 'P')
 perz=(geocode /= 'F')

!this initialise the array to zero, it will work only if bigdft library is enabled
 pot_ion(1:n1i*n2i*n3d)=zero

 do iat=1,nat
   ityp=iatype(iat)
   rx=rxyz(1,iat)
   ry=rxyz(2,iat)
   rz=rxyz(3,iat)

   rloc=psppar(0,0,ityp)
   charge=real(nelpsp(ityp),kind=dp)/(2._dp*pi*sqrt(2._dp*pi)*rloc**3)
   cutoff=10._dp*rloc

   isx=floor((rx-cutoff)/hxh)
   isy=floor((ry-cutoff)/hyh)
   isz=floor((rz-cutoff)/hzh)

   iex=ceiling((rx+cutoff)/hxh)
   iey=ceiling((ry+cutoff)/hyh)
   iez=ceiling((rz+cutoff)/hzh)

!  Calculate Ionic Density
!  using HGH parameters.
!  Eq. 1.104, T. Deutsch and L. Genovese, JDN. 12, 2011
   do i3=isz,iez
     z=real(i3,kind=dp)*hzh-rz
     call ind_positions_mklocl(perz,i3,n3i,j3,goz)
     if(fftn3_distrib(j3)==iproc) then
       i3loc=ffti3_local(j3)
       do i2=isy,iey
         y=real(i2,kind=dp)*hyh-ry
         call ind_positions_mklocl(pery,i2,n2i,j2,goy)
         do i1=isx,iex
           x=real(i1,kind=dp)*hxh-rx
           call ind_positions_mklocl(perx,i1,n1i,j1,gox)
           r2=x**2+y**2+z**2
           if (goz  .and. goy  .and. gox ) then
             ind=j1+(j2-1)*n1i+(i3loc-1)*n1i*n2i
             r2=(gridcart(1,ind)-rx)**2+(gridcart(2,ind)-ry)**2+(gridcart(3,ind)-rz)**2
           end if
           arg=r2/rloc**2
           xp=exp(-.5d0*arg)
           if (goz  .and. goy  .and. gox ) then
             pot_ion(ind)=pot_ion(ind)-xp*charge
           else
             rholeaked=rholeaked+xp*charge
           end if
         end do
       end do
     end if
   end do

 end do

!Check
 tt=0._dp
 do j3= 1,n3d
   do i2= 1,n2i
     do i1= 1,n1i
       ind=i1+(i2-1)*n1i+(j3-1)*n1i*n2i
       tt=tt+pot_ion(ind)
     end do
   end do
 end do

 tt=tt*hxh*hyh*hzh
 rholeaked=rholeaked*hxh*hyh*hzh

 call xmpi_sum(tt,tt_tot,spaceworld,ierr)
 call xmpi_sum(rholeaked,rholeaked_tot,spaceworld,ierr)

 if (iproc.eq.0) then
   write(std_out,'(1x,a,f26.12,2x,1pe10.3)') &
&   'total ionic charge, leaked charge ',tt_tot,rholeaked_tot
 end if

!Here the value of the datacode must be kept fixed
!there can be some problems when running this stuff in parallel,
!  if the ionic potential distribution does not agree with the
!  plane distribution which is supposed to hold for the Poisson Solver
 call psolver(geocode,'D',iproc,nproc,n1i,n2i,n3i,0,hxh,hyh,hzh,&
& pot_ion,kernel%kernel,pot_ion,ehart,eexcu,vexcu,0._dp,.false.,1)

!Add the remaining short-range local terms
 do iat=1,nat
   ityp=iatype(iat)

   rx=rxyz(1,iat)
   ry=rxyz(2,iat)
   rz=rxyz(3,iat)

!  determine number of local terms
   rloc=psppar(0,0,ityp)
   cutoff=10._dp*rloc
   charge=real(nelpsp(ityp),kind=dp)

!  determine number of local terms (HGH pot)
   nloc=0
   do iloc=1,4
     if (psppar(0,iloc,ityp).ne.0._dp) nloc=iloc
   end do

!  PAW specifics
   if (usepaw==1) then
     msz=pawtab(ityp)%wvl%rholoc%msz
     rad    => pawtab(ityp)%wvl%rholoc%rad(1:msz)
     vloc   => pawtab(ityp)%wvl%rholoc%d(1:msz,3)
     d2vloc => pawtab(ityp)%wvl%rholoc%d(1:msz,4)
     rzero=rad(1);if (rzero<=1.d-10) rzero=rad(2)
   end if

   isx=floor((rx-cutoff)/hxh)
   isy=floor((ry-cutoff)/hyh)
   isz=floor((rz-cutoff)/hzh)

   iex=ceiling((rx+cutoff)/hxh)
   iey=ceiling((ry+cutoff)/hyh)
   iez=ceiling((rz+cutoff)/hzh)

   do i3=isz,iez
     z=real(i3,kind=dp)*hzh-rz
     call ind_positions_mklocl(perz,i3,n3i,j3,goz)
     if(fftn3_distrib(j3) == iproc .and. goz) then !MPI
       i3loc=ffti3_local(j3)
       if (goz) then
         do i2=isy,iey
           y=real(i2,kind=dp)*hyh-ry
           call ind_positions_mklocl(pery,i2,n2i,j2,goy)
           if (goy) then
             do i1=isx,iex
               x=real(i1,kind=dp)*hxh-rx
               call ind_positions_mklocl(perx,i1,n1i,j1,gox)
               if (gox) then
                 ind=j1+(j2-1)*n1i+(i3loc-1)*n1i*n2i
                 r2=(gridcart(1,ind)-rx)**2+(gridcart(2,ind)-ry)**2+(gridcart(3,ind)-rz)**2
!                r2=x**2+y**2+z**2

!                HGH: V_S=gaussian potential of Eq. (9) in JCP 129, 014109(2008)
                 if (usepaw==0) then
                   if (nloc /= 0) then
                     arg=r2/rloc**2
                     xp=exp(-.5d0*arg)
                     tt=psppar(0,nloc,ityp)
                     do iloc=nloc-1,1,-1
                       tt=arg*tt+psppar(0,iloc,ityp)
                     end do
                     pot_ion(ind)=pot_ion(ind)+xp*tt
                   end if

!                PAW: V_PAW-V_L^HGH
                 else
                   rr(1)=sqrt(r2)
                   if (rr(1)>=rzero) then
                     call paw_splint(msz,rad,vloc,d2vloc,1,rr,vpaw)
                     call calcVloc_mklocl(vhgh,rr(1),rloc,charge)
                     pot_ion(ind)=pot_ion(ind)+vpaw(1)-vhgh
                   else
                     pot_ion(ind)=pot_ion(ind)+vloc_zero_mklocl(charge,rloc,msz,rad,vloc,d2vloc)
                   end if
                 end if

               end if
             end do
           end if
         end do
       end if
     end if
   end do

 end do  !iat
#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) geocode,iproc,nproc,ntypes,nat,n1i,n2i,n3i,n3d,spaceworld,usepaw,&
& hxh,hyh,hzh,iatype(1),nelpsp(1),fftn3_distrib(1),ffti3_local(1),gridcart(1,1),psppar(1,1,1),&
& rxyz(1,1),pot_ion(1),pawtab(1)%mesh_size,kernel%co
#endif

 CONTAINS
!!***

!----------------------------------------------------------------------

!!****f* mklocl_realspace/calcVloc_mklocl
!! NAME
!!  calcVloc_mklocl
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      mklocl_realspace
!!
!! CHILDREN
!!
!! SOURCE

subroutine calcVloc_mklocl(yy,xx,rloc,Z)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in)  :: xx,rloc,Z
 real(dp),intent(out) :: yy

!Local variables-------------------------------
 !scalars
 real(dp):: arg,tt

! *************************************************************************

   arg=xx/(sqrt(2.0)*rloc)
   call derf_ab(tt,arg)
   yy=-Z/xx*tt

 end subroutine calcVloc_mklocl
!!***

!----------------------------------------------------------------------

!!****f* mklocl_realspace/vloc_zero_mklocl
!! NAME
!!  vloc_zero_mklocl
!!
!! FUNCTION
!!  Use a quadratic interpolation to get limit of Vloc(x) at x->0
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

function vloc_zero_mklocl(charge,rloc,msz,rad,vloc,d2vloc)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: msz
 real(dp) :: vloc_zero_mklocl
 real(dp),intent(in)  :: charge,rloc
!arrays
 real(dp) :: rad(msz),vloc(msz),d2vloc(msz)

!Local variables-------------------------------
!scalars
 real(dp) :: y1,y2,y3,zz=0._dp
!arrays
 real(dp) :: ll(3),xx(3),yy(3)

! *************************************************************************

!Select 3 points x1,x2,x3 near 0
   if (rad(1)>1.d-10) then
     xx(1:3)=rad(1:3)
   else
     xx(1:3)=rad(2:4)
   end if

!Find the corresponding values of y=(V^PAW(x)-V^HGH(x))/x
   call paw_splint(msz,rad,vloc,d2vloc,3,xx,yy)
   call calcVloc_mklocl(y1,xx(1),rloc,charge)
   call calcVloc_mklocl(y2,xx(2),rloc,charge)
   call calcVloc_mklocl(y3,xx(3),rloc,charge)
   yy(1)= yy(1)-y1
   yy(2)= yy(2)-y2
   yy(3)= yy(3)-y3

!Find a polynomial of the form (z=0):
!P(z) = y1.L1(z) + y2.L2(z) + y3.L3(z)

!L1(z) = (z-x2)(z-x3)/((x1-x2)(x1-x3))
   ll(1)=(zz-xx(2))*(zz-xx(3))/((xx(1)-xx(2))*(xx(1)-xx(3)))
!L2(z) = (z-x1)(z-x3)/((x2-x1)(x2-x3))
   ll(2)=(zz-xx(1))*(zz-xx(3))/((xx(2)-xx(1))*(xx(2)-xx(3)))
!L3(z) = (z-x1)(z-x2)/((x3-x1)(x3-x2))
   ll(3)=(zz-xx(1))*(zz-xx(2))/((xx(3)-xx(1))*(xx(3)-xx(2)))

   vloc_zero_mklocl=yy(1)*ll(1)+yy(2)*ll(2)+yy(3)*ll(3)

 end function vloc_zero_mklocl
!!***

end subroutine createIonicPotential_new
!!***

!----------------------------------------------------------------------

!!****f* mklocl_realspace/local_forces_new
!! NAME
!!  local_forces_new
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      mklocl_realspace
!!
!! CHILDREN
!!
!! SOURCE

subroutine local_forces_new(fftn3_distrib,ffti3_local,&
     geocode,iproc,ntypes,nat,iatype,rxyz,gridcart,psppar,nelpsp,hxh,hyh,hzh,&
     n1,n2,n3,n3d,rho,pot,floc,pawtab,usepaw)


!Arguments -------------------------------
!scalars
 integer, intent(in) :: iproc,ntypes,nat,n1,n2,n3,n3d,usepaw
 character(len=1), intent(in) :: geocode
 real(dp), intent(in) :: hxh,hyh,hzh
!arrays
 integer, dimension(*), intent(in) ::fftn3_distrib,ffti3_local
 integer, dimension(nat), intent(in) :: iatype
 integer, dimension(ntypes), intent(in) :: nelpsp
 real(dp), dimension(3,n1*n2*n3d), intent(in) :: gridcart
 real(dp), dimension(0:4,0:6,ntypes), intent(in) :: psppar
 real(dp), dimension(3,nat), intent(in) :: rxyz
 real(dp), dimension(*), intent(in) :: rho,pot
 real(dp), dimension(3,nat), intent(out) :: floc
 type(pawtab_type),intent(in)  :: pawtab(ntypes*usepaw)

!Local variables -------------------------
!scalars
 integer :: isx,isy,isz,iex,iey,iez,i1,i2,i3,j1,j2,j3,ind,iat,ityp,iloc,i3loc,msz,nloc
 logical :: perx,pery,perz,gox,goy,goz
 real(dp) :: arg,charge,cutoff,dvhgh,fxerf,fyerf,fzerf,fxgau,fygau,fzgau,forceleaked
 real(dp) :: forceloc,prefactor,rloc,rloc2,rhoel,rx,ry,rz,rzero,r2,tt,x,xp,y,z,Vel
 real(dp), dimension(4) :: cprime
!arrays
 real(dp) :: dvpawdr(1),rr(1)
 real(dp),pointer :: rad(:),vloc(:),d2vloc(:)

! *********************************************************************

   if (iproc == 0) write(std_out,'(1x,a)',advance='no')'Calculate local forces...'

!Conditions for periodicity in the three directions
   perx=(geocode /= 'F')
   pery=(geocode == 'P')
   perz=(geocode /= 'F')

   forceleaked=zero

   do iat=1,nat
     ityp=iatype(iat)
!  Coordinates of the center
     rx=rxyz(1,iat)
     ry=rxyz(2,iat)
     rz=rxyz(3,iat)

!  Initialization of the forces
!  ion-electron term, error function part
     fxerf=zero
     fyerf=zero
     fzerf=zero
!  ion-electron term, gaussian part
     fxgau=zero
     fygau=zero
     fzgau=zero

!  Building array of coefficients of the derivative of the gaussian part
     cprime(1)=2._dp*psppar(0,2,ityp)-psppar(0,1,ityp)
     cprime(2)=4._dp*psppar(0,3,ityp)-psppar(0,2,ityp)
     cprime(3)=6._dp*psppar(0,4,ityp)-psppar(0,3,ityp)
     cprime(4)=-psppar(0,4,ityp)

!  Determine number of local terms (HGH pot)
     nloc=0
     do iloc=1,4
       if (psppar(0,iloc,ityp).ne.zero) nloc=iloc
     end do

!  Some constants depending on the atom type
     rloc=psppar(0,0,ityp) ; rloc2=rloc**2
     charge=real(nelpsp(ityp),kind=dp)
     prefactor=charge/(2._dp*pi*sqrt(2._dp*pi)*rloc**5)

!  PAW specifics
     if (usepaw==1) then
       msz=pawtab(ityp)%wvl%rholoc%msz
       rad    => pawtab(ityp)%wvl%rholoc%rad(1:msz)
       vloc   => pawtab(ityp)%wvl%rholoc%d(1:msz,3)
       d2vloc => pawtab(ityp)%wvl%rholoc%d(1:msz,4)
       rzero=rad(1);if (rad(1)<=1.d-10) rzero=rad(2)
     end if

!  Maximum extension of the gaussian
     cutoff=10._dp*rloc
     isx=floor((rx-cutoff)/hxh)
     isy=floor((ry-cutoff)/hyh)
     isz=floor((rz-cutoff)/hzh)
     iex=ceiling((rx+cutoff)/hxh)
     iey=ceiling((ry+cutoff)/hyh)
     iez=ceiling((rz+cutoff)/hzh)

!  Calculate the forces near the atom due to the gaussian
!  and error function parts of the potential
     do i3=isz,iez
       z=real(i3,kind=dp)*hzh-rz
       call ind_positions_mklocl(perz,i3,n3,j3,goz)
       if(fftn3_distrib(j3)==iproc) then
         i3loc=ffti3_local(j3)
         do i2=isy,iey
           y=real(i2,kind=dp)*hyh-ry
           call ind_positions_mklocl(pery,i2,n2,j2,goy)
           do i1=isx,iex
             x=real(i1,kind=dp)*hxh-rx
             call ind_positions_mklocl(perx,i1,n1,j1,gox)

             if (goz.and.goy.and.gox) then
               ind=j1+(j2-1)*n1+(i3loc-1)*n1*n2
               x=(gridcart(1,ind)-rx)
               y=(gridcart(2,ind)-ry)
               z=(gridcart(3,ind)-rz)
               r2=x**2+y**2+z**2
               xp=exp(-0.5_dp*r2/rloc**2)

!            Short range part
               rhoel=rho(ind)
!            HGH: V_S^prime=gaussian
               if (usepaw==0) then
                 if (nloc/=0) then
                   arg=r2/rloc**2
                   tt=cprime(nloc)
                   do iloc=nloc-1,1,-1
                     tt=arg*tt+cprime(iloc)
                   end do
                   forceloc=xp*tt*rhoel
                 else
                   forceloc=zero
                 end if
!            PAW: V_PAW^prime-V_L^prime
               else
                 rr(1)=sqrt(r2)
                 if (rr(1)>=rzero) then
                   call paw_splint_der(msz,rad,vloc,d2vloc,1,rr,dvpawdr)
                   call calcdVloc_mklocl(dvhgh,rr(1),rloc,charge)
                   forceloc=rhoel*rloc2*(dvpawdr(1)-dvhgh)/rr(1)
                 else
                   forceloc=rhoel*rloc2*dvloc_zero_mklocl(charge,rloc,msz,rad,vloc,d2vloc)
                 end if
               end if

               fxgau=fxgau+forceloc*x
               fygau=fygau+forceloc*y
               fzgau=fzgau+forceloc*z

!            Long range part: error function
               Vel=pot(ind)
               fxerf=fxerf+xp*Vel*x
               fyerf=fyerf+xp*Vel*y
               fzerf=fzerf+xp*Vel*z

             else if (nloc>0) then
               r2=x**2+y**2+z**2
               arg=r2/rloc**2
               xp=exp(-0.5_dp*arg)
               tt=cprime(nloc)
               do iloc=nloc-1,1,-1
                 tt=arg*tt+cprime(iloc)
               end do
               forceleaked=forceleaked+xp*(1._dp+tt)
             end if
           end do
         end do
       end if
     end do

!  Final result of the forces
     floc(1,iat)=(hxh*hyh*hzh*prefactor)*fxerf+(hxh*hyh*hzh/rloc**2)*fxgau
     floc(2,iat)=(hxh*hyh*hzh*prefactor)*fyerf+(hxh*hyh*hzh/rloc**2)*fygau
     floc(3,iat)=(hxh*hyh*hzh*prefactor)*fzerf+(hxh*hyh*hzh/rloc**2)*fzgau

   end do

   forceleaked=forceleaked*prefactor*hxh*hyh*hzh
   if (iproc.eq.0) write(std_out,'(a,1pe12.5)') 'done. Leaked force: ',forceleaked

   CONTAINS
!!***

!----------------------------------------------------------------------

!!****f* mklocl_realspace/calcdVloc_mklocl
!! NAME
!!  calcdVloc_mklocl
!!
!! FUNCTION
!!  Compute 1st-derivative of long-range HGH local ionic potential (derf)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      mklocl_realspace
!!
!! CHILDREN
!!
!! SOURCE

subroutine calcdVloc_mklocl(yy,xx,rloc,Z)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in)  :: xx,rloc,Z
 real(dp),intent(out) :: yy

!Local variables-------------------------------
 !scalars
 real(dp):: arg,tt

! *************************************************************************

     arg=xx/(sqrt(2._dp)*rloc)
     call derf_ab(tt,arg)
     yy=(Z/(xx**2))* ( tt - 2._dp/sqrt(pi)*arg*exp(-arg**2) )

   end subroutine calcdVloc_mklocl
!!***

!----------------------------------------------------------------------

!!****f* mklocl_wavelets/dvloc_zero_mklocl
!! NAME
!!  dvloc_zero_mklocl
!!
!! FUNCTION
!!  Use a quadratic interpolation to get limit of (1/x).dVloc(x)/dx at x->0
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

function dvloc_zero_mklocl(charge,rloc,msz,rad,vloc,d2vloc)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: msz
 real(dp) :: dvloc_zero_mklocl
 real(dp),intent(in)  :: charge,rloc
!arrays
 real(dp) :: rad(msz),vloc(msz),d2vloc(msz)

!Local variables-------------------------------
!scalars
 real(dp) :: y1,y2,y3,zz=0._dp
!arrays
 real(dp) :: ll(3),xx(3),yy(3)

! *************************************************************************

!Select 3 points x1,x2,x3 near 0
     if (rad(1)>1.d-10) then
       xx(1:3)=rad(1:3)
     else
       xx(1:3)=rad(2:4)
     end if

!Find the corresponding values of y=(V^PAW(x)-V^HGH(x))/x
     call paw_splint_der(msz,rad,vloc,d2vloc,3,xx,yy)
     call calcdVloc_mklocl(y1,xx(1),rloc,charge)
     call calcdVloc_mklocl(y2,xx(2),rloc,charge)
     call calcdVloc_mklocl(y3,xx(3),rloc,charge)
     yy(1)=(yy(1)-y1)/xx(1)
     yy(2)=(yy(2)-y2)/xx(2)
     yy(3)=(yy(3)-y3)/xx(3)

!Find a polynomial of the form (z=0):
!P(z) = y1.L1(z) + y2.L2(z) + y3.L3(z)

!L1(z) = (z-x2)(z-x3)/((x1-x2)(x1-x3))
     ll(1)=(zz-xx(2))*(zz-xx(3))/((xx(1)-xx(2))*(xx(1)-xx(3)))
!L2(z) = (z-x1)(z-x3)/((x2-x1)(x2-x3))
     ll(2)=(zz-xx(1))*(zz-xx(3))/((xx(2)-xx(1))*(xx(2)-xx(3)))
!L3(z) = (z-x1)(z-x2)/((x3-x1)(x3-x2))
     ll(3)=(zz-xx(1))*(zz-xx(2))/((xx(3)-xx(1))*(xx(3)-xx(2)))

     dvloc_zero_mklocl=yy(1)*ll(1)+yy(2)*ll(2)+yy(3)*ll(3)

   end function dvloc_zero_mklocl
!!***

end subroutine local_forces_new
!!***

!----------------------------------------------------------------------

!!****f* mklocl_realspace/ind_positions_mklocl
!! NAME
!!  ind_positions_mklocl
!!
!! FUNCTION
!! determine the index in which the potential must be inserted, following the BC
!! determine also whether the index is inside or outside the box for free BC
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      mklocl_realspace
!!
!! CHILDREN
!!
!! SOURCE

subroutine ind_positions_mklocl(periodic,i,n,j,go)


!Arguments -------------------------------
 logical, intent(in) :: periodic
 integer, intent(in) :: i,n
 logical, intent(out) :: go
 integer, intent(out) :: j

! *********************************************************************

     if (periodic) then
       go=.true.
       j=modulo(i-1,n)+1
     else
       j=i
       go=(i >= 1 .and. i <= n)
     end if

   end subroutine ind_positions_mklocl
!!***

!!****f* ABINIT/mklocl_wavelets
!!
!! NAME
!! mklocl_wavelets
!!
!! FUNCTION
!! Compute the ionic local potential when the pseudo-potentials are GTH, using
!! the special decomposition of these pseudo. The resulting potential is computed with
!! free boundary conditions. It gives the same result than mklocl_realspace for the
!! GTH pseudo only with a different way to compute the potential.
!!
!! Optionally compute :
!!  option=1 : local ionic potential throughout unit cell
!!  option=2 : contribution of local ionic potential to E gradient wrt xred
!!
!! INPUTS
!!  efield (3)=external electric field
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms
!!  nfft=size of vpsp (local potential)
!!  nspden=number of spin-density components
!!  option=type of calculation (potential, forces, ...)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  xcart(3,natom)=cartesian atomic coordinates.
!!  wvl_den=density-potential BigDFT object
!!  wvl_descr=wavelet BigDFT object
!!
!! OUTPUT
!!  (if option==1) vpsp(nfft)=the potential resulting from the ionic
!!                 density of charge.
!!  (if option==2) grtn(3,natom)=grads of Etot wrt tn. These gradients are in
!!                 reduced coordinates. Multiply them by rprimd to get
!!                 gradients in cartesian coordinates.
!!
!! PARENTS
!!      mklocl,wvl_wfsinp_scratch
!!
!! CHILDREN
!!      calcdvloc_wvl,derf_ab,paw_splint_der
!!
!! SOURCE

subroutine mklocl_wavelets(efield, grtn, mpi_enreg, natom, nfft, &
     & nspden, option, rprimd, vpsp, wvl_den, wvl_descr, xcart)

 use defs_wvltypes
#if defined HAVE_BIGDFT
 use BigDFT_API, only : ELECTRONIC_DENSITY,createIonicPotential,local_forces
 use poisson_solver, only : H_potential
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option, natom, nfft, nspden
 type(MPI_type),intent(in) :: mpi_enreg
 type(wvl_denspot_type), intent(inout) :: wvl_den
 type(wvl_internal_type), intent(in) :: wvl_descr
!arrays
 real(dp),intent(in) :: rprimd(3,3),efield(3)
 real(dp),intent(inout) :: grtn(3,natom)
 real(dp), intent(inout) :: vpsp(nfft)
 real(dp),intent(inout) :: xcart(3,natom)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
!scalars
 integer :: i,i1,i2,i3,ia,ierr,igeo,me,nproc,shift,comm
 real(dp) :: energ
 character(len=500) :: message
!arrays
 real(dp) :: epot(3)
 real(dp) :: elecfield(3)=(/zero,zero,zero/) ! Not used here
 real(dp),allocatable :: gxyz(:,:),vhartr(:),rhov(:,:)
#endif

! *********************************************************************

#if defined HAVE_BIGDFT

 elecfield=zero !not used here

!Manage parallelism
 comm=mpi_enreg%comm_wvl
 nproc=xmpi_comm_size(comm)
 me=xmpi_comm_rank(comm)

!----------------------------------------------------------------------
! ----- Option 1: compute local ionic potential                   -----
!----------------------------------------------------------------------
 if (option == 1) then

!  We get the kernel for the Poisson solver (used to go from the ionic
!  charge to the potential).If the kernel is uncomputed, it does it now.
!  call psolver_kernel(wvl_den%denspot%dpbox%hgrids, 2, icoulomb, me, kernel, &
!&                     comm, wvl_den%denspot%dpbox%ndims, nproc, nscforder)
!  if (.not.associated(kernel%co%kernel)) then
!    call psolver_kernel(wvl_den%denspot%dpbox%hgrids, 1, icoulomb, me, kernel, &
!&                       comm, wvl_den%denspot%dpbox%ndims, nproc, nscforder)
!  end if

   message=ch10//' mklocl_wavelets: Create local potential from ions.'
   call wrtout(std_out,message,'COLL')

   shift = 1 + wvl_den%denspot%dpbox%ndims(1) * wvl_den%denspot%dpbox%ndims(2) &
&   * wvl_den%denspot%dpbox%i3xcsh

!  Call the BigDFT routine
   call createIonicPotential(wvl_descr%atoms%astruct%geocode, me, nproc, (me == 0), wvl_descr%atoms, &
&   xcart, wvl_den%denspot%dpbox%hgrids(1), wvl_den%denspot%dpbox%hgrids(2), &
&   wvl_den%denspot%dpbox%hgrids(3), &
&   elecfield, wvl_descr%Glr%d%n1, wvl_descr%Glr%d%n2, wvl_descr%Glr%d%n3, &
&   wvl_den%denspot%dpbox%n3pi, wvl_den%denspot%dpbox%i3s + wvl_den%denspot%dpbox%i3xcsh, &
&   wvl_den%denspot%dpbox%ndims(1), wvl_den%denspot%dpbox%ndims(2), &
&   wvl_den%denspot%dpbox%ndims(3), wvl_den%denspot%pkernel, vpsp(shift), 0.d0,wvl_descr%rholoc)

!  Copy vpsp into bigdft object:
   call wvl_rhov_abi2big(1,vpsp,wvl_den%denspot%v_ext,shift=(shift-1))

!  Eventually add the electric field
   if (maxval(efield) > tol12) then
     message=ch10//'mklocl_wavelets: Add the electric field.'
     call wrtout(std_out,message,'COLL')
!    We add here the electric field since in BigDFT, the field must be on x...
     epot(:) = real(0.5, dp) * efield(:) * wvl_den%denspot%dpbox%hgrids(:)
     do i3 = 1, wvl_den%denspot%dpbox%n3pi, 1
       ia = (i3 - 1) * wvl_den%denspot%dpbox%ndims(1) * wvl_den%denspot%dpbox%ndims(2)
       do i2 = -14, 2 * wvl_descr%Glr%d%n2 + 16, 1
         i = ia + (i2 + 14) * wvl_den%denspot%dpbox%ndims(1)
         do i1 = -14, 2 * wvl_descr%Glr%d%n1 + 16, 1
           i = i + 1
           vpsp(shift + i) = vpsp(shift + i) + &
&           epot(1) * real(i1 - wvl_descr%Glr%d%n1, dp) + &
&           epot(2) * real(i2 - wvl_descr%Glr%d%n2, dp) + &
&           epot(3) * real(i3 - wvl_descr%Glr%d%n3, dp)
         end do
       end do
     end do
   end if

!----------------------------------------------------------------------
! ----- Option 2: compute forces induced by local ionic potential -----
!----------------------------------------------------------------------
 else if (option == 2) then

   message=ch10//' mklocl_wavelets: compute local forces.'
   call wrtout(std_out,message,'COLL')

   if (wvl_den%denspot%rhov_is/=ELECTRONIC_DENSITY) then
     message='denspot bigdft datstructure should contain rhor!'
     MSG_BUG(message)
   end if

!  Extract density rhor from bigDFT datastructure
   ABI_ALLOCATE(rhov,(nfft, nspden))
   ABI_ALLOCATE(vhartr,(nfft))
   shift = wvl_den%denspot%dpbox%ndims(1) * wvl_den%denspot%dpbox%ndims(2) &
&   * wvl_den%denspot%dpbox%i3xcsh
   do i = 1, nfft
     rhov(i, 1) = wvl_den%denspot%rhov(i + shift)
     vhartr(i)  = wvl_den%denspot%rhov(i + shift)
   end do
   if (nspden == 2) then
     shift = shift + wvl_den%denspot%dpbox%ndims(1) * wvl_den%denspot%dpbox%ndims(2) &
&     * wvl_den%denspot%dpbox%n3d
     do i = 1, nfft
       rhov(i, 2) =             wvl_den%denspot%rhov(i + shift)
       vhartr(i)  = vhartr(i) + wvl_den%denspot%rhov(i + shift)
     end do
   end if

!  Compute Hartree potential from rhor
   call H_potential('D',wvl_den%denspot%pkernel,vhartr,vhartr,energ,zero,.false.)

!  Allocate temporary array for forces
   ABI_ALLOCATE(gxyz,(3, natom))

!  Calculate local part of the forces grtn (modified BigDFT routine)
   call local_forces_wvl(me,natom,xcart,&
&   wvl_den%denspot%dpbox%hgrids(1),&
&   wvl_den%denspot%dpbox%hgrids(2),&
&   wvl_den%denspot%dpbox%hgrids(3),&
&   wvl_descr%Glr%d%n1,wvl_descr%Glr%d%n2,wvl_descr%Glr%d%n3,&
&   wvl_den%denspot%dpbox%n3p,&
&   wvl_den%denspot%dpbox%i3s+wvl_den%denspot%dpbox%i3xcsh,&
&   wvl_den%denspot%dpbox%ndims(1),wvl_den%denspot%dpbox%ndims(2),&
&   rhov,vhartr,gxyz,wvl_descr)
!    call local_forces(me, wvl_descr%atoms, xcart, &
! &   wvl_den%denspot%dpbox%hgrids(1), wvl_den%denspot%dpbox%hgrids(2), &
! &   wvl_den%denspot%dpbox%hgrids(3), &
! &   wvl_descr%Glr%d%n1, wvl_descr%Glr%d%n2, wvl_descr%Glr%d%n3, &
! &   wvl_den%denspot%dpbox%n3p, wvl_den%denspot%dpbox%i3s + wvl_den%denspot%dpbox%i3xcsh, &
! &   wvl_den%denspot%dpbox%ndims(1), wvl_den%denspot%dpbox%ndims(2), &
! &   rhov, vhartr,gxyz,locstrten,charge)

!    Pending: floc,locstrten and charge are not used here.
!    Pending: put mpi_enreg%nscatterarr... in object denspot, initialize object, etc.

   if (nproc > 1) then
     call xmpi_sum(gxyz, comm, ierr)
   end if

!  Forces should be in reduced coordinates.
   do ia = 1, natom, 1
     do igeo = 1, 3, 1
       grtn(igeo, ia) = - rprimd(1, igeo) * gxyz(1, ia) &
&       - rprimd(2, igeo) * gxyz(2, ia) &
&       - rprimd(3, igeo) * gxyz(3, ia)
     end do
   end do

!  Deallocate local variables
   ABI_DEALLOCATE(vhartr)
   ABI_DEALLOCATE(rhov)
   ABI_DEALLOCATE(gxyz)

!----------------------------------------------------------------------

 else ! option switch
   message = 'Internal error, option should be 1 or 2!'
   MSG_ERROR(message)
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) option,natom,nfft,nspden,mpi_enreg%me,&
& wvl_den%symObj,wvl_descr%h(1),rprimd(1,1),efield(1),grtn(1,1),vpsp(1),xcart(1,1)
#endif

end subroutine mklocl_wavelets
!!***

!----------------------------------------------------------------------

!!****f* mklocl_wavelets/local_forces_wvl
!! NAME
!!  local_forces_wvl
!!
!! FUNCTION
!!
!! INPUTS
!!  hxh,hyh,hzh=wavelet grid spacings
!!  iproc=current MPI process number
!!  n1,n2,n3=number of wavelet points in each direction
!!  n1i,n2i=size of intermediate xy wvl grid
!!  n3pi=number of xy wvl planes handled by current MPI process
!!  i3s=starting index of local potential for current MPI process
!!  natom=number of atoms
!!  pot(*)=Hartree ionic potential
!!  rho(*)=electronic density
!!  rxyz(3,natom)=cartesian coordinates of atoms
!!  wvl=wavelet BigDFT object
!!
!! OUTPUT
!!  floc(3,natom)=local ionic potential contribution to forces
!!
!! PARENTS
!!      mklocl_wavelets
!!
!! CHILDREN
!!      calcdvloc_wvl,derf_ab,paw_splint_der
!!
!! SOURCE

subroutine local_forces_wvl(iproc,natom,rxyz,hxh,hyh,hzh,n1,n2,n3,n3pi,i3s,n1i,n2i,&
&                           rho,pot,floc,wvl)

 use defs_wvltypes
#if defined HAVE_BIGDFT
 use BigDFT_API, only : PSPCODE_PAW,ind_positions
#endif

!Arguments -------------------------------
!scalars
 integer,intent(in) :: i3s,iproc,n1,n1i,n2,n2i,n3,n3pi,natom
 real(dp),intent(in) :: hxh,hyh,hzh
 type(wvl_internal_type),intent(in) :: wvl
!arrays
 real(dp),intent(in) :: rxyz(3,natom)
 real(dp),dimension(*),intent(in) :: rho,pot
 real(dp),intent(out) :: floc(3,natom)

!Local variables -------------------------
#if defined HAVE_BIGDFT
!scalars
 integer :: i1,i2,i3,iat,iex,iey,iez,iloc,ind,isx,isy,isz,ityp
 integer :: j1,j2,j3,msz=0
 integer :: nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,nloc
 logical :: perx,pery,perz,gox,goy,goz,USE_PAW
 real(dp) :: arg,charge,cutoff,dvhgh,forceloc,fxerf,fyerf,fzerf,fxgau,fygau,fzgau
 real(dp) :: forceleaked,prefactor,r2,rhoel,rloc,rloc2,rx,ry,rz,rzero,tt,vel,x,xp,y,z
!arrays
 real(dp) :: cprime(4),dvpawdr(1),rr(1)
 real(dp),pointer :: psppar(:,:),rad(:),vloc(:),d2vloc(:)
#endif

! *********************************************************************

#if defined HAVE_BIGDFT

 if (iproc==0) write(std_out,'(a)',advance='no') ' Calculate local forces...'

!PAW or NCPP ?
 USE_PAW=any(wvl%atoms%npspcode==PSPCODE_PAW)

!Conditions for periodicity in the three directions
 perx=(wvl%atoms%astruct%geocode /= 'F')
 pery=(wvl%atoms%astruct%geocode == 'P')
 perz=(wvl%atoms%astruct%geocode /= 'F')
 call ext_buffers(perx,nbl1,nbr1)
 call ext_buffers(pery,nbl2,nbr2)
 call ext_buffers(perz,nbl3,nbr3)

 forceleaked=zero

 do iat=1,natom
   ityp=wvl%atoms%astruct%iatype(iat)

!  Coordinates of the center
   rx=rxyz(1,iat)
   ry=rxyz(2,iat)
   rz=rxyz(3,iat)

!  Initialization of the forces
!  ion-electron term, error function part
   fxerf=zero
   fyerf=zero
   fzerf=zero
!  ion-electron term, gaussian part
   fxgau=zero
   fygau=zero
   fzgau=zero

!  Building array of coefficients of the derivative of the gaussian part
   psppar => wvl%atoms%psppar(:,:,ityp)
   cprime(1)=2._dp*wvl%atoms%psppar(0,2,ityp)-wvl%atoms%psppar(0,1,ityp)
   cprime(2)=4._dp*wvl%atoms%psppar(0,3,ityp)-wvl%atoms%psppar(0,2,ityp)
   cprime(3)=6._dp*wvl%atoms%psppar(0,4,ityp)-wvl%atoms%psppar(0,3,ityp)
   cprime(4)=-wvl%atoms%psppar(0,4,ityp)

!  Determine number of local terms (HGH pot)
   nloc=0
   do iloc=1,4
     if (wvl%atoms%psppar(0,iloc,ityp).ne.zero) nloc=iloc
   end do

!  Some constants depending on the atom type
   rloc=wvl%atoms%psppar(0,0,ityp) ; rloc2=rloc**2
   charge=real(wvl%atoms%nelpsp(ityp),kind=dp)
   prefactor=charge/(2._dp*pi*sqrt(2._dp*pi)*rloc**5)

!  PAW specifics
   if (USE_PAW) then
     msz=wvl%rholoc%msz(ityp)
     rad    => wvl%rholoc%rad(1:msz,ityp)
     vloc   => wvl%rholoc%d(1:msz,3,ityp)
     d2vloc => wvl%rholoc%d(1:msz,4,ityp)
     rzero=rad(1);if (rzero<=1.d-10) rzero=rad(2)
   end if

!  Maximum extension of the gaussian
   cutoff=10._dp*rloc
   isx=floor((rx-cutoff)/hxh)
   isy=floor((ry-cutoff)/hyh)
   isz=floor((rz-cutoff)/hzh)
   iex=ceiling((rx+cutoff)/hxh)
   iey=ceiling((ry+cutoff)/hyh)
   iez=ceiling((rz+cutoff)/hzh)

!  Calculate the forces near the atom due to the gaussian
!  and error function parts of the potential
   if (n3pi>0) then
     do i3=isz,iez
       z=real(i3,kind=dp)*hzh-rz
       call ind_positions(perz,i3,n3,j3,goz)
       j3=j3+nbl3+1
       do i2=isy,iey
         y=real(i2,kind=dp)*hyh-ry
         call ind_positions(pery,i2,n2,j2,goy)
         do i1=isx,iex
           x=real(i1,kind=dp)*hxh-rx
           call ind_positions(perx,i1,n1,j1,gox)

           r2=x**2+y**2+z**2
           xp=exp(-0.5_dp*r2/rloc2)

           if ((j3>=i3s.and.j3<=i3s+n3pi-1) .and. (gox.and.goy)) then
             ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i

!            Short range part
             rhoel=rho(ind)
!            HGH: V_S^prime=gaussian
             if (.not.USE_PAW) then
               if (nloc/=0) then
                 arg=r2/rloc2
                 tt=cprime(nloc)
                 do iloc=nloc-1,1,-1
                   tt=arg*tt+cprime(iloc)
                 end do
                 forceloc=xp*tt*rhoel
               else
                 forceloc=zero
               end if
!            PAW: V_PAW^prime-V_L^prime
             else
               rr(1)=sqrt(r2)

               if (rr(1)>=rzero) then
                 call paw_splint_der(msz,rad,vloc,d2vloc,1,rr,dvpawdr)
                 call calcdVloc_wvl(dvhgh,rr(1),rloc,charge)
                 forceloc=rhoel*rloc2*(dvpawdr(1)-dvhgh)/rr(1)
               else
                 forceloc=rhoel*rloc2*dvloc_zero_wvl(charge,rloc,msz,rad,vloc,d2vloc)
               end if
             end if

             fxgau=fxgau+forceloc*x
             fygau=fygau+forceloc*y
             fzgau=fzgau+forceloc*z

!            Long range part: error function
             vel=pot(ind)
             fxerf=fxerf+xp*vel*x
             fyerf=fyerf+xp*vel*y
             fzerf=fzerf+xp*vel*z

           else if ((.not.goz).and.(nloc>0)) then
             arg=r2/rloc2
             tt=cprime(nloc)
             do iloc=nloc-1,1,-1
               tt=arg*tt+cprime(iloc)
             end do
             forceleaked=forceleaked+prefactor*xp*tt*rho(1)
           end if

         end do ! i1
       end do ! i2
     end do ! i3
   end if ! n3pi>0

!  Final result of the forces
   floc(1,iat)=(hxh*hyh*hzh*prefactor)*fxerf+(hxh*hyh*hzh/rloc2)*fxgau
   floc(2,iat)=(hxh*hyh*hzh*prefactor)*fyerf+(hxh*hyh*hzh/rloc2)*fygau
   floc(3,iat)=(hxh*hyh*hzh*prefactor)*fzerf+(hxh*hyh*hzh/rloc2)*fzgau

 end do ! iat

 forceleaked=forceleaked*hxh*hyh*hzh
 if (iproc.eq.0) write(std_out,'(a,1pe12.5)') 'done. Leaked force: ',forceleaked

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) i3s,iproc,n1,n1i,n2,n2i,n3,n3pi,natom,hxh,hyh,hzh,&
& rxyz(1,1),floc(1,1),rho(1),pot(1),wvl%h(1)
#endif

 CONTAINS
!!***

!----------------------------------------------------------------------

!!****f* mklocl_wavelets/calcdVloc_wvl
!! NAME
!!  calcdVloc_wvl
!!
!! FUNCTION
!!  Compute 1st-derivative of long-range HGH local ionic potential (derf)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      mklocl_wavelets
!!
!! CHILDREN
!!      calcdvloc_wvl,derf_ab,paw_splint_der
!!
!! SOURCE

subroutine calcdVloc_wvl(yy,xx,rloc,Z)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in)  :: xx,rloc,Z
 real(dp),intent(out) :: yy

!Local variables-------------------------------
 !scalars
 real(dp):: arg,tt

! *************************************************************************

   arg=xx/(sqrt(2._dp)*rloc)
   call derf_ab(tt,arg)
   yy=(Z/(xx**2))* ( tt - 2._dp/sqrt(pi)*arg*exp(-arg**2) )

 end subroutine calcdVloc_wvl
!!***

!----------------------------------------------------------------------

!!****f* mklocl_wavelets/dvloc_zero_wvl
!! NAME
!!  dvloc_zero_wvl
!!
!! FUNCTION
!!  Use a quadratic interpolation to get limit of (1/x).dVloc(x)/dx at x->0
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

function dvloc_zero_wvl(charge,rloc,msz,rad,vloc,d2vloc)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: msz
 real(dp) :: dvloc_zero_wvl
 real(dp),intent(in)  :: charge,rloc
!arrays
 real(dp) :: rad(msz),vloc(msz),d2vloc(msz)

!Local variables-------------------------------
!scalars
 real(dp) :: y1,y2,y3,zz=0._dp
!arrays
 real(dp) :: ll(3),xx(3),yy(3)

! *************************************************************************

!Select 3 points x1,x2,x3 near 0
   if (rad(1)>1.d-10) then
     xx(1:3)=rad(1:3)
   else
     xx(1:3)=rad(2:4)
   end if

!Find the corresponding values of y=(V^PAW(x)-V^HGH(x))/x
   call paw_splint_der(msz,rad,vloc,d2vloc,3,xx,yy)
   call calcdVloc_wvl(y1,xx(1),rloc,charge)
   call calcdVloc_wvl(y2,xx(2),rloc,charge)
   call calcdVloc_wvl(y3,xx(3),rloc,charge)
   yy(1)=(yy(1)-y1)/xx(1)
   yy(2)=(yy(2)-y2)/xx(2)
   yy(3)=(yy(3)-y3)/xx(3)

!Find a polynomial of the form (z=0):
!P(z) = y1.L1(z) + y2.L2(z) + y3.L3(z)

!L1(z) = (z-x2)(z-x3)/((x1-x2)(x1-x3))
   ll(1)=(zz-xx(2))*(zz-xx(3))/((xx(1)-xx(2))*(xx(1)-xx(3)))
!L2(z) = (z-x1)(z-x3)/((x2-x1)(x2-x3))
   ll(2)=(zz-xx(1))*(zz-xx(3))/((xx(2)-xx(1))*(xx(2)-xx(3)))
!L3(z) = (z-x1)(z-x2)/((x3-x1)(x3-x2))
   ll(3)=(zz-xx(1))*(zz-xx(2))/((xx(3)-xx(1))*(xx(3)-xx(2)))

   dvloc_zero_wvl=yy(1)*ll(1)+yy(2)*ll(2)+yy(3)*ll(3)

 end function dvloc_zero_wvl
!!***

end subroutine local_forces_wvl
!!***

end module m_mklocl_realspace
!!***
