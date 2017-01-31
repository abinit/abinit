!{\src2tex{textfont=tt}}
!!****f* ABINIT/mklocl_realspace
!! NAME
!!  mklocl_realspace
!!
!! FUNCTION
!! This method is equivalent to mklocl_recipspace except that
!! it uses real space pseudo-potentials. It is usefull for isolated
!! systems. Then the option 3 and 4 are not available for this
!! implementation.
!!
!! Optionally compute :
!!  option=1 : local ionic potential throughout unit cell
!!  option=2 : contribution of local ionic potential to E gradient wrt xred
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR,TRangel)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in unit cell.
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms.
!!  option= (see above)
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
!! SIDE EFFECTS
!!
!!
!! PARENTS
!!      mklocl
!!
!! CHILDREN
!!      ind_positions_
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mklocl_realspace(grtn,icoulomb,mpi_enreg, natom, nattyp, nfft, ngfft, &
                          & nscforder,nspden,ntypat,option,paral_kgb,pawtab,psps, rhog, rhor, &
                          & rprimd, typat,ucvol,usewvl,vpsp, xred)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_splines, only : splint
 use m_xmpi
 use m_profiling_abi
 use m_errors

 use m_mpinfo,      only : ptabs_fourdp
 use m_pawtab,      only : pawtab_type

#if defined HAVE_BIGDFT
 use BigDFT_API,    only : coulomb_operator,deallocate_coulomb_operator
 use defs_PSolver
#else
 use defs_wvltypes, only : coulomb_operator
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mklocl_realspace'
 use interfaces_18_timing
 use interfaces_41_geometry
 use interfaces_53_ffts
 use interfaces_62_poisson
 use interfaces_67_common, except_this_one => mklocl_realspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,nspden,ntypat,option
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
 type(pawtab_type),intent(in)  :: pawtab(ntypat*psps%usepaw)
!arrays
 integer,intent(in)  :: icoulomb,nscforder,paral_kgb,usewvl
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

 if (icoulomb == 1) then
   geocode='F'
 else if (icoulomb == 2) then
   geocode='S'
 end if

!Keep track of total time spent in mklocl
 if(option==2)then
   call timab(72,1,tsec)
 end if

 if (testing) then
   call system_clock(count_rate = countParSeconde)
   call system_clock(tpsStart, count_rate = countParSeconde)
 end if

 n1 = ngfft(1) ; n2 = ngfft(2) ; n3 = ngfft(3)
 nproc_fft = ngfft(10) ;  me_fft = ngfft(11)
 n3d = ngfft(13) !for parallel runs
 if(nproc_fft==1) n3d=n3  !for serial runs
 comm_fft=mpi_enreg%comm_fft

 if(me_fft /= mpi_enreg%me_fft .or. nproc_fft /= mpi_enreg%nproc_fft) then
   MSG_BUG("mpi_enreg%x_fft not equal to the corresponding values in ngfft")
 end if

!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)


!Store xcart for each atom
 ABI_ALLOCATE(xcart,(3, natom))
 call xred2xcart(natom, rprimd, xcart, xred)
!Store cartesian coordinates for each grid points
 ABI_ALLOCATE(gridcart,(3, nfft))
 call mkgrid_fft(ffti3_local,fftn3_distrib,gridcart,nfft,ngfft,rprimd)


!the branch with the HGH treatment of the PSP will presumably start here
!here we need to put the if statement for the PSP code =2,3,10 for GTH-HGH

!see whether all the PSP considered are of type GTH-HGH or PAW
 doIt=.true.
!doIt=.false.
 do ii=1,psps%npsp
   doIt=doIt .and.&
   (psps%pspcod(ii)==2 .or.psps%pspcod(ii)==3 .or. psps%pspcod(ii)==10 .or. psps%pspcod(ii)==7)
 end do

 if (doIt) then

!  definition of the grid spacings as in the kernel routine
   hgx = rprimd(1,1)/(n1)
   hgy = rprimd(2,2)/(n2)
   hgz = rprimd(3,3)/(n3)

   call psolver_kernel( (/ hgx, hgy, hgz /), 2, icoulomb, me_fft, kernel, comm_fft, (/n1,n2,n3/), &
&   nproc_fft, nscforder)

   if (option==1) then

     call createIonicPotential_new(fftn3_distrib,ffti3_local,&
&     geocode,me_fft, nproc_fft, natom, &
&     ntypat, typat, psps%gth_params%psppar, &
&     int(psps%ziontypat), xcart,gridcart, hgx,hgy,hgz, &
&     n1,n2,n3d,n3, kernel, vpsp, comm_fft,pawtab,psps%usepaw)

   else if (option ==2) then

!    the local forces with this formalism are calculated differently

!    Compute Hartree's potential from rhor.
     ABI_ALLOCATE(vhartr,(nfft))
     call psolver_hartree(entmp, (/ hgx, hgy, hgz /), icoulomb, me_fft, comm_fft, nfft, &
&     (/n1,n2,n3/), nproc_fft, nscforder, nspden, rhor, vhartr, &
&     usewvl)

     ABI_ALLOCATE(gxyz,(3, natom))
!    calculate local part of the forces grtn (inspired from BigDFT routine)
     call local_forces_new(fftn3_distrib,ffti3_local,geocode,me_fft, ntypat, natom, &
&     typat, xcart, gridcart, psps%gth_params%psppar, &
&     int(psps%ziontypat), hgx,hgy,hgz, n1,n2,n3,n3d,&
&     rhor,vhartr, gxyz)
     ABI_DEALLOCATE(vhartr)

!    Forces should be in reduced coordinates.
     do ia = 1, natom, 1
       do igeo = 1, 3, 1
         grtn(igeo, ia) = - rprimd(1, igeo) * gxyz(1, ia) - &
&         rprimd(2, igeo) * gxyz(2, ia) - &
&         rprimd(3, igeo) * gxyz(3, ia)
       end do
     end do

!    Deallocate local variables
     ABI_DEALLOCATE(gxyz)
   end if


   ABI_DEALLOCATE(xcart)
   ABI_DEALLOCATE(gridcart)

!  else statement for the non GTH-HGH PSP
 else

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
       call fourdp(1, rhog_testing, rhor_testing, -1, mpi_enreg, nfft, ngfft, paral_kgb,0)
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
     call fourdp(1, rhog_interpol, rhor_work, 1, mpi_enreg, nfft * n_interpol, ngfft_interpol, paral_kgb,0)

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

   ABI_DEALLOCATE(xcart)
   ABI_DEALLOCATE(gridcart)

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
 end if

 if(option==2)then
   call timab(72,2,tsec)
 end if

end subroutine mklocl_realspace
!!***

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
!!      ind_positions_
!!
!! SOURCE
subroutine createIonicPotential_new(fftn3_distrib,ffti3_local,geocode,iproc,&
&  nproc,nat,ntypes,iatype,psppar,nelpsp,rxyz,gridcart,&
&  hxh,hyh,hzh,n1i,n2i,n3d,n3i,kernel,pot_ion,spaceworld,pawtab,usepaw)

 use defs_datatypes
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use defs_basis,    only : std_out,std_out_default
 use defs_wvltypes, only : coulomb_operator
 use m_pawtab,      only : pawtab_type
 use m_splines,     only : splint

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'createIonicPotential_new'
 use interfaces_43_wvl_wrappers
 use interfaces_67_common, except_this_one => createIonicPotential_new
!End of the abilint section

implicit none

!Arguments -------------------------------

character(len=1), intent(in) :: geocode
integer, intent(in) :: iproc,nproc,ntypes,nat,n1i,n2i,n3i,n3d,spaceworld
integer, intent(in) :: usepaw
real(kind=8), intent(in) :: hxh,hyh,hzh
integer, dimension(nat), intent(in) :: iatype
integer, dimension(ntypes), intent(in) :: nelpsp
integer, dimension(*), intent(in) ::fftn3_distrib,ffti3_local
real(kind=8), dimension(3,n1i*n2i*n3d), intent(in) :: gridcart
real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
real(kind=8), dimension(3,nat), intent(in) :: rxyz
real(kind=8), dimension(*), intent(inout) :: pot_ion
type(pawtab_type),intent(in)  :: pawtab(ntypes*usepaw)
type(coulomb_operator), intent(in) :: kernel

!Local variables -------------------------
#if defined HAVE_BIGDFT
logical :: perx,pery,perz,gox,goy,goz
integer :: iat,i1,i2,i3,j1,j2,j3,isx,isy,isz,iex,iey,iez,ierr,ityp
integer :: ind,nloc,iloc,i3loc
real(kind=8) :: rholeaked,rloc,charge,cutoff,x,y,z,r2,arg,xp,tt
real(kind=8) :: raux,raux2,rx,ry,rz,rr
real(kind=8) :: tt_tot,rholeaked_tot
real(kind=8) :: ehart,eexcu,vexcu!,r2paw
real(kind=8) :: raux1(1),rr1(1)
#endif

! *********************************************************************

#if defined HAVE_BIGDFT

 if(nproc<0)then
   MSG_ERROR('nproc should not be negative')
 end if

!Ionic charge (must be calculated for the PS active processes)
 rholeaked=0.d0
!Ionic energy (can be calculated for all the processors)

!here we should insert the calculation of the ewald energy for the periodic BC case
!!!  eion=0.d0
!!!  do iat=1,nat
!!!     ityp=iatype(iat)
!!!     rx=rxyz(1,iat)
!!!     ry=rxyz(2,iat)
!!!     rz=rxyz(3,iat)
!!!     !    ion-ion interaction
!!!     do jat=1,iat-1
!!!        dist=sqrt( (rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2 )
!!!        jtyp=iatype(jat)
!!!        eion=eion+real(nelpsp(jtyp)*nelpsp(ityp),kind=8)/dist
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
   charge=real(nelpsp(ityp),kind=8)/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
   cutoff=10.d0*rloc

   isx=floor((rx-cutoff)/hxh)
   isy=floor((ry-cutoff)/hyh)
   isz=floor((rz-cutoff)/hzh)

   iex=ceiling((rx+cutoff)/hxh)
   iey=ceiling((ry+cutoff)/hyh)
   iez=ceiling((rz+cutoff)/hzh)

!  Calculate Ionic Density
!  using HGH parameters.
!  Eq. 1.104, T. Deutsch and L. Genovese, JDN. 12, 2011
   if(usepaw==0) then

     do i3=isz,iez
       z=real(i3,kind=8)*hzh-rz
       call ind_positions_(perz,i3,n3i,j3,goz)
       if(fftn3_distrib(j3)==iproc) then
         i3loc=ffti3_local(j3) 
         do i2=isy,iey
           y=real(i2,kind=8)*hyh-ry
           call ind_positions_(pery,i2,n2i,j2,goy)
           do i1=isx,iex
             x=real(i1,kind=8)*hxh-rx
             call ind_positions_(perx,i1,n1i,j1,gox)
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
!    Calculate Ionic Density using splines, 
!    PAW case
   else
!    r2paw=(pawtab(ityp)%rpaw)**2
     do i3=isz,iez
       z=real(i3,kind=8)*hzh-rz
       call ind_positions_(perz,i3,n3i,j3,goz)
       if(fftn3_distrib(j3)==iproc) then
         i3loc=ffti3_local(j3) 
         do i2=isy,iey
           y=real(i2,kind=8)*hyh-ry
           call ind_positions_(pery,i2,n2i,j2,goy)
           do i1=isx,iex
             x=real(i1,kind=8)*hxh-rx
             call ind_positions_(perx,i1,n1i,j1,gox)
             r2=x**2+y**2+z**2
             if (goz  .and. goy  .and. gox ) then
               ind=j1+(j2-1)*n1i+(i3loc-1)*n1i*n2i
               r2=(gridcart(1,ind)-rx)**2+(gridcart(2,ind)-ry)**2+(gridcart(3,ind)-rz)**2
             end if
!            if(r2>r2paw) cycle
             rr=sqrt(r2)
!            This converges very slow with hh
!            call splint(pawtab(ityp)%wvl%rholoc%msz,pawtab(ityp)%wvl%rholoc%rad(:),&
!            &             pawtab(ityp)%wvl%rholoc%d(:,1),pawtab(ityp)%wvl%rholoc%d(:,2),1,rr,raux,ierr)
!            Take the HGH form for rho_L (long range)
             arg=r2/rloc**2
             xp=exp(-.5d0*arg)
             raux=-xp*charge
             if (goz  .and. goy  .and. gox ) then
               pot_ion(ind)=pot_ion(ind)+raux
             else
               rholeaked=rholeaked-raux
             end if
           end do
         end do
       end if
     end do
   end if !if usepaw
 end do
!end if
!Check
 tt=0.d0
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

!here the value of the datacode must be kept fixed
!there can be some problems when running this stuff in parallel, if the ionic potential distribution does not agree with the
!plane distribution which is supposed to hold for the Poisson Solver
 call psolver(geocode,'D',iproc,nproc,n1i,n2i,n3i,0,hxh,hyh,hzh,&
& pot_ion,kernel%kernel,pot_ion,ehart,eexcu,vexcu,0.d0,.false.,1)



!The following is just for HGH/GTH pseudos
!write(std_out,*) 'ehartree',ehart
!if (n3i > 0) then
 do iat=1,nat
   ityp=iatype(iat)

   rx=rxyz(1,iat)
   ry=rxyz(2,iat)
   rz=rxyz(3,iat)

!  if (iat==1) then
!  write(std_out,*) rx/hxh,ry/hyh,rz/hzh
!  write(std_out,*) rx,ry,rz
!  write(std_out,*) hxh,hyh,hzh
!  end if

!  stop
!  determine number of local terms
   rloc=psppar(0,0,ityp)
   cutoff=10.d0*rloc

   isx=floor((rx-cutoff)/hxh)
   isy=floor((ry-cutoff)/hyh)
   isz=floor((rz-cutoff)/hzh)

   iex=ceiling((rx+cutoff)/hxh)
   iey=ceiling((ry+cutoff)/hyh)
   iez=ceiling((rz+cutoff)/hzh)

   if (usepaw==0) then
!    Only for HGH pseudos
!    Add the remaining local terms of Eq. (9)
!    in JCP 129, 014109(2008)
     nloc=0
     do iloc=1,4
       if (psppar(0,iloc,ityp).ne.0.d0) nloc=iloc
     end do



     if (nloc /= 0) then

!      write(std_out,*) 'nloc=',nloc

       do i3=isz,iez
         z=real(i3,kind=8)*hzh-rz
         call ind_positions_(perz,i3,n3i,j3,goz)
         if(fftn3_distrib(j3) == iproc .and. goz) then !MPI
           i3loc=ffti3_local(j3)
           if (goz) then
             do i2=isy,iey
               y=real(i2,kind=8)*hyh-ry
               call ind_positions_(pery,i2,n2i,j2,goy)
               if (goy) then
                 do i1=isx,iex
                   x=real(i1,kind=8)*hxh-rx
                   call ind_positions_(perx,i1,n1i,j1,gox)
                   if (gox) then
                     ind=j1+(j2-1)*n1i+(i3loc-1)*n1i*n2i
                     r2=(gridcart(1,ind)-rx)**2+(gridcart(2,ind)-ry)**2+(gridcart(3,ind)-rz)**2
!                    r2=x**2+y**2+z**2
                     arg=r2/rloc**2
                     xp=exp(-.5d0*arg)
                     tt=psppar(0,nloc,ityp)
                     do iloc=nloc-1,1,-1
                       tt=arg*tt+psppar(0,iloc,ityp)
                     end do
                     pot_ion(ind)=pot_ion(ind)+xp*tt
                   end if
                 end do
               end if
             end do
           end if
         end if
       end do
     end if !nloc
   else !HGH or PAW
!    For PAW, add V^PAW-V_L^HGH
     charge=real(nelpsp(ityp),kind=8)
     do i3=isz,iez
       z=real(i3,kind=8)*hzh-rz
       call ind_positions_(perz,i3,n3i,j3,goz)
       if(fftn3_distrib(j3) == iproc .and. goz) then !MPI
         i3loc=ffti3_local(j3)
         do i2=isy,iey
           y=real(i2,kind=8)*hyh-ry
           call ind_positions_(pery,i2,n2i,j2,goy)
           if (goy) then
             do i1=isx,iex
               x=real(i1,kind=8)*hxh-rx
               call ind_positions_(perx,i1,n1i,j1,gox)
               if (gox) then
                 ind=j1+(j2-1)*n1i+(i3loc-1)*n1i*n2i
                 r2=(gridcart(1,ind)-rx)**2+(gridcart(2,ind)-ry)**2+(gridcart(3,ind)-rz)**2
!                r2=x**2+y**2+z**2
                 rr=sqrt(r2)
!                1) V_L^HGH
                 if(rr>0.01d0) then
                   arg=rr/(sqrt(2.0)*rloc)
                   call derf_ab(tt,arg)
                   raux2=-charge/rr*tt
                 else
!                  In this case we deduce the values
!                  from a quadratic interpolation (due to 1/rr factor)
                   call interpol_vloc(rr,rloc,charge,raux2)
                 end if
!                2) V^PAW from splines
                 rr1(1)=rr
                 call splint(pawtab(ityp)%wvl%rholoc%msz,pawtab(ityp)%wvl%rholoc%rad(:),&
&                 pawtab(ityp)%wvl%rholoc%d(:,3),pawtab(ityp)%wvl%rholoc%d(:,4),&
&                 1,rr1,raux1,ierr)
                 raux=raux1(1)
                 pot_ion(ind)=pot_ion(ind)+raux-raux2
               end if
             end do
           end if
         end do
       end if
     end do
   end if !usepaw
 end do  !iat

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) geocode,iproc,nproc,ntypes,nat,n1i,n2i,n3i,n3d,spaceworld,usepaw,&
& hxh,hyh,hzh,iatype(1),nelpsp(1),fftn3_distrib(1),ffti3_local(1),gridcart(1,1),psppar(1,1,1),&
& rxyz(1,1),pot_ion(1),pawtab(1)%mesh_size,kernel%co
#endif

 contains
!!***

!!****f* mklocl_realspace/interpol_vloc
!! NAME
!!  interpol_vloc
!!
!! FUNCTION
!! We use a quadratic interpolation to get vloc(x)
!! useful for small values of x
!!
!! COPYRIGHT
!! Copyright (C) 2013-2016 ABINIT group (TRangel)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      mklocl_realspace
!!
!! CHILDREN
!!      ind_positions_
!!
!! SOURCE
subroutine interpol_vloc(xx,rloc,charge,yy)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'interpol_vloc'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  real(dp),intent(in)  :: xx,rloc,charge
  real(dp),intent(out) :: yy

!Local variables-------------------------------
 !scalars
  real(dp)::l0,l1,l2,x0,x1,x2,y0,y1,y2
 
! *************************************************************************

!  Find 3 points (x0,y0), (x1,y1), (x2,y2).
   x0=0.01d0; x1=0.02d0; x2=0.03d0
   call calcVloc(y0,x0,rloc,charge)
   call calcVloc(y1,x1,rloc,charge)
   call calcVloc(y2,x2,rloc,charge)

!  Find a polynomial of the form:
!  P(x)=y0L0(x) + y1L1(x) + y2L2(x)

!  L0(x) = (x-x1)(x-x2)/((x0-x1)(x0-x2))
   l0=(xx-x1)*(xx-x2)/((x0-x1)*(x0-x2))
!  L1(x) = (x-x0)(x-x2)/((x1-x0)(x1-x2))
   l1=(xx-x0)*(xx-x2)/((x1-x0)*(x1-x2))
!  L2(x) = (x-x0)(x-x1)/((x2-x0)(x2-x1))
   l2=(xx-x0)*(xx-x1)/((x2-x0)*(x2-x1))

   yy=y0*l0+y1*l1+y2*l2

 end subroutine interpol_vloc
!!***

!!****f* mklocl_realspace/calcVloc
!! NAME
!!  calcVloc
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2013-2016 ABINIT group (TRangel)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      mklocl_realspace
!!
!! CHILDREN
!!      ind_positions_
!!
!! SOURCE
subroutine calcVloc(yy,xx,rloc,Z)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calcVloc'
 use interfaces_43_wvl_wrappers
!End of the abilint section

 implicit none
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

 end subroutine calcVloc
!!***

end subroutine createIonicPotential_new
!!***

!!****f* mklocl_realspace/ind_positions_
!! NAME
!!  ind_positions_
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
!!      mkcore_paw,mklocl_realspace
!!
!! CHILDREN
!!      ind_positions_
!!
!! SOURCE

subroutine ind_positions_(periodic,i,n,j,go)

 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ind_positions_'
!End of the abilint section

  implicit none

!Arguments -------------------------------

  logical, intent(in) :: periodic
  integer, intent(in) :: i,n
  logical, intent(out) :: go
  integer, intent(out) :: j

!Local variables -------------------------

! *********************************************************************

   if (periodic) then
     go=.true.
     j=modulo(i-1,n)+1
   else
     j=i
     if (i >= 1 .and. i <= n) then
       go=.true.
     else
       go=.false.
     end if
   end if

 end subroutine ind_positions_
!!***

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
!!      ind_positions_
!!
!! SOURCE
subroutine local_forces_new(fftn3_distrib,ffti3_local,&
     geocode,iproc,ntypes,nat,iatype,rxyz,gridcart,psppar,nelpsp,hxh,hyh,hzh,&
     n1,n2,n3,n3d,rho,pot,floc)

 use m_profiling_abi
! Calculates the local forces acting on the atoms belonging to iproc
  use defs_basis,only: std_out,std_out_default

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'local_forces_new'
 use interfaces_67_common, except_this_one => local_forces_new
!End of the abilint section

  implicit none

!Arguments -------------------------------

  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,ntypes,nat,n1,n2,n3,n3d
  real(kind=8), intent(in) :: hxh,hyh,hzh
  real(kind=8), dimension(3,n1*n2*n3d), intent(in) :: gridcart
  real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(*), intent(in) :: rho,pot
  integer, dimension(*), intent(in) ::fftn3_distrib,ffti3_local
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: nelpsp
  real(kind=8), dimension(3,nat), intent(out) :: floc

!Local variables -------------------------

  logical :: perx,pery,perz,gox,goy,goz
  real(kind=8) :: pi,prefactor,cutoff,rloc,Vel,rhoel
  real(kind=8) :: fxerf,fyerf,fzerf,fxgau,fygau,fzgau,forceleaked,forceloc
  real(kind=8) :: rx,ry,rz,x,y,z,arg,r2,xp,tt
  integer :: isx,isy,isz,iex,iey,iez,i1,i2,i3,j1,j2,j3,ind,iat,ityp,nloc,iloc,i3loc
  !array of coefficients of the derivative
  real(kind=8), dimension(4) :: cprime

! *********************************************************************

   pi=4.d0*atan(1.d0)

   if (iproc == 0) write(std_out,'(1x,a)',advance='no')'Calculate local forces...'
   forceleaked=0.d0

!conditions for periodicity in the three directions
   perx=(geocode /= 'F')
   pery=(geocode == 'P')
   perz=(geocode /= 'F')

   do iat=1,nat
     ityp=iatype(iat)
!  coordinates of the center
     rx=rxyz(1,iat)
     ry=rxyz(2,iat)
     rz=rxyz(3,iat)

!  inizialization of the forces
!  ion-electron term, error function part
     fxerf=0.d0
     fyerf=0.d0
     fzerf=0.d0
!  ion-electron term, gaussian part
     fxgau=0.d0
     fygau=0.d0
     fzgau=0.d0

!  building array of coefficients of the derivative of the gaussian part
     cprime(1)=2.d0*psppar(0,2,ityp)-psppar(0,1,ityp)
     cprime(2)=4.d0*psppar(0,3,ityp)-psppar(0,2,ityp)
     cprime(3)=6.d0*psppar(0,4,ityp)-psppar(0,3,ityp)
     cprime(4)=-psppar(0,4,ityp)

!  determine number of local terms
     nloc=0
     do iloc=1,4
       if (psppar(0,iloc,ityp).ne.0.d0) nloc=iloc
     end do

!  local part
     rloc=psppar(0,0,ityp)
     prefactor=real(nelpsp(ityp),kind=8)/(2.d0*pi*sqrt(2.d0*pi)*rloc**5)
!  maximum extension of the gaussian
     cutoff=10.d0*rloc
     isx=floor((rx-cutoff)/hxh)
     isy=floor((ry-cutoff)/hyh)
     isz=floor((rz-cutoff)/hzh)

     iex=ceiling((rx+cutoff)/hxh)
     iey=ceiling((ry+cutoff)/hyh)
     iez=ceiling((rz+cutoff)/hzh)

!  calculate the forces near the atom due to the error function part of the potential
     do i3=isz,iez
       z=real(i3,kind=8)*hzh-rz
       call ind_positions_(perz,i3,n3,j3,goz)
       if(fftn3_distrib(j3)==iproc) then
         i3loc=ffti3_local(j3)
         do i2=isy,iey
           y=real(i2,kind=8)*hyh-ry
           call ind_positions_(pery,i2,n2,j2,goy)
           do i1=isx,iex
             x=real(i1,kind=8)*hxh-rx
             call ind_positions_(perx,i1,n1,j1,gox)
             r2=x**2+y**2+z**2
             if (goz  .and. goy  .and. gox ) then
               ind=j1+(j2-1)*n1+(i3loc-1)*n1*n2
               x=(gridcart(1,ind)-rx)
               y=(gridcart(2,ind)-ry)
               z=(gridcart(3,ind)-rz)
               r2=x**2+y**2+z**2
             end if
             arg=r2/rloc**2
             xp=exp(-.5d0*arg)
             if (goz  .and. goy  .and. gox ) then
!            gaussian part
               if (nloc /= 0) then
                 tt=cprime(nloc)
                 do iloc=nloc-1,1,-1
                   tt=arg*tt+cprime(iloc)
                 end do
                 rhoel=rho(ind)
                 forceloc=xp*tt*rhoel
                 fxgau=fxgau+forceloc*x
                 fygau=fygau+forceloc*y
                 fzgau=fzgau+forceloc*z
               end if
!            error function part
               Vel=pot(ind)
               fxerf=fxerf+xp*Vel*x
               fyerf=fyerf+xp*Vel*y
               fzerf=fzerf+xp*Vel*z
             else
               forceleaked=forceleaked+xp*(1.d0+tt)
             end if
           end do
         end do
       end if
     end do

!  final result of the forces

     floc(1,iat)=(hxh*hyh*hzh*prefactor)*fxerf+(hxh*hyh*hzh/rloc**2)*fxgau
     floc(2,iat)=(hxh*hyh*hzh*prefactor)*fyerf+(hxh*hyh*hzh/rloc**2)*fygau
     floc(3,iat)=(hxh*hyh*hzh*prefactor)*fzerf+(hxh*hyh*hzh/rloc**2)*fzgau

   end do

   forceleaked=forceleaked*prefactor*hxh*hyh*hzh
   if (iproc.eq.0) write(std_out,'(a,1pe12.5)') 'done. Leaked force: ',forceleaked

 end subroutine local_forces_new
!!***

