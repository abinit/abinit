!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dens
!! NAME
!!  m_dens
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2008-2016 ABINIT group ()
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_dens

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_splines

 implicit none

 private

! public procedures.
 public :: dens_hirsh      ! Compute the Hirshfeld charges
!!***

contains

!----------------------------------------------------------------------

!!****f* m_dens/dens_hirsh
!! NAME
!! dens_hirsh
!!
!! FUNCTION
!! Compute the Hirshfeld charges
!!
!! INPUTS
!!  mpoint=Maximum number of points in radial meshes.
!!  radii(mpoint, ntypat)=Radial meshes for each type
!!  aeden(mpoint, nytpat)=All-electron densities.
!!  npoint(ntypat)=The number of the last point with significant density is stored in npoint(itypat)
!!  minimal_den=Tolerance on the minum value of the density
!!  grid_den(nrx,nry,nrz)= density on the grid
!!  natom = number of atoms in the unit cell
!!  nrx,nry,nrz= number of points in the grid for the three directions
!!  ntypat=number of types of atoms in unit cell.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  typat(natom)=type of each atom
!!  xcart(3,natom) = different positions of the atoms in the unit cell
!!  zion=(ntypat)gives the ionic charge for each type of atom
!!  znucl(ntypat)=gives the nuclear number for each type of atom
!!  prtcharge=1 to write the Hirshfeld charge decomposition
!!
!! OUTPUT
!!  hcharge(natom), hden(natom), hweight(natom)= Hirshfeld charges, densities, weights.
!!
!! PARENTS
!!      m_cut3d
!!
!! CHILDREN
!!      metric,spline,xcart2xred
!!
!! SOURCE

subroutine dens_hirsh(mpoint,radii,aeden,npoint,minimal_den,grid_den, &
  natom,nrx,nry,nrz,ntypat,rprimd,xcart,typat,zion,znucl,prtcharge,hcharge,hden,hweight)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dens_hirsh'
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nrx,nry,nrz,ntypat,prtcharge,mpoint
 real(dp),intent(in) :: minimal_den
!arrays
 integer,intent(in) :: typat(natom),npoint(ntypat)
 real(dp),intent(in) :: grid_den(nrx,nry,nrz),rprimd(3,3),zion(ntypat)
 real(dp),intent(in) :: znucl(ntypat)
 real(dp),intent(in) :: xcart(3,natom)
 real(dp),intent(in) :: radii(mpoint,ntypat),aeden(mpoint,ntypat)
 real(dp),intent(out) :: hcharge(natom),hden(natom),hweight(natom)

!Local variables -------------------------
!scalars
 integer :: i1,i2,i3,iatom,icell,igrid,ii,inmax,inmin,istep,itypat
 integer :: k1,k2,k3,mcells,nfftot,ngoodpoints,npt
 real(dp) :: aa,bb,coeff1,coeff2,coeff3,den,factor,h_inv,hh,maxrad
 real(dp) :: rr,rr2,total_charge,total_weight,total_zion,ucvol
 real(dp) :: yp1,ypn
 !character(len=500) :: msg
!arrays
 integer :: highest(3),lowest(3)
 integer,allocatable :: ncells(:)
 real(dp) :: coordat(3),coord23_1,coord23_2,coord23_3,diff1,diff2,diff3,gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp) :: vperp(3),width(3)
 real(dp),allocatable :: coord1(:,:),local_den(:,:,:,:)
 real(dp),allocatable :: step(:,:),sum_den(:,:,:),work(:)
 real(dp),allocatable :: xcartcells(:,:,:),xred(:,:),yder2(:)

! *********************************************************************

!1. Read the 1D all-electron atomic files
!Store the radii in radii(:,itypat), and the all-electron
!densities in aeden(:,itypat). The number of the last
!point with significant density is stored in npoint(itypat)

!2. Compute the list of atoms that are sufficiently close to the cell

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 nfftot=nrx*nry*nrz
!DEBUG
!do k3=1,nrz
!do k2=1,nry
!do k1=1,nrx
!total_charge=total_charge+grid_den(k1,k2,k3)
!end do
!end do
!end do
!write(std_out,*)' total_charge=',total_charge*ucvol/dble(nfftot)
!ENDDEBUG

 ABI_ALLOCATE(xred,(3,natom))
 call xcart2xred(natom,rprimd,xcart,xred)

!Compute the widths of the cell
!First width : perpendicular vector length
 vperp(:)=rprimd(:,1)-rprimd(:,2)*rmet(1,2)/rmet(2,2) -rprimd(:,3)*rmet(1,3)/rmet(3,3)
 width(1)=sqrt(dot_product(vperp,vperp))
!Second width
 vperp(:)=rprimd(:,2)-rprimd(:,1)*rmet(2,1)/rmet(1,1) -rprimd(:,3)*rmet(2,3)/rmet(3,3)
 width(2)=sqrt(dot_product(vperp,vperp))
!Third width
 vperp(:)=rprimd(:,3)-rprimd(:,1)*rmet(3,1)/rmet(1,1) -rprimd(:,2)*rmet(3,2)/rmet(2,2)
 width(3)=sqrt(dot_product(vperp,vperp))

!Compute the number of cells that will make up the supercell
 ABI_ALLOCATE(ncells,(natom))
 mcells=0
 do iatom=1,natom
   itypat=typat(iatom)
   maxrad=radii(npoint(itypat),itypat)
!  Compute the lower and higher indices of the supercell
!  for this atom
   do ii=1,3
     lowest(ii)=floor(-xred(ii,iatom)-maxrad/width(ii))
     highest(ii)=ceiling(-xred(ii,iatom)+maxrad/width(ii)+1)
!    Next coding, still incorrect
!    lowest(ii)=floor(xred(ii,iatom)-maxrad/width(ii))-1
!    highest(ii)=ceiling(xred(ii,iatom)+maxrad/width(ii))+1
!    Old coding incorrect
!    lowest(ii)=ceiling(-xred(ii,iatom)-maxrad/width(ii))
!    highest(ii)=floor(-xred(ii,iatom)+maxrad/width(ii)+1)
   end do
   ncells(iatom)=(highest(1)-lowest(1)+1)* &
&   (highest(2)-lowest(2)+1)* &
&   (highest(3)-lowest(3)+1)
!  DEBUG
!  write(std_out,*)' maxrad=',maxrad
!  write(std_out,*)' lowest(:)=',lowest(:)
!  write(std_out,*)' highest(:)=',highest(:)
!  write(std_out,*)' ncells(iatom)=',ncells(iatom)
!  ENDDEBUG
 end do
 mcells=maxval(ncells(:))

!Compute, for each atom, the set of image atoms in the whole supercell
 ABI_ALLOCATE(xcartcells,(3,mcells,natom))
 do iatom=1,natom
   itypat=typat(iatom)
   maxrad=radii(npoint(itypat),itypat)
!  Compute the lower and higher indices of the supercell
!  for this atom

   do ii=1,3
     lowest(ii)=floor(-xred(ii,iatom)-maxrad/width(ii))
     highest(ii)=ceiling(-xred(ii,iatom)+maxrad/width(ii)+1)
   end do
   icell=0
   do i1=lowest(1),highest(1)
     do i2=lowest(2),highest(2)
       do i3=lowest(3),highest(3)
         icell=icell+1
         xcartcells(:,icell,iatom)=xcart(:,iatom)+i1*rprimd(:,1)+i2*rprimd(:,2)+i3*rprimd(:,3)
       end do
     end do
   end do
 end do

!Compute, for each atom, the all-electron pro-atom density
!at each point in the primitive cell
 ABI_ALLOCATE(local_den,(nrx,nry,nrz,natom))
 ABI_ALLOCATE(step,(2,mpoint))
 ABI_ALLOCATE(work,(mpoint))
 ABI_ALLOCATE(yder2,(mpoint))
 ABI_ALLOCATE(coord1,(3,nrx))
 coeff1=one/nrx
 coeff2=one/nry
 coeff3=one/nrz

 do iatom=1,natom
   itypat=typat(iatom)
   npt=npoint(itypat)
   maxrad=radii(npt,itypat)
!   write(std_out,*)
!   write(std_out,'(a,i4)' )' hirsh : accumulating density for atom ',iatom
!  write(std_out,*)' ncells(iatom)=',ncells(iatom)
   do istep=1,npt-1
     step(1,istep)=radii(istep+1,itypat) - radii(istep,itypat)
     step(2,istep)=one/step(1,istep)
   end do
!  Approximate first derivative for small radii
   yp1=(aeden(2,itypat)-aeden(1,itypat))/(radii(2,itypat)-radii(1,itypat))
   ypn=zero
   call spline(radii(1:npt,itypat),aeden(1:npt,itypat),npt,yp1,ypn,yder2)

   local_den(:,:,:,iatom)=zero

!  Big loop on the cells
   do icell=1,ncells(iatom)
!    write(std_out,*)' icell=',icell
     coordat(:)=xcartcells(:,icell,iatom)

!    Big loop on the grid points
     do k1 = 1,nrx
       coord1(:,k1)=rprimd(:,1)*(k1-1)*coeff1
     end do
     do k3 = 1, nrz
       do k2 = 1, nry
         coord23_1=rprimd(1,2)*(k2-1)*coeff2+rprimd(1,3)*(k3-1)*coeff3-coordat(1)
         coord23_2=rprimd(2,2)*(k2-1)*coeff2+rprimd(2,3)*(k3-1)*coeff3-coordat(2)
         coord23_3=rprimd(3,2)*(k2-1)*coeff2+rprimd(3,3)*(k3-1)*coeff3-coordat(3)
         do k1 = 1, nrx
           diff1=coord1(1,k1)+coord23_1
           diff2=coord1(2,k1)+coord23_2
           diff3=coord1(3,k1)+coord23_3
           rr2=diff1**2+diff2**2+diff3**2
           if(rr2<maxrad**2)then

             rr=sqrt(rr2)
!            Find the index of the radius by bissection
             if (rr < radii(1,itypat)) then
!              Linear extrapolation
               den=aeden(1,itypat)+(rr-radii(1,itypat))/(radii(2,itypat)-radii(1,itypat))&
&               *(aeden(2,itypat)-aeden(1,itypat))
             else
!              Use the spline interpolation
!              Find the index of the radius by bissection
               inmin=1
               inmax=npt
               igrid=1
               do
                 if(inmax-inmin==1)exit
                 igrid=(inmin+inmax)/2
                 if(rr>=radii(igrid,itypat))then
                   inmin=igrid
                 else
                   inmax=igrid
                 end if
               end do
               igrid=inmin
!              write(std_out,*)' igrid',igrid

               hh=step(1,igrid)
               h_inv=step(2,igrid)
               aa= (radii(igrid+1,itypat)-rr)*h_inv
               bb= (rr-radii(igrid,itypat))*h_inv
               den = aa*aeden(igrid,itypat) + bb*aeden(igrid+1,itypat)  &
&               +( (aa*aa*aa-aa)*yder2(igrid)         &
&               +(bb*bb*bb-bb)*yder2(igrid+1) ) *hh*hh*sixth
             end if ! Select small radius or spline

             local_den(k1,k2,k3,iatom)=local_den(k1,k2,k3,iatom)+den
           end if ! dist2<maxrad

         end do ! k1
       end do ! k2
     end do ! k3

   end do ! icell
 end do ! iatom

!Compute, the total all-electron density at each point in the primitive cell
 ABI_ALLOCATE(sum_den,(nrx,nry,nrz))
 sum_den(:,:,:)=zero
 do iatom=1,natom
   sum_den(:,:,:)=sum_den(:,:,:)+local_den(:,:,:,iatom)
 end do

!DEBUG
!do k3=1,nrz
!do k2=1,nry
!do k1=1,nrx
!write(std_out,'(3i4,3es16.6)' )k1,k2,k3,local_den(k1,k2,k3,1:2),sum_den(k1,k2,k3)
!end do
!end do
!end do
!write(std_out,*)' hirsh : before accumulate the integral of the density'
!ENDDEBUG


!Accumulate the integral of the density, to get Hirshfeld charges
!There is a minus sign because the electron has a negative charge
 ngoodpoints = 0
 hcharge(:)=zero
 hweight(:)=zero
 do k3=1,nrz
   do k2=1,nry
     do k1=1,nrx
!      Use minimal_den in order to avoid divide by zero
       if (abs(sum_den(k1,k2,k3)) > minimal_den) then
         ngoodpoints = ngoodpoints+1
         factor=grid_den(k1,k2,k3)/(sum_den(k1,k2,k3)+minimal_den)
         do iatom=1,natom
           hden(iatom)=hden(iatom)+local_den(k1,k2,k3,iatom)
           hcharge(iatom)=hcharge(iatom)-local_den(k1,k2,k3,iatom)*factor
           hweight(iatom)=hweight(iatom)+local_den(k1,k2,k3,iatom)/(sum_den(k1,k2,k3)+minimal_den)
         end do
       end if
     end do
   end do
 end do

!DEBUG
!do iatom=1,natom
!write(std_out,'(i9,3es17.6)' )iatom,hden(iatom),hcharge(iatom),hweight(iatom)
!end do
!ENDDEBUG

 hcharge(:)=hcharge(:)*ucvol/dble(nfftot)
 hweight(:)=hweight(:)/dble(nfftot)

!Check on the total charge
 total_zion=sum(zion(typat(1:natom)))
 total_charge=sum(hcharge(1:natom))
 total_weight=sum(hweight(1:natom))

!DEBUG
!write(std_out,*)' ngoodpoints = ', ngoodpoints, ' out of ', nfftot
!write(std_out,*)' total_weight=',total_weight
!write(std_out,*)' total_weight=',total_weight
!ENDDEBUG

!Output
 if (prtcharge == 1) then
     write(std_out,*)
     write(std_out,*)'    Hirshfeld analysis'
     write(std_out,*)'    Atom       Zion       Electron  Charge       Net charge '
     write(std_out,*)
     do iatom=1,natom
       write(std_out,'(i9,3es17.6)' )&
&       iatom,zion(typat(iatom)),hcharge(iatom),hcharge(iatom)+zion(typat(iatom))
     end do
     write(std_out,*)
     write(std_out,'(a,3es17.6)')'    Total',total_zion,total_charge,total_charge+total_zion
     write(std_out,*)
 end if

 ABI_DEALLOCATE(coord1)
 ABI_DEALLOCATE(local_den)
 ABI_DEALLOCATE(ncells)
 ABI_DEALLOCATE(step)
 ABI_DEALLOCATE(sum_den)
 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(xcartcells)
 ABI_DEALLOCATE(xred)
 ABI_DEALLOCATE(yder2)

end subroutine dens_hirsh
!!***

end module m_dens
!!***
