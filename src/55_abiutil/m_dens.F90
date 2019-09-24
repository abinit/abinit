!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dens
!! NAME
!!  m_dens
!!
!! FUNCTION
!! Module containing the definition of the constrained_dft_t data type and methods used to handle it,
!! and also includes the computation of integrated atomic charge and magnetization, as well as Hirshfeld charges.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (MT,ILuk,MVer,EB,SPr)
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
 use m_abicore
 use m_xmpi
 use m_splines

 use defs_abitypes,   only : MPI_type
 use m_fft,           only : fourdp
 use m_time,          only : timab
 use m_numeric_tools, only : wrap2_zero_one
 use m_io_tools,      only : open_file
 use m_geometry,      only : dist2, xcart2xred, metric
 use m_mpinfo,        only : ptabs_fourdp

 implicit none

 private

 public :: dens_hirsh              ! Compute the Hirshfeld charges
 public :: add_atomic_fcts         ! Add atomic functions to real space function
 public :: constrained_dft_ini     ! Initialize the constrained_dft datastructure
 public :: constrained_dft_free    ! Free the constrained_dft datastructure
 public :: constrained_residual    ! Recompute the potential residual, to account for constraints
 public :: mag_penalty             ! Compute the potential corresponding to constrained magnetic moments (using add_atomic_fcts) with the penalty function.
 public :: mag_penalty_e           ! Compute the energy corresponding to constrained magnetic moments.
 public :: calcdensph              ! Compute and print integral of total density inside spheres around atoms.
!!***

!----------------------------------------------------------------------

!!****t* m_dens/constrained_dft_t
!! NAME
!! constrained_dft_t
!!
!! FUNCTION
!! Structure gathering the relevant information for constrained DFT calculations
!!
!! SOURCE

 type,public :: constrained_dft_t

!scalars
  integer :: natom                           ! Number of atoms
  integer :: nfftf                           ! Number of FFT grid points (for this processor) for the "fine" grid
  integer :: nspden                          ! Number of spin-density components
  integer :: ntypat                          ! Number of type of atoms

  real(dp) :: magcon_lambda                  ! Strength of the atomic spherical constraint
  real(dp) :: ucvol                          ! Unit cell volume

!arrays

  integer :: ngfftf(18)                      ! Number of FFT grid points (for this processor) for the "fine" grid

  integer,allocatable :: typat(:)
  ! typat(natom)
  ! Type of each natom 

  integer,allocatable :: constraint_kind(:)
  ! constraint_kind(ntypat)
  ! Constraint kind to be applied to each type of atom. See corresponding input variable

  real(dp) :: gmet(3,3)
  ! Reciprocal space metric tensor, Bohr^2 units.

  real(dp) :: rprimd(3,3)
  ! Direct lattice vectors, Bohr units.

  real(dp),allocatable :: chrgat(:)
  ! chrgat(natom)
  ! Target charge for each atom. Not always used, it depends on the value of constraint_kind

  real(dp),allocatable :: intgf2(:)
  ! intgf2(natom)
  ! Integral of the square of the spherical integrating function, for each atom. 
  ! Initialized using some xred values, will not change during the SCF cycles, except for exotic algorithms, not in production,
  ! and even in this case, the change for different xred should be very small.

  real(dp),allocatable :: ratsph(:)
  ! ratsph(ntypat)
  ! Radius of the atomic sphere for each type of atom

  real(dp),allocatable :: spinat(:,:)
  ! spinat(3,natom)
  ! Target magnetization for each atom. Possibly only the direction or the magnitude, depending on constraint_kind

  real(dp),allocatable :: ziontypat(:)
  ! ziontypat(ntypat)
  ! Ionic charge, per type of atom
 
 end type constrained_dft_t

!!***

CONTAINS

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
  natom,nrx,nry,nrz,ntypat,rprimd,xcart,typat,zion,prtcharge,hcharge,hden,hweight)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nrx,nry,nrz,ntypat,prtcharge,mpoint
 real(dp),intent(in) :: minimal_den
!arrays
 integer,intent(in) :: typat(natom),npoint(ntypat)
 real(dp),intent(in) :: grid_den(nrx,nry,nrz),rprimd(3,3),zion(ntypat)
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

!!****f* m_dens/add_atomic_fcts
!! NAME
!! add_atomic_fcts
!!
!! FUNCTION
!! This routine is called to add the atomic spherical functions to some input function (usually an input potential residual).
!! The contributions from each atomic sphere are governed by parameters coeff_constr_dft, input to the present routine.
!!
!! INPUTS
!!  natom=number of atoms
!!  nspden = number of spin densities (1 2 or 4)
!!  rprimd=lattice vectors (dimensionful)
!!  mpi_enreg=mpi structure with communicator info
!!  nfft=number of points in standard fft grid
!!  ngfft=FFT grid dimensions
!!  ntypat=number of types of atoms
!!  ratsph(ntypat)=radii for muffin tin spheres of each atom
!!  typat(natom)=types of atoms
!!  xred(3,natom)=reduced atomic positions
!!
!! SIDE EFFECTS
!!  nv_constr_dft_r=the constrained potential or density in real space
!!
!! PARENTS
!!
!! CHILDREN
!!      metric,ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

subroutine add_atomic_fcts(natom,nspden,rprimd,mpi_enreg,nfft,ngfft,ntypat,ratsph, &
  typat,coeffs_constr_dft,nv_constr_dft_r,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,nspden,ntypat
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in)  :: typat(natom)
 integer,intent(in)  :: ngfft(18)
 real(dp),intent(in) :: coeffs_constr_dft(nspden,natom)
 real(dp),intent(inout) :: nv_constr_dft_r(nfft,nspden)
 real(dp),intent(in) :: ratsph(ntypat)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: xred(3,natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: ishift=5
 integer :: iatom, ierr
 integer :: n1a, n1b, n3a, n3b, n2a, n2b
 integer :: n1, n2, n3
 integer :: ifft_local
 integer ::  i1,i2,i3,ix,iy,iz,izloc
 real(dp) :: dify,difz,fsm,r2atsph,rr1,rr2,rr3,ratsm,ratsm2,rx23,ry23,rz23
 real(dp) :: r2,r2_11,r2_123,r2_23
 real(dp) :: ucvol
 real(dp),parameter :: delta=0.99_dp
!arrays
 real(dp), allocatable :: difx(:)
 real(dp) :: gprimd(3,3),rmet(3,3),gmet(3,3)
 real(dp) :: tsec(2)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)

! ***********************************************************************************************

!We need the metric because it is needed to compute the "box" around each atom
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 n1 = ngfft(1)
 n2 = ngfft(2)
 n3 = ngfft(3)

 ratsm = 0.05_dp ! default value for the smearing region radius - may become input variable later

 nv_constr_dft_r = zero

!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!Loop over atoms
!-------------------------------------------
 do iatom=1,natom

   if(sum(coeffs_constr_dft(1:nspden,iatom)**2)<tol12)then
     cycle
   endif

!  Define a "box" around the atom
   r2atsph=1.0000001_dp*ratsph(typat(iatom))**2
   rr1=sqrt(r2atsph*gmet(1,1))
   rr2=sqrt(r2atsph*gmet(2,2))
   rr3=sqrt(r2atsph*gmet(3,3))

   n1a=int((xred(1,iatom)-rr1+ishift)*n1+delta)-ishift*n1
   n1b=int((xred(1,iatom)+rr1+ishift)*n1      )-ishift*n1
   n2a=int((xred(2,iatom)-rr2+ishift)*n2+delta)-ishift*n2
   n2b=int((xred(2,iatom)+rr2+ishift)*n2      )-ishift*n2
   n3a=int((xred(3,iatom)-rr3+ishift)*n3+delta)-ishift*n3
   n3b=int((xred(3,iatom)+rr3+ishift)*n3      )-ishift*n3

   ratsm2 = -(ratsm**2 - 2*ratsph(typat(iatom))*ratsm)

   ABI_ALLOCATE(difx,(n1a:n1b))
   do i1=n1a,n1b
     difx(i1)=dble(i1)/dble(n1)-xred(1,iatom)
   enddo ! i1

   do i3=n3a,n3b
     iz=mod(i3+ishift*n3,n3)
     if(fftn3_distrib(iz+1)==mpi_enreg%me_fft) then
       izloc = ffti3_local(iz+1) - 1
       difz=dble(i3)/dble(n3)-xred(3,iatom)
       do i2=n2a,n2b
         iy=mod(i2+ishift*n2,n2)
         dify=dble(i2)/dble(n2)-xred(2,iatom)
         rx23=dify*rprimd(1,2)+difz*rprimd(1,3)
         ry23=dify*rprimd(2,2)+difz*rprimd(2,3)
         rz23=dify*rprimd(3,2)+difz*rprimd(3,3)
         r2_23=rx23**2+ry23**2+rz23**2
         r2_11=rprimd(1,1)**2+rprimd(2,1)**2+rprimd(3,1)**2
         r2_123=2*(rprimd(1,1)*rx23+rprimd(2,1)*ry23+rprimd(3,1)*rz23)
         do i1=n1a,n1b
           r2=(difx(i1)*r2_11+r2_123)*difx(i1)+r2_23
           if (r2 > r2atsph) cycle
           fsm = radsmear(r2, r2atsph, ratsm2)
           ix=mod(i1+ishift*n1,n1)
!          Identify the fft indexes of the rectangular grid around the atom
           ifft_local=1+ix+n1*(iy+n2*izloc)
           nv_constr_dft_r(ifft_local,1:nspden)=nv_constr_dft_r(ifft_local,1:nspden) + fsm*coeffs_constr_dft(1:nspden,iatom)

         end do  ! i1
       end do  ! i2
     end if  ! if this is my fft slice
   end do ! i3
   ABI_DEALLOCATE(difx)

!  end loop over atoms
 end do

!MPI parallelization
!TODO: test if xmpi_sum does the correct stuff for a slice of nv_constr_dft_r
 if(mpi_enreg%nproc_fft>1)then
   call timab(48,1,tsec)
   call xmpi_sum(nv_constr_dft_r,mpi_enreg%comm_fft,ierr)
   call timab(48,2,tsec)
 end if

! write (201,*) '# potential 1'
! write (201,*) nv_constr_dft_r(:,1)

! write (202,*) '# potential 2'
! write (202,*) nv_constr_dft_r(:,2)

! if (nspden > 2) then
!   write (203,*) '# potential 3'
!   write (203,*) nv_constr_dft_r(:,3)

!   write (204,*) '# potential 4'
!   write (204,*) nv_constr_dft_r(:,4)
! end if

end subroutine add_atomic_fcts
!!***

!!****f* m_dens/constrained_dft_ini
!! NAME
!! constrained_dft_ini
!!
!! FUNCTION
!! Initialize the constrained_dft datastructure.
!! Mostly copying already available (dtset) information, but also computing intgf2

!!
!! INPUTS
!!  constraint_kinds(ntypat)=for each type of atom, 0=no constraint,
!!    1=fix only the magnetization direction, following spinat direction,
!!    2=fix the magnetization vector to be the spinat one,
!!    3=fix the magnetization amplitude to be the spinat one, but does not fix its direction
!!    other future values will constrain the local atomic charge and possibly mix constraints if needed.
!!  magcon_lambda=strength of the atomic spherical constraint
!!  mpi_enreg=mpi structure with communicator info 
!!  natom=number of atoms
!!  nfft=number of points in standard fft grid
!!  ngfft=FFT grid dimensions
!!  nspden = number of spin densities (1 2 or 4)
!!  ntypat=number of types of atoms
!!  ratsph(ntypat)=radii for muffin tin spheres of each atom
!!  rprimd=lattice vectors (dimensioned)
!!  spinat(3,natom)=magnetic moments vectors, possible targets according to the value of constraint_kinds
!!  typat(natom)=types of atoms
!!  xred(3,natom)=reduced atomic positions
!!
!! OUTPUT
!!  constrained_dft=datastructure that contain the needed information to enforce the density and magnetization constraints
!!    Most of the data are simply copied from dtset, but also constrained_dft%intgf2(1:natom) is computed from the available data.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 subroutine constrained_dft_ini(constrained_dft,constraint_kind,magcon_lambda,mpi_enreg,natom,nfftf,ngfftf,nspden,ntypat,&
& ratsph,rprimd,spinat,typat,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: natom,nfftf,nspden,ntypat
 real(dp),intent(in) :: magcon_lambda
 type(MPI_type),intent(in) :: mpi_enreg
 type(constrained_dft_t),intent(out):: constrained_dft
!arrays
 integer,intent(in)  :: constraint_kind(ntypat)
 integer,intent(in)  :: ngfftf(18)
 integer,intent(in)  :: typat(natom)
 real(dp),intent(in) :: ratsph(ntypat)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: spinat(3,natom)
 real(dp),intent(in) :: xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: cplex1=1
 real(dp) :: ucvol
!arrays
 real(dp), allocatable :: intgf2(:) ! natom
 real(dp), allocatable :: rhor_dum(:,:) ! nfftf,nspden
 real(dp) :: gprimd(3,3),rmet(3,3),gmet(3,3)

! ***********************************************************************************************

 ABI_ALLOCATE(intgf2,(natom))
 ABI_ALLOCATE(rhor_dum,(nfftf,nspden))

!We need the metric because it is needed in calcdensph.F90
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!We need the integrated magnetic moments and the smoothing function
 call calcdensph(gmet,mpi_enreg,natom,nfftf,ngfftf,nspden,ntypat,std_out,&
&  ratsph,rhor_dum,rprimd,typat,ucvol,xred,0,cplex1,intgf2=intgf2)

 constrained_dft%gmet            =gmet
 constrained_dft%magcon_lambda   =magcon_lambda
 constrained_dft%natom           =natom
 constrained_dft%nfftf           =nfftf
 constrained_dft%ngfftf          =ngfftf
 constrained_dft%nspden          =nspden
 constrained_dft%ntypat          =ntypat
 constrained_dft%rprimd          =rprimd
 constrained_dft%ucvol           =ucvol

 ABI_ALLOCATE(constrained_dft%constraint_kind,(ntypat))
 ABI_ALLOCATE(constrained_dft%intgf2,(natom))
 ABI_ALLOCATE(constrained_dft%ratsph,(ntypat))
 ABI_ALLOCATE(constrained_dft%spinat,(3,natom))
 ABI_ALLOCATE(constrained_dft%typat,(natom))

 constrained_dft%constraint_kind=constraint_kind
 constrained_dft%intgf2=intgf2
 constrained_dft%ratsph=ratsph
 constrained_dft%spinat=spinat
 constrained_dft%typat=typat

 ABI_DEALLOCATE(intgf2)
 ABI_DEALLOCATE(rhor_dum)

end subroutine constrained_dft_ini
!!***


!!****f* m_dens/constrained_dft_free
!! NAME
!! constrained_dft_free
!!
!! FUNCTION
!! Free the constrained_dft datastructure.
!! Mostly copying already available (dtset) information, but also computing intgf2

!!
!! INPUTS
!!  constrained_dft=datastructure that contain the needed information to enforce the density and magnetization constraints
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 subroutine constrained_dft_free(constrained_dft)

!Arguments ------------------------------------
!scalars
 type(constrained_dft_t),intent(out):: constrained_dft

!Local variables-------------------------------

! ***********************************************************************************************

 ABI_SFREE(constrained_dft%constraint_kind)
 ABI_SFREE(constrained_dft%intgf2)
 ABI_SFREE(constrained_dft%ratsph)
 ABI_SFREE(constrained_dft%spinat)
 ABI_SFREE(constrained_dft%typat)

end subroutine constrained_dft_free
!!***


!!****f* m_dens/constrained_residual
!! NAME
!! constrained_residual
!!
!! FUNCTION
!! Recompute the residual to take into account the constraints, within constrained DFT.
!! The kind of constraint is given by constraint_kind, and the target values are given by spinat, for the local atomic magnetization,
!! and another future argument (possibly chrgat), for the local atomic charge.

!!
!! INPUTS
!!  c_dft <type(constrained_dft_t)>=datastructure for the information related to constrained DFT
!!   ! constraint_kind(ntypat)=for each type of atom, 0=no constraint, 
!!   !  1=fix only the magnetization direction, following spinat direction, 
!!   !  2=fix the magnetization vector to be the spinat one,
!!   !  3=fix the magnetization amplitude to be the spinat one, but does not fix its direction
!!   !  other future values will constrain the local atomic charge and possibly mix constraints if needed.
!!   ! intgf2(natom)=(precomputed) integral of the square of the spherical integration function for each atom in a sphere of radius ratsph.
!!   ! magcon_lambda=strength of the atomic spherical constraint
!!   ! natom=number of atoms
!!   ! nfftf=number of points in fine fft grid
!!   ! ngfftf=FFT grid dimensions
!!   ! nspden = number of spin densities (1 2 or 4)
!!   ! ntypat=number of types of atoms
!!   ! ratsph(ntypat)=radii for muffin tin spheres of each atom
!!   ! rprimd=lattice vectors (dimensioned)
!!   ! spinat(3,natom)=magnetic moments vectors, possible targets according to the value of constraint_kind
!!   ! typat(natom)=types of atoms
!!  mpi_enreg=mpi structure with communicator info
!!  rhor(nfft,nspden)=array for electron density in el./bohr**3. At output it will be constrained.
!!  xred(3,natom)=reduced atomic positions
!!
!! SIDE EFFECTS
!!  vresid(nfft,nspden)==array for potential residual in real space
!!    At output it will be modified: projected onto the space orthogonal to the atomic spherical functions (if there is a related
!!    constrained, and augmented by such atomic spherical functions multiplied by the difference between the actual
!!    integrated charge or magnetization and the target ones.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 subroutine constrained_residual(c_dft,mpi_enreg,rhor,vresid,xred)

!Arguments ------------------------------------
!scalars
 type(constrained_dft_t),intent(in) :: c_dft
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 real(dp),intent(in) :: rhor(c_dft%nfftf,c_dft%nspden)
 real(dp),intent(inout) :: vresid(c_dft%nfftf,c_dft%nspden)
 real(dp),intent(in) :: xred(3,c_dft%natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,natom,nfftf,nspden,ntypat
 integer :: cplex1=1
 real(dp) :: conkind,intgden_norm,intgden_proj,norm
!arrays
 real(dp), allocatable :: coeffs_constr_dft(:,:) ! nspden,natom
 real(dp), allocatable :: intgden(:,:) ! nspden,natom
 real(dp), allocatable :: intgres(:,:) ! nspden,natom
 real(dp) :: spinat_normed(3)

! ***********************************************************************************************

 natom=c_dft%natom
 nfftf=c_dft%nfftf
 nspden=c_dft%nspden
 ntypat=c_dft%ntypat

!We need the integrated magnetic moments 
 ABI_ALLOCATE(intgden,(nspden,natom))
 call calcdensph(c_dft%gmet,mpi_enreg,natom,nfftf,c_dft%ngfftf,nspden,ntypat,std_out,&
&  c_dft%ratsph,rhor,c_dft%rprimd,c_dft%typat,c_dft%ucvol,xred,1,cplex1,intgden=intgden)

!We need the integrated residuals
 ABI_ALLOCATE(intgres,(nspden,natom))
 call calcdensph(c_dft%gmet,mpi_enreg,natom,nfftf,c_dft%ngfftf,nspden,ntypat,std_out,&
&  c_dft%ratsph,vresid,c_dft%rprimd,c_dft%typat,c_dft%ucvol,xred,1,cplex1,intgden=intgres)

!The proper combination of intgden and intgres is stored in intgden: it is an effective charge/magnetization, to be compared to the target one.
 do iatom=1,natom

   !The total charge, conjugate to the average potential
   intgden(1,iatom)=c_dft%magcon_lambda * intgden(1,iatom) + half*(intgres(1,iatom)+intgres(2,iatom))/c_dft%intgf2(iatom)

   if(nspden==2)then

     !The magnetization along z
     intgden(2,iatom)=c_dft%magcon_lambda * (intgden(2,iatom)-half*intgden(1,iatom)) + (intgres(1,iatom)-intgres(2,iatom))/c_dft%intgf2(iatom)

   else if(nspden==4)then

     !The magnetization along z
     intgden(2,iatom)=c_dft%magcon_lambda * intgden(2,iatom) + (intgres(1,iatom)-intgres(2,iatom))/c_dft%intgf2(iatom)
     intgden(3,iatom)=c_dft%magcon_lambda * intgden(3,iatom) + intgres(3,iatom) / c_dft%intgf2(iatom)
     intgden(4,iatom)=c_dft%magcon_lambda * intgden(4,iatom) - intgres(4,iatom) / c_dft%intgf2(iatom)

   endif
 enddo

 ABI_ALLOCATE(coeffs_constr_dft,(nspden,natom))
 coeffs_constr_dft=zero

!Loop over atoms, to obtain the charge/magnetization correction.
 do iatom=1,natom

   conkind=c_dft%constraint_kind(c_dft%typat(iatom))

   if(conkind >=10)then

!    coeffs_constr_dft(1,iatom)=intgden(1,iatom)-c_dft%chrgat(iatom)

   endif

!HERE

   if(conkind==1 .and. nspden>1)then

     coeffs_constr_dft(2:nspden,iatom)=intgden(2:nspden,iatom)-c_dft%spinat(2:nspden,iatom)

   else if( (conkind==2 .or. conkind==3) .and. nspden>1)then

     norm = sqrt(sum(c_dft%spinat(:,iatom)**2))
     if (norm > tol10) then

       if( conkind==2 )then
         if(nspden==4)then
           !Fix the direction
           spinat_normed(:) = c_dft%spinat(:,iatom) / norm
           !Calculate the scalar product of the fixed mag. mom. vector and calculated mag. mom. vector
           !This is actually the size of the projection of the calc. mag. mom. vector on the fixed mag. mom. vector
           intgden_proj=spinat_normed(1)*intgden(2,iatom)+ &
&            spinat_normed(2)*intgden(3,iatom)+ &
&            spinat_normed(3)*intgden(4,iatom)
           coeffs_constr_dft(2:nspden,iatom)=intgden(2:nspden,iatom)-spinat_normed(2:nspden)*intgden_proj 
         else if(nspden==2)then
           !The direction must be correct, collinear, so no change.
           coeffs_constr_dft(2,iatom)=zero
         endif
       else if( conkind==3 )then
         intgden_norm = sqrt(sum(intgden(2:nspden,iatom)**2))
         coeffs_constr_dft(2:nspden,iatom)=(one-norm/intgden_norm)*intgden(2:nspden,iatom)
       endif

     else 
       !In this case, we set the atomic magnetization to zero.
       coeffs_constr_dft(2:nspden,iatom)=intgden(2:nspden,iatom)
     endif
       
   end if

 enddo ! iatom

!THERE IS LIKELY A MODIFICATION TO BE DONE IN ORDER TO TREAT PROPERLY THE DENSITY VS POTENTIAL UPDATE ...

!Now compute the new residual
 call add_atomic_fcts(natom,nspden,c_dft%rprimd,mpi_enreg,nfftf,c_dft%ngfftf,ntypat,&
&  c_dft%ratsph,c_dft%typat,coeffs_constr_dft,vresid,xred)

 ABI_DEALLOCATE(coeffs_constr_dft)
 ABI_DEALLOCATE(intgden)
 ABI_DEALLOCATE(intgres)

 end subroutine constrained_residual
!!***

!!****f* m_dens/mag_penalty
!! NAME
!! mag_penalty
!!
!! FUNCTION
!! This routine is called to compute the potential corresponding to constrained magnetic moments using the penalty function algorithm.
!!
!! INPUTS
!!  natom=number of atoms
!!  spinat=fixed magnetic moments vectors
!!  nspden = number of spin densities (1 2 or 4)
!!  magconon=constraining option (on/off); 1=fix only the direction, 2=fix the direction and size
!!  magcon_lambda=the size of the penalty terms
!!  rprimd=lattice vectors (dimensionful)
!!  mpi_enreg=mpi structure with communicator info
!!  nfft=number of points in standard fft grid
!!  ngfft=FFT grid dimensions
!!  ntypat=number of types of atoms
!!  ratsph=radii for muffin tin spheres of each atom
!!  rhor=density in real space
!!  typat=types of atoms
!!  xred=reduced atomic positions
!!
!! OUTPUT
!!  nv_constr_dft_r=the constrained potential
!!
!! PARENTS
!!      energy,rhotov,setvtr
!!
!! CHILDREN
!!      calcdensph,metric,ptabs_fourdp,timab,xmpi_sum
!!
!! NOTES
!!  based on html notes for the VASP implementation at
!!  http://cms.mpi.univie.ac.at/vasp/vasp/Constraining_direction_magnetic_moments.html
!!
!! SOURCE

subroutine mag_penalty(natom,spinat,nspden,magconon,magcon_lambda,rprimd, &
                      mpi_enreg,nfft,ngfft,ntypat,ratsph,rhor, &
                      typat,nv_constr_dft_r,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,magconon,nfft,nspden
 integer,intent(in) :: ntypat
 real(dp),intent(in) :: magcon_lambda
 real(dp),intent(out) :: nv_constr_dft_r(nfft,nspden)
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in)  :: typat(natom)
 integer,intent(in)  :: ngfft(18)
 real(dp),intent(in) :: ratsph(ntypat)
 real(dp),intent(in) :: rhor(nfft,nspden)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: spinat(3,natom)
 real(dp),intent(in) :: xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: iatom
 integer :: cplex1=1
 real(dp):: cmm_x,cmm_y,cmm_z
 real(dp) :: intgden_proj,norm,ucvol
!arrays
 real(dp), allocatable :: coeffs_constr_dft(:,:) ! nspden,natom
 real(dp), allocatable :: intgden(:,:) ! nspden,natom
 real(dp) :: spinat_norm(3,natom)
 real(dp) :: gprimd(3,3),rmet(3,3),gmet(3,3)

! ***********************************************************************************************

 ABI_ALLOCATE(coeffs_constr_dft,(nspden,natom))
 ABI_ALLOCATE(intgden,(nspden,natom))

!We need the metric because it is needed in calcdensph.F90
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!We need the integrated magnetic moments and the smoothing function
 call calcdensph(gmet,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,std_out,ratsph,rhor,rprimd,typat,ucvol,xred,1,cplex1,intgden=intgden)

!Loop over atoms
!-------------------------------------------
 do iatom=1,natom

   norm = sqrt(sum(spinat(:,iatom)**2))
   spinat_norm(:,iatom) = zero
   if (norm > tol10) then
     spinat_norm(:,iatom) = spinat(:,iatom) / norm
   else if (magconon == 1) then
!    if spinat = 0 and we are imposing the direction only, skip this atom
     cycle
   end if

!  Calculate the x- and y-components of the square bracket term
   cmm_x = zero
   cmm_y = zero
   cmm_z = zero
   intgden_proj = zero
   if (nspden == 4) then
     if (magconon==1) then
!      Calculate the scalar product of the fixed mag. mom. vector and calculated mag. mom. vector
!      This is actually the size of the projection of the calc. mag. mom. vector on the fixed mag. mom. vector
       intgden_proj=spinat_norm(1,iatom)*intgden(2,iatom)+ &
&        spinat_norm(2,iatom)*intgden(3,iatom)+ &
&        spinat_norm(3,iatom)*intgden(4,iatom)

       cmm_x=intgden(2,iatom)
       cmm_x=cmm_x-spinat_norm(1,iatom)*intgden_proj

       cmm_y=intgden(3,iatom)
       cmm_y=cmm_y-spinat_norm(2,iatom)*intgden_proj

     else if (magconon==2 .and. nspden == 4) then
       cmm_x=intgden(2,iatom)-spinat(1,iatom)
       cmm_y=intgden(3,iatom)-spinat(2,iatom)
     end if

!    Calculate the constraining potential for x- and y- components of the mag. mom. vector
!    Eric Bousquet has derived the relationship between spin components and potential spin matrix elements:
!    1 = up up     = +z
!    2 = down down = -z
!    3 = up down   = +x
!    4 = down up   = -y
     coeffs_constr_dft(3,iatom)= 2*magcon_lambda*cmm_x
     coeffs_constr_dft(4,iatom)=-2*magcon_lambda*cmm_y
   end if ! nspden 4

!  Calculate the z-component of the square bracket term
   if (magconon==1) then
     if (nspden == 4) then
       ! m_z - spinat_z * <m | spinat>
       cmm_z = intgden(4,iatom) - spinat_norm(3,iatom)*intgden_proj
     else if (nspden == 2) then
       ! this will be just a sign +/- : are we in the same direction as spinat_z?
       !    need something more continuous??? To make sure the gradient pushes the state towards FM/AFM?
       cmm_z = -sign(one, (intgden(1,iatom)-intgden(2,iatom))*spinat_norm(3,iatom))
     end if
   else if (magconon==2) then
     if (nspden == 4) then
       cmm_z=intgden(4,iatom)-spinat(3,iatom)
     else if (nspden == 2) then
       ! this is up spins - down spins - requested moment ~ 0
       ! EB: note that intgden comes from calcdensph, which, in nspden=2 case, returns
       ! intgden(1)=rho_up=n+m
       ! intgden(2)=rho_dn=n-m
       ! Then, is the following line be
       ! cmm_z=half*(intgden(1,iatom)-intgden(2,iatom)) - spinat(3,iatom)
       ! ??
       cmm_z=intgden(1,iatom)-intgden(2,iatom) - spinat(3,iatom)
     end if
   endif

!  Calculate the constraining potential for z-component of the mag. mom. vector
   coeffs_constr_dft(1,iatom)= 2*magcon_lambda*cmm_z
   coeffs_constr_dft(2,iatom)=-2*magcon_lambda*cmm_z

 enddo ! iatom
 
!Now compute the potential in real space
 call add_atomic_fcts(natom,nspden,rprimd,mpi_enreg,nfft,ngfft,ntypat,ratsph, &
   typat,coeffs_constr_dft,nv_constr_dft_r,xred)

 ABI_DEALLOCATE(coeffs_constr_dft)
 ABI_DEALLOCATE(intgden)

end subroutine mag_penalty
!!***

!!****f* m_dens/mag_penalty_e
!! NAME
!! mag_penalty_e
!!
!! FUNCTION
!! Compute the energy corresponding to constrained magnetic moments.
!!
!! INPUTS
!!  magconon=constraining option (on/off); 1=fix only the direction, 2=fix the direction and size
!!  spinat=fixed magnetic moments vectors
!!  magcon_lambda=the size of the penalty terms
!!
!! OUTPUT
!!  Epen=penalty contribution to the total energy corresponding to the constrained potential
!!  Econstr=???
!!  Eexp=???
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      calcdensph,metric,wrtout
!!
!! SOURCE


subroutine mag_penalty_e(magconon,magcon_lambda,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,ratsph,rhor,rprimd,spinat,typat,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,magconon,nspden,nfft,ntypat
 real(dp),intent(in) :: magcon_lambda
!arrays
 integer, intent(in) :: ngfft(18),typat(natom)
 real(dp),intent(in) :: spinat(3,natom), rprimd(3,3)
 real(dp),intent(in) :: ratsph(ntypat),rhor(nfft,nspden),xred(3,natom)
 type(MPI_type),intent(in) :: mpi_enreg

!Local variables-------------------------------
!scalars
 integer :: iatom,ii
 integer :: cplex1=1    ! dummy argument for calcdensphere
 real(dp) :: intgden_proj, Epen,Econstr,lVp, norm
!arrays
 real(dp) :: intmm(3), mag_1atom(3)
 real(dp), allocatable :: intgden(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),ucvol
 real(dp) :: spinat_norm(3,natom)
 character(len=500) :: msg

! *********************************************************************

!We need the metric because it is needed in calcdensph.F90
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 ABI_ALLOCATE (intgden, (nspden,natom))

!We need the integrated magnetic moments
 cplex1=1
 call calcdensph(gmet,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,std_out,ratsph,rhor,rprimd,typat,ucvol,xred,&
& 1,cplex1,intgden=intgden)

 Epen=0
 Econstr=0
 lVp=0

!Loop over atoms
!-------------------------------------------
 do iatom=1,natom

   norm = sqrt(sum(spinat(:,iatom)**2))
   spinat_norm(:,iatom) = zero
   if (norm > tol10) then
     spinat_norm(:,iatom) = spinat(:,iatom) / norm
   else if (magconon == 1) then
!    if spinat = 0 and we are imposing the direction only, skip this atom
     cycle
   end if
!  Calculate the scalar product of the fixed mag. mom. vector and calculated mag. mom. vector
!  This is actually the size of the projection of the calc. mag. mom. vector on the fixed mag. mom. vector

! for the collinear spin case, set up a fictitious 3D vector along z
   if (nspden == 4) then
     mag_1atom(1:3) = intgden(2:4,iatom)
   else if (nspden == 2) then
     mag_1atom = zero
     mag_1atom(3) = intgden(1,iatom)-intgden(2,iatom)
   end if

   intgden_proj = zero
   intmm = zero
!  Calculate the square bracket term
   if (magconon==1) then
     intgden_proj=spinat_norm(1,iatom)*mag_1atom(1)+ &
&     spinat_norm(2,iatom)*mag_1atom(2)+ &
&     spinat_norm(3,iatom)*mag_1atom(3)

     do ii=1,3
       intmm(ii)=mag_1atom(ii)-spinat_norm(ii,iatom)*intgden_proj
     end do

!    Calculate the energy Epen corresponding to the constraining potential
!    Econstr and lVp do not have a clear meaning (yet)
     Epen=Epen+magcon_lambda*(intmm(1)*intmm(1)+intmm(2)*intmm(2)+intmm(3)*intmm(3))
     Econstr=Econstr-magcon_lambda*(intmm(1)*mag_1atom(1)+intmm(2)*mag_1atom(2)+intmm(3)*mag_1atom(3))
     lVp=lVp+2*magcon_lambda*(intmm(1)*mag_1atom(1)+intmm(2)*mag_1atom(2)+intmm(3)*mag_1atom(3))

   else if (magconon==2) then
     do ii=1,3
       intmm(ii)=mag_1atom(ii)-spinat(ii,iatom)
     end do

!    Calculate the energy Epen corresponding to the constraining potential
!    Epen = -Econstr - lVp
!    Econstr = -M**2 + spinat**2
!    lVp = +2 M \cdot spinat
     Epen=Epen+magcon_lambda*(intmm(1)*intmm(1)+intmm(2)*intmm(2)+intmm(3)*intmm(3))
     Econstr=Econstr-magcon_lambda*(mag_1atom(1)*mag_1atom(1)+&
&     mag_1atom(2)*mag_1atom(2)+&
&     mag_1atom(3)*mag_1atom(3)) &
&     +magcon_lambda*(spinat(1,iatom)*spinat(1,iatom)+&
&     spinat(2,iatom)*spinat(2,iatom)+&
&     spinat(3,iatom)*spinat(3,iatom))
     lVp=lVp+2*magcon_lambda*(intmm(1)*mag_1atom(1)+intmm(2)*mag_1atom(2)+intmm(3)*mag_1atom(3))
   end if

   write(msg, *) 'atom             constraining magnetic field'
   call wrtout(std_out,msg,'COLL')
   write(msg, '(I3,A2,E12.5,A2,E12.5,A2,E12.5)') &
   iatom,'  ',magcon_lambda*intmm(1),'  ',magcon_lambda*intmm(2),'  ',magcon_lambda*intmm(3)
   call wrtout(std_out,msg,'COLL')

!  End loop over atoms
!  -------------------------------------------
 end do

!Printing
 write(msg, '(A17,E10.3)' ) ' magcon_lambda    = ',magcon_lambda
 call wrtout(std_out,msg,'COLL')
 write(msg, '(A17,E12.5)' ) ' Lagrange penalty = ',Epen
 call wrtout(std_out,msg,'COLL')
 write(msg, '(A17,E12.5)' ) ' E_constraint     = ',Econstr
 call wrtout(std_out,msg,'COLL')
 write(msg, '(A17,E12.5)' ) ' lVp = ',lVp
 call wrtout(std_out,msg,'COLL')

 ABI_DEALLOCATE (intgden)

end subroutine mag_penalty_e
!!***

!!****f* m_dens/calcdensph
!! NAME
!! calcdensph
!!
!! FUNCTION
!! Compute and print integral of total density inside spheres around atoms.
!!
!! INPUTS
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  ntypat=number of atom types
!!  nunit=number of the unit for printing
!!  prtopt = if 1, the default printing is on (to unit nunit), if -1, 2, 3, 4, special printing options, if 0 no printing.
!!  ratsph(ntypat)=radius of spheres around atoms
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!!   (total in first half and spin-up in second half if nspden=2)
!!   (total in first comp. and magnetization in comp. 2 to 4 if nspden=4)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  typat(natom)=type of each atom
!!  ucvol=unit cell volume in bohr**3
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  dentot(nspden)=integrated density (magnetization...) over full u.c. vol, optional argument
!!  intgden(nspden, natom)=integrated density (magnetization...) for each atom in a sphere of radius ratsph. Optional arg
!!    Note that when intgden is present, the definition of the spherical integration function changes, as it is smoothed.
!!  intgf2(natom)=integral of the square of the spherical integration function for each atom in a sphere of radius ratsph. Optional arg
!!    if present, the routine also checks that the spheres do not overlap
!!  Rest is printing
!!
!! PARENTS
!!      dfpt_scfcv,mag_penalty,mag_constr_e,outscfcv
!!
!! CHILDREN
!!      timab,wrtout,xmpi_sum
!!
!! SOURCE

subroutine calcdensph(gmet,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,nunit,ratsph,rhor,rprimd,typat,ucvol,xred,&
&    prtopt,cplex,intgden,dentot,intgf2)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in)        :: natom,nfft,nspden,ntypat,nunit
 real(dp),intent(in)       :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
 integer ,intent(in)       :: prtopt
 integer, intent(in)       :: cplex
!arrays
 integer,intent(in)  :: ngfft(18),typat(natom)
 real(dp),intent(in) :: gmet(3,3),ratsph(ntypat),rhor(cplex*nfft,nspden),rprimd(3,3)
 real(dp),intent(in) :: xred(3,natom)
!integer,intent(out),optional   :: atgridpts(nfft)
 real(dp),intent(out),optional  :: intgden(nspden,natom)
 real(dp),intent(out),optional  :: intgf2(natom)
 real(dp),intent(out),optional  :: dentot(nspden)
!Local variables ------------------------------
!scalars
 integer,parameter :: ishift=5
 integer :: ii,i1,i2,i3,iatom,ierr,ifft_local,ix,iy,iz,izloc,jatom,n1,n1a,n1b,n2,ifft
 integer :: n2a,n2b,n3,n3a,n3b,nfftot
 integer :: cmplex_den,jfft
 real(dp),parameter :: delta=0.99_dp
 real(dp) :: difx,dify,difz,dist_ij,r2,r2atsph,rr1,rr2,rr3,rx,ry,rz
 real(dp) :: fsm, ratsm, ratsm2
 real(dp) :: mag_coll   , mag_x, mag_y, mag_z ! EB
 real(dp) :: mag_coll_im, mag_x_im, mag_y_im, mag_z_im ! SPr
 real(dp) :: sum_mag, sum_mag_x,sum_mag_y,sum_mag_z,sum_rho_up,sum_rho_dn,sum_rho_tot ! EB
 real(dp) :: rho_tot, rho_tot_im
! real(dp) :: rho_up,rho_dn,rho_tot !EB
 logical   :: grid_found
 character(len=500) :: msg
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: intgden_(nspden,natom),tsec(2), my_xred(3, natom), xshift(3, natom)

! *************************************************************************

 n1=ngfft(1);n2=ngfft(2);n3=ngfft(3)
 nfftot=n1*n2*n3
 intgden_=zero

 ! This routine is not able to handle xred positions that are "far" from the
 ! first unit cell so wrap xred into [0, 1[ interval here.
 call wrap2_zero_one(xred, my_xred, xshift)

!Default value for the smearing region radius  - may become input variable later
 ratsm = zero
 if(present(intgden).or.present(intgf2)) then
   ratsm = 0.05_dp ! So, the presence of intgden as argument changes the default smearing
 end if

!If intgf2 present, check that the spheres do not overlap. 
 if(present(intgf2))then
   do iatom=1,natom
     do jatom=1,natom
       dist_ij=dist2(xred(:,iatom),xred(:,jatom),rprimd,1)
       if(dist_ij-ratsph(typat(iatom))-ratsph(typat(jatom))<tol10)then
         write(msg,'(a,a,a,i5,a,i5,a,f12.4,a,a,a,f12.4,a,f12.4)')&
         'In constrained DFT, the spheres around different atoms cannot overlap. See input var ratsph.',ch10,&
&        'It is found that for atoms ',iatom,' and ',jatom,', the distance is ',dist_ij,' Bohr, while ',ch10,&
&        'the radii of the spheres are ',ratsph(typat(iatom)),' and ',ratsph(typat(jatom))
         MSG_ERROR(msg)
       endif
     enddo
   enddo
 endif

!Get the distrib associated with this fft_grid
 grid_found=.false.

 if(n2 == mpi_enreg%distribfft%n2_coarse ) then

   if(n3== size(mpi_enreg%distribfft%tab_fftdp3_distrib)) then

     fftn3_distrib => mpi_enreg%distribfft%tab_fftdp3_distrib
     ffti3_local => mpi_enreg%distribfft%tab_fftdp3_local
     grid_found=.true.

   end if

 end if

 if(n2 == mpi_enreg%distribfft%n2_fine ) then

   if(n3 == size(mpi_enreg%distribfft%tab_fftdp3dg_distrib)) then

     fftn3_distrib => mpi_enreg%distribfft%tab_fftdp3dg_distrib
     ffti3_local => mpi_enreg%distribfft%tab_fftdp3dg_local
     grid_found = .true.

   end if

 end if
 if(.not.(grid_found)) then

   MSG_BUG("Unable to find an allocated distrib for this fft grid")

 end if

!Loop over atoms
!-------------------------------------------
 ii=0

 do iatom=1,natom

   !if (present(atgridpts)) then
   !  npts(iatom)=0           !SPr: initialize the number of grid points within atomic sphere around atom i
   !  ii=ii+1                 !SPr: initialize running index for constructing an array of grid point indexes
   !end if                    !     within atomic spheres

!  Define a "box" around the atom
   r2atsph=1.0000001_dp*ratsph(typat(iatom))**2
   rr1=sqrt(r2atsph*gmet(1,1))
   rr2=sqrt(r2atsph*gmet(2,2))
   rr3=sqrt(r2atsph*gmet(3,3))

   n1a=int((my_xred(1,iatom)-rr1+ishift)*n1+delta)-ishift*n1
   n1b=int((my_xred(1,iatom)+rr1+ishift)*n1      )-ishift*n1
   n2a=int((my_xred(2,iatom)-rr2+ishift)*n2+delta)-ishift*n2
   n2b=int((my_xred(2,iatom)+rr2+ishift)*n2      )-ishift*n2
   n3a=int((my_xred(3,iatom)-rr3+ishift)*n3+delta)-ishift*n3
   n3b=int((my_xred(3,iatom)+rr3+ishift)*n3      )-ishift*n3

   !This is the "width" of the zone of smearing, in term of the square of radius
   ratsm2 = (2*ratsph(typat(iatom))-ratsm)*ratsm

   do i3=n3a,n3b
     iz=mod(i3+ishift*n3,n3)

     if(fftn3_distrib(iz+1)==mpi_enreg%me_fft) then

       izloc = ffti3_local(iz+1) - 1
       difz=dble(i3)/dble(n3)-my_xred(3,iatom)
       do i2=n2a,n2b
         iy=mod(i2+ishift*n2,n2)
         dify=dble(i2)/dble(n2)-my_xred(2,iatom)
         do i1=n1a,n1b
           ix=mod(i1+ishift*n1,n1)

           difx=dble(i1)/dble(n1)-my_xred(1,iatom)
           rx=difx*rprimd(1,1)+dify*rprimd(1,2)+difz*rprimd(1,3)
           ry=difx*rprimd(2,1)+dify*rprimd(2,2)+difz*rprimd(2,3)
           rz=difx*rprimd(3,1)+dify*rprimd(3,2)+difz*rprimd(3,3)
           r2=rx**2+ry**2+rz**2


!          Identify the fft indexes of the rectangular grid around the atom
           if(r2 > r2atsph) then
             cycle
           end if

           fsm = radsmear(r2, r2atsph, ratsm2)

           ifft_local=1+ix+n1*(iy+n2*izloc)

           !if(present(atgridpts)) then
           !  ii=ii+1
           !  atgridpts(ii)=ifft_local    !SPr: save grid point index (dbg: to check whether ifft_local is a valid "global" index )
           !end if

           if(present(intgf2))then
!            intgden_(1,iatom)= integral of the square of the spherical integrating function
             intgf2(iatom)=intgf2(iatom)+fsm*fsm
           endif

           if(nspden==1) then
!            intgden_(1,iatom)= integral of total density
             intgden_(1,iatom)=intgden_(1,iatom)+fsm*rhor(ifft_local,1)
           elseif(nspden==2) then
!            intgden_(1,iatom)= integral of up density
!            intgden_(2,iatom)= integral of dn density
             intgden_(1,iatom)=intgden_(1,iatom)+fsm*rhor(ifft_local,2)
             intgden_(2,iatom)=intgden_(2,iatom)+fsm*rhor(ifft_local,1)-rhor(ifft_local,2)
           else
!            intgden_(1,iatom)= integral of total density
!            intgden_(2,iatom)= integral of magnetization, x-component
!            intgden_(3,iatom)= integral of magnetization, y-component
!            intgden_(4,iatom)= integral of magnetization, z-component
             intgden_(1,iatom)=intgden_(1,iatom)+fsm*rhor(ifft_local,1)
             intgden_(2,iatom)=intgden_(2,iatom)+fsm*rhor(ifft_local,2)
             intgden_(3,iatom)=intgden_(3,iatom)+fsm*rhor(ifft_local,3)
             intgden_(4,iatom)=intgden_(4,iatom)+fsm*rhor(ifft_local,4)
           end if

         end do
       end do
     end if
   end do

   if(present(intgf2))then
     intgf2(iatom)=intgf2(iatom)*ucvol/dble(nfftot)
   endif

   intgden_(:,iatom)=intgden_(:,iatom)*ucvol/dble(nfftot)

   !if (present(atgridpts)) then
   !  npts(iatom)=ii-1
   !  do is=1,iatom-1,1
   !    npts(iatom)=npts(iatom)-npts(is)-1
   !  end do
   !  atgridpts(ii-npts(iatom))=npts(iatom)    !SPr: save number of grid points around atom i
   !end if


!  End loop over atoms
!  -------------------------------------------

 end do

! EB  - Compute magnetization of the whole cell
! SPr - in case of complex density array set cmplex_den to one
 cmplex_den=0

 if(cplex==2) then

   cmplex_den=1

 end if

 if(nspden==2) then
   mag_coll=zero
   mag_coll_im=zero
   rho_tot=zero
   rho_tot_im=zero
   do ifft=1,nfft
     jfft=(cmplex_den+1)*ifft
!    rho_up=rho_up+rhor(ifft,2)
!    rho_dn=rho_dn+(rhor(ifft,2)-rhor(ifft,1))
     rho_tot=rho_tot+rhor(jfft-cmplex_den,1)             ! real parts of density and magnetization
     mag_coll=mag_coll+2*rhor(jfft-cmplex_den,2)-rhor(jfft-cmplex_den,1)

     rho_tot_im=rho_tot_im+rhor(jfft,1)                  ! imaginary parts of density and magnetization
     mag_coll_im=mag_coll_im+2*rhor(jfft,2)-rhor(jfft,1) ! in case of real array both are equal to corresponding real parts
   end do
   mag_coll=mag_coll*ucvol/dble(nfftot)
   rho_tot =rho_tot *ucvol/dble(nfftot)

   mag_coll_im=mag_coll_im*ucvol/dble(nfftot)
   rho_tot_im =rho_tot_im *ucvol/dble(nfftot)
!  rho_up=rho_up*ucvol/dble(nfftot)
!  rho_dn=rho_dn*ucvol/dble(nfftot)
!  rho_tot=rho_tot*ucvol/dble(nfftot)
 else if(nspden==4) then
   rho_tot=0
   rho_tot_im=0
   mag_x=0
   mag_y=0
   mag_z=0
   mag_x_im=0
   mag_y_im=0
   mag_z_im=0
   do ifft=1,nfft
     jfft=(cmplex_den+1)*ifft
     rho_tot=rho_tot+rhor(jfft-cmplex_den,1)
     mag_x=mag_x+rhor(jfft-cmplex_den,2)
     mag_y=mag_y+rhor(jfft-cmplex_den,3)
     mag_z=mag_z+rhor(jfft-cmplex_den,4)
     rho_tot_im=rho_tot_im+rhor(jfft,1)
     mag_x_im=mag_x_im+rhor(jfft,2)
     mag_y_im=mag_y_im+rhor(jfft,3)
     mag_z_im=mag_z_im+rhor(jfft,4)
   end do
   rho_tot=rho_tot*ucvol/dble(nfftot)
   mag_x=mag_x*ucvol/dble(nfftot)
   mag_y=mag_y*ucvol/dble(nfftot)
   mag_z=mag_z*ucvol/dble(nfftot)
   rho_tot_im=rho_tot_im*ucvol/dble(nfftot)
   mag_x_im=mag_x_im*ucvol/dble(nfftot)
   mag_y_im=mag_y_im*ucvol/dble(nfftot)
   mag_z_im=mag_z_im*ucvol/dble(nfftot)
 end if

!MPI parallelization
 if(mpi_enreg%nproc_fft>1)then
   call timab(48,1,tsec)
   call xmpi_sum(intgden_,mpi_enreg%comm_fft,ierr)
   call xmpi_sum(rho_tot,mpi_enreg%comm_fft,ierr)  ! EB
   call xmpi_sum(mag_coll,mpi_enreg%comm_fft,ierr) ! EB

   call xmpi_sum(rho_tot_im,mpi_enreg%comm_fft,ierr)
   call xmpi_sum(mag_coll_im,mpi_enreg%comm_fft,ierr)
!  call xmpi_sum(rho_up,mpi_enreg%comm_fft,ierr)  ! EB
!  call xmpi_sum(rho_dn,mpi_enreg%comm_fft,ierr)  ! EB
   call xmpi_sum(mag_x,mpi_enreg%comm_fft,ierr)    ! EB
   call xmpi_sum(mag_y,mpi_enreg%comm_fft,ierr)    ! EB
   call xmpi_sum(mag_z,mpi_enreg%comm_fft,ierr)    ! EB
   call xmpi_sum(mag_x_im,mpi_enreg%comm_fft,ierr)    ! EB
   call xmpi_sum(mag_y_im,mpi_enreg%comm_fft,ierr)    ! EB
   call xmpi_sum(mag_z_im,mpi_enreg%comm_fft,ierr)    ! EB
   call timab(48,2,tsec)
 end if

!Printing
 sum_mag=zero
 sum_mag_x=zero
 sum_mag_y=zero
 sum_mag_z=zero
 sum_rho_up=zero
 sum_rho_dn=zero
 sum_rho_tot=zero

 if(prtopt==1) then

   if(nspden==1) then
     write(msg, '(4a)' ) &
&     ' Integrated electronic density in atomic spheres:',ch10,&
&     ' ------------------------------------------------'
     call wrtout(nunit,msg,'COLL')
     write(msg, '(a)' ) ' Atom  Sphere_radius  Integrated_density'
     call wrtout(nunit,msg,'COLL')
     do iatom=1,natom
       write(msg, '(i5,f15.5,f20.8)' ) iatom,ratsph(typat(iatom)),intgden_(1,iatom)
       call wrtout(nunit,msg,'COLL')
     end do
   elseif(nspden==2) then
     write(msg, '(4a)' ) &
&     ' Integrated electronic and magnetization densities in atomic spheres:',ch10,&
&     ' ---------------------------------------------------------------------'
     call wrtout(nunit,msg,'COLL')
     write(msg, '(3a)' ) ' Note: Diff(up-dn) is a rough ',&
&     'approximation of local magnetic moment'
     call wrtout(nunit,msg,'COLL')
     write(msg, '(a)' ) ' Atom    Radius    up_density   dn_density  Total(up+dn)  Diff(up-dn)'
     call wrtout(nunit,msg,'COLL')
     do iatom=1,natom
       write(msg, '(i5,f10.5,2f13.6,a,f12.6,a,f12.6)' ) iatom,ratsph(typat(iatom)),intgden_(1,iatom),intgden_(2,iatom),&
&       '  ',(intgden_(1,iatom)+intgden_(2,iatom)),' ',(intgden_(1,iatom)-intgden_(2,iatom))
       call wrtout(nunit,msg,'COLL')
       ! Compute the sum of the magnetization
       sum_mag=sum_mag+intgden_(1,iatom)-intgden_(2,iatom)
       sum_rho_up=sum_rho_up+intgden_(1,iatom)
       sum_rho_dn=sum_rho_dn+intgden_(2,iatom)
       sum_rho_tot=sum_rho_tot+intgden_(1,iatom)+intgden_(2,iatom)
     end do
     write(msg, '(a)') ' ---------------------------------------------------------------------'
     call wrtout(nunit,msg,'COLL')
     write(msg, '(a,2f13.6,a,f12.6,a,f12.6)') '  Sum:         ', sum_rho_up,sum_rho_dn,'  ',sum_rho_tot,' ',sum_mag
     call wrtout(nunit,msg,'COLL')
     write(msg, '(a,f14.6)') ' Total magnetization (from the atomic spheres):       ', sum_mag
     call wrtout(nunit,msg,'COLL')
     write(msg, '(a,f14.6)') ' Total magnetization (exact up - dn):                 ', mag_coll
     call wrtout(nunit,msg,'COLL')
   ! EB for testing purpose print rho_up, rho_dn and rho_tot
!    write(msg, '(a,3f14.4,2i8)') ' rho_up, rho_dn, rho_tot, nfftot, nfft: ', rho_up,rho_dn,rho_tot,nfft,nfft
!   call wrtout(nunit,msg,'COLL')

   elseif(nspden==4) then

     write(msg, '(4a)' ) &
&     ' Integrated electronic and magnetization densities in atomic spheres:',ch10,&
&     ' ---------------------------------------------------------------------'
     call wrtout(nunit,msg,'COLL')
     write(msg, '(3a)' ) ' Note:      this is a rough approximation of local magnetic moments'
     call wrtout(nunit,msg,'COLL')
     write(msg, '(a)' ) ' Atom   Radius      Total density     mag(x)      mag(y)      mag(z)  '
     call wrtout(nunit,msg,'COLL')
     do iatom=1,natom
       write(msg, '(i5,f10.5,f16.6,a,3f12.6)' ) iatom,ratsph(typat(iatom)),intgden_(1,iatom),'  ',(intgden_(ix,iatom),ix=2,4)
       call wrtout(nunit,msg,'COLL')
       ! Compute the sum of the magnetization in x, y and z directions
       sum_mag_x=sum_mag_x+intgden_(2,iatom)
       sum_mag_y=sum_mag_y+intgden_(3,iatom)
       sum_mag_z=sum_mag_z+intgden_(4,iatom)
     end do
     write(msg, '(a)') ' ---------------------------------------------------------------------'
     call wrtout(nunit,msg,'COLL')
!    write(msg, '(a,f12.6,f12.6,f12.6)') ' Total magnetization :           ', sum_mag_x,sum_mag_y,sum_mag_z
     write(msg, '(a,f12.6,f12.6,f12.6)') ' Total magnetization (spheres)   ', sum_mag_x,sum_mag_y,sum_mag_z
     call wrtout(nunit,msg,'COLL')
     write(msg, '(a,f12.6,f12.6,f12.6)') ' Total magnetization (exact)     ', mag_x,mag_y,mag_z
     call wrtout(nunit,msg,'COLL')
!    SPr for dfpt debug
!    write(msg, '(a,f12.6)') ' Total density (exact)           ', rho_tot
!   call wrtout(nunit,msg,'COLL')
   end if

 elseif(prtopt==-1) then

   write(msg, '(2a)') ch10,' ------------------------------------------------------------------------'
   call wrtout(nunit,msg,'COLL')

   if(nspden==1) then
     write(msg, '(4a)' ) &
&     ' Fermi level charge density n_f:',ch10,&
&     ' ------------------------------------------------------------------------',ch10
   else
     write(msg, '(4a)' ) &
&     ' Fermi level charge density n_f and magnetization m_f:',ch10,&
&     ' ------------------------------------------------------------------------',ch10
   end if
   call wrtout(nunit,msg,'COLL')

   if(cmplex_den==0) then
     write(msg, '(a,f13.8)') '     n_f   = ',rho_tot
   else
     write(msg, '(a,f13.8,a,f13.8)') '  Re[n_f]= ', rho_tot,"   Im[n_f]= ",rho_tot_im
   end if
   call wrtout(nunit,msg,'COLL')
   if(nspden==2) then
     if(cmplex_den==0) then
       write(msg, '(a,f13.8)') '     m_f    = ', mag_coll
     else
       write(msg, '(a,f13.8,a,f13.8)') '  Re[m_f]= ', mag_coll,"   Im[m_f]= ",mag_coll_im
     end if
     call wrtout(nunit,msg,'COLL')
   elseif (nspden==4) then
     write(msg, '(a,f13.8)') '     mx_f  = ',mag_x
     call wrtout(nunit,msg,'COLL')
     write(msg, '(a,f13.8)') '     my_f  = ',mag_y
     call wrtout(nunit,msg,'COLL')
     write(msg, '(a,f13.8)') '     mz_f  = ',mag_z
     call wrtout(nunit,msg,'COLL')
   end if

   write(msg, '(3a)') ch10,' ------------------------------------------------------------------------',ch10
   call wrtout(nunit,msg,'COLL')


 else if (prtopt==2 .or. prtopt==3 .or. prtopt==4) then
   ! Used in the DFPT case, prtopt=idir+1

   if(abs(rho_tot)<1.0d-10) then
     rho_tot=0
   end if

   write(msg, '(2a)') ch10,' ------------------------------------------------------------------------'
   call wrtout(nunit,msg,'COLL')

   if(nspden==1) then
     write(msg, '(4a)' ) &
&     ' Integral of the first order density n^(1):',ch10,&
&     ' ------------------------------------------------------------------------',ch10
   else
     write(msg, '(4a)' ) &
&     ' Integrals of the first order density n^(1) and magnetization m^(1):',ch10,&
&     ' ------------------------------------------------------------------------',ch10
   end if
   call wrtout(nunit,msg,'COLL')

   if(cmplex_den==0) then
     write(msg, '(a,e16.8)') '     n^(1)    = ', rho_tot
   else
     write(msg, '(a,e16.8,a,e16.8)') '  Re[n^(1)] = ', rho_tot,"   Im[n^(1)] = ",rho_tot_im
   end if
   call wrtout(nunit,msg,'COLL')

   if(nspden==2) then

     if(cmplex_den==0) then
       write(msg, '(a,e16.8)') '     m^(1)    = ', mag_coll
     else
       write(msg, '(a,e16.8,a,e16.8)') '  Re[m^(1)] = ', mag_coll,"   Im[m^(1)] = ",mag_coll_im
     end if
     call wrtout(nunit,msg,'COLL')

   elseif (nspden==4) then
     if(cmplex_den==0) then
       write(msg, '(a,e16.8)') '     mx^(1)   = ', mag_x
       call wrtout(nunit,msg,'COLL')
       write(msg, '(a,e16.8)') '     my^(1)   = ', mag_y
       call wrtout(nunit,msg,'COLL')
       write(msg, '(a,e16.8)') '     mz^(1)   = ', mag_z
       call wrtout(nunit,msg,'COLL')
     else
       write(msg, '(a,e16.8,a,e16.8)') '  Re[mx^(1)]= ',  mag_x, "   Im[mx^(1)]= ", mag_x_im
       call wrtout(nunit,msg,'COLL')
       write(msg, '(a,e16.8,a,e16.8)') '  Re[my^(1)]= ',  mag_y, "   Im[my^(1)]= ", mag_y_im
       call wrtout(nunit,msg,'COLL')
       write(msg, '(a,e16.8,a,e16.8)') '  Re[mz^(1)]= ',  mag_z, "   Im[mz^(1)]= ", mag_z_im
       call wrtout(nunit,msg,'COLL')
     end if
   end if

   write(msg, '(3a)') ch10,' ------------------------------------------------------------------------',ch10
   call wrtout(nunit,msg,'COLL')

 end if

 if(present(intgden)) then
   intgden = intgden_
 end if

 if(present(dentot)) then
   if(nspden==2) then
     dentot(1)=rho_tot
     dentot(2)=mag_coll
   elseif(nspden==4) then
     dentot(1)=rho_tot
     dentot(2)=mag_x
     dentot(3)=mag_y
     dentot(4)=mag_z
   end if
 end if

end subroutine calcdensph
!!***

!!****f* m_dens/radsmear
!! NAME
!! radsmear
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function radsmear(r, rsph, rsm)

!Arguments ------------------------------------
!scalars
 real(dp) :: radsmear
 real(dp), intent(in) :: r, rsph, rsm

!Local variables ------------------------------
!scalars
 real(dp) :: xx

!******************************************************************

 radsmear = zero
 if (r < rsph - rsm - tol12) then
   radsmear = one
 else if (r < rsph - tol12) then
   xx = (rsph - r) / rsm
   radsmear = xx**2*(3+xx*(1+xx*(-6+3*xx)))
 end if

end function radsmear
!!***

!!****f* ABINIT/printmagvtk
!! NAME
!!  printmagvtk
!!
!! FUNCTION
!!  Auxiliary routine for printing out magnetization density in VTK format.
!!  Output file name is DEN.vtk
!!
!! INPUTS
!!  mpi_enreg = information about adopted parallelization strategy
!!  nspden    = number of components of density matrix (possible values re 1,2, or 4)
!!              nspden:   1 -> rho
!!              nspden:   2 -> rho_up,rho_dwn
!!              nspden:   4 -> rho,mx,my,mz
!!  nfft      = number of fft points per FFT processor
!!  ngfft     = full information about FFT mesh
!!  rhor      = density array stored in the memory of current FFT processor
!!  rprimd    = array of lattice vectors
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  At the moment this routine is mainly used for development and debugging
!!  of gs and dfpt calculations with non-collinear spins. If needed, can be used
!!  to print final density in vtk format.
!!  IMPORTANT: implementation is thoroughly checked only for npspinor = 1,
!!             for other case might need to change the part gathering
!!             the FFT mesh info
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_gather
!!
!! SOURCE

subroutine printmagvtk(mpi_enreg,cplex,nspden,nfft,ngfft,rhor,rprimd,fname)

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(in)   :: mpi_enreg
 integer,intent(in)          :: nfft,nspden,cplex
!arrays
 integer,intent(in)          :: ngfft(18)
 real(dp),intent(in)         :: rhor(cplex*nfft,nspden),rprimd(3,3)
 character(len=*),intent(in) :: fname

!Local variables-------------------------------
!scalars
 integer :: denvtk,denxyz,denxyz_im,nfields
 integer :: nx,ny,nz,nfft_tot
 integer :: ii,jj,kk,ind,ispden
 integer :: mpi_comm,mpi_head,mpi_rank,ierr
 real(dp)    :: rx,ry,rz
 integer :: nproc_fft,ir
 character(len=500) :: msg
 character(len=10)  :: outformat
 character(len=50)   :: fname_vtk
 character(len=50)   :: fname_xyz
 character(len=50)   :: fname_xyz_re
 character(len=50)   :: fname_xyz_im
!arrays
 real(dp),allocatable :: rhorfull(:,:)


! *************************************************************************

 DBG_ENTER("COLL")

! if (option/=1 .and. option/=2 ) then
!   write(msg,'(3a,i0)')&
!&   'The argument option should be 1 or 2,',ch10,&
!&   'however, option=',option
!   MSG_BUG(msg)
! end if
!
! if (sizein<1) then
!   write(msg,'(3a,i0)')&
!&   'The argument sizein should be a positive number,',ch10,&
!&   'however, sizein=',sizein
!   MSG_ERROR(msg)
! end if

 DBG_EXIT("COLL")

 fname_vtk=adjustl(adjustr(fname)//".vtk")
 fname_xyz=adjustl(adjustr(fname)//".xyz")
 fname_xyz_re=adjustl(adjustr(fname)//"_re.xyz")
 fname_xyz_im=adjustl(adjustr(fname)//"_im.xyz")
 !write(std_out,*) ' Writing out .vtk file: ',fname_vtk
 !write(std_out,*) ' Writing out .xyz file: ',fname_xyz

  !if 1 or two component density then write out either 1 or 2 scalar density fields
  !if 4, then write one scalar field (density) and one vector field (magnetization density)
 if(nspden/=4)then
   nfields=nspden
 else
   nfields=2
 end if

 nfields=nfields*cplex

  ! FFT mesh specifications: full grid
 nx=ngfft(1)        ! number of points along 1st lattice vector
 ny=ngfft(2)        ! number of points along 2nd lattice vector
 nz=ngfft(3)        ! number of points along 3rd lattice vector
 nfft_tot=nx*ny*nz  ! total number of fft mesh points (can be different from nfft in case of distributed memory of nproc_fft processors)


  ! Gather information about memory distribution
 mpi_head=0
 mpi_comm = mpi_enreg%comm_fft
 mpi_rank = xmpi_comm_rank(mpi_comm)
 nproc_fft=ngfft(10)

  ! Create array to host full FFT mesh
 if(mpi_rank==mpi_head)then
   ABI_ALLOCATE(rhorfull,(cplex*nfft_tot,nspden))
 end if

  ! Fill in the full mesh
 if(nproc_fft==1)then
   rhorfull=rhor
 else
   do ir=1,nspden
     call xmpi_gather(rhor(:,ir),cplex*nfft,rhorfull(:,ir),cplex*nfft,mpi_head,mpi_comm,ierr)
   end do
 end if

 if(mpi_rank==mpi_head)then

    ! Open the output vtk file
   if (open_file(fname_vtk,msg,newunit=denvtk,status='replace',form='formatted') /=0) then
     MSG_WARNING(msg)
     RETURN
   end if

   if(cplex==1) then
     if (open_file(fname_xyz,msg,newunit=denxyz,status='replace',form='formatted') /=0) then
       MSG_WARNING(msg)
       RETURN
     end if
   else if (cplex==2) then
     if (open_file(fname_xyz_re,msg,newunit=denxyz,status='replace',form='formatted') /=0) then
       MSG_WARNING(msg)
       RETURN
     end if
     if (open_file(fname_xyz_im,msg,newunit=denxyz_im,status='replace',form='formatted') /=0) then
       MSG_WARNING(msg)
       RETURN
     end if
   end if

    ! Write the header of the output vtk file
   write(denvtk,"(a)") '# vtk DataFile Version 2.0'
   write(denvtk,"(a)") 'Electron density components'
   write(denvtk,"(a)") 'ASCII'
   write(denvtk,"(a)") 'DATASET STRUCTURED_GRID'
   write(denvtk,"(a,3i6)") 'DIMENSIONS ', nx,ny,nz
   write(denvtk,"(a,i18,a)") 'POINTS ',nfft_tot,' double'

   if (nspden==1) then
     outformat="(4e16.8)"
   else if (nspden==2) then
     outformat="(5e16.8)"
   else
     outformat="(7e16.8)"
   end if

    ! Write out information about grid points
   do kk=0,nz-1
     do jj=0,ny-1
       do ii=0,nx-1

         rx=(dble(ii)/nx)*rprimd(1,1)+(dble(jj)/ny)*rprimd(1,2)+(dble(kk)/nz)*rprimd(1,3)
         ry=(dble(ii)/nx)*rprimd(2,1)+(dble(jj)/ny)*rprimd(2,2)+(dble(kk)/nz)*rprimd(2,3)
         rz=(dble(ii)/nx)*rprimd(3,1)+(dble(jj)/ny)*rprimd(3,2)+(dble(kk)/nz)*rprimd(3,3)
         write(denvtk,'(3f16.8)') rx,ry,rz  !coordinates of the grid point
         ind=1+ii+nx*(jj+ny*kk)
         if (cplex==1) then
           write(denxyz,outformat) rx,ry,rz,(rhorfull(ind,ispden),ispden=1,nspden)
         else
           write(denxyz,outformat)    rx,ry,rz,(rhorfull(2*ind-1,ispden),ispden=1,nspden)
           write(denxyz_im,outformat) rx,ry,rz,(rhorfull(2*ind  ,ispden),ispden=1,nspden)
         end if
       end do
     end do
   end do

   if(cplex==1) then
     close(denxyz)
   else
     close(denxyz)
     close(denxyz_im)
   end if

    ! Write out information about field defined on the FFT mesh
   write(denvtk,"(a,i18)") 'POINT_DATA ',nfft_tot
   write(denvtk,"(a,i6)")  'FIELD Densities ',nfields


    ! Write out different fields depending on the number of density matrix components
   if(nspden==1)then

      !single component, so just write out the density
     if(cplex==1) then
       write(denvtk,"(a,i18,a)") 'rho 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=1+ii+nx*(jj+ny*kk)
             write(denvtk,'(f16.8)') rhorfull(ind,1)
           end do
         end do
       end do
     else
       write(denvtk,"(a,i18,a)") 'Re_rho 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))-1
             write(denvtk,'(f16.8)') rhorfull(ind,1)
           end do
         end do
       end do
       write(denvtk,"(a,i18,a)") 'Im_rho 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))
             write(denvtk,'(f16.8)') rhorfull(ind,1)
           end do
         end do
       end do
     end if

   else if(nspden==2)then

      !two component, write the density for spin_up and spin_down channels
     if(cplex==1) then

       write(denvtk,"(a,i18,a)") 'rho 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=1+ii+nx*(jj+ny*kk)
             write(denvtk,'(f16.8)') rhorfull(ind,1)
           end do
         end do
       end do
       write(denvtk,"(a,i18,a)") 'mag 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=1+ii+nx*(jj+ny*kk)
             write(denvtk,'(f16.8)') 2*rhorfull(ind,2)-rhorfull(ind,1)
           end do
         end do
       end do

     else

       write(denvtk,"(a,i18,a)") 'Re_rho 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))-1
             write(denvtk,'(f16.8)') rhorfull(ind,1)
           end do
         end do
       end do
       write(denvtk,"(a,i18,a)") 'Im_rho 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))
             write(denvtk,'(f16.8)') rhorfull(ind,1)
           end do
         end do
       end do
       write(denvtk,"(a,i18,a)") 'Re_mag 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))-1
             write(denvtk,'(f16.8)') 2*rhorfull(ind,2)-rhorfull(ind,1)
           end do
         end do
       end do
       write(denvtk,"(a,i18,a)") 'Im_mag 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))
             write(denvtk,'(f16.8)') 2*rhorfull(ind,2)-rhorfull(ind,1)
           end do
         end do
       end do

     end if

   else  !here is the last option: nspden==4

     if(cplex==1) then

        !four component, write the density (scalar field) and magnetization density (vector field)
       write(denvtk,"(a,i18,a)") 'rho 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=1+ii+nx*(jj+ny*kk)
             write(denvtk,'(f16.8)') rhorfull(ind,1)
           end do
         end do
       end do
       write(denvtk,"(a,i18,a)") 'mag 3 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=1+ii+nx*(jj+ny*kk)
             write(denvtk,'(3f16.8)') rhorfull(ind,2),rhorfull(ind,3),rhorfull(ind,4)
           end do
         end do
       end do

     else

       write(denvtk,"(a,i18,a)") 'Re_rho 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))-1
             write(denvtk,'(f16.8)') rhorfull(ind,1)
           end do
         end do
       end do
       write(denvtk,"(a,i18,a)") 'Im_rho 1 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))
             write(denvtk,'(f16.8)') rhorfull(ind,1)
           end do
         end do
       end do
       write(denvtk,"(a,i18,a)") 'Re_mag 3 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))-1
             write(denvtk,'(3f16.8)') rhorfull(ind,2),rhorfull(ind,3),rhorfull(ind,4)
           end do
         end do
       end do
       write(denvtk,"(a,i18,a)") 'Im_mag 3 ',nfft_tot,' double'
       do kk=0,nz-1
         do jj=0,ny-1
           do ii=0,nx-1
             ind=2*(1+ii+nx*(jj+ny*kk))
             write(denvtk,'(3f16.8)') rhorfull(ind,2),rhorfull(ind,3),rhorfull(ind,4)
           end do
         end do
       end do

     end if

   end if ! nspden options condition

   close (denvtk)

    !clean up the gathered FFT mesh
   ABI_DEALLOCATE(rhorfull)

 end if

end subroutine printmagvtk
!!***

end module m_dens
!!***
