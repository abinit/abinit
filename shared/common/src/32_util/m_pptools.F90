!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_pptools
!! NAME
!! m_pptools
!!
!! FUNCTION
!!  Helper functions used for post-processing.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2019 ABINIT group (MG, ZL, MJV, BXu)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

MODULE m_pptools

 use defs_basis
 use m_errors
 use m_abicore
 use m_krank

 use m_io_tools,        only : open_file
 use m_fstrings,        only : sjoin, itoa
 use m_numeric_tools,   only : wrap2_pmhalf

 implicit none

 private

 public :: prmat                ! print real(dp) matrices in an attractive format.
 public :: printxsf             ! Write a generic array in the XSF format (XCrysden format)
 public :: print_fofr_ri        ! Print the [real, imaginary] part of an array
 public :: print_fofr_xyzri     ! Print the Cartesian coordinates and the [real,imaginary] part of an array
 public :: print_fofr_cube      ! Print ||fofr|| in CUBE format.
 public :: printbxsf            ! Print band structure energies in XCrysDen format.
 public :: printvtk             ! Print band structure energies and velocities in VTK format.

CONTAINS  !===========================================================
!!***

!!****f* m_pptools/prmat
!! NAME
!! prmat
!!
!! FUNCTION
!! This subroutine prints real*8 matrices in an attractive format.
!!
!! INPUTS
!!  mat(mi,nj)= matrix to be printed
!!  mi        = no rows of mat
!!  ni        = no rows to print
!!  nj        = no colums of mat
!!  unitm     = unit to print to, if not provided std_out is chosen
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      chiscwrt,ioniondist,newkpt,pawuj_det,pawuj_red,pawuj_utils,shellstruct
!!
!! CHILDREN
!!      wrap2_pmhalf
!!
!! SOURCE

subroutine prmat (mat, ni, nj, mi, unitm)


!Arguments ------------------------------------
!scalars
 integer,intent(in)           :: mi,ni,nj
 integer,intent(in), optional :: unitm
!arrays
 real(dp),intent(in)          :: mat(mi,nj)

!Local variables-------------------------------
!scalars
 character(len=1000)    :: message
 integer,parameter      :: nline=10
 integer                :: ii,jj,jstart,jstop,unitn

! *************************************************************************

 if (present(unitm)) then ! standard printing to std_out
   unitn=unitm
 else
   unitn=std_out
 end if

 do  jstart = 1, nj, nline
   jstop = min(nj, jstart+nline-1)
   write(message, '(3x,10(i4,8x))' ) (jj,jj=jstart,jstop)
   call wrtout(unitn,message,'COLL')
 end do

 do ii = 1,ni
   do jstart= 1, nj, nline
     jstop = min(nj, jstart+nline-1)
     if (jstart==1) then
       write(message, '(i3,1p,10e12.4)' ) ii, (mat(ii,jj),jj=jstart,jstop)
       call wrtout(unitn,message,'COLL')
     else
       write(message, '(3x,1p,10e12.4)' )    (mat(ii,jj),jj=jstart,jstop)
       call wrtout(unitn,message,'COLL')
     end if
   end do
 end do

end subroutine prmat
!!***

!----------------------------------------------------------------------

!!****f* m_pptools/printxsf
!! NAME
!! printxsf
!!
!! FUNCTION
!! Write a generic array in the XSF format (XCrysden format)
!!
!! INPUTS
!! basis(3,3) = basis vectors of the direct real lattice or of the reciprocal lattice (fortran convention)
!!              (Bohr units if realrecip=0, Bohr^-1 if realrecip=1, see below)
!! realrecip = 0  for a plot in real space
!!             1  for a plot in reciprocal space
!! nunit   = unit number of the output file (already open by the caller, not closed here!)
!! n1=grid size along x
!! n2=grid size along y
!! n3=grid size along z
!! origin(3) = origin of the grid
!! datagrid(n1*n2*n3) = datagrid values stored using the fortran convention
!!
!! OUTPUT
!! Only write
!!
!! PARENTS
!!      denfgr,exc_plot,m_kxc,m_nesting,m_wfd,spin_current
!!
!! CHILDREN
!!      wrap2_pmhalf
!!
!! SOURCE

subroutine printxsf(n1,n2,n3,datagrid,basis,origin,natom,ntypat,typat,xcart,znucl,nunit,realrecip)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1,n2,n3,nunit,realrecip
 integer,intent(in) :: natom,ntypat
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: basis(3,3),datagrid(n1*n2*n3),origin(3)
 real(dp),intent(in) :: xcart(3,natom), znucl(ntypat)

!Local variables-------------------------------
!scalars
 integer :: ix,iy,iz,nslice,nsym,iatom
 real(dp) :: fact
!arrays
 real(dp) :: tau(3,natom)

! *************************************************************************

 DBG_ENTER("COLL")

 if (all(realrecip/= [0,1])) then
   MSG_BUG(sjoin('The argument realrecip should be 0 or 1, received:', itoa(realrecip)))
 end if

!conversion between ABINIT default units and XCrysden units
 fact=Bohr_Ang; if (realrecip ==1) fact=one/fact  !since we are in reciprocal space

!TODO insert crystalline structure and dummy atoms in case of reciprocal space
!need to convert basis too

 write(nunit,'(1X,A)')  'DIM-GROUP'
 write(nunit,*) '3  1'
 write(nunit,'(1X,A)') 'PRIMVEC'
 do iy = 1,3
   write(nunit,'(3(ES17.10,2X))') (Bohr_Ang*basis(ix,iy), ix=1,3)
 end do
!
!generate translated coordinates to fit origin shift
!
 do iatom = 1,natom
   tau (:,iatom) = xcart(:,iatom) - origin(:)
 end do

 write(nunit,'(1X,A)') 'PRIMCOORD'
 write(nunit,*) natom, ' 1'
 do iatom = 1,natom
   write(nunit,'(i9,3(3X,ES17.10))') NINT(znucl(typat(iatom))), &  ! WARNING alchemy not supported by XCrysden
     Bohr_Ang*tau(1,iatom), &
     Bohr_Ang*tau(2,iatom), &
     Bohr_Ang*tau(3,iatom)
 end do
 write(nunit,'(1X,A)') 'ATOMS'
 do iatom = 1,natom
   write(nunit,'(i9,3(3X,ES17.10))') NINT(znucl(typat(iatom))), & ! WARNING alchemy not supported by XCrysden
     Bohr_Ang*tau(1,iatom), &
     Bohr_Ang*tau(2,iatom), &
     Bohr_Ang*tau(3,iatom)
 end do

 write(nunit,'(a)')' BEGIN_BLOCK_DATAGRID3D'
 write(nunit,'(a)')' datagrid'
 write(nunit,'(a)')' DATAGRID_3D_DENSITY'
!NOTE: XCrysden uses aperiodical data grid
 write(nunit,*)n1+1,n2+1,n3+1
 write(nunit,*)origin
 write(nunit,*)basis(:,1)*fact
 write(nunit,*)basis(:,2)*fact
 write(nunit,*)basis(:,3)*fact

 nslice=1
 do iz=1,n3
   do iy=1,n2
     write(nunit,'(8es16.8)') datagrid(1+n1*(nslice-1):n1+n1*(nslice-1)),datagrid(1+n1*(nslice-1))
     nslice=nslice+1
   end do
   nsym=nslice-n2
   write (nunit,'(8es16.8)') datagrid(1+n1*(nsym-1):n1+n1*(nsym-1)),datagrid(1+n1*(nsym-1))
 end do

!Now write upper plane
 nslice=1
 do iy=1,n2
   write (nunit,'(8es16.8)') datagrid(1+n1*(nslice-1):n1+n1*(nslice-1)),datagrid(1+n1*(nslice-1))
   nslice=nslice+1
 end do

 nsym=nslice-n2
 write (nunit,'(8es16.8)') datagrid(1+n1*(nsym-1):n1+n1*(nsym-1)),datagrid(1+n1*(nsym-1))

 write (nunit,'(a)')' END_DATAGRID_3D'
 write (nunit,'(a)')' END_BLOCK_DATAGRID3D'

 DBG_EXIT("COLL")

end subroutine printxsf
!!***

!----------------------------------------------------------------------

!!****f* m_pptools/print_fofr_ri
!! NAME
!!  print_fofr_ri
!!
!! FUNCTION
!!  Print the [real,imaginary] part of fofr on unit unit
!!
!! INPUTS
!!  ri_mode =
!!    "RI" if both real and imag part are wanted
!!    "R"  for real part
!!    "I"  for imaginary part
!!  nx,ny,nz,ldx,ldy,ldz = Logical and physical dimensions of the array.
!!  fofr(2,ldx,ldy,ldz) = Input data
!!  [unit] = Fortran unit number. Default: std_out
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_cut3d
!!
!! CHILDREN
!!      wrap2_pmhalf
!!
!! SOURCE


subroutine print_fofr_ri(ri_mode,nx,ny,nz,ldx,ldy,ldz,fofr,unit)


!Arguments -----------------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz
 integer,optional,intent(in) :: unit
 character(len=*),intent(in) :: ri_mode
!arrays
 real(dp),intent(in) :: fofr(2,ldx,ldy,ldz)

!Local variables-------------------------------
!scalars
 integer :: ount,ix,iy,iz
!arrays

! *************************************************************************

 ount = std_out; if (PRESENT(unit)) ount = unit

 SELECT CASE (ri_mode)
 CASE ("RI","ri")
   do iz=1,nz
     do iy=1,ny
       do ix=1,nx
         write(ount,'(2f20.16)') fofr(:,ix,iy,iz)
       end do
     end do
   end do

 CASE ("R","r")
   do iz=1,nz
     do iy=1,ny
       do ix=1,nx
         write(ount,'(f20.16)') fofr(1,ix,iy,iz)
       end do
     end do
   end do

 CASE ("I","i")
   do iz=1,nz
     do iy=1,ny
       do ix=1,nx
         write(ount,'(f20.16)') fofr(2,ix,iy,iz)
       end do
     end do
   end do

 CASE DEFAULT
   MSG_ERROR("Wrong ri_mode")
 END SELECT

end subroutine print_fofr_ri
!!***

!----------------------------------------------------------------------

!!****f* m_pptools/print_fofr_xyzri
!! NAME
!!  print_fofr_xyzri
!!
!! FUNCTION
!!  Print the Cartesian coordinates and the [real,imaginary] part of fofr on unit unit
!!
!! INPUTS
!!  ri_mode =
!!    "RI" if both real and imag part are wanted
!!    "R"  for real part
!!    "I"  for imaginary part
!!  nx,ny,nz,ldx,ldy,ldz = Logical and physical dimensions of the array.
!!  fofr(2,ldx,ldy,ldz) = Input data
!!  rprimd(3,3)=Lattive vectors in Bohr
!!  [conv_fact] = Conversion factor for rprimd (rprimd is multiplied by conv_fact). Default is one
!!  [unit] = Fortran unit number. Default: std_out
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_cut3d
!!
!! CHILDREN
!!      wrap2_pmhalf
!!
!! SOURCE


subroutine print_fofr_xyzri(ri_mode,nx,ny,nz,ldx,ldy,ldz,fofr,rprimd,conv_fact,unit)


!Arguments -----------------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz
 integer,optional,intent(in) :: unit
 real(dp),optional,intent(in) :: conv_fact
 character(len=*),intent(in) :: ri_mode
!arrays
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: fofr(2,ldx,ldy,ldz)

!Local variables-------------------------------
!scalars
 integer :: ount,ix,iy,iz
 real(dp) :: xnow,ynow,znow,my_cfact

! *************************************************************************

 ount = std_out; if (PRESENT(unit)) ount = unit
 my_cfact = one; if (PRESENT(conv_fact)) my_cfact = conv_fact

 SELECT CASE (ri_mode)
 CASE ("RI","ri")
   do iz=1,nz
     do iy=1,ny
       do ix=1,nx
         xnow = rprimd(1,1)*(ix-1)/nx + rprimd(1,2)*(iy-1)/ny + rprimd(1,3)*(iz-1)/nz
         ynow = rprimd(2,1)*(ix-1)/nx + rprimd(2,2)*(iy-1)/ny + rprimd(2,3)*(iz-1)/nz
         znow = rprimd(3,1)*(ix-1)/nx + rprimd(3,2)*(iy-1)/ny + rprimd(3,3)*(iz-1)/nz
         write(ount,'(3f16.10,2f20.16)') my_cfact*xnow, my_cfact*ynow, my_cfact*znow,fofr(:,ix,iy,iz)
       end do
     end do
   end do

 CASE ("R","r")
   do iz=1,nz
     do iy=1,ny
       do ix=1,nx
         xnow = rprimd(1,1)*(ix-1)/nx + rprimd(1,2)*(iy-1)/ny + rprimd(1,3)*(iz-1)/nz
         ynow = rprimd(2,1)*(ix-1)/nx + rprimd(2,2)*(iy-1)/ny + rprimd(2,3)*(iz-1)/nz
         znow = rprimd(3,1)*(ix-1)/nx + rprimd(3,2)*(iy-1)/ny + rprimd(3,3)*(iz-1)/nz
         write(ount,'(3f16.10,f20.16)') my_cfact*xnow, my_cfact*ynow, my_cfact*znow,fofr(1,ix,iy,iz)
       end do
     end do
   end do

 CASE ("I","i")
   do iz=1,nz
     do iy=1,ny
       do ix=1,nx
         xnow = rprimd(1,1)*(ix-1)/nx + rprimd(1,2)*(iy-1)/ny + rprimd(1,3)*(iz-1)/nz
         ynow = rprimd(2,1)*(ix-1)/nx + rprimd(2,2)*(iy-1)/ny + rprimd(2,3)*(iz-1)/nz
         znow = rprimd(3,1)*(ix-1)/nx + rprimd(3,2)*(iy-1)/ny + rprimd(3,3)*(iz-1)/nz
         write(ount,'(3f16.10,f20.16)') my_cfact*xnow, my_cfact*ynow, my_cfact*znow,fofr(2,ix,iy,iz)
       end do
     end do
   end do

 CASE DEFAULT
   MSG_ERROR("Wrong ri_mode")
 END SELECT

end subroutine print_fofr_xyzri
!!***

!----------------------------------------------------------------------

!!****f* m_pptools/print_fofr_cube
!! NAME
!!  print_fofr_cube
!!
!! FUNCTION
!!  Print array fofr in the cube file format
!!
!! INPUTS
!!  nx,ny,nz,ldx,ldy,ldz = Logical and physical dimensions of the array.
!!  fofr(2,ldx,ldy,ldz) = Input data
!!  rprimd(3,3)=Lattive vectors in Bohr
!!  [unit] = Fortran unit number. Default: std_out
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_cut3d
!!
!! CHILDREN
!!      wrap2_pmhalf
!!
!! SOURCE


subroutine print_fofr_cube(nx,ny,nz,ldx,ldy,ldz,fofr,rprimd,natom,znucl_atom,xcart,unit)


!Arguments -----------------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,natom
 integer,optional,intent(in) :: unit
!arrays
 integer,intent(in) :: znucl_atom(natom)
 real(dp),intent(in) :: rprimd(3,3),xcart(3,natom)
 real(dp),intent(in) :: fofr(2,ldx,ldy,ldz)

!Local variables-------------------------------
!scalars
 integer,parameter :: cplx=2
 integer :: ount,ix,iy,iz,iatom
!arrays

! *************************************************************************

 ount = std_out; if (PRESENT(unit)) ount = unit

! EXAMPLE FROM THE WEB
! CPMD CUBE FILE.
! OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z
! 3    0.000000    0.000000    0.000000
! 40    0.283459    0.000000    0.000000
! 40    0.000000    0.283459    0.000000
! 40    0.000000    0.000000    0.283459
! 8    0.000000    5.570575    5.669178    5.593517
! 1    0.000000    5.562867    5.669178    7.428055
! 1    0.000000    7.340606    5.669178    5.111259
! -0.25568E-04  0.59213E-05  0.81068E-05  0.10868E-04  0.11313E-04  0.35999E-05

 write(ount,'(a)') 'ABINIT generated cube file'
 write(ount,'(a)') 'from cut3d tool'

 write(ount,'(i9,3(1x,f12.6))') natom,0.,0.,0.
 write(ount,'(i9,3(1x,f12.6))') nx,(rprimd(iy,1)/nx, iy=1,3)
 write(ount,'(i9,3(1x,f12.6))') ny,(rprimd(iy,2)/ny, iy=1,3)
 write(ount,'(i9,3(1x,f12.6))') nz,(rprimd(iy,3)/nz, iy=1,3)

 do iatom=1,natom
   write(ount,'(i9,4(3X,ES17.10))') znucl_atom(iatom),0.d0,xcart(1:3,iatom)
 end do

! Note C ordering of the indexes
 if (cplx==2) then
   do ix=1,nx
     do iy=1,ny
       do iz=1,nz
         write(ount,'(6(f12.6,2x))') sqrt(fofr(1,ix,iy,iz)**2 + fofr(2,ix,iy,iz)**2 )
       end do
     end do
   end do
 else
   do ix=1,nx
     do iy=1,ny
       do iz=1,nz
         write(ount,'(6(f12.6,2x))') fofr(1,ix,iy,iz)
       end do
     end do
   end do
 end if

end subroutine print_fofr_cube
!!***

!----------------------------------------------------------------------

!!****f* m_pptools/printbxsf
!! NAME
!! printbxsf
!!
!! FUNCTION
!!  Print band structure energies in XCrysDen format.
!!
!! INPUTS
!!  eigen(mband,nkpt,nsppol) = eigenvalues in hartree
!!  ewind = energy window around the fermi level.
!!          if ewind /= 0 ==> a band is considered in the plot of FSurf
!!                            only if it is inside [ ef-ewind, ef+ewind ] for some k point
!!          if ewind == 0 ==> all bands will be keept in the _BXSF file
!!  fermie = Fermi energy (Hartree)
!!  gprimd(3,3) = dimensional primitive translations for reciprocal space (bohr^-1)
!!  kptrlatt(3,3) = reciprocal of lattice vectors for full kpoint grid
!!  mband = maximum number of bands
!!  nsppol = 1 for unpolarized, 2 for spin-polarized
!!  shiftk(3,nshiftk) =shift vector for k point grid
!!  fname = filename for the fortran file
!!  symafm(nsym)=(Anti)ferromagnetic symmetries.
!!  use_afm=.TRUE. if (anti)ferromagnetic symmetries are used.
!!
!! OUTPUT
!!  ierr=Status error.
!!  BXSF file.
!!
!! PARENTS
!!      m_ebands,m_ifc
!!
!! CHILDREN
!!      wrap2_pmhalf
!!
!! SOURCE

subroutine printbxsf(eigen,ewind,fermie,gprimd,kptrlatt,mband,&
& nkptirred,kptirred,nsym,use_afm,symrec,symafm,use_tr,nsppol,shiftk,nshiftk,fname,ierr)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,nkptirred,nshiftk,nsppol,nsym
 integer,intent(out) :: ierr
 real(dp),intent(in) :: ewind,fermie
 logical,intent(in) :: use_afm,use_tr
 character(len=*),intent(in) :: fname
!arrays
 integer,intent(in) :: kptrlatt(3,3),symafm(nsym),symrec(3,3,nsym)
 real(dp),intent(in) :: eigen(mband,nkptirred,nsppol),gprimd(3,3)
 real(dp),intent(in) :: kptirred(3,nkptirred),shiftk(3,nshiftk)

!Local variables-------------------------------
!scalars
 integer,parameter :: enough = 50
 integer :: iband,ik1,ik2,ik3,ikgrid,ikpt,indx
 integer :: isppol,isym,maxband,minband,nk1,nk2,nk3,nkptfull,ubxsf,timrev
 integer :: symkptrank, nsymfm, isymfm
 real(dp) :: ene
 character(len=500) :: msg
 type(krank_t) :: krank
!arrays
 integer,allocatable :: fulltoirred(:),symrecfm(:,:,:)
 real(dp) :: kptgrid(3),gmet(3,3)

! *************************************************************************

 ierr = 0

 ! Error if klatt is no simple orthogonal lattice (in red space)
 ! for generalization to MP grids, need new version of XCrysDen
 if (kptrlatt(1,2)/=0 .or. kptrlatt(1,3)/=0 .or. kptrlatt(2,1)/=0 .or. &
     kptrlatt(2,3)/=0 .or. kptrlatt(3,1)/=0 .or. kptrlatt(3,2)/=0 ) then
   write(msg,'(3a)')&
    'kptrlatt should be diagonal, for the FS calculation ',ch10,&
    'Action: use an orthogonal k-grid for the GS calculation '
   MSG_COMMENT(msg)
   ierr = ierr + 1
 end if

 ! Error if there are not at least 2 kpts in each direction:
 ! kptrank will fail for the intermediate points below
 if (abs(kptrlatt(1,1)) < 2 .or. abs(kptrlatt(2,2)) < 2 .or. abs(kptrlatt(3,3)) < 2) then
   write(msg,'(3a)')&
    'You need at least 2 points in each direction in k space to output BXSF files ',ch10,&
    'Action: use an augmented k-grid for the GS calculation (at least 2x2x2) '
   MSG_COMMENT(msg)
   ierr = ierr + 1
 end if

 if (ANY(ABS(shiftk) > tol10)) then
   write(msg,'(3a)')&
    'Origin of the k-grid should be (0,0,0) for the FS calculation ',ch10,&
    'Action: use a non-shifted k-grid for the GS calculation. Returning '
   MSG_COMMENT(msg)
   ierr = ierr + 1
 end if

 if (ierr /= 0) return

 ! Compute reciprocal space metric.
 gmet = MATMUL(TRANSPOSE(gprimd), gprimd)

 if (use_afm) then
   nsymfm = 0
   do isym = 1, nsym
     if (symafm(isym) == 1) nsymfm = nsymfm+1
   end do
   ABI_MALLOC(symrecfm, (3,3,nsymfm))
   isymfm = 0
   do isym = 1, nsym
     if (symafm(isym) == 1) then
       isymfm = isymfm + 1
       symrecfm(:,:,isymfm) = symrec(:,:,isym)
     end if
   end do
 else
   nsymfm = nsym
   ABI_MALLOC(symrecfm, (3,3,nsymfm))
   symrecfm = symrec
 end if

 ! Xcrysden uses aperiodic data-grid (images are included in the grid)
 nk1 = kptrlatt(1,1); nk2 = kptrlatt(2,2); nk3 = kptrlatt(3,3)
 nkptfull = (nk1+1) * (nk2+1) * (nk3+1)

 ABI_MALLOC(fulltoirred, (nkptfull))
 timrev = 0; if (use_tr) timrev=1

 krank = krank_new(nkptirred, kptirred, nsym=nsymfm, symrec=symrecfm, time_reversal=use_tr)

 ! Xcrysden employs the C-ordering for the Fermi Surface (x-y-z)
 ikgrid=0
 do ik1=0,nk1
   do ik2=0,nk2
     do ik3=0,nk3

       ikgrid = ikgrid+1
       kptgrid(1) = DBLE(ik1)/kptrlatt(1,1)
       kptgrid(2) = DBLE(ik2)/kptrlatt(2,2)
       kptgrid(3) = DBLE(ik3)/kptrlatt(3,3)

       ! Find correspondence between the Xcrysden grid and the IBZ
       symkptrank = krank%get_rank(kptgrid)
       fulltoirred(ikgrid) = krank%invrank(symkptrank)

       if (fulltoirred(ikgrid) < 1) then
         if (ierr <= enough) then
           write(msg,'(a,3es16.8,2a,i0,2a)')&
            'kpt = ',kptgrid,ch10,' with rank ', symkptrank, ch10,&
            'has no symmetric among the k-points used in the GS calculation '
           MSG_WARNING(msg)
         end if
         ierr = ierr + 1
       end if

     end do !ik1
   end do !ik2
 end do !ik3

 call krank%free()

 ABI_CHECK(ierr == 0, "See above warnings")

 if (abs(ewind) < tol12 ) then
   ! Keep all bands.
   minband=1
   maxband=mband
 else
   ! Select a subset of bands.
   minband = mband
   maxband = 0
   ene=abs(ewind)
   do isppol=1,nsppol
     do iband=1,mband
       if(minval(eigen(iband,:,isppol))-fermie < -ene) minband = iband
     end do
     do iband=mband,1,-1
       if (maxval(eigen(iband,:,isppol))-fermie > ene) maxband = iband
     end do
   end do ! isppol

 end if ! abs(energy_window)

 ! Dump results to file
 if (open_file(fname,msg, newunit=ubxsf, status='unknown', action="write", form='formatted') /= 0 ) then
   MSG_WARNING(msg)
   ierr=ierr +1; RETURN
 end if

 ! Write header
 write(ubxsf,*)' BEGIN_INFO'
 write(ubxsf,*)'   #'
 write(ubxsf,*)'   # this is a Band-XCRYSDEN-Structure-File for Visualization of Fermi Surface'
 write(ubxsf,*)'   # generated by the ABINIT package'
 write(ubxsf,*)'   #'
 write(ubxsf,*)'   #  bands between ',minband,' and ',maxband
 write(ubxsf,*)'   #'
 if (nsppol == 2 ) then
   write(ubxsf,*)'   # NOTE: the first band is relative to spin-up electrons,'
   write(ubxsf,*)'   # the second band to spin-down and so on .. '
   write(ubxsf,*)'   #'
 end if
 write(ubxsf,*)'   # Launch as: xcrysden --bxsf '
 write(ubxsf,*)'   #'
 write(ubxsf,'(a,es16.8)')'   Fermi Energy: ',fermie
 write(ubxsf,*)' END_INFO'
 write(ubxsf,*)' '
 write(ubxsf,*)' BEGIN_BLOCK_BANDGRID_3D'
 write(ubxsf,*)' band_energies'
 write(ubxsf,*)' BEGIN_BANDGRID_3D'

 write(ubxsf,*)' ',(maxband-minband+1)*nsppol
 write(ubxsf,*)' ',nk1+1,nk2+1,nk3+1
 write(ubxsf,*)' ',shiftk(:,1)
 ! Angstrom units are used in the BXSF format
 write(ubxsf,*)' ',gprimd(:,1)/Bohr_Ang
 write(ubxsf,*)' ',gprimd(:,2)/Bohr_Ang
 write(ubxsf,*)' ',gprimd(:,3)/Bohr_Ang

 ! print out data for all relevant bands and full kpt grid (redundant, yes)
 ! for each kpt in full zone, find equivalent irred kpt and print eigenval
 indx = 0
 do iband=minband,maxband
   do isppol=1,nsppol
     write(ubxsf,*)' BAND: ',indx+minband
     write(ubxsf,'(7(es16.8))')(eigen(iband,fulltoirred(ikpt),isppol),ikpt=1,nkptfull)
     indx=indx+1
   end do
 end do

 write(ubxsf,*)'  END_BANDGRID_3D'
 write(ubxsf,*)' END_BLOCK_BANDGRID_3D'
 close(ubxsf)

 ABI_FREE(fulltoirred)
 ABI_FREE(symrecfm)

end subroutine printbxsf
!!***

!!****f* m_pptools/printvtk
!! NAME
!! printvtk
!!
!! FUNCTION
!!  Print band structure energies and velocities in VTK format.
!!
!! INPUTS
!!  eigen(mband,nkpt,nsppol) = eigenvalues in hartree
!!  ewind = energy window around the fermi level.
!!          if ewind /= 0 ==> a band is considered in the plot of FSurf
!!                            only if it is inside [ ef-ewind, ef+ewind ] for some k point
!!          if ewind == 0 ==> all bands will be keept in the _BXSF file
!!  fermie = Fermi energy (Hartree)
!!  gprimd(3,3) = dimensional primitive translations for reciprocal space (bohr^-1)
!!  kptrlatt(3,3) = reciprocal of lattice vectors for full kpoint grid
!!  mband = maximum number of bands
!!  nsppol = 1 for unpolarized, 2 for spin-polarized
!!  shiftk(3,nshiftk) =shift vector for k point grid
!!  fname = filename for the fortran file
!!  symafm(nsym)=(Anti)ferromagnetic symmetries.
!!  use_afm=.TRUE. if (anti)ferromagnetic symmetries are used.
!!
!! OUTPUT
!!  ierr=Status error.
!!  BXSF file.
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      wrap2_pmhalf
!!
!! SOURCE

subroutine printvtk(eigen,v_surf,ewind,fermie,gprimd,kptrlatt,mband,&
& nkptirred,kptirred,nsym,use_afm,symrec,symafm,use_tr,nsppol,shiftk,nshiftk,fname,ierr)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,nkptirred,nshiftk,nsppol,nsym
 integer,intent(out) :: ierr
 real(dp),intent(in) :: ewind,fermie
 logical,intent(in) :: use_afm,use_tr
 character(len=*),intent(in) :: fname
!arrays
 integer,intent(in) :: kptrlatt(3,3),symafm(nsym),symrec(3,3,nsym)
 real(dp),intent(in) :: eigen(mband,nkptirred,nsppol),gprimd(3,3)
 real(dp),intent(in) :: v_surf(mband,kptrlatt(1,1)+1,kptrlatt(2,2)+1,kptrlatt(3,3)+1,3,nsppol)
 real(dp),intent(in) :: kptirred(3,nkptirred),shiftk(3,nshiftk)

!Local variables-------------------------------
!scalars
 integer :: iband,ikgrid,ikpt1,indx
 integer :: ikpt,jkpt,kkpt,ikpt_fine, ik1, ik2, ik3
 integer :: isppol,isym,itim,maxband,minband,nk1,nk2,nk3,nkptfull,uvtk,timrev
 real(dp) :: ene,res,ss,timsign
 logical :: found
 character(len=500) :: msg, format_str
!arrays
 integer,allocatable :: fulltoirred(:)
 real(dp) :: kconv(3),kpt(3),kptgrid(3),kptsym(3)

! *************************************************************************

 ierr=0

!Error if klatt is no simple orthogonal lattice (in red space)
!for generalization to MP grids, need new version of XCrysDen

 if (kptrlatt(1,2)/=0 .or. kptrlatt(1,3)/=0 .or. kptrlatt(2,1)/=0 .or. &
     kptrlatt(2,3)/=0 .or. kptrlatt(3,1)/=0 .or. kptrlatt(3,2)/=0 ) then
   write(msg,'(3a)')&
    'kptrlatt should be diagonal, for the FS calculation ',ch10,&
    'action: use an orthogonal k-grid for the GS calculation '
   MSG_COMMENT(msg)
   ierr=ierr+1
 end if

 if (ANY(ABS(shiftk(:,:))>tol10)) then
   write(msg,'(3a)')&
    'Origin of the k-grid should be (0,0,0) for the FS calculation ',ch10,&
    'Action: use a non-shifted k-grid for the GS calculation. Returning '
   MSG_COMMENT(msg)
   ierr=ierr+1
 end if

 if (ierr/=0) RETURN

 ! Xcrysden uses aperiodical data-grid
 nk1 = kptrlatt(1,1)
 nk2 = kptrlatt(2,2)
 nk3 = kptrlatt(3,3)
 nkptfull=(nk1+1)*(nk2+1)*(nk3+1)

 ABI_ALLOCATE(fulltoirred,(nkptfull))
 timrev=0; if (use_tr) timrev=1

 !Xcrysden employs C-ordering for the Fermi Surface.
 ierr = 0
 ikgrid=0
 do ik1=0,nk1
   do ik2=0,nk2
     do ik3=0,nk3

       ikgrid=ikgrid+1
       kptgrid(1)=DBLE(ik1)/kptrlatt(1,1)
       kptgrid(2)=DBLE(ik2)/kptrlatt(2,2)
       kptgrid(3)=DBLE(ik3)/kptrlatt(3,3)
       call wrap2_pmhalf(kptgrid(1),kpt(1),res)
       call wrap2_pmhalf(kptgrid(2),kpt(2),res)
       call wrap2_pmhalf(kptgrid(3),kpt(3),res)

       ! === Find correspondence between the Xcrysden grid and the IBZ ===
       ! If AFM case, use only Ferromagetic symmetries.
       found=.FALSE.
       irred: do ikpt1=1,nkptirred
         do itim=0,timrev
           do isym=1,nsym
             if (use_afm.and.symafm(isym)==-1) CYCLE
             timsign = one-two*itim
             kptsym(:) = timsign*(symrec(:,1,isym)*kptirred(1,ikpt1) + &
                                  symrec(:,2,isym)*kptirred(2,ikpt1) + &
                                  symrec(:,3,isym)*kptirred(3,ikpt1))
             call wrap2_pmhalf(kptsym(1),kconv(1),res)
             call wrap2_pmhalf(kptsym(2),kconv(2),res)
             call wrap2_pmhalf(kptsym(3),kconv(3),res)
             ! is kconv equivalent to kpt?
             ss= (kpt(1)-kconv(1))**2 + (kpt(2)-kconv(2))**2 + (kpt(3)-kconv(3))**2
             if (ss < tol6) then
               found=.TRUE.
               fulltoirred(ikgrid)=ikpt1
               exit irred
             end if

           end do !itim
         end do !isym
       end do irred

       if (.not.found) then
         write(msg,'(a,3es16.8,2a)')&
          ' kpt = ',kpt,ch10,' has no symmetric among the irred k-points used in the GS calculation '
         ierr=ierr+1
         MSG_ERROR(msg)
       end if

     end do !ik1
   end do !ik2
 end do !ik3


 if (ierr/=0) then
   ABI_DEALLOCATE(fulltoirred)
   RETURN
 end if

 if (abs(ewind) < tol12 ) then
   ! Keep all bands.
   minband=1
   maxband=mband
 else
   ! Select a subset of bands.
   minband = mband
   maxband = 0
   ene=abs(ewind)
   do isppol=1,nsppol
     do iband=1,mband
       if(minval(eigen(iband,:,isppol))-fermie < -ene) then
         minband = iband
       end if
     end do
     do iband=mband,1,-1
       if (maxval(eigen(iband,:,isppol))-fermie > ene) then
         maxband = iband
       end if
     end do
   end do ! isppol

 end if ! abs(energy_window)

 ! Dump the results on file ===
 if (open_file(fname,msg,newunit=uvtk,status='unknown',form='formatted') /= 0) then
   ABI_DEALLOCATE(fulltoirred)
   MSG_WARNING(msg)
   ierr=ierr +1; RETURN
 end if

 ! write header
 write(uvtk,"(a)") '# vtk DataFile Version 2.0'
 write(uvtk,"(a)") 'Eigen values for the Fermi surface'
 write(uvtk,"(a)") 'ASCII'
 write(uvtk,*) ''
 write(uvtk,"(a)") 'DATASET STRUCTURED_GRID'
 write(uvtk,"(a,3i6)") 'DIMENSIONS', nk1+1,nk2+1,nk3+1
 write(uvtk,"(a,i6,a)") 'POINTS',nkptfull,' float'

 do ik3 = 0, nk3
   do ik2 = 0, nk2
     do ik1 = 0, nk1
       write(uvtk,'(3es16.8)') dble(ik1)/nk1*gprimd(1,1)+ &
                               dble(ik2)/nk2*gprimd(1,2)+ &
                               dble(ik3)/nk3*gprimd(1,3), &
                               dble(ik1)/nk1*gprimd(2,1)+ &
                               dble(ik2)/nk2*gprimd(2,2)+ &
                               dble(ik3)/nk3*gprimd(2,3), &
                               dble(ik1)/nk1*gprimd(3,1)+ &
                               dble(ik2)/nk2*gprimd(3,2)+ &
                               dble(ik3)/nk3*gprimd(3,3)
     end do
   end do
 end do

!print out data for all relevant bands and full kpt grid (redundant, yes)
!for each kpt in full zone, find equivalent irred kpt and print eigenval
 write(uvtk,*) ''
 write(uvtk,"(a,i6)") 'POINT_DATA',nkptfull
 indx=0
 do iband=minband,maxband
   do isppol=1,nsppol
     if (minband+indx < 10) then
       format_str="(a14,i1,1X,a)"
     else
       format_str="(a14,i2,1X,a)"
     end if
     write(uvtk,format_str) 'SCALARS eigval', minband+indx, 'float 1'
     write(uvtk,"(a)") 'LOOKUP_TABLE default'
     write(uvtk,*) ' '
     do kkpt = nk3/2+1, nk3+nk3/2+1
       do jkpt = nk2/2+1, nk2+nk2/2+1
         do ikpt = nk1/2+1, nk1+nk1/2+1
           ik1 = ikpt
           ik2 = jkpt
           ik3 = kkpt
           if (ikpt > nk1+1) ik1 = ikpt - nk1
           if (jkpt > nk2+1) ik2 = jkpt - nk2
           if (kkpt > nk3+1) ik3 = kkpt - nk3
!          get the index with zyx order
           ikpt_fine = (ik1-1)*(nk2+1)*(nk3+1) + (ik2-1)*(nk3+1) + ik3
           write(uvtk,'(es16.8)') eigen(iband,fulltoirred(ikpt_fine),isppol)
         end do
       end do
     end do
     indx=indx+1
   end do
 end do

 write(uvtk,*) ''
 indx=0
 do iband=minband,maxband
   do isppol=1,nsppol
     if (minband+indx < 10) then
       format_str="(a10,i1,1X,a)"
     else
       format_str="(a10,i2,1X,a)"
     end if
     write(uvtk,format_str) 'SCALARS ve', minband+indx, 'float'
     write(uvtk,"(a)") 'LOOKUP_TABLE default'
     write(uvtk,*) ' '
     do kkpt = nk3/2+1, nk3+nk3/2+1
       do jkpt = nk2/2+1, nk2+nk2/2+1
         do ikpt = nk1/2+1, nk1+nk1/2+1
           ik1 = ikpt
           ik2 = jkpt
           ik3 = kkpt
           if (ikpt > nk1+1) ik1 = ikpt - nk1
           if (jkpt > nk2+1) ik2 = jkpt - nk2
           if (kkpt > nk3+1) ik3 = kkpt - nk3
!          write(uvtk,'(3i6,3es16.8)') ik1,ik2,ik3,v_surf(iband,ik1,ik2,ik3,1,isppol), &
!          &                                                 v_surf(iband,ik1,ik2,ik3,2,isppol), &
!          &                                                 v_surf(iband,ik1,ik2,ik3,3,isppol)
           write(uvtk,'(es16.8)') sqrt(v_surf(iband,ik1,ik2,ik3,1,isppol)*v_surf(iband,ik1,ik2,ik3,1,isppol)+ &
&           v_surf(iband,ik1,ik2,ik3,2,isppol)*v_surf(iband,ik1,ik2,ik3,2,isppol)+ &
&           v_surf(iband,ik1,ik2,ik3,3,isppol)*v_surf(iband,ik1,ik2,ik3,3,isppol))
         end do
       end do
     end do
     write(uvtk,format_str) 'SCALARS vx', minband+indx, 'float'
     write(uvtk,"(a)") 'LOOKUP_TABLE default'
     write(uvtk,*) ' '
     do kkpt = nk3/2+1, nk3+nk3/2+1
       do jkpt = nk2/2+1, nk2+nk2/2+1
         do ikpt = nk1/2+1, nk1+nk1/2+1
           ik1 = ikpt
           ik2 = jkpt
           ik3 = kkpt
           if (ikpt > nk1+1) ik1 = ikpt - nk1
           if (jkpt > nk2+1) ik2 = jkpt - nk2
           if (kkpt > nk3+1) ik3 = kkpt - nk3
           write(uvtk,'(es16.8)') v_surf(iband,ik1,ik2,ik3,1,isppol)
         end do
       end do
     end do
     write(uvtk,format_str) 'SCALARS vy', minband+indx, 'float'
     write(uvtk,"(a)") 'LOOKUP_TABLE default'
     write(uvtk,*) ' '
     do kkpt = nk3/2+1, nk3+nk3/2+1
       do jkpt = nk2/2+1, nk2+nk2/2+1
         do ikpt = nk1/2+1, nk1+nk1/2+1
           ik1 = ikpt
           ik2 = jkpt
           ik3 = kkpt
           if (ikpt > nk1+1) ik1 = ikpt - nk1
           if (jkpt > nk2+1) ik2 = jkpt - nk2
           if (kkpt > nk3+1) ik3 = kkpt - nk3
           write(uvtk,'(es16.8)') v_surf(iband,ik1,ik2,ik3,2,isppol)
         end do
       end do
     end do
     write(uvtk,format_str) 'SCALARS vz', minband+indx, 'float'
     write(uvtk,"(a)") 'LOOKUP_TABLE default'
     write(uvtk,*) ' '
     do kkpt = nk3/2+1, nk3+nk3/2+1
       do jkpt = nk2/2+1, nk2+nk2/2+1
         do ikpt = nk1/2+1, nk1+nk1/2+1
           ik1 = ikpt
           ik2 = jkpt
           ik3 = kkpt
           if (ikpt > nk1+1) ik1 = ikpt - nk1
           if (jkpt > nk2+1) ik2 = jkpt - nk2
           if (kkpt > nk3+1) ik3 = kkpt - nk3
           write(uvtk,'(es16.8)') v_surf(iband,ik1,ik2,ik3,3,isppol)
         end do
       end do
     end do
     indx=indx+1
   end do
 end do

 close (uvtk)
 ABI_DEALLOCATE(fulltoirred)

end subroutine printvtk
!!***

END MODULE m_pptools
