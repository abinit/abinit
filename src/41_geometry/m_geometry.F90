!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_geometry
!! NAME
!!  m_geometry
!!
!! FUNCTION
!!  This module contains basic tools to operate on vectors expressed in reduced coordinates.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2019 ABINIT group (MG, MT, FJ, TRangel, DCA, XG, AHR, DJA, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_geometry

 use defs_basis
 use m_abicore
 use m_errors
 use m_atomdata
 use m_sort

 use m_io_tools,       only : open_file
 use m_numeric_tools,  only : uniformrandom, isinteger, set2unit
 use m_symtk,          only : mati3inv, mati3det, matr3inv, symdet
 use m_hide_lapack,    only : matr3eigval
 use m_pptools,        only : prmat
 use m_numeric_tools,  only : wrap2_pmhalf
 use m_hide_lapack,    only : matrginv

 implicit none

 private

 public :: normv              ! Norm of vector(s) in reduced coordinates either in real or reciprocal space.
 public :: vdotw              ! Scalar product between two reduced vectors either in real or reciprocal space.
 public :: acrossb            ! Cross product of two 3-vectors.
 public :: wigner_seitz       ! Find the grid of points falling inside the Wigner-Seitz cell.
 public :: phdispl_cart2red   ! Calculate the displacement vectors for all branches in reduced coordinates.
 public :: getspinrot         ! Compute the components of the spinor rotation matrix
 public :: spinrot_cmat       ! Construct 2x2 complex matrix representing rotation operator in spin-space.
 public :: rotmat             ! Finds the rotation matrix.
 public :: fixsym             ! Check that iatfix does not break symmetry.
 public :: metric             ! Compute metric matrices.
 public :: mkradim            ! Make rprim and acell from rprimd
 public :: mkrdim             ! Make rprimd from acell from rprim
 public :: chkrprimd          ! Test if {rprim,acell,rprimd} are consistent
 public :: chkdilatmx         ! check if dilatation of unit cell is consistent with initial G-sphere
 public :: xcart2xred         ! From cart coords to reduced
 public :: xred2xcart         ! From reduced coords to cart.
 public :: fred2fcart         ! Convert reduced forces into cartesian forces
 public :: fcart2fred         ! Convert cartesian forces into reduced forces
 public :: bonds_lgth_angles  ! Write GEO file
 public :: randomcellpos      ! Creates unit cell with random atomic positions.
 public :: ioniondist         ! Compute ion-ion distances
 public :: dist2              ! Calculates the distance of v1 and v2 in a crystal by epeating the unit cell
 public :: shellstruct        ! Calculates shell structure (multiplicities, radii)
 public :: remove_inversion   ! Remove the inversion symmetry and improper rotations
 public :: symredcart         ! Convert a symmetry operation from reduced coordinates (integers) to cart coords (reals)
 public :: strainsym          ! Symmetrize the strain tensor.
 public :: stresssym          ! Symmetrize the stress tensor.
 public :: strconv            ! Convert from symmetric storage mode in reduced coords to cart coords.
 public :: littlegroup_pert   ! Determines the set of symmetries that leaves a perturbation invariant.
 public :: irreducible_set_pert  ! Determines a set of perturbations that form a basis

 interface normv
  module procedure normv_rdp_vector
  module procedure normv_int_vector
  !module procedure normv_int_vector_array  ! WARNING for the time being, do not use these 2 procedures,
  !module procedure normv_rdp_vector_array  ! sunstudio12 is not able to resolve which sub should be called.
 end interface normv

 interface vdotw
  module procedure vdotw_rr_vector
  module procedure vdotw_rc_vector
 end interface vdotw

CONTAINS  !===========================================================
!!***

!!****f* m_geometry/normv_rdp_vector
!! NAME
!! normv_rdp_vector
!!
!! FUNCTION
!! Compute the norm of a vector expressed in reduced coordinates using the metric met.
!! The result is multiplied by 2pi in case of a vector in reciprocal space
!! to take into account the correct normalisation of the reciprocal lattice vectors
!!
!! INPUTS
!!  xv(3)=Vector in reduced coordinates
!!  met(3,3)=Metric tensor
!!  space=Character defining whether we are working in real (r|R) or reciprocal space (g|G)
!!
!! OUTPUT
!!  normv_rdp_vector=norm of xv
!!
!! NOTES
!!  The routine is able to deal both with a single vector as well as arrays of vectors.
!!  Versions for integer and real vectors are provided.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function normv_rdp_vector(xv,met,space) result(res)


!Arguments ------------------------------------
!scalars
 real(dp) :: res
 character(len=1),intent(in) :: space
!arrays
 real(dp),intent(in) :: met(3,3),xv(3)

! *************************************************************************

 res =  (xv(1)*met(1,1)*xv(1) + xv(2)*met(2,2)*xv(2) + xv(3)*met(3,3)*xv(3)  &
&  +two*(xv(1)*met(1,2)*xv(2) + xv(1)*met(1,3)*xv(3) + xv(2)*met(2,3)*xv(3)) )

 select case (space)
 case ('r','R')
   res=SQRT(res)
 case ('g','G')
   res=two_pi*SQRT(res)
 case default
   MSG_BUG('Wrong value for space')
 end select

end function normv_rdp_vector
!!***

!----------------------------------------------------------------------

!!****f* m_geometry/normv_int_vector
!! NAME
!!  normv_int_vector
!!
!! FUNCTION
!!  Returns the norm of an integer 3D vector expressed in reduced coordinates.
!!  either in real or reciprocal space. In the later case the factor 2pi has
!!  to be included, due to the conventions used in abinit to define the reciprocal lattice.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

function normv_int_vector(xv,met,space) result(res)


!Arguments ------------------------------------
!scalars
 real(dp) :: res
 character(len=1),intent(in) :: space
!arrays
 real(dp),intent(in) :: met(3,3)
 integer,intent(in) :: xv(3)

! *************************************************************************

 res =  ( xv(1)*met(1,1)*xv(1) + xv(2)*met(2,2)*xv(2) + xv(3)*met(3,3)*xv(3)  &
&  +two*( xv(1)*met(1,2)*xv(2) + xv(1)*met(1,3)*xv(3) + xv(2)*met(2,3)*xv(3)) )

 select case (space)
 case ('r','R')
   res=SQRT(res)
 case ('g','G')
   res=two_pi*SQRT(res)
 case default
   MSG_BUG('Wrong value for space')
 end select

end function normv_int_vector
!!***

!----------------------------------------------------------------------

!!****f* m_geometry/normv_int_vector_array
!! NAME
!!  normv_int_vector_array
!!
!! FUNCTION
!!  Returns the norm of an array of integer 3D vectors expressed in reduced coordinates.
!!  either in real or reciprocal space. In the later case the factor 2pi has
!!  to be included, due to the conventions used in abinit to define the reciprocal lattice.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

function normv_int_vector_array(xv,met,space) result(res)


!Arguments ------------------------------------
!scalars
 character(len=1),intent(in) :: space
!arrays
 real(dp),intent(in) :: met(3,3)
 integer,intent(in) :: xv(:,:)
 !this awful trick is needed to avoid problems with abilint
 real(dp) :: res(SIZE(xv(1,:)))
 !real(dp) :: res(SIZE(xv,DIM=2))

! *************************************************************************

 res(:) = ( xv(1,:)*met(1,1)*xv(1,:) + xv(2,:)*met(2,2)*xv(2,:) + xv(3,:)*met(3,3)*xv(3,:)  &
&     +two*(xv(1,:)*met(1,2)*xv(2,:) + xv(1,:)*met(1,3)*xv(3,:) + xv(2,:)*met(2,3)*xv(3,:)) )

 select case (space)
 case ('r','R')
   res(:)=SQRT(res(:))
 case ('g','G')
   res(:)=two_pi*SQRT(res(:))
 case default
   MSG_BUG('Wrong value for space')
 end select

end function normv_int_vector_array
!!***

!----------------------------------------------------------------------

!!****f* m_geometry/normv_rdp_vector_array
!! NAME
!!  normv_rdp_vector_array
!!
!! FUNCTION
!!  Returns the norm of an array of real 3D vectors expressed in reduced coordinates.
!!  either in real or reciprocal space. In the later case the factor 2pi has
!!  to be included, due to the conventions used in abinit to define the reciprocal lattice.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

function normv_rdp_vector_array(xv,met,space) result(res)


!Arguments ------------------------------------
!scalars
 character(len=1),intent(in) :: space
!arrays
 real(dp),intent(in) :: met(3,3)
 real(dp),intent(in) :: xv(:,:)
 !this awful trick is needed to avoid problems with abilint
 real(dp) :: res(SIZE(xv(1,:)))
 !real(dp) :: res(SIZE(xv,DIM=2))

! *************************************************************************

 res(:) = ( xv(1,:)*met(1,1)*xv(1,:) + xv(2,:)*met(2,2)*xv(2,:) + xv(3,:)*met(3,3)*xv(3,:)  &
&     +two*(xv(1,:)*met(1,2)*xv(2,:) + xv(1,:)*met(1,3)*xv(3,:) + xv(2,:)*met(2,3)*xv(3,:)) )

 select case (space)
 case ('r','R')
   res(:)=SQRT(res(:))
 case ('g','G')
   res(:)=two_pi*SQRT(res)
 case default
   MSG_BUG('Wrong value for space')
 end select

end function normv_rdp_vector_array
!!***

!----------------------------------------------------------------------

!!****f* m_geometry/vdotw_rr_vector
!! NAME
!! vdotw_rr_vector
!!
!! FUNCTION
!! Compute the scalar product between two vectors expressed in reduced coordinates
!! The result is multiplied by (2pi)**2 in case of vectors in reciprocal space
!! to take into account the correct normalisation of the reciprocal lattice vectors
!!
!! INPUTS
!!  xv(3),xw(3)=Vectors in reduced coordinates
!!  met(3,3)=Metric tensor
!!  space=Character defining whether we are working in real (r) or reciprocal space (g)
!!
!! OUTPUT
!!  res=scalar product of xv and xw
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function vdotw_rr_vector(xv,xw,met,space) result(res)


!Arguments ------------------------------------
!scalars
 real(dp) :: res
 character(len=1),intent(in) :: space
!arrays
 real(dp),intent(in) :: met(3,3),xv(3),xw(3)

! *************************************************************************

 res = (  met(1,1)* xv(1)*xw(1)                &
&        +met(2,2)* xv(2)*xw(2)                &
&        +met(3,3)* xv(3)*xw(3)                &
&        +met(1,2)*(xv(1)*xw(2) + xv(2)*xw(1)) &
&        +met(1,3)*(xv(1)*xw(3) + xv(3)*xw(1)) &
&        +met(2,3)*(xv(2)*xw(3) + xv(3)*xw(2)) )

 select case (space)
 case ('r','R')
   return
 case ('g','G')
   res= res * (two_pi**2)
 case default
   MSG_BUG('Wrong value for space')
 end select

end function vdotw_rr_vector
!!***

!----------------------------------------------------------------------

!!****f* m_geometry/vdotw_rc_vector
!! NAME
!! vdotw_rc_vector
!!
!! FUNCTION
!! Compute the scalar product between two vectors expressed in reduced coordinates
!! First vector is real, the second one is complex.
!! The result is multiplied by (2pi)**2 in case of vectors in reciprocal space
!! to take into account the correct normalisation of the reciprocal lattice vectors
!!
!! INPUTS
!!  xv(3),xw(3)=Vectors in reduced coordinates
!!  met(3,3)=Metric tensor
!!  space=Character defining whether we are working in real (r) or reciprocal space (g)
!!
!! OUTPUT
!!  res=complex scalar product of xv and xw
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function vdotw_rc_vector(xv,xw,met,space) result(res)


!Arguments ------------------------------------
!scalars
 complex(dpc) :: res
 character(len=1),intent(in) :: space
!arrays
 real(dp),intent(in) :: met(3,3),xv(3)
 complex(dpc),intent(in) :: xw(3)

! *************************************************************************

 res = (  met(1,1)* xv(1)*xw(1)                &
&        +met(2,2)* xv(2)*xw(2)                &
&        +met(3,3)* xv(3)*xw(3)                &
&        +met(1,2)*(xv(1)*xw(2) + xv(2)*xw(1)) &
&        +met(1,3)*(xv(1)*xw(3) + xv(3)*xw(1)) &
&        +met(2,3)*(xv(2)*xw(3) + xv(3)*xw(2)) )

 select case (space)
 case ('r','R')
   return
 case ('g','G')
   res= res * (two_pi**2)
 case default
   MSG_BUG('Wrong value for space')
 end select

end function vdotw_rc_vector
!!***

!----------------------------------------------------------------------

!!****f* m_geometry/acrossb
!! NAME
!! acrossb
!!
!! FUNCTION
!! Calculates the cross product of two 3-vectors
!!
!! INPUTS
!!   a(3): real(dp) vector
!!   b(3): real(dp) vector
!!
!! OUTPUT
!!   c(3): real(dp) vector = a X b
!!
!! PARENTS
!!      calc_b_matrix,m_abimover,simple_j_dia
!!
!! CHILDREN
!!
!! SOURCE

subroutine acrossb(a,b,c)


!Arguments ---------------------------------------------
!arrays
 real(dp),intent(in) :: a(3),b(3)
 real(dp),intent(out) :: c(3)

! *********************************************************************

 c(1) =  a(2)*b(3) - a(3)*b(2)
 c(2) = -a(1)*b(3) + a(3)*b(1)
 c(3) =  a(1)*b(2) - b(1)*a(2)

end subroutine acrossb
!!***

!!****f* m_geometry/wigner_seitz
!! NAME
!! wigner_seitz
!!
!! FUNCTION
!! Calculates a grid of points that falls inside of (and eventually on the surface of)
!! the Wigner-Seitz supercell centered on the origin of the B lattice with primitive
!! translations nmonkh(1)*a_1+nmonkh(2)*a_2+nmonkh(3)*a_3.
!! Subroutine taken from the Wannier90 code. Modified by MG to fulfil abinit coding rules.
!! API slightly changed wrt the wannier90 version.
!!
!! COPYRIGHT
!! Copyright (C) 2007 Jonathan Yates, Arash Mostofi,
!! Young-Su Lee, Nicola Marzari, Ivo Souza, David Vanderbilt.
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  center(3)=The Wigner-Seitz cell is centered on this point in reduced coordinates.
!!  rmet(3,3)=Real space metric ($\textrm{bohr}^{2}$).
!!  kptrlatt(3)=Values defining the supercell.
!!  prtvol=If different from 0 print out the points falling inside the W-S cell and the correponding weights.
!!  lmax(3)=see Notes below.
!!
!! OUTPUT
!!  npts=number of points falling inside the Wigner-Seitz cell
!!  irvec(3,npts)=Reduced coordinated of the points inside the W-S cell
!!  ndegen(npts)=Weigths associated to each point.
!!
!! SIDE EFFECTS
!!  irvec and ndegen are are allocated with the correct
!!  size inside the routine and returned to the caller.
!!
!! NOTES
!! The Wannier functions live in a supercell of the real space unit cell.
!! This supercell is mp_grid unit cells long in each direction
!! The algorithm loops over grid points r on a unit cell that is 8 times larger than this
!! primitive supercell.
!! One of these points is in the W-S cell if it is closer to center(:)
!! than any of the other points R where R are the translation vectors of the supercell.
!! In the end npts contains the total number of grid points that have been found in the Wigner-Seitz cell
!! The number of lattice vectors R along each direction of the supercell is defined by lmax.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wigner_seitz(center, lmax, kptrlatt, rmet, npts, irvec, ndegen, prtvol)


!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol
 integer,intent(out) :: npts
!arrays
 integer,intent(in) :: kptrlatt(3,3),lmax(3)
 integer,allocatable,intent(out) :: irvec(:,:),ndegen(:)
 real(dp),intent(in) :: center(3),rmet(3,3)

!Local variables-------------------------------
!scalars
 integer :: in1,in2,in3,l1,l2,l3,ii,icount,n1,n2,n3
 integer :: l0,l1_max,l2_max,l3_max,nl,verbose,mm1,mm2,mm3
 real(dp) :: tot,dist_min
 real(dp),parameter :: TOL_DIST=tol7
 character(len=500) :: msg
!arrays
 real(dp) :: diff(3)
 real(dp),allocatable :: dist(:),swap2(:,:),swap1(:)

! *************************************************************************

 verbose = 0; if (PRESENT(prtvol)) verbose = prtvol

 if (kptrlatt(1,2)/=0 .or. kptrlatt(2,1)/=0 .or. &
&    kptrlatt(1,3)/=0 .or. kptrlatt(3,1)/=0 .or. &
&    kptrlatt(2,3)/=0 .or. kptrlatt(3,2)/=0 ) then
   MSG_ERROR('Off-diagonal elements of kptrlatt must be zero')
 end if

 n1 = kptrlatt(1,1)
 n2 = kptrlatt(2,2)
 n3 = kptrlatt(3,3)

 l1_max = lmax(1)
 l2_max = lmax(2)
 l3_max = lmax(3)

 nl=(2*l1_max+1)*(2*l2_max+1)*(2*l3_max+1)
 l0=1+l1_max*(1+(2*l2_max+1)**2+(2*l3_max+1)) ! Index of the origin.
 ABI_MALLOC(dist,(nl))

 ! Allocate with maximum size
 mm1 = 2 * n1 + 1
 mm2 = 2 * n2 + 1
 mm3 = 2 * n3 + 1
 ABI_MALLOC(irvec, (3, mm1*mm2*mm3))
 ABI_MALLOC(ndegen, (mm1*mm2*mm3))

 npts = 0
 do in1=-n1,n1
   do in2=-n2,n2
     do in3=-n3,n3

      ! Loop over the nl points R. R=0 corresponds to l1=l2=l3=1, or icount=l0
      icount = 0
      do l1=-l1_max,l1_max
        do l2=-l2_max,l2_max
          do l3=-l3_max,l3_max
            ! Calculate |r - R -r0|^2.
            diff(1) = in1 - l1 * n1 - center(1)
            diff(2) = in2 - l2 * n2 - center(2)
            diff(3) = in3 - l3 * n3 - center(3)
            icount = icount+1
            dist(icount) = DOT_PRODUCT(diff, MATMUL(rmet, diff))
          end do
        end do
      end do

      dist_min = MINVAL(dist)

      if (ABS(dist(l0) - dist_min) < TOL_DIST) then
        npts = npts + 1
        ndegen (npts) = 0
        do ii=1,nl
          if (ABS(dist(ii) - dist_min) < TOL_DIST) ndegen(npts) = ndegen(npts) + 1
        end do
        irvec(1, npts) = in1
        irvec(2, npts) = in2
        irvec(3, npts) = in3
      end if
     end do !in3
   end do !in2
 end do !in1

 if (verbose>=1) then
   write(msg,'(a,i0)')' lattice points in Wigner-Seitz supercell: ',npts
   call wrtout(std_out,msg,'COLL')
   do ii=1,npts
     write(msg,'(a,3(i3),a,i4)')'  vector ', irvec(:,ii),' degeneracy: ', ndegen(ii)
     call wrtout(std_out,msg,'COLL')
   end do
 end if

 ! Check the "sum rule"
 tot = zero
 do ii=1,npts
   tot = tot + one/ndegen(ii)
 end do
 if (ABS(tot-(n1*n2*n3))>tol8) then
   write(msg,'(a,es16.8,a,i0)')'Something wrong in the generation of WS mesh: tot ',tot,' /= ',n1*n2*n3
   MSG_ERROR(msg)
 end if

 ABI_FREE(dist)

 ! Reallocate the output with correct size.
 ABI_MALLOC(swap1,(npts))
 swap1(:)=ndegen(1:npts)
 ABI_FREE(ndegen)
 ABI_MALLOC(ndegen,(npts))
 ndegen=swap1
 ABI_FREE(swap1)

 ABI_MALLOC(swap2,(3,npts))
 swap2(:,:)=irvec(1:3,1:npts)
 ABI_FREE(irvec)
 ABI_MALLOC(irvec,(3,npts))
 irvec=swap2
 ABI_FREE(swap2)

end subroutine wigner_seitz
!!***

!----------------------------------------------------------------------

!!****f* m_geometry/phdispl_cart2red
!! NAME
!!  phdispl_cart2red
!!
!! FUNCTION
!!  Calculates the displacement vectors for all branches in reduced coordinates.
!!  $ displ_red = displ_cart \cdot gprimd $ for each phonon branch.
!!
!! INPUTS
!!  natom=Number of atoms.
!!  gprimd(3,3)Dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  displ_cart(2,3*natom,3*natom)=Phonon displacement in Cartesian coordinates.
!!
!! OUTPUT
!!  displ_red(2,3*natom,3*natom)=Phonon displacement in reduded coordinates.
!!
!! PARENTS
!!      get_tau_k,m_ddb,m_ifc,mka2f,mkph_linwid,read_gkk
!!
!! CHILDREN
!!
!! SOURCE

subroutine phdispl_cart2red(natom,gprimd,displ_cart,displ_red)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
!arrays
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: displ_cart(2,3*natom,3*natom)
 real(dp),intent(out) :: displ_red(2,3*natom,3*natom)

!Local variables-------------------------
!scalars
 integer :: nbranch,jbranch,iatom,idir,ibranch,kdir,k1

! *************************************************************************

 displ_red = zero

 nbranch=3*natom

 do jbranch=1,nbranch
   !
   do iatom=1,natom
     do idir=1,3
       ibranch=idir+3*(iatom-1)
       do kdir=1,3
         k1 = kdir+3*(iatom-1)
         ! WARNING: could be non-transpose of rprimd matrix : to be checked.
         ! 23 june 2004: rprimd becomes gprimd
         ! could be gprim and then multiply by acell...
         ! Nope, checked and ok with gprimd 24 jun 2004
         displ_red(1,ibranch,jbranch) = displ_red(1,ibranch,jbranch) + &
&         gprimd(kdir,idir) * displ_cart(1,k1,jbranch)

         displ_red(2,ibranch,jbranch) = displ_red(2,ibranch,jbranch) + &
&         gprimd(kdir,idir) * displ_cart(2,k1,jbranch)

       end do !kdir
     end do !idir
   end do !iatom
   !
 end do !jbranch

end subroutine phdispl_cart2red
!!***

!----------------------------------------------------------------------

!!****f* m_geometry/getspinrot
!! NAME
!! getspinrot
!!
!! FUNCTION
!! From the symmetry matrix symrel_conv expressed in the coordinate system rprimd,
!! compute the components of the spinor rotation matrix
!!
!! INPUTS
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!! symrel_conv(3,3)=symmetry operation in real space in terms
!!  of primitive translations rprimd
!!
!! OUTPUT
!! spinrot(4)=components of the spinor rotation matrix :
!!  spinrot(1)=$\cos \phi / 2$
!!  spinrot(2)=$\sin \phi / 2 \times u_x$
!!  spinrot(3)=$\sin \phi / 2 \times u_y$
!!  spinrot(4)=$\sin \phi / 2 \times u_z$
!!  where $\phi$ is the angle of rotation, and
!!  $(u_x,u_y,u_z)$ is the normalized direction of the rotation axis
!!
!! NOTES
!! Only the proper part of the symmetry operation is taken into account:
!! pure rotations, while the inversion part is taken away, if present.
!!
!! The whole collection of symmetry matrices is call symrel(3,3,nsym)
!! symrel1 contains just one of those matrices symrel1(3,3)
!!
!! PARENTS
!!      cg_rotate,m_crystal,wfconv
!!
!! CHILDREN
!!      mati3det,matr3inv
!!
!! SOURCE

subroutine getspinrot(rprimd,spinrot,symrel_conv)


!Arguments ------------------------------------
!arrays
 integer,intent(in) :: symrel_conv(3,3)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(out) :: spinrot(4)

!Local variables-------------------------------
!scalars
 integer :: det,ii
 real(dp) :: cos_phi,norminv,phi,scprod,sin_phi
 !character(len=500) :: message
!arrays
 integer :: identity(3,3),symrel1(3,3)
 real(dp) :: axis(3),coord(3,3),coordinvt(3,3),matr1(3,3),matr2(3,3)
 real(dp) :: rprimd_invt(3,3),vecta(3),vectb(3),vectc(3)

!**************************************************************************

 symrel1(:,:)=symrel_conv(:,:)

!Compute determinant of the matrix
 call mati3det(symrel1,det)

!Produce a rotation from an improper symmetry
 if(det==-1)symrel1(:,:)=-symrel1(:,:)

!Test the possibility of the unit matrix
 identity(:,:)=0
 identity(1,1)=1 ; identity(2,2)=1 ; identity(3,3)=1

 if( sum((symrel1(:,:)-identity(:,:))**2)/=0)then

!  Transform symmetry matrix in the system defined by rprimd
   call matr3inv(rprimd,rprimd_invt)
   do ii=1,3
     coord(:,ii)=rprimd_invt(ii,:)
   end do
   call matr3inv(coord,coordinvt)
   do ii=1,3
     matr1(:,ii)=symrel1(:,1)*coord(1,ii)+&
&     symrel1(:,2)*coord(2,ii)+&
&     symrel1(:,3)*coord(3,ii)
   end do
   do ii=1,3
     matr2(:,ii)=coordinvt(1,:)*matr1(1,ii)+&
&     coordinvt(2,:)*matr1(2,ii)+&
&     coordinvt(3,:)*matr1(3,ii)
   end do

!  Find the eigenvector with unit eigenvalue of the
!  rotation matrix in cartesian coordinate, matr2

   matr1(:,:)=matr2(:,:)
   matr1(1,1)=matr1(1,1)-one
   matr1(2,2)=matr1(2,2)-one
   matr1(3,3)=matr1(3,3)-one

!  Compute the axis of rotation and the cos and sin of rotation angle
   if(matr1(1,1)**2 + matr1(2,1)**2 + matr1(3,1)**2 < tol8 )then
!    The first direction is the axis
     axis(1)=one ; axis(2)=zero ; axis(3)=zero
     cos_phi=matr2(2,2)
     sin_phi=matr2(3,2)
   else if(matr1(1,2)**2 + matr1(2,2)**2 + matr1(3,2)**2 < tol8 )then
!    The second direction is the axis
     axis(1)=zero ; axis(2)=one ; axis(3)=zero
     cos_phi=matr2(3,3)
     sin_phi=matr2(1,3)
   else
!    In this case, try use the first and second vector to build the
!    rotation axis : compute their cross product
     axis(1)=matr1(2,1)*matr1(3,2)-matr1(2,2)*matr1(3,1)
     axis(2)=matr1(3,1)*matr1(1,2)-matr1(3,2)*matr1(1,1)
     axis(3)=matr1(1,1)*matr1(2,2)-matr1(1,2)*matr1(2,1)
!    Then, try to normalize it
     scprod=sum(axis(:)**2)
     if(scprod<tol8)then
!      The first and second vectors were linearly dependent
!      Thus, use the first and third vectors
       axis(1)=matr1(2,1)*matr1(3,3)-matr1(2,3)*matr1(3,1)
       axis(2)=matr1(3,1)*matr1(1,3)-matr1(3,3)*matr1(1,1)
       axis(3)=matr1(1,1)*matr1(2,3)-matr1(1,3)*matr1(2,1)
!      Normalize the vector
       scprod=sum(axis(:)**2)
       if(scprod<tol8)then
         MSG_BUG('Cannot find the rotation axis.')
       end if
     end if
     norminv=one/sqrt(scprod)
     axis(:)=axis(:)*norminv

!    Project the axis vector out of the first unit vector,
!    and renormalize the projected vector
!    (the first vector cannot be the axis, as tested before)
     vecta(1)=one-axis(1)**2
     vecta(2)=-axis(1)*axis(2)
     vecta(3)=-axis(1)*axis(3)
     scprod=sum(vecta(:)**2)
     norminv=one/sqrt(scprod)
     vecta(:)=vecta(:)*norminv
!    Rotate the vector A, to get vector B
     vectb(:)=matr2(:,1)*vecta(1)+matr2(:,2)*vecta(2)+matr2(:,3)*vecta(3)
!    Get dot product of vectors A and B, giving cos of the rotation angle
     cos_phi=sum(vecta(:)*vectb(:))
!    Compute the cross product of the axis and vector A
     vectc(1)=axis(2)*vecta(3)-axis(3)*vecta(2)
     vectc(2)=axis(3)*vecta(1)-axis(1)*vecta(3)
     vectc(3)=axis(1)*vecta(2)-axis(2)*vecta(1)
!    Get dot product of vectors B and C, giving sin of the rotation angle
     sin_phi=sum(vectb(:)*vectc(:))
   end if

!  Get the rotation angle, then the parameters of the spinor rotation
!  Here, treat possible inaccurate values of the cosine of phi
   if(cos_phi>  one-tol8 )cos_phi=  one-tol8
   if(cos_phi<-(one-tol8))cos_phi=-(one-tol8)
   phi=acos(cos_phi)
   if(sin_phi<zero)phi=-phi
!  Rectify the angle, such that its absolute values corresponds to
!  180, 120, 90, 60, or 0 degrees
   phi=(nint(six*phi/pi))/six*pi
!  Compute components of the spinor matrix
   spinrot(1)=cos(half*phi)
   spinrot(2)=axis(1)*sin(half*phi)
   spinrot(3)=axis(2)*sin(half*phi)
   spinrot(4)=axis(3)*sin(half*phi)

 else

!  Here, the case of the unit matrix
   axis(:)=zero
   phi=zero
   spinrot(1)=one
   spinrot(2)=zero
   spinrot(3)=zero
   spinrot(4)=zero

 end if ! the case of the identity matrix

!DEBUG
!write(std_out,*)' getspinrot :'
!write(std_out,*)' symrel_conv =',symrel_conv(:,:)
!write(std_out,*)' symrel =',symrel1(:,:)
!write(std_out,*)' rprimd =',rprimd(:,:)
!write(std_out,*)' matr2 =',matr2(:,:)
!write(std_out,*)' matr1 =',matr1(:,:)
!write(std_out,*)' phi (degree)=',phi*180._dp/pi
!write(std_out,'(a,3d16.6)' )' axis=',axis(:)
!write(std_out,*)' vecta=',vecta(:)
!stop
!ENDDEBUG

end subroutine getspinrot
!!***

!!****f* m_geometry/spinrot_cmat
!! NAME
!!  spinrot_cmat
!!
!! FUNCTION
!!  Construct 2x2 complex matrix representing rotation operator in spin-space.
!!
!! INPUTS
!!  spinrot(4)=components of the spinor rotation matrix computed by getspinrot
!!
!! OUTPUT
!!  spinrot(2,2)=Rotation matrix (complex array)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure function spinrot_cmat(spinrot)


!Arguments ------------------------------------
 real(dp),intent(in) :: spinrot(4)
 complex(dpc) :: spinrot_cmat(2,2)

! *************************************************************************

 ! Rotation in spinor space (same equations as in wfconv)
 spinrot_cmat(1,1) = spinrot(1) + j_dpc*spinrot(4)
 spinrot_cmat(1,2) = spinrot(3) + j_dpc*spinrot(2)
 spinrot_cmat(2,1) =-spinrot(3) + j_dpc*spinrot(2)
 spinrot_cmat(2,2) = spinrot(1) - j_dpc*spinrot(4)

 ! My equation
 !spinrot_cmat(1,1) = spinrot(1) - j_dpc*spinrot(4)
 !spinrot_cmat(1,2) =-spinrot(3) - j_dpc*spinrot(2)
 !spinrot_cmat(2,1) = spinrot(3) - j_dpc*spinrot(2)
 !spinrot_cmat(2,2) = spinrot(1) + j_dpc*spinrot(4)

end function spinrot_cmat
!!***

!----------------------------------------------------------------------

!!****f* m_geometry/rotmat
!! NAME
!! rotmat
!!
!! FUNCTION
!! Finds the rotation matrix.
!!
!! INPUTS
!!  xaxis(3)= vectors defining the x axis
!!  zaxis(3)= vectors defining the z axis

!! OUTPUT
!!  inversion_flag = flag that indicates that an inversion operation
!!   on the coordinate system should be done
!!  umat(3,3)= matrix that rotates the x=(1 0 0) and z=(0 0 1) to the new
!!   values defined in xaxis and zaxis
!!
!! NOTES
!! Here I set that the axe x is originally at the 1 0 0 direction and z is originally 0 0 1.
!! So calling rotmat(x',z') will find the rotation
!! matrix for the case in which we rotate the x and z
!! axes from their default values to x' and z'.
!!
!! PARENTS
!!      mlwfovlp_ylmfac,mlwfovlp_ylmfar
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine rotmat(xaxis,zaxis,inversion_flag,umat)


!Arguments ------------------------------------
!scalars
 integer,intent(out) :: inversion_flag
!arrays
 real(dp),intent(in) :: xaxis(3),zaxis(3)
 real(dp),intent(out) :: umat(3,3)

!Local variables-------------------------------
!scalars
 real(dp) :: cosine,xmod,zmod
 character(len=500) :: message
!arrays
 real(dp) :: yaxis(3)

! *************************************************************************

 xmod = sqrt(xaxis(1)**2 + xaxis(2)**2 + xaxis(3)**2)
 zmod = sqrt(zaxis(1)**2 + zaxis(2)**2 + zaxis(3)**2)

 if(xmod < 1.d-8)then
   write(message,'(a,a,a,i6)')&
&   'The module of the xaxis should be greater than 1.d-8,',ch10,&
&   'however, |xaxis|=',xmod
   MSG_BUG(message)
 end if

 if(zmod < 1.d-8)then
   write(message,'(a,a,a,i6)')&
&   'The module of the zaxis should be greater than 1.d-8,',ch10,&
&   'however, |zaxis|=',zmod
   MSG_ERROR(message)
 end if

!verify that both axis are perpendicular
 cosine = (xaxis(1)*zaxis(1) + xaxis(2)*zaxis(2) &
& + xaxis(3)*zaxis(3))/(xmod*zmod)

 if(abs(cosine) > 1.d-8)then
   write(message,'(a,a,a,i6)')&
&   'xaxis and zaxis should be perpendicular,',ch10,&
&   'however, cosine=',cosine
   MSG_BUG(message)
 end if

!new y axis as cross product
 yaxis(1) = (zaxis(2)*xaxis(3) - xaxis(2)*zaxis(3))/(xmod*zmod)
 yaxis(2) = (zaxis(3)*xaxis(1) - xaxis(3)*zaxis(1))/(xmod*zmod)
 yaxis(3) = (zaxis(1)*xaxis(2) - xaxis(1)*zaxis(2))/(xmod*zmod)

!hack to allow inversion operation on coordinate transformation
!uses unlikely large but legal values of proj_x and/or proj_z
!to flag inversion
 inversion_flag=0
 if(xmod>10._dp .or. zmod>10._dp) then
   inversion_flag=1
   write(message, '(4a)' )&
&   'inversion operation will be appended to axis transformation',ch10,&
&   'Action: If you did not intend this, make |z|<10 and |x|<10 ',ch10
   call wrtout(std_out,message,'COLL')
 end if

 umat(1,:) = xaxis(:)/xmod
 umat(2,:) = yaxis(:)
 umat(3,:) = zaxis(:)/zmod

end subroutine rotmat
!!***

!!****f* m_geometry/fixsym
!! NAME
!! fixsym
!!
!! FUNCTION
!! Using input indsym which tells which atoms are related by symmetry,
!! check that iatfix consistently fixes (freezes) all atoms which are
!! related by symmetry, i.e. that iatfix does not break symmetry.
!!
!! INPUTS
!! iatfix(3,natom)=integer array with 1 in every position for which
!!  the atom is to be kept fixed
!!  NOTE that this is not the input data structure for iatfix but it is
!!  the internal data structure used through most of the subroutines
!! indsym(4,nsym,natom)=indirect indexing array for symmetrically related
!!  atoms; 4th element is label of symmetrically related atom
!! natom=number of atoms
!! nsym=number of symmetries (should be > 1 when this is called)
!!
!! OUTPUT
!!  (only checking)
!!
!! NOTE
!!  Stops execution with an error message if iatfix breaks symmetry.
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!
!! SOURCE

subroutine fixsym(iatfix,indsym,natom,nsym)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsym
!arrays
 integer,intent(in) :: iatfix(3,natom),indsym(4,nsym,natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,isym,jatom
 character(len=500) :: message

! *************************************************************************

 if (nsym > 1) then
   do iatom=1,natom
     do isym=1,nsym
!      jatom is the label of a symmetrically related atom
       jatom=indsym(4,isym,iatom)
!      Thus the atoms jatom and iatom must be fixed along the same directions
       if ( iatfix(1,jatom) /=  iatfix(1,iatom) .or. &
&       iatfix(2,jatom) /=  iatfix(2,iatom) .or. &
&       iatfix(3,jatom) /=  iatfix(3,iatom)       ) then
         write(message, '(a,i0,a,i0,7a)' )&
&         'Atom number ',jatom,' is symmetrically  equivalent to atom number ',iatom,',',ch10,&
&         'but according to iatfix, iatfixx, iatfixy and iatfixz, they',ch10,&
&         'are not fixed along the same directions, which is forbidden.',ch10,&
&         'Action: modify either the symmetry or iatfix(x,y,z) and resubmit.'
         MSG_ERROR(message)
       end if
     end do
   end do
 end if

end subroutine fixsym
!!***

!!****f* m_geometry/metric
!! NAME
!! metric
!!
!! FUNCTION
!! Compute first dimensional primitive translation vectors in reciprocal space
!! gprimd from rprimd, and eventually writes out.
!! Then, computes metrics for real and recip space rmet and gmet using length
!! dimensional primitive translation vectors in columns of rprimd(3,3) and gprimd(3,3).
!!  gprimd is the inverse transpose of rprimd.
!!  i.e. $ rmet_{i,j}= \sum_k ( rprimd_{k,i}*rprimd_{k,j} )  $
!!       $ gmet_{i,j}= \sum_k ( gprimd_{k,i}*gprimd_{k,j} )  $
!! Also computes unit cell volume ucvol in $\textrm{bohr}^3$
!!
!! INPUTS
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  iout=unit number of output file.  If iout<0, do not write output.
!!
!! OUTPUT
!!  gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{-2}$).
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  rmet(3,3)=real space metric ($\textrm{bohr}^{2}$).
!!  ucvol=unit cell volume ($\textrm{bohr}^{3}$).
!!
!! PARENTS
!!      afterscfloop,bethe_salpeter,chkinp,clnup1,conducti_nc,conducti_paw
!!      conducti_paw_core,cut3d,d2frnl,dfpt_eltfrhar,dfpt_eltfrkin,dfpt_eltfrxc
!!      dfpt_looppert,dfpt_newvtr,dfpt_scfcv,dist2,elpolariz,emispec_paw,energy
!!      extrapwf,fftprof,finddistrproc,forces,forstr,get_npert_rbz,getkgrid
!!      hartre,hartrestr,ingeo,initaim,initberry,inkpts,inqpt,invacuum,invars2m
!!      ks_ddiago,linear_optics_paw,m_ab7_symmetry,m_crystal,m_cut3d,m_ddb
!!      m_dens,m_effective_potential,m_effective_potential_file,m_fft
!!      m_fft_prof,m_fit_data,m_hamiltonian,m_io_kss,m_ioarr,m_mep,m_pawpwij
!!      m_screening,m_tdep_latt,m_use_ga,m_vcoul,m_wfk,mag_constr,mag_constr_e
!!      memory_eval,mkcore_wvl,mlwfovlp_qp,moddiel,mpi_setup,mrgscr,newrho
!!      newvtr,nres2vres,odamix,optic,pawgrnl,prcref,prcref_PMA,pred_bfgs
!!      pred_delocint,pred_isothermal,pred_langevin,pred_lbfgs,pred_nose
!!      pred_srkna14,pred_verlet,prt_cif,prtefield,prtimg,psolver_rhohxc
!!      rhotoxc,scfcv,screening,setup1,setup_bse,setup_screening,setup_sigma
!!      sigma,smallprim,stress,strhar,symmetrize_rprimd,testkgrid,thmeig
!!      vdw_dftd2,vdw_dftd3,wrt_moldyn_netcdf,wvl_initro,xchybrid_ncpp_cc
!!      xfpack_vin2x,xfpack_x2vin
!!
!! CHILDREN
!!      matr3inv,wrtout
!!
!! SOURCE

subroutine metric(gmet,gprimd,iout,rmet,rprimd,ucvol)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout
 real(dp),intent(out) :: ucvol
!arrays
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(out) :: gmet(3,3),gprimd(3,3),rmet(3,3)

!Local variables-------------------------------
!scalars
 integer :: nu
 character(len=500) :: message
!arrays
 real(dp) :: angle(3)

! *************************************************************************

!Compute unit cell volume
 ucvol=rprimd(1,1)*(rprimd(2,2)*rprimd(3,3)-rprimd(3,2)*rprimd(2,3))+&
& rprimd(2,1)*(rprimd(3,2)*rprimd(1,3)-rprimd(1,2)*rprimd(3,3))+&
& rprimd(3,1)*(rprimd(1,2)*rprimd(2,3)-rprimd(2,2)*rprimd(1,3))

!Check that the input primitive translations are not linearly dependent (and none is zero); i.e. ucvol~=0
!Also ask that the mixed product is positive.
 if (abs(ucvol)<tol12) then
!  write(std_out,*)"rprimd",rprimd,"ucvol",ucvol
   write(message,'(5a)')&
&   'Input rprim and acell gives vanishing unit cell volume.',ch10,&
&   'This indicates linear dependency between primitive lattice vectors',ch10,&
&   'Action: correct either rprim or acell in input file.'
   MSG_ERROR(message)
 end if
 if (ucvol<zero)then
   write(message,'(2a,3(a,3es16.6,a),7a)')&
&   'Current rprimd gives negative (R1xR2).R3 . ',ch10,&
&   'Rprimd =',rprimd(:,1),ch10,&
&   '        ',rprimd(:,2),ch10,&
&   '        ',rprimd(:,3),ch10,&
&   'Action: if the cell size and shape are fixed (optcell==0),',ch10,&
&   '        exchange two of the input rprim vectors;',ch10,&
&   '        if you are optimizing the cell size and shape (optcell/=0),',ch10,&
&   '        maybe the move was too large, and you might try to decrease strprecon.'
   MSG_ERROR(message)
 end if

!Generates gprimd
 call matr3inv(rprimd,gprimd)

!Write out rprimd, gprimd and ucvol
 if (iout>=0) then
   write(message,'(2a)')' Real(R)+Recip(G) ','space primitive vectors, cartesian coordinates (Bohr,Bohr^-1):'
   call wrtout(iout,message,'COLL')
   do nu=1,3
     write(message, '(1x,a,i1,a,3f11.7,2x,a,i1,a,3f11.7)' ) &
&     'R(',nu,')=',rprimd(:,nu)+tol10,&
&     'G(',nu,')=',gprimd(:,nu)+tol10
     call wrtout(iout,message,'COLL')
   end do
   write(message,'(a,1p,e15.7,a)') ' Unit cell volume ucvol=',ucvol+tol10,' bohr^3'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

!Compute real space metric.
 rmet = MATMUL(TRANSPOSE(rprimd),rprimd)

!Compute reciprocal space metric.
 gmet = MATMUL(TRANSPOSE(gprimd),gprimd)

!Write out the angles
 if (iout>=0) then
   angle(1)=acos(rmet(2,3)/sqrt(rmet(2,2)*rmet(3,3)))/two_pi*360.0d0
   angle(2)=acos(rmet(1,3)/sqrt(rmet(1,1)*rmet(3,3)))/two_pi*360.0d0
   angle(3)=acos(rmet(1,2)/sqrt(rmet(1,1)*rmet(2,2)))/two_pi*360.0d0
   write(message, '(a,3es16.8,a)' )' Angles (23,13,12)=',angle(1:3),' degrees'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

end subroutine metric
!!***

!!****f* m_geometry/mkradim
!! NAME
!! mkradim
!!
!! FUNCTION
!!  Not so trivial subroutine to make dimensionless real space
!!  primitive translations rprim(3,3) from dimensional rprimd(3).
!!  also make acell(3).
!!
!! INPUTS
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!              where: rprimd(i,j)=rprim(i,j)*acell(j)
!!
!! OUTPUT
!!  acell(3)=unit cell length scales (bohr)
!!  rprim(3,3)=dimensionless real space primitive translations
!!
!! PARENTS
!!      gstate,gstateimg,ingeo,m_ddk,m_pimd,m_use_ga,pred_steepdesc
!!      predict_pimd,wvl_memory,xfpack_vin2x
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkradim(acell,rprim,rprimd)


!Arguments ------------------------------------
!arrays
 real(dp),intent(out) :: acell(3),rprim(3,3)
 real(dp),intent(in) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii

! *************************************************************************

!Use a representation based on normalised rprim vectors
 do ii=1,3
   acell(ii)=sqrt(rprimd(1,ii)**2+rprimd(2,ii)**2+rprimd(3,ii)**2)
   rprim(:,ii)=rprimd(:,ii)/acell(ii)
 end do

end subroutine mkradim
!!***

!!****f* m_geometry/chkrprimd
!!
!! NAME
!! chkrprimd
!!
!! FUNCTION
!! Test if {rprim,acell,rprimd} are consistent
!! It means that rprimd can be reconstructed from the rprim and acell
!! Output a message if is not the case
!!
!! INPUTS
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine chkrprimd(acell,rprim,rprimd,iout)


!Arguments ------------------------------------
!scalars
integer,intent(in) :: iout
!arrays
real(dp),intent(in) :: rprim(3,3)
real(dp),intent(in) :: rprimd(3,3)
real(dp),intent(in) :: acell(3)

!Local variables-------------------------------
!scalars
integer :: ii,jj
!arrays
real(dp) :: rprimd_test(3,3)
logical :: equal

! ***********************************************************

!###########################################################
!### 1. Compute rprimd from rprim and acell
 do ii=1,3
   do jj=1,3
     rprimd_test(ii,jj)=rprim(ii,jj)*acell(jj)
   end do
 end do


!###########################################################
!### 2. Compare rprimd and rprimd_test

 equal=.TRUE.
 do ii=1,3
   do jj=1,3
     if (abs(rprimd_test(ii,jj)-rprimd(ii,jj))>1.E-12) then
       equal=.FALSE.
     end if
   end do
 end do

 if (equal)then
   write(iout,*) 'chkrprimd: rprimd is consistent'
 else
   write(iout,*) 'chkrprimd: rprimd is NOT consistent ERROR'
 end if

end subroutine chkrprimd
!!***

!!****f* m_geometry/chkdilatmx
!! NAME
!! chkdilatmx
!!
!! FUNCTION
!! Check whether the new rprimd does not give a too large number
!! of plane waves, compared to the one booked for rprimd, taking
!! into account the maximal dilatation dilatmx. Actually check whether
!! the new Fermi sphere is inside the old one, dilated.
!!
!! INPUTS
!!  chkdilatmx_ = if 1, will prevent to have any vector outside the Fermi sphere, possibly
!!       by rescaling (three times at most), and then stopping the execution
!!                if 0, simply send a warning, but continues execution
!!  dilatmx     = maximal dilatation factor (usually the input variable)
!!  rprimd      = new primitive vectors
!!  rprimd_orig = original primitive vectors (usually the input variable)
!!
!! OUTPUT
!!  dilatmx_errmsg=Emptry string if calculation can continue.
!!            If the calculation cannot continue, dilatmx_errmsg will contain
!!            the message that should be reported in the output file.
!!
!!            Client code should handle a possible problem with the following test:
!!
!!              if (LEN_TRIM(dilatmx_errmsg) then
!!                dump dilatmx_errmsg to the main output file.
!!                handle_error
!!              end if
!!
!!
!! PARENTS
!!      driver,mover
!!
!! CHILDREN
!!      matr3eigval,matr3inv
!!
!! SOURCE

subroutine chkdilatmx(chkdilatmx_,dilatmx,rprimd,rprimd_orig,dilatmx_errmsg)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: chkdilatmx_
 real(dp),intent(in) :: dilatmx
 character(len=500),intent(out) :: dilatmx_errmsg
!arrays
 real(dp),intent(inout) :: rprimd(3,3)
 real(dp),intent(in) :: rprimd_orig(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,mu
 real(dp) :: alpha,dilatmx_new
!arrays
 real(dp) :: eigval(3),gprimd_orig(3,3),met(3,3),old_to_new(3,3)
 character(len=500) :: message

! *************************************************************************

!Generates gprimd
 call matr3inv(rprimd_orig,gprimd_orig)

!Find the matrix that transform an original xcart to xred, then to the new xcart
 do mu=1,3
   old_to_new(mu,:)=rprimd(mu,1)*gprimd_orig(:,1)+&
&   rprimd(mu,2)*gprimd_orig(:,2)+&
&   rprimd(mu,3)*gprimd_orig(:,3)
 end do

!The largest increase in length will be obtained thanks
!to the diagonalization of the corresponding metric matrix :
!it is the square root of its largest eigenvalue.
 do ii=1,3
   do jj=1,3
     met(ii,jj)=old_to_new(1,ii)*old_to_new(1,jj)+&
&     old_to_new(2,ii)*old_to_new(2,jj)+&
&     old_to_new(3,ii)*old_to_new(3,jj)
   end do
 end do

 call matr3eigval(eigval,met)

 dilatmx_new=sqrt(maxval(eigval(:)))

 dilatmx_errmsg = ""
 if(dilatmx_new>dilatmx+tol6)then

! MJV 2014 07 22: correct rprim to maximum jump allowed by dilatmx
! XG 20171011 : eigenvalues of "old_to_old" tensor are of course the unity !

   if(chkdilatmx_/=0)then
     alpha = (dilatmx - one) / (dilatmx_new - one)
!    for safety, only 90 percent of max jump
     alpha = 0.9_dp * alpha

     rprimd = alpha * rprimd + (one - alpha) * rprimd_orig

     write(dilatmx_errmsg,'(3a,es16.6,4a,es16.6,2a,es16.6,a)')&
&     'The new primitive vectors rprimd (an evolving quantity)',ch10,&
&     'are too large with respect to the old rprimd and the accompanying dilatmx:',dilatmx,ch10,&
&     'This large change of unit cell parameters is not allowed by the present value of dilatmx.',ch10,&
&     'An adequate value would have been dilatmx_new=',dilatmx_new,ch10,&
&     'Calculation continues with limited jump, by rescaling the projected move by the factor',alpha,'.'
   else
     write(message, '(3a,es16.6,2a,es16.6,2a)' )&
&     'The new primitive vectors rprimd (an evolving quantity)',ch10,&
&     'are too large, given the initial rprimd and the accompanying dilatmx:',dilatmx,ch10,&
&     'An adequate value would have been dilatmx_new=',dilatmx_new,ch10,&
&     'As chkdilatmx=0, assume experienced user. Execution will continue.'
     MSG_WARNING(message)
   end if

 end if

end subroutine chkdilatmx
!!***

!!****f* m_geometry/mkrdim
!! NAME
!! mkrdim
!!
!! FUNCTION
!!  Trivial subroutine to make dimensional real space
!!  primitive translations from length scales acell(3)
!!  and dimensionless translations rprim(3,3).
!!
!! INPUTS
!!  acell(3)=unit cell length scales (bohr)
!!  rprim(3,3)=dimensionless real space primitive translations
!!
!! OUTPUT
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!              where: rprimd(i,j)=rprim(i,j)*acell(j)
!!
!! PARENTS
!!      bethe_salpeter,dfpt_looppert,dfpt_symph,driver,finddistrproc
!!      get_npert_rbz,gstateimg,harmonic_thermo,ingeo,invars1,invars2m,m_ddb
!!      m_ifc,m_results_img,m_use_ga,memory_eval,mpi_setup,outvar_o_z,pred_bfgs
!!      pred_isothermal,pred_lbfgs,pred_steepdesc,pred_verlet,predict_pimd
!!      randomcellpos,screening,setup1,setup_bse,setup_screening,setup_sigma
!!      sigma,thmeig,wvl_setboxgeometry,xfpack_x2vin
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkrdim(acell,rprim,rprimd)

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: acell(3),rprim(3,3)
 real(dp),intent(out) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,jj

! *************************************************************************

 do ii=1,3
   do jj=1,3
     rprimd(ii,jj)=rprim(ii,jj)*acell(jj)
   end do
 end do

end subroutine mkrdim
!!***

!!****f* m_geometry/xcart2xred
!! NAME
!! xcart2xred
!!
!! FUNCTION
!! Convert from cartesian coordinates xcart(3,natom) in bohr to
!! dimensionless reduced coordinates xred(3,natom) by using
!! xred(mu,ia)=gprimd(1,mu)*xcart(1,ia)
!!            +gprimd(2,mu)*xcart(2,ia)
!!            +gprimd(3,mu)*xcart(3,ia)
!! where gprimd is the inverse of rprimd
!! Note that the reverse operation is deon by xred2xcart
!!
!! INPUTS
!!  natom=number of atoms in unit cell
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  xcart(3,natom)=cartesian coordinates of atoms (bohr)
!!
!! OUTPUT
!!  xred(3,natom)=dimensionless reduced coordinates of atoms
!!
!! PARENTS
!!      driver,evdw_wannier,ingeo,m_cut3d,m_dens,m_effective_potential
!!      m_effective_potential_file,m_mep,m_paw_pwaves_lmn,m_pred_lotf
!!      mkcore_paw,mkcore_wvl,mover_effpot,pawmkaewf,pimd_langevin_npt
!!      pimd_langevin_nvt,pimd_nosehoover_npt,pimd_nosehoover_nvt,prcref
!!      prcref_PMA,pred_delocint,pred_diisrelax,pred_isokinetic,pred_isothermal
!!      pred_langevin,pred_moldyn,pred_nose,pred_srkna14,pred_steepdesc
!!      pred_velverlet,pred_verlet,relaxpol,wrt_moldyn_netcdf
!!      wvl_setboxgeometry
!!
!! CHILDREN
!!      matr3inv
!!
!! SOURCE

subroutine xcart2xred(natom,rprimd,xcart,xred)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
!arrays
 real(dp),intent(in) :: rprimd(3,3),xcart(3,natom)
 real(dp),intent(out) :: xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,mu
!arrays
 real(dp) :: gprimd(3,3)

! *************************************************************************

 call matr3inv(rprimd,gprimd)
 do iatom=1,natom
   do mu=1,3
     xred(mu,iatom)= gprimd(1,mu)*xcart(1,iatom)+gprimd(2,mu)*xcart(2,iatom)+gprimd(3,mu)*xcart(3,iatom)
   end do
 end do

end subroutine xcart2xred
!!***

!!****f* m_geometry/xred2xcart
!! NAME
!! xred2xcart
!!
!! FUNCTION
!! Convert from dimensionless reduced coordinates xred(3,natom)
!! to cartesian coordinates xcart(3,natom) in bohr by using
!! xcart(mu,ia)=rprimd(mu,1)*xred(1,ia)
!!             +rprimd(mu,2)*xred(2,ia)
!!             +rprimd(mu,3)*xred(3,ia)
!! Note that the reverse operation is done by xcart2xred.F90
!!
!! INPUTS
!!  natom=number of atoms in unit cell
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  xred(3,natom)=dimensionless reduced coordinates of atoms
!!
!! OUTPUT
!!  xcart(3,natom)=cartesian coordinates of atoms (bohr)
!!
!! PARENTS
!!      afterscfloop,berryphase,berryphase_new,bonds_lgth_angles,constrf,cut3d
!!      denfgr,driver,evdw_wannier,forstr,ingeo,ionion_realspace,ionion_surface
!!      m_abihist,m_crystal,m_ddb,m_effective_potential,m_fit_polynomial_coeff
!!      m_mep,m_pred_lotf,m_results_img,m_tdep_abitypes,make_efg_el
!!      make_efg_ion,mkcore_paw,mkcore_wvl,mkgrid_fft,mklocl,mklocl_realspace
!!      mlwfovlp_projpaw,mover_effpot,out1dm,outqmc,outvar_o_z,outxml
!!      pimd_langevin_npt,pimd_langevin_nvt,pimd_nosehoover_npt
!!      pimd_nosehoover_nvt,prec_simple,pred_delocint,pred_diisrelax,pred_hmc
!!      pred_isokinetic,pred_isothermal,pred_langevin,pred_moldyn,pred_nose
!!      pred_srkna14,pred_steepdesc,pred_velverlet,pred_verlet,prtimg
!!      prtspgroup,prtxfase,randomcellpos,rhotov,setvtr,spin_current,symspgr
!!      thmeig,vso_realspace_local,vtorho,wrt_moldyn_netcdf,wvl_denspot_set
!!      wvl_initro,wvl_memory,wvl_nhatgrid,wvl_projectors_set,wvl_rwwf
!!      wvl_setboxgeometry,wvl_wfs_set,wvl_wfsinp_reformat,wvl_wfsinp_scratch
!!      xfh_recover_deloc
!!
!! CHILDREN
!!
!! SOURCE

subroutine xred2xcart(natom,rprimd,xcart,xred)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
!arrays
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom)
 real(dp),intent(out) :: xcart(3,natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,mu

! *************************************************************************

 do iatom=1,natom
   do mu=1,3
     xcart(mu,iatom)=rprimd(mu,1)*xred(1,iatom)+rprimd(mu,2)*xred(2,iatom)+rprimd(mu,3)*xred(3,iatom)
   end do
 end do

end subroutine xred2xcart
!!***

!!****f* m_geometry/fred2fcart
!! NAME
!! fred2fcart
!!
!! FUNCTION
!! Convert reduced forces into cartesian forces
!!
!! INPUTS
!!  fred(3,natom)=symmetrized grtn = d(etotal)/d(xred)
!!  natom=Number of atoms in the unitary cell
!!  Favgz_null=TRUE if the average cartesian force has to be set to zero
!!             FALSE if it is set to zero only in x,y directions (not z)
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1)
!!
!! OUTPUT
!!  fcart(3,natom)=forces in cartesian coordinates (Ha/Bohr)
!!
!! NOTES
!!    Unlike fred, fcart has been corrected by enforcing
!!    the translational symmetry, namely that the sum of force
!!    on all atoms is zero (except is a slab is used)
!!
!! PARENTS
!!      forces,m_mep
!!
!! CHILDREN
!!
!! SOURCE

subroutine fred2fcart(favg,Favgz_null,fcart,fred,gprimd,natom)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
 logical :: Favgz_null
!arrays
 real(dp),intent(out) :: fcart(3,natom)
 real(dp),intent(in) :: fred(3,natom)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(out) :: favg(3)

!Local variables-------------------------------
!scalars
 integer :: iatom,mu

! *************************************************************************

!Note conversion to cartesian coordinates (bohr) AND
!negation to make a force out of a gradient
 favg(:)=zero
 do iatom=1,natom
   do mu=1,3
     fcart(mu,iatom)= - (gprimd(mu,1)*fred(1,iatom)+&
&     gprimd(mu,2)*fred(2,iatom)+&
&     gprimd(mu,3)*fred(3,iatom))
     favg(mu)=favg(mu)+fcart(mu,iatom)
   end do
 end do

!Subtract off average force from each force component
!to avoid spurious drifting of atoms across cell.
 favg(:)=favg(:)/dble(natom)
 if(.not.Favgz_null) favg(3)=zero
 do iatom=1,natom
   fcart(:,iatom)=fcart(:,iatom)-favg(:)
 end do

end subroutine fred2fcart
!!***

!!****f* m_geometry/fcart2fred
!!
!! NAME
!! fcart2fred
!!
!! FUNCTION
!! Convert cartesian forces into reduced forces
!!
!! INPUTS
!!  fcart(3,natom)=forces in cartesian coordinates (Ha/Bohr)
!!  natom=Number of atoms in the unitary cell
!!  rprimd(3,3)=dimensional primitive
!!
!! OUTPUT
!!  fred(3,natom)=symmetrized grtn = d(etotal)/d(xred)
!!
!! NOTES
!!  Unlike fred, fcart has been corrected by enforcing
!!  the translational symmetry, namely that the sum of force
!!  on all atoms is zero.
!!
!! PARENTS
!!      gstateimg,m_abihist,m_effective_potential,m_mep,mover,prec_simple
!!      pred_bfgs,pred_delocint,pred_lbfgs,pred_verlet,prtxfase
!!
!! CHILDREN
!!
!! SOURCE

subroutine fcart2fred(fcart,fred,rprimd,natom)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
!arrays
 real(dp),intent(in) :: fcart(3,natom)
 real(dp),intent(out) :: fred(3,natom)
 real(dp),intent(in) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: iatom,mu

! *************************************************************************

!MT, april 2012: the coding was not consistent with fred2fcart
 do iatom=1,natom
   do mu=1,3
     fred(mu,iatom)= - (rprimd(1,mu)*fcart(1,iatom)+&
&     rprimd(2,mu)*fcart(2,iatom)+&
&     rprimd(3,mu)*fcart(3,iatom))
   end do
 end do

!Previous version
!do iatom=1,natom
!do mu=1,3
!fred(mu,iatom)= - (rprimd(mu,1)*fcart(1,iatom)+&
!&     rprimd(mu,2)*fcart(2,iatom)+&
!&     rprimd(mu,3)*fcart(3,iatom))
!end do
!end do

end subroutine fcart2fred
!!***

!!****f* m_geometry/bonds_lgth_angles
!! NAME
!! bonds_lgth_angles
!!
!! FUNCTION
!! From list of coordinates and primitive translations, output
!! a list of bonds lengths and bond angles.
!!
!! INPUTS
!!  coordn = maximum coordination number to be taken into account
!!  fnameabo_app_geo=name of file for _GEO data
!!  natom  = number of atoms in unit cell
!!  ntypat = number of types of atoms in unit cell.
!!  rprimd(3,3)  = real space dimensional primitive translations (bohr)
!!  typat(natom) = type integer for each atom in cell
!!  znucl(ntypat)= real(dp), atomic number of atom type
!!  xred(3,natom)= reduced coordinates of atoms
!!
!! OUTPUT
!! data written in file fnameabo_app_geo
!!
!! NOTES
!!  The tolerance tol8 aims at giving a machine-independent ordering.
!!  (this trick is used in bonds.f, listkk.f, prtrhomxmn.f and rsiaf9.f)
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      atomdata_from_znucl,wrtout,xred2xcart
!!
!! SOURCE

subroutine bonds_lgth_angles(coordn,fnameabo_app_geo,natom,ntypat,rprimd,typat,xred,znucl)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: coordn,natom,ntypat
 character(len=*),intent(in) :: fnameabo_app_geo
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: rprimd(3,3),znucl(ntypat)
 real(dp),intent(inout) :: xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: done,ia,ib,ic,ii,ineighb,jneighb,mneighb,mu,ndig,nu,t1,t2,t3,tmax,temp_unit
 real(dp) :: adotb,asq,bsq,co,length,sq,thdeg
!real(dp)u1,u2,u3,v1,v2,v3
 character(len=500) :: message
 type(atomdata_t) :: atom
!arrays
 integer,allocatable :: list_neighb(:,:,:)
 real(dp) :: bab(3),bac(3),dif(3),rmet(3,3)
 real(dp),allocatable :: sqrlength(:),xangst(:,:),xcart(:,:)
 character(len=8),allocatable :: iden(:)

! *************************************************************************

!Initialize the file
 write(message, '(a,a,a)' )' bonds_lgth_angles : about to open file ',trim(fnameabo_app_geo),ch10
 call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')

 if (open_file(fnameabo_app_geo,message,newunit=temp_unit,status='unknown',form='formatted') /= 0) then
   MSG_ERROR(message)
 end if
 rewind(temp_unit)

 write(message, '(a,a)' ) ch10,' ABINIT package : GEO file '
 call wrtout(temp_unit,message,'COLL')

!Compute maximum number of neighbors is the neighbor list,
!from the indicative coordination number
!Note : the following formula includes next nearest neighbors, but not others
 mneighb=1+coordn+coordn*(coordn-1)

 write(message, '(a,a,i2,a,a,i4,a,a,a,i4,a)' ) ch10,&
& ' Maximal coordination number, as estimated by the user : ',coordn,ch10,&
& '  giving a maximum of ',coordn*coordn,&
& ' nearest neighbors and next nearest neighbors, ',ch10,&
& '                  and ',(coordn*(coordn-1))/2,&
& ' distinct angles between nearest neighbors'
 call wrtout(temp_unit,message,'COLL')

!Compute metric tensor in real space rmet
 do nu=1,3
   do mu=1,3
     rmet(mu,nu)=rprimd(1,mu)*rprimd(1,nu)+&
&     rprimd(2,mu)*rprimd(2,nu)+&
&     rprimd(3,mu)*rprimd(3,nu)
   end do
 end do

 write(message, '(a,a)' )ch10,' Primitive vectors of the periodic cell (bohr)'
 call wrtout(temp_unit,message,'COLL')
 do nu=1,3
   write(message, '(1x,a,i1,a,3f10.5)' ) '  R(',nu,')=',rprimd(:,nu)
   call wrtout(temp_unit,message,'COLL')
 end do

 write(message, '(a,a)' ) ch10,&
& ' Atom list        Reduced coordinates          Cartesian coordinates (bohr)'
 call wrtout(temp_unit,message,'COLL')

!Set up a list of character identifiers for all atoms : iden(ia)
 ABI_ALLOCATE(iden,(natom))
 iden(:)='        '
 do ia=1,natom
   ndig=int(log10(dble(ia)+0.5d0))+1
   call atomdata_from_znucl(atom,znucl(typat(ia)))
   if(ndig==1) write(iden(ia), '(a,a,i1,a)' )  atom%symbol,'(',ia,')   '
   if(ndig==2) write(iden(ia), '(a,a,i2,a)' )  atom%symbol,'(',ia,')  '
   if(ndig==3) write(iden(ia), '(a,a,i3,a)' )  atom%symbol,'(',ia,') '
   if(ndig==4) write(iden(ia), '(a,a,i4,a)' )  atom%symbol,'(',ia,')'
   if(ndig>4)then
     close(temp_unit)
     write(message, '(a,i8,a,a)' )&
&     'bonds_lgth_angles cannot handle more than 9999 atoms, while natom=',natom,ch10,&
&     'Action: decrease natom, or contact ABINIT group.'
     MSG_BUG(message)
   end if
 end do

!Compute cartesian coordinates, and print reduced and cartesian coordinates
!then print coordinates in angstrom, with the format neede for xmol
 ABI_ALLOCATE(xangst,(3,natom))
 ABI_ALLOCATE(xcart,(3,natom))
 call xred2xcart(natom,rprimd,xcart,xred)
 xangst(:,:)=xcart(:,:)*Bohr_Ang

 do ia=1,natom
   write(message, '(a,a,3f10.5,a,3f10.5)' ) &
&   '   ',iden(ia),(xred(ii,ia)+tol10,ii=1,3),&
&   '    ',(xcart(ii,ia)+tol10,ii=1,3)
   call wrtout(temp_unit,message,'COLL')
 end do

 write(message, '(a,a,a,a,i4,a)' )ch10,&
& ' XMOL data : natom, followed by cartesian coordinates in Angstrom',&
& ch10,ch10,natom,ch10
 call wrtout(temp_unit,message,'COLL')

 do ia=1,natom
   call atomdata_from_znucl(atom,znucl(typat(ia)))
   write(message, '(a,a,3f10.5)' )'   ',atom%symbol,xangst(1:3,ia)
   call wrtout(temp_unit,message,'COLL')
 end do

 ABI_DEALLOCATE(xangst)
 ABI_DEALLOCATE(xcart)

 ABI_ALLOCATE(list_neighb,(0:mneighb+1,4,2))
 ABI_ALLOCATE(sqrlength,(0:mneighb+1))

!Compute list of neighbors
 do ia=1,natom

   write(message, '(a,a,a,a,a,a,a,a,a)' ) ch10,'===========',&
&   '=====================================================================',&
&   ch10,' ',iden(ia),ch10,ch10,' Bond lengths '
   call wrtout(temp_unit,message,'COLL')

!  Search other atoms for bonds, but must proceed
!  in such a way to consider a search box sufficiently large,
!  so increase the size of the search box until the
!  final bond length list do not change
   do tmax=0,5

!    Set initial list of neighbors to zero,
!    and initial square of bond lengths to a very large number.
!    Note that the dimension is larger than neighb to ease
!    the later sorting : neighbors 0 and neighb+1 are non-existent, while
!    neighbor 1 will be the atom itself ...
     list_neighb(0:mneighb+1,1:4,1)=0
     sqrlength(1:mneighb+1)=huge(0.0d0)
     sqrlength(0)=-1.0d0

!    Here search on all atoms inside the box defined by tmax
     do ib=1,natom
       do t3=-tmax,tmax
         do t2=-tmax,tmax
           do t1=-tmax,tmax
             dif(1)=xred(1,ia)-(xred(1,ib)+dble(t1))
             dif(2)=xred(2,ia)-(xred(2,ib)+dble(t2))
             dif(3)=xred(3,ia)-(xred(3,ib)+dble(t3))
             sq=rsdot(dif(1),dif(2),dif(3),dif(1),dif(2),dif(3),rmet)

!            Insert the atom at the proper place in the neighbor list.
             do ineighb=mneighb,0,-1
!              Note the tolerance
               if(sq+tol8>sqrlength(ineighb))then
                 sqrlength(ineighb+1)=sq
                 list_neighb(ineighb+1,1,1)=ib
                 list_neighb(ineighb+1,2,1)=t1
                 list_neighb(ineighb+1,3,1)=t2
                 list_neighb(ineighb+1,4,1)=t3
!                DEBUG
!                if(ineighb/=mneighb)then
!                write(std_out,*)' '
!                do ii=1,mneighb
!                write(std_out,*)ii,sqrlength(ii)
!                end do
!                end if
!                ENDDEBUG
                 exit
               else
                 sqrlength(ineighb+1)=sqrlength(ineighb)
                 list_neighb(ineighb+1,1:4,1)=list_neighb(ineighb,1:4,1)
               end if
             end do

           end do
         end do
       end do
!      end ib loop:
     end do

!    Now, check that the box defined by tmax was large enough :
!    require the present and old lists to be the same
     done=0

     if(tmax>0)then
       done=1
       do ineighb=1,mneighb
!        DEBUG
!        write(std_out,'(5i5,f12.5)' )ineighb,list_neighb(ineighb,1:4,1),&
!        &                                    sqrlength(ineighb)
!        write(std_out,'(5i5)' )ineighb,list_neighb(ineighb,1:4,2)
!        ENDDEBUG
         if( list_neighb(ineighb,1,1)/=list_neighb(ineighb,1,2) .or. &
&         list_neighb(ineighb,2,1)/=list_neighb(ineighb,2,2) .or. &
&         list_neighb(ineighb,3,1)/=list_neighb(ineighb,3,2) .or. &
&         list_neighb(ineighb,4,1)/=list_neighb(ineighb,4,2)       )then
           done=0
         end if
       end do
     end if

!    If done==1, then one can exit the loop : the correct list of
!    neighbors is contained in list_neighb(1:neighb,1:4,1),
!    with the first neighbor being the atom itself
     if(done==1)exit

!    If the work is not done, while tmax==5, then there is a problem .
     if(tmax==5)then
       close(temp_unit)
       write(message, '(2a)' )&
&       'Did not succeed to generate a reliable list of bonds ',&
&       'since tmax is exceeded.'
       MSG_BUG(message)
     end if

!    Copy the new list into the old list.
     list_neighb(1:mneighb,1:4,2)=list_neighb(1:mneighb,1:4,1)

!    Loop on tmax (note that there are exit instruction inside the loop)
   end do



!  Output the bond list
   do ineighb=2,mneighb
     ib=list_neighb(ineighb,1,1)
     length=sqrt(sqrlength(ineighb))
     write(message, '(a,a,a,a,3i2,t27,a,f10.5,a,f9.5,a)' )&
&     '  ',trim(iden(ia)),' - ',trim(iden(ib)),&
&     list_neighb(ineighb,2:4,1),'bond length is ',&
&     length,' bohr  ( or ',Bohr_Ang*length,' Angst.)'
     call wrtout(temp_unit,message,'COLL')
   end do

!  Output the angle list
   if(coordn>1)then

     write(message, '(a,a)' ) ch10,' Bond angles '
     call wrtout(temp_unit,message,'COLL')

     do ineighb=2,coordn
       do jneighb=ineighb+1,coordn+1

         ib=list_neighb(ineighb,1,1)
         ic=list_neighb(jneighb,1,1)
         do mu=1,3
           bab(mu)=xred(mu,ib)+dble(list_neighb(ineighb,1+mu,1))-xred(mu,ia)
           bac(mu)=xred(mu,ic)+dble(list_neighb(jneighb,1+mu,1))-xred(mu,ia)
         end do
         asq=rsdot(bab(1),bab(2),bab(3),bab(1),bab(2),bab(3),rmet)
         bsq=rsdot(bac(1),bac(2),bac(3),bac(1),bac(2),bac(3),rmet)
         adotb=rsdot(bab(1),bab(2),bab(3),bac(1),bac(2),bac(3),rmet)
         co=adotb/sqrt(asq*bsq)
         if( abs(co)-1.0d0 >= 0.0d0 )then
           if( abs(co)-1.0d0 <= 1.0d-12 )then
!            Allows for a small numerical inaccuracy
             thdeg=0.0d0
             if(co < 0.0d0) thdeg=180.0d0
           else
             MSG_BUG('the evaluation of the angle is wrong.')
           end if
         else
           thdeg=acos(co)*180.d0*piinv
         end if

         write(message, '(a,a,3i2,a,a,a,a,3i2,t44,a,f13.5,a)' )&
&         '  ',trim(iden(ib)),list_neighb(ineighb,2:4,1),' - ',&
&         trim(iden(ia)),' - ',trim(iden(ic)),&
&         list_neighb(jneighb,2:4,1),'bond angle is ',thdeg,' degrees '
         call wrtout(temp_unit,message,'COLL')
       end do
     end do

   end if
 end do !  End big ia loop:

 ABI_DEALLOCATE(iden)
 ABI_DEALLOCATE(list_neighb)
 ABI_DEALLOCATE(sqrlength)

 close(temp_unit)

 contains

   function rsdot(u1,u2,u3,v1,v2,v3,rmet)

   real(dp) :: rsdot
   real(dp),intent(in) :: u1,u2,u3,v1,v2,v3
   real(dp),intent(in) :: rmet(3,3)
   rsdot=rmet(1,1)*u1*v1+rmet(2,1)*u2*v1+&
&   rmet(3,1)*u3*v1+rmet(1,2)*u1*v2+rmet(2,2)*u2*v2+&
&   rmet(3,2)*u3*v2+rmet(1,3)*u1*v3+rmet(2,3)*u2*v3+rmet(3,3)*u3*v3
 end function rsdot

end subroutine bonds_lgth_angles
!!***

!!****f* m_geometry/randomcellpos
!! NAME
!!  randomcellpos
!!
!! FUNCTION
!!  This subroutine creates a unit cell with random atomic positions. It is
!!  assumed that the cell parameters are given and fixed. Several methods are
!!  used to generate the cell.
!!
!! INPUTS
!! natom=number of atoms
!! npsp=number of pseudopotentials (needed for the dimension of znucl)
!! ntypat=number of type of atoms
!! random_atpos=input variable
!!   0 no generation of random atomic potision
!!   1 completely random atomic potisions
!!   2 random atomic positions, avoiding too close atoms
!!     (prevent coming closer than a fraction of the sum of covalent radii)
!!   3 same than 2 but also generates the rprim and acell randomly
!!    within some given ranges (angles between 50 and 130)
!! ratsph(1:ntypat)=radius of the atomic sphere
!! rprimd(3,3)=dimensional primitive translations in real space (bohr)
!! typat(1:natom)= input variable giving the type of each atom
!! znucl(1:npsp)=nuclear number of atom as specified in psp file
!!
!! OUTPUT
!! xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      ingeo
!!
!! CHILDREN
!!      atomdata_from_znucl,mkrdim,xred2xcart
!!
!! SOURCE

subroutine randomcellpos(natom,npsp,ntypat,random_atpos,ratsph,rprim,rprimd,typat,xred,znucl,acell)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,npsp,ntypat,random_atpos
!arrays
 integer, intent(in)   :: typat(natom)
 real(dp),intent(in)   :: ratsph(ntypat)
 real(dp), intent(inout)  :: rprim(3,3)
 real(dp), intent(inout)  :: rprimd(3,3)
 real(dp), intent(inout) :: xred(3,natom)
 real(dp), intent(in) :: znucl(npsp)
 real(dp), intent(inout) :: acell(3)

!Local variables-------------------------------
 integer ::   iatom=0,ii,idum=-20
 real(dp) ::  rij(3), rijd(3), radiuscovi, radiuscovj, dist, rati, ratj, angdeg(3)
 real(dp) ::  cosang,aa,cc,a2
 character(len=500) :: message
 type(atomdata_t) :: atom

! *************************************************************************

!DEBUG
!For the time being, print rprimd to keep it as an argument, in spite of abirule checking.
!write (std_out,*) ' randomcellpos : enter'
!write(std_out,*)' rprimd=',rprimd
!write(std_out,*)' znucl=',znucl
!write(std_out,*)' typat=',typat
!write(std_out,*)' random_atpos=',random_atpos
!ENDDEBUG

 if(random_atpos==2 .and. npsp/=ntypat)then
   write(message, '(a,i5,2a,i5,a,i5,4a)' )&
&   'Input variable random_atpos= ',random_atpos,ch10,&
&   'However, the number of pseudopotentials ',npsp,', is not equal to the number of type of atoms ',ntypat,ch10,&
&   'The use of alchemical mixing cannot be combined with the constraint based on the mixing of covalent radii.',ch10,&
&   'Action: switch to another value of random_atpos.'
   MSG_ERROR(message)
 end if

!random_atpos = 0   Default value, no random initialisation
!random_atpos = 1   Fully random (Is it really useful ???)
!random_atpos = 2   Random, but the sum of the two covalent radii is
!less than the interatomic distance
!random_atpos = 3   Random, but the sum of the two (other type of)
!radii is less than the interatomic distance
!random_atpos = 4   Random, but the sum of the two pseudopotential
!radii is less than the interatomic distance
!random_atpos = 5   Random, but the interatomic distance must be bigger
!than the sum of
!some input variable (well, instead of defining a new variable, why
!not use ratsph ?)
!Right now we are not using a factor for the tested distance.. something to be done, after a new variable has been defined

 if (random_atpos /= 0) then
   select case (random_atpos)
   case (1)
     do ii=1,natom
       xred(1,ii)=uniformrandom(idum)
       xred(2,ii)=uniformrandom(idum)
       xred(3,ii)=uniformrandom(idum)
     end do
   case (2)
     iatom=0
     do
       iatom=iatom+1
       xred(1,iatom)=uniformrandom(idum)
       xred(2,iatom)=uniformrandom(idum)
       xred(3,iatom)=uniformrandom(idum)
       call atomdata_from_znucl(atom,znucl(typat(iatom)))
       radiuscovi = atom%rcov
       do ii=1,iatom-1
         rij=xred(:,iatom)-xred(:,ii)
!          periodic boundary conditions
         rij = rij - 0.5
         rij = rij - anint (rij)
!          coming back to cube between (0,1)
         rij = rij + 0.5
!          convert reduced coordinates to cartesian coordinates
         call xred2xcart(1,rprimd,rijd,rij)
         dist=dot_product(rijd,rijd)
         call atomdata_from_znucl(atom,znucl(typat(ii)))
         radiuscovj = atom%rcov
         if (dist<(radiuscovj+radiuscovi)) then
           iatom = iatom -1
           EXIT
         end if
       end do
       if (iatom>=natom) EXIT
     end do
   case(3)
     iatom=0
     do
       iatom=iatom+1
       xred(1,iatom)=uniformrandom(idum)
       xred(2,iatom)=uniformrandom(idum)
       xred(3,iatom)=uniformrandom(idum)
       call atomdata_from_znucl(atom,znucl(typat(iatom)))
       radiuscovi = atom%rcov
       do ii=1,iatom-1
         rij=xred(:,iatom)-xred(:,ii)
!          periodic boundary conditions
         rij = rij - 0.5
         rij = rij - anint (rij)
!          coming back to cube between (0,1)
         rij = rij + 0.5
!          convert reduced coordinates to cartesian coordinates
         call xred2xcart(1,rprimd,rijd,rij)
         dist=dot_product(rijd,rijd)
         call atomdata_from_znucl(atom,znucl(typat(ii)))
         radiuscovj = atom%rcov
         if (dist<(radiuscovj+radiuscovi)) then
           iatom = iatom -1
           EXIT
         end if
       end do
       if (iatom>=natom) EXIT
     end do
     do ii=1,3
!        generates cells with angles between 60 and 120 degrees
       angdeg(ii)=60_dp+uniformrandom(idum)*60.0_dp
     end do
     if (angdeg(1)+angdeg(2)+angdeg(3)>360._dp) then
       angdeg(3)=360._dp-angdeg(1)-angdeg(2)
     end if
!      check if angles are between the limits and create rprim
     if( abs(angdeg(1)-angdeg(2))<tol12 .and. &
&     abs(angdeg(2)-angdeg(3))<tol12 .and. &
&     abs(angdeg(1)-90._dp)+abs(angdeg(2)-90._dp)+abs(angdeg(3)-90._dp)>tol12 )then
!        Treat the case of equal angles (except all right angles) :
!        generates trigonal symmetry wrt third axis
       cosang=cos(pi*angdeg(1)/180.0_dp)
       a2=2.0_dp/3.0_dp*(1.0_dp-cosang)
       aa=sqrt(a2)
       cc=sqrt(1.0_dp-a2)
       rprim(1,1)=aa        ; rprim(2,1)=0.0_dp                 ; rprim(3,1)=cc
       rprim(1,2)=-0.5_dp*aa ; rprim(2,2)= sqrt(3.0_dp)*0.5_dp*aa ; rprim(3,2)=cc
       rprim(1,3)=-0.5_dp*aa ; rprim(2,3)=-sqrt(3.0_dp)*0.5_dp*aa ; rprim(3,3)=cc
!        DEBUG
!        write(std_out,*)' ingeo : angdeg=',angdeg(1:3)
!        write(std_out,*)' ingeo : aa,cc=',aa,cc
!        ENDDEBUG
     else
!        Treat all the other cases
       rprim(:,:)=0.0_dp
       rprim(1,1)=1.0_dp
       rprim(1,2)=cos(pi*angdeg(3)/180.0_dp)
       rprim(2,2)=sin(pi*angdeg(3)/180.0_dp)
       rprim(1,3)=cos(pi*angdeg(2)/180.0_dp)
       rprim(2,3)=(cos(pi*angdeg(1)/180.0_dp)-rprim(1,2)*rprim(1,3))/rprim(2,2)
       rprim(3,3)=sqrt(1.0_dp-rprim(1,3)**2-rprim(2,3)**2)
     end if
!      generate acell
     aa=zero
     do ii=1,npsp
       aa=znucl(ii)
     end do
     do ii=1,3
       acell(ii)=aa+uniformrandom(idum)*4.0
     end do
     call mkrdim(acell,rprim,rprimd)
   case(4)
     write(std_out,*) 'Not implemented yet'
   case(5)
     iatom=0
     do
       iatom=iatom+1
       xred(1,iatom)=uniformrandom(idum)
       xred(2,iatom)=uniformrandom(idum)
       xred(3,iatom)=uniformrandom(idum)
       rati=ratsph(typat(iatom))
       do ii=1,iatom-1
         ratj=ratsph(typat(ii))
!          apply periodic boundary conditions
         rij=(xred(:,iatom)-xred(:,ii))-0.5
         rij = rij - ANINT ( rij )
         rij = rij + 0.5
         call xred2xcart(natom,rprimd,rijd,rij)
         dist=dot_product(rijd,rijd)
         if (dist<(rati+ratj)) EXIT
       end do
       if (iatom==natom) EXIT
       if (ii<(iatom-1)) iatom=iatom-1
     end do
   end select
 end if

end subroutine randomcellpos
!!***

!!****f* m_geometry/shellstruct
!! NAME
!!  shellstruct
!!
!! FUNCTION
!!  Calculates shell structure (multiplicities, radii) of an atomic configuration
!!
!! INPUTS
!!  natom=number of atoms in unit cell
!!  xred=reduced coordinates of atoms
!!  rprimd=unit cell vectors
!!  magv = magnetic ordering of atoms given as 1 and -1, if not given fm is assumed
!!  atp = atom on which the perturbation was done
!!
!! OUTPUT
!!  sdisv(nat)= distance of each shell to central atom (only the first nsh entries are relevant)
!!  nsh= number of shells
!!  mult(nat) = number of atoms on shell (only the first nsh entries are relevant)
!!
!! PARENTS
!!      pawuj_det
!!
!! CHILDREN
!!      ioniondist,prmat,sort_dp,sort_int,wrtout
!!
!! SOURCE

subroutine shellstruct(xred,rprimd,natom,magv,distv,smult,sdisv,nsh,atp,prtvol)


!Arguments ------------------------------------
!scalars
 integer,intent(in)              :: natom
 integer,intent(in),optional     :: atp
 integer,intent(in),optional     :: prtvol
 integer,intent(out)             :: nsh
!arrays
 real(dp),intent(in)             :: rprimd(3,3)
 real(dp),intent(in)             :: xred(3,natom)
 integer,intent(out)             :: smult(natom)
 integer,intent(in),optional     :: magv(natom)
 real(dp),intent(out)            :: sdisv(natom)
 real(dp),intent(out)            :: distv(natom)

!Local variables-------------------------------
!scalars
 integer                      :: iatom,atpp,ish,prtvoll
 character(len=500)           :: message
 real(dp),parameter           :: rndfact=10000_dp
!arrays
 integer                      :: iperm(natom),jperm(natom)
 real(dp)                     :: distvh(natom,natom)
 real(dp)                     :: magvv(natom)

! *************************************************************************

 if (present(magv)) then
   magvv=magv
 else
   magvv=(/ (1, iatom=1,natom)  /)
 end if

 if (present(atp)) then
   atpp=atp
 else
   atpp=1
 end if

 if (present(prtvol)) then
   prtvoll=prtvol
 else
   prtvoll=1
 end if

!DEBUB
 write(std_out,*)'shellstruct start'
!END DEBUG

!Calculate ionic distances
 call ioniondist(natom,rprimd,xred,distvh,1,magv=int(magvv),atp=atpp)
 distv=distvh(1,:)

 if (prtvol>2) then
   write(std_out,'(a)')' shellstruct ionic distances in cell (distv) : '
   call prmat(distv(1:natom),1,natom,1,std_out)
 end if

 iperm=(/ (iatom, iatom=1,natom ) /)
 jperm=iperm
 distv=anint(distv*rndfact)/rndfact
!Sort distances
 call sort_dp(natom,distv,iperm,10d-5)
 call sort_int(natom,iperm,jperm)

 smult=0
 sdisv=dot_product(rprimd(1,:),rprimd(1,:))+dot_product(rprimd(2,:),rprimd(2,:))+dot_product(rprimd(3,:),rprimd(3,:))

 nsh=1
 smult(1)=1
 sdisv(1)=distv(1)

 do iatom=2,natom
   do ish=1,natom
     if (distv(iatom)>sdisv(ish)) then
       cycle
     else if (distv(iatom)==sdisv(ish)) then
       smult(ish)=smult(ish)+1
       exit
     else if (distv(iatom)<sdisv(ish)) then
       smult(ish+1:natom)=smult(ish:natom-1)
       sdisv(ish+1:natom)=sdisv(ish:natom-1)
       smult(ish)=1
       sdisv(ish)=distv(iatom)
       nsh=nsh+1
       exit
     end if
   end do
 end do

 distv=(/ ( distv(jperm(iatom)),iatom=1,natom ) /)

 if (prtvoll>2) then
   write(message,'(a,i4,a)')' shellstruct found ',nsh,' shells at distances (sdisv) '
   call wrtout(std_out,message,'COLL')
   call prmat(sdisv(1:nsh),1,nsh,1,std_out)
   write(message,fmt='(a,150i4)')' and multiplicities (smult) ', smult(1:nsh)
   call wrtout(std_out,message,'COLL')
 end if

!DEBUB
 write(std_out,*)'shellstruct leave'
!END DEBUG

end subroutine shellstruct
!!***

!!****f* m_geometry/ioniondist
!! NAME
!! ioniondist
!!
!! FUNCTION
!!  Compute ion-ion distances
!!
!! INPUTS
!!  natom= number of atoms in unit cell
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  xred(3,natom)=dimensionless reduced coordinates of atoms
!!  inm(natom,natom)=index (m,n) of the atom
!!  option= 1 output ion-ion distances / 2 output ordering of ion-ion
!!          distances / 3 output variables in varlist
!!          according to ion-ion distances * magnetic ordering
!!          magv magnetic ordering of atoms given als 1 and -1, if not
!!          given fm is assumed
!!  varlist=List of variables
!!  magv(natom)= magnetic ordering of atoms
!!  atp=atom on which the perturbation was done
!!
!! OUTPUT
!!
!! PARENTS
!!      pawuj_utils,shellstruct
!!
!! CHILDREN
!!      prmat,wrtout
!!
!! SOURCE

subroutine ioniondist(natom,rprimd,xred,inm,option,varlist,magv,atp,prtvol)


!Arguments ------------------------------------
!scalars
 integer,intent(in)              :: natom,option
 integer,intent(in),optional     :: atp                   !atom on which the perturbation was done
!arrays
 real(dp),intent(in)             :: rprimd(3,3)
 real(dp),intent(in)             :: xred(3,natom)
 real(dp),intent(out)            :: inm(natom,natom)
 integer,intent(in),optional     :: magv(natom)
 real(dp),intent(in),optional    :: varlist(natom)
 integer,intent(in),optional     :: prtvol

!Local variables-------------------------------
!scalars
 integer                      :: iatom,jatom,katom,kdum,atpp,prtvoll
 !character(len=500)           :: message
!arrays
 integer                      :: interq(natom)
 real(dp)                     :: hxcart(3,natom),distm(natom,natom)
 real(dp)                     :: magvv(natom)

! *************************************************************************

 hxcart=matmul(rprimd,xred)
 interq=(/(iatom,iatom=1,natom)/)
 inm=0

 if (present(magv)) then
   magvv=magv
 else
   magvv=(/ (1, iatom=1,natom)  /)
 end if

 if (present(atp)) then
   atpp=atp
 else
   atpp=1
 end if

 if (present(prtvol)) then
   prtvoll=prtvol
 else
   prtvoll=1
 end if

 if (option==3.and.(.not.present(varlist))) then
   call  wrtout(std_out,'ioniondist error: option=3 but no variable list provided for symmetrization','COLL')
   return
 end if


!DEBUG
!write(message, '(a,a)' ) ch10,' ioniondist start '
!call wrtout(std_out,message,'COLL')
!END DEBUG

 distm=0
 katom=atpp-1
 do iatom=1,natom
   katom=katom+1
   if (katom > natom) katom=1
   distm(iatom,iatom)=0
   do jatom=iatom,natom
     distm(iatom,jatom)=dist2(xred(:,katom),xred(:,jatom),rprimd,1)*magvv(katom)*magvv(jatom)
     distm(jatom,iatom)=distm(iatom,jatom)
   end do
 end do

 if (prtvoll>=3) then
   call  wrtout(std_out,'ioniondist: ionic distances:','COLL')
   call prmat(distm,natom,natom,natom,std_out)
 end if

 distm=anint(distm*10000_dp)/10000_dp           ! rounding needed else distm(iatom,jatom)/= distm(1,kdum) sometimes fails

 do iatom=1,natom
   if (option==1) then
     inm(iatom,:)=distm(iatom,:)
   else
     do jatom=iatom,natom
       kdum=1
       do while ( (kdum <= natom) .and. (distm(iatom,jatom)/= distm(1,kdum)) )
         kdum=kdum+1
       end do
       if (option==2) then
         inm(iatom,jatom)=interq(kdum)
       else if (option==3) then
         inm(iatom,jatom)=varlist(kdum)
       end if
       inm(jatom,iatom)=inm(iatom,jatom)
     end do
   end if
 end do

 if (prtvoll==2) then
   call wrtout(std_out,'ioniondist: symmetrized matrix:','COLL')
   call prmat(distm,1,natom,natom,std_out)
 else if (prtvoll>=3) then
   call wrtout(std_out,'ioniondist: symmetrized matrix:','COLL')
   call prmat(distm,natom,natom,natom,std_out)
 end if

end subroutine ioniondist
!!***

!!****f* m_geometry/dist2
!! NAME
!!  dist2
!!
!! FUNCTION
!!  Calculates the distance of v1 and v2 in a crystal by epeating the unit cell
!!
!! INPUTS
!!  v1,v2
!!  rprimd: dimensions of the unit cell. if not given 1,0,0/0,1,0/0,0,1 is assumed
!!  option: 0 v1, v2 given in cartesian coordinates (default) / 1 v1,v2 given in reduced coordinates
!!
!! OUTPUT
!!  dist2
!!
!! PARENTS
!!  ioniondist
!!
!! CHILDREN
!!
!! SOURCE

function dist2(v1,v2,rprimd,option)


!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: option
 real(dp) :: dist2
!arrays
 real(dp),intent(in),optional :: rprimd(3,3)
 real(dp),intent(in) :: v1(3),v2(3)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,opt,s1,s2,s3
 real(dp):: min2,norm2,ucvol
!arrays
 integer :: limits(3)
 real(dp) :: corner(3),dred(3),dtot(3),dv(3),dwrap(3),sh(3)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp) :: vprimd(3,3)

! *************************************************************************

 if (.not.PRESENT(rprimd)) then
   vprimd=reshape((/1,0,0,  0,1,0,  0,0,1/),(/3,3/))
 else
   vprimd=rprimd
 end if

 call metric(gmet,gprimd,-1,rmet,vprimd,ucvol)

 dv(:)=v2(:)-v1(:)

!If in cartesian coordinates, need to be transformed to reduced coordinates.
 opt=0
 if(present(option))then
   opt=option
 end if
 if(opt==0)then
   dred(:)=gprimd(1,:)*dv(1)+gprimd(2,:)*dv(2)+gprimd(3,:)*dv(3)
 else
   dred(:)=dv(:)
 end if

!Wrap in the ]-1/2,1/2] interval
 call wrap2_pmhalf(dred(1),dwrap(1),sh(1))
 call wrap2_pmhalf(dred(2),dwrap(2),sh(2))
 call wrap2_pmhalf(dred(3),dwrap(3),sh(3))

!Compute the limits of the parallelipipedic box that contains the Wigner-Seitz cell
!The reduced coordinates of the corners of the Wigner-Seitz cell are computed (multiplied by two)
!Then, the maximal values of these reduced coordinates are stored.
 limits(:)=0
 do s1=-1,1,2
   do s2=-1,1,2
     do s3=-1,1,2
       corner(:)=gmet(:,1)*s1*rmet(1,1)+gmet(:,2)*s2*rmet(2,2)+gmet(:,3)*s3*rmet(3,3)
       limits(1)=max(limits(1),ceiling(abs(corner(1))+tol14))
       limits(2)=max(limits(2),ceiling(abs(corner(2))+tol14))
       limits(3)=max(limits(3),ceiling(abs(corner(3))+tol14))
     end do
   end do
 end do

!Use all relevant primitive real space lattice vectors to find the minimal difference vector
 min2=huge(zero)
 do i1=-limits(1),limits(1)
   do i2=-limits(2),limits(2)
     do i3=-limits(3),limits(3)
       dtot(1)=dwrap(1)+i1
       dtot(2)=dwrap(2)+i2
       dtot(3)=dwrap(3)+i3
       norm2=dtot(1)*rmet(1,1)*dtot(1)+dtot(2)*rmet(2,2)*dtot(2)+dtot(3)*rmet(3,3)*dtot(3)+&
&       2*(dtot(1)*rmet(1,2)*dtot(2)+dtot(2)*rmet(2,3)*dtot(3)+dtot(3)*rmet(3,1)*dtot(1))
       min2=min(norm2,min2)
     end do
   end do
 end do
 dist2=sqrt(min2)

end function dist2
!!***

!!****f* m_geometry/remove_inversion
!! NAME
!! remove_inversion
!!
!! FUNCTION
!!  Remove the inversion symmetry from a symmetry set as well
!!  all the improper rotations (if present)
!!
!! INPUTS
!!  nsym=initial number of symmetries
!!  symrel(3,3,nsym)=Initial set of symmetry operarations in real space
!!  tnons(3,nsym)=Initial fractional translations
!!
!! OUTPUT
!!  nsym_out=Number of symmetries in the set without improper rotation
!!  symrel_out(:,:) [pointer] = output symmetries without improper rotations
!!  tnons_out(:) [pointer] = fractional translations associated to symrel_out
!!  pinv=-1 if the inversion has been removed, 1 otherwise
!!
!! NOTES
!!  Note the use of pointers, memory is allocated inside the procedure and passed back
!!  to the caller. Thus memory deallocation is relegated to the caller. To be on the safe side
!!  the pointers should be nullified before entering.
!!
!! PARENTS
!!      m_crystal,m_io_kss,outkss
!!
!! CHILDREN
!!      set2unit,symdet,wrtout
!!
!! SOURCE

subroutine remove_inversion(nsym,symrel,tnons,nsym_out,symrel_out,tnons_out,pinv)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(out) :: nsym_out,pinv
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 integer,pointer :: symrel_out(:,:,:)
 real(dp),intent(in) :: tnons(3,nsym)
 real(dp),pointer :: tnons_out(:,:)

!Local variables-------------------------------
!scalars
 integer :: is,is2,is_discarded,is_inv,is_retained,nsym2
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: determinant(nsym),inversion(3,3),symrel2(3,3,nsym)
 real(dp) :: dtnons(3),tnons2(3,nsym)

! *********************************************************************

 MSG_WARNING('Removing inversion related symmetrie from initial set')

 ! Find the occurence of the inversion symmetry.
 call set2unit(inversion) ; inversion=-inversion

 is_inv=0; found=.FALSE.
 do while (is_inv<nsym .and. .not.found)
   is_inv=is_inv+1; found=ALL(symrel(:,:,is_inv)==inversion)
 end do
 if (found) then
   write(msg,'(a,i3)')' The inversion is symmetry operation no. ',is_inv
 else
   write(msg,'(a)')' The inversion was not found in the symmetries list.'
 end if
 call wrtout(std_out,msg,'COLL')

 ! Find the symmetries that are related through the inversion symmetry
 call symdet(determinant,nsym,symrel)
 nsym2=0
 do is=1,nsym-1
   do is2=is+1,nsym

     dtnons(:)=tnons(:,is2)-tnons(:,is)-tnons(:,is_inv)
     found=ALL(symrel(:,:,is)==-symrel(:,:,is2)).and.isinteger(dtnons,tol8)

     if (found) then
       nsym2=nsym2+1
       ! Retain symmetries with positive determinant
       if (ALL(tnons(:,is2)<tol8).and.ALL(tnons(:,is)<tol8)) then
         is_retained=is2 ; is_discarded=is
         if (determinant(is)==1) then
           is_retained=is  ; is_discarded=is2
         end if
       else if (ALL(tnons(:,is2)<tol8)) then
         is_retained=is2 ; is_discarded=is
       else
         is_retained=is ;  is_discarded=is2
       end if

       symrel2(:,:,nsym2)=symrel(:,:,is_retained)
       tnons2   (:,nsym2)=tnons   (:,is_retained)
       write(msg,'(a,i3,a,i3,3a,i3,a)')&
&       ' Symmetry operations no. ',is,' and no. ',is2,&
&       ' are related through the inversion.',ch10,&
&       ' Symmetry operation no. ',is_discarded,' will be suppressed.'
       call wrtout(std_out,msg,'COLL')
     end if ! found

   end do !is2
 end do !is

 if (nsym2/=(nsym/2).or.nsym==1) then
   call wrtout(std_out, ' Program uses the original set of symmetries ', 'COLL')
   nsym_out=nsym
   ABI_ALLOCATE(symrel_out,(3,3,nsym))
   ABI_ALLOCATE(tnons_out,(3,nsym))
   symrel_out(:,:,:)=symrel(:,:,1:nsym)
   tnons_out(:,:)=tnons(:,1:nsym)
   pinv=1
 else
   write(msg,'(a)')' Inversion related operations have been suppressed from symmetries list.'
   call wrtout(std_out,msg,'COLL')
   nsym_out=nsym2
   ABI_ALLOCATE(symrel_out,(3,3,nsym2))
   ABI_ALLOCATE(tnons_out,(3,nsym2))
   symrel_out(:,:,:)=symrel2(:,:,1:nsym2)
   tnons_out(:,:)=tnons(:,1:nsym2)
   pinv=-1
 end if

end subroutine remove_inversion
!!***

!!****f* m_geometry/symredcart
!! NAME
!! symredcart
!!
!! FUNCTION
!! Convert a symmetry operation from reduced coordinates (integers)
!! to cartesian coordinates (reals). Can operate in real or reciprocal space
!!
!! INPUTS
!! symred(3,3)=symmetry matrice in reduced coordinates (integers) (real or reciprocal space)
!! aprim(3,3)=real or reciprocal space dimensional primitive translations (see below)
!! bprim(3,3)=real or reciprocal space dimensional primitive translations (see below)
!!
!! OUTPUT
!! symcart(3,3)=symmetry matrice in cartesian coordinates (reals)
!!
!! NOTES
!! When aprim=rprimd and bprim=gprimd, the routine operates in real space (on a real space symmetry)
!! When aprim=gprimd and bprim=rprimd, the routine operates in reciprocal space (on a real space symmetry)
!!
!! PARENTS
!!      m_matlu,m_phonons,symrhg
!!
!! CHILDREN
!!
!! SOURCE

subroutine symredcart(aprim,bprim,symcart,symred)

!Arguments ------------------------------------
!arrays
 integer,intent(in) :: symred(3,3)
 real(dp),intent(in) :: aprim(3,3),bprim(3,3)
 real(dp),intent(out) :: symcart(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,kk
 real(dp) :: symtmp
!arrays
 real(dp) :: work(3,3)

! *************************************************************************

 work=zero
 do kk=1,3
   do jj=1,3
     symtmp=dble(symred(jj,kk))
     do ii=1,3
       work(ii,jj)=work(ii,jj)+bprim(ii,kk)*symtmp
     end do
   end do
 end do

 ! work = bprim * symred^T

 symcart=zero
 do kk=1,3
   do jj=1,3
     symtmp=work(jj,kk)
     do ii=1,3
       ! symcart = aprim * work^T = aprim * symred * bprim^T
       symcart(ii,jj)=symcart(ii,jj)+aprim(ii,kk)*symtmp
     end do
   end do
 end do

end subroutine symredcart
!!***

!!****f* m_geometry/strainsym
!! NAME
!! strainsym
!!
!! FUNCTION
!! For given order of point group, symmetrizes the strain tensor,
!! then produce primitive vectors based on the symmetrized strain.
!!
!! INPUTS
!! nsym=order of group.
!! rprimd(3,3)= primitive vectors, to be symmetrized
!! rprimd0(3,3)= reference primitive vectors, already symmetrized
!! symrel(3,3,nsym)=symmetry operators in terms of action on primitive translations
!!
!! OUTPUT
!! rprimd_symm(3,3)= symmetrized primitive vectors
!!
!! PARENTS
!!      xfpack_vin2x,xfpack_x2vin
!!
!! CHILDREN
!!      dgemm,mati3inv,matrginv
!!
!! SOURCE

subroutine strainsym(nsym,rprimd0,rprimd,rprimd_symm,symrel)

 use m_linalg_interfaces

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: rprimd(3,3),rprimd0(3,3)
 real(dp),intent(out) :: rprimd_symm(3,3)

!Local variables-------------------------------
!scalars
 integer :: isym
!arrays
 integer :: symrel_it(3,3)
 real(dp) :: rprimd0_inv(3,3),strain(3,3),strain_symm(3,3),tmp_mat(3,3)

!**************************************************************************

!copy initial rprimd input and construct inverse
 rprimd0_inv = rprimd0
 call matrginv(rprimd0_inv,3,3)

!define strain as rprimd = strain * rprimd0 (in cartesian frame)
!so strain = rprimd * rprimd0^{-1}
!transform to triclinic frame with rprimd0^{-1} * strain * rprimd0
!giving strain as rprimd0^{-1} * rprimd
 call dgemm('N','N',3,3,3,one,rprimd0_inv,3,rprimd,3,zero,strain,3)

!loop over symmetry elements to obtain symmetrized strain matrix
 strain_symm = zero
 do isym = 1, nsym

!  this loop accumulates symrel^{-1}*strain*symrel into strain_symm

!  mati3inv gives the inverse transpose of symrel
   call mati3inv(symrel(:,:,isym),symrel_it)
   call dgemm('N','N',3,3,3,one,strain,3,dble(symrel(:,:,isym)),3,zero,tmp_mat,3)
   call dgemm('T','N',3,3,3,one,dble(symrel_it),3,tmp_mat,3,one,strain_symm,3)

 end do

!normalize by number of symmetry operations
 strain_symm = strain_symm/dble(nsym)

!this step is equivalent to r_new = r_old * strain * r_old^{-1} * r_old,
!that is, convert strain back to cartesian frame and then multipy by r_old,
!to get the r_new primitive vectors

 call dgemm('N','N',3,3,3,one,rprimd0,3,strain_symm,3,zero,rprimd_symm,3)

end subroutine strainsym
!!***

!!****f* m_geometry/stresssym
!! NAME
!! stresssym
!!
!! FUNCTION
!! For given order of point group, symmetrizes the stress tensor,
!! in symmetrized storage mode and cartesian coordinates, using input
!! 3x3 symmetry operators in reduced coordinates.
!! symmetrized tensor replaces input tensor.
!!
!! INPUTS
!! gprimd(3,3)=dimensional primitive translations for reciprocal space (bohr**-1)
!! nsym=order of group.
!! sym(3,3,nsym)=symmetry operators (usually symrec=expressed in terms
!!               of action on reciprocal lattice primitive translations);
!!               integers.
!!
!! OUTPUT
!! stress(6)=stress tensor, in cartesian coordinates, in symmetric storage mode
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      dfpt_nselt,dfpt_nstpaw,forstrnps,littlegroup_pert,pawgrnl,stress
!!
!! CHILDREN
!!      matr3inv,strconv
!!
!! SOURCE

subroutine stresssym(gprimd,nsym,stress,sym)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
!arrays
 integer,intent(in) :: sym(3,3,nsym)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(inout) :: stress(6)

!Local variables-------------------------------
!scalars
 integer :: ii,isym,mu,nu
 real(dp) :: summ,tmp
!arrays
 real(dp) :: rprimd(3,3),rprimdt(3,3),strfrac(6),tensor(3,3),tt(3,3)

!*************************************************************************

!Obtain matrix of real space dimensional primitive translations
!(inverse tranpose of gprimd), and its transpose
 call matr3inv(gprimd,rprimd)
 rprimdt=transpose(rprimd)

!Compute stress tensor in reduced coordinates
 call strconv(stress,rprimdt,strfrac)

!Switch to full storage mode
 tensor(1,1)=strfrac(1)
 tensor(2,2)=strfrac(2)
 tensor(3,3)=strfrac(3)
 tensor(3,2)=strfrac(4)
 tensor(3,1)=strfrac(5)
 tensor(2,1)=strfrac(6)
 tensor(2,3)=tensor(3,2)
 tensor(1,3)=tensor(3,1)
 tensor(1,2)=tensor(2,1)

 do nu=1,3
   do mu=1,3
     tt(mu,nu)=tensor(mu,nu)/dble(nsym)
     tensor(mu,nu)=0.0_dp
   end do
 end do

!loop over all symmetry operations:
 do isym=1,nsym
   do mu=1,3
     do nu=1,3
       summ=0._dp
       do ii=1,3
         tmp=tt(ii,1)*sym(nu,1,isym)+tt(ii,2)*sym(nu,2,isym)+&
&         tt(ii,3)*sym(nu,3,isym)
         summ=summ+sym(mu,ii,isym)*tmp
       end do
       tensor(mu,nu)=tensor(mu,nu)+summ
     end do
   end do
 end do

!Switch back to symmetric storage mode
 strfrac(1)=tensor(1,1)
 strfrac(2)=tensor(2,2)
 strfrac(3)=tensor(3,3)
 strfrac(4)=tensor(3,2)
 strfrac(5)=tensor(3,1)
 strfrac(6)=tensor(2,1)

!Convert back stress tensor (symmetrized) in cartesian coordinates
 call strconv(strfrac,gprimd,stress)

end subroutine stresssym
!!***

!!****f* m_geometry/strconv
!! NAME
!! strconv
!!
!! FUNCTION
!! If original gprimd is input, convert from symmetric storage mode
!! 3x3 tensor in reduced coordinates "frac" to symmetric storage mode
!! symmetric tensor in cartesian coordinates "cart".
!!
!! INPUTS
!!  frac(6)=3x3 tensor in symmetric storage mode, reduced coordinates
!!  gprimd(3,3)=reciprocal space dimensional primitive translations (bohr^-1)
!!
!! OUTPUT
!!  cart(6)=symmetric storage mode for symmetric 3x3 tensor in cartesian coords.
!!
!! NOTES
!! $cart(i,j)=G(i,a) G(j,b) frac(a,b)$
!! "Symmetric" storage mode for 3x3 tensor is 6 element array with
!! elements 11, 22, 33, 32, 31, and 21.
!! "cart" may be same array as "frac".
!! If rprimd transpose is input instead of gprimd, then convert tensor
!! in cartesian coordinates to reduced coordinates
!!
!! PARENTS
!!      ctocprj,d2frnl,mkcore,mkcore_paw,mkcore_wvl,nonlop_pl,nonlop_ylm
!!      stresssym
!!
!! CHILDREN
!!
!! SOURCE

subroutine strconv(frac,gprimd,cart)


!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: frac(6),gprimd(3,3)
 real(dp),intent(inout) :: cart(6) ! alias of frac   !vz_i

!Local variables-------------------------------
!scalars
 integer :: ii,jj
!arrays
 real(dp) :: work1(3,3),work2(3,3)

! *************************************************************************

 work1(1,1)=frac(1)
 work1(2,2)=frac(2)
 work1(3,3)=frac(3)
 work1(3,2)=frac(4) ; work1(2,3)=frac(4)
 work1(3,1)=frac(5) ; work1(1,3)=frac(5)
 work1(2,1)=frac(6) ; work1(1,2)=frac(6)

 do ii=1,3
   work2(:,ii)=zero
   do jj=1,3
     work2(:,ii)=work2(:,ii)+gprimd(ii,jj)*work1(:,jj)
   end do
 end do

 do ii=1,3
   work1(ii,:)=zero
   do jj=1,3
     work1(ii,:)=work1(ii,:)+gprimd(ii,jj)*work2(jj,:)
   end do
 end do

 cart(1)=work1(1,1)
 cart(2)=work1(2,2)
 cart(3)=work1(3,3)
 cart(4)=work1(2,3)
 cart(5)=work1(1,3)
 cart(6)=work1(1,2)

end subroutine strconv
!!***

!!****f* m_geometry/littlegroup_pert
!!
!! NAME
!! littlegroup_pert
!!
!! FUNCTION
!! If syuse==0 and abs(rfmeth)==2, determines the set of symmetries that leaves a perturbation invariant.
!! (Actually, all symmetries that leaves a q-wavevector invariant should be used to reduce the number
!! of k-points for all perturbations. Unfortunately, one has to take into account the sign reversal of the
!! perturbation under the symmetry operations, which makes GS routines not usable for the respfn code.
!! The intermediate choice was to select only those that keep also the perturbation invariant.
!! Note that the wavevector of the perturbation must also be invariant,
!! a translation vector in real space is NOT allowed ).
!!
!! INPUTS
!! gprimd(3,3)=dimensional primitive translations for reciprocal space (bohr**-1)
!! idir=direction of the perturbation
!! indsym(4,nsym,natom)=indirect indexing of atom labels--see subroutine symatm for definition (if nsym>1)
!! iout=if non-zero, output on unit iout
!! ipert=characteristics of the perturbation
!! natom= number of atoms
!! nsym=number of space group symmetries
!! rfmeth =
!!   1 or -1 if non-stationary block
!!   2 or -2 if stationary block
!!   3 or -3 if third order derivatives
!!   positive if symmetries are used to set elements to zero whenever possible, negative to prevent this to happen.
!! symq(4,2,nsym)= Table computed by littlegroup_q.
!!   three first numbers define the G vector;
!!   fourth number is zero if the q-vector is not preserved, is 1 otherwise
!!   second index is one without time-reversal symmetry, two with time-reversal symmetry
!! symafm(nsym)=(anti)ferromagnetic part of the symmetry operations
!! symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!! symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!! syuse= flag to use the symmetries or not. If 0 usei it, if 1 do not use it.
!! tnons(3,nsym)=nonsymmorphic translations of space group in terms
!!  of real space primitive translations (may be 0)
!! [unit]=By default the routine writes to std_out and this is very annoying if we are inside a big loop.
!!   Use unit=dev_null or a negative integer to disable writing.
!!
!! OUTPUT
!! nsym1 =number of space group symmetries that leaves the perturbation invariant
!! symaf1(nsym1)=(anti)ferromagnetic part of the corresponding symmetry operations
!! symrl1(3,3,nsym1)=corresponding 3x3 matrices of the group symmetries (real space)
!! tnons1(3,nsym1)=corresponding nonsymmorphic translations of space group in terms
!!   of real space primitive translations (may be 0)!!
!!
!! PARENTS
!!      dfpt_looppert,get_npert_rbz,m_dvdb,read_gkk
!!
!! CHILDREN
!!      stresssym,wrtout
!!
!! SOURCE

subroutine littlegroup_pert(gprimd,idir,indsym,iout,ipert,natom,nsym,nsym1, &
&    rfmeth,symafm,symaf1,symq,symrec,symrel,symrl1,syuse,tnons,tnons1, &
&    unit) ! Optional


!Arguments -------------------------------
!scalars
 integer,intent(in) :: idir,iout,ipert,natom,nsym,rfmeth,syuse
 integer,intent(in),optional :: unit
 integer,intent(out) :: nsym1
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symafm(nsym),symq(4,2,nsym)
 integer,intent(in) :: symrec(3,3,nsym),symrel(3,3,nsym)
 integer,intent(out) :: symaf1(nsym),symrl1(3,3,nsym)
 real(dp),intent(in) :: gprimd(3,3),tnons(3,nsym)
 real(dp),intent(out) :: tnons1(3,nsym)

!Local variables -------------------------
!scalars
 integer :: idir1,ii,istr,isym,jj,nsym_test,tok,ount
 character(len=500) :: msg
!arrays
 integer :: sym_test(3,3,2)
 real(dp) :: str_test(6)

! *********************************************************************

 ount = std_out; if (present(unit)) ount = unit

 nsym1=0
 if((ipert==natom+3 .or. ipert==natom+4) .and. syuse==0 .and. abs(rfmeth)==2) then
!  Strain perturbation section
!  Use ground state routine which symmetrizes cartesian stress as a quick
!  and dirty test for the invariance of the strain (ipert,idir) under
!  each candidate symmetry
!  I am presently assuming that translations are acceptable because I dont
!  see why not.

   istr=3*(ipert-natom-3)+idir
   nsym_test=2
!  Store identity as first element for test
   sym_test(:,:,1)=0
   sym_test(1,1,1)=1; sym_test(2,2,1)=1; sym_test(3,3,1)=1
   do isym=1,nsym
     sym_test(:,:,2)=symrec(:,:,isym)
     str_test(:)=0.0_dp
     str_test(istr)=1.0_dp
     call stresssym(gprimd,nsym_test,str_test,sym_test)
     if(abs(str_test(istr)-1.0_dp)<tol8)then
!      The test has been successful !
       nsym1=nsym1+1
       symaf1(nsym1)=symafm(isym)
       do ii=1,3
         tnons1(ii,nsym1)=tnons(ii,isym)
         do jj=1,3
           symrl1(ii,jj,nsym1)=symrel(ii,jj,isym)
         end do
       end do
     end if
   end do

 else if(ipert>natom .or. syuse/=0 .or. abs(rfmeth)/=2)then

!  Not yet coded for d/dk or electric field perturbations
   nsym1=1
   do ii=1,3
     tnons1(ii,1)=0._dp
     symaf1(1)=1
     do jj=1,3
       symrl1(ii,jj,1)=0
       if(ii==jj)symrl1(ii,jj,1)=1
     end do
   end do

 else

   do isym=1,nsym
!    Check that the symmetry operation preserves the wavevector
!    (a translation is NOT allowed)
     if(symq(4,1,isym)==1 .and.&
&     symq(1,1,isym)==0 .and.&
&     symq(2,1,isym)==0 .and.&
&     symq(3,1,isym)==0          )then
!      Check that the symmetry operation preserves the atom
       if(ipert==indsym(4,isym,ipert))then
!        Check if the direction is preserved
         tok=1
         do idir1=1,3
           if((idir1==idir.and.symrec(idir,idir1,isym)/=1) .or.&
&           (idir1/=idir.and.symrec(idir,idir1,isym)/=0))then
             tok=0
           end if
         end do
         if(tok==1)then
!          All the tests have been successful !
           nsym1=nsym1+1
           symaf1(nsym1)=symafm(isym)
           do ii=1,3
             tnons1(ii,nsym1)=tnons(ii,isym)
             do jj=1,3
               symrl1(ii,jj,nsym1)=symrel(ii,jj,isym)
             end do
           end do
         end if

       end if
     end if
   end do
 end if

 if (nsym1<1) then
   write(msg,'(a,i0,a)')' The number of selected symmetries should be > 0, while it is nsym= ',nsym1,'.'
   MSG_BUG(msg)
 end if

 if (nsym1 /= 1) then
   if (iout /= ount .and. iout > 0) then
     write(msg,'(a,i5,a)')' Found ',nsym1,' symmetries that leave the perturbation invariant.'
     call wrtout(iout,msg,'COLL')
   end if
   write(msg,'(a,i5,a)')' littlegroup_pert: found ',nsym1,' symmetries that leave the perturbation invariant: '
   call wrtout(ount,msg,'COLL')
 else
   if (iout /= ount .and. iout > 0) then
     write(msg,'(a,a)')' The set of symmetries contains',' only one element for this perturbation.'
     call wrtout(iout,msg,'COLL')
   end if
   write(msg,'(a)')' littlegroup_pert: only one element in the set of symmetries for this perturbation:'
   call wrtout(ount,msg,'COLL')
 end if

 if (ount > 0) then
   do isym=1,nsym1
     write(msg, '(9i4)' )((symrl1(ii,jj,isym),ii=1,3),jj=1,3)
     call wrtout(ount,msg,'COLL')
   end do
 end if

end subroutine littlegroup_pert
!!***

!!****f* ABINIT/irreducible_set_pert
!! NAME
!! irreducible_set_pert
!!
!! FUNCTION
!! Determines a set of perturbations that form a basis
!! in that, using symmetry, they can be used to generate
!! all other perturbations that are asked to be calculated (target).
!!
!! INPUTS
!!  indsym(4,nsym,natom)=indirect indexing array described above: for each
!!   isym,iatom, fourth element is label of atom into which iatom is sent by
!!   INVERSE of symmetry operation isym; first three elements are the primitive
!!   translations which must be subtracted after the transformation to get back
!!   to the original unit cell.
!!  mpert =maximum number of iper
!!  natom= number of atoms
!!  nsym=number of space group symmetries
!!  rfdir(3)=direction for the perturbations
!!  rfpert(mpert)=information on the perturbations
!!  symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!!  symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!!  symq(4,2,nsym)= (integer) three first numbers define the G vector;
!!   fourth number is 0 if the q-vector is not preserved, is 1 otherwise
!!   second index is one without time-reversal symmetry, two with time-reversal symmetry
!!
!! OUTPUT
!!   pertsy(3,mpert)= the target perturbation is described by the two last indices (idir, and ipert),
!!                    the value is 0, 1 or -1, see notes.
!!
!! NOTES
!! Output will be in the pertsy array,
!!   0 for non-target perturbations
!!   1 for basis perturbations
!!  -1 for perturbations that can be found from basis perturbations
!!
!! PARENTS
!!      get_npert_rbz,m_dvdb,respfn
!!
!! CHILDREN
!!
!! SOURCE

subroutine irreducible_set_pert(indsym,mpert,natom,nsym,pertsy,rfdir,rfpert,symq,symrec,symrel)


!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,nsym
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),rfdir(3),rfpert(mpert)
 integer,intent(in) :: symq(4,2,nsym),symrec(3,3,nsym),symrel(3,3,nsym)
 integer,intent(out) :: pertsy(3,mpert)

!Local variables -------------------------
!scalars
 integer :: found,idir1,idisy1,ii,ipert1,ipesy1,isign,isym,itirev,jj
!arrays
 integer :: sym1(3,3)

! *********************************************************************

!Zero pertsy
 pertsy(:,:)=0

 do ipert1=1,mpert
   do idir1=1,3
     if(rfpert(ipert1)==1.and.rfdir(idir1)==1)then
!      write(std_out,*)' for candidate idir =',idir1,' ipert = ',ipert1

!      Loop on all symmetries, including time-reversal
       do isym=1,nsym
         do itirev=1,2
           isign=3-2*itirev

           if(symq(4,itirev,isym)/=0)then

             found=1

!            Here select the symmetric of ipert1
             if(ipert1<=natom)then
               ipesy1=indsym(4,isym,ipert1)
               do ii=1,3
                 do jj=1,3
                   sym1(ii,jj)=symrec(ii,jj,isym)
                 end do
               end do
             else if(ipert1==(natom+2) .or. ipert1==(natom+6))then
               ipesy1=ipert1
               do ii=1,3
                 do jj=1,3
                   sym1(ii,jj)=symrel(ii,jj,isym)
                 end do
               end do
             else
               found=0
             end if

!            Now that a symmetric perturbation has been obtained,
!            including the expression of the symmetry matrix, see
!            if the symmetric perturbations are available
             if( found==1 ) then

               do idisy1=1,3
                 if(sym1(idir1,idisy1)/=0)then
                   if(pertsy(idisy1,ipesy1)==0)then
                     found=0
                     exit
                   end if
                 end if
               end do
             end if

!            Now, if still found, then it is a symmetric
!            of some linear combination of existing perturbations
             if(found==1)then

!              DEBUG
!              write(std_out,*)' all found !  isym, isign= ',isym,isign
!              write(std_out,1010)((sym1(ii,jj),ii=1,3),jj=1,3)
!              write(std_out,1010)((sym2(ii,jj),ii=1,3),jj=1,3)
!              write(std_out,*)sumr,sumi
!              1010    format(9i4)
!              ENDDEBUG

               pertsy(idir1,ipert1)=-1
               exit ! Exit loop on symmetry operations

             end if

           end if !  End loop on all symmetries + time-reversal
         end do
       end do

!      Now that all symmetries have been examined,
!      if still not symmetric of a linear combination
!      of basis perturbations, then it is a basis perturbation
       if(pertsy(idir1,ipert1)/=-1) pertsy(idir1,ipert1)=1
!      write(std_out,'(a,3i5)' ) ' irreducible_set_pert :',idir1,ipert1,pertsy(idir1,ipert1)

     end if ! End big loop on all elements
   end do
 end do

end subroutine irreducible_set_pert
!!***

end module  m_geometry
!!***
