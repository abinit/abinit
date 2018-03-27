!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_geometry
!! NAME
!!  m_geometry
!!
!! FUNCTION
!!  This module contains basic tools to operate on vectors expressed in reduced coordinates.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2018 ABINIT group (MG, MT, FJ, TRangel, DCA, XG, AHR)
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
 use m_profiling_abi
 use m_errors
 use m_atomdata
 use m_sort

 use m_io_tools,       only : open_file
 use m_numeric_tools,  only : uniformrandom
 use m_pptools,        only : prmat

 implicit none

 private

 public :: normv              ! Norm of vector(s) in reduced coordinates either in real or reciprocal space.
 public :: vdotw              ! Scalar product between two reduced vectors either in real or reciprocal space.
 public :: acrossb            ! Cross product of two 3-vectors.
 public :: wigner_seitz       ! Find the grid of points falling inside the Wigner-Seitz cell.
 public :: phdispl_cart2red   ! Calculate the displacement vectors for all branches in reduced coordinates.
 public :: spinrot_cmat       ! Construct 2x2 complex matrix representing rotation operator in spin-space.
 public :: rotmat             ! Finds the rotation matrix.
 public :: fixsym             ! Check that iatfix does not break symmetry.
 public :: mkradim            ! Make rprim and acell from rprimd
 public :: mkrdim             ! Make rprimd from acell from rprim
 public :: xcart2xred         ! From cart coords to reduced
 public :: xred2xcart         ! From reduced coords to cart.
 public :: fred2fcart         ! Convert reduced forces into cartesian forces
 public :: fcart2fred         ! Convert cartesian forces into reduced forces
 public :: bonds_lgth_angles  ! Write GEO file
 public :: randomcellpos      ! Creates unit cell with random atomic positions.
 public :: ioniondist         ! Compute ion-ion distances
 public :: shellstruct        ! Calculates shell structure (multiplicities, radii)

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'normv_rdp_vector'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'normv_int_vector'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'normv_int_vector_array'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'normv_rdp_vector_array'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vdotw_rr_vector'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vdotw_rc_vector'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'acrossb'
!End of the abilint section

 implicit none

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
!!  In input irvec and ndegen are NULL pointers. They are allocated with the correct
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

subroutine wigner_seitz(center,lmax,kptrlatt,rmet,npts,irvec,ndegen,prtvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wigner_seitz'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol
 integer,intent(out) :: npts
!arrays
 integer,intent(in) :: kptrlatt(3,3),lmax(3)
 integer,pointer :: irvec(:,:),ndegen(:)
 real(dp),intent(in) :: center(3),rmet(3,3)

!Local variables-------------------------------
!scalars
 integer :: in1,in2,in3,l1,l2,l3,ii,icount,n1,n2,n3
 integer :: l0,l1_max,l2_max,l3_max,nl,verbose,mm1,mm2,mm3
 real(dp) :: tot,dist_min
 real(dp),parameter :: TOL_DIST=tol6
 character(len=500) :: msg
!arrays
 real(dp) :: diff(3)
 real(dp),allocatable :: dist(:)
 real(dp),allocatable :: swap2(:,:),swap1(:)

! *************************************************************************

 verbose=0; if (PRESENT(prtvol)) verbose=prtvol

 if (kptrlatt(1,2)/=0 .or. kptrlatt(2,1)/=0 .or. &
&    kptrlatt(1,3)/=0 .or. kptrlatt(3,1)/=0 .or. &
&    kptrlatt(2,3)/=0 .or. kptrlatt(3,2)/=0 ) then
   MSG_ERROR('Off-diagonal elements of kptrlatt must be zero')
 end if

 n1=kptrlatt(1,1)
 n2=kptrlatt(2,2)
 n3=kptrlatt(3,3)

 l1_max=lmax(1)
 l2_max=lmax(2)
 l3_max=lmax(3)

 nl=(2*l1_max+1)*(2*l2_max+1)*(2*l3_max+1)
 l0=1+l1_max*(1+(2*l2_max+1)**2+(2*l3_max+1)) ! Index of the origin.
 ABI_MALLOC(dist,(nl))

 ! Allocate with maximum size
 mm1=2*n1+1
 mm2=2*n2+1
 mm3=2*n3+1
 ABI_MALLOC(irvec,(3,mm1*mm2*mm3))
 ABI_MALLOC(ndegen,(mm1*mm2*mm3))

 npts=0
 do in1=-n1,n1
   do in2=-n2,n2
     do in3=-n3,n3
      !
      ! Loop over the nl points R. R=0 corresponds to l1=l2=l3=1, or icount=l0
      icount=0
      do l1=-l1_max,l1_max
        do l2=-l2_max,l2_max
          do l3=-l3_max,l3_max
            ! * Calculate |r-R-r_0|^2.
            diff(1)= in1 -l1*n1 -center(1)
            diff(2)= in2 -l2*n2 -center(2)
            diff(3)= in3 -l3*n3 -center(3)
            icount=icount+1
            dist(icount)=DOT_PRODUCT(diff,MATMUL(rmet,diff))
          end do
        end do
      end do

      dist_min=MINVAL(dist)

      if (ABS(dist(l0)-dist_min)<TOL_DIST) then
        npts=npts+1
        ndegen(npts)=0
        do ii=1,nl
          if (ABS(dist(ii)-dist_min)<TOL_DIST) ndegen(npts)=ndegen(npts)+1
        end do
        irvec(1,npts)=in1
        irvec(2,npts)=in2
        irvec(3,npts)=in3
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

 ! === Check the "sum rule" ===
 tot=zero
 do ii=1,npts
   tot=tot+one/ndegen(ii)
 end do
 if (ABS(tot-(n1*n2*n3))>tol8) then
   write(msg,'(a,es16.8,a,i5)')'Something wrong in the generation of the mesh ',tot,' /= ',n1*n2*n3
   MSG_ERROR(msg)
 end if

 ABI_FREE(dist)

 ! === Reallocate the output with correct size ===
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phdispl_cart2red'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spinrot_cmat'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rotmat'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

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
&   'Action : If you did not intend this, make |z|<10 and |x|<10 ',ch10
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fixsym'
!End of the abilint section

 implicit none

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
         write(message, '(a,i6,a,a,i6,a,a,a,a,a,a,a)' )&
&         '  Atom number ',jatom,' is symmetrically',&
&         ' equivalent to atom number ',iatom,',',ch10,&
&         '  but according to iatfix, iatfixx, iatfixy and iatfixz, they',ch10,&
&         '  are not fixed along the same directions, which is forbidden.',ch10,&
&         '  Action : modify either the symmetry or iatfix(x,y,z) and resubmit.'
         MSG_ERROR(message)
       end if
     end do
   end do
 end if

end subroutine fixsym
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkradim'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkrdim'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcart2xred'
 use interfaces_32_util
!End of the abilint section

 implicit none

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
     xred(mu,iatom)= gprimd(1,mu)*xcart(1,iatom)+gprimd(2,mu)*xcart(2,iatom)+&
&     gprimd(3,mu)*xcart(3,iatom)
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xred2xcart'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fred2fcart'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fcart2fred'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bonds_lgth_angles'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

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
&     'Action : decrease natom, or contact ABINIT group.'
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rsdot'
!End of the abilint section

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'randomcellpos'
!End of the abilint section

 implicit none

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
&   'Action : switch to another value of random_atpos.'
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

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'shellstruct'
 use interfaces_14_hidewrite
 use interfaces_41_geometry, except_this_one => shellstruct
!End of the abilint section

 implicit none

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

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ioniondist'
 use interfaces_14_hidewrite
 use interfaces_41_geometry, except_this_one => ioniondist
!End of the abilint section

 implicit none

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

end module  m_geometry
!!***
