!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_geometry
!! NAME
!!  m_geometry
!!
!! FUNCTION
!!  This module contains basic tools to operate on vectors expressed in reduced coordinates.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2018 ABINIT group (MG, MT, FJ, TRangel)
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
!! API slightly changed the wrt wannier90 version.
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

end module  m_geometry
!!***
