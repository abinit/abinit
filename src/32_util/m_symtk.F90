!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_symtk
!! NAME
!!  m_symtk
!!
!! FUNCTION
!!  Low-level tools related to symmetries
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR, MG)
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

module m_symtk

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_numeric_tools,  only : isinteger, wrap2_pmhalf

 implicit none

 private
!!***

 public :: mati3inv             ! Invert and transpose orthogonal 3x3 matrix of INTEGER elements.
 public :: mati3det             ! Compute the determinant of a 3x3 matrix of INTEGER elements.
 public :: matr3inv             ! Invert and TRANSPOSE general 3x3 matrix of real*8 elements.
 public :: symdet               ! Compute determinant of each input symmetry matrix sym(3,3,i)
 public :: chkgrp               ! Checks that a set of input symmetries constitutes a group.
 public :: sg_multable          ! Checks that a set of input symmetries constitutes a group.
                                ! TODO: This improved version should replace chkgrp.
 public :: chkorthsy            ! Check the orthogonality of the symmetry operations
 public :: chkprimit            ! Check whether the cell is primitive or not.
 public :: symrelrot            ! Transform symmetry matrices to new new coordinate system.
 public :: littlegroup_q        ! Determines the symmetry operations by which reciprocal vector q is preserved.
!!***

contains
!!***

!!****f* m_symtk/mati3inv
!! NAME
!! mati3inv
!!
!! FUNCTION
!! Invert and transpose orthogonal 3x3 matrix of INTEGER elements.
!!
!! INPUTS
!! mm = integer matrix to be inverted
!!
!! OUTPUT
!! mit = inverse of mm input matrix
!!
!! NOTES
!! Used for symmetry operations.
!! This routine applies to ORTHOGONAL matrices only.
!! Since these form a group, inverses are also integer
!! arrays.  Returned array is TRANSPOSE of inverse, as needed.
!! Note use of integer arithmetic.
!!
!! PARENTS
!!      cg_rotate,chkgrp,classify_bands,debug_tools,dfpt_nstdy,get_full_kgrid
!!      get_npert_rbz,getkgrid,ingeo,m_ab7_symmetry,m_crystal,m_ddb,m_dvdb
!!      m_dynmat,m_fft_mesh,m_fock,m_ptgroups,m_tdep_sym,matpointsym
!!      memory_eval,optic,read_gkk,setsym,strainsym,thmeig,wfconv
!!
!! CHILDREN
!!
!! SOURCE

subroutine mati3inv(mm, mit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mati3inv'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 integer,intent(in) :: mm(3,3)
 integer,intent(out) :: mit(3,3)

!Local variables-------------------------------
!scalars
 integer :: dd
 character(len=500) :: message
!arrays
 integer :: tt(3,3)

! *************************************************************************

 tt(1,1) = mm(2,2) * mm(3,3) - mm(3,2) * mm(2,3)
 tt(2,1) = mm(3,2) * mm(1,3) - mm(1,2) * mm(3,3)
 tt(3,1) = mm(1,2) * mm(2,3) - mm(2,2) * mm(1,3)
 tt(1,2) = mm(3,1) * mm(2,3) - mm(2,1) * mm(3,3)
 tt(2,2) = mm(1,1) * mm(3,3) - mm(3,1) * mm(1,3)
 tt(3,2) = mm(2,1) * mm(1,3) - mm(1,1) * mm(2,3)
 tt(1,3) = mm(2,1) * mm(3,2) - mm(3,1) * mm(2,2)
 tt(2,3) = mm(3,1) * mm(1,2) - mm(1,1) * mm(3,2)
 tt(3,3) = mm(1,1) * mm(2,2) - mm(2,1) * mm(1,2)
 dd  = mm(1,1) * tt(1,1) + mm(2,1) * tt(2,1) + mm(3,1) * tt(3,1)

!Make sure matrix is not singular
 if (dd/=0) then
   mit(:,:)=tt(:,:)/dd
 else
   write(message, '(2a,2x,9i5,a)' )&
&   'Attempting to invert integer array',ch10,mm(:,:),'   ==> determinant is zero.'
   MSG_BUG(message)
 end if

!If matrix is orthogonal, determinant must be 1 or -1
 if (abs(dd)/=1) then
   write(message, '(2a,2x,9i5,a)' )&
&   'Absolute value of determinant should be one',ch10,'but determinant= ',dd
   MSG_BUG(message)
 end if

end subroutine mati3inv
!!***

!!****f* m_symtk/mati3det
!! NAME
!! mati3det
!!
!! FUNCTION
!! Compute the determinant of a 3x3 matrix of INTEGER elements.
!!
!! INPUTS
!! mm = integer matrix
!!
!! OUTPUT
!! det = determinant of the matrix
!!
!! PARENTS
!!      get_kpt_fullbz,getspinrot,m_ab7_symmetry,m_anaddb_dataset,symdet
!!
!! CHILDREN
!!
!! SOURCE

subroutine mati3det(mm, det)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mati3det'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 integer,intent(in) :: mm(3,3)
 integer,intent(out) :: det

! *************************************************************************
 det=mm(1,1)*(mm(2,2) * mm(3,3) - mm(3,2) * mm(2,3)) &
& + mm(2,1)*(mm(3,2) * mm(1,3) - mm(1,2) * mm(3,3)) &
& + mm(3,1)*(mm(1,2) * mm(2,3) - mm(2,2) * mm(1,3))

end subroutine mati3det
!!***

!!****f* m_symtk/matr3inv
!! NAME
!! matr3inv
!!
!! FUNCTION
!! Invert and transpose general 3x3 matrix of real*8 elements.
!!
!! INPUTS
!! aa = 3x3 matrix to be inverted
!!
!! OUTPUT
!! ait = inverse of aa input matrix
!!
!! NOTES
!! Returned array is TRANSPOSE of inverse, as needed to get g from r.
!!
!! PARENTS
!!      berryphase,chkdilatmx,conducti_nc,ddb_hybrid,dfpt_mkvxc,dfpt_mkvxcstr
!!      dfpt_symph,electrooptic,ep_el_weights,ep_fs_weights,ep_ph_weights
!!      fock_getghc,get_kpt_fullbz,getkgrid,getspinrot,gstate,harmonic_thermo
!!      invars2,inwffil,m_cut3d,m_ddb,m_ddk,m_double_grid,m_dynmat
!!      m_effective_potential,m_esymm,m_ewald,m_fock,m_fstab,m_ifc,m_pimd
!!      m_psps,m_strain,m_supercell,m_tdep_latt,make_efg_el,make_efg_ion,metric
!!      mover,optic,outwant,pimd_langevin_npt,prtxf,relaxpol,respfn,smpbz
!!      stresssym,symbrav,symlatt,symmetrize_rprimd,symrelrot,symrhg,tddft
!!      testkgrid,thmeig,uderiv,xcart2xred,xfpack_x2vin
!!
!! CHILDREN
!!
!! SOURCE

subroutine matr3inv(aa, ait)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'matr3inv'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: aa(3,3)
 real(dp),intent(out) :: ait(3,3)

!Local variables-------------------------------
!scalars
 real(dp) :: dd,det,t1,t2,t3
 character(len=500) :: message

! *************************************************************************

 t1 = aa(2,2) * aa(3,3) - aa(3,2) * aa(2,3)
 t2 = aa(3,2) * aa(1,3) - aa(1,2) * aa(3,3)
 t3 = aa(1,2) * aa(2,3) - aa(2,2) * aa(1,3)
 det  = aa(1,1) * t1 + aa(2,1) * t2 + aa(3,1) * t3

!Make sure matrix is not singular
 if (abs(det)>tol16) then
   dd=one/det
 else
   write(message, '(2a,2x,9es16.8,a,a,es16.8,a)' )&
&   'Attempting to invert real(8) 3x3 array',ch10,aa(:,:),ch10,'   ==> determinant=',det,' is zero.'
   MSG_BUG(message)
 end if

 ait(1,1) = t1 * dd
 ait(2,1) = t2 * dd
 ait(3,1) = t3 * dd
 ait(1,2) = (aa(3,1)*aa(2,3)-aa(2,1)*aa(3,3)) * dd
 ait(2,2) = (aa(1,1)*aa(3,3)-aa(3,1)*aa(1,3)) * dd
 ait(3,2) = (aa(2,1)*aa(1,3)-aa(1,1)*aa(2,3)) * dd
 ait(1,3) = (aa(2,1)*aa(3,2)-aa(3,1)*aa(2,2)) * dd
 ait(2,3) = (aa(3,1)*aa(1,2)-aa(1,1)*aa(3,2)) * dd
 ait(3,3) = (aa(1,1)*aa(2,2)-aa(2,1)*aa(1,2)) * dd

end subroutine matr3inv
!!***

!!****f* m_symtk/symdet
!! NAME
!! symdet
!!
!! FUNCTION
!! Compute determinant of each input symmetry matrix sym(3,3,i)
!! and check that the determinant is always +/- 1.  Integer arithmetic.
!!
!! INPUTS
!! nsym=number of symmetry operations
!! sym(3,3,nsym)=integer symmetry array
!!
!! OUTPUT
!! determinant(nsym)=determinant of each symmetry operation
!!
!! PARENTS
!!      remove_inversion,setsym,symptgroup,symspgr
!!
!! CHILDREN
!!      mati3det
!!
!! SOURCE

subroutine symdet(determinant, nsym, sym)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symdet'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
!arrays
 integer,intent(in) :: sym(3,3,nsym)
 integer,intent(out) :: determinant(nsym)

!Local variables-------------------------------
!scalars
 integer :: det,isym
 character(len=500) :: message

! *************************************************************************

 do isym=1,nsym
   call mati3det(sym(:,:,isym),det)
   determinant(isym)=det
   if (abs(det)/=1) then
     write(message,'(a,i5,a,i10,a,a,a,a,a)')&
&     'Abs(determinant) for symmetry number',isym,' is',det,' .',ch10,&
&     'For a legitimate symmetry, abs(determinant) must be 1.',ch10,&
&     'Action: check your symmetry operations (symrel) in input file.'
     MSG_ERROR(message)
   end if
 end do

end subroutine symdet
!!***

!!****f* m_symtk/chkgrp
!! NAME
!! chkgrp
!!
!! FUNCTION
!! Checks that a set of input symmetries constitutes a group.
!!
!! INPUTS
!! nsym = number of symmetry operations
!! symafm = (anti)ferromagnetic part of symmetry operations
!! symrel = 3D matrix containg symmetry operations
!!
!! OUTPUT
!!  ierr=Status error.
!!
!! TODO
!! SHOULD ALSO CHECK THE tnons !
!!
!! PARENTS
!!      chkinp,gensymspgr,m_bz_mesh,m_esymm,m_sigmaph,setsym
!!
!! CHILDREN
!!
!! SOURCE

subroutine chkgrp(nsym,symafm,symrel,ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkgrp'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(out) :: ierr
!arrays
 integer,intent(in) :: symafm(nsym),symrel(3,3,nsym)

!Local variables-------------------------------
!scalars
 integer :: isym,jsym,ksym,symafmchk,testeq=1
 logical :: found_inv
 character(len=500) :: msg
!arrays
 integer :: chk(3,3)

! *************************************************************************

!DEBUG
!write(std_out,*)' chkgrp : enter'
!write(std_out,*)'     isym         symrel            symafm '
!do isym=1,nsym
!write(std_out,'(i3,a,9i3,a,i3)' )isym,'   ',symrel(:,:,isym),'   ',symafm(isym)
!end do
!ENDDEBUG

 ierr = 0

!1) Identity must be the first symmetry.
 if (ANY(symrel(:,:,1) /= identity_3d .or. symafm(1)/=1 )) then
   MSG_WARNING("First operation must be the identity operator")
   ierr = ierr+1
 end if
!
!2) The inverse of each element must belong to the group.
 do isym=1,nsym
   call mati3inv(symrel(:,:,isym),chk)
   chk = TRANSPOSE(chk)
   found_inv = .FALSE.
   do jsym=1,nsym
     if ( ALL(symrel(:,:,jsym) == chk) .and. (symafm(jsym)*symafm(isym) == 1 )) then
       found_inv = .TRUE.; EXIT
     end if
   end do

   if (.not.found_inv) then
     write(msg,'(a,i0,2a)')&
&     "Cannot find the inverse of symmetry operation ",isym,ch10,&
&     "Input symmetries do not form a group "
     MSG_WARNING(msg)
     ierr = ierr+1
   end if

 end do
!
!Closure relation under composition.
 do isym=1,nsym
   do jsym=1,nsym
!
!    Compute the product of the two symmetries
     chk = MATMUL(symrel(:,:,jsym), symrel(:,:,isym))
     symafmchk=symafm(jsym)*symafm(isym)
!
!    Check that product array is one of the original symmetries.
     do ksym=1,nsym
       testeq=1
       if ( ANY(chk/=symrel(:,:,ksym) )) testeq=0
#if 0
!      FIXME this check make v4/t26 and v4/t27 fails.
!      The rotational part is in the group but with different magnetic part!
       if (symafmchk/=symafm(ksym))testeq=0
#endif
       if (testeq==1) exit ! The test is positive
     end do
!
     if(testeq==0) then ! The test is negative
       write(msg, '(a,2i3,a,7a)' )&
&       'product of symmetries',isym,jsym,' is not in group.',ch10,&
&       'This indicates that the input symmetry elements',ch10,&
&       'do not possess closure under group composition.',ch10,&
&       'Action: check symrel, symafm and fix them.'
       MSG_WARNING(msg)
       ierr = ierr+1
     end if

   end do ! jsym
 end do ! isym

end subroutine chkgrp
!!***

!!****f* m_symtk/sg_multable
!! NAME
!! sg_multable
!!
!! FUNCTION
!! Checks that a set of input symmetries constitutes a group.
!!
!! INPUTS
!! nsym=number of symmetry operations
!! symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym)=symmetry operations in real space.
!! tnons(3,nsym)=Fractional translations.
!!
!! OUTPUT
!!  ierr=Status error. A non-zero value signals a failure.
!!  [multable(4,nsym,nsym)]= Optional output.
!!    multable(1,sym1,sym2) gives the index of the symmetry product S1 * S2 in the symrel array. 0 if not found.
!!    multable(2:4,sym1,sym2)= the lattice vector that has to added to the fractional translation
!!      of the operation of index multable(1,sym1,sym2) to obtain the fractional traslation of the product S1 * S2.
!!  [toinv(4,nsym)]= Optional output.
!!    toinv(1,sym1)=Gives the index of the inverse of the symmetry operation.
!!     S1 * S1^{-1} = {E, L} with E identity and L a lattice vector
!!    toinv(2:4,sym1)=The lattice vector L
!!      Note that toinv can be easily obtained from multable but sometimes we do not need the full table.
!!
!! TODO
!!  This improved version should replace chkgrp.
!!
!! PARENTS
!!      m_crystal,m_shirley
!!
!! CHILDREN
!!
!! SOURCE

subroutine sg_multable(nsym,symafm,symrel,tnons,tnons_tol,ierr,multable,toinv)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sg_multable'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(out) :: ierr
 real(dp),intent(in) :: tnons_tol
!arrays
 integer,intent(in) :: symafm(nsym),symrel(3,3,nsym)
 integer,optional,intent(out) :: multable(4,nsym,nsym)
 integer,optional,intent(out) :: toinv(4,nsym)
 real(dp),intent(in) :: tnons(3,nsym)

!Local variables-------------------------------
!scalars
 integer :: sym1,sym2,sym3,prd_symafm
 logical :: found_inv,iseq
 character(len=500) :: msg
!arrays
 integer :: prd_symrel(3,3)
 real(dp) :: prd_tnons(3)

! *************************************************************************

 ierr = 0

!1) Identity must be the first symmetry. Do not check tnons, cell might not be primitive.
 if (ANY(symrel(:,:,1) /= identity_3d .or. symafm(1)/=1 )) then
   MSG_WARNING("First operation must be the identity operator")
   ierr = ierr+1
 end if
!
!2) The inverse of each element must belong to the group.
 do sym1=1,nsym
   found_inv = .FALSE.
   do sym2=1,nsym
     prd_symrel = MATMUL(symrel(:,:,sym1), symrel(:,:,sym2))
     prd_tnons = tnons(:,sym1) + MATMUL(symrel(:,:,sym1),tnons(:,sym2))
     prd_symafm = symafm(sym1)*symafm(sym2)
     if ( ALL(prd_symrel == identity_3d) .and. isinteger(prd_tnons,tnons_tol) .and. prd_symafm == 1 ) then
       found_inv = .TRUE.
       if (PRESENT(toinv)) then
         toinv(1,sym1) = sym2
         toinv(2:4,sym1) = NINT(prd_tnons)
       end if
       EXIT
     end if
   end do

   if (.not.found_inv) then
     write(msg,'(a,i0,2a)')&
&     "Cannot find the inverse of symmetry operation ",sym1,ch10,&
&     "Input symmetries do not form a group "
     MSG_WARNING(msg)
     ierr = ierr+1
   end if
 end do
!
!Check closure relation under composition and construct multiplication table.
 do sym1=1,nsym
   do sym2=1,nsym
!
!    Compute the product of the two symmetries. Convention {A,a} {B,b} = {AB, a+Ab}
     prd_symrel = MATMUL(symrel(:,:,sym1), symrel(:,:,sym2))
     prd_symafm = symafm(sym1)*symafm(sym2)
     prd_tnons = tnons(:,sym1) + MATMUL(symrel(:,:,sym1),tnons(:,sym2))
!
     iseq=.FALSE.
     do sym3=1,nsym ! Check that product array is one of the original symmetries.
       iseq = ( ALL(prd_symrel==symrel(:,:,sym3) )           .and. &
&       isinteger(prd_tnons-tnons(:,sym3),tnons_tol) .and. &
&       prd_symafm==symafm(sym3) )  ! Here v4/t26 and v4/t27 will fail.
!      The rotational part is in the group but with different magnetic part!

       if (iseq) then ! The test is positive
         if (PRESENT(multable)) then
           multable(1,sym1,sym2) = sym3
           multable(2:4,sym1,sym2) = NINT(prd_tnons-tnons(:,sym3))
         end if
         EXIT
       end if
     end do
!
     if (.not.iseq) then ! The test is negative
       write(msg, '(a,2(i0,1x),a,7a)' )&
&       'Product of symmetries:',sym1,sym2,' is not in group.',ch10,&
&       'This indicates that the input symmetry elements',ch10,&
&       'do not possess closure under group composition.',ch10,&
&       'Action: check symrel, symafm and fix them.'
       MSG_WARNING(msg)
       ierr = ierr+1
       if (PRESENT(multable)) then
         multable(1,sym1,sym2) = 0
         multable(2:4,sym1,sym2) = HUGE(0)
       end if
     end if

   end do ! sym2
 end do ! sym1

end subroutine sg_multable
!!***

!!****f* m_symtk/chkorthsy
!! NAME
!! chkorthsy
!!
!! FUNCTION
!! Check the orthogonality of the symmetry operations
!! (lengths and absolute values of scalar products should be preserved)
!!
!! INPUTS
!! gprimd(3,3)=dimensional primitive transl. for reciprocal space (bohr**-1)
!! rmet=Real space metric.
!! nsym=actual number of symmetries
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!! symrel(3,3,1:nsym)=symmetry operations in real space in terms of primitive translations
!!
!! SIDE EFFECTS
!! iexit= if 0 at input, will do the check, and stop if there is a problem, return 0 if no problem
!!        if 1 at input, will always input, return 0 if no problem, -1 if there is a problem,
!!                       also, suppresses printing of problem
!!
!! PARENTS
!!      chkinp,ingeo,symmetrize_rprimd
!!
!! CHILDREN
!!
!! SOURCE

subroutine chkorthsy(gprimd,iexit,nsym,rmet,rprimd,symrel)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkorthsy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(inout) :: iexit
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: gprimd(3,3),rmet(3,3),rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,isym,jj
 real(dp),parameter :: tol=2.0d-8
 real(dp) :: residual,rmet2
 character(len=500) :: message
!arrays
 real(dp) :: prods(3,3),rmet_sym(3,3),rprimd_sym(3,3)

! *************************************************************************

 rmet2=zero
 do ii=1,3
   do jj=1,3
     rmet2=rmet2+rmet(ii,jj)*2
   end do
 end do

!Loop over all symmetry operations
 do isym=1,nsym

!  Compute symmetric of primitive vectors under point symmetry operations
   do ii=1,3
     rprimd_sym(:,ii)=symrel(1,ii,isym)*rprimd(:,1)+&
&     symrel(2,ii,isym)*rprimd(:,2)+&
&     symrel(3,ii,isym)*rprimd(:,3)
   end do

!  If the new lattice is the same as the original one,
!  the lengths and angles are preserved
   do ii=1,3
     rmet_sym(ii,:)=rprimd_sym(1,ii)*rprimd_sym(1,:)+&
&     rprimd_sym(2,ii)*rprimd_sym(2,:)+&
&     rprimd_sym(3,ii)*rprimd_sym(3,:)
   end do

   residual=zero
   do ii=1,3
     do jj=1,3
       residual=residual+(rmet_sym(ii,jj)-rmet(ii,jj))**2
     end do
   end do

   if(sqrt(residual)>tol*sqrt(rmet2))then
     if(iexit==0)then
       write(message, '(a,i5,a,a,a,a,a,es12.4,a,a,a,a,a,a,a)' )&
&       'The symmetry operation number ',isym,' does not preserve',ch10,&
&       'vector lengths and angles.',ch10,&
&       'The value of the residual is ',residual,'.',ch10,&
&       'Action: modify rprim, acell and/or symrel so that',ch10,&
&       'vector lengths and angles are preserved.',ch10,&
&       'Beware, the tolerance on symmetry operations is very small.'
       MSG_ERROR(message)
     else
       iexit=-1
     end if
   end if

!  Also, the scalar product of rprimd_sym and gprimd must give integer numbers
   do ii=1,3
     prods(ii,:)=rprimd_sym(1,ii)*gprimd(1,:)+ &
&     rprimd_sym(2,ii)*gprimd(2,:)+ &
&     rprimd_sym(3,ii)*gprimd(3,:)
   end do

   do ii=1,3
     do jj=1,3
       residual=prods(ii,jj)-anint(prods(ii,jj))
       if(abs(residual)>tol)then
         if(iexit==0)then
           write(message, '(a,i0,a,a,a,a,a,a,a)' )&
&           'The symmetry operation number ',isym,' generates',ch10,&
&           'a different lattice.',ch10,&
&           'Action: modify rprim, acell and/or symrel so that',ch10,&
&           'the lattice is preserved.'
           MSG_ERROR(message)
         else
           iexit=-1
         end if
       end if
     end do
   end do

   if(iexit==-1) exit
 end do ! isym

 if(iexit==1)iexit=0

end subroutine chkorthsy
!!***

!!****f* m_symtk/chkprimit
!! NAME
!! chkprimit
!!
!! FUNCTION
!! Check whether the cell is primitive or not.
!! If chkprim/=0 and the cell is non-primitive, stops.
!!
!! INPUTS
!! chkprim= if non-zero, check that the unit cell is primitive.
!! nsym=actual number of symmetries
!! symafm(nsym)= (anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym)= nsym symmetry operations in real space in terms
!!   of primitive translations
!!
!! OUTPUT
!!  multi=multiplicity of the unit cell
!!
!! PARENTS
!!      symanal
!!
!! CHILDREN
!!
!! SOURCE

subroutine chkprimit(chkprim, multi, nsym, symafm, symrel)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkprimit'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: chkprim,nsym
 integer,intent(out) :: multi
!arrays
 integer,intent(in) :: symafm(nsym),symrel(3,3,nsym)

!Local variables-------------------------------
!scalars
 integer :: isym
 character(len=500) :: message

!**************************************************************************

!Loop over each symmetry operation of the Bravais lattice
!Find whether it is the identity, or a pure translation,
!without change of sign of the spin
 multi=0
 do isym=1,nsym
   if( abs(symrel(1,1,isym)-1)+&
&   abs(symrel(2,2,isym)-1)+&
&   abs(symrel(3,3,isym)-1)+&
&   abs(symrel(1,2,isym))+abs(symrel(2,1,isym))+&
&   abs(symrel(2,3,isym))+abs(symrel(3,2,isym))+&
&   abs(symrel(3,1,isym))+abs(symrel(1,3,isym))+&
&   abs(symafm(isym)-1) == 0 )then
     multi=multi+1
   end if
 end do

!Check whether the cell is primitive
 if(multi>1)then
   if(chkprim/=0)then
     write(message,'(a,a,a,i0,a,a,a,a,a,a,a,a,a)')&
&     'According to the symmetry finder, the unit cell is',ch10,&
&     'NOT primitive. The multiplicity is ',multi,' .',ch10,&
&     'The use of non-primitive unit cells is allowed',ch10,&
&     'only when the input variable chkprim is 0.',ch10,&
&     'Action: either change your unit cell (rprim or angdeg),',ch10,&
&     'or set chkprim to 0.'
     MSG_ERROR(message)
   else
     write(message,'(3a,i0,a,a,a)')&
&     'According to the symmetry finder, the unit cell is',ch10,&
&     'not primitive, with multiplicity= ',multi,'.',ch10,&
&     'This is allowed, as the input variable chkprim is 0.'
     MSG_COMMENT(message)
   end if
 end if

end subroutine chkprimit
!!***

!!****f* m_symtk/symrelrot
!! NAME
!! symrelrot
!!
!! FUNCTION
!! Transform the symmetry matrices symrel expressed in the coordinate system rprimd,
!! to symmetry matrices symrel expressed in the new coordinate system rprimd_new
!!
!! INPUTS
!! nsym=number of symmetries
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!! rprimd_new(3,3)=new dimensional primitive translations for real space (bohr)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! symrel(3,3,nsym)=symmetry operations in real space in terms
!! of primitive translations rprimd at input and rprimd_new at output
!!
!! PARENTS
!!      ingeo,m_esymm,symbrav,symlatt,symspgr
!!
!! CHILDREN
!!      matr3inv
!!
!! SOURCE

subroutine symrelrot(nsym,rprimd,rprimd_new,symrel,tolsym)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symrelrot'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 real(dp),intent(in) :: tolsym
!arrays
 integer,intent(inout) :: symrel(3,3,nsym)
 real(dp),intent(in) :: rprimd(3,3),rprimd_new(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,isym,jj
 real(dp) :: val
 character(len=500) :: message
!arrays
 real(dp) :: coord(3,3),coordinvt(3,3),matr1(3,3),matr2(3,3),rprimd_invt(3,3)

!**************************************************************************

!Compute the coordinates of rprimd_new in the system defined by rprimd(:,:)
 call matr3inv(rprimd,rprimd_invt)
 do ii=1,3
   coord(:,ii)=rprimd_new(1,ii)*rprimd_invt(1,:)+ &
&   rprimd_new(2,ii)*rprimd_invt(2,:)+ &
&   rprimd_new(3,ii)*rprimd_invt(3,:)
 end do

!Transform symmetry matrices in the system defined by rprimd_new
 call matr3inv(coord,coordinvt)
 do isym=1,nsym
   do ii=1,3
     matr1(:,ii)=symrel(:,1,isym)*coord(1,ii)+&
&     symrel(:,2,isym)*coord(2,ii)+&
&     symrel(:,3,isym)*coord(3,ii)
   end do
   do ii=1,3
     matr2(:,ii)=coordinvt(1,:)*matr1(1,ii)+&
&     coordinvt(2,:)*matr1(2,ii)+&
&     coordinvt(3,:)*matr1(3,ii)
   end do

!  Check that the new symmetry matrices are made of integers, and store them
   do ii=1,3
     do jj=1,3
       val=matr2(ii,jj)
!      Need to allow for twice tolsym, in case of centered Bravais lattices (but do it for all lattices ...)
       if(abs(val-nint(val))>two*tolsym)then
         write(message,'(2a,a,i3,a,a,3es14.6,a,a,3es14.6,a,a,3es14.6)')&
&         'One of the components of symrel is non-integer,',ch10,&
&         '  for isym=',isym,ch10,&
&         '  symrel=',matr2(:,1),ch10,&
&         '         ',matr2(:,2),ch10,&
&         '         ',matr2(:,3)
         MSG_ERROR_CLASS(message, "TolSymError")
       end if
       symrel(ii,jj,isym)=nint(val)
     end do
   end do
!  End loop isym
 end do

end subroutine symrelrot
!!***

!!****f* m_symtk/littlegroup_q
!! NAME
!! littlegroup_q
!!
!! FUNCTION
!! Determines the symmetry operations by which reciprocal vector q is preserved,
!! modulo a primitive reciprocal lattice vector, and the time-reversal symmetry.
!!
!! INPUTS
!! nsym=number of space group symmetries
!! qpt(3)= vector in reciprocal space
!! symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!! [prtvol]=integer flag defining the verbosity of output. =0 if no output is provided.
!! prtgkk= integer flag. If 1 provide output of electron-phonon "gkk" matrix elements, for further
!!     treatment by mrggkk utility or anaddb utility. If 0 no output is provided.
!!
!! OUTPUT
!! symq(4,2,nsym)= three first numbers define the G vector;
!!     fourth number is zero if the q-vector is not preserved, is 1 otherwise
!!     second index is one without time-reversal symmetry, two with time-reversal symmetry
!! timrev=1 if the time-reversal symmetry preserves the wavevector,
!!   modulo a reciprocal lattice vector (in principle, see below).
!!
!! NOTES
!! The condition is:
!!
!!    $q =  O  S(q) - G$
!!
!! with O being either the identity or the time reversal symmetry (= inversion in reciprocal space)
!! and G being a primitive vector of the reciprocal lattice.
!! If the time-reversal (alone) also preserves q, modulo a lattice vector, then timrev is set to 1, otherwise 0.
!!
!! TODO
!! timrev is put to 1 only for Gamma.
!! Better handling should be provided in further version.
!!
!! PARENTS
!!      get_npert_rbz,m_bz_mesh,m_ddb,m_dvdb,m_dynmat,m_esymm,m_gkk,m_phgamma
!!      m_sigmaph,memory_eval,read_gkk,respfn
!!
!! CHILDREN
!!      wrap2_pmhalf,wrtout
!!
!! SOURCE

subroutine littlegroup_q(nsym,qpt,symq,symrec,symafm,timrev,prtvol,use_sym)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'littlegroup_q'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(in),optional :: prtvol,use_sym
 integer,intent(out) :: timrev
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 integer,intent(in) :: symafm(nsym)
 integer,intent(out) :: symq(4,2,nsym)
 real(dp),intent(in) :: qpt(3)

!Local variables -------------------------
!scalars
 integer :: ii,isign,isym,itirev,my_prtvol
 real(dp),parameter :: tol=2.d-8
 real(dp) :: reduce
 character(len=500) :: message
!arrays
 real(dp) :: difq(3),qsym(3),shift(3)

! *********************************************************************

 my_prtvol=0 ; if (PRESENT(prtvol)) my_prtvol=prtvol

! Initialise the array symq
 symq = 0

 isym = symafm(1) ! just to fool abirules and use symafm for the moment

 do isym=1,nsym
!   if (symafm(isym) /= 1) cycle ! skip afm symops
! TODO: check how much of the afm syms are coded in the rf part of the code. cf
! test v3 / 12
   do itirev=1,2
     isign=3-2*itirev  ! isign is 1 without time-reversal, -1 with time-reversal

!    Get the symmetric of the vector
     do ii=1,3
       qsym(ii)=qpt(1)*isign*symrec(ii,1,isym)&
&       +qpt(2)*isign*symrec(ii,2,isym)&
&       +qpt(3)*isign*symrec(ii,3,isym)
     end do

!    Get the difference between the symmetric and the original vector

     symq(4,itirev,isym)=1
     do ii=1,3
       difq(ii)=qsym(ii)-qpt(ii)
!      Project modulo 1 in the interval ]-1/2,1/2] such that difq = reduce + shift
       call wrap2_pmhalf(difq(ii),reduce,shift(ii))
       if(abs(reduce)>tol)symq(4,itirev,isym)=0
     end do

!    SP: When prtgkk is asked (GKK matrix element will be output), one has to
!    disable symmetries. There is otherwise a jauge problem with the unperturbed
!    and the perturbed wavefunctions. This leads to a +- 5% increase in computational
!    cost but provide the correct GKKs (i.e. the same as without the use of
!    symmerties.)

     if (PRESENT(use_sym)) then
       if (use_sym == 0) then
         symq(4,itirev,isym)=0
         symq(4,itirev,1)=1
       end if
     end if

!    If the operation succeded, change shift from real(dp) to integer, then exit loop
     if(symq(4,itirev,isym)/=0)then
       if (my_prtvol>0) then
         if(itirev==1)write(message,'(a,i4,a)')' littlegroup_q : found symmetry',isym,' preserves q '
         if(itirev==2)write(message,'(a,i4,a)')' littlegroup_q : found symmetry ',isym,' + TimeReversal preserves q '
         call wrtout(std_out,message,'COLL')
       end if
!      Uses the mathematical function NINT = nearest integer
       do ii=1,3
         symq(ii,itirev,isym)=nint(shift(ii))
       end do
     end if

   end do !itirev
 end do !isym

!Test time-reversal symmetry
 timrev=1
 do ii=1,3
!  Unfortunately, this version does not work yet ...
!  call wrap2_pmhalf(2*qpt(ii),reduce,shift(ii))
!  if(abs(reduce)>tol)timrev=0
!  So, this is left ...
   if(abs(qpt(ii))>tol)timrev=0
 end do

 if(timrev==1.and.my_prtvol>0)then
   write(message, '(a,a,a)' )&
&   ' littlegroup_q : able to use time-reversal symmetry. ',ch10,&
&   '  (except for gamma, not yet able to use time-reversal symmetry)'
   call wrtout(std_out,message,'COLL')
 end if

end subroutine littlegroup_q
!!***

end module m_symtk
!!***
