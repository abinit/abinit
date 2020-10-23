!!****m* ABINIT/m_symtk
!! NAME
!!  m_symtk
!!
!! FUNCTION
!!  Low-level tools related to symmetries
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (RC, XG, GMR, MG, JWZ)
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
 use m_abicore

 use m_numeric_tools,  only : isinteger, wrap2_pmhalf
 use m_hide_lapack,    only : matrginv

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
 public :: matpointsym          ! Symmetrizes a 3x3 input matrix using the point symmetry of the input atom
 public :: holocell             ! Examine whether the trial conventional cell described by cell_base
                                ! is coherent with the required holohedral group.
 public :: symmetrize_xred      ! Symmetrize atomic coordinates using input symmetry matrices symrel
 public :: symmetrize_rprimd    ! Generate new rprimd on the basis of the expected characteristics of the conventional cell
 public :: symchk               ! Symmetry checker for atomic coordinates.
 public :: symatm               ! Build indsym table describing the action of the symmetry operations on the atomic positions.
 public :: symcharac            ! Get the type of axis for the symmetry.
 public :: smallprim            ! Find the smallest possible primitive vectors for an input lattice
 public :: print_symmetries     ! Helper function to print symmetries in a nice format.
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
!!      m_ab7_symmetry,m_cgtk,m_classify_bands,m_crystal,m_ddb,m_dfpt_scfcv
!!      m_dtset,m_dvdb,m_dynmat,m_fft_mesh,m_fock,m_geometry,m_ingeo,m_inwffil
!!      m_iogkk,m_kpts,m_memeval,m_mover_effpot,m_ptgroups,m_spacepar,m_symtk
!!      m_tdep_sym,m_thmeig,optic
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine mati3inv(mm, mit)

!Arguments ------------------------------------
!arrays
 integer,intent(in) :: mm(3,3)
 integer,intent(out) :: mit(3,3)

!Local variables-------------------------------
!scalars
 integer :: dd
 character(len=500) :: msg
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
   write(msg, '(2a,2x,9i5,a)' )'Attempting to invert integer array',ch10,mm(:,:),'   ==> determinant is zero.'
   MSG_BUG(msg)
 end if

!If matrix is orthogonal, determinant must be 1 or -1
 if (abs(dd)/=1) then
   write(msg, '(2a,2x,9i5,a)' )'Absolute value of determinant should be one',ch10,'but determinant= ',dd
   MSG_BUG(msg)
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
!!      m_ab7_symmetry,m_anaddb_dataset,m_geometry,m_kpts,m_symtk
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine mati3det(mm, det)

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
!!      m_berryphase,m_conducti,m_cut3d,m_ddb,m_dfpt_mkvxc,m_dfpt_mkvxcstr
!!      m_double_grid,m_dvdb,m_dynmat,m_effective_potential,m_elpolariz
!!      m_epweights,m_esymm,m_ewald,m_fock,m_fock_getghc,m_fstab,m_geometry
!!      m_gstate,m_harmonic_thermo,m_ifc,m_inwffil,m_kpts,m_longwave,m_mover
!!      m_mover_effpot,m_nucprop,m_outwant,m_phonons,m_pimd,m_pimd_langevin
!!      m_psps,m_raman,m_relaxpol,m_respfn_driver,m_scup_dataset,m_spacepar
!!      m_strain,m_supercell,m_supercell_maker,m_symfind,m_symtk,m_tddft
!!      m_tdep_latt,m_thmeig,m_xfpack,optic
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine matr3inv(aa, ait)

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: aa(3,3)
 real(dp),intent(out) :: ait(3,3)

!Local variables-------------------------------
!scalars
 real(dp) :: dd,det,t1,t2,t3
 character(len=500) :: msg

! *************************************************************************

 t1 = aa(2,2) * aa(3,3) - aa(3,2) * aa(2,3)
 t2 = aa(3,2) * aa(1,3) - aa(1,2) * aa(3,3)
 t3 = aa(1,2) * aa(2,3) - aa(2,2) * aa(1,3)
 det  = aa(1,1) * t1 + aa(2,1) * t2 + aa(3,1) * t3

!Make sure matrix is not singular
 if (abs(det)>tol16) then
   dd=one/det
 else
   write(msg, '(2a,2x,9es16.8,a,a,es16.8,a)' )&
     'Attempting to invert real(8) 3x3 array',ch10,aa(:,:),ch10,'   ==> determinant=',det,' is zero.'
   MSG_BUG(msg)
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
!!      m_geometry,m_spacepar,m_spgdata,m_symfind
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine symdet(determinant, nsym, sym)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
!arrays
 integer,intent(in) :: sym(3,3,nsym)
 integer,intent(out) :: determinant(nsym)

!Local variables-------------------------------
!scalars
 integer :: det,isym
 character(len=500) :: msg

! *************************************************************************

 do isym=1,nsym
   call mati3det(sym(:,:,isym),det)
   determinant(isym)=det
   if (abs(det)/=1) then
     write(msg,'(a,i0,a,i0,a,a,a,a,a)')&
      'Abs(determinant) for symmetry number ',isym,' is ',det,' .',ch10,&
      'For a legitimate symmetry, abs(determinant) must be 1.',ch10,&
      'Action: check your symmetry operations (symrel) in input file.'
     MSG_ERROR(msg)
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
!!      m_bz_mesh,m_chkinp,m_esymm,m_lgroup,m_spacepar,m_spgbuilder
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine chkgrp(nsym, symafm, symrel, ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(out) :: ierr
!arrays
 integer,intent(in) :: symafm(nsym),symrel(3,3,nsym)

!Local variables-------------------------------
!scalars
 integer :: isym,jsym,ksym,print_warning=1,symafmchk,testeq=1
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
 ! 2) The inverse of each element must belong to the group.
 do isym=1,nsym
   call mati3inv(symrel(:,:,isym), chk)
   chk = TRANSPOSE(chk)
   found_inv = .FALSE.
   do jsym=1,nsym
     if ( ALL(symrel(:,:,jsym) == chk) .and. (symafm(jsym) * symafm(isym) == 1 )) then
       found_inv = .TRUE.; EXIT
     end if
   end do

   if (.not.found_inv) then
     write(msg,'(a,i0,2a)')&
      "Cannot find the inverse of symmetry operation ",isym,ch10,&
      "Input symmetries do not form a group "
     MSG_WARNING(msg)
     ierr = ierr+1
   end if

 end do
 !
 ! Closure relation under composition.
 do isym=1,nsym
   do jsym=1,nsym
     !
     ! Compute the product of the two symmetries
     chk = MATMUL(symrel(:,:,jsym), symrel(:,:,isym))
     symafmchk=symafm(jsym)*symafm(isym)
     !
     ! Check that product array is one of the original symmetries.
     do ksym=1,nsym
       testeq=1
       if ( ANY(chk/=symrel(:,:,ksym) )) testeq=0
#if 0
!      FIXME this check make v4/t26 and v4/t27 fails.
!      The rotational part is in the group but with different magnetic part!
       if (symafmchk /= symafm(ksym)) testeq=0
#endif
       if (testeq==1) exit ! The test is positive
     end do
!
     if(testeq==0 .and. print_warning==1) then
       ! The test is negative
       write(msg, '(a,2i3,a,9a)' )&
        'Product of symmetries',isym,jsym,' is not in group.',ch10,&
        'This indicates that the input symmetry elements',ch10,&
        'do not possess closure under group composition.',ch10,&
        'ABINIT might stop with an ERROR after trying to correct and making a few more checks.',ch10,&
        'Action: check symrel, symafm and possibly atomic positions, and fix them.'
       MSG_WARNING(msg)
       ierr = ierr+1
       print_warning=0
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
!!     S1 * S1^{-1} = {E, L} with E the identity and L a real-space lattice vector.
!!    toinv(2:4,sym1)=The lattice vector L
!!      Note that toinv can be easily obtained from multable but sometimes we do not need the full table.
!!
!! TODO
!!  This improved version should replace chkgrp.
!!
!! PARENTS
!!      m_crystal
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine sg_multable(nsym, symafm, symrel, tnons, tnons_tol, ierr, multable, toinv)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(out) :: ierr
 real(dp),intent(in) :: tnons_tol
!arrays
 integer,intent(in) :: symafm(nsym),symrel(3,3,nsym)
 integer,optional,intent(out) :: multable(4,nsym,nsym), toinv(4,nsym)
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

 ! 1) Identity must be the first symmetry. Do not check if tnon == 0 as cell might not be primitive.
 if (any(symrel(:,:,1) /= identity_3d .or. symafm(1) /= 1)) then
   MSG_WARNING("First operation must be the identity operator")
   ierr = ierr + 1
 end if

 ! 2) The inverse of each element must belong to the group.
 do sym1=1,nsym
   found_inv = .FALSE.
   do sym2=1,nsym
     prd_symrel = matmul(symrel(:,:,sym1), symrel(:,:,sym2))
     prd_tnons = tnons(:,sym1) + matmul(symrel(:,:,sym1), tnons(:,sym2))
     prd_symafm = symafm(sym1)*symafm(sym2)
     if ( all(prd_symrel == identity_3d) .and. isinteger(prd_tnons, tnons_tol) .and. prd_symafm == 1 ) then
       found_inv = .TRUE.
       if (present(toinv)) then
         toinv(1, sym1) = sym2
         toinv(2:4, sym1) = nint(prd_tnons)
       end if
       exit
     end if
   end do

   if (.not. found_inv) then
     write(msg,'(a,i0,2a)')&
      "Cannot find the inverse of symmetry operation ",sym1,ch10,&
      "Input symmetries do not form a group "
     MSG_WARNING(msg)
     ierr = ierr + 1
   end if
 end do

 ! Check closure relation under composition and construct multiplication table.
 do sym1=1,nsym
   do sym2=1,nsym

     ! Compute the product of the two symmetries. Convention {A,a} {B,b} = {AB, a+Ab}
     prd_symrel = matmul(symrel(:,:,sym1), symrel(:,:,sym2))
     prd_symafm = symafm(sym1) * symafm(sym2)
     prd_tnons = tnons(:, sym1) + matmul(symrel(:,:,sym1), tnons(:,sym2))

     ! Check that product array is one of the original symmetries.
     iseq = .False.
     do sym3=1,nsym
       iseq = (all(prd_symrel == symrel(:,:,sym3) ) .and. &
               isinteger(prd_tnons - tnons(:,sym3), tnons_tol) .and. &
               prd_symafm == symafm(sym3) )  ! Here v4/t26 and v4/t27 will fail.

       ! The rotational part is in the group but with different magnetic part!
       if (iseq) then
         ! The test is positive
         if (present(multable)) then
           multable(1,sym1,sym2) = sym3
           multable(2:4,sym1,sym2) = nint(prd_tnons - tnons(:,sym3))
         end if
         exit
       end if
     end do

     if (.not. iseq) then
       ! The test is negative
       write(msg, '(a,2(i0,1x),a,7a)' )&
         'Product of symmetries:',sym1,sym2,' is not in group.',ch10,&
         'This indicates that the input symmetry elements',ch10,&
         'do not possess closure under group composition.',ch10,&
         'Action: check symrel, symafm and fix them.'
       MSG_WARNING(msg)
       ierr = ierr + 1
       if (present(multable)) then
         multable(1, sym1, sym2) = 0
         multable(2:4, sym1, sym2) = huge(0)
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
!! tolsym=defines the tolerance on the orthogonality, after multiplication by 2.
!!
!! SIDE EFFECTS
!! iexit= if 0 at input, will do the check, and stop if there is a problem, return 0 if no problem
!!        if 1 at input, will always output, return 0 if no problem, -1 if there is a problem,
!!                       also, suppresses printing of problem
!!
!! PARENTS
!!      m_chkinp,m_ingeo,m_symtk
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine chkorthsy(gprimd,iexit,nsym,rmet,rprimd,symrel,tolsym)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(inout) :: iexit
 real(dp),intent(in) :: tolsym
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: gprimd(3,3),rmet(3,3),rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,isym,jj
 real(dp) :: residual,rmet2
 character(len=500) :: msg
!arrays
 real(dp) :: prods(3,3),rmet_sym(3,3),rprimd_sym(3,3)

! *************************************************************************

!DEBUG
!write(std_out,'(a)') ' chkorthsy : enter '
!write(std_out,'(a,i3)') ' nsym=',nsym
!do isym=1,nsym
!  write(std_out,'(9i4)')symrel(:,:,isym)
!enddo
!write(std_out, '(a)') ' Matrix rprimd :'
!do ii=1,3
!  write(std_out, '(3es16.8)')rprimd(:,ii)
!enddo
!write(std_out, '(a)') ' Matrix rmet :'
!do ii=1,3
!  write(std_out, '(3es16.8)')rmet(:,ii)
!enddo
!ENDDEBUG

 rmet2=zero
 do ii=1,3
   do jj=1,3
     rmet2=rmet2+rmet(ii,jj)*2
   end do
 end do

!Loop over all symmetry operations
 do isym=1,nsym

!DEBUG
!write(std_out,'(a,a,i4)') ch10,' Check for isym=',isym
!ENDDEBUG

!  Compute symmetric of primitive vectors under point symmetry operations
   do ii=1,3
     rprimd_sym(:,ii)=symrel(1,ii,isym)*rprimd(:,1)+&
                      symrel(2,ii,isym)*rprimd(:,2)+&
                      symrel(3,ii,isym)*rprimd(:,3)
   end do

!  If the new lattice is the same as the original one,
!  the lengths and angles are preserved
   do ii=1,3
     rmet_sym(ii,:)=rprimd_sym(1,ii)*rprimd_sym(1,:)+&
                    rprimd_sym(2,ii)*rprimd_sym(2,:)+&
                    rprimd_sym(3,ii)*rprimd_sym(3,:)
   end do

   residual=zero
   do ii=1,3
     do jj=1,3
       residual=residual+(rmet_sym(ii,jj)-rmet(ii,jj))**2
     end do
   end do

   if(sqrt(residual) > two*tolsym*sqrt(rmet2))then
     if(iexit==0)then
       write(std_out, '(a)') ' Matrix rprimd :'
       do ii=1,3
         write(std_out, '(3es16.8)')rprimd(:,ii)
       enddo
       write(std_out, '(a)') ' Matrix rmet :'
       do ii=1,3
         write(std_out, '(3es16.8)')rmet(:,ii)
       enddo
       write(std_out, '(a)') ' Matrix rprimd_sym :'
       do ii=1,3
         write(std_out, '(3es16.8)')rprimd_sym(:,ii)
       enddo
       write(std_out, '(a)') ' Matrix rmet_sym :'
       do ii=1,3
         write(std_out, '(3es16.8)')rmet_sym(:,ii)
       enddo
       write(std_out, '(a)') ' Matrix rmet_sym-rmet :'
       do ii=1,3
         write(std_out, '(3es16.8)')(rmet_sym(:,ii)-rmet(:,ii))
       enddo
       write(msg, '(a,i0,5a,es12.4,a,es12.4,6a)' )&
        'The symmetry operation number ',isym,' does not preserve',ch10,&
        'vector lengths and angles.',ch10,&
        'The value of the square root of residual is: ',sqrt(residual),&
&       '  that is greater than threshold:', two*tolsym*sqrt(rmet2),ch10,&
        'Action: modify rprim, acell and/or symrel so that',ch10,&
        'vector lengths and angles are preserved.',ch10,&
        'Beware, the tolerance on symmetry operations is very small.'
       MSG_ERROR(msg)
     else
       iexit=-1
     end if
   end if

!  Also, the scalar product of rprimd_sym and gprimd must give integer numbers
   do ii=1,3
     prods(ii,:)=rprimd_sym(1,ii)*gprimd(1,:)+ &
                 rprimd_sym(2,ii)*gprimd(2,:)+ &
                 rprimd_sym(3,ii)*gprimd(3,:)
   end do

   do ii=1,3
     do jj=1,3
       residual=prods(ii,jj)-anint(prods(ii,jj))
       if(abs(residual)>two*tolsym)then
         if(iexit==0)then
           write(msg, '(a,i0,5a,es12.4,a,es12.4,4a)' )&
            'The symmetry operation number ',isym,' generates',ch10,&
            'a different lattice.',ch10,&
            'The value of the residual is: ',residual, 'that is greater than the threshold:', two*tolsym, ch10,&
            'Action: modify rprim, acell and/or symrel so that',ch10,&
            'the lattice is preserved.'
           MSG_ERROR(msg)
         else
           iexit=-1
         end if
       end if
     end do
   end do

   if(iexit==-1) exit
 end do ! isym

 if(iexit==1)iexit=0

!DEBUG
!write(std_out,'(a)') ' chkorthsy : exit '
!ENDDEBUG

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
!!      m_symfind
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine chkprimit(chkprim, multi, nsym, symafm, symrel)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: chkprim,nsym
 integer,intent(out) :: multi
!arrays
 integer,intent(in) :: symafm(nsym),symrel(3,3,nsym)

!Local variables-------------------------------
!scalars
 integer :: isym
 character(len=500) :: msg

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
     write(msg,'(a,a,a,i0,a,a,a,a,a,a,a,a,a)')&
     'According to the symmetry finder, the unit cell is',ch10,&
     'NOT primitive. The multiplicity is ',multi,' .',ch10,&
     'The use of non-primitive unit cells is allowed',ch10,&
     'only when the input variable chkprim is 0.',ch10,&
     'Action: either change your unit cell (rprim or angdeg),',ch10,&
     'or set chkprim to 0.'
     MSG_ERROR(msg)
   else
     write(msg,'(3a,i0,a,a,a)')&
      'According to the symmetry finder, the unit cell is',ch10,&
      'not primitive, with multiplicity= ',multi,'.',ch10,&
      'This is allowed, as the input variable chkprim is 0.'
     MSG_COMMENT(msg)
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
!!      m_esymm,m_ingeo,m_symfind
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine symrelrot(nsym,rprimd,rprimd_new,symrel,tolsym)

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
 character(len=500) :: msg
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
!      Need to allow for four times tolsym, in case of centered Bravais lattices (but do it for all lattices ...)
       if(abs(val-nint(val))>four*tolsym)then
         write(msg,'(2a,a,i3,a,a,3es14.6,a,a,3es14.6,a,a,3es14.6)')&
         'One of the components of symrel is non-integer within 4*tolsym,',ch10,&
         '  for isym=',isym,ch10,&
         '  symrel=',matr2(:,1),ch10,&
         '         ',matr2(:,2),ch10,&
         '         ',matr2(:,3)
         MSG_ERROR_CLASS(msg, "TolSymError")
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
!!     fourth number is zero if the q-vector is not preserved, 1 otherwise
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
!!      m_bz_mesh,m_ddb,m_dtset,m_dvdb,m_dynmat,m_esymm,m_gkk,m_iogkk,m_lgroup
!!      m_memeval,m_nonlinear,m_phgamma,m_respfn_driver
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine littlegroup_q(nsym,qpt,symq,symrec,symafm,timrev,prtvol,use_sym)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(in),optional :: prtvol,use_sym
 integer,intent(out) :: timrev
!arrays
 integer,intent(in) :: symrec(3,3,nsym), symafm(nsym)
 integer,intent(out) :: symq(4,2,nsym)
 real(dp),intent(in) :: qpt(3)

!Local variables -------------------------
!scalars
 integer :: ii,isign,isym,itirev,my_prtvol
 real(dp),parameter :: tol=2.d-8
 real(dp) :: reduce
 character(len=500) :: msg
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

     ! Get the symmetric of the vector
     do ii=1,3
       qsym(ii)=qpt(1)*isign*symrec(ii,1,isym)&
         +qpt(2)*isign*symrec(ii,2,isym)&
         +qpt(3)*isign*symrec(ii,3,isym)
     end do

     ! Get the difference between the symmetric and the original vector
     symq(4,itirev,isym)=1
     do ii=1,3
       difq(ii)=qsym(ii)-qpt(ii)
       ! Project modulo 1 in the interval ]-1/2,1/2] such that difq = reduce + shift
       call wrap2_pmhalf(difq(ii),reduce,shift(ii))
       if(abs(reduce)>tol)symq(4,itirev,isym)=0
     end do

     ! SP: When prtgkk is asked (GKK matrix element will be output), one has to
     ! disable symmetries. There is otherwise a jauge problem with the unperturbed
     ! and the perturbed wavefunctions. This leads to a +- 5% increase in computational
     ! cost but provide the correct GKKs (i.e. the same as without the use of symmetries.)

     if (PRESENT(use_sym)) then
       if (use_sym == 0) then
         symq(4,itirev,isym)=0
         symq(4,itirev,1)=1
       end if
     end if

     ! If the operation succeded, change shift from real(dp) to integer, then exit loop
     if(symq(4,itirev,isym)/=0)then
       if (my_prtvol>0) then
         if(itirev==1)write(msg,'(a,i4,a)')' littlegroup_q : found symmetry',isym,' preserves q '
         if(itirev==2)write(msg,'(a,i4,a)')' littlegroup_q : found symmetry ',isym,' + TimeReversal preserves q '
         call wrtout(std_out,msg)
       end if
       ! Uses the mathematical function NINT = nearest integer
       do ii=1,3
         symq(ii,itirev,isym)=nint(shift(ii))
       end do
     end if

   end do !itirev
 end do !isym

 ! Test time-reversal symmetry
 timrev=1
 do ii=1,3
   ! Unfortunately, this version does not work yet ...
   ! call wrap2_pmhalf(2*qpt(ii),reduce,shift(ii))
   ! if(abs(reduce)>tol)timrev=0
   ! So, this is left ...
   if(abs(qpt(ii))>tol)timrev=0
 end do

 if(timrev==1.and.my_prtvol>0)then
   write(msg, '(3a)' )&
   ' littlegroup_q : able to use time-reversal symmetry. ',ch10,&
   '  (except for gamma, not yet able to use time-reversal symmetry)'
   call wrtout(std_out,msg)
 end if

end subroutine littlegroup_q
!!***

!!****f* m_symtk/matpointsym
!! NAME
!! matpointsym
!!
!! FUNCTION
!! For given order of point group, symmetrizes a 3x3 input matrix using the
!! point symmetry of the input atom
!!
!! INPUTS
!! iatom=index of atom to symmetrize around
!! natom=number of atoms in cell
!! nsym=order of group
!! rprimd(3,3)= real space primitive vectors
!! symrel(3,3,nsym)=symmetry operators in terms of action on primitive translations
!! tnons(3,nsym) = nonsymmorphic translations
!! xred(3,natom)=locations of atoms in reduced coordinates
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! mat3(3,3) = matrix to be symmetrized, in cartesian frame
!!
!! PARENTS
!!      m_nucprop,m_paw_nmr
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine matpointsym(iatom,mat3,natom,nsym,rprimd,symrel,tnons,xred)

 use m_linalg_interfaces

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iatom,natom,nsym
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: rprimd(3,3),tnons(3,nsym),xred(3,natom)
 real(dp),intent(inout) :: mat3(3,3)

!Local variables-------------------------------
!scalars
 integer :: cell_index,cell_indexp,ii,isym,nsym_point
 real(dp) :: xreddiff
!arrays
 integer :: symrel_it(3,3)
 real(dp) :: mat3_tri(3,3),mat3_tri_sym(3,3),rprimd_inv(3,3),tmp_mat(3,3)
 real(dp) :: xredp(3)

!**************************************************************************

!copy rprimd input and construct inverse
 rprimd_inv = rprimd
 call matrginv(rprimd_inv,3,3)

!transform input mat3 to triclinic frame with rprimd^{-1} * mat3 * rprimd
 call dgemm('N','N',3,3,3,one,rprimd_inv,3,mat3,3,zero,tmp_mat,3)
 call dgemm('N','N',3,3,3,one,tmp_mat,3,rprimd,3,zero,mat3_tri,3)

!loop over symmetry elements to obtain symmetrized input matrix
 mat3_tri_sym = zero
 nsym_point = 0
 do isym = 1, nsym

! skip any nonsymmorphic symmetry elements, want to consider point elements only
   if(dot_product(tnons(:,isym),tnons(:,isym))>tol8) cycle

! for current symmetry element, find transformed reduced coordinates of target atom
! via xredp = symrel * xred
   call dgemv('N',3,3,one,dble(symrel(:,:,isym)),3,xred(:,iatom),1,zero,xredp,1)


! shift xredp into the same unit cell as xred, for comparison
! label cells as 0..1:0 1..2:1 2..3:2 and -1..0:-1 -2..-1:-2 and so forth
   do ii = 1, 3

     cell_index = int(xred(ii,iatom))
     if(xred(ii,iatom) < zero) cell_index = cell_index - 1
     cell_indexp = int(xredp(ii))
     if(xredp(ii) < zero) cell_indexp = cell_indexp - 1

     do while (cell_indexp < cell_index)
       xredp(ii) = xredp(ii)+one
       cell_indexp = cell_indexp + 1
     end do
     do while (cell_indexp > cell_index)
       xredp(ii) = xredp(ii)-one
       cell_indexp = cell_indexp - 1
     end do

   end do

! now compare xredp to xred
   xreddiff = dot_product(xredp-xred(:,iatom),xredp-xred(:,iatom))

   if (xreddiff < tol8) then

!  accumulate symrel^{-1}*mat3_tri*symrel into mat3_tri_sym iff xredp = xred + L,
!  where is a lattice vector, so symrel leaves the target atom invariant

!  mati3inv gives the inverse transpose of symrel
     call mati3inv(symrel(:,:,isym),symrel_it)
     call dgemm('N','N',3,3,3,one,mat3_tri,3,dble(symrel(:,:,isym)),3,zero,tmp_mat,3)
     call dgemm('T','N',3,3,3,one,dble(symrel_it),3,tmp_mat,3,one,mat3_tri_sym,3)
     nsym_point = nsym_point + 1
   end if

 end do

!normalize by number of point symmetry operations
 mat3_tri_sym = mat3_tri_sym/dble(nsym_point)

!transform mat3_tri_sym to cartesian frame with rprimd * mat3_tri_sym * rprimd^{-1}

 call dgemm('N','N',3,3,3,one,mat3_tri_sym,3,rprimd_inv,3,zero,tmp_mat,3)
 call dgemm('N','N',3,3,3,one,rprimd,3,tmp_mat,3,zero,mat3,3)

end subroutine matpointsym
!!***

!!****f* m_symtk/holocell
!! NAME
!! holocell
!!
!! FUNCTION
!! Examine whether the trial conventional cell described by cell_base
!! is coherent with the required holohedral group.
!! Possibly enforce the holohedry and modify the basis vectors.
!! Note: for iholohedry=4, the tetragonal axis is not required to be along the C axis.
!!
!! INPUTS
!!  enforce= if 0, only check; if =1, enforce exactly the holohedry
!!  iholohedry=required holohegral group (uses its absolute value, since when the multiplicity of the cell is 
!!   more than one, the sign of iholohedry is changed).
!!  iholohedry=1   triclinic      1bar
!!  iholohedry=2   monoclinic     2/m
!!  iholohedry=3   orthorhombic   mmm
!!  iholohedry=4   tetragonal     4/mmm
!!  iholohedry=5   trigonal       3bar m
!!  iholohedry=6   hexagonal      6/mmm
!!  iholohedry=7   cubic          m3bar m
!!  tolsym=tolerance for the symmetry operations
!!
!! OUTPUT
!!  foundc=1 if the basis vectors supports the required holohedry ; =0 otherwise
!!
!! SIDE EFFECTS
!!  cell_base(3,3)=basis vectors of the conventional cell  (changed if enforce==1, otherwise unchanged)
!!
!! PARENTS
!!      m_symfind,m_symtk
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine holocell(cell_base,enforce,foundc,iholohedry,tolsym)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: enforce,iholohedry
 integer,intent(out) :: foundc
 real(dp),intent(in) :: tolsym
!arrays
 real(dp),intent(inout) :: cell_base(3,3)

!Local variables ------------------------------
!scalars
 integer :: allequal,ii,orth
 real(dp):: aa,scprod1
 character(len=500) :: msg
!arrays
 integer :: ang90(3),equal(3)
 real(dp) :: length(3),metric(3,3),norm(3),rbasis(3,3),rconv(3,3),rconv_new(3,3)
 real(dp) :: rnormalized(3,3),symmetrized_length(3)

!**************************************************************************

 if(abs(iholohedry)<1 .or. abs(iholohedry)>7)then
   write(msg, '(a,i0)' )&
&    'Abs(iholohedry) should be between 1 and 7, while iholohedry=',iholohedry 
   MSG_BUG(msg)
 end if

 do ii=1,3
   metric(:,ii)=cell_base(1,:)*cell_base(1,ii)+&
&   cell_base(2,:)*cell_base(2,ii)+&
&   cell_base(3,:)*cell_base(3,ii)
 end do

!Examine the angles and vector lengths
 ang90(:)=0
 if(metric(1,2)**2<tolsym**2*metric(1,1)*metric(2,2))ang90(3)=1
 if(metric(1,3)**2<tolsym**2*metric(1,1)*metric(3,3))ang90(2)=1
 if(metric(2,3)**2<tolsym**2*metric(2,2)*metric(3,3))ang90(1)=1
 orth=0
 if(ang90(1)==1 .and. ang90(2)==1 .and. ang90(3)==1) orth=1
 equal(:)=0
 if(abs(metric(1,1)-metric(2,2))<tolsym*half*(metric(1,1)+metric(2,2)))equal(3)=1
 if(abs(metric(1,1)-metric(3,3))<tolsym*half*(metric(1,1)+metric(3,3)))equal(2)=1
 if(abs(metric(2,2)-metric(3,3))<tolsym*half*(metric(2,2)+metric(3,3)))equal(1)=1
 allequal=0
 if(equal(1)==1 .and. equal(2)==1 .and. equal(3)==1) allequal=1

!DEBUG
!write(std_out, '(a,i4)' )' holocell : iholohedry=',iholohedry
!write(std_out, '(a,3i4)' )' holocell : ang90=',ang90
!write(std_out, '(a,3i4)' )' holocell : equal=',equal
!ENDDEBUG

 foundc=0
 if(abs(iholohedry)==1)                                      foundc=1
 if(abs(iholohedry)==2 .and. ang90(1)+ang90(3)==2 )          foundc=1
 if(abs(iholohedry)==3 .and. orth==1)                        foundc=1
 if(abs(iholohedry)==4 .and. orth==1 .and.          &
&  (equal(3)==1 .or. equal(2)==1 .or. equal(1)==1) ) foundc=1
 if(abs(iholohedry)==5 .and. allequal==1 .and. &
&  (abs(metric(1,2)-metric(2,3))<tolsym*metric(2,2)) .and. &
&  (abs(metric(1,2)-metric(1,3))<tolsym*metric(1,1))         )      foundc=1
 if(abs(iholohedry)==6 .and. equal(3)==1 .and. &
&   ang90(1)==1 .and. ang90(2)==1 .and. &
&   (2*metric(1,2)-metric(1,1))<tolsym*metric(1,1) )      foundc=1
 if(abs(iholohedry)==7 .and. orth==1 .and. allequal==1)      foundc=1

!DEBUG
!write(std_out, '(a,i4)' )' holocell : foundc=',foundc
!ENDDEBUG

!-------------------------------------------------------------------------------------
!Possibly enforce the holohedry (if it is to be enforced !)

 if(foundc==0.and.enforce==1.and.abs(iholohedry)/=1)then

!  Copy the cell_base vectors, and possibly fix the tetragonal axis to be the c-axis
!  XG20201016 WARNING : in principle, one should NOT use the 'equal' information, since precisely this enforcement
!  has the aim to reinstall the symmetries while they are broken !!
   if(abs(iholohedry)==4.and.equal(1)==1)then
     rconv(:,3)=cell_base(:,1) ; rconv(:,1)=cell_base(:,2) ; rconv(:,2)=cell_base(:,3)
   else if (abs(iholohedry)==4.and.equal(2)==1)then
     rconv(:,3)=cell_base(:,2) ; rconv(:,2)=cell_base(:,1) ; rconv(:,1)=cell_base(:,3)
   else
     rconv(:,:)=cell_base(:,:)
   end if

!  Compute the length of the three conventional vectors
   length(1)=sqrt(sum(rconv(:,1)**2))
   length(2)=sqrt(sum(rconv(:,2)**2))
   length(3)=sqrt(sum(rconv(:,3)**2))

!  Take care of the first conventional vector aligned with rbasis(:,3) (or aligned with the trigonal axis if rhombohedral)
!  and choice of the first normalized direction
   if(abs(iholohedry)==5)then
     rbasis(:,3)=third*(rconv(:,1)+rconv(:,2)+rconv(:,3))
   else
     rbasis(:,3)=rconv(:,3)
   end if
   norm(3)=sqrt(sum(rbasis(:,3)**2))
   rnormalized(:,3)=rbasis(:,3)/norm(3)

!  Projection of the first conventional vector perpendicular to rbasis(:,3)
!  and choice of the first normalized direction
   scprod1=sum(rnormalized(:,3)*rconv(:,1))
   rbasis(:,1)=rconv(:,1)-rnormalized(:,3)*scprod1
   norm(1)=sqrt(sum(rbasis(:,1)**2))
   rnormalized(:,1)=rbasis(:,1)/norm(1)

!  Generation of the second vector, perpendicular to the third and first
   rnormalized(1,2)=rnormalized(2,3)*rnormalized(3,1)-rnormalized(3,3)*rnormalized(2,1)
   rnormalized(2,2)=rnormalized(3,3)*rnormalized(1,1)-rnormalized(1,3)*rnormalized(3,1)
   rnormalized(3,2)=rnormalized(1,3)*rnormalized(2,1)-rnormalized(2,3)*rnormalized(1,1)

!  Compute the vectors of the conventional cell, on the basis of iholohedry
   if(abs(iholohedry)==2)then
     rconv_new(:,3)=rconv(:,3)
     rconv_new(:,1)=rconv(:,1)
     rconv_new(:,2)=rnormalized(:,2)*length(2) ! Now, the y axis is perpendicular to the two others, that have not been changed
   else if(abs(iholohedry)==3.or.abs(iholohedry)==4.or.abs(iholohedry)==7)then
     if(abs(iholohedry)==7)then
       symmetrized_length(1:3)=sum(length(:))*third
     else if(abs(iholohedry)==4)then
       symmetrized_length(3)=length(3)
       symmetrized_length(1:2)=half*(length(1)+length(2))
     else if(abs(iholohedry)==3)then
       symmetrized_length(:)=length(:)
     end if
     do ii=1,3
       rconv_new(:,ii)=rnormalized(:,ii)*symmetrized_length(ii)
     end do
   else if(abs(iholohedry)==5)then
!    In the normalized basis, they have coordinates (a,0,c), and (-a/2,+-sqrt(3)/2*a,c)
!    c is known, but a is computed from the knowledge of the average length of the initial vectors
     aa=sqrt(sum(length(:)**2)*third-norm(3)**2)
     rconv_new(:,1)=aa*rnormalized(:,1)+rbasis(:,3)
     rconv_new(:,2)=aa*half*(-rnormalized(:,1)+sqrt(three)*rnormalized(:,2))+rbasis(:,3)
     rconv_new(:,3)=aa*half*(-rnormalized(:,1)-sqrt(three)*rnormalized(:,2))+rbasis(:,3)
   else if(abs(iholohedry)==6)then

!    In the normalized basis, they have coordinates (a,0,0), (-a/2,+-sqrt(3)/2*a,0), and (0,0,c)
!    c is known, but a is computed from the knowledge of the average length of the initial vectors
     aa=half*(length(1)+length(2))
     rconv_new(:,1)=aa*rnormalized(:,1)
     rconv_new(:,2)=aa*half*(-rnormalized(:,1)+sqrt(three)*rnormalized(:,2))
     rconv_new(:,3)=rconv(:,3)
   end if

!  Copy back the cell_base vectors
   if(abs(iholohedry)==4.and.equal(1)==1)then
     cell_base(:,3)=rconv_new(:,2) ; cell_base(:,2)=rconv_new(:,1) ; cell_base(:,1)=rconv_new(:,3)
   else if (abs(iholohedry)==4.and.equal(2)==1)then
     cell_base(:,3)=rconv_new(:,1) ; cell_base(:,1)=rconv_new(:,2) ; cell_base(:,2)=rconv_new(:,3)
   else
     cell_base(:,:)=rconv_new(:,:)
   end if

 end if

end subroutine holocell
!!***

!!****f* m_symtk/symmetrize_rprimd
!! NAME
!! symmetrize_rprimd
!!
!! FUNCTION
!! Supposing the input rprimd does not preserve the length and angles
!! following the symmetries, will generates a new set rprimd,
!! on the basis of the expected characteristics of the conventional cell,
!! as specified in bravais(:)
!!
!! INPUTS
!! bravais(11): bravais(1)=iholohedry
!!              bravais(2)=center
!!              bravais(3:11)=coordinates of rprimd in the axes
!!              of the conventional bravais lattice (*2 if center/=0)
!! nsym=actual number of symmetries
!! symrel(3,3,1:nsym)=symmetry operations in real space in terms of primitive translations
!! tolsym=tolerance for the symmetry operations (only for checking purposes, the new set rprimd will
!!     be coherent with the symmetry operations at a much accurate level).
!!
!! SIDE EFFECTS
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!
!! PARENTS
!!      m_ingeo
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine symmetrize_rprimd(bravais,nsym,rprimd,symrel,tolsym)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 real(dp),intent(in) :: tolsym
!arrays
 integer,intent(in) :: bravais(11),symrel(3,3,nsym)
 real(dp),intent(inout) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: foundc,iexit,ii,jj
 real(dp):: rprimd_maxabs
!character(len=500) :: msg
!arrays
 real(dp):: aa(3,3),ait(3,3),cell_base(3,3),gprimd(3,3),rmet(3,3),rprimd_new(3,3)

! *************************************************************************

!DEBUG
!write(std_out,'(a)') ' symmetrize_rprimd : enter '
!ENDDEBUG

!Build the conventional cell basis vectors in cartesian coordinates
 aa(:,1)=bravais(3:5)
 aa(:,2)=bravais(6:8)
 aa(:,3)=bravais(9:11)
!Inverse transpose
 call matr3inv(aa,ait)
 do ii=1,3
   cell_base(:,ii)=ait(ii,1)*rprimd(:,1)+ait(ii,2)*rprimd(:,2)+ait(ii,3)*rprimd(:,3)
 end do

!DEBUG
!write(std_out,'(a)') ' before holocell, cell_base ='
!do ii=1,3
!  write(std_out,'(3es16.8)') cell_base(:,ii)
!enddo
!ENDDEBUG

!Enforce the proper holohedry on the conventional cell vectors.
 call holocell(cell_base,1,foundc,bravais(1),tolsym)

!DEBUG
!write(std_out,'(a)') ' after holocell, cell_base ='
!do ii=1,3
!  write(std_out,'(3es16.8)') cell_base(:,ii)
!enddo
!ENDDEBUG

!Reconstruct the dimensional primitive vectors
 do ii=1,3
   rprimd_new(:,ii)=aa(1,ii)*cell_base(:,1)+aa(2,ii)*cell_base(:,2)+aa(3,ii)*cell_base(:,3)
 end do

!Suppress meaningless values
 rprimd_maxabs=maxval(abs(rprimd_new))
 do ii=1,3
   do jj=1,3
     if(abs(rprimd(ii,jj))<tol12*rprimd_maxabs)rprimd(ii,jj)=zero
   enddo
 enddo

!! WRONG TEST
!Check whether the modification make sense
! do ii=1,3
!   do jj=1,3
!     reldiff=(rprimd_new(ii,jj)-rprimd(ii,jj))/sqrt(sum(rprimd(:,jj)**2))
!!    Allow for twice tolsym
!     if(abs(reldiff)>two*tolsym)then
!       write(msg,'(a,6(2a,3es14.6))')&
!!!         This is CRAZY : one detects symmetry problems above tolsym, and then requires the lattice vectors
!!!         not to be modify by more than 2 tolsym !!!
!&       'Failed rectification of lattice vectors to comply with Bravais lattice identification, modifs are too large',ch10,&
!&       '  rprimd    =',rprimd(:,1),ch10,&
!&       '             ',rprimd(:,2),ch10,&
!&       '             ',rprimd(:,3),ch10,&
!&       '  rprimd_new=',rprimd_new(:,1),ch10,&
!&       '             ',rprimd_new(:,2),ch10,&
!&       '             ',rprimd_new(:,3)
!       MSG_ERROR_CLASS(msg, "TolSymError")
!     end if
!   end do
! end do

 rprimd(:,:)=rprimd_new(:,:)

!Check whether the symmetry operations are consistent with the lattice vectors
 rmet = MATMUL(TRANSPOSE(rprimd), rprimd)
 call matr3inv(rprimd, gprimd)
 !call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 iexit=0

 call chkorthsy(gprimd,iexit,nsym,rmet,rprimd,symrel,tolsym)

!DEBUG
!write(std_out,'(a)') ' symmetrize_rprimd : exit '
!ENDDEBUG

end subroutine symmetrize_rprimd
!!***

!!****f* m_symtk/symmetrize_xred
!! NAME
!! symmetrize_xred
!!
!! FUNCTION
!! Symmetrize atomic coordinates using input symmetry matrices symrel
!! which are expressed in terms of the basis of real space primitive
!! translations (array elements are integers).
!! Input array indsym(4,isym,iatom) gives label of atom into which iatom
!! is rotated by INVERSE of symmetry element isym and also gives primitive
!! translation to get back to unit cell.
!! This version uses improvement in algorithm suggested by Andrew
!! Horsfield (see symatm.f).
!!
!! INPUTS
!! indsym(4,nsym,natom)=indirect indexing array giving label of atom
!!   into which iatom is rotated by symmetry element isym
!! natom=number of atoms
!! nsym=number of symmetries in group
!! symrel(3,3,nsym)=symmetry matrices in terms of real space
!!   primitive translations
!! tnons(3,nsym)=nonsymmorphic translations for symmetries
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! xred(3,natom)=
!!  (input) atomic coordinates in terms of real space translations
!!  (output) symmetrized atomic coordinates in terms
!!    of real space translations
!!
!! PARENTS
!!      m_ingeo,m_longwave,m_mover,m_nonlinear,m_respfn_driver,m_scfcv_core
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine symmetrize_xred(indsym,natom,nsym,symrel,tnons,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsym
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrel(3,3,nsym)
 real(dp),intent(in) :: tnons(3,nsym)
 real(dp),intent(inout) :: xred(3,natom)

!Local variables-------------------------------
!scalars
 integer  :: iatom,ib,isym
 integer  :: ii,jj
 real(dp) :: fc1,fc2,fc3
 real(dp) :: diff
 logical  :: dissimilar
!arrays
 real(dp) :: tsum(3),tt(3)
 real(dp),allocatable :: xredsym(:,:)
 real(dp) :: transl(3) ! translation vector

! *************************************************************************
!
!Check whether group contains more than identity;
!if not then simply return
 if (nsym>1) then

!  loop over atoms
   ABI_ALLOCATE(xredsym,(3,natom))
   do iatom=1,natom
     tsum(1)=0.0d0
     tsum(2)=0.0d0
     tsum(3)=0.0d0
!
!    loop over symmetries
     do isym=1,nsym
!      atom ib is atom into which iatom is rotated by inverse of
!      symmetry isym (inverse of symrel(mu,nu,isym))
       ib=indsym(4,isym,iatom)
!      Find the reduced coordinates after translation=t(indsym)+transl
       fc1=xred(1,ib)+dble(indsym(1,isym,iatom))
       fc2=xred(2,ib)+dble(indsym(2,isym,iatom))
       fc3=xred(3,ib)+dble(indsym(3,isym,iatom))
!      Compute [S * (x(indsym)+transl) ] + tnonsymmorphic
       tt(:)=dble(symrel(:,1,isym))*fc1+&
&       dble(symrel(:,2,isym))*fc2+&
&       dble(symrel(:,3,isym))*fc3+ tnons(:,isym)

!      Average over nominally equivalent atomic positions
       tsum(:)=tsum(:)+tt(:)
     end do
!
!    Set symmetrized result to sum over number of terms
     xredsym(:,iatom)=tsum(:)/dble(nsym)

!    End loop over iatom
   end do

   transl(:)=xredsym(:,1)-nint(xredsym(:,1))

!  Compute the smallest translation to an integer
   do jj=2,natom
     do ii=1,3
       diff=xredsym(ii,jj)-nint(xredsym(ii,jj))
       if (diff<transl(ii)) transl(ii)=diff
     end do
   end do

!  Test if the translation on each direction is small
!  Tolerance 1E-13
   do ii=1,3
     if (abs(transl(ii))>1e-13) transl(ii)=0.0
   end do

!  Execute translation
   do jj=1,natom
     do ii=1,3
       xredsym(ii,jj)=xredsym(ii,jj)-transl(ii)
     end do
   end do

!  Test if xredsym is too similar to xred
!  Tolerance 1E-15
   dissimilar=.FALSE.
   do jj=1,natom
     do ii=1,3
       if (abs(xredsym(ii,jj)-xred(ii,jj))>1E-15) dissimilar=.TRUE.
     end do
   end do

   if (dissimilar) xred(:,:)=xredsym(:,:)
   ABI_DEALLOCATE(xredsym)

!  End condition of nsym/=1
 end if

end subroutine symmetrize_xred
!!***

!!****f* m_symtk/symchk
!! NAME
!! symchk
!!
!! FUNCTION
!! Symmetry checker for atomic coordinates.
!! Checks for translated atomic coordinate tratom(3) to agree
!! with some coordinate xred(3,iatom) where atomic types agree too.
!! All coordinates are "reduced", i.e. given in terms of primitive
!! reciprocal translations.
!!
!! INPUTS
!! natom=number of atoms in unit cell
!! tratom(3)=reduced coordinates for a single atom which presumably
!!   result from the application of a symmetry operation to an atomic
!!   coordinate
!! trtypat=type of atom (integer) translated to tratom
!! typat(natom)=types of all atoms in unit cell (integer)
!! xred(3,natom)=reduced coordinates for all atoms in unit cell
!!
!! OUTPUT
!! difmin(3)=minimum difference between apparently equivalent atoms
!!   (give value separately for each coordinate)--note that value
!!   may be NEGATIVE so take abs later if needed
!! eatom=atom label of atom which is SAME as tratom to within a primitive
!!   cell translation ("equivalent atom")
!! transl(3)=primitive cell translation to make iatom same as tratom (integers)
!!
!! PARENTS
!!      m_polynomial_coeff,m_symtk
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine symchk(difmin,eatom,natom,tratom,transl,trtypat,typat,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,trtypat
 integer,intent(out) :: eatom
!arrays
 integer,intent(in) :: typat(natom)
 integer,intent(out) :: transl(3)
 real(dp),intent(in) :: tratom(3),xred(3,natom)
 real(dp),intent(out) :: difmin(3)

!Local variables-------------------------------
!scalars
 integer :: iatom,jatom,trans1,trans2,trans3
 real(dp) :: test,test1,test2,test3,testmn

! *************************************************************************

!DEBUG
! write(std_out,'(a,a,i4,3f18.12)') ch10,' symchk : enter, trtypat,tratom=',trtypat,tratom
!ENDDEBUG

!Start testmn out at large value
 testmn=1000000.d0

!Loop through atoms--
!when types agree, check for agreement after primitive translation
 jatom=1
 do iatom=1,natom
   if (trtypat/=typat(iatom)) cycle

!  Check all three components
   test1=tratom(1)-xred(1,iatom)
   test2=tratom(2)-xred(2,iatom)
   test3=tratom(3)-xred(3,iatom)
!  Find nearest integer part of difference
   trans1=nint(test1)
   trans2=nint(test2)
   trans3=nint(test3)
!  Check whether, after translation, they agree
   test1=test1-dble(trans1)
   test2=test2-dble(trans2)
   test3=test3-dble(trans3)
   test=abs(test1)+abs(test2)+abs(test3)
   if (test<tol10) then
!    Note that abs() is not taken here
     difmin(1)=test1
     difmin(2)=test2
     difmin(3)=test3
     jatom=iatom
     transl(1)=trans1
     transl(2)=trans2
     transl(3)=trans3
!    Break out of loop when agreement is within tolerance
     exit
   else
!    Keep track of smallest difference if greater than tol10
     if (test<testmn) then
       testmn=test
!      Note that abs() is not taken here
       difmin(1)=test1
       difmin(2)=test2
       difmin(3)=test3
       jatom=iatom
       transl(1)=trans1
       transl(2)=trans2
       transl(3)=trans3
     end if
   end if

!  End loop over iatom. Note a "cycle" and an "exit" inside the loop
 end do

 eatom=jatom

end subroutine symchk
!!***

!!****f* m_symtk/symatm
!! NAME
!! symatm
!!
!! FUNCTION
!! For each symmetry operation, find the number of the position to
!! which each atom is sent in the unit cell by the INVERSE of the
!! symmetry operation inv(symrel); i.e. this is the atom which, when acted
!! upon by the given symmetry element isym, gets transformed into atom iatom.
!!
!! This routine uses the fact that inv(symrel)=trans(symrec),
!! the inverse of the symmetry operation expressed in the basis of real
!! space primitive translations equals the transpose of the same symmetry
!! operation expressed in the basis of reciprocal space primitive transl:
!!
!!      $ xred(nu,indsym(4,isym,ia)) = symrec(mu,nu,isym)*(xred(mu,ia)-tnons(mu,isym)) - transl(mu)$
!!
!! where $transl$ is also a set of integers and
!! where translation transl places coordinates within unit cell (note sign).
!! Note that symrec is the set of arrays which are actually input here.
!! These arrays have integer elements.
!! tnons is the nonsymmorphic translation or else is zero.
!! If nsym=1 (i.e. only the identity symmetry is present) then
!! indsym merely takes each atom into itself.
!! The array of integer translations "transl" gets included within array "indsym" as seen below.
!! This routine has been improved using ideas of p. 649 of notes,
!! implementing suggestion of Andrew Horsfield: replace search for
!! equivalent atoms using direct primitive cell translations by
!! use of dot product relation which must produce an integer.
!! Relation:
!!
!!      $[inv(S(i))*(x(a)-tnons(i)) - x(inv(S)(i,a))] = integer$
!!
!! where $S(i) =$ symmetry matrix in real space, tnons=nonsymmorphic translation
!! (may be 0 0 0), and $x(inv(S)(i,a))$ is sought atom into which $x(a)$ gets
!! rotated by $inv(S)$.  Integer gives primitive translation coordinates to get
!! back to original unit cell.
!! Equivalent to $S*t(b)+tnons-x(a)=another$ $integer$ for $x(b)=x(inv(S))$.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! natom=number of atoms in cell.
!! nsym=number of space group symmetries.
!! symrec(3,3,nsym)=symmetries expressed in terms of their action on
!!                  reciprocal space primitive translations (integer).
!! tnons(3,nsym)=nonsymmorphic translations for each symmetry (would
!!               be 0 0 0 each for a symmorphic space group)
!! typat(natom)=integer identifying type of atom.
!! xred(3,natom)=reduced coordinates of atoms in terms of real space
!!               primitive translations
!! tolsym=tolerance for the symmetries
!! [print_indsym]: Print indsym table to std_out if the number of atoms is smaller that print_indsym
!!  Default: -1 i.e. no output is provided.
!!
!! OUTPUT
!! indsym(4,nsym,natom)=indirect indexing array described above: for each
!!                      isym,iatom, fourth element is label of atom into
!!                      which iatom is sent by INVERSE of symmetry operation
!!                      isym; first three elements are the primitive
!!                      translations which must be subtracted after the
!!                      transformation to get back to the original unit cell.
!!
!! PARENTS
!!      m_ab7_symmetry,m_berryphase_new,m_crystal,m_ddb,m_dtset,m_ingeo
!!      m_mover_effpot,m_orbmag,m_polynomial_coeff,m_spacepar,m_tdep_sym
!!      m_thmeig
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine symatm(indsym, natom, nsym, symrec, tnons, tolsym, typat, xred, print_indsym)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsym
 integer,optional,intent(in) :: print_indsym
 real(dp), intent(in) :: tolsym
!arrays
 integer,intent(in) :: symrec(3,3,nsym),typat(natom)
 integer,intent(out) :: indsym(4,nsym,natom)
 real(dp),intent(in) :: tnons(3,nsym),xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: eatom,errout,iatom,ii,isym,mu,print_indsym_
 real(dp) :: difmax,err
 character(len=500) :: msg
!arrays
 integer :: transl(3)
 real(dp) :: difmin(3),tratom(3)

! *************************************************************************

!DEBUG
!write(std_out,'(a,es12.4)')' symatm : enter, tolsym=',tolsym
!write(std_out,'(a,es12.4)')' symatm : xred='
!do ii=1,natom
!  write(std_out,'(i4,3es18.10)')ii,xred(1:3,ii)
!enddo
!write(std_out,'(a,es12.4)')' symatm : isym,symrec,tnons='
!do isym=1,nsym
!  write(std_out,'(i6,9i4,3es18.10)')isym,symrec(:,:,isym),tnons(1:3,isym)
!enddo
!ENDDEBUG

 err=zero
 errout=0

 do isym=1,nsym
   do iatom=1,natom

     do mu=1,3 ! Apply inverse transformation to original coordinates. Note transpose of symrec.
       tratom(mu) = dble(symrec(1,mu,isym))*(xred(1,iatom)-tnons(1,isym))&
&       +dble(symrec(2,mu,isym))*(xred(2,iatom)-tnons(2,isym))&
&       +dble(symrec(3,mu,isym))*(xred(3,iatom)-tnons(3,isym))
     end do
!
!    Find symmetrically equivalent atom
     call symchk(difmin,eatom,natom,tratom,transl,typat(iatom),typat,xred)
!
!    Put information into array indsym: translations and label
     indsym(1,isym,iatom)=transl(1)
     indsym(2,isym,iatom)=transl(2)
     indsym(3,isym,iatom)=transl(3)
     indsym(4,isym,iatom)=eatom
!
!    Keep track of maximum difference between transformed coordinates and
!    nearest "target" coordinate
     difmax=max(abs(difmin(1)),abs(difmin(2)),abs(difmin(3)))
     err=max(err,difmax)

     if(errout==3)then
       write(msg, '(a)' )&
       ' Suppress warning about finding symmetrically equivalent atoms, as mentioned already three times.'
       MSG_WARNING(msg)
       errout=errout+1
     endif

     if (difmax>tolsym .and. errout<3) then ! Print warnings if differences exceed tolerance
       write(msg, '(3a,i3,a,i6,a,i3,a,a,3f18.12,3a,es12.4)' )&
       ' Trouble finding symmetrically equivalent atoms',ch10,&
       ' Applying inv of symm number',isym,' to atom number',iatom,'  of typat',typat(iatom),ch10,&
       ' gives tratom=',tratom(1:3),'.',ch10,&
       ' This is further away from every atom in crystal than the allowed tolerance, tolsym=',tolsym
       MSG_WARNING(msg)

       write(msg, '(a,3i3,a,a,3i3,a,a,3i3)' ) &
       '  The inverse symmetry matrix is',symrec(1,1:3,isym),ch10,&
       '                                ',symrec(2,1:3,isym),ch10,&
       '                                ',symrec(3,1:3,isym)
       call wrtout(std_out,msg)
       write(msg, '(a,3f18.12)' )'  and the nonsymmorphic transl. tnons =',(tnons(mu,isym),mu=1,3)

       call wrtout(std_out,msg)
       write(msg, '(a,1p,3es12.4,a,a,i5)' ) &
        '  The nearest coordinate differs by',difmin(1:3),ch10,&
        '  for indsym(nearest atom)=',indsym(4,isym,iatom)
       call wrtout(std_out,msg)
!
!      Use errout to reduce volume of error diagnostic output
       if (errout==0) then
         write(msg,'(6a)') ch10,&
          '  This indicates that when symatm attempts to find atoms symmetrically',ch10, &
          '  related to a given atom, the nearest candidate is further away than some',ch10,&
          '  tolerance.  Should check atomic coordinates and symmetry group input data.'
         call wrtout(std_out,msg)
       end if
       errout=errout+1

     end if !difmax>tol
   end do !iatom
 end do !isym

 ! MG: Do not change this behaviour. symatm is called many times in the EPH code in which we have tons of q-points
 ! and it's really annoying to see this output repeated over and over again.
 ! If you need to print the indsym table at the beginning of the calculation, find the call to symatm
 ! and pass the optional argument print_indsym_ or use `abitk crystal_print FILE --prtvol 1`
 print_indsym_ = -1; if (present(print_indsym)) print_indsym_ = print_indsym
 if (natom <= print_indsym_) then
   do iatom=1,natom
     write(msg, '(a,i0,a)' )' symatm: atom number ',iatom,' is reached starting at atom'
     call wrtout(std_out,msg)
     do ii=1,(nsym-1)/24+1
       if(natom<100)then
         write(msg, '(1x,24i3)' ) (indsym(4,isym,iatom),isym=1+(ii-1)*24,min(nsym,ii*24))
       else
         write(msg, '(1x,24i6)' ) (indsym(4,isym,iatom),isym=1+(ii-1)*24,min(nsym,ii*24))
       end if
       call wrtout(std_out,msg)
     end do
   end do
 end if

 if (err>tolsym) then
   write(msg, '(1x,a,1p,e14.5,a,e12.4)' )'symatm: maximum (delta t)=',err,' is larger than tol=',tolsym
   MSG_WARNING(msg)
 end if

!Stop execution if error is really big
 if (err>0.01d0) then
   write(msg,'(5a)')&
    'Largest error (above) is so large (0.01) that either input atomic coordinates (xred)',ch10,&
    'are wrong or space group symmetry data is wrong.',ch10,&
    'Action: correct your input file.'
   MSG_ERROR(msg)
 end if

end subroutine symatm
!!***

!!****f* m_symtk/symcharac
!! NAME
!! symcharac
!!
!! FUNCTION
!! Get the type of axis for the symmetry.
!!
!! INPUTS
!! center=bravais(2)
!! determinant=the value of the determinant of sym
!! iholohedry=bravais(1)
!! isym=number of the symmetry operation that is currently analyzed
!! order=the order of the symmetry
!! symrel(3,3)= the symmetry matrix
!! tnons(3)=nonsymmorphic translations
!!
!! OUTPUT
!! label=a human readable text for the characteristic of the symmetry
!! type_axis=an identifier for the type of symmetry
!!
!! PARENTS
!!      m_ab7_symmetry,m_symfind
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine symcharac(center, determinant, iholohedry, isym, label, symrel, tnons, type_axis)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: determinant, center, iholohedry, isym
 integer, intent(out) :: type_axis
 character(len = 128), intent(out) :: label
 !arrays
 integer,intent(in) :: symrel(3,3)
 real(dp),intent(in) :: tnons(3)

 !Local variables-------------------------------
 !scalars
 logical,parameter :: verbose=.FALSE.
 integer :: tnons_order, identified, ii, order, iorder
 character(len=500) :: msg
 !arrays
 integer :: identity(3,3),matrix(3,3),trial(3,3)
 real(dp) :: reduced(3),trialt(3)

 !**************************************************************************

 identity(:,:)=0
 identity(1,1)=1 ; identity(2,2)=1 ; identity(3,3)=1
 trial(:,:)=identity(:,:)
 matrix(:,:)=symrel(:,:)

 order=0
 do iorder=1,6
   trial=matmul(matrix,trial)
   if(sum((trial-identity)**2)==0)then
     order=iorder
     exit
   end if
   if(sum((trial+identity)**2)==0)then
     order=iorder
     exit
   end if
 end do

 if(order==0)then
   type_axis = -2
   return
 end if

!Determination of the characteristics of proper symmetries (rotations)
 if(determinant==1)then

!  Determine the translation vector associated to the rotations
!  and its order : apply the symmetry operation
!  then analyse the resulting vector.
   identified=0
   trialt(:)=zero
   do ii=1,order
     trialt(:)=matmul(symrel(:,:),trialt(:))+tnons(:)
   end do
!  Gives the associated translation, with components in the
!  interval [-0.5,0.5] .
   reduced(:)=trialt(:)-nint(trialt(:)-tol6)

   if(sum(abs(reduced(:)))<tol6)identified=1
   if( (center==1 .or. center==-3) .and. &
&   sum(abs(reduced(:)-(/zero,half,half/)))<tol6 )identified=2
   if( (center==2 .or. center==-3) .and. &
&   sum(abs(reduced(:)-(/half,zero,half/)))<tol6 )identified=3
   if( (center==3 .or. center==-3) .and. &
&   sum(abs(reduced(:)-(/half,half,zero/)))<tol6 )identified=4
   if(center==-1.and. sum(abs(reduced(:)-(/half,half,half/)))<tol6 )identified=5

!  If the symmetry operation has not been identified, there is a problem ...
   if(identified==0) then
     type_axis = -1
     return
   end if

!  Compute the translation vector associated with one rotation
   trialt(:)=trialt(:)/order
   trialt(:)=trialt(:)-nint(trialt(:)-tol6)

!  Analyse the resulting vector.
   identified=0
   do ii=1,order
     reduced(:)=ii*trialt(:)-nint(ii*trialt(:)-tol6)
     if(sum(abs(reduced(:)))<tol6)identified=1
     if( (center==1 .or. center==-3) .and. &
&     sum(abs(reduced(:)-(/zero,half,half/)))<tol6 )identified=2
     if( (center==2 .or. center==-3) .and. &
&     sum(abs(reduced(:)-(/half,zero,half/)))<tol6 )identified=3
     if( (center==3 .or. center==-3) .and. &
&     sum(abs(reduced(:)-(/half,half,zero/)))<tol6 )identified=4
     if(center==-1.and. sum(abs(reduced(:)-(/half,half,half/)))<tol6 )identified=5

     if(identified/=0)then
       tnons_order=ii
       exit
     end if
   end do ! ii

!  Determinant (here=+1, as we are dealing with proper symmetry operations),
!  order, tnons_order and identified are enough to
!  determine the kind of symmetry operation

   select case(order)
   case(1)                       ! point symmetry 1
     if(identified==1) then
       type_axis=8                 ! 1
       write(label,'(a)') 'the identity'
     else
       type_axis=7                 ! t
       write(label,'(a)') 'a pure translation '
     end if

     if (verbose) then
       write(msg,'(a,i3,2a)')' symspgr : the symmetry operation no. ',isym,' is ',trim(label)
       call wrtout(std_out,msg)
     end if

   case(2,3,4,6)                 ! point symmetry 2,3,4,6 - rotations
     call symaxes(center,iholohedry,isym,symrel,label,order,tnons_order,trialt,type_axis)
   end select

 else if (determinant==-1)then

!  Now, take care of the improper symmetry operations.
!  Their treatment is relatively easy, except for the mirror planes
   select case(order)
   case(1)                       ! point symmetry 1
     type_axis=5                  ! -1
     write(label,'(a)') 'an inversion'
   case(2)                       ! point symmetry 2 - planes
     call symplanes(center,iholohedry,isym,symrel,tnons,label,type_axis)
   case(3)                       ! point symmetry 3
     type_axis=3                  ! -3
     write(label,'(a)') 'a -3 axis '
   case(4)                       ! point symmetry 1
     type_axis=2                  ! -4
     write(label,'(a)') 'a -4 axis '
   case(6)                       ! point symmetry 1
     type_axis=1                  ! -6
     write(label,'(a)') 'a -6 axis '
   end select

   if (order /= 2 .and. verbose) then
     write(msg,'(a,i3,2a)')' symspgr : the symmetry operation no. ',isym,' is ',trim(label)
     call wrtout(std_out,msg)
   end if

 end if ! determinant==1 or -1

end subroutine symcharac
!!***

!!****f* m_symtk/symaxes
!! NAME
!! symaxes
!!
!! FUNCTION
!! Determines the type of symmetry operation, for
!! the proper symmetries 2,2_1,3,3_1,3_2,4,4_1,4_2,4_3,6,6_1,...6_5
!!
!! INPUTS
!! center=type of bravais lattice centering
!!   	  center=0        no centering
!!        center=-1       body-centered
!!        center=-3       face-centered
!!        center=1        A-face centered
!!        center=2        B-face centered
!!        center=3        C-face centered
!! iholohedry=type of holohedry
!!            iholohedry=1   triclinic      1bar
!!            iholohedry=2   monoclinic     2/m
!!            iholohedry=3   orthorhombic   mmm
!!            iholohedry=4   tetragonal     4/mmm
!!            iholohedry=5   trigonal       3bar m  (rhombohedral Bravais latt)
!!            iholohedry=6   hexagonal      6/mmm
!!            iholohedry=7   cubic          m3bar m
!! isym=number of the symmetry operation that is currently analyzed
!! isymrelconv=symrel matrix for the particular operation, in conv. axes
!! ordersym=order of the symmetry operation
!! tnons_order=order of the screw translation
!! trialt(3)=screw translation associated with the symmetry operation
!!           in conventional axes (all components in the range ]-1/2,1/2] )
!!
!! OUTPUT
!! label=a user friendly label for the rotation
!! type_axis=type of the symmetry operation
!!
!! NOTES
!! It is assumed that the symmetry operations will be entered in the
!! symrel tnonsconv arrays, for the CONVENTIONAL cell.
!! For proper symmetries (rotations), the
!! associated translation is determined.
!!
!! There is a subtlety with translations associated with rotations :
!! all the rotations with axis
!! parallel to the one analysed do not all have the
!! same translation characteristics. This is clearly seen
!! in the extended Hermann-Mauguin symbols, see the international
!! table for crystallography, chapter 4.
!! In the treatment that we adopt, one will distinguish
!! the cases of primitive Bravais lattices, and centered
!! bravais lattices. In the latter case, in the present routine,
!! at the exception of the trigonal axis for the
!! cubic system, we explicitely generate the correct ratio of different
!! translations, so that their type can be explicitely assigned,
!! without confusion. By contrast, for primitive lattices,
!! the "tnons" that has been transmitted to the present routine
!! might be one of the few possible translations vectors,
!! nearly at random. We deal with this case by the explicit
!! examination of the system classes, and the identification
!! of such a possibility. In particular:
!! (1) for the trigonal axis in the rhombohedral Bravais lattice,
!! or in the cubic system, there is an equal number of 3, 3_1,
!! and 3_2 axes parallel to each other, in a cell that
!! is primitive (as well as conventional). In this particular case,
!! in the present
!! routine, all 3, 3_1 and 3_2 axes are assigned to be 3 axes,
!! independently of the centering.
!! (2) for the 4- or 6- axes, no confusion is possible :
!! in the primitive cell, there is only one possible translation,
!! while in the centered cells, the correct ratio of translation
!! vectors will be generated
!! (3) for the binary axes, there is no problem when the cell
!! is centered, but there are problems
!! (3a) for the tP Bravais lattice, for an axis in a tertiary direction,
!! (see the description of the lattice symmetry directions
!!  table 2.4.1 of the international tables for crystallography),
!!  where the family of axes is made equally of 2 and 2_1 axis.
!!  In this case, we attribute the binary axis to the specific class
!!  of "tertiary 2-axis". We keep track of the 2 or 2_1
!!  characteristics of all other binary axes
!! (3b) for the tI Bravais lattice, in all the directions,
!!  there is an equal number of 2 and 2_1 axes. We distinguish
!!  the primary and secondary family from the tertiary family.
!! (3c) for the hP Bravais lattice, each binary axis can present
!!  no translation or be a screw axis (in the same direction).
!!  For primary axes, one need the "2" and "2_1" classification,
!!  while for secondary and tertiary axes, the associated
!!  translation vector will have not importance.
!!  However, one will need to distinguish secondary from
!!  tertiary, and these from primary axes.
!!  So, this is the most complicated case, for binary axes,
!!  with the following sets of binary axes : "2", "2_1",
!!  "secondary 2" and "tertiary 2".
!! (3d) for the hR Bravais lattice, each binary axis can present
!!  no translation or be a screw axis (in the same direction).
!!  There is no distinction between tertiary axes and other, so that
!!  we simply assign a binary axis to "2-axis"
!! (3e) for the cP lattice, the binary axes along tertiary directions
!!  can also have different translation vectors, while for the primary
!!  direction, there is no such ambiguity. So, we will attribute
!!  tertiary 2 axis to the "tertiary 2-axis" set (there are always 6),
!!  and attribute 2 and 2_1 primary axes to the corresponding sets.
!!
!! PARENTS
!!      m_symtk
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine symaxes(center,iholohedry,isym,isymrelconv,label,ordersym,tnons_order,trialt,type_axis)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: center,iholohedry,isym,ordersym,tnons_order
 integer,intent(out) :: type_axis
 character(len=128),intent(out) :: label
!arrays
 integer,intent(in) :: isymrelconv(3,3)
 real(dp),intent(in) :: trialt(3)

!Local variables-------------------------------
!scalars
 logical,parameter :: verbose=.FALSE.
 character(len=500) :: msg
 integer :: direction,directiontype
 real(dp),parameter :: nzero=1.0d-6

!**************************************************************************

!write(std_out,*)' symaxes : enter, isym=',isym
!write(std_out,*)' symaxes : iholohedry, ',iholohedry
!write(std_out,*)' symaxes : center, ',center

 select case(ordersym)

 case(2)                       ! point symmetry 2
!    Must characterize directiontype for cP, tP, tI, and hP Bravais lattices
   directiontype=1
   if( iholohedry==4 .or. iholohedry==7) then ! tP or cP Bravais lattices
     if(abs(isymrelconv(1,1))+ &
&     abs(isymrelconv(2,2))+ &
&     abs(isymrelconv(3,3))  ==1) directiontype=3
   else if(iholohedry==6)then   ! hP Bravais lattice
     if(sum(isymrelconv(:,:))/=-1 )directiontype=2
     if(sum(isymrelconv(:,:))==0 .or. sum(isymrelconv(:,:))==-3 )&
&     directiontype=3
!      directiontype=1 corresponds to a primary axis
!      directiontype=2 corresponds to a tertiary axis
!      directiontype=3 corresponds to a secondary axis
   end if

!    DEBUG
!    write(std_out,*)' directiontype=',directiontype
!    write(std_out,'(a,3i6)' )' isymrelconv(1:3)=',isymrelconv(:,1)
!    write(std_out,'(a,3i6)' )' isymrelconv(4:6)=',isymrelconv(:,2)
!    write(std_out,'(a,3i6)' )' isymrelconv(7:9)=',isymrelconv(:,3)
!    write(std_out,'(a,i)' )' tnons_order=',tnons_order
!    ENDDEBUG

!    Now, classify the 2 axes
   if(directiontype==2)then
     type_axis=4                 ! secondary 2  (only in the hP Bravais latt case)
     write(label,'(a)') 'a secondary 2-axis '

   else if(directiontype==3 .and. iholohedry==4)then
     type_axis=21                ! tertiary 2
     write(label,'(a)') 'a tertiary 2-axis '
   else if(directiontype==3 .and. &
&     center==0 .and. (iholohedry==6.or.iholohedry==7) )then
     type_axis=21                ! tertiary 2
     write(label,'(a)') 'a tertiary 2-axis '
   else if(tnons_order==1 .or. (iholohedry==4 .and. center==-1) .or. &
&     iholohedry==5)then
     type_axis=9                 ! 2
     write(label,'(a)') 'a 2-axis '
   else
     type_axis=20                ! 2_1
     write(label,'(a)') 'a 2_1-axis '
   end if

 case(3)                       ! point symmetry 3
   if(tnons_order==1)then
     type_axis=10                ! 3
     write(label,'(a)') 'a 3-axis '
   else if(iholohedry==5 .or. iholohedry==7)then
!      This is a special situation : in the same family of parallel 3-axis,
!      one will have an equal number of 3, 3_1 and 3_2 axes, so that
!      it is non-sense to try to classify one of them.
     type_axis=10                ! 3, 3_1 or 3_2, undistinguishable
     write(label,'(a)') 'a 3, 3_1 or 3_2 axis '
   else
!      DEBUG
!      write(std_out,*)'isymrelconv=',isymrelconv(:,:)
!      write(std_out,*)'trialt=',trialt(:)
!      ENDDEBUG
!      Must recognize 3_1 or 3_2
     if(isymrelconv(1,1)==0)then  ! 3+
       if(abs(trialt(3)-third)<nzero)type_axis=22   ! 3_1
       if(abs(trialt(3)+third)<nzero)type_axis=23   ! 3_2
     else if(isymrelconv(1,1)==-1)then  ! 3-
       if(abs(trialt(3)-third)<nzero)type_axis=23   ! 3_2
       if(abs(trialt(3)+third)<nzero)type_axis=22   ! 3_1
     end if
     write(label,'(a)') 'a 3_1 or 3_2-axis '
   end if

 case(4)                       ! point symmetry 4
   if(tnons_order==1)then
     type_axis=12                ! 4
     write(label,'(a)') 'a 4-axis '
   else if(tnons_order==2)then
     type_axis=25                ! 4_2
     write(label,'(a)') 'a 4_2-axis '
   else if(center/=0)then
     type_axis=24                ! 4_1 or 4_3
     write(label,'(a)') 'a 4_1 or 4_3-axis '
   else
!      DEBUG
!      write(std_out,*)'isymrelconv=',isymrelconv(:,:)
!      write(std_out,*)'trialt=',trialt(:)
!      ENDDEBUG
!      Must recognize 4_1 or 4_3, along the three primary directions
     do direction=1,3
       if(isymrelconv(direction,direction)==1)then  !
         if( (direction==1 .and. isymrelconv(2,3)==-1) .or. &
&         (direction==2 .and. isymrelconv(3,1)==-1) .or. &
&         (direction==3 .and. isymrelconv(1,2)==-1)       )then ! 4+
           if(abs(trialt(direction)-quarter)<nzero)type_axis=24    ! 4_1
           if(abs(trialt(direction)+quarter)<nzero)type_axis=26    ! 4_3
         else if( (direction==1 .and. isymrelconv(2,3)==1) .or. &
&           (direction==2 .and. isymrelconv(3,1)==1) .or. &
&           (direction==3 .and. isymrelconv(1,2)==1)       )then ! 4-
           if(abs(trialt(direction)-quarter)<nzero)type_axis=26    ! 4_3
           if(abs(trialt(direction)+quarter)<nzero)type_axis=24    ! 4_1
         end if
       end if
     end do
     write(label,'(a)') 'a 4_1 or 4_3-axis '
   end if

 case(6)                       ! point symmetry 6
   if(tnons_order==1)then
     type_axis=14                ! 6
     write(label,'(a)') 'a 6-axis '
   else if(tnons_order==2)then
     type_axis=29                ! 6_3
     write(label,'(a)') 'a 6_3-axis '
   else if(tnons_order==3)then
     !write(std_out,*)'isymrelconv=',isymrelconv(:,:)
     !write(std_out,*)'trialt=',trialt(:)
     !Must recognize 6_2 or 6_4
     if(isymrelconv(1,1)==1)then  ! 6+
       if(abs(trialt(3)-third)<nzero)type_axis=28   ! 6_2
       if(abs(trialt(3)+third)<nzero)type_axis=30   ! 6_4
     else if(isymrelconv(1,1)==0)then  ! 6-
       if(abs(trialt(3)-third)<nzero)type_axis=30   ! 6_4
       if(abs(trialt(3)+third)<nzero)type_axis=28   ! 6_2
     end if
     write(label,'(a)') 'a 6_2 or 6_4-axis '
   else
     !write(std_out,*)'isymrelconv=',isymrelconv(:,:)
     !write(std_out,*)'trialt=',trialt(:)
     !Must recognize 6_1 or 6_5
     if(isymrelconv(1,1)==1)then  ! 6+
       if(abs(trialt(3)-sixth)<nzero)type_axis=27   ! 6_1
       if(abs(trialt(3)+sixth)<nzero)type_axis=31   ! 6_5
     else if(isymrelconv(1,1)==0)then  ! 6-
       if(abs(trialt(3)-sixth)<nzero)type_axis=31   ! 6_5
       if(abs(trialt(3)+sixth)<nzero)type_axis=27   ! 6_1
     end if
     write(label,'(a)') 'a 6_1 or 6_5-axis '
   end if

 end select

 if (verbose) then
   write(msg,'(a,i3,a,a)')' symaxes : the symmetry operation no. ',isym,' is ', trim(label)
   call wrtout(std_out,msg)
 end if

end subroutine symaxes
!!***

!!****f* m_symtk/symplanes
!! NAME
!! symplanes
!!
!! FUNCTION
!! Determines the type of symmetry mirror planes: m,a,b,c,d,n,g.
!! This is used (see symlist.f) to identify the space group.
!!
!! INPUTS
!! center=type of bravais lattice centering
!!   center=0        no centering
!!   center=-1       body-centered
!!   center=-3       face-centered
!!   center=1        A-face centered
!!   center=2        B-face centered
!!   center=3        C-face centered
!! iholohedry=type of holohedry
!!   iholohedry=1   triclinic      1bar
!!   iholohedry=2   monoclinic     2/m
!!   iholohedry=3   orthorhombic   mmm
!!   iholohedry=4   tetragonal     4/mmm
!!   iholohedry=5   trigonal       3bar m
!!   iholohedry=6   hexagonal      6/mmm
!!   iholohedry=7   cubic          m3bar m
!! isym=number of the symmetry operation that is currently analyzed
!! isymrelconv=symrel matrix for the particular operation, in conv. coord.
!! itnonsconv=tnons vector for the particular operation, in conv. coord
!!
!! OUTPUT
!! label=user friendly label of the plane
!! type_axis=type of the symmetry operation
!!
!! NOTES
!! One follows the
!! conventions explained in table 1.3 of the international tables for
!! crystallography. In the case of the rhombohedral system,
!! one takes into account the first footnote of this table 1.3 .
!! In general, we will assign the different symmetries to
!! the following numbers :  m -> 15 , (a, b or c) -> 16,
!!  d -> 17, n -> 18 , g -> 19
!! However, there is the same problem as for binary axes,
!! namely, for parallel mirror planes, one can find different
!! translation vectors, and these might be found at random,
!! depending on the input tnons.
!! (1) In the tP case, one will distinguish tertiary
!!  mirror plane, for which it is important to know whether they are
!!  m or c (for tertiary planes in tP, g is equivalent to m and n is equivalent to c).
!!  On the other hand, it is important to distinguish among
!!  primary and secondary mirror planes, those that are m,(a or b),c, or n.
!!  To summarize, the number of the symmetry will be :
!!  m (primary, secondary or tertiary) -> 15 ,
!!  secondary (a or b) -> 16, secondary c -> 17,
!!  primary or secondary n -> 18 , tertiary c -> 19
!! (2) In the tI case, one will distinguish tertiary
!!  mirror plane, for which it is important to know whether they are
!!  m or d (for tertiary planes in tI, c is equivalent to m.
!!  On the other hand, it is important to distinguish among
!!  primary and secondary mirror planes, those that are m (equivalent to n),
!!  or a,b or c.
!!  To summarize, the number of the symmetry will be :
!!  m (primary, secondary, tertiary) -> 15 ,
!!  a,b or c (primary or secondary) -> 16, tertiary d -> 17
!! (3) For hP and hR, a m plane is always coupled to a a or b plane,
!!  while a c plane is always coupled to an n plane. On the other
!!  hand, it is important to distinguish between primary or secondary
!!  mirror planes, and tertiary mirror planes. So we will keep the
!!  following sets : m non-tertiary (that includes a or b non-tertiary) -> 15,
!!  c non-tertiary (that includes n non-tertiary) -> 16,
!!  m tertiary (that includes a or b non-tertiary) -> 17,
!!  c tertiary (that includes n non-tertiary) -> 18.
!!  For hR, all mirror planes are secondary.
!! (4) For the cP lattice, in the same spirit, one can see that
!!  the tertiary m and g mirror planes are to be classified as "m" -> 15,
!!  while n, a and c are to be classified as "n" -> 18. There is no need
!!  to distinguish between primary, secondary or tertiary axes.
!!
!! PARENTS
!!      m_symtk
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine symplanes(center,iholohedry,isym,isymrelconv,itnonsconv,label,type_axis)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: center,iholohedry,isym
 integer,intent(out) :: type_axis
 character(len = 128), intent(out) :: label
!arrays
 integer,intent(in) :: isymrelconv(3,3)
 real(dp),intent(in) :: itnonsconv(3)

!Local variables-------------------------------
!scalars
 logical,parameter :: verbose=.FALSE.
 character(len=500) :: msg
 integer :: directiontype,sum_elements
 real(dp),parameter :: nzero=1.0d-6
!arrays
 integer :: identity(3,3),mirrormxy(3,3),mirrormyz(3,3),mirrormzx(3,3)
 integer :: mirrorx(3,3),mirrorxy(3,3),mirrory(3,3),mirroryz(3,3),mirrorz(3,3)
 integer :: mirrorzx(3,3)
 real(dp) :: trialt(3)
! real(dp) :: itnonsconv2(3),trialt2(3)

!**************************************************************************

!write(std_out,*)' symplanes : enter'
!write(std_out,*)' center,iholohedry,isym,isymrelconv,itnonsconv=',center,iholohedry,isym,isymrelconv,itnonsconv

 identity(:,:)=0
 identity(1,1)=1 ; identity(2,2)=1 ; identity(3,3)=1

!Will be a mirror plane, but one must characterize
!(1) the type of plane (primary, secondary or tertiary)
!(2) the gliding vector. One now defines a few matrices.
 mirrorx(:,:)=identity(:,:) ; mirrorx(1,1)=-1
 mirrory(:,:)=identity(:,:) ; mirrory(2,2)=-1
 mirrorz(:,:)=identity(:,:) ; mirrorz(3,3)=-1
 mirrorxy(:,:)=0 ; mirrorxy(1,2)=1 ; mirrorxy(2,1)=1 ; mirrorxy(3,3)=1
 mirrorzx(:,:)=0 ; mirrorzx(1,3)=1 ; mirrorzx(3,1)=1 ; mirrorzx(2,2)=1
 mirroryz(:,:)=0 ; mirroryz(2,3)=1 ; mirroryz(3,2)=1 ; mirroryz(1,1)=1
 mirrormxy(:,:)=0 ; mirrormxy(1,2)=-1 ; mirrormxy(2,1)=-1 ; mirrormxy(3,3)=1
 mirrormzx(:,:)=0 ; mirrormzx(1,3)=-1 ; mirrormzx(3,1)=-1 ; mirrormzx(2,2)=1
 mirrormyz(:,:)=0 ; mirrormyz(2,3)=-1 ; mirrormyz(3,2)=-1 ; mirrormyz(1,1)=1

!Determine the type of plane. At the end,
!directiontype=1 will correspond to a primary axis (or equivalent
!axes for orthorhombic)
!directiontype=2 will correspond to a secondary axis
!directiontype=3 will correspond to a tertiary axis
!See table 2.4.1, 11.2 and 11.3 of the international tables for crystallography
 directiontype=0
!The sum of elements of the matrices allow to characterize them
 sum_elements=sum(isymrelconv(:,:))

 if(sum_elements==1)then
!  The mirror plane perpendicular to the c axis is always primary
   if( sum(abs(isymrelconv(:,:)-mirrorz(:,:)))==0 )then
     directiontype=1
!    All the other planes with a symrel matrix whose sum of elements is 1
!    are a or b planes. They are primary or
!    secondary planes, depending the holohedry.
   else if(sum(isymrelconv(:,:))==1)then
     if( iholohedry==2 .or. iholohedry==3 .or. iholohedry==7 )then
       directiontype=1
     else if(iholohedry==4 .or. iholohedry==6)then
       directiontype=2
     end if
   end if
 end if

!All the planes with a symrel matrix whose sum of elements
!is 2 are secondary planes (table 11.3).
 if( sum_elements==2 ) directiontype=2

!The planes with a symrel matrix whose sum of elements
!is 3 or 0 are tertiary planes
 if( sum_elements==3 .or. sum_elements==0 )directiontype=3

!One is left with sum_elements=-1, tertiary for tetragonal
!or cubic, secondary for hexagonal
 if( sum_elements==-1)then
   if(iholohedry==4 .or. iholohedry==7)directiontype=3
   if(iholohedry==6)directiontype=2
 end if


!Now, determine the gliding vector
!First, apply the symmetry operation
!to itnonsconv, in order to get the translation vector
!under the application of twice the symmetry operation
 trialt(:)=matmul(isymrelconv(:,:),itnonsconv(:)) +itnonsconv(:)
!Get the translation associated with one application,
!and force its components to be in the interval ]-0.5,0.5] .
 trialt(:)=trialt(:)*half
 trialt(:)=trialt(:)-nint(trialt(:)-nzero)

!If there is a glide vector for the initial choice of itnonsconv,
!it might be that it disappears if itnonsconv is translated by a
!lattice vector of the conventional cell
!if(trialt(1)**2+trialt(2)**2+trialt(3)**2>tol5)then
!do ii=1,3
!itnonsconv2(:)=itnonsconv(:)
!itnonsconv2(ii)=itnonsconv(ii)+one
!trialt2(:)=matmul(isymrelconv(:,:),itnonsconv2(:)) +itnonsconv2(:)
!trialt2(:)=trialt2(:)*half
!trialt2(:)=trialt2(:)-nint(trialt2(:)-nzero)
!if(trialt2(1)**2+trialt2(2)**2+trialt2(3)**2<tol5)then
!trialt(:)=trialt2(:)
!endif
!enddo
!endif

 write(msg,'(a)') ' symplanes...'

!Must use the convention of table 1.3 of the international
!tables for crystallography, see also pp 788 and 789.
!Often, one needs to specialize the selection according
!to the Bravais lattice or the system.

 if(sum(abs(trialt(:)))<nzero .and. iholohedry/=6)then
   type_axis=15  ! m
   write(label,'(a)') 'a mirror plane'
 else if(iholohedry==4 .and. center==0)then    ! primitive tetragonal

   if(directiontype==1)then
     type_axis=18  ! primary n
     write(label,'(a)') 'a primary n plane'
   else if(directiontype==2)then
     if(sum(abs(trialt(:)-(/half,zero,zero/)))<nzero .or. sum(abs(trialt(:)-(/zero,half,zero/)))<nzero)then
       type_axis=16  ! secondary a or b
       write(label,'(a)') 'a secondary a or b plane'
     else if(sum(abs(trialt(:)-(/zero,zero,half/)))<nzero)then
       type_axis=17    ! secondary c
       write(label,'(a)') 'a secondary c plane'
     else
       type_axis=18    ! secondary n
       write(label,'(a)') 'a secondary n plane'
     end if ! directiontype==2
   else if(directiontype==3)then
     if( abs(trialt(3))<nzero )then
       type_axis=15    ! tertiary m
       write(label,'(a)') 'a tertiary m plane'
     else if( abs(trialt(3)-half)<nzero )then
       type_axis=19    ! tertiary c
       write(label,'(a)') 'a tertiary c plane'
     end if
   end if

 else if(iholohedry==4 .and. center==-1)then    ! inner tetragonal

   if(directiontype==1 .or. directiontype==2)then
     if(sum(abs(trialt(:)-(/half,zero,zero/)))<nzero .or. &
        sum(abs(trialt(:)-(/zero,half,zero/)))<nzero .or. &
        sum(abs(trialt(:)-(/zero,zero,half/)))<nzero      )then
       type_axis=16    ! a, b, or c
       write(label,'(a)') 'an a, b or c plane'
     else if(sum(abs(trialt(:)-(/half,half,zero/)))<nzero .or. &
             sum(abs(trialt(:)-(/zero,half,half/)))<nzero .or. &
             sum(abs(trialt(:)-(/half,zero,half/)))<nzero       )then
       type_axis=15    ! n plane, equivalent to m
       write(label,'(a)') 'a m plane'
     end if ! directiontype==1 or 2
   else if(directiontype==3)then
     if( abs(trialt(3))<nzero .or. abs(trialt(3)-half)<nzero )then
       type_axis=15    ! tertiary c, equivalent to m
       write(label,'(a)') 'a tertiary m plane'
     else
       type_axis=17    ! tertiary d
       write(label,'(a)') 'a tertiary d plane'
     end if
   end if

 else if(iholohedry==5)then    ! hR

   if( abs(sum(abs(trialt(:)))-one) < nzero) then
     type_axis=15    ! secondary m
     write(label,'(a)') 'a secondary m plane'
   else if( abs(sum(abs(trialt(:)))-half) < nzero .or. abs(sum(abs(trialt(:)))-three*half) < nzero )then
     type_axis=16    ! secondary c
     write(label,'(a)') 'a secondary c plane'
   end if

 else if(iholohedry==6)then    ! hP

   if(directiontype==1)then
     if( abs(trialt(3)) <nzero )then
       type_axis=15    ! primary m
       write(label,'(a)') 'a primary m plane'
     end if
   else if(directiontype==2)then
     if( abs(trialt(3)) <nzero )then
       type_axis=15    ! secondary m
       write(label,'(a)') 'a secondary m plane'
     else if( abs(trialt(3)-half) < nzero ) then
       type_axis=16    ! secondary c
       write(label,'(a)') 'a secondary c plane'
     end if
   else if(directiontype==3)then
     if( abs(trialt(3)) <nzero )then
       type_axis=17    ! tertiary m
       write(label,'(a)') 'a tertiary m plane'
     else if( abs(trialt(3)-half) < nzero ) then
       type_axis=18    ! tertiary c
       write(label,'(a)') 'a tertiary c plane'
     end if
   end if ! directiontype

!  else if(iholohedry==7 .and. center==0)then    ! cP
 else if(iholohedry==7)then    ! cP

   if(directiontype==1)then
     if((sum(abs(isymrelconv(:,:)-mirrorx(:,:)))==0 .and.  &
         sum(abs(two*abs(trialt(:))-(/zero,half,half/)))<nzero   ).or. &
         (sum(abs(isymrelconv(:,:)-mirrory(:,:)))==0 .and.  &
         sum(abs(two*abs(trialt(:))-(/half,zero,half/)))<nzero   ).or. &
         (sum(abs(isymrelconv(:,:)-mirrorz(:,:)))==0 .and.  &
         sum(abs(two*abs(trialt(:))-(/half,half,zero/)))<nzero   )    ) then
       type_axis=17     ! d
       write(label,'(a)') 'a d plane'
     else
       type_axis=18    ! primary n
       write(label,'(a)') 'a primary n plane'
     end if
   else if(directiontype==3)then
     if(sum(abs(two*abs(trialt(:))-(/half,half,half/)))<nzero       )then
       type_axis=17     ! d
       write(label,'(a)') 'a d plane'
     else if( abs(sum(abs(trialt(:)))-half) < nzero .or. abs(sum(abs(trialt(:)))-three*half) < nzero ) then
       type_axis=18    ! tertiary n
       write(label,'(a)') 'a tertiary n plane'
     else if( abs(sum(abs(trialt(:)))-one) < nzero )then
       type_axis=15    ! tertiary m
       write(label,'(a)') 'a tertiary m plane'
     end if
   end if

!  Now, treat all other cases (including other centered Bravais lattices)
 else if(sum(abs(trialt(:)-(/half,zero,zero/)))<nzero .or. &
         sum(abs(trialt(:)-(/zero,half,zero/)))<nzero .or. &
         sum(abs(trialt(:)-(/zero,zero,half/)))<nzero       )then
   type_axis=16     ! a, b or c
   write(label,'(a)') 'an a,b, or c plane'
 else if( (directiontype==1 .or. directiontype==2) .and. &
          (sum(abs(trialt(:)-(/half,half,zero/)))<nzero .or. &
           sum(abs(trialt(:)-(/zero,half,half/)))<nzero .or. &
           sum(abs(trialt(:)-(/half,zero,half/)))<nzero     ) )then
   type_axis=18     ! n
   write(label,'(a)') 'an n plane'
 else if( directiontype==3 .and. sum(abs(trialt(:)-(/half,half,half/)))<nzero )then
   type_axis=18     ! n
   write(label,'(a)') 'an n plane'
 else if((sum(abs(isymrelconv(:,:)-mirrorx(:,:)))==0 .and.  &
          sum(abs(two*abs(trialt(:))-(/zero,half,half/)))<nzero   ).or. &
          (sum(abs(isymrelconv(:,:)-mirrory(:,:)))==0 .and.  &
          sum(abs(two*abs(trialt(:))-(/half,zero,half/)))<nzero   ).or. &
          (sum(abs(isymrelconv(:,:)-mirrorz(:,:)))==0 .and.  &
          sum(abs(two*abs(trialt(:))-(/half,half,zero/)))<nzero   )    ) then
   type_axis=17     ! d
   write(label,'(a)') 'a d plane'
 else if( directiontype==3 .and. sum(abs(two*abs(trialt(:))-(/half,half,half/)))<nzero)then
   type_axis=17     ! d
   write(label,'(a)') 'a d plane'
 else
   type_axis=19     ! g (all other planes with
!  unconventional glide vector)
   write(label,'(a)') 'a g plane'
 end if

 if (verbose) then
   write(msg,'(a,i3,a,a)')' symplanes : the symmetry operation no. ',isym,' is ', trim(label)
   call wrtout(std_out,msg)
 end if

end subroutine symplanes
!!***

!!****f* m_symtk/smallprim
!!
!! NAME
!! smallprim
!!
!! FUNCTION
!! Find the smallest possible primitive vectors for an input lattice
!! This algorithm is not as restrictive as the conditions mentioned at p.740
!! of the international tables for crystallography (1983).
!! The final vectors form a right-handed basis, while their
!! sign and ordering is chosen such as to maximize the overlap
!! with the original vectors in order.
!!
!! INPUTS
!!  rprimd(3,3)=primitive vectors
!!
!! OUTPUT
!!  metmin(3,3)=metric for the new (minimal) primitive vectors
!!  minim(3,3)=minimal primitive translations
!!
!! NOTES
!! The routine might as well be defined without
!! metmin as argument, but it is more convenient to have it
!!
!! PARENTS
!!      m_kpts,m_symfind
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE


subroutine smallprim(metmin,minim,rprimd)

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(out) :: metmin(3,3),minim(3,3)

!Local variables-------------------------------
!scalars
 integer :: ia,ib,ii,itrial,minimal
 integer :: iiter, maxiter = 100000
 real(dp) :: determinant,length2,metsum
 character(len=500) :: msg
!arrays
 integer :: nvecta(3),nvectb(3)
 real(dp) :: rmet(3,3),scprod(3),tmpvect(3)

!**************************************************************************

 !call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 rmet = MATMUL(TRANSPOSE(rprimd),rprimd)

 nvecta(1)=2 ; nvectb(1)=3
 nvecta(2)=1 ; nvectb(2)=3
 nvecta(3)=1 ; nvectb(3)=2

 minim(:,:)=rprimd(:,:)
 metmin(:,:)=rmet(:,:)

!DEBUG
!write(std_out,*)' smallprim : starting values, rprim '
!write(std_out,'(3f16.8)' )rprimd(:,1)
!write(std_out,'(3f16.8)' )rprimd(:,2)
!write(std_out,'(3f16.8)' )rprimd(:,3)
!write(std_out,*)' smallprim : starting values, rmet '
!write(std_out,'(3f16.8)' )rmet(:,1)
!write(std_out,'(3f16.8)' )rmet(:,2)
!write(std_out,'(3f16.8)' )rmet(:,3)
!ENDDEBUG

!Note this loop without index
 do iiter = 1, maxiter

!  Will exit if minimal=1 is still valid after a trial
!  to reduce the vectors of each of the three pairs
   minimal=1

   do itrial=1,3

     ia=nvecta(itrial) ; ib=nvectb(itrial)
!    Make sure the scalar product is negative
     if(metmin(ia,ib)>tol8)then
       minim(:,ia)=-minim(:,ia)
       metmin(ia,ib)=-metmin(ia,ib) ; metmin(ib,ia)=-metmin(ib,ia)
       metmin(ia,itrial)=-metmin(ia,itrial)
       metmin(itrial,ia)=-metmin(itrial,ia)
     end if
!    Compute the length of the sum vector
     length2=metmin(ia,ia)+2*metmin(ia,ib)+metmin(ib,ib)
!    Replace the first vector by the sum vector if the latter is smaller
     if(length2/metmin(ia,ia) < one-tol8)then
       minim(:,ia)=minim(:,ia)+minim(:,ib)
       metmin(ia,ia)=length2
       metmin(ia,ib)=metmin(ia,ib)+metmin(ib,ib)
       metmin(ia,itrial)=metmin(ia,itrial)+metmin(ib,itrial)
       metmin(ib,ia)=metmin(ia,ib)
       metmin(itrial,ia)=metmin(ia,itrial)
       minimal=0
!      Replace the second vector by the sum vector if the latter is smaller
     else if(length2/metmin(ib,ib) < one-tol8)then
       minim(:,ib)=minim(:,ia)+minim(:,ib)
       metmin(ib,ib)=length2
       metmin(ia,ib)=metmin(ia,ib)+metmin(ia,ia)
       metmin(itrial,ib)=metmin(itrial,ib)+metmin(itrial,ia)
       metmin(ib,ia)=metmin(ia,ib)
       metmin(ib,itrial)=metmin(itrial,ib)
       minimal=0
     end if

   end do

   if(minimal==1)exit
 end do

 if (iiter >= maxiter) then
   write(msg,'(a,i0,a)') 'the loop has failed to find a set of minimal vectors in ',maxiter,' iterations.'
   MSG_BUG(msg)
 end if

!At this stage, the three vectors have angles between each other that are
!comprised between 90 and 120 degrees. It might still be that minus the vector
!that is the sum of the three vectors is smaller than the longest of these vectors
 do iiter = 1, maxiter

!  Will exit if minimal=1 is still valid after a trial
!  to replace each of the three vectors by minus the summ of the three vectors
   minimal=1
   metsum=sum(metmin(:,:))
   do itrial=1,3
     ia=nvecta(itrial) ; ib=nvectb(itrial)
     if(metmin(ia,ia)/metsum > one + tol8)then
       minim(:,ia)=-minim(:,1)-minim(:,2)-minim(:,3)
       metmin(ia,ib)=-sum(metmin(:,ib))
       metmin(ia,itrial)=-sum(metmin(:,itrial))
       metmin(ia,ia)=metsum
       metmin(ib,ia)=metmin(ia,ib)
       metmin(itrial,ia)=metmin(ia,itrial)
       minimal=0
     end if
   end do

   if(minimal==1)exit

 end do

 if (iiter >= maxiter) then
   write(msg, '(a,i0,a)') 'the second loop has failed to find a set of minimal vectors in ',maxiter, 'iterations.'
   MSG_BUG(msg)
 end if

!DEBUG
!write(std_out,'(a,3es14.6,a,3es14.6,a,3es14.6)')' rprimd=',rprimd(:,1),ch10,rprimd(:,2),ch10,rprimd(:,3)
!write(std_out,'(a,3es14.6,a,3es14.6,a,3es14.6)')' minim =',minim(:,1),ch10,minim(:,2),ch10,minim(:,3)
!ENDDEBUG

!DEBUG
!Change sign of the third vector if not right-handed basis
!determinant=minim(1,1)*(minim(2,2)*minim(3,3)-minim(3,2)*minim(2,3))+&
!&            minim(2,1)*(minim(3,2)*minim(1,3)-minim(1,2)*minim(3,3))+&
!&            minim(3,1)*(minim(1,2)*minim(2,3)-minim(2,2)*minim(1,3))
!write(std_out,*)' smallprim: determinant=',determinant
!ENDDEBUG

!Choose the first vector
!Compute the scalar product of the three minimal vectors with the first original vector
 scprod(:)=zero
 do ii=1,3
   scprod(:)=scprod(:)+minim(ii,:)*rprimd(ii,1)
 end do
!Determine the vector with the maximal absolute overlap
 itrial=1
 if(abs(scprod(2))>abs(scprod(1))+tol8)itrial=2
 if(abs(scprod(3))>abs(scprod(itrial))+tol8)itrial=3
!Switch the vectors if needed
 if(itrial/=1)then
   tmpvect(:)=minim(:,1)
   minim(:,1)=minim(:,itrial)
   minim(:,itrial)=tmpvect(:)
 end if
!Choose the sign
 if(scprod(itrial)<tol8)minim(:,1)=-minim(:,1)

!DEBUG
!Change sign of the third vector if not right-handed basis
!determinant=minim(1,1)*(minim(2,2)*minim(3,3)-minim(3,2)*minim(2,3))+&
!&            minim(2,1)*(minim(3,2)*minim(1,3)-minim(1,2)*minim(3,3))+&
!&            minim(3,1)*(minim(1,2)*minim(2,3)-minim(2,2)*minim(1,3))
!write(std_out,*)' smallprim: determinant=',determinant
!ENDDEBUG

!Choose the second vector
!Compute the scalar product of the second and third minimal vectors with the second original vector
 scprod(2:3)=zero
 do ii=1,3
   scprod(2:3)=scprod(2:3)+minim(ii,2:3)*rprimd(ii,2)
 end do
!Determine the vector with the maximal absolute overlap
 itrial=2
 if(abs(scprod(3))>abs(scprod(2))+tol8)itrial=3
!Switch the vectors if needed
 if(itrial/=2)then
   tmpvect(:)=minim(:,2)
   minim(:,2)=minim(:,itrial)
   minim(:,itrial)=tmpvect(:)
 end if
!Choose the sign
 if(scprod(itrial)<tol8)minim(:,2)=-minim(:,2)

!Change sign of the third vector if not right-handed basis
 determinant=minim(1,1)*(minim(2,2)*minim(3,3)-minim(3,2)*minim(2,3))+&
& minim(2,1)*(minim(3,2)*minim(1,3)-minim(1,2)*minim(3,3))+&
& minim(3,1)*(minim(1,2)*minim(2,3)-minim(2,2)*minim(1,3))
 if(determinant<-tol8)minim(:,3)=-minim(:,3)
 if(abs(determinant)<tol8)then
   MSG_BUG('minim gives vanishing unit cell volume.')
 end if

!Final computation of metmin
 do ii=1,3
   metmin(ii,:)=minim(1,ii)*minim(1,:)+ minim(2,ii)*minim(2,:)+ minim(3,ii)*minim(3,:)
 end do

!DEBUG
!write(std_out,'(a,3es14.6,a,3es14.6,a,3es14.6)')' rprimd=',rprimd(:,1),ch10,rprimd(:,2),ch10,rprimd(:,3)
!write(std_out,'(a,3es16.8,a,3es16.8,a,3es16.8)')' minim =',minim(:,1),ch10,minim(:,2),ch10,minim(:,3)
!write(std_out,'(a,3es16.8,a,3es16.8,a,3es16.8)')' metmin =',metmin(:,1),ch10,metmin(:,2),ch10,metmin(:,3)
!ENDDEBUG

end subroutine smallprim
!!***

!!****f* m_symtk/print_symmetries
!! NAME
!! print_symmetries
!!
!! FUNCTION
!!  Helper function to print the set of symmetries.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_crystal,m_hdr,m_spgbuilder
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine print_symmetries(nsym, symrel, tnons, symafm, unit, mode_paral)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,optional,intent(in) :: unit
 character(len=4),optional,intent(in) :: mode_paral
!arrays
 integer,intent(in) :: symrel(3,3,nsym),symafm(nsym)
 real(dp),intent(in) :: tnons(3,nsym)

!Local variables-------------------------------
 integer :: my_unt,isym,isymin,isymend,ii,jj
 character(len=500) :: msg
 character(len=4) :: my_mode
! *********************************************************************

 my_unt =std_out; if (PRESENT(unit      )) my_unt =unit
 my_mode='COLL' ; if (PRESENT(mode_paral)) my_mode=mode_paral

 !write(msg,'(2a)')ch10,' Rotations                           Translations     Symafm '
 !do isym=1,nsym
 ! write(msg,'(1x,3(3i3,1x),4x,3(f11.7,1x),6x,i2)')symrel(:,:,isym),tnons(:,isym),symafm(isym)
 ! call wrtout(my_unt,msg,my_mode)
 !end do

 write(msg,'(2a)')ch10,' Symmetry operations in real space (Rotation tnons AFM)'
 call wrtout(my_unt,msg,my_mode)
 do isymin=1,nsym,4
   isymend=isymin+3
   if (isymend>nsym) isymend=nsym
   do ii=1,3
     write(msg,'(4(3i3,f8.3,i3,3x))')((symrel(ii,jj,isym),jj=1,3),tnons(ii,isym),symafm(isym),isym=isymin,isymend)
     call wrtout(my_unt,msg,my_mode)
   end do
   write(msg,'(a)')ch10
   call wrtout(my_unt,msg,my_mode)
 end do

end subroutine print_symmetries
!!***

end module m_symtk
!!***
