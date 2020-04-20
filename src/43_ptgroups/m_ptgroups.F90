!!****m* ABINIT/m_ptgroups
!! NAME
!! m_ptgroups
!!
!! FUNCTION
!!  This module contains the irreducible representations and the
!!  character tables of the 32 point groups.
!!
!! COPYRIGHT
!! Copyright (C) 2010-2020 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_ptgroups

 use defs_basis
 use m_defs_ptgroups
 use m_errors
 use m_abicore

 use m_io_tools,       only : open_file
 use m_fstrings,       only : sjoin
 use m_numeric_tools,  only : get_trace, cmplx_sphcart
 use m_symtk,          only : mati3inv

! Import group tables
 use m_ptg_C1
 use m_ptg_Ci
 use m_ptg_C2
 use m_ptg_Cs
 use m_ptg_C2h
 use m_ptg_D2
 use m_ptg_C2v
 use m_ptg_D2h
 use m_ptg_C4
 use m_ptg_S4
 use m_ptg_C4h
 use m_ptg_D4
 use m_ptg_C4v
 use m_ptg_D2d
 use m_ptg_D4h
 use m_ptg_C3
 use m_ptg_C3i
 use m_ptg_D3
 use m_ptg_C3v
 use m_ptg_D3d
 use m_ptg_C6
 use m_ptg_C3h
 use m_ptg_C6h
 use m_ptg_D6
 use m_ptg_C6v
 use m_ptg_D3h
 use m_ptg_D6h
 use m_ptg_T
 use m_ptg_Th
 use m_ptg_O
 use m_ptg_Td
 use m_ptg_Oh

 implicit none

 private

 public :: get_point_group
 public :: get_classes
 public :: show_character_tables
 public :: point_group_free
 public :: point_group_init
 public :: point_group_print
 public :: locate_sym
 public :: irrep_free
 public :: copy_irrep
 public :: init_irrep
 public :: sum_irreps
 public :: mult_table
 !public :: polish_irreps

 interface irrep_free
   module procedure irrep_free_0d
   module procedure irrep_free_1d
 end interface irrep_free

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_ptgroups/get_point_group
!! NAME
!!  get_point_group
!!
!! FUNCTION
!!
!! INPUTS
!! ptg_name=point group name as returned by symptgroup.
!!
!! OUTPUT
!! nsym=Number of symmetries in the point group.
!! nclass=Number of classes.
!! sym(3,3,nsym)=Elements of the point group ordered by classe.
!! class_ids(2,nclass)=Initial and final index in sym, for each
!! Irreps(nclass)=Datatype gathering data on the different irreducible representations.
!!
!! PARENTS
!!      m_ptgroups
!!
!! CHILDREN
!!      irrep_free
!!
!! SOURCE

subroutine get_point_group(ptg_name,nsym,nclass,sym,class_ids,class_names,Irreps)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: nclass,nsym
 character(len=*),intent(in) :: ptg_name
!arrays
 integer,allocatable,intent(out) :: sym(:,:,:),class_ids(:,:)
 character(len=5),allocatable,intent(out) :: class_names(:)
 type(irrep_t),allocatable,intent(out) :: Irreps(:)

!Local variables-------------------------------
!scalars
 integer :: irp,isym
 !character(len=500) :: msg

! *************************************************************************

 SELECT CASE (TRIM(ADJUSTL(ptg_name)))
 CASE ('1')
   call ptg_C1  (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('-1')
   call ptg_Ci  (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('2')
   call ptg_C2  (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('m',"-2") ! Abinit uses "-2"
   call ptg_Cs  (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('2/m')
   call ptg_C2h (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('222')
   call ptg_D2  (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('mm2')
   call ptg_C2v (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('mmm')
   call ptg_D2h (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('4')
   call ptg_C4  (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('-4')
   call ptg_S4  (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('4/m')
   call ptg_C4h (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('422')
   call ptg_D4  (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('4mm')
   call ptg_C4v (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('-42m')
   call ptg_D2d (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('4/mmm')
   call ptg_D4h (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('3')
   call ptg_C3  (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('-3')
   call ptg_C3i (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('32')
   call ptg_D3  (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('3m')
   call ptg_C3v (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('-3m')
   call ptg_D3d (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('6')
   call ptg_C6  (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('-6')
   call ptg_C3h (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('6/m')
   call ptg_C6h (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('622')
   call ptg_D6  (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('6mm')
   call ptg_C6v (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('-62m')
   call ptg_D3h (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('6/mmm')
   call ptg_D6h (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('23')
   call ptg_T   (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('m-3')
   call ptg_Th  (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('432')
   call ptg_O   (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('-43m')
   call ptg_Td  (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE ('m-3m')
   call ptg_Oh  (nsym,nclass,sym,class_ids,class_names,Irreps)
 CASE DEFAULT
   MSG_BUG(sjoin("Unknown value for ptg_name:", ptg_name))
 END SELECT

 ! Calculate the trace of each irreducible representation in order to have the character at hand.
 do irp=1,SIZE(Irreps)
   ABI_MALLOC(Irreps(irp)%trace, (nsym))
   do isym=1,nsym
     Irreps(irp)%trace(isym) = get_trace(Irreps(irp)%mat(:,:,isym))
   end do
 end do

end subroutine get_point_group
!!***

!----------------------------------------------------------------------

!!****f* m_ptgroups/get_classes
!! NAME
!! get_classes
!!
!! FUNCTION
!!  Given a set of nsym 3x3 operations in reciprocal or real space,
!!  which are supposed to form a group, this routine divides the group into classes.
!!
!! INPUTS
!! nsym=number of symmetry operation.
!! sym(3,3,nsym)=The symmetry operations.
!!
!! OUTPUT
!! nclass=The number of classes
!! nelements(1:nclass)=For each class, the number of elements
!! elements_idx(ii,1:nclass)=For each class, this table gives the index
!!   of its elements (ii=1,..,nelements(iclass))
!!
!! NOTES
!!  * A class is defined as the set of distinct elements obtained by
!!    considering for each element, S, of the group all its conjugate
!!    elements X^-1 S X where X range over all the elements of the group.
!!
!!  * It does not work in case of non-collinear magnetism.
!!
!!  * The routine assumes that anti-ferromagnetic symmetries (if any) have been removed by the caller.
!!
!! PARENTS
!!      m_esymm,m_ptgroups
!!
!! CHILDREN
!!      irrep_free
!!
!! SOURCE

subroutine get_classes(nsym,sym,nclass,nelements,elements_idx)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(out) :: nclass
!arrays
 integer,intent(in) :: sym(3,3,nsym)
 integer,intent(out) :: nelements(nsym),elements_idx(nsym,nsym)

!Local variables-------------------------------
!scalars
 integer :: isym,jsym,ksym,identity_idx !,ierr
 character(len=500) :: msg
!arrays
 integer :: cjg(3,3),ss(3,3),xx(3,3),xxm1(3,3),test(3,3)
 integer :: identity(3,3)
 logical :: found(nsym),found_identity

!************************************************************************

 ! === Check if identity is present in the first position ===
 identity=RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/)); found_identity=.FALSE.

 do isym=1,nsym
   if (ALL(sym(:,:,isym)==identity)) then
     found_identity=.TRUE.; identity_idx=isym; EXIT
   end if
 end do
 if (.not.found_identity.or.identity_idx/=1) then
  write(msg,'(3a)')&
&  'Either identity is not present or it is not the first operation ',ch10,&
&  'check set of symmetry operations '
  MSG_ERROR(msg)
 end if
 !
 ! Is it a group? Note that I assume that AFM sym.op (if any) have been pruned in the caller.
 !dummy_symafm=1
 !call chkgrp(nsym,dummy_symafm,sym,ierr)
 !ABI_CHECK(ierr==0,"Error in group closure")

 nclass=0; nelements(:)=0; elements_idx(:,:)=0; found(:)=.FALSE.
 do isym=1,nsym
   if (.not.found(isym)) then
     nclass=nclass+1
     ss(:,:)=sym(:,:,isym)

     do jsym=1,nsym ! Form conjugate.
       xx(:,:)=sym(:,:,jsym)
       call mati3inv(xx,xxm1) ; xxm1=TRANSPOSE(xxm1)
       cjg(:,:)=MATMUL(xxm1,MATMUL(ss,xx))
       do ksym=1,nsym ! Is it already found?
         test(:,:)=sym(:,:,ksym)
         if (.not.found(ksym).and.(ALL((test-cjg)==0))) then
           found(ksym)=.TRUE.
           nelements(nclass)=nelements(nclass)+1
           elements_idx(nelements(nclass),nclass)=ksym
         end if
       end do
     end do

   end if
 end do

end subroutine get_classes
!!***

!----------------------------------------------------------------------

!!****f* m_ptgroups/show_character_tables
!! NAME
!!  show_character_tables
!!
!! FUNCTION
!!   Printout of the caracter tables of the 32 point groups.
!!
!! INPUTS
!!  [unit]=Unit number of output file. Defaults to std_out
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!
!! CHILDREN
!!      irrep_free
!!
!! SOURCE

subroutine show_character_tables(unit)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit

!Local variables-------------------------------
 integer :: igrp,my_unt
 character(len=5) :: ptg_name
 type(point_group_t) :: Ptg
!arrays
 !integer,allocatable :: elements_idx(:,:),nelements(:)

! *********************************************************************

 my_unt = std_out; if (PRESENT(unit)) my_unt=unit

 do igrp=1,SIZE(ptgroup_names)
   ptg_name = ptgroup_names(igrp)
   call point_group_init(Ptg,ptg_name)
   call point_group_print(Ptg,unit=my_unt)
   !allocate(nelements(Ptg%nsym),elements_idx(Ptg%nsym,Ptg%nsym))
   !call get_classes(Ptg%nsym,Ptg%sym,nclass,nelements,elements_idx)
   !deallocate(nelements,elements_idx)
   call point_group_free(Ptg)
 end do

end subroutine show_character_tables
!!***

!----------------------------------------------------------------------

!!****f* m_ptgroups/point_group_free
!! NAME
!! point_group_free
!!
!! FUNCTION
!!  Deallocate all memory allocated in the point_group_t datatype.
!!
!! PARENTS
!!      m_ptgroups
!!
!! CHILDREN
!!      irrep_free
!!
!! SOURCE

subroutine point_group_free(Ptg)

 implicit none

!Arguments ------------------------------------
 type(point_group_t),intent(inout) :: Ptg

!Local variables-------------------------------
! *********************************************************************

 !@point_group_t
 if (allocated(Ptg%class_ids))   then
   ABI_FREE(Ptg%class_ids)
 end if
 if (allocated(Ptg%sym)) then
   ABI_FREE(Ptg%sym)
 end if
 if (allocated(Ptg%class_names)) then
   ABI_FREE(Ptg%class_names)
 end if

 if (allocated(Ptg%Irreps)) then
   call irrep_free(Ptg%Irreps)
   ABI_DT_FREE(Ptg%Irreps)
 end if

end subroutine point_group_free
!!***

!----------------------------------------------------------------------

!!****f* m_ptgroups/point_group_init
!! NAME
!! point_group_init
!!
!! FUNCTION
!!  Creation method for the point_group_t datatype.
!!
!! INPUTS
!!  ptg_name=The name of the point group (International conventions).
!!
!! OUTPUT
!!  The datatype completely initialized.
!!
!! PARENTS
!!      m_esymm,m_ptgroups
!!
!! CHILDREN
!!      irrep_free
!!
!! SOURCE


subroutine point_group_init(Ptg,ptg_name)

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=5),intent(in) :: ptg_name
 type(point_group_t),intent(inout) :: Ptg
! *********************************************************************

 !@point_group_t
 !call wrtout(std_out," Retrieving point group data for: "//TRIM(ptg_name),"COLL")

 Ptg%gname = ptg_name
 call get_point_group(Ptg%gname,Ptg%nsym,Ptg%nclass,Ptg%sym,Ptg%class_ids,Ptg%class_names,Ptg%Irreps)

end subroutine point_group_init
!!***

!----------------------------------------------------------------------

!!****f* m_ptgroups/point_group_print
!! NAME
!! point_group_print
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_ptgroups
!!
!! CHILDREN
!!      irrep_free
!!
!! SOURCE

subroutine point_group_print(Ptg,header,unit,mode_paral,prtvol)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header
 type(point_group_t),target,intent(in) :: Ptg

!Local variables-------------------------------
 integer :: my_unt,my_prtvol,irp,icls,sidx
 complex(dpc) :: trace
 character(len=4) :: my_mode
 character(len=500) :: msg
 type(irrep_t),pointer :: Row
! *********************************************************************

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=' ==== Point Group Table ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 write(std_out,*)REPEAT("=",80)
 write(std_out,*)" Point group : ",TRIM(Ptg%gname)," Number of symmetries ",Ptg%nsym," Number of classes    ",Ptg%nclass

 write(std_out,"(a6)",advance="no")"Class "
 do icls=1,Ptg%nclass
   write(std_out,"('|',a10)",advance="no")Ptg%class_names(icls)
 end do
 write(std_out,"('|')",advance="no")
 write(std_out,*)" "

 write(std_out,"(a6)",advance="no")"Mult  "
 do icls=1,Ptg%nclass
   write(std_out,"('|',i10)",advance="no")Ptg%class_ids(2,icls)-Ptg%class_ids(1,icls) + 1
 end do
 write(std_out,"('|')",advance="no")
 write(std_out,*)" "

 do irp=1,SIZE(Ptg%Irreps)
   Row =>  Ptg%Irreps(irp)
   write(std_out,'(a6)',advance="no")TRIM(Row%name)

   do icls=1,Ptg%nclass
     sidx = Ptg%class_ids(1,icls)
     trace = Row%trace(sidx)
     if (ABS(AIMAG(trace)) > tol6) then
        write(std_out,"('|',(2f5.2))",advance="no")trace
      else
        write(std_out,"('|',(f10.2))",advance="no")REAL(trace)
      end if
   end do

   write(std_out,"('|')",advance="no")
   write(std_out,*)" "
 end do

 write(std_out,*)REPEAT("=",80)

end subroutine point_group_print
!!***

!----------------------------------------------------------------------

!!****f* m_ptgroups/locate_sym
!! NAME
!!  locate_sym
!!
!! FUNCTION
!!  Given a symmetry operation asym, this routine returns its index in the Ptg%sym
!!  array and the index of the class it belongs to.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_esymm
!!
!! CHILDREN
!!      irrep_free
!!
!! SOURCE

subroutine locate_sym(Ptg,asym,sym_idx,cls_idx,ierr)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: sym_idx,cls_idx
 integer,optional,intent(out) :: ierr
 type(point_group_t),intent(in) :: Ptg
!arrays
 integer,intent(in) :: asym(3,3)

!Local variables-------------------------------
!scalars
 integer :: isym,icls
 character(len=500) :: msg

! *********************************************************************

 sym_idx = 0
 do isym=1,Ptg%nsym
   if (ALL(asym == Ptg%sym(:,:,isym) )) then
     sym_idx = isym
     EXIT
   end if
 end do

 cls_idx = 0
 do icls=1,Ptg%nclass
   if (sym_idx >= Ptg%class_ids(1,icls) .and. &
&      sym_idx <= Ptg%class_ids(2,icls) ) then
     cls_idx = icls
     EXIT
   end if
 end do

 if (PRESENT(ierr)) ierr=0
 if (sym_idx==0 .or. cls_idx==0) then
   write(msg,'(a,9(i0,1x),3a,i1,a,i1)')&
&    " Symmetry: ",asym," not found in point group table ",ch10,&
&    " sym_idx= ",sym_idx, " and cls_idx= ",cls_idx
   if (PRESENT(ierr)) then
     ierr=1
     MSG_WARNING(msg)
   else
     MSG_ERROR(msg)
   end if
 end if

end subroutine locate_sym
!!***

!----------------------------------------------------------------------

!!****f* m_ptgroups/mult_table
!! NAME
!! mult_table
!!
!! FUNCTION
!!  Given a set of nsym 3x3 operations which are supposed to form a group,
!!  this routine constructs the multiplication table of the group.
!!
!! INPUTS
!! nsym=number of symmetry operation
!! sym(3,3,nsym)=the operations
!!
!! OUTPUT
!!  mtab(nsym,nsym)=The index of the product S_i * S_j in the input set sym.
!!
!! PARENTS
!!
!! CHILDREN
!!      irrep_free
!!
!! SOURCE

subroutine mult_table(nsym,sym,mtab)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
!arrays
 integer,intent(in) :: sym(3,3,nsym)
 integer,intent(out) :: mtab(nsym,nsym)

!Local variables-------------------------------
!scalars
 integer :: isym,jsym,ksym
 !character(len=500) :: msg
!arrays
 integer :: prod_ij(3,3),found(nsym)

!************************************************************************

 do jsym=1,nsym
   found(:)=0 ! Each symmetry should compare only once in a given (row|col).

   do isym=1,nsym
     prod_ij = MATMUL(sym(:,:,isym),sym(:,:,jsym))
     do ksym=1,nsym
       if ( ALL(prod_ij == sym(:,:,ksym)) ) then
         found(ksym)=found(ksym)+1
         mtab(isym,jsym) = ksym
       end if
     end do
   end do ! jsym

   if ( ANY(found /= 1)) then
     write(std_out,*)"found = ",found
     MSG_ERROR("Input elements do not form a group")
   end if
 end do ! isym

end subroutine mult_table
!!***

!----------------------------------------------------------------------

!!****f* m_ptgroups/groupk_from_file
!! NAME
!!  groupk_from_file
!!
!! FUNCTION
!!  Initialize the group_k_t datatype from an external database retrieved from
!!  the Bilbao server via the ptg.py script.
!!
!! INPUTS
!!  fname(len=*)=file name
!!
!! OUTPUT
!!  ierr=Status error
!!  Lgrps<group_k_t>=The structure completely initialized.
!!
!! TODO
!!   This is a stub. I still have to complete the fileformat for the Bilbao database.
!!
!! PARENTS
!!
!! CHILDREN
!!      irrep_free
!!
!! SOURCE

subroutine groupk_from_file(Lgrps,spgroup,fname,nkpt,klist,ierr)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spgroup
 integer,intent(out) :: ierr,nkpt
 character(len=fnlen),intent(in) :: fname
!arrays
 type(group_k_t),target,allocatable :: Lgrps(:)
 real(dp),pointer :: klist(:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: last_file_version=1
 integer :: unt,fvers,ik,nsym_ltgk,isym,nirreps_k,irp,icls
 integer :: irrep_idx,irrep_dim,sym_idx,ita_spgroup,nels,prev,now
 character(len=IRREPNAME_LEN) :: irrep_name
 character(len=1) :: basis
 character(len=500) :: msg
 type(group_k_t),pointer :: Gk
!arrays
 integer,allocatable :: nelements(:),elements_idx(:,:)
 real(dp) :: kpt(3)
 character(len=10),allocatable :: kname(:)
 type(irrep_t),pointer :: OneIrr

! *************************************************************************

 ierr=0
 if (open_file(fname,msg,newunit=unt,form="formatted") /=0) then
   MSG_ERROR(msg)
 end if

 read(unt,*,ERR=10)          ! Skip the header.
 read(unt,*,ERR=10) fvers    ! File version.
 if (fvers > last_file_version) then
   write(msg,"(2(a,i0))")" Found file format= ",fvers," but the latest supported version is: ",last_file_version
   MSG_ERROR(msg)
 end if
 read(unt,*,ERR=10) ita_spgroup
 read(unt,*,ERR=10) basis

 if (spgroup/=ita_spgroup) then
   write(msg,'(a,2i0)')&
&   " Input space group does not match with the value reported on file: ",spgroup,ita_spgroup
   MSG_ERROR(msg)
 end if

 if (basis /= "b") then
   msg=" Wrong value for basis: "//TRIM(basis)
   MSG_ERROR(msg)
 end if

 ! * Read the list of the k-points.
 read(unt,*,ERR=10) nkpt

 ABI_DT_MALLOC(Lgrps,(nkpt))

 ABI_MALLOC(klist,(3,nkpt))
 ABI_MALLOC(kname,(nkpt))
 do ik=1,nkpt
   read(unt,*,ERR=10) klist(:,ik), kname(ik)
 end do

 ! * Read tables for each k-point
 do ik=1,nkpt

   read(unt,*,ERR=10) kpt
   read(unt,*,ERR=10) nsym_ltgk
   Gk  => Lgrps(ik)

   Gk%spgroup = ita_spgroup
   Gk%nsym    = nsym_ltgk
   Gk%point   = kpt
   ABI_MALLOC(Gk%sym,(3,3,nsym_ltgk))
   ABI_MALLOC(Gk%tnons,(3,nsym_ltgk))

   do isym=1,nsym_ltgk ! Read symmetries of the little group.
     read(unt,*,ERR=10) Gk%sym(:,:,isym)
     read(unt,*,ERR=10) Gk%tnons(:,isym)
   end do

   ABI_MALLOC(nelements,(nsym_ltgk))
   ABI_MALLOC(elements_idx,(nsym_ltgk,nsym_ltgk))

   call get_classes(nsym_ltgk,Gk%sym,Gk%nclass,nelements,elements_idx)

   ! The operations reported on the file are supposed to be packed in classes
   ! otherwise one should perform a rearrangement of the indices.
   prev = 0
   do icls=1,Gk%nclass
     do isym=1,nelements(icls)
       now = elements_idx(isym,icls)
       if ( (now-prev) /= 1 ) then
         write(msg,"(2(a,i0))")&
&          " Symmetries on file are not ordered in classes. icls= ",icls,", isym= ",isym
         MSG_ERROR(msg)
       else
         prev = now
       end if
     end do
   end do

   ABI_MALLOC(Gk%class_ids,(2,Gk%nclass))
   do icls=1,Gk%nclass
     nels = nelements(icls)
     Gk%class_ids(1,icls) = elements_idx(1,   icls)
     Gk%class_ids(2,icls) = elements_idx(nels,icls)
   end do

   ABI_FREE(nelements)
   ABI_FREE(elements_idx)

   ! Read the irreducible representations.
   read(unt,*,ERR=10) nirreps_k
   ABI_CHECK(Gk%nclass == nirreps_k,"Gk%nclass /= nirreps_k")

   !$$ allocate(Gk%class_names(Gk%nclass))
   ABI_DT_MALLOC(Gk%Irreps,(nirreps_k))

   do irp=1,nirreps_k
     OneIrr =>  Gk%Irreps(irp)
     read(unt,*,ERR=10) irrep_idx, irrep_dim, irrep_name
     call init_irrep(OneIrr,nsym_ltgk,irrep_dim,irrep_name)
     do isym=1,nsym_ltgk
       read(unt,*,ERR=10) sym_idx, OneIrr%mat(:,:,isym)
       ABI_CHECK(sym_idx==irp,"sym_idx/=irp!")
       ! Matrix elements on file are in the form (rho, theta) with theta given in degrees.
       call cmplx_sphcart(OneIrr%mat(:,:,isym),from="Sphere",units="Degrees")
       OneIrr%trace(isym) = get_trace(OneIrr%mat(:,:,isym))
     end do
   end do

 end do

 close(unt)
 RETURN
 !
 ! Handle IO-error.
10 ierr=1
 close(unt)
 RETURN

end subroutine groupk_from_file
!!***

!----------------------------------------------------------------------

!!****f* m_ptgroups/irrep_free_0d
!! NAME
!! irrep_free_0d
!!
!! FUNCTION
!!  Deallocate all memory allocated in the irrep_t datatype.
!!
!! PARENTS
!!      m_ptgroups
!!
!! CHILDREN
!!      irrep_free
!!
!! SOURCE

subroutine irrep_free_0d(Irrep)

 implicit none

!Arguments ------------------------------------
 type(irrep_t),intent(inout) :: Irrep

! *********************************************************************

 !@irrep_t
 if (allocated(Irrep%trace))  then
   ABI_FREE(Irrep%trace)
 end if
 if (allocated(Irrep%mat))  then
   ABI_FREE(Irrep%mat)
 end if

end subroutine irrep_free_0d
!!***

!----------------------------------------------------------------------

!!****f* m_ptgroups/irrep_free_1d
!! NAME
!! irrep_free_1d
!!
!! FUNCTION
!!  Deallocate all memory allocated in the irrep_t datatype.
!!
!! PARENTS
!!
!! CHILDREN
!!      irrep_free
!!
!! SOURCE

subroutine irrep_free_1d(Irrep)

 implicit none

!Arguments ------------------------------------
 type(irrep_t),intent(inout) :: Irrep(:)

!Local variables-------------------------------
 integer :: irp
! *********************************************************************

 do irp=1,SIZE(Irrep)
   call irrep_free_0d(Irrep(irp))
 end do

end subroutine irrep_free_1d
!!***

!----------------------------------------------------------------------

!!****f* m_ptgroups/copy_irrep
!! NAME
!!  copy_irrep
!!
!! FUNCTION
!!  Perform a copy of a set of irrep_t datatypes. Optionally one can multiply
!!  by a phase factor.
!!
!! PARENTS
!!      m_esymm
!!
!! CHILDREN
!!      irrep_free
!!
!! SOURCE

subroutine copy_irrep(In_irreps,Out_irreps,phase_fact)

 implicit none

!Arguments ------------------------------------
 type(irrep_t),intent(in) :: In_irreps(:)
 type(irrep_t),intent(inout) :: Out_irreps(:)
 complex(dpc),optional,intent(in) :: phase_fact(:)

!Local variables-------------------------------
!scalars
 integer :: irp,dim1,dim2,in_nsym,in_dim,isym
 character(len=500) :: msg
!arrays
 complex(dpc) :: my_phase_fact(In_irreps(1)%nsym)
! *********************************************************************

 !@irrep_t
 dim1 = SIZE( In_irreps)
 dim2 = SIZE(Out_irreps)
 if (dim1/=dim2) then
   msg = " irreps to be copied have different dimension"
   MSG_ERROR(msg)
 end if

 my_phase_fact=cone
 if (PRESENT(phase_fact)) then
   my_phase_fact=phase_fact
   if (SIZE(phase_fact) /= In_irreps(1)%nsym) then
     msg = " irreps to be copied have different dimension"
     MSG_ERROR(msg)
   end if
 end if

 do irp=1,dim1
   in_dim  = In_irreps(irp)%dim
   in_nsym = In_irreps(irp)%nsym
   call init_irrep(Out_irreps(irp),in_nsym,in_dim)
   Out_irreps(irp)%name = In_irreps(irp)%name
   do isym=1,in_nsym
     Out_irreps(irp)%mat(:,:,isym) = In_irreps(irp)%mat(:,:,isym) * my_phase_fact(isym)
     Out_irreps(irp)%trace(isym) = get_trace(Out_irreps(irp)%mat(:,:,isym))
   end do
 end do

end subroutine copy_irrep
!!***

!----------------------------------------------------------------------

!!****f* m_ptgroups/init_irrep
!! NAME
!!  alloc_irrep
!!
!! FUNCTION
!!  Initialize an instance of the irrep_t datatype.
!!
!! INPUTS
!!  nsym=The number of symmetries.
!!  irr_dim=The dimension of the irrep.
!!  [irr_name]=The name of theirrep. "???" is used if not given
!!
!! OUTPUT
!!  Irrep<irrep_t>=
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_irrep(Irrep,nsym,irr_dim,irr_name)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
!arrays
 integer,intent(in) :: irr_dim
 character(len=*),optional,intent(in) :: irr_name
 type(irrep_t),intent(inout) :: Irrep

!Local variables-------------------------------
 !character(len=500) :: msg
! *********************************************************************

 !@irrep_t
 Irrep%dim      = irr_dim
 Irrep%nsym     = nsym
 Irrep%name     = "???"
 if (PRESENT(irr_name)) Irrep%name = irr_name

 ABI_MALLOC(Irrep%mat,(irr_dim,irr_dim,nsym))
 Irrep%mat = czero

 ABI_MALLOC(Irrep%trace,(nsym))
 Irrep%trace = czero

end subroutine init_irrep
!!***

!----------------------------------------------------------------------

!!****f* m_ptgroups/sum_irreps
!! NAME
!!  sum_irreps
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

function sum_irreps(Irrep1,Irrep2,ii,jj,kk,ll) result(res)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ii,jj,kk,ll
!arrays
 type(irrep_t),intent(in) :: Irrep1,Irrep2
 complex(dpc) :: res

!Local variables-------------------------------
 integer :: isym,nsym,ierr
 !character(len=500) :: msg
! *********************************************************************

 ierr=0; res = czero

 nsym = Irrep1%nsym
 if (nsym /= Irrep2%nsym) then
   MSG_WARNING("Irreps have different nsym")
   ierr=ierr+1
 end if

 if (Irrep1%dim /= Irrep2%dim) then
   MSG_WARNING("Irreps have different dimensions")
   write(std_out,*)Irrep1%dim,Irrep2%dim
   ierr=ierr+1
 end if

 if (ii>Irrep2%dim .or. jj>Irrep2%dim .or. &
&    kk>Irrep1%dim .or. ll>Irrep1%dim) then
   MSG_WARNING("Wrong indeces")
   write(std_out,*)ii,Irrep2%dim,jj,Irrep2%dim,kk>Irrep1%dim,ll,Irrep1%dim
   ierr=ierr+1
 end if

 if (ierr/=0) RETURN

 do isym=1,nsym
   res = res + DCONJG(Irrep1%mat(ii,jj,isym)) * Irrep2%mat(kk,ll,isym)
 end do

end function sum_irreps
!!***

!----------------------------------------------------------------------

!!****f* m_ptgroups/groupk_free
!! NAME
!!  groupk_free
!!
!! FUNCTION
!!  Deallocate all memory allocate in the group_k_t.
!!
!! PARENTS
!!
!! CHILDREN
!!      irrep_free
!!
!! SOURCE

subroutine groupk_free(Gk)

 implicit none

!Arguments ------------------------------------
 type(group_k_t),intent(inout) :: Gk

! *************************************************************************

! integer
 if (allocated(Gk%class_ids))   then
   ABI_FREE(Gk%class_ids)
 end if
 if (allocated(Gk%sym))   then
   ABI_FREE(Gk%sym)
 end if

!real
 if (allocated(Gk%tnons)) then
   ABI_FREE(Gk%tnons)
 end if

!character
 if (allocated(Gk%class_names)) then
   ABI_FREE(Gk%class_names)
 end if

!type
 if (allocated(Gk%Irreps)) then
   call irrep_free(Gk%Irreps)
   ABI_DT_FREE(Gk%Irreps)
 end if

end subroutine groupk_free
!!***

END MODULE m_ptgroups
!!***
