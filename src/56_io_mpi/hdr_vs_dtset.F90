!{\src2tex{textfont=tt}}
!!****f* ABINIT/hdr_vs_dtset
!! NAME
!! hdr_vs_dtset
!!
!! FUNCTION
!!  Check the compatibility of the Abinit header with respect to the
!!  input variables defined in the input file. 
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  Dtset<type(dataset_type)>=all input variables for this dataset
!!  Hdr <type(hdr_type)>=the header structured variable
!!
!! OUTPUT
!!  Only check
!!
!! PARENTS
!!      eph,setup_bse,setup_screening,setup_sigma,wfk_analyze
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine hdr_vs_dtset(Hdr,Dtset)
    
 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi

 use m_fstrings, only : ltoa
 use m_crystal,  only : print_symmetries

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_vs_dtset'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(Hdr_type),intent(in) :: Hdr
 type(Dataset_type),intent(in) :: Dtset

!Local variables-------------------------------
 integer :: ik,jj,ierr
 logical :: test
 logical :: tsymrel,ttnons,tsymafm
 character(len=500) :: msg      
! *************************************************************************

!=== Check basic dimensions ===
 ierr=0
 call compare_int('natom',  Hdr%natom,  Dtset%natom,  ierr)
 call compare_int('nkpt',   Hdr%nkpt,   Dtset%nkpt,   ierr)
 call compare_int('npsp',   Hdr%npsp,   Dtset%npsp,   ierr)
 call compare_int('nspden', Hdr%nspden, Dtset%nspden, ierr)
 call compare_int('nspinor',Hdr%nspinor,Dtset%nspinor,ierr)
 call compare_int('nsppol', Hdr%nsppol, Dtset%nsppol, ierr)
 call compare_int('nsym',   Hdr%nsym,   Dtset%nsym,   ierr)
 call compare_int('ntypat', Hdr%ntypat, Dtset%ntypat, ierr)
 call compare_int('usepaw', Hdr%usepaw, Dtset%usepaw, ierr)
 call compare_int('usewvl', Hdr%usewvl, Dtset%usewvl, ierr)
 call compare_int('kptopt', Hdr%kptopt, Dtset%kptopt, ierr)
 call compare_int('pawcpxocc', Hdr%pawcpxocc, Dtset%pawcpxocc, ierr)
 call compare_int('nshiftk_orig', Hdr%nshiftk_orig, Dtset%nshiftk_orig, ierr)
 call compare_int('nshiftk', Hdr%nshiftk, Dtset%nshiftk, ierr)

!=== The number of fatal errors must be zero ===
 if (ierr/=0) then 
   write(msg,'(3a)')&
&   'Cannot continue, basic dimensions reported in the header do not agree with input file. ',ch10,&
&   'Check consistency between the content of the external file and the input file. '
   MSG_ERROR(msg)
 end if

 test=ALL(ABS(Hdr%xred-Dtset%xred_orig(:,1:Dtset%natom,1))<tol6)
 ABI_CHECK(test,'Mismatch in xred')

 test=ALL(Hdr%typat==Dtset%typat(1:Dtset%natom)) 
 ABI_CHECK(test,'Mismatch in typat')
!
!* Check if the lattice from the input file agrees with that read from the KSS file
 if ( (ANY(ABS(Hdr%rprimd-Dtset%rprimd_orig(1:3,1:3,1))>tol6)) ) then
   write(msg,'(6a)')ch10,&
&   ' hdr_vs_dtset : ERROR - ',ch10,&
&   ' real lattice vectors read from Header ',ch10,&
&   ' differ from the values specified in the input file'
   call wrtout(std_out,msg,'COLL')
   write(msg,'(3a,3(3es16.6),3a,3(3es16.6),3a)')ch10,&
&   ' rprimd from Hdr file   = ',ch10,(Hdr%rprimd(:,jj),jj=1,3),ch10,&
&   ' rprimd from input file = ',ch10,(Dtset%rprimd_orig(:,jj,1),jj=1,3),ch10,ch10,&
&   '  Modify the lattice vectors in the input file '
   call wrtout(std_out,msg,'COLL') 
   MSG_ERROR("")
 end if 

!=== Check symmetry operations ===
 tsymrel=(ALL(Hdr%symrel==Dtset%symrel(:,:,1:Dtset%nsym)))
 if (.not.tsymrel) then
   write(msg,'(6a)')ch10,&
&   ' hdr_vs_dtset : ERROR - ',ch10,&
&   ' real space symmetries read from Header ',ch10,&
&   ' differ from the values inferred from the input file'
   call wrtout(std_out,msg,'COLL')
   tsymrel=.FALSE.
 end if 

 ttnons=ALL(ABS(Hdr%tnons-Dtset%tnons(:,1:Dtset%nsym))<tol6)
 if (.not.ttnons) then
   write(msg,'(6a)')ch10,&
&   ' hdr_vs_dtset : ERROR - ',ch10,&
&   ' fractional translations read from Header ',ch10,&
&   ' differ from the values inferred from the input file'
   call wrtout(std_out,msg,'COLL')
   ttnons=.FALSE.
 end if 

 tsymafm=ALL(Hdr%symafm==Dtset%symafm(1:Dtset%nsym))
 if (.not.tsymafm) then
   write(msg,'(6a)')ch10,&
&   ' hdr_vs_dtset : ERROR - ',ch10,&
&   ' AFM symmetries read from Header ',ch10,&
&   ' differ from the values inferred from the input file'
   call wrtout(std_out,msg,'COLL')
   tsymafm=.FALSE.
 end if

 if (.not.(tsymrel.and.ttnons.and.tsymafm)) then
   write(msg,'(a)')' Header ' 
   call wrtout(std_out,msg,'COLL') 
   call print_symmetries(Hdr%nsym,Hdr%symrel,Hdr%tnons,Hdr%symafm)
   write(msg,'(a)')' Dtset  ' 
   call wrtout(std_out,msg,'COLL') 
   call print_symmetries(Dtset%nsym,Dtset%symrel,Dtset%tnons,Dtset%symafm)
   MSG_ERROR('Check symmetry operations')
 end if

 if (abs(Dtset%nelect-hdr%nelect)>tol6) then
   write(msg,'(2(a,f8.2))')&
&   "File contains ", hdr%nelect," electrons but nelect initialized from input is ",Dtset%nelect
   MSG_ERROR(msg)
 end if
 if (abs(Dtset%charge-hdr%charge)>tol6) then
   write(msg,'(2(a,f8.2))')&
&   "File contains charge ", hdr%charge," but charge from input is ",Dtset%charge
   MSG_ERROR(msg)
 end if

 if (any(hdr%kptrlatt_orig /= dtset%kptrlatt_orig)) then
   write(msg,"(5a)")&
   "hdr%kptrlatt_orig: ",trim(ltoa(reshape(hdr%kptrlatt_orig,[9]))),ch10,&
   "dtset%kptrlatt_orig: ",trim(ltoa(reshape(dtset%kptrlatt_orig, [9])))
   MSG_ERROR(msg)
 end if

 if (any(hdr%kptrlatt /= dtset%kptrlatt)) then
   write(msg,"(5a)")&
   "hdr%kptrlatt: ",trim(ltoa(reshape(hdr%kptrlatt, [9]))),ch10,&
   "dtset%kptrlatt: ",trim(ltoa(reshape(dtset%kptrlatt, [9])))
   MSG_ERROR(msg)
 end if

 if (any(abs(hdr%shiftk_orig - dtset%shiftk_orig(:,1:dtset%nshiftk_orig)) > tol6)) then
   write(msg,"(5a)")&
   "hdr%shiftk_orig: ",trim(ltoa(reshape(hdr%shiftk_orig, [3*hdr%nshiftk_orig]))),ch10,&
   "dtset%shiftk_orig: ",trim(ltoa(reshape(dtset%shiftk_orig, [3*dtset%nshiftk_orig])))
   MSG_ERROR(msg)
 end if

 if (any(abs(hdr%shiftk - dtset%shiftk(:,1:dtset%nshiftk)) > tol6)) then
   write(msg,"(5a)")&
   "hdr%shiftk: ",trim(ltoa(reshape(hdr%shiftk, [3*hdr%nshiftk]))),ch10,&
   "dtset%shiftk: ",trim(ltoa(reshape(dtset%shiftk, [3*dtset%nshiftk])))
   MSG_ERROR(msg)
 end if

!* Check if the k-points from the input file agrees with that read from the WFK file
 if ( (ANY(ABS(Hdr%kptns(:,:)-Dtset%kpt(:,1:Dtset%nkpt))>tol6)) ) then
   write(msg,'(9a)')ch10,&
&   ' hdr_vs_dtset : ERROR - ',ch10,&
&   '  k-points read from Header ',ch10,&
&   '  differ from the values specified in the input file',ch10,&
&   '  k-points from Hdr file                        | k-points from input file ',ch10
   call wrtout(std_out,msg,'COLL') 
   do ik=1,Dtset%nkpt
     write(msg,'(3(3es16.6,3x))')Hdr%kptns(:,ik),Dtset%kpt(:,ik)
     call wrtout(std_out,msg,'COLL') 
   end do
   MSG_ERROR('Modify the k-mesh in the input file')
 end if 

 if (ANY(ABS(Hdr%wtk(:)-Dtset%wtk(1:Dtset%nkpt))>tol6)) then
   write(msg,'(9a)')ch10,&
&   ' hdr_vs_dtset : ERROR - ',ch10,&
&   '  k-point weights read from Header ',ch10,&
&   '  differ from the values specified in the input file',ch10,&
&   '  Hdr file  |  File ',ch10
   call wrtout(std_out,msg,'COLL') 
   do ik=1,Dtset%nkpt
     write(msg,'(2(f11.5,1x))')Hdr%wtk(ik),Dtset%wtk(ik)
     call wrtout(std_out,msg,'COLL') 
   end do
   MSG_ERROR('Check the k-mesh and the symmetries of the system. ')
 end if 

!Check istwfk storage
 if ( (ANY(Hdr%istwfk(:)/=Dtset%istwfk(1:Dtset%nkpt))) ) then
   write(msg,'(9a)')ch10,&
&   ' hdr_vs_dtset : ERROR - ',ch10,&
&   '  istwfk read from Header ',ch10,&
&   '  differ from the values specified in the input file',ch10,&
&   '  Hdr | input ',ch10
   call wrtout(std_out,msg,'COLL') 
   do ik=1,Dtset%nkpt
     write(msg,'(i5,3x,i5)')Hdr%istwfk(ik),Dtset%istwfk(ik)
     call wrtout(std_out,msg,'COLL') 
   end do
   MSG_ERROR('Modify istwfk in the input file')
 end if 

 CONTAINS  !===========================================================
!!***

!!****f* hdr_vs_dtset/compare_int
!! NAME
!! compare_int
!!
!! FUNCTION
!!  Compare two int value and may raise an exception on error.
!!
!! INPUTS
!!  name=Name of the variable
!!  iexp= expected value.
!!  ifound=the actuval value
!!
!! SIDE EFFECTS
!!  ierr=increased by one if values differ 
!!
!! PARENTS
!!      hdr_vs_dtset
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

 subroutine compare_int(name,iexp,ifound,ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'compare_int'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: iexp,ifound
 integer,intent(inout) :: ierr
 character(len=*),intent(in) :: name

!Local variables-------------------------------
 logical :: leq                       
 character(len=500) :: msg                                              
! *************************************************************************

   leq=(iexp==ifound)

   if (.not.leq) then 
     write(msg,'(4a,i6,a,i6)')ch10,&
     ' hdr_vs_dtset : WARNING - Mismatch in '//TRIM(name),ch10,&
     '  Expected = ',iexp,' Found = ',ifound
     call wrtout(std_out,msg,'COLL') 
!      Increase ierr to signal we should stop in the caller.
     ierr=ierr+1 
   end if

 end subroutine compare_int
!!***

end subroutine hdr_vs_dtset
