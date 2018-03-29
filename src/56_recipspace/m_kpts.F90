!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_kpts
!! NAME
!!  m_kpts
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2008-2018 ABINIT group (MG, DRH, XG)
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

MODULE m_kpts

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_crystal

 use m_fstrings,       only : sjoin, itoa, ltoa
 use m_tetrahedron,    only : t_tetrahedron, init_tetra, destroy_tetra

 implicit none

 private

 public :: kpts_timrev_from_kptopt     ! Returns the value of timrev from kptopt
 public :: kpts_ibz_from_kptrlatt      ! Determines the IBZ, the weights and the BZ from kptrlatt
 public :: tetra_from_kptrlatt         ! Create an instance of `t_tetrahedron` from kptrlatt and shiftk
 public :: symkchk                     ! Checks that the set of k points has the full space group symmetry,
                                       ! modulo time reversal if appropriate.
!!***

!----------------------------------------------------------------------

contains  !============================================================
!!***

!!****f* m_kpts/kpts_timrev_from_kptopt
!! NAME
!!  kpts_timrev_from_kptopt
!!
!! FUNCTION
!!  Returns the value of timrev from kptopt
!!  1 if the use of time-reversal is allowed; 0 otherwise
!!
!! INPUTS
!!  kptopt=option for the generation of k points
!!    (defines whether spatical symmetries and/or time-reversal can be used)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer pure function kpts_timrev_from_kptopt(kptopt) result(timrev)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kpts_timrev_from_kptopt'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: kptopt

! *********************************************************************

 timrev = 1; if (any(kptopt == [3, 4])) timrev = 0

end function kpts_timrev_from_kptopt
!!***

!!****f* m_kpts/kpts_ibz_from_kptrlatt
!! NAME
!!  kpts_ibz_from_kptrlatt
!!
!! FUNCTION
!!  Determines the irreducible wedge, the corresponding weights and the list
!!  of k-points in the Brillouin Zone starting from kptrlatt and the set shifts.
!!
!! INPUTS
!!  cryst<crystal_t> = crystalline structure with info on symmetries and time-reversal.
!!  kptopt=option for the generation of k points
!!    (defines whether spatical symmetries and/or time-reversal can be used)
!!  kptrlatt(3,3)=integer coordinates of the primitive vectors of the
!!   lattice reciprocal to the k point lattice to be generated here
!!   If diagonal, the three values are the Monkhorst-Pack usual values, in case of simple cubic.
!!  nshiftk= number of shift vectors in the repeated cell
!!  shiftk(3,nshiftk) = vectors that will be used to determine the shifts from (0. 0. 0.).
!!
!! OUTPUT
!!  nkibz,nkbz = Number of points in IBZ and BZ, respectively.
!!  The following arrays are allocated and returned by the routine:
!!  kibz(3,nkibz) = k-points in the IBZ.
!!  wtk(nkibz) = weights of the k-points in the IBZ (normalized to one).
!!  kbz(3,nkbz) = k-points in the BZ.
!!  [new_kptrlatt] = New value of kptrlatt returned by getkgrid
!!  [new_shiftk(3,new_nshiftk)] = New set of shifts returned by getkgrid
!!
!! PARENTS
!!      m_dvdb,m_ebands,m_gruneisen,m_ifc,m_kpts,m_phgamma,m_phonons,m_sigmaph
!!
!! CHILDREN
!!      getkgrid,init_tetra,kpts_ibz_from_kptrlatt,listkk
!!
!! SOURCE

subroutine kpts_ibz_from_kptrlatt(cryst, kptrlatt, kptopt, nshiftk, shiftk, nkibz, kibz, wtk, nkbz, kbz, &
  new_kptrlatt, new_shiftk)  ! Optional


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kpts_ibz_from_kptrlatt'
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nshiftk,kptopt
 integer,intent(out) :: nkibz,nkbz
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 integer,optional,intent(out) :: new_kptrlatt(3,3)
 real(dp),intent(in) :: shiftk(3,nshiftk)
 real(dp),allocatable,intent(out) :: wtk(:),kibz(:,:),kbz(:,:)
 real(dp),optional,allocatable,intent(out) :: new_shiftk(:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: iout0=0,chksymbreak0=0,iscf2=2
 integer :: my_nshiftk,nkpt_computed
 real(dp) :: kptrlen
!arrays
 integer,parameter :: vacuum0(3)=[0,0,0]
 integer :: my_kptrlatt(3,3)
 real(dp) :: my_shiftk(3,210)

! *********************************************************************

 ! First call to getkgrid to obtain the number of points in the BZ.
 ABI_MALLOC(kibz, (3,0))
 ABI_MALLOC(wtk, (0))

 ! Copy kptrlatt and shifts because getkgrid can change them
 ! Be careful as getkgrid expects shiftk(3,210).
 ABI_CHECK(nshiftk > 0 .and. nshiftk <= 210, "nshiftk must be in [1,210]")
 my_nshiftk = nshiftk; my_shiftk = zero; my_shiftk(:,1:nshiftk) = shiftk
 my_kptrlatt = kptrlatt

 call getkgrid(chksymbreak0,iout0,iscf2,kibz,kptopt,my_kptrlatt,kptrlen,&
   cryst%nsym,0,nkibz,my_nshiftk,cryst%nsym,cryst%rprimd,my_shiftk,cryst%symafm,cryst%symrel,vacuum0,wtk)

 ABI_FREE(kibz)
 ABI_FREE(wtk)

 ! Recall getkgrid to get kibz and wtk.
 ABI_MALLOC(kibz, (3, nkibz))
 ABI_MALLOC(wtk, (nkibz))

 call getkgrid(chksymbreak0,iout0,iscf2,kibz,kptopt,my_kptrlatt,kptrlen,&
   cryst%nsym,nkibz,nkpt_computed,my_nshiftk,cryst%nsym,cryst%rprimd,my_shiftk,&
   cryst%symafm,cryst%symrel,vacuum0,wtk,fullbz=kbz)

 nkbz = size(kbz, dim=2)

 ! Optionally, return new shifts and new_kptrlatt
 if (present(new_shiftk)) then
   ABI_MALLOC(new_shiftk, (3, my_nshiftk))
   new_shiftk = my_shiftk(:, 1:my_nshiftk)
 end if
 if (present(new_kptrlatt)) new_kptrlatt = my_kptrlatt

 DBG_CHECK(abs(sum(wtk) - one) < tol10, "sum(wtk) != one")

end subroutine kpts_ibz_from_kptrlatt
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/tetra_from_kptrlatt
!! NAME
!! tetra_from_kptrlatt
!!
!! FUNCTION
!!  Create an instance of `t_tetrahedron` from kptrlatt and shiftk
!!
!! INPUTS
!!  cryst<cryst_t>=Crystalline structure.
!!  kptopt=Option for the k-point generation.
!!  kptrlatt(3,3)=k-point lattice specification
!!  nshiftk= number of shift vectors.
!!  shiftk(3,nshiftk)=shift vectors for k point generation
!!  nkibz=Number of points in the IBZ
!!  kibz(3,nkibz)=Reduced coordinates of the k-points in the IBZ.
!!
!! OUTPUT
!!  tetra<t_tetrahedron>=Tetrahedron object, fully initialized if ierr == 0.
!!  msg=Error message if ierr /= 0
!!  ierr=Exit status
!!
!! PARENTS
!!      gstate,wfk_analyze
!!
!! CHILDREN
!!      init_tetra,listkk,smpbz
!!
!! SOURCE

type(t_tetrahedron) function tetra_from_kptrlatt( &
&  cryst, kptopt, kptrlatt, nshiftk, shiftk, nkibz, kibz, msg, ierr) result (tetra)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'tetra_from_kptrlatt'
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: kptopt,nshiftk,nkibz
 integer,intent(out) :: ierr
 character(len=*),intent(out) :: msg
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 real(dp),intent(in) :: shiftk(3,nshiftk),kibz(3,nkibz)

!Local variables-------------------------------
!scalars
 integer :: nkfull,timrev,sppoldbl,my_nkibz,new_nshiftk
 real(dp) :: dksqmax
 character(len=80) :: errorstring
!arrays
 integer :: new_kptrlatt(3,3)
 integer,allocatable :: indkk(:,:)
 real(dp) :: rlatt(3,3),klatt(3,3)
 real(dp),allocatable :: kfull(:,:),my_kibz(:,:),my_wtk(:),new_shiftk(:,:)

! *************************************************************************

 ierr = 0

 ! Refuse only 1 kpoint: the algorithms are no longer valid. DOH!
 if (nkibz == 1) then
   msg = 'You need at least 2 kpoints to use the tetrahedron method.'
   ierr = 1; goto 10
 end if
 if (all(kptrlatt == 0)) then
   msg = 'Cannot generate tetrahedron because input kptrlatt == 0.'
   ierr = 1; goto 10
 end if
 if (kptopt <= 0) then
   msg = sjoin("Cannot generate tetrahedron because input kptopt:", itoa(kptopt))
   ierr = 1; goto 10
 end if

 call kpts_ibz_from_kptrlatt(cryst, kptrlatt, kptopt, nshiftk, shiftk, &
   my_nkibz, my_kibz, my_wtk, nkfull, kfull, new_kptrlatt=new_kptrlatt, new_shiftk=new_shiftk)

 ABI_FREE(my_kibz)
 ABI_FREE(my_wtk)
 new_nshiftk = size(new_shiftk, dim=2)

 if (my_nkibz /= nkibz) then
   msg = sjoin("Input nkibz:", itoa(nkibz), "does not agree with computed value:", itoa(my_nkibz))
   ierr = 1; goto 10
 end if

 ! Do not support new_nshiftk > 1: lattice must be decomposed into boxes
 ! and this is not always possible (I think) with bizzare shiftks
 ! normally at this point we have incorporated everything into
 ! new_kptrlatt, and only 1 shift is needed (in particular for MP grids).
 if (new_nshiftk > 1) then
   write(msg, "(9a)") &
     'Cannot create tetrahedron object...',ch10, &
     'Only simple lattices are supported. Action: use nshiftk=1.',ch10, &
     'new_shiftk: ', trim(ltoa(reshape(new_shiftk, [3*new_nshiftk]))),ch10, &
     'new_kptrlatt: ', trim(ltoa(reshape(new_kptrlatt, [9])))
   ierr = 2; goto 10
 end if

 ! Costruct full BZ and create mapping BZ --> IBZ
 ! Note:
 !   - we don't change the value of nsppol hence sppoldbl is set to 1
 !   - we use symrec (operations in reciprocal space)
 !
 sppoldbl = 1; timrev = kpts_timrev_from_kptopt(kptopt)
 ABI_MALLOC(indkk, (nkfull*sppoldbl,6))

 ! Compute k points from input file closest to the output file
 call listkk(dksqmax,cryst%gmet,indkk,kibz,kfull,nkibz,nkfull,cryst%nsym,&
    sppoldbl,cryst%symafm,cryst%symrec,timrev,use_symrec=.True.)

 if (dksqmax > tol12) then
   write(msg, '(3a,es16.6,6a)' )&
   'At least one of the k points could not be generated from a symmetrical one.',ch10,&
   'dksqmax=',dksqmax,ch10,&
   'new_kptrkatt= ',trim(ltoa(reshape(new_kptrlatt, [9]))),ch10,&
   'new_shiftk= ',trim(ltoa(reshape(new_shiftk, [3*new_nshiftk])))
   ierr = 2; goto 10
 end if

 rlatt = new_kptrlatt; call matr3inv(rlatt,klatt)

 call init_tetra(indkk(:,1), cryst%gprimd, klatt, kfull, nkfull, tetra, ierr, errorstring)
 if (ierr /= 0) msg = errorstring

 10 continue
 if (allocated(indkk)) then
   ABI_FREE(indkk)
 end if
 if (allocated(kfull)) then
   ABI_FREE(kfull)
 end if
 if (allocated(new_shiftk)) then
   ABI_FREE(new_shiftk)
 end if

end function tetra_from_kptrlatt
!!***

!!****f* m_kpts/symkchk
!! NAME
!! symkchk
!!
!! FUNCTION
!! Checks that the set of k points chosen for a response function
!! calculation has the full space group symmetry, modulo time reversal if appropriate.
!! Returns ierr/=0 with error message if not satisfied
!! Currently used only when strain perturbation is treated. Based on symkpt.
!!
!! INPUTS
!! kptns(3,nkpt)= k vectors in reciprocal space
!! nkpt = number of k-points whose weights are wtk
!! nsym=number of space group symmetries
!! symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!! timrev: if 1, the time reversal operation has to be taken into account
!! if 0, no time reversal symmetry.
!!
!! OUTPUT
!!  msg=Error message if ierr /= 0
!!
!! TODO
!!  This version should scale badly with the number of k-points. Replace loops with listkk
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

integer function symkchk(kptns,nkpt,nsym,symrec,timrev,errmsg) result(ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symkchk'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nkpt,nsym,timrev
 character(len=*),intent(out) :: errmsg
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 real(dp),intent(in) :: kptns(3,nkpt)

!Local variables -------------------------
!scalars
 integer :: identi,ii,ikpt,ikpt2,imatch,isym,jj,tident
 real(dp) :: difk,reduce
 character(len=500) :: message
!arrays
 real(dp) :: ksym(3)

! *********************************************************************
 ierr = 0

 if(timrev/=1 .and. timrev/=0)then
   write(errmsg, '(3a,i0,a)' )&
&   'timrev should be 0 or 1, while',ch10,&
&   'it is equal to ',timrev,'.'
   ierr = 1; return
 end if

 if(nsym/=1)then
!  Find the identity symmetry operation
   do isym=1,nsym
     tident=1
     do jj=1,3
       if(symrec(jj,jj,isym)/=1)tident=0
       do ii=1,3
         if( ii/=jj .and.&
&         symrec(ii,jj,isym)/=0)tident=0
       end do
     end do
     if(tident==1)then
       identi=isym
       call wrtout(std_out,sjoin(' symkchk: found identity with number:', itoa(identi)))
       exit
     end if
   end do
   if(tident==0)then
     errmsg = 'Did not found the identity operation.'
     ierr = 1; return
   end if
 end if

!Here begins the serious business
!The length sorting, etc. of symkpt have been dropped because the
!computational cost is estimated to be negligible.

 if(nsym>1 .or. timrev==1)then

!  Outer loop over kpts
   do ikpt=1,nkpt-1

!    Loop on the symmetries
!    For each k-point and each symmetry transformation, a matching
!    k-pointpt must be found, modulo time reversal if appropriate
     do isym=1,nsym

!      Get the symmetric of the vector
       do ii=1,3
         ksym(ii)= kptns(1,ikpt)*symrec(ii,1,isym)&
&         +kptns(2,ikpt)*symrec(ii,2,isym)&
&         +kptns(3,ikpt)*symrec(ii,3,isym)
       end do

!      Second loop k-points
       do ikpt2=1,nkpt

!        Test for match of symmetric and any vector (including original)
         imatch=1
         do ii=1,3
           difk= ksym(ii)-kptns(ii,ikpt2)
           reduce=difk-anint(difk)
           if(abs(reduce)>tol8)imatch=0
         end do
         if(imatch==1)exit

!        Test for match with time reversal
         if(timrev==1)then
           imatch=1
           do ii=1,3
             difk= ksym(ii)+kptns(ii,ikpt2)
             reduce=difk-anint(difk)
             if(abs(reduce)>tol8)imatch=0
           end do
           if(imatch==1)exit
         end if

       end do ! End secondary loop over k-points
       if (imatch/=1) then
         write(errmsg, '(a,a,a,i4,a,i4,a,a,a,a)' )&
&         'k-point set must have full space-group symmetry',ch10,&
&         'there is no match for kpt',ikpt,' transformed by symmetry',isym,ch10,&
&         'Action: change kptopt to 2 or 3 and/or change or use shiftk',ch10,&
&         'shiftk = 0 0 0 is always a safe choice.'
         ierr = 2; return
       end if

     end do ! End loop on isym
   end do ! End primary loop over k-points

   write(message,'(a)')' symkchk : k-point set has full space-group symmetry.'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
 end if

end function symkchk
!!***

end module m_kpts
!!***
