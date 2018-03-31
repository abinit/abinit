!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_kpts
!! NAME
!!  m_kpts
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2008-2018 ABINIT group (XG, MG, DRH)
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
 use m_sort

 use m_fstrings,       only : sjoin, itoa, ltoa
 use m_tetrahedron,    only : t_tetrahedron, init_tetra, destroy_tetra

 implicit none

 private

 public :: kpts_timrev_from_kptopt     ! Returns the value of timrev from kptopt
 public :: kpts_ibz_from_kptrlatt      ! Determines the IBZ, the weights and the BZ from kptrlatt
 public :: tetra_from_kptrlatt         ! Create an instance of `t_tetrahedron` from kptrlatt and shiftk
 public :: symkchk                     ! Checks that the set of k points has the full space group symmetry,
                                       ! modulo time reversal if appropriate.
 public :: listkk                      ! Find correspondence between two set of k-points.
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

!!****f* m_kpts/listkk
!! NAME
!! listkk
!!
!! FUNCTION
!! Given a list of nkpt1 initial k points kptns1 and a list of nkpt2
!! final k points kptns2, associates each final k pt with a "closest"
!! initial k point (or symmetric thereof, also taking possible umklapp)
!! as determined by a metric gmet, that commutes with the symmetry operations.
!! The algorithm does not scale as nkpt1 times nkpt2, thanks
!! to the ordering of the kptns1 and kptns2 vectors according to their
!! lengths, and comparison first between vectors of similar lengths.
!! Returns indirect indexing list indkk.
!!
!! INPUTS
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!  kptns1(3,nkpt1)=list of initial k points (reduced coordinates)
!!  kptns2(3,nkpt2)=list of final k points
!!  nkpt1=number of initial k points
!!  nkpt2=number of final k points
!!  nsym=number of symmetry elements in space group
!!  sppoldbl=if 1, no spin-polarisation doubling
!!           if 2, spin-polarisation doubling using symafm
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symmat(3,3,nsym)=symmetry operations (symrel or symrec, depending on
!!                   value of use_symrec
!!  timrev=1 if the use of time-reversal is allowed; 0 otherwise
!!  use_symrec: if present and true, symmat assumed to be symrec, otherwise assumed to be symrel (default)
!!
!! OUTPUT
!!  dksqmax=maximal value of the norm**2 of the difference between
!!    a kpt2 vector and the closest k-point found from the kptns1 set, using symmetries.
!!  indkk(nkpt2*sppoldbl,6)=describe k point number of kpt1 that allows to
!!    generate wavefunctions closest to given kpt2
!!    if sppoldbl=2, use symafm to generate spin down wfs from spin up wfs
!!
!!    indkk(:,1)=k point number of kptns1
!!    indkk(:,2)=symmetry operation to be applied to kpt1, to give kpt1a
!!      (if 0, means no symmetry operation, equivalent to identity )
!!    indkk(:,3:5)=shift in reciprocal space to be given to kpt1a,
!!      to give kpt1b, that is the closest to kpt2.
!!    indkk(:,6)=1 if time-reversal was used to generate kpt1a from kpt1, 0 otherwise
!!
!! NOTES
!!  The tolerances tol12 and tol8 aims at giving a machine-independent ordering.
!!  (this trick is used in bonds.f, listkk.f, prtrhomxmn.f and rsiaf9.f)
!!  The tolerance tol12 is used for each component of the k vectors,
!!  and for the length of the vectors
!!  while the tolerance tol8 is used for the comparison of the squared lengths
!!  of the separate vectors.
!!
!! PARENTS
!!      initberry,initorbmag,inwffil,m_dvdb,m_ebands,m_eprenorms,m_exc_diago
!!      m_fock,m_fstab,m_haydock,m_ifc,m_kpts,m_phgamma,m_sigmaph,mlwfovlp_qp
!!
!! CHILDREN
!!      sort_dp,timab
!!
!! SOURCE

subroutine listkk(dksqmax,gmet,indkk,kptns1,kptns2,nkpt1,nkpt2,nsym,&
& sppoldbl,symafm,symmat,timrev,use_symrec)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'listkk'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt1,nkpt2,nsym,sppoldbl,timrev
 real(dp),intent(out) :: dksqmax
 logical,optional,intent(in) :: use_symrec
!arrays
 integer,intent(in) :: symafm(nsym),symmat(3,3,nsym)
 integer,intent(out) :: indkk(nkpt2*sppoldbl,6)
 real(dp),intent(in) :: gmet(3,3),kptns1(3,nkpt1),kptns2(3,nkpt2)

!Local variables-------------------------------
!scalars
 integer :: l3,ig1,ig2,ig3,ii,ikpg1,ikpt1,ikpt2,ikpt2_done
 integer :: ilarger,ismaller,itrial
 integer :: isppol,isym,itimrev,jkpt1,jsym,jtime,limit
 integer :: nsym_used,timrev_used,usesym
 real(dp) :: dksq,dksqmn,lk2,llarger,ldiff,lsmaller,ltrial,min_l
 character(len=500) :: message
!arrays
 integer :: dkint(3),jdkint(3),k1int(3),k2int(3)
 integer, allocatable :: isort(:)
 real(dp) :: tsec(2)
 real(dp) :: dk(3),kpg1(3),kpt1a(3),k1(3),k2(3)
!real(dp) :: kasq,ka(3)
 real(dp),allocatable :: lkpg1(:),lkpg1_sorted(:)

! *************************************************************************

!write(std_out,*)' listkk : nkpt1,nkpt2,nsym=',nkpt1,nkpt2,nsym
 call timab(1021,1,tsec)

 if(sppoldbl<1 .or. sppoldbl>2)then
   write(message, '(a,i4,3a)' )&
&   'The value of sppoldbl is',sppoldbl,',',ch10,&
&   'but it should be either 1 or 2.'
   MSG_BUG(message)
 end if

!When usesym=0, the old way of converting the wavefunctions (without
!using the symmetries), is recovered.
 usesym=1

 nsym_used=nsym
 timrev_used=timrev
 if(usesym==0)nsym_used=1
 if(usesym==0)timrev_used=0

!Precompute the length of the kpt1 vectors, also taking into account
!possible umpklapp vectors
 limit=1 ; l3 = (2*limit+1)**3
 ABI_ALLOCATE(lkpg1,(l3*nkpt1))
 ABI_ALLOCATE(lkpg1_sorted,(l3*nkpt1))
 ABI_ALLOCATE(isort,(l3*nkpt1))
!write(std_out,*)' List of kpt1 vectors '
!write(std_out,*)' Length of the kpt1 vectors :'

 do ikpt1=1,nkpt1
   k1(:)=kptns1(:,ikpt1)
!  write(std_out,*)ikpt1,k1(:)
   k1int(:)=nint(k1(:)+tol12)
   k1(:)=k1(:)-k1int(:)
   do ig1=-limit,limit
     kpg1(1)=k1(1)+ig1
     do ig2=-limit,limit
       kpg1(2)=k1(2)+ig2
       do ig3=-limit,limit
         kpg1(3)=k1(3)+ig3

         ikpg1=ig1+limit+1 + (2*limit+1)*(ig2+limit) + (2*limit+1)**2*(ig3+limit) + l3*(ikpt1-1)
!        Compute the norm of the vector (also taking into account possible umklapp)
         lkpg1(ikpg1)=sqrt(gmet(1,1)*kpg1(1)**2+gmet(2,2)*kpg1(2)**2+&
&         gmet(3,3)*kpg1(3)**2+two*(gmet(2,1)*kpg1(2)*kpg1(1)+&
&         gmet(3,2)*kpg1(3)*kpg1(2)+gmet(3,1)*kpg1(3)*kpg1(1)))
         lkpg1_sorted(ikpg1)=lkpg1(ikpg1)
         isort(ikpg1)=ikpg1
!        write(std_out,*)' ikpt1,ig1,ig2,ig3,lkpg1=',ikpt1,ig1,ig2,ig3,lkpg1(ikpg1)
       end do
     end do
   end do
 end do

 call sort_dp( l3*nkpt1,lkpg1_sorted,isort,tol12)

!DEBUG
!write(std_out,*)' listkk : output list of kpt1 for checking purposes '
!write(std_out,*)' ii,ikpt1,isort(ii)-l3*(ikpt1-1),lkpg1_sorted(ii),lkpg1(isort(ii)) '
!do ii=1,l3*nkpt1
!ikpt1=(isort(ii)-1)/l3+1
!write(std_out,*)ii,ikpt1,isort(ii)-l3*(ikpt1-1),lkpg1_sorted(ii),lkpg1(isort(ii))
!enddo
!stop
!ENDDEBUG

 dksqmax=zero
 do isppol=1,sppoldbl
   do ikpt2=1,nkpt2

     ikpt2_done=0
!    Precompute the length of the kpt2 vector, with the Umklapp vector such that it is the closest to the Gamma point
     k2(:)=kptns2(:,ikpt2)
     k2int(:)=nint(k2(:)+tol12)
     k2(:)=k2(:)-k2int(:)
     lk2=sqrt(gmet(1,1)*k2(1)**2+gmet(2,2)*k2(2)**2+&
&     gmet(3,3)*k2(3)**2+two*(gmet(2,1)*k2(2)*k2(1)+&
&     gmet(3,2)*k2(3)*k2(2)+gmet(3,1)*k2(3)*k2(1)))

!    DEBUG
!    write(std_out, '(a,i4,7es16.6)' )' listkk : ikpt2,kptns2(:,ikpt2),k2(:),lk2=',ikpt2,kptns2(:,ikpt2),k2(:),lk2
!    if(ikpt2/=17)cycle
!    ENDDEBUG

!    Find the kpt1 vector whose length is the most similar to the length of lk2
!    up to a tolerance. Use a bissection algorithm.
     ismaller=0       ; lsmaller=zero
     ilarger=l3*nkpt1+1 ; llarger=huge(one)
!    This loop should never reach l3*nkpt1, since this is a bissection algorithm
     do ii=1,l3*nkpt1
       if((ilarger-ismaller)<2 .or. (llarger-lsmaller)<2*tol12)exit
       itrial=(ilarger+ismaller)/2 ; ltrial=lkpg1_sorted(itrial)
       if((ltrial-lk2)>tol12)then
         ilarger=itrial ; llarger=ltrial
       else if((ltrial-lk2)<-tol12)then
         ismaller=itrial ; lsmaller=ltrial
       else
         ismaller=itrial ; lsmaller=ltrial
         ilarger=itrial ; llarger=ltrial
       end if
     end do
     itrial=ismaller
     if(abs(llarger-lk2)<abs(lsmaller-lk2)-tol12)itrial=ilarger
     if(itrial==0)itrial=ilarger
     ismaller=itrial ; ilarger=itrial
!    write(std_out,*)' listkk : starting search at itrial=',itrial

     dksqmn=huge(one)

!    The ii index is dummy. This avoids an infinite loop.
     do ii=1,l3*nkpt1
!      do ikpt1=1,nkpt1

!      If the difference in length between the trial vector and the target vector is bigger
!      than the already achieved distance, the search is finished ...
       ldiff=abs(lkpg1_sorted(itrial)-lk2)


!      DEBUG
!      write(std_out,*)' listkk : ii,itrial,lkpg1_sorted(itrial),lk2,ldiff,dksqmn=',ii,itrial,lkpg1_sorted(itrial),lk2,ldiff,dksqmn
!      ENDDEBUG
       if(ldiff**2>dksqmn+tol8)exit

!      If this k-point has already been examined in a previous batch, skip it
!      First, compute the minimum of the difference of length of the sets of associated vectors thanks to Umklapp vectors
!      with the target vector
       ikpt1=(isort(itrial)-1)/l3+1
       min_l=minval(abs(lkpg1((ikpt1-1)*l3+1:(ikpt1-1)*l3+l3)-lk2))
!      Then compare with the current ldiff

!      DEBUG
!      write(std_out,*)' listkk : ikpt1,min_l,ldiff=',ikpt1,min_l,ldiff
!      ENDDEBUG

       if(min_l > ldiff-tol12)then

!        Now, will examine the trial vector, and the symmetric ones
!MG FIXME:
! Here there's a possible problem with the order of symmetries because
! in symkpt, time-reversal is the innermost loop. This can create inconsistencies in the symmetry tables.
         do itimrev=0,timrev_used
           do isym=1,nsym_used

!            Select magnetic characteristic of symmetries
             if(isppol==1 .and. symafm(isym)==-1)cycle
             if(isppol==2 .and. symafm(isym)==1)cycle

!            Compute symmetric point to kpt1
             if(usesym==1)then
!              original code only used transpose(symrel)
!              kpt1a(:)=symrel(1,:,isym)*kptns1(1,ikpt1)+&
!              &             symrel(2,:,isym)*kptns1(2,ikpt1)+&
!              &             symrel(3,:,isym)*kptns1(3,ikpt1)
               if (present(use_symrec)) then
                 if (use_symrec) then
                   kpt1a(:) = MATMUL(symmat(:,:,isym),kptns1(:,ikpt1))
                 else
                   kpt1a(:) = MATMUL(TRANSPOSE(symmat(:,:,isym)),kptns1(:,ikpt1))
                 end if
               else
                 kpt1a(:) = MATMUL(TRANSPOSE(symmat(:,:,isym)),kptns1(:,ikpt1))
               end if
               kpt1a(:)=(1-2*itimrev)*kpt1a(:)
             else
               kpt1a(:)=kptns1(:,ikpt1)
             end if

!            Compute difference with respect to kpt2, modulo a lattice vector
             dk(:)=kptns2(:,ikpt2)-kpt1a(:)
             if(usesym==1)then
!              The tolerance insure similar behaviour on different platforms
!              XG120418 : Actually, *assumes* that the closest point will have reduced
!              coordinates differing by less than 1/2 . There might be elongated
!              cells where this is not correct ...
               dkint(:)=nint(dk(:)+tol12)
               dk(:)=dk(:)-dkint(:)
             else
               dkint(:)=0
             end if

!            Compute norm of the difference vector, and update kpt1 if better.
             dksq=gmet(1,1)*dk(1)**2+gmet(2,2)*dk(2)**2+&
&             gmet(3,3)*dk(3)**2+two*(gmet(2,1)*dk(2)*dk(1)+&
&             gmet(3,2)*dk(3)*dk(2)+gmet(3,1)*dk(3)*dk(1))

             if (dksq<dksqmn+tol8) then

!              If exactly the right point (without using symmetries neither umklapp vector), will exit the search
!              Note that in this condition, each coordinate is tested separately, without squaring. So, it is a much stronger
!              condition than dksqmn<tol12
               if(sum(abs(kptns2(:,ikpt2)-kptns1(:,ikpt1)))<3*tol12)then
                 ikpt2_done=1
               end if

!              Update in three cases : either if succeeded to have exactly the vector, or the distance is better,
!              or the distance is only slightly worsened so select the lowest itimrev, isym or ikpt1, in order to respect previous ordering
               if(  ikpt2_done==1 .or. &
&               dksq+tol12<dksqmn .or. &
&               ( abs(dksq-dksqmn)<tol12 .and. &
&               ((itimrev<jtime) .or. &
&               (itimrev==jtime .and. isym<jsym) .or. &
&               (itimrev==jtime .and. isym==jsym .and. ikpt1<jkpt1))))then

                 dksqmn=dksq
                 jkpt1=ikpt1
                 jsym=isym
                 jtime=itimrev
                 jdkint(:)=dkint(:)

!                DEBUG
!                write(std_out,*)' ikpt1,ikpt2=',ikpt1,ikpt2
!                write(std_out,*)' timrev_used=',timrev_used
!                write(std_out,*)' Succeeded to lower dskmn,ikpt2_done=',dksqmn,ikpt2_done
!                write(std_out,*)' ikpt1,isym,dkint(:),itimrev=',ikpt1,isym,dkint(:),itimrev
!                ka(:)=kpt1a(:)+dkint(:)
!                write(std_out,*)'        k1=',kpt1a(:)
!                write(std_out,*)'     dkint=',dkint(:)
!                write(std_out,*)' Actual k1=',ka(:)
!                write(std_out,*)'        k2=',kptns2(:,ikpt2)
!                kasq=gmet(1,1)*ka(1)**2+gmet(2,2)*ka(2)**2+&
!                &                  gmet(3,3)*ka(3)**2+two*(gmet(2,1)*ka(2)*ka(1)+&
!                &                  gmet(3,2)*ka(3)*ka(2)+gmet(3,1)*ka(3)*ka(1))
!                write(std_out,*)' Actual k1sq=',kasq
!                ENDDEBUG
               end if

             end if
             if(ikpt2_done==1)exit

           end do ! isym
           if(ikpt2_done==1)exit

         end do ! itimrev
         if(ikpt2_done==1)exit

       end if

!      Update the interval that has been explored
       if(itrial<ismaller)ismaller=itrial
       if(itrial>ilarger)ilarger=itrial

!      Select the next index to be tried (preferably the smaller indices, but this is a bit arbitrary).

!      DEBUG
!      write(std_out,*)' before choosing the next index :'
!      write(std_out,*)' ismaller,itrial,ilarger=',ismaller,itrial,ilarger
!      write(std_out,*)' lkpg1_sorted(ismaller-1),lk2,lkpg1_sorted(ilarger+1)=',lkpg1_sorted(ismaller-1),lk2,lkpg1_sorted(ilarger+1)
!      ENDDEBUG
       if(ismaller>1 .and. ilarger<l3*nkpt1)then
         if(abs(lkpg1_sorted(ismaller-1)-lk2)<abs(lkpg1_sorted(ilarger+1)-lk2)+tol12)then
           itrial=ismaller-1
         else
           itrial=ilarger+1
         end if
       end if
       if(ismaller==1 .and. ilarger<l3*nkpt1)itrial=ilarger+1
       if(ismaller>1 .and. ilarger==l3*nkpt1)itrial=ismaller-1
!      if(ismaller==1 .and. ilarger==l3*nkpt1), we are done with the loop !

     end do ! ikpt1

     indkk(ikpt2+(isppol-1)*nkpt2,1)=jkpt1
     indkk(ikpt2+(isppol-1)*nkpt2,2)=jsym
     indkk(ikpt2+(isppol-1)*nkpt2,3:5)=jdkint(:)
     indkk(ikpt2+(isppol-1)*nkpt2,6)=jtime
     dksqmax=max(dksqmax,dksqmn)

     if(dksqmn<-tol12)then
       write(message, '(a,es16.6)' )'  The minimum square of dk has negative norm: dksqmn=',dksqmn
       MSG_BUG(message)
     end if

!    DEBUG
!    write(std_out,'(a,i6,i2,2x,i6,5i3,es24.14)' )' listkk: ikpt2,isppol,indkk(ikpt2+(isppol-1)*nkpt2,:)=',ikpt2,isppol,indkk(ikpt2+(isppol-1)*nkpt2,:),dksqmn
!    if(nkpt1==17)stop
!    ENDDEBUG

   end do ! ikpt2
 end do ! isppol

 ABI_DEALLOCATE(isort)
 ABI_DEALLOCATE(lkpg1)
 ABI_DEALLOCATE(lkpg1_sorted)

 call timab(1021,2,tsec)

end subroutine listkk
!!***

end module m_kpts
!!***
