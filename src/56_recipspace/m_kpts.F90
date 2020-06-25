!!****m* ABINIT/m_kpts
!! NAME
!!  m_kpts
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (XG, MG, MJV, DRH, DCA, JCC, MM)
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

module m_kpts

 use defs_basis
 use m_errors
 use m_abicore
 use m_crystal
 use m_sort
 use m_krank
 use m_htetra
 use m_xmpi

 use m_time,           only : timab
 use m_time,           only : cwtime, cwtime_report
 use m_copy,           only : alloc_copy
 use m_symtk,          only : mati3inv, mati3det, matr3inv, smallprim
 use m_fstrings,       only : sjoin, itoa, ltoa
 use m_numeric_tools,  only : wrap2_pmhalf
 use m_geometry,       only : metric
 use m_symkpt,         only : symkpt, symkpt_new

 implicit none

 private

 public :: kpts_timrev_from_kptopt   ! Returns the value of timrev from kptopt
 public :: kpts_ibz_from_kptrlatt    ! Determines the IBZ, the weights and the BZ from kptrlatt
 public :: tetra_from_kptrlatt       ! Create an instance from kptrlatt and shiftk
 public :: symkchk                   ! Checks that the set of k points has the full space group symmetry,
                                     ! modulo time reversal if appropriate.
 public :: listkk                    ! Find correspondence between two set of k-points.
 public :: getkgrid                  ! Compute the grid of k points in the irreducible Brillouin zone.
 !FIXME: Deprecated
 public :: get_full_kgrid            ! Create full grid of kpoints and find equivalent irred ones.
                                     ! Duplicates work in getkgrid, but need all outputs of kpt_fullbz, and indkpt

 private :: get_kpt_fullbz           ! Create full grid of kpoints from kptrlatt and shiftk

 public :: smpbz                     ! Generate a set of special k (or q) points which samples in a homogeneous way the BZ
 public :: testkgrid                 ! Test different grids of k points.

 ! FIXME: deprecated
 public :: mknormpath
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
!!    (defines whether spatial symmetries and/or time-reversal can be used)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer pure function kpts_timrev_from_kptopt(kptopt) result(timrev)

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
!!    (defines whether spatial symmetries and/or time-reversal can be used)
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
!!  [bz2ibz(6,nkbz)]=Mapping BZ --> IBZ
!!
!! PARENTS
!!      m_dvdb,m_ebands,m_gruneisen,m_ifc,m_kpts,m_phgamma,m_phonons,m_sigmaph
!!
!! CHILDREN
!!      getkgrid,kpts_ibz_from_kptrlatt,listkk
!!
!! SOURCE

subroutine kpts_ibz_from_kptrlatt(cryst, kptrlatt, kptopt, nshiftk, shiftk, nkibz, kibz, wtk, nkbz, kbz, &
  new_kptrlatt, new_shiftk, bz2ibz)  ! Optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nshiftk,kptopt
 integer,intent(out) :: nkibz,nkbz
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 integer,optional,allocatable,intent(out) :: bz2ibz(:,:)
 integer,optional,intent(out) :: new_kptrlatt(3,3)
 real(dp),intent(in) :: shiftk(3,nshiftk)
 real(dp),allocatable,intent(out) :: wtk(:),kibz(:,:),kbz(:,:)
 real(dp),optional,allocatable,intent(out) :: new_shiftk(:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: iout0 = 0, chksymbreak0 = 0, iscf2 = 2
 integer :: my_nshiftk
 real(dp) :: kptrlen
!arrays
 integer,parameter :: vacuum0(3) = [0, 0, 0]
 integer :: my_kptrlatt(3,3)
 integer,allocatable :: indkpt(:),bz2ibz_smap(:,:)
 real(dp) :: my_shiftk(3,MAX_NSHIFTK)

! *********************************************************************

 ! Copy kptrlatt and shifts because getkgrid can change them
 ! Be careful as getkgrid expects shiftk(3,MAX_NSHIFTK).
 ABI_CHECK(nshiftk > 0 .and. nshiftk <= MAX_NSHIFTK, sjoin("nshiftk must be between 1 and", itoa(MAX_NSHIFTK)))
 my_nshiftk = nshiftk; my_shiftk = zero; my_shiftk(:,1:nshiftk) = shiftk
 my_kptrlatt = kptrlatt

 call getkgrid_low(chksymbreak0,iout0,iscf2,kibz,kptopt,my_kptrlatt,kptrlen,&
   cryst%nsym,-1,nkibz,my_nshiftk,cryst%nsym,cryst%rprimd,my_shiftk,cryst%symafm,&
   cryst%symrel,vacuum0,wtk,indkpt,bz2ibz_smap,fullbz=kbz)

 if (present(bz2ibz)) then
   ABI_MOVE_ALLOC(bz2ibz_smap, bz2ibz)
 else
   ABI_SFREE(bz2ibz_smap)
 endif
 ABI_SFREE(indkpt)

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

!!****f* m_kpts/tetra_from_kptrlatt
!! NAME
!! tetra_from_kptrlatt
!!
!! FUNCTION
!!  Create an instance from kptrlatt and shiftk
!!
!! INPUTS
!!  cryst<cryst_t>=Crystalline structure.
!!  kptopt=Option for the k-point generation.
!!  kptrlatt(3,3)=k-point lattice specification
!!  nshiftk= number of shift vectors.
!!  shiftk(3,nshiftk)=shift vectors for k point generation
!!  nkibz=Number of points in the IBZ
!!  kibz(3,nkibz)=Reduced coordinates of the k-points in the IBZ.
!!  comm= MPI communicator
!!
!! OUTPUT
!!  tetra<htetra_t>=Tetrahedron object, fully initialized if ierr == 0.
!!  msg=Error message if ierr /= 0
!!  ierr=Exit status
!!
!! PARENTS
!!      gstate,wfk_analyze
!!
!! CHILDREN
!!      listkk,smpbz
!!
!! SOURCE

type(htetra_t) function tetra_from_kptrlatt( &
  cryst, kptopt, kptrlatt, nshiftk, shiftk, nkibz, kibz, comm, msg, ierr) result (htetra)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: kptopt,nshiftk,nkibz,comm
 integer,intent(out) :: ierr
 character(len=*),intent(out) :: msg
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 real(dp),intent(in) :: shiftk(3,nshiftk),kibz(3,nkibz)

!Local variables-------------------------------
!scalars
 integer :: nkfull,my_nkibz,new_nshiftk
 character(len=80) :: errorstring
!arrays
 integer :: new_kptrlatt(3,3)
 integer,allocatable :: indkk(:)
 integer,allocatable :: bz2ibz(:,:)
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
   my_nkibz, my_kibz, my_wtk, nkfull, kfull, new_kptrlatt=new_kptrlatt, new_shiftk=new_shiftk, bz2ibz=bz2ibz)

 ABI_FREE(my_wtk)
 new_nshiftk = size(new_shiftk, dim=2)

 if (my_nkibz /= nkibz .or. all(my_kibz /= kibz) ) then
   msg = sjoin("Input nkibz:", itoa(nkibz), "does not agree with computed value:", itoa(my_nkibz))
   ierr = 1; goto 10
 end if

 ! Do not support new_nshiftk > 1: lattice must be decomposed into boxes
 ! and this is not always possible (I think) with bizarre shifts
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

 rlatt = new_kptrlatt; call matr3inv(rlatt, klatt)

 ABI_MALLOC(indkk,(nkfull))
 indkk(:) = bz2ibz(1,:)
 ABI_SFREE(bz2ibz)
 call htetra_init(htetra, indkk, cryst%gprimd, klatt, kfull, nkfull, my_kibz, my_nkibz, ierr, errorstring, comm)
 if (ierr /= 0) msg = errorstring

 10 continue
 ABI_SFREE(my_kibz)
 ABI_SFREE(indkk)
 ABI_SFREE(kfull)
 ABI_SFREE(new_shiftk)

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
 character(len=500) :: msg
!arrays
 real(dp) :: ksym(3)

! *********************************************************************
 ierr = 0

 if(timrev/=1 .and. timrev/=0)then
   write(errmsg, '(3a,i0,a)' )&
    'timrev should be 0 or 1, while',ch10,&
    'it is equal to ',timrev,'.'
   ierr = 1; return
 end if

 if(nsym/=1)then
   ! Find the identity symmetry operation
   do isym=1,nsym
     tident=1
     do jj=1,3
       if(symrec(jj,jj,isym)/=1)tident=0
       do ii=1,3
         if( ii/=jj .and.symrec(ii,jj,isym)/=0)tident=0
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
!    k-point must be found, modulo time reversal if appropriate
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
         write(errmsg, '(a,a,a,i0,a,i0,a,a,a,a)' )&
          'k-point set must have full space-group symmetry',ch10,&
          'there is no match for kpt: ',ikpt,' transformed by symmetry: ',isym,ch10,&
          'Action: change kptopt to 2 or 3 and/or change or use shiftk',ch10,&
          'shiftk = 0 0 0 is always a safe choice.'
         ierr = 2; return
       end if

     end do ! End loop on isym
   end do ! End primary loop over k-points

   write(msg,'(a)')' symkchk : k-point set has full space-group symmetry.'
   call wrtout([std_out, ab_out], msg, 'COLL')
 end if

end function symkchk
!!***

!!****f* m_kpts/listkk
!! NAME
!! listkk
!!
!! FUNCTION
!! Given a list of nkpt1 initial k points kptns1 and a list of nkpt2
!! final k points kptns2, associates each final kpt with a "closest"
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
!!  symmat(3,3,nsym)=symmetry operations (symrel or symrec, depending on value of use_symrec
!!  timrev=1 if the use of time-reversal is allowed; 0 otherwise
!!  comm=MPI communicator.
!!  [exit_loop]: if present and True, exit the loop over k-points in the sphere as soon as the lenght**2 of the
!!    difference vector is smaller than tol12. Default: False
!!  [use_symrec]: if present and true, symmat assumed to be symrec, otherwise assumed to be symrel (default)
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
!!  and for the length of the vectors while the tolerance tol8 is used for
!!  the comparison of the squared lengths of the separate vectors.
!!
!! PARENTS
!!      initberry,initorbmag,inwffil,m_dvdb,m_ebands,m_eprenorms,m_exc_diago
!!      m_fock,m_fstab,m_haydock,m_ifc,m_kpts,m_phgamma,m_sigmaph,mlwfovlp_qp
!!
!! CHILDREN
!!      sort_dp,timab
!!
!! SOURCE

subroutine listkk(dksqmax,gmet,indkk,kptns1,kptns2,nkpt1,nkpt2,nsym,sppoldbl,symafm,symmat,timrev,comm, &
                  exit_loop, use_symrec) ! optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt1,nkpt2,nsym,sppoldbl,timrev,comm
 real(dp),intent(out) :: dksqmax
 logical,optional,intent(in) :: use_symrec, exit_loop
!arrays
 integer,intent(in) :: symafm(nsym),symmat(3,3,nsym)
 integer,intent(out) :: indkk(nkpt2*sppoldbl,6)
 real(dp),intent(in) :: gmet(3,3),kptns1(3,nkpt1),kptns2(3,nkpt2)

!Local variables-------------------------------
!scalars
 integer,parameter :: usesym=1, limit=1
 integer :: nprocs, my_rank, ierr, isk_start, isk_stop
 integer :: l3,ig1,ig2,ig3,ii,ikpg1,ikpt1,ikpt2,ikpt2_done, isk
 integer :: ilarger,ismaller,itrial
 integer :: isppol,isym,itimrev,jkpt1,jsym,jtime
 integer :: nsym_used,timrev_used
 real(dp) :: dksq,dksqmn,lk2,llarger,ldiff,lsmaller,ltrial,min_l
 real(dp) :: cpu,wall,gflops
 character(len=500) :: msg
!arrays
 integer :: dkint(3),jdkint(3),k1int(3),k2int(3)
 integer, allocatable :: isort(:), tmp_indkk(:,:)
 real(dp) :: tsec(2)
 real(dp) :: dk(3),kpg1(3),kpt1a(3),k1(3),k2(3)
 !real(dp) :: kasq,ka(3)
 real(dp),allocatable :: lkpg1(:),lkpg1_sorted(:)

! *************************************************************************

 call timab(1021, 1, tsec)
 call cwtime(cpu, wall, gflops, "start")

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 if (sppoldbl<1 .or. sppoldbl>2) then
   write(msg, '(a,i0,a)' )'The value of sppoldbl is: ',sppoldbl,', but it should be either 1 or 2.'
   MSG_BUG(msg)
 end if

 ! When usesym=0, the old way of converting the wavefunctions (without using the symmetries), is recovered.
 nsym_used=nsym
 timrev_used=timrev
 if(usesym==0)nsym_used=1
 if(usesym==0)timrev_used=0

 ! Precompute the length of the kpt1 vectors, also taking into account possible umklapp vectors
 l3 = (2*limit+1)**3
 ABI_CALLOC(lkpg1, (l3*nkpt1))
 ABI_CALLOC(lkpg1_sorted, (l3*nkpt1))
 ABI_MALLOC(isort, (l3*nkpt1))
 isort = 0

 call xmpi_split_work(nkpt1, comm, isk_start, isk_stop)
 !write(std_out,*)' List of kpt1 vectors'; write(std_out,*)' Length of the kpt1 vectors:'

!$OMP PARALLEL DO PRIVATE(k1, k1int, kpg1, ikpg1)
 do ikpt1=isk_start,isk_stop
   k1(:) = kptns1(:,ikpt1)  !; write(std_out,*)ikpt1,k1(:)
   k1int(:) = nint(k1(:) + tol12)
   k1(:) = k1(:) - k1int(:)
   do ig3=-limit,limit
     kpg1(3) = k1(3) + ig3
     do ig2=-limit,limit
       kpg1(2) = k1(2) + ig2
       do ig1=-limit,limit
         kpg1(1) = k1(1) + ig1

         ikpg1 = ig1 + limit + 1 + (2*limit+1)*(ig2+limit) + (2*limit+1)**2*(ig3+limit) + l3*(ikpt1-1)
         ! Compute the norm of the vector (also taking into account possible umklapp)
         lkpg1(ikpg1) = sqrt(gmet(1,1)*kpg1(1)**2+gmet(2,2)*kpg1(2)**2 + &
                             gmet(3,3)*kpg1(3)**2+two*(gmet(2,1)*kpg1(2)*kpg1(1) + &
                             gmet(3,2)*kpg1(3)*kpg1(2)+gmet(3,1)*kpg1(3)*kpg1(1)))
         lkpg1_sorted(ikpg1) = lkpg1(ikpg1)
         isort(ikpg1) = ikpg1
         !write(std_out,*)' ikpt1,ig1,ig2,ig3,lkpg1=',ikpt1,ig1,ig2,ig3,lkpg1(ikpg1)
       end do
     end do
   end do
 end do

 if (nprocs > 1) then
   call xmpi_sum(lkpg1_sorted, comm, ierr)
   call xmpi_sum(lkpg1, comm, ierr)
   call xmpi_sum(isort, comm, ierr)
 end if
 call cwtime_report(" listkk_loop1", cpu, wall, gflops)

 call sort_dp(l3*nkpt1, lkpg1_sorted, isort, tol12)
 ! From "precompute" to "sort_dp" represents more than 50% of the overall wall time for large meshes.
 call cwtime_report(" listkk_sort", cpu, wall, gflops)

 !write(std_out,*)' listkk : output list of kpt1 for checking purposes '
 !write(std_out,*)' ii,ikpt1,isort(ii)-l3*(ikpt1-1),lkpg1_sorted(ii),lkpg1(isort(ii)) '
 !do ii=1,l3*nkpt1
 !  ikpt1=(isort(ii)-1)/l3+1
 !  write(std_out,*)ii,ikpt1,isort(ii)-l3*(ikpt1-1),lkpg1_sorted(ii),lkpg1(isort(ii))
 !enddo

 dksqmax = zero
 indkk = 0
 ! TODO: Should change API to use this shape.
 ! workspace array for improved memory access.
 ABI_MALLOC(tmp_indkk, (6, nkpt2*sppoldbl))
 tmp_indkk = 0

 ! Split loop in contiguous blocks
 call xmpi_split_work(sppoldbl * nkpt2, comm, isk_start, isk_stop)

 do isppol=1,sppoldbl
   do ikpt2=1,nkpt2
     isk = ikpt2 + (isppol-1)*nkpt2
     if (isk < isk_start .or. isk > isk_stop) cycle

     ikpt2_done=0
     ! Precompute the length of the kpt2 vector, with the Umklapp vector such that it is the closest to the Gamma point
     k2(:)=kptns2(:,ikpt2)
     k2int(:)=nint(k2(:)+tol12)
     k2(:)=k2(:)-k2int(:)
     lk2=sqrt(gmet(1,1)*k2(1)**2+gmet(2,2)*k2(2)**2+&
              gmet(3,3)*k2(3)**2+two*(gmet(2,1)*k2(2)*k2(1)+&
              gmet(3,2)*k2(3)*k2(2)+gmet(3,1)*k2(3)*k2(1)))
     ! write(std_out, '(a,i4,7es16.6)' )' listkk : ikpt2,kptns2(:,ikpt2),k2(:),lk2=',ikpt2,kptns2(:,ikpt2),k2(:),lk2

     ! Find the kpt1 vector whose length is the most similar to the length of lk2 up to a tolerance.
     ! Use a bisection algorithm.
     ismaller=0; lsmaller=zero
     ilarger=l3*nkpt1+1; llarger=huge(one)

     ! This loop should never reach l3*nkpt1, since this is a bisection algorithm
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
     !write(std_out,*)' listkk : starting search at itrial=',itrial

     dksqmn=huge(one)

     ! The ii index is dummy. This avoids an infinite loop.
     do ii=1,l3*nkpt1
       ! If the difference in length between the trial vector and the target vector is bigger
       ! than the already achieved distance, the search is finished ...
       ldiff = abs(lkpg1_sorted(itrial) - lk2)
       ! write(std_out,*)' listkk : ii,itrial,lkpg1_sorted(itrial),lk2,ldiff,&
       ! dksqmn=',ii,itrial,lkpg1_sorted(itrial),lk2,ldiff,dksqmn

       if (ldiff**2 > dksqmn+tol8) exit

       ! If this k-point has already been examined in a previous batch, skip it
       ! First, compute the minimum of the difference of length of the sets of
       ! associated vectors thanks to Umklapp vectors with the target vector
       ikpt1 = (isort(itrial)-1) /l3 + 1
       min_l = minval(abs(lkpg1((ikpt1-1)*l3+1:(ikpt1-1)*l3+l3)-lk2))

       ! Then compare with the current ldiff
       ! write(std_out,*)' listkk : ikpt1,min_l,ldiff=',ikpt1,min_l,ldiff
       if (min_l > ldiff-tol12) then

         ! Now, will examine the trial vector, and the symmetric ones
         ! MG FIXME: Here there's a possible problem with the order of symmetries because
         ! in symkpt, time-reversal is the innermost loop. This can create inconsistencies in the symmetry tables.
         ! Besides, one should use symrel^{-1 T} to keep the correspondence between isym -> R or S
         do itimrev=0,timrev_used
           do isym=1,nsym_used
           !do itimrev=0,timrev_used

             ! Select magnetic characteristic of symmetries
             if (isppol == 1 .and. symafm(isym) == -1) cycle
             if (isppol == 2 .and. symafm(isym) == 1) cycle

             ! Compute symmetric point to kpt1
             if (usesym==1) then
               ! original code only used transpose(symrel)
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

             ! Compute difference with respect to kpt2, modulo a lattice vector
             dk(:)=kptns2(:,ikpt2)-kpt1a(:)
             if (usesym==1) then
               ! The tolerance insure similar behaviour on different platforms
               ! XG120418: Actually, *assumes* that the closest point will have reduced
               ! coordinates differing by less than 1/2. There might be elongated cells where this is not correct ...
               dkint(:)=nint(dk(:)+tol12)
               dk(:)=dk(:)-dkint(:)
             else
               dkint(:)=0
             end if

             ! Compute norm of the difference vector, and update kpt1 if better.
             dksq=gmet(1,1)*dk(1)**2+gmet(2,2)*dk(2)**2+ &
                  gmet(3,3)*dk(3)**2+two*(gmet(2,1)*dk(2)*dk(1)+ &
                  gmet(3,2)*dk(3)*dk(2)+gmet(3,1)*dk(3)*dk(1))

             if (dksq < dksqmn+tol8) then
               ! If exactly the right point (without using symmetries neither umklapp vector), will exit the search
               ! Note that in this condition, each coordinate is tested separately, without squaring.
               ! So, it is a much stronger condition than dksqmn < tol12
               if (sum(abs(kptns2(:,ikpt2)-kptns1(:,ikpt1)))<3*tol12) ikpt2_done = 1

               ! This line leads to a significant speedup for dense meshes but ~30 tests fail after this change.
               if (present(exit_loop)) then
                 if (exit_loop) then
                   if (dksq < tol12) ikpt2_done = 1
                 end if
               end if

               ! Update in three cases: either if succeeded to have exactly the vector, or the distance is better,
               ! or the distance is only slightly worsened so select the lowest itimrev, isym or ikpt1,
               ! in order to respect previous ordering
               if (ikpt2_done==1 .or. &
                  dksq+tol12<dksqmn .or. &
                  ( abs(dksq-dksqmn)<tol12 .and. &
                   ((itimrev<jtime) .or. &
                   (itimrev==jtime .and. isym<jsym) .or. &
                   (itimrev==jtime .and. isym==jsym .and. ikpt1<jkpt1))))then

                 dksqmn = dksq
                 jkpt1 = ikpt1
                 jsym = isym
                 jtime = itimrev
                 jdkint(:) = dkint(:)

                 !if (ikpt2_done == 1) then
                 !  write(std_out,*)'Succeeded to lower dskmn,ikpt2_done=',dksqmn,ikpt2_done
                 !  write(std_out,*)'  ikpt1,ikpt2=',ikpt1, ikpt2
                 !  write(std_out,*)'  ikpt1,isym,dkint(:),itimrev=',ikpt1,isym,dkint(:),itimrev
                 !  ka(:) = kpt1a(:) + dkint(:)
                 !  kasq=gmet(1,1)*ka(1)**2+gmet(2,2)*ka(2)**2+&
                 !       gmet(3,3)*ka(3)**2+two*(gmet(2,1)*ka(2)*ka(1)+&
                 !       gmet(3,2)*ka(3)*ka(2)+gmet(3,1)*ka(3)*ka(1))
                 !  write(std_out,*)'             k1 = ',kpt1a(:)
                 !  write(std_out,*)'          dkint = ',dkint(:)
                 !  write(std_out,*)'      Actual k1 = ',ka(:)
                 !  write(std_out,*)'             k2 = ',kptns2(:,ikpt2)
                 !  write(std_out,*)'      Actual k1sq = ',kasq
                 !end if
               end if
             end if

             if (ikpt2_done==1) exit
           end do ! isym
           if (ikpt2_done==1) exit
         end do ! itimrev
         if (ikpt2_done==1) exit
       end if

       ! Update the interval that has been explored
       if (itrial < ismaller) ismaller = itrial
       if (itrial > ilarger) ilarger = itrial

       ! Select the next index to be tried (preferably the smaller indices, but this is a bit arbitrary).
       ! write(std_out,*)' before choosing the next index :'
       ! write(std_out,*)' ismaller,itrial,ilarger=',ismaller,itrial,ilarger
       ! write(std_out,*)' lkpg1_sorted(ismaller-1),lk2,lkpg1_sorted(ilarger+1)=',&
       ! lkpg1_sorted(ismaller-1),lk2,lkpg1_sorted(ilarger+1)

       if (ismaller>1 .and. ilarger<l3*nkpt1) then
         if (abs(lkpg1_sorted(ismaller-1)-lk2) < abs(lkpg1_sorted(ilarger+1)-lk2)+tol12) then
           itrial = ismaller-1
         else
           itrial = ilarger+1
         end if
       end if
       if (ismaller==1 .and. ilarger<l3*nkpt1) itrial = ilarger+1
       if (ismaller>1 .and. ilarger==l3*nkpt1) itrial = ismaller-1
       !if(ismaller==1 .and. ilarger==l3*nkpt1), we are done with the loop !
     end do ! ikpt1

     ! Store indices (lots of cache miss here)
     !indkk(isk, 1) = jkpt1
     !indkk(isk, 2) = jsym
     !indkk(isk, 3:5) = jdkint(:)
     !indkk(isk, 6) = jtime

     tmp_indkk(1, isk) = jkpt1
     tmp_indkk(2, isk) = jsym
     tmp_indkk(3:5, isk) = jdkint(:)
     tmp_indkk(6, isk) = jtime

     dksqmax = max(dksqmax, dksqmn)

     if (dksqmn < -tol12) then
       write(msg, '(a,es16.6)' )'The minimum square of dk has negative norm: dksqmn= ',dksqmn
       MSG_BUG(msg)
     end if

     !write(std_out,'(a,i6,i2,2x,i6,5i3,es24.14)' )' listkk: ikpt2,isppol,indkk(isk,:)=',ikpt2,isppol,indkk(isk,:),dksqmn
   end do ! ikpt2
 end do ! isppol

 ABI_FREE(isort)
 ABI_FREE(lkpg1)
 ABI_FREE(lkpg1_sorted)

 indkk = transpose(tmp_indkk)
 ABI_FREE(tmp_indkk)
 if (nprocs > 1) then
   call xmpi_sum(indkk, comm, ierr)
   dksqmn = dksqmax
   call xmpi_max(dksqmn, dksqmax, comm, ierr)
 end if

 call timab(1021, 2, tsec)
 call cwtime_report(" listkk_end", cpu, wall, gflops)

end subroutine listkk
!!***

!!****f* m_kpts/getkgrid
!! NAME
!! getkgrid
!!
!! FUNCTION
!! Compute the grid of k points in the irreducible Brillouin zone.
!! Note that nkpt (and nkpthf) can be computed by calling this routine with nkpt=0, provided that kptopt/=0.
!! If downsampling is present, also compute a downsampled k grid.
!!
!! INPUTS
!! chksymbreak= if 1, will check whether the k point grid is symmetric (for kptopt=1,2 and 4), and stop if not.
!! iout=unit number for echoed output . 0 if no output is wished.
!! iscf= ( <= 0 =>non-SCF), >0 => SCF)  MG: FIXME I don't understand why we have to pass the value iscf.
!! kptopt=option for the generation of k points
!!   (defines whether spatial symmetries and/or time-reversal can be used)
!! msym=default maximal number of symmetries
!! nsym=number of symmetries
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym)=symmetry operations in real space in terms of primitive translations
!! vacuum(3)=for each direction, 0 if no vacuum, 1 if vacuum
!! [downsampling(3) = input variable that governs the downsampling]
!!
!! OUTPUT
!! kptrlen=length of the smallest real space supercell vector associated with the lattice of k points.
!! nkpt_computed=number of k-points in the IBZ computed in the present routine
!! If nkpt/=0  the following are also output:
!!   kpt(3,nkpt)=reduced coordinates of k points.
!!   wtk(nkpt)=weight assigned to each k point.
!! [fullbz(3,nkpt_fullbz)]=k-points generated in the full Brillouin zone.
!!   In output: allocated array with the list of k-points in the BZ.
!! [kpthf(3,nkpthf)]=k-points generated in the full Brillouin zone, possibly downsampled (for Fock).
!!
!! NOTES
!!  msym not needed since nsym is the last index.
!!
!! SIDE EFFECTS
!! Input/Output
!! nkpt=number of k points (might be zero, see output description)
!! kptrlatt(3,3)=k-point lattice specification
!! nshiftk=actual number of k-point shifts in shiftk
!! shiftk(3,MAX_NSHIFTK)=shift vectors for k point generation
!! [nkpthf] = number of k points in the full BZ, for the Fock operator.
!!
!! PARENTS
!!      ep_setupqpt,getshell,inkpts,inqpt,m_ab7_kpoints,m_bz_mesh,m_kpts
!!      nonlinear,testkgrid,thmeig
!!
!! CHILDREN
!!      mati3inv,matr3inv,metric,smallprim,smpbz,symkpt
!!
!! SOURCE

subroutine getkgrid(chksymbreak,iout,iscf,kpt,kptopt,kptrlatt,kptrlen,&
& msym,nkpt,nkpt_computed,nshiftk,nsym,rprimd,shiftk,symafm,symrel,vacuum,wtk,&
& fullbz,nkpthf,kpthf,downsampling) ! optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: chksymbreak,iout,iscf,kptopt,msym,nkpt,nsym
 integer,intent(inout),optional :: nkpthf
 integer,intent(inout) :: nshiftk
 integer,intent(inout) :: nkpt_computed !vz_i
 real(dp),intent(out) :: kptrlen
!arrays
 integer,intent(in) :: symafm(msym),symrel(3,3,msym),vacuum(3)
 integer,optional,intent(in) :: downsampling(3)
 integer,intent(inout) :: kptrlatt(3,3)
 integer,allocatable :: indkpt(:)
 integer,allocatable :: bz2ibz_smap(:,:)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: shiftk(3,MAX_NSHIFTK)
 real(dp),intent(inout) :: kpt(3,nkpt) !vz_i
 real(dp),intent(inout) :: wtk(nkpt)
 real(dp),optional,allocatable,intent(out) :: fullbz(:,:)
 real(dp),optional,intent(out) :: kpthf(:,:)

!Local variables-------------------------------
 real(dp),allocatable :: kpt_tmp(:,:), wtk_tmp(:)

 call getkgrid_low(chksymbreak,iout,iscf,kpt_tmp,kptopt,kptrlatt,kptrlen,&
   msym,nkpt,nkpt_computed,nshiftk,nsym,rprimd,shiftk,symafm,symrel,vacuum,wtk_tmp,indkpt,bz2ibz_smap,&
   fullbz,nkpthf,kpthf,downsampling)

 if (nkpt > 0) then
   kpt(:,1:nkpt) = kpt_tmp(:,1:nkpt)
   wtk(1:nkpt)   = wtk_tmp(1:nkpt)
 end if

 ABI_SFREE(kpt_tmp)
 ABI_SFREE(wtk_tmp)
 ABI_SFREE(indkpt)
 ABI_SFREE(bz2ibz_smap)

end subroutine getkgrid
!!***

!!****f* m_kpts/getkgrid_low
!! NAME
!! getkgrid_low
!!
!! FUNCTION
!! Compute the grid of k points in the irreducible Brillouin zone.
!! Note that nkpt (and nkpthf) can be computed by calling this routine with nkpt=0, provided that kptopt/=0.
!! If downsampling is present, also compute a downsampled k grid.
!!
!! INPUTS
!! chksymbreak= if 1, will check whether the k point grid is symmetric (for kptopt=1,2 and 4), and stop if not.
!! iout=unit number for echoed output . 0 if no output is wished.
!! iscf= ( <= 0 =>non-SCF), >0 => SCF)  MG: FIXME I don't understand why we have to pass the value iscf.
!! kptopt=option for the generation of k points
!!   (defines whether spatial symmetries and/or time-reversal can be used)
!! msym=default maximal number of symmetries
!! nsym=number of symmetries
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym)=symmetry operations in real space in terms of primitive translations
!! vacuum(3)=for each direction, 0 if no vacuum, 1 if vacuum
!! [downsampling(3) = input variable that governs the downsampling]
!!
!! OUTPUT
!! kptrlen=length of the smallest real space supercell vector associated with the lattice of k points.
!! nkpt_computed=number of k-points in the IBZ computed in the present routine
!! If nkpt/=0  the following are also output:
!!   kpt(3,nkpt)=reduced coordinates of k points.
!!   wtk(nkpt)=weight assigned to each k point.
!! bz2ibz_smap(nkbz, 6)= Mapping BZ --> IBZ.
!! [fullbz(3,nkpt_fullbz)]=k-points generated in the full Brillouin zone.
!!   In output: allocated array with the list of k-points in the BZ.
!! [kpthf(3,nkpthf)]=k-points generated in the full Brillouin zone, possibly downsampled (for Fock).
!!
!! NOTES
!!  msym not needed since nsym is the last index.
!!
!! SIDE EFFECTS
!! Input/Output
!! nkpt=number of k points (might be zero, see output description)
!! kptrlatt(3,3)=k-point lattice specification
!! nshiftk=actual number of k-point shifts in shiftk
!! shiftk(3,MAX_NSHIFTK)=shift vectors for k point generation
!! [nkpthf] = number of k points in the full BZ, for the Fock operator.
!!
!! PARENTS
!!      ep_setupqpt,getshell,inkpts,inqpt,m_ab7_kpoints,m_bz_mesh,m_kpts
!!      nonlinear,testkgrid,thmeig
!!
!! CHILDREN
!!      mati3inv,matr3inv,metric,smallprim,smpbz,symkpt
!!
!! SOURCE

subroutine getkgrid_low(chksymbreak,iout,iscf,kpt,kptopt,kptrlatt,kptrlen,&
& msym,nkpt,nkpt_computed,nshiftk,nsym,rprimd,shiftk,symafm,symrel,vacuum,wtk,indkpt,bz2ibz_smap,&
& fullbz,nkpthf,kpthf,downsampling) ! optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: chksymbreak,iout,iscf,kptopt,msym,nkpt,nsym
 integer,intent(inout),optional :: nkpthf
 integer,intent(inout) :: nshiftk
 integer,intent(inout) :: nkpt_computed !vz_i
 real(dp),intent(out) :: kptrlen
!arrays
 integer,intent(in) :: symafm(msym),symrel(3,3,msym),vacuum(3)
 integer,optional,intent(in) :: downsampling(3)
 integer,intent(inout) :: kptrlatt(3,3)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: shiftk(3,MAX_NSHIFTK)
 integer,allocatable,intent(out) :: indkpt(:)
 integer,allocatable,intent(out) :: bz2ibz_smap(:,:)
 real(dp),allocatable,intent(out) :: kpt(:,:) !vz_i
 real(dp),allocatable,intent(out) :: wtk(:)
 real(dp),optional,allocatable,intent(out) :: fullbz(:,:)
 real(dp),optional,intent(out) :: kpthf(:,:)

!Local variables-------------------------------
!scalars
 integer, parameter :: max_number_of_prime=47
 integer :: brav,decreased,found,ii,ikpt,iprime,ishiftk,isym,jshiftk,kshiftk,mkpt,mult
 integer :: nkpthf_computed,nkpt_fullbz,nkptlatt,nshiftk2,nsym_used,option
 integer :: test_prime,timrev
 integer :: nkpt_use
 real(dp) :: length2,ucvol,ucvol_super
 character(len=500) :: msg
!arrays
 integer, parameter :: prime_factor(max_number_of_prime)=(/2,3,5,7,9, 11,13,17,19,23,&
&  29,31,37,41,43, 47,53,59,61,67,&
&  71,73,79,83,89, 97,101,103,107,109,&
&  113,127,131,137,139, 149,151,157,163,167,&
&  173,179,181,191,193, 197,199/)
 integer :: kptrlatt2(3,3)
 integer,allocatable :: belong_chain(:),generator(:),number_in_chain(:)
 integer,allocatable :: repetition_factor(:),symrec(:,:,:)
! real(dp) :: cart(3,3)
 real(dp) :: dijk(3),delta_dmult(3),dmult(3),fact_vacuum(3),gmet(3,3)
 real(dp) :: gmet_super(3,3),gprimd(3,3),gprimd_super(3,3),klatt2(3,3)
 real(dp) :: klatt3(3,3),kptrlattr(3,3),ktransf(3,3),ktransf_invt(3,3)
 real(dp) :: metmin(3,3),minim(3,3),rmet(3,3),rmet_super(3,3),rprimd_super(3,3)
 real(dp),allocatable :: deltak(:,:),kpt_fullbz(:,:),shiftk2(:,:),shiftk3(:,:),spkpt(:,:),wtk_folded(:),wtk_fullbz(:)

! *************************************************************************

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 !call cwtime(cpu, wall, gflops, "start")
 if (kptopt==1.or.kptopt==4) then
! Cannot use antiferromagnetic symmetry operations to decrease the number of k points
!XG20191123 : now, antiferromagnetic symmetry operations can be used to decrease the number of k points for kptopt==4
   nsym_used=0
   do isym=1,nsym
     if(symafm(isym)==1 .or. kptopt==4)nsym_used=nsym_used+1
   end do
   ABI_MALLOC(symrec,(3,3,nsym_used))
   nsym_used=0
   do isym=1,nsym ! Get the symmetry matrices in terms of reciprocal basis
     if(symafm(isym)==1 .or. kptopt==4)then
       nsym_used=nsym_used+1
       call mati3inv(symrel(:,:,isym),symrec(:,:,nsym_used))
     end if
   end do
 else if (kptopt==2) then
   !Use only the time-reversal
   nsym_used=1
   ABI_MALLOC(symrec,(3,3,1))
   symrec(1:3,1:3,1)=0
   do ii=1,3
     symrec(ii,ii,1)=1
   end do
 end if

 kptrlatt2(:,:)=kptrlatt(:,:)
 nshiftk2=nshiftk
 ABI_MALLOC(shiftk2,(3,MAX_NSHIFTK))
 ABI_MALLOC(shiftk3,(3,MAX_NSHIFTK))
 shiftk2(:,:)=shiftk(:,:)

!Find a primitive k point lattice, if possible, by decreasing the number of shifts.
 if(nshiftk2/=1)then

   do
     ! Loop to be repeated if there has been a successful reduction of nshiftk2
     ABI_MALLOC(deltak,(3,nshiftk2))
     ABI_MALLOC(repetition_factor,(nshiftk2))
     ABI_MALLOC(generator,(nshiftk2))
     ABI_MALLOC(belong_chain,(nshiftk2))
     ABI_MALLOC(number_in_chain,(nshiftk2))

     decreased=0
     deltak(1,1:nshiftk2)=shiftk2(1,1:nshiftk2)-shiftk2(1,1)
     deltak(2,1:nshiftk2)=shiftk2(2,1:nshiftk2)-shiftk2(2,1)
     deltak(3,1:nshiftk2)=shiftk2(3,1:nshiftk2)-shiftk2(3,1)
     deltak(:,:)=deltak(:,:)-floor(deltak(:,:)+tol8)

!    Identify for each shift, the smallest repetition prime factor that yields a reciprocal lattice vector.
     repetition_factor(:)=0
     repetition_factor(1)=1
     do ishiftk=2,nshiftk2
       do iprime=1,max_number_of_prime
         test_prime=prime_factor(iprime)
         dmult(:)=test_prime*deltak(:,ishiftk)
         if(sum(abs( dmult(:)-nint(dmult(:)) ))<tol8)then
           repetition_factor(ishiftk)=test_prime
           exit
         end if
       end do
     end do

!    Initialize the selection of tentative generators
     generator(:)=1
     do ishiftk=1,nshiftk2
       if(repetition_factor(ishiftk)==0 .or. repetition_factor(ishiftk)==1)generator(ishiftk)=0
     end do

!    Try different shifts as generators, by order of increasing repetition factor,
!    provided they are equal or bigger than 2
     do iprime=1,max_number_of_prime
       do ishiftk=2,nshiftk2
         ! Note that ishiftk=1 is never a generator. It is the reference starting point.
         if(generator(ishiftk)==1 .and. repetition_factor(ishiftk)==prime_factor(iprime))then
!          Test the generator : is it indeed closed ?
           if(prime_factor(iprime)/=2)then
             do mult=2,prime_factor(iprime)-1
               dmult(:)=mult*deltak(:,ishiftk)
               found=0
               do jshiftk=1,nshiftk2
                 delta_dmult(:)=deltak(:,jshiftk)-dmult(:)
                 if(sum(abs(delta_dmult(:)-nint(delta_dmult(:)) ))<tol8)then
                   found=1
                   exit
                 end if
               end do
               if(found==0)exit
             end do
             if(found==0)generator(ishiftk)=0
           end if
           if(generator(ishiftk)==0)cycle
         else
           cycle
         end if
!        Now, test whether all k points can be found in all possible chains
         belong_chain(:)=0
         do jshiftk=1,nshiftk2
!          Initialize a chain starting from a k point not yet in a chain
           if(belong_chain(jshiftk)==0)then
             number_in_chain(:)=0   ! Not a member of the chain (yet)
             number_in_chain(jshiftk)=1   ! The first point in chain
             do mult=1,prime_factor(iprime)-1
               dmult(:)=mult*deltak(:,ishiftk)
               found=0
               do kshiftk=jshiftk+1,nshiftk2
                 delta_dmult(:)=deltak(:,kshiftk)-deltak(:,jshiftk)-dmult(:)
                 if(sum(abs(delta_dmult(:)-nint(delta_dmult(:)) ))<tol8)then
                   found=1
                   number_in_chain(kshiftk)=mult+1
                   exit
                 end if
               end do
               if(found==0)then
                 generator(ishiftk)=0
                 exit
               end if
             end do
             if(generator(ishiftk)==1)then
!              Store the chain
               do kshiftk=1,nshiftk2
                 if(number_in_chain(kshiftk)/=0)belong_chain(kshiftk)=number_in_chain(kshiftk)
               end do
             else
               exit
             end if
           end if
         end do

         if(generator(ishiftk)==0)cycle

!        For the generator based on ishiftk, all the k points have been found to belong to one chain.
!        All the initializing k points in the different chains have belong_chain(:)=1 .
!        They must be kept, and the others thrown away.
         ktransf(:,:)=0.0_dp
         ktransf(1,1)=1.0_dp
         ktransf(2,2)=1.0_dp
         ktransf(3,3)=1.0_dp
!        Replace one of the unit vectors by the shift vector deltak(:,ishiftk).
!        However, must pay attention not to make linear combinations.
!        Also, choose positive sign for first-non-zero value.
         if(abs(deltak(1,ishiftk)-nint(deltak(1,ishiftk)))>tol8)then
           if(deltak(1,ishiftk)>0)ktransf(:,1)= deltak(:,ishiftk)
           if(deltak(1,ishiftk)<0)ktransf(:,1)=-deltak(:,ishiftk)
         else if(abs(deltak(2,ishiftk)-nint(deltak(2,ishiftk)))>tol8)then
           if(deltak(2,ishiftk)>0)ktransf(:,2)= deltak(:,ishiftk)
           if(deltak(2,ishiftk)<0)ktransf(:,2)=-deltak(:,ishiftk)
         else if(abs(deltak(3,ishiftk)-nint(deltak(3,ishiftk)))>tol8)then
           if(deltak(3,ishiftk)>0)ktransf(:,3)= deltak(:,ishiftk)
           if(deltak(3,ishiftk)<0)ktransf(:,3)=-deltak(:,ishiftk)
         end if
!        Copy the integers to real(dp)
         kptrlattr(:,:)=kptrlatt2(:,:)
!        Go to reciprocal space
         call matr3inv(kptrlattr,klatt2)
!        Make the transformation
         do ii=1,3
           klatt3(:,ii)=ktransf(1,ii)*klatt2(:,1)+ktransf(2,ii)*klatt2(:,2)+ktransf(3,ii)*klatt2(:,3)
         end do
!        Back to real space
         call matr3inv(klatt3,kptrlattr)
!        real(dp) to integer
         kptrlatt2(:,:)=nint(kptrlattr(:,:))
!        Prepare the transformation of the shifts
         call matr3inv(ktransf,ktransf_invt)
         decreased=1
         kshiftk=0
         do jshiftk=1,nshiftk2
           if(belong_chain(jshiftk)==1)then
             kshiftk=kshiftk+1
!            Place the shift with index jshiftk in place of the one in kshiftk,
!            also transform the shift from the old to the new coordinate system
             shiftk3(:,kshiftk)=ktransf_invt(1,:)*shiftk2(1,jshiftk)+&
&             ktransf_invt(2,:)*shiftk2(2,jshiftk)+&
&             ktransf_invt(3,:)*shiftk2(3,jshiftk)
           end if
         end do
         nshiftk2=nshiftk2/prime_factor(iprime)
         shiftk2(:,1:nshiftk2)=shiftk3(:,1:nshiftk2)-floor(shiftk3(:,1:nshiftk2)+tol8)
         if(kshiftk/=nshiftk2)then
           MSG_BUG('The search for a primitive k point lattice contains a bug.')
         end if

!        If this trial shift was successful, must exit the loop on trial ishiftk,
!        and reinitialize the global loop
         if(decreased==1)exit
       end do ! ishiftk
       if(decreased==1)exit
     end do ! iprime

     ABI_FREE(belong_chain)
     ABI_FREE(deltak)
     ABI_FREE(number_in_chain)
     ABI_FREE(repetition_factor)
     ABI_FREE(generator)

     if(decreased==0 .or. nshiftk2==1)exit

   end do ! Infinite loop

 end if !  End nshiftk being 1 or larger

!Impose shiftk coordinates to be in [0,1[
 do ishiftk=1,nshiftk2
   do ii=1,3
     if(shiftk2(ii,ishiftk)>one-tol8) shiftk2(ii,ishiftk)=shiftk2(ii,ishiftk)-1.0_dp
     if(shiftk2(ii,ishiftk)<-tol8)    shiftk2(ii,ishiftk)=shiftk2(ii,ishiftk)+1.0_dp
   end do
 end do

!Compute the number of k points in the G-space unit cell
 nkptlatt=kptrlatt2(1,1)*kptrlatt2(2,2)*kptrlatt2(3,3) &
& +kptrlatt2(1,2)*kptrlatt2(2,3)*kptrlatt2(3,1) &
& +kptrlatt2(1,3)*kptrlatt2(2,1)*kptrlatt2(3,2) &
& -kptrlatt2(1,2)*kptrlatt2(2,1)*kptrlatt2(3,3) &
& -kptrlatt2(1,3)*kptrlatt2(2,2)*kptrlatt2(3,1) &
& -kptrlatt2(1,1)*kptrlatt2(2,3)*kptrlatt2(3,2)

!Check whether the number of k points is positive, otherwise, change the handedness of kptrlatt2
 if(nkptlatt<=0)then
   ! write(std_out,*)' getkgrid : nkptlatt is negative !'
   kptrlatt2(:,3)=-kptrlatt2(:,3)
   nkptlatt=-nkptlatt
   do ishiftk=1,nshiftk2
     shiftk2(3,ishiftk)=-shiftk2(3,ishiftk)
   end do
 end if

!Determine the smallest supercell R-vector whose contribution
!is not taken correctly into account in the k point integration.
!Increase enormously the size of the cell when vacuum is present.
 fact_vacuum(:)=1
 if(vacuum(1)==1)fact_vacuum(1)=1000.0_dp
 if(vacuum(2)==1)fact_vacuum(2)=1000.0_dp
 if(vacuum(3)==1)fact_vacuum(3)=1000.0_dp
 do ii=1,3
   rprimd_super(:,ii)=fact_vacuum(1)*rprimd(:,1)*kptrlatt2(1,ii)+&
&   fact_vacuum(2)*rprimd(:,2)*kptrlatt2(2,ii)+&
&   fact_vacuum(3)*rprimd(:,3)*kptrlatt2(3,ii)
 end do

 call metric(gmet_super,gprimd_super,-1,rmet_super,rprimd_super,ucvol_super)
 call smallprim(metmin,minim,rprimd_super)
 length2=min(metmin(1,1),metmin(2,2),metmin(3,3))
 kptrlen=sqrt(length2)

 !write(msg,'(a,es16.6)' )' getkgrid : length of smallest supercell vector (bohr)=',kptrlen
 !call wrtout(std_out,msg,'COLL')
! If the number of shifts has been decreased, determine the set of kptrlatt2 vectors
! with minimal length (without using fact_vacuum)
! It is worth to determine the minimal set of vectors so that the kptrlatt that is output
! does not seem screwy, although correct but surprising.
 if(nshiftk/=nshiftk2)then
   do ii=1,3
     rprimd_super(:,ii)=rprimd(:,1)*kptrlatt2(1,ii)+rprimd(:,2)*kptrlatt2(2,ii)+rprimd(:,3)*kptrlatt2(3,ii)
   end do
   call metric(gmet_super,gprimd_super,-1,rmet_super,rprimd_super,ucvol_super)
!  Shift vectors in cartesian coordinates (reciprocal space)
   do ishiftk=1,nshiftk2
     shiftk3(:,ishiftk)=gprimd_super(:,1)*shiftk2(1,ishiftk)+&
&     gprimd_super(:,2)*shiftk2(2,ishiftk)+&
&     gprimd_super(:,3)*shiftk2(3,ishiftk)
   end do
   call smallprim(metmin,minim,rprimd_super)
   call metric(gmet_super,gprimd_super,-1,rmet_super,minim,ucvol_super)
!  This is the new kptrlatt2
   do ii=1,3
     dijk(:)=gprimd(1,:)*minim(1,ii)+&
&     gprimd(2,:)*minim(2,ii)+&
&     gprimd(3,:)*minim(3,ii)
     kptrlatt2(:,ii)=nint(dijk(:))
   end do
!  Shifts in the new set of kptrlatt vectors
   do ishiftk=1,nshiftk2
     shiftk2(:,ishiftk)=minim(1,:)*shiftk3(1,ishiftk)+&
&     minim(2,:)*shiftk3(2,ishiftk)+&
&     minim(3,:)*shiftk3(3,ishiftk)
   end do
 end if

!brav=1 is able to treat all bravais lattices.
 brav=1
 mkpt=nkptlatt*nshiftk2

 ABI_MALLOC(spkpt,(3,mkpt))
 option=0
 if(iout/=0)option=1

 !call cwtime_report(' shifts', cpu, wall, gflops)

 if (present(downsampling))then
   call smpbz(brav,iout,kptrlatt2,mkpt,nkpthf_computed,nshiftk2,option,shiftk2,spkpt,downsampling=downsampling)
   if (present(kpthf) .and. nkpthf/=0) then
     ! Returns list of k-points in the Full BZ, possibly downsampled for Fock
     kpthf = spkpt(:,1:nkpthf)
   end if
   nkpthf=nkpthf_computed
 end if

 call smpbz(brav,iout,kptrlatt2,mkpt,nkpt_fullbz,nshiftk2,option,shiftk2,spkpt)
 !call cwtime_report(' smpbz', cpu, wall, gflops)

 if(kptopt==1 .or. kptopt==2 .or. kptopt==4)then

   ABI_MALLOC(indkpt,(nkpt_fullbz))
   ABI_MALLOC(kpt_fullbz,(3,nkpt_fullbz))
   ABI_MALLOC(bz2ibz_smap, (6, nkpt_fullbz))
#if 1
   ABI_MALLOC(wtk_fullbz,(nkpt_fullbz))
   ABI_MALLOC(wtk_folded,(nkpt_fullbz))

   kpt_fullbz(:,:)=spkpt(:,1:nkpt_fullbz)
   wtk_fullbz(1:nkpt_fullbz)=1.0_dp/dble(nkpt_fullbz)

   timrev=1;if (kptopt==4) timrev=0

   indkpt = 0
   call symkpt(chksymbreak,gmet,indkpt,iout,kpt_fullbz,nkpt_fullbz,&
&   nkpt_computed,nsym_used,symrec,timrev,wtk_fullbz,wtk_folded,bz2ibz_smap,xmpi_comm_self)

   ABI_FREE(symrec)
   ABI_FREE(wtk_fullbz)

   !do ikpt=1,nkpt_fullbz
   !  write(*,*) ikpt, indkpt(ikpt), bz2ibz_smap(1,ikpt), indkpt(bz2ibz_smap(1,ikpt))
   !end do
#else
   kpt_fullbz(:,:)=spkpt(:,1:nkpt_fullbz)

   timrev=1;if (kptopt==4) timrev=0

   call symkpt_new(chksymbreak,gmet,indkpt,iout,kpt_fullbz,nkpt_fullbz,&
&   nkpt_computed,nsym_used,symrec,timrev,bz2ibz_smap,xmpi_comm_self)

   ABI_FREE(symrec)
   ABI_CALLOC(wtk_folded,(nkpt_fullbz))
   do ii=1,nkpt_fullbz
    ikpt = indkpt(bz2ibz_smap(1,ii))
    wtk_folded(ikpt) = wtk_folded(ikpt) + one
   end do
   wtk_folded = wtk_folded / nkpt_fullbz
#endif

 else if(kptopt==3)then
   ABI_CALLOC(bz2ibz_smap, (6, nkpt_fullbz))
   bz2ibz_smap(1,:) = [(ii,ii=1,nkpt_fullbz)]
   bz2ibz_smap(2,:) = 1 !isym
   nkpt_computed=nkpt_fullbz
 end if
 !call cwtime_report(' symkpt', cpu, wall, gflops)

!The number of k points has been computed from kptopt, kptrlatt, nshiftk, shiftk,
!and the eventual symmetries, it is presently called nkpt_computed.
 nkpt_use = nkpt
 if (nkpt<0) nkpt_use = nkpt_computed

!Check that the argument nkpt is coherent with nkpt_computed, if nkpt/=0.
 if(nkpt_use/=nkpt_computed .and. nkpt/=0)then
   write(msg, '(a,i0,5a,i0,7a)') &
&   'The argument nkpt = ',nkpt_use,', does not match',ch10,&
&   'the number of k points generated by kptopt, kptrlatt, shiftk,',ch10,&
&   'and the eventual symmetries, that is, nkpt= ',nkpt_computed,'.',ch10,&
&   'However, note that it might be due to the user,',ch10,&
&   'if nkpt is explicitely defined in the input file.',ch10,&
&   'In this case, please check your input file.'
   MSG_BUG(msg)
 end if

 ABI_MALLOC(kpt,(3,nkpt_use))
 ABI_MALLOC(wtk,(nkpt_use))

 if(kptopt==1 .or. kptopt==2 .or. kptopt==4)then

   if(nkpt_use/=0)then
     do ikpt=1,nkpt_use
       kpt(:,ikpt)=kpt_fullbz(:,indkpt(ikpt))
       if(iscf>=0 .or. iscf==-3 .or. iscf==-1.or.iscf==-2)wtk(ikpt)=wtk_folded(indkpt(ikpt))
     end do
   end if

   if (present(fullbz)) then
     ! Returns list of k-points in the Full BZ.
     ABI_MOVE_ALLOC(kpt_fullbz,fullbz)
   else
     ABI_FREE(kpt_fullbz)
   end if

   ABI_FREE(wtk_folded)

 else if(kptopt==3)then

   if(nkpt_use/=0)then
     kpt(:,1:nkpt_use)=spkpt(:,1:nkpt_use)
     if(iscf>1 .or. iscf==-3 .or. iscf==-1.or.iscf==-2)wtk(1:nkpt_use)=1.0_dp/dble(nkpt_use)
   end if

   if (present(fullbz)) then
     ! Returns list of k-points in the Full BZ.
     ABI_MALLOC(fullbz,(3,nkpt_fullbz))
     fullbz = spkpt(:,1:nkpt_fullbz)
   end if

 end if

 ABI_FREE(spkpt)
 kptrlatt(:,:)=kptrlatt2(:,:)
 nshiftk=nshiftk2
 shiftk(:,1:nshiftk)=shiftk2(:,1:nshiftk)
 ABI_FREE(shiftk2)
 ABI_FREE(shiftk3)

end subroutine getkgrid_low
!!***

!!****f* m_kpts/get_full_kgrid
!! NAME
!! get_full_kgrid
!!
!! FUNCTION
!! Create full grid of kpoints and find equivalent
!! irred ones. Duplicates work in getkgrid, but need all outputs of kpt_fullbz, and indkpt
!!
!! INPUTS
!!  kpt(3,nkpt)=irreducible kpoints
!!  kptrlatt(3,3)=lattice vectors for full kpoint grid
!!  nkpt=number of irreducible kpoints
!!  nkpt_fullbz=number of kpoints in full brillouin zone
!!  nshiftk=number of kpoint grid shifts
!!  nsym=number of symmetries
!!  shiftk(3,nshiftk)=kpoint shifts
!!  symrel(3,3,nsym)=symmetry matrices in real space
!!
!! OUTPUT
!!  indkpt(nkpt_fullbz)=non-symmetrized indices of the k-points (see symkpt.f)
!!  kpt_fullbz(3,nkpt_fullbz)=kpoints in full brillouin zone
!!
!! NOTES
!!  MG: The present implementation always assumes kptopt==1 !!!!
!!
!! TODO: This routine should be removed
!!
!! PARENTS
!!      m_phonons
!!
!! CHILDREN
!!      destroy_kptrank,get_kpt_fullbz,get_rank_1kpt,mati3inv,mkkptrank
!!
!! SOURCE

subroutine get_full_kgrid(indkpt,kpt,kpt_fullbz,kptrlatt,nkpt,nkpt_fullbz,nshiftk,nsym,shiftk,symrel)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,nkpt_fullbz,nshiftk,nsym
!arrays
 integer,intent(in) :: kptrlatt(3,3),symrel(3,3,nsym)
 integer,intent(out) :: indkpt(nkpt_fullbz)
 real(dp),intent(in) :: kpt(3,nkpt),shiftk(3,nshiftk)
 real(dp),intent(out) :: kpt_fullbz(3,nkpt_fullbz)

!Local variables-------------------------------
!scalars
 integer :: ikpt,isym,itim,timrev
 integer :: symrankkpt
 character(len=500) :: msg
 type(krank_t) :: krank
!arrays
 integer :: inv_symrel(3,3,nsym)
 real(dp) :: k2(3)

! *********************************************************************

!Invert symrels => gives symrels for kpoints

 do isym=1,nsym
   call mati3inv (symrel(:,:,isym),inv_symrel(:,:,isym))
 end do

 call get_kpt_fullbz(kpt_fullbz,kptrlatt,nkpt_fullbz,nshiftk,shiftk)

 ! make full k-point rank arrays
 krank = krank_new(nkpt, kpt)

 !find equivalence to irred kpoints in kpt
 indkpt(:) = 0
 timrev=1 ! includes the time inversion symmetry
 do ikpt=1,nkpt_fullbz
   do isym=1,nsym
     do itim=1,(1-2*timrev),-2

       k2(:) = itim*(inv_symrel(:,1,isym)*kpt_fullbz(1,ikpt) + &
                     inv_symrel(:,2,isym)*kpt_fullbz(2,ikpt) + &
                     inv_symrel(:,3,isym)*kpt_fullbz(3,ikpt))

       symrankkpt = krank%get_rank(k2)
       if (krank%invrank(symrankkpt) /= -1) indkpt(ikpt) = krank%invrank(symrankkpt)

     end do ! loop time reversal symmetry
   end do !  loop sym ops

   if (indkpt(ikpt) == 0) then
     write(msg,'(a,i0)')' indkpt(ikpt) is still 0: no irred kpoint is equiv to ikpt ',ikpt
     MSG_BUG(msg)
   end if
 end do !  loop full kpts

 call krank%free()

end subroutine get_full_kgrid
!!***

!!****f* m_kpts/get_kpt_fullbz
!! NAME
!! get_kpt_fullbz
!!
!! FUNCTION
!! Create full grid of kpoints from kptrlatt and shiftk
!!
!! INPUTS
!!  kptrlatt(3,3)=lattice vectors for full kpoint grid
!!  nkpt_fullbz=number of kpoints in full brillouin zone
!!  nshiftk=number of kpoint grid shifts
!!  shiftk(3,nshiftk)=kpoint shifts
!!
!! OUTPUT
!!  kpt_fullbz(3,nkpt_fullbz)=kpoints in full brillouin zone
!!
!! PARENTS
!!      get_full_kgrid,invars2
!!
!! CHILDREN
!!      mati3det,matr3inv,wrap2_pmhalf
!!
!! SOURCE

subroutine get_kpt_fullbz(kpt_fullbz,kptrlatt,nkpt_fullbz,nshiftk,shiftk)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt_fullbz,nshiftk
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 real(dp),intent(in) :: shiftk(3,nshiftk)
 real(dp),intent(out) :: kpt_fullbz(3,nkpt_fullbz)

!Local variables-------------------------------
!scalars
 integer, parameter :: max_number_of_prime=47
 integer :: det,ii,ikshft,iprim,jj,kk,nn
 character(len=500) :: msg
!arrays
 integer :: boundmax(3),boundmin(3),common_factor(3)
 integer, parameter :: prime_factor(max_number_of_prime)=(/2,3,5,7,9, 11,13,17,19,23,&
&  29,31,37,41,43, 47,53,59,61,67,&
&  71,73,79,83,89, 97,101,103,107,109,&
&  113,127,131,137,139, 149,151,157,163,167,&
&  173,179,181,191,193, 197,199/)
 real(dp) :: k1(3),k2(3),klatt(3,3),rlatt(3,3),shift(3),test_rlatt(3,3)

! *********************************************************************

!Identify first factors that can be used to rescale the three kptrlatt vectors
!Only test a large set of prime factors, though ...
 do jj=1,3
   common_factor(jj)=1
   rlatt(:,jj)=kptrlatt(:,jj)
   do iprim=1,max_number_of_prime
     test_rlatt(:,jj)=rlatt(:,jj)/dble(prime_factor(iprim))
!    If one of the components is lower than 1 in absolute value, then it is not worth to continue the search.
     if(minval(abs(abs(test_rlatt(:,jj))-half))<half-tol8)exit
     do
       if(sum(abs(test_rlatt(:,jj)-nint(test_rlatt(:,jj)) ))<tol8)then
         common_factor(jj)=prime_factor(iprim)*common_factor(jj)
         rlatt(:,jj)=rlatt(:,jj)/dble(prime_factor(iprim))
         test_rlatt(:,jj)=test_rlatt(:,jj)/dble(prime_factor(iprim))
       else
         exit
       end if
     end do
   end do
 end do
 call mati3det(kptrlatt,det)
 det=det/(common_factor(1)*common_factor(2)*common_factor(3))

 rlatt(:,:)=kptrlatt(:,:)
 call matr3inv(rlatt,klatt)
!Now, klatt contains the three primitive vectors of the k lattice,
!in reduced coordinates. One builds all k vectors that
!are contained in the first Brillouin zone, with coordinates
!in the interval [0,1[ . First generate boundaries of a big box.
!In order to generate all possible vectors in the reciprocal space,
!one must consider all multiples of the primitive ones, until a vector with only integers is found.
!The maximum bound is the scale of the corresponding kptrlatt vector, times the determinant of kptrlatt. Also consider negative vectors.
!On this basis, compute the bounds.
 do jj=1,3
!  To accomodate the shifts, boundmin starts from -1
!  Well, this is not a complete solution ...
   boundmin(jj)=-1-common_factor(jj)*abs(det)
   boundmax(jj)=common_factor(jj)*abs(det)
 end do

 nn=1
 do kk=boundmin(3),boundmax(3)
   do jj=boundmin(2),boundmax(2)
     do ii=boundmin(1),boundmax(1)
       do ikshft=1,nshiftk

!        Coordinates of the trial k point with respect to the k primitive lattice
         k1(1)=ii+shiftk(1,ikshft)
         k1(2)=jj+shiftk(2,ikshft)
         k1(3)=kk+shiftk(3,ikshft)

!        Reduced coordinates of the trial k point
         k2(:)=k1(1)*klatt(:,1)+k1(2)*klatt(:,2)+k1(3)*klatt(:,3)

!        Eliminate the point if outside [0,1[
         if(k2(1)<-tol10)cycle ; if(k2(1)>one-tol10)cycle
         if(k2(2)<-tol10)cycle ; if(k2(2)>one-tol10)cycle
         if(k2(3)<-tol10)cycle ; if(k2(3)>one-tol10)cycle

!        Wrap the trial values in the interval ]-1/2,1/2] .
         call wrap2_pmhalf(k2(1),k1(1),shift(1))
         call wrap2_pmhalf(k2(2),k1(2),shift(2))
         call wrap2_pmhalf(k2(3),k1(3),shift(3))
         if(nn > nkpt_fullbz) then
           write (msg,'(a,i0)')' nkpt_fullbz mis-estimated, exceed nn=',nn
           MSG_BUG(msg)
         end if
         kpt_fullbz(:,nn)=k1(:)
         nn=nn+1
       end do
     end do
   end do
 end do
 nn = nn-1

 if (nn /= nkpt_fullbz) then
   write (msg,'(2(a,i0),a,a)')' nkpt_fullbz= ',nkpt_fullbz,' underestimated  nn=',nn,&
&   ch10, "Perhaps your k grid or shifts do not correspond to the symmetry?"
   MSG_BUG(msg)
 end if

end subroutine get_kpt_fullbz
!!***

!!****f* m_kpts/smpbz
!! NAME
!! smpbz
!!
!! FUNCTION
!! Generate a set of special k (or q) points which samples in a homogeneous way
!! the entire Brillouin zone of a simple lattice, face-centered cubic,
!! body-centered lattice and hexagonal lattice.
!! If kptrlatt is diagonal, the algorithm used here reduces to the usual
!! Monkhorst-Pack set of k points.
!!
!! INPUTS
!!  brav = 1 or -1 -> simple lattice; 2 -> face-centered cubic;
!!   3 -> body-centered lattice; 4 -> hexagonal lattice (D6h)
!!  downsampling(3) [optional, for brav=1 only]
!!    Three integer numbers, describing the downsampling of the k grid
!!    If present, in any case, only the first shiftk is taken into account
!!    The absolute value of one number gives, for the corresponding k-coordinate, the factor of decrease of the sampling
!!    If zero, only one point is used to sample along this direction
!!    The sign has also a meaning :
!!    - if three numbers are negative, perform a face-centered sampling
!!    - if two numbers are negative, perform a body-centered sampling
!!    - if one number is negative, perform a face-centered sampling for the two-dimensional lattice of the other directions
!!    - if one number is zero and at least one number is negative, perform face-centered sampling for the non-zero directions.
!!  iout = unit number for output
!!  kptrlatt(3,3)=integer coordinates of the primitive vectors of the
!!   lattice reciprocal to the k point lattice to be generated here
!!   If diagonal, the three values are the Monkhorst-Pack usual values, in case of simple cubic.
!!  mkpt = maximum number of k points
!!  nshiftk= number of shift vectors in the repeated cell
!!  option= Flag defining what will be printed of iout: 0 for k points, anything else for q points.
!!    Also, for q points, if the Gamma point is present, place it first in the list.
!!  shiftk(3,nshiftk) = vectors that will be used to determine the shifts from (0. 0. 0.).
!!
!! OUTPUT
!!  nkpt = number of k points
!!  spkpt(3,mkpt) = the nkpt first values contain the special k points
!!   obtained by the Monkhorst & Pack method, in reduced coordinates.
!!   These vectors have to be multiplied by the reciprocal basis vectors
!!   gprimd(3,3) (in cartesian coordinates) to obtain the special k points
!!   set in cartesian coordinates.
!!
!! NOTES
!!  also allows for more than one vector in repeated cell.
!!  this routine should be rewritten, to use the Wigner-Seitz cell,
!!  and thus unify the different treatments.
!!  References :
!!  H.J. Monkhorst and J.D. Pack, Phys. Rev. B 13, 5188 (1976) [[cite:Monkhorst1976]]
!!  J.D. Pack and H.J. Monkhorst, Phys. Rev. B 16, 1748 (1977) [[cite:Pack1977]]
!!  A.H. MacDonald, Phys. Rev. B 18, 5897 (1978) [[cite:MacDonald1978]]
!!  R.A. Evarestov and V.P. Smirnov, Phys. Stat. Sol. (b) 119, 9 (1983) [[cite:Evarestov1983]]
!!
!! PARENTS
!!      ep_setupqpt,getkgrid,harmonic_thermo,initberry,initorbmag,m_fstab,m_ifc
!!      m_tdep_abitypes
!!
!! CHILDREN
!!      matr3inv,wrap2_pmhalf,wrtout
!!
!! SOURCE

subroutine smpbz(brav,iout,kptrlatt,mkpt,nkpt,nshiftk,option,shiftk,spkpt,downsampling)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: brav,iout,mkpt,nshiftk,option
 integer,intent(out) :: nkpt
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 integer,optional,intent(in) :: downsampling(3)
 real(dp),intent(in) :: shiftk(3,nshiftk)
 real(dp),intent(out) :: spkpt(3,mkpt)

!Local variables -------------------------
!scalars
 integer,parameter :: prtvol=0
 integer :: dividedown,ii,ikshft,jj,kk,nkpout,nkptlatt,nn,proddown
 real(dp) :: shift
 character(len=500) :: msg
!arrays
 integer :: ads(3),boundmax(3),boundmin(3),cds(3),coord(3),ngkpt(3)
 integer, allocatable :: found1(:,:),found2(:,:),found3(:,:)
 real(dp) :: k1(3),k2(3),kcar(3),klatt(3,3),ktest(3),rlatt(3,3)

! *********************************************************************

!DEBUG
!write(std_out,*)' smpbz : brav,iout,mkpt,nkpt,option=',brav,iout,mkpt,nkpt,option
!write(std_out,*)' smpbz : kptrlatt(:,:)=',kptrlatt(:,:)
!write(std_out,*)' smpbz : nshiftk=',nshiftk
!write(std_out,*)' smpbz : shiftk(:,:)=',shiftk(:,:)
!write(std_out,*)' smpbz : downsampling(:)=',downsampling(:)
!ENDDEBUG

 if(option/=0) call wrtout(iout,'       Homogeneous q point set in the B.Z.  ','COLL')

 if(abs(brav)/=1)then
!  Only generate Monkhorst-Pack lattices
   if(kptrlatt(1,2)/=0 .or. kptrlatt(2,1)/=0 .or. &
&   kptrlatt(1,3)/=0 .or. kptrlatt(3,1)/=0 .or. &
&   kptrlatt(2,3)/=0 .or. kptrlatt(3,2)/=0     ) then
     write(msg, '(2a,a,3i0,a,a,3i4,a,a,3i4)' )&
&     'When abs(brav)/=1, kptrlatt must be diagonal, while it is',ch10,&
&     'kptrlatt(:,1)= ',kptrlatt(:,1),ch10,&
&     'kptrlatt(:,2)= ',kptrlatt(:,2),ch10,&
&     'kptrlatt(:,3)= ',kptrlatt(:,3)
     MSG_BUG(msg)
   end if

   ngkpt(1)=kptrlatt(1,1)
   ngkpt(2)=kptrlatt(2,2)
   ngkpt(3)=kptrlatt(3,3)
!
   if( (ngkpt(1)<=0.or.ngkpt(2)<=0.or.ngkpt(3)<=0) .and. (ngkpt(1)/=0.or.ngkpt(2)/=0.or.ngkpt(3)/=0) ) then
     write(msg, '(5a,i4,a,a,i0,a,a,i0,a,a)' )&
&     'All ngkpt (or ngqpt) must be strictly positive',ch10,&
&     'or all ngk(q)pt must be zero (for Gamma sampling), but :',ch10,&
&     'ngk(q)pt(1) = ',ngkpt(1),ch10,&
&     'ngk(q)pt(2) = ',ngkpt(2),ch10,&
&     'ngk(q)pt(3) = ',ngkpt(3),ch10,&
&     'Action: correct ngkpt or ngqpt in the input file.'
     MSG_BUG(msg)
   end if
 end if

!Just in case the user wants the grid downsampled to the Gamma point, checks that it is present, and possibly exits
 if(present(downsampling))then
   if(sum(abs(downsampling(:)))==0)then
     do ikshft=1,nshiftk
       if(sum(abs(shiftk(:,ikshft)))>tol12)cycle
       nkpt=1
       spkpt(:,1)=zero
       return
     end do
   end if
 end if

!*********************************************************************

 if(abs(brav)==1)then

!  Compute the number of k points in the G-space unit cell
!  (will be multiplied by nshiftk later).
   nkptlatt=kptrlatt(1,1)*kptrlatt(2,2)*kptrlatt(3,3) &
&   +kptrlatt(1,2)*kptrlatt(2,3)*kptrlatt(3,1) &
&   +kptrlatt(1,3)*kptrlatt(2,1)*kptrlatt(3,2) &
&   -kptrlatt(1,2)*kptrlatt(2,1)*kptrlatt(3,3) &
&   -kptrlatt(1,3)*kptrlatt(2,2)*kptrlatt(3,1) &
&   -kptrlatt(1,1)*kptrlatt(2,3)*kptrlatt(3,2)

   if(present(downsampling))then
     if(.not.(downsampling(1)==1 .and. downsampling(2)==1 .and. downsampling(3)==1))then
       if(nshiftk>1)then
         write(msg, '(a,3i4,2a,i4,4a)' )&
&         'Real downsampling is activated, with downsampling(1:3)=',downsampling(1:3),ch10,&
&         'However, nshiftk must be 1 in this case, while the input nshiftk=',nshiftk,ch10,&
&         'Action: either choose not to downsample the k point grid (e.g. fockdownsampling=1),',ch10,&
&         'or set nshiftk=1.'
         MSG_ERROR(msg)
       end if
       proddown=downsampling(1)*downsampling(2)*downsampling(3)
       if(proddown/=0)then
         dividedown=abs(proddown)
         if(minval(downsampling(:))<0)then   ! If there is at least one negative number
           dividedown=dividedown*2
           if(proddown>0)dividedown=dividedown*2 ! If there are two negative numbers
         end if
       end if
       if(mod(nkptlatt,dividedown)==0)then
         nkptlatt=nkptlatt/dividedown
       else
         write(msg, '(a,3i4,2a,i4,4a)' )&
&         'The requested downsampling, with downsampling(1:3)=',downsampling(1:3),ch10,&
&         'is not compatible with kptrlatt=',ch10,&
&         kptrlatt(:,:),ch10,&
&         'that gives nkptlatt=',nkptlatt,ch10,&
&         'Action: either choose not to downsample the k point grid (e.g. fockdownsampling=1),',ch10,&
&         'or modify your k-point grid and/or your downsampling in order for them to be compatible.'
         MSG_ERROR(msg)
       end if
     end if
   end if

!  Simple Lattice
   if (prtvol > 0) call wrtout(std_out,'       Simple Lattice Grid ','COLL')
   if (mkpt<nkptlatt*nshiftk) then
     write(msg, '(a,a,a,i8,a,a,a,a,a)' )&
&     'The value of mkpt is not large enough. It should be',ch10,&
&     'at least',nkptlatt*nshiftk,',',ch10,&
&     'Action: set mkpt to that value in the main routine,',ch10,&
&     'and recompile the code.'
     MSG_BUG(msg)
   end if

!  Build primitive vectors of the k lattice
   rlatt(:,:)=kptrlatt(:,:)
   call matr3inv(rlatt,klatt)

!  write(std_out,*)' First primitive vector of the k lattice :',klatt(:,1)
!  write(std_out,*)' Second primitive vector of the k lattice :',klatt(:,2)
!  write(std_out,*)' Third primitive vector of the k lattice :',klatt(:,3)

!  Now, klatt contains the three primitive vectors of the k lattice,
!  in reduced coordinates. One builds all k vectors that
!  are contained in the first Brillouin zone, with coordinates
!  in the interval [0,1[ . First generate boundaries of a big box.

   do jj=1,3

!    Mathematically, one has to find the coordinates of the corners of a
!    rectangular paralleliped with integer coordinates, that multiplies the klatt primitive cell and allows
!    it to incorporate completely the [0,1]^3 box. Then take the minimum and maximum
!    of these coordinates, and round them negatively and positively to the next integer.
!    This can be done easily using kptrlatt, considering each coordinate in turn
!    and boils down to enlarging the boundaries for jj by the value of kptrlatt(:,jj),
!    acting on boundmin or boundmax depending on the sign ot kptrlatt(:,jj).
!    XG171020 The coding before 171020 was correct, despite being very simple.
     boundmin(jj)=0 ; boundmax(jj)=0
     do ii=1,3
       if(kptrlatt(ii,jj)<0)boundmin(jj)=boundmin(jj)+kptrlatt(ii,jj)
       if(kptrlatt(ii,jj)>0)boundmax(jj)=boundmax(jj)+kptrlatt(ii,jj)
     end do

!    To accomodate the shifts, boundmin and boundmax don't start from 0, but are enlarged by one
!    positively and/or negatively.
!    XG171020 Coding in v8.6.0 and before was not correct. This one is even simpler actually.
     boundmin(jj)=boundmin(jj)-ceiling(maxval(shiftk(jj,:))+tol14)
     boundmax(jj)=boundmax(jj)-floor(minval(shiftk(jj,:))-tol14)

   end do

   if(present(downsampling))then
     ABI_MALLOC(found1,(boundmin(2):boundmax(2),boundmin(3):boundmax(3)))
     ABI_MALLOC(found2,(boundmin(1):boundmax(1),boundmin(3):boundmax(3)))
     ABI_MALLOC(found3,(boundmin(1):boundmax(1),boundmin(2):boundmax(2)))
     found1=0 ; found2=0 ; found3=0
   end if

   nn=1
   do kk=boundmin(3),boundmax(3)
     coord(3)=kk
     do jj=boundmin(2),boundmax(2)
       coord(2)=jj
       do ii=boundmin(1),boundmax(1)
         coord(1)=ii

!        Here, apply the downsampling : skip some of the trials
         if(present(downsampling))then

           if(downsampling(1)==0 .and. found1(coord(2),coord(3))==1)cycle
           if(downsampling(2)==0 .and. found2(coord(1),coord(3))==1)cycle
           if(downsampling(3)==0 .and. found3(coord(1),coord(2))==1)cycle

           ads(:)=abs(downsampling(:))
           if(ads(1)>0 .and. mod(coord(1),ads(1))/=0)cycle
           if(ads(2)>0 .and. mod(coord(2),ads(2))/=0)cycle
           if(ads(3)>0 .and. mod(coord(3),ads(2))/=0)cycle
           cds(:)=coord(:)/ads(:)
           if(minval(downsampling(:))<0)then   ! If there is at least one negative number

             if(downsampling(1)*downsampling(2)*downsampling(3)/=0)then  ! If there is no zero number
!              Face-centered case
               if(downsampling(1)<0 .and. downsampling(2)<0 .and. downsampling(3)<0)then ! All three are negative
                 if(mod(sum(cds(:)),2)/=0)cycle
!              One-face-centered case
               else if(downsampling(1)*downsampling(2)*downsampling(3)<0)then  ! Only one is negative
                 if(downsampling(1)<0 .and. mod(cds(2)+cds(3),2)/=0)cycle
                 if(downsampling(2)<0 .and. mod(cds(1)+cds(3),2)/=0)cycle
                 if(downsampling(3)<0 .and. mod(cds(1)+cds(2),2)/=0)cycle
!              Body-centered case ! What is left : two are negative
               else
                 ! Either all are zero, or all are one, so skip when sum is 1 or 2.
                 if(sum(mod(cds(:),2))==1 .or. sum(mod(cds(:),2))==2)cycle
               end if
             else
               if(downsampling(1)==0 .and. mod(cds(2)+cds(3),2)/=0)cycle
               if(downsampling(2)==0 .and. mod(cds(1)+cds(3),2)/=0)cycle
               if(downsampling(3)==0 .and. mod(cds(1)+cds(2),2)/=0)cycle
             end if
           end if
         end if

         do ikshft=1,nshiftk

!          Only the first shiftk is taken into account if downsampling
!          if(.false.)then
           if(present(downsampling))then
             if(.not.(downsampling(1)==1 .and. downsampling(2)==1 .and. downsampling(3)==1))then
               if(ikshft>1)cycle
             end if
           end if

!          Coordinates of the trial k point with respect to the k primitive lattice
           k1(1)=ii+shiftk(1,ikshft)
           k1(2)=jj+shiftk(2,ikshft)
           k1(3)=kk+shiftk(3,ikshft)
!          Reduced coordinates of the trial k point
           k2(:)=k1(1)*klatt(:,1)+k1(2)*klatt(:,2)+k1(3)*klatt(:,3)
!          Eliminate the point if outside [0,1[
           if(k2(1)<-tol10)cycle ; if(k2(1)>one-tol10)cycle
           if(k2(2)<-tol10)cycle ; if(k2(2)>one-tol10)cycle
           if(k2(3)<-tol10)cycle ; if(k2(3)>one-tol10)cycle
!          Wrap the trial values in the interval ]-1/2,1/2] .
           call wrap2_pmhalf(k2(1),k1(1),shift)
           call wrap2_pmhalf(k2(2),k1(2),shift)
           call wrap2_pmhalf(k2(3),k1(3),shift)
           spkpt(:,nn)=k1(:)
           nn=nn+1

           if(present(downsampling))then
             found1(coord(2),coord(3))=1
             found2(coord(1),coord(3))=1
             found3(coord(1),coord(2))=1
           end if

         end do
       end do
     end do
   end do
   nkpt=nn-1

   if(present(downsampling))then
     ABI_FREE(found1)
     ABI_FREE(found2)
     ABI_FREE(found3)
   end if

   if(nkpt/=nkptlatt*nshiftk)then
     write(msg, '(a,i0,3a,i0,a)' )&
     'The number of k points ',nkpt,' is not equal to',ch10,&
     'nkptlatt*nshiftk which is ',nkptlatt*nshiftk,'.'
     MSG_BUG(msg)
   end if

 else if(brav==2)then

!  Face-Centered Lattice
   if (prtvol > 0) call wrtout(std_out,'       Face-Centered Lattice Grid ','COLL')
   if (mkpt<ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk/2) then
     write(msg, '(a,a,a,i0,a,a,a,a,a)' )&
&     'The value of mkpt is not large enough. It should be',ch10,&
&     'at least',(ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk)/2,',',ch10,&
&     'Action: set mkpt to that value in the main routine,',ch10,&
&     'and recompile the code.'
     MSG_BUG(msg)
   end if
   nn=1
   if (ngkpt(1)/=ngkpt(2).or.ngkpt(1)/=ngkpt(3)) then
     write(msg, '(4a,3(a,i0,a),a)' )&
&     'For face-centered lattices, the numbers ngqpt(1:3)',ch10,&
&     'must be equal, while they are :',ch10,&
&     'ngqpt(1) = ',ngkpt(1),ch10,&
&     'ngqpt(2) = ',ngkpt(2),ch10,&
&     'ngqpt(3) = ',ngkpt(3),ch10,&
&     'Action: modify ngqpt(1:3) in the input file.'
     MSG_BUG(msg)
   end if
   if ((ngkpt(1)*nshiftk)/=(((ngkpt(1)*nshiftk)/2)*2)) then
     write(msg, '(4a,3(a,i0,a),a)' )&
&     'For face-centered lattices, the numbers ngqpt(1:3)*nshiftk',ch10,&
&     'must be even, while they are :',ch10,&
&     'ngqpt(1)*nshiftk = ',ngkpt(1)*nshiftk,ch10,&
&     'ngqpt(2)*nshiftk = ',ngkpt(2)*nshiftk,ch10,&
&     'ngqpt(3)*nshiftk = ',ngkpt(3)*nshiftk,ch10,&
&     'Action: modify ngqpt(1:3)*nshiftk in the input file.'
     MSG_ERROR(msg)
   end if
   if (ngkpt(1)==0.or.ngkpt(2)==0.or.ngkpt(3)==0) then
     spkpt(1,1)=0.0_dp
     spkpt(2,1)=0.0_dp
     spkpt(3,1)=0.0_dp
     nkpt=1
   else
     do kk=1,ngkpt(3)
       do jj=1,ngkpt(2)
         do ii=1,ngkpt(1)
           do ikshft=1,nshiftk
             k1(1)=(ii-1+shiftk(1,ikshft))/ngkpt(1)
             k1(2)=(jj-1+shiftk(2,ikshft))/ngkpt(2)
             k1(3)=(kk-1+shiftk(3,ikshft))/ngkpt(3)
!            Wrap the trial values in the interval ]-1/2,1/2] .
             call wrap2_pmhalf(k1(1),k2(1),shift)
             call wrap2_pmhalf(k1(2),k2(2),shift)
             call wrap2_pmhalf(k1(3),k2(3),shift)
!            Test whether it is inside the FCC BZ.
             ktest(1)=2*k2(1)-1.0d-10
             ktest(2)=2*k2(2)-2.0d-10
             ktest(3)=2*k2(3)-5.0d-10
             if (abs(ktest(1))+abs(ktest(2))+abs(ktest(3))<1.5_dp) then
               kcar(1)=ktest(1)+1.0d-10
               kcar(2)=ktest(2)+2.0d-10
               kcar(3)=ktest(3)+5.0d-10
               spkpt(1,nn)=0.5_dp*kcar(2)+0.5_dp*kcar(3)
               spkpt(2,nn)=0.5_dp*kcar(1)+0.5_dp*kcar(3)
               spkpt(3,nn)=0.5_dp*kcar(1)+0.5_dp*kcar(2)
               nn=nn+1
             end if
           end do
         end do
       end do
     end do
     nkpt=nn-1
     if(nkpt/=ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk/2)then
       write(msg, '(a,i8,a,a,a,i8,a)' )&
&       'The number of k points ',nkpt,'  is not equal to',ch10,&
&       '(ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk)/2 which is',&
&       (ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk)/2,'.'
       MSG_BUG(msg)
     end if
   end if

 else if(brav==3)then

!  Body-Centered Lattice (not mandatory cubic !)
   if (prtvol > 0) call wrtout(std_out,'       Body-Centered Lattice Grid ','COLL')
   if (mkpt<ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk/4) then
     write(msg, '(a,a,a,i8,a,a,a,a,a)' )&
&     'The value of mkpt is not large enough. It should be',ch10,&
&     'at least',(ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk)/4,',',ch10,&
&     'Action: set mkpt to that value in the main routine,',ch10,&
&     'and recompile the code.'
     MSG_BUG(msg)
   end if
   nn=1
   if ((ngkpt(1)*nshiftk)/=(((ngkpt(1)*nshiftk)/2)*2) .or.&
&   (ngkpt(2)*nshiftk)/=(((ngkpt(2)*nshiftk)/2)*2) .or.&
&   (ngkpt(3)*nshiftk)/=(((ngkpt(3)*nshiftk)/2)*2) ) then
     write(msg, '(4a,3(a,i6,a),a)' )&
&     'For body-centered lattices, the numbers ngqpt(1:3)',ch10,&
&     'must be even, while they are :',ch10,&
&     'ngqpt(1)*nshiftk = ',ngkpt(1)*nshiftk,ch10,&
&     'ngqpt(2)*nshiftk = ',ngkpt(2)*nshiftk,ch10,&
&     'ngqpt(3)*nshiftk = ',ngkpt(3)*nshiftk,ch10,&
&     'Action: modify ngqpt(1:3) in the input file.'
     MSG_ERROR(msg)
   end if
   if (ngkpt(1)==0.or.ngkpt(2)==0.or.ngkpt(3)==0) then
     spkpt(1,1)=0.0_dp
     spkpt(2,1)=0.0_dp
     spkpt(3,1)=0.0_dp
     nkpt=1
   else
     do kk=1,ngkpt(3)
       do jj=1,ngkpt(2)
         do ii=1,ngkpt(1)
           do ikshft=1,nshiftk
             k1(1)=(ii-1+shiftk(1,ikshft))/ngkpt(1)
             k1(2)=(jj-1+shiftk(2,ikshft))/ngkpt(2)
             k1(3)=(kk-1+shiftk(3,ikshft))/ngkpt(3)
!            Wrap the trial values in the interval ]-1/2,1/2] .
             call wrap2_pmhalf(k1(1),k2(1),shift)
             call wrap2_pmhalf(k1(2),k2(2),shift)
             call wrap2_pmhalf(k1(3),k2(3),shift)
!            Test whether it is inside the BCC BZ.
             ktest(1)=2*k2(1)-1.0d-10
             ktest(2)=2*k2(2)-2.0d-10
             ktest(3)=2*k2(3)-5.0d-10
             if (abs(ktest(1))+abs(ktest(2))<1._dp) then
               if (abs(ktest(1))+abs(ktest(3))<1._dp) then
                 if (abs(ktest(2))+abs(ktest(3))<1._dp) then
                   kcar(1)=ktest(1)+1.0d-10
                   kcar(2)=ktest(2)+2.0d-10
                   kcar(3)=ktest(3)+5.0d-10
                   spkpt(1,nn)=-0.5*kcar(1)+0.5*kcar(2)+0.5*kcar(3)
                   spkpt(2,nn)=0.5*kcar(1)-0.5*kcar(2)+0.5*kcar(3)
                   spkpt(3,nn)=0.5*kcar(1)+0.5*kcar(2)-0.5*kcar(3)
                   nn=nn+1
                 end if
               end if
             end if
           end do
         end do
       end do
     end do
     nkpt=nn-1
     if(nkpt==0)then
       write(msg, '(3a)' )&
&       'BCC lattice, input ngqpt=0, so no kpt is generated.',ch10,&
&       'Action: modify ngqpt(1:3) in the input file.'
       MSG_ERROR(msg)
     end if
     if(nkpt/=(ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk)/4)then
       write(msg, '(a,i0,3a,i0,a)' )&
&       'The number of k points ',nkpt,' is not equal to',ch10,&
&       '(ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk)/4 which is',(ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk)/4,'.'
       MSG_BUG(msg)
     end if
   end if

 else if(brav==4)then

!  Hexagonal Lattice  (D6h)
   if (prtvol > 0) call wrtout(std_out,'       Hexagonal Lattice Grid ','COLL')
   if (mkpt<ngkpt(1)*ngkpt(2)*ngkpt(3)) then
     write(msg, '(a,a,a,i0,a,a,a,a,a)' )&
&     'The value of mkpt is not large enough. It should be',ch10,&
&     'at least',ngkpt(1)*ngkpt(2)*ngkpt(3),',',ch10,&
&     'Action: set mkpt to that value in the main routine,',ch10,&
&     'and recompile the code.'
     MSG_BUG(msg)
   end if
   nn=1
   if (ngkpt(1)/=ngkpt(2)) then
     write(msg, '(4a,2(a,i0,a),a)' )&
&     'For hexagonal lattices, the numbers ngqpt(1:2)',ch10,&
&     'must be equal, while they are:',ch10,&
&     'ngqpt(1) = ',ngkpt(1),ch10,&
&     'ngqpt(2) = ',ngkpt(2),ch10,&
&     'Action: modify ngqpt(1:3) in the input file.'
     MSG_ERROR(msg)
   end if
   if (ngkpt(1)==0.or.ngkpt(2)==0.or.ngkpt(3)==0) then
     write(msg, '(3a)' )&
&     'For hexagonal lattices, ngqpt(1:3)=0 is not permitted',ch10,&
&     'Action: modify ngqpt(1:3) in the input file.'
     MSG_ERROR(msg)
   else
     do kk=1,ngkpt(3)
       do jj=1,ngkpt(2)
         do ii=1,ngkpt(1)
           do ikshft=1,nshiftk
             k1(1)=(ii-1+shiftk(1,ikshft))/ngkpt(1)
             k1(2)=(jj-1+shiftk(2,ikshft))/ngkpt(2)
             k1(3)=(kk-1+shiftk(3,ikshft))/ngkpt(3)
!            Wrap the trial values in the interval ]-1/2,1/2] .
             call wrap2_pmhalf(k1(1),k2(1),shift)
             call wrap2_pmhalf(k1(2),k2(2),shift)
             call wrap2_pmhalf(k1(3),k2(3),shift)
             spkpt(:,nn)=k2(:)
             nn=nn+1
           end do
         end do
       end do
     end do
     nkpt=nn-1
     if(nkpt/=ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk)then
       write(msg, '(a,i0,3a,i0,a)' )&
&       'The number of k points ',nkpt,'  is not equal to',ch10,&
&       'ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk which is',ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk,'.'
       MSG_BUG(msg)
     end if
   end if

 else

   write(msg, '(a,i0,a,a,a)' )&
&   'The calling routine asks brav= ',brav,'.',ch10,&
&   'but only brav=1 or -1,2,3 or 4 are allowed.'
   MSG_BUG(msg)
 end if

 if (option/=0) then
!  Put the Gamma point first
   if(nkpt>1)then
     do ii=1,nkpt
       if(sum(abs(spkpt(:,ii)))<tol8)then
         spkpt(:,ii)=spkpt(:,1)
         spkpt(:,1)=zero
         exit
       end if
     end do
   end if

   write(msg,'(a,i8)')' Grid q points  : ',nkpt
   call wrtout(iout,msg,'COLL')
   nkpout=nkpt
   if(nkpt>80)then
     call wrtout(iout,' greater than 80, so only write 20 of them ','COLL')
     nkpout=20
   end if
   do ii=1,nkpout
     write(msg, '(1x,i2,a2,3es16.8)' )ii,') ',spkpt(1,ii),spkpt(2,ii),spkpt(3,ii)
     call wrtout(iout,msg,'COLL')
   end do
 end if

end subroutine smpbz
!!***

!!****f* m_kpts/testkgrid
!! NAME
!! testkgrid
!!
!! FUNCTION
!! Test different grids of k points. The algorithm used is based on the idea of testing different
!! one-dimensional sets of possible k point grids. It is not exhaustive (other families could be included),
!! but should do a respectable job in all cases. The Monkhorst-Pack set of grids (defined with respect to
!! symmetry axes, and not primitive axes) is always tested.
!!
!! INPUTS
!!  bravais(11): bravais(1)=iholohedry
!!               bravais(2)=center
!!               bravais(3:11)=coordinates of rprim in the axes of the conventional bravais lattice (*2 if center/=0)
!!  iout=unit number for echoed output
!!  msym=default maximal number of symmetries
!!  nsym=number of symmetries
!!  prtkpt=if non-zero, will write the characteristics of k grids, then stop
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrel(3,3,nsym)=symmetry operations in real space in terms of primitive translations
!!  vacuum(3)=for each direction, 0 if no vacuum, 1 if vacuum
!!
!! OUTPUT
!!  kptrlatt(3,3)=k-point lattice specification
!!  nshiftk=number of k-point shifts in shiftk (always 1 from this routine)
!!  shiftk(3,MAX_NSHIFTK)=shift vectors for k point generation
!!
!! SIDE EFFECTS
!!  kptrlen=length of the smallest real space supercell vector associated with the lattice of k points.
!!
!! NOTES
!! Note that nkpt can be computed by calling this routine with input value nkpt=0
!! Note that kptopt is always =1 in this routine.
!!
!! PARENTS
!!      inkpts,m_ab7_kpoints
!!
!! CHILDREN
!!      getkgrid,abi_abort,matr3inv,metric,smallprim,wrtout,xmpi_abort
!!
!! SOURCE

subroutine testkgrid(bravais,iout,kptrlatt,kptrlen,&
& msym,nshiftk,nsym,prtkpt,rprimd,shiftk,symafm,symrel,vacuum)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,msym,nsym,prtkpt
 integer,intent(out) :: nshiftk
 real(dp),intent(inout) :: kptrlen
!arrays
 integer,intent(in) :: bravais(11),symafm(msym),symrel(3,3,msym),vacuum(3)
 integer,intent(out) :: kptrlatt(3,3)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: shiftk(3,MAX_NSHIFTK) !vz_i

!Local variables-------------------------------
!scalars
 integer,parameter :: kptopt=1,mkpt_list=100000
 integer :: ang90,center,dirvacuum,equal,igrid,igrid_current,iholohedry,ii,init_mult,iscale,iscf
 integer :: iset,mult1,mult2,mult3,ndims,nkpt,nkpt_current,nkpt_trial,nset
 real(dp) :: buffer_scale,determinant,fact,factor,kptrlen_current,kptrlen_max,kptrlen_target
 real(dp) :: kptrlen_trial,length1,length2,length3,length_axis1,length_axis2
 real(dp) :: length_axis3,merit_factor,mult1h,mult2h,mult3h,reduceda,reducedb
 real(dp) :: sca,scb,scc,surface,ucvol
 character(len=500) :: msg
!arrays
 integer :: kptrlatt_current(3,3),kptrlatt_trial(3,3)
 integer,allocatable :: grid_list(:)
 real(dp) :: axes(3,3),gmet(3,3),gprimd(3,3),matrix1(3,3),matrix2(3,3)
 real(dp) :: metmin(3,3),minim(3,3),r2d(3,3),rmet(3,3),rsuper(3,3)
 real(dp) :: shiftk_current(3,MAX_NSHIFTK),shiftk_trial(3,MAX_NSHIFTK)
 real(dp),allocatable :: kpt(:,:),kptrlen_list(:),wtk(:)

! *************************************************************************

 kptrlen_target=kptrlen

!The vacuum array must be made of 0 or 1
 do ii=1,3
   if(vacuum(ii)/=0 .and. vacuum(ii)/=1)then
     write(msg,'(a,a,a,i1,a,i3,a,a)')&
&     'The values of vacuum must be 0 or 1.',ch10,&
&     'However, the input vacuum(',ii,') is',vacuum(ii),ch10,&
&     'Action: correct vacuum in your input file.'
     MSG_ERROR(msg)
   end if
 end do

!Specific preparation for 2-dimensional system
 if(sum(vacuum(:))==1)then

!  Make the non-active vector orthogonal to the active vectors,
!  and take it along the z direction
   if(vacuum(1)==1)then
     r2d(1,3)=rprimd(2,2)*rprimd(3,3)-rprimd(3,2)*rprimd(2,3)
     r2d(2,3)=rprimd(3,2)*rprimd(1,3)-rprimd(1,2)*rprimd(3,3)
     r2d(3,3)=rprimd(1,2)*rprimd(2,3)-rprimd(2,2)*rprimd(1,3)
     r2d(:,1)=rprimd(:,2)
     r2d(:,2)=rprimd(:,3)
     dirvacuum=1
   else if(vacuum(2)==1)then
     r2d(1,3)=rprimd(2,3)*rprimd(3,1)-rprimd(3,3)*rprimd(2,1)
     r2d(2,3)=rprimd(3,3)*rprimd(1,1)-rprimd(1,3)*rprimd(3,1)
     r2d(3,3)=rprimd(1,3)*rprimd(2,1)-rprimd(2,3)*rprimd(1,1)
     r2d(:,1)=rprimd(:,3)
     r2d(:,2)=rprimd(:,1)
     dirvacuum=2
   else if(vacuum(3)==1)then
     r2d(1,3)=rprimd(2,1)*rprimd(3,2)-rprimd(3,1)*rprimd(2,2)
     r2d(2,3)=rprimd(3,1)*rprimd(1,2)-rprimd(1,1)*rprimd(3,2)
     r2d(3,3)=rprimd(1,1)*rprimd(2,2)-rprimd(2,1)*rprimd(1,2)
     r2d(:,1)=rprimd(:,1)
     r2d(:,2)=rprimd(:,2)
     dirvacuum=3
   end if
   surface=sqrt(sum(r2d(:,3)**2))
!  Identify the 2-D Bravais lattice
!  DEBUG
!  write(std_out,*)' r2d=',r2d(:,:)
!  ENDDEBUG
   call metric(gmet,gprimd,-1,rmet,r2d,ucvol)
   call smallprim(metmin,minim,r2d)
!  DEBUG
!  write(std_out,*)' minim=',minim(:,:)
!  ENDDEBUG
   ang90=0 ; equal=0 ; center=0
   axes(:,:)=minim(:,:)
   if(abs(metmin(1,2))<tol8)ang90=1
   if(abs(metmin(1,1)-metmin(2,2))<tol8)equal=1
   if(ang90==1)then
     if(equal==1)iholohedry=4
     if(equal==0)iholohedry=2
   else if(equal==1)then
     reduceda=metmin(1,2)/metmin(1,1)
     if(abs(reduceda+0.5_dp)<tol8)then
       iholohedry=3
     else if(abs(reduceda-0.5_dp)<tol8)then
       iholohedry=3
!      Use conventional axes
       axes(:,2)=minim(:,2)-minim(:,1)
     else
       iholohedry=2 ; center=1
       axes(:,1)=minim(:,1)+minim(:,2)
       axes(:,2)=minim(:,2)-minim(:,1)
     end if
   else
     reduceda=metmin(1,2)/metmin(1,1)
     reducedb=metmin(1,2)/metmin(2,2)
     if(abs(reduceda+0.5_dp)<tol8)then
       iholohedry=2 ; center=1
       axes(:,2)=2.0_dp*minim(:,2)+minim(:,1)
     else if(abs(reduceda-0.5_dp)<tol8)then
       iholohedry=2 ; center=1
       axes(:,2)=2.0_dp*minim(:,2)-minim(:,1)
     else if(abs(reducedb+0.5_dp)<tol8)then
       iholohedry=2 ; center=1
       axes(:,1)=2.0_dp*minim(:,1)+minim(:,2)
     else if(abs(reducedb-0.5_dp)<tol8)then
       iholohedry=2 ; center=1
       axes(:,1)=2.0_dp*minim(:,1)-minim(:,2)
     else
       iholohedry=1
     end if
   end if
!  Make sure that axes form a right-handed coordinate system
   determinant=axes(1,1)*axes(2,2)*axes(3,3) &
&   +axes(1,2)*axes(2,3)*axes(3,1) &
&   +axes(1,3)*axes(3,2)*axes(2,1) &
&   -axes(1,1)*axes(3,2)*axes(2,3) &
&   -axes(1,3)*axes(2,2)*axes(3,1) &
&   -axes(1,2)*axes(2,1)*axes(3,3)
   if(determinant<zero)then
     axes(:,1)=-axes(:,1)
   end if
!  Prefer symmetry axes on the same side as the primitive axes
   sca=axes(1,1)*r2d(1,1)+axes(2,1)*r2d(2,1)+axes(3,1)*r2d(3,1)
   scb=axes(1,2)*r2d(1,2)+axes(2,2)*r2d(2,2)+axes(3,2)*r2d(3,2)
   scc=axes(1,3)*rprimd(1,dirvacuum)&
&   +axes(2,3)*rprimd(2,dirvacuum)&
&   +axes(3,3)*rprimd(3,dirvacuum)
   if(sca<-tol8 .and. scb<-tol8)then
     axes(:,1)=-axes(:,1) ; sca=-sca
     axes(:,2)=-axes(:,2) ; scb=-scb
   end if
!  Doing this might change the angle between vectors, so that
!  the cell is not conventional anymore
!  if(sca<-tol8 .and. scc<-tol8)then
!  axes(:,1)=-axes(:,1) ; sca=-sca
!  axes(:,3)=-axes(:,3) ; scc=-scc
!  end if
!  if(scb<-tol8 .and. scc<-tol8)then
!  axes(:,2)=-axes(:,2) ; scb=-scb
!  axes(:,3)=-axes(:,3) ; scc=-scc
!  end if
   length_axis1=sqrt(axes(1,1)**2+axes(2,1)**2+axes(3,1)**2)
   length_axis2=sqrt(axes(1,2)**2+axes(2,2)**2+axes(3,2)**2)

!  DEBUG
!  write(std_out,*)' testkgrid: iholohedry, center =',iholohedry,center
!  write(std_out,*)' testkgrid: axis 1=',axes(:,1)
!  write(std_out,*)' testkgrid: axis 2=',axes(:,2)
!  write(std_out,*)' testkgrid: axis 3=',axes(:,3)
!  write(std_out,*)' testkgrid: length_axis=',length_axis1,length_axis2
!  ENDDEBUG

!  End special treatment of 2-D case
 end if

!3-dimensional system
 if(sum(vacuum(:))==0)then
   iholohedry=bravais(1)
   center=bravais(2)
   fact=1.0_dp
   if(center/=0)fact=0.5_dp
   matrix1(:,1)=bravais(3:5)*fact
   matrix1(:,2)=bravais(6:8)*fact
   matrix1(:,3)=bravais(9:11)*fact
   call matr3inv(matrix1,matrix2)
   do ii=1,3
     axes(:,ii)=rprimd(:,1)*matrix2(ii,1)+rprimd(:,2)*matrix2(ii,2)+rprimd(:,3)*matrix2(ii,3)
   end do
   length_axis1=sqrt(axes(1,1)**2+axes(2,1)**2+axes(3,1)**2)
   length_axis2=sqrt(axes(1,2)**2+axes(2,2)**2+axes(3,2)**2)
   length_axis3=sqrt(axes(1,3)**2+axes(2,3)**2+axes(3,3)**2)
!  DEBUG
!  write(std_out,*)' testkgrid: axes=',axes(:,:)
!  write(std_out,*)' length_axis=',length_axis1,length_axis2,length_axis3
!  ENDDEBUG
 end if

!This routine examine only primitive k lattices.
 nshiftk=1

!If prtkpt/=0, will examine more grids than strictly needed
 buffer_scale=one
 if(prtkpt/=0)buffer_scale=two

 if(prtkpt/=0)then
   write(msg,'(a,a,a,a,a,a,a,a)' )ch10,&
     ' testkgrid : will perform the analysis of a series of k-grids.',ch10,&
     '  Note that kptopt=1 in this analysis, irrespective of its input value.',ch10,ch10,&
     ' Grid#    kptrlatt         shiftk         kptrlen       nkpt  iset',ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(iout,msg,'COLL')
   ABI_MALLOC(grid_list,(mkpt_list))
   ABI_MALLOC(kptrlen_list,(mkpt_list))
   grid_list(:)=0
   kptrlen_list(:)=0.0_dp
 end if

 if(sum(vacuum(:))==3)then

   kptrlatt(:,:)=0
   kptrlatt(1,1)=1
   kptrlatt(2,2)=1
   kptrlatt(3,3)=1
   shiftk(:,1)=0.0_dp
   kptrlen=1000.0_dp
   nkpt_current=1
   igrid_current=1

   if(prtkpt/=0)then
     write(msg,&
&     '(a,3i4,a,es14.4,a,es14.4,i8,i6,a,a,3i4,a,es14.4,a,a,3i4,a,es14.4,a)' )&
&     '    1  ',kptrlatt(:,1),'  ',shiftk(1,1),'  ',kptrlen,1,1,ch10,&
&     '       ',kptrlatt(:,2),'  ',shiftk(2,1),ch10,&
&     '       ',kptrlatt(:,3),'  ',shiftk(3,1),ch10
     call wrtout(std_out,msg,'COLL')
     call wrtout(iout,msg,'COLL')
!    The unit cell volume is fake
     ucvol=kptrlen**3
   end if

 else

   nkpt=0 ; nkpt_current=0 ; iscf=1 ; iset=1
   kptrlen_current=0.0_dp
   mult1=0 ; mult2=0 ; mult3=0 ; init_mult=1
   ABI_MALLOC(kpt,(3,nkpt))
   ABI_MALLOC(wtk,(nkpt))
   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!  Loop on different grids, the upper limit is only to avoid an infinite loop
   do igrid=1,1000

     kptrlatt_trial(:,:)=0
     kptrlatt_trial(1,1)=1
     kptrlatt_trial(2,2)=1
     kptrlatt_trial(3,3)=1
     shiftk_trial(:,1)=0.0_dp

!    1-dimensional system
     if(sum(vacuum(:))==2)then
       if(vacuum(1)==0)then
         kptrlatt_trial(1,1)=2*igrid ; shiftk_trial(1,1)=0.5_dp
       else if(vacuum(2)==0)then
         kptrlatt_trial(2,2)=2*igrid ; shiftk_trial(2,1)=0.5_dp
       else if(vacuum(3)==0)then
         kptrlatt_trial(3,3)=2*igrid ; shiftk_trial(3,1)=0.5_dp
       end if
     end if

!    2-dimensional system
     if(sum(vacuum(:))==1)then

!      Treat hexagonal holohedries separately
       if(iholohedry==3)then

!        write(std_out,*)' testkgrid: 2D, hexagonal'

         mult1=mult1+1
         nset=4
         if(iset==1)then
           rsuper(:,1)=axes(:,1)*mult1
           rsuper(:,2)=axes(:,2)*mult1
           shiftk_trial(:,1)=0.0_dp
         else if(iset==2)then
           rsuper(:,1)=(axes(:,1)-axes(:,2))  *mult1
           rsuper(:,2)=(axes(:,1)+2*axes(:,2))*mult1
           shiftk_trial(1,1)=1.0_dp/3.0_dp
           shiftk_trial(2,1)=1.0_dp/3.0_dp
         else if(iset==3)then
           rsuper(:,1)=(axes(:,1)-axes(:,2))  *mult1
           rsuper(:,2)=(axes(:,1)+2*axes(:,2))*mult1
           shiftk_trial(:,1)=0.0_dp
         else if(iset==4)then
           rsuper(:,1)=axes(:,1)*mult1
           rsuper(:,2)=axes(:,2)*mult1
           shiftk_trial(1,1)=0.5_dp
           shiftk_trial(2,1)=0.5_dp
         end if

       else
!        Now treat all other holohedries
         length1=length_axis1*mult1
         length2=length_axis2*mult2
!        write(std_out,*)' testkgrid: (2d) length=',length1,length2
         if(abs(length1-length2)<tol8)then
           mult1=mult1+1
           mult2=mult2+1
         else if(length1>length2)then
           mult2=mult2+1
         else if(length2>length1)then
           mult1=mult1+1
         end if
         nset=4
!        iset==5 and 6 are allowed only for centered lattice
         if(center==1)nset=6
         if(iset==1 .or. iset==2)then
           rsuper(:,1)=axes(:,1)*mult1
           rsuper(:,2)=axes(:,2)*mult2
         else if(iset==3 .or. iset==4)then
           rsuper(:,1)=axes(:,1)*mult1-axes(:,2)*mult2
           rsuper(:,2)=axes(:,1)*mult1+axes(:,2)*mult2
         else if(iset==5 .or. iset==6)then
           rsuper(:,1)=axes(:,1)*(mult1-0.5_dp)-axes(:,2)*(mult2-0.5_dp)
           rsuper(:,2)=axes(:,1)*(mult1-0.5_dp)+axes(:,2)*(mult2-0.5_dp)
         end if
!        This was the easiest way to code all even mult1 and mult2 pairs:
!        make separate series for this possibility.
         if(iset==2 .or. iset==4 .or. iset==6)then
           rsuper(:,1)=2.0_dp*rsuper(:,1)
           rsuper(:,2)=2.0_dp*rsuper(:,2)
         end if
         shiftk_trial(1,1)=0.5_dp
         shiftk_trial(2,1)=0.5_dp

       end if

!      Put back the inactive direction
       if(dirvacuum==1)then
         rsuper(:,3)=rsuper(:,1)
         shiftk_trial(3,1)=shiftk_trial(1,1)
         rsuper(:,1)=rprimd(:,1)
         shiftk_trial(1,1)=0.0_dp
       else if(dirvacuum==2)then
         rsuper(:,3)=rsuper(:,1)
         shiftk_trial(3,1)=shiftk_trial(1,1)
         rsuper(:,1)=rsuper(:,2)
         shiftk_trial(1,1)=shiftk_trial(2,1)
         rsuper(:,2)=rprimd(:,2)
         shiftk_trial(2,1)=0.0_dp
       else if(dirvacuum==3)then
         rsuper(:,3)=rprimd(:,3)
         shiftk_trial(3,1)=0.0_dp
       end if

!      The supercell and the corresponding shift have been generated !
!      Convert cartesian coordinates into kptrlatt_trial
       do ii=1,3
         kptrlatt_trial(:,ii)=nint( gprimd(1,:)*rsuper(1,ii)+&
&         gprimd(2,:)*rsuper(2,ii)+&
&         gprimd(3,:)*rsuper(3,ii)  )
       end do

!      End of 2-dimensional system
     end if

!    3-dimensional system
     if(sum(vacuum(:))==0)then
!      Treat hexagonal holohedries separately
       if(iholohedry==6)then
         length1=length_axis1*mult1
         length3=length_axis3*mult3
!        write(std_out,*)' testkgrid: (hex) lengths=',length1,length2
         if(abs(length1-length3)<tol8)then
           mult1=mult1+1
           mult3=mult3+1
         else if(length1>length3)then
           mult3=mult3+1
         else if(length3>length1)then
           mult1=mult1+1
         end if
         nset=4
         if(iset==1)then
           rsuper(:,1)=axes(:,1)*mult1
           rsuper(:,2)=axes(:,2)*mult1
           rsuper(:,3)=axes(:,3)*mult3
           shiftk_trial(:,1)=0.0_dp
           shiftk_trial(3,1)=0.5_dp
         else if(iset==2)then
           rsuper(:,1)=(axes(:,1)-axes(:,2))  *mult1
           rsuper(:,2)=(axes(:,1)+2*axes(:,2))*mult1
           rsuper(:,3)=axes(:,3)*mult3
           shiftk_trial(1,1)=1.0_dp/3.0_dp
           shiftk_trial(2,1)=1.0_dp/3.0_dp
           shiftk_trial(3,1)=0.5_dp
         else if(iset==3)then
           rsuper(:,1)=(axes(:,1)-axes(:,2))  *mult1
           rsuper(:,2)=(axes(:,1)+2*axes(:,2))*mult1
           rsuper(:,3)=axes(:,3)*mult3
           shiftk_trial(:,1)=0.0_dp
           shiftk_trial(3,1)=0.5_dp
         else if(iset==4)then
           rsuper(:,1)=axes(:,1)*mult1
           rsuper(:,2)=axes(:,2)*mult1
           rsuper(:,3)=axes(:,3)*mult3
           shiftk_trial(:,1)=0.5_dp
         end if

       else
!        Now treat all other holohedries
         length1=length_axis1*mult1
         length2=length_axis2*mult2
         length3=length_axis3*mult3
!        write(std_out,*)' testkgrid: length=',length1,length2,length3
         if(length2>length1+tol8 .and. length3>length1+tol8)then
           mult1=mult1+1
         else if(length1>length2+tol8 .and. length3>length2+tol8)then
           mult2=mult2+1
         else if(length1>length3+tol8 .and. length2>length3+tol8)then
           mult3=mult3+1
         else if(abs(length2-length3)<tol8 .and. &
&           abs(length1-length3)<tol8 .and. &
&           abs(length1-length2)<tol8        )then
           mult1=mult1+1 ; mult2=mult2+1 ; mult3=mult3+1
         else if(abs(length1-length2)<tol8)then
           mult1=mult1+1 ; mult2=mult2+1
         else if(abs(length1-length3)<tol8)then
           mult1=mult1+1 ; mult3=mult3+1
         else if(abs(length2-length3)<tol8)then
           mult2=mult2+1 ; mult3=mult3+1
         end if
         nset=6
         if(center==-1 .or. center==-3)nset=8
         if(iset==1 .or. iset==2)then
!          Simple lattice of k points
           rsuper(:,1)=axes(:,1)*mult1
           rsuper(:,2)=axes(:,2)*mult2
           rsuper(:,3)=axes(:,3)*mult3
           shiftk_trial(:,1)=0.5_dp
         else if(iset==3 .or. iset==4)then
!          FCC lattice of k points = BCC lattice in real space
           rsuper(:,1)=-axes(:,1)*mult1+axes(:,2)*mult2+axes(:,3)*mult3
           rsuper(:,2)= axes(:,1)*mult1-axes(:,2)*mult2+axes(:,3)*mult3
           rsuper(:,3)= axes(:,1)*mult1+axes(:,2)*mult2-axes(:,3)*mult3
           shiftk_trial(:,1)=0.5_dp
         else if(iset==5 .or. iset==6)then
!          BCC lattice of k points = FCC lattice in real space
           rsuper(:,1)=                 axes(:,2)*mult2+axes(:,3)*mult3
           rsuper(:,2)= axes(:,1)*mult1                +axes(:,3)*mult3
           rsuper(:,3)= axes(:,1)*mult1+axes(:,2)*mult2
!          The BCC lattice has no empty site with full symmetry
           shiftk_trial(:,1)=0.0_dp
         else if(iset==7 .or. iset==8)then
!          iset==7 and 8 are allowed only for centered lattice
           mult1h=mult1-0.5_dp
           mult2h=mult2-0.5_dp
           mult3h=mult3-0.5_dp
           if(center==-1)then
!            FCC lattice of k points = BCC lattice in real space
             rsuper(:,1)=-axes(:,1)*mult1h+axes(:,2)*mult2h+axes(:,3)*mult3h
             rsuper(:,2)= axes(:,1)*mult1h-axes(:,2)*mult2h+axes(:,3)*mult3h
             rsuper(:,3)= axes(:,1)*mult1h+axes(:,2)*mult2h-axes(:,3)*mult3h
             shiftk_trial(:,1)=0.5_dp
           else if(center==-3)then
!            BCC lattice of k points = FCC lattice in real space
             rsuper(:,1)=                  axes(:,2)*mult2h+axes(:,3)*mult3h
             rsuper(:,2)= axes(:,1)*mult1h                 +axes(:,3)*mult3h
             rsuper(:,3)= axes(:,1)*mult1h+axes(:,2)*mult2h
!            The BCC lattice has no empty site with full symmetry
             shiftk_trial(:,1)=0.0_dp
           end if
         end if
!        This was the easiest way to code all even mult1, mult2, mult3 triplets:
!        make separate series for this possibility.
         if(2*(iset/2)==iset)then
           rsuper(:,1)=2.0_dp*rsuper(:,1)
           rsuper(:,2)=2.0_dp*rsuper(:,2)
           rsuper(:,3)=2.0_dp*rsuper(:,3)
         end if
       end if

!      write(std_out,*)' testkgrid: gprimd=',gprimd(:,:)
!      write(std_out,*)' testkgrid: rsuper=',rsuper(:,:)
!      write(std_out,*)' testkgrid: iset  =',iset

!      The supercell and the corresponding shift have been generated!
!      Convert cartesian coordinates into kptrlatt_trial
       do ii=1,3
         kptrlatt_trial(:,ii)=nint( gprimd(1,:)*rsuper(1,ii)+&
&         gprimd(2,:)*rsuper(2,ii)+&
&         gprimd(3,:)*rsuper(3,ii)  )
       end do

!      End of 3-dimensional system
     end if

!    write(std_out,*)' testkgrid: before getkgrid'
!    write(std_out,*)' testkgrid: rprimd=',rprimd(:,:)
!    write(std_out,*)' testkgrid: kptrlatt_trial=',kptrlatt_trial(:,:)

     call getkgrid(0,0,iscf,kpt,&
&     kptopt,kptrlatt_trial,kptrlen_trial,&
&     msym,nkpt,nkpt_trial,nshiftk,nsym,rprimd,&
&     shiftk_trial,symafm,symrel,vacuum,wtk)

!    write(std_out,*)' testkgrid: after getkgrid'

!    In case one does not need the full list of grids, will take a shortcut, and go to one of the last grids of the series,
!    that generates a kptrlen_trial that is just below kptrlen.
     if(prtkpt==0 .and. init_mult==1 .and. kptrlen_trial<(half-tol8)*kptrlen )then
       iscale=int((one-tol8)*kptrlen/kptrlen_trial)
       mult1=mult1*iscale
       mult2=mult2*iscale
       mult3=mult3*iscale
       init_mult=0
!       write(std_out,*)' testkgrid: iscale=',iscale
       kptrlatt_trial(:,:)=kptrlatt_trial(:,:)*iscale
       call getkgrid(0,0,iscf,kpt,&
&       kptopt,kptrlatt_trial,kptrlen_trial,&
&       msym,nkpt,nkpt_trial,nshiftk,nsym,rprimd,&
&       shiftk_trial,symafm,symrel,vacuum,wtk)
     end if

     if( (kptrlen_trial+tol8>kptrlen*(1.0_dp+tol8) .and. nkpt_current==0) .or. &
&     (kptrlen_trial+tol8>kptrlen*(1.0_dp+tol8) .and. nkpt_trial<nkpt_current) .or. &
&     (nkpt_trial==nkpt_current  .and. kptrlen_trial>kptrlen_current*(1.0_dp+tol8)))then

       kptrlatt_current(:,:)=kptrlatt_trial(:,:)
       nkpt_current=nkpt_trial
       shiftk_current(:,:)=shiftk_trial(:,:)
       kptrlen_current=kptrlen_trial
       igrid_current=igrid
     end if

     if(prtkpt/=0)then
       write(msg,'(i5,a,3i4,a,es14.4,a,es14.4,i8,i6,a,a,3i4,a,es14.4,a,a,3i4,a,es14.4,a)' )&
&       igrid,'  ',kptrlatt_trial(:,1),'  ',shiftk_trial(1,1),&
&       '  ',kptrlen_trial,nkpt_trial,iset,ch10,&
&       '       ',kptrlatt_trial(:,2),'  ',shiftk_trial(2,1),ch10,&
&       '       ',kptrlatt_trial(:,3),'  ',shiftk_trial(3,1),ch10
       call wrtout(std_out,msg,'COLL')
       call wrtout(iout,msg,'COLL')

!      Keep track of this grid, if it is worth
       if(kptrlen_trial > kptrlen_list(nkpt_trial)*(1.0_dp+tol8))then
         grid_list(nkpt_trial)=igrid
         kptrlen_list(nkpt_trial)=kptrlen_trial
       end if
     end if

!    Treat 1-D case
     if( sum(vacuum(:))==2 .and. kptrlen_trial>buffer_scale*(1.0_dp+tol8)*kptrlen )exit

!    Treat 2-D case or 3-D case
     if( sum(vacuum(:))<=1 .and. kptrlen_trial>buffer_scale*(1.0_dp+tol8)*kptrlen )then
!      The present set of sets of k points is finished:
!      either it was the last, or one has to go to the next one
       if(iset==nset)exit
       iset=iset+1
       mult1=0 ; mult2=0 ; mult3=0 ; init_mult=1
     end if

   end do ! igrid=1,1000

   ABI_FREE(kpt)
   ABI_FREE(wtk)

   kptrlatt(:,:)=kptrlatt_current(:,:)
   shiftk(:,:)=shiftk_current(:,:)
   kptrlen=kptrlen_current

 end if ! test on the number of dimensions

 if(prtkpt/=0)then

!  sqrt(1/2) comes from the FCC packing, the best one
   factor=sqrt(0.5_dp)/ucvol/dble(nsym)
   ndims=3
   if(sum(vacuum(:))/=0)then
     if(sum(vacuum(:))==1)then
!      sqrt(3/4) comes from the hex packing, the best one
!      one multiplies by 2 because nsym is likely twice the number
!      of symmetries that can be effectively used in 2D
       ndims=2 ; factor=sqrt(0.75_dp)/surface/dble(nsym)*2
       write(msg,'(2a)' )ch10,' Note that the system is bi-dimensional.'
     else if(sum(vacuum(:))==2)then
       ndims=1 ; factor=1/ucvol
       write(msg,'(2a)' )ch10,' Note that the system is uni-dimensional.'
     else if(sum(vacuum(:))==3)then
       ndims=0
       write(msg,'(2a)' )ch10,' Note that the system is zero-dimensional.'
     end if
     call wrtout(std_out,msg,'COLL')
     call wrtout(iout,msg,'COLL')
   end if

!  The asymptotic value of the merit factor is determined
!  by the set of symmetries: in 3D, if it includes the
!  inversion symmetry, the limit will be 1, if not, it
!  will be two. In 2D, if it includes the inversion symmetry
!  and an operation that maps z on -z, it will tend to one,
!  while if only one of these operations is present,
!  it will tend to two, and if none is present, it will tend to four.
   write(msg,'(11a)' )ch10,&
&   ' List of best grids, ordered by nkpt.',ch10,&
&   '  (stop at a value of kptrlen 20% larger than the target value).',ch10,&
&   '  (the merit factor will tend to one or two in 3 dimensions)',ch10,&
&   '  (and to one, two or four in 2 dimensions)',ch10,ch10,&
&   '    nkpt   kptrlen    grid#  merit_factor'
   call wrtout(std_out,msg,'COLL')
   call wrtout(iout,msg,'COLL')

   kptrlen_max=0.0_dp
   do ii=1,mkpt_list
     if(kptrlen_list(ii)>kptrlen_max*(1.0_dp+tol8))then
       kptrlen_max=kptrlen_list(ii)
       merit_factor=kptrlen_max**ndims/dble(ii)*factor
       write(msg, '(i6,es14.4,i6,f12.4)' )ii,kptrlen_max,grid_list(ii),merit_factor
       call wrtout(std_out,msg,'COLL')
       call wrtout(iout,msg,'COLL')
     end if
     if(kptrlen_max>1.2_dp*(1.0_dp-tol8)*kptrlen_target)exit
   end do

   write(msg,'(a,a,es14.4,a,a,i6,a,a,a,es14.4,a,i6)' )ch10,&
&   ' For target kptrlen=',kptrlen_target,',',&
&   ' the selected grid is number',igrid_current,',',ch10,&
&   '     giving kptrlen=',kptrlen_current,' with nkpt=',nkpt_current
   call wrtout(std_out,msg,'COLL')
   call wrtout(iout,msg,'COLL')

   write(msg,'(a,a,a,a)' )ch10,&
&   ' testkgrid : stop after analysis of a series of k-grids.',ch10,&
&   '  For usual production runs, set prtkpt back to 0 (the default).'
   call wrtout(std_out,msg,'COLL',do_flush=.True.)
   call wrtout(iout,msg,'COLL',do_flush=.True.)

   call abi_abort('PERS',exit_status=0,print_config=.false.)
 end if

end subroutine testkgrid
!!***

!!****f* m_kpts/mknormpath
!! NAME
!! mknormpath
!!
!! FUNCTION
!! Please do not use this  routine, use make_normpath instead.
!! mknormpath should be removed
!!
!!  This simple routine generates a normalized path that can be used to plot a band
!!  structures in an easy way. For normalized path we mean a path where the number
!!  of division on each segment is proportional to the length of the segment itself.
!!  To generate the above mentioned path, the subroutine must be called twice.
!!  The first call reports the total number of divisions in the normalized path, dimension
!!  that is required to correctly allocate the array.
!!  The second call calculates the reduced coordinates of the circuit.
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2020 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! nbounds=number of points defining the path
!! ndiv_small=number of points to be used to sample the smallest
!!  segment defined by bounds(:,1:nbounds)
!! bounds(3,nbounds)=points defining the path
!! gmet(3,3)=metric
!!
!! OUTPUT
!! ndiv(nbounds-1)= number of divisions for each segment
!! npt_tot=total number of points sampled along the circuit
!! path(3,npt_tot)= normalized path in reciprocal space
!!
!! TODO
!!  Do not use this routine, it is obsolete and should be replaced by make_path in m_bz_mesh.
!!
!! PARENTS
!!      inkpts
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE


subroutine mknormpath(nbounds,bounds,gmet,ndiv_small,ndiv,npt_tot,path)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbounds,ndiv_small
 integer,intent(inout) :: npt_tot
!arrays
 integer,intent(inout) :: ndiv(nbounds-1)
 real(dp),intent(in) :: bounds(3,nbounds),gmet(3,3)
 real(dp),intent(out),optional :: path(3,npt_tot)

!Local variables-------------------------------
!scalars
 integer :: idx,ii,jp
 real(dp) :: fct
 character(len=500) :: msg
!arrays
 real(dp) :: dd(3),lng(nbounds-1)

! *************************************************************************

 if (ndiv_small<=0) then
   write(msg,'(3a,i0)')&
&   'The argument ndiv_small should be a positive number,',ch10,&
&   'however, ndiv_small=',ndiv_small
   MSG_ERROR(msg)
 end if

 do ii=1,nbounds-1
   dd(:)=bounds(:,ii+1)-bounds(:,ii)
   lng(ii)= sqrt( dd(1)*gmet(1,1)*dd(1)+ &
&   dd(2)*gmet(2,2)*dd(2)+ &
&   dd(3)*gmet(3,3)*dd(3)+ &
&   2.0d0*(dd(1)*gmet(1,2)*dd(2)+ &
&   dd(1)*gmet(1,3)*dd(3)+ &
&   dd(2)*gmet(2,3)*dd(3)) &
&   )
 end do
 write(std_out,*)lng
 fct=minval(lng)

!Avoid division by zero if k(:,i+1)=k(:,i)
 if (abs(fct)<tol6) then
   write(msg,'(3a)')&
&   'found two consecutive points in the path which are equal',ch10,&
&   'This is not allowed, please modify the path in your input file'
   MSG_ERROR(msg)
 end if

 fct=fct/ndiv_small
 ndiv(:)=nint(lng(:)/fct)
!The 1 stand for the first point
 npt_tot=sum(ndiv)+1

 if (.not.present(path)) then
   write(msg,'(2a,i8)')ch10,' mknormpath : total number of points on the path: ',npt_tot
   call wrtout(std_out,msg,'COLL')
   write(msg,'(2a)')ch10,' Number of divisions for each segment of the normalized path: '
   call wrtout(std_out,msg,'COLL')
   do ii=1,nbounds-1
     write(msg,'(2(3f8.5,a),i5,a)')&
     bounds(:,ii),' ==> ',bounds(:,ii+1),' ( ndiv: ',ndiv(ii),' )'
     call wrtout(std_out,msg,'COLL')
   end do
   write(msg,'(a)')ch10
   call wrtout(std_out,msg,'COLL')
 else
   write(msg,'(2a)')ch10,' Normalized Path: '
   call wrtout(std_out,msg,'COLL')
   idx=1
   do ii=1,nbounds-1
     do jp=1,ndiv(ii)
       path(:,idx)=bounds(:,ii)+(jp-1)*(path(:,ii+1)-path(:,ii))/ndiv(ii)
       write(msg,'(i4,4x,3(f8.5,1x))')idx,path(:,idx)
       call wrtout(std_out,msg,'COLL')
       idx=idx+1
     end do
   end do
 end if

end subroutine mknormpath
!!***

end module m_kpts
!!***
