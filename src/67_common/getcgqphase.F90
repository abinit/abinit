!{\src2tex{textfont=tt}}
!!****f* ABINIT/getcgqphase
!! NAME
!! getcgqphase
!!
!! FUNCTION
!! extract phases from wave functions, to cancel contributions to gkk matrix elements
!!
!! COPYRIGHT
!! Copyright (C) 2011-2016 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  timrev = flag for use of time reversal symmetry
!!  cg = input wavefunctions
!!  mcg = dimension of cg = nspinor*mband*mpw*mkmem
!!  cgq = input wavefunctions at k+q
!!  mcgq = dimension of cgq = nspinor*mband*mpw*mkmem
!!  mpi_enreg = datastructure for mpi communication
!!  nkpt_rbz = number of k-points in reduced zone for present q point
!!  npwarr = array of numbers of plane waves for each k-point
!!  npwar1 = array of numbers of plane waves for each k+q point 
!!
!! OUTPUT
!!  phasecg = phase of different wavefunction products <k,n | k+q,n'>
!!
!! PARENTS
!!      dfpt_looppert
!!
!! CHILDREN
!!      smatrix,wrtout,xmpi_barrier,xmpi_sum_master
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"



subroutine getcgqphase(dtset, timrev, cg,  mcg,  cgq, mcgq, mpi_enreg, &
&    nkpt_rbz, npwarr, npwar1, phasecg)

 use m_profiling_abi
 use defs_basis
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getcgqphase'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments -------------------------------
 ! scalars 
 integer, intent(in) :: mcg, mcgq, timrev
 integer, intent(in) :: nkpt_rbz
 type(dataset_type), intent(in) :: dtset

 ! arrays
 integer, intent(in) :: npwarr(nkpt_rbz)
 integer, intent(in) :: npwar1(nkpt_rbz)
 real(dp), intent(in) :: cg(2,mcg)
 real(dp), intent(in) :: cgq(2,mcgq)
 type(MPI_type), intent(in) :: mpi_enreg
 real(dp),intent(out) :: phasecg(2, dtset%mband*dtset%mband*nkpt_rbz&
&     *dtset%nsppol)

!Local variables -------------------------
 ! local vars
 integer :: icg, icgq, isppol, ikpt, ipw
 integer :: istate, iband1, iband2
 integer :: npw_k, npw_q
 integer :: me, ierr, master, spaceComm, nprocs
 integer :: usepaw
 integer :: ddkflag, itrs, job, maxbd, mcg1_k, minbd, shiftbd
 real(dp) :: normsmat

 integer, allocatable :: sflag_k(:)
 integer, allocatable :: pwind_k(:)

 real(dp) :: cg1_dummy(1,1)
 real(dp) :: smat_inv_dummy(1,1,1)
 real(dp) :: smat_k_paw_dummy(1,1,1)
 real(dp) :: dtm_k_dummy(2)
 real(dp), allocatable :: smat_k(:,:,:)
 real(dp), allocatable :: pwnsfac_k(:,:)
 logical, allocatable :: my_kpt(:,:)
 character(len=500) :: message

! *********************************************************************

 ABI_ALLOCATE(smat_k,(2,dtset%mband,dtset%mband))
 ABI_ALLOCATE(sflag_k,(dtset%mband))

!dummy use of timrev so abirules stops complaining.
 icg = timrev

!!MPI data for future use
 spaceComm=mpi_enreg%comm_cell
 nprocs=xmpi_comm_size(spaceComm)
 master=0
 me=mpi_enreg%me_kpt

 ABI_ALLOCATE(my_kpt, (nkpt_rbz, dtset%nsppol))
 my_kpt = .true.
 if (mpi_enreg%nproc_kpt > 1) then
   do isppol = 1, dtset%nsppol
     do ikpt = 1, nkpt_rbz
       my_kpt(ikpt, isppol) = .not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,&
&       dtset%nband(ikpt),isppol,me))
     end do
   end do
 end if


!make trivial association of G vectors: we just want <psi_k| psi_k+q>
!TODO: check this is correct wrt arrangement of kg vectors for k+q
!looks ok : usually made in initberry, from the translations associated
!to the symops, scalar product with the G vectors. The symop is the one
!used to go from the irreducible k to the full zone k. In present context
!we should be using only the reduced zone, and anyhow have the same k-grid
!for the gkk matrix elements and for the cg here...
 ABI_ALLOCATE(pwind_k,(dtset%mpw))
 ABI_ALLOCATE(pwnsfac_k,(4,dtset%mpw))
 do ipw = 1, dtset%mpw
   pwind_k(ipw) = ipw
   pwnsfac_k(1,ipw) = one
   pwnsfac_k(2,ipw) = zero
   pwnsfac_k(3,ipw) = one
   pwnsfac_k(4,ipw) = zero
 end do

!flags for call to smatrix
 usepaw = 0 ! for now
 ddkflag = 0
 itrs = 0
 job = 0
 maxbd = 1
 mcg1_k = 1
 minbd = 1
 shiftbd = 1

!from overlap matrix for each wavefunction, extract phase
 icg = 0
 icgq = 0
 istate = 0

 phasecg = zero
 do isppol = 1, dtset%nsppol
   do ikpt = 1, nkpt_rbz
!    each proc only has certain k
     if (.not. my_kpt(ikpt, isppol)) then
       istate = istate +  dtset%nband(ikpt)*dtset%nband(ikpt)
       cycle
     end if

     npw_k = npwarr(ikpt)
     npw_q= npwar1(ikpt)

!    TODO: question: are the k-points in the ibz correctly ordered in cg and cgq? if not the icg below have to be adapted.
     sflag_k = 0 ! make sure all elements are calculated
     smat_k = zero

     call smatrix(cg, cgq, cg1_dummy, ddkflag, dtm_k_dummy, icg, icgq,&
&     itrs, job, maxbd, mcg, mcgq, mcg1_k, minbd,dtset%mpw, dtset%mband, dtset%mband,&
&     npw_k, npw_q, dtset%nspinor, pwind_k, pwnsfac_k, sflag_k, shiftbd,&
&     smat_inv_dummy, smat_k, smat_k_paw_dummy, usepaw)

     icg  = icg  + npw_k*dtset%nspinor*dtset%nband(ikpt)
     icgq = icgq + npw_q*dtset%nspinor*dtset%nband(ikpt)

     do iband1 = 1, dtset%nband(ikpt)
       do iband2 = 1, dtset%nband(ikpt)
         istate = istate + 1
!        normalise the overlap matrix element to get just the phase difference phi_k - phi_k+q
         normsmat = sqrt(smat_k(1,iband2, iband1)**2 &
&         + smat_k(2,iband2, iband1)**2)
         if (normsmat > tol12) then
           phasecg(:,istate) = smat_k(:,iband2, iband1) / normsmat
!          NOTE: 21/9/2011 these appear to be always 1, i, or -i, to within 1.e-5 at worst!
         end if
       end do
     end do
   end do
 end do

!eventually do an mpi allreduce over the k-points for phasecg
 if (nprocs>1) then
   call xmpi_barrier(spaceComm)
   call xmpi_sum_master(phasecg,master,spaceComm,ierr)
   call xmpi_barrier(spaceComm)
   if (1==1) then
     write(message,'(a)') '  In getcgqphase - contributions to phasecg collected'
     call wrtout(std_out,message,'PERS')
   end if
 end if

 ABI_DEALLOCATE(sflag_k)
 ABI_DEALLOCATE(smat_k)
 ABI_DEALLOCATE(pwind_k)
 ABI_DEALLOCATE(pwnsfac_k)
 ABI_DEALLOCATE(my_kpt)

end subroutine getcgqphase
!!***
