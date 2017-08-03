!{\src2tex{textfont=tt}}
!!****f* ABINIT/pead_nl_mv
!! NAME
!! pead_nl_mv
!!
!! FUNCTION
!! Compute the finite difference expression of the k-point derivative
!! using the PEAD formulation of the third-order energy
!! (see Nunes and Gonze PRB 63, 155107 (2001) Eq. 102)
!! and the finite difference formula of Marzari and Vanderbilt
!! (see Marzari and Vanderbilt, PRB 56, 12847 (1997), Appendix B)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol) = array for planewave
!!                                         coefficients of wavefunctions
!!  cgindex(nkpt2,nsppol) = for each k-point, cgindex stores the location
!!                          of the WF in the cg array
!!  cg1 = first-order wavefunction relative to the perturbations i1pert
!!  cg3 = first-order wavefunction relative to the perturbations i3pert
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  i1pert,i3pert = type of perturbation that has to be computed
!!  i1dir,i3dir=directions of the corresponding perturbations
!!  kneigh(30,nkpt2) = index of the neighbours of each k-point
!!  kg_neigh(30,nkpt2,3) = necessary to construct the vector joining a k-point 
!!                         to its nearest neighbour in case of a single k-point, 
!!                         a line of k-points or a plane of k-points. 
!!                         See getshell.F90 for details
!!  kptindex(2,nkpt3)= index of the k-points in the reduced BZ
!!                     related to a k-point in the full BZ
!!  kpt3(3,nkpt3) = reduced coordinates of k-points in the full BZ
!!  mband = maximum number of bands
!!  mkmem = maximum number of k points which can fit in core memory
!!  mkmem_max = maximal number of k-points on each processor (MPI //)
!!  mk1mem = maximum number of k points for first-order WF
!!           which can fit in core memory
!!  mpi_enreg=MPI-parallelisation information
!!  mpw   = maximum number of planewaves in basis sphere (large number)
!!  mvwtk(30,nkpt) = weights to compute the finite difference ddk
!!  natom = number of atoms in unit cell
!!  nkpt2 = number of k-points in the reduced part of the BZ
!!          nkpt2 = nkpt/2 in case of time-reversal symmetry (kptopt = 2)
!!  nkpt3 = number of k-points in the full BZ
!!  nneigh = total number of neighbours required to evaluate the finite
!!           difference formula
!!  npwarr(nkpt) = array holding npw for each k point
!!  nspinor = number of spinorial components of the wavefunctions
!!  nsppol = number of channels for spin-polarization (1 or 2)
!!  pwind(mpw,nneigh,mkmem) = array used to compute the overlap matrix smat
!!                           between k-points
!!
!! OUTPUT
!!  d3_berry(2,3) = Berry-phase part of the third-order energy
!!
!! SIDE EFFECTS
!!  mpi_enreg=MPI-parallelisation information
!!
!! NOTES
!! For a given set of values of i1pert,i3pert,i1dir and
!! i3dir, the routine computes the k-point derivatives for
!! 12dir = 1,2,3
!!
!! TODO
!!
!! PARENTS
!!      pead_nl_loop
!!
!! CHILDREN
!!      dzgedi,dzgefa,mpi_recv,mpi_send,status,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine pead_nl_mv(cg,cgindex,cg1,cg3,dtset,dtfil,d3_berry,gmet,&
&                   i1pert,i3pert,i1dir,i3dir,kneigh,kg_neigh,kptindex,&
&                   kpt3,mband,mkmem,mkmem_max,mk1mem,mpi_enreg,&
&                   mpw,mvwtk,natom,nkpt2,nkpt3,nneigh,npwarr,nspinor,&
&                   nsppol,pwind)

 use m_profiling_abi

 use defs_basis
 use defs_abitypes
 use m_xmpi

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pead_nl_mv'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!
!---  Arguments : integer scalars
 integer, intent(in) :: i1dir,i1pert,i3dir,i3pert,mband,mk1mem
 integer, intent(in) :: mkmem,mkmem_max,mpw,natom
 integer, intent(in) :: nkpt2,nkpt3,nneigh,nspinor,nsppol
!
!---  Arguments : integer arrays
 integer, intent(in) :: cgindex(nkpt2,nsppol)
 integer, intent(in) :: kneigh(30,nkpt2),kg_neigh(30,nkpt2,3),kptindex(2,nkpt3)
 integer, intent(in) :: npwarr(nkpt2),pwind(mpw,nneigh,mkmem)
!
!---  Arguments : real(dp) scalars
!
!---  Arguments : real(dp) arrays
 real(dp), intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp), intent(in) :: cg1(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp), intent(in) :: cg3(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp), intent(in) :: gmet(3,3),kpt3(3,nkpt3)
 real(dp), intent(in) :: mvwtk(30,nkpt2)
 real(dp), intent(out) :: d3_berry(2,3)
!
!---  Arguments : structured datatypes
 type(MPI_type), intent(in) :: mpi_enreg
 type(datafiles_type), intent(in) :: dtfil
 type(dataset_type), intent(in) :: dtset

!Local variables-------------------------------
!
!---- Local variables : integer scalars
 integer :: count,counter,count1,iband,icg
 integer :: ierr,iexit,ii,ikpt,ikpt_loc,ikpt2
 integer :: ikpt_rbz,ineigh,info,ipw,isppol,jband,jcg,jj,jkpt,job,jpw, jkpt2, jkpt_rbz
 integer :: lband,lpband,nband_occ,npw_k,npw_k1,my_source,his_source,dest,tag
 integer :: spaceComm
 integer,parameter :: level=52
 integer :: bdtot_index
!
!---- Local variables : integer arrays
 integer,allocatable :: ipvt(:)
 integer, allocatable :: bd_index(:,:)
!
!---- Local variables : real(dp) scalars
 real(dp) :: dotnegi,dotnegr,dotposi,dotposr
! real(dp) :: c1,c2 ! appear commented out below
!
!---- Local variables : real(dp) arrays
 real(dp) :: d3_aux(2,3),det(2,2),dk(3),dk_(3)
 real(dp) :: z1(2),z2(2)
 real(dp),allocatable :: buffer(:,:),cgq(:,:),cg1q(:,:),cg3q(:,:)
 real(dp),allocatable :: qmat(:,:,:),s13mat(:,:,:),s1mat(:,:,:),s3mat(:,:,:)
 real(dp),allocatable :: smat(:,:,:),zgwork(:,:)
!
!---- Local variables : character variables
 character(len=500) :: message
!
!---- Local variables : structured datatypes

#if defined HAVE_MPI
             integer :: status1(MPI_STATUS_SIZE)
!BEGIN TF_CHANGES
             spaceComm=mpi_enreg%comm_cell
!END TF_CHANGES
#endif


! ***********************************************************************

!DEBUG
!write(std_out,*)' pead_nl_mv : enter '
!flush(6)
!stop
!ENDDEBUG


 call status(0,dtfil%filstat,iexit,level,'enter         ')

 write(message,'(8a)') ch10,&
& ' pead_nl_mv : finite difference expression of the k-point derivative',ch10,&
& '           is performed using the PEAD formulation of ',&
& 'the third-order energy',ch10,&
& '           (see Nunes and Gonze PRB 63, 155107 (2001) Eq. 102)',ch10
!call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')


!fab: I think that the following restriction must be eliminated: 
!isppol = 1

 ikpt_loc = 0
 d3_aux(:,:) = 0_dp

 ABI_ALLOCATE(s13mat,(2,mband,mband))
 ABI_ALLOCATE(smat,(2,mband,mband))
 ABI_ALLOCATE(s1mat,(2,mband,mband))
 ABI_ALLOCATE(qmat,(2,mband,mband))
 ABI_ALLOCATE(ipvt,(mband))
 ABI_ALLOCATE(s3mat,(2,mband,mband))
 ABI_ALLOCATE(zgwork,(2,mband))
 ABI_ALLOCATE(bd_index, (nkpt2, nsppol))

 bdtot_index = 0
 do isppol = 1, nsppol
   do ikpt_rbz = 1, nkpt2
     bd_index(ikpt_rbz,isppol) = bdtot_index
     bdtot_index = bdtot_index + dtset%nband(ikpt_rbz+nkpt2*(isppol-1))
   end do
 end do

!fab: I think here I have to add the loop over spin
 
 do isppol = 1, nsppol

!  Loop over k-points
!  COMMENT: Every processor has to make mkmem_max iterations
!  even if mkmem < mkemem_max. This is due to the fact
!  that it still has to communicate its wavefunctions
!  to other processors even if it has no more overlap
!  matrices to compute.

   ikpt_loc = 0 ; ikpt = 0

   do while (ikpt_loc < mkmem_max)

     if (ikpt_loc < mkmem) ikpt = ikpt + 1

     if (xmpi_paral == 1) then
!      if ((minval(abs(mpi_enreg%proc_distrb(ikpt,1:mband,1:dtset%nsppol) &
!      &       - mpi_enreg%me)) /= 0).and.(ikpt_loc < mkmem)) cycle
       if(ikpt>nkpt2)then
         ikpt_loc=mkmem_max
         cycle
       end if
       if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,mband,-1,mpi_enreg%me)) then
         if(ikpt==nkpt2) ikpt_loc=mkmem_max
         cycle
       end if
     end if

     ikpt_loc = ikpt_loc + 1
     npw_k = npwarr(ikpt)
     counter = 100*ikpt
     call status(counter,dtfil%filstat,iexit,level,'loop over k   ')

     ii = cgindex(ikpt,isppol)

!    Loop on the  neighbours

     do ineigh = 1,nneigh

       s13mat(:,:,:) = zero
       smat(:,:,:) = zero
       s1mat(:,:,:) = zero
       s3mat(:,:,:) = zero
       qmat(:,:,:) = zero

       ikpt2  = kneigh(ineigh,ikpt)
       ikpt_rbz = kptindex(1,ikpt2)   ! index of the k-point in the reduced BZ
       jj = cgindex(ikpt_rbz,isppol)
       ! previous fixed value for nband_k now called nband_occ:
       !nband_occ = dtset%nband(ikpt_rbz+nkpt2*(isppol-1))
       ! TODO: check if all these bands are occupied in nsppol = 2 case
       nband_occ = 0
       do iband = 1, dtset%nband(ikpt_rbz+nkpt2*(isppol-1))
         if (dtset%occ_orig(bd_index(ikpt_rbz,isppol) + iband) > tol10) nband_occ = nband_occ + 1
       end do
       npw_k1 = npwarr(ikpt_rbz)
       dk_(:) = kpt3(:,ikpt2) - dtset%kptns(:,ikpt)
       dk(:)  = dk_(:) - nint(dk_(:)) + real(kg_neigh(ineigh,ikpt,:),dp)

       count = nspinor*mband*npw_k1
       ABI_ALLOCATE(cgq,(2,count))
       ABI_ALLOCATE(cg1q,(2,count))
       ABI_ALLOCATE(cg3q,(2,count))

#if defined HAVE_MPI

       my_source = mpi_enreg%proc_distrb(ikpt_rbz,1,1)

!      do dest = 0, mpi_enreg%nproc-1
       do dest = 0, maxval(mpi_enreg%proc_distrb(1:nkpt2,1:mband,1:dtset%nsppol))

         if ((dest==mpi_enreg%me).and.(ikpt_loc <= mkmem)) then
!          I am dest and have something to do

           if (my_source == mpi_enreg%me) then
!            I am destination and source
             jcg = cgindex(ikpt_rbz,isppol)

             cgq(:,1:count)  = cg(:,jcg+1:jcg+count)
             cg1q(:,1:count) = cg1(:,jcg+1:jcg+count)
             cg3q(:,1:count) = cg3(:,jcg+1:jcg+count)

           else
!            I am the destination but not the source -> receive

             tag = ikpt_rbz

             ABI_ALLOCATE(buffer,(2,3*count))

             call MPI_RECV(buffer,2*3*count,MPI_DOUBLE_PRECISION,my_source,tag,spaceComm,status1,ierr)

             cgq(:,1:count)  = buffer(:,1:count)
             cg1q(:,1:count) = buffer(:,count+1:2*count)
             cg3q(:,1:count) = buffer(:,2*count+1:3*count)
             ABI_DEALLOCATE(buffer)

           end if
           
         else if (ikpt_loc <= mpi_enreg%mkmem(dest)) then  ! dest != me and the dest has a k-point to treat

           jkpt=mpi_enreg%kpt_loc2ibz_sp(dest, ikpt_loc,1)
           jkpt2  = kneigh(ineigh,jkpt)
           jkpt_rbz = kptindex(1,jkpt2)   ! index of the k-point in the reduced BZ

           his_source = mpi_enreg%proc_distrb(jkpt_rbz,1,1)
           
           if (his_source == mpi_enreg%me) then

             jcg = cgindex(jkpt_rbz,isppol)

             tag = jkpt_rbz
             count1 = npwarr(jkpt_rbz)*mband*nspinor
             ABI_ALLOCATE(buffer,(2,3*count1))
             buffer(:,1:count1)            = cg(:,jcg+1:jcg+count1)
             buffer(:,count1+1:2*count1)   = cg1(:,jcg+1:jcg+count1)
             buffer(:,2*count1+1:3*count1) = cg3(:,jcg+1:jcg+count1)

             call MPI_SEND(buffer,2*3*count1,MPI_DOUBLE_PRECISION,dest,tag,spaceComm,ierr)

             ABI_DEALLOCATE(buffer)

           end if

         end if

       end do
!      
!      do jkpt = 1, nkpt2
!      
!      if ((jkpt == ikpt_rbz).and.(source /= mpi_enreg%me).and.&
!      &         (ikpt_loc <= mkmem)) then
!      
!      tag = jkpt
!      
!      allocate(buffer(2,3*count))
!      call MPI_RECV(buffer,2*3*count,MPI_DOUBLE_PRECISION,&
!      source,tag,spaceComm,status1,ierr)
!      
!      cgq(:,1:count)  = buffer(:,1:count)
!      cg1q(:,1:count) = buffer(:,count+1:2*count)
!      cg3q(:,1:count) = buffer(:,2*count+1:3*count)
!      deallocate(buffer)
!      
!      end if
!      
!      !        ----------------------------------------------------------------------------
!      !        --------------- Here: send the WF to all the cpus that need it -------------
!      !        ----------------------------------------------------------------------------
!      
!      do dest = 1, mpi_enreg%nproc
!      
!      if ((minval(abs(mpi_enreg%proc_distrb(jkpt,1:mband,1:dtset%nsppol) &
!      &           - mpi_enreg%me)) == 0).and.&
!      &           (mpi_enreg%kptdstrb(dest,ineigh,ikpt_loc) == jkpt)) then
!      
!      
!      
!      jcg = cgindex(jkpt,isppol)
!      
!      if (((dest-1) == mpi_enreg%me)) then
!      
!      cgq(:,1:count)  = cg(:,jcg+1:jcg+count)
!      cg1q(:,1:count) = cg1(:,jcg+1:jcg+count)
!      cg3q(:,1:count) = cg3(:,jcg+1:jcg+count)
!      
!      else
!      
!      tag = jkpt
!      count1 = npwarr(jkpt)*mband*nspinor
!      allocate(buffer(2,3*count1))
!      buffer(:,1:count1)            = cg(:,jcg+1:jcg+count1)
!      buffer(:,count1+1:2*count1)   = cg1(:,jcg+1:jcg+count1)
!      buffer(:,2*count1+1:3*count1) = cg3(:,jcg+1:jcg+count1)
!      
!      call MPI_SEND(buffer,2*3*count1,MPI_DOUBLE_PRECISION,(dest-1),tag,spaceComm,ierr)
!      
!      deallocate(buffer)
!      
!      end if
!      
!      end if
!      
!      end do          ! loop over dest
!      
!      end do          ! loop over jkpt

       if (ikpt_loc > mkmem) then
         ABI_DEALLOCATE(cgq)
         ABI_DEALLOCATE(cg1q)
         ABI_DEALLOCATE(cg3q)
         cycle
       end if

#else
!      no // over k-points

       cgq(:,1:count)  = cg(:,jj+1:jj+count)
       cg1q(:,1:count) = cg1(:,jj+1:jj+count)
       cg3q(:,1:count) = cg3(:,jj+1:jj+count)

#endif

!      Compute overlap matrices

       if (kptindex(2,ikpt2) == 0) then  ! no time-reversal symmetry

         do ipw = 1, npw_k

           jpw = pwind(ipw,ineigh,ikpt_loc)
           if (jpw /= 0) then

             do iband = 1, nband_occ
               do jband = 1, nband_occ

                 icg = ii + (iband-1)*npw_k + ipw
                 jcg = (jband-1)*npw_k1 + jpw

                 smat(1,iband,jband) = smat(1,iband,jband) + &
&                 cg(1,icg)*cgq(1,jcg) + cg(2,icg)*cgq(2,jcg)
                 smat(2,iband,jband) = smat(2,iband,jband) + &
&                 cg(1,icg)*cgq(2,jcg) - cg(2,icg)*cgq(1,jcg)

                 s13mat(1,iband,jband) = s13mat(1,iband,jband) + &
&                 cg1(1,icg)*cg3q(1,jcg) + cg1(2,icg)*cg3q(2,jcg)
                 s13mat(2,iband,jband) = s13mat(2,iband,jband) + &
&                 cg1(1,icg)*cg3q(2,jcg) - cg1(2,icg)*cg3q(1,jcg)

                 s1mat(1,iband,jband) = s1mat(1,iband,jband) + &
&                 cg1(1,icg)*cgq(1,jcg) + cg1(2,icg)*cgq(2,jcg) + &
&                 cg(1,icg)*cg1q(1,jcg) + cg(2,icg)*cg1q(2,jcg)
                 s1mat(2,iband,jband) = s1mat(2,iband,jband) + &
&                 cg1(1,icg)*cgq(2,jcg) - cg1(2,icg)*cgq(1,jcg) + &
&                 cg(1,icg)*cg1q(2,jcg) - cg(2,icg)*cg1q(1,jcg)

                 s3mat(1,iband,jband) = s3mat(1,iband,jband) + &
&                 cg3(1,icg)*cgq(1,jcg) + cg3(2,icg)*cgq(2,jcg) + &
&                 cg(1,icg)*cg3q(1,jcg) + cg(2,icg)*cg3q(2,jcg)
                 s3mat(2,iband,jband) = s3mat(2,iband,jband) + &
&                 cg3(1,icg)*cgq(2,jcg) - cg3(2,icg)*cgq(1,jcg) + &
&                 cg(1,icg)*cg3q(2,jcg) - cg(2,icg)*cg3q(1,jcg)

               end do
             end do

           end if

         end do   ! ipw

       else                              ! use time-reversal symmetry

         do ipw = 1,npw_k

           jpw = pwind(ipw,ineigh,ikpt_loc)
           if (jpw /= 0) then

             do iband = 1, nband_occ
               do jband = 1, nband_occ

                 icg = ii + (iband-1)*npw_k + ipw
                 jcg = (jband-1)*npw_k1 + jpw

                 smat(1,iband,jband) = smat(1,iband,jband) + &
&                 cg(1,icg)*cgq(1,jcg) - cg(2,icg)*cgq(2,jcg)
                 smat(2,iband,jband) = smat(2,iband,jband) - &
&                 cg(1,icg)*cgq(2,jcg) - cg(2,icg)*cgq(1,jcg)

                 s13mat(1,iband,jband) = s13mat(1,iband,jband) + &
&                 cg1(1,icg)*cg3q(1,jcg) - cg1(2,icg)*cg3q(2,jcg)
                 s13mat(2,iband,jband) = s13mat(2,iband,jband) - &
&                 cg1(1,icg)*cg3q(2,jcg) - cg1(2,icg)*cg3q(1,jcg)

                 s1mat(1,iband,jband) = s1mat(1,iband,jband) + &
&                 cg1(1,icg)*cgq(1,jcg) - cg1(2,icg)*cgq(2,jcg) + &
&                 cg(1,icg)*cg1q(1,jcg) - cg(2,icg)*cg1q(2,jcg)
                 s1mat(2,iband,jband) = s1mat(2,iband,jband) - &
&                 cg1(1,icg)*cgq(2,jcg) - cg1(2,icg)*cgq(1,jcg) - &
&                 cg(1,icg)*cg1q(2,jcg) - cg(2,icg)*cg1q(1,jcg)

                 s3mat(1,iband,jband) = s3mat(1,iband,jband) + &
&                 cg3(1,icg)*cgq(1,jcg) - cg3(2,icg)*cgq(2,jcg) + &
&                 cg(1,icg)*cg3q(1,jcg) - cg(2,icg)*cg3q(2,jcg)
                 s3mat(2,iband,jband) = s3mat(2,iband,jband) - &
&                 cg3(1,icg)*cgq(2,jcg) - cg3(2,icg)*cgq(1,jcg) - &
&                 cg(1,icg)*cg3q(2,jcg) - cg(2,icg)*cg3q(1,jcg)

               end do
             end do

           end if

         end do   ! ipw

       end if

       ABI_DEALLOCATE(cgq)
       ABI_DEALLOCATE(cg1q)
       ABI_DEALLOCATE(cg3q)

!      Compute qmat, the inverse of smat

       job = 1  ! compute inverse only
       qmat(:,:,:) = smat(:,:,:)

       call dzgefa(qmat,mband,nband_occ,ipvt,info)
       call dzgedi(qmat,mband,nband_occ,ipvt,det,zgwork,job)

!      DEBUG
!      write(100,*)
!      write(100,*)'ikpt = ',ikpt,'ineigh = ',ineigh
!      do iband = 1,nband_occ
!      do jband = 1,nband_occ
!      c1 = 0_dp ; c2 = 0_dp
!      do lband = 1,nband_occ
!      c1 = c1 + smat(1,iband,lband)*qmat(1,lband,jband) - &
!      &           smat(2,iband,lband)*qmat(2,lband,jband)
!      c2 = c2 + smat(1,iband,lband)*qmat(2,lband,jband) + &
!      &           smat(2,iband,lband)*qmat(1,lband,jband)
!      end do
!      write(100,'(2(2x,i2),2(2x,f16.9))')iband,jband,&
!      & c1,c2
!      end do
!      end do
!      ENDDEBUG



!      Accumulate sum over bands

       dotposr = 0_dp ; dotposi = 0_dp
       dotnegr = 0_dp ; dotnegi = 0_dp
       do iband = 1, nband_occ
         do jband = 1, nband_occ

           dotposr = dotposr + &
&           s13mat(1,iband,jband)*qmat(1,jband,iband) - &
&           s13mat(2,iband,jband)*qmat(2,jband,iband)
           dotposi = dotposi + &
&           s13mat(1,iband,jband)*qmat(2,jband,iband) + &
&           s13mat(2,iband,jband)*qmat(1,jband,iband)


           do lband = 1, nband_occ
             do lpband= 1, nband_occ

               z1(1) = s1mat(1,iband,jband)*qmat(1,jband,lband) - &
&               s1mat(2,iband,jband)*qmat(2,jband,lband)
               z1(2) = s1mat(1,iband,jband)*qmat(2,jband,lband) + &
&               s1mat(2,iband,jband)*qmat(1,jband,lband)

               z2(1) = s3mat(1,lband,lpband)*qmat(1,lpband,iband) - &
&               s3mat(2,lband,lpband)*qmat(2,lpband,iband)
               z2(2) = s3mat(1,lband,lpband)*qmat(2,lpband,iband) + &
&               s3mat(2,lband,lpband)*qmat(1,lpband,iband)

               dotnegr = dotnegr + &
&               z1(1)*z2(1) - z1(2)*z2(2)
               dotnegi = dotnegi + &
&               z1(1)*z2(2) + z1(2)*z2(1)

             end do   ! lpband
           end do   ! lband

         end do   ! jband
       end do   ! iband

       d3_aux(1,:) = d3_aux(1,:) + &
&       dk(:)*mvwtk(ineigh,ikpt)*dtset%wtk(ikpt)*(2_dp*dotposr-dotnegr)
       d3_aux(2,:) = d3_aux(2,:) + &
&       dk(:)*mvwtk(ineigh,ikpt)*dtset%wtk(ikpt)*(2_dp*dotposi-dotnegi)

     end do        ! End loop over neighbours


   end do      ! End loop over k-points

 end do  ! fab: end loop over spin 



 
 call xmpi_sum(d3_aux,spaceComm,ierr)


 ABI_DEALLOCATE(s13mat)
 ABI_DEALLOCATE(smat)
 ABI_DEALLOCATE(s1mat)
 ABI_DEALLOCATE(qmat)
 ABI_DEALLOCATE(ipvt)
 ABI_DEALLOCATE(s3mat)
 ABI_DEALLOCATE(zgwork)
 ABI_DEALLOCATE(bd_index)


!fab: I think that in the following we have to make a distinction:
!for the spin unpolarized case we leave the PEAD expression as it is, while
!in the spin polarized case we have simply to divide by 2 
!(see eq.19 di PRB 63,155107, eq. 7 di PRB 71,125107 and eq 13 di PRB 71, 125107...
!in this latter equation the 2 must be simply replaced by the sum over the spin components...
!and indeed we have inserted the loop over the spin,
!but there was a factor 2 already present in the routine due to spin degenracy that had to be removed)


 if (nsppol==1) then

!  Take minus the imaginary part

   d3_berry(1,:) = -1_dp*d3_aux(2,:)
   d3_berry(2,:) = d3_aux(1,:)

   d3_berry(2,:) = 0_dp

 else 

   d3_berry(1,:) = -1_dp*d3_aux(2,:)/2._dp
   d3_berry(2,:) = d3_aux(1,:)/2._dp

   d3_berry(2,:) = 0_dp/2._dp
   
 end if




!DEBUG
!write(100,*)'pead_nl_mv.f : d3_berry'
!write(100,*)'Perturbation',i1dir,i3dir
!write(100,*)
!write(100,*)'before transformation'
!write(100,*)'real part'
!write(100,'(3(2x,f20.9))')d3_berry(1,:)
!write(100,*)
!write(100,*)'imaginary part'
!write(100,'(3(2x,f20.9))')d3_berry(2,:)
!write(100,*)
!write(100,*)'after transformation'
!ENDDEBUG

!Compute the projection on the basis vectors of
!reciprocal space

 d3_aux(:,:) = 0_dp
 do ii = 1,3
   do jj = 1,3
     d3_aux(:,ii) = d3_aux(:,ii) + gmet(ii,jj)*d3_berry(:,jj)
   end do
 end do
 d3_berry(:,:) = d3_aux(:,:)

!Write out the berryphase part of the third order energy

 if (mpi_enreg%me == 0) then

   write(message,'(a,a,a)')ch10,&
&   ' Berryphase part of the third-order energy:',ch10
!  call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

   if (i1pert < natom + 1) then
     write(message,'(a,i3,a,i3)')&
&     '            j1: Displacement of atom ',i1pert,&
&     ' along direction ',i1dir
   else if (i1pert == natom + 2) then
     write(message,'(a,i3)')&
&     '            j1: homogenous electric field along direction ',i1dir
   end if
!  call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a)')&
&   '            j2: k-point derivative along direction i2dir '
!  call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

   if (i3pert < natom + 1) then
     write(message,'(a,i3,a,i3,a)')&
&     '            j3: Displacement of atom ',i3pert,&
&     ' along direction ',i3dir,ch10
   else if (i3pert == natom + 2) then
     write(message,'(a,i3,a)')&
&     '            j3: homogenous electric field along direction ',i3dir,ch10
   end if
!  call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

!  write(ab_out,'(5x,a5,8x,a9,5x,a14)')'i2dir','real part','imaginary part'
   write(std_out,'(5x,a5,8x,a9,5x,a14)')'i2dir','real part','imaginary part'
   do ii = 1,3
     write(std_out,'(7x,i1,3x,f16.9,3x,f16.9)')ii,&
&     d3_berry(1,ii),d3_berry(2,ii)
     write(std_out,'(7x,i1,3x,f16.9,3x,f16.9)')ii,&
&     d3_berry(1,ii),d3_berry(2,ii)
   end do

 end if    ! mpi_enreg%me == 0

!DEBUG
!write(100,*)'real part'
!write(100,'(3(2x,f20.9))')d3_berry(1,:)
!write(100,*)
!write(100,*)'imaginary part'
!write(100,'(3(2x,f20.9))')d3_berry(2,:)
!ENDDEBUG



!close(dtfil%unwff1)
!close(dtfil%unwff2)

 call status(0,dtfil%filstat,iexit,level,'exit          ')

!DEBUG
!write(std_out,*)' pead_nl_mv : exit '
!ENDDEBUG

end subroutine pead_nl_mv
!!***
