!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmv
!! NAME
!! initmv
!!
!! FUNCTION
!! Initialize finite difference calculation of the ddk im dfptnl_mv.f
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)> = all input variables in this dataset
!!  gmet(3,3) = reciprocal space metric tensor in bohr**-2
!!  kg(3,mpw*mkmem) = reduced (integer) coordinates of G vecs in basis sphere
!!  kneigh(30,nkpt2) = index of the neighbours of each k-point
!!  kg_neigh(30,nkpt2,3) = necessary to construct the vector joining a k-point 
!!                         to its nearest neighbour in case of a single k-point, 
!!                         a line of k-points or a plane of k-points. 
!!                         See getshell.F90 for details
!!  kptindex(2,nkpt3)= index of the k-points in the reduced BZ
!!                     related to a k-point in the full BZ
!!  kpt3(3,nkpt3) = reduced coordinates of k-points in the full BZ
!!  mband = maximum number of bands
!!  mkmem = number of k points which can fit in memory
!!  mpi_enreg = informations about MPI parallelization
!!  mpw = maximum number of plane waves
!!  nband(nkpt*nsppol)=number of bands at each k point, for each polarization
!!  nkpt2 = number of k-points in the reduced BZ
!!  nkpt3 = number of k-points in the full BZ
!!  nneigh = total number of neighbours required to evaluate the finite
!!          difference formula
!!  npwarr(nkpt2)=number of planewaves at each k point
!!  nsppol = number of spin polarizations
!!  occ(mband*nkpt*nsppol) = occupation number for each band for each k
!!
!! OUTPUT
!! cgindex(nkpt2,nsppol) = for each k-point, cgindex tores the location
!!                         of the WF in the cg array
!!        me = index of the current processor
!!        ineigh = index of a neighbour
!!        ikpt_loc = index of the iteration on ikpt on the current processor
!!        ikpt_rbz = index of a k-point in the reduced BZ
!! pwind(mpw,nneigh,mkmem) = array used to compute the overlap matrix smat
!!                           between k-points
!!                           (see initberry.f for more explanations)
!!
!! COMMENTS
!!
!! PARENTS
!!      nonlinear
!!
!! CHILDREN
!!      kpgio,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initmv(cgindex,dtset,gmet,kg,kneigh,kg_neigh,kptindex,&
&  kpt3,mband,mkmem,mpi_enreg,mpw,nband,nkpt2,&
&  nkpt3,nneigh,npwarr,nsppol,occ,pwind)

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initmv'
 use interfaces_32_util
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mkmem,mpw,nkpt2,nkpt3,nneigh,nsppol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: kg(3,mpw*mkmem),kneigh(30,nkpt2),kg_neigh(30,nkpt2,3)
 integer,intent(in) :: nband(nkpt2*nsppol),npwarr(nkpt2),kptindex(2,nkpt3)
 integer,intent(out) :: cgindex(nkpt2,nsppol),pwind(mpw,nneigh,mkmem)
 real(dp),intent(in) :: gmet(3,3),kpt3(3,nkpt3),occ(mband*nkpt2*nsppol)

!Local variables-------------------------------
!scalars
 integer :: flag,iband,icg,ierr,ikg,ikg1,ikpt,ikpt2,ikpt_loc,ikpt_rbz
 integer :: index,ineigh,ipw,isppol,jpw,nband_k,mband_occ,mband_occ_k,npw_k
 integer :: npw_k1,orig,spaceComm
 real(dp) :: ecut_eff,sdeg
 character(len=500) :: message
!arrays
 integer :: dg(3)
 integer,allocatable :: kg1(:,:),kg1_k(:,:),npwar1(:),npwtot(:)
 real(dp) :: dk(3),dk_(3)
 real(dp),allocatable :: kpt1(:,:)

!************************************************************************

!DEBUG
!write(std_out,*)' initmv : enter '
!ENDDEBUG

 if (xmpi_paral== 1) then
!  BEGIN TF_CHANGES
   spaceComm=mpi_enreg%comm_cell
!  END TF_CHANGES
   mpi_enreg%kpt_loc2ibz_sp(:,:,:) = 0
   mpi_enreg%mkmem(:) = 0
 end if

 ecut_eff = dtset%ecut*(dtset%dilatmx)**2
 ABI_ALLOCATE(kg1_k,(3,mpw))
 ABI_ALLOCATE(kg1,(3,mkmem*mpw))
 ABI_ALLOCATE(kpt1,(3,nkpt2))
 ABI_ALLOCATE(npwar1,(nkpt2))
 ABI_ALLOCATE(npwtot,(nkpt2))
 kg1_k(:,:) = 0
 pwind(:,:,:) = 0
 cgindex(:,:) = 0

!Compute the number of occupied bands.
!Check that it is the same for every k-point and that
!nband(ikpt) is equal to this value

 if (nsppol == 1) then
   sdeg = two
 else if (nsppol == 2) then
   sdeg = one
 end if

!DEBUG
!write(std_out,*)' list of nband '
!do isppol = 1, nsppol
!do ikpt = 1, nkpt2
!
!nband_k = nband(ikpt + (isppol - 1)*nkpt2)
!write(std_out,*)' isppol, ikpt, nband_k=',isppol, ikpt, nband_k
!end do
!end do
!ENDDEBUG


 




 index = 0
 do isppol = 1, nsppol
   do ikpt = 1, nkpt2

     mband_occ_k = 0
     nband_k = nband(ikpt + (isppol - 1)*nkpt2)

     do iband = 1, nband_k
       index = index + 1
       if (abs(occ(index) - sdeg) < tol8) mband_occ_k = mband_occ_k + 1
     end do

     if (nband_k /= mband_occ_k) then
       write(message,'(a,a,a)')&
&       '  In a non-linear response calculation, nband must be equal ',ch10,&
&       '  to the number of valence bands.'
       MSG_ERROR(message)
     end if

!    Note that the number of bands can be different for spin up and spin down
     if (ikpt > 1) then
       if (mband_occ /= mband_occ_k) then
         message = 'The number of valence bands is not the same for every k-point'
         MSG_ERROR(message)
       end if
     else
       mband_occ = mband_occ_k
     end if

   end do                ! close loop over ikpt
 end do                ! close loop over isppol

!Find the location of each wavefunction

 icg = 0
 do isppol = 1, nsppol
   do ikpt = 1, nkpt2


!    fab: inserted the shift due to the spin...

     nband_k = dtset%nband(ikpt+(isppol - 1)*nkpt2)
     npw_k = npwarr(ikpt)

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,mpi_enreg%me)) cycle

     cgindex(ikpt,isppol) = icg
     icg = icg + dtset%nspinor*npw_k*nband_k

   end do
 end do


!Build pwind 

 do ineigh = 1, nneigh

   do ikpt = 1, nkpt2
     ikpt2  = kneigh(ineigh,ikpt)
     ikpt_rbz = kptindex(1,ikpt2)   ! index of the k-point in the reduced BZ
     kpt1(:,ikpt) = dtset%kptns(:,ikpt_rbz)
   end do

!  Set up the basis sphere of plane waves at kpt1
   kg1(:,:) = 0
   call kpgio(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,kg1,&
&   kpt1,mkmem,dtset%nband,nkpt2,'PERS',mpi_enreg,mpw,&
&   npwar1,npwtot,dtset%nsppol)

   ikg = 0 ; ikg1 = 0 ; ikpt_loc = 0

   if(dtset%nsppol/=1)then
     if(mpi_enreg%nproc/=1)then
       message = ' At present, non-linear response calculations for spin-polarized system cannot be done in parallel.'
       MSG_ERROR(message)
     else
       isppol=1
     end if
   else
     isppol=1
   end if

   do ikpt = 1, nkpt2

     nband_k = dtset%nband(ikpt+(isppol - 1)*nkpt2)
     ikpt2  = kneigh(ineigh,ikpt)
     ikpt_rbz = kptindex(1,ikpt2)   ! index of the k-point in the reduced BZ

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,-1,mpi_enreg%me))  cycle

     ikpt_loc = ikpt_loc + 1

     mpi_enreg%kpt_loc2ibz_sp(mpi_enreg%me, ikpt_loc, 1) = ikpt

     flag = 0
     npw_k = npwarr(ikpt)
     npw_k1 = npwarr(ikpt_rbz)
     dk_(:) = kpt3(:,ikpt2) - dtset%kptns(:,ikpt)
     dk(:)  = dk_(:) - nint(dk_(:)) + real(kg_neigh(ineigh,ikpt,:),dp)
     dg(:)  = nint(dk(:) - dk_(:))


     if (kptindex(2,ikpt2) == 0) then
       kg1_k(:,1:npw_k1) = kg1(:,ikg1+1:ikg1+npw_k1)
       if (dg(1)==0.and.dg(2)==0.and.dg(3)==0) flag = 1
     else
       kg1_k(:,1:npw_k1) = -1*kg1(:,ikg1+1:ikg1+npw_k1)
     end if

     orig = 1
     do ipw = 1, npw_k
       do jpw = orig, npw_k1

         if ((kg(1,ikg + ipw) == kg1_k(1,jpw) - dg(1)).and. &
&         (kg(2,ikg + ipw) == kg1_k(2,jpw) - dg(2)).and. &
&         (kg(3,ikg + ipw) == kg1_k(3,jpw) - dg(3)))  then

           pwind(ipw,ineigh,ikpt_loc) = jpw
           if (flag == 1)  orig = jpw + 1
           exit

         end if

       end do
     end do

     ikg = ikg + npw_k
     ikg1 = ikg1 + npw_k1

   end do     ! close loop over k-points

 end do    ! close loop over ineigh
 mpi_enreg%mkmem(mpi_enreg%me) = mkmem

 call xmpi_sum(mpi_enreg%kpt_loc2ibz_sp,spaceComm,ierr)
 call xmpi_sum(mpi_enreg%mkmem,spaceComm,ierr)
 

!----------------------------------------------------------------------------

 ABI_DEALLOCATE(kg1)
 ABI_DEALLOCATE(kg1_k)
 ABI_DEALLOCATE(kpt1)
 ABI_DEALLOCATE(npwar1)
 ABI_DEALLOCATE(npwtot)

end subroutine initmv
!!***
