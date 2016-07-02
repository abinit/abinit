!{\src2tex{textfont=tt}}
!!****f* ABINIT/cgq_builder
!! NAME
!! cgq_builder
!!
!! FUNCTION
!! This routine locates cgq for efield calculations, especially
!! for // case
!!
!! COPYRIGHT
!! Copyright (C) 2003-2016 ABINIT  group 
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  berryflag = logical flag determining use of electric field variables
!!  cg(2,mcg)=planewave coefficients of wavefunctions.
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  ikpt=index of current k kpt 
!!  ikpt_loc=index of k point on current processor (see vtorho.F90)
!!  isspol=value of spin polarization currently treated
!!  me_distrb=current value from spaceComm_distrb (see vtorho.F90)
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcgq=size of cgq array (see vtorho.F90)
!!  mkgq=size of pwnsfacq array (see vtorho.F90)
!!  my_nspinor=nspinor value determined by current // set up
!!  nband_k=number of bands at each k point
!!  nproc_distrb=nproc from spaceComm_distrb (see vtorho.F90)
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  spaceComm_distrb=comm_cell from mpi_enreg
!!
!! OUTPUT
!!  cgq(2,mcgq)=planewave coefficients of wavenfunctions adjacent to cg at ikpt
!!  pwnsfacq(2,mkgq)=phase factors for non-symmorphic translations for cg's adjacent to cg(ikpt)
!!
!! SIDE EFFECTS
!! Input/Output
!!   dtefield <type(efield_type)> = efield variables
!!   mpi_enreg=information about MPI parallelization
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      timab,xmpi_recv,xmpi_send
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine cgq_builder(berryflag,cg,cgq,dtefield,dtset,ikpt,ikpt_loc,isppol,mcg,mcgq,&
&                      me_distrb,mkgq,mpi_enreg,my_nspinor,nband_k,nproc_distrb,&
&                      npwarr,pwnsfac,pwnsfacq,pwind_alloc,spaceComm_distrb)

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_efield
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cgq_builder'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: ikpt,ikpt_loc,isppol,me_distrb,mcg,mcgq,mkgq,my_nspinor,nband_k
 integer,intent(in) :: nproc_distrb,pwind_alloc,spaceComm_distrb
 logical,intent(in) :: berryflag
 type(dataset_type), intent(in) :: dtset
 type(efield_type), intent(inout) :: dtefield
 type(MPI_type), intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: npwarr(dtset%nkpt) 
 real(dp),intent(in) :: cg(2,mcg),pwnsfac(2,pwind_alloc)
 real(dp),intent(out) :: cgq(2,mcgq),pwnsfacq(2,mkgq)

!Local variables -------------------------
!scalars
 integer :: count,count1,icg1,icg2,dest,his_source
 integer :: idir,ierr,ifor,ikg1,ikg2,ikptf,ikpt1f,ikpt1i
 integer :: jkpt,jkpt1i,jkptf,jkpt1f,jsppol,my_source,npw_k1,tag
!arrays
 integer,allocatable :: flag_send(:,:), flag_receive(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: buffer(:,:)

! *************************************************************************

!DEBUG
!write(std_out,'(a)')'cgq_builder enter'
!DEBUG

 if (mcgq==0.or.mkgq==0) return
 
 call timab(983,1,tsec)

!Test compatbility of berryflag
 if (berryflag) then
   ABI_ALLOCATE(flag_send,(0:nproc_distrb-1,dtefield%fnkpt))
 end if
 ABI_ALLOCATE(flag_receive,(dtset%nkpt))
 flag_send(:,:) = 0
 flag_receive(:) = 0

 if (berryflag) ikptf = dtefield%i2fbz(ikpt)

 do idir = 1, 3

!  skip idir values for which efield_dot(idir) = 0
   if (berryflag .and. abs(dtefield%efield_dot(idir)) < tol12 ) cycle 

   do ifor = 1, 2

     if(berryflag) then
       dtefield%sflag(:,ikpt + dtset%nkpt*(isppol - 1),ifor,idir) = 0
       ikpt1f = dtefield%ikpt_dk(ikptf,ifor,idir)
       ikpt1i = dtefield%indkk_f2ibz(ikpt1f,1)
     end if

     npw_k1 = npwarr(ikpt1i)
     count = npw_k1*my_nspinor*nband_k
     my_source = mpi_enreg%proc_distrb(ikpt1i,1,isppol)

     do dest = 0, nproc_distrb-1

       if ((dest==me_distrb).and.(ikpt_loc <= dtset%mkmem)) then
!        I am dest and have something to do

         if ( my_source == me_distrb ) then
!          I am destination and source

           if(berryflag) then 
             ikg1 = dtefield%fkgindex(ikpt1f)
             ikg2 = dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*dtset%nkpt)
             icg1 = dtefield%cgindex(ikpt1i,isppol)
             icg2 = dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*dtset%nkpt)
           end if

           pwnsfacq(:,ikg2 + 1:ikg2 + npw_k1) = pwnsfac(:,ikg1 + 1:ikg1 + npw_k1)
           cgq(:,icg2 + 1:icg2 + count) = cg(:,icg1 + 1:icg1 + count)

         else !  I am the destination but not the source -> receive
!          receive pwnsfacq
           if(berryflag) then
             tag = ikpt1f + (isppol - 1)*dtefield%fnkpt
             ikg1 = dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*dtset%nkpt)
           end if
           ABI_ALLOCATE(buffer,(2,npw_k1))
           call xmpi_recv(buffer,my_source,tag,spaceComm_distrb,ierr)
           pwnsfacq(:,ikg1+1:ikg1+npw_k1) = buffer(:,1:npw_k1)
           ABI_DEALLOCATE(buffer)

!          receive cgq if necessary
           if(flag_receive(ikpt1i) == 0) then
             ABI_ALLOCATE(buffer,(2,count))
             tag = ikpt1i + (isppol - 1)*dtset%nkpt
             call xmpi_recv(buffer,my_source,tag,spaceComm_distrb,ierr)
             if(berryflag) icg1 = dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*dtset%nkpt)
             cgq(:,icg1+1:icg1+count) = buffer(:,1:count)
             ABI_DEALLOCATE(buffer)
             flag_receive(ikpt1i) = 1
           end if ! end if flag_receive == 0
         end if ! end tasks if I am the destination

       else if (ikpt_loc <= mpi_enreg%mkmem(dest)) then  ! dest != me and the dest has a k-point to treat

!        jkpt is the kpt which is being treated by dest (in ibz)
!        jsppol is his isppol
         jkpt = mpi_enreg%kpt_loc2ibz_sp(dest, ikpt_loc,1)
         jsppol = mpi_enreg%kpt_loc2ibz_sp(dest, ikpt_loc,2)

         if(jkpt > 0 .and. jsppol > 0) then

           if(berryflag) then 
             jkptf = dtefield%i2fbz(jkpt)
             jkpt1f = dtefield%ikpt_dk(jkptf,ifor,idir)
             jkpt1i = dtefield%indkk_f2ibz(jkpt1f,1)
           end if
           his_source = mpi_enreg%proc_distrb(jkpt1i,1,jsppol)

           if (his_source == me_distrb) then

!            send
!            pwnsfacq
             if(berryflag) then
               ikg1 = dtefield%fkgindex(jkpt1f)
               tag = jkpt1f + (jsppol - 1)*dtefield%fnkpt
             end if
             count1 = npwarr(jkpt1i)
             ABI_ALLOCATE(buffer,(2,count1))
             buffer(:,1:count1)  = pwnsfac(:,ikg1+1:ikg1+count1)
             call xmpi_send(buffer,dest,tag,spaceComm_distrb,ierr)
             ABI_DEALLOCATE(buffer)

!            send cgq if necessary
             if(flag_send(dest, jkpt1i)==0) then
               if(berryflag) icg1 = dtefield%cgindex(jkpt1i,jsppol)
               tag = jkpt1i + (jsppol - 1)*dtset%nkpt
               count1 = npwarr(jkpt1i)*nband_k*my_nspinor
               ABI_ALLOCATE(buffer,(2,count1))
               buffer(:,1:count1)  = cg(:,icg1+1:icg1+count1)
               call xmpi_send(buffer,dest,tag,spaceComm_distrb,ierr)
               ABI_DEALLOCATE(buffer)
               flag_send(dest, jkpt1i)=1
             end if ! if send cgq

           end if ! end check that his_source == me

         end if ! end check on jkpt > 0 and jsppol > 0

       end if ! end check on me = dest else if me != dest

     end do ! end loop over dest = 0, nproc-1

   end do !end loop over ifor

 end do !end loop over idir

 call timab(983,2,tsec)

 ABI_DEALLOCATE(flag_send)
 ABI_DEALLOCATE(flag_receive)

!DEBUG
!write(std_out,'(a)')'cgq_builder exit'
!END_DEBUG

end subroutine cgq_builder
!!***
