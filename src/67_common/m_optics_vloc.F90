!!****m* ABINIT/m_optics_vloc
!! NAME
!!  m_optics_vloc
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2010-2020 ABINIT group (SM,VR,FJ,MT)
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

module m_optics_vloc

 use defs_basis
 use m_abicore
 use m_errors
 use m_wffile
 use m_xmpi
 use m_hdr
 use m_dtset
 use m_dtfil

use defs_abitypes,   only : MPI_type
 use m_time,         only : timab
 use m_io_tools,     only : get_unit
 use m_mpinfo,       only : proc_distrb_cycle

 implicit none

 private
!!***

 public :: optics_vloc
!!***

contains
!!***

!!****f* ABINIT/optics_vloc
!! NAME
!! optics_vloc
!!
!! FUNCTION
!! Compute matrix elements need for optical conductivity in a LOCAL potential
!! and store them in a file.
!! Matrix elements = <Phi_i|Nabla|Phi_j>
!!
!! INPUTS
!!  cg(2,mcg)=planewave coefficients of wavefunctions.
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  fildata= name of the output file
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mkmem =number of k points treated by this node.
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimensioned size of npw.
!!  nkpt=number of k points.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!
!! OUTPUT
!!  (only writing in a file)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_outscfcv
!!
!! CHILDREN
!!      hdr_io,timab,wffclose,wffopen,xmpi_exch,xmpi_sum_master
!!
!! SOURCE

 subroutine optics_vloc(cg,dtfil,dtset,eigen0,gprimd,hdr,kg,mband,mcg,mkmem,mpi_enreg,mpw,&
&                       nkpt,npwarr,nsppol)
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcg,mkmem,mpw,nkpt,nsppol
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(in) :: mpi_enreg
 type(hdr_type),intent(inout) :: hdr
!arrays
 integer,intent(in) :: kg(3,mpw*mkmem),npwarr(nkpt)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: eigen0(mband*nkpt*nsppol)
 real(dp),intent(inout) :: cg(2,mcg)

!Local variables-------------------------------
!scalars
 integer :: iomode,bdtot_index,cplex,etiq,fformopt,ib,icg,ierr,ikg,ikpt
 integer :: ipw,isppol,istwf_k,iwavef,jb,jwavef
 integer :: me,me_kpt,my_nspinor,nband_k,npw_k,sender,ount
 integer :: spaceComm_band,spaceComm_bandfftspin,spaceComm_fft,spaceComm_k
 logical :: mykpt
 real(dp) :: cgnm1,cgnm2
 character(len=500) :: msg
!arrays
 integer :: tmp_shape(3)
 integer,allocatable :: kg_k(:,:)
 real(dp) :: kpoint(3),tsec(2)
 real(dp),allocatable :: eig0_k(:),kpg_k(:,:)
 real(dp),allocatable :: psinablapsi(:,:,:,:),tnm(:,:,:,:)
 type(wffile_type) :: wff1
! ************************************************************************

 DBG_ENTER("COLL")

!Compatibility tests

!----------------------------------------------------------------------------------
!2- Computation of <psi_n|-i.nabla|psi_m> for each k
!----------------------------------------------------------------------------------

!Init parallelism
 if (mpi_enreg%paral_kgb==1) then
   spaceComm_k=mpi_enreg%comm_kpt
   spaceComm_fft=mpi_enreg%comm_fft
   spaceComm_band=mpi_enreg%comm_band
   spaceComm_bandfftspin=mpi_enreg%comm_bandspinorfft
   me_kpt=mpi_enreg%me_kpt
 else
   spaceComm_band=0;spaceComm_fft=0;spaceComm_bandfftspin=0
   spaceComm_k=mpi_enreg%comm_cell
   me=xmpi_comm_rank(spaceComm_k)
   me_kpt=me
 end if
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

!Initialize main variables
 ABI_MALLOC(psinablapsi,(2,3,mband,mband))
 psinablapsi=zero

 iomode= IO_MODE_FORTRAN_MASTER
 fformopt=612
 ount = get_unit()
 call WffOpen(iomode,spaceComm_k,dtfil%fnameabo_app_opt,ierr,wff1,0,me,ount)
 call hdr_io(fformopt,hdr,2,wff1)
 !if (me == 0) then
 !call hdr_fort_write(hdr, wff1%unwff, fformopt,ierr)
 !ABI_CHECK(ierr /= 0, "hdr_fort_write returned ierr = 0")
 !end if

!LOOP OVER SPINS
 icg=0
 do isppol=1,nsppol

!  LOOP OVER k POINTS
   ikg=0
   do ikpt=1,nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
     etiq=ikpt+(isppol-1)*nkpt
     if (me==0) then
       ABI_MALLOC(eig0_k,(nband_k))
       eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
     end if

!    Select my k-points
     mykpt=.true.
     mykpt=(.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_kpt)))
     if (mykpt) then

!      Allocations depending on k-point
       kpoint(:)=dtset%kptns(:,ikpt)
       istwf_k=dtset%istwfk(ikpt)
       npw_k=npwarr(ikpt)
       cplex=2;if (istwf_k>1) cplex=1
       ABI_MALLOC(kg_k,(3,npw_k))
       ABI_MALLOC(kpg_k,(npw_k*dtset%nspinor,3))

!      Get G-vectors for this k-point
       kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
       ikg=ikg+npw_k

!      Calculation of k+G in cartesian coordinates
       do ipw=1,npw_k
         kpg_k(ipw,1)=(kpoint(1)+kg_k(1,ipw))*gprimd(1,1)&
&         +(kpoint(2)+kg_k(2,ipw))*gprimd(1,2)&
&         +(kpoint(3)+kg_k(3,ipw))*gprimd(1,3)
         kpg_k(ipw,2)=(kpoint(1)+kg_k(1,ipw))*gprimd(2,1)&
&         +(kpoint(2)+kg_k(2,ipw))*gprimd(2,2)&
&         +(kpoint(3)+kg_k(3,ipw))*gprimd(2,3)
         kpg_k(ipw,3)=(kpoint(1)+kg_k(1,ipw))*gprimd(3,1)&
&         +(kpoint(2)+kg_k(2,ipw))*gprimd(3,2)&
&         +(kpoint(3)+kg_k(3,ipw))*gprimd(3,3)
       end do !ipw
       kpg_k=two_pi*kpg_k
       if (dtset%nspinor==2) kpg_k(npw_k+1:2*npw_k,1:3)=kpg_k(1:npw_k,1:3)
       ABI_FREE(kg_k)

!      2-A Computation of <psi_tild_n|-i.nabla|psi_tild_m>
!      ----------------------------------------------------------------------------------
!      Computation of (C_nk^*)*C_mk*(k+g) in cartesian coordinates

       ABI_MALLOC(tnm,(2,3,nband_k,nband_k))
       tnm=zero

!      Loops on bands
       do jb=1,nband_k
         jwavef=(jb-1)*npw_k*my_nspinor+icg
         if (mpi_enreg%paral_kgb/=1) then
           tmp_shape = shape(mpi_enreg%proc_distrb)
           if (ikpt > tmp_shape(1)) then
             msg='  ikpt out of bounds '
             ABI_BUG(msg)
           end if
           if (abs(mpi_enreg%proc_distrb(ikpt,jb,isppol)-me_kpt)/=0) cycle
         end if
         do ib=1,jb
           iwavef=(ib-1)*npw_k*my_nspinor+icg

!          Computation of (C_nk^*)*C_mk*(k+g) in cartesian coordinates
           if (cplex==1) then
             do ipw=1,npw_k*my_nspinor
               cgnm1=cg(1,ipw+iwavef)*cg(1,ipw+jwavef)
               tnm(1,1:3,ib,jb)=tnm(1,1:3,ib,jb)+cgnm1*kpg_k(ipw,1:3)
             end do
           else
             do ipw=1,npw_k*my_nspinor
               cgnm1=cg(1,ipw+iwavef)*cg(1,ipw+jwavef)+cg(2,ipw+iwavef)*cg(2,ipw+jwavef)
               cgnm2=cg(1,ipw+iwavef)*cg(2,ipw+jwavef)-cg(2,ipw+iwavef)*cg(1,ipw+jwavef)
               tnm(1,1:3,ib,jb)=tnm(1,1:3,ib,jb)+cgnm1*kpg_k(ipw,1:3)
               tnm(2,1:3,ib,jb)=tnm(2,1:3,ib,jb)+cgnm2*kpg_k(ipw,1:3)
             end do
           end if

!          Second half of the (n,m) matrix
           if (ib/=jb) then
             tnm(1,1:3,jb,ib)= tnm(1,1:3,ib,jb)
             tnm(2,1:3,jb,ib)=-tnm(2,1:3,ib,jb)
           end if

         end do ! ib
       end do ! jb

!      Reduction in case of parallelism
       if (mpi_enreg%paral_kgb == 1) then
         call timab(48,1,tsec)
         call xmpi_sum_master(tnm,0,spaceComm_bandfftspin,ierr)
         call timab(48,2,tsec)
       end if

       psinablapsi(:,:,:,:)=tnm(:,:,:,:)

       ABI_FREE(tnm)

       if (mkmem/=0) then
         icg = icg + npw_k*my_nspinor*nband_k
       end if

       ABI_FREE(kpg_k)

       if (me==0) then
         write(ount)(eig0_k(ib),ib=1,nband_k)
         write(ount)((psinablapsi(1:2,1,ib,jb),ib=1,nband_k),jb=1,nband_k)
         write(ount)((psinablapsi(1:2,2,ib,jb),ib=1,nband_k),jb=1,nband_k)
         write(ount)((psinablapsi(1:2,3,ib,jb),ib=1,nband_k),jb=1,nband_k)
       elseif (mpi_enreg%me_band==0.and.mpi_enreg%me_fft==0) then
         call xmpi_exch(psinablapsi,etiq,me_kpt,psinablapsi,0,spaceComm_k,ierr)
       end if

     elseif (me==0) then
       sender=minval(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol))
       call xmpi_exch(psinablapsi,etiq,sender,psinablapsi,0,spaceComm_k,ierr)
       write(ount)(eig0_k(ib),ib=1,nband_k)
       write(ount)((psinablapsi(1:2,1,ib,jb),ib=1,nband_k),jb=1,nband_k)
       write(ount)((psinablapsi(1:2,2,ib,jb),ib=1,nband_k),jb=1,nband_k)
       write(ount)((psinablapsi(1:2,3,ib,jb),ib=1,nband_k),jb=1,nband_k)
     end if ! mykpt

     bdtot_index=bdtot_index+nband_k
     if (me==0)  then
       ABI_FREE(eig0_k)
     end if
!    End loop on spin,kpt
   end do ! ikpt
 end do !isppol

!Close file
 call WffClose(wff1,ierr)

!Datastructures deallocations
 ABI_FREE(psinablapsi)

 DBG_EXIT("COLL")

end subroutine optics_vloc
!!***

end module m_optics_vloc
!!***
