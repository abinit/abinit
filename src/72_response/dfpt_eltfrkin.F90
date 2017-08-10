!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_eltfrkin
!! NAME
!! dfpt_eltfrkin
!!
!! FUNCTION
!! Compute the frozen-wavefunction kinetic enegy contribution to the
!! elastic tensor
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DRH, DCA, XG, GM, AR, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=<G|Cnk>=Fourier coefficients of wavefunction
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha) (NOT NEEDED !)
!!  effmass=effective mass for electrons (1. in common case)
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=work array for coordinates of G vectors in basis
!!  kptns(3,nkpt)=coordinates of k points in terms of reciprocal space
!!   primitive translations
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem=number of k points treated by this node.
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimension for number of planewaves
!!  nband(nkpt*nsppol)=number of bands being considered per k point
!!  nkpt=number of k points
!!  ngfft(18)=contain all needed information about 3D FFT, i
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  npwarr(nkpt)=number of planewaves at each k point, and boundary
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for polarized
!!  occ(mband*nkpt*nsppol)=occupation numbers of bands (usually 2)
!!    at each k point
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  wtk(nkpt)=k point weights
!!
!! OUTPUT
!!  eltfrkin(6,6)=non-symmetrized kinetic energy contribution to the
!!                    elastic tensor
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      d2kindstr2,metric,sphereboundary,timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_eltfrkin(cg,eltfrkin,ecut,ecutsm,effmass,&
&  istwfk,kg,kptns,mband,mgfft,mkmem,mpi_enreg,&
&  mpw,nband,nkpt,ngfft,npwarr,nspinor,nsppol,occ,rprimd,wtk)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_xmpi
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_eltfrkin'
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_52_fft_mpi_noabirule
 use interfaces_72_response, except_this_one => dfpt_eltfrkin
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mgfft,mkmem,mpw,nkpt,nspinor,nsppol
 real(dp),intent(in) :: ecut,ecutsm,effmass
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: istwfk(nkpt),kg(3,mpw*mkmem),nband(nkpt*nsppol)
 integer,intent(in) :: ngfft(18),npwarr(nkpt)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),kptns(3,nkpt)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),rprimd(3,3),wtk(nkpt)
 real(dp),intent(out) :: eltfrkin(6,6)

!Local variables-------------------------------
!scalars
 integer :: bdtot_index,iband,icg,ierr,ii,ikg
 integer :: ikpt,index,ipw,isppol,istwf_k,jj,master,me,n1,n2
 integer :: n3,nband_k,nkinout,npw_k,spaceComm
 real(dp) :: ucvol
!arrays
 integer,allocatable :: gbound(:,:),kg_k(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),kpoint(3),rmet(3,3),tsec(2)
 real(dp),allocatable :: cwavef(:,:),ekinout(:)
 real(dp),allocatable :: eltfrkink(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

!Default for sequential use
 master=0
!Init mpi_comm
 spaceComm=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt

!Compute gmet, gprimd and ucvol from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 eltfrkin(:,:)=0.0_dp
 bdtot_index=0
 icg=0

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(cwavef,(2,mpw*nspinor))
 ABI_ALLOCATE(eltfrkink,(6,6))

!Define k-points distribution

!LOOP OVER SPINS
 do isppol=1,nsppol
   ikg=0

!  Loop over k points
   do ikpt=1,nkpt

     nband_k=nband(ikpt+(isppol-1)*nkpt)
     istwf_k=istwfk(ikpt)
     npw_k=npwarr(ikpt)

!    Skip this k-point if not the proper processor
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
       bdtot_index=bdtot_index+nband_k
       cycle
     end if

     ABI_ALLOCATE(gbound,(2*mgfft+8,2))
     kpoint(:)=kptns(:,ikpt)

     kg_k(:,:) = 0

!$OMP PARALLEL DO PRIVATE(ipw) SHARED(ikg,kg,kg_k,npw_k)
     do ipw=1,npw_k
       kg_k(1,ipw)=kg(1,ipw+ikg)
       kg_k(2,ipw)=kg(2,ipw+ikg)
       kg_k(3,ipw)=kg(3,ipw+ikg)
     end do

     call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)

     index=1+icg

     eltfrkink(:,:)=0.0_dp

     nkinout=6*6
     ABI_ALLOCATE(ekinout,(nkinout))
     ekinout(:)=zero

     do iband=1,nband_k

       if(mpi_enreg%proc_distrb(ikpt,iband,isppol) /= me) cycle

       cwavef(:,1:npw_k*nspinor)=cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)

       call d2kindstr2(cwavef,ecut,ecutsm,effmass,ekinout,gmet,gprimd,&
&       istwf_k,kg_k,kpoint,npw_k,nspinor)

       eltfrkink(:,:)=eltfrkink(:,:)+ occ(iband+bdtot_index)* reshape(ekinout(:), (/6,6/) )

     end do !iband

     ABI_DEALLOCATE(ekinout)

     eltfrkin(:,:)=eltfrkin(:,:)+wtk(ikpt)*eltfrkink(:,:)

     ABI_DEALLOCATE(gbound)

     bdtot_index=bdtot_index+nband_k

     if (mkmem/=0) then
!      Handle case in which kg, cg, are kept in core
       icg=icg+npw_k*nspinor*nband_k
       ikg=ikg+npw_k
     end if

   end do
 end do  ! End loops on isppol and ikpt

!Fill in lower triangle
 do jj=2,6
   do ii=1,jj-1
     eltfrkin(jj,ii)=eltfrkin(ii,jj)
   end do
 end do

!Accumulate eltfrkin on all proc.
 call timab(48,1,tsec)
 call xmpi_sum(eltfrkin,spaceComm,ierr)
 call timab(48,2,tsec)

 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(eltfrkink)
 ABI_DEALLOCATE(kg_k)

 DBG_EXIT("COLL")

end subroutine dfpt_eltfrkin
!!***
