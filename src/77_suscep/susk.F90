!{\src2tex{textfont=tt}}
!!****f* ABINIT/susk
!! NAME
!! susk
!!
!! FUNCTION
!! Compute the contribution of one k point to the susceptibility matrix
!! from input wavefunctions, band occupations, and k point wts.
!! Include the usual sum-over-state terms, but also the
!! corrections due to the change of the Fermi level in the metallic
!! case, as well as implicit sum over higher lying conduction
!! states, thanks to the closure relation (referred to as an extrapolation).
!! Compared to the routine suskmm, there is no particular attention
!! to the use of the memory, so the code is simpler.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms
!!  bdtot_index=index for the number of the band
!!  cg(2,mcg)=wfs in G space
!!  cprj_k(natom,nspinor*nband_k)= wave functions projected with non-local projectors:
!!                                 cprj_k=<p_i|Cnk> where p_i is a non-local projector.
!!  doccde(mband*nkpt*nsppol)=derivative of occupancies wrt
!!           the energy for each band and k point
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  extrap: if==1, the closure relation (an extrapolation) must be used
!!  gbound(2*mgfftdiel+8,2)=G sphere boundary for going from WF sphere to
!!      medium size FFT grid
!!  gbound_diel(2*mgfftdiel+8,2)=G sphere boundary for going from medium size
!!      FFT grid to small sphere.
!!  gylmg_diel(npwdiel,lmax_diel,ntypat*usepaw)= -PAW only- Fourier transform of g_l(r).Y_ml(r) shape functions
!!                                               for dielectric matrix
!!  icg=index for cg
!!  ikpt=number of the k point
!!  isp=number of the current spin
!!  istwf_k=input option parameter that describes the storage of wfs
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  kg_k(3,npw_k)=coordinates of planewaves in basis sphere.
!!  lmax_diel=1+max. value of l angular momentum used for dielectric matrix
!!  mband=maximum number of bands
!!  mcg=dimension of cg
!!  mgfftdiel=maximum size of 1D FFTs, for the computation of
!!     the dielectric matrix
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell
!!  nband_k=number of bands at this k point for that spin polarization
!!  ndiel4,ndiel5,ndiel6= FFT dimensions, modified to avoid cache trashing
!!  neglect_pawhat=1 if PAW contribution from hat density (compensation charge)
!!                 has to be neglected (to be used when only an estimation of
!!                 suscep. matrix has to be evaluated, i.e. for SCF precondictioning)
!!  ngfftdiel(18)=contain all needed information about 3D FFT, for dielectric matrix,
!!    see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkpt=number of k points
!!  npwdiel=third and fifth dimension of the susmat array.
!!  npw_k=number of plane waves at this k point
!!  nspden=number of spin-density components
!!  nspden_eff=number of spin-density components actually computed in sussceptibility
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  occ(mband*nkpt*nsppol)=
!!          occupation numbers for each band (usually 2.0) at each k point
!!  occopt=option for occupancies
!!  occ_deavg(mband)=factor for extrapolation (occup. divided by an energy gap)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ph3d_diel(2,npwdiel,natom*usepaw)=3-dim structure factors, for each atom and plane wave, for dielectric matrix
!!  typat(natom)=type (integer) for each atom
!!  ucvol=unit cell volume (Bohr**3)
!!  usepaw=flag for PAW
!!  wtk(nkpt)=k point weights (they sum to 1.0)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! These quantities are accumulated in this routine:
!!  drhode(2,npwdiel,nspden_eff)=weighted density, needed to compute the
!!   effect of change of fermi energy
!!  rhoextrap(ndiel4,ndiel5,ndiel6,nspinor)=density-like array, needed for the
!!   extrapolation procedure.
!!  sumdocc=sum of weighted occupation numbers, needed to compute the
!!   effect of change of fermi energy
!!  susmat(2,npwdiel,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!
!! NOTES
!! Band-fft parallel treatment: Each processor will treat his own band, but susmat will be known by all.
!! This means that cg will not have the same meaning in sequential or parallel mode.
!! In parallel mode, it will contain the set of all bands treated by the currrent processor.
!! To achieve this, the argument cg has been replaced by cg_mpi, with the "target" attribute.
!! In sequential mode, the pointer cg will point towards cg_mpi. In parallel mode, cg will point
!! to a new array cg_local, containing the bands treated by the currrent processor.
!! This allows to minimize the overhead incurred by the parallelization  of the sequential version.
!! A similar treatment is performed on kg_k, npw_k.
!! A future version might have objects like kg_k_gather as arguments, instead of computing them.
!! This is in slight violation of programming rules, but I think it is safe, since the pointers remain local
!! GZ
!! PARENTS
!!      suscep_stat
!!
!! CHILDREN
!!      destroy_mpi_enreg,fourwf,init_distribfft_seq,initmpi_seq,pawsushat
!!      sphereboundary,timab,xmpi_allgather,xmpi_allgatherv,xmpi_alltoallv
!!      xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine susk(atindx,bdtot_index,cg_mpi,cprj_k,doccde,drhode,eigen,extrap,gbound,&
&  gbound_diel,gylmg_diel,icg_mpi,ikpt,isp,istwf_k,kg_diel,kg_k_mpi,&
&  lmax_diel,mband,mcg,mgfftdiel,mpi_enreg,&
&  natom,nband_k,ndiel4,ndiel5,ndiel6,neglect_pawhat,ngfftdiel,nkpt,&
&  npwdiel,npw_k_mpi,nspden,nspden_eff,nspinor,nsppol,ntypat,occ,occopt,occ_deavg,&
&  pawang,pawtab,ph3d_diel,rhoextrap,sumdocc,&
&  susmat,typat,ucvol,usepaw,wtk)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi

 use m_pawang,  only : pawang_type
 use m_pawtab,  only : pawtab_type
 use m_pawcprj, only : pawcprj_type
 use m_mpinfo,  only : destroy_mpi_enreg

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'susk'
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_52_fft_mpi_noabirule
 use interfaces_53_ffts
 use interfaces_65_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!This type is defined in defs_mpi
!scalars
 integer,intent(in) :: bdtot_index,extrap,ikpt,isp,istwf_k,lmax_diel,mband,mcg
 integer,intent(in) :: mgfftdiel,natom,nband_k,ndiel4,ndiel5,ndiel6,neglect_pawhat
 integer,intent(in) :: nkpt,npwdiel,nspden,nspden_eff,nspinor,nsppol
 integer,intent(in) :: ntypat,occopt,usepaw
 integer,intent(in),target :: icg_mpi,npw_k_mpi
 real(dp),intent(in) :: ucvol
 real(dp),intent(inout) :: sumdocc
 type(MPI_type),intent(in) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: atindx(natom),gbound_diel(2*mgfftdiel+8,2)
 integer,intent(in) :: kg_diel(3,npwdiel),ngfftdiel(18),typat(natom)
 integer,intent(in),target :: kg_k_mpi(3,npw_k_mpi)
 integer,intent(inout) :: gbound(2*mgfftdiel+8,2)
 integer,pointer :: kg_k(:,:)
 real(dp),intent(in) :: doccde(mband*nkpt*nsppol),eigen(mband*nkpt*nsppol)
 real(dp),intent(in) :: gylmg_diel(npwdiel,lmax_diel**2,ntypat*usepaw)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),occ_deavg(mband)
 real(dp),intent(in) :: ph3d_diel(2,npwdiel,natom*usepaw),wtk(nkpt)
 real(dp),intent(in),target :: cg_mpi(2,mcg)
 real(dp),intent(inout) :: drhode(2,npwdiel,nspden_eff)
 real(dp),intent(inout) :: rhoextrap(ndiel4,ndiel5,ndiel6,nspinor)
 real(dp),intent(inout) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 type(pawcprj_type) :: cprj_k(natom,nspinor*nband_k*usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
! real(dp), allocatable :: cg_disk(:,:)
!Local variables for MPI
!scalars
 integer :: blocksize,i1,i2,i3,iband,iband_loc,ibd1,ibd2,ibdblock,ier
 integer :: iproc,iproc_fft,ipw,ipw1,ipw2,isp1,isp2,ispinor,iwf,jsp,me_bandfft
 integer :: nbdblock,ndatarecv,ndiel1,ndiel2,ndiel3
 integer :: sizemax_per_proc,spaceComm,testocc,tim_fourwf
 integer,pointer :: icg,npw_k
 integer,target :: icg_loc=0,npw_k_loc,npw_tot
 real(dp) :: eigdiff,occdiff,tolocc,weight,wght1,wght2
 type(MPI_type) :: mpi_enreg_diel
!arrays
 integer,allocatable :: band_loc(:),kg_k_gather(:,:),npw_per_proc(:),rdispls(:)
 integer,allocatable :: rdispls_all(:),rdisplsloc(:),recvcounts(:)
 integer,allocatable :: recvcountsloc(:),sdispls(:),sdisplsloc(:),sendcounts(:)
 integer,allocatable :: sendcountsloc(:),susmat_mpi(:,:,:)
 integer,allocatable,target :: kg_k_gather_all(:,:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: cwavef(:,:),cwavef_alltoall(:,:)
 real(dp),allocatable :: cwavef_alltoall_gather(:,:),dummy(:,:),rhoaug(:,:,:)
 real(dp),allocatable :: wfprod(:,:),wfraug(:,:,:,:),wfrspa(:,:,:,:,:,:)
 real(dp),allocatable,target :: cg_local(:,:)
 real(dp),pointer :: cg(:,:)
 logical,allocatable :: treat_band(:)

! *************************************************************************

!DEBUG
!write(std_out,*)' susk : enter '
!if(.true.)stop
!ENDDEBUG

 call timab(750,1,tsec)
 call timab(751,1,tsec)

 ndiel1=ngfftdiel(1) ; ndiel2=ngfftdiel(2) ; ndiel3=ngfftdiel(3)

!The dielectric stuff is performed in sequential mode.
!Set mpi_enreg_diel accordingly
 call initmpi_seq(mpi_enreg_diel)
 call init_distribfft_seq(mpi_enreg_diel%distribfft,'c',ngfftdiel(2),ngfftdiel(3),'all')
 me_bandfft=xmpi_comm_rank(mpi_enreg%comm_bandfft)

 testocc=1
!DEBUG
!write(std_out,*)' susk : set testocc to 0 '
!testocc=0
!write(std_out,*)' susk : set extrap to 0 '
!extrap=0
!ENDDEBUG

!Allocations, initializations
 ABI_ALLOCATE(rhoaug,(ndiel4,ndiel5,ndiel6))
 ABI_ALLOCATE(wfraug,(2,ndiel4,ndiel5,ndiel6))
 ABI_ALLOCATE(wfprod,(2,npwdiel))
 ABI_ALLOCATE(wfrspa,(2,ndiel4,ndiel5,ndiel6,nspinor,mband))
 ABI_ALLOCATE(dummy,(2,1))
 wfrspa(:,:,:,:,:,:)=zero
 ABI_ALLOCATE(treat_band,(nband_k))
 treat_band(:)=.true.
 isp1=isp;isp2=isp
 if (nspden_eff==2.and.nspinor==2) isp2=isp+1

!BAND-FFT parallelism
 if (mpi_enreg%paral_kgb==1) then
   treat_band(:)=.false.
!  We gather the wavefunctions treated by this proc in cg_local
   spaceComm=mpi_enreg%comm_band
   blocksize=mpi_enreg%nproc_band
   nbdblock=nband_k/blocksize
   ABI_ALLOCATE(sdispls,(blocksize))
   ABI_ALLOCATE(sdisplsloc,(blocksize))
   ABI_ALLOCATE(sendcounts,(blocksize))
   ABI_ALLOCATE(sendcountsloc,(blocksize))
   ABI_ALLOCATE(rdispls,(blocksize))
   ABI_ALLOCATE(rdisplsloc,(blocksize))
   ABI_ALLOCATE(recvcounts,(blocksize))
   ABI_ALLOCATE(recvcountsloc,(blocksize))
!  First gather the kg_k in kg_k_gather_all
   npw_k_loc=npw_k_mpi
   call xmpi_allgather(npw_k_loc,recvcounts,spaceComm,ier)
   rdispls(1)=0
   do iproc=2,blocksize
     rdispls(iproc)=rdispls(iproc-1)+recvcounts(iproc-1)
   end do
   ndatarecv=rdispls(blocksize)+recvcounts(blocksize)
   ABI_ALLOCATE(kg_k_gather,(3,ndatarecv))
   recvcountsloc(:)=recvcounts(:)*3
   rdisplsloc(:)=rdispls(:)*3
   call xmpi_allgatherv(kg_k_mpi,3*npw_k_loc,kg_k_gather,recvcountsloc(:),rdisplsloc,spaceComm,ier)
   ABI_ALLOCATE(npw_per_proc,(mpi_enreg%nproc_fft))
   ABI_ALLOCATE(rdispls_all,(mpi_enreg%nproc_fft))
   spaceComm=mpi_enreg%comm_fft
   call xmpi_allgather(ndatarecv,npw_per_proc,spaceComm,ier)
   rdispls_all(1)=0
   do iproc=2,mpi_enreg%nproc_fft
     rdispls_all(iproc)=rdispls_all(iproc-1)+npw_per_proc(iproc-1)
   end do
   npw_tot=rdispls_all(mpi_enreg%nproc_fft)+npw_per_proc(mpi_enreg%nproc_fft)
   ABI_ALLOCATE(kg_k_gather_all,(3,npw_tot))
   call xmpi_allgatherv(kg_k_gather,3*ndatarecv,kg_k_gather_all,3*npw_per_proc(:),3*rdispls_all,spaceComm,ier)
!  At this point kg_k_gather_all contains all the kg
   if(allocated(cwavef))  then
     ABI_DEALLOCATE(cwavef)
   end if
   ABI_ALLOCATE(cwavef,(2,npw_k_loc*nspinor*blocksize))
   sizemax_per_proc=nband_k/(mpi_enreg%nproc_band*mpi_enreg%nproc_fft)+1
   ABI_ALLOCATE(band_loc,(sizemax_per_proc))
   ABI_ALLOCATE(cg_local,(2,sizemax_per_proc*npw_tot*nspinor))
   iband_loc=0
   do ibdblock=1,nbdblock
     cwavef(:,1:npw_k_loc*nspinor*blocksize)=&
&     cg_mpi(:,1+(ibdblock-1)*npw_k_loc*nspinor*blocksize+icg_mpi:ibdblock*npw_k_loc*nspinor*blocksize+icg_mpi)
     sendcounts(:)=npw_k_loc
     do iproc=1,blocksize
       sdispls(iproc)=(iproc-1)*npw_k_loc
     end do
     ABI_ALLOCATE(cwavef_alltoall,(2,ndatarecv*nspinor))
     recvcountsloc(:)=recvcounts(:)*2*nspinor
     rdisplsloc(:)=rdispls(:)*2*nspinor
     sendcountsloc(:)=sendcounts(:)*2*nspinor
     sdisplsloc(:)=sdispls(:)*2*nspinor
     call timab(547,1,tsec)
     spaceComm=mpi_enreg%comm_band
     call xmpi_alltoallv(cwavef,sendcountsloc,sdisplsloc,cwavef_alltoall,recvcountsloc,rdisplsloc,spaceComm,ier)
     call timab(547,2,tsec)
     ABI_ALLOCATE(cwavef_alltoall_gather,(2,npw_tot*nspinor))
     blocksize=mpi_enreg%nproc_band
     spaceComm=mpi_enreg%comm_fft
     call xmpi_allgatherv(cwavef_alltoall,2*nspinor*ndatarecv,cwavef_alltoall_gather,&
&     2*nspinor*npw_per_proc,2*nspinor*rdispls_all,spaceComm,ier)
     iproc_fft=modulo(ibdblock-1,mpi_enreg%nproc_fft)
     if(mpi_enreg%me_fft==iproc_fft) then !All nproc_band procs of index me_fft will treat these bands
       iband_loc=iband_loc+1
       iband=1+mpi_enreg%me_band+mpi_enreg%nproc_band*mpi_enreg%me_fft+(iband_loc-1)*mpi_enreg%nproc_fft*mpi_enreg%nproc_band
       treat_band(iband)=.true.
       band_loc(iband_loc)=iband
       cg_local(:,1+(iband_loc-1)*npw_tot*nspinor:iband_loc*npw_tot*nspinor)=cwavef_alltoall_gather(:,1:npw_tot*nspinor)
     end if
     ABI_DEALLOCATE(cwavef_alltoall_gather)
     ABI_DEALLOCATE(cwavef_alltoall)
   end do
!  On exit:
!  npw_tot will be npw
!  kg_k_gather_all will be kg_k
!  cg_local will be cg
!  icg will be zero
   npw_k=>npw_tot
   kg_k=>kg_k_gather_all(:,:)
   cg=>cg_local(:,:)
   icg=>icg_loc
   call sphereboundary(gbound,istwf_k,kg_k,mgfftdiel,npw_k)
   ABI_DEALLOCATE(npw_per_proc)
   ABI_DEALLOCATE(rdispls_all)
   ABI_DEALLOCATE(sendcounts)
   ABI_DEALLOCATE(recvcounts)
   ABI_DEALLOCATE(sdispls)
   ABI_DEALLOCATE(rdispls)
   ABI_DEALLOCATE(sendcountsloc)
   ABI_DEALLOCATE(sdisplsloc)
   ABI_DEALLOCATE(recvcountsloc)
   ABI_DEALLOCATE(rdisplsloc)
   ABI_DEALLOCATE(kg_k_gather)
   ABI_DEALLOCATE(cwavef)
!  Because they will be summed over all procs, and arrive on input, rescale drhode and rhoextrap
   if(occopt>=3)drhode(:,:,isp1:isp2)=drhode(:,:,isp1:isp2)/real(mpi_enreg%nproc_fft*mpi_enreg%nproc_band,dp)
   if(extrap==1)rhoextrap(:,:,:,:)=rhoextrap(:,:,:,:)/real(mpi_enreg%nproc_fft*mpi_enreg%nproc_band,dp)
   do i1=isp1,isp2
     susmat(:,:,i1,:,i1)=susmat(:,:,i1,:,i1)/real(mpi_enreg%nproc_fft*mpi_enreg%nproc_band,dp)
   end do

!  No BAND-FFT parallelism
 else ! use argument variables
   cg=>cg_mpi
   kg_k=>kg_k_mpi
   npw_k=>npw_k_mpi
   icg=>icg_mpi
 end if
 iband_loc=0

 call timab(751,2,tsec)
 call timab(752,1,tsec)

!Loop over bands to fft and store Fourier transform of wavefunction
 ABI_ALLOCATE(cwavef,(2,npw_k))
 do iband=1,nband_k
   if(.not. treat_band(iband))  cycle ! I am not treating this band (only for the parallel case)
   iband_loc=iband_loc+1

!  Loop on spinorial components
   do ispinor=1,nspinor
     iwf=(ispinor-1)*npw_k+(iband_loc-1)*npw_k*nspinor+icg
     jsp=isp+ispinor-1;if (nspden_eff==1) jsp=isp

!    Obtain Fourier transform in fft box
     tim_fourwf=8
     cwavef(:,1:npw_k)=cg(:,1+iwf:npw_k+iwf)
     call fourwf(1,rhoaug,cwavef,dummy,wfraug,gbound,gbound,&
&     istwf_k,kg_k,kg_k,mgfftdiel,mpi_enreg_diel,1,ngfftdiel,npw_k,1,ndiel4,ndiel5,ndiel6,&
&     0,mpi_enreg_diel%paral_kgb,tim_fourwf,weight,weight)

     wfrspa(:,:,:,:,ispinor,iband)=wfraug(:,:,:,:)

     if( (occopt>=3 .and. testocc==1) .or. extrap==1 )then
!      In the case of metallic occupation, or if the extrapolation
!      over higher bands is included, must compute the
!      Fourier transform of the density of each band, then
!      generate the part of the susceptibility matrix due
!      varying occupation numbers.

       weight=-two*occ_deavg(iband)*wtk(ikpt)/ucvol
       do i3=1,ndiel3
         do i2=1,ndiel2
           do i1=1,ndiel1
             wfraug(1,i1,i2,i3)=wfraug(1,i1,i2,i3)**2+wfraug(2,i1,i2,i3)**2
             wfraug(2,i1,i2,i3)=zero
           end do
         end do
!        If extrapolation, accumulate density in real space
         if(extrap==1.and.usepaw==0)then
           do i2=1,ndiel2
             do i1=1,ndiel1
               rhoextrap(i1,i2,i3,ispinor)=rhoextrap(i1,i2,i3,ispinor)+weight*wfraug(1,i1,i2,i3)
             end do
           end do
         end if
       end do

!      In case of PAW, add compensation charge contribution
       if (usepaw==1.and.extrap==1.and.neglect_pawhat==0) then
         call pawsushat(atindx,cprj_k,gbound_diel,gylmg_diel,iband,iband,ispinor,ispinor,1,kg_diel,&
&         lmax_diel,mgfftdiel,natom,nband_k,ndiel4,ndiel5,ndiel6,&
&         ngfftdiel,npwdiel,nspinor,ntypat,1,&
&         pawang,pawtab,ph3d_diel,typat,dummy,wfraug,&
&         mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
         rhoextrap(:,:,:,ispinor)=rhoextrap(:,:,:,ispinor)+weight*wfraug(1,:,:,:)
       end if

!      Performs the Fourier Transform of the density of the band,
!      and store it in wfprod
       tim_fourwf=9
       call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&       1,kg_diel,kg_diel,mgfftdiel,mpi_enreg_diel,1,ngfftdiel,1,npwdiel,&
&       ndiel4,ndiel5,ndiel6,3,mpi_enreg_diel%paral_kgb,tim_fourwf,weight,weight)
!      In case of PAW, add compensation charge contribution if not already done
       if (usepaw==1.and.extrap==0.and.neglect_pawhat==0) then
         call pawsushat(atindx,cprj_k,gbound_diel,gylmg_diel,ibd1,ibd2,ispinor,ispinor,1,kg_diel,&
&         lmax_diel,mgfftdiel,natom,nband_k,ndiel4,ndiel5,ndiel6,&
&         ngfftdiel,npwdiel,nspinor,ntypat,0,&
&         pawang,pawtab,ph3d_diel,typat,wfprod,dummy,&
&         mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
       end if

!      Perform now the summation of terms related to direct change of eigenvalues
!      or extrapolation over higher bands
       wght1=zero ; wght2=zero
       if(occopt>=3 .and. testocc==1) wght1=doccde(iband+bdtot_index)*wtk(ikpt)/ucvol
       if(extrap==1) wght2=two*occ_deavg(iband)*wtk(ikpt)/ucvol
       weight=wght1+wght2

       if (abs(weight)>tol12) then
         do ipw2=1,npwdiel
!          Only fills lower half of the matrix (here, the susceptibility matrix)
!          Note that wfprod of the first index must behave like a density,
!          so that it is used as generated by fourwf, while wfprod of the
!          second index will be implicitely used to make a scalar product
!          with a potential change, meaning that its complex conjugate must be
!          used. This explains the following signs...
           do ipw1=ipw2,npwdiel
             susmat(1,ipw1,jsp,ipw2,jsp)=susmat(1,ipw1,jsp,ipw2,jsp)+&
&             weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
             susmat(2,ipw1,jsp,ipw2,jsp)=susmat(2,ipw1,jsp,ipw2,jsp)+&
&             weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
           end do
         end do
       end if

       if( occopt>=3 .and. testocc==1 .and. abs(wght1)>tol12) then
!        Accumulate product of band densities by their doccde, for the
!        computation of the effect of change of Fermi level.
         do ipw=1,npwdiel
           drhode(1,ipw,jsp)=drhode(1,ipw,jsp)+wfprod(1,ipw)*wght1
           drhode(2,ipw,jsp)=drhode(2,ipw,jsp)+wfprod(2,ipw)*wght1
         end do
!        Also accumulate weighted sum of doccde
         sumdocc=sumdocc+wght1
       end if

!      End condition of metallic occupancies or extrapolation
     end if

!    End loop on spinorial components
   end do
!  End loop on iband
 end do

 call timab(752,2,tsec)
 call timab(753,1,tsec)

 ABI_DEALLOCATE(cwavef)

!Stuff for parallelism (bands-FFT)
 if(mpi_enreg%paral_kgb==1) then
   call xmpi_sum(wfrspa,mpi_enreg%comm_bandfft,ier)
   if(occopt>=3) then
     call xmpi_sum(drhode(:,:,isp1:isp2),mpi_enreg%comm_bandfft,ier)
   end if
   if(extrap==1) then
     call xmpi_sum(rhoextrap,mpi_enreg%comm_bandfft,ier)
   end if
   if(occopt>=3) then
     call xmpi_sum(sumdocc,mpi_enreg%comm_bandfft,ier)
   end if
   ABI_ALLOCATE(susmat_mpi,(2,npwdiel,npwdiel))
   do i1=isp1,isp2
     susmat_mpi(:,:,:)=susmat(:,:,i1,:,i1)
     call xmpi_sum(susmat_mpi,mpi_enreg%comm_bandfft,ier)
     susmat(:,:,i1,:,i1)=susmat_mpi(:,:,:)/real(mpi_enreg%nproc_fft*mpi_enreg%nproc_band,dp)
   end do
   ABI_DEALLOCATE(susmat_mpi)
 end if
 call timab(753,2,tsec)

!-- Wavefunctions have been generated in real space ------------------------
!-- Now, compute product of wavefunctions for different bands --------------
 call timab(754,1,tsec)
!if (occopt<3) then
 tolocc=1.0d-3
!else
!tolocc=1.0d-8
!end if
 iproc=-1

 if(nband_k>1)then
   do ibd1=1,nband_k-1
     do ibd2=ibd1+1,nband_k
       iproc=iproc+1
       if(modulo(iproc,mpi_enreg%nproc_fft*mpi_enreg%nproc_band) /= me_bandfft) cycle
!      If the occupation numbers are sufficiently different, or
!      if extrapolation is used and the corresponding factor is not zero,
!      then there is a contribution
       occdiff=occ(ibd1+bdtot_index)-occ(ibd2+bdtot_index)
       if( abs(occdiff)>tolocc      .or. &
&       ( extrap==1 .and.            &
&       ( abs(occ_deavg(ibd1)) + abs(occ_deavg(ibd2)) ) >tolocc ) &
&       ) then

         eigdiff=eigen(ibd1+bdtot_index)-eigen(ibd2+bdtot_index)
!        DEBUG
!        write(std_out,*)' susk : contribution from bands',ibd1,ibd2
!        write(std_out,*)'   occ diff =',occdiff
!        write(std_out,*)'   eig diff =',eigdiff
!        ENDDEBUG

!        Loop on spinorial components
         do ispinor=1,nspinor
           jsp=isp+ispinor-1;if (nspden_eff==1) jsp=isp

!          Store the contribution in wfraug
           do i3=1,ndiel3
             do i2=1,ndiel2
               do i1=1,ndiel1
                 wfraug(1,i1,i2,i3)=wfrspa(1,i1,i2,i3,ispinor,ibd1)*wfrspa(1,i1,i2,i3,ispinor,ibd2)&
&                 +wfrspa(2,i1,i2,i3,ispinor,ibd1)*wfrspa(2,i1,i2,i3,ispinor,ibd2)
                 wfraug(2,i1,i2,i3)=wfrspa(2,i1,i2,i3,ispinor,ibd1)*wfrspa(1,i1,i2,i3,ispinor,ibd2)&
&                 -wfrspa(1,i1,i2,i3,ispinor,ibd1)*wfrspa(2,i1,i2,i3,ispinor,ibd2)
               end do
             end do
           end do

!          Performs the Fourier Transform of the product, and store it in wfprod
           tim_fourwf=19
           call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&           1,kg_diel,kg_diel,mgfftdiel,mpi_enreg_diel,1,ngfftdiel,1,npwdiel,&
&           ndiel4,ndiel5,ndiel6,3,mpi_enreg_diel%paral_kgb,tim_fourwf,weight,weight)

!          In case of PAW, add compensation charge contribution
           if (usepaw==1.and.neglect_pawhat==0) then
             call pawsushat(atindx,cprj_k,gbound_diel,gylmg_diel,ibd1,ibd2,ispinor,ispinor,1,kg_diel,&
&             lmax_diel,mgfftdiel,natom,nband_k,ndiel4,ndiel5,ndiel6,&
&             ngfftdiel,npwdiel,nspinor,ntypat,0,&
&             pawang,pawtab,ph3d_diel,typat,wfprod,dummy,&
&             mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
           end if

!          Perform now the summation
           wght1=zero ; wght2=zero
           if(abs(occdiff)>tolocc) wght1= occdiff/eigdiff * two*wtk(ikpt)/ucvol
           if(extrap==1) wght2=(occ_deavg(ibd1)+occ_deavg(ibd2)) * two*wtk(ikpt)/ucvol
           weight=wght1+wght2

!          DEBUG
!          write(std_out,*)' weight =',weight
!          norm=zero
!          do ipw=1,npwdiel
!          norm=norm+wfprod(1,ipw)**2+wfprod(2,ipw)**2
!          end do
!          write(std_out,*)' norm in reciprocal space  =',norm
!          ENDDEBUG

           if (abs(weight)>tol12) then
             do ipw2=1,npwdiel
!              Only fills lower half of the matrix (here, the susceptibility matrix)
!              Note that wfprod of the first index must behave like a density,
!              so that it is used as generated by fourwf, while wfprod of the
!              second index will be implicitely used to make a scalar product
!              with a potential change, meaning that its complex conjugate must be
!              used. This explains the following signs...
               do ipw1=ipw2,npwdiel
                 susmat(1,ipw1,jsp,ipw2,jsp)=susmat(1,ipw1,jsp,ipw2,jsp)+&
&                 weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
                 susmat(2,ipw1,jsp,ipw2,jsp)=susmat(2,ipw1,jsp,ipw2,jsp)+&
&                 weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
               end do
             end do
           end if

!          End loop on spinorial components
         end do
!        End condition of different occupation numbers or extrapolation
       end if
!      End internal loop over bands
     end do
!    End external loop over bands
   end do
!  End condition of having more than one band
 end if

 call timab(754,2,tsec)
 call timab(755,1,tsec)

 if(mpi_enreg%paral_kgb==1) then
   ABI_ALLOCATE(susmat_mpi,(2,npwdiel,npwdiel))
   do i1=isp1,isp2
     susmat_mpi(:,:,:)=susmat(:,:,i1,:,i1)
     call xmpi_sum(susmat_mpi,mpi_enreg%comm_bandfft,ier)
     susmat(:,:,i1,:,i1)=susmat_mpi(:,:,:)
   end do
   ABI_DEALLOCATE(susmat_mpi)
   ABI_DEALLOCATE(band_loc)
   ABI_DEALLOCATE(treat_band)
   ABI_DEALLOCATE(cg_local)
   ABI_DEALLOCATE(kg_k_gather_all)
 end if

 call destroy_mpi_enreg(mpi_enreg_diel)
 ABI_DEALLOCATE(dummy)
 ABI_DEALLOCATE(rhoaug)
 ABI_DEALLOCATE(wfprod)
 ABI_DEALLOCATE(wfraug)
 ABI_DEALLOCATE(wfrspa)

 call timab(755,2,tsec)
 call timab(750,2,tsec)

end subroutine susk
!!***
