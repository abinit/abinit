!{\src2tex{textfont=tt}}
!!****f* ABINIT/suskmm
!! NAME
!! suskmm
!!
!! FUNCTION
!! Compute the contribution of one k point to the susceptibility matrix
!! from input wavefunctions, band occupations, and k point wts.
!! Include the usual sum-over-state terms, but also the
!! corrections due to the change of the Fermi level in the metallic
!! case, as well as implicit sum over higher lying conduction
!! states, thanks to the closure relation (referred to as an extrapolation).
!!
!! This routine is similar to susk, but use blocking on wavefunctions
!! to decrease memory requirements, at the expense of CPU time.
!!
!! NOTES
!! There is still room for optimization !!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms
!!  bdtot_index=index for the number of the band
!!  cg(2,mcg)=wf in G space
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
!!  gylmg_diel(npwdiel,lmax_diel**2,ntypat*usepaw)= -PAW only- Fourier transform of g_l(r).Y_ml(r) shape functions
!!                                                   for dielectric matrix
!!  icg=index for cg
!!  ikpt=number of the k point
!!  isp=number of the current spin
!!  istwf_k=input option parameter that describes the storage of wfs
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  kg_k(3,npw)=coordinates of planewaves in basis sphere.
!!  lmax_diel=1+max. value of l angular momentum used for dielectric matrix
!!  mband=maximum number of bands
!!  mcg=dimension of cg
!!  mgfftdiel=maximum size of 1D FFTs, for the computation of
!!     the dielectric matrix
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  nband_k=number of bands at this k point for that spin polarization
!!  ndiel4,ndiel5,ndiel6= FFT dimensions, modified to avoid cache trashing
!!  neglect_pawhat=1 if PAW contribution from hat density (compensation charge)
!!                 has to be neglected (to be used when only an estimation of
!!                 suscep. matrix has to be evaluated, i.e. for SCF precondictioning)
!!  ngfftdiel(18)=contain all needed information about 3D FFT, for dielectric matrix,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
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
!!  drhode(2,npwdiel,nspden_eff)=weighted density, needed to compute the
!!   effect of change of fermi energy
!!  rhoextrap(ndiel4,ndiel5,ndiel6,nspinor)=density-like array, needed for the
!!   extrapolation procedure.
!!  sumdocc=sum of weighted occupation numbers, needed to compute the
!!   effect of change of fermi energy
!!  susmat(2,npwdiel,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!
!! PARENTS
!!      suscep_stat
!!
!! CHILDREN
!!      destroy_mpi_enreg,fourwf,init_distribfft_seq,initmpi_seq,pawsushat
!!      timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine suskmm(atindx,bdtot_index,cg,cprj_k,doccde,drhode,eigen,extrap,gbound,&
&  gbound_diel,gylmg_diel,icg,ikpt,isp,istwf_k,kg_diel,kg_k,&
&  lmax_diel,mband,mcg,mgfftdiel,mpi_enreg,&
&  natom,nband_k,ndiel4,ndiel5,ndiel6,neglect_pawhat,ngfftdiel,nkpt,&
&  npwdiel,npw_k,nspden,nspden_eff,nspinor,nsppol,ntypat,occ,occopt,occ_deavg,paral_kgb,&
&  pawang,pawtab,ph3d_diel,rhoextrap,sumdocc,&
&  susmat,typat,ucvol,usepaw,wtk)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_mpinfo,  only : destroy_mpi_enreg

 use m_pawang,  only : pawang_type
 use m_pawtab,  only : pawtab_type
 use m_pawcprj, only : pawcprj_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'suskmm'
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_65_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!This type is defined in defs_mpi
!scalars
 integer,intent(in) :: bdtot_index,extrap,icg,ikpt,isp,istwf_k,lmax_diel,mband,mcg
 integer,intent(in) :: mgfftdiel,natom,nband_k,ndiel4,ndiel5,ndiel6,neglect_pawhat
 integer,intent(in) :: nkpt,npw_k,npwdiel,nspden,nspden_eff,nspinor
 integer,intent(in) :: nsppol,ntypat,occopt,paral_kgb,usepaw
 real(dp),intent(in) :: ucvol
 real(dp),intent(inout) :: sumdocc
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: atindx(natom),gbound(2*mgfftdiel+8,2)
 integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
 integer,intent(in) :: kg_diel(3,npwdiel),kg_k(3,npw_k),ngfftdiel(18)
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: cg(2,mcg),doccde(mband*nkpt*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
 real(dp),intent(in) :: gylmg_diel(npwdiel,lmax_diel**2,ntypat*usepaw)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),occ_deavg(mband)
 real(dp),intent(in) :: ph3d_diel(2,npwdiel,natom*usepaw),wtk(nkpt)
 real(dp),intent(inout) :: drhode(2,npwdiel,nspden_eff)
 real(dp),intent(inout) :: rhoextrap(ndiel4,ndiel5,ndiel6,nspinor)
 real(dp),intent(inout) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 type(pawcprj_type) :: cprj_k(natom,nspinor*nband_k*usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: comm_fft,i1,i2,i3,iband,iband_shift,iband_shift2,ibd1,ibd2,ibdshft1,ibdshft2
 integer :: iblk1,iblk2,ipw,ipw1,ipw2,ispinor,iwf,jsp,mblk
 integer :: nblk,nbnd_current,nbnd_in_blk,nbnd_in_blk1,ndiel1,ndiel2,ndiel3
 integer :: testocc,tim_fourwf
 real(dp) :: eigdiff,occdiff,tolocc,weight,wght1,wght2
 character(len=500) :: message
 type(MPI_type) :: mpi_enreg_diel
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: cwavef(:,:),dummy(:,:),rhoaug(:,:,:),wfprod(:,:)
 real(dp),allocatable :: wfraug(:,:,:,:),wfrspa1(:,:,:,:,:,:)
 real(dp),allocatable :: wfrspa2(:,:,:,:,:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' suskmm : enter '
!if(.true.)stop
!ENDDEBUG

 call timab(760,1,tsec)
 call timab(761,1,tsec)

!Allocations, initializations
 ndiel1=ngfftdiel(1) ; ndiel2=ngfftdiel(2) ; ndiel3=ngfftdiel(3)
 testocc=1
 ABI_ALLOCATE(rhoaug,(ndiel4,ndiel5,ndiel6))
 ABI_ALLOCATE(wfraug,(2,ndiel4,ndiel5,ndiel6))
 ABI_ALLOCATE(wfprod,(2,npwdiel))
 ABI_ALLOCATE(dummy,(2,1))

!The dielectric stuff is performed in sequential mode.
!Set mpi_enreg_diel accordingly
 call initmpi_seq(mpi_enreg_diel)
 call init_distribfft_seq(mpi_enreg_diel%distribfft,'c',ngfftdiel(2),ngfftdiel(3),'all')

 comm_fft=mpi_enreg%comm_fft
 
!Prepare the blocking : compute the number of blocks,
!the number of bands in each normal block,
!and the number in the first one, usually smaller.

!Consider that if the number of bands is large, there are at most 8 blocks
 nbnd_in_blk=0
 if(nband_k>=48)then
   mblk=8
   nbnd_in_blk=(nband_k-1)/mblk+1
!  If the number of bands is medium, place 6 bands per block
 else if(nband_k>=12)then
   nbnd_in_blk=6
!  Otherwise, must have at least 2 blocks
 else if(nband_k>=2)then
   mblk=2
   nbnd_in_blk=(nband_k-1)/mblk+1
 else
   write(message, '(a,a,a,i2,a,a,a)')&
&   '  The number of bands must be larger or equal to 2, in suskmm.',ch10,&
&   '  It is equal to ',nband_k,'.',ch10,&
&   '  Action : choose another preconditioner.'
   MSG_ERROR(message)
 end if

!Compute the effective number of blocks, and the number of bands in
!the first block.
 nblk=(nband_k-1)/nbnd_in_blk+1
 nbnd_in_blk1=nband_k-(nblk-1)*nbnd_in_blk

!DEBUG
!write(std_out,*)' suskmm : nband_k,nblk,nbnd_in_blk,nbnd_in_blk1 '
!write(std_out,*)nband_k,nblk,nbnd_in_blk,nbnd_in_blk1
!stop
!ENDDEBUG

!wfrspa1 will contain the wavefunctions of the slow sampling (iblk1)
 ABI_ALLOCATE(wfrspa1,(2,ndiel4,ndiel5,ndiel6,nspinor,nbnd_in_blk))
!wfrspa2 will contain the wavefunctions of the rapid sampling (iblk2)
 ABI_ALLOCATE(wfrspa2,(2,ndiel4,ndiel5,ndiel6,nspinor,nbnd_in_blk))

 ABI_ALLOCATE(cwavef,(2,npw_k))

 call timab(761,2,tsec)

!First loop over blocks
 do iblk1=1,nblk

   call timab(762,1,tsec)

!  Initialisation
   if(iblk1==1)then

     nbnd_current=nbnd_in_blk1
     iband_shift=0
!    Loop over bands to fft and store Fourier transform of wavefunction
     do iband=1,nbnd_current
!      Loop on spinorial components
       do ispinor=1,nspinor
         iwf=(ispinor-1)*npw_k+(iband-1)*npw_k*nspinor+icg
!        Obtain Fourier transform in fft box
         tim_fourwf=21
         cwavef(:,1:npw_k)=cg(:,1+iwf:npw_k+iwf)
         call fourwf(1,rhoaug,cwavef,dummy,wfraug,gbound,gbound,&
&         istwf_k,kg_k,kg_k,mgfftdiel,mpi_enreg_diel,1,ngfftdiel,npw_k,1,ndiel4,ndiel5,ndiel6,&
&         0,paral_kgb,tim_fourwf,weight,weight)
         wfrspa1(:,:,:,:,ispinor,iband)=wfraug(:,:,:,:)
       end do
     end do

   else

!    The Fourier transform of wavefunctions have already been obtained
     nbnd_current=nbnd_in_blk
     iband_shift=nbnd_in_blk1+(iblk1-2)*nbnd_in_blk

   end if

!  Loop over bands of this block, to generate band-diagonal
   do iband=1,nbnd_current

!    Loop on spinorial components
     do ispinor=1,nspinor
       jsp=isp+ispinor-1;if (nspden_eff==1) jsp=isp

       if( (occopt>=3 .and. testocc==1) .or. extrap==1 )then
!        In the case of metallic occupation, or if the extrapolation
!        over higher bands is included, must compute the
!        Fourier transform of the density of each band, then
!        generate the part of the susceptibility matrix due
!        varying occupation numbers.
         weight=-two*occ_deavg(iband+iband_shift)*wtk(ikpt)/ucvol
         do i3=1,ndiel3
           do i2=1,ndiel2
             do i1=1,ndiel1
               wfraug(1,i1,i2,i3)=wfrspa1(1,i1,i2,i3,ispinor,iband)**2&
&               +wfrspa1(2,i1,i2,i3,ispinor,iband)**2
               wfraug(2,i1,i2,i3)=zero
             end do
           end do
!          If extrapolation, accumulate density in real space
           if(extrap==1.and.usepaw==0)then
             do i2=1,ndiel2
               do i1=1,ndiel1
                 rhoextrap(i1,i2,i3,ispinor)=rhoextrap(i1,i2,i3,ispinor)+weight*wfraug(1,i1,i2,i3)
               end do
             end do
           end if
         end do

!        In case of PAW, add compensation charge contribution
         if (usepaw==1.and.extrap==1.and.neglect_pawhat==0) then
           call pawsushat(atindx,cprj_k,gbound_diel,gylmg_diel,iband,iband,ispinor,ispinor,1,kg_diel,&
&           lmax_diel,mgfftdiel,natom,nband_k,ndiel4,ndiel5,ndiel6,&
&           ngfftdiel,npwdiel,nspinor,ntypat,1,&
&           pawang,pawtab,ph3d_diel,typat,dummy,wfraug,&
&           mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
           rhoextrap(:,:,:,ispinor)=rhoextrap(:,:,:,ispinor)+weight*wfraug(1,:,:,:)
         end if

!        Performs the Fourier Transform of the density of the band,
!        and store it in wfprod
         tim_fourwf=31
         call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&         1,kg_diel,kg_diel,&
&         mgfftdiel,mpi_enreg_diel,1,ngfftdiel,1,npwdiel,ndiel4,ndiel5,ndiel6,3,paral_kgb,tim_fourwf,weight,weight)
!        In case of PAW, add compensation charge contribution if not already done
         if (usepaw==1.and.extrap==0.and.neglect_pawhat==0) then
           call pawsushat(atindx,cprj_k,gbound_diel,gylmg_diel,iband,iband,1,1,1,kg_diel,&
&           lmax_diel,mgfftdiel,natom,nband_k,ndiel4,ndiel5,ndiel6,&
&           ngfftdiel,npwdiel,nspinor,ntypat,0,&
&           pawang,pawtab,ph3d_diel,typat,wfprod,dummy,&
&           mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
         end if

!        Perform now the summation of terms related to direct change of eigenvalues
!        or extrapolation over higher bands
         wght1=zero ; wght2=zero
         if(occopt>=3 .and. testocc==1) wght1=doccde(iband+iband_shift+bdtot_index)*wtk(ikpt)/ucvol
         if(extrap==1) wght2=two*occ_deavg(iband+iband_shift)*wtk(ikpt)/ucvol
         weight=wght1+wght2

         if (abs(weight)>tol12) then
           do ipw2=1,npwdiel
!            Only fills lower half of the matrix (here, the susceptibility matrix)
!            Note that wfprod of the first index must behave like a density,
!            so that it is used as generated by fourwf, while wfprod of the
!            second index will be implicitely used to make a scalar product
!            with a potential change, meaning that its complex conjugate must be
!            used. This explains the following signs...
             do ipw1=ipw2,npwdiel
               susmat(1,ipw1,jsp,ipw2,jsp)=susmat(1,ipw1,jsp,ipw2,jsp)+&
&               weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
               susmat(2,ipw1,jsp,ipw2,jsp)=susmat(2,ipw1,jsp,ipw2,jsp)+&
&               weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
             end do
           end do
         end if

         if( occopt>=3 .and. testocc==1 .and. abs(wght1)>tol12) then
!          Accumulate product of band densities by their doccde, for the
!          computation of the effect of change of Fermi level.
           do ipw=1,npwdiel
             drhode(1,ipw,jsp)=drhode(1,ipw,jsp)+wfprod(1,ipw)*wght1
             drhode(2,ipw,jsp)=drhode(2,ipw,jsp)+wfprod(2,ipw)*wght1
           end do
!          Also accumulate weighted sum of doccde
           sumdocc=sumdocc+wght1
         end if

!        End condition of metallic occupancies or extrapolation
       end if

!      End loop on spinorial components
     end do
!    End loop on iband
   end do

   call timab(762,2,tsec)

!  -- Compute now off-band-diagonal terms ------------------------------------
!  -- Compute product of wavefunctions for different bands, inside the block -

   call timab(763,1,tsec)

!  if (occopt<3) then
   tolocc=1.0d-3
!  else
!  tolocc=1.0d-8
!  end if

   if(nbnd_current>1)then
     do ibd1=1,nbnd_current-1
       ibdshft1=ibd1+iband_shift
       do ibd2=ibd1+1,nbnd_current
         ibdshft2=ibd2+iband_shift

!        If the occupation numbers are sufficiently different, or
!        if extrapolation is used and the corresponding factor is not zero,
!        then there is a contribution
         occdiff=occ(ibdshft1+bdtot_index)-occ(ibdshft2+bdtot_index)
         if( abs(occdiff)>tolocc      .or. &
&         ( extrap==1 .and.            &
&         ( abs(occ_deavg(ibdshft1)) + abs(occ_deavg(ibdshft2)) ) >tolocc ) &
&         ) then

           eigdiff=eigen(ibdshft1+bdtot_index) - eigen(ibdshft2+bdtot_index)

!          Loop on spinorial components
           do ispinor=1,nspinor
             jsp=isp+ispinor-1;if (nspden_eff==1) jsp=isp

!            Store the contribution in wfraug
             do i3=1,ndiel3
               do i2=1,ndiel2
                 do i1=1,ndiel1
                   wfraug(1,i1,i2,i3)=wfrspa1(1,i1,i2,i3,ispinor,ibd1)*wfrspa1(1,i1,i2,i3,ispinor,ibd2)&
&                   +wfrspa1(2,i1,i2,i3,ispinor,ibd1)*wfrspa1(2,i1,i2,i3,ispinor,ibd2)
                   wfraug(2,i1,i2,i3)=wfrspa1(2,i1,i2,i3,ispinor,ibd1)*wfrspa1(1,i1,i2,i3,ispinor,ibd2)&
&                   -wfrspa1(1,i1,i2,i3,ispinor,ibd1)*wfrspa1(2,i1,i2,i3,ispinor,ibd2)
                 end do
               end do
             end do

!            Performs the Fourier Transform of the product, and store it in wfprod
             tim_fourwf=32
             call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&             1,kg_diel,kg_diel, mgfftdiel,mpi_enreg_diel,1,ngfftdiel,1,npwdiel,&
&             ndiel4,ndiel5,ndiel6,3,paral_kgb,tim_fourwf,weight,weight)

!            In case of PAW, add compensation charge contribution
             if (usepaw==1.and.neglect_pawhat==0) then
               call pawsushat(atindx,cprj_k,gbound_diel,gylmg_diel,ibd1,ibd2,ispinor,ispinor,1,kg_diel,&
&               lmax_diel,mgfftdiel,natom,nband_k,ndiel4,ndiel5,ndiel6,&
&               ngfftdiel,npwdiel,nspinor,ntypat,0,&
&               pawang,pawtab,ph3d_diel,typat,wfprod,dummy,&
&               mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
             end if

!            Perform now the summation
             wght1=zero ; wght2=zero
             if(abs(occdiff)>tolocc) wght1= occdiff/eigdiff * two*wtk(ikpt)/ucvol
             if(extrap==1) wght2=(occ_deavg(ibdshft1)+occ_deavg(ibdshft2)) * two*wtk(ikpt)/ucvol
             weight=wght1+wght2

             if (abs(weight)>tol12) then
               do ipw2=1,npwdiel
!                Only fills lower half of the matrix (here, the susceptibility matrix)
!                Note that wfprod of the first index must behave like a density,
!                so that it is used as generated by fourwf, while wfprod of the
!                second index will be implicitely used to make a scalar product
!                with a potential change, meaning that its complex conjugate must be
!                used. This explains the following signs...
                 do ipw1=ipw2,npwdiel
                   susmat(1,ipw1,jsp,ipw2,jsp)=susmat(1,ipw1,jsp,ipw2,jsp)+&
&                   weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
                   susmat(2,ipw1,jsp,ipw2,jsp)=susmat(2,ipw1,jsp,ipw2,jsp)+&
&                   weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
                 end do
               end do
             end if

!            End loop on spinorial components
           end do
!          End condition of different occupation numbers or extrapolation
         end if
!        End internal loop over bands
       end do
!      End external loop over bands
     end do
!    End condition of having more than one band
   end if

!  Loop on secondary block, with fast varying index, in decreasing order.
   if(iblk1/=nblk)then
     do iblk2=nblk,iblk1+1,-1
       iband_shift2=nbnd_in_blk1+(iblk2-2)*nbnd_in_blk

!      Loop over bands to fft and store Fourier transform of wavefunction
       iband_shift2=nbnd_in_blk1+(iblk2-2)*nbnd_in_blk
       do iband=1,nbnd_in_blk
!        Loop on spinorial components
         do ispinor=1,nspinor
           iwf=(ispinor-1)*npw_k+(iband+iband_shift2-1)*npw_k*nspinor+icg

!          Obtain Fourier transform in fft box
           tim_fourwf=22
           cwavef(:,1:npw_k)=cg(:,1+iwf:npw_k+iwf)
           call fourwf(1,rhoaug,cwavef,dummy,wfraug,gbound,gbound,&
&           istwf_k,kg_k,kg_k,mgfftdiel,mpi_enreg_diel,1,ngfftdiel,npw_k,1,&
&           ndiel4,ndiel5,ndiel6,0,paral_kgb,tim_fourwf,weight,weight)
           wfrspa2(:,:,:,:,ispinor,iband)=wfraug(:,:,:,:)
         end do
       end do

       do ibd1=1,nbnd_current
         ibdshft1=ibd1+iband_shift
         do ibd2=1,nbnd_in_blk
           ibdshft2=ibd2+iband_shift2

!          If the occupation numbers are sufficiently different, or
!          if extrapolation is used and the corresponding factor is not zero,
!          then there is a contribution
           occdiff=occ(ibdshft1+bdtot_index)-occ(ibdshft2+bdtot_index)
           if( abs(occdiff)>tolocc      .or. &
&           ( extrap==1 .and.            &
&           ( abs(occ_deavg(ibdshft1)) + abs(occ_deavg(ibdshft2)) ) >tolocc ) &
&           ) then

             eigdiff=eigen(ibdshft1+bdtot_index) - eigen(ibdshft2+bdtot_index)

!            Loop on spinorial components
             do ispinor=1,nspinor
               jsp=isp+ispinor-1;if (nspden_eff==1) jsp=isp

!              Store the contribution in wfraug
               do i3=1,ndiel3
                 do i2=1,ndiel2
                   do i1=1,ndiel1
                     wfraug(1,i1,i2,i3)=wfrspa1(1,i1,i2,i3,ispinor,ibd1)*wfrspa2(1,i1,i2,i3,ispinor,ibd2)&
&                     +wfrspa1(2,i1,i2,i3,ispinor,ibd1)*wfrspa2(2,i1,i2,i3,ispinor,ibd2)
                     wfraug(2,i1,i2,i3)=wfrspa1(2,i1,i2,i3,ispinor,ibd1)*wfrspa2(1,i1,i2,i3,ispinor,ibd2)&
&                     -wfrspa1(1,i1,i2,i3,ispinor,ibd1)*wfrspa2(2,i1,i2,i3,ispinor,ibd2)
                   end do
                 end do
               end do

!              Performs the Fourier Transform of the product, and store it in wfprod
               tim_fourwf=32
               call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&               1,kg_diel,kg_diel,mgfftdiel,mpi_enreg_diel,1,ngfftdiel,1,npwdiel,&
&               ndiel4,ndiel5,ndiel6,3,paral_kgb,tim_fourwf,weight,weight)

!              In case of PAW, add compensation charge contribution
               if (usepaw==1.and.neglect_pawhat==0) then
                 call pawsushat(atindx,cprj_k,gbound_diel,gylmg_diel,ibd1,ibdshft2,ispinor,ispinor,1,kg_diel,&
&                 lmax_diel,mgfftdiel,natom,nband_k,ndiel4,ndiel5,ndiel6,&
&                 ngfftdiel,npwdiel,nspinor,ntypat,0,&
&                 pawang,pawtab,ph3d_diel,typat,wfprod,dummy,&
&                 mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
               end if

!              Perform now the summation
               wght1=zero ; wght2=zero
               if(abs(occdiff)>tolocc) wght1= occdiff/eigdiff * two*wtk(ikpt)/ucvol
               if(extrap==1) wght2=(occ_deavg(ibdshft1)+occ_deavg(ibdshft2)) * two*wtk(ikpt)/ucvol
               weight=wght1+wght2

               if (abs(weight)>tol12) then
                 do ipw2=1,npwdiel
!                  Only fills lower half of the matrix (here, the susceptibility matrix)
!                  Note that wfprod of the first index must behave like a density,
!                  so that it is used as generated by fourwf, while wfprod of the
!                  second index will be implicitely used to make a scalar product
!                  with a potential change, meaning that its complex conjugate must be
!                  used. This explains the following signs...
                   do ipw1=ipw2,npwdiel
                     susmat(1,ipw1,jsp,ipw2,jsp)=susmat(1,ipw1,jsp,ipw2,jsp)+&
&                     weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
                     susmat(2,ipw1,jsp,ipw2,jsp)=susmat(2,ipw1,jsp,ipw2,jsp)+&
&                     weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
                   end do
                 end do
               end if

!              End loop on spinorial components
             end do
!            End condition of different occupation numbers or extrapolation
           end if
!          End internal loop over bands
         end do
!        End external loop over bands
       end do
!      End loop on bloks
     end do

!    Finish the loop on blok with iblk2=iblk1+1, so can use the
!    FFTd wavefunctions for the next iblk1.
     do iband=1,nbnd_in_blk
       wfrspa1(:,:,:,:,1:nspinor,iband)=wfrspa2(:,:,:,:,1:nspinor,iband)
     end do

!    End condition of iblk1/=nblk
   end if

   call timab(763,2,tsec)

!  End loop on iblk1
 end do

!DEBUG
!write(std_out,*)' suskmm : exit '
!do ipw1=1,npwdiel
!write(std_out,*)ipw1,susmat(1,ipw1,1,ipw1,1),susmat(2,ipw1,1,ipw1,1)
!end do
!write(std_out,*)' suskmm : end of susmat '
!stop
!ENDDEBUG

 call destroy_mpi_enreg(mpi_enreg_diel)
 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(dummy)
 ABI_DEALLOCATE(rhoaug)
 ABI_DEALLOCATE(wfprod)
 ABI_DEALLOCATE(wfraug)
 ABI_DEALLOCATE(wfrspa1)
 ABI_DEALLOCATE(wfrspa2)

 call timab(760,2,tsec)

end subroutine suskmm
!!***
