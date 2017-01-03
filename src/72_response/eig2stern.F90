!{\src2tex{textfont=tt}}
!!****f* ABINIT/eig2stern
!! NAME
!! eig2stern
!!
!! FUNCTION
!! This routine calculates the second-order eigenvalues.
!! The output eig2nkq is this quantity for the input k points.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (PB, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!!  bdeigrf = number of bands for which to calculate the second-order eigenvalues.
!!  clflg(3,mpert)= array on calculated perturbations for eig2rf.
!!  dim_eig2nkq = 1 if eig2nkq is to be computed.
!!  cg1_pert(2,mpw1*nspinor*mband*mk1mem*nsppol,3,mpert) = first-order wf in G 
!!            space for each perturbation. The wavefunction is orthogonal to the
!!            active space.
!!  gh0c1_pert(2,mpw1*nspinor*mband*mk1mem*nsppol,3,mpert) = matrix containing the
!!            vector:  <G|H(0)|psi(1)>, for each perturbation.
!!  gh1c_pert(2,mpw1*nspinor*mband*mk1mem*nsppol,3,mpert)) = matrix containing the
!!            vector:  <G|H(1)|n,k>, for each perturbation. The wavefunction is 
!!            orthogonal to the active space. 
!!  eigbrd(2,mband*nsppol,nkpt,3,natom,3,natom) = broadening factors for the 
!!            electronic eigenvalues (optional).
!!  eigen0(nkpt_rbz*mband*nsppol) = 0-order eigenvalues at all K-points: 
!!            <k,n'|H(0)|k,n'> (hartree).
!!  eigenq(nkpt_rbz*mband*nsppol) = 0-order eigenvalues at all shifted K-points:
!!            <k+Q,n'|H(0)|k+Q,n'> (hartree).
!!  eigen1(nkpt_rbz*2*nsppol*mband**2,3,mpert) = matrix of first-order: 
!!            <k+Q,n'|H(1)|k,n> (hartree) (calculated in dfpt_cgwf).
!!  eig2nkq(2,mband*nsppol,nkpt,3,natom,3,natom*dim_eig2nkq) = second derivatives of
!!            the electronic eigenvalues.
!!  elph2_imagden = imaginary part of the denominator of the sum-over-state expression
!!            for the electronic eigenenergy shift due to second-order electron-phonon
!!            interation.
!!  ieig2rf = integer for calculation type.
!!  indsym(4,nsym,natom) = indirect indexing array for atom labels
!!            (not used yet, but will be used with symmetries).
!!  istwfk_pert(nkpt_rbz,3,mpert) = integer for choice of storage of wavefunction at
!!            each k point for each perturbation.
!!  mband = maximum number of bands.
!!  mk1mem = maximum number of k points which can fit in memory (RF data);
!!            0 if use disk.
!!  mpert = maximum number of perturbations.
!!  natom = number of atoms in the unit cell.
!!  npert = number of phonon perturbations, without taking into account directions:
!!            natom. 
!!  nsym = number of symmetries (not used yet).
!!  mpi_enreg = informations about MPI parallelization.
!!  mpw1 = maximum number of planewaves used to represent first-order wavefunctions.
!!  nkpt_rbz = number of k-points for each perturbation.
!!  npwar1(nkpt_rbz,mpert) = number of planewaves at k-point for first-order.
!!  nspinor = number of spinorial components of the wavefunctions.
!!  nsppol = 1 for unpolarized, 2 for spin-polarized.
!!  occ(mband*nkpt*nsppol)=occup number for each band (often 2) at each k point
!!  smdelta = integer controling the calculation of electron lifetimes.
!!  symq(4,2,nsym) = 1 if symmetry preserves present qpoint. From littlegroup_q (not used yet).
!!  symrec(3,3,nsym) = 3x3 matrices of the group symmetries (reciprocal space)
!!            (not used yet).
!!  symrel(3,3,nsym) = array containing the symmetries in real space (not used yet).
!!  timrev = 1 if time-reversal preserves the q wavevector; 0 otherwise 
!!            (not in use yet).
!!  dtset = OPTIONAL, dataset structure containing the input variable of the
!!            calculation. This is required to use the k-interpolation routine.
!!  eigenq_fine(mband_fine,mkpt_fine,nsppol_fine) = OPTIONAL, 0-order eigenvalues
!!            at all shifted K-points: <k+Q,n'|H(0)|k+Q,n'> (hartree) of the
!!            fine grid. This information is read from the WF dense k-grid file.  
!!  hdr_fine = OPTIONAL, header of the WF file of the fine k-point grid. This
!!            variable is required for the k-interpolation routine.  
!!  hdr0     = OPTIONAL, header of the GS WF file of the corse k-point grid. This
!!            variable is required for the k-interpolation routine.  
!!
!! OUTPUT
!!  eig2nkq(2,mband*nsppol,nkpt_rbz,3,npert,3,npert)= diagonal part of the 
!!            second-order eigenvalues: E^{(2),diag}_{k,q,j}.
!!  eigbrd(2,mband*nsppol,nkpt_rbz,3,npert,3,npert)= OPTIONAL, array containing the
!!            electron lifetimes.
!!
!! PARENTS
!!      dfpt_looppert
!!
!! CHILDREN
!!      distrb2,dotprod_g,kptfine_av,smeared_delta,timab,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine eig2stern(occ,bdeigrf,clflg,cg1_pert,dim_eig2nkq,dim_eig2rf,eigen0,eigenq,&
&  eigen1,eig2nkq,elph2_imagden,esmear,gh0c1_pert,gh1c_pert,ieig2rf,istwfk_pert,&
&  mband,mk1mem,mpert,npert,mpi_enreg,mpw1,nkpt_rbz,npwar1,nspinor,nsppol,smdelta,&
&  dtset,eigbrd,eigenq_fine,hdr_fine,hdr0) 

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_cgtools
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eig2stern'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_51_manage_mpi
 use interfaces_72_response, except_this_one => eig2stern
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bdeigrf,dim_eig2nkq,dim_eig2rf,ieig2rf,mband,mk1mem,mpert,mpw1,nkpt_rbz
 integer,intent(in) :: npert,nspinor,nsppol,smdelta
 real(dp),intent(in) :: elph2_imagden,esmear
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: clflg(3,mpert)
 integer,intent(in) :: istwfk_pert(nkpt_rbz,3,mpert)
 integer,intent(in) :: npwar1(nkpt_rbz,mpert)
 real(dp),intent(in) :: cg1_pert(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf,3,mpert)
 real(dp),intent(in) :: gh0c1_pert(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf,3,mpert)
 real(dp),intent(in) :: gh1c_pert(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf,3,mpert)
 real(dp),intent(inout) :: eigen0(nkpt_rbz*mband*nsppol)
 real(dp),intent(in) :: eigen1(nkpt_rbz*2*nsppol*mband**2,3,mpert)
 real(dp),intent(inout) :: eigenq(nkpt_rbz*mband*nsppol)
 real(dp),intent(out) :: eig2nkq(2,mband*nsppol,nkpt_rbz,3,npert,3,npert*dim_eig2nkq)
 real(dp),intent(out),optional :: eigbrd(2,mband*nsppol,nkpt_rbz,3,npert,3,npert)
 real(dp),intent(in),pointer,optional :: eigenq_fine(:,:,:)
 real(dp), intent(in) :: occ(mband*nkpt_rbz*nsppol)
 type(dataset_type), intent(in) :: dtset
 type(hdr_type),intent(in),optional :: hdr_fine,hdr0

!Local variables-------------------------------
!tolerance for non degenerated levels
!scalars
 integer :: band2tot_index,band_index,bandtot_index,iband,icg2,idir1,idir2
 integer :: ikpt,ipert1,ipert2,isppol,istwf_k,jband,npw1_k,nkpt_sub,ikpt2
!integer :: ipw
 integer :: master,me,spaceworld,ierr
!real(dp),parameter :: etol=1.0d-3
 real(dp),parameter :: etol=1.0d-6
!real(dp),parameter :: etol=zero 
 real(dp) :: ar,ai,deltae,den,dot2i,dot2r,dot3i,dot3r,doti,dotr,eig1_i1,eig1_i2
 real(dp) :: eig1_r1,eig1_r2,eig2_diai,den_av
 real(dp) :: wgt_int
 real(dp) :: eig2_diar,eigbrd_i,eigbrd_r
 character(len=500) :: message
 character(len=500) :: msg
!DBSP
! character(len=300000) :: message2
!END
 logical :: test_do_band
!arrays
 integer, allocatable :: nband_rbz(:),icg2_rbz(:,:)
 integer,pointer      :: kpt_fine_sub(:)
 real(dp)             :: tsec(2)
 real(dp),allocatable :: cwavef(:,:),cwavef2(:,:),center(:),eigen0tmp(:),eigenqtmp(:)
 real(dp) :: eigen(mband*nsppol),eigen_prime(mband*nsppol)
 real(dp),allocatable :: gh(:,:),gh1(:,:),ghc(:,:)
 real(dp),allocatable :: smdfun(:,:)
 real(dp),pointer     :: wgt_sub(:)

! *********************************************************************

!Init parallelism
 master =0
 spaceworld=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt
!DEBUG
!write(std_out,*)' eig2stern : enter '
!write(std_out,*)' mpw1=',mpw1
!write(std_out,*)' mband=',mband
!write(std_out,*)' nsppol=',nsppol
!write(std_out,*)' nkpt_rbz=',nkpt_rbz
!write(std_out,*)' npert=',npert
!ENDDEBUG

!Init interpolation method
 if(present(eigenq_fine))then
   ABI_ALLOCATE(center,(3))
 end if

 call timab(148,1,tsec)

 if(nsppol==2)then
   message = 'nsppol=2 is still under development. Be careful when using it ...'
   MSG_COMMENT(message)
 end if        

 band2tot_index =0
 bandtot_index=0
 band_index=0

!Add scissor shift to eigenenergies
 if (dtset%dfpt_sciss > tol6 ) then
   write(msg,'(a,f7.3,2a)')&
&   ' A scissor operator of ',dtset%dfpt_sciss*Ha_eV,' [eV] has been applied to the eigenenergies',ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   ABI_ALLOCATE(eigen0tmp,(nkpt_rbz*mband*nsppol))
   ABI_ALLOCATE(eigenqtmp,(nkpt_rbz*mband*nsppol))
   eigen0tmp =   eigen0(:)
   eigenqtmp =   eigenq(:)
   eigen0 = zero
   eigenq = zero
 end if

 if(ieig2rf > 0) then
   eig2nkq(:,:,:,:,:,:,:) = zero
 end if
 if(present(eigbrd))then
   eigbrd(:,:,:,:,:,:,:) = zero
 end if

 if(xmpi_paral==1) then
   ABI_ALLOCATE(mpi_enreg%proc_distrb,(nkpt_rbz,mband,nsppol))
   ABI_ALLOCATE(nband_rbz,(nkpt_rbz*nsppol))
   if (allocated(mpi_enreg%my_kpttab)) then
     ABI_DEALLOCATE(mpi_enreg%my_kpttab)
   end if
   ABI_ALLOCATE(mpi_enreg%my_kpttab,(nkpt_rbz))
!  Assume the number of bands is the same for all k points.
   nband_rbz(:)=mband
   call distrb2(mband,nband_rbz,nkpt_rbz,mpi_enreg%nproc_cell,nsppol,mpi_enreg)
 end if

 icg2=0
 ipert1=1 ! Suppose that the situation is the same for all perturbations
 ABI_ALLOCATE(icg2_rbz,(nkpt_rbz,nsppol))
 do isppol=1,nsppol
   do ikpt=1,nkpt_rbz
     icg2_rbz(ikpt,isppol)=icg2
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,mband,isppol,me)) cycle
     icg2 = icg2 + npwar1(ikpt,ipert1)*nspinor*mband
   end do
 end do

 do isppol=1,nsppol
   do ikpt =1,nkpt_rbz

     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,mband,isppol,me)) then
       band2tot_index = band2tot_index + 2*mband**2
       bandtot_index = bandtot_index + mband
       cycle
     end if

     if(present(eigenq_fine))then
       write(std_out,*) 'Start of the energy denominator interpolation method.'
       nkpt_sub = 0
!      center is the k+q point around which we will average the kpt_fine
       center = hdr0%kptns(:,ikpt)+ dtset%qptn(:) 

       call kptfine_av(center,dtset%qptrlatt,hdr_fine%kptns,hdr_fine%nkpt,&
&       kpt_fine_sub,nkpt_sub,wgt_sub)
       write(std_out,'(a,3f8.4,a,i3)') 'Number of k-points of the fine grid &
&       around the k+Q point ',center,' is:',nkpt_sub
       write(std_out,'(a,f10.5)') 'The sum of the weights of the k-points is: ',SUM(wgt_sub)
     end if

!    Add scissor shift to eigenenergies
     if (dtset%dfpt_sciss > tol6 ) then
       do iband=1,mband
         if (occ(iband+bandtot_index) < tol6) then
           eigen0(iband+bandtot_index) = eigen0tmp(iband+bandtot_index) + dtset%dfpt_sciss
           eigenq(iband+bandtot_index) = eigenqtmp(iband+bandtot_index) + dtset%dfpt_sciss
         else
           eigen0(iband+bandtot_index) = eigen0tmp(iband+bandtot_index)
           eigenq(iband+bandtot_index) = eigenqtmp(iband+bandtot_index)
         end if
       end do
     end if


     if(smdelta >0) then   !broadening
       if(.not.allocated(smdfun))  then
         ABI_ALLOCATE(smdfun,(mband,mband))
       end if
       smdfun(:,:) = zero
       do iband=1,mband
         eigen(iband) = eigen0(iband+bandtot_index)
         eigen_prime(iband) =eigenq(iband+bandtot_index)
       end do
       if(esmear>tol6) then
         call smeared_delta(eigen,eigen_prime,esmear,mband,smdelta,smdfun)
       end if
     end if
     icg2=icg2_rbz(ikpt,isppol)

     ipert1=1 ! Suppose all perturbations lead to the same number of planewaves
     npw1_k = npwar1(ikpt,ipert1)
     ABI_ALLOCATE(cwavef,(2,npw1_k*nspinor))
     ABI_ALLOCATE(cwavef2,(2,npw1_k*nspinor))
     ABI_ALLOCATE(gh,(2,npw1_k*nspinor))
     ABI_ALLOCATE(gh1,(2,npw1_k*nspinor))
     ABI_ALLOCATE(ghc,(2,npw1_k*nspinor))

     do iband=1,bdeigrf

!      If the k point and band belong to me, compute the contribution
       test_do_band=.true.
       if(mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me)test_do_band=.false.
       
       if(test_do_band)then

         do ipert1=1,npert

           do idir1=1,3
             if(clflg(idir1,ipert1)==0)cycle
             istwf_k = istwfk_pert(ikpt,idir1,ipert1)

             do ipert2=1,npert
               do idir2=1,3
                 if(clflg(idir2,ipert2)==0)cycle

                 eig2_diar = zero ; eig2_diai = zero ; eigbrd_r = zero ; eigbrd_i = zero

                 do jband=1,mband
                   eig1_r1 = eigen1(2*jband-1+(iband-1)*2*mband+band2tot_index,idir1,ipert1)
                   eig1_r2 = eigen1(2*jband-1+(iband-1)*2*mband+band2tot_index,idir2,ipert2)
                   eig1_i1 = eigen1(2*jband+(iband-1)*2*mband+band2tot_index,idir1,ipert1)
                   eig1_i2 = - eigen1(2*jband+(iband-1)*2*mband+band2tot_index,idir2,ipert2) !the negative sign is from the CC
!                  If no interpolation, fallback on to the previous
!                  implementation
                   if(.not. present(eigenq_fine))then
                     deltae=eigenq(jband+bandtot_index)-eigen0(iband+bandtot_index)
                   end if
                   ar=eig1_r1*eig1_r2-eig1_i1*eig1_i2
                   ai=eig1_r1*eig1_i2+eig1_i1*eig1_r2

!                  Sum over all active space to retrieve the diagonal gauge
                   if(ieig2rf == 1 .or. ieig2rf ==2 ) then
!                    if(abs(deltae)>etol) then ! This is commented because
!                    there is no problem with divergencies with elph2_imag != 0
                     if( present(eigenq_fine))then
                       den_av = zero
                       wgt_int = zero
                       do ikpt2=1,nkpt_sub
                         deltae=eigenq_fine(jband,kpt_fine_sub(ikpt2),1)&
&                         -eigen0(iband+bandtot_index)
                         den_av = den_av-(wgt_sub(ikpt2)*deltae)/(deltae**2+elph2_imagden**2)
                         wgt_int = wgt_int+wgt_sub(ikpt2)
                       end do
                       den = den_av/wgt_int
                     else
                       if(abs(elph2_imagden) < etol) then  
                         if(abs(deltae)>etol) then
                           den=-one/(deltae**2+elph2_imagden**2)
                         else
                           den= zero
                         end if
                       else
                         den=-one/(deltae**2+elph2_imagden**2)
                       end if
                     end if
                     
!                    The following should be the most general implementation of the presence of elph2_imagden
!                    eig2_diar=eig2_diar+(ar*deltae+ai*elph2_imagden)*den
!                    eig2_diai=eig2_diai+(ai*deltae-ar*elph2_imagden)*den
!                    This gives back the implementation without elph2_imagden
!                    eig2_diar=eig2_diar+ar*deltae*den
!                    eig2_diai=eig2_diai+ai*deltae*den
!                    This is what Samuel had implemented
!                    eig2_diar=eig2_diar+ar*deltae*den
!                    eig2_diai=eig2_diai+ai*elph2_imagden*den
!                    Other possibility : throw away the broadening part, that is actually treated separately.
                     if( present(eigenq_fine))then
                       eig2_diar=eig2_diar+ar*den
                       eig2_diai=eig2_diai+ai*den
                     else
                       eig2_diar=eig2_diar+ar*deltae*den
                       eig2_diai=eig2_diai+ai*deltae*den
!DBSP              
!                       if (iband+band_index==2 .and. ikpt==1 .and. idir1==1 .and. ipert1==1 .and. idir2==1 .and. ipert2==1) then
!                         write(message2,*) 'eig2_diar1=',eig2_diar,' ar=',ar,' deltae=',deltae,' den=',den
!                         call wrtout(std_out,message2,'PERS')
!                       endif
!END

                     end if
                   end if ! ieig2rf==1 or 2

                   if(present(eigbrd))then
                     if(smdelta >0) then   !broadening
                       eigbrd_r = eigbrd_r + ar*smdfun(iband,jband)
                       eigbrd_i = eigbrd_i + ai*smdfun(iband,jband)
                     end if
                   end if

                 end do !jband

!                Add the contribution of non-active bands, if DFPT calculation (= Sternheimer)
                 if(ieig2rf == 1 .or. ieig2rf ==3 .or. ieig2rf ==4 .or. ieig2rf==5 ) then
!                  if(ieig2rf == 1   ) then

                   dotr=zero ; doti=zero
                   dot2r=zero ; dot2i=zero
                   dot3r=zero ; dot3i=zero


                   cwavef(:,:) = cg1_pert(:,1+(iband-1)*npw1_k*nspinor+icg2:iband*npw1_k*nspinor+icg2,idir2,ipert2)
                   cwavef2(:,:)= cg1_pert(:,1+(iband-1)*npw1_k*nspinor+icg2:iband*npw1_k*nspinor+icg2,idir1,ipert1)
                   gh1(:,:)    = gh1c_pert(:,1+(iband-1)*npw1_k*nspinor+icg2:iband*npw1_k*nspinor+icg2,idir1,ipert1)
                   gh(:,:)     = gh1c_pert(:,1+(iband-1)*npw1_k*nspinor+icg2:iband*npw1_k*nspinor+icg2,idir2,ipert2)
                   ghc(:,:)    = gh0c1_pert(:,1+(iband-1)*npw1_k*nspinor+icg2:iband*npw1_k*nspinor+icg2,idir1,ipert1)

!                  The first two dotprod corresponds to:  <Psi(1)nkq|H(1)k+q,k|Psi(0)nk> and <Psi(0)nk|H(1)k,k+q|Psi(1)nkq>
!                  They are calculated using wavefunctions <Psi(1)| that are orthogonal to the active space.
                   call dotprod_g(dotr,doti,istwf_k,npw1_k*nspinor,2,cwavef,gh1,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
                   call dotprod_g(dot2r,dot2i,istwf_k,npw1_k*nspinor,2,gh,cwavef2,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

!                  This dotprod corresponds to : <Psi(1)nkq|H(0)k+q- E(0)nk|Psi(1)nkq>
!                  It is calculated using wavefunctions that are orthogonal to the active space.
!                  Should work for metals. (But adiabatic approximation is bad in this case...)
                   call dotprod_g(dot3r,dot3i,istwf_k,npw1_k*nspinor,2,cwavef,ghc,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

                   eig2_diar= eig2_diar + dotr + dot2r + dot3r
                   eig2_diai= eig2_diai + doti + dot2i + dot3i

                 end if

!                Store the contribution
                 if(ieig2rf > 0) then
                   eig2nkq(1,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = eig2_diar
                   eig2nkq(2,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = eig2_diai 
                 end if

                 if(present(eigbrd))then
                   if(smdelta >0) then   !broadening
                     eigbrd(1,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = eigbrd_r
                     eigbrd(2,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = eigbrd_i
                   end if
                 end if

               end do !idir2
             end do !ipert2
           end do  !idir1
         end do   !ipert1

       end if ! Selection of processor
       
     end do !iband

     ABI_DEALLOCATE(cwavef)
     ABI_DEALLOCATE(cwavef2)
     ABI_DEALLOCATE(gh)
     ABI_DEALLOCATE(gh1)
     ABI_DEALLOCATE(ghc)
     band2tot_index = band2tot_index + 2*mband**2
     bandtot_index = bandtot_index + mband

     if(present(eigenq_fine))then
       ABI_DEALLOCATE(kpt_fine_sub) ! Deallocate the variable
       ABI_DEALLOCATE(wgt_sub)
     end if

   end do    !ikpt
   band_index = band_index + mband
 end do !isppol

!Accumulate eig2nkq and/or eigbrd
 if(xmpi_paral==1) then
   if(ieig2rf == 1 .or. ieig2rf == 2) then
     call xmpi_sum(eig2nkq,spaceworld,ierr)
     if (dtset%dfpt_sciss > tol6 ) then
       call xmpi_sum(eigen0,spaceworld,ierr)
       call xmpi_sum(eigenq,spaceworld,ierr)
     end if
   end if
   if(present(eigbrd) .and. (ieig2rf == 1 .or. ieig2rf == 2))then
     if(smdelta >0) then
       call xmpi_sum(eigbrd,spaceworld,ierr)
     end if
   end if
   ABI_DEALLOCATE(nband_rbz)
   ABI_DEALLOCATE(mpi_enreg%proc_distrb)
   ABI_DEALLOCATE(mpi_enreg%my_kpttab)
 end if 

 if(ieig2rf==1 .or. ieig2rf==2 ) then
   write(ab_out,'(a)')' Components of second-order derivatives of the electronic energy, EIGR2D.'
   write(ab_out,'(a)')' For automatic tests, printing the matrix for the first k-point, first band, first atom.'
   do idir1=1,3
     do idir2=1,3
       ar=eig2nkq(1,1,1,idir1,1,idir2,1) ; if(abs(ar)<tol10)ar=zero
       ai=eig2nkq(2,1,1,idir1,1,idir2,1) ; if(abs(ai)<tol10)ai=zero
       write (ab_out,'(4i4,2es20.10)') idir1,1,idir2,1,ar,ai 
     end do ! idir2
   end do ! idir1
 end if 

 if(present(eigbrd))then
   if(smdelta >0) then   !broadening
     write(ab_out,'(a)')' '
     write(ab_out,'(a)')' Components of second-order derivatives of the electronic energy, EIGI2D.'
     write(ab_out,'(a)')' For automatic tests, printing the matrix for the first k-point, first band, first atom.'
     do idir1=1,3
       do idir2=1,3
         ar=eigbrd(1,1,1,idir1,1,idir2,1) ; if(abs(ar)<tol10)ar=zero
         ai=eigbrd(2,1,1,idir1,1,idir2,1) ; if(abs(ai)<tol10)ai=zero
         write (ab_out,'(4i4,2es20.10)') idir1,1,idir2,1,ar,ai
       end do
     end do !nband
   end if
 end if

 if(allocated(smdfun))  then
   ABI_DEALLOCATE(smdfun)
 end if
 ABI_DEALLOCATE(icg2_rbz)
 if(present(eigenq_fine))then
   ABI_DEALLOCATE(center)
 end if
 if (dtset%dfpt_sciss > tol6 ) then
   ABI_DEALLOCATE(eigen0tmp)
   ABI_DEALLOCATE(eigenqtmp)
 end if

 call timab(148,2,tsec)

!DEBUG
!write(std_out,*)' eig2stern: exit'
!ENDDEBUG

end subroutine eig2stern
!!***

