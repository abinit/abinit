!{\src2tex{textfont=tt}}
!!****f* ABINIT/eig2tot
!! NAME
!! eig2tot
!!
!! FUNCTION
!! This routine calculates the second-order eigenvalues.
!! The output eig2nkq is this quantity for the input k points.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (SP,PB,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!!  bdeigrf = number of bands for which to calculate the second-order eigenvalues.
!!  clflg(3,mpert)= array on calculated perturbations for eig2rf.
!!  dim_eig2nkq = 1 if eig2nkq is to be computed.
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
!!  mband = maximum number of bands.
!!  mpert = maximum number of perturbations.
!!  natom = number of atoms in the unit cell.
!!  npert = number of phonon perturbations, without taking into account directions:
!!            natom. 
!!  nsym = number of symmetries (not used yet).
!!  mpi_enreg = informations about MPI parallelization.
!!  nkpt_rbz = number of k-points for each perturbation.
!!  nsppol = 1 for unpolarized, 2 for spin-polarized.
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
!!  hdr0     = header of the GS WF file of the corse k-point grid. 
!!            
!!
!! OUTPUT
!!  eig2nkq(2,mband*nsppol,nkpt_rbz,3,npert,3,npert)= diagonal part of the 
!!            second-order eigenvalues: E^{(2),diag}_{k,q,j}.
!!  eigbrd(2,mband*nsppol,nkpt_rbz,3,npert,3,npert)= OPTIONAL, array containing the
!!            electron lifetimes.
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      crystal_free,crystal_init,ddb_io_out,distrb2,ebands_free,ebands_init
!!      eigr2d_free,eigr2d_init,eigr2d_ncwrite,fan_free,fan_init,fan_ncwrite
!!      gkk_free,gkk_init,gkk_ncwrite,kptfine_av,outbsd,psddb8,smeared_delta
!!      timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine eig2tot(dtfil,xred,psps,pawtab,natom,bdeigrf,clflg,dim_eig2nkq,eigen0,eigenq,eigen1,eig2nkq,&
&  elph2_imagden,esmear,ieig2rf,mband,mpert,npert,mpi_enreg,doccde,&
&  nkpt_rbz,nsppol,smdelta,rprimd,dtset,occ_rbz,hdr0,eigbrd,eigenq_fine,hdr_fine) 


 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_nctk
 use m_errors
#ifdef HAVE_TRIO_NETCDF
 use netcdf
#endif
 use m_ebands

 use m_fstrings,   only : strcat
 use m_eig2d,      only : eigr2d_init,eigr2d_t,eigr2d_ncwrite,eigr2d_free,fan_t,&
                          & fan_init,fan_ncwrite,fan_free, gkk_t, gkk_init, &
                          & gkk_ncwrite,gkk_free
 use m_crystal,    only : crystal_init, crystal_free, crystal_t
 use m_crystal_io, only : crystal_ncwrite
 use m_pawtab,   only : pawtab_type
 use m_ddb,      only : psddb8

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eig2tot'
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_51_manage_mpi
 use interfaces_72_response, except_this_one => eig2tot
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bdeigrf,dim_eig2nkq,ieig2rf,mband,mpert,natom,nkpt_rbz
 integer,intent(in) :: npert,nsppol,smdelta
 real(dp),intent(in) :: elph2_imagden,esmear
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type), intent(in) :: dtfil
 type(pseudopotential_type), intent(inout) :: psps
!arrays
 type(dataset_type), intent(in) :: dtset
 integer,intent(in) :: clflg(3,mpert)
 real(dp),intent(in) :: doccde(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: eigen0(nkpt_rbz*mband*nsppol)
 real(dp),intent(in) :: eigen1(nkpt_rbz*2*nsppol*mband**2,3,mpert)
 real(dp),intent(in) :: eigenq(nkpt_rbz*mband*nsppol)
 real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol)
 real(dp),intent(inout) :: eig2nkq(2,mband*nsppol,nkpt_rbz,3,npert,3,npert*dim_eig2nkq)
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom)
 real(dp),intent(inout),optional :: eigbrd(2,mband*nsppol,nkpt_rbz,3,npert,3,npert)
 real(dp),intent(in),pointer,optional :: eigenq_fine(:,:,:)
 type(pawtab_type), intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 type(hdr_type),intent(in) :: hdr0
 type(hdr_type),intent(in),optional :: hdr_fine

!Local variables-------------------------------
!tolerance for non degenerated levels
!scalars
 integer :: band2tot_index,band_index,bantot,bandtot_index,choice,fullinit,iband,idir1,idir2
 integer :: ikpt,ipert1,ipert2,isppol,jband,nkpt_sub,ikpt2,nblok,unitout,vrsddb,ncid
!integer :: ipw
 character(len=fnlen) :: dscrpt,fname
 integer :: master,me,spaceworld,ierr
! real(dp),parameter :: etol=1.0d-6
 real(dp),parameter :: etol=1.0d-7
!real(dp),parameter :: etol=zero 
 real(dp) :: ar,ai,deltae,den,eig1_i1,eig1_i2,eigen_corr
 real(dp) :: eig1_r1,eig1_r2,eig2_diai,den_av
 real(dp) :: eig2_diar,eigbrd_i,eigbrd_r,tolwfr,wgt_int
 character(len=500) :: message
 logical :: remove_inv,test_do_band
 type(crystal_t) :: Crystal
 type(ebands_t)  :: Bands
 type(eigr2d_t)  :: eigr2d,eigi2d
 type(fan_t)     :: fan2d
 type(gkk_t)     :: gkk2d
!arrays
 integer, allocatable :: nband_rbz(:)
 integer,pointer      :: kpt_fine_sub(:)
 real(dp)             :: tsec(2)
 real(dp),allocatable :: center(:)
 real(dp) :: eigen(mband*nsppol),eigen_prime(mband*nsppol)
 real(dp),allocatable :: fan(:,:,:,:,:,:,:)
 real(dp),allocatable :: gkk(:,:,:,:,:)
 real(dp),allocatable :: eig2nkq_tmp(:,:,:,:,:,:,:)
 real(dp),allocatable :: smdfun(:,:)
 real(dp),pointer     :: wgt_sub(:)

! *********************************************************************

!Init parallelism
 master =0
 spaceworld=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt
!DEBUG
!write(std_out,*)' eig2tot : enter '
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

 if(ieig2rf == 4 ) then
   ABI_STAT_ALLOCATE(fan,(2*mband*nsppol,dtset%nkpt,3,natom,3,natom*dim_eig2nkq,mband), ierr)
   ABI_CHECK(ierr==0, "out-of-memory in fan")
   fan(:,:,:,:,:,:,:) = zero
   ABI_STAT_ALLOCATE(eig2nkq_tmp,(2,mband*nsppol,dtset%nkpt,3,natom,3,natom*dim_eig2nkq), ierr)
   ABI_CHECK(ierr==0, "out-of-memory in eig2nkq_tmp")
   eig2nkq_tmp(:,:,:,:,:,:,:) = zero
!  This is not efficient because double the memory. Alternative: use buffer and
!  print part by part.
   eig2nkq_tmp = eig2nkq
   if(present(eigbrd))then
     eigbrd(:,:,:,:,:,:,:)=zero
   end if
   eigen_corr = 0
 end if

 if(ieig2rf == 5 ) then
   ABI_STAT_ALLOCATE(gkk,(2*mband*nsppol,dtset%nkpt,3,natom,mband), ierr)
   ABI_CHECK(ierr==0, "out-of-memory in gkk")
   gkk(:,:,:,:,:) = zero
   ABI_STAT_ALLOCATE(eig2nkq_tmp,(2,mband*nsppol,dtset%nkpt,3,natom,3,natom*dim_eig2nkq), ierr)
   ABI_CHECK(ierr==0, "out-of-memory in eig2nkq_tmp")
   eig2nkq_tmp(:,:,:,:,:,:,:) = zero
!  This is not efficient because double the memory. Alternative: use buffer and
!  print part by part.
   eig2nkq_tmp = eig2nkq
   if(present(eigbrd))then
     eigbrd(:,:,:,:,:,:,:)=zero
   end if
   eigen_corr = 0
 end if
 
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

     ipert1=1 ! Suppose all perturbations lead to the same number of planewaves

     do iband=1,bdeigrf

!      If the k point and band belong to me, compute the contribution
       test_do_band=.true.
       if(mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me)test_do_band=.false.
       
       if(test_do_band)then
!        ------------------------------------------------------------------------------------------------------!
!        ------- ieig2rf ==3 : Non dynamic traditional AHC theory with Sternheimer (computed in eig2stern.F90)-!
!        ------------------------------------------------------------------------------------------------------!
!        Note that ieig2rf==4 and ieig2rf==5 also goes into that part only for later printing of the ZPR in the ouput of abinit
!        later in the code
         if(ieig2rf==3 .or. ieig2rf==4 .or. ieig2rf==5) then
           do ipert1=1,npert
             do idir1=1,3
               if(clflg(idir1,ipert1)==0) cycle
               do ipert2=1,npert
                 do idir2=1,3
                   if(clflg(idir2,ipert2)==0)cycle
                   eig2_diar = zero ; eig2_diai = zero ; eigbrd_r = zero ; eigbrd_i = zero
                   do jband=1,mband
                     eig1_r1 = eigen1(2*jband-1+(iband-1)*2*mband+band2tot_index,idir1,ipert1)
                     eig1_r2 = eigen1(2*jband-1+(iband-1)*2*mband+band2tot_index,idir2,ipert2)
                     eig1_i1 = eigen1(2*jband+(iband-1)*2*mband+band2tot_index,idir1,ipert1)
                     eig1_i2 = - eigen1(2*jband+(iband-1)*2*mband+band2tot_index,idir2,ipert2) !the negative sign is from the CC
!                    If no interpolation, fallback on to the previous
!                    implementation
                     if(.not. present(eigenq_fine))then
                       deltae=eigenq(jband+bandtot_index)-eigen0(iband+bandtot_index)
                     end if
                     ar=eig1_r1*eig1_r2-eig1_i1*eig1_i2
                     ai=eig1_r1*eig1_i2+eig1_i1*eig1_r2

!                    Sum over all active space to retrieve the diagonal gauge
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
                     
                     if( present(eigenq_fine))then
                       eig2_diar=eig2_diar+ar*den
                       eig2_diai=eig2_diai+ai*den
                     else
                       eig2_diar=eig2_diar+ar*deltae*den
                       eig2_diai=eig2_diai+ai*deltae*den
                     end if

                     if(present(eigbrd))then
                       if(smdelta >0) then   !broadening
                         eigbrd_r = eigbrd_r + ar*smdfun(iband,jband)
                         eigbrd_i = eigbrd_i + ai*smdfun(iband,jband)
                       end if
                     end if
                   end do !jband

!                  Store the contribution
                   eig2nkq(1,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = &
&                   eig2nkq(1,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) + eig2_diar
                   eig2nkq(2,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = &
&                   eig2nkq(2,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) + eig2_diai 

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
         end if !ieig2rf 3

!        -------------------------------------------------------------------------------------------!
!        ------- ieig2rf ==4  Dynamic AHC using second quantization and Sternheimer from eig2stern -!
!        -------------------------------------------------------------------------------------------!   
         if(ieig2rf ==4 ) then
           do ipert1=1,npert
             do idir1=1,3
               if(clflg(idir1,ipert1)==0) cycle
               do ipert2=1,npert
                 do idir2=1,3
                   if(clflg(idir2,ipert2)==0)cycle
                   do jband=1,mband
                     eig1_r1 = eigen1(2*jband-1+(iband-1)*2*mband+band2tot_index,idir1,ipert1)
                     eig1_r2 = eigen1(2*jband-1+(iband-1)*2*mband+band2tot_index,idir2,ipert2)
                     eig1_i1 = eigen1(2*jband+(iband-1)*2*mband+band2tot_index,idir1,ipert1)
                     eig1_i2 = - eigen1(2*jband+(iband-1)*2*mband+band2tot_index,idir2,ipert2) !the negative sign is from the CC
                     ar=eig1_r1*eig1_r2-eig1_i1*eig1_i2
                     ai=eig1_r1*eig1_i2+eig1_i1*eig1_r2
!                  Store the contribution
                     fan(2*iband-1+2*band_index,ikpt,idir1,ipert1,idir2,ipert2,jband) = &
&                     fan(2*iband-1+2*band_index,ikpt,idir1,ipert1,idir2,ipert2,jband) + ar
                     fan(2*iband+2*band_index,ikpt,idir1,ipert1,idir2,ipert2,jband) = &
&                     fan(2*iband+2*band_index,ikpt,idir1,ipert1,idir2,ipert2,jband) + ai
                   end do !jband
                 end do !idir2
               end do !ipert2
             end do  !idir1
           end do   !ipert1
         end if !ieig2rf 4
!        --------------------------------------------------------------------------------!
!        ------- ieig2rf ==5  Dynamic AHC with Sternheimer from eig2stern but print GKK -!
!        --------------------------------------------------------------------------------!   
         if(ieig2rf ==5 ) then
           do ipert1=1,npert
             do idir1=1,3
               if(clflg(idir1,ipert1)==0) cycle
               do jband=1,mband
                 eig1_r1 = eigen1(2*jband-1+(iband-1)*2*mband+band2tot_index,idir1,ipert1)
                 eig1_i1 = eigen1(2*jband+(iband-1)*2*mband+band2tot_index,idir1,ipert1)
!              Store the contribution
                 gkk(2*iband-1+2*band_index,ikpt,idir1,ipert1,jband) = &
&                 gkk(2*iband-1+2*band_index,ikpt,idir1,ipert1,jband) + eig1_r1
                 gkk(2*iband+2*band_index,ikpt,idir1,ipert1,jband) = &
&                 gkk(2*iband+2*band_index,ikpt,idir1,ipert1,jband) + eig1_i1
               end do !jband
             end do  !idir1
           end do   !ipert1
         end if !ieig2rf 5
       end if ! Selection of processor
     end do !iband

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
   if(ieig2rf == 3) then
     call xmpi_sum(eig2nkq,spaceworld,ierr)
   end if
   if(ieig2rf == 4) then
     call xmpi_sum(eig2nkq,spaceworld,ierr)
     call xmpi_sum(eig2nkq_tmp,spaceworld,ierr)
     call xmpi_sum(fan,spaceworld,ierr)
   end if
   if(ieig2rf == 5) then
     call xmpi_sum(eig2nkq,spaceworld,ierr)
     call xmpi_sum(eig2nkq_tmp,spaceworld,ierr)
     call xmpi_sum(gkk,spaceworld,ierr)
   end if
   if(present(eigbrd) .and. (ieig2rf == 3 .or. ieig2rf == 4 .or. ieig2rf == 5))then
     if(smdelta >0) then
       call xmpi_sum(eigbrd,spaceworld,ierr)
     end if
   end if
   ABI_DEALLOCATE(nband_rbz)
   ABI_DEALLOCATE(mpi_enreg%proc_distrb)
   ABI_DEALLOCATE(mpi_enreg%my_kpttab)
 end if

 if(ieig2rf > 2) then
   write(ab_out,'(a)')' Components of second-order derivatives of the electronic energy, EIGR2D.'
   write(ab_out,'(a)')' For automatic tests, printing the matrix for the first k-point, first band, first atom.'
   band_index = 0
   do isppol=1,dtset%nsppol
     do idir1=1,3
       do idir2=1,3
         ar=eig2nkq(1,1+band_index,1,idir1,1,idir2,1) ; if(abs(ar)<tol10)ar=zero
         ai=eig2nkq(2,1+band_index,1,idir1,1,idir2,1) ; if(abs(ai)<tol10)ai=zero
         write (ab_out,'(4i4,2es20.10)') idir1,1,idir2,1,ar,ai 
       end do ! idir2
     end do ! idir1
     band_index = band_index + mband
     write(ab_out,'(a)')' '
   end do
 end if 

 if(present(eigbrd))then
   if(smdelta >0) then   !broadening
     write(ab_out,'(a)')' Components of second-order derivatives of the electronic energy, EIGI2D.'
     write(ab_out,'(a)')' For automatic tests, printing the matrix for the first k-point, first band, first atom.'
     band_index = 0
     do isppol=1,dtset%nsppol
       do idir1=1,3
         do idir2=1,3
           ar=eigbrd(1,1+band_index,1,idir1,1,idir2,1) ; if(abs(ar)<tol10)ar=zero
           ai=eigbrd(2,1+band_index,1,idir1,1,idir2,1) ; if(abs(ai)<tol10)ai=zero
           write (ab_out,'(4i4,2es20.10)') idir1,1,idir2,1,ar,ai
         end do
       end do
       band_index = band_index + mband
       write(ab_out,'(a)')' '
     end do
   end if
 end if

 if(allocated(smdfun))  then
   ABI_DEALLOCATE(smdfun)
 end if
 if(present(eigenq_fine))then
   ABI_DEALLOCATE(center)
 end if

 master=0
 if (me==master) then
!  print _EIGR2D file for this perturbation in the case of ieig2rf 3 or 4 or 5
   if (ieig2rf == 3 .or. ieig2rf == 4 .or. ieig2rf == 5) then 
     unitout = dtfil%unddb
     vrsddb=100401
     dscrpt=' Note : temporary (transfer) database '
!    tolwfr must be initialized here, but it is a dummy value
     tolwfr=1.0_dp
     call ddb_io_out (dscrpt,dtfil%fnameabo_eigr2d,dtset%natom,dtset%mband,&
&     dtset%nkpt,dtset%nsym,dtset%ntypat,dtfil%unddb,vrsddb,&
&     dtset%acell_orig(1:3,1),dtset%amu_orig(:,1),dtset%dilatmx,dtset%ecut,dtset%ecutsm,&
&     dtset%intxc,dtset%iscf,dtset%ixc,dtset%kpt,dtset%kptnrm,&
&     dtset%natom,dtset%nband,dtset%ngfft,dtset%nkpt,dtset%nspden,dtset%nspinor,&
&     dtset%nsppol,dtset%nsym,dtset%ntypat,occ_rbz,dtset%occopt,dtset%pawecutdg,&
&     dtset%rprim_orig(1:3,1:3,1),dtset%dfpt_sciss,dtset%spinat,dtset%symafm,dtset%symrel,&
&     dtset%tnons,tolwfr,dtset%tphysel,dtset%tsmear,&
&     dtset%typat,dtset%usepaw,dtset%wtk,xred,psps%ziontypat,dtset%znucl)
     nblok=1 ; fullinit=1 ; choice=2
     call psddb8 (choice,psps%dimekb,psps%ekb,fullinit,psps%indlmn,&
&     psps%lmnmax,nblok,dtset%ntypat,dtfil%unddb,pawtab,&
&     psps%pspso,psps%usepaw,psps%useylm,vrsddb)
   end if
   if(ieig2rf == 3 ) then
     call outbsd(bdeigrf,dtset,eig2nkq,dtset%natom,nkpt_rbz,unitout)
   end if
   if(ieig2rf == 4 .or. ieig2rf == 5 ) then
     call outbsd(bdeigrf,dtset,eig2nkq_tmp,dtset%natom,nkpt_rbz,unitout)
   end if
!  Output of the EIGR2D.nc file.
   fname = strcat(dtfil%filnam_ds(4),"_EIGR2D.nc")
!  Crystalline structure.
   remove_inv=.false.
   if(dtset%nspden==4 .and. dtset%usedmft==1) remove_inv=.true.
   call crystal_init(Crystal,dtset%spgroup,dtset%natom,dtset%npsp,psps%ntypat, &
&   dtset%nsym,rprimd,dtset%typat,xred,dtset%ziontypat,dtset%znucl,1,&
&   dtset%nspden==2.and.dtset%nsppol==1,remove_inv,hdr0%title,&
&   dtset%symrel,dtset%tnons,dtset%symafm)
!  Electronic band energies.
   bantot= dtset%mband*dtset%nkpt*dtset%nsppol
   call ebands_init(bantot,Bands,dtset%nelect,doccde,eigen0,hdr0%istwfk,hdr0%kptns,&
&   hdr0%nband, hdr0%nkpt,hdr0%npwarr,hdr0%nsppol,hdr0%nspinor,&
&   hdr0%tphysel,hdr0%tsmear,hdr0%occopt,hdr0%occ,hdr0%wtk,&
&   hdr0%charge, hdr0%kptopt, hdr0%kptrlatt_orig, hdr0%nshiftk_orig, hdr0%shiftk_orig, &
&   hdr0%kptrlatt, hdr0%nshiftk, hdr0%shiftk)

!  Second order derivative EIGR2D (real and Im)
   if(ieig2rf == 3 ) then
     call eigr2d_init(eig2nkq,eigr2d,dtset%mband,hdr0%nsppol,nkpt_rbz,dtset%natom)
   end if
   if(ieig2rf == 4 .or. ieig2rf == 5 ) then
     call eigr2d_init(eig2nkq_tmp,eigr2d,dtset%mband,hdr0%nsppol,nkpt_rbz,dtset%natom)
   end if
#ifdef HAVE_TRIO_NETCDF
   NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating EIGR2D file")
   NCF_CHECK(crystal_ncwrite(Crystal,ncid))
   NCF_CHECK(ebands_ncwrite(Bands, ncid))
   call eigr2d_ncwrite(eigr2d,dtset%qptn(:),dtset%wtq,ncid)
   NCF_CHECK(nf90_close(ncid))
#else
   ABI_UNUSED(ncid)
#endif

!  print _FAN file for this perturbation. Note that the Fan file will only be produced if
!  abinit is compiled with netcdf.
   if(ieig2rf == 4 ) then
!    Output of the Fan.nc file.
#ifdef HAVE_TRIO_NETCDF
     fname = strcat(dtfil%filnam_ds(4),"_FAN.nc")
     call fan_init(fan,fan2d,dtset%mband,hdr0%nsppol,nkpt_rbz,dtset%natom)
     NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating FAN file")
     NCF_CHECK(crystal_ncwrite(Crystal, ncid))
     NCF_CHECK(ebands_ncwrite(Bands, ncid))
     call fan_ncwrite(fan2d,dtset%qptn(:),dtset%wtq, ncid)
     NCF_CHECK(nf90_close(ncid))
#else
     MSG_ERROR("Dynamical calculation with ieig2rf 4 only work with NETCDF support.")
     ABI_UNUSED(ncid)
#endif
     ABI_DEALLOCATE(fan)
     ABI_DEALLOCATE(eig2nkq_tmp)
   end if
!  print _GKK.nc file for this perturbation. Note that the GKK file will only be produced if
!  abinit is compiled with netcdf.
   if(ieig2rf == 5 ) then
!    Output of the GKK.nc file.
#ifdef HAVE_TRIO_NETCDF
     fname = strcat(dtfil%filnam_ds(4),"_GKK.nc")
     call gkk_init(gkk,gkk2d,dtset%mband,hdr0%nsppol,nkpt_rbz,dtset%natom,3)
     NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating GKK file")
     NCF_CHECK(crystal_ncwrite(Crystal, ncid))
     NCF_CHECK(ebands_ncwrite(Bands, ncid))
     call gkk_ncwrite(gkk2d,dtset%qptn(:),dtset%wtq, ncid)
     NCF_CHECK(nf90_close(ncid))
#else
     MSG_ERROR("Dynamical calculation with ieig2rf 5 only work with NETCDF support.")
     ABI_UNUSED(ncid)
#endif
     ABI_DEALLOCATE(gkk)
     ABI_DEALLOCATE(eig2nkq_tmp)
   end if
!  print _EIGI2D file for this perturbation
   if (ieig2rf /= 5 ) then
     if(smdelta>0) then
       unitout = dtfil%unddb
       vrsddb=100401
       dscrpt=' Note : temporary (transfer) database '
!      tolwfr must be initialized here, but it is a dummy value
       tolwfr=1.0_dp
       call ddb_io_out (dscrpt,dtfil%fnameabo_eigi2d,dtset%natom,dtset%mband,&
&       dtset%nkpt,dtset%nsym,dtset%ntypat,dtfil%unddb,vrsddb,&
&       dtset%acell_orig(1:3,1),dtset%amu_orig(:,1),dtset%dilatmx,dtset%ecut,dtset%ecutsm,&
&       dtset%intxc,dtset%iscf,dtset%ixc,dtset%kpt,dtset%kptnrm,&
&       dtset%natom,dtset%nband,dtset%ngfft,dtset%nkpt,dtset%nspden,dtset%nspinor,&
&       dtset%nsppol,dtset%nsym,dtset%ntypat,occ_rbz,dtset%occopt,dtset%pawecutdg,&
&       dtset%rprim_orig(1:3,1:3,1),dtset%dfpt_sciss,dtset%spinat,dtset%symafm,dtset%symrel,&
&       dtset%tnons,tolwfr,dtset%tphysel,dtset%tsmear,&
&       dtset%typat,dtset%usepaw,dtset%wtk,xred,psps%ziontypat,dtset%znucl)

       nblok=1 ; fullinit=1 ; choice=2
       call psddb8 (choice,psps%dimekb,psps%ekb,fullinit,psps%indlmn,&
&       psps%lmnmax,nblok,dtset%ntypat,dtfil%unddb,pawtab,&
&       psps%pspso,psps%usepaw,psps%useylm,vrsddb)

       call outbsd(bdeigrf,dtset,eigbrd,dtset%natom,nkpt_rbz,unitout)

!      Output of the EIGI2D.nc file.
       fname = strcat(dtfil%filnam_ds(4),"_EIGI2D.nc")
!      Broadening EIGI2D (real and Im)
       call eigr2d_init(eigbrd,eigi2d,dtset%mband,hdr0%nsppol,nkpt_rbz,dtset%natom)
#ifdef HAVE_TRIO_NETCDF
       NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating EIGI2D file")
       NCF_CHECK(crystal_ncwrite(Crystal, ncid))
       NCF_CHECK(ebands_ncwrite(Bands, ncid))
       call eigr2d_ncwrite(eigi2d,dtset%qptn(:),dtset%wtq,ncid)
       NCF_CHECK(nf90_close(ncid))
#else
       ABI_UNUSED(ncid)
#endif
     end if !smdelta
   end if 
 end if

 if (allocated(fan)) then
   ABI_DEALLOCATE(fan)
 end if
 if (allocated(eig2nkq_tmp)) then
   ABI_DEALLOCATE(eig2nkq_tmp)
 end if
 if (allocated(gkk)) then
   ABI_DEALLOCATE(gkk)
 end if

 call crystal_free(Crystal)
 call ebands_free(Bands)
 call eigr2d_free(eigr2d)
 call eigr2d_free(eigi2d)
 call fan_free(fan2d)
 call gkk_free(gkk2d)


 call timab(148,2,tsec)
!DEBUG
!write(std_out,*)' eig2tot: exit'
!ENDDEBUG

end subroutine eig2tot
!!***

