!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_mkrho
!! NAME
!!  m_mkrho
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR, LSI, AR, MB, MT)
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

module m_mkrho

 use defs_basis
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_xmpi
 use m_errors

 use m_time,         only : timab
 use m_fftcore,      only : sphereboundary
 use m_fft,          only : fftpac, zerosym
 use m_hamiltonian,  only : gs_hamiltonian_type
 use m_bandfft_kpt,  only : bandfft_kpt_set_ikpt
 use m_paw_dmft,     only : paw_dmft_type
 use m_spacepar,     only : symrhg
 use defs_datatypes, only : pseudopotential_type
 use m_atomdata,     only : atom_length
 use m_mpinfo,       only : ptabs_fourdp
 use m_pawtab,       only : pawtab_type

 implicit none

 private
!!***

 public :: mkrho
 public :: initro
!!***

contains
!!***

!!****f* m_mkrho/mkrho
!! NAME
!! mkrho
!!
!! FUNCTION
!! Depending on option argument value:
!! --Compute charge density rho(r) and rho(G) in electrons/bohr**3
!!   from input wavefunctions, band occupations, and k point wts.
!! --Compute kinetic energy density tau(r) and tau(G) in bohr**-5
!!   from input wavefunctions, band occupations, and k point wts.
!! --Compute a given element of the kinetic energy density tensor
!!   tau_{alpha,beta}(r) and tau_{alpha,beta}(G) in bohr**-5
!!   from input wavefunctions, band occupations, and k point wts.
!!
!! INPUTS
!!  cg(2,mcg)=wf in G space
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | istwfk(nkpt)=input option parameter that describes the storage of wfs
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs
!!   | mkmem=Number of k points treated by this node
!!   | mpw=maximum allowed value for npw
!!   | nband(nkpt*nsppol)=number of bands to be included in summation
!!   |  at each k point for each spin channel
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   | ngfft(18)=contain all needed information about 3D FFT,
!!   |  see ~abinit/doc/variables/vargs.htm#ngfft
!!   | nkpt=number of k points
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in group (at least 1 for identity)
!!   | symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!   | symrel(3,3,nsym)=symmetry matrices in real space (integers)
!!   | wtk(nkpt)=k point weights (they sum to 1.0)
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  kg(3,mpw*mkmem)=reduced planewave coordinates
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mpi_enreg=informations about MPI parallelization
!!  npwarr(nkpt)=number of planewaves and boundary planewaves at each k
!!  occ(mband*nkpt*nsppol)=
!!          occupation numbers for each band (usually 2.0) at each k point
!!  option if 0: compute rhor (electron density)
!!         if 1: compute taur (kinetic energy density)
!!               (i.e. Trace over the kinetic energy density tensor)
!!         if 2: compute taur_{alpha,beta} !!NOT YET IMPLEMENTED
!!               (a given element of the kinetic energy density tensor)
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  rprimd(3,3)=dimensional real space primitive translations
!!  tim_mkrho=timing code of the calling routine(can be set to 0 if not attributed)
!!  ucvol=unit cell volume (Bohr**3)
!!  wvl_den <type(wvl_denspot_type)>=density informations for wavelets
!!  wvl_wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets
!!
!! OUTPUT
!! rhog(2,nfft)=total electron density in G space
!! rhor(nfft,nspden)=electron density in r space
!!   (if spin polarized, array contains total density in first half and spin-up density in second half)
!!   (for non-collinear magnetism, first element: total density, 3 next ones: mx,my,mz in units of hbar/2)
!!
!! PARENTS
!!      afterscfloop,energy,gstate,respfn,scfcv,vtorho
!!
!! CHILDREN
!!      bandfft_kpt_set_ikpt,fftpac,fourwf,prep_fourwf,prtrhomxmn
!!      sphereboundary,symrhg,timab,wrtout,wvl_mkrho,xmpi_sum
!!
!! SOURCE

subroutine mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,npwarr,occ,paw_dmft,phnons,&
&                rhog,rhor,rprimd,tim_mkrho,ucvol,wvl_den,wvl_wfs,&
&                option) !optional

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkrho'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_53_ffts
 use interfaces_66_wfs
 use interfaces_67_common, except_this_one => mkrho
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mcg,tim_mkrho
 integer,intent(in),optional :: option
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(paw_dmft_type), intent(in)  :: paw_dmft
 type(wvl_wf_type),intent(inout) :: wvl_wfs
 type(wvl_denspot_type), intent(inout) :: wvl_den
!no_abirules
!nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
 integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,  &
&               (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem),npwarr(dtset%nkpt)
 real(dp), intent(in) :: gprimd(3,3)
 real(dp), intent(in) :: cg(2,mcg)
 real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
!nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
 real(dp), intent(in) :: phnons(2,(dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3))**(1-1/dtset%nsym),  &
&                                 (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 real(dp), intent(in) :: rprimd(3,3)
 real(dp), intent(out) :: rhor(dtset%nfft,dtset%nspden),rhog(2,dtset%nfft)

!Local variables-------------------------------
!scalars
 integer,save :: nskip=0
 integer :: alpha,use_nondiag_occup_dmft,bdtot_index,beta,blocksize,iband,iband1,ibandc1,ib,iblock,icg,ierr
 integer :: ifft,ikg,ikpt,ioption,ipw,ipwsp,ishf,ispden,ispinor,ispinor_index
 integer :: isppol,istwf_k,jspinor_index
 integer :: me,my_nspinor,n1,n2,n3,n4,n5,n6,nalpha,nband_k,nbandc1,nbdblock,nbeta
 integer :: ndat,nfftot,npw_k,spaceComm,tim_fourwf
 real(dp) :: kpt_cart,kg_k_cart,gp2pi1,gp2pi2,gp2pi3,cwftmp
 real(dp) :: weight,weight_i
 character(len=500) :: message
!arrays
 integer,allocatable :: gbound(:,:),kg_k(:,:)
 logical :: locc_test,nspinor1TreatedByThisProc,nspinor2TreatedByThisProc
 real(dp) :: dummy(2,1),tsec(2)
 real(dp),allocatable :: cwavef(:,:,:),cwavefb(:,:,:),cwavef_x(:,:)
 real(dp),allocatable :: cwavef_y(:,:),cwavefb_x(:,:),cwavefb_y(:,:),kg_k_cart_block(:)
 real(dp),allocatable :: occ_k(:),rhoaug(:,:,:),rhoaug_down(:,:,:)
 real(dp),allocatable :: rhoaug_mx(:,:,:),rhoaug_my(:,:,:),rhoaug_up(:,:,:)
 real(dp),allocatable :: taur_alphabeta(:,:,:,:),wfraug(:,:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(790+tim_mkrho,1,tsec)
 call timab(799,1,tsec)

 if(.not.(present(option))) then
   ioption=0
 else
   ioption=option
 end if

 if(ioption/=0.and.paw_dmft%use_sc_dmft==1) then
   MSG_ERROR('option argument value of this routines should be 0 if usedmft=1.')
 end if
 if(paw_dmft%use_sc_dmft/=0) then
   nbandc1=(paw_dmft%mbandc-1)*paw_dmft%use_sc_dmft+1
 else
   nbandc1=1
 end if
 use_nondiag_occup_dmft=0

!if(dtset%nspinor==2.and.paw_dmft%use_sc_dmft==1) then
!write(message, '(a,a,a,a)' )ch10,&
!&   ' mkrho : ERROR -',ch10,&
!&   '  nspinor argument value of this routines should be 1 if usedmft=1. '
!call wrtout(std_out,message,'COLL')
!end if

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 if (mpi_enreg%paral_spinor==0) then
   ispinor_index=1;jspinor_index=1
   nspinor1TreatedByThisProc=.true.
   nspinor2TreatedByThisProc=(dtset%nspinor==2)
 else
   ispinor_index=mpi_enreg%me_spinor+1;jspinor_index=3-ispinor_index
   nspinor1TreatedByThisProc=(mpi_enreg%me_spinor==0)
   nspinor2TreatedByThisProc=(mpi_enreg%me_spinor==1)
 end if

!Set local variable which depend on option argument

!nalpha*nbeta is the number of element of the kinetic energy density tensor
!to be computed in the irreducible Brillouin Zone (BZ) to get the result in the full BZ.
!In case of electron density calculation, nalpha=nbeta=1
 select case (ioption)
 case (0)
   nalpha = 1
   nbeta = 1
 case (1)
   nalpha = 3
   nbeta = 1
   ABI_ALLOCATE(taur_alphabeta,(dtset%nfft,dtset%nspden,3,1))
 case (2)
   nalpha = 3
   nbeta = 3
   ABI_ALLOCATE(taur_alphabeta,(dtset%nfft,dtset%nspden,3,3))
 case default
   MSG_BUG('ioption argument value should be 0,1 or 2.')
 end select

!Init me
 me=mpi_enreg%me_kpt
!zero the charge density array in real space
!$OMP PARALLEL DO COLLAPSE(2)
 do ispden=1,dtset%nspden
   do ifft=1,dtset%nfft
     rhor(ifft,ispden)=zero
   end do
 end do

!WVL - Branching with a separate mkrho procedure in wavelet.
 if (dtset%usewvl == 1) then
   select case(ioption)
   case (0)
     call wvl_mkrho(dtset, irrzon, mpi_enreg, phnons, rhor, wvl_wfs, wvl_den)
     return
   case (1)
     !call wvl_mkrho(dtset, mpi_enreg, occ, rhor, wvl_wfs, wvl_den)
     MSG_ERROR("kinetic energy density (taur) is not yet implemented in wavelet formalism.")
   case (2)
     !call wvl_mkrho(dtset, mpi_enreg, occ, rhor, wvl_wfs, wvl_den)
     MSG_BUG('kinetic energy density tensor (taur_(alpha,beta)) is not yet implemented in wavelet formalism.')
   end select
 end if
!WVL - Following is done in plane waves.

!start loop over alpha and beta

 do alpha=1,nalpha
   do beta=1,nbeta

!    start loop over spin and k points
     bdtot_index=0
     icg=0

!    n4,n5,n6 are FFT dimensions, modified to avoir cache trashing
     n1 = dtset%ngfft(1) ; n2 = dtset%ngfft(2) ; n3 = dtset%ngfft(3)
     n4 = dtset%ngfft(4) ; n5 = dtset%ngfft(5) ; n6 = dtset%ngfft(6)
     ndat = 1 ; if (mpi_enreg%paral_kgb==1) ndat = mpi_enreg%bandpp
     ABI_ALLOCATE(cwavef,(2,dtset%mpw,my_nspinor))
     ABI_ALLOCATE(rhoaug,(n4,n5,n6))
     ABI_ALLOCATE(wfraug,(2,n4,n5,n6*ndat))
     ABI_ALLOCATE(cwavefb,(2,dtset%mpw*paw_dmft%use_sc_dmft,my_nspinor))
     if(dtset%nspden==4) then
       ABI_ALLOCATE(rhoaug_up,(n4,n5,n6))
       ABI_ALLOCATE(rhoaug_down,(n4,n5,n6))
       ABI_ALLOCATE(rhoaug_mx,(n4,n5,n6))
       ABI_ALLOCATE(rhoaug_my,(n4,n5,n6))
       rhoaug_up(:,:,:)=zero
       rhoaug_down(:,:,:)=zero
       rhoaug_mx(:,:,:)=zero
       rhoaug_my(:,:,:)=zero
     end if

     do isppol=1,dtset%nsppol
       ikg=0

       rhoaug(:,:,:)=zero
       do ikpt=1,dtset%nkpt

         nband_k = dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
         npw_k=npwarr(ikpt)
         istwf_k = dtset%istwfk(ikpt)

         if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
           bdtot_index=bdtot_index+nband_k
           cycle
         end if

         ABI_ALLOCATE(gbound,(2*dtset%mgfft+8,2))
         ABI_ALLOCATE(kg_k,(3,npw_k))

         kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
         call sphereboundary(gbound,istwf_k,kg_k,dtset%mgfft,npw_k)

!        Loop over bands to fft and square for rho(r)
!        Shoulb be changed to treat bands by batch always

         if(mpi_enreg%paral_kgb /= 1) then  ! Not yet parallelized on spinors
           do iband=1,nband_k
!            if(paw_dmft%use_sc_dmft==1) then
!            write(std_out,*) 'iband  ',iband,occ(iband+bdtot_index),paw_dmft%occnd(iband,iband,ikpt,isppol)
!            else
!            write(std_out,*) 'iband  ',iband,occ(iband+bdtot_index)
!            endif
             if(mpi_enreg%paralbd==1)then
               if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband,iband,isppol,me)) cycle
             end if
             do ibandc1=1,nbandc1 ! in case of DMFT
!              Check if DMFT and only treat occupied states (check on occ.)
               if(paw_dmft%use_sc_dmft == 1) then
                 iband1 = paw_dmft%include_bands(ibandc1)
                 if(paw_dmft%band_in(iband)) then
                   if(.not. paw_dmft%band_in(iband1))  stop
                   use_nondiag_occup_dmft = 1
                   locc_test = abs(paw_dmft%occnd(1,iband,iband1,ikpt,isppol)) +&
&                   abs(paw_dmft%occnd(2,iband,iband1,ikpt,isppol))>tol8
!                  write(std_out,*) "mkrho,ikpt,iband,use_occnd",ikpt,iband
                 else
                   use_nondiag_occup_dmft = 0
                   locc_test = abs(occ(iband+bdtot_index))>tol8
                   if(ibandc1 /=1 .and. .not. paw_dmft%band_in(iband)) cycle
                 end if
               else
                 use_nondiag_occup_dmft = 0
                 locc_test = abs(occ(iband+bdtot_index))>tol8
               end if

               if (locc_test) then
!                Obtain Fourier transform in fft box and accumulate the density or the kinetic energy density
!                Not yet parallise on nspinor if paral_kgb non equal to 1
                 ipwsp=(iband-1)*npw_k*my_nspinor +icg
                 cwavef(:,1:npw_k,1)=cg(:,1+ipwsp:ipwsp+npw_k)
                 if (my_nspinor==2) cwavef(:,1:npw_k,2)=cg(:,ipwsp+npw_k+1:ipwsp+2*npw_k)

                 if(ioption==1)then
!                  Multiplication by 2pi i (k+G)_alpha
                   gp2pi1=gprimd(alpha,1)*two_pi ; gp2pi2=gprimd(alpha,2)*two_pi ; gp2pi3=gprimd(alpha,3)*two_pi
                   kpt_cart=gp2pi1*dtset%kptns(1,ikpt)+gp2pi2*dtset%kptns(2,ikpt)+gp2pi3*dtset%kptns(3,ikpt)
                   do ispinor=1,my_nspinor
                     do ipw=1,npw_k
                       kg_k_cart=gp2pi1*kg_k(1,ipw)+gp2pi2*kg_k(2,ipw)+gp2pi3*kg_k(3,ipw)+kpt_cart
                       cwftmp=-cwavef(2,ipw,ispinor)*kg_k_cart
                       cwavef(2,ipw,ispinor)=cwavef(1,ipw,ispinor)*kg_k_cart
                       cwavef(1,ipw,ispinor)=cwftmp
                     end do
                   end do
                 else if(ioption==2)then
                   MSG_ERROR('kinetic energy density tensor (taur_(alpha,beta)) is not yet implemented.')
                 end if

!                Non diag occupation in DMFT.
                 if(use_nondiag_occup_dmft==1) then
                   ipwsp=(iband1-1)*npw_k*my_nspinor +icg
                   cwavefb(:,1:npw_k,1)=cg(:,1+ipwsp:ipwsp+npw_k)
                   if (my_nspinor==2) cwavefb(:,1:npw_k,2)=cg(:,ipwsp+npw_k+1:ipwsp+2*npw_k)
                   weight  =paw_dmft%occnd(1,iband,iband1,ikpt,isppol)*dtset%wtk(ikpt)/ucvol
                   if(dtset%nspinor==1) weight_i=zero
                   if(dtset%nspinor==2) weight_i=paw_dmft%occnd(2,iband,iband1,ikpt,isppol)*dtset%wtk(ikpt)/ucvol

                 else
                   weight=occ(iband+bdtot_index)*dtset%wtk(ikpt)/ucvol
                   weight_i=weight
                 end if

                 if(mpi_enreg%paralbd==0) tim_fourwf=3
                 if(mpi_enreg%paralbd==1) tim_fourwf=6

!                The same section of code is also found in vtowfk.F90 : should be rationalized !

                 call fourwf(1,rhoaug,cwavef(:,:,1),dummy,wfraug,gbound,gbound,&
&                 istwf_k,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,&
&                 npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight,weight_i,&
&                 use_ndo=use_nondiag_occup_dmft,fofginb=cwavefb(:,:,1),&
&                 use_gpu_cuda=dtset%use_gpu_cuda)


                 if(dtset%nspinor==2)then
!                  DEBUG GZ !To obtain a x-directed magnetization(test)
!                  cwavef1(1,1:npw_k)=-cwavef(2,1:npw_k)
!                  cwavef1(2,1:npw_k)= cwavef(1,1:npw_k)
!                  ENDDEBUG

                   if(dtset%nspden==1) then
!                    We need only the total density : accumulation continues on top of rhoaug

                     call fourwf(1,rhoaug,cwavef(:,:,2),dummy,wfraug,gbound,gbound,&
&                     istwf_k,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,&
&                     npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight,weight_i,&
&                     use_ndo=use_nondiag_occup_dmft,fofginb=cwavefb(:,:,2),use_gpu_cuda=dtset%use_gpu_cuda)


                   else if(dtset%nspden==4) then
!                    Build the four components of rho. We use only norm quantities and, so fourwf.
!$\sum_{n} f_n \Psi^{* \alpha}_n \Psi^{\alpha}_n =\rho^{\alpha \alpha}$
!$\sum_{n} f_n (\Psi^{1}+\Psi^{2})^*_n (\Psi^{1}+\Psi^{2})_n=rho+m_x$
!$\sum_{n} f_n (\Psi^{1}-i \Psi^{2})^*_n (\Psi^{1}-i \Psi^{2})_n=rho+m_y$
                     ABI_ALLOCATE(cwavef_x,(2,npw_k))
                     ABI_ALLOCATE(cwavef_y,(2,npw_k))
                     ABI_ALLOCATE(cwavefb_x,(2,npw_k*paw_dmft%use_sc_dmft))
                     ABI_ALLOCATE(cwavefb_y,(2,npw_k*paw_dmft%use_sc_dmft))
!$(\Psi^{1}+\Psi^{2})$
                     cwavef_x(:,:)=cwavef(:,1:npw_k,1)+cwavef(:,1:npw_k,2)
!$(\Psi^{1}-i \Psi^{2})$
                     cwavef_y(1,:)=cwavef(1,1:npw_k,1)+cwavef(2,1:npw_k,2)
                     cwavef_y(2,:)=cwavef(2,1:npw_k,1)-cwavef(1,1:npw_k,2)
                     if(use_nondiag_occup_dmft==1) then
                       cwavefb_x(:,:)=cwavefb(:,1:npw_k,1)+cwavefb(:,1:npw_k,2)
                       cwavefb_y(1,:)=cwavefb(1,1:npw_k,1)+cwavefb(2,1:npw_k,2)
                       cwavefb_y(2,:)=cwavefb(2,1:npw_k,1)-cwavefb(1,1:npw_k,2)
                     end if
                     rhoaug_up(:,:,:)=rhoaug(:,:,:) !Already computed

                     call fourwf(1,rhoaug_down,cwavef(:,:,2),dummy,wfraug,gbound,gbound,&
&                     istwf_k,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,&
&                     npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight,weight_i,&
&                     use_ndo=use_nondiag_occup_dmft,fofginb=cwavefb(:,:,2),use_gpu_cuda=dtset%use_gpu_cuda)

                     call fourwf(1,rhoaug_mx,cwavef_x,dummy,wfraug,gbound,gbound,&
&                     istwf_k,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,&
&                     npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight,weight_i,&
&                     use_ndo=use_nondiag_occup_dmft,fofginb=cwavefb_x,use_gpu_cuda=dtset%use_gpu_cuda)

                     call fourwf(1,rhoaug_my,cwavef_y,dummy,wfraug,gbound,gbound,&
&                     istwf_k,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,&
&                     npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight,weight_i,&
&                     use_ndo=use_nondiag_occup_dmft,fofginb=cwavefb_y,use_gpu_cuda=dtset%use_gpu_cuda)

                     ABI_DEALLOCATE(cwavef_x)
                     ABI_DEALLOCATE(cwavef_y)
                     ABI_DEALLOCATE(cwavefb_x)
                     ABI_DEALLOCATE(cwavefb_y)

                   end if ! dtset%nspden/=4
                 end if

               else
!                Accumulate the number of one-way 3D ffts skipped
                 nskip=nskip+1
               end if ! abs(occ(iband+bdtot_index))>tol8
!              End loop on iband
             end do ! iband1=1,(nband_k-1)*paw_dmft%use_sc_dmft+1
           end do ! iband=1,nband_k

         else !paral_kgb==1
           if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) cycle

           call bandfft_kpt_set_ikpt(ikpt,mpi_enreg)

           nbdblock=nband_k/(mpi_enreg%nproc_band * mpi_enreg%bandpp)
           blocksize=nband_k/nbdblock
           if(allocated(cwavef))  then
             ABI_DEALLOCATE(cwavef)
           end if
           ABI_ALLOCATE(cwavef,(2,npw_k*blocksize,dtset%nspinor))
           if(ioption==1)  then
             ABI_ALLOCATE(kg_k_cart_block,(npw_k))
           end if
           ABI_ALLOCATE(occ_k,(nband_k))
           occ_k(:)=occ(bdtot_index+1:bdtot_index+nband_k)

           do iblock=1,nbdblock
             if (dtset%nspinor==1) then
               cwavef(:,1:npw_k*blocksize,1)=cg(:,1+(iblock-1)*npw_k*blocksize+icg:iblock*npw_k*blocksize+icg)
             else
               if (mpi_enreg%paral_spinor==0) then
                 ishf=(iblock-1)*npw_k*my_nspinor*blocksize+icg
                 do ib=1,blocksize
                   cwavef(:,(ib-1)*npw_k+1:ib*npw_k,1)=cg(:,1+(2*ib-2)*npw_k+ishf:(2*ib-1)*npw_k+ishf)
                   cwavef(:,(ib-1)*npw_k+1:ib*npw_k,2)=cg(:,1+(2*ib-1)*npw_k+ishf:ib*2*npw_k+ishf)
                 end do
               else
                 ishf=(iblock-1)*npw_k*my_nspinor*blocksize+icg
                 do ib=1,blocksize
                   cwavef(:,(ib-1)*npw_k+1:ib*npw_k,ispinor_index)=&
&                   cg(:,1+(ib-1)*npw_k+ishf:ib*npw_k+ishf)
                   cwavef(:,(ib-1)*npw_k+1:ib*npw_k,jspinor_index)=zero
                 end do
                 call xmpi_sum(cwavef,mpi_enreg%comm_spinor,ierr)
               end if
             end if
             if(ioption==1)then
!              Multiplication by 2pi i (k+G)_alpha
               gp2pi1=gprimd(alpha,1)*two_pi ; gp2pi2=gprimd(alpha,2)*two_pi ; gp2pi3=gprimd(alpha,3)*two_pi
               kpt_cart=gp2pi1*dtset%kptns(1,ikpt)+gp2pi2*dtset%kptns(2,ikpt)+gp2pi3*dtset%kptns(3,ikpt)
               kg_k_cart_block(1:npw_k)=gp2pi1*kg_k(1,1:npw_k)+gp2pi2*kg_k(2,1:npw_k)+gp2pi3*kg_k(3,1:npw_k)+kpt_cart
               do ib=1,blocksize
                 do ipw=1,npw_k
                   cwftmp=-cwavef(2,ipw+(ib-1)*npw_k,1)*kg_k_cart_block(ipw)
                   cwavef(2,ipw,1)=cwavef(1,ipw+(ib-1)*npw_k,1)*kg_k_cart_block(ipw)
                   cwavef(1,ipw,1)=cwftmp
                   if (my_nspinor==2) then
                     cwftmp=-cwavef(2,ipw+(ib-1)*npw_k,2)*kg_k_cart_block(ipw)
                     cwavef(2,ipw,2)=cwavef(1,ipw+(ib-1)*npw_k,2)*kg_k_cart_block(ipw)
                     cwavef(1,ipw,2)=cwftmp
                   end if
                 end do
               end do
             else if(ioption==2)then
               MSG_ERROR("kinetic energy density tensor (taur_(alpha,beta)) is not yet implemented.")
             end if

             call timab(538,1,tsec)
             if (nspinor1TreatedByThisProc) then
               call prep_fourwf(rhoaug,blocksize,cwavef(:,:,1),wfraug,iblock,istwf_k,dtset%mgfft,mpi_enreg,&
&               nband_k,ndat,dtset%ngfft,npw_k,n4,n5,n6,occ_k,1,ucvol,&
&               dtset%wtk(ikpt),use_gpu_cuda=dtset%use_gpu_cuda)
             end if
             call timab(538,2,tsec)
             if(dtset%nspinor==2)then
               if (dtset%nspden==1) then
                 if (nspinor2TreatedByThisProc) then
                   call prep_fourwf(rhoaug,blocksize,cwavef(:,:,2),wfraug,&
&                   iblock,istwf_k,dtset%mgfft,mpi_enreg,&
&                   nband_k,ndat,dtset%ngfft,npw_k,n4,n5,n6,occ_k,1,ucvol,&
&                   dtset%wtk(ikpt),use_gpu_cuda=dtset%use_gpu_cuda)
                 end if
               else if(dtset%nspden==4 ) then
                 ABI_ALLOCATE(cwavef_x,(2,npw_k*blocksize))
                 ABI_ALLOCATE(cwavef_y,(2,npw_k*blocksize))
                 cwavef_x(:,:)=cwavef(:,:,1)+cwavef(:,:,2)
                 cwavef_y(1,:)=cwavef(1,:,1)+cwavef(2,:,2)
                 cwavef_y(2,:)=cwavef(2,:,1)-cwavef(1,:,2)
                 call timab(538,1,tsec)
                 if (nspinor1TreatedByThisProc) then
                   call prep_fourwf(rhoaug_down,blocksize,cwavef(:,:,2),wfraug,&
&                   iblock,istwf_k,dtset%mgfft,mpi_enreg,&
&                   nband_k,ndat,dtset%ngfft,npw_k,n4,n5,n6,occ_k,1,ucvol,&
&                   dtset%wtk(ikpt),use_gpu_cuda=dtset%use_gpu_cuda)
                 end if
                 if (nspinor2TreatedByThisProc) then
                   call prep_fourwf(rhoaug_mx,blocksize,cwavef_x,wfraug,&
&                   iblock,istwf_k,dtset%mgfft,mpi_enreg,&
&                   nband_k,ndat,dtset%ngfft,npw_k,n4,n5,n6,occ_k,1,ucvol,&
&                   dtset%wtk(ikpt),use_gpu_cuda=dtset%use_gpu_cuda)
                   call prep_fourwf(rhoaug_my,blocksize,cwavef_y,wfraug,&
&                   iblock,istwf_k,dtset%mgfft,mpi_enreg,&
&                   nband_k,ndat,dtset%ngfft,npw_k,n4,n5,n6,occ_k,1,ucvol,&
&                   dtset%wtk(ikpt),use_gpu_cuda=dtset%use_gpu_cuda)
                 end if
                 call timab(538,2,tsec)
                 ABI_DEALLOCATE(cwavef_x)
                 ABI_DEALLOCATE(cwavef_y)
               end if
             end if
           end do !iblock
           if(ioption==1)  then
             ABI_DEALLOCATE(kg_k_cart_block)
           end if
           if (allocated(cwavef))  then
             ABI_DEALLOCATE(cwavef)
           end if
           ABI_DEALLOCATE(occ_k)
         end if

         ABI_DEALLOCATE(gbound)
         ABI_DEALLOCATE(kg_k)

         bdtot_index=bdtot_index+nband_k

         if (dtset%mkmem/=0) then
           icg=icg+npw_k*my_nspinor*nband_k
           ikg=ikg+npw_k
         end if

       end do ! ikpt

       if(mpi_enreg%paral_kgb == 1) then
         call bandfft_kpt_set_ikpt(-1,mpi_enreg)
         if (dtset%nspden==4) then
!          Sum the contribution of the band and of the FFT
           call xmpi_sum(rhoaug     ,mpi_enreg%comm_bandspinorfft, ierr)
           call xmpi_sum(rhoaug_down,mpi_enreg%comm_bandspinorfft, ierr)
           call xmpi_sum(rhoaug_mx ,mpi_enreg%comm_bandspinorfft, ierr)
           call xmpi_sum(rhoaug_my ,mpi_enreg%comm_bandspinorfft, ierr)
           rhoaug_up(:,:,:) = rhoaug(:,:,:)
         else
           call xmpi_sum(rhoaug,mpi_enreg%comm_bandspinorfft,ierr)
         end if
       end if

!      Transfer density on augmented fft grid to normal fft grid in real space
!      Take also into account the spin, to place it correctly in rhor.
       if(dtset%nspden==1 .or. dtset%nspden==2) then
         call fftpac(isppol,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug,1)
       else if(dtset%nspden==4) then
         ispden=1
         call fftpac(ispden,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug_up,1)
         ispden=2
         call fftpac(ispden,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug_mx,1)
         ispden=3
         call fftpac(ispden,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug_my,1)
         ispden=4
         call fftpac(ispden,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug_down,1)
         ABI_DEALLOCATE(rhoaug_up)
         ABI_DEALLOCATE(rhoaug_down)
         ABI_DEALLOCATE(rhoaug_mx)
         ABI_DEALLOCATE(rhoaug_my)
       end if

     end do ! isppol

     if(allocated(cwavef))  then
       ABI_DEALLOCATE(cwavef)
     end if
     ABI_DEALLOCATE(cwavefb)
     ABI_DEALLOCATE(rhoaug)
     ABI_DEALLOCATE(wfraug)

!    Recreate full rhor on all proc.
     call timab(48,1,tsec)
     call timab(71,1,tsec)
     spaceComm=mpi_enreg%comm_cell
     if (mpi_enreg%paral_hf==1)spaceComm=mpi_enreg%comm_kpt
     if(mpi_enreg%paral_kgb==1)spaceComm=mpi_enreg%comm_kpt
     call xmpi_sum(rhor,spaceComm,ierr)
     call timab(71,2,tsec)
     call timab(48,2,tsec)

     call timab(799,2,tsec)
     call timab(549,1,tsec)

     if(ioption==1 .or. ioption==2) then
!$OMP PARALLEL DO COLLAPSE(2)
       do ispden=1,dtset%nspden
         do ifft=1,dtset%nfft
           taur_alphabeta(ifft,ispden,alpha,beta) = rhor(ifft,ispden)
         end do
       end do
     end if

   end do !  beta=1,nbeta
 end do !  alpha=1,nalpha

!Compute the trace over the kinetic energy density tensor. i.e. Sum of the 3 diagonal elements.
 if(ioption==1)then
!  zero rhor array in real space
   do ispden=1,dtset%nspden
!$OMP PARALLEL DO
     do ifft=1,dtset%nfft
       rhor(ifft,ispden)=zero
     end do
   end do
   do alpha = 1, nalpha
!$OMP PARALLEL DO COLLAPSE(2)
     do ispden=1,dtset%nspden
       do ifft=1,dtset%nfft
         rhor(ifft,ispden) = rhor(ifft,ispden) + taur_alphabeta(ifft,ispden,alpha,1)
       end do
     end do
   end do
 end if

 nfftot=dtset%ngfft(1) * dtset%ngfft(2) * dtset%ngfft(3)

 select case (ioption)
 case(0, 1)
   call symrhg(1,gprimd,irrzon,mpi_enreg,dtset%nfft,nfftot,dtset%ngfft,dtset%nspden,dtset%nsppol,dtset%nsym,&
   dtset%paral_kgb,phnons,rhog,rhor,rprimd,dtset%symafm,dtset%symrel)
   if(ioption==1)then
!$OMP PARALLEL DO
     do ifft=1,dtset%nfft
       do ispden=1,dtset%nspden
         rhor(ifft,ispden)=1.0d0/2.0d0*rhor(ifft,ispden)
       end do
       rhog(:,ifft)=1.0d0/2.0d0*rhog(:,ifft)
     end do
   end if
 case(2)
   MSG_BUG('kinetic energy density tensor (taur_(alpha,beta)) is not yet implemented.')
   !call symtaug(1,gprimd,irrzon,mpi_enreg,dtset%nfft,nfftot,dtset%ngfft,dtset%nspden,dtset%nsppol,dtset%nsym,&
   !dtset%paral_kgb,phnons,rhog,rhor,rprimd,dtset%symafm,dtset%symrel)
 end select

 call timab(549,2,tsec)

!We now have both rho(r) and rho(G), symmetrized, and if dtset%nsppol=2
!we also have the spin-up density, symmetrized, in rhor(:,2).
!In case of non collinear magnetism, we have rho,mx,my,mz. No symmetry is applied

 call timab(799,1,tsec)

 if(ioption==1 .or. ioption==2)  then
   ABI_DEALLOCATE(taur_alphabeta)
 end if

!Find and print minimum and maximum total electron density
!(or total kinetic energy density, or total element of kinetic energy density tensor) and locations
 call wrtout(std_out,'mkrho: echo density (plane-wave part only)','COLL')
 call prtrhomxmn(std_out,mpi_enreg,dtset%nfft,dtset%ngfft,dtset%nspden,1,rhor,optrhor=ioption,ucvol=ucvol)

 call timab(799,2,tsec)
 call timab(790+tim_mkrho,2,tsec)

 DBG_EXIT("COLL")

end subroutine mkrho
!!***

!!****f* m_mkrho/initro
!!
!! NAME
!! initro
!!
!! FUNCTION
!! Initialize the density using either:
!!  - a gaussian of adjustable decay length (norm-conserving psp)
!!  - PS atomic valence density from psp file (PAW or NC psps with valence change in the pp file)
!!
!! INPUTS
!! atindx(natom)=index table for atoms (see gstate.f)
!! densty(ntypat,4)=parameters for initialisation of the density of each atom type
!! gmet(3,3)=reciprocal space metric (Bohr**-2)
!! gsqcut=cutoff G**2 for included G s in fft box (larger sphere).
!! izero=if 1, unbalanced components of rho(g) have to be set to zero
!! mgfft=maximum size of 1D FFTs
!! mpi_enreg=informations about mpi parallelization
!! mqgrid=number of grid pts in q array for n^AT(q) spline.
!! natom=number of atoms in cell.
!! nattyp(ntypat)=number of atoms of each type in cell.
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! ntypat=number of types of atoms in cell.
!! nspden=number of spin-density components
!! psps<type(pseudopotential_type)>=variables related to pseudopotentials
!! pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! ph1d(2,3*(2*mgfft+1)*natom)=1-dim phase information for given atom coordinates.
!! qgrid(mqgrid)=q grid for spline atomic valence density n^AT(q) from 0 to qmax.
!! spinat(3,natom)=initial spin of each atom, in unit of hbar/2.
!! ucvol=unit cell volume (Bohr**3).
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!! zion(ntypat)=charge on each type of atom (real number)
!! znucl(ntypat)=atomic number, for each type of atom
!!
!! OUTPUT
!! rhog(2,nfft)=initialized total density in reciprocal space
!! rhor(nfft,nspden)=initialized total density in real space.
!!         as well as spin-up part if spin-polarized
!!
!! PARENTS
!!      gstate,setup_positron
!!
!! CHILDREN
!!      fourdp,ptabs_fourdp,wrtout,zerosym
!!
!! SOURCE

subroutine initro(atindx,densty,gmet,gsqcut,izero,mgfft,mpi_enreg,mqgrid,natom,nattyp,&
&  nfft,ngfft,nspden,ntypat,paral_kgb,psps,pawtab,ph1d,qgrid,rhog,rhor,spinat,ucvol,usepaw,zion,znucl)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initro'
 use interfaces_14_hidewrite
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: izero,mgfft,mqgrid,natom,nfft,nspden,ntypat,paral_kgb
 integer,intent(in) :: usepaw
 real(dp),intent(in) :: gsqcut,ucvol
 type(mpi_type),intent(in) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx(natom),nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: densty(ntypat,4),gmet(3,3),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: qgrid(mqgrid),spinat(3,natom),zion(ntypat)
 real(dp),intent(in) :: znucl(ntypat)
 real(dp),intent(out) :: rhog(2,nfft),rhor(nfft,nspden)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!The decay lengths should be optimized element by element, and even pseudopotential by pseudopotential.
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i3,ia,ia1,ia2,id1,id2,id3,ig1,ig2,ig3,ii,ispden
 integer :: itypat,jj,jtemp,me_fft,n1,n2,n3,nproc_fft
 real(dp),parameter :: tolfix=1.000000001_dp
 real(dp) :: aa,alf2pi2,bb,cc,cutoff,dd,diff,dq,dq2div6,dqm1,fact,fact0,gmag
 real(dp) :: gsquar,rhoat,sfi,sfr
 real(dp) :: xnorm
 character(len=500) :: message
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:),fftn3_distrib(:),ffti3_local(:)
 real(dp),allocatable :: length(:),spinat_indx(:,:),work(:)
 logical,allocatable :: use_gaussian(:)

! *************************************************************************

 if(nspden==4)then
   write(std_out,*)' initro : might work yet for nspden=4 (not checked)'
   write(std_out,*)' spinat',spinat(1:3,1:natom)
!  stop
 end if

 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)
 me_fft=ngfft(11)
 nproc_fft=ngfft(10)
 ABI_ALLOCATE(work,(nfft))
 ABI_ALLOCATE(spinat_indx,(3,natom))

 ! Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!Transfer the spinat array to an array in which the atoms have the proper order, type by type.
 do ia=1,natom
   spinat_indx(:,atindx(ia))=spinat(:,ia)
 end do

!Check whether the values of spinat are acceptable
 if(nspden==2)then
   ia1=1
   do itypat=1,ntypat
!    ia1,ia2 sets range of loop over atoms:
     ia2=ia1+nattyp(itypat)-1
     do ia=ia1,ia2
       if( sqrt(spinat_indx(1,ia)**2+spinat_indx(2,ia)**2+spinat_indx(3,ia)**2) &
&       > abs(zion(itypat))*(1.0_dp + epsilon(0.0_dp)) ) then
         write(message, '(a,a,a,a,i4,a,a,3es11.4,a,a,a,es11.4)' ) ch10,&
&         ' initro : WARNING - ',ch10,&
&         '  For type-ordered atom number ',ia,ch10,&
&         '  input spinat=',spinat_indx(:,ia),'  is larger, in magnitude,',ch10,&
&         '  than zion(ia)=',zion(itypat)
         call wrtout(std_out,message,'COLL')
         call wrtout(ab_out,message,'COLL')
       end if
     end do
     ia1=ia2+1
   end do
 end if

!Compute the decay length of each type of atom
 ABI_ALLOCATE(length,(ntypat))
 ABI_ALLOCATE(use_gaussian,(ntypat))
 jtemp=0
 do itypat=1,ntypat

   use_gaussian(itypat)=.true.
   if (usepaw==0) use_gaussian(itypat) = .not. psps%nctab(itypat)%has_tvale
   if (usepaw==1) use_gaussian(itypat)=(pawtab(itypat)%has_tvale==0)
   if (.not.use_gaussian(itypat)) jtemp=jtemp+1

   if (use_gaussian(itypat)) then
     length(itypat) = atom_length(densty(itypat,1),zion(itypat),znucl(itypat))
     write(message,'(a,i3,a,f12.4,a,a,a,f12.4,a,i3,a,es12.4,a)' )&
&     ' initro: for itypat=',itypat,', take decay length=',length(itypat),',',ch10,&
&     ' initro: indeed, coreel=',znucl(itypat)-zion(itypat),', nval=',int(zion(itypat)),' and densty=',densty(itypat,1),'.'
     call wrtout(std_out,message,'COLL')
   else
     write(message,"(a,i3,a)")' initro: for itypat=',itypat,", take pseudo charge density from pp file"
     call wrtout(std_out,message,"COLL")
   end if

 end do

 if (jtemp>0) then
   dq=(qgrid(mqgrid)-qgrid(1))/dble(mqgrid-1)
   dqm1=1.0_dp/dq
   dq2div6=dq**2/6.0_dp
 end if

 cutoff=gsqcut*tolfix
 xnorm=1.0_dp/ucvol

 id1=n1/2+2
 id2=n2/2+2
 id3=n3/2+2

 if(nspden /= 4) then

   do ispden=nspden,1,-1
!    This loop overs spins will actually be as follows :
!    ispden=2 for spin up
!    ispden=1 for total spin (also valid for non-spin-polarized calculations)
!    The reverse ispden order is chosen, in order to end up with
!    rhog containing the proper total density.

     rhog(:,:)=zero

     ia1=1
     do itypat=1,ntypat

       if (use_gaussian(itypat)) alf2pi2=(two_pi*length(itypat))**2

!      ia1,ia2 sets range of loop over atoms:
       ia2=ia1+nattyp(itypat)-1
       ii=0
       jtemp=0

       do i3=1,n3
         ig3=i3-(i3/id3)*n3-1
         do i2=1,n2
           ig2=i2-(i2/id2)*n2-1
           if (fftn2_distrib(i2)==me_fft) then
             do i1=1,n1

               ig1=i1-(i1/id1)*n1-1
               ii=ii+1
!              gsquar=gsq_ini(ig1,ig2,ig3)
               gsquar=dble(ig1*ig1)*gmet(1,1)+dble(ig2*ig2)*gmet(2,2)+&
&               dble(ig3*ig3)*gmet(3,3)+dble(2*ig1*ig2)*gmet(1,2)+&
&               dble(2*ig2*ig3)*gmet(2,3)+dble(2*ig3*ig1)*gmet(3,1)

!              Skip G**2 outside cutoff:
               if (gsquar<=cutoff) then

!                Assemble structure factor over all atoms of given type,
!                also taking into account the spin-charge on each atom:
                 sfr=zero;sfi=zero
                 if(ispden==1)then
                   do ia=ia1,ia2
                     sfr=sfr+phre_ini(ig1,ig2,ig3,ia)
                     sfi=sfi-phimag_ini(ig1,ig2,ig3,ia)
                   end do
                   if (use_gaussian(itypat)) then
                     sfr=sfr*zion(itypat)
                     sfi=sfi*zion(itypat)
                   end if
                 else
                   fact0=half;if (.not.use_gaussian(itypat)) fact0=half/zion(itypat)
                   do ia=ia1,ia2
!                    Here, take care only of the z component
                     fact=fact0*(zion(itypat)+spinat_indx(3,ia))
                     sfr=sfr+phre_ini(ig1,ig2,ig3,ia)*fact
                     sfi=sfi-phimag_ini(ig1,ig2,ig3,ia)*fact
                   end do
                 end if

!                Charge density integrating to one
                 if (use_gaussian(itypat)) then
                   rhoat=xnorm*exp(-gsquar*alf2pi2)
!                  Multiply structure factor times rhoat (atomic density in reciprocal space)
                   rhog(re,ii)=rhog(re,ii)+sfr*rhoat
                   rhog(im,ii)=rhog(im,ii)+sfi*rhoat
                 else
                   gmag=sqrt(gsquar)
                   jj=1+int(gmag*dqm1)
                   diff=gmag-qgrid(jj)
                   bb = diff*dqm1
                   aa = one-bb
                   cc = aa*(aa**2-one)*dq2div6
                   dd = bb*(bb**2-one)*dq2div6
                   if (usepaw == 1) then
                     rhoat=(aa*pawtab(itypat)%tvalespl(jj,1)+bb*pawtab(itypat)%tvalespl(jj+1,1)+&
&                     cc*pawtab(itypat)%tvalespl(jj,2)+dd*pawtab(itypat)%tvalespl(jj+1,2)) *xnorm
                   else if (usepaw == 0) then
                     rhoat=(aa*psps%nctab(itypat)%tvalespl(jj,1)+bb*psps%nctab(itypat)%tvalespl(jj+1,1)+&
                     cc*psps%nctab(itypat)%tvalespl(jj,2)+dd*psps%nctab(itypat)%tvalespl(jj+1,2))*xnorm
                   else
                     MSG_BUG('Initialization of density is non consistent.')
                   end if
!                  Multiply structure factor times rhoat (atomic density in reciprocal space)
                   rhog(re,ii)=rhog(re,ii)+sfr*rhoat
                   rhog(im,ii)=rhog(im,ii)+sfi*rhoat
                 end if

               else
                 jtemp=jtemp+1
               end if

             end do ! End loop on i1
           end if
         end do ! End loop on i2
       end do ! End loop on i3
       ia1=ia2+1

     end do ! End loop on type of atoms

!    Set contribution of unbalanced components to zero
     if (izero==1) then
       call zerosym(rhog,2,n1,n2,n3,comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
     end if
     !write(std_out,*)"initro: ispden, ucvol * rhog(:2,1)",ispden, ucvol * rhog(:2,1)

!    Note, we end with ispden=1, so that rhog contains the total density
     call fourdp(1,rhog,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
     rhor(:,ispden)=work(:)
   end do ! End loop on spins

 else if(nspden==4) then
   do ispden=nspden,1,-1
!    This loop overs spins will actually be as follows :
!    ispden=2,3,4 for mx,my,mz
!    ispden=1 for total spin (also valid for non-spin-polarized calculations)
!    The reverse ispden order is chosen, in order to end up with
!    rhog containing the proper total density.

     rhog(:,:)=zero

     ia1=1
     do itypat=1,ntypat

       if (use_gaussian(itypat)) alf2pi2=(two_pi*length(itypat))**2

!      ia1,ia2 sets range of loop over atoms:
       ia2=ia1+nattyp(itypat)-1
       ii=0
       jtemp=0
       do i3=1,n3
         ig3=i3-(i3/id3)*n3-1
         do i2=1,n2
           ig2=i2-(i2/id2)*n2-1
           if (fftn2_distrib(i2)==me_fft) then
             do i1=1,n1

               ig1=i1-(i1/id1)*n1-1
               ii=ii+1
!              gsquar=gsq_ini(ig1,ig2,ig3)
               gsquar=dble(ig1*ig1)*gmet(1,1)+dble(ig2*ig2)*gmet(2,2)+&
&               dble(ig3*ig3)*gmet(3,3)+dble(2*ig1*ig2)*gmet(1,2)+&
&               dble(2*ig2*ig3)*gmet(2,3)+dble(2*ig3*ig1)*gmet(3,1)

!              Skip G**2 outside cutoff:
               if (gsquar<=cutoff) then

!                Assemble structure factor over all atoms of given type,
!                also taking into account the spin-charge on each atom:
                 sfr=zero;sfi=zero
                 if(ispden==1)then
                   do ia=ia1,ia2
                     sfr=sfr+phre_ini(ig1,ig2,ig3,ia)
                     sfi=sfi-phimag_ini(ig1,ig2,ig3,ia)
                   end do
                   if (use_gaussian(itypat)) then
                     sfr=sfr*zion(itypat)
                     sfi=sfi*zion(itypat)
                   end if
                 else
                   fact0=one;if (.not.use_gaussian(itypat)) fact0=one/zion(itypat)
                   do ia=ia1,ia2
!                    Here, take care of the components of m
                     fact=fact0*spinat_indx(ispden-1,ia)
                     sfr=sfr+phre_ini(ig1,ig2,ig3,ia)*fact
                     sfi=sfi-phimag_ini(ig1,ig2,ig3,ia)*fact
                   end do
                 end if

!                Charge density integrating to one
                 if (use_gaussian(itypat)) then
                   rhoat=xnorm*exp(-gsquar*alf2pi2)
                 else
                   gmag=sqrt(gsquar)
                   jj=1+int(gmag*dqm1)
                   diff=gmag-qgrid(jj)
                   bb = diff*dqm1
                   aa = one-bb
                   cc = aa*(aa**2-one)*dq2div6
                   dd = bb*(bb**2-one)*dq2div6
                   if (usepaw == 1) then
                     rhoat=(aa*pawtab(itypat)%tvalespl(jj,1)+bb*pawtab(itypat)%tvalespl(jj+1,1)+&
&                     cc*pawtab(itypat)%tvalespl(jj,2)+dd*pawtab(itypat)%tvalespl(jj+1,2)) *xnorm
                   else if (usepaw == 0) then
                     rhoat=(aa*psps%nctab(itypat)%tvalespl(jj,1)+bb*psps%nctab(itypat)%tvalespl(jj+1,1)+&
                     cc*psps%nctab(itypat)%tvalespl(jj,2)+dd*psps%nctab(itypat)%tvalespl(jj+1,2))*xnorm
                   else
                     MSG_BUG('Initialization of density is non consistent.')
                   end if
                 end if

!                Multiply structure factor times rhoat (atomic density in reciprocal space)
                 rhog(re,ii)=rhog(re,ii)+sfr*rhoat
                 rhog(im,ii)=rhog(im,ii)+sfi*rhoat
               else
                 jtemp=jtemp+1
               end if

             end do ! End loop on i1
           end if
         end do ! End loop on i2
       end do ! End loop on i3
       ia1=ia2+1
     end do ! End loop on type of atoms

!    Set contribution of unbalanced components to zero
     if (izero==1) then
       call zerosym(rhog,2,n1,n2,n3,comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
     end if
     !write(std_out,*)"initro: ispden, ucvol * rhog(:2,1)",ispden, ucvol * rhog(:2,1)

!    Note, we end with ispden=1, so that rhog contains the total density
     call fourdp(1,rhog,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
     rhor(:,ispden)=work(:)

   end do ! End loop on spins

!  Non-collinear magnetism: avoid zero magnetization, because it produces numerical instabilities
!    Add a small real to the magnetization
   if (all(abs(spinat(:,:))<tol10)) rhor(:,4)=rhor(:,4)+tol14

 end if ! nspden==4

 ABI_DEALLOCATE(length)
 ABI_DEALLOCATE(use_gaussian)
 ABI_DEALLOCATE(spinat_indx)
 ABI_DEALLOCATE(work)

 contains

!Real and imaginary parts of phase.
   function phr_ini(x1,y1,x2,y2,x3,y3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phr_ini'
!End of the abilint section

   real(dp) :: phr_ini
   real(dp),intent(in) :: x1,x2,x3,y1,y2,y3
   phr_ini=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
 end function phr_ini

   function phi_ini(x1,y1,x2,y2,x3,y3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phi_ini'
!End of the abilint section

   real(dp) :: phi_ini
   real(dp),intent(in) :: x1,x2,x3,y1,y2,y3
   phi_ini=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
 end function phi_ini

   function ph1_ini(nri,ig1,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ph1_ini'
!End of the abilint section

   real(dp) :: ph1_ini
   integer,intent(in) :: nri,ig1,ia
   ph1_ini=ph1d(nri,ig1+1+n1+(ia-1)*(2*n1+1))
 end function ph1_ini

   function ph2_ini(nri,ig2,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ph2_ini'
!End of the abilint section

   real(dp) :: ph2_ini
   integer,intent(in) :: nri,ig2,ia
   ph2_ini=ph1d(nri,ig2+1+n2+(ia-1)*(2*n2+1)+natom*(2*n1+1))
 end function ph2_ini

   function ph3_ini(nri,ig3,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ph3_ini'
!End of the abilint section

   real(dp) :: ph3_ini
   integer,intent(in) :: nri,ig3,ia
   ph3_ini=ph1d(nri,ig3+1+n3+(ia-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1))
 end function ph3_ini

   function phre_ini(ig1,ig2,ig3,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phre_ini'
!End of the abilint section

   real(dp) :: phre_ini
   integer,intent(in) :: ig1,ig2,ig3,ia
   phre_ini=phr_ini(ph1_ini(re,ig1,ia),ph1_ini(im,ig1,ia),&
&   ph2_ini(re,ig2,ia),ph2_ini(im,ig2,ia),ph3_ini(re,ig3,ia),ph3_ini(im,ig3,ia))
 end function phre_ini

   function phimag_ini(ig1,ig2,ig3,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phimag_ini'
!End of the abilint section

   real(dp) :: phimag_ini
   integer,intent(in) :: ig1,ig2,ig3,ia
   phimag_ini=phi_ini(ph1_ini(re,ig1,ia),ph1_ini(im,ig1,ia),&
&   ph2_ini(re,ig2,ia),ph2_ini(im,ig2,ia),ph3_ini(re,ig3,ia),ph3_ini(im,ig3,ia))
 end function phimag_ini

end subroutine initro
!!***

end module m_mkrho
!!***
