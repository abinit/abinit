!{\src2tex{textfont=tt}}
!!****m* m_mkrho/m_mkrho
!! NAME
!!  m_mkrho
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA, XG, GMR, LSI, AR, MB, MT)
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
 use m_abicore
 use m_xmpi
 use m_errors

 use m_time,         only : timab
 use m_fftcore,      only : sphereboundary
 use m_fft,          only : fftpac, zerosym, fourwf, fourdp
 use m_hamiltonian,  only : gs_hamiltonian_type
 use m_bandfft_kpt,  only : bandfft_kpt_set_ikpt
 use m_paw_dmft,     only : paw_dmft_type
 use m_spacepar,     only : symrhg
 use defs_datatypes, only : pseudopotential_type
 use m_atomdata,     only : atom_length
 use m_mpinfo,       only : ptabs_fourdp, proc_distrb_cycle
 use m_pawtab,       only : pawtab_type
 use m_io_tools,     only : open_file
 use m_splines,      only : spline,splint
 use m_sort,         only : sort_dp
 use m_prep_kgb,     only : prep_fourwf
 use m_wvl_rho,      only : wvl_mkrho
 use m_rot_cg,       only : rot_cg

 implicit none

 private
!!***

 public :: mkrho
 public :: initro
 public :: prtrhomxmn
 public :: read_atomden
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
!!      sphereboundary,symrhg,timab,wrtout,wvl_mkrho,xmpi_sum,diag_occ_rot_cg
!!
!! SOURCE

subroutine mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,npwarr,occ,paw_dmft,phnons,&
&                rhog,rhor,rprimd,tim_mkrho,ucvol,wvl_den,wvl_wfs,&
&                option) !optional

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
 !character(len=500) :: message
!arrays
 integer,allocatable :: gbound(:,:),kg_k(:,:)
 logical :: locc_test,nspinor1TreatedByThisProc,nspinor2TreatedByThisProc
 real(dp) :: dummy(2,1),tsec(2)
 real(dp),allocatable :: cwavef(:,:,:),cwavefb(:,:,:),cwavef_x(:,:)
 real(dp),allocatable :: cwavef_y(:,:),cwavefb_x(:,:),cwavefb_y(:,:),kg_k_cart_block(:)
 real(dp),allocatable :: occ_k(:),rhoaug(:,:,:),rhoaug_down(:,:,:)
 real(dp),allocatable :: rhoaug_mx(:,:,:),rhoaug_my(:,:,:),rhoaug_up(:,:,:)
 real(dp),allocatable :: taur_alphabeta(:,:,:,:),wfraug(:,:,:,:)
 real(dp),allocatable :: occ_diag(:)
! real(dp),allocatable :: occ_nd(2, :, :)
 real(dp),allocatable :: cwavef_rot(:,:,:,:)

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
                   weight_i=paw_dmft%occnd(2,iband,iband1,ikpt,isppol)*dtset%wtk(ikpt)/ucvol

                 else
                   weight=occ(iband+bdtot_index)*dtset%wtk(ikpt)/ucvol
                   weight_i=weight
                 end if

                 if(mpi_enreg%paralbd==0) tim_fourwf=3
                 if(mpi_enreg%paralbd==1) tim_fourwf=6

!                The same section of code is also found in vtowfk.F90 : should be rationalized !

                 call fourwf(1,rhoaug,cwavef(:,:,1),dummy,wfraug,gbound,gbound,&
&                 istwf_k,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,&
&                 npw_k,1,n4,n5,n6,1,tim_fourwf,weight,weight_i,&
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
&                     npw_k,1,n4,n5,n6,1,tim_fourwf,weight,weight_i,&
&                     use_ndo=use_nondiag_occup_dmft,fofginb=cwavefb(:,:,2),use_gpu_cuda=dtset%use_gpu_cuda)


                   else if(dtset%nspden==4) then
                     ! Build the four components of rho. We use only norm quantities and, so fourwf.
                     ! $\sum_{n} f_n \Psi^{* \alpha}_n \Psi^{\alpha}_n =\rho^{\alpha \alpha}$
                     ! $\sum_{n} f_n (\Psi^{1}+\Psi^{2})^*_n (\Psi^{1}+\Psi^{2})_n=rho+m_x$
                     ! $\sum_{n} f_n (\Psi^{1}-i \Psi^{2})^*_n (\Psi^{1}-i \Psi^{2})_n=rho+m_y$
                     ABI_ALLOCATE(cwavef_x,(2,npw_k))
                     ABI_ALLOCATE(cwavef_y,(2,npw_k))
                     ABI_ALLOCATE(cwavefb_x,(2,npw_k*paw_dmft%use_sc_dmft))
                     ABI_ALLOCATE(cwavefb_y,(2,npw_k*paw_dmft%use_sc_dmft))
                     ! $(\Psi^{1}+\Psi^{2})$
                     cwavef_x(:,:)=cwavef(:,1:npw_k,1)+cwavef(:,1:npw_k,2)
                     ! $(\Psi^{1}-i \Psi^{2})$
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
&                     npw_k,1,n4,n5,n6,1,tim_fourwf,weight,weight_i,&
&                     use_ndo=use_nondiag_occup_dmft,fofginb=cwavefb(:,:,2),use_gpu_cuda=dtset%use_gpu_cuda)

                     call fourwf(1,rhoaug_mx,cwavef_x,dummy,wfraug,gbound,gbound,&
&                     istwf_k,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,&
&                     npw_k,1,n4,n5,n6,1,tim_fourwf,weight,weight_i,&
&                     use_ndo=use_nondiag_occup_dmft,fofginb=cwavefb_x,use_gpu_cuda=dtset%use_gpu_cuda)

                     call fourwf(1,rhoaug_my,cwavef_y,dummy,wfraug,gbound,gbound,&
&                     istwf_k,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,&
&                     npw_k,1,n4,n5,n6,1,tim_fourwf,weight,weight_i,&
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

! ---------- DMFT
           if(allocated(cwavef_rot))  then
             ABI_DEALLOCATE(cwavef_rot)
             ABI_DEALLOCATE(occ_diag)
             ! ABI_DEALLOCATE(occ_nd)
           end if
           if(paw_dmft%use_sc_dmft==1) then
             ! Allocation of DMFT temporaries arrays
             ABI_ALLOCATE(cwavef_rot,(2,npw_k,blocksize,dtset%nspinor))
             ABI_ALLOCATE(occ_diag,(blocksize))
             ! ABI_ALLOCATE(occ_nd,(2, blocksize, blocksize, dtset%nspinor))
           end if
! ---------- END DMFT

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

! ---------- DMFT
             if(paw_dmft%use_sc_dmft==1) then
               ! initialisation of DMFT arrays
               cwavef_rot(:,:,:,:) = zero
               occ_diag(:) = zero
               ! occ_nd(:,:,:,:) = paw_dmft%occnd(:,:,:,ikpt,:)

               do ib=1,blocksize
                 cwavef_rot(:, :, ib, :) = cwavef(:, 1+(ib-1)*npw_k:ib*npw_k, :)
               end do

               call rot_cg(paw_dmft%occnd(:,:,:,ikpt,isppol), cwavef_rot, npw_k, nband_k, blocksize,&
&                          dtset%nspinor, paw_dmft%include_bands(1), paw_dmft%mbandc, occ_diag)
               do ib=1,blocksize
                 cwavef(:, 1+(ib-1)*npw_k:ib*npw_k, :) = cwavef_rot(:, :, ib, :)
               end do

               occ_k(:) = occ_diag(:)
             end if
! ---------- END DMFT

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

     if(allocated(cwavef_rot))  then
       ABI_DEALLOCATE(cwavef_rot)
       ABI_DEALLOCATE(occ_diag)
       ! ABI_DEALLOCATE(occ_nd)
     end if

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
&   phnons,rhog,rhor,rprimd,dtset%symafm,dtset%symrel)
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
 call wrtout(std_out,' mkrho: echo density (plane-wave part only)','COLL')
 call prtrhomxmn(std_out,mpi_enreg,dtset%nfft,dtset%ngfft,dtset%nspden,1,rhor,optrhor=ioption,ucvol=ucvol)

 call timab(799,2,tsec)
 call timab(790+tim_mkrho,2,tsec)

 DBG_EXIT("COLL")

end subroutine mkrho
!!***

!!****f* m_mkrho/initro
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
&  nfft,ngfft,nspden,ntypat,psps,pawtab,ph1d,qgrid,rhog,rhor,spinat,ucvol,usepaw,zion,znucl)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: izero,mgfft,mqgrid,natom,nfft,nspden,ntypat
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
     call fourdp(1,rhog,work,1,mpi_enreg,nfft,1,ngfft,0)
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
     call fourdp(1,rhog,work,1,mpi_enreg,nfft,1,ngfft,0)
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

   real(dp) :: phr_ini
   real(dp),intent(in) :: x1,x2,x3,y1,y2,y3
   phr_ini=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
 end function phr_ini

   function phi_ini(x1,y1,x2,y2,x3,y3)

   real(dp) :: phi_ini
   real(dp),intent(in) :: x1,x2,x3,y1,y2,y3
   phi_ini=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
 end function phi_ini

   function ph1_ini(nri,ig1,ia)

   real(dp) :: ph1_ini
   integer,intent(in) :: nri,ig1,ia
   ph1_ini=ph1d(nri,ig1+1+n1+(ia-1)*(2*n1+1))
 end function ph1_ini

   function ph2_ini(nri,ig2,ia)

   real(dp) :: ph2_ini
   integer,intent(in) :: nri,ig2,ia
   ph2_ini=ph1d(nri,ig2+1+n2+(ia-1)*(2*n2+1)+natom*(2*n1+1))
 end function ph2_ini

   function ph3_ini(nri,ig3,ia)

   real(dp) :: ph3_ini
   integer,intent(in) :: nri,ig3,ia
   ph3_ini=ph1d(nri,ig3+1+n3+(ia-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1))
 end function ph3_ini

   function phre_ini(ig1,ig2,ig3,ia)

   real(dp) :: phre_ini
   integer,intent(in) :: ig1,ig2,ig3,ia
   phre_ini=phr_ini(ph1_ini(re,ig1,ia),ph1_ini(im,ig1,ia),&
&   ph2_ini(re,ig2,ia),ph2_ini(im,ig2,ia),ph3_ini(re,ig3,ia),ph3_ini(im,ig3,ia))
 end function phre_ini

   function phimag_ini(ig1,ig2,ig3,ia)

   real(dp) :: phimag_ini
   integer,intent(in) :: ig1,ig2,ig3,ia
   phimag_ini=phi_ini(ph1_ini(re,ig1,ia),ph1_ini(im,ig1,ia),&
&   ph2_ini(re,ig2,ia),ph2_ini(im,ig2,ia),ph3_ini(re,ig3,ia),ph3_ini(im,ig3,ia))
 end function phimag_ini

end subroutine initro
!!***

!!****f* m_mkrho/prtrhomxmn
!! NAME
!! prtrhomxmn
!!
!! FUNCTION
!! If option==1, compute the maximum and minimum of the density (and spin-polarization
!! if nspden==2), and print it.
!! If option==2, also compute and print the second maximum or minimum
!!
!! INPUTS
!!  iout=unit for output file
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  option, see above
!!  optrhor=option for rhor (If optrhor==0, rhor is expected to be the electron density)
!!                          (If optrhor==1, rhor is expected to be the kinetic energy density (taur))
!!                          (If optrhor==2, rhor is expected to be the gradient of the electron density (grhor))
!!                          (If optrhor==3, rhor is expected to be the laplacian of the electron density (lrhor))
!!                          (If optrhor==4, rhor is expected to be the ELF (elfr))
!!  rhor(nfft,nspden)=electron density (electrons/bohr^3)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  The tolerance tol12 aims at giving a machine-independent ordering.
!!  (this trick is used in bonds.f, listkk.f, prtrhomxmn.f and rsiaf9.f)
!!
!! PARENTS
!!      afterscfloop,bethe_salpeter,clnup1,mkrho,screening,sigma,vtorho
!!
!! CHILDREN
!!      wrtout,xmpi_sum
!!
!! SOURCE

subroutine prtrhomxmn(iout,mpi_enreg,nfft,ngfft,nspden,option,rhor,optrhor,ucvol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,nfft,nspden,option
 integer,intent(in),optional :: optrhor
 real(dp),intent(in),optional :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhor(nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,ierr,ifft,ii,iisign,iitems,index1,ioptrhor
 integer :: index2,indsign,iproc,istart,me,n1,n2,n3,nitems
 integer :: nfft_,nfftot,nproc,spaceComm
 real(dp) :: temp,value1,value2
 character(len=500) :: message,txt1_in_mssg,txt2_in_mssg,txt3_in_mssg
 logical :: reduce=.false.
!arrays
 integer,allocatable :: iindex(:,:,:),index_fft(:,:,:,:)
 real(dp) :: rhomn1(4),rhomn2(4),rhomx1(4),rhomx2(4),ri_rhomn1(3,4)
 real(dp) :: ri_rhomn2(3,4),ri_rhomx1(3,4),ri_rhomx2(3,4),ri_zetmn1(3,2)
 real(dp) :: ri_zetmn2(3,2),ri_zetmx1(3,2),ri_zetmx2(3,2),zetmn1(2)
 real(dp) :: zetmn2(2),zetmx1(2),zetmx2(2)
 real(dp),allocatable :: array(:),coord(:,:,:,:),value(:,:,:),integrated(:)
 real(dp),allocatable :: value_fft(:,:,:)

! *************************************************************************

 if(.not.(present(optrhor))) then
   ioptrhor=0
 else
   ioptrhor=optrhor
 end if

 if(option/=1 .and. option/=2)then
   write(message, '(a,i0)' )' Option must be 1 or 2, while it is ',option
   MSG_BUG(message)
 end if

 if (mpi_enreg%nproc_wvl>1) then
!  nfft is always the potential size (in GGA, the density has buffers).
   nfft_ = ngfft(1) * ngfft(2) * mpi_enreg%nscatterarr(mpi_enreg%me_wvl, 2)
   n1 = ngfft(1)
   n2 = ngfft(2)
   n3 = sum(mpi_enreg%nscatterarr(:, 2))
   istart = mpi_enreg%nscatterarr(mpi_enreg%me_wvl, 4)
 else
   nfft_ = nfft
   n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
   istart = 0
 end if

!--------------------------------------------------------------------------
!One has to determine the maximum and minimum (etc...) values
!over all space, and then output it, as well as to identify
!the point at which it occurs ...
!This will require a bit of data exchange, and correct indirect indexing ...

!For the local processor, find different items :
!maximum and minimum total electron density and locations
!and also spin-polarisation and magnetization
!also keep the second maximal or minimal value
 if(nspden==1)nitems=1   ! Simply the total density
 if(nspden==2)nitems=5   ! Total density, spin up, spin down, magnetization, zeta
 if(nspden==4)nitems=6   ! Total density, x, y, z, magnetization, zeta

 ABI_ALLOCATE(value,(2,2,nitems))
 ABI_ALLOCATE(iindex,(2,2,nitems))
 ABI_ALLOCATE(array,(nfft))
 ABI_ALLOCATE(integrated,(nitems))

 do iitems=1,nitems

!  Copy the correct values into the array
!  First set of items : the density, for each spin component
   if(iitems<=nspden)then
     array(:)=rhor(:,iitems)
   end if
!  Case nspden==2, some computation to be done
   if(nspden==2)then
     if(iitems==3)then ! Spin down
       array(:)=rhor(:,1)-rhor(:,2)
     else if(iitems==4)then  ! Magnetization
       array(:)=2*rhor(:,2)-rhor(:,1)
     else if(iitems==5)then  ! zeta = relative magnetization
       ! Avoid 0/0: the limit of (x - y) / (x+ y) depends on the direction.
       array(:)=zero
       where (abs(rhor(:,1)) > tol12) array(:)=(2*rhor(:,2)-rhor(:,1))/rhor(:,1)
     end if
!    Case nspden==4, some other computation to be done
   else if(nspden==4)then
     if(iitems==5)then ! Magnetization
       array(:)=sqrt(rhor(:,2)**2+rhor(:,3)**2+rhor(:,4)**2)
     else if(iitems==6)then ! zeta = relative magnetization
       array(:)=(sqrt(rhor(:,2)**2+rhor(:,3)**2+rhor(:,4)**2))/rhor(:,1)
     end if
   end if

!  Zero all the absolute values that are lower than tol8, for portability reasons.
   do ifft = 1, nfft_
     if(abs(array(ifft))<tol8)array(ifft)=zero
   end do

!  DEBUG
!  write(std_out,*) ' iitems,array(1:2)=',iitems,array(1:2)
!  ENDDEBUG

   do indsign=1,2 ! Find alternatively the maximum and the minimum
     iisign=3-2*indsign

     if (nfft_ > 1) then
!      Initialize the two first values
       value1=array(istart + 1) ; value2=array(istart + 2)
       index1=1 ; index2=2

!      Ordering, if needed
       if( iisign*(value2+tol12) > iisign*(value1)) then
         temp=value2 ; value2=value1 ; value1=temp
         index1=2 ; index2=1
       end if

!      Integration, if relevant
       if(present(ucvol).and. indsign==1)then
         integrated(iitems) = array(istart + 1)+array(istart + 2)
       end if
     else
       value1 = zero; value2 = zero
       index1 = 0;    index2 = 0
     end if

!    DEBUG
!    write(std_out,*) ' value1,value2,index1,index2=',value1,value2,index1,index2
!    ENDDEBUG

!    Loop over all points
     do ifft = 3, nfft_

       temp=array(istart + ifft)
       if(present(ucvol).and. indsign==1)integrated(iitems) = integrated(iitems)+temp
!      Compares it to the second value
       if( iisign*(temp+tol12) > iisign*value2 ) then
!        Compare it to the first value
         if( iisign*(temp+tol12) > iisign*value1 ) then
           value2=value1 ; index2=index1
           value1=temp   ; index1=ifft
         else
           value2=temp   ; index2=ifft
         end if
       end if

     end do ! ifft

     value(1,indsign,iitems)=value1
     value(2,indsign,iitems)=value2
     iindex(1,indsign,iitems)=index1
     iindex(2,indsign,iitems)=index2

!    DEBUG
!    write(std_out,*) ' it,v1,i1=',iitems, value1,index1
!    write(std_out,*) ' it,v2,i2=',iitems, value2,index2
!    ENDDEBUG

   end do ! indsign

   if(present(ucvol))then
     nfftot=ngfft(1) * ngfft(2) * ngfft(3)
     integrated(iitems)=integrated(iitems)*ucvol/nfftot
   end if

!  Integrate the array
!  integrated(iitems)=zero
!  do ifft=1,nfft_
!  integrated(iitems) = integrated(iitems) + array(istart + ifft)
!  enddo
!  if(present(ucvol))integrated(iitems) = integrated(iitems)*ucvol/nfft_
!  write(std_err,*)present(ucvol)
!  if(present(ucvol))then
!  write(std_err,*)ucvol
!  endif

 end do ! iitems

 ABI_DEALLOCATE(array)

!-------------------------------------------------------------------
!Enter section for FFT parallel case
!if(mpi_enreg%paral_kgb>1) spaceComm=mpi_enreg%comm_fft; reduce=.true.
 spaceComm=mpi_enreg%comm_fft; reduce=.false.
 if(mpi_enreg%nproc_fft>1) then
   spaceComm=mpi_enreg%comm_fft; reduce=.true.
 else if(mpi_enreg%nproc_wvl>1) then
   spaceComm=mpi_enreg%comm_wvl; reduce=.true.
 end if
 nproc=xmpi_comm_size(spaceComm)
 me=xmpi_comm_rank(spaceComm)

 if (reduce) then

!  Communicate all data to all processors with only two global communications
   ABI_ALLOCATE(value_fft,(5,nitems,nproc))
   ABI_ALLOCATE(index_fft,(2,2,nitems,nproc))
   value_fft(:,:,:)=zero
   index_fft(:,:,:,:)=zero
   value_fft(1,:,me + 1)=value(1,1,:)
   value_fft(2,:,me + 1)=value(2,1,:)
   value_fft(3,:,me + 1)=value(1,2,:)
   value_fft(4,:,me + 1)=value(2,2,:)
   if(present(ucvol))value_fft(5,:,me + 1)=integrated(:)
   index_fft(:,:,:,me + 1)=iindex(:,:,:)
   call xmpi_sum(value_fft,spaceComm,ierr)
   call xmpi_sum(index_fft,spaceComm,ierr)

!  Determine the global optimum and second optimum for each item
!  Also, the integrated quantities, if relevant.
   do iitems=1,nitems

     if(present(ucvol))integrated(iitems)=sum(value_fft(5,iitems,1:nproc))

     do indsign=1,2 ! Find alternatively the maximum and the minimum
       iisign=3-2*indsign

!      Initialisation
       value1=value_fft(2*indsign-1,iitems,1)
       value2=value_fft(2*indsign  ,iitems,1)
       index1=index_fft(1,indsign,iitems,1)
       index2=index_fft(2,indsign,iitems,1)

!      Loop
       do iproc=1, nproc, 1
         do ii=1,2
           if(iproc>1 .or. ii==2)then

             temp=value_fft(ii+2*(indsign-1),iitems,iproc)
!            Compares it to the second value
             if( iisign*(temp+tol12) > iisign*value2 ) then
!              Compare it to the first value
               if( iisign*(temp+tol12) > iisign*value1 ) then
                 value2=value1 ; index2=index1
                 value1=temp   ; index1=index_fft(ii,indsign,iitems,iproc)
               else
                 value2=temp   ; index2=index_fft(ii,indsign,iitems,iproc)
               end if
             end if

           end if ! if(iproc>1 .or. ii==2)
         end do ! ii
       end do ! iproc

       value(1,indsign,iitems)=value1
       value(2,indsign,iitems)=value2
       iindex(1,indsign,iitems)=index1
       iindex(2,indsign,iitems)=index2

     end do ! iisign
   end do ! iitems

   ABI_DEALLOCATE(value_fft)
   ABI_DEALLOCATE(index_fft)

 end if !if(reduce)

!-------------------------------------------------------------------

!Determines the reduced coordinates of the min and max for each item
 ABI_ALLOCATE(coord,(3,2,2,nitems))
 do iitems=1,nitems
   do indsign=1,2
     do ii=1,2
       index1=iindex(ii,indsign,iitems)
       i3=(index1-1)/n1/n2
       i2=(index1-1-i3*n1*n2)/n1
       i1=index1-1-i3*n1*n2-i2*n1
       coord(1,ii,indsign,iitems)=dble(i1)/dble(n1)+tol12
       coord(2,ii,indsign,iitems)=dble(i2)/dble(n2)+tol12
       coord(3,ii,indsign,iitems)=dble(i3)/dble(n3)+tol12
!      DEBUG
!      write(std_out,*)' ii,indsign,iitems,coord(1:3)=',ii,indsign,iitems,coord(:,ii,indsign,iitems)
!      write(std_out,*)' value ', value(ii, indsign, iitems)
!      ENDDEBUG
     end do
   end do
 end do

!-------------------------------------------------------------------------
!Output
 if (mpi_enreg%paral_kgb==0.or.mpi_enreg%me_fft==0) then
   if(.true.)then
     do iitems=1,nitems

       if(ioptrhor==4 .and. iitems>2)exit

       select case (ioptrhor)
       case(0)

         if(iitems==1) write(message,'(a)')' Total charge density [el/Bohr^3]'
         if(nspden==2)then
           if(iitems==2) write(message,'(a)')' Spin up density      [el/Bohr^3]'
           if(iitems==3) write(message,'(a)')' Spin down density    [el/Bohr^3]'
           if(iitems==4) write(message,'(a)')' Magnetization (spin up - spin down) [el/Bohr^3]'
           if(iitems==5) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         else if(nspden==4)then
           if(iitems==2) write(message,'(a)')' x component of magnetization [el/Bohr^3]'
           if(iitems==3) write(message,'(a)')' y component of magnetization [el/Bohr^3]'
           if(iitems==4) write(message,'(a)')' z component of magnetization [el/Bohr^3]'
           if(iitems==5) write(message,'(a)')' Magnetization (absolute value) [el/Bohr^3]'
           if(iitems==6) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         end if

       case(1)

         if(iitems==1) write(message,'(a)')' Total kinetic energy density [Ha/Bohr^3]'
         if(nspden==2)then
           if(iitems==2) write(message,'(a)')' Spin up density      [Ha/Bohr^3]'
           if(iitems==3) write(message,'(a)')' Spin down density    [Ha/Bohr^3]'
           if(iitems==4) write(message,'(a)')' Magnetization (spin up - spin down) [Ha/Bohr^3]'
           if(iitems==5) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         else if(nspden==4)then
           if(iitems==2) write(message,'(a)')' x component of magnetization [Ha/Bohr^3]'
           if(iitems==3) write(message,'(a)')' y component of magnetization [Ha/Bohr^3]'
           if(iitems==4) write(message,'(a)')' z component of magnetization [Ha/Bohr^3]'
           if(iitems==5) write(message,'(a)')' Magnetization (absolute value) [Ha/Bohr^3]'
           if(iitems==6) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         end if

       case(2)

         if(iitems==1) write(message,'(a)')' Gradient of the electronic density [el/Bohr^4]'
         if(nspden==2)then
           if(iitems==2) write(message,'(a)')' Spin up density      [el/Bohr^4]'
           if(iitems==3) write(message,'(a)')' Spin down density    [el/Bohr^4]'
           if(iitems==4) write(message,'(a)')' Magnetization (spin up - spin down) [el/Bohr^4]'
           if(iitems==5) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         else if(nspden==4)then
           if(iitems==2) write(message,'(a)')' x component of magnetization [el/Bohr^4]'
           if(iitems==3) write(message,'(a)')' y component of magnetization [el/Bohr^4]'
           if(iitems==4) write(message,'(a)')' z component of magnetization [el/Bohr^4]'
           if(iitems==5) write(message,'(a)')' Magnetization (absolute value) [el/Bohr^4]'
           if(iitems==6) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         end if

       case(3)

         if(iitems==1) write(message,'(a)')' Laplacian of the electronic density [el/Bohr^5]'
         if(nspden==2)then
           if(iitems==2) write(message,'(a)')' Spin up density      [el/Bohr^5]'
           if(iitems==3) write(message,'(a)')' Spin down density    [el/Bohr^5]'
           if(iitems==4) write(message,'(a)')' Magnetization (spin up - spin down) [el/Bohr^5]'
           if(iitems==5) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         else if(nspden==4)then
           if(iitems==2) write(message,'(a)')' x component of magnetization [el/Bohr^5]'
           if(iitems==3) write(message,'(a)')' y component of magnetization [el/Bohr^5]'
           if(iitems==4) write(message,'(a)')' z component of magnetization [el/Bohr^5]'
           if(iitems==5) write(message,'(a)')' Magnetization (absolute value) [el/Bohr^5]'
           if(iitems==6) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         end if

       case(4)

         if(iitems==1) write(message,'(a)')' Electron Localization Function (ELF) [min:0;max:1]'
         if(nspden==2)then
           if(iitems==2) write(message,'(a)')' Spin up ELF      [min:0;max:1]'
!            if(iitems==3) write(message,'(a)')' Spin down ELF    [min:0;max:1]'
!            if(iitems==4) write(message,'(a)')' Magnetization (spin up - spin down) [el/Bohr^4]'
!            if(iitems==5) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         else if(nspden==4)then
!            if(iitems==2) write(message,'(a)')' x component of magnetization [el/Bohr^4]'
!            if(iitems==3) write(message,'(a)')' y component of magnetization [el/Bohr^4]'
!            if(iitems==4) write(message,'(a)')' z component of magnetization [el/Bohr^4]'
!            if(iitems==5) write(message,'(a)')' Magnetization (spin up - spin down) [el/Bohr^4]'
!            if(iitems==6) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         end if


       end select

       call wrtout(iout,message,'COLL')

       write(message,'(a,es13.4,a,3f10.4)') '      Maximum= ',&
&       value(1,1,iitems),'  at reduced coord.',coord(:,1,1,iitems)
       call wrtout(iout,message,'COLL')
       if(option==2)then
         write(message,'(a,es13.4,a,3f10.4)')' Next maximum= ',&
&         value(2,1,iitems),'  at reduced coord.',coord(:,2,1,iitems)
         call wrtout(iout,message,'COLL')
       end if
       write(message,'(a,es13.4,a,3f10.4)') '      Minimum= ',&
&       value(1,2,iitems),'  at reduced coord.',coord(:,1,2,iitems)
       call wrtout(iout,message,'COLL')
       if(option==2)then
         write(message,'(a,es13.4,a,3f10.4)')' Next minimum= ',&
&         value(2,2,iitems),'  at reduced coord.',coord(:,2,2,iitems)
         call wrtout(iout,message,'COLL')
       end if
       if(present(ucvol))then
         if(.not.(nspden==2.and.iitems==5) .and. .not.(nspden==4.and.iitems==6))then
           if(abs(integrated(iitems))<tol10)integrated(iitems)=zero
           write(message,'(a,es13.4)')'   Integrated= ',integrated(iitems)
           call wrtout(iout,message,'COLL')
         end if
       end if

     end do ! iitems
   end if

   if(.false.)then

     select case(optrhor)
     case(0)
       write(txt1_in_mssg, '(a)')" Min el dens="
       write(txt2_in_mssg, '(a)')" el/bohr^3 at reduced coord."
       write(txt3_in_mssg, '(a)')" Max el dens="
     case(1)
       write(txt1_in_mssg, '(a)')" Min kin energy dens="
       write(txt2_in_mssg, '(a)')" bohr^(-5) at reduced coord."
       write(txt3_in_mssg, '(a)')" Max kin energy dens="
     end select

     write(message, '(a,a,1p,e12.4,a,0p,3f8.4)' ) ch10,&
&     trim(txt1_in_mssg),value(1,2,1),&
&     trim(txt2_in_mssg),coord(:,1,2,1)
     call wrtout(iout,message,'COLL')
     if(option==2)then
       write(message, '(a,1p,e12.4,a,0p,3f8.4)' ) &
&       ',   next min=',value(2,2,1),&
&       trim(txt2_in_mssg),coord(:,2,2,1)
       call wrtout(iout,message,'COLL')
     end if
     write(message, '(a,1p,e12.4,a,0p,3f8.4)' )&
&     trim(txt3_in_mssg),value(1,1,1),&
&     trim(txt2_in_mssg),coord(:,1,1,1)
     call wrtout(iout,message,'COLL')
     if(option==2)then
       write(message, '(a,1p,e12.4,a,0p,3f8.4)' )&
&       ',   next max=',value(2,1,1),&
&       trim(txt2_in_mssg),coord(:,2,1,1)
       call wrtout(iout,message,'COLL')
     end if

     if(nspden>=2)then
       write(message, '(a,a,1p,e12.4,a,0p,3f8.4)' ) ch10,&
&       ',Min spin pol zeta=',value(1,2,4+nspden/2),&
&       ' at reduced coord.',coord(:,1,2,4+nspden/2)
       call wrtout(iout,message,'COLL')
       if(option==2)then
         write(message, '(a,1p,e12.4,a,0p,3f8.4)' )&
&         ',         next min=',value(2,2,4+nspden/2),&
&         ' at reduced coord.',coord(:,2,2,4+nspden/2)
         call wrtout(iout,message,'COLL')
       end if
       write(message, '(a,1p,e12.4,a,0p,3f8.4)' )&
&       ',Max spin pol zeta=',value(1,1,4+nspden/2),&
&       ' at reduced coord.',coord(:,1,1,4+nspden/2)
       call wrtout(iout,message,'COLL')
       if(option==2)then
         write(message, '(a,1p,e12.4,a,0p,3f8.4)' )&
&         ',         next max=',value(2,1,4+nspden/2),&
&         ' at reduced coord.',coord(:,2,1,4+nspden/2)
         call wrtout(iout,message,'COLL')
       end if
     end if ! nspden

   end if ! second section always true

   if(nspden==2 .and. .false.)then
     write(message,'(a)')&
&     '                               Position in reduced coord.       (  x         y         z )'
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Total  el-den) : [el/Bohr^3]',&
&     rhomn1(1),'  at',ri_rhomn1(1,1),ri_rhomn1(2,1),ri_rhomn1(3,1)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Spin-up   den) : [el/Bohr^3]',&
&     rhomn1(2),'  at',ri_rhomn1(1,2),ri_rhomn1(2,2),ri_rhomn1(3,2)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Spin-down den) : [el/Bohr^3]',&
&     zetmn1(1),'  at',ri_zetmn1(1,1),ri_zetmn1(2,1),ri_zetmn1(3,1)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Spin pol zeta) :   [m/|m|]  ',&
&     zetmn1(2),'  at',ri_zetmn1(1,2),ri_zetmn1(2,2),ri_zetmn1(3,2)
     call wrtout(iout,message,'COLL')
     if(option==2)then
       write(message,'(a,es13.4,a,3f10.4)')' Next minimum (Total  el-den) : [el/Bohr^3]',&
&       rhomn2(1),'  at',ri_rhomn2(1,1),ri_rhomn2(2,1),ri_rhomn2(3,1)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next minimum (Spin-up   den) : [el/Bohr^3]',&
&       rhomn2(2),'  at',ri_rhomn2(1,2),ri_rhomn2(2,2),ri_rhomn2(3,2)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next minimum (Spin-down den) : [el/Bohr^3]',&
&       zetmn2(1),'  at',ri_zetmn2(1,1),ri_zetmn2(2,1),ri_zetmn2(3,1)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next minimum (Spin pol zeta) :   [m/|m|]  ',&
&       zetmn2(2),'  at',ri_zetmn2(1,2),ri_zetmn2(2,2),ri_zetmn2(3,2)
       call wrtout(iout,message,'COLL')
     end if
     write(message,*)' '
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Total  el-den) : [el/Bohr^3]',&
&     rhomx1(1),'  at',ri_rhomx1(1,1),ri_rhomx1(2,1),ri_rhomx1(3,1)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Spin-up   den) : [el/Bohr^3]',&
&     rhomx1(2),'  at',ri_rhomx1(1,2),ri_rhomx1(2,2),ri_rhomx1(3,2)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Spin-down den) : [el/Bohr^3]',&
&     zetmx1(1),'  at',ri_zetmx1(1,1),ri_zetmx1(2,1),ri_zetmx1(3,1)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Spin pol zeta) :   [m/|m|]  ',&
&     zetmx1(2),'  at',ri_zetmx1(1,2),ri_zetmx1(2,2),ri_zetmx1(3,2)
     call wrtout(iout,message,'COLL')
     if(option==2)then
       write(message,'(a,es13.4,a,3f10.4)')' Next maximum (Total  el-den) : [el/Bohr^3]',&
&       rhomx2(1),'  at',ri_rhomx2(1,1),ri_rhomx2(2,1),ri_rhomx2(3,1)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next maximum (Spin-up   den) : [el/Bohr^3]',&
&       rhomx2(2),'  at',ri_rhomx2(1,2),ri_rhomx2(2,2),ri_rhomx2(3,2)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next maximum (Spin-down den) : [el/Bohr^3]',&
&       zetmx2(1),'  at',ri_zetmx2(1,1),ri_zetmx2(2,1),ri_zetmx2(3,1)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next maximum (Spin pol zeta) :   [m/|m|]  ',&
&       zetmx2(2),'  at',ri_zetmx2(1,2),ri_zetmx2(2,2),ri_zetmx2(3,2)
       call wrtout(iout,message,'COLL')
     end if
   end if

   if(nspden==4 .and. .false.)then
     write(message,'(a)')&
&     '                               Position in reduced coord.       (  x         y         z )'
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Total  el-den) : [el/Bohr^3]',&
&     rhomn1(1),'  at',ri_rhomn1(1,1),ri_rhomn1(2,1),ri_rhomn1(3,1)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Magnetizat.-x) :   [m/|m|]  ',&
&     rhomn1(2),'  at',ri_rhomn1(1,2),ri_rhomn1(2,2),ri_rhomn1(3,2)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Magnetizat.-y) :   [m/|m|]  ',&
&     rhomn1(3),'  at',ri_rhomn1(1,3),ri_rhomn1(2,3),ri_rhomn1(3,3)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Magnetizat.-z) :   [m/|m|]  ',&
&     rhomn1(4),'  at',ri_rhomn1(1,4),ri_rhomn1(2,4),ri_rhomn1(3,4)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Spin pol zeta) :   [m/|m|]  ',&
&     zetmn1(1),'  at',ri_zetmn1(1,1),ri_zetmn1(2,1),ri_zetmn1(3,1)
     call wrtout(iout,message,'COLL')
     if(option==2)then
       write(message,'(a,es13.4,a,3f10.4)')' Next-Minimum (Total  el-den) : [el/Bohr^3]',&
&       rhomn2(1),'  at',ri_rhomn2(1,1),ri_rhomn2(2,1),ri_rhomn2(3,1)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next-Minimum (Magnetizat.-x) :   [m/|m|]  ',&
&       rhomn2(2),'  at',ri_rhomn2(1,2),ri_rhomn2(2,2),ri_rhomn2(3,2)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next-Minimum (Magnetizat.-y) :   [m/|m|]  ',&
&       rhomn2(3),'  at',ri_rhomn2(1,3),ri_rhomn2(2,3),ri_rhomn2(3,3)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next-Minimum (Magnetizat.-z) :   [m/|m|]  ',&
&       rhomn2(4),'  at',ri_rhomn2(1,4),ri_rhomn2(2,4),ri_rhomn2(3,4)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next-Minimum (Spin pol zeta) :   [m/|m|]  ',&
&       zetmn2(1),'  at',ri_zetmn2(1,1),ri_zetmn2(2,1),ri_zetmn2(3,1)
       call wrtout(iout,message,'COLL')
     end if
     write(message,*)' '
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Total  el-den) : [el/Bohr^3]',&
&     rhomx1(1),'  at',ri_rhomx1(1,1),ri_rhomx1(2,1),ri_rhomx1(3,1)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Magnetizat.-x) :   [m/|m|]  ',&
&     rhomx1(2),'  at',ri_rhomx1(1,2),ri_rhomx1(2,2),ri_rhomx1(3,2)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Magnetizat.-y) :   [m/|m|]  ',&
&     rhomx1(3),'  at',ri_rhomx1(1,3),ri_rhomx1(2,3),ri_rhomx1(3,3)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Magnetizat.-z) :   [m/|m|]  ',&
&     rhomx1(4),'  at',ri_rhomx1(1,4),ri_rhomx1(2,4),ri_rhomx1(3,4)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Spin pol zeta) :   [m/|m|]  ',&
&     zetmx1(1),'  at',ri_zetmx1(1,1),ri_zetmx1(2,1),ri_zetmx1(3,1)
     call wrtout(iout,message,'COLL')
     if(option==2)then
       write(message,'(a,es13.4,a,3f10.4)')' Next-Maximum (Total  el-den) : [el/Bohr^3]',&
&       rhomx2(1),'  at',ri_rhomx2(1,1),ri_rhomx2(2,1),ri_rhomx2(3,1)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next-Maximum (Magnetizat.-x) :   [m/|m|]  ',&
&       rhomx2(2),'  at',ri_rhomx2(1,2),ri_rhomx2(2,2),ri_rhomx2(3,2)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next-Maximum (Magnetizat.-y) :   [m/|m|]  ',&
&       rhomx2(3),'  at',ri_rhomx2(1,3),ri_rhomx2(2,3),ri_rhomx2(3,3)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next-Maximum (Magnetizat.-z) :   [m/|m|]  ',&
&       rhomx2(4),'  at',ri_rhomx2(1,4),ri_rhomx2(2,4),ri_rhomx2(3,4)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next-Maximum (Spin pol zeta) :   [m/|m|]  ',&
&       zetmx2(1),'  at',ri_zetmx2(1,1),ri_zetmx2(2,1),ri_zetmx2(3,1)
       call wrtout(iout,message,'COLL')
     end if
   end if
 end if

 ABI_DEALLOCATE(coord)
 ABI_DEALLOCATE(value)
 ABI_DEALLOCATE(iindex)
 ABI_DEALLOCATE(integrated)

end subroutine prtrhomxmn
!!***

!!****f* m_mkrho/read_atomden
!! NAME
!! read_atomden
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (SM,VR,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! natom : number of atoms in cell
!! nfft=(effective) number of FFT grid points (for this processor) - fine grid
!! ngfft(18)=contain all needed information about 3D FFT,
!! nspden : number of spin densities
!! ntypat : number of types of atoms in the cell
!! typat(natom) : list of atom types
!!
!! OUTPUT
!! rhor_atm(nfft,nspden) : full electron density on the (fine) grid
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      atomden
!!
!! SOURCE

subroutine read_atomden(MPI_enreg,natom,nfft,ngfft,nspden,ntypat, &
&                       rhor_atm,typat,rprimd,xred,prtvol,file_prefix)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,nspden,ntypat,prtvol
!arrays
 type(MPI_type),intent(in) :: MPI_enreg
 integer,intent(in) :: ngfft(18),typat(natom)
 real(dp), intent(in) :: rprimd(3,3),xred(3,natom)
 real(dp),intent(inout) :: rhor_atm(nfft,nspden)
 character(len=7), intent(in) :: file_prefix

!Local variables-------------------------------
!scalars
 character(len=500) :: message
 character(len=120) :: filename
 character(len=7) :: calctype='replace'
 integer :: igrid,i,i1,i2,i3,io_err,itypat,unt
 integer :: natomgrmax,nlines,ngrid,n1,n2,n3
 real(dp) :: difx,dify,difz,ucvol!,norm
!arrays
 integer :: natomgr(ntypat)
 real(dp) :: a(3),b(3),c(3)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: atomrgrid(:,:),r_vec_grid(:,:),density(:,:)
 real(dp),allocatable :: rho(:)

! ************************************************************************

!Initialise various variables
 ngrid = nfft
 a(:) = rprimd(:,1)
 b(:) = rprimd(:,2)
 c(:) = rprimd(:,3)
 ABI_ALLOCATE(rho,(ngrid))
 if (nspden/=1) then
   MSG_ERROR('read_atomden: Only nspden=1 allowed.')
 end if
 rho = rhor_atm(:,1)
 gmet=zero;gprimd=zero;rmet=zero;ucvol=zero


!Calculate the r vector (reduced coord.) of the fine gridpoints
 ABI_ALLOCATE(r_vec_grid,(3,ngrid))
 igrid = 0
 n1 = ngfft(1)
 n2 = ngfft(2)
 n3 = ngfft(3)
 do i3=0,n3-1
   difz=dble(i3)/dble(n3)
   do i2=0,n2-1
     dify=dble(i2)/dble(n2)
     do i1=0,n1-1
       difx=dble(i1)/dble(n1)
       igrid = igrid + 1
       r_vec_grid(1,igrid)=difx*rprimd(1,1)+dify*rprimd(1,2)+difz*rprimd(1,3)
       r_vec_grid(2,igrid)=difx*rprimd(2,1)+dify*rprimd(2,2)+difz*rprimd(2,3)
       r_vec_grid(3,igrid)=difx*rprimd(3,1)+dify*rprimd(3,2)+difz*rprimd(3,3)
     end do
   end do
 end do
 if (igrid/=ngrid) then
   MSG_ERROR('read_atomden: igrid not equal to ngrid')
 end if

!Read in atomic density data for each atom type
!first check how many datapoints are in each file
 do itypat=1,ntypat
   filename='';io_err=0;
   if (itypat>0)  write(filename,'(a,a,i1,a)') trim(file_prefix), &
&   '_density_atom_type',itypat,'.dat'
   if (itypat>10) write(filename,'(a,a,i2,a)') trim(file_prefix), &
&   '_density_atom_type',itypat,'.dat'
   if (open_file(filename, message, newunit=unt, status='old',action='read') /= 0) then
     write(std_out,*) 'ERROR in read_atomden: Could not open file: ',filename
     write(std_out,*) ' Current implementation requires this file to be present'
     write(std_out,*) ' for each type of atom.'
     write(std_out,*)trim(message)
     MSG_ERROR("Cannot continue")
   end if
!  Check number of lines in file
   nlines = 1;io_err=0;
   do
     read(unt,*,iostat=io_err)
     if (io_err<0) exit
     nlines = nlines + 1
   end do
   close(unt)
   natomgr(itypat) = nlines - 2
 end do ! Atom type
!Allocate arrays and read in data
 natomgrmax = maxval(natomgr)
 ABI_ALLOCATE(atomrgrid,(natomgrmax,ntypat))
 ABI_ALLOCATE(density,(natomgrmax,ntypat))
 atomrgrid = zero ; density = zero
 do itypat=1,ntypat
   filename='';io_err=0;
   if (itypat>0)  write(filename,'(a,a,i1,a)') trim(file_prefix), &
&   '_density_atom_type',itypat,'.dat'
   if (itypat>10) write(filename,'(a,a,i2,a)') trim(file_prefix), &
&   '_density_atom_type',itypat,'.dat'
   if (open_file(filename,message,newunit=unt,status='old',action='read') /= 0) then
     MSG_ERROR(message)
   end if
   read(unt,*) ! Skip comment line
   do i=1,natomgr(itypat)
     read(unt,*) atomrgrid(i,itypat),density(i,itypat)
   end do
   close(unt)
   if (atomrgrid(1,itypat)/=zero) then
     write(std_out,*) 'ERROR in read_atomden, in file: ',filename
     write(std_out,*) ' First gridpoint has to be the origin.'
     MSG_ERROR("Cannot continue")
   end if
 end do ! Atom type

!write(std_out,*) '*** --- In read_atomden before call--- ***'
!write(std_out,*) '  calctype:',calctype,' natom:',natom
!write(std_out,*) '    ntypat:',ntypat,' typat:',typat
!write(std_out,*) '     ngrid:',ngrid
!write(std_out,*) '         a:',a
!write(std_out,*) '         b:',b
!write(std_out,*) '         c:',c
!write(std_out,*) '      xred:',xred
!write(std_out,*) '   natomgr:',natomgr
!write(std_out,*) 'natomgrmax:',natomgrmax
!write(std_out,*) ' atomrgrid:',atomrgrid
!write(std_out,*) '   density:',density
!write(std_out,*) 'r_vec_grid:'
!write(std_out,*) r_vec_grid

!Call atomden
 call atomden(MPI_enreg,natom,ntypat,typat,ngrid,r_vec_grid,rho,a,b,c,xred, &
& natomgr,natomgrmax,atomrgrid,density,prtvol,calctype)

!if (prtvol>9) then ! calculate norm
!call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
!norm = SUM(rho(:))*ucvol/dble(n1*n2*n3)
!write(message,'(a,F8.4)') '  In read_atomden - NORM OF DENSITY: ',norm
!call wrtout(std_out,message,'COLL')
!end if

 rhor_atm(:,1) = rho

 if (allocated(atomrgrid))  then
   ABI_DEALLOCATE(atomrgrid)
 end if
 if (allocated(density))  then
   ABI_DEALLOCATE(density)
 end if
 if (allocated(r_vec_grid))  then
   ABI_DEALLOCATE(r_vec_grid)
 end if
 if (allocated(rho))  then
   ABI_DEALLOCATE(rho)
 end if

 return

end subroutine read_atomden
!!***

!!****f* m_mkrho/atomden
!! NAME
!! atomden
!!
!! FUNCTION
!! Construct atomic proto-bulk density (i.e. the superposed density
!! from neutral, isolated atoms at the bulk atomic positions).
!! This is useful if one wants to construct the bonding density:
!!
!! rho^{bnd} = rho^{bulk}(r)
!!                 - \sum_{\alpha}\rho^{atm}_{\alpha}(r-R_{\alpha})
!!
!! Where rho^{bulk} is the bulk density, rho^{atm} the atomic density
!! and the index \alpha sums over all atoms. the R_{\alpha} are the
!! atomic positions in the bulk. This routine calculates the sum over
!! rho^{atm}_{\alpha}(r-R_{\alpha}) on a grid.
!!
!! Units are atomic.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (SM,VR,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! calctype : type of calculation
!!          'replace' zero the input/output density array
!!          'add'     add to the input/output density array
!! natom : number of atoms in cell
!! ntypat : number of different types of atoms in cell
!! typat(natom) : type of each atom
!! ngrid : number of gridpoints
!! r_vec_grid(3,ngrid) : real (non-reduced) coordinates for grid points
!! rho(ngrid) : input/output density array
!! a(3),b(3),c(3) : real-space basis vectors
!! atom_pos(3,natom) : reduced coordinates for atomic positions
!! natomgr(ntypat) : number of gridpoints for each atomic density grid
!! natomgrmax : max(natomgr(ntypat))
!! atomrgrid(natomgrmax,ntypat)
!! density(natomgrmax,ntypat)
!!
!! OUTPUT
!! rho(ngrid) : input/output density array
!!
!! SIDE EFFECTS
!!
!! NOTES
!! There are two ways to compile the proto density in real space
!! for a solid. One alternative is that the density is calculated
!! for an extended grid encompassing the sphere of points around
!! one atom, and the results are folded back into the unit cell.
!! On the other hand one can, around each grid point, identify the
!! number of atoms in a sphere equivalent to the length of the radial
!! grid for each type of atom.
!! The second approach, with some modification, is taken here. The
!! numer of atoms in a supercell cell are listed such that the supercell
!! encompasses the atoms which could contribute to any point in the grid.
!! That list is kept and cycled through, to avoid recalculating it at
!! each point.
!! Note that the density calculated from the atom is the spherical
!! average, since there is no preferred direction without any
!! external field (and it's simpler)
!!
!!
!! PARENTS
!!      read_atomden
!!
!! CHILDREN
!!      sort_dp,spline,splint,wrtout,xmpi_barrier,xmpi_sum_master
!!
!! SOURCE

subroutine atomden(MPI_enreg,natom,ntypat,typat,ngrid,r_vec_grid,rho,a,b,c,atom_pos, &
&                  natomgr,natomgrmax,atomrgrid,density,prtvol,calctype)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat,ngrid,natomgrmax,prtvol
 character(len=7),intent(in) :: calctype
!arrays
 type(MPI_type),intent(in) :: MPI_enreg
 integer,intent(in) :: typat(natom),natomgr(ntypat)
 real(dp),intent(in) :: r_vec_grid(3,ngrid),a(3),b(3),c(3)
 real(dp),intent(in) :: atom_pos(3,natom),atomrgrid(natomgrmax,ntypat)
 real(dp),intent(in) :: density(natomgrmax,ntypat)
 real(dp),intent(inout) :: rho(ngrid)

!Local variables-------------------------------
!scalars
 character(len=500) :: message
 integer :: cnt,delta,i,l,m,n,iatom,itypat,igrid,ncells,n_grid_p
 integer :: ierr,spaceComm,nprocs,master,rank,remainder
 real(dp) :: a_norm,b_norm,c_norm
 real(dp) :: r_max,R_sphere_max,dp_dummy,ybcbeg,ybcend
!arrays
 integer :: n_equiv_atoms(ntypat),grid_index(ngrid)
 integer :: my_start_equiv_atoms(ntypat)
 integer :: my_end_equiv_atoms(ntypat)
 integer :: l_min(ntypat),m_min(ntypat),n_min(ntypat)
 integer :: l_max(ntypat),m_max(ntypat),n_max(ntypat)
 real(dp) :: center(3),dp_vec_dummy(3),delta_a(3),delta_b(3),delta_c(3)
 real(dp) :: r_atom(3),grid_distances(ngrid)
 integer, allocatable :: new_index(:),i_1d_dummy(:)
 real(dp),allocatable :: equiv_atom_dist(:,:),equiv_atom_pos(:,:,:),rho_temp(:,:)
 real(dp),allocatable :: dp_1d_dummy(:),dp_2d_dummy(:,:),ypp(:)
 real(dp),allocatable :: x_fit(:),y_fit(:)


! ************************************************************************

!initialise and check parallel execution
 spaceComm=MPI_enreg%comm_cell
 nprocs=xmpi_comm_size(spaceComm)
 rank=MPI_enreg%me_kpt

 master=0

!initialise variables and vectors
 a_norm = sqrt(dot_product(a,a))
 b_norm = sqrt(dot_product(b,b))
 c_norm = sqrt(dot_product(c,c))
 center = (a+b+c)*half
 dp_dummy = dot_product(a,b)/(b_norm*b_norm)
 dp_vec_dummy = dp_dummy*b
 delta_a = a - dp_vec_dummy
 dp_dummy = dot_product(b,a)/(a_norm*a_norm)
 dp_vec_dummy = dp_dummy*a
 delta_b = b - dp_vec_dummy
 dp_dummy = dot_product(c,(a+b))/(dot_product((a+b),(a+b)))
 dp_vec_dummy = dp_dummy*(a+b)
 delta_c = c - dp_vec_dummy
 ABI_ALLOCATE(rho_temp,(ngrid,ntypat))
 rho_temp = zero

!write(std_out,*) '*** --- In atomden --- ***'
!write(std_out,*) ' a_norm:',a_norm,' b_norm:',b_norm,' c_norm:',c_norm
!write(std_out,*) 'delta_a:',delta_a,'delta_b:',delta_b,'delta_c:',delta_c
!write(std_out,*) ' center:',center

!Find supercell which will contain all possible contributions
!for all atoms, and enumerate positions for all atoms
!TODO list of atoms can be "pruned", i.e identify all atoms
!that can't possibly contribute and remove from list.
!Should be most important for very oblique cells
 do itypat=1,ntypat
   R_sphere_max = atomrgrid(natomgr(itypat),itypat)
   l_min(itypat) = -ceiling(R_sphere_max/sqrt(dot_product(delta_a,delta_a)))
   l_max(itypat) = -l_min(itypat)
   m_min(itypat) = -ceiling(R_sphere_max/sqrt(dot_product(delta_b,delta_b)))
   m_max(itypat) = -m_min(itypat)
   n_min(itypat) = -ceiling(R_sphere_max/sqrt(dot_product(delta_c,delta_c)))
   n_max(itypat) = -n_min(itypat)
   ncells = (l_max(itypat)-l_min(itypat)+1) &
&   *(m_max(itypat)-m_min(itypat)+1) &
&   *(n_max(itypat)-n_min(itypat)+1)
   n_equiv_atoms(itypat) = 0
   do iatom=1,natom
     if (typat(iatom)==itypat) then
       n_equiv_atoms(itypat) = n_equiv_atoms(itypat) + ncells
     end if ! if type=itypat
   end do ! number of atoms per cell
   if ((rank==master).and.(prtvol>9)) then
     write(message,'(a)') '*** --- In atomden --- find box ***'
     call wrtout(std_out,message,'COLL')
     write(message,'(a,I4)') ' itypat:',itypat
     call wrtout(std_out,message,'COLL')
     write(message,'(2(a,I4))') ' l_min:',l_min(itypat),' l_max:',l_max(itypat)
     call wrtout(std_out,message,'COLL')
     write(message,'(2(a,I4))') ' m_min:',m_min(itypat),' m_max:',m_max(itypat)
     call wrtout(std_out,message,'COLL')
     write(message,'(2(a,I4))') ' n_min:',n_min(itypat),' n_max:',n_max(itypat)
     call wrtout(std_out,message,'COLL')
     write(message,'(2(a,I4))') ' n_equiv_atoms:',n_equiv_atoms(itypat)
     call wrtout(std_out,message,'COLL')
   end if
 end do !atom type

!allocate arrays
 n = maxval(n_equiv_atoms)
 ABI_ALLOCATE(equiv_atom_pos,(3,n,ntypat))
 ABI_ALLOCATE(equiv_atom_dist,(n,ntypat))
 equiv_atom_pos = zero
 equiv_atom_dist = zero

!Find positions and distance of atoms from center of cell
 do itypat=1,ntypat
   i = 1
   do l=l_min(itypat),l_max(itypat)
     do m=m_min(itypat),m_max(itypat)
       do n=n_min(itypat),n_max(itypat)
         do iatom=1,natom
           if (typat(iatom)==itypat) then
             if (i>n_equiv_atoms(itypat)) then
               MSG_ERROR('atomden: i>n_equiv_atoms')
             end if
             equiv_atom_pos(:,i,itypat) = (atom_pos(1,iatom)+dble(l))*a &
&             + (atom_pos(2,iatom)+dble(m))*b &
&             + (atom_pos(3,iatom)+dble(n))*c
             dp_vec_dummy = equiv_atom_pos(:,i,itypat)-center
             equiv_atom_dist(i,itypat) = &
&             sqrt(dot_product(dp_vec_dummy,dp_vec_dummy))
             i = i + 1
           end if
         end do
       end do !n
     end do !m
   end do !l
!  write(std_out,*) '*** --- In atomden --- find equiv ***'
!  write(std_out,*) ' itypat:',itypat
!  write(std_out,*) ' equiv_atom_pos:'
!  write(std_out,*) equiv_atom_pos(:,:,itypat)
!  write(std_out,*) ' equiv_atom_dist:',equiv_atom_dist(:,itypat)
 end do !atom type

!Sort the atoms after distance so that the density from the ones
!furthest away can be added first. This is to prevent truncation error.
 do itypat=1,ntypat
   n = n_equiv_atoms(itypat)
   ABI_ALLOCATE(dp_1d_dummy,(n))
   ABI_ALLOCATE(new_index,(n))
   ABI_ALLOCATE(dp_2d_dummy,(3,n))
   dp_1d_dummy = equiv_atom_dist(1:n,itypat)
   dp_2d_dummy = equiv_atom_pos(1:3,1:n,itypat)
   do i=1,n
     new_index(i) = i
   end do
   call sort_dp(n,dp_1d_dummy,new_index,tol14)
   do i=1,n
!    write(std_out,*) i,' -> ',new_index(i)
     equiv_atom_pos(1:3,n+1-i,itypat) = dp_2d_dummy(1:3,new_index(i))
     equiv_atom_dist(1:n,itypat) = dp_1d_dummy
   end do
   ABI_DEALLOCATE(dp_1d_dummy)
   ABI_DEALLOCATE(new_index)
   ABI_DEALLOCATE(dp_2d_dummy)
!  write(std_out,*) '*** --- In atomden ---  sorting atoms ***'
!  write(std_out,*) ' itypat:',itypat
!  write(std_out,*) ' equiv_atom_pos:'
!  write(std_out,*) equiv_atom_pos(:,:,itypat)
!  write(std_out,*) ' equiv_atom_dist:',equiv_atom_dist(:,itypat)
 end do ! atom type

!Divide the work in case of parallel execution
 if (nprocs==1) then ! Make sure everything runs with one proc
   if (prtvol>9) then
     write(message,'(a)') '  In atomden - number of processors:     1'
     call wrtout(std_out,message,'COLL')
     write(message,'(a)') '  Calculation of proto-atomic density done in serial'
     call wrtout(std_out,message,'COLL')
   end if
   do itypat=1,ntypat
     if (prtvol>9) then
       write(message,'(a,I6)') '  Number of equivalent atoms:',n_equiv_atoms(itypat)
       call wrtout(std_out,message,'COLL')
     end if
     my_start_equiv_atoms(itypat) = 1
     my_end_equiv_atoms(itypat) = n_equiv_atoms(itypat)
   end do
 else
   if (rank==master.and.prtvol>9) then
     write(message,'(a,I5)') '  In atomden - number of processors:',nprocs
     call wrtout(std_out,message,'COLL')
     write(message,'(a)') '  Calculation of proto-atomic density done in parallel'
     call wrtout(std_out,message,'COLL')
   end if
   do itypat=1,ntypat
     if (rank==master.and.prtvol>9) then
       write(message,'(a,I6)') '  Number of equivalent atoms:',n_equiv_atoms(itypat)
       call wrtout(std_out,message,'COLL')
     end if
!    Divide the atoms among the processors by shuffling indices
     delta = int(floor(real(n_equiv_atoms(itypat))/real(nprocs)))
     remainder = n_equiv_atoms(itypat)-nprocs*delta
     my_start_equiv_atoms(itypat) = 1+rank*delta
     my_end_equiv_atoms(itypat) = (rank+1)*delta
!    Divide the remainder points among the processors
!    by shuffling indices
     if ((rank+1)>remainder) then
       my_start_equiv_atoms(itypat) = my_start_equiv_atoms(itypat) + remainder
       my_end_equiv_atoms(itypat) = my_end_equiv_atoms(itypat) + remainder
     else
       my_start_equiv_atoms(itypat) = my_start_equiv_atoms(itypat) + rank
       my_end_equiv_atoms(itypat) = my_end_equiv_atoms(itypat) + rank + 1
     end if
     if (prtvol>9) then
       write(message,'(a,I3)') '          For atom type: ',itypat
       call wrtout(std_out,message,'PERS')
!      write(message,'(a,I6)') '  I''ll take atoms from: ',my_start_equiv_atoms(itypat)
!      call wrtout(std_out,message,'PERS')
!      write(message,'(a,I6)') '           total for me: ',my_end_equiv_atoms(itypat)
!      call wrtout(std_out,message,'PERS')
       write(message,'(a,I6)') '            total for me: ', &
&       my_end_equiv_atoms(itypat)+1-my_start_equiv_atoms(itypat)
       call wrtout(std_out,message,'PERS')
     end if
   end do
 end if

!Loop over types of atoms and equivalent atoms and
!interpolate density onto grid
 do itypat=1,ntypat
!  do iatom=my_start_equiv_atoms(itypat),my_end_equiv_atoms(itypat)

   cnt = 0
   iatom = rank+1 - nprocs
!  Parallel execution of loop
   do
     cnt = cnt + 1
     iatom = iatom + nprocs
     if (iatom>n_equiv_atoms(itypat)) exit ! Exit if index is too large

     if (mod(cnt,100)==0.and.prtvol>0) then
       write(message,'(2(a,I6))') ' atoms so far',cnt,' of: ',n_equiv_atoms(itypat)/nprocs
       call wrtout(std_out,message,'PERS')
     end if

     r_max = atomrgrid(natomgr(itypat),itypat)
     r_atom = equiv_atom_pos(:,iatom,itypat)

!    Set up an array with the gridpoint distances
     i = 1
     grid_distances = zero
     grid_index = 0
     do igrid=1,ngrid
       dp_vec_dummy(:) = r_vec_grid(:,igrid) - r_atom(:)
       dp_dummy = sqrt(dot_product(dp_vec_dummy,dp_vec_dummy))
       if (dp_dummy <= r_max) then
         grid_distances(i) = dp_dummy
         grid_index(i) = igrid
         i = i + 1
       else
         cycle ! cycle if point is too far away
       end if
     end do
     n_grid_p = i - 1

     if (n_grid_p==0) cycle ! Cycle if no point needs
!    to be interpolated

!    Sort points to be interpolated in ascending order
     ABI_ALLOCATE(dp_1d_dummy,(n_grid_p))
     ABI_ALLOCATE(new_index,(n_grid_p))
     ABI_ALLOCATE(i_1d_dummy,(n_grid_p))
     dp_1d_dummy = grid_distances(1:n_grid_p)
     do i=1,n_grid_p
       new_index(i) = i
     end do
     call sort_dp(n_grid_p,dp_1d_dummy,new_index,tol16)
     grid_distances(1:n_grid_p) = dp_1d_dummy
     i_1d_dummy = grid_index(1:n_grid_p)
     do i=1,n_grid_p
!      write(std_out,*) i_1d_dummy(i),' -> ',i_1d_dummy(new_index(i))
       grid_index(i) = i_1d_dummy(new_index(i))
     end do
     ABI_DEALLOCATE(dp_1d_dummy)
     ABI_DEALLOCATE(new_index)
     ABI_DEALLOCATE(i_1d_dummy)

!    Interpolate density onto all grid points
     ABI_ALLOCATE(ypp,(natomgr(itypat)))
     ABI_ALLOCATE(x_fit,(n_grid_p))
     ABI_ALLOCATE(y_fit,(n_grid_p))
     ypp = zero; y_fit = zero
     ybcbeg = zero; ybcend = zero
     x_fit = grid_distances(1:n_grid_p)
     call spline(atomrgrid(1:natomgr(itypat),itypat), &
&     density(1:natomgr(itypat),itypat), &
&     natomgr(itypat),ybcbeg,ybcend,ypp)
     call splint(natomgr(itypat),atomrgrid(1:natomgr(itypat),itypat), &
&     density(1:natomgr(itypat),itypat),ypp,n_grid_p, &
&     x_fit,y_fit)

!    Save the interpolated points to grid
     do i=1,n_grid_p
       rho_temp(grid_index(i),itypat) = rho_temp(grid_index(i),itypat) + y_fit(i)
     end do
     ABI_DEALLOCATE(ypp)
     ABI_DEALLOCATE(x_fit)
     ABI_DEALLOCATE(y_fit)

   end do ! n equiv atoms
 end do ! type of atom

!Collect all contributions to rho_temp if
!we are running in parallel
 if (nprocs>1) then
   call xmpi_barrier(spaceComm)
   call xmpi_sum_master(rho_temp,master,spaceComm,ierr)
   call xmpi_barrier(spaceComm)
   if (prtvol>9) then
     write(message,'(a)') '  In atomden - contributions to rho_temp collected'
     call wrtout(std_out,message,'PERS')
   end if
 end if

!Now rho_temp contains the atomic protodensity for each atom.
!Check whether this is to replace or be added to the input/output array
!and sum up contributions
 if (trim(calctype)=='replace') rho = zero
 do itypat=1,ntypat
   rho(:) = rho(:) + rho_temp(:,itypat)
 end do

!deallocations
 if (allocated(rho_temp))  then
   ABI_DEALLOCATE(rho_temp)
 end if
 if (allocated(equiv_atom_pos))  then
   ABI_DEALLOCATE(equiv_atom_pos)
 end if
 if (allocated(equiv_atom_dist))  then
   ABI_DEALLOCATE(equiv_atom_dist)
 end if
!if (allocated()) deallocate()

 return

 end subroutine atomden
!!***

end module m_mkrho
!!***
