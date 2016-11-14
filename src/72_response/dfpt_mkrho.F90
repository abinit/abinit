!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_mkrho
!! NAME
!! dfpt_mkrho
!!
!! FUNCTION
!! Compute RF charge density rho1(r) and rho1(G) in electrons/bohr**3
!! from input RF and GS wavefunctions, band occupations, and k point weights.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR, LSI, AR, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=wf in G space
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=first-order wf in G space
!!  cplex=1 if rhor1 is real, 2 if rhor1 is complex
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the
!!    storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates, GS data.
!!  kg1(3,mpw1*mkmem1)=reduced planewave coordinates, RF data.
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem=Number of k points treated by this node (GS data)
!!  mk1mem=Number of k points treated by this node (RF data)
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum allowed value for npw (GS wfs)
!!  mpw1=maximum allowed value for npw1 (RF data)
!!  nband_rbz(nkpt_rbz*nsppol)=number of bands to be included in summation
!!   at each k point for each spin channel.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt_rbz=number of k points in the reduced Brillouin zone
!!  npwarr(nkpt_rbz)=number of planewaves and boundary planewaves at k points
!!  npwar1(nkpt_rbz)=number of planewaves and boundary planewaves at k+q points
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of symmetry elements in group (at least 1 for identity)
!!  occ_rbz(mband*nkpt_rbz*nsppol)=occupation numbers for each band
!!   (usually 2.0) at each k point of the reduced Brillouin zone
!!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  rprimd(3,3)=dimensional real space primitive translations
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrel(3,3,nsym)=symmetry matrices in real space (integers)
!!  ucvol=unit cell volume (Bohr**3)
!!  wtk_rbz(nkpt_rbz)=k point weights (they sum to 1.0).
!!
!! OUTPUT
!!  rhog1(2,nfft)=total electron density in G space
!!  rhor1(cplex*nfft,nspden)=electron density in r space
!!   (if spin polarized, array contains total density in first half and
!!    spin-up density in second half)
!!
!! PARENTS
!!      dfpt_looppert
!!
!! CHILDREN
!!      cg_zcopy,fftpac,fourwf,sphereboundary,symrhg,timab,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_mkrho(cg,cg1,cplex,gprimd,irrzon,istwfk_rbz,&
& kg,kg1,mband,mgfft,mkmem,mk1mem,mpi_enreg,mpw,mpw1,nband_rbz,&
& nfft,ngfft,nkpt_rbz,npwarr,npwar1,nspden,nspinor,nsppol,nsym,&
& occ_rbz,paral_kgb,phnons,rhog1,rhor1,rprimd,symafm,symrel,ucvol,wtk_rbz)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_cgtools
 use m_xmpi

 use m_io_tools,  only : get_unit, iomode_from_fname

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_mkrho'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_52_fft_mpi_noabirule
 use interfaces_53_ffts
 use interfaces_67_common
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,mband,mgfft,mk1mem,mkmem,mpw,mpw1,nfft,nkpt_rbz
 integer,intent(in) :: nspden,nspinor,nsppol,nsym,paral_kgb
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))
 integer,intent(in) :: istwfk_rbz(nkpt_rbz),kg(3,mpw*mkmem),kg1(3,mpw1*mk1mem)
 integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol),ngfft(18),npwar1(nkpt_rbz)
 integer,intent(in) :: npwarr(nkpt_rbz),symafm(nsym),symrel(3,3,nsym)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol),gprimd(3,3)
 real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))
 real(dp),intent(in) :: rprimd(3,3),wtk_rbz(nkpt_rbz)
 real(dp),intent(out) :: rhog1(2,nfft),rhor1(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourwf7=7,tim_rwwf15=15
 integer,save :: nskip=0
 integer :: bdtot_index,i1,i2,i3,iband,icg,icg1,ierr,ifft,ikg,ptr
 integer :: ikg1,ikpt,ispden,ispinor,isppol,istwf_k,ptr1,ptr2
 integer :: me,n1,n2,n3,n4,n5,n6,nband_k,npw1_k
 integer :: npw_k,spaceworld
 real(dp) :: im0,im1,re0,re1,weight
 real(dp) :: im0_up,im1_up,re0_up,re1_up,im0_down,im1_down,re0_down,re1_down
 character(len=500) :: message
!arrays
 integer,allocatable :: gbound(:,:),gbound1(:,:),kg1_k(:,:)
 integer,allocatable :: kg_k(:,:)
 real(dp) :: tsec(2)
 real(dp),allocatable,target :: cwavef(:,:),cwavef1(:,:)
 real(dp),allocatable :: dummy(:,:),rhoaug(:,:,:,:)
 real(dp),allocatable :: rhoaug1(:,:,:,:),wfraug(:,:,:,:),wfraug1(:,:,:,:)
 real(dp),allocatable :: wfraug1_up(:,:,:,:),wfraug1_down(:,:,:,:)
 real(dp),allocatable :: wfraug_up(:,:,:,:),wfraug_down(:,:,:,:)
 real(dp),allocatable :: cwave0_up(:,:),cwave0_down(:,:),cwave1_up(:,:),cwave1_down(:,:)

! *************************************************************************

!DBG_ENTER("COLL")

 if(nspden==4)then
!  NOTE : see mkrho for the modifications needed for non-collinear treatment
   write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&   ' dfpt_mkrho : WARNING -',ch10,&
&   '  Linear-response calculations are under construction with nspden=4',ch10,&
&   ' Action : modify value of nspden in input file unless you know what you are doing.'
!   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

!Init spaceworld
 spaceworld=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt

!zero the charge density array in real space
!$OMP PARALLEL DO 
 do ispden=1,nspden
   do ifft=1,cplex*nfft
     rhor1(ifft,ispden)=zero
   end do
 end do

!start loop over spin and k points
 bdtot_index=0; icg=0; icg1=0

 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 n4=ngfft(4); n5=ngfft(5); n6=ngfft(6) !n4,n5,n6 are FFT dimensions, modified to avoid cache trashing

!Note that the dimensioning of cwavef and cwavef1 does not include nspinor
 ABI_ALLOCATE(cwavef,(2,mpw))
 ABI_ALLOCATE(cwavef1,(2,mpw1))
!Actually, rhoaug is not needed, except for strong dimensioning requirement
 ABI_ALLOCATE(dummy,(2,1))
 ABI_ALLOCATE(rhoaug,(n4,n5,n6,nspinor**2))
 ABI_ALLOCATE(rhoaug1,(cplex*n4,n5,n6,nspinor**2))
 ABI_ALLOCATE(wfraug,(2,n4,n5,n6))
 ABI_ALLOCATE(wfraug1,(2,n4,n5,n6))

! EB FR Separate collinear and non-collinear magnetism
 if (nspden /= 4) then  ! EB FR nspden check
   do isppol=1,nsppol

     ikg=0; ikg1=0

     rhoaug1(:,:,:,:)=zero

     do ikpt=1,nkpt_rbz

       nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
       istwf_k=istwfk_rbz(ikpt)
       npw_k=npwarr(ikpt)
       npw1_k=npwar1(ikpt)

       if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
         bdtot_index=bdtot_index+nband_k
         cycle
       end if
       
       ABI_ALLOCATE(gbound,(2*mgfft+8,2))
       ABI_ALLOCATE(kg_k,(3,npw_k))
       ABI_ALLOCATE(gbound1,(2*mgfft+8,2))
       ABI_ALLOCATE(kg1_k,(3,npw1_k))

       kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
       call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)

       kg1_k(:,1:npw1_k)=kg1(:,1+ikg1:npw1_k+ikg1)
       call sphereboundary(gbound1,istwf_k,kg1_k,mgfft,npw1_k)

!    Loop over bands to fft and square for rho(r)
       do iband=1,nband_k
         if (mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me) cycle
!      Only treat occupied states
         if (abs(occ_rbz(iband+bdtot_index))>tol8) then
!        Treat separately the two spinor components
           do ispinor=1,nspinor
!          Obtain Fourier transform in fft box and accumulate the density
             ptr = 1 + (ispinor-1)*npw_k + (iband-1)*npw_k*nspinor + icg
             call cg_zcopy(npw_k, cg(1,ptr), cwavef)

!      In these two calls, rhoaug, rhoaug1 and weight are dummy variables, and are not modified
             call fourwf(1,rhoaug,cwavef,dummy,wfraug,gbound,gbound,&
&             istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,ngfft,npw_k,1,n4,n5,n6,0,paral_kgb,tim_fourwf7,weight,weight)

! TODO: here ispinor should be ispinorp to get full matrix and nspden 4
             ptr = 1 + (ispinor-1)*npw1_k + (iband-1)*npw1_k*nspinor + icg1
             call cg_zcopy(npw1_k, cg1(1,ptr), cwavef1)

             call fourwf(cplex,rhoaug1,cwavef1,dummy,wfraug1,gbound1,gbound1,&
&             istwf_k,kg1_k,kg1_k,mgfft,mpi_enreg,1,ngfft,npw1_k,1,n4,n5,n6,0,&
&             paral_kgb,tim_fourwf7,weight,weight)

!          Compute the weight, note that the factor 2 is
!          not the spin factor (see Eq.44 of PRB55,10337 (1997))
             weight=two*occ_rbz(iband+bdtot_index)*wtk_rbz(ikpt)/ucvol

!          Accumulate density
             if(cplex==2)then
!$OMP PARALLEL DO PRIVATE(im0,im1,re0,re1) 
               do i3=1,n3
                 do i2=1,n2
                   do i1=1,n1
                     re0=wfraug(1,i1,i2,i3) ; im0=wfraug(2,i1,i2,i3)
                     re1=wfraug1(1,i1,i2,i3); im1=wfraug1(2,i1,i2,i3)
                     rhoaug1(2*i1-1,i2,i3,1)=rhoaug1(2*i1-1,i2,i3,1)+weight*(re0*re1+im0*im1)
                     rhoaug1(2*i1  ,i2,i3,1)=rhoaug1(2*i1  ,i2,i3,1)+weight*(re0*im1-im0*re1)
                   end do
                 end do
               end do
             else
!$OMP PARALLEL DO 
               do i3=1,n3
                 do i2=1,n2
                   do i1=1,n1
                     rhoaug1(i1,i2,i3,1)=rhoaug1(i1,i2,i3,1)+&
&                     weight*( wfraug(1,i1,i2,i3)*wfraug1(1,i1,i2,i3) + wfraug(2,i1,i2,i3)*wfraug1(2,i1,i2,i3)  )
                   end do
                 end do
               end do
             end if ! cplex
           end do ! ispinor
         else !abs(occ_rbz(iband+bdtot_index))>tol8  
           nskip=nskip+1  ! if the state is not occupied. Accumulate the number of one-way 3D ffts skipped
         end if ! abs(occ_rbz(iband+bdtot_index))>tol8
       end do ! iband

       ABI_DEALLOCATE(gbound)
       ABI_DEALLOCATE(kg_k)
       ABI_DEALLOCATE(gbound1)
       ABI_DEALLOCATE(kg1_k)

       bdtot_index=bdtot_index+nband_k

       icg=icg+npw_k*nband_k
       ikg=ikg+npw_k
       icg1=icg1+npw1_k*nband_k
       ikg1=ikg1+npw1_k

     end do ! ikpt

     if (xmpi_paral==0) then !  Write the number of one-way 3D ffts skipped until now
       write(message,'(a,i8)')' mkrho3 : number of one-way 3D ffts skipped in mkrho3 until now =',nskip
       call wrtout(std_out,message,'PERS')
     end if

! TODO: if n+magnetization basis is used for the density, need to rotate rhoaug1 to that now, before packing into rhor1

!  Transfer density on augmented fft grid to normal fft grid in real space
!  Take also into account the spin, to place it correctly in rhor1.
!  Note the use of cplex
     call fftpac(isppol,mpi_enreg,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,ngfft,rhor1,rhoaug1,1)

   end do ! loop over isppol spins

 else ! nspden = 4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Part added for the non collinear magnetism
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The same lines of code are in 72_response/accrho3.F90
! TODO: merge these lines in a single routine??!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ikg=0; ikg1=0

   rhoaug1(:,:,:,:)=zero

   do ikpt=1,nkpt_rbz

     nband_k=nband_rbz(ikpt)
     istwf_k=istwfk_rbz(ikpt)
     npw_k=npwarr(ikpt)
     npw1_k=npwar1(ikpt)

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,1,me)) then
       bdtot_index=bdtot_index+nband_k
       cycle
     end if
     
     ABI_ALLOCATE(gbound,(2*mgfft+8,2))
     ABI_ALLOCATE(kg_k,(3,npw_k))
     ABI_ALLOCATE(gbound1,(2*mgfft+8,2))
     ABI_ALLOCATE(kg1_k,(3,npw1_k))

     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)

     kg1_k(:,1:npw1_k)=kg1(:,1+ikg1:npw1_k+ikg1)
     call sphereboundary(gbound1,istwf_k,kg1_k,mgfft,npw1_k)

!    Loop over bands to fft and square for rho(r)
     do iband=1,nband_k

       if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband,iband,1,me)) cycle

!      Only treat occupied states
       if (abs(occ_rbz(iband+bdtot_index))>tol8) then

! Build the four components of rho. We use only norm quantities and, so fourwf.

         ABI_ALLOCATE(wfraug_up,(2,n4,n5,n6))
         ABI_ALLOCATE(wfraug_down,(2,n4,n5,n6))
         ABI_ALLOCATE(wfraug1_up,(2,n4,n5,n6))
         ABI_ALLOCATE(wfraug1_down,(2,n4,n5,n6))
         ABI_ALLOCATE(cwave0_up,(2,mpw))
         ABI_ALLOCATE(cwave0_down,(2,mpw))
         ABI_ALLOCATE(cwave1_up,(2,mpw1))
         ABI_ALLOCATE(cwave1_down,(2,mpw1))

! EB FR build spinorial wavefunctions
! Obtain Fourier transform in fft box and accumulate the density
! EB FR How do we manage the following lines for non-collinear????
! zero order up and down spins
         ptr1 = 1 + (iband-1)*npw_k*nspinor + icg
         call cg_zcopy(npw_k, cg(1,ptr1), cwave0_up)
         ptr2 = 1 + npw_k + (iband-1)*npw_k*nspinor + icg
         call cg_zcopy(npw_k, cg(1,ptr2), cwave0_down)
! first order up and down spins
         ptr1 = 1 + (iband-1)*npw1_k*nspinor + icg1
         call cg_zcopy(npw1_k, cg1(1,ptr1), cwave1_up)
         ptr2 = 1 + npw1_k + (iband-1)*npw1_k*nspinor + icg1
         call cg_zcopy(npw1_k, cg1(1,ptr2), cwave1_down)

! TODO: here ispinor should be ispinorp to get full matrix and nspden 4
!        ptr = 1 + (ispinor-1)*npw1_k + (iband-1)*npw1_k*nspinor + icg1
!        call cg_zcopy(npw1_k, cg1(1,ptr), cwavef1)

! EB FR lines to be managed (?????)

! zero order
!        cwave0_up => cwavef(:,1:npw_k)
!        cwave0_down => cwavef(:,1+npw_k:2*npw_k)
! first order
!        cwave1_up => cwavef1(:,1:npw_k)
!        cwave1_down => cwavef1(:,1+npw_k:2*npw_k)

!    The factor 2 is not the spin factor (see Eq.44 of PRB55,10337 (1997)??)
         weight=two*occ_rbz(iband+bdtot_index)*wtk_rbz(ikpt)/ucvol
!density components 
!GS wfk Fourrier Tranform 
! EB FR in the fourwf calls rhoaug(:,:,:,2) is a dummy argument
         call fourwf(1,rhoaug(:,:,:,2),cwave0_up,dummy,wfraug_up,gbound,gbound,istwf_k,kg_k,kg_k,&
&         mgfft,mpi_enreg,1,ngfft,npw_k,1,n4,n5,n6,&
&         0,paral_kgb,tim_fourwf7,weight,weight)
         call fourwf(1,rhoaug(:,:,:,2),cwave0_down,dummy,wfraug_down,gbound,gbound,istwf_k,kg_k,kg_k,&
&         mgfft,mpi_enreg,1,ngfft,npw_k,1,n4,n5,n6,&
&         0,paral_kgb,tim_fourwf7,weight,weight)
 !1st order wfk Fourrier Transform
         call fourwf(1,rhoaug1(:,:,:,2),cwave1_up,dummy,wfraug1_up,gbound,gbound,istwf_k,kg_k,kg_k,&
&         mgfft,mpi_enreg,1,ngfft,npw_k,1,n4,n5,n6,&
&         0,paral_kgb,tim_fourwf7,weight,weight)
         call fourwf(1,rhoaug1(:,:,:,2),cwave1_down,dummy,wfraug1_down,gbound,gbound,istwf_k,kg_k,kg_k,&
&         mgfft,mpi_enreg,1,ngfft,npw_k,1,n4,n5,n6,&
&         0,paral_kgb,tim_fourwf7,weight,weight)

!    Accumulate 1st-order density (x component)
         if (cplex==2) then
           do i3=1,n3
             do i2=1,n2
               do i1=1,n1
                 re0_up=wfraug_up(1,i1,i2,i3)  ;     im0_up=wfraug_up(2,i1,i2,i3)
                 re1_up=wfraug1_up(1,i1,i2,i3) ;     im1_up=wfraug1_up(2,i1,i2,i3)
                 re0_down=wfraug_down(1,i1,i2,i3)  ; im0_down=wfraug_down(2,i1,i2,i3)
                 re1_down=wfraug1_down(1,i1,i2,i3) ; im1_down=wfraug1_down(2,i1,i2,i3)
             rhoaug1(2*i1-1,i2,i3,1)=rhoaug1(2*i1-1,i2,i3,1)+weight*(re0_up*re1_up+im0_up*im1_up) !n_upup
             rhoaug1(2*i1  ,i2,i3,1)=zero ! imag part of n_upup at k
             rhoaug1(2*i1-1,i2,i3,4)=rhoaug1(2*i1-1,i2,i3,4)+weight*(re0_down*re1_down+im0_down*im1_down) ! n_dndn
             rhoaug1(2*i1  ,i2,i3,4)=zero ! imag part of n_dndn at k
             rhoaug1(2*i1-1,i2,i3,2)=rhoaug1(2*i1-1,i2,i3,2)+weight*(re1_up*re0_down+re0_up*re1_down &
&                             +im0_up*im1_down+im0_down*im1_up)! mx; the factor two is inside weight
             rhoaug1(2*i1  ,i2,i3,2)=zero ! imag part of mx
             rhoaug1(2*i1-1,i2,i3,3)=rhoaug1(2*i1-1,i2,i3,3)+weight*(re1_up*im0_down-im1_up*re0_down &
&                             +re0_up*im1_down-im0_up*re1_down) ! my; the factor two is inside weight
             rhoaug1(2*i1  ,i2,i3,3)=zero ! imag part of my at k
               end do
             end do
           end do
         else
           do i3=1,n3
             do i2=1,n2
               do i1=1,n1
                 re0_up=wfraug_up(1,i1,i2,i3)  ;     im0_up=wfraug_up(2,i1,i2,i3)
                 re1_up=wfraug1_up(1,i1,i2,i3) ;     im1_up=wfraug1_up(2,i1,i2,i3)
                 re0_down=wfraug_down(1,i1,i2,i3)  ; im0_down=wfraug_down(2,i1,i2,i3)
                 re1_down=wfraug1_down(1,i1,i2,i3) ; im1_down=wfraug1_down(2,i1,i2,i3)
             rhoaug1(i1,i2,i3,1)=rhoaug1(i1,i2,i3,1)+weight*(re0_up*re1_up+im0_up*im1_up) ! n_upup
             rhoaug1(i1,i2,i3,4)=rhoaug1(i1,i2,i3,4)+weight*(re0_down*re1_down+im0_down*im1_down) ! n_dndn
             rhoaug1(i1,i2,i3,2)=rhoaug1(i1,i2,i3,2)+weight*(re1_up*re0_down+re0_up*re1_down &
&                             +im0_up*im1_down+im0_down*im1_up) !mx; the factor two is inside weight
             rhoaug1(i1,i2,i3,3)=rhoaug1(i1,i2,i3,3)+weight*(re1_up*im0_down-im1_up*re0_down &
&                             +re0_up*im1_down-im0_up*re1_down) !my; the factor two is inside weight
               end do
             end do
           end do
         end if
         ABI_DEALLOCATE(wfraug_up)
         ABI_DEALLOCATE(wfraug_down)
         ABI_DEALLOCATE(wfraug1_up)
         ABI_DEALLOCATE(wfraug1_down)
         ABI_DEALLOCATE(cwave0_up)
         ABI_DEALLOCATE(cwave0_down)
         ABI_DEALLOCATE(cwave1_up)
         ABI_DEALLOCATE(cwave1_down)

       end if ! occupied states
     end do ! End loop on iband

     ABI_DEALLOCATE(gbound)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(gbound1)
     ABI_DEALLOCATE(kg1_k)

     bdtot_index=bdtot_index+nband_k

     icg=icg+npw_k*nband_k
     ikg=ikg+npw_k

     icg1=icg1+npw1_k*nband_k
     ikg1=ikg1+npw1_k

   end do ! End loop on ikpt


   if (xmpi_paral==0) then !  Write the number of one-way 3D ffts skipped until now
     write(message,'(a,i8)')' dfpt_mkrho : number of one-way 3D ffts skipped in mkrho3 until now =',nskip
     call wrtout(std_out,message,'PERS')
   end if

! TODO: if n+magnetization basis is used for the density, need to rotate rhoaug1 to that now, before packing into rhor1

!  Transfer density on augmented fft grid to normal fft grid in real space
!  Take also into account the spin, to place it correctly in rhor1.
!  Note the use of cplex
   call fftpac(1,mpi_enreg,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,ngfft,rhor1,rhoaug1,1)
 end if ! nspden /= 4

!if (xmpi_paral==1) then
!call timab(63,1,tsec)
!call wrtout(std_out,'dfpt_mkrho: loop on k-points and spins done in parallel','COLL')
!call xmpi_barrier(spaceworld)
!call timab(63,2,tsec)
!end if

 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(cwavef1)
 ABI_DEALLOCATE(dummy)
 ABI_DEALLOCATE(rhoaug)
 ABI_DEALLOCATE(rhoaug1)
 ABI_DEALLOCATE(wfraug)
 ABI_DEALLOCATE(wfraug1)

!Recreate full rhor1 on all proc.
 call timab(48,1,tsec)
 call timab(71,1,tsec)
 call xmpi_sum(rhor1,spaceworld,ierr)
 call timab(71,2,tsec)
 call timab(48,2,tsec)

 call symrhg(cplex,gprimd,irrzon,mpi_enreg,nfft,nfft,ngfft,nspden,nsppol,nsym,paral_kgb,phnons,&
& rhog1,rhor1,rprimd,symafm,symrel)

!We now have both rho(r) and rho(G), symmetrized, and if nsppol=2
!we also have the spin-up density, symmetrized, in rhor1(:,2).

!DBG_EXIT("COLL")

end subroutine dfpt_mkrho
!!***
