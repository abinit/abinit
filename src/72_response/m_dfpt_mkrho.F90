!!****m* ABINIT/m_dfpt_mkrho
!! NAME
!!  m_dfpt_mkrho
!!
!! FUNCTION
!! Compute RF charge density rho1(r) and rho1(G) in electrons/bohr**3
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, LSI, AR, MB, MT, SPr)
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

module m_dfpt_mkrho

 use defs_basis
 use m_abicore
 use m_errors
 use m_cgtools
 use m_xmpi


 use defs_abitypes, only : MPI_type
 use m_time,            only : timab
 use m_io_tools,        only : get_unit, iomode_from_fname
 use m_fftcore,         only : sphereboundary
 use m_fft,             only : fftpac, fourwf
 use m_spacepar,        only : symrhg
 use m_hamiltonian,     only : gs_hamiltonian_type
 use m_pawrhoij,        only : pawrhoij_type
 use m_pawcprj,         only : pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_paw_occupancies, only : pawaccrhoij
 use m_paral_atom,      only : get_my_atmtab
 use m_mpinfo,          only : proc_distrb_cycle
 use m_cgprj,           only : getcprj

 implicit none

 private
!!***

 public :: dfpt_mkrho
 public :: dfpt_accrho
!!***

contains
!!***

!!****f* ABINIT/dfpt_mkrho
!! NAME
!! dfpt_mkrho
!!
!! FUNCTION
!! Compute RF charge density rho1(r) and rho1(G) in electrons/bohr**3
!! from input RF and GS wavefunctions, band occupations, and k point weights.
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=wf in G space
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=first-order wf in G space
!!  cplex=1 if rhor1 is real, 2 if rhor1 is complex
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates, GS data.
!!  kg1(3,mpw1*mkmem1)=reduced planewave coordinates, RF data.
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem=Number of k points treated by this node (GS data)
!!  mk1mem=Number of k points treated by this node (RF data)
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum allowed value for npw (GS wfs)
!!  mpw1=maximum allowed value for npw1 (RF data)
!!  nband_rbz(nkpt_rbz*nsppol)=number of bands to be included in summation
!!   at each k point for each spin channel.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/variables/vargs.htm#ngfft
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
!!      m_dfpt_looppert
!!
!! CHILDREN
!!      fourwf,get_my_atmtab,getcprj,pawaccrhoij,pawcprj_alloc,pawcprj_free
!!
!! SOURCE

subroutine dfpt_mkrho(cg,cg1,cplex,gprimd,irrzon,istwfk_rbz,&
& kg,kg1,mband,mgfft,mkmem,mk1mem,mpi_enreg,mpw,mpw1,nband_rbz,&
& nfft,ngfft,nkpt_rbz,npwarr,npwar1,nspden,nspinor,nsppol,nsym,&
& occ_rbz,phnons,rhog1,rhor1,rprimd,symafm,symrel,tnons,ucvol,wtk_rbz)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,mband,mgfft,mk1mem,mkmem,mpw,mpw1,nfft,nkpt_rbz
 integer,intent(in) :: nspden,nspinor,nsppol,nsym
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))
 integer,intent(in) :: istwfk_rbz(nkpt_rbz),kg(3,mpw*mkmem),kg1(3,mpw1*mk1mem)
 integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol),ngfft(18),npwar1(nkpt_rbz)
 integer,intent(in) :: npwarr(nkpt_rbz),symafm(nsym),symrel(3,3,nsym)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol),gprimd(3,3)
 real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))
 real(dp),intent(in) :: rprimd(3,3),tnons(3,nsym)
 real(dp),intent(in) :: wtk_rbz(nkpt_rbz)
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
!  NOTE: see mkrho for the modifications needed for non-collinear treatment
   write(message, '(3a)' )&
    ' Linear-response calculations are under construction with nspden=4',ch10,&
    ' Action: modify value of nspden in input file unless you know what you are doing.'
   ABI_WARNING(message)
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
&             istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,ngfft,npw_k,1,n4,n5,n6,0,tim_fourwf7,weight,weight)

! TODO: here ispinor should be ispinorp to get full matrix and nspden 4
             ptr = 1 + (ispinor-1)*npw1_k + (iband-1)*npw1_k*nspinor + icg1
             call cg_zcopy(npw1_k, cg1(1,ptr), cwavef1)

             call fourwf(cplex,rhoaug1,cwavef1,dummy,wfraug1,gbound1,gbound1,&
&             istwf_k,kg1_k,kg1_k,mgfft,mpi_enreg,1,ngfft,npw1_k,1,n4,n5,n6,0,&
&             tim_fourwf7,weight,weight)

!          Compute the weight, note that the factor 2 is
!          not the spin factor (see Eq.44 of PRB55,10337 (1997) [[cite:Gonze1997]])
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

!    The factor 2 is not the spin factor (see Eq.44 of PRB55,10337 (1997) ?? [[cite:Gonze1997]])
         weight=two*occ_rbz(iband+bdtot_index)*wtk_rbz(ikpt)/ucvol
!density components
!GS wfk Fourrier Tranform
! EB FR in the fourwf calls rhoaug(:,:,:,2) is a dummy argument
         call fourwf(1,rhoaug(:,:,:,2),cwave0_up,dummy,wfraug_up,gbound,gbound,istwf_k,kg_k,kg_k,&
&         mgfft,mpi_enreg,1,ngfft,npw_k,1,n4,n5,n6,&
&         0,tim_fourwf7,weight,weight)
         call fourwf(1,rhoaug(:,:,:,2),cwave0_down,dummy,wfraug_down,gbound,gbound,istwf_k,kg_k,kg_k,&
&         mgfft,mpi_enreg,1,ngfft,npw_k,1,n4,n5,n6,&
&         0,tim_fourwf7,weight,weight)
 !1st order wfk Fourrier Transform
         call fourwf(1,rhoaug1(:,:,:,2),cwave1_up,dummy,wfraug1_up,gbound,gbound,istwf_k,kg_k,kg_k,&
&         mgfft,mpi_enreg,1,ngfft,npw_k,1,n4,n5,n6,&
&         0,tim_fourwf7,weight,weight)
         call fourwf(1,rhoaug1(:,:,:,2),cwave1_down,dummy,wfraug1_down,gbound,gbound,istwf_k,kg_k,kg_k,&
&         mgfft,mpi_enreg,1,ngfft,npw_k,1,n4,n5,n6,&
&         0,tim_fourwf7,weight,weight)

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
&                 +im0_up*im1_down+im0_down*im1_up) ! mx; the factor two is inside weight
                 rhoaug1(2*i1  ,i2,i3,2)=zero ! imag part of mx
                 rhoaug1(2*i1-1,i2,i3,3)=rhoaug1(2*i1-1,i2,i3,3)+weight*(re1_up*im0_down-im1_up*re0_down &
&                 +re0_up*im1_down-im0_up*re1_down) ! my; the factor two is inside weight
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
&                 +im0_up*im1_down+im0_down*im1_up) !mx; the factor two is inside weight
                 rhoaug1(i1,i2,i3,3)=rhoaug1(i1,i2,i3,3)+weight*(re1_up*im0_down-im1_up*re0_down &
&                 +re0_up*im1_down-im0_up*re1_down) !my; the factor two is inside weight
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

 call symrhg(cplex,gprimd,irrzon,mpi_enreg,nfft,nfft,ngfft,nspden,nsppol,nsym,phnons,&
             rhog1,rhor1,rprimd,symafm,symrel,tnons)

!We now have both rho(r) and rho(G), symmetrized, and if nsppol=2
!we also have the spin-up density, symmetrized, in rhor1(:,2).

!DBG_EXIT("COLL")

end subroutine dfpt_mkrho
!!***

!!****f* ABINIT/dfpt_accrho
!! NAME
!! dfpt_accrho
!!
!! FUNCTION
!! Response function calculation only:
!!  Accumulate contribution to first-order density due to current (k,band)
!!  Also accumulate zero-order potential part of the 2nd-order total energy (if needed)
!!
!! INPUTS
!!  cplex=1 if 1st-order density is real, 2 if 1st-order density is complex
!!  cwave0(2,npw*nspinor)=GS wavefunction at k, in reciprocal space
!!  cwave1(2,npw1*nspinor)=1st-order wavefunction at k,q, in reciprocal space
!!  cwavef(2,npw1*nspinor)=1st-order wavefunction at k,q, in reciprocal space, without correction due to occupation change
!!  cwaveprj0(natom,nspinor*usecprj)= GS wave function at k projected with nl projectors
!!  cwaveprj1(natom,nspinor*usecprj)= 1st-order wave function at k,q projected with nl projectors
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  iband=index of current band
!!  idir=direction of the current perturbation
!!  ipert=type of the perturbation
!!  isppol=1 index of current spin component
!!  kptopt=option for the generation of k points
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell
!!  nband_k=number of bands at this k point for that spin polarization
!!  ncpgr=number of gradients stored in cprj array (cprj=<p_i|Cnk>)
!!  npw_k=number of planewaves in basis sphere at k
!!  npw1_k=number of planewaves in basis sphere at k+q
!!  nspinor=number of spinorial components of the wavefunctions
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  option= 1: accumulate 1st-order density,
!!          2: accumulate 0-order potential part of the 2nd-order total energy
!!          3: accumulate both
!!  tim_fourwf= timing code for fourwf (5 from dfpt_vtowfk, 18 from dfpt_nstwf)
!!  wf_corrected=flag put to 1 if cwave1 is different from cwavef (if there is a contribution from occ. change)
!!  wtk_k=weight assigned to the k point.
!!
!! OUTPUT
!!  ====== if option=2 or option=3 =====
!!  eloc0_k=zero-order local contribution to 2nd-order total energy for current band and k
!!
!! SIDE EFFECTS
!!  ====== if option=1 or option=3 =====
!!    rhoaug1(cplex*n4,n5,n6,nvloc)= density in electrons/bohr**3,
!!    ==== if gs_hamkq%usepaw=1 =====
!!    pawrhoij1(natom) <type(pawrhoij_type)>= 1st-order paw rhoij occupancies and related data
!!                                            (cumulative, so input as well as output)
!!
!! NOTES
!!  In this part of the treatment of one band, one has to
!!  perform Fourier transforms, and to treat separately the
!!  two spinorial components of the wavefunction.
!!  Was part of dfpt_vtowfk before.
!!
!! PARENTS
!!      m_dfpt_nstwf,m_dfpt_scfcv,m_dfpt_vtowfk
!!
!! CHILDREN
!!      fourwf,get_my_atmtab,getcprj,pawaccrhoij,pawcprj_alloc,pawcprj_free
!!
!! SOURCE

subroutine dfpt_accrho(cplex,cwave0,cwave1,cwavef,cwaveprj0,cwaveprj1,&
&                  eloc0_k,gs_hamkq,iband,idir,ipert,isppol,kptopt,&
&                  mpi_enreg,natom,nband_k,ncpgr,npw_k,npw1_k,nspinor,occ_k,&
&                  option,pawrhoij1,rhoaug1,tim_fourwf,wf_corrected,&
&                  wtk_k,comm_atom,mpi_atmtab)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,iband,idir,ipert,isppol,kptopt,natom,nband_k
 integer,intent(in) :: ncpgr,npw_k,npw1_k,nspinor,option,tim_fourwf,wf_corrected
 integer,optional,intent(in) :: comm_atom
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: wtk_k
 real(dp),intent(out) :: eloc0_k
 type(gs_hamiltonian_type),intent(inout),target :: gs_hamkq
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 real(dp),intent(in),target :: cwave0(2,npw_k*nspinor),cwave1(2,npw1_k*nspinor),cwavef(2,npw1_k*nspinor)
 real(dp),intent(in) :: occ_k(nband_k)
 real(dp),intent(inout) :: rhoaug1(cplex*gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,gs_hamkq%nvloc)
 type(pawcprj_type),intent(in) :: cwaveprj0(natom,nspinor*gs_hamkq%usecprj)
 type(pawcprj_type),intent(in) :: cwaveprj1(natom,nspinor*gs_hamkq%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij1(:)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=14
 integer :: choice,cplex_cprj,i1,i2,i3,ispinor,my_comm_atom,my_natom,n1,n2,n3,option_rhoij
 logical :: my_atmtab_allocated,paral_atom
 logical :: usetimerev
 real(dp) :: im0,im1,re0,re1,valuer,diag,offdiag,weight
 real(dp) :: im0_up,im1_up,re0_up,re1_up,im0_down,im1_down,re0_down,re1_down
!arrays
 integer,pointer :: my_atmtab(:)
 real(dp) :: dummy(2,1)
 real(dp),allocatable :: rhoaug(:,:,:,:),wfraug(:,:,:,:),wfraug1(:,:,:,:)
 real(dp),allocatable :: wfraug1_up(:,:,:,:),wfraug1_down(:,:,:,:)
 real(dp),allocatable :: wfraug_up(:,:,:,:),wfraug_down(:,:,:,:)
 real(dp),pointer :: cwavef_sp(:,:),cwavef_up(:,:),cwavef_down(:,:)
 real(dp),pointer :: cwave0_up(:,:),cwave0_down(:,:),cwave1_up(:,:),cwave1_down(:,:)
 real(dp),pointer :: vlocal(:,:,:,:)=>null()
 type(pawcprj_type),allocatable :: cwaveprj_tmp(:,:)

! *********************************************************************
 DBG_ENTER("COLL")

 if (option/=1.and.option/=2.and.option/=3) return

!Initializations
 ABI_ALLOCATE(rhoaug,(gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,gs_hamkq%nvloc))
 n1=gs_hamkq%ngfft(1);n2=gs_hamkq%ngfft(2);n3=gs_hamkq%ngfft(3)
 if (option==2.or.option==3) eloc0_k=zero
 if (option==2.or.option==3) vlocal => gs_hamkq%vlocal

!Loop on spinorial components
! TODO : double loop on spinors for full rhoaug1 matrix if nspden =4
 if (gs_hamkq%nvloc/=4) then  ! see later EB FR
   ABI_ALLOCATE(wfraug1,(2,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))

   do ispinor=1,nspinor

!  Part devoted to the accumulation of the 0-order potential part of the 2nd-order total energy
!  --------------------------------------------------------------------------------------------

!  Fourier transform of cwavef. Here, rhoaug is a dummy variable.
     if (wf_corrected==0.or.option==2.or.option==3) then
       if (ispinor==1) then
         cwavef_sp => cwavef(:,1:npw1_k)
       else
         cwavef_sp => cwavef(:,1+npw1_k:2*npw1_k)
       end if
       !make an inverse FFT from cwavef_sp to wfraug1
       call fourwf(cplex,rhoaug,cwavef_sp,dummy,wfraug1,gs_hamkq%gbound_kp,gs_hamkq%gbound_kp,&
&       gs_hamkq%istwf_k,gs_hamkq%kg_kp,gs_hamkq%kg_kp,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&       gs_hamkq%npw_kp,1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,tim_fourwf,&
&       weight,weight,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
       nullify(cwavef_sp)
     end if

!  Compute contribution of this band to zero-order potential part of the 2nd-order total energy
!  NB: this is spinor diagonal
     if (option==2.or.option==3) then
       valuer=zero
       do i3=1,n3
         do i2=1,n2
           do i1=1,n1
             valuer=valuer+vlocal(i1,i2,i3,1)*(wfraug1(1,i1,i2,i3)**2+wfraug1(2,i1,i2,i3)**2)
           end do
         end do
       end do

!    Local potential energy of this band
       eloc0_k=eloc0_k+two*valuer/dble(gs_hamkq%nfft)
     end if ! option

!  Part devoted to the accumulation of the 1st-order density
!  ---------------------------------------------------------
     if (option==1.or.option==3) then

!    Compute 1st-order WF in real space
!    One needs the Fourier transform of cwave1. However, only the one of
!    cwavef is available. If cwavef and cwave1 differs, this Fourier
!    transform must be computed. In both case the result is in wfraug1.
       if (wf_corrected==1) then
         if (ispinor==1) then
           cwavef_sp => cwave1(:,1:npw1_k)
         else
           cwavef_sp => cwave1(:,1+npw1_k:2*npw1_k)
         end if
         call fourwf(cplex,rhoaug,cwavef_sp,dummy,wfraug1,gs_hamkq%gbound_kp,gs_hamkq%gbound_kp,&
&         gs_hamkq%istwf_k,gs_hamkq%kg_kp,gs_hamkq%kg_kp,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&         gs_hamkq%npw_kp,1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,tim_fourwf,&
&         weight,weight,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
         nullify(cwavef_sp)
       end if

!    Compute 0-order WF in real space
! TODO: add loop over ispinor_prime here
       ABI_ALLOCATE(wfraug,(2,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))
       if (ispinor==1) then
         cwavef_sp => cwave0(:,1:npw_k)
       else
         cwavef_sp => cwave0(:,1+npw_k:2*npw_k)
       end if
       call fourwf(1,rhoaug,cwavef_sp,dummy,wfraug,gs_hamkq%gbound_k,gs_hamkq%gbound_k,&
&       gs_hamkq%istwf_k,gs_hamkq%kg_k,gs_hamkq%kg_k,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&       gs_hamkq%npw_k,1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,tim_fourwf,&
&       weight,weight,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
       nullify(cwavef_sp)

!    The factor 2 is not the spin factor (see Eq.44 of PRB55,10337 (1997) [[cite:Gonze1997]])
       weight=two*occ_k(iband)*wtk_k/gs_hamkq%ucvol
!    Accumulate 1st-order density
       if (cplex==2) then
         do i3=1,n3
           do i2=1,n2
             do i1=1,n1
               re0=wfraug(1,i1,i2,i3)  ; im0=wfraug(2,i1,i2,i3)
               re1=wfraug1(1,i1,i2,i3) ; im1=wfraug1(2,i1,i2,i3)
! TODO: check which terms (ispinor ispinorp) enter a given element of rhoaug1
               rhoaug1(2*i1-1,i2,i3,1)=rhoaug1(2*i1-1,i2,i3,1)+weight*(re0*re1+im0*im1)
               rhoaug1(2*i1  ,i2,i3,1)=rhoaug1(2*i1  ,i2,i3,1)+weight*(re0*im1-im0*re1)
             end do
           end do
         end do
       else
         do i3=1,n3
           do i2=1,n2
             do i1=1,n1
               rhoaug1(i1,i2,i3,1)=rhoaug1(i1,i2,i3,1) &
&               +weight*(wfraug(1,i1,i2,i3)*wfraug1(1,i1,i2,i3) &
&               +wfraug(2,i1,i2,i3)*wfraug1(2,i1,i2,i3))
             end do
           end do
         end do
       end if
       ABI_DEALLOCATE(wfraug)
     end if ! option

   end do ! Loop on spinorial components if nspden=1 or 2

   ABI_DEALLOCATE(wfraug1)
 else ! nvloc = 4
! The same lines of code are in 72_response/dfpt_mkrho.F90
! TODO merge these lines in a single routine??!!
   ABI_ALLOCATE(wfraug1_up,(2,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))
   ABI_ALLOCATE(wfraug1_down,(2,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))

   wfraug1_up(:,:,:,:)=zero
   wfraug1_down(:,:,:,:)=zero

!  Part devoted to the accumulation of the 0-order potential part of the 2nd-order total energy
!  --------------------------------------------------------------------------------------------

!  Fourier transform of cwavef. Here, rhoaug is a dummy variable.
   if (wf_corrected==0.or.option==2.or.option==3) then
     cwavef_up => cwavef(:,1:npw1_k) ! wfs up spin-polarized
     call fourwf(cplex,rhoaug(:,:,:,1),cwavef_up,dummy,wfraug1_up,gs_hamkq%gbound_kp,gs_hamkq%gbound_kp,&
&     gs_hamkq%istwf_k,gs_hamkq%kg_kp,gs_hamkq%kg_kp,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&     gs_hamkq%npw_kp,1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,tim_fourwf,&
&     weight,weight,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
     nullify(cwavef_up)

     cwavef_down => cwavef(:,1+npw1_k:2*npw1_k) ! wfs down spin-polarized
     call fourwf(cplex,rhoaug(:,:,:,1),cwavef_down,dummy,wfraug1_down,gs_hamkq%gbound_kp,gs_hamkq%gbound_kp,&
&     gs_hamkq%istwf_k,gs_hamkq%kg_kp,gs_hamkq%kg_kp,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&     gs_hamkq%npw_kp,1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,tim_fourwf,&
&     weight,weight,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
     nullify(cwavef_down)
   end if
   if (option==2.or.option==3) then
     valuer=zero
     diag=zero
     offdiag=zero
   ! EB FR 2nd term in Eq. 91 PRB52,1096 [[cite:Gonze1995]] for non-collinear magnetism
     do i3=1,n3
       do i2=1,n2
         do i1=1,n1
           diag=vlocal(i1,i2,i3,1)*(wfraug1_up(1,i1,i2,i3)**2+wfraug1_up(2,i1,i2,i3)**2)&
&           +vlocal(i1,i2,i3,2)*(wfraug1_down(1,i1,i2,i3)**2+wfraug1_down(2,i1,i2,i3)**2)
           offdiag=(two*vlocal(i1,i2,i3,3)*((wfraug1_up(1,i1,i2,i3)*wfraug1_down(1,i1,i2,i3))+&
&           (wfraug1_up(2,i1,i2,i3)*wfraug1_down(2,i1,i2,i3))))+&
&           (two*vlocal(i1,i2,i3,4)*((-wfraug1_down(2,i1,i2,i3)*wfraug1_up(1,i1,i2,i3))+&
&           (wfraug1_down(1,i1,i2,i3)*wfraug1_up(2,i1,i2,i3))))
           valuer=valuer+diag+offdiag
         end do
       end do
     end do
!    Local potential energy of this band
     eloc0_k=eloc0_k+two*valuer/dble(gs_hamkq%nfft)
   end if ! option

!  Part devoted to the accumulation of the 1st-order density
!  ---------------------------------------------------------

! first order
!
   if (option==1.or.option==3) then

     !SPr: condition on wf_corrected not to do FFTs of the same Bloch functions
     if (wf_corrected==1) then
       cwave1_up => cwave1(:,1:npw1_k)
       cwave1_down => cwave1(:,1+npw1_k:2*npw1_k)
       wfraug1_up(:,:,:,:)=zero
       wfraug1_down(:,:,:,:)=zero

       call fourwf(cplex,rhoaug(:,:,:,1),cwave1_up,dummy,wfraug1_up,gs_hamkq%gbound_kp,gs_hamkq%gbound_kp,&
&       gs_hamkq%istwf_k,gs_hamkq%kg_kp,gs_hamkq%kg_kp,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&       gs_hamkq%npw_kp,1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,tim_fourwf,&
&       weight,weight,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
       nullify(cwave1_up)

       call fourwf(cplex,rhoaug(:,:,:,1),cwave1_down,dummy,wfraug1_down,gs_hamkq%gbound_kp,gs_hamkq%gbound_kp,&
&       gs_hamkq%istwf_k,gs_hamkq%kg_kp,gs_hamkq%kg_kp,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&       gs_hamkq%npw_kp,1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,tim_fourwf,&
&       weight,weight,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
       nullify(cwave1_down)
     end if


! EB FR build spinorial wavefunctions
! zero order
     cwave0_up => cwave0(:,1:npw_k)
     cwave0_down => cwave0(:,1+npw_k:2*npw_k)
     ABI_ALLOCATE(wfraug_up,(2,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))
     ABI_ALLOCATE(wfraug_down,(2,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))
     wfraug_up(:,:,:,:)=zero
     wfraug_down(:,:,:,:)=zero
!
     !density components
     !GS wfk Fourrier Tranform
     ! EB FR in the fourwf calls rhoaug(:,:,:,2) is a dummy argument
     call fourwf(1,rhoaug(:,:,:,2),cwave0_up,dummy,wfraug_up,gs_hamkq%gbound_k,gs_hamkq%gbound_k,&
&     gs_hamkq%istwf_k,gs_hamkq%kg_k,gs_hamkq%kg_k,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&     gs_hamkq%npw_k,1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,tim_fourwf,&
&     weight,weight,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
     nullify(cwave0_up)
     call fourwf(1,rhoaug(:,:,:,2),cwave0_down,dummy,wfraug_down,gs_hamkq%gbound_k,gs_hamkq%gbound_k,&
&     gs_hamkq%istwf_k,gs_hamkq%kg_k,gs_hamkq%kg_k,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&     gs_hamkq%npw_k,1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,0,tim_fourwf,&
&     weight,weight,use_gpu_cuda=gs_hamkq%use_gpu_cuda)
     nullify(cwave0_down)
!    Accumulate 1st-order density (x component)
     re0_up=zero;im0_up=zero;re1_up=zero;im1_up=zero;re0_down=zero;im0_down=zero
     re1_down=zero;im1_down=zero
!    The factor 2 is not the spin factor (see Eq.44 of PRB55,10337 (1997) [[cite:Gonze1997]])
!    SPr: the following treatment with factor=2 is ok for perturbations not breaking the
!         time reversal symmetry of the Hamiltonian (due to Kramer's degeneracy) hence
!         not applicable for magnetic field perturbation (for phonons with SOC, H^(0) has
!         time reversal symmetry though). The formulas below are rectified in dfpt_scfcv
!         in case of broken time-reversal upon reconstructing rhor1_pq and rhor1_mq.
     weight=two*occ_k(iband)*wtk_k/gs_hamkq%ucvol
     if (cplex==2) then
       do i3=1,n3
         do i2=1,n2
           do i1=1,n1
             re0_up=wfraug_up(1,i1,i2,i3)  ;     im0_up=wfraug_up(2,i1,i2,i3)
             re1_up=wfraug1_up(1,i1,i2,i3) ;     im1_up=wfraug1_up(2,i1,i2,i3)
             re0_down=wfraug_down(1,i1,i2,i3)  ; im0_down=wfraug_down(2,i1,i2,i3)
             re1_down=wfraug1_down(1,i1,i2,i3) ; im1_down=wfraug1_down(2,i1,i2,i3)
             !SPr: in case of +q/-q calculation, the factor will be corrected later from dfpt_scfcv level
             !     along with the reconstruction of correct rhor1_{+q} and rhor1_{-q}
             !     here, rhoaug1_{sigma,sigma'} = \sum_{n,k} u1_{sigma} u0*_{sigma'} independent of the sign of q
             rhoaug1(2*i1-1,i2,i3,1)=rhoaug1(2*i1-1,i2,i3,1)+weight*(re0_up*re1_up+im0_up*im1_up) !n_upup
             rhoaug1(2*i1  ,i2,i3,1)=rhoaug1(2*i1  ,i2,i3,1)+weight*(re0_up*im1_up-im0_up*re1_up)
             rhoaug1(2*i1-1,i2,i3,4)=rhoaug1(2*i1-1,i2,i3,4)+weight*(re0_down*re1_down+im0_down*im1_down) ! n_dndn
             rhoaug1(2*i1  ,i2,i3,4)=rhoaug1(2*i1  ,i2,i3,4)+weight*(re0_down*im1_down-im0_down*re1_down)

             rhoaug1(2*i1-1,i2,i3,2)=rhoaug1(2*i1-1,i2,i3,2)+weight*(re1_up*re0_down+im1_up*im0_down)& !Re[m1x]
&            +weight*(re1_down*re0_up+im1_down*im0_up)
             rhoaug1(2*i1  ,i2,i3,2)=rhoaug1(2*i1  ,i2,i3,2)+weight*(-re1_up*im0_down+im1_up*re0_down)& !Im[m1x]
&            +weight*(-re1_down*im0_up+im1_down*re0_up)

             rhoaug1(2*i1-1,i2,i3,3)=rhoaug1(2*i1-1,i2,i3,3)+weight*(+re1_up*im0_down-im1_up*re0_down)& !Re[m1y]
&            +weight*(-re1_down*im0_up+im1_down*re0_up)
             rhoaug1(2*i1  ,i2,i3,3)=rhoaug1(2*i1  ,i2,i3,3)+weight*(+re1_up*re0_down+im1_up*im0_down)& !Im[m1y]
&            +weight*(-re1_down*re0_up-im1_down*im0_up)
           end do
         end do
       end do
     else !cplex
       re0_up=zero;im0_up=zero;re1_up=zero;im1_up=zero;re0_down=zero;im0_down=zero
       re1_down=zero;im1_down=zero
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
&             +im0_up*im1_down+im0_down*im1_up) !mx; the factor two is inside weight
             rhoaug1(i1,i2,i3,3)=rhoaug1(i1,i2,i3,3)+weight*(re1_up*im0_down-im1_up*re0_down &
&             +re0_up*im1_down-im0_up*re1_down) !my; the factor two is inside weight
           end do
         end do
       end do
     end if !cplex

     ABI_DEALLOCATE(wfraug_up)
     ABI_DEALLOCATE(wfraug_down)

   end if ! option

   ABI_DEALLOCATE(wfraug1_up)
   ABI_DEALLOCATE(wfraug1_down)

 end if ! nvloc /= 4

 ABI_DEALLOCATE(rhoaug)

!Part devoted to the accumulation of the 1st-order occupation matrix in PAW case
! TODO: parse for more nspden 4 dependencies on spinors
! EB FR CHECK: to be modified for non-collinear?????
!-------------------------------------------------------------------------------

 if ((option==1.or.option==3).and.gs_hamkq%usepaw==1) then

!  Set up parallelism over atoms
   my_natom=natom; if(gs_hamkq%usepaw==1) my_natom=size(pawrhoij1)
   paral_atom=(present(comm_atom).and.(my_natom/=natom))
   my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
   nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
   call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

   cplex_cprj=2;if (gs_hamkq%istwf_k>1) cplex_cprj=1
   option_rhoij=2;usetimerev=(kptopt>0.and.kptopt<3)

   if (gs_hamkq%usecprj==1) then
     call pawaccrhoij(gs_hamkq%atindx,cplex_cprj,cwaveprj0,cwaveprj1,ipert,isppol,&
&     my_natom,natom,nspinor,occ_k(iband),option_rhoij,pawrhoij1,usetimerev,wtk_k,&
&     comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
   else
     ABI_DATATYPE_ALLOCATE(cwaveprj_tmp,(natom,nspinor))
     call pawcprj_alloc(cwaveprj_tmp,ncpgr,gs_hamkq%dimcprj)
     choice=2
     call getcprj(choice,0,cwave0,cwaveprj_tmp,&
&     gs_hamkq%ffnl_k,idir,gs_hamkq%indlmn,gs_hamkq%istwf_k,&
&     gs_hamkq%kg_k,gs_hamkq%kpg_k,gs_hamkq%kpt_k,gs_hamkq%lmnmax,&
&     gs_hamkq%mgfft,mpi_enreg,gs_hamkq%natom,gs_hamkq%nattyp,gs_hamkq%ngfft,&
&     gs_hamkq%nloalg,gs_hamkq%npw_k,gs_hamkq%nspinor,gs_hamkq%ntypat,gs_hamkq%phkxred,&
&     gs_hamkq%ph1d,gs_hamkq%ph3d_k,gs_hamkq%ucvol,gs_hamkq%useylm)
     call pawaccrhoij(gs_hamkq%atindx,cplex_cprj,cwaveprj_tmp,cwaveprj1,ipert,isppol,&
&     my_natom,gs_hamkq%natom,nspinor,occ_k(iband),option_rhoij,pawrhoij1,usetimerev,wtk_k, &
&     comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
     call pawcprj_free(cwaveprj_tmp)
     ABI_DATATYPE_DEALLOCATE(cwaveprj_tmp)
   end if

 end if

 DBG_EXIT("COLL")

end subroutine dfpt_accrho
!!***

end module m_dfpt_mkrho
!!***
