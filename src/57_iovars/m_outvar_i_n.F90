!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_outvar_i_n
!! NAME
!!  m_outvar_i_n
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA, XG, GMR, MM)
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

module m_outvar_i_n

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_abicore
 use m_errors
 use m_results_out
 use m_errors

 use m_parser,  only : prttagm, prttagm_images

 implicit none

 private
!!***

 public :: outvar_i_n
!!***

contains
!!***

!!****f* ABINIT/outvar_i_n
!! NAME
!! outvar_i_n
!!
!! FUNCTION
!! Echo variables between acell and natom (by alphabetic order) for the ABINIT code.
!!
!! INPUTS
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables
!!  iout=unit number for echoed output
!!  jdtset_(0:ndtset_alloc)=actual index of the dataset (equal to dtsets(:)%jdtset)
!!  marr=maximum number of numbers in an array (might need to be increased ... !)
!!  multivals= <type ab_dimensions>  either 0 or 1 , depending whether the
!!     dimension has different values for different datasets
!!  mxvals= <type ab_dimensions>
!!     maximum size of some arrays along all datasets, including
!!         lpawu      =maximal value of input lpawu for all the datasets
!!         gw_nqlwl   =maximal value of input gw_nqlwl for all the datasets
!!         mband      =maximum number of bands
!!         natom      =maximal value of input natom for all the datasets
!!         natpawu    =maximal value of number of atoms on which +U is applied for all the datasets
!!         natsph     =maximal value of input natsph for all the datasets
!!         natvshift  =maximal value of input natvshift for all the datasets
!!         nconeq     =maximal value of input nconeq for all the datasets
!!         nimage     =maximal value of input nimage for all the datasets
!!         nkptgw     =maximal value of input nkptgw for all the datasets
!!         nkpthf     =maximal value of input nkpthf for all the datasets
!!         nkpt       =maximal value of input nkpt for all the datasets
!!         nnos       =maximal value of input nnos for all the datasets
!!         nqptdm     =maximal value of input nqptdm for all the datasets
!!         nspinor    =maximal value of input nspinor for all the datasets
!!         nsppol     =maximal value of input nsppol for all the datasets
!!         nsym       =maximum number of symmetries
!!         ntypat     =maximum number of type of atoms
!!         nzchempot  =maximal value of input nzchempot for all the datasets
!!  ndtset=number of datasets
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!      one data set. Use for most dimensioned arrays.
!!  npsp=number of pseudopotentials
!!  nqptdm=the number of q vectors provided by the user to calculate DM in GW
!!  prtvol_glob= if 0, minimal output volume, if 1, no restriction.
!!  response_(0:ndtset_alloc)= 1 if response variables must be output, 0 otherwise,
!!   for different datasets
!!  results_out(0:ndtset_alloc)=<type results_out_type>contains the results
!!   needed for outvars, including evolving variables
!!
!! OUTPUT
!!
!! NOTES
!! Note that this routine is called only by the processor me==0 .
!! In consequence, no use of message and wrtout routine.
!! The lines of code needed to output the defaults are preserved
!! (see last section of the routine, but are presently disabled)
!!
!! PARENTS
!!      outvars
!!
!! CHILDREN
!!      prttagm,prttagm_images
!!
!! SOURCE

subroutine outvar_i_n (dtsets,iout,&
& jdtset_,marr,multivals,mxvals,ncid,ndtset,ndtset_alloc,npsp,prtvol_glob,&
& response_,results_out,strimg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,marr,ndtset
 integer,intent(in) :: ndtset_alloc,prtvol_glob,ncid,npsp
!arrays
 integer,intent(in) :: jdtset_(0:ndtset_alloc),response_(ndtset_alloc)
 type(ab_dimensions),intent(in) :: multivals,mxvals
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
 type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)
 character(len=8),intent(in) :: strimg(mxvals%nimage)

!Local variables-------------------------------
!scalars
 integer,parameter :: nkpt_max=50
 integer :: allowed,allowed_sum,iatom,idtset,ii,iimage,ikpt,kptopt,narr
 integer :: multival,multi_natfix,multi_natfixx,multi_natfixy,multi_natfixz
 integer :: multi_atsph,multi_occopt
 integer :: natfix,natfixx,natfixy,natfixz,natom
 integer :: ndtset_kptopt,nimage,nqpt,nkpt_eff
 integer :: ntypalch,ntypat,size1,size2,tnkpt
 real(dp) :: kpoint
 character(len=1) :: firstchar_gpu
!arrays
 integer,allocatable :: iatfixio_(:,:),iatfixx_(:,:),iatfixy_(:,:)
 integer,allocatable :: iatfixz_(:,:),intarr(:,:),istwfk_2(:,:)
 integer,allocatable :: jdtset_kptopt(:),natfix_(:),natfixx_(:),natfixy_(:)
 integer,allocatable :: natfixz_(:)
 integer,allocatable :: narrm(:)
 integer,allocatable :: nimagem(:),prtimg(:,:)
 real(dp),allocatable :: dprarr(:,:),dprarr_images(:,:,:)

! *************************************************************************

!###########################################################
!### 01. Initial allocations and initialisations.

 DBG_ENTER("COLL")

 ABI_ALLOCATE(dprarr,(marr,0:ndtset_alloc))
 ABI_ALLOCATE(dprarr_images,(marr,mxvals%nimage,0:ndtset_alloc))
 ABI_ALLOCATE(intarr,(marr,0:ndtset_alloc))
 ABI_ALLOCATE(narrm,(0:ndtset_alloc))
 ABI_ALLOCATE(nimagem,(0:ndtset_alloc))
 ABI_ALLOCATE(prtimg,(mxvals%nimage,0:ndtset_alloc))

 do idtset=0,ndtset_alloc
   nimagem(idtset)=dtsets(idtset)%nimage
 end do

 firstchar_gpu=' '; if (maxval(dtsets(1:ndtset_alloc)%use_gpu_cuda)>0) firstchar_gpu='-'

!if(multivals%natom==0)natom=dtsets(1)%natom
 natom=dtsets(1)%natom
!if(multivals%nimage==0)nimage=dtsets(1)%nimage
 nimage=dtsets(1)%nimage
!if(multivals%ntypalch==0)ntypalch=dtsets(1)%ntypalch
 ntypalch=dtsets(1)%ntypalch
!if(multivals%ntypat==0)ntypat=dtsets(1)%ntypat
 ntypat=dtsets(1)%ntypat

!###########################################################
!### 02. Specific treatment for partially fixed atoms. Also compute multi_occopt for nband

!Must treat separately the translation of iatfix from the internal
!representation to the input/output representation
 ABI_ALLOCATE(natfix_,(0:ndtset_alloc))
 ABI_ALLOCATE(iatfixio_,(mxvals%natom,0:ndtset_alloc))
 ABI_ALLOCATE(natfixx_,(0:ndtset_alloc))
 ABI_ALLOCATE(iatfixx_,(mxvals%natom,0:ndtset_alloc))
 ABI_ALLOCATE(natfixy_,(0:ndtset_alloc))
 ABI_ALLOCATE(iatfixy_,(mxvals%natom,0:ndtset_alloc))
 ABI_ALLOCATE(natfixz_,(0:ndtset_alloc))
 ABI_ALLOCATE(iatfixz_,(mxvals%natom,0:ndtset_alloc))
 natfix_(0:ndtset_alloc)=0 ; iatfixio_(:,0:ndtset_alloc)=0
 natfixx_(0:ndtset_alloc)=0 ; iatfixx_(:,0:ndtset_alloc)=0
 natfixy_(0:ndtset_alloc)=0 ; iatfixy_(:,0:ndtset_alloc)=0
 natfixz_(0:ndtset_alloc)=0 ; iatfixz_(:,0:ndtset_alloc)=0
 natfixx=-1; natfixy=-1; natfixz=-1;nimage=-1;
 do idtset=1,ndtset_alloc
!  DEBUG
!  write(std_out,*)' outvar_i_n : iatfix_ for idtset= ',idtset
!  ENDDEBUG
   do iatom=1,dtsets(idtset)%natom
!    First look whether the atom is fixed along the three directions
     if( dtsets(idtset)%iatfix(1,iatom)+ &
&     dtsets(idtset)%iatfix(2,iatom)+ &
&     dtsets(idtset)%iatfix(3,iatom)   ==3 )then
       natfix_(idtset)=natfix_(idtset)+1
!      DEBUG
!      write(std_out,*)' outvar_i_n: iatom,natfix_(idtset)=',iatom,natfix_(idtset)
!      ENDDEBUG
       iatfixio_(natfix_(idtset),idtset)=iatom
     else
!      Now examine each direction, one at a time
       if( dtsets(idtset)%iatfix(1,iatom) ==1)then
         natfixx_(idtset)=natfixx_(idtset)+1
         iatfixx_(natfixx_(idtset),idtset)=iatom
       end if
       if( dtsets(idtset)%iatfix(2,iatom) ==1)then
         natfixy_(idtset)=natfixy_(idtset)+1
         iatfixy_(natfixy_(idtset),idtset)=iatom
       end if
       if( dtsets(idtset)%iatfix(3,iatom) ==1)then
         natfixz_(idtset)=natfixz_(idtset)+1
         iatfixz_(natfixz_(idtset),idtset)=iatom
       end if
     end if
   end do
!  DEBUG
!  write(std_out,*)' natfix ...'
!  write(std_out,*)natfix_(idtset),natfixx_(idtset),natfixy_(idtset),natfixz_(idtset)
!  ENDDEBUG
 end do

 multi_natfix=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(natfix_(1)/=natfix_(idtset))multi_natfix=1
   end do
 end if
 if(multi_natfix==0)natfix=natfix_(1)

 multi_natfixx=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(natfixx_(1)/=natfixx_(idtset))multi_natfixx=1
   end do
 end if
 if(multi_natfixx==0)natfixx=natfixx_(1)

 multi_natfixy=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(natfixy_(1)/=natfixy_(idtset))multi_natfixy=1
   end do
 end if
 if(multi_natfixy==0)natfixy=natfixy_(1)

 multi_natfixz=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(natfixz_(1)/=natfixz_(idtset))multi_natfixz=1
   end do
 end if
 if(multi_natfixz==0)natfixz=natfixz_(1)


 multi_occopt=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%occopt/=dtsets(idtset)%occopt)multi_occopt=1
   end do
 end if

!write(ab_out,*)' outvar_i_n : I '
!call flush(ab_out)
!###########################################################
!### 03. Print all the input variables (I)
!##

!iatfix
 narr=natfix                    ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=natfix_(idtset)
   if(idtset==0)narrm(idtset)=mxvals%natom
   intarr(1:narrm(idtset),idtset)=iatfixio_(1:narrm(idtset),idtset)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'iatfix','INT',multi_natfix)

!iatfixx
 narr=natfixx                   ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=natfixx_(idtset)
   if(idtset==0)narrm(idtset)=mxvals%natom
   intarr(1:narrm(idtset),idtset)=iatfixx_(1:narrm(idtset),idtset)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'iatfixx','INT',multi_natfixx)

!iatfixy
 narr=natfixy                   ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=natfixy_(idtset)
   if(idtset==0)narrm(idtset)=mxvals%natom
   intarr(1:narrm(idtset),idtset)=iatfixy_(1:narrm(idtset),idtset)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'iatfixy','INT',multi_natfixy)

!iatfixz
 narr=natfixz                   ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=natfixz_(idtset)
   if(idtset==0)narrm(idtset)=mxvals%natom
   intarr(1:narrm(idtset),idtset)=iatfixz_(1:narrm(idtset),idtset)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'iatfixz','INT',multi_natfixz)

!iatsph
 multi_atsph=1
 narr=dtsets(1)%natsph          ! default size for all datasets
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=dtsets(idtset)%natsph
   if(idtset==0)narrm(idtset)=mxvals%natsph
!  Need to be printed only if there is some occurence of prtdos==3 or pawfatbnd
   if (narrm(idtset)>0) then
     intarr(1:narrm(idtset),idtset)=dtsets(idtset)%iatsph(1:narrm(idtset))
   end if
   if(dtsets(idtset)%prtdos==3.or.dtsets(idtset)%pawfatbnd>0)then
     narrm(idtset)=dtsets(idtset)%natsph
   else
     narrm(idtset)=0
   end if
 end do
 if (ndtset_alloc==1.and.sum(narrm(1:ndtset_alloc))==1) multi_atsph=0

 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'iatsph','INT',multi_atsph) ! Emulating the case of multiple narr


 intarr(1,:)=dtsets(:)%iboxcut
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'iboxcut','INT',0)

 intarr(1,:)=dtsets(:)%icoulomb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'icoulomb','INT',0)

 intarr(1,:)=dtsets(:)%icutcoul
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'icutcoul','INT',0)

 intarr(1,:)=dtsets(:)%ieig2rf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ieig2rf','INT',0)

 intarr(1,:)=dtsets(:)%imgmov
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'imgmov','INT',0)

 intarr(1,:)=dtsets(:)%imgwfstor
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'imgwfstor','INT',0)

 intarr(1,:)=dtsets(:)%inclvkb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'inclvkb','INT',0)

 intarr(1,:)=dtsets(:)%intxc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'intxc','INT',0)

 intarr(1,:)=dtsets(:)%ionmov
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ionmov','INT',0)

 intarr(1,:)=dtsets(:)%iprcel
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'iprcel','INT',0)

 intarr(1,:)=dtsets(:)%iprcfc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'iprcfc','INT',0)

 intarr(1,:)=dtsets(:)%irandom
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irandom','INT',0)

 intarr(1,:)=dtsets(:)%irdbscoup
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdbscoup','INT',0)

 intarr(1,:)=dtsets(:)%irdbseig
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdbseig','INT',0)

 intarr(1,:)=dtsets(:)%irdbsreso
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdbsreso','INT',0)

 intarr(1,:)=dtsets(:)%irdddk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdddk','INT',0)

 intarr(1,:)=dtsets(:)%irdden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdden','INT',0)

 intarr(1,:)=dtsets(:)%irddvdb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irddvdb','INT',0)

 intarr(1,:)=dtsets(:)%irdefmas
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdefmas','INT',0)

 intarr(1,:)=dtsets(:)%irdhaydock
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdhaydock','INT',0)

 intarr(1,:)=dtsets(:)%irdpawden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdpawden','INT',0)

 intarr(1,:)=dtsets(:)%irdqps
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdqps','INT',0)

 intarr(1,:)=dtsets(:)%irdscr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdscr','INT',0)

 intarr(1,:)=dtsets(:)%irdsuscep
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdsuscep','INT',0)

 intarr(1,:)=dtsets(:)%irdvdw
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdvdw','INT',0)

 intarr(1,:)=dtsets(:)%irdwfk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdwfk','INT',0)

 intarr(1,:)=dtsets(:)%irdwfkfine
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdwfkfine','INT',0)

 intarr(1,:)=dtsets(:)%irdwfq
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdwfq','INT',0)

 intarr(1,:)=dtsets(:)%ird1den
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ird1den','INT',0)

 intarr(1,:)=dtsets(:)%ird1wf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ird1wf','INT',0)

 intarr(1,:)=dtsets(:)%iscf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'iscf','INT',0)

 intarr(1,:)=dtsets(:)%isecur
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'isecur','INT',0)

 intarr(1,:)=dtsets(:)%istatimg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'istatimg','INT',0)

 intarr(1,:)=dtsets(:)%istatr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'istatr','INT',0)

 intarr(1,:)=dtsets(:)%istatshft
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'istatshft','INT',0)

 if (allocated(dtsets(0)%istwfk)) then
   ! istwfk (must first restore the default istwf=0 for non-allowed k points)
   ABI_ALLOCATE(istwfk_2,(mxvals%nkpt,0:ndtset_alloc))
   istwfk_2=0;allowed_sum=0
   do idtset=1,ndtset_alloc
     nqpt=dtsets(idtset)%nqpt
     do ikpt=1,dtsets(idtset)%nkpt
       allowed=1
       do ii=1,3
         kpoint=dtsets(idtset)%kpt(ii,ikpt)/dtsets(idtset)%kptnrm
         if(nqpt/=0.and.response_(idtset)==0)kpoint=kpoint+dtsets(idtset)%qptn(ii)
         if(abs(kpoint)>1.d-10.and.abs(kpoint-0.5_dp)>1.e-10_dp )allowed=0
       end do
       allowed_sum=allowed_sum+allowed
       if(allowed==1)istwfk_2(ikpt,idtset)=dtsets(idtset)%istwfk(ikpt)
     end do
   end do

   !istwfk
   tnkpt=0
   intarr(1:marr,0:ndtset_alloc)=0 ! default value
   do idtset=1,ndtset_alloc
     nkpt_eff=dtsets(idtset)%nkpt
     if(prtvol_glob==0 .and. nkpt_eff>nkpt_max)then
       nkpt_eff=nkpt_max
       tnkpt=1
     end if
     if((multivals%nkpt/=0).and.(sum(istwfk_2(1:nkpt_eff,idtset))==0)) nkpt_eff=0
     narrm(idtset)=nkpt_eff
     intarr(1:narrm(idtset),idtset)=istwfk_2(1:narrm(idtset),idtset)
   end do

   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,nkpt_eff,narrm,ncid,ndtset_alloc,'istwfk','INT',multivals%nkpt)

   if(tnkpt==1 .and. sum(istwfk_2(1:nkpt_eff,1:ndtset_alloc))/=0 ) &
     write(iout,'(23x,a,i3,a)' ) 'outvar_i_n : Printing only first ',nkpt_max,' k-points.'
   ABI_DEALLOCATE(istwfk_2)
 end if

!ixc
 intarr(1,:)=dtsets(:)%ixc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ixc','INT',0)

!ixcpositron
 intarr(1,:)=dtsets(:)%ixcpositron
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ixcpositron','INT',0)

!ixcrot
 intarr(1,:)=dtsets(:)%ixcrot
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ixcrot','INT',0)

!ixc_sigma
 intarr(1,:)=dtsets(:)%ixc_sigma
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ixc_sigma','INT',0)

!write(ab_out,*)' outvar_i_n : J '
!call flush(ab_out)
!###########################################################
!### 03. Print all the input variables (J)
!##

 if (ndtset > 0) write(iout,"(1x,a16,1x,(t22,10i5))") 'jdtset',jdtset_(1:ndtset)

 intarr(1,:)=dtsets(:)%jellslab
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'jellslab','INT',0)

 intarr(1,:)=dtsets(:)%jfielddir(1)
 intarr(2,:)=dtsets(:)%jfielddir(2)
 intarr(3,:)=dtsets(:)%jfielddir(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'jfielddir','INT',0)

!jpawu
 prtimg(:,:)=1
 do idtset=0,ndtset_alloc
   narrm(idtset)=dtsets(idtset)%ntypat
   if (idtset==0) narrm(idtset)=mxvals%ntypat
   do iimage=1,nimagem(idtset)
     if (narrm(idtset)>0) then
       dprarr_images(1:narrm(idtset),iimage,idtset)=dtsets(idtset)%jpawu(1:narrm(idtset),iimage)
     end if
   end do
 end do
 call prttagm_images(dprarr_images,iout,jdtset_,1,marr,narrm,&
& ncid,ndtset_alloc,'jpawu','ENE',mxvals%nimage,nimagem,ndtset,prtimg,strimg)


!write(ab_out,*)' outvar_i_n : K '
!call flush(ab_out)
!###########################################################
!### 03. Print all the input variables (K)
!##

!kberry
 narr=3*dtsets(1)%nberry ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   if(idtset/=0)then
     narrm(idtset)=3*dtsets(idtset)%nberry
     if (narrm(idtset)>0)&
&     intarr(1:narrm(idtset),idtset)= reshape(dtsets(idtset)%kberry(1:3,1:dtsets(idtset)%nberry), [narrm(idtset)] )
   else
     narrm(idtset)=3*mxvals%nberry
     if (narrm(idtset)>0)&
&     intarr(1:narrm(idtset),idtset)=reshape(dtsets(idtset)%kberry(1:3,1:mxvals%nberry), [narrm(idtset)] )
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'kberry','INT',multivals%nberry)

 ! kpt
 if (allocated(dtsets(0)%kpt)) then
   tnkpt=0
   dprarr(:,0)=0
   narr=3*dtsets(1)%nkpt            ! default size for all datasets
   if(prtvol_glob==0 .and. narr>3*nkpt_max)then
     narr=3*nkpt_max
     tnkpt=1
   end if

   do idtset=1,ndtset_alloc       ! specific size for each dataset
     narrm(idtset)=3*dtsets(idtset)%nkpt
     if (narrm(idtset)>0) then
       dprarr(1:narrm(idtset),idtset)=reshape(dtsets(idtset)%kpt(1:3,1:dtsets(idtset)%nkpt), [narrm(idtset)] )
     end if

     if(prtvol_glob==0 .and. narrm(idtset)>3*nkpt_max)then
       narrm(idtset)=3*nkpt_max
       tnkpt=1
     end if

   end do
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr, narrm,ncid,ndtset_alloc,'kpt','DPR',multivals%nkpt)
   if(tnkpt==1) write(iout,'(23x,a,i3,a)' ) 'outvar_i_n : Printing only first ',nkpt_max,' k-points.'
 end if

 !intarr(1,:)=dtsets(:)%kptbounds
 !call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'kptbounds','INT',0)

!kptgw
 narr=3*dtsets(1)%nkptgw ! default size for all datasets
 dprarr(:,0)=zero
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   if(idtset/=0)then
     narrm(idtset)=3*dtsets(idtset)%nkptgw
     if (narrm(idtset)>0)&
&     dprarr(1:narrm(idtset),idtset) = reshape(dtsets(idtset)%kptgw(1:3,1:dtsets(idtset)%nkptgw), [narrm(idtset)])
   else
     narrm(idtset)=mxvals%nkptgw
     if (narrm(idtset)>0)&
&     dprarr(1:narrm(idtset),idtset) = reshape(dtsets(idtset)%kptgw(1:3,1:mxvals%nkptgw), [narrm(idtset)] )
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'kptgw','DPR',multivals%nkptgw)


 ! kptns_hf
 if (sum(dtsets(1:ndtset_alloc)%usefock) /=0 .and. allocated(dtsets(0)%kptns_hf)) then
   tnkpt=0
   dprarr(:,0)=0
   do idtset=1,ndtset_alloc       ! specific size for each dataset
     if(dtsets(idtset)%usefock/=0)then
       narrm(idtset)=3*dtsets(idtset)%nkpthf
       narr=narrm(idtset)
       if (narrm(idtset)>0) then
         dprarr(1:narrm(idtset),idtset)=reshape(dtsets(idtset)%kptns_hf(1:3,1:dtsets(idtset)%nkpthf), [narrm(idtset)] )
       end if
     else
       narrm(idtset)=0
     end if
     if(prtvol_glob==0 .and. narrm(idtset)>3*nkpt_max)then
       narrm(idtset)=3*nkpt_max
       narr=narrm(idtset)
       tnkpt=1
     end if
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'kptns_hf','DPR',multivals%nkpthf)
   if(tnkpt==1) write(iout,'(23x,a,i3,a)' ) 'outvar_i_n : Printing only first ',nkpt_max,' k-points.'
 end if

 dprarr(1,:)=dtsets(:)%kptnrm
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'kptnrm','DPR',0)

 intarr(1,:)=dtsets(:)%kptopt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'kptopt','INT',0)

!kptrlatt
 if(sum((dtsets(1:ndtset_alloc)%kptopt)**2)/=0)then
   ndtset_kptopt=0
   intarr(1:9,0)=reshape( dtsets(0)%kptrlatt, [9] )
   ABI_ALLOCATE(jdtset_kptopt,(0:ndtset_alloc))
!  Define the set of datasets for which kptopt>0
   do idtset=1,ndtset_alloc
     kptopt=dtsets(idtset)%kptopt
     if(kptopt>0)then
       ndtset_kptopt=ndtset_kptopt+1
       jdtset_kptopt(ndtset_kptopt)=jdtset_(idtset)
       intarr(1:9,ndtset_kptopt)=reshape( dtsets(idtset)%kptrlatt , [9] )
     end if
   end do
   if(ndtset_kptopt>0)then
     call prttagm(dprarr,intarr,iout,jdtset_kptopt,6,marr,9,narrm,ncid,ndtset_kptopt,'kptrlatt','INT',0)
   end if
   ABI_DEALLOCATE(jdtset_kptopt)
 end if

 dprarr(1,:)=dtsets(:)%kptrlen
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'kptrlen','DPR',0)

 intarr(1,:)=dtsets(:)%kssform
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'kssform','INT',0)

!write(ab_out,*)' outvar_i_n : L '
!call flush(ab_out)
!###########################################################
!### 03. Print all the input variables (L)
!##

!lexexch
 narr=mxvals%ntypat             ! default size for all datasets
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=dtsets(idtset)%ntypat
   if(idtset==0)narrm(idtset)=mxvals%ntypat
   if (narrm(idtset)>0) intarr(1:narrm(idtset),idtset)=dtsets(idtset)%lexexch(1:narrm(idtset))
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,narr, narrm,ncid,ndtset_alloc,'lexexch','INT',multivals%ntypat)

!ldaminushalf
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=dtsets(idtset)%ntypat
   if(idtset==0)narrm(idtset)=mxvals%ntypat
   if (narrm(idtset)>0) then
     intarr(1:narrm(idtset),idtset)=dtsets(idtset)%ldaminushalf(1:narrm(idtset))
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,narr, narrm,ncid,ndtset_alloc,'ldaminushalf','INT',multivals%ntypat)

!localrdwf
 intarr(1,:)=dtsets(:)%localrdwf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'localrdwf','INT',0)

!lpawu
 narr=mxvals%ntypat                    ! default size for all datasets
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=dtsets(idtset)%ntypat
   if(idtset==0)narrm(idtset)=mxvals%ntypat
   if (narrm(idtset)>0) then
     intarr(1:narrm(idtset),idtset)=dtsets(idtset)%lpawu(1:narrm(idtset))
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,narr,narrm,ncid,ndtset_alloc,'lpawu','INT',multivals%ntypat)

#if defined HAVE_LOTF
 if (any(dtsets(:)%ionmov==23)) then
   intarr(1,:)=dtsets(:)%lotf_classic
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'lotf_classic','INT',0)
   intarr(1,:)=dtsets(:)%lotf_nitex
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'lotf_nitex','INT',0)
   intarr(1,:)=dtsets(:)%lotf_nneigx
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'lotf_nneigx','INT',0)
   intarr(1,:)=dtsets(:)%lotf_version
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'lotf_version','INT',0)
 end if
#endif


!write(ab_out,*)' outvar_i_n : M '
!call flush(ab_out)
!###########################################################
!### 03. Print all the input variables (M)
!##

 intarr(1,:)=dtsets(:)%macro_uj
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'macro_uj','INT',0)

 dprarr(1,:)=dtsets(:)%magcon_lambda
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'magcon_lambda','DPR',0)

 intarr(1,:)=dtsets(:)%magconon
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'magconon','INT',0)

 dprarr(1,:)=dtsets(:)%maxestep
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'maxestep','ENE',0)

 intarr(1,:)=dtsets(:)%max_ncpus
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'max_ncpus','INT',0)

 intarr(1,:)=dtsets(:)%maxnsym
 call prttagm(dprarr,intarr,iout,jdtset_,4,marr,1,narrm,ncid,ndtset_alloc,'maxnsym','INT',0)

 dprarr(1,:)=dtsets(:)%mbpt_sciss
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'mbpt_sciss','ENE',0)

 dprarr(1,:)=dtsets(:)%mdf_epsinf
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'mdf_epsinf','DPR',0)

 dprarr(1,:)=dtsets(:)%mdtemp(1)
 dprarr(2,:)=dtsets(:)%mdtemp(2)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,2,narrm,ncid,ndtset_alloc,'mdtemp','DPR',0)

 dprarr(1,:)=dtsets(:)%mdwall
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'mdwall','LEN',0)

 intarr(1,:)=dtsets(:)%mem_test
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'mem_test','INT',0)

 dprarr(1,:)=dtsets(:)%mep_mxstep
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'mep_mxstep','LEN',0)

 intarr(1,:)=dtsets(:)%mep_solver
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'mep_solver','INT',0)

 intarr(1,:)=dtsets(:)%mffmem
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'mffmem','INT',0)

!mixalch
 prtimg(:,:)=1
 do idtset=0,ndtset_alloc
   if(idtset/=0)then
     size1=dtsets(idtset)%npspalch ; size2=dtsets(idtset)%ntypalch
   else
     size1=npsp ; size2=mxvals%ntypat
   end if
   narrm(idtset)=size1*size2
   do iimage=1,nimagem(idtset)
     if (narrm(idtset)>0) then
       dprarr_images(1:narrm(idtset),iimage,idtset)=&
&       reshape(results_out(idtset)%mixalch(1:size1,1:size2,iimage), [narrm(idtset)] )
     end if
   end do
 end do
 call prttagm_images(dprarr_images,iout,jdtset_,1,marr,narrm,ncid,ndtset_alloc,'mixalch','DPR',&
& mxvals%nimage,nimagem,ndtset,prtimg,strimg)

 intarr(1,:)=dtsets(:)%mixprec
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'mixprec','INT',0)

!mixesimgf
 dprarr(1:marr,0)=zero              ! default value
 narr=mxvals%nimage                 ! default size for all datasets
 if(any(abs(dtsets(1:ndtset_alloc)%imgmov-6)==0))then
   multival=multivals%nimage
   do idtset=1,ndtset_alloc           ! specific size and array for each dataset
     narrm(idtset)=dtsets(idtset)%nimage
     dprarr(1:narrm(idtset),idtset)=dtsets(idtset)%mixesimgf(1:narrm(idtset))
   end do
 else
   multival=0
   narrm(1:ndtset_alloc)=narr
   dprarr(1:marr,1:ndtset_alloc)=zero
 endif
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'mixesimgf','DPR',multival)

 intarr(1,:)=dtsets(:)%mkmem
 call prttagm(dprarr,intarr,iout,jdtset_,5,marr,1,narrm,ncid,ndtset_alloc,'mkmem','INT',0)

 if(sum(response_(:))/=0)then
   intarr(1,:)=dtsets(:)%mkqmem
   call prttagm(dprarr,intarr,iout,jdtset_,5,marr,1,narrm,ncid,ndtset_alloc,'mkqmem','INT',0)
 end if

 if(sum(response_(:))/=0)then
   intarr(1,:)=dtsets(:)%mk1mem
   call prttagm(dprarr,intarr,iout,jdtset_,5,marr,1,narrm,ncid,ndtset_alloc,'mk1mem','INT',0)
 end if

 intarr(1,:)=dtsets(:)%mqgrid
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'mqgrid','INT',0)

 intarr(1,:)=dtsets(:)%mqgriddg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'mqgriddg','INT',0)

!###########################################################
!### 03. Print all the input variables (N)
!##

 intarr(1,:)=natfix_(:)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'natfix','INT',0)

 intarr(1,:)=natfixx_(:)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'natfixx','INT',0)

 intarr(1,:)=natfixy_(:)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'natfixy','INT',0)

 intarr(1,:)=natfixz_(:)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'natfixz','INT',0)

 intarr(1,:)=dtsets(0:ndtset_alloc)%natom
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'natom','INT',0,forceprint=2)

!natsph
!Need to be printed only if there is some occurence of prtdos==3 or
!pawfatbnd>0
 narr=1                      ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=1
   intarr(1,idtset)=dtsets(idtset)%natsph

   if(dtsets(idtset)%prtdos==3.or.dtsets(idtset)%pawfatbnd>0)then
     narrm(idtset)=1
   else
     narrm(idtset)=0
   end if
 end do
 if (ndtset_alloc==1.and.sum(narrm(1:ndtset_alloc))==1) multi_atsph=0
 ! Emulating multiple size for narrm
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'natsph','INT',multi_atsph)

!natsph_extra
 intarr(1,:)=dtsets(0:ndtset_alloc)%natsph_extra
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'natsph_extra','INT',0)

 if(dtsets(1)%occopt==2)then
   narr=dtsets(1)%nkpt*dtsets(1)%nsppol                      ! default size for all datasets
 else
   narr=1
 end if
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   if(dtsets(idtset)%occopt==2)then
     narrm(idtset)=dtsets(idtset)%nkpt*dtsets(idtset)%nsppol
   else
     narrm(idtset)=1
   end if

   if (narrm(idtset)>0) then
     intarr(1:narrm(idtset),idtset)=dtsets(idtset)%nband(1:narrm(idtset))
   end if
 end do

 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,narr,&
& narrm,ncid,ndtset_alloc,'nband','INT',multivals%nkpt+multivals%nsppol+multi_occopt)

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%natvshift
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'natvshift','INT',0)

 if(sum(dtsets(1:ndtset_alloc)%usefock)/=0)then
   intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%nbandhf
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nbandhf','INT',0)
 end if

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%nbandkss
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nbandkss','INT',0)

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%nbdblock
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nbdblock','INT',0)

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%nbdbuf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nbdbuf','INT',0)

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%nberry
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nberry','INT',0)

 intarr(1,:)=dtsets(:)%nc_xccc_gspace
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nc_xccc_gspace','INT',0) !, firstchar="-")

 intarr(1,:)=dtsets(:)%nconeq
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nconeq','INT',0)

 intarr(1,:)=dtsets(:)%nctime
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nctime','INT',0)

!ndtset
 if(ndtset>0)then
   intarr(1,:)=ndtset
   intarr(1,0)=0
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ndtset','INT',0)
 end if

 intarr(1,:)=dtsets(:)%ndivsm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ndivsm','INT',0)

 !intarr(1,:)=dtsets(:)%nkpath
 !call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nkpath','INT',0)

 intarr(1,:)=dtsets(:)%ndynimage
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ndynimage','INT',0)

 intarr(1,:)=dtsets(:)%neb_algo
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'neb_algo','INT',0)

 dprarr(1,:)=dtsets(:)%neb_spring(1)
 dprarr(2,:)=dtsets(:)%neb_spring(2)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,2,narrm,ncid,ndtset_alloc,'neb_spring','DPR',0)

 intarr(1,:)=dtsets(:)%nfreqim
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nfreqim','INT',0)

 intarr(1,:)=dtsets(:)%nfreqre
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nfreqre','INT',0)

 intarr(1,:)=dtsets(:)%nfreqsp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nfreqsp','INT',0)

 intarr(1,:)=dtsets(:)%ngfft(1)
 intarr(2,:)=dtsets(:)%ngfft(2)
 intarr(3,:)=dtsets(:)%ngfft(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'ngfft','INT',0)

 intarr(1,:)=dtsets(:)%ngfftdg(1)
 intarr(2,:)=dtsets(:)%ngfftdg(2)
 intarr(3,:)=dtsets(:)%ngfftdg(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'ngfftdg','INT',0)

 intarr(1,:)=dtsets(:)%nimage
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nimage','INT',0)

 intarr(1,:)=dtsets(:)%nkpt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nkpt','INT',0)

 intarr(1,:)=dtsets(:)%nkptgw
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nkptgw','INT',0)

 if(sum(dtsets(1:ndtset_alloc)%usefock)/=0)then
   intarr(1,:)=dtsets(:)%nkpthf
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nkpthf','INT',0)
 end if

 intarr(1,:)=dtsets(:)%nline
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nline','INT',0)

 intarr(1,:)=dtsets(:)%nloalg(1)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nloc_alg','INT',0)

 intarr(1,:)=dtsets(:)%nloalg(2)*(dtsets(:)%nloalg(3)+1)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nloc_mem','INT',0)

 intarr(1,:)=dtsets(:)%nnos
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nnos','INT',0)

 intarr(1,:)=dtsets(:)%nnsclo
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nnsclo','INT',0)

 intarr(1,:)=dtsets(:)%nnsclohf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nnsclohf','INT',0)

 intarr(1,:)=dtsets(:)%nomegasf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nomegasf','INT',0)

 intarr(1,:)=dtsets(:)%nomegasi
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nomegasi','INT',0)

 intarr(1,:)=dtsets(:)%nomegasrd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nomegasrd','INT',0)

 intarr(1,:)=dtsets(:)%nonlinear_info
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nonlinear_info','INT',0)

 dprarr(1,:)=dtsets(:)%noseinert
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'noseinert','DPR',0)

 intarr(1,:)=dtsets(:)%npband
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npband','INT',0,firstchar="-")

 intarr(1,:)=dtsets(:)%npfft
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npfft','INT',0,firstchar="-")

 intarr(1,:)=dtsets(:)%nphf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nphf','INT',0,firstchar="-")

 intarr(1,:)=dtsets(:)%npimage
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npimage','INT',0,firstchar="-")

 intarr(1,:)=dtsets(:)%npkpt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npkpt','INT',0,firstchar='-')

 intarr(1,:)=dtsets(:)%nppert
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nppert','INT',0,firstchar="-")

 if(multivals%ntypat/=0 .or. (multivals%ntypat==0 .and. ntypat/=npsp) )then
   intarr(1,:)=dtsets(:)%npsp
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npsp','INT',0)
 end if

 intarr(1,:)=dtsets(:)%npspinor
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npspinor','INT',0,firstchar="-")

 intarr(1,:)=dtsets(0:ndtset_alloc)%npulayit
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npulayit','INT',0)

 intarr(1,:)=dtsets(:)%npweps
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npweps','INT',0)

 intarr(1,:)=dtsets(:)%npwkss
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npwkss','INT',0)

 intarr(1,:)=dtsets(:)%npwsigx
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npwsigx','INT',0)

 intarr(1,:)=dtsets(:)%npwwfn
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npwwfn','INT',0)

 intarr(1,:)=dtsets(:)%np_slk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'np_slk','INT',0,firstchar="-")

 intarr(1,:)=dtsets(:)%nqpt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nqpt','INT',0)

 intarr(1,:)=dtsets(:)%nqptdm
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'nqptdm','INT',0)

 intarr(1,:)=dtsets(:)%npvel
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npvel','INT',0)

 intarr(1,:)=dtsets(:)%nscforder
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nscforder','INT',0)

!nshiftk
 if(sum((dtsets(1:ndtset_alloc)%kptopt)**2)/=0)then
   ndtset_kptopt=0
   intarr(1:1,0)=dtsets(0)%nshiftk
   ABI_ALLOCATE(jdtset_kptopt,(0:ndtset_alloc))
!  Define the set of datasets for which kptopt>0
   do idtset=1,ndtset_alloc
     kptopt=dtsets(idtset)%kptopt
     if(kptopt>0)then
       ndtset_kptopt=ndtset_kptopt+1
       jdtset_kptopt(ndtset_kptopt)=jdtset_(idtset)
       intarr(1:1,ndtset_kptopt)=dtsets(idtset)%nshiftk
     end if
   end do
   if(ndtset_kptopt>0)then
     call prttagm(dprarr,intarr,iout,jdtset_kptopt,2,marr,1,narrm,ncid,ndtset_kptopt,'nshiftk','INT',0)
   end if
   ABI_DEALLOCATE(jdtset_kptopt)
 end if

 intarr(1,:)=dtsets(:)%nspden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nspden','INT',0)

 intarr(1,:)=dtsets(:)%nspinor
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nspinor','INT',0)

 intarr(1,:)=dtsets(:)%nsppol
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nsppol','INT',0)

 intarr(1,:)=dtsets(:)%nstep
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nstep','INT',0)

 intarr(1,:)=dtsets(:)%nsym
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nsym','INT',0)

 intarr(1,:)=dtsets(:)%ntime
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ntime','INT',0)

 intarr(1,:)=dtsets(:)%ntimimage
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ntimimage','INT',0)

 intarr(1,:)=dtsets(:)%ntypalch
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ntypalch','INT',0)

 intarr(1,:)=dtsets(:)%ntypat
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ntypat','INT',0,forceprint=2)

!nucdipmom
 dprarr(:,0)=0.0_dp
 narr=3*natom ! default size for all datasets
 do idtset=1,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=3*dtsets(idtset)%natom
   if (narrm(idtset)>0) then
     dprarr(1:narrm(idtset),idtset)=&
&     reshape(dtsets(idtset)%nucdipmom(1:3,1:dtsets(idtset)%natom), [narrm(idtset)])
   end if
   if(sum(abs( dtsets(idtset)%nucdipmom(1:3,1:dtsets(idtset)%natom))) < tol12 ) narrm(idtset)=0
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,narr,narrm,ncid,ndtset_alloc,'nucdipmom','DPR',multivals%natom)

 intarr(1,:)=dtsets(:)%nwfshist
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nwfshist','INT',0)

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%nzchempot
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'nzchempot','INT',0)

!###########################################################
!## Deallocation for generic arrays, and for i-n variables

 ABI_DEALLOCATE(dprarr)
 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(narrm)
 ABI_DEALLOCATE(nimagem)
 ABI_DEALLOCATE(dprarr_images)
 ABI_DEALLOCATE(prtimg)

 ABI_DEALLOCATE(natfix_)
 ABI_DEALLOCATE(iatfixio_)
 ABI_DEALLOCATE(natfixx_)
 ABI_DEALLOCATE(iatfixx_)
 ABI_DEALLOCATE(natfixy_)
 ABI_DEALLOCATE(iatfixy_)
 ABI_DEALLOCATE(natfixz_)
 ABI_DEALLOCATE(iatfixz_)

 DBG_EXIT("COLL")

end subroutine outvar_i_n
!!***

end module m_outvar_i_n
!!***
