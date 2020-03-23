!!****m* ABINIT/m_outvar_o_z
!! NAME
!!  m_outvar_o_z
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, MM)
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

module m_outvar_o_z

 use defs_basis

 use m_errors
 use m_results_out
 use m_abicore
 use m_xmpi
 use m_dtset

 use m_geometry,     only : mkrdim, xred2xcart
 use m_parser,       only : prttagm, prttagm_images, ab_dimensions

 implicit none

 private
!!***

 public :: outvar_o_z
!!***

contains
!!***

!!****f* ABINIT/outvar_o_z
!! NAME
!! outvar_o_z
!!
!! FUNCTION
!! Echo variables between acell and gw_ ... (by alphabetic order) for the ABINIT code.
!!
!! INPUTS
!!  choice= 1 if echo of preprocessed variables, 2 if echo after call driver
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
!!         nkpt       =maximal value of input nkpt for all the datasets
!!         nnos       =maximal value of input nnos for all the datasets
!!         nqptdm     =maximal value of input nqptdm for all the datasets
!!         nspinor    =maximal value of input nspinor for all the datasets
!!         nsppol     =maximal value of input nsppol for all the datasets
!!         nsym       =maximum number of symmetries
!!         ntypat     =maximum number of type of atoms
!!         nzchempot  =maximal value of input nzchempot for all the datasets
!!  ncid= NetCDF handler
!!  ndtset=number of datasets
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!      one data set. Use for most dimensioned arrays.
!!  npsp=number of pseudopotentials
!!  prtvol_glob= if 0, minimal output volume, if 1, no restriction.
!!  results_out(0:ndtset_alloc)=<type results_out_type>contains the results
!!   needed for outvars, including evolving variables
!!  timopt=input variable to modulate the timing
!!
!! OUTPUT
!!
!! NOTES
!! Note that this routine is called only by the processor me==0 .
!! In consequence, no use of message and wrtout routine.
!! The lines of code needed to output the defaults are preserved
!! (see last section of the routine, but are presently disabled)
!!
!!  Note that acell, occ, rprim, xred and vel might have been modified by the
!!  computation, so that their values if choice=1 or choice=2 will differ.
!!
!! PARENTS
!!      outvars
!!
!! CHILDREN
!!      mkrdim,prtocc,prttagm,prttagm_images,xred2xcart
!!
!! SOURCE

 subroutine outvar_o_z(choice,dtsets,iout,&
& jdtset_,marr,multivals,mxvals,ncid,ndtset,ndtset_alloc,npsp,prtvol_glob,&
& results_out,strimg,timopt)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,iout,marr,ndtset
 integer,intent(in) :: ndtset_alloc,prtvol_glob,ncid,npsp,timopt
!arrays
 integer,intent(in) :: jdtset_(0:ndtset_alloc)
 type(ab_dimensions),intent(in) :: multivals,mxvals
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
 type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)
 character(len=8),intent(in) :: strimg(mxvals%nimage)

!Local variables-------------------------------
!scalars
 integer,parameter :: nkpt_max=50
 integer :: iat,icount,idtset,ii,iimage,ndtset_alloc_tmp
 integer :: narr
 integer :: multi_kptopt
 integer :: natom
 integer :: nimage,nnos,nsym
 integer :: ntypalch,ntypat,size1,size2,tnkpt,timopt_default,tmpimg0
 logical :: compute_static_images
 character(len=1) :: firstchar_gpu
!arrays
 integer,allocatable :: narrm(:)
 integer,allocatable :: nimagem(:),prtimg(:,:)
 integer,allocatable :: intarr(:,:)
 real(dp) :: rprimd(3,3)
 real(dp),allocatable :: dprarr(:,:),dprarr_images(:,:,:)
 real(dp),allocatable :: xangst(:,:),xcart(:,:),xred(:,:)
 real(dp),allocatable :: xangst_(:,:,:,:),xcart_(:,:,:,:)

! *************************************************************************

!###########################################################
!### 01. Initial allocations and initialisations.

 ABI_ALLOCATE(dprarr,(marr,0:ndtset_alloc))
 ABI_ALLOCATE(dprarr_images,(marr,mxvals%nimage,0:ndtset_alloc))
 ABI_ALLOCATE(intarr,(marr,0:ndtset_alloc))
 ABI_ALLOCATE(narrm,(0:ndtset_alloc))
 ABI_ALLOCATE(nimagem,(0:ndtset_alloc))
 ABI_ALLOCATE(prtimg,(mxvals%nimage,0:ndtset_alloc))

 do idtset=0,ndtset_alloc
   nimagem(idtset)=dtsets(idtset)%nimage
 end do

 firstchar_gpu=' ';if (maxval(dtsets(1:ndtset_alloc)%use_gpu_cuda)>0) firstchar_gpu='-'

 natom=dtsets(1)%natom
 nimage=dtsets(1)%nimage
 nnos=dtsets(1)%nnos
 nsym=-1;nsym=dtsets(1)%nsym
 ntypalch=dtsets(1)%ntypalch
 ntypat=dtsets(1)%ntypat

!###########################################################
!### 02. Specific treatment for occopt, xangst, xcart, xred

!Must compute xangst and xcart
 ABI_ALLOCATE(xangst_,(3,mxvals%natom,mxvals%nimage,0:ndtset_alloc))
 ABI_ALLOCATE(xcart_,(3,mxvals%natom,mxvals%nimage,0:ndtset_alloc))
 xangst_(:,:,:,:)=0.0_dp ; xcart_(:,:,:,:)=0.0_dp

 do idtset=1,ndtset_alloc
   natom=dtsets(idtset)%natom
   ABI_ALLOCATE(xred,(3,natom))
   ABI_ALLOCATE(xangst,(3,natom))
   ABI_ALLOCATE(xcart,(3,natom))
   do iimage=1,dtsets(idtset)%nimage
     xred(:,1:natom)=results_out(idtset)%xred(:,1:natom,iimage)
     call mkrdim(results_out(idtset)%acell(:,iimage),results_out(idtset)%rprim(:,:,iimage),rprimd)
!    Compute xcart from xred and rprimd
     call xred2xcart(natom,rprimd,xcart,xred)
!    Compute xangst from xcart
     xangst(:,:)=xcart(:,:)*Bohr_Ang
!    Save the data
     xangst_(1:3,1:natom,iimage,idtset)=xangst(:,:)
     xcart_(1:3,1:natom,iimage,idtset)=xcart(:,:)
   end do
   if(dtsets(idtset)%nimage/=mxvals%nimage)then
     xangst_(1:3,1:natom,dtsets(idtset)%nimage+1:mxvals%nimage,idtset)=zero
     xcart_(1:3,1:natom,dtsets(idtset)%nimage+1:mxvals%nimage,idtset)=zero
   end if
   ABI_DEALLOCATE(xred)
   ABI_DEALLOCATE(xangst)
   ABI_DEALLOCATE(xcart)
 end do

!###########################################################
!### 03. Print all the input variables (O)
!##

!occ
!The use of prttagm for occ if occopt>=2 is not possible because
!the different k-point and spins must be separated on different lines ...
!Also, if prtvol=-1 and NC file has been created, only print one dataset
 ndtset_alloc_tmp=ndtset_alloc
 if(ncid<0)ndtset_alloc_tmp=1
 call prtocc(dtsets,iout,jdtset_,mxvals,ndtset_alloc_tmp,nimagem,prtvol_glob,results_out,strimg)

 intarr(1,:)=dtsets(:)%occopt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'occopt','INT',0)

 dprarr(1,:)=dtsets(:)%omegasimax
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'omegasimax','ENE',0)

 dprarr(1,:)=dtsets(:)%omegasrdmax
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'omegasrdmax','ENE',0)

 intarr(1,:)=dtsets(:)%optcell
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'optcell','INT',0)

 intarr(1,:)=dtsets(:)%optdriver
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'optdriver','INT',0)

 intarr(1,:)=dtsets(:)%optforces
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'optforces','INT',0)

 intarr(1,:)=dtsets(:)%optnlxccc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'optnlxccc','INT',0)

 intarr(1,:)=dtsets(:)%optstress
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'optstress','INT',0)

 intarr(1,:)=dtsets(:)%orbmag
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'orbmag','INT',0)

 intarr(1,:)=dtsets(:)%ortalg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ortalg','INT',0,firstchar=firstchar_gpu)

!###########################################################
!### 03. Print all the input variables (P)
!##

 intarr(1,:)=dtsets(:)%paral_atom
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'paral_atom','INT',0, firstchar="-")
 !call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'paral_atom','INT',0)

 intarr(1,:)=dtsets(:)%paral_kgb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'paral_kgb','INT',0)

 intarr(1,:)=dtsets(:)%paral_rf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'paral_rf','INT',0)

 intarr(1,:)=dtsets(:)%pawcpxocc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawcpxocc','INT',0)

 intarr(1,:)=dtsets(:)%pawcross
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawcross','INT',0)

 dprarr(1,:)=dtsets(:)%pawecutdg
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'pawecutdg','ENE',0)

 intarr(1,:)=dtsets(:)%pawfatbnd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawfatbnd','INT',0)

 intarr(1,:)=dtsets(:)%pawlcutd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawlcutd','INT',0)

 intarr(1,:)=dtsets(:)%pawlmix
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawlmix','INT',0)

 intarr(1,:)=dtsets(:)%pawmixdg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawmixdg','INT',0)

 intarr(1,:)=dtsets(:)%pawnhatxc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawnhatxc','INT',0)

 intarr(1,:)=dtsets(:)%pawnphi
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawnphi','INT',0)

 intarr(1,:)=dtsets(:)%pawntheta
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawntheta','INT',0)

 intarr(1,:)=dtsets(:)%pawnzlm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawnzlm','INT',0)

 intarr(1,:)=dtsets(:)%pawoptmix
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawoptmix','INT',0)

 intarr(1,:)=dtsets(:)%pawoptosc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawoptosc','INT',0)

 dprarr(1,:)=dtsets(:)%pawovlp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawovlp','DPR',0)

 intarr(1,:)=dtsets(:)%pawprtdos
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawprtdos','INT',0)

 intarr(1,:)=dtsets(:)%pawprtvol
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawprtvol','INT',0)

 intarr(1,:)=dtsets(:)%pawprtwf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawprtwf','INT',0)

 intarr(1,:)=dtsets(:)%pawprt_b
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawprt_b','INT',0)

 intarr(1,:)=dtsets(:)%pawprt_k
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawprt_k','INT',0)

 intarr(1,:)=dtsets(:)%pawspnorb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawspnorb','INT',0)

 intarr(1,:)=dtsets(:)%pawstgylm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawstgylm','INT',0)

 intarr(1,:)=dtsets(:)%pawsushat
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawsushat','INT',0)

 intarr(1,:)=dtsets(:)%pawujat
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawujat','INT',0)

 dprarr(1,:)=dtsets(:)%pawujrad
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawujrad','LEN',0)

 dprarr(1,:)=dtsets(:)%pawujv
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawujv','ENE',0)

 intarr(1,:)=dtsets(:)%pawusecp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawusecp','INT',0)

 intarr(1,:)=dtsets(:)%pawxcdev
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawxcdev','INT',0)

 intarr(1,:)=dtsets(:)%ph_intmeth
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'ph_intmeth','INT',0)

 intarr(1,:)=dtsets(:)%ph_ndivsm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ph_ndivsm','INT',0)

 do idtset=0,ndtset_alloc
   intarr(1:3,idtset)=dtsets(idtset)%ph_ngqpt
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'ph_ngqpt','INT',0)

 intarr(1,:)=dtsets(:)%ph_nqpath
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ph_nqpath','INT',0)

 intarr(1,:)=dtsets(:)%ph_nqshift
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'ph_nqshift','INT',0)

 dprarr(1,:)=dtsets(:)%ph_smear
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'ph_smear','ENE',0)

 dprarr(1,:)=dtsets(:)%ph_wstep
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'ph_wstep','ENE',0)

!pimass
 icount=0
 do idtset=0, ndtset_alloc
   do ii = 1, ntypat
     dprarr(ii,idtset) = dtsets(idtset)%pimass(ii)
     if (dtsets(idtset)%pimass(ii)/=dtsets(idtset)%amu_orig(ii,1)) icount=1
   end do ! end loop over ntypat
 end do ! end loop over datasets
 if (icount/=0) then
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,narrm,ncid,ndtset_alloc,'pimass','DPR',0)
 end if

 intarr(1,:)=dtsets(:)%pimd_constraint
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pimd_constraint','INT',0)

 intarr(1,:)=dtsets(:)%pitransform
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pitransform','INT',0)

 intarr(1,:)=dtsets(:)%plowan_bandi
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'plowan_bandi','INT',0)

 intarr(1,:)=dtsets(:)%plowan_bandf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'plowan_bandf','INT',0)

 intarr(1,:)=dtsets(:)%plowan_compute
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'plowan_compute','INT',0)

 intarr(1,:)=dtsets(:)%plowan_natom
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'plowan_natom','INT',0)

 intarr(1,:)=dtsets(:)%plowan_nt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'plowan_nt','INT',0)

 intarr(1,:)=dtsets(:)%plowan_realspace
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'plowan_realspace','INT',0)


!plowan_it
 narr=100
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=3*dtsets(idtset)%plowan_nt
   if(idtset==0)narrm(idtset)=100
   if (narrm(idtset)>0.and.dtsets(idtset)%plowan_compute>0) then
     intarr(1:narrm(idtset),idtset)=dtsets(idtset)%plowan_it(1:narrm(idtset))
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'plowan_it','INT',1)

!plowan_iatom
 narr=mxvals%natom
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=dtsets(idtset)%plowan_natom
   if(idtset==0)narrm(idtset)=mxvals%natom
   if (narrm(idtset)>0.and.dtsets(idtset)%plowan_compute>0) then
     intarr(1:narrm(idtset),idtset)=dtsets(idtset)%plowan_iatom(1:narrm(idtset))
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'plowan_iatom','INT',1)

!plowan_nbl
 narr=mxvals%natom
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=dtsets(idtset)%plowan_natom
   if(idtset==0)narrm(idtset)=mxvals%natom
   if (narrm(idtset)>0.and.dtsets(idtset)%plowan_compute>0) then
     intarr(1:narrm(idtset),idtset)=dtsets(idtset)%plowan_nbl(1:narrm(idtset))
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'plowan_nbl','INT',1)

!plowan_lcalc
 narr=12*mxvals%natom
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=sum(dtsets(idtset)%plowan_nbl(1:dtsets(idtset)%plowan_natom))
   if(idtset==0)narrm(idtset)=12*mxvals%natom
   if (narrm(idtset)>0.and.dtsets(idtset)%plowan_compute>0) then
     intarr(1:narrm(idtset),idtset)=dtsets(idtset)%plowan_lcalc(1:narrm(idtset))
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'plowan_lcalc','INT',1)

!plowan_projcalc
 narr=12*mxvals%natom
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=sum(dtsets(idtset)%plowan_nbl(1:dtsets(idtset)%plowan_natom))
   if(idtset==0)narrm(idtset)=12*mxvals%natom
   if (narrm(idtset)>0.and.dtsets(idtset)%plowan_compute>0) then
     intarr(1:narrm(idtset),idtset)=dtsets(idtset)%plowan_projcalc(1:narrm(idtset))
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'plowan_projcalc','INT',1)

 dprarr(1,:)=dtsets(:)%polcen(1)
 dprarr(2,:)=dtsets(:)%polcen(2)
 dprarr(3,:)=dtsets(:)%polcen(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'polcen','DPR',0)

 intarr(1,:)=dtsets(:)%posdoppler
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'posdoppler','INT',0)

 intarr(1,:)=dtsets(:)%positron
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'positron','INT',0)

 intarr(1,:)=dtsets(:)%posnstep
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'posnstep','INT',0)

 dprarr(1,:)=dtsets(:)%posocc
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'posocc','DPR',0)

 dprarr(1,:)=dtsets(:)%postoldfe
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'postoldfe','ENE',0)

 dprarr(1,:)=dtsets(:)%postoldff
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'postoldff','DPR',0)

 dprarr(1,:)=dtsets(:)%ppmfrq
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'ppmfrq','ENE',0)

 intarr(1,:)=dtsets(:)%ppmodel
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ppmodel','INT',0)


 intarr(1,:)=dtsets(:)%prepanl
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prepanl','INT',0)

 intarr(1,:)=dtsets(:)%prepgkk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prepgkk','INT',0)

!prtatlist
 if(multivals%natom==0)then
   do idtset=0,ndtset_alloc
     intarr(1:natom,idtset)=dtsets(idtset)%prtatlist(1:natom)
   end do
   intarr(1:mxvals%natom,0)=(/ (ii,ii=1,mxvals%natom) /)
   call prttagm(dprarr,intarr,iout,jdtset_,4,marr,natom,narrm,ncid,ndtset_alloc,'prtatlist','INT',0)
 else
!  This thing will disapear with new generalized prttagm
 end if

 intarr(1,:)=dtsets(:)%prtbbb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtbbb','INT',0)

 intarr(1,:)=dtsets(:)%prtbltztrp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtbltztrp','INT',0)

 intarr(1,:)=dtsets(:)%prtcif
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtcif','INT',0)

 intarr(1,:)=dtsets(:)%prtden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtden','INT',0)

 intarr(1,:)=dtsets(:)%prtdensph
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtdensph','INT',0)

 intarr(1,:)=dtsets(:)%prtdos
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtdos','INT',0)

 intarr(1,:)=dtsets(:)%prtdosm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtdosm','INT',0)

 intarr(1,:)=dtsets(:)%prtebands
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtebands','INT',0)

 intarr(1,:)=dtsets(:)%prtefg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtefg','INT',0)

 intarr(1,:)=dtsets(:)%prtefmas
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtefmas','INT',0)

 intarr(1,:)=dtsets(:)%prteig
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prteig','INT',0)

 intarr(1,:)=dtsets(:)%prtelf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtelf','INT',0)

 intarr(1,:)=dtsets(:)%prteliash
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prteliash','INT',0)

 intarr(1,:)=dtsets(:)%prtfc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtfc','INT',0)

 intarr(1,:)=dtsets(:)%prtfull1wf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtfull1wf','INT',0)

 intarr(1,:)=dtsets(:)%prtfsurf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtfsurf','INT',0)

 intarr(1,:)=dtsets(:)%prtgden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtgden','INT',0)

 intarr(1,:)=dtsets(:)%prtgeo
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtgeo','INT',0)

 intarr(1,:)=dtsets(:)%prtgkk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtgkk','INT',0)

 intarr(1,:)=dtsets(:)%prtgsr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtgsr','INT',0)

 intarr(1,:)=dtsets(:)%prtkden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtkden','INT',0)

 intarr(1,:)=dtsets(:)%prtlden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtlden','INT',0)

 intarr(1,:)=dtsets(:)%prtnabla
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtnabla','INT',0)

 intarr(1,:)=dtsets(:)%prtphbands
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtphbands','INT',0)

 intarr(1,:)=dtsets(:)%prtphdos
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtphdos','INT',0)

 intarr(1,:)=dtsets(:)%prtphsurf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtphsurf','INT',0)

 intarr(1,:)=dtsets(:)%prtposcar
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtposcar','INT',0)

 intarr(1,:)=dtsets(:)%prtprocar
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtprocar','INT',0)

 intarr(1,:)=dtsets(:)%prtpot
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtpot','INT',0)

 intarr(1,:)=dtsets(:)%prtpsps
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtpsps','INT',0)

 intarr(1,:)=dtsets(:)%prtspcur
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtspcur','INT',0)

 intarr(1,:)=dtsets(:)%prtstm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtstm','INT',0)

 intarr(1,:)=dtsets(:)%prtsuscep
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtsuscep','INT',0)

 intarr(1,:)=dtsets(:)%prtvclmb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtvclmb','INT',0)

 intarr(1,:)=dtsets(:)%prtvha
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtvha','INT',0)

 intarr(1,:)=dtsets(:)%prtvhxc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtvhxc','INT',0)

 intarr(1,:)=dtsets(:)%prtkbff
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtkbff','INT',0)

 intarr(1,:)=dtsets(:)%prtvol
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtvol','INT',0)

 intarr(1,:)=dtsets(:)%prtvolimg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtvolimg','INT',0)

 intarr(1,:)=dtsets(:)%prtvpsp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtvpsp','INT',0)

 intarr(1,:)=dtsets(:)%prtvxc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtvxc','INT',0)

 intarr(1,:)=dtsets(:)%prtwant
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtwant','INT',0)

 intarr(1,:)=dtsets(:)%prtwf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtwf','INT',0)

 intarr(1,:)=dtsets(:)%prtwf_full
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtwf_full','INT',0)

 intarr(1,:)=dtsets(:)%prtxml
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtxml','INT',0)

!prt1dm
 intarr(1,:)=dtsets(:)%prt1dm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prt1dm','INT',0)

 do idtset=0, ndtset_alloc
   do ii = 1, ntypat
     dprarr(ii,idtset) = dtsets(idtset)%ptcharge(ii)
   end do ! end loop over ntypat
 end do ! end loop over datasets
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,narrm,ncid,ndtset_alloc,'ptcharge','DPR',0)

 intarr(1,:)=dtsets(:)%ptgroupma
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ptgroupma','INT',0)

 dprarr(1,:)=dtsets(:)%pvelmax(1)
 dprarr(2,:)=dtsets(:)%pvelmax(2)
 dprarr(3,:)=dtsets(:)%pvelmax(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'pvelmax','DPR',0)

 dprarr(1,:)=dtsets(:)%pw_unbal_thresh
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'pw_unbal_thresh','DPR',0)

!###########################################################
!### 03. Print all the input variables (Q)
!##

!qmass
 narr=nnos ! default size for all datasets
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=dtsets(idtset)%nnos
   if(idtset==0)narrm(idtset)=mxvals%nnos
   if (narrm(idtset)>0) then
     dprarr(1:narrm(idtset),idtset)=dtsets(idtset)%qmass(1:narrm(idtset))
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'qmass','DPR',multivals%nnos)

 intarr(1,:)=dtsets(:)%qprtrb(1)
 intarr(2,:)=dtsets(:)%qprtrb(2)
 intarr(3,:)=dtsets(:)%qprtrb(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'qprtrb','INT',0)

 dprarr(1,:)=dtsets(:)%qptn(1)
 dprarr(2,:)=dtsets(:)%qptn(2)
 dprarr(3,:)=dtsets(:)%qptn(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'qpt','DPR',0)

!qptdm
 narr=3*dtsets(1)%nqptdm ! default size for all datasets
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   if(idtset/=0)then
     narrm(idtset)=3*dtsets(idtset)%nqptdm
     if (narrm(idtset)>0)&
&     dprarr(1:narrm(idtset),idtset)=&
&     reshape(dtsets(idtset)%qptdm(1:3,&
&     1:dtsets(idtset)%nqptdm),&
&     (/ narrm(idtset) /) )
   else
     narrm(idtset)=3*mxvals%nqptdm
     if (narrm(idtset)>0)&
&     dprarr(1:narrm(idtset),idtset)=&
&     reshape(dtsets(idtset)%qptdm(1:3,&
&     1:mxvals%nqptdm),&
&     (/ narrm(idtset) /) )
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'qptdm','DPR',multivals%nqptdm)

 do idtset=0, ndtset_alloc
   do ii = 1, ntypat
     dprarr(ii,idtset) = dtsets(idtset)%quadmom(ii)
   end do ! end loop over ntypat
 end do ! end loop over datasets
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,narrm,ncid,ndtset_alloc,'quadmom','DPR',0)

!###########################################################
!### 03. Print all the input variables (R)
!##

!variables used for the random positions in unit cell
 intarr(1,:)=dtsets(:)%random_atpos
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'random_atpos','INT',0)

 dprarr(1,:)=dtsets(:)%ratsm
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'ratsm','LEN',0)

 do idtset=0, ndtset_alloc
   do ii = 1, ntypat
     dprarr(ii,idtset) = dtsets(idtset)%ratsph(ii)
   end do ! end loop over ntypat
 end do ! end loop over datasets
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,narrm,ncid,ndtset_alloc,'ratsph','LEN',0)

 dprarr = zero
 dprarr(1,:) = dtsets(:)%ratsph_extra
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'ratsph_extra','LEN',0)

 dprarr(1,:)=dtsets(:)%rcut
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'rcut','LEN',0)

!Variables used for recursion method
 dprarr(1,:)=dtsets(:)%recefermi
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'recefermi','ENE',0)

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%recgratio
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'recgratio','INT',0)

 intarr(1,:)=dtsets(:)%recnpath
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'recnpath','INT',0)

 intarr(1,:)=dtsets(:)%recnrec
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'recnrec','INT',0)

 intarr(1,:)=dtsets(:)%recptrott
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'recptrott','INT',0)

 dprarr(1,:)=dtsets(:)%recrcut
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'recrcut','LEN',0)

 intarr(1,:)=dtsets(:)%rectesteg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rectesteg','INT',0)

 dprarr(1,:)=dtsets(:)%rectolden
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'rectolden','DPR',0)

 dprarr(1,:)=dtsets(:)%red_dfield(1)    !!HONG
 dprarr(2,:)=dtsets(:)%red_dfield(2)
 dprarr(3,:)=dtsets(:)%red_dfield(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'red_dfield','DPR',0)

 dprarr(1,:)=dtsets(:)%red_efield(1)    !!HONG
 dprarr(2,:)=dtsets(:)%red_efield(2)
 dprarr(3,:)=dtsets(:)%red_efield(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'red_efield','DPR',0)

 dprarr(1,:)=dtsets(:)%red_efieldbar(1)   !!HONG
 dprarr(2,:)=dtsets(:)%red_efieldbar(2)
 dprarr(3,:)=dtsets(:)%red_efieldbar(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'red_efieldbar','DPR',0)

 intarr(1,:)=dtsets(:)%restartxf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'restartxf','INT',0)

 intarr(1,:)=dtsets(:)%rfasr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rfasr','INT',0)

 intarr(1,:)=dtsets(:)%rfatpol(1)
 intarr(2,:)=dtsets(:)%rfatpol(2)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,2,narrm,ncid,ndtset_alloc,'rfatpol','INT',0)

 intarr(1,:)=dtsets(:)%rfddk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rfddk','INT',0)

 intarr(1,:)=dtsets(:)%rfdir(1)
 intarr(2,:)=dtsets(:)%rfdir(2)
 intarr(3,:)=dtsets(:)%rfdir(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'rfdir','INT',0)

 intarr(1,:)=dtsets(:)%rfelfd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rfelfd','INT',0)

 intarr(1,:)=dtsets(:)%rfmagn
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rfmagn','INT',0)

 intarr(1,:)=dtsets(:)%rfmeth
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rfmeth','INT',0)

 intarr(1,:)=dtsets(:)%rfphon
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rfphon','INT',0)

 intarr(1,:)=dtsets(:)%rfstrs
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rfstrs','INT',0)

 intarr(1,:)=dtsets(:)%rfuser
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rfuser','INT',0)

 intarr(1,:)=dtsets(:)%rf2_dkdk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rf2_dkdk','INT',0)

 intarr(1,:)=dtsets(:)%rf2_dkde
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rf2_dkde','INT',0)

 intarr(1,:)=dtsets(:)%rf2_pert1_dir(1)
 intarr(2,:)=dtsets(:)%rf2_pert1_dir(2)
 intarr(3,:)=dtsets(:)%rf2_pert1_dir(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'rf2_pert1_dir','INT',0)

 intarr(1,:)=dtsets(:)%rf2_pert2_dir(1)
 intarr(2,:)=dtsets(:)%rf2_pert2_dir(2)
 intarr(3,:)=dtsets(:)%rf2_pert2_dir(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'rf2_pert2_dir','INT',0)

 dprarr(1,:)=dtsets(:)%rhoqpmix
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'rhoqpmix','DPR',0)

!rprim
 prtimg(:,:)=1
 do idtset=0,ndtset_alloc
   narrm(idtset)=9
   do iimage=1,nimagem(idtset)
     if (narrm(idtset)>0) then
       dprarr_images(1:narrm(idtset),iimage,idtset)=&
&       reshape(results_out(idtset)%rprim(1:3,1:3,iimage), (/ narrm(idtset) /) )
     end if
   end do
 end do
 call prttagm_images(dprarr_images,iout,jdtset_,-2,marr,narrm,ncid,ndtset_alloc,'rprim','DPR',&
& mxvals%nimage,nimagem,ndtset,prtimg,strimg,forceprint=2)


!###########################################################
!### 03. Print all the input variables (S)
!##

!shiftk (printed only when kptopt>0)
 if(sum((dtsets(1:ndtset_alloc)%kptopt)**2)/=0)then
   multi_kptopt=0
   dprarr(:,0)=0.0_dp
   narr=3*dtsets(1)%nshiftk ! default size for all datasets
   do idtset=1,ndtset_alloc       ! specific size for each dataset
     narrm(idtset)=3*dtsets(idtset)%nshiftk
     if (narrm(idtset)>0) then
       dprarr(1:narrm(idtset),idtset)=&
&       reshape(dtsets(idtset)%shiftk(1:3,1:dtsets(idtset)%nshiftk),(/ narrm(idtset) /) )
     end if
     if(dtsets(idtset)%kptopt<=0)then
       narrm(idtset)=0
       multi_kptopt=1
     end if
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'shiftk','DPR',multivals%nshiftk)
!  End of test to see whether kptopt/=0 for some dataset
 end if

 intarr(1,:)=dtsets(:)%sigma_bsum_range(1)
 intarr(2,:)=dtsets(:)%sigma_bsum_range(2)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,2,narrm,ncid,ndtset_alloc,'sigma_bsum_range','INT',0)

 dprarr(1,:)=dtsets(:)%sigma_erange(1)
 dprarr(2,:)=dtsets(:)%sigma_erange(2)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,2,narrm,ncid,ndtset_alloc,'sigma_erange','ENE',0)

 intarr(1,:)=dtsets(:)%transport_ngkpt(1)
 intarr(2,:)=dtsets(:)%transport_ngkpt(2)
 intarr(3,:)=dtsets(:)%transport_ngkpt(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'transport_ngkpt','INT',0)

 intarr(1,:)=dtsets(:)%sigma_ngkpt(1)
 intarr(2,:)=dtsets(:)%sigma_ngkpt(2)
 intarr(3,:)=dtsets(:)%sigma_ngkpt(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'sigma_ngkpt','INT',0)

 intarr(1,:)=dtsets(:)%signperm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'signperm','INT',0)

 dprarr(1,:)=dtsets(:)%slabwsrad
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'slabwsrad','DPR',0)

 dprarr(1,:)=dtsets(:)%slabzbeg
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'slabzbeg','DPR',0)

 dprarr(1,:)=dtsets(:)%slabzend
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'slabzend','DPR',0)

 intarr(1,:)=dtsets(:)%slk_rankpp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'slk_rankpp','INT',0)

 intarr(1,:)=dtsets(:)%smdelta
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'smdelta','INT',0)

 do idtset=0,ndtset_alloc
   intarr(1:npsp,idtset)=dtsets(idtset)%so_psp(1:npsp)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,npsp,narrm,ncid,ndtset_alloc,'so_psp','INT',0)

 dprarr(1,:)=dtsets(:)%spbroad
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'spbroad','ENE',0)

 intarr(1,:)=dtsets(:)%spgroup
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'spgroup','INT',0)

!spinat
 dprarr(:,0)=0.0_dp
 narr=3*natom ! default size for all datasets
 do idtset=1,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=3*dtsets(idtset)%natom
   if (narrm(idtset)>0) then
     dprarr(1:narrm(idtset),idtset)=reshape(dtsets(idtset)%spinat(1:3,1:dtsets(idtset)%natom), (/narrm(idtset)/))
   end if
   if(sum(abs( dtsets(idtset)%spinat(1:3,1:dtsets(idtset)%natom))) < tol12 ) narrm(idtset)=0
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,narr,narrm,ncid,ndtset_alloc,'spinat','DPR',multivals%natom)

 dprarr(1,:)=dtsets(:)%spinmagntarget
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'spinmagntarget','DPR',0)

 intarr(1,:)=dtsets(:)%spmeth
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'spmeth','INT',0)

 dprarr(1,:)=dtsets(:)%spnorbscl
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'spnorbscl','DPR',0)

 dprarr(1,:)=dtsets(:)%stmbias
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'stmbias','DPR',0)

 dprarr(1,:)=dtsets(:)%strfact
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'strfact','DPR',0)

 intarr(1,:)=dtsets(:)%string_algo
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'string_algo','INT',0)

 do ii=1,6
   dprarr(ii,:)=dtsets(:)%strtarget(ii)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,6,narrm,ncid,ndtset_alloc,'strtarget','DPR',0)

!strten
 if(choice==2)then
   prtimg(:,:)=1
   do idtset=0,ndtset_alloc       ! specific size for each dataset
     compute_static_images=(dtsets(idtset)%istatimg>0)
     narrm(idtset)=6
     if(dtsets(idtset)%iscf>=0)then
       do iimage=1,dtsets(idtset)%nimage
         if (narrm(idtset)>0) then
           dprarr_images(1:narrm(idtset),iimage,idtset)=results_out(idtset)%strten(:,iimage)
         end if
         if(.not.(dtsets(idtset)%dynimage(iimage)==1.or.compute_static_images))then
           prtimg(iimage,idtset)=0
         end if
       end do
     else
       narrm(idtset)=0
     end if
   end do
!  This is a trick to force printing of strten even if zero, still not destroying the value of nimagem(0).
   tmpimg0=nimagem(0)
   nimagem(0)=0
   call prttagm_images(dprarr_images,iout,jdtset_,2,marr,narrm,ncid,ndtset_alloc,'strten','DPR',&
&   mxvals%nimage,nimagem,ndtset,prtimg,strimg)
   nimagem(0)=tmpimg0
 end if

!symafm
 intarr(:,0)=1
 narr=nsym ! default size for all datasets
 do idtset=1,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=dtsets(idtset)%nsym
   if (narrm(idtset)>0) then
     intarr(1:narrm(idtset),idtset)=dtsets(idtset)%symafm(1:narrm(idtset))
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'symafm','INT', multivals%nsym)

 intarr(1,:)=dtsets(:)%symchi
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'symchi','INT',0)

 intarr(1,:)=dtsets(:)%symdynmat
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'symdynmat','INT',0)

 intarr(1,:)=dtsets(:)%symmorphi
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'symmorphi','INT',0)

!symrel
 intarr(1:9,0)=(/ 1,0,0, 0,1,0, 0,0,1 /)
 narr=9*nsym ! default size for all datasets
 do idtset=1,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=9*dtsets(idtset)%nsym
   if (narrm(idtset)>0) then
     intarr(1:narrm(idtset),idtset)=&
&     reshape(dtsets(idtset)%symrel(1:3,1:3,1:dtsets(idtset)%nsym), [narrm(idtset)] )
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,3,marr,narr,narrm,ncid,ndtset_alloc,'symrel','INT', multivals%nsym)

 intarr(1,:)=dtsets(:)%symsigma
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'symsigma','INT',0)

 intarr(1,:)=dtsets(:)%symv1scf
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'symv1scf','INT',0)

!###########################################################
!### 03. Print all the input variables (T)
!##

 dprarr(1,:)=dtsets(:)%td_maxene
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'td_maxene','DPR',0)

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%td_mexcit
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'td_mexcit','INT',0)

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%tfkinfunc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'tfkinfunc','INT',0)

 dprarr(1,:)=dtsets(:)%tfw_toldfe
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tfw_toldfe','ENE',0)

 intarr(1,:)=dtsets(:)%tim1rev
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'tim1rev','INT',0)


!timopt
 timopt_default=1; if(xmpi_paral==1) timopt_default=0

 if(timopt/=timopt_default)then
   intarr(1,:)=timopt
   intarr(1,0)=timopt_default
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'timopt','INT',0)
 end if

!WVL - tails related variables
 intarr(1,:)=dtsets(:)%tl_nprccg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'tl_nprccg','INT',0)
 dprarr(1,:)=dtsets(:)%tl_radius
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tl_radius','DPR',0)

!tnons
 dprarr(:,0)=0.0_dp
 narr=3*nsym ! default size for all datasets
 do idtset=1,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=3*dtsets(idtset)%nsym
   if (narrm(idtset)>0) then
     dprarr(1:narrm(idtset),idtset)=reshape(dtsets(idtset)%tnons(1:3,1:dtsets(idtset)%nsym), [narrm(idtset)])
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,-3,marr,narr,narrm,ncid,ndtset_alloc,'tnons','DPR',multivals%nsym)

 dprarr(1,:)=dtsets(:)%toldfe
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'toldfe','ENE',0)

 dprarr(1,:)=dtsets(:)%tolmxde
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tolmxde','ENE',0)

 dprarr(1,:)=dtsets(:)%toldff
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'toldff','DPR',0)

 dprarr(1,:)=dtsets(:)%tolimg
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tolimg','ENE',0)

 dprarr(1,:)=dtsets(:)%tolmxf
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tolmxf','DPR',0)

 dprarr(1,:)=dtsets(:)%tolrde
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tolrde','DPR',0)

 dprarr(1,:)=dtsets(:)%tolrff
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tolrff','DPR',0)

 dprarr(1,:)=dtsets(:)%tolsym
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tolsym','DPR',0)

 dprarr(1,:)=dtsets(:)%tolvrs
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tolvrs','DPR',0)

 dprarr(1,:)=dtsets(:)%tolwfr
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tolwfr','DPR',0)

 dprarr(1,:)=dtsets(:)%tphysel
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tphysel','ENE',0)

 dprarr(1,:) = dtsets(:)%tmesh(1); dprarr(2,:) = dtsets(:)%tmesh(2); dprarr(3,:) = dtsets(:)%tmesh(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'tmesh','DPR',0)

 dprarr(1,:)=dtsets(:)%tsmear
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tsmear','ENE',0)

!typat
 narr=natom                      ! default size for all datasets
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=dtsets(idtset)%natom
   if(idtset==0)narrm(idtset)=mxvals%natom
   if (narrm(idtset)>0) then
     intarr(1:narrm(idtset),idtset)=dtsets(idtset)%typat(1:narrm(idtset))
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,4,marr,narr,narrm,ncid,ndtset_alloc,'typat','INT',multivals%natom,forceprint=2)

!###########################################################
!### 03. Print all the input variables (U)
!##
 intarr(1,:)=dtsets(:)%ucrpa
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ucrpa','INT',0)

 intarr(1,:)=dtsets(:)%ucrpa_bands(1)
 intarr(2,:)=dtsets(:)%ucrpa_bands(2)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,2,narrm,ncid,ndtset_alloc,'ucrpa_bands','INT',0)

 dprarr(1,:)=dtsets(:)%ucrpa_window(1)
 dprarr(2,:)=dtsets(:)%ucrpa_window(2)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,2,narrm,ncid,ndtset_alloc,'ucrpa_window','ENE',0)

!upawu
 prtimg(:,:)=1
 do idtset=0,ndtset_alloc
   narrm(idtset)=dtsets(idtset)%ntypat
   if (idtset==0) narrm(idtset)=mxvals%ntypat
   do iimage=1,nimagem(idtset)
     if (narrm(idtset)>0) then
       dprarr_images(1:narrm(idtset),iimage,idtset)=dtsets(idtset)%upawu(1:narrm(idtset),iimage)
     end if
   end do
 end do
 call prttagm_images(dprarr_images,iout,jdtset_,1,marr,narrm,&
& ncid,ndtset_alloc,'upawu','ENE',mxvals%nimage,nimagem,ndtset,prtimg,strimg)

 intarr(1,:)=dtsets(:)%usedmatpu
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'usedmatpu','INT',0)

 intarr(1,:)=dtsets(:)%usedmft
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'usedmft','INT',0)

 intarr(1,:)=dtsets(:)%useexexch
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'useexexch','INT',0)

 intarr(1,:)=dtsets(:)%usefock
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'usefock','INT',0)

 intarr(1,:)=dtsets(:)%usepotzero
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'usepotzero','INT',0)

 intarr(1,:)=dtsets(:)%usekden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'usekden','INT',0)

 intarr(1,:)=dtsets(:)%use_gemm_nonlop
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'use_gemm_nonlop','INT',0)

 intarr(1,:)=dtsets(:)%use_yaml
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'use_yaml','INT',0)

 intarr(1,:)=dtsets(:)%use_nonscf_gkk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'use_nonscf_gkk','INT',0)

 intarr(1,:)=dtsets(:)%usepawu
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'usepawu','INT',0)

 intarr(1,:)=dtsets(:)%usepead
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'usepead','INT',0)

 intarr(1,:)=dtsets(:)%useria
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'useria','INT',0)

 intarr(1,:)=dtsets(:)%userib
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'userib','INT',0)

 intarr(1,:)=dtsets(:)%useric
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'useric','INT',0)

 intarr(1,:)=dtsets(:)%userid
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'userid','INT',0)

 intarr(1,:)=dtsets(:)%userie
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'userie','INT',0)

 dprarr(1,:)=dtsets(:)%userra
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'userra','DPR',0)

 dprarr(1,:)=dtsets(:)%userrb
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'userrb','DPR',0)

 dprarr(1,:)=dtsets(:)%userrc
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'userrc','DPR',0)

 dprarr(1,:)=dtsets(:)%userrd
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'userrd','DPR',0)

 dprarr(1,:)=dtsets(:)%userre
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'userre','DPR',0)

 intarr(1,:)=dtsets(:)%usewvl
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'usewvl','INT',0)

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%usexcnhat_orig
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'usexcnhat','INT',0)

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%useylm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'useylm','INT',0,firstchar=firstchar_gpu)

 intarr(1,:)=dtsets(:)%use_gpu_cuda
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'use_gpu_cuda','INT',0,firstchar=firstchar_gpu)

 intarr(1,:)=dtsets(:)%use_slk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'use_slk','INT',0, firstchar="-")


!###########################################################
!### 03. Print all the input variables (V)
!##

 dprarr(1,:)=dtsets(:)%vcutgeo(1)
 dprarr(2,:)=dtsets(:)%vcutgeo(2)
 dprarr(3,:)=dtsets(:)%vcutgeo(3)
 call prttagm(dprarr,intarr,iout,jdtset_,3,marr,3,narrm,ncid,ndtset_alloc,'vcutgeo','DPR',0)

 if(sum(dtsets(1:ndtset_alloc)%prtwant) >1)then
!  van der Waals correction with MLWFs related variables
   if(any(dtsets(1:ndtset_alloc)%vdw_xc==10).or.any(dtsets(1:ndtset_alloc)%vdw_xc==11).or.&
&   any(dtsets(1:ndtset_alloc)%vdw_xc==14))then
     intarr(1,:)=dtsets(:)%vdw_nfrag
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'vdw_nfrag','INT',0)
   end if !vdw_xc==10,11,14
   if(any(dtsets(1:ndtset_alloc)%vdw_xc==10).or.any(dtsets(1:ndtset_alloc)%vdw_xc==11).or.&
&   any(dtsets(1:ndtset_alloc)%vdw_xc==14))then
     intarr(1,:)=dtsets(:)%vdw_supercell(1)
     intarr(2,:)=dtsets(:)%vdw_supercell(2)
     intarr(3,:)=dtsets(:)%vdw_supercell(3)
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'vdw_supercell','INT',0)
   end if !vdw_xc==10,11,14
 end if !prtwant>1

 dprarr(1,:)=dtsets(:)%vdw_tol
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'vdw_tol','DPR',0)
 dprarr(1,:)=dtsets(:)%vdw_tol_3bt
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'vdw_tol_3bt','DPR',0)

 if(sum(dtsets(1:ndtset_alloc)%prtwant) >1)then
!  van der Waals correction with MLWFs related variables
   if(any(dtsets(1:ndtset_alloc)%vdw_xc==10).or.any(dtsets(1:ndtset_alloc)%vdw_xc==11))then
     do iat=1,mxvals%natom
       intarr(iat,:)=dtsets(:)%vdw_typfrag(iat)
     end do
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,mxvals%natom,narrm,ncid,ndtset_alloc,'vdw_typfrag','INT',0)
   end if !vdw_xc==10 or xc==11
 end if !prtwant>1

 intarr(1,:)=dtsets(:)%vdw_xc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'vdw_xc','INT',0)

 if(sum(dtsets(1:ndtset_alloc)%prtvdw) >1)then
   if(any(dtsets(1:ndtset_alloc)%vdw_xc<10))then
     dprarr(1,:)=dtsets(:)%vdw_df_threshold
     call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'vdw_df_threshold','ENE',0)
   end if
 end if

!vel
 prtimg(:,:)=1
 do idtset=0,ndtset_alloc
   if(idtset/=0)then
     size1=dtsets(idtset)%natom
   else
     size1=mxvals%natom
   end if
   narrm(idtset)=3*size1
   do iimage=1,nimagem(idtset)
     if (narrm(idtset)>0) then
       dprarr_images(1:narrm(idtset),iimage,idtset)=&
&       reshape(results_out(idtset)%vel(1:3,1:size1,iimage), (/ narrm(idtset) /) )
     end if
   end do
 end do
 call prttagm_images(dprarr_images,iout,jdtset_,2,marr,narrm,ncid,ndtset_alloc,'vel','DPR',&
& mxvals%nimage,nimagem,ndtset,prtimg,strimg)

!vel_cell
!At present, vel_cell does not depend on image... but this might change in the future.
 prtimg(:,:)=1
 if (.true.) then
!  if(mxvals%nimage==1)then
   do idtset=0,ndtset_alloc
     dprarr(1:9,idtset)= reshape(results_out(idtset)%vel_cell(:,:,1),(/9/))
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,9,narrm,ncid,ndtset_alloc,'vel_cell','DPR',0)
!  else
!  do idtset=1,ndtset_alloc       ! specific size for each dataset
!  nimagem(idtset)=dtsets(idtset)%nimage
!  narrm(idtset)=9
!  do iimage=1,dtsets(idtset)%nimage
!  if (narrm(idtset)>0) then
!  dprarr_images(1:narrm(idtset),iimage,idtset)=&
!  &         reshape(results_out(idtset)%vel_cell(1:3,1:3,iimage),&
!  &         (/ narrm(idtset) /) )
!  end if
!  end do
!  end do
!  call prttagm_images(dprarr_images,iout,jdtset_,&
!  &   marr,narrm,ncid,ndtset_alloc,'vel_cell',&
!  &   mxvals%nimage,nimagem,ndtset,prtimg,strimg)
 end if

 dprarr(1,:)=dtsets(:)%vis
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'vis','DPR',0)

 dprarr(1,:)=dtsets(:)%vprtrb(1)
 dprarr(2,:)=dtsets(:)%vprtrb(2)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,2,narrm,ncid,ndtset_alloc,'vprtrb','ENE',0)


!###########################################################
!### 03. Print all the input variables (W)
!##

 dprarr(1,:)=dtsets(:)%wfmix
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'wfmix','DPR',0)

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%wfk_task
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'wfk_task','INT',0)

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%wfoptalg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'wfoptalg','INT',0,firstchar=firstchar_gpu)

!wtatcon
 narr=3*natom*dtsets(1)%nconeq ! default size for all datasets
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   if(idtset/=0)then
     narrm(idtset)=3*dtsets(idtset)%natom*dtsets(idtset)%nconeq
     if (narrm(idtset)>0)&
&     dprarr(1:narrm(idtset),idtset)=&
&     reshape(dtsets(idtset)%wtatcon(1:3,1:dtsets(idtset)%natom,1:dtsets(idtset)%nconeq),(/ narrm(idtset) /) )
   else
     narrm(idtset)=3*mxvals%natom*mxvals%nconeq
     if (narrm(idtset)>0)&
&     dprarr(1:narrm(idtset),idtset)=&
&     reshape(dtsets(idtset)%wtatcon(1:3,1:mxvals%natom,1:mxvals%nconeq),(/ narrm(idtset) /) )
   end if
 end do

 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'wtatcon','DPR',&
   multivals%natom+multivals%nconeq)

!wtk
 if (allocated(dtsets(0)%wtk)) then
   tnkpt=0
   dprarr(:,0)=1
   narr=dtsets(1)%nkpt ! default size for all datasets
   if(prtvol_glob==0 .and. narr>nkpt_max)then
     narr=nkpt_max
     tnkpt=1
   end if
   do idtset=1,ndtset_alloc       ! specific size for each dataset
     narrm(idtset)=dtsets(idtset)%nkpt
     if (narrm(idtset)>0) dprarr(1:narrm(idtset),idtset)=dtsets(idtset)%wtk(1:narrm(idtset))+tol12

     if(prtvol_glob==0 .and. narrm(idtset)>nkpt_max)then
       narrm(idtset)=nkpt_max
       tnkpt=1
     end if
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,4,marr,narr,narrm,ncid,ndtset_alloc,'wtk','DPR',multivals%nkpt)
   if(tnkpt==1) write(iout,'(23x,a,i3,a)' ) 'outvars : Printing only first ',nkpt_max,' k-points.'
 end if

!WVL - wavelets variables
 if (any(dtsets(:)%usewvl==1)) then
   intarr(1,:)=dtsets(:)%wvl_bigdft_comp
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'wvl_bigdft_comp','INT',0)
   dprarr(1,:)=dtsets(:)%wvl_crmult
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'wvl_crmult','DPR',0)
   dprarr(1,:)=dtsets(:)%wvl_frmult
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'wvl_frmult','DPR',0)
   dprarr(1,:)=dtsets(:)%wvl_hgrid
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'wvl_hgrid','DPR',0)
   intarr(1,:)=dtsets(:)%wvl_nprccg
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'wvl_nprccg','INT',0)
 end if

!Wannier90 interface related variables
 if(sum(dtsets(1:ndtset_alloc)%prtwant) >1)then
   intarr(1,:)=dtsets(:)%w90iniprj
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'w90iniprj','INT',0)
   intarr(1,:)=dtsets(:)%w90prtunk
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'w90prtunk','INT',0)
 end if !prtwant>1

!###########################################################
!### 03. Print all the input variables (X)
!##

!xangst
 prtimg(:,:)=1
 do idtset=0,ndtset_alloc
   if(idtset/=0)then
     size1=dtsets(idtset)%natom
   else
     size1=mxvals%natom
   end if
   narrm(idtset)=3*size1
   do iimage=1,nimagem(idtset)
     if (narrm(idtset)>0) then
       dprarr_images(1:narrm(idtset),iimage,idtset)=reshape(xangst_(1:3,1:size1,iimage,idtset), (/narrm(idtset)/))
     end if
   end do
 end do

 call prttagm_images(dprarr_images,iout,jdtset_,-2,marr,narrm,ncid,ndtset_alloc,'xangst','DPR',&
& mxvals%nimage,nimagem,ndtset,prtimg,strimg)

!xcart
 prtimg(:,:)=1
 do idtset=0,ndtset_alloc
   if(idtset/=0)then
     size1=dtsets(idtset)%natom
   else
     size1=mxvals%natom
   end if
   narrm(idtset)=3*size1
   do iimage=1,nimagem(idtset)
     if (narrm(idtset)>0) then
       dprarr_images(1:narrm(idtset),iimage,idtset)=reshape(xcart_(1:3,1:size1,iimage,idtset), (/ narrm(idtset) /) )
     end if
   end do
 end do

 call prttagm_images(dprarr_images,iout,jdtset_,-2,marr,narrm,ncid,ndtset_alloc,'xcart','DPR',&
& mxvals%nimage,nimagem,ndtset,prtimg,strimg)

 dprarr(1,:)=dtsets(:)%xc_denpos
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'xc_denpos','DPR',0)

 dprarr(1,:)=dtsets(:)%xc_tb09_c
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'xc_tb09_c','DPR',0)

!xred
 prtimg(:,:)=1
 do idtset=0,ndtset_alloc
   if(idtset/=0)then
     size2=dtsets(idtset)%natom
   else
     size2=mxvals%natom
   end if
   narrm(idtset)=3*size2
   do iimage=1,nimagem(idtset)
     if (narrm(idtset)>0) then
       dprarr_images(1:narrm(idtset),iimage,idtset)=&
&       reshape(results_out(idtset)%xred(:,1:size2,iimage), (/ narrm(idtset) /) )
     end if
   end do
 end do
 call prttagm_images(dprarr_images,iout,jdtset_,-2,marr,narrm,ncid,ndtset_alloc,'xred','DPR',&
& mxvals%nimage,nimagem,ndtset,prtimg,strimg,forceprint=2)

!xredsph_extra
 do idtset=0,ndtset_alloc
   if(idtset/=0)then
     size2=dtsets(idtset)%natsph_extra
   else
     size2=0
   end if
   narrm(idtset)=3*size2
   if (narrm(idtset)>0) then
     dprarr(1:narrm(idtset),idtset)= reshape(dtsets(idtset)%xredsph_extra(:,1:size2), (/ narrm(idtset) /) )
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'xredsph_extra','DPR',1)


!###########################################################
!### 03. Print all the input variables (Y)
!##

!###########################################################
!### 03. Print all the input variables (Z)
!##

 dprarr(1,:)=dtsets(:)%zcut
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'zcut','ENE',0)

!zeemanfield
 dprarr(1,:)=dtsets(:)%zeemanfield(1)
 dprarr(2,:)=dtsets(:)%zeemanfield(2)
 dprarr(3,:)=dtsets(:)%zeemanfield(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'zeemanfield','BFI',0)

!ziontypat   ! After all, should always echo this value
 if(sum(dtsets(:)%ntypalch)>0)then
   narr=ntypat                    ! default size for all datasets
   do idtset=0,ndtset_alloc       ! specific size for each dataset
     narrm(idtset)=dtsets(idtset)%ntypat
     if(idtset==0)narrm(idtset)=mxvals%ntypat
     if (narrm(idtset)>0) then
       dprarr(1:narrm(idtset),idtset)=dtsets(idtset)%ziontypat(1:narrm(idtset))
     end if
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
&   narrm,ncid,ndtset_alloc,'ziontypat','DPR',multivals%ntypat,forceprint=2)
 end if

 do idtset=0,ndtset_alloc
   dprarr(1:npsp,idtset)=dtsets(idtset)%znucl(1:npsp)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,4,marr,npsp,narrm,ncid,ndtset_alloc,'znucl','DPR',0,forceprint=2)

!###########################################################
!## Deallocation for generic arrays, and for n-z variables

 ABI_DEALLOCATE(dprarr)
 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(narrm)
 ABI_DEALLOCATE(nimagem)
 ABI_DEALLOCATE(dprarr_images)
 ABI_DEALLOCATE(prtimg)
 ABI_DEALLOCATE(xangst_)
 ABI_DEALLOCATE(xcart_)

contains
!!***

!!****f* ABINIT/prtocc
!!
!! NAME
!! prtocc
!!
!! FUNCTION
!! Print the content of occ.
!! Due to the need to distinguish between different k-points and
!! different spin polarisations, prttagm.f cannot be used.
!! So, need a dedicated routine.
!!
!! INPUTS
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables
!!  iout=unit number for echoed output
!!  jdtset_(0:ndtset_alloc)=list of dataset indices.
!!  ndtset_alloc=govern second dimension of intarr and dprarr
!!  prtvol_glob= if 0, minimal output volume, if 1, no restriction.
!!  results_out(0:ndtset_alloc)=<type results_out_type>contains the results
!!   needed for outvars, including occ, an evolving variable
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      outvar_o_z
!!
!! CHILDREN
!!      appdig
!!
!! SOURCE

subroutine prtocc(dtsets,iout,jdtset_,mxvals,ndtset_alloc,nimagem,prtvol_glob,results_out,strimg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,ndtset_alloc,prtvol_glob
!arrays
 integer,intent(in) :: jdtset_(0:ndtset_alloc)
 integer,intent(in) :: nimagem(0:ndtset_alloc)
 type(ab_dimensions),intent(in) :: mxvals
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
 type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)
 character(len=8),intent(in) :: strimg(mxvals%nimage)

!Local variables-------------------------------
 character(len=*), parameter :: f_occ    ="(1x,a16,1x,(t22,6f10.6))"
 character(len=*), parameter :: f_occa   ="(1x,a16,a,1x,(t22,6f10.6))"
 character(len=*), parameter :: token='occ'
!scalars
 integer,parameter :: nkpt_max=50
 integer :: generic,iban,idtset,ikpsp,ikpt,isppol,jdtset,multi,multi_nband
 integer :: multi_nimage
 integer :: multi_nkpt,multi_nsppol,multi_occopt,nban,nkpt,nkpt_eff
 integer :: multi_tsmear
 integer :: print,tnkpt
 logical, allocatable :: test_multiimages(:)
 character(len=4) :: appen
 character(len=16) :: keywd
 character(len=500) :: message

! *************************************************************************

 if(ndtset_alloc<1)then
   write(message, '(a,i0,a)' )' ndtset_alloc=',ndtset_alloc,', while it should be >= 1.'
   MSG_BUG(message)
 end if

 if(ndtset_alloc>9999)then
   write(message, '(a,i0,a)' )' ndtset_alloc=',ndtset_alloc,', while it must be lower than 100.'
   MSG_BUG(message)
 end if

!It is important to take iscf into account, since when it is -2, occupation numbers must be ignored

 multi_occopt=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%occopt/=dtsets(idtset)%occopt .and. dtsets(idtset)%iscf/=-2 )multi_occopt=1
   end do
 end if

 multi_tsmear=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%tsmear/=dtsets(idtset)%tsmear .and. dtsets(idtset)%iscf/=-2 )multi_tsmear=1
   end do
 end if

 multi_nkpt=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nkpt/=dtsets(idtset)%nkpt .and. dtsets(idtset)%iscf/=-2 )multi_nkpt=1
   end do
 end if
 if(multi_nkpt==0)nkpt=dtsets(1)%nkpt

 multi_nsppol=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nsppol/=dtsets(idtset)%nsppol .and. dtsets(idtset)%iscf/=-2 )multi_nsppol=1
   end do
 end if

 if(multi_nsppol==0 .and. multi_nkpt==0)then
   multi_nband=0
   if(ndtset_alloc>1)then
     do idtset=1,ndtset_alloc
       if(dtsets(idtset)%iscf/=-2)then
         do ikpsp=1,dtsets(1)%nkpt*dtsets(1)%nsppol
           if(dtsets(1)%nband(ikpsp)/=dtsets(idtset)%nband(ikpsp))multi_nband=1
         end do
       end if
     end do
   end if
 else
   multi_nband=1
 end if

 multi_nimage=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nimage/=dtsets(idtset)%nimage .and. dtsets(idtset)%iscf/=-2 )multi_nimage=1
   end do
 end if

!DEBUG
! write(std_out,*)' prtocc : 2, multi_nimage= ',multi_nimage
!ENDDEBUG

!Test whether for this variable, the content of different images differ.
!test_multiimages(idtset)=.false. if, for that dataset, the content for different
!images is identical.
 ABI_MALLOC(test_multiimages,(0:ndtset_alloc))
 test_multiimages=.false.
 do idtset=1,ndtset_alloc
   if(nimagem(idtset)>1)then
     nban=sum(dtsets(idtset)%nband(1:dtsets(idtset)%nsppol*dtsets(idtset)%nkpt))
     do iban=1,nban
       if(sum(abs( results_out(idtset)%occ(iban,2:nimagem(idtset))- results_out(idtset)%occ(iban,1)))>tol12)then
         test_multiimages(idtset)=.true.
       end if
     end do
   end if
 end do
 if(nimagem(0)==0)test_multiimages(0)=.true.

!DEBUG
! write(std_out,*)' prtocc : 3, test_multiimages= ',test_multiimages
! write(std_out,*)' prtocc : multi_occopt, multi_nband, multi_nimage=',multi_occopt, multi_nband, multi_nimage
! write(std_out,*)' prtocc : test_multiimages(1:ndtset_alloc)=',test_multiimages(1:ndtset_alloc)
! write(std_out,*)' prtocc : any(test_multiimages(1:ndtset_alloc))=',any(test_multiimages(1:ndtset_alloc))
!ENDDEBUG

!There is a possibility of a single generic occupation-number set (common to all datasets) if
!multi_occopt==0 and multi_nband==0  and (multi_nimage==0  or the content of the different images is always the same)
!This might occur even if occupation numbers differ for different images.
 multi=1
 if(multi_occopt==0 .and. multi_nband==0 .and. (multi_nimage==0 .or. .not. any(test_multiimages(1:ndtset_alloc)))) then
   nban=sum(dtsets(1)%nband(1:dtsets(1)%nsppol*dtsets(1)%nkpt))
   multi=0
   if(ndtset_alloc>1)then
     do idtset=1,ndtset_alloc
       if(dtsets(idtset)%iscf/=-2)then
!        nban counts all bands and kpoints and spins: see above
         do iimage=1,nimagem(idtset)
           if(iimage==1 .or. test_multiimages(idtset))then
             do iban=1,nban
!              Use of tol8, because the format for multi=1 is f16.6, so will not
!              discriminate between relative values, or absolute values that
!              agree within more than 6 digits
               if( abs(results_out(1)%occ(iban,iimage)-results_out(idtset)%occ(iban,iimage)) > tol8) multi=1
             end do
           end if
         end do
       end if
     end do
   end if
 end if

! write(std_out,*)' prtocc : 4, multi= ',multi

!At this stage, if multi==1, the occ must be printed
!if multi==0, then it might be that we have the default values.
!Since the default is all zeros, it only happens when iscf=-2
!Also initialize the number of a idtset that can be used as generic
!(this might not be the case for idtset=1 !)

 generic=0
 print=0
 do idtset=1,ndtset_alloc
   if(dtsets(idtset)%iscf/=-2)then
     print=1
     generic=idtset
   end if
 end do

! write(std_out,*)' prtocc : 5, print= ',print

!Now, print occ in the generic occupation-number set case (occ is independent of the dtset).
 if(print==1 .and. multi==0)then
!  Might restrict the number of k points to be printed
   tnkpt=0
   nkpt_eff=dtsets(1)%nkpt
   if(prtvol_glob==0 .and. nkpt_eff>nkpt_max)then
     nkpt_eff=nkpt_max
     tnkpt=1
   end if

! write(std_out,*)' prtocc : 6, do-loop over iimage '

   do iimage=1,nimagem(generic)
     if(iimage==1 .or. test_multiimages(generic) )then
       keywd=token//trim(strimg(iimage))
!      The quantity of data to be output vary with occopt
       if(dtsets(generic)%occopt>=2)then
         iban=1
         do isppol=1,dtsets(generic)%nsppol
           do ikpt=1,nkpt_eff
             ikpsp=ikpt+dtsets(generic)%nkpt*(isppol-1)
             nban=dtsets(generic)%nband(ikpsp)
             if(ikpsp==1)then
               write(iout, '(1x,a16,1x,(t22,6f10.6))' )&
&               trim(keywd),results_out(generic)%occ(iban:iban+nban-1,iimage)
             else
               write(iout, '((t22,6f10.6))' )results_out(generic)%occ(iban:iban+nban-1,iimage)
             end if
             iban=iban+nban
           end do
           if(tnkpt==1) write(iout,'(23x,a)' ) 'prtocc : prtvol=0, do not print more k-points.'
         end do
       else
!        The number of bands is identical for all k points and spin
         nban=dtsets(generic)%nband(1)
         write(iout, '(1x,a16,1x,(t22,6f10.6))' )trim(keywd),results_out(generic)%occ(1:nban,iimage)
!        if occopt==1, the occ might differ with the spin
         if(dtsets(generic)%nsppol/=1)then
           write(iout,'((t22,6f10.6))')results_out(generic)%occ(nban*dtsets(generic)%nkpt+1:&
&           nban*dtsets(generic)%nkpt+nban,iimage)
         end if
       end if
     end if
   end do
 end if

! write(std_out,*)' prtocc : 7, finished do-loop over iimage '

!Now, print occ in the other cases (occ depends on the dataset)
 if(print==1 .and. multi==1)then
   do idtset=1,ndtset_alloc
!    Might restrict the number of k points to be printed
     tnkpt=0
     nkpt_eff=dtsets(idtset)%nkpt
     if(prtvol_glob==0 .and. nkpt_eff>nkpt_max)then
       nkpt_eff=nkpt_max
       tnkpt=1
     end if
     if(dtsets(idtset)%iscf/=-2)then
       jdtset=jdtset_(idtset)
       call appdig(jdtset,'',appen)
       do iimage=1,nimagem(idtset)
         if(iimage==1 .or. test_multiimages(idtset) )then
           keywd=trim(token)//trim(strimg(iimage))
!          The quantity of data to be output vary with occopt
           if(dtsets(idtset)%occopt>=2)then
             iban=1
             do isppol=1,dtsets(idtset)%nsppol
               do ikpt=1,nkpt_eff
                 ikpsp=ikpt+dtsets(idtset)%nkpt*(isppol-1)
                 nban=dtsets(idtset)%nband(ikpsp)
                 if(ikpsp==1)then
                   write(iout, '(1x,a16,a,1x,(t22,6f10.6))' )&
&                   trim(keywd),appen,results_out(idtset)%occ(iban:iban+nban-1,iimage)
                 else
                   write(iout, '((t22,6f10.6))' )results_out(idtset)%occ(iban:iban+nban-1,iimage)
                 end if
                 iban=iban+nban
               end do
               if(tnkpt==1) write(iout,'(23x,a)' ) 'prtocc : prtvol=0, do not print more k-points.'
             end do
           else
!            The number of bands is identical for all k points and spin
             nban=dtsets(idtset)%nband(1)
             write(iout, '(1x,a16,a,1x,(t22,6f10.6))' )&
&             trim(keywd),appen,results_out(idtset)%occ(1:nban,iimage)
!            if occopt==1, the occ might differ with the spin
             if(dtsets(idtset)%nsppol/=1)then
               write(iout, '((t22,6f10.6))' ) &
&               results_out(idtset)%occ(nban*dtsets(idtset)%nkpt+1:nban*dtsets(idtset)%nkpt+nban,iimage)
             end if
           end if
         end if
       enddo
     end if
!    Endloop on idtset
   end do
 end if

 ABI_DEALLOCATE(test_multiimages)

end subroutine prtocc
!!***

end subroutine outvar_o_z
!!***

end module m_outvar_o_z
!!***
