!!****m* ABINIT/m_outvar_a_h
!! NAME
!!  m_outvar_a_h
!!
!! FUNCTION
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

module m_outvar_a_h

 use defs_basis
 use m_abicore
 use m_results_out
 use m_dtset

 use m_parser,  only : prttagm, prttagm_images, ab_dimensions

 implicit none

 private
!!***

 public :: outvar_a_h
!!***

contains
!!***

!!****f* ABINIT/outvar_a_h
!! NAME
!! outvar_a_h
!!
!! FUNCTION
!! Echo variables between acell and gw_ ... (by alphabetic order) for the ABINIT code.
!!
!! INPUTS
!!  choice= 1 if echo of preprocessed variables, 2 if echo after call driver
!!  dmatpuflag=flag controlling the use of an initial density matrix in PAW+U (max. value over datasets)
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
!!         natvshift  =maximal value of input natvshift for all the datasets
!!         nconeq     =maximal value of input nconeq for all the datasets
!!         nimage     =maximal value of input nimage for all the datasets
!!         nkpt       =maximal value of input nkpt for all the datasets
!!         nkptgw     =maximal value of input nkptgw for all the datasets
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
!!  Note that acell, occ, rprim, xred and vel might have been modified by the
!!  computation, so that their values if choice=1 or choice=2 will differ.
!!
!! PARENTS
!!      m_outvars
!!
!! CHILDREN
!!      prttagm,prttagm_images
!!
!! SOURCE

subroutine outvar_a_h (choice,dmatpuflag,dtsets,iout,&
& jdtset_,marr,multivals,mxvals,ncid,ndtset,ndtset_alloc,&
& results_out,strimg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,dmatpuflag,iout,marr,ndtset
 integer,intent(in) :: ndtset_alloc,ncid
!arrays
 integer,intent(in) :: jdtset_(0:ndtset_alloc)
 type(ab_dimensions),intent(in) :: multivals,mxvals
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
 type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)
 character(len=8),intent(in) :: strimg(mxvals%nimage)

!Local variables-------------------------------
!scalars
 integer,parameter :: nkpt_max=50
 integer :: defo,idtset,ii,iimage,ga_n_rules,nn
 integer :: lpawu1,narr,mxnsp
 integer :: natom,nimfrqs,nimage
 integer :: ntypalch,ntypat,print_constraint,size1,size2,tmpimg0
 logical :: compute_static_images
 real(dp) :: cpus
 character(len=1) :: firstchar_fftalg,firstchar_gpu
 character(len=14) :: str_hyb
!arrays
 integer,allocatable :: narrm(:)
 integer,allocatable :: nimagem(:),prtimg(:,:)
 integer,allocatable :: intarr(:,:)
 real(dp),allocatable :: dprarr(:,:),dprarr_images(:,:,:)

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

!if(multivals%ga_n_rules==0)ga_n_rules=dtsets(1)%ga_n_rules
 ga_n_rules=dtsets(1)%ga_n_rules
!if(multivals%natom==0)natom=dtsets(1)%natom
 natom=dtsets(1)%natom
!if(multivals%nimage==0)nimage=dtsets(1)%nimage
 nimage=dtsets(1)%nimage

 nimfrqs=dtsets(1)%cd_customnimfrqs
!if(multivals%ntypalch==0)ntypalch=dtsets(1)%ntypalch
 ntypalch=dtsets(1)%ntypalch
!if(multivals%ntypat==0)ntypat=dtsets(1)%ntypat
 ntypat=dtsets(1)%ntypat

!###########################################################
!### 03. Print all the input variables (A)
!##

 intarr(1,:)=dtsets(:)%iomode
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'iomode','INT',0,firstchar="-")

 intarr(1,:)=dtsets(:)%accuracy
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'accuracy','INT',0)

!acell
 prtimg(:,:)=1
 do idtset=0,ndtset_alloc
   narrm(idtset)=3
   do iimage=1,nimagem(idtset)
     if (narrm(idtset)>0) then
       dprarr_images(1:narrm(idtset),iimage,idtset)=results_out(idtset)%acell(1:3,iimage)
     end if
   end do
 end do
 call prttagm_images(dprarr_images,iout,jdtset_,2,marr,narrm,ncid,ndtset_alloc,'acell','LEN',&
   mxvals%nimage,nimagem,ndtset,prtimg,strimg)

!adpimd and adpimd_gamma
 intarr(1,:)=dtsets(:)%adpimd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'adpimd','INT',0)

 dprarr(1,:)=dtsets(:)%adpimd_gamma
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'adpimd_gamma','DPR',0)

!algalch
 narr=ntypalch                      ! default size for all datasets
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=dtsets(idtset)%ntypalch
   if(idtset==0)narrm(idtset)=mxvals%ntypat
   intarr(1:narrm(idtset),idtset)=dtsets(idtset)%algalch(1:narrm(idtset))
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
   narrm,ncid,ndtset_alloc,'algalch','INT',multivals%ntypalch)

!amu
 prtimg(:,:)=1
 do idtset=0,ndtset_alloc
   if(idtset/=0)then
     size1=dtsets(idtset)%ntypat
   else
     size1=mxvals%ntypat
   end if
   narrm(idtset)=size1
   do iimage=1,nimagem(idtset)
     if (narrm(idtset)>0) then
       dprarr_images(1:narrm(idtset),iimage,idtset)=results_out(idtset)%amu(1:size1,iimage)
     end if
   end do
 end do
 call prttagm_images(dprarr_images,iout,jdtset_,1,marr,narrm,ncid,ndtset_alloc,'amu','DPR',&
& mxvals%nimage,nimagem,ndtset,prtimg,strimg,forceprint=2)

 intarr(1,:)=dtsets(:)%asr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'asr','INT',0)

!atvshift
 if(mxvals%natpawu>0)then
   narr=dtsets(1)%natvshift*dtsets(1)%nsppol*mxvals%natom ! default size for all datasets
   do idtset=0,ndtset_alloc       ! specific size for each dataset
     if(idtset/=0)then
       narrm(idtset)=dtsets(idtset)%natvshift*dtsets(idtset)%nsppol*mxvals%natom
       if(narrm(idtset)/=0)&
&       dprarr(1:narrm(idtset),idtset)=&
&       reshape(dtsets(idtset)%atvshift(1:dtsets(idtset)%natvshift,&
&       1:dtsets(idtset)%nsppol,1:mxvals%natom),&
&       (/ narrm(idtset) /) )
     else
       narrm(idtset)=mxvals%natvshift*mxvals%nsppol*mxvals%natom
       if(narrm(idtset)/=0)&
&       dprarr(1:narrm(idtset),idtset)=&
&       reshape(dtsets(idtset)%atvshift(1:mxvals%natvshift,&
&       1:mxvals%nsppol,1:mxvals%natom),&
&       (/ narrm(idtset) /) )
     end if
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,5,marr,narr,&
&   narrm,ncid,ndtset_alloc,'atvshift','DPR',&
&   multivals%natvshift+multivals%nsppol+multivals%natom)
 end if

 intarr(1,:)=dtsets(:)%autoparal
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'autoparal','INT',0)

 intarr(1,:)=dtsets(:)%auxc_ixc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'auxc_ixc','INT',0)

 dprarr(1,:)=dtsets(:)%auxc_scal
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'auxc_scal','DPR',0)

 intarr(1,:)=dtsets(:)%awtr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'awtr','INT',0)

!###########################################################
!### 03. Print all the input variables (B)
!##

 intarr(1,:)=dtsets(:)%bandpp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bandpp','INT',0)

 intarr(1,:)=dtsets(:)%bdberry(1)
 intarr(2,:)=dtsets(:)%bdberry(2)
 intarr(3,:)=dtsets(:)%bdberry(3)
 intarr(4,:)=dtsets(:)%bdberry(4)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,4,narrm,ncid,ndtset_alloc,'bdberry','INT',0)

 intarr(1,:)=dtsets(:)%bdeigrf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bdeigrf','INT',0)

!bdgw
 narr=2*dtsets(1)%nkptgw*dtsets(1)%nsppol ! default size for all datasets
 do idtset=0,ndtset_alloc        ! specific size for each dataset
   if(idtset/=0)then
     narrm(idtset)=2*dtsets(idtset)%nkptgw*dtsets(idtset)%nsppol
     if (narrm(idtset)>0)&
&     intarr(1:narrm(idtset),idtset)=&
&     reshape(dtsets(idtset)%bdgw(1:2,1:dtsets(idtset)%nkptgw,1:dtsets(idtset)%nsppol),(/narrm(idtset)/))
   else
     narrm(idtset)=2*mxvals%nkptgw*mxvals%nsppol
     if (narrm(idtset)>0)&
&     intarr(1:narrm(idtset),idtset)=&
&     reshape(dtsets(idtset)%bdgw(1:2,1:mxvals%nkptgw,1:mxvals%nsppol),(/ narrm(idtset) /) )
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,narr,&
& narrm,ncid,ndtset_alloc,'bdgw','INT',multivals%nkptgw+multivals%nsppol)

 intarr(1,:)=dtsets(:)%berryopt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'berryopt','INT',0)

 intarr(1,:)=dtsets(:)%berrysav
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'berrysav','INT',0)

 intarr(1,:)=dtsets(:)%berrystep
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'berrystep','INT',0)

 dprarr(1,:)=dtsets(:)%bfield(1)
 dprarr(2,:)=dtsets(:)%bfield(2)
 dprarr(3,:)=dtsets(:)%bfield(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'bfield','DPR',0)

 dprarr(1,:)=dtsets(:)%bmass
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'bmass','DPR',0)

 dprarr(1,:)=dtsets(:)%boxcenter(1)
 dprarr(2,:)=dtsets(:)%boxcenter(2)
 dprarr(3,:)=dtsets(:)%boxcenter(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'boxcenter','DPR',0)

 dprarr(1,:)=dtsets(:)%boxcutmin
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'boxcutmin','DPR',0)

 intarr(1,:)=dtsets(:)%brav
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'brav','INT',0)

 intarr(1,:)=dtsets(:)%bs_algorithm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_algorithm','INT',0)

 intarr(1,:)=dtsets(:)%bs_calctype
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_calctype','INT',0)

 intarr(1,:)=dtsets(:)%bs_coulomb_term
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_coulomb_term','INT',0)

 intarr(1,:)=dtsets(:)%bs_coupling
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_coupling','INT',0)

 do idtset=0,ndtset_alloc
   dprarr(1:2,idtset)=dtsets(idtset)%bs_eh_cutoff(1:2)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,2,narrm,ncid,ndtset_alloc,'bs_eh_cutoff','ENE',0)

 intarr(1,:)=dtsets(:)%bs_exchange_term
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_exchange_term','INT',0)

 do idtset=0,ndtset_alloc
   dprarr(1:3,idtset)=dtsets(idtset)%bs_freq_mesh(1:3)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'bs_freq_mesh','ENE',0)

 intarr(1,:)=dtsets(:)%bs_haydock_niter
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_haydock_niter','INT',0)

 do idtset=0,ndtset_alloc
   dprarr(1:2,idtset)=dtsets(idtset)%bs_haydock_tol(:)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,2,narrm,ncid,ndtset_alloc,'bs_haydock_tol','DPR',0)

 intarr(1,:)=dtsets(:)%bs_hayd_term
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_hayd_term','INT',0)

 do idtset=0,ndtset_alloc
   intarr(1:3,idtset)=dtsets(idtset)%bs_interp_kmult(1:3)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'bs_interp_kmult','INT',0)

 intarr(1,:)=dtsets(:)%bs_interp_method
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_interp_method','INT',0)

 intarr(1,:)=dtsets(:)%bs_interp_mode
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_interp_mode','INT',0)

 intarr(1,:)=dtsets(:)%bs_interp_prep
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_interp_prep','INT',0)

 intarr(1,:)=dtsets(:)%bs_interp_rl_nb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_interp_rl_nb','INT',0)

!bs_loband
 narr=dtsets(1)%nsppol ! default size for all datasets
 intarr = 0
 do idtset=0,ndtset_alloc        ! specific size for each dataset
   if(idtset/=0)then
     narrm(idtset)=dtsets(idtset)%nsppol
   else
     narrm(idtset)=mxvals%nsppol
   end if
   intarr(1:narrm(idtset),idtset)=dtsets(idtset)%bs_loband(1:narrm(idtset))
 end do

 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,narr,narrm,ncid,ndtset_alloc,'bs_loband','INT',multivals%nsppol)

 intarr(1,:)=dtsets(:)%bs_nstates
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_nstates','INT',0)

 intarr(1,:)=dtsets(:)%builtintest
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'builtintest','INT',0)

 dprarr(1,:)=dtsets(:)%bxctmindg
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'bxctmindg','DPR',0)

!###########################################################
!### 03. Print all the input variables (C)
!##

 if (ANY(dtsets(:)%cd_customnimfrqs/=0)) then
   intarr(1,:)=dtsets(:)%cd_customnimfrqs
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'cd_customnimfrqs','INT',0)
 end if

 intarr(1,:)=dtsets(:)%cd_frqim_method
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'cd_frqim_method','INT',0)

 intarr(1,:)=dtsets(:)%cd_full_grid
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'cd_full_grid','INT',0)

 dprarr(1,:)=dtsets(:)%cd_halfway_freq
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'cd_halfway_freq','ENE',0)

!cd_imfrqs
 narr=mxvals%nimfrqs            ! default size for all datasets
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=dtsets(idtset)%cd_customnimfrqs
   if(idtset==0)narrm(idtset)=mxvals%nimfrqs
   if (narrm(idtset)>0) then
     dprarr(1:narrm(idtset),idtset)=dtsets(idtset)%cd_imfrqs(1:narrm(idtset))
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,6,marr,narr,narrm,ncid,ndtset_alloc,'cd_imfrqs','ENE',multivals%nimfrqs)

 dprarr(1,:)=dtsets(:)%cd_max_freq
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'cd_max_freq','ENE',0)

 if (ANY(dtsets(:)%cd_subset_freq(1)/=0)) then
   intarr(1,:)=dtsets(:)%cd_subset_freq(1)
   intarr(2,:)=dtsets(:)%cd_subset_freq(2)
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,2,narrm,ncid,ndtset_alloc,'cd_subset_freq','INT',0)
 end if

 dprarr(1,:)=dtsets(:)%charge
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'charge','DPR',0)

!chempot
 narr=3*mxvals%nzchempot*mxvals%ntypat ! default size for all datasets
 if(narr/=0)then
   do idtset=0,ndtset_alloc       ! specific size for each dataset
     if(idtset/=0)then
       narrm(idtset)=3*dtsets(idtset)%nzchempot*dtsets(idtset)%ntypat
       if(narrm(idtset)/=0)&
&       dprarr(1:narrm(idtset),idtset)=&
&       reshape(dtsets(idtset)%chempot(1:3,1:dtsets(idtset)%nzchempot,&
&       1:dtsets(idtset)%ntypat),&
&       (/ narrm(idtset) /) )
     else
       narrm(idtset)=3*mxvals%nzchempot*mxvals%ntypat
       if(narrm(idtset)/=0)&
&       dprarr(1:narrm(idtset),idtset)=&
&       reshape(dtsets(idtset)%chempot(1:3,1:mxvals%nzchempot,1:mxvals%ntypat),(/ narrm(idtset) /) )
     end if
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'chempot','DPR',1)
 end if

 intarr(1,:)=dtsets(:)%chkdilatmx
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'chkdilatmx','INT',0)

 intarr(1,:)=dtsets(:)%chkexit
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'chkexit','INT',0)

 intarr(1,:)=dtsets(:)%chkprim
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'chkprim','INT',0)

 intarr(1,:)=dtsets(:)%chksymbreak
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'chksymbreak','INT',0)

 intarr(1,:)=dtsets(:)%chneut
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'chneut','INT',0)

!chrgat
 narr=mxvals%natom              ! default size for all datasets
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=dtsets(idtset)%natom
   if(idtset==0)narrm(idtset)=mxvals%natom
   if (narrm(idtset)>0) dprarr(1:narrm(idtset),idtset)=dtsets(idtset)%chrgat(1:narrm(idtset))
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'chrgat','DPR',multivals%natom)

 intarr(1,:)=dtsets(:)%cineb_start
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'cineb_start','INT',0)

 if(dtsets(1)%cpus>one)then
   cpus=dtsets(1)%cpus
   write(iout,'(1x,a16,1x,1p,t22,g10.2,t25,a)') 'cpus',cpus,'(seconds)'
   write(iout,'(1x,a16,1x,1p,t22,g10.2,t25,a)') 'cpum',cpus/60.0_dp,'(minutes)'
   write(iout,'(1x,a16,1x,1p,t22,g10.2,t25,a)') 'cpuh',cpus/3600.0_dp,'(hours)'
 end if

!constraint_kind
 narr=mxvals%ntypat             ! default size for all datasets
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=dtsets(idtset)%ntypat
   if(idtset==0)narrm(idtset)=mxvals%ntypat
   if (narrm(idtset)>0) intarr(1:narrm(idtset),idtset)=dtsets(idtset)%constraint_kind(1:narrm(idtset))
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'constraint_kind','INT',multivals%ntypat)

!corecs
 narr=mxvals%ntypat             ! default size for all datasets
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=dtsets(idtset)%ntypat
   if(idtset==0)narrm(idtset)=mxvals%ntypat
   if (narrm(idtset)>0) dprarr(1:narrm(idtset),idtset)=dtsets(idtset)%corecs(1:narrm(idtset))
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'corecs','DPR',multivals%ntypat)

!###########################################################
!### 03. Print all the input variables (D)
!##

 dprarr(1,:)=dtsets(:)%ddamp
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'ddamp','DPR',0)

 do idtset=0,ndtset_alloc
   intarr(1:3,idtset)=dtsets(idtset)%ddb_ngqpt
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'ddb_ngqpt','INT',0)

 dprarr(1,:)=dtsets(:)%ddb_shiftq(1)
 dprarr(2,:)=dtsets(:)%ddb_shiftq(2)
 dprarr(3,:)=dtsets(:)%ddb_shiftq(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'ddb_shiftq','DPR',0)

 intarr(1,:)=dtsets(:)%delayperm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'delayperm','INT',0)

 intarr(1,:)=dtsets(:)%densfor_pred
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'densfor_pred','INT',0)

!densty
 narr=mxvals%ntypat              ! default size for all datasets
 do idtset=0,ndtset_alloc        ! specific size for each dataset
   narrm(idtset)=dtsets(idtset)%ntypat
   if(idtset==0)narrm(idtset)=mxvals%ntypat
!  Only one component of densty is used until now
   dprarr(1:narrm(idtset),idtset)=dtsets(idtset)%densty(1:narrm(idtset),1)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'densty','DPR',multivals%ntypat)

 dprarr(1,:)=dtsets(:)%dfield(1)
 dprarr(2,:)=dtsets(:)%dfield(2)
 dprarr(3,:)=dtsets(:)%dfield(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'dfield','DPR',0)

 dprarr(1,:)=dtsets(:)%dfpt_sciss
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'dfpt_sciss','ENE',0)

 dprarr(1,:)=dtsets(:)%diecut
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'diecut','ENE',0)

 dprarr(1,:)=dtsets(:)%diegap
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'diegap','ENE',0)

 dprarr(1,:)=dtsets(:)%dielam
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'dielam','DPR',0)

 dprarr(1,:)=dtsets(:)%dielng
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'dielng','LEN',0)

 dprarr(1,:)=dtsets(:)%diemac
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'diemac','DPR',0)

 dprarr(1,:)=dtsets(:)%diemix
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'diemix','DPR',0)

 if (any(dtsets(1:ndtset_alloc)%diemixmag/=dtsets(1:ndtset_alloc)%diemix)) then
   dprarr(1,:)=dtsets(:)%diemixmag
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'diemixmag','DPR',0)
 end if

 intarr(1,:)=dtsets(:)%diismemory
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'diismemory','INT',0)

 dprarr(1,:)=dtsets(:)%dilatmx
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'dilatmx','DPR',0)

 intarr(1,:)=dtsets(:)%dipdip
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dipdip','INT',0)

!dmatpawu
 if (dmatpuflag==1.and.mxvals%natpawu>0) then
   prtimg(:,:)=1
   do idtset=0,ndtset_alloc
     mxnsp=max(dtsets(idtset)%nsppol,dtsets(idtset)%nspinor)
     lpawu1=maxval(dtsets(idtset)%lpawu(:))
     narrm(idtset)=((2*lpawu1+1)**2)*mxnsp*dtsets(idtset)%natpawu
     do iimage=1,nimagem(idtset)
       if (narrm(idtset)>0) then
         dprarr_images(1:narrm(idtset),iimage,idtset)= &
&         reshape(dtsets(idtset)%dmatpawu(&
&         1:2*lpawu1+1,1:2*lpawu1+1,1:mxnsp,1:dtsets(idtset)%natpawu,iimage),&
&         (/narrm(idtset)/))
       end if
     end do
   end do
   call prttagm_images(dprarr_images,iout,jdtset_,5,marr,narrm,&
     ncid,ndtset_alloc,'dmatpawu','DPR',mxvals%nimage,nimagem,ndtset,prtimg,strimg)
 end if

 intarr(1,:)=dtsets(:)%dmatpuopt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmatpuopt','INT',0)

 intarr(1,:)=dtsets(:)%dmatudiag
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmatudiag','INT',0)

 intarr(1,:)=dtsets(:)%dmftbandf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmftbandf','INT',0)

 intarr(1,:)=dtsets(:)%dmftbandi
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmftbandi','INT',0)

 intarr(1,:)=dtsets(:)%dmftcheck
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmftcheck','INT',0)

 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_dc','INT',0)

 intarr(1,:)=dtsets(:)%dmft_iter
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_iter','INT',0)

! intarr(1,:)=dtsets(:)%dmft_kspectralfunc
! call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_kspectralfunc','INT',0)

 dprarr(1,:)=dtsets(:)%dmft_mxsf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_mxsf','DPR',0)

 intarr(1,:)=dtsets(:)%dmft_nwli
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_nwli','INT',0)

 intarr(1,:)=dtsets(:)%dmft_nwlo
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_nwlo','INT',0)

 intarr(1,:)=dtsets(:)%dmft_occnd_imag
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_occnd_imag','INT',0)

 intarr(1,:)=dtsets(:)%dmft_read_occnd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_read_occnd','INT',0)

 intarr(1,:)=dtsets(:)%dmft_rslf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_rslf','INT',0)

 intarr(1,:)=dtsets(:)%dmft_solv
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_solv','INT',0)

 dprarr(1,:)=dtsets(:)%dmft_tolfreq
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_tolfreq','DPR',0)

 dprarr(1,:)=dtsets(:)%dmft_tollc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_tollc','DPR',0)

 dprarr(1,:)=dtsets(:)%dmft_charge_prec
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_charge_prec','DPR',0)

 dprarr(1,:)=dtsets(:)%dosdeltae
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'dosdeltae','ENE',0)

 dprarr(1,:)=dtsets(:)%dtion
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'dtion','DPR',0,forceprint=2)

 intarr(1,:)=dtsets(:)%dvdb_add_lr
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'dvdb_add_lr','INT',0)

 dprarr(1,:)=dtsets(:)%dvdb_qcache_mb
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'dvdb_qcache_mb','DPR',0)

 dprarr(1,:)=dtsets(:)%dvdb_qdamp
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'dvdb_qdamp','DPR',0)

 intarr(1,:)=dtsets(:)%dvdb_rspace_cell
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'dvdb_rspace_cell','INT',0)

!dynimage
 intarr(1:marr,0)=1                 ! default value
 narr=nimage                        ! default size for all datasets
 do idtset=1,ndtset_alloc           ! specific size and array for each dataset
   narrm(idtset)=dtsets(idtset)%nimage
   intarr(1:narrm(idtset),idtset)=dtsets(idtset)%dynimage(1:narrm(idtset))
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'dynimage','INT',multivals%nimage)

!Variables for nonlinear response
 intarr(1,:)=dtsets(:)%d3e_pert1_atpol(1)
 intarr(2,:)=dtsets(:)%d3e_pert1_atpol(2)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,2,narrm,ncid,ndtset_alloc,'d3e_pert1_atpol','INT',0)

 intarr(1,:)=dtsets(:)%d3e_pert1_dir(1)
 intarr(2,:)=dtsets(:)%d3e_pert1_dir(2)
 intarr(3,:)=dtsets(:)%d3e_pert1_dir(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'d3e_pert1_dir','INT',0)

 intarr(1,:)=dtsets(:)%d3e_pert1_elfd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'d3e_pert1_elfd','INT',0)

 intarr(1,:)=dtsets(:)%d3e_pert1_phon
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'d3e_pert1_phon','INT',0)

 intarr(1,:)=dtsets(:)%d3e_pert2_atpol(1)
 intarr(2,:)=dtsets(:)%d3e_pert2_atpol(2)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,2,narrm,ncid,ndtset_alloc,'d3e_pert2_atpol','INT',0)

 intarr(1,:)=dtsets(:)%d3e_pert2_dir(1)
 intarr(2,:)=dtsets(:)%d3e_pert2_dir(2)
 intarr(3,:)=dtsets(:)%d3e_pert2_dir(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'d3e_pert2_dir','INT',0)

 intarr(1,:)=dtsets(:)%d3e_pert2_elfd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'d3e_pert2_elfd','INT',0)

 intarr(1,:)=dtsets(:)%d3e_pert2_phon
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'d3e_pert2_phon','INT',0)

 intarr(1,:)=dtsets(:)%d3e_pert2_strs
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'d3e_pert2_strs','INT',0)

 intarr(1,:)=dtsets(:)%d3e_pert3_atpol(1)
 intarr(2,:)=dtsets(:)%d3e_pert3_atpol(2)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,2,narrm,ncid,ndtset_alloc,'d3e_pert3_atpol','INT',0)

 intarr(1,:)=dtsets(:)%d3e_pert3_dir(1)
 intarr(2,:)=dtsets(:)%d3e_pert3_dir(2)
 intarr(3,:)=dtsets(:)%d3e_pert3_dir(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'d3e_pert3_dir','INT',0)

 intarr(1,:)=dtsets(:)%d3e_pert3_elfd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'d3e_pert3_elfd','INT',0)

 intarr(1,:)=dtsets(:)%d3e_pert3_phon
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'d3e_pert3_phon','INT',0)

!###########################################################
!### 03. Print all the input variables (E)
!##

 dprarr(1,:)=dtsets(:)%ecut
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'ecut','ENE',0)

 dprarr(1,:)=dtsets(:)%ecuteps
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'ecuteps','ENE',0)

 dprarr(1,:)=dtsets(:)%ecutsigx
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'ecutsigx','ENE',0)

 dprarr(1,:)=dtsets(:)%ecutsm
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'ecutsm','ENE',0)

 dprarr(1,:)=dtsets(:)%ecutwfn
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'ecutwfn','ENE',0)

 dprarr(1,:)=dtsets(:)%effmass_free
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'effmass_free','DPR',0)

 dprarr(1,:)=dtsets(:)%efield(1)
 dprarr(2,:)=dtsets(:)%efield(2)
 dprarr(3,:)=dtsets(:)%efield(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'efield','DPR',0)

 nn = size(dtsets(0)%einterp)
 do ii=1,nn
   dprarr(ii,:)=dtsets(:)%einterp(ii)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,nn,narrm,ncid,ndtset_alloc,'einterp','DPR',0)

 dprarr(1,:)=dtsets(:)%elph2_imagden
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'elph2_imagden','ENE',0)

 intarr(1,:)=dtsets(:)%enunit
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'enunit','INT',0)

 dprarr(1,:)=dtsets(:)%eph_ecutosc
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'eph_ecutosc','ENE',0)

 dprarr(1,:)=dtsets(:)%eph_phwinfact
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'eph_phwinfact','DPR',0)

 dprarr(1,:)=dtsets(:)%eph_extrael
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'eph_extrael','DPR',0)

 dprarr(1,:)=dtsets(:)%eph_fermie
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'eph_fermie','ENE',0)

 intarr(1,:)=dtsets(:)%eph_frohlichm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'eph_frohlichm','INT',0)

 dprarr(1,:)=dtsets(:)%eph_fsewin
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'eph_fsewin','ENE',0)

 dprarr(1,:)=dtsets(:)%eph_fsmear
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'eph_fsmear','ENE',0)

 intarr(1,:)=dtsets(:)%eph_intmeth
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'eph_intmeth','INT',0)

 dprarr(1,:)=dtsets(:)%eph_mustar
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'eph_mustar','DPR',0)

 do idtset=0,ndtset_alloc
   intarr(1:3,idtset)=dtsets(idtset)%eph_ngqpt_fine
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'eph_ngqpt_fine','INT',0)

 narr = size(dtsets(0)%eph_np_pqbks)
 do idtset=0,ndtset_alloc
   intarr(1:narr,idtset) = dtsets(idtset)%eph_np_pqbks
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'eph_np_pqbks','INT',0, firstchar="-")

 intarr(1,:)=dtsets(:)%eph_phrange(1)
 intarr(2,:)=dtsets(:)%eph_phrange(2)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,2,narrm,ncid,ndtset_alloc,'eph_phrange','INT',0)

 intarr(1,:)=dtsets(:)%eph_restart
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'eph_restart','INT',0)

 intarr(1,:)=dtsets(:)%eph_stern
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'eph_stern','INT',0)

 intarr(1,:)=dtsets(:)%eph_task
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'eph_task','INT',0)

 dprarr(1,:)=dtsets(:)%eph_tols_idelta(1)
 dprarr(2,:)=dtsets(:)%eph_tols_idelta(2)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,2,narrm,ncid,ndtset_alloc,'eph_tols_idelta','DPR',0)

 intarr(1,:)=dtsets(:)%eph_transport
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'eph_transport','INT',0)

 intarr(1,:)=dtsets(:)%eph_use_ftinterp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'eph_use_ftinterp','INT',0)


 dprarr(1,:)=dtsets(:)%eshift
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'eshift','ENE',0)

 dprarr(1,:)=dtsets(:)%esmear
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'esmear','ENE',0)

!etotal
 if(choice==2)then
   prtimg(:,:)=1
   do idtset=0,ndtset_alloc       ! specific size for each dataset
     compute_static_images=(dtsets(idtset)%istatimg>0)
     narrm(idtset)=1

     if(dtsets(idtset)%iscf>=0 .or. dtsets(idtset)%iscf==-3)then
       do iimage=1,dtsets(idtset)%nimage
         if (narrm(idtset)>0) then
           dprarr_images(1:narrm(idtset),iimage,idtset)=results_out(idtset)%etotal(iimage)
         end if
         if(.not.(dtsets(idtset)%dynimage(iimage)==1.or.compute_static_images))then
           prtimg(iimage,idtset)=0
         end if
       end do
     else
       narrm(idtset)=0
     end if
   end do
!  This is a trick to force printing of etotal even if zero, still not destroying the value of nimagem(0).
   tmpimg0=nimagem(0)
   nimagem(0)=0
   call prttagm_images(dprarr_images,iout,jdtset_,2,marr,narrm,ncid,ndtset_alloc,'etotal','DPR',&
     mxvals%nimage,nimagem,ndtset,prtimg,strimg)
   nimagem(0)=tmpimg0
 end if

 dprarr(1,:)=dtsets(:)%exchmix
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'exchmix','DPR',0)

 intarr(1,:)=dtsets(:)%exchn2n3d
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'exchn2n3d','INT',0)

 intarr(1,:)=dtsets(:)%extrapwf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'extrapwf','INT',0)

!###########################################################
!### 03. Print all the input variables (F)
!##

!fcart
 if(choice==2)then
   prtimg(:,:)=1
   do idtset=0,ndtset_alloc       ! specific size for each dataset
     compute_static_images=(dtsets(idtset)%istatimg>0)
     size2=dtsets(idtset)%natom
     if(idtset==0)size2=0
     narrm(idtset)=3*size2
     if(dtsets(idtset)%iscf>=0 .or. idtset==0)then
       do iimage=1,dtsets(idtset)%nimage
         if (narrm(idtset)>0) then
           dprarr_images(1:narrm(idtset),iimage,idtset)=&
&           reshape(results_out(idtset)%fcart(1:3,1:size2,iimage),(/ narrm(idtset) /) )
         end if
         if(.not.(dtsets(idtset)%dynimage(iimage)==1.or.compute_static_images))then
           prtimg(iimage,idtset)=0
         end if
       end do
     else
       narrm(idtset)=0
     end if
   end do
!  This is a trick to force printing of fcart even if zero, still not destroying the value of nimagem(0).
   tmpimg0=nimagem(0)
   nimagem(0)=0
   call prttagm_images(dprarr_images,iout,jdtset_,2,marr,narrm,ncid,ndtset_alloc,'fcart','DPR',&
     mxvals%nimage,nimagem,ndtset,prtimg,strimg)
   nimagem(0)=tmpimg0
 end if

 dprarr(1,:)=dtsets(:)%fermie_nest
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'fermie_nest','DPR',0)

 firstchar_fftalg = "_"
 intarr(1,:)=dtsets(:)%ngfft(7)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'fftalg','INT',0,firstchar="-",forceprint=3)

 intarr(1,:)=dtsets(:)%ngfft(8)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'fftcache','INT',0)

 intarr(1,:)=dtsets(:)%fftgw
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'fftgw','INT',0)

 intarr(1,:)=dtsets(:)%fockoptmix
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'fockoptmix','INT',0)

 dprarr(1,:)=dtsets(:)%focktoldfe
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'focktoldfe','DPR',0)

 intarr(1,:)=dtsets(:)%fockdownsampling(1)
 intarr(2,:)=dtsets(:)%fockdownsampling(2)
 intarr(3,:)=dtsets(:)%fockdownsampling(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'fockdownsampling','INT',0)

 intarr(1,:)=dtsets(:)%fock_icutcoul
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'fock_icsing','INT',0)

 dprarr(1,:)=dtsets(:)%freqim_alpha
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'freqim_alpha','DPR',0)

 dprarr(1,:)=dtsets(:)%freqremax
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'freqremax','ENE',0)

 dprarr(1,:)=dtsets(:)%freqremin
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'freqremin','ENE',0)

 dprarr(1,:)=dtsets(:)%freqspmax
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'freqspmax','ENE',0)

 dprarr(1,:)=dtsets(:)%freqspmin
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'freqspmin','ENE',0)

 dprarr(1,:)=dtsets(:)%friction
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'friction','DPR',0)

 intarr(1,:)=dtsets(:)%frzfermi
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'frzfermi','INT',0)

 dprarr(1,:)=dtsets(:)%fxcartfactor
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'fxcartfactor','DPR',0)

!f4of2_sla
 narr=mxvals%ntypat                    ! default size for all datasets
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=dtsets(idtset)%ntypat
   if(idtset==0)narrm(idtset)=mxvals%ntypat
   if (narrm(idtset)>0) then
     dprarr(1:narrm(idtset),idtset)=dtsets(idtset)%f4of2_sla(1:narrm(idtset))
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'f4of2_sla','DPR',multivals%ntypat)

!f6of2_sla
 narr=mxvals%ntypat                    ! default size for all datasets
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=dtsets(idtset)%ntypat
   if(idtset==0)narrm(idtset)=mxvals%ntypat
   if (narrm(idtset)>0) then
     dprarr(1:narrm(idtset),idtset)=dtsets(idtset)%f6of2_sla(1:narrm(idtset))
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'f6of2_sla','DPR',multivals%ntypat)

!###########################################################
!### 03. Print all the input variables (G)
!##

 intarr(1,:)=dtsets(:)%ga_algor
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ga_algor','INT',0)

 intarr(1,:)=dtsets(:)%ga_fitness
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ga_fitness','INT',0)

 intarr(1,:)=dtsets(:)%ga_n_rules
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ga_n_rules','INT',0)

 dprarr(1,:)=dtsets(:)%ga_opt_percent
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ga_opt_percent','DPR',0)

!ga_rules
 narr=ga_n_rules                    ! default size for all datasets
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   narrm(idtset)=dtsets(idtset)%ga_n_rules
   if(idtset==0)narrm(idtset)=mxvals%ga_n_rules
   intarr(1:narrm(idtset),idtset)=dtsets(idtset)%ga_rules(1:narrm(idtset))
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'ga_rules','INT',multivals%ga_n_rules)

 intarr(1,:)=dtsets(:)%getbscoup
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getbscoup','INT',0)

 intarr(1,:)=dtsets(:)%getbseig
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getbseig','INT',0)

 intarr(1,:)=dtsets(:)%getbsreso
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getbsreso','INT',0)

 intarr(1,:)=dtsets(:)%getcell
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getcell','INT',0)

 intarr(1,:)=dtsets(:)%getddb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getddb','INT',0)

 intarr(1,:)=dtsets(:)%getddk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getddk','INT',0)

 intarr(1,:)=dtsets(:)%getdelfd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getdelfd','INT',0)

 intarr(1,:)=dtsets(:)%getdkdk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getdkdk','INT',0)

 intarr(1,:)=dtsets(:)%getdkde
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getdkde','INT',0)

 intarr(1,:)=dtsets(:)%getden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getden','INT',0)

 intarr(1,:)=dtsets(:)%getdvdb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getdvdb','INT',0)

 intarr(1,:)=dtsets(:)%getefmas
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getefmas','INT',0)

 intarr(1,:)=dtsets(:)%getgam_eig2nkq
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getgam_eig2nkq','INT',0)

 intarr(1,:)=dtsets(:)%gethaydock
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gethaydock','INT',0)

 intarr(1,:)=dtsets(:)%getocc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getocc','INT',0)

 intarr(1,:)=dtsets(:)%getpawden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getpawden','INT',0)

 intarr(1,:)=dtsets(:)%getqps
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getqps','INT',0)

 intarr(1,:)=dtsets(:)%getscr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getscr','INT',0)

 intarr(1,:)=dtsets(:)%getsuscep
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getsuscep','INT',0)

 intarr(1,:)=dtsets(:)%getvel
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getvel','INT',0)

 intarr(1,:)=dtsets(:)%getwfk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getwfk','INT',0)

 intarr(1,:)=dtsets(:)%getwfkfine
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getwfkfine','INT',0)

 intarr(1,:)=dtsets(:)%getwfq
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getwfq','INT',0)

 intarr(1,:)=dtsets(:)%getxcart
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getxcart','INT',0)

 intarr(1,:)=dtsets(:)%getxred
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getxred','INT',0)

 intarr(1,:)=dtsets(:)%get1den
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'get1den','INT',0)

 intarr(1,:)=dtsets(:)%get1wf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'get1wf','INT',0)

 intarr(1,:)=dtsets(:)%goprecon
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'goprecon','INT',0)

 dprarr(1,:)=dtsets(:)%goprecprm(1)
 dprarr(2,:)=dtsets(:)%goprecprm(2)
 dprarr(3,:)=dtsets(:)%goprecprm(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'goprecprm','DPR',0)

 intarr(1,:)=dtsets(:)%gpu_devices(1) ; intarr(2,:)=dtsets(:)%gpu_devices(2)
 intarr(3,:)=dtsets(:)%gpu_devices(3) ; intarr(4,:)=dtsets(:)%gpu_devices(4)
 intarr(5,:)=dtsets(:)%gpu_devices(5)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,5,narrm,ncid,ndtset_alloc,'gpu_devices','INT',0)

 intarr(1,:)=dtsets(:)%gpu_linalg_limit
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gpu_linalg_limit','INT',0)

!grchrg
 print_constraint=0
 do idtset=1,ndtset_alloc
   if(any(dtsets(idtset)%constraint_kind(:)>=10))print_constraint=1
 enddo
 if(print_constraint==1)then
!if(any(dtsets(1:ndtset_alloc)%constraint_kind(:)>=10))then
   if(choice==2)then
     prtimg(:,:)=1
     do idtset=0,ndtset_alloc       ! specific size for each dataset
       compute_static_images=(dtsets(idtset)%istatimg>0)
       size2=dtsets(idtset)%natom
       if(idtset==0)size2=0
       narrm(idtset)=size2
       if(dtsets(idtset)%iscf>=0 .or. idtset==0)then
         do iimage=1,dtsets(idtset)%nimage
           if (narrm(idtset)>0) then
             dprarr_images(1:narrm(idtset),iimage,idtset)=&
&             results_out(idtset)%intgres(1,1:size2,iimage)
           end if
           if(.not.(dtsets(idtset)%dynimage(iimage)==1.or.compute_static_images))then
             prtimg(iimage,idtset)=0
           end if
         end do
       else
         narrm(idtset)=0
       end if
     end do
!    This is a trick to force printing of fcart even if zero, still not destroying the value of nimagem(0).
     tmpimg0=nimagem(0)
     nimagem(0)=0
     call prttagm_images(dprarr_images,iout,jdtset_,1,marr,narrm,ncid,ndtset_alloc,'grchrg','DPR',&
       mxvals%nimage,nimagem,ndtset,prtimg,strimg)
     nimagem(0)=tmpimg0
   end if
 endif

!grspin
 print_constraint=0
 do idtset=1,ndtset_alloc
   if(any(mod(dtsets(idtset)%constraint_kind(:),10)>0))print_constraint=1
 enddo
 if(print_constraint==1)then
!if(any(mod(dtsets(1:ndtset_alloc)%constraint_kind(:),10)/=0))then
   if(choice==2)then
     prtimg(:,:)=1
     do idtset=0,ndtset_alloc       ! specific size for each dataset
       compute_static_images=(dtsets(idtset)%istatimg>0)
       size2=dtsets(idtset)%natom
       if(idtset==0)size2=0
       narrm(idtset)=3*size2
       if(dtsets(idtset)%iscf>=0 .or. idtset==0)then
         do iimage=1,dtsets(idtset)%nimage
           if (narrm(idtset)>0) then
             dprarr_images(1:narrm(idtset),iimage,idtset)=&
&             reshape(results_out(idtset)%intgres(2:4,1:size2,iimage),(/ narrm(idtset) /) )
           end if
           if(.not.(dtsets(idtset)%dynimage(iimage)==1.or.compute_static_images))then
             prtimg(iimage,idtset)=0
           end if
         end do
       else
         narrm(idtset)=0
       end if
     end do
!    This is a trick to force printing of fcart even if zero, still not destroying the value of nimagem(0).
     tmpimg0=nimagem(0)
     nimagem(0)=0
     call prttagm_images(dprarr_images,iout,jdtset_,2,marr,narrm,ncid,ndtset_alloc,'grspin','DPR',&
       mxvals%nimage,nimagem,ndtset,prtimg,strimg)
     nimagem(0)=tmpimg0
   end if
 endif

 intarr(1,:)=dtsets(:)%gwcalctyp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gwcalctyp','INT',0)

 intarr(1,:)=dtsets(:)%gwcomp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gwcomp','INT',0)

 dprarr(1,:)=dtsets(:)%gwencomp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gwencomp','ENE',0)

 intarr(1,:)=dtsets(:)%gwgamma
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gwgamma','INT',0)

 intarr(1,:)=dtsets(:)%gwmem
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gwmem','INT',0)

 intarr(1,:)=dtsets(:)%gwpara
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gwpara','INT',0, firstchar="-")

 intarr(1,:)=dtsets(:)%gwrpacorr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gwrpacorr','INT',0)

!gw_customnfreqsp
!It actually overrides the content of nfreqsp (which is forbidden !) in dtset.
!This is to be cleaned ...
 if (ANY(dtsets(:)%gw_customnfreqsp/=0)) then
   intarr(1,:)=dtsets(:)%gw_customnfreqsp
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_customnfreqsp','INT',0)
 end if

!gw_freqsp
!This is to be cleaned ... See above ...
 narr=mxvals%nfreqsp ! default size for all datasets
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   dprarr(1:narr,idtset)=zero
   narrm(idtset)=dtsets(idtset)%gw_customnfreqsp
   if(idtset==0)narrm(idtset)=mxvals%nfreqsp
   if (narrm(idtset)>0) then
     dprarr(1:narrm(idtset),idtset)=dtsets(idtset)%gw_freqsp(1:narrm(idtset))
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,6,marr,narr,narrm,ncid,ndtset_alloc,'gw_freqsp','ENE',multivals%nfreqsp)

 intarr(1,:)=dtsets(:)%gw_frqim_inzgrid
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_frqim_inzgrid','INT',0)

 intarr(1,:)=dtsets(:)%gw_frqre_inzgrid
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_frqre_inzgrid','INT',0)

 intarr(1,:)=dtsets(:)%gw_frqre_tangrid
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_frqre_tangrid','INT',0)


 intarr(1,:)=dtsets(:)%gw_invalid_freq
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_invalid_freq','INT',0)

 intarr(1,:)=dtsets(:)%gw_qprange
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_qprange','INT',0)

 intarr(1,:)=dtsets(:)%gw_nqlwl
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_nqlwl','INT',0)

 intarr(1,:)=dtsets(:)%gw_nstep
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_nstep','INT',0)

!gw_qlwl
 narr=3*dtsets(1)%gw_nqlwl ! default size for all datasets
 do idtset=0,ndtset_alloc       ! specific size for each dataset
   if(idtset/=0)then
     narrm(idtset)=3*dtsets(idtset)%gw_nqlwl
     if (narrm(idtset)>0)then
       dprarr(1:narrm(idtset),idtset)=&
&       reshape(dtsets(idtset)%gw_qlwl(1:3,1:dtsets(idtset)%gw_nqlwl),(/ narrm(idtset) /) )
     end if
   else
     narrm(idtset)=3*mxvals%gw_nqlwl
     if (narrm(idtset)>0)then
       dprarr(1:narrm(idtset),idtset)=zero
       dprarr(1:3,idtset)=(/0.00001_dp, 0.00002_dp, 0.00003_dp/)
     end if
   end if
 end do

 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,narrm,ncid,ndtset_alloc,'gw_qlwl','DPR',multivals%gw_nqlwl)

 intarr(1,:)=dtsets(:)%gw_sctype
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_sctype','INT',0)

 intarr(1,:)=dtsets(:)%gw_sigxcore
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_sigxcore','INT',0)

 intarr(1,:)=dtsets(:)%gw_icutcoul
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_icutcoul','INT',0)

 dprarr(1,:)=dtsets(:)%gw_toldfeig
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'gw_toldfeig','ENE',0)

 intarr(1,:)=dtsets(:)%hmcsst
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'hmcsst','INT',0)

 intarr(1,:)=dtsets(:)%hmctt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'hmctt','INT',0)

!Special treatment of the default values for the hybrid functional parameters.
 do ii=1,4
   if(ii==1)dprarr(1,:)=dtsets(:)%hyb_mixing
   if(ii==2)dprarr(1,:)=dtsets(:)%hyb_mixing_sr
   if(ii==3)dprarr(1,:)=dtsets(:)%hyb_range_dft
   if(ii==4)dprarr(1,:)=dtsets(:)%hyb_range_fock
   defo=1
   do idtset=1,ndtset_alloc
     if(dprarr(1,idtset)<-tol8 .and. abs(dprarr(1,idtset)+999.0_dp)>tol8)defo=0
   end do
   if(defo==0)then
     do idtset=1,ndtset_alloc
!      Change the sign of user defined input value
       if(dprarr(1,idtset)<-tol8 .and. abs(dprarr(1,idtset)+999.0_dp)>tol8)then
         dprarr(1,idtset)=abs(dprarr(1,idtset))
       end if
     end do
     if(ii==1)str_hyb='hyb_mixing'
     if(ii==2)str_hyb='hyb_mixing_sr'
     if(ii==3)str_hyb='hyb_range_dft'
     if(ii==4)str_hyb='hyb_range_fock'
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,str_hyb,'DPR',0)
   end if
 end do

!###########################################################
!## Deallocation for generic arrays, and for n-z variables

 ABI_DEALLOCATE(dprarr)
 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(narrm)
 ABI_DEALLOCATE(nimagem)
 ABI_DEALLOCATE(dprarr_images)
 ABI_DEALLOCATE(prtimg)

end subroutine outvar_a_h
!!***

end module m_outvar_a_h
!!***
