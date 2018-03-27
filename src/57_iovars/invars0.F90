!{\src2tex{textfont=tt}}
!!****f* ABINIT/invars0
!! NAME
!! invars0
!!
!! FUNCTION
!! Initialisation phase : prepare the main input subroutine call by
!! reading most of the NO MULTI variables, as well as natom, nimage, and ntypat,
!! needed for allocating some input arrays in abinit, and also useri
!! and userr. The variable usewvl is also read here for later reading
!! of input path for the atomic orbital file (if required).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  lenstr=actual length of string
!!  ndtset= number of datasets to be read; if 0, no multi-dataset mode
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!               one data set.
!!  string*(*)=string of characters containing all input variables and data
!!
!! OUTPUT
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables,
!!   some of which are initialized here :
!!   cpus,jdtset,natom,nimage,npsp,ntypat,useri*,userr*
!!  istatr=repetition rate for status file
!!  istatshft=shift of the repetition rate for status file
!!  msym=maximal value of input msym for all the datasets
!!  mxnatom=maximal value of input natom for all the datasets
!!  mxnimage=maximal value of input nimage for all the datasets
!!  mxntypat=maximal value of input ntypat for all the datasets
!!  npsp=number of pseudopotentials
!!
!! PARENTS
!!      m_ab7_invars_f90
!!
!! CHILDREN
!!      get_ndevice,intagm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine invars0(dtsets,istatr,istatshft,lenstr,&
& msym,mxnatom,mxnimage,mxntypat,ndtset,ndtset_alloc,npsp,papiopt,timopt,string)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
#if defined HAVE_GPU_CUDA
 use m_initcuda, only : Get_ndevice
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'invars0'
 use interfaces_32_util
 use interfaces_42_parser
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lenstr,ndtset,ndtset_alloc
 integer,intent(out) :: istatr,istatshft,msym,mxnatom,mxnimage,mxntypat,npsp,papiopt
 integer,intent(inout) :: timopt
 character(len=*),intent(in) :: string
!arrays
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc) !vz_i

!Local variables-------------------------------
!scalars
 integer :: i1,i2,idtset,ii,jdtset,marr,multiplicity,tjdtset,tread,treadh,treadm
 integer :: treads,use_gpu_cuda
 real(dp) :: cpus
 character(len=500) :: message
!arrays
 integer :: supercell_latt(3,3)
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:)

!******************************************************************

!Set ii to avoid warning of uninitialised variable
 ii = 0

 marr=max(9,ndtset_alloc,2)
 ABI_ALLOCATE(dprarr,(marr))
 ABI_ALLOCATE(intarr,(marr))

!Set up jdtset
 if(ndtset/=0)then

!  Default values
   dtsets(0)%jdtset = -1 ! unused value
   dtsets(1:ndtset_alloc)%jdtset=(/ (ii,ii=1,ndtset_alloc) /)

!  Read explicitly the jdtset array
   call intagm(dprarr,intarr,0,marr,ndtset,string(1:lenstr),'jdtset',tjdtset,'INT')
   if(tjdtset==1) dtsets(1:ndtset)%jdtset=intarr(1:ndtset)

!  Read the udtset array
   call intagm(dprarr,intarr,0,marr,2,string(1:lenstr),'udtset',tread,'INT')

!  jdtset and udtset cannot be defined together
   if(tjdtset==1 .and. tread==1)then
     write(message, '(3a)' )&
&     'jdtset and udtset cannot be defined both in the input file.',ch10,&
&     'Action: remove one of them from your input file.'
     MSG_ERROR(message)
   end if

!  Check values of udtset
   if(tread==1)then
     if(intarr(1)<1 .or. intarr(1)>999)then
       write(message, '(a,i0,3a)' )&
&       'udtset(1) must be between 1 and 999, but it is ',intarr(1),'.',ch10,&
&       'Action: change the value of udtset(1) in your input file.'
       MSG_ERROR(message)
     end if
     if(intarr(2)<1 .or. intarr(2)>9)then
       write(message, '(a,i0,3a)' )&
&       'udtset(2) must be between 1 and 9, but it is ',intarr(2),'.',ch10,&
&       'Action: change the value of udtset(2) in your input file.'
       MSG_ERROR(message)
     end if
     if(intarr(1)*intarr(2) /= ndtset)then
       write(message, '(3a,i0,3a,i0,a,i0,3a,i0,3a)' )&
&       'udtset(1)*udtset(2) must be equal to ndtset,',ch10,&
&       'but it is observed that udtset(1) = ',intarr(1),',',ch10,&
&       'and udtset(2) = ',intarr(2),' so that their product is ',intarr(1)*intarr(2),',',ch10,&
&       'while ndtset is ',ndtset,'.',ch10,&
&       'Action: change udtset or ndtset in your input file.'
       MSG_ERROR(message)
     end if
     idtset=0
     do i1=1,intarr(1)
       do i2=1,intarr(2)
         idtset=idtset+1
         dtsets(idtset)%jdtset=i1*10+i2
       end do
     end do
   end if

!  Final check on the jdtset values
   do idtset=1,ndtset
     if(dtsets(idtset)%jdtset<1 .or. dtsets(idtset)%jdtset>9999)then
       write(message, '(3a,i0,a,i0,a,a)' )&
&       'The components of jdtset must be between 1 and 9999.',ch10,&
&       'However, the input value of the component ',idtset,' of jdtset is ',dtsets(idtset)%jdtset,ch10,&
&       'Action : correct jdtset in your input file.'
       MSG_ERROR(message)
     end if
   end do

 else
   dtsets(1)%jdtset=0
 end if

 papiopt = 0
 call intagm(dprarr,intarr,0,1,1,string(1:lenstr),'papiopt',tread,'INT')
 if(tread==1) papiopt=intarr(1)

!Read timopt and pass it to timab
 call intagm(dprarr,intarr,0,1,1,string(1:lenstr),'timopt',tread,'INT')
 if(tread==1) timopt=intarr(1)

 istatr=0
 dtsets(0)%istatr=istatr
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),'istatr',tread,'INT')
 if(tread==1) istatr=intarr(1)
 dtsets(1:)%istatr=istatr

 istatshft=1
 dtsets(0)%istatshft=istatshft
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),'istatshft',tread,'INT')
 if(tread==1) istatshft=intarr(1)
 dtsets(1:)%istatshft=istatshft

 cpus=zero
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),'cpus ',treads,'DPR')
 if(treads==1) cpus=dprarr(1)
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),'cpum ',treadm,'DPR')
 if(treadm==1) cpus=dprarr(1)*60.0_dp
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),'cpuh ',treadh,'DPR')

 if(treadh==1) cpus=dprarr(1)*3600.0_dp
 if(treads+treadm+treadh>1)then
   write(message, '(5a)' )&
&   'More than one input variable is used to defined the CPU time limit.',ch10,&
&   'This is not allowed.',ch10,&
&   'Action: in the input file, suppress either cpus, cpum or cpuh.'
   MSG_ERROR(message)
 end if
 dtsets(:)%cpus=cpus

!Default for natom, nimage, ntypat, useri and userr
 dtsets(:)%natom=1
 dtsets(:)%nimage=1
 dtsets(:)%ntypat=1 ; dtsets(0)%ntypat=0    ! Will always echo ntypat
 dtsets(:)%macro_uj=0
 dtsets(:)%maxnsym=384
 dtsets(:)%useria=0
 dtsets(:)%userib=0
 dtsets(:)%useric=0
 dtsets(:)%userid=0
 dtsets(:)%userie=0
 dtsets(:)%userra=zero
 dtsets(:)%userrb=zero
 dtsets(:)%userrc=zero
 dtsets(:)%userrd=zero
 dtsets(:)%userre=zero
 dtsets(:)%usewvl = 0
 dtsets(:)%plowan_compute=0

!Loop on datasets, to find natom and mxnatom, as well as useri and userr
 do idtset=1,ndtset_alloc
   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0

! proposal: supercell generation in input string before it is read in
! call expand_supercell_input(jdtset, lenstr, string)
!  find supercell, else exit
!  determinant = ncells
!  copy rprim,    acell,    xred,    xcart,    xangst,    vel,    typat,   to
!       rprim_uc, acell_uc, xred_uc, xcart_uc, xangst_uc, vel_uc, typat_uc
!     NB: also rprim and angdeg need to be updated in non diagonal case!!!
!  generate supercell info for each of these copying out with translation vectors etc...
!  set chkprim to 0
!  done!

!  Generate the supercell if supercell_latt is specified and update string
   supercell_latt(:,:) = 0
   do ii=1,3
     supercell_latt(ii,ii) = 1
   end do
   call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),"supercell_latt",tread,'INT')
   if (tread==1) dtsets(idtset)%supercell_latt(:,:)=reshape(intarr(1:marr),(/3,3/))
   !This test should be update if in the future we allow non-diagonal supercell
   if (any(supercell_latt(:,:) < zero) .or. ( supercell_latt(1,1) < tol10 .or.&
&          supercell_latt(2,2) <tol10 .or. supercell_latt(3,3) < tol10 )) then
     write(message, '(5a)' )&
&     'supercell_latt must have positif parameters and diagonal part',ch10,&
&     'This is not allowed.  ',ch10,&
&     'Action : modify supercell_lat in the input file.'
     MSG_ERROR(message)     
   end if
!  Compute the multiplicity of the supercell   
   call mati3det(dtsets(idtset)%supercell_latt,multiplicity)

!  Read natom from string
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natom',tread,'INT')
!  Might initialize natom from XYZ file
   if(tread==0)then
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'_natom',tread,'INT')
   end if

   if(tread==1)then
     dtsets(idtset)%natom=intarr(1)
   else
     write(message, '(a,i0,2a)' )&
&     'Input natom must be defined, but was absent for dataset ',jdtset,ch10,&
&     'Action: check the input file.'
     MSG_ERROR(message)
   end if
!  Check that natom is greater than 0
   if (dtsets(idtset)%natom<=0) then
     write(message, '(a,i0,2a,i0,3a)' )&
&     'Input natom must be > 0, but was ',dtsets(idtset)%natom,ch10,&
&     'for dataset ',jdtset,'. This is not allowed.',ch10,&
&     'Action: check the input file.'
     MSG_ERROR(message)
   end if

   if(multiplicity > 1)then
     dtsets(idtset)%natom = dtsets(idtset)%natom * multiplicity
   end if
   
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nimage',tread,'INT')
   if(tread==1) dtsets(idtset)%nimage=intarr(1)

!  Check that nimage is greater than 0
   if (dtsets(idtset)%nimage<=0) then
     write(message, '(a,i0,4a)' )&
&     'nimage must be > 0, but was ',dtsets(idtset)%nimage,ch10,&
&     'This is not allowed.',ch10,&
&     'Action: check the input file.'
     MSG_ERROR(message)
   end if

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ntypat',tread,'INT')
   if(tread==1)dtsets(idtset)%ntypat=intarr(1)
!  Check that ntypat is greater than 0
   if (dtsets(idtset)%ntypat<=0) then
     write(message, '(a,i0,2a,i0,3a)' )&
&     'Input ntypat must be > 0, but was ',dtsets(idtset)%ntypat,ch10,&
&     'for dataset ',jdtset,'. This is not allowed.',ch10,&
&     'Action: check the input file.'
     MSG_ERROR(message)
   end if

!  Read msym from string
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'maxnsym',tread,'INT')
   if(tread==1)dtsets(idtset)%maxnsym=intarr(1)
!  Check that maxnsym is greater than 1
   if (dtsets(idtset)%maxnsym<1) then
     write(message, '(a,i0,2a,i0,3a)' )&
&     'Input maxnsym must be > 1, but was ',dtsets(idtset)%maxnsym,ch10,&
&     'for dataset ',jdtset,'. This is not allowed.',ch10,&
&     'Action: check the input file.'
     MSG_ERROR(message)
   end if

! Read plowan_compute
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'plowan_compute',tread,'INT')
   if(tread==1) dtsets(idtset)%plowan_compute=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'useria',tread,'INT')
   if(tread==1) dtsets(idtset)%useria=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'userib',tread,'INT')
   if(tread==1) dtsets(idtset)%userib=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'useric',tread,'INT')
   if(tread==1) dtsets(idtset)%useric=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'userid',tread,'INT')
   if(tread==1) dtsets(idtset)%userid=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'userie',tread,'INT')
   if(tread==1) dtsets(idtset)%userie=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'userra',tread,'DPR')
   if(tread==1) dtsets(idtset)%userra=dprarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'userrb',tread,'DPR')
   if(tread==1) dtsets(idtset)%userrb=dprarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'userrc',tread,'DPR')
   if(tread==1) dtsets(idtset)%userrc=dprarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'userrd',tread,'DPR')
   if(tread==1) dtsets(idtset)%userrd=dprarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'userre',tread,'DPR')
   if(tread==1) dtsets(idtset)%userre=dprarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'usewvl',tread,'INT')
   if(tread==1) dtsets(idtset)%usewvl=intarr(1)

 end do

!mxnatom =maxval(dtsets(1:ndtset_alloc)%natom)
!mxntypat =maxval(dtsets(1:ndtset_alloc)%ntypat)
!msym =maxval(dtsets(1:ndtset_alloc)%maxnsym)
!There is a bug in the HP compiler, the following should execute properly
 mxnatom=dtsets(1)%natom ; mxnimage=dtsets(1)%nimage
 mxntypat=dtsets(1)%ntypat ; msym=dtsets(1)%maxnsym
 if(ndtset_alloc>1)then
   do idtset=2,ndtset_alloc
     mxnatom =max(dtsets(idtset)%natom,mxnatom)
     mxnimage=max(dtsets(idtset)%nimage,mxnimage)
     mxntypat=max(dtsets(idtset)%ntypat,mxntypat)
     msym    =max(dtsets(idtset)%maxnsym,msym)
   end do
 end if

 if(mxnimage>1)then
   do idtset=2,ndtset_alloc
     if(mxnatom/=dtsets(idtset)%natom)then
       write(message,'(5a,i0,a,i0,3a,i0,a)')&
&       'When there exist one dataset with more than one image,',ch10,&
&       'the number of atoms in each dataset must be the same.',ch10,&
&       'However, it has been found that for dataset= ',idtset,ch10,&
&       'natom= ',dtsets(idtset)%natom,' differs from the maximum number',ch10,&
&       'of atoms, mxnatom= ',mxnatom,&
&       'Action: check the input variables natom for different datasets.'
       MSG_ERROR(message)
     end if
   end do
 end if

!Set up npsp
 npsp=mxntypat   ! Default value
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),'npsp',tread,'INT')
 if(tread==1)then
   npsp=intarr(1)
 else
   if(ndtset_alloc>1)then
     do idtset=1,ndtset_alloc
       if(dtsets(idtset)%ntypat/=mxntypat)then
         write(message, '(5a,i0,a,i0,2a,i0,2a)' )&
&         'When npsp is not defined, the input variable ntypat must be',ch10,&
&         'the same for all datasets. However, it has been found that for',ch10,&
&         'jdtset: ',dtsets(idtset)%jdtset,', ntypat= ',dtsets(idtset)%ntypat,ch10,&
&         'differs from the maximum value of ntypat= ',mxntypat,ch10,&
&         'Action: check the input variables npsp and ntypat.'
         MSG_ERROR(message)
       end if
     end do
   end if
 end if
 dtsets(0)%npsp = mxntypat   ! Default value
 dtsets(1:ndtset_alloc)%npsp = npsp

!KGB parallelism information (needed at this stage)
 dtsets(:)%paral_kgb=0
 do idtset=1,ndtset_alloc
   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'paral_kgb',tread,'INT')
   if(tread==1)dtsets(idtset)%paral_kgb=intarr(1)

   if (dtsets(idtset)%paral_kgb<0 .or. dtsets(idtset)%paral_kgb>1) then
     write(message,'(a,i0,2a,i0,3a)')&
&     'Input paral_kgb must be 0 or 1, but was ',dtsets(idtset)%paral_kgb,ch10,&
&     'for dataset',jdtset,'. This is not allowed.',ch10,&
&     'Action: check the input file.'
     MSG_ERROR(message)
   end if
 end do

!GPU information
 use_gpu_cuda=0
 dtsets(:)%use_gpu_cuda=0
#if defined HAVE_GPU_CUDA && defined HAVE_GPU_CUDA_DP
 call Get_ndevice(ii)
 if (ii>0) then
   do i1=1,ndtset_alloc
     dtsets(i1)%use_gpu_cuda=-1
   end do
 end if
#endif
 do idtset=1,ndtset_alloc
   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'use_gpu_cuda',tread,'INT')
   if(tread==1)dtsets(idtset)%use_gpu_cuda=intarr(1)
   if (dtsets(idtset)%use_gpu_cuda==1) use_gpu_cuda=1
 end do
 if (use_gpu_cuda==1) then
#if defined HAVE_GPU_CUDA && defined HAVE_GPU_CUDA_DP
   if (ii<=0) then
     write(message,'(3a)')&
&     'Input variables use_gpu_cuda is on',ch10,&
&     'but no available GPU device has been detected !',ch10,&
&     'Action: change the input variable use_gpu_cuda.'
     MSG_ERROR(message)
   end if
#else
   write(message,'(7a)')&
&   'Input variables use_gpu_cuda is on but abinit hasn''t been built',ch10,&
&   'with (double precision) gpu mode enabled !',ch10,&
&   'Action: change the input variable use_gpu_cuda',ch10,&
&   '        or re-compile ABINIT with double-precision Cuda enabled.'
   MSG_ERROR(message)
#endif
 end if

 ABI_DEALLOCATE(dprarr)
 ABI_DEALLOCATE(intarr)

!We allocate the internal array, depending on the computed values.
!WARNING : do not forget to deallocate these arrays in the routine dtset_free
!(should make a separate subroutine for allocating/deallocating these records)
 do idtset=0,ndtset_alloc
   ABI_ALLOCATE(dtsets(idtset)%acell_orig,(3,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%algalch,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%amu_orig,(mxntypat,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%corecs,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%densty,(mxntypat,4))
   ABI_ALLOCATE(dtsets(idtset)%dynimage,(mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%iatfix,(3,mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%f4of2_sla,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%f6of2_sla,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%jpawu,(mxntypat,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%kberry,(3,20))
   ABI_ALLOCATE(dtsets(idtset)%lexexch,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%ldaminushalf,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%lpawu,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%mixalch_orig,(npsp,mxntypat,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%nucdipmom,(3,mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%pimass,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%ptcharge,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%prtatlist,(mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%quadmom,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%ratsph,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%rprim_orig,(3,3,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%rprimd_orig,(3,3,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%so_psp,(npsp))
   ABI_ALLOCATE(dtsets(idtset)%spinat,(3,mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%shiftk,(3,210))
   ABI_ALLOCATE(dtsets(idtset)%typat,(mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%upawu,(mxntypat,mxnimage))
!   if (dtsets(idtset)%plowan_compute>0) then
   ABI_ALLOCATE(dtsets(idtset)%plowan_iatom,(mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%plowan_it,(100*3))
   ABI_ALLOCATE(dtsets(idtset)%plowan_nbl,(mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%plowan_lcalc,(12*mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%plowan_projcalc,(12*mxnatom))
!   endif
   ABI_ALLOCATE(dtsets(idtset)%vel_orig,(3,mxnatom,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%vel_cell_orig,(3,3,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%xred_orig,(3,mxnatom,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%ziontypat,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%znucl,(npsp))
 end do

!DEBUG
!write(std_out,*)' invars0 : nimage, mxnimage = ',dtsets(:)%nimage, mxnimage
!write(std_out,*)' invars0 : natom = ',dtsets(:)%natom
!write(std_out,*)' invars0 : mxnatom = ',mxnatom
!ENDDEBUG

end subroutine invars0
!!***
