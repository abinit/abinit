!{\src2tex{textfont=tt}}
!!****f* ABINIT/invars1m
!! NAME
!! invars1m
!!
!! FUNCTION
!! Initialisation phase : prepare the main input subroutine call by
!! reading all the NO MULTI variables, as well as the dimensions
!! needed for allocating the input arrays in abinit.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (XG, MKV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  iout=unit number of output file
!!  lenstr=actual length of string
!!  msym=default maximal number of symmetries
!!  mxnatom=maximal value of input natom for all the datasets
!!  mxnimage=maximal value of input nimage for all the datasets
!!  ndtset= number of datasets to be read; if 0, no multi-dataset mode
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!               one data set.
!!  npsp= number of pseudopotential files
!!  string*(*)=string of characters containing all input variables and data
!!  zionpsp(npsp)= valence charge over all psps
!!
!! OUTPUT
!!  dmatpuflag=flag controlling the use of an initial density matrix in PAW+U (max. value over datasets)
!!  mband_upper_(0:ndtset_alloc)=list of mband_upper values
!!  mxga_n_rules=maximal value of input ga_n_rules for all the datasets
!!  mxgw_nqlwl=maximal value of input gw_nqlwl for all the datasets
!!  mxlpawu=maximal value of input lpawu for all the datasets
!!  mxmband_upper=maximal value of input nband for all the datasets
!!  mxnatpawu=maximal value of number of atoms on which +U is applied for all the datasets
!!  mxnatsph=maximal value of input natsph for all the datasets
!!  mxnatsph_extra=maximal value of input natsph_extra for all the datasets
!!  mxnatvshift=maximal value of input natvshift for all the datasets
!!  mxnconeq=maximal value of input nconeq for all the datasets
!!  mxnkptgw=maximal value of input nkptgw for all the datasets
!!  mxnkpthf=maximal value of input nkpthf for all the datasets
!!  mxnkpt=maximal value of input nkpt for all the datasets
!!  mxnnos=maximal value of input nnos for all the datasets
!!  mxnqptdm=maximal value of input nqptdm for all the datasets
!!  mxnspinor=maximal value of input nspinor for all the datasets
!!  mxnsppol=maximal value of input nsppol for all the datasets
!!  mxnsym=maximum number of symmetries
!!  mxntypat=maximum number of types of atoms
!!  mxnzchempot=maximal value of input nzchempot for all the datasets
!!
!! SIDE EFFECTS
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables,
!!   some of which are initialized here (see invars1.f for more details on the
!!   initialized records)
!!
!! PARENTS
!!      m_ab7_invars_f90
!!
!! CHILDREN
!!      indefo1,invars1
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine invars1m(dmatpuflag,dtsets,iout,lenstr,mband_upper_,&
& msym,mxga_n_rules,mxgw_nqlwl,mxlpawu,mxmband_upper,mxnatom,&
& mxnatpawu,mxnatsph,mxnatsph_extra,mxnatvshift,mxnconeq,&
& mxnimage,mxn_efmas_dirs,mxnkpt,mxnkptgw,mxnkpthf,mxnnos,mxnqptdm,mxnspinor, &
& mxnsppol,mxnsym,mxntypat,mxnimfrqs,mxnfreqsp,mxnzchempot,&
& mxn_projection_frequencies,ndtset,ndtset_alloc,string,npsp,zionpsp)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'invars1m'
 use interfaces_57_iovars, except_this_one => invars1m
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,lenstr,msym,mxnatom,mxnimage,ndtset,ndtset_alloc,npsp
 integer,intent(out) :: dmatpuflag,mxga_n_rules,mxgw_nqlwl,mxlpawu,mxmband_upper,mxnatpawu
 integer,intent(out) :: mxnatsph, mxnatsph_extra
 integer,intent(out) :: mxnatvshift,mxnconeq,mxn_efmas_dirs,mxnkpt,mxnkptgw,mxnkpthf,mxnnos
 integer,intent(out) :: mxnqptdm,mxnspinor,mxnsppol,mxnsym,mxntypat
 integer,intent(out) :: mxnimfrqs,mxnfreqsp,mxnzchempot,mxn_projection_frequencies
 character(len=*),intent(inout) :: string
!arrays
 integer,intent(out) :: mband_upper_(0:ndtset_alloc)
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
 real(dp),intent(in) :: zionpsp(npsp)

!Local variables-------------------------------
!scalars
 integer :: idtset,ii,jdtset,lpawu,mband_upper,iatom,nat,nsp
!arrays
 integer,allocatable :: symafm_(:,:),symrel_(:,:,:,:)
 integer,allocatable :: symafm(:),symrel(:,:,:)
 real(dp),allocatable :: tnons_(:,:,:),tnons(:,:)

!******************************************************************

!Here, allocation of the arrays that depend on msym.
 ABI_ALLOCATE(symrel_,(3,3,msym,0:ndtset_alloc))
 ABI_ALLOCATE(symafm_,(msym,0:ndtset_alloc))
 ABI_ALLOCATE(tnons_,(3,msym,0:ndtset_alloc))
 ABI_ALLOCATE(symafm,(msym))
 ABI_ALLOCATE(symrel,(3,3,msym))
 ABI_ALLOCATE(tnons,(3,msym))

!Set up default values (note that the default acell, amu
!mkmem, mkmem1,mkqmem, and nkpt must be overcome

 do idtset=0,ndtset_alloc
   call indefo1(dtsets(idtset))
 end do

!natom and nimage are already initialized in invars0
 dtsets(0)%natom=-1
 dtsets(0)%nimage=1

!Initialization for parallelization data has changed
!these lines aim to keep old original default values
 dtsets(0)%npimage=1
 dtsets(0)%npkpt=1
 dtsets(0)%npspinor=1
 dtsets(0)%npfft=1
 dtsets(0)%npband=1
 dtsets(0)%bandpp=1

 symafm_(:,0)=1
 symrel_(:,:,:,0)=0
 symrel_(1,1,:,0)=1 ; symrel_(2,2,:,0)=1 ; symrel_(3,3,:,0)=1
 tnons_(:,:,0)=0.0_dp

!Loop on datasets
 do idtset=1,ndtset_alloc
   !write(std_out,'(2a,i0)') ch10,' invars1m : enter jdtset= ',jdtset
   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0

!  Input default values
   dtsets(idtset)%bravais(:)=0
   symafm(:)=symafm_(:,0)
   symrel(:,:,:)=symrel_(:,:,:,0)
   tnons(:,:)=tnons_(:,:,0)

   call invars1(dtsets(idtset)%bravais,dtsets(idtset),iout,jdtset,lenstr,&
&   mband_upper,msym,npsp,string,symafm,symrel,tnons,zionpsp)

   mband_upper_ (idtset)=mband_upper
   symafm_(:,idtset)=symafm(:)
   symrel_(:,:,:,idtset)=symrel(:,:,:)
   tnons_(:,:,idtset)=tnons(:,:)
 end do

 mxmband_upper =maxval(mband_upper_ (1:ndtset_alloc))

 dmatpuflag=0;mxnatpawu=0;mxlpawu=0
 mxnatsph=dtsets(1)%natsph
 mxnatsph_extra=dtsets(1)%natsph_extra
 mxnatvshift=dtsets(1)%natvshift
 mxnconeq=dtsets(1)%nconeq
 mxn_efmas_dirs=0
 mxga_n_rules = dtsets(1)%ga_n_rules
 mxgw_nqlwl = dtsets(1)%gw_nqlwl
 mxnimfrqs = 0
 mxnfreqsp = 0
 mxn_projection_frequencies=0
 mxnkpt  =dtsets(1)%nkpt
 mxnkptgw=dtsets(1)%nkptgw
 mxnkpthf=dtsets(1)%nkpthf
 mxnnos  =dtsets(1)%nnos
 mxnqptdm=dtsets(1)%nqptdm
 mxnspinor=dtsets(1)%nspinor
 mxnsppol=dtsets(1)%nsppol
 mxntypat=dtsets(1)%ntypat
 mxnzchempot=dtsets(1)%nzchempot

!Get MAX dimension over datasets
 do ii=1,ndtset_alloc
   mxnatsph=max(dtsets(ii)%natsph,mxnatsph)
   mxnatsph_extra=max(dtsets(ii)%natsph_extra,mxnatsph_extra)
   mxnconeq=max(dtsets(ii)%nconeq,mxnconeq)
   mxn_efmas_dirs=max(dtsets(ii)%efmas_n_dirs,mxn_efmas_dirs)
   mxga_n_rules = max(dtsets(ii)%ga_n_rules,mxga_n_rules)
   mxgw_nqlwl = max(dtsets(ii)%gw_nqlwl,mxgw_nqlwl)
   mxnimfrqs = max(dtsets(ii)%cd_customnimfrqs,mxnimfrqs)
   mxnfreqsp = max(dtsets(ii)%gw_customnfreqsp,mxnfreqsp)
   mxn_projection_frequencies = max(dtsets(ii)%gwls_n_proj_freq,mxn_projection_frequencies)
   mxnkpt  =max(dtsets(ii)%nkpt,mxnkpt)
   mxnkptgw=max(dtsets(ii)%nkptgw,mxnkptgw)
   mxnkpthf=max(dtsets(ii)%nkpthf,mxnkpthf)
   mxnnos  =max(dtsets(ii)%nnos,mxnnos)
   mxnqptdm=max(dtsets(ii)%nqptdm,mxnqptdm)
   mxnspinor=max(dtsets(ii)%nspinor,mxnspinor)
   mxnsppol=max(dtsets(ii)%nsppol,mxnsppol)
   mxntypat=max(dtsets(ii)%ntypat,mxntypat)
   mxnzchempot=max(dtsets(ii)%nzchempot,mxnzchempot)
   if (dtsets(ii)%usepawu>0) then
     if (dtsets(ii)%usedmatpu/=0) dmatpuflag=1
     lpawu=maxval(dtsets(ii)%lpawu(:))
     mxlpawu=max(lpawu,mxlpawu)
     !dtsets(ii)%natpawu=count(dtsets(ii)%lpawu(dtsets(ii)%typat((/(i1,i1=1,dtsets(ii)%natom)/)))/=-1)
     ! Old fashon way that should do fine
     dtsets(ii)%natpawu = 0
     do iatom=1, dtsets(ii)%natom
       if (dtsets(ii)%lpawu(dtsets(ii)%typat(iatom)) /= -1 ) dtsets(ii)%natpawu = dtsets(ii)%natpawu + 1
     end do
     mxnatpawu=max(dtsets(ii)%natpawu,mxnatpawu)
     if (dtsets(ii)%macro_uj/=0) dtsets(ii)%natvshift=lpawu*2+1
   end if
   mxnatvshift=max(dtsets(ii)%natvshift,mxnatvshift)
 end do

!mxnsym=maxval(dtsets(1:ndtset_alloc)%nsym) ! This might not work properly with HP compiler
 mxnsym=dtsets(1)%nsym
 do idtset=1,ndtset_alloc
   mxnsym=max(dtsets(idtset)%nsym,mxnsym)
 end do

 do idtset=0,ndtset_alloc
   ABI_ALLOCATE(dtsets(idtset)%atvshift,(mxnatvshift,mxnsppol,mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%bs_loband,(mxnsppol))
   ABI_ALLOCATE(dtsets(idtset)%bdgw,(2,mxnkptgw,mxnsppol))
   ABI_ALLOCATE(dtsets(idtset)%cd_imfrqs,(mxnimfrqs))
   ABI_ALLOCATE(dtsets(idtset)%chempot,(3,mxnzchempot,mxntypat))
   nsp=max(mxnsppol,mxnspinor);nat=mxnatpawu*dmatpuflag
   ABI_ALLOCATE(dtsets(idtset)%dmatpawu,(2*mxlpawu+1,2*mxlpawu+1,nsp,nat,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%efmas_bands,(2,mxnkpt))
   ABI_ALLOCATE(dtsets(idtset)%efmas_dirs,(3,mxn_efmas_dirs))
   ABI_ALLOCATE(dtsets(idtset)%gw_freqsp,(mxnfreqsp))
   ABI_ALLOCATE(dtsets(idtset)%gwls_list_proj_freq,(mxn_projection_frequencies))
   ABI_ALLOCATE(dtsets(idtset)%gw_qlwl,(3,mxgw_nqlwl))
   ABI_ALLOCATE(dtsets(idtset)%kpt,(3,mxnkpt))
   ABI_ALLOCATE(dtsets(idtset)%kptgw,(3,mxnkptgw))
   ABI_ALLOCATE(dtsets(idtset)%kptns,(3,mxnkpt))
   ABI_ALLOCATE(dtsets(idtset)%kptns_hf,(3,mxnkpthf))
   ABI_ALLOCATE(dtsets(idtset)%iatsph,(mxnatsph))
   ABI_ALLOCATE(dtsets(idtset)%istwfk,(mxnkpt))
   ABI_ALLOCATE(dtsets(idtset)%nband,(mxnkpt*mxnsppol))
   ABI_ALLOCATE(dtsets(idtset)%occ_orig,(mxmband_upper*mxnkpt*mxnsppol))
   ABI_ALLOCATE(dtsets(idtset)%qmass,(mxnnos))
   ABI_ALLOCATE(dtsets(idtset)%qptdm,(3,mxnqptdm))
   ABI_ALLOCATE(dtsets(idtset)%symafm,(mxnsym))
   ABI_ALLOCATE(dtsets(idtset)%symrel,(3,3,mxnsym))
   ABI_ALLOCATE(dtsets(idtset)%tnons,(3,mxnsym))
   ABI_ALLOCATE(dtsets(idtset)%wtatcon,(3,mxnatom,mxnconeq))
   ABI_ALLOCATE(dtsets(idtset)%wtk,(mxnkpt))
   ABI_ALLOCATE(dtsets(idtset)%xredsph_extra,(3,mxnatsph_extra))
   dtsets(idtset)%symrel(:,:,:)=symrel_(:,:,1:mxnsym,idtset)
   dtsets(idtset)%symafm(:)    =symafm_(1:mxnsym,idtset)
   dtsets(idtset)%tnons (:,:)  =tnons_ (:,1:mxnsym,idtset)
 end do

 ABI_DEALLOCATE(symafm_)
 ABI_DEALLOCATE(symrel_)
 ABI_DEALLOCATE(tnons_)
 ABI_DEALLOCATE(symafm)
 ABI_DEALLOCATE(symrel)
 ABI_DEALLOCATE(tnons)

end subroutine invars1m
!!***
