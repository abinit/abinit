!!****m* ABINIT/m_invars2
!! NAME
!!  m_invars2
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1999-2020 ABINIT group (XG)
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

module m_invars2

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_nctk
 use m_sort
 use m_dtset
 use libxc_functionals
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_datatypes, only : pspheader_type
 use m_time,      only : timab
 use m_fstrings,  only : sjoin, itoa, ltoa, tolower
 use m_symtk,     only : matr3inv
 use m_parser,    only : intagm, intagm_img
 use m_geometry,   only : mkrdim, metric
 use m_gsphere,    only : setshells
 use m_xcdata,    only : get_auxc_ixc, get_xclevel
 use m_inkpts,    only : inkpts
 use m_ingeo,     only : invacuum

 implicit none

 private
!!***

 public :: invars2m
 public :: invars2
!!***

contains
!!***

!!****f* ABINIT/invars2m
!! NAME
!! invars2m
!!
!! FUNCTION
!! Initialisation phase - main input routine.
!! Big loop on the datasets :
!! - for each of the datasets, write one line about the crystallographic data
!! - call invars2, that read the eventual single dataset input values;
!! - compute mgfft,mpw,nfft,... for this data set;
!! - compute quantities for the susceptibility computation
!! - compute the memory needs for this data set.
!!
!! INPUTS
!!  bravais_(11,0:ndtset_alloc)=characteristics of Bravais lattice
!!  iout=unit number of output file
!!  lenstr=actual length of string
!!  mband_upper_(0:ndtset_alloc)=list of mband_upper values
!!  msym=default maximal number of symmetries
!!  ndtset= number of datasets to be read; if 0, no multi-dataset mode
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!      one data set.
!!  npsp=number of pseudopotentials
!!  pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!  string*(*)=character string containing all the input data.
!!   Initialized previously in instrng.
!! comm=MPI communicator
!!
!! OUTPUT
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables,
!!   some of which are initialized here, while other were already
!!   initialized previously.
!!
!! NOTES
!! The outputs of this routine are the values of input variables,
!! their default value is stored at the index 0 of the last dimension
!! of their multi-dataset representation.
!!
!! PARENTS
!!      m_common
!!
!! CHILDREN
!!      dtset%chkneu,get_auxc_ixc,get_xclevel,inkpts,intagm,intagm_img,invacuum
!!      libxc_functionals_end,libxc_functionals_get_hybridparams
!!      libxc_functionals_init,sort_int,timab,wrtout
!!
!! SOURCE

subroutine invars2m(dtsets,iout,lenstr,mband_upper_,msym,ndtset,ndtset_alloc,npsp,pspheads,string, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,lenstr,msym,ndtset,ndtset_alloc,npsp, comm
 character(len=*),intent(in) :: string
!arrays
 integer,intent(in) :: mband_upper_(0:ndtset_alloc)
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
 type(pspheader_type),intent(in) :: pspheads(npsp)

!Local variables -------------------------------
!scalars
 integer :: idtset,jdtset,mband_upper,nsheps,nshsigx,nshwfn,usepaw
 real(dp) :: ucvol
!arrays
 integer :: bravais(11)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: zionpsp(:)

!*************************************************************************

 do idtset=1,ndtset_alloc
   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0
   bravais(:)=dtsets(idtset)%bravais(:)
   mband_upper  =mband_upper_(idtset)
   usepaw=dtsets(idtset)%usepaw
   ! Allocate arrays
   ABI_ALLOCATE(zionpsp,(npsp))
   zionpsp(:)=pspheads(1:npsp)%zionpsp

   call mkrdim(dtsets(idtset)%acell_orig(1:3,1),dtsets(idtset)%rprim_orig(1:3,1:3,1),rprimd)
   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

   ! Here, nearly all the remaining input variables are initialized
   call invars2(bravais,dtsets(idtset),iout,jdtset,lenstr,mband_upper,msym,npsp,string,usepaw,zionpsp,ucvol,comm)

   ABI_DEALLOCATE(zionpsp)

   if (ANY(dtsets(idtset)%optdriver == [RUNL_SCREENING,RUNL_SIGMA,RUNL_BSE])) then
    ! For GW or BSE calculations, we only use (npwwfn|ecutwfn) G-vectors read from the KSS file,
    ! therefore the FFT box for the density should be defined according to ecut=ecutwfn.

     nshwfn=0
     call setshells(dtsets(idtset)%ecutwfn,dtsets(idtset)%npwwfn,nshwfn,&
       dtsets(idtset)%nsym,gmet,gprimd,dtsets(idtset)%symrel,'wfn',ucvol)

     ! MG: Hack to avoid portability problems under gfortran:
     ! getng and getmpw are indeed quite sensitive if ecut is small
     ! and, in the GW tests, mpw and ngfft might depend on the compiler used.
     ! the problem shows up if we use npwwfn instead of ecutwfn, a good reason for removing npwwfn!
     dtsets(idtset)%ecutwfn=dtsets(idtset)%ecutwfn-tol14
     ! MG: This is a kind of a hack, but the problem is ecutwfn that is too much redundant!
     dtsets(idtset)%ecut=dtsets(idtset)%ecutwfn

     ! Close the shell for (W|chi0)
     nsheps=0
     call setshells(dtsets(idtset)%ecuteps,dtsets(idtset)%npweps,nsheps,&
      dtsets(idtset)%nsym,gmet,gprimd,dtsets(idtset)%symrel,'eps',ucvol)

     ! Close the shell for the exchange term.
     nshsigx=0
     call setshells(dtsets(idtset)%ecutsigx,dtsets(idtset)%npwsigx,nshsigx,&
      dtsets(idtset)%nsym,gmet,gprimd,dtsets(idtset)%symrel,'sigx',ucvol)
   end if ! (SIGMA|SCREENING|SCGW|BSE|RDFTM)

 end do

end subroutine invars2m
!!***

!!****f* ABINIT/invars2
!! NAME
!! invars2
!!
!! FUNCTION
!! Initialize variables for the ABINIT code, for one particular dataset, characterized by jdtset.
!! Note: some parameters have already been read in invars0 and invars1,
!! and were used to dimension the arrays needed here.
!!
!! INPUTS
!! bravais(11): bravais(1)=iholohedry
!!              bravais(2)=center
!!              bravais(3:11)=coordinates of rprim in the axes
!!               of the conventional bravais lattice (*2 if center/=0)
!! iout=unit number for echoed output
!! jdtset=number of the dataset looked for
!! lenstr=actual length of the string
!! mband=maximum number of bands for any k-point
!! msym=default maximal number of symmetries
!! npsp=number of pseudopotentials
!! string*(*)=character string containing all the input data.
!!  Initialized previously in instrng.
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!! zionpsp(npsp)=valence charge of each type of atom (coming from the psp files)
!! ucvol: Unit cell volume
!! comm=MPI communicator
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output:
!!  dtset=<type datafiles_type>contains all input variables,
!!   some of which are initialized here, while other were already
!! All rest of arguments given alphabetically from acell (length scales)
!! to wtk (k-point weights), used to control running of the main routine.
!! See abinit_help for definitions of these variables.
!! These values become defined by being read from string,
!! that contains all information from the input file,
!! in a compressed, standardized, format
!! At the input, they already contain a default value.
!!
!! PARENTS
!!      m_invars2
!!
!! CHILDREN
!!      dtset%chkneu,get_auxc_ixc,get_xclevel,inkpts,intagm,intagm_img,invacuum
!!      libxc_functionals_end,libxc_functionals_get_hybridparams
!!      libxc_functionals_init,sort_int,timab,wrtout
!!
!! SOURCE

subroutine invars2(bravais,dtset,iout,jdtset,lenstr,mband,msym,npsp,string,usepaw,zionpsp,ucvol,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,jdtset,lenstr,mband,msym,npsp,usepaw, comm
 character(len=*),intent(in) :: string
 type(dataset_type),intent(inout) :: dtset
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: bravais(11)
 real(dp),intent(in) :: zionpsp(npsp)

!Local variables-------------------------------
!scalars
 integer :: bantot,berryopt,dmatsize,ndim,getocc,narr,nprocs
 integer :: iat,iatom,iband,ii,iimage,ikpt,intimage,ionmov,isppol,ixc_current
 integer :: densfor_pred,ipsp,iscf,isiz,itypat,jj,kptopt,lpawu,marr,natom,nband1,nberry
 integer :: niatcon,nimage,nkpt,nkpthf,npspalch,nqpt,nsp,nspinor,nsppol,nsym,ntypalch,ntypat,ntyppure
 integer :: occopt,occopt_tmp,response,sumnbl,tfband,tnband,tread,tread_alt,tread_dft,tread_fock,tread_key, tread_extrael
 integer :: itol, itol_gen, ds_input, ifreq,ncerr
 real(dp) :: areaxy,charge,fband,kptrlen,nelectjell,sum_spinat
 real(dp) :: rhoavg,zelect,zval
 real(dp) :: toldfe_, tolrff_, toldff_, tolwfr_, tolvrs_
 real(dp) :: tolmxde_, tolmxf_
 character(len=500) :: msg
 character(len=fnlen) :: key_value
!arrays
 integer :: vacuum(3)
 integer,allocatable :: iatcon(:),natcon(:), intarr(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: dmatpawu_tmp(:), dprarr(:)

! *************************************************************************

 call timab(191,1,tsec)

 nprocs = xmpi_comm_size(comm)

 ! Compute the maximum size of arrays intarr and dprarr
 natom=dtset%natom
 nimage=dtset%nimage
 nkpt=dtset%nkpt
 nkpthf=dtset%nkpthf
 npspalch=dtset%npspalch
 nspinor=dtset%nspinor
 nsppol=dtset%nsppol
 ntypat=dtset%ntypat
 ntypalch=dtset%ntypalch
 ntyppure=dtset%ntyppure

 dmatsize=0 ! dmatpu not available for usepawu<0
 if (dtset%usepawu>0.and.dtset%usedmatpu/=0) then
   do iatom=1,natom
     lpawu=dtset%lpawu(dtset%typat(iatom))
     if (lpawu/=-1) dmatsize=dmatsize+nsppol*nspinor*(2*lpawu+1)**2
   end do
 end if

 marr=max(3*natom,&
   nkpt*nsppol*mband,&
   2*dtset%nkptgw*nsppol,&
   dmatsize,&
   3*nkpt,&
   npsp,&
   ntypat,&
   9*msym,&
   60,100,&
   3*dtset%nconeq*natom,&
   nimage,&
   3*dtset%nqptdm,&
   3*dtset%natsph_extra,&
   dtset%natvshift*nsppol*natom,&
   3*dtset%nzchempot*ntypat)

 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

 !----------------------------------------------------------------------------

 !****   Read parameters which set remaining array dimensions ****

 !Note: some parameters have already been read in invars0 and invars1
 !Also, some checking is needed here.

 !Read ngfft(1), ngfft(2), and ngfft(3), then ngfft(7)=fftalg and ngfft(8)=fftcache.
 !Read ngfftdg(1:3)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'ngfft',tread,'INT')
 if(tread==1) dtset%ngfft(1:3)=intarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fftalg',tread,'INT')
 if(tread==1) then
   dtset%ngfft(7)=intarr(1)
   if (usepaw==1) dtset%ngfftdg(7)=intarr(1)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fftcache',tread,'INT')

 if(tread==1) then
   dtset%ngfft(8)=intarr(1)
   if (usepaw==1) dtset%ngfftdg(8)=intarr(1)
 end if

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'ngfftdg',tread,'INT')
 if(tread==1) dtset%ngfftdg(1:3)=intarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'mqgrid',tread,'INT')
 if(tread==1) dtset%mqgrid=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'mqgriddg',tread,'INT')
 if(tread==1) dtset%mqgriddg=intarr(1)

 ! Make special arrangements to check nband: may be a scalar
 ! (for occopt=0, 1 or 3, 4, 5, 6, 7, 8) or may be an array (for occopt=2)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'occopt',tread,'INT')
 if(tread==1) dtset%occopt=intarr(1)
 occopt=dtset%occopt

 ! check for variables related to genetic algorithm. ga_n_rules has been already read

 if (dtset%imgmov==4) then
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ga_algor',tread,'INT')
   if(tread==1) dtset%ga_algor=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ga_fitness',tread,'INT')
   if(tread==1) dtset%ga_fitness=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ga_opt_percent',tread,'ENE')
   if(tread==1) dtset%ga_opt_percent=dprarr(1)

   call intagm(dprarr,intarr,jdtset,marr,dtset%ga_n_rules,string(1:lenstr),'ga_rules',tread,'INT')
   if(tread==1)then
     dtset%ga_rules(1:dtset%ga_n_rules)=intarr(1:dtset%ga_n_rules)
     do ii=1,dtset%ga_n_rules
       if(dtset%ga_rules(ii)<0)then
         write(msg, '(3a)' )&
         'All values for Genetic rules must be greater than 0.',ch10,&
         'Action: check the values of ga_rules.'
         MSG_ERROR(msg)
       end if
     end do
   end if
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwaclowrank',tread,'INT')
 if(tread==1) dtset%gwaclowrank=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwcalctyp',tread,'INT')
 if(tread==1) dtset%gwcalctyp=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwcomp',tread,'INT')
 if(tread==1) dtset%gwcomp=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwencomp',tread,'ENE')
 if(tread==1) dtset%gwencomp=dprarr(1)

 if (usepaw==1) then
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gw_sigxcore',tread,'INT')
   if(tread==1) dtset%gw_sigxcore=intarr(1)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwmem',tread,'INT')
 if(tread==1) dtset%gwmem=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gw_sctype',tread,'INT')
 if(tread==1) dtset%gw_sctype=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gw_nstep',tread,'INT')
 if(tread==1) dtset%gw_nstep=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gw_toldfeig',tread,'ENE')
 if(tread==1) dtset%gw_toldfeig=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwls_stern_kmax',tread,'INT')
 if(tread==1) dtset%gwls_stern_kmax=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwls_npt_gauss_quad',tread,'INT')
 if(tread==1) dtset%gwls_npt_gauss_quad=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwls_diel_model',tread,'INT')
 if(tread==1) dtset%gwls_diel_model=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwls_model_parameter',tread,'ENE')
 if(tread==1) dtset%gwls_model_parameter=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwls_recycle',tread,'INT')
 if(tread==1) dtset%gwls_recycle=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwls_print_debug',tread,'INT')
 if(tread==1) dtset%gwls_print_debug=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwls_nseeds',tread,'INT')
 if(tread==1) dtset%gwls_nseeds=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwls_kmax_complement',tread,'INT')
 if(tread==1) dtset%gwls_kmax_complement=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwls_kmax_poles',tread,'INT')
 if(tread==1) dtset%gwls_kmax_poles=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwls_kmax_analytic',tread,'INT')
 if(tread==1) dtset%gwls_kmax_analytic=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwls_kmax_numeric',tread,'INT')
 if(tread==1) dtset%gwls_kmax_numeric=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwls_exchange',tread,'INT')
 if(tread==1) dtset%gwls_exchange=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwls_correlation',tread,'INT')
 if(tread==1) dtset%gwls_correlation=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwls_band_index',tread,'INT')
 if(tread==1) dtset%gwls_band_index=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gw1rdm',tread,'INT')
 if(tread==1) dtset%gw1rdm=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gw1rdm_energies',tread,'INT')
 if(tread==1) dtset%gw1rdm_energies=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'x1rdm',tread,'INT')
 if(tread==1) dtset%x1rdm=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'chkp_rdm',tread,'INT')
 if(tread==1) dtset%chkp_rdm=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'chkp_rdm_lo',tread,'INT')
 if(tread==1) dtset%chkp_rdm_lo=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwls_first_seed',tread,'INT')
 if(tread==1) then
   dtset%gwls_first_seed=intarr(1)
 else
   dtset%gwls_first_seed=dtset%gwls_band_index
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rhoqpmix',tread,'DPR')
 if(tread==1) dtset%rhoqpmix=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rifcsph',tread,'DPR')
 if(tread==1) dtset%rifcsph=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nfreqim',tread,'INT')
 if(tread==1) dtset%nfreqim=intarr(1)
 if (dtset%cd_customnimfrqs/=0) then
   write(msg, '(3a)' )&
   'cd_customnimfrqs not equal to zero and not equal to nfreqim',ch10,&
   'setting nfreqim = cd_customnimfrqs'
   MSG_WARNING(msg)
   dtset%nfreqim = dtset%cd_customnimfrqs
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'freqim_alpha',tread,'DPR')
 if(tread==1) dtset%freqim_alpha=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'freqremin',tread,'ENE')
 if(tread==1) dtset%freqremin=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'freqremax',tread,'ENE')
 if(tread==1) dtset%freqremax=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nfreqre',tread,'INT')
 if(tread==1) dtset%nfreqre=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nfreqsp',tread,'INT')
 if(tread==1) dtset%nfreqsp=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'freqspmax',tread,'ENE')
 if(tread==1) dtset%freqspmax=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'freqspmin',tread,'ENE')
 if(tread==1) then
   ! If found, set it
   dtset%freqspmin=dprarr(1)
 else
   ! Else give it the value -freqspmax
   dtset%freqspmin=-dtset%freqspmax
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gw_frqim_inzgrid',tread,'INT')
 if(tread==1) dtset%gw_frqim_inzgrid=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gw_frqre_inzgrid',tread,'INT')
 if(tread==1) dtset%gw_frqre_inzgrid=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gw_frqre_tangrid',tread,'INT')
 if(tread==1) dtset%gw_frqre_tangrid=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gw_invalid_freq',tread,'INT')
 if(tread==1) dtset%gw_invalid_freq=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gw_icutcoul',tread,'INT')
 if(tread==1) dtset%gw_icutcoul=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gw_qprange',tread,'INT')
 if(tread==1) dtset%gw_qprange=intarr(1)

 if(dtset%cd_customnimfrqs/=0) then
   call intagm(dprarr,intarr,jdtset,marr,dtset%cd_customnimfrqs,string(1:lenstr),'cd_imfrqs',tread,'ENE')
   if(tread==1) then
     dtset%cd_imfrqs(1:dtset%cd_customnimfrqs)=dprarr(1:dtset%cd_customnimfrqs)
     do ifreq=1,dtset%cd_customnimfrqs-1
       if (dtset%cd_imfrqs(ifreq+1)<dtset%cd_imfrqs(ifreq)) then
         write(msg, '(3a)' )&
         'The frequencies specified in cd_imfrqs must be strictly increasing',ch10,&
         'Action: Correct this in your input file'
         MSG_ERROR(msg)
       end if
     end do
   end if
 end if

 if(dtset%gw_customnfreqsp/=0) then
   call intagm(dprarr,intarr,jdtset,marr,dtset%gw_customnfreqsp,string(1:lenstr),'gw_freqsp',tread,'ENE')
   if(tread==1) then
     dtset%gw_freqsp(1:dtset%gw_customnfreqsp)=dprarr(1:dtset%gw_customnfreqsp)
     do ifreq=1,dtset%gw_customnfreqsp-1
       if (dtset%gw_freqsp(ifreq+1)<dtset%gw_freqsp(ifreq)) then
         write(msg, '(3a)' )&
         'The frequencies specified in gw_freqsp must be strictly increasing',ch10,&
         'Action: Correct this in your input file'
         MSG_ERROR(msg)
       end if
     end do
   end if
 end if

 if(dtset%gwls_n_proj_freq/=0) then
   call intagm(dprarr,intarr,jdtset,marr,dtset%gwls_n_proj_freq,string(1:lenstr),'gwls_list_proj_freq',tread,'ENE')
   if(tread==1) then
     dtset%gwls_list_proj_freq(1:dtset%gwls_n_proj_freq)=dprarr(1:dtset%gwls_n_proj_freq)
   end if
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'cd_full_grid',tread,'INT')
 if(tread==1) dtset%cd_full_grid=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'cd_frqim_method',tread,'INT')
 if(tread==1) dtset%cd_frqim_method=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'cd_halfway_freq',tread,'ENE')
 if(tread==1) dtset%cd_halfway_freq=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'cd_max_freq',tread,'ENE')
 if(tread==1) dtset%cd_max_freq=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'cd_subset_freq',tread,'INT')
 if(tread==1) dtset%cd_subset_freq(1:2)=intarr(1:2)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwrpacorr',tread,'INT')
 if(tread==1) dtset%gwrpacorr=intarr(1)

 ! RESPFN integer input variables (needed here to get the value of response)
 ! Warning: rfddk,rfelfd,rfmagn,rfphon,rfstrs,rfuser,rf2_dkdk and rf2_dkde are also read in invars1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfasr',tread,'INT')
 if(tread==1) dtset%rfasr=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'rfatpol',tread,'INT')
 if(tread==1) dtset%rfatpol(1:2)=intarr(1:2)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'rfdir',tread,'INT')
 if(tread==1) dtset%rfdir(1:3)=intarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfddk',tread,'INT')
 if(tread==1) dtset%rfddk=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfelfd',tread,'INT')
 if(tread==1) dtset%rfelfd=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfmagn',tread,'INT')
 if(tread==1) dtset%rfmagn=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfmeth',tread,'INT')
 if(tread==1) dtset%rfmeth=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfphon',tread,'INT')
 if(tread==1) dtset%rfphon=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfstrs',tread,'INT')
 if(tread==1) dtset%rfstrs=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfuser',tread,'INT')
 if(tread==1) dtset%rfuser=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rf2_dkdk',tread,'INT')
 if(tread==1) dtset%rf2_dkdk=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rf2_dkde',tread,'INT')
 if(tread==1) dtset%rf2_dkde=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'rf2_pert1_dir',tread,'INT')
 if(tread==1) dtset%rf2_pert1_dir(1:3)=intarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'rf2_pert2_dir',tread,'INT')
 if(tread==1) dtset%rf2_pert2_dir(1:3)=intarr(1:3)

 ! Set value of response to 1 and also set rfdir to 1 1 1 if we are doing
 ! response calculation but rfdir was left at default 0 0 0 value.
 ! For rf2_dkdk and rf2_dkde, we do the same for rf2_pert1_dir and rf2_pert2_dir
 response=0
 if(dtset%rfddk/=0 .or. dtset%rf2_dkdk/=0 .or. dtset%rf2_dkde/=0 .or. dtset%rfelfd/=0 .or. &
   dtset%rfphon/=0 .or. dtset%rfstrs/=0 .or. dtset%rfuser/=0 ) then
   response=1
   if( (dtset%rfdir(1) == 0) .and. (dtset%rfdir(2) == 0) .and. (dtset%rfdir(3) == 0) ) dtset%rfdir(1:3) = 1
   if (dtset%rf2_dkdk/=0 .or. dtset%rf2_dkde/=0) then
     if (sum(abs(dtset%rf2_pert1_dir)) == 0) dtset%rf2_pert1_dir(:) = 1
     if (sum(abs(dtset%rf2_pert2_dir)) == 0) dtset%rf2_pert2_dir(:) = 1
   end if
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nonlinear_info',tread,'INT')
 if(tread==1) dtset%nonlinear_info=intarr(1)

 ! NONLINEAR integer input variables (same definition as for rfarr)
 ! Presently, rf?asr, rf?meth,rf?strs and rf?thrd are not used
 ! --Keep the old input variables for backward compatibility
 if(dtset%optdriver==RUNL_NONLINEAR) then
   tread_key=0

   call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'d3e_pert1_atpol',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'rf1atpol',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert1_atpol(1:2)=intarr(1:2)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'d3e_pert1_dir',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'rf1dir',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert1_dir(1:3)=intarr(1:3)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'d3e_pert1_elfd',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rf1elfd',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert1_elfd=intarr(1)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'d3e_pert1_phon',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rf1phon',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert1_phon=intarr(1)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'d3e_pert2_atpol',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'rf2atpol',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert2_atpol(1:2)=intarr(1:2)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'d3e_pert2_dir',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'rf2dir',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert2_dir(1:3)=intarr(1:3)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'d3e_pert2_elfd',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rf2elfd',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert2_elfd=intarr(1)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'d3e_pert2_phon',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rf2phon',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert2_phon=intarr(1)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'d3e_pert3_atpol',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'rf3atpol',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert3_atpol(1:2)=intarr(1:2)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'d3e_pert3_dir',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'rf3dir',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert3_dir(1:3)=intarr(1:3)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'d3e_pert3_elfd',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rf3elfd',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert3_elfd=intarr(1)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'d3e_pert3_phon',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rf3phon',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert3_phon=intarr(1)
   if (tread_alt==1) tread_key=1

   if (tread_key==1) then
     msg='The following input keywords are obsolete:'//ch10//&
     '  rfxatpol, rfxdir, rfxrlfd, rfxphon (with x=1,2,3)'//ch10//&
     'Action: change to the d3e_pertx_*** input parameters!'
     MSG_WARNING(msg)
   end if

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'usepead',tread,'INT')
   if(tread==1) dtset%usepead=intarr(1)
 end if

 response=0
 if(dtset%rfddk/=0 .or. dtset%rfphon/=0 .or. dtset%rfelfd/=0 .or. &
    dtset%rfstrs/=0 .or. dtset%rfuser/=0 .or. &
    dtset%rf2_dkdk/=0 .or. dtset%rf2_dkde/=0 .or. &
    dtset%d3e_pert1_elfd/=0 .or. dtset%d3e_pert1_phon/=0 .or. &
    dtset%d3e_pert2_elfd/=0 .or. dtset%d3e_pert2_phon/=0 .or. &
    dtset%d3e_pert3_elfd/=0 .or. dtset%d3e_pert3_phon/=0 ) response=1

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prepanl',tread,'INT')
 if(tread==1) dtset%prepanl=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prepgkk',tread,'INT')
 if(tread==1) dtset%prepgkk=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'use_nonscf_gkk',tread,'INT')
 if(tread==1) dtset%use_nonscf_gkk=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'use_gemm_nonlop',tread,'INT')
 if(tread==1) dtset%use_gemm_nonlop=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'use_yaml',tread,'INT')
 if(tread==1) dtset%use_yaml=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'brav',tread,'INT')
 if(tread==1) dtset%brav=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'boxcutmin',tread,'DPR')
 if(tread==1) dtset%boxcutmin=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'charge',tread,'DPR')
 if(tread==1) dtset%charge=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dosdeltae',tread,'ENE')
 if(tread==1) dtset%dosdeltae=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dtion',tread,'DPR')
 if(tread==1) dtset%dtion=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ecut',tread,'ENE')
 if(tread==1) dtset%ecut=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,size(dtset%einterp),string(1:lenstr),'einterp',tread,'DPR')
 if(tread==1) dtset%einterp=dprarr(1:size(dtset%einterp))

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'elph2_imagden',tread,'ENE')
 if(tread==1) dtset%elph2_imagden=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'esmear',tread,'ENE')
 if(tread==1) dtset%esmear=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fermie_nest',tread,'ENE')
 if(tread==1) dtset%fermie_nest=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dfpt_sciss',tread,'ENE')
 if(tread==1) dtset%dfpt_sciss=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'tmesh',tread,'DPR')
 if(tread==1) dtset%tmesh=dprarr(1:3)
 ABI_CHECK(all(dtset%tmesh >= zero), sjoin("Invalid tmesh containg T < 0:", ltoa(dtset%tmesh)))

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'tsmear',tread,'ENE')
 if(tread==1) dtset%tsmear=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vis',tread,'DPR')
 if(tread==1) dtset%vis=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ecutsm',tread,'ENE')
 if(tread==1) dtset%ecutsm=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'exchmix',tread,'DPR')
 if(tread==1) dtset%exchmix=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dilatmx',tread,'DPR')
 if(tread==1) dtset%dilatmx=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fxcartfactor',tread,'DPR')
 if(tread==1) dtset%fxcartfactor=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'strfact',tread,'DPR')
 if(tread==1) dtset%strfact=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'effmass_free',tread,'DPR')
 if(tread==1) dtset%effmass_free=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'mdtemp',tread,'DPR')
 if(tread==1) dtset%mdtemp(1:2)=dprarr(1:2)

!LONG WAVE integer input variables
!FIXME
! if(dtset%optdriver==RUNL_LONGWAVE) then
   tread_key=0

   call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'d3e_pert1_atpol',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'rf1atpol',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert1_atpol(1:2)=intarr(1:2)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'d3e_pert1_dir',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'rf1dir',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert1_dir(1:3)=intarr(1:3)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'d3e_pert1_elfd',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rf1elfd',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert1_elfd=intarr(1)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'d3e_pert1_phon',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rf1phon',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert1_phon=intarr(1)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'d3e_pert2_atpol',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'rf2atpol',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert2_atpol(1:2)=intarr(1:2)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'d3e_pert2_dir',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'rf2dir',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert2_dir(1:3)=intarr(1:3)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'d3e_pert2_elfd',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rf2elfd',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert2_elfd=intarr(1)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'d3e_pert2_phon',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rf2phon',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert2_phon=intarr(1)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'d3e_pert2_strs',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rf2strs',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert2_strs=intarr(1)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'d3e_pert3_atpol',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'rf3atpol',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert3_atpol(1:2)=intarr(1:2)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'d3e_pert3_dir',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'rf3dir',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert3_dir(1:3)=intarr(1:3)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'d3e_pert3_elfd',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rf3elfd',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert3_elfd=intarr(1)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'d3e_pert3_phon',tread,'INT')
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rf3phon',tread_alt,'INT')
   if(tread==1.or.tread_alt==1) dtset%d3e_pert3_phon=intarr(1)
   if (tread_alt==1) tread_key=1

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'lw_flexo',tread,'INT')
   if(tread==1) dtset%lw_flexo=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'lw_qdrpl',tread,'INT')
   if(tread==1) dtset%lw_qdrpl=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prepalw',tread,'INT')
   if(tread==1) dtset%prepalw=intarr(1)
! end if

 ! Recursion input variables
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'tfkinfunc',tread,'INT')
 if(tread==1) dtset%tfkinfunc=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'recgratio',tread,'INT')
 if(tread==1) dtset%recgratio=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'recefermi',tread,'ENE')
 if(tread==1) dtset%recefermi=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'recnpath',tread,'INT')
 if(tread==1) dtset%recnpath=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'recnrec',tread,'INT')
 if(tread==1) dtset%recnrec=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'recrcut',tread,'LEN')
 if(tread==1) dtset%recrcut=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'recptrott',tread,'INT')
 if(tread==1) dtset%recptrott=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rectesteg',tread,'INT')
 if(tread==1) dtset%rectesteg=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rectolden',tread,'DPR')
 if(tread==1) dtset%rectolden=dprarr(1)

 ! Constant NPT Molecular Dynamics Input variables
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'noseinert',tread,'DPR')
 if(tread==1) dtset%noseinert=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bmass',tread,'DPR')
 if(tread==1) dtset%bmass=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,dtset%nnos,string(1:lenstr),'qmass',tread,'DPR')
 if(tread==1) dtset%qmass(:)=dprarr(1:dtset%nnos)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'tphysel',tread,'ENE')
 if(tread==1) dtset%tphysel=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,6,string(1:lenstr),'strtarget',tread,'DPR')
 if(tread==1) dtset%strtarget(1:6)=dprarr(1:6)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'vcutgeo',tread,'DPR')
 if(tread==1) dtset%vcutgeo(1:3)=dprarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'friction',tread,'DPR')
 if(tread==1) dtset%friction=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'mdwall',tread,'LEN')
 if(tread==1) dtset%mdwall=dprarr(1)

 ! Path-Integral Molecular Dynamics Input variables
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'adpimd',tread,'INT')
 if(tread==1) dtset%adpimd=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'adpimd_gamma',tread,'DPR')
 if(tread==1) dtset%adpimd_gamma=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pimd_constraint',tread,'INT')
 if(tread==1) dtset%pimd_constraint=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pitransform',tread,'INT')
 if(tread==1) dtset%pitransform=intarr(1)

 ! Default for pimass is amu
 dtset%pimass(1:ntypat)=dtset%amu_orig(1:ntypat,1)
 ! NOTE: initialisation with the first image only. TO BE MODIFIED ....
 call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),'pimass',tread,'DPR')
 if(tread==1) dtset%pimass(1:ntypat)=dprarr(1:ntypat)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spinmagntarget',tread,'DPR')
 if(tread==1) dtset%spinmagntarget=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'eshift',tread,'ENE')
 if(tread==1) dtset%eshift=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'boxcenter',tread,'DPR')
 if(tread==1) dtset%boxcenter(1:3)=dprarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ecuteps',tread,'ENE')
 if(tread==1) dtset%ecuteps=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ecutsigx',tread,'ENE')
 if(tread==1) dtset%ecutsigx=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ecutwfn',tread,'ENE')
 if(tread==1) then
   dtset%ecutwfn=dprarr(1)
 else
   if(dtset%optdriver==RUNL_SCREENING .or. dtset%optdriver==RUNL_SIGMA) dtset%ecutwfn=dtset%ecut
 end if
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'omegasimax',tread,'ENE')
 if(tread==1) dtset%omegasimax=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'omegasrdmax',tread,'ENE')
 if(tread==1) dtset%omegasrdmax=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'mbpt_sciss',tread,'ENE')
 if(tread==1) dtset%mbpt_sciss=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spbroad',tread,'ENE')
 if(tread==1) dtset%spbroad=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'stmbias',tread,'ENE')
 if(tread==1) dtset%stmbias=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'awtr',tread,'INT')
 if(tread==1) dtset%awtr=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'inclvkb',tread,'INT')
 if(tread==1) dtset%inclvkb=intarr(1)
 if (dtset%inclvkb == 1) then
   MSG_ERROR("inclvkb == 1 is not allowed anymore. Choose between 0 and 1.")
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nomegasf',tread,'INT')
 if(tread==1) dtset%nomegasf=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spmeth',tread,'INT')
 if(tread==1) dtset%spmeth=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'symchi',tread,'INT')
 if(tread==1) dtset%symchi=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getefmas',tread,'INT')
 if(tread==1) dtset%getefmas=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getscr',tread,'INT')
 if(tread==1) dtset%getscr=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getscr_filepath',tread,'KEY', key_value=key_value)
 if(tread==1) dtset%getscr_filepath = key_value

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwgamma',tread,'INT')
 if(tread==1) dtset%gwgamma=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'irdefmas',tread,'INT')
 if(tread==1) dtset%irdefmas=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'irdsuscep',tread,'INT')
 if(tread==1) dtset%irdsuscep=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'irdscr',tread,'INT')
 if(tread==1) dtset%irdscr=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nomegasi',tread,'INT')
 if(tread==1) dtset%nomegasi=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ppmodel',tread,'INT')
 if(tread==1) dtset%ppmodel=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'symsigma',tread,'INT')
 if(tread==1) dtset%symsigma=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fftgw',tread,'INT')
 if(tread==1) dtset%fftgw=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getsuscep',tread,'INT')
 if(tread==1) dtset%getsuscep=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getqps',tread,'INT')
 if(tread==1) dtset%getqps=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwpara',tread,'INT')
 if(tread==1) dtset%gwpara=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'irdqps',tread,'INT')
 if(tread==1) dtset%irdqps=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ppmfrq',tread,'ENE')
 if(tread==1) dtset%ppmfrq=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rcut',tread,'LEN')
 if(tread==1) dtset%rcut=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'zcut',tread,'ENE')
 if(tread==1) dtset%zcut = dprarr(1)
 !else
 !  ! Change default value in EPH calculations
 !  !if (dtset%optdriver == RUNL_EPH) then
 !  !  dtset%zcut = 0.001_dp * eV_Ha
 !  !end if
 !end if

 ! q-points for long wave-length limit.
 if (dtset%gw_nqlwl>0) then
   call intagm(dprarr,intarr,jdtset,marr,3*dtset%gw_nqlwl,string(1:lenstr),'gw_qlwl',tread,'DPR')
   if(tread==1) dtset%gw_qlwl(1:3,1:dtset%gw_nqlwl) = reshape(dprarr(1:3*dtset%gw_nqlwl),(/3,dtset%gw_nqlwl/))
 end if

 !@bethe_salpeter
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bs_nstates',tread,'INT')
 if(tread==1) dtset%bs_nstates=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bs_algorithm',tread,'INT')
 if(tread==1) dtset%bs_algorithm=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bs_haydock_niter',tread,'INT')
 if(tread==1) dtset%bs_haydock_niter=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bs_hayd_term',tread,'INT')
 if(tread==1) dtset%bs_hayd_term=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bs_exchange_term',tread,'INT')
 if(tread==1) dtset%bs_exchange_term=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bs_calctype',tread,'INT')
 if(tread==1) dtset%bs_calctype=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bs_coulomb_term',tread,'INT')
 if(tread==1) dtset%bs_coulomb_term=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bs_coupling',tread,'INT')
 if(tread==1) dtset%bs_coupling=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'irdbseig',tread,'INT')
 if(tread==1) dtset%irdbseig=intarr(1)

 !call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'irdbsr',tread,'INT')
 !if(tread==1) dtset%irdbsr=intarr(1)

 !call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'irdbsc',tread,'INT')
 !if(tread==1) dtset%irdbsc=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getbseig',tread,'INT')
 if(tread==1) dtset%getbseig=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getbsreso',tread,'INT')
 if(tread==1) dtset%getbsreso=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getbscoup',tread,'INT')
 if(tread==1) dtset%getbscoup=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gethaydock',tread,'INT')
 if(tread==1) dtset%gethaydock=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bs_interp_method',tread,'INT')
 if(tread==1) dtset%bs_interp_method=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bs_interp_mode',tread,'INT')
 if(tread==1) dtset%bs_interp_mode=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bs_interp_prep',tread,'INT')
 if(tread==1) dtset%bs_interp_prep=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bs_interp_m3_width',tread,'DPR')
 if(tread==1) dtset%bs_interp_m3_width=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,dtset%nsppol,string(1:lenstr),'bs_loband',tread,'INT')
 if(tread==1) dtset%bs_loband=intarr(1:dtset%nsppol)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'bs_interp_kmult',tread,'INT')
 if(tread==1) dtset%bs_interp_kmult(1:3) = intarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'bs_eh_cutoff',tread,'ENE')
 if(tread==1) dtset%bs_eh_cutoff(1:2)=dprarr(1:2)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'bs_freq_mesh',tread,'ENE')
 if(tread==1) dtset%bs_freq_mesh(1:3)=dprarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'bs_haydock_tol',tread,'DPR')
 if(tread==1) dtset%bs_haydock_tol=dprarr(1:2)

 ntypalch=dtset%ntypalch
 npspalch=dtset%npspalch

 ! Compute ziontypat. When the pseudo-atom is pure, simple copy
 if(ntyppure>0)then
   do itypat=1,ntyppure
     dtset%ziontypat(itypat)=zionpsp(itypat)
   end do
 end if

 ! When the pseudo-atom is alchemical, must make mixing
 if(ntypalch>0)then
   do itypat=ntyppure+1,ntypat
     dtset%ziontypat(itypat)=zero
     do ipsp=ntyppure+1,npsp
       dtset%ziontypat(itypat)=dtset%ziontypat(itypat) &
        + dtset%mixalch_orig(ipsp-ntyppure,itypat-ntyppure,1)*zionpsp(ipsp)
     end do
   end do
 end if

 charge=dtset%charge

 if (occopt==0 .or. occopt==1 .or. (occopt>=3 .and. occopt<=8) ) then
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nband',tnband,'INT')
   if(tnband==1) then
     nband1=intarr(1)
   else
     ! Default value in the metallic case, or in the insulating case
     fband=0.5_dp
     if(occopt==1)fband=0.125_dp
     if((occopt/=1).and.(dtset%accuracy==5.or.dtset%accuracy==6)) fband =0.75_dp
     if (dtset%usewvl == 1) fband = zero
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fband',tfband,'DPR')
     if(tfband==1)then
       fband=dprarr(1)
       write(msg, '(a,es16.8,a)' )' invars2: read the value of fband=',fband,' from input file.'
     else
       write(msg, '(a,es16.8)' )' invars2: take the default value of fband=',fband
     end if
     dtset%fband=fband
     call wrtout(std_out,msg,'COLL')
     ! First compute the total valence charge
     zval=zero
     sum_spinat=zero
     do iatom=1,natom
       zval=zval+dtset%ziontypat(dtset%typat(iatom))
       sum_spinat=sum_spinat+dtset%spinat(3,dtset%typat(iatom))
     end do
     zelect=zval-charge
     ! Then select the minimum number of bands, and add the required number.
     ! Note that this number might be smaller than the one computed
     ! by a slightly different formula in invars1 (difference in fband).
     nband1=dtset%nspinor * ((ceiling(zelect-tol10)+1)/2 + ceiling( fband*natom - tol10 )) &
&     + (nsppol-1)*(ceiling(half*(sum_spinat -tol10)))
   end if

   ! Set nband to same input number for each k point and spin
   ! where nband1 is the eventual input, computed value, or default
   do ikpt=1,nkpt*nsppol
     dtset%nband(ikpt)=nband1
   end do

 else if (occopt==2) then
   ! Give nband explicitly for each k point and spin
   call intagm(dprarr,intarr,jdtset,nkpt*nsppol,nkpt*nsppol,string(1:lenstr),'nband',tnband,'INT')
   if(tnband==1) dtset%nband(1:nkpt*nsppol)=intarr(1:nkpt*nsppol)

 else
   write(msg, '(a,i0,3a)' )'occopt=',occopt,' not allowed.',ch10,'Action: correct your input file.'
   MSG_ERROR(msg)
 end if

 !----------------------------------------------------------------------------
 ! Read other parameters
 ! All checking should be done in chkinp.f

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'auxc_scal',tread,'DPR')
 if(tread==1) dtset%auxc_scal=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'builtintest',tread,'INT')
 if(tread==1) dtset%builtintest=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'expert_user',tread,'INT')
 if(tread==1) dtset%expert_user=intarr(1)

 if(dtset%expert_user==0)then

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'chkdilatmx',tread,'INT')
   if(tread==1) dtset%chkdilatmx=intarr(1)
  
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'chkprim',tread,'INT')
   if(tread==1) dtset%chkprim=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'chksymbreak',tread,'INT')
   if(tread==1) dtset%chksymbreak=intarr(1)
  
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'chksymtnons',tread,'INT')
   if(tread==1) dtset%chksymtnons=intarr(1)

 else

   dtset%chkdilatmx=0
   dtset%chkprim=0
   dtset%chksymbreak=0
   dtset%chksymtnons=0

 endif

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'expert_user',tread,'INT')
 if(tread==1) dtset%expert_user=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fockoptmix',tread,'INT')
 if(tread==1) dtset%fockoptmix=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getocc',tread,'INT')
 if(tread==1) dtset%getocc=intarr(1)
 getocc=dtset%getocc

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getwfk',tread,'INT')
 if(tread==1) dtset%getwfk=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getwfk_filepath',tread,'KEY', key_value=key_value)
 if(tread==1) dtset%getwfk_filepath = key_value

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getwfkfine',tread,'INT')
 if(tread==1) dtset%getwfkfine=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getwfkfine_filepath',tread,'KEY', key_value=key_value)
 if(tread==1) dtset%getwfkfine_filepath = key_value

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getxcart',tread,'INT')
 if(tread==1) dtset%getxcart=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getxred',tread,'INT')
 if(tread==1) dtset%getxred=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getddb',tread,'INT')
 if(tread==1) dtset%getddb=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getdvdb',tread,'INT')
 if(tread==1) dtset%getdvdb=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getden',tread,'INT')
 if(tread==1) dtset%getden=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getden_filepath',tread,'KEY', key_value=key_value)
 if(tread==1) dtset%getden_filepath = key_value

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getpawden',tread,'INT')
 if(tread==1) dtset%getpawden=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getddb_filepath',tread,'KEY', key_value=key_value)
 if(tread==1) dtset%getddb_filepath = key_value

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getdvdb_filepath',tread,'KEY', key_value=key_value)
 if(tread==1) dtset%getdvdb_filepath = key_value

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getpot_filepath',tread,'KEY', key_value=key_value)
 if(tread==1) dtset%getpot_filepath = key_value

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getsigeph_filepath',tread,'KEY', key_value=key_value)
 if(tread==1) dtset%getsigeph_filepath = key_value

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getcell',tread,'INT')
 if(tread==1) dtset%getcell=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getgam_eig2nkq',tread,'INT')
 if(tread==1) dtset%getgam_eig2nkq=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getwfq',tread,'INT')
 if(tread==1) dtset%getwfq=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getwfq_filepath',tread,'KEY', key_value=key_value)
 if(tread==1) dtset%getwfq_filepath = key_value

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'get1wf',tread,'INT')
 if(tread==1) dtset%get1wf=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getddk',tread,'INT')
 if(tread==1) dtset%getddk=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getdelfd',tread,'INT')
 if(tread==1) dtset%getdelfd=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getdkdk',tread,'INT')
 if(tread==1) dtset%getdkdk=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getdkde',tread,'INT')
 if(tread==1) dtset%getdkde=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getvel',tread,'INT')
 if(tread==1) dtset%getvel=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'get1den',tread,'INT')
 if(tread==1) dtset%get1den=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'iomode',tread,'INT')
 if(tread==1) dtset%iomode=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'enunit',tread,'INT')
 if(tread==1) dtset%enunit=intarr(1)

 ! eph variables
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'asr',tread,'INT')
 if(tread==1) dtset%asr=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dipdip',tread,'INT')
 if(tread==1) dtset%dipdip=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'chneut',tread,'INT')
 if(tread==1) dtset%chneut=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'symdynmat',tread,'INT')
 if(tread==1) dtset%symdynmat=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'symv1scf',tread,'INT')
 if(tread==1) dtset%symv1scf = intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'eph_restart',tread,'INT')
 if(tread==1) dtset%eph_restart=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'eph_task',tread,'INT')
 if(tread==1) dtset%eph_task=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'eph_mustar',tread,'DPR')
 if(tread==1) dtset%eph_mustar=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'eph_tols_idelta',tread,'DPR')
 if (tread == 1) dtset%eph_tols_idelta = dprarr(1:2)

 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'eph_phrange',tread,'INT')
 if (tread == 1) dtset%eph_phrange = intarr(1:2)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'eph_intmeth',tread,'INT')
 if (tread == 1) then
   dtset%eph_intmeth = intarr(1)
 else
   ! eph_intmeth depends on eph_task.
   !if (abs(dtset%eph_task) == 4) dtset%eph_intmeth = 1
   if (dtset%eph_task == +4) dtset%eph_intmeth = 1
   if (dtset%eph_task == -4 .and. dtset%symsigma == 0) dtset%eph_intmeth = 1
 end if

 ! Allow use to dope the system or shift artificially the Fermi level
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'eph_extrael',tread_extrael,'DPR')
 if (tread_extrael == 1) dtset%eph_extrael = dprarr(1)

 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'eph_doping', tread, 'DPR')
 if (tread == 1) then
   ABI_CHECK(tread_extrael == 0, "eph_extrael and eph_doping are mutually exclusive!")
   ! Units of eph_doping is e_charge / cm^3
   dtset%eph_extrael = - dprarr(1) * ucvol * (Bohr_meter * 100) ** 3
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'eph_fermie',tread,'ENE')
 if(tread==1) dtset%eph_fermie=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'eph_frohlichm',tread,'INT')
 if(tread==1) dtset%eph_frohlichm=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'eph_fsmear',tread,'ENE')
 if(tread==1) dtset%eph_fsmear=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'eph_fsewin',tread,'ENE')
 if(tread==1) dtset%eph_fsewin=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'eph_ecutosc',tread,'ENE')
 if(tread==1) dtset%eph_ecutosc=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'eph_phwinfact',tread,'DPR')
 if(tread==1) dtset%eph_phwinfact = dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'eph_ngqpt_fine',tread,'INT')
 if(tread==1) dtset%eph_ngqpt_fine=intarr(1:3)

 narr = size(dtset%eph_np_pqbks)
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),'eph_np_pqbks',tread,'INT')
 if (tread==1) dtset%eph_np_pqbks = intarr(1:narr)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'eph_stern',tread,'INT')
 if(tread==1) dtset%eph_stern = intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'eph_use_ftinterp',tread,'INT')
 if(tread==1) dtset%eph_use_ftinterp = intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'eph_transport',tread,'INT')
 if(tread==1) dtset%eph_transport=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ph_wstep',tread,'ENE')
 if(tread==1) dtset%ph_wstep=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ph_intmeth',tread,'INT')
 if(tread==1) dtset%ph_intmeth=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ph_smear',tread,'ENE')
 if(tread==1) dtset%ph_smear=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'ddb_ngqpt',tread,'INT')
 if(tread==1) dtset%ddb_ngqpt=intarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'ddb_shiftq',tread,'DPR')
 if(tread==1) dtset%ddb_shiftq=dprarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dvdb_qcache_mb',tread,'DPR')
 if(tread==1) dtset%dvdb_qcache_mb = dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dvdb_qdamp',tread,'DPR')
 if(tread==1) dtset%dvdb_qdamp = dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dvdb_rspace_cell',tread, 'INT')
 if(tread==1) dtset%dvdb_rspace_cell = intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dvdb_add_lr',tread,'INT')
 if(tread==1) dtset%dvdb_add_lr = intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ph_ndivsm',tread,'INT')
 if(tread==1) dtset%ph_ndivsm=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ph_nqpath',tread,'INT')
 if(tread==1) dtset%ph_nqpath=intarr(1)

 if (dtset%ph_nqpath > 0) then
   ! Read qpath for phonons
   ABI_MALLOC(dtset%ph_qpath, (3, dtset%ph_nqpath))
   ABI_CHECK(3 * dtset%ph_nqpath <= marr, "3 * dtset%ph_nqpath > marr!")
   call intagm(dprarr,intarr,jdtset,marr,3*dtset%ph_nqpath,string(1:lenstr),'ph_qpath',tread,'DPR')
   if (tread == 0) then
     MSG_ERROR("When ph_nqpath > 0, ph_qpath should be specified")
   end if
   dtset%ph_qpath = reshape(dprarr(1:3*dtset%ph_nqpath), [3, dtset%ph_nqpath])
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ph_nqshift',tread,'INT')
 if(tread==1) dtset%ph_nqshift=intarr(1)

 if (dtset%ph_nqshift > 0) then
   ! Read ph_qshift for phonons, default is [0,0,0]
   ABI_CALLOC(dtset%ph_qshift, (3, dtset%ph_nqshift))
   if (tread == 1) then
     ABI_CHECK(3 * dtset%ph_nqshift <= marr, "3 * dtset%ph_nqshift > marr!")
     call intagm(dprarr,intarr,jdtset,marr,3*dtset%ph_nqshift,string(1:lenstr),'ph_qshift',tread,'DPR')
     if (tread == 0) then
       MSG_ERROR("When ph_nqshift > 0, ph_qshift should be specified")
     end if
     dtset%ph_qshift = reshape(dprarr(1:3*dtset%ph_nqshift), [3, dtset%ph_nqshift])
   end if
 end if

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'ph_ngqpt',tread,'INT')
 if(tread==1) dtset%ph_ngqpt=intarr(1:3)
!end e-ph variables

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'exchn2n3d',tread,'INT')
 if(tread==1) dtset%exchn2n3d=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'extrapwf',tread,'INT')
 if(tread==1) dtset%extrapwf=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'iboxcut',tread,'INT')
 if(tread==1) dtset%iboxcut=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'icutcoul',tread,'INT')
 if(tread==1) dtset%icutcoul=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'imgwfstor',tread,'INT')
 if(tread==1) dtset%imgwfstor=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'mdf_epsinf',tread,'DPR')
 if(tread==1) dtset%mdf_epsinf=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'mem_test',tread,'INT')
 if(tread==1) dtset%mem_test=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'mixprec',tread,'INT')
 if (tread==1) dtset%mixprec = intarr(1)

 if (dtset%imgmov==0.or.dtset%imgmov==2.or.dtset%imgmov==5) then
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'mep_mxstep',tread,'LEN')
   if(tread==1) then
     dtset%mep_mxstep=dprarr(1)
   else if(dtset%imgmov==5) then
     dtset%mep_mxstep=0.4_dp
   end if
 end if
 if (dtset%imgmov==2.or.dtset%imgmov==5) then
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'mep_solver',tread,'INT')
   if(tread==1) then
     dtset%mep_solver=intarr(1)
   else if(dtset%imgmov==2) then
     dtset%mep_solver=0
   else if(dtset%imgmov==5) then
     dtset%mep_solver=0
   end if
   if (dtset%imgmov==2) then
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'string_algo',tread,'INT')
     if(tread==1) dtset%string_algo=intarr(1)
   end if
   if (dtset%imgmov==5) then
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'neb_algo',tread,'INT')
     if(tread==1) dtset%neb_algo=intarr(1)
     call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'neb_spring',tread,'DPR')
     if(tread==1) then
       dtset%neb_spring(1:2)=dprarr(1:2)
     else if (dtset%neb_algo==2) then
       dtset%neb_spring(1:2)=(/0.02_dp,0.15_dp/)
     end if
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'cineb_start',tread,'INT')
     if(tread==1) dtset%cineb_start=intarr(1)
   end if
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'istatimg',tread,'INT')
 if(tread==1) dtset%istatimg=intarr(1)

 ! variables for random positions in unit cell
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'random_atpos',tread,'INT')
 if(tread==1) dtset%random_atpos=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ionmov',tread,'INT')
 if(tread==1) dtset%ionmov=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'intxc',tread,'INT')
 if(tread==1) dtset%intxc=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'densfor_pred',tread,'INT')
 if(tread==1) then
   dtset%densfor_pred=intarr(1)
 else
   if (dtset%paral_kgb==1) dtset%densfor_pred=6
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'hmcsst',tread,'INT')
 if(tread==1) dtset%hmcsst=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'hmctt',tread,'INT')
 if(tread==1) dtset%hmctt=intarr(1)

 ! For the moment LOTF does not use different from 2
 if(dtset%ionmov==23) then
#ifdef HAVE_LOTF
   if(dtset%densfor_pred/=2)then
     densfor_pred=2
     dtset%densfor_pred=densfor_pred
     write(msg, '(4a)' )&
       'When ionmov==23, densfor_pred must be 2.',ch10,&
       'Set densfor_pred to 2.',ch10
     MSG_COMMENT(msg)
   end if
#else
   dtset%ionmov=12
   write(msg, '(4a)' )&
     'LOTF is disabled, ionmov can not be 23.',ch10,&
    'Set ionmov to 12.',ch10
   MSG_COMMENT(msg)
#endif
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'iprcel',tread,'INT')
 if(tread==1) dtset%iprcel=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'iprcfc',tread,'INT')
 if(tread==1) dtset%iprcfc=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'irandom',tread,'INT')
 if(tread==1) dtset%irandom=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'irdddb',tread,'INT')
 if(tread==1) dtset%irdddb=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'irdden',tread,'INT')
 if(tread==1) dtset%irdden=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'irddvdb',tread,'INT')
 if(tread==1) dtset%irddvdb = intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'irdpawden',tread,'INT')
 if(tread==1) dtset%irdpawden=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'irdhaydock',tread,'INT')
 if(tread==1) dtset%irdhaydock=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'irdvdw',tread,'INT')
 if(tread==1) dtset%irdvdw=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'irdwfk',tread,'INT')
 if(tread==1) dtset%irdwfk=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'irdwfkfine',tread,'INT')
 if(tread==1) dtset%irdwfkfine=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'iscf',tread,'INT')
 if(tread==1) then
   dtset%iscf=intarr(1)
   if (dtset%usewvl==1) then
     !wvl_bigdft_comp should be 1 for iscf=0, 0 for iscf>0, except if it is set by the user
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'wvl_bigdft_comp',tread,'INT')
     if (tread==0.and.dtset%iscf==0) dtset%wvl_bigdft_comp=1
     if (tread==0.and.dtset%iscf >0) dtset%wvl_bigdft_comp=0
   end if
 else if (dtset%optdriver==RUNL_RESPFN.and.dtset%iscf>=10) then
   dtset%iscf=dtset%iscf-10
 else if (dtset%optdriver==RUNL_GWLS) then
   dtset%iscf=-2
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'isecur',tread,'INT')
 if(tread==1) dtset%isecur=intarr(1)

 ! Reading ixc must be immediately followed by reading xcname
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ixc',tread,'INT')
 if(tread==1) dtset%ixc=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),"xcname",tread_key,'KEY',key_value=key_value)
 if(tread_key==1)then
   if(tread==1)then
     write(msg, '(5a)' )&
      'ixc and xcname cannot be specified simultaneously',ch10,&
      'for the same dataset.',ch10,'Action: check the input file.'
     MSG_ERROR(msg)
   else
     !Note that xcname is a 'key' variable : its value is stored in keyw at output of intagm
     if(trim(key_value) == 'PW92') then
        dtset%ixc=7
     else
       MSG_ERROR(sjoin("Don't know how to convert xcname", key_value, "to ixc"))
     end if
     tread=1
   end if
 end if
 ixc_current=dtset%ixc

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ixcrot',tread,'INT')
 if(tread==1) dtset%ixcrot=intarr(1)

 ! Read the ixc for an advanced functional
 ! If present, and relevant (only specific values for gcalctyp),
 ! the other internal variable will be adjusted to this other functional)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ixc_sigma',tread,'INT')
 if(tread==1)then
   dtset%ixc_sigma=intarr(1)
   if( dtset%optdriver==RUNL_SIGMA .and. mod(dtset%gwcalctyp,10)==5)ixc_current=dtset%ixc_sigma
 end if

 ! Initialize xclevel and usefock
 call get_xclevel(ixc_current,dtset%xclevel,dtset%usefock)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'auxc_ixc',tread,'INT')
 if(tread==1) dtset%auxc_ixc=intarr(1)
 ! If the default value had been given, possibly switch on
 ! the auxc_ixc corresponding to ixc, if the latter is an hybrid
 if(dtset%auxc_ixc==0)then
   call get_auxc_ixc(dtset%auxc_ixc,ixc_current)
 end if

 ! Now take care of the parameters for hybrid functionals
 if(dtset%usefock==1)then

   if(ixc_current ==40 .or. ixc_current ==41 .or. ixc_current ==42)then
     dtset%hyb_mixing_sr=zero
     dtset%hyb_range_dft=zero ; dtset%hyb_range_fock=zero
     if(ixc_current==40)dtset%hyb_mixing=one
     if(ixc_current==41)dtset%hyb_mixing=quarter
     if(ixc_current==42)dtset%hyb_mixing=third
   else if(ixc_current==-427)then
     ! Special case of HSE03
     dtset%hyb_mixing=zero  ; dtset%hyb_mixing_sr=quarter
     dtset%hyb_range_dft=0.15_dp*two**third  ; dtset%hyb_range_fock=0.15_dp*sqrt(half)
   else if (ixc_current<0) then
     call libxc_functionals_init(ixc_current,dtset%nspden,xc_tb09_c=dtset%xc_tb09_c)
     call libxc_functionals_get_hybridparams(hyb_mixing=dtset%hyb_mixing,hyb_mixing_sr=dtset%hyb_mixing_sr,&
       hyb_range=dtset%hyb_range_dft)
     call libxc_functionals_end()
     dtset%hyb_range_fock=dtset%hyb_range_dft
   end if

   ! Warning: the user-defined parameters for hybrids are by convention stored as negative numbers
   ! This trick will allow to echo them, and only them.
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'hyb_mixing',tread,'DPR')
   if(tread==1)then
     if(dprarr(1)<-tol14)then
       write(msg, '(a,es16.8,2a)' )&
       'A negative value for hyb_mixing is not allowed, while at input hyb_mixing=',dprarr(1),ch10,&
       'Action: modify hyb_mixing in the input file.'
       MSG_ERROR(msg)
     end if
     dtset%hyb_mixing=-dprarr(1) ! Note the minus sign
   end if

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'hyb_mixing_sr',tread,'DPR')
   if(tread==1)then
     if(dprarr(1)<-tol14)then
       write(msg, '(a,es16.8,2a)' )&
       'A negative value for hyb_mixing_sr is not allowed, while at input hyb_mixing_sr=',dprarr(1),ch10,&
       'Action: modify hyb_mixing_sr in the input file.'
       MSG_ERROR(msg)
     end if
     dtset%hyb_mixing_sr=-dprarr(1) ! Note the minus sign
   end if

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'hyb_range_dft',tread_dft,'DPR')
   if(tread_dft==1)then
     if(dprarr(1)<-tol14)then
       write(msg, '(a,es16.8,2a)' )&
       'A negative value for hyb_range_dft is not allowed, while at input hyb_range_dft=',dprarr(1),ch10,&
       'Action: modify hyb_range_dft in the input file.'
       MSG_ERROR(msg)
     end if
     dtset%hyb_range_dft=-dprarr(1) ! Note the minus sign
   end if

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'hyb_range_fock',tread_fock,'DPR')
   if(tread_fock==1)then
     if(dprarr(1)<-tol14)then
       write(msg, '(a,es16.8,2a)' )&
       'A negative value for hyb_range_fock is not allowed, while at input hyb_range_fock=',dprarr(1),ch10,&
       'Action: modify hyb_range_fock in the input file.'
       MSG_ERROR(msg)
     end if
     dtset%hyb_range_fock=-dprarr(1) ! Note the minus sign
   end if

   if(tread_fock==1 .and. tread_dft==0)dtset%hyb_range_dft=dtset%hyb_range_fock
   if(tread_fock==0 .and. tread_dft==1)dtset%hyb_range_fock=dtset%hyb_range_dft

 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_acutmin',tread,'DPR')
 if(tread==1) dtset%vdw_df_acutmin=dprarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_aratio',tread,'DPR')
 if(tread==1) dtset%vdw_df_aratio=dprarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_damax',tread,'DPR')
 if(tread==1) dtset%vdw_df_damax=dprarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_damin',tread,'DPR')
 if(tread==1) dtset%vdw_df_damin=dprarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_dcut',tread,'DPR')
 if(tread==1) dtset%vdw_df_dcut=dprarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_dratio',tread,'DPR')
 if(tread==1) dtset%vdw_df_dratio=dprarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_dsoft',tread,'DPR')
 if(tread==1) dtset%vdw_df_dsoft=dprarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_gcut',tread,'DPR')
 if(tread==1) dtset%vdw_df_gcut=dprarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_ndpts',tread,'INT')
 if(tread==1) dtset%vdw_df_ndpts=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_ngpts',tread,'INT')
 if(tread==1) dtset%vdw_df_ngpts=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_nqpts',tread,'INT')
 if(tread==1) dtset%vdw_df_nqpts=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_nrpts',tread,'INT')
 if(tread==1) dtset%vdw_df_nrpts=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_nsmooth',tread,'INT')
 if(tread==1) dtset%vdw_df_nsmooth=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_phisoft',tread,'DPR')
 if(tread==1) dtset%vdw_df_phisoft=dprarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_qcut',tread,'DPR')
 if(tread==1) dtset%vdw_df_qcut=dprarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_qratio',tread,'DPR')
 if(tread==1) dtset%vdw_df_qratio=dprarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_rcut',tread,'DPR')
 if(tread==1) dtset%vdw_df_rcut=dprarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_rsoft',tread,'DPR')
 if(tread==1) dtset%vdw_df_rsoft=dprarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_tolerance',tread,'DPR')
 if(tread==1) dtset%vdw_df_tolerance=dprarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_threshold',tread,'DPR',ds_input)
 if(tread==1) dtset%vdw_df_threshold=dprarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_tweaks',tread,'INT')
 if(tread==1) dtset%vdw_df_tweaks=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_df_zab',tread,'DPR')
 if(tread==1) dtset%vdw_df_zab=dprarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_xc',tread,'INT')
 if(tread==1) dtset%vdw_xc=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'posdoppler',tread,'INT')
 if(tread==1) dtset%posdoppler=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'positron',tread,'INT')
 if(tread==1) dtset%positron=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ixcpositron',tread,'INT')
 if(tread==1) dtset%ixcpositron=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'posnstep',tread,'INT')
 if(tread==1) dtset%posnstep=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'posocc',tread,'DPR')
 if(tread==1) dtset%posocc=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'frzfermi',tread,'INT')
 if(tread==1) dtset%frzfermi=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nqpt',tread,'INT')
 if(tread==1) dtset%nqpt=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'efmas',tread,'INT')
 if(tread==1) dtset%efmas=intarr(1)

 if(dtset%efmas>0) then
   call intagm(dprarr,intarr,jdtset,marr,2*nkpt,string(1:lenstr),'efmas_bands',tread,'INT')
   if(tread==1) then
     dtset%efmas_bands(1:2,1:nkpt)=reshape(intarr(1:2*nkpt),(/2,nkpt/))
   else
     dtset%efmas_bands(1,:)=1
     dtset%efmas_bands(2,:)=dtset%nband(:)
   end if
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'efmas_deg',tread,'INT')
 if(tread==1) dtset%efmas_deg=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'efmas_deg_tol',tread,'ENE')
 if(tread==1) dtset%efmas_deg_tol=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'efmas_dim',tread,'INT')
 if(tread==1) dtset%efmas_dim=intarr(1)

 if(dtset%efmas_calc_dirs/=0 .and. dtset%efmas_n_dirs>0) then
   ndim=2+MERGE(1,0,ABS(dtset%efmas_calc_dirs)<3)
   call intagm(dprarr,intarr,jdtset,marr,ndim*dtset%efmas_n_dirs,string(1:lenstr),'efmas_dirs',tread,'DPR')
   if(tread==1) then
     dtset%efmas_dirs(1:ndim,1:dtset%efmas_n_dirs)=reshape(dprarr(1:ndim*dtset%efmas_n_dirs),(/ndim,dtset%efmas_n_dirs/))
   end if
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'efmas_ntheta',tread,'INT')
 if(tread==1) dtset%efmas_ntheta=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ieig2rf',tread,'INT')
 if(tread==1) dtset%ieig2rf=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'smdelta',tread,'INT')
 if(tread==1) dtset%smdelta=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bdeigrf',tread,'INT')
 if(tread==1) dtset%bdeigrf=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'restartxf',tread,'INT')
 if(tread==1) dtset%restartxf=intarr(1)
 if (dtset%restartxf == 1) then
   MSG_ERROR("restartxf == 1 has been removed in Abinit8. Use 0,-1,-2-3")
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'optcell',tread,'INT')
 if(tread==1) dtset%optcell=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'irdwfq',tread,'INT')
 if(tread==1) dtset%irdwfq=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ird1den',tread,'INT')
 if(tread==1) dtset%ird1den=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ird1wf',tread,'INT')
 if(tread==1) dtset%ird1wf=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'irdddk',tread,'INT')
 if(tread==1) dtset%irdddk=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'kptopt',tread,'INT')
 if(tread==1) dtset%kptopt=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'chkexit',tread,'INT')
 if(tread==1) dtset%chkexit=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nbdblock',tread,'INT')
 if(tread==1) dtset%nbdblock=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nbdbuf',tread,'INT')
 if(tread==1)then
   dtset%nbdbuf=intarr(1)
   ! A negative value is interpreted as percentage of nband
   if (dtset%nbdbuf < 0) then
     ABI_CHECK(abs(dtset%nbdbuf) < 100, "abs(nbdbuf) should be < 100")
     dtset%nbdbuf = max(nint(abs(dtset%nbdbuf) / 100.0_dp * maxval(dtset%nband)), 1)
   end if
 else
   if(response/=1 .and. dtset%iscf<0)dtset%nbdbuf=2*dtset%nspinor
   if(response==1 .and. 3<=occopt .and. occopt<=8 )dtset%nbdbuf=2*dtset%nspinor
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'localrdwf',tread,'INT')
 if(tread==1) dtset%localrdwf=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'magconon',tread,'INT')
 if(tread==1) dtset%magconon=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'magcon_lambda',tread,'DPR')
 if(tread==1) dtset%magcon_lambda=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ratsm',tread,'DPR')
 if(tread==1) dtset%ratsm=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'optforces',tread,'INT')
 if(tread==1) dtset%optforces=intarr(1)
 if(dtset%usedmft>0.and.dtset%optforces/=0) then
   write(msg, '(3a,i0)' )&
    'When DFT+DMFT is activated ', ch10, &
    'optforces must be equal to 0 instead of ',dtset%optforces
   MSG_ERROR(msg)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'optstress',tread,'INT')
 if(tread==1) dtset%optstress=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'optnlxccc',tread,'INT')
 if(tread==1) dtset%optnlxccc=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nberry',tread,'INT')
 if(tread==1) dtset%nberry=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nc_xccc_gspace',tread,'INT')
 if(tread==1) dtset%nc_xccc_gspace=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,2*nsppol,string(1:lenstr),'bdberry',tread,'INT')
 if(tread==1) then
   dtset%bdberry(1)=intarr(1); dtset%bdberry(2)=intarr(2)
   if(nsppol==2)then
     dtset%bdberry(3)=intarr(3); dtset%bdberry(4)=intarr(4)
   end if
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'berrystep',tread,'INT')
 if(tread==1) dtset%berrystep=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'delayperm',tread,'INT')
 if(tread==1) dtset%delayperm=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'signperm',tread,'INT')
 if(tread==1) dtset%signperm=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nbandkss',tread,'INT')
 if(tread==1) dtset%nbandkss=intarr(1)
 if ( dtset%usedmft > 0  .and. (dtset%nbandkss==0.or.dtset%iscf<0)) then
   if (dtset%usepawu==4.or.dtset%usepawu==14)  dtset%usepawu=14
   if (dtset%usepawu/=4.and.dtset%usepawu/=14) dtset%usepawu=10
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'npwkss',tread,'INT')
 if(tread==1) dtset%npwkss=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'kssform',tread,'INT')
 if(tread==1) dtset%kssform=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'td_mexcit',tread,'INT')
 if(tread==1) dtset%td_mexcit=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'npvel',tread,'INT')
 if(tread==1) dtset%npvel=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'npulayit',tread,'INT')
 if(tread==1) dtset%npulayit=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'diismemory',tread,'INT')
 if(tread==1) dtset%diismemory=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'goprecon',tread,'INT')
 if(tread==1) dtset%goprecon=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'goprecprm',tread,'DPR')
 if(tread==1) dtset%goprecprm(1:3)=dprarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nwfshist',tread,'INT')
 if(tread==1) dtset%nwfshist=intarr(1)

 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'nscforder', tread, 'INT')
 if (tread == 1) dtset%nscforder = intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nomegasrd',tread,'INT')
 if(tread==1) dtset%nomegasrd=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),'corecs',tread,'DPR')
 if(tread==1)then
   dtset%corecs(1:ntypat)=dprarr(1:ntypat)
 end if

 call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),'ptcharge',tread,'DPR')
 if(tread==1)then
   dtset%ptcharge(1:ntypat)=dprarr(1:ntypat)
 end if

 call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),'quadmom',tread,'DPR')
 if(tread==1)then
   dtset%quadmom(1:ntypat)=dprarr(1:ntypat)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawcpxocc',tread,'INT')
 if(tread==1) then
   dtset%pawcpxocc=intarr(1)
 else if (dtset%nspinor==2.and.(dtset%usepawu/=0.or.dtset%usedmft>0)) then
   dtset%pawcpxocc=2
 else if (dtset%pawspnorb>0.and.(dtset%kptopt<=0.or.dtset%kptopt>=3)) then
   if (dtset%optdriver/=RUNL_GSTATE.or.dtset%ionmov<6.or.dtset%iscf<10) dtset%pawcpxocc=2
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawcross',tread,'INT')
 if(tread==1) dtset%pawcross=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawecutdg',tread,'ENE')
 if(tread==1) then
   dtset%pawecutdg=dprarr(1)
 else if (dtset%accuracy==1.or.dtset%accuracy==2) then
   dtset%pawecutdg=dtset%ecut
 else if (dtset%accuracy==3.) then
   dtset%pawecutdg=1.2_dp*dtset%ecut
 else if (dtset%accuracy==4) then
   dtset%pawecutdg=1.5_dp*dtset%ecut
 else if (dtset%accuracy==5.or.dtset%accuracy==6) then
   dtset%pawecutdg=two*dtset%ecut
! NOT YET ACTIVATED
! else if (dtset%accuracy==0.and.usepaw==1.and.dtset%usewvl==0) then
!  dtset%pawecutdg=1.5_dp*dtset%ecut
 end if
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawlcutd',tread,'INT')
 if(tread==1) dtset%pawlcutd=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawlmix',tread,'INT')
 if(tread==1) dtset%pawlmix=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawmixdg',tread,'INT')
 if(tread==1) then
   dtset%pawmixdg=intarr(1)
 else if (dtset%npfft>1.and.usepaw==1) then
   dtset%pawmixdg=1
 else if (dtset%usewvl==1.and.usepaw==1) then
   dtset%pawmixdg=1
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawnhatxc',tread,'INT')
 if(tread==1) dtset%pawnhatxc=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawntheta',tread,'INT')
 if(tread==1) dtset%pawntheta=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawnphi',tread,'INT')
 if(tread==1) dtset%pawnphi=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawnzlm',tread,'INT')
 if(tread==1) dtset%pawnzlm=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawoptmix',tread,'INT')
 if(tread==1) then
   dtset%pawoptmix=intarr(1)
 else
   if (usepaw==1.and.dtset%iscf<10.and.(dtset%positron==1.or.dtset%positron<0)) dtset%pawoptmix=1
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawoptosc',tread,'INT')
 if(tread==1) dtset%pawoptosc=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawovlp',tread,'DPR')
 if(tread==1) dtset%pawovlp=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawprtdos',tread,'INT')
 if(tread==1) dtset%pawprtdos=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawprtvol',tread,'INT')
 if(tread==1) dtset%pawprtvol=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawprtwf',tread,'INT')
 if(tread==1) dtset%pawprtwf=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawprt_k',tread,'INT')
 if(tread==1) dtset%pawprt_k=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawprt_b',tread,'INT')
 if(tread==1) dtset%pawprt_b=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawstgylm',tread,'INT')
 if(tread==1) dtset%pawstgylm=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawsushat',tread,'INT')
 if(tread==1) dtset%pawsushat=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawusecp',tread,'INT')
 if(tread==1) dtset%pawusecp=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawxcdev',tread,'INT')
 if(tread==1) dtset%pawxcdev=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'spnorbscl',tread,'DPR')
 if(tread==1) dtset%spnorbscl=dprarr(1)

 if (dtset%usedmft>0) then
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmft_dc',tread,'INT')
   if(tread==1) dtset%dmft_dc=intarr(1)
   if (dtset%usepawu==14.and.dtset%dmft_dc/=5) then
     write(msg, '(a,a,a)' )&
      'usepawu == 4 and usedmft == 1, dmft_dc should be equal to 5 ',ch10,&
      'imposing dmft_dc = 5'
     MSG_WARNING(msg)
     dtset%dmft_dc=5
   end if
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmft_iter',tread,'INT')
   if(tread==1) dtset%dmft_iter=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmft_kspectralfunc',tread,'INT')
   if(tread==1) dtset%dmft_kspectralfunc=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmft_mxsf',tread,'DPR')
   if(tread==1) dtset%dmft_mxsf=dprarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmft_nwli',tread,'INT')
   if(tread==1) dtset%dmft_nwli=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmft_nwlo',tread,'INT')
   if(tread==1) dtset%dmft_nwlo=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmft_occnd_imag',tread,'INT')
   if(tread==1) dtset%dmft_occnd_imag=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmft_read_occnd',tread,'INT')
   if(tread==1) dtset%dmft_read_occnd=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmft_rslf',tread,'INT')
   if(tread==1) dtset%dmft_rslf=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmft_solv',tread,'INT')
   if(tread==1) dtset%dmft_solv=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmft_t2g',tread,'INT')
   if(tread==1) dtset%dmft_t2g=intarr(1)
!  call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmft_x2my2d',tread,'INT')
!  if(tread==1) dtset%dmft_x2my2d=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmft_tolfreq',tread,'DPR')
   if(tread==1) dtset%dmft_tolfreq=dprarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmft_tollc',tread,'DPR')
   if(tread==1) dtset%dmft_tollc=dprarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmft_charge_prec',tread,'DPR')
   if(tread==1) dtset%dmft_charge_prec=dprarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmftbandi',tread,'INT')
   if(tread==1) dtset%dmftbandi=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmftbandf',tread,'INT')
   if(tread==1) dtset%dmftbandf=intarr(1)
   if((dtset%dmftbandf-dtset%dmftbandi+1)<2*maxval(dtset%lpawu(:))+1.and.&
!     ((dtset%dmft_t2g==0).and.(dtset%dmft_x2my2d==0))) then
      (dtset%dmft_t2g==0) ) then
     write(msg, '(4a,i2,2a)' )&
     '   dmftbandf-dmftbandi+1)<2*max(lpawu(:))+1)',ch10, &
     '   Number of bands to construct Wannier functions is not', &
     ' sufficient to build Wannier functions for l=',maxval(dtset%lpawu(:)),ch10, &
     '   Action: select a correct number of KS bands with dmftbandi and dmftbandf.'
     MSG_ERROR(msg)
   end if

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmftcheck',tread,'INT')
   if(tread==1) dtset%dmftcheck=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmft_entropy',tread,'INT')
   if(tread==1) dtset%dmft_entropy=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmft_nlambda',tread,'INT')
   if(tread==1) dtset%dmft_nlambda=intarr(1)
   if(dtset%dmft_solv>=4) then
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmftqmc_n',tread,'DPR')
     if(tread==1) then
       dtset%dmftqmc_n=dprarr(1)
     else if(dtset%ucrpa==0) then
     else
       write(msg, '(5a)' )&
        'When DFT+DMFT is activated and one of QMC solvers is used,', ch10, &
        'dmftqmc_n MUST be defined.',ch10,&
        'Action: add dmftqmc_n keyword in input file.'
       MSG_ERROR(msg)
     end if
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmftqmc_l',tread,'INT')
     if(tread==1) then
       dtset%dmftqmc_l=intarr(1)
     else if(dtset%ucrpa==0.and.dtset%dmft_solv/=9) then
       write(msg, '(5a)' )&
        'When DFT+DMFT is activated and one of QMC solvers is used,', ch10, &
        'dmftqmc_l MUST be defined.',ch10,&
        'Action: add dmftqmc_l keyword in input file.'
       MSG_ERROR(msg)
     end if
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmftqmc_seed',tread,'INT')
     if(tread==1) dtset%dmftqmc_seed=intarr(1)
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmftqmc_therm',tread,'INT')
     if(tread==1) then
       dtset%dmftqmc_therm=intarr(1)
     else if(dtset%ucrpa==0.and.dtset%dmft_solv/=9) then
       write(msg, '(5a)' )&
        'When DFT+DMFT is activated and one of QMC solvers is used,', ch10, &
        'dmftqmc_therm MUST be defined.',ch10,&
        'Action: add dmftqmc_therm keyword in input file.'
       MSG_ERROR(msg)
     end if
     if(dtset%dmft_solv==5.or.dtset%dmft_solv==8.or.dtset%dmft_solv==9) then
    ! if(dtset%dmft_solv==5) then
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmftctqmc_basis',tread,'INT')
       if(tread==1) dtset%dmftctqmc_basis  =intarr(1)
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmftctqmc_check',tread,'INT')
       if(tread==1) dtset%dmftctqmc_check  =intarr(1)
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmftctqmc_correl',tread,'INT')
       if(tread==1) dtset%dmftctqmc_correl =intarr(1)
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmftctqmc_gmove',tread,'INT')
       if(tread==1) dtset%dmftctqmc_gmove  =intarr(1)
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmftctqmc_grnns',tread,'INT')
       if(tread==1) dtset%dmftctqmc_grnns  =intarr(1)
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmftctqmc_meas',tread,'INT')
       if(tread==1) dtset%dmftctqmc_meas   =intarr(1)
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmftctqmc_mrka',tread,'INT')
       if(tread==1) dtset%dmftctqmc_mrka   =intarr(1)
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmftctqmc_mov',tread,'INT')
       if(tread==1) dtset%dmftctqmc_mov    =intarr(1)
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmftctqmc_order',tread,'INT')
       if(tread==1) dtset%dmftctqmc_order  =intarr(1)
     end if
     if(dtset%dmft_solv>=6.and.dtset%dmft_solv<=7) then
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmftctqmc_triqs_nleg',tread,'INT')
       if(tread==1) dtset%dmftctqmc_triqs_nleg  =intarr(1)
     end if
   end if
 end if

 if (dtset%usepawu/=0.or.dtset%usedmft>0) then
   call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),'f4of2_sla',tread,'ENE')
   if(tread==1) dtset%f4of2_sla(1:ntypat)=dprarr(1:ntypat)
   call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),'f6of2_sla',tread,'ENE')
   if(tread==1) dtset%f6of2_sla(1:ntypat)=dprarr(1:ntypat)
 end if
 if (dtset%usepawu>0.or.dtset%usedmft>0) then
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmatpuopt',tread,'INT')
   if(tread==1) dtset%dmatpuopt=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dmatudiag',tread,'INT')
   if(tread==1) dtset%dmatudiag=intarr(1)
 end if

 if (dtset%macro_uj>0) then
   dtset%pawujat=minval(pack((/ (ii,ii=1,dtset%natom) /) ,dtset%lpawu(dtset%typat(:))>0))
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawujat',tread,'INT')
   if(tread==1) dtset%pawujat=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawujrad',tread,'LEN')
   if(tread==1) dtset%pawujrad=dprarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawujv',tread,'ENE')
   if(tread==1) dtset%pawujv=dprarr(1)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawfatbnd',tread,'INT')
 if(tread==1) dtset%pawfatbnd=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bxctmindg',tread,'DPR')
 if(tread==1) then
   dtset%bxctmindg=dprarr(1)
 else if(usepaw==1.and.dtset%accuracy==0)then
   dtset%bxctmindg=dtset%boxcutmin
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'usexcnhat',tread,'INT')
 if(tread==1) dtset%usexcnhat_orig=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'useylm',tread,'INT')
 if(tread==1) then
   dtset%useylm=intarr(1)
   if ((usepaw==1).and.(dtset%useylm==0)) then
     write(msg, '(5a)' )&
      'Pseudopotential file is PAW format (pspcod=7 or 17) while',ch10,&
      'input variable "useylm" has the incompatible value 0 !',ch10,&
      'Action: change psp format or "useylm" value in your input file.'
     MSG_ERROR(msg)
   end if
   if ((dtset%tfkinfunc==2).and.(dtset%useylm==0)) then
     write(msg, '(a,a,a,a,a)' )&
      'You are using recursion method (tfkinfunc=2)  while',ch10,&
      'input variable "useylm" has the incompatible value 0 !',ch10,&
      'Action: change  "useylm" value in your input file.'
     MSG_ERROR(msg)
   end if
   if ((dtset%efmas==1).and.(dtset%useylm==0)) then
     write(msg, '(a,a,a,a,a)' )&
      'The calculation of effective masses requires the input variable, ',ch10,&
      '"useylm" to be 1, while in your input file, useylm=0',ch10,&
      'Action: change "useylm" value in your input file.'
     MSG_ERROR(msg)
   end if
   if ((dtset%rf2_dkdk/=0).and.(dtset%useylm==0)) then
     write(msg, '(a,a,a,a,a)' )&
      'The calculation of 2nd order k perturbation requires the input variable, ',ch10,&
      '"useylm" to be 1, while in your input file, useylm=0',ch10,&
      'Action: change "useylm" value in your input file.'
     MSG_ERROR(msg)
   end if
   if ((dtset%rf2_dkde/=0).and.(dtset%useylm==0)) then
     write(msg, '(a,a,a,a,a)' )&
      'The calculation of the 2nd order k/Efield perturbation requires the input variable, ',ch10,&
      '"useylm" to be 1, while in your input file, useylm=0',ch10,&
      'Action: change "useylm" value in your input file.'
     MSG_ERROR(msg)
   end if
 end if
 if (usepaw==1) dtset%useylm=1
 if (usepaw==1 .and. dtset%usewvl==1) then
   dtset%useylm=0
 end if
 if (dtset%efmas==1.or.dtset%use_gpu_cuda==1.or.dtset%rf2_dkdk/=0.or.dtset%rf2_dkde/=0) dtset%useylm=1

 if(dtset%tfkinfunc==2 .and. dtset%usewvl==0 ) then
   dtset%useylm=1
   dtset%userec=1
 end if

 ionmov=dtset%ionmov ; densfor_pred=dtset%densfor_pred ; iscf=dtset%iscf ; nqpt=dtset%nqpt
 kptopt=dtset%kptopt; nberry=dtset%nberry ; berryopt=dtset%berryopt

 ! Dielectric real(dp) input variables
 ! Reading of diemix/diemixmag must be inserted after iprcel
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'diecut',tread,'ENE')
 if(tread==1) dtset%diecut=dprarr(1)
 ! Special treatment if iscf==-1
 if(iscf==-1) dtset%diecut=four*dtset%ecut
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dielng',tread,'LEN')
 if(tread==1) dtset%dielng=dprarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'diemac',tread,'DPR')
 if(tread==1) dtset%diemac=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'diemix',tread,'DPR')
 if(tread==1) then
   dtset%diemix=dprarr(1)
 else
   if (mod(dtset%iprcel,100)>19) dtset%diemix=one
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'diemixmag',tread,'DPR')
 if(tread==1) then
   dtset%diemixmag=dprarr(1)
 else
   if (dtset%iscf<10.or.dtset%nspden==1.or.dtset%iprcel==0.or.(dtset%iprcel>70.and.dtset%iprcel<80)) then
     dtset%diemixmag=dtset%diemix
   else
     dtset%diemixmag=-dtset%diemix
   end if
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'diegap',tread,'ENE')
 if(tread==1) dtset%diegap=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dielam',tread,'DPR')
 if(tread==1) dtset%dielam=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'maxestep',tread,'DPR')
 if(tread==1) dtset%maxestep=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'td_maxene',tread,'DPR')
 if(tread==1) dtset%td_maxene=dprarr(1)

 if((iscf==5.or.iscf==6) .and. ionmov==4 .and. densfor_pred/=3 )then
   densfor_pred=3
   dtset%densfor_pred=densfor_pred
   write(msg, '(3a)' )&
   'When ionmov is [4, 5, 6], densfor_pred must be 3.',ch10,&
   'Set densfor_pred to 3.'
   MSG_COMMENT(msg)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'mffmem',tread,'INT')
 if(tread==1) dtset%mffmem=intarr(1)

 ! Set default values of 0 for occupation numbers and k pt wts
 bantot = 0
 ! nkpt and nband must be defined to execute following loop
 if ( tnband == 1 ) then
   do ikpt=1,nkpt*nsppol
     do ii=1,dtset%nband(ikpt)
       bantot=bantot+1
       dtset%occ_orig(bantot,:)=zero
     end do
   end do
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nloc_alg',tread,'INT')
 if(tread==1) dtset%nloalg(1)=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nloc_mem',tread,'INT')
 if(tread==1) then
   dtset%nloalg(2)=1 ; if(intarr(1)<0)dtset%nloalg(2)=-1
   dtset%nloalg(3)=abs(intarr(1))-1
 end if

 ! LOOP variables
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nline',tread,'INT')
 if(tread==1) dtset%nline=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nnsclo',tread,'INT')
 if(tread==1) dtset%nnsclo=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nnsclohf',tread,'INT')
 if(tread==1) dtset%nnsclohf=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nstep',tread,'INT')
 if(tread==1) dtset%nstep=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ntime',tread,'INT')
 if(tread==1) dtset%ntime=intarr(1)

 if (tread == 0) then
   ! if ntime is not given but ionmov > 0 and imgmv is != 0, use a reasonable number of iterations!
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ionmov',tread,'INT')
   if (tread == 1 .and. intarr(1) /= 0 .and. dtset%imgmov == 0) then
     dtset%ntime = 1000
     MSG_COMMENT("Found ionmov /= 0 without ntime in the input. ntime has been set automatically to 1000")
   end if
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nctime',tread,'INT')
 if(tread==1) dtset%nctime=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'orbmag',tread,'INT')
 if(tread==1) dtset%orbmag=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ortalg',tread,'INT')
 if(tread==1) then
   dtset%ortalg=intarr(1)
 else if (dtset%wfoptalg>=10 .and. dtset%ortalg>0) then
   dtset%ortalg=-dtset%ortalg
 end if

 ! Print variables
 call intagm(dprarr,intarr,jdtset,marr,natom,string(1:lenstr),'prtatlist',tread,'INT')
 if(tread==1) then
   dtset%prtatlist=intarr(1:natom)
   call sort_int(natom,dtset%prtatlist,intarr)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtbbb',tread,'INT')
 if(tread==1) dtset%prtbbb=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtbltztrp',tread,'INT')
 if(tread==1) dtset%prtbltztrp=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtcif',tread,'INT')
 if(tread==1) dtset%prtcif=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtden',tread,'INT')
 if(tread==1) dtset%prtden=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtdipole',tread,'INT')
 if(tread==1) dtset%prtdipole=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtdos',tread,'INT')
 if(tread==1) dtset%prtdos=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtdosm',tread,'INT')
 if(tread==1) dtset%prtdosm=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtebands',tread,'INT')
 if(tread==1) dtset%prtebands=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtefg',tread,'INT')
 if(tread==1) dtset%prtefg=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtefmas',tread,'INT')
 if(tread==1) dtset%prtefmas=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prteig',tread,'INT')
 if(tread==1) dtset%prteig=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtelf',tread,'INT')
 if(tread==1) dtset%prtelf=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prteliash',tread,'INT')
 if(tread==1) dtset%prteliash=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtfc',tread,'INT')
 if(tread==1) dtset%prtfc=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtfull1wf',tread,'INT')
 if(tread==1) dtset%prtfull1wf=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtfsurf',tread,'INT')
 if(tread==1) dtset%prtfsurf=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtgsr',tread,'INT')
 if(tread==1) dtset%prtgsr=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtgden',tread,'INT')
 if(tread==1) dtset%prtgden=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtgeo',tread,'INT')
 if(tread==1) dtset%prtgeo=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtgkk',tread,'INT')
 if(tread==1) dtset%prtgkk=intarr(1)

 ! Read usekden and set prtkden to 1 by default.
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'usekden',tread,'INT')
 if(tread==1) then
   dtset%usekden=intarr(1)
 else
   dtset%usekden=merge(1,0,libxc_functionals_ismgga().or.dtset%ixc==31.or.dtset%ixc==34.or.dtset%ixc==35)
 end if
 if (dtset%usekden == 1 .and. dtset%nimage == 1) dtset%prtkden = 1

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtkden',tread,'INT')
 if(tread==1) dtset%prtkden=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtlden',tread,'INT')
 if(tread==1) dtset%prtlden=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtnabla',tread,'INT')
 if(tread==1) dtset%prtnabla=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtnest',tread,'INT')
 if(tread==1) dtset%prtnest=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtphbands',tread,'INT')
 if(tread==1) dtset%prtphbands=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtphdos',tread,'INT')
 if(tread==1) dtset%prtphdos=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtphsurf',tread,'INT')
 if(tread==1) dtset%prtphsurf=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtposcar',tread,'INT')
 if(tread==1) dtset%prtposcar=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtprocar',tread,'INT')
 if(tread==1) dtset%prtprocar=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtpot',tread,'INT')
 if(tread==1) dtset%prtpot=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtpsps',tread,'INT')
 if(tread==1) dtset%prtpsps=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtspcur',tread,'INT')
 if(tread==1) dtset%prtspcur=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtsuscep',tread,'INT')
 if(tread==1) dtset%prtsuscep=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtstm',tread,'INT')
 if(tread==1) dtset%prtstm=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prt1dm',tread,'INT')
 if(tread==1) dtset%prt1dm=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtvclmb',tread,'INT')
 if(tread==1) dtset%prtvclmb=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtvdw',tread,'INT')
 if(tread==1) dtset%prtvdw=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtvha',tread,'INT')
 if(tread==1) dtset%prtvha=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtvhxc',tread,'INT')
 if(tread==1) dtset%prtvhxc=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtkbff',tread,'INT')
 if(tread==1) dtset%prtkbff=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtvol',tread,'INT')
 if(tread==1) dtset%prtvol=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtvolimg',tread,'INT')
 if(tread==1) dtset%prtvolimg=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtvpsp',tread,'INT')
 if(tread==1) dtset%prtvpsp=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtvxc',tread,'INT')
 if(tread==1) dtset%prtvxc=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtwant',tread,'INT')
 if(tread==1) dtset%prtwant=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtwf',tread,'INT')
 if(tread==1) then
   dtset%prtwf=intarr(1)
   if ((dtset%tfkinfunc==2).and.(dtset%prtwf==1)) then
     write(msg, '(5a)' )&
      'You are using recursion method (tfkinfunc=2)  while',ch10,&
      'input variable "prtwf" has the incompatible value 1 !',ch10,&
      'Action: change  "prtwf=0" value in your input file.'
     MSG_ERROR(msg)
   end if
 end if
 if (dtset%tfkinfunc==2) dtset%prtwf=0

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtwf_full',tread,'INT')
 if (tread==1) dtset%prtwf_full=intarr(1)

 if (dtset%prtwf_full == 1 .and. dtset%prtwf == 0) then
   write(msg, '(3a)' )&
   'You are using prtwf_full=1 while prtwf=0',ch10,&
   'Action: set "prtwf=0" in your input file.'
   MSG_ERROR(msg)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtxml',tread,'INT')
 if(tread==1) dtset%prtxml=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ratsph_extra',tread,'DPR')
 if(tread==1) dtset%ratsph_extra=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'qprtrb',tread,'INT')
 if(tread==1) dtset%qprtrb(:)=intarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'strprecon',tread,'DPR')
 if(tread==1) dtset%strprecon=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'tim1rev',tread,'INT')
 if (tread==1) dtset%tim1rev=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'xc_denpos',tread,'DPR')
 if(tread==1) dtset%xc_denpos=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'xc_tb09_c',tread,'DPR')
 if(tread==1) dtset%xc_tb09_c=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,3*dtset%natsph_extra,string(1:lenstr),'xredsph_extra',tread,'DPR')
 if(tread==1) dtset%xredsph_extra=reshape(dprarr(1:3*dtset%natsph_extra), (/3,dtset%natsph_extra/))

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'wfmix',tread,'DPR')
 if(tread==1) dtset%wfmix=dprarr(1)

 ! WVL - Wavelets related values
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'wvl_hgrid',tread,'DPR')
 if(tread==1) dtset%wvl_hgrid=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'wvl_crmult',tread,'DPR')
 if(tread==1) dtset%wvl_crmult=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'wvl_frmult',tread,'DPR')
 if(tread==1) dtset%wvl_frmult=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'wvl_ngauss',tread,'INT')
 if(tread==1) dtset%wvl_ngauss(1:2)=intarr(1:2)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'wvl_nprccg',tread,'INT')
 if(tread==1) dtset%wvl_nprccg=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'tl_nprccg',tread,'INT')
 if(tread==1) dtset%tl_nprccg=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'tl_radius',tread,'DPR')
 if(tread==1) dtset%tl_radius=dprarr(1)

 ! Wannier90 interface related variables
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'w90iniprj',tread,'INT')
 if(tread==1) then
   dtset%w90iniprj=intarr(1)
   if ( usepaw == 0 .and. ( dtset%w90iniprj /= 2  .and. dtset%w90iniprj/=0 .and. dtset%w90iniprj /= 1 )) then
     write(msg, '(a,i0,2a)' )&
      'w90iniprj should be set to 0, 1 or 2, however, it was ',dtset%w90iniprj,ch10,&
      'Action: check the values of w90iniprj.'
     MSG_ERROR(msg)
   end if
   if ( usepaw == 1 .and. ( dtset%w90iniprj < 2 .or. dtset%w90iniprj>6 ) &
      .and. ( dtset%w90iniprj /= 1 .and. dtset%w90iniprj/=0 )) then
     write(msg, '(a,i0,2a)' )&
      'In the PAW case, the only valid values for w90iniprj are 0, 1, 2, 5 and 6 however, it was ',dtset%w90iniprj,ch10,&
      'Action: check the values of w90iniprj.'
     MSG_ERROR(msg)
   end if
 end if

 ! w90prtunk
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'w90prtunk',tread,'INT')
 if(tread==1) then
   dtset%w90prtunk=intarr(1)
   if ( dtset%w90prtunk < 0 ) then
     write(msg, '(4a)' )&
      'w90prtunk should be greater or equal to zero, however, it was ',dtset%w90prtunk,ch10,&
      'Action: check the values of w90prtunk.'
     MSG_ERROR(msg)
   end if
 end if

 ! Wannier90 - GW quasiparticle interface
 if(dtset%prtwant==3) then
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'getqps',tread,'INT')
   if(tread==1) dtset%getqps=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'irdqps',tread,'INT')
   if(tread==1) dtset%irdqps=intarr(1)
 end if

  !van der Waals with DFT-D approach
 if(dtset%vdw_xc==5.or.dtset%vdw_xc==6.or.dtset%vdw_xc==7) then
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_tol',tread,'DPR')
   if(tread==1) dtset%vdw_tol=dprarr(1)
   if (dtset%vdw_xc==6.or.dtset%vdw_xc==7) then
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_tol_3bt',tread,'DPR')
     if(tread==1) dtset%vdw_tol_3bt=dprarr(1)
   end if
 end if

 ! van der Waals with Wannier functions (Silvestrelli's approach versions 1,2 and 3)
 if(dtset%vdw_xc==10.or.dtset%vdw_xc==11.or.dtset%vdw_xc==12.or.dtset%vdw_xc==14) then
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vdw_nfrag',tread,'INT')
   if(tread==1) dtset%vdw_nfrag=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,dtset%natom,string(1:lenstr),'vdw_typfrag',tread,'INT')
   if(tread==1)then
     dtset%vdw_typfrag(1:dtset%natom)=intarr(1:dtset%natom)
     do ii=1,dtset%natom
       if(dtset%vdw_typfrag(ii)<0)then
         write(msg, '(a,a,a,i0,a,i0,a,a)' )&
          'All the components of vdw_typfrag must be greater than 0.',ch10,&
          'The component',ii,' is equal to ',dtset%vdw_typfrag(ii),ch10,&
          'Action: check the values of vdw_typfrag.'
         MSG_ERROR(msg)
       end if
     end do
   end if
 end if
 if(dtset%vdw_xc==10.or.dtset%vdw_xc==11.or.dtset%vdw_xc==12.or.dtset%vdw_xc==14) then
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'vdw_supercell',tread,'INT')
   if(tread==1)then
     dtset%vdw_supercell(1:3)=intarr(1:3)
     if (dtset%vdw_supercell(1)<zero.and.dtset%vdw_supercell(2)<zero) then
       write(msg, '(7a)' )&
        ' only one component of vdw_supercell could be < 0, however, it was ',ch10,&
        dtset%vdw_supercell(1),dtset%vdw_supercell(2),dtset%vdw_supercell(3),ch10,&
        'Action: check the components of vdw_supercell.'
       MSG_ERROR(msg)
     end if
     if (dtset%vdw_supercell(2)<zero.and.dtset%vdw_supercell(3)<zero) then
       write(msg, '(7a)' )&
        'only one component of vdw_supercell could be < 0, however, it was ',ch10,&
        dtset%vdw_supercell(1),dtset%vdw_supercell(2),dtset%vdw_supercell(3),ch10,&
        'Action: check the components of vdw_supercell.'
       MSG_ERROR(msg)
     end if
     if (dtset%vdw_supercell(1)<zero.and.dtset%vdw_supercell(3)<zero) then
       write(msg, '(7a)' )&
        'only one component of vdw_supercell could be < 0, however, it was ',ch10,&
        dtset%vdw_supercell(1),dtset%vdw_supercell(2),dtset%vdw_supercell(3),ch10,&
        'Action: check the components of vdw_supercell.'
       MSG_ERROR(msg)
     end if
   end if

 end if !vdw_xc==10 or vdw_xc==11 or vdw_xc==12 or vdw_xc==14

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ucrpa',tread,'INT')
 if(tread==1) dtset%ucrpa=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'ucrpa_bands',tread_key,'INT')
 if(tread_key==1) then
   if(dtset%nsym/=1) then
     write(msg,"(a,i0,3a)")&
      'Value of nsym different from 1 when ucrpa_bands is used is under test ',dtset%nsym,&
      ' (because symmetry is not yet used)',ch10,&
      'Action: check your calculation  with nsym=1'
     MSG_WARNING(msg)
   end if
   dtset%ucrpa_bands(1:2)=intarr(1:2)
 end if

 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'ucrpa_window',tread,'ENE')
 if(tread==1) then
   if(dtset%nsym/=1) then
     write(msg,*)&
      'Value of nsym different from 1 when ucrpa_windows is used is under test ',dtset%nsym,&
      ' (because symmetry is not yet used)',ch10,&
      'Action: check your calculation  with nsym=1'
     MSG_WARNING(msg)
   end if
   dtset%ucrpa_window(1:2)=dprarr(1:2)
 end if

 if(tread==1.and.tread_key==1) then
   write(msg, '(a,a,a,a,a)' )&
    'ucrpa_bands and ucrpa_window cannot be specified simultaneously',ch10,&
    'for the same dataset.',ch10,&
    'Action: check the input file.'
   MSG_ERROR(msg)
 end if

 ! Tolerance variables
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'tolimg',tread,'ENE')
 if(tread==1) dtset%tolimg=dprarr(1)

 ! find which tolerance is used for the ionic relaxations
 tolmxde_=dtset%tolmxde
 tolmxf_=dtset%tolmxf
 itol=0
 itol_gen=0

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'tolmxde',tread,'ENE',ds_input)
 if(tread==1) then
   if (ds_input == 0) then
     tolmxde_=dprarr(1)
     if(abs(dprarr(1))>tiny(0._dp))itol_gen=itol_gen+1
   else
     dtset%tolmxde=dprarr(1)
     if(abs(dprarr(1))>tiny(0._dp))itol=itol+1
   end if
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'tolmxf',tread,'DPR',ds_input)
 if(tread==1) then
   if (ds_input == 0) then
     tolmxf_=dprarr(1)
     if(abs(dprarr(1))>tiny(0._dp))itol_gen=itol_gen+1
   else
     dtset%tolmxf=dprarr(1)
     if(abs(dprarr(1))>tiny(0._dp))itol=itol+1
   end if
 else
   itol_gen=itol_gen+1
 end if

 ! check for definitions of tolmxf and tolmxde for the present dataset
 if (itol > 1 .or. itol_gen > 1) then
   write(msg, '(5a)' )&
    'Only one of the tolmxf/tolmxde variables may be defined at once.',ch10,&
    'Action: check values of tolmxf, tolmxde. If you want to use ',ch10,&
    'tolmxde, you should explicitly put tolmxf to 0.0.'
   MSG_ERROR(msg)
 end if

 ! if no value is given for jdtset, use defaults
 if (itol == 0 .and. itol_gen == 1) then
   dtset%tolmxde=tolmxde_
   dtset%tolmxf=tolmxf_
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'focktoldfe',tread,'DPR')
 if(tread==1) dtset%focktoldfe=dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'tolrde',tread,'DPR')
 if(tread==1) dtset%tolrde=dprarr(1)

 ! find which tolXXX are defined generically and for this jdtset
 tolwfr_=zero
 toldfe_=zero
 toldff_=zero
 tolrff_=zero
 tolvrs_=zero
 itol=0
 itol_gen=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'tolwfr',tread,'DPR',ds_input)
 if(tread==1) then
   if (ds_input == 0) then
     tolwfr_=dprarr(1)
     if(abs(dprarr(1))>tiny(0._dp))itol_gen=itol_gen+1
   else
     dtset%tolwfr=dprarr(1)
     if(abs(dprarr(1))>tiny(0._dp))itol=itol+1
   end if
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'toldff',tread,'DPR',ds_input)
 if(tread==1) then
   if (ds_input == 0) then
     toldff_=dprarr(1)
     if(abs(dprarr(1))>tiny(0._dp))itol_gen=itol_gen+1
   else
     dtset%toldff=dprarr(1)
     if(abs(dprarr(1))>tiny(0._dp))itol=itol+1
   end if
   dtset%optforces=1
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'tolrff',tread,'DPR',ds_input)
 if(tread==1) then
   if (ds_input == 0) then
     tolrff_=dprarr(1)
     if(abs(dprarr(1))>tiny(0._dp))itol_gen=itol_gen+1
   else
     dtset%tolrff=dprarr(1)
     if(abs(dprarr(1))>tiny(0._dp))itol=itol+1
   end if
   dtset%optforces=1
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'toldfe',tread,'ENE',ds_input)
 if(tread==1) then
   if (ds_input == 0) then
     toldfe_=dprarr(1)
     if(abs(dprarr(1))>tiny(0._dp))itol_gen=itol_gen+1
   else
     dtset%toldfe=dprarr(1)
     if(abs(dprarr(1))>tiny(0._dp))itol=itol+1
   end if
 end if

 if((dtset%accuracy>0).and.(itol>0.or.itol_gen>0)) tolvrs_=zero

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'tolvrs',tread,'DPR',ds_input)
 if(tread==1) then
   if (ds_input == 0) then
     tolvrs_=dprarr(1)
     if(abs(dprarr(1))>tiny(0._dp))itol_gen=itol_gen+1
   else
     dtset%tolvrs=dprarr(1)
     if(abs(dprarr(1))>tiny(0._dp))itol=itol+1
   end if
 end if

 ! check for multiple definitions of tolXXX for the present dataset
 if (itol > 1 .or. itol_gen > 1) then
   write(msg, '(3a)' )&
   'Only one of the tolXXX variables may be defined at once.',ch10,&
   'Action: check values of tolvrs, toldfe, tolrff, tolwfr, and toldff.'
   MSG_ERROR(msg)
 end if

 ! if no value is given for jdtset, use defaults
 if (itol == 0 .and. itol_gen == 1) then
   dtset%tolwfr=tolwfr_
   dtset%toldfe=toldfe_
   dtset%toldff=toldff_
   dtset%tolrff=tolrff_
   dtset%tolvrs=tolvrs_
 end if

 ! Tolerance variable for TFW initialization step
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'tfw_toldfe',tread,'ENE',ds_input)
 if(tread==1) then
   dtset%tfw_toldfe=dprarr(1)
 else if (dtset%tfkinfunc>10.and.dtset%toldfe>tiny(0._dp)) then
   dtset%tfw_toldfe=dtset%toldfe
 end if

 ! Tolerance variables for electrons-positron tc-dft
 toldfe_=zero;toldff_=zero
 itol=0;itol_gen=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'postoldff',tread,'DPR',ds_input)
 if(tread==1) then
   if (ds_input == 0) then
     toldff_=dprarr(1)
     if(abs(dprarr(1))>tiny(0._dp))itol_gen=itol_gen+1
   else
     dtset%postoldff=dprarr(1)
     if(abs(dprarr(1))>tiny(0._dp))itol=itol+1
   end if
   dtset%postoldfe=zero
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'postoldfe',tread,'ENE',ds_input)
 if(tread==1) then
   if (ds_input == 0) then
     toldfe_=dprarr(1)
     if(abs(dprarr(1))>tiny(0._dp))itol_gen=itol_gen+1
   else
     dtset%postoldfe=dprarr(1)
     if(abs(dprarr(1))>tiny(0._dp))itol=itol+1
   end if
 end if
 if (itol > 1.or.itol_gen >1) then
   write(msg, '(3a)' )&
   'Only one of the postolXXX variables may be defined at once.',ch10,&
   'Action: check values of postoldfe and postoldff.'
   MSG_ERROR(msg)
 end if
 if (itol==0.and.itol_gen==1) then
   dtset%postoldfe=toldfe_
   dtset%postoldff=toldff_
 end if

 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'vprtrb',tread,'ENE')
 if(tread==1) dtset%vprtrb(:)=dprarr(1:2)

!TODO: should this test be just on natsph without caring about the others? natsph might be used for other purposes
 if((dtset%pawfatbnd>0.or.dtset%prtdos==3) .and. dtset%natsph>0)then
   call intagm(dprarr,intarr,jdtset,marr,dtset%natsph,string(1:lenstr),'iatsph',tread,'INT')
   if(tread==1) dtset%iatsph(1:dtset%natsph)=intarr(1:dtset%natsph)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtdensph',tread,'INT')
 if(tread==1) then
   dtset%prtdensph=intarr(1)
 else if (dtset%nsppol==1.and.dtset%nspden==2) then
   dtset%prtdensph=1
 end if

 ! Initialize atvshift
 if(dtset%natvshift>0)then
   call intagm(dprarr,intarr,jdtset,marr,dtset%natvshift*nsppol*dtset%natom,string(1:lenstr),'atvshift',tread,'ENE')
   if(tread==1) dtset%atvshift(1:dtset%natvshift,1:nsppol,1:dtset%natom)=&
&   reshape(dprarr(1:dtset%natvshift*nsppol*dtset%natom),(/dtset%natvshift,nsppol,dtset%natom/))
 end if

 ! Initialize wtatcon
 if(dtset%nconeq>0)then

   ! Read and check natcon
   ABI_ALLOCATE(natcon,(dtset%nconeq))
   call intagm(dprarr,intarr,jdtset,marr,dtset%nconeq,string(1:lenstr),'natcon',tread,'INT')
   if(tread==1)then
     natcon(:)=intarr(1:dtset%nconeq)
   else
     write(msg, '(3a)' )&
      'When nconeq is positive, natcon MUST be defined.',ch10,&
      'Action: check the values of nconeq and natcon.'
     MSG_ERROR(msg)
   end if
   do ii=1,dtset%nconeq
     if(natcon(ii)<0)then
       write(msg, '(3a,i0,a,i0,2a)' )&
        'All the components of natcon must be greater than 0.',ch10,&
        'The component',ii,' is equal to ',natcon(ii),ch10,&
        'Action: check the values of natcon.'
       MSG_ERROR(msg)
     end if
   end do
   niatcon=sum(natcon(:))

   ! Read and check iatcon
   ABI_ALLOCATE(iatcon,(niatcon))
   call intagm(dprarr,intarr,jdtset,marr,niatcon,string(1:lenstr),'iatcon',tread,'INT')
   if(tread==1)then
     iatcon(:)=intarr(1:niatcon)
   else
     write(msg, '(3a)' )&
     'When nconeq is positive, natcon MUST be defined.',ch10,&
     'Action: check the values of nconeq and natcon.'
     MSG_ERROR(msg)
   end if
   do ii=1,niatcon
     if(iatcon(ii)<0)then
       write(msg, '(a,a,a,i4,a,i4,a,a)' )&
        'All the components of iatcon must be greater than 0.',ch10,&
        'The component',ii,' is equal to ',iatcon(ii),ch10,&
        'Action: check the values of iatcon.'
       MSG_ERROR(msg)
     end if
   end do

   ! Read wtatcon, and unfold it.
   call intagm(dprarr,intarr,jdtset,marr,3*niatcon,string(1:lenstr),'wtatcon',tread,'DPR')
   if(tread/=1)then
     write(msg, '(3a)' )&
     'When nconeq is positive, wtatcon MUST be defined.',ch10,&
     'Action: check the values of nconeq and wtatcon.'
     MSG_ERROR(msg)
   end if
   iat=0
   do ii=1,dtset%nconeq
     do jj=1,natcon(ii)
       dtset%wtatcon(1:3,iatcon(jj+iat),ii)=dprarr(1+3*(jj+iat-1):3+3*(jj+iat-1))
     end do
     iat=iat+natcon(ii)
   end do

   ABI_DEALLOCATE(iatcon)
   ABI_DEALLOCATE(natcon)
 end if

  ! Initialize chempot
 if(dtset%nzchempot>0)then
   call intagm(dprarr,intarr,jdtset,marr,3*dtset%nzchempot*ntypat,string(1:lenstr),'chempot',tread,'DPR')
   if(tread==1) dtset%chempot(1:3,1:dtset%nzchempot,1:ntypat)=&
&   reshape(dprarr(1:3*dtset%nzchempot*ntypat),(/3,dtset%nzchempot,ntypat/))
 end if

 ! Initialize the list of k points, as well as wtk and istwfk
 intimage=1 ; if(nimage>2)intimage=(1+nimage)/2
 call invacuum(jdtset,lenstr,natom,dtset%rprimd_orig(1:3,1:3,intimage),string,vacuum,dtset%xred_orig(1:3,1:natom,intimage))

 ! reset rfdir to 1 1 1 in case it was left at the default of 0 0 0 for berryopt -1
 if( (dtset%berryopt == -1) .AND. &
& (dtset%rfdir(1) == 0) .AND. &
& (dtset%rfdir(2) == 0) .AND. &
& (dtset%rfdir(3) == 0) ) dtset%rfdir(1:3) = 1

!In case of a Berryphase calculation put response = 1
!in order to set istwfk = 1 at all k-points
 if ( dtset%berryopt/=0 ) then
   response = 1
 end if

!In case of longwave calculation put response = 1
!in order to set istwfk = 1 at all k-points
 if ( dtset%optdriver==RUNL_LONGWAVE ) then
   response = 1
 end if

 nsym=dtset%nsym
 ii=0;if (mod(dtset%wfoptalg,10)==4) ii=2
 if(dtset%ngfft(7)==314)ii=1
 if(dtset%usefock==1.and.dtset%optdriver/=RUNL_SIGMA.and.mod(dtset%wfoptalg,10)/=5) ii=1

 call inkpts(bravais,dtset%chksymbreak,dtset%fockdownsampling,iout,iscf,dtset%istwfk(1:nkpt),jdtset,&
   dtset%kpt(:,1:nkpt),dtset%kptns_hf(:,1:nkpthf),kptopt,dtset%kptnrm,&
   dtset%kptrlatt_orig,dtset%kptrlatt,kptrlen,lenstr,nsym, dtset%getkerange_filepath, &
   nkpt,nkpthf,nqpt,dtset%ngkpt,dtset%nshiftk,dtset%nshiftk_orig,dtset%shiftk_orig,nsym,&
   occopt,dtset%qptn,response,dtset%rprimd_orig(1:3,1:3,intimage),dtset%shiftk,string,&
   dtset%symafm(1:nsym),dtset%symrel(:,:,1:nsym),vacuum,dtset%wtk(1:nkpt), comm, impose_istwf_1=ii)

 dtset%kptrlen=kptrlen

 dtset%kptns(:,1:nkpt)=dtset%kpt(:,1:nkpt)/dtset%kptnrm
 if(nqpt>=1 .and. dtset%optdriver/=RUNL_RESPFN)then
   dtset%kptns(1,1:nkpt)=dtset%kptns(1,1:nkpt)+dtset%qptn(1)
   dtset%kptns(2,1:nkpt)=dtset%kptns(2,1:nkpt)+dtset%qptn(2)
   dtset%kptns(3,1:nkpt)=dtset%kptns(3,1:nkpt)+dtset%qptn(3)
 end if

 if(nkpthf/=0)then
   dtset%kptns_hf(:,1:nkpthf)=dtset%kptns_hf(:,1:nkpthf)/dtset%kptnrm
 end if

 ! Read variables defining the k-path
 ! If kptopt < 0  --> Band structure and kptbounds size is given by abs(kptopt)
 ! If kptopt >= 0 --> We may have a k-path specified by nkpath and kptbounds (used by post-processing tools)
 ! TODO: ndivk?

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ndivsm',tread,'INT')
 if (tread == 1) dtset%ndivsm = intarr(1)

 if (dtset%kptopt < 0) then
   dtset%nkpath = abs(dtset%kptopt)
 else
   dtset%nkpath = 0
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nkpath',tread,'INT')
   if (tread==1) dtset%nkpath = intarr(1)
 end if

 if (dtset%nkpath /= 0) then
   call intagm(dprarr,intarr,jdtset,marr,3*dtset%nkpath,string(1:lenstr),'kptbounds',tread,'DPR')

   if (tread == 1) then
     ABI_MALLOC(dtset%kptbounds, (3, dtset%nkpath))
     dtset%kptbounds = reshape(dprarr(1:3*dtset%nkpath), [3, dtset%nkpath])
   else
     MSG_ERROR("When nkpath /= 0 or kptopt < 0, kptbounds must be defined in the input file.")
   end if
 else
   ABI_MALLOC(dtset%kptbounds, (0,0))
 end if

 ! if prtkpt==-2, write the k-points in netcdf format and exit here so that AbiPy can read the data.
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtkpt',tread,'INT')
 if (tread == 1 .and. intarr(1) == -2) then
#ifdef HAVE_NETCDF
   ncerr= nctk_write_ibz("kpts.nc", dtset%kptns(:,1:nkpt), dtset%wtk(1:nkpt))
   NCF_CHECK(ncerr)
#endif
   MSG_ERROR_NODUMP("kpts.nc file written. Aborting now")
 end if

 if (dtset%nkptgw>0) then
   ! Read bdgw.
   call intagm(dprarr,intarr,jdtset,marr,2*dtset%nkptgw*dtset%nsppol,string(1:lenstr),'bdgw',tread,'INT')
   if(tread==1) then
     dtset%bdgw(1:2,1:dtset%nkptgw,1:dtset%nsppol) =  &
       reshape(intarr(1:2*dtset%nkptgw*dtset%nsppol),[2,dtset%nkptgw,dtset%nsppol])
   end if

   ! Test bdgw values.
   if (dtset%optdriver == RUNL_SIGMA) then  
     if (any(dtset%bdgw(1:2,1:dtset%nkptgw,1:dtset%nsppol) <= 0)) then
       MSG_ERROR("bdgw entries cannot be <= 0. Check input file")
     end if
     if (any(dtset%bdgw(1,1:dtset%nkptgw,1:dtset%nsppol) > dtset%bdgw(2,1:dtset%nkptgw,1:dtset%nsppol))) then
       MSG_ERROR("First band index in bdgw cannot be greater than the second index")
     end if
   end if

   call intagm(dprarr,intarr,jdtset,marr,3*dtset%nkptgw,string(1:lenstr),'kptgw',tread,'DPR')
   if(tread==1) dtset%kptgw(1:3,1:dtset%nkptgw) = reshape(dprarr(1:3*dtset%nkptgw), [3, dtset%nkptgw])
 end if

 if (dtset%nqptdm > 0) then
   call intagm(dprarr,intarr,jdtset,marr,3*dtset%nqptdm,string(1:lenstr),'qptdm',tread,'DPR')
   if(tread==1) dtset%qptdm(1:3,1:dtset%nqptdm)=reshape(dprarr(1:3*dtset%nqptdm), [3, dtset%nqptdm])
 end if

 if(dtset%npvel>0) then
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'pvelmax',tread,'DPR')
   if(tread==1) dtset%pvelmax(1:3)=dprarr(1:3)
 end if

 if(dtset%imgmov==6)then
   call intagm(dprarr,intarr,jdtset,marr,dtset%nimage,string(1:lenstr),'mixesimgf',tread,'DPR')
   if(tread==1)dtset%mixesimgf(1:dtset%nimage)=dprarr(1:dtset%nimage)
 endif

 ! Read remaining input variables depending on images : occ_orig, upawu, jpawu.
 ! MT oct 15: these lines are inspired from invars1 but are not really understood
 intimage=2 ; if(dtset%nimage==1)intimage=1
 do ii=1,dtset%nimage+1
   iimage=ii
   if(dtset%nimage==1 .and. ii==2)exit
   if(dtset%nimage==2 .and. ii==3)exit
   if(dtset%nimage> 2 .and. ii==intimage)cycle ! Will do the intermediate reference image at the last reading
   if(dtset%nimage>=2 .and. ii==dtset%nimage+1)iimage=intimage

   ! Only read occ if (iscf >0 or iscf=-1 or iscf=-3) and (occopt==0 or occopt==2)
   if  (iscf>=0.or.iscf==-1.or.iscf==-3)  then
     if (occopt==2 .and. getocc==0) then
       ! Read occ(nband(kpt)*nkpt*nsppol) explicitly
       call wrtout(std_out,' invars2: reading occ(nband*nkpt*nsppol) explicitly','COLL')
       call intagm(dprarr,intarr,jdtset,marr,bantot,string(1:lenstr),'occ',tread,'DPR')
       if(tread==1) dtset%occ_orig(1:bantot,iimage)=dprarr(1:bantot)
       call intagm_img(dprarr(1:bantot),iimage,jdtset,lenstr,dtset%nimage,bantot,string,'occ',tread_alt,'DPR')
       if(tread_alt==1) dtset%occ_orig(1:bantot,iimage)=dprarr(1:bantot)
     else if(occopt==0) then
       nband1=dtset%nband(1)
       ! Read usual occupancy--same for all k points but might differ for spins
       call intagm(dprarr,intarr,jdtset,marr,nband1*nsppol,string(1:lenstr),'occ',tread,'DPR')
       if(tread==1) dtset%occ_orig(1:nband1*nsppol,iimage)=dprarr(1:nband1*nsppol)
       call intagm_img(dprarr(1:nband1*nsppol),iimage,jdtset,lenstr,dtset%nimage,nband1*nsppol,string,'occ',tread_alt,'DPR')
       if(tread_alt==1) dtset%occ_orig(1:nband1*nsppol,iimage)=dprarr(1:nband1*nsppol)
       ! Fill in full occ array using input values for each k and spin
       ! (make a separate copy for each k point and spin)
       if(nkpt>1)then
         do isppol=nsppol,1,-1
           do ikpt=2,nkpt
             dtset%occ_orig(1+(ikpt-1)*nband1+nkpt*nband1*(isppol-1):ikpt*nband1+nkpt*nband1*(isppol-1),iimage)=&
&             dtset%occ_orig(1+nband1*(isppol-1):nband1*isppol,iimage)
           end do
         end do
       endif
     end if
   end if

   if (dtset%usepawu/=0.or.dtset%usedmft>0) then

     dprarr(1:ntypat)=dtset%upawu(1:ntypat,iimage)
     call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),'upawu',tread,'ENE')
     if(tread==1) dtset%upawu(1:ntypat,iimage)=dprarr(1:ntypat)
     call intagm_img(dprarr(1:ntypat),iimage,jdtset,lenstr,dtset%nimage,ntypat,string,'upawu',tread_alt,'ENE')
     if(tread_alt==1) dtset%upawu(1:ntypat,iimage)=dprarr(1:ntypat)

     dprarr(1:ntypat)=dtset%jpawu(1:ntypat,iimage)
     call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),'jpawu',tread,'ENE')
     if(tread==1) dtset%jpawu(1:ntypat,iimage)=dprarr(1:ntypat)
     call intagm_img(dprarr(1:ntypat),iimage,jdtset,lenstr,dtset%nimage,ntypat,string,'jpawu',tread_alt,'ENE')
     if(tread_alt==1) dtset%jpawu(1:ntypat,iimage)=dprarr(1:ntypat)

     if (dtset%usedmatpu/=0.and.dmatsize>0) then
       nsp=nsppol*nspinor
       call intagm(dprarr,intarr,jdtset,marr,dmatsize,string(1:lenstr),'dmatpawu',tread,'DPR')
       if(tread==1) then
         iat=1;jj=1
         do iatom=1,natom
           lpawu=dtset%lpawu(dtset%typat(iatom))
           if (lpawu/=-1) then
             isiz=nsp*(2*lpawu+1)**2
             dtset%dmatpawu(1:2*lpawu+1,1:2*lpawu+1,1:nsp,iat,iimage)= &
&             reshape(dprarr(jj:jj+isiz),(/2*lpawu+1,2*lpawu+1,nsp/))
             iat=iat+1;jj=jj+isiz
           end if
         end do
       end if

       ABI_ALLOCATE(dmatpawu_tmp,(dmatsize))
       iat=1;jj=1
       do iatom=1,natom
         lpawu=dtset%lpawu(dtset%typat(iatom))
         if (lpawu/=-1) then
           isiz=nsp*(2*lpawu+1)**2
           dmatpawu_tmp(jj:jj+isiz-1)= &
&           reshape(dtset%dmatpawu(1:2*lpawu+1,1:2*lpawu+1,1:nsp,iat,iimage),(/isiz/))
           iat=iat+1;jj=jj+isiz
         end if
       end do
       call intagm_img(dmatpawu_tmp,iimage,jdtset,lenstr,dtset%nimage,dmatsize, string,'dmatpawu',tread_alt,'DPR')
       if(tread_alt==1) then
         iat=1;jj=1
         do iatom=1,natom
           lpawu=dtset%lpawu(dtset%typat(iatom))
           if (lpawu/=-1) then
             isiz=nsp*(2*lpawu+1)**2
             dtset%dmatpawu(1:2*lpawu+1,1:2*lpawu+1,1:nsp,iat,iimage)= &
&             reshape(dmatpawu_tmp(jj:jj+isiz-1),(/2*lpawu+1,2*lpawu+1,nsp/))
             iat=iat+1;jj=jj+isiz
           end if
         end do
       end if
       ABI_DEALLOCATE(dmatpawu_tmp)

       if (tread/=1.and.tread_alt/=1) then
         write(msg, '(3a)' )&
         'When LDA/GGA+U is activated and usedmatpu/=0, dmatpawu MUST be defined.',ch10,&
         'Action: add dmatpawu keyword in input file.'
         MSG_ERROR(msg)
       end if
     end if

   end if

 end do

 if(dtset%jellslab/=0)then
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'slabwsrad',tread,'LEN')
   if(tread==1) dtset%slabwsrad=dprarr(1)

   ! Update number of electrons taking into account jellium
   ! Suppose that the cell has z axis perpendicular to x and y axes.
   ! This will be checked later
   areaxy=abs(dtset%rprimd_orig(1,1,1)*dtset%rprimd_orig(2,2,1)-dtset%rprimd_orig(1,2,1)*dtset%rprimd_orig(2,1,1))
   rhoavg=three/(four_pi*dtset%slabwsrad**3)
   nelectjell=areaxy*(dtset%slabzend-dtset%slabzbeg)*rhoavg
   charge=charge-nelectjell
 end if

 ! Initialize occ if occopt==1 or 3 ... 8,
 ! while if getocc/=0, make a fake initialization
 ! If iscf>0, check the charge of the system, and compute nelect.
 occopt_tmp=occopt
 if(getocc/=0)occopt_tmp=1
 call dtset%chkneu(charge, occopt_tmp)

 ! Now that the occupation numbers have been initialized, can meaningfully define nbandhf.
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nbandhf',tread,'INT')
 if(tread==1) then
   dtset%nbandhf=intarr(1)
 else
   ! If the occupation numbers might change, must keep the maximum number of bands
   if(occopt>=3 .and. occopt<=8)then
     dtset%nbandhf=maxval(dtset%nband(1:nkpt*nsppol))
   else if(occopt==0 .or. occopt==1 .or. occopt==2) then
     ! Eliminate all the bands that are never occupied
     nband1=0 ; bantot=0
     do isppol=1,dtset%nsppol
       do ikpt=1,dtset%nkpt
         do iband=1,dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
           bantot=bantot+1
           if(maxval(abs(dtset%occ_orig(bantot,:)))>tol8)then
             if(iband>nband1)nband1=iband
           end if
         end do
       end do
     end do
     dtset%nbandhf=nband1
   else
     write(msg, '(a,i0,3a)' )'occopt=',occopt,' not allowed.',ch10,'Action: correct your input file.'
     MSG_ERROR(msg)
   end if
 end if

 ! Initialize Berry phase vectors
 ! Should check that nberry is smaller than 20
 if(berryopt>0 .and. nberry>0)then
   call intagm(dprarr,intarr,jdtset,marr,3*nberry,string(1:lenstr),'kberry',tread,'INT')
   if(tread==1) dtset%kberry(1:3,1:nberry)=reshape(intarr(1:3*nberry), [3,nberry])
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ddamp',tread,'DPR')
 if (tread==1) dtset%ddamp = dprarr(1)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'polcen',tread,'DPR')
 if (tread==1) dtset%polcen(1:3) = dprarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),'densty',tread,'DPR')
 if(tread==1) dtset%densty(1:ntypat,1)=dprarr(1:ntypat)

 call intagm(dprarr,intarr,jdtset,marr,npsp,string(1:lenstr),'so_psp',tread,'INT')
 if (tread==1.and.dtset%usepaw==0) dtset%so_psp(1:npsp) = intarr(1:npsp)

 ! LOTF variables
#if defined HAVE_LOTF
 if(dtset%ionmov==23) then
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'lotf_classic',tread,'INT')
   if(tread==1) dtset%lotf_classic=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'lotf_nitex',tread,'INT')
   if(tread==1) dtset%lotf_nitex=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'lotf_nneigx',tread,'INT')
   if(tread==1) dtset%lotf_nneigx=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'lotf_version',tread,'INT')
   if(tread==1) dtset%lotf_version=intarr(1)
 end if
#endif

!Read input variables related to Projected Local Orbitals Wannier functions (plowan)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'plowan_bandi',tread,'INT')
 if(tread==1) dtset%plowan_bandi=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'plowan_bandf',tread,'INT')
 if(tread==1) dtset%plowan_bandf=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'plowan_realspace',tread,'INT')
 if(tread==1) dtset%plowan_realspace=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'plowan_natom',tread,'INT')
 if(tread==1) dtset%plowan_natom=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'plowan_nt',tread,'INT')
 if(tread==1) dtset%plowan_nt=intarr(1)

! if (dtset%plowan_compute>0) then
   call intagm(dprarr,intarr,jdtset,marr,3*dtset%plowan_nt,string(1:lenstr),'plowan_it',tread,'INT')
   if(tread==1) dtset%plowan_it(1:3*dtset%plowan_nt)=intarr(1:3*dtset%plowan_nt)

   call intagm(dprarr,intarr,jdtset,marr,dtset%plowan_natom,string(1:lenstr),'plowan_iatom',tread,'INT')
   if(tread==1) dtset%plowan_iatom(1:dtset%plowan_natom)=intarr(1:dtset%plowan_natom)

   call intagm(dprarr,intarr,jdtset,marr,dtset%plowan_natom,string(1:lenstr),'plowan_nbl',tread,'INT')
   if(tread==1) dtset%plowan_nbl(1:dtset%plowan_natom)=intarr(1:dtset%plowan_natom)
   sumnbl=sum(dtset%plowan_nbl(:))

   call intagm(dprarr,intarr,jdtset,marr,sumnbl,string(1:lenstr),'plowan_lcalc',tread,'INT')
   if(tread==1) dtset%plowan_lcalc(1:sumnbl)=intarr(1:sumnbl)

   call intagm(dprarr,intarr,jdtset,marr,sumnbl,string(1:lenstr),'plowan_projcalc',tread,'INT')
   if(tread==1) dtset%plowan_projcalc(1:sumnbl)=intarr(1:sumnbl)
! end if


   if ((dtset%ucrpa > 0 .and. dtset%plowan_natom == 0).or.(dtset%nbandkss /= 0 .and. dtset%usedmft/=0)) then
     dtset%plowan_natom=1
     dtset%plowan_nbl(:)=1
     dtset%plowan_nt=1
     dtset%plowan_it(:)=0
     dtset%plowan_realspace=1
     do iatom=1,dtset%natom
       lpawu=dtset%lpawu(dtset%typat(iatom))
       if (lpawu/=-1) then
         dtset%plowan_lcalc(:)=lpawu
         dtset%plowan_iatom(:)=iatom
         dtset%plowan_projcalc(:)=-2
       end if
     end do
     dtset%plowan_bandi=dtset%dmftbandi
     dtset%plowan_bandf=dtset%dmftbandf
     if (dtset%nbandkss /= 0 .and. dtset%usedmft/=0) then
       dtset%plowan_compute=1
       dtset%usedmft=0
     else if (dtset%optdriver==3) then
       dtset%plowan_compute=10
     else if(dtset%optdriver==4) then
       dtset%plowan_compute=10
     end if
   end if

 ! band range for self-energy sum
 call intagm(dprarr, intarr, jdtset, marr, 2, string(1:lenstr), 'sigma_bsum_range', tread, 'INT')
 if (tread == 1) then
    dtset%sigma_bsum_range = intarr(1:2)
    ABI_CHECK(all(dtset%sigma_bsum_range > 0), "sigma_bsum_range cannot be negative")
    ABI_CHECK(dtset%sigma_bsum_range(2) >= dtset%sigma_bsum_range(1), "sigma_bsum_range(2) must be >= (1)")
 end if

 ! band range for self-energy corrections.
 call intagm(dprarr, intarr, jdtset, marr, 2, string(1:lenstr), 'sigma_erange', tread, 'ENE')
 if (tread == 1) dtset%sigma_erange = dprarr(1:2)

 ! IBZ k-points for transport calculation in terms of transport_ngkpt
 call intagm(dprarr, intarr, jdtset, marr, 3, string(1:lenstr), 'transport_ngkpt', tread, 'INT')
 if (tread == 1) dtset%transport_ngkpt = intarr(1:3)

 ! IBZ k-points for electron self-energy given in terms of sigma_ngkpt
 call intagm(dprarr, intarr, jdtset, marr, 3, string(1:lenstr), 'sigma_ngkpt', tread, 'INT')

 if (tread == 1) then
   ! sigma_ngkpt mode --> initialize shifts, provide default if not given in input
   ! Consistency check: nkptgw must be zero, sigma_erange should not be given.
   ABI_CHECK(dtset%nkptgw == 0, "nkptgw and sigma_ngkpt are mutually exclusive.")

   dtset%sigma_ngkpt = intarr(1:3)
   call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'sigma_nshiftk', tread, 'INT')
   if (tread == 1) dtset%sigma_nshiftk = intarr(1)
   if (dtset%sigma_nshiftk < 1 .or. dtset%sigma_nshiftk > MAX_NSHIFTK ) then
     write(msg,  '(a,i0,2a,i0,3a)' )&
     'The only allowed values of nshiftk are between 1 and ',MAX_NSHIFTK,ch10,&
     'while it is found to be ',dtset%sigma_nshiftk,'.',ch10,&
     'Action: change the value of sigma_nshiftk in your input file, or change kptopt.'
     MSG_ERROR(msg)
   end if

   call intagm(dprarr, intarr, jdtset, marr, 3*dtset%sigma_nshiftk, string(1:lenstr), 'sigma_shiftk', tread, 'DPR')
   ! Yes, I know that multidatasets will likely crash if sigma_nshiftk changes...
   ABI_CALLOC(dtset%sigma_shiftk, (3, dtset%sigma_nshiftk))
   if (tread == 1) then
     dtset%shiftk(:,1:dtset%sigma_nshiftk) = reshape(dprarr(1:3*dtset%sigma_nshiftk), [3,dtset%sigma_nshiftk])
   end if
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),"wfk_task",tread,'KEY',key_value=key_value)
 if (tread==1) dtset%wfk_task = str2wfktask(tolower(key_value))
 if (dtset%optdriver == RUNL_WFK .and. dtset%wfk_task == WFK_TASK_NONE) then
   MSG_ERROR(sjoin("A valid wfk_task must be specified when optdriver= ", itoa(dtset%optdriver), ", Received:", key_value))
 end if

 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)

 call timab(191,2,tsec)

end subroutine invars2
!!***

end module m_invars2
!!***
