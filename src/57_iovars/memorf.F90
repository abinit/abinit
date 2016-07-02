!{\src2tex{textfont=tt}}
!!****f* ABINIT/memorf
!! NAME
!! memorf
!!
!! FUNCTION
!! Estimation of the memory needed for a response-function job.
!! According to the value of the option variable,
!! might also try to allocate this amount of memory, and if it fails,
!! might estimate the available memory.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex=1 or 2, indicate whether the den and pot functions are real or complex
!!  getcell=if non-zero, the values of acell and rprim are taken from
!!   the output of another dataset
!!  idtset=number of the current dataset
!!  intxc=control xc quadrature
!!  iout=unit number for output of formatted data.
!!  iprcel=govern the choice of preconditioner for the SCF cycle
!!  iscf=governs the choice of SCF algorithm, or non-SCF calculation.
!!  jdtset=index of the current dataset
!!  lmnmax=max. number of (l,m,n) components over all type of psps
!!  lnmax =max. number of (l,n)   components over all type of psps
!!  mband =maximum number of bands
!!  mffmem =governs the number of FFT arrays which are fit in core memory
!!  mgfft =maximum single fft dimension
!!  mkmems=number of k points which can fit in memory; set to 0 if use disk
!!    the three values correspond to mkmem, mkqmem and mk1mem
!!  mpi_enreg=informations about MPI parallelization
!!  mpssang is 1+maximum angular momentum for nonlocal pseudopotential
!!  mpssoang is 1+maximum (spin*angular momentum) for nonlocal pseudopotential
!!  mpw   =maximum number of planewaves in basis sphere (large number)
!!  mqgrid=maximum dimension of grid of q values for psp representations
!!  natom =number of atoms in unit cell
!!  nband(nkpt*nsppol)=number of bands at each k point, for each polarization
!!  nfft  =(effective) number of FFT grid points (for one processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt  =number of k points
!!  nloalg(3)=governs the choice of the algorithm for non-local operator.
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=number of channels for spin-polarization (1 or 2)
!!  nsym  =number of symmetry elements in space group
!!  ntypat=number of types of atoms
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  occopt=option for occupation numbers. If 3<=occopt<=7, varying occupation
!!  optddk=1 if ddk is computed during run
!!  optphon=1 if phonons are computed during run
!!  option : if 0 , no test of available memory
!!           if 1 , the routine tries to allocate the estimated memory, for testing
!!                    purposes, and if a failure occurs, the routine stops.
!!           if 2 , like 1, but before stopping, the routine will provide
!!                    an estimation of the available memory.
!!  optstrs=1 if strain perturbation is computing during run
!!  prtvol=control print volume
!!  useylm=governs the way the nonlocal operator is to be applied:
!!         1=using Ylm, 0=using Legendre polynomials
!!  use_gpu_cuda=1 if Cuda (GPU) is on
!!  xclevel= level of the XC functional
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! for the estimation, it is only taken into account those
!! arrays that have some probability of being larger than 1000*8 bytes :
!! - All the arrays that have large numbers as one of their dimensions
!! (mqgrid, mpw, nfft, ngfft(4)*ngfft(5)*ngfft(6),n1xccc
!!                                      or a constant larger than 1000)
!! - All the arrays that have a product of two moderately large numbers
!! (potential size above 30  : mband, mgfft, mkmems, natom, nkpt, nsym,
!!  or a constant larger than 30)
!! After this estimation, an amount of (176 + 55 + 6*natom) Kbytes is added
!! to take into account the static arrays declared
!! in rhohxc and daughter routines (at maximum 22*1000 dp numbers),
!! as well as other arrays like
!! character(len=500) :: message (present in about 100 routines), or the different
!! arrays allocated in move.f, brdmin.f, gstate.f (xf array) or pspini.f
!! In the case 3<=occopt<=7 this amount is increased by 760 Kbytes
!! to take into account the arrays smdfun, occfun, entfun, workfun and xgrid,
!! declared in getnel
!!
!! The current version takes into account only :
!! 1) and 2) the "main chain" in its two slightly different versions :
!! driver - respfn - dfpt_looppert - dfpt_scfcv - dfpt_vtorho - dfpt_vtowfk -
!!     dfpt_cgwf - getghc - fourwf or (nonlop+opernl)
!!
!! Also, it is assumed that the potentials are non-local, even if there
!!     are local ! It would be necessary to update this routine
!!     now that the beginning of psp files is read before
!!     the present call (XG 980502)
!!
!! Some BIG approximations, not present in the GS corresponding routine
!!  have been done : nsym=nsym1, nkpt=nkpt_rbz, mpw=mpw1 ...
!!
!! PARENTS
!!      memory_eval
!!
!! CHILDREN
!!      memana,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine memorf(cplex,n1xccc,getcell,idtset,intxc,iout,iprcel,&
& iscf,jdtset,lmnmax,lnmax,mband,mffmem,mgfft,&
& mkmems,mpi_enreg,mpsang,mpssoang,mpw,mqgrid,&
& natom,nband,nfft,ngfft,&
& nkpt,nloalg,nspden,nspinor,nsppol,nsym,ntypat,&
& occopt,optddk,optphon,option,optstrs,prtvol,useylm,use_gpu_cuda,xclevel)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memorf'
 use interfaces_14_hidewrite
 use interfaces_57_iovars, except_this_one => memorf
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,getcell,idtset,intxc,iout,iprcel,iscf
 integer,intent(in) :: jdtset,lmnmax,lnmax,mband,mffmem,mgfft,mpsang
 integer,intent(in) :: mpssoang,mpw,mqgrid,n1xccc,natom,nfft,nkpt
 integer,intent(in) :: nspden,nspinor,nsppol,nsym,ntypat,occopt
 integer,intent(in) :: optddk,option,optphon,optstrs,prtvol,useylm
 integer,intent(in) :: use_gpu_cuda,xclevel
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: mkmems(3),nband(nkpt*nsppol),ngfft(18)
 integer,intent(in) :: nloalg(3)

!Local variables-------------------------------
!marrays= maximal number of arrays to be monitored (or group of arrays)
!cmpw(marrays)=count of blocks of size mpw bytes
!cfft(marrays)=number of blocks of size nfft bytes
!cadd(marrays)=additional storage needed (in bytes)
!dttyp(marrays)=datatype of the array : 4 for integers, 8 for real(dp)
!nchain= number of different chains of routines
!chain(marrays,nchain)=different chains of routines
!scalars
 integer,parameter :: marrays=150,nchain=2
 integer :: fftalgb,matblk,maxmkmem,mincat,mk1mem,mkmem,mkqmem,mu,n_fftgr
 integer :: narr_fourdp,ngrad,nprocwf
 integer :: my_natom
 real(dp) :: mbcg,mbdiskpd,mbdiskwf,mbf_fftgr,mbgylm
 character(len=500) :: message
 character(len=1) :: firstchar
!arrays
 integer :: dttyp(marrays)
 real(dp) :: cadd(marrays),cfft(marrays),cmpw(marrays)
 real(dp),allocatable :: cfft_dum(:)
 logical :: chain(marrays,nchain)

! **************************************************************************

 if(option<0 .or. option>2)then
   write(message, '(a,i0,a)')'option= ',option,' while the only allowed values are 0, 1, or 2.'
   MSG_BUG(message)
 end if

 firstchar=' ';if (use_gpu_cuda==1) firstchar='_'
 cmpw(:)=zero ; cfft(:)=zero ; cadd(:)=zero 
 dttyp(:)=0

 call wrtout(std_out,' memorf : analysis of memory needs ','COLL')

 if(jdtset>=100)then
   write(message,'(80a,a,a,i5,a)')('=',mu=1,80),ch10,&
&   ' Values of the parameters that define the memory need for DATASET',jdtset,&
&   ' (RF).'
 else if(jdtset/=0)then
   write(message,'(80a,a,a,i3,a)')('=',mu=1,80),ch10,&
&   ' Values of the parameters that define the memory need for DATASET',jdtset,&
&   ' (RF).'
 else
   write(message,'(80a,a,a,a)')('=',mu=1,80),ch10,&
&   ' Values of the parameters that define the memory need of the present run',&
&   ' (RF).'
 end if
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 mkmem=mkmems(1)
 mkqmem=mkmems(2)
 mk1mem=mkmems(3)
 my_natom=natom;if (mpi_enreg%nproc_atom>1) my_natom=mpi_enreg%my_natom

 write(message,'( 4(a,i8),a,4(a,i8) )' ) &
& '     intxc =',intxc   ,'      iscf =',iscf,&
& '    lmnmax =',lmnmax  ,'     lnmax =',lnmax,ch10,&
& '     mgfft =',mgfft,'  mpssoang =',mpssoang,&
& '    mqgrid =',mqgrid,'     natom =',natom
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 write(message,'( 4(a,i8),a,4(a,i8),a,4(a,i8) )' ) &
& '  nloc_mem =',nloalg(2)*(nloalg(3)+1),'    nspden =',nspden ,&
& '   nspinor =',nspinor,'    nsppol =',nsppol ,ch10,&
& '      nsym =',nsym,'    n1xccc =',n1xccc ,&
& '    ntypat =',ntypat,'    occopt =',occopt ,ch10,&
& '   xclevel =',xclevel
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 write(message,'(4(3(a,i12),a))') &
& '-    mband =',mband  ,'        mffmem =',mffmem,&
& '         mkmem =',mkmem  ,ch10,&
& '-   mkqmem =',mkqmem ,'        mk1mem =',mk1mem,&
& '           mpw =',mpw  ,ch10,&
& '      nfft =',nfft ,'          nkpt =',nkpt
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 if (my_natom/=natom)then
   write(message,'(a,i10)') 'Pmy_natom=',my_natom
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

 write(message,'(80a)') ('=',mu=1,80)
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 if(getcell>0 .or. (getcell<0 .and. idtset+getcell>0) )then
   write(message,'(a,a,a,a,a,a,i3,a,i3,a,a,a,a,a,a)' )ch10,&
&   ' memorf : COMMENT -',ch10,&
&   '  The determination of memory needs at this stage is meaningless,',ch10,&
&   '  since getcell = ',getcell,' is non-zero, while idtset=',idtset,'.',ch10,&
&   '  The following numbers are obtained by supposing that acell and rprim',ch10,&
&   '  are NOT taken from a previous dataset. You cannot rely on them.',ch10
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

 n_fftgr=1
 if(iscf==1)            n_fftgr=5
 if(iscf==2.or.iscf==3) n_fftgr=4
 if(iscf==5.or.iscf==6) n_fftgr=10

!work1 and work2 in fourdp : take into account approximately fftalgb
 fftalgb=mod(ngfft(7),100)/10
 if(fftalgb==0)narr_fourdp=2*2
 if(fftalgb==1)narr_fourdp=2

 ngrad=1
 if(xclevel==2)ngrad=2

!(0)                     in main, driver, and respfn -------------------
!indsym (respfn)
 cadd(1)=4*nsym*natom          ; dttyp(1)=4
!rhor,rhog (respfn)
 cfft(2)=nspden+2              ; dttyp(2)=8
!occ (driver), doccde (respfn)
 cadd(3)=2*mband*nkpt*nsppol   ; dttyp(3)=8
!qgrid,vlspl,ffspl (driver)
 cadd(4)=mqgrid*(1+2*ntypat*(1+lnmax))   &
& ; dttyp(4)=8
!xccc1d (driver)
 cadd(5)=n1xccc*6*ntypat       ; dttyp(5)=8
!vtrial (respfn)
 cfft(6)=nspden                ; dttyp(6)=8
!kxc (respfn)
 cfft(7)=2*nspden-1            ; dttyp(7)=8

!(1-2)                   in dfpt_looppert --------------------------------------
!ph1d
 cadd(11)=2*3*(2*mgfft+1)*natom ; dttyp(11)=8
!vpsp1
 cfft(12)=cplex                ; dttyp(12)=8
!indsy1  assume that nsym=nsym1
 cadd(13)=4*nsym*natom         ; dttyp(13)=4
!irrzonr1 and phnons1  assume that nsym=nsym1
 if(nsym/=1)then
   cfft(14)=(2+(nspden/4))*((nspden/nsppol)-3*nspden/3)     ; dttyp(14)=4
   cfft(15)=2*((nspden/nsppol)-3*nspden/3)                  ; dttyp(15)=8
 end if
!doccde_rbz, eigen0, eigenq, occ_rbz, docckqde, occkq, resid
!assume than nkpt=nkpt_rbz
 cadd(16)=7*mband*nkpt*nsppol  ; dttyp(16)=8
!kg
 cmpw(18)=3*mkmem              ; dttyp(18)=4
!cg
 cmpw(19)=2*nspinor*mband*mkmem*nsppol  ; dttyp(19)=8
!kg1
 cmpw(21)=3*mk1mem             ; dttyp(21)=4
!cgq
 cmpw(22)=2*nspinor*mband*mkqmem*nsppol  ; dttyp(22)=8
!cg1
 cmpw(23)=2*nspinor*mband*mk1mem*nsppol  ; dttyp(23)=8
!rhor1,rhog1
 cfft(24)=cplex*nspden+2       ; dttyp(24)=8
!eigen1
!assume than nkpt=nkpt_rbz
 cadd(25)=2*mband*mband*nkpt*nsppol      ; dttyp(25)=8
!ylm
 cmpw(26)=mkmem*mpsang*mpsang*useylm     ; dttyp(26)=8

!(3)                     in dfpt_scfcv --------------------------------------

!vhartr1,vtrial1,vxc
 cfft(31)=cplex+cplex*nspden+nspden      ; dttyp(31)=8
 if(iscf>0)then
!  f_fftgr
   cfft(32)=cplex*nspden*n_fftgr*mffmem    ; dttyp(32)=8
 end if

!(4)                   in dfpt_vtorho----------------------------------------

!proc_distrb
 cadd(41)=nkpt*mband*nsppol    ; dttyp(41)=4
!kg_k,kg1_k
 cmpw(42)=6                    ; dttyp(42)=4
!rhoaug1, vlocal, vlocal1
 cfft(43)=2*cplex+1            ; dttyp(43)=8
 cadd(43)=(2*cplex+1)*(ngfft(4)*ngfft(5)*ngfft(6)-nfft)

 if(mkqmem==0)then
!  cgq_disk
   cmpw(45)=2*nspinor*mband      ; dttyp(45)=8
 end if
!doccde_k,doccde_kq,eig0_k, ..., eig1_k, rocceig
 cadd(47)=(14+3*mband)*mband   ; dttyp(47)=8
!ylm_k,ylm1_k
 cmpw(49)=2*mpsang*mpsang*useylm  ; dttyp(49)=8

!(5)                     in dfpt_vtowfk --------------------------------------

!dkinpw,kinpw1
 cmpw(51)=2                    ; dttyp(51)=8
!ffnlk,ffnl1,ffnlkq
 cmpw(52)=2*(ntypat+2)*lmnmax  ; dttyp(52)=8
!ghc,gvnlc,gvnl1
 cmpw(53)=6*nspinor            ; dttyp(53)=8
!ph3d
 matblk=NLO_MINCAT
 if(nloalg(2)<=0)matblk=natom
 cmpw(54)=2*matblk             ; dttyp(54)=8
!wfraug,wfraug1,rhoaug
 cfft(55)=5                    ; dttyp(55)=8
 cadd(55)=5*(ngfft(4)*ngfft(5)*ngfft(6)-nfft)
!cwavef,cwave0,cwave1
 cmpw(56)=6*nspinor            ; dttyp(56)=8

!(6)                     in dfpt_cgwf ----------------------------------------

!gh1, gh_direc, gvnl_direc, conjgr, direc, vresid, cwaveq
 cmpw(61)=14*nspinor            ; dttyp(61)=8

!(9a)                    in getghc and fourwf----------------------------

!work (in getghc)
 cfft(91)=2                    ; dttyp(91)=8
 cadd(92)=2*(ngfft(4)*ngfft(5)*ngfft(6)-nfft)
!work1 (in fourwf)
 cfft(92)=2                    ; dttyp(92)=8
 cadd(92)=2*(ngfft(4)*ngfft(5)*ngfft(6)-nfft)

!(9b)                    in getghc, nonlop and opernl--------------------
 mincat=min(NLO_MINCAT,natom-ntypat+1)
 if (useylm==0) then                          ! ===== nonlop_pl
!  gxa  (in nonlop)
   cadd(94)=2*20*mincat*2       ; dttyp(94)=8
!  dgxdt  (in nonlop)
   cadd(95)=2*3*20*mincat*2    ; dttyp(95)=8
!  dgxds  (in nonlop)
   cadd(96)=2*56*mincat*2      ; dttyp(96)=8
!  teffv (in opernl4 - no distinction is made for opernl, opernl2 or opernl3)
!  kpgx, ffkg
!  here, evaluate an upper value, with nproj=2, p,d and f orbitals, but not
!  considering the stress, since it will be called outside of the main chain
   cadd(97)=NLO_MBLKPW*40        ; dttyp(97)=8
!  kpg if nloalg(3)=1
   cadd(98)=3*mpw*nloalg(3)     ; dttyp(98)=8
 else                                        ! ===== nonlop_ylm
!  gx + gxfac
   cadd(94)=2* 2*mpw*lmnmax*mincat    ; dttyp(94)=8
!  dgxdt + dgxdtfac + d2gxdt
   if (optddk>0.and.optphon==0.and.optstrs==0) cadd(95)=2*2*mpw*lmnmax*mincat
   if (optphon>0) cadd(95)=12*2*mpw*lmnmax*mincat
   if (optstrs>0) cadd(95)=72*2*mpw*lmnmax*mincat
   dttyp(95)=8
!  kpg
   cadd(96)=2*3*mpw       ; dttyp(96)=8
   if (optphon>0) cadd(96)=cadd(96)+2*6*mpw
!  miscelaneous: indlmn_typ, ffnl_typ
   cadd(97)=lmnmax*(6+mpw*(2+optstrs)); dttyp(97)=8
!  opernla_ylm: scalar,scali,scalarr,scalari
   cadd(98)=2*mpw+2*mpw
   if (optddk>0.and.optstrs==0) cadd(98)=cadd(98)+2*mpw
   if (optstrs>0) cadd(98)=cadd(98)+9*2*mpw
   dttyp(98)=8
 end if

!--------------------------------------------------------------------------

 chain(:,:)=.true.

!Define the main chain version a (fourwf)
 chain(93:100,1)=.false.

!Define the main chain version b (nonlop+opernl)
 chain(91:92,2)=.false.

!The memory needed for each chain has been computed
!-------------------------------------------------------------------------
!Still need some auxiliary data : estimate the disk space
!or the maximum segment size.

!XG030513 : MPIWF need to multiply mbdiskwf by the number of processors
!in the WF group. For the time being, nprocwf=1
 nprocwf=mpi_enreg%nproc_fft

 mbdiskwf=(8*2*mpw*nprocwf*sum(nband(1:nkpt*nsppol)))/1024._dp**2 + 0.002_dp
 mbdiskpd=(8*nfft*nsppol)/1024._dp**2 + 0.002_dp

!Determine the largest array out of cg,cg1,cgq, cg_disk or f_fftgr (f_fftgr_disk)
 if(mkmem==0 .and. mk1mem==0 .and. mkqmem==0)then
   mbcg=(8*2*mpw*nspinor*mband)/1024._dp**2 + 0.002_dp
 else
   maxmkmem=maxval(mkmems(:))
   mbcg=(8*2*mpw*nspinor*mband*maxmkmem*nsppol)/1024._dp**2 + 0.002_dp
 end if
 if(mffmem==0)then
   mbf_fftgr=(8*cplex*nfft*n_fftgr)/1024._dp**2 + 0.002_dp
 else
   mbf_fftgr=(8*cplex*nfft*n_fftgr*nspden*mffmem)/1024._dp**2 + 0.002_dp
 end if

!---------------------------------------------------------------------
!Now, analyze the data

!DEBUG
!write(std_out,*)' memorf : nchain=',nchain
!ENDDEBUG

 ABI_ALLOCATE(cfft_dum,(marrays))
 cfft_dum=zero
 mbgylm=zero
 call memana(cadd,cfft,cfft_dum,chain,cmpw,dttyp,iout,iprcel,iscf,&
& marrays,mbcg,mbdiskpd,mbdiskwf,mbf_fftgr,mbgylm,mffmem,&
& mpw,natom,nchain,nfft,nfft,occopt,option,prtvol)
 ABI_DEALLOCATE(cfft_dum)

end subroutine memorf
!!***
