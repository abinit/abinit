!!****m* ABINIT/m_memeval
!! NAME
!! m_memeval
!!
!! FUNCTION
!!  Functions to estimate memory requirements from the calculation parameters.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (XG, DC, DW)
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

MODULE m_memeval

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_dtset

 use defs_datatypes, only : pspheader_type
 use defs_abitypes,   only : MPI_type
 use m_geometry,      only : mkradim, mkrdim, xred2xcart, metric
 use m_symtk,         only : mati3inv, littlegroup_q
 use m_spgdata,       only : prtspgroup
 use m_fftcore,       only : getng
 use m_kg,            only : getmpw
 use m_libpaw_tools,  only : libpaw_write_comm_set

 implicit none

 private
!!***

 public :: memory_eval   ! Main entry point
 public :: getdim_nloc   ! Determine the dimensions of arrays with non-local projectors: ekb, ffspl, indlmn
 public :: setmqgrid     ! Sets the number of points needed to represent the pseudopotentials in q-space

contains

!!****f* m_memeval/memory_eval
!! NAME
!! memory_eval
!!
!! FUNCTION
!! Big loop on the datasets:
!! - for each of the datasets, write one line about the crystallographic data
!! - compute the memory needs for this data set.
!!
!! INPUTS
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables
!!  iout=unit number of output file
!!  mpi_enregs=information about MPI parallelization
!!  ndtset= number of datasets to be read; if 0, no multi-dataset mode
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!      one data set.
!!  npsp=number of pseudopotentials
!!  pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!
!! OUTPUT
!!   printing only
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      atomic_info,createwavefunctionsdescriptors,deallocate_lr
!!      memoryestimator,mkradim,wrtout,wvl_descr_atoms_set,wvl_descr_free
!!      wvl_setboxgeometry,xred2xcart
!!
!! SOURCE

subroutine memory_eval(dtsets,iout,mpi_enregs,ndtset,ndtset_alloc,npsp,pspheads)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,ndtset,ndtset_alloc,npsp
 type(MPI_type),intent(inout) :: mpi_enregs(0:ndtset_alloc)
!arrays
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
 type(pspheader_type),intent(in) :: pspheads(npsp)

!Local variables -------------------------------
!scalars
 integer :: cplex,exchn2n3d,extrapwf,getcell,idtset,ii,intxc,densfor_pred,iprcel
 integer :: iscf,isym,jdtset,lmnmax,mem_test
 integer :: lmnmax_eff,lmnmaxso,lnmax,lnmax_eff,lnmaxso,mband
 integer :: me_fft,mffmem,mgfftdiel,mgfftf,mkmem,mpsang,mpspso
 integer :: mpssoang,mpw,mqgrid,mqgriddg,mqgrid_ff,mqgrid_vl,n1xccc,natom
 integer :: nfftdiel,nfftf,nkpt,nproc_fft,nptsgvec,npulayit,npwdiel,nspden,nspinor
 integer :: nsppol,nsym,ntypat,occopt,optddk,optforces,optphon,optstress
 integer :: optstrs,paral_fft,pawcpxocc,pawmixdg,pawnhatxc,pawspnorb,pawstgylm,prtvol,ptgroupma,response
 integer :: spgroup,timrev,usepaw,useylm,use_gpu_cuda,xclevel
 real(dp) :: diecut,dilatmx,ecut,ecut_eff,ecutdg_eff,ecutsus,ucvol
!arrays
 integer :: bravais(11),mkmems(3),ngfftdiel(18)
 integer :: ngfftf(18),nloalg(3)
 integer,allocatable :: nband(:),symq(:,:,:),symrec(:,:,:),symrel(:,:,:)
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp) :: genafm(3),gmet(3,3),gprimd(3,3),kpt_diel(3),qphon(3),rmet(3,3),rprimd(3,3)

!*************************************************************************

 do idtset=1,ndtset_alloc
   if(mpi_enregs(idtset)%me<0) cycle
   call abi_io_redirect(new_io_comm=mpi_enregs(idtset)%comm_world)
   call libpaw_write_comm_set(mpi_enregs(idtset)%comm_world)

!  Initialisations
   bravais(:)=dtsets(idtset)%bravais(:)
   exchn2n3d=dtsets(idtset)%exchn2n3d
   extrapwf=dtsets(idtset)%extrapwf
   genafm(:) =dtsets(idtset)%genafm(:)
   getcell=dtsets(idtset)%getcell
   intxc=dtsets(idtset)%intxc
   densfor_pred=dtsets(idtset)%densfor_pred
   iprcel=dtsets(idtset)%iprcel
   iscf=dtsets(idtset)%iscf
   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0
   me_fft=mpi_enregs(idtset)%me_fft
   mffmem=dtsets(idtset)%mffmem
   mpw=dtsets(idtset)%mpw
   mqgrid=dtsets(idtset)%mqgrid
   mqgriddg=dtsets(idtset)%mqgriddg
   natom=dtsets(idtset)%natom
   nkpt  =dtsets(idtset)%nkpt
   nloalg(:)=dtsets(idtset)%nloalg(:)
   nproc_fft=mpi_enregs(idtset)%nproc_fft
   npulayit=dtsets(idtset)%npulayit
   nspden=dtsets(idtset)%nspden
   nspinor=dtsets(idtset)%nspinor
   nsppol=dtsets(idtset)%nsppol
   nsym     =dtsets(idtset)%nsym
   ntypat=dtsets(idtset)%ntypat
   occopt=dtsets(idtset)%occopt
   optforces=dtsets(idtset)%optforces
   paral_fft=mpi_enregs(idtset)%paral_kgb
   pawcpxocc=dtsets(idtset)%pawcpxocc
   pawmixdg=dtsets(idtset)%pawmixdg
   pawnhatxc=dtsets(idtset)%pawnhatxc
   pawspnorb=dtsets(idtset)%pawspnorb
   pawstgylm=dtsets(idtset)%pawstgylm
   prtvol=dtsets(idtset)%prtvol
   ptgroupma =dtsets(idtset)%ptgroupma
   qphon(:)=dtsets(idtset)%qptn(:)
   spgroup   =dtsets(idtset)%spgroup
   usepaw=dtsets(idtset)%usepaw
   useylm=dtsets(idtset)%useylm
   use_gpu_cuda=dtsets(idtset)%use_gpu_cuda
   xclevel=dtsets(idtset)%xclevel

   ABI_MALLOC(symrel,(3,3,nsym))
   symrel(:,:,1:nsym)=dtsets(idtset)%symrel(:,:,1:nsym)

!  Space group output
   call prtspgroup(bravais,genafm,iout,jdtset,ptgroupma,spgroup)

   if (dtsets(idtset)%toldff>tol16.and.optforces==0) optforces=1
   if (dtsets(idtset)%tolrff>tol16.and.optforces==0) optforces=1
   if (dtsets(idtset)%ionmov>tol16.and.optforces==0) optforces=1
   if (dtsets(idtset)%imgmov>tol16.and.optforces==0) optforces=1
   optstress=dtsets(idtset)%optstress
   optddk=0;optphon=0;optstrs=0
   if (dtsets(idtset)%rfddk>0.or.dtsets(idtset)%rf2_dkdk>0.or.dtsets(idtset)%rf2_dkde>0) optddk=1
   if (dtsets(idtset)%rfelfd>0.or.dtsets(idtset)%d3e_pert1_elfd>0.or.&
&   dtsets(idtset)%d3e_pert2_elfd>0.or.dtsets(idtset)%d3e_pert3_elfd>0) optddk=1
   if (dtsets(idtset)%rfphon>0.or.dtsets(idtset)%d3e_pert1_phon>0.or.&
&   dtsets(idtset)%d3e_pert2_phon>0.or.dtsets(idtset)%d3e_pert3_phon>0) optphon=1
   if (dtsets(idtset)%rfstrs>0) optstrs=1

   ABI_MALLOC(nband,(nkpt*nsppol))
   nband(1:nkpt*nsppol)=dtsets(idtset)%nband(1:nkpt*nsppol)
   mband=maxval(nband(1:nkpt*nsppol))
   dtsets(idtset)%mband=mband

!  mpsang=max(maxval(pspheads(1:npsp)%lmax)+1,1) ! Likely problems with the HP compiler
!  n1xccc=maxval(pspheads(1:npsp)%xccc)
   mpsang=1
   n1xccc=pspheads(1)%xccc
   do ii=1,npsp
     mpsang=max(pspheads(ii)%lmax+1,mpsang)
     n1xccc=max(pspheads(ii)%xccc,n1xccc)
   end do

!  Determine the maximum number of projectors, for the set of pseudo atom
   call getdim_nloc(lmnmax,lmnmaxso,lnmax,lnmaxso,dtsets(idtset)%mixalch_orig,&
&   dtsets(idtset)%nimage,npsp,dtsets(idtset)%npspalch,ntypat,dtsets(idtset)%ntypalch,pspheads)

!  Treatment of the effect of using a spin-orbit part
!  Warning : mpspso is different for each dataset; not relevant for PAW
   mpspso=1
   if (dtsets(idtset)%usepaw==0) then
     do ii=1,npsp
       if(nspinor/=1)then
         if(pspheads(ii)%pspso/=0)then
           if(dtsets(idtset)%so_psp(ii)/=0)then
             mpspso=2
           end if
         end if
       end if
     end do
   end if
!  In case of no spin-orbit
   if(mpspso==1)then
     mpssoang=mpsang ; lmnmax_eff =lmnmax; lnmax_eff =lnmax
   else ! spin-orbit will be used
     mpssoang=2*mpsang-1 ; lmnmax_eff =lmnmaxso ; lnmax_eff =lnmaxso
   end if
!  lmnmax is not used if the Ylm are not used
   if (useylm==0) lmnmax_eff =lnmax_eff

   ecut     =dtsets(idtset)%ecut
   dilatmx  =dtsets(idtset)%dilatmx
   ecut_eff=ecut*dilatmx**2
   ecutdg_eff=dtsets(idtset)%pawecutdg*dtsets(idtset)%dilatmx**2

!  Compute mgfft,mpw,nfft for this data set
   call mkrdim(dtsets(idtset)%acell_orig(1:3,1),dtsets(idtset)%rprim_orig(1:3,1:3,1),rprimd)
   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

   if (usepaw==0) then
     mgfftf=dtsets(idtset)%mgfft;nfftf=dtsets(idtset)%nfft;ngfftf(:)=dtsets(idtset)%ngfft(:)
   else
     mgfftf=dtsets(idtset)%mgfftdg;nfftf=dtsets(idtset)%nfftdg;ngfftf(:)=dtsets(idtset)%ngfftdg(:)
   end if
   response=0
   if(dtsets(idtset)%rfddk/=0  .or. dtsets(idtset)%rf2_dkdk/=0 .or. dtsets(idtset)%rf2_dkde/=0 .or. &
&   dtsets(idtset)%rfphon/=0 .or. dtsets(idtset)%rfelfd/=0 .or. &
&   dtsets(idtset)%rfstrs/=0 .or. dtsets(idtset)%rfuser/=0 .or. &
&   dtsets(idtset)%rfmagn/=0    ) response=1

!  Compute mgfftdiel,npwdiel,nfftdiel for this data set
   if((modulo(iprcel,100)>=20 .and.modulo(iprcel,100) < 71).or. iscf==-1)then
!    Get diecut, and the fft grid to be used for the susceptibility computation
     diecut=abs(dtsets(idtset)%diecut)
     if( dtsets(idtset)%diecut < zero )then
       ecutsus=ecut
     else
       ecutsus= ( sqrt(ecut) *0.5_dp + sqrt(diecut) *0.25_dp )**2
     end if
!    Beware, for the dielectric matrix fftalg=ngfftdiel(7) is default here
     ngfftdiel(1:3)=0 ; ngfftdiel(7)=101 ; ngfftdiel(8:18)=dtsets(idtset)%ngfft(8:18)
     if(iscf==-1)ngfftdiel(7)=102
     ecut_eff=ecutsus*dilatmx**2
     call getng(dtsets(idtset)%boxcutmin,ecut_eff,gmet,k0,me_fft,mgfftdiel,nfftdiel,&
&     ngfftdiel,nproc_fft,nsym,paral_fft,symrel,&
&     use_gpu_cuda=dtsets(idtset)%use_gpu_cuda)
!    Compute the size of the dielectric matrix : npwdiel
     kpt_diel(1:3)=(/ 0.0_dp, 0.0_dp, 0.0_dp /)
     ecut_eff=diecut*dilatmx**2
     call getmpw(ecut_eff,exchn2n3d,gmet,(/1/),kpt_diel,mpi_enregs(idtset),npwdiel,1)
   else
     npwdiel=1 ; mgfftdiel=1 ; nfftdiel=1 ; ngfftdiel(1:8)=1
   end if

!  Special treatment for the value of mqgrid to be fed in memory.F90

   nptsgvec         = 200 ! At present, this has to be chosen once and for all ...
   if ( dtsets(idtset)%usewvl == 0) then
     call setmqgrid(mqgrid,mqgriddg,ecut_eff,ecutdg_eff,gprimd,nptsgvec,usepaw)
   else
     call setmqgrid(mqgrid,mqgriddg,one,one,gprimd,nptsgvec,usepaw)
   end if
   mqgrid_ff=mqgrid
   if (usepaw==0) mqgrid_vl=mqgrid
   if (usepaw==1) mqgrid_vl=mqgriddg

!  Compute the memory needs for this data set.
   if(response==0)then

     if (dtsets(idtset)%usewvl == 0) then
       mkmem=dtsets(idtset)%mkmem
       mband=maxval(dtsets(idtset)%nband(1:nkpt*nsppol))

       ! Don't perform memory tests if MBPT.
       mem_test = dtsets(idtset)%mem_test
       if (any(dtsets(idtset)%optdriver == [RUNL_SIGMA, RUNL_SCREENING, RUNL_BSE, RUNL_EPH])) mem_test = 0

       call memory(n1xccc,extrapwf,getcell,idtset,dtsets(idtset)%icoulomb,&
&       intxc,dtsets(idtset)%ionmov,iout,densfor_pred,&
&       iprcel,iscf,jdtset,lmnmax_eff,lnmax_eff,mband,mffmem,dtsets(idtset)%mgfft,mgfftdiel,mgfftf,mkmem,&
&       mpi_enregs(idtset),mpsang,mpssoang,mpw,mqgrid_ff,mqgrid_vl,natom,nband,dtsets(idtset)%nfft,nfftdiel,nfftf,&
&       dtsets(idtset)%ngfft,ngfftdiel,ngfftf,dtsets(idtset)%nimage,nkpt,nloalg,npsp,npulayit,npwdiel,nspden,nspinor,&
&       nsppol,nsym,ntypat,occopt,optforces,mem_test,optstress,pawcpxocc,pawmixdg,&
&       pawnhatxc,pawspnorb,pawstgylm,prtvol,pspheads,qphon,dtsets(idtset)%tfkinfunc,&
&       dtsets(idtset)%typat,ucvol,usepaw,useylm,use_gpu_cuda,xclevel)
     else if( dtsets(idtset)%usepaw==0) then
       if (mpi_enregs(idtset)%me == 0) then
         call wvl_memory(dtsets(idtset), idtset, mpi_enregs(idtset), npsp, 1, pspheads)
       end if
     end if

   else
!    Compute the value of cplex, for which one needs symrec
     ABI_MALLOC(symq,(4,2,nsym))
     ABI_MALLOC(symrec,(3,3,nsym))
     do isym=1,nsym
       call mati3inv(symrel(:,:,isym),symrec(:,:,isym))
     end do
     call littlegroup_q(nsym,qphon,symq,symrec,dtsets(idtset)%symafm,timrev)
     cplex=2-timrev
     ABI_FREE(symq)
     ABI_FREE(symrec)
     mkmems(1)=dtsets(idtset)%mkmem
     mkmems(2)=dtsets(idtset)%mkqmem
     mkmems(3)=dtsets(idtset)%mk1mem

     mem_test = dtsets(idtset)%mem_test

     call memorf(cplex,n1xccc,getcell,idtset,intxc,iout,iprcel,&
&     iscf,jdtset,lmnmax_eff,lnmax_eff,mband,mffmem,dtsets(idtset)%mgfft,&
&     mkmems,mpi_enregs(idtset),mpsang,mpssoang,mpw,mqgrid_ff,natom,nband,dtsets(idtset)%nfft,&
&     dtsets(idtset)%ngfft,nkpt,nloalg,nspden,nspinor,nsppol,nsym,&
&     ntypat,occopt,optddk,optphon,mem_test,optstrs,prtvol,useylm,use_gpu_cuda,xclevel)
   end if

!  Deallocate temporary arrays (when they will really be temporary !)
   ABI_FREE(nband)
   ABI_FREE(symrel)

 end do ! idtset

end subroutine memory_eval
!!***

!!****f* m_memeval/memory
!! NAME
!! memory
!!
!! FUNCTION
!! Estimation of the memory needed for a ground-state job.
!! According to the value of the option variable,
!! might also try to allocate this amount of memory, and if it fails,
!! might estimate the available memory.
!!
!! INPUTS
!!  extrapwf=flag controlling the extrapolation of wave functions during MD or relaxation
!!  getcell=if non-zero, the values of acell and rprim are taken from
!!   the output of another dataset
!!  idtset=number of the current dataset
!!  icoulomb=0 for periodic Fourier calculation of Hartree potential; 1 for isolated system using Poisson solver.
!!  intxc=control xc quadrature
!!  ionmov=control force calculations
!!  iout=unit number for output of formatted data.
!!  densfor_pred=govern the choice of density prediction and/or forces correction
!!  iprcel=govern the choice of preconditioner for the SCF cycle
!!  iscf=governs the choice of SCF algorithm, or non-SCF calculation.
!!  jdtset=index of the current dataset
!!  lmnmax=max. number of (l,m,n) components over all type of psps
!!  lnmax =max. number of (l,n)   components over all type of psps
!!  mband =maximum number of bands
!!  mffmem =governs the number of FFT arrays which are fit in core memory
!!  mgfftf =maximum single fft dimension (fine grid, if PAW)
!!  mgfft  =maximum single fft dimension (coarse grid, if PAW)
!!  mgfftdiel =maximum single fft dimension for susceptibility and dielectric
!!   matrices.
!!  mkmem =maximum number of k points which can fit in core memory
!!  mpi_enreg=information about MPI parallelization
!!  mpssang is 1+maximum angular momentum for nonlocal pseudopotential
!!  mpssoang is 1+maximum (spin*angular momentum) for nonlocal pseudopotential
!!  mpw   =maximum number of planewaves in basis sphere (large number)
!!  mqgrid_ff=dimension of q (or G) grid for nl form factors (array ffspl)
!!  mqgrid_vl=dimension of q (or G) grid for Vloc (array vlspl)
!!  natom =number of atoms in unit cell
!!  nband(nkpt*nsppol)=number of bands at each k point, for each polarization
!!  nfftf =number of fft grid points for density        (fine grid, if PAW)
!!  nfft  =number of fft grid points for wavefunctions  (coarse grid, if PAW)
!!  nfftdiel  =maximum number of fft grid points for susceptibility
!!    and dielectric matrices
!!  ngfftf(18)=contain all needed information about 3D FFT (fine grid, if PAW)
!!  ngfft(18) =contain all needed information about 3D FFT (coarse grid, if PAW)
!!  ngfftdiel(18)=contain all needed information about 3D FFT, dielectric case,
!!                 see ~abinit/doc/variables/vargs.htm#ngfft
!!    for susceptibility and dielectric matrices
!!  nimage=number of images (replicas) of the cell
!!  nkpt  =number of k points
!!  npsp=number of different pseudopotentials
!!  npwdiel=number of plane wave for susceptibility and dielectric matrix
!!  npulayit=number of iterations used in Pulay SCF mixing
!!  nloalg(3)=governs the choice of the algorithm for non-local operator.
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=number of channels for spin-polarization (1 or 2)
!!  nsym  =number of symmetry elements in space group
!!  ntypat =number of types of atoms
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  occopt=option for occupation numbers. If 3<=occopt<=8, varying occupation
!!  optforces=1 if forces are computed during run
!!  option : if 0 , no test of available memory
!!           if 1 , the routine tries to allocate the estimated memory, for testing
!!                    purposes, and if a failure occurs, the routine stops.
!!           if 2 , like 1, but before stopping, the routine will provide
!!                    an estimation of the available memory.
!!  optstress=1 if stresses are computed during run
!!  pawcpxocc=2 if PAW occupancies (rhoij) are complex
!!  pawmixdg=1 if mixing (in PAW) is done on the fine grid
!!  pawnhatxc=1 if nhat PAW density has to be analytically included in XC
!!  pawspnorb=1 when spin-orbit is activated within PAW
!!  pawstgylm=1 if g_l(r).Y_lm(r) factors are stored in memory (PAW)
!!  prtvol=control print volume
!!  pspheads(npsp)=<type pspheader_type>all the important information from the header
!!  tfkinfun=flag controling the use of Thomas-Fermi algorithme (without WF)
!!  typat(natom)=type of each atom
!!  ucvol= unit cell volume
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  useylm=governs the way the nonlocal operator is to be applied:
!!         1=using Ylm, 0=using Legendre polynomials
!!  use_gpu_cuda=1 if Cuda (GPU) is on
!!  xclevel=XC functional level
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! for the estimation, it is only taken into account those
!! arrays that have some probability of being larger than 1000*8 bytes :
!! - All the arrays that have large numbers as one of their dimensions
!! (mqgrid, mpw, nfft, ngfft(4)*ngfft(5)*ngfft(6),
!!                     ngfftdiel(4)*ngfftdiel(5)*ngfftdiel(6), n1xccc
!!                                      or a constant larger than 1000)
!! - All the arrays that have a product of two moderately large numbers
!! (potential size above 30  : mband, mgfft, mkmem, natom, nkpt, nsym,
!!  or a constant larger than 30)
!! After this estimation, an amount of (176 + 55 + 6*natom) Kbytes is added
!! to take into account the static arrays declared
!! in rhotoxc and daughter routines (at maximum 22*1000 dp numbers),
!! as well as other arrays like
!! character(len=500) :: message (present in about 100 routines), or the different
!! arrays allocated in move.f, brdmin.f, gstate.f (xf array) or pspini.f
!! In the case 3<=occopt<=8 this amount is increased by 760 Kbytes
!! to take into account the arrays smdfun, occfun, entfun, workfun and xgrid,
!! declared in getnel.
!!
!! The current version takes into account
!! 1) and 2) the "main chain" in its two slightly different versions :
!! driver - gstate - (move or brdmin) - scfcv - vtorho - vtowfk -
! !!     cgwf - getghc - fourwf or (nonlop+opernl)
!! 3) the xc chain :
!! driver - gstate - (move or brdmin) - scfcv - (vresfo) - rhotoxc - xcden
!! 4) the mkrho chain :
!! driver - gstate - (move or brdmin) - scfcv - vtorho - mkrho
!! 5) the fourdp chain :
!! driver - gstate - (move or brdmin) - scfcv - vtorho
!!         ( + ftofr - fourdp - symrhg )
!! 6) the parallel k-point chain :
!! driver - gstate - (move or brdmin) - scfcv - vtorho - MPI_ALLREDUCE
!! 7) the newvtr chain :
!! driver - gstate - (move or brdmin) - scfcv - newvtr
!! 8) the susceptibility chain :
!! driver - gstate - (move or brdmin) - scfcv - vtorho - suscep - suskmm
!! 9) the dielectric chain :
!! driver - gstate - (move or brdmin) - scfcv - vtorho - dielmt
!! 10) the tddft chain :
!! driver - gstate - (move or brdmin) - scfcv - vtorho - tddft
!!
!! It is valid for all values of iscf, but not for nstep=0 (when the chain
!! goes through energy instead of vtorho).
!!
!! Also, it is assumed that the potentials are non-local, even if there
!! are local ! It would be necessary to update this routine
!! now that the beginning of psp files is read before the present call (XG 980502)
!!
!! One might also estimate if there must be a chain arriving at:
!!  strnps , mkffnl, mkcore, mklocl, mkrho, prcpot, irrzg, initro, clnup1.
!! This is because there are allocated arrays in these routines.
!!
!! PARENTS
!!      m_memeval
!!
!! CHILDREN
!!      atomic_info,createwavefunctionsdescriptors,deallocate_lr
!!      memoryestimator,mkradim,wrtout,wvl_descr_atoms_set,wvl_descr_free
!!      wvl_setboxgeometry,xred2xcart
!!
!! SOURCE

subroutine memory(n1xccc,extrapwf,getcell,idtset,icoulomb,intxc,ionmov,iout,densfor_pred,iprcel,&
& iscf,jdtset,lmnmax,lnmax,&
& mband,mffmem,mgfft,mgfftdiel,mgfftf,mkmem,mpi_enreg,mpsang,mpssoang,mpw,mqgrid_ff,mqgrid_vl,&
& natom,nband,nfft,nfftdiel,nfftf,ngfft,ngfftdiel,ngfftf,nimage,&
& nkpt,nloalg,npsp,npulayit,npwdiel,nspden,nspinor,nsppol,nsym,ntypat,&
& occopt,optforces,option,optstress,pawcpxocc,pawmixdg,pawnhatxc,pawspnorb,pawstgylm,&
& prtvol,pspheads,qphon,tfkinfunc,typat,ucvol,usepaw,useylm,use_gpu_cuda,xclevel)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: extrapwf,getcell,icoulomb,idtset,intxc,ionmov,iout,densfor_pred
 integer,intent(in) :: iprcel,iscf,jdtset,lmnmax,lnmax,mband,mffmem,mgfft
 integer,intent(in) :: mgfftdiel,mgfftf,mkmem,mpsang,mpssoang,mpw,mqgrid_ff
 integer,intent(in) :: mqgrid_vl,n1xccc,natom,nfft,nfftdiel,nfftf,nimage,nkpt,npsp
 integer,intent(in) :: npulayit,npwdiel,nspden,nspinor,nsppol,nsym,ntypat
 integer,intent(in) :: occopt,optforces,option,optstress
 integer,intent(in) :: pawcpxocc,pawmixdg,pawnhatxc,pawspnorb,pawstgylm
 integer,intent(in) :: prtvol,tfkinfunc,usepaw,useylm,use_gpu_cuda,xclevel
 real(dp) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: nband(nkpt*nsppol),ngfft(18),ngfftdiel(18),ngfftf(18)
 integer,intent(in) :: nloalg(3),typat(natom)
 real(dp),intent(in) :: qphon(3)
 type(pspheader_type) :: pspheads(npsp)

!Local variables-------------------------------
!marrays=maximal number of arrays to be monitored (or group of arrays)
!cmpw(marrays)=count of blocks of size mpw bytes
!cfft(marrays) =count of blocks of size nfft bytes (coarse grid, if PAW)
!cfftf(marrays)=count of blocks of size nfft bytes (fine grid, if PAW)
!cadd(marrays)=count of additional storage needed (in bytes)
!dttyp(marrays)=datatype of the array : 4 for integers, 8 for real(dp)
!nchain=number of different chains of routines
!chain(marrays,nchain)=different chains of routines
 ! The cfoo arrays are used to store the allocated memory in the different
 ! routines of the program. Each stack of the program can allocate some
 ! memory and the amount is estimated and stored in cfoo(i). The lower i,
 ! the higher routine. cfft is memory used by FFT handling, cmpw for
 ! plane waves storage and cadd for miscellaneous memory occupation.
 ! The unit is the multiplier of the size of nfft for cfft, the multiplier
 ! of mpw for cmpw and the actually allocated memory for cadd.
 ! This array stores the size of each chunk of memory (8 for double
 ! floating point precision, 4 for integers and so on).
 ! This array defines if the chain defined above allocate or not the
 ! memory (depending on options).
!scalars
 integer,parameter :: marrays=150,nchain=10
 integer :: fftalgb,histsz,ii,iscf10,jj,l_max,l_size_max,matblk,mblk,mincat,mu
 integer :: my_natom,n_fftgr,narr_fourdp,nbnd_in_blk,ndiel4,ndiel456,ndiel5,ndiel6
 integer :: ngrad,nprocwf,nspgrad,qphase_rhoij,rhoij_nspden
 real(dp) :: mbcg,mbdiskpd,mbdiskwf,mbf_fftgr,mbgylm
 character(len=500) :: msg
! character(len=1) :: firstchar
!arrays
 integer :: dttyp(marrays),nattyp(ntypat)
 integer,allocatable :: basis_size(:),l_size(:),lmn2_size(:),lmn_size(:)
 integer,allocatable :: mesh_size(:),my_nattyp(:),pawver(:),shape_type(:)
 real(dp) :: cadd(marrays),cfft(marrays),cfftf(marrays),cmpw(marrays)
 real(dp),allocatable :: rshp(:)
 logical :: chain(marrays,nchain)

! **************************************************************************

 if(option<0 .or. option>2)then
   write(msg,'(A,I0,A)')'option=',option,' while the only allowed values are 0, 1, or 2.'
   ABI_BUG(msg)
 end if

!firstchar=' ';if (use_gpu_cuda==1) firstchar='_'
 cmpw(:)=zero ; cfft(:)=zero ; cfftf(:)=zero ; cadd(:)=zero
 dttyp(:)=0

 my_natom=natom;if (mpi_enreg%nproc_atom>1) my_natom=mpi_enreg%my_natom

 call wrtout(std_out,'memory: analysis of memory needs ','COLL')

 if(jdtset>=100)then
   write(msg,'(80a,a,a,i5,a)')('=',mu=1,80),ch10,&
    ' Values of the parameters that define the memory need for DATASET',jdtset,'.'
 else if(jdtset/=0)then
   write(msg,'(80a,a,a,i3,a)')('=',mu=1,80),ch10,&
    ' Values of the parameters that define the memory need for DATASET',jdtset,'.'
 else
   write(msg,'(80a,a,a)')('=',mu=1,80),ch10,&
    ' Values of the parameters that define the memory need of the present run '
 end if
 call wrtout(iout,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 write(msg,'( 4(a,i8),a,4(a,i8) )' ) &
& '     intxc =',intxc   ,'    ionmov =',ionmov,&
& '      iscf =',iscf    ,'    lmnmax =',lmnmax,ch10,&
& '     lnmax =',lnmax   ,'     mgfft =',mgfft,&
& '  mpssoang =',mpssoang,'    mqgrid =',mqgrid_vl
 call wrtout(iout,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 write(msg,'( 4(a,i8),a,4(a,i8),a,4(a,i8) )' ) &
& '     natom =',natom  ,'  nloc_mem =',nloalg(2)*(nloalg(3)+1),&
& '    nspden =',nspden ,'   nspinor =',nspinor,ch10,&
& '    nsppol =',nsppol ,'      nsym =',nsym,&
& '    n1xccc =',n1xccc ,'    ntypat =',ntypat,ch10,&
& '    occopt =',occopt ,'   xclevel =',xclevel
 call wrtout(iout,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 write(msg,'(4(3(a,i12),a))') &
& '-    mband =',mband  ,'        mffmem =',mffmem,&
& '         mkmem =',mkmem  ,ch10,&
& '       mpw =',mpw    ,'          nfft =',nfft ,&
& '          nkpt =',nkpt
 call wrtout(iout,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 if (my_natom/=natom)then
   write(msg,'(a,i10)') 'Pmy_natom=',my_natom
   call wrtout(iout,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
 end if

!Additional information if imgmov is activated (use of replicas of the cell)
 if (nimage>1) then
   write(msg,'(1(a,i10))' ) '  nimage =',nimage
   call wrtout(iout,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
 end if

!Additional information on FFT grids if PAW
 if (usepaw==1) then
   write(msg, '(a,a,a,i10,a,i10)' )&
&   ' PAW method is used; the additional fine FFT grid is defined by:',ch10,&
&   '   mgfftf=',mgfftf,'    nfftf =',nfftf
   call wrtout(iout,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
 end if

!Additional information if GPU
 if (use_gpu_cuda==1) then
!  write(msg, '(a)' )' GPU method is used'
!  call wrtout(iout,msg,'COLL')
!  call wrtout(std_out,msg,'COLL')
 end if

!Additional information needed for the susceptibility and dielectric matrices
 if((modulo(iprcel,100)>=20.and.modulo(iprcel,100)<70) .or. iscf==-1)then

!  Compute the number of bands in blocks (nbnd_in_blk) from mband (see suskmm.f)
!  Consider that if the number of bands is large, there are at most 8 blocks
   if(mband>=48)then
     mblk=8
     nbnd_in_blk=(mband-1)/mblk+1
!    If the number of bands is medium, place 6 bands per block
   else if(mband>=12)then
     nbnd_in_blk=6
!    Otherwise, must have at least 2 blocks
   else
     mblk=2
     nbnd_in_blk=(mband-1)/mblk+1
   end if

   write(msg, '(a,a,a,i10,a,i6,a,i10,a,i10)' )&
&   ' For the susceptibility and dielectric matrices, or tddft :',ch10,&
&   '   mgfft =',mgfftdiel,'  nbnd_in_blk=',nbnd_in_blk,'    nfft =',nfftdiel,&
&   '     npw =',npwdiel
   call wrtout(iout,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
   ndiel4=ngfftdiel(4) ; ndiel5=ngfftdiel(5) ; ndiel6=ngfftdiel(6)
   ndiel456=ndiel4*ndiel5*ndiel6
 else
!  To be sure of initialisation.
   ndiel456 = 1
 end if

 write(msg,'(80a)') ('=',mu=1,80)
 call wrtout(iout,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 if(getcell>0 .or. (getcell<0 .and. idtset+getcell>0) )then
   write(msg,'(a,a,a,a,a,a,i3,a,i3,a,a,a,a,a,a)' )ch10,&
&   ' memory : COMMENT -',ch10,&
&   '  The determination of memory needs at this stage is meaningless,',ch10,&
&   '  since getcell = ',getcell,' is non-zero, while idtset=',idtset,'.',ch10,&
&   '  The following numbers are obtained by supposing that acell and rprim',ch10,&
&   '  are NOT taken from a previous dataset. You cannot rely on them.',ch10
   call wrtout(iout,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
 end if

!Compute number of atoms per type for current proc
 nattyp(:)=0
 do ii=1,natom
   nattyp(typat(ii))=nattyp(typat(ii))+1
 end do

!PAW: store useful dims
 if (usepaw==1) then
   ABI_MALLOC(basis_size,(npsp))
   ABI_MALLOC(l_size,(npsp))
   ABI_MALLOC(lmn_size,(npsp))
   ABI_MALLOC(lmn2_size,(npsp))
   ABI_MALLOC(mesh_size,(npsp))
   ABI_MALLOC(shape_type,(npsp))
   ABI_MALLOC(pawver,(npsp))
   ABI_MALLOC(rshp,(npsp))
   do ii=1,npsp
     basis_size(ii)=pspheads(ii)%pawheader%basis_size
     mesh_size(ii)=pspheads(ii)%pawheader%mesh_size
     l_size(ii)=pspheads(ii)%pawheader%l_size
     lmn_size(ii)=pspheads(ii)%pawheader%lmn_size
     lmn2_size(ii)=lmn_size(ii)*(lmn_size(ii)+1)/2
     pawver(ii)=pspheads(ii)%pawheader%pawver
     rshp(ii)=pspheads(ii)%pawheader%rshp
     shape_type(ii)=pspheads(ii)%pawheader%shape_type
   end do
   l_max=maxval(pspheads(:)%lmax)
   l_size_max=maxval(pspheads(:)%pawheader%l_size)
   rhoij_nspden=nspden;if (pawspnorb>0) rhoij_nspden=4
   ABI_MALLOC(my_nattyp,(ntypat))
   if ((mpi_enreg%nproc_atom<=1).or.(.not.associated(mpi_enreg%my_atmtab))) then
     my_nattyp=nattyp
   else
     my_nattyp=0
     do ii=1,my_natom
       jj=typat(mpi_enreg%my_atmtab(ii))
       my_nattyp(jj)=my_nattyp(jj)+1
     end do
   end if
   qphase_rhoij=merge(2,1,any(qphon(:)>tol8))
 else
!  Do the allocation to avoid uninitialised variables.
   ABI_MALLOC(my_nattyp,(1))
   ABI_MALLOC(basis_size,(1))
   ABI_MALLOC(l_size,(1))
   ABI_MALLOC(lmn_size,(1))
   ABI_MALLOC(lmn2_size,(1))
   ABI_MALLOC(mesh_size,(1))
   ABI_MALLOC(shape_type,(1))
   ABI_MALLOC(pawver,(1))
   ABI_MALLOC(rshp,(1))
   rhoij_nspden=nspden
   l_size_max=1
   l_max=1
   qphase_rhoij=1
 end if

 n_fftgr=1;iscf10=mod(iscf,10)
 if(iscf10==1)              n_fftgr=5
 if(iscf10==2)              n_fftgr=3
 if(iscf10==3)              n_fftgr=4
 if(iscf10==4)              n_fftgr=6
 if(iscf10==5.or.iscf10==6) n_fftgr=10
 if(iscf10==7)              n_fftgr=2+2*npulayit

!work1 and work2 in fourdp : take into account approximately fftalgb
 fftalgb=mod(ngfft(7),100)/10
 if(fftalgb==0)narr_fourdp=2*2
 if(fftalgb==1)narr_fourdp=2

 ngrad=1;if(xclevel==2.or.tfkinfunc>10)ngrad=2

!(1)                     in main, driver, gstate and brdmin ----------------
!in move, nothing interesting is allocated.
!kg (gstate)
 cmpw(1)=3*mkmem               ; dttyp(1)=4
!indsym (gstate)
 cadd(3)=4*nsym*natom          ; dttyp(3)=4
!irrzon  (gstate)
 if(nsym/=1)then
   cfft(4)=2*((nspden/nsppol)-3*(nspden/4))    ; dttyp(4)=4
 end if
!ylm (gstate)
 cmpw(5)=mkmem*mpsang*mpsang*useylm ; dttyp(5)=8
!
!rhor,rhog (gstate)
 cfftf(5)=nspden+2              ; dttyp(5)=8
!cg (gstate)
 cmpw(6)=2*nspinor*mband*mkmem*nsppol  ; dttyp(6)=8
!eigen,resid,occ (occ is initialized in abinit, and not in driver)
 cadd(7)=3*mband*nkpt*nsppol   ; dttyp(7)=8
!qgrid_vl,qgrid_ff,vlspl,ffspl
 cadd(8)=mqgrid_vl*(1+2*ntypat)   &
& +mqgrid_ff*(1+2*ntypat*lnmax)   &
& ; dttyp(8)=8
!ph1d (actually allocated in scfcv !!)
 cadd(9)=2*3*(2*mgfft+1)*natom ; dttyp(9)=8
 cadd(9)=cadd(9)+2*3*(2*mgfftf+1)*natom*usepaw  !Additional ph1df for PAW
!phnons (in gstate)
 if(nsym/=1)then
   cfft(10)=2*((nspden/nsppol)-3*(nspden/4))    ; dttyp(10)=8
 end if
!xccc1d (in driver)
 cadd(11)=n1xccc*6*ntypat      ; dttyp(11)=8

!hessin in brdmin
 if(ionmov==2)then
   cadd(15)=3*natom*3*natom      ; dttyp(15)=8
 end if

!Additional PAW arrays
!PAW datasets (pawtab)
 if (usepaw==1) then
   dttyp(16)=8 ; dttyp(17)=4
   do ii=1,npsp
     cadd(16)=cadd(16)+2*mesh_size(ii)*basis_size(ii)   !phi,tphi
     cadd(16)=cadd(16)+2*mesh_size(ii)*basis_size(ii)&  !phiphj,tphiphj
&    *(basis_size(ii)+1)/2
     cadd(16)=cadd(16)+mesh_size(ii)*l_size(ii)         !shapefunc
     cadd(16)=cadd(16)+lmn2_size(ii)*l_size(ii)**2      !qijl
     cadd(16)=cadd(16)+l_size(ii)*5                     !gnorm,shape_a,shape_q
     cadd(16)=cadd(16)+lmn2_size(ii)*(4+lmn2_size(ii))  !eijkl,dltij,dij0,rhoij0,sij
     cadd(17)=cadd(17)+lmn2_size(ii)*8                  !indklmn
     cadd(16)=cadd(16)+mesh_size(ii)*5                  !coreden,tcoreden,rad,radfact,simfact
     if (shape_type(ii)==-1) cadd(16)=cadd(16)+4*mesh_size(ii)*l_size(ii)  !dshpfunc
     cadd(16)=cadd(16)+mqgrid_vl*2                      !tncorespl
     if (pawver(ii)>=4) cadd(16)=cadd(16)+mqgrid_vl*2   !tnvalespl
   end do
!  additional arrays
   cadd(16)=cadd(16)+l_size_max*2*l_max*nsym                 !zarot
   cadd(16)=cadd(16)+(2*l_max-1)**2*l_max**2*(l_max**2+1)/2  !realgnt
   cadd(17)=cadd(17)+nfft+nfftf                              ! fintocoa,coatofin
   do ii=1,ntypat
     cadd(16)=cadd(16)+my_nattyp(ii)*lmn2_size(ii)*rhoij_nspden*pawcpxocc ! Rhoij and related data
     cadd(17)=cadd(17)+my_nattyp(ii)*(2+lmn2_size(ii))    ! (rhoijselect, ...)
   end do
 end if

!SCF history (if selected)
 if (abs(densfor_pred)==5.or.abs(densfor_pred)==6) then          ! scf_history...
   histsz=2
   cfftf(18)=nspden*(histsz+1)+1      ; dttyp(18)=8  ! %deltarhor, %atmrho_last, %rhor_last
   cadd(19)=3*natom*2*histsz          ; dttyp(19)=8  ! %xreddiff,xred_last
   dttyp(20)=4
   if (usepaw==1) then
     do ii=1,ntypat
       cadd(19)=cadd(19)+histsz*2*my_nattyp(ii)*lmn2_size(ii)*rhoij_nspden*qphase_rhoij*pawcpxocc ! %pawrhoij()%rhoijp
       cadd(20)=cadd(20)+histsz*2*my_nattyp(ii)*(2+lmn2_size(ii))*nspden ! %pawrhoij()%rhoijselect
     end do
   end if
   if (extrapwf>0) then
     cadd(19)=cadd(19)+histsz*2*nspinor*mband*mkmem*nsppol  ; dttyp(19)=8  ! %cg
   end if
 end if

!(2)                     in scfcv----------------------------------------

!vhartr,vpsp,vtrial,vxc
 cfftf(21)=2+2*nspden           ; dttyp(21)=8
!kxc
 if (abs(densfor_pred)>0.and.iscf>=10) then
   cfftf(21)=cfftf(21)+3*nspden
   if (densfor_pred<0.and.xclevel==2) cfftf(21)=cfftf(21)+20*nspden
 end if
 if(iscf>0)then
!  f_fftgr
   if (pawmixdg==1) then
     cfftf(22)=nspden*n_fftgr*mffmem; dttyp(22)=8
   else
     cfft(22)=nspden*n_fftgr*mffmem; dttyp(22)=8
   end if
 end if
 if( iscf>0 .and. (modulo(iprcel,100)>=20.and.modulo(iprcel,100)<70))then
!  dielinv, susmat
   cadd(23)=4*(npwdiel*min(nspden,2))**2; dttyp(23)=8
 end if
!Kernel of Poisson's solver
 if (icoulomb == 1) then
   cadd(24) = ngfft(4)*ngfft(5)*ngfft(6) ; dttyp(24) = 8
 end if
 if( (iscf>0 .and. modulo(iprcel,100)>=20 .and. modulo(iprcel,100)<70) .or. iscf==-1 )then
!  kg_diel
   cadd(27)=3*npwdiel             ; dttyp(27)=4
   if(nsym/=1)then
!    irrzondiel
     cadd(27)=cadd(27)+2*nfftdiel*(nspden/nsppol)
!    phnonsdiel
     cadd(28)=2*nfftdiel*(nspden/nsppol)   ; dttyp(28)=8
   end if
 end if
 if(n1xccc/=0)then
!  xccc3d
   cfftf(29)=1                    ; dttyp(29)=8
 end if

!Additional PAW arrays
 dttyp(25)=8 ; dttyp(26)=4
 if (usepaw==1) then
   do ii=1,ntypat
     jj=(1+int(nfftf*four_pi*rshp(ii)**3/(three*ucvol)))        ! pawfgrtab
     cadd(26)=cadd(26)+my_nattyp(ii)*jj                         !   %ifftsph
     cadd(25)=cadd(25)+my_nattyp(ii)*jj*(1-pawstgylm)*3         !   %rfgd (if pawstgylm=0)
     cadd(25)=cadd(25)+my_nattyp(ii)*jj*pawstgylm*l_size(ii)**2 !   %gylm (if pawstgylm=1)
     if (optforces==1) cadd(25)=cadd(25)+my_nattyp(ii)*jj&      !   %gylmgr,%rfgd (if pawstgylm=1)
&    *pawstgylm*(3*l_size(ii)**2+3*optstress)
     cadd(26)=cadd(26)+my_nattyp(ii)*l_size(ii)**2/32           ! lmselect  !now a boolean
     cadd(25)=cadd(25)+my_nattyp(ii)*lmn2_size(ii)*nspinor**3   ! dij
     if (iscf>0) then
       cadd(25)=cadd(25)+my_nattyp(ii)*lmn2_size(ii)*rhoij_nspden*pawcpxocc                ! rhoijres
       cadd(25)=cadd(25)+my_nattyp(ii)*lmn2_size(ii)*rhoij_nspden*pawcpxocc*n_fftgr*mffmem ! f_paw
     end if
   end do
   cadd(25)=cadd(25)+(1+3*pawnhatxc*(ngrad/2))*nspden*nfftf       !nhat,nhatgr
 end if

!(3)                     in rhotoxc, xcden -------------------------------

 if(xclevel/=0)then
   if(n1xccc/=0)then
!    rhocorval
     cfftf(31)=nspden               ; dttyp(31)=8
   end if
!  dnexcdn, rhonow
   nspgrad=nspden*ngrad
   if(nspden==2 .and. ngrad==2)nspgrad=5
   cfftf(32)=nspden*ngrad*ngrad+nspgrad  ; dttyp(32)=8
   if(intxc==1 .or. ngrad==2)then
!    wkcmpx,work in xcden +work1,work2 in fourdp
     cfftf(33)=3+narr_fourdp        ; dttyp(33)=8
     cadd(33)=narr_fourdp*(ngfftf(4)*ngfftf(5)*ngfftf(6)-nfftf)
   end if
   if(ngrad==2)then
!    workgr in xcden
     cfftf(34)=2                    ; dttyp(34)=8
   end if
 end if
 if(iscf>0)then
!  In this case, rhotoxc is called from rhotov also,
!  for which vresid was allocated in scfcv
!  vresid
   cfftf(35)=nspden               ; dttyp(35)=8
 end if
!Poisson's solver with zero padding
 if (icoulomb == 1) then
   cfft(36) = 8                   ; dttyp(36) = 8
   cadd(36) = ngfft(4) * ngfft(5) * ngfft(6) - nfft
 end if

!Note : in hartre, called by rhotoxc, one uses
!2 dp arrays of total size 3*nfft,
!and 2 arrays of total size 4*n4*n5*n6 for fourdp
!This will be smaller than the total use for symrhg

!(4)                     in newvtr/newrho --------------------------------------

 if(iscf>0)then
!  vresid (allocated in scfcv) and vrespc
   if (pawmixdg==1) then
     cfftf(41)=2*nspden             ; dttyp(41)=8
   else
     cfft(41)=2*nspden             ; dttyp(41)=8
   end if
   if(mffmem==0)then
!    f_fftgr_disk
     if (pawmixdg==1) then
       cfftf(42)=nspden*n_fftgr       ; dttyp(42)=8
     else
       cfft(42)=nspden*n_fftgr       ; dttyp(42)=8
     end if
!    f_paw_disk
     if (usepaw==1) then
       dttyp(43)=8
       do ii=1,ntypat
         cadd(43)=cadd(43)+my_nattyp(ii)*lmn2_size(ii)*nspden*n_fftgr
       end do
     end if
   end if
!  rhoupdn, n(v)resid0, vtrialg, rhog2, magng
   if (pawmixdg==1) then
     cfftf(43)=2*nspden       ; dttyp(43)=8
   else
     cfft(43)=2*nspden       ; dttyp(43)=8
     if (nspden>1) cfftf(43)=2*(nspden-1)
   end if
 end if

!(5-6)                   in vtorho-----------------------------------------

!Note : (5) is for the arrays inside the spin and k-point loop
!they belong to the main chain
!(6) is for the arrays after the spin and k-point loop
!(6a) is for the arrays after that loop, for the parallel k-point chain
!(6b) is for the arrays in mkrho, for the mkrho chain
!(6c) is for the arrays in symrhg, for the fourdp chain
!(6d) is for the arrays in suscep, for the suscep chain, see (10)
!(6e) is for the arrays in dielmt, for the dielmt chain, see (11)
!(6f) is for the arrays in pawmkrhoij

!eknlk, enlxnk, grnlnk
 cadd(51)=(11+3*natom)*mband*nkpt*nsppol &
& ; dttyp(51)=8
!kg_k
 cmpw(52)=3                    ; dttyp(52)=4
!rhoaug,vlocal
 cfft(53)=2                    ; dttyp(53)=8
 cadd(53)=2*(ngfft(4)*ngfft(5)*ngfft(6)-nfft)
!rhowfr,rhowfg
 cfft(53)=cfft(53)+2+nspden
 if(mkmem==0)then
!  cg_disk
   cmpw(54)=2*nspinor*mband      ; dttyp(54)=8
 end if
!eig_k, ek_k, enlx_k, grnl_k, occ_k, resid_k
 cadd(56)=(14+3*natom)*mband   ; dttyp(56)=8
!ylm_k
 cmpw(57)=mpsang*mpsang*useylm ; dttyp(57)=8
!PAW:cprj
 if (usepaw==1) then
   dttyp(58)=8
   do ii=1,ntypat
     cadd(58)=cadd(58)+2*nattyp(ii)*nkpt*nspinor*mband*nsppol*lmn_size(ii)/max(mpi_enreg%nproc_band,1)
   end do
 end if

!(6)                     in vtorho----------------------------------------

!doccde
 cadd(60)=mband*nkpt*nsppol    ; dttyp(60)=8

!(6a)                    in vtorho----------------------------------------
 if(xmpi_paral==1)then
!  Parallel case
!  buffer1
!  buffer2
   if(occopt>=3 .and. occopt <=8) then
     dttyp(61)=8
     if(nsppol*nfft >= (13+3*natom)*mband*nkpt*nspden)then
       cfft(61)=2*nspden
     else
       cadd(61)=(13+3*natom)*mband*nkpt*nspden
     end if
   else
     cfft(61)=2*nspden             ; dttyp(61)=8
     cadd(61)=9+3*natom+2+2*mband*nkpt*nspden
   end if
 end if


!(6b)                    in mkrho, called by vtorho--------------------------
 if(occopt>=3 .and. occopt <=8)then
   if(mkmem==0)then
!    cg_disk
     cmpw(62)=2*nspinor*mband      ; dttyp(62)=8
   end if
!  cwavef
   cmpw(65)=2*nspinor            ; dttyp(65)=8

!  rhoaug, wfraug, work1 in fourwf
   cfft(66)=5                    ; dttyp(66)=8
   cadd(66)=5*(ngfft(4)*ngfft(5)*ngfft(6)-nfft)
 end if

!(6c)                    in symrhg, called by vtorho--------------------------
 if(iscf>0)then
   cfft(67)=narr_fourdp          ; dttyp(67)=8
   cadd(67)=narr_fourdp*(ngfft(4)*ngfft(5)*ngfft(6)-nfft)
   if(nsym>1)then
!    work1  in symrhg
     cfft(68)=2                    ; dttyp(68)=8
     cadd(68)=2*(ngfft(4)*ngfft(5)*ngfft(6)-nfft)
   end if
 end if


!(6d) and (6e)           in suscep and dielmt, called by vtorho,
!see (10) and (11) -------------------------------

!(6f)  in pawmkrhoij or pawrhoij_symrhoij called by pawmkrho, called by vtorho--------
!only when paralellim over atoms is activated
 dttyp(63)=8
 if((usepaw==1) .and. ((iscf>0) .or. (iscf == -3) .and. mpi_enreg%nproc_atom>1 ))then
   do ii=1,ntypat
     cadd(63)=cadd(63)+nattyp(ii)*lmn2_size(ii)*rhoij_nspden*pawcpxocc*qphase_rhoij ! Rhoij_gather and related data
     cadd(63)=cadd(63)+nattyp(ii)*(2+lmn2_size(ii)) ! Rhoij_gather (rhoijselect, ...)
   end do
 end if

!(7)                     in vtowfk----------------------------------------

!evec
 cadd(71)=2*mband*mband        ; dttyp(71)=8
!subham, subvnlx(if not PAW or if usefock_ACE)
 cadd(72)=(1+usepaw)*mband*(mband+1)    ; dttyp(72)=8
!gkpsq
 cmpw(73)=1                    ; dttyp(73)=8
!ffnl
 cmpw(74)=2*ntypat*lmnmax      ; dttyp(74)=8
!ph3d
 matblk=min(NLO_MINCAT,maxval(nattyp))
 if(nloalg(2)<=0)matblk=natom
 cmpw(75)=2*matblk             ; dttyp(75)=8
!gsc(if PAW)
 cmpw(76)=2*mband*nspinor*usepaw          ; dttyp(76)=8
!Note : matvnl and mat1 do not belong to a chain defined until now
!
 if(occopt<3 .and. iscf>0)then
!  cwavef
   cmpw(77)=2*nspinor            ; dttyp(77)=8
!  wfraug
   cfft(78)=2                    ; dttyp(78)=8
   cadd(78)=2*(ngfft(4)*ngfft(5)*ngfft(6)-nfft)
!  work1 in fourwf
   cfft(79)=2                    ; dttyp(79)=8
   cadd(79)=2*(ngfft(4)*ngfft(5)*ngfft(6)-nfft)
 end if

!(8)                     in cgwf------------------------------------------

!conjgr, cwavef, direc, gh_direc, gvnlx_direc
 cmpw(81)=2*5*nspinor          ; dttyp(81)=8
!ghc,gvnlxc
 cmpw(82)=2*2*nspinor          ; dttyp(82)=8
!PAW: scwavef,direc_tmp,ghc_all
 cmpw(83)=2*(2+mband)*nspinor*usepaw  ; dttyp(83)=8


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
!  dgxdt  (in nonlop)            !MT20072002: not allocated in getghc !!
   if (optforces==1) then
     cadd(95)=2*3*20*mincat*2    ; dttyp(95)=8
   end if
!  teffv (in opernl4 - no distinction is made for opernl, opernl2 or opernl3)
!  kpgx, ffkg
!  here, evaluate an upper value, with nproj=2, p,d and f orbitals, but not
!  considering the stress, since it will be called outside of the main chain
   cadd(97)=NLO_MBLKPW*40        ; dttyp(97)=8
!  kpg if nloalg(3)=1
   cadd(98)=3*mpw*nloalg(3)      ; dttyp(98)=8
 else                                        ! ===== nonlop_ylm
!  gx + gxfac + gxfac_sij
   cadd(94)=2*lmnmax*mincat*(mpw+1+usepaw)    ; dttyp(94)=8
!  kpg
   cadd(95)=3*mpw       ; dttyp(95)=8
!  indlmn_typ, ffnl_typ
   cadd(96)=lmnmax*6; dttyp(96)=4
!  ffnl_typ
   cadd(97)=lmnmax*mpw; dttyp(97)=8
!  opernla_ylm: scalar,scali
   cadd(98)=2*mpw; dttyp(98)=8
 end if

!(10)                    in suscep and suskmm ----------------------------

 if(modulo(iprcel,100)>=20.and.modulo(iprcel,100)<70)then
!  Variables allocated in suscep
   if(mkmem==0)then
!    cg_disk
     cmpw(101)=2*mband             ; dttyp(101)=8
   end if
   if(occopt>=3)then
!    drhode
     cadd(103)=2*npwdiel*nsppol    ; dttyp(103)=8
   end if
!  rhoextrap (always included, although it appears only when extrap==1)
   cadd(104)=ndiel456            ; dttyp(104)=8

!  Variables allocated in suskmm
!  cwavef
   cmpw(106)=2                   ; dttyp(106)=8
!  rhoaug, wfraug
   cadd(107)=3*ndiel456          ; dttyp(107)=8
!  wfprod
   cadd(108)=2*npwdiel           ; dttyp(108)=8
!  wfrspa1, wfrspa2
   cadd(109)=4*ndiel456*nbnd_in_blk ; dttyp(109)=8

 end if

!(11)                    in dielmt ---------------------------------------

 if(modulo(iprcel,100)>=20.and.modulo(iprcel,100)<70)then
!  dielh,dielvec,eig_diel,zhpev1,zhpev2
   cadd(111)=3*npwdiel*npwdiel                   &
&   +9*npwdiel           ; dttyp(111)=8
 end if

!(12)                    in tddft  ---------------------------------------

 if(iscf==-1)then
   if(mkmem/=0)then
!    cg_disk
     cmpw(121)=2*mband            ; dttyp(121)=8
   end if
!  cwavef
   cmpw(124)=2*mband             ; dttyp(124)=8
!  rhoaug,wfraug,wfrspa
   cadd(125)=(2+mband)*ndiel456  ; dttyp(125)=8
 end if

!--------------------------------------------------------------------------

 chain(:,:)=.true.

!Define the main chain version a (fourwf)
 chain(31:50,1)=.false.
 chain(60:70,1)=.false.
 chain(77:80,1)=.false.
 chain(93:100,1)=.false.
 chain(101:marrays,1)=.false.

!Define the main chain version b (nonlop+opernl)
 chain(31:50,2)=.false.
 chain(60:70,2)=.false.
 chain(77:80,2)=.false.
 chain(91:92,2)=.false.
 chain(101:marrays,2)=.false.

!Define the XC chain ( 31:40 belong only to this chain)
 chain(41:marrays,3)=.false.

!Define the mkrho chain ( 62:66 and 76:77 belong only to this chain)
!is it sure that they have to be summed ?)
 chain(31:50,4)=.false.
 chain(51:59,4)=.false.
 chain(61   ,4)=.false.
 chain(67:70,4)=.false.
 chain(71:marrays,4)=.false.
 chain(77:80,4)=.true.

!Define the fourdp chain ( 67:70 belong only to this chain)
 chain(31:50,5)=.false.
 chain(51:66,5)=.false.
 chain(60   ,5)=.true.
 chain(71:marrays,5)=.false.

!Define the parallel k-point chain ( 61 belong only to this chain )
 chain(31:50,6)=.false.
 chain(51:59,6)=.false.
 chain(62:70,6)=.false.
 chain(71:marrays,6)=.false.

!Define the newvtr chain ( 41:50 belong only to this chain)
 chain(31:40,7)=.false.
 chain(51:marrays,7)=.false.

!Define the suscep chain ( 101:110 belong only to this chain)
 chain(31:marrays,8)=.false.
 chain(60    ,8)=.true.
 chain(101:110,8)=.true.

!Define the dielmt chain ( 111:120 belong only to this chain)
 chain(31:marrays,9)=.false.
 chain(60    ,9)=.true.
 chain(111:120,9)=.true.

!Define the tddft chain ( 121:130 belong only to this chain)
 chain(31:marrays,10)=.false.
 chain(60    ,10)=.true.
 chain(121:130,10)=.true.

!The memory needed for each chain has been computed
!-------------------------------------------------------------------------
!Still need some auxiliary data : estimate the disk space
!or the maximum segment size.

!XG030513 : MPIWF need to multiply mbdiskwf by the number of processors
!in the WF group. For the time being, nprocwf=1
 nprocwf=mpi_enreg%nproc_fft

 mbdiskwf=(8*two*mpw*nprocwf*sum(nband(1:nkpt*nsppol)))/1024._dp**2 + 0.002_dp
 mbdiskpd=(8*nfftf*nsppol)/1024._dp**2 + 0.002_dp

!Determine the largest array out of cg (cg_disk), f_fftgr (f_fftgr_disk), or pawfgrtab%gylm
 if(mkmem==0)then
   mbcg=(8*2*mpw*mband)/1024._dp**2 + 0.002_dp
 else
   mbcg=(8*2*mpw*mband*mkmem*nsppol)/1024._dp**2 + 0.002_dp
 end if
 if(mffmem==0)then
   if (pawmixdg==1) then
     mbf_fftgr=(8*nfftf*n_fftgr)/1024._dp**2 + 0.002_dp
   else
     mbf_fftgr=(8*nfft*n_fftgr)/1024._dp**2 + 0.002_dp
   end if
 else
   if (pawmixdg==1) then
     mbf_fftgr=(8*nfftf*n_fftgr*nsppol*mffmem)/1024._dp**2 + 0.002_dp
   else
     mbf_fftgr=(8*nfft*n_fftgr*nsppol*mffmem)/1024._dp**2 + 0.002_dp
   end if
 end if
 if(usepaw==1)then
   mbgylm=0
   do ii=1,ntypat                                        ! pawfgrtab
     jj=(1+int(nfftf*four_pi/(three*ucvol)*rshp(ii)**3))
     mbgylm=mbgylm+my_nattyp(ii)*jj &
&     *( l_size(ii)**2*pawstgylm &                              !   %gylm   (if pawstgylm=1)
&    +3*max((optforces+1)/2,optstress)*l_size(ii)**2*pawstgylm& !   %gylmgr (if pawstgylm=1)
&    +3*optstress*pawstgylm&                                    !   %rfgd   (if pawstgylm=1)
&    +3*(1-pawstgylm) )                                         !   %rfgd   (if pawstgylm=0)
   end do
   mbgylm=8*mbgylm/1024._dp**2 + 0.002_dp
 else
   mbgylm=0
 end if

!-------------------------------------------------------------------------
 ABI_FREE(my_nattyp)
 ABI_FREE(basis_size)
 ABI_FREE(l_size)
 ABI_FREE(lmn_size)
 ABI_FREE(lmn2_size)
 ABI_FREE(mesh_size)
 ABI_FREE(pawver)
 ABI_FREE(shape_type)
 ABI_FREE(rshp)

!---------------------------------------------------------------------
!Now, analyze the data

 call memana(cadd,cfft,cfftf,chain,cmpw,dttyp,iout,iprcel,iscf,&
& marrays,mbcg,mbdiskpd,mbdiskwf,mbf_fftgr,mbgylm,mffmem,&
& mpw,natom,nchain,nfft,nfftf,occopt,option,prtvol)

end subroutine memory
!!***

!!****f* m_memeval/memana
!! NAME
!! memana
!!
!! FUNCTION
!! Analysis of the memory and disk space needed for the job,
!! thanks to the data computed in the calling routine: for each
!! array, the number of blocks of size mpw or nfft bytes, and the
!! additional memory occupation;
!! the list of arrays that are used for each chain.
!!
!! According to the value of the option variable,
!! the routine will eventually try to allocate this amount of memory,
!! and if it fails, estimate the maximum value nfft compatible with
!! the available memory.
!!
!! INPUTS
!!  cadd(marrays)= count of bytes needed in addition of cmpw, cfftc and cfft.
!!  cfft(marrays) =for each array, count of blocks of size nfft bytes (coarse grid, if PAW)
!!  cfftf(marrays)=for each array, count of blocks of size nfft bytes (fine grid, if PAW)
!!  chain(marrays,nchain)=logical variable, that informs whether an array
!!    belongs to a given chain.
!!  cmpw(marrays)=for each array, count of blocks of size mpw bytes.
!!  dttyp(marrays)=datatype of the array : 4 for integers, 8 for real(dp)
!!  iout=unit number for output of formatted data.
!!  iprcel=govern the choice of preconditioner for the SCF cycle
!!  iscf=governs the choice of SCF algorithm, or non-SCF calculation.
!!  marrays=maximal number of arrays (or group of arrays) to be monitored.
!!  mbcg=number of MB needed for the cg array.
!!  mbdiskpd=number of MB needed to store a density or potential file on disk
!!  mbdiskwf=number of MB needed to store a wavefunction file on disk
!!  mbf_fftgr=number of MB needed for the f_fftgr array.
!!  mbgylm=number of MB needed for the pawfgrtab%gylm array (paw only)
!!  mffmem =governs the number of FFT arrays which are fit in core memory
!!  mpw   =maximum number of planewaves in basis sphere (large number)
!!  natom =number of atoms in unit cell
!!  nchain=number of chains to be used in the estimation of memory.
!!  nfft =(effective) number of FFT grid points (for one processor) (coarse grid, if PAW)
!!  nfftf=(effective) number of FFT grid points (for one processor) (fine grid, if PAW)
!!  occopt=option for occupation numbers. If 3<=occopt<=8, varying occupation
!!  option : if 0 , no test of available memory
!!           if 1 , the routine tries to allocate the estimated memory, for testing
!!                    purposes, and if a failure occurs, the routine stops.
!!           if 2 , like 1, but before stopping, the routine will provide
!!                    an estimation of the available memory.
!!  prtvol=control print volume
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      m_memeval
!!
!! CHILDREN
!!      atomic_info,createwavefunctionsdescriptors,deallocate_lr
!!      memoryestimator,mkradim,wrtout,wvl_descr_atoms_set,wvl_descr_free
!!      wvl_setboxgeometry,xred2xcart
!!
!! SOURCE

subroutine memana(cadd,cfft,cfftf,chain,cmpw,dttyp,iout,iprcel,iscf,&
& marrays,mbcg,mbdiskpd,mbdiskwf,mbf_fftgr,mbgylm,mffmem,&
& mpw,natom,nchain,nfft,nfftf,occopt,option,prtvol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,iprcel,iscf,marrays,mffmem,mpw,natom,nchain
 integer,intent(in) :: nfft,nfftf,occopt,option,prtvol
 real(dp),intent(in) :: mbcg,mbdiskpd,mbdiskwf,mbf_fftgr,mbgylm
!arrays
 integer,intent(in) :: dttyp(marrays)
 logical,intent(in) :: chain(marrays,nchain)
 real(dp),intent(in) :: cadd(marrays),cfft(marrays),cfftf(marrays),cmpw(marrays)

!Local variables-------------------------------
!scalars
 integer :: biggest,ichain,ier,ier1,ier2,ier3,ier4,ier5,ier6,ier7,ier8,ii
!integer :: jj,kk
 integer :: mu,nmbytes,nquarter_mbytes,quit
 real(dp) :: mbbigarr,mbbiggest
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: bigarray(:,:),bigarray1(:,:),bigarray2(:,:)
 real(dp),allocatable :: bigarray3(:,:),bigarray4(:,:),bigarray5(:,:)
 real(dp),allocatable :: bigarray6(:,:),bigarray7(:,:),bigarray8(:,:)
 real(dp),allocatable :: cdpadd(:),cdpfft(:),cdpfftf(:),cdpmpw(:)
 real(dp),allocatable :: cintfft(:),cintfftf(:),cintmpw(:),cintadd(:)
 real(dp),allocatable :: mbdpadd(:),mbdpfft(:),mbdpfftf(:)
 real(dp),allocatable :: mbdpmpw(:),mbintadd(:),mbintfft(:),mbintfftf(:)
 real(dp),allocatable :: mbintmpw(:),mbother(:),mbtot(:)

! **************************************************************************

!write(std_out,*)' memana : nchain=',nchain

 ABI_MALLOC(cdpfftf,(nchain))
 ABI_MALLOC(cdpfft,(nchain))
 ABI_MALLOC(cdpmpw,(nchain))
 ABI_MALLOC(cintfftf,(nchain))
 ABI_MALLOC(cintfft,(nchain))
 ABI_MALLOC(cintmpw,(nchain))
 ABI_MALLOC(cdpadd,(nchain))
 ABI_MALLOC(cintadd,(nchain))
 ABI_MALLOC(mbdpadd,(nchain))
 ABI_MALLOC(mbdpfftf,(nchain))
 ABI_MALLOC(mbdpfft,(nchain))
 ABI_MALLOC(mbdpmpw,(nchain))
 ABI_MALLOC(mbintadd,(nchain))
 ABI_MALLOC(mbintfftf,(nchain))
 ABI_MALLOC(mbintfft,(nchain))
 ABI_MALLOC(mbintmpw,(nchain))
 ABI_MALLOC(mbother,(nchain))
 ABI_MALLOC(mbtot,(nchain))

 biggest=0
 mbbiggest=0.0_dp

!For each chain, compute the number of bytes
 do ichain=1,nchain

!  First, the number of integer or real(dp), fft, mpw or add blocks
   cdpmpw(ichain) =sum(cmpw(:),MASK=(dttyp(:)==8).and.chain(:,ichain))
   cintmpw(ichain)=sum(cmpw(:),MASK=(dttyp(:)==4).and.chain(:,ichain))
   cdpfftf(ichain) =sum(cfftf(:),MASK=(dttyp(:)==8).and.chain(:,ichain))
   cintfftf(ichain)=sum(cfftf(:),MASK=(dttyp(:)==4).and.chain(:,ichain))
   cdpfft(ichain) =sum(cfft(:),MASK=(dttyp(:)==8).and.chain(:,ichain))
   cintfft(ichain)=sum(cfft(:),MASK=(dttyp(:)==4).and.chain(:,ichain))
   cdpadd(ichain) =sum(cadd(:),MASK=(dttyp(:)==8).and.chain(:,ichain))
   cintadd(ichain)=sum(cadd(:),MASK=(dttyp(:)==4).and.chain(:,ichain))

!  Compute the corresponding number of Mbytes
   mbdpmpw(ichain) =8*cdpmpw(ichain) *dble(mpw) /1024._dp**2
   mbintmpw(ichain)=4*cintmpw(ichain)*dble(mpw) /1024._dp**2
   mbdpfftf(ichain) =8*cdpfftf(ichain) *dble(nfftf)/1024._dp**2
   mbintfftf(ichain)=4*cintfftf(ichain)*dble(nfftf)/1024._dp**2
   mbdpfft(ichain) =8*cdpfft(ichain) *dble(nfft)/1024._dp**2
   mbintfft(ichain)=4*cintfft(ichain)*dble(nfft)/1024._dp**2
   mbdpadd(ichain) =8*cdpadd(ichain)              /1024._dp**2
   mbintadd(ichain)=4*cintadd(ichain)             /1024._dp**2
   mbother(ichain) =dble(231+6*natom)/1024._dp
   if(3<=occopt .and. occopt<=8)mbother(ichain)=dble(991+natom)/1024._dp

!  Compute the total number of Mbytes
   mbtot(ichain)=mbdpmpw(ichain)+mbintmpw(ichain)&
&   +mbdpfftf(ichain)+mbintfftf(ichain)&
&   +mbdpfft(ichain)+mbintfft(ichain)&
&   +mbdpadd(ichain)+mbintadd(ichain)+mbother(ichain)

!  Select the biggest chain
   if(mbtot(ichain)>mbbiggest)then
     mbbiggest=mbtot(ichain)
     biggest=ichain
   end if
 end do
!When iprcel<20, the biggest chains cannot be number 8 or 9 ...
 if(modulo(iprcel,100)<20 .and. (biggest==8 .or. biggest==9))then
   write(msg,'(a,a,a,a,i3,a,a,a)') ch10,&
&   ' memana: BUG -',ch10,&
&   '  The biggest chain is number',biggest,' while iprcel==20.',ch10,&
&   '  This is not allowed.'
   call wrtout(std_out,msg,'COLL')
 end if

 write(msg, '(a,f11.3,a)' ) &
& 'P This job should need less than                 ',&
& mbbiggest+tol10,' Mbytes of memory. '
 call wrtout(std_out,msg,'COLL')
 call wrtout(iout,msg,'COLL')

 if(prtvol>=10)then
   if(biggest==1)write(msg,'(a)')'P Max. in main chain + fourwf.f '
   if(biggest==2)write(msg,'(a)')'P Max. in main chain + nonlop.f + opernl.f '
   if(biggest==3)write(msg,'(a)')'P Max. in XC chain '
   if(biggest==4)write(msg,'(a)')'P Max. in mkrho chain '
   if(biggest==5)write(msg,'(a)')'P Max. in fourdp chain '
   if(biggest==6)write(msg,'(a)')'P Max. in parallel k-point chain '
   if(biggest==7)write(msg,'(a)')'P Max. in newvtr chain '
   if(biggest==8)write(msg,'(a)')'P Max. in suscep chain '
   if(biggest==9)write(msg,'(a)')'P Max. in dielmt chain '
   if(biggest==10)write(msg,'(a)')'P Max. in tddft chain '
   call wrtout(iout,msg,'COLL')

   write(msg, '(a,i13,a,f11.3,a)' )&
&   'P',nint(cintmpw(biggest)),' blocks of mpw  integer numbers, for',&
&   mbintmpw(biggest)+tol10,' Mbytes. '
   call wrtout(iout,msg,'COLL')
   write(msg, '(a,i13,a,f11.3,a)' )&
&   'P',nint(cdpmpw(biggest)),' blocks of mpw  real(dp)  numbers, for',&
&   mbdpmpw(biggest)+tol10,' Mbytes. '
   call wrtout(iout,msg,'COLL')
   if (nfft==nfftf) then
     if(mbintfft(biggest)+mbintfftf(biggest)>0.001)then
       write(msg, '(a,i13,a,f11.3,a)' )&
&       'P',nint(cintfft(biggest)+cintfftf(biggest)),' blocks of nfft integer numbers, for',&
&       mbintfft(biggest)+mbintfftf(biggest)+tol10,' Mbytes. '
       call wrtout(iout,msg,'COLL')
     end if
     write(msg, '(a,i13,a,f11.3,a)' )&
&     'P',nint(cdpfft(biggest)+cdpfftf(biggest)),' blocks of nfft real(dp)  numbers, for',&
&     mbdpfft(biggest)+mbdpfftf(biggest)+tol10,' Mbytes. '
     call wrtout(iout,msg,'COLL')
   else
     if(mbintfftf(biggest)>0.001)then
       write(msg, '(a,i13,a,f11.3,a)' )&
&       'P',nint(cintfftf(biggest)),' blocks of nfft (fine grid) integer numbers, for',&
&       mbintfftf(biggest)+tol10,' Mbytes. '
       call wrtout(iout,msg,'COLL')
     end if
     write(msg, '(a,i13,a,f11.3,a)' )&
&     'P',nint(cdpfftf(biggest)),' blocks of nfft (fine grid) real(dp)  numbers, for',&
&     mbdpfftf(biggest)+tol10,' Mbytes. '
     call wrtout(iout,msg,'COLL')
     if(mbintfft(biggest)>0.001)then
       write(msg, '(a,i13,a,f11.3,a)' )&
&       'P',nint(cintfft(biggest)),' blocks of nfft (coarse grid) integer numbers, for',&
&       mbintfft(biggest)+tol10,' Mbytes. '
       call wrtout(iout,msg,'COLL')
     end if
     write(msg, '(a,i13,a,f11.3,a)' )&
&     'P',nint(cdpfft(biggest)),' blocks of nfft (coarse grid) real(dp)  numbers, for',&
&     mbdpfft(biggest)+tol10,' Mbytes. '
     call wrtout(iout,msg,'COLL')
   end if
   if(mbintadd(biggest)>0.001)then
     write(msg, '(a,13x,a,f11.3,a)' )'P',' Additional     integer numbers, for',mbintadd(biggest)+tol10,' Mbytes. '
     call wrtout(iout,msg,'COLL')
   end if
   write(msg, '(a,13x,a,f11.3,a)' )'P',' Additional     real(dp)  numbers, for',mbdpadd(biggest)+tol10,' Mbytes. '
   call wrtout(iout,msg,'COLL')
   write(msg, '(a,13x,a,f11.3,a)' )'P',' With residue estimated to be       ',mbother(biggest)+tol10,' Mbytes. '
   call wrtout(iout,msg,'COLL')
   write(msg, '(a)' )'P'
   call wrtout(iout,msg,'COLL')
   write(msg, '(a)' )'P Comparison of the memory needs of different chains'
   call wrtout(iout,msg,'COLL')

   write(msg, '(a,f11.3,a)' )'P Main chain + fourwf.f           ',mbtot(1)+tol10,' Mbytes. '
   call wrtout(iout,msg,'COLL')
   write(msg, '(a,f11.3,a)' )'P Main chain + nonlop.f + opernl.f',mbtot(2)+tol10,' Mbytes. '
   call wrtout(iout,msg,'COLL')

!  The next chains are not defined in the RF case.
   if(nchain>2)then
     write(msg, '(a,f11.3,a)' )'P XC chain                        ',mbtot(3)+tol10,' Mbytes. '
     call wrtout(iout,msg,'COLL')
     write(msg, '(a,f11.3,a)' )&
&     'P mkrho chain                     ',mbtot(4)+tol10,' Mbytes. '
     call wrtout(iout,msg,'COLL')
     write(msg, '(a,f11.3,a)' )&
&     'P fourdp chain                    ',mbtot(5)+tol10,' Mbytes. '
     call wrtout(iout,msg,'COLL')
     if(xmpi_paral==1)then
       write(msg, '(a,f11.3,a)' )&
&       '- parallel k-point chain          ',mbtot(6)+tol10,' Mbytes. '
       call wrtout(iout,msg,'COLL')
     end if
     write(msg, '(a,f11.3,a)' )&
&     'P newvtr chain                    ',mbtot(7)+tol10,' Mbytes. '
     call wrtout(iout,msg,'COLL')
     if(modulo(iprcel,100)>=20.and.modulo(iprcel,100)<70)then
       write(msg, '(a,f11.3,a)' )&
&       'P suscep chain                    ',mbtot(8)+tol10,' Mbytes. '
       call wrtout(iout,msg,'COLL')
       write(msg, '(a,f11.3,a)' )&
&       'P dielmt chain                    ',mbtot(9)+tol10,' Mbytes. '
       call wrtout(iout,msg,'COLL')
     end if
     if(iscf==-1)then
       write(msg, '(a,f11.3,a)' )&
&       'P tddft  chain                    ',mbtot(10)+tol10,' Mbytes. '
     end if
   end if ! nchain>2

 end if

!--------------------------------------------------------------------

 write(msg, '(a)' ) '  Rough estimation (10% accuracy) of disk space for files :'
 call wrtout(iout,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 write(msg, '(a,f11.3,a,a,f11.3,a)' ) &
& '_ WF disk file :',mbdiskwf+tol10,' Mbytes ;',&
& ' DEN or POT disk file :',mbdiskpd+tol10,' Mbytes.'
 call wrtout(iout,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 if(mffmem==0 .and. iscf>0)then
   if(iscf==1)then
     write(msg, '(a,a,a)' )&
&     '  mffmem==0, iscf==1 => use of 1 FFT temporary disk file,',ch10,&
&     '                       5 times bigger than a DEN file.'
   else if(iscf==2.or.iscf==12)then
     write(msg, '(a,a,a)' )&
&     '  mffmem==0, iscf==2 => use of 1 FFT temporary disk file,',ch10,&
&     '                       3 times bigger than a DEN file.'
   else if(iscf==3.or.iscf==13)then
     write(msg, '(a,a,a)' )&
&     '  mffmem==0, iscf==3 => use of 1 FFT temporary disk file,',ch10,&
&     '                       4 times bigger than a DEN file.'
   else if(iscf==4.or.iscf==14)then
     write(msg, '(a,a,a)' )&
&     '  mffmem==0, iscf==4 => use of 1 FFT temporary disk file,',ch10,&
&     '                       6 times bigger than a DEN file.'
   else if(iscf==5)then
     write(msg, '(a,a,a)' )&
&     '  mffmem==0, iscf==5 => use of 1 FFT temporary disk file,',ch10,&
&     '                       10 times bigger than a DEN file.'
   else if(iscf==6)then
     write(msg, '(a,a,a)' )&
&     '  mffmem==0, iscf==6 => use of 1 FFT temporary disk file,',ch10,&
&     '                       10 times bigger than a DEN file.'
   else if(iscf==7.or.iscf==17)then
     write(msg, '(a,a,a)' )&
&     '  mffmem==0, iscf==7 => use of 1 FFT temporary disk file,',ch10,&
&     '                       (2+2*npulayit) times bigger than a DEN file.'
   end if
   call wrtout(iout,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
 end if

!Temporary msg - estimation of PAW specific data has to be done...
!Have to add the usepaw argument to use this.
!if (usepaw==1) then
!write(msg,'(5a)') '  WARNING: You are using PAW formalism;',ch10,&
!&       '           Above estimations do not take PAW',ch10,&
!&       '           specific data into account !'
!call wrtout(iout,msg,'COLL')
!call wrtout(std_out,msg,'COLL')
!end if

 write(msg,'(80a,a)') ('=',mu=1,80),ch10
 call wrtout(iout,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

!--------------------------------------------------------------------
!Here, each processor must test its memory, so use
!the PERS mode for error msgs, followed by synchronisation

 mbbigarr=max(mbf_fftgr,mbcg,mbgylm)
 if(mbbigarr==mbcg) then
   write(msg, '(a,f12.4,a)' ) ' Biggest array : cg(disk), with',mbcg+tol10,' MBytes.'
 else if (mbbigarr==mbf_fftgr) then
   write(msg, '(a,f12.4,a)' ) ' Biggest array : f_fftgr(disk), with',mbf_fftgr+tol10,' MBytes.'
 else if (mbbigarr==mbgylm)then
   write(msg, '(a,f12.4,a)' ) ' Biggest array : pawfgrtab%gylm(gr), with',mbgylm+tol10,' MBytes.'
 end if
 call wrtout(std_out,msg,'COLL')

!if (mpi_enreg%my_nimage>1) then
!write(msg, '(a,f12.4,a)' ) &
!&   ' These estimations take the distribution over replicas (images) of the cell into account.'
!call wrtout(std_out,msg,'COLL')
!end if

 quit=0

 if(option>=1)then

!  Test the ability to allocate the biggest array
   nquarter_mbytes=4.0_dp*mbbigarr+1.0_dp
   ABI_STAT_MALLOC(bigarray,(32*1024,nquarter_mbytes), ier)
   if(ier/=0)then
     write(msg,'(a,f11.3,a,a,a,a,a,a,a)')&
&     'Test failed to allocate an array of',mbbigarr,' Mbytes',ch10,&
&     'It is not worth to continue ',ch10,&
&     'Action: modify input variable to fit the available memory,',ch10,&
&     'increase limit on maximal array size or set mem_test to 0 to disable this test.'
     call wrtout(std_out,msg,'PERS')
     if(option==1)then
       ABI_ERROR_CLASS(msg, "MemanaError")
     else
       ABI_WARNING(msg)
       quit=1
     end if
   end if
   if(quit==0)then
     write(msg,'(a,f11.3,a)')' memana : allocated an array of',mbbigarr+tol10,' Mbytes, for testing purposes. '
     call wrtout(std_out,msg,'COLL')
   end if
   if(allocated(bigarray)) then
     ABI_FREE(bigarray)
   end if

!  Test the ability to allocate the needed total memory : use 8 segments,
!  hoping that the maximal segment size is not so much smaller than the
!  total memory
   nquarter_mbytes=0.5_dp*mbbiggest+1.0_dp
   ABI_STAT_MALLOC(bigarray1,(32*1024,nquarter_mbytes), ier1)
   ABI_STAT_MALLOC(bigarray2,(32*1024,nquarter_mbytes), ier2)
   ABI_STAT_MALLOC(bigarray3,(32*1024,nquarter_mbytes), ier3)
   ABI_STAT_MALLOC(bigarray4,(32*1024,nquarter_mbytes), ier4)
   ABI_STAT_MALLOC(bigarray5,(32*1024,nquarter_mbytes), ier5)
   ABI_STAT_MALLOC(bigarray6,(32*1024,nquarter_mbytes), ier6)
   ABI_STAT_MALLOC(bigarray7,(32*1024,nquarter_mbytes), ier7)
   ABI_STAT_MALLOC(bigarray8,(32*1024,nquarter_mbytes), ier8)

   if(ier1/=0 .or. ier2/=0 .or. ier3/=0 .or. ier4/=0 .or. ier5/=0 .or. ier6/=0 .or. ier7/=0 .or. ier8/=0) then
     write(msg,'(a,f11.3,a,a,a,a,a,a,a)')&
&     'Test failed to allocate ',mbbiggest,' Mbytes',ch10,&
&     'It is not worth to continue ',ch10,&
&     'Action: modify input variables or submission parameters to fit the available memory,',ch10,&
&     'increase limit on available memory or set mem_test to 0 to disable this test.'
     if(option==1)then
       ABI_ERROR_CLASS(msg, "MemanaError")
     else
       ABI_WARNING(msg)
       quit=1
     end if
   end if

   if(quit==0)then
     write(msg,'(a,f11.3,a,a,a)')&
&     ' memana: allocated ',mbbiggest,'Mbytes, for testing purposes. ',ch10,&
&     ' The job will continue.'
     call wrtout(std_out,msg,'COLL')
   end if
   if(allocated(bigarray1)) then
     ABI_FREE(bigarray1)
   end if
   if(allocated(bigarray2)) then
     ABI_FREE(bigarray2)
   end if
   if(allocated(bigarray3)) then
     ABI_FREE(bigarray3)
   end if
   if(allocated(bigarray4)) then
     ABI_FREE(bigarray4)
   end if
   if(allocated(bigarray5)) then
     ABI_FREE(bigarray5)
   end if
   if(allocated(bigarray6)) then
     ABI_FREE(bigarray6)
   end if
   if(allocated(bigarray7)) then
     ABI_FREE(bigarray7)
   end if
   if(allocated(bigarray8)) then
     ABI_FREE(bigarray8)
   end if

 end if

!--------------------------------------------------------------------

 if(option==2 .and. quit==1 )then

!  Estimation of the available memory
!
!  A quarter of Mbyte is 256*1024/8 real(dp) numbers,
!  that is 32*1024 dp numbers.
!  One begins with the allocation of 4 Mbytes. If successful,
!  one increases that number, until the allocation is not successfull
!  any more. Unfortunately, on a P6 with the pghpf compiler, the
!  allocate instruction generate a core dump, instead of returning
!  an error code, so that this part of code has been made optional.

   nquarter_mbytes=16
   nmbytes=nquarter_mbytes/4.0_dp

!  With an increase ratio of 1.25_dp (see below), ii=5 leads to 9 MB,
!  ii=10 leads to 28 MB, ii=15 leads to 85 MB, ii=18 leads to 165 MB,
!  ii=30 is over 2 GB
   do ii=1,30
     ABI_STAT_MALLOC(bigarray,(32*1024,nquarter_mbytes), ier)
     if(ier/=0)then
       write(msg,'(a,i0,a)')' memana : failed to allocate ',nmbytes,' Mbytes'
       call wrtout(std_out,msg,'PERS')
       exit
     end if
     write(msg,'(a,i0,a)')' memana : succeeded to allocate ',nmbytes,' Mbytes'
     call wrtout(std_out,msg,'PERS')
!    Here really test the space
!    do kk=1,nquarter_mbytes
!    do jj=1,32*1024,37
!    bigarray(jj,kk)=0.0_dp
!    end do
!    write(std_out,*)' memana : wrote ',kk,' quarter of mbytes'
!    end do
     ABI_FREE(bigarray)
     nquarter_mbytes=dble(nquarter_mbytes)*1.25_dp
     nmbytes=nquarter_mbytes/4.0_dp
   end do
   if(allocated(bigarray)) then
     ABI_FREE(bigarray)
   end if

   ABI_ERROR_CLASS("in memana with option==2 .and. quit==1", "MemanaError")
 end if !  End the test of the available memory

!--------------------------------------------------------------------

 ABI_FREE(cdpfftf)
 ABI_FREE(cdpfft)
 ABI_FREE(cdpmpw)
 ABI_FREE(cintfftf)
 ABI_FREE(cintfft)
 ABI_FREE(cintmpw)
 ABI_FREE(cdpadd)
 ABI_FREE(cintadd)
 ABI_FREE(mbdpadd)
 ABI_FREE(mbdpfftf)
 ABI_FREE(mbdpfft)
 ABI_FREE(mbdpmpw)
 ABI_FREE(mbintadd)
 ABI_FREE(mbintfftf)
 ABI_FREE(mbintfft)
 ABI_FREE(mbintmpw)
 ABI_FREE(mbother)
 ABI_FREE(mbtot)

end subroutine memana
!!***

!!****f* m_memeval/memorf
!! NAME
!! memorf
!!
!! FUNCTION
!! Estimation of the memory needed for a response-function job.
!! According to the value of the option variable,
!! might also try to allocate this amount of memory, and if it fails,
!! might estimate the available memory.
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
!!  mpi_enreg=information about MPI parallelization
!!  mpssang is 1+maximum angular momentum for nonlocal pseudopotential
!!  mpssoang is 1+maximum (spin*angular momentum) for nonlocal pseudopotential
!!  mpw   =maximum number of planewaves in basis sphere (large number)
!!  mqgrid=maximum dimension of grid of q values for psp representations
!!  natom =number of atoms in unit cell
!!  nband(nkpt*nsppol)=number of bands at each k point, for each polarization
!!  nfft  =(effective) number of FFT grid points (for one processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
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
!! in rhotoxc and daughter routines (at maximum 22*1000 dp numbers),
!! as well as other arrays like
!! character(len=500) :: msg (present in about 100 routines), or the different
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
!!      m_memeval
!!
!! CHILDREN
!!      atomic_info,createwavefunctionsdescriptors,deallocate_lr
!!      memoryestimator,mkradim,wrtout,wvl_descr_atoms_set,wvl_descr_free
!!      wvl_setboxgeometry,xred2xcart
!!
!! SOURCE

subroutine memorf(cplex,n1xccc,getcell,idtset,intxc,iout,iprcel,&
& iscf,jdtset,lmnmax,lnmax,mband,mffmem,mgfft,&
& mkmems,mpi_enreg,mpsang,mpssoang,mpw,mqgrid,&
& natom,nband,nfft,ngfft,&
& nkpt,nloalg,nspden,nspinor,nsppol,nsym,ntypat,&
& occopt,optddk,optphon,option,optstrs,prtvol,useylm,use_gpu_cuda,xclevel)

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
 character(len=500) :: msg
 character(len=1) :: firstchar
!arrays
 integer :: dttyp(marrays)
 real(dp) :: cadd(marrays),cfft(marrays),cmpw(marrays)
 real(dp),allocatable :: cfft_dum(:)
 logical :: chain(marrays,nchain)

! **************************************************************************

 if(option<0 .or. option>2)then
   write(msg, '(a,i0,a)')'option= ',option,' while the only allowed values are 0, 1, or 2.'
   ABI_BUG(msg)
 end if

 firstchar=' ';if (use_gpu_cuda==1) firstchar='_'
 cmpw(:)=zero ; cfft(:)=zero ; cadd(:)=zero
 dttyp(:)=0

 call wrtout(std_out,' memorf : analysis of memory needs ','COLL')

 if(jdtset>=100)then
   write(msg,'(80a,a,a,i5,a)')('=',mu=1,80),ch10,&
   ' Values of the parameters that define the memory need for DATASET',jdtset,' (RF).'
 else if(jdtset/=0)then
   write(msg,'(80a,a,a,i3,a)')('=',mu=1,80),ch10,&
   ' Values of the parameters that define the memory need for DATASET',jdtset,' (RF).'
 else
   write(msg,'(80a,a,a,a)')('=',mu=1,80),ch10,&
   ' Values of the parameters that define the memory need of the present run',' (RF).'
 end if
 call wrtout(iout,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 mkmem=mkmems(1)
 mkqmem=mkmems(2)
 mk1mem=mkmems(3)
 my_natom=natom;if (mpi_enreg%nproc_atom>1) my_natom=mpi_enreg%my_natom

 write(msg,'( 4(a,i8),a,4(a,i8) )' ) &
& '     intxc =',intxc   ,'      iscf =',iscf,&
& '    lmnmax =',lmnmax  ,'     lnmax =',lnmax,ch10,&
& '     mgfft =',mgfft,'  mpssoang =',mpssoang,&
& '    mqgrid =',mqgrid,'     natom =',natom
 call wrtout(iout,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 write(msg,'( 4(a,i8),a,4(a,i8),a,4(a,i8) )' ) &
& '  nloc_mem =',nloalg(2)*(nloalg(3)+1),'    nspden =',nspden ,&
& '   nspinor =',nspinor,'    nsppol =',nsppol ,ch10,&
& '      nsym =',nsym,'    n1xccc =',n1xccc ,&
& '    ntypat =',ntypat,'    occopt =',occopt ,ch10,&
& '   xclevel =',xclevel
 call wrtout(iout,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 write(msg,'(4(3(a,i12),a))') &
& '-    mband =',mband  ,'        mffmem =',mffmem,&
& '         mkmem =',mkmem  ,ch10,&
& '-   mkqmem =',mkqmem ,'        mk1mem =',mk1mem,&
& '           mpw =',mpw  ,ch10,&
& '      nfft =',nfft ,'          nkpt =',nkpt
 call wrtout(iout,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 if (my_natom/=natom)then
   write(msg,'(a,i10)') 'Pmy_natom=',my_natom
   call wrtout(iout,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
 end if

 write(msg,'(80a)') ('=',mu=1,80)
 call wrtout(iout,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 if(getcell>0 .or. (getcell<0 .and. idtset+getcell>0) )then
   write(msg,'(a,a,a,a,a,a,i3,a,i3,a,a,a,a,a,a)' )ch10,&
&   ' memorf : COMMENT -',ch10,&
&   '  The determination of memory needs at this stage is meaningless,',ch10,&
&   '  since getcell = ',getcell,' is non-zero, while idtset=',idtset,'.',ch10,&
&   '  The following numbers are obtained by supposing that acell and rprim',ch10,&
&   '  are NOT taken from a previous dataset. You cannot rely on them.',ch10
   call wrtout(iout,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
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
!ghc,gvnlxc,gvnlx1
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

!gh1, gh_direc, gvnlx_direc, conjgr, direc, vresid, cwaveq
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

 ABI_MALLOC(cfft_dum,(marrays))
 cfft_dum=zero
 mbgylm=zero
 call memana(cadd,cfft,cfft_dum,chain,cmpw,dttyp,iout,iprcel,iscf,&
& marrays,mbcg,mbdiskpd,mbdiskwf,mbf_fftgr,mbgylm,mffmem,&
& mpw,natom,nchain,nfft,nfft,occopt,option,prtvol)
 ABI_FREE(cfft_dum)

end subroutine memorf
!!***

!!****f* m_memeval/getdim_nloc
!! NAME
!! getdim_nloc
!!
!! FUNCTION
!! Determine the dimensions of arrays that contain
!! the definition of non-local projectors : ekb, ffspl, indlmn
!!
!! INPUTS
!!  mixalch(npspalch,ntypalch,nimage)=alchemical mixing coefficients
!!  nimage=number of images
!!  npsp=number of pseudopotentials
!!  npspalch=number of pseudopotentials for alchemical purposes
!!  ntypat=number of types of pseudo atoms
!!  ntypalch=number of types of alchemical pseudo atoms
!!  pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file headers, as well as the psp file names
!!
!! OUTPUT
!!  lmnmax=maximum number of l,m,n projectors, not taking into account the spin-orbit
!!  lmnmaxso=maximum number of l,m,n projectors, taking into account the spin-orbit
!!  lnmax=maximum number of l,n projectors, not taking into account the spin-orbit
!!  lnmaxso=maximum number of l,n projectors, taking into account the spin-orbit
!!
!! PARENTS
!!      m_memeval,m_psps
!!
!! CHILDREN
!!      atomic_info,createwavefunctionsdescriptors,deallocate_lr
!!      memoryestimator,mkradim,wrtout,wvl_descr_atoms_set,wvl_descr_free
!!      wvl_setboxgeometry,xred2xcart
!!
!! SOURCE

subroutine getdim_nloc(lmnmax,lmnmaxso,lnmax,lnmaxso,mixalch,nimage,npsp,npspalch,&
& ntypat,ntypalch,pspheads)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nimage,npsp,npspalch,ntypalch,ntypat
 integer,intent(out) :: lmnmax,lmnmaxso,lnmax,lnmaxso
!arrays
 real(dp),intent(in) :: mixalch(npspalch,ntypalch,nimage)
 type(pspheader_type),intent(in) :: pspheads(npsp)

!Local variables-------------------------------
!scalars
 integer :: ilang,ipsp,ipspalch,itypalch,itypat,ntyppure
!integer :: llmax
 character(len=500) :: msg
!arrays
 integer,allocatable :: lmnproj_typat(:),lmnprojso_typat(:),lnproj_typat(:)
 integer,allocatable :: lnprojso_typat(:),nproj_typat(:,:),nprojso_typat(:,:)

! *************************************************************************

!write(std_out,*)' getdim_nloc: 'pspheads(1)%nproj(0:3)=',pspheads(1)%nproj(0:3)

 ABI_MALLOC(lmnproj_typat,(ntypat))
 ABI_MALLOC(lmnprojso_typat,(ntypat))
 ABI_MALLOC(lnproj_typat,(ntypat))
 ABI_MALLOC(lnprojso_typat,(ntypat))
 ABI_MALLOC(nproj_typat,(0:3,ntypat))
 ABI_MALLOC(nprojso_typat,(3,ntypat))
 lmnproj_typat(:)=0 ; lmnprojso_typat(:)=0
 lnproj_typat(:)=0 ; lnprojso_typat(:)=0
 nproj_typat(:,:)=0 ; nprojso_typat(:,:)=0

 ntyppure=ntypat-ntypalch

!For each type of pseudo atom, compute the number of projectors
!First, pure pseudo atoms
 if(ntyppure>0)then
   do itypat=1,ntyppure
     nproj_typat(0:3,itypat)=pspheads(itypat)%nproj(0:3)
     nprojso_typat(:,itypat)=pspheads(itypat)%nprojso(:)
   end do
 end if

!Then, alchemical pseudo atoms
 if(ntypalch>0)then
   do itypat=ntyppure+1,ntypat
     itypalch=itypat-ntyppure
     do ipsp=ntyppure+1,npsp
       ipspalch=ipsp-ntyppure
!      If there is some mixing, must accumulate the projectors
       if(sum(abs(mixalch(ipspalch,itypalch,:)))>tol10)then
         nproj_typat(0:3,itypat)=nproj_typat(0:3,itypat)+pspheads(ipsp)%nproj(0:3)
         nprojso_typat(:,itypat)=nprojso_typat(:,itypat)+pspheads(ipsp)%nprojso(:)
       end if
     end do
   end do
 end if

!Now that the number of projectors is known, accumulate the dimensions
 do itypat=1,ntypat
   do ilang=0,3
     lnproj_typat(itypat)=lnproj_typat(itypat)+nproj_typat(ilang,itypat)
     lmnproj_typat(itypat)=lmnproj_typat(itypat)+nproj_typat(ilang,itypat)*(2*ilang+1)
   end do
   lnprojso_typat(itypat)=lnproj_typat(itypat)
   lmnprojso_typat(itypat)=lmnproj_typat(itypat)
   do ilang=1,3
     lnprojso_typat(itypat)=lnprojso_typat(itypat)+nprojso_typat(ilang,itypat)
     lmnprojso_typat(itypat)=lmnprojso_typat(itypat)+nprojso_typat(ilang,itypat)*(2*ilang+1)
   end do
 end do

!Compute the maximal bounds, at least equal to 1, even for local psps
 lmnmax=1;lmnmaxso=1;lnmax=1;lnmaxso=1
 do itypat=1,ntypat
   lmnmax  =max(lmnmax  ,lmnproj_typat  (itypat))
   lmnmaxso=max(lmnmaxso,lmnprojso_typat(itypat))
   lnmax   =max(lnmax   ,lnproj_typat   (itypat))
   lnmaxso =max(lnmaxso ,lnprojso_typat (itypat))
 end do
!The initial coding (below) was not totally portable (MT 110215)
!lmnmax=max(maxval(lmnproj_typat(1:ntypat)),1)
!lmnmaxso=max(maxval(lmnprojso_typat(1:ntypat)),1)
!lnmax=max(maxval(lnproj_typat(1:ntypat)),1)
!lnmaxso=max(maxval(lnprojso_typat(1:ntypat)),1)

 if(maxval(lmnproj_typat(1:ntypat))==0)then
   write(msg, '(3a)' )&
    'Despite there is only a local part to pseudopotential(s),',ch10,&
    'lmnmax and lnmax are set to 1.'
   ABI_COMMENT(msg)
 end if

!XG040806 : These lines make modifications of lnmax and lmnmax
!that are unjustified in many cases, according to the many tests cases
!where they produce a changes, while the test case was working properly.
!One should understand better the needs, and code more appropriate changes ...
!lnmax/lmnmax has to be bigger than 1+lmax (for compatibility reasons)
!llmax=maxval(pspheads(1:ntypat)%lmax)+1 ! And this line might have trouble with HP compiler
!if (lnmax   <llmax) lnmax=llmax
!if (lnmaxso <llmax) lnmaxso=llmax
!if (lmnmax  <llmax) lmnmax=llmax
!if (lmnmaxso<llmax) lmnmaxso=llmax

 write(msg, '(a,a,i4,a,i4,3a,i4,a,i4,a)' ) ch10,&
& ' getdim_nloc : deduce lmnmax  =',lmnmax,', lnmax  =',lnmax,',',ch10,&
& '                      lmnmaxso=',lmnmaxso,', lnmaxso=',lnmaxso,'.'
 call wrtout(std_out,msg,'COLL')

 ABI_FREE(lmnproj_typat)
 ABI_FREE(lmnprojso_typat)
 ABI_FREE(lnproj_typat)
 ABI_FREE(lnprojso_typat)
 ABI_FREE(nproj_typat)
 ABI_FREE(nprojso_typat)

end subroutine getdim_nloc
!!***

!!****f* m_memeval/setmqgrid
!! NAME
!!  setmqgrid
!!
!! FUNCTION
!!  Sets the number of points needed to represent the pseudopotentials in
!!  reciprocal space for a specified resolution.
!!
!! INPUTS
!!  ecut=cutoff energy for the wavefunctions
!!  ecutdg=cutoff energy for the fine grid in case usepaw==1
!!  gprimd=primitive translation vectors for reciprocal space
!!  nptsgvec=number of points along the smallest primitive translation vector
!!    of the reciprocal space
!!  usepaw=1 if PAW is used, 0 otherwise
!!
!! OUTPUT
!!
!! PARENTS
!!      m_memeval,m_psps
!!
!! CHILDREN
!!      atomic_info,createwavefunctionsdescriptors,deallocate_lr
!!      memoryestimator,mkradim,wrtout,wvl_descr_atoms_set,wvl_descr_free
!!      wvl_setboxgeometry,xred2xcart
!!
!! SOURCE

subroutine setmqgrid(mqgrid,mqgriddg,ecut,ecutdg,gprimd,nptsgvec,usepaw)

!Arguments ------------------------------------
 integer , intent(inout)  :: mqgrid,mqgriddg
 integer , intent(in)  :: nptsgvec,usepaw
 real(dp), intent(in) :: ecut,ecutdg
 real(dp), intent(in) :: gprimd(3,3)

!Local variables-------------------------------
 integer :: mqgrid2,mqgriddg2
 real(dp) :: gmax,gmaxdg,gvecnorm
 character(len=500) :: msg

! *************************************************************************

 gvecnorm=sqrt(min(dot_product(gprimd(:,1),gprimd(:,1)), &
& dot_product(gprimd(:,2),gprimd(:,2)), &
& dot_product(gprimd(:,3),gprimd(:,3))))
 gmax=one/(sqrt2*pi)*sqrt(ecut)

 if (mqgrid == 0) then
   mqgrid2=ceiling(gmax/gvecnorm*nptsgvec)
   mqgrid=max(mqgrid2,3001)
   write(msg, '(5a,i0,a)' )&
&   'The number of points "mqgrid" in reciprocal space used for the',ch10,&
&   'description of the pseudopotentials has been set automatically',ch10,&
&   'by abinit to: ',mqgrid,'.'
   !ABI_COMMENT(msg)
 else
   mqgrid2=ceiling(gmax/gvecnorm*nptsgvec)
   if (mqgrid2>mqgrid) then
     write(msg, '(3a,i8,3a,i8,3a)' )&
&     'The number of points "mqgrid" in reciprocal space used for the',ch10,&
&     'description of the pseudopotentials is : ',mqgrid,'.',ch10,&
&     'It would be better to increase it to at least ',mqgrid2,', or',ch10,&
&     'let abinit choose it automatically by setting mqgrid = 0.'
     ABI_WARNING(msg)
   end if
 end if

 if (usepaw==1) then
   if(ecutdg<tol6)then
     write(msg,'(a)')'The value of (paw)ecutdg is zero or negative, which is forbidden.'
     ABI_ERROR(msg)
   end if
   gmaxdg=one/(sqrt2*pi)*sqrt(ecutdg)
   if (mqgriddg == 0) then
     mqgriddg2=ceiling(gmaxdg/gvecnorm*nptsgvec)
     mqgriddg=max(mqgriddg2,3001)
     write(msg, '(5a,i0,a)' )&
&     'The number of points "mqgriddg" in reciprocal space used for the',ch10,&
&     'description of the pseudopotentials has been set automatically',ch10,&
&     'by abinit to: ',mqgriddg,'.'
     !ABI_COMMENT(msg)
   else
     mqgriddg2=ceiling(gmax/gvecnorm*nptsgvec)
     if (mqgriddg2>mqgriddg) then
       write(msg, '(3a,i8,3a,i8,3a)' )&
&       'The number of points "mqgriddg" in reciprocal space used for the',ch10,&
&       'description of the pseudopotentials (fine grid) is :',mqgriddg,'.',ch10,&
&       'It would be better to increase it to at least ',mqgriddg2,', or',ch10,&
&       'let abinit choose it automatically by setting mqgrid = 0.'
       ABI_WARNING(msg)
     end if
   end if
 end if

end subroutine setmqgrid
!!***

!!****f* m_memeval/wvl_memory
!! NAME
!! wvl_memory
!!
!! FUNCTION
!! Estimation of the memory needed for waelet based computation job.
!! According to the value of the option variable,
!! might also try to allocate this amount of memory, and if it fails,
!! might estimate the available memory.
!!
!! INPUTS
!!  dtset=<type datafiles_type>contains all input variables.
!!  idtset=number of the current dataset
!!  mpi_enreg=information about MPI parallelization
!!  npsp=number of pseudopotentials
!!  option: if 0, no test of available memory
!!          if 1, the routine tries to allocate the estimated memory, for testing
!!                purposes, and if a failure occurs, the routine stops.
!!          if 2, like 1, but before stopping, the routine will provide
!!                an estimation of the available memory.
!!  pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! The estimator is the one provided by BigDFT.
!!
!! PARENTS
!!      m_memeval
!!
!! CHILDREN
!!      atomic_info,createwavefunctionsdescriptors,deallocate_lr
!!      memoryestimator,mkradim,wrtout,wvl_descr_atoms_set,wvl_descr_free
!!      wvl_setboxgeometry,xred2xcart
!!
!! SOURCE

subroutine wvl_memory(dtset, idtset, mpi_enreg, npsp, option, pspheads)

 use defs_wvltypes
 use m_abi2big, only : wvl_setBoxGeometry
 use m_wvl_descr_psp,    only : wvl_descr_free, wvl_descr_atoms_set

#if defined HAVE_BIGDFT
 use BigDFT_API, only: MemoryEstimator, createWavefunctionsDescriptors, deallocate_lr, &
      & atomic_info, memory_estimation
#endif

!Arguments ------------------------------------
  !scalars
  integer,intent(in) :: idtset, npsp, option
  type(dataset_type),intent(in) :: dtset
  type(MPI_type),intent(in) :: mpi_enreg
  !arrays
  type(pspheader_type),intent(in) :: pspheads(npsp)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
  !scalars
  integer :: ityp, i, mu, nstates, me, nproc, comm
  character(len=500) :: msg
  real(dp) :: ehomo, radfine
  type(wvl_internal_type) :: wvl
  type(memory_estimation) :: peakmem
  !arrays
  real(dp) :: acell(3), rprimd(3,3), rprim(3,3)
  real(dp), allocatable :: radii_cf(:,:)
  real(dp), allocatable :: xred(:,:), xcart(:,:)
#endif

! **************************************************************************

#if defined HAVE_BIGDFT

 comm=mpi_enreg%comm_wvl
 me=xmpi_comm_rank(comm)
 nproc=xmpi_comm_size(comm)

 if(option<0 .or. option>2)then
   write(msg, '(A,A,A,A,I0,A)') ch10,&
&   ' wvl_memory : BUG -',ch10,&
&   '  option=',option,' while the only allowed values are 0, 1, or 2.'
   call wrtout(std_out,msg,'COLL')
 end if

 wvl%paw%usepaw=0 !no PAW here
 nullify(wvl%rholoc%d)
 nullify(wvl%rholoc%msz)
 nullify(wvl%rholoc%rad)
 nullify(wvl%rholoc%radius)
 nullify(wvl%paw%spsi)
 nullify(wvl%paw%indlmn)
 nullify(wvl%paw%spsi)
 nullify(wvl%paw%indlmn)

 write(msg,*)' wvl_memory : analysis of memory needs '
 call wrtout(std_out,msg,'COLL')

 if(idtset>=100)then
   write(msg,'(80a,a,a,i5,a)')('=',mu=1,80),ch10,&
&   ' Values of the parameters that define the memory need for DATASET', idtset,&
&   ' (WVL).'
 else if(idtset/=0)then
   write(msg,'(80a,a,a,i3,a)')('=',mu=1,80),ch10,&
&   ' Values of the parameters that define the memory need for DATASET', idtset,&
&   ' (WVL).'
 else
   write(msg,'(80a,a,a,a)')('=',mu=1,80),ch10,&
&   ' Values of the parameters that define the memory need of the present run',&
&   ' (WVL).'
 end if
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 write(msg,'( a,f7.3,a,i7,2(a,F7.3),a,a,f7.3,a,i7 )' ) &
& '  wvl_hgrid =', dtset%wvl_hgrid , '   nwfshist =', dtset%nwfshist, &
& ' wvl_crmult =', dtset%wvl_crmult, ' wvl_frmult =', dtset%wvl_frmult, ch10,&
& '  tl_radius =', dtset%tl_radius , '  tl_nprccg =', dtset%tl_nprccg
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 if (dtset%nsppol == 2) then
   nstates = dtset%nelect
 else
   nstates = dtset%mband
 end if
 write(msg,'(4(a,i7))')&
& '      natom =', dtset%natom, '     ntypat =', dtset%ntypat, &
& '    nstates =', nstates,     '     nsppol =', dtset%nsppol
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 write(msg,'(80a)') ('=',mu=1,80)
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

!First, use eleconf to get radii_cf().
 ABI_MALLOC(radii_cf,(npsp, 3))
 do ityp = 1, npsp, 1
   call atomic_info(int(pspheads(ityp)%znuclpsp), int(pspheads(ityp)%zionpsp), ehomo = ehomo)

!  new method for assigning the radii
   radii_cf(ityp, 1) = one / sqrt(abs(two * ehomo))
   radfine = 100.d0
   do i = 0, 4, 1
     if (pspheads(ityp)%GTHradii(i) /= zero) then
       radfine = min(radfine, pspheads(ityp)%GTHradii(i))
     end if
   end do
   radii_cf(ityp,2) = radfine
 end do

!Compute the shifted positions and acell
 acell = dtset%acell_orig(1:3,1)
 call wvl_descr_atoms_set(acell, dtset%icoulomb, dtset%natom, dtset%ntypat, dtset%typat, wvl)
 ABI_MALLOC(xred,(3, dtset%natom))
 xred = dtset%xred_orig(:,:,1)
 rprimd = dtset%rprimd_orig(1:3,1:3,1)
 wvl%h(:) = dtset%wvl_hgrid
 call wvl_setBoxGeometry(1, radii_cf, rprimd, xred, &
& wvl, dtset%wvl_crmult, dtset%wvl_frmult)
!Compute acell and rprim from rprimd
 call mkradim(acell,rprim,rprimd)
 ABI_MALLOC(xcart,(3, dtset%natom))
 call xred2xcart(dtset%natom, rprimd, xcart, xred)
 call createWavefunctionsDescriptors(me, wvl%h(1), wvl%h(2), wvl%h(3), &
& wvl%atoms, xcart, radii_cf, dtset%wvl_crmult, dtset%wvl_frmult, wvl%Glr)
 call MemoryEstimator(nproc, dtset%nwfshist, wvl%Glr, &
& dtset%mband, dtset%nspinor, dtset%nkpt, 0, dtset%nsppol, &
& 0, dtset%iscf, peakmem)

 call deallocate_lr(wvl%Glr)
 call wvl_descr_free(wvl)
 ABI_FREE(radii_cf)
 ABI_FREE(xred)
 ABI_FREE(xcart)

 write(msg,'(80a,a)') ('=',mu=1,80), ch10
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) idtset,npsp,option,dtset%nstep,mpi_enreg%nproc,pspheads(1)%zionpsp
#endif

end subroutine wvl_memory
!!***

end module m_memeval
!!***
