!{\src2tex{textfont=tt}}
!!****f* ABINIT/memory
!! NAME
!! memory
!!
!! FUNCTION
!! Estimation of the memory needed for a ground-state job.
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
!!  mpi_enreg=informations about MPI parallelization
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
!!                 see ~abinit/doc/input_variables/vargs.htm#ngfft
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
!! in rhohxc and daughter routines (at maximum 22*1000 dp numbers),
!! as well as other arrays like
!! character(len=500) :: message (present in about 100 routines), or the different
!! arrays allocated in move.f, brdmin.f, gstate.f (xf array) or pspini.f
!! In the case 3<=occopt<=8 this amount is increased by 760 Kbytes
!! to take into account the arrays smdfun, occfun, entfun, workfun and xgrid,
!! declared in getnel
!!
!! The current version takes into account
!! 1) and 2) the "main chain" in its two slightly different versions :
!! driver - gstate - (move or brdmin) - scfcv - vtorho - vtowfk -
! !!     cgwf - getghc - fourwf or (nonlop+opernl)
!! 3) the xc chain :
!! driver - gstate - (move or brdmin) - scfcv - (vresfo) - rhohxc - xcden
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
!!     goes through energy instead of vtorho).
!!
!! Also, it is assumed that the potentials are non-local, even if there
!!     are local ! It would be necessary to update this routine
!!     now that the beginning of psp files is read before
!!     the present call (XG 980502)
!!
!! One might also estimate if there must be a chain arriving at :
!!  strnps , mkffnl, mkcore, mklocl, mkrho, prcpot, irrzg, initro,
!!  clnup1.
!! This is because there are allocated arrays in these routines.
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

subroutine memory(n1xccc,extrapwf,getcell,idtset,icoulomb,intxc,ionmov,iout,densfor_pred,iprcel,&
& iscf,jdtset,lmnmax,lnmax,&
& mband,mffmem,mgfft,mgfftdiel,mgfftf,mkmem,mpi_enreg,mpsang,mpssoang,mpw,mqgrid_ff,mqgrid_vl,&
& natom,nband,nfft,nfftdiel,nfftf,ngfft,ngfftdiel,ngfftf,nimage,&
& nkpt,nloalg,npsp,npulayit,npwdiel,nspden,nspinor,nsppol,nsym,ntypat,&
& occopt,optforces,option,optstress,pawcpxocc,pawmixdg,pawnhatxc,pawspnorb,pawstgylm,&
& prtvol,pspheads,tfkinfunc,typat,ucvol,usepaw,useylm,use_gpu_cuda,xclevel)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memory'
 use interfaces_14_hidewrite
 use interfaces_57_iovars, except_this_one => memory
!End of the abilint section

 implicit none

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
 integer :: ngrad,nprocwf,nspgrad,rhoij_nspden
 real(dp) :: mbcg,mbdiskpd,mbdiskwf,mbf_fftgr,mbgylm
 character(len=500) :: message
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
   write(message,'(A,I0,A)')'option=',option,' while the only allowed values are 0, 1, or 2.'
   MSG_BUG(message)
 end if

!firstchar=' ';if (use_gpu_cuda==1) firstchar='_'
 cmpw(:)=zero ; cfft(:)=zero ; cfftf(:)=zero ; cadd(:)=zero  
 dttyp(:)=0

 my_natom=natom;if (mpi_enreg%nproc_atom>1) my_natom=mpi_enreg%my_natom

 call wrtout(std_out,'memory : analysis of memory needs ','COLL')

 if(jdtset>=100)then
   write(message,'(80a,a,a,i5,a)')('=',mu=1,80),ch10,&
&   ' Values of the parameters that define the memory need for DATASET',jdtset,'.'
 else if(jdtset/=0)then
   write(message,'(80a,a,a,i3,a)')('=',mu=1,80),ch10,&
&   ' Values of the parameters that define the memory need for DATASET',jdtset,'.'
 else
   write(message,'(80a,a,a)')('=',mu=1,80),ch10,&
&   ' Values of the parameters that define the memory need of the present run '
 end if
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 write(message,'( 4(a,i8),a,4(a,i8) )' ) &
& '     intxc =',intxc   ,'    ionmov =',ionmov,&
& '      iscf =',iscf    ,'    lmnmax =',lmnmax,ch10,&
& '     lnmax =',lnmax   ,'     mgfft =',mgfft,&
& '  mpssoang =',mpssoang,'    mqgrid =',mqgrid_vl
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 write(message,'( 4(a,i8),a,4(a,i8),a,4(a,i8) )' ) &
& '     natom =',natom  ,'  nloc_mem =',nloalg(2)*(nloalg(3)+1),&
& '    nspden =',nspden ,'   nspinor =',nspinor,ch10,&
& '    nsppol =',nsppol ,'      nsym =',nsym,&
& '    n1xccc =',n1xccc ,'    ntypat =',ntypat,ch10,&
& '    occopt =',occopt ,'   xclevel =',xclevel
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 write(message,'(4(3(a,i12),a))') &
& '-    mband =',mband  ,'        mffmem =',mffmem,&
& '         mkmem =',mkmem  ,ch10,&
& '       mpw =',mpw    ,'          nfft =',nfft ,&
& '          nkpt =',nkpt   
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 if (my_natom/=natom)then
   write(message,'(a,i10)') 'Pmy_natom=',my_natom
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

!Additional information if imgmov is activated (use of replicas of the cell)
 if (nimage>1) then
   write(message,'(1(a,i10))' ) '  nimage =',nimage
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

!Additional information on FFT grids if PAW
 if (usepaw==1) then
   write(message, '(a,a,a,i10,a,i10)' )&
&   ' PAW method is used; the additional fine FFT grid is defined by:',ch10,&
&   '   mgfftf=',mgfftf,'    nfftf =',nfftf
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

!Additional information if GPU
 if (use_gpu_cuda==1) then
!  write(message, '(a)' )' GPU method is used'
!  call wrtout(iout,message,'COLL')
!  call wrtout(std_out,message,'COLL')
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

   write(message, '(a,a,a,i10,a,i6,a,i10,a,i10)' )&
&   ' For the susceptibility and dielectric matrices, or tddft :',ch10,&
&   '   mgfft =',mgfftdiel,'  nbnd_in_blk=',nbnd_in_blk,'    nfft =',nfftdiel,&
&   '     npw =',npwdiel
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   ndiel4=ngfftdiel(4) ; ndiel5=ngfftdiel(5) ; ndiel6=ngfftdiel(6)
   ndiel456=ndiel4*ndiel5*ndiel6
 else
!  To be sure of initialisation.
   ndiel456 = 1
 end if

 write(message,'(80a)') ('=',mu=1,80)
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 if(getcell>0 .or. (getcell<0 .and. idtset+getcell>0) )then
   write(message,'(a,a,a,a,a,a,i3,a,i3,a,a,a,a,a,a)' )ch10,&
&   ' memory : COMMENT -',ch10,&
&   '  The determination of memory needs at this stage is meaningless,',ch10,&
&   '  since getcell = ',getcell,' is non-zero, while idtset=',idtset,'.',ch10,&
&   '  The following numbers are obtained by supposing that acell and rprim',ch10,&
&   '  are NOT taken from a previous dataset. You cannot rely on them.',ch10
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

!Compute number of atoms per type for current proc
 nattyp(:)=0
 do ii=1,natom
   nattyp(typat(ii))=nattyp(typat(ii))+1
 end do

!PAW: store useful dims
 if (usepaw==1) then
   ABI_ALLOCATE(basis_size,(npsp))
   ABI_ALLOCATE(l_size,(npsp))
   ABI_ALLOCATE(lmn_size,(npsp))
   ABI_ALLOCATE(lmn2_size,(npsp))
   ABI_ALLOCATE(mesh_size,(npsp))
   ABI_ALLOCATE(shape_type,(npsp))
   ABI_ALLOCATE(pawver,(npsp))
   ABI_ALLOCATE(rshp,(npsp))
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
   ABI_ALLOCATE(my_nattyp,(ntypat))
   if ((mpi_enreg%nproc_atom<=1).or.(.not.associated(mpi_enreg%my_atmtab))) then
     my_nattyp=nattyp
   else
     my_nattyp=0
     do ii=1,my_natom
       jj=typat(mpi_enreg%my_atmtab(ii))
       my_nattyp(jj)=my_nattyp(jj)+1
     end do
   end if
 else
!  Do the allocation to avoid uninitialised variables.
   ABI_ALLOCATE(my_nattyp,(1))
   ABI_ALLOCATE(basis_size,(1))
   ABI_ALLOCATE(l_size,(1))
   ABI_ALLOCATE(lmn_size,(1))
   ABI_ALLOCATE(lmn2_size,(1))
   ABI_ALLOCATE(mesh_size,(1))
   ABI_ALLOCATE(shape_type,(1))
   ABI_ALLOCATE(pawver,(1))
   ABI_ALLOCATE(rshp,(1))
   rhoij_nspden=nspden
   l_size_max=1
   l_max=1
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
       cadd(19)=cadd(19)+histsz*2*my_nattyp(ii)*lmn2_size(ii)*rhoij_nspden*pawcpxocc  ! %pawrhoij()%rhoijp
       cadd(20)=cadd(20)+histsz*2*my_nattyp(ii)*(2+lmn2_size(ii))*nspden              ! %pawrhoij()%rhoijselect
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

!(3)                     in rhohxc, xcden -------------------------------

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
!  In this case, rhohxc is called from rhotov also,
!  for which vresid was allocated in scfcv
!  vresid
   cfftf(35)=nspden               ; dttyp(35)=8
 end if
!Poisson's solver with zero padding
 if (icoulomb == 1) then
   cfft(36) = 8                   ; dttyp(36) = 8
   cadd(36) = ngfft(4) * ngfft(5) * ngfft(6) - nfft
 end if

!Note : in hartre, called by rhohxc, one uses
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

!eknlk, enlnk, grnlnk
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
!eig_k, ek_k, enl_k, grnl_k, occ_k, resid_k
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

!(6f)  in pawmkrhoij or symrhoij called by pawmkrho, called by vtorho--------
!only when paralellim over atoms is activated
 dttyp(63)=8
 if((usepaw==1) .and. ((iscf>0) .or. (iscf == -3) .and. mpi_enreg%nproc_atom>1 ))then
   do ii=1,ntypat
     cadd(63)=cadd(63)+nattyp(ii)*lmn2_size(ii)*rhoij_nspden*pawcpxocc   ! Rhoij_gather and related data
     cadd(63)=cadd(63)+nattyp(ii)*(2+lmn2_size(ii))    ! Rhoij_gather (rhoijselect, ...)
   end do
 end if

!(7)                     in vtowfk----------------------------------------

!evec
 cadd(71)=2*mband*mband        ; dttyp(71)=8
!subham, subvnl(if not PAW)
 cadd(72)=(1+usepaw)*mband*(mband+1)    ; dttyp(72)=8
!gkpsq
 cmpw(73)=1                    ; dttyp(73)=8
!ffnl
 cmpw(74)=2*ntypat*lmnmax      ; dttyp(74)=8
!ph3d
 matblk=NLO_MINCAT
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

!conjgr, cwavef, direc, gh_direc, gvnl_direc
 cmpw(81)=2*5*nspinor          ; dttyp(81)=8
!ghc,gvnlc
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
 ABI_DEALLOCATE(my_nattyp)
 ABI_DEALLOCATE(basis_size)
 ABI_DEALLOCATE(l_size)
 ABI_DEALLOCATE(lmn_size)
 ABI_DEALLOCATE(lmn2_size)
 ABI_DEALLOCATE(mesh_size)
 ABI_DEALLOCATE(pawver)
 ABI_DEALLOCATE(shape_type)
 ABI_DEALLOCATE(rshp)

!---------------------------------------------------------------------
!Now, analyze the data

 call memana(cadd,cfft,cfftf,chain,cmpw,dttyp,iout,iprcel,iscf,&
& marrays,mbcg,mbdiskpd,mbdiskwf,mbf_fftgr,mbgylm,mffmem,&
& mpw,natom,nchain,nfft,nfftf,occopt,option,prtvol)

end subroutine memory
!!***
