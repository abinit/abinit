!!****m* ABINIT/m_lobpcgwf_old
!! NAME
!!   m_lobpcgwf_old
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group ()
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

module m_lobpcgwf_old

 implicit none

 private
!!***

 public :: lobpcgwf
!!***

contains
!!***

!!****f* ABINIT/lobpcgwf
!! NAME
!! lobpcgwf
!!
!! FUNCTION
!! this routine updates the whole wave functions at a given k-point,
!! using the lobpcg method
!! for a given spin-polarization, from a fixed hamiltonian
!! but might also simply compute eigenvectors and eigenvalues at this k point.
!! it will also update the matrix elements of the hamiltonian.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group (FBottin,GZ,AR,MT,FDahm)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variales for this dataset
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cg
!!  igsc=shift to be applied on the location of data in the array gsc
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (hartree)
!!  mcg=second dimension of the cg array
!!  mgsc=second dimension of the gsc array
!!  mpi_enreg=informations about MPI parallelization
!!  nband_k=number of bands at this k point for that spin polarization
!!  nbdblock : number of blocks
!!  npw_k=number of plane waves at this k point
!!  prtvol=control print volume and debugging output
!!  use_totvnlx=1 if one has to compute totvnlx
!!
!! OUTPUT
!!  resid_k(nband_k)=residuals for each states
!!  subham(nband_k*(nband_k+1))=the matrix elements of h
!!  If gs_hamk%usepaw==0:
!!    gsc(2,mgsc)=<g|s|c> matrix elements (s=overlap)
!!    totvnlx(nband_k*use_totvnlx,nband_k*use_totvnlx)=the matrix elements of vnl+vfockACE
!!
!! SIDE EFFECTS
!!  cg(2,mcg)=updated wavefunctions
!!
!! PARENTS
!!      vtowfk
!!
!! CHILDREN
!!      abi_xcopy,abi_xgemm,abi_xheev,abi_xhegv,abi_xorthonormalize,abi_xtrsm
!!      alloc_on_gpu,copy_from_gpu,copy_on_gpu,dealloc_on_gpu,getghc,gpu_xgemm
!!      gpu_xorthonormalize,gpu_xtrsm,prep_getghc,setwfparameter,timab,wfcopy
!!      wrtout,xmpi_sum,xprecon
!!
!! SOURCE

subroutine lobpcgwf(cg,dtset,gs_hamk,gsc,icg,igsc,kinpw,mcg,mgsc,mpi_enreg,&
&                   nband_k,nbdblock,npw_k,prtvol,resid_k,subham,totvnlx,use_totvnlx)


 use defs_basis
 use m_abicore
 use m_lobpcg
 use m_abi_linalg
 use m_wfutils
 use m_xmpi
 use m_errors
 use iso_c_binding
 use m_dtset

 use defs_abitypes, only : mpi_type
 use m_time,        only : timab
 use m_hamiltonian, only : gs_hamiltonian_type
 use m_pawcprj,     only : pawcprj_type
 use m_getghc,      only : getghc
 use m_prep_kgb,    only : prep_getghc

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: icg,igsc,mcg,mgsc,nband_k,nbdblock,npw_k,prtvol,use_totvnlx
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(dataset_type),intent(in) :: dtset
 type(mpi_type),intent(inout) :: mpi_enreg
 real(dp),intent(inout) :: cg(2,mcg),gsc(2,mgsc)
 real(dp),intent(in) :: kinpw(npw_k)
 real(dp),intent(out) :: resid_k(nband_k)
 real(dp),intent(inout) :: subham(nband_k*(nband_k+1))
 real(dp),intent(inout) :: totvnlx((3-gs_hamk%istwf_k)*nband_k*use_totvnlx,nband_k*use_totvnlx)

!Local variables-------------------------------
 integer, parameter :: tim_getghc=5
 integer :: activepsize,activersize,bblocksize,bigorder,blocksize,cpopt
 integer :: cond_try
 integer :: iblocksize,iblock,ierr,ii,info,istwf_k,isubh
 integer :: iterationnumber
 integer :: iwavef,i1,i2,i3,i4,maxiterations,my_nspinor
 integer :: nrestart,optekin,optpcon,restart
 integer :: sij_opt,timopt,tim_wfcopy,tim_xeigen
 integer :: tim_xortho,tim_xprecon,use_lapack_gpu,use_linalg_gpu,vectsize
 logical :: gen_eigenpb
 integer :: cplx
 real(dp) :: condestgramb,deltae,deold,dum
 complex(dpc) :: cminusone
 real(dp) :: zvar(2)
 logical :: havetoprecon
 real(dp) :: tsec(2)
 real(dp), allocatable :: gwavef(:,:),cwavef(:,:),gvnlxc(:,:)
 real(dp), allocatable :: swavef(:,:)
 real(dp), allocatable :: residualnorms(:),eigen(:)
 real(dp), allocatable :: tmpeigen(:)
 real(dp), allocatable :: pcon(:,:)
 real(dp), allocatable :: blockvectorx(:,:),blockvectorvx(:,:),blockvectorax(:,:),blockvectorbx(:,:)
 real(dp),allocatable  :: blockvectorr(:,:),blockvectorvr(:,:),blockvectorar(:,:),blockvectorbr(:,:)
 real(dp),allocatable  :: blockvectorp(:,:),blockvectorvp(:,:),blockvectorap(:,:),blockvectorbp(:,:),blockvectordumm(:,:)
 real(dp),allocatable  :: blockvectory(:,:),blockvectorby(:,:),blockvectorz(:,:)
 real(dp),allocatable  :: gramxax(:,:),gramxar(:,:),gramxap(:,:),gramrar(:,:),gramrap(:,:),grampap(:,:)
 real(dp),allocatable  :: gramxbx(:,:),gramxbr(:,:),gramxbp(:,:),gramrbr(:,:),gramrbp(:,:),grampbp(:,:)
 real(dp),allocatable  :: coordx1(:,:),coordx2(:,:),coordx3(:,:),lambda(:,:),grama(:,:),gramb(:,:),gramyx(:,:)
 real(dp),allocatable  :: tmpgramb(:,:),transf3(:,:,:),transf5(:,:,:)
 real(dp), allocatable :: tsubham(:,:)
 type(pawcprj_type) :: cprj_dum(gs_hamk%natom,0)
 character(len=500) :: message
 character, dimension(2) :: cparam
 type(c_ptr) :: A_gpu,C_gpu,coordx2_gpu,coordx3_gpu,bblockvector_gpu,gram_gpu
 type(c_ptr) :: blockvectorr_gpu,blockvectorar_gpu,blockvectorbr_gpu

!Index of a given band
!gramindex(iblocksize)=(iblocksize-1)*cplx+1

! *********************************************************************

 DBG_ENTER("COLL")

 call timab(530,1,tsec)
 if(abs(dtset%timopt)==4) then
   call timab(520,1,tsec)
 end if

!###########################################################################
!################ INITIALISATION  ##########################################
!###########################################################################

!For timing
 timopt=dtset%timopt
 tim_wfcopy=584
 !tim_xcopy=584
 tim_xeigen=587
 !tim_xgemm=532
 tim_xortho=535
 tim_xprecon=536
 !tim_xtrsm=535

!Variables
 maxiterations=dtset%nline
 gen_eigenpb=(gs_hamk%usepaw==1)
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 cminusone=-cone
 istwf_k=gs_hamk%istwf_k
 info = 0
 cparam(1)='t'
 cparam(2)='c'

!Depends on istwfk
 if ( istwf_k == 2 ) then
   cplx=1
   if (mpi_enreg%me_g0 == 1) then
     vectsize=2*npw_k*my_nspinor-1
   else
     vectsize=2*npw_k*my_nspinor
   end if
 else
   cplx=2
   vectsize=npw_k*my_nspinor
 end if

!For preconditionning
 optekin=0;if (dtset%wfoptalg>10) optekin=0
 optpcon=1;if (dtset%wfoptalg>10) optpcon=0

!For communication
 blocksize=mpi_enreg%nproc_fft
 if(mpi_enreg%paral_kgb==1) blocksize=mpi_enreg%nproc_band*mpi_enreg%bandpp
 !IF you want to compare with new lobpcg in sequential uncomment the following
 !line
 !blocksize=mpi_enreg%nproc_band*mpi_enreg%bandpp

!Iniitializations/allocations of GPU parallelism
 use_linalg_gpu=0;use_lapack_gpu=0
 if ((dtset%use_gpu_cuda==1).and. &
& (vectsize*blocksize*blocksize>dtset%gpu_linalg_limit)) use_linalg_gpu=1
#if defined HAVE_LINALG_MAGMA
 use_lapack_gpu=use_linalg_gpu
#endif
 if(use_linalg_gpu==1) then
   call alloc_on_gpu(A_gpu,cplx*dp*vectsize*blocksize)
   call alloc_on_gpu(C_gpu,cplx*dp*vectsize*blocksize)
   call alloc_on_gpu(blockvectorr_gpu,cplx*dp*vectsize*blocksize)
   call alloc_on_gpu(blockvectorar_gpu,cplx*dp*vectsize*blocksize)
   call alloc_on_gpu(blockvectorbr_gpu,cplx*dp*vectsize*blocksize)
   call alloc_on_gpu(coordx2_gpu,cplx*dp*blocksize*blocksize)
   call alloc_on_gpu(coordx3_gpu,cplx*dp*blocksize*blocksize)
 end if

 if(abs(dtset%timopt)==4) then
   call timab(520,2,tsec)
 end if

!###########################################################################
!################ BIG LOOP OVER BLOCKS  ####################################
!###########################################################################

 do iblock=1,nbdblock

   if(abs(dtset%timopt)==4) then
     call timab(521,1,tsec)
   end if

   havetoprecon=.true.
   nrestart=0
   bblocksize=(iblock-1)*blocksize

!  allocations
   ABI_ALLOCATE(pcon,(npw_k,blocksize))
   ABI_ALLOCATE(blockvectorx,(cplx*vectsize,blocksize))
   ABI_ALLOCATE(blockvectorax,(cplx*vectsize,blocksize))
   ABI_ALLOCATE(blockvectorbx,(cplx*vectsize,blocksize))
   ABI_ALLOCATE(blockvectorr,(cplx*vectsize,blocksize))
   ABI_ALLOCATE(blockvectorar,(cplx*vectsize,blocksize))
   ABI_ALLOCATE(blockvectorbr,(cplx*vectsize,blocksize))
   ABI_ALLOCATE(blockvectorp,(cplx*vectsize,blocksize))
   ABI_ALLOCATE(blockvectorap,(cplx*vectsize,blocksize))
   ABI_ALLOCATE(blockvectorbp,(cplx*vectsize,blocksize))
   ABI_ALLOCATE(blockvectordumm,(cplx*vectsize,blocksize))
   ABI_ALLOCATE(gramxax,(cplx*blocksize,blocksize))
   ABI_ALLOCATE(gramxar,(cplx*blocksize,blocksize))
   ABI_ALLOCATE(gramxap,(cplx*blocksize,blocksize))
   ABI_ALLOCATE(gramrar,(cplx*blocksize,blocksize))
   ABI_ALLOCATE(gramrap,(cplx*blocksize,blocksize))
   ABI_ALLOCATE(grampap,(cplx*blocksize,blocksize))
   ABI_ALLOCATE(gramxbx,(cplx*blocksize,blocksize))
   ABI_ALLOCATE(gramxbr,(cplx*blocksize,blocksize))
   ABI_ALLOCATE(gramxbp,(cplx*blocksize,blocksize))
   ABI_ALLOCATE(gramrbr,(cplx*blocksize,blocksize))
   ABI_ALLOCATE(gramrbp,(cplx*blocksize,blocksize))
   ABI_ALLOCATE(grampbp,(cplx*blocksize,blocksize))
   ABI_ALLOCATE(transf3,(cplx*blocksize,blocksize,3))
   ABI_ALLOCATE(transf5,(cplx*blocksize,blocksize,5))
   ABI_ALLOCATE(lambda,(cplx*blocksize,blocksize))
   ABI_ALLOCATE(residualnorms,(blocksize))

   ABI_ALLOCATE(blockvectory,(cplx*vectsize,bblocksize))
   ABI_ALLOCATE(blockvectorby,(cplx*vectsize,bblocksize))
   ABI_ALLOCATE(gramyx,(cplx*bblocksize,blocksize))
   if (gs_hamk%usepaw==0) then
     ABI_ALLOCATE(blockvectorvx,(cplx*vectsize,blocksize))
     ABI_ALLOCATE(blockvectorvr,(cplx*vectsize,blocksize))
     ABI_ALLOCATE(blockvectorvp,(cplx*vectsize,blocksize))
   end if

   if(use_linalg_gpu==1) then
     if(iblock/=1) then
       call alloc_on_gpu(bblockvector_gpu,cplx*dp*vectsize*bblocksize)
       call alloc_on_gpu(gram_gpu,cplx*dp*bblocksize*blocksize)
     else
       call alloc_on_gpu(bblockvector_gpu,cplx*dp*vectsize*blocksize)
       call alloc_on_gpu(gram_gpu,cplx*dp*blocksize*blocksize)
     end if
   end if

! Initialize global variables in m_wfutils.
   call setWFParameter(cplx,mpi_enreg%me_g0,npw_k,my_nspinor,icg,igsc,blocksize)

!  transfer array of wf coeff in iblock to blockvectorx
   call wfcopy('D',blocksize*vectsize,cg,1,blockvectorx,1,blocksize,iblock,'C',withbbloc=.true.,&
&   timopt=timopt,tim_wfcopy=tim_wfcopy)

!  !!!!!!!!!!!!!!!!!!!!!!!! Begin if iblock /=1 !!!!!!!!!!!!!!!!!!!!!!!!!!
!  transfer array of wf coeff less than iblock to blockvectory
   if(iblock /=1) then
     call wfcopy('D',bblocksize*vectsize,cg,1,blockvectory,1,bblocksize,iblock,'C',withbbloc=.false.,&
&     timopt=timopt,tim_wfcopy=tim_wfcopy)

     if(gen_eigenpb) then
       call wfcopy('D',bblocksize*vectsize,gsc,1,blockvectorby,1,bblocksize,iblock,'S',withbbloc=.false.,&
&       timopt=timopt,tim_wfcopy=tim_wfcopy)
     else
       call abi_xcopy(vectsize*bblocksize,blockvectory,1,blockvectorby,1,x_cplx=x_cplx)
     end if

!    b-orthogonalize x to the constraint y (supposed b-orthonormal)
!    blockvectorx=blockvectorx-matmul(blockvectory,matmul((blockvectorby)^T,blockvectorx))

     call abi_xgemm(cparam(cplx),'n',bblocksize,blocksize,vectsize,cone,blockvectorby,&
&     vectsize,blockvectorx,vectsize,czero,gramyx,bblocksize,x_cplx=x_cplx)

     if(abs(dtset%timopt)==3) then
       call timab(533,1,tsec)
     end if
     call xmpi_sum(gramyx,mpi_enreg%comm_bandspinorfft,ierr)
     if(abs(dtset%timopt)==3) then
       call timab(533,2,tsec)
     end if

     call abi_xgemm('n','n',vectsize,blocksize,bblocksize,cminusone,blockvectory,&
&     vectsize,gramyx,bblocksize,cone,blockvectorx,vectsize,x_cplx=x_cplx)

   end if
!  !!!!!!!!!!!!!!!!!!!!!!!! End if iblock /=1 !!!!!!!!!!!!!!!!!!!!!!!!!!!

   ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor*blocksize))
   ABI_ALLOCATE(gwavef,(2,npw_k*my_nspinor*blocksize))
   ABI_ALLOCATE(gvnlxc,(2,npw_k*my_nspinor*blocksize))
   ABI_ALLOCATE(swavef,(2,npw_k*my_nspinor*blocksize))

   call wfcopy('I',vectsize*blocksize,blockvectorx,1,cwavef,1,blocksize,iblock,'W',withbbloc=.false.,&
&   timopt=timopt,tim_wfcopy=tim_wfcopy)

   if(abs(dtset%timopt)==4) then
     call timab(521,2,tsec)
   end if
   if(abs(dtset%timopt)==4) then
     call timab(526,1,tsec)
   end if

   cpopt=-1;sij_opt=0;if (gen_eigenpb) sij_opt=1

   if (mpi_enreg%paral_kgb==0) then
     call getghc(cpopt,cwavef,cprj_dum,gwavef,swavef,gs_hamk,gvnlxc,dum,&
&     mpi_enreg,blocksize,prtvol,sij_opt,tim_getghc,0)
   else
     call prep_getghc(cwavef,gs_hamk,gvnlxc,gwavef,swavef,dum,blocksize,mpi_enreg,&
&     prtvol,sij_opt,cpopt,cprj_dum,already_transposed=.false.)
   end if
   if(abs(dtset%timopt)==4) then
     call timab(526,2,tsec)
   end if
   if(abs(dtset%timopt)==4) then
     call timab(522,1,tsec)
   end if

   if ( gen_eigenpb ) then
     call wfcopy('D',vectsize*blocksize,swavef,1,blockvectorbx,1,blocksize,iblock,'W',withbbloc=.false.,&
&     timopt=timopt,tim_wfcopy=tim_wfcopy)
   else
     call wfcopy('D',vectsize*blocksize,gvnlxc,1,blockvectorvx,1,blocksize,iblock,'W',withbbloc=.false.,&
&     timopt=timopt,tim_wfcopy=tim_wfcopy)
     call abi_xcopy(vectsize*blocksize,blockvectorx,1,blockvectorbx,1,x_cplx=x_cplx)
   end if

   call wfcopy('D',vectsize*blocksize,gwavef,1,blockvectorax,1,blocksize,iblock,'W',withbbloc=.false.,&
&   timopt=timopt,tim_wfcopy=tim_wfcopy)

   ABI_DEALLOCATE(cwavef)
   ABI_DEALLOCATE(gwavef)
   ABI_DEALLOCATE(gvnlxc)
   ABI_DEALLOCATE(swavef)

   call abi_xorthonormalize(blockvectorx,blockvectorbx,blocksize,mpi_enreg%comm_bandspinorfft,gramxbx,vectsize,&
&   x_cplx,timopt=timopt,tim_xortho=tim_xortho)

   call abi_xtrsm('r','u','n','n',vectsize,blocksize,cone,gramxbx,blocksize,blockvectorbx,vectsize,x_cplx=x_cplx)
   call abi_xtrsm('r','u','n','n',vectsize,blocksize,cone,gramxbx,blocksize,blockvectorax,vectsize,x_cplx=x_cplx)

   if (gs_hamk%usepaw==0) then
     call abi_xtrsm('r','u','n','n',vectsize,blocksize,cone,gramxbx,blocksize,blockvectorvx,vectsize,x_cplx=x_cplx)
   end if

!  Do rayleigh ritz on a in space x
!  gramxax=matmul(transpose(blockvectorx),blockvectorax)
   call abi_xgemm(cparam(cplx),'n',blocksize,blocksize,vectsize,cone,blockvectorx,&
&   vectsize,blockvectorax,vectsize,czero,gramxax,blocksize,x_cplx=x_cplx)

   if(abs(dtset%timopt)==3) then
     call timab(533,1,tsec)
   end if
   call xmpi_sum(gramxax,mpi_enreg%comm_bandspinorfft,ierr)
   if(abs(dtset%timopt)==3) then
     call timab(533,2,tsec)
   end if
   ABI_ALLOCATE(eigen,(blocksize))

   call abi_xheev('v','u',blocksize,gramxax,blocksize,eigen,x_cplx=cplx,istwf_k=istwf_k, &
   timopt=timopt,tim_xeigen=tim_xeigen,use_slk=dtset%use_slk,use_gpu=use_lapack_gpu)

!  blockvectorx=matmul(blockvectorx,gramxax)
   call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
&   vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize,x_cplx=x_cplx)
   call abi_xcopy(vectsize*blocksize,blockvectordumm,1,blockvectorx,1,x_cplx=x_cplx)

!  blockvectorax=matmul(blockvectorax,gramxax)
   call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorax,&
&   vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize,x_cplx=x_cplx)
   call abi_xcopy(vectsize*blocksize,blockvectordumm,1,blockvectorax,1,x_cplx=x_cplx)

!  blockvectorvx=matmul(blockvectorvx,gramxax)
   if (gs_hamk%usepaw==0) then
     call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorvx,&
&     vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize,x_cplx=x_cplx)
     call abi_xcopy(vectsize*blocksize,blockvectordumm,1,blockvectorvx,1,x_cplx=x_cplx)
   end if

!  blockvectorbx=matmul(blockvectorbx,gramxax)
   call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbx,&
&   vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize,x_cplx=x_cplx)
   call abi_xcopy(vectsize*blocksize,blockvectordumm,1,blockvectorbx,1,x_cplx=x_cplx)

   do iblocksize=1,blocksize
     zvar=(/eigen(iblocksize),zero/)
     call abi_xcopy(1,zvar,1,lambda(cplx*(iblocksize-1)+1:cplx*iblocksize,iblocksize),1,x_cplx=x_cplx)
   end do
   ABI_DEALLOCATE(eigen)

   if(abs(dtset%timopt)==4) then
     call timab(522,2,tsec)
   end if

!  ###########################################################################
!  ################ PERFORM LOOP ON NLINE ####################################
!  ###########################################################################
!  now the main alogrithm
   iter: do iterationnumber=1,maxiterations

     if(abs(dtset%timopt)==4) then
       call timab(523,1,tsec)
     end if

!    Build residual
!    blockvectorr=blockvectorax-matmul(blockvectorx,lambda)
     call xprecon(blockvectorbx,lambda,blocksize,&
&     iterationnumber,kinpw,mpi_enreg,npw_k,my_nspinor,&
&     optekin,optpcon,pcon,blockvectorax,blockvectorr,vectsize,timopt=timopt,tim_xprecon=tim_xprecon)

     residualnorms=sum(blockvectorr**2,dim=1)

     if(abs(dtset%timopt)==3) then
       call timab(533,1,tsec)
     end if
     call xmpi_sum(residualnorms,mpi_enreg%comm_bandspinorfft,ierr)
     if(abs(dtset%timopt)==3) then
       call timab(533,2,tsec)
     end if

     resid_k(bblocksize+1:bblocksize+blocksize)=residualnorms(1:blocksize)

!    If residual sufficiently small stop line minimizations
     if (abs(maxval(residualnorms(1:blocksize)))<dtset%tolwfr) then
       if (prtvol > 0) then
         write(message, '(a,i0,a,i0,a,es12.4)' ) &
&         ' lobpcgwf: block ',iblock,' converged after ',iterationnumber,&
&         ' line minimizations: maxval(resid(1:blocksize)) =',maxval(residualnorms(1:blocksize))
         call wrtout(std_out,message,'PERS')
       end if
       havetoprecon=.false.
       exit
     end if

     if(use_linalg_gpu==1) then
       call copy_on_gpu(blockvectorr,blockvectorr_gpu,cplx*dp*vectsize*blocksize)
     end if

     if(iblock /=1) then
!      Residuals orthogonal to blockvectorby
!      blockvectorr=blockvectorr-matmul(blockvectory,matmul((blockvectorby)^T,blockvectorr))

       if(use_linalg_gpu==1) then
         call copy_on_gpu(blockvectorby,bblockvector_gpu,cplx*dp*vectsize*bblocksize)
         call gpu_xgemm(cplx,cparam(cplx),'n',bblocksize,blocksize,vectsize,cone,bblockvector_gpu,&
&         vectsize,blockvectorr_gpu,vectsize,czero,gram_gpu,bblocksize)
         call copy_from_gpu(gramyx,gram_gpu,cplx*dp*bblocksize*blocksize)
       else
         call abi_xgemm(cparam(cplx),'n',bblocksize,blocksize,vectsize,cone,blockvectorby,&
&         vectsize,blockvectorr,vectsize,czero,gramyx,bblocksize,x_cplx=x_cplx)
       end if

       if(abs(dtset%timopt)==3) then
         call timab(533,1,tsec)
       end if
       call xmpi_sum(gramyx,mpi_enreg%comm_bandspinorfft,ierr)
       if(abs(dtset%timopt)==3) then
         call timab(533,2,tsec)
       end if

       if(use_linalg_gpu==1) then
         call copy_on_gpu(gramyx,gram_gpu,cplx*dp*bblocksize*blocksize)
         call copy_on_gpu(blockvectory,bblockvector_gpu,cplx*dp*vectsize*bblocksize)
         call gpu_xgemm(cplx,'n','n',vectsize,blocksize,bblocksize,cminusone,bblockvector_gpu,&
&         vectsize,gram_gpu,bblocksize,cone,blockvectorr_gpu,vectsize)
       else
         call abi_xgemm('n','n',vectsize,blocksize,bblocksize,cminusone,blockvectory,&
&         vectsize,gramyx,bblocksize,cone,blockvectorr,vectsize,x_cplx=x_cplx)
       end if

     end if

!    Residuals orthogonal to blockvectorx
!    blockvectorr=blockvectorr-matmul(blockvectorx,matmul((blockvectorbx)^T,blockvectorr))
     if(use_linalg_gpu==1) then
       call copy_on_gpu(blockvectorbx,C_gpu,cplx*dp*vectsize*blocksize)
       call gpu_xgemm(cplx,cparam(cplx),'n',blocksize,blocksize,vectsize,cone,C_gpu,&
&       vectsize,blockvectorr_gpu,vectsize,czero,gram_gpu,blocksize)
       call copy_from_gpu(gramxax,gram_gpu,cplx*dp*blocksize*blocksize)
     else
       call abi_xgemm(cparam(cplx),'n',blocksize,blocksize,vectsize,cone,blockvectorbx,&
&       vectsize,blockvectorr,vectsize,czero,gramxax,blocksize,x_cplx=x_cplx)
     end if

     if(abs(dtset%timopt)==3) then
       call timab(533,1,tsec)
     end if
     call xmpi_sum(gramxax,mpi_enreg%comm_bandspinorfft,ierr)
     if(abs(dtset%timopt)==3) then
       call timab(533,2,tsec)
     end if

     if(use_linalg_gpu==1) then
       call copy_on_gpu(gramxax,gram_gpu,cplx*dp*blocksize*blocksize)
       call copy_on_gpu(blockvectorx,C_gpu,cplx*dp*vectsize*blocksize)
       call gpu_xgemm(cplx,'n','n',vectsize,blocksize,blocksize,cminusone,C_gpu,&
&       vectsize,gram_gpu,blocksize,cone,blockvectorr_gpu,vectsize)
       call copy_from_gpu(blockvectorr,blockvectorr_gpu,cplx*dp*vectsize*blocksize)
     else
       call abi_xgemm('n','n',vectsize,blocksize,blocksize,cminusone,blockvectorx,&
&       vectsize,gramxax,blocksize,cone,blockvectorr,vectsize,x_cplx=x_cplx)
     end if

     ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor*blocksize))
     ABI_ALLOCATE(gwavef,(2,npw_k*my_nspinor*blocksize))
     ABI_ALLOCATE(gvnlxc,(2,npw_k*my_nspinor*blocksize))
     ABI_ALLOCATE(swavef,(2,npw_k*my_nspinor*blocksize))

     call wfcopy('I',vectsize*blocksize,blockvectorr,1,cwavef,1,blocksize,iblock,'W',withbbloc=.false.,&
&     timopt=timopt,tim_wfcopy=tim_wfcopy)

     cpopt=-1;sij_opt=0;if (gen_eigenpb) sij_opt=1

     if(abs(dtset%timopt)==4) then
       call timab(523,2,tsec)
     end if
     if(abs(dtset%timopt)==4) then
       call timab(526,1,tsec)
     end if

     if (mpi_enreg%paral_kgb==0) then
       call getghc(cpopt,cwavef,cprj_dum,gwavef,swavef,gs_hamk,gvnlxc,dum,&
&       mpi_enreg,blocksize,prtvol,sij_opt,tim_getghc,0)
     else
       call prep_getghc(cwavef,gs_hamk,gvnlxc,gwavef,swavef,dum,blocksize,mpi_enreg,&
&       prtvol,sij_opt,cpopt,cprj_dum,already_transposed=.false.)
     end if

     if(abs(dtset%timopt)==4) then
       call timab(526,2,tsec)
     end if
     if(abs(dtset%timopt)==4) then
       call timab(524,1,tsec)
     end if

     if (gen_eigenpb) then
       call wfcopy('D',vectsize*blocksize,swavef,1,blockvectorbr,1,blocksize,iblock,'W',withbbloc=.false.,&
&       timopt=timopt,tim_wfcopy=tim_wfcopy)
     else
       call abi_xcopy(vectsize*blocksize,blockvectorr,1,blockvectorbr,1,x_cplx=x_cplx)
       call wfcopy('D',vectsize*blocksize,gvnlxc,1,blockvectorvr,1,blocksize,iblock,'W',withbbloc=.false.,&
&       timopt=timopt,tim_wfcopy=tim_wfcopy)
     end if

     call wfcopy('D',vectsize*blocksize,gwavef,1,blockvectorar,1,blocksize,iblock,'W',withbbloc=.false.,&
&     timopt=timopt,tim_wfcopy=tim_wfcopy)

     ABI_DEALLOCATE(cwavef)
     ABI_DEALLOCATE(gwavef)
     ABI_DEALLOCATE(gvnlxc)
     ABI_DEALLOCATE(swavef)

     if(use_linalg_gpu==1) then
       call copy_on_gpu(blockvectorbr,blockvectorbr_gpu,cplx*dp*vectsize*blocksize)
       call gpu_xorthonormalize(blockvectorr_gpu,blockvectorbr_gpu,blocksize,mpi_enreg%comm_bandspinorfft,gram_gpu,vectsize,&
&       x_cplx,timopt=timopt,tim_xortho=tim_xortho)
       call copy_from_gpu(blockvectorr,blockvectorr_gpu,cplx*dp*vectsize*blocksize)
       call gpu_xtrsm(cplx,'r','u','n','n',vectsize,blocksize,cone,gram_gpu,blocksize,blockvectorbr_gpu,vectsize)
       call copy_from_gpu(blockvectorbr,blockvectorbr_gpu,cplx*dp*vectsize*blocksize)
       call copy_on_gpu(blockvectorar,blockvectorar_gpu,cplx*dp*vectsize*blocksize)
       call gpu_xtrsm(cplx,'r','u','n','n',vectsize,blocksize,cone,gram_gpu,blocksize,blockvectorar_gpu,vectsize)
       call copy_from_gpu(blockvectorar,blockvectorar_gpu,cplx*dp*vectsize*blocksize)
       if (gs_hamk%usepaw==0) then
         call copy_on_gpu(blockvectorvr,A_gpu,cplx*dp*vectsize*blocksize)
         call gpu_xtrsm(cplx,'r','u','n','n',vectsize,blocksize,cone,gram_gpu,blocksize,A_gpu,vectsize)
         call copy_from_gpu(blockvectorvr,A_gpu,cplx*dp*vectsize*blocksize)
       end if
       call copy_from_gpu(gramrbr,gram_gpu,cplx*dp*blocksize*blocksize)
     else
       call abi_xorthonormalize(blockvectorr,blockvectorbr,blocksize,mpi_enreg%comm_bandspinorfft,gramrbr,vectsize,&
&       x_cplx,timopt=timopt,tim_xortho=tim_xortho)
       call abi_xtrsm('r','u','n','n',vectsize,blocksize,cone,gramrbr,blocksize,blockvectorbr,vectsize,x_cplx=x_cplx)
       call abi_xtrsm('r','u','n','n',vectsize,blocksize,cone,gramrbr,blocksize,blockvectorar,vectsize,x_cplx=x_cplx)
       if (gs_hamk%usepaw==0) then
         call abi_xtrsm('r','u','n','n',vectsize,blocksize,cone,gramrbr,blocksize,blockvectorvr,vectsize,x_cplx=x_cplx)
       end if
     end if

     if(iterationnumber>1) then
       if(use_linalg_gpu==1) then
         call copy_on_gpu(blockvectorp,A_gpu,cplx*dp*vectsize*blocksize)
         call copy_on_gpu(blockvectorbp,C_gpu,cplx*dp*vectsize*blocksize)
         call gpu_xorthonormalize(A_gpu,C_gpu,blocksize,mpi_enreg%comm_bandspinorfft,gram_gpu,vectsize,&
&         x_cplx,timopt=timopt,tim_xortho=tim_xortho)
         call copy_from_gpu(blockvectorp,A_gpu,cplx*dp*vectsize*blocksize)
         call gpu_xtrsm(cplx,'r','u','n','n',vectsize,blocksize,cone,gram_gpu,blocksize,C_gpu,vectsize)
         call copy_from_gpu(blockvectorbp,C_gpu,cplx*dp*vectsize*blocksize)
         call copy_on_gpu(blockvectorap,A_gpu,cplx*dp*vectsize*blocksize)
         call gpu_xtrsm(cplx,'r','u','n','n',vectsize,blocksize,cone,gram_gpu,blocksize,A_gpu,vectsize)
         call copy_from_gpu(blockvectorap,A_gpu,cplx*dp*vectsize*blocksize)
         if (gs_hamk%usepaw==0) then
           call copy_on_gpu(blockvectorvp,A_gpu,cplx*dp*vectsize*blocksize)
           call gpu_xtrsm(cplx,'r','u','n','n',vectsize,blocksize,cone,gram_gpu,blocksize,A_gpu,vectsize)
           call copy_from_gpu(blockvectorvp,A_gpu,cplx*dp*vectsize*blocksize)
         end if
         call copy_from_gpu(grampbp,gram_gpu,cplx*dp*blocksize*blocksize)
       else
!        call orthonormalize(blockvectorp,blockvectorbp,blockvectorap)
         call abi_xorthonormalize(blockvectorp,blockvectorbp,blocksize,mpi_enreg%comm_bandspinorfft,grampbp,vectsize,&
&         x_cplx,timopt=timopt,tim_xortho=tim_xortho)
!        blockvectorap=matmul(blockvectorap,grampbp)
         call abi_xtrsm('r','u','n','n',vectsize,blocksize,cone,grampbp,blocksize,blockvectorbp,vectsize,x_cplx=x_cplx)
         call abi_xtrsm('r','u','n','n',vectsize,blocksize,cone,grampbp,blocksize,blockvectorap,vectsize,x_cplx=x_cplx)
         if (gs_hamk%usepaw==0) then
           call abi_xtrsm('r','u','n','n',vectsize,blocksize,cone,grampbp,blocksize,blockvectorvp,vectsize,x_cplx=x_cplx)
         end if
       end if
     end if

     activersize=blocksize
     if (iterationnumber==1) then
       activepsize=0
       restart=1
     else
       activepsize=blocksize
       restart=0
     end if

!    gramxar=matmul((blockvectorax)^T,blockvectorr)
!    gramrar=matmul((blockvectorar)^T,blockvectorr)
!    gramxax=matmul((blockvectorax)^T,blockvectorx)
     if(use_linalg_gpu==1) then
       call copy_on_gpu(blockvectorax,A_gpu,cplx*dp*vectsize*blocksize)
       call gpu_xgemm(cplx,cparam(cplx),'n',blocksize,blocksize,vectsize,cone,A_gpu,&
&       vectsize,blockvectorr_gpu,vectsize,czero,gram_gpu,blocksize)
       call copy_from_gpu(gramxar,gram_gpu,cplx*dp*blocksize*blocksize)
       call gpu_xgemm(cplx,cparam(cplx),'n',blocksize,blocksize,vectsize,cone,blockvectorar_gpu,&
&       vectsize,blockvectorr_gpu,vectsize,czero,gram_gpu,blocksize)
       call copy_from_gpu(gramrar,gram_gpu,cplx*dp*blocksize*blocksize)
       call copy_on_gpu(blockvectorx,C_gpu,cplx*dp*vectsize*blocksize)
       call gpu_xgemm(cplx,cparam(cplx),'n',blocksize,blocksize,vectsize,cone,A_gpu,&
&       vectsize,C_gpu,vectsize,czero,gram_gpu,blocksize)
       call copy_from_gpu(gramxax,gram_gpu,cplx*dp*blocksize*blocksize)
     else
       call abi_xgemm(cparam(cplx),'n',blocksize,blocksize,vectsize,cone,blockvectorax,&
&       vectsize,blockvectorr,vectsize,czero,gramxar,blocksize,x_cplx=x_cplx)
       call abi_xgemm(cparam(cplx),'n',blocksize,blocksize,vectsize,cone,blockvectorar,&
&       vectsize,blockvectorr,vectsize,czero,gramrar,blocksize,x_cplx=x_cplx)
       call abi_xgemm(cparam(cplx),'n',blocksize,blocksize,vectsize,cone,blockvectorax,&
&       vectsize,blockvectorx,vectsize,czero,gramxax,blocksize,x_cplx=x_cplx)
     end if

     call abi_xcopy(blocksize*blocksize,gramxar,1,transf3(:,:,1),1,x_cplx=x_cplx)
     call abi_xcopy(blocksize*blocksize,gramrar,1,transf3(:,:,2),1,x_cplx=x_cplx)
     call abi_xcopy(blocksize*blocksize,gramxax,1,transf3(:,:,3),1,x_cplx=x_cplx)
     if(abs(dtset%timopt)==3) then
       call timab(533,1,tsec)
     end if
     call xmpi_sum(transf3,mpi_enreg%comm_bandspinorfft,ierr)
     if(abs(dtset%timopt)==3) then
       call timab(533,2,tsec)
     end if

     call abi_xcopy(blocksize*blocksize,transf3(:,:,1),1,gramxar,1,x_cplx=x_cplx)
     call abi_xcopy(blocksize*blocksize,transf3(:,:,2),1,gramrar,1,x_cplx=x_cplx)
     call abi_xcopy(blocksize*blocksize,transf3(:,:,3),1,gramxax,1,x_cplx=x_cplx)

!    gramxbx=matmul((blockvectorbx)^T,blockvectorx)
!    gramrbr=matmul((blockvectorbr)^T,blockvectorr)
!    gramxbr=matmul((blockvectorbx)^T,blockvectorr)
!    Note that the gramb matrix is more easier to construct than grama:
!    i) <x|B|x>=<r|B|r>=<p|B|p>=(1;0)
!    since the x, r and p blockvector are normalized
!    ii) <r|B|x>=(0;0)
!    since the x and r blockvector are orthogonalized
!    iii) The <p|B|r> and <p|B|x> have to be computed.
     gramxbx(:,:)=zero
     gramrbr(:,:)=zero
     gramxbr(:,:)=zero
     do iblocksize=1,blocksize
       zvar=(/one,zero/)
       call abi_xcopy(1,zvar,1,gramxbx(cplx*(iblocksize-1)+1:cplx*iblocksize,iblocksize),1,x_cplx=x_cplx)
       call abi_xcopy(1,zvar,1,gramrbr(cplx*(iblocksize-1)+1:cplx*iblocksize,iblocksize),1,x_cplx=x_cplx)
     end do

!    ###########################################################################
!    ################ PERFORM LOOP ON COND #####################################
!    ###########################################################################

     i1=0;i2=blocksize;i3=2*blocksize;i4=3*blocksize
     cond: do cond_try=1,2 !2 when restart
       if (restart==0) then

!        gramxap=matmul((blockvectorax)^T,blockvectorp)
!        gramrap=matmul((blockvectorar)^T,blockvectorp)
!        grampap=matmul((blockvectorap)^T,blockvectorp)
!        gramxbp=matmul((blockvectorbx)^T,blockvectorp)
!        gramrbp=matmul((blockvectorbr)^T,blockvectorp)
!        grampbp=matmul((blockvectorbp)^T,blockvectorp)
         if(use_linalg_gpu==1) then
           call copy_on_gpu(blockvectorp,C_gpu,cplx*dp*vectsize*blocksize)
           call copy_on_gpu(blockvectorax,A_gpu,cplx*dp*vectsize*blocksize)
           call gpu_xgemm(cplx,cparam(cplx),'n',blocksize,blocksize,vectsize,&
&           cone,A_gpu,vectsize,C_gpu,vectsize,czero,gram_gpu,blocksize)
           call copy_from_gpu(gramxap,gram_gpu,cplx*dp*blocksize*blocksize)
           call gpu_xgemm(cplx,cparam(cplx),'n',blocksize,blocksize,vectsize,&
&           cone,blockvectorar_gpu,vectsize,C_gpu,vectsize,czero,gram_gpu,blocksize)
           call copy_from_gpu(gramrap,gram_gpu,cplx*dp*blocksize*blocksize)
           call copy_on_gpu(blockvectorap,A_gpu,cplx*dp*vectsize*blocksize)
           call gpu_xgemm(cplx,cparam(cplx),'n',blocksize,blocksize,vectsize,&
&           cone,A_gpu,vectsize,C_gpu,vectsize,czero,gram_gpu,blocksize)
           call copy_from_gpu(grampap,gram_gpu,cplx*dp*blocksize*blocksize)
           call copy_on_gpu(blockvectorbx,A_gpu,cplx*dp*vectsize*blocksize)
           call gpu_xgemm(cplx,cparam(cplx),'n',blocksize,blocksize,vectsize,&
&           cone,A_gpu,vectsize,C_gpu,vectsize,czero,gram_gpu,blocksize)
           call copy_from_gpu(gramxbp,gram_gpu,cplx*dp*blocksize*blocksize)
           call gpu_xgemm(cplx,cparam(cplx),'n',blocksize,blocksize,vectsize,&
&           cone,blockvectorbr_gpu,vectsize,C_gpu,vectsize,czero,gram_gpu,blocksize)
           call copy_from_gpu(gramrbp,gram_gpu,cplx*dp*blocksize*blocksize)
         else
           call abi_xgemm(cparam(cplx),'n',blocksize,blocksize,vectsize,cone,blockvectorax,&
&           vectsize,blockvectorp,vectsize,czero,gramxap,blocksize,x_cplx=x_cplx)
           call abi_xgemm(cparam(cplx),'n',blocksize,blocksize,vectsize,cone,blockvectorar,&
&           vectsize,blockvectorp,vectsize,czero,gramrap,blocksize,x_cplx=x_cplx)
           call abi_xgemm(cparam(cplx),'n',blocksize,blocksize,vectsize,cone,blockvectorap,&
&           vectsize,blockvectorp,vectsize,czero,grampap,blocksize,x_cplx=x_cplx)
           call abi_xgemm(cparam(cplx),'n',blocksize,blocksize,vectsize,cone,blockvectorbx,&
&           vectsize,blockvectorp,vectsize,czero,gramxbp,blocksize,x_cplx=x_cplx)
           call abi_xgemm(cparam(cplx),'n',blocksize,blocksize,vectsize,cone,blockvectorbr,&
&           vectsize,blockvectorp,vectsize,czero,gramrbp,blocksize,x_cplx=x_cplx)
         end if
!        It's not necessary to compute the last one: <p|B|p>=(1;0) (see above)
         transf5(:,:,1)=gramxap(:,:)
         transf5(:,:,2)=gramrap(:,:)
         transf5(:,:,3)=grampap(:,:)
         transf5(:,:,4)=gramxbp(:,:)
         transf5(:,:,5)=gramrbp(:,:)
         if(abs(dtset%timopt)==3) then
           call timab(533,1,tsec)
         end if
         call xmpi_sum(transf5,mpi_enreg%comm_bandspinorfft,ierr)
         if(abs(dtset%timopt)==3) then
           call timab(533,2,tsec)
         end if
         gramxap(:,:)=transf5(:,:,1)
         gramrap(:,:)=transf5(:,:,2)
         grampap(:,:)=transf5(:,:,3)
         gramxbp(:,:)=transf5(:,:,4)
         gramrbp(:,:)=transf5(:,:,5)
         grampbp(:,:)=zero
         do iblocksize=1,blocksize
           zvar=(/one,zero/)
           call abi_xcopy(1,zvar,1,grampbp(cplx*(iblocksize-1)+1:cplx*iblocksize,iblocksize),1,x_cplx=x_cplx)
         end do
         bigorder=i4
         ABI_ALLOCATE(grama,(cplx*i4,i4))
         ABI_ALLOCATE(gramb,(cplx*i4,i4))
         ABI_ALLOCATE(eigen,(i4))
!        ABI_ALLOCATE(coordx,(cplx*i4,blocksize))
         ABI_ALLOCATE(coordx1,(cplx*blocksize,blocksize))
         ABI_ALLOCATE(coordx2,(cplx*blocksize,blocksize))
         ABI_ALLOCATE(coordx3,(cplx*blocksize,blocksize))
         grama(:,:)=zero;gramb(:,:)=zero
         grama(gramindex(i1+1):gramindex(i2)+cplx-1,i1+1:i2)=gramxax
         grama(gramindex(i1+1):gramindex(i2)+cplx-1,i2+1:i3)=gramxar
         grama(gramindex(i1+1):gramindex(i2)+cplx-1,i3+1:i4)=gramxap
         grama(gramindex(i2+1):gramindex(i3)+cplx-1,i2+1:i3)=gramrar
         grama(gramindex(i2+1):gramindex(i3)+cplx-1,i3+1:i4)=gramrap
         grama(gramindex(i3+1):gramindex(i4)+cplx-1,i3+1:i4)=grampap
         gramb(gramindex(i1+1):gramindex(i2)+cplx-1,i1+1:i2)=gramxbx
         gramb(gramindex(i1+1):gramindex(i2)+cplx-1,i2+1:i3)=gramxbr
         gramb(gramindex(i1+1):gramindex(i2)+cplx-1,i3+1:i4)=gramxbp
         gramb(gramindex(i2+1):gramindex(i3)+cplx-1,i2+1:i3)=gramrbr
         gramb(gramindex(i2+1):gramindex(i3)+cplx-1,i3+1:i4)=gramrbp
         gramb(gramindex(i3+1):gramindex(i4)+cplx-1,i3+1:i4)=grampbp
       else
         bigorder=i3
         ABI_ALLOCATE(grama,(cplx*i3,i3))
         ABI_ALLOCATE(gramb,(cplx*i3,i3))
         ABI_ALLOCATE(eigen,(i3))
!        ABI_ALLOCATE(coordx,(cplx*i3,blocksize))
         ABI_ALLOCATE(coordx1,(cplx*blocksize,blocksize))
         ABI_ALLOCATE(coordx2,(cplx*blocksize,blocksize))
         grama(:,:)=zero;gramb(:,:)=zero
         grama(gramindex(i1+1):gramindex(i2)+cplx-1,i1+1:i2)=gramxax
         grama(gramindex(i1+1):gramindex(i2)+cplx-1,i2+1:i3)=gramxar
         grama(gramindex(i2+1):gramindex(i3)+cplx-1,i2+1:i3)=gramrar
         gramb(gramindex(i1+1):gramindex(i2)+cplx-1,i1+1:i2)=gramxbx
         gramb(gramindex(i1+1):gramindex(i2)+cplx-1,i2+1:i3)=gramxbr
         gramb(gramindex(i2+1):gramindex(i3)+cplx-1,i2+1:i3)=gramrbr
       end if

       ABI_ALLOCATE(tmpgramb,(cplx*bigorder,bigorder))
       ABI_ALLOCATE(tmpeigen,(bigorder))
       tmpgramb=gramb

       call abi_xheev('v','u',bigorder,tmpgramb,bigorder,tmpeigen,x_cplx=cplx,istwf_k=istwf_k, &
&       timopt=timopt,tim_xeigen=tim_xeigen,use_slk=dtset%use_slk,use_gpu=use_lapack_gpu)

       condestgramb=tmpeigen(bigorder)/tmpeigen(1)
       ABI_DEALLOCATE(tmpgramb)
       ABI_DEALLOCATE(tmpeigen)

       if (condestgramb.gt.1d+5.or.condestgramb.lt.0.d0.or.info/=0) then
         write(std_out,*)'condition number of the Gram matrix = ',condestgramb
         if (cond_try==1.and.restart==0) then
           ABI_DEALLOCATE(grama)
           ABI_DEALLOCATE(gramb)
           ABI_DEALLOCATE(eigen)
!          ABI_DEALLOCATE(coordx)
           ABI_DEALLOCATE(coordx1)
           ABI_DEALLOCATE(coordx2)
           if(bigorder==i4) then
             ABI_DEALLOCATE(coordx3)
           end if
           if (nrestart.gt.1) then
             MSG_WARNING('the minimization is stopped for this block')
             exit iter
           else
             restart=1
             nrestart=nrestart+1
             call wrtout(std_out,'Lobpcgwf: restart performed',"PERS")
           end if
         else
           MSG_WARNING('Gramm matrix ill-conditionned: results may be unpredictable')
         end if
       else
         exit cond
       end if
     end do cond

!    ###########################################################################
!    ################ END LOOP ON COND #########################################
!    ###########################################################################

     call abi_xhegv(1,'v','u',bigorder,grama,bigorder,gramb,bigorder,eigen,x_cplx=cplx,istwf_k=istwf_k, &
     timopt=timopt,tim_xeigen=tim_xeigen,use_slk=dtset%use_slk,use_gpu=use_lapack_gpu)

     deltae=-one
     do iblocksize=1,blocksize
       call abi_xcopy(1,lambda(cplx*(iblocksize-1)+1,iblocksize),1,zvar,1,x_cplx=x_cplx)
       deltae=max(deltae,abs(cmplx(zvar(1),zvar(2))-eigen(iblocksize)))
       zvar=(/eigen(iblocksize),zero/)
       call abi_xcopy(1,zvar,1,lambda(cplx*(iblocksize-1)+1,iblocksize),1,x_cplx=x_cplx)
     end do

!    DEBUG
!    write(std_out,*)'eigen',eigen(1:blocksize)
!    ENDDEBUG

!    coordx(1:bigorder*cplx,1:blocksize)=grama(1:bigorder*cplx,1:blocksize)
     coordx1(:,:) =  grama(1+i1*cplx : i2*cplx,1:blocksize)
     coordx2(:,:) =  grama(1+i2*cplx : i3*cplx,1:blocksize)
     if(bigorder==i4) then
       coordx3(:,:) =  grama(1+i3*cplx : i4*cplx,1:blocksize)
     end if


     if(use_linalg_gpu==1) then
       call copy_on_gpu(coordx2,coordx2_gpu,cplx*dp*blocksize*blocksize)
       if(bigorder==i4) then
         call copy_on_gpu(coordx3,coordx3_gpu,cplx*dp*blocksize*blocksize)
       end if
     end if

     ABI_DEALLOCATE(grama)
     ABI_DEALLOCATE(gramb)
     ABI_DEALLOCATE(eigen)
     if (restart==0 .and. iterationnumber >1) then

!      blockvectorp=matmul(blockvectorr,coordx(i2+1:i3,:))+&
!      &               matmul(blockvectorp,coordx(i3+1:i4,:))
       if(use_linalg_gpu==1) then
!        call copy_on_gpu(blockvectorr,blockvectorr_gpu,cplx*dp*vectsize*blocksize)
         call gpu_xgemm(cplx,'n','n',vectsize,blocksize,blocksize,cone,blockvectorr_gpu,&
&         vectsize,coordx2_gpu,blocksize,czero,C_gpu,vectsize)
         call copy_on_gpu(blockvectorp,A_gpu,cplx*dp*vectsize*blocksize)
         call gpu_xgemm(cplx,'n','n',vectsize,blocksize,blocksize,cone,A_gpu,vectsize,&
&         coordx3_gpu,blocksize,cone,C_gpu,vectsize)
         call copy_from_gpu(blockvectorp,C_gpu,cplx*dp*vectsize*blocksize)
       else
         call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorr,&
&         vectsize,coordx2,blocksize,czero,blockvectordumm,vectsize,x_cplx=x_cplx)
         call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorp,&
&         vectsize,coordx3,blocksize,cone,blockvectordumm,vectsize,x_cplx=x_cplx)
         call abi_xcopy(vectsize*blocksize,blockvectordumm,1,blockvectorp,1,x_cplx=x_cplx)
       end if

!      blockvectorap=matmul(blockvectorar,coordx(i2+1:i3,:))+&
!      &                matmul(blockvectorap,coordx(i3+1:i4,:))
       if(use_linalg_gpu==1) then
!        call copy_on_gpu(blockvectorar,blockvectorar_gpu,cplx*dp*vectsize*blocksize)
         call gpu_xgemm(cplx,'n','n',vectsize,blocksize,blocksize,cone,blockvectorar_gpu,&
&         vectsize,coordx2_gpu,blocksize,czero,C_gpu,vectsize)
         call copy_on_gpu(blockvectorap,A_gpu,cplx*dp*vectsize*blocksize)
         call gpu_xgemm(cplx,'n','n',vectsize,blocksize,blocksize,cone,A_gpu,vectsize,&
&         coordx3_gpu,blocksize,cone,C_gpu,vectsize)
         call copy_from_gpu(blockvectorap,C_gpu,cplx*dp*vectsize*blocksize)
       else
         call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorar,&
&         vectsize,coordx2,blocksize,czero,blockvectordumm,vectsize,x_cplx=x_cplx)
         call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorap,&
&         vectsize,coordx3,blocksize,cone,blockvectordumm,vectsize,x_cplx=x_cplx)
         call abi_xcopy(vectsize*blocksize,blockvectordumm,1,blockvectorap,1,x_cplx=x_cplx)
       end if


!      blockvectorvp=matmul(blockvectorvr,coordx(i2+1:i3,:))+&
!      &                matmul(blockvectorvp,coordx(i3+1:i4,:))
       if (gs_hamk%usepaw==0) then
         if(use_linalg_gpu==1) then
           call copy_on_gpu(blockvectorvr,A_gpu,cplx*dp*vectsize*blocksize)
           call gpu_xgemm(cplx,'n','n',vectsize,blocksize,blocksize,cone,A_gpu,&
&           vectsize,coordx2_gpu,blocksize,czero,C_gpu,vectsize)
           call copy_on_gpu(blockvectorvp,A_gpu,cplx*dp*vectsize*blocksize)
           call gpu_xgemm(cplx,'n','n',vectsize,blocksize,blocksize,cone,A_gpu,&
&           vectsize,coordx3_gpu,blocksize,cone,C_gpu,vectsize)
           call copy_from_gpu(blockvectorvp,C_gpu,cplx*dp*vectsize*blocksize)
         else
           call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorvr,&
&           vectsize,coordx2,blocksize,czero,blockvectordumm,vectsize,x_cplx=x_cplx)
           call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorvp,&
&           vectsize,coordx3,blocksize,cone,blockvectordumm,vectsize,x_cplx=x_cplx)
           call abi_xcopy(vectsize*blocksize,blockvectordumm,1,blockvectorvp,1,x_cplx=x_cplx)
         end if
       end if

!      blockvectorbp=matmul(blockvectorbr,coordx(i2+1:i3,:))+&
!      &                matmul(blockvectorbp,coordx(i3+1:i4,:))
       if(use_linalg_gpu==1) then
!        call copy_on_gpu(blockvectorbr,blockvectorbr_gpu,cplx*dp*vectsize*blocksize)
         call gpu_xgemm(cplx,'n','n',vectsize,blocksize,blocksize,cone,blockvectorbr_gpu,&
&         vectsize,coordx2_gpu,blocksize,czero,C_gpu,vectsize)
         call copy_on_gpu(blockvectorbp,A_gpu,cplx*dp*vectsize*blocksize)
         call gpu_xgemm(cplx,'n','n',vectsize,blocksize,blocksize,cone,A_gpu,vectsize,&
&         coordx3_gpu,blocksize,cone,C_gpu,vectsize)
         call copy_from_gpu(blockvectorbp,C_gpu,cplx*dp*vectsize*blocksize)
       else
         call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbr,&
&         vectsize,coordx2,blocksize,czero,blockvectordumm,vectsize,x_cplx=x_cplx)
         call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbp,&
&         vectsize,coordx3,blocksize,cone,blockvectordumm,vectsize,x_cplx=x_cplx)
         call abi_xcopy(vectsize*blocksize,blockvectordumm,1,blockvectorbp,1,x_cplx=x_cplx)
       end if

     else

!      blockvectoSz =matmul(blockvectorr,coordx(i2+1:i3,:))
       if(use_linalg_gpu==1) then
!        call copy_on_gpu(blockvectorr,blockvectorr_gpu,cplx*dp*vectsize*blocksize)
         call gpu_xgemm(cplx,'n','n',vectsize,blocksize,blocksize,cone,blockvectorr_gpu,&
&         vectsize,coordx2_gpu,blocksize,czero,C_gpu,vectsize)
         call copy_from_gpu(blockvectorp,C_gpu,cplx*dp*vectsize*blocksize)
       else
         call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorr,&
&         vectsize,coordx2,blocksize,czero,blockvectorp,vectsize,x_cplx=x_cplx)
       end if

!      blockvectorap=matmul(blockvectorar,coordx(i2+1:i3,:))
       if(use_linalg_gpu==1) then
!        call copy_on_gpu(blockvectorar,blockvectorar_gpu,cplx*dp*vectsize*blocksize)
         call gpu_xgemm(cplx,'n','n',vectsize,blocksize,blocksize,cone,blockvectorar_gpu,&
&         vectsize,coordx2_gpu,blocksize,czero,C_gpu,vectsize)
         call copy_from_gpu(blockvectorap,C_gpu,cplx*dp*vectsize*blocksize)
       else
         call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorar,&
&         vectsize,coordx2,blocksize,czero,blockvectorap,vectsize,x_cplx=x_cplx)
       end if
!      blockvectorvp=matmul(blockvectorvr,coordx(i2+1:i3,:))
       if (gs_hamk%usepaw==0) then
         if(use_linalg_gpu==1) then
           call copy_on_gpu(blockvectorvr,A_gpu,cplx*dp*vectsize*blocksize)
           call gpu_xgemm(cplx,'n','n',vectsize,blocksize,blocksize,cone,A_gpu,&
&           vectsize,coordx2_gpu,blocksize,czero,C_gpu,vectsize)
           call copy_from_gpu(blockvectorvp,C_gpu,cplx*dp*vectsize*blocksize)
         else
           call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorvr,&
&           vectsize,coordx2,blocksize,czero,blockvectorvp,vectsize,x_cplx=x_cplx)
         end if
       end if

!      blockvectorbp=matmul(blockvectorbr,coordx(i2+1:i3,:))
       if(use_linalg_gpu==1) then
!        call copy_on_gpu(blockvectorbr,blockvectorbr_gpu,cplx*dp*vectsize*blocksize)
         call gpu_xgemm(cplx,'n','n',vectsize,blocksize,blocksize,cone,blockvectorbr_gpu,&
&         vectsize,coordx2_gpu,blocksize,czero,C_gpu,vectsize)
         call copy_from_gpu(blockvectorbp,C_gpu,cplx*dp*vectsize*blocksize)
       else
         call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbr,&
&         vectsize,coordx2,blocksize,czero,blockvectorbp,vectsize,x_cplx=x_cplx)
       end if
     end if

     if(use_linalg_gpu==1) then
       call copy_on_gpu(coordx1,coordx2_gpu,cplx*dp*blocksize*blocksize)
     end if

!    blockvectorx = matmul(blockvectorx,coordx(i1+1:i2,:))+blockvectorp
     if(use_linalg_gpu==1) then
       call copy_on_gpu(blockvectorx,A_gpu,cplx*dp*vectsize*blocksize)
       call gpu_xgemm(cplx,'n','n',vectsize,blocksize,blocksize,cone,A_gpu,&
&       vectsize,coordx2_gpu,blocksize,czero,C_gpu,vectsize)
       call copy_from_gpu(blockvectordumm,C_gpu,cplx*dp*vectsize*blocksize)
     else
       call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
&       vectsize,coordx1,blocksize,czero,blockvectordumm,vectsize,x_cplx=x_cplx)
     end if
     blockvectorx = blockvectordumm+blockvectorp

!    blockvectorax= matmul(blockvectorax,coordx(i1+1:i2,:))+blockvectorap
     if(use_linalg_gpu==1) then
       call copy_on_gpu(blockvectorax,A_gpu,cplx*dp*vectsize*blocksize)
       call gpu_xgemm(cplx,'n','n',vectsize,blocksize,blocksize,cone,A_gpu,&
&       vectsize,coordx2_gpu,blocksize,czero,C_gpu,vectsize)
       call copy_from_gpu(blockvectordumm,C_gpu,cplx*dp*vectsize*blocksize)
     else
       call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorax,&
&       vectsize,coordx1,blocksize,czero,blockvectordumm,vectsize,x_cplx=x_cplx)
     end if
     blockvectorax = blockvectordumm+blockvectorap

!    blockvectorvx= matmul(blockvectorvx,coordx(i1+1:i2,:))+blockvectorvp
     if (gs_hamk%usepaw==0) then
       if(use_linalg_gpu==1) then
         call copy_on_gpu(blockvectorvx,A_gpu,cplx*dp*vectsize*blocksize)
         call gpu_xgemm(cplx,'n','n',vectsize,blocksize,blocksize,cone,A_gpu,&
&         vectsize,coordx2_gpu,blocksize,czero,C_gpu,vectsize)
         call copy_from_gpu(blockvectordumm,C_gpu,cplx*dp*vectsize*blocksize)
       else
         call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorvx,&
&         vectsize,coordx1,blocksize,czero,blockvectordumm,vectsize,x_cplx=x_cplx)
       end if
       blockvectorvx = blockvectordumm+blockvectorvp
     end if

!    blockvectorbx= matmul(blockvectorbx,coordx(i1+1:i2,:))+blockvectorbp
     if(use_linalg_gpu==1) then
       call copy_on_gpu(blockvectorbx,A_gpu,cplx*dp*vectsize*blocksize)
       call gpu_xgemm(cplx,'n','n',vectsize,blocksize,blocksize,cone,A_gpu,&
&       vectsize,coordx2_gpu,blocksize,czero,C_gpu,vectsize)
       call copy_from_gpu(blockvectordumm,C_gpu,cplx*dp*vectsize*blocksize)
     else
       call abi_xgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbx,&
&       vectsize,coordx1,blocksize,czero,blockvectordumm,vectsize,x_cplx=x_cplx)
     end if
     blockvectorbx = blockvectordumm+blockvectorbp

!    ABI_DEALLOCATE(coordx)
     ABI_DEALLOCATE(coordx1)
     ABI_DEALLOCATE(coordx2)
     if(bigorder==i4) then
       ABI_DEALLOCATE(coordx3)
     end if

!    Check convergence on energy and eventually exit
     if (iterationnumber==1) then
       deold=deltae
     else if (iterationnumber>1) then
       if ((abs(deltae)<0.005*abs(deold)).and.(iterationnumber/=maxiterations))then
         if(prtvol>=10)then
           write(message, '(2(a,i4),1x,a,1p,e12.4,a,e12.4,a)' ) &
&           ' lobpcgwf: block',iblock,', line',iterationnumber,&
&           ', deltae=',deltae,' < 0.005*',deold,' =>skip lines !'
           call wrtout(std_out,message,'PERS')
         end if
         exit
       else if (abs(deltae)>0.005*abs(deold)) then
         if(prtvol>=10)then
           write(message, '(2(a,i4),1x,a,1p,e12.4,a,e12.4,a)' ) &
&           ' lobpcgwf: block',iblock,', line',iterationnumber,&
&           ', deltae=',deltae,' > 0.005*',deold,' =>keep on working !'
           call wrtout(std_out,message,'PERS')
         end if
       end if
     end if

     if(abs(dtset%timopt)==4) then
       call timab(524,2,tsec)
     end if

   end do iter

!  ###########################################################################
!  ################## END LOOP ON NLINE ######################################
!  ###########################################################################

   if(abs(dtset%timopt)==4) then
     call timab(525,1,tsec)
   end if

   if (havetoprecon) then
     call xprecon(blockvectorbx,lambda,blocksize,&
&     iterationnumber,kinpw,mpi_enreg,npw_k,my_nspinor,&
&     optekin,optpcon,pcon,blockvectorax,blockvectorr,vectsize,timopt=timopt,tim_xprecon=tim_xprecon)

     residualnorms=sum(abs(blockvectorr)**2,dim=1)

     if(abs(dtset%timopt)==3) then
       call timab(533,1,tsec)
     end if
     call xmpi_sum(residualnorms,mpi_enreg%comm_bandspinorfft,ierr)
     if(abs(dtset%timopt)==3) then
       call timab(533,2,tsec)
     end if

     resid_k(bblocksize+1:bblocksize+blocksize)=residualnorms(1:blocksize)
   end if

   call wfcopy('I',vectsize*blocksize,blockvectorx,1,cg,1,blocksize,iblock,'C',withbbloc=.true.,&
&   timopt=timopt,tim_wfcopy=tim_wfcopy)

   if(gen_eigenpb) then
     call wfcopy('I',vectsize*blocksize,blockvectorbx,1,gsc,1,blocksize,iblock,'S',withbbloc=.true.,&
&     timopt=timopt,tim_wfcopy=tim_wfcopy)
   end if

!  The Vnl+VFockACE part of the Hamiltonian is no more stored in the packed form such as it was the case for subvnlx(:).
!  Now, the full matrix is stored in totvnlx(:,:). This trick permits:
!  1) to avoid the reconstruction of the total matrix in vtowfk.F90 (double loop over bands)
!  2) to use two optimized matrix-matrix blas routine for general (in lobpcgccwf.F90) or hermitian (in vtowfk.F90)
!  operators, zgemm.f and zhemm.f respectively, rather than a triple loop in both cases.
   iwavef=iblock*blocksize
   isubh=1+2*bblocksize*(bblocksize+1)/2

   ABI_ALLOCATE(blockvectorz,(cplx*vectsize,iwavef))
   if(bblocksize > 0 ) then
     call abi_xcopy(bblocksize*vectsize,blockvectory(:,1:bblocksize),1,blockvectorz(:,1:bblocksize),1,x_cplx=x_cplx)
   end if
   call abi_xcopy( blocksize*vectsize,blockvectorx(:,1:blocksize) ,1,blockvectorz(:,bblocksize+1:iwavef),1,x_cplx=x_cplx)

   ABI_ALLOCATE(tsubham,(cplx*iwavef,blocksize))
   tsubham(:,:)=zero
   call abi_xgemm(cparam(cplx),'n',iwavef,blocksize,vectsize,cone,blockvectorz,vectsize,&
&   blockvectorax,vectsize,czero,tsubham,iwavef,x_cplx=x_cplx)

   if (gs_hamk%usepaw==0) then
     ! MG FIXME: Here gfortran4.9 allocates temporary array for C in abi_d2zgemm.
     call abi_xgemm(cparam(cplx),'n',blocksize,iwavef,vectsize,cone,blockvectorvx,vectsize,&
&     blockvectorz,vectsize,czero,totvnlx(cplx*bblocksize+1:cplx*iwavef,1:iwavef),blocksize,x_cplx=x_cplx)
   end if

   do iblocksize=1,blocksize
     do ii=1,bblocksize+iblocksize
       if ( cplx == 1 ) then
         subham(isubh)  = tsubham(ii,iblocksize)
         subham(isubh+1)= zero
       else
         subham(isubh)  = tsubham(2*ii-1,iblocksize)
         subham(isubh+1)= tsubham(2*ii  ,iblocksize)
       end if
       isubh=isubh+2
     end do
   end do
   ABI_DEALLOCATE(tsubham)
   ABI_DEALLOCATE(blockvectorz)
!  comm for subham and subvnlx are made in vtowfk

   ABI_DEALLOCATE(pcon)
   ABI_DEALLOCATE(blockvectory)
   ABI_DEALLOCATE(blockvectorby)
   ABI_DEALLOCATE(gramyx)
   ABI_DEALLOCATE(blockvectorx)
   ABI_DEALLOCATE(blockvectorax)
   ABI_DEALLOCATE(blockvectorbx)
   ABI_DEALLOCATE(blockvectorr)
   ABI_DEALLOCATE(blockvectorar)
   ABI_DEALLOCATE(blockvectorbr)
   ABI_DEALLOCATE(blockvectorp)
   ABI_DEALLOCATE(blockvectorap)
   ABI_DEALLOCATE(blockvectorbp)
   if (gs_hamk%usepaw==0) then
     ABI_DEALLOCATE(blockvectorvx)
     ABI_DEALLOCATE(blockvectorvp)
     ABI_DEALLOCATE(blockvectorvr)
   end if
   ABI_DEALLOCATE(blockvectordumm)
   ABI_DEALLOCATE(gramxax)
   ABI_DEALLOCATE(gramxar)
   ABI_DEALLOCATE(gramxap)
   ABI_DEALLOCATE(gramrar)
   ABI_DEALLOCATE(gramrap)
   ABI_DEALLOCATE(grampap)
   ABI_DEALLOCATE(gramxbx)
   ABI_DEALLOCATE(gramxbr)
   ABI_DEALLOCATE(gramxbp)
   ABI_DEALLOCATE(gramrbr)
   ABI_DEALLOCATE(gramrbp)
   ABI_DEALLOCATE(grampbp)
   ABI_DEALLOCATE(transf3)
   ABI_DEALLOCATE(transf5)
   ABI_DEALLOCATE(lambda)
   ABI_DEALLOCATE(residualnorms)
   if(use_linalg_gpu==1) then
     call dealloc_on_gpu(bblockvector_gpu)
     call dealloc_on_gpu(gram_gpu)
   end if

 end do  ! End big loop over bands inside blocks

 if(use_linalg_gpu==1) then
   call dealloc_on_gpu(blockvectorr_gpu)
   call dealloc_on_gpu(blockvectorar_gpu)
   call dealloc_on_gpu(blockvectorbr_gpu)
   call dealloc_on_gpu(A_gpu)
   call dealloc_on_gpu(C_gpu)
   call dealloc_on_gpu(coordx2_gpu)
   call dealloc_on_gpu(coordx3_gpu)
   !call gpu_linalg_shutdown()
 end if

 if(abs(dtset%timopt)==4) then
   call timab(525,2,tsec)
 end if
 call timab(530,2,tsec)

 DBG_ENTER("COLL")

 contains

   function gramindex(iblocksize)

   integer :: gramindex,iblocksize
   gramindex=(iblocksize-1)*cplx+1
 end function gramindex

end subroutine lobpcgwf
!!***

end module m_lobpcgwf_old
!!***
