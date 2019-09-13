!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_chebfiwf
!! NAME
!! m_chebfiwf
!!
!! FUNCTION
!! this routine updates the whole wave functions at a given k-point,
!! using the lobpcg method
!! for a given spin-polarization, from a fixed hamiltonian
!! but might also simply compute eigenvectors and eigenvalues at this k point.
!! it will also update the matrix elements of the hamiltonian.
!!
!! COPYRIGHT
!! Copyright (C) 2018-2019 ABINIT group (BS)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! PARENTS
!!      vtowfk
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_chebfiwf

  use defs_abitypes
  use defs_basis
  use m_abicore
  use m_chebfi
  use m_xmpi
  use m_errors
  use m_time
  use m_xomp
  use m_fstrings
  use m_xg
  use m_chebfi2
 
  use m_hamiltonian, only : gs_hamiltonian_type
  use m_pawcprj,     only : pawcprj_type
  use m_nonlop,      only : nonlop
  use m_prep_kgb,    only : prep_getghc, prep_nonlop
  use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_free
  use m_getghc,      only : multithreaded_getghc
  use m_invovl
  
  use iso_c_binding, only: c_associated, c_loc, c_ptr

  implicit none

  private
 
  integer, parameter :: l_tim_getghc=7
  double precision, parameter :: inv_sqrt2 = 1/sqrt2
  ! Use in getghc_gsc1
  integer,save  :: l_cpopt
  integer,save  :: l_icplx
  integer,save  :: l_istwf
  integer,save  :: l_npw
  integer,save  :: l_nband_filter
  integer,save  :: l_nspinor
  logical,save  :: l_paw
  integer,save  :: l_prtvol
  integer,save  :: l_sij_opt
  integer, save :: l_paral_kgb
  integer, save :: l_useria
  real(dp), allocatable,save :: l_gvnlc(:,:)
  real(dp), allocatable,save ::  l_pcon(:)
  type(mpi_type),pointer,save :: l_mpi_enreg
  type(gs_hamiltonian_type),pointer,save :: l_gs_hamk
  
  integer, parameter :: DEBUG_ROWS = 5
  integer, parameter :: DEBUG_COLUMNS = 5
  
  type(c_ptr) :: cptr
 
  public :: chebfiwf2
 
  contains
 
  subroutine chebfiwf2(cg,dtset,eig,enl_out,gs_hamk,kinpw,mpi_enreg,&
&                   nband,npw,nspinor,prtvol,resid,counter) 


  use m_cgtools, only : dotprod_g
  use iso_c_binding

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chebfiwf2'
!End of the abilint section

  implicit none
 
  !Arguments ------------------------------------
  integer,intent(in) :: nband,npw,prtvol,nspinor
  type(gs_hamiltonian_type),target,intent(inout) :: gs_hamk
  type(dataset_type)              ,intent(in   ) :: dtset
  type(mpi_type)           ,target,intent(inout) :: mpi_enreg
  real(dp)                 ,target,intent(inout) :: cg(2,npw*nspinor*nband)!,gsc(2,nspinor*nband*npw)
  real(dp)                        ,intent(in   ) :: kinpw(npw)
  real(dp)                 ,target,intent(  out) :: resid(nband)
  real(dp)                        ,intent(  out) :: enl_out(nband) !*(1-gs_hamk%usepaw)
  real(dp)                 ,target,intent(  out) :: eig(nband)
  integer                         ,intent(in   ) :: counter            
 
  !Local variables-------------------------------
 
  type(xgBlock_t) :: xgx0
  type(xgBlock_t) :: xgeigen
  type(xgBlock_t) :: xgresidu
  type(chebfi_t) :: chebfi
  
  type(xgBlock_t) :: HELPER
 
  integer :: space, blockdim, nline !blockdim
  integer :: ipw

  integer :: nthreads
 
  double precision :: chebfiMem(2)
  double precision :: localmem
 
  integer, parameter :: tim_chebfiwf2 = 1750
  double precision :: cputime, walltime
  double precision :: tsec(2)
 
  type(c_ptr) :: cptr
  real(dp), pointer :: eig_ptr(:,:) => NULL()
  real(dp), pointer :: resid_ptr(:,:) => NULL()
 
  ! Stupid things for NC
  integer,parameter :: choice=1, paw_opt=0, signs=1
  type(pawcprj_type) :: cprj_dum(gs_hamk%natom,0)
  integer :: iband, shift
  real(dp) :: gsc_dummy(0,0)
  
! *********************************************************************


!###########################################################################
!################ INITIALIZATION  ##########################################
!###########################################################################

  call timab(tim_chebfiwf2,1,tsec)
  cputime = abi_cpu_time()
  walltime = abi_wtime()
 
  ! Set module variables
  l_paw = (gs_hamk%usepaw==1) 
  l_cpopt=-1;l_sij_opt=0;if (l_paw) l_sij_opt=1 
  l_istwf=gs_hamk%istwf_k 
  l_npw = npw
  !print * , "l_npw", l_npw
  !stop
  l_nspinor = nspinor
  l_prtvol = prtvol
  l_mpi_enreg => mpi_enreg
  l_gs_hamk => gs_hamk
  l_nband_filter = nband
  l_paral_kgb = dtset%paral_kgb
 
  !Variables
  nline=dtset%nline
  blockdim=l_mpi_enreg%nproc_band*l_mpi_enreg%bandpp
  !for debug
  l_useria=dtset%useria
  !rint *, "BANDPP", l_mpi_enreg%nproc_band*l_mpi_enreg%bandpp
  !top  

  !Depends on istwfk
  if ( l_istwf == 2 ) then ! Real only
    ! SPACE_CR mean that we have complex numbers but no re*im terms only re*re
    ! and im*im so that a vector of complex is consider as a long vector of real
    ! therefore the number of data is (2*npw*nspinor)*nband
    ! This space is completely equivalent to SPACE_R but will correctly set and
    ! get the array data into the xgBlock
    space = SPACE_CR
    l_icplx = 2
    
  else ! complex
    space = SPACE_C
    l_icplx = 1
  end if
    
  !print *, "SPACEEEEEEEEEEEEEEE", l_icplx
  !stop
 
  ! Memory info
  if ( prtvol >= 3 ) then
    chebfiMem = chebfi_memInfo(nband, l_icplx*l_npw*l_nspinor, space) !blockdim
    localMem = (l_npw+2*l_npw*l_nspinor+2*nband)*kind(1.d0) !blockdim
    write(std_out,'(1x,A,F10.6,1x,A)') "Each MPI process calling chebfi should need around ", &
    (localMem+sum(chebfiMem))/1e9, &
    "GB of peak memory as follows :"
    write(std_out,'(4x,A,F10.6,1x,A)') "Permanent memory in chebfiwf : ", &
    (localMem)/1e9, "GB"
    write(std_out,'(4x,A,F10.6,1x,A)') "Permanent memory in m_chebfi : ", &
    (chebfiMem(1))/1e9, "GB"
    write(std_out,'(4x,A,F10.6,1x,A)') "Temporary memory in m_chebfi : ", &
    (chebfiMem(2))/1e9, "GB"
  end if
 
  !For preconditionning 
  ABI_MALLOC(l_pcon,(1:l_icplx*npw))

  !$omp parallel do schedule(static), shared(l_pcon,kinpw)
  do ipw=1-1,l_icplx*npw-1
    if(kinpw(ipw/l_icplx+1)>huge(0.0_dp)*1.d-11) then
      l_pcon(ipw+1)=0.d0
    else
      l_pcon(ipw+1) = (27+kinpw(ipw/l_icplx+1)*(18+kinpw(ipw/l_icplx+1)*(12+8*kinpw(ipw/l_icplx+1)))) &
&    / (27+kinpw(ipw/l_icplx+1)*(18+kinpw(ipw/l_icplx+1)*(12+8*kinpw(ipw/l_icplx+1))) + 16*kinpw(ipw/l_icplx+1)**4)
    end if
  end do
       
  
  !Local variables for chebfi
  !call xg_init(xgx0,space,icplx*npw*nspinor,nband)
!  print *, "l_icplx", l_icplx
!  print *, "l_npw", l_npw
!  print *, "l_nspinor", l_nspinor
!  print *, "l_icplx*l_npw*l_nspinor", l_icplx*l_npw*l_nspinor
!  print *, "nband", nband
  !stop
  if (l_useria == 121212) then
    !call random_number(cg) !to test validity of algo
    !print *, "USAO RANDOM"
  end if
  
  call xgBlock_map(xgx0,cg,space,l_icplx*l_npw*l_nspinor,nband,l_mpi_enreg%comm_bandspinorfft) 
  
  !print *, "l_icplx*l_npw*l_nspinor", l_icplx*l_npw*l_nspinor
  !stop
  !print *, "l_npw", l_npw
  !print *, "nband", nband
  !stop
    
  if ( l_istwf == 2 ) then ! Real only
    ! Scale cg
    call xgBlock_scale(xgx0,sqrt2,1)  !ALL MPI processes do this
   
    ! This is possible since the memory in cg and xgx0 is the same
    ! Don't know yet how to deal with this with xgBlock
    !print *, "npw", npw
    !print *, "nspinor", nspinor
    !print *, "l_mpi_enreg%me_g0", l_mpi_enreg%me_g0
    !stop
    !MPI HANDLES THIS AUTOMATICALLY (only proc 0 is me_g0)
    if(l_mpi_enreg%me_g0 == 1) cg(:, 1:npw*nspinor*nband:npw) = cg(:, 1:npw*nspinor*nband:npw) * inv_sqrt2 
       
  end if
  
  !if (xmpi_comm_rank(l_mpi_enreg%comm_bandspinorfft) == 1) then   
    !print *, "EIG", eig
    !print *, "nband", nband
    !print *, "npw", npw
  !end if
  !stop
  !call xmpi_barrier(l_mpi_enreg%comm_bandspinorfft)
  
  !stop
  !call xg_init(xgeigen,SPACE_R,nband,1,l_mpi_enreg%comm_bandspinorfft)
  ! Trick with C is to change rank of arrays (:) to (:,:)
    
  cptr = c_loc(eig)			
  call c_f_pointer(cptr,eig_ptr,(/ nband,1 /))
  call xgBlock_map(xgeigen,eig_ptr,SPACE_R,nband,1,l_mpi_enreg%comm_bandspinorfft)
  !if (xmpi_comm_rank(l_mpi_enreg%comm_bandspinorfft) == 1) then  
    !print *, "EIG FROM OUT", l_nband_filter
    !print *, eig
    !print *, "XGEIGEN print"
    !call xgBlock_print(xgeigen, 200+xmpi_comm_rank(xmpi_world))
  !end if
  !print *, "nband", nband
  !stop
 
!  if (xmpi_comm_rank(l_mpi_enreg%comm_bandspinorfft) == 0) then  
!      !print *, "BEGIN" 
!      stop
!      print *, "eig", eig
!      !print *, "eig_ptr", eig_ptr
!      !call xgBlock_print(xgeigen, 6)
!      print *, "rank", xmpi_comm_rank(l_mpi_enreg%comm_bandspinorfft)
!      print *, "nband", nband
!      call xgBlock_print(xgeigen, 6)
!  end if
  
  !call xmpi_barrier(l_mpi_enreg%comm_bandspinorfft)
  
  !stop
  
  !call xg_init(xgresidu,SPACE_R,nband,1,l_mpi_enreg%comm_bandspinorfft)
  ! Trick the with C to change rank of arrays (:) to (:,:)
  cptr = c_loc(resid)		
  call c_f_pointer(cptr,resid_ptr,(/ nband,1 /))
  call xgBlock_map(xgresidu,resid_ptr,SPACE_R,nband,1,l_mpi_enreg%comm_bandspinorfft)
 
  ABI_MALLOC(l_gvnlc,(2,l_npw*l_nspinor*l_nband_filter)) !*blockdim
 
  !print *, "mpi_enreg%nproc_band", mpi_enreg%nproc_band
  !print *, "xmpi_comm_size(xmpi_world)", xmpi_comm_size(xmpi_world)
  !print *, "l_mpi_enreg%bandpp", l_mpi_enreg%bandpp
  !stop
  
  !print *, "l_icplx*l_npw*l_nspinor", l_icplx*l_npw*l_nspinor
  !print *, "l_mpi_enreg%nproc_band", l_mpi_enreg%nproc_band
  !stop
  !print *, "INIT SPACE", space
  !stop


  call chebfi_init(chebfi, nband, l_icplx*l_npw*l_nspinor, dtset%tolwfr, dtset%ecut, dtset%paral_kgb, l_mpi_enreg%nproc_band, &
                   l_mpi_enreg%bandpp, l_mpi_enreg%nproc_fft, nline, space, 1, l_gs_hamk%istwf_k, l_mpi_enreg%comm_bandspinorfft, l_mpi_enreg%me_g0, l_paw)
  !chebfi_init(chebfi, neigenpairs, spacedim, tolerance, ecut, nline, space, eigenProblem, istwf_k, spacecom, me_g0, paw)                
  !call chebfi_init(chebfi, nband, l_icplx*l_npw*l_nspinor, dtset%tolwfr, dtset%ecut, nline, space, 1, l_gs_hamk%istwf_k, l_mpi_enreg%comm_bandspinorfft, l_mpi_enreg%me_g0, l_paw) 
 
  !###########################################################################
  !################    RUUUUUUUN    ##########################################
  !###########################################################################
  !print *, "BEFORE"
  !stop
  !call xg_getPointer(xgx0)   
  !stop
  ! Run chebfi
  print *, "CB2 RUN STARTED!"
  !stop
  print *, "mpi_enreg%comm_fft BEFORE CB2, rank", l_mpi_enreg%comm_fft, xmpi_comm_rank(l_mpi_enreg%comm_bandspinorfft)
  print *, "mpi_enreg%comm_band BEFORE CB2, rank", l_mpi_enreg%comm_band, xmpi_comm_rank(l_mpi_enreg%comm_bandspinorfft)
  !print *, "mpi_enreg%comm_kpt BEFORE CB2", l_mpi_enreg%comm_kpt
  !print *, "mpi_enreg%comm_spinor BEFORE CB2", l_mpi_enreg%comm_spinor 
  !print *, "mpi_enreg%comm_bandspinor BEFORE CB2", l_mpi_enreg%comm_bandspinor 
  !print *, "mpi_enreg%comm_kptband BEFORE CB2", l_mpi_enreg%comm_kptband 
  !print *, "mpi_enreg%comm_spinorfft BEFORE CB2", l_mpi_enreg%comm_spinorfft 
  !print *, "mpi_enreg%comm_bandfft BEFORE CB2", l_mpi_enreg%comm_bandfft 
  !print *, "mpi_enreg%comm_bandspinorfft BEFORE CB2", l_mpi_enreg%comm_bandspinorfft 
  
  call chebfi_run(chebfi, xgx0, getghc_gsc1, getBm1X, precond1, xgeigen, xgresidu, counter, l_mpi_enreg) 
  !stop  

  print *, "mpi_enreg%comm_fft AFTER CB2, rank", l_mpi_enreg%comm_fft, xmpi_comm_rank(l_mpi_enreg%comm_bandspinorfft)
  print *, "mpi_enreg%comm_band AFTER CB2, rank", l_mpi_enreg%comm_band, xmpi_comm_rank(l_mpi_enreg%comm_bandspinorfft)
  !print *, "mpi_enreg%comm_kpt AFTER CB2", l_mpi_enreg%comm_kpt
  !print *, "mpi_enreg%comm_spinor AFTER CB2", l_mpi_enreg%comm_spinor 
  !print *, "mpi_enreg%comm_bandspinor AFTER CB2", l_mpi_enreg%comm_bandspinor 
  !print *, "mpi_enreg%comm_kptband AFTER CB2", l_mpi_enreg%comm_kptband 
  !print *, "mpi_enreg%comm_spinorfft AFTER CB2", l_mpi_enreg%comm_spinorfft 
  !print *, "mpi_enreg%comm_bandfft AFTER CB2", l_mpi_enreg%comm_bandfft 
  !print *, "mpi_enreg%comm_bandspinorfft AFTER CB2", l_mpi_enreg%comm_bandspinorfft 
  
  !if (xmpi_comm_rank(l_mpi_enreg%comm_bandspinorfft) == 0) then
  !print *, "CB2 RUN FINISHED!"
  !stop
  print *, "CB2WF PREPNONLOP pi_enreg%comm_fft, rank", mpi_enreg%comm_fft, xmpi_comm_rank(l_mpi_enreg%comm_bandspinorfft)
  print *, "CB2WF PREPNONLOP pi_enreg%comm_band, rank", mpi_enreg%comm_band, xmpi_comm_rank(l_mpi_enreg%comm_bandspinorfft)
  !stop
    !call xg_getPointer(xgx0) 
  !end if  
  !stop
  
  !call debug_helper_linalg(xgx0, l_mpi_enreg%comm_bandspinorfft, l_icplx*l_npw*l_nspinor, nband, 1, 1) !MPI 1,2 OK
  !stop
  
  !call xgBlock_setBlock(xgx0, HELPER, 1, 2, 5) 
  !call xgBlock_print(HELPER, 201) 
  !stop
  !stop
  ! Free preconditionning since not needed anymore
  ABI_FREE(l_pcon)
    
  ! Scale back
!  if(l_istwf == 2) then
!    call xgBlock_scale(xgx0,inv_sqrt2,1)
!    if(l_mpi_enreg%me_g0 == 1) cg(:, 1:npw*nspinor*nband:npw) = cg(:, 1:npw*nspinor*nband:npw) * sqrt2 !not sure what this does exactly
!  end if
  
  print *, "BEFORE CBWF2 NONLOP"
  !stop  
  ! Compute enlout (nonlocal energy for each band if necessary) This is the best
  ! quick and dirty trick to compute this part in NC. gvnlc cannot be part of
  ! chebfi algorithm 
  if ( .not. l_paw ) then
    !Check l_gvnlc size
    !if ( size(l_gvnlc) < 2*nband*l_npw*l_nspinor ) then
    if ( size(l_gvnlc) /= 0 ) then
      ABI_FREE(l_gvnlc)
      !ABI_MALLOC(l_gvnlc,(2,nband*l_npw*l_nspinor))
      ABI_MALLOC(l_gvnlc,(0,0))
    end if
    !Call nonlop
    if (mpi_enreg%paral_kgb==0) then
      print *, "USAO"
      call nonlop(choice,l_cpopt,cprj_dum,enl_out,l_gs_hamk,0,eig,mpi_enreg,nband,1,paw_opt,&
&                signs,gsc_dummy,l_tim_getghc,cg,l_gvnlc)

      print *, "IZASAO"
      stop
    else
      !do iband=1,nband
      !  shift = (iband-1)*l_npw*l_nspinor
      !  print *, "IBAND", iband
      !  call prep_nonlop(choice,l_cpopt,cprj_dum, &
      !&       enl_out((iband-1)+1:iband),l_gs_hamk,0,&
      !&       eig((iband-1)+1:iband),1,mpi_enreg,1,paw_opt,signs,& !1=blockdim?
      !&       gsc_dummy,l_tim_getghc, &
      !&       cg(:,shift+1:shift+l_npw*l_nspinor),&
      !&       l_gvnlc(:,shift+1:shift*l_npw*l_nspinor),&
      !!&       l_gvnlc(:,:),&
      !&       already_transposed=.false.)
      !end do
     do iband=1,nband/blockdim
       shift = (iband-1)*blockdim*l_npw*l_nspinor
       call prep_nonlop(choice,l_cpopt,cprj_dum, &
&       enl_out((iband-1)*blockdim+1:iband*blockdim),l_gs_hamk,0,&
&       eig((iband-1)*blockdim+1:iband*blockdim),blockdim,mpi_enreg,1,paw_opt,signs,&
&       gsc_dummy,l_tim_getghc, &
&       cg(:,shift+1:shift+blockdim*l_npw*l_nspinor),&
!&      l_gvnlc(:,shift+1:shift+blockdim*l_npw*l_nspinor),&
&       l_gvnlc(:,:),&
&       already_transposed=.false.)
      end do
    end if
    !Compute enlout
!   do iband=1,nband
!     shift = npw*nspinor*(iband-1)
!       call dotprod_g(enl_out(iband),dprod_i,l_gs_hamk%istwf_k,npw*nspinor,1,cg(:, shift+1:shift+npw*nspinor),&
!  &     l_gvnlc(:, shift+1:shift+npw*nspinor),mpi_enreg%me_g0,mpi_enreg%comm_bandspinorfft)
!   end do
  end if

  print *, "AFTER CBWF2 NONLOP"
  stop
  ABI_FREE(l_gvnlc)

  ! Free chebfi
  call chebfi_free(chebfi)
  
  !###########################################################################
  !################    SORRY IT'S ALREADY FINISHED : )  ######################
  !###########################################################################
 
  call timab(tim_chebfiwf2,2,tsec)
  cputime = abi_cpu_time() - cputime
  walltime = abi_wtime() - walltime
  nthreads = xomp_get_num_threads(open_parallel = .true.)
  if ( cputime/walltime/dble(nthreads) < 0.75 .and. (int(cputime/walltime)+1) /= nthreads) then
    if ( prtvol >= 3 ) then
      write(std_out,'(a)',advance='no') sjoin(" Chebfi took", sec2str(cputime), "of cpu time")
      write(std_out,*) sjoin("for a wall time of", sec2str(walltime))
      write(std_out,'(a,f6.2)') " -> Ratio of ", cputime/walltime
    end if
    MSG_COMMENT(sjoin("You should set the number of threads to something close to",itoa(int(cputime/walltime)+1)))
  end if
 
  print *, "END OF CBWF2"
  !stop
  DBG_EXIT("COLL")
 
  end subroutine chebfiwf2
 
  subroutine getghc_gsc1(X,AX,BX,transposer)	
    use m_xg, only : xg_t, xgBlock_get, xgBlock_set, xgBlock_getSize, xgBlock_t
    use m_xgTransposer !, only xgTransposer, xgTransposer_getCPURow
#ifdef HAVE_OPENMP
    use omp_lib
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getghc_gsc1'
!End of the abilint section

    type(xgBlock_t), intent(inout) :: X
    type(xgBlock_t), intent(inout) :: AX
    type(xgBlock_t), intent(inout) :: BX
    type(xgTransposer_t), optional, intent(inout) :: transposer
    integer         :: blockdim
    integer         :: spacedim
    type(pawcprj_type) :: cprj_dum(l_gs_hamk%natom,0) 
    !l_nspinor*l_nband_filter
    double precision :: eval
    double precision, parameter :: inv_sqrt2 = 1/sqrt2
    double precision, pointer :: cg(:,:)
    double precision, pointer :: ghc(:,:)
    double precision, pointer :: gsc(:,:)
    
    integer :: cpuRow
    
    call xgBlock_getSize(X,spacedim,blockdim)

    !print *, "spacedim X", spacedim
    !print *, "blockdim X", blockdim
    
    !call xgBlock_getSize(AX,spacedim,blockdim)
    
    !print *, "spacedim AX", spacedim
    !print *, "blockdim AX", blockdim
    
    !stop
    
    !call xgBlock_getSize(BX,spacedim,blockdim)
    
   ! print *, "spacedim BX", spacedim
    !print *, "blockdim BX", blockdim
    
    !stop

    spacedim = spacedim/l_icplx
    
    !call xgBlock_get(X,cg(:,1:blockdim*spacedim),0,spacedim)
    call xgBlock_reverseMap(X,cg,l_icplx,spacedim*blockdim)
    call xgBlock_reverseMap(AX,ghc,l_icplx,spacedim*blockdim)
    call xgBlock_reverseMap(BX,gsc,l_icplx,spacedim*blockdim)
    
    !call debug_helper_colrows(X, blockdim)
    !stop
    
    !print *, "spacedim", spacedim
    !print *, "blockdim", blockdim
    !print *, "l_npw", l_npw
    !stop
    
    ! scale back cg
    if(l_istwf == 2) then
      !cg(:,1:spacedim*blockdim) = cg(:,1:spacedim*blockdim) * inv_sqrt2
      call xgBlock_scale(X,inv_sqrt2,1)
      !suppress scaling to debug MPI
      if (l_paral_kgb == 0) then
        if(l_mpi_enreg%me_g0 == 1) cg(:, 1:spacedim*blockdim:l_npw) = cg(:, 1:spacedim*blockdim:l_npw) * sqrt2
      else
        !call xgTransposer_getCPURow(transposer, cpuRow)
        cpuRow = xgTransposer_getRank(transposer, 2)
        !print *, "CPUROW", cpurow
        !l_npw is not same for all processes so the lower formula cannot be like
        !that
        if (cpuRow == 0) then
          cg(:, 1:spacedim*blockdim:spacedim) = cg(:, 1:spacedim*blockdim:spacedim) * sqrt2
        end if
      end if
    end if
      
    !call debug_helper_colrows(X, blockdim)
    !stop 

    if ( size(l_gvnlc) < 2*blockdim*spacedim ) then
      ABI_FREE(l_gvnlc)
      ABI_MALLOC(l_gvnlc,(2,blockdim*spacedim))
    end if
      
    !print *, "before getghc"  
      
    if (l_useria == 121212) then
      !call random_number(ghc)
      ghc = 1
      gsc = 2
      !call random_number(gsc)
      !call xgBlock_copy(X,AX)
      !call xgBlock_copy(X,BX)
    else
      if (l_mpi_enreg%paral_kgb==0) then
        !l_cpopt = -1
        call multithreaded_getghc(l_cpopt,cg,cprj_dum,ghc,gsc,& 
          l_gs_hamk,l_gvnlc,eval,l_mpi_enreg,blockdim,l_prtvol,l_sij_opt,l_tim_getghc,0) 
      else
        !print *, "blockdim", blockdim
        !print *, "spacedim", spacedim
        !print *, "before prep_getghc"
        print *, "l_mpi_enreg%comm_fft INSIDE getAX_BX", l_mpi_enreg%comm_fft
        print *, "l_mpi_enreg%comm_band INSIDE getAX_BX", l_mpi_enreg%comm_band
        call prep_getghc(cg(:,1:blockdim*spacedim),l_gs_hamk,l_gvnlc,ghc,gsc(:,1:blockdim*spacedim),eval,blockdim,l_mpi_enreg,&
&                       l_prtvol,l_sij_opt,l_cpopt,cprj_dum,already_transposed=.true.)  !already_transposed = true (previous)
        !print *, "after prep_getghc"
      end if
    end if
    
    !print *, "after getghc"  
    
    !call debug_helper_colrows(AX, blockdim)
    !stop
  
    ! scale cg, ghc, gsc
    if ( l_istwf == 2 ) then
      !cg(:,1:spacedim*blockdim) = cg(:,1:spacedim*blockdim) * sqrt2
      !ghc(:,1:spacedim*blockdim) = ghc(:,1:spacedim*blockdim) * sqrt2
      call xgBlock_scale(X,sqrt2,1)
      call xgBlock_scale(AX,sqrt2,1)
      !suppress scaling to debug MPI
      if (l_paral_kgb == 0) then
        if(l_mpi_enreg%me_g0 == 1) then
          !print *, "RANK", xmpi_comm_rank(l_mpi_enreg%comm_bandspinorfft)
          cg(:, 1:spacedim*blockdim:l_npw) = cg(:, 1:spacedim*blockdim:l_npw) * inv_sqrt2
          ghc(:, 1:spacedim*blockdim:l_npw) = ghc(:, 1:spacedim*blockdim:l_npw) * inv_sqrt2
        endif
      else 
        if (cpuRow == 0) then
          !print *, "RANK", xmpi_comm_rank(l_mpi_enreg%comm_bandspinorfft)
          cg(:, 1:spacedim*blockdim:spacedim) = cg(:, 1:spacedim*blockdim:spacedim) * inv_sqrt2
          ghc(:, 1:spacedim*blockdim:spacedim) = ghc(:, 1:spacedim*blockdim:spacedim) * inv_sqrt2
        end if
      end if
      if(l_paw) then
        !gsc(:,1:spacedim*blockdim) = gsc(:,1:spacedim*blockdim) * sqrt2
        call xgBlock_scale(BX,sqrt2,1)
        if (l_paral_kgb == 0) then
          if(l_mpi_enreg%me_g0 == 1) gsc(:, 1:spacedim*blockdim:l_npw) = gsc(:, 1:spacedim*blockdim:l_npw) * inv_sqrt2
        else
          if (cpuRow == 0) then
            !print *, "USAO CPUROW0"
            !stop
            gsc(:, 1:spacedim*blockdim:spacedim) = gsc(:, 1:spacedim*blockdim:spacedim) * inv_sqrt2
          end if
        end if
      end if
    end if
    
    !print *, "BAJA DO JAJA"
    !call debug_helper_colrows(AX, blockdim)
    !stop
    
    if ( .not. l_paw ) call xgBlock_copy(X,BX)
  
    !call xgBlock_set(AX,ghc,0,spacedim)
    !call xgBlock_set(BX,gsc(:,1:blockdim*spacedim),0,spacedim)
  end subroutine getghc_gsc1
 
  subroutine getBm1X(X,Bm1X,iline_t,xXColsRows,X1,transposer)
    use m_xg, only : xgBlock_t, xgBlock_gemm
    use m_xgTransposer !, only xgTransposer, xgTransposer_getCPURow
!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getBm1X'
!End of the abilint section

    type(xgBlock_t), intent(inout) :: X
    type(xgBlock_t), intent(inout) :: Bm1X
    type(integer), intent(inout) :: iline_t
    type(xgBlock_t), intent(inout) :: xXColsRows
    type(xgBlock_t), intent(inout) :: X1
    type(xgTransposer_t), optional, intent(inout) :: transposer
    integer         :: blockdim
    integer         :: spacedim   
    
    integer :: cpuRow

    double precision, pointer :: ghc_filter(:,:)
    double precision, pointer :: gsm1hc_filter(:,:)
    type(pawcprj_type), allocatable :: cwaveprj_next(:,:) !dummy
    
    call xgBlock_getSize(X,spacedim,blockdim)
    !call xgBlock_getSize(Bm1X,spacedim,blockdim)
    spacedim = spacedim/l_icplx
    
    !print *, "spacedim", spacedim
    !print *, "BLOCKDIM", blockdim
    !print *, "l_icplx", l_icplx
    
    !print *, "spacedim*blockdim 2", spacedim*blockdim
    !stop
!    if (iline_t == 2) then
!      print *, "USAO U getBm1X"
!      call xgBlock_setBlock(Bm1X, xXColsRows, 1, 3888, 96)  
!      call xgTransposer_transpose(transposer,STATE_LINALG) !all_to_all
!      !call xmpi_barrier(chebfi%spacecom)
!      call debug_helper_linalg(X1)
!      stop
!    end if
    
    !stop
    call xgBlock_reverseMap(X,ghc_filter,l_icplx,spacedim*blockdim)
    !stop
    call xgBlock_reverseMap(Bm1X,gsm1hc_filter,l_icplx,spacedim*blockdim)
    
    if (l_useria == 121212) then
      !call random_number(ghc_filter)
      !call random_number(gsm1hc_filter)
      gsm1hc_filter = 3
      !call xgBlock_copy(X,Bm1X)
    end if
    !stop
!        print *, "spacedim", spacedim
!          print *, "blockdim", blockdim
!          print *, "l_npw", l_npw
!          stop     
    ! scale back cg
    if(l_istwf == 2) then 
      !cg(:,1:spacedim*blockdim) = cg(:,1:spacedim*blockdim) * inv_sqrt2
      call xgBlock_scale(X,inv_sqrt2,1)
      if (l_paral_kgb == 0) then
        if(l_mpi_enreg%me_g0 == 1) ghc_filter(:, 1:spacedim*blockdim:l_npw) = ghc_filter(:, 1:spacedim*blockdim:l_npw) * sqrt2
      else
        !call xgTransposer_getCPURow(transposer, cpuRow)
        cpuRow = xgTransposer_getRank(transposer, 2)
        if (cpuRow == 0) then
          ghc_filter(:, 1:spacedim*blockdim:spacedim) = ghc_filter(:, 1:spacedim*blockdim:spacedim) * sqrt2
        end if
      end if
      if(l_paw) then
        if (iline_t == 2) then
          !gsc(:,1:spacedim*blockdim) = gsc(:,1:spacedim*blockdim) * sqrt2
          !print *, "PRINTING BM1X"
          !call debug_helper_colrows(Bm1X)
          !stop
!          call xgBlock_setBlock(Bm1X, xXColsRows, 1, 3888, 96)  
!          call xgTransposer_transpose(transposer,STATE_LINALG) !all_to_all
!          call debug_helper_linalg(X1)
!          stop
        end if
        call xgBlock_scale(Bm1X,inv_sqrt2,1)
        if (l_paral_kgb == 0) then
          if(l_mpi_enreg%me_g0 == 1) gsm1hc_filter(:, 1:spacedim*blockdim:l_npw) = gsm1hc_filter(:, 1:spacedim*blockdim:l_npw) * sqrt2 
        else
          if (cpuRow == 0) then
            print *, "CPUROW0"
            !!TODO OVDE JE PROBLEM NE ZNAM ZASTO
            gsm1hc_filter(:, 1:spacedim*blockdim:spacedim) = gsm1hc_filter(:, 1:spacedim*blockdim:spacedim) * sqrt2
          end if
         if (iline_t == 2) then
            !print *, "spacedim", spacedim
            !print *, "blockdim", blockdim
            !print *, "l_npw", l_npw
            !stop
            !print *, "USAO U getBm1X"
            !call xgBlock_setBlock(Bm1X, xXColsRows, 1, 3888, 96)  
            !call xgTransposer_transpose(transposer,STATE_LINALG) !all_to_all
            !call debug_helper_linalg(X1)
            !stop
          end if
        end if
      end if
    end if
    
!    if (iline_t == 2) then
!      print *, "USAO U getBm1X"
!      call xgBlock_setBlock(Bm1X, xXColsRows, 1, 3888, 96)  
!      call xgTransposer_transpose(transposer,STATE_LINALG) !all_to_all
!      !call xmpi_barrier(chebfi%spacecom)
!      call debug_helper_linalg(X1)
!      stop
!    end if
    
    !print *, "PROSAO 1"
    !stop
    !cwaveprj dummy allocate

    if(l_paw) then
      !print *, "blockdim", blockdim
      !stop
      if (l_useria /= 121212) then
        ABI_DATATYPE_ALLOCATE(cwaveprj_next, (l_gs_hamk%natom,l_nspinor*blockdim)) !nband_filter
        !print *, "ABI_DATATYPE_ALLOCATE PASS"
        !stop
        !print *, "l_gs_hamk%dimcprj", l_gs_hamk%dimcprj
        call pawcprj_alloc(cwaveprj_next,0,l_gs_hamk%dimcprj) 
        !print *, "DUMMY ALLOCATE PASS"
        !stop
        call apply_invovl(l_gs_hamk, ghc_filter(:,:), gsm1hc_filter(:,:), cwaveprj_next(:,:), & 
&       spacedim, blockdim, l_mpi_enreg, l_nspinor)   !!ghc_filter size problem if transposed
        !print *, "APPLY INVOVL PASS"
        !stop
      end if
    else
      gsm1hc_filter(:,:) = ghc_filter(:,:)
    end if
    
    !print *, "OVDEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
    ! scale cg, ghc, gsc
    if ( l_istwf == 2 ) then
      !cg(:,1:spacedim*blockdim) = cg(:,precond1:spacedim*blockdim) * sqrt2
      !ghc(:,1:spacedim*blockdim) = ghc(:,1:spacedim*blockdim) * sqrt2
      call xgBlock_scale(X,sqrt2,1)
      if (l_paral_kgb == 0) then
        if(l_mpi_enreg%me_g0 == 1) then
          ghc_filter(:, 1:spacedim*blockdim:l_npw) = ghc_filter(:, 1:spacedim*blockdim:l_npw) * inv_sqrt2
        endif
      else
        if (cpuRow == 0) then
          ghc_filter(:, 1:spacedim*blockdim:spacedim) = ghc_filter(:, 1:spacedim*blockdim:spacedim) * inv_sqrt2
        end if
      end if
      if(l_paw) then
        !gsc(:,1:spacedim*blockdim) = gsc(:,1:spacedim*blockdim) * sqrt2
        call xgBlock_scale(Bm1X,sqrt2,1)
        if (l_paral_kgb == 0) then
          if(l_mpi_enreg%me_g0 == 1) gsm1hc_filter(:, 1:spacedim*blockdim:l_npw) = gsm1hc_filter(:, 1:spacedim*blockdim:l_npw) * inv_sqrt2
        else
          if (cpuRow == 0) then
            gsm1hc_filter(:, 1:spacedim*blockdim:spacedim) = gsm1hc_filter(:, 1:spacedim*blockdim:spacedim) * inv_sqrt2
          end if
        end if
      end if
    end if

    if (l_paw) then
      if (l_useria /= 121212) then
        call pawcprj_free(cwaveprj_next)
        ABI_DATATYPE_DEALLOCATE(cwaveprj_next)
      end if
    end if
    
    !print *, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    
  end subroutine getBm1X
 
  subroutine precond1(W)
    use m_xg, only : xg_t, xgBlock_colwiseMul

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'precond1'
!End of the abilint section

    type(xgBlock_t), intent(inout) :: W
    integer :: ispinor
    !integer :: cplx
    
    ! precondition resid_vec
    do ispinor = 1,l_nspinor
      !do cplx = 1, l_icplx
        call xgBlock_colwiseMul(W,l_pcon,l_icplx*l_npw*(ispinor-1))
      !end do
    end do

  end subroutine precond1
  
!  subroutine debug_me_g0_COLROWS(debugBlock, chebfi)
!  
!    type(xgBlock_t) , intent(inout) :: debugBlock
!    type(chebfi_t) , intent(inout) :: chebfi
!    type(xgBlock_t) :: HELPER
!    
!    call xmpi_barrier(chebfi%spacecom)
!    
!    if (xmpi_comm_size(xmpi_world) == 1) then !only one MPI proc
!      if (xmpi_comm_rank(chebfi%spacecom) == 0) then
!        call xgBlock_setBlock(debugBlock, HELPER, 1, 2, 5) 
!        call xgBlock_print(HELPER, 100) 
!      
!        call xgBlock_setBlock(debugBlock, HELPER, chebfi%bandpp/2+1, 2, 5) 
!        call xgBlock_print(HELPER, 101) 
!      end if
!    else 
!      call xgBlock_setBlock(debugBlock, HELPER, 1, 2, 5) 
!      call xgBlock_print(HELPER, 200 + xmpi_comm_rank(chebfi%spacecom)) 
!    end if
!    
!    call xmpi_barrier(chebfi%spacecom)   
!    
!  end subroutine debug_me_g0_COLROWS
  
!  subroutine debug_helper_colrows(debugBlock, chebfi)
!      
!    type(xgBlock_t) , intent(inout) :: debugBlock
!    type(chebfi_t) , intent(inout) :: chebfi
!    type(xgBlock_t) :: HELPER
!    
!    call xmpi_barrier(chebfi%spacecom)
!    
!    if (xmpi_comm_size(xmpi_world) == 1) then !only one MPI proc
!      !print *, "USAOOSOOSOSOSO"
!      call xgBlock_setBlock(debugBlock, HELPER, 1, DEBUG_ROWS, DEBUG_COLUMNS) 
!      call xgBlock_print(HELPER, 100+xmpi_comm_rank(chebfi%spacecom)) 
!      
!      call xgBlock_setBlock(debugBlock, HELPER, chebfi%bandpp/2+1, DEBUG_ROWS, DEBUG_COLUMNS) 
!      call xgBlock_print(HELPER, 100+xmpi_comm_rank(chebfi%spacecom)+1) 
!    else
!      !print *, "2222222222222222222"
!      call xgBlock_setBlock(debugBlock, HELPER, 1, DEBUG_ROWS, DEBUG_COLUMNS) 
!      call xgBlock_print(HELPER, 200+xmpi_comm_rank(chebfi%spacecom)) 
!    end if
!    
!    !print *, "debugBlock%rows", rows(debugBlock)
!    !print *, "debugBlock%cols", cols(debugBlock)
!    call xmpi_barrier(chebfi%spacecom)

!  end subroutine debug_helper_colrows
!  
!  subroutine debug_helper_linalg(debugBlock)
!      
!    type(xgBlock_t) , intent(inout) :: debugBlock
!    type(xgBlock_t) :: HELPER
!    
!    integer, parameter :: FROW = 1, FCOL = 1, DROWS = 1, DCOLS = 10
!     
!      if (xmpi_comm_size(xmpi_world) == 1) then !only one MPI proc
!        call xgBlock_setBlock1(debugBlock, HELPER, 1, FCOL, DROWS, DCOLS) 
!        call xgBlock_print(HELPER, 100) 
!        
!        call xgBlock_setBlock1(debugBlock, HELPER, 3888/2+1, FCOL, DROWS, DCOLS) 
!        call xgBlock_print(HELPER, 101)    
!      else
!        if (xmpi_comm_rank(xmpi_world) == 0) then
!          call xgBlock_setBlock1(debugBlock, HELPER, 1, FCOL, DROWS, DCOLS) 
!          call xgBlock_print(HELPER, 200)
!        else 
!          call xgBlock_setBlock1(debugBlock, HELPER, 1, FCOL, DROWS, DCOLS) 
!          call xgBlock_print(HELPER, 201)
!        end if
!      end if

!  end subroutine debug_helper_linalg
  
end module m_chebfiwf
