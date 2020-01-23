!!****f* ABINIT/m_lobpcgwf
!! NAME
!! m_lobpcgwf
!!
!! FUNCTION
!! this routine updates the whole wave functions at a given k-point,
!! using the lobpcg method
!! for a given spin-polarization, from a fixed hamiltonian
!! but might also simply compute eigenvectors and eigenvalues at this k point.
!! it will also update the matrix elements of the hamiltonian.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group (JB)
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

module m_lobpcgwf

 use defs_basis
 use m_abicore
 use m_lobpcg
 use m_xmpi
 use m_errors
 use m_time
 use m_xomp
 use m_fstrings
 use m_xg
 use m_lobpcg2
 use m_dtset

 use defs_abitypes, only : mpi_type
 use m_hamiltonian, only : gs_hamiltonian_type
 use m_pawcprj,     only : pawcprj_type
 use m_nonlop,      only : nonlop
 use m_prep_kgb,    only : prep_getghc, prep_nonlop
 use m_getghc,      only : multithreaded_getghc

 private

 integer, parameter :: l_tim_getghc=5
 double precision, parameter :: inv_sqrt2 = 1/sqrt2
 ! Use in getghc_gsc
 integer,save  :: l_cpopt
 integer,save  :: l_icplx
 integer,save  :: l_istwf
 integer,save  :: l_npw
 integer,save  :: l_nspinor
 logical,save  :: l_paw
 integer,save  :: l_prtvol
 integer,save  :: l_sij_opt
 real(dp), allocatable,save :: l_gvnlxc(:,:)
 real(dp), allocatable,save ::  l_pcon(:)
 type(mpi_type),pointer,save :: l_mpi_enreg
 type(gs_hamiltonian_type),pointer,save :: l_gs_hamk

 public :: lobpcgwf2

 contains

subroutine lobpcgwf2(cg,dtset,eig,enl_out,gs_hamk,kinpw,mpi_enreg,&
&                   nband,npw,nspinor,prtvol,resid)


 use m_cgtools, only : dotprod_g
 use iso_c_binding
 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: nband,npw,prtvol,nspinor
 type(gs_hamiltonian_type),target,intent(inout) :: gs_hamk
 type(dataset_type)              ,intent(in   ) :: dtset
 type(mpi_type)           ,target,intent(inout) :: mpi_enreg
 real(dp)                 ,target,intent(inout) :: cg(2,nspinor*nband*npw)!,gsc(2,nspinor*nband*npw)
 real(dp)                        ,intent(in   ) :: kinpw(npw)
 real(dp)                 ,target,intent(  out) :: resid(nband)
 real(dp)                        ,intent(  out) :: enl_out(nband)
 real(dp)                 ,target,intent(  out) :: eig(nband)

!Local variables-------------------------------

 type(xgBlock_t) :: xgx0
 type(xgBlock_t) :: xgeigen
 type(xgBlock_t) :: xgresidu
 type(lobpcg_t) :: lobpcg

 integer :: space, blockdim,  nline
 integer :: ipw

 integer :: nthreads

 double precision :: lobpcgMem(2)
 double precision :: localmem

 integer, parameter :: tim_lobpcgwf2 = 1650
 double precision :: cputime, walltime
 double precision :: tsec(2)

 type(c_ptr) :: cptr
 real(dp), pointer :: eig_ptr(:,:) => NULL()
 real(dp), pointer :: resid_ptr(:,:) => NULL()

 ! Important things for NC
 integer,parameter :: choice=1, paw_opt=0, signs=1
 type(pawcprj_type) :: cprj_dum(gs_hamk%natom,0)
 integer :: iband, shift
 real(dp) :: gsc_dummy(0,0)

! *********************************************************************


!###########################################################################
!################ INITIALISATION  ##########################################
!###########################################################################

 call timab(tim_lobpcgwf2,1,tsec)
 cputime = abi_cpu_time()
 walltime = abi_wtime()

 ! Set module variables
 l_paw = (gs_hamk%usepaw==1)
 l_cpopt=-1;l_sij_opt=0;if (l_paw) l_sij_opt=1
 l_istwf=gs_hamk%istwf_k
 l_npw = npw
 l_nspinor = nspinor
 l_prtvol = prtvol
 l_mpi_enreg => mpi_enreg
 l_gs_hamk => gs_hamk

!Variables
 nline=dtset%nline
 blockdim=l_mpi_enreg%nproc_band*l_mpi_enreg%bandpp

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

 ! Memory info
 if ( prtvol >= 3 ) then
   lobpcgMem = lobpcg_memInfo(nband,l_icplx*l_npw*l_nspinor,blockdim,space)
   localMem = (l_npw+2*l_npw*l_nspinor*blockdim+2*nband)*kind(1.d0)
   write(std_out,'(1x,A,F10.6,1x,A)') "Each MPI process calling lobpcg should need around ", &
   (localMem+sum(lobpcgMem))/1e9, &
   "GB of peak memory as follows :"
   write(std_out,'(4x,A,F10.6,1x,A)') "Permanent memory in lobpcgwf : ", &
   (localMem)/1e9, "GB"
   write(std_out,'(4x,A,F10.6,1x,A)') "Permanent memory in m_lobpcg : ", &
   (lobpcgMem(1))/1e9, "GB"
   write(std_out,'(4x,A,F10.6,1x,A)') "Temporary memory in m_lobpcg : ", &
   (lobpcgMem(2))/1e9, "GB"
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

 ! Local variables for lobpcg
 !call xg_init(xgx0,space,icplx*npw*nspinor,nband)
 call xgBlock_map(xgx0,cg,space,l_icplx*l_npw*l_nspinor,nband,l_mpi_enreg%comm_bandspinorfft)
 if ( l_istwf == 2 ) then ! Real only
   ! Scale cg
   call xgBlock_scale(xgx0,sqrt2,1)
   ! This is possible since the memory in cg and xgx0 is the same
   ! Don't know yet how to deal with this with xgBlock
   if(l_mpi_enreg%me_g0 == 1) cg(:, 1:npw*nspinor*nband:npw) = cg(:, 1:npw*nspinor*nband:npw) * inv_sqrt2
 end if

 !call xg_init(xgeigen,SPACE_R,nband,1,l_mpi_enreg%comm_bandspinorfft)
 ! Trick the with C to change rank of arrays (:) to (:,:)
 cptr = c_loc(eig)
 call c_f_pointer(cptr,eig_ptr,(/ nband,1 /))
 call xgBlock_map(xgeigen,eig_ptr,SPACE_R,nband,1)

 !call xg_init(xgresidu,SPACE_R,nband,1,l_mpi_enreg%comm_bandspinorfft)
 ! Trick the with C to change rank of arrays (:) to (:,:)
 cptr = c_loc(resid)
 call c_f_pointer(cptr,resid_ptr,(/ nband,1 /))
 call xgBlock_map(xgresidu,resid_ptr,SPACE_R,nband,1)

 ABI_MALLOC(l_gvnlxc,(2,l_npw*l_nspinor*blockdim))

 call lobpcg_init(lobpcg,nband, l_icplx*l_npw*l_nspinor, blockdim,dtset%tolwfr,nline,space, l_mpi_enreg%comm_bandspinorfft)

!###########################################################################
!################    RUUUUUUUN    ##########################################
!###########################################################################

 ! Run lobpcg
 call lobpcg_run(lobpcg,xgx0,getghc_gsc,precond,xgeigen,xgresidu,prtvol)

 ! Free preconditionning since not needed anymore
 ABI_FREE(l_pcon)

 ! Scale back
 if(l_istwf == 2) then
   call xgBlock_scale(xgx0,inv_sqrt2,1)
   if(l_mpi_enreg%me_g0 == 1) cg(:, 1:npw*nspinor*nband:npw) = cg(:, 1:npw*nspinor*nband:npw) * sqrt2
 end if

 ! Compute enlout (nonlocal energy for each band if necessary) This is the best
 ! quick and dirty trick to compute this part in NC. gvnlxc cannot be part of
 ! lobpcg algorithm
 if ( .not. l_paw ) then
   !Check l_gvnlxc size
   !if ( size(l_gvnlxc) < 2*nband*l_npw*l_nspinor ) then
   if ( size(l_gvnlxc) /= 0 ) then
     ABI_FREE(l_gvnlxc)
     !ABI_MALLOC(l_gvnlxc,(2,nband*l_npw*l_nspinor))
     ABI_MALLOC(l_gvnlxc,(0,0))
   end if
   !Call nonlop
   if (mpi_enreg%paral_kgb==0) then

     call nonlop(choice,l_cpopt,cprj_dum,enl_out,l_gs_hamk,0,eig,mpi_enreg,nband,1,paw_opt,&
&                signs,gsc_dummy,l_tim_getghc,cg,l_gvnlxc)

   else
     do iband=1,nband/blockdim
       shift = (iband-1)*blockdim*l_npw*l_nspinor
      call prep_nonlop(choice,l_cpopt,cprj_dum, &
&       enl_out((iband-1)*blockdim+1:iband*blockdim),l_gs_hamk,0,&
&       eig((iband-1)*blockdim+1:iband*blockdim),blockdim,mpi_enreg,1,paw_opt,signs,&
&       gsc_dummy,l_tim_getghc, &
&       cg(:,shift+1:shift+blockdim*l_npw*l_nspinor),&
!&       l_gvnlxc(:,shift+1:shift+blockdim*l_npw*l_nspinor),&
&       l_gvnlxc(:,:),&
&       already_transposed=.false.)
     end do
   end if
   !Compute enlout
!   do iband=1,nband
!     shift = npw*nspinor*(iband-1)
!       call dotprod_g(enl_out(iband),dprod_i,l_gs_hamk%istwf_k,npw*nspinor,1,cg(:, shift+1:shift+npw*nspinor),&
!  &     l_gvnlxc(:, shift+1:shift+npw*nspinor),mpi_enreg%me_g0,mpi_enreg%comm_bandspinorfft)
!   end do
 end if

 ABI_FREE(l_gvnlxc)

 ! Free lobpcg
 call lobpcg_free(lobpcg)

!###########################################################################
!################    SORRY IT'S ALREADY FINISHED : )  ######################
!###########################################################################


 call timab(tim_lobpcgwf2,2,tsec)
 cputime = abi_cpu_time() - cputime
 walltime = abi_wtime() - walltime
 nthreads = xomp_get_num_threads(open_parallel = .true.)
 if ( cputime/walltime/dble(nthreads) < 0.75 .and. (int(cputime/walltime)+1) /= nthreads) then
   if ( prtvol >= 3 ) then
     write(std_out,'(a)',advance='no') sjoin(" Lobpcg took", sec2str(cputime), "of cpu time")
     write(std_out,*) sjoin("for a wall time of", sec2str(walltime))
     write(std_out,'(a,f6.2)') " -> Ratio of ", cputime/walltime
   end if
   MSG_COMMENT(sjoin("You should set the number of threads to something close to",itoa(int(cputime/walltime)+1)))
 end if


 DBG_EXIT("COLL")

end subroutine lobpcgwf2


 subroutine getghc_gsc(X,AX,BX)
   use m_xg, only : xg_t, xgBlock_get, xgBlock_set, xgBlock_getSize, xgBlock_t
#ifdef HAVE_OPENMP
   use omp_lib
#endif
  type(xgBlock_t), intent(inout) :: X
  type(xgBlock_t), intent(inout) :: AX
  type(xgBlock_t), intent(inout) :: BX
  integer         :: blockdim
  integer         :: spacedim
  type(pawcprj_type) :: cprj_dum(l_gs_hamk%natom,0)
  double precision :: dum
  double precision, parameter :: inv_sqrt2 = 1/sqrt2
  double precision, pointer :: cg(:,:)
  double precision, pointer :: ghc(:,:)
  double precision, pointer :: gsc(:,:)

  call xgBlock_getSize(X,spacedim,blockdim)
  spacedim = spacedim/l_icplx


  !call xgBlock_get(X,cg(:,1:blockdim*spacedim),0,spacedim)
  call xgBlock_reverseMap(X,cg,l_icplx,spacedim*blockdim)
  call xgBlock_reverseMap(AX,ghc,l_icplx,spacedim*blockdim)
  call xgBlock_reverseMap(BX,gsc,l_icplx,spacedim*blockdim)

  ! scale back cg
 if(l_istwf == 2) then
   !cg(:,1:spacedim*blockdim) = cg(:,1:spacedim*blockdim) * inv_sqrt2
   call xgBlock_scale(X,inv_sqrt2,1)
   if(l_mpi_enreg%me_g0 == 1) cg(:, 1:spacedim*blockdim:l_npw) = cg(:, 1:spacedim*blockdim:l_npw) * sqrt2
 end if

 if ( size(l_gvnlxc) < 2*blockdim*spacedim ) then
   ABI_FREE(l_gvnlxc)
   ABI_MALLOC(l_gvnlxc,(2,blockdim*spacedim))
 end if

  if (l_mpi_enreg%paral_kgb==0) then

    call multithreaded_getghc(l_cpopt,cg(:,1:blockdim*spacedim),cprj_dum,ghc,gsc(:,1:blockdim*spacedim),&
      l_gs_hamk,l_gvnlxc,dum, l_mpi_enreg,blockdim,l_prtvol,l_sij_opt,l_tim_getghc,0)

  else
    call prep_getghc(cg(:,1:blockdim*spacedim),l_gs_hamk,l_gvnlxc,ghc,gsc(:,1:blockdim*spacedim),dum,blockdim,l_mpi_enreg,&
&                     l_prtvol,l_sij_opt,l_cpopt,cprj_dum,already_transposed=.false.)
  end if

  ! scale cg, ghc, gsc
  if ( l_istwf == 2 ) then
    !cg(:,1:spacedim*blockdim) = cg(:,1:spacedim*blockdim) * sqrt2
    !ghc(:,1:spacedim*blockdim) = ghc(:,1:spacedim*blockdim) * sqrt2
    call xgBlock_scale(X,sqrt2,1)
    call xgBlock_scale(AX,sqrt2,1)
    if(l_mpi_enreg%me_g0 == 1) then
      cg(:, 1:spacedim*blockdim:l_npw) = cg(:, 1:spacedim*blockdim:l_npw) * inv_sqrt2
      ghc(:, 1:spacedim*blockdim:l_npw) = ghc(:, 1:spacedim*blockdim:l_npw) * inv_sqrt2
    endif
    if(l_paw) then
      !gsc(:,1:spacedim*blockdim) = gsc(:,1:spacedim*blockdim) * sqrt2
      call xgBlock_scale(BX,sqrt2,1)
      if(l_mpi_enreg%me_g0 == 1) gsc(:, 1:spacedim*blockdim:l_npw) = gsc(:, 1:spacedim*blockdim:l_npw) * inv_sqrt2
    end if
  end if

  if ( .not. l_paw ) call xgBlock_copy(X,BX)

  !call xgBlock_set(AX,ghc,0,spacedim)
  !call xgBlock_set(BX,gsc(:,1:blockdim*spacedim),0,spacedim)
 end subroutine getghc_gsc

 subroutine precond(W)
   use m_xg, only : xg_t, xgBlock_colwiseMul
   type(xgBlock_t), intent(inout) :: W
   integer :: ispinor
   !integer :: cplx

   ! precondition resid_vec
   do ispinor = 1,l_nspinor
     !do cplx = 1, l_icplx
     call xgBlock_colwiseMul(W,l_pcon,l_icplx*l_npw*(ispinor-1))
      !end do
   end do

 end subroutine precond

 end module m_lobpcgwf
!!***
