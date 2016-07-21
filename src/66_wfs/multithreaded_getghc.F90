!{\src2tex{textfont=tt}}
!!****f* ABINIT/getghc
!!
!! NAME
!! multithreaded_getghc
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2016 ABINIT group (JB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cpopt=flag defining the status of cwaveprj%cp(:)=<Proj_i|Cnk> scalars (PAW only)
!!       (same meaning as in nonlop.F90 routine)
!!       if cpopt=-1, <p_lmn|in> (and derivatives) are computed here (and not saved)
!!       if cpopt= 0, <p_lmn|in> are computed here and saved
!!       if cpopt= 1, <p_lmn|in> and first derivatives are computed here and saved
!!       if cpopt= 2  <p_lmn|in> are already in memory;
!!       if cpopt= 3  <p_lmn|in> are already in memory; first derivatives are computed here and saved
!!       if cpopt= 4  <p_lmn|in> and first derivatives are already in memory;
!! cwavef(2,npw*my_nspinor*ndat)=planewave coefficients of wavefunction.
!! gs_ham <type(gs_hamiltonian_type)>=all data for the Hamiltonian to be applied
!! lambda=factor to be used when computing <G|H-lambda.S|C> - only for sij_opt=-1
!!        Typically lambda is the eigenvalue (or its guess)
!! mpi_enreg=informations about MPI parallelization
!! ndat=number of FFT to do in parallel
!! prtvol=control print volume and debugging output
!! sij_opt= -PAW ONLY-  if  0, only matrix elements <G|H|C> have to be computed
!!    (S=overlap)       if  1, matrix elements <G|S|C> have to be computed in gsc in addition to ghc
!!                      if -1, matrix elements <G|H-lambda.S|C> have to be computed in ghc (gsc not used)
!! tim_getghc=timing code of the calling subroutine(can be set to 0 if not attributed)
!! type_calc= option governing which part of Hamitonian is to be applied:
!             0: whole Hamiltonian
!!            1: local part only
!!            2: non-local+kinetic only (added to the exixting Hamiltonian)
!!            3: local + kinetic only (added to the existing Hamiltonian)
!! ===== Optional inputs =====
!!   [kg_fft_k(3,:)]=optional, (k+G) vector coordinates to be used for the FFT tranformation
!!                   instead of the one contained in gs_ham datastructure.
!!                   Typically used for real WF (in parallel) which are FFT-transformed 2 by 2.
!!   [kg_fft_kp(3,:)]=optional, (k^prime+G) vector coordinates to be used for the FFT tranformation
!!   [select_k]=optional, option governing the choice of k points to be used.
!!             gs_ham datastructure contains quantities needed to apply Hamiltonian
!!             in reciprocal space between 2 kpoints, k and k^prime (equal in most cases);
!!             if select_k=1, <k^prime|H|k>       is applied [default]
!!             if select_k=2, <k|H|k^prime>       is applied
!!             if select_k=3, <k|H|k>             is applied
!!             if select_k=4, <k^prime|H|k^prime> is applied
!!
!! OUTPUT
!!  ghc(2,npw*my_nspinor*ndat)=matrix elements <G|H|C> (if sij_opt>=0)
!!                                          or <G|H-lambda.S|C> (if sij_opt=-1)
!!  gvnlc(2,npw*my_nspinor*ndat)=matrix elements <G|Vnonlocal|C> (if sij_opt>=0)
!!                                            or <G|Vnonlocal-lambda.S|C> (if sij_opt=-1)
!!  if (sij_opt=1)
!!    gsc(2,npw*my_nspinor*ndat)=matrix elements <G|S|C> (S=overlap).
!!
!! SIDE EFFECTS
!!  cwaveprj(natom,my_nspinor*(1+cpopt)*ndat)= wave function projected on nl projectors (PAW only)
!!
!! PARENTS
!!      cgwf,chebfi,dfpt_cgwf,gwls_hamiltonian,ks_ddiago,lobpcgwf,m_io_kss
!!      mkresi,prep_getghc
!!
!! CHILDREN
!!      fourwf
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine multithreaded_getghc(cpopt,cwavef,cwaveprj,ghc,gsc,gs_ham,gvnlc,lambda,mpi_enreg,ndat,&
&                 prtvol,sij_opt,tim_getghc,type_calc,&
&                 kg_fft_k,kg_fft_kp,select_k) ! optional arguments

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

 use m_pawcprj,     only : pawcprj_type,pawcprj_alloc,pawcprj_free,pawcprj_getdim
 use m_bandfft_kpt, only : bandfft_kpt,bandfft_kpt_get_ikpt
 use m_hamiltonian, only : gs_hamiltonian_type,KPRIME_H_K,K_H_KPRIME,K_H_K,KPRIME_H_KPRIME
 use m_fock,        only : fock_type,fock_get_getghc_call

#ifdef HAVE_OPENMP
   use omp_lib
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'multithreaded_getghc'
 use interfaces_66_wfs, except_this_one => multithreaded_getghc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cpopt,ndat, prtvol
 integer,intent(in) :: sij_opt,tim_getghc,type_calc
 integer,intent(in),optional :: select_k
 real(dp),intent(in) :: lambda
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout),target :: gs_ham
!arrays
 integer,intent(in),optional,target :: kg_fft_k(:,:),kg_fft_kp(:,:)
 real(dp),intent(out),target :: gsc(:,:)
 real(dp),intent(inout) :: cwavef(:,:)
 real(dp),intent(out) :: ghc(:,:),gvnlc(:,:)
 type(pawcprj_type),intent(inout),target :: cwaveprj(:,:)

!Local variables-------------------------------
!scalars
 integer :: firstelt, lastelt
 integer :: nthreads
 integer :: ithread
 integer :: chunk
 integer :: residuchunk
 integer :: firstband
 integer :: lastband
 integer :: spacedim
#ifdef HAVE_OPENMP
 logical :: is_nested
#endif

 integer :: select_k_default

 ! *************************************************************************

   select_k_default = 1; if ( present(select_k) ) select_k_default = select_k

   spacedim = size(cwavef,dim=2)/ndat

    !$omp parallel default (none) private(ithread,nthreads,chunk,firstband,lastband,residuchunk,firstelt,lastelt, is_nested), &
    !$omp& shared(cwavef,ghc,gsc, gvnlc,spacedim,ndat,kg_fft_k,kg_fft_kp,gs_ham,cwaveprj,mpi_enreg), &
    !$omp& firstprivate(cpopt,lambda,prtvol,sij_opt,tim_getghc,type_calc,select_k_default)
#ifdef HAVE_OPENMP
    ithread = omp_get_thread_num()
    nthreads = omp_get_num_threads()
    is_nested = omp_get_nested()
    call omp_set_nested(.false.)
#ifdef HAVE_LINALG_MKL_THREADS
    call mkl_set_num_threads(1)
#endif
#else
    ithread = 0
    nthreads = 1
#endif
    chunk = ndat/nthreads ! Divide by 2 to construct chunk of even number of bands
    residuchunk = ndat - nthreads*chunk
    if ( ithread < nthreads-residuchunk ) then
      firstband = ithread*chunk+1
      lastband = (ithread+1)*chunk
    else
      firstband = (nthreads-residuchunk)*chunk + ( ithread -(nthreads-residuchunk) )*(chunk+1) +1
      lastband = firstband+chunk
    endif

    if ( lastband /= 0 ) then
      firstelt = (firstband-1)*spacedim+1
      lastelt = lastband*spacedim
      ! Don't know how to manage optional arguments .... :(
      if ( present(kg_fft_k) ) then
        if (present(kg_fft_kp)) then
          call getghc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj,ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*gs_ham%usepaw),&
            gs_ham,gvnlc(:,firstelt:lastelt),lambda, mpi_enreg,lastband-firstband+1,prtvol,sij_opt,tim_getghc,type_calc,&
            select_k=select_k_default,kg_fft_k=kg_fft_k,kg_fft_kp=kg_fft_kp)
        else
          call getghc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj,ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*gs_ham%usepaw),&
            gs_ham,gvnlc(:,firstelt:lastelt),lambda, mpi_enreg,lastband-firstband+1,prtvol,sij_opt,tim_getghc,type_calc,&
            select_k=select_k_default,kg_fft_k=kg_fft_k)
        endif
      else
        if (present(kg_fft_kp)) then
          call getghc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj,ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*gs_ham%usepaw),&
            gs_ham,gvnlc(:,firstelt:lastelt),lambda, mpi_enreg,lastband-firstband+1,prtvol,sij_opt,tim_getghc,type_calc,&
            select_k=select_k_default,kg_fft_kp=kg_fft_kp)
        else
          call getghc(cpopt,cwavef(:,firstelt:lastelt),cwaveprj,ghc(:,firstelt:lastelt),gsc(:,firstelt:lastelt*gs_ham%usepaw),&
            gs_ham,gvnlc(:,firstelt:lastelt),lambda, mpi_enreg,lastband-firstband+1,prtvol,sij_opt,tim_getghc,type_calc,&
            select_k=select_k_default)
        endif
      end if
    end if
#ifdef HAVE_OPENMP
    call omp_set_nested(is_nested)
#ifdef HAVE_LINALG_MKL_THREADS
    call mkl_set_num_threads(nthreads)
#endif
#endif
    !$omp end parallel

end subroutine multithreaded_getghc
!!***

