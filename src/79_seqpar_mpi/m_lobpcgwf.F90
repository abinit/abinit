!{\src2tex{textfont=tt}}
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
!! Copyright (C) 1998-2018 ABINIT group (JB)
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

 use defs_abitypes
 use defs_basis
 use m_profiling_abi
 use m_lobpcg
 use m_xmpi
 use m_errors
 use m_time
 use m_xomp
 use m_fstrings
 use m_xg
 use m_xgTransposer
 use m_lobpcg2

 use m_hamiltonian, only : gs_hamiltonian_type
 use m_pawcprj,     only : pawcprj_type
 use m_nonlop,      only : nonlop

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
 real(dp), allocatable,save :: l_gvnlc(:,:)
 real(dp), allocatable,save ::  l_pcon(:)
 type(mpi_type),pointer,save :: l_mpi_enreg
 type(gs_hamiltonian_type),pointer,save :: l_gs_hamk

 public :: lobpcgwf2
  contains

subroutine lobpcgwf2(cg,dtset,eig,enl_out,gs_hamk,kinpw,mpi_enreg,&
&                   nband,npw,nspinor,prtvol,resid)


 use m_cgtools, only : dotprod_g
 use iso_c_binding

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lobpcgwf2'
 use interfaces_66_wfs
!End of the abilint section

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

 ! Stupid things for NC
 integer,parameter :: choice=1, paw_opt=0, signs=1
 type(pawcprj_type) :: cprj_dum(gs_hamk%natom,0)
 integer :: iband, shift
 real(dp) :: gsc_dummy(0,0)

 real(dp), allocatable :: first(:,:)
 integer :: ncpuCols, ncpuRows
 real(dp) :: errmax, maxt
 type(xgTransposer_t) :: xgTransposer

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
! cptr = c_loc(eig)
! call c_f_pointer(cptr,eig_ptr,(/ nband,1 /))
! call xgBlock_map(xgeigen,eig_ptr,SPACE_R,nband,1)

 !call xg_init(xgresidu,SPACE_R,nband,1,l_mpi_enreg%comm_bandspinorfft)
 ! Trick the with C to change rank of arrays (:) to (:,:)
! cptr = c_loc(resid)
! call c_f_pointer(cptr,resid_ptr,(/ nband,1 /))
! call xgBlock_map(xgresidu,resid_ptr,SPACE_R,nband,1)

 !ABI_MALLOC(l_gvnlc,(2,l_npw*l_nspinor*blockdim))

 !call lobpcg_init(lobpcg,nband, l_icplx*l_npw*l_nspinor, blockdim,dtset%tolwfr,nline,space, l_mpi_enreg%comm_bandspinorfft)

!------------------------ TEST --------------------------!

! Copy initial xgx0
ABI_MALLOC(first,(2,l_npw*l_nspinor*nband))
first = cg
!call xgBlock_get(xgx0,first,0,l_npw*l_nspinor)

! Get full number of npw
call xmpi_sum(l_npw,mpi_enreg%comm_bandspinorfft,iband)
! Allocate memory
!ABI_MALLOC(l_gvnlc,(2,l_npw*l_nspinor*nband/dtset%npband))
! Map a xgBlock onto this memory
!call xgBlock_map(xgeigen,l_gvnlc,space,l_icplx*l_npw*l_nspinor,nband/dtset%npband,l_mpi_enreg%comm_fft)

cputime = 0
ncpuCols = dtset%npband
ncpuRows = dtset%npfft*dtset%npspinor

call xgTransposer_init(xgTransposer,xgx0,xgeigen,ncpuRows,ncpuCols,STATE_LINALG,dtset%useric)

do iband =1 , 10
  walltime = abi_wtime()
  call xgTransposer_transpose(xgTransposer,STATE_COLSROWS)
  call random_number(cg)
  !call xgBlock_scale(xgx0,0.d0,1)
  !call xgBlock_print(xgeigen,6)
  call xgTransposer_transpose(xgTransposer,STATE_LINALG)
  !call xgBlock_print(xgx0,6)
  call xmpi_barrier(l_mpi_enreg%comm_bandspinorfft)
  walltime = abi_wtime() - walltime
  cputime = cputime + walltime
  call xmpi_max(walltime,maxt,mpi_enreg%comm_bandspinorfft,nthreads)
  write(std_out,*) "walltime", maxt
end do
!ABI_FREE(l_gvnlc)
call xmpi_max(cputime,maxt,mpi_enreg%comm_bandspinorfft,nthreads)
call xgTransposer_free(xgTransposer)
write(std_out,*) "mean", maxt/l_istwf
errmax = (sum(first-cg))/nband
call xmpi_sum(errmax,l_mpi_enreg%comm_bandspinorfft,nthreads)
write(std_out,*) "difference:",errmax
call flush(std_out)
ABI_FREE(first)
call xmpi_abort()

 ! Free lobpcg
 call lobpcg_free(lobpcg)

 DBG_EXIT("COLL")

end subroutine lobpcgwf2



 end module m_lobpcgwf

!!***
