!{\src2tex{textfont=tt}}
!!****p* ABINIT/testTransposer
!! NAME
!! testTransposer
!!
!! FUNCTION
!! test the xgTransposer module with 8 MPI. No more no less.
!! It includes testing of complex and real numbers, and all2all and gatherv
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (JB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
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
program testTransposer
  use m_xg
  use m_xgTransposer
  use m_xmpi
  use m_time
  use defs_basis
  use m_profiling_abi
  use m_errors

  implicit none

  integer :: npw
  integer :: nband
  integer :: ncycle
  integer :: i
  integer :: ierr
  double precision :: errmax
  double precision :: walltime
  double precision :: cputime
  double precision :: maxt
  integer :: nCpuCols, nCpuRows
  double precision, allocatable :: cg(:,:)
  double precision, allocatable :: cg0(:,:)
  double precision, allocatable :: gh(:,:)
  double precision, allocatable :: ghc(:,:)
  character(len=40) :: names(8)

  double precision :: nflops, ftimes(2)
  integer :: ncount
  double precision :: times(2)

  type(xgBlock_t) :: xcgLinalg
  type(xgBlock_t) :: xcgColsRows
  type(xgBlock_t) :: xghLinalg
  type(xgBlock_t) :: xghColsRows
  type(xgBlock_t) :: xghcLinalg
  type(xgBlock_t) :: xghcColsRows
  type(xgTransposer_t) :: xgTransposer


  names(1662-1661) = 'xgTransposer_transpose@ColsRows'
  names(1663-1661) = 'xgTransposer_transpose@Linalg  '
  names(1664-1661) = 'xgTransposer_*@all2all         '
  names(1665-1661) = 'xgTransposer_*@gatherv         '
  names(1666-1661) = 'xgTransposer_@reorganize       '
  names(1667-1661) = 'xgTransposer_*constructor      '
  names(1668-1661) = 'xgTransposer_free              '
  names(1669-1661) = 'xgTransposer_transpose         '

  call xmpi_init()

  npw = 4000+xmpi_comm_rank(xmpi_world)
  nband = 2000
  ncycle = 20
  if ( xmpi_comm_size(xmpi_world) > 1 ) then
    nCpuRows = 2
    nCpuCols = xmpi_comm_size(xmpi_world)/nCpuRows
  else
    nCpuRows = 1
    nCpuCols = 1
  end if

  std_out = 6+xmpi_comm_rank(xmpi_world)

 ! Initialize memory profiling if it is activated
 ! if a full memocc.prc report is desired, set the argument of abimem_init to "2" instead of "0"
 ! note that memocc.prc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

  call test1()
  
  call test2()

  call xg_finalize()


 ! Writes information on file about the memory before ending mpi module, if memory profiling is enabled
 call abinit_doctor("__testtransposer")

  call xmpi_end()

  contains
!!***

!!****f* testTransposer/test1
!!
!! NAME
!! test1
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (JB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  subroutine test1()
    ABI_MALLOC(cg, (2,npw*nband))
    ABI_MALLOC(cg0, (2,npw*nband))
  
    call random_number(cg)
    cg0(:,:) = cg(:,:)
  
    call xgBlock_map(xcgLinalg,cg,SPACE_C,npw,nband,xmpi_world)
  
    write(std_out,*) " Complex all2all"
    call xgTransposer_constructor(xgTransposer,xcgLinalg,xcgColsRows,nCpuRows,nCpuCols,STATE_LINALG,TRANS_ALL2ALL)
    call backAndForth()
    call xgTransposer_free(xgTransposer)
    call printTimes()
  
    write(std_out,*) " Complex gatherv"
    call xgTransposer_constructor(xgTransposer,xcgLinalg,xcgColsRows,nCpuRows,nCpuCols,STATE_LINALG,TRANS_GATHER)
    call backAndForth()
    call xgTransposer_free(xgTransposer)
    call printTimes()
  
    call xgBlock_map(xcgLinalg,cg,SPACE_CR,2*npw,nband,xmpi_world)
  
    write(std_out,*) " Real all2all"
    call xgTransposer_constructor(xgTransposer,xcgLinalg,xcgColsRows,nCpuRows,nCpuCols,STATE_LINALG,TRANS_ALL2ALL)
    call backAndForth()
    call xgTransposer_free(xgTransposer)
    call printTimes()
  
    write(std_out,*) " Real gatherv"
    call xgTransposer_constructor(xgTransposer,xcgLinalg,xcgColsRows,nCpuRows,nCpuCols,STATE_LINALG,TRANS_GATHER)
    call backAndForth()
    call xgTransposer_free(xgTransposer)
    call printTimes()
  
    ABI_FREE(cg)
    ABI_FREE(cg0)
  end subroutine test1
!!***


!!****f* testTransposer/test2
!!
!! NAME
!! test2
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (JB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  subroutine test2()
    type(xgTransposer_t) :: xgTransposerGh
    type(xgTransposer_t) :: xgTransposerGhc
    type(xg_t) :: dotLinalg
    type(xg_t) :: dotColsRows
    double precision :: maxdiff

    write(std_out,*) "Allocation"
    ABI_MALLOC(cg,(2,npw*nband))
    ABI_MALLOC(gh,(2,npw*nband))
    ABI_MALLOC(ghc,(2,npw*nband))

    write(std_out,*) "Mapping"
    call xgBlock_map(xcgLinalg,cg,SPACE_C,npw,nband,xmpi_world)
    call xgBlock_map(xghLinalg,gh,SPACE_C,npw,nband,xmpi_world)
    call xgBlock_map(xghcLinalg,ghc,SPACE_C,npw,nband,xmpi_world)

    write(std_out,*) "Transposer constructor"
    call xgTransposer_constructor(xgTransposer,xcgLinalg,xcgColsRows,nCpuRows,nCpuCols,STATE_LINALG,TRANS_ALL2ALL)
    call xgTransposer_copyConstructor(xgTransposerGh,xgTransposer,xghLinalg,xghColsRows,STATE_LINALG)
    call xgTransposer_copyConstructor(xgTransposerGhc,xgTransposer,xghcLinalg,xghcColsRows,STATE_LINALG)
    
    write(std_out,*) "Init data"
    call random_number(cg)
    call random_number(gh)
    !call initVectors()

    write(std_out,*) "Linalg division"
    call xgBlock_colwiseDivision(xcgLinalg,xghLinalg,xghcLinalg)

    write(std_out,*) "Linalg norm2"
    call xg_init(dotLinalg,SPACE_R,nband,1,xmpi_world)
    call xgBlock_colwiseNorm2(xghcLinalg,dotLinalg%self)
    !call xgBlock_print(dotLinalg%self,std_out)

    write(std_out,*) "Transposer transpose"
    call xgTransposer_transpose(xgTransposer,STATE_COLSROWS)
    call xgTransposer_transpose(xgTransposerGh,STATE_COLSROWS)
    call xgTransposer_transpose(xgTransposerGhc,STATE_COLSROWS)

    write(std_out,*) "ColsRows divisions"
    call xgBlock_colwiseDivision(xcgColsRows,xghColsRows,xghcColsRows)

    write(std_out,*) "Transposer transpose back"
    call xgTransposer_transpose(xgTransposerGhc,STATE_LINALG)

    write(std_out,*) "ColsRows norm2"
    call xg_init(dotColsRows,SPACE_R,nband,1,xmpi_world)
    call xgBlock_colwiseNorm2(xghcLinalg,dotColsRows%self)
    !call xgBlock_print(dotColsRows%self,std_out)

    write(std_out,*) "Compare"
    call xgBlock_saxpy(dotLinalg%self, -1.0d0, dotColsRows%self)
    call xgBlock_reshape(dotLinalg%self, (/1,nband/))
    call xgBlock_colwiseNorm2(dotLinalg%self,dotColsRows%self,max_val=maxdiff)
    write(std_out,"(a,f20.4)") " Difference: ",sqrt(maxdiff)

    write(std_out,*) "Free everything"
    call xg_free(dotLinalg)
    call xg_free(dotColsRows)

    call xgTransposer_free(xgTransposer)
    call xgTransposer_free(xgTransposerGh)
    call xgTransposer_free(xgTransposerGhc)

    ABI_FREE(cg)
    ABI_FREE(gh)
    ABI_FREE(ghc)

  end subroutine test2

  subroutine initVectors()
    integer, allocatable :: seed(:)
    integer :: n, iseed
    integer :: icol, irow

    call random_seed(size=n)
    ABI_MALLOC(seed,(n))
    do icol = 1, nband
      do iseed = 1, n
        seed(iseed) = (xmpi_comm_rank(xmpi_world)*nband+icol)*n+iseed
      end do
      call random_seed(put=seed)
      do irow = 1, npw
        call random_number(cg(:,(icol-1)*npw+1:icol*npw))
        call random_number(gh(:,(icol-1)*npw+1:icol*npw))
      end do
    end do
    ABI_FREE(seed)
  end subroutine initVectors

!!****f* testTransposer/backAndForth
!!
!! NAME
!! backAndForth
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (JB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

    subroutine backAndForth()

      maxt = 0
      cputime = 0
      do i=1,ncycle
        walltime = abi_wtime()
        call xgTransposer_transpose(xgTransposer,STATE_COLSROWS)
        if ( ncpucols > 1 ) then ! for 1 both states are aliased !!
          call random_number(cg)
        end if
        !call xgBlock_scale(xcgLinalg,0.d0,1)
        !call xgBlock_print(xgeigen,6)
        call xgTransposer_transpose(xgTransposer,STATE_LINALG)
        !call xgBlock_print(xgx0,6)
        call xmpi_barrier(xmpi_world)
        walltime = abi_wtime() - walltime
        cputime = cputime + walltime
        call xmpi_max(walltime,maxt,xmpi_world,ierr)
      end do
      call xmpi_max(cputime,maxt,xmpi_world,ierr)
      write(std_out,"(a,f20.4)") ".Mean time: ", maxt/ncycle
      errmax = (sum(cg0-cg))/nband
      call xmpi_sum(errmax,xmpi_world,ierr)
      write(std_out,"(a,f20.4)") " Difference: ",errmax
      call xmpi_barrier(xmpi_world)
    end subroutine backAndForth
!!***

!!****f* testTransposer/printTimes
!!
!! NAME
!! printTimes
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (JB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

    subroutine printTimes()

      double precision :: total(2)
      integer :: ntot
      write(std_out,'(1x,a30,a8,a17,a17)') "counter", "calls", "cpu_time", "wall_time"
      ntot = 0
      total(:) = 0.d0
      do i=1662,1669
        call time_accu(i,ncount,times,nflops,ftimes)
        total(1) = total(1) + times(1)
        total(2) = total(2) + times(2)
        ntot = ntot + ncount
        write(std_out,'(a,a30,i8,2F17.3)') "-",trim(names(i-1661)), ncount, times(1), times(2)
      end do
      write(std_out,'(a,a30,i8,2F17.3)') "-","total", ntot, total(1), total(2)
      call timab(1,0,times)

    end subroutine printTimes


  end program testTransposer
!!***

