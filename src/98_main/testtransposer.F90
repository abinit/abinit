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
  character(len=40) :: names(8)

  double precision :: nflops, ftimes(2)
  integer :: ncount
  double precision :: times(2)
  
  integer :: ncols
  integer :: nrows
  
  integer :: counter
  integer :: debug_rank = 0

  type(xgBlock_t) :: xcgLinalg
  type(xgBlock_t) :: xcgColsRows
  type(xgTransposer_t) :: xgTransposer


  names(1662-1661) = 'xgTransposer_transpose@ColsRows'
  names(1663-1661) = 'xgTransposer_transpose@Linalg  '
  names(1664-1661) = 'xgTransposer_*@all2all         '
  names(1665-1661) = 'xgTransposer_*@gatherv         '
  names(1666-1661) = 'xgTransposer_@reorganize       '
  names(1667-1661) = 'xgTransposer_init              '
  names(1668-1661) = 'xgTransposer_free              '
  names(1669-1661) = 'xgTransposer_transpose         '

  call xmpi_init()

  npw = 6+xmpi_comm_rank(xmpi_world) !4000 big !6 small
  
  !print *, "npw", npw
  !stop
  npw = 3888/2
  nband = 192 !2000 big !4 small
  !npw = 
  ncycle = 20
  if ( xmpi_comm_size(xmpi_world) > 1 ) then
    nCpuRows = 2
    nCpuCols = xmpi_comm_size(xmpi_world)/nCpuRows
    
    !override
    nCpuRows = 1
    nCpuCols = 2
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

  ABI_MALLOC(cg, (2,npw*nband))
  ABI_MALLOC(cg0, (2,npw*nband))

  call random_number(cg)
  cg0(:,:) = cg(:,:) !this is actually copy, not just pointer assignment

  call xgBlock_map(xcgLinalg,cg,SPACE_C,npw,nband,xmpi_world)
  !print *, "NPW", npw
  !if (xmpi_comm_rank(xmpi_world) == debug_rank) then
  !  call xgBlock_print(xcgLinalg,6)
  !end if
  !print *, cg
  !stop
  
  call xgBlock_getSize(xcgLinalg,nrows,ncols)
  
  !print *, "NROWS 1", nrows
  !print *, "NCLOS 1", ncols
  
  write(std_out,*) " Complex all2all"
!  if (xmpi_comm_rank(xmpi_world) == debug_rank) then
!    print *, "xcgColsRows PRVI PUT:"
!    call xgBlock_print(xcgColsRows,6)
!  end if


  print *, "nCpuRows", nCpuRows
  print *, "nCpuCols", nCpuCols
 
  call xgTransposer_init(xgTransposer,xcgLinalg,xcgColsRows,nCpuRows,nCpuCols,STATE_LINALG,1,debug_rank)
  call tester(xcgColsRows)
  call xgTransposer_free(xgTransposer)
  call printTimes()

  write(std_out,*) " Complex gatherv"
  call xgTransposer_init(xgTransposer,xcgLinalg,xcgColsRows,nCpuRows,nCpuCols,STATE_LINALG,2,debug_rank)
  call tester(xcgColsRows)
  call xgTransposer_free(xgTransposer)
  call printTimes()

  call xgBlock_map(xcgLinalg,cg,SPACE_CR,2*npw,nband,xmpi_world)

  write(std_out,*) " Real all2all"
  call xgTransposer_init(xgTransposer,xcgLinalg,xcgColsRows,nCpuRows,nCpuCols,STATE_LINALG,1,debug_rank)
  call tester(xcgColsRows)
  call xgTransposer_free(xgTransposer)
  call printTimes()

  write(std_out,*) " Real gatherv"
  call xgTransposer_init(xgTransposer,xcgLinalg,xcgColsRows,nCpuRows,nCpuCols,STATE_LINALG,2,debug_rank)
  call tester(xcgColsRows)
  call xgTransposer_free(xgTransposer)
  call printTimes()

  ABI_FREE(cg)
  ABI_FREE(cg0)

  call xg_finalize()


 ! Writes information on file about the memory before ending mpi module, if memory profiling is enabled
 call abinit_doctor("__testtransposer")

  call xmpi_end()

  contains
!!***

!!****f* testTransposer/tester
!!
!! NAME
!! tester
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

    subroutine tester(xcgColsRows)


      type(xgBlock_t) , intent(inout) :: xcgColsRows
      
      maxt = 0
      cputime = 0
      do i=1,ncycle
        walltime = abi_wtime()
        !if (xmpi_comm_rank(xmpi_world) == debug_rank) then
        !  call xgBlock_print(xcgLinalg,6)
        !end if
        !stop     
!        if (xmpi_comm_rank(xmpi_world) == debug_rank) then
!          print *, "RANK (before transpose): ", debug_rank
!          call xgBlock_print(xcgLinalg,6)
!        end if
        call xgTransposer_transpose(xgTransposer,STATE_COLSROWS)
        !if (xmpi_comm_rank(xmpi_world) == debug_rank) then
          !print *, "RANK (after transpose): ", debug_rank
          !call xgBlock_print(xcgLinalg,6)
        !end if
        !stop
        call debug_helper(xcgColsRows) 
        stop
        
        if ( ncpucols > 1 ) then ! for 1 both states are aliased !!
          call random_number(cg)
        end if
        !call xgBlock_scale(xcgLinalg,0.d0,1)
        !call xgBlock_print(xcgLinalg,6)
        call xgTransposer_transpose(xgTransposer,STATE_LINALG)
        !call xgBlock_print(xcgLinalg,6)
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
    end subroutine tester
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
    
    subroutine debug_helper(debugBlock)
        
      type(xgBlock_t) , intent(inout) :: debugBlock
      type(xgBlock_t) :: HELPER
      
      integer :: DEBUG_ROWS = 20
      integer :: DEBUG_COLUMNS = 2

      !call xgBlock_setBlock(debugBlock, HELPER, 1, DEBUG_ROWS, DEBUG_COLUMNS) 
      !call xgBlock_print(HELPER, 100+xmpi_comm_rank(xmpi_world)) 
   
      if (xmpi_comm_size(xmpi_world) == 1) then !only one MPI proc
        call xgBlock_setBlock(debugBlock, HELPER, 1, DEBUG_ROWS, DEBUG_COLUMNS) 
        call xgBlock_print(HELPER, 200+xmpi_comm_rank(xmpi_world)) 
      
        call xgBlock_setBlock(debugBlock, HELPER, 192/2+1, DEBUG_ROWS, DEBUG_COLUMNS) 
        call xgBlock_print(HELPER, 200+xmpi_comm_rank(xmpi_world)+1) 
      else
        call xgBlock_setBlock(debugBlock, HELPER, 1, DEBUG_ROWS, DEBUG_COLUMNS) 
        call xgBlock_print(HELPER, 100+xmpi_comm_rank(xmpi_world)) 
      end if
      
      !print *, "debugBlock%rows", rows(debugBlock)
      !print *, "debugBlock%cols", cols(debugBlock)

    end subroutine debug_helper


  end program testTransposer
!!***
