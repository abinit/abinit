program testTransposer
  use m_xg
  use m_xgTransposer
  use m_xmpi
  use m_time
  use defs_basis

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
  double precision, allocatable :: cg(:,:)
  double precision, allocatable :: cg0(:,:)

  type(xgBlock_t) :: xcgLinalg
  type(xgBlock_t) :: xcgColsRows
  type(xgTransposer_t) :: xgTransposer

  call xmpi_init()

  npw = 2000+xmpi_comm_rank(xmpi_world)
  nband = 200
  ncycle = 100

  std_out = 6+xmpi_comm_rank(xmpi_world)

  allocate(cg(2,npw*nband))
  allocate(cg0(2,npw*nband))

  call random_number(cg)
  cg0(:,:) = cg(:,:)

  call xgBlock_map(xcgLinalg,cg,SPACE_C,npw,nband,xmpi_world)

  write(std_out,*) "Complex iall2all"
  call xgTransposer_init(xgTransposer,xcgLinalg,xcgColsRows,2,4,STATE_LINALG,1)
  call tester()
  call xgTransposer_free(xgTransposer)

  write(std_out,*) "Complex igatherv"
  call xgTransposer_init(xgTransposer,xcgLinalg,xcgColsRows,2,4,STATE_LINALG,2)
  call tester()
  call xgTransposer_free(xgTransposer)

  call xgBlock_map(xcgLinalg,cg,SPACE_CR,2*npw,nband,xmpi_world)

  write(std_out,*) "Real iall2all"
  call xgTransposer_init(xgTransposer,xcgLinalg,xcgColsRows,2,4,STATE_LINALG,1)
  call tester()
  call xgTransposer_free(xgTransposer)

  write(std_out,*) "Real igatherv"
  call xgTransposer_init(xgTransposer,xcgLinalg,xcgColsRows,2,4,STATE_LINALG,2)
  call tester()
  call xgTransposer_free(xgTransposer)

  deallocate(cg)
  deallocate(cg0)

  call xg_finalize()

  call xmpi_end()

  contains

    subroutine tester()
      maxt = 0
      do i=1,ncycle
        walltime = abi_wtime()
        call xgTransposer_transpose(xgTransposer,STATE_COLSROWS)
        call random_number(cg)
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
      write(std_out,*) "Mean time: ", maxt/ncycle
      errmax = (sum(cg0-cg))/nband
      call xmpi_sum(errmax,xmpi_world,ierr)
      write(std_out,*) "Difference: ",errmax
      call flush(std_out)
      call xmpi_barrier(xmpi_world)
    end subroutine tester



  end program testTransposer
