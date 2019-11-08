
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

#if defined HAVE_FC_SHIFTLR
#  define SHIFTL_ shiftl
#  define SHIFTR_ shiftr
#else
#  define SHIFTL_ lshift
#  define SHIFTR_ rshift
#endif

! xoroshiro128plus method random number generator
! adapted by hexu for usage in Abinit (downgrade to Fortran 90 and added
! some functions )
! TODO: move this to 28_numeric_noabirule?

!** License for the xoroshiro128plus random number generator **
!
!Written in 2016 by David Blackman and Sebastiano Vigna (vigna@acm.org)
!Translated to Fortran 2008 by Jannis Teunissen
!
!To the extent possible under law, the author has dedicated all copyright
!and related and neighboring rights to this software to the public domain
!worldwide. This software is distributed without any warranty.
!
!See <http://creativecommons.org/publicdomain/zero/1.0/>.
!
! << This is the original documentation of Blackman and Vigna
!This is the successor to xorshift128+. It is the fastest full-period
!generator passing BigCrush without systematic failures, but due to the
!relatively short period it is acceptable only for applications with a
!mild amount of parallelism; otherwise, use a xorshift1024* generator.
!
!Beside passing BigCrush, this generator passes the PractRand test suite
!up to (and included) 16TB, with the exception of binary rank tests,
!which fail due to the lowest bit being an LFSR; all other bits pass all
!tests. We suggest to use a sign test to extract a random Boolean value.
!
!Note that the generator uses a simulated rotate operation, which most C
!compilers will turn into a single instruction. In Java, you can use
!Long.rotateLeft(). In languages that do not make low-level rotation
!instructions accessible xorshift128+ could be faster.
!
!The state must be seeded so that it is not everywhere zero. If you have
!a 64-bit seed, we suggest to seed a splitmix64 generator and use its
!output to fill s.
! =================
!

!> Module for pseudo random number generation. The internal pseudo random
!> generator is the xoroshiro128plus method.

module m_random_xoroshiro128plus

  use defs_basis
  use m_profiling_abi

  implicit none
  private

  ! A 64 bit floating point type
  !integer, parameter :: dp = kind(0.0d0)

  ! A 32 bit integer type
  integer, parameter :: i4 = selected_int_kind(9)

  ! A 64 bit integer type
  integer, parameter :: i8 = selected_int_kind(18)

  !> Random number generator type, which contains the state
  type rng_t
     !> The rng state (always use your own seed)
     integer(i8), private       :: s(2) = [123456789_i8, 987654321_i8]
     integer(i8), private       :: separator(32) ! Separate cache lines (parallel use)

     real(dp) :: residual =0.0d0 ! for saving residual in normal function.
     logical :: has_residual = .False.
  ! contains
  !   procedure, non_overridable :: set_seed    ! Seed the generator
  !   procedure, non_overridable :: jump        ! Jump function (see below)
  !   procedure, non_overridable :: rand_int4       ! 4-byte random integer
  !   procedure, non_overridable :: rand_int8       ! 8-byte random integer
  !   procedure, non_overridable :: rand_unif_01     ! Uniform (0,1] real
  !   procedure, non_overridable :: rand_unif_01_array     ! Uniform (0,1] real
  !   procedure, non_overridable :: rand_two_normals ! Two normal(0,1) samples
  !   procedure, non_overridable :: rand_normal ! Two normal(0,1) samples
  !   procedure, non_overridable :: rand_normal_array ! Two normal(0,1) samples
  !   procedure, non_overridable :: rand_poisson     ! Sample from Poisson-dist.
  !   procedure, non_overridable :: rand_circle      ! Sample on a rand_circle
  !   procedure, non_overridable :: rand_sphere      ! Sample on a rand_sphere
  !   procedure, non_overridable :: next        ! Internal method
  end type rng_t

  !> Parallel random number generator type
  type prng_t
     type(rng_t), allocatable :: rngs(:)
   !contains
   !  procedure, non_overridable :: init_parallel
   ! MG: free method not defined!!
  end type prng_t


  public :: rng_t
  public :: prng_t
  public :: set_seed
  public :: jump
  public :: rand_int4
  public :: rand_int8
  public :: rand_unif_01
  public :: rand_two_normals
  public :: rand_normal
  public :: rand_normal_array
  public :: rand_poisson
  public :: rand_circle
  public :: rand_sphere


contains

  !> Initialize a collection of rng's for parallel use
  subroutine init_parallel(self, n_proc, rng)

    class(prng_t), intent(inout) :: self
    type(rng_t), intent(inout)   :: rng
    integer, intent(in)          :: n_proc
    integer                      :: n

    ABI_MALLOC(self%rngs, (n_proc))
    self%rngs(1) = rng

    do n = 2, n_proc
       self%rngs(n) = self%rngs(n-1)
       !call self%rngs(n)%jump()
       call jump(self%rngs(n))
    end do
  end subroutine init_parallel

  !> Set a seed for the rng
  subroutine set_seed(self, the_seed)

    class(rng_t), intent(inout) :: self
    integer(i8), intent(in)     :: the_seed(2)

    self%s = the_seed

    ! Simulate calls to next() to improve randomness of first number
    !call self%jump()
    call jump(self)
  end subroutine set_seed

  ! This is the jump function for the generator. It is equivalent
  ! to 2^64 calls to next(); it can be used to generate 2^64
  ! non-overlapping subsequences for parallel computations.
  subroutine jump(self)

    class(rng_t), intent(inout) :: self
    integer                     :: i, b
    integer(i8)                 :: t(2), dummy

    ! The signed equivalent of the unsigned constants
    integer(i8), parameter      :: jmp_c(2) = &
         (/-4707382666127344949_i8, -2852180941702784734_i8/)

    t = 0
    do i = 1, 2
       do b = 0, 63
          if (iand(jmp_c(i), SHIFTL_(1_i8, b)) /= 0) then
             t = ieor(t, self%s)
          end if
          !dummy = self%next()
          dummy = next(self)
       end do
    end do

    self%s = t
  end subroutine jump

  !> Return 4-byte integer
  integer(i4) function rand_int4(self)

    class(rng_t), intent(inout) :: self
    !rand_int4 = int(self%next(), i4)
    rand_int4 = int(next(self), i4)
  end function rand_int4

  !> Return 8-byte integer
  integer(i8) function rand_int8(self)

    class(rng_t), intent(inout) :: self
    !rand_int8 = self%next()
    rand_int8 = next(self)
  end function rand_int8

  !> Get a uniform [0,1) random real (double precision)
  real(dp) function rand_unif_01(self)

    class(rng_t), intent(inout) :: self
    integer(i8)                 :: x
    real(dp)                    :: tmp

    !x   = self%next()
    x   = next(self)
    x   = ior(SHIFTL_(1023_i8, 52), SHIFTR_(x, 12))
    rand_unif_01 = transfer(x, tmp) - 1.0_dp
  end function rand_unif_01


  subroutine rand_unif_01_array(self, output, size_array)

    class(rng_t), intent(inout) :: self
    integer, intent(in) :: size_array
    real(dp), intent(inout) :: output(size_array)
    integer(i8)                 :: x, i
    real(dp)                    :: tmp

    !x   = self%next()
    do i=1, size_array
       x   = next(self)
       x   = ior(SHIFTL_(1023_i8, 52), SHIFTR_(x, 12))
       output(i) = transfer(x, tmp) - 1.0_dp
    enddo
  end subroutine rand_unif_01_array

  !> Return two normal random variates with mean 0 and variance 1.
  !> http://en.wikipedia.org/wiki/Marsaglia_polar_method
  function rand_two_normals(self) result(rands)

    class(rng_t), intent(inout) :: self
    real(dp)                    :: rands(2), sum_sq
    do
       !rands(1) = 2 * self%rand_unif_01() - 1
       rands(1) = 2 * rand_unif_01(self) - 1
       !rands(2) = 2 * self%rand_unif_01() - 1
       rands(2) = 2 * rand_unif_01(self) - 1
       sum_sq = sum(rands**2)
       if (sum_sq < 1.0_dp .and. sum_sq > 0.0_dp) exit
    end do
    rands = rands * sqrt(-2 * log(sum_sq) / sum_sq)
  end function rand_two_normals

  !> Return one normal random variates with mean 0 and variance 1.
  !> http://en.wikipedia.org/wiki/Marsaglia_polar_method
  function rand_normal(self) result(r)

    class(rng_t), intent(inout) :: self
    real(dp) :: r, rands(2), sum_sq
    if(self%has_residual) then
        self%has_residual=.False.
        r=self%residual
    else
        do
        !rands(1) = 2 * self%rand_unif_01() - 1
        rands(1) = 2 * rand_unif_01(self) - 1
        !rands(2) = 2 * self%rand_unif_01() - 1
        rands(2) = 2 * rand_unif_01(self) - 1
        sum_sq = sum(rands**2)
        if (sum_sq < 1.0_dp .and. sum_sq > 0.0_dp) exit
        end do
    rands = rands * sqrt(-2 * log(sum_sq) / sum_sq)
    r=rands(1)
    self%has_residual=.True.
    self%residual=rands(2)
    endif
  end function rand_normal

  subroutine rand_normal_array(self, a, asize)

    class(rng_t), intent(inout) :: self
    integer ,intent(in) :: asize
    real(dp), intent(inout) :: a(asize)
    integer :: i
    do i=1, asize
        a(i)=rand_normal(self)
    enddo
  end subroutine rand_normal_array


  !> Return Poisson random variate with rate lambda. Works well for lambda < 30
  !> or so. For lambda >> 1 it can produce wrong results due to roundoff error.
  function rand_poisson(self, lambda) result(rr)

    class(rng_t), intent(inout) :: self
    real(dp), intent(in)        :: lambda
    integer(i4)                 :: rr
    real(dp)                    :: expl, p

    expl = exp(-lambda)
    rr   = 0
    !p    = self%rand_unif_01()
    p    = rand_unif_01(self)

    do while (p > expl)
       rr = rr + 1
       !p = p * self%rand_unif_01()
       p = p * rand_unif_01(self)
    end do
  end function rand_poisson

  !> Sample point on a rand_circle with given radius
  function rand_circle(self, radius) result(xy)

    class(rng_t), intent(inout) :: self
    real(dp), intent(in)        :: radius
    real(dp)                    :: rands(2), xy(2)
    real(dp)                    :: sum_sq

    ! Method for uniform sampling on rand_circle
    do
       !rands(1) = 2 * self%rand_unif_01() - 1
       rands(1) = 2 * rand_unif_01(self) - 1
       !rands(2) = 2 * self%rand_unif_01() - 1
       rands(2) = 2 * rand_unif_01(self) - 1
       sum_sq   = sum(rands**2)
       if (sum_sq <= 1) exit
    end do

    xy(1) = (rands(1)**2 - rands(2)**2) / sum_sq
    xy(2) = 2 * rands(1) * rands(2) / sum_sq
    xy    = xy * radius
  end function rand_circle

  !> Sample point on a rand_sphere with given radius
  function rand_sphere(self, radius) result(xyz)

    class(rng_t), intent(inout) :: self
    real(dp), intent(in)        :: radius
    real(dp)                    :: rands(2), xyz(3)
    real(dp)                    :: sum_sq, tmp_sqrt

    ! Marsaglia method for uniform sampling on rand_sphere
    do
       !rands(1) = 2 * self%rand_unif_01() - 1
       !rands(2) = 2 * self%rand_unif_01() - 1
       rands(1) = 2 * rand_unif_01(self) - 1
       rands(2) = 2 * rand_unif_01(self) - 1
       sum_sq   = sum(rands**2)
       if (sum_sq <= 1) exit
    end do

    tmp_sqrt = sqrt(1 - sum_sq)
    xyz(1:2) = 2 * rands(1:2) * tmp_sqrt
    xyz(3)   = 1 - 2 * sum_sq
    xyz      = xyz * radius
  end function rand_sphere

  !> Interal routine: get the next value (returned as 64 bit signed integer)
  function next(self) result(res)

    class(rng_t), intent(inout) :: self
    integer(i8)                 :: res
    integer(i8)                 :: t(2)

    t         = self%s
    res       = t(1) + t(2)
    t(2)      = ieor(t(1), t(2))
    self%s(1) = ieor(ieor(rotl(t(1), 55), t(2)), SHIFTL_(t(2), 14))
    self%s(2) = rotl(t(2), 36)
  end function next

  !> Helper function for next()
  pure function rotl(x, k) result(res)

    integer(i8), intent(in) :: x
    integer, intent(in)     :: k
    integer(i8)             :: res

    res = ior(SHIFTL_(x, k), SHIFTR_(x, 64 - k))
  end function rotl

end module m_random_xoroshiro128plus
