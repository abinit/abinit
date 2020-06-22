include "fftw3.f03"
integer, parameter :: N = 10
complex :: a1(N), a2(N)
integer :: plan

call sfftw_plan_dft_1d(plan, N, a1, a2, FFTW_FORWARD, FFTW_ESTIMATE)
call sfftw_execute_dft(plan, a1, a2)
call sfftw_destroy_plan(plan)
