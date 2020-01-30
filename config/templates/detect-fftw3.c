#include <fftw3.h>

fftw_plan *plan;
fftw_complex *a1, *a2;
dfftw_execute_dft(plan, a1, a2);
