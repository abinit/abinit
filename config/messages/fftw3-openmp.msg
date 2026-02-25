OpenMP support is enabled, but the FFT flavor is NOT
set to 'fftw3-threads'. This means the FFTW3 library
will be only linked in its single-threaded version,
and FFT operations will NOT benefit from multi-threads
parallelization.

Certain implementations \\(e.g., cprj_in_memory\\) may
suffer from this configuration.

If you intended to use multithreaded FFT, please
re-run configure with:

--with-fft-flavor=fftw3-threads
