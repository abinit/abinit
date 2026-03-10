OpenMP support is enabled, but the FFT flavor is NOT
set to 'fftw3-threads'. This means the FFTW3 library
will be only linked in its single-threaded version,
and FFT operations will NOT benefit from multi-threads
parallelization.

If you intended to use multithreaded FFT \\(recommended\\),
please re-run configure with:

--with-fft-flavor=fftw3-threads
