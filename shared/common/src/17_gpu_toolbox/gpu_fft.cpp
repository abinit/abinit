/* gpu_fft.cpp */

/*
 * Copyright (C) 2008-2024 ABINIT Group (MSarraute)
 * this file is distributed under the terms of the
 * gnu general public license, see ~abinit/COPYING
 * or http://www.gnu.org/copyleft/gpl.txt .
 * for the initials of contributors, see ~abinit/doc/developers/contributors.txt.
 *
 * The main goal of this file is to contain GPU linear algebra encapsulation routines,
 * that will be callable from fortran routines, by pointing to the relevant
 * GPU runtime libraries (CUDA for NVIDIA targets, HIP for AMD targets).
 *
 */

#include <gpu_fft.h>

#ifdef HAVE_GPU_CUDA
#include "gpu_fft_cuda.cpp"
#endif
#ifdef HAVE_GPU_HIP
#include "gpu_fft_hip.cpp"
#endif
