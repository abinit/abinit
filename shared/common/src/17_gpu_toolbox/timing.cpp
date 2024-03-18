/* timing_cuda.cpp */

/*
 * Copyright (C) 2008-2024 ABINIT Group (MMancini)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_GPU_CUDA
#include "timing_cuda.cpp"
#endif
#ifdef HAVE_GPU_HIP
#include "timing_hip.cpp"
#endif
