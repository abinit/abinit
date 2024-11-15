/* dev_spec.cpp*/

/*
 * Copyright (C) 2008-2024 ABINIT Group (MMancini,FDahm)
 * this file is distributed under the terms of the
 * gnu general public license, see ~abinit/COPYING
 * or http://www.gnu.org/copyleft/gpl.txt.
 *
 */

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_GPU_CUDA
#include "dev_spec_cuda.cpp"
#endif
#ifdef HAVE_GPU_HIP
#include "dev_spec_hip.cpp"
#endif


