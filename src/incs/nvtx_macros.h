/*
 * Copyright (C) 2008-2025 ABINIT Group
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */
#ifndef ABINIT_NVTX_MACRO_H
#define ABINIT_NVTX_MACRO_H

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

/*
 * Note:
 * nvtx_activated is a boolean variable defined in module
 * m_nvtx_data (44_abitools/m_nvtx_data.F90).
 *
 * It can only be true if GPU (NVIDIA CUDA > v10 or AMD ROCm) is enabled.
 *
 * We need these macro because subroutine abi_nvtx_start_range and abi_nvtx_end_range
 * only exists when GPU markers are enabled.
 */

#if defined(HAVE_GPU_MARKERS)
#define ABI_NVTX_START_RANGE(id) call abi_nvtx_start_range(id)
#define ABI_NVTX_END_RANGE() call abi_nvtx_end_range()
#define NVTX_INIT() call nvtx_init()
#define NVTX_PROFILER_START() call nvtxProfilerStart()
#define NVTX_PROFILER_STOP() call nvtxProfilerStop()
#else
#define ABI_NVTX_START_RANGE(id)
#define ABI_NVTX_END_RANGE()
#define NVTX_INIT()
#define NVTX_PROFILER_START()
#define NVTX_PROFILER_STOP()
#endif

#endif /* ABINIT_NVTX_MACRO_H */
