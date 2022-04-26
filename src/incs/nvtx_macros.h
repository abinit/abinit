i/*
 * Copyright (C) 2008-2021 ABINIT Group
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */
#ifndef ABINIT_52_MANAGE_CUDA_NVTX_MACRO_H
#define ABINIT_52_MANAGE_CUDA_NVTX_MACRO_H

#include "config.h"

/*
 * Note:
 * nvtx_activated is a boolean variable defined in module
 * m_nvtx_data (52_manage_gpu/m_nvtx_data.F90).
 *
 * It can only be true if GPU/Cuda is enabled and Cuda version >= 10.
 *
 * We need these macro because subroutine abi_nvtx_start_range and abi_nvtx_end_range
 * only exists when GPU is enabled.
 */

#if defined(HAVE_GPU_CUDA) && defined(HAVE_GPU_NVTX_V3)
#define ABI_NVTX_START_RANGE(id) call abi_nvtx_start_range(id)
#define ABI_NVTX_END_RANGE() call abi_nvtx_end_range()
#define NVTX_INIT(value) call nvtx_init(value)
#else
#define ABI_NVTX_START_RANGE(id)
#define ABI_NVTX_END_RANGE()
#define NVTX_INIT(value)
#endif

#endif /* ABINIT_52_MANAGE_CUDA_NVTX_MACRO_H */
