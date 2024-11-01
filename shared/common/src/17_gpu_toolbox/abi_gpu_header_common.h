/*
 * Copyright (C) 2008-2024 ABINIT Group
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#ifndef ABI_GPU_HEADER_COMMON_H
#define ABI_GPU_HEADER_COMMON_H

//Interfaces
#ifdef __cplusplus
extern "C" {
#endif

    void abi_cabort();
    void check_gpu_mem_(const char* str);

#ifdef __cplusplus
}
#endif

#endif /* ABI_GPU_HEADER_COMMON_H */
