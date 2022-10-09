// abi_common_kokkos.h

/*
 * Copyright (C) 2022 ABINIT Group
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#ifndef _ABINIT_COMMON_KOKKOS_H
#define _ABINIT_COMMON_KOKKOS_H

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include <Kokkos_Core.hpp>
#include <typeinfo>

#include <cstdio> // for printf
#include <cstdlib> // for EXIT_SUCCESS
#include <iostream>
#include <cstdint>
#include <cassert>

// Use (void) to silence unused warnings.
#define assertm(exp, msg) assert(((void)msg, exp))

#ifndef STR
#define STR(x) #x
#endif

#define MY_ASSERT(x)

#ifndef ABI_CHECK
#define ABI_CHECK(expr, msg) if (!(expr)) { printf("Abinit failed: (%s), in file %s, line %d : %s.\n", STR(expr), __FILE__, __LINE__, msg); abort(); }
#endif

using memory_space = Kokkos::DefaultExecutionSpace::memory_space;
using memory_trait = Kokkos::MemoryTraits<Kokkos::Unmanaged>;

//! device views types for array of doubles, unmanaged memory (allocated on fortran side)
using AbiView_r64_1d = Kokkos::View<double*, memory_space, memory_trait >;
using AbiView_r64_2d = Kokkos::View<double**, Kokkos::LayoutLeft,   memory_space, memory_trait >;
using AbiView_r64_3d = Kokkos::View<double***, Kokkos::LayoutLeft,  memory_space, memory_trait >;
using AbiView_r64_4d = Kokkos::View<double****, Kokkos::LayoutLeft, memory_space, memory_trait >;

//! device views types for array of doubles, unmanaged memory (allocated on fortran side), const
using AbiView_r64_1d_const = Kokkos::View<const double*, memory_space, memory_trait >;
using AbiView_r64_2d_const = Kokkos::View<const double**, Kokkos::LayoutLeft,   memory_space, memory_trait >;
using AbiView_r64_3d_const = Kokkos::View<const double***, Kokkos::LayoutLeft,  memory_space, memory_trait >;
using AbiView_r64_4d_const = Kokkos::View<const double****, Kokkos::LayoutLeft, memory_space, memory_trait >;

//! device views types for array of doubles, managed memory
using AbiView_r64_3d_managed = Kokkos::View<double***, Kokkos::LayoutLeft,  memory_space >;


//! device views types for array of integers, unmanaged memory (allocated on fortran side)
using AbiView_i32_1d = Kokkos::View<int32_t*, memory_space, memory_trait >;
using AbiView_i32_2d = Kokkos::View<int32_t**, Kokkos::LayoutLeft,   memory_space, memory_trait >;
using AbiView_i32_3d = Kokkos::View<int32_t***, Kokkos::LayoutLeft,   memory_space, memory_trait >;

//! device views types for array of complex doubles, unmanaged memory (allocated on fortran side)
using cplx_t = Kokkos::complex<double>;
using AbiView_c64_1d = Kokkos::View<cplx_t*, memory_space, memory_trait >;
using AbiView_c64_2d = Kokkos::View<cplx_t**, Kokkos::LayoutLeft,   memory_space, memory_trait >;
using AbiView_c64_3d = Kokkos::View<cplx_t***, Kokkos::LayoutLeft,  memory_space, memory_trait >;
using AbiView_c64_4d = Kokkos::View<cplx_t****, Kokkos::LayoutLeft, memory_space, memory_trait >;

//! device views types for array of complex doubles, unmanaged memory (allocated on fortran side), const
using AbiView_c64_1d_const = Kokkos::View<const cplx_t*, memory_space, memory_trait >;
using AbiView_c64_2d_const = Kokkos::View<const cplx_t**, Kokkos::LayoutLeft,   memory_space, memory_trait >;
using AbiView_c64_3d_const = Kokkos::View<const cplx_t***, Kokkos::LayoutLeft,  memory_space, memory_trait >;
using AbiView_c64_4d_const = Kokkos::View<const cplx_t****, Kokkos::LayoutLeft, memory_space, memory_trait >;


#endif // _ABINIT_COMMON_KOKKOS_H
