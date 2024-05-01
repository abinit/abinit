// abi_common_kokkos.h

/*
 * Copyright (C) 2022-2024 ABINIT Group
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

// generic template aliases for multidimensional views
template<class value_t>
using AbiView_1d = Kokkos::View<value_t*,                        memory_space, memory_trait >;
template<class value_t>
using AbiView_2d = Kokkos::View<value_t**,   Kokkos::LayoutLeft, memory_space, memory_trait >;
template<class value_t>
using AbiView_3d = Kokkos::View<value_t***,  Kokkos::LayoutLeft, memory_space, memory_trait >;
template<class value_t>
using AbiView_4d = Kokkos::View<value_t****, Kokkos::LayoutLeft, memory_space, memory_trait >;

// generic template aliases for multidimensional const views
template<class value_t>
using AbiView_1d_const = Kokkos::View<const value_t*,                        memory_space, memory_trait >;
template<class value_t>
using AbiView_2d_const = Kokkos::View<const value_t**,   Kokkos::LayoutLeft, memory_space, memory_trait >;
template<class value_t>
using AbiView_3d_const = Kokkos::View<const value_t***,  Kokkos::LayoutLeft, memory_space, memory_trait >;
template<class value_t>
using AbiView_4d_const = Kokkos::View<const value_t****, Kokkos::LayoutLeft, memory_space, memory_trait >;

//! device views types for array of doubles, unmanaged memory (allocated on fortran side)
using AbiView_r64_1d = AbiView_1d<double>;
using AbiView_r64_2d = AbiView_2d<double>;
using AbiView_r64_3d = AbiView_3d<double>;
using AbiView_r64_4d = AbiView_4d<double>;

//! device views types for array of doubles, unmanaged memory (allocated on fortran side), const
using AbiView_r64_1d_const = AbiView_1d_const<double>;
using AbiView_r64_2d_const = AbiView_2d_const<double>;
using AbiView_r64_3d_const = AbiView_3d_const<double>;
using AbiView_r64_4d_const = AbiView_4d_const<double>;

//! device views types for array of doubles, managed memory
using AbiView_r64_3d_managed = Kokkos::View<double***, Kokkos::LayoutLeft,  memory_space >;


//! device views types for array of integers, unmanaged memory (allocated on fortran side)
using AbiView_i32_1d = AbiView_1d<int32_t>;
using AbiView_i32_2d = AbiView_2d<int32_t>;
using AbiView_i32_3d = AbiView_3d<int32_t>;
using AbiView_i32_4d = AbiView_4d<int32_t>;

//! device views types for array of complex doubles, unmanaged memory (allocated on fortran side)
using cplx_t = Kokkos::complex<double>;
using AbiView_c64_1d = AbiView_1d<cplx_t>;
using AbiView_c64_2d = AbiView_2d<cplx_t>;
using AbiView_c64_3d = AbiView_3d<cplx_t>;
using AbiView_c64_4d = AbiView_4d<cplx_t>;

//! device views types for array of complex doubles, unmanaged memory (allocated on fortran side), const
using AbiView_c64_1d_const = AbiView_1d_const<cplx_t>;
using AbiView_c64_2d_const = AbiView_2d_const<cplx_t>;
using AbiView_c64_3d_const = AbiView_3d_const<cplx_t>;
using AbiView_c64_4d_const = AbiView_4d_const<cplx_t>;

#endif // _ABINIT_COMMON_KOKKOS_H
