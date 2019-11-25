#!/usr/bin/env bash
# Copyright (C) 2019-2020 ABINIT group (XG)
# 
# The purpose of this script is to update all reference files, by brute force.

#Issue this script on the main reference machine, in the directory tests.
#  ../dev*/main*
#
# Normal directories
cp TestBot_MPI1/v1_t*/*out TestBot_MPI1/v1_t*/*DOS* TestBot_MPI1/v1_t*/*cml v1/Refs
cp TestBot_MPI1/v2_t*/*out                                                  v2/Refs
cp TestBot_MPI1/v3_t*/*out TestBot_MPI1/v3_t*/*DOS*                         v3/Refs
cp TestBot_MPI1/v3_t*/*_k*_b* TestBot_MPI1/v3_t*/t79_*                      v3/Refs
cp TestBot_MPI1/v4_t*/*out TestBot_MPI1/v4_t*/*DOS* TestBot_MPI1/v4_t*/*AV TestBot_MPI1/v4_t*/*DDB  v4/Refs
cp TestBot_MPI1/v5_t*/*out TestBot_MPI1/v5_t*/*DOS* TestBot_MPI1/v5_t*/*xml                         v5/Refs
cp TestBot_MPI1/v5_t*/*PROCAR* TestBot_MPI1/v5_t*/t10_*  TestBot_MPI1/v5_t*/*FATBANDS*              v5/Refs
cp TestBot_MPI1/v6_t*/*out TestBot_MPI1/v6_t*/*BLZTRP*  v6/Refs
cp TestBot_MPI1/v7_t*/*out TestBot_MPI1/v7_t*/*DOS* TestBot_MPI1/v7_t*/*dat v7/Refs
cp TestBot_MPI1/v7_t*/*DMFT TestBot_MPI1/v7_t*/*TETRA TestBot_MPI1/v7_t*/t56_* v7/Refs
cp TestBot_MPI1/v7_t*/t84_* TestBot_MPI1/v7_t*/t88_* v7/Refs
cp TestBot_MPI1/v8_t*/*out TestBot_MPI1/v8_t*/*dat* TestBot_MPI1/v8_t*/*xml v8/Refs
cp TestBot_MPI1/v9_t*/*out TestBot_MPI1/v9_t*/*KERANGE   v9/Refs
cp TestBot_MPI1/v67mbpt_t*/*out  TestBot_MPI1/v67mbpt_t*/*MDF TestBot_MPI1/v67mbpt_t*/*_k*_b*       v67mbpt/Refs
cp TestBot_MPI1/atompaw_t*/*out  atompaw/Refs
cp TestBot_MPI1/bigdft_t*/*out   bigdft/Refs
cp TestBot_MPI1/etsf_io_t*/*out  etsf_io/Refs
cp TestBot_MPI1/fast_t*/*out  TestBot_MPI1/fast_t*/*GEO     fast/Refs
cp TestBot_MPI1/libxc_t*/*out TestBot_MPI1/libxc_t*/*data   libxc/Refs
cp TestBot_MPI1/psml_t*/*out     psml/Refs
cp TestBot_MPI1/tutomultibinit_t*/*out TestBot_MPI1/tutomultibinit_t*/*xml tutomultibinit/Refs
cp TestBot_MPI1/tutoplugs_t*/*out      tutoplugs/Refs
cp TestBot_MPI1/tutorespfn_t*/*out     tutorespfn/Refs
cp TestBot_MPI1/tutorial_t*/*out  TestBot_MPI1/tutorial_t*/*f2b   tutorial/Refs
cp TestBot_MPI1/unitary_t*/*out        unitary/Refs
cp TestBot_MPI1/vdwxc_t*/*out          vdwxc/Refs
cp TestBot_MPI1/wannier90_t*/*out      wannier90/Refs

# Directories for parallel tests
cp TestBot_MPI*/bigdft_paral_t*/*out   bigdft_paral/Refs
cp TestBot_MPI*/mpiio_t*/*out TestBot_MPI*/mpiio_t*/*DOS*  TestBot_MPI*/mpiio_t*/*VCMLB*       mpiio/Refs
cp TestBot_MPI*/paral_t*/*out TestBot_MPI*/paral_t*/*MDF   TestBot_MPI*/paral_t*/*xml    TestBot_MPI*/paral_t*/*dat  TestBot_MPI*/paral_t*/*DDB paral/Refs
cp TestBot_MPI*/tutoparal_t*/*out TestBot_MPI*/tutoparal_t*/*PRESS*     tutoparal/Refs
cp TestBot_MPI*/unitary_t*/*out        unitary/Refs

# For the serial reference machine
#cp TestBot_MPI1/seq_t*/*out      seq/Refs

# For the gpu reference machine
#cp TestBot_MPI*/gpu_t*/*out      gpu/Refs

