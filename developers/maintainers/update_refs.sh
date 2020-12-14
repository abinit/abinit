#!/usr/bin/env bash
# Copyright (C) 2019-2020 ABINIT group (XG)
# 
# The purpose of this script is to update all reference files, by brute force.

#Issue this script on the main reference machine, in the directory tests.
#  ../dev*/main*/update_refs.sh
#
# Then check whether all fldiff files have been erased (provided the correspondence below has been correctly implemented).
#
# Normal directories
cp TestBot_MPI1/v1_t*/*out  v1/Refs
cp TestBot_MPI1/v2_t*/*out  v2/Refs
cp TestBot_MPI1/v3_t*/*out TestBot_MPI1/v3_t*/*DOS*                         v3/Refs
cp TestBot_MPI1/v3_t*/*_k*_b* TestBot_MPI1/v3_t*/t79_*                      v3/Refs
cp TestBot_MPI1/v4_t*/*out* TestBot_MPI1/v4_t*/*DOS* TestBot_MPI1/v4_t*/*AV  v4/Refs
cp TestBot_MPI1/v5_t*/*out TestBot_MPI1/v5_t*/*DOS* TestBot_MPI1/v5_t*/*xml                         v5/Refs
cp TestBot_MPI1/v5_t*/t10_*  TestBot_MPI1/v5_t*/*FATBANDS*              v5/Refs
cp TestBot_MPI1/v6_t*/*out TestBot_MPI1/v6_t*/*B*TR*P*  v6/Refs
cp TestBot_MPI1/v7_t*/*out TestBot_MPI1/v7_t*/*DOS*     v7/Refs
cp TestBot_MPI1/v7_t*/*DMFT TestBot_MPI1/v7_t*/*SBK  TestBot_MPI1/v7_t*/t88* v7/Refs
cp TestBot_MPI1/v8_t*/*out* TestBot_MPI1/v8_t*/*dat* TestBot_MPI1/v8_t*/*xml TestBot_MPI1/v8_t*/t58* TestBot_MPI1/v8_t*/*DDB  v8/Refs
cp TestBot_MPI1/v9_t*/*out* TestBot_MPI1/v9_t*/*KERANGE   v9/Refs
cp TestBot_MPI1/v67mbpt_t*/*out  TestBot_MPI1/v67mbpt_t*/*MDF TestBot_MPI1/v67mbpt_t*/*_k*_b*       v67mbpt/Refs
cp TestBot_MPI1/atompaw_t*/*out  atompaw/Refs
cp TestBot_MPI1/bigdft_t*/*out   bigdft/Refs
cp TestBot_MPI1/etsf_io_t*/*out  etsf_io/Refs
cp TestBot_MPI1/fast_t*/*out  TestBot_MPI1/fast_t*/*GEO     fast/Refs
cp TestBot_MPI1/libxc_t*/*out TestBot_MPI1/libxc_t*/*data   libxc/Refs
cp TestBot_MPI1/psml_t*/*out     psml/Refs
cp TestBot_MPI1/tutomultibinit_t*/*out TestBot_MPI1/tutomultibinit_t*/*xml tutomultibinit/Refs
cp TestBot_MPI1/tutoplugs_t*/*out      tutoplugs/Refs
cp TestBot_MPI1/tutorespfn_t*/*out TestBot_MPI1/tutorespfn_t*/*MRTA* TestBot_MPI1/tutorespfn_t*/*SERTA*    tutorespfn/Refs
cp TestBot_MPI1/tutorial_t*/*out  TestBot_MPI1/tutorial_t*/*f2b  TestBot_MPI1/tutorial_t*/*abo  tutorial/Refs
cp TestBot_MPI1/unitary_t*/*out        unitary/Refs
cp TestBot_MPI1/vdwxc_t*/*out          vdwxc/Refs
cp TestBot_MPI1/wannier90_t*/*out      wannier90/Refs

# Directories for parallel tests
cp TestBot_MPI*/bigdft_paral_t*/*out   bigdft_paral/Refs
cp TestBot_MPI*/mpiio_t*/*out TestBot_MPI*/mpiio_t*/*DOS*        mpiio/Refs
cp TestBot_MPI*/paral_t*/*out TestBot_MPI*/paral_t*/*MDF   TestBot_MPI*/paral_t*/*xml    TestBot_MPI*/paral_t*/*dat  TestBot_MPI*/paral_t*/*DDB paral/Refs
cp TestBot_MPI*/tutoparal_t*/*out      tutoparal/Refs
cp TestBot_MPI*/unitary_t*/*out        unitary/Refs

# For the serial reference machine
#p TestBot_MPI1/seq_t*/*out      seq/Refs

# For the gpu reference machine
#cp TestBot_MPI*/gpu_t*/*out      gpu/Refs

# Now, suppress the corresponding fldiff files, to see whether some relevant files have not been copied
# Normal directories
rm TestBot_MPI1/v1_t*/*out.fldiff
rm TestBot_MPI1/v2_t*/*out.fldiff  
rm TestBot_MPI1/v3_t*/*out.fldiff TestBot_MPI1/v3_t*/*DOS*.fldiff
rm TestBot_MPI1/v3_t*/*_k*_b*.fldiff TestBot_MPI1/v3_t*/t79_*.fldiff
rm TestBot_MPI1/v4_t*/*out*.fldiff TestBot_MPI1/v4_t*/*DOS*.fldiff TestBot_MPI1/v4_t*/*AV.fldiff 
rm TestBot_MPI1/v5_t*/*out.fldiff TestBot_MPI1/v5_t*/*DOS*.fldiff TestBot_MPI1/v5_t*/*xml.fldiff                   
rm TestBot_MPI1/v5_t*/t10_*.fldiff  TestBot_MPI1/v5_t*/*FATBANDS*.fldiff   
rm TestBot_MPI1/v6_t*/*out.fldiff TestBot_MPI1/v6_t*/*B*TR*P*.fldiff 
rm TestBot_MPI1/v7_t*/*out.fldiff TestBot_MPI1/v7_t*/*DOS*.fldiff
rm TestBot_MPI1/v7_t*/*DMFT.fldiff TestBot_MPI1/v7_t*/*SBK.fldiff  TestBot_MPI1/v7_t*/t88*.fldiff
rm TestBot_MPI1/v8_t*/*out*.fldiff TestBot_MPI1/v8_t*/*dat*.fldiff TestBot_MPI1/v8_t*/*xml.fldiff TestBot_MPI1/v8_t*/t58*.fldiff TestBot_MPI1/v8_t*/*DDB.fldiff
rm TestBot_MPI1/v9_t*/*out*.fldiff TestBot_MPI1/v9_t*/*KERANGE.fldiff  
rm TestBot_MPI1/v67mbpt_t*/*out.fldiff  TestBot_MPI1/v67mbpt_t*/*MDF.fldiff TestBot_MPI1/v67mbpt_t*/*_k*_b*.fldiff
rm TestBot_MPI1/atompaw_t*/*out.fldiff  
rm TestBot_MPI1/bigdft_t*/*out.fldiff
rm TestBot_MPI1/etsf_io_t*/*out.fldiff
rm TestBot_MPI1/fast_t*/*out.fldiff  TestBot_MPI1/fast_t*/*GEO.fldiff
rm TestBot_MPI1/libxc_t*/*out.fldiff TestBot_MPI1/libxc_t*/*data.fldiff
rm TestBot_MPI1/psml_t*/*out.fldiff   
rm TestBot_MPI1/tutomultibinit_t*/*out.fldiff TestBot_MPI1/tutomultibinit_t*/*xml.fldiff
rm TestBot_MPI1/tutoplugs_t*/*out.fldiff
rm TestBot_MPI1/tutorespfn_t*/*out.fldiff TestBot_MPI1/tutorespfn_t*/*MRTA*.fldiff TestBot_MPI1/tutorespfn_t*/*SERTA*.fldiff
rm TestBot_MPI1/tutorial_t*/*out.fldiff  TestBot_MPI1/tutorial_t*/*f2b.fldiff TestBot_MPI1/tutorial_t*/*abo.fldiff
rm TestBot_MPI1/unitary_t*/*out.fldiff
rm TestBot_MPI1/vdwxc_t*/*out.fldiff
rm TestBot_MPI1/wannier90_t*/*out.fldiff

# Directories for parallel tests
rm TestBot_MPI*/bigdft_paral_t*/*out.fldiff
rm TestBot_MPI*/mpiio_t*/*out.fldiff TestBot_MPI*/mpiio_t*/*DOS*.fldiff
rm TestBot_MPI*/paral_t*/*out.fldiff TestBot_MPI*/paral_t*/*MDF.fldiff   TestBot_MPI*/paral_t*/*xml.fldiff    TestBot_MPI*/paral_t*/*dat.fldiff  TestBot_MPI*/paral_t*/*DDB.fldiff 
rm TestBot_MPI*/tutoparal_t*/*out.fldiff  
rm TestBot_MPI*/unitary_t*/*out.fldiff

# For the serial reference machine
#rm TestBot_MPI1/seq_t*/*out.fldiff

# For the gpu reference machine
#rm TestBot_MPI*/gpu_t*/*out.fldiff      

