# FIXME CMake doesn't manage OpenMP target flags yet (v3.26, 2023/06/09)
# This function is adapted from a proposed solution from related issue on CMake Gitlab:
# https://gitlab.kitware.com/cmake/cmake/-/issues/21784#note_1354500


# Get the OpenMP device offload flags for the current Fortran compiler.
#
#   TARGET <offload target>
#       names the target for offloading (optional).
#
#   ARCH <target architecture>
#       names the architcure to compiler for (optional).
#
#   ADD_FLAGS <other flags>
#       additional flags that may be needed (optional).
#
#   RESULT <avr name>
#       the flags are stored in this variable.
#
function(get_offload_flags)
    set(opts "")
    set(nvpo ARCH TARGET ADD_FLAGS RESULT)
    set(mvo)
    cmake_parse_arguments(PARSE_ARGV 0 OMP_DO "${opts}" "${nvpo}" "${mvo}")
    set(tmp)
    if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "LLVMFlang")
        set(tmp "${tmp} -fopenmp")
        if (OMP_DO_TARGET)
            set(tmp "${tmp} -fopenmp-targets=${OMP_DO_TARGET}")
        endif()
        if (OMP_DO_ARCH)
            set(tmp "${tmp} --offload-arch=${OMP_DO_ARCH}")
        endif()
    elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "IntelLLVM")
        set(tmp "${tmp} -qopenmp")
        if (OMP_DO_TARGET)
            set(tmp "${tmp} -fopenmp-targets=${OMP_DO_TARGET}")
        endif()
    elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
        set(tmp "${tmp} -fopenmp")
        if (OMP_DO_TARGET)
            set(tmp "${tmp} -foffload=${OMP_DO_TARGET}")
        endif()
        if (OMP_DO_ARCH)
            set(tmp "${tmp} --offload-options=${OMP_DO_TARGET}=-march=${OMP_DO_ARCH}")
        endif()
    elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "NVHPC")
        set(tmp "-mp=gpu")
        if (OMP_DO_ARCH)
            set(tmp "${tmp}" "-gpu=${OMP_DO_ARCH}")
        endif()
        if (OMP_DO_TARGET)
          set(tmp "${tmp}" "-gpu=cc${OMP_DO_TARGET}")
        endif()
    endif()
    if (OMP_DO_ADD_FLAGS)
        set(tmp "${tmp}" "${OMP_DO_ADD_FLAGS}")
    endif()
    if ("${tmp}" STREQUAL "")
        message(WARNING "OpenMP offload flags not known for ${CMAKE_Fortran_COMPILER_ID}")
    else()
        message(STATUS "OpenMP offload flags for ${CMAKE_Fortran_COMPILER_ID} are ${tmp}")
    endif()
    set(${OMP_DO_RESULT} ${tmp} PARENT_SCOPE)
endfunction()

