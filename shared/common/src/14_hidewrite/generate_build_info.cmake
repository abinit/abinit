macro(generate_build_info)

  set(TARGET_CPU ${CMAKE_HOST_SYSTEM_PROCESSOR})
  cmake_host_system_information(RESULT TARGET_CPU_LONG QUERY PROCESSOR_DESCRIPTION)
  cmake_host_system_information(RESULT OS_NAME QUERY OS_NAME)

  set(abi_target_os   ${OS_NAME})
  set(abi_cc_vendor   ${CMAKE_C_COMPILER_ID})
  set(abi_cc_version  ${CMAKE_C_COMPILER_VERSION})
  set(abi_cxx_vendor  ${CMAKE_CXX_COMPILER_ID})
  set(abi_cxx_version ${CMAKE_CXX_COMPILER_VERSION})
  set(abi_fc_vendor   ${CMAKE_Fortran_COMPILER_ID})
  set(abi_fc_version  ${CMAKE_Fortran_COMPILER_VERSION})

  set(ABINIT_TARGET ${TARGET_CPU}_${abi_target_os}_${abi_fc_vendor}_${abi_fc_version})

  message(STATUS "Generating m_build_info.F90...")
  configure_file(m_build_info.F90.cmake.in m_build_info.F90 @ONLY)

endmacro()
