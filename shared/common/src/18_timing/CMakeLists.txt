#
# 18_timing
#
add_library(18_timing STATIC
  m_xpapi.F90
  m_time.F90
  )

set_target_properties(18_timing
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(18_timing
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  )

target_include_directories(18_timing
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
  $<INSTALL_INTERFACE:include/libpaw>
  )

target_compile_definitions(18_timing
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(18_timing
  PUBLIC
  abinit::16_hideleave
  )

if(ABINIT_PAPI_FOUND)
  target_link_libraries(18_timing
    PUBLIC
    abinit::papi)
endif()

add_library(abinit::18_timing ALIAS 18_timing)
