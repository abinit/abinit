#
# 14_hidewrite
#

# generate m_build_info.F90
include(${CMAKE_CURRENT_SOURCE_DIR}/generate_build_info.cmake)
generate_build_info()

# TODO : refactor, there is currently no systematic way to
# to override compilation flags; we should write a cmake macro
# to todo that using set/get_target_properties for property
# "COMPILE_FLAGS"
configure_file(m_optim_dumper.F90.cmake.in m_optim_dumper.F90 @ONLY)

add_library(14_hidewrite STATIC
  m_abicore.F90
  m_cppopts_dumper.F90
  m_fstrings.F90
  m_io_tools.F90
  m_specialmsg.F90
  # generated sources
  ${CMAKE_CURRENT_BINARY_DIR}/m_build_info.F90
  ${CMAKE_CURRENT_BINARY_DIR}/m_optim_dumper.F90)

set_target_properties(14_hidewrite
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(14_hidewrite
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  )

target_include_directories(14_hidewrite
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
  $<INSTALL_INTERFACE:include/libpaw>
  )

target_compile_definitions(14_hidewrite
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(14_hidewrite
  PUBLIC
  abinit::10_defs
  abinit::11_memory_mpi
  abinit::12_hide_mpi
  MPI::MPI_Fortran)

add_library(abinit::14_hidewrite ALIAS 14_hidewrite)
