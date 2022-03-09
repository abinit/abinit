#
# 51_manage_mpi
#
add_library(51_manage_mpi STATIC
  m_io_redirect.F90
  m_mpinfo.F90
  m_wffile.F90
  m_xredistribute.F90
  )

set_target_properties(51_manage_mpi
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(51_manage_mpi
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(51_manage_mpi
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(51_manage_mpi
  PUBLIC
  abinit::44_abitypes_defs)

add_library(abinit::51_manage_mpi ALIAS 51_manage_mpi)
