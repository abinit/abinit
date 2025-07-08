#
# 17_minimax
#
add_library(17_minimax STATIC
  minimax_grids.F90
  minimax_omega.F90
  minimax_tau.F90
  minimax_utils.F90
  )

set_target_properties(17_minimax
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

target_include_directories(17_minimax
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  )

target_include_directories(17_minimax
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
  $<INSTALL_INTERFACE:include/libpaw>
  )

target_compile_definitions(17_minimax
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(17_minimax
  PUBLIC
  abinit::16_hideleave
  )

add_library(abinit::17_minimax ALIAS 17_minimax)
