#
# 67_triqs_ext
#
add_library(67_triqs_ext STATIC
  triqs_cthyb_qmc.cpp
  triqs_interface_ctqmc.F90
  )

set_target_properties(67_triqs_ext
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(67_triqs_ext
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(67_triqs_ext
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(67_triqs_ext
 PUBLIC
 abinit::10_defs
 abinit::16_hideleave
 )

if(TRIQS_FOUND)
  target_link_libraries(67_triqs_ext
    PUBLIC
    abinit::triqs)
endif()

add_library(abinit::67_triqs_ext ALIAS 67_triqs_ext)
