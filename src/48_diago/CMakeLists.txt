#
# 48_diago
#
add_library(48_diago STATIC
  m_chebfi2.F90
  m_lobpcg2.F90
  )

set_target_properties(48_diago
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(48_diago
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  ${CMAKE_SOURCE_DIR}/src/incs
  )

target_compile_definitions(48_diago
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(48_diago
  PUBLIC
  abinit::44_abitools
  abinit::44_abitypes_defs
  abinit::45_xgTools)

if(ABINIT_ENABLE_GPU_CUDA)
  target_link_libraries(48_diago
    PUBLIC
    abinit::46_manage_cuda
    )
endif()

add_library(abinit::48_diago ALIAS 48_diago)
