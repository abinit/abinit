#
# 67_python_invocation_ext
#
add_library(67_python_invocation_ext STATIC
  execute_python.c
  m_invocation_tools.F90
  m_invoke_python.F90
  )

set_target_properties(67_python_invocation_ext
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

target_include_directories(67_python_invocation_ext
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_BINARY_DIR}/modules
  )

target_compile_definitions(67_python_invocation_ext
  PRIVATE
  HAVE_CONFIG_H)

#target_link_libraries(67_python_invocation_ext
#  PUBLIC
#  abinit::10_defs
#  abinit::16_hideleave
#  )

add_library(abinit::67_python_invocation_ext ALIAS 67_python_invocation_ext)
