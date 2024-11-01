#
# 57_iopsp_parser
#
add_library(57_iopsp_parser STATIC
  flib_pwscf.F90
  functionals_pwscf.F90
  m_pseudo_pwscf.F90
  m_pspheads.F90
  m_psxml2ab.F90
  m_read_upf_pwscf.F90
  # modules taken from QE to support UPF2
  pseudo_types.F90
  read_upf_new.F90
  upf_kinds.F90
  xmltools.F90
  upf_utils.F90
  )

set_target_properties(57_iopsp_parser
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(57_iopsp_parser
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(57_iopsp_parser
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(57_iopsp_parser
  PUBLIC
  abinit::21_hashfuncs
  abinit::44_abitypes_defs
  )

add_library(abinit::57_iopsp_parser ALIAS 57_iopsp_parser)
