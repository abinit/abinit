#
# 02_clib
#
add_library(02_clib STATIC
  calling_levmar.c
  cclock.c
  etime.c
  sockets.c
  fsi_posix.c
  gnu_tools.c
  mallinfo.c
  md5.c
  md5.h
  rlimit.c
  xexit.c
  xmalloc.c
  xmalloc.h
  m_fsockets.F90
  m_clib.F90)

set_target_properties(02_clib
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(02_clib
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  )

target_include_directories(02_clib
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
  $<INSTALL_INTERFACE:include/libpaw>
  )

target_compile_definitions(02_clib
  PRIVATE
  HAVE_CONFIG_H)

if(levmar_FOUND)
  target_link_libraries(02_clib
    PUBLIC
    abinit::levmar)
endif()

add_library(abinit::02_clib ALIAS 02_clib)
