#
# 52_fft_mpi_noabirule
#
add_library(52_fft_mpi_noabirule STATIC
  defs_fftdata.F90
  m_fftw3.F90
  m_dfti.F90
  m_fftcore.F90
  m_fft_mesh.F90
  m_sgfft.F90
  m_sg2002.F90
  )

set_target_properties(52_fft_mpi_noabirule
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(52_fft_mpi_noabirule
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(52_fft_mpi_noabirule
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(52_fft_mpi_noabirule
  PUBLIC
  abinit::44_abitools
  abinit::51_manage_mpi)

if(OpenMP_FOUND)
  target_link_libraries(52_fft_mpi_noabirule
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

if (ABINIT_FFT_FLAVOR MATCHES "GOEDECKER")
  # no linking required
else()
  target_link_libraries(52_fft_mpi_noabirule
    PUBLIC
    abinit::fft)
endif()

add_library(abinit::52_fft_mpi_noabirule ALIAS 52_fft_mpi_noabirule)
