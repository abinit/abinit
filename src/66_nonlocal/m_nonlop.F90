#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_nonlop
  use defs_basis
  use defs_abitypes
  use m_errors
  use m_profiling_abi
  use m_gemm_nonlop
  use m_hamiltonian, only :
  use m_pawcprj,     only : pawcprj_type,pawcprj_alloc,pawcprj_free,pawcprj_copy
  use m_xmpi
  use m_cgtools
  use m_opernl
  use m_opernl_ylm, only : opernla_ylm, opernlb_ylm, opernlc_ylm, opernld_ylm

  implicit none

contains

#include "nonlop.finc"
#include "nonlop_gpu.finc"
#include "nonlop_pl.finc"
#include "nonlop_ylm.finc"
end module m_nonlop
