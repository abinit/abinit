#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_opernl_ylm
  use m_profiling_abi
  use m_errors
  use m_xmpi

  implicit none

contains 

#include "opernla_ylm.finc"
#include "opernlb_ylm.finc"
#include "opernlc_ylm.finc"
#include "opernld_ylm.finc"
end module m_opernl_ylm
