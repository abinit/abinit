#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_opernl
  use m_errors
  use m_profiling_abi

  implicit none

contains 

#include "opernl2.finc"
#include "opernl3.finc"
#include "opernl4a.finc"
#include "opernl4b.finc"
end module m_opernl
