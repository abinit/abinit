
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBTETRA_ABINIT
#include "abi_common.h"
#  define TETRA_ALLOCATE(ARR,SIZE) ABI_ALLOCATE(ARR,SIZE)
#  define TETRA_DEALLOCATE(ARR)  ABI_DEALLOCATE(ARR)
#  define TETRA_ERROR(MSG) MSG_ERROR(MSG)
#  define USE_MEMORY_PROFILING use m_profiling_abi
#  define USE_MSG_HANDLING use m_errors, only : msg_hndl

#else
#  define TETRA_ALLOCATE(ARR,SIZE) allocate(ARR SIZE)
#  define TETRA_DEALLOCATE(ARR)  deallocate(ARR)
#  define TETRA_ERROR(MSG) write (stdout,'(a)') MSG ; stop "ERROR"
#  define USE_MEMORY_PROFILING 
#  define USE_MSG_HANDLING 
#endif

