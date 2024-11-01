#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>

/*
 * Set stack size limit to maximum allowed value. Return soft and hard limit and exit status.
 * */

void ulimit_stack(long int *rlim_cur, long int *rlim_max, int *ierr) {

  int result;
  struct rlimit rlim = {RLIM_INFINITY, RLIM_INFINITY};
  /* From https://linux.die.net/man/2/getrlimit
  struct rlimit {
      rlim_t rlim_cur;  Soft limit
      rlim_t rlim_max;  Hard limit (ceiling for rlim_cur)
  }; */

  /* Try to set both soft and hard limits to INFINITY */
  result = setrlimit(RLIMIT_STACK, &rlim);
  if (result == -1) {
    /* Try to set soft == sys hard */
    if (getrlimit(RLIMIT_STACK, &rlim) == 0) {
        rlim.rlim_cur = rlim.rlim_max;
        result = setrlimit(RLIMIT_STACK, &rlim);
    }
  }

  *ierr = result;
  *rlim_cur = -1;
  *rlim_max = -1;
  if (getrlimit(RLIMIT_STACK, &rlim) == 0) {
    *rlim_cur = rlim.rlim_cur;
    *rlim_max = rlim.rlim_max;
  }

}
