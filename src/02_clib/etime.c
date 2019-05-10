#include "abi_clib.h"

#if defined HAVE_SYS_TIME_H && defined HAVE_SYS_RESOURCE_H
#include <sys/time.h>
#include <sys/resource.h>
#endif

double etime(tt)
float tt[2];
{
#if defined HAVE_SYS_TIME_H && defined HAVE_SYS_RESOURCE_H
  int who;
  struct rusage used;
  who = 0;
  getrusage(who,&used);
  tt[0] = used.ru_utime.tv_sec+((used.ru_utime.tv_usec)/1000000.);
  tt[1] = used.ru_stime.tv_sec+((used.ru_stime.tv_usec)/1000000.);
  return(tt[0]+tt[1]);
#else
  tt[0]=-1;
  tt[1]=-1;
  return(tt[0]+tt[1]);
#endif
}
