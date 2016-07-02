#include <stdlib.h>
#include <stdio.h>

int main(void)
{
#if defined __GNUC__
#if defined __i386
#define arch "x86_32"
#elif defined __ia64
#define arch "ia64"
#else
#define arch "unknown"
#endif
#if defined __linux__
#define system "linux"
#elif defined __macosx__
#define system "macosx"
#else
#define system "unknown"
#endif
 printf("gnu %d.%d %s %s\n",__GNUC__,__GNUC_MINOR__,arch,system);
#endif

 return(0);
}
