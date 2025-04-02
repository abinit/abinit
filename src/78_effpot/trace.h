#define DMSG(msg) msg=''; write(msg,'(2A,I0)') TRIM(__FILE__)//' line: ', '', __LINE__
#define DMSG2(msg) msg=TRIM(__FILE__); print *, trim(msg), __LINE__

