#ifndef __time___
#define __time___

#include        <sys/types.h>
#include        <sys/times.h>
#include        <unistd.h>

extern struct tms      tmsstart, tmsend;

double timeDiff();

#endif
