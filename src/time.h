#ifndef __time___
#define __time___

#include <sys/times.h>
#include <sys/types.h>
#include <unistd.h>
#include <cstdio>
#include <ctime>

extern struct tms tmsstart, tmsend;

double timeDiff();

extern time_t timestart, timeend;
extern double timediff(time_t &s, time_t &e);
extern char prev_info[500];
extern int timecount;
extern void print_time_diff(const char *info);
extern void print_section_info(const char *info);

#endif
