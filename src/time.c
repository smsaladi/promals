#include "time.h"
#include "all.h"

struct tms      tmsstart, tmsend;

//double pr_times(struct tms *tmsstart, struct tms *tmsend)
double timeDiff()
{
        static long             clktck = 0;
	double CPU_time;

        if (clktck == 0)        /* fetch clock ticks per second first time */
                if ( (clktck = sysconf(_SC_CLK_TCK)) < 0)
                        fprintf(stderr,"sysconf error");
	CPU_time = (tmsend.tms_utime - tmsstart.tms_utime)/(double) clktck + (tmsend.tms_stime - tmsstart.tms_stime) / (double) clktck;
        //fprintf(stdout, "CPU time: %7.2f\n", (tmsend->tms_cutime - tmsstart->tms_cutime) / (double) clktck + (tmsend->tms_cstime - tmsstart->tms_cstime) / (double) clktck );
	fprintf(stdout, "CPU time elapsed: %f\n", CPU_time); fflush(stdout);
	return CPU_time;
}
