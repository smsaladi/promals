#include "all.h"
#include "time.h"

struct tms tmsstart, tmsend;

// double pr_times(struct tms *tmsstart, struct tms *tmsend)
double timeDiff() {
  static long clktck = 0;
  double CPU_time;

  if (clktck == 0) /* fetch clock ticks per second first time */
    if ((clktck = sysconf(_SC_CLK_TCK)) < 0) fprintf(stderr, "sysconf error");
  CPU_time = (tmsend.tms_utime - tmsstart.tms_utime) / (double)clktck +
             (tmsend.tms_stime - tmsstart.tms_stime) / (double)clktck;
  // fprintf(stdout, "CPU time: %7.2f\n", (tmsend->tms_cutime -
  // tmsstart->tms_cutime) / (double) clktck + (tmsend->tms_cstime -
  // tmsstart->tms_cstime) / (double) clktck );
  fprintf(stdout, "CPU time elapsed: %f\n", CPU_time);
  fflush(stdout);
  return CPU_time;
}

time_t timestart, timeend;
double timediff(time_t &s, time_t &e) {
  double diff = difftime(e, s);
  return diff;
}

char prev_info[500];
int timecount;
void print_time_diff(const char *info) {
  time(&timeend);
  struct tm *timeinfo = localtime(&timeend);
  double diff = timediff(timestart, timeend);
  if (timecount == 0) strcpy(prev_info, "beginning");
  cout << "------------" << endl;
  printf("Current time and date: %s", asctime(timeinfo));
  cout << "Time between '" << info << "' and '" << prev_info << "': ";
  cout << diff << endl;
  cout << "----" << endl;
  timestart = timeend;
  timecount++;
  strcpy(prev_info, info);
}
void print_section_info(const char *info) {
  cout << "=====================" << endl;
  cout << info << endl;
  cout << endl;
}
