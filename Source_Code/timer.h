//From Paul Bauman as provided in HPC 1 at UB

#include <time.h>
#include <sys/time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

void diff( struct timespec* start, struct timespec* end, struct timespec* diff);

void get_time( struct timespec* time );
