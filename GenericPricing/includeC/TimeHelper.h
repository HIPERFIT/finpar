#ifndef TIME_HELPER_CHR
#define TIME_HELPER_CHR

#include <sys/time.h>
#include <time.h>

#define TIME_RESOLUTION_MICROSECOND

#ifdef TIME_RESOLUTION_MICROSECOND
// CHR: added helper function for computing time differences at microsecond resolution
int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1)
{
    unsigned int resolution=1000000;
    long int diff = (t2->tv_usec + resolution * t2->tv_sec) - (t1->tv_usec + resolution * t1->tv_sec);
    result->tv_sec = diff / resolution;
    result->tv_usec = diff % resolution;
    return (diff<0);
}
#endif

#endif
