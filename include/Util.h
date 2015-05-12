#ifndef GENERIC_UTILITIES
#define GENERIC_UTILITIES

/*******************************************************/
/*****  Utilities Related to Time Instrumentation  *****/
/*******************************************************/
//#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
//#include <sys/timeb.h>
#include <assert.h>

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif


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

#if 0
typedef struct timeb mlfi_timeb;
#define mlfi_ftime ftime

#define mlfi_diff_time(t1,t2) \
  (t1.time - t2.time) * 1000 + (t1.millitm - t2.millitm)
#endif
/*******************************************************/
/*****   Utilities Related to Read/Write to File   *****/
/*******************************************************/

int get_CPU_num_threads() {
    int procs;

#ifdef ENABLE_OPENMP
#pragma omp parallel shared(procs)
    {
        int th_id = omp_get_thread_num();
        if(th_id == 0) { procs = omp_get_num_threads(); }
    }

    bool valid_procs = (procs > 0) && (procs <= 1024);
    assert(valid_procs && "Number of threads NOT in {1, ..., 1024}");
#else
    procs = 1;
#endif
    return procs;
}

//int get_GPU_num_threads() {
//    return 1024;
//}

bool is_pow2(int atr_val) {
    int x = 1;

    for(int i = 0; i < 31; i++) {
        if(x == atr_val) return true;
        x = (x << 1);
    }
    return false;
}

#endif //GENERIC_UTILITIES
