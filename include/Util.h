#ifndef GENERIC_UTILITIES
#define GENERIC_UTILITIES

/*******************************************************/
/*****  Utilities Related to Time Instrumentation  *****/
/*******************************************************/
//#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/timeb.h>
#include <assert.h>
#include <omp.h>

typedef struct timeb mlfi_timeb;
#define mlfi_ftime ftime

#define mlfi_diff_time(t1,t2) \
  (t1.time - t2.time) * 1000 + (t1.millitm - t2.millitm)

/*******************************************************/
/*****   Utilities Related to Read/Write to File   *****/
/*******************************************************/

int get_tot_num_threads() {
    int procs;

#pragma omp parallel shared(procs)
    {
        int th_id = omp_get_thread_num();
        if(th_id == 0) { procs = omp_get_num_threads(); }
    }

    bool valid_procs = (procs > 0) && (procs <= 1024);
    assert(valid_procs && "Number of threads NOT in {1, ..., 1024}");
    return procs;
}

#endif //GENERIC_UTILITIES
