#ifndef TIME_HELPER_CHR
#define TIME_HELPER_CHR

#include <sys/time.h>
#include <time.h>



#if 0
    unsigned long int elapsed;
    mlfi_timeb  t_start, t_end;
    mlfi_ftime(&t_start);
    // ...
    mlfi_ftime(&t_end);
    elapsed = mlfi_diff_time(t_end,t_start);
#endif


#endif
