/****************************************/
/*** Translated to C++ & parallelized ***/
/*** in OpenCL by Cosmin E. Oancea    ***/
/***   from an original OCaml code    ***/
/***    (sequential) by LexiFi        ***/
/***   and a sequential Python implem ***/
/***    using a different genetic alg.***/
/***    by Christian Andreetta        ***/  
/****************************************/ 


#include <stdio.h>
#include <stdlib.h> 
#include <math.h>   
#include <algorithm>  
#include <assert.h>

#include "Constants.h"
#include "Date.h"
#include "MathModule.h"
#include "G2ppUtil.h"
#include "G2PP.h"

#include "GenAlgUtil.h"
#include "Candidate.h"
#include "GenAlgFlat.h"

#include "FindBestKer.h"

void test_all() {
    test_dates             ();
    test_math              ();
    test_g2ppUtil          ();
    test_pricer_of_swaption();
}


int main() {
    //test_all();

    //pricer  ();

//    int  N;
//    bool adjustable;
//    getIregShapeAdjusted(  64, N );
//    getIregShapeAdjusted( 128, N );
//    getIregShapeAdjusted( 256, N );
//    getIregShapeAdjusted( 512, N );
//    return 1;

    printf("\n\n\n\nNow Running the Genetic Algorithm to Find the Best Estimation!\n\n");

    {
//        mlfi_timeb  t_start, t_end;
//        unsigned long int elapsed;
//        mlfi_ftime(&t_start);

        mainKernel (  );

//        mlfi_ftime(&t_end);
//        elapsed = mlfi_diff_time(t_end,t_start);
//        printf("\n\nTotal Application Calibration Time: %lu !\n\n", elapsed);

    }
    return 1;
}
