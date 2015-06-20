#ifndef MAKE_IRREG_SHAPE
#define MAKE_IRREG_SHAPE

#include "Constants.h"

short* getIregShapeAdjusted( const int   LWG, 
                             UINT&       N 
) {
    bool adjustable = true;
    N               = 0;

    /**
     * Compute the segmented sum:
     *    sgmscan     = scan  ( op +, 0, sizes )
     *    adjustable  = reduce( op &&, True ) . map(\sz->sz <= LWG) sizes 
     */
    for( int i = 0; i < NUM_SWAP_QUOTES; i++ ) {
        int sz = 12 * SwaptionQuotes[4*i+2] / SwaptionQuotes[4*i+1];
        N += sz;
        adjustable = adjustable && (sz <= LWG) && (sz > 0);
    }

    assert(adjustable && "Irregular Array Not Adjustable!");

    /**
     * Compute how many local workgroups need adjustements,
     *   and the total number of iddle threads
     */
    {
        // num of adjusted segments: NUM_SWAP_QUOTES + ( (N % LWG == 0) ? N/LWG : N/LWG + 1  
        int  Nbeg           = 0;
        int  Ncur           = 0;
        int  tot_iddle_thds = 0;
        for( int i = 0; i < NUM_SWAP_QUOTES; i++ ) {
            int sz = 12 * SwaptionQuotes[4*i+2] / SwaptionQuotes[4*i+1];
            if( Ncur + sz < Nbeg+LWG ) {
                Ncur += sz;
            } else if (Ncur + sz == Nbeg + LWG ) {
                Nbeg += LWG; 
                Ncur  = Nbeg;
            } else {
                tot_iddle_thds += Nbeg + LWG - Ncur;
                Nbeg += LWG;
                Ncur = Nbeg + sz;
            }
        }

        // adjust for the last workgroup.
        int mmod = Ncur % LWG;
        if( mmod != 0 ) {
            Ncur           += (LWG - mmod);
            tot_iddle_thds += (LWG - mmod); 
        }
        N = Ncur;

        REAL overhead = ((REAL)Ncur)/(Ncur-tot_iddle_thds);
        overhead = (overhead - 1.0) * 100.0;

//        printf("Summary: Ncurr: %d, tot_iddle_thds: %d, overhead: %f percent\n", 
//                Ncur, tot_iddle_thds,  overhead );

        if ( overhead > 20.0 ) {
            adjustable = false;
            return NULL;
        } 
    }    

    short* flags = new short[4*N];
    for( int i = 0; i < 4*N; i++ ) {
        flags[i] = 0;
    }

    {
        // I need:
        //    an array that records the indexes of the iddle threads,
        //    an array that records the indexes of the start of the segments.
        int  Nbeg           = 0;
        int  Ncur           = 0;
        for( int i = 0; i < NUM_SWAP_QUOTES; i++ ) {
            int sz = 12 * SwaptionQuotes[4*i+2] / SwaptionQuotes[4*i+1];

            int diff = Nbeg + LWG - (Ncur + sz);
            if( diff >= 0 ) {
                flags[Ncur] = 1;  flags[N+Ncur] = sz;  flags[N+N+Ncur] = i;
                Ncur += sz;
                if(diff == 0) { Nbeg += LWG; }
            } else {
                Nbeg += LWG;
                int ids = Nbeg - Ncur;
                flags[Ncur] = -1; flags[N+Ncur] = -ids; flags[N+N+Ncur] = i-1; //-1; //i; 
                flags[Nbeg] = 1;  flags[N+Nbeg] = sz;  flags[N+N+Nbeg] = i;
                Ncur = Nbeg + sz;
            }
        }

        // mark iddle threads for the last workgroup
        if( (Ncur % LWG) != 0 ) {
            flags[Ncur] = -1; flags[N+Ncur] = -(Nbeg+LWG-Ncur); flags[N+N+Ncur] = NUM_SWAP_QUOTES-1;//-1;
        }

        { // segmented scan with (+) to compute the indexes in
          // iota(size) and the corresponding swap indexes.
            int sz = 0, ind = 0, count = 0;
            for( int i = 0; i < N; i++ ) {
                if(flags[i] != 0) {
                    sz    = flags[  N+i]; 
                    ind   = flags[N+N+i]; 
                    count = 0;
                } else {
                    //flags[  N+i] = sz;
                    flags[N+N+i] = ind; // swaption global index
                    count++;
                }
                flags[  N+i] = count; // iota
                flags[3*N+i] = sz;    // size
            }
        }

#if 0
        printf("Flag Array Is: \n { \n ");
        for( int i = N; i < 2*N; i += LWG ) {
            printf("Row %d: { ", i / LWG);
            for( int j = 0; j < LWG; j++ ) {
                printf("%d, ", flags[i + j]);
            }
            printf(" }\n");
        }
        printf(" \n } \n\n\n");
#endif

    }

    return flags;
}

/**
 * Input  array of size N,
 * Result array of size M.
 */
int* getStartInd(const uint& N, short* flags, const uint& M) {
    int   count      = 0;
    int  *start_inds = new int[M];
    
    // compute the start index into expanded arrays 
    // for all NUM_SWAP_QUOTES, i.e., M, iterations
    for( UINT k = 0; k < N; k++ ) {
        if( flags[k] == 1 ) {
            start_inds[count++] = k;
        }
    }
    assert(count == M && "Invariant failed in getStartInd!");

    return start_inds;
}


#endif // ifndef MAKE_IRREG_SHAPE

