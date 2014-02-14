#ifndef PAR_PREFIX_UTIL
#define PAR_PREFIX_UTIL

#include <time.h>
#include <sys/timeb.h>

// TIMING

typedef struct timeb mlfi_timeb;
#define mlfi_ftime ftime

#define mlfi_diff_time(t1,t2) \
  (t1.time - t2.time) * 1000 + (t1.millitm - t2.millitm)

/////////////////////////////////////////////////////////////


/**
 * Multiplies 2x2 matrixes `a' and `b' and
 * stores the result in(-place in) `a'.
 */
void matmult2(REAL* a, REAL* b) {
    REAL a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3],
         val = 1.0/(a0*b[0]);
            //val = 1.0/(a0*b[0]+a1*b[2]);
    a[0] = (a0*b[0] + a1*b[2])*val;
    a[1] = (a0*b[1] + a1*b[3])*val;
    a[2] = (a2*b[0] + a3*b[2])*val;
    a[3] = (a2*b[1] + a3*b[3])*val;
}


void scan_matmult_2by2(REAL* tmp, int N) {
    REAL* prev = NULL;

    for(int i=1; i<N; i++) {
        prev = tmp;
        tmp += 4;

        matmult2(tmp, prev);

    }
}

REAL map_matmult(REAL* tmp, REAL val1, REAL val2) {
#if 1
    REAL nom   = tmp[0]*val1+tmp[1]*val2;
    REAL denom = tmp[2]*val1+tmp[3]*val2;

    return nom/denom;
#else
    REAL denom = tmp[2]*val1+tmp[3]*val2;
    REAL nom = (tmp[0]/denom)*val1 + (tmp[1]/denom)*val2;

    return nom;
#endif
}

// Linear Function Composition and Scan
void linfuncomp(REAL* ab1, REAL* ab2) {
    //ab2[0] = ab1[0]+ab1[1]*ab2[0];
    //ab2[1] = ab1[1]*ab2[1];

    ab2[0] = ab2[0]+ab2[1]*ab1[0];
    ab2[1] = ab1[1]*ab2[1];
}

void scan_linfuncomp(REAL* tmp, int N) {
    REAL* prev = NULL;

    for(int i=1; i<N; i++) {
        prev = tmp;
        tmp += 2;

        linfuncomp(prev, tmp);
    }
}

REAL map_linfuncomp(REAL* tmp, REAL val) {
    return tmp[0] + tmp[1]*val;
}


//  tridag
void tridag_scan_array(
    const REAL*   a,
    const REAL*   b,
    const REAL*   c,
    const REAL*   d,
    int           n,
    REAL*         y,
    REAL*         u,
    REAL*         tmp)
{
    int    i, k;
    REAL beta;

    // initialize u[0]
    u[0] = b[0];
#if 0
    for(i=1; i<n; i++)
    {
        beta = a[i] / u[i-1];

        u[i] = b[i] - beta*c[i-1];
        //y[i] = d[i] - beta*y[i-1];
    }
#else
    { // 1. solve u[i] = b[i] - (a[i]/u[i-1])*c[i-1]

        // initialize tmp[4..4*n-1]

        tmp[0] = 1.0; tmp[1] = 0.0; tmp[2] = 0.0; tmp[3] = 1.0; // NEW

        for(i=1,k=4; i<n; i++, k+=4) {   // NEW: k = 4;
            tmp[k]   = b[i];
            //tmp[k+1] = 0.0 - a[i]*c[i-1];
            tmp[k+1] = 0.0 - a[i]*c[i-1];   // NEW
            tmp[k+2] = 1.0;
            tmp[k+3] = 0.0;
        }

        // scan with matrix multiplication on tmp
        //scan_matmult_2by2(tmp, n-1);
        scan_matmult_2by2(tmp, n);  // NEW

        //printf("tmp[0]: %f, %f, %f, %f\n", tmp[0], tmp[1], tmp[2], tmp[3]);
        //printf("b[1]: %f, a[1]: %f, c[0]: %f, a[1]*c[0]: %f\n", b[1], a[1], c[0], a[1]*c[0]);

        // map the result back:
        for(i=1,k=4; i<n; i++, k+=4) {   // NEW k=4
            u[i] = map_matmult(tmp+k, b[0], 1.0);

#if 0
            REAL beta = a[i] / u[i-1];

            beta = b[i] - beta*c[i-1];

            REAL err = beta - u[i];
            if(err < 0.0) err = 0.0 - err;
            if(err > 0.001) {
                printf("ERROR: i: %d, u[%d]= %.16f != %.16f\n\n", i, i, u[i], beta);

                printf("elem %d: %.8f, %.8f, %.8f, %.8f, b: %.8f, a: %.8f, c: %.8f, b0: %.8f, u: %.8f\n",
                        i, tmp[k], tmp[k+1], tmp[k+2], tmp[k+3], b[i], a[i], c[i-1], b[0], u[i]);
                exit(0);
            }
#endif
        }

    }
#endif

    // init y[0]
    y[0] = d[0];
#if 0
    for(i=1; i<n; i++)
    {
        beta = a[i] / u[i-1];
        y[i] = d[i] - beta*y[i-1];
    }
#else
    { // 2. compute the forward recurrence y[i] = d[i] - (a[i]/u[i-1])*y[i-1]
        //init tmp[2..2*n-1]
        tmp[0] = 0.0; tmp[1] = 1.0; // NEW
        for(i=1,k=2; i<n; i++, k+=2) {  // NEW: k=2
            tmp[k]   = d[i];
            tmp[k+1] = 0.0 - a[i]/u[i-1];
        }

        // forward scan with linear function composition on tmp
        scan_linfuncomp(tmp, n); // NEW: n

        // map the result back:
        for(i=1,k=2; i<n; i++, k+=2) {  // NEW: k = 2;
            y[i] = map_linfuncomp(tmp+k, y[0]);
        }
    }
#endif

    y[n-1] = y[n-1]/u[n-1];
#if 0
    for(i=n-2; i>=0; i--) {
        y[i] = (y[i] - c[i]*y[i+1]) / u[i];
    }
#else
    { // 2. compute the backward recurrence y[i] = (y[i] - c[i]*y[i+1]) / u[i];
        //init tmp[2..2*n-1]
        tmp[0] = 0.0; tmp[1] = 1.0; // NEW
        for(i=n-2,k=2; i>=0; i--, k+=2) {  // NEW k=2
            tmp[k]   = y[i]/u[i];
            tmp[k+1] = 0.0 - c[i]/u[i];
        }

        // backward scan with linear function composition on tmp
        scan_linfuncomp(tmp, n); // NEW: n

        // map the result back:
        for(i=n-2,k=2; i>=0; i--, k+=2) {   // NEW k=2
            y[i] = map_linfuncomp(tmp+k, y[n-1]);
        }
    }
#endif
}

void test_matmult(){
    int N = 16;
    REAL* arr = new REAL[N*4];
    REAL* tmp = arr;
    for(int i=0; i<N; i++){
        tmp[0] = 1.0+i; tmp[1] = 0.0;
        tmp[2] = 1.0+i; tmp[3] = 0.0;
        tmp += 4;
    }

    scan_matmult_2by2(arr, N);

    tmp = arr;
    for(int i=0; i<N; i++){
        printf("(%f, %f, %f, %f), ", tmp[0], tmp[1], tmp[2], tmp[3]);
        tmp += 4;
    }
    printf("\n\n");

    delete[] arr;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void tridag_scan_array_womatmult(
    const REAL*   a,
    const REAL*   b,
    const REAL*   c,
    const REAL*   d,
    int           n,
    REAL*         y,
    REAL*         u,
    REAL*         tmp)
{
    int    i, k;
    REAL beta;

#if 0
    // initialize u[0]
    u[0] = b[0];
    { // 1. solve u[i] = b[i] - (a[i]/u[i-1])*c[i-1]

        // initialize tmp[4..4*n-1]

        for(i=1,k=0; i<n; i++, k+=4) {
            tmp[k]   = b[i];
            tmp[k+1] = 0.0 - a[i]*c[i-1];
            tmp[k+2] = 1.0;
            tmp[k+3] = 0.0;
        }

        // scan with matrix multiplication on tmp
        scan_matmult_2by2(tmp, n-1);

        //printf("tmp[0]: %f, %f, %f, %f\n", tmp[0], tmp[1], tmp[2], tmp[3]);
        //printf("b[1]: %f, a[1]: %f, c[0]: %f, a[1]*c[0]: %f\n", b[1], a[1], c[0], a[1]*c[0]);

        // map the result back:
        for(i=1,k=4; i<n; i++, k+=4) {
            u[i] = map_matmult(tmp+k, b[0], 1.0);

        }

    }
#else
    // scan with matrix multiplication on tmp
    scan_matmult_2by2(tmp, n-1);

    //printf("tmp[0]: %f, %f, %f, %f\n", tmp[0], tmp[1], tmp[2], tmp[3]);
    //printf("b[1]: %f, a[1]: %f, c[0]: %f, a[1]*c[0]: %f\n", b[1], a[1], c[0], a[1]*c[0]);

    // map the result back:
    for(i=1,k=4; i<n; i++, k+=4) {
        u[i] = map_matmult(tmp+k, b[0], 1.0);

    }
#endif

#if 0
    // init y[0]
    y[0] = d[0];
    { // 2. compute the forward recurrence y[i] = d[i] - (a[i]/u[i-1])*y[i-1]
        tmp[0] = 0.0; tmp[1] = 1.0;
        //init tmp[2..2*n-1]
        for(i=1,k=2; i<n; i++, k+=2) {
            tmp[k]   = d[i];
            tmp[k+1] = 0.0 - a[i]/u[i-1];
        }

        // forward scan with linear function composition on tmp
        scan_linfuncomp(tmp, n); // n-1

        // map the result back:
        for(i=1,k=2; i<n; i++, k+=2) {
            y[i] = map_linfuncomp(tmp+k, y[0]);
        }
    }
#else

    // forward scan with linear function composition on tmp
    //scan_linfuncomp(tmp, n); // n-1

    // map the result back:
    //for(i=1,k=2; i<n; i++, k+=2) {
    //    y[i] = map_linfuncomp(tmp+k, y[0]);
    //}
#endif

#if 1
    y[n-1] = y[n-1]/u[n-1];
    { // 2. compute the forward recurrence y[i] = (y[i] - c[i]*y[i+1]) / u[i];
        //init tmp[2..2*n-1]
        for(i=n-2,k=0; i>=0; i--, k+=2) {
            tmp[k]   = y[i]/u[i];
            tmp[k+1] = 0.0 - c[i]/u[i];
        }

        // forward scan with linear function composition on tmp
        scan_linfuncomp(tmp, n-1);

        // map the result back:
        for(i=n-2,k=0; i>=0; i--, k+=2) {
            y[i] = map_linfuncomp(tmp+k, y[n-1]);
        }
    }
#else
    {
        //scan_linfuncomp(tmp, n-1);

        // map the result back:
        for(i=n-2,k=0; i>=0; i--, k+=2) {
            y[i] = map_linfuncomp(tmp+k, y[n-1]);
        }
    }
#endif


}

inline void tridag_seq_32(
    REAL*   a,
    REAL*   b,
    REAL*   c,
    REAL*   d,
    int     n,
    REAL*   y,
    REAL*   u)
{
    int    i, UB, offset;
    REAL   beta;

    y[0] = d[0];
    u[0] = b[0];

    UB = n << 5;
    for(i=32; i<UB; i+=32) {
        beta = a[i] / u[i-32];

        u[i] = b[i] - beta*c[i-32];
        y[i] = d[i] - beta*y[i-32];
    }

    UB -= 32;
    y[UB] = y[UB] / u[UB];

    UB = (n-2) * 32;
    for(i=UB; i>=0; i-=32) {
        y[i] = (y[i] - c[i]*y[i+32]) / u[i];
    }
}

inline void tridag_seq_1(
    REAL*   a,
    REAL*   b,
    REAL*   c,
    REAL*   d,
    int     n,
    REAL*   y,
    REAL*   u)
{
    int    i, UB, offset;
    REAL   beta;

    y[0] = d[0];
    u[0] = b[0];

    UB = n << 5;
    for(i=1; i<n; i++) {
        beta = a[i] / u[i-1];

        u[i] = b[i] - beta*c[i-1];
        y[i] = d[i] - beta*y[i-1];
    }

    y[n-1] = y[n-1] / u[n-1];

    for(i=n-2; i>=0; i--) {
        y[i] = (y[i] - c[i]*y[i+1]) / u[i];
    }
}



#endif // end include PAR_PREFIX_UTIL

