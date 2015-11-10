#ifndef NORDEA_CPU
#define NORDEA_CPU

#include "../includeC/ParPrefixUtil.h"

////////////////////////////////////////////////////////////
//// FINALLY THE MEAT: The CPU Version of the Main Loop!
////////////////////////////////////////////////////////////

void iteration_expanded_CPU (
        const unsigned int time_ind,
        const REAL         alpha,
        const REAL         beta,
        const REAL         nu,
        REAL*              a, // [OUTER_LOOP_COUNT*NUM_XY]
        REAL*              b, // [OUTER_LOOP_COUNT*NUM_XY]
        REAL*              c, // [OUTER_LOOP_COUNT*NUM_XY]
        REAL*              y, // [OUTER_LOOP_COUNT*NUM_XY]
        REAL*              u, // [OUTER_LOOP_COUNT*NUM_XY]
        REAL*              v, // [OUTER_LOOP_COUNT*NUM_XY]
        REAL*        scan_tmp // [4*OUTER_LOOP_COUNT*NUM_XY]
) {
    unsigned int i, j, k;

    REAL dtInv = 1.0/(myTimeline[time_ind+1]-myTimeline[time_ind]);

    REAL *res_arr = &myResArr[0];

#pragma omp parallel for default(shared) private(k,j,i) if(OUTER_LOOP_COUNT>4)
    for( k=0; k<OUTER_LOOP_COUNT; ++k ) {

        const unsigned int glob_ind = k*NUM_X*NUM_Y;

        //  explicit x and y
        for(j=0; j<NUM_Y; j++)  {
            for(i=0; i<NUM_X; i++) {

                // decls (constants)
                const unsigned int ind = glob_ind + j*NUM_X + i;
                unsigned int mul3;
                REAL tmp; //tmp2, tmp1 = dtInv*res_arr[j*NUM_X+i]; //

                REAL cur_myMuX  = 0.0, cur_myMuY = 0.0;
                REAL cur_myVarY = nu*nu;
                REAL cur_myVarX = exp(2*(beta*log(myX[i]) + myY[j] - 0.5*nu*nu*myTimeline[time_ind]));

                // second loop
                mul3 = REAL3_CT * j;
                tmp = 0.0;

                if(j!=0)
                    tmp += (cur_myMuY*myDy[mul3] + 0.5*cur_myVarY*myDyy[mul3])*res_arr[ind-NUM_X]; //[(j-1)*NUM_X+i];

                tmp += (cur_myMuY*myDy[mul3+1] + 0.5*cur_myVarY*myDyy[mul3+1])*res_arr[ind]; // [j*NUM_X+i];

                if(j!=NUM_Y-1)
                    tmp += (cur_myMuY*myDy[mul3+2] + 0.5*cur_myVarY*myDyy[mul3+2])*res_arr[ind+NUM_X]; //[(j+1)*NUM_X+i];

                v[glob_ind + i*NUM_Y + j] = tmp;   // v[i*NUM_Y + j] = tmp;     //v[i][j] = tmp2;

                // first loop
                mul3 = REAL3_CT * i;
                tmp += dtInv*res_arr[ind]; // [j*NUM_X+i];

                if(i!=0)
                    tmp += 0.5*(cur_myMuX*myDx[mul3] + 0.5*cur_myVarX*myDxx[mul3])*res_arr[ind-1]; //[j*NUM_X+i-1];

                tmp += 0.5*(cur_myMuX*myDx[mul3+1] + 0.5*cur_myVarX*myDxx[mul3+1])*res_arr[ind]; //[j*NUM_X+i];

                if(i!=NUM_X-1)
                    tmp += 0.5*(cur_myMuX*myDx[mul3+2] + 0.5*cur_myVarX*myDxx[mul3+2])*res_arr[ind+1]; //[j*NUM_X+i+1];

                u[ind] = tmp; //u[j*NUM_X + i] = tmp1 + tmp2; //u[j][i] = tmp1 + tmp2;

                // third loop (computing the implicit x parameters)
                a[ind] =        - 0.5*(cur_myMuX*myDx[mul3  ] + 0.5*cur_myVarX*myDxx[mul3  ]); // a[j*NUM_X+i] = ...
                b[ind] =  dtInv - 0.5*(cur_myMuX*myDx[mul3+1] + 0.5*cur_myVarX*myDxx[mul3+1]); // b[j*NUM_X+i] = ...
                c[ind] =        - 0.5*(cur_myMuX*myDx[mul3+2] + 0.5*cur_myVarX*myDxx[mul3+2]); // c[j*NUM_X+i] = ...
            }
        }
    }    // end OUTER_LOOP

#pragma omp parallel for default(shared) private(k,j,i) if(OUTER_LOOP_COUNT>4)
    for( k=0; k<OUTER_LOOP_COUNT; ++k ) {
        // implicit x tridag
        for(j=0;j<NUM_Y;j++) {

            const unsigned int ind = k*NUM_X*NUM_Y + j*NUM_X;

#if 1
            // As a (vectorized) scan (the vector is dimension "y")
            tridag_scan_array(
                    a+ind, b+ind, c+ind,
                    u+ind, NUM_X, u+ind,
                    y+ind, scan_tmp+4*ind
                );
#else
            tridag( a+ind, b+ind, c+ind, 
                    u+ind, NUM_X, u+ind, 
                    y+ind );
#endif
        }
    }    // end OUTER_LOOP

#pragma omp parallel for default(shared) private(k,j,i) if(OUTER_LOOP_COUNT>4)
    for( k=0; k<OUTER_LOOP_COUNT; ++k ) {
        //  implicit y
        for(i=0;i<NUM_X;i++) {
            for(j=0;j<NUM_Y;j++) {

                const unsigned int ind = k*NUM_X*NUM_Y + i*NUM_Y + j;

                unsigned int mul3j = REAL3_CT*j;
                double cur_myMuY = 0.0;
                double cur_myVarY = nu*nu;

                a[ind] =        - 0.5*(cur_myMuY*myDy[mul3j  ] + 0.5*cur_myVarY*myDyy[mul3j  ]); // a[i*NUM_Y+j] = ...
                b[ind] =  dtInv - 0.5*(cur_myMuY*myDy[mul3j+1] + 0.5*cur_myVarY*myDyy[mul3j+1]); // b[i*NUM_Y+j] = ...
                c[ind] =        - 0.5*(cur_myMuY*myDy[mul3j+2] + 0.5*cur_myVarY*myDyy[mul3j+2]); // c[i*NUM_Y+j] = ...

                v[ind] =  dtInv * u[k*NUM_X*NUM_Y + j*NUM_X+i] - 0.5*v[ind];  // y[j] = dtInv*u[j][i] - 0.5*v[i][j];
            }
        }
    }    // end OUTER_LOOP

#pragma omp parallel for default(shared) private(k,j,i) if(OUTER_LOOP_COUNT>4)
    for( k=0; k<OUTER_LOOP_COUNT; ++k ) {
        for(i=0; i<NUM_X; i++) {

            const unsigned int ind = k*NUM_X*NUM_Y + i*NUM_Y;

#if 1
            // As a (vectorized) scan (the vector is dimension "x")
            tridag_scan_array(
                    a+ind, b+ind, c+ind,
                    v+ind, NUM_Y, y+ind,
                    u+ind, scan_tmp+2*ind
                );
#else
            tridag( a+ind, b+ind, c+ind, 
                    v+ind, NUM_Y, y+ind, 
                    u+ind );
#endif
        }
    }    // end OUTER_LOOP

#pragma omp parallel for default(shared) private(k,j,i) if(OUTER_LOOP_COUNT>4)
    for( k=0; k<OUTER_LOOP_COUNT; ++k ) {  // segmented transpose followed by copy out!
        const unsigned int glb_ind = k*NUM_X*NUM_Y;

        for(i=0; i<NUM_X; i++) {
            for(int j=0; j<NUM_Y; j++) {
               res_arr[glb_ind + j*NUM_X+i] = y[glb_ind + i*NUM_Y+j];
            }
        }
    }    // end OUTER_LOOP
}

#endif // end include NORDEA_CPU
