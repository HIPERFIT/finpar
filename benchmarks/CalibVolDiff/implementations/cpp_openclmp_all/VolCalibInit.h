#ifndef NORDEA_INIT
#define NORDEA_INIT

#include "DataStructConst.h"

///////////////////////////////////////////////////////
//// INITIALISATION !
///////////////////////////////////////////////////////

void initGrid(
        const real_t s0,
        const real_t alpha,
        const real_t nu,
        const real_t t,
        const unsigned int numX,
        const unsigned int numY,
        const unsigned int numT
) {
    unsigned int i;

    for( i=0; i<numT; ++i )
        myTimeline[i] = t*i/(numT-1);

    const real_t stdX = 20*alpha*s0*sqrt(t);
    const real_t dx   = stdX/numX;
    myXindex = static_cast<unsigned int>(s0/dx);

    for( i=0; i<numX; ++i )
        myX[i] = i*dx - myXindex*dx + s0;

    const real_t stdY     = 10*nu*sqrt(t);
    const real_t dy       = stdY/numY;
    const real_t logAlpha = log(alpha);

    myYindex = numY/2;

    for( i=0; i<numY; ++i )
        myY[i] = i*dy - myYindex*dy + logAlpha;
}

void initOperator(
        const unsigned int N,
        const real_t*        x,   // size: N
        real_t*              Dx,  // size: N x 3
        real_t*              Dxx  // size: N x 3
) {
    real_t dxl, dxu;
    unsigned int i, k;

    //  lower boundary
    dxl      =  0.0;
    dxu      =  x[1] - x[0];

    Dx[0]  =  0.0;     // Dx[0][0]
    Dx[1]  = -1.0/dxu; // Dx[0][1]
    Dx[2]  =  1.0/dxu; // Dx[0][2]

    Dxx[0] =  0.0;  // Dxx[0][0]
    Dxx[1] =  0.0;  // Dxx[0][1]
    Dxx[2] =  0.0;  // Dxx[0][2]

#if (REAL3_CT == 4)
    Dx[3] = 1.0; Dxx[3] = 1.0;
#endif

    //  standard case
    for(i=1, k=REAL3_CT; i<N-1; i++, k+=REAL3_CT) {
        dxl      = x[i]   - x[i-1];
        dxu      = x[i+1] - x[i];

        Dx[k  ]  = -dxu/dxl/(dxl+dxu);                 // Dx[i][0]
        Dx[k+1]  = (dxu/dxl - dxl/dxu)/(dxl+dxu);      // Dx[i][1]
        Dx[k+2]  =  dxl/dxu/(dxl+dxu);                 // Dx[i][2]

        Dxx[k  ] =  2.0/dxl/(dxl+dxu);                 // Dxx[i][0]
        Dxx[k+1] = -2.0*(1.0/dxl + 1.0/dxu)/(dxl+dxu); // Dxx[i][1]
        Dxx[k+2] =  2.0/dxu/(dxl+dxu);                 // Dxx[i][2]

#if (REAL3_CT == 4)
        Dx[k+3] = 1.0; Dxx[k+3] = 1.0;
#endif
    }

    //  upper boundary
    dxl        =  x[N-1] - x[N-2];
    dxu        =  0.0;

    k = (N-1)*REAL3_CT;
    Dx[k  ]  = -1.0/dxl; // Dx[N-1][0]
    Dx[k+1]  =  1.0/dxl; // Dx[N-1][1]
    Dx[k+2]  =  0.0;     // Dx[N-1][2]

    Dxx[k  ] = 0.0;      // Dxx[N-1][0]
    Dxx[k+1] = 0.0;      // Dxx[N-1][1]
    Dxx[k+2] = 0.0;      // Dxx[N-1][2]

#if (REAL3_CT == 4)
    Dx[k+3] = 1.0; Dxx[k+3] = 1.0;
#endif
}

void whole_loop_nest_init (
        const real_t  s0,
        const real_t  t,
        const real_t  alpha,
        const real_t  nu,
        const real_t  beta
) {
    initGrid(s0, alpha, nu, t, NUM_X, NUM_Y, NUM_T);
    initOperator(NUM_X, myX, myDx, myDxx);
    initOperator(NUM_Y, myY, myDy, myDyy);
}

///////////////////////////////////////////////////////
//// PAYOFF UPDATE !
///////////////////////////////////////////////////////

void setPayoff_expanded(const real_t strike, real_t* myResArr) {
    for(unsigned int i=0; i<NUM_X; i++) {
        real_t payoff = max( myX[i]-strike, (real_t)0.0 );
        for(unsigned int j=0; j<NUM_Y; j++) {
            myResArr[j*NUM_X+i] = payoff;
        }
    }
}

#endif // end include Init functionality!
