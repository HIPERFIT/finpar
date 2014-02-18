#ifndef NORDEA_INIT
#define NORDEA_INIT

#include "../includeC/DataStructConst.h"

///////////////////////////////////////////////////////
//// INITIALISATION !
///////////////////////////////////////////////////////

void initGrid(
        const REAL s0,
        const REAL alpha,
        const REAL nu,
        const REAL t,
        const unsigned int numX,
        const unsigned int numY,
        const unsigned int numT
) {
    unsigned int i;

    for( i=0; i<numT; ++i )
        myTimeline[i] = t*i/(numT-1);

    const REAL stdX = 20*alpha*s0*sqrt(t);
    const REAL dx   = stdX/numX;
    myXindex = static_cast<unsigned int>(s0/dx);

    for( i=0; i<numX; ++i )
        myX[i] = i*dx - myXindex*dx + s0;

    const REAL stdY     = 10*nu*sqrt(t);
    const REAL dy       = stdY/numY;
    const REAL logAlpha = log(alpha);

    myYindex = numY/2;

    for( i=0; i<numY; ++i )
        myY[i] = i*dy - myYindex*dy + logAlpha;
}

void initOperator(
        const unsigned int N,
        const REAL*        x,   // size: N
        REAL*              Dx,  // size: N x 3
        REAL*              Dxx  // size: N x 3
) {
    REAL dxl, dxu;
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
        const REAL  s0,
        const REAL  t,
        const REAL  alpha,
        const REAL  nu,
        const REAL  beta
) {
    initGrid(s0, alpha, nu, t, NUM_X, NUM_Y, NUM_T);
    initOperator(NUM_X, myX, myDx, myDxx);
    initOperator(NUM_Y, myY, myDy, myDyy);
}

///////////////////////////////////////////////////////
//// PAYOFF UPDATE !
///////////////////////////////////////////////////////

void setPayoff_expanded(const REAL strike, REAL* myResArr) {
    for(unsigned int i=0; i<NUM_X; i++) {
        REAL payoff = max( myX[i]-strike, (REAL)0.0 );
        for(unsigned int j=0; j<NUM_Y; j++) {
            myResArr[j*NUM_X+i] = payoff;
        }
    }
}

#endif // end include Init functionality!
