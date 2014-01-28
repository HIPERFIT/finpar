#ifndef NORDEA_DATA_STRUCT
#define NORDEA_DATA_STRUCT

#include <vector>
#include <cmath>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include <cmath>

using namespace std;

///////////////////////////////////////////////////////
//// FLAGS and CONSTANT SCALARS!
///////////////////////////////////////////////////////

#define DEBUG           1

//#define DEBUG_PRINT_GPU_INFO

//typedef double       REAL;
//typedef unsigned int UINT;

const unsigned int OUTER_LOOP_COUNT = 128; //128; //1024; //100;
const unsigned int NUM_X            = 128; //256; //256; //64; //256;
const unsigned int NUM_Y            = 128; //32; //32; //32;
const unsigned int NUM_T            = 64; //64; //64; //64;
const unsigned int NUM_XY           = NUM_X*NUM_Y;

///////////////////////////////////////////////////////
//// GLOBAL ARRAYS !
///////////////////////////////////////////////////////

///  grid  ///
REAL myX       [NUM_X]; // 1-dim, size: NUM_X
REAL myY       [NUM_Y]; // 1-dim, size: NUM_Y
REAL myTimeline[NUM_T]; // 1-dim, size: NUM_T

unsigned int  myXindex, myYindex;

///  variable  ///
REAL myResArr[OUTER_LOOP_COUNT * NUM_Y * NUM_X]; // 3-dim, size: OUTER_LOOP_COUNT x NUM_Y x NUM_X


///  coeffs  ///
REAL myMuX [NUM_X * NUM_Y], myVarX[NUM_X * NUM_Y]; // 2-dim, size: NUM_X x NUM_Y
REAL myMuY [NUM_X * NUM_Y], myVarY[NUM_X * NUM_Y]; // 2-dim, size: NUM_X x NUM_Y

//  operators
REAL myDx[NUM_X * REAL3_CT], myDxx[NUM_X * REAL3_CT]; // 2-dim, size: NUM_X x 3
REAL myDy[NUM_Y * REAL3_CT], myDyy[NUM_Y * REAL3_CT]; // 2-dim, size: NUM_Y x 3

#endif // end include NORDEA_DATA_STRUCT
