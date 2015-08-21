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

#define TRIDAG_ALL_OPT_ON
#define MOST_OPTIMISED_ON

//#define IS_GPU
//#define DEBUG_PRINT_GPU_INFO

unsigned int OUTER_LOOP_COUNT; //128; //1024; //100;
unsigned int NUM_X           ; //256; //256; //64; //256;
unsigned int NUM_Y           ; //32; //32; //32;
unsigned int NUM_T           ; //64; //64; //64;
unsigned int NUM_XY          ; //NUM_X*NUM_Y;

///////////////////////////////////////////////////////
//// GLOBAL ARRAYS !
///////////////////////////////////////////////////////

///  grid  ///
real_t* myX       ; // 1-dim, size: NUM_X
real_t* myY       ; // 1-dim, size: NUM_Y
real_t* myTimeline; // 1-dim, size: NUM_T

unsigned int myXindex, myYindex;

///  variable  ///
real_t* myResArr; //[OUTER_LOOP_COUNT * NUM_Y * NUM_X]; // 3-dim, size: OUTER_LOOP_COUNT x NUM_Y x NUM_X


///  coeffs  ///
real_t *myMuX, *myVarX; // 2-dim, size: NUM_X x NUM_Y
real_t *myMuY, *myVarY; // 2-dim, size: NUM_X x NUM_Y

//  operators
real_t *myDx, *myDxx; // 2-dim, size: NUM_X x REAL3_CT
real_t *myDy, *myDyy; // 2-dim, size: NUM_Y x REAL3_CT

void allocGlobArrs() {
    myX = new real_t[NUM_X];
    myY = new real_t[NUM_Y];
    myTimeline = new real_t[NUM_T];

    myMuX = new real_t[NUM_X * NUM_Y];
    myVarX= new real_t[NUM_X * NUM_Y];
    myMuY = new real_t[NUM_X * NUM_Y];
    myVarY= new real_t[NUM_X * NUM_Y];

    myDx  = new real_t[NUM_X * REAL3_CT];
    myDxx = new real_t[NUM_X * REAL3_CT];
    myDy  = new real_t[NUM_Y * REAL3_CT];
    myDyy = new real_t[NUM_Y * REAL3_CT];

    myResArr = new real_t[OUTER_LOOP_COUNT * NUM_Y * NUM_X];
}

void deallocGlobArrs() {
    delete[] myX;
    delete[] myY;
    delete[] myTimeline;

    delete[] myMuX;
    delete[] myVarX;
    delete[] myMuY;
    delete[] myVarY;

    delete[] myDx;
    delete[] myDxx;
    delete[] myDy;
    delete[] myDyy;

    delete[] myResArr;
}

#endif // end include NORDEA_DATA_STRUCT
