#ifndef PARSE_INPUT
#define PARSE_INPUT

#include "ParserC.h"
#include "Util.h"
#include <assert.h>
#include <math.h>
#include <algorithm>

using namespace std;

#include <iostream>

#if REAL_IS_FLOAT
    #define read_real read_float
#else
    #define read_real read_double
#endif

const float EPS = 0.00001;

/***********************************/
/********** READ DATA SET **********/
/***********************************/

void readDataSet(   unsigned int& outer, 
                    unsigned int& num_X,
                    unsigned int& num_Y,
                    unsigned int& num_T
) {
    if(     read_int( static_cast<unsigned int*>( &outer ) ) ||
            read_int( static_cast<unsigned int*>( &num_X ) ) ||
            read_int( static_cast<unsigned int*>( &num_Y ) ) ||
            read_int( static_cast<unsigned int*>( &num_T ) )  ) {
        
        fprintf(stderr, "Syntax error when reading the dataset, i.e., four ints.\n");
        exit(1);
    }

    {   // check dataset invariants:
        bool atr_ok = true;

        atr_ok  = outer > 0;
        assert(atr_ok && "Outer loop count less than 0!");

        atr_ok  = (num_X > 0) && (num_X <= WORKGROUP_SIZE) && is_pow2(num_X); 
        assert(atr_ok && "Illegal NUM_X value!");

        atr_ok  = (num_Y > 0) && (num_Y <= WORKGROUP_SIZE) && is_pow2(num_Y); 
        assert(atr_ok && "Illegal NUM_X value!");

        atr_ok  = num_T > 0;
        assert(atr_ok && "NUM_T value less or equal to zero!!");
    }
}

#endif // PARSE_INPUT

