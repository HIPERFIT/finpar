#ifndef PARSE_INPUT
#define PARSE_INPUT

#include "ParserC.h"
#include "Constants.h"
#include <assert.h>
#include <math.h>
#include <algorithm>


#if _OPTIMIZATION_USE_FLOATS
    #define read_real read_float
#else
    #define read_real read_double
#endif


void readDataSet(   LoopROScalars& scals, int*& sobol_dirvcts, 
                    ModelArrays& md_arrs, BrowBridgeArrays& bb_arrs 
) {
    if( read_int( &scals.contract   ) ||
        read_int( &scals.num_mcits  ) ||
        read_int( &scals.num_dates  ) ||
        read_int( &scals.num_under  ) ||
        read_int( &scals.num_models ) ||
        read_int( &scals.sobol_bits )   ) {
        
        fprintf(stderr, "Syntax error when reading the first four ints params.\n");
        exit(1);
    }

    int64_t shape[3];

    { // reading sobol's direction vectors
        if (read_array(sizeof(int), read_int, (void**)&sobol_dirvcts, shape, 2) ) {
            fprintf(stderr, "Syntax error when reading sobol direction-vector matrix of dim [%d,%d].\n",
                            scals.num_under*scals.num_dates, scals.sobol_bits);
            exit(1);
        }

        bool ok = ( shape[1] == scals.sobol_bits ) && 
                  ( shape[0] == scals.num_under*scals.num_dates );
        assert(ok && "Incorrect shape of sobol direction vectors!");
    }

    { // reading the market (models) data
        bool ok;

        // md_c
        if (read_array(sizeof(REAL), read_real, (void**)&md_arrs.md_c, shape, 3) ) {
            fprintf(stderr, "Syntax error when reading md_c [%d,%d,%d].\n",
                            scals.num_models, scals.num_under, scals.num_under);
            exit(1);
        }
        ok = ( shape[0] == scals.num_models) && 
             ( shape[1] == scals.num_under ) && 
             ( shape[2] == scals.num_under )  ;
        assert(ok && "Incorrect shape of md_c!");

        // md_vols (volatility)
        if (read_array(sizeof(REAL), read_real, (void**)&md_arrs.md_vols, shape, 3) ) {
            fprintf(stderr, "Syntax error when reading md_vols [%d,%d,%d].\n",
                            scals.num_models, scals.num_dates, scals.num_under);
            exit(1);
        }
        ok = ( shape[0] == scals.num_models) && 
             ( shape[1] == scals.num_dates ) && 
             ( shape[2] == scals.num_under )  ;
        assert(ok && "Incorrect shape of md_vols!");

        // md_drifts
        if (read_array(sizeof(REAL), read_real, (void**)&md_arrs.md_drifts, shape, 3) ) {
            fprintf(stderr, "Syntax error when reading md_drifts [%d,%d,%d].\n",
                            scals.num_models, scals.num_dates, scals.num_under);
            exit(1);
        }
        ok = ( shape[0] == scals.num_models) && 
             ( shape[1] == scals.num_dates ) && 
             ( shape[2] == scals.num_under )  ;
        assert(ok && "Incorrect shape of md_drifts!");

        // md_starts
        if (read_array(sizeof(REAL), read_real, (void**)&md_arrs.md_starts, shape, 2) ) {
            fprintf(stderr, "Syntax error when reading md_drifts [%d,%d].\n",
                            scals.num_models, scals.num_under);
            exit(1);
        }
        ok = ( shape[0] == scals.num_models) && 
             ( shape[1] == scals.num_under )  ;
        assert(ok && "Incorrect shape of md_starts!");

        // md_detvals (Model deterministic values)
        if (read_array(sizeof(REAL), read_real, (void**)&md_arrs.md_detvals, shape, 2) ) {
            fprintf(stderr, "Syntax error when reading md_detvals [%d,..].\n",
                            scals.num_models);
            exit(1);
        }
        ok = ( shape[0] == scals.num_models );
        assert(ok && "Incorrect shape of md_detvals!");
        scals.num_det_pricers = shape[1];

        // md_discts (Model discounts)
        if (read_array(sizeof(REAL), read_real, (void**)&md_arrs.md_discts, shape, 2) ) {
            fprintf(stderr, "Syntax error when reading md_discts [%d,..].\n",
                            scals.num_models);
            exit(1);
        }
        ok = ( shape[0] == scals.num_models );
        assert(ok && "Incorrect shape of md_discts!");
        scals.num_cash_flows  = shape[1];
        printf("discount[2]: %f\n", md_arrs.md_discts[2]);
    }

    { // Brownian Bridge Indirect Arrays and Data
        if (read_array(sizeof(int), read_int, (void**)&bb_arrs.bb_inds, shape, 2) ) {
            fprintf(stderr, "Syntax error when reading Brownian Bridge indirect arrays [3,%d].\n",
                            scals.num_dates);
            exit(1);
        }
        bool ok = ( shape[0] == 3 ) && ( shape[1] == scals.num_dates );
        assert(ok && "Incorrect shape of bb_inds (brownian bridge indirect arrays)!");

                // md_starts
        if (read_array(sizeof(REAL), read_real, (void**)&bb_arrs.bb_data, shape, 2) ) {
            fprintf(stderr, "Syntax error when reading md_drifts [3,%d].\n",
                            scals.num_dates );
            exit(1);
        }
        ok = ( shape[0] == 3 ) && ( shape[1] == scals.num_dates )  ;
        assert(ok && "Incorrect shape of bb_data (brownian bridge data arrays)!");
    }
}

#endif // PARSE_INPUT
