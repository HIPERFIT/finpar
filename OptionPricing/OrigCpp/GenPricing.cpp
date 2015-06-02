#define FAST_BB

#include "ParseInput.h"
#include "SobolGaussBB.h"
#include "Contracts.h"

double* run_CPUkernel(  const LoopROScalars   &   scals,
                        const SobolArrays     &   sob_arrs,
                        const ModelArrays     &   md_arrs,
                        const BrowBridgeArrays&   bb_arrs
) {
    const UINT sob_dim      = scals.num_under * scals.num_dates;
    const UINT num_under_sq = scals.num_under * scals.num_under;
    const int* sob_matT     = sob_arrs.sobol_dirvcts;

    // allocate for Ps threads
    const UINT Sd = do_padding( sob_dim          );
    UINT* sob_vct = static_cast<UINT*>( malloc( Sd * sizeof(UINT)) );
    REAL*  md_vct = static_cast<REAL*>( malloc( Sd * sizeof(REAL)) );
    REAL* trj_vct = static_cast<REAL*>( malloc( Sd * sizeof(REAL)) );

    const UINT Sv = do_padding( scals.num_models );
    double* vhat  = static_cast<double*>( malloc( Sv * sizeof(double)) );

    // init the prices to 0
    for(int i = 0; i < scals.num_models; i ++) {
        vhat[ i ] = 0.0;
    }

    // Sobol number 0 is a vector of 0s. 
    for(int i = 0; i < sob_dim; i++) {
        sob_vct[i] = 0;
    }

    for( int i = 0; i < scals.num_mcits; i ++ ) {

        // Sobol strength-reduced (recurrence) formula.
        sobolRec( sob_dim, scals.sobol_bits, sob_matT, i+1, sob_vct );

        // transform the normal [0,1) to gaussian distribution [-inf, +inf]
        uGaussian( scals.sob_norm_fact, sob_dim, sob_vct, md_vct );

        // correlate the dates on each path using a Brownian Bridge
        brownianBridge( scals.num_under,    scals.num_dates, 
                        bb_arrs.bb_inds,    bb_arrs.bb_data, 
                        md_vct,             trj_vct         );

#ifdef FAST_BB
        REAL* traj = md_vct;
#else
        REAL* traj = trj_vct; //md_vct;
#endif

        // compute trajectory
        for ( int m = 0; m < scals.num_models; m ++ ) {
            REAL* md_c      = md_arrs.md_c      + ( m * num_under_sq    );
            REAL* md_vols   = md_arrs.md_vols   + ( m * sob_dim         );
            REAL* md_drifts = md_arrs.md_drifts + ( m * sob_dim         );
            REAL* md_starts = md_arrs.md_starts + ( m * scals.num_under );
    
            for( int i = 0; i < scals.num_dates; i ++ )  {
                for( int j = 0; j < scals.num_under; j ++ ) {
                    REAL temp = 0.0;
                    int k = i*scals.num_under + j;

                    for ( int l = 0; l <= j; l ++ ) {
#ifdef FAST_BB
                        REAL md_val = trj_vct[i*scals.num_under + l];
#else
                        REAL md_val = md_vct [i*scals.num_under + l];
#endif
                        temp += md_c[j*scals.num_under + l] * md_val;
                    }

                    temp = exp( temp*md_vols[k] + md_drifts[k] );

                    traj[k] = (k < scals.num_under) ? 
                                    md_starts[ k ] * temp            :
                                    traj[k - scals.num_under] * temp ;
                }
            }

            aggregDiscountedPayoff( m,  scals.contract,  
                                        scals.num_under,      
                                        scals.num_cash_flows, 
                                        scals.num_det_pricers,   
                                        md_arrs.md_discts,  
                                        md_arrs.md_detvals,   
                                        traj,   vhat      );                
        }        
    } 

    
    for(int i = 0; i < scals.num_models; i ++) {
        vhat[i] = vhat[i] / scals.num_mcits;
    }

    // clean-up!
    free(sob_vct);
    free( md_vct);
    free(trj_vct);

    return vhat;
}


int main() {
    LoopROScalars    scals;
    SobolArrays      sob_arrs;
    ModelArrays      md_arrs;
    BrowBridgeArrays bb_arrs;

    fprintf(stdout, "\n// Generic Pricing, Original Sequential Benchmark:\n");

    readDataSet(scals, sob_arrs, md_arrs, bb_arrs);

    fprintf(stdout, "// Contract: %d, MC Its#: %d, #Underlyings: %d, #Path Dates: %d, chunk: %d\n\n", 
            scals.contract, scals.num_mcits, scals.num_under, scals.num_dates, scals.chunk      );

    double* prices;
    unsigned long int elapsed;
    { // run kernel
        struct timeval t_start, t_end, t_diff;
        gettimeofday(&t_start, NULL);

        computeSobolFixIndex( sob_arrs, scals.chunk );
        { // do work and cleanup
            prices = run_CPUkernel( scals, sob_arrs, md_arrs, bb_arrs );
        }

        gettimeofday(&t_end, NULL);
        timeval_subtract(&t_diff, &t_end, &t_start);
        elapsed = t_diff.tv_sec*1e6+t_diff.tv_usec;

        md_arrs .cleanup();
        bb_arrs .cleanup();
        sob_arrs.cleanup();
    }

    {   // validation and writeback of the result
        bool is_valid = validate   ( scals.num_models, prices );
        writeStatsAndResult( is_valid, scals.num_models, prices, false, 1, elapsed ); 
        free(prices);       
    }

    return 0;
}

