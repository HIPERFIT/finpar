#define FAST_BB                     1
#define _OPTIMIZATION_USE_FLOATS    1


#include "ParseInput.h"
#include "TimeHelper.h"
#include "SobolGaussBB.h"
#include "Contracts.h"
#include <omp.h>

UINT do_padding(const UINT& n) {
    return (((n / 64) * 64) + 64);
}

void run_CPUkernel(    
                    LoopROScalars   &   scals,
                    int*                sob_mat,       
                    ModelArrays     &   md_arrs,
                    BrowBridgeArrays&   bb_arrs
) {
    const UINT sob_dim      = scals.num_under * scals.num_dates;
    const UINT num_under_sq = scals.num_under * scals.num_under;
    int Ps = -1;

    // get the number of OMP threads
#pragma omp parallel shared(Ps)
    {
            UINT th_id = omp_get_thread_num();
            if(th_id == 0) { Ps = omp_get_num_threads(); }
    }
    assert(Ps > 0 && "Number of OMP threads <= 0!");


    // allocate for Ps threads
    const UINT Sd = do_padding( sob_dim          );
    UINT* sob_glb_vct = static_cast<UINT*>( malloc( Ps * Sd * sizeof(UINT)) );
    REAL*  md_glb_vct = static_cast<REAL*>( malloc( Ps * Sd * sizeof(REAL)) );
    REAL* trj_glb_vct = static_cast<REAL*>( malloc( Ps * Sd * sizeof(REAL)) );

    const UINT Sv = do_padding( scals.num_models );
    double* vhat_glb  = static_cast<double*>( malloc( Ps * Sv * sizeof(double)) );

    // Transpose the sobol direction vectors!
    int* sob_matT = static_cast<int*> ( malloc( do_padding(sob_dim * scals.sobol_bits) * sizeof(int) ) );
    for( int j = 0; j < scals.sobol_bits; j++ ) {
        for( int i = 0; i< sob_dim; i++ ) {
            sob_matT[j*sob_dim + i] = sob_mat[i*scals.sobol_bits + j];
        }
    }

    printf("CHUNK: %d\n\n", scals.chunk);

    struct timeval t_start, t_end, t_diff;
    gettimeofday(&t_start, NULL);
    unsigned long int elapsed;

#pragma omp parallel
    {
        UINT  th_id   = omp_get_thread_num();
        UINT* sob_vct = sob_glb_vct + ( th_id * Sd );
        REAL*  md_vct =  md_glb_vct + ( th_id * Sd );
        REAL* trj_vct = trj_glb_vct + ( th_id * Sd );
        double* vhat  = vhat_glb    + ( th_id * Sv ); 

        for(int i = 0; i < scals.num_models; i ++) 
            vhat[ i ] = 0.0;

#pragma omp for schedule(dynamic) 
        for( int i = 0; i < scals.num_mcits; i += scals.chunk ) {

            unsigned bound = std::min(i + scals.chunk, scals.num_mcits);

            for( int k = i; k < bound; k ++ ) {
//            new_trajectoryGPU(k, k==i, ro_scal, ro_arr, priv_arr);

                // generate the next Sobol quasi-random vector: 1st iter uses
                // the independent formula; the others the recurrent formula.
                if ( k == i)    sobolInd( sob_dim, scals.sobol_bits, sob_matT, k, sob_vct );
                else            sobolRec( sob_dim, scals.sobol_bits, sob_matT, k, sob_vct );

                // transform the normal [0,1) to gaussian distribution [-inf, +inf]
                uGaussian( scals.sob_norm_fact, sob_dim, sob_vct, md_vct );

                // correlate the dates on each path using a Brownian Bridge
                brownianBridge( scals.num_under,    scals.num_dates, 
                                bb_arrs.bb_inds,    bb_arrs.bb_data, 
                                md_vct,             trj_vct         );
#if FAST_BB
                REAL* traj = md_vct;
#else
                REAL* traj = trj_vct; //md_vct;
#endif
#if 0
                if (k==1) {
                    printf("\n\nBB vector: [ ");
                    for ( int ii = 0; ii < scals.num_dates * scals.num_under; ii ++ ) {
                        printf("%8f, ",  md_vct[ii]);
                    }
                    printf(" ]\n\n");
                }
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
#if FAST_BB
                                REAL md_val = trj_vct[i*scals.num_under + l];
#else
                                REAL md_val = md_vct [i*scals.num_under + l];
#endif
                                temp += md_c[j*scals.num_under + l] * md_val;
                            }

                            temp = exp( temp*md_vols[k] + md_drifts[k] );

                            traj[k] = (k < scals.num_under) ? 
                                            md_starts[ k ] * temp               :
                                            traj[k - scals.num_under] * temp ;
                        }
                    }
#if 0
                if (k==1) {
                    printf("\n\nTrajectory vector: [ ");
                    for ( int ii = 0; ii < scals.num_dates * scals.num_under; ii ++ ) {
                        printf("%8f, ", traj[ii]);
                    }
                    printf(" ]\n\n");
                }
#endif
                    aggregDiscountedPayoff( m,  scals.contract,  
                                                scals.num_under,      
                                                scals.num_cash_flows, 
                                                scals.num_det_pricers,   
                                                md_arrs.md_discts,  
                                                md_arrs.md_detvals,   
                                                traj,   vhat      );

                }            
            }
        } 
    }

    for( int p = 1; p < Ps; p ++ ) {
        for(int i = 0; i < scals.num_models; i ++) { 
            vhat_glb[i] += vhat_glb[ i + p*Sv ];
        }
    }

    gettimeofday(&t_end, NULL);
    timeval_subtract(&t_diff, &t_end, &t_start);
    elapsed = t_diff.tv_sec*1e6+t_diff.tv_usec;

    printf("CPU Run Time: %lu !!! \n\n\n", elapsed);

    printf("Price for vhat[0] is: %.16f !!!\n\n", (vhat_glb[0] / scals.num_mcits) );

//  updateState(&priv_arr, seq_count, data);                          // UPDATE IMPLEM!
//  updateModels(models, &ro_scal, &priv_arr);

    md_arrs.cleanup();
    bb_arrs.cleanup();
    free(sob_glb_vct);
    free( md_glb_vct);
    free(trj_glb_vct);
    free( vhat_glb  );
}


int main() {
    int*             sobol_dirvcts;
    LoopROScalars    scals;
    ModelArrays      md_arrs;
    BrowBridgeArrays bb_arrs;

    readDataSet(scals,sobol_dirvcts, md_arrs, bb_arrs);
    scals.init();

    printf("Contract: %d, MC Its#: %d, #Underlyings: %d, #Path Dates: %d, #Sobol Bits: %d, #Models: %d, chunk: %d\n\n", 
            scals.contract, scals.num_mcits, scals.num_under, scals.num_dates, scals.sobol_bits, scals.num_models, scals.chunk);

    printf("#det pricers: %d, #cash flows: %d\n\n", scals.num_det_pricers, scals.num_cash_flows);

    run_CPUkernel( scals, sobol_dirvcts, md_arrs, bb_arrs );
#if 0
    printf("Row 0 of sobol direction vectors:\n[ ");
    for(int i=0; i<scals.sobol_bits; i++) {
        printf("%d, ", sobol_dirvcts[i]);
    }
    printf(" ]\n\n");
#endif

    return 1;
}

