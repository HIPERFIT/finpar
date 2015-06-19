#ifndef EVAL_GENOME_INLINED
#define EVAL_GENOME_INLINED

#include "Constants.h"
#include "GenAlgUtil.h"
#include "G2PP.h"


/**
 * MOST IMPORTANTLY: GENOME EVALUATION By Pricer of Swaption & BLACK PRICE
 */
REAL eval_genome(   const REAL a, 
                    const REAL b, 
                    const REAL rho, 
                    const REAL nu, 
                    const REAL sigma, 

                    REAL* anew_quote, // [NUM_SWAP_QUOTES]
                    REAL* anew_price, // [NUM_SWAP_QUOTES]

                    const int N, 
                    short* flags,      // [N]
 
                    int  * start_inds  // [NUM_SWAP_QUOTES]
) {
    REAL rms = 0.0;
    bool sanity = true;

    // expand the irregular arrays
    REAL *irreg_arrays = new REAL[8 * N];
    REAL *gbai, *gbbi, *gaici, *glog_aici, *gci, *gt1_cs, *gscale, *ghat_scale;
    { // map intermediary arrays
        REAL* ptr = irreg_arrays;
        gbai      = ptr; ptr += N;
        gbbi      = ptr; ptr += N;
        gaici     = ptr; ptr += N;
        glog_aici = ptr; ptr += N; 
        gci       = ptr; ptr += N;
        gt1_cs    = ptr; ptr += N;
        gscale    = ptr; ptr += N;
        ghat_scale= ptr; ptr += N;
    }

    // BIG KERNEL
    for( UINT ttt = 0; ttt < NUM_SWAP_QUOTES; ttt++ ) {
        const REAL mat_year  = SwaptionQuotes[4*ttt + 0];
        const REAL swap_freq = SwaptionQuotes[4*ttt + 1];
        const REAL term_year = SwaptionQuotes[4*ttt + 2];
        const REAL quote     = SwaptionQuotes[4*ttt + 3];

        // adjust temporary array
        const int  beg_ind = start_inds[ttt];
        REAL *ci = gci + beg_ind, *bai = gbai + beg_ind, *bbi = gbbi + beg_ind, *aici = gaici + beg_ind;
        REAL *log_aici = glog_aici + beg_ind, *t1_cs = gt1_cs + beg_ind, *scale = gscale + beg_ind;
        REAL *hat_scale = ghat_scale + beg_ind;

        // new_quote does not have to be computed all the time, does it?
        // it can just be computed once and indexed into an array ... 
        //const REAL new_quote = black_price ( TODAY, mat_year, swap_freq, term_year, quote );

        const REAL maturity   = add_years( TODAY, mat_year );
        const UINT n_schedi   = static_cast<UINT>(12.0 * term_year / swap_freq);
        const REAL tmat0      = date_act_365( maturity, TODAY );

        REAL strike;   // new_quote, 
        { // BLACK PRICE computation
            REAL lvl = 0.0, t0 = MAX_DATE, tn = MIN_DATE; 

            for(UINT i = 0; i < n_schedi; i++) {  // reduce o map => in local memory
                const REAL a1 = add_months( maturity, swap_freq*i );
                const REAL a2 = add_months( a1,       swap_freq   );

                // Reduction( lvl: +, t0 : min, tn : max )
                lvl += zc(a2) * date_act_365(a2, a1);
                t0   = std::min(t0, a1);
                tn   = std::max(tn, a2);
            }

            strike = ( zc(t0) - zc(tn) ) / lvl;
            const REAL d1     = 0.5 * quote * tmat0;
            anew_quote[ttt] = lvl * strike * ( uGaussian_P(d1) - uGaussian_P(-d1) );
        } // END BLACK PRICE


        //const REAL new_price = pricer_of_swaption( TODAY, mat_year, swap_freq, term_year, g_a, g_b, g_rho, g_nu, g_sigma);
        //REAL new_price;
        { // PRICER OF SWAPTION COMPUTATION
            REAL v0_mat, dummy1, dummy2;
            bigv( a, b, rho, nu, sigma, tmat0, v0_mat, dummy1, dummy2);
//
            const REAL mux     = - bigmx( a, b, rho, nu, sigma, TODAY, maturity, TODAY, maturity );
            const REAL muy     = - bigmy( a, b, rho, nu, sigma, TODAY, maturity, TODAY, maturity );
//
            const REAL zc_mat   = zc(maturity);
//
            const REAL sqrt_bfun_a = sqrt( b_fun(2.0*a, tmat0) );
            const REAL sqrt_bfun_b = sqrt( b_fun(2.0*b, tmat0) );
            const REAL rhoxy  = rho * b_fun(a+b, tmat0) / (sqrt_bfun_a * sqrt_bfun_b);
            const REAL sigmax = sigma * sqrt_bfun_a;
            const REAL sigmay = nu    * sqrt_bfun_b;

            const REAL rhoxyc  = 1.0 - rhoxy * rhoxy;  // used in reduction kernel
            const REAL rhoxycs = sqrt( rhoxyc );       // used in reduction kernel
            const REAL sigmay_rhoxycs = sigmay * rhoxycs;
            const REAL t4      = (rhoxy * sigmay) / sigmax;
//

            for( UINT i = 0; i < n_schedi; i++ ) {
                const REAL beg_date = add_months( maturity, swap_freq*i ); //scheduleix[i];
                const REAL end_date = add_months( beg_date, swap_freq   ); //scheduleiy[i];
                const REAL res      = date_act_365( end_date, beg_date ) * strike;

                const REAL cii = ( i == n_schedi-1 ) ?  1.0 + res : res;
        
                REAL v0_end, vt_end, baii, bbii, date_tod1, date_tod2;
                date_tod1 = date_act_365(end_date, TODAY);
                bigv( a, b, rho, nu, sigma, date_tod1, v0_end, dummy1, dummy2 );
                date_tod2 = date_act_365(end_date, maturity);
                bigv( a, b, rho, nu, sigma, date_tod2, vt_end, baii, bbii );

                const REAL expo_aici = 0.5 * (vt_end - v0_end + v0_mat);
                const REAL fact_aici = cii * zc(end_date) / zc_mat;
                ci      [i] = fact_aici; // reuse the space to hold the factor of t1_cs
                t1_cs   [i] = bbii * (mux * t4 - (muy - 0.5*rhoxyc*sigmay*sigmay*bbii) ) + expo_aici; 
                // hold only the exponent of the original t1_cs;

                bai     [i] = baii;
                bbi     [i] = bbii;
                aici    [i] = fact_aici * exp( expo_aici );
                log_aici[i] = log( fact_aici ) + expo_aici;
                scale   [i] = - ( baii + bbii * t4 );  

                sanity = ! ( isinf(aici[i]) || isnan(aici[i]) );
                assert(sanity && "Nan aici in pricer of swaption. Exiting!\n"); 
            }
//
            const REAL eps = 0.5 * sigmax;

            const REAL f   = exactYhat( n_schedi, b, sigmax, sigmay, rhoxy, rhoxyc, rhoxycs, 
                                        mux, muy, bai, bbi, aici, log_aici, hat_scale, mux       );
            const REAL g   = exactYhat( n_schedi, b, sigmax, sigmay, rhoxy, rhoxyc, rhoxycs, 
                                        mux, muy, bai, bbi, aici, log_aici, hat_scale, mux + eps );
            const REAL h   = exactYhat( n_schedi, b, sigmax, sigmay, rhoxy, rhoxyc, rhoxycs, 
                                        mux, muy, bai, bbi, aici, log_aici, hat_scale, mux - eps );
            const REAL df  = 0.5 * ( g - h ) / eps;

            const REAL sqrt2sigmax    = sqrt(2.0) * sigmax;;
            const REAL t2             = rhoxy / (sigmax*rhoxycs);;

            REAL accum = 0.0;
            
            for( UINT j = 0; j < NUM_HERMITE; j++ ) {
                const REAL x_quad = HermiteCoeffs [j];
                const REAL w_quad = HermiteWeights[j];

                const REAL      x = sqrt2sigmax * x_quad + mux;        
                const REAL yhat_x = f + df*(x - mux);
                const REAL h1     = ( (yhat_x - muy) / sigmay_rhoxycs ) - t2*( x - mux );

                REAL accum1 = 0.0;
                for( UINT i = 0; i < n_schedi; i++ ) {
                    const REAL h2  = h1 + bbi[i] * sigmay_rhoxycs;
            
                    const REAL expo_aici = t1_cs[i] + scale[i]*x;
                    const REAL fact_aici = ci[i];
    //              accum1 += fact_aici * exp(expo_aici) * uGaussian_P(-h2);
                    const REAL expo_part = uGaussian_P_withExpFactor( -h2, expo_aici );
                    accum1 += fact_aici * expo_part;
                }

                sanity = ! ( isnan(accum1) || isinf(accum1) );
                assert(sanity && "Nan accum1 in pricer of swaption. Exiting!\n"); 

                REAL tmp = sqrt(2.0) * x_quad; //(x - mux) / sigmax;
                const REAL t1     = exp( - 0.5 * tmp * tmp );
        
                accum += w_quad * t1 * ( uGaussian_P(-h1) - accum1 );   
            }

            sanity = ! ( isnan(accum) || isinf(accum) );
            assert(sanity && "Nan accum1 in pricer of swaption. Exiting!\n"); 

            anew_price[ttt] = zc_mat * ( accum / sqrt( PI ) );

//        const REAL tmp       = (anew_price[ttt] - anew_quote[ttt]) / anew_quote[ttt]; 
//        rms += tmp * tmp;
            const REAL lik = logLikelihood( anew_quote[ttt], anew_price[ttt] );
            rms += lik;
        }
    }

    delete[] irreg_arrays; 

    return rms;
//    return 0.0-rms;
}

#endif // end ifndef EVAL_GENOME_INLINED

