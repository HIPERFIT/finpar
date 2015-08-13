#ifndef EVAL_GENOME_INLINED
#define EVAL_GENOME_INLINED

#include "Constants.h"
#include "GenAlgUtil.h"
#include "G2PP.h"
#include "G2PPorig.h"


/**
 * MOST IMPORTANTLY: GENOME EVALUATION By Pricer of Swaption & BLACK PRICE
 */
void eval_genome_new (   
                    const REAL&   a, 
                    const REAL&   b, 
                    const REAL&   rho, 
                    const REAL&   nu, 
                    const REAL&   sigma,
                    const REAL*   swaption,
                    IntermElem*   tmp_arrs,
                          REAL&   new_quote, // output
                          REAL&   new_price  // output
) {
    bool sanity = true;

    const REAL swap_freq  = swaption[1];
    const REAL maturity   = add_years( TODAY, swaption[0] );
    const UINT n_schedi   = static_cast<UINT>(12.0 * swaption[2] / swap_freq);
    const REAL tmat0      = date_act_365( maturity, TODAY );

    REAL strike;   
    { // BLACK PRICE computation
        REAL lvl = 0.0, t0 = MAX_DATE, tn = MIN_DATE; 

        for(UINT i = 0; i < n_schedi; i++) {  // reduce o map 
            // Map computes a1 and a2 (depends on i)
            const REAL a1 = add_months( maturity, swap_freq*i );
            const REAL a2 = add_months( a1,       swap_freq   );

            // Reduction( lvl: +, t0 : min, tn : max )
            lvl += zc(a2) * date_act_365(a2, a1);
            t0   = std::min(t0, a1);
            tn   = std::max(tn, a2);
        }

        strike = ( zc(t0) - zc(tn) ) / lvl;
        const REAL d1 = 0.5 * swaption[3] * tmat0;
        new_quote = lvl * strike * ( uGaussian_P(d1) - uGaussian_P(-d1) );
    } // END BLACK PRICE

    { // PRICER OF SWAPTION COMPUTATION
        REAL v0_mat, dummy1, dummy2;
        bigv( a, b, rho, nu, sigma, tmat0, v0_mat, dummy1, dummy2);
//
        const REAL mux = - bigmx( a, b, rho, nu, sigma, TODAY, maturity, TODAY, maturity );
        const REAL muy = - bigmy( a, b, rho, nu, sigma, TODAY, maturity, TODAY, maturity );
//
        const REAL zc_mat = zc(maturity);
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
            tmp_arrs[i].ci       = fact_aici; // reuse the space to hold the factor of t1_cs
            tmp_arrs[i].t1_cs    = bbii * (mux * t4 - (muy - 0.5*rhoxyc*sigmay*sigmay*bbii) ) + expo_aici; 
                // hold only the exponent of the original t1_cs;

            tmp_arrs[i].bai      = baii;
            tmp_arrs[i].bbi      = bbii;
            tmp_arrs[i].aici     = fact_aici * exp( expo_aici );
            tmp_arrs[i].log_aici = log( fact_aici ) + expo_aici;
            tmp_arrs[i].scale    = - ( baii + bbii * t4 );

            sanity = ! ( isinf(tmp_arrs[i].aici) || isnan(tmp_arrs[i].aici) );
            assert(sanity && "Nan aici in pricer of swaption. Exiting!\n"); 
        }

        const REAL eps = 0.5 * sigmax;

        const REAL f   = exactYhatOrig( n_schedi, b, sigmax, sigmay, rhoxy, rhoxyc, rhoxycs, 
                                        mux, muy, tmp_arrs, mux       );

        const REAL g   = exactYhatOrig( n_schedi, b, sigmax, sigmay, rhoxy, rhoxyc, rhoxycs, 
                                        mux, muy, tmp_arrs, mux + eps );

        const REAL h   = exactYhatOrig( n_schedi, b, sigmax, sigmay, rhoxy, rhoxyc, rhoxycs, 
                                        mux, muy, tmp_arrs, mux - eps );

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
                const REAL h2  = h1 + tmp_arrs[i].bbi * sigmay_rhoxycs;
            
                const REAL expo_aici = tmp_arrs[i].t1_cs + tmp_arrs[i].scale*x;
                const REAL fact_aici = tmp_arrs[i].ci;
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

        new_price = zc_mat * ( accum / sqrt( PI ) );
    }
}



#endif // end ifndef EVAL_GENOME_INLINED

