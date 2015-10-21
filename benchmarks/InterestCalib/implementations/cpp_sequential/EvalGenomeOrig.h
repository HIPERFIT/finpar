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
                    const real_t&   a, 
                    const real_t&   b, 
                    const real_t&   rho, 
                    const real_t&   nu, 
                    const real_t&   sigma,
                    const real_t*   swaption,
                    IntermElem*   tmp_arrs,
                          real_t&   new_quote, // output
                          real_t&   new_price  // output
) {
    bool sanity = true;

    const real_t swap_freq  = swaption[1];
    const real_t maturity   = add_years( TODAY, swaption[0] );
    const UINT n_schedi   = static_cast<UINT>(12.0 * swaption[2] / swap_freq);
    const real_t tmat0      = date_act_365( maturity, TODAY );

    real_t strike;   
    { // BLACK PRICE computation
        real_t lvl = 0.0, t0 = MAX_DATE, tn = MIN_DATE; 

        for(UINT i = 0; i < n_schedi; i++) {  // reduce o map 
            // Map computes a1 and a2 (depends on i)
            const real_t a1 = add_months( maturity, swap_freq*i );
            const real_t a2 = add_months( a1,       swap_freq   );

            // Reduction( lvl: +, t0 : min, tn : max )
            lvl += zc(a2) * date_act_365(a2, a1);
            t0   = std::min(t0, a1);
            tn   = std::max(tn, a2);
        }

        strike = ( zc(t0) - zc(tn) ) / lvl;
        const real_t d1 = 0.5 * swaption[3] * tmat0;
        new_quote = lvl * strike * ( uGaussian_P(d1) - uGaussian_P(-d1) );
    } // END BLACK PRICE

    { // PRICER OF SWAPTION COMPUTATION
        real_t v0_mat, dummy1, dummy2;
        bigv( a, b, rho, nu, sigma, tmat0, v0_mat, dummy1, dummy2);
//
        const real_t mux = - bigmx( a, b, rho, nu, sigma, TODAY, maturity, TODAY, maturity );
        const real_t muy = - bigmy( a, b, rho, nu, sigma, TODAY, maturity, TODAY, maturity );
//
        const real_t zc_mat = zc(maturity);
//
        const real_t sqrt_bfun_a = sqrt( b_fun(2.0*a, tmat0) );
        const real_t sqrt_bfun_b = sqrt( b_fun(2.0*b, tmat0) );
        const real_t rhoxy  = rho * b_fun(a+b, tmat0) / (sqrt_bfun_a * sqrt_bfun_b);
        const real_t sigmax = sigma * sqrt_bfun_a;
        const real_t sigmay = nu    * sqrt_bfun_b;

        const real_t rhoxyc  = 1.0 - rhoxy * rhoxy;  // used in reduction kernel
        const real_t rhoxycs = sqrt( rhoxyc );       // used in reduction kernel
        const real_t sigmay_rhoxycs = sigmay * rhoxycs;
        const real_t t4      = (rhoxy * sigmay) / sigmax;

        for( UINT i = 0; i < n_schedi; i++ ) {
            const real_t beg_date = add_months( maturity, swap_freq*i ); //scheduleix[i];
            const real_t end_date = add_months( beg_date, swap_freq   ); //scheduleiy[i];
            const real_t res      = date_act_365( end_date, beg_date ) * strike;

            const real_t cii = ( i == n_schedi-1 ) ?  1.0 + res : res;
        
            real_t v0_end, vt_end, baii, bbii, date_tod1, date_tod2;
            date_tod1 = date_act_365(end_date, TODAY);
            bigv( a, b, rho, nu, sigma, date_tod1, v0_end, dummy1, dummy2 );
            date_tod2 = date_act_365(end_date, maturity);
            bigv( a, b, rho, nu, sigma, date_tod2, vt_end, baii, bbii );

            const real_t expo_aici = 0.5 * (vt_end - v0_end + v0_mat);
            const real_t fact_aici = cii * zc(end_date) / zc_mat;
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

        const real_t eps = 0.5 * sigmax;

        const real_t f   = exactYhatOrig( n_schedi, b, sigmax, sigmay, rhoxy, rhoxyc, rhoxycs, 
                                        mux, muy, tmp_arrs, mux       );

        const real_t g   = exactYhatOrig( n_schedi, b, sigmax, sigmay, rhoxy, rhoxyc, rhoxycs, 
                                        mux, muy, tmp_arrs, mux + eps );

        const real_t h   = exactYhatOrig( n_schedi, b, sigmax, sigmay, rhoxy, rhoxyc, rhoxycs, 
                                        mux, muy, tmp_arrs, mux - eps );

        const real_t df  = 0.5 * ( g - h ) / eps;

        const real_t sqrt2sigmax    = sqrt(2.0) * sigmax;;
        const real_t t2             = rhoxy / (sigmax*rhoxycs);;

        real_t accum = 0.0;

        for( UINT j = 0; j < NUM_HERMITE; j++ ) {
            const real_t x_quad = HermiteCoeffs [j];
            const real_t w_quad = HermiteWeights[j];

            const real_t      x = sqrt2sigmax * x_quad + mux;        
            const real_t yhat_x = f + df*(x - mux);
            const real_t h1     = ( (yhat_x - muy) / sigmay_rhoxycs ) - t2*( x - mux );

            real_t accum1 = 0.0;
            for( UINT i = 0; i < n_schedi; i++ ) {
                const real_t h2  = h1 + tmp_arrs[i].bbi * sigmay_rhoxycs;
            
                const real_t expo_aici = tmp_arrs[i].t1_cs + tmp_arrs[i].scale*x;
                const real_t fact_aici = tmp_arrs[i].ci;
                const real_t expo_part = uGaussian_P_withExpFactor( -h2, expo_aici );
                accum1 += fact_aici * expo_part;
            }

            sanity = ! ( isnan(accum1) || isinf(accum1) );
            assert(sanity && "Nan accum1 in pricer of swaption. Exiting!\n"); 

            real_t tmp = sqrt(2.0) * x_quad; //(x - mux) / sigmax;
            const real_t t1     = exp( - 0.5 * tmp * tmp );
        
            accum += w_quad * t1 * ( uGaussian_P(-h1) - accum1 );   
        }

        sanity = ! ( isnan(accum) || isinf(accum) );
        assert(sanity && "Nan accum1 in pricer of swaption. Exiting!\n"); 

        new_price = zc_mat * ( accum / sqrt( PI ) );
    }
}



#endif // end ifndef EVAL_GENOME_INLINED

