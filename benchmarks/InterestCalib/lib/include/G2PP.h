#ifndef G2PP
#define G2PP

#include "Constants.h"
#include "Date.h"
#include "MathModule.h"
#include "G2ppUtil.h"

#include "GenAlgUtil.h"

//////////////////////////
// Root finder
//////////////////////////

real_t exactYhat( const UINT& n_schedi,

                const real_t& b,         // scals begins
                const real_t& sigmax,
                const real_t& sigmay,
                const real_t& rhoxy,
                const real_t& rhoxyc,
                const real_t& rhoxycs,
                const real_t& mux,
                const real_t& muy,      // scals ends

                const real_t* bai,      // babaicis begins
                const real_t* bbi,
                const real_t* aici,
                const real_t* log_aici, // babaicis ends

                      real_t* scales,
                const real_t& x
) {
    // ugaussian_Pinv(k)=1.0e~4
    const real_t k = - 3.71901648545568;

    real_t up = 0.0, lo = -INFTY;
    for( UINT i = 0; i < n_schedi; i++ ) {
        real_t baix = bai[i] * x;

        real_t up_term = aici[i] * exp( -baix );
        scales[i]    = up_term;
        up += up_term;
        lo  = std::max( lo, ( log_aici[i] - baix ) / bbi[i] );
    }

//    if ( n_schedi == 1 ) { return lo; }    // delete[] scales; 

//    return up;

// CHECKING uplo!!!!

    const real_t log_s = log(up);
    real_t tmp   = log_s / bbi[n_schedi-1];

    if ( tmp <= 0.0 ) {
        up = tmp;
    } else {
        tmp = log_s / bbi[0];
        if ( 0.0 <= tmp ) up = tmp;
        else              up = - INFTY;
    }

    const real_t yl = lo - EPS;
    const real_t yu = up + EPS;

    const real_t y0 = sigmay * ( rhoxy * (x-mux) / sigmax + k * rhoxycs ) - rhoxyc/b + muy;
    const real_t y1 = sigmay * ( rhoxy * (x-mux) / sigmax - k * rhoxycs )            + muy;

    real_t res;
    if      ( y1 <= yl ) res = y1 + 1.0;  // yhat is greater than y1 => 1 - phi(h_i(x, yhat)) < EPS
    else if ( yu <= y0 ) res = y0 - 1.0;  // yhat is lower than y0 => phi(h_i(x, yhat)) < EPS)
    else {        
        const real_t root_lb = std::max( yl, y0 );
        const real_t root_ub = std::min( yu, y1 );

        real_t root, error; UINT iter;
        rootFinding_Brent(1, n_schedi, scales, bbi, root_lb, root_ub, 1.0e-4, 1000, root, iter, error);
        //rootBisection(1, n_schedi, scales, bbi, root_lb, root_ub, 1.0e-4, 1000, root, iter, error);

        res = ( error == -INFTY ) ?  y0 - 1.0 : ( error ==  INFTY ) ? y1 + 1.0 : root;
    }

    return res;
}



///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
/// Main function of Module G2PP: pricer_of_swaption    ///
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

//////////////
// def a_fun(end_date,a,b,rho,nu,sigma,today,maturity,zc_mat,v0_mat):
//   # Brigo and Mercurio: defined top p. 148
//   v0_end,dummyA,dummyB=bigv(a,b,rho,nu,sigma,tau=Date.act_365(end_date,today))
//   vt_end,ba,bb=bigv(a,b,rho,nu,sigma,tau=Date.act_365(end_date,maturity))
//   res=zc(end_date)/zc_mat*N.exp(0.5*(vt_end-v0_end+v0_mat))
//   return res,ba,bb
//////////////

real_t pricer_of_swaption(    const real_t& today,
                            
                            const real_t& sw_mat,   // swaptionQuote begins ...
                            const real_t& sw_freq,
                            const real_t& sw_ty,    // swaptionQuote ends.

                            const real_t& a,        // genome begins ...
                            const real_t& b, 
                            const real_t& rho, 
                            const real_t& nu, 
                            const real_t& sigma     // genome ends.
) {
    //SwapOfSwap sos = extended_swaption_of_swaption(sw_mat, sw_freq, sw_ty);
    const real_t maturity  = add_years( TODAY, sw_mat );
    const UINT n_schedi  = static_cast<UINT>(12.0 * sw_ty / sw_freq);
//
    const real_t tmat0    = date_act_365( maturity, today ); // BIG BUG -- was TODAY
    real_t v0_mat, dummy1, dummy2;
    bigv( a, b, rho, nu, sigma, tmat0, v0_mat, dummy1, dummy2);
//
    const real_t mux     = - bigmx( a, b, rho, nu, sigma, today, maturity, today, maturity );
    const real_t muy     = - bigmy( a, b, rho, nu, sigma, today, maturity, today, maturity );
//
    const real_t zc_mat   = zc(maturity);
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
//
    real_t lvl = 0.0, t0 = MAX_DATE, tn = MIN_DATE; 
    for( UINT i = 0; i < n_schedi; i++ ) {
        real_t a1 = add_months( maturity, sw_freq*i );
        //sos.swap_sched1[i] = a1;
        real_t a2 = add_months( a1, sw_freq );
        //sos.swap_sched2[i] = a2;

        // Reduction( lvl: +, t0 : min, tn : max )
        lvl += zc(a2) * date_act_365(a2, a1);
        t0   = std::min(t0, a1);
        tn   = std::max(tn, a2);
    }

    const real_t strike   = ( zc(t0) - zc(tn) ) / lvl;
//
    real_t* ci       = new real_t[n_schedi];
    real_t* bai      = new real_t[n_schedi];
    real_t* bbi      = new real_t[n_schedi];
    real_t* aici     = new real_t[n_schedi];
    real_t* log_aici = new real_t[n_schedi];
    real_t* t1_cs    = new real_t[n_schedi];
    real_t* scale    = new real_t[n_schedi];
    real_t* hat_scale= new real_t[n_schedi];

    for( UINT i = 0; i < n_schedi; i++ ) {
        const real_t beg_date = add_months( maturity, sw_freq*i ); //scheduleix[i];
        const real_t end_date = add_months( beg_date, sw_freq   ); //scheduleiy[i];
        const real_t res      = date_act_365( end_date, beg_date ) * strike;

        const real_t cii = ( i == n_schedi-1 ) ?  1.0 + res : res;

        real_t v0_end, vt_end, baii, bbii, date_tod1, date_tod2;
        date_tod1 = date_act_365(end_date, today);
        bigv( a, b, rho, nu, sigma, date_tod1, v0_end, dummy1, dummy2 );
        date_tod2 = date_act_365(end_date, maturity);
        bigv( a, b, rho, nu, sigma, date_tod2, vt_end, baii, bbii );

        const real_t expo_aici = 0.5 * (vt_end - v0_end + v0_mat);
        const real_t fact_aici = cii * zc(end_date) / zc_mat;
        ci      [i] = fact_aici; // reuse the space to hold the factor of t1_cs
        t1_cs   [i] = bbii * (mux * t4 - (muy - 0.5*rhoxyc*sigmay*sigmay*bbii) ) + expo_aici; // hold only the exponent of the original t1_cs;

        bai     [i] = baii;
        bbi     [i] = bbii;
        aici    [i] = fact_aici * exp( expo_aici );
        log_aici[i] = log( fact_aici ) + expo_aici;
        scale   [i] = - ( baii + bbii * t4 );  

        if( isinf(aici[i]) || isnan(aici[i]) ) {
            fprintf(stderr, "NaN aici in pricer of swaption. Exiting!\n"); exit(0);
        }
    }
//

    const real_t eps = 0.5 * sigmax;

    const real_t f   = exactYhat( n_schedi, b, sigmax, sigmay, rhoxy, rhoxyc, rhoxycs, mux, muy, bai, bbi, aici, log_aici, hat_scale, mux       );
    const real_t g   = exactYhat( n_schedi, b, sigmax, sigmay, rhoxy, rhoxyc, rhoxycs, mux, muy, bai, bbi, aici, log_aici, hat_scale, mux + eps );
    const real_t h   = exactYhat( n_schedi, b, sigmax, sigmay, rhoxy, rhoxyc, rhoxycs, mux, muy, bai, bbi, aici, log_aici, hat_scale, mux - eps );
    const real_t df  = 0.5 * ( g - h ) / eps;

    // AT this point we need:
    //    scalars: sqrt2sigmax, t2, f, df, mux, muy, sigmay_rhoxycs, zc_mat
    //    arrays : t1_cs, bbi, scale, ci 

    const real_t sqrt2sigmax = sqrt(2.0) * sigmax;
    const real_t t2      = rhoxy / (sigmax*rhoxycs);
    // t2 * (sqrt2sigmax * x_quad) -> sqrt(2.0) * rhoxy / rhoxycs * xquad

    real_t accum0 = 0.0;
    for( UINT j = 0; j < NUM_HERMITE; j++ ) {
        const real_t x_quad = HermiteCoeffs [j];
        const real_t w_quad = HermiteWeights[j];

        const real_t x = sqrt2sigmax * x_quad + mux;

        const real_t yhat_x = f + df*(x - mux);
        const real_t h1     = ( (yhat_x - muy) / sigmay_rhoxycs ) - t2*( x - mux );

        real_t accum1 = 0.0;
        for( UINT i = 0; i < n_schedi; i++ ) {
            const real_t h2  = h1 + bbi[i] * sigmay_rhoxycs;

            const real_t expo_aici = t1_cs[i] + scale[i]*x;
            const real_t fact_aici = ci[i];
//            accum1 += fact_aici * exp(expo_aici) * uGaussian_P(-h2);
            const real_t expo_part = uGaussian_P_withExpFactor( -h2, expo_aici );
            accum1 += fact_aici * expo_part;

            if( isnan(accum1) || isinf(accum1) ) {
                fprintf(stderr, "accum is Nan -- expo_aici: %f, fact_aici: %f, h2: %f, t1_cs: %f, scale: %f, x: %f, bbi: %f, mux: %f, muy: %f\n", 
                        expo_aici, fact_aici, h2, t1_cs[i], scale[i], x, bbi[i], mux, muy );
                fprintf(stderr, "accum is Nan -- t4: %f, rhoxyc: %f, v0mat: %f, b: %f, aici: %f\n", 
                            t4, rhoxyc, v0_mat, b, aici[i] );
                exit(0);
            }
        }

        real_t tmp = sqrt(2.0) * x_quad; //(x - mux) / sigmax;
        const real_t t1     = exp( - 0.5 * tmp * tmp );

        accum0 += w_quad * t1 * ( uGaussian_P(-h1) - accum1 );     
    }

    if( isnan(accum0) || isinf(accum0) ) {
        fprintf(stderr,"Numerically Instable Implem: end of pricer_of_swaption: accum0 in NaN or Inf!\n");
        exit(0);
    }

    delete[] ci;  delete[] hat_scale;
    delete[] bai; delete[] aici;  delete[] log_aici;
    delete[] bbi; delete[] t1_cs; delete[] scale;

    return zc_mat * ( accum0 / sqrt( PI ) );
}


////////////////////////////////
//// Pricer for one genome: ////
////////////////////////////////

real_t pricer (const real_t a = 0.02453, const real_t b = 0.98376, const real_t rho = -0.82400, const real_t nu = 0.11830, const real_t sigma = 0.02398) {
    //let genome = {0.02453, 0.98376, ~0.82400, 0.11830, 0.02398} in

    real_t rms = 0.0, logLik = 0.0;
    for( UINT i = 0; i < NUM_SWAP_QUOTES; i++ ) {
        const real_t mat_year  = SwaptionQuotes[4*i+0];
        const real_t swap_freq = SwaptionQuotes[4*i+1];
        const real_t term_year = SwaptionQuotes[4*i+2];
        const real_t quote     = SwaptionQuotes[4*i+3];

        const real_t g2pp_price   = pricer_of_swaption( TODAY, mat_year, swap_freq, term_year, a, b, rho, nu, sigma);
        const real_t market_price = black_price       ( TODAY, mat_year, swap_freq, term_year, quote               );

        const real_t tmp          = (g2pp_price - market_price) / market_price; 
        rms    += tmp * tmp;
        logLik += logLikelihood( market_price, g2pp_price );
    }

    rms = 100.0 * sqrt ( rms / NUM_SWAP_QUOTES );

    fprintf(stderr, "\n\n!!!!!!!!!!!!!Computed RMS is: %f, Computed LogLik is: %f\n\n", rms, logLik); 

    return rms;
}


//////////////////////////////////////////////////////
//// Test pricer_of_swaption:
////
//// params=params2dict(a = 0.02453, b = 0.98376, sigma = 0.02398, nu = 0.11830, rho = -0.82400)
//// swaption=['swaption_maturity_in_year': 10, 'swap_term_in_year': 4, 'swap_frequency': 6]
//// assert "%.3f" % (1e4*pricer_of_swaption(today=today,zc=zc,swaption=swaption,params=params)) == "657.822"
//// swaption=['swaption_maturity_in_year': 30, 'swap_term_in_year': 30, 'swap_frequency': 6]
//// assert "%.3f" % (1e4*pricer_of_swaption(today=today,zc=zc,swaption=swaption,params=params)) == "1902.976"
////
//////////////////////////////////////////////////////


void test_pricer_of_swaption() {
    const real_t a        = 0.02453;
    const real_t b        = 0.98376;
    const real_t rho      = -0.82400;
    const real_t nu       = 0.11830;
    const real_t sigma    = 0.02398;

    // (maturity, frequency, term) = swaption
    real_t maturity = 10.0;
    real_t freq     =  6.0;
    real_t term     =  4.0;

    const real_t price1   = 1.0e4 * pricer_of_swaption( TODAY, maturity, freq, term, a, b, rho, nu, sigma);
    printf("Pricer_of_swaption test I = 657.82158867845 :   ");
    if( equalEps(price1, 657.82158867845) ) printf(" SUCCESS!\n\n");
    else                                    printf("%f FAILS!\n\n", price1);

    // let swaption = {30.0, 6.0, 30.0}
    maturity = 30.0;
    freq     =  6.0;
    term     = 30.0;

    const real_t price2   = 1.0e4 * pricer_of_swaption( TODAY, maturity, freq, term, a, b, rho, nu, sigma );
    printf("Pricer_of_swaption test II = 1902.97628191498 :   ");
    if( equalEps(price2, 1902.97628191498) ) printf(" SUCCESS!\n\n");
    else                                     printf("%f FAILS!\n\n", price2);


    //let swaption = {30.0, 6.0, 25.0}                                   in
    maturity = 30.0;
    freq     =  6.0;
    term     = 25.0;

    const real_t price3   = 1.0e4 * pricer_of_swaption( TODAY, maturity, freq, term, a, b, rho, nu, sigma );
    printf("Pricer_of_swaption test III = 1840.859126408099 :   ");
    if( equalEps(price3, 1840.859126408099) ) printf(" SUCCESS!\n\n");
    else                                      printf("%f FAILS!\n\n", price3);
}

#endif // ifndef G2PP
