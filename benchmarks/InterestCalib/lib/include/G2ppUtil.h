#ifndef G2PP_UTIL
#define G2PP_UTIL


#include "Constants.h"
#include "MathModule.h"

real_t zc( const real_t& t ) {
    return exp( - R * date_act_365( t, TODAY ) );
}

/**
Triple<real_t> accumSched( Triple<real_t> xx, Triple<real_t> yy ) {
    real_t x = xx.fst, d1 = xx.snd, d2 = xx.thrd, y = yy.fst, tf = yy.snd, tp = yy.thrd;
    return Triple<real_t>( x + zc(tp) * date_act_365(tp, tf), std::min(d1,tf), std::max(d2,tp) );
}
**/


SwapOfSwap extended_swaption_of_swaption( 
            const real_t& sw_mat, 
            const real_t& freq, 
            const real_t& sw_ty  
) { 
    const real_t maturity  = add_years( TODAY, sw_mat );
    const UINT nschedule  = static_cast<UINT>(12.0 * sw_ty / freq);

    SwapOfSwap sos(nschedule);

    real_t lvl = 0.0, t0 = MAX_DATE, tn = MIN_DATE; 
    for(UINT i=0; i<nschedule; i++) {
        real_t a1 = add_months( maturity, freq*i );
        sos.swap_sched1[i] = a1;

        real_t a2 = add_months( a1, freq );
        sos.swap_sched2[i] = a2;

        // Reduction( lvl: +, t0 : min, tn : max )
        lvl += zc(a2) * date_act_365(a2, a1);
        t0   = std::min(t0, a1);
        tn   = std::max(tn, a2);
    }

    sos.maturity = maturity;
    sos.strike   = ( zc(t0) - zc(tn) ) / lvl;

    return sos;
}


real_t b_fun( const real_t& z0, const real_t& tau ) {
    //const real_t z0 = (z == 0.0) ? EPS0 : z; 
    return ( 1.0 - exp(-z0*tau) ) / z0;
}

real_t t_fun( const real_t& sigma, const real_t& x0, const real_t& tau ) {
    //const real_t x0       = ( fabs(x) >= EPS ) ? x : ( x >= 0.0 ) ? EPS : -EPS; 
    const real_t expxtau  = exp( -x0*tau ) ;
    const real_t exp2xtau = expxtau*expxtau;

    return sigma*sigma/(x0*x0) * ( tau + 2.0/x0*expxtau-1.0/(2.0*x0)*exp2xtau-3.0/(2.0*x0) );
}


///////////////////////////////////////////////////////////////////
// the first parameter `genome' is the five-genes genome used in
//     the genetic algorithms that models the interest rate
// the second parameter is the time
// the result is V in Brigo and Mercurio's book page 148.
//     \var(\int_t^T [x(u)+y(u)]du)
///////////////////////////////////////////////////////////////////
void bigv(  const real_t& g_a, 
            const real_t& g_b, 
            const real_t& g_rho, 
            const real_t& g_nu, 
            const real_t& g_sig,
            const real_t& tau,
                  real_t& v, 
                  real_t& bai, 
                  real_t& bbi
) {
    const real_t g_sigma = g_sig; //(g_sig == 0.0) ? 1.0e-10 : g_sig;

    bai = b_fun(g_a,        tau);
    bbi = b_fun(g_b,        tau);
    v   = t_fun(g_sigma,g_a,tau);
    v  += t_fun(g_nu,   g_b,tau);
    v  += 2.0 * g_rho * g_nu * g_sigma / (g_a * g_b)*
             ( tau - bai - bbi + b_fun(g_a+g_b, tau) );
}


///////////////////////////////////////////////////////////////////
// the first parameter `genome' is the five-genes genome used in
//     the genetic algorithms that models the interest rate
// the other parameter are times: today, maturity, and the
//      lower and upper bound of the considered time interval
//
// the result is: x drift term in tmat-forward measure
///////////////////////////////////////////////////////////////////
real_t bigmx( const real_t& a, 
            const real_t& b, 
            const real_t& rho, 
            const real_t& nu, 
            const real_t& sigma,  // ends genome
            const real_t& today,
            const real_t& tmat,
            const real_t& s,
            const real_t& t
) {
    const real_t ts    = date_act_365(t,    s)    ;
    const real_t tmatt = date_act_365(tmat, t)    ;

    const real_t tmat0 = date_act_365(tmat, today);
    const real_t tmats = date_act_365(tmat, s)    ;
    const real_t t0    = date_act_365(t,    today);
    const real_t s0    = date_act_365(s,    today);

    real_t tmp, res = 0.0;

    tmp   = (sigma*sigma)/(a*a)+(sigma*rho*nu)/(a*b);
    tmp  *= ( 1.0 - exp(- a * ts) );
    res  += tmp;

    tmp   = sigma * sigma / (2.0 * a * a);
    tmp  *= exp(- a * tmatt) - exp(- a * (tmats + ts));
    res  -= tmp;

    tmp   = rho * sigma * nu / (b * (a + b));
    tmp  *= exp(-b * tmatt) - exp(-b*tmat0 - a*t0 + (a+b)*s0);
    res  -= tmp;

    return res;
}


///////////////////////////////////////////////////////////////////
// the first parameter `genome' is the five-genes genome used in
//     the genetic algorithms that models the interest rate
// the other parameter are times: today, maturity, and the
//      lower and upper bound of the considered time interval
//
// the result is: y drift term in tmat-forward measure
///////////////////////////////////////////////////////////////////

real_t bigmy( const real_t& a, 
            const real_t& b, 
            const real_t& rho, 
            const real_t& nu, 
            const real_t& sigma,  // ends genome
            const real_t& today,
            const real_t& tmat,
            const real_t& s,
            const real_t& t
) {
    const real_t ts    = date_act_365(t,    s)    ;
    const real_t tmatt = date_act_365(tmat, t)    ;
    const real_t tmat0 = date_act_365(tmat, today);
    const real_t tmats = date_act_365(tmat, s)    ;
    const real_t t0    = date_act_365(t,    today);
    const real_t s0    = date_act_365(s,    today);

    real_t tmp, res = 0.0;

    tmp  = nu*nu/(b*b)+sigma*rho*nu/(a*b);
    tmp *= 1.0 - exp(-b * ts);
    res += tmp;

    tmp  = nu * nu / (2.0 * b * b);
    tmp *= exp(-b * tmatt) - exp(-b * (tmats + ts));
    res -= tmp;

    tmp  = sigma * rho * nu / (a * (a + b));
    tmp *= exp(-a * tmatt) - exp(-a*tmat0 - b*t0 + (a+b)*s0);
    res -= tmp;

    return res;
}


real_t black_price (  const real_t& today, 
                    const real_t& sw_mat, // swaption.1
                    const real_t& freq,   // swaption.2
                    const real_t& sw_ty,  // swaption.3
                    const real_t& vol
) {  
    // inlined extended_swaption_of_swaptions
    const real_t maturity   = add_years( TODAY, sw_mat );
    const UINT nschedule  = static_cast<UINT>(12.0 * sw_ty / freq);
    const real_t sqrtt      = date_act_365( maturity, today );

    // morally equivalent to `swap_schedule2lvl(swap_schedule)' but in map-reduce form!!

    real_t lvl = 0.0, t0 = MAX_DATE, tn = MIN_DATE; 
    for(UINT i=0; i<nschedule; i++) {
        const real_t a1 = add_months( maturity, freq*i );
        const real_t a2 = add_months( a1, freq );

        // Reduction( lvl: +, t0 : min, tn : max )
        lvl += zc(a2) * date_act_365(a2, a1);
        t0   = std::min(t0, a1);
        tn   = std::max(tn, a2);
    }

    const real_t strike = ( zc(t0) - zc(tn) ) / lvl;
//    const real_t d1     = log( strike / strike ) / ( vol * sqrtt ) + 0.5 * vol * sqrtt;
    const real_t d1     = 0.5 * vol * sqrtt;
    const real_t d2     = 0.0 - d1; //d1 - vol * sqrtt;

    return ( lvl * ( strike * uGaussian_P(d1) - strike * uGaussian_P(d2) ) );
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
/// Testing g2pp minus the main function,
///    i.e., pricer_of_swaption
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void test_g2ppUtil() {
    real_t today = 9000.0;
    real_t tmat  = 18000.0;
    real_t s     = 400000.0;
    real_t t     = 9000000.0;

    ///////////////////////////////////////////
    // testing b_fun, bigv, bigmx, bigmy
    ///////////////////////////////////////////

    real_t res_b_fun = b_fun(3.24, 1.362);
    real_t res_bigmx = bigmx(0.02, 0.02, 0.0, 0.01, 0.04, today, tmat, s, t);
    real_t res_bigmy = bigmy(0.02, 0.02, 0.0, 0.01, 0.04, today, tmat, s, t);

    real_t res_bigv1, res_bigv2, res_bigv3;
    bigv (0.02, 0.02, 0.0, 0.01, 0.04, 1.12, res_bigv1, res_bigv2, res_bigv3);

    printf("b_fun test = 0.30490117 : ");
    if( equalEps(res_b_fun, 0.30490117) )           printf(" SUCCESS!\n\n");
    else                                            printf(" %f FAILS!\n\n", res_b_fun);


    printf("bigv test  = { 7.8288965347e-4, 1.107549139, 1.107549139 } : ");
    if( equalEps(res_bigv1, 7.8288965347e-4) &&
        equalEps(res_bigv2, 1.107549139    ) && 
        equalEps(res_bigv3, 1.107549139    )   )    printf(" SUCCESS!\n\n");
    else                                            printf(" { %f %f %f } FAILS!\n\n", res_bigv1, res_bigv2, res_bigv3);

    printf("bigmx test = -0.2356067470979 : ");
    if( equalEps(res_bigmx, -0.2356067470979) )     printf(" SUCCESS!\n\n");
    else                                            printf(" %f FAILS!\n\n", res_bigmx);

    printf("bigmy test = -0.01472542169362 : ");
    if( equalEps(res_bigmy, -0.01472542169362) )    printf(" SUCCESS!\n\n");
    else                                            printf(" %f FAILS!\n\n", res_bigmy);


    //////////////////////////////////
    // testing extended_swaption_of_swaption
    // The right value to test against is "654.142965".
    // However, because "today is not tomorrow", i.e., the date module
    // is very approximatively implemented, we test against ``655.250458''
    //////////////////////////////////

    //let swaption = {10.0, 6.0, 4.0}     in
    const real_t sw_mat = 10.0; 
    const real_t freq   =  6.0;
    const real_t sw_ty  =  4.0;

    real_t maturity          = 22094640.0;
    real_t strike            = 0.030226283149239714;
    real_t swap_schedule1[8] = { 22094640.0, 22355280.0, 22620240.0, 22880880.0, 23145840.0, 23407920.0, 23672880.0, 23933520.0 };
    real_t swap_schedule2[8] = { 22355280.0, 22620240.0, 22880880.0, 23145840.0, 23407920.0, 23672880.0, 23933520.0, 24198480.0 };
    //SwapOfSwap sos(8, 22094640.0, 0.030226283149239714, swap_schedule1, swap_schedule2);

    SwapOfSwap sos = extended_swaption_of_swaption(sw_mat, freq, sw_ty);

    bool mat_ok    = equalEps(maturity, sos.maturity );
    bool strike_ok = equalEps(strike,   sos.strike   );
    bool sched_ok  = true;
    for(UINT i=0; i<sos.n; i++) {
        sched_ok = sched_ok && equalEps( swap_schedule1[i], sos.swap_sched1[i] );
        sched_ok = sched_ok && equalEps( swap_schedule2[i], sos.swap_sched2[i] );
    }

    sos.cleanUp();

    printf("Testing extended_swaption_of_swaption: ");
    if(mat_ok && strike_ok && sched_ok) printf("SUCCESS!\n");
    else                                printf("FAILS!  \n");
    printf("\t\tmaturity: %f, strike: %f\n\n\n", sos.maturity, sos.strike);


    //////////////////////////////////
    // testing black_price
    // The right value to test against is "654.142965".
    // However, because "today is not tomorrow", i.e., the date module
    // is very approximatively implemented, we test against ``654.1689526995502''
    //////////////////////////////////

    const real_t vol       = 0.2454;
    real_t black_price_res = black_price( TODAY, sw_mat, freq, sw_ty, vol) * 10000.0;

    printf("\n\nTesting Black Price = 654.1429648454:  ");
    if( equalEps(black_price_res, 654.1429648454) ) printf(" SUCCESS!\n\n");
    else                                            printf("%f FAILS!\n\n", black_price_res);

    printf("\n\n\n");
}

#endif // end G2PP_UTIL
