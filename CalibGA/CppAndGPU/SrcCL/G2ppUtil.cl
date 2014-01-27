////////////////////////////////////////////////////////////////////
///// Trivial (w.r.t. parallelism) Helper Functions for G2PP
////////////////////////////////////////////////////////////////////

#ifndef G2PP_UTIL
#define G2PP_UTIL

////////////////////////////////////
///   MATH RELATED
////////////////////////////////////

//-------------------------------------------------------------------------
// polynomial expansion of the erf() function, with error<=1.5e-7
//   formula 7.1.26 (page 300), Handbook of Mathematical Functions, Abramowitz and Stegun
//   http://people.math.sfu.ca/~cbm/aands/frameindex.htm

inline
REAL erff1( REAL x ) {
    REAL poly = 0.0, t, t2;

    t       = 1.0 / (1.0 + 0.3275911*x);
    poly   +=   0.254829592 * t;

    t2      = t * t;
    poly   += (-0.284496736) * t2;
    poly   +=   1.421413741 * t * t2;

    t2      = t2 * t2; 
    poly   += (-1.453152027) * t2;
    poly   +=   1.061405429  * t * t2;

    return ( 1.0 - poly * exp(-(x*x)) );

//    REAL p = 0.3275911;
//    REAL a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429;
//
//    REAL t  = 1.0/(1.0+p*x);
//    REAL t2 = t  * t;
//    REAL t3 = t  * t2;
//    REAL t4 = t2 * t2;
//    REAL t5 = t2 * t3;
//    return ( 1.0 - (a1*t + a2*t2 + a3*t3 + a4*t4 + a5*t5) * exp(-(x*x)) );

}

//-------------------------------------------------------------------------
// Cumulative Distribution Function for a standard normal distribution

inline
REAL uGaussian_P( REAL x ) {
    REAL u = x / sqrt(2.0);
    REAL e = ( u < 0.0 ) ? -erff1(-u) : erff1(u);

    return ( 0.5 * (1.0 + e) );
}

inline
REAL erff_poly_only( REAL x ) {
    REAL poly = 0.0, t, t2;

    t     = 1.0 / (1.0 + 0.3275911*x);
    poly +=   0.254829592 * t;

    t2    = t * t;
    poly += (-0.284496736) * t2;
    poly +=   1.421413741 * t * t2;

    t2    = t2 * t2; 
    poly += (-1.453152027) * t2;
    poly +=   1.061405429  * t * t2;
    
    return poly;
}

inline
REAL uGaussian_P_withExpFactor( REAL x, REAL exp_factor ) {
    REAL u = fabs( x / sqrt(2.0) );
    REAL e = erff_poly_only(u);

    REAL res = 0.5 * e * exp(exp_factor-u*u);

    if( x >= 0.0 ) {
        res = exp(exp_factor) - res;
    }

    return res;
}

/////////////////////////
//// G2PP Utilities
/////////////////////////

inline
REAL zc( REAL t ) {
    return exp( - R * date_act_365( t, TODAY ) );
}

inline
REAL b_fun( REAL z0, REAL tau ) {
    return ( 1.0 - exp(-z0*tau) ) / z0;
}

inline
REAL t_fun( REAL sigma, REAL x0, REAL tau ) {
    REAL expxtau  = exp( -x0*tau ) ;
    REAL exp2xtau = expxtau*expxtau;

    return sigma*sigma/(x0*x0) * ( tau + 2.0/x0*expxtau-1.0/(2.0*x0)*exp2xtau-3.0/(2.0*x0) );
}


///////////////////////////////////////////////////////////////////
// the first parameter `genome' is the five-genes genome used in
//     the genetic algorithms that models the interest rate
// the second parameter is the time
// the result is V in Brigo and Mercurio's book page 148.
//     \var(\int_t^T [x(u)+y(u)]du)
///////////////////////////////////////////////////////////////////
//res = (v, bai, bbi, _)
inline
REAL4 bigv(  REAL g_a, 
             REAL g_b, 
             REAL g_rho, 
             REAL g_nu, 
             REAL g_sig,
             REAL tau
) {
    REAL4 res4;

    res4.w = g_sig; //(g_sig == 0.0) ? 1.0e-10 : g_sig; // g_sigma
    res4.y = b_fun(g_a,        tau);
    res4.z = b_fun(g_b,        tau);
    res4.x = t_fun(res4.w, g_a,tau);
    res4.x+= t_fun(g_nu,   g_b,tau);
    res4.x+= 2.0 * g_rho * g_nu * res4.w / (g_a * g_b) *
             ( tau - res4.y - res4.z + b_fun(g_a+g_b, tau) );

    return res4;
}


///////////////////////////////////////////////////////////////////
// the first parameter `genome' is the five-genes genome used in
//     the genetic algorithms that models the interest rate
// the other parameter are times: today, maturity, and the
//      lower and upper bound of the considered time interval
//
// the result is: x drift term in tmat-forward measure
///////////////////////////////////////////////////////////////////
inline
REAL bigmx( REAL a, 
            REAL b, 
            REAL rho, 
            REAL nu, 
            REAL sigma,  // ends genome
            REAL today,
            REAL tmat,
            REAL s,
            REAL t
) {
    REAL ts    = date_act_365(t,    s)    ;
    REAL tmatt = date_act_365(tmat, t)    ;

    REAL tmat0 = date_act_365(tmat, today);
    REAL tmats = date_act_365(tmat, s)    ;
    REAL t0    = date_act_365(t,    today);
    REAL s0    = date_act_365(s,    today);

    REAL tmp, res = 0.0;

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
inline
REAL bigmy( REAL a, 
            REAL b, 
            REAL rho, 
            REAL nu, 
            REAL sigma,  // ends genome
            REAL today,
            REAL tmat,
            REAL s,
            REAL t
) {
    REAL ts    = date_act_365(t,    s    );
    REAL tmatt = date_act_365(tmat, t    );
    REAL tmat0 = date_act_365(tmat, today);
    REAL tmats = date_act_365(tmat, s    );
    REAL t0    = date_act_365(t,    today);
    REAL s0    = date_act_365(s,    today);

    REAL tmp, res = 0.0;

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

#endif // ifndef G2PP_UTIL

