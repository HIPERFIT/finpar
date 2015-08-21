#ifndef MATH_MODULE
#define MATH_MODULE

#include "Constants.h"

////////////////////////////////////////////////////////////////
///// MATH MODULE
////////////////////////////////////////////////////////////////

//-------------------------------------------------------------------------
// polynomial expansion of the erf() function, with error<=1.5e-7
//   formula 7.1.26 (page 300), Handbook of Mathematical Functions, Abramowitz and Stegun
//   http://people.math.sfu.ca/~cbm/aands/frameindex.htm

real_t erff1( const real_t& x ) {
    const real_t p = 0.3275911;
    const real_t a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429;

    const real_t t  = 1.0/(1.0+p*x);
    const real_t t2 = t  * t;
    const real_t t3 = t  * t2;
    const real_t t4 = t2 * t2;
    const real_t t5 = t2 * t3;

    return ( 1.0 - (a1*t + a2*t2 + a3*t3 + a4*t4 + a5*t5) * exp(-(x*x)) );
}

//-------------------------------------------------------------------------
// Cumulative Distribution Function for a standard normal distribution

real_t uGaussian_P( const real_t& x ) {
    const real_t u = x / sqrt(2.0);
    const real_t e = ( u < 0.0 ) ? -erff1(-u) : erff1(u);

    return ( 0.5 * (1.0 + e) );
}


real_t erff_poly_only( const real_t& x ) {
    const real_t p = 0.3275911;
    const real_t a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429;

    const real_t t  = 1.0/(1.0+p*x);
    const real_t t2 = t  * t;
    const real_t t3 = t  * t2;
    const real_t t4 = t2 * t2;
    const real_t t5 = t2 * t3;

    const real_t res = (a1*t + a2*t2 + a3*t3 + a4*t4 + a5*t5);
    return res;
    // erff was: ( 1.0 - res * exp(-(x*x)) );
}


real_t uGaussian_P_withExpFactor( const real_t& x, const real_t exp_factor ) {
    const real_t u = fabs( x / sqrt(2.0) );
    const real_t e = erff_poly_only(u);

    real_t res = 0.5 * e * exp(exp_factor-u*u);

    if( x >= 0.0 ) {
        res = exp(exp_factor) - res;
    }

    return res;
}


//-------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////
// THE FUNction-PARAMeter to rootFinding_Brent
// if fid = 33 then: to_solve(real x) = (x+3.0)*(x-1.0)*(x-1.0)
// otherwise follows the real implementation
////////////////////////////////////////////////////////////////////////
// N is the size of scales (bbi)
real_t to_solve( const UINT& fid, const UINT& N, const real_t* scales, const real_t* bbi, const real_t& yhat ) { 
    if ( fid == 33 ) return (yhat+3.0)*(yhat-1.0)*(yhat-1.0);
    else {
        real_t accum = 0.0;
        for( UINT i = 0; i < N; i++ ) {
            accum += scales[i] * exp( - bbi[i] * yhat );
        }
        return accum - 1.0;
    }
}

/**
def bisection(f,lb,ub,N=50):
  for j in xrange(N):
    x = (lb + ub) / 2.0
    if f(lb)*f(x)>0: lb = x
    else: ub = x
    print "bisection: j: %d, x: %.3f, f(x): %.3f" % (j,x,f(x))
  return x
**/
void rootBisection( 
                            const UINT&   fid, 
                            const UINT&   N, 
                            const real_t*   scale, 
                            const real_t*   bbi, 
                            const real_t&   lb, 
                            const real_t&   ub, 
                            const real_t&   toll, 
                            const UINT&   it_mx,
                                  real_t&   root,  // result
                                  UINT&   it,
                                  real_t&   fb
) {
   // IMPLEMENT THIS HERE!!!
    const real_t tol      = (toll  <= 0.0) ?  1.0e-9 : toll;
    const real_t iter_max = (it_mx <= 0  ) ?  IT_MAX : it_mx;

    real_t a = lb, b = ub;

    real_t x, fx;
    it         = 0;

    real_t fa = to_solve(fid, N, scale, bbi, a);
         fb = to_solve(fid, N, scale, bbi, b);

    if( fa*fb >= 0.0 ) {
        root = 0.0;  it   = 0;
        if ( a >= 0.0 ) fb =  INFTY;
        else            fb = -INFTY;
        return;
    } 

    for( UINT i = 0; i < iter_max; i++ ) {
        x=(a+b)/2.0;
        //fa = to_solve(fid, N, scale, bbi, a);
        fx = to_solve(fid, N, scale, bbi, x);
        
        if ( fa*fx > 0.0 ) { a = x; fa = fx; }
        else               { b = x; }
        

        // reporting non-convergence!
        if(i == iter_max-1) {
             printf("# ERROR: Bisection method not converged, error: %f %f %d\n\n", b, fb, i);
        }

        it = i;
        if ( fx == 0.0 || fabs(b-a) <= tol ) { break; }
    }

    root = x;
    fb = fx;
}



void rootFinding_Brent ( 
                            const UINT&   fid, 
                            const UINT&   N, 
                            const real_t*   scale, 
                            const real_t*   bbi, 
                            const real_t&   lb, 
                            const real_t&   ub, 
                            const real_t&   toll, 
                            const UINT&   it_mx,
                                  real_t&   root,  // result
                                  UINT&   it,
                                  real_t&   fb
) { 
    const real_t tol      = (toll  <= 0.0) ?  1.0e-9 : toll;
    const real_t iter_max = (it_mx <= 0  ) ?  IT_MAX : it_mx;

    real_t a = lb, b = ub;

    real_t fa = to_solve(fid, N, scale, bbi, a);
         fb = to_solve(fid, N, scale, bbi, b);

    if( fa*fb >= 0.0 ) {
        root = 0.0;  it   = 0;
        if ( a >= 0.0 ) fb =  INFTY;
        else            fb = -INFTY;
        return;
    } 

    if( fabs(fa) < fabs(fb) ) { real_t tmp = fa; fa = fb; fb = tmp; tmp = a; a = b; b = tmp; }
    
    real_t c = a, fc = fa;
    bool mflag = true;
    real_t d     = 0.0;
    it         = 0;

    for( UINT i = 0; i < iter_max; i++ ) {
        if ( fb != 0.0 && fabs(b-a) >= tol ) {
            real_t s;
            if( fa == fc || fb == fc ) {
                s = b - fb * (b - a) / (fb - fa);
            } else {
                real_t s1 = (a*fb*fc)/( (fa-fb)*(fa-fc) );
                real_t s2 = (b*fa*fc)/( (fb-fa)*(fb-fc) );
                real_t s3 = (c*fa*fb)/( (fc-fa)*(fc-fb) );
                s = s1 + s2 + s3;
            }

            if ( ( (3.0 * a + b) /4.0 > s || s > b)        ||
                 (  mflag  && fabs(b-c)/2.0 <= fabs(s-b) ) ||
                 ( !mflag  && fabs(c-d)/2.0 <= fabs(s-b) ) ||
                 (  mflag  && fabs(b-c)     <= fabs(tol) ) ||
                 ( !mflag  && fabs(c-d)     <= fabs(tol) )    ) {
                mflag = true;
                s     = (a + b) / 2.0;
            } else {
                mflag = false;
            }
    
            real_t fs = to_solve(fid, N, scale, bbi, s);

            // d is assigned for the first time here:
            // it's not used above because mflag is set
            d = c;
            c = b; fc = fb;
            
            if( fa*fs < 0.0 ) { b = s; fb = fs; }
            else              { a = s; fa = fs; }

            if( fabs(fa) < fabs(fb) ) { 
                real_t tmp;
                tmp = a;   a =  b;  b = tmp;
                tmp = fa; fa = fb; fb = tmp;
            }

            // reporting non-convergence!
            if(i == iter_max-1) {
                    printf("# ERROR: Brent method not converged, error: %f %f %d\n\n", b, fb, i);
            }

            it = i;
        }
    }

    root = b;
    
    //return BrentRes(b, it, fb);
}




int test_math() {
    // Rootfinder.brent (-4.) (4./.3.) (fun x -> (x+.3.)*.(x-.1.)**2.) 1e-4 == -3
    printf("# Brent test: %f => ", -3.0);
    real_t scale[1] = { 0.0 }, bbi[1] = { 0.0 };

    real_t root, err; UINT it;
    rootFinding_Brent(33, 1, scale, bbi, -4.0, 4.0/3.0, 0.0, 0, root, it, err); 
    if( equalEps(root, -3.0) ) printf(" SUCCESS!\n\n");
    else                       printf(" %f FAILS!\n\n", root);

    // erff 0. == 0. ;; 100. *. erff (1./.sqrt 2.)
    printf("# Erf test: ");
    if( equalEps( erff1(0.0), 0.0 ) && equalEps( floor( 100.0 * erff1( 1.0 / sqrt(2.0) ) ), 68.0 ) )
         printf(" SUCCESS!\n\n");
    else printf(" FAILS  !\n\n");

    // ugaussian_P 0. ;; ugaussian_P 1. +. ugaussian_P (-1.)
    printf("# Gaussian test: ");
    if( equalEps( uGaussian_P(0.0), 0.5 ) && equalEps( uGaussian_P(-1.0)+uGaussian_P(1.0), 1.0 ) )
         printf(" SUCCESS!\n\n");
    else printf(" FAILS  !\n\n");
}

#endif // ifndef MATH_MODULE

