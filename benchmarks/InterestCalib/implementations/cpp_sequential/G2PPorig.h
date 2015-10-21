#ifndef G2PPORIG
#define G2PPORIG

#include "Constants.h"
#include "Date.h"
#include "MathModule.h"
#include "G2ppUtil.h"

#include "GenAlgUtil.h"

real_t to_solve_orig( const UINT&       N, 
                    const IntermElem* tmp_arrs, 
                    const real_t&       yhat  ) {     
    real_t accum = 0.0;
    for( UINT i = 0; i < N; i++ ) {
        accum += tmp_arrs[i].hat_scale * exp( - tmp_arrs[i].bbi * yhat );
    }
    return accum - 1.0;
}

//////////////////////////
// Root finder
//////////////////////////

void rootFinding_Brent_orig ( 
                            const UINT&   N, 
                            IntermElem*   tmp_arrs, 
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

    real_t fa = to_solve_orig(N, tmp_arrs, a);
         fb = to_solve_orig(N, tmp_arrs, b);

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
    
            real_t fs = to_solve_orig(N, tmp_arrs, s);

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

real_t exactYhatOrig( 
                const UINT& n_schedi,

                const real_t& b,         // scals begins
                const real_t& sigmax,
                const real_t& sigmay,
                const real_t& rhoxy,
                const real_t& rhoxyc,
                const real_t& rhoxycs,
                const real_t& mux,
                const real_t& muy,      // scals ends
                IntermElem* tmp_arrs, 

                const real_t& x   // output
) {
    // ugaussian_Pinv(k)=1.0e~4
    const real_t k = - 3.71901648545568;

    real_t up = 0.0, lo = -INFTY;
    for( UINT i = 0; i < n_schedi; i++ ) {
        real_t baix = tmp_arrs[i].bai * x;

        real_t up_term = tmp_arrs[i].aici * exp( -baix );
        tmp_arrs[i].hat_scale = up_term;
        up += up_term;
        lo  = std::max( lo, ( tmp_arrs[i].log_aici - baix ) / tmp_arrs[i].bbi );
    }

// CHECKING uplo!!!!

    const real_t log_s = log(up);
    real_t tmp   = log_s / tmp_arrs[n_schedi-1].bbi;

    if ( tmp <= 0.0 ) {
        up = tmp;
    } else {
        tmp = log_s / tmp_arrs[0].bbi;
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
        rootFinding_Brent_orig(n_schedi, tmp_arrs, root_lb, root_ub, 1.0e-4, 1000, root, iter, error);
        //rootBisection(1, n_schedi, scales, bbi, root_lb, root_ub, 1.0e-4, 1000, root, iter, error);

        res = ( error == -INFTY ) ?  y0 - 1.0 : ( error ==  INFTY ) ? y1 + 1.0 : root;
    }

    return res;
}


#endif // ifndef G2PP
