#ifndef G2PPORIG
#define G2PPORIG

#include "Constants.h"
#include "Date.h"
#include "MathModule.h"
#include "G2ppUtil.h"

#include "GenAlgUtil.h"

REAL to_solve_orig( const UINT&       N, 
                    const IntermElem* tmp_arrs, 
                    const REAL&       yhat  ) {     
    REAL accum = 0.0;
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
                            const REAL&   lb, 
                            const REAL&   ub, 
                            const REAL&   toll, 
                            const UINT&   it_mx,
                                  REAL&   root,  // result
                                  UINT&   it,
                                  REAL&   fb
) { 
    const REAL tol      = (toll  <= 0.0) ?  1.0e-9 : toll;
    const REAL iter_max = (it_mx <= 0  ) ?  IT_MAX : it_mx;

    REAL a = lb, b = ub;

    REAL fa = to_solve_orig(N, tmp_arrs, a);
         fb = to_solve_orig(N, tmp_arrs, b);

    if( fa*fb >= 0.0 ) {
        root = 0.0;  it   = 0;
        if ( a >= 0.0 ) fb =  INFTY;
        else            fb = -INFTY;
        return;
    } 

    if( fabs(fa) < fabs(fb) ) { REAL tmp = fa; fa = fb; fb = tmp; tmp = a; a = b; b = tmp; }
    
    REAL c = a, fc = fa;
    bool mflag = true;
    REAL d     = 0.0;
    it         = 0;

    for( UINT i = 0; i < iter_max; i++ ) {
        if ( fb != 0.0 && fabs(b-a) >= tol ) {
            REAL s;
            if( fa == fc || fb == fc ) {
                s = b - fb * (b - a) / (fb - fa);
            } else {
                REAL s1 = (a*fb*fc)/( (fa-fb)*(fa-fc) );
                REAL s2 = (b*fa*fc)/( (fb-fa)*(fb-fc) );
                REAL s3 = (c*fa*fb)/( (fc-fa)*(fc-fb) );
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
    
            REAL fs = to_solve_orig(N, tmp_arrs, s);

            // d is assigned for the first time here:
            // it's not used above because mflag is set
            d = c;
            c = b; fc = fb;
            
            if( fa*fs < 0.0 ) { b = s; fb = fs; }
            else              { a = s; fa = fs; }

            if( fabs(fa) < fabs(fb) ) { 
                REAL tmp;
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

REAL exactYhatOrig( 
                const UINT& n_schedi,

                const REAL& b,         // scals begins
                const REAL& sigmax,
                const REAL& sigmay,
                const REAL& rhoxy,
                const REAL& rhoxyc,
                const REAL& rhoxycs,
                const REAL& mux,
                const REAL& muy,      // scals ends
                IntermElem* tmp_arrs, 

                const REAL& x   // output
) {
    // ugaussian_Pinv(k)=1.0e~4
    const REAL k = - 3.71901648545568;

    REAL up = 0.0, lo = -INFTY;
    for( UINT i = 0; i < n_schedi; i++ ) {
        REAL baix = tmp_arrs[i].bai * x;

        REAL up_term = tmp_arrs[i].aici * exp( -baix );
        tmp_arrs[i].hat_scale = up_term;
        up += up_term;
        lo  = std::max( lo, ( tmp_arrs[i].log_aici - baix ) / tmp_arrs[i].bbi );
    }

// CHECKING uplo!!!!

    const REAL log_s = log(up);
    REAL tmp   = log_s / tmp_arrs[n_schedi-1].bbi;

    if ( tmp <= 0.0 ) {
        up = tmp;
    } else {
        tmp = log_s / tmp_arrs[0].bbi;
        if ( 0.0 <= tmp ) up = tmp;
        else              up = - INFTY;
    }

    const REAL yl = lo - EPS;
    const REAL yu = up + EPS;

    const REAL y0 = sigmay * ( rhoxy * (x-mux) / sigmax + k * rhoxycs ) - rhoxyc/b + muy;
    const REAL y1 = sigmay * ( rhoxy * (x-mux) / sigmax - k * rhoxycs )            + muy;

    REAL res;
    if      ( y1 <= yl ) res = y1 + 1.0;  // yhat is greater than y1 => 1 - phi(h_i(x, yhat)) < EPS
    else if ( yu <= y0 ) res = y0 - 1.0;  // yhat is lower than y0 => phi(h_i(x, yhat)) < EPS)
    else {        
        const REAL root_lb = std::max( yl, y0 );
        const REAL root_ub = std::min( yu, y1 );

        REAL root, error; UINT iter;
        rootFinding_Brent_orig(n_schedi, tmp_arrs, root_lb, root_ub, 1.0e-4, 1000, root, iter, error);
        //rootBisection(1, n_schedi, scales, bbi, root_lb, root_ub, 1.0e-4, 1000, root, iter, error);

        res = ( error == -INFTY ) ?  y0 - 1.0 : ( error ==  INFTY ) ? y1 + 1.0 : root;
    }

    return res;
}


#endif // ifndef G2PP
