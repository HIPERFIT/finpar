#ifndef SOBOL_GAUSS
#define SOBOL_GAUSS

/**
 * Computes the x's Sobol quasirandom vector 
 *      of size `sob_dim' in result `res'
 *      using the Sobol-independent formula, 
 *      which does not depend on the previous 
 *      Sobol vector.
 */
inline static
void sobolInd(  const int&  sobol_dim,  // size of the quasi-random vector
                const int&  bit_count,  // number of bits of Sobol numbers
                const int*  dir_matrix, // [bit_count, sobol_dim]: Sobol's
                                        //    direction vectors
                const UINT& x,          // x'th sobol number
                UINT*       res         // result (quasi-random vector)
) {
    UINT gs, y = x + 1;
    gs = (y >> 1);
    gs = (y ^ gs);

    for( int i = 0; i < sobol_dim; i ++ ) {
        res[i] = 0;
    }
    for( int k = 0; k < bit_count; k ++ ) {
        const int* dir_row = dir_matrix + k*sobol_dim;
		if( gs & 1 ) {
            for( int i = 0; i < sobol_dim; i++ ) {
                res[i] ^= dir_row[i];
            }
        }
        gs = gs >> 1;
    }
}

/**
 * Computes the x's Sobol quasirandom vector 
 *      of size `sob_dim' in result `res'
 *      using the Sobol-reccurent formula, 
 *      which depends on the previous Sobol
 *      vector (stored in `res').
 */
inline static
void sobolRec(  const int&  sobol_dim,  // size of the quasi-random vector
                const int&  bit_count,  // number of bits of Sobol numbers
                const int*  dir_matrix, // [bit_count, sobol_dim]: Sobol's
                                        //    direction vectors
                const UINT& x,          // x'th sobol number
                UINT*       res         // in-out param-result    
) {
    UINT y = 0;
    UINT c = x;
	while(c & 1) { 
        y++;
        c >>= 1;
    }
	//assert( (y < bit_count) && "In sobolRec ell > bit_count!!!");

	for( int i = 0; i < sobol_dim; i ++ ) {
        res[i] ^= dir_matrix[y*sobol_dim + i];
		//md_zd[i] = res[i] * ro_scal->sobol_last_den_inv;
	}
}

inline static void
sobolRecOpt(    const int&  sobol_dim,  // size of the quasi-random vector
                const int&  bit_count,  // number of bits of Sobol numbers
                const int*  dir_matrix, // [bit_count, sobol_dim]: Sobol's
                                        //    direction vectors
                const UINT& f_ind,      // lookup-table index
                UINT*       res         // in-out param-result    
) {
    const UINT f_ind_dim = f_ind * sobol_dim;
    for( int i = 0; i < sobol_dim; i ++ )
	    res[i] ^= dir_matrix[ f_ind_dim + i];
}


/******************************************************/
/*** Normal-to-Gaussian Distribution Transformation ***/
/******************************************************/

#define rat_eval(a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7) \
    (x*(x*(x*(x*(x*(x*(x*a7+a6)+a5)+a4)+a3)+a2)+a1)+a0)/ \
    (x*(x*(x*(x*(x*(x*(x*b7+b6)+b5)+b4)+b3)+b2)+b1)+b0)

inline static REAL small_case(const REAL& q) {
  REAL x = 0.180625 - q * q;
  return q * rat_eval(
                      3.387132872796366608,
                      133.14166789178437745,
                      1971.5909503065514427,
                      13731.693765509461125,
                      45921.953931549871457,
                      67265.770927008700853,
                      33430.575583588128105,
                      2509.0809287301226727,

                      1.0,
                      42.313330701600911252,
                      687.1870074920579083,
                      5394.1960214247511077,
                      21213.794301586595867,
                      39307.89580009271061,
                      28729.085735721942674,
                      5226.495278852854561);
}

inline static REAL intermediate(const REAL& r) {
  REAL x = r - 1.6;
  return rat_eval(
                  1.42343711074968357734,
                  4.6303378461565452959,
                  5.7694972214606914055,
                  3.64784832476320460504,
                  1.27045825245236838258,
                  0.24178072517745061177,
                  0.0227238449892691845833,
                  7.7454501427834140764e-4,

                  1.0,
                  2.05319162663775882187,
                  1.6763848301838038494,
                  0.68976733498510000455,
                  0.14810397642748007459,
                  0.0151986665636164571966,
                  5.475938084995344946e-4,
                  1.05075007164441684324e-9);
}

inline static REAL tail(const REAL& r) {
  REAL x = r - 5.0;
  return rat_eval(
                  6.6579046435011037772,
                  5.4637849111641143699,
                  1.7848265399172913358,
                  0.29656057182850489123,
                  0.026532189526576123093,
                  0.0012426609473880784386,
                  2.71155556874348757815e-5,
                  2.01033439929228813265e-7,

                  1.0,
                  0.59983220655588793769,
                  0.13692988092273580531,
                  0.0148753612908506148525,
                  7.868691311456132591e-4,
                  1.8463183175100546818e-5,
                  1.4215117583164458887e-7,
                  2.04426310338993978564e-15);
}

/**
 * Transforms the Sobol quasi-random (int) vector `sov_vct'
 * to a gaussian distribution, uniformely distributed in [-inf, +inf].
 */ 
inline 
void uGaussian(  const REAL& fract,    // 2 ^ sobol_bits
                 const UINT& dim,      // size of the Sobol quasi-random vector
                 const UINT* sob_vct,  // [dim] Sobol quasi-random (int) vector
                       REAL* gauss_vct // result
) {
    for ( int i = 0; i < dim; i ++ ) {
        // sobol int number -> sobol real number in [0,1) 
        REAL sob_real = static_cast<REAL>(sob_vct[i]) * fract;

        // sobol real number [0,1) -> gaussian number in (-inf, +inf)
        REAL dp  = sob_real - 0.5;
        if ( fabs(dp) <= 0.425 ) {
            gauss_vct[i]= small_case(dp);
        } else {
            REAL pp     = (dp < 0.0) ? (0.5 + dp) : (0.5 - dp);
            REAL r      = sqrt (- log(pp));
            REAL x      = (r <= 5.0) ? intermediate(r) : tail(r);
            gauss_vct[i]= (dp < 0.0) ? (0.0 - x) : x;
        }
    }
}

/******************************************************/
/*** Brownian Bridge: Correlate each path's dates   ***/
/******************************************************/

/** 
 * Brownian Bridge entry point.
 */
inline 
void brownianBridge(
        const  UINT&    num_under,
        const  UINT&    num_dates,
        const  int*     bb_inds, // [3, num_dates] Brownian Bridge's indirect indexing
        const  REAL*    bb_data, // [3, num_dates] Brownian Bridge's data
               REAL*    md_zd,   // [num_dates, num_under] the gaussian vector
                                 //                        also holds the result!
               REAL*    res      // [num_dates, num_under] temporary array
) {
    const int * bb_bi = bb_inds;
    const int * bb_li = bb_inds + num_dates;
    const int * bb_ri = bb_inds + 2*num_dates;
    const REAL* bb_sd = bb_data;
    const REAL* bb_lw = bb_data + num_dates;
    const REAL* bb_rw = bb_data + 2*num_dates;

    for ( int m = 0; m < num_under; m ++ ) {
        res[ (bb_bi[0]-1) * num_under + m ] = bb_sd[0] * md_zd[m];
        
        for( int i = 1; i < num_dates; i ++) {
            int j = bb_li[i] - 1;
            int k = bb_ri[i] - 1;
            int l = bb_bi[i] - 1;

            double wk = res  [k*num_under + m];
            double zi = md_zd[i*num_under + m];

            res[l*num_under + m] = (j == -1) ?
                    bb_rw[i] * wk + bb_sd[i] * zi :
                    bb_rw[i] * wk + bb_sd[i] * zi + bb_lw[i] * res[j*num_under + m];
		}
    }

    const UINT S = num_under * num_dates;
#ifdef FAST_BB
	for ( int i = S - 1; i >= num_under; i -- ) {
        res[i] -= res[i - num_under];
	}
#else
    for( int i = 0; i < S; i ++ ) {
        if (i < num_under) 
                md_zd[i] = res[i];
        else    md_zd[i] = res[i] - res[i - num_under];
    }
#endif


}

#endif //SOBOL_GAUSS

