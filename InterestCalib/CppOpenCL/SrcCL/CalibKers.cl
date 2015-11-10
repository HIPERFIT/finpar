// BOP is a generic binary associative operator, 
//     e.g., +, *, AND, OR
//#define BOP(a,b) ( (a)+(b) ) 
 
 
#include "KerConsts.h"

#define TH_ID    (get_local_id(0))     
#define WARP_ID  (TH_ID    >> lgWARP)
#define WARP_FST (WARP_ID  << lgWARP)
#define WARP_LST (WARP_FST + (WARP-1))


#if WITH_FLOAT
    typedef float2        REAL2; 
    typedef float4        REAL4; 
#else
    typedef double2       REAL2; 
    typedef double4       REAL4; 
#endif 

inline REAL minR(REAL a, REAL b) { return ( (a <= b) ? a : b ); }
inline REAL maxR(REAL a, REAL b) { return ( (a <= b) ? b : a ); }

#include "GenAlg.cl" 
#include "BestIndValRedKer.cl" 
  
/******************************************************/
/********* EVAL_GENOME KERNELS and HELPERS ************/
/******************************************************/ 

/** Assumes global shape has the following layout:
 *    1. flags, i.e., 1 if start of segment, -1 if iddle thread, and 0 otherwise
 *    2. iota(n_schedi)
 *    3. the swaption index
 *    4. the size of the current segment
 *  Local shape has the following structure:
 *    1. flags as uchar, where 1 denote the start of a segment,
 *    2. cond  as uchar, where 0 denotes thread is active,
 *    3. index as short in iota(n_schedi),
 *    4. swaption's global index as short
 *    5. size fo the current segment as short
 */ 
inline
bool fillShapeMeta(__global short* shape,
                    __local short* sh_mem,
                            uint   Nshp
) {
    bool iddle = false;
    uint  loc_ind = TH_ID, glb_ind = get_global_id(0); // ignore global dimension 1
    short tmp = shape[ glb_ind ];

    { // fill the boolean part!
        __local uchar* shape_meta = (__local uchar*) sh_mem;

        shape_meta[ loc_ind ] = ( tmp != 0 ) ? (uchar) 1 : (uchar) 0;
        loc_ind += get_local_size(0);  
        shape_meta[ loc_ind ] = (uchar) 0;
    }

    glb_ind += Nshp;
    sh_mem[ loc_ind ] = shape[ glb_ind ];          // setting iota(n_schedi) 

    loc_ind += get_local_size(0);
    glb_ind += Nshp;
    sh_mem[ loc_ind ] = shape[ glb_ind ];          // setting the swaption index

    loc_ind += get_local_size(0);
    glb_ind += Nshp;
    tmp      = shape[ glb_ind ];
    if( tmp < 0 ) { iddle = true; tmp = -tmp; }   
    sh_mem[loc_ind] = tmp;

    barrier(CLK_LOCAL_MEM_FENCE);
    return iddle;
}

inline
short getIotaInd  (__local short* shape) { 
    return shape[ TH_ID + get_local_size(0) ];
} 
inline
short getSwapGlbInd(__local short* shape) {
    return shape[ (get_local_size(0) << 1) + TH_ID ]; 
}
inline 
short getSwapLocInd(__local short* shape) {
    short res = getSwapGlbInd(shape) - shape[ get_local_size(0) << 1 ];
    return res;
}
inline 
short getSgmSize   (__local short* shape ) {
    return shape[ (get_local_size(0) << 1) + get_local_size(0) + TH_ID ];
}

#include "Date.cl"
#include "G2ppUtil.cl"
#include "Reductions.cl"
#include "ExactYhat.cl"


/**
 * Two-dimensional kernel: 
 *    1. the outer dimension is POP_SIZE, 
 *    2. the inner dimension is Nshp, i.e., an approximation  
 *         (with idde threads) of the flattened parallelism  
 *         of the loop of count `NUM_SWAP_QUOTES'.
 */

__kernel __attribute__((reqd_work_group_size(LWG_EG, 1, 1)))     
void eval_genome_main (
    __global short *glb_shape,
    __global REAL  *SwaptionQuotes, 
    __global REAL  *genomes,
    __global REAL  *glb_arrs,          // ci ++ t1_cs ++ scale ++ bbi
    __global REAL  *new_quote_price,   // [POP_SIZE * NUM_SWAP_QUOTES * 2]
    __global REAL  *interm_scalars,    // [POP_SIZE * NUM_SWAP_QUOTES * 8]
                                       // { mux, muy, sqrt_sigmax = sqrt(2.0) * sigmax, 
                                       //   t2 = rhoxy / (sigmax*rhoxycs), sigmay_rhoxycs, zc_mat, f, df } 
    uint           Nswap,
    uint           Npop,
    uint           Nshp
) { 
    REAL4 tmp4;
 
    __local  short  shape_meta [ (LWG_EG << 2) ];
    __local  REAL4  sh_mem4    [ LWG_EG + (LWG_EG >> 1) ];

    bool iddle = fillShapeMeta( glb_shape, shape_meta, Nshp );

    REAL swap_freq = 1.0, maturity = 1.0, tmat0 = 1.0; uint n_schedi;
    if( !iddle ) {
        int ttt   = getSwapGlbInd(shape_meta);
        swap_freq = SwaptionQuotes[ (ttt<<2) + 1 ];
        maturity  = add_years( TODAY, SwaptionQuotes[ (ttt<<2) ] );
        n_schedi  = (uint) (12.0 * SwaptionQuotes[ (ttt<<2) + 2 ] / swap_freq);
        tmat0     = date_act_365( maturity, TODAY );
    }

    REAL strike, beg_date, end_date; 
    { // BLACK PRICE COMPUTATION.
        __local uchar* flags   = (__local uchar*) shape_meta;

        tmp4 = (REAL4) (0.0, MAX_DATE, MIN_DATE, 1.0);
        if(!iddle) {
            tmp4.w  = (REAL) getIotaInd(shape_meta);            // i
            tmp4.y  = add_months( maturity, swap_freq*tmp4.w ); // t0
            tmp4.z  = add_months( tmp4.y,   swap_freq        ); // tn
            tmp4.x  = zc(tmp4.z) * date_act_365(tmp4.z,tmp4.y); // lvl
        }
        beg_date = tmp4.y; end_date = tmp4.z;

        sh_mem4[ TH_ID ] = tmp4;
        barrier(CLK_LOCAL_MEM_FENCE);
        segm_reduce_plusminmax( sh_mem4, flags ); // ToDo 2

        int last_ind  = TH_ID - getIotaInd( shape_meta );
        if(!iddle) last_ind += n_schedi - 1;
        tmp4 = sh_mem4[last_ind];      // the reduction result!
        barrier(CLK_LOCAL_MEM_FENCE); 
        
        strike = ( zc(tmp4.y) - zc(tmp4.z) ) / tmp4.x;
        if( flags[TH_ID] && (!iddle) ) {
            // if thread is the start of a segment
            //     write new_quote to global memory!
            int ttt   = getSwapGlbInd(shape_meta);
            tmp4.y = 0.5 * SwaptionQuotes[ (ttt<<2) + 3 ] * tmat0;
            // was get_group_id(1)
            new_quote_price[ get_global_id(1)*Nswap + ttt ] = 
                tmp4.x * strike * ( uGaussian_P(tmp4.y) - uGaussian_P(-tmp4.y) );
        }
    } // END BLACK PRICE.

    { // PRICER OF SWAPTION COMPUTATION

        REAL y0, y1, eps, eps_contrib, mux, t4, t1_cs_first;
        {
            REAL a, b, rho, nu, sigma;
            {
                uint tmp_ind = get_global_id(1) + Npop;
                uint Npop2   = (Npop << 1);
                a      = genomes[ tmp_ind ]; tmp_ind += Npop2;
                b      = genomes[ tmp_ind ]; tmp_ind += Npop2;
                rho    = genomes[ tmp_ind ]; tmp_ind += Npop2;
                nu     = genomes[ tmp_ind ]; tmp_ind += Npop2;
                sigma  = genomes[ tmp_ind ]; tmp_ind += Npop2;
            }
//
            { 
                REAL v0_mat, v0_end, vt_end, baii, bbii, tmp;
                tmp4   = bigv( a, b, rho, nu, sigma, tmat0 );
                v0_mat = tmp4.x;
//
                tmp    = date_act_365(end_date, TODAY);
                tmp4   = bigv( a, b, rho, nu, sigma, tmp );
                v0_end = tmp4.x;
//
                tmp    = date_act_365(end_date, maturity);
                tmp4   = bigv( a, b, rho, nu, sigma, tmp );
                vt_end = tmp4.x; 
                //baii   = tmp4.y;
                //bbii   = tmp4.z;
                tmp4.w = 0.5 * ( vt_end - v0_end + v0_mat );  // expo_aici
                tmp4.x = date_act_365( end_date, beg_date ) * strike; // res
            }
//
            {  
                //bool is_first = ((__local uchar*)shape_meta)[TH_ID] && (!iddle);
                bool is_first = (!iddle) && ( (TH_ID == get_local_size(0)-1) || (((__local uchar*)shape_meta)[TH_ID+1]==1) ); // actually is_last!
                REAL  muy, rhoxy, sigmax, sigmay, rhoxyc, rhoxycs, zc_mat, tmp1, tmp2;
                uint gind = get_global_id(1)*Nswap + getSwapGlbInd(shape_meta), offset = Npop*Nswap;   // was get_group_id(1)
//
                mux    = - bigmx( a, b, rho, nu, sigma, TODAY, maturity, TODAY, maturity );
                if(is_first) interm_scalars[ gind ] = mux; gind += offset;
//
                muy    = - bigmy( a, b, rho, nu, sigma, TODAY, maturity, TODAY, maturity );
                if(is_first) interm_scalars[ gind ] = muy; gind += offset;
//
                tmp1   =  sqrt( b_fun(2.0*a, tmat0) );   // sqrt_bfun_a
                sigmax =  sigma * tmp1;
                if(is_first) interm_scalars[ gind ] = sqrt(2.0) * sigmax; gind += offset;
//
                tmp2   =  sqrt( b_fun(2.0*b, tmat0) );   // sqrt_bfun_b
                sigmay =  nu    * tmp2;
                rhoxy  =  rho * b_fun(a+b, tmat0) / (tmp1 * tmp2);                
                t4     = (rhoxy * sigmay) / sigmax;
//
                rhoxyc  = 1.0 - rhoxy * rhoxy;  // used in reduction kernel
                rhoxycs = sqrt( rhoxyc );       // used in reduction kernel
                if(is_first) interm_scalars[ gind ] = rhoxy / (sigmax*rhoxycs); gind += offset;
                if(is_first) interm_scalars[ gind ] = sigmay * rhoxycs;         gind += offset;
//
                zc_mat  = zc(maturity);
                if(is_first) tmp4.x += 1.0;     //cii
                tmp4.x *= zc(end_date) / zc_mat; // fact_aici
                if(is_first) interm_scalars[ gind ] = zc_mat;                   gind += offset;                
//
                t1_cs_first = tmp4.z * (mux * t4 - (muy - 0.5*rhoxyc*sigmay*sigmay*tmp4.z) );
                eps         = 0.5 * sigmax;
                eps_contrib = sigmay * (rhoxy * eps / sigmax);
                y0          = sigmay * KKK * rhoxycs;
                y1          = muy - y0;
                y0         += muy - ( rhoxyc / b );
            }
        } 

        {
            { // save to global arrays, ci, t1_cs, scale, bbi, and 
              // store (aici, baii, bbii, log_aici) in tmp4!
                uint gind = get_global_id(1) * Nshp + get_global_id(0), offset = Nshp * Npop;       
                REAL fact_aici = tmp4.x;
                 
                if( !iddle ) glb_arrs[gind] = fact_aici;                  gind += offset; // ci   [i] 
                if( !iddle ) glb_arrs[gind] = t1_cs_first + tmp4.w;       gind += offset; // t1_cs[i]
                if( !iddle ) glb_arrs[gind] = - ( tmp4.y + tmp4.z * t4 ); gind += offset; // scale[i]
                if( !iddle ) glb_arrs[gind] = tmp4.z;                                     // bbi  [i]

                tmp4.x *= exp(tmp4.w);    // aici 
                tmp4.w += log(fact_aici); // log_aici    
 
                sh_mem4[ TH_ID ] = tmp4;  
                barrier(CLK_LOCAL_MEM_FENCE);    
            } 

            { // finally call `exactYhat', and store `f' and `df' to global memory.        
                REAL f, g, h, df; 
                uint gind = getSwapGlbInd(shape_meta);   
                __local uchar* flags = (__local uchar*)shape_meta; 
                //uint n_schedi = (iddle) ? 1 : (uint)(12.0 * SwaptionQuotes[ (gind<<2) + 2 ] / swap_freq);    

                f = exactYhat( iddle, flags, n_schedi, y0, y1, sh_mem4, mux ); 

                if( flags[TH_ID] && (!iddle) ) { // start of segment
                    uint offset = Npop*Nswap;   
                    gind += get_global_id(1)*Nswap + 6*offset;   // was get_group_id(1)          
                    interm_scalars[ gind ] = f;        
                    gind += offset;           
                } 

                y0 += eps_contrib;        
                y1 += eps_contrib;   
                mux+= eps;
                g = exactYhat( iddle, flags, n_schedi, y0, y1, sh_mem4, mux );

                y0 -= 2.0*eps_contrib; 
                y1 -= 2.0*eps_contrib;
                mux-= 2.0*eps;
                h = exactYhat( iddle, flags, n_schedi, y0, y1, sh_mem4, mux );

                if( flags[TH_ID] && (!iddle) ) { // start of segment
                    interm_scalars[ gind ] = 0.5 * ( g - h ) / eps; // i.e., df
                }  
            }
        }
    }
}

/********************************************************/
/*************** Main Reduction kernel 1 ****************/
/********************************************************/
/**
 * Two-dimensional kernel: 
 *    1. the outermost dimension is POP_SIZE,
 *    2. the middle    dimension is Gauss_DIM
 *    2. the inner dimension is Nshp, i.e., an approximation
 *         (with idde threads) of the flattened parallelism
 *         of the loop of count `NUM_SWAP_QUOTES'.
 */

__kernel __attribute__((reqd_work_group_size(LWG_EG, 1, 1)))
void eval_genome_red1 (
    __global   short *glb_shape,
    __constant REAL  *gauss_coefs,    // [Gauss_DIM]
    __constant REAL  *gauss_weights,  // [Gauss_DIM]

    __global REAL  *glb_arrs,          // ci ++ t1_cs ++ scale ++ bbi
    __global REAL  *interm_scalars,    // [POP_SIZE * NUM_SWAP_QUOTES * 8]
                                       // { mux, muy, sqrt_sigmax = sqrt(2.0) * sigmax, 
                                       //   t2 = rhoxy / (sigmax*rhoxycs), sigmay_rhoxycs, zc_mat, f, df } 

    __global REAL  *accum0,            // [ POP_SIZE * NUM_SWAP_QUOTES * Gauss_DIM ]

    uint           Nswap,
    uint           Nshp
) { 
    __local  short  shape_meta [ LWG_EG << 2 ];
    __local  REAL   sh_mem     [ LWG_EG ];

    bool iddle;
    {
        iddle = fillShapeMeta( glb_shape, shape_meta, Nshp ); 
    }

    REAL x_quad = gauss_coefs  [get_global_id(1)];

    REAL sigmay_rhoxycs, x, h1;
    if(!iddle){
        REAL mux, muy, sqrt2sigmax, t2, f, df, yhat_x;  // Npop == POP_SIZE = get_global_size(2)
        uint gind = get_global_id(2) * Nswap + getSwapGlbInd(shape_meta), offset = get_global_size(2) * Nswap;

        mux            = interm_scalars[ gind ]; gind += offset;
        muy            = interm_scalars[ gind ]; gind += offset;

        // x = sqrt2sigmax * x_quad + mux;
        sqrt2sigmax    = interm_scalars[ gind ]; gind += offset;   // 
        x  = sqrt2sigmax * x_quad + mux;

        t2             = interm_scalars[ gind ]; gind += offset;

        sigmay_rhoxycs = interm_scalars[ gind ]; gind += offset + offset;
        f              = interm_scalars[ gind ]; gind += offset;
        df             = interm_scalars[ gind ];

        yhat_x = f + df*(x - mux);
        h1     = ( (yhat_x - muy) / sigmay_rhoxycs ) - t2*( x - mux );
    }

    if(!iddle){
        uint gind  = get_global_id(2) * Nshp + get_global_id(0), offset = Nshp * get_global_size(2);
             gind += 3 * offset;   
        REAL bbi   = glb_arrs[gind]; gind -= offset;
        REAL h2    = h1 + bbi * sigmay_rhoxycs;
        
        REAL scale = glb_arrs[gind]; gind -= offset;
        REAL t1_cs = glb_arrs[gind]; gind -= offset;
        REAL expo_aici = t1_cs + scale*x;

        REAL expo_part = uGaussian_P_withExpFactor( -h2, expo_aici );
        
        x = glb_arrs[gind] * expo_part;
    }

    { // segmented reduce with operator plus the x's

        __local uchar* flags = (__local uchar*) shape_meta;
        bool is_last = (!iddle) && ( (TH_ID == get_local_size(0)-1) || ( flags[TH_ID+1] == 1) );
        
        sh_mem[TH_ID] = (iddle) ?  0.0 : x; 
        barrier(CLK_LOCAL_MEM_FENCE);

        segm_reduce_plus( sh_mem, flags );

        if( is_last ) {
            REAL accum = sh_mem[TH_ID];
            REAL tmp = sqrt(2.0) * x_quad; //(x - mux) / sigmax;
                 tmp = exp( - 0.5 * tmp * tmp );

            REAL w_quad = gauss_weights[get_global_id(1)]; 

            REAL res    = w_quad * tmp * ( uGaussian_P(-h1) - accum );

            //uint GaussDIM = get_global_size(1);
#if 0
            uint gind = get_global_id(2) * get_global_size(1) * Nswap  +  
                        getSwapGlbInd(shape_meta) * get_global_size(1) + get_global_id(1) ;
            accum0[gind] = res; 
#else
            // accum0 is put it in transposed form, i.e.,
            // accum0[ j * Nswap + ttt ] and is used in 
            // this way by the next kernel
            uint gind = get_global_id(2) * get_global_size(1) * Nswap  +  
                        get_global_id(1) * Nswap + getSwapGlbInd(shape_meta);
            accum0[gind] = res;             
#endif
        }
    }
}



/********************************************************/
/*************** Main Reduction kernel 2 ****************/
/********************************************************/
inline REAL cauchy_pdf( REAL z, REAL mu, REAL gamma ) { 
    REAL x = (z-mu) / gamma;
    return 1.0 / ( PI * gamma * (1.0+x*x) );
}
inline REAL logLikelihood_cauchy( REAL y_ref, REAL y ) {
    REAL LLHOOD_CAUCHY_OFFS = 5.0;
    REAL gamma = ( fabs(y_ref) / 50.0 ) * LLHOOD_CAUCHY_OFFS + 0.01;
    REAL pdfs  = cauchy_pdf( y, y_ref, gamma );
    pdfs += 1.0e-20; // avoid NaNs
    return log(pdfs);
}

/**
 * Two-dimensional kernel: 
 *    1. the outermost global dimension is POP_SIZE,
 *    2. the innermost global dimension is equal to the local 
 *         dimension, i.e., LWG_EG, which is an overapproximation 
 *         of NUM_SWAP_QUOTES.
 * Each local groups process exactly NUM_SWAP_QUOTES elements,
 *   and the rest of the htreads are iddle.
 * This is taking a shortcut to implement reduce in a regular fashion.
 */

__kernel __attribute__((reqd_work_group_size(LWG_EG, 1, 1)))
void eval_genome_red2 (
    __global REAL  *genomes,
    __global REAL  *interm_scalars,    // [POP_SIZE * NUM_SWAP_QUOTES * 8]
                                       // { mux, muy, sqrt_sigmax = sqrt(2.0) * sigmax, 
                                       //   t2 = rhoxy / (sigmax*rhoxycs), sigmay_rhoxycs, zc_mat, f, df } 
    __global REAL  *accum0,            // [ POP_SIZE * NUM_SWAP_QUOTES * Gauss_DIM ]
    __global REAL  *new_quote_price,   // [POP_SIZE * NUM_SWAP_QUOTES * 2]

    uint           Nswap,
    uint           GaussDIM
) { 
    __local  REAL   sh_mem     [ LWG_EG ];

    bool iddle = ( TH_ID >= Nswap );     

    REAL lik = 0.0;
    if (!iddle) { 
        REAL accum = 0.0;
        uint gid = get_global_id(1) * GaussDIM * Nswap + get_global_id(0);
        for( uint j = 0; j < GaussDIM; j++, gid += Nswap ) {
            accum += accum0[ gid ]; 
        } 

        gid = get_global_id(1) * Nswap + TH_ID;  
        REAL zc_mat = interm_scalars[ gid + (5 * get_global_size(1) * Nswap) ];

        REAL nprice = zc_mat * ( accum / sqrt( PI ) );
        new_quote_price[ gid + get_global_size(1)*Nswap ] = nprice;

        lik = logLikelihood_cauchy( new_quote_price[gid], nprice );
    }

    sh_mem[ TH_ID ] = lik;
    barrier(CLK_LOCAL_MEM_FENCE);

    REAL rms = reduce_reg_plus( sh_mem );

    if(TH_ID == 0) genomes[ 11*get_global_size(1) + get_global_id(1) ] = rms;
}

