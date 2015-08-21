// Various Kinds of Regular and Segmented Reduction
#ifndef GEN_ALG_KERNELS
#define GEN_ALG_KERNELS

real_t getSobolNum( uint n, __constant int* SOBOL_DIR_VCT ) {
    int n_gray = (n >> 1) ^ n;
    int res = 0;
    for( int i=0; i < SOBOL_BITS_NUM; i++ ) {
        int  t    = (1 << i);
        bool cond = ( (n_gray & t) == t ); 
        if ( cond ) {
            res = res ^ SOBOL_DIR_VCT[i];
        }
    }
    real_t rres = ((real_t)res) / ( (1<<SOBOL_BITS_NUM) + 1.0F);
    return rres;
}

#define SOBOL_TRANSP 1

/***************************************/
/******** Initialization Kernel ********/
/***************************************/
/**
 * Two-dimensional kernel: 
 *   outer dimension is 6, i.e., {a, b, rho, nu, sigma, bf_rat}, 
 *   inner dim is population size, i.e., Npop = POP_SIZE,
 *        and LWG == GWG!  
 */
__kernel void init_genome(
        __global   real_t *genomes,
        __constant real_t *lub,            //[2, 5]. i.e., [lower,upper] gene range 
        __constant int  *SOBOL_DIR_VCT,  //[SOBOL_BITS_NUM == 30]   
                 uint   sob_offset,
                 uint   Npop
){
    uint gid     = get_global_id(0);
    uint gene_num= get_global_id(1);

    if( gid < Npop ) {
        if( gene_num == 5 ) { // bf_rat 
            //bf_rat[i] = 1.0; 
            genomes[12 * Npop + gid] = 1.0F;
        } else {

#if (SOBOL_TRANSP == 1)
            uint sobnum = sob_offset + gid * 5 + gene_num;
#else
            uint sobnum = sob_offset + gene_num * Npop + gid;
#endif
            real_t r01    = getSobolNum( sobnum, SOBOL_DIR_VCT );

            real_t lb = lub[gene_num], ub = lub[gene_num + 5];
            r01 = r01 * ( ub - lb ) + lb;

            gid += ((gene_num * Npop) << 1);
            genomes[gid       ] = r01;
            genomes[gid + Npop] = r01;
        }
    }
}

/**
 * One-dimensional kernel: the global size is Npop. 
 * Always accepts the proposed likelihood. 
 */
__kernel void copyLogLik (
        __global real_t  *genomes,
        uint           Npop
) { 
    uint gid = get_global_id(0);
    genomes += 10 * Npop;
    if ( gid < Npop ) {
        // logLik[i] = logLik[i+POP_SIZE]; 
        genomes[ gid ] = genomes[ gid + Npop ];
    }
}


/*******************************************/
/*********** ACCEPT NEW GENOME   ***********/
/*******************************************/

/**
 * One-dimensional kernel: the global size is Npop = LWG = GWG 
 * Computes the acceptance condition.
 * Receives Npop random numbers!
 */
__kernel void accept_cond(
        __global real_t  *genomes,
        __constant int *SOBOL_DIR_VCT,  //[SOBOL_BITS_NUM == 30]   
        __global uchar *to_accept,
                 uint   sob_offset,
                 uint   Npop
){
    uint gid = get_global_id(0);
    genomes += 11 * Npop;
    if ( gid < Npop ) {
        real_t acceptance = minR( 1.0F, exp( genomes[gid] - genomes[gid-Npop] ) * genomes[gid+Npop] );
        real_t r01        = getSobolNum( gid + sob_offset, SOBOL_DIR_VCT );
        to_accept[gid]  = ( r01 < acceptance ) ? 1 : 0;
    }
}
/**
 * Two-dimensional kernel: outer dimension is 6, 
 * inner dim is Npop and LWG == GWG! 
 * Always accepts the proposed likelihood. 
 */
__kernel void accept_prop(
        __global real_t  *genomes,
        __global uchar *to_accept,
        uint            Npop
){
    uint gid = get_global_id(0);
    genomes += get_global_id(1) * 2 * Npop;
    if( gid < Npop && to_accept[gid] ) {
        genomes[gid] = genomes[gid+Npop];
    }
}


/*******************************************/
/***********    MUTATE HELPERS   ***********/
/*******************************************/

inline real2_t 
mutate_helper( 
        real_t gene, 
        real_t gene_prop, 
        real_t g_min, 
        real_t g_max,
        real_t r01,
        real_t amplitude_ratio 
) {
    real_t forward_range, backward_range;
    real_t tmp_min_max,   tmp_max_min;

    real_t amplitude     = fabs( (g_max - g_min) * amplitude_ratio );
    real_t semiamplitude = amplitude / 2.0F;

    tmp_min_max = minR( g_max, gene + semiamplitude );
    tmp_max_min = maxR( g_min, gene - semiamplitude );
    forward_range = tmp_min_max - tmp_max_min;

    tmp_min_max = minR( g_max, gene_prop + semiamplitude );
    tmp_max_min = maxR( g_min, gene_prop - semiamplitude );
    backward_range = tmp_min_max - tmp_max_min;

    real_t bf_fact = ( semiamplitude > 0.0F ) ? (backward_range / forward_range) : 1.0F;

    // assign p'
    real_t diff = amplitude * r01 - semiamplitude;
    return (real2_t) (gene + diff, bf_fact);
}

inline 
real_t constrain_gene( real_t gene, real_t lb, real_t ub ) {
    return maxR( lb, minR( ub, gene ) );
}

/*******************************************/
/**********  MUTATE DIMs KERNEL  ***********/
/*******************************************/

/**
 * One-Dimensional Kernel:
 *    Dimension is POP_SIZE (LWG = GWG)
 */
__kernel void 
mutate_dims (   __global   real_t* genomes,    //[13, POP_SIZE]
                __constant real_t* lub,        //[2, 5]. i.e., [lower,upper] gene range 
                __constant int * SOBOL_DIR_VCT,  //[SOBOL_BITS_NUM == 30]   
                           uint  sob_offset,
                           uint  dim_j,
                           real_t  amplitude_ratio,  // = MOVES_UNIF_AMPL_RATIO 
                           uint  Npop              // POP_SIZE
) {
    uint gid    = get_global_id(0);
#if (SOBOL_TRANSP == 1)
    uint gidr   = gid*5 + sob_offset;
#else
    uint gidr   = gid + sob_offset;
#endif
    uint offset = ( Npop << 1 );
    real_t fb_rat = 1.0F;
    real2_t tmp;

    if( gid < Npop ) {
        gid += Npop;
        for(int i=0; i < 5; i++) {
            real_t lb  = lub[i], ub = lub[i+5];
            real_t r01 = getSobolNum(gidr, SOBOL_DIR_VCT);

            real_t ampl = (dim_j > 4 || dim_j == i) ? amplitude_ratio : 0.0F;
            tmp = mutate_helper( genomes[gid-Npop], genomes[gid], 
                                 lb, ub, r01, ampl );

            real_t val = constrain_gene(tmp.x, lb, ub);
            genomes[gid] = val; //(i == 2 && val == 0.0) ? EPS0 : val;

            fb_rat      *= tmp.y;
            gid         += offset;
#if (SOBOL_TRANSP == 1)
            gidr        += 1;
#else
            gidr        += Npop;
#endif
        }
        gid += Npop;
        genomes[gid] = fb_rat;
    }
}


/**********************************************/
/************* Crossover MCMC_DE **************/
/**********************************************/

inline real_t getRandUnifNorm(real_t r) {
    return r;
}
//Returns a (pseudo) random integer in [0, ub)
inline uint getRandIntNorm(uint ub, real_t r01) {
    return (uint) (r01 * ub);
}

inline
real_t perturbation(  real_t  gene,
                    real_t  gene_k,
                    real_t  gene_l,
                    real_t  r01,
                    real_t  lb,
                    real_t  ub, 
                    real_t  gamma1,
                    real_t  amplitude_ratio //= MOVES_UNIF_AMPL_RATIO 
) {
    real_t amplitude     = fabs( (ub - lb) * amplitude_ratio );
    real_t semiamplitude = amplitude / 2.0F;
    real_t perturb       = ( amplitude * r01 - semiamplitude );
    return ( gene + perturb + gamma1 * ( gene_k - gene_l ) );
}

/**
 * Crossover
 * One-Dimensional of size POP_SIZE = LWG = GWG
 */
__kernel void 
mcmc_DE( __global   real_t *genomes,
         __constant real_t *lub,            // [2, 5]. i.e., [lower,upper] gene range 
         __constant int  *SOBOL_DIR_VCT,  // [SOBOL_BITS_NUM == 30]   
                    uint  sob_offset,
                    real_t  gamma_avg,      // 2.38 / sqrt(2.0*GENOME_DIM),
                    real_t  ampl_ratio,     // 0.1 * MOVES_UNIF_AMPL_RATIO
                    uint  Npop
) {
        real_t r01;
        uint gid    = get_global_id(0);

#if (SOBOL_TRANSP == 1)
        uint gidr   = gid*8 + sob_offset;
        uint rinc   = 1;
#else
        uint gidr   = gid + sob_offset;
        uint rinc   = Npop;
#endif

        uint offset = ( Npop << 1 ); 
        
        if( gid < Npop ) { 
            uint cand_UB = Npop - 1;
            r01 = getSobolNum(gidr, SOBOL_DIR_VCT );  gidr += rinc;
            uint k = getRandIntNorm( cand_UB, r01 );
            if ( k == gid ) {
                k = cand_UB;
                cand_UB -= 1;
            }

            r01 = getSobolNum(gidr, SOBOL_DIR_VCT );  gidr += rinc;
            uint l = getRandIntNorm( cand_UB, r01 );
            if ( l == gid || l == k ) {
                l = cand_UB;
            }            

            // proposal
            //   gamma: integrated out from the adviced Multivariate Gaussian with Gaussian target (Braak, 2006)
            //   random.uniform: is it supposed to be Sobol or such ? Does C++ rand() works?
            // gamma = N.random.uniform( gamma_avg-0.5, gamma_avg+0.5 )
            r01 = getSobolNum(gidr, SOBOL_DIR_VCT );  gidr += rinc;
            real_t gamma1 = gamma_avg - 0.5F + getRandUnifNorm( r01 );

            for( int j = 0; j < 5; j ++ ) {
                real_t lb  = lub[j], ub = lub[j+5];

                r01 = getSobolNum(gidr, SOBOL_DIR_VCT );
                real_t val = perturbation  ( genomes[gid], genomes[k], genomes[l], r01, lb, ub, gamma1, ampl_ratio );
                val      = constrain_gene( val, lb, ub );

                genomes[ gid + Npop ] = val; //( j == 2 && val == 0.0 ) ? EPS0 : val;

                gidr += rinc;
                gid  += offset;
                k    += offset;
                l    += offset;
            }

            gid += offset;
            genomes[gid] = 1.0F; // TODO: fix
        }
}
 

#endif   // GEN_ALG_KERNELS

