// Various Kinds of Regular and Segmented Reduction
#ifndef GEN_ALG_KERNELS
#define GEN_ALG_KERNELS

/**
 * One-dimensional kernel: the global size is Npop. 
 * Always accepts the proposed likelihood. 
 */
__kernel void copyLogLik (
        __global REAL  *genomes,
        uint           Npop
) { 
    uint gid = get_global_id(0);
    genomes += 10 * Npop;
    if ( gid < Npop ) {
        // logLik[i] = logLik[i+POP_SIZE]; 
        genomes[ gid ] = genomes[ gid + Npop ];
    }
}

/**
 * One-dimensional kernel: the global size is Npop. 
 * Computes the acceptance condition.
 * Receives Npop random numbers!
 */
__kernel void accept_cond(
        __global REAL  *genomes,
        __global REAL  *rand_nums,
        __global uchar *to_accept,
                 uint   Npop
){
    uint gid = get_global_id(0);
    genomes += 11 * Npop;
    if ( gid < Npop ) {
        // logLik[i] = logLik[i+POP_SIZE]; 
        REAL acceptance = min( 1.0, exp( genomes[gid] - genomes[gid-Npop] ) * genomes[gid+Npop] );
        to_accept[gid] = ( rand_nums[gid] < acceptance );
    }
}
/**
 * Two-dimensional kernel: outer dimension is 6, inner dim is Npop. 
 * Always accepts the proposed likelihood. 
 */
__kernel void accept_prop(
        __global REAL  *genomes,
        __global uchar *to_accept,
        uint           Npop
){
    uint gid = get_global_id(0);
    genomes += get_global_id(1) * 2 * Npop;
    if( gid < Npop && to_accept[gid] ) {
        genomes[gid] = genomes[gid+Npop];
    }
}

/*****************************************/
/*****************************************/
/*****************************************/

inline REAL2 
mutate_helper( 
        REAL gene, 
        REAL gene_prop, 
        REAL g_max,
        REAL g_min,
        REAL r01,
        REAL amplitude_ratio 
) {
    REAL forward_range, backward_range;
    REAL tmp_min_max,   tmp_max_min;

    REAL amplitude     = fabs( (g_max - g_min) * amplitude_ratio );
    REAL semiamplitude = amplitude / 2.0;

    tmp_min_max = min( g_max, gene + semiamplitude );
    tmp_max_min = max( g_min, gene - semiamplitude );
    forward_range = tmp_min_max - tmp_max_min;

    tmp_min_max = min( g_max, gene_prop + semiamplitude );
    tmp_max_min = max( g_min, gene_prop - semiamplitude );
    backward_range = tmp_min_max - tmp_max_min;

    REAL bf_fact = ( semiamplitude > 0 ) ? (backward_range / forward_range) : 1.0;

    // assign p'
    REAL diff = amplitude * r01 - semiamplitude;
    return (REAL2) (gene + diff, bf_fact);
}

inline 
REAL constrain_gene( REAL gene, REAL lb, REAL ub ) {
    return max( lb, min( ub, gene ) );
}

/**
 * One-Dimensional Kernel:
 *    Dimension is POP_SIZE
 */
__kernel void 
mutate_dims_all(__global   REAL* genomes,    //[13, POP_SIZE]
                __global   REAL* rand_nums,  //[5, POP_SIZE]
                __constant REAL* lub,        //[2, 5]. i.e., [lower,upper] gene range 
                           REAL  amplitude_ratio,  // = MOVES_UNIF_AMPL_RATIO 
                           uint  Npop              // POP_SIZE
) {
    uint gidr = get_global_id(0);
    uint gid  = gidr + Npop, offset = 2 * Npop;
    REAL fb_rat = 1.0;    
    REAL2 tmp;

    if(gidr > Npop) {
        for(int i=0; i < 5; i++) {
            REAL lb = lub[i], ub = lub[i+5];

            tmp = mutate_helper( genomes[gid-Npop], genomes[gid], 
                                 lb, ub, rand_nums[gidr], amplitude_ratio );

            REAL val = constrain_gene(tmp.x, lb, ub);
            genomes[gid] = (i == 2 && val == 0) ? EPS0 : val;

            fb_rat      *= tmp.y;

            gid         += offset;
            gidr        += Npop;
        }
        gid += Npop;
        genomes[gid] = fb_rat;
    }
}

/**
 * One-Dimensional Kernel:
 *    Dimension is POP_SIZE
 */
__kernel void 
mutate_dims_one(__global   REAL* genomes,    //[13, POP_SIZE]
                __global   REAL* rand_nums,  //[5, POP_SIZE]
                __constant REAL* lub,        //[2, 5]. i.e., [lower,upper] gene range
                           uint  dim_j, 
                           REAL  amplitude_ratio,  // MOVES_UNIF_AMPL_RATIO 
                           uint  Npop              // POP_SIZE
) {
    uint gidr = get_global_id(0);
    uint gid  = gidr + Npop, offset = 2 * Npop;
    REAL fb_rat = 1.0;    
    REAL2 tmp;

    if(gidr > Npop) {
        for(int i=0; i < 5; i++) {
            REAL lb = lub[i], ub = lub[i+5];

            tmp = mutate_helper( genomes[gid-Npop], genomes[gid], 
                                 lb, ub, rand_nums[gidr], (dim_j == i) ? amplitude_ratio : 0.0 );

            REAL val = constrain_gene(tmp.x, lb, ub);
            genomes[gid] = (i == 2 && val == 0) ? EPS0 : val;

            fb_rat      *= tmp.y;

            gid         += offset;
            gidr        += Npop;
        }
        gid += Npop;
        genomes[gid] = fb_rat;
    }
}


inline REAL getRandUnifNorm(REAL r) {
    return r;
}
//Returns a (pseudo) random integer in [0, ub)
inline uint getRandIntNorm(uint ub, REAL r01) {
    return (uint) (r01 * ub);
}

inline
REAL perturbation(  REAL  gene,
                    REAL  gene_k,
                    REAL  gene_l,
                    REAL  r01,
                    REAL  lb,
                    REAL  ub, 
                    REAL  gamma1,
                    REAL  amplitude_ratio //= MOVES_UNIF_AMPL_RATIO 
) {
    REAL amplitude     = fabs( (ub - lb) * amplitude_ratio );
    REAL semiamplitude = amplitude / 2.0;
    REAL perturb       = ( amplitude * r01 - semiamplitude );
    return ( gene + perturb + gamma1 * ( gene_k - gene_l ) );
}

/**
 * Crossover
 * Requires a vector of size [ POP_SIZE*8 ] of random uniform numbers.
 */
__kernel void 
mcmc_DE( __global   REAL *genomes,
         __global   REAL *rand_nums,  // [8, POP_SIZE]
         __constant REAL* lub,        //[2, 5]. i.e., [lower,upper] gene range
                    REAL  gamma_avg,  // 2.38 / sqrt(2.0*GENOME_DIM),
                    REAL  ampl_ratio, // 0.1 * MOVES_UNIF_AMPL_RATIO
                    uint  Npop
) {
        uint gidr   = get_global_id(0);
        uint gid    = gidr;
        uint offset = 2 * Npop; 
        
        if( gidr < Npop ) { 
            uint cand_UB = Npop - 1;
            uint k = getRandIntNorm( cand_UB, rand_nums[gidr] ); gidr += Npop;
            if ( k == gid ) {
                k == cand_UB;
                cand_UB -= 1;
            }
            uint l = getRandIntNorm( cand_UB, rand_nums[gidr] ); gidr += Npop;
            if ( l == gid || l == k ) {
                l = cand_UB;
            }            

            // proposal
            //   gamma: integrated out from the adviced Multivariate Gaussian with Gaussian target (Braak, 2006)
            //   random.uniform: is it supposed to be Sobol or such ? Does C++ rand() works?
            // gamma = N.random.uniform( gamma_avg-0.5, gamma_avg+0.5 )
            REAL gamma1 = gamma_avg - 0.5 + getRandUnifNorm( rand_nums[gidr] ); gidr += Npop;

            for( int j = 0; j < 5; j ++ ) {
                REAL lb  = lub[j], ub = lub[j+5];

                REAL val = perturbation  ( genomes[gid], genomes[k], genomes[l], rand_nums[gidr], lb, ub, gamma1, ampl_ratio );
                val      = constrain_gene( val, lb, ub );

                genomes[ gid + Npop ] = ( j == 2 && val == 0 ) ? EPS0 : val;

                gidr += Npop;
                gid  += offset;
                k    += offset;
                l    += offset;
            }

            gid += offset;
            genomes[gid] = 1.0; // TODO: fix
        }
}


#endif   // GEN_ALG_KERNELS

