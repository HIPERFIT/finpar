#ifndef CPU_UTILITIES
#define CPU_UTILITIES

/************************/
/*** Helper Functions ***/
/************************/
bool is_pow2(uint n) {
    uint c = 1;
    while(c < n) {
        c = c << 1;
    }
    return (c == n);
}


struct CpuArrays {
    // size of the shape array
    const uint SS;

    // shape of the irregular arrays for one genome
    short* shape;

    // [13 * POP_SIZE] = 
    // { a, a_p, b, b_p, rho, rho_p, nu, nu_p, 
    //   sigma, sigma_p, logLik, logLik_p, bf_rat }
    REAL* genomes; 

    // [ 4 * SS * POP_SIZE ]
    REAL* ci_t1cs_scale;  

    // [ 2 * NUM_SWAP_QUOTES * POP_SIZE ]
    REAL* new_quote_price;

    // [ NUM_SWAP_QUOTES*NUM_HERMITE ]
    REAL* accum0;

    // [ 8 * NUM_SWAP_QUOTES * POP_SIZE ]
    // { mux, muy, sqrt_sigmax = sqrt(2.0) * sigmax, 
    //   t2 = rhoxy / (sigmax*rhoxycs), sigmay_rhoxycs, zc_mat, f, df } 
    REAL* scalars; 

    REAL  gene_ranges[10];

    CpuArrays (const uint n, short* shp) : SS(n) {
        shape           = shp;
    
        genomes         = new REAL[ 13 * POP_SIZE ];
        ci_t1cs_scale   = new REAL[  4 * SS * POP_SIZE ];
        scalars         = new REAL[  8 * NUM_SWAP_QUOTES * POP_SIZE ];

        new_quote_price = new REAL [ 2 * NUM_SWAP_QUOTES * POP_SIZE ];
        accum0          = new REAL [ NUM_SWAP_QUOTES * NUM_HERMITE * POP_SIZE ];

        for(int i=0; i<5; i++) { gene_ranges[i]   = g_mins[i]; }
        for(int i=0; i<5; i++) { gene_ranges[i+5] = g_maxs[i]; }
    }

    void releaseResources() {
        delete[] shape;
        delete[] genomes;
        delete[] ci_t1cs_scale;
        delete[] scalars;
        delete[] new_quote_price;
        delete[] accum0;
    }

    // genome helpers
    REAL* get_a     () { return genomes; }
    REAL* get_b     () { return genomes + POP_SIZE*2; }
    REAL* get_rho   () { return genomes + POP_SIZE*4; }
    REAL* get_nu    () { return genomes + POP_SIZE*6; }
    REAL* get_sigma () { return genomes + POP_SIZE*8; }
    REAL* get_logLik() { return genomes + POP_SIZE*10;}
    REAL* get_bf_rat() { return genomes + POP_SIZE*12;}

    // get the shape of the irregular array (for one genome)
    short* get_shape (){ return shape; }
    
    // get the start iterator into arrays ci, t1cs, scale for genome i
    REAL* get_ci    (int i) { return ci_t1cs_scale + SS*i;                 }
    REAL* get_t1cs  (int i) { return ci_t1cs_scale + SS*i +   SS*POP_SIZE; }
    REAL* get_scale (int i) { return ci_t1cs_scale + SS*i + 2*SS*POP_SIZE; }
    REAL* get_bbi   (int i) { return ci_t1cs_scale + SS*i + 3*SS*POP_SIZE; }

    // get the start iterator into array-expanded scalars for genome i
    REAL* get_mux   (int i) { return scalars + i * NUM_SWAP_QUOTES ; }
    REAL* get_muy   (int i) { return scalars + i * NUM_SWAP_QUOTES +     NUM_SWAP_QUOTES * POP_SIZE; }
    REAL* get_sqsigx(int i) { return scalars + i * NUM_SWAP_QUOTES + 2 * NUM_SWAP_QUOTES * POP_SIZE; }
    REAL* get_t2    (int i) { return scalars + i * NUM_SWAP_QUOTES + 3 * NUM_SWAP_QUOTES * POP_SIZE; }
    REAL* get_sigrho(int i) { return scalars + i * NUM_SWAP_QUOTES + 4 * NUM_SWAP_QUOTES * POP_SIZE; }
    REAL* get_zcmat (int i) { return scalars + i * NUM_SWAP_QUOTES + 5 * NUM_SWAP_QUOTES * POP_SIZE; }
    REAL* get_f     (int i) { return scalars + i * NUM_SWAP_QUOTES + 6 * NUM_SWAP_QUOTES * POP_SIZE; }
    REAL* get_df    (int i) { return scalars + i * NUM_SWAP_QUOTES + 7 * NUM_SWAP_QUOTES * POP_SIZE; }

    // get the quotes prices for genome i
    REAL* get_quote (int i) { return new_quote_price + i * NUM_SWAP_QUOTES; }
    REAL* get_price (int i) { return new_quote_price + i * NUM_SWAP_QUOTES + NUM_SWAP_QUOTES * POP_SIZE; }

    // get the accumulator for the irregular array for genome i 
    REAL* get_accum (int i) { return accum0 + i * NUM_SWAP_QUOTES * NUM_HERMITE; }
};
#endif // CPU_UTILITIES
