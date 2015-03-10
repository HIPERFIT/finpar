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


struct IntermArrays {
    REAL ci;
    REAL t1_cs;
    REAL bai;
    REAL bbi;
    REAL aici;
    REAL log_aici;
    REAL scale;
    REAL hat_scale;  

    IntermArrays() { }
};

struct CpuArrays {

    // max size of temp arrays
    const uint    MaxSzTmpArrs; 

    // [MaxSzTmpArrs]
    IntermArrays* tmp_arrs;

    // [2*POP_SIZE]     
    Genome* genomes;   

    // [ 2 * NUM_SWAP_QUOTES * POP_SIZE ]
    REAL* new_quote_price;

    REAL  gene_ranges[10];

    CpuArrays (const uint m) : MaxSzTmpArrs(m) {
        genomes         = new Genome[ 2 * POP_SIZE ];
        tmp_arrs        = new IntermArrays[MaxSzTmpArrs];

        new_quote_price = new REAL [ 2 * NUM_SWAP_QUOTES * POP_SIZE ];

        for(int i=0; i<5; i++) { gene_ranges[i]   = g_mins[i]; }
        for(int i=0; i<5; i++) { gene_ranges[i+5] = g_maxs[i]; }
    }

    void releaseResources() {
        delete[] genomes;
        delete[] tmp_arrs;
        delete[] new_quote_price;
    }

    // get the quotes prices for genome i
    REAL* get_quote (int i) { return new_quote_price + i * NUM_SWAP_QUOTES; }
    REAL* get_price (int i) { return new_quote_price + i * NUM_SWAP_QUOTES + NUM_SWAP_QUOTES * POP_SIZE; }
};
#endif // CPU_UTILITIES
