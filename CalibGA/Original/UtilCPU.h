#ifndef CPU_UTILITIES
#define CPU_UTILITIES

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

    // [ MaxSzTmpArrs ]
    IntermArrays* tmp_arrs;

    // [ 2 * POP_SIZE ]     
    Genome*       genomes;   

    // [ 2 * NUM_SWAP_QUOTES ]
    REAL*         new_quote_price;

    REAL          gene_ranges[10];

    CpuArrays (const uint m) : MaxSzTmpArrs(m) {
        genomes         = new Genome[ 2 * POP_SIZE ];
        tmp_arrs        = new IntermArrays[MaxSzTmpArrs];

        new_quote_price = new REAL [ 2 * NUM_SWAP_QUOTES ];

        for(int i=0; i<5; i++) { gene_ranges[i]   = g_mins[i]; }
        for(int i=0; i<5; i++) { gene_ranges[i+5] = g_maxs[i]; }
    }

    void releaseResources() {
        delete[] genomes;
        delete[] tmp_arrs;
        delete[] new_quote_price;
    }

    // get the quotes prices for genome i
    REAL* get_quote () { return new_quote_price;                   }
    REAL* get_price () { return new_quote_price + NUM_SWAP_QUOTES; }
};
#endif // CPU_UTILITIES
