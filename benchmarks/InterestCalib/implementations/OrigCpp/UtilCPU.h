#ifndef CPU_UTILITIES
#define CPU_UTILITIES

struct IntermElem {
    REAL ci;
    REAL t1_cs;
    REAL bai;
    REAL bbi;
    REAL aici;
    REAL log_aici;
    REAL scale;
    REAL hat_scale;  

    IntermElem() { }
};

struct SeqArrays {

    // max size of temp arrays
    const uint    MaxSzTmpArrs; 

    // [ MaxSzTmpArrs ]
    IntermElem* tmp_arrs;

    // [ 2 * POP_SIZE ]     
    Genome*       genomes;   

    // [ 2 * NUM_SWAP_QUOTES ]
    REAL*         new_quote_price;

    SeqArrays (const uint m) : MaxSzTmpArrs(m) {
        genomes         = new Genome[ 2 * POP_SIZE ];
        tmp_arrs        = new IntermElem[MaxSzTmpArrs];

        new_quote_price = new REAL [ 2 * NUM_SWAP_QUOTES ];
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
