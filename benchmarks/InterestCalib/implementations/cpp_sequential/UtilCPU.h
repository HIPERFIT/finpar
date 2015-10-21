#ifndef CPU_UTILITIES
#define CPU_UTILITIES

struct IntermElem {
    real_t ci;
    real_t t1_cs;
    real_t bai;
    real_t bbi;
    real_t aici;
    real_t log_aici;
    real_t scale;
    real_t hat_scale;  

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
    real_t*         new_quote_price;

    SeqArrays (const uint m) : MaxSzTmpArrs(m) {
        genomes         = new Genome[ 2 * POP_SIZE ];
        tmp_arrs        = new IntermElem[MaxSzTmpArrs];

        new_quote_price = new real_t [ 2 * NUM_SWAP_QUOTES ];
    }

    void releaseResources() {
        delete[] genomes;
        delete[] tmp_arrs;
        delete[] new_quote_price;
    }

    // get the quotes prices for genome i
    real_t* get_quote () { return new_quote_price;                   }
    real_t* get_price () { return new_quote_price + NUM_SWAP_QUOTES; }
};
#endif // CPU_UTILITIES
