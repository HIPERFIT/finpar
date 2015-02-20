/*************************/
/**** types & tunning ****/
/*************************/

//#define WARP   32
//#define lgWARP 5

// max number of elements to be processed 
// locally by a thread, e.g., for scan.
#define SEQ_ELMS 4




// BOP is a generic binary associative operator, 
//     e.g., +, *, AND, OR
#define BOP(a,b) ((a)+(b))

#define TH_ID    (get_local_id(0))
#define WARP_ID  (TH_ID    >> lgWARP)
#define WARP_FST (WARP_ID  << lgWARP)
#define WARP_LST (WARP_FST + (WARP-1))

/*********************************/
/******** CAS-like Update ********/
/*********************************/

inline 
void raw_cas_update(__local volatile TYPE* sh_hist, const ulong ind, const TYPE val) {
    TYPE acc;
    TYPE id  = (TYPE)get_local_id(0);
    if (val == 0.0) return;
    bool repeat = true;
    while(repeat) {
        acc = sh_hist[ind];
        sh_hist[ind] = id; 
        if( sh_hist[ ind ] == id ) {
            sh_hist[ ind ] = acc + val;
            repeat = false;
        }
    }
}

/**
 * With while is INCORRECT!!!!
 * With for is correct, but it is prohibitively expensive!   
 */
inline 
void raw_cas_update_bar(__local TYPE* sh_hist, const ulong ind, const TYPE val) {
    TYPE acc;
    TYPE id  = (TYPE)get_local_id(0);
    if (val == 0.0) return;
    bool repeat = true;
    //while(repeat) {
    for(uint i=0; i<get_local_size(0); i++) {
        barrier(CLK_LOCAL_MEM_FENCE);
        if(repeat) acc = sh_hist[ind];
        barrier(CLK_LOCAL_MEM_FENCE);
        if(repeat) sh_hist[ind] = id; 
        barrier(CLK_LOCAL_MEM_FENCE);
        if(repeat) {
            if( sh_hist[ ind ] == id ) {
                sh_hist[ ind ] = acc + val;
                repeat = false;
            }
        }
    }
}

inline 
void raw_cas_update_glb(
            __global volatile TYPE*  hist, 
            const             ulong  ind, 
            const             TYPE   val, 
            __local volatile  uchar* sh_flags, 
            const             ulong  M
) {
    TYPE id   = (TYPE)get_local_id(0);
    uint pos  = ind % M; 
    bool repeat = (val != 0.0);
    while(repeat) {
        sh_flags[ pos ] = id;
        if( sh_flags[ pos ] == id ) {
            hist[ ind ] += val;
            repeat = false;
        }
    }
}

/********************************/
/******* Segmented Reduce *******/
/********************************/

/**
 * Inclusive segmeneted scan.
 * ASSERTS: 
 *    1. `size' is a power of two smaller than LWG.
 *    2. `sh_data' pointer is NOT adjusted per warp.
 */
inline
TYPE segm_scan_reg_warp(
        __local volatile TYPE* sh_data,
                         uint  size
) {
    const uint th_id = get_local_id(0) & (WARP-1);

    if( th_id >= 1  && size >= 2  ) sh_data[th_id] = BOP( sh_data[th_id], sh_data[th_id-1 ] );
    if( th_id >= 2  && size >= 4  ) sh_data[th_id] = BOP( sh_data[th_id], sh_data[th_id-2 ] );
    if( th_id >= 4  && size >= 8  ) sh_data[th_id] = BOP( sh_data[th_id], sh_data[th_id-4 ] );
    if( th_id >= 8  && size >= 16 ) sh_data[th_id] = BOP( sh_data[th_id], sh_data[th_id-8 ] );
#if WARP == 32
    if( th_id >= 16 && size >= 32 ) sh_data[th_id] = BOP( sh_data[th_id], sh_data[th_id-16] );
#endif
    return sh_data[th_id];
}

/*
 * `sh_size'   is the size of the shared memory to be scanned.
 * `sgm_size' is the size of the segment. 
 * ASSERTS that:
 *    1. `sgm_size' is a power of two, and
 *    2. `sgm_size' divides the size of the
 *       local group, in particular is less or eq. (<= LWG), and
 *    3. the size of the local group divides `sh_size'. 
**/
inline
void segm_scan_reg_block (
            __local volatile REAL* sh_data,
                             uint  sh_size,
                             uint  sgm_size
) {
    TYPE res;

    for( uint i=TH_ID; i < sh_size; i += get_local_size(0) )
    {
        // perform scan at warp level
        res = segm_scan_reg_warp(sh_data + WARP_FST, sgm_size);
        barrier(CLK_LOCAL_MEM_FENCE);

        if( sgm_size > WARP ) {
            // if last thread in a warp, record it 
            //    at the beginning of sh_data
            if ( TH_ID == WARP_LST ) {
                sh_data[ WARP_ID ] = res; 
            }
            barrier(CLK_LOCAL_MEM_FENCE);

            // first warp scans the per warp results (again)
            //sgm_size = sgm_size >> lgWARP;
            if( TH_ID < WARP ) {   
                segm_scan_reg_warp( sh_data + WARP_FST, (sgm_size >> lgWARP) );
            }

            barrier(CLK_LOCAL_MEM_FENCE);

            // accumulate results from previous step IF `sgm_size'
            // does NOT divide the current thread group number
            if ( WARP_ID % (sgm_size >> lgWARP) ) {
                res = BOP( res, sh_data[ WARP_ID - 1 ] );
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            
            // finally, put the data back in shared memory!
            sh_data[ TH_ID ] = res;

            barrier(CLK_LOCAL_MEM_FENCE);
        }
        
        sh_data += get_local_size(0);
    }
}

/**
 * Accumulates `K' elements locally, then reduces the segments. 
 * In addition to the ASSERTS of segm_scan_reg_block,
 * `K' must divide `sgm_size' (and be a power of two greater than 1!)
 * ASSUMES: K <= sgm_size!
 */
inline
void segm_red_reg_block_K (
            __local volatile TYPE* sh_data,
                             uint  sh_size,
                             uint  sgm_size,
                             char  K
) {

    if (K > 1) {
        // accumulate locally `K' consecutive values,
        //    and update the local memory accordingly
        const uint KKK = K * get_local_size(0); 
        uint i   = K * TH_ID;
        for( ; i < sh_size; i += KKK ) {        
            TYPE res = sh_data[i];
            for( uint j = 1; j < K; j ++ ) {
                res = BOP( res, sh_data[i + j] );
            }
    
            barrier(CLK_LOCAL_MEM_FENCE);
    
            sh_data[i / K] = res;
        }
    }

    { // do regular segmented scan at block level
        sh_size  =  sh_size / K;
        sgm_size = sgm_size / K;

        barrier(CLK_LOCAL_MEM_FENCE);
        if(sgm_size > 1) 
            segm_scan_reg_block(sh_data, sh_size, sgm_size);
    }

    // write back to the beginning of the histogram

    if(sgm_size > 1)
    for(uint i = TH_ID; i < sh_size/sgm_size; i += get_local_size(0) ) {
        TYPE val = sh_data[ i*sgm_size + sgm_size - 1 ];
        barrier(CLK_LOCAL_MEM_FENCE);
        sh_data[ i ] = val;
    }

/*
    if( TH_ID < (sh_size / sgm_size) ) {
        sh_data[ TH_ID ] = sh_data[ (TH_ID + 1)*sgm_size - 1 ];
    }
*/
    barrier(CLK_LOCAL_MEM_FENCE);
}

/**********************************/
/**** Segmented Irregular Scan ****/
/**********************************/

/**
 * ASSERTS: 
 *    1. ptr & hd are adjusted for each warp
 */
inline REAL segm_scan_warp( 
            __local volatile REAL* ptr,  
            __local volatile FLAG* hd ) {
    const uint th_id = TH_ID & (WARP-1);

    if( th_id >= 1 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : BOP(ptr[th_id-1], ptr[th_id]);
        hd [th_id] = hd[th_id - 1] | hd[th_id];
    }
    if( th_id >= 2 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : BOP(ptr[th_id-2], ptr[th_id]);
        hd [th_id] = hd[th_id - 2] | hd[th_id];
    }
    if( th_id >= 4 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : BOP(ptr[th_id-4], ptr[th_id]);
        hd [th_id] = hd[th_id - 4] | hd[th_id];
    }
    if( th_id >= 8 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : BOP(ptr[th_id-8], ptr[th_id]);
        hd [th_id] = hd[th_id - 8] | hd[th_id];
    }
#if WARP == 32
    if( th_id >= 16 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : BOP(ptr[th_id-16], ptr[th_id]);
        hd [th_id] = hd[th_id - 16] | hd[th_id];
    }
#endif
    return ptr[th_id];
}

inline REAL segm_scan_block( 
            __local volatile REAL *ptr, 
            __local volatile FLAG *hd  
) {
    // 1a: record whether this warp begins 
    //      with an ``open'' segment.
    bool warp_is_open = (hd[WARP_FST] == 0);

    barrier(CLK_LOCAL_MEM_FENCE);

    // 1b: intra-warp segmented scan for each warp
    REAL val = segm_scan_warp( ptr + WARP_FST, hd + WARP_FST );

    // 2a: the last value is the correct partial result
    REAL warp_total = ptr[WARP_LST];
    
    // 2b: warp_flag is the OR-reduction of the flags 
    //     in a warp, and is computed indirectly from
    //     the mindex in hd[]
    bool warp_flag = hd[WARP_LST]!=0 || !warp_is_open;
    bool will_accum= warp_is_open && (hd[TH_ID] == 0);
    
    barrier(CLK_LOCAL_MEM_FENCE);

    // 2c: the last thread in a warp writes partial results
    if( TH_ID == WARP_LST ) {
        ptr[WARP_ID] = warp_total;
        hd [WARP_ID] = warp_flag;
    } 

    barrier(CLK_LOCAL_MEM_FENCE);

    // 3: one warp scans the per-warp results
    if( WARP_ID == 0 ) {
        segm_scan_warp( ptr + WARP_FST, hd + WARP_FST );
    } 

    barrier(CLK_LOCAL_MEM_FENCE);

    // 4: accumulate results from step 3:
    if ( WARP_ID != 0 && will_accum ) {
        val = BOP( val, ptr[ WARP_ID-1 ] );
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    ptr[TH_ID] = val;

    barrier(CLK_LOCAL_MEM_FENCE);

    return val;
}

/**
 * This is `segm_scan_block' when `K==1'. 
 * Otherwise each thread processes `K' 
 *    elements sequentially.
 * IF `warp_level == true' then each warp
 *    (always) begins a new segment, 
 *    i.e., not a real scan.
 * ASSERTS:
 *     1. K <= SEQ_ELMS
 *     2. size of ptr and hd is [K * local_group_size(0)]
 */
inline void segm_scan_block_K( 
                __local volatile REAL *ptr, 
                __local volatile FLAG *hd ,   
                                 uint  K  ,
                                 bool  warp_level
) { 
    if(K == 1 && !warp_level) {
        segm_scan_block(ptr, hd);
    } else {
        uint th_idmK  = TH_ID * K;
        REAL reg_ptr[SEQ_ELMS];
        bool reg_hd [SEQ_ELMS];

        { // 1) each threads copies its elems in regs
            for( uint i = 0; i< K; i++ ) {
                reg_ptr[i] = ptr[th_idmK + i];
                reg_hd [i] = hd [th_idmK + i];
            }
        }

        { // 2) now write to shared memory local_group_size(0) elements
            FLAG flag;
            REAL data; 
            data = reg_ptr[K-1];     
    
            if ( reg_hd[K-1] == 0 ) { // the sum of the last segm
                uint i = K-1;
                do { 
                    i --;
                    data += reg_ptr[i]; 
                } while (i > 0 && reg_hd[i] == 0);

                flag = (i != 0 || reg_hd[0]);
            } else {                  // the last element
                flag = 1;
            }

            barrier(CLK_LOCAL_MEM_FENCE);

            ptr[ TH_ID ] = data;
            hd [ TH_ID ] = flag;

            barrier(CLK_LOCAL_MEM_FENCE);
        }

        { // 3) perform the segmented scan for get_local_size(0) elements
            REAL accum;

            if(warp_level) {
                segm_scan_warp ( ptr + WARP_FST, hd + WARP_FST );
                barrier(CLK_LOCAL_MEM_FENCE);
                accum = (TH_ID & (WARP-1)) ? ptr[TH_ID-1] : 0.0;
            } else {
                segm_scan_block( ptr, hd );
                accum = (TH_ID           ) ? ptr[TH_ID-1] : 0.0;
            }

          // 4) add the contributions of the K-elem to the result
            reg_ptr[K-1] = ptr[TH_ID];

            barrier(CLK_LOCAL_MEM_FENCE);

            for(int i=0; i<K-1; i++) {
                accum = ( reg_hd[i] == 0 ) ? accum + reg_ptr[i] : reg_ptr[i];
                ptr[th_idmK + i] = accum;
            }
            ptr[th_idmK + K-1] = reg_ptr[K-1]; // update last
        }

        barrier(CLK_LOCAL_MEM_FENCE);
    }
}



/*****************************************/
/********** LOCAL BITONIC SORT ***********/
/*****************************************/

inline void ComparatorPrivate(
    uint *keyA,
    REAL *valA,
    uint *keyB,
    REAL *valB,
    uint arrowDir
){
    if( (*keyA > *keyB) == arrowDir ){
        uint t;
        t = *keyA; *keyA = *keyB; *keyB = t;
        REAL r;
        r = *valA; *valA = *valB; *valB = r;
    }
}

// volatile
inline void ComparatorLocal(
    __local  uint *keyA,      
    __local  REAL *valA,
    __local  uint *keyB,
    __local  REAL *valB,
    uint arrowDir
){
    if( (*keyA > *keyB) == arrowDir ){
        uint t;
        t = *keyA; *keyA = *keyB; *keyB = t;
        REAL r;
        r = *valA; *valA = *valB; *valB = r;
    }
}


/**
 * `sort_len' is the length of the periodic sorted array subset, default is 1.
 */
void bitonicSortLocal(
    __local volatile uint *l_key,
    __local volatile REAL *l_val,
    uint arrayLength,
    uint sort_len
    //uint sortDir
){
    // if sorting spans more than one warp than barrier is necessary!
    bool with_barrier = ( arrayLength > (WARP << 1) );

    for(uint size = sort_len*2; size < arrayLength; size <<= 1) {  // size = sort_len*2 //size = 2
        //Bitonic merge
        uint dir = ( (get_local_id(0) & (size / 2)) != 0 );
        for(uint stride = size / 2; stride > 0; stride >>= 1){
            if(with_barrier) 
                barrier(CLK_LOCAL_MEM_FENCE);
            uint pos = 2 * get_local_id(0) - (get_local_id(0) & (stride - 1));
            ComparatorLocal(
                &l_key[pos +      0], &l_val[pos +      0],
                &l_key[pos + stride], &l_val[pos + stride],
                dir
            );
        }
    }

    //dir == sortDir for the last bitonic merge step
    {
        for(uint stride = arrayLength / 2; stride > 0; stride >>= 1){
            if(with_barrier) 
                barrier(CLK_LOCAL_MEM_FENCE);
            uint pos = 2 * get_local_id(0) - (get_local_id(0) & (stride - 1));
            ComparatorLocal(
                &l_key[pos +      0], &l_val[pos +      0],
                &l_key[pos + stride], &l_val[pos + stride],
                1//sortDir
            );
        }
    }

    barrier(CLK_LOCAL_MEM_FENCE);
}

/*********************************************/
/*********************************************/
/*********** Matrix Transposition ************/
/*********************************************/
/*********************************************/

/** 
 * This kernel is optimized to ensure all global reads and writes are coalesced,
 * and to avoid bank conflicts in shared memory.  This kernel is up to 11x faster
 * than the naive kernel below.  Note that the shared memory array is sized to 
 * (BLOCK_DIM+1)*BLOCK_DIM.  This pads each row of the 2D block in shared memory 
 * so that bank conflicts do not occur when threads address the array column-wise.
 */
inline void transpMatrix(
        __global TYPE *idata, 
        __global TYPE *odata, 
                 uint  width, 
                 uint  height,
        __local  TYPE *block,
          const  uint  BLOCK_DIM 
) {
    uint xIndex, yIndex;

    // adjust the input and output arrays for the third dimension!
    yIndex = get_global_id(2) * width * height;
    odata = odata + yIndex;
    idata = idata + yIndex;

    // read the matrix tile into shared memory
	xIndex = get_global_id(0);
	yIndex = get_global_id(1); 
  
	if((xIndex < width) && (yIndex < height))
	{
		uint index_in = yIndex * width + xIndex;
		block[get_local_id(1)*(BLOCK_DIM+1)+get_local_id(0)] = idata[index_in];
	} 
 
	barrier(CLK_LOCAL_MEM_FENCE);

	// write the transposed matrix tile to global memory
	xIndex = get_group_id(1) * BLOCK_DIM + get_local_id(0);
	yIndex = get_group_id(0) * BLOCK_DIM + get_local_id(1);
	if((xIndex < height) && (yIndex < width))
    {
		uint index_out = yIndex * height + xIndex;
		odata[index_out] = block[get_local_id(0)*(BLOCK_DIM+1)+get_local_id(1)];
	} 
}

