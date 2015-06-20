// Various Kinds of Regular and Segmented Reduction
#ifndef REDUCTION_KERNELS
#define REDUCTION_KERNELS

/*************************************************/
/*** Regular Total Reduction With Operator AND ***/
/*************************************************/

inline
uchar reduce_reg_and_warp(
        __local volatile uchar* sh_data
) {
    const uint th_id = get_local_id(0) & (WARP-1);

    if( th_id >= 1 ) sh_data[th_id] = sh_data[th_id] && sh_data[th_id-1 ];
    if( th_id >= 2 ) sh_data[th_id] = sh_data[th_id] && sh_data[th_id-2 ];
    if( th_id >= 4 ) sh_data[th_id] = sh_data[th_id] && sh_data[th_id-4 ];
    if( th_id >= 8 ) sh_data[th_id] = sh_data[th_id] && sh_data[th_id-8 ];
#if (WARP == 32)
    if( th_id >= 16) sh_data[th_id] = sh_data[th_id] && sh_data[th_id-16];
#endif
    return sh_data[th_id];
}

inline 
uchar reduce_reg_and ( __local volatile uchar* sh_data ) {
    // perform scan at warp level
    uchar res = reduce_reg_and_warp( sh_data + WARP_FST );
    barrier( CLK_LOCAL_MEM_FENCE );

    // gather per-warp results in the first WARP elements
    if ( TH_ID == WARP_LST ) {
        sh_data[ WARP_ID ] = res; 
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    
    if( TH_ID < WARP ) {   
        reduce_reg_and_warp( sh_data );
    }
    barrier(CLK_LOCAL_MEM_FENCE);
#if 0
    if( WARP_ID != 0 ) {
        res = res && sh_data[ WARP_ID-1 ];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    sh_data[TH_ID] = res;

    barrier(CLK_LOCAL_MEM_FENCE);

    res = sh_data[ get_local_size(0) - 1 ];
    barrier(CLK_LOCAL_MEM_FENCE);
    return res;
#else
    res = sh_data[ (get_local_size(0) >> lgWARP) - 1 ]; //sh_data[WARP-1];
    barrier(CLK_LOCAL_MEM_FENCE);
    return res;
#endif
}

/*************************************************/
/*** Regular Total Reduction With Operator PLUS***/
/*************************************************/

inline
REAL reduce_reg_plus_warp(
        __local volatile REAL* sh_data
) {
    const uint th_id = get_local_id(0) & (WARP-1);

    if( th_id >= 1 ) sh_data[th_id] = sh_data[th_id] + sh_data[th_id-1 ];
    if( th_id >= 2 ) sh_data[th_id] = sh_data[th_id] + sh_data[th_id-2 ];
    if( th_id >= 4 ) sh_data[th_id] = sh_data[th_id] + sh_data[th_id-4 ];
    if( th_id >= 8 ) sh_data[th_id] = sh_data[th_id] + sh_data[th_id-8 ];
#if (WARP == 32)
    if( th_id >= 16) sh_data[th_id] = sh_data[th_id] + sh_data[th_id-16];
#endif
    return sh_data[th_id];
}

inline 
REAL reduce_reg_plus ( __local volatile REAL* sh_data ) {
    // perform scan at warp level
    REAL res = reduce_reg_plus_warp( sh_data + WARP_FST );
    barrier( CLK_LOCAL_MEM_FENCE );

    // gather per-warp results in the first WARP elements
    if ( TH_ID == WARP_LST ) {
        sh_data[ WARP_ID ] = res; 
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    
    if( TH_ID < WARP ) {   
        reduce_reg_plus_warp( sh_data );
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if( WARP_ID != 0 ) {
        res = res + sh_data[ WARP_ID-1 ];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    sh_data[TH_ID] = res;
    barrier(CLK_LOCAL_MEM_FENCE);    

//    res = sh_data[ get_local_size(0) - 1 ];
//    barrier(CLK_LOCAL_MEM_FENCE);
//    return res;

    return sh_data[ get_local_size(0) - 1 ];
}

/*************************************************/
/*** Irregular, Segmented Reduction With Op. + ***/
/*************************************************/

inline REAL segm_reduce_plus_warp( 
            __local volatile REAL*  ptr,  
            __local volatile uchar* hd ) {
    const uint th_id = TH_ID & (WARP-1);

    if( th_id >= 1 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : ptr[th_id-1] + ptr[th_id];
        hd [th_id] = hd[th_id - 1] | hd[th_id];
    }
    if( th_id >= 2 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : ptr[th_id-2] + ptr[th_id];
        hd [th_id] = hd[th_id - 2] | hd[th_id];
    }
    if( th_id >= 4 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : ptr[th_id-4] + ptr[th_id];
        hd [th_id] = hd[th_id - 4] | hd[th_id];
    }
    if( th_id >= 8 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : ptr[th_id-8] + ptr[th_id];
        hd [th_id] = hd[th_id - 8] | hd[th_id];
    }
#if (WARP == 32)
    if( th_id >= 16 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : ptr[th_id-16]+ ptr[th_id];
        hd [th_id] = hd[th_id - 16] | hd[th_id];
    }
#endif
    return ptr[th_id];
}

inline void segm_reduce_plus( 
            __local volatile REAL  *ptr, 
            __local volatile uchar *hd  
) {
    // 1a: record whether this warp begins 
    //      with an ``open'' segment.
    uchar th_flag      = hd[ TH_ID ];
    bool  warp_is_open = ( hd[WARP_FST] == 0 );
    barrier(CLK_LOCAL_MEM_FENCE);

    // 1b: intra-warp segmented scan for each warp
    REAL val = segm_reduce_plus_warp( ptr + WARP_FST, hd + WARP_FST );

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
        segm_reduce_plus_warp( ptr + WARP_FST, hd + WARP_FST );
    } 

    barrier(CLK_LOCAL_MEM_FENCE);

    // 4: accumulate results from step 3:
    if ( WARP_ID != 0 && will_accum ) {
        val = val + ptr[ WARP_ID-1 ];
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    ptr[TH_ID] = val;
    hd [TH_ID] = th_flag;

    barrier(CLK_LOCAL_MEM_FENCE);
}

/****************************************************/
/*** Irregular, Segmented Reduce With Op. (+,max) ***/
/****************************************************/

inline REAL2 OP_PLUS_MAX(REAL2 a, REAL2 b) {
    return (REAL2) ( a.x + b.x, maxR(a.y, b.y) );
}

inline REAL2 segm_reduce_plusmax_warp( 
            __local volatile REAL2* ptr,  
            __local volatile uchar* hd ) {
    const uint th_id = TH_ID & (WARP-1);

    if( th_id >= 1 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : OP_PLUS_MAX(ptr[th_id-1], ptr[th_id]);
        hd [th_id] = hd[th_id - 1] | hd[th_id];
    }
    if( th_id >= 2 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : OP_PLUS_MAX(ptr[th_id-2], ptr[th_id]);
        hd [th_id] = hd[th_id - 2] | hd[th_id];
    }
    if( th_id >= 4 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : OP_PLUS_MAX(ptr[th_id-4], ptr[th_id]);
        hd [th_id] = hd[th_id - 4] | hd[th_id];
    }
    if( th_id >= 8 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : OP_PLUS_MAX(ptr[th_id-8], ptr[th_id]);
        hd [th_id] = hd[th_id - 8] | hd[th_id];
    }
#if (WARP == 32)
    if( th_id >= 16 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : OP_PLUS_MAX(ptr[th_id-16], ptr[th_id]);
        hd [th_id] = hd[th_id - 16] | hd[th_id];
    }
#endif
    return ptr[th_id];
}

inline void segm_reduce_plusmax( 
            __local volatile REAL2 *ptr, 
            __local volatile uchar *hd  
) {
    // 1a: record whether this warp begins 
    //      with an ``open'' segment.
    uchar th_flag      = hd[ TH_ID ];
    bool  warp_is_open = ( hd[WARP_FST] == 0 );
    barrier(CLK_LOCAL_MEM_FENCE);

    // 1b: intra-warp segmented scan for each warp
    REAL2 val = segm_reduce_plusmax_warp( ptr + WARP_FST, hd + WARP_FST );

    // 2a: the last value is the correct partial result
    REAL2 warp_total = ptr[WARP_LST];
    
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
        segm_reduce_plusmax_warp( ptr + WARP_FST, hd + WARP_FST );
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    // 4: accumulate results from step 3:
    if ( WARP_ID != 0 && will_accum ) {
        val = OP_PLUS_MAX( val, ptr[ WARP_ID-1 ] );
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    ptr[TH_ID] = val;
    hd [TH_ID] = th_flag;

    barrier(CLK_LOCAL_MEM_FENCE);
}


/****************************************************/
/*** Irregular, Segmented Reduce With Op. (+,max) ***/
/****************************************************/

inline REAL4 OP_PLUS_MIN_MAX( REAL4 a, REAL4 b ) {
    return ( REAL4 ) ( a.x + b.x, minR(a.y, b.y), maxR(a.z, b.z), 1.0 );
}

inline REAL4 segm_reduce_plusminmax_warp( 
            __local volatile REAL4* ptr,  
            __local volatile uchar* hd ) {
    const uint th_id = TH_ID & (WARP-1);

    if( th_id >= 1 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : OP_PLUS_MIN_MAX(ptr[th_id-1], ptr[th_id]);
        hd [th_id] = hd[th_id - 1] | hd[th_id];
    }
    if( th_id >= 2 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : OP_PLUS_MIN_MAX(ptr[th_id-2], ptr[th_id]);
        hd [th_id] = hd[th_id - 2] | hd[th_id];
    }
    if( th_id >= 4 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : OP_PLUS_MIN_MAX(ptr[th_id-4], ptr[th_id]);
        hd [th_id] = hd[th_id - 4] | hd[th_id];
    }
    if( th_id >= 8 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : OP_PLUS_MIN_MAX(ptr[th_id-8], ptr[th_id]);
        hd [th_id] = hd[th_id - 8] | hd[th_id];
    }
#if (WARP == 32)
    if( th_id >= 16 ) {
        ptr[th_id] = hd[th_id] ? ptr[th_id] : OP_PLUS_MIN_MAX(ptr[th_id-16], ptr[th_id]);
        hd [th_id] = hd[th_id - 16] | hd[th_id];
    }
#endif
    return ptr[th_id];
}

inline void segm_reduce_plusminmax( 
            __local volatile REAL4 *ptr, 
            __local volatile uchar *hd  
) {
    // 1a: record whether this warp begins 
    //      with an ``open'' segment.
    uchar th_flag      = hd[ TH_ID ];
    bool  warp_is_open = ( hd[WARP_FST] == 0 );
    barrier(CLK_LOCAL_MEM_FENCE);

    // 1b: intra-warp segmented scan for each warp
    REAL4 val = segm_reduce_plusminmax_warp( ptr + WARP_FST, hd + WARP_FST );

    // 2a: the last value is the correct partial result
    REAL4 warp_total = ptr[WARP_LST];
    
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
        segm_reduce_plusminmax_warp( ptr + WARP_FST, hd + WARP_FST );
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    // 4: accumulate results from step 3:
    if ( WARP_ID != 0 && will_accum ) {
        val = OP_PLUS_MIN_MAX( val, ptr[ WARP_ID-1 ] );
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    ptr[TH_ID] = val;
    hd [TH_ID] = th_flag;

    barrier(CLK_LOCAL_MEM_FENCE);
}

#endif // REDUCTION_KERNELS
