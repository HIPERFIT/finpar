/******************************************************/
/********* FIND_BEST KERNEL and HELPERS ***************/
/******************************************************/

inline REAL2 OP_BEST(REAL2 iv1, REAL2 iv2) {
    return ( ( iv1.y > iv2.y ) ? iv1 : iv2 ); 
}

inline REAL2 reduce_best_warp( __local volatile REAL2* sh_data ) {
    const uint th_id = get_local_id(0) & (WARP-1);

    if( th_id >= 1  ) sh_data[th_id] = OP_BEST( sh_data[th_id], sh_data[th_id-1 ] );
    if( th_id >= 2  ) sh_data[th_id] = OP_BEST( sh_data[th_id], sh_data[th_id-2 ] );
    if( th_id >= 4  ) sh_data[th_id] = OP_BEST( sh_data[th_id], sh_data[th_id-4 ] );
    if( th_id >= 8  ) sh_data[th_id] = OP_BEST( sh_data[th_id], sh_data[th_id-8 ] );
#if (WARP == 32)
    if( th_id >= 16 ) sh_data[th_id] = OP_BEST( sh_data[th_id], sh_data[th_id-16] );
#endif
    return sh_data[th_id];
}

inline void reduce_best_local ( __local volatile REAL2* sh_data ) {
    // perform scan at warp level
    REAL2 res = reduce_best_warp( sh_data + WARP_FST );
    barrier( CLK_LOCAL_MEM_FENCE );

    // gather per-warp results in the first WARP elements
    if ( TH_ID == WARP_LST ) {
        sh_data[ WARP_ID ] = res; 
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    
    if( TH_ID < WARP ) {   
        reduce_best_warp( sh_data );
    }
    barrier(CLK_LOCAL_MEM_FENCE);
}


__kernel void 
redbest_ker1 (
    __global REAL *arr,
    __global uint *res_ind,
    __global REAL *res_val,
    uint           offset,
    uint           N,
    uint           U
) {
    __local REAL2 locred[ LWG_FB ];

    arr += offset;

    REAL2 seqbest = (REAL2) ( 0.0, -INFTY );

    uint  i = get_group_id(0) * (U * get_local_size(0)) + get_local_id(0);

    
    //uint  B = min( i + (U * get_local_size(0)), N );
    uint  B = i + (U * get_local_size(0));
    if( B > N ) { B = N; }

    for( ; i < B; i += get_local_size(0) ) {
        REAL2 cur = (REAL2) ( (REAL)i, arr[i] );
        seqbest = OP_BEST(seqbest, cur);
    }

    locred[ get_local_id(0) ] = seqbest;
    barrier(CLK_LOCAL_MEM_FENCE);
    
    reduce_best_local ( locred );

    if ( get_local_id(0) == get_local_size(0) - 1 ) {
        REAL2 res = locred[WARP-1];
        res_ind[ get_group_id(0) ] = (uint) res.x;
        res_val[ get_group_id(0) ] = res.y;
    }
}

/**
 * Assumes we have to reduce exactly one local workgroup.
 */
__kernel void 
redbest_ker2 (
    __global uint *arr_ind,
    __global REAL *arr_val,
    uint           N
) {
    __local REAL2 locred[ LWG_FB ];

    uint th_id = get_local_id(0);

    if ( th_id >= N ) {
        REAL2 tmp = (REAL2) ( 0.0, -INFTY );
        locred[ th_id ] = tmp;
    } else {
        REAL2 tmp = (REAL2) ( (REAL)arr_ind[th_id], arr_val[th_id] );
        locred[ th_id ] = tmp;
    }

    barrier(CLK_LOCAL_MEM_FENCE);    
    reduce_best_local( locred );

    if ( th_id == get_local_size(0) - 1 ) {
        REAL2 res = locred[ WARP-1 ];
        arr_ind[ 0 ] = (uint) res.x;
        arr_val[ 0 ] = res.y;
    }
}

