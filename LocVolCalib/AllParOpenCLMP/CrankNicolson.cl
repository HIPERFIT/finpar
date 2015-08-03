#include "../includeC/Constants.h"
 
    
#if (WITH_FLOATS)
    typedef float4 REAL4;
    typedef float2 REAL2;
    typedef float3 REAL3;
#else
    typedef double4 REAL4;
    typedef double2 REAL2;
    typedef double3 REAL3;
#endif

/**************************************************************************/
/**************** SHARED-MEMORY ACCESSOR FUNCTIONS ************************/
/**************************************************************************/

#define WITH_INTERLEAVED_BANKS4 1 
#define WITH_INTERLEAVED_BANKS2 1

#if WITH_INTERLEAVED_BANKS2
inline void store2(REAL2 val, unsigned int ind,  __local REAL* cache) {
    ind = (ind & (WORKGROUP_SIZE-1)) + ( (ind >> logWORKGROUP_SIZE) << (logWORKGROUP_SIZE+1) );
    cache[ind] = val.x;
    ind       += WORKGROUP_SIZE;
    cache[ind] = val.y;
}

inline REAL2 load2(unsigned int ind,  __local REAL* cache) {
    ind = (ind & (WORKGROUP_SIZE-1)) + ( (ind >> logWORKGROUP_SIZE) << (logWORKGROUP_SIZE+1) );
    REAL2 val;
    val.x = cache[ind];
    ind  += WORKGROUP_SIZE;
    val.y = cache[ind];
    return val;
}
#else
inline void store2(REAL2 val, unsigned int ind, volatile __local REAL2* cache) {
    cache[ind] = val;
}

inline REAL2 load2(unsigned int ind, volatile __local REAL2* cache) {
    return cache[ind];
}
#endif


#if WITH_INTERLEAVED_BANKS4
inline void store4(REAL4 val, unsigned int ind,  __local REAL* cache) {
    //ind = (ind / WORKGROUP_SIZE)*(WORKGROUP_SIZE) * 4 + (ind % WORKGROUP_SIZE);
    ind = (ind & (WORKGROUP_SIZE-1)) + ( (ind >> logWORKGROUP_SIZE) << (logWORKGROUP_SIZE+2) );
    cache[ind] = val.x;
    ind       += WORKGROUP_SIZE;
    cache[ind] = val.y;
    ind       += WORKGROUP_SIZE;
    cache[ind] = val.z;
    ind       += WORKGROUP_SIZE;
    cache[ind] = val.w;
}

inline REAL4 load4(unsigned int ind,  __local REAL* cache) {
    REAL4 val;
    //ind = (ind / WORKGROUP_SIZE)*(WORKGROUP_SIZE) * 4 + (ind % WORKGROUP_SIZE);
    ind = (ind & (WORKGROUP_SIZE-1)) + ( (ind >> logWORKGROUP_SIZE) << (logWORKGROUP_SIZE+2) );
    val.x = cache[ind];
    ind  += WORKGROUP_SIZE;
    val.y = cache[ind];
    ind  += WORKGROUP_SIZE;
    val.z = cache[ind];
    ind  += WORKGROUP_SIZE;
    val.w = cache[ind];
    return val;
}
#else
inline void store4(REAL4 val, unsigned int ind, volatile __local REAL4* cache) {
    cache[ind] = val;
}

inline REAL4 load4(unsigned int ind, volatile __local REAL4* cache) {
    return cache[ind];
}
#endif
/**************************************************************************/
/*********** PREPARE FOR TRIDAG X *****************************************/
/**************************************************************************/
__kernel void prepare_tridag_x (
        /*** Read-Only Scalars ***/
        __global RWScalars*   ro_scals,   // __constant
        /*** Read-Only Arrays  ***/
        __global REAL*         myX,
        __global REAL3*        myDx,
        __global REAL3*        myDxx,
        __constant REAL*         myY,
        __constant REAL*        myDy,
        __constant REAL*        myDyy,
        /*** Read-Write Temporary Arrays ***/
        __global   REAL*        a,
        __global   REAL*        b,
        __global   REAL*        c,
        __global   REAL*        y,
        __global   REAL*        u,
        __global   REAL*        v,
        __global   REAL4*       scan_tmp,
        /*** Temporary Array in Shared Space ***/   
        __local    REAL*        cache_tmp,
        /*** The Result Array ***/
        __global   REAL*        res_arr        
) {
    REAL3 myDy_elem, myDyy_elem, myres_elem;
    REAL cur_myVar = ro_scals->nu; 
    REAL cur_myMu  = 0.0;

    unsigned int glob_ind = get_global_id(2) * get_global_size(1) * get_global_size(0);
    const unsigned int ind = glob_ind + get_global_id(1)*get_global_size(0) + get_global_id(0);

    REAL tmp = 0.0;
 
    cur_myVar *= cur_myVar;
 
    // second loop    
    {
        myres_elem = (REAL3)(   (get_global_id(1)!=0) ? res_arr[ind - get_global_size(0)] : 0.0, 
                                res_arr[ind], 
                                (get_global_id(1) != get_global_size(1)-1) ? res_arr[ind + get_global_size(0)] : 0.0
                            );
 
        myDy_elem  = (REAL3)( myDy[get_global_id(1)<<2], myDy[ (get_global_id(1)<<2)+1], myDy[(get_global_id(1)<<2)+2] );
        myDyy_elem = (REAL3)( myDyy[get_global_id(1)<<2], myDyy[(get_global_id(1)<<2)+1], myDyy[(get_global_id(1)<<2)+2] );

        myDy_elem = cur_myMu*myDy_elem + 0.5f*cur_myVar*myDyy_elem;
        myDy_elem *= myres_elem;
        tmp += myDy_elem.x + myDy_elem.y + myDy_elem.z;

        v[glob_ind + get_global_id(0)*get_global_size(1) + get_global_id(1)] = tmp;
    }

    glob_ind = get_local_id(2)*get_local_size(1)*get_local_size(0) + get_local_id(1)*get_local_size(0) + get_local_id(0);

    // first loop 
    cur_myMu  = 0.0; // X
    cur_myVar = exp( 2 * ( ro_scals->beta*log(myX[get_global_id(0)]) + 
                                 myY[get_global_id(1)] - 0.5f*cur_myVar*ro_scals->timeline_i  // cur_myVar == nu*nu
                          )
                    ); 
    { 
        // CACHING ASSUMES ALL ELEMENTS IN X DIMENSION fit in the LOCALGROUP!
        cache_tmp[glob_ind] = res_arr[ ind ];

        myDy_elem  = myDx[ get_global_id(0) ]; 
        myDyy_elem = myDxx[ get_global_id(0) ];     

        barrier(CLK_LOCAL_MEM_FENCE);   


        myres_elem = (REAL3)(     (get_global_id(0)!=0) ? cache_tmp[ glob_ind - 1 ] : 0.0, 
                                  cache_tmp[ glob_ind ], 
                                  (get_global_id(0) != get_global_size(0)-1) ? cache_tmp[ glob_ind + 1 ] : 0.0
                            );

        tmp += ro_scals->dtInv*myres_elem.y;
        myDy_elem = 0.5f*( cur_myMu*myDy_elem + 0.5f*cur_myVar*myDyy_elem );

        c[ind] = - myDy_elem.z;   

        barrier(CLK_GLOBAL_MEM_FENCE);

        a[ind] = - myDy_elem.x; // b[ind] = ro_scals->dtInv - myDy_elem.y;

        myDyy_elem = myDy_elem*myres_elem;
        tmp += myDyy_elem.x + myDyy_elem.y + myDyy_elem.z;
    }    

    u[ind] = tmp; 

    

    tmp = ro_scals->dtInv - myDy_elem.y;  

    // third loop 
    if( get_global_id(0) > 0 ) { 
        scan_tmp[ind] = (REAL4) ( tmp, myDy_elem.x * c[ind-1], 1.0, 0.0 );
    } else { // i == get_global_id(0) == 0
        y[ind] = tmp;    // = b [ind]
        scan_tmp[ind] = (REAL4) ( tmp, 0.0, 0.0, tmp ); 
    }
} 


/**************************************************************************/
/*********** IN-BETWEEN TRIDAGS -- i.e., prepare for TRIDAG Y *************/
/**************************************************************************/
//#define USE_LOCAL_Y 

__kernel void prepare_tridag_y (
        /*** Read-Only Scalars ***/
        __global RWScalars*   ro_scals,  // __constant
        /*** Read-Only Arrays  ***/
        __global REAL3*        myDy,
        __global REAL3*        myDyy,
        /*** Read-Write Temporary Arrays ***/
        __global   REAL*        a,
        __global   REAL*        b,
        __global   REAL*        c,
        __global   REAL*        u,
        __global   REAL*        v,
        __global   REAL*        y,
        __global   REAL4*       scan_tmp,
        /*** Temporary Array in Shared Space ***/   
        __local    REAL*        cache_tmp
) {
 
    // !!!0 <-> j; 1<-> i; 2 <-> k!!!
    const unsigned int ind       = get_global_id(2) * get_global_size(1) * get_global_size(0)
                                    + get_global_id(1) * get_global_size(0) + get_global_id(0); 
                                    // + get_global_id(1)*NUM_Y + get_global_id(0); 
    const unsigned int glob_ind = get_local_id(2)*get_local_size(1)*get_local_size(0) 
                                    + get_local_id(1)*get_local_size(0) + get_local_id(0);
    REAL3 myDy_elem;
    REAL dt_inv = ro_scals->dtInv;

    REAL cur_myMuY = 0.0;
    REAL cur_myVarY = ro_scals->nu; cur_myVarY *= cur_myVarY*0.5f; //nu*nu;

    myDy_elem  = 0.5f*( cur_myMuY*myDy [get_global_id(0)] + cur_myVarY*myDyy[get_global_id(0)] );

#ifdef USE_LOCAL_Y
    cache_tmp[glob_ind] = - myDy_elem.z;
#endif

    a[ind] = - myDy_elem.x; c[ind] = - myDy_elem.z;

    // NUM_X = get_global_size(1)
    v[ind] =    dt_inv * 
                u[  get_global_id(2)* get_global_size(1) * get_global_size(0) + 
                    get_global_id(0)*get_global_size(1) + get_global_id(1)     ] 
                - 0.5f*v[ind];  // y[j] = dtInv*u[j][i] - 0.5*v[i][j];

#ifdef USE_LOCAL_Y
    barrier(CLK_LOCAL_MEM_FENCE);
#endif

    // third loop 
    if( get_global_id(0) > 0 ) { 
#ifdef USE_LOCAL_Y
        scan_tmp[ind] = (REAL4) ( dt_inv - myDy_elem.y, myDy_elem.x * cache_tmp[glob_ind-1], 1.0, 0.0 );
#else
        const REAL elem = - 0.5f*(cur_myMuY*myDy[get_global_id(0)-1].z + cur_myVarY*myDyy[get_global_id(0)-1].z);
        scan_tmp[ind] = (REAL4) ( dt_inv - myDy_elem.y, myDy_elem.x * elem, 1.0, 0.0 );   
#endif
        //scan_tmp[ind] = (REAL4) ( dt_inv - myDy_elem.y, myDy_elem.x*c[ind-1], 1.0, 0.0 );

    } else { // j == get_global_id(0) == 0
        y[ind] = dt_inv - myDy_elem.y; //b[ind];
        scan_tmp[ind] = (REAL4) ( dt_inv - myDy_elem.y, 0.0, 0.0, dt_inv - myDy_elem.y );
    }
}


/**********************************************************************/
/************* TRIDAG KERNELS *****************************************/
/**********************************************************************/

/*************************************/
/*** 1. Scan with Matrix Multiply ****/
/*************************************/
 
/**
 * Multiplies 2x2 matrixes `a' and `b' and
 * stores the result in(-place in) `a'.
 */ 
inline REAL4 matmult2(REAL4 a, REAL4 b) {
    REAL val = (a.x*b.x + a.y*b.z);

    a = (REAL4) (val, a.x*b.y + a.y*b.w, a.z*b.x + a.w*b.z, a.z*b.y + a.w*b.w);
    val = 1.0/val;
    return a * val;
}

 

//#define LOG2_WARP_SIZE 5U
//#define WARP_SIZE (1U << LOG2_WARP_SIZE)
 
#if WITH_INTERLEAVED_BANKS4
inline REAL4 warpScanInclMatMult(REAL4 idata, volatile __local REAL4 *l_Data4, uint size){
    volatile __local REAL* l_Data = (volatile __local REAL*) l_Data4;

    uint pos = 2 * get_local_id(0) - (get_local_id(0) & (size - 1));
    store4( (REAL4)(1.0, 0.0, 0.0, 1.0), pos, l_Data );
    pos += size;
    store4( idata, pos, l_Data );  

    if(size >=  2) { REAL4 tmp = matmult2( load4(pos, l_Data),load4(pos-1, l_Data) ); store4(tmp, pos, l_Data); }
    if(size >=  4) { REAL4 tmp = matmult2( load4(pos, l_Data),load4(pos-2, l_Data) ); store4(tmp, pos, l_Data); }
    if(size >=  8) { REAL4 tmp = matmult2( load4(pos, l_Data),load4(pos-4, l_Data) ); store4(tmp, pos, l_Data); }
    if(size >= 16) { REAL4 tmp = matmult2( load4(pos, l_Data),load4(pos-8, l_Data) ); store4(tmp, pos, l_Data); }
#if WARP == 32
    if(size >= 32) { REAL4 tmp = matmult2( load4(pos, l_Data),load4(pos-16,l_Data) ); store4(tmp, pos, l_Data); }
#endif
    return load4(pos, l_Data);
}
#else
inline REAL4 warpScanInclMatMult(REAL4 idata, volatile __local REAL4 *l_Data, uint size){
    uint pos = 2 * get_local_id(0) - (get_local_id(0) & (size - 1));
    l_Data[pos] = (REAL4)(1.0, 0.0, 0.0, 1.0); 
    pos += size;
    l_Data[pos] = idata;  

    if(size >=  2) l_Data[pos] = matmult2( l_Data[pos], l_Data[pos-1] );  
    if(size >=  4) l_Data[pos] = matmult2( l_Data[pos], l_Data[pos-2] );  
    if(size >=  8) l_Data[pos] = matmult2( l_Data[pos], l_Data[pos-4] );  
    if(size >= 16) l_Data[pos] = matmult2( l_Data[pos], l_Data[pos-8] );  
#if WARP == 32
    if(size >= 32) l_Data[pos] = matmult2( l_Data[pos], l_Data[pos-16]);  
#endif
    return l_Data[pos];
}
#endif
  
//Vector scan: the array to be scanned is stored
//in work-item private memory as uint4
inline REAL4 scanMatMultInclLocal(REAL4 idata4, __local REAL4 *l_Data, uint size){
    if(size > WARP){
        //Bottom-level inclusive warp scan
        REAL4 warpResult = warpScanInclMatMult(idata4, l_Data, WARP);
 
        //Save top elements of each warp for exclusive warp scan
        //sync to wait for warp scans to complete (because l_Data is being overwritten)
        barrier(CLK_LOCAL_MEM_FENCE);
        if( (get_local_id(0) & (WARP - 1)) == (WARP - 1) ) {
#if WITH_INTERLEAVED_BANKS4
            store4(warpResult, get_local_id(0) >> lgWARP, (__local REAL*)l_Data);
#else
            l_Data[get_local_id(0) >> lgWARP] = warpResult;
#endif
        }
        //wait for warp scans to complete
        barrier(CLK_LOCAL_MEM_FENCE); 
     
        if( get_local_id(0) < (WORKGROUP_SIZE >> lgWARP) ){
#if WITH_INTERLEAVED_BANKS4
            REAL4 val = load4(get_local_id(0), (__local REAL*)l_Data);    
            store4( warpScanInclMatMult(val, l_Data, size >> lgWARP), get_local_id(0), (__local REAL*)l_Data);
#else
            //grab top warp elements
            REAL4 val = l_Data[get_local_id(0)];
            //calculate inclusive scan and write back to shared memory
            l_Data[get_local_id(0)] = warpScanInclMatMult(val, l_Data, size >> lgWARP);
#endif
        }
    
        //return updated warp scans with exclusive scan results
        barrier(CLK_LOCAL_MEM_FENCE);

        //if( (get_local_id(0) >> LOG2_WARP_SIZE) > 0 )
        if( ( (get_local_id(0) & (size-1)) >> lgWARP ) > 0 ) {
#if WITH_INTERLEAVED_BANKS4
            REAL4 tmp = load4( (get_local_id(0) >> lgWARP) - 1,   (__local REAL*)l_Data );   
            warpResult = matmult2( warpResult, tmp );
#else
            warpResult = matmult2( warpResult, l_Data[ (get_local_id(0) >> lgWARP) - 1] );
#endif
        }

        return warpResult; 
    }else{
        return warpScanInclMatMult(idata4, l_Data, size);
    }
} 

__kernel __attribute__((reqd_work_group_size(WORKGROUP_SIZE, 1, 1)))
void scanMatMultIncl(
    __global REAL4 *data,
    __local  REAL4 *l_Data,
    uint size
){
    //Load data
    REAL4 idata4 = data[get_global_id(0)];

    //Calculate exclusive scan
    REAL4 odata4  = scanMatMultInclLocal( idata4, l_Data, size );

    //Write back
    data[get_global_id(0)] = odata4; 
}


/*************************************************/
/*** 1. Scan with Linear Function Composition ****/
/*************************************************/

// Linear Function Composition and Scan
inline REAL2 linfuncomp(REAL2 ab2, REAL2 ab1) {
    ab2.x = ab2.y * ab1.x + ab2.x;
    ab2.y = ab2.y * ab1.y;
    return ab2;
}

#if WITH_INTERLEAVED_BANKS2
inline REAL2 warpScanInclLinFunComp(REAL2 idata, volatile __local REAL2 *l_Data2, uint size){
    volatile __local REAL* l_Data = (volatile __local REAL*) l_Data2;

    uint pos = 2 * get_local_id(0) - (get_local_id(0) & (size - 1));
    store2( (REAL2)(0.0, 1.0), pos, l_Data );
    pos += size;
    store2( idata, pos, l_Data );  

    if(size >=  2) { REAL2 tmp = linfuncomp( load2(pos, l_Data),load2(pos-1, l_Data) ); store2(tmp, pos, l_Data); }
    if(size >=  4) { REAL2 tmp = linfuncomp( load2(pos, l_Data),load2(pos-2, l_Data) ); store2(tmp, pos, l_Data); }
    if(size >=  8) { REAL2 tmp = linfuncomp( load2(pos, l_Data),load2(pos-4, l_Data) ); store2(tmp, pos, l_Data); }
    if(size >= 16) { REAL2 tmp = linfuncomp( load2(pos, l_Data),load2(pos-8, l_Data) ); store2(tmp, pos, l_Data); }
#if WARP == 32
    if(size >= 32) { REAL2 tmp = linfuncomp( load2(pos, l_Data),load2(pos-16,l_Data) ); store2(tmp, pos, l_Data); }
#endif
    return load2(pos, l_Data);
}
#else
inline REAL2 warpScanInclLinFunComp(REAL2 idata, volatile __local REAL2 *l_Data, uint size){
    uint pos = 2 * get_local_id(0) - (get_local_id(0) & (size - 1));
    l_Data[pos] = (REAL2)(0.0, 1.0);
    pos += size;
    l_Data[pos] = idata;

    if(size >=  2) l_Data[pos] = linfuncomp( l_Data[pos], l_Data[pos-1] );  
    if(size >=  4) l_Data[pos] = linfuncomp( l_Data[pos], l_Data[pos-2] );  
    if(size >=  8) l_Data[pos] = linfuncomp( l_Data[pos], l_Data[pos-4] );   
    if(size >= 16) l_Data[pos] = linfuncomp( l_Data[pos], l_Data[pos-8] );
#if WARP == 32
    if(size >= 32) l_Data[pos] = linfuncomp( l_Data[pos], l_Data[pos-16]);
#endif
    return l_Data[pos];
}
#endif

inline REAL2 scanLinFunCompInclLocal(REAL2 idata4, __local REAL2 *l_Data, uint size){
    if(size > WARP){
        //Bottom-level inclusive warp scan
        REAL2 warpResult = warpScanInclLinFunComp(idata4, l_Data, WARP);

        //Save top elements of each warp for exclusive warp scan
        //sync to wait for warp scans to complete (because l_Data is being overwritten)
        barrier(CLK_LOCAL_MEM_FENCE);
        if( (get_local_id(0) & (WARP - 1)) == (WARP - 1) ) { 
#if WITH_INTERLEAVED_BANKS2
            store2(warpResult, get_local_id(0) >> lgWARP, (__local REAL*)l_Data);
#else
            l_Data[get_local_id(0) >> lgWARP] = warpResult;
#endif
        }
 
        //wait for warp scans to complete
        barrier(CLK_LOCAL_MEM_FENCE);
  
        if( get_local_id(0) < (WORKGROUP_SIZE / WARP) ){
#if WITH_INTERLEAVED_BANKS2
            REAL2 val = load2(get_local_id(0), (__local REAL*)l_Data);    
            store2( warpScanInclLinFunComp(val, l_Data, size >> lgWARP), get_local_id(0), (__local REAL*)l_Data);
#else
            //grab top warp elements
            REAL2 val = l_Data[get_local_id(0)];
            //calculate inclusive scan and write back to shared memory
            l_Data[get_local_id(0)] = warpScanInclLinFunComp(val, l_Data, size >> lgWARP);
#endif
        }

        //return updated warp scans with exclusive scan results
        barrier(CLK_LOCAL_MEM_FENCE);

        if( ( (get_local_id(0) & (size-1)) >> lgWARP ) > 0 ) {
#if WITH_INTERLEAVED_BANKS2
            REAL2 tmp = load2( (get_local_id(0) >> lgWARP) - 1,   (__local REAL*)l_Data );   
            warpResult = linfuncomp( warpResult, tmp );
#else
            warpResult = linfuncomp( warpResult, l_Data[ (get_local_id(0) >> lgWARP) - 1] );
#endif
        }
        return warpResult;
    }else{
        return warpScanInclLinFunComp(idata4, l_Data, size);
    }
} 
  

__kernel __attribute__((reqd_work_group_size(WORKGROUP_SIZE, 1, 1)))
void scanLinFunCompIncl(
    __global REAL2 *data,
    __local  REAL2 *l_Data,
    uint size
){
    //Load data
    REAL2 idata4 = data[get_global_id(0)];

    //Calculate exclusive scan
    REAL2 odata4  = scanLinFunCompInclLocal( idata4, l_Data, size );

    //Write back
    data[get_global_id(0)] = odata4;
}


/*************************************/
/*** 2. TRIDAG - in between scans ****/
/*************************************/

inline REAL map_matmult(REAL4 tmp, REAL val1) {

    REAL nom   = tmp.x*val1+tmp.y;
    REAL denom = tmp.z*val1+tmp.w;

    return nom/denom;

    //REAL denom = tmp.z*val1+tmp.w*val2;
    //REAL nom   = (tmp.x/denom)*val1 + (tmp.y/denom)*val2;
    //return nom;
}

__kernel void conclude_matmult (
        /*** Read-Write Temporary Arrays ***/
        __global   REAL*        b,
        __global   REAL*        d,
        __global   REAL*        y,
        __global   REAL*        u,
        __global   REAL4*       scan_tmp
) {
    unsigned int ind = get_global_id(2) * get_global_size(1) * get_global_size(0) 
                                + get_global_id(1) * get_global_size(0);
    REAL b0 = u[ind]; //b[ ind ];
    ind += get_global_id(0);

    if(get_global_id(0) > 0) {
        u[ind] = map_matmult( scan_tmp[ind], b0 );
    } else {
        y[ind] = d[ind];
    }
}

 
__kernel void prelude_fwd_fun_comp (
        /*** Read-Write Temporary Arrays ***/
        __global   REAL*        a,
        __global   REAL*        d,
        __global   REAL*        u,
        __global   REAL2*       scan_tmp
) {
    const unsigned int ind = get_global_id(2) * get_global_size(1) * get_global_size(0)  
                                + get_global_id(1) * get_global_size(0) + get_global_id(0);

    if( get_global_id(0) > 0 ) { 
        scan_tmp[ind] = (REAL2) ( d[ind], (0.0 - (a[ind]/u[ind-1])) );
    } else {
        scan_tmp[ind] = (REAL2) (0.0, 1.0);
    }
}

/////   

__kernel void post_fwd_fun_comp (
        /*** Read-Write Temporary Arrays ***/
        __global   REAL*   y,
        __global   REAL2*  scan_tmp
) {
    unsigned int ind = get_global_id(2) * get_global_size(1) * get_global_size(0) 
                                + get_global_id(1) * get_global_size(0);

    REAL y0 = y[ ind ];

    ind += get_global_id(0);

    if( get_global_id(0) > 0 ) { 
        REAL2 fun = scan_tmp[ind];
        y[ind] = fun.x + y0 * fun.y;
    }
}

__kernel void prelude_bwd_fun_comp (
        /*** Read-Write Temporary Arrays ***/
        __global   REAL*        c,
        __global   REAL*        y,
        __global   REAL*        u,
        __global   REAL2*       scan_tmp
) {
    const unsigned int ind     = get_global_id(2) * get_global_size(1) * get_global_size(0) 
                                    + get_global_id(1) * get_global_size(0) + get_global_id(0);
    const unsigned int inv_ind = ind - 2*get_global_id(0) + get_global_size(0) - 2; 

    if( get_global_id(0) <  get_global_size(0) - 1 ) {
        scan_tmp[ind] = (REAL2)(y[inv_ind]/u[inv_ind], 0.0 - c[inv_ind]/u[inv_ind]);
    } else {    // get_global_id(0) == get_global_size(0) - 1
        scan_tmp[ind] = (REAL2)(0.0, 1.0);
        y[ind] = y[ind]/u[ind];
    }
}
 
/**
 * 1. map back the result from the (bacward) parallel prefix 
 *      with linear function composition
 * 2. update the RW scalars rw_scals->dtInv and rw_scals->timeline_i
 */
__kernel void post_bwd_fun_comp (
        __global   REAL*        y,
        __global   REAL2*       scan_tmp
) {
    // map back the result from the (bacward) parallel prefix 
    //    with linear function composition

    const unsigned int ind     = get_global_id(2) * get_global_size(1) * get_global_size(0)  
                                    + get_global_id(1) * get_global_size(0) + get_global_id(0);
    const unsigned int inv_ind = ind - 2*get_global_id(0) + get_global_size(0) - 2; 

    REAL ynm1 = y[inv_ind + get_global_id(0) + 1];

    if( get_global_id(0) < get_global_size(0) - 1 ) { 
        REAL2 fun = scan_tmp[ind];
        y[inv_ind] = fun.x + ynm1 * fun.y;
    }
}


/*********************************************/
/*********** Matrix Transposition ************/
/*********************************************/

// This kernel is optimized to ensure all global reads and writes are coalesced,
// and to avoid bank conflicts in shared memory.  The shared memory array is sized  
// to (BLOCK_DIM+1)*BLOCK_DIM.  This pads each row of the 2D block in shared memory 
// so that bank conflicts do not occur when threads address the array column-wise.
inline void transposeMatrix(
        __global REAL *odata, 
        __global REAL *idata, 
        unsigned int width, 
        unsigned int height, 
        __local REAL* block)
{
    unsigned int xIndex, yIndex;

    // adjust the input and output arrays for the third dimension!
    yIndex = get_global_id(2) * width * height;
    odata = odata + yIndex;
    idata = idata + yIndex;

    // read the matrix tile into shared memory
	xIndex = get_global_id(0);
	yIndex = get_global_id(1); 
  
	if((xIndex < width) && (yIndex < height))
	{
		unsigned int index_in = yIndex * width + xIndex;
		block[get_local_id(1)*(BLOCK_DIM+1)+get_local_id(0)] = idata[index_in];
	} 
 
	barrier(CLK_LOCAL_MEM_FENCE);

	// write the transposed matrix tile to global memory
	xIndex = get_group_id(1) * BLOCK_DIM + get_local_id(0);
	yIndex = get_group_id(0) * BLOCK_DIM + get_local_id(1);
	if((xIndex < height) && (yIndex < width))
    {
		unsigned int index_out = yIndex * height + xIndex;
		odata[index_out] = block[get_local_id(0)*(BLOCK_DIM+1)+get_local_id(1)];
	} 
}
   
__kernel void transpose(
        __global REAL *odata, 
        __global REAL *idata, 
        //int offset, 
        unsigned int width, 
        unsigned int height, 
        __local REAL* block)
{
    transposeMatrix( odata, idata, width, height, block );
}
 
__kernel void transposeUpdateScalars(
        __global REAL *odata, 
        __global REAL *idata, 
        //int offset, 
        unsigned int width, 
        unsigned int height, 
        __local REAL* block,
        __global   RWScalars*   ro_scals,  // __constant
        __constant REAL*        timeline
) {
    transposeMatrix( odata, idata, width, height, block );
    //odata[get_global_id(2)*get_global_size(1)*get_global_size(0) + get_global_id(1)*get_global_size(0) + get_global_id(0)] 
    //    = 33.33;

    // update time-loop-variant scalars
    if(get_global_id(2) + get_global_id(1) + get_global_id(0) == 0) {
        int t_ind =  --ro_scals->t_ind;
        ro_scals->dtInv      = 1.0 / ( timeline[t_ind+1] - timeline[t_ind] );
        ro_scals->timeline_i = timeline[ t_ind ];
    }
}
 

/***************************************************************/
/**************** INLINED TRIDAG VERSION ***********************/
/***************************************************************/

inline void tridag_inline_local (
        /*** Read-Write Temporary Arrays ***/
//        __global   REAL*        a,
        REAL                    a_elem,
//        __global   REAL*        c,
        REAL                    c_elem,
//        __global   REAL*         d,
        REAL                    d_elem,
        __global   REAL*        y,
//        __global   REAL4*       scan_tmp4,
        REAL4                   data4,
        const unsigned int      SIZE,
/*** Temporary Array in Shared Space ***/   
        __local    REAL*        cache_tmp,   // size: [8*WORKSIZE*sizeof(REAL)]
        const size_t            first_id,
        REAL                    b0,
        REAL                    d0            
) { 
    //const size_t first_id = get_global_id(0) - (get_global_id(0) & (SIZE-1));
    //REAL  b0;

    { // SCAN with matrix multiplication
        __local REAL4*  l_data = (__local REAL4*)cache_tmp;

        //Calculate exclusive scan
        data4  = scanMatMultInclLocal( data4, l_data, SIZE );
        barrier(CLK_LOCAL_MEM_FENCE); // IMPORTANT!
        
        REAL  data = map_matmult( data4, b0 );
        cache_tmp[get_local_id(0)] = data;
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    { // Forward scan with linear function composition
        REAL2 data2;
        // prepare for scan
        if(get_global_id(0) == first_id) {
            //y[get_global_id(0)] = d0; //d[get_global_id(0)];  // y[ind] = d[ind];
            data2 = (REAL2) (0.0, 1.0);
        } else {
            data2 = (REAL2) ( d_elem, (0.0 - (a_elem/cache_tmp[get_local_id(0)-1])) );
        }

        {//Calculate exclusive scan
            __local    REAL2*  l_Data = ((__local REAL2*)cache_tmp) + get_local_size(0);
            data2  = scanLinFunCompInclLocal( data2, l_Data, SIZE );
        }

        //POST FWD MAP: y[ind] = fun.x + y0 * fun.y;
        //y[get_global_id(0)] = data2.x + d0 * data2.y;
        {
            cache_tmp[get_local_size(0) + get_local_id(0)] = data2.x + d0 * data2.y;
            barrier(CLK_LOCAL_MEM_FENCE);
        }
    }

    { // Backward scan with linear function composition
        const unsigned int inv_ind = (SIZE-1) - (get_global_id(0) & (SIZE-1));
        REAL2 data2;

        unsigned int local_offset = get_local_id(0) - (get_local_id(0) & (SIZE-1));
        b0 = cache_tmp[local_offset + get_local_size(0) + (SIZE-1)]/cache_tmp[local_offset + (SIZE-1)];

        // prepare for scan!
        if(get_global_id(0) == first_id) {
            data2 = (REAL2) (0.0, 1.0);
        } else {
            REAL myu = cache_tmp[local_offset + inv_ind];
            data2 = (REAL2) (   cache_tmp[local_offset + inv_ind + get_local_size(0)] / myu, 
                                0.0 - c_elem / myu );
                                //0.0 - c[first_id + inv_ind]/myu );
        }

        {// BACKWARD SCAN with linear function compsition
            __local    REAL2*  l_Data = ((__local REAL2*)cache_tmp) + get_local_size(0);
            data2  = scanLinFunCompInclLocal( data2, l_Data, SIZE );
        }

        // //POST BWD MAP: y[ind] = fun.x + y0 * fun.y;
        y[first_id + inv_ind] = data2.x + b0 * data2.y;
    }

}

__kernel __attribute__((reqd_work_group_size(WORKGROUP_SIZE, 1, 1)))
void tridag_inlined (
        /*** Read-Write Temporary Arrays ***/
        __global   REAL*        a,
        __global   REAL*        c,
        __global   REAL*        d,
        __global   REAL*        y,
        __global   REAL*        u,
        __global   REAL4*       scan_tmp4,
        const unsigned int      SIZE,
/*** Temporary Array in Shared Space ***/   
        __local    REAL*        cache_tmp   // size: [8*WORKSIZE*sizeof(REAL)]
) {
    const size_t first_id = get_global_id(0) - (get_global_id(0) & (SIZE-1));
    //const unsigned int inv_ind = (SIZE-1) - (get_global_id(0) & (SIZE-1));

    tridag_inline_local ( 
        a[get_global_id(0)], 
        c[first_id + (SIZE-1) - (get_global_id(0) & (SIZE-1))], 
            //(get_global_id(0) != first_id) ? c[first_id + (SIZE-1) - (get_global_id(0) & (SIZE-1))] : 0.0, //c
        d[get_global_id(0)], 
        y, 
        scan_tmp4[get_global_id(0)], 
        SIZE, 
        cache_tmp, 
        first_id, 
        scan_tmp4[first_id].x, 
        d[first_id] 
    );
}

  

/***************************************************************/
/************* TWO KERNELS VERSION: X & Y **********************/
/***************************************************************/

__kernel __attribute__((reqd_work_group_size(WORKGROUP_SIZE, 1, 1)))
void nordea_kernel_x (
        /*** Read-Only Scalars ***/
        __constant RWScalars*   ro_scals,   // __constant
        /*** Read-Only Arrays  ***/
        __global REAL*         myX,
        __global REAL3*        myDx,
        __global REAL3*        myDxx,
        __constant REAL*       myY,
        __constant REAL*       myDy,
        __constant REAL*       myDyy,
        /*** Read-Write Temporary Arrays ***/
        __global   REAL*       u,
        __global   REAL*       v,
        /*** Temporary Array in Shared Space ***/   
        __local    REAL*       cache_tmp,
        /*** The Result Array ***/
        __global   REAL*       res_arr        
) { 
   
    REAL3 myD_elem;
    
    unsigned int ind_x = get_global_id(0) & (ro_scals->NUM_X-1); 
    unsigned int ind_y = (get_global_id(0) & (ro_scals->NUM_XY-1)) / ro_scals->NUM_X;

    REAL tmp = 0.0;
    
    { 
        // I. second loop
        REAL3 myres_elem;
        REAL  cur_myMu  = 0.0;
        REAL  cur_myVar = ro_scals->nu; 
        cur_myVar *= cur_myVar;

        myres_elem = (REAL3)( myDyy[ind_y<<2], myDyy[(ind_y<<2)+1], myDyy[(ind_y<<2)+2] );
        myD_elem   = (REAL3)( myDy [ind_y<<2], myDy [(ind_y<<2)+1], myDy [(ind_y<<2)+2] );

        myD_elem   = cur_myMu*myD_elem + 0.5f*cur_myVar*myres_elem; 

        myres_elem = (REAL3)(   (ind_y != 0) ? res_arr[get_global_id(0) - ro_scals->NUM_X] : 0.0, 
                                 res_arr[get_global_id(0)], 
                                 (ind_y != ro_scals->NUM_Y-1) ? res_arr[get_global_id(0) + ro_scals->NUM_X] : 0.0
                            ); 

        myD_elem *= myres_elem;
        tmp += myD_elem.x + myD_elem.y + myD_elem.z;

#if TRANSPOSE_UV
        v[ get_global_id(0) ] = tmp;
#else
        v[ (get_global_id(0)/ro_scals->NUM_XY)*ro_scals->NUM_XY + ind_x*ro_scals->NUM_Y + ind_y] = tmp;
#endif

        // II. first loop 
        cur_myMu  = 0.0; // X
        cur_myVar = exp( 2 * ( ro_scals->beta*log(myX[ind_x]) + 
                                 myY[ind_y] - 0.5f*cur_myVar*ro_scals->timeline_i  // cur_myVar == nu*nu
                             )
                       ); 
        // CACHING ASSUMES ALL ELEMENTS IN X DIMENSION fit in the LOCALGROUP!
        cache_tmp[get_local_id(0)] = res_arr[ get_global_id(0) ];

        myD_elem  = myDx[ ind_x ];    

        barrier(CLK_LOCAL_MEM_FENCE);   


        myres_elem = (REAL3)(     (ind_x!=0) ? cache_tmp[ get_local_id(0) - 1 ] : 0.0, 
                                  cache_tmp[ get_local_id(0) ], 
                                  (ind_x != ro_scals->NUM_X-1) ? cache_tmp[ get_local_id(0) + 1 ] : 0.0
                            );

        tmp += ro_scals->dtInv*myres_elem.y;
        myD_elem = 0.5f*( cur_myMu*myD_elem + 0.5f*cur_myVar*myDxx[ ind_x ] );

        // write in cache the values of c!
        cache_tmp[get_local_size(0) + get_local_id(0)] = - myD_elem.z;      // c[get_global_id(0)] = - myD_elem.z;    

        barrier(CLK_LOCAL_MEM_FENCE);

        myres_elem = myD_elem*myres_elem;
        tmp += myres_elem.x + myres_elem.y + myres_elem.z;

        //a[get_global_id(0)] = - myD_elem.x; // b[ind] = ro_scals->dtInv - myDy_elem.y;
        myD_elem.x = - myD_elem.x; // holds `a[glb_ind]'
    }    


    { 
        // prepare for and call TRIDAG
        REAL4 scan_elem;
        //u[get_global_id(0)] = tmp; //u[j*NUM_X + i] = tmp1 + tmp2; //u[j][i] = tmp1 + tmp2;
        myD_elem.z = tmp;
    
        tmp = ro_scals->dtInv - myD_elem.y; // b[ind] = tmp;

        //barrier(CLK_LOCAL_MEM_FENCE);
        

        // third loop 
        if( ind_x > 0 ) { 
            myD_elem.y = cache_tmp[get_local_size(0) + get_local_id(0) - 1]; // i.e., c[glob_ind-1]
            scan_elem  = (REAL4) ( tmp, 0.0 - myD_elem.x * myD_elem.y, 1.0, 0.0 );
            myD_elem.y = cache_tmp[get_local_size(0) + get_local_id(0) + (ro_scals->NUM_X-1) - 2*ind_x]; // i.e., c[glb_ind_inv]
        } else { // i == get_global_id(0) == 0
            scan_elem = (REAL4) ( tmp, 0.0, 0.0, tmp ); //(REAL4) ( 1.0, 0.0, 0.0, 1.0 );
            myD_elem.y = cache_tmp[get_local_size(0) + get_local_id(0) + (ro_scals->NUM_X-1) ];

            cache_tmp[get_local_id(0)    ] = tmp;
            cache_tmp[get_local_id(0) + 1] = myD_elem.z;
        }

        barrier(CLK_LOCAL_MEM_FENCE);

        tmp     = cache_tmp[get_local_id(0) - ind_x    ];
        REAL d0 = cache_tmp[get_local_id(0) - ind_x + 1];

        barrier(CLK_LOCAL_MEM_FENCE);

        // CALL TRIDAG
        tridag_inline_local ( 
            myD_elem.x,                // a[glb_ind]
            myD_elem.y,                // c[glb_ind_inv]
            myD_elem.z,                // d[glb_ind] == u[glb_ind]
            u,                         // y
            scan_elem,                 // scan_tmp[glb_ind]
            ro_scals->NUM_X,           // SIZE
            cache_tmp,                 // cache
            get_global_id(0) - ind_x,  // first_id 
            tmp,                       // b0
            d0                         // d0
        );
    }
} 

////////////////////////////////////////////////////////////////
////////// AND KERNEL FOR THE Y-INNERMOST DIMENSION ////////////
////////////////////////////////////////////////////////////////

__kernel __attribute__((reqd_work_group_size(WORKGROUP_SIZE, 1, 1)))
void nordea_kernel_y (
        /*** Read-Only Scalars ***/
        __constant RWScalars*   ro_scals,  // __constant
        /*** Read-Only Arrays  ***/
        __global REAL3*        myDy,
        __global REAL3*        myDyy,
        /*** Read-Write Temporary Arrays ***/
        __global   REAL*        u,
        __global   REAL*        v,
        /*** Temporary Array in Shared Space ***/   
        __local    REAL*        cache_tmp    // SIZE: 8 * WORKGROUP_SIZE * sizeof(REAL)
) {

    REAL3 myDy_elem; 
    
    unsigned int ind_y = get_global_id(0) & (ro_scals->NUM_Y-1); 
    unsigned int ind_x = (get_global_id(0) & (ro_scals->NUM_XY-1)) / ro_scals->NUM_Y;

    {
        REAL dt_inv = ro_scals->dtInv;
#if TRANSPOSE_UV
        REAL my_u   = u[get_global_id(0)];
#else
        REAL my_u   = u[    (get_global_id(0)/ro_scals->NUM_XY)*ro_scals->NUM_XY + 
                            ind_y * ro_scals->NUM_X + ind_x ];
#endif

        REAL tmp = 0.0;
        REAL cur_myVarY = ro_scals->nu; cur_myVarY *= cur_myVarY*0.5f; //nu*nu;
        
        myDy_elem  = 0.0f - 0.5f*( tmp*myDy[ind_y] + cur_myVarY*myDyy[ind_y] );
 
        // a[ind] = myDy_elem.x; c[ind] = myDy_elem.z;   
        cache_tmp[get_local_size(0) + get_local_id(0)] = myDy_elem.z;  // c[ind]
        myDy_elem.y = dt_inv + myDy_elem.y; // b[ind]

        
        //v[get_global_id(0)] = dt_inv * my_u - 0.5 * v[get_global_id(0)];
        myDy_elem.z = dt_inv * my_u - 0.5f * v[get_global_id(0)]; // i.e., v[ind]

        if(ind_y == 0) {
            cache_tmp[get_local_id(0)    ] = myDy_elem.y; // b0
            cache_tmp[get_local_id(0) + 1] = myDy_elem.z; // d0, i.e., v[0]
        }

    }

    barrier(CLK_LOCAL_MEM_FENCE);

    {
        REAL4 scan_elem; 
        REAL b0, d0;

        // third loop 
        if( ind_y > 0 ) { 
            scan_elem = (REAL4) ( myDy_elem.y, 0.0 - myDy_elem.x * cache_tmp[get_local_size(0)+get_local_id(0)-1], 1.0, 0.0 );
        } else { 
            scan_elem = (REAL4) ( 1.0, 0.0, 0.0, 1.0 );
        }

        myDy_elem.y = cache_tmp[get_local_size(0) + get_local_id(0) + (ro_scals->NUM_Y-1) - 2*ind_y]; // i.e., c[ind_inv]

        b0 = cache_tmp[get_local_id(0) - ind_y    ];
        d0 = cache_tmp[get_local_id(0) - ind_y + 1]; // d0

        barrier(CLK_LOCAL_MEM_FENCE);

        // CALL TRIDAG
        tridag_inline_local ( 
            myDy_elem.x,               // a[glb_ind]
            myDy_elem.y,               // c[glb_ind_inv]
            myDy_elem.z,               // d[glb_ind] == v[glb_ind]
            v,                         // y
            scan_elem,                 // scan_tmp[glb_ind]
            ro_scals->NUM_Y,           // SIZE
            cache_tmp,                 // cache
            get_global_id(0) - ind_y,  // first_id 
            b0,                        // b0
            d0                         // d0
        );
    }
}
