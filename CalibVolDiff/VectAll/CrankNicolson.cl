#include "Constants.h"


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
/*********** PREPARE FOR TRIDAG X *****************************************/
/**************************************************************************/

#if 0

__kernel void prepare_tridag_x (
        /*** Read-Only Scalars ***/
        __global RWScalars*   ro_scals,   // __constant
        /*** Read-Only Arrays  ***/
        __constant REAL*        myX,
        __constant REAL*        myDx,
        __constant REAL*        myDxx,
        __constant REAL*        myY,
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
    REAL cur_myVar = ro_scals->nu; 
    REAL cur_myMu  = 0.0;

    const unsigned int glob_ind = get_global_id(2) * get_global_size(1) * get_global_size(0);
    const unsigned int ind = glob_ind + get_global_id(1)*get_global_size(0) + get_global_id(0);

    unsigned int mul3;
    REAL tmp;

    cur_myVar *= cur_myVar;

    // second loop    
    mul3 = REAL3_CT * get_global_id(1);
    tmp = 0.0;
    
    if(get_global_id(1)!=0)
        tmp += (cur_myMu*myDy[mul3] + 0.5*cur_myVar*myDyy[mul3])*res_arr[ind - get_global_size(0)]; //[(j-1)*NUM_X+i];

        tmp += (cur_myMu*myDy[mul3+1] + 0.5*cur_myVar*myDyy[mul3+1])*res_arr[ind]; // [j*NUM_X+i];

    if(get_global_id(1) != get_global_size(1)-1) //ro_scals->NUM_Y-1)
        tmp += (cur_myMu*myDy[mul3+2] + 0.5*cur_myVar*myDyy[mul3+2])*res_arr[ind + get_global_size(0)]; //[(j+1)*NUM_X+i];

    v[glob_ind + get_global_id(0)*get_global_size(1) + get_global_id(1)] = tmp;   // v[i*NUM_Y + j] = tmp;     //v[i][j] = tmp2;

    // first loop 
    cur_myMu  = 0.0; // X
    cur_myVar = exp( 2 * ( ro_scals->beta*log(myX[get_global_id(0)]) + 
                                 myY[get_global_id(1)] - 0.5*cur_myVar*ro_scals->timeline_i  // cur_myVar == nu*nu
                          )
                    ); 

    mul3 = REAL3_CT * get_global_id(0);
    tmp += ro_scals->dtInv*res_arr[ind]; // [j*NUM_X+i];

    if(get_global_id(0) != 0) //if(i!=0)
        tmp += 0.5*(cur_myMu*myDx[mul3] + 0.5*cur_myVar*myDxx[mul3])*res_arr[ind-1]; //[j*NUM_X+i-1];

        tmp += 0.5*(cur_myMu*myDx[mul3+1] + 0.5*cur_myVar*myDxx[mul3+1])*res_arr[ind]; //[j*NUM_X+i];

    if(get_global_id(0) != get_global_size(0) - 1) //if(i!=NUM_X-1)
        tmp += 0.5*(cur_myMu*myDx[mul3+2] + 0.5*cur_myVar*myDxx[mul3+2])*res_arr[ind+1]; //[j*NUM_X+i+1];

    u[ind] = tmp; //u[j*NUM_X + i] = tmp1 + tmp2; //u[j][i] = tmp1 + tmp2;
 
 
    // third loop (computing the implicit x parameters)
    a[ind] =                  - 0.5*(cur_myMu*myDx[mul3  ] + 0.5*cur_myVar*myDxx[mul3  ]); // a[j*NUM_X+i] = ...
    b[ind] =  ro_scals->dtInv - 0.5*(cur_myMu*myDx[mul3+1] + 0.5*cur_myVar*myDxx[mul3+1]); // b[j*NUM_X+i] = ...
    c[ind] =                  - 0.5*(cur_myMu*myDx[mul3+2] + 0.5*cur_myVar*myDxx[mul3+2]); // c[j*NUM_X+i] = ...

    barrier(CLK_GLOBAL_MEM_FENCE);

    // third loop 
    if( get_global_id(0) > 0 ) { 
        scan_tmp[ind] = (REAL4) ( b[ind], 0.0 - a[ind]*c[ind-1], 1.0, 0.0 );
    } else { // i == get_global_id(0) == 0
        y[ind] = b [ind];
        scan_tmp[ind] = (REAL4) ( 1.0, 0.0, 0.0, 1.0 );
    }
}

#else

__kernel void prepare_tridag_x (
        /*** Read-Only Scalars ***/
        __global RWScalars*   ro_scals,   // __constant
        /*** Read-Only Arrays  ***/
        __global REAL*         myX,
        __global REAL3*        myDx,
        __global REAL3*        myDxx,
        __constant REAL*         myY,
        __constant REAL3*        myDy,
        __constant REAL3*        myDyy,
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

    const unsigned int glob_ind = get_global_id(2) * get_global_size(1) * get_global_size(0);
    const unsigned int ind = glob_ind + get_global_id(1)*get_global_size(0) + get_global_id(0);

    REAL tmp = 0.0;

    cur_myVar *= cur_myVar;

    // second loop    
    {
        myres_elem = (REAL3)(   (get_global_id(1)!=0) ? res_arr[ind - get_global_size(0)] : 0.0, 
                                res_arr[ind], 
                                (get_global_id(1) != get_global_size(1)-1) ? res_arr[ind + get_global_size(0)] : 0.0
                            );
 
        myDy_elem  = myDy[ get_global_id(1) ]; 
        myDyy_elem = myDyy[ get_global_id(1) ];

        myDy_elem = cur_myMu*myDy_elem + 0.5*cur_myVar*myDyy_elem;
        myDy_elem *= myres_elem;
        tmp += myDy_elem.x + myDy_elem.y + myDy_elem.z;

        v[glob_ind + get_global_id(0)*get_global_size(1) + get_global_id(1)] = tmp;
        //v[ind] = tmp;
    }

    // first loop 
    cur_myMu  = 0.0; // X
    cur_myVar = exp( 2 * ( ro_scals->beta*log(myX[get_global_id(0)]) + 
                                 myY[get_global_id(1)] - 0.5*cur_myVar*ro_scals->timeline_i  // cur_myVar == nu*nu
                          )
                    ); 
    {
        // CACHING ASSUMES ALL ELEMENTS IN X DIMENSION fit in the LOCALGROUP!
        cache_tmp[get_local_id(0)] = res_arr[ ind ];
        
        myDy_elem  = myDx[ get_global_id(0) ]; 
        myDyy_elem = myDxx[ get_global_id(0) ];

        barrier(CLK_LOCAL_MEM_FENCE);


        myres_elem = (REAL3)(     (get_global_id(0)!=0) ? cache_tmp[ get_local_id(0) - 1 ] : 0.0, 
                                  cache_tmp[ get_local_id(0) ], 
                                  (get_global_id(0) != get_global_size(0)-1) ? cache_tmp[ get_local_id(0) + 1 ] : 0.0
                            );

        tmp += ro_scals->dtInv*myres_elem.y;
        myDy_elem = 0.5*( cur_myMu*myDy_elem + 0.5*cur_myVar*myDyy_elem );

        a[ind] = - myDy_elem.x; c[ind] = - myDy_elem.z;   // b[ind] = ro_scals->dtInv - myDy_elem.y;

        myDyy_elem = myDy_elem*myres_elem;
        tmp += myDyy_elem.x + myDyy_elem.y + myDyy_elem.z;
    }    

    u[ind] = tmp; //u[j*NUM_X + i] = tmp1 + tmp2; //u[j][i] = tmp1 + tmp2;

    barrier(CLK_GLOBAL_MEM_FENCE);

    tmp = ro_scals->dtInv - myDy_elem.y; // b[ind] = tmp;

    // third loop 
    if( get_global_id(0) > 0 ) { 
#if 0 
        cur_myVar = exp( 2 * ( ro_scals->beta*log(myX[get_global_id(0)-1]) + 
                                 myY[get_global_id(1)] - 0.5*ro_scals->nu*ro_scals->nu*ro_scals->timeline_i  // cur_myVar == nu*nu
                          )
                    ); 
        tmp = 0.0 - 0.5*(cur_myMu*myDx[get_global_id(0)-1].z + 0.5*cur_myVar*myDxx[get_global_id(0)-1].z); // c[j*NUM_X+i] = ...
        scan_tmp[ind] = (REAL4) ( b[ind], 0.0 - a[ind]*tmp, 1.0, 0.0 );
        c[ind-1] = tmp;
#else
        scan_tmp[ind] = (REAL4) ( tmp, myDy_elem.x * c[ind-1], 1.0, 0.0 );
        //scan_tmp[ind] = (REAL4) ( tmp, 0.0 - a[ind]*c[ind-1], 1.0, 0.0 );
#endif
    } else { // i == get_global_id(0) == 0
        y[ind] = tmp;    // = b [ind]
        scan_tmp[ind] = (REAL4) ( 1.0, 0.0, 0.0, 1.0 );
        //c[ind+get_global_size(0)-1] = 0.0 - 0.5*(cur_myMu*myDx[3*(get_global_size(0)-1)+2] + 0.5*cur_myVar*myDxx[3*(get_global_size(0)-1)+2]);
    }
} 

#endif // versions for prepare_for_tridag_x

/**************************************************************************/
/*********** IN-BETWEEN TRIDAGS -- i.e., prepare for TRIDAG Y *************/
/**************************************************************************/

__kernel void prepare_tridag_y (
        /*** Read-Only Scalars ***/
        __global RWScalars*   ro_scals,  // __constant
        /*** Read-Only Arrays  ***/
        __constant REAL*        myDy,
        __constant REAL*        myDyy,
        /*** Read-Write Temporary Arrays ***/
        __global   REAL*        a,
        __global   REAL*        b,
        __global   REAL*        c,
        __global   REAL*        u,
        __global   REAL*        v,
        __global   REAL*        y,
        __global   REAL4*       scan_tmp
        /*** Temporary Array in Shared Space ***/   
        //__local    REAL*        cache_tmp
) {
 
    // !!!0 <-> j; 1<-> i; 2 <-> k!!!
    const unsigned int ind      = get_global_id(2) * get_global_size(1) * get_global_size(0)
                                    + get_global_id(1) * get_global_size(0) + get_global_id(0); 
                                    // + get_global_id(1)*NUM_Y + get_global_id(0); 
    unsigned int mul3j = REAL3_CT * get_global_id(0);
    REAL cur_myMuY = 0.0;
    REAL cur_myVarY = ro_scals->nu; cur_myVarY *= cur_myVarY; //nu*nu;
    REAL dt_inv = ro_scals->dtInv;

    a[ind] =         - 0.5*(cur_myMuY*myDy[mul3j  ] + 0.5*cur_myVarY*myDyy[mul3j  ]); // a[i*NUM_Y+j] = ...
    b[ind] =  dt_inv - 0.5*(cur_myMuY*myDy[mul3j+1] + 0.5*cur_myVarY*myDyy[mul3j+1]); // b[i*NUM_Y+j] = ...
    //c[ind] =         - 0.5*(cur_myMuY*myDy[mul3j+2] + 0.5*cur_myVarY*myDyy[mul3j+2]); // c[i*NUM_Y+j] = ...

    // NUM_X = get_global_size(1)
    v[ind] =    dt_inv * 
                u[  get_global_id(2)* get_global_size(1) * get_global_size(0) + 
                    get_global_id(0)*get_global_size(1) + get_global_id(1)     ] 
                - 0.5*v[ind];  // y[j] = dtInv*u[j][i] - 0.5*v[i][j];

    //barrier(CLK_GLOBAL_MEM_FENCE);

    // third loop 
    if( get_global_id(0) > 0 ) { 
        const REAL elem = - 0.5*(cur_myMuY*myDy[mul3j-(REAL3_CT-2)] + 0.5*cur_myVarY*myDyy[mul3j-(REAL3_CT-2)]);
        scan_tmp[ind] = (REAL4) ( b[ind], 0.0 - a[ind]*elem, 1.0, 0.0 );   
        c[ind-1] = elem;

        //scan_tmp[ind] = (REAL4) ( b[ind], 0.0 - a[ind]*c[ind-1], 1.0, 0.0 );
    } else { // j == get_global_id(0) == 0
        y[ind] = b[ind];
        scan_tmp[ind] = (REAL4) ( 1.0, 0.0, 0.0, 1.0 );
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
#if 1
inline REAL4 matmult2(REAL4 a, REAL4 b) {
    REAL ax = a.x, az = a.z;

     
    REAL val = 1.0/(ax*b.x);
    //REAL val = (ax*b.x + a.y*b.z); if(val>=1.0) val = 1.0/val; else val = 1.0;

    a.x = (ax*b.x + a.y*b.z)*val;
    a.y = (ax*b.y + a.y*b.w)*val;
    a.z = (az*b.x + a.w*b.z)*val;
    a.w = (az*b.y + a.w*b.w)*val;
 
    //ax = 1.0/(ax*b.x);
    //b = (REAL4)(ax,ax,ax,ax);
    //a *= ax;
    return a;
}
#else
inline REAL4 matmult2(REAL4 a, REAL4 b) {
    REAL bx = b.x, bz = b.z, val = 1.0/(bx*a.x);

    b.x = (bx*a.x + b.y*a.z)*val;
    b.y = (bx*a.y + b.y*a.w)*val;
    b.z = (bz*a.x + b.w*a.z)*val;
    b.w = (bz*a.y + b.w*a.w)*val;

    //bx = 1.0/(bx*a.x);
    //b *= ((REAL4)(ax,ax,ax,ax));
    return b;
}

#endif


#define LOG2_WARP_SIZE 5U
#define WARP_SIZE (1U << LOG2_WARP_SIZE)

inline REAL4 warpScanInclMatMult(REAL4 idata, volatile __local REAL4 *l_Data, uint size){
    uint pos = 2 * get_local_id(0) - (get_local_id(0) & (size - 1));
    l_Data[pos] = (REAL4)(1.0, 0.0, 0.0, 1.0); 
    pos += size;
    l_Data[pos] = idata;

    if(size >=  2) l_Data[pos] = matmult2( l_Data[pos], l_Data[pos-1] );  //l_Data[pos] += l_Data[pos -  1];
    if(size >=  4) l_Data[pos] = matmult2( l_Data[pos], l_Data[pos-2] );  //l_Data[pos] += l_Data[pos -  2];
    if(size >=  8) l_Data[pos] = matmult2( l_Data[pos], l_Data[pos-4] );   //l_Data[pos] += l_Data[pos -  4];
    if(size >= 16) l_Data[pos] = matmult2( l_Data[pos], l_Data[pos-8] );   //l_Data[pos] += l_Data[pos -  8];
    if(size >= 32) l_Data[pos] = matmult2( l_Data[pos], l_Data[pos-16]);  //l_Data[pos] += l_Data[pos - 16];

    return l_Data[pos];
}


//Vector scan: the array to be scanned is stored
//in work-item private memory as uint4
inline REAL4 scanMatMultInclLocal(REAL4 idata4, __local REAL4 *l_Data, uint size){
    if(size > WARP_SIZE){
        //Bottom-level inclusive warp scan
        REAL4 warpResult = warpScanInclMatMult(idata4, l_Data, WARP_SIZE);

        //Save top elements of each warp for exclusive warp scan
        //sync to wait for warp scans to complete (because l_Data is being overwritten)
        barrier(CLK_LOCAL_MEM_FENCE);
        if( (get_local_id(0) & (WARP_SIZE - 1)) == (WARP_SIZE - 1) )
            l_Data[get_local_id(0) >> LOG2_WARP_SIZE] = warpResult;

        //wait for warp scans to complete
        barrier(CLK_LOCAL_MEM_FENCE); 
     
        if( get_local_id(0) < (WORKGROUP_SIZE / WARP_SIZE) ){
            //grab top warp elements
            REAL4 val = l_Data[get_local_id(0)];
            //calculate inclusive scan and write back to shared memory
            l_Data[get_local_id(0)] = warpScanInclMatMult(val, l_Data, size >> LOG2_WARP_SIZE);
        }

        //return updated warp scans with exclusive scan results
        barrier(CLK_LOCAL_MEM_FENCE);

        //if( (get_local_id(0) >> LOG2_WARP_SIZE) > 0 )
        if( ( (get_local_id(0) & (size-1)) >> LOG2_WARP_SIZE ) > 0 )
            warpResult = matmult2( warpResult, l_Data[ (get_local_id(0) >> LOG2_WARP_SIZE) - 1] );

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
#if 0
    REAL4 idata4           = data[get_global_id(0)]; //(REAL4)(1.0+get_global_id(0), 0.0, 1.0, 0);
    REAL4 odata4           = warpScanInclMatMult(idata4, l_Data, size); //scanMatMultInclLocal( idata4, l_Data, size );
    data[get_global_id(0)] = odata4;
#endif
#if 1
    //Load data
    REAL4 idata4 = data[get_global_id(0)];

    //Calculate exclusive scan
    REAL4 odata4  = scanMatMultInclLocal( idata4, l_Data, size );

    //Write back
    data[get_global_id(0)] = odata4; //(REAL4)(33.33,33.33,33.33,33.33);//odata4;
#endif
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

inline REAL2 warpScanInclLinFunComp(REAL2 idata, volatile __local REAL2 *l_Data, uint size){
    uint pos = 2 * get_local_id(0) - (get_local_id(0) & (size - 1));
    l_Data[pos] = (REAL2)(0.0, 1.0);
    pos += size;
    l_Data[pos] = idata;

    if(size >=  2) l_Data[pos] = linfuncomp( l_Data[pos], l_Data[pos-1] );  
    if(size >=  4) l_Data[pos] = linfuncomp( l_Data[pos], l_Data[pos-2] );  
    if(size >=  8) l_Data[pos] = linfuncomp( l_Data[pos], l_Data[pos-4] );   
    if(size >= 16) l_Data[pos] = linfuncomp( l_Data[pos], l_Data[pos-8] );
    if(size >= 32) l_Data[pos] = linfuncomp( l_Data[pos], l_Data[pos-16]);

    return l_Data[pos];
}

inline REAL2 scanLinFunCompInclLocal(REAL2 idata4, __local REAL2 *l_Data, uint size){
    if(size > WARP_SIZE){
        //Bottom-level inclusive warp scan
        REAL2 warpResult = warpScanInclLinFunComp(idata4, l_Data, WARP_SIZE);

        //Save top elements of each warp for exclusive warp scan
        //sync to wait for warp scans to complete (because l_Data is being overwritten)
        barrier(CLK_LOCAL_MEM_FENCE);
        if( (get_local_id(0) & (WARP_SIZE - 1)) == (WARP_SIZE - 1) )
            l_Data[get_local_id(0) >> LOG2_WARP_SIZE] = warpResult;

        //wait for warp scans to complete
        barrier(CLK_LOCAL_MEM_FENCE);

        if( get_local_id(0) < (WORKGROUP_SIZE / WARP_SIZE) ){
            //grab top warp elements
            REAL2 val = l_Data[get_local_id(0)];
            //calculate inclusive scan and write back to shared memory
            l_Data[get_local_id(0)] = warpScanInclLinFunComp(val, l_Data, size >> LOG2_WARP_SIZE);
        }

        //return updated warp scans with exclusive scan results
        barrier(CLK_LOCAL_MEM_FENCE);

        if( ( (get_local_id(0) & (size-1)) >> LOG2_WARP_SIZE ) > 0 )
            warpResult = linfuncomp( warpResult, l_Data[ (get_local_id(0) >> LOG2_WARP_SIZE) - 1] );

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

inline REAL map_matmult(REAL4 tmp, REAL val1, REAL val2) {

    REAL nom   = tmp.x*val1+tmp.y*val2;
    REAL denom = tmp.z*val1+tmp.w*val2;

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
        u[ind] = map_matmult( scan_tmp[ind], b0, 1.0 );
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
// and to avoid bank conflicts in shared memory.  This kernel is up to 11x faster
// than the naive kernel below.  Note that the shared memory array is sized to 
// (BLOCK_DIM+1)*BLOCK_DIM.  This pads each row of the 2D block in shared memory 
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
 

