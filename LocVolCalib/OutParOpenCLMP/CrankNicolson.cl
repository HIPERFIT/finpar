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



/***********************************************/
/************ sequential TRIDAG   **************/
/***********************************************/

inline void tridag_seq(
    __global REAL*   a,
    __global REAL*   b,
    __global REAL*   c,
    __global REAL*   d,
    int              n,
    __global REAL*   y,
    __global REAL*   u)
{
    int    i, UB;
    REAL   beta;

    y[0] = d[0];
    u[0] = b[0];

    UB = n << lgWARP;
    for(i=WARP; i<UB; i+=WARP) {
        beta = a[i] / u[i-WARP];

        u[i] = b[i] - beta*c[i-WARP];
        y[i] = d[i] - beta*y[i-WARP];
    }

    UB -= WARP;
    y[UB] = y[UB] / u[UB];

    UB = (n-2) << lgWARP;
    for(i=UB; i>=0; i-=WARP) {
        y[i] = (y[i] - c[i]*y[i+WARP]) / u[i];
    }
}

/**************************************************************************/
/*********** PREPARE FOR TRIDAG X *****************************************/
/**************************************************************************/

__kernel void nordea_kernel_x (
        /*** Read-Only Scalars ***/
        __constant RWScalars*  ro_scals,
        /*** Read-Only Arrays  ***/
        __constant REAL*       myX,
        __constant REAL*       myDx,
        __constant REAL*       myDxx,
        __global   REAL*       myY,
        __global   REAL4*      myDy,
        __global   REAL4*      myDyy,
        /*** Read-Write Temporary Arrays ***/
        __global   REAL*        a,
        __global   REAL*        b,
        __global   REAL*        c,
        __global   REAL*        y,
        __global   REAL*        u,
        __global   REAL*        v,
        /*** The Result Array ***/
        __global   REAL*        res_arr        
) {
          unsigned int ind, i, ind_rev;
    //const unsigned int step_rev = ro_scals->NUM_Y; //(ro_scals->NUM_Y << LOG2_WARP_SIZE); 
    REAL  cur_myVar, cur_myMu, exp_add;
    REAL4 myDy_el;

    { // cur_myVarY = nu*nu; cur_myMuY = 0.0; compute ind_y and exp_add
        cur_myVar = ro_scals->nu;  cur_myVar *= cur_myVar; 
        cur_myMu  = 0.0;

        //ind_y = get_global_id(0) & (ro_scals->NUM_Y - 1);
        exp_add = myY[get_global_id(0) & (ro_scals->NUM_Y - 1)] - 0.5*cur_myVar*ro_scals->timeline_i;
    }

    { // adjust arrays: except for `res_arr' and `v' which are accesses as `v[i,j]', hence memory is already coalesced!
        unsigned int offset =  (get_global_id(0)>>lgWARP); 
        offset *= ( ro_scals->NUM_X << lgWARP );
        offset += ( get_global_id(0) & (WARP-1) );

        a       += offset;
        b       += offset;
        c       += offset;
        y       += offset;
        u       += offset;
    }

    { // computing myDy_el = (cur_myMuY*myDy[j] + 0.5*cur_myVarY*myDyy[j])
        i = get_global_id(0) & (ro_scals->NUM_Y - 1);
        myDy_el = cur_myMu*myDy[i] + (0.5f*cur_myVar)*myDyy[i];
        //i = REAL3_CT*(get_global_id(0) & (ro_scals->NUM_Y - 1));
        //myDy_el = (REAL4)(myDyy[i], myDyy[i+1], myDyy[i+2], 1.0);
        //myDy_el = myDy_el * (0.5*cur_myVar);
        //myDy_el += cur_myMu*((REAL4)(myDy[i], myDy[i+1], myDy[i+2], 1.0));
    }


    ind_rev  = (get_global_id(0)/ro_scals->NUM_Y)*ro_scals->NUM_XY + (get_global_id(0) & (ro_scals->NUM_Y - 1)); 
    for( i=0, ind=0; i<ro_scals->NUM_X; i++, ind+=WARP, ind_rev+=ro_scals->NUM_Y ) {
        myDy_el.w = 0.0;

        { // comput cur_myMy and cur_myVar
            cur_myMu  = 0.0; // X
            cur_myVar = exp( 2 * ( ro_scals->beta*log(myX[i]) + exp_add ) ); 
            //cur_myVar = exp(2*(ro_scals->beta*log(myX[i]) + myY[get_global_id(0) & (ro_scals->NUM_Y - 1)] - 0.5*ro_scals->nu*ro_scals->nu*ro_scals->timeline_i));
        }
        
        { // compute v[i, j]
            if( (ind_rev & (ro_scals->NUM_Y-1)) != 0 )                 myDy_el.w += myDy_el.x*res_arr[ind_rev-1]; //[(j-1)*NUM_X+i];
                                                                       myDy_el.w += myDy_el.y*res_arr[ind_rev  ]; // [j*NUM_X+i];
            if( (ind_rev & (ro_scals->NUM_Y-1)) != ro_scals->NUM_Y-1 ) myDy_el.w += myDy_el.z*res_arr[ind_rev+1]; //[(j+1)*NUM_X+i];

            //v[k*NUM_X*NUM_Y + i*NUM_Y + j] = tmp;   // v[i*NUM_Y + j] = tmp;     //v[i][j] = tmp2;

            v[ ind_rev ] = myDy_el.w;
        }

        { // compute u, a, b, c
            const unsigned int im3 = REAL3_CT*i;

            if(i!=0)
                myDy_el.w += 0.5*(cur_myMu*myDx[im3  ] + 0.5*cur_myVar*myDxx[im3  ])*res_arr[ind_rev-ro_scals->NUM_Y]; //[j*NUM_X+i-1];

            // tmp += ro_scals->dtInv * res_arr[ind]; WAS included in the next line!
            myDy_el.w += ( 0.5*(cur_myMu*myDx[im3+1] + 0.5*cur_myVar*myDxx[im3+1]) + ro_scals->dtInv)*res_arr[ind_rev  ]; //[j*NUM_X+i];

            if(i!=ro_scals->NUM_X-1)
                myDy_el.w += 0.5*(cur_myMu*myDx[im3+2] + 0.5*cur_myVar*myDxx[im3+2])*res_arr[ind_rev+ro_scals->NUM_Y]; //[j*NUM_X+i+1];

            u[ind] = myDy_el.w;

            // compute a,b,c
            a[ind] =              0.0 - 0.5*(cur_myMu*myDx[im3  ] + 0.5*cur_myVar*myDxx[im3  ]); // a[j*NUM_X+i] = ...
            b[ind] =  ro_scals->dtInv - 0.5*(cur_myMu*myDx[im3+1] + 0.5*cur_myVar*myDxx[im3+1]); // b[j*NUM_X+i] = ...
            c[ind] =              0.0 - 0.5*(cur_myMu*myDx[im3+2] + 0.5*cur_myVar*myDxx[im3+2]); // c[j*NUM_X+i] = ...
        }
    }

    tridag_seq(a, b, c, u, ro_scals->NUM_X, u, y);
} 


/**************************************************************************/
/*********** IN-BETWEEN TRIDAGS -- i.e., prepare for TRIDAG Y *************/
/**************************************************************************/

__kernel void nordea_kernel_y (
        /*** Read-Only Scalars ***/
        __constant RWScalars*   ro_scals,  // __constant
        /*** Read-Only Arrays  ***/
        __constant REAL*        myDy,
        __constant REAL*        myDyy,
        /*** Read-Write Temporary Arrays ***/
        __global   REAL*        a,
        __global   REAL*        b,
        __global   REAL*        c,
        __global   REAL*        u,
        __global   REAL*        v,
        __global   REAL*        y
) {
    //const unsigned int step_rev = ro_scals->NUM_X; //(ro_scals->NUM_X << LOG2_WARP_SIZE);
    unsigned int ind, ind_rev, i;
    REAL  cur_myVar, cur_myMu;

    { // cur_myVarY = nu*nu; cur_myMuY = 0.0; compute ind_y and exp_add
        cur_myVar = ro_scals->nu;  cur_myVar *= cur_myVar; 
        cur_myMu  = 0.0;
    }

    { // adjust arrays! except for `u' whose access pattern is `u[j,i]' which is already coalesced!
        unsigned int offset = (get_global_id(0)>>lgWARP); 
        offset *= ( ro_scals->NUM_Y << lgWARP );
        offset += ( get_global_id(0) & (WARP-1) );

        a       += offset;
        b       += offset;
        c       += offset;
        y       += offset;
        v       += offset;
    }

    ind_rev  = (get_global_id(0)/ro_scals->NUM_X)*ro_scals->NUM_XY + (get_global_id(0) & (ro_scals->NUM_X - 1)); // should step with `NUM_Y*32'
    for( i=0, ind=0; i<ro_scals->NUM_Y; i++, ind+=WARP, ind_rev+=ro_scals->NUM_X) {
        unsigned int im3 = REAL3_CT*i;

        a[ind] =              0.0 - 0.5*(cur_myMu*myDy[im3  ] + 0.5*cur_myVar*myDyy[im3  ]); // a[i*NUM_Y+j] = ...
        b[ind] =  ro_scals->dtInv - 0.5*(cur_myMu*myDy[im3+1] + 0.5*cur_myVar*myDyy[im3+1]); // b[i*NUM_Y+j] = ...
        c[ind] =              0.0 - 0.5*(cur_myMu*myDy[im3+2] + 0.5*cur_myVar*myDyy[im3+2]); // c[i*NUM_Y+j] = ...

        v[ind] =  ro_scals->dtInv * u[ind_rev] - 0.5 * v[ind];
                  //u[ (get_global_id(0)/ro_scals->NUM_XY)*ro_scals->NUM_XY + i*ro_scals->NUM_X + ind_x] -
    }

    tridag_seq(a, b, c, v, ro_scals->NUM_Y, v, y);
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

    // update time-loop-variant scalars
    if(get_global_id(2) + get_global_id(1) + get_global_id(0) == 0) {
        int t_ind =  --ro_scals->t_ind;
        ro_scals->dtInv      = 1.0 / ( timeline[t_ind+1] - timeline[t_ind] );
        ro_scals->timeline_i = timeline[ t_ind ];
    }
}
 


