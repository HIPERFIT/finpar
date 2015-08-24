#include "Constants.h"

#ifdef VECT_KERNEL
    #define underlyings(i,j)   \
        (inst_traj[(i*num_under + j)<<lgWARP])
#else
    #define underlyings(i,j)   \
        (inst_traj[i*num_under + j])
#endif

inline
void trajectory_inner(
        const       UINT    num_cash_flows, // number of discounts
        const       UINT    model_ind,      // the model index
        const       UINT    disct_index,    // index of the discount 
        const       real_t    amount,         // update with amount
        __constant real_t*    md_discts,  
        __local    real_t*    local_vhat
) {
    UINT index = model_ind * get_local_size(0); //block_size;
    local_vhat[index] += ( amount * md_discts[model_ind*num_cash_flows+disct_index] );
}

// Medium Dataset Contract Definition
inline 
void payoffFunction( 
        const      UINT     model_num,  // the index of the current model
        const      UINT     num_under,  // the number of underlyings
        const      UINT     num_cash_flows, // the number of discounts
        const      UINT     num_pricers,// the number of deterministic procers 
        __constant real_t*    md_discts,  // [num_models][num_cash_flows] discounts
        __constant real_t*    md_detvals, // [num_models, num_det_pricers]  pricers 
#ifdef VECT_KERNEL
        __global   real_t*    inst_traj,
#else
        const      real_t*    inst_traj,  // [num_dates, num_under] current trajectory
#endif
        __local    real_t*    vhat        // [model_num] Accumulated per-model price        
) {
	real_t x50;
    UINT date;
    real_t amount;
 
    if        ( 1. <= fmin( underlyings(0,1) / 11840., fmin( underlyings(0,2) / 1200., underlyings(0,0) / 3758.05 ) ) ) {
        //goto L21;
        date = 0; amount = 1150.;
    } else if ( 1. <= fmin( underlyings(1,1) / 11840., fmin( underlyings(1,2) / 1200., underlyings(1,0) / 3758.05 ) ) ) {
        //goto L22;
        date = 1; amount = 1300.;
    } else if ( 1. <= fmin( underlyings(2,1) / 11840., fmin( underlyings(2,2) / 1200., underlyings(2,0) / 3758.05 ) ) ) {
        //goto L23;
        date = 2; amount = 1450.;
    } else if ( 1. <= fmin( underlyings(3,1) / 11840., fmin( underlyings(3,2) / 1200., underlyings(3,0) / 3758.05 ) ) ) {
        //goto L24;
        date = 3; amount = 1600.;
    } else {
        date = 4; 
        x50 = fmin( underlyings(4,1) / 11840., fmin( underlyings(4,2) / 1200., underlyings(4,0) / 3758.05 ) );
        if( 1. <= x50 )        { amount = 1750.;     }
        else if ( 0.75 < x50 ) { amount = 1000.;     }
        else                   { amount = x50*1000.; }
    }

    trajectory_inner( num_cash_flows, model_num, date, amount, md_discts, vhat );
}

