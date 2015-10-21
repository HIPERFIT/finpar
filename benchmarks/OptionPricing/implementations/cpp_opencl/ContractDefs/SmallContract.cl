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

// Small Dataset Contract Definition
inline 
void payoffFunction( // SIMPLE
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
    real_t amount = ((underlyings(0,0)) - 4000.0) * md_detvals[model_num*num_pricers];
    amount = fmax((real_t)0.0, amount);
	trajectory_inner( num_cash_flows, model_num, 0, amount, md_discts, vhat );
}

