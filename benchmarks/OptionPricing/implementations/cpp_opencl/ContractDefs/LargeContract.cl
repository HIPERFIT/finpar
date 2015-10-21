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

// Large Dataset Contract Definition
inline 
void payoffFunction(                    // BARRIER
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
    bool goto_L40, x3309 = false; 
    real_t mult  = 0.0;
    for(int i=0; (i<367 && (!x3309)); i++) {
        x3309 = x3309 || 
                    (underlyings(i,0) <= 2630.6349999999998) ||
                    (underlyings(i,1) <= 8288.) ||
                    (underlyings(i,2) <= 840. ) ;
    }

    trajectory_inner( num_cash_flows, model_num, 0, 100., md_discts, vhat );

    goto_L40 = x3309 && 
                ( (underlyings(366,0) < 3758.05) ||
                  (underlyings(366,1) < 11840. ) ||
                  (underlyings(366,2) < 1200.  ) );

    mult = (goto_L40) ? 1.0 : 0.0;

    const real_t min_amount = fmin( (underlyings(366,1) / 11840.) - 1., 
                                  fmin( (underlyings(366,2) / 1200.  ) - 1., 
                                        (underlyings(366,0) / 3758.05) - 1. ) );

    const real_t amount = 1000. * ( 1. + mult*min_amount ) ;

    trajectory_inner( num_cash_flows, model_num, 1, amount, md_discts, vhat );
}

