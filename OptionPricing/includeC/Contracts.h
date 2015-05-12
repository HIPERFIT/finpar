#ifndef CONTRACTS_CODE
#define CONTRACTS_CODE

#define underlyings(i,j)   \
		(inst_traj[i*num_under + j])

inline
void trajectory_inner(
        const UINT&     num_cash_flows,
        const UINT&     model_ind,  // the model index
        const UINT&     disct_index,// index of the discount 
        const REAL&     amount,     // update with amount
        const REAL*     md_discts,  // [model_num][num_cash_flow] discounts
              double*   vhat        // [model_num] Accumulated per-model price
) {
    UINT ind = disct_index + model_ind*num_cash_flows;
    vhat[model_ind] += (amount * md_discts[ind]);
}

inline 
void trajectory_contract1( // SMALL
        const UINT&     model_num,  // the index of the current model
        const UINT&     num_under,  // the number of underlyings
        const UINT&     num_cash_flows,
        const UINT&     num_pricers,// the number of deterministic procers 
        const REAL*     md_discts,  // [num_models][num_cash_flows] discounts
        const REAL*     md_detvals, // [num_models, num_det_pricers]  pricers 
        const REAL*     inst_traj,  // [num_dates, num_under] current trajectory
              double*   vhat        // [model_num] Accumulated per-model price
        
) {
    REAL amount = ((underlyings(0,0)) - 4000.0) * md_detvals[model_num*num_pricers];
    amount = fmax(0.0, amount);
	trajectory_inner( num_cash_flows, model_num, 0, amount, md_discts, vhat );
}

inline 
void trajectory_contract2( // MEDIUM
        const UINT&     model_num,  // the index of the current model
        const UINT&     num_under,  // the number of underlyings
        const UINT&     num_cash_flows,
        const UINT&     num_pricers,// the number of deterministic procers 
        const REAL*     md_discts,  // [num_models][num_cash_flows] discounts
        const REAL*     md_detvals, // [num_models, num_det_pricers]  pricers 
        const REAL*     inst_traj,  // [num_dates, num_under] current trajectory
              double*   vhat        // [model_num] Accumulated per-model price
) {
	double x50;
    int    date;
    double amount;
 
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


inline 
void trajectory_contract3( // MARGE (Barrier Contract)
        const UINT&     model_num,  // the index of the current model
        const UINT&     num_under,  // the number of underlyings
        const UINT&     num_cash_flows,
        const UINT&     num_pricers,// the number of deterministic procers 
        const REAL*     md_discts,  // [num_models][num_scash_flow] discounts
        const REAL*     md_detvals, // [num_models, num_det_pricers]  pricers 
        const REAL*     inst_traj,  // [num_dates, num_under] current trajectory
              double*   vhat        // [model_num] Accumulated per-model price
) {
    bool goto_L40, x3309 = false; 
    REAL mult  = 0.0;
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

    const REAL min_amount = fmin( (underlyings(366,1) / 11840.) - 1., 
                                  fmin( (underlyings(366,2) / 1200.  ) - 1., 
                                        (underlyings(366,0) / 3758.05) - 1. ) );

    const REAL amount = 1000. * ( 1. + mult*min_amount ) ;

    trajectory_inner( num_cash_flows, model_num, 1, amount, md_discts, vhat );
}


inline
void aggregDiscountedPayoff(
        const UINT&     model_num,  // the index of the current model
        const UINT&     contr_num,
        const UINT&     num_under,  // the number of underlyings
        const UINT&     num_cash_flows,
        const UINT&     num_pricers,// the number of deterministic procers 
        const REAL*     md_discts,  // [num_models][num_scash_flow] discounts
        const REAL*     md_detvals, // [num_models, num_det_pricers]  pricers 
        const REAL*     inst_traj,  // [num_dates, num_under] current trajectory
              double*   vhat        // [model_num] Accumulated per-model price
) {
    if        (contr_num == 1) {
       	trajectory_contract1(model_num, num_under, num_cash_flows, num_pricers, md_discts, md_detvals, inst_traj, vhat);
    } else if (contr_num == 2) {
        trajectory_contract2(model_num, num_under, num_cash_flows, num_pricers, md_discts, md_detvals, inst_traj, vhat);
    } else if (contr_num == 3) {
        trajectory_contract3(model_num, num_under, num_cash_flows, num_pricers, md_discts, md_detvals, inst_traj, vhat);
    } else {
        trajectory_contract1(model_num, num_under, num_cash_flows, num_pricers, md_discts, md_detvals, inst_traj, vhat);
    }
}

#endif //CONTRACTS_CODE

