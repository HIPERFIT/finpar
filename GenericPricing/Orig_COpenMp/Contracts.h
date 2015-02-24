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
#if 0
    UINT prod = ro_scal->num_contracts * ro_scal->inst_num;
    UINT index = model_num*ro_scal->num_contracts + contract_number;
    index += cur_iter*prod;
#endif
    UINT ind = disct_index + model_ind*num_cash_flows;
    vhat[model_ind] += (amount * md_discts[ind]);

    //printf("Updating vhat with: %f %f \n", amount, md_discts[ind] );

	//instance *inst = (instance*) model->model_data;
	//inst->vhat[contract_number] += amount * inst->discounts[date_index];
}

inline 
void trajectory_contract1( // SIMPLE
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
void trajectory_contract2( // WORST OFF
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

    if ((1. <= fmin((underlyings(0,1) / 11840.), fmin((underlyings(0,2) / 1200.), (underlyings(0,0) / 3758.05))))) goto L21;
    if ((1. <= fmin((underlyings(1,1) / 11840.), fmin((underlyings(1,2) / 1200.), (underlyings(1,0) / 3758.05))))) goto L22;
    if ((1. <= fmin((underlyings(2,1) / 11840.), fmin((underlyings(2,2) / 1200.), (underlyings(2,0) / 3758.05))))) goto L23;
    if ((1. <= fmin((underlyings(3,1) / 11840.), fmin((underlyings(3,2) / 1200.), (underlyings(3,0) / 3758.05))))) goto L24;
    x50=fmin((underlyings(4,1) / 11840.), fmin((underlyings(4,2) / 1200.), (underlyings(4,0) / 3758.05)));

    //model->notify_cash_flow(model, 0, 1000., 4 /*2017-02-03, 2017-01-27, EUR*/);
    trajectory_inner( num_cash_flows, model_num, 4, 1000., md_discts, vhat );

	if ((1. <= x50)) goto L25;
	if ((0.75 < x50)) return;

	//model->notify_cash_flow(model, 0, (-(1000. * (1. - x50))), 4 /*2017-02-03, 2017-01-27, EUR*/); return;
    trajectory_inner( num_cash_flows, model_num, 4, (-(1000. * (1. - x50))), md_discts, vhat ); return;
L25:
	//model->notify_cash_flow(model, 0, 750., 4 /*2017-02-03, 2017-01-27, EUR*/); return;
    trajectory_inner( num_cash_flows, model_num, 4, 750., md_discts, vhat ); return;
L24:
	//model->notify_cash_flow(model, 0, 1600., 3 /*2016-02-03, 2016-01-27, EUR*/); return;
    trajectory_inner( num_cash_flows, model_num, 3, 1600., md_discts, vhat ); return;
L23:
	//model->notify_cash_flow(model, 0, 1450., 2 /*2015-02-03, 2015-01-27, EUR*/); return;
    trajectory_inner( num_cash_flows, model_num, 2, 1450., md_discts, vhat ); return;
L22:
	//model->notify_cash_flow(model, 0, 1300., 1 /*2014-02-03, 2014-01-27, EUR*/); return;
    trajectory_inner( num_cash_flows, model_num, 1, 1300., md_discts, vhat ); return;
L21:
	//model->notify_cash_flow(model, 0, 1150., 0 /*2013-02-01, 2013-01-27, EUR*/); return;
    trajectory_inner( num_cash_flows, model_num, 0, 1150., md_discts, vhat ); return;
}




// LENGTHY CONTRACT
inline 
void trajectory_contract3( // BARRIER
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

#if 0
inline 
void trajectory_contract3_OLD( // BARRIER
        const UINT&     model_num,  // the index of the current model
        const UINT&     num_under,  // the number of underlyings
        const UINT&     num_cash_flows,
        const UINT&     num_pricers,// the number of deterministic procers 
        const REAL*     md_discts,  // [num_models][num_scash_flow] discounts
        const REAL*     md_detvals, // [num_models, num_det_pricers]  pricers 
        const REAL*     inst_traj,  // [num_dates, num_under] current trajectory
              double*   vhat        // [model_num] Accumulated per-model price
) {
	//assert(model_num==0 && "NUM MODELS > 1");
	UINT x3309;

	  if ((underlyings(0,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(0,2) <= 840.)) goto L50;
	  if ((underlyings(0,1) <= 8288.)) goto L50;
	  if ((underlyings(1,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(1,2) <= 840.)) goto L50;
	  if ((underlyings(1,1) <= 8288.)) goto L50;
	  if ((underlyings(2,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(2,2) <= 840.)) goto L50;
	  if ((underlyings(2,1) <= 8288.)) goto L50;
	  if ((underlyings(3,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(3,2) <= 840.)) goto L50;
	  if ((underlyings(3,1) <= 8288.)) goto L50;
	  if ((underlyings(4,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(4,2) <= 840.)) goto L50;
	  if ((underlyings(4,1) <= 8288.)) goto L50;
	  if ((underlyings(5,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(5,2) <= 840.)) goto L50;
	  if ((underlyings(5,1) <= 8288.)) goto L50;
	  if ((underlyings(6,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(6,2) <= 840.)) goto L50;
	  if ((underlyings(6,1) <= 8288.)) goto L50;
	  if ((underlyings(7,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(7,2) <= 840.)) goto L50;
	  if ((underlyings(7,1) <= 8288.)) goto L50;
	  if ((underlyings(8,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(8,2) <= 840.)) goto L50;
	  if ((underlyings(8,1) <= 8288.)) goto L50;
	  if ((underlyings(9,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(9,2) <= 840.)) goto L50;
	  if ((underlyings(9,1) <= 8288.)) goto L50;
	  if ((underlyings(10,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(10,2) <= 840.)) goto L50;
	  if ((underlyings(10,1) <= 8288.)) goto L50;
	  if ((underlyings(11,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(11,2) <= 840.)) goto L50;
	  if ((underlyings(11,1) <= 8288.)) goto L50;
	  if ((underlyings(12,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(12,2) <= 840.)) goto L50;
	  if ((underlyings(12,1) <= 8288.)) goto L50;
	  if ((underlyings(13,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(13,2) <= 840.)) goto L50;
	  if ((underlyings(13,1) <= 8288.)) goto L50;
	  if ((underlyings(14,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(14,2) <= 840.)) goto L50;
	  if ((underlyings(14,1) <= 8288.)) goto L50;
	  if ((underlyings(15,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(15,2) <= 840.)) goto L50;
	  if ((underlyings(15,1) <= 8288.)) goto L50;
	  if ((underlyings(16,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(16,2) <= 840.)) goto L50;
	  if ((underlyings(16,1) <= 8288.)) goto L50;
	  if ((underlyings(17,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(17,2) <= 840.)) goto L50;
	  if ((underlyings(17,1) <= 8288.)) goto L50;
	  if ((underlyings(18,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(18,2) <= 840.)) goto L50;
	  if ((underlyings(18,1) <= 8288.)) goto L50;
	  if ((underlyings(19,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(19,2) <= 840.)) goto L50;
	  if ((underlyings(19,1) <= 8288.)) goto L50;
	  if ((underlyings(20,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(20,2) <= 840.)) goto L50;
	  if ((underlyings(20,1) <= 8288.)) goto L50;
	  if ((underlyings(21,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(21,2) <= 840.)) goto L50;
	  if ((underlyings(21,1) <= 8288.)) goto L50;
	  if ((underlyings(22,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(22,2) <= 840.)) goto L50;
	  if ((underlyings(22,1) <= 8288.)) goto L50;
	  if ((underlyings(23,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(23,2) <= 840.)) goto L50;
	  if ((underlyings(23,1) <= 8288.)) goto L50;
	  if ((underlyings(24,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(24,2) <= 840.)) goto L50;
	  if ((underlyings(24,1) <= 8288.)) goto L50;
	  if ((underlyings(25,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(25,2) <= 840.)) goto L50;
	  if ((underlyings(25,1) <= 8288.)) goto L50;
	  if ((underlyings(26,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(26,2) <= 840.)) goto L50;
	  if ((underlyings(26,1) <= 8288.)) goto L50;
	  if ((underlyings(27,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(27,2) <= 840.)) goto L50;
	  if ((underlyings(27,1) <= 8288.)) goto L50;
	  if ((underlyings(28,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(28,2) <= 840.)) goto L50;
	  if ((underlyings(28,1) <= 8288.)) goto L50;
	  if ((underlyings(29,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(29,2) <= 840.)) goto L50;
	  if ((underlyings(29,1) <= 8288.)) goto L50;
	  if ((underlyings(30,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(30,2) <= 840.)) goto L50;
	  if ((underlyings(30,1) <= 8288.)) goto L50;
	  if ((underlyings(31,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(31,2) <= 840.)) goto L50;
	  if ((underlyings(31,1) <= 8288.)) goto L50;
	  if ((underlyings(32,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(32,2) <= 840.)) goto L50;
	  if ((underlyings(32,1) <= 8288.)) goto L50;
	  if ((underlyings(33,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(33,2) <= 840.)) goto L50;
	  if ((underlyings(33,1) <= 8288.)) goto L50;
	  if ((underlyings(34,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(34,2) <= 840.)) goto L50;
	  if ((underlyings(34,1) <= 8288.)) goto L50;
	  if ((underlyings(35,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(35,2) <= 840.)) goto L50;
	  if ((underlyings(35,1) <= 8288.)) goto L50;
	  if ((underlyings(36,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(36,2) <= 840.)) goto L50;
	  if ((underlyings(36,1) <= 8288.)) goto L50;
	  if ((underlyings(37,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(37,2) <= 840.)) goto L50;
	  if ((underlyings(37,1) <= 8288.)) goto L50;
	  if ((underlyings(38,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(38,2) <= 840.)) goto L50;
	  if ((underlyings(38,1) <= 8288.)) goto L50;
	  if ((underlyings(39,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(39,2) <= 840.)) goto L50;
	  if ((underlyings(39,1) <= 8288.)) goto L50;
	  if ((underlyings(40,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(40,2) <= 840.)) goto L50;
	  if ((underlyings(40,1) <= 8288.)) goto L50;
	  if ((underlyings(41,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(41,2) <= 840.)) goto L50;
	  if ((underlyings(41,1) <= 8288.)) goto L50;
	  if ((underlyings(42,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(42,2) <= 840.)) goto L50;
	  if ((underlyings(42,1) <= 8288.)) goto L50;
	  if ((underlyings(43,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(43,2) <= 840.)) goto L50;
	  if ((underlyings(43,1) <= 8288.)) goto L50;
	  if ((underlyings(44,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(44,2) <= 840.)) goto L50;
	  if ((underlyings(44,1) <= 8288.)) goto L50;
	  if ((underlyings(45,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(45,2) <= 840.)) goto L50;
	  if ((underlyings(45,1) <= 8288.)) goto L50;
	  if ((underlyings(46,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(46,2) <= 840.)) goto L50;
	  if ((underlyings(46,1) <= 8288.)) goto L50;
	  if ((underlyings(47,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(47,2) <= 840.)) goto L50;
	  if ((underlyings(47,1) <= 8288.)) goto L50;
	  if ((underlyings(48,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(48,2) <= 840.)) goto L50;
	  if ((underlyings(48,1) <= 8288.)) goto L50;
	  if ((underlyings(49,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(49,2) <= 840.)) goto L50;
	  if ((underlyings(49,1) <= 8288.)) goto L50;
	  if ((underlyings(50,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(50,2) <= 840.)) goto L50;
	  if ((underlyings(50,1) <= 8288.)) goto L50;
	  if ((underlyings(51,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(51,2) <= 840.)) goto L50;
	  if ((underlyings(51,1) <= 8288.)) goto L50;
	  if ((underlyings(52,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(52,2) <= 840.)) goto L50;
	  if ((underlyings(52,1) <= 8288.)) goto L50;
	  if ((underlyings(53,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(53,2) <= 840.)) goto L50;
	  if ((underlyings(53,1) <= 8288.)) goto L50;
	  if ((underlyings(54,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(54,2) <= 840.)) goto L50;
	  if ((underlyings(54,1) <= 8288.)) goto L50;
	  if ((underlyings(55,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(55,2) <= 840.)) goto L50;
	  if ((underlyings(55,1) <= 8288.)) goto L50;
	  if ((underlyings(56,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(56,2) <= 840.)) goto L50;
	  if ((underlyings(56,1) <= 8288.)) goto L50;
	  if ((underlyings(57,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(57,2) <= 840.)) goto L50;
	  if ((underlyings(57,1) <= 8288.)) goto L50;
	  if ((underlyings(58,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(58,2) <= 840.)) goto L50;
	  if ((underlyings(58,1) <= 8288.)) goto L50;
	  if ((underlyings(59,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(59,2) <= 840.)) goto L50;
	  if ((underlyings(59,1) <= 8288.)) goto L50;
	  if ((underlyings(60,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(60,2) <= 840.)) goto L50;
	  if ((underlyings(60,1) <= 8288.)) goto L50;
	  if ((underlyings(61,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(61,2) <= 840.)) goto L50;
	  if ((underlyings(61,1) <= 8288.)) goto L50;
	  if ((underlyings(62,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(62,2) <= 840.)) goto L50;
	  if ((underlyings(62,1) <= 8288.)) goto L50;
	  if ((underlyings(63,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(63,2) <= 840.)) goto L50;
	  if ((underlyings(63,1) <= 8288.)) goto L50;
	  if ((underlyings(64,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(64,2) <= 840.)) goto L50;
	  if ((underlyings(64,1) <= 8288.)) goto L50;
	  if ((underlyings(65,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(65,2) <= 840.)) goto L50;
	  if ((underlyings(65,1) <= 8288.)) goto L50;
	  if ((underlyings(66,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(66,2) <= 840.)) goto L50;
	  if ((underlyings(66,1) <= 8288.)) goto L50;
	  if ((underlyings(67,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(67,2) <= 840.)) goto L50;
	  if ((underlyings(67,1) <= 8288.)) goto L50;
	  if ((underlyings(68,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(68,2) <= 840.)) goto L50;
	  if ((underlyings(68,1) <= 8288.)) goto L50;
	  if ((underlyings(69,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(69,2) <= 840.)) goto L50;
	  if ((underlyings(69,1) <= 8288.)) goto L50;
	  if ((underlyings(70,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(70,2) <= 840.)) goto L50;
	  if ((underlyings(70,1) <= 8288.)) goto L50;
	  if ((underlyings(71,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(71,2) <= 840.)) goto L50;
	  if ((underlyings(71,1) <= 8288.)) goto L50;
	  if ((underlyings(72,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(72,2) <= 840.)) goto L50;
	  if ((underlyings(72,1) <= 8288.)) goto L50;
	  if ((underlyings(73,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(73,2) <= 840.)) goto L50;
	  if ((underlyings(73,1) <= 8288.)) goto L50;
	  if ((underlyings(74,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(74,2) <= 840.)) goto L50;
	  if ((underlyings(74,1) <= 8288.)) goto L50;
	  if ((underlyings(75,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(75,2) <= 840.)) goto L50;
	  if ((underlyings(75,1) <= 8288.)) goto L50;
	  if ((underlyings(76,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(76,2) <= 840.)) goto L50;
	  if ((underlyings(76,1) <= 8288.)) goto L50;
	  if ((underlyings(77,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(77,2) <= 840.)) goto L50;
	  if ((underlyings(77,1) <= 8288.)) goto L50;
	  if ((underlyings(78,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(78,2) <= 840.)) goto L50;
	  if ((underlyings(78,1) <= 8288.)) goto L50;
	  if ((underlyings(79,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(79,2) <= 840.)) goto L50;
	  if ((underlyings(79,1) <= 8288.)) goto L50;
	  if ((underlyings(80,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(80,2) <= 840.)) goto L50;
	  if ((underlyings(80,1) <= 8288.)) goto L50;
	  if ((underlyings(81,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(81,2) <= 840.)) goto L50;
	  if ((underlyings(81,1) <= 8288.)) goto L50;
	  if ((underlyings(82,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(82,2) <= 840.)) goto L50;
	  if ((underlyings(82,1) <= 8288.)) goto L50;
	  if ((underlyings(83,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(83,2) <= 840.)) goto L50;
	  if ((underlyings(83,1) <= 8288.)) goto L50;
	  if ((underlyings(84,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(84,2) <= 840.)) goto L50;
	  if ((underlyings(84,1) <= 8288.)) goto L50;
	  if ((underlyings(85,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(85,2) <= 840.)) goto L50;
	  if ((underlyings(85,1) <= 8288.)) goto L50;
	  if ((underlyings(86,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(86,2) <= 840.)) goto L50;
	  if ((underlyings(86,1) <= 8288.)) goto L50;
	  if ((underlyings(87,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(87,2) <= 840.)) goto L50;
	  if ((underlyings(87,1) <= 8288.)) goto L50;
	  if ((underlyings(88,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(88,2) <= 840.)) goto L50;
	  if ((underlyings(88,1) <= 8288.)) goto L50;
	  if ((underlyings(89,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(89,2) <= 840.)) goto L50;
	  if ((underlyings(89,1) <= 8288.)) goto L50;
	  if ((underlyings(90,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(90,2) <= 840.)) goto L50;
	  if ((underlyings(90,1) <= 8288.)) goto L50;
	  if ((underlyings(91,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(91,2) <= 840.)) goto L50;
	  if ((underlyings(91,1) <= 8288.)) goto L50;
	  if ((underlyings(92,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(92,2) <= 840.)) goto L50;
	  if ((underlyings(92,1) <= 8288.)) goto L50;
	  if ((underlyings(93,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(93,2) <= 840.)) goto L50;
	  if ((underlyings(93,1) <= 8288.)) goto L50;
	  if ((underlyings(94,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(94,2) <= 840.)) goto L50;
	  if ((underlyings(94,1) <= 8288.)) goto L50;
	  if ((underlyings(95,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(95,2) <= 840.)) goto L50;
	  if ((underlyings(95,1) <= 8288.)) goto L50;
	  if ((underlyings(96,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(96,2) <= 840.)) goto L50;
	  if ((underlyings(96,1) <= 8288.)) goto L50;
	  if ((underlyings(97,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(97,2) <= 840.)) goto L50;
	  if ((underlyings(97,1) <= 8288.)) goto L50;
	  if ((underlyings(98,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(98,2) <= 840.)) goto L50;
	  if ((underlyings(98,1) <= 8288.)) goto L50;
	  if ((underlyings(99,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(99,2) <= 840.)) goto L50;
	  if ((underlyings(99,1) <= 8288.)) goto L50;
	  if ((underlyings(100,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(100,2) <= 840.)) goto L50;
	  if ((underlyings(100,1) <= 8288.)) goto L50;
	  if ((underlyings(101,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(101,2) <= 840.)) goto L50;
	  if ((underlyings(101,1) <= 8288.)) goto L50;
	  if ((underlyings(102,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(102,2) <= 840.)) goto L50;
	  if ((underlyings(102,1) <= 8288.)) goto L50;
	  if ((underlyings(103,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(103,2) <= 840.)) goto L50;
	  if ((underlyings(103,1) <= 8288.)) goto L50;
	  if ((underlyings(104,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(104,2) <= 840.)) goto L50;
	  if ((underlyings(104,1) <= 8288.)) goto L50;
	  if ((underlyings(105,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(105,2) <= 840.)) goto L50;
	  if ((underlyings(105,1) <= 8288.)) goto L50;
	  if ((underlyings(106,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(106,2) <= 840.)) goto L50;
	  if ((underlyings(106,1) <= 8288.)) goto L50;
	  if ((underlyings(107,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(107,2) <= 840.)) goto L50;
	  if ((underlyings(107,1) <= 8288.)) goto L50;
	  if ((underlyings(108,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(108,2) <= 840.)) goto L50;
	  if ((underlyings(108,1) <= 8288.)) goto L50;
	  if ((underlyings(109,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(109,2) <= 840.)) goto L50;
	  if ((underlyings(109,1) <= 8288.)) goto L50;
	  if ((underlyings(110,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(110,2) <= 840.)) goto L50;
	  if ((underlyings(110,1) <= 8288.)) goto L50;
	  if ((underlyings(111,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(111,2) <= 840.)) goto L50;
	  if ((underlyings(111,1) <= 8288.)) goto L50;
	  if ((underlyings(112,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(112,2) <= 840.)) goto L50;
	  if ((underlyings(112,1) <= 8288.)) goto L50;
	  if ((underlyings(113,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(113,2) <= 840.)) goto L50;
	  if ((underlyings(113,1) <= 8288.)) goto L50;
	  if ((underlyings(114,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(114,2) <= 840.)) goto L50;
	  if ((underlyings(114,1) <= 8288.)) goto L50;
	  if ((underlyings(115,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(115,2) <= 840.)) goto L50;
	  if ((underlyings(115,1) <= 8288.)) goto L50;
	  if ((underlyings(116,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(116,2) <= 840.)) goto L50;
	  if ((underlyings(116,1) <= 8288.)) goto L50;
	  if ((underlyings(117,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(117,2) <= 840.)) goto L50;
	  if ((underlyings(117,1) <= 8288.)) goto L50;
	  if ((underlyings(118,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(118,2) <= 840.)) goto L50;
	  if ((underlyings(118,1) <= 8288.)) goto L50;
	  if ((underlyings(119,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(119,2) <= 840.)) goto L50;
	  if ((underlyings(119,1) <= 8288.)) goto L50;
	  if ((underlyings(120,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(120,2) <= 840.)) goto L50;
	  if ((underlyings(120,1) <= 8288.)) goto L50;
	  if ((underlyings(121,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(121,2) <= 840.)) goto L50;
	  if ((underlyings(121,1) <= 8288.)) goto L50;
	  if ((underlyings(122,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(122,2) <= 840.)) goto L50;
	  if ((underlyings(122,1) <= 8288.)) goto L50;
	  if ((underlyings(123,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(123,2) <= 840.)) goto L50;
	  if ((underlyings(123,1) <= 8288.)) goto L50;
	  if ((underlyings(124,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(124,2) <= 840.)) goto L50;
	  if ((underlyings(124,1) <= 8288.)) goto L50;
	  if ((underlyings(125,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(125,2) <= 840.)) goto L50;
	  if ((underlyings(125,1) <= 8288.)) goto L50;
	  if ((underlyings(126,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(126,2) <= 840.)) goto L50;
	  if ((underlyings(126,1) <= 8288.)) goto L50;
	  if ((underlyings(127,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(127,2) <= 840.)) goto L50;
	  if ((underlyings(127,1) <= 8288.)) goto L50;
	  if ((underlyings(128,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(128,2) <= 840.)) goto L50;
	  if ((underlyings(128,1) <= 8288.)) goto L50;
	  if ((underlyings(129,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(129,2) <= 840.)) goto L50;
	  if ((underlyings(129,1) <= 8288.)) goto L50;
	  if ((underlyings(130,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(130,2) <= 840.)) goto L50;
	  if ((underlyings(130,1) <= 8288.)) goto L50;
	  if ((underlyings(131,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(131,2) <= 840.)) goto L50;
	  if ((underlyings(131,1) <= 8288.)) goto L50;
	  if ((underlyings(132,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(132,2) <= 840.)) goto L50;
	  if ((underlyings(132,1) <= 8288.)) goto L50;
	  if ((underlyings(133,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(133,2) <= 840.)) goto L50;
	  if ((underlyings(133,1) <= 8288.)) goto L50;
	  if ((underlyings(134,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(134,2) <= 840.)) goto L50;
	  if ((underlyings(134,1) <= 8288.)) goto L50;
	  if ((underlyings(135,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(135,2) <= 840.)) goto L50;
	  if ((underlyings(135,1) <= 8288.)) goto L50;
	  if ((underlyings(136,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(136,2) <= 840.)) goto L50;
	  if ((underlyings(136,1) <= 8288.)) goto L50;
	  if ((underlyings(137,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(137,2) <= 840.)) goto L50;
	  if ((underlyings(137,1) <= 8288.)) goto L50;
	  if ((underlyings(138,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(138,2) <= 840.)) goto L50;
	  if ((underlyings(138,1) <= 8288.)) goto L50;
	  if ((underlyings(139,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(139,2) <= 840.)) goto L50;
	  if ((underlyings(139,1) <= 8288.)) goto L50;
	  if ((underlyings(140,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(140,2) <= 840.)) goto L50;
	  if ((underlyings(140,1) <= 8288.)) goto L50;
	  if ((underlyings(141,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(141,2) <= 840.)) goto L50;
	  if ((underlyings(141,1) <= 8288.)) goto L50;
	  if ((underlyings(142,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(142,2) <= 840.)) goto L50;
	  if ((underlyings(142,1) <= 8288.)) goto L50;
	  if ((underlyings(143,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(143,2) <= 840.)) goto L50;
	  if ((underlyings(143,1) <= 8288.)) goto L50;
	  if ((underlyings(144,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(144,2) <= 840.)) goto L50;
	  if ((underlyings(144,1) <= 8288.)) goto L50;
	  if ((underlyings(145,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(145,2) <= 840.)) goto L50;
	  if ((underlyings(145,1) <= 8288.)) goto L50;
	  if ((underlyings(146,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(146,2) <= 840.)) goto L50;
	  if ((underlyings(146,1) <= 8288.)) goto L50;
	  if ((underlyings(147,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(147,2) <= 840.)) goto L50;
	  if ((underlyings(147,1) <= 8288.)) goto L50;
	  if ((underlyings(148,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(148,2) <= 840.)) goto L50;
	  if ((underlyings(148,1) <= 8288.)) goto L50;
	  if ((underlyings(149,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(149,2) <= 840.)) goto L50;
	  if ((underlyings(149,1) <= 8288.)) goto L50;
	  if ((underlyings(150,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(150,2) <= 840.)) goto L50;
	  if ((underlyings(150,1) <= 8288.)) goto L50;
	  if ((underlyings(151,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(151,2) <= 840.)) goto L50;
	  if ((underlyings(151,1) <= 8288.)) goto L50;
	  if ((underlyings(152,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(152,2) <= 840.)) goto L50;
	  if ((underlyings(152,1) <= 8288.)) goto L50;
	  if ((underlyings(153,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(153,2) <= 840.)) goto L50;
	  if ((underlyings(153,1) <= 8288.)) goto L50;
	  if ((underlyings(154,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(154,2) <= 840.)) goto L50;
	  if ((underlyings(154,1) <= 8288.)) goto L50;
	  if ((underlyings(155,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(155,2) <= 840.)) goto L50;
	  if ((underlyings(155,1) <= 8288.)) goto L50;
	  if ((underlyings(156,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(156,2) <= 840.)) goto L50;
	  if ((underlyings(156,1) <= 8288.)) goto L50;
	  if ((underlyings(157,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(157,2) <= 840.)) goto L50;
	  if ((underlyings(157,1) <= 8288.)) goto L50;
	  if ((underlyings(158,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(158,2) <= 840.)) goto L50;
	  if ((underlyings(158,1) <= 8288.)) goto L50;
	  if ((underlyings(159,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(159,2) <= 840.)) goto L50;
	  if ((underlyings(159,1) <= 8288.)) goto L50;
	  if ((underlyings(160,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(160,2) <= 840.)) goto L50;
	  if ((underlyings(160,1) <= 8288.)) goto L50;
	  if ((underlyings(161,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(161,2) <= 840.)) goto L50;
	  if ((underlyings(161,1) <= 8288.)) goto L50;
	  if ((underlyings(162,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(162,2) <= 840.)) goto L50;
	  if ((underlyings(162,1) <= 8288.)) goto L50;
	  if ((underlyings(163,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(163,2) <= 840.)) goto L50;
	  if ((underlyings(163,1) <= 8288.)) goto L50;
	  if ((underlyings(164,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(164,2) <= 840.)) goto L50;
	  if ((underlyings(164,1) <= 8288.)) goto L50;
	  if ((underlyings(165,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(165,2) <= 840.)) goto L50;
	  if ((underlyings(165,1) <= 8288.)) goto L50;
	  if ((underlyings(166,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(166,2) <= 840.)) goto L50;
	  if ((underlyings(166,1) <= 8288.)) goto L50;
	  if ((underlyings(167,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(167,2) <= 840.)) goto L50;
	  if ((underlyings(167,1) <= 8288.)) goto L50;
	  if ((underlyings(168,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(168,2) <= 840.)) goto L50;
	  if ((underlyings(168,1) <= 8288.)) goto L50;
	  if ((underlyings(169,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(169,2) <= 840.)) goto L50;
	  if ((underlyings(169,1) <= 8288.)) goto L50;
	  if ((underlyings(170,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(170,2) <= 840.)) goto L50;
	  if ((underlyings(170,1) <= 8288.)) goto L50;
	  if ((underlyings(171,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(171,2) <= 840.)) goto L50;
	  if ((underlyings(171,1) <= 8288.)) goto L50;
	  if ((underlyings(172,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(172,2) <= 840.)) goto L50;
	  if ((underlyings(172,1) <= 8288.)) goto L50;
	  if ((underlyings(173,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(173,2) <= 840.)) goto L50;
	  if ((underlyings(173,1) <= 8288.)) goto L50;
	  if ((underlyings(174,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(174,2) <= 840.)) goto L50;
	  if ((underlyings(174,1) <= 8288.)) goto L50;
	  if ((underlyings(175,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(175,2) <= 840.)) goto L50;
	  if ((underlyings(175,1) <= 8288.)) goto L50;
	  if ((underlyings(176,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(176,2) <= 840.)) goto L50;
	  if ((underlyings(176,1) <= 8288.)) goto L50;
	  if ((underlyings(177,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(177,2) <= 840.)) goto L50;
	  if ((underlyings(177,1) <= 8288.)) goto L50;
	  if ((underlyings(178,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(178,2) <= 840.)) goto L50;
	  if ((underlyings(178,1) <= 8288.)) goto L50;
	  if ((underlyings(179,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(179,2) <= 840.)) goto L50;
	  if ((underlyings(179,1) <= 8288.)) goto L50;
	  if ((underlyings(180,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(180,2) <= 840.)) goto L50;
	  if ((underlyings(180,1) <= 8288.)) goto L50;
	  if ((underlyings(181,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(181,2) <= 840.)) goto L50;
	  if ((underlyings(181,1) <= 8288.)) goto L50;
	  if ((underlyings(182,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(182,2) <= 840.)) goto L50;
	  if ((underlyings(182,1) <= 8288.)) goto L50;
	  if ((underlyings(183,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(183,2) <= 840.)) goto L50;
	  if ((underlyings(183,1) <= 8288.)) goto L50;
	  if ((underlyings(184,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(184,2) <= 840.)) goto L50;
	  if ((underlyings(184,1) <= 8288.)) goto L50;
	  if ((underlyings(185,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(185,2) <= 840.)) goto L50;
	  if ((underlyings(185,1) <= 8288.)) goto L50;
	  if ((underlyings(186,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(186,2) <= 840.)) goto L50;
	  if ((underlyings(186,1) <= 8288.)) goto L50;
	  if ((underlyings(187,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(187,2) <= 840.)) goto L50;
	  if ((underlyings(187,1) <= 8288.)) goto L50;
	  if ((underlyings(188,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(188,2) <= 840.)) goto L50;
	  if ((underlyings(188,1) <= 8288.)) goto L50;
	  if ((underlyings(189,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(189,2) <= 840.)) goto L50;
	  if ((underlyings(189,1) <= 8288.)) goto L50;
	  if ((underlyings(190,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(190,2) <= 840.)) goto L50;
	  if ((underlyings(190,1) <= 8288.)) goto L50;
	  if ((underlyings(191,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(191,2) <= 840.)) goto L50;
	  if ((underlyings(191,1) <= 8288.)) goto L50;
	  if ((underlyings(192,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(192,2) <= 840.)) goto L50;
	  if ((underlyings(192,1) <= 8288.)) goto L50;
	  if ((underlyings(193,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(193,2) <= 840.)) goto L50;
	  if ((underlyings(193,1) <= 8288.)) goto L50;
	  if ((underlyings(194,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(194,2) <= 840.)) goto L50;
	  if ((underlyings(194,1) <= 8288.)) goto L50;
	  if ((underlyings(195,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(195,2) <= 840.)) goto L50;
	  if ((underlyings(195,1) <= 8288.)) goto L50;
	  if ((underlyings(196,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(196,2) <= 840.)) goto L50;
	  if ((underlyings(196,1) <= 8288.)) goto L50;
	  if ((underlyings(197,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(197,2) <= 840.)) goto L50;
	  if ((underlyings(197,1) <= 8288.)) goto L50;
	  if ((underlyings(198,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(198,2) <= 840.)) goto L50;
	  if ((underlyings(198,1) <= 8288.)) goto L50;
	  if ((underlyings(199,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(199,2) <= 840.)) goto L50;
	  if ((underlyings(199,1) <= 8288.)) goto L50;
	  if ((underlyings(200,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(200,2) <= 840.)) goto L50;
	  if ((underlyings(200,1) <= 8288.)) goto L50;
	  if ((underlyings(201,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(201,2) <= 840.)) goto L50;
	  if ((underlyings(201,1) <= 8288.)) goto L50;
	  if ((underlyings(202,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(202,2) <= 840.)) goto L50;
	  if ((underlyings(202,1) <= 8288.)) goto L50;
	  if ((underlyings(203,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(203,2) <= 840.)) goto L50;
	  if ((underlyings(203,1) <= 8288.)) goto L50;
	  if ((underlyings(204,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(204,2) <= 840.)) goto L50;
	  if ((underlyings(204,1) <= 8288.)) goto L50;
	  if ((underlyings(205,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(205,2) <= 840.)) goto L50;
	  if ((underlyings(205,1) <= 8288.)) goto L50;
	  if ((underlyings(206,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(206,2) <= 840.)) goto L50;
	  if ((underlyings(206,1) <= 8288.)) goto L50;
	  if ((underlyings(207,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(207,2) <= 840.)) goto L50;
	  if ((underlyings(207,1) <= 8288.)) goto L50;
	  if ((underlyings(208,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(208,2) <= 840.)) goto L50;
	  if ((underlyings(208,1) <= 8288.)) goto L50;
	  if ((underlyings(209,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(209,2) <= 840.)) goto L50;
	  if ((underlyings(209,1) <= 8288.)) goto L50;
	  if ((underlyings(210,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(210,2) <= 840.)) goto L50;
	  if ((underlyings(210,1) <= 8288.)) goto L50;
	  if ((underlyings(211,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(211,2) <= 840.)) goto L50;
	  if ((underlyings(211,1) <= 8288.)) goto L50;
	  if ((underlyings(212,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(212,2) <= 840.)) goto L50;
	  if ((underlyings(212,1) <= 8288.)) goto L50;
	  if ((underlyings(213,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(213,2) <= 840.)) goto L50;
	  if ((underlyings(213,1) <= 8288.)) goto L50;
	  if ((underlyings(214,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(214,2) <= 840.)) goto L50;
	  if ((underlyings(214,1) <= 8288.)) goto L50;
	  if ((underlyings(215,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(215,2) <= 840.)) goto L50;
	  if ((underlyings(215,1) <= 8288.)) goto L50;
	  if ((underlyings(216,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(216,2) <= 840.)) goto L50;
	  if ((underlyings(216,1) <= 8288.)) goto L50;
	  if ((underlyings(217,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(217,2) <= 840.)) goto L50;
	  if ((underlyings(217,1) <= 8288.)) goto L50;
	  if ((underlyings(218,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(218,2) <= 840.)) goto L50;
	  if ((underlyings(218,1) <= 8288.)) goto L50;
	  if ((underlyings(219,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(219,2) <= 840.)) goto L50;
	  if ((underlyings(219,1) <= 8288.)) goto L50;
	  if ((underlyings(220,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(220,2) <= 840.)) goto L50;
	  if ((underlyings(220,1) <= 8288.)) goto L50;
	  if ((underlyings(221,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(221,2) <= 840.)) goto L50;
	  if ((underlyings(221,1) <= 8288.)) goto L50;
	  if ((underlyings(222,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(222,2) <= 840.)) goto L50;
	  if ((underlyings(222,1) <= 8288.)) goto L50;
	  if ((underlyings(223,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(223,2) <= 840.)) goto L50;
	  if ((underlyings(223,1) <= 8288.)) goto L50;
	  if ((underlyings(224,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(224,2) <= 840.)) goto L50;
	  if ((underlyings(224,1) <= 8288.)) goto L50;
	  if ((underlyings(225,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(225,2) <= 840.)) goto L50;
	  if ((underlyings(225,1) <= 8288.)) goto L50;
	  if ((underlyings(226,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(226,2) <= 840.)) goto L50;
	  if ((underlyings(226,1) <= 8288.)) goto L50;
	  if ((underlyings(227,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(227,2) <= 840.)) goto L50;
	  if ((underlyings(227,1) <= 8288.)) goto L50;
	  if ((underlyings(228,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(228,2) <= 840.)) goto L50;
	  if ((underlyings(228,1) <= 8288.)) goto L50;
	  if ((underlyings(229,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(229,2) <= 840.)) goto L50;
	  if ((underlyings(229,1) <= 8288.)) goto L50;
	  if ((underlyings(230,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(230,2) <= 840.)) goto L50;
	  if ((underlyings(230,1) <= 8288.)) goto L50;
	  if ((underlyings(231,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(231,2) <= 840.)) goto L50;
	  if ((underlyings(231,1) <= 8288.)) goto L50;
	  if ((underlyings(232,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(232,2) <= 840.)) goto L50;
	  if ((underlyings(232,1) <= 8288.)) goto L50;
	  if ((underlyings(233,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(233,2) <= 840.)) goto L50;
	  if ((underlyings(233,1) <= 8288.)) goto L50;
	  if ((underlyings(234,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(234,2) <= 840.)) goto L50;
	  if ((underlyings(234,1) <= 8288.)) goto L50;
	  if ((underlyings(235,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(235,2) <= 840.)) goto L50;
	  if ((underlyings(235,1) <= 8288.)) goto L50;
	  if ((underlyings(236,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(236,2) <= 840.)) goto L50;
	  if ((underlyings(236,1) <= 8288.)) goto L50;
	  if ((underlyings(237,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(237,2) <= 840.)) goto L50;
	  if ((underlyings(237,1) <= 8288.)) goto L50;
	  if ((underlyings(238,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(238,2) <= 840.)) goto L50;
	  if ((underlyings(238,1) <= 8288.)) goto L50;
	  if ((underlyings(239,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(239,2) <= 840.)) goto L50;
	  if ((underlyings(239,1) <= 8288.)) goto L50;
	  if ((underlyings(240,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(240,2) <= 840.)) goto L50;
	  if ((underlyings(240,1) <= 8288.)) goto L50;
	  if ((underlyings(241,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(241,2) <= 840.)) goto L50;
	  if ((underlyings(241,1) <= 8288.)) goto L50;
	  if ((underlyings(242,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(242,2) <= 840.)) goto L50;
	  if ((underlyings(242,1) <= 8288.)) goto L50;
	  if ((underlyings(243,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(243,2) <= 840.)) goto L50;
	  if ((underlyings(243,1) <= 8288.)) goto L50;
	  if ((underlyings(244,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(244,2) <= 840.)) goto L50;
	  if ((underlyings(244,1) <= 8288.)) goto L50;
	  if ((underlyings(245,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(245,2) <= 840.)) goto L50;
	  if ((underlyings(245,1) <= 8288.)) goto L50;
	  if ((underlyings(246,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(246,2) <= 840.)) goto L50;
	  if ((underlyings(246,1) <= 8288.)) goto L50;
	  if ((underlyings(247,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(247,2) <= 840.)) goto L50;
	  if ((underlyings(247,1) <= 8288.)) goto L50;
	  if ((underlyings(248,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(248,2) <= 840.)) goto L50;
	  if ((underlyings(248,1) <= 8288.)) goto L50;
	  if ((underlyings(249,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(249,2) <= 840.)) goto L50;
	  if ((underlyings(249,1) <= 8288.)) goto L50;
	  if ((underlyings(250,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(250,2) <= 840.)) goto L50;
	  if ((underlyings(250,1) <= 8288.)) goto L50;
	  if ((underlyings(251,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(251,2) <= 840.)) goto L50;
	  if ((underlyings(251,1) <= 8288.)) goto L50;
	  if ((underlyings(252,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(252,2) <= 840.)) goto L50;
	  if ((underlyings(252,1) <= 8288.)) goto L50;
	  if ((underlyings(253,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(253,2) <= 840.)) goto L50;
	  if ((underlyings(253,1) <= 8288.)) goto L50;
	  if ((underlyings(254,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(254,2) <= 840.)) goto L50;
	  if ((underlyings(254,1) <= 8288.)) goto L50;
	  if ((underlyings(255,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(255,2) <= 840.)) goto L50;
	  if ((underlyings(255,1) <= 8288.)) goto L50;
	  if ((underlyings(256,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(256,2) <= 840.)) goto L50;
	  if ((underlyings(256,1) <= 8288.)) goto L50;
	  if ((underlyings(257,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(257,2) <= 840.)) goto L50;
	  if ((underlyings(257,1) <= 8288.)) goto L50;
	  if ((underlyings(258,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(258,2) <= 840.)) goto L50;
	  if ((underlyings(258,1) <= 8288.)) goto L50;
	  if ((underlyings(259,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(259,2) <= 840.)) goto L50;
	  if ((underlyings(259,1) <= 8288.)) goto L50;
	  if ((underlyings(260,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(260,2) <= 840.)) goto L50;
	  if ((underlyings(260,1) <= 8288.)) goto L50;
	  if ((underlyings(261,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(261,2) <= 840.)) goto L50;
	  if ((underlyings(261,1) <= 8288.)) goto L50;
	  if ((underlyings(262,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(262,2) <= 840.)) goto L50;
	  if ((underlyings(262,1) <= 8288.)) goto L50;
	  if ((underlyings(263,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(263,2) <= 840.)) goto L50;
	  if ((underlyings(263,1) <= 8288.)) goto L50;
	  if ((underlyings(264,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(264,2) <= 840.)) goto L50;
	  if ((underlyings(264,1) <= 8288.)) goto L50;
	  if ((underlyings(265,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(265,2) <= 840.)) goto L50;
	  if ((underlyings(265,1) <= 8288.)) goto L50;
	  if ((underlyings(266,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(266,2) <= 840.)) goto L50;
	  if ((underlyings(266,1) <= 8288.)) goto L50;
	  if ((underlyings(267,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(267,2) <= 840.)) goto L50;
	  if ((underlyings(267,1) <= 8288.)) goto L50;
	  if ((underlyings(268,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(268,2) <= 840.)) goto L50;
	  if ((underlyings(268,1) <= 8288.)) goto L50;
	  if ((underlyings(269,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(269,2) <= 840.)) goto L50;
	  if ((underlyings(269,1) <= 8288.)) goto L50;
	  if ((underlyings(270,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(270,2) <= 840.)) goto L50;
	  if ((underlyings(270,1) <= 8288.)) goto L50;
	  if ((underlyings(271,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(271,2) <= 840.)) goto L50;
	  if ((underlyings(271,1) <= 8288.)) goto L50;
	  if ((underlyings(272,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(272,2) <= 840.)) goto L50;
	  if ((underlyings(272,1) <= 8288.)) goto L50;
	  if ((underlyings(273,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(273,2) <= 840.)) goto L50;
	  if ((underlyings(273,1) <= 8288.)) goto L50;
	  if ((underlyings(274,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(274,2) <= 840.)) goto L50;
	  if ((underlyings(274,1) <= 8288.)) goto L50;
	  if ((underlyings(275,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(275,2) <= 840.)) goto L50;
	  if ((underlyings(275,1) <= 8288.)) goto L50;
	  if ((underlyings(276,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(276,2) <= 840.)) goto L50;
	  if ((underlyings(276,1) <= 8288.)) goto L50;
	  if ((underlyings(277,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(277,2) <= 840.)) goto L50;
	  if ((underlyings(277,1) <= 8288.)) goto L50;
	  if ((underlyings(278,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(278,2) <= 840.)) goto L50;
	  if ((underlyings(278,1) <= 8288.)) goto L50;
	  if ((underlyings(279,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(279,2) <= 840.)) goto L50;
	  if ((underlyings(279,1) <= 8288.)) goto L50;
	  if ((underlyings(280,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(280,2) <= 840.)) goto L50;
	  if ((underlyings(280,1) <= 8288.)) goto L50;
	  if ((underlyings(281,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(281,2) <= 840.)) goto L50;
	  if ((underlyings(281,1) <= 8288.)) goto L50;
	  if ((underlyings(282,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(282,2) <= 840.)) goto L50;
	  if ((underlyings(282,1) <= 8288.)) goto L50;
	  if ((underlyings(283,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(283,2) <= 840.)) goto L50;
	  if ((underlyings(283,1) <= 8288.)) goto L50;
	  if ((underlyings(284,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(284,2) <= 840.)) goto L50;
	  if ((underlyings(284,1) <= 8288.)) goto L50;
	  if ((underlyings(285,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(285,2) <= 840.)) goto L50;
	  if ((underlyings(285,1) <= 8288.)) goto L50;
	  if ((underlyings(286,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(286,2) <= 840.)) goto L50;
	  if ((underlyings(286,1) <= 8288.)) goto L50;
	  if ((underlyings(287,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(287,2) <= 840.)) goto L50;
	  if ((underlyings(287,1) <= 8288.)) goto L50;
	  if ((underlyings(288,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(288,2) <= 840.)) goto L50;
	  if ((underlyings(288,1) <= 8288.)) goto L50;
	  if ((underlyings(289,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(289,2) <= 840.)) goto L50;
	  if ((underlyings(289,1) <= 8288.)) goto L50;
	  if ((underlyings(290,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(290,2) <= 840.)) goto L50;
	  if ((underlyings(290,1) <= 8288.)) goto L50;
	  if ((underlyings(291,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(291,2) <= 840.)) goto L50;
	  if ((underlyings(291,1) <= 8288.)) goto L50;
	  if ((underlyings(292,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(292,2) <= 840.)) goto L50;
	  if ((underlyings(292,1) <= 8288.)) goto L50;
	  if ((underlyings(293,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(293,2) <= 840.)) goto L50;
	  if ((underlyings(293,1) <= 8288.)) goto L50;
	  if ((underlyings(294,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(294,2) <= 840.)) goto L50;
	  if ((underlyings(294,1) <= 8288.)) goto L50;
	  if ((underlyings(295,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(295,2) <= 840.)) goto L50;
	  if ((underlyings(295,1) <= 8288.)) goto L50;
	  if ((underlyings(296,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(296,2) <= 840.)) goto L50;
	  if ((underlyings(296,1) <= 8288.)) goto L50;
	  if ((underlyings(297,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(297,2) <= 840.)) goto L50;
	  if ((underlyings(297,1) <= 8288.)) goto L50;
	  if ((underlyings(298,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(298,2) <= 840.)) goto L50;
	  if ((underlyings(298,1) <= 8288.)) goto L50;
	  if ((underlyings(299,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(299,2) <= 840.)) goto L50;
	  if ((underlyings(299,1) <= 8288.)) goto L50;
	  if ((underlyings(300,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(300,2) <= 840.)) goto L50;
	  if ((underlyings(300,1) <= 8288.)) goto L50;
	  if ((underlyings(301,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(301,2) <= 840.)) goto L50;
	  if ((underlyings(301,1) <= 8288.)) goto L50;
	  if ((underlyings(302,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(302,2) <= 840.)) goto L50;
	  if ((underlyings(302,1) <= 8288.)) goto L50;
	  if ((underlyings(303,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(303,2) <= 840.)) goto L50;
	  if ((underlyings(303,1) <= 8288.)) goto L50;
	  if ((underlyings(304,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(304,2) <= 840.)) goto L50;
	  if ((underlyings(304,1) <= 8288.)) goto L50;
	  if ((underlyings(305,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(305,2) <= 840.)) goto L50;
	  if ((underlyings(305,1) <= 8288.)) goto L50;
	  if ((underlyings(306,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(306,2) <= 840.)) goto L50;
	  if ((underlyings(306,1) <= 8288.)) goto L50;
	  if ((underlyings(307,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(307,2) <= 840.)) goto L50;
	  if ((underlyings(307,1) <= 8288.)) goto L50;
	  if ((underlyings(308,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(308,2) <= 840.)) goto L50;
	  if ((underlyings(308,1) <= 8288.)) goto L50;
	  if ((underlyings(309,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(309,2) <= 840.)) goto L50;
	  if ((underlyings(309,1) <= 8288.)) goto L50;
	  if ((underlyings(310,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(310,2) <= 840.)) goto L50;
	  if ((underlyings(310,1) <= 8288.)) goto L50;
	  if ((underlyings(311,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(311,2) <= 840.)) goto L50;
	  if ((underlyings(311,1) <= 8288.)) goto L50;
	  if ((underlyings(312,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(312,2) <= 840.)) goto L50;
	  if ((underlyings(312,1) <= 8288.)) goto L50;
	  if ((underlyings(313,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(313,2) <= 840.)) goto L50;
	  if ((underlyings(313,1) <= 8288.)) goto L50;
	  if ((underlyings(314,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(314,2) <= 840.)) goto L50;
	  if ((underlyings(314,1) <= 8288.)) goto L50;
	  if ((underlyings(315,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(315,2) <= 840.)) goto L50;
	  if ((underlyings(315,1) <= 8288.)) goto L50;
	  if ((underlyings(316,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(316,2) <= 840.)) goto L50;
	  if ((underlyings(316,1) <= 8288.)) goto L50;
	  if ((underlyings(317,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(317,2) <= 840.)) goto L50;
	  if ((underlyings(317,1) <= 8288.)) goto L50;
	  if ((underlyings(318,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(318,2) <= 840.)) goto L50;
	  if ((underlyings(318,1) <= 8288.)) goto L50;
	  if ((underlyings(319,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(319,2) <= 840.)) goto L50;
	  if ((underlyings(319,1) <= 8288.)) goto L50;
	  if ((underlyings(320,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(320,2) <= 840.)) goto L50;
	  if ((underlyings(320,1) <= 8288.)) goto L50;
	  if ((underlyings(321,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(321,2) <= 840.)) goto L50;
	  if ((underlyings(321,1) <= 8288.)) goto L50;
	  if ((underlyings(322,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(322,2) <= 840.)) goto L50;
	  if ((underlyings(322,1) <= 8288.)) goto L50;
	  if ((underlyings(323,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(323,2) <= 840.)) goto L50;
	  if ((underlyings(323,1) <= 8288.)) goto L50;
	  if ((underlyings(324,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(324,2) <= 840.)) goto L50;
	  if ((underlyings(324,1) <= 8288.)) goto L50;
	  if ((underlyings(325,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(325,2) <= 840.)) goto L50;
	  if ((underlyings(325,1) <= 8288.)) goto L50;
	  if ((underlyings(326,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(326,2) <= 840.)) goto L50;
	  if ((underlyings(326,1) <= 8288.)) goto L50;
	  if ((underlyings(327,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(327,2) <= 840.)) goto L50;
	  if ((underlyings(327,1) <= 8288.)) goto L50;
	  if ((underlyings(328,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(328,2) <= 840.)) goto L50;
	  if ((underlyings(328,1) <= 8288.)) goto L50;
	  if ((underlyings(329,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(329,2) <= 840.)) goto L50;
	  if ((underlyings(329,1) <= 8288.)) goto L50;
	  if ((underlyings(330,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(330,2) <= 840.)) goto L50;
	  if ((underlyings(330,1) <= 8288.)) goto L50;
	  if ((underlyings(331,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(331,2) <= 840.)) goto L50;
	  if ((underlyings(331,1) <= 8288.)) goto L50;
	  if ((underlyings(332,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(332,2) <= 840.)) goto L50;
	  if ((underlyings(332,1) <= 8288.)) goto L50;
	  if ((underlyings(333,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(333,2) <= 840.)) goto L50;
	  if ((underlyings(333,1) <= 8288.)) goto L50;
	  if ((underlyings(334,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(334,2) <= 840.)) goto L50;
	  if ((underlyings(334,1) <= 8288.)) goto L50;
	  if ((underlyings(335,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(335,2) <= 840.)) goto L50;
	  if ((underlyings(335,1) <= 8288.)) goto L50;
	  if ((underlyings(336,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(336,2) <= 840.)) goto L50;
	  if ((underlyings(336,1) <= 8288.)) goto L50;
	  if ((underlyings(337,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(337,2) <= 840.)) goto L50;
	  if ((underlyings(337,1) <= 8288.)) goto L50;
	  if ((underlyings(338,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(338,2) <= 840.)) goto L50;
	  if ((underlyings(338,1) <= 8288.)) goto L50;
	  if ((underlyings(339,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(339,2) <= 840.)) goto L50;
	  if ((underlyings(339,1) <= 8288.)) goto L50;
	  if ((underlyings(340,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(340,2) <= 840.)) goto L50;
	  if ((underlyings(340,1) <= 8288.)) goto L50;
	  if ((underlyings(341,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(341,2) <= 840.)) goto L50;
	  if ((underlyings(341,1) <= 8288.)) goto L50;
	  if ((underlyings(342,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(342,2) <= 840.)) goto L50;
	  if ((underlyings(342,1) <= 8288.)) goto L50;
	  if ((underlyings(343,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(343,2) <= 840.)) goto L50;
	  if ((underlyings(343,1) <= 8288.)) goto L50;
	  if ((underlyings(344,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(344,2) <= 840.)) goto L50;
	  if ((underlyings(344,1) <= 8288.)) goto L50;
	  if ((underlyings(345,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(345,2) <= 840.)) goto L50;
	  if ((underlyings(345,1) <= 8288.)) goto L50;
	  if ((underlyings(346,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(346,2) <= 840.)) goto L50;
	  if ((underlyings(346,1) <= 8288.)) goto L50;
	  if ((underlyings(347,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(347,2) <= 840.)) goto L50;
	  if ((underlyings(347,1) <= 8288.)) goto L50;
	  if ((underlyings(348,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(348,2) <= 840.)) goto L50;
	  if ((underlyings(348,1) <= 8288.)) goto L50;
	  if ((underlyings(349,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(349,2) <= 840.)) goto L50;
	  if ((underlyings(349,1) <= 8288.)) goto L50;
	  if ((underlyings(350,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(350,2) <= 840.)) goto L50;
	  if ((underlyings(350,1) <= 8288.)) goto L50;
	  if ((underlyings(351,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(351,2) <= 840.)) goto L50;
	  if ((underlyings(351,1) <= 8288.)) goto L50;
	  if ((underlyings(352,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(352,2) <= 840.)) goto L50;
	  if ((underlyings(352,1) <= 8288.)) goto L50;
	  if ((underlyings(353,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(353,2) <= 840.)) goto L50;
	  if ((underlyings(353,1) <= 8288.)) goto L50;
	  if ((underlyings(354,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(354,2) <= 840.)) goto L50;
	  if ((underlyings(354,1) <= 8288.)) goto L50;
	  if ((underlyings(355,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(355,2) <= 840.)) goto L50;
	  if ((underlyings(355,1) <= 8288.)) goto L50;
	  if ((underlyings(356,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(356,2) <= 840.)) goto L50;
	  if ((underlyings(356,1) <= 8288.)) goto L50;
	  if ((underlyings(357,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(357,2) <= 840.)) goto L50;
	  if ((underlyings(357,1) <= 8288.)) goto L50;
	  if ((underlyings(358,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(358,2) <= 840.)) goto L50;
	  if ((underlyings(358,1) <= 8288.)) goto L50;
	  if ((underlyings(359,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(359,2) <= 840.)) goto L50;
	  if ((underlyings(359,1) <= 8288.)) goto L50;
	  if ((underlyings(360,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(360,2) <= 840.)) goto L50;
	  if ((underlyings(360,1) <= 8288.)) goto L50;
	  if ((underlyings(361,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(361,2) <= 840.)) goto L50;
	  if ((underlyings(361,1) <= 8288.)) goto L50;
	  if ((underlyings(362,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(362,2) <= 840.)) goto L50;
	  if ((underlyings(362,1) <= 8288.)) goto L50;
	  if ((underlyings(363,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(363,2) <= 840.)) goto L50;
	  if ((underlyings(363,1) <= 8288.)) goto L50;
	  if ((underlyings(364,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(364,2) <= 840.)) goto L50;
	  if ((underlyings(364,1) <= 8288.)) goto L50;
	  if ((underlyings(365,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(365,2) <= 840.)) goto L50;
	  if ((underlyings(365,1) <= 8288.)) goto L50;
	  if ((underlyings(366,0) <= 2630.6349999999998)) goto L50;
	  if ((underlyings(366,2) <= 840.)) goto L50;
	  x3309=(underlyings(366,1) <= 8288.);
	  goto L48;

L50:
    x3309=1;

L48:
    //model->notify_cash_flow(model, 0, 100., 0 /*2013-01-27, 1980-01-01, EUR*/);
    trajectory_inner( num_cash_flows, model_num, 0, 100., md_discts, vhat );

    if (((underlyings(366,0) < 3758.05) && x3309)) goto L40;
    if (((underlyings(366,2) < 1200.) && x3309)) goto L40;
    if (((underlyings(366,1) < 11840.) && x3309)) goto L40;

    //model->notify_cash_flow(model, 0, 1000., 1 /*2013-01-27, 2013-01-27, EUR*/); return;
    trajectory_inner( num_cash_flows, model_num, 1, 1000., md_discts, vhat ); return;

L40:
	  //model->notify_cash_flow(model, 0, (1000. * (1. + fmin(((underlyings(366,1) / 11840.) - 1.), fmin(((underlyings(366,2) / 1200.) - 1.), ((underlyings(366,0) / 3758.05) - 1.))))), 1 /*2013-01-27, 2013-01-27, EUR*/);
    const REAL amount = (1000. * (1. + fmin(((underlyings(366,1) / 11840.) - 1.), fmin(((underlyings(366,2) / 1200.) - 1.), ((underlyings(366,0) / 3758.05) - 1.)))));
	trajectory_inner( num_cash_flows, model_num, 1, amount, md_discts, vhat );

	return;
}

#endif

#define CONTRACT_NUM 2

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
//#if   (CONTRACT_NUM == 1)
       	trajectory_contract1(model_num, num_under, num_cash_flows, num_pricers, md_discts, md_detvals, inst_traj, vhat);
    } else if (contr_num == 2) {
//#elif (CONTRACT_NUM == 2)
        trajectory_contract2(model_num, num_under, num_cash_flows, num_pricers, md_discts, md_detvals, inst_traj, vhat);
    } else if (contr_num == 3) {
//#elif (CONTRACT_NUM == 3)
        trajectory_contract3(model_num, num_under, num_cash_flows, num_pricers, md_discts, md_detvals, inst_traj, vhat);
    } else {
//#else
        trajectory_contract1(model_num, num_under, num_cash_flows, num_pricers, md_discts, md_detvals, inst_traj, vhat);
    }
//#endif

}

#endif //CONTRACTS_CODE

