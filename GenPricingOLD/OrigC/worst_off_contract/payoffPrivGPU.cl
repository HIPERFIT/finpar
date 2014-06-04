#include "../GenericPricingPrivOpt.cl"      
#include "../Optimizations.h"
                   
                 
inline void trajectory_contract( // SIMPLE              
		UINT                        model_num,                
		__constant LoopROScalars*   ro_scal,    
        __constant REAL*            pc_coefs, 
        UINT                        block_size,   
        REAL*                       inst_traj, 
        __local    REAL*            local_vhat    
) {                                    
	__constant REAL*  model_deter_val = pc_coefs + ro_scal->pc_deter_val_beg;     
	__constant REAL*  model_discounts = pc_coefs + ro_scal->pc_discounts_beg;  
      
	REAL x50;                              
        
	if ((1. <= fmin((underlyings(0,1) / 11840.), fmin((underlyings(0,2) / 1200.), (underlyings(0,0) / 3758.05))))) goto L21; 
	if ((1. <= fmin((underlyings(1,1) / 11840.), fmin((underlyings(1,2) / 1200.), (underlyings(1,0) / 3758.05))))) goto L22;
	if ((1. <= fmin((underlyings(2,1) / 11840.), fmin((underlyings(2,2) / 1200.), (underlyings(2,0) / 3758.05))))) goto L23;
	if ((1. <= fmin((underlyings(3,1) / 11840.), fmin((underlyings(3,2) / 1200.), (underlyings(3,0) / 3758.05))))) goto L24;
	x50=fmin((underlyings(4,1) / 11840.), fmin((underlyings(4,2) / 1200.), (underlyings(4,0) / 3758.05)));
      
	//model->notify_cash_flow(model, 0, 1000., 4 /*2017-02-03, 2017-01-27, EUR*/); 
	trajectory_inner( model_num, block_size, ro_scal, model_discounts, local_vhat, 0, 1000., 4 );  

	if ((1. <= x50)) goto L25; 
	if ((0.75 < x50)) return;             
            
    
	//model->notify_cash_flow(model, 0, (-(1000. * (1. - x50))), 4 /*2017-02-03, 2017-01-27, EUR*/); return;
    trajectory_inner( model_num, block_size, ro_scal, model_discounts, local_vhat, 0, (-(1000. * (1. - x50))), 4 ); return;
L25:    
	//model->notify_cash_flow(model, 0, 750., 4 /*2017-02-03, 2017-01-27, EUR*/); return;
    trajectory_inner( model_num, block_size, ro_scal, model_discounts, local_vhat,  0, 750., 4 ); return;
L24:  
	//model->notify_cash_flow(model, 0, 1600., 3 /*2016-02-03, 2016-01-27, EUR*/); return;
    trajectory_inner( model_num, block_size, ro_scal, model_discounts, local_vhat,  0, 1600., 3 ); return;
L23:
	//model->notify_cash_flow(model, 0, 1450., 2 /*2015-02-03, 2015-01-27, EUR*/); return;
    trajectory_inner( model_num, block_size, ro_scal, model_discounts, local_vhat,  0, 1450., 2 ); return;
L22:
	//model->notify_cash_flow(model, 0, 1300., 1 /*2014-02-03, 2014-01-27, EUR*/); return;
    trajectory_inner( model_num, block_size, ro_scal, model_discounts, local_vhat,  0, 1300., 1 ); return;
L21:
	//model->notify_cash_flow(model, 0, 1150., 0 /*2013-02-01, 2013-01-27, EUR*/); return;
    trajectory_inner( model_num, block_size, ro_scal, model_discounts, local_vhat, 0, 1150., 0  ); return;
}




