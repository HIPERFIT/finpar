#include "../GenericPricingPrivOpt.cl"
   
inline void trajectory_contract( // SIMPLE                  
		UINT                      model_num,   
		__constant LoopROScalars* ro_scal,
        __constant REAL*          pc_coefs,
        UINT                      block_size,    
        REAL*                     inst_traj,
        __local    REAL*          local_vhat
) {     
    __constant REAL*  model_deter_val = pc_coefs + ro_scal->pc_deter_val_beg;
    __constant REAL*  model_discounts = pc_coefs + ro_scal->pc_discounts_beg;

	//__local REAL* inst_traj = inst_trajWF + model_num * (ro_scal->md_dim * ro_scal->md_nb_path_dates) * block_size;

	REAL amount = fmax(0.0, ((underlyings(0,0)) - 4000.0) * model_deter_val[model_num*ro_scal->num_det_pricers + 0]);

	trajectory_inner( model_num, block_size, ro_scal, model_discounts, local_vhat, 0, amount, 0 ); 
} 
  
