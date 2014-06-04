#include "../Optimizations.h" 

#ifdef _OPTIMIZATION_MEM_COALES_ON
  #include "../GenericPricingVectOpt.cl"
#else
  #include "../GenericPricingVectUncoalesced.cl"
#endif
      
inline void trajectory_contract( // SIMPLE                     
		UINT                      model_num,     
		__constant LoopROScalars* ro_scal, 
        __constant REAL*          pc_coefs, 
        __global    REAL*         inst_traj,
        __global    REAL*         vhat
) {  
    __constant REAL*  model_deter_val = pc_coefs + ro_scal->pc_deter_val_beg;
    __constant REAL*  model_discounts = pc_coefs + ro_scal->pc_discounts_beg;

	REAL amount = fmax(0.0, ((underlyings(0,0)) - 4000.0) * model_deter_val[model_num*ro_scal->num_det_pricers + 0]);

	trajectory_inner( model_num, ro_scal, model_discounts, vhat, 0, amount, 0 );
}
                
