#ifndef GEN_ALG_KERNELS
#define GEN_ALG_KERNELS

#include "SDK_stub.h"
#include "Constants.h"

class GenAlgKers {
    private:
        size_t  LWG[3];
        size_t  GWG[3];

        const uint lwg;
//              uint sobol_offset;

        CpuArrays    cpu_arrs;
        OclBuffers   gpu_buffs;
        OclObjects   ocl_objs;

        // resources owned by this kernel
        cl_kernel     init_ker;
        cl_kernel    cplik_ker;
        cl_kernel    accnd_ker;
        cl_kernel    acprp_ker;
        cl_kernel    mutat_ker;
        cl_kernel     mcde_ker;

    public:
        GenAlgKers( const CpuArrays &   arrs,
                    const OclBuffers&   buffs,
                    const OclObjects&   objs
        ) : cpu_arrs(arrs), gpu_buffs(buffs), ocl_objs(objs), lwg( mkLWG(POP_SIZE) )
        { 
            //sobol_offset = 1;

            bool valid = (POP_SIZE >= 32) && (POP_SIZE <= 512); 
            assert( valid && "Population size NOT in range [32, 512]!" );
             init_ker = mkInitKernel       ();
            cplik_ker = mkCpLiKernel       ();
            accnd_ker = mkAcceptCondKernel ();
            acprp_ker = mkAcceptPropKernel ();
            mutat_ker = mkMutateDimKernel  ();
             mcde_ker = mkMcMcDcKernel     ();
        }

        virtual ~GenAlgKers() {
            clReleaseKernel( init_ker);
            clReleaseKernel(cplik_ker);
            clReleaseKernel(accnd_ker);
            clReleaseKernel(acprp_ker);
            clReleaseKernel(mutat_ker);
            clReleaseKernel( mcde_ker);
        }

        void run_init( ) {
            cl_int ciErr1;

            LWG[0] = lwg;    LWG[1] = 1;        LWG[2] = 1;
            GWG[0] = lwg;    GWG[1] = 6;        GWG[2] = 1;

            // set the index of the next random number!
            ciErr1 |= clSetKernelArg( init_ker, 3, sizeof(uint), &sobol_offset );
            oclCheckError(ciErr1, CL_SUCCESS);
            ciErr1 = clEnqueueNDRangeKernel (  
                                ocl_objs.getCommandQueue(), 
                                init_ker, 2, NULL, this->GWG, this->LWG, 
                                0, NULL, NULL          
                            );
            oclCheckError(ciErr1, CL_SUCCESS);
            ciErr1 |= clFinish( ocl_objs.getCommandQueue() );            
            oclCheckError(ciErr1, CL_SUCCESS);
            sobol_offset += (5 * POP_SIZE);
        }

        void run_cpLik( ) {
            cl_int ciErr1;

            LWG[0] = lwg;    LWG[1] = 1;        LWG[2] = 1;
            GWG[0] = lwg;    GWG[1] = 1;        GWG[2] = 1;

            ciErr1 = clEnqueueNDRangeKernel (  
                                ocl_objs.getCommandQueue(), 
                                cplik_ker, 1, NULL, this->GWG, this->LWG, 
                                0, NULL, NULL          
                            );
            ciErr1 |= clFinish( ocl_objs.getCommandQueue() );            
            oclCheckError(ciErr1, CL_SUCCESS);
        }

        void run_accept_cond( ) {
            cl_int ciErr1;

            LWG[0] = lwg;    LWG[1] = 1;        LWG[2] = 1;
            GWG[0] = lwg;    GWG[1] = 1;        GWG[2] = 1;

            // set the index of the next random number!
            ciErr1 |= clSetKernelArg( accnd_ker, 3, sizeof(uint), &sobol_offset);
            ciErr1 = clEnqueueNDRangeKernel (  
                                ocl_objs.getCommandQueue(), 
                                accnd_ker, 1, NULL, this->GWG, this->LWG, 
                                0, NULL, NULL          
                            );
            ciErr1 |= clFinish( ocl_objs.getCommandQueue() );            
            oclCheckError(ciErr1, CL_SUCCESS);

            sobol_offset += (1 * POP_SIZE);
        }

        void run_accept_prop(  ) {
            cl_int ciErr1;

            LWG[0] = lwg;    LWG[1] = 1;        LWG[2] = 1;
            GWG[0] = lwg;    GWG[1] = 6;        GWG[2] = 1;

            ciErr1 = clEnqueueNDRangeKernel (  
                                ocl_objs.getCommandQueue(), 
                                acprp_ker, 2, NULL, this->GWG, this->LWG, 
                                0, NULL, NULL          
                            );
            ciErr1 |= clFinish( ocl_objs.getCommandQueue() );            
            oclCheckError(ciErr1, CL_SUCCESS);
        }

        // IF dim_j > 4 then is a mutate_all!
        void run_mutate( const uint dim_j ) {
            cl_int ciErr1;

            LWG[0] = lwg;    LWG[1] = 1;        LWG[2] = 1;
            GWG[0] = lwg;    GWG[1] = 1;        GWG[2] = 1;

            // set the index of the next random number!
            ciErr1 |= clSetKernelArg( mutat_ker, 3, sizeof(uint), &sobol_offset);
            ciErr1 |= clSetKernelArg( mutat_ker, 4, sizeof(uint), &dim_j);

            ciErr1 = clEnqueueNDRangeKernel (  
                                ocl_objs.getCommandQueue(), 
                                mutat_ker, 1, NULL, this->GWG, this->LWG, 
                                0, NULL, NULL          
                            );
            ciErr1 |= clFinish( ocl_objs.getCommandQueue() );            
            oclCheckError(ciErr1, CL_SUCCESS);

            sobol_offset += (5 * POP_SIZE);
        }

        void run_McMcDc( ) {
            cl_int ciErr1;

            LWG[0] = lwg;    LWG[1] = 1;        LWG[2] = 1;
            GWG[0] = lwg;    GWG[1] = 1;        GWG[2] = 1;

            // set the index of the next random number!
            ciErr1 |= clSetKernelArg( mcde_ker, 3, sizeof(uint), &sobol_offset);

            ciErr1 = clEnqueueNDRangeKernel (  
                                ocl_objs.getCommandQueue(), 
                                mcde_ker, 1, NULL, this->GWG, this->LWG, 
                                0, NULL, NULL          
                            );
            ciErr1 |= clFinish( ocl_objs.getCommandQueue() );            
            oclCheckError(ciErr1, CL_SUCCESS);

            sobol_offset += (8 * POP_SIZE);
        }


    private:

        uint mkLWG( const uint Npop ) {
            if( (Npop % WARP) == 0 ) return Npop;
            else return (Npop / WARP + 1) * WARP;
        }

        cl_kernel mkInitKernel ( ) {
            cl_kernel kernel;
            UINT      counter = 0;
            cl_int    ciErr1;

            kernel = clCreateKernel( ocl_objs.program, "init_genome", &ciErr1 );
            oclCheckError( ciErr1, CL_SUCCESS );

            counter = 0; 
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.genomes        );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.gene_ranges    );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.sobol_dir_vct  );

            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &sobol_offset             );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &POP_SIZE                 );
            oclCheckError( ciErr1, CL_SUCCESS );

            return kernel;
        }

        cl_kernel mkCpLiKernel ( ) {
            cl_kernel kernel;
            UINT      counter = 0;
            cl_int    ciErr1;

            kernel = clCreateKernel( ocl_objs.program, "copyLogLik", &ciErr1 );
            oclCheckError( ciErr1, CL_SUCCESS );

            counter = 0; 
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.genomes        );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &POP_SIZE                 );
            oclCheckError( ciErr1, CL_SUCCESS );

            return kernel;
        }

        cl_kernel mkAcceptCondKernel ( ) {
            cl_kernel kernel;
            UINT      counter = 0;
            cl_int    ciErr1;

            kernel = clCreateKernel( ocl_objs.program, "accept_cond", &ciErr1 );
            oclCheckError( ciErr1, CL_SUCCESS );

            counter = 0; 
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.genomes        );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.sobol_dir_vct  );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.accept_cond    );

            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &sobol_offset             );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &POP_SIZE                 );
            oclCheckError( ciErr1, CL_SUCCESS );

            return kernel;
        }

        cl_kernel mkAcceptPropKernel ( ) {
            cl_kernel kernel;
            UINT      counter = 0;
            cl_int    ciErr1;

            kernel = clCreateKernel( ocl_objs.program, "accept_prop", &ciErr1 );
            oclCheckError( ciErr1, CL_SUCCESS );

            counter = 0; 
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.genomes        );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.accept_cond    );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &POP_SIZE                 );
            oclCheckError( ciErr1, CL_SUCCESS );

            return kernel;
        }

        cl_kernel mkMutateDimKernel ( ) {
            cl_kernel kernel;
            UINT      counter = 0;
            cl_int    ciErr1;
            uint      dim_j   = 7;
            

            kernel = clCreateKernel( ocl_objs.program, "mutate_dims", &ciErr1 );
            oclCheckError( ciErr1, CL_SUCCESS );

            counter = 0; 
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.genomes        );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.gene_ranges    );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.sobol_dir_vct  );

            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &sobol_offset             );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &dim_j                    );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(REAL),          &MOVES_UNIF_AMPL_RATIO    );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &POP_SIZE                 );
            oclCheckError( ciErr1, CL_SUCCESS );

            return kernel;
        }

        cl_kernel mkMcMcDcKernel ( ) {
            cl_kernel kernel;
            UINT      counter = 0;
            cl_int    ciErr1;

            const REAL  gamma_avg = 2.38 / sqrt(2.0*GENOME_DIM);
            const REAL  ampl_ratio = 0.1 * MOVES_UNIF_AMPL_RATIO;


            kernel = clCreateKernel( ocl_objs.program, "mcmc_DE", &ciErr1 );
            oclCheckError( ciErr1, CL_SUCCESS );

            counter = 0; 
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.genomes        );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.gene_ranges    );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.sobol_dir_vct  );

            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &sobol_offset             );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(REAL),          &gamma_avg                );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(REAL),          &ampl_ratio               );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &POP_SIZE                 );
            oclCheckError( ciErr1, CL_SUCCESS );

            return kernel;
        }
};

#endif
