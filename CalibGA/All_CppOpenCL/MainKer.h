#ifndef MAIN_KERNEL
#define MAIN_KERNEL

#include "SDK_stub.h"
#include "Constants.h"

class MainKer {
    private:
        size_t  LWG[3];
        size_t  GWG[3];

        CpuArrays          cpu_arrs;
        OclBuffers         gpu_buffs;
        OclObjects         ocl_objs;


        // resources owned by this kernel
        cl_kernel    ker;
        cl_kernel    ker_red1;
        cl_kernel    ker_red2;

    public:
        MainKer ( const CpuArrays &   arrs,
                  const OclBuffers&   buffs,
                  const OclObjects&   objs
        ) : cpu_arrs(arrs), gpu_buffs(buffs), ocl_objs(objs)
        { 
//            LWG[0] = LWG_EG;  LWG[1] = 1;        LWG[2] = 1;
//            GWG[0] = arrs.SS; GWG[1] = POP_SIZE; GWG[2] = 1;
            
            ker      = mkKernel1 ();  
            ker_red1 = mkKernel2 ();  
            ker_red2 = mkKernel3 ();       
        }

        virtual ~MainKer() {
            clReleaseKernel(ker);
            clReleaseKernel(ker_red1);
            clReleaseKernel(ker_red2);
        }

        void run ( ) {
            cl_int ciErr1;

            //printf("Before run kernel 1 (LWG,GWG,UNROLL): (%ld,%ld,%d)\n", LWG, GWG, UNROLL);
#if (GPU_VERSION == 1)
            copy_in_buffers1();
#endif
            { // run the main, big kernel with local segmented reduction for each local group
                LWG[0] = LWG_EG;      LWG[1] = 1;        LWG[2] = 1;
                GWG[0] = cpu_arrs.SS; GWG[1] = POP_SIZE; GWG[2] = 1;

                ciErr1 = clEnqueueNDRangeKernel (  
                                ocl_objs.getCommandQueue(), 
                                ker, 3, NULL, this->GWG, this->LWG, 
                                0, NULL, NULL          
                            );
                ciErr1 |= clFinish( ocl_objs.getCommandQueue() );            
                oclCheckError(ciErr1, CL_SUCCESS);
            }
    
            { // run the second segmented reduction kernel
                LWG[0] = LWG_EG;      LWG[1] = 1;           LWG[2] = 1;
                GWG[0] = cpu_arrs.SS; GWG[1] = NUM_HERMITE; GWG[2] = POP_SIZE;
                
                ciErr1 = clEnqueueNDRangeKernel (  
                                ocl_objs.getCommandQueue(), 
                                ker_red1, 3, NULL, this->GWG, this->LWG, 
                                0, NULL, NULL          
                            );
                ciErr1 |= clFinish( ocl_objs.getCommandQueue() );            
                oclCheckError(ciErr1, CL_SUCCESS);
            }

            { // run the third regular reduction kernel
                LWG[0] = LWG_EG;      LWG[1] = 1;         LWG[2] = 1;
                GWG[0] = LWG_EG;      GWG[1] = POP_SIZE;  GWG[2] = 1;
                
                ciErr1 = clEnqueueNDRangeKernel (  
                                ocl_objs.getCommandQueue(), 
                                ker_red2, 3, NULL, this->GWG, this->LWG, 
                                0, NULL, NULL          
                            );
                ciErr1 |= clFinish( ocl_objs.getCommandQueue() );            
                oclCheckError(ciErr1, CL_SUCCESS);
            }

#if (GPU_VERSION == 1)
//            copy_out_buffers1();
//            copy_out_buffers2();
            copy_out_buffers3();
#endif
        }

    private:
        cl_kernel mkKernel1 ( ) {
            cl_kernel kernel;
            UINT      counter = 0; 
            cl_int    ciErr1;

            kernel = clCreateKernel( ocl_objs.program, "eval_genome_main", &ciErr1 );
            oclCheckError( ciErr1, CL_SUCCESS );

            counter = 0; 
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.shape          );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.swap_quotes    );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.genomes        );

            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.ci_t1cs_scale  );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.new_quote_price);
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.scalars        );

            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &NUM_SWAP_QUOTES         );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &POP_SIZE                );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &cpu_arrs.SS             );
            oclCheckError( ciErr1, CL_SUCCESS );

            return kernel;
        }

        cl_kernel mkKernel2 ( ) {
            cl_kernel kernel;
            UINT      counter = 0; 
            cl_int    ciErr1;

            kernel = clCreateKernel( ocl_objs.program, "eval_genome_red1", &ciErr1 );
            oclCheckError( ciErr1, CL_SUCCESS );

            counter = 0; 
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.shape          );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.gauss_coefs    );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.gauss_weights  );


            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.ci_t1cs_scale  );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.scalars        );

            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.accum0         );

            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &NUM_SWAP_QUOTES         );
//            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &POP_SIZE                );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &cpu_arrs.SS             );
            oclCheckError( ciErr1, CL_SUCCESS );

            return kernel;
        }

        cl_kernel mkKernel3 ( ) {
            cl_kernel kernel;
            UINT      counter = 0; 
            cl_int    ciErr1;

            kernel = clCreateKernel( ocl_objs.program, "eval_genome_red2", &ciErr1 );
            oclCheckError( ciErr1, CL_SUCCESS );

            counter = 0; 
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.genomes        );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.scalars        );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.accum0         );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&gpu_buffs.new_quote_price);

            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &NUM_SWAP_QUOTES          );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &NUM_HERMITE              );
            oclCheckError( ciErr1, CL_SUCCESS );

            return kernel;
        }


        void copy_in_buffers1() {
            cl_int ciErr;
            size_t cur_size;
            cl_command_queue& cmd_queue = ocl_objs.getCommandQueue();

            cur_size = 13 * POP_SIZE * sizeof(REAL);
            ciErr |= clEnqueueWriteBuffer(cmd_queue, gpu_buffs.genomes, CL_TRUE, 0, cur_size, cpu_arrs.genomes, 0, NULL, NULL);

            oclCheckError(ciErr, CL_SUCCESS);
        }

        void copy_out_buffers1() {
            cl_int ciErr;
            size_t cur_size;
            cl_command_queue& cmd_queue = ocl_objs.getCommandQueue();

            cur_size = 4 * cpu_arrs.SS * POP_SIZE * sizeof(REAL);
            ciErr |= clEnqueueReadBuffer(cmd_queue, gpu_buffs.ci_t1cs_scale, CL_TRUE, 0, cur_size, cpu_arrs.ci_t1cs_scale, 0, NULL, NULL);

            cur_size = NUM_SWAP_QUOTES * POP_SIZE * sizeof(REAL); // need only to write back the quotes (not the prices)
            ciErr |= clEnqueueReadBuffer(cmd_queue, gpu_buffs.new_quote_price, CL_TRUE, 0, cur_size, cpu_arrs.new_quote_price, 0, NULL, NULL);

            cur_size = 8 * NUM_SWAP_QUOTES * POP_SIZE * sizeof(REAL);
            ciErr |= clEnqueueReadBuffer(cmd_queue, gpu_buffs.scalars, CL_TRUE, 0, cur_size, cpu_arrs.scalars, 0, NULL, NULL);
        }

        void copy_out_buffers2() {
            cl_int ciErr;
            size_t cur_size;
            cl_command_queue& cmd_queue = ocl_objs.getCommandQueue();

            cur_size = NUM_SWAP_QUOTES * NUM_HERMITE * POP_SIZE * sizeof(REAL);
            ciErr |= clEnqueueReadBuffer(cmd_queue, gpu_buffs.accum0, CL_TRUE, 0, cur_size, cpu_arrs.accum0, 0, NULL, NULL);
        }

        void copy_out_buffers3() {
            cl_int ciErr;
            size_t cur_size, offset;
            cl_command_queue& cmd_queue = ocl_objs.getCommandQueue();

            cur_size =      POP_SIZE * sizeof(REAL);
            offset   = 11 * POP_SIZE * sizeof(REAL); 
            ciErr |= clEnqueueReadBuffer(cmd_queue, gpu_buffs.genomes, CL_TRUE, offset, cur_size, cpu_arrs.get_logLik() + POP_SIZE, 0, NULL, NULL);

//            cur_size = 13 * POP_SIZE * sizeof(REAL);
//            ciErr |= clEnqueueReadBuffer(cmd_queue, gpu_buffs.genomes, CL_TRUE, 0, cur_size, cpu_arrs.genomes, 0, NULL, NULL);

//            cur_size = 2 * NUM_SWAP_QUOTES * POP_SIZE * sizeof(REAL); // need only to write back the quotes (not the prices)
//            ciErr |= clEnqueueReadBuffer(cmd_queue, gpu_buffs.new_quote_price, CL_TRUE, 0, cur_size, cpu_arrs.new_quote_price, 0, NULL, NULL);
        }
};

#endif // MAIN_KERNEL

