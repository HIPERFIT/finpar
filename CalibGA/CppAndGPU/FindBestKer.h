#ifndef FIND_BEST_KERNEL
#define FIND_BEST_KERNEL

#include "SDK_stub.h"
#include "Constants.h"

class FindBestKer {
    // Assumptions: 
    //   1. the input  array, `arr' has `size(arr) >= N',
    //   2. the result arrays, has `size(r_ind)=size(r_val) >= LWG'
    //
    // Computation is split in two stages:
    //   - one that compute at most LWG elements (denoting partial reduces)
    //   - a second one that reduces the resulting LWG elements in one local workgroup.
    private:
        const size_t       LWG;
              size_t       GWG;
              uint         UNROLL;

        // inputs, i.e., resources not owned by this kernel!
        const cl_mem       arr;  // of size `N', and `offset' into it!
        const uint         N;
        const uint         offset;


        OclObjects         ocl;

        // res_ind and res_val have at least size LWG_FB
        const cl_mem       res_ind;
        const cl_mem       res_val;


        // resources owned by this kernel
        cl_kernel    ker1; 
        cl_kernel    ker2;

    public:
        FindBestKer ( const uint        sz, 
                      const uint        off,
                      const cl_mem&     a,
                      const cl_mem&     r_ind,
                      const cl_mem&     r_val,
                      const OclObjects& o,
                      const int&        lwg
        ) : N(sz), offset(off), LWG(lwg),  
            arr(a), res_ind(r_ind), res_val(r_val), ocl(o)
        { 
            UNROLL = mkUnrollFactor();
            GWG    = mkGlobalWork  ();
            ker1   = mkKernel1     ();
            ker2   = mkKernel2     ();
            
            bool sanity = is_pow2(lwg) && lwg <= LWG_FB;
            assert( sanity && "In FindBestKer: Workgroup size not a pow of two or too large!" );
        }

        virtual ~FindBestKer() {
            clReleaseKernel(ker1);
            clReleaseKernel(ker2);
        }

        void run ( int& best_ind, REAL& best_val ) {
            cl_int ciErr1;

            //printf("Before run kernel 1 (LWG,GWG,UNROLL): (%ld,%ld,%d)\n", LWG, GWG, UNROLL);

            // run the reduction for each local group
            ciErr1 = clEnqueueNDRangeKernel (  
                            ocl.getCommandQueue(), 
                            ker1, 1, NULL, & GWG, & LWG, 
                            0, NULL, NULL          
                        );
            ciErr1 |= clFinish( ocl.getCommandQueue() );            
            oclCheckError(ciErr1, CL_SUCCESS);

            if( GWG > LWG ) {
                // run the reduction for the resulted local group
                ciErr1 = clEnqueueNDRangeKernel (  
                                ocl.getCommandQueue(), 
                                ker2, 1, NULL, & LWG, & LWG,
                                0, NULL, NULL          
                            );
                ciErr1 |= clFinish( ocl.getCommandQueue() );            
                oclCheckError(ciErr1, CL_SUCCESS);
            }

            { // Finally, read-back the histogram result!
                ciErr1 |= clEnqueueReadBuffer (
                               ocl.getCommandQueue(), res_val, CL_TRUE, 0, 
                               sizeof(REAL), &best_val, 0, NULL, NULL
                            );
                ciErr1 |= clEnqueueReadBuffer (
                               ocl.getCommandQueue(), res_ind, CL_TRUE, 0, 
                               sizeof(uint), &best_ind, 0, NULL, NULL
                            );
            }
            oclCheckError(ciErr1, CL_SUCCESS);
//            printf("After enqueue Best (Ind,Val): (%d,%f)\n\n", best_ind, best_val);
        }

    private:
        int mkUnrollFactor() {
            const int  sqLWG  = LWG * LWG;
            const bool exact1 = (N % sqLWG) == 0;
            const int  unroll = ( exact1 ) ? ( N / sqLWG ) : ( (N / sqLWG)   + 1 );
            return unroll;
        }

        int mkGlobalWork() {         
            const int  urolLWG= UNROLL * LWG;
            const bool exact2 = (N % urolLWG) == 0;
            const int  numGWG = ( exact2 ) ? (N / urolLWG) : ( (N / urolLWG) + 1 );
            return numGWG * LWG;
        }

        cl_kernel mkKernel1( ) {
            cl_kernel kernel;
            UINT      counter = 0; 
            cl_int    ciErr1;

            kernel = clCreateKernel( ocl.program, "redbest_ker1", &ciErr1 );
            oclCheckError( ciErr1, CL_SUCCESS );

            counter = 0; 
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&this->arr     );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&this->res_ind );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&this->res_val );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &this->offset  );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &this->N       );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &this->UNROLL  );

            oclCheckError( ciErr1, CL_SUCCESS );

            return kernel;
        }

        cl_kernel mkKernel2( ) {
            cl_kernel kernel;
            UINT      counter = 0; 
            cl_int    ciErr1;

            kernel = clCreateKernel( ocl.program, "redbest_ker2", &ciErr1 );
            oclCheckError( ciErr1, CL_SUCCESS );

            uint newN = GWG / LWG;

            counter = 0; 
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&this->res_ind );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(cl_mem), (void*)&this->res_val );
            ciErr1 |= clSetKernelArg( kernel, counter++, sizeof(uint),          &newN    );

            oclCheckError( ciErr1, CL_SUCCESS );

            return kernel;
        }

};

#endif // FIND_BEST_KERNEL

