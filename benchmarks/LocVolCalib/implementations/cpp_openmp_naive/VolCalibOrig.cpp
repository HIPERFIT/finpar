#include <vector>
#include <cmath>

#define WORKGROUP_SIZE  512 

#include "real.h"
#include "Util.h"
#include "ParseInput.h"

using namespace std;

struct PrivGlobs {

    //	grid
    vector<real_t>			myX;
    vector<real_t>          myY; 
    vector<real_t>          myTimeline;
    unsigned				myXindex; 
    unsigned                myYindex;

    //	variable
    vector<vector<real_t> > myResult;

    //	coeffs
    vector<vector<real_t> > myMuX;
    vector<vector<real_t> > myVarX;
    vector<vector<real_t> > myMuY;
    vector<vector<real_t> > myVarY;

    //	operators
    vector<vector<real_t> >	myDx;
    vector<vector<real_t> > myDxx;

    vector<vector<real_t> > myDy; 
    vector<vector<real_t> > myDyy;
} __attribute__ ((aligned (128)));





/***********************************/

void updateParams(const unsigned g, const real_t alpha, const real_t beta, const real_t nu, PrivGlobs& globs)
{
    for(unsigned i=0;i<globs.myX.size();++i)
        for(unsigned j=0;j<globs.myY.size();++j) {
            globs.myMuX[i][j]  = 0.0;
            globs.myVarX[i][j] = exp(2*(    beta*log(globs.myX[i])   
                                          + globs.myY[j]             
                                          - 0.5*nu*nu*globs.myTimeline[g] )
                                    );
            globs.myMuY[i][j]  = 0.0;
            globs.myVarY[i][j] = nu*nu;
        }
}


void initGrid(  const real_t s0, const real_t alpha, const real_t nu,const real_t t, 
                const unsigned numX, const unsigned numY, const unsigned numT, PrivGlobs& globs   
) {
    globs.myX.resize(numX);
    globs.myY.resize(numY);
    globs.myTimeline.resize(numT);

    for(unsigned i=0;i<numT;++i)
        globs.myTimeline[i] = t*i/(numT-1);

    const real_t stdX = 20*alpha*s0*sqrt(t);
    const real_t dx = stdX/numX;
    globs.myXindex = static_cast<unsigned>(s0/dx);

    for(unsigned i=0;i<numX;++i)
        globs.myX[i] = i*dx - globs.myXindex*dx + s0;

    const real_t stdY = 10*nu*sqrt(t);
    const real_t dy = stdY/numY;
    const real_t logAlpha = log(alpha);
    globs.myYindex = numY/2;

    for(unsigned i=0;i<numY;++i)
        globs.myY[i] = i*dy - globs.myYindex*dy + logAlpha;


    globs.myMuX.resize(numX);
    globs.myVarX.resize(numX);
    globs.myMuY.resize(numX);
    globs.myVarY.resize(numX);
    for(unsigned i=0;i<numX;++i) {
        globs.myMuX[i].resize(numY);
        globs.myVarX[i].resize(numY);
        globs.myMuY[i].resize(numY);
        globs.myVarY[i].resize(numY);
    }
}

void initOperator(const vector<real_t>& x, vector<vector<real_t> >& Dx, vector<vector<real_t> >& Dxx)
{
	const unsigned n = x.size();

	Dx.resize(n);
	Dxx.resize(n);

	for(unsigned i=0;i<n;++i)
	{
		Dx[i].resize(3);
		Dxx[i].resize(3);
	}

	real_t dxl, dxu;

	//	lower boundary
	dxl		 =  0.0;
	dxu		 =  x[1] - x[0];

	Dx[0][0]  =  0.0;
	Dx[0][1]  = -1.0/dxu;
	Dx[0][2]  =  1.0/dxu;
	
	Dxx[0][0] =  0.0;
	Dxx[0][1] =  0.0;
	Dxx[0][2] =  0.0;
	
	//	standard case
	for(unsigned i=1;i<n-1;i++)
	{
		dxl      = x[i]   - x[i-1];
		dxu      = x[i+1] - x[i];

		Dx[i][0]  = -dxu/dxl/(dxl+dxu);
		Dx[i][1]  = (dxu/dxl - dxl/dxu)/(dxl+dxu);
		Dx[i][2]  =  dxl/dxu/(dxl+dxu);

		Dxx[i][0] =  2.0/dxl/(dxl+dxu);
		Dxx[i][1] = -2.0*(1.0/dxl + 1.0/dxu)/(dxl+dxu);
		Dxx[i][2] =  2.0/dxu/(dxl+dxu); 
	}

	//	upper boundary
	dxl		   =  x[n-1] - x[n-2];
	dxu		   =  0.0;

	Dx[n-1][0]  = -1.0/dxl;
	Dx[n-1][1]  =  1.0/dxl;
	Dx[n-1][2]  =  0.0;

	Dxx[n-1][0] = 0.0;
	Dxx[n-1][1] = 0.0;
	Dxx[n-1][2] = 0.0;
}


void setPayoff(const real_t strike, PrivGlobs& globs )
{
    globs.myResult.resize(globs.myX.size());
	for(unsigned i=0;i<globs.myX.size();++i)
	{
		globs.myResult[i].resize(globs.myY.size());
		real_t payoff = max(globs.myX[i]-strike,(real_t)0.0);
		for(unsigned j=0;j<globs.myY.size();++j)
			globs.myResult[i][j] = payoff;
	}
}

inline void tridag(
    const vector<real_t>&   a,
    const vector<real_t>&   b,
    const vector<real_t>&   c,
    const vector<real_t>&   r,
    const int               n,
          vector<real_t>&   u,
          vector<real_t>&   uu)
{
    int    i, offset;
    real_t beta;

    u[0]  = r[0];
    uu[0] = b[0];

    for(i=1; i<n; i++) {
        beta  = a[i] / uu[i-1];

        uu[i] = b[i] - beta*c[i-1];
        u[i]  = r[i] - beta*u[i-1];
    }

    u[n-1] = u[n-1] / uu[n-1];

    for(i=n-2; i>=0; i--) {
        u[i] = (u[i] - c[i]*u[i+1]) / uu[i];
    }
}


void
rollback( const unsigned g, PrivGlobs& globs ) {
    unsigned numX = globs.myX.size(),
             numY = globs.myY.size();

    unsigned numZ = max(numX,numY);

    int k, l;
    unsigned i, j;

    int kl, ku, ll, lu;

    real_t dtInv = 1.0/(globs.myTimeline[g+1]-globs.myTimeline[g]);

    vector<vector<real_t> > u(numY,vector<real_t>(numX)), v(numX,vector<real_t>(numY));
    vector<real_t> a(numZ), b(numZ), c(numZ), y(numZ), yy(numZ);

    //	explicit x
    for(i=0;i<numX;i++) {
        kl =	 1*(i==0);
        ku = 2 - 1*(i==numX-1);
        for(j=0;j<numY;j++) {
            u[j][i] = dtInv*globs.myResult[i][j];
            for(k=kl;k<=ku;k++)
                u[j][i] += 0.5*( globs.myMuX[i][j]*globs.myDx[i][k] + 0.5*globs.myVarX[i][j]*globs.myDxx[i][k] )
                              * globs.myResult[i+k-1][j];
        }
    }

    //	explicit y
    for(j=0;j<numY;j++)
    {
        ll =	 1*(j==0);
        lu = 2 - 1*(j==numY-1);
        for(i=0;i<numX;i++) {
            v[i][j] = 0.0;
            for(l=ll;l<=lu;l++) {
                v[i][j] +=   (globs.myMuY[i][j]*globs.myDy[j][l] + 0.5*globs.myVarY[i][j]*globs.myDyy[j][l])
                            *  globs.myResult[i][j+l-1];
            }
            u[j][i] += v[i][j]; 
        }
    }

    //	implicit x
    for(j=0;j<numY;j++) {
        for(i=0;i<numX;i++) {
            a[i] =		 - 0.5*(globs.myMuX[i][j]*globs.myDx[i][0] + 0.5*globs.myVarX[i][j]*globs.myDxx[i][0]);
            b[i] = dtInv - 0.5*(globs.myMuX[i][j]*globs.myDx[i][1] + 0.5*globs.myVarX[i][j]*globs.myDxx[i][1]);
            c[i] =		 - 0.5*(globs.myMuX[i][j]*globs.myDx[i][2] + 0.5*globs.myVarX[i][j]*globs.myDxx[i][2]);
        }

        tridag(a,b,c,u[j],numX,u[j],yy);
    }

    //	implicit y
    for(i=0;i<numX;i++) {
        for(j=0;j<numY;j++) {
            a[j] =		 - 0.5*(globs.myMuY[i][j]*globs.myDy[j][0] + 0.5*globs.myVarY[i][j]*globs.myDyy[j][0]);
            b[j] = dtInv - 0.5*(globs.myMuY[i][j]*globs.myDy[j][1] + 0.5*globs.myVarY[i][j]*globs.myDyy[j][1]);
            c[j] =		 - 0.5*(globs.myMuY[i][j]*globs.myDy[j][2] + 0.5*globs.myVarY[i][j]*globs.myDyy[j][2]);
        }

        for(j=0;j<numY;j++)
            y[j] = dtInv*u[j][i] - 0.5*v[i][j];

        tridag(a,b,c,y,numY,globs.myResult[i],yy);
    }
}

real_t value(   const real_t s0,
                const real_t strike, 
                const real_t t, 
                const real_t alpha, 
                const real_t nu, 
                const real_t beta,
                const unsigned int numX,
                const unsigned int numY,
                const unsigned int numT
) {	
    PrivGlobs   globs;
	initGrid(s0,alpha,nu,t, numX, numY, numT, globs);
	initOperator(globs.myX,globs.myDx,globs.myDxx);
	initOperator(globs.myY,globs.myDy,globs.myDyy);

	setPayoff(strike, globs);
	for(int i = globs.myTimeline.size()-2;i>=0;--i)
	{
		updateParams(i,alpha,beta,nu,globs);
		rollback(i, globs);
	}

	return globs.myResult[globs.myXindex][globs.myYindex];
}

real_t* run_CPUkernel(  
                const unsigned int&     outer,
                const unsigned int&     numX,
                const unsigned int&     numY,
                const unsigned int&     numT,
                const real_t&           s0,
                const real_t&           t, 
                const real_t&           alpha, 
                const real_t&           nu, 
                const real_t&           beta,
                const vector<real_t>&   strikes,
                      vector<real_t>&   res
) {
#pragma omp parallel for default(shared) schedule(static) if(outer>4)
    for( unsigned i = 0; i < outer; ++ i ) {
        res[i] = value( s0, strikes[i], t, alpha, nu, beta,
                        numX, numY, numT );
    }
}


int main()
{
    unsigned int OUTER_LOOP_COUNT, NUM_X, NUM_Y, NUM_T; 
	const real_t s0 = 0.03, strike = 0.03, t = 5.0, alpha = 0.2, nu = 0.6, beta = 0.5;

    cout<<"\n// Running Original (Parallel) Volatility-Calibration Benchmark"<<endl;

    readDataSet( OUTER_LOOP_COUNT, NUM_X, NUM_Y, NUM_T ); 

	vector<real_t> strikes(OUTER_LOOP_COUNT),res(OUTER_LOOP_COUNT);

	for(unsigned i=0;i<OUTER_LOOP_COUNT;++i)
		strikes[i] = 0.001*i;

    const int Ps = get_CPU_num_threads();

        unsigned long int elapsed_usec = 0;
    {   // Main Computational Kernel
        struct timeval t_start, t_end, t_diff;
        gettimeofday(&t_start, NULL);

        run_CPUkernel( OUTER_LOOP_COUNT, NUM_X, NUM_Y, NUM_T, s0, t, alpha, nu, beta, strikes, res );

        gettimeofday(&t_end, NULL);
        timeval_subtract(&t_diff, &t_end, &t_start);
        elapsed_usec = t_diff.tv_sec*1e6+t_diff.tv_usec;
    }

    {
      FILE* runtime = fopen(getenv("HIPERMARK_RUNTIME"), "w");
      FILE* result = fopen(getenv("HIPERMARK_RESULT"), "w");
      fprintf(runtime, "%d\n", elapsed_usec / 1000);
      fclose(runtime);
      write_1Darr(result, res.data(), OUTER_LOOP_COUNT);
      fclose(result);
    }

    return 0;
}

