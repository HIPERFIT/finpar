#include <vector>
#include <cmath>

#define WORKGROUP_SIZE  512

#ifdef REAL_TYPE
typedef REAL_TYPE real_t;
#else
typedef double real_t;
#endif

#include "Util.h"
#include "ParseInput.h"

using namespace std;

//	grid
vector<real_t>			myX, myY, myTimeline;
unsigned				myXindex, myYindex;

//	variable
vector<vector<real_t> > myResult;

//	coeffs
vector<vector<real_t> > myMuX, myVarX, myMuY, myVarY;

//	operators
vector<vector<real_t> >	myDx, myDxx, myDy, myDyy;




/***********************************/

void updateParams(const unsigned g, const real_t alpha, const real_t beta, const real_t nu)
{
	for(unsigned i=0;i<myX.size();++i)
		for(unsigned j=0;j<myY.size();++j)
		{
			myMuX[i][j] = 0.0;
			myVarX[i][j] = exp(2*(beta*log(myX[i]) + myY[j] - 0.5*nu*nu*myTimeline[g]));
			myMuY[i][j] = 0.0;
			myVarY[i][j] = nu*nu;
		}
}


void initGrid(const real_t s0, const real_t alpha, const real_t nu,const real_t t, const unsigned numX, const unsigned numY, const unsigned numT)
{
	myX.resize(numX);
	myY.resize(numY);
	myTimeline.resize(numT);

	for(unsigned i=0;i<numT;++i)
		myTimeline[i] = t*i/(numT-1);

	const real_t stdX = 20*alpha*s0*sqrt(t);
	const real_t dx = stdX/numX;
	myXindex = static_cast<unsigned>(s0/dx);

	for(unsigned i=0;i<numX;++i)
		myX[i] = i*dx - myXindex*dx + s0;

	const real_t stdY = 10*nu*sqrt(t);
	const real_t dy = stdY/numY;
	const real_t logAlpha = log(alpha);
	myYindex = numY/2;

	for(unsigned i=0;i<numY;++i)
		myY[i] = i*dy - myYindex*dy + logAlpha;


	myMuX.resize(numX);
	myVarX.resize(numX);
	myMuY.resize(numX);
	myVarY.resize(numX);
	for(unsigned i=0;i<numX;++i)
	{
		myMuX[i].resize(numY);
		myVarX[i].resize(numY);
		myMuY[i].resize(numY);
		myVarY[i].resize(numY);
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


void setPayoff(const real_t strike)
{
	myResult.resize(myX.size());
	for(unsigned i=0;i<myX.size();++i)
	{
		myResult[i].resize(myY.size());
		real_t payoff = max(myX[i]-strike,(real_t)0.0);
		for(unsigned j=0;j<myY.size();++j)
			myResult[i][j] = payoff;
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
    real_t   beta;

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
rollback(const unsigned g)
{
	unsigned numX = myX.size(),
			 numY = myY.size();

	unsigned numZ = max(numX,numY);

	int k, l;
	unsigned i, j;

	int kl, ku, ll, lu;

	real_t dtInv = 1.0/(myTimeline[g+1]-myTimeline[g]);

	vector<vector<real_t> > u(numY,vector<real_t>(numX)), v(numX,vector<real_t>(numY));

	vector<real_t> a(numZ), b(numZ), c(numZ), y(numZ), yy(numZ);

	//	explicit x
	for(i=0;i<numX;i++)
	{
		kl =	 1*(i==0);
		ku = 2 - 1*(i==numX-1);
		for(j=0;j<numY;j++)
		{
			u[j][i] = dtInv*myResult[i][j];
			for(k=kl;k<=ku;k++)
				u[j][i] += 0.5*(myMuX[i][j]*myDx[i][k] + 0.5*myVarX[i][j]*myDxx[i][k])*myResult[i+k-1][j];
		}
	}

	//	explicit y
	for(j=0;j<numY;j++)
	{
		ll =	 1*(j==0);
		lu = 2 - 1*(j==numY-1);
		for(i=0;i<numX;i++)
		{
			v[i][j] = 0.0;
			for(l=ll;l<=lu;l++)
				v[i][j] +=     (myMuY[i][j]*myDy[j][l] + 0.5*myVarY[i][j]*myDyy[j][l])*myResult[i][j+l-1];

			u[j][i] += v[i][j]; 
		}
	}

	//	implicit x
	for(j=0;j<numY;j++)
	{
		for(i=0;i<numX;i++)
		{
			a[i] =		 - 0.5*(myMuX[i][j]*myDx[i][0] + 0.5*myVarX[i][j]*myDxx[i][0]);
			b[i] = dtInv - 0.5*(myMuX[i][j]*myDx[i][1] + 0.5*myVarX[i][j]*myDxx[i][1]);
			c[i] =		 - 0.5*(myMuX[i][j]*myDx[i][2] + 0.5*myVarX[i][j]*myDxx[i][2]);
		}

		tridag(a,b,c,u[j],numX,u[j],yy);
	}

	//	implicit y
	for(i=0;i<numX;i++)
	{
		for(j=0;j<numY;j++)
		{
			a[j] =		 - 0.5*(myMuY[i][j]*myDy[j][0] + 0.5*myVarY[i][j]*myDyy[j][0]);
			b[j] = dtInv - 0.5*(myMuY[i][j]*myDy[j][1] + 0.5*myVarY[i][j]*myDyy[j][1]);
			c[j] =		 - 0.5*(myMuY[i][j]*myDy[j][2] + 0.5*myVarY[i][j]*myDyy[j][2]);
		}

		for(j=0;j<numY;j++)
			y[j] = dtInv*u[j][i] - 0.5*v[i][j];

		tridag(a,b,c,y,numY,myResult[i],yy);
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
	initGrid(s0,alpha,nu,t, numX, numY, numT);
	initOperator(myX,myDx,myDxx);
	initOperator(myY,myDy,myDyy);

	setPayoff(strike);
	for(int i = myTimeline.size()-2;i>=0;--i)
	{
		updateParams(i,alpha,beta,nu);
		rollback(i);
	}

	return myResult[myXindex][myYindex];
}

int main()
{
    unsigned int OUTER_LOOP_COUNT, NUM_X, NUM_Y, NUM_T; 
	const real_t s0 = 0.03, strike = 0.03, t = 5.0, alpha = 0.2, nu = 0.6, beta = 0.5;

    fprintf(stdout, "\n// Original (Sequential) Volatility Calibration Benchmark:\n");
    readDataSet( OUTER_LOOP_COUNT, NUM_X, NUM_Y, NUM_T ); 

	vector<real_t> strikes(OUTER_LOOP_COUNT),res(OUTER_LOOP_COUNT);

	for(unsigned i=0;i<OUTER_LOOP_COUNT;++i)
		strikes[i] = 0.001*i;

    unsigned long int elapsed_usec = 0;
    {   // Main Computational Kernel
        struct timeval t_start, t_end, t_diff;
        gettimeofday(&t_start, NULL);

    	for(unsigned i=0;i<OUTER_LOOP_COUNT;++i) {
	    	res[i] = value( s0, strikes[i], t, alpha, nu, beta,
                            NUM_X, NUM_Y, NUM_T );
        }

        gettimeofday(&t_end, NULL);
        timeval_subtract(&t_diff, &t_end, &t_start);
        elapsed_usec = t_diff.tv_sec*1e6+t_diff.tv_usec;
    }

    {   FILE* runtime = fopen(getenv("HIPERMARK_RUNTIME"), "w");
      FILE* result = fopen(getenv("HIPERMARK_RESULT"), "w");
        fprintf(runtime, "%d\n", elapsed_usec / 1000);
        fclose(runtime);
        write_1Darr(result, res.data(), OUTER_LOOP_COUNT);
        fclose(result);
    }
	return 0;
}