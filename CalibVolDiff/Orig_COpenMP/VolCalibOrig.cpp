#include <vector>
#include <cmath>

typedef double REAL;
#define WORKGROUP_SIZE  512 

#include "Util.h"
#include "../includeC/ParseInput.h"

using namespace std;

struct PrivGlobs {

    //	grid
    vector<double>			myX;
    vector<double>          myY; 
    vector<double>          myTimeline;
    unsigned				myXindex; 
    unsigned                myYindex;

    //	variable
    vector<vector<double> > myResult;

    //	coeffs
    vector<vector<double> > myMuX;
    vector<vector<double> > myVarX;
    vector<vector<double> > myMuY;
    vector<vector<double> > myVarY;

    //	operators
    vector<vector<double> >	myDx;
    vector<vector<double> > myDxx;

    vector<vector<double> > myDy; 
    vector<vector<double> > myDyy;
} __attribute__ ((aligned (128)));





/***********************************/

void updateParams(const unsigned g, const double alpha, const double beta, const double nu, PrivGlobs& globs)
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


void initGrid(  const double s0, const double alpha, const double nu,const double t, 
                const unsigned numX, const unsigned numY, const unsigned numT, PrivGlobs& globs   
) {
    globs.myX.resize(numX);
    globs.myY.resize(numY);
    globs.myTimeline.resize(numT);

    for(unsigned i=0;i<numT;++i)
        globs.myTimeline[i] = t*i/(numT-1);

    const double stdX = 20*alpha*s0*sqrt(t);
    const double dx = stdX/numX;
    globs.myXindex = static_cast<unsigned>(s0/dx);

    for(unsigned i=0;i<numX;++i)
        globs.myX[i] = i*dx - globs.myXindex*dx + s0;

    const double stdY = 10*nu*sqrt(t);
    const double dy = stdY/numY;
    const double logAlpha = log(alpha);
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

void initOperator(const vector<double>& x, vector<vector<double> >& Dx, vector<vector<double> >& Dxx)
{
	const unsigned n = x.size();

	Dx.resize(n);
	Dxx.resize(n);

	for(unsigned i=0;i<n;++i)
	{
		Dx[i].resize(3);
		Dxx[i].resize(3);
	}

	double dxl, dxu;

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


void setPayoff(const double strike, PrivGlobs& globs )
{
    globs.myResult.resize(globs.myX.size());
	for(unsigned i=0;i<globs.myX.size();++i)
	{
		globs.myResult[i].resize(globs.myY.size());
		double payoff = max(globs.myX[i]-strike,0.0);
		for(unsigned j=0;j<globs.myY.size();++j)
			globs.myResult[i][j] = payoff;
	}
}

inline void tridag(
    const vector<double>&   a,
    const vector<double>&   b,
    const vector<double>&   c,
    const vector<double>&   r,
    const int               n,
          vector<double>&   u,
          vector<double>&   uu)
{
    int    i, offset;
    REAL   beta;

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

    double dtInv = 1.0/(globs.myTimeline[g+1]-globs.myTimeline[g]);

    vector<vector<double> > u(numY,vector<double>(numX)), v(numX,vector<double>(numY));
    vector<double> a(numZ), b(numZ), c(numZ), y(numZ), yy(numZ);

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

double value(   const double s0,
                const double strike, 
                const double t, 
                const double alpha, 
                const double nu, 
                const double beta,
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

double* run_CPUkernel(  
                const unsigned int&     outer,
                const unsigned int&     numX,
                const unsigned int&     numY,
                const unsigned int&     numT,
                const double&           s0,
                const double&           t, 
                const double&           alpha, 
                const double&           nu, 
                const double&           beta,
                const vector<double>&   strikes,
                      vector<double>&   res
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
	REAL s0, t, alpha, nu, beta;

    cout<<"\n// Running Original (Parallel) Volatility-Calibration Benchmark"<<endl;

    readDataSet( OUTER_LOOP_COUNT, NUM_X, NUM_Y, NUM_T, s0, t, alpha, nu, beta ); 

	vector<double> strikes(OUTER_LOOP_COUNT),res(OUTER_LOOP_COUNT);

	for(unsigned i=0;i<OUTER_LOOP_COUNT;++i)
		strikes[i] = 0.001*i;

    const int Ps = get_CPU_num_threads();

    unsigned long int elapsed = 0;
    {   // Main Computational Kernel
        struct timeval t_start, t_end, t_diff;
        gettimeofday(&t_start, NULL);

        run_CPUkernel( OUTER_LOOP_COUNT, NUM_X, NUM_Y, NUM_T, s0, t, alpha, nu, beta, strikes, res );

        gettimeofday(&t_end, NULL);
        timeval_subtract(&t_diff, &t_end, &t_start);
        elapsed = t_diff.tv_sec*1e6+t_diff.tv_usec;
    }

    {   // validation and writeback of the result
        bool is_valid = validate   ( res.data(), OUTER_LOOP_COUNT );
        writeStatsAndResult( is_valid, res.data(), OUTER_LOOP_COUNT, 
                             NUM_X, NUM_Y, NUM_T, false, Ps, elapsed );        
    }

    return 0;
}

