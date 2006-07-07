/** \file WHighJac.cpp
 * \author Andre Oeckerath
 * \brief WHighJac.cpp contains the implementation of the class WHighJac.
 * \see WHighJac.h
 */
 
#include "WHighJac.h"
#include "../functions/residuum.h"
#include <iostream>
#include <cmath>


namespace mg
{

// calculation of weighting operator W in gridpoint (px/py)	
void WHighJac::choleskySolve(NumericArray& W, const NumericArray& phiTable,
		const Stencil& stencil, const PositionArray& phiX, const PositionArray& phiY,
		const int weightingSize, const Index px, const Index py, const Index nx,
		const Index ny) const
{
	Precision temp = 0.0;
	Precision innertemp = 0.0;
	Precision outertemp = 0.0;
	Precision sum = 0.0;

	Precision Matrix[9][9];
	NumericArray rhs(0.0,weightingSize);
	NumericArray V(0.0,weightingSize);

	NumericArray operatorL=stencil.getL(C,px,py,nx,ny);
    PositionArray jX=stencil.getJx(C,nx,ny);
    PositionArray jY=stencil.getJy(C,nx,ny);

	// create right side
	for(int i=0; i<weightingSize; i++)
	{
		temp = 0.0;
		for(Index nue=0; nue<operatorL.size(); nue++)
		{
			temp += operatorL[nue] 
				* phiTable[6-jX[nue]-phiX[i]+13*(6-jY[nue]-phiY[i])];
		}		
		rhs[i] = temp*operatorL[0];
	}

	//create Matrix
	for(int i=0; i<weightingSize; i++)
	{
		for(int j=0; j<weightingSize; j++)
		{
			outertemp = 0.0;
			for(Index mue=0; mue<operatorL.size(); mue++)
			{
				innertemp = 0.0;
				for(Index nue=0; nue<operatorL.size(); nue++)
				{
					innertemp += operatorL[nue] 
					* phiTable[6+jX[mue]-jX[nue]+phiX[i]-phiX[j]+
					       13*(6+jY[mue]-jY[nue]+phiY[i]-phiY[j])];
				}
				outertemp += operatorL[mue] * innertemp;
			}
			Matrix[i][j] = outertemp;			
		}
	}

	// Cholesky decomposition for Matrix
	Matrix[0][0] = sqrt(Matrix[0][0]);
	for(int j=1; j<weightingSize; j++)
	{
		Matrix[0][j] = Matrix[0][j] / Matrix[0][0];
	}
	for(int i=1; i<weightingSize; i++)
	{
		sum = 0.0;
		for(int k=0; k < i; k++)
		{
			sum += Matrix[k][i] * Matrix[k][i];
		}
		Matrix[i][i] = sqrt(Matrix[i][i] - sum);

		for(int j=i+1 ; j<weightingSize; j++)
		{  
			sum = 0.0;
			for(int k=0; k < i; k++)
			{
				sum += Matrix[k][i] * Matrix[k][j];
			}
			Matrix[i][j] = (Matrix[i][j] - sum)/Matrix[i][i];
		}
	}

	// Solve G^T * V = rhs
	for(int i=0; i<weightingSize; i++)
	{
		sum= 0.0;
		for(int j=0; j<i;j++)
		{
			sum += Matrix[j][i]*V[j];  
		}
		V[i] = (rhs[i] - sum)/Matrix[i][i];
	}

	// Solve G * W = V
	for(int i=0; i<weightingSize; i++)
	{
		sum=0.0;
		for(int j=weightingSize-i; j<weightingSize ; j++)
		{
			sum += Matrix[weightingSize-1-i][j] * W[j];
		}
		W[weightingSize-1-i] = (V[weightingSize-1-i] - sum)/Matrix[weightingSize-1-i][weightingSize-1-i];
	}
}

///////////////////////////////////////////////////////////////////////////

void WHighJac::relax(
    NumericArray &u,
    const NumericArray &f, 
    const Stencil &stencil,
    const Index nx,
    const Index ny) const
{	
    NumericArray W(0.0,9);
	NumericArray resid = residuum(u,f,stencil,nx,ny);
	Precision sum=0.0;
	NumericArray operatorL=stencil.getL(C,2,2,nx,ny);

	//relax for stencil's with constant coefficients
	if( stencil.isConstant() == true )
	{
		choleskySolve(W, phiTable, stencil, phiX, phiY, weightingSize_, 2, 2, nx, ny);

		for(Index i=1; i<nx; i++)
		{
			for(Index j=1; j<ny; j++)
			{
				sum=0.0;
				for(int k=0; k<weightingSize_; k++)
				{
					sum += W[k] * resid[i+jXW[k]+(j+jYW[k])*(ny+1)];
				}
				u[i+j*(ny+1)] += sum/operatorL[0];
			}
		}
	}
	//relax for stencil's with nonconstant coefficients
	else
	{
		for(Index i=1; i<nx; i++)
		{
			for(Index j=1; j<ny; j++)
			{
				operatorL=stencil.getL(C,i,j,nx,ny);
				sum=0.0;
				choleskySolve(W, phiTable, stencil, phiX, phiY, weightingSize_, i, j, nx, ny);
				for(int k=0; k<weightingSize_; k++)
				{
					sum += W[k] * resid[i+jXW[k]+(j+jYW[k])*(ny+1)];
				}
				u[i+j*(ny+1)] += sum;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////

// command line output of the weighting operator W in gridpoint (i/j)
void WHighJac::OutputW(const Stencil& stencil, const Index i, const Index j, const Index nx,
		const Index ny) const
{	

	NumericArray W(0.0,9);
	choleskySolve(W, phiTable, stencil, phiX, phiY, weightingSize_, i, j, nx, ny);
	
	std::cout << "weighting operator W in point (" << i << "/" << j << ")" << std::endl;
	std::cout<< W[5] << " \t " << W[2] << " \t " << W[6] << std::endl;
	std::cout<< W[1] << " \t " << W[0] << " \t " << W[3] << std::endl;
	std::cout<< W[8] << " \t " << W[4] << " \t " << W[7] << std::endl;
	std::cout << std::endl;
}

}

