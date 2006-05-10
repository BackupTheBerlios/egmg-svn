/** \file SPAI.cpp
 * \author Andre Oeckerath
 * \brief SPAI.cpp contains the implementation of the class SPAI.
 * \see SPAI.h
 */
 
#include "SPAI.h"
#include "../functions/residuum.h"
#include <iostream>
#include <cmath>


namespace mg
{
	
void choleskySolve(NumericArray& W, const NumericArray& phiTable,
		const Stencil& stencil, const PositionArray& phiX, const PositionArray& phiY,
		const int weightingSize, const Index px, const Index py, const Index nx,
		const Index ny) 
{
	Precision temp = 0.0;
	Precision innertemp = 0.0;
	Precision outertemp = 0.0;
	Precision sum = 0.0;

	Precision Matrix[9][9];
	NumericArray rhs(0.0,weightingSize);
	NumericArray V(0.0,weightingSize);

	NumericArray operatorL=stencil.getL(C,px,py,nx,ny);
    PositionArray jX=stencil.getJx(C);
    PositionArray jY=stencil.getJy(C);

	//erzeuge rechte Seite
	for(int i=0; i<weightingSize; i++)
	{
		temp = 0.0;
		for(Index nue=0; nue<operatorL.size(); nue++)
		{
			temp += operatorL[nue] 
				* phiTable[6-jX[nue]-phiX[i]+13*(6-jY[nue]-phiY[i])];
		}		
		rhs[i] = temp;
	}

	//erzeuge Matrix
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

	//  Matrix Check
/*	std::cout << "Ausgangsmatrix:" << std::endl;
	for(int i=0; i<9;i++)
	{
		for(int j=6; j<9; j++)
		{
			std::cout<< Matrix[i][j] << "\t  ";
		}
		std::cout<< std::endl << std::endl;
	}
*/
	// Cholesky Zerlegung für Matrix * W = rhs
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

/*	//  Cholesky Check
	std::cout << "Choleskyzerlegung:" << std::endl;
	for(int i=0; i<9;i++)
	{
		for(int j=i; j<9; j++)
		{
			std::cout<< Matrix[i][j] << "\t  ";
		}
		std::cout<< std::endl << std::endl;
	}
*/

/*	// rhs Check
	std::cout << "rhs:" << std::endl;
	for(int i=0; i<weightingSize;i++)
	{		
		std::cout<< rhs[i] << std::endl;
	}
	std::cout << std::endl;
*/
	// Solve G^T V = rhs
	for(int i=0; i<weightingSize; i++)
	{
		sum= 0.0;
		for(int j=0; j<i;j++)
		{
			sum += Matrix[j][i]*V[j];  //eigentlich Matrix[i][j] aber betrachte transponierte Matrix
		}
		V[i] = (rhs[i] - sum)/Matrix[i][i];
	}

/*	// V Check
	std::cout << "V:" << std::endl;
	for(int i=0; i<weightingSize;i++)
	{
		
		std::cout<< V[i] << std::endl;
	}
	std::cout << std::endl;
*/
	// Solve G W = V
	for(int i=0; i<weightingSize; i++)
	{
		sum=0.0;
		for(int j=weightingSize-i; j<weightingSize ; j++)
		{
			sum += Matrix[weightingSize-1-i][j] * W[j];
		}
		W[weightingSize-1-i] = (V[weightingSize-1-i] - sum)/Matrix[weightingSize-1-i][weightingSize-1-i];
	}

/*	// W Check
	std::cout << "Gewichtungsoperator W:" << std::endl;
	for(int i=0; i<weightingSize;i++)
	{
		
		std::cout<< W[i] << std::endl;
	}
	std::cout << std::endl;
*/			
}

///////////////////////////////////////////////////////////////////////////

void SPAI::relax(
    NumericArray &u,
    const NumericArray &f, 
    const Stencil &stencil,
    const Index nx,
    const Index ny) const
{	
    NumericArray W(0.0,9);
	NumericArray resid = residuum(u,f,stencil,nx,ny);
	Precision sum=0.0;

//  Table Check
/*	for(int i=3; i<10;i++)
	{
		for(int j=3; j<10; j++)
		{
			std::cout<< phiTable[i+j*13] << "   \t   ";
		}
		std::cout<< std::endl << std::endl;
	}
*/

//	choleskySolve(W, phiTable, stencil, phiX, phiY, weightingSize_, 1, 1, nx, ny);

/*	for(Index px=1; px<nx; px++)
	{
		for(Index py=1; py<ny; py++)
		{
		   choleskySolve(W, phiTable, stencil, phiX, phiY, px, py, nx, ny);
		}
	}
*/
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
				u[i+j*(ny+1)] += sum;
			}
		}
	}
	else
	{
		for(Index i=1; i<nx; i++)
		{
			for(Index j=1; j<ny; j++)
			{
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
}

