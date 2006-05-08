/** \file SPAI.cpp
 * \author Andre Oeckerath
 * \brief SPAI.cpp contains the implementation of the class SPAI.
 * \see SPAI.h
 */
 
#include "SPAI.h"
#include "../functions/residuum.h"
#include <iostream>


namespace mg
{
	
void choleskySolve(NumericArray& W, const NumericArray& phiTable,
		const Stencil& stencil, const PositionArray& phiX, const PositionArray& phiY,
		const Index px, const Index py, const Index nx, const Index ny) 
{
	Precision temp = 0.0;
	Precision innertemp = 0.0;
	Precision outertemp = 0.0;

	Precision Matrix[9][9];
	NumericArray rhs(0.0,9);



	NumericArray operatorL=stencil.getL(C,px,py,nx,ny);
    PositionArray jX=stencil.getJx(C);
    PositionArray jY=stencil.getJy(C);

	//erzeuge rechte Seite
	for(int i=0; i<9; i++)
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
	for(int i=0; i<9; i++)
	{
		for(int j=0; j<9; j++)
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
/*	for(int i=4; i<9;i++)
	{
		for(int j=4; j<9; j++)
		{
			std::cout<< Matrix[i][j] << "\t  ";
		}
		std::cout<< std::endl << std::endl;
	}
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

//  Table Check
	for(int i=3; i<10;i++)
	{
		for(int j=3; j<10; j++)
		{
			std::cout<< phiTable[i+j*13] << "   \t   ";
		}
		std::cout<< std::endl << std::endl;
	}


/*	for(Index px=1; px<nx; px++)
	{
		for(Index py=1; py<ny; py++)
		{
		   choleskySolve(W, phiTable, stencil, phiX, phiY, px, py, nx, ny);
		}
	}
*/

	
}
}

