/** \file WHighJac.h
 * \author Andre Oeckerath
 * \brief WHighJac.h contains the interface of the class WHighJac.
 * \see Relaxation.h
 */
#ifndef WHighJac_H_
#define WHighJac_H_


#include "Relaxation.h"
#include "../functions/WHighJacScalarProduct.h"

namespace mg
{
/**
 * \brief WHighJac is an abstract class for weighted relaxation using
 *  the WHighJac approach for minimising the L_{high}^2 Norm
 */
class WHighJac : public mg::Relaxation
{
private: 
	const int weightingSize_;
	NumericArray phiTable;
	
	void initTable() 
	{
		for (int i=-6; i<7; i++)
		{
			for (int j=-6; j<7; j++)
			{
				phiTable[6+i+13*(6+j)] = 0.0;	
			}
		}
		// only tableelements with uneven index (exept (0/0))
		// have a value != 0		
		for (int i=-5; i<6; i+=2)
		{
			for (int j=-5; j<6; j+=2)
			{
				phiTable[6+i+13*(6+j)] = WHighJacScalarProduct(i,j);
			}
		}	
        for (int i=-5; i<6; i+=2)
		{
			phiTable[6+i+13*6] = WHighJacScalarProduct(i,0);
		}
        for (int j=-5; j<6; j+=2)
		{
			phiTable[6+13*(6+j)] = WHighJacScalarProduct(0,j);
		}
		phiTable[6+13*6] = WHighJacScalarProduct(0,0);

	}	

	PositionArray jXW;
	PositionArray jYW;
	PositionArray phiX;
	PositionArray phiY;

	void initArrays() 
	{
		jXW[0] = phiX[0] = 0; jYW[0] = phiY[0] = 0;
		jXW[1] = phiX[1] = -1; jYW[1] = phiY[1] = 0;
		jXW[2] = phiX[2] = 0; jYW[2] = phiY[2] = 1;
		jXW[3] = phiX[3] = 1; jYW[3] = phiY[3] = 0;
		jXW[4] = phiX[4] = 0; jYW[4] = phiY[4] = -1;
		jXW[5] = phiX[5] = -1; jYW[5] = phiY[5] = 1;
		jXW[6] = phiX[6] = 1; jYW[6] = phiY[6] = 1;
		jXW[7] = phiX[7] = 1; jYW[7] = phiY[7] = -1;
		jXW[8] = phiX[8] = -1; jYW[8] = phiY[8] = -1;
	}

	void choleskySolve(NumericArray& W, const NumericArray& phiTable,
		const Stencil& stencil, const PositionArray& phiX, const PositionArray& phiY,
		const int weightingSize, const Index px, const Index py, const Index nx,
		const Index ny) const;
    //We don't want auto generated copy constructors and assignment operators
	WHighJac( const WHighJac& );
	WHighJac& operator=(const WHighJac&);
public:
    /**
     * \brief The constructor of a WHighJac object
     * 
     * WHighJac constructs a WHighJac object with:
     * \param[in]  weightingSize    number of elements in the weighting operator
     */
    WHighJac(const int weightingSize = 9)
        : weightingSize_(weightingSize),phiTable(169),jXW(9), jYW(9),
		phiX(9),phiY(9)	{ initTable(); initArrays(); }

    virtual ~WHighJac() {}

    /**
     * \brief relax() executes one relaxation step on the input vector
     * 
     * relax() exectues one weighted SPAI relaxation step on the input vector
     * one a rectangular 2D grid with lexicographic ordering and the
     * discretazation given by stencil for a pde.
     * 
     * \param[in,out] u     the vector representation of the 2D grid to perform
     *                      the relaxation on
     * \param[in] f         the right hand side of the pde
     * \param[in] stencil   the stencil rep. of the pde
     * \param[in] nx        number of steps in x direction
     * \param[in] ny        number of steps in y direction
     */
    void relax(
        NumericArray& u,
        const NumericArray& f,
        const Stencil& stencil,
        const Index nx,
        const Index ny) const;
	
	void OutputW(const Stencil& stencil, const Index i, const Index j, const Index nx,
		const Index ny) const; 
};


}
#endif /*WHighJac_H_*/
