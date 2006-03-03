/** \file TransferOperator.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface of the abstract class TransferOperator
 */
#ifndef TRANSFEROPERATOR_H_
#define TRANSFEROPERATOR_H_


#include "parameters.h"
#include "../Stencil/Stencil.h"

namespace mg
{

/**
 * \brief TransferOperator is the interface of a 2D transfer operator
 * \todo Fix/Implement Matrix Dependent Transfer Operators
 */
class TransferOperator
{
public:
	/**
	 * \brief getI returns the stencil representation of the TransferOperator
	 * 
	 * Retruns the coefficents of the stencil representation as valarray. For a
	 * 9 Point Stencil the entries are numbered like this:
	 * \f[
	 * \left[
	 * \begin{array}{ccc}
	 * 5 & 2 & 6 \\
	 * 1 & 0 & 3 \\
	 * 8 & 4 & 7 \\
	 * \end{array}
	 * \right]
	 * =
	 * \left[
	 * \begin{array}{ccc}
	 * nw & n & ne \\
	 * w & c & e \\
	 * sw & s & se \\
	 * \end{array}
	 * \right]
	 * \f]
	 * The enum Postion has the same ordering, please use the enum values for
	 * addresing. It is guaranteed that all 9 Values are set. To iterate over
	 * the entries you can also use the valarrays returnd by getJx and
	 * getJy. This is analog to the getL and getJx getJy Methods of the
	 * Stencil interface.
	 * 
	 * \see Position
	 * \see Stencil.getL
	 * 
	 * \param[in] sx		the x position
	 * \param[in] sy		the y position
	 * \param[in] nx		number of points in x direction
	 * \param[in] ny		number of points in y direction
	 * \param[in] stencil	the stencil rep. of the PDE needed for matrix
	 * 						dependend TransferOperators
	 * \return				the coefficients of the TransferOperator
	 */
	virtual const NumericArray& getI(
		const Index sx,
		const Index sy,
		const Index nx,
		const Index ny, 
		const Stencil& stencil) const =0;
		
	/**
	 * \brief returns the coordinate vector in x dir. for the coefficient vec. I
	 * \see getI
	 * 
	 * \return			the coordinates vector in y direction
	 */
	virtual const PositionArray& getJx() const =0;
	
	/**
	 * \brief returns the coordinate vector in y dir. for the coefficient vec. I
	 * \see getI
	 * 
	 * \return			the coordinates vector in y direction
	 */
	virtual const PositionArray& getJy() const =0;
};

}

#endif /*TRANSFEROPERATOR_H_*/
