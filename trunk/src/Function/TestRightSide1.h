#ifndef TESTRIGHTSIDE1_H_
#define TESTRIGHTSIDE1_H_

#include "Function.h"

namespace mg
{

/**
 * \brief TestRightSide1 defines a function \f$ f: K^2 \rightarrow K \f$
 * 
 * f() calculates the function:
 * \f[
 * f(x,y) = -6(x^2y(x^2+2y^2))
 * \f]
 * 
 * \return      the value of f(x,y)
 */

 // * \param[in] x     the x coordinate to evaluate f on
 // * \param[in] y     the y coordinate to evaluate f on
class TestRightSide1: public Function
{
protected:
	virtual Precision apply(Precision x, Precision y) const;
};

}

#endif
