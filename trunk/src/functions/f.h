/** \file f.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the function f
 */
#ifndef F_H_
#define F_H_
#include "../general/parameters.h"

namespace mg
{
	/**
	 * \brief f() defines a function \f$ f: K^2 \rightarrow K \f$
	 * 
	 * f() calculates the function:
	 * \f[
	 * f(x,y) = -6(x^2y(x^2+2y^2))
	 * \f]
	 * \see function2D
	 * 
	 * \param[in] x		the x coordinate to evaluate f on
	 * \param[in] y		the y coordinate to evaluate f on
	 * \return 		the value of f(x,y)
	 */
	inline precision f(const precision x, const precision y) throw()
	{
		return -6*(x*x*y*(x*x+2*y*y));
	}
}

#endif /*F_H_*/
