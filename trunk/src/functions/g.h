/** \file g.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the function g
 */
#ifndef G_H_
#define G_H_

#include "../general/parameters.h"

namespace mg
{
    /**
     * \brief g() defines a function \f$ g: K^2 \rightarrow K \f$
     * 
     * g() calculates the function:
     * \f[
     * g(x,y) = x^4y^3
     * \f]
     * \see function2D
     *  
     * \param[in] x     the x coordinate to evaluate g on
     * \param[in] y     the y coordinate to evaluate g on
     * \return      the value of g(x,y)
     */
    inline Precision g(
        const Precision x,
        const Precision y) throw()
    {
        return x*x*x*x*y*y*y;
    }
}

#endif /*G_H_*/
