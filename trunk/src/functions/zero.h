/** \file zero.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the function zero
 */
#ifndef ZERO_H_
#define ZERO_H_

#include "../general/parameters.h"

namespace mg
{
    /**
     * \brief zero() defines the function \f$ zero: K^2 \rightarrow \{0\} \f$
     * 
     * zero always returns 0.0 and is needed for evaluation purpose.
     * \see function2D
     * 
     * \param[in]   the x coordinate
     * \param[in]   the y coordinate
     * \return      0.0
     */
    inline Precision zero(
        const Precision,
        const Precision) throw()
    {
        return 0.0;
    }
}

#endif /*ZERO_H_*/
