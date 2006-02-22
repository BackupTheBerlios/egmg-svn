/** \file error.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface of the function error
 */
#ifndef ERROR_H_
#define ERROR_H_

#include <valarray>
#include "../general/parameters.h"

namespace mg
{
    /**
     * \brief error() calculates the error of the calcuated solution
     * 
     * error() calculates maximum norm of the diffrence between the given
     * discrete solution of a pde and the also given exact solution. This
     * function is of cause only used for evaluation purposes.
     * 
     * \param[in] u         the discrete solution of the discrete pde
     * \param[in] nx        number of steps in x direction
     * \param[in] ny        number of steps in y direction
     * \param[in] solution  the exact solution of the contious pde
     * \return              the maximum norm of the error
     */
    Precision error(
        const std::valarray<Precision>& u,
        const size_t nx,
        const size_t ny,
        const function2D solution);
}

#endif /*ERROR_H_*/
