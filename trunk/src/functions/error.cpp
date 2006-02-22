/** \file error.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the function error
 */

#include <cmath>
#include <utility>
#include "error.h"

namespace mg
{
    Precision error(
        const std::valarray<Precision>& u,
        const size_t nx,
        const size_t ny,
        const function2D solution)
    {
        Precision result=0.0;
        const Precision hx=1.0/nx;
        const Precision hy=1.0/ny;
        for (size_t j=1; j<nx; j++)
            for(size_t i=1; i<ny; i++)
            {
                Precision temp=std::fabs(solution(i*hx,j*hy)-u[j*(nx+1)+i]);
                result=std::max(result,temp);
            }
            
        return result;
    }
}
