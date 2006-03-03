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
        const NumericArray& u,
        const size_t nx,
        const size_t ny,
        const function2D solution)
    {
        Precision result=0.0;
        const Precision hx=1.0/nx;
        const Precision hy=1.0/ny;
        for (size_t sx=1; sx<nx; sx++)
            for(size_t sy=1; sy<ny; sy++)
            {
                Precision temp=std::fabs(solution(sy*hx,sx*hy)-u[sx*(nx+1)+sy]);
                result=std::max(result,temp);
            }
            
        return result;
    }
}
