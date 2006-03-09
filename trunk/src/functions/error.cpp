/** \file error.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the function error
 */

#include <cmath>
#include <algorithm>
#include "error.h"

namespace mg
{
    Precision error(
        const NumericArray& u,
        const Index nx,
        const Index ny,
        const function2D solution)
    {
        Precision result=0.0;
        const Precision hx=1.0/nx;
        const Precision hy=1.0/ny;
        for (Index sx=1; sx<nx; sx++)
            for(Index sy=1; sy<ny; sy++)
            {
                Precision temp=std::fabs(solution(sy*hx,sx*hy)-u[sx*(nx+1)+sy]);
                result=std::max(result,temp);
            }
            
        return result;
    }
}
