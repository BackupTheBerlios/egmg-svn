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
        const Problem& problem,
        const Function& solution)
    {
        Precision result=0.0;
        const Index nx = problem.getNx();
        const Index ny = problem.getNy();
        const DiscreteFunction& u = problem.getSolution();
        const Precision hx=1.0/nx;
        const Precision hy=1.0/ny;
        for (Index sx=0; sx<=nx; sx++)
            for(Index sy=0; sy<=ny; sy++)
            {
                Precision temp=std::fabs(solution(sx*hx,sy*hy)-u(sx,sy));
                result=std::max(result,temp);
            }
            
        return result;
    }
}
