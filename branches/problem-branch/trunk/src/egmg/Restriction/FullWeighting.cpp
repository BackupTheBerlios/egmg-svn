/** \file FullWeighting.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the class FullWeighting.
 */

#include "FullWeighting.h"
#include <stdexcept>

namespace mg
{

DiscreteFunction FullWeighting::restriction(
        const Problem& problem, const DiscreteFunction& u) const
{
    const Index nx = problem.getNx();
    const Index ny = problem.getNy();
    //if it is not possible to do standart coarsening throw an exeption
    if ((nx%2 != 0) || (ny%2 != 0))
        throw std::domain_error("u");
    const Index nxNew=nx/2;
    const Index nyNew=ny/2;
    DiscreteFunction result(0.0,nxNew,nyNew);
    //do injection on the boarders
    for (Index sy=0;sy<=nyNew;sy++)
    {
            result(0,sy)=u(0,2*sy);
            result(nxNew,sy)=u(nx,2*sy);
    }
    for (Index sx=0;sx<=nxNew;sx++)
    {
            result(sx,0)=u(2*sx,0);
            result(sx,nyNew)=u(2*sx,ny);
    }
    
    for (Index sy=1;sy<nyNew;sy++)
    {
        for(Index sx=1;sx<nxNew;sx++)
            result(sx,sy)=
                weight_*(4*u(2*sx,2*sy)
                        +2*u(2*sx-1,2*sy)
                        +2*u(2*sx+1,2*sy)
                        +2*u(2*sx,2*sy-1)
                        +2*u(2*sx,2*sy+1)
                          +u(2*sx-1,2*sy+1)
                          +u(2*sx-1,2*sy+1)
                          +u(2*sx+1,2*sy-1)
                          +u(2*sx+1,2*sy-1)
                        )/16.0;
    }
    return result;
}

}

