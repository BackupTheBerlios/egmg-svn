/** \file FullWeighting.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the class FullWeighting.
 */

#include "FullWeighting.h"
#include <stdexcept>

namespace mg
{
NumericArray FullWeighting::restriction(
    const NumericArray& u,
    const Stencil&,
    const Index nx,const Index ny) const
{
    //if it is not possible to do standart coarsening throw an exeption
    if ((nx%2 != 0) || (ny%2 != 0))
        throw std::domain_error("u");
    const Index nxNew=nx/2;
    const Index nyNew=ny/2;
    NumericArray result(0.0,(nxNew+1)*(nyNew+1));
    //do injection on the boarders
    for (Index sy=0;sy<=nyNew;sy++)
    {
        result[sy*(nxNew+1)]=u[2*sy*(nx+1)];
        result[sy*(nxNew+1)+nxNew]=u[2*sy*(nx+1)+nx];
    }
    for (Index sx=0;sx<=nxNew;sx++)
    {
        result[sx]=u[2*sx];
        result[nyNew*(nxNew+1)+sx]=u[ny*(nx+1)+2*sx];
    }
    for (Index sy=1;sy<nyNew;sy++)
    {
        for(Index sx=1;sx<nxNew;sx++)
            result[sy*(nxNew+1)+sx]=
                weight_*(4*u[2*sy*(nx+1)+2*sx]
                        +2*u[2*sy*(nx+1)+2*sx+1]+2*u[2*sy*(nx+1)+2*sx-1]
                        +2*u[(2*sy+1)*(nx+1)+2*sx]+2*u[(2*sy-1)*(nx+1)+2*sx]
                          +u[(2*sy+1)*(nx+1)+2*sx+1]+u[(2*sy+1)*(nx+1)+2*sx-1]
                          +u[(2*sy-1)*(nx+1)+2*sx+1]+u[(2*sy-1)*(nx+1)+2*sx-1]
                        )/16.0;
    }
    return result;
}
}

