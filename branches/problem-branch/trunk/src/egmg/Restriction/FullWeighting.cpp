/** \file FullWeighting.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the class FullWeighting.
 */

#include "FullWeighting.h"
#include <stdexcept>

namespace mg
{
NumericArray FullWeighting::restriction(
        const Problem& problem) const
{
    const Index nx = problem.getNx();
    const Index ny = problem.getNy();
    const Point llc = problem.getLowerLeftCorner();
    const Point fp = problem.getFirstPoint();
    const Point lp = problem.getLastPoint();
    const NumericArray& u = problem.getSolution();
    //if it is not possible to do standart coarsening throw an exeption
    if ((nx%2 != 0) || (ny%2 != 0))
        throw std::domain_error("u");
    const Index nxNew=nx/2;
    const Index nyNew=ny/2;
    NumericArray result(0.0,(nxNew+1)*(nyNew+1));
    //do injection on the boarders
    for (Index sy=fp.sy;sy<=lp.sy/2;sy++)
    {
        result[sy*(nxNew+1)]=u[2*sy*(nx+1)];
        result[sy*(nxNew+1)+lp.sy/2]=u[2*sy*(nx+1)+lp.sy];
    }
    for (Index sx=fp.sx;sx<=lp.sx;sx++)
    {
        result[sx]=u[2*sx];
        result[nyNew*(nxNew+1)+sx]=u[ny*(nx+1)+2*sx];
    }
    for (Index sy=fp.sy;sy<lp.sy/2;sy++)
    {
        for(Index sx=fp.sx;sx<lp.sx/2;sx++)
            result[sy*(nxNew+1)+sx]=
                weight_*(4*u[2*sy*(nx+1)+2*sx]
                        +2*u[2*sy*(nx+1)+2*sx+1]+2*u[2*sy*(nx+1)+2*sx-1]
                        +2*u[(2*sy+1)*(nx+1)+2*sx]+2*u[(2*sy-1)*(nx+1)+2*sx]
                          +u[(2*sy+1)*(nx+1)+2*sx+1]+u[(2*sy+1)*(nx+1)+2*sx-1]
                          +u[(2*sy-1)*(nx+1)+2*sx+1]+u[(2*sy-1)*(nx+1)+2*sx-1]
                        )/16.0;
    }
    problem.applyBoundaryConstraint( result );
    return result;
}
}

