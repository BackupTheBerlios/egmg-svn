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
    const Precision hx = problem.getHx();
    const Precision hy = problem.getHy();
    const Point origin = problem.getOrigin();
    //if it is not possible to do standart coarsening throw an exeption
    if ((nx%2 != 0) || (ny%2 != 0))
        throw std::domain_error("u");
    const Index nxNew=nx/2;
    const Index nyNew=ny/2;
    const Precision hxNew = 2.0*hx;
    const Precision hyNew = 2.0*hy;
    DiscreteFunction result(0.0,origin,nxNew,nyNew,hxNew,hyNew);
    
    const IndexPair fp = problem.getFirstPoint();
    const IndexPair lp = problem.getLastPoint();
    
    if ( fp.sx > 0 )
    {
        //do injection on the south boarder
        for (Index sy=0;sy<=nyNew;sy++)
            result(0,sy)=u(0,2*sy);
    }
    if ( lp.sx < nx )
    {
        //do injection on the north boarder
        for (Index sy=0;sy<=nyNew;sy++)
            result(nxNew,sy)=u(nx,2*sy);
    }
    if ( fp.sy > 0 )
    {
        // do injection on the west boarder
        for (Index sx=0;sx<=nxNew;sx++)
            result(sx,0)=u(2*sx,0);
    }
    if ( lp.sy < ny )
    {
        //do injection on the east boarder
        for (Index sx=0;sx<=nxNew;sx++)
            result(sx,nyNew)=u(2*sx,ny);
    }
    
    IndexPair nfp = problem.getFirstPoint(nxNew,nyNew);
    IndexPair nlp = problem.getLastPoint(nxNew,nyNew);

    for (Index sy=nfp.sy;sy<=nlp.sy;sy++)
    {
        for(Index sx=nfp.sx;sx<=nlp.sx;sx++)
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

const NumericArray& FullWeighting::getI(
    const Position,
    const Index,
    const Index,
    const Index,
    const Index,
    const Precision,
    const Precision,
    const Point,
    const Stencil&) const
{
    return i_;  
}

const PositionArray& FullWeighting::getJx( const Position ) const
{
    return jx_; 
}

const PositionArray& FullWeighting::getJy( const Position ) const
{
    return jy_;
}

}

