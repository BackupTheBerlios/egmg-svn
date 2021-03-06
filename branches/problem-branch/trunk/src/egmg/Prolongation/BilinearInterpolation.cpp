/** \file BilinearInterpolation.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the class BilinearInterpolation
 */
#include "BilinearInterpolation.h"

namespace mg
{
    
DiscreteFunction BilinearInterpolation::prolongate(
    const Problem& problem) const
{
    const Index nx = problem.getNx();
    const Index ny = problem.getNy();
    const Index nxNew = 2*nx;
    const Index nyNew = 2*ny;
    const Precision hx = problem.getHx();
    const Precision hy = problem.getHy();
    const Precision hxNew = hx/2.0;
    const Precision hyNew = hy/2.0;
    const Point origin = problem.getOrigin();
    const DiscreteFunction& u = problem.getSolution();
    DiscreteFunction result(0.0,origin,nxNew,nyNew,hxNew,hyNew);

    //"interpolation" of coarse grid points
    for (Index sy=0; sy<=ny; ++sy)
        for (Index sx=0; sx<=nx; ++sx)
            result(2*sx,2*sy)=u(sx,sy);

    //interpolation of fine grid points on coarse grid lines
    for (Index sy=0; sy<=nyNew; sy+=2)
        for (Index sx=1; sx<=nxNew; sx+=2)
            result(sx,sy)=
                1./2*(result(sx-1,sy)
                     +result(sx+1,sy));
    //interpolation of fine grid points on fine grid lines
    for (Index sy=1; sy<=nyNew; sy+=2)
        for (Index sx=0; sx<=nxNew; ++sx)
            result(sx,sy)=
                1./2*(result(sx,sy-1)
                     +result(sx,sy+1));
    return result;
}

const NumericArray& BilinearInterpolation::getI(
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
    
const PositionArray& BilinearInterpolation::getJx( const Position ) const
{
    return jx_; 
}

const PositionArray& BilinearInterpolation::getJy( const Position ) const
{
    return jy_;
}

}
