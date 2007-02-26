#include "ChannelProblem.h"

namespace mg
{

ChannelProblem::ChannelProblem(
    Stencil& stencil,
    Point origin,
    Index nx, Index ny,
    Precision hx, Precision hy)
    : Problem(stencil, origin, nx, ny, hx, hy)
{
}

ChannelProblem::~ChannelProblem()
{
}

void ChannelProblem::setBoundaryConstraint( const Function& boundaryConstraint )
{
    const Index nx = getNx();
    const Index ny = getNy();
    const Precision hx = getHy();
    const Precision hy = getHx();
    const Point origin = getOrigin();
    for (Index sx=0; sx<=nx; sx++)
    {
        solution_(sx,-1)=boundaryConstraint(origin.x+sx*hx,origin.y-hy);
        solution_(sx,0)=boundaryConstraint(origin.x+sx*hx,origin.y+0.0);       //top border
        solution_(sx,ny)=boundaryConstraint(origin.x+sx*hx,origin.y+ny*hy);     //bottom border
        solution_(sx,ny+1)=boundaryConstraint(origin.x+sx*hx,origin.y+(ny+1)*hy);
    }
}

void ChannelProblem::applyBoundaryConstraint()
{
    applyBoundaryConstraint(solution_);
}

void ChannelProblem::applyBoundaryConstraint( DiscreteFunction& u ) const
{
    const Index nx = getNx();
    const Index ny = getNy();
    for(Index sy=0;sy<=ny;++sy)
    {
        u(0,sy)=u(nx,sy);
        u(nx+1,sy)=u(1,sy);
    }
}

DiscreteFunction ChannelProblem::residuum()
{
    const Index nx = getNx();
    const Index ny = getNy();
    const Precision hx=getHx();
    const Precision hy=getHy();
    const Point origin = getOrigin();
    DiscreteFunction result(0.0,origin,nx,ny,hx,hy);
    if (stencil_.size()<2)
    {
        for (Index sy=1; sy<ny; sy++)
            for(Index sx=1; sx<=nx; sx++)
                result(sx,sy)=rightHandSide_(sx,sy)
                            -stencil_.apply(solution_,C,sx,sy);
    }
    else
    {
        //TODO: This needs to be fixed
        /*//south west corner
        result(1,1)=rightHandSide_(1,1)
                    -stencil_.apply(solution_,SW,1,1,nx_,ny_);
        //south east corner
        result(nx_-1,1)=rightHandSide_(nx_-1,1)
                    -stencil_.apply(solution_,SE,nx_-1,1,nx_,ny_);
        //north west corner
        result(1,ny_-1)=rightHandSide_(1,ny_-1)
                    -stencil_.apply(solution_,NW,1,ny_-1,nx_,ny_);
        //north east corner
        result(nx_-1,ny_-1)=rightHandSide_(nx_-1,ny_-1)
                    -stencil_.apply(solution_,NE,nx_-1,ny_-1,nx_,ny_);
        //south boarder
        for (Index sx=2; sx<nx_-1; sx++)
            result(sx,1)=rightHandSide_(sx,1)
                    -stencil_.apply(solution_,S,sx,1,nx_,ny_);
        //north boarder
        for (Index sx=2; sx<nx_-1; sx++)
            result(sx,ny_-1)=rightHandSide_(sx,ny_-1)
                    -stencil_.apply(solution_,N,sx,ny_-1,nx_,ny_);
        //west boarder
        for (Index sy=2; sy<ny_-1; sy++)
            result(1,sy)=rightHandSide_(1,sy)
                    -stencil_.apply(solution_,W,1,sy,nx_,ny_);
        //east boarder;
        for (Index sy=2; sy<ny_-1; sy++)
            result(nx_-1,sy)=rightHandSide_(nx_-1,sy)
                    -stencil_.apply(solution_,E,nx_-1,sy,nx_,ny_);
        //the center
        for (Index sy=2; sy<ny_-1; sy++)
            for (Index sx=2; sx<nx_-1; sx++)
                result(sx,sy)=rightHandSide_(sx,sy)
                        -stencil_.apply(solution_,C,sx,sy,nx_,ny_);*/
    }
	applyBoundaryConstraint(result);
    return result;
}

IndexPair ChannelProblem::getFirstPoint(Index, Index) const
{
    return IndexPair(1,1);
}

IndexPair ChannelProblem::getLastPoint(Index nx, Index ny) const
{
    return IndexPair(nx,ny-1);
}

ProblemPtr ChannelProblem::getCoarsGridProblem(
    Index nxNew,
    Index nyNew,
    Precision hxNew,
    Precision hyNew) const
{
    return ProblemPtr(
        new ChannelProblem(stencil_,getOrigin(),nxNew,nyNew,hxNew,hyNew));
}

void ChannelProblem::fillBorderValues(
        NumericArray& matrix,
        NumericArray& rightSide,
        const Index dimension) const
{
    const Index nx = getNx();
    const Index ny = getNy();
    for (Index sx=0; sx<=nx; ++sx)
    {
        //lower border
        matrix[sx*dimension+sx]=1;
        rightSide[sx]=solution_(sx,0);
        //upper border
        matrix[(dimension-1-sx)*dimension+(dimension-1-sx)]=1;
        rightSide[dimension-1-sx]=solution_(nx-sx,ny);
    }
    for (Index sy=1; sy<ny; ++sy)
    {
        matrix[(sy*(nx+1)+0)*dimension+(sy*(nx+1)+0)]=1.0;
        matrix[(sy*(nx+1)+0)*dimension+(sy*(nx+1)+nx)]=-1.0;
        rightSide[(sy*(nx+1)+0)]=0.0;
        
        matrix[(sy*(nx+1)+nx)*dimension+(sy*(nx+1)+nx)]=2.0;
        matrix[(sy*(nx+1)+nx)*dimension+(sy*(nx+1)+nx-1)]=-1.0;
        matrix[(sy*(nx+1)+nx)*dimension+(sy*(nx+1)+1)]+=-1.0;
        rightSide[(sy*(nx+1)+nx)]=0.0;
    }
}

}
