#include "ChannelProblem.h"

namespace mg
{

ChannelProblem::ChannelProblem(
    Stencil& stencil,
    Index nx,
    Index ny)
    : Problem(stencil, nx, ny)
{
}

ChannelProblem::~ChannelProblem()
{
}

void ChannelProblem::setBoundaryConstraint( const Function& boundaryConstraint )
{
    const Precision hx = 1.0/nx_;
    const Precision hy = 1.0/ny_;
    for (Index sx=0; sx<=nx_; sx++)
    {
        solution_(sx,-1)=boundaryConstraint(sx*hx,-hy);
        solution_(sx,0)=boundaryConstraint(sx*hx,0.0);       //top border
        solution_(sx,ny_)=boundaryConstraint(sx*hx,1.0);     //bottom border
        solution_(sx,ny_+1)=boundaryConstraint(sx*hx,1.0+hy);
    }
}

void ChannelProblem::applyBoundaryConstraint()
{
    applyBoundaryConstraint(solution_);
}

void ChannelProblem::applyBoundaryConstraint( DiscreteFunction& u ) const
{
    for(Index sy=0;sy<=ny_;++sy)
    {
        u(0,sy)=u(nx_,sy);
        u(nx_+1,sy)=u(1,sy);
    }
}

DiscreteFunction ChannelProblem::residuum()
{
    DiscreteFunction result(0.0,nx_,ny_);
    if (stencil_.size()<2)
    {
        for (Index sy=1; sy<ny_; sy++)
            for(Index sx=1; sx<=nx_; sx++)
                result(sx,sy)=rightHandSide_(sx,sy)
                            -stencil_.apply(solution_,C,sx,sy,nx_,ny_);
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
    return result;
}

Point ChannelProblem::getFirstPoint() const
{
    return Point(1,1);
}

Point ChannelProblem::getLastPoint() const
{
    return Point(nx_,ny_-1);
}

ChannelProblem* ChannelProblem::getCoarsGridProblem(
    Index nxNew,
    Index nyNew) const
{
    return new ChannelProblem(stencil_,nxNew,nyNew);
}

void ChannelProblem::fillBorderValues(
        NumericArray& matrix,
        NumericArray& rightSide,
        const Index dimension) const
{
    for (Index sx=0; sx<=nx_; ++sx)
    {
        //lower border
        matrix[sx*dimension+sx]=1;
        rightSide[sx]=solution_[sx];
        //upper border
        matrix[(dimension-1-sx)*dimension+(dimension-1-sx)]=1;
        rightSide[dimension-1-sx]=solution_[dimension-1-sx];
    }
    for (Index sy=1; sy<ny_; ++sy)
    {
        matrix[(sy*(nx_+1)+0)*dimension+(sy*(nx_+1)+0)]=1.0;
        matrix[(sy*(nx_+1)+0)*dimension+(sy*(nx_+1)+nx_)]=-1.0;
        rightSide[(sy*(nx_+1)+0)]=0.0;
        
        matrix[(sy*(nx_+1)+nx_)*dimension+(sy*(nx_+1)+nx_)]=2.0;
        matrix[(sy*(nx_+1)+nx_)*dimension+(sy*(nx_+1)+nx_-1)]=-1.0;
        matrix[(sy*(nx_+1)+nx_)*dimension+(sy*(nx_+1)+1)]+=-1.0;
        rightSide[(sy*(nx_+1)+nx_)]=0.0;
    }
}

}
