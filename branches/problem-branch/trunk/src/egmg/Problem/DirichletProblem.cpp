#include "DirichletProblem.h"

namespace mg
{

DirichletProblem::DirichletProblem(
    Stencil& stencil,
    Index nx,
    Index ny)
    : Problem(stencil, nx, ny)
{
}

DirichletProblem::~DirichletProblem()
{
}

void DirichletProblem::setBoundaryConstraint( const Function& boundaryConstraint )
{
    const Precision hx = 1.0/nx_;
    const Precision hy = 1.0/ny_;
    solution_(-1,-1)=boundaryConstraint(-hx,-hy);
    for (Index sx=0; sx<=nx_; sx++)
    {
        solution_(sx,-1)=boundaryConstraint(sx*hx,-hy);
        solution_(sx,0)=boundaryConstraint(sx*hx,0.0);       //top border
        solution_(sx,ny_)=boundaryConstraint(sx*hx,1.0);     //bottom border
        solution_(sx,ny_+1)=boundaryConstraint(sx*hx,1.0+hy);
    }
    solution_(nx_+1,-1)=boundaryConstraint(1.0+hx,-hy);
    solution_(-1,ny_+1)=boundaryConstraint(-hx,1.0+hy);
    for (Index sy=1; sy<ny_; sy++)
    {
        solution_(-1,sy)=boundaryConstraint(-hx,sy*hy);
        solution_(0,sy)=boundaryConstraint(0.0,sy*hy);        //left border
        solution_(nx_,sy)=boundaryConstraint(1.0,sy*hy);      //right border
        solution_(nx_+1,sy)=boundaryConstraint(1.0+hx,sy*hy);
    }
    solution_(nx_+1,ny_+1)=boundaryConstraint(1.0+hx,1.0+hy);
}

void DirichletProblem::applyBoundaryConstraint()
{
}

void DirichletProblem::applyBoundaryConstraint( DiscreteFunction& ) const
{
}

DiscreteFunction DirichletProblem::residuum()
{
    DiscreteFunction result(0.0,nx_,ny_);
    if (stencil_.size()<2)
    {
        for (Index sy=1; sy<ny_; sy++)
        {
            for(Index sx=1; sx<nx_; sx++)
            {
                result(sx,sy)=rightHandSide_(sx,sy)
                            -stencil_.apply(solution_,C,sx,sy,nx_,ny_);
            }
        }
    }
    else
    {
        //south west corner
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
                        -stencil_.apply(solution_,C,sx,sy,nx_,ny_);
    }
    return result;
}

Point DirichletProblem::getFirstPoint() const
{
    return Point(1,1);
}

Point DirichletProblem::getLastPoint() const
{
    return Point(nx_-1,ny_-1);
}

DirichletProblem* DirichletProblem::getCoarsGridProblem(
    Index nxNew,
    Index nyNew) const
{
    return new DirichletProblem(stencil_,nxNew,nyNew);
}

void DirichletProblem::fillBorderValues(
        NumericArray& matrix,
        NumericArray& rightSide,
        const Index dimension) const
{
    //border values XDIR
    for (Index sx=0; sx<=nx_; ++sx)
    {
        //lower border
        matrix[sx*dimension+sx]=1;
        rightSide[sx]=solution_[sx];
        //upper border
        matrix[(dimension-1-sx)*dimension+(dimension-1-sx)]=1;
        rightSide[dimension-1-sx]=solution_[dimension-1-sx];
    }
    //border values YDIR
    //corners have been process in XDIR (y=1..ny-1 instead of y=0..ny)
    for (Index sy=1; sy<ny_; ++sy)
    {
        //left border
        matrix[sy*(nx_+1)*dimension+sy*(nx_+1)]=1;
        rightSide[sy*(nx_+1)]=solution_[sy*(nx_+1)];
        //right border
        matrix[(sy*(nx_+1)+nx_)*dimension+(sy*(nx_+1)+nx_)]=1;
        rightSide[sy*(nx_+1)+nx_]=solution_[sy*(nx_+1)+nx_];
    }
}

}
