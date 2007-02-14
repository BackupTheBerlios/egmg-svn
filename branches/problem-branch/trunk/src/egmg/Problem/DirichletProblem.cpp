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
    for (Index sx=0; sx<=nx_; sx++)
    {
        solution_(sx,0)=boundaryConstraint(sx*hx,0.0);       //top border
        solution_(sx,ny_)=boundaryConstraint(sx*hx,1.0);     //bottom border
    }
    for (Index sy=1; sy<ny_; sy++)
    {
        solution_(0,sy)=boundaryConstraint(0.0,sy*hy);        //left border
        solution_(nx_,sy)=boundaryConstraint(1.0,sy*hy);      //right border
    }
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

}
