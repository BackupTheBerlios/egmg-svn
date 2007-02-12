#include "DirichletProblem.h"

namespace mg
{

DirichletProblem::DirichletProblem(
    Stencil& stencil,
    Index nx,
    Index ny)
    : Problem(stencil, nx, ny), solution_(0.0,(nx+1)*(ny+1))
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
        solution_[sx]=boundaryConstraint(sx*hx,0);                     //top border
        solution_[solution_.size()-1-sx]=boundaryConstraint((nx_-sx)*hx,1);     //bottom border
    }
    for (Index sy=0; sy<=ny_; sy++)
    {
        solution_[sy*(nx_+1)]=boundaryConstraint(0,sy*hy);              //left border
        solution_[(sy+1)*(nx_+1)-1]=boundaryConstraint(1,sy*hy);        //right border
    }
}

void DirichletProblem::applyBoundaryConstraint()
{
}

void DirichletProblem::applyBoundaryConstraint( NumericArray& ) const
{
}

NumericArray& DirichletProblem::getSolution()
{
    return solution_;
}

const NumericArray& DirichletProblem::getSolution() const
{
    return solution_;
}

NumericArray DirichletProblem::residuum()
{
    NumericArray result(0.0,solution_.size());
    if (stencil_.size()<2)
    {
        for (Index sy=1; sy<ny_; sy++)
            for(Index sx=1; sx<nx_; sx++)
                result[sy*(nx_+1)+sx]=
                        rightHandSide_[sy*(nx_+1)+sx]
                            -stencil_.apply(solution_,C,sx,sy,nx_,ny_);
    }
    else
    {
        //south west corner
        result[1*(nx_+1)+1]=
                rightHandSide_[1*(nx_+1)+1]
                    -stencil_.apply(solution_,SW,1,1,nx_,ny_);
        //south east corner
        result[1*(nx_+1)+(nx_-1)]=
                rightHandSide_[1*(nx_+1)+(nx_-1)]
                    -stencil_.apply(solution_,SE,nx_-1,1,nx_,ny_);
        //north west corner
        result[(nx_-1)*(nx_+1)+1]=
                rightHandSide_[(nx_-1)*(nx_+1)+1]
                    -stencil_.apply(solution_,NW,1,ny_-1,nx_,ny_);
        //north east corner
        result[(nx_-1)*(nx_+1)+(nx_-1)]=
                rightHandSide_[(nx_-1)*(nx_+1)+(nx_-1)]
                    -stencil_.apply(solution_,NE,nx_-1,ny_-1,nx_,ny_);
        //south boarder
        for (Index sx=2; sx<nx_-1; sx++)
            result[1*(nx_+1)+sx]=
                rightHandSide_[1*(nx_+1)+sx]
                    -stencil_.apply(solution_,S,sx,1,nx_,ny_);
        //north boarder
        for (Index sx=2; sx<nx_-1; sx++)
            result[(nx_-1)*(nx_+1)+sx]=
                rightHandSide_[(nx_-1)*(nx_+1)+sx]
                    -stencil_.apply(solution_,N,sx,ny_-1,nx_,ny_);
        //west boarder
        for (Index sy=2; sy<ny_-1; sy++)
            result[sy*(nx_+1)+1]=
                rightHandSide_[sy*(nx_+1)+1]
                    -stencil_.apply(solution_,W,1,sy,nx_,ny_);
        //east boarder;
        for (Index sy=2; sy<ny_-1; sy++)
            result[sy*(nx_+1)+(nx_-1)]=
                rightHandSide_[sy*(nx_+1)+(nx_-1)]
                    -stencil_.apply(solution_,E,nx_-1,sy,nx_,ny_);
        //the center
        for (Index sy=2; sy<ny_-1; sy++)
            for (Index sx=2; sx<nx_-1; sx++)
                result[sy*(nx_+1)+sx]=
                    rightHandSide_[sy*(nx_+1)+sx]
                        -stencil_.apply(solution_,C,sx,sy,nx_,ny_);
    }
    return result;
}

Point DirichletProblem::getLowerLeftCorner() const
{
    return Point(0,0);
}

Point DirichletProblem::getUpperRightCorner() const
{
    return Point(nx_,ny_);
}

Point DirichletProblem::getFirstPoint() const
{
    return Point(1,1);
}

Point DirichletProblem::getLastPoint() const
{
    return Point(nx_-1,ny_-1);
}

}
