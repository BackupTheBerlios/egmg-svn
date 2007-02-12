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
    const Precision hx = 1.0/nx;
    const Precision hy = 1.0/ny;
    for (Index sx=0; sx<=nx; sx++)
    {
        solution_[sx]=boundaryConstraint(sx*hx,0);                     //top border
        solution_[u.size()-1-sx]=boundaryConstraint((nx-sx)*hx,1);     //bottom border
    }
    for (Index sy=0; sy<=ny; sy++)
    {
        solution_[sy*(nx+1)]=boundaryConstraint(0,sy*hy);              //left border
        solution_[(sy+1)*(nx+1)-1]=boundaryConstraint(1,sy*hy);        //right border
    }
}

void DirichletProblem::applyBoundaryConstraint()
{
}

NumericArray& DirichletProblem::getSolution()
{
    return solution_;
}

NumericArray DirichletProblem::residuum()
{
    NumericArray result(0.0,solution_.size());
    if (stencil_.size()<2)
    {
        for (Index sy=1; sy<ny; sy++)
            for(Index sx=1; sx<nx; sx++)
                result[sy*(nx+1)+sx]=
                        rightHandSide_[sy*(nx+1)+sx]
                            -stencil_.apply(solution_,C,sx,sy,nx,ny);
    }
    else
    {
        //south west corner
        result[1*(nx+1)+1]=
                rightHandSide_[1*(nx+1)+1]
                    -stencil.apply(solution_,SW,1,1,nx,ny);
        //south east corner
        result[1*(nx+1)+(nx-1)]=
                rightHandSide_[1*(nx+1)+(nx-1)]
                    -stencil.apply(solution_,SE,nx-1,1,nx,ny);
        //north west corner
        result[(nx-1)*(nx+1)+1]=
                rightHandSide_[(nx-1)*(nx+1)+1]
                    -stencil.apply(solution_,NW,1,ny-1,nx,ny);
        //north east corner
        result[(nx-1)*(nx+1)+(nx-1)]=
                rightHandSide_[(nx-1)*(nx+1)+(nx-1)]
                    -stencil.apply(solution_,NE,nx-1,ny-1,nx,ny);
        //south boarder
        for (Index sx=2; sx<nx-1; sx++)
            result[1*(nx+1)+sx]=
                rightHandSide_[1*(nx+1)+sx]
                    -stencil.apply(solution_,S,sx,1,nx,ny);
        //north boarder
        for (Index sx=2; sx<nx-1; sx++)
            result[(nx-1)*(nx+1)+sx]=
                rightHandSide_[(nx-1)*(nx+1)+sx]
                    -stencil.apply(solution_,N,sx,ny-1,nx,ny);
        //west boarder
        for (Index sy=2; sy<ny-1; sy++)
            result[sy*(nx+1)+1]=
                rightHandSide_[sy*(nx+1)+1]
                    -stencil.apply(solution_,W,1,sy,nx,ny);
        //east boarder;
        for (Index sy=2; sy<ny-1; sy++)
            result[sy*(nx+1)+(nx-1)]=
                rightHandSide_[sy*(nx+1)+(nx-1)]
                    -stencil.apply(solution_,E,nx-1,sy,nx,ny);
        //the center
        for (Index sy=2; sy<ny-1; sy++)
            for (Index sx=2; sx<nx-1; sx++)
                result[sy*(nx+1)+sx]=
                    rightHandSide_[sy*(nx+1)+sx]
                        -stencil.apply(solution_,C,sx,sy,nx,ny);
    }
    return result;
}

}
