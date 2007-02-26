#include "DirichletProblem.h"

namespace mg
{

DirichletProblem::DirichletProblem(
    Stencil& stencil,
        Point origin,
        Index nx, Index ny,
        Precision hx, Precision hy)
    : Problem(stencil, origin, nx, ny, hx, hy)
{
}

DirichletProblem::~DirichletProblem()
{
}

void DirichletProblem::setBoundaryConstraint( const Function& boundaryConstraint )
{
    const Index nx = getNx();
    const Index ny = getNy();
    const Precision hx=getHx();
    const Precision hy=getHy();
    const Point origin = getOrigin();
    solution_(-1,-1)=boundaryConstraint(origin.x-hx,origin.y-hy);
    for (Index sx=0; sx<=nx; sx++)
    {
        solution_(sx,-1)=boundaryConstraint(origin.x+sx*hx,origin.y-hy);
        solution_(sx,0)=boundaryConstraint(origin.x+sx*hx,origin.y+0.0);       //top border
        solution_(sx,ny)=boundaryConstraint(origin.x+sx*hx,origin.y+ny*hy);     //bottom border
        solution_(sx,ny+1)=boundaryConstraint(origin.x+sx*hx,origin.y+(ny+1)*hy);
    }
    solution_(nx+1,-1)=boundaryConstraint(origin.x+(nx+1)*hx,origin.y-hy);
    solution_(-1,ny+1)=boundaryConstraint(origin.x-hx,origin.y+(ny+1)*hy);
    for (Index sy=1; sy<ny; sy++)
    {
        solution_(-1,sy)=boundaryConstraint(origin.x-hx,origin.y+sy*hy);
        solution_(0,sy)=boundaryConstraint(origin.x+0.0,origin.y+sy*hy);        //left border
        solution_(nx,sy)=boundaryConstraint(origin.x+nx*hx,origin.y+sy*hy);      //right border
        solution_(nx+1,sy)=boundaryConstraint(origin.x+(nx+1)*hx,origin.y+sy*hy);
    }
    solution_(nx+1,ny+1)=boundaryConstraint(origin.x+(nx+1)*hx,origin.y+(ny+1)*hy);
}

void DirichletProblem::applyBoundaryConstraint()
{
}

void DirichletProblem::applyBoundaryConstraint( DiscreteFunction& ) const
{
}

DiscreteFunction DirichletProblem::residuum()
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
        {
            for(Index sx=1; sx<nx; sx++)
            {
                result(sx,sy)=rightHandSide_(sx,sy)
                            -stencil_.apply(solution_,C,sx,sy);
            }
        }
    }
    else
    {
        //south west corner
        result(1,1)=rightHandSide_(1,1)
                    -stencil_.apply(solution_,SW,1,1);
        //south east corner
        result(nx-1,1)=rightHandSide_(nx-1,1)
                    -stencil_.apply(solution_,SE,nx-1,1);
        //north west corner
        result(1,ny-1)=rightHandSide_(1,ny-1)
                    -stencil_.apply(solution_,NW,1,ny-1);
        //north east corner
        result(nx-1,ny-1)=rightHandSide_(nx-1,ny-1)
                    -stencil_.apply(solution_,NE,nx-1,ny-1);
        //south boarder
        for (Index sx=2; sx<nx-1; sx++)
            result(sx,1)=rightHandSide_(sx,1)
                    -stencil_.apply(solution_,S,sx,1);
        //north boarder
        for (Index sx=2; sx<nx-1; sx++)
            result(sx,ny-1)=rightHandSide_(sx,ny-1)
                    -stencil_.apply(solution_,N,sx,ny-1);
        //west boarder
        for (Index sy=2; sy<ny-1; sy++)
            result(1,sy)=rightHandSide_(1,sy)
                    -stencil_.apply(solution_,W,1,sy);
        //east boarder;
        for (Index sy=2; sy<ny-1; sy++)
            result(nx-1,sy)=rightHandSide_(nx-1,sy)
                    -stencil_.apply(solution_,E,nx-1,sy);
        //the center
        for (Index sy=2; sy<ny-1; sy++)
            for (Index sx=2; sx<nx-1; sx++)
                result(sx,sy)=rightHandSide_(sx,sy)
                        -stencil_.apply(solution_,C,sx,sy);
    }
    return result;
}

IndexPair DirichletProblem::getFirstPoint(Index, Index) const
{
    return IndexPair(1,1);
}

IndexPair DirichletProblem::getLastPoint(Index nx, Index ny) const
{
    return IndexPair(nx-1,ny-1);
}

ProblemPtr DirichletProblem::getCoarsGridProblem(
    Index nxNew,
    Index nyNew,
    Precision hxNew,
    Precision hyNew) const
{
    return ProblemPtr(
        new DirichletProblem(stencil_,getOrigin(),nxNew,nyNew,hxNew,hyNew));
}

void DirichletProblem::fillBorderValues(
        NumericArray& matrix,
        NumericArray& rightSide,
        const Index dimension) const
{
    const Index nx = getNx();
    const Index ny = getNy();
    //border values XDIR
    for (Index sx=0; sx<=nx; ++sx)
    {
        //lower border
        matrix[sx*dimension+sx]=1;
        rightSide[sx]=solution_(sx,0);
        //upper border
        matrix[(dimension-1-sx)*dimension+(dimension-1-sx)]=1;
        rightSide[dimension-1-sx]=solution_(nx-sx,ny);
    }
    //border values YDIR
    //corners have been process in XDIR (y=1..ny-1 instead of y=0..ny)
    for (Index sy=1; sy<ny; ++sy)
    {
        //left border
        matrix[sy*(nx+1)*dimension+sy*(nx+1)]=1;
        rightSide[sy*(nx+1)]=solution_(0,sy);
        //right border
        matrix[(sy*(nx+1)+nx)*dimension+(sy*(nx+1)+nx)]=1;
        rightSide[sy*(nx+1)+nx]=solution_(nx,sy);
    }
}

}
