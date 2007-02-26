#include "PeriodicProblem.h"

namespace mg
{

PeriodicProblem::PeriodicProblem(
    Stencil& stencil,
    Point origin,
    Index nx, Index ny,
    Precision hx, Precision hy)
    : Problem(stencil, origin, nx, ny, hx, hy)
{
}

PeriodicProblem::~PeriodicProblem()
{
}

void PeriodicProblem::setBoundaryConstraint( const Function& )
{
}

void PeriodicProblem::applyBoundaryConstraint()
{
    applyBoundaryConstraint(solution_);
}

void PeriodicProblem::applyBoundaryConstraint( DiscreteFunction& u ) const
{
    const Index nx = getNx();
    const Index ny = getNy();
    for(Index sx=0;sx<=nx+1;++sx)
    {
        u(sx,0)=u(sx,ny);
        u(sx,ny+1)=u(sx,1);
    }
    for(Index sy=0;sy<=ny+1;++sy)
    {
        u(0,sy)=u(nx,sy);
        u(nx+1,sy)=u(1,sy);
    }
}

DiscreteFunction PeriodicProblem::residuum()
{
    const Index nx = getNx();
    const Index ny = getNy();
    const Precision hx=getHx();
    const Precision hy=getHy();
    const Point origin = getOrigin();
    DiscreteFunction result(0.0,origin,nx,ny,hx,hy);
    if (stencil_.size()<2)
    {
        for (Index sy=1; sy<=ny; sy++)
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

IndexPair PeriodicProblem::getFirstPoint(Index, Index) const
{
    return IndexPair(1,1);
}

IndexPair PeriodicProblem::getLastPoint(Index nx, Index ny) const
{
    return IndexPair(nx,ny);
}

ProblemPtr PeriodicProblem::getCoarsGridProblem(
    Index nxNew,
    Index nyNew,
    Precision hxNew,
    Precision hyNew) const
{
    return ProblemPtr(
        new PeriodicProblem(stencil_,getOrigin(),nxNew,nyNew,hxNew,hyNew));
}

void PeriodicProblem::fillBorderValues(
        NumericArray& matrix,
        NumericArray& rightSide,
        const Index dimension) const
{
    const Index nx = getNx();
    const Index ny = getNy();
    for (Index sx=0; sx<=nx; ++sx)
    {
        matrix[(0*(nx+1)+sx)*dimension+(0*(nx+1)+sx)]=1.0;
        matrix[(0*(nx+1)+sx)*dimension+(ny*(nx+1)+sx)]=-1.0;
        rightSide[(0*(nx+1)+sx)]=0.0;
        
        matrix[(ny*(nx+1)+sx)*dimension+(ny*(nx+1)+sx)]=-2.0;
        matrix[(ny*(nx+1)+sx)*dimension+(1*(nx+1)+sx)]=1.0;
        matrix[(ny*(nx+1)+sx)*dimension+((ny-1)*(nx+1)+sx)]+=1.0;
        rightSide[(ny*(nx+1)+sx)]=0.0;
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
        
        //matrix[(sy*(nx_+1)+nx_)*dimension+sy*(nx_+1)]=1;
        //rightSide[(sy*(nx_+1)+nx_)]=solution_(nx_,sy);
    }
}

}
