#include "Problem.h"

namespace mg
{

Problem::Problem(
        Stencil& stencil,
        Point origin,
        Index nx, Index ny,
        Precision hx, Precision hy)
    : stencil_(stencil),
      rightHandSide_(0.0,origin,nx,ny,hx,hy),
      solution_(0.0,origin,nx,ny,hx,hy)
{
}

Problem::~Problem()
{
}

Stencil& Problem::getStencil()
{
    return stencil_;
}

void Problem::setRightHandSide( const Function& rightHandSide )
{
    const Index nx = getNx();
    const Index ny = getNy();
    const Precision hx=getHx();
    const Precision hy=getHy();
    const Point origin = getOrigin();
    for (Index sy=0; sy<=ny; sy++)
        for (Index sx=0; sx<=nx; sx++)
            rightHandSide_(sx,sy)=rightHandSide(origin.x+sx*hx,origin.y+sy*hy);
}

void Problem::setRightHandSide( const DiscreteFunction& rightHandSide )
{
    ASSERT(solution_.checkSimilarity(rightHandSide));
    rightHandSide_ = rightHandSide;
}

const DiscreteFunction& Problem::getRightHandSide()
{
    return rightHandSide_;
}

DiscreteFunction& Problem::getSolution()
{
    return solution_;
}

const DiscreteFunction& Problem::getSolution() const
{
    return solution_;
}

IndexPair Problem::getFirstPoint() const
{
    return getFirstPoint(getNx(),getNy());
}

IndexPair Problem::getLastPoint() const
{
    return getLastPoint(getNx(),getNy());
}

Index Problem::getNx() const
{
    return solution_.getNx();
}

Index Problem::getNy() const
{
    return solution_.getNy();
}

Precision Problem::getHx() const
{
    return solution_.getHx();
}

Precision Problem::getHy() const
{
    return solution_.getHy();
}

Point Problem::getOrigin() const
{
    return solution_.getOrigin();
}

void Problem::fillRightHandSide(NumericArray& rightSide) const
{
    const Index nx = getNx();
    const Index ny = getNy();
    for(Index sy=0;sy<=ny;++sy)
        for(Index sx=0;sx<=nx;++sx)
            rightSide[sy*(nx+1)+sx]=rightHandSide_(sx,sy);
}

void Problem::pointFillMatrix(
    NumericArray& matrix,
    const Index dimension,
    const Position postion,
    const Index sx,
    const Index sy) const
{
    const Index nx = getNx();
    const Index ny = getNy();
    const Precision hx = getHx();
    const Precision hy = getHy();
    const Point origin = getOrigin();
    PositionArray jX=stencil_.getJx(postion,nx,ny);
    PositionArray jY=stencil_.getJy(postion,nx,ny);
    NumericArray operatorL=stencil_.getL(postion,sx,sy,nx,ny,hx,hy,origin);
    for (Index i=0; i<operatorL.size(); ++i)
        matrix[(sy*(nx+1)+sx)*dimension+((sy+jY[i])*(nx+1)+sx+jX[i])]
            =operatorL[i];
}

void Problem::fillMatrix(NumericArray& matrix, Index dimension) const
{
    const Index nx = getNx();
    const Index ny = getNy();
    //corner points
    pointFillMatrix(matrix,dimension,NW,1,ny-1);
    pointFillMatrix(matrix,dimension,NE,nx-1,ny-1);
    pointFillMatrix(matrix,dimension,SE,nx-1,1);
    pointFillMatrix(matrix,dimension,SW,1,1);
    //border points
    for (Index sy=2; sy<(ny-1); ++sy)
        pointFillMatrix(matrix,dimension,W,1,sy);
    for (Index sy=2; sy<(ny-1); ++sy)
        pointFillMatrix(matrix,dimension,E,nx-1,sy);
    for (Index sx=2; sx<(nx-1); ++sx)
        pointFillMatrix(matrix,dimension,N,sx,ny-1);
    for (Index sx=2; sx<(nx-1); ++sx)
        pointFillMatrix(matrix,dimension,S,sx,1);
    //center Points
    for (Index sy=2; sy<(ny-1); ++sy)
        for (Index sx=2; sx<(nx-1); ++sx)
            pointFillMatrix(matrix,dimension,C,sx,sy);
}

LinearEquationSystem Problem::getLinearEquationSystem() const
{
    const Index nx = getNx();
    const Index ny = getNy();
    const Index dimension=(nx+1)*(ny+1);
    NumericArray matrix(0.0,dimension*dimension);
    NumericArray rightSide(0.0,dimension);
    
    fillRightHandSide(rightSide);
    
    fillBorderValues(matrix,rightSide,dimension);

    fillMatrix(matrix,dimension);
        
    return std::make_pair(matrix,rightSide);
}

void Problem::setSolution(NumericArray& solution)
{
    const Index nx = getNx();
    const Index ny = getNy();
    for(Index sy=0;sy<=ny;++sy)
        for(Index sx=0;sx<=nx;++sx)
            solution_(sx,sy)=solution[sy*(nx+1)+sx];
    applyBoundaryConstraint();
}

}
