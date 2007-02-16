#include "Problem.h"

namespace mg
{

Problem::Problem(Stencil& stencil,
                 Index nx,
                 Index ny)
    : stencil_(stencil),
      nx_(nx),
      ny_(ny),
      rightHandSide_(0.0,nx,ny),
      solution_(0.0,nx,ny)
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
    const Precision hx=1.0/nx_;
    const Precision hy=1.0/ny_;
    for (Index sy=0; sy<=ny_; sy++)
        for (Index sx=0; sx<=nx_; sx++)
            rightHandSide_(sx,sy)=rightHandSide(sx*hx,sy*hy);
}

void Problem::setRightHandSide( const DiscreteFunction& rightHandSide )
{
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

Index Problem::getNx() const
{
    return nx_;
}

Index Problem::getNy() const
{
    return ny_;
}

void Problem::fillRightHandSide(NumericArray& rightSide) const
{
    for(Index sy=0;sy<=ny_;++sy)
        for(Index sx=0;sx<=nx_;++sx)
            rightSide[sy*(nx_+1)+sx]=rightHandSide_(sx,sy);
}

void Problem::pointFillMatrix(
    NumericArray& matrix,
    const Index dimension,
    const Position postion,
    const Index sx,
    const Index sy) const
{
    PositionArray jX=stencil_.getJx(postion,nx_,ny_);
    PositionArray jY=stencil_.getJy(postion,nx_,ny_);
    NumericArray operatorL=stencil_.getL(postion,sx,sy,nx_,ny_);
    for (Index i=0; i<operatorL.size(); ++i)
        matrix[(sy*(nx_+1)+sx)*dimension+((sy+jY[i])*(nx_+1)+sx+jX[i])]
            =operatorL[i];
}

void Problem::fillMatrix(NumericArray& matrix, Index dimension) const
{
    //corner points
    pointFillMatrix(matrix,dimension,NW,1,ny_-1);
    pointFillMatrix(matrix,dimension,NE,nx_-1,ny_-1);
    pointFillMatrix(matrix,dimension,SE,nx_-1,1);
    pointFillMatrix(matrix,dimension,SW,1,1);
    //border points
    for (Index sy=2; sy<(ny_-1); ++sy)
        pointFillMatrix(matrix,dimension,W,1,sy);
    for (Index sy=2; sy<(ny_-1); ++sy)
        pointFillMatrix(matrix,dimension,E,nx_-1,sy);
    for (Index sx=2; sx<(nx_-1); ++sx)
        pointFillMatrix(matrix,dimension,N,sx,ny_-1);
    for (Index sx=2; sx<(nx_-1); ++sx)
        pointFillMatrix(matrix,dimension,S,sx,1);
    //center Points
    for (Index sy=2; sy<(ny_-1); ++sy)
        for (Index sx=2; sx<(nx_-1); ++sx)
            pointFillMatrix(matrix,dimension,C,sx,sy);
}

LinearEquationSystem Problem::getLinearEquationSystem() const
{
    const Index dimension=(nx_+1)*(ny_+1);
    NumericArray matrix(0.0,dimension*dimension);
    NumericArray rightSide(0.0,dimension);
    
    fillRightHandSide(rightSide);
    
    fillBorderValues(matrix,rightSide,dimension);

    fillMatrix(matrix,dimension);
        
    return std::make_pair(matrix,rightSide);
}

void Problem::setSolution(NumericArray& solution)
{
    for(Index sy=0;sy<=ny_;++sy)
        for(Index sx=0;sx<=nx_;++sx)
            solution_(sx,sy)=solution[sy*(nx_+1)+sx];
    applyBoundaryConstraint();
}

}
