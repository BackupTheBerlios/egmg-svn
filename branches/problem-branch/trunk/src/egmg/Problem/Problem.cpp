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

}
