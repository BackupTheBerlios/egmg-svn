#include "Problem.h"

namespace mg
{

Problem::Problem(Stencil& stencil,
                 Index nx,
                 Index ny)
    : stencil_(stencil),
      nx_(nx),
      ny_(ny)
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
    rightHandSide_.resize((nx_+1)*(ny_+1));
    const Precision hx=1.0/nx_;
    const Precision hy=1.0/ny_;
    for (Index sy=0; sy<=ny_; sy++)
        for (Index sx=0; sx<=nx_; sx++)
            rightHandSide_[sy*(nx_+1)+sx]=rightHandSide(sx*hx,sy*hy);
}

void Problem::setRightHandSide( const NumericArray& rightHandSide )
{
    rightHandSide_.resize(rightHandSide.size());
    rightHandSide_ = rightHandSide;
}

const NumericArray& Problem::getRightHandSide()
{
    return rightHandSide_;
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
