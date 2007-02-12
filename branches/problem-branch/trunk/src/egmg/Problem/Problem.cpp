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
    rightHandSide_.resize((nx+1)*(ny+1));
    const Precision hx=1.0/nx;
    const Precision hy=1.0/ny;
    for (Index sy=0; sy<=ny; sy++)
        for (Index sx=0; sx<=nx; sx++)
            rightHandSide_[sy*(nx+1)+sx]=rightHandSide(sx*hx,sy*hy);
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

}
