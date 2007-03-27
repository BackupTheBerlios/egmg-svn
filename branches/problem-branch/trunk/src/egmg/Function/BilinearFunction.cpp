#include "BilinearFunction.h"

namespace mg
{

BilinearFunction::BilinearFunction(): a_(initA())
{}

BilinearFunction::BilinearFunction(const NumericArray& a): a_(a)
{}

Precision BilinearFunction::apply(Precision x, Precision y) const
{
	return a_[0]+a_[1]*x+a_[2]*y+a_[3]*x*y;
}

NumericArray BilinearFunction::initA() const
{
    const Precision t[]={1.0,1.0,1.0,1.0};
    return NumericArray(t,4);
}

}
