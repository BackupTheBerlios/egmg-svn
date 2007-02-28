#include "RossbyHaurwitzWave.h"
#include <cmath>

namespace mg
{
RossbyHaurwitzWave::RossbyHaurwitzWave(
    Precision t /*=0.0*/,
    Precision c /*=1.0*/,
    Precision beta /*=1.0*/,
    Precision k /*=1*/,
    Precision l /*=1*/)
        : t_(t), c_(c), beta_(beta), k_(k), l_(l)
{}

void RossbyHaurwitzWave::setT(Precision t)
{
    t_=t;
}

Precision RossbyHaurwitzWave::getT()
{
    return t_;
}

Precision RossbyHaurwitzWave::operator() (
    Precision x, Precision y, Precision t) const
{
    return c_*std::sin(k_*x+(beta_*k_)/(k_*k_+l_*l_)*t)
                *std::sin(l_*y);
}

Precision RossbyHaurwitzWave::apply(Precision x, Precision y) const
{
    return operator()(x,y,t_);
}

}
