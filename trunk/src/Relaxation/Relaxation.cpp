/** \file Relaxation.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementaion of the class Relaxation.
 * \see Relaxation.h
 */
#include "Relaxation.h"

namespace mg
{

virtual void Relaxation::preSmooth(
    std::valarray<Precision>& u,
    const std::valarray<Precision>& f,
    const Stencil& stencil,
    const size_t nx,
    const size_t ny)
{
    for(int i=0; i<preSmoothingSteps_; ++i)
        relax(u,f,stencil,nx,ny);
}
        
virtual void Relaxation::postSmooth(
    std::valarray<Precision>& u,
    const std::valarray<Precision>& f,
    const Stencil& stencil,
    const size_t nx,
    const size_t ny)
{
    for(int i=0; i<postSmoothingSteps_; ++i)
        relax(u,f,stencil,nx,ny);
}
}
