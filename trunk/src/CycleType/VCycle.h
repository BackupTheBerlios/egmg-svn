/** \file VCycle.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>,
 *  <a href="mailto:kai.lists.berlios.de@tragetaschen.dyndns.org">Kai Ruhnau</a>
 * \brief Contains the class VCycle.
 * \see cycle.h
 */
#ifndef VCYCLE_H_
#define VCYCLE_H_

#include "WCycle.h"

namespace mg
{

/**
 * \brief VCycle is a V Cylce
 */
class VCycle : public mg::WCycle
{
public:
    /**
     * \brief The constructor of a VCycle object
     * 
     * VCycle constructs a VCycle object with:
     * \param[in] preSmoothingSteps     number of pre smoothing steps  (def. 1)
     * \param[in] postSmoothingSteps    number of post smoothing steps (def. 1)
     * \param[in] maximalDepth  Number of Grid Levels
     */
    VCycle(
        const int preSmoothingSteps,
        const int postSmoothingSteps,
        const Index maximalDepth)
        : WCycle(preSmoothingSteps, postSmoothingSteps,maximalDepth,1)
    {}    
};

}

#endif /*VCYCLE_H_*/
