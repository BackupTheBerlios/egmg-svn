/** \file WCycle.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>,
 *  <a href="mailto:kai.lists.berlios.de@tragetaschen.dyndns.org">Kai Ruhnau</a>
 * \brief Contains the class WCycle.
 * \see cycle.h
 */
#ifndef WCYCLE_H_
#define WCYCLE_H_

#include <vector>
#include "LevelCycle.h"

namespace mg
{

/**
 * \brief WCycle is a W Cylce
 */
class WCycle : public mg::LevelCycle
{
private:
    //we don't want the autogenerated copy constructor and assignment operator
    WCycle(const WCycle&);
    WCycle& operator=(const WCycle&);
public:
    /**
     * \brief The constructor of a WCycle object
     * 
     * WCycle constructs a WCycle object with:
     * \param[in] preSmoothingSteps     number of pre smoothing steps  (def. 1)
     * \param[in] postSmoothingSteps    number of post smoothing steps (def. 1)
     * \param[in] maximalDepth  Number of Grid Levels
     * \param[in] gamma         Number of repeats on each grid level
     */
    WCycle(
        const int preSmoothingSteps,
        const int postSmoothingSteps,
        const Index maximalDepth,
        const Index gamma)
        : LevelCycle(
                preSmoothingSteps,
                postSmoothingSteps,
                std::vector<Index>(maximalDepth,gamma))
    {}
};

}

#endif /*WCYCLE_H_*/