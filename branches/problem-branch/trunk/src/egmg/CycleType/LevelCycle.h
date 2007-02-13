/** \file LevelCycle.h
 * \author 
 *  <a href="mailto:kai.lists.berlios.de@tragetaschen.dyndns.org">Kai Ruhnau</a>
 * \brief Contains the class LevelCycle.
 * \see cycle.h
 */
#ifndef LEVELCYCLE_H_
#define LEVELCYCLE_H_

#include <vector>
#include "CycleType.h"
#include <iostream>

namespace mg
{

/**
 * \brief LevelCycle is a WCylce with different gamma for each level
 */
class LevelCycle : public mg::CycleType
{
private:
    std::vector<Index> gamma_;
    std::vector<Index> repeats_;
    Index currentDepth_;
public:
    /**
     * \brief The constructor of a LevelCycle object
     * 
     * LevelCycle constructs a LevelCycle object with:
     * \param[in] preSmoothingSteps     number of pre smoothing steps  (def. 1)
     * \param[in] postSmoothingSteps    number of post smoothing steps (def. 1)
     * \param[in] gamma         Number of repeats on different grid levels
     */
    LevelCycle(
        const int preSmoothingSteps,
        const int postSmoothingSteps,
        const std::vector<Index>& gamma)
        : CycleType(preSmoothingSteps, postSmoothingSteps),
          gamma_(gamma),
          repeats_(gamma.size()+1,0),
          currentDepth_(0)
    {}
    
    virtual ~LevelCycle() {}
    
    virtual void incrementGridLevel()
    {
        repeats_[++currentDepth_]=0;
    }
    
    virtual void decrementGridLevel()
    {
        --currentDepth_;
    }
    
    /**
     * \brief solve() on this grid level directly?
     * 
     * solve() checks if we are on the coarsest grid level.
     * 
     * \return  true if we are on the coarsest grid
     */
    virtual bool solve() const
    {
        return currentDepth_==gamma_.size();
    }
    
    /**
     * \brief repeat() the iteration on this level?
     * 
     * repeat() check if appropriate iterations have been done on 
     *          this grid level
     * 
     * \return  true if a another iteration should be done on this grid level
     */
    virtual bool repeat()
    {
        ++repeats_[currentDepth_];
        if (currentDepth_==1)
            return repeats_[currentDepth_]<=1;
        else
            return repeats_[currentDepth_]<=gamma_[currentDepth_];
    }
    
    /**
     * \brief accelerate() the current solution
     * 
     * does nothing for LevelCycle
     */
    virtual void accelerate(
        Problem&)
    {}
};

}

#endif /*LevelCYCLE_H_*/
