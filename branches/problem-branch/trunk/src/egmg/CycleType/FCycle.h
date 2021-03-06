/** \file FCycle.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>,
 *  <a href="mailto:kai.lists.berlios.de@tragetaschen.dyndns.org">Kai Ruhnau</a>
 * \brief Contains the class FCycle.
 * \see cycle.h
 */
#ifndef FCYCLE_H_
#define FCYCLE_H_

#include <vector>
#include "CycleType.h"

namespace mg
{

/**
 * \brief FCycle is a F Cylce
 */    
class FCycle : public mg::CycleType
{
private:
    const Index maximalDepth_;
    Index currentDepth_;
    enum State_{INIT,FCYCLE,VCYCLE};
    std::vector<State_ > currentState_;
    //we don't want the autogenerated copy constructor and assignment operator
    FCycle(const FCycle&);
    FCycle& operator=(const FCycle&);

public:
    /**
     * \brief The constructor of a FCycle object
     * 
     * FCycle constructs a FCycle object with:
     * \param[in] preSmoothingSteps     number of pre smoothing steps  (def. 1)
     * \param[in] postSmoothingSteps    number of post smoothing steps (def. 1)
     * \param[in] maximalDepth  Number of Grid Levels
     */
    FCycle(
        const int preSmoothingSteps,
        const int postSmoothingSteps,
        Index maximalDepth)
        : CycleType(preSmoothingSteps, postSmoothingSteps),
          maximalDepth_(maximalDepth),
          currentDepth_(0),
          currentState_(maximalDepth+1,INIT)
    {}
    
    virtual void incrementGridLevel()
    {
        ++currentDepth_;
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
        return maximalDepth_==currentDepth_;
    }
    
    /**
     * \brief repeat() the iteration on this level?
     * 
     * \return  true if a another iteration should be done on this grid level
     */
    virtual bool repeat()
    {
        if (currentState_[currentDepth_]==INIT)
            currentState_[currentDepth_]=FCYCLE;
        else if (currentState_[currentDepth_]==FCYCLE)
            currentState_[currentDepth_]=VCYCLE;
        else if (currentState_[currentDepth_]==VCYCLE)
        {
            //reinitilize for next cycle
            currentState_[currentDepth_]=INIT;
            return false;
        }
        return true;
    }
    
    /**
     * \brief accelerate() the current solution
     * 
     * does nothing for FCycle
     */
    virtual void accelerate(
        Problem&)
    {}
};
}

#endif /*FCYCLE_H_*/
