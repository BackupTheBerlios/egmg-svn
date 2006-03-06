/** \file FCycle.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>,
 *  <a href="mailto:kai.lists.berlios.de@tragetaschen.dyndns.org">Kai Ruhnau</a>
 * 
 * \brief Contains the class FCycle.
 * \see cycle.h
 */
#ifndef FCYCLE_H_
#define FCYCLE_H_

#include "../general/parameters.h"
#include "../Stencil/Stencil.h"

namespace mg
{
    
class FCycle //: interface mg::CycleType
{
private:
    const Index maximumDepth_;
    const Index currentDepth_;
    enum State_{INIT,FCYCLE,VCYCLE};
    State_ currentState_;
    State_ lastState_;

public:
    /**
     * \brief The constructor of a FCycle object
     * 
     * FCycle constructs a FCycle object with:
     * \param[in] maximalDepth  Number of Grid Levels
     */
    FCycle(Index maximumDepth)
        : maximumDepth_(maximumDepth),
          currentDepth_(0),
          currentState_(INIT),
          lastState_(INIT)
    {}
    
    /**
     * \brief The copy constructor of a FCycle object
     * 
     * FCycle constructs a FCycle object with the same parameters than rhs,
     * but a current level counter is increased by one.
     * \param[in] rhs   the FCycle to copy
     */
    FCycle(const FCycle& rhs)
        : maximumDepth_(rhs.maximumDepth_),
          currentDepth_(rhs.currentDepth_+1),
          currentState_(rhs.lastState_),
          lastState_(rhs.lastState_)
    {}
    
    /**
     * \brief solve() on this grid level directly?
     * 
     * solve() checks if we are on the coarsest grid level.
     * 
     * \return  true if we are on the coarsest grid
     */
    bool solve() const
    {
        return maximumDepth_==currentDepth_;
    }
    
    /**
     * \brief repeat() the iteration on this level?
     * 
     * \return  true if a another iteration should be done on this grid level
     */
    bool repeat()
    {
        lastState_=currentState_;
        if (currentState_==INIT)
            currentState_=FCYCLE;
        else if (currentState_==FCYCLE)
            currentState_=VCYCLE;
        return lastState_!=VCYCLE;
    }
    
    /**
     * \brief accelerate() the current solution
     * 
     * does nothing for FCycle
     */
    void accelerate(
        NumericArray&,
        const NumericArray&,
        const Stencil&,
        const Index,
        const Index)
    {}
};
}

#endif /*FCYCLE_H_*/
