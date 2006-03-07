/** \file WCycle.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>,
 *  <a href="mailto:kai.lists.berlios.de@tragetaschen.dyndns.org">Kai Ruhnau</a>
 * \brief Contains the class WCycle.
 * \see cycle.h
 */
#ifndef WCYCLE_H_
#define WCYCLE_H_

#include <vector>
#include "CycleType.h"

namespace mg
{

/**
 * \brief WCycle is a W Cylce
 */
class WCycle : public mg::CycleType
{
private:
    const Index gamma_;
    std::vector<Index> repeats_;
    const Index maximalDepth_;
    Index currentDepth_;
public:
    /**
     * \brief The constructor of a WCycle object
     * 
     * WCycle constructs a WCycle object with:
     * \param[in] maximalDepth  Number of Grid Levels
     * \param[in] gamma         Number of repeats on each grid level
     */
    WCycle(const Index maximalDepth,const Index gamma)
        : gamma_(gamma),
          repeats_(maximalDepth+1,0),
          maximalDepth_(maximalDepth),
          currentDepth_(0)
    {}
    
    virtual ~WCycle() {}
    
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
        return currentDepth_==maximalDepth_;
    }
    
    /**
     * \brief repeat() the iteration on this level?
     * 
     * repeat() check if gamma iterations have been done on this grid level
     * 
     * \return  true if a another iteration should be done on this grid level
     */
    virtual bool repeat()
    {
        ++repeats_[currentDepth_];
        return repeats_[currentDepth_]<=gamma_;
    }
    
    /**
     * \brief accelerate() the current solution
     * 
     * does nothing for WCycle
     */
    virtual void accelerate(
        NumericArray&,
        const NumericArray&,
        const Stencil&,
        const Index,
        const Index)
    {}
};

}

#endif /*WCYCLE_H_*/
