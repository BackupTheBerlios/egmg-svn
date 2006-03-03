/** \file WCycle.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>,
 *  <a href="mailto:kai.lists.berlios.de@tragetaschen.dyndns.org">Kai Ruhnau</a>
 * 
 * \brief Contains the class WCycle.
 * \see cycle.h
 */
#ifndef WCYCLE_H_
#define WCYCLE_H_

#include "../general/parameters.h"
#include "../Stencil/Stencil.h"

namespace mg
{

/**
 * \brief VCycle is a W Cylce
 */
class WCycle //: interface mg::CycleType
{
private:
    const Index gamma_;
    Index repeats_;
    const Index maximalDepth_;
    const Index currentDepth_;
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
          repeats_(0),
          maximalDepth_(maximalDepth),
          currentDepth_(0)
    {}
    
    /**
     * \brief The copy constructor of a WCycle object
     * 
     * WCycle constructs a WCycle object with the same parameters than rhs,
     * but a current level counter is increased by one.
     * \param[in] rhs   the WCycle to copy
     */
    WCycle(const WCycle& rhs)
        : gamma_(rhs.gamma_),
          repeats_(0),
          maximalDepth_(rhs.maximalDepth_),
          currentDepth_(rhs.currentDepth_+1)
    {}

    virtual ~WCycle() {}
    
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
        ++repeats_;
        return repeats_<=gamma_;
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
