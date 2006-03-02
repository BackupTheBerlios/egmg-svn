#ifndef WCYCLE_H_
#define WCYCLE_H_

#include "../general/parameters.h"
#include "../Stencil/Stencil.h"

namespace mg
{

class WCycle //: interface mg::CycleType
{
private:
    const size_t gamma_;
    size_t repeats_;
    const size_t maximalDepth_;
    const size_t currentDepth_;
public:
    WCycle(const size_t maximalDepth,const size_t gamma)
        : gamma_(gamma),
          repeats_(0),
          maximalDepth_(maximalDepth),
          currentDepth_(0)
    {}
    
    WCycle(const WCycle& rhs)
        : gamma_(rhs.gamma_),
          repeats_(0),
          maximalDepth_(rhs.maximalDepth_),
          currentDepth_(rhs.currentDepth_+1)
    {}

    virtual ~WCycle() {}
    
    virtual bool solve() const
    {
        return currentDepth_==maximalDepth_;
    }
    
    virtual bool repeat()
    {
        ++repeats_;
        return repeats_<=gamma_;
    }
    
    virtual void accelerate(
        NumericArray&,
        const NumericArray&,
        const Stencil&,
        const size_t,
        const size_t)
    {}
};

}

#endif /*WCYCLE_H_*/
