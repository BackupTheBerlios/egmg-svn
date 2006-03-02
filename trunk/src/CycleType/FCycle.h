#ifndef FCYCLE_H_
#define FCYCLE_H_

#include "../general/parameters.h"
#include "../Stencil/Stencil.h"

namespace mg
{
    
class FCycle //: interface mg::CycleType
{
private:
    const size_t maximumDepth_;
    const size_t currentDepth_;
    enum State_{INIT,FCYCLE,VCYCLE};
    State_ currentState_;

public:
    FCycle(size_t maximumDepth)
        : maximumDepth_(maximumDepth),
          currentDepth_(0),
          currentState_(INIT)
    {}
    
    FCycle(const FCycle& rhs)
        : maximumDepth_(rhs.maximumDepth_),
          currentDepth_(rhs.currentDepth_+1),
          currentState_(rhs.currentState_)
    {}
    
    bool solve() const
    {
        return maximumDepth_==currentDepth_;
    }
    
    bool repeat()
    {
        State_ result=currentState_;
        if (currentState_==INIT)
            currentState_=FCYCLE;
        else if (currentState_==FCYCLE)
            currentState_=VCYCLE;
        return result!=VCYCLE;
    }
    
    void accelerate(
        NumericArray&,
        const NumericArray&,
        const Stencil&,
        const size_t,
        const size_t)
    {}
};
}

#endif /*FCYCLE_H_*/
