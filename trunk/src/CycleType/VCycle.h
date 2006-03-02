#ifndef VCYCLE_H_
#define VCYCLE_H_

#include "WCycle.h"

namespace mg
{

class VCycle : public mg::WCycle
{
public:
	explicit VCycle(const size_t maximalDepth)
        : WCycle(maximalDepth,1)
    {}    
};

}

#endif /*VCYCLE_H_*/
