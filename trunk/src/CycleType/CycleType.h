#ifndef CYCLETYPE_H_
#define CYCLETYPE_H_

#include "../general/parameters.h"
#include "../Stencil/Stencil.h"

namespace mg
{

class CycleType
{
public:
    CycleType();
    
    CycleType(const CycleType& rhs);
    
    ~CycleType();
    
    bool solve() const;
    
    bool repeat();
    
    void accelerate(
        NumericArray& u,
        const NumericArray& f,
        const Stencil& stencil,
        const size_t nx,
        const size_t ny);
};

}

#endif /*CYCLETYPE_H_*/
