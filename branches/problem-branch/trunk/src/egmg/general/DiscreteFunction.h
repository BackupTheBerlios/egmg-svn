#ifndef DISCRETEFUNCTION_H_
#define DISCRETEFUNCTION_H_

#include "parameters.h"
#include "../Function/Function.h"
#include <ostream>

namespace mg
{

class DiscreteFunction : public NumericArray
{
public:
    DiscreteFunction();
	DiscreteFunction(Precision initialValue, Index nx, Index ny);
    DiscreteFunction(const DiscreteFunction& rhs);
    DiscreteFunction(const Function& function, Index nx, Index ny);
    const DiscreteFunction& operator =(const DiscreteFunction& rhs);
    Precision& operator()(Integer sx, Integer sy);
    const Precision& operator()(Integer sx, Integer sy) const;
    void write(std::ostream& out) const;
    Precision twoNorm() const;
private:
    Index calculateIndex(Integer sx, Integer sy) const;
    Index nx_;
    Index ny_;
};

std::ostream& operator<<(std::ostream& stream, const DiscreteFunction& function);

}

#endif /*DISCRETEFUNCTION_H_*/
