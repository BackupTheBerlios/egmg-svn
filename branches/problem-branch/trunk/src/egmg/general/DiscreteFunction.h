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
    const DiscreteFunction abs() const;
    Index getNx() const;
    Index getNy() const;
private:
    Index calculateIndex(Integer sx, Integer sy) const;
    Index nx_;
    Index ny_;
};

std::ostream& operator<<(std::ostream& stream, const DiscreteFunction& function);
const DiscreteFunction operator -(const DiscreteFunction& lhs, const DiscreteFunction& rhs);
const DiscreteFunction operator +(const DiscreteFunction& lhs, const DiscreteFunction& rhs);
const DiscreteFunction operator *(const DiscreteFunction& lhs, const DiscreteFunction& rhs);
const DiscreteFunction operator *(const Precision lhs, const DiscreteFunction& rhs);
const DiscreteFunction operator *(const DiscreteFunction& rhs, const Precision lhs);


}

#endif /*DISCRETEFUNCTION_H_*/
