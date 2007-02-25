#ifndef DISCRETEFUNCTION_H_
#define DISCRETEFUNCTION_H_

#include "parameters.h"
#include "../Function/Function.h"
#include <ostream>

namespace mg
{

class DiscreteFunction
{
public:
    DiscreteFunction();
	DiscreteFunction(Precision initialValue, Index nx, Index ny);
    DiscreteFunction(
        Precision initialValue,
        Point origin,
        Precision hx,
        Precision hy,
        Index nx,
        Index ny);
    DiscreteFunction(const DiscreteFunction& rhs);
    DiscreteFunction(const Function& function, Index nx, Index ny);
    DiscreteFunction(
        const Function& function,
        Point origin,
        Precision hx,
        Precision hy,
        Index nx,
        Index ny);
    const DiscreteFunction& operator =(const DiscreteFunction& rhs);
    Precision& operator()(Integer sx, Integer sy);
    const Precision& operator()(Integer sx, Integer sy) const;
    void write(std::ostream& out) const;
    Precision twoNorm() const;
    const DiscreteFunction abs() const;
    Index getNx() const;
    Index getNy() const;
    Precision getHx() const;
    Precision getHy() const;
    Point getOrigin() const;
    
    const DiscreteFunction operator +=(const DiscreteFunction rhs);
    const DiscreteFunction operator -=(const DiscreteFunction rhs);
    const DiscreteFunction operator +=(Precision rhs);
    const DiscreteFunction operator -=(Precision rhs);
    const DiscreteFunction operator *=(const DiscreteFunction rhs);
    const DiscreteFunction operator /=(const DiscreteFunction rhs);
    const DiscreteFunction operator *=(Precision rhs);
    const DiscreteFunction operator /=(Precision rhs);
    
private:
    Index calculateIndex(Integer sx, Integer sy) const;
    bool checkSimilarity( const DiscreteFunction& rhs) const;
    Index nx_;
    Index ny_;
    Precision hx_;
    Precision hy_;
    Point origin_;
    NumericArray data_;
};

std::ostream& operator<<(std::ostream& stream, const DiscreteFunction& function);
const DiscreteFunction operator -(const DiscreteFunction& lhs, const DiscreteFunction& rhs);
const DiscreteFunction operator +(const DiscreteFunction& lhs, const DiscreteFunction& rhs);
const DiscreteFunction operator *(const DiscreteFunction& lhs, const DiscreteFunction& rhs);
const DiscreteFunction operator *(const Precision lhs, const DiscreteFunction& rhs);
const DiscreteFunction operator *(const DiscreteFunction& rhs, const Precision lhs);


}

#endif /*DISCRETEFUNCTION_H_*/
