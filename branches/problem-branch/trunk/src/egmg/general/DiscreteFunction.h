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
        Index nx,
        Index ny,
        Precision hx,
        Precision hy);
    DiscreteFunction(const DiscreteFunction& rhs);
    DiscreteFunction(const Function& function, Index nx, Index ny);
    DiscreteFunction(
        const Function& function,
        Point origin,
        Index nx,
        Index ny,
        Precision hx,
        Precision hy);
    const DiscreteFunction& operator =(const DiscreteFunction& rhs);
    const DiscreteFunction& operator =(const Function& rhs);
    const DiscreteFunction& operator =(Precision rhs);
    Precision& operator()(Integer sx, Integer sy);
    const Precision& operator()(Integer sx, Integer sy) const;
    void write(std::ostream& out) const;
    Precision twoNorm() const;
    Precision oneNorm() const;
    Precision maxNorm() const;
    const DiscreteFunction abs() const;
    Index getNx() const;
    Index getNy() const;
    Precision getHx() const;
    Precision getHy() const;
    Point getOrigin() const;
    bool checkSimilarity( const DiscreteFunction& rhs) const;

    const DiscreteFunction& operator +=(const DiscreteFunction rhs);
    const DiscreteFunction& operator -=(const DiscreteFunction rhs);
    const DiscreteFunction& operator +=(Precision rhs);
    const DiscreteFunction& operator -=(Precision rhs);
    const DiscreteFunction& operator *=(const DiscreteFunction rhs);
    const DiscreteFunction& operator /=(const DiscreteFunction rhs);
    const DiscreteFunction& operator *=(Precision rhs);
    const DiscreteFunction& operator /=(Precision rhs);

private:
    Index calculateIndex(Integer sx, Integer sy) const;
    Index nx_;
    Index ny_;
    Precision hx_;
    Precision hy_;
    Point origin_;
    NumericArray data_;
};

std::ostream& operator<<(
    std::ostream& stream,
    const DiscreteFunction& function);

const DiscreteFunction operator -(
    const DiscreteFunction& lhs,
    const DiscreteFunction& rhs);

const DiscreteFunction operator +(
    const DiscreteFunction& lhs,
    const DiscreteFunction& rhs);

const DiscreteFunction operator *(
    const DiscreteFunction& lhs,
    const DiscreteFunction& rhs);

const DiscreteFunction operator /(
    const DiscreteFunction& lhs,
    const DiscreteFunction& rhs);

const DiscreteFunction operator +(
    const Precision lhs,
    const DiscreteFunction& rhs);

const DiscreteFunction operator +(
    const DiscreteFunction& rhs,
    const Precision lhs);

const DiscreteFunction operator -(
    const Precision lhs,
    const DiscreteFunction& rhs);

const DiscreteFunction operator -(
    const DiscreteFunction& rhs,
    const Precision lhs);

const DiscreteFunction operator *(
    const Precision lhs,
    const DiscreteFunction& rhs);

const DiscreteFunction operator *(
    const DiscreteFunction& rhs,
    const Precision lhs);

const DiscreteFunction operator /(
    const Precision lhs,
    const DiscreteFunction& rhs);

}

#endif /*DISCRETEFUNCTION_H_*/
