#ifndef NUMERICARRAY_H_
#define NUMERICARRAY_H_

#include <valarray>
#include "typedefs.h"

namespace mg
{

/**
 * \brief vector holding numeric values
 */
class NumericArray : public std::valarray<Precision>
{
public:
    NumericArray();
    NumericArray(Precision initialValue, Index size);
	NumericArray(Precision initialValue, Index nx, Index ny);
    NumericArray(const NumericArray& rhs);
    const NumericArray& operator =(const NumericArray& rhs);
    const NumericArray& operator =(Precision rhs);
    Precision& operator()(Integer sx, Integer sy);
    const Precision& operator()(Integer sx, Integer sy) const;
private:
    Index calculateIndex(Integer sx, Integer sy) const;
    const Index nx_;
    const Index ny_;
};

}

#endif /*NUMERICARRAY_H_*/
