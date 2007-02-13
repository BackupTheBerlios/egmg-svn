#include "NumericArray.h"
#include <cassert>

namespace mg
{
    
NumericArray::NumericArray()
    : std::valarray<Precision>()
{
}
    
NumericArray::NumericArray(Precision initialValue, Index size)
    : std::valarray<Precision>(initialValue,size)
{
}

NumericArray::NumericArray(const NumericArray& rhs)
    : std::valarray<Precision>(rhs), nx_(rhs.nx_), ny_(rhs.ny_)
{
}

NumericArray::NumericArray(Precision initialValue, Index nx, Index ny)
    : std::valarray<Precision>(initialValue,(nx+3)*(ny+3)), nx_(nx), ny_(ny)
{
}

const NumericArray& NumericArray::operator =(const NumericArray& rhs)
{
    return std::valarray<Precision>::operator =(rhs);
}

const NumericArray& NumericArray::operator =(Precision rhs)
{
    return std::valarray<Precision>::operator =(rhs);
}

Precision& NumericArray::operator()(Integer sx, Integer sy)
{
    return operator[](calculateIndex(sx,sy));
}

const Precision& NumericArray::operator()(Integer sx, Integer sy) const
{
    
    return operator[](calculateIndex(sx,sy));
}

Index NumericArray::calculateIndex(Integer sx, Integer sy) const
{
    assert( sx+1 >= 0 && sy+1 >= 0 );
    Index sx_ = static_cast<Index>(sx+1);
    Index sy_ = static_cast<Index>(sy+1);
    assert( sx_ < nx_+3 && sy_ < ny_+3 );
    return sy_*(nx_+3)+sx_;
}

}
