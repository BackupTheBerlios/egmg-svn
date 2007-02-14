#include "DiscreteFunction.h"
#include <cassert>
#include <iomanip>
#include <cmath>

namespace mg
{
    
DiscreteFunction::DiscreteFunction()
    : NumericArray(), nx_(0), ny_(0)
{
}

DiscreteFunction::DiscreteFunction(const DiscreteFunction& rhs)
    : NumericArray(rhs), nx_(rhs.nx_), ny_(rhs.ny_)
{
}

DiscreteFunction::DiscreteFunction(Precision initialValue, Index nx, Index ny)
    : NumericArray(initialValue,(nx+3)*(ny+3)), nx_(nx), ny_(ny)
{
}

const DiscreteFunction& DiscreteFunction::operator =(const DiscreteFunction& rhs)
{
    NumericArray::resize(rhs.size());
    NumericArray::operator =(rhs);
    nx_ = rhs.nx_;
    ny_ = rhs.ny_;
    return *this;
}

Precision& DiscreteFunction::operator()(Integer sx, Integer sy)
{
    return operator[](calculateIndex(sx,sy));
}

const Precision& DiscreteFunction::operator()(Integer sx, Integer sy) const
{
    
    return operator[](calculateIndex(sx,sy));
}

Index DiscreteFunction::calculateIndex(Integer sx, Integer sy) const
{
    assert( sx+1 >= 0 && sy+1 >= 0 );
    Index sx_ = static_cast<Index>(sx+1);
    Index sy_ = static_cast<Index>(sy+1);
    assert( sx_ < nx_+3 && sy_ < ny_+3 );
    return sy_*(nx_+3)+sx_;
}

void DiscreteFunction::write(std::ostream& out) const
{
    out<<"#begin points"<<std::endl;
    out<<std::setw(10)<<std::left<<"#x"<<" "
       <<std::setw(10)<<std::left<<"y"<<" "
       <<std::setw(10)<<std::left<<"value"<<std::endl;
    Precision hx=1.0/nx_;
    Precision hy=1.0/ny_;
    for (Index sy=0; sy<=ny_; ++sy)
        for (Index sx=0; sx<=nx_; ++sx)
            out<<std::setw(10)<<std::left<<sx*hx<<" "
               <<std::setw(10)<<std::left<<sy*hy<<" "
               <<std::setw(10)<<std::left<<operator[](calculateIndex(sx,sy))
               <<std::endl;
    out<<"#end points"<<std::endl;
}

Precision DiscreteFunction::twoNorm() const
{
    Precision result = 0.0;
    for (Index sy=1; sy<ny_; ++sy)
        for (Index sx=1; sx<nx_; ++sx)
        {
            Precision temp = operator[](calculateIndex(sx,sy));
            result+=temp*temp;
        }
    return std::sqrt(result);
}

std::ostream& operator<<(std::ostream& stream, const DiscreteFunction& function)
{
    function.write(stream);
    return stream;
}

}
