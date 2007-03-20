#include "DiscreteFunction.h"
#include <iomanip>
#include <cmath>

namespace mg
{

DiscreteFunction::DiscreteFunction()
    : nx_(0), ny_(0), hx_(0.0), hy_(0.0), origin_(0.0,0.0)
{
}

DiscreteFunction::DiscreteFunction(const DiscreteFunction& rhs)
    : nx_(rhs.nx_), ny_(rhs.ny_),
      hx_(rhs.hx_), hy_(rhs.hy_),
      origin_(rhs.origin_), data_(rhs.data_)
{
}

DiscreteFunction::DiscreteFunction(Precision initialValue, Index nx, Index ny)
    : nx_(nx), ny_(ny),
      hx_(1.0/nx), hy_(1.0/ny),
      origin_(0.0,0.0),
      data_(initialValue,(nx+3)*(ny+3))
{
}

DiscreteFunction::DiscreteFunction(
        Precision initialValue,
        Point origin,
        Index nx,
        Index ny,
        Precision hx,
        Precision hy)
    : nx_(nx), ny_(ny),
      hx_(hx), hy_(hy),
      origin_(origin),
      data_(initialValue,(nx+3)*(ny+3))
{
}

DiscreteFunction::DiscreteFunction(const Function& function, Index nx, Index ny)
    : nx_(nx), ny_(ny),
      hx_(1.0/nx), hy_(1.0/ny),
      origin_(0.0,0.0), data_(0.0,(nx+3)*(ny+3))
{
    for (Integer sy=-1; sy<=static_cast<Integer>(ny_+1); ++sy)
        for (Integer sx=-1; sx<=static_cast<Integer>(nx_+1); ++sx)
               data_[calculateIndex(sx,sy)]=
                    function(origin_.x+sx*hx_,origin_.y+sy*hy_);
}

DiscreteFunction::DiscreteFunction(
        const Function& function,
        Point origin,
        Index nx,
        Index ny,
        Precision hx,
        Precision hy)
    : nx_(nx), ny_(ny),
      hx_(hx), hy_(hy),
      origin_(origin), data_(0.0,(nx+3)*(ny+3))
{
    for (Integer sy=-1; sy<=static_cast<Integer>(ny_+1); ++sy)
        for (Integer sx=-1; sx<=static_cast<Integer>(nx_+1); ++sx)
               data_[calculateIndex(sx,sy)]=
                    function(origin_.x+sx*hx_,origin_.y+sy*hy_);
}

const DiscreteFunction& DiscreteFunction::operator =(
    const DiscreteFunction& rhs)
{
    data_.resize(rhs.data_.size());
    data_ = rhs.data_;
    nx_ = rhs.nx_;
    ny_ = rhs.ny_;
    hx_ = rhs.hx_;
    hy_ = rhs.hy_;
    origin_ = rhs.origin_;
    return *this;
}

const DiscreteFunction& DiscreteFunction::operator =(const Function& rhs)
{
    for (Integer sy=-1; sy<=static_cast<Integer>(ny_+1); ++sy)
        for (Integer sx=-1; sx<=static_cast<Integer>(nx_+1); ++sx)
               data_[calculateIndex(sx,sy)]=
                    rhs(origin_.x+sx*hx_,origin_.y+sy*hy_);
    return *this;
}

Precision& DiscreteFunction::operator()(Integer sx, Integer sy)
{
    return data_[calculateIndex(sx,sy)];
}

const Precision& DiscreteFunction::operator()(Integer sx, Integer sy) const
{
    return data_[calculateIndex(sx,sy)];
}

Index DiscreteFunction::calculateIndex(Integer sx, Integer sy) const
{
    ASSERT( sx+1 >= 0 && sy+1 >= 0 );
    Index sx_ = static_cast<Index>(sx+1);
    Index sy_ = static_cast<Index>(sy+1);
    ASSERT( sx_ < nx_+3 && sy_ < ny_+3 );
    return sy_*(nx_+3)+sx_;
}

bool DiscreteFunction::checkSimilarity( const DiscreteFunction& rhs) const
{
    return nx_ == rhs.nx_ && ny_ == rhs.ny_ &&
           hx_ == rhs.hx_ && hy_ == rhs.hy_ &&
           origin_ == rhs.origin_;
}

void DiscreteFunction::write(std::ostream& out) const
{
    out<<"#begin points"<<std::endl;
    out<<std::setw(10)<<std::left<<"#x"<<" "
       <<std::setw(10)<<std::left<<"y"<<" "
       <<std::setw(10)<<std::left<<"value"<<std::endl;
    //for (Integer sy=-1; sy<=static_cast<Integer>(ny_+1); ++sy)
    for (Index sy=0; sy<=ny_; ++sy)
        //for (Integer sx=-1; sx<=static_cast<Integer>(nx_+1); ++sx)
        for (Index sx=0; sx<=nx_; ++sx)
            out<<std::setw(10)<<std::left<<origin_.x + sx*hx_<<" "
               <<std::setw(10)<<std::left<<origin_.y + sy*hy_<<" "
               <<std::setw(10)<<std::left<<data_[calculateIndex(sx,sy)]
               <<std::endl;
    out<<"#end points"<<std::endl;
}

Precision DiscreteFunction::twoNorm() const
{
    Precision result = 0.0;
    for (Index sy=0; sy<=ny_; ++sy)
        for (Index sx=0; sx<=nx_; ++sx)
        {
            Precision temp = data_[calculateIndex(sx,sy)];
            result+=temp*temp;
        }
    return std::sqrt(result/((nx_+1)*(ny_+1)));
}

const DiscreteFunction DiscreteFunction::abs() const
{
    DiscreteFunction result(*this);
    for (Index sy=0; sy<=ny_; ++sy)
        for (Index sx=0; sx<=nx_; ++sx)
            result(sx,sy) = std::abs(data_[calculateIndex(sx,sy)]);
    return result;
}

Index DiscreteFunction::getNx() const
{
    return nx_;
}

Index DiscreteFunction::getNy() const
{
    return ny_;
}

Precision DiscreteFunction::getHx() const
{
    return hx_;
}

Precision DiscreteFunction::getHy() const
{
    return hy_;
}

Point DiscreteFunction::getOrigin() const
{
    return origin_;
}

const DiscreteFunction& DiscreteFunction::operator +=(const DiscreteFunction rhs)
{
    ASSERT( checkSimilarity(rhs) );
    data_+=rhs.data_;
    return *this;

}
const DiscreteFunction& DiscreteFunction::operator -=(const DiscreteFunction rhs)
{
    ASSERT( checkSimilarity(rhs) );
    data_-=rhs.data_;
    return *this;
}
const DiscreteFunction& DiscreteFunction::operator +=(Precision rhs)
{
    data_+=rhs;
    return *this;
}
const DiscreteFunction& DiscreteFunction::operator -=(Precision rhs)
{
    data_-=rhs;
    return *this;
}
const DiscreteFunction& DiscreteFunction::operator *=(const DiscreteFunction rhs)
{
    ASSERT( checkSimilarity(rhs) );
    data_*=rhs.data_;
    return *this;
}
const DiscreteFunction& DiscreteFunction::operator /=(const DiscreteFunction rhs)
{
    ASSERT( checkSimilarity(rhs) );
    data_/=rhs.data_;
    return *this;
}
const DiscreteFunction& DiscreteFunction::operator *=(Precision rhs)
{
    data_+=rhs;
    return *this;
}
const DiscreteFunction& DiscreteFunction::operator /=(Precision rhs)
{
    data_/=rhs;
    return *this;
}


std::ostream& operator<<(std::ostream& stream, const DiscreteFunction& function)
{
    function.write(stream);
    return stream;
}

const DiscreteFunction operator -(
    const DiscreteFunction& lhs,
    const DiscreteFunction& rhs)
{
    DiscreteFunction result(lhs);
    result-=rhs;
    return result;
}
const DiscreteFunction operator +(
    const DiscreteFunction& lhs,
    const DiscreteFunction& rhs)
{
    DiscreteFunction result(lhs);
    result+=rhs;
    return result;
}
const DiscreteFunction operator *(
    const DiscreteFunction& lhs,
    const DiscreteFunction& rhs)
{
    DiscreteFunction result(lhs);
    result*=rhs;
    return result;
}
const DiscreteFunction operator /(
    const DiscreteFunction& lhs,
    const DiscreteFunction& rhs)
{
    DiscreteFunction result(lhs);
    result/=rhs;
    return result;
}

const DiscreteFunction operator +(
    const Precision lhs,
    const DiscreteFunction& rhs)
{
    DiscreteFunction result(rhs);
    result+=lhs;
    return result;
}
const DiscreteFunction operator +(
    const DiscreteFunction& rhs,
    const Precision lhs)
{
    return lhs+rhs;
}

const DiscreteFunction operator -(
    const Precision lhs,
    const DiscreteFunction& rhs)
{
    DiscreteFunction result(rhs);
    result-=lhs;
    return result;
}
const DiscreteFunction operator -(
    const DiscreteFunction& rhs,
    const Precision lhs)
{
    return lhs-rhs;
}

const DiscreteFunction operator *(
    const Precision lhs,
    const DiscreteFunction& rhs)
{
    DiscreteFunction result(rhs);
    result*=lhs;
    return result;
}
const DiscreteFunction operator *(
    const DiscreteFunction& rhs,
    const Precision lhs)
{
    return lhs*rhs;
}
const DiscreteFunction operator /(
    const Precision lhs,
    const DiscreteFunction& rhs)
{
    DiscreteFunction result(rhs);
    result/=lhs;
    return result;
}


}
