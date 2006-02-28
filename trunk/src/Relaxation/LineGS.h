#ifndef LINEGS_H_
#define LINEGS_H_

#include "Relaxation.h"

namespace mg
{

class LineGS : public mg::Relaxation
{
private:
    Direction direction_;
    void ninepointxline(
        std::valarray<Precision> &u,
        const std::valarray<Precision> &f, 
        std::valarray<Precision> &rhs,
        const Stencil &stencil,
        const size_t nx, 
        const size_t ny) const;
    void ninepointyline(
        std::valarray<Precision> &u,
        const std::valarray<Precision> &f, 
        std::valarray<Precision> &rhs,
        const Stencil &stencil,
        const size_t nx, 
        const size_t ny) const;
    void xline(
        std::valarray<Precision> &u,
        const std::valarray<Precision> &f, 
        std::valarray<Precision> &rhs,
        const Stencil &stencil,
        const size_t nx, 
        const size_t ny) const;
    void yline(
        std::valarray<Precision> &u,
        const std::valarray<Precision> &f, 
        std::valarray<Precision> &rhs,
        const Stencil &stencil,
        const size_t nx, 
        const size_t ny) const;
public:
	LineGS();
	virtual ~LineGS();
};

}

#endif /*LINEGS_H_*/
