#ifndef GALERKIN_H_
#define GALERKIN_H_

#include <map>
#include <vector>
#include <valarray>
#include "Stencil.h"
#include "../general/parameters.h"
#include "../Prolongation/Prolongation.h"
#include "../Restriction/Restriction.h"

namespace mg
{

class Galerkin : public Stencil
{
private:
    class Quadruple
    {
    private:
        const size_t i_, j_, k_, l_;
    public:
        Quadruple(const size_t i,const size_t j,const size_t k,const size_t l)
            :i_(i),j_(j),k_(k),l_(l) {}

        bool operator==(const Quadruple& rhs) const
        {
            if (this == &rhs)
                return true;
            return i_==rhs.i_ && j_==rhs.j_ && k_== rhs.k_ && l_==rhs.l_;
        }

        bool operator!=(const Quadruple& rhs) const
        {
            if (this == &rhs)
                return false;
            return !(*this==rhs);
        }

        bool operator>(const Quadruple& rhs) const
        {
            return !(*this<rhs);
        }
        bool operator<(const Quadruple& rhs) const

        {
            return 
                   i_<rhs.i_ 
                || i_==rhs.i_ && j_<rhs.j_
                || i_==rhs.i_ && j_==rhs.j_ && k_<rhs.k_
                || i_==rhs.i_ && j_==rhs.j_ && k_==rhs.k_ && l_<rhs.l_;
        }

    };

    std::vector<std::map<Quadruple,std::valarray<Precision> > > data_;
    std::vector<std::valarray<int> > jx_;
    std::vector<std::valarray<int> > jy_;
    const Stencil& stencil_;
    size_t size_;
    std::vector<const Prolongation*> prolongations_;
    std::vector<const Restriction*> restrictions_;
    std::vector<std::valarray<int> > initJx_(const Stencil&);
    std::vector<std::valarray<int> > initJy_(const Stencil&);
    void updateSize_();
    void updateJxJy_();
public:
    Galerkin(const Stencil& fineGridOperator)
        : jx_(initJx_(fineGridOperator)), jy_(initJy_(fineGridOperator)),
          stencil_(fineGridOperator), size_(fineGridOperator.size()) {}
    virtual ~Galerkin() {}

    virtual Precision apply(
        const std::valarray<Precision>&,
        const Position,
        const size_t,
        const size_t,
        const size_t,
        const size_t) const;

    virtual Precision getCenter(
        const Position,
        const size_t,
        const size_t,
        const size_t,
        const size_t) const;

    virtual const std::valarray<Precision>& getL(
        const Position,
        const size_t,
        const size_t,
        const size_t,
        const size_t) const;

    inline const std::valarray<int>& getJx(const Position p) const
    {
        return jx_[p];
    }

    inline const std::valarray<int>& getJy(const Position p) const
    {
        return jy_[p];
    }

    void pushProlongation(const Prolongation& prolongation)
    {
        prolongations_.push_back(&prolongation);
        updateSize_();
        updateJxJy_();
    }

    void popProlongation()
    {
        prolongations_.pop_back();
        updateSize_();
        updateJxJy_();
    }

    void pushRestriction(const Restriction& restriction)
    {
        restrictions_.push_back(&restriction);
        updateSize_();
        updateJxJy_();
    }

    void popRestriction()
    {
        restrictions_.pop_back();
        updateSize_();
        updateJxJy_();
    }

    inline size_t size() const
    {
        return size_;
    }

    inline bool isConstant() const
    {
        return stencil_.isConstant();
    }
};

}

#endif /*GALERKIN_H_*/
