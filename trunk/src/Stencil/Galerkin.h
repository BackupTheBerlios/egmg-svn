/** \file Galerkin.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface the class Galerkin.
 */
#ifndef GALERKIN_H_
#define GALERKIN_H_

#include <map>
#include <vector>

#include "Stencil.h"
#include "../general/parameters.h"
#include "../Prolongation/Prolongation.h"
#include "../Restriction/Restriction.h"

namespace mg
{

namespace
{

class GalerkinImpl
{
private:
	class Quadruple
    {
    private:
        const Index v1_, v2_, v3_, v4_;
    public:
        Quadruple(const Index v1,const Index v2,const Index v3,const Index v4)
            :v1_(v1),v2_(v2),v3_(v3),v4_(v4) {}

        bool operator==(const Quadruple& rhs) const
        {
            if (this == &rhs)
                return true;
            return v1_==rhs.v1_ && v2_==rhs.v2_ && v3_== rhs.v3_ && v4_==rhs.v4_;
        }

        bool operator!=(const Quadruple& rhs) const
        {
            if (this == &rhs)
                return false;
            return !(*this==rhs);
        }
        bool operator<(const Quadruple& rhs) const
        {
            return 
                   v1_<rhs.v1_ 
                || v1_==rhs.v1_ && v2_<rhs.v2_
                || v1_==rhs.v1_ && v2_==rhs.v2_ && v3_<rhs.v3_
                || v1_==rhs.v1_ && v2_==rhs.v2_ && v3_==rhs.v3_ && v4_<rhs.v4_;
        }
		bool operator>(const Quadruple& rhs) const
		{
			return (*this!=rhs) && !(*this<rhs);
		}
		bool operator>=(const Quadruple& rhs) const
        {
            return !(*this<rhs);
        }
		bool operator<=(const Quadruple& rhs) const
		{
			return !(*this>rhs);
		}
    };
	std::map<Quadruple,NumericArray > operatorValues_;
	const std::vector<PositionArray > jX_;
	const std::vector<PositionArray > jY_;
	const Restriction* const restriction_;
	const Stencil* const stencil_;
	const Prolongation* const prolongation_;
	const Index size_;
	initJX_(
		const Restriction* const restriction,
		const Stencil* const stencil,
		const Prolongation* const prolongation);
	initJY_(
		const Restriction* const restriction,
		const Stencil* const stencil,
		const Prolongation* const prolongation);
	initJX_(const Stencil* const stencil);
	initJy_(const Stencil* const stencil);
public:
	GalerkinImpl(
		const Restriction* const restriction,
		const Stencil* const stencil,
		const Prolongation* const prolongation)
		: jX_(initJX_(restriction,stencil,prolongation)),
		  jY_(initJY_((restriction,stencil,prolongation)),
		  restriction_(restriction),
		  stencil_(stencil),
		  prolongation_(prolongation),
		  size_(initSize_(restriction,stencil,prolongation))
	{}
	GalerkinImpl(const Stencil* const stencil)
		: jX_(initJX_(stencil)),jY_(initJY_(stencil)),
		  restriction_(0),
		  stencil_(stencil),
		  prolongation_(0),
		  size_(stencil.size())
	{}
}

}

class Galerkin : public Stencil
{
private:
    

    std::vector< > data_;
	std::vector<std::vector<PositionArray > > jX_;
	std::vector<std::vector<PositionArray > > jY_;
    const Stencil& stencil_;
    Index size_;
    std::vector<const Prolongation&> prolongations_;
    std::vector<const Restriction&> restrictions_;
	Index currentDepth_;
    std::vector<std::vector<PositionArray > > initJx_(const Stencil&);
    std::vector<std::vector<PositionArray > > initJy_(const Stencil&);
    void updateSize_();
    void updateJxJy_();
	NumericArray computeL(
        const Position position,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny) const;
public:
    Galerkin(const Stencil& fineGridOperator)
        : data_(9),
		  jX_(initJx_(fineGridOperator)),jY_(initJy_(fineGridOperator)),
          stencil_(fineGridOperator),
		  size_(fineGridOperator.size()),
		  prolongations_(0),
		  restrictions_(0),
		  currentDepth_(prolongations_.size()){}
    virtual ~Galerkin() {}

    virtual Precision apply(
        const NumericArray& u,
        const Position position,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny) const;

    virtual Precision getCenter(
		const Position position,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny) const;

    virtual const NumericArray& getL(
        const Position position,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny) const;

    inline const PositionArray& getJx(const Position p) const
    {
		return jX_[currentDepth_][p];
    }

    inline const PositionArray& getJy(const Position p) const
    {
        return jY_[currentDepth_][p];
    }

    void pushProlongation(const Prolongation& prolongation)
    {
        prolongations_.push_back(prolongation);
		currentDepth_=prolongations_.size();
        updateJxJy_();
		updateSize_();
    }

    void popProlongation()
    {
        prolongations_.pop_back();
		currentDepth_=prolongations_.size();
        updateJxJy_();
		updateSize_();
    }

    void pushRestriction(const Restriction& restriction)
    {
        restrictions_.push_back(restriction);
		currentDepth_=restrictions_.size();
        updateJxJy_();
		updateSize_();
    }

    void popRestriction()
    {
        restrictions_.pop_back();
		currentDepth_=restrictions_.size();
        updateJxJy_();
		updateSize_();
    }

    inline Index size() const
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
