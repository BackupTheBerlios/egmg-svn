/** \file Galerkin.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface the class Galerkin.
 */
#ifndef GALERKIN_H_
#define GALERKIN_H_


#include <algorithm>
#include <cmath>
#include <iostream>

#include "Stencil.h"
#include "StencilBoard.h"
#include "../general/parameters.h"
#include "../Prolongation/Prolongation.h"
#include "../Restriction/Restriction.h"

namespace mg
{

void computeGalerkin(
    NumericArray& resultL,
    PositionArray& resultJx,
    PositionArray& resultJy,
    const Position pos,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny,
    const Restriction& restriction,
    const Stencil& stencil,
    const Prolongation& prolongation);

class Galerkin : public Stencil
{
private:
    mutable StencilBoard data_;
    const   Stencil& fineGridOperator_;
    Index size_;
    std::vector<const Prolongation*> prolongations_;
    std::vector<const Restriction*> restrictions_;
    Index currentDepth_;
    void update( const Index nx, const Index ny );
    //we don't want the autogenerated copy constructor and assignment operator
    Galerkin(const Galerkin&);
    Galerkin& operator=(const Galerkin&);
public:
    Galerkin( const Stencil& fineGridOperator )
        : fineGridOperator_(fineGridOperator),
          size_(fineGridOperator.size()),
          prolongations_(0),
          restrictions_(0),
          currentDepth_(prolongations_.size()) {}
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
        const Index ny) const
    {
        return getL(position,sx,sy,nx,ny)[0];
    }

    virtual const NumericArray& getL(
        const Position position,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny) const
    {
        if ( currentDepth_ == 0 )
            return fineGridOperator_.getL( position, sx, sy, nx, ny );
        //ASSERT( data_.find( currentDepth_, sx, sy, nx, ny ) );
        return data_.getL( currentDepth_, position, sx, sy, nx, ny );
    }

    inline const PositionArray& getJx(const Position p) const
    {
        if ( currentDepth_ == 0 )
            return fineGridOperator_.getJx( p );
        return data_.getJx( currentDepth_, p );
    }

    inline const PositionArray& getJy(const Position p) const
    {
        if ( currentDepth_ == 0 )
            return fineGridOperator_.getJy( p );
        return data_.getJy( currentDepth_, p );
    }
    
    virtual void pushTransferOperators(
        const Restriction& restriction,
        const Prolongation& prolongation,
        const Index nx,
        const Index ny )
    {
        prolongations_.push_back(&prolongation);
        restrictions_.push_back(&restriction);
        update(nx,ny);
    }

    virtual void popTransferOperators()
    {
        prolongations_.pop_back();
        restrictions_.pop_back();
        currentDepth_ = prolongations_.size();
    }
    
    inline Index size() const
    {
        return size_;
    }

    inline bool isConstant() const
    {
        return fineGridOperator_.isConstant();
    }
};

}

#endif /*GALERKIN_H_*/
