#include "PreCalculatedStencil.h"

namespace mg
{
    
Index PreCalculatedStencil::getLevel( const Index nx, const Index ny )const
{
    LevelMap::const_iterator result =
        levelMap_.find( std::make_pair( nx, ny ) );
    if ( result == levelMap_.end() )
        return 0;
    return result->second;
}

Precision PreCalculatedStencil::apply(
  const NumericArray& u,
    const Position position,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny) const
{
    Precision result=0;
    NumericArray operatorL = getL( position, sx, sy, nx, ny );
    PositionArray jX = getJx( position, nx, ny );
    PositionArray jY = getJy( position, nx, ny );
    for ( Index i = 0; i < operatorL.size(); ++i )
        result += operatorL[ i ] *
                  u[ ( sy + jY[ i ] ) * ( nx + 1 ) + sx + jX[ i ] ];
    return result;
}

Precision PreCalculatedStencil::getCenter(
    const Position position,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny) const
{
    return getL( position, sx, sy, nx, ny )[ 0 ];
}

const NumericArray& PreCalculatedStencil::getL(
    const Position position,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny) const
{
    const Index level = getLevel( nx, ny );
    if ( level == 0 )
        return fineGridOperator_.getL( position, sx, sy, nx, ny );
    //ASSERT( data_.contains( currentDepth_, sx, sy, nx, ny ) );
    return data_.getL( level, position, sx, sy );
}

const PositionArray& PreCalculatedStencil::getJx(
    const Position p,
    const Index nx,
    const Index ny ) const
{
    const Index level = getLevel( nx, ny );
    if ( level == 0 )
        return fineGridOperator_.getJx( p, nx, ny );
    return data_.getJx( level, p );
}

const PositionArray& PreCalculatedStencil::getJy(
    const Position p,
    const Index nx,
    const Index ny ) const
{
    const Index level = getLevel( nx, ny );
    if ( level == 0 )
        return fineGridOperator_.getJy( p, nx, ny );
    return data_.getJy( level, p );
}

void PreCalculatedStencil::pushTransferOperators(
        const Restriction& restriction,
        const Prolongation& prolongation,
        const Index nx,
        const Index ny )
{
    ++currentDepth_;
    LevelMap::const_iterator entry = levelMap_.find( std::make_pair( nx, ny ) );
    if ( entry == levelMap_.end() )
        levelMap_[ std::make_pair( nx, ny ) ] = currentDepth_;
    if ( currentDepth_ == levelMap_.size() )
    update( restriction, prolongation, nx, ny );
    size_ = std::max(
        std::max( std::abs( data_.getJx( currentDepth_, C ).max() ),
                  std::abs( data_.getJx( currentDepth_, C ).min() ) ),
        std::max( std::abs( data_.getJy( currentDepth_, C ).max() ),
                  std::abs( data_.getJy( currentDepth_, C ).min() ) ) );
}

void PreCalculatedStencil::popTransferOperators()
{
    --currentDepth_;
    update();
    size_ = std::max(
        std::max( std::abs( data_.getJx( currentDepth_, C ).max() ),
                  std::abs( data_.getJx( currentDepth_, C ).min() ) ),
        std::max( std::abs( data_.getJy( currentDepth_, C ).max() ),
                  std::abs( data_.getJy( currentDepth_, C ).min() ) ) );
}
        
Index PreCalculatedStencil::size() const
{
    return size_;
}

bool PreCalculatedStencil::isConstant() const
{
    return fineGridOperator_.isConstant();
}

}
