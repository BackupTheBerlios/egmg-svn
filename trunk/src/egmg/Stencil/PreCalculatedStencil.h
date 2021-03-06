#ifndef PRECALCULATEDSTENCIL_H_
#define PRECALCULATEDSTENCIL_H_

#include "Stencil.h"
#include "StencilBoard.h"

namespace mg
{

class PreCalculatedStencil : public Stencil
{

public:
	PreCalculatedStencil( const Stencil& fineGridOperator );
	virtual ~PreCalculatedStencil();

    virtual Precision apply(
        const NumericArray& u,
        const Position position,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny ) const;

    virtual Precision getCenter(
        const Position position,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny ) const;
        
    virtual NumericArray getL(
        const Position position,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny ) const;
        
    virtual PositionArray getJx(
        const Position p,
        const Index nx,
        const Index ny ) const;

    virtual PositionArray getJy(
        const Position p,
        const Index nx,
        const Index ny ) const;
    
    virtual void pushTransferOperators(
        const Restriction&,
        const Prolongation&,
        const Index,
        const Index );

    virtual void popTransferOperators();
    
    virtual Index size() const;

    virtual bool isConstant() const;
protected:
    Index getLevel( const Index nx, const Index ny ) const;
    virtual void update(
        const Restriction& restriction,
        const Prolongation& prolongation,
        const Index nx,
        const Index ny ) =0;
    virtual void update() =0;
    mutable StencilBoard data_;
    const   Stencil& fineGridOperator_;
    Index size_;
    Index currentDepth_;
private:
    typedef std::map<std::pair<Index,Index>, Index>   LevelMap;
    LevelMap                        levelMap_;
};

}

#endif /*PRECALCULATEDSTENCIL_H_*/
