#ifndef STENCILBOARD_H_
#define STENCILBOARD_H_

#include "../general/parameters.h"
#include <vector>


namespace mg
{

class StencilBoard
{
public:
	StencilBoard();
	virtual ~StencilBoard();
    void insert( const Index level,
                 const Position pos,
                 const Index sx,
                 const Index sy,
                 const Index nx,
                 const Index ny,
                 const NumericArray& opL );
    void insert( const Index level,
                 const Position pos,
                 const PositionArray& jX,
                 const PositionArray& jY );
    bool find( const Index level,
               const Position pos,
               const Index sx,
               const Index sy,
               const Index nx,
               const Index ny ) const;
    bool find( const Index level,
               const Position pos ) const;
    const NumericArray& getL( const Index level,
                              const Position pos,
                              const Index sx,
                              const Index sy,
                              const Index nx,
                              const Index ny ) const;
    const PositionArray& getJx( const Index level,
                                const Position pos ) const;
    const PositionArray& getJy( const Index level,
                                const Position pos ) const;
private:
    std::vector<std::vector<NumericArray> >     m_coefficients;
    std::vector<std::vector<NumericArray> >     m_coefficients_b[4];
    std::vector<std::vector<NumericArray> >     m_coefficients_c;
    std::vector<std::vector<PositionArray> >    m_posJx;
    std::vector<std::vector<PositionArray> >    m_posJy;
};

}

#endif /*STENCILBOARD_H_*/
