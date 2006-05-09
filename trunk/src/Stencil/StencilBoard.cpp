#include "StencilBoard.h"

namespace mg
{

StencilBoard::StencilBoard()
{
}

StencilBoard::~StencilBoard()
{
}

void StencilBoard::insert( const Index level,
                           const Index sx,
                           const Index sy,
                           const Index nx,
                           const Index ny,
                           const NumericArray& opL )
{
    //ASSERT( level <= m_coefficients.size() )
    if ( level - 1 == m_coefficients.size() )
    {
        m_coefficients.push_back( std::vector<NumericArray>( nx*ny ) );
    }
    m_coefficients[ level - 1 ][ sy*nx+sx ].resize( opL.size() );
    m_coefficients[ level - 1 ][ sy*nx+sx ] = opL;
}
void StencilBoard::insert( const Index level,
             const Position pos,
             const PositionArray& jX,
             const PositionArray& jY )
{
    //ASSERT( level <= m_posJx.size() && m_posJx.size() == m_posJy.size() )
    //ASSERT( jX.size() == jY.size() )
    if ( level - 1 == m_posJx.size() )
    {
        m_posJx.push_back( std::vector<PositionArray>( 9 ) );
        m_posJy.push_back( std::vector<PositionArray>( 9 ) );
    }
    m_posJx[ level - 1 ][ pos ].resize( jX.size() );
    m_posJy[ level - 1 ][ pos ].resize( jY.size() );
    m_posJx[ level - 1 ][ pos ] = jX;
    m_posJy[ level - 1 ][ pos ] = jY;
}
bool StencilBoard::find( const Index level,
                         const Index sx,
                         const Index sy,
                         const Index nx,
                         const Index ny ) const
{
    if ( m_coefficients.size() < level )
        return false;
    if ( m_coefficients[ level - 1 ].size() != nx*ny )
        return false;
    if ( m_coefficients[ level - 1 ][ sy*nx+sx ].size() == 0 )
        return false;
    return true;
}
bool StencilBoard::find( const Index level,
                         const Position pos ) const
{
    if ( m_posJx.size() < level || m_posJy.size() < level  )
        return false;
    if ( m_posJx[ level - 1 ][ pos ].size() == 0 || 
         m_posJy[ level - 1 ][ pos ].size() == 0 )
        return false;
    return true;
}
const NumericArray& StencilBoard::getL(
    const Index level,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ) const
{
    //ASSERT( level - 1 < m_coefficients.size() )
    //ASSERT( m_coefficients[ level - 1 ].size() == nx*ny )
    //ASSERT( m_coefficients[ level - 1 ][ sy*nx+sx ].size() > 0 )
    return m_coefficients[ level - 1 ][ sy*nx+sx ];
}
const PositionArray& StencilBoard::getJx( const Index level,
                                          const Position pos ) const
{
    //ASSERT( level - 1 < m_posJx.size() && m_posJx.size() == m_posJy.size() )
    //ASSERT( m_posJx[ level - 1 ][ pos ].size() == 
    //        m_posJy[ level - 1 ][ pos ].size() )
    //ASSERT( m_posJx[ level - 1 ][ pos ].size() > 0 )
    return m_posJx[ level - 1 ][ pos ];
}
const PositionArray& StencilBoard::getJy( const Index level,
                                          const Position pos ) const
{
    //ASSERT( level - 1 <= m_posJy.size() && m_posJx.size() == m_posJy.size() )
    //ASSERT( m_posJx[ level - 1 ][ pos ].size() ==
    //        m_posJy[ level - 1 ][ pos ].size() )
    //ASSERT( m_posJy[ level - 1 ][ pos ].size() > 0 )
    return m_posJy[ level - 1 ][ pos ];
}

}
