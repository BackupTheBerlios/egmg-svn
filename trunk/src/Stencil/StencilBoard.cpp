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
                           const Position pos,
                           const Index sx,
                           const Index sy,
                           const Index nx,
                           const Index ny,
                           const NumericArray& opL )
{
    if ( pos == C )
    {
        //ASSERT( level <= m_coefficients.size() )
        if ( level - 1 == m_coefficients.size() )
        {
            m_coefficients.push_back( 
                std::vector<NumericArray>( (nx+1)*(ny+1) ) );
        }
        m_coefficients[ level - 1 ][ sy*(nx+1)+sx ].resize( opL.size() );
        m_coefficients[ level - 1 ][ sy*(nx+1)+sx ] = opL;
    }
    else if ( pos == W || pos == N || pos == E || pos == S )
    {
        //ASSERT( level <= m_coefficients_b[ pos - 1 ].size() )
        if ( level - 1 == m_coefficients_b[ pos - 1 ].size() )
        {
            m_coefficients_b[ pos - 1 ].push_back(
                std::vector<NumericArray>( std::max( nx+1, ny+1 ) ) );
        }
        m_coefficients_b[ pos-1 ][ level - 1 ][ pos%2 == 0 ? sx : sy ].resize(
            opL.size() );
        m_coefficients_b[ pos-1 ][ level - 1 ][ pos%2 == 0 ? sx : sy ] = opL;
    }
    else
    {
        //ASSERT( level <= m_coefficients_c.size() )
        if ( level - 1 == m_coefficients_c.size() )
        {
            m_coefficients_c.push_back( std::vector<NumericArray>( 4 ) );
        }
        m_coefficients_c[ level - 1 ][ pos - 5 ].resize( opL.size() );
        m_coefficients_c[ level - 1 ][ pos - 5 ] = opL;
    }
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
                         const Position pos,
                         const Index sx,
                         const Index sy,
                         const Index nx,
                         const Index ny ) const
{
    if ( pos == C )
    {
        if ( m_coefficients.size() < level )
            return false;
        if ( m_coefficients[ level - 1 ].size() != (nx+1)*(ny+1) )
            return false;
        if ( m_coefficients[ level - 1 ][ sy*(nx+1)+sx ].size() == 0 )
            return false;
    }
    else if ( pos == W || pos == N || pos == E || pos == S )
    {
        if ( m_coefficients_b[ pos - 1 ].size() < level )
            return false;
        if ( m_coefficients_b[ pos - 1 ][ level - 1 ].size() != 
                std::max( (nx+1), (ny+1) ) )
            return false;
        if ( m_coefficients_b[ pos - 1 ]
             [ level - 1 ][ pos%2 == 0 ? sx : sy ].size() == 0 )
            return false;
    }
    else
    {
        if ( m_coefficients_c.size() < level )
            return false;
        if ( m_coefficients_c[ level - 1 ].size() != 4 )
            return false;
        if ( m_coefficients_c[ level - 1 ][ pos - 5 ].size() == 0 )
            return false;
    }
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
    const Position pos,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ) const
{
    if ( pos == C )
    {
        //ASSERT( level - 1 < m_coefficients.size() )
        //ASSERT( m_coefficients[ level - 1 ].size() == (nx+1)*(ny+1) )
        //ASSERT( m_coefficients[ level - 1 ][ sy*(nx+1)+sx ].size() > 0 )
        return m_coefficients[ level - 1 ][ sy*(nx+1)+sx ];
    }
    else if ( pos == W || pos == N || pos == E || pos == S )
    {
        //ASSERT( level - 1 < m_coefficients_b[ pos - 1 ].size() )
        //ASSERT( m_coefficients_b[ pos - 1 ][ level - 1 ].size() == 
        //        std::max( nx+1, ny+1 ) )
        //ASSERT( m_coefficients_b[ pos - 1 ][ level - 1 ]
        //        [ pos%2 == 0 ? sx : sy ].size() > 0 )
        return m_coefficients_b[ pos - 1 ][ level - 1 ][ pos%2 == 0 ? sx : sy ];
    }
    else
    {
        //ASSERT( level - 1 < m_coefficients_c.size() )
        //ASSERT( m_coefficients_c[ level - 1 ].size() == 4 )
        //ASSERT( m_coefficients_c[ level - 1 ][ pos - 5 ].size() > 0 )
        return m_coefficients_c[ level - 1 ][ pos - 5 ];
    }
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
