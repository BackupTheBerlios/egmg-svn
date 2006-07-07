/** \file StencilBoard.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the class StencilBoard.
 */
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
                           const NumericArray& opL )
{
    if ( level - 1 == m_coefficients[ pos ].size() )
    {
        m_coefficients[ pos ].push_back( 
            std::map<std::pair<Index,Index>,NumericArray>() );
    }
    m_coefficients[ pos ]
                  [ level - 1 ]
                  [ std::make_pair( sx, sy ) ].resize( opL.size() );
    m_coefficients[ pos ]
                  [ level - 1 ]
                  [ std::make_pair( sx, sy ) ] = opL;
}

void StencilBoard::insert( const Index level,
                           const Position pos,
                           const PositionArray& jX,
                           const PositionArray& jY )
{
   
    //ASSERT( level <= m_posJx.size() && m_posJx.size() == m_posJy.size() )
    //ASSERT( jX.size() == jY.size() )
    if ( level - 1 == m_posJx[ pos ].size() )
    {
        m_posJx[ pos ].push_back( PositionArray( jX.size() ) );
        m_posJy[ pos ].push_back( PositionArray( jY.size() ) );
    }
    m_posJx[ pos ][ level - 1 ].resize( jX.size() );
    m_posJx[ pos ][ level - 1 ] = jX;
    m_posJy[ pos ][ level - 1 ].resize( jY.size() );
    m_posJy[ pos ][ level - 1 ] = jY;
}
bool StencilBoard::contains( const Index level,
                             const Position pos,
                             const Index sx,
                             const Index sy ) const
{
    if ( m_coefficients[ pos ].size() < level )
        return false;
    if ( m_posJx[ pos ].size() < level || m_posJy[ pos ].size() < level  )
        return false;
        
    if ( m_coefficients[ pos ][ level - 1 ].find( std::make_pair( sx, sy ) ) ==
         m_coefficients[ pos ][ level - 1 ].end() )
        return false;
    if ( m_posJx[ pos ][ level - 1 ].size() == 0 || 
         m_posJy[ pos ][ level - 1 ].size() == 0 )
        return false;
    
    return true;
}

const NumericArray& StencilBoard::getL(
    const Index level,
    const Position pos,
    const Index sx,
    const Index sy ) const
{
    CoefficientMap::const_iterator result =
        m_coefficients[ pos ][ level -1 ].find( std::make_pair( sx, sy ) );
    return result->second;
}
const PositionArray& StencilBoard::getJx( const Index level,
                                          const Position pos  ) const
{
//    assert( level - 1 <= m_posJy[ pos ].size() &&
//            m_posJx[ pos ].size() == m_posJy[ pos ].size() );
//    assert( m_posJx[ pos ][ level - 1 ][ pos ].size() ==
//            m_posJy[ pos ][ level - 1 ].size() );
//    assert( m_posJy[ pos ][ level - 1 ][ pos ].size() > 0 );
    return m_posJx[ pos ][ level - 1 ];
}
const PositionArray& StencilBoard::getJy( const Index level,
                                          const Position pos ) const
{
//    assert( level - 1 <= m_posJy[ pos ].size() &&
//            m_posJx[ pos ].size() == m_posJy[ pos ].size() );
//    assert( m_posJx[ pos ][ level - 1 ][ pos ].size() ==
//            m_posJy[ pos ][ level - 1 ].size() );
//    assert( m_posJy[ pos ][ level - 1 ][ pos ].size() > 0 );
    return m_posJy[ pos ][ level - 1 ];
}

}
