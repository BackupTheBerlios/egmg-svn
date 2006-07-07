/** \file expansion.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains function expansion.
 */
#ifndef EXPANSION_H_
#define EXPANSION_H_

#include "../general/parameters.h"
#include <algorithm>

namespace mg
{

/**
 * \brief expansion() calculates the expansion of a stencil
 * 
 * expansion calculates the expansion of a stencil with the given
 * PositionArrays according to the size of a stencil. /see Stencil
 * 
 * \param[in] jX    PositionArray in x-direction
 * \param[in] jY    PositionArray in y-direction
 * \return  the expansion of the stencil
 */
inline Index expansion(
    const PositionArray& jX,
    const PositionArray& jY)
{
	assert( jX.size() > 0 );
	assert( jY.size() > 0 );
    const Index expansionJx1 = static_cast<Index>(std::abs(jX.max()));
    const Index expansionJx2 = static_cast<Index>(std::abs(jX.min()));
    const Index expansionJy1 = static_cast<Index>(std::abs(jY.max()));
    const Index expansionJy2 = static_cast<Index>(std::abs(jY.min()));
    const Index expansionJx = std::max(expansionJx1,expansionJx2);
    const Index expansionJy = std::max(expansionJy1,expansionJy2);
    return std::max(expansionJx,expansionJy);
}

}
#endif /*EXPANSION_H_*/
