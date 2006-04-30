/** \file generatePositionArrays.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains interface of the function generatePositionArrays.
 */
#ifndef GENERATEPOSITIONARRAYS_H_
#define GENERATEPOSITIONARRAYS_H_

#include "../general/parameters.h"

namespace mg
{

/**
 * \brief generatePositionArrays() generate position arrays for stencils
 * 
 * generatePositionArrays() generate position arrays for stencils according
 * to the numbering as discribed in the class \see stencil. This function 
 * resizes the passed PositionArray to the correct size.
 * 
 * \param[in,out] jX        PositionArray in x-direction
 * \param[in,out] jY        PositionArray in y-direction
 * \param[in] pos           the relative position of the center element
 * \param[in] size          the size of the stencil expansion(jX,jY) will
 *                          return this (\see expansion).
 * \param[in] sizeToBorder  the number of points to the closest boarder
 * \return  the expansion of the stencil
 */
void generatePositionArrays(
    PositionArray& jX,
    PositionArray& jY,
    const Position pos,
    const Index size,
    Index sizeToBorder);
    
}

#endif /*GENERATEPOSITIONARRAYS_H_*/
