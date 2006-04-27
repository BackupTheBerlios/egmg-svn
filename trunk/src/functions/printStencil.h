/** \file printStencil.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief contains the interface of the function printStencil
 */
#ifndef PRINTSTENCIL_H_
#define PRINTSTENCIL_H_

#include "../general/parameters.h"

namespace mg
{

/**
 * \brief prints the given stencil to out
 * 
 * printStencil prints a stencil given by a coeficient array and two positon
 * arrays according to the discription of Stencil. /see Stencil
 * 
 * \param[in] L     the coeficients of the stencil
 * \param[in] jX    the position array of the stencil in x dir
 * \param[in] jY    the position array of the stencil in x dir
 * \param[out] out  the ostream to write to
 */
void printStencil(
    const NumericArray& L,
    const PositionArray& jX,
    const PositionArray& jY,
    std::ostream& out);

}

#endif /*PRINTSTENCIL_H_*/
