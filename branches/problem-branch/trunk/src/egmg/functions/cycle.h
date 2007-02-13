/** \file cycle.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the template function cycle
 */
#ifndef CYCLE_H_
#define CYCLE_H_


#include "../general/parameters.h"
#include "../Prolongation/Prolongation.h"
#include "../Relaxation/Relaxation.h"
#include "../Restriction/Restriction.h"
#include "../CycleType/CycleType.h"
#include "../Stencil/Stencil.h"
#include "../Problem/Problem.h"

namespace mg
{

/**
 * \brief cycle() does one iterative multigridcycle
 * 
 * cycle does one iterative multigridcycle to solve a partial differential
 * equation (pde) on a discrete rectangular grid. The discretization of the
 * pde is given by stencil.\n
 * CycleType controlls the CylceType \see CycleType.
 * 
 * \param[in] cycleType         the cycle type
 * \param[in,out] u             the vector represantation of the unknown
 *                              function to solve on a rectangular grid.
 * \param[in] f                 the right hand side of the pde
 * \param[in] stencil           the stencil representaion of the discrete
 *                              pde
 * \param[in] prolongation      the prolongation to use
 * \param[in] restriction       the restriction to use
 * \param[in] relaxation        the relaxation to use
 * \param[in] nx                number of steps in x direction
 * \param[in] ny                number of steps in y direction
 * \throw std::domain_error     if nx or ny is not divedable by 2
 */
void cycle(
    CycleType& cycleType,
    Problem& problem,
    const Prolongation& prolongation,
    const Restriction& restriction,
    const Relaxation& relaxation);
}

#endif /*CYCLE_H_*/
