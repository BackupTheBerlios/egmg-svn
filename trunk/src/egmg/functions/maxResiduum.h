/** \file maxResiduum.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface of the function maxResiduum
 */
#ifndef MAXRESIDUUM_H_
#define MAXRESIDUUM_H_


#include "../general/parameters.h"
#include "../Stencil/Stencil.h"

namespace mg
{
    
/**
 * \brief maxResiduum calculates the maximum norm of the residuum
 * 
 * maxResiduum calcutates the maximum norm of the residuum (defect) of the
 * discrete PDE given by stencil, i.e
 * \f[
 * \parallel d_h\parallel_\infty =\parallel f_h-L_h\cdot u_h\parallel_\infty
 * \f]
 * 
 * \param[in] u         the acctual approximation of the solution
 * \param[in] f         the right hand side of the pde
 * \param[in] stencil   the stencil rep. of the pde
 * \param[in] nx        number of steps in x direction
 * \param[in] ny        number of steps in y direction
 * \return              the maximum norm of the residuum
 */
Precision maxResiduum(
    const NumericArray& u,
    const NumericArray& f,
    const Stencil& stencil,
    const Index nx,
    const Index ny);
}

#endif /*MAXRESIDUUM_H_*/
