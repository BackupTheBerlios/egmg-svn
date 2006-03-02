/** \file max_residuum.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface of the function max_residuum
 */
#ifndef MAX_RESIDUUM_H_
#define MAX_RESIDUUM_H_


#include "../general/parameters.h"
#include "../Stencil/Stencil.h"

namespace mg
{
    
    /**
     * \brief max_residuum calculates the maximum norm of the residuum
     * 
     * max_residuum calcutates the maximum norm of the residuum (defect) of the
     * discrete PDE given by stencil, i.e
     * \f[
     * \parallel d_h\parallel_\infty =\parallel f_h-L_h\cdot u_h\parallel_\infty
     * \f]
     * 
     * \param[in] u         the acctual approximation of the solution
     * \param[in] fv        the right hand side of the pde
     * \param[in] stencil   the stencil rep. of the pde
     * \param[in] nx        number of steps in x direction
     * \param[in] ny        number of steps in y direction
     * \return              the maximum norm of the residuum
     */
    Precision max_residuum(
        const NumericArray& u,
        const NumericArray& fv,
        const Stencil& stencil,
        const size_t nx,
        const size_t ny);
}

#endif /*MAX_RESIDUUM_H_*/
