/** \file putrhs.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface of the function putrhs
 */
#ifndef PUTRHS_H_
#define PUTRHS_H_


#include "../general/parameters.h"
#include "../Function/Function.h"

namespace mg
{
    /**
     * \brief putrhs() fills the given vector with the rhs of the pde
     * 
     * \param[in,out] fv    the vector of the right hand side of the pde
     * \param[in] nx        number of steps in x direction
     * \param[in] ny        number of steps in y direction
     * \param[in] f         the right hand side function of the pde 
     */
    void putrhs(
        NumericArray& fv,
        const Index nx,
        const Index ny,
        const Function& f);
}

#endif /*PUTRHS_H_*/
