/** \file Restriction.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface Restriction
 */
#ifndef RESTRICTION_H_
#define RESTRICTION_H_

#include "../Problem/Problem.h"
#include "../general/TransferOperator.h"

namespace mg
{

/**
 * \brief Restriction is a 2D restriction operator
 */
class Restriction : public TransferOperator
{
public:
    virtual ~Restriction() {}
    
    /**
     * \brief restriction() restricts the given vector to a smaller grid
     * 
     * restriction() restricts the given vector which represents a rectangular
     * grid to a smaller grid.
     * \param[in] u                 the vector to restrict
     * \param[in] stencil           the stencil rep. of the pde needed
     *                              for matrix dep. Restrictions
     * \param[in] prolongation      Prolongation used, needed for matrix dep.
     *                              Restrictions
     * \param[in] nx                number of steps in x direction
     * \param[in] ny                number of steps in y direction
     * \throw std::domain_error     if nx or ny is not divedable by 2
     * \return                      a vector with the values on the restricted
     *                              grid
     */
    virtual NumericArray restriction(
        const Problem& problem) const =0;
};

}

#endif /*RESTRICTION_H_*/
