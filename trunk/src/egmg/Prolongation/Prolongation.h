/** \file Prolongation.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface of the abstract class Prolongation
 */
#ifndef PROLONGATION_H_
#define PROLONGATION_H_


#include "../general/parameters.h"
#include "../general/TransferOperator.h"
#include "../Stencil/Stencil.h"

namespace mg
{

/**
 * \brief Prolongation is the interface of a 2D Prolongation operator
 */
class Prolongation : public TransferOperator
{
public:
    virtual ~Prolongation() {}
    
    /**
     * \brief prolongate does an interpolation on the input vector
     * 
     * \param[in] u         the vector that represents a rectangel to prolongate
     * \param[in] stencil   the stencil rep. of the pde needed for matrix dep.
     *                      Prolongations.
     * \param[in] nx        Number of steps in x direction
     * \param[in] ny        Number of steps in y direction
     * \return              a vector representing the prolongated rectangel
     */
    virtual NumericArray prolongate(
        const NumericArray& u,
        const Stencil& stencil,
        const Index nx,
        const Index ny) const =0;
};

}

#endif /*PROLONGATION_H_*/
