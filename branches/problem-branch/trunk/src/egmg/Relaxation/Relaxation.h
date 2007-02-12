/** \file Relaxation.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface Relaxation
 */
#ifndef RELAXATION_H_
#define RELAXATION_H_


#include "../Problem/Problem.h"

namespace mg
{

/**
 * \brief Relaxation is a 2D Relaxation operator
 * \todo Implement ILU Smoother
 */
class Relaxation
{
public:
    virtual ~Relaxation() {}
    
    /**
     * \brief relax() executes one relaxation step on the input vector
     * 
     * \param[in,out] u         the vector representation of the 2D grid to
     *                          perform the relaxation on
     * \param[in] f             the right hand side of the pde
     * \param[in] stencil       the stencil rep. of the pde
     * \param[in] nx            number of steps in x direction
     * \param[in] ny            number of steps in y direction
     */
    virtual void relax(
        Problem& problem) const=0;
};

}

#endif /*RELAXATION_H_*/
