/** \file Jacobi.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Jacobi.h contains the interface of the class Jacobi.
 * \see Relaxation.h
 */
#ifndef JACOBI_H_
#define JACOBI_H_

#include "Relaxation.h"

namespace mg
{

/**
 * \brief Jacobi is a class for a Jacobi relaxation
 * 
 * Jacobi represents a Jacobi relaxation with the discretazation
 * given by Stencil.
 */
class Jacobi : public mg::Relaxation
{
private:
    const Precision omega_;
public:
    /**
     * \brief The constructor of a Jacobi object
     * 
     * Jacobi constructs a Jacobi object with:
     * \param[in] preSmoothingSteps     number of pre smoothing steps  (def. 1)
     * \param[in] postSmoothingSteps    number of post smoothing steps (def. 1)
     * \param[in] omega                 relaxation parameter (def. 1.0)
     */
    Jacobi(
        const int preSmoothingSteps =1,
        const int postSmoothingSteps =1,
        const Precision omega =1.0)
        : Relaxation(preSmoothingSteps,postSmoothingSteps),
          omega_(omega) {}
    virtual ~Jacobi() {}
    
    /**
     * \brief relax() executes one relaxation step on the input vector
     * 
     * relax() exectues one Jacobi relaxation step on the input vector
     * one a rectangular 2D grid with lexicographic ordering and the
     * discretazation given by stencil for a pde.
     * 
     * \param[in,out] u     the vector representation of the 2D grid to perform
     *                      the relaxation on
     * \param[in] f         the right hand side of the pde
     * \param[in] stencil   the stencil rep. of the pde
     * \param[in] nx        number of steps in x direction
     * \param[in] ny        number of steps in y direction
     */
    void relax(
        NumericArray& u,
        const NumericArray& f,
        const Stencil& stencil,
        const size_t nx,
        const size_t ny) const;
};
}

#endif /*JACOBI_H_*/
