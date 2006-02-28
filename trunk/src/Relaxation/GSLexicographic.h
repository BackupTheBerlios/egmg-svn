/** \file GSLexicographic.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief GSLexicographic.h contains the interface of the class GSLexicographic.
 * \see Relaxation.h
 */
#ifndef GSLEXICOGRAPHIC_H_
#define GSLEXICOGRAPHIC_H_

#include<valarray>
#include "Relaxation.h"

namespace mg
{

/**
 * \brief GSLexicographic is a class for a Gauss Seidel relaxation
 * 
 * GSLexicographic represents a Gauss Seidel relaxation with the discretazation
 * given by Stencil and lexicographic ordering.
 */
class GSLexicographic : public mg::Relaxation
{
public:
	/**
	 * \brief The constructor of a GSLexicographic object
	 * 
	 * GSLexicographic constructs a GSLexicographic object with:
	 * \param[in] preSmoothingSteps     number of pre smoothing steps  (def. 1)
     * \param[in] postSmoothingSteps    number of post smoothing steps (def. 1)
	 */
	GSLexicographic(
        const int preSmoothingSteps =1,
        const int postSmoothingSteps =1)
        : Relaxation(preSmoothingSteps,postSmoothingSteps) {}
	virtual ~GSLexicographic() {}
    
	/**
	 * \brief relax() executes one relaxation step on the input vector
	 * 
	 * relax() exectues one Gaus seidel relaxation step on the input vector
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
        std::valarray<Precision>& u,
        const std::valarray<Precision>& f,
        const Stencil& stencil,
        const size_t nx,
        const size_t ny) const;
};

}
#endif /*GSLEXICOGRAPHIC_H_*/
