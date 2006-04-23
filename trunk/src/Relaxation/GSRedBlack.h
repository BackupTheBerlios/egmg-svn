/** \file GSRedBlack.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief GSRedBlack.h contains the interface of the class GSRedBlack.
 * \see Relaxation.h
 */
#ifndef GSREDBLACK_H_
#define GSREDBLACK_H_


#include "Relaxation.h"

namespace mg
{

/**
 * \brief GSRedBlack is a class for a Gauss Seidel relaxation
 * 
 * GSRedBlack represents a Gauss Seidel relaxation with red-black ordering.
 */
class GSRedBlack : public mg::Relaxation
{
public:
    /**
     * \brief The constructor of a GSRedBlack object
     * 
     * GSRedBlack constructs a GSRedBlack object with:
     * \param[in] preSmoothingSteps     number of pre smoothing steps  (def. 1)
     * \param[in] postSmoothingSteps    number of post smoothing steps (def. 1)
     */
    GSRedBlack(const int preSmoothingSteps =1,const int postSmoothingSteps =1)
        : Relaxation(preSmoothingSteps,postSmoothingSteps) {}
        
    virtual ~GSRedBlack() {}
    
    /**
     * \brief relax() executes one relaxation step on the input vector
     * 
     * relax() exectues one Gaus seidel relaxation step on the input vector
     * one a rectangular 2D gird with red-black ordering and the
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
        const Index nx,
        const Index ny) const;
};

}

#endif /*GSREDBLACK_H_*/
