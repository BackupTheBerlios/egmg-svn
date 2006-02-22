/** \file Relaxation.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface Relaxation
 */
#ifndef RELAXATION_H_
#define RELAXATION_H_

#include <valarray>
#include "../general/parameters.h"
#include "../Stencil/Stencil.h"

namespace mg
{

/**
 * \brief Relaxation is a 2D Relaxation operator
 * \todo    move implementaion of getPreSmoothingSteps, setPreSmoothingSteps,
 *          getPostSmoothingSteps and setPostSmoothingSteps into Relaxation
 * \todo    let the member function relax do the correct number of relaxations
 *          by itself
 */
class Relaxation
{
public:
    virtual ~Relaxation() {}
    
    /**
     * \brief getPreSmoothingSteps() returns the number of pre smothing steps
     * \return the number of pre somthing steps
     */
    virtual int getPreSmoothingSteps() const =0;
    /**
     * \brief setPreSmoothingSteps() sets the number of pre smothing steps
     * \param[in] preSmoothingSteps     the number of pre somthing steps
     */
    virtual void setPreSmoothingSteps(int preSmoothingSteps) =0;
    /**
     * \brief getPostSmoothingSteps() returns the number of post somthing steps
     * \return the number of post somthing steps
     */
    virtual int getPostSmoothingSteps() const =0;
    /**
     * \brief setPostSmoothingSteps() sets the number of post smothing steps
     * \param[in] postSmoothingSteps    the number of post somthing steps
     */
    virtual void setPostSmoothingSteps(int postSmoothingSteps) =0;
    
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
        std::valarray<Precision>& u,
        const std::valarray<Precision>& f,
        const Stencil& stencil,
        const size_t nx,
        const size_t ny) const=0;
};

}

#endif /*RELAXATION_H_*/
