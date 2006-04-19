/** \file LineGS.h
 * \author Andre Oeckerath
 * \brief LineGS.h contains the interface of the class LineGS.
 * \see Relaxation.h
 */
#ifndef LINEGS_H_
#define LINEGS_H_


#include "LineRelaxation.h"
#include "GSLexicographic.h"

namespace mg
{
/**
 * \brief LineGS is a class for a Gauss Seidel line relaxation
 */
class LineGS : public mg::LineRelaxation
{
private:
    const GSLexicographic gsLexicographic_;
    void ninepointxline(
        NumericArray &u,
        const NumericArray &f, 
        NumericArray &rhs,
        const Stencil &stencil,
        const Index nx, 
        const Index ny) const;
    void ninepointyline(
        NumericArray &u,
        const NumericArray &f, 
        NumericArray &rhs,
        const Stencil &stencil,
        const Index nx, 
        const Index ny) const;
    void xline(
        NumericArray &u,
        const NumericArray &f, 
        NumericArray &rhs,
        const Stencil &stencil,
        const Index nx, 
        const Index ny) const;
    void yline(
        NumericArray &u,
        const NumericArray &f, 
        NumericArray &rhs,
        const Stencil &stencil,
        const Index nx, 
        const Index ny) const;
public:
    /**
     * \brief The constructor of a LineGS object
     * 
     * LineGS constructs a LineGS object with:
     * \param[in] preSmoothingSteps     number of pre smoothing steps  (def. 1)
     * \param[in] postSmoothingSteps    number of post smoothing steps (def. 1)
     * \param[in] direction             direction of the line relaxation
     *                                  (def. alternating directions)
     * \see Direction
     */ 
    LineGS(
        const int preSmoothingSteps =1,
        const int postSmoothingSteps =1,
        const Direction direction =ALTDIR)
        : LineRelaxation(preSmoothingSteps,postSmoothingSteps,direction),
          gsLexicographic_() {}
    virtual ~LineGS() {}

    /**
     * \brief relax() executes one relaxation step on the input vector
     * 
     * relax() exectues one Gauss Seidel line relaxation step on the 
     * input vector on a rectangular 2D gird with lexicographic ordering and the
     * discretazation Stencil for a pde
     * 
     * \param[in,out] u     the vector representation of the 2D grid to perform
     *                      the relaxation on
     * \param[in] f         the right hand side of the pde
     * \param[in] stencil   the stencil rep. of the pde
     * \param[in] nx        number of steps in x direction
     * \param[in] ny        number of steps in y direction
     */
    void relax(
        NumericArray &u,
        const NumericArray &f, 
        const Stencil &stencil,
        const Index nx,
        const Index ny) const;
};

}

#endif /*LINEGS_H_*/
