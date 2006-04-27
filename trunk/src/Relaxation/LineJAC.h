/** \file LineJAC.h
 * \author Andre Oeckerath
 * \brief LineJAC.h contains the interface of the class LineJAC.
 * \see LineRelaxation.h
 */
#ifndef LINEJAC_H_
#define LINEJAC_H_

#include "LineRelaxation.h"
#include "Jacobi.h"

namespace mg
{
/**
 * \brief LineJAC is a class for a Jacobi line relaxation
 */
class LineJAC : public mg::LineRelaxation
{
private:
    const Precision omega_;
    const Jacobi jacobi_;
    void ninepointxline(
        NumericArray &u,
        const NumericArray &f, 
        NumericArray resid,
        const Stencil &stencil,
        const Index nx, 
        const Index ny) const;
    void ninepointyline(
        NumericArray &u,
        const NumericArray &f, 
        NumericArray resid,
        const Stencil &stencil,
        const Index nx, 
        const Index ny) const;
    void xline(
        NumericArray &u,
        const NumericArray &f, 
        NumericArray resid,
        const Stencil &stencil,
        const Index nx, 
        const Index ny) const;
    void yline(
        NumericArray &u,
        const NumericArray &f, 
        NumericArray resid,
        const Stencil &stencil,
        const Index nx, 
        const Index Ny) const;
public:
    /**
     * \brief The constructor of a LineJAC object
     * 
     * LineJAC constructs a LineJAC object with:
     * \param[in] direction             direction of the line relaxation
     *                                  (def. alternating directions)
     * \param[in] omega                 relaxation parameter (def. 1.0)
     * \see Direction
     */ 
    LineJAC(
        const Direction direction =ALTDIR,
        const Precision omega =1.0)
        : LineRelaxation(direction),
          omega_(omega), jacobi_(omega) {}
    virtual ~LineJAC() {}
    
    /**
     * \brief relax() executes one relaxation step on the input vector
     * 
     * relax() exectues one Jacobi line relaxation step on the 
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

#endif /*LINEJAC_H_*/
