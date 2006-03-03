/** \file ZebraLine.h
 * \author Andre Oeckerath
 * \brief ZebraLine.h contains the interface of the class ZebraLine.
 * \see LineRelaxation.h
 */
#ifndef ZEBRALINE_H_
#define ZEBRALINE_H_

#include "LineRelaxation.h"
#include "GSRedBlack.h"

namespace mg
{
/**
 * \brief ZebraLine is a class for a zebra line relaxation
 */
class ZebraLine : public mg::LineRelaxation
{
private:
    const Precision omega_;
    const GSRedBlack gsRedBlack_;
    void ninepointxzebra(
        NumericArray &u,
        const NumericArray &f, 
        NumericArray resid,
        const Stencil &stencil,
        const size_t nx, 
        const size_t ny) const;
    void ninepointyzebra(
        NumericArray &u,
        const NumericArray &f, 
        NumericArray resid,
        const Stencil &stencil,
        const size_t nx, 
        const size_t ny) const;
    void xzebra(
        NumericArray &u,
        const NumericArray &f, 
        NumericArray resid,
        const Stencil &stencil,
        const size_t nx, 
        const size_t ny) const;
    void yzebra(
        NumericArray &u,
        const NumericArray &f, 
        NumericArray resid,
        const Stencil &stencil,
        const size_t nx, 
        const size_t ny) const;
public:
    /**
     * \brief The constructor of a ZebraLine object
     * 
     * ZebraLine constructs a ZebraLine object with:
     * \param[in] preSmoothingSteps     number of pre smoothing steps  (def. 1)
     * \param[in] postSmoothingSteps    number of post smoothing steps (def. 1)
     * \param[in] direction             direction of the line relaxation
     *                                  (def. alternating directions)
     * \param[in] omega                 relaxation parameter (def. 1.0)
     * \see Direction
     */ 
    ZebraLine(
        const int preSmoothingSteps =1,
        const int postSmoothingSteps =1,
        const Direction direction =ALTDIR,
        const Precision omega =1.0)
        : LineRelaxation(preSmoothingSteps,postSmoothingSteps,direction),
          omega_(omega),gsRedBlack_() {}
    virtual ~ZebraLine() {}
    
	/**
      * \brief relax() executes one relaxation step on the input vector
      * 
      * relax() exectues one zebra line relaxation step on the input
	  * vector on a rectangular 2D gird with lexicographic ordering and the
      * discretazation Stencil for a pde
      * 
     * \param u     the vector representation of the 2D grid to perform the
     *              relaxation on this vector will be changed
     * \param f     the right hand side of the pde
     * \param nx    number of steps in x direction
     * \param ny    number of steps in y direction
     */
    void relax(
        NumericArray &u,
        const NumericArray &f, 
        const Stencil &stencil,
        const size_t nx,
        const size_t ny) const;
};
}
#endif /* ZEBRALINE_H_ */
