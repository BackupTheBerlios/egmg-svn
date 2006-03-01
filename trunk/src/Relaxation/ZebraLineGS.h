/** \file ZebraLineGS.h
 * \author Andre Oeckerath
 * \brief ZebraLineGS.h contains the interface of the class ZebraLineGS.
 * \see Relaxation.h
 */
#ifndef ZEBRALINEGS_H_
#define ZEBRALINEGS_H_

#include<valarray>
#include "LineRelaxation.h"
#include "GSRedBlack.h"


namespace mg
{
/**
 * \brief ZebraLineGS is a class for a Gauss Seidel zebra line relaxation
 */
class ZebraLineGS : public mg::LineRelaxation
{
private:
    const GSRedBlack gsRedBlack_;
    void ninepointxzebra(
        std::valarray<Precision> &u,
        const std::valarray<Precision> &f, 
        std::valarray<Precision> &rhs,
        const Stencil &stencil,
        const size_t nx, 
        const size_t ny) const;
    void ninepointyzebra(
        std::valarray<Precision> &u,
        const std::valarray<Precision> &f, 
        std::valarray<Precision> &rhs,
        const Stencil &stencil,
        const size_t nx, 
        const size_t ny) const;
    void xzebra(
        std::valarray<Precision> &u,
        const std::valarray<Precision> &f, 
        std::valarray<Precision> &rhs,
        const Stencil &stencil,
        const size_t nx, 
        const size_t ny) const;
    void yzebra(
        std::valarray<Precision> &u,
        const std::valarray<Precision> &f, 
        std::valarray<Precision> &rhs,
        const Stencil &stencil,
        const size_t nx, 
        const size_t ny) const;
public:
    /**
     * \brief The constructor of a ZebraLineGS object
     * 
     * ZebraLineGS constructs a ZebraLineGS object with:
     * \param[in] preSmoothingSteps     number of pre smoothing steps  (def. 1)
     * \param[in] postSmoothingSteps    number of post smoothing steps (def. 1)
     * \param[in] direction             direction of the line relaxation
     *                                  (def. alternating directions)
     * \see Direction
     */ 
    ZebraLineGS(
        const int preSmoothingSteps =1,
        const int postSmoothingSteps =1,
        const Direction direction =ALTDIR)
        : LineRelaxation(preSmoothingSteps,postSmoothingSteps,direction),
          gsRedBlack_() {}
    virtual ~ZebraLineGS() {}

    /**
     * \brief relax() executes one relaxation step on the input vector
     * 
     * relax() exectues one Gauss Seidel zebra line relaxation step on the 
     * input vector on a rectangular 2D gird with lexicographic ordering and the
     * discretazation Stencil for a pde
     * 
     * \param u		the vector representation of the 2D grid to perform the
     * 				relaxation on this vector will be changed
     * \param f		the right hand side of the pde
     * \param nx	number of steps in x direction
     * \param ny	number of steps in y direction
     */
    void relax(
        std::valarray<Precision> &u,
        const std::valarray<Precision> &f, 
		const Stencil &stencil,
        const size_t nx,
        const size_t ny) const;	
};
}
#endif /*ZEBRALINEGS_H_*/
