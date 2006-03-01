/** \file ZebraLineGS.h
 * \author Andre Oeckerath
 * \brief ZebraLineGS.h contains the interface of the class ZebraLineGS.
 * \see Relaxation.h
 */
#ifndef ZEBRALINEGS_H_
#define ZEBRALINEGS_H_

#include<valarray>
#include "Relaxation.h"
#include "GSRedBlack.h"


namespace mg
{
/**
 * \brief ZebraLineGS is a class for a Gauss Seidel zebra line relaxation
 */
class ZebraLineGS : public mg::Relaxation
{
private:
    const Direction direction_;
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
    void xLRSolver(
        std::valarray<Precision>& u,
        const size_t sy,
        const size_t nx,
        std::valarray<Precision>& rhs,
        std::valarray<Precision>& ndiagL,
        std::valarray<Precision>& diagR,
        const std::valarray<Precision>& ndiagR) const;
    void yLRSolver(
        std::valarray<Precision>& u,
        const size_t sx,
        const size_t nx,
        const size_t ny,
        std::valarray<Precision>& rhs,
        std::valarray<Precision>& ndiagL,
        std::valarray<Precision>& diagR,
        const std::valarray<Precision>& ndiagR) const;
    void xLRSolver(
        std::valarray<Precision>& u,
        const size_t sy,
        const size_t nx,
        std::valarray<Precision>& rhs,
        std::valarray<Precision>& ndiagL1,
        std::valarray<Precision>& ndiagL2,
        std::valarray<Precision>& diagR,
        std::valarray<Precision>& ndiagR1,
        const std::valarray<Precision>& ndiagR2) const;
    void yLRSolver(
        std::valarray<Precision>& u,
        const size_t sx,
        const size_t nx,
        const size_t ny,
        std::valarray<Precision>& rhs,
        std::valarray<Precision>& ndiagL1,
        std::valarray<Precision>& ndiagL2,
        std::valarray<Precision>& diagR,
        std::valarray<Precision>& ndiagR1,
        const std::valarray<Precision>& ndiagR2) const;
public:
    /**
     * \brief The constructor of a ZebraLineGS object
     * 
     * ZebraLineGS constructs a ZebraLineGS object with:
     * \param[in] preSmoothingSteps     number of pre smoothing steps  (def. 1)
     * \param[in] postSmoothingSteps    number of post smoothing steps (def. 1)
     * \param[in] direction             direction of the line relaxation
     * \see Direction
     */ 
    ZebraLineGS(
        const int preSmoothingSteps =1,
        const int postSmoothingSteps =1,
        const Direction direction =ALTDIR)
        : Relaxation(preSmoothingSteps,postSmoothingSteps),
          direction_(direction), gsRedBlack_() {}
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
