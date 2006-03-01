/** \file LineRelaxation.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief LineRelaxation.h contains the interface of the class LineRelaxation.
 * \see Relaxation.h
 */
#ifndef LINERELAXATION_H_
#define LINERELAXATION_H_

#include<valarray>
#include "Relaxation.h"

namespace mg
{
/**
 * \brief LineRelaxation is an abstract class for generale line relaxations
 */
class LineRelaxation : public mg::Relaxation
{
protected:
    const Direction direction_;
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
     * \brief The constructor of a LineRelaxation object
     * 
     * LineRelaxation constructs a LineRelaxation object with:
     * \param[in] preSmoothingSteps     number of pre smoothing steps  (def. 1)
     * \param[in] postSmoothingSteps    number of post smoothing steps (def. 1)
     * \param[in] direction             direction of the line relaxation
     *                                  (def. alternating directions)
     * \see Direction
     */ 
    LineRelaxation(
        const int preSmoothingSteps =1,
        const int postSmoothingSteps =1,
        const Direction direction =ALTDIR)
        : Relaxation(preSmoothingSteps,postSmoothingSteps),
          direction_(direction) {}
};

}

#endif /*LINERELAXATION_H_*/
