/** \file LineRelaxation.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief LineRelaxation.h contains the interface of the class LineRelaxation.
 * \see Relaxation.h
 */
#ifndef LINERELAXATION_H_
#define LINERELAXATION_H_


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
        NumericArray& u,
        const Index sy,
        const Index nx,
        NumericArray& rhs,
        NumericArray& ndiagL,
        NumericArray& diagR,
        const NumericArray& ndiagR) const;
    void yLRSolver(
        NumericArray& u,
        const Index sx,
        const Index nx,
        const Index ny,
        NumericArray& rhs,
        NumericArray& ndiagL,
        NumericArray& diagR,
        const NumericArray& ndiagR) const;
    void xLRSolver(
        NumericArray& u,
        const Index sy,
        const Index nx,
        NumericArray& rhs,
        NumericArray& ndiagL1,
        NumericArray& ndiagL2,
        NumericArray& diagR,
        NumericArray& ndiagR1,
        const NumericArray& ndiagR2) const;
    void yLRSolver(
        NumericArray& u,
        const Index sx,
        const Index nx,
        const Index ny,
        NumericArray& rhs,
        NumericArray& ndiagL1,
        NumericArray& ndiagL2,
        NumericArray& diagR,
        NumericArray& ndiagR1,
        const NumericArray& ndiagR2) const;
public:
    /**
     * \brief The constructor of a LineRelaxation object
     * 
     * LineRelaxation constructs a LineRelaxation object with:
     * \param[in] direction             direction of the line relaxation
     *                                  (def. alternating directions)
     * \see Direction
     */ 
    LineRelaxation(const Direction direction=ALTDIR):direction_(direction) {}
};

}

#endif /*LINERELAXATION_H_*/
