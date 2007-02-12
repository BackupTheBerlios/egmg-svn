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
    void ninePointX(
        NumericArray &u,
        const NumericArray &f, 
        const Stencil &stencil,
        const Index nx, 
        const Index ny) const;
    void ninePointY(
        NumericArray &u,
        const NumericArray &f, 
        const Stencil &stencil,
        const Index nx, 
        const Index ny) const;
    void fullX(
        NumericArray &u,
        const NumericArray &f, 
        const Stencil &stencil,
        const Index nx, 
        const Index ny) const;
    void fullY(
        NumericArray &u,
        const NumericArray &f, 
        const Stencil &stencil,
        const Index nx, 
        const Index ny) const;
public:
    /**
     * \brief The constructor of a ZebraLine object
     * 
     * ZebraLine constructs a ZebraLine object with:
     * \param[in] direction             direction of the line relaxation
     *                                  (def. alternating directions)
     * \param[in] omega                 relaxation parameter (def. 1.0)
     * \see Direction
     */ 
    ZebraLine(
        const Direction direction =ALTDIR,
        const Precision omega =1.0)
        : LineRelaxation(direction),
          omega_(omega),gsRedBlack_() {}
    virtual ~ZebraLine() {}
};
}
#endif /* ZEBRALINE_H_ */
