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
    //we don't want the autogenerated copy constructor and assignment operator
    LineGS(const LineGS&);
    LineGS& operator=(const LineGS&);
public:
    /**
     * \brief The constructor of a LineGS object
     * 
     * LineGS constructs a LineGS object with:
     * \param[in] direction             direction of the line relaxation
     *                                  (def. alternating directions)
     * \see Direction
     */ 
    LineGS(const Direction direction =ALTDIR)
        : LineRelaxation(direction),
          gsLexicographic_() {}
    virtual ~LineGS() {}
};

}

#endif /*LINEGS_H_*/
