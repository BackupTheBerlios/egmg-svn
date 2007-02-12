/** \file BilinearInterpolation.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface of the class BilinearInterpolation
 */
#ifndef BILINEARINTERPOLATION_H_
#define BILINEARINTERPOLATION_H_

#include "Prolongation.h"

namespace mg
{

/**
 * \brief BilinearInterpolation is a 2D Prolongation Operator
 * 
 * BilinearInterpolation represents a 2D Prolongation Operator that uses
 * bilinear interpolation to do its job.
 */
class BilinearInterpolation : public Prolongation
{
private:
    const NumericArray i_;
    const PositionArray jx_;
    const PositionArray jy_;

    //initialize i_, makes it possible to make i_ const
    NumericArray initI_() const
    {
        const Precision t[] = {
            1.0,   1.0/2, 1.0/2,
            1.0/2, 1.0/2, 1.0/4,
            1.0/4, 1.0/4, 1.0/4};
        return NumericArray(t,9);
    }

    //initialize jx_, makes it possible to make jx_ const
    PositionArray initJx_() const
    {
        const int t[] = {
            0, -1,  0,
            1,  0, -1,
            1,  1, -1};
        return PositionArray(t,9);
    }

    //initialize jy_, makes it possible to make jy_ const
    PositionArray initJy_() const
    {
        const int t[] = {
            0,  0,  1,
            0, -1,  1,
            1, -1, -1};
        return PositionArray(t,9);
    }
    //we don't want these autogenerated contrs and operators
    BilinearInterpolation(const BilinearInterpolation&);
    BilinearInterpolation& operator=(const BilinearInterpolation&);
public:
    BilinearInterpolation()
        : i_(initI_()), jx_(initJx_()), jy_(initJy_()) {}
    virtual ~BilinearInterpolation() {}
    /**
     * \brief prolongate does a bilinear interpolation on the input vector
     * 
     * \param[in] u         the vector representing a rectangel to prolongate
     * \param[in] stencil   the stencil rep. of the pde needed for matrix dep.
     *                      Prolongations (not used for Bilinear int.).
     * \param[in] nx        Number of steps in x direction
     * \param[in] ny        Number of steps in y direction
     * \return              a vector representing the prolongated rectangel of 
     *                      size 2*(nx+1)*2*(ny+1)
     */
    NumericArray prolongate(
        const NumericArray& u,
        const Stencil&,
        const Index nx,
        const Index ny) const;

    const NumericArray& getI(
        const Position,
        const Index,
        const Index,
        const Index,
        const Index,
        const Stencil&) const
    {
        return i_;  
    }

    const PositionArray& getJx( const Position ) const
    {
        return jx_; 
    }

    const PositionArray& getJy( const Position ) const
    {
        return jy_;
    }
};

}

#endif /*BILINEARINTERPOLATION_H_*/