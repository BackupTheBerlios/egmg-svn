/** \file BicubicInterpolation.h
 * \brief Contains the interface of the class BicubicInterpolation
 * \author Benedikt Engbroks
 * 
 * This file contains the interface of BicubicInterpolation. The
 * implementation is in BicubicInterpolation.cpp
 */
#ifndef BICUBICINTERPOLATION_H_
#define BICUBICINTERPOLATION_H_

#include "Prolongation.h"

namespace mg
{

/**
 * \brief BicubicInterpolation is a 2D Prolongation Operator
 * 
 * BicubicInterpolation represents a 2D Prolongation Operator that uses
 * bicubic interpolation to do its job.
 */
class BicubicInterpolation : public Prolongation
{
private:
    const NumericArray i_;
    mutable NumericArray iBorder_;
    const PositionArray jx_;
    const PositionArray jy_;

    //initialize i_, makes it possible to make i_ const
    NumericArray initI_() const
    {   // use the saved i_ only for inner points
        const Precision t[] = {
             1.0,      9.0/16,   9.0/16,   9.0/16,    9.0/16,
             81.0/256, 81.0/256, 81.0/256, 81.0/256, -1.0/16,
            -9.0/256,  1.0/256, -9.0/256, -1.0/16,   -9.0/256,
             1.0/256, -9.0/256, -1.0/16,  -9.0/256,   1.0/256,
            -9.0/256, -1.0/16,  -9.0/256,  1.0/256,  -9.0/256};
        return NumericArray(t,25);
    }

    //initialize jx_, makes it possible to make jx_ const
    PositionArray initJx_() const
    {
        const int t[] = {
             0, -1,  0,  1,  0,
            -1,  1,  1, -1, -3,
            -3, -3, -1,  0,  1,
             3,  3,  3,  3,  3,
             1,  0, -1, -3, -3};
        return PositionArray(t,25);
    }

    //initialize jy_, makes it possible to make jy_ const
    PositionArray initJy_() const
    {
        const int t[] = {
             0,  0,  1,  0, -1,
             1,  1, -1, -1,  0,
             1,  3,  3,  3,  3,
             3,  1,  0, -1, -3,
            -3, -3, -3, -3, -1};
        return PositionArray(t,25);
    }
    //we don't want the autogenerated copy constructor and assignment operator
    BicubicInterpolation(const BicubicInterpolation&);
    BicubicInterpolation& operator=(const BicubicInterpolation&);
public:
    BicubicInterpolation()
        : i_(initI_()), iBorder_(25), jx_(initJx_()), jy_(initJy_()) {}

    virtual ~BicubicInterpolation() {}

    /**
     * \brief prolongate does a bicubic interpolation on the input vector
     * 
     * \param[in] u         the vector representing a rectangel to prolongate
     * \param[in] stencil   the stencil rep. of the pde needed for matrix dep.
     *                      Prolongations (not used).
     * \param[in] nx        Number of steps in x direction
     * \param[in] ny        Number of steps in y direction
     * \return              a vector representing the prolongated rectangel of 
     *                      size 2*(nx+1)*2*(ny+1)
     */
    NumericArray prolongate(
        const NumericArray& u,
        const Stencil& stencil,
        const Index nx,
        const Index ny) const;

    const NumericArray& getI(
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny,
        const Stencil&) const
    {
        if (sy==1)        // unterer Rand
        {
            if (sx==1)   // untere linke Ecke
            {
                const Precision t[] = {
                    1.0,     3.0/8,   3.0/4,   3.0/4,  3.0/8,
                    9.0/32,  9.0/16,  9.0/32,  9.0/64, 0.0,
                    0.0,     0.0,    -3.0/64, -1.0/8, -3.0/32,
                    1.0/64, -3.0/32, -1.0/8,  -3.0/64, 0.0,
                    0.0,     0.0,     0.0,     0.0,    0.0};
                iBorder_ = NumericArray(t,25);
                return iBorder_;                
            }
            else if (sx==nx-1)   // untere rechte Ecke
            {
                const Precision t[] = {
                     1.0,    3.0/4,   3.0/4,   3.0/8,   3.0/8,
                     9.0/16, 9.0/32,  9.0/64,  9.0/32, -1.0/8,
                    -3.0/32, 1.0/64, -3.0/32, -1.0/8,  -3.0/64,
                     0.0,    0.0,     0.0,     0.0,     0.0,
                     0.0,    0.0,     0.0,     0.0,    -3.0/64};
                iBorder_ = NumericArray(t,25);
                return iBorder_;                
            }
            else
            {
                const Precision t[] = {
                     1.0,      9.0/16,   3.0/4,    9.0/16,    3.0/8,
                     27.0/64,  27.0/64,  27.0/128, 27.0/128, -1.0/16,
                    -3.0/64,   1.0/128, -9.0/128, -1.0/8,    -9.0/128,
                     1.0/128, -3.0/64,  -1.0/16,  -3.0/128,   0.0,
                     0.0,      0.0,      0.0,      0.0,      -3.0/128};
                iBorder_ = NumericArray(t,25);
                return iBorder_;                
            }
        }
        else if (sy==ny-1)  // oberer Rand
        {
            if (sx==1)   // obere linke Ecke
            {
                const Precision t[] = {
                     1.0,     3.0/8,   3.0/8,  3.0/4,  3.0/4,
                     9.0/64,  9.0/32,  9.0/16, 9.0/32, 0.0,
                     0.0,     0.0,     0.0,    0.0,    0.0,
                     0.0,    -3.0/64, -1.0/8, -3.0/32, 1.0/64,
                    -3.0/32, -1.0/8,  -3.0/64, 0.0,    0.0};
                iBorder_ = NumericArray(t,25);
                return iBorder_;                
            }
            else if (sx==nx-1)   // obere rechte Ecke
            {
                const Precision t[] = {
                     1.0,     3.0/4,  3.0/8,  3.0/8,   3.0/4,
                     9.0/32,  9.0/64, 9.0/32, 9.0/16, -1.0/8,
                    -3.0/64,  0.0,    0.0,    0.0,     0.0,
                     0.0,     0.0,    0.0,    0.0,     0.0,
                    -3.0/64, -1.0/8, -3.0/32, 1.0/64, -3.0/32};
                iBorder_ = NumericArray(t,25);
                return iBorder_;                
            }
            else
            {
                const Precision t[] = {
                     1.0,      9.0/16,   3.0/8,   9.0/16,   3.0/4,
                     27.0/128, 27.0/128, 27.0/64, 27.0/64, -1.0/16,
                    -3.0/128,  0.0,      0.0,     0.0,      0.0,
                     0.0,     -3.0/128, -1.0/16, -3.0/64,   1.0/128,
                    -9.0/128, -1.0/8,   -9.0/128, 1.0/128, -3.0/64};
                iBorder_ = NumericArray(t,25);
                return iBorder_;                
            }
        }
        else if (sx==1)   // linker Rand
        {
            const Precision t[] = {
                1.0,      3.0/8,    9.0/16,   3.0/4,    9.0/16,
                27.0/128, 27.0/64,  27.0/64,  27.0/128, 0.0,
                0.0,      0.0,     -3.0/128, -1.0/16,  -3.0/64,
                1.0/128, -9.0/128, -1.0/8,   -9.0/128,  1.0/128,
                -3.0/64, -1.0/16,  -3.0/128,  0.0,      0.0};
            iBorder_ = NumericArray(t,25);
            return iBorder_;                
        }
        else if (sx==nx-1)  // rechter Rand
        {
            const Precision t[] = {
                 1.0,      3.0/4,    9.0/16,   3.0/8,    9.0/16,
                 27.0/64,  27.0/128, 27.0/128, 27.0/64, -1.0/8,
                -9.0/128,  1.0/128, -3.0/64,  -1.0/16,  -3.0/128,
                 0.0,      0.0,      0.0,      0.0,      0.0,
                -3.0/128, -1.0/16,  -3.0/64,   1.0/128, -9.0/128};
            iBorder_ = NumericArray(t,25);
            return iBorder_;                
        }
        else
            return i_;
    }

    const PositionArray& getJx() const
    {
        return jx_; 
    }

    const PositionArray& getJy() const
    {
        return jy_;
    }
};

}

#endif /*BICUBICINTERPOLATION_H_*/
