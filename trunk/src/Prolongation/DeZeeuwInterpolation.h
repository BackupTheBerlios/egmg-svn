/** \file DeZeeuwInterpolation.h
 * \author Benedikt Engbroks
 * \brief Contains the interface of the class DeZeeuwInterpolation
 * 
 * This file contains the interface of DeZeeuwInterpolation. The
 * implementation is in DeZeeuwInterpolation.cpp
 * \todo This file needs cleanup
 * \todo make this work with all type of stencils
 */
#ifndef DEZEEUWINTERPOLATION_H_
#define DEZEEUWINTERPOLATION_H_


#include "Prolongation.h"
#include "../Stencil/Stencil.h"

namespace mg
{

/**
 * \brief DeZeeuwInterpolation is a matrix-dependent 2D Prolongation Operator
 * 
 * DeZeeuwInterpolation represents a 2D Prolongation Operator that uses
 * matrix-dependent interpolation by De Zeeuw to do its job.
 */
class DeZeeuwInterpolation : public mg::Prolongation
{
private:
    NumericArray t_;
    const PositionArray jx_;
    const PositionArray jy_;

    //initialize jx_, makes it possible to make jx_ const
    PositionArray initJx_() const
    {
        const int t[]={0,-1,0,1,0,-1,1,1,-1};
        return PositionArray(t,9);
    }

    //initialize jy_, makes it possible to make jy_ const
    PositionArray initJy_() const
    {
        const int t[]={0,0,1,0,-1,1,1,-1,-1};
        return PositionArray(t,9);
    }

    //we don't want the autogenerated copy constructor and assignment operator
    DeZeeuwInterpolation(const DeZeeuwInterpolation&);
    DeZeeuwInterpolation& operator=(const DeZeeuwInterpolation&);
    
    /**
     * \todo for this the enum Position should be used
    static const int SW=0;
    static const int S=1;
    static const int SE=2;
    static const int W=3;
    static const int C=4;
    static const int E=5;
    static const int NW=6;
    static const int N=7;
    static const int NE=8;
     */
public:
    DeZeeuwInterpolation() : t_(9), jx_(initJx_()), jy_(initJy_()) {}

    virtual ~DeZeeuwInterpolation() {}

    /**
     * \brief prolongate does a matrix-dependent interpolation on the input 
     *        vector
     * 
     * \param u     the vector representing a rectangle to prolongate
     * \param nx    Number of steps in x direction
     * \param ny    Number of steps in y direction
     * \return      a vector representing the prolongated rectangle of 
     *              size 2*(nx+1)*2*(ny+1)
     */
    NumericArray prolongate(
        const NumericArray& u,
        const Stencil& stencil, 
        const size_t nx,
        const size_t ny) const;

    const NumericArray& get_I(
        const size_t sx,
        const size_t sy, 
        const size_t nx,
        const size_t ny,
        const Stencil& stencil);

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

#endif /*DEZEEUWINTERPOLATION_H_*/
