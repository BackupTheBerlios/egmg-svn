/** \file DendyInterpolation.h
 * \author Benedikt Engbroks
 * \brief Contains the interface of the class DendyInterpolation
 * 
 * This file contains the interface of DendyInterpolation. The
 * implementation is in DendyInterpolation.cpp
 */

//TODO: This file needs cleanup

#ifndef DENDY_INTERPOLATION_H_
#define DENDY_INTERPOLATION_H_

#include "Prolongation.h"
#include "../general/parameters.h"
#include "../Stencil/Stencil.h"

namespace mg
{

/**
 * \brief DendyInterpolation is a matrix-dependent 2D Prolongation Operator
 * 
 * DendyInterpolation represents a 2D Prolongation Operator that uses
 * matrix-dependent interpolation by Dendy to do its job.
 */
class DendyInterpolation : public mg::Prolongation
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
    DendyInterpolation(const DendyInterpolation& rhs);
    DendyInterpolation& operator=(const DendyInterpolation& rhs);

public:                                                              
    DendyInterpolation() : t_(9), jx_(initJx_()), jy_(initJy_()) {}

    virtual ~DendyInterpolation() {}

    /**
     * \brief prolongate does a matrix-dependent interpolation on the input vector
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
                                    
    const NumericArray& getI(
        const size_t sx,
        const size_t sy, 
        const size_t nx,
        const size_t ny,
        const Stencil& stencil)
    {
        PositionArray jx=stencil.getJx(C);
        PositionArray jy=stencil.getJy(C);
        std::valarray<size_t> position(9);

        for (size_t i=0; i<jx.size(); ++i)
            position[(jx[i]+1)+3*(jy[i]+1)]=i;

        Precision scale=0;
        Precision weight1=0;
        Precision weight2=0;
        Precision erg=0;

        // C
        t_[0]=1.0;

        // W
        NumericArray stencilL=stencil.getL(C,sx-1,sy,nx,ny);
        scale=-stencilL[0];
        weight2=0;
        if (position[S]!=0)
            scale-=stencilL[position[S]];
        if (position[N]!=0)
            scale-=stencilL[position[N]];
        if (position[E]!=0)
            weight2+=stencilL[position[E]];
        if (position[SE]!=0)
            weight2+=stencilL[position[SE]];
        if (position[NE]!=0)
            weight2+=stencilL[position[NE]];
        t_[1]=weight2/scale; 

        // N
        stencilL=stencil.getL(C,sx,sy+1,nx,ny);
        scale=-stencilL[0];
        weight1=0;
        if (position[W]!=0)
            scale-=stencilL[position[W]];
        if (position[E]!=0)
            scale-=stencilL[position[E]];
        if (position[S]!=0)
            weight1+=stencilL[position[S]];
        if (position[SW]!=0)
            weight1+=stencilL[position[SW]];
        if (position[SE]!=0)
            weight1+=stencilL[position[SE]];
        t_[2]=weight1/scale;

        // E
        stencilL=stencil.getL(C,sx+1,sy,nx,ny);
        scale=-stencilL[0];
        weight1=0;
        if (position[S]!=0)
            scale-=stencilL[position[S]];
        if (position[N]!=0)
            scale-=stencilL[position[N]];
        if (position[W]!=0)
            weight1+=stencilL[position[W]];
        if (position[SW]!=0)
            weight1+=stencilL[position[SW]];
        if (position[NW]!=0)
            weight1+=stencilL[position[NW]];
        t_[3]=weight1/scale;      

        // S
        stencilL=stencil.getL(C,sx,sy-1,nx,ny);
        scale=-stencilL[0];
        weight2=0;
        if (position[W]!=0)
            scale-=stencilL[position[W]];
        if (position[E]!=0)
            scale-=stencilL[position[E]];
        if (position[N]!=0)
            weight2+=stencilL[position[N]];
        if (position[NW]!=0)
            weight2+=stencilL[position[NW]];
        if (position[NE]!=0)
            weight2+=stencilL[position[NE]];
        t_[4]=weight2/scale;

        // NW
        stencilL=stencil.getL(C,sx-1,sy+1,nx,ny);
        scale=-stencilL[0];
        erg=0;
        if (position[E]!=0)
            erg+=stencilL[position[E]]*t_[2];
        if (position[S]!=0)
            erg+=stencilL[position[S]]*t_[1];
        if (position[SE]!=0)
            erg+=stencilL[position[SE]];
        t_[5]=erg/scale;

        // NE
        stencilL=stencil.getL(C,sx+1,sy+1,nx,ny);
        scale=-stencilL[0];
        erg=0;
        if (position[W]!=0)
            erg+=stencilL[position[W]]*t_[2];
        if (position[S]!=0)
            erg+=stencilL[position[S]]*t_[3];
        if (position[SW]!=0)
            erg+=stencilL[position[SW]];
        t_[6]=erg/scale;

        // SE
        stencilL=stencil.getL(C,sx+1,sy-1,nx,ny);
        scale=-stencilL[0];
        erg=0;
        if (position[W]!=0)
            erg+=stencilL[position[W]]*t_[4];
        if (position[N]!=0)
            erg+=stencilL[position[N]]*t_[3];
        if (position[NW]!=0)
            erg+=stencilL[position[NW]];
        t_[7]=erg/scale;

        // SW
        stencilL=stencil.getL(C,sx-1,sy-1,nx,ny);
        scale=-stencilL[0];
        erg=0;
        if (position[E]!=0)
            erg+=stencilL[position[E]]*t_[4];
        if (position[N]!=0)
            erg+=stencilL[position[N]]*t_[1];
        if (position[NE]!=0)
            erg+=stencilL[position[NE]];
        t_[8]=erg/scale;      

        return t_;          
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

#endif /*DENDY_INTERPOLATION_H_*/
