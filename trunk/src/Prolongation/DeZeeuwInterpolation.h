/** \file DeZeeuwInterpolation.h
 * \author Benedikt Engbroks
 * \brief Contains the interface of the class DeZeeuwInterpolation
 * 
 * This file contains the interface of DeZeeuwInterpolation. The
 * implementation is in DeZeeuwInterpolation.cpp
 */

//TODO: This file needs cleanup

#ifndef DEZEEUWINTERPOLATION_H_
#define DEZEEUWINTERPOLATION_H_

#include <valarray>
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
    std::valarray<Precision> t_;
    const std::valarray<int> jx_;
    const std::valarray<int> jy_;

    //initialize jx_, makes it possible to make jx_ const
    std::valarray<int> initJx_() const
    {
        const int t[]={0,-1,0,1,0,-1,1,1,-1};
        return std::valarray<int>(t,9);
    }

    //initialize jy_, makes it possible to make jy_ const
    std::valarray<int> initJy_() const
    {
        const int t[]={0,0,1,0,-1,1,1,-1,-1};
        return std::valarray<int>(t,9);
    }

    //we don't want the autogenerated copy constructor and assignment operator
    DeZeeuwInterpolation(const DeZeeuwInterpolation&);
    DeZeeuwInterpolation& operator=(const DeZeeuwInterpolation&);
    
    static const int SW=0;
    static const int S=1;
    static const int SE=2;
    static const int W=3;
    static const int C=4;
    static const int E=5;
    static const int NW=6;
    static const int N=7;
    static const int NE=8;
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
    std::valarray<Precision> prolongate(
        const std::valarray<Precision>& u,
        const Stencil& stencil, 
        const size_t nx,
        const size_t ny) const;

    const std::valarray<Precision>& get_I(
        const size_t sx,
        const size_t sy, 
        const size_t nx,
        const size_t ny,
        const Stencil& stencil)
    {
        std::valarray<Precision> stencilL=stencil.getL(C,0,0,nx,ny);
        std::valarray<int> jx=stencil.getJx(C);
        std::valarray<int> jy=stencil.getJy(C);
        std::valarray<size_t> position(9);
        for (size_t jj=0;jj<jx.size();jj++)
        {
            position[(jx[jj]+1)+3*(jy[jj]+1)]=jj;
        }
        std::valarray<Precision> mS(9);
        std::valarray<Precision> mT(9);

        Precision scale=0;
        Precision weight1=0;
        Precision weight2=0;
        Precision erg=0;
        Precision symsum=0;
        Precision d_w=0;
        Precision d_e=0;
        Precision d_n=0;
        Precision d_s=0;
        Precision sigma1=0;
        Precision c_1=0;
        Precision w_w=0;
        Precision w_e=0;
        Precision sigma2=0;
        Precision c_2=0;
        Precision w_s=0;
        Precision w_n=0;

        // C
        t_[0]=1.0;

        // W
        stencilL=stencil.getL(C,sx-1,sy,nx,ny);
        symsum=0;
        
        // Divide the stencil defined by stencilL und position into a symmetric 
        // and an antisymmetric part.
        for (size_t k=0; k<9; ++k)
        {
            mS[k]=0.5*(stencilL[position[k]]+stencilL[position[8-k]]);
            symsum+=mS[k];
            mT[k]=0.5*(stencilL[position[k]]-stencilL[position[8-k]]);
        }
        d_w=std::max(std::fabs(mS[SW]+mS[W]+mS[NW]),
                     std::max(std::fabs(mS[SW]),
                              std::fabs(mS[NW])));
        d_e=std::max(std::fabs(mS[SE]+mS[E]+mS[NE]),
                     std::max(std::fabs(mS[SE]),
                              std::fabs(mS[NE])));
        d_n=std::max(std::fabs(mS[NW]+mS[N]+mS[NE]),
                     std::max(std::fabs(mS[NW]),
                              std::fabs(mS[NE])));
        d_s=std::max(std::fabs(mS[SW]+mS[S]+mS[SE]),
                     std::max(std::fabs(mS[SW]),
                              std::fabs(mS[SE])));
        sigma1=0.5*std::min(1.0, std::fabs(1-symsum/stencilL[0]));
        c_1=mT[SE]+mT[E]+mT[NE]-mT[SW]-mT[W]-mT[NW];
        w_w=sigma1*(1+(d_w-d_e)/(d_w+d_e)+c_1/(d_w+d_e+d_n+d_s));
        w_e=2*sigma1-w_w;
        weight2=std::min(2*sigma1, std::max(w_e, 0.0));
        t_[1]=weight2; 

        // N
        stencilL=stencil.getL(C,sx,sy+1,nx,ny);
        symsum=0;
        
        // Divide the stencil defined by stencilL und position into a symmetric 
        // and an antisymmetric part.
        for (int k=0; k<9; ++k)
        {
            mS[k]=0.5*(stencilL[position[k]]+stencilL[position[8-k]]);
            symsum+=mS[k];
            mT[k]=0.5*(stencilL[position[k]]-stencilL[position[8-k]]);
        }
        d_w=std::max(std::fabs(mS[SW]+mS[W]+mS[NW]),
                     std::max(std::fabs(mS[SW]),
                              std::fabs(mS[NW])));
        d_e=std::max(std::fabs(mS[SE]+mS[E]+mS[NE]),
                     std::max(std::fabs(mS[SE]),
                              std::fabs(mS[NE])));
        d_n=std::max(std::fabs(mS[NW]+mS[N]+mS[NE]),
                     std::max(std::fabs(mS[NW]),
                              std::fabs(mS[NE])));
        d_s=std::max(std::fabs(mS[SW]+mS[S]+mS[SE]),
                     std::max(std::fabs(mS[SW]),
                              std::fabs(mS[SE])));
        sigma2=0.5*std::min(1.0, std::fabs(1-symsum/stencilL[0]));
        c_2=mT[NW]+mT[N]+mT[NE]-mT[SW]-mT[S]-mT[SE];
        w_n=sigma2*(1+(d_s-d_n)/(d_s+d_n)+c_2/(d_w+d_e+d_n+d_s));
        w_s=2*sigma2-w_n;
        weight1=std::min(2*sigma2, std::max(w_s, 0.0));
        weight2=std::min(2*sigma2, std::max(w_n, 0.0));
        t_[2]=weight1;

        // E
        stencilL=stencil.getL(C,sx+1,sy,nx,ny);
        symsum=0;
        
        // Divide the stencil defined by stencilL und position into a symmetric 
        // and an antisymmetric part.
        for (size_t k=0; k<9; ++k)
        {
            mS[k]=0.5*(stencilL[position[k]]+stencilL[position[8-k]]);
            symsum+=mS[k];
            mT[k]=0.5*(stencilL[position[k]]-stencilL[position[8-k]]);
        }
        d_w=std::max(std::fabs(mS[SW]+mS[W]+mS[NW]),
                     std::max(std::fabs(mS[SW]),
                              std::fabs(mS[NW])));
        d_e=std::max(std::fabs(mS[SE]+mS[E]+mS[NE]),
                     std::max(std::fabs(mS[SE]),
                              std::fabs(mS[NE])));
        d_n=std::max(std::fabs(mS[NW]+mS[N]+mS[NE]),
                     std::max(std::fabs(mS[NW]),
                              std::fabs(mS[NE])));
        d_s=std::max(std::fabs(mS[SW]+mS[S]+mS[SE]),
                     std::max(std::fabs(mS[SW]),
                              std::fabs(mS[SE])));
        sigma1=0.5*std::min(1.0, std::fabs(1-symsum/stencilL[0]));
        c_1=mT[SE]+mT[E]+mT[NE]-mT[SW]-mT[W]-mT[NW];
        w_w=sigma1*(1+(d_w-d_e)/(d_w+d_e)+c_1/(d_w+d_e+d_n+d_s));
        w_e=2*sigma1-w_w;
        weight1=std::min(2*sigma1, std::max(w_w, 0.0));
        t_[3]=weight1; 

        // S
        stencilL=stencil.getL(C,sx,sy-1,nx,ny);
        symsum=0;
        
        // Divide the stencil defined by stencilL und position into a symmetric 
        // and an antisymmetric part.
        for (int k=0; k<9; ++k)
        {
            mS[k]=0.5*(stencilL[position[k]]+stencilL[position[8-k]]);
            symsum+=mS[k];
            mT[k]=0.5*(stencilL[position[k]]-stencilL[position[8-k]]);
        }
        d_w=std::max(std::fabs(mS[SW]+mS[W]+mS[NW]),
                     std::max(std::fabs(mS[SW]),
                              std::fabs(mS[NW])));
        d_e=std::max(std::fabs(mS[SE]+mS[E]+mS[NE]),
                     std::max(std::fabs(mS[SE]),
                              std::fabs(mS[NE])));
        d_n=std::max(std::fabs(mS[NW]+mS[N]+mS[NE]),
                     std::max(std::fabs(mS[NW]),
                              std::fabs(mS[NE])));
        d_s=std::max(std::fabs(mS[SW]+mS[S]+mS[SE]),
                     std::max(std::fabs(mS[SW]),
                              std::fabs(mS[SE])));
        sigma2=0.5*std::min(1.0, std::fabs(1-symsum/stencilL[0]));
        c_2=mT[NW]+mT[N]+mT[NE]-mT[SW]-mT[S]-mT[SE];
        w_n=sigma2*(1+(d_s-d_n)/(d_s+d_n)+c_2/(d_w+d_e+d_n+d_s));
        w_s=2*sigma2-w_n;
        weight2=std::min(2*sigma2, std::max(w_n, 0.0));
        t_[4]=weight2;

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

    const std::valarray<int>& getJx() const
    {
        return jx_; 
    }

    const std::valarray<int>& getJy() const
    {
        return jy_;
    }
    
};

}

#endif /*DEZEEUWINTERPOLATION_H_*/