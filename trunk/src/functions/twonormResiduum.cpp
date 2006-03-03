/** \file twonormResiduum.cpp
 * \author Matthias Rettenmeier
 * \brief contains the implimentaion of the function twonormResiduum
 * \see twonormResiduum.h
 */
#include "twonormResiduum.h"

namespace mg
{
    Precision twonormResiduum(
        const NumericArray& u,
        const NumericArray& f,
        const Stencil& stencil,
        const Index nx,
        const Index ny)
    {
        Precision result=0;
        if (stencil.size()<2)
            for (Index sy=1; sy<ny; ++sy)
                for (Index sx=1; sx<nx; ++sx)
                {
                    Precision temp_res=
                            f[sy*(nx+1)+sx]-stencil.apply(u,C,sx,sy,nx,ny);
                    result+=temp_res*temp_res;
                }
        else
        {
            //south west corner
            Precision temp_res=f[1*(nx+1)+1]-stencil.apply(u,SW,1,1,nx,ny);
            result+=temp_res*temp_res;
            //south east corner
            temp_res=f[1*(nx+1)+(nx-1)]-stencil.apply(u,SE,nx-1,1,nx,ny);
            result+=temp_res*temp_res;
            //north west corner
            temp_res=f[(nx-1)*(nx+1)+1]-stencil.apply(u,NW,1,ny-1,nx,ny);
            result+=temp_res*temp_res;
            //north east corner
            temp_res=
                    f[(nx-1)*(nx+1)+(nx-1)]
                    -stencil.apply(u,NE,nx-1,ny-1,nx,ny);
            result+=temp_res*temp_res;
            //south boarder
            for (Index sx=2; sx<nx-1; sx++)
            {
                temp_res=f[1*(nx+1)+sx]-stencil.apply(u,S,sx,1,nx,ny);
                result+=temp_res*temp_res;
            }
            //north boarder
            for (Index sx=2; sx<nx-1; sx++)
            {
                temp_res=f[(nx-1)*(nx+1)+sx]-stencil.apply(u,N,sx,ny-1,nx,ny);
                result+=temp_res*temp_res;
            }
            //west boarder
            for (Index sy=2; sy<ny-1; sy++)
            {
                temp_res=f[sy*(nx+1)+1]-stencil.apply(u,W,1,sy,nx,ny);
                result+=temp_res*temp_res;
            }
            //east boarder
            for (Index sy=2; sy<ny-1; sy++)
            {
                temp_res=f[sy*(nx+1)+(nx-1)]-stencil.apply(u,E,nx-1,sy,nx,ny);
                result+=temp_res*temp_res;
            }
            //the center
            for (Index sy=2; sy<ny-1; sy++)
                for (Index sx=2; sx<nx-1; sx++)
                {
                    temp_res=f[sy*(nx+1)+sx]-stencil.apply(u,C,sx,sy,nx,ny);
                    result+=temp_res*temp_res;
                }
        }
        result=std::sqrt(result);
        return result;
    }
}
