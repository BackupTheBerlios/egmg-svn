/** \file maxResiduum.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementaion of the function maxResiduum
 * \see maxResiduum.h
 */

#include <algorithm>
#include "maxResiduum.h"

namespace mg
{
Precision maxResiduum(
    const NumericArray& u,
    const NumericArray& f,
    const Stencil& stencil,
    const Index nx,
    const Index ny)
{
    Precision result=0;

    if (stencil.size()<2)
        for (Index sy=1; sy<ny; sy++)
            for (Index sx=1; sx<nx; sx++)
            {
                Precision temp_res=
                        f[sy*(nx+1)+sx]-stencil.apply(u,C,sx,sy,nx,ny);
                result=std::max(std::fabs(temp_res),result);
            }
    else
    {
        //south west corner
        Precision tempRes=f[1*(nx+1)+1]-stencil.apply(u,SW,1,1,nx,ny);
        result=std::max(std::fabs(tempRes),result);
        //south east corner
        tempRes=f[1*(nx+1)+(nx-1)]-stencil.apply(u,SE,nx-1,1,nx,ny);
        result=std::max(std::fabs(tempRes),result);
        //north west corner
        tempRes=f[(nx-1)*(nx+1)+1]-stencil.apply(u,NW,1,ny-1,nx,ny);
        result=std::max(std::fabs(tempRes),result);
        //north east corner
        tempRes=
                f[(nx-1)*(nx+1)+(nx-1)]
                -stencil.apply(u,NE,nx-1,ny-1,nx,ny);
        result=std::max(std::fabs(tempRes),result);
        //south boarder
        for (Index sx=2; sx<nx-1; sx++)
        {
            tempRes=f[1*(nx+1)+sx]-stencil.apply(u,S,sx,1,nx,ny);
            result=std::max(std::fabs(tempRes),result);
        }
        //north boarder
        for (Index sx=2; sx<nx-1; sx++)
        {
            tempRes=f[(nx-1)*(nx+1)+sx]-stencil.apply(u,N,sx,ny-1,nx,ny);
            result=std::max(std::fabs(tempRes),result);
        }
        //west boarder
        for (Index sy=2; sy<ny-1; sy++)
        {
            tempRes=f[sy*(nx+1)+1]-stencil.apply(u,W,1,sy,nx,ny);
            result=std::max(std::fabs(tempRes),result);
        }
        //east boarder
        for (Index sy=2; sy<ny-1; sy++)
        {
            tempRes=f[sy*(nx+1)+(nx-1)]-stencil.apply(u,E,nx-1,sy,nx,ny);
            result=std::max(std::fabs(tempRes),result);
        }
        //the center
        for (Index sy=2; sy<ny-1; sy++)
            for (Index sx=2; sx<nx-1; sx++)
            {
                tempRes=f[sy*(nx+1)+sx]-stencil.apply(u,C,sx,sy,nx,ny);
                result=std::max(std::fabs(tempRes),result);
            }
    }
    return result;
}
}
