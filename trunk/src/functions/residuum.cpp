/** \file residuum.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementaion of the function residuum.
 * \see residuum.h
 */
#include "residuum.h"

namespace mg
{
    NumericArray residuum(
        const NumericArray& u,
        const NumericArray& f,
        const Stencil& stencil,
        const size_t nx,
        const size_t ny)
    {
        NumericArray result(0.0,u.size());
        if (stencil.size()<2)
            for (size_t sy=1; sy<ny; sy++)
                for(size_t sx=1; sx<nx; sx++)
                    result[sy*(nx+1)+sx]=
                            f[sy*(nx+1)+sx]-stencil.apply(u,C,sx,sy,nx,ny);
        else
        {
            //south west corner
            result[1*(nx+1)+1]=
                    f[1*(nx+1)+1]-stencil.apply(u,SW,1,1,nx,ny);
            //south east corner
            result[1*(nx+1)+(nx-1)]=
                    f[1*(nx+1)+(nx-1)]-stencil.apply(u,SE,nx-1,1,nx,ny);
            //north west corner
            result[(nx-1)*(nx+1)+1]=
                    f[(nx-1)*(nx+1)+1]-stencil.apply(u,NW,1,ny-1,nx,ny);
            //north east corner
            result[(nx-1)*(nx+1)+(nx-1)]=
                    f[(nx-1)*(nx+1)+(nx-1)]
                    -stencil.apply(u,NE,nx-1,ny-1,nx,ny);
            //south boarder
            for (size_t sx=2; sx<nx-1; sx++)
                result[1*(nx+1)+sx]=
                        f[1*(nx+1)+sx]-stencil.apply(u,S,sx,1,nx,ny);
            //north boarder
            for (size_t sx=2; sx<nx-1; sx++)
                result[(nx-1)*(nx+1)+sx]=
                        f[(nx-1)*(nx+1)+sx]-stencil.apply(u,N,sx,ny-1,nx,ny);
            //west boarder
            for (size_t sy=2; sy<ny-1; sy++)
                result[sy*(nx+1)+1]=
                        f[sy*(nx+1)+1]-stencil.apply(u,W,1,sy,nx,ny);
            //east boarder;
            for (size_t sy=2; sy<ny-1; sy++)
                result[sy*(nx+1)+(nx-1)]=
                        f[sy*(nx+1)+(nx-1)]-stencil.apply(u,E,nx-1,sy,nx,ny);
            //the center
            for (size_t sy=2; sy<ny-1; sy++)
                for (size_t sx=2; sx<nx-1; sx++)
                    result[sy*(nx+1)+sx]=
                            f[sy*(nx+1)+sx]-stencil.apply(u,C,sx,sy,nx,ny);
        }
        return result;
    }
}
