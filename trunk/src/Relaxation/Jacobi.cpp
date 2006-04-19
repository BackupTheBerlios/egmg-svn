/** \file Jacobi.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Jacobi.cpp contains the implementaion of the class Jacobi.
 * \see Jacobi.h
 */

#include "Jacobi.h"

namespace mg
{
void Jacobi::relax(
    NumericArray& u,
    const NumericArray& f,
    const Stencil& stencil,
    const Index nx,
    const Index ny) const
{
    Precision factor = 1.0;
    NumericArray result(u);
    if (stencil.size() < 2)
    {
        for (Index sy=1;sy<ny;sy++)
            for (Index sx=1;sx<nx;sx++)
            {
                factor = 1.0/stencil.getCenter(C,sx,sy,nx,ny);
                result[sy*(nx+1)+sx]+=factor*(f[sy*(nx+1)+sx]
                        -stencil.apply(u,C,sx,sy,nx,ny));
            }
    }
    else
    {
        //corners
        factor = 1.0/stencil.getCenter(SW,1,1,nx,ny);
        result[1*(nx+1)+1]+=factor*(f[1*(nx+1)+1]
                -stencil.apply(u,SW,1,1,nx,ny));
        factor = 1.0/stencil.getCenter(SE,(nx-1),1,nx,ny);
        result[1*(nx+1)+(nx-1)]+=factor*(f[1*(nx+1)+(nx-1)]
                -stencil.apply(u,SE,(nx-1),1,nx,ny));
        factor = 1.0/stencil.getCenter(NW,1,(nx-1),nx,ny);
        result[(nx-1)*(nx+1)+1]+=factor*(f[(nx-1)*(nx+1)+1]
                -stencil.apply(u,NW,1,(ny-1),nx,ny));
        factor = 1.0/stencil.getCenter(NE,(nx-1),(nx-1),nx,ny);
        result[(nx-1)*(nx+1)+(nx-1)]+=factor*(f[(nx-1)*(nx+1)+(nx-1)]
            -stencil.apply(u,NE,(nx-1),(ny-1),nx,ny));
        //boarders
        for (Index sx=2;sx<(nx-1);sx++)
        {
            factor = 1.0/stencil.getCenter(S,sx,1,nx,ny);
            result[1*(nx+1)+sx]+=factor*(f[1*(nx+1)+sx]
                -stencil.apply(u,S,sx,1,nx,ny));
            factor = 1.0/stencil.getCenter(N,sx,(nx-1),nx,ny);
            result[(nx-1)*(nx+1)+sx]+=factor*(f[(nx-1)*(nx+1)+sx]
                -stencil.apply(u,N,sx,(ny-1),nx,ny));
        }
        for (Index sy=2;sy<(ny-1);sy++)
        {
            factor = 1.0/stencil.getCenter(W,1,sy,nx,ny);
            result[sy*(nx+1)+1]+=factor*(f[sy*(nx+1)+1]
                -stencil.apply(u,W,1,sy,nx,ny));
            factor = 1.0/stencil.getCenter(E,(nx-1),sy,nx,ny);
            result[sy*(nx+1)+(nx-1)]+=factor*(f[sy*(nx+1)+(nx-1)]
                -stencil.apply(u,E,(nx-1),sy,nx,ny));
        }
        //center
        for (Index sy=2;sy<(ny-1);sy++)
        {
            for (Index sx=2;sx<(nx-1);sx++)
            {
                factor = 1.0/stencil.getCenter(C,sx,sy,nx,ny);
                result[sy*(nx+1)+sx]+=factor*(f[sy*(nx+1)+sx]
                    -stencil.apply(u,C,sx,sy,nx,ny));
            }
        }
    }
    u=(1-omega_)*u+omega_*result;
}
}
