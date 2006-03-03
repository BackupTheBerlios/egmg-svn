/** \file BilinearInterpolation.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the class BilinearInterpolation
 */
#include "BilinearInterpolation.h"

namespace mg
{
    
NumericArray BilinearInterpolation::prolongate(
    const NumericArray& u,
    const Stencil&,
    const size_t nx,
    const size_t ny) const
{
    const size_t nxNew = 2*nx;
    const size_t nyNew = 2*ny;
    NumericArray result(0.0,(nxNew+1)*(nxNew+1));

    //"interpolation" of coarse grid points
    for (size_t sy=0; sy<=ny; ++sy)
        for (size_t sx=0; sx<=nx; ++sx)
            result[2*sy*(nxNew+1)+2*sx]=u[sy*(nx+1)+sx];

    //interpolation of fine grid points on coarse grid lines
    for (size_t sy=0; sy<=nyNew; sy+=2)
        for (size_t sx=1; sx<=nxNew; sx+=2)
            result[sy*(nxNew+1)+sx]=
                1./2*(result[sy*(nxNew+1)+sx-1]
                     +result[sy*(nxNew+1)+sx+1]);
    //interpolation of fine grid points on fine grid lines
    for (size_t sy=1; sy<=nyNew; sy+=2)
        for (size_t sx=0; sx<=nxNew; ++sx)
            result[sy*(nxNew+1)+sx]=
                1./2*(result[(sy-1)*(nxNew+1)+sx]
                     +result[(sy+1)*(nxNew+1)+sx]);
    return result;
}

}
