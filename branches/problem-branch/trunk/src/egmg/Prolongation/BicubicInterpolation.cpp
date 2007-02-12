/** \file BicubicInterpolation.cpp
 * \author Benedikt Engbroks
 * \brief Contains the implementaion of the class BicubicInterpolation
 * \see BicubicInterpolation.h
 */
#include "BicubicInterpolation.h"

namespace mg
{
    
NumericArray BicubicInterpolation::prolongate(
    const NumericArray& u,
    const Stencil&,
    const Index nx,
    const Index ny) const
{
    const Index nxNew=2*nx;
    const Index nyNew=2*ny;
    NumericArray result((nxNew+1)*(nyNew+1));
    
    //"interpolation" of coarse grid points
    for (Index sy=0; sy<=ny; ++sy)
        for (Index sx=0; sx<=nx; ++sx)
            result[2*sy*(nxNew+1)+2*sx]=u[sy*(nx+1)+sx];
        
    if (nxNew>3 && nyNew>3)
    {
        //interpolation of fine grid points on coarse grid lines
        for (Index sy=0; sy<=nyNew; sy+=2)
        {
            result[sy*(nxNew+1)+1]=
                1./8*(3*result[sy*(nxNew+1)]
                     +6*result[sy*(nxNew+1)+2]
                     -1*result[sy*(nxNew+1)+4]);

            for (Index sx=3; sx<=nxNew-2; sx+=2)
                result[sy*(nxNew+1)+sx]=
                    1./16*(-1*result[sy*(nxNew+1)+sx-3]
                           +9*result[sy*(nxNew+1)+sx-1]
                           +9*result[sy*(nxNew+1)+sx+1]
                           -1*result[sy*(nxNew+1)+sx+3]);

            result[sy*(nxNew+1)+nxNew-1]=
                1./8*(-1*result[sy*(nxNew+1)+nxNew-4]
                      +6*result[sy*(nxNew+1)+nxNew-2]
                      +3*result[sy*(nxNew+1)+nxNew]);
        }
        
        //interpolation of fine grid points on fine grid lines
        for (Index sx=0; sx<=nxNew; ++sx)
        {
            result[nxNew+1+sx]=
                1./8*(3*result[sx]+
                      6*result[2*(nxNew+1)+sx]
                     -1*result[4*(nxNew+1)+sx]); 

            for (Index j=3; j<=nyNew-2; j+=2)
                result[j*(nxNew+1)+sx]=
                    1./16*(-1*result[(j-3)*(nxNew+1)+sx]
                           +9*result[(j-1)*(nxNew+1)+sx]
                           +9*result[(j+1)*(nxNew+1)+sx]
                           -1*result[(j+3)*(nxNew+1)+sx]);

            result[(nyNew-1)*(nxNew+1)+sx]=
                1./8*(-1*result[(nyNew-4)*(nxNew+1)+sx]
                      +6*result[(nyNew-2)*(nxNew+1)+sx]
                      +3*result[nyNew*(nxNew+1)+sx]);        
        }
    }
    else
    {
        //bilinear interpolation on a grid with less than 4 fine grid points in 
        //each direction

        //interpolation of fine grid points on coarse grid lines
        for (Index sy=0; sy<=nyNew; sy+=2)
            for (Index sx=1; sx<=nxNew; sx+=2)
                result[sy*(nxNew+1)+sx]=
                    1./2*(result[sy*(nxNew+1)+sx-1]
                         +result[sy*(nxNew+1)+sx+1]);
        
        //interpolation of fine grid points on fine grid lines
        for (Index sy=1; sy<=nyNew; sy+=2)
            for (Index sx=0; sx<=nxNew; ++sx)
                result[sy*(nxNew+1)+sx]=
                    1./2*(result[(sy-1)*(nxNew+1)+sx]
                         +result[(sy+1)*(nxNew+1)+sx]);
    }
    return result;
}
}
