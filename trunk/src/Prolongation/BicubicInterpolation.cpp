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
    const size_t nx,
    const size_t ny) const
{
    const size_t nxNew=2*nx;
    const size_t nyNew=2*ny;
    NumericArray result((nxNew+1)*(nyNew+1));
    
    //"interpolation" of coarse grid points
    for (size_t j=0; j<=ny; ++j)
        for (size_t i=0; i<=nx; ++i)
            result[2*j*(nxNew+1)+2*i]=u[j*(nx+1)+i];
        
    if (nxNew>3 && nyNew>3)
    {
        //interpolation of fine grid points on coarse grid lines
        for (size_t j=0; j<=nyNew; j+=2)
        {
            result[j*(nxNew+1)+1]=
                1./8*(3*result[j*(nxNew+1)]
                     +6*result[j*(nxNew+1)+2]
                     -1*result[j*(nxNew+1)+4]);

            for (size_t i=3; i<=nxNew-2; i+=2)
                result[j*(nxNew+1)+i]=
                    1./16*(-1*result[j*(nxNew+1)+i-3]
                           +9*result[j*(nxNew+1)+i-1]
                           +9*result[j*(nxNew+1)+i+1]
                           -1*result[j*(nxNew+1)+i+3]);

            result[j*(nxNew+1)+nxNew-1]=
                1./8*(-1*result[j*(nxNew+1)+nxNew-4]
                      +6*result[j*(nxNew+1)+nxNew-2]
                      +3*result[j*(nxNew+1)+nxNew]);
        }
        
        //interpolation of fine grid points on fine grid lines
        for (size_t i=0; i<=nxNew; ++i)
        {
            result[nxNew+1+i]=
                1./8*(3*result[i]+
                      6*result[2*(nxNew+1)+i]
                     -1*result[4*(nxNew+1)+i]); 

            for (size_t j=3; j<=nyNew-2; j+=2)
                result[j*(nxNew+1)+i]=
                    1./16*(-1*result[(j-3)*(nxNew+1)+i]
                           +9*result[(j-1)*(nxNew+1)+i]
                           +9*result[(j+1)*(nxNew+1)+i]
                           -1*result[(j+3)*(nxNew+1)+i]);

            result[(nyNew-1)*(nxNew+1)+i]=
                1./8*(-1*result[(nyNew-4)*(nxNew+1)+i]
                      +6*result[(nyNew-2)*(nxNew+1)+i]
                      +3*result[nyNew*(nxNew+1)+i]);        
        }
    }
    else
    {
        //bilinear interpolation on a grid with less than 4 fine grid points in each direction

        //interpolation of fine grid points on coarse grid lines
        for (size_t j=0; j<=nyNew; j+=2)
            for (size_t i=1; i<=nxNew; i+=2)
                result[j*(nxNew+1)+i]=
                    1./2*(result[j*(nxNew+1)+i-1]
                         +result[j*(nxNew+1)+i+1]);
        
        //interpolation of fine grid points on fine grid lines
        for (size_t j=1; j<=nyNew; j+=2)
            for (size_t i=0; i<=nxNew; ++i)
                result[j*(nxNew+1)+i]=
                    1./2*(result[(j-1)*(nxNew+1)+i]
                         +result[(j+1)*(nxNew+1)+i]);
    }
    return result;
}

}
