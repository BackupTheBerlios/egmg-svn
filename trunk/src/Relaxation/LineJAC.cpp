/** \file LineJAC.cpp
 * \author Andre Oeckerath
 * \brief LineJAC.cpp contains the implementation of the class LineJAC.
 * \see LineJAC.h
 * \todo check LRSovler calls for corretness
 * \todo clean up ninepointxline, ninepointyline, xline and yline by
 *       doing redudant things in seperate functions
 */
 
#include <iostream>
#include "LineJAC.h"

namespace mg
{

void LineJAC::ninePointX(
    NumericArray &u,
    const NumericArray &f, 
    const Stencil &stencil,
    const Index nx, 
    const Index ny) const
{ 
    if (nx > 2)
    {
        NumericArray temp(u);

        solveNinePointXLine(u, f, temp, stencil, nx, ny, 1, 1, SW, S, SE);

        for(Index sy=2; sy<ny-1 ; sy++)
            solveNinePointXLine(u, f, temp, stencil, nx, ny, 1, sy, W, C, E);

        solveNinePointXLine(u, f, temp, stencil, nx, ny, 1, ny-1, NW, N, NE);

        u = (1-omega_) * u + omega_ * temp;
    }
    else
    {
        jacobi_.relax(u,f,stencil,nx,ny);
    }
}

void LineJAC::ninePointY(
    NumericArray &u,
    const NumericArray &f, 
    const Stencil &stencil,
    const Index nx, 
    const Index ny) const
{ 
    if(ny > 2)
    {
        NumericArray temp(u);

        solveNinePointYLine(u, f, temp, stencil, nx, ny, 1, 1, SW, W, NW);

        for(Index sx=2; sx<nx-1 ; sx++)
            solveNinePointYLine(u, f, temp, stencil, nx, ny, sx, 1, S, C, N);

        solveNinePointYLine(u, f, temp, stencil, nx, ny, nx-1, 1, SE, E, NE);

        u = (1-omega_) * u + omega_ * temp;
    }
    else
    {
        jacobi_.relax(u,f,stencil,nx,ny);
    }
}

void LineJAC::fullX(
    NumericArray &u,
    const NumericArray &f, 
    const Stencil &stencil,
    const Index nx, 
    const Index ny) const

{
    if((ny > 4) && (nx > 4))
    {
        NumericArray temp(u);
        solveFullXLine(u, f, temp, stencil, nx, ny, 1, 1, SW, S, SE);
        
        for(Index sy=2; sy<ny-1; sy++)
            solveFullXLine(u, f, temp, stencil, nx, ny, 1, sy, W, C, E);

        solveFullXLine(u, f, temp,  stencil, nx, ny, 1, ny-1, NW, N, NE);

        u = (1-omega_) * u + omega_ * temp;
    }
    else 
    {
        for(int k=0; k<2; k++)
        {
            jacobi_.relax(u,f,stencil,nx,ny);
        }
    }
}

void LineJAC::fullY(
    NumericArray &u,
    const NumericArray &f, 
    const Stencil &stencil,
    const Index nx, 
    const Index ny) const
{
    if ((ny > 4) && (nx > 4))
    {
        NumericArray temp(u);
        solveFullYLine(u, f, temp, stencil, nx, ny, 1, 1, SW, W, NW);

        for(Index sx=2; sx<nx-1; sx++)
            solveFullYLine(u, f, temp, stencil, nx, ny, sx, 1, S, C, N);

        solveFullYLine(u, f, temp, stencil, nx, ny, nx-1, 1, SE, E, NE);

        u = (1-omega_) * u + omega_ * temp;
    }
    else
    {     
        for(int k=0; k<2; k++)
        {
            jacobi_.relax(u,f,stencil,nx,ny);
        }
    }
}
}
