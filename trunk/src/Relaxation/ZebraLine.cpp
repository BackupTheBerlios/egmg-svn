/** \file ZebraLine.cpp
 * \author Andre Oeckerath
 * \brief ZebraLine.cpp contains the implementation of the class ZebraLine.
 * \see ZebraLine.h
 * \todo check LRSovler calls for corretness
 * \todo clean up ninepointxzebra, ninepointyzebra, xzebra and yzebra by
 *       doing redudant things in seperate functions
 */
#include <iostream>
#include "ZebraLine.h"

namespace mg
{
void ZebraLine::ninePointX(
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

        for(Index sy=3; sy<ny-1 ; sy+=2)
        {
            solveNinePointXLine(u, f, temp, stencil, nx, ny, 1, sy, W, C, E);
        }

        solveNinePointXLine(u, f, temp, stencil, nx, ny, 1, ny-1, NW, N, NE);

        u = (1 - omega_) * u + omega_ * temp;

        for(Index sy=2; sy<ny ; sy+=2) 
        {                   
            solveNinePointXLine(u, f, temp, stencil, nx, ny, 1, sy, W, C, E);

            for (Index sx=1; sx<=nx; ++sx)
                u[sy*(ny+1)+sx] =
                         (1 - omega_) * u[sy*(ny+1)+sx] 
                        + omega_ * temp[sy*(ny+1)+sx];
        }

    }       
    else
    {
        gsRedBlack_.relax(u,f,stencil,nx,ny);
    }
}

void ZebraLine::ninePointY(
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

        for(Index sx=3; sx<nx-1 ; sx+=2)
        {
            solveNinePointYLine(u, f, temp, stencil, nx, ny, sx, 1, S, C, N);
        }

        solveNinePointYLine(u, f, temp, stencil, nx, ny, nx-1, 1, SE, E, NE);

        u = (1 - omega_) * u + omega_ * temp;

        for(Index sx=2; sx<nx-1 ; sx+=2)
        {
            solveNinePointYLine(u, f, temp, stencil, nx, ny, sx, 1, S, C, N);

            for (Index sy=1; sy<=ny; ++sy)
                u[sy*(ny+1)+sx] =
                         (1 - omega_) * u[sy*(ny+1)+sx] 
                        + omega_ * temp[sy*(ny+1)+sx];
        }
    }
    else
    {
        gsRedBlack_.relax(u,f,stencil,nx,ny);
    }
}

void ZebraLine::fullX(
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

        for(Index sy=3; sy<ny-2; sy+=2)
            solveFullXLine(u, f, temp, stencil, nx, ny, 1, sy, W, C, E);
        
        solveFullXLine(u, f, temp, stencil, nx, ny, 1, ny-1, NW, N, NE);

        u = (1 - omega_) * u + omega_ * temp;

        for(Index sy=2; sy<ny-1; sy+=2)
        {
            solveFullXLine(u, f, temp, stencil, nx, ny, 1, sy, W, C, E);

            for (Index sx=1; sx<=nx; ++sx)
                u[sy*(ny+1)+sx] =
                         (1 - omega_) * u[sy*(ny+1)+sx] 
                        + omega_ * temp[sy*(ny+1)+sx];
        }
     }
     else 
     {
        for(int k=0; k<2; k++)
        {
            gsRedBlack_.relax(u,f,stencil,nx,ny);
        }
    }
}

void ZebraLine::fullY(
    NumericArray &u,
    const NumericArray &f, 
    const Stencil &stencil,
    const Index nx, 
    const Index ny) const
{
    if((ny > 4) && (nx > 4))
    {
        NumericArray temp(u);

        solveFullYLine(u, f, temp, stencil, nx, ny, 1, 1, SW, W, NW);

        for(Index sx=3; sx<nx-2; sx+=2)
            solveFullYLine(u, f, temp, stencil, nx, ny, sx, 1, S, C, N);

        solveFullYLine(u, f, temp, stencil, nx, ny, nx-1, 1, SE, E, NE);

        u = (1 - omega_) * u + omega_ * temp;

        for(Index sx=2; sx<nx-1; sx+=2)
        {
            solveFullYLine(u, f, temp, stencil, nx, ny, sx, 1, S, C, N);

            for (Index sy=1; sy<=ny; ++sy)
                u[sy*(ny+1)+sx] =
                         (1 - omega_) * u[sy*(ny+1)+sx] 
                        + omega_ * temp[sy*(ny+1)+sx];
        }
    }
    else
    {
        for(int k=0; k<2; k++)
        {
            gsRedBlack_.relax(u,f,stencil,nx,ny);
        }
    }
}
}
