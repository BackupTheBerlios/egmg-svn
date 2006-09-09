/** \file ZebraLineGS.cpp
 * \author Andre Oeckerath
 * \brief ZebraLineGS.cpp contains the implementation of the class ZebraLineGS.
 * \see ZebraLineGS.h
 * \todo clean up ninepointxzebra, ninepointyzebra, xzebra and yzebra by
 *       doing redudant things in seperate functions
 */

#include <iostream>
#include "ZebraLineGS.h"

namespace mg
{
void ZebraLineGS::ninePointX(
    NumericArray &u,
    const NumericArray &f, 
    const Stencil &stencil,
    const Index nx, 
    const Index ny) const
{ 
    if (nx > 2)
    {
        solveNinePointXLine(u,f,u,stencil,nx,ny,1,1,SW,S,SE);

        for(Index sy=3; sy<ny-1 ; sy+=2)
            solveNinePointXLine(u,f,u,stencil,nx,ny,1,sy,W,C,E);

        solveNinePointXLine(u,f,u,stencil,nx,ny,1,ny-1,NW,N,NE);

        for(Index sy=2; sy<ny ; sy+=2) 
            solveNinePointXLine(u,f,u,stencil,nx,ny,1,sy,W,C,E);
    }       
    else
    {
        gsRedBlack_.relax(u,f,stencil,nx,ny);
    }
}

void ZebraLineGS::ninePointY(
    NumericArray &u,
    const NumericArray &f, 
    const Stencil &stencil,
    const Index nx, 
    const Index ny) const
{ 
    if (ny > 2)
    {
        solveNinePointYLine(u, f, u, stencil, nx, ny, 1, 1, SW, W, NW);

        for(Index sx=3; sx<nx-1 ; sx+=2)
            solveNinePointYLine(u, f, u, stencil, nx, ny, sx, 1, S, C, N);

        solveNinePointYLine(u, f, u, stencil, nx, ny, nx-1, 1, SE, E, NE);

        for(Index sx=2; sx<nx ; sx+=2)
            solveNinePointYLine(u, f, u, stencil, nx, ny, sx, 1, S, C, N);
    }
    else
    {
        gsRedBlack_.relax(u,f,stencil,nx,ny);
    }                       
}

void ZebraLineGS::fullX(
    NumericArray &u,
    const NumericArray &f, 
    const Stencil &stencil,
    const Index nx, 
    const Index ny) const
{
    if((ny > 4) && (nx > 4))
    {   
        solveFullXLine(u, f, u, stencil, nx, ny, 1, 1, SW, S, SE);

        for(Index sy=3; sy < ny-2; sy+=2)
            solveFullXLine(u, f, u, stencil, nx, ny, 1, sy, W, C, E);

        solveFullXLine(u, f, u, stencil, nx, ny, 1, ny-1, NW, N, NE);

        for(Index sy=2; sy<ny-1; sy+=2)
            solveFullXLine(u, f, u, stencil, nx, ny, 1, sy, W, C, E);
    }
    else //nx,ny to small
    {
        for(int i=0; i<2; i++)
        {
            gsRedBlack_.relax(u,f,stencil,nx,ny);
        }
    }
}

void ZebraLineGS::fullY(
    NumericArray &u,
    const NumericArray &f, 
    const Stencil &stencil,
    const Index nx, 
    const Index ny) const
{
    if ((ny > 4) && (nx > 4))
    {
        solveFullYLine(u, f, u, stencil, nx, ny, 1, 1, SW, W, NW);

        for(Index sx=3; sx<nx-2; sx+=2)
            solveFullYLine(u, f, u, stencil, nx, ny, sx, 1, S, C, N);

        solveFullYLine(u, f, u, stencil, nx, ny, nx-1, 1, SE, E, NE);

        for(Index sx=2; sx<nx-1; sx+=2)
            solveFullYLine(u, f, u, stencil, nx, ny, sx, 1, S, C, N);
    }
    else
    {
        for(int i=0; i<2; i++)
        {
            gsRedBlack_.relax(u,f,stencil,nx,ny);
        }
    }
}
}
